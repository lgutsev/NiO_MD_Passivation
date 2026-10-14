"""Unit systems and energy conventions, treated as first-class metadata.

Two independent things can make two correct engines disagree about the same
model and the same geometry:

1. **Units.** ASE reports eV and eV/Angstrom. LAMMPS ``metal`` units happen to
   agree; LAMMPS ``real`` units do not (kcal/mol). OpenMM reports kJ/mol and
   nm. Comparing a raw OpenMM number with a raw ASE number is meaningless.
2. **Energy convention.** A machine-learned potential's total energy usually
   includes per-element atomic self-energies (the isolated-atom references,
   ``E0``). OpenMM-ML's MACE implementation defaults to the *interaction*
   energy, with those self-energies removed. For a few thousand atoms the
   difference is thousands of eV, so an ASE/OpenMM comparison can look
   catastrophically wrong while the forces agree perfectly -- forces are
   unaffected, because a per-element constant has zero gradient.

Everything in this subsystem is normalised to the canonical system (eV,
eV/Angstrom, Angstrom, eV/Angstrom^3) and to a declared convention before any
two numbers are allowed to meet. Engines still run in their own native units.

Conversion factors are built from the exact SI-defining constants (the
elementary charge and the Avogadro constant have had exact values since the
2019 SI redefinition) rather than copied from a table, so they cannot drift.
``tests/test_mlip_units.py`` cross-checks them against ``ase.units``.
"""
from __future__ import annotations

from dataclasses import dataclass

from .errors import EnergyConventionError, UnitError

# Exact SI-defining constants.
ELEMENTARY_CHARGE_C = 1.602176634e-19
AVOGADRO_PER_MOL = 6.02214076e23

#: 1 kJ/mol expressed in eV.
KJ_PER_MOL_IN_EV = 1000.0 / AVOGADRO_PER_MOL / ELEMENTARY_CHARGE_C
#: 1 kcal/mol expressed in eV (thermochemical calorie, exactly 4.184 J).
KCAL_PER_MOL_IN_EV = 4.184 * KJ_PER_MOL_IN_EV
#: 1 nm expressed in Angstrom.
NM_IN_ANGSTROM = 10.0
#: 1 eV/Angstrom^3 expressed in GPa and in bar (the LAMMPS metal pressure unit).
EV_PER_ANGSTROM3_IN_GPA = ELEMENTARY_CHARGE_C * 1e30 / 1e9
EV_PER_ANGSTROM3_IN_BAR = ELEMENTARY_CHARGE_C * 1e30 / 1e5

QUANTITIES = ("energy", "length", "force", "stress")


@dataclass(frozen=True)
class UnitSystem:
    """A named (energy, length) pair plus the derived force and stress units.

    ``energy_in_eV`` and ``length_in_angstrom`` are the multiplicative factors
    that take a value *in this system* to the canonical system.
    """

    name: str
    energy: str
    length: str
    energy_in_eV: float
    length_in_angstrom: float

    @property
    def force(self) -> str:
        return f"{self.energy}/{self.length}"

    @property
    def stress(self) -> str:
        return f"{self.energy}/{self.length}^3"

    @property
    def force_in_eV_per_angstrom(self) -> float:
        return self.energy_in_eV / self.length_in_angstrom

    @property
    def stress_in_eV_per_angstrom3(self) -> float:
        return self.energy_in_eV / self.length_in_angstrom**3

    def factor(self, quantity: str) -> float:
        """Factor taking ``quantity`` from this system to the canonical one."""
        if quantity == "energy":
            return self.energy_in_eV
        if quantity == "length":
            return self.length_in_angstrom
        if quantity == "force":
            return self.force_in_eV_per_angstrom
        if quantity == "stress":
            return self.stress_in_eV_per_angstrom3
        raise UnitError(
            f"unknown quantity {quantity!r}; expected one of {', '.join(QUANTITIES)}"
        )

    def label(self, quantity: str) -> str:
        return {
            "energy": self.energy,
            "length": self.length,
            "force": self.force,
            "stress": self.stress,
        }[quantity]


#: The comparison system. Every result object this subsystem hands out is in
#: these units, whatever the engine used internally.
CANONICAL = UnitSystem("canonical", "eV", "Angstrom", 1.0, 1.0)
ASE = UnitSystem("ase", "eV", "Angstrom", 1.0, 1.0)
LAMMPS_METAL = UnitSystem("lammps_metal", "eV", "Angstrom", 1.0, 1.0)
LAMMPS_REAL = UnitSystem("lammps_real", "kcal/mol", "Angstrom", KCAL_PER_MOL_IN_EV, 1.0)
OPENMM = UnitSystem("openmm", "kJ/mol", "nm", KJ_PER_MOL_IN_EV, NM_IN_ANGSTROM)

UNIT_SYSTEMS: dict[str, UnitSystem] = {
    system.name: system for system in (CANONICAL, ASE, LAMMPS_METAL, LAMMPS_REAL, OPENMM)
}

#: LAMMPS ``units`` keyword -> the unit system it implies.
LAMMPS_UNIT_SYSTEMS = {"metal": LAMMPS_METAL, "real": LAMMPS_REAL}


def unit_system(name) -> UnitSystem:
    if isinstance(name, UnitSystem):
        return name
    try:
        return UNIT_SYSTEMS[name]
    except KeyError:
        raise UnitError(
            f"unknown unit system {name!r}; known: {', '.join(sorted(UNIT_SYSTEMS))}"
        ) from None


def lammps_unit_system(units: str) -> UnitSystem:
    """Map a LAMMPS ``units`` keyword onto a unit system.

    Only the styles this subsystem can honestly convert are accepted; an
    ``lj``-units model, for instance, has no absolute energy scale at all and
    must not be quietly treated as eV.
    """
    try:
        return LAMMPS_UNIT_SYSTEMS[units]
    except KeyError:
        raise UnitError(
            f"LAMMPS units {units!r} have no canonical (eV/Angstrom) mapping in this "
            f"subsystem; supported: {', '.join(sorted(LAMMPS_UNIT_SYSTEMS))}"
        ) from None


def convert(value, quantity: str, source, target=CANONICAL):
    """Convert a scalar or (nested) sequence between unit systems."""
    factor = unit_system(source).factor(quantity) / unit_system(target).factor(quantity)
    return _scale(value, factor)


def to_canonical(value, quantity: str, source):
    """Convert into eV / Angstrom / eV-per-Angstrom / eV-per-Angstrom^3."""
    return convert(value, quantity, source, CANONICAL)


def from_canonical(value, quantity: str, target):
    """Convert canonical values into an engine's native units."""
    return convert(value, quantity, CANONICAL, target)


def _scale(value, factor: float):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return value * factor
    if hasattr(value, "shape") and hasattr(value, "__mul__"):  # numpy array
        return value * factor
    return [_scale(item, factor) for item in value]


# --------------------------------------------------------------------------
# Energy conventions
# --------------------------------------------------------------------------

#: Energy includes the per-element atomic self-energies (isolated-atom E0s).
#: This is what a stock MACE ASE calculator and a MACE LAMMPS pair style report.
TOTAL = "total"
#: Atomic self-energies removed. This is OpenMM-ML's MACE default.
INTERACTION = "interaction"
ENERGY_CONVENTIONS = (TOTAL, INTERACTION)


def check_energy_convention(value: str, *, field: str = "energy_convention") -> str:
    if value not in ENERGY_CONVENTIONS:
        raise EnergyConventionError(
            f"{field} must be one of {', '.join(ENERGY_CONVENTIONS)}; got {value!r}"
        )
    return value


def self_energy_eV(symbols, atomic_reference_energies) -> float:
    """Sum of isolated-atom reference energies (E0) over a structure, in eV.

    This is exactly the offset between the two conventions. It is a
    composition-dependent constant, so it shifts energies while leaving forces
    and stresses untouched.
    """
    if not atomic_reference_energies:
        raise EnergyConventionError(
            "converting between the 'total' and 'interaction' energy conventions "
            "needs per-element atomic reference energies (E0), and none are known. "
            "Declare them as potential.atomic_reference_energies, or compare the two "
            "engines under their shared native convention instead."
        )
    missing = sorted({s for s in symbols if s not in atomic_reference_energies})
    if missing:
        raise EnergyConventionError(
            "no atomic reference energy (E0) for element(s) "
            f"{', '.join(missing)}; cannot convert energy convention"
        )
    return float(sum(atomic_reference_energies[s] for s in symbols))


def convert_energy_convention(
    energy_eV: float,
    symbols,
    *,
    source: str,
    target: str,
    atomic_reference_energies=None,
) -> float:
    """Shift an energy between the ``total`` and ``interaction`` conventions."""
    check_energy_convention(source, field="source")
    check_energy_convention(target, field="target")
    if source == target:
        return float(energy_eV)
    offset = self_energy_eV(symbols, atomic_reference_energies)
    if source == INTERACTION and target == TOTAL:
        return float(energy_eV) + offset
    return float(energy_eV) - offset


def describe_conventions() -> dict[str, str]:
    """Human-readable notes, reused by ``mlip inspect`` and by the manifest."""
    return {
        TOTAL: (
            "total energy including per-element atomic self-energies (E0); the "
            "convention of a stock MACE ASE calculator and of MACE LAMMPS pair styles"
        ),
        INTERACTION: (
            "interaction energy with atomic self-energies removed; the default of "
            "OpenMM-ML's MACE potential (returnEnergyType='interaction_energy')"
        ),
    }


__all__ = [
    "UnitSystem",
    "CANONICAL",
    "ASE",
    "LAMMPS_METAL",
    "LAMMPS_REAL",
    "OPENMM",
    "UNIT_SYSTEMS",
    "LAMMPS_UNIT_SYSTEMS",
    "QUANTITIES",
    "KJ_PER_MOL_IN_EV",
    "KCAL_PER_MOL_IN_EV",
    "NM_IN_ANGSTROM",
    "EV_PER_ANGSTROM3_IN_GPA",
    "EV_PER_ANGSTROM3_IN_BAR",
    "unit_system",
    "lammps_unit_system",
    "convert",
    "to_canonical",
    "from_canonical",
    "TOTAL",
    "INTERACTION",
    "ENERGY_CONVENTIONS",
    "check_energy_convention",
    "self_energy_eV",
    "convert_energy_convention",
    "describe_conventions",
]
