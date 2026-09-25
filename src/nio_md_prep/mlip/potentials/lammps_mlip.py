"""LAMMPS-native machine-learned potentials, as a family rather than a framework.

"LAMMPS MLIP" here means *any* pair style that evaluates a machine-learned
model inside LAMMPS. DeepMD is one such style; so are ML-IAP, PACE, and a
MACE pair style. Nothing in this module branches on which one it is: the
pair style string is the datum, and :data:`KNOWN_FRAMEWORKS` exists only to
suggest which LAMMPS package a route needs and to label provenance.

Adding a new framework therefore means adding a row to
:data:`KNOWN_FRAMEWORKS` -- or nothing at all, if the user is willing to
declare ``required_packages`` themselves.
"""
from __future__ import annotations

from collections.abc import Mapping

from ..capabilities import CapabilitySet
from ..environment import Availability, executable_available, module_available
from ..errors import ConfigError
from ..specs import LammpsMlipPotentialSpec
from ..units import lammps_unit_system
from .base import PotentialAdapter


#: Hints keyed by the first token of ``pair_style``. Only hints: an unknown
#: pair style is accepted, and simply carries no inferred package.
KNOWN_FRAMEWORKS: dict[str, dict] = {
    "deepmd": {
        "framework": "deepmd",
        "packages": ("USER-DEEPMD / DEEPMD",),
        "per_atom_energy": True,
        "stress": True,
    },
    "mliap": {
        "framework": "mliap",
        "packages": ("ML-IAP", "ML-SNAP"),
        "per_atom_energy": True,
        "stress": True,
        "note": "the ML-IAP interface also carries MACE models; atomic virials available",
    },
    "pace": {
        "framework": "pace",
        "packages": ("ML-PACE",),
        "per_atom_energy": True,
        "stress": True,
    },
    "mace": {
        "framework": "mace",
        "packages": ("a LAMMPS build with the MACE pair style",),
        # pair_mace tallies eatom, but MACE site energies through LAMMPS have
        # not been demonstrated against a reference: not claimed.
        "per_atom_energy": False,
        "stress": True,
        "note": "per-atom energies are not claimed for pair_style mace (not demonstrated)",
    },
    "nequip": {
        "framework": "nequip",
        "packages": ("pair_nequip",),
        "per_atom_energy": True,
        "stress": False,
        "note": "stress support depends on the pair_nequip build; declared conservatively",
    },
    "quip": {
        "framework": "quip",
        "packages": ("ML-QUIP",),
        "per_atom_energy": True,
        "stress": True,
    },
}


class LammpsMlipAdapter(PotentialAdapter):
    kind = "lammps"
    requires = ()

    spec: LammpsMlipPotentialSpec

    @property
    def framework(self) -> str | None:
        """The declared framework, else the first token of the pair style."""
        if self.spec.framework:
            return self.spec.framework
        token = self.spec.pair_style.split()[0] if self.spec.pair_style else ""
        return KNOWN_FRAMEWORKS.get(token, {}).get("framework") or token or None

    @property
    def hints(self) -> Mapping:
        tokens = self.spec.pair_style.split() if self.spec.pair_style else []
        token = tokens[0] if tokens else ""
        hints = dict(KNOWN_FRAMEWORKS.get(token, {}))
        if token == "mliap" and len(tokens) > 1 and tokens[1] == "unified":
            # A python-coupled ML-IAP model (a MACE export, typically): whether
            # its site energies are right has not been demonstrated here.
            hints.update(
                per_atom_energy=False,
                note=(
                    "pair_style mliap unified runs a python-coupled model (e.g. a MACE "
                    "export); per-atom energies are not claimed for it (not demonstrated)"
                ),
            )
        return hints

    def capabilities(self) -> CapabilitySet:
        """Advertise conservatively: a pair style cannot be interrogated from here.

        Energy, forces and the global virial are guaranteed by LAMMPS itself
        for any pair style. Per-atom energy is claimed only when the framework
        hint says the style computes ``eatom``, because a job asking for a
        per-atom decomposition must fail loudly rather than receive zeros.
        """
        hints = self.hints
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=bool(hints.get("stress", True)),
            per_atom_energy=bool(hints.get("per_atom_energy", False)),
            periodic=True,
            gpu=False,
            elements=frozenset(self.spec.elements),
            # Fixed by the pair style's package and model build and not readable
            # from here: no precision is guaranteed (an empty set, not None,
            # which would mean "any").
            precisions=frozenset(),
            engines=frozenset({"lammps", "ase"}),
            native_energy_convention=self.spec.energy_convention,
            native_units=lammps_unit_system(self.spec.units).name,
            notes=(
                f"pair_style {self.spec.pair_style}",
                *( (hints["note"],) if "note" in hints else () ),
            ),
        )

    def availability(self) -> Availability:
        """A LAMMPS build cannot be interrogated for pair styles from here.

        So this reports whether *some* LAMMPS is reachable (the Python module
        or an executable) and lists the packages the pair style needs, leaving
        the final word to LAMMPS itself at run time. Claiming more would be
        dishonest; claiming less would make ``mlip inspect`` useless.
        """
        missing = [] if module_available("lammps") else ["lammps"]
        detail = "the pair style's own package is verified by LAMMPS at run time"
        if missing and not executable_available("lmp"):
            return Availability(
                available=False,
                missing=tuple(missing),
                detail="neither the LAMMPS python module nor an 'lmp' executable was found",
            )
        return Availability(available=True, missing=(), detail=detail)

    def required_packages(self) -> tuple[str, ...]:
        if self.spec.required_packages:
            return self.spec.required_packages
        return tuple(self.hints.get("packages", ()))

    def missing_model_files(self) -> tuple[str, ...]:
        return tuple(str(p) for p in self.spec.model_paths if not p.exists())

    def inspect(self) -> dict:
        report = super().inspect()
        report.update(
            {
                "pair_style": self.spec.pair_style,
                "pair_coeff": list(self.spec.pair_coeff),
                "framework": self.framework,
                "units": self.spec.units,
                "unit_system": self.spec.unit_system_name,
                "atom_style": self.spec.atom_style,
                "type_map": {str(k): v for k, v in sorted(self.spec.type_map.items())},
                "required_packages": list(self.required_packages()),
                "model_paths": [str(p) for p in self.spec.model_paths],
                "model_files_missing": list(self.missing_model_files()),
                "rendered_commands": list(self.spec.render_pair_commands()),
                "availability": self.availability().as_dict(),
            }
        )
        return report


def build_adapter(spec: LammpsMlipPotentialSpec) -> LammpsMlipAdapter:
    if not isinstance(spec, LammpsMlipPotentialSpec):
        raise ConfigError(f"expected a LAMMPS MLIP potential spec; got {type(spec).__name__}")
    return LammpsMlipAdapter(spec)


__all__ = ["KNOWN_FRAMEWORKS", "LammpsMlipAdapter", "build_adapter"]
