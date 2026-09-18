"""Unit systems and energy conventions.

The conversion factors are built from exact SI constants rather than copied
from a table, so these tests cross-check them against an independent source
(``ase.units``) and pin the behaviour that matters scientifically: energies
shift between conventions, forces do not.
"""
import math

import pytest

from nio_md_prep.mlip import units as u
from nio_md_prep.mlip.errors import EnergyConventionError, UnitError

ase_units = pytest.importorskip("ase.units")

#: ASE's module-level units are CODATA 2014 by default, while this subsystem
#: builds its factors from the exact post-2019 SI constants. The two differ in
#: the ninth significant figure, which is far below anything that matters
#: physically but far above a float comparison -- so the cross-check is made
#: against ASE's own CODATA 2018 set, and the vintage gap is asserted
#: separately rather than papered over with a loose tolerance.
CODATA_2018 = ase_units.create_units("2018")


def test_energy_factors_agree_with_ase():
    assert u.KJ_PER_MOL_IN_EV == pytest.approx(
        CODATA_2018["kJ"] / CODATA_2018["mol"], rel=1e-12
    )
    assert u.KCAL_PER_MOL_IN_EV == pytest.approx(
        CODATA_2018["kcal"] / CODATA_2018["mol"], rel=1e-12
    )


def test_stress_factors_agree_with_ase():
    assert u.EV_PER_ANGSTROM3_IN_GPA == pytest.approx(1.0 / CODATA_2018["GPa"], rel=1e-12)
    # 1 GPa is 1e4 bar, so the bar factor must be exactly 1e4 times the GPa one.
    assert u.EV_PER_ANGSTROM3_IN_BAR == pytest.approx(
        u.EV_PER_ANGSTROM3_IN_GPA * 1e4, rel=1e-12
    )


def test_the_codata_vintage_gap_is_negligible_but_real():
    """Documents why the cross-check pins a vintage instead of a tolerance."""
    relative = abs(u.EV_PER_ANGSTROM3_IN_GPA - 1.0 / ase_units.GPa) / (
        1.0 / ase_units.GPa
    )
    assert 0 < relative < 1e-7


def test_openmm_round_trip_is_lossless():
    """kJ/mol/nm out and back must not drift; this is the ASE/OpenMM boundary."""
    forces = [[1.0, -2.0, 3.0], [0.5, 0.0, -0.25]]
    canonical = u.to_canonical(forces, "force", u.OPENMM)
    back = u.from_canonical(canonical, "force", u.OPENMM)
    for original_row, returned_row in zip(forces, back):
        for original, returned in zip(original_row, returned_row):
            assert returned == pytest.approx(original, rel=1e-14)


def test_openmm_force_conversion_includes_the_length_unit():
    """A kJ/mol/nm force is not a kJ/mol energy: the nm must be converted too."""
    energy_factor = u.to_canonical(1.0, "energy", u.OPENMM)
    force_factor = u.to_canonical(1.0, "force", u.OPENMM)
    assert force_factor == pytest.approx(energy_factor / 10.0, rel=1e-14)


def test_lammps_metal_is_already_canonical():
    assert u.LAMMPS_METAL.energy_in_eV == 1.0
    assert u.LAMMPS_METAL.length_in_angstrom == 1.0
    assert u.to_canonical(7.5, "energy", u.LAMMPS_METAL) == 7.5


def test_lammps_real_is_not_canonical():
    """The mistake this catches: assuming every LAMMPS run reports eV."""
    assert u.to_canonical(1.0, "energy", u.LAMMPS_REAL) == pytest.approx(
        u.KCAL_PER_MOL_IN_EV
    )


def test_lj_units_are_refused_rather_than_assumed():
    with pytest.raises(UnitError, match="no canonical"):
        u.lammps_unit_system("lj")


def test_unknown_quantity_is_refused():
    with pytest.raises(UnitError, match="unknown quantity"):
        u.CANONICAL.factor("temperature")


# --- energy conventions ---------------------------------------------------

E0 = {"Ni": -5.78, "O": -4.95, "H": -1.12}


def test_self_energy_is_the_offset_between_conventions():
    symbols = ["Ni", "Ni", "O", "H"]
    offset = u.self_energy_eV(symbols, E0)
    assert offset == pytest.approx(2 * E0["Ni"] + E0["O"] + E0["H"])


def test_interaction_to_total_adds_the_self_energy():
    symbols = ["Ni", "O"]
    interaction = -1.5
    total = u.convert_energy_convention(
        interaction,
        symbols,
        source=u.INTERACTION,
        target=u.TOTAL,
        atomic_reference_energies=E0,
    )
    assert total == pytest.approx(interaction + E0["Ni"] + E0["O"])
    back = u.convert_energy_convention(
        total,
        symbols,
        source=u.TOTAL,
        target=u.INTERACTION,
        atomic_reference_energies=E0,
    )
    assert back == pytest.approx(interaction)


def test_conversion_without_reference_energies_refuses_rather_than_guesses():
    with pytest.raises(EnergyConventionError, match="atomic reference energies"):
        u.convert_energy_convention(
            1.0, ["Ni"], source=u.INTERACTION, target=u.TOTAL
        )


def test_conversion_refuses_an_element_with_no_reference_energy():
    with pytest.raises(EnergyConventionError, match="no atomic reference energy"):
        u.convert_energy_convention(
            1.0,
            ["Ni", "P"],
            source=u.INTERACTION,
            target=u.TOTAL,
            atomic_reference_energies=E0,
        )


def test_same_convention_needs_no_reference_energies():
    assert u.convert_energy_convention(
        3.0, ["Ni"], source=u.TOTAL, target=u.TOTAL
    ) == 3.0


def test_unknown_convention_is_rejected():
    with pytest.raises(EnergyConventionError, match="must be one of"):
        u.check_energy_convention("binding")


def test_convention_offset_is_large_enough_to_look_like_a_bug():
    """Why this matters: on a realistic cell the offset dwarfs any real error.

    A 3000-atom NiO slab carries thousands of eV of atomic self-energy. An
    ASE/OpenMM comparison that ignores the convention therefore looks
    catastrophically broken even when the forces agree exactly.
    """
    symbols = ["Ni", "O"] * 1500
    offset = u.self_energy_eV(symbols, E0)
    assert abs(offset) > 1000.0
    assert math.isfinite(offset)
