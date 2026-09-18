"""LAMMPS deck rendering, which is pure and therefore testable everywhere.

The rendered ``pair_style``/``pair_coeff`` strings are what a LAMMPS run
actually is, and they are what the manifest preserves verbatim. Rendering them
without LAMMPS installed is deliberate: ``mlip validate`` must be able to show
the deck it would submit from a laptop.
"""
import pytest

from nio_md_prep.mlip.engines import lammps_engine
from nio_md_prep.mlip.registry import build_bridge
from nio_md_prep.mlip.specs import (
    EngineSpec,
    LammpsMlipPotentialSpec,
    MacePotentialSpec,
    SimulationSpec,
)

DEEPMD = LammpsMlipPotentialSpec(
    label="nio-deepmd",
    pair_style="deepmd nio.pb",
    pair_coeff=("* *",),
    type_map={1: "Ni", 2: "O"},
    units="metal",
    framework="deepmd",
)

SINGLEPOINT = SimulationSpec(task="singlepoint", compute_per_atom_energy=True)
NVT = SimulationSpec(
    task="md",
    ensemble="nvt",
    temperature_K=400,
    timestep_fs=0.5,
    steps=100,
    seed=12345,
    thermostat_damping_fs=100.0,
    trajectory_interval=10,
)


def render(potential=DEEPMD, simulation=SINGLEPOINT, **kwargs):
    return lammps_engine.render_deck(
        potential, simulation, n_types=len(potential.type_map), **kwargs
    )


def test_the_deck_sets_units_and_atom_style_before_reading_data():
    deck = render()
    assert "units metal" in deck
    assert "atom_style atomic" in deck
    assert deck.index("units metal") < deck.index(f"read_data {lammps_engine.DATA_FILE}")


def test_pair_commands_appear_verbatim_after_the_structure():
    deck = render()
    assert "pair_style deepmd nio.pb" in deck
    assert "pair_coeff * *" in deck
    assert deck.index(f"read_data {lammps_engine.DATA_FILE}") < deck.index(
        "pair_style deepmd nio.pb"
    )


def test_a_singlepoint_deck_runs_zero_steps_and_dumps_forces():
    deck = render()
    assert "run 0" in deck
    dump = next(line for line in deck if line.startswith("dump mlip_forces"))
    assert "fx fy fz" in dump and "c_mlip_pe_atom" in dump
    assert "dump_modify mlip_forces sort id" in deck


def test_metal_timestep_is_picoseconds_not_femtoseconds():
    """A factor of 1000 here would silently ruin every trajectory."""
    deck = render(simulation=NVT)
    assert "timestep 0.0005" in deck
    assert lammps_engine.TIMESTEP_PER_FS["metal"] == 1e-3
    assert lammps_engine.TIMESTEP_PER_FS["real"] == 1.0


def test_real_units_use_femtoseconds():
    real = LammpsMlipPotentialSpec(
        pair_style="x", pair_coeff=("* *",), type_map={1: "Ni"}, units="real"
    )
    assert "timestep 0.5" in render(potential=real, simulation=NVT)


def test_nvt_renders_a_thermostat_and_a_seeded_velocity_command():
    deck = render(simulation=NVT)
    assert any(line.startswith("fix mlip_integrate all nvt temp 400 400") for line in deck)
    assert any("velocity all create 400 12345" in line for line in deck)
    assert "run 100" in deck


def test_npt_renders_a_barostat():
    npt = SimulationSpec(
        task="md",
        ensemble="npt",
        temperature_K=400,
        pressure_bar=1.0,
        timestep_fs=0.5,
        steps=50,
    )
    deck = render(simulation=npt)
    fix = next(line for line in deck if line.startswith("fix mlip_integrate"))
    assert "npt temp 400 400" in fix and "iso 1 1" in fix


def test_nve_renders_no_thermostat():
    nve = SimulationSpec(task="md", ensemble="nve", timestep_fs=0.5, steps=10)
    deck = render(simulation=nve)
    assert "fix mlip_integrate all nve" in deck
    assert not any("nvt" in line for line in deck)


def test_a_non_periodic_structure_gets_fixed_boundaries():
    assert "boundary f f f" in render(periodic=False)
    assert "boundary p p p" in render(periodic=True)


def test_thermo_values_are_printed_for_the_executable_route():
    deck = render()
    assert any('print "pe ${mlip_pe}" file' in line for line in deck)
    assert any('print "pxx ${mlip_pxx}" append' in line for line in deck)


def test_stress_conversion_flips_the_sign_and_the_unit():
    """LAMMPS reports pressure in bars; canonical stress is eV/Angstrom^3."""
    from nio_md_prep.mlip.units import EV_PER_ANGSTROM3_IN_BAR

    thermo = {"pxx": 10000.0, "pyy": 0.0, "pzz": 0.0, "pyz": 0.0, "pxz": 0.0, "pxy": 0.0}
    stress = lammps_engine.stress_from_pressure(thermo, "metal")
    assert stress[0] == pytest.approx(-10000.0 / EV_PER_ANGSTROM3_IN_BAR)
    assert stress[0] < 0  # positive pressure is compressive, i.e. negative stress


def test_stress_is_not_invented_for_units_with_a_different_pressure_scale():
    thermo = {k: 1.0 for k in ("pxx", "pyy", "pzz", "pyz", "pxz", "pxy")}
    assert lammps_engine.stress_from_pressure(thermo, "real") is None


def test_incomplete_thermo_yields_no_stress():
    assert lammps_engine.stress_from_pressure({"pxx": 1.0}, "metal") is None


# --- MACE projected onto the LAMMPS representation ------------------------


def mace_spec(tmp_path):
    return MacePotentialSpec(
        label="nio",
        model_path=tmp_path / "nio.model",
        declared_elements=("Ni", "O", "P", "C", "H"),
    )


def test_mace_on_lammps_renders_the_mliap_pair_style(tmp_path):
    bridge = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps"))
    commands = bridge.engine_parameters()["pair_commands"]
    assert commands[0].startswith("pair_style mliap unified nio-mliap.pt 0")
    assert commands[1] == "pair_coeff * * Ni O P C H"


def test_mace_on_lammps_can_use_the_older_pair_style(tmp_path):
    bridge = build_bridge(
        mace_spec(tmp_path), EngineSpec(kind="lammps"), implementation="pair-mace"
    )
    commands = bridge.engine_parameters()["pair_commands"]
    assert commands[0] == "pair_style mace no_domain_decomposition"
    assert commands[1] == "pair_coeff * * nio-lammps.pt Ni O P C H"


def test_only_the_mliap_route_claims_per_atom_energies(tmp_path):
    """ML-IAP brings atomic virials and site energies; the older route is not assumed to."""
    spec = mace_spec(tmp_path)
    mliap = build_bridge(spec, EngineSpec(kind="lammps"))
    pair_mace = build_bridge(
        spec, EngineSpec(kind="lammps"), implementation="pair-mace"
    )
    assert mliap.capabilities().per_atom_energy
    assert not pair_mace.capabilities().per_atom_energy


def test_the_exported_model_path_can_be_overridden(tmp_path):
    engine = EngineSpec(kind="lammps", options={"model_path": "/models/custom.pt"})
    bridge = build_bridge(mace_spec(tmp_path), engine)
    assert bridge.exported_model_path().name == "custom.pt"


def test_a_mace_lammps_route_renders_a_complete_deck(tmp_path):
    bridge = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps"))
    spec = bridge.as_lammps_spec()
    deck = lammps_engine.render_deck(spec, SINGLEPOINT, n_types=5)
    assert "units metal" in deck
    assert any(line.startswith("pair_style mliap unified") for line in deck)
    assert spec.type_map == {1: "Ni", 2: "O", 3: "P", 4: "C", 5: "H"}
