"""Capability negotiation and element coverage.

These are the checks that must reject a job *before* a batch script exists:
NPT without a trustworthy virial, per-atom energies from a route that does not
decompose them, and a structure carrying an element the model never saw.
"""
import pytest

from nio_md_prep.mlip.capabilities import (
    CapabilitySet,
    RequirementSet,
    check_elements,
    negotiate,
)
from nio_md_prep.mlip.errors import CapabilityError, ConfigError, ElementCoverageError
from nio_md_prep.mlip.specs import EngineSpec, MacePotentialSpec, SimulationSpec

FULL = CapabilitySet(
    energy=True,
    forces=True,
    stress=True,
    per_atom_energy=True,
    periodic=True,
    gpu=True,
    elements=frozenset({"Ni", "O", "P", "C", "H"}),
)


def test_intersection_keeps_only_what_both_sides_have():
    potential = CapabilitySet(energy=True, forces=True, stress=True, per_atom_energy=True)
    engine = CapabilitySet(energy=True, forces=True, stress=False, per_atom_energy=False)
    combined = potential.intersect(engine)
    assert combined.forces and not combined.stress and not combined.per_atom_energy


def test_unrestricted_element_set_is_the_identity():
    """An engine does not restrict elements; a potential does."""
    potential = CapabilitySet(elements=frozenset({"Ni", "O"}))
    engine = CapabilitySet(elements=None)
    assert potential.intersect(engine).elements == frozenset({"Ni", "O"})


def test_npt_requires_a_stress_route():
    simulation = SimulationSpec(
        task="md",
        ensemble="npt",
        temperature_K=300,
        pressure_bar=1.0,
        timestep_fs=0.5,
        steps=10,
    )
    requirements = simulation.required_capabilities()
    assert requirements.stress
    assert "virial" in requirements.reasons["stress"]


def test_npt_on_a_stressless_route_is_rejected_with_a_reason():
    """The OpenMM case: no virial, so constant pressure must fail here."""
    simulation = SimulationSpec(
        task="md",
        ensemble="npt",
        temperature_K=300,
        pressure_bar=1.0,
        timestep_fs=0.5,
        steps=10,
    )
    openmm_like = CapabilitySet(energy=True, forces=True, stress=False, periodic=True)
    with pytest.raises(CapabilityError) as excinfo:
        negotiate(openmm_like, simulation.required_capabilities(), label="mace -> openmm")
    assert "stress" in str(excinfo.value)
    assert "integrates the simulation cell" in str(excinfo.value)


def test_optimisation_requires_forces():
    requirements = SimulationSpec(task="optimize").required_capabilities()
    assert requirements.forces
    energy_only = CapabilitySet(energy=True, forces=False)
    with pytest.raises(CapabilityError, match="forces"):
        negotiate(energy_only, requirements, label="energy-only route")


def test_per_atom_energy_request_is_rejected_where_unavailable():
    simulation = SimulationSpec(task="singlepoint", compute_per_atom_energy=True)
    no_decomposition = CapabilitySet(energy=True, forces=True, per_atom_energy=False)
    with pytest.raises(CapabilityError, match="per_atom_energy"):
        negotiate(no_decomposition, simulation.required_capabilities(), label="route")


def test_a_periodic_structure_requires_a_periodic_route():
    requirements = SimulationSpec().required_capabilities(periodic=True)
    aperiodic = CapabilitySet(energy=True, forces=True, periodic=False)
    with pytest.raises(CapabilityError, match="periodic"):
        negotiate(aperiodic, requirements, label="route")


def test_element_coverage_is_checked_before_anything_runs():
    """A Ni/O/P/C/H model must not silently receive an extra atom type."""
    with pytest.raises(ElementCoverageError) as excinfo:
        check_elements(FULL, ["Ni", "O", "P", "C", "H", "F"], label="nio_phosphonate")
    message = str(excinfo.value)
    assert "F" in message
    assert "nio_phosphonate" in message
    assert "Ni" in message  # the message must say what IS covered


def test_element_coverage_accepts_a_subset():
    check_elements(FULL, ["Ni", "O"], label="nio_phosphonate")


def test_precision_mismatch_is_reported():
    requirements = RequirementSet(precision="float64")
    float32_only = CapabilitySet(precisions=frozenset({"float32"}))
    assert any("precision" in p for p in requirements.unmet(float32_only))


def test_an_unreachable_energy_convention_is_reported():
    requirements = SimulationSpec(
        task="singlepoint", energy_convention="total"
    ).required_capabilities()
    interaction_only = CapabilitySet(
        energy=True, forces=True, native_energy_convention="interaction"
    )
    problems = requirements.unmet(interaction_only)
    assert any("energy convention" in p for p in problems)


def test_a_route_with_reference_energies_can_reach_both_conventions():
    convertible = CapabilitySet(
        energy=True,
        forces=True,
        native_energy_convention="interaction",
        convertible_energy_conventions=frozenset({"total", "interaction"}),
    )
    requirements = SimulationSpec(
        task="singlepoint", energy_convention="total"
    ).required_capabilities()
    assert not requirements.unmet(convertible)


def test_mace_advertises_gpu_only_when_the_device_asks_for_one():
    from nio_md_prep.mlip.potentials.mace import MaceAdapter

    cpu = MaceAdapter(
        MacePotentialSpec(declared_elements=("Ni",), device="cpu")
    ).capabilities()
    cuda = MaceAdapter(
        MacePotentialSpec(declared_elements=("Ni",), device="cuda:1")
    ).capabilities()
    assert not cpu.gpu and cuda.gpu


def test_openmm_engine_declares_no_stress():
    """Pinned deliberately: this is what makes the NPT rejection happen."""
    from nio_md_prep.mlip.engines.openmm_engine import OpenMMEngine

    capabilities = OpenMMEngine(EngineSpec(kind="openmm")).capabilities()
    assert capabilities.energy and capabilities.forces
    assert not capabilities.stress
    assert not capabilities.per_atom_energy


def test_capability_notes_are_not_repeated():
    combined = CapabilitySet(notes=("same",)).intersect(CapabilitySet(notes=("same",)))
    assert combined.notes == ("same",)


# --- per-axis periodicity, frozen atoms and geometry refusals --------------

SLAB = (True, True, False)


def npt(**kwargs):
    base = dict(task="md", ensemble="npt", temperature_K=300, pressure_bar=1.0,
                timestep_fs=0.5, steps=10)
    base.update(kwargs)
    return SimulationSpec(**base)


def test_a_slab_requires_per_axis_periodicity_not_just_periodic():
    requirements = SimulationSpec().required_capabilities(pbc=SLAB)
    assert requirements.periodic and requirements.partial_periodic
    assert requirements.periodic_axes == SLAB
    assert requirements.as_dict()["periodic_axes"] == [True, True, False]
    all_or_nothing = CapabilitySet(energy=True, forces=True, periodic=True)
    with pytest.raises(CapabilityError, match="partial_periodic"):
        negotiate(all_or_nothing, requirements, label="openmm-like")


def test_bulk_and_cluster_need_no_per_axis_support():
    bulk = SimulationSpec().required_capabilities(pbc=(True, True, True))
    cluster = SimulationSpec().required_capabilities(pbc=(False, False, False))
    assert bulk.periodic and not bulk.partial_periodic
    assert not cluster.periodic and not cluster.partial_periodic
    # The older boolean spelling still means "fully periodic".
    assert SimulationSpec().required_capabilities(periodic=True).periodic_axes == (
        True, True, True
    )


def test_engine_flags_come_from_the_engine_side_of_an_intersection():
    """A potential never declares slab/FixAtoms support; the engine does."""
    potential = CapabilitySet(energy=True, forces=True, periodic=True)
    engine = CapabilitySet(
        energy=True, forces=True, periodic=True, partial_periodic=True, fixed_atoms=True
    )
    combined = potential.intersect(engine)
    assert combined.partial_periodic and combined.fixed_atoms
    assert not engine.intersect(CapabilitySet()).fixed_atoms


@pytest.mark.parametrize(
    "simulation",
    [SimulationSpec(compute_stress=True), npt(), npt(barostat_coupling="isotropic")],
    ids=["stress", "npt-default", "npt-isotropic"],
)
def test_stress_and_npt_need_a_fully_periodic_cell(simulation):
    with pytest.raises(ConfigError, match="periodic along all three axes"):
        simulation.required_capabilities(pbc=SLAB)


def test_in_plane_npt_on_a_slab_is_left_to_the_engine():
    requirements = npt(barostat_coupling="in-plane").required_capabilities(pbc=SLAB)
    assert requirements.stress and requirements.partial_periodic
    with pytest.raises(ConfigError, match="at least two periodic axes"):
        npt(barostat_coupling="in-plane").required_capabilities(pbc=(True, False, False))


def test_npt_on_a_vacuum_slab_written_as_bulk_is_refused():
    with pytest.raises(ConfigError, match="vacuum slab") as excinfo:
        npt().required_capabilities(pbc=(True, True, True), vacuum_gaps={2: 15.0})
    assert "15.00 Angstrom along lattice vector c" in str(excinfo.value)
    # Declaring in-plane coupling is the way through.
    npt(barostat_coupling="in-plane").required_capabilities(
        pbc=(True, True, True), vacuum_gaps={2: 15.0}
    )


def test_npt_with_frozen_atoms_is_refused_everywhere():
    with pytest.raises(ConfigError, match="frozen atom"):
        npt().required_capabilities(pbc=(True, True, True), fixed_atoms=4)


def test_frozen_atoms_need_an_engine_that_honours_them_in_dynamics():
    nvt = SimulationSpec(task="md", ensemble="nvt", temperature_K=300, timestep_fs=0.5,
                         steps=10)
    requirements = nvt.required_capabilities(pbc=(True, True, True), fixed_atoms=4)
    assert requirements.fixed_atoms and "4 atom(s)" in requirements.reasons["fixed_atoms"]
    # A single point evaluates raw forces on every atom; no engine support needed.
    assert not SimulationSpec().required_capabilities(fixed_atoms=4).fixed_atoms


def mock_bridge(engine=None):
    from nio_md_prep.mlip.registry import build_bridge
    from nio_md_prep.mlip.specs import MockPotentialSpec

    return build_bridge(
        MockPotentialSpec(declared_elements=("Ni", "O")), engine or EngineSpec(kind="ase")
    )


def test_bridge_requirements_read_the_geometry_from_the_structure(nio_structure):
    from ase.constraints import FixAtoms

    atoms = nio_structure.copy()
    atoms.pbc = SLAB
    atoms.set_constraint(FixAtoms(indices=[0, 1]))
    nvt = SimulationSpec(task="md", ensemble="nvt", temperature_K=300, timestep_fs=0.5,
                         steps=10)
    requirements = mock_bridge().requirements(nvt, atoms)
    assert requirements.periodic_axes == SLAB
    assert requirements.partial_periodic and requirements.fixed_atoms
    # Engine kind and engine.precision are negotiated, not dead fields.
    assert requirements.engine == "ase"
    assert requirements.precision is None


def test_the_bridge_refuses_an_unsupported_constraint_at_validate(nio_structure):
    from ase.constraints import FixCom

    atoms = nio_structure.copy()
    atoms.set_constraint(FixCom())
    with pytest.raises(ConfigError, match="FixCom"):
        mock_bridge().validate(SimulationSpec(), atoms)


def test_engine_precision_is_a_negotiated_cross_check(nio_structure):
    """The mock evaluates in float64; engine.precision = float32 contradicts it."""
    bridge = mock_bridge(EngineSpec(kind="ase", precision="float32"))
    with pytest.raises(CapabilityError, match="precision 'float32' is not available"):
        bridge.validate(SimulationSpec(), nio_structure)
    mock_bridge(EngineSpec(kind="ase", precision="float64")).validate(
        SimulationSpec(), nio_structure
    )


def test_engine_capabilities_declare_per_axis_and_frozen_atom_support():
    """Pinned: engine agents flip these only when the engine really honours them."""
    from nio_md_prep.mlip.engines.ase_engine import AseEngine
    from nio_md_prep.mlip.engines.openmm_engine import OpenMMEngine

    ase = AseEngine(EngineSpec(kind="ase")).capabilities()
    assert ase.partial_periodic and ase.fixed_atoms
    openmm = OpenMMEngine(EngineSpec(kind="openmm")).capabilities()
    assert not openmm.partial_periodic  # openmm-ml periodicity is all-or-nothing
