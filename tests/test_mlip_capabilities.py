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
from nio_md_prep.mlip.errors import CapabilityError, ElementCoverageError
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
