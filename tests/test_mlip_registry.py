"""Bridge resolution: the compatibility matrix, enforced.

These tests are the executable statement of the matrix in
``docs/mlip_architecture.md``. The one that matters most is
:func:`test_lammps_potential_on_openmm_is_an_intentional_error` -- asking for
an impossible combination must produce a typed, explanatory error rather than
an import traceback from a half-written adapter.
"""
import sys

import pytest

from nio_md_prep.mlip import registry
from nio_md_prep.mlip.errors import (
    UnknownEngineError,
    UnknownPotentialError,
    UnsupportedCombinationError,
)
from nio_md_prep.mlip.specs import EngineSpec, MacePotentialSpec

SUPPORTED_CELLS = [
    ("mace", "ase", "mace-ase-calculator"),
    ("mace", "lammps", "mliap"),
    ("mace", "openmm", "openmm-ml"),
    ("lammps", "lammps", "native"),
    ("lammps", "ase", "ase-lammps"),
]


@pytest.mark.parametrize("potential,engine,implementation", SUPPORTED_CELLS)
def test_every_supported_cell_resolves(potential, engine, implementation):
    assert registry.resolve_bridge(potential, engine).implementation == implementation


def test_lammps_potential_on_openmm_is_an_intentional_error():
    """The single unsupported cell, and the reason it is unsupported."""
    with pytest.raises(UnsupportedCombinationError) as excinfo:
        registry.resolve_bridge("lammps", "openmm")
    message = str(excinfo.value)
    assert "compiled C++ inside LAMMPS" in message
    assert "own potential kind" in message
    # The error must name the routes that do work, so the message is actionable.
    assert "lammps -> lammps" in message
    assert "lammps -> ase" in message


def test_the_unsupported_cell_is_registered_rather_than_absent():
    """Registered-with-a-reason, not merely missing: the distinction is the point."""
    entries = registry.REGISTRY.for_pair("lammps", "openmm")
    assert entries and all(not e.supported for e in entries)
    assert entries[0].reason


def test_resolution_imports_no_backend():
    """Resolving must stay cheap: no torch, no OpenMM, no LAMMPS."""
    for potential, engine, _ in SUPPORTED_CELLS:
        registry.resolve_bridge(potential, engine)
    assert "torch" not in sys.modules
    assert "openmm" not in sys.modules
    assert "mace" not in sys.modules


def test_unknown_potential_and_engine_are_distinguished():
    with pytest.raises(UnknownPotentialError, match="unknown potential kind"):
        registry.resolve_bridge("dpa3", "ase")
    with pytest.raises(UnknownEngineError, match="unknown engine kind"):
        registry.resolve_bridge("mace", "gromacs")


def test_mliap_is_preferred_over_the_older_pair_style():
    """ML-IAP first: GPU acceleration, multi-GPU inference and atomic virials."""
    assert registry.resolve_bridge("mace", "lammps").implementation == "mliap"


def test_implementation_preference_selects_the_older_route():
    chosen = registry.resolve_bridge("mace", "lammps", preferences=("pair-mace",))
    assert chosen.implementation == "pair-mace"


def test_an_unknown_preference_falls_back_to_the_default():
    chosen = registry.resolve_bridge("mace", "lammps", preferences=("nonexistent",))
    assert chosen.implementation == "mliap"


def test_an_explicitly_requested_missing_implementation_is_an_error():
    with pytest.raises(UnsupportedCombinationError, match="no implementation named"):
        registry.resolve_bridge("mace", "ase", implementation="torchscript")


def test_build_bridge_reads_the_preference_from_the_spec(tmp_path):
    """The spec's implementation preference is honoured without a separate argument."""
    spec = MacePotentialSpec(
        model_path=tmp_path / "model.pt",
        declared_elements=("Ni", "O"),
        implementation=("pair-mace",),
    )
    bridge = registry.build_bridge(spec, EngineSpec(kind="lammps"))
    assert bridge.implementation == "pair-mace"


def test_matrix_covers_every_potential_and_engine():
    matrix = registry.compatibility_matrix()
    assert set(matrix) >= {"mace", "lammps"}
    for row in matrix.values():
        assert set(row) == {"ase", "lammps", "openmm"}
    # MACE reaches all three engines; a LAMMPS-native MLIP reaches two.
    assert all(matrix["mace"][engine] for engine in ("ase", "lammps", "openmm"))
    assert not any(e["status"] == "supported" for e in matrix["lammps"]["openmm"])


def test_duplicate_registration_is_refused():
    entry = registry.REGISTRY.for_pair("mace", "ase")[0]
    with pytest.raises(ValueError, match="already registered"):
        registry.REGISTRY.register(entry)
