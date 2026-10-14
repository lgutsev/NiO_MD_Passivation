"""ASE ``Atoms`` as the MLIP interchange structure, and the LAMMPS converter.

The converter exists so a prepared classical build can be read by the MLIP
layer. It is one-way on purpose: the topology-rich representation in
:mod:`nio_md_prep.lammps` stays where the classical workflows need it, and the
MLIP layer never pretends to reconstruct it.
"""
from pathlib import Path

import pytest

from nio_md_prep.mlip.errors import ConfigError
from nio_md_prep.mlip.structures import (
    describe_structure,
    from_lammps_data,
    is_periodic,
    load_structure,
    structure_digest,
    to_lammps_type_order,
)

pytest.importorskip("ase")

ROOT = Path(__file__).parents[1]
SURFACE = ROOT / "inputs/surfaces/corrugated-nio-110/surface.lmp"


def test_digest_is_stable_across_a_text_round_trip(tmp_path, nio_structure):
    from ase.io import read, write

    path = tmp_path / "cell.xyz"
    write(str(path), nio_structure, format="extxyz")
    assert structure_digest(read(str(path))) == structure_digest(nio_structure)


def test_digest_changes_when_an_atom_moves(nio_structure):
    moved = nio_structure.copy()
    positions = moved.get_positions()
    positions[0][0] += 0.01
    moved.set_positions(positions)
    assert structure_digest(moved) != structure_digest(nio_structure)


def test_digest_changes_when_the_cell_changes(nio_structure):
    strained = nio_structure.copy()
    strained.set_cell(strained.get_cell() * 1.01, scale_atoms=True)
    assert structure_digest(strained) != structure_digest(nio_structure)


def test_description_is_json_friendly_and_complete(nio_structure):
    import json

    report = describe_structure(nio_structure)
    assert report["n_atoms"] == len(nio_structure)
    assert report["elements"] == ["Ni", "O"]
    assert report["composition"] == {"Ni": 8, "O": 8}
    assert report["periodic"] is True
    assert len(report["sha256"]) == 64
    json.dumps(report)


def test_load_structure_reads_an_extxyz_file(tmp_path, nio_structure):
    from ase.io import write

    from nio_md_prep.mlip.specs import StructureSpec

    path = tmp_path / "cell.xyz"
    write(str(path), nio_structure, format="extxyz")
    loaded = load_structure(StructureSpec(path=path))
    assert len(loaded) == len(nio_structure)
    assert is_periodic(loaded)


def test_unknown_format_is_refused(tmp_path):
    path = tmp_path / "cell.weird"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ConfigError, match="structure.format"):
        load_structure(path, format="gromacs")


def test_missing_structure_reports_the_path(tmp_path):
    with pytest.raises(FileNotFoundError, match="structure not found"):
        load_structure(tmp_path / "absent.xyz")


@pytest.mark.skipif(not SURFACE.exists(), reason="the NiO surface fixture is absent")
def test_the_repository_lammps_representation_converts_without_being_replaced():
    """A prepared classical build becomes a valid MLIP interchange structure."""
    from nio_md_prep.lammps import parse

    atoms = from_lammps_data(SURFACE)
    data = parse(SURFACE)
    assert len(atoms) == data.count("Atoms")
    assert set(atoms.get_chemical_symbols()) == {"Ni", "O"}
    assert is_periodic(atoms)
    # Topology is deliberately dropped; it stays in the LAMMPS representation.
    assert not hasattr(atoms, "bonds")


@pytest.mark.skipif(not SURFACE.exists(), reason="the NiO surface fixture is absent")
def test_an_explicit_type_map_overrides_mass_inference():
    atoms = from_lammps_data(SURFACE, type_map={1: "Ni", 2: "O"})
    assert set(atoms.get_chemical_symbols()) <= {"Ni", "O"}


def test_type_ordering_matches_the_potential_type_map(nio_structure):
    types = to_lammps_type_order(nio_structure, {1: "Ni", 2: "O"})
    assert set(types) == {1, 2}
    assert types == [1 if s == "Ni" else 2 for s in nio_structure.get_chemical_symbols()]


def test_an_element_with_no_lammps_type_is_refused(nio_structure):
    """The LAMMPS analogue of element-coverage validation."""
    with pytest.raises(ConfigError, match="no LAMMPS type"):
        to_lammps_type_order(nio_structure, {1: "Ni"})
