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
    constraint_summary,
    describe_structure,
    fixed_atom_indices,
    from_lammps_data,
    is_partially_periodic,
    is_periodic,
    largest_empty_gaps,
    load_structure,
    outside_nonperiodic_faces,
    periodic_axes,
    structure_digest,
    to_lammps_type_order,
    vacuum_gaps,
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
    assert report["pbc"] == [True, True, True]
    assert report["constraints"] == {"fixed_atoms": [], "n_fixed": 0, "unsupported": []}
    assert len(report["largest_empty_gap_angstrom"]) == 3
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
    """A prepared classical build becomes a valid MLIP interchange structure.

    The classical slab decks use ``boundary p p f``; the data file cannot say
    so, and used to be read as fully periodic. Its box starts at zlo = -15,
    which must become the cell origin rather than being dropped.
    """
    from nio_md_prep.lammps import parse

    atoms = from_lammps_data(SURFACE, pbc=(True, True, False))
    data = parse(SURFACE)
    assert len(atoms) == data.count("Atoms")
    assert set(atoms.get_chemical_symbols()) == {"Ni", "O"}
    assert periodic_axes(atoms) == (True, True, False)
    assert not is_periodic(atoms) and is_partially_periodic(atoms)
    assert atoms.cell.lengths() == pytest.approx([125.1, 41.7, 285.0])
    first = data.sections["Atoms"][0].fields
    assert atoms.positions[0][2] == pytest.approx(float(first[6]) + 15.0)
    assert outside_nonperiodic_faces(atoms) == {}
    # Topology is deliberately dropped; it stays in the LAMMPS representation.
    assert not hasattr(atoms, "bonds")


@pytest.mark.skipif(not SURFACE.exists(), reason="the NiO surface fixture is absent")
def test_an_explicit_type_map_overrides_mass_inference():
    atoms = from_lammps_data(SURFACE, pbc=(True, True, False), type_map={1: "Ni", 2: "O"})
    assert set(atoms.get_chemical_symbols()) <= {"Ni", "O"}


@pytest.mark.skipif(not SURFACE.exists(), reason="the NiO surface fixture is absent")
def test_a_lammps_data_file_needs_an_explicit_pbc():
    with pytest.raises(ConfigError, match="explicit pbc"):
        from_lammps_data(SURFACE)
    with pytest.raises(ConfigError, match="explicit pbc"):
        load_structure(SURFACE, format="lammps-data-nio")


TRICLINIC_DATA = """triclinic test box

2 atoms
1 atom types

1.0 5.0 xlo xhi
-2.0 2.0 ylo yhi
0.5 6.5 zlo zhi
1.5 -0.5 0.25 xy xz yz

Masses

1 58.6934 # Ni

Atoms # full

1 1 1 0.0 1.0 -2.0 0.5
2 1 1 0.0 4.0 1.0 3.5
"""


def test_a_triclinic_box_keeps_its_tilt_and_origin(tmp_path):
    """The tilt used to be dropped (a 60-degree cell became 90 degrees)."""
    path = tmp_path / "tri.lmp"
    path.write_text(TRICLINIC_DATA, encoding="utf-8")
    atoms = from_lammps_data(path, pbc=(True, True, True))
    import numpy as np

    assert np.allclose(atoms.get_cell(), [[4.0, 0.0, 0.0], [1.5, 4.0, 0.0], [-0.5, 0.25, 6.0]])
    # Positions are relative to the (xlo, ylo, zlo) corner.
    assert np.allclose(atoms.positions, [[0.0, 0.0, 0.0], [3.0, 3.0, 3.0]])
    assert atoms.get_chemical_symbols() == ["Ni", "Ni"]


def test_an_atomic_style_data_file_is_refused_rather_than_misread(tmp_path):
    path = tmp_path / "atomic.lmp"
    path.write_text(
        TRICLINIC_DATA.replace("Atoms # full", "Atoms # atomic")
        .replace("1 1 1 0.0 1.0 -2.0 0.5", "1 1 1.0 -2.0 0.5")
        .replace("2 1 1 0.0 4.0 1.0 3.5", "2 1 4.0 1.0 3.5"),
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="atom_style full"):
        from_lammps_data(path, pbc=(True, True, True))


def test_a_declared_pbc_replaces_the_files_periodicity(tmp_path, nio_structure):
    """A POSCAR slab is always read fully periodic; structure.pbc says otherwise."""
    from ase.io import write

    from nio_md_prep.mlip.specs import StructureSpec

    path = tmp_path / "POSCAR"
    write(str(path), nio_structure, format="vasp")
    assert periodic_axes(load_structure(StructureSpec(path=path, format="vasp"))) == (
        True, True, True
    )
    slab = load_structure(StructureSpec(path=path, format="vasp", pbc=(True, False, True)))
    assert periodic_axes(slab) == (True, False, True)


def test_a_periodic_axis_without_a_cell_vector_is_refused(tmp_path):
    from ase import Atoms
    from ase.io import write

    path = tmp_path / "cluster.xyz"
    write(str(path), Atoms("Ni2", positions=[[0, 0, 0], [0, 0, 2.5]]), format="extxyz")
    with pytest.raises(ConfigError, match="zero length"):
        load_structure(path, pbc=(True, True, True))


# --- per-axis periodicity helpers -----------------------------------------


def nio_slab(vacuum: float = 7.5, axis: int = 2):
    from ase.build import bulk

    atoms = bulk("NiO", "rocksalt", a=4.17, cubic=True).repeat((2, 2, 2))
    atoms.center(vacuum=vacuum, axis=axis)
    return atoms


def test_bulk_has_no_vacuum_gap(nio_structure):
    gaps = largest_empty_gaps(nio_structure)
    # Rocksalt fcc-primitive cell: widest empty slab is one (111) spacing, ~2.4 A.
    assert all(g is not None and g < 3.0 for g in gaps)
    assert vacuum_gaps(nio_structure) == {}


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_a_slab_written_as_a_periodic_cell_is_recognised(axis):
    """Vacuum normal to x, y or z is found along that lattice direction only."""
    slab = nio_slab(vacuum=7.5, axis=axis)
    gaps = vacuum_gaps(slab, threshold_angstrom=5.0)
    assert set(gaps) == {axis}
    # center(vacuum=7.5) leaves 7.5 A beyond the outermost layer on each side.
    assert gaps[axis] == pytest.approx(15.0, abs=1e-9)
    assert vacuum_gaps(slab, threshold_angstrom=20.0) == {}


def test_gaps_are_measured_normal_to_the_lattice_plane():
    """A tilted c vector must not inflate the gap: distance is 1/|b_c|, not |c|."""
    slab = nio_slab(vacuum=7.5, axis=2)
    cell = slab.get_cell().array.copy()
    cell[2, 0] += 6.0  # tilt c along x; the planes of constant c-fraction do not move apart
    slab.set_cell(cell, scale_atoms=False)
    assert largest_empty_gaps(slab)[2] == pytest.approx(15.0, abs=1e-9)


def test_non_periodic_axes_report_no_gap_and_bounds_are_checked():
    slab = nio_slab(vacuum=7.5, axis=2)
    slab.pbc = (True, True, False)
    assert largest_empty_gaps(slab)[2] is None
    assert outside_nonperiodic_faces(slab) == {}
    slab.positions[3, 2] = -0.1
    slab.positions[5, 2] = slab.cell[2, 2]  # exactly on the upper face: outside [0, 1)
    assert outside_nonperiodic_faces(slab) == {2: (3, 5)}


def test_bounds_need_a_real_box():
    from ase import Atoms

    cluster = Atoms("Ni2", positions=[[0, 0, 0], [0, 0, 2.5]])
    with pytest.raises(ConfigError, match="rank-3 cell"):
        outside_nonperiodic_faces(cluster)


# --- constraints ------------------------------------------------------------


def test_fixed_atoms_are_extracted_sorted(nio_structure):
    from ase.constraints import FixAtoms

    atoms = nio_structure.copy()
    atoms.set_constraint([FixAtoms(indices=[5, 1]), FixAtoms(mask=[i == 3 for i in range(16)])])
    assert fixed_atom_indices(atoms) == (1, 3, 5)
    assert constraint_summary(atoms) == {
        "fixed_atoms": [1, 3, 5], "n_fixed": 3, "unsupported": []
    }


def test_any_other_constraint_is_refused_by_name(nio_structure):
    from ase.constraints import FixAtoms, FixBondLength

    atoms = nio_structure.copy()
    atoms.set_constraint([FixAtoms(indices=[0]), FixBondLength(0, 1)])
    # ASE's FixBondLength factory builds a FixBondLengths constraint.
    with pytest.raises(ConfigError, match="FixBondLengths"):
        fixed_atom_indices(atoms)
    # The description records it without raising.
    assert constraint_summary(atoms)["unsupported"] == ["FixBondLengths"]
    assert describe_structure(atoms)["constraints"]["n_fixed"] == 1


def test_selective_dynamics_from_a_poscar_become_fixed_atoms(tmp_path, nio_structure):
    from ase.constraints import FixAtoms
    from ase.io import write

    atoms = nio_structure.copy()
    atoms.set_constraint(FixAtoms(indices=[0, 2]))
    path = tmp_path / "POSCAR"
    write(str(path), atoms, format="vasp")
    assert fixed_atom_indices(load_structure(path, format="vasp")) == (0, 2)


# --- identity -------------------------------------------------------------


def test_digest_distinguishes_magnetic_initialisations(nio_structure):
    """FM and AFM NiO of the same geometry are different electronic states."""
    ferro = nio_structure.copy()
    ferro.set_initial_magnetic_moments([2.0 if s == "Ni" else 0.0 for s in ferro.symbols])
    anti = ferro.copy()
    moments = anti.get_initial_magnetic_moments()
    ni = [i for i, s in enumerate(anti.symbols) if s == "Ni"]
    moments[ni[::2]] *= -1
    anti.set_initial_magnetic_moments(moments)
    assert structure_digest(ferro) != structure_digest(anti)


def test_digest_treats_all_zero_arrays_as_absent(nio_structure):
    """ASE reads an absent array as zeros; the digest must agree."""
    bare = nio_structure.copy()
    for name in ("initial_magmoms", "initial_charges", "tags"):
        if bare.has(name):
            del bare.arrays[name]
    zeros = bare.copy()
    zeros.set_initial_magnetic_moments([0.0] * len(zeros))
    zeros.set_initial_charges([0.0] * len(zeros))
    zeros.set_tags([0] * len(zeros))
    assert structure_digest(zeros) == structure_digest(bare)


@pytest.mark.parametrize("change", ["charges", "tags", "fixed", "pbc"])
def test_digest_changes_with_charges_tags_constraints_and_pbc(nio_structure, change):
    from ase.constraints import FixAtoms

    changed = nio_structure.copy()
    if change == "charges":
        changed.set_initial_charges([0.5] + [0.0] * (len(changed) - 1))
    elif change == "tags":
        changed.set_tags([1] + [0] * (len(changed) - 1))
    elif change == "fixed":
        changed.set_constraint(FixAtoms(indices=[0]))
    else:
        changed.pbc = (True, True, False)
    assert structure_digest(changed) != structure_digest(nio_structure)


def test_digest_survives_an_extxyz_round_trip_with_magmoms_and_constraints(
    tmp_path, nio_structure
):
    from ase.constraints import FixAtoms
    from ase.io import read, write

    atoms = nio_structure.copy()
    atoms.set_initial_magnetic_moments([1.7 if s == "Ni" else 0.0 for s in atoms.symbols])
    atoms.set_constraint(FixAtoms(indices=[0, 1]))
    path = tmp_path / "cell.xyz"
    write(str(path), atoms, format="extxyz")
    assert structure_digest(read(str(path))) == structure_digest(atoms)


def test_type_ordering_matches_the_potential_type_map(nio_structure):
    types = to_lammps_type_order(nio_structure, {1: "Ni", 2: "O"})
    assert set(types) == {1, 2}
    assert types == [1 if s == "Ni" else 2 for s in nio_structure.get_chemical_symbols()]


def test_an_element_with_no_lammps_type_is_refused(nio_structure):
    """The LAMMPS analogue of element-coverage validation."""
    with pytest.raises(ConfigError, match="no LAMMPS type"):
        to_lammps_type_order(nio_structure, {1: "Ni"})
