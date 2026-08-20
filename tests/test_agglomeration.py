from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from nio_md_prep.agglomeration import prepare_agglomeration
from nio_md_prep.geometry import elements
from nio_md_prep.lammps import atom_coordinates, parse

ROOT = Path(__file__).parents[1]


def _config(path: Path) -> Path:
    path.write_text(
        """[agglomeration]
replicas = 1
base_seed = 12345
radius_angstrom = 15.0
packmol_tolerance_angstrom = 2.0
vacuum_angstrom = 10.0
minimum_distance_angstrom = 1.0
center_scales = [1.0, 1.5]

[[molecules]]
slug = "me-4pacz"
count = 2
""",
        encoding="utf-8",
    )
    return path


def _packed(path: Path) -> Path:
    template = parse(ROOT / "inputs/molecules/me-4pacz/ligpargen.lmp")
    symbols = elements(template)
    coordinates = atom_coordinates(template)
    rows = []
    for shift in (-12.0, 12.0):
        rows.extend(
            f"{symbol} {x + shift:.10f} {y:.10f} {z:.10f}"
            for symbol, (x, y, z) in zip(symbols, coordinates)
        )
    path.write_text(
        f"{len(rows)}\nordered two-molecule fixture\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def _read_ordered_xyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()[2:]
    return [
        (fields[0], tuple(float(value) for value in fields[1:4]))
        for fields in (line.split() for line in lines)
    ]


def _molecule_centers(case: Path):
    atoms = _read_ordered_xyz(case / "structure.xyz")
    manifest = json.loads((case / "agglomeration_manifest.json").read_text())
    grouped = {}
    for (_, coordinate), identity in zip(atoms, manifest["atom_map"]):
        grouped.setdefault(identity["molecule_id"], []).append(coordinate)
    return {
        molecule_id: tuple(sum(point[axis] for point in points) / len(points) for axis in range(3))
        for molecule_id, points in grouped.items()
    }


def _distance(first, second):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def _potcar(path: Path, elements: list[str]) -> Path:
    path.write_text(
        "\n".join(f"   TITEL  = PAW_PBE {element} 01Jan2001" for element in elements)
        + "\n",
        encoding="utf-8",
    )
    return path


def test_agglomeration_builds_fixed_cell_center_scaled_vasp_cases(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
    (reference / "KPOINTS").write_text("Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8")
    launcher = reference / "runvasp.sh"
    launcher.write_text("#!/bin/bash\nvasp_std\n", encoding="utf-8")
    launcher.chmod(0o755)
    (reference / "POSCAR").write_text("must not be copied\n", encoding="utf-8")
    _potcar(reference / "POTCAR", ["C", "H", "N", "O", "P"])
    output = tmp_path / "agglomeration"
    manifest_path = prepare_agglomeration(
        _config(tmp_path / "config.toml"),
        output,
        reference_dir=reference,
        packed_xyz=_packed(tmp_path / "packed.xyz"),
    )
    root_manifest = json.loads(manifest_path.read_text())
    assert root_manifest["status"] == "COMPLETE"
    assert root_manifest["sanity_status"] == "PASS"
    assert root_manifest["vasp_reference_sanity"]["potcar_status"] == "PASS"
    assert root_manifest["output_mode"] == "vasp_runs"
    cases = sorted((output / "vasp_runs").iterdir())
    assert len(cases) == 2
    assert (cases[0] / "INCAR").read_text() == "ENCUT = 520\n"
    assert (cases[0] / "KPOINTS").is_file()
    assert (cases[0] / "POTCAR").is_file()
    assert (cases[0] / "runvasp.sh").is_file()
    assert (cases[0] / "runvasp.sh").stat().st_mode & 0o111
    assert not (cases[0] / "CONTCAR").exists()
    poscars = [(case / "POSCAR").read_text().splitlines() for case in cases]
    assert poscars[0][5].split() == sorted(poscars[0][5].split())
    assert poscars[0][2:5] == poscars[1][2:5]
    centers = [_molecule_centers(case) for case in cases]
    base = _distance(centers[0][1], centers[0][2])
    separated = _distance(centers[1][1], centers[1][2])
    assert separated == pytest.approx(1.5 * base, rel=1e-9)
    assert len(json.loads((cases[0] / "agglomeration_manifest.json").read_text())["atom_map"]) == 90
    assert len((output / "agglomeration_index.csv").read_text().splitlines()) == 3
    assert len((output / "agglomeration_sanity.csv").read_text().splitlines()) == 3
    sanity = json.loads((cases[0] / "sanity_report.json").read_text())
    assert sanity["status"] == "PASS"
    assert all(sanity["checks"].values())
    assert sanity["minimum_face_clearance_angstrom"] >= 10.0 - 1.0e-7
    assert sanity["periodic_image_separation_lower_bound_angstrom"] >= 20.0 - 1.0e-7
    assert sanity["maximum_intramolecular_distance_deviation_angstrom"] < 1.0e-7
    assert (cases[0] / "sanity_report.txt").read_text().startswith("PASS\n")


def test_agglomeration_packmol_inputs_can_be_resumed(tmp_path, monkeypatch):
    monkeypatch.setattr("nio_md_prep.agglomeration.shutil.which", lambda _: None)
    config = _config(tmp_path / "config.toml")
    output = tmp_path / "agglomeration"
    first = prepare_agglomeration(config, output, structures_only=True)
    assert json.loads(first.read_text())["status"] == "PACKMOL_REQUIRED"
    work = output / "packmol/replica_000"
    assert "inside sphere" in (work / "packmol.inp").read_text()
    _packed(work / "packed.xyz")
    second = prepare_agglomeration(config, output, structures_only=True)
    assert json.loads(second.read_text())["status"] == "COMPLETE"
    assert len(list((output / "structures").iterdir())) == 2


def test_agglomeration_requires_complete_vasp_reference_by_default(tmp_path):
    with pytest.raises(ValueError, match="complete VASP --reference-dir is required"):
        prepare_agglomeration(
            _config(tmp_path / "config.toml"),
            tmp_path / "output",
            packed_xyz=_packed(tmp_path / "packed.xyz"),
        )


def test_agglomeration_rejects_incomplete_vasp_reference(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing KPOINTS, POTCAR"):
        prepare_agglomeration(
            _config(tmp_path / "config.toml"),
            tmp_path / "output",
            reference_dir=reference,
            packed_xyz=_packed(tmp_path / "packed.xyz"),
        )


def test_agglomeration_rejects_a_single_molecule(tmp_path):
    config = tmp_path / "one.toml"
    config.write_text(
        """[agglomeration]
replicas = 1
[[molecules]]
slug = "me-4pacz"
count = 1
""",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="at least two molecules"):
        prepare_agglomeration(config, tmp_path / "output")


def test_agglomeration_rejects_wrong_potcar_order(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
    (reference / "KPOINTS").write_text("Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8")
    _potcar(reference / "POTCAR", ["H", "C", "N", "O", "P"])
    with pytest.raises(ValueError, match="POTCAR element order.*does not match"):
        prepare_agglomeration(
            _config(tmp_path / "config.toml"),
            tmp_path / "output",
            reference_dir=reference,
            packed_xyz=_packed(tmp_path / "packed.xyz"),
        )
