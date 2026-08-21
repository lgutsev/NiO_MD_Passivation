from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import pytest

from nio_md_prep.agglomeration import (
    MoleculeInstance,
    _compact_molecule_centers,
    prepare_agglomeration,
)
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


def test_adaptive_compaction_connects_an_isolated_molecule_without_overlap():
    instances = [
        MoleculeInstance("toy", index + 1, index, index + 1, (1.0,))
        for index in range(3)
    ]
    compacted, _ = _compact_molecule_centers(
        [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (20.0, 0.0, 0.0)],
        ["H", "H", "H"],
        instances,
        3.0,
    )
    distances = [
        min(
            _distance(compacted[index], compacted[other])
            for other in range(3)
            if other != index
        )
        for index in range(3)
    ]
    assert min(distances) >= 3.0 - 1.0e-7
    assert max(distances) <= 3.0 + 1.0e-7


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
    assert sanity["minimum_intermolecular_oxygen_oxygen_distance_angstrom"] > 2.5
    assert sanity["closest_intermolecular_pair"]["left"]["molecule_id"] != sanity[
        "closest_intermolecular_pair"
    ]["right"]["molecule_id"]
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


def test_xtb_stage_resumes_into_300_and_400_kelvin_vasp_runs(tmp_path):
    config = tmp_path / "xtb.toml"
    config.write_text(
        """[agglomeration]
replicas = 1
base_seed = 12345
radius_angstrom = 15.0
packmol_tolerance_angstrom = 3.0
vacuum_angstrom = 8.0
minimum_distance_angstrom = 2.5
compact_to_distance_angstrom = 3.0
center_scales = [1.0]

[xtb]
enabled = true
minimum_oxygen_oxygen_distance_angstrom = 2.2

[vasp_md]
heating_steps = 111
hold_steps = 222

[[molecules]]
slug = "me-4pacz"
count = 2
""",
        encoding="utf-8",
    )
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "INCAR").write_text(
        "ENCUT = 520\nNSW = 999\nTEBEG = 1\nTEEND = 2\n",
        encoding="utf-8",
    )
    (reference / "KPOINTS").write_text(
        "Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8"
    )
    launcher = reference / "runvasp.sh"
    launcher.write_text(
        "#!/bin/bash\n#SBATCH -J old\nSECONDS=0\nsrun -n128 vasp_gam\n",
        encoding="utf-8",
    )
    launcher.chmod(0o755)
    _potcar(reference / "POTCAR", ["C", "H", "N", "O", "P"])
    output = tmp_path / "campaign"
    first = prepare_agglomeration(
        config,
        output,
        reference_dir=reference,
        packed_xyz=_packed(tmp_path / "packed.xyz"),
    )
    assert json.loads(first.read_text())["status"] == "XTB_REQUIRED"
    xtb_launcher = (output / "run_xtb_array.sbatch").read_text()
    assert "#SBATCH -p single" in xtb_launcher
    assert "#SBATCH -c 1" in xtb_launcher
    assert "XTB_ENV:=/project/lgutsev/env/xtb_env" in xtb_launcher
    assert not (output / "vasp_runs/r000_s00_1p000/POSCAR").exists()
    xtb_work = output / "xtb/r000_s00_1p000"
    shutil.copy2(xtb_work / "input.xyz", xtb_work / "xtbopt.xyz")

    second = prepare_agglomeration(
        config,
        output,
        reference_dir=reference,
        packed_xyz=_packed(tmp_path / "packed-again.xyz"),
    )
    assert json.loads(second.read_text())["status"] == "COMPLETE"
    case = output / "vasp_runs/r000_s00_1p000"
    incar300 = (case / "300K/INCAR").read_text()
    assert "NSW = 222" in incar300
    assert "TEBEG = 300" in incar300 and "TEEND = 300" in incar300
    incar_heat = (case / "400K/01_heat_300_to_400K/INCAR").read_text()
    assert "NSW = 111" in incar_heat
    assert "TEBEG = 300" in incar_heat and "TEEND = 400" in incar_heat
    incar400 = (case / "400K/02_hold_400K/INCAR").read_text()
    assert "NSW = 222" in incar400
    assert "TEBEG = 400" in incar400 and "TEEND = 400" in incar400
    assert "cp ../01_heat_300_to_400K/CONTCAR POSCAR" in (
        case / "400K/02_hold_400K/runvasp.sh"
    ).read_text()
    assert (case / "submit_temperature_jobs.sh").stat().st_mode & 0o111


def test_agglomeration_uses_loose_cost_conscious_defaults(tmp_path):
    config = tmp_path / "defaults.toml"
    config.write_text(
        """[agglomeration]
replicas = 1

[[molecules]]
slug = "me-4pacz"
count = 2
""",
        encoding="utf-8",
    )
    manifest = prepare_agglomeration(
        config,
        tmp_path / "output",
        packed_xyz=_packed(tmp_path / "packed.xyz"),
        structures_only=True,
    )
    settings = json.loads(manifest.read_text())["settings"]
    assert settings["radius_angstrom"] == 14.0
    assert settings["packmol_tolerance_angstrom"] == 3.0
    assert settings["vacuum_angstrom"] == 8.0
    assert settings["minimum_distance_angstrom"] == 2.5
    assert settings["center_scales"] == [1.0]


def test_size_stratified_campaign_allocates_more_small_replicas(tmp_path, monkeypatch):
    monkeypatch.setattr("nio_md_prep.agglomeration.shutil.which", lambda _: None)
    output = tmp_path / "campaign"
    manifest_path = prepare_agglomeration(
        ROOT / "examples/agglomeration/me-4pacz.toml",
        output,
        structures_only=True,
    )
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "PACKMOL_REQUIRED"
    expected = {"n02": 3, "n03": 3, "n04": 3, "n06": 2, "n08": 1}
    assert {
        family["name"]: family["replicas"]
        for family in manifest["agglomerates"]
    } == expected
    assert len(manifest["replicas"]) == 12
    for name, replicas in expected.items():
        work = output / "packmol" / name
        assert len(list(work.glob("replica_*"))) == replicas
    assert "number 2" in (output / "packmol/n02/replica_000/packmol.inp").read_text()
    assert "number 8" in (output / "packmol/n08/replica_000/packmol.inp").read_text()


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
