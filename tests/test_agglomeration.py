from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from nio_md_prep.agglomeration import (
    MoleculeInstance,
    _compact_molecule_centers,
    _head_bias_for_replica,
    _nearest_phosphonate_pairs,
    _potcar_elements,
    prepare_agglomeration,
)
from nio_md_prep.cli import (
    _agglomeration_completion_message,
    _agglomeration_xtb_pending_message,
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


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _potcar(path: Path, elements: list[str]) -> Path:
    path.write_text(
        "\n".join(f"   TITEL  = PAW_PBE {element} 01Jan2001" for element in elements)
        + "\n",
        encoding="utf-8",
    )
    return path


def _molecular_poscar(
    path: Path,
    *,
    coordinate_offset: float = 0.0,
    reverse_within_elements: bool = False,
) -> Path:
    return _molecular_poscar_for_slug(
        path,
        "me-4pacz",
        coordinate_offset=coordinate_offset,
        reverse_within_elements=reverse_within_elements,
    )


def _molecular_poscar_for_slug(
    path: Path,
    slug: str,
    *,
    coordinate_offset: float = 0.0,
    reverse_within_elements: bool = False,
) -> Path:
    data = parse(ROOT / f"inputs/molecules/{slug}/ligpargen.lmp")
    symbols = elements(data)
    coordinates = atom_coordinates(data)
    element_order = sorted(set(symbols))
    rows = []
    for element in element_order:
        indexed = list(enumerate(zip(symbols, coordinates)))
        if reverse_within_elements:
            indexed.reverse()
        for index, (symbol, coordinate) in indexed:
            if symbol != element:
                continue
            x, y, z = coordinate
            if index == 0:
                x += coordinate_offset
            rows.append(f"{x + 25:.10f} {y + 25:.10f} {z + 25:.10f}")
    counts = [symbols.count(element) for element in element_order]
    path.write_text(
        f"{slug} molecular reference\n"
        "1.0\n"
        "50 0 0\n0 50 0\n0 0 50\n"
        + " ".join(element_order)
        + "\n"
        + " ".join(str(count) for count in counts)
        + "\nCartesian\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )
    return path


def _mixed_config(path: Path) -> Path:
    path.write_text(
        """[agglomeration]
replicas = 1
base_seed = 22345
radius_angstrom = 20.0
packmol_tolerance_angstrom = 2.0
vacuum_angstrom = 10.0
minimum_distance_angstrom = 1.0
center_scales = [1.0]

[[molecules]]
slug = "me-4pacz"
count = 1

[[molecules]]
slug = "dcz-4p"
count = 1
""",
        encoding="utf-8",
    )
    return path


def _mixed_packed(path: Path) -> Path:
    rows = []
    for slug, shift in (("me-4pacz", -15.0), ("dcz-4p", 15.0)):
        data = parse(ROOT / f"inputs/molecules/{slug}/ligpargen.lmp")
        rows.extend(
            f"{symbol} {x + shift:.10f} {y:.10f} {z:.10f}"
            for symbol, (x, y, z) in zip(elements(data), atom_coordinates(data))
        )
    path.write_text(
        f"{len(rows)}\nordered mixed-molecule fixture\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def test_completion_message_reports_created_vasp_folders(tmp_path):
    manifest = {
        "case_root": "vasp_runs",
        "xtb": {"enabled": True},
        "vasp_md": {"schedules": ["300K", "300K_to_400K", "400K"]},
        "replicas": [
            {
                "structures": [
                    {"path": "vasp_runs/n02/r000_s00_1p000"},
                    {"path": "vasp_runs/n02/r001_s00_1p000"},
                ]
            }
        ],
    }
    message = _agglomeration_completion_message(
        manifest, tmp_path / "campaign", structures_only=False
    )
    assert f"created under {tmp_path / 'campaign/vasp_runs'}" in message
    assert "2 structure case(s)" in message
    assert "6 temperature-stage run folder(s)" in message
    assert "300K, 300K_to_400K, 400K" in message
    assert "submit_temperature_jobs.sh" in message
    assert "launch_vasp_runs.sh --dry-run" in message
    regenerated = _agglomeration_completion_message(
        manifest,
        tmp_path / "campaign",
        structures_only=False,
        regenerated_vasp=True,
    )
    assert "VASP folders regenerated under" in regenerated


def test_xtb_pending_message_points_to_audit_and_safe_launcher(tmp_path):
    manifest = {
        "replicas": [
            {
                "structures": [
                    {"path": "vasp_runs/finished"},
                    {"status": "XTB_REQUIRED"},
                    {"status": "XTB_REQUIRED"},
                ]
            }
        ]
    }
    message = _agglomeration_xtb_pending_message(manifest, tmp_path / "campaign")
    assert "1 completed case(s), 2 pending case(s)" in message
    assert str(tmp_path / "campaign/audit_xtb_runs.py") in message
    assert str(tmp_path / "campaign/launch_xtb.sh") in message
    assert "--dry-run" in message


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


def test_phosphonate_pairing_never_reuses_a_head_or_pairs_one_molecule():
    instances = [
        MoleculeInstance("toy", 1, 0, 2, (1.0, 1.0), (0, 1), (1, 0)),
        MoleculeInstance("toy", 2, 2, 4, (1.0, 1.0), (2, 3), (3, 2)),
        MoleculeInstance("toy", 3, 4, 5, (1.0,), (4,), (4,)),
    ]
    pairs = _nearest_phosphonate_pairs(
        [(0.0, 0.0, 0.0), (50.0, 0.0, 0.0), (1.0, 0.0, 0.0),
         (51.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
        instances,
    )
    used = [index for pair in pairs for index in (
        pair["left_atom_index"], pair["right_atom_index"]
    )]
    assert len(used) == len(set(used))
    assert all(
        pair["left_molecule_id"] != pair["right_molecule_id"]
        for pair in pairs
    )
    assert len(pairs) == 2
    assert all("left_carbon_atom_index" in pair for pair in pairs)
    assert all("right_carbon_atom_index" in pair for pair in pairs)


def test_only_requested_replicas_receive_head_steering():
    settings = {"head_bias_replicas_per_family": 1}
    assert [_head_bias_for_replica(index, 3, settings) for index in range(3)] == [
        True, False, False
    ]
    assert not _head_bias_for_replica(0, 3, {"head_bias_enabled": False})


def test_agglomeration_builds_fixed_cell_center_scaled_vasp_cases(tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
    (reference / "KPOINTS").write_text("Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8")
    launcher = reference / "runvasp.sh"
    launcher.write_text("#!/bin/bash\nvasp_std\n", encoding="utf-8")
    launcher.chmod(0o755)
    _molecular_poscar(reference / "POSCAR")
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


def test_reference_poscar_is_the_authoritative_molecular_geometry(
    tmp_path, monkeypatch
):
    monkeypatch.setattr("nio_md_prep.agglomeration.shutil.which", lambda _: None)
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
    (reference / "KPOINTS").write_text(
        "Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8"
    )
    _molecular_poscar(
        reference / "POSCAR",
        coordinate_offset=0.1,
        reverse_within_elements=True,
    )
    _potcar(reference / "POTCAR", ["C", "H", "N", "O", "P"])

    output = tmp_path / "campaign"
    manifest_path = prepare_agglomeration(
        _config(tmp_path / "config.toml"),
        output,
        reference_dir=reference,
    )
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "PACKMOL_REQUIRED"
    assert manifest["molecular_geometry"]["source"] == "reference_poscar"
    assert manifest["molecular_geometry"]["sha256"]
    assert manifest["molecular_geometry"]["mapping"].startswith(
        "bond-graph isomorphism"
    )

    written = _read_ordered_xyz(output / "templates/me-4pacz.xyz")
    original = atom_coordinates(
        parse(ROOT / "inputs/molecules/me-4pacz/ligpargen.lmp")
    )
    modified = list(original)
    modified[0] = (original[0][0] + 0.1, original[0][1], original[0][2])

    def distance_inventory(coordinates):
        return sorted(
            _distance(coordinates[left], coordinates[right])
            for left in range(len(coordinates))
            for right in range(left + 1, len(coordinates))
        )

    written_coordinates = [coordinate for _, coordinate in written]
    assert distance_inventory(written_coordinates) == pytest.approx(
        distance_inventory(modified), abs=1.0e-8
    )
    assert distance_inventory(written_coordinates) != pytest.approx(
        distance_inventory(original), abs=1.0e-8
    )


def test_mixed_agglomeration_uses_two_slug_qualified_vasp_references(tmp_path):
    references = {}
    for slug, encut in (("me-4pacz", 400), ("dcz-4p", 520)):
        reference = tmp_path / slug
        reference.mkdir()
        (reference / "INCAR").write_text(f"ENCUT = {encut}\n", encoding="utf-8")
        (reference / "KPOINTS").write_text(
            "Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8"
        )
        _molecular_poscar_for_slug(reference / "POSCAR", slug)
        _potcar(reference / "POTCAR", ["C", "H", "N", "O", "P"])
        references[slug] = reference

    output = tmp_path / "mixed"
    manifest_path = prepare_agglomeration(
        _mixed_config(tmp_path / "mixed.toml"),
        output,
        reference_dirs=references,
        vasp_template_slug="dcz-4p",
        packed_xyz=_mixed_packed(tmp_path / "mixed.xyz"),
    )
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "COMPLETE"
    assert manifest["molecular_geometry"]["source"] == "reference_poscars"
    assert set(manifest["molecular_geometry"]["by_slug"]) == {
        "me-4pacz", "dcz-4p"
    }
    sanity = manifest["vasp_reference_sanity"]
    assert sanity["mode"] == "mixed_species"
    assert sanity["calculation_template_slug"] == "dcz-4p"
    assert sanity["combined_potcar_elements"] == ["C", "H", "N", "O", "P"]

    case = output / "vasp_runs/r000_s00_1p000"
    assert (case / "INCAR").read_text() == "ENCUT = 520\n"
    assert _potcar_elements(case / "POTCAR") == ["C", "H", "N", "O", "P"]
    case_manifest = json.loads((case / "agglomeration_manifest.json").read_text())
    assert {row["molecule_slug"] for row in case_manifest["atom_map"]} == {
        "me-4pacz", "dcz-4p"
    }
    assert len(case_manifest["atom_map"]) == 121


def test_mixed_agglomeration_rejects_conflicting_potcar_datasets(tmp_path):
    references = {}
    for slug in ("me-4pacz", "dcz-4p"):
        reference = tmp_path / slug
        reference.mkdir()
        (reference / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
        (reference / "KPOINTS").write_text(
            "Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8"
        )
        _molecular_poscar_for_slug(reference / "POSCAR", slug)
        _potcar(reference / "POTCAR", ["C", "H", "N", "O", "P"])
        references[slug] = reference
    references["dcz-4p"].joinpath("POTCAR").write_text(
        references["dcz-4p"].joinpath("POTCAR").read_text().replace(
            "PAW_PBE C", "PAW_PBE C_hard"
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="conflicting POTCAR datasets for C"):
        prepare_agglomeration(
            _mixed_config(tmp_path / "mixed.toml"),
            tmp_path / "mixed",
            reference_dirs=references,
            packed_xyz=_mixed_packed(tmp_path / "mixed.xyz"),
        )


def test_mixed_agglomeration_resumes_from_xtb_into_vasp_runs(tmp_path):
    config = _mixed_config(tmp_path / "mixed-xtb.toml")
    config.write_text(
        config.read_text(encoding="utf-8")
        + (
            "\n[xtb]\nenabled = true\nhead_bias_enabled = false\n"
            "maximum_nearest_molecule_distance_angstrom = 40.0\n"
        ),
        encoding="utf-8",
    )
    references = {}
    for slug, encut in (("me-4pacz", 400), ("dcz-4p", 520)):
        reference = tmp_path / slug
        reference.mkdir()
        (reference / "INCAR").write_text(f"ENCUT = {encut}\n", encoding="utf-8")
        (reference / "KPOINTS").write_text(
            "Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8"
        )
        _molecular_poscar_for_slug(reference / "POSCAR", slug)
        _potcar(reference / "POTCAR", ["C", "H", "N", "O", "P"])
        references[slug] = reference

    output = tmp_path / "mixed"
    first = prepare_agglomeration(
        config,
        output,
        reference_dirs=references,
        vasp_template_slug="dcz-4p",
        packed_xyz=_mixed_packed(tmp_path / "mixed.xyz"),
    )
    assert json.loads(first.read_text())["status"] == "XTB_REQUIRED"
    xtb_work = output / "xtb/r000_s00_1p000"
    assert int((xtb_work / "input.xyz").read_text().splitlines()[0]) == 121
    shutil.copy2(xtb_work / "input.xyz", xtb_work / "xtbfinal.xyz")
    original_xtb_input = (xtb_work / "input.xyz").read_text(encoding="utf-8")
    original_xtb_protocol = (xtb_work / "xtb_protocol.json").read_text(
        encoding="utf-8"
    )
    assert not (xtb_work / "xtb_protocol.complete").exists()
    config.write_text(
        config.read_text(encoding="utf-8") + "\n# harmless post-run comment\n",
        encoding="utf-8",
    )
    assert json.loads(first.read_text())["config_sha256"] != _file_sha256(config)

    incomplete = prepare_agglomeration(
        config,
        output,
        reference_dirs=references,
        vasp_template_slug="dcz-4p",
        packed_xyz=_mixed_packed(tmp_path / "mixed-incomplete.xyz"),
        regenerate_vasp=True,
    )
    assert json.loads(incomplete.read_text())["status"] == "XTB_REQUIRED"
    assert not (output / "launch_vasp_runs.sh").exists()
    regenerated_xtb_launch = prepare_agglomeration(
        config,
        output,
        regenerate_xtb_launcher=True,
    )
    regenerated_xtb_manifest = json.loads(regenerated_xtb_launch.read_text())
    assert regenerated_xtb_manifest["status"] == "XTB_REQUIRED"
    assert regenerated_xtb_manifest["xtb_launcher_regenerated"] is True
    assert regenerated_xtb_manifest["xtb_launcher"] == "launch_xtb.sh"
    assert regenerated_xtb_manifest["xtb_pool_launcher"] == "run_xtb_pool.sbatch"
    assert (xtb_work / "input.xyz").read_text(encoding="utf-8") == original_xtb_input
    assert (
        xtb_work / "xtb_protocol.json"
    ).read_text(encoding="utf-8") == original_xtb_protocol
    assert not (output / "launch_vasp_runs.sh").exists()

    completed_protocol = xtb_work / "xtb_protocol.complete"
    completed_protocol.write_text("legacy-protocol-hash\n", encoding="utf-8")
    legacy_audit = subprocess.run(
        [str(output / "audit_xtb_runs.py")],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "HASH_MISMATCH=1" in legacy_audit.stdout
    assert "--regenerate-vasp" in legacy_audit.stdout
    second = prepare_agglomeration(
        config,
        output,
        reference_dirs=references,
        vasp_template_slug="dcz-4p",
        packed_xyz=_mixed_packed(tmp_path / "mixed-again.xyz"),
        regenerate_vasp=True,
    )
    assert json.loads(second.read_text())["status"] == "COMPLETE"
    case = output / "vasp_runs/r000_s00_1p000"
    assert (case / "300K/INCAR").read_text().startswith("ENCUT = 520\n")
    assert sum(map(int, (case / "POSCAR").read_text().splitlines()[6].split())) == 121
    campaign_manifest = json.loads(second.read_text())
    assert campaign_manifest["vasp_launcher"] == "launch_vasp_runs.sh"
    assert (xtb_work / "input.xyz").read_text(encoding="utf-8") == original_xtb_input
    assert (
        xtb_work / "xtb_protocol.json"
    ).read_text(encoding="utf-8") == original_xtb_protocol
    case_manifest = json.loads((case / "agglomeration_manifest.json").read_text())
    assert case_manifest["xtb"]["protocol"][
        "regeneration_reused_completed_xtb"
    ] is True
    assert case_manifest["xtb"]["protocol"][
        "legacy_completion_marker_present"
    ] is True
    assert case_manifest["xtb"]["protocol"][
        "legacy_completion_marker_migrated"
    ] is True
    assert completed_protocol.read_text(encoding="utf-8") == (
        xtb_work / "xtb_protocol.sha256"
    ).read_text(encoding="utf-8")
    migrated_audit = subprocess.run(
        [str(output / "audit_xtb_runs.py")],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "xTB progress: 1/1 safely complete (100.0%); 0 pending" in (
        migrated_audit.stdout
    )
    campaign_launcher = output / "launch_vasp_runs.sh"
    assert campaign_launcher.stat().st_mode & 0o111
    subprocess.run(["bash", "-n", str(campaign_launcher)], check=True)
    dry_run = subprocess.run(
        [str(campaign_launcher), "--dry-run"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "DRY RUN: no jobs will be submitted" in dry_run.stdout
    assert dry_run.stdout.count("sbatch --parsable") == 3
    cancelled = subprocess.run(
        [str(campaign_launcher)],
        input="n\n",
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Submission cancelled; no jobs were launched." in cancelled.stdout
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()
    fake_sbatch = fake_bin / "sbatch"
    fake_sbatch.write_text(
        "#!/bin/bash\nprintf '%s\\n' \"$*\" >> \"$SBATCH_LOG\"\nprintf '12345\\n'\n",
        encoding="utf-8",
    )
    fake_sbatch.chmod(0o755)
    sbatch_log = tmp_path / "sbatch.log"
    launched = subprocess.run(
        [str(campaign_launcher), "--yes"],
        check=True,
        capture_output=True,
        text=True,
        env=os.environ
        | {
            "PATH": f"{fake_bin}:{os.environ['PATH']}",
            "SBATCH_LOG": str(sbatch_log),
        },
    )
    assert "All VASP stage jobs were submitted." in launched.stdout
    assert len(sbatch_log.read_text(encoding="utf-8").splitlines()) == 3

    campaign_launcher.unlink()
    (case / "300K/INCAR").write_text("STALE INPUT\n", encoding="utf-8")
    regenerated = prepare_agglomeration(
        config,
        output,
        reference_dirs=references,
        vasp_template_slug="dcz-4p",
        regenerate_vasp=True,
    )
    regenerated_manifest = json.loads(regenerated.read_text())
    assert regenerated_manifest["status"] == "COMPLETE"
    assert regenerated_manifest["regenerated_vasp"] is True
    assert campaign_launcher.is_file()
    assert (case / "300K/INCAR").read_text().startswith("ENCUT = 520\n")

    (case / "300K/OUTCAR").write_text("started VASP run\n", encoding="utf-8")
    with pytest.raises(ValueError, match="VASP runtime outputs.*OUTCAR"):
        prepare_agglomeration(
            config,
            output,
            reference_dirs=references,
            vasp_template_slug="dcz-4p",
            regenerate_vasp=True,
        )


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
    _molecular_poscar(reference / "POSCAR")
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
    assert "#SBATCH -p workq" in xtb_launcher
    assert "#SBATCH --cpus-per-task=8" in xtb_launcher
    assert f"root={output.resolve()}" in xtb_launcher
    assert "SLURM_SUBMIT_DIR" not in xtb_launcher
    assert 'dirname "$0"' not in xtb_launcher
    assert "XTB_ENV:=/project/lgutsev/env/xtb_env" in xtb_launcher
    assert 'OMP_NUM_THREADS="${threads},1"' in xtb_launcher
    assert 'MKL_NUM_THREADS="$threads"' in xtb_launcher
    assert "OMP_MAX_ACTIVE_LEVELS=1" in xtb_launcher
    assert "OMP_STACKSIZE=${XTB_OMP_STACKSIZE:-1G}" in xtb_launcher
    assert "--md --input md.inp" in xtb_launcher
    assert "--md --input md_steered.inp" in xtb_launcher
    assert "xtbopt_initial.xyz" in xtb_launcher
    assert "xtbsteered_last.xyz" in xtb_launcher
    assert "xtbmd_unbiased_last.xyz" in xtb_launcher
    assert "xtbfinal.xyz" in xtb_launcher
    assert "xtb_protocol.complete" in xtb_launcher
    assert "xtb_stages.protocol" in xtb_launcher
    assert "xtb_input.complete" in xtb_launcher
    assert xtb_launcher.index("input_changed=false") < xtb_launcher.index(
        "if [[ -s xtbopt.xyz && ! -s xtbopt_initial.xyz ]]"
    )
    xtb_audit = output / "audit_xtb_runs.py"
    xtb_submitter = output / "launch_xtb.sh"
    xtb_pool = output / "run_xtb_pool.sbatch"
    assert xtb_audit.stat().st_mode & 0o111
    assert xtb_submitter.stat().st_mode & 0o111
    assert xtb_pool.stat().st_mode & 0o111
    subprocess.run([sys.executable, "-m", "py_compile", str(xtb_audit)], check=True)
    subprocess.run(["bash", "-n", str(xtb_submitter)], check=True)
    subprocess.run(["bash", "-n", str(xtb_pool)], check=True)
    pool_text = xtb_pool.read_text(encoding="utf-8")
    assert "#SBATCH --ntasks-per-node=8" in pool_text
    assert "#SBATCH --cpus-per-task=8" in pool_text
    assert "#SBATCH --exclusive" in pool_text
    assert "srun --exclusive" in pool_text
    assert 'SLURM_ARRAY_TASK_ID="$index"' in pool_text
    audit_result = subprocess.run(
        [str(xtb_audit)], check=True, capture_output=True, text=True
    )
    assert "xTB progress: 0/1 safely complete (0.0%); 1 pending" in audit_result.stdout
    assert "Resources: one 64-core node carousel" in audit_result.stdout
    assert "8 simultaneous worker(s) x 8 OpenMP cores" in audit_result.stdout
    assert "INITIAL_OPT_COMPLETE=0" not in audit_result.stdout
    assert "[PENDING] xtb/r000_s00_1p000" in audit_result.stdout
    pending_indices = subprocess.run(
        [str(xtb_audit), "--pending-indices"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert pending_indices.stdout.strip() == "0"
    case_list = output / "xtb_cases.tsv"
    original_case_list = case_list.read_text(encoding="utf-8")
    synthetic_complete = output / "xtb/synthetic-complete"
    synthetic_complete.mkdir()
    (synthetic_complete / "xtbfinal.xyz").write_text("done\n", encoding="utf-8")
    (synthetic_complete / "xtb_protocol.sha256").write_text(
        "current\n", encoding="utf-8"
    )
    (synthetic_complete / "xtb_protocol.complete").write_text(
        "current\n", encoding="utf-8"
    )
    case_list.write_text(
        original_case_list
        + "xtb/synthetic-complete\t0\t0\n"
        + "xtb/synthetic-pending-a\t0\t0\n"
        + "xtb/synthetic-pending-b\t0\t0\n",
        encoding="utf-8",
    )
    compact_pending = subprocess.run(
        [str(xtb_audit), "--pending-indices"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert compact_pending.stdout.strip() == "0,2-3"
    case_list.write_text(original_case_list, encoding="utf-8")
    assert json.loads((output / "xtb_progress.json").read_text())["pending"] == 1
    assert "array_index,case,status,charge,uhf,detail" in (
        output / "xtb_progress.csv"
    ).read_text()
    submit_dry_run = subprocess.run(
        [str(xtb_submitter), "--dry-run"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "1 xTB case(s) remain pending" in submit_dry_run.stdout
    assert "Selected pending case indices: 0" in submit_dry_run.stdout
    assert "workq carousel: 8 workers x 8 cores = 64 cores" in (
        submit_dry_run.stdout
    )
    assert "run_xtb_pool.sbatch" in submit_dry_run.stdout
    assert "run_xtb_array.sbatch" not in submit_dry_run.stdout
    assert "DRY RUN: no jobs will be submitted" in submit_dry_run.stdout
    fake_srun_bin = tmp_path / "fake-srun-bin"
    fake_srun_bin.mkdir()
    fake_srun = fake_srun_bin / "srun"
    fake_srun.write_text(
        "#!/bin/bash\nprintf '%s\\n' \"$*\" >> \"$SRUN_LOG\"\n",
        encoding="utf-8",
    )
    fake_srun.chmod(0o755)
    srun_log = tmp_path / "srun.log"
    carousel = subprocess.run(
        [str(xtb_pool)],
        check=True,
        capture_output=True,
        text=True,
        env=os.environ
        | {
            "PATH": f"{fake_srun_bin}:{os.environ['PATH']}",
            "SRUN_LOG": str(srun_log),
            "XTB_ARRAY_INDICES": "0,2-3",
            "SLURM_JOB_ID": "123456",
            "TMPDIR": str(tmp_path),
        },
    )
    assert "3 selected case(s)" in carousel.stdout
    assert "8 simultaneous worker(s), 8 cores per worker" in carousel.stdout
    srun_lines = srun_log.read_text(encoding="utf-8").splitlines()
    assert len(srun_lines) == 3
    assert all("--cpus-per-task=8" in line for line in srun_lines)
    assert {
        next(field for field in line.split() if field.startswith("SLURM_ARRAY_TASK_ID="))
        for line in srun_lines
    } == {
        "SLURM_ARRAY_TASK_ID=0",
        "SLURM_ARRAY_TASK_ID=2",
        "SLURM_ARRAY_TASK_ID=3",
    }
    cancelled_xtb = subprocess.run(
        [str(xtb_submitter)],
        input="n\n",
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Submission cancelled; no jobs were launched" in cancelled_xtb.stdout
    fake_xtb_bin = tmp_path / "fake-xtb-bin"
    fake_xtb_bin.mkdir()
    fake_xtb_sbatch = fake_xtb_bin / "sbatch"
    fake_xtb_sbatch.write_text(
        "#!/bin/bash\nprintf '%s\\n' \"$*\" >> \"$SBATCH_LOG\"\nprintf '98765\\n'\n",
        encoding="utf-8",
    )
    fake_xtb_sbatch.chmod(0o755)
    xtb_sbatch_log = tmp_path / "xtb-sbatch.log"
    submitted_xtb = subprocess.run(
        [str(xtb_submitter), "--yes"],
        check=True,
        capture_output=True,
        text=True,
        env=os.environ
        | {
            "PATH": f"{fake_xtb_bin}:{os.environ['PATH']}",
            "SBATCH_LOG": str(xtb_sbatch_log),
        },
    )
    assert "xTB work submitted" in submitted_xtb.stdout
    assert xtb_sbatch_log.read_text(encoding="utf-8").strip() == (
        "--export=ALL,XTB_ARRAY_INDICES=0 run_xtb_pool.sbatch"
    )
    assert not (output / "vasp_runs/r000_s00_1p000/POSCAR").exists()
    xtb_work = output / "xtb/r000_s00_1p000"
    md_input = (xtb_work / "md.inp").read_text()
    assert "$seed 117074" in md_input
    assert "temp=450" in md_input
    assert "time=14" in md_input
    assert (xtb_work / "steering_plan.json").is_file()
    subprocess.run(
        [
            sys.executable,
            str(output / "write_xtb_steering_input.py"),
            str(xtb_work / "input.xyz"),
            str(xtb_work / "steering_plan.json"),
            str(xtb_work / "md_steered.inp"),
            str(xtb_work / "steering_restraints.json"),
        ],
        check=True,
    )
    steered_input = (xtb_work / "md_steered.inp").read_text()
    assert "$seed 12345" in steered_input
    assert "force constant=0.005" in steered_input
    assert "angle: 18,19,64,170" in steered_input
    assert "angle: 63,64,19,170" in steered_input
    restraints = json.loads((xtb_work / "steering_restraints.json").read_text())
    restraint = restraints["pairs"][0][
        "restraint_distance_angstrom"
    ]
    assert f"distance: 19,64,{restraint:.8f}" in steered_input
    assert restraint <= restraints["pairs"][0][
        "post_optimization_distance_angstrom"
    ]
    assert "time=6" in steered_input
    first_frame = (xtb_work / "input.xyz").read_text()
    second_frame = first_frame.replace(
        "Packmol seed 12345", "final synthetic MD frame", 1
    )
    trajectory = xtb_work / "xtb.trj"
    trajectory.write_text(first_frame + second_frame, encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(output / "extract_last_xtb_frame.py"),
            str(trajectory),
            str(xtb_work / "xtbmd_last.xyz"),
        ],
        check=True,
    )
    assert "final synthetic MD frame" in (
        xtb_work / "xtbmd_last.xyz"
    ).read_text()
    shutil.copy2(xtb_work / "input.xyz", xtb_work / "xtbfinal.xyz")
    shutil.copy2(
        xtb_work / "xtb_protocol.sha256",
        xtb_work / "xtb_protocol.complete",
    )
    completed_audit = subprocess.run(
        [str(xtb_audit)], check=True, capture_output=True, text=True
    )
    assert "xTB progress: 1/1 safely complete (100.0%); 0 pending" in (
        completed_audit.stdout
    )
    no_submit = subprocess.run(
        [str(xtb_submitter), "--yes"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "All xTB cases are safely complete; no work was submitted" in (
        no_submit.stdout
    )

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
    case_manifest = json.loads((case / "agglomeration_manifest.json").read_text())
    assert case_manifest["xtb"]["vasp_source_xyz"] == "xtbfinal.xyz"
    assert case_manifest["xtb"]["protocol"]["total_md_time_ps"] == 20.0
    assert case_manifest["xtb"]["protocol"]["input_xyz_sha256"] == _file_sha256(
        xtb_work / "input.xyz"
    )
    assert case_manifest["xtb"]["protocol"]["head_bias"]["enabled"] is True
    assert case_manifest["xtb"]["protocol"]["head_bias"]["unbiased_time_ps"] == 14.0
    assert case_manifest["xtb"]["protocol"]["head_bias"][
        "target_axis_angle_degrees"
    ] == 170.0


def test_xtb_input_geometry_change_invalidates_completed_protocol(tmp_path):
    config = tmp_path / "xtb-unbiased.toml"
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
head_bias_enabled = false

[[molecules]]
slug = "me-4pacz"
count = 2
""",
        encoding="utf-8",
    )
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "INCAR").write_text("ENCUT = 520\n", encoding="utf-8")
    (reference / "KPOINTS").write_text(
        "Gamma\n0\nGamma\n1 1 1\n0 0 0\n", encoding="utf-8"
    )
    _molecular_poscar(reference / "POSCAR")
    _potcar(reference / "POTCAR", ["C", "H", "N", "O", "P"])

    first_packed = _packed(tmp_path / "packed-first.xyz")
    output = tmp_path / "campaign"
    first = prepare_agglomeration(
        config, output, reference_dir=reference, packed_xyz=first_packed
    )
    assert json.loads(first.read_text())["status"] == "XTB_REQUIRED"
    xtb_work = output / "xtb/r000_s00_1p000"
    first_protocol = (xtb_work / "xtb_protocol.sha256").read_text()
    shutil.copy2(xtb_work / "input.xyz", xtb_work / "xtbfinal.xyz")
    shutil.copy2(
        xtb_work / "xtb_protocol.sha256", xtb_work / "xtb_protocol.complete"
    )

    changed_packed = _packed(tmp_path / "packed-changed.xyz")
    lines = changed_packed.read_text(encoding="utf-8").splitlines()
    atoms_per_molecule = parse(
        ROOT / "inputs/molecules/me-4pacz/ligpargen.lmp"
    ).count("Atoms")
    for index in range(2, 2 + atoms_per_molecule):
        fields = lines[index].split()
        fields[1] = f"{float(fields[1]) + 1.0:.10f}"
        lines[index] = " ".join(fields)
    changed_packed.write_text("\n".join(lines) + "\n", encoding="utf-8")

    second = prepare_agglomeration(
        config, output, reference_dir=reference, packed_xyz=changed_packed
    )
    assert json.loads(second.read_text())["status"] == "XTB_REQUIRED"
    assert (xtb_work / "xtb_protocol.sha256").read_text() != first_protocol
    assert (xtb_work / "xtb_protocol.complete").read_text() == first_protocol


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
    _molecular_poscar(reference / "POSCAR")
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
    _molecular_poscar(reference / "POSCAR")
    _potcar(reference / "POTCAR", ["H", "C", "N", "O", "P"])
    with pytest.raises(ValueError, match="POTCAR element order.*does not match"):
        prepare_agglomeration(
            _config(tmp_path / "config.toml"),
            tmp_path / "output",
            reference_dir=reference,
            packed_xyz=_packed(tmp_path / "packed.xyz"),
        )
