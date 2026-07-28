from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pytest

from nio_md_prep.build import build
from nio_md_prep.geometry import elements
from nio_md_prep.lammps import Record, parse, write
from nio_md_prep.lego import (
    identify_x_tunnel,
    prepare_lego_continuation,
    prepare_lego_stage2,
)


ROOT = Path(__file__).parents[1]


def packed_fixture(tmp_path: Path, shift=(15.0, 10.0, 62.0)) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    template = parse(ROOT / "inputs/molecules/me-4pacz/ligpargen.lmp")
    rows = []
    for symbol, atom in zip(elements(template), template.sections["Atoms"]):
        rows.append(
            f"{symbol} "
            f"{float(atom.fields[4]) + shift[0]} "
            f"{float(atom.fields[5]) + shift[1]} "
            f"{float(atom.fields[6]) + shift[2]}"
        )
    path = tmp_path / "packed.xyz"
    path.write_text(
        f"{len(rows)}\nmock Packmol\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def primary_fixture(tmp_path: Path) -> Path:
    packed = packed_fixture(tmp_path / "primary-packed")
    stage1 = build(
        ROOT / "tests/data/small-study.toml",
        tmp_path / "stage1",
        packed_xyz=packed,
    )
    data = parse(stage1 / "topology_output.lmp")
    data.sections["Velocities"] = [
        Record([atom.fields[0], "1.0", "2.0", "3.0"])
        for atom in data.sections["Atoms"]
    ]
    primary = stage1 / "held-300K.data"
    write(data, primary)
    return primary


def coverage_fixture(tmp_path: Path) -> Path:
    """Make a 40%-wide gap that wraps across the periodic x boundary."""
    probability = np.ones((6, 20), dtype=float)
    probability[:, :4] = 0.0
    probability[:, 16:] = 0.0
    path = tmp_path / "coverage_probability.npz"
    np.savez_compressed(
        path,
        probability=probability,
        fractional_x=(np.arange(20) + 0.5) / 20,
        fractional_y=(np.arange(6) + 0.5) / 6,
    )
    return path


def test_identify_x_tunnel_centers_periodic_gap(tmp_path):
    primary = primary_fixture(tmp_path)
    coverage = coverage_fixture(tmp_path)
    plan = identify_x_tunnel(coverage, primary)

    assert plan["gap_wraps_periodic_x"] is True
    assert plan["empty_column_count"] == 8
    assert plan["gap_width_fraction"] == pytest.approx(0.4)
    assert plan["original_gap_center_fraction"] == pytest.approx(0.0)
    assert abs(plan["coordinate_shift_fraction_x"]) == pytest.approx(0.5)
    assert plan["shifted_wall_x_fraction"] == pytest.approx([0.3, 0.7])
    assert plan["shifted_packing_x_fraction"][0] > 0.3
    assert plan["shifted_packing_x_fraction"][1] < 0.7


def test_identify_x_tunnel_rejects_nonlocalized_maps(tmp_path):
    primary = primary_fixture(tmp_path)
    for name, value, message in (
        ("occupied", 1.0, "no x interval"),
        ("empty", 0.0, "entire x direction"),
    ):
        path = tmp_path / f"{name}.npz"
        np.savez_compressed(path, probability=np.full((4, 8), value))
        with pytest.raises(ValueError, match=message):
            identify_x_tunnel(path, primary)


def test_lego_build_seeds_gap_without_lateral_confinement(tmp_path):
    primary = primary_fixture(tmp_path)
    coverage = coverage_fixture(tmp_path)
    secondary = packed_fixture(
        tmp_path / "secondary-packed",
        shift=(50.0, 20.0, 100.0),
    )
    original = parse(primary)

    output = prepare_lego_stage2(
        ROOT / "tests/data/small-study.toml",
        primary,
        coverage,
        tmp_path / "lego",
        packed_xyz=secondary,
        packmol_seed=1234567,
        velocity_seed=7654321,
    )

    plan = json.loads((output / "lego_plan.json").read_text())
    profile = (output / "lego_tunnel_profile.csv").read_text()
    assert "mean_occupancy_probability" in profile
    assert profile.count(",1\n") == 8
    shifted = parse(output / "lego-stage1-shifted.data")
    topology = parse(output / "topology_output.lmp")
    lx = original.bounds["x"][1] - original.bounds["x"][0]
    dx = plan["coordinate_shift_fraction_x"] * lx
    primary_max_z = max(
        float(atom.fields[6]) for atom in original.sections["Atoms"]
    )
    expected_endpoint = primary_max_z + 30.0
    expected_final_endpoint = primary_max_z + 15.0
    assert plan["stage1_max_z_angstrom"] == pytest.approx(primary_max_z)
    assert plan["deposition_wall_endpoint_angstrom"] == pytest.approx(
        expected_endpoint
    )
    assert plan["final_deposition_wall_endpoint_angstrom"] == pytest.approx(
        expected_final_endpoint
    )
    assert plan["deposition_continuation_steps"] == 300000
    assert plan["lateral_confinement"] is False
    for before, after in zip(
        original.sections["Atoms"], shifted.sections["Atoms"]
    ):
        before_unwrapped = float(before.fields[4])
        after_unwrapped = float(after.fields[4]) + int(after.fields[7]) * lx
        assert after_unwrapped == pytest.approx(before_unwrapped + dx, abs=1e-7)
    for section in ("Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"):
        assert [
            row.fields for row in topology.sections[section]
        ][: shifted.count(section)] == [
            row.fields for row in shifted.sections[section]
        ]

    packmol = (output / "packmol.inp").read_text()
    assert "seed 1234567" in packmol
    assert "inside box 41.530000 2.000000" in packmol
    assert "83.570000 39.700000" in packmol

    deposition = (output / "deposition.in").read_text()
    assert "region lego_tunnel block" not in deposition
    assert "fix lego_wall" not in deposition
    assert "unfix lego_wall" not in deposition
    assert f"variable zend equal {expected_endpoint}" in deposition
    assert "compute deposition_zmax all reduce max z" in deposition
    assert (
        'variable safe_zstart equal "max(v_requested_zstart,'
        'v_measured_deposition_zmax+v_cutoff+1.0)"'
    ) in deposition
    assert "variable zstart equal ${safe_zstart}" in deposition
    assert "variable required_deposition_box_zhi equal" in deposition
    assert "change_box all z final" in deposition
    continuation = (output / "continue-deposition.in").read_text()
    assert "read_data deposited.data" in continuation
    assert f"variable continuation_zstart equal {expected_endpoint}" in continuation
    assert (
        f"variable continuation_zend equal {expected_final_endpoint}"
        in continuation
    )
    assert "run 300000 start 0 stop 300000" in continuation
    assert "velocity " not in continuation
    assert "write_data deposited-continued.data nocoeff" in continuation
    assert (
        f"variable hold_wall_hi equal {expected_final_endpoint}"
        in (output / "hold-300K.in").read_text()
    )
    for name in (
        "hold-300K.in",
        "hold-400K.in",
        "decompress-300K.in",
        "decompress-400K.in",
    ):
        assert "lego_" not in (output / name).read_text()

    manifest = json.loads((output / "assembly_manifest.json").read_text())
    assert (
        manifest["deposition_guidance"]["method"]
        == "periodic_x_gap_seeded"
    )
    assert "laterally unconstrained" in plan["bias_warning"]


def test_prepare_legacy_lego_continuation_preserves_deposited_data(tmp_path):
    primary = primary_fixture(tmp_path)
    coverage = coverage_fixture(tmp_path)
    secondary = packed_fixture(
        tmp_path / "secondary-packed",
        shift=(50.0, 20.0, 100.0),
    )
    output = prepare_lego_stage2(
        ROOT / "tests/data/small-study.toml",
        primary,
        coverage,
        tmp_path / "legacy-lego",
        packed_xyz=secondary,
    )
    deposited = output / "deposited.data"
    shutil.copyfile(output / "topology_output.lmp", deposited)
    deposited_hash = hashlib.sha256(deposited.read_bytes()).hexdigest()

    plan_path = output / "lego_plan.json"
    plan = json.loads(plan_path.read_text())
    for key in (
        "final_deposition_clearance_above_stage1_angstrom",
        "final_deposition_wall_endpoint_angstrom",
        "deposition_continuation_steps",
    ):
        plan.pop(key)
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    manifest_path = output / "assembly_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["deposition_guidance"] = plan
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    prepare_lego_continuation(
        ROOT / "tests/data/small-study.toml",
        output,
        additional_drop_angstrom=15.0,
        continuation_steps=300000,
    )

    assert hashlib.sha256(deposited.read_bytes()).hexdigest() == deposited_hash
    updated = json.loads(plan_path.read_text())
    expected_final = updated["deposition_wall_endpoint_angstrom"] - 15.0
    assert updated["final_deposition_wall_endpoint_angstrom"] == pytest.approx(
        expected_final
    )
    assert (
        f"variable continuation_zend equal {expected_final}"
        in (output / "continue-deposition.in").read_text()
    )
    assert (
        f"variable hold_wall_hi equal {expected_final}"
        in (output / "hold-300K.in").read_text()
    )


def test_lego_can_reproduce_confined_tunnel_control(tmp_path):
    primary = primary_fixture(tmp_path)
    coverage = coverage_fixture(tmp_path)
    secondary = packed_fixture(
        tmp_path / "secondary-packed",
        shift=(50.0, 20.0, 100.0),
    )

    output = prepare_lego_stage2(
        ROOT / "tests/data/small-study.toml",
        primary,
        coverage,
        tmp_path / "lego-confined",
        lateral_confinement=True,
        packed_xyz=secondary,
    )

    plan = json.loads((output / "lego_plan.json").read_text())
    deposition = (output / "deposition.in").read_text()
    assert plan["method"] == "periodic_x_gap_tunnel"
    assert plan["lateral_confinement"] is True
    assert "region lego_tunnel block" in deposition
    assert "fix lego_wall stage2 wall/region lego_tunnel harmonic" in deposition
    assert deposition.index("unfix lego_wall") < deposition.index(
        "write_data deposited.data"
    )
