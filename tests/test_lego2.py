from pathlib import Path
import json
import re

import numpy as np
import pytest

from nio_md_prep.build import build
from nio_md_prep.geometry import elements
from nio_md_prep.lammps import Record, parse, write
from nio_md_prep.lego2 import (
    _largest_true_rectangle,
    _periodic_components,
    identify_2d_void,
    prepare_lego2_stage2,
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
    """Localized rectangular void wrapping in both periodic directions."""
    probability = np.ones((10, 20), dtype=float)
    rows = [0, 1, 2, 7, 8, 9]
    columns = [0, 1, 2, 3, 16, 17, 18, 19]
    probability[np.ix_(rows, columns)] = 0.0
    path = tmp_path / "coverage_probability.npz"
    np.savez_compressed(
        path,
        probability=probability,
        fractional_x=(np.arange(20) + 0.5) / 20,
        fractional_y=(np.arange(10) + 0.5) / 10,
    )
    return path


def test_periodic_components_merge_across_x_and_y():
    mask = np.zeros((6, 8), dtype=bool)
    mask[np.ix_([0, 5], [0, 7])] = True
    components = _periodic_components(mask)
    assert len(components) == 1
    assert len(components[0]) == 4


def test_largest_true_rectangle_stays_inside_component():
    mask = np.array(
        [
            [True, True, False, False],
            [True, True, True, False],
            [False, True, True, False],
        ],
        dtype=bool,
    )
    top, bottom, left, right = _largest_true_rectangle(mask)
    rectangle = np.zeros_like(mask)
    rectangle[top:bottom, left:right] = True
    assert rectangle.sum() == 4
    assert bool(np.all(mask[rectangle]))


def test_identify_2d_void_centers_periodic_patch(tmp_path):
    primary = primary_fixture(tmp_path)
    coverage = coverage_fixture(tmp_path)
    plan = identify_2d_void(coverage, primary)

    assert plan["method"] == "periodic_2d_void_seeded"
    assert plan["largest_void_cell_count"] == 48
    assert plan["largest_void_area_fraction"] == pytest.approx(0.24)
    assert abs(plan["coordinate_shift_fraction_x"]) == pytest.approx(0.5)
    assert abs(plan["coordinate_shift_fraction_y"]) == pytest.approx(0.5)
    assert plan["shifted_void_rectangle_x_fraction"] == pytest.approx(
        [0.3, 0.7]
    )
    assert plan["shifted_void_rectangle_y_fraction"] == pytest.approx(
        [0.2, 0.8]
    )
    assert plan["packing_width_angstrom"] >= 12.0
    assert plan["packing_height_angstrom"] >= 12.0


def test_identify_2d_void_rejects_nonlocalized_maps(tmp_path):
    primary = primary_fixture(tmp_path)
    for name, value, message in (
        ("occupied", 1.0, "no 2D cells"),
        ("empty", 0.0, "entire surface"),
    ):
        path = tmp_path / f"{name}.npz"
        np.savez_compressed(path, probability=np.full((6, 8), value))
        with pytest.raises(ValueError, match=message):
            identify_2d_void(path, primary)


def test_lego2_build_seeds_largest_2d_void_without_confinement(tmp_path):
    primary = primary_fixture(tmp_path)
    coverage = coverage_fixture(tmp_path)
    secondary = packed_fixture(
        tmp_path / "secondary-packed",
        shift=(50.0, 20.0, 100.0),
    )
    original = parse(primary)

    output = prepare_lego2_stage2(
        ROOT / "tests/data/small-study.toml",
        primary,
        coverage,
        tmp_path / "lego2",
        packed_xyz=secondary,
        packmol_seed=1234567,
        velocity_seed=7654321,
    )

    plan = json.loads((output / "lego2_plan.json").read_text())
    shifted = parse(output / "lego2-stage1-shifted.data")
    topology = parse(output / "topology_output.lmp")
    assert (output / "lego2_void_map.csv").is_file()
    assert (output / "lego2_void_map.npz").is_file()
    assert plan["method"] == "periodic_2d_void_seeded"
    assert plan["lateral_confinement"] is False

    lx = original.bounds["x"][1] - original.bounds["x"][0]
    ly = original.bounds["y"][1] - original.bounds["y"][0]
    dx = plan["coordinate_shift_fraction_x"] * lx
    dy = plan["coordinate_shift_fraction_y"] * ly
    for before, after in zip(
        original.sections["Atoms"], shifted.sections["Atoms"]
    ):
        before_x = float(before.fields[4])
        before_y = float(before.fields[5])
        after_x = float(after.fields[4]) + int(after.fields[7]) * lx
        after_y = float(after.fields[5]) + int(after.fields[8]) * ly
        assert after_x == pytest.approx(before_x + dx, abs=1e-7)
        assert after_y == pytest.approx(before_y + dy, abs=1e-7)
    for section in ("Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"):
        assert [
            row.fields for row in topology.sections[section]
        ][: shifted.count(section)] == [
            row.fields for row in shifted.sections[section]
        ]

    packmol = (output / "packmol.inp").read_text()
    match = re.search(
        r"inside box ([0-9.]+) ([0-9.]+) ([0-9.]+) "
        r"([0-9.]+) ([0-9.]+) ([0-9.]+)",
        packmol,
    )
    assert match
    xlo, ylo, _zlo, xhi, yhi, _zhi = map(float, match.groups())
    assert xhi - xlo == pytest.approx(plan["packing_width_angstrom"], abs=1e-5)
    assert yhi - ylo == pytest.approx(plan["packing_height_angstrom"], abs=1e-5)

    deposition = (output / "deposition.in").read_text()
    assert "region lego_tunnel block" not in deposition
    assert "fix lego_wall" not in deposition
    assert "2D-void lego2" in deposition
    assert "write_data deposited.data nocoeff" in deposition
    continuation = (output / "continue-deposition.in").read_text()
    assert "read_data deposited.data" in continuation
    assert "write_data deposited-continued.data nocoeff" in continuation
    manifest = json.loads((output / "assembly_manifest.json").read_text())
    assert (
        manifest["deposition_guidance"]["method"]
        == "periodic_2d_void_seeded"
    )
