from pathlib import Path
import csv
import json

import numpy as np
import pytest

from nio_md_prep.analysis.coverage import (
    _label_dependency,
    _periodic_patch_sizes,
    analyze_coverage,
    count_dump_frames,
    iter_dump_frames,
    rasterize_periodic,
)


def write_topology(path: Path) -> None:
    path.write_text(
        """LAMMPS data

3 atoms
3 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 12.011
2 15.999
3 58.6934

Atoms # full

1 1 1 0.0 0.1 5.0 10.0
2 2 2 0.0 5.0 5.0 10.0
3 0 3 2.0 5.0 5.0 0.0
""",
        encoding="utf-8",
    )


def dump_frame(step: int, first_x: float, second_x: float) -> str:
    return f"""ITEM: TIMESTEP
{step}
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp ff
0 10
0 10
0 20
ITEM: ATOMS id mol type q x y z
1 1 1 0.0 {first_x} 5.0 10.0
2 2 2 0.0 {second_x} 5.0 10.0
3 0 3 2.0 5.0 5.0 0.0
"""


def build_fixture(tmp_path: Path) -> Path:
    build = tmp_path / "build"
    build.mkdir()
    write_topology(build / "topology_output.lmp")
    (build / "assembly_manifest.json").write_text(
        json.dumps(
            {
                "components": [
                    {"component": "primary", "molecule_ids": [1, 1]},
                    {"component": "secondary", "molecule_ids": [2, 2]},
                ]
            }
        ),
        encoding="utf-8",
    )
    (build / "equilibration-300K.lammpstrj").write_text(
        dump_frame(0, 0.1, 5.0)
        + dump_frame(1000, 1.0, 5.0)
        + dump_frame(2000, 4.8, 5.0),
        encoding="utf-8",
    )
    return build


def test_periodic_disk_wraps_across_x_boundary():
    x = np.asarray([0.05])
    y = np.asarray([5.0])
    radii = np.asarray([1.0])
    mask = rasterize_periodic(
        x,
        y,
        radii,
        ((0.0, 10.0), (0.0, 10.0), (0.0, 20.0)),
        (100, 100),
    )
    assert mask[:, 0].any()
    assert mask[:, -1].any()
    assert mask.mean() * 100.0 == pytest.approx(3.14, abs=0.2)


def test_periodic_patch_sizes_merges_across_boundaries():
    label = _label_dependency()

    wrapping = np.zeros((6, 6), dtype=bool)
    wrapping[:, 0] = True
    wrapping[:, 5] = True
    assert _periodic_patch_sizes(np, label, wrapping) == [12]

    wrapping_y = np.zeros((6, 6), dtype=bool)
    wrapping_y[0, :] = True
    wrapping_y[5, :] = True
    assert _periodic_patch_sizes(np, label, wrapping_y) == [12]

    isolated = np.zeros((6, 6), dtype=bool)
    isolated[0, 0] = True
    isolated[3, 3] = True
    assert sorted(_periodic_patch_sizes(np, label, isolated)) == [1, 1]

    assert _periodic_patch_sizes(np, label, np.zeros((4, 4), dtype=bool)) == []


def test_dump_stream_and_count(tmp_path):
    path = tmp_path / "trajectory.lammpstrj"
    path.write_text(
        dump_frame(0, 0.1, 5.0) + dump_frame(1000, 0.2, 5.0),
        encoding="utf-8",
    )
    assert count_dump_frames(path) == 2
    frames = list(iter_dump_frames(path))
    assert [frame.step for frame in frames] == [0, 1000]
    assert frames[0].molecule_ids.tolist() == [1, 2, 0]


def test_analysis_uses_last_frames_and_reports_component_overlap(tmp_path):
    build = build_fixture(tmp_path)
    summary_path = analyze_coverage(
        build,
        last_frames=2,
        grid_spacing=0.1,
        blocks=2,
    )
    summary = json.loads(summary_path.read_text())
    assert summary["method"] == "periodic_projected_vdw_union"
    assert summary["frames_in_trajectory"] == 3
    assert summary["frames_analyzed"] == 2
    assert summary["analyzed_frame_indices"] == [1, 2]
    assert set(summary["components"]) == {"primary", "secondary"}
    assert summary["metrics"]["overlap"]["mean_percent"] > 0.0
    assert summary["metrics"]["void_patch_count"]["mean"] >= 0.0
    assert summary["metrics"]["void_largest_patch_percent"]["mean_percent"] >= 0.0
    assert summary["metrics"]["void_mean_patch_percent"]["mean_percent"] >= 0.0
    assert (summary_path.parent / "coverage_probability.npz").is_file()
    with (summary_path.parent / "coverage_timeseries.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert int(rows[0]["step"]) == 1000
    assert "void_patch_count" in rows[0]
    assert "void_largest_patch_percent" in rows[0]
    total = float(rows[-1]["coverage_total_percent"])
    primary = float(rows[-1]["coverage_primary_percent"])
    secondary = float(rows[-1]["coverage_secondary_percent"])
    overlap = float(rows[-1]["coverage_overlap_percent"])
    assert total == pytest.approx(primary + secondary - overlap)


def test_unmapped_positive_molecule_id_is_rejected(tmp_path):
    build = build_fixture(tmp_path)
    manifest = build / "assembly_manifest.json"
    manifest.write_text(
        json.dumps(
            {"components": [{"component": "primary", "molecule_ids": [1, 1]}]}
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="does not map all positive molecule IDs"):
        analyze_coverage(build, last_frames=1)
