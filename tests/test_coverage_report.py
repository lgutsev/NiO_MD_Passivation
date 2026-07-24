from pathlib import Path
import json

import pytest

from nio_md_prep.analysis.report import (
    COMPONENT_HEADERS,
    RESULT_HEADERS,
    collect_coverage_results,
    create_coverage_workbook,
)


def write_summary(
    prepared: Path,
    system: str,
    temperature: int,
    *,
    primary: float,
    secondary: float,
    overlap: float,
    stage: str = "hold",
) -> None:
    output = prepared / system / f"coverage-analysis-{stage}-{temperature}K"
    output.mkdir(parents=True)
    total = primary + secondary - overlap
    summary = {
        "method": "periodic_projected_vdw_union",
        "trajectory": str(
            prepared / system / f"{stage}-{temperature}K.lammpstrj"
        ),
        "frames_analyzed": 100,
        "grid_target_spacing_angstrom": 0.2,
        "radius_scale": 1.0,
        "exclude_hydrogen": False,
        "metrics": {
            "total": {
                "mean_percent": total,
                "block_sem_percent": 0.8,
                "frame_std_percent": 2.1,
            },
            "uncovered": {
                "mean_percent": 100.0 - total,
                "block_sem_percent": 0.8,
                "frame_std_percent": 2.1,
            },
            "overlap": {
                "mean_percent": overlap,
                "block_sem_percent": 0.2,
                "frame_std_percent": 0.6,
            },
        },
        "components": {
            "primary-sam": {
                "molecule_id_range": [1, 180],
                "mean_percent": primary,
                "block_sem_percent": 0.7,
                "frame_std_percent": 1.9,
            },
            "secondary-sam": {
                "molecule_id_range": [181, 240],
                "mean_percent": secondary,
                "block_sem_percent": 0.3,
                "frame_std_percent": 0.8,
            },
        },
    }
    (output / "coverage_summary.json").write_text(
        json.dumps(summary),
        encoding="utf-8",
    )


def test_collect_results_sorts_holds_and_checks_union(tmp_path):
    prepared = tmp_path / "prepared"
    write_summary(
        prepared,
        "mixed-sam",
        400,
        primary=58.0,
        secondary=17.0,
        overlap=4.0,
    )
    write_summary(
        prepared,
        "mixed-sam",
        300,
        primary=55.0,
        secondary=15.0,
        overlap=3.0,
    )

    results, components = collect_coverage_results(prepared)

    assert [row["temperature"] for row in results] == [300, 400]
    assert results[0]["total"] == pytest.approx(67.0)
    assert results[0]["balance_qc"] == pytest.approx(0.0)
    assert results[0]["union_qc"] == pytest.approx(0.0)
    assert results[0]["status"] == "OK"
    assert [row["name"] for row in components[:2]] == [
        "primary-sam",
        "secondary-sam",
    ]


def test_create_workbook_contains_summary_components_and_chart(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    prepared = tmp_path / "prepared"
    write_summary(
        prepared,
        "mixed-sam",
        300,
        primary=55.0,
        secondary=15.0,
        overlap=3.0,
    )
    output = prepared / "coverage_summary.xlsx"

    assert create_coverage_workbook(prepared, output) == output

    workbook = openpyxl.load_workbook(output)
    assert workbook.sheetnames == ["Results", "Components", "Methods"]
    results = workbook["Results"]
    components = workbook["Components"]
    assert [cell.value for cell in results[1]] == RESULT_HEADERS
    assert [cell.value for cell in components[1]] == COMPONENT_HEADERS
    assert results["F2"].value == pytest.approx(67.0)
    assert results["T2"].value == "OK"
    assert len(results._charts) == 1
    assert "CoverageResultsTable" in results.tables
    assert components.max_row == 3
    assert "CoverageComponentsTable" in components.tables


def test_collect_results_includes_relaxed_holds(tmp_path):
    prepared = tmp_path / "prepared"
    write_summary(
        prepared,
        "mixed-sam",
        300,
        primary=52.0,
        secondary=14.0,
        overlap=2.0,
        stage="relax",
    )

    results, _ = collect_coverage_results(prepared)

    assert len(results) == 1
    assert results[0]["temperature"] == 300
    assert results[0]["stage"] == "relax-300K"


def test_workbook_requires_completed_summaries(tmp_path):
    with pytest.raises(FileNotFoundError, match="no compressed-hold or relaxed-hold"):
        collect_coverage_results(tmp_path / "prepared")
