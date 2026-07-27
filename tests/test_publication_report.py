from pathlib import Path

import pytest

from nio_md_prep.analysis import publication_report
from nio_md_prep.analysis.publication_report import (
    ALL_METRICS_HEADERS,
    COMPARISON_HEADERS,
    SUMMARY_HEADERS,
    build_publication_rows,
    create_publication_workbook,
)


def _fixture_rows():
    coverage = []
    interface = []
    components = []
    terminals = []
    systems = (
        (
            "me-4pacz-dcz-4p-cosam",
            71.0,
            29.0,
            42.0,
            35.0,
            14.0,
            62.0,
        ),
        (
            "me-4pacz-then-dcz-4p",
            76.0,
            24.0,
            27.0,
            18.0,
            9.0,
            74.0,
        ),
    )
    for index, (
        system,
        total,
        uncovered,
        secondary_bound,
        persistent,
        p_height,
        tilt,
    ) in enumerate(systems):
        coverage.append(
            {
                "system": system,
                "temperature": 300,
                "stage": "hold-300K",
                "total": total,
                "sem": 0.7,
                "uncovered": 100.0 - total,
                "primary_mean": total - 15.0,
                "secondary_name": "dcz-4p",
                "secondary_mean": 18.0,
                "overlap": 3.0,
                "status": "OK",
            }
        )
        interface.append(
            {
                "system": system,
                "temperature": 300,
                "stage": "hold-300K",
                "primary_name": "me-4pacz",
                "secondary_name": "dcz-4p",
                "empty": uncovered,
                "primary_bound": 82.0,
                "secondary_bound": secondary_bound,
                "secondary_persistent": persistent,
                "secondary_p_height": p_height,
                "secondary_tilt": tilt,
                "packmol_seed": 202607242 + index,
                "velocity_seed": 31415927 + index,
                "run_directory": system,
                "status": "OK",
            }
        )
        for component, sem in (("me-4pacz", 1.1), ("dcz-4p", 1.7)):
            components.append(
                {
                    "system": system,
                    "temperature": 300,
                    "stage": "hold-300K",
                    "name": component,
                    "bound_sem": sem,
                }
            )
        for count, population in ((0, 20.0), (1, 65.0), (2, 15.0)):
            terminals.append(
                {
                    "system": system,
                    "temperature": 300,
                    "stage": "hold-300K",
                    "component": "dcz-4p",
                    "count": count,
                    "population": population,
                }
            )
    return coverage, interface, components, terminals


def test_build_publication_rows_pairs_sequential_and_cosam():
    summary, comparison = build_publication_rows(*_fixture_rows())

    assert len(summary) == 2
    assert [row["assembly"] for row in summary] == ["CoSAM", "Sequential"]
    assert summary[0]["secondary"] == "DCZ-4P"
    assert summary[0]["dcz_1"] == pytest.approx(65.0)
    assert len(comparison) == 1
    assert comparison[0]["delta_coverage"] == pytest.approx(5.0)
    assert comparison[0]["delta_empty"] == pytest.approx(-5.0)
    assert comparison[0]["delta_secondary_bound"] == pytest.approx(-15.0)
    assert comparison[0]["delta_secondary_p_height"] == pytest.approx(-5.0)


def test_create_publication_workbook_has_draft_tables(tmp_path, monkeypatch):
    openpyxl = pytest.importorskip("openpyxl")
    coverage, interface, components, terminals = _fixture_rows()
    monkeypatch.setattr(
        publication_report,
        "collect_coverage_results",
        lambda _root: (coverage, []),
    )
    monkeypatch.setattr(
        publication_report,
        "collect_interfacial_results",
        lambda _root: (interface, components, terminals, [], [], []),
    )

    output = tmp_path / "publication_summary.xlsx"
    assert create_publication_workbook(Path("prepared"), output) == output

    workbook = openpyxl.load_workbook(output)
    assert workbook.sheetnames == [
        "Draft Summary",
        "CoSAM vs Sequential",
        "All Metrics",
        "Methods",
    ]
    assert [cell.value for cell in workbook["Draft Summary"][1]] == SUMMARY_HEADERS
    assert [
        cell.value for cell in workbook["All Metrics"][1]
    ] == ALL_METRICS_HEADERS
    assert [
        cell.value for cell in workbook["CoSAM vs Sequential"][1]
    ] == COMPARISON_HEADERS
    assert "DraftSummaryTable" in workbook["Draft Summary"].tables
    assert (
        "CoSAMSequentialComparisonTable"
        in workbook["CoSAM vs Sequential"].tables
    )
    assert "AllPublicationMetricsTable" in workbook["All Metrics"].tables
    assert workbook["Draft Summary"]["D2"].value == "71.00 ± 0.70"
    assert workbook["Draft Summary"]["L2"].value == "20.00 / 65.00 / 15.00"
    assert workbook["Draft Summary"]["M2"].value == "OK"
    assert workbook["CoSAM vs Sequential"]["F2"].value == pytest.approx(5.0)
    assert "single assembly seed" in workbook["Methods"]["B3"].value
