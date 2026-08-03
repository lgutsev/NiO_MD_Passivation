from pathlib import Path

import pytest

from nio_md_prep.analysis import publication_report
from nio_md_prep.analysis.publication_report import (
    ALL_METRICS_HEADERS,
    COMPARISON_HEADERS,
    CORRELATION_HEADERS,
    SUMMARY_HEADERS,
    build_publication_rows,
    build_correlation_rows,
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
        trajectory_sha256 = f"trajectory-{index}"
        topology_sha256 = f"topology-{index}"
        manifest_sha256 = f"manifest-{index}"
        shared_provenance = {
            "frames": 501,
            "first_step": 0,
            "last_step": 500000,
            "stride": 1,
            "schema_version": "2",
            "git_commit": "abc123",
            "trajectory_sha256": trajectory_sha256,
            "trajectory_size_bytes": 123456,
            "trajectory_mtime_utc": "2026-07-30T12:00:00+00:00",
            "topology_sha256": topology_sha256,
            "manifest_sha256": manifest_sha256,
            "requested_blocks": 5,
        }
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
                "void_patch_count": 3 if index == 0 else 1,
                "void_largest_patch": 12.0 if index == 0 else 4.0,
                "void_mean_patch": 8.0 if index == 0 else 4.0,
                "roughness_rms": 3.0 if index == 0 else 2.0,
                "roughness_mean_absolute": 2.0 if index == 0 else 1.0,
                "roughness_peak_to_valley": 10.0 if index == 0 else 7.0,
                "near_surface": 61.0 if index == 0 else 70.0,
                "anchor_conditioned": 45.0 if index == 0 else 56.0,
                "near_surface_void_largest_patch": 18.0 if index == 0 else 8.0,
                "p_near_surface_void_largest_patch": 25.0 if index == 0 else 12.0,
                "grid": 0.2,
                "radius_scale": 1.0,
                "hydrogen_included": True,
                "near_surface_height": 5.0,
                "surface_reference_depth": 4.0,
                "run_directory": system,
                "status": "OK",
                **shared_provenance,
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
                "primary_unbound_anchor_density": 0.0,
                "secondary_bound": secondary_bound,
                "secondary_unbound_anchor_density": 1.5 if index == 0 else 0.4,
                "secondary_persistent": persistent,
                "secondary_p_height": p_height,
                "secondary_tilt": tilt,
                "z_dipole_moment": -20.0 + index,
                "z_dipole_potential_step_volts": -0.4 + 0.1 * index,
                "packmol_seed": 202607242 + index,
                "velocity_seed": 31415927 + index,
                "cutoff": 3.25,
                "cutoff_method": "Ni-O distance",
                "surface_sites": 600,
                "persistence_threshold": 0.5,
                "exchange_min_dwell_frames": 2,
                "exchange_max_vacancy_gap_frames": 3,
                "run_directory": system,
                "status": "OK",
                **shared_provenance,
            }
        )
        interface.append(
            {
                "system": system,
                "stage": "deposition",
                "run_directory": system,
                "cross_component_exchange_events": 20 + index,
                "exchange_rate_per_site_per_ns": 1.2 + 0.3 * index,
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
                    "run_directory": system,
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
                    "run_directory": system,
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
    assert comparison[0]["delta_roughness_rms"] == pytest.approx(-1.0)
    assert comparison[0]["delta_deposition_exchange_rate"] == pytest.approx(0.3)
    assert comparison[0]["delta_near_surface_coverage"] == pytest.approx(9.0)
    assert comparison[0]["delta_p_near_surface_coverage"] == pytest.approx(11.0)


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
    assert workbook["Draft Summary"]["M2"].value == "20.00 / 65.00 / 15.00"
    assert workbook["Draft Summary"]["P2"].value == "OK"
    assert workbook["Draft Summary"]["Q2"].value == pytest.approx(61.0)
    assert workbook["CoSAM vs Sequential"]["G2"].value == pytest.approx(5.0)
    assert "available MD seed" in workbook["Methods"]["B3"].value


def test_create_publication_workbook_adds_correlation_sheet_when_experimental_csv_given(
    tmp_path, monkeypatch
):
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

    # Synthetic placeholder JV numbers -- not real device data.
    experimental_csv = tmp_path / "jv.csv"
    experimental_csv.write_text(
        "secondary,assembly,batch,scan_direction,voc_v,jsc_ma_cm2,ff_percent,pce_percent\n"
        "DCZ-4P,CoSAM,synthetic,R,1.100,10.000,70.0,7.000\n"
        "DCZ-4P,CoSAM,synthetic,F,1.100,10.000,72.0,7.500\n"
        "DCZ-4P,Sequential,synthetic,R,1.200,12.000,80.0,9.500\n"
        "DCZ-4P,Sequential,synthetic,F,1.200,12.000,82.0,10.000\n",
        encoding="utf-8",
    )

    output = tmp_path / "publication_summary.xlsx"
    assert (
        create_publication_workbook(
            Path("prepared"), output, experimental_csv=experimental_csv
        )
        == output
    )

    workbook = openpyxl.load_workbook(output)
    assert workbook.sheetnames == [
        "Draft Summary",
        "CoSAM vs Sequential",
        "All Metrics",
        "Structure-Property Correlation",
        "Methods",
    ]
    sheet = workbook["Structure-Property Correlation"]
    assert [cell.value for cell in sheet[1]] == CORRELATION_HEADERS
    assert "StructurePropertyCorrelationTable" in sheet.tables
    rows = {row[1].value: row for row in sheet.iter_rows(min_row=2)}
    assert set(rows) == {"CoSAM", "Sequential"}
    assert rows["CoSAM"][0].value == "DCZ-4P"
    assert rows["CoSAM"][21].value == pytest.approx(7.25)
    assert rows["CoSAM"][22].value == pytest.approx(0.0)
    assert rows["CoSAM"][23].value == 1
    assert rows["CoSAM"][24].value == 1
    assert rows["CoSAM"][25].value == 2
    assert rows["Sequential"][21].value == pytest.approx(9.75)
    assert rows["CoSAM"][26].value == pytest.approx(61.0)
    assert len(sheet._charts) == 3
    assert "Correlated rows" in {row[0].value for row in workbook["Methods"].iter_rows()}


def test_build_publication_rows_marks_mismatched_frame_windows_incomparable():
    coverage, interface, components, terminals = _fixture_rows()
    # Coverage says 501 frames over step 0-500000; interface (unchanged
    # fixture) implicitly has no provenance at all for the first system, so
    # give system 0 explicit, MISMATCHED provenance on both sides.
    coverage[0]["frames"] = 501
    coverage[0]["first_step"] = 0
    coverage[0]["last_step"] = 500000
    coverage[0]["stride"] = 1
    hold_interface = next(
        row for row in interface if row.get("stage") == "hold-300K" and row["system"] == coverage[0]["system"]
    )
    hold_interface["frames"] = 101
    hold_interface["first_step"] = 0
    hold_interface["last_step"] = 100000
    hold_interface["stride"] = 1

    summary, comparison = build_publication_rows(
        coverage, interface, components, terminals
    )

    mismatched = next(
        row for row in summary if row["system_key"] == coverage[0]["system"]
    )
    assert mismatched["status"] == "INCOMPARABLE"
    assert "coverage/interface provenance mismatch" in mismatched["comparability_note"]
    assert "frames=501/101" in mismatched["comparability_note"]
    # The mismatched pair's comparison row must be excluded outright, not
    # merely blanked -- only the still-comparable system remains.
    assert len(comparison) == 0

    other = next(
        row for row in summary if row["system_key"] != coverage[0]["system"]
    )
    assert other["status"] == "OK"
    assert other["comparability_note"] is None


def test_build_publication_rows_rejects_cross_assembly_window_mismatch():
    coverage, interface, components, terminals = _fixture_rows()
    sequential_system = "me-4pacz-then-dcz-4p"
    for rows in (coverage, interface):
        for row in rows:
            if row.get("system") == sequential_system and row.get("stage") == "hold-300K":
                row["frames"] = 101
                row["last_step"] = 100000

    summary, comparison = build_publication_rows(
        coverage, interface, components, terminals
    )

    assert comparison == []
    assert {row["status"] for row in summary} == {"INCOMPARABLE"}
    assert all(
        "cross-assembly coverage mismatch" in row["comparability_note"]
        for row in summary
    )


def test_build_correlation_rows_keeps_states_and_averages_md_seeds():
    summary, _ = build_publication_rows(*_fixture_rows())
    first = dict(summary[0])
    replicate = dict(first)
    replicate["run_directory"] = f"replicates/seed-02/{first['system_key']}"
    replicate["coverage"] = first["coverage"] + 2.0
    relaxed = dict(first)
    relaxed["stage"] = "relax-300K"
    relaxed["film_state"] = "Relaxed"
    relaxed["coverage"] = first["coverage"] - 4.0
    experimental = {
        ("DCZ-4P", "CoSAM"): {
            "n": 2,
            "n_batches": 1,
            "n_independent_units": 2,
            "n_measurements": 4,
            **{
                field: {"mean": mean, "std": 0.1}
                for field, mean in (
                    ("voc_v", 1.2),
                    ("jsc_ma_cm2", 18.0),
                    ("ff_percent", 80.0),
                    ("pce_percent", 17.0),
                )
            },
        }
    }

    rows = build_correlation_rows(
        [first, replicate, relaxed], experimental
    )

    assert {(row["film_state"], row["md_seed_count"]) for row in rows} == {
        ("Compressed", 2),
        ("Relaxed", 1),
    }
    compressed = next(row for row in rows if row["film_state"] == "Compressed")
    assert compressed["coverage"] == pytest.approx(first["coverage"] + 1.0)
