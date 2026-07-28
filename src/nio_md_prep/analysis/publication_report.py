"""Compact, publication-facing tables joined from the detailed analyses."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import math

from .experimental import aggregate_device_results, load_device_results
from .interfacial_report import collect_interfacial_results
from .report import (
    _openpyxl,
    _style_table_sheet,
    collect_coverage_results,
)


DRAFT_HEADERS = [
    "System",
    "Temperature (K)",
    "Film State",
    "Projected Coverage, Mean ± Block SEM (%)",
    "Uncovered Area (%)",
    "Empty Exposed-Ni Sites (%)",
    "Primary Bound, Mean ± Block SEM (%)",
    "Secondary Bound, Mean ± Block SEM (%)",
    "Secondary Persistent Bound (%)",
    "Secondary P Height (Å)",
    "Secondary Tilt (deg)",
    "DCZ-4P Terminals 0 / 1 / 2 (%)",
    "QC Status",
]

ALL_METRICS_HEADERS = [
    "System",
    "Assembly",
    "Secondary Molecule",
    "Temperature (K)",
    "Film State",
    "Total Projected Coverage (%)",
    "Coverage Block SEM (%)",
    "Uncovered Area (%)",
    "Primary Projected Coverage (%)",
    "Secondary Projected Coverage (%)",
    "Component Overlap (%)",
    "Void Patch Count",
    "Largest Void Patch (% of cell)",
    "Mean Void Patch (% of cell)",
    "Empty Exposed-Ni Sites (%)",
    "Primary Bound Molecules (%)",
    "Primary Bound Block SEM (%)",
    "Primary Unbound Anchor Density (per nm²)",
    "Secondary Bound Molecules (%)",
    "Secondary Bound Block SEM (%)",
    "Secondary Unbound Anchor Density (per nm²)",
    "Secondary Persistent Bound (%)",
    "Secondary P Height (Å)",
    "Secondary Tilt (deg)",
    "DCZ-4P 0-Terminal (%)",
    "DCZ-4P 1-Terminal (%)",
    "DCZ-4P 2-Terminal (%)",
    "Packmol Seed",
    "Velocity Seed",
    "QC Status",
    "Run Directory",
]

# Public name retained for callers that imported the original constant.
SUMMARY_HEADERS = DRAFT_HEADERS

COMPARISON_HEADERS = [
    "Secondary Molecule",
    "Temperature (K)",
    "Film State",
    "CoSAM Coverage (%)",
    "Sequential Coverage (%)",
    "Δ Coverage: Seq-CoSAM (%)",
    "CoSAM Empty Ni Sites (%)",
    "Sequential Empty Ni Sites (%)",
    "Δ Empty Sites: Seq-CoSAM (%)",
    "CoSAM Secondary Bound (%)",
    "Sequential Secondary Bound (%)",
    "Δ Secondary Bound: Seq-CoSAM (%)",
    "CoSAM Secondary Persistent (%)",
    "Sequential Secondary Persistent (%)",
    "Δ Secondary Persistent: Seq-CoSAM (%)",
    "CoSAM Secondary P Height (Å)",
    "Sequential Secondary P Height (Å)",
    "Δ Secondary P Height: Seq-CoSAM (Å)",
    "CoSAM Secondary Tilt (deg)",
    "Sequential Secondary Tilt (deg)",
    "Δ Secondary Tilt: Seq-CoSAM (deg)",
    "CoSAM Secondary Unbound Anchor Density (per nm²)",
    "Sequential Secondary Unbound Anchor Density (per nm²)",
    "Δ Unbound Anchor Density: Seq-CoSAM (per nm²)",
    "CoSAM Largest Void Patch (% of cell)",
    "Sequential Largest Void Patch (% of cell)",
    "Δ Largest Void Patch: Seq-CoSAM (% of cell)",
    "QC Status",
]

CORRELATION_HEADERS = [
    "Secondary Molecule",
    "Assembly",
    "Projected Coverage (%)",
    "Secondary Bound (%)",
    "Secondary Persistent Bound (%)",
    "Secondary Tilt (deg)",
    "Secondary Unbound Anchor Density (per nm²)",
    "Largest Void Patch (% of cell)",
    "Void Patch Count",
    "Voc, Mean ± SD (V)",
    "Jsc, Mean ± SD (mA/cm²)",
    "FF, Mean ± SD (%)",
    "PCE, Mean ± SD (%)",
    "N Devices",
]


def _display_component(name: str | None) -> str | None:
    labels = {
        "me-4pacz": "Me-4PACz",
        "stage1_primary": "Me-4PACz",
        "meo-2pacz": "MeO-2PACz",
        "meo-4padbc": "MeO-4PADBC",
        "dcz-4p": "DCZ-4P",
    }
    if name is None:
        return None
    return labels.get(name, name)


def _assembly(system: str) -> str:
    if "-then-" in system:
        return "Sequential"
    if system.endswith("-cosam"):
        return "CoSAM"
    if system.endswith("-high-dose"):
        return "High-dose control"
    return "Pure SAM"


def _system_label(
    system: str,
    assembly: str,
    secondary: str | None,
) -> str:
    if assembly == "Pure SAM":
        return "Me-4PACz"
    if assembly == "High-dose control":
        return "Me-4PACz (high dose)"
    separator = " → " if assembly == "Sequential" else " + "
    return f"Me-4PACz{separator}{secondary or system}"


def _film_state(stage: str) -> str:
    if stage.startswith("hold-"):
        return "Compressed"
    if stage.startswith("relax-"):
        return "Relaxed"
    return stage


def _delta(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return right - left


def _mean_sem(mean: float | None, sem: float | None) -> str | None:
    if mean is None:
        return None
    if sem is None:
        return f"{mean:.2f}"
    return f"{mean:.2f} ± {sem:.2f}"


def _terminal_populations(row: dict[str, Any]) -> str | None:
    populations = (row["dcz_0"], row["dcz_1"], row["dcz_2"])
    if all(value is None for value in populations):
        return None
    return " / ".join(
        "—" if value is None else f"{value:.2f}" for value in populations
    )


def build_publication_rows(
    coverage_results: list[dict[str, Any]],
    interface_results: list[dict[str, Any]],
    interface_components: list[dict[str, Any]],
    terminal_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Join detailed outputs into draft-summary and paired-comparison rows."""
    coverage_by_key = {
        (row["system"], row["temperature"], row["stage"]): row
        for row in coverage_results
    }
    interface_by_key = {
        (row["system"], row["temperature"], row["stage"]): row
        for row in interface_results
    }
    component_by_key = {
        (row["system"], row["temperature"], row["stage"], row["name"]): row
        for row in interface_components
    }
    terminal_by_key = {
        (
            row["system"],
            row["temperature"],
            row["stage"],
            row["component"],
            int(row["count"]),
        ): row
        for row in terminal_rows
    }

    keys = sorted(
        set(coverage_by_key) | set(interface_by_key),
        key=lambda key: (
            key[0],
            math.inf if key[1] is None else key[1],
            0 if key[2].startswith("hold-") else 1,
            key[2],
        ),
    )
    summary_rows: list[dict[str, Any]] = []
    for key in keys:
        coverage = coverage_by_key.get(key)
        interface = interface_by_key.get(key)
        system, temperature, stage = key
        secondary_raw = (
            interface.get("secondary_name")
            if interface is not None
            else coverage.get("secondary_name") if coverage is not None else None
        )
        secondary = _display_component(secondary_raw)
        assembly = _assembly(system)
        primary_name = interface.get("primary_name") if interface else None
        primary_component = component_by_key.get(
            (system, temperature, stage, primary_name)
        )
        secondary_component = component_by_key.get(
            (system, temperature, stage, secondary_raw)
        )

        terminal_populations = {0: None, 1: None, 2: None}
        if secondary_raw == "dcz-4p":
            for count in terminal_populations:
                terminal = terminal_by_key.get(
                    (system, temperature, stage, secondary_raw, count)
                )
                if terminal is not None:
                    terminal_populations[count] = terminal["population"]

        if coverage is None:
            status = "MISSING COVERAGE"
        elif interface is None:
            status = "MISSING INTERFACE"
        elif coverage["status"] != "OK" or interface["status"] != "OK":
            status = "CHECK"
        else:
            status = "OK"

        summary_rows.append(
            {
                "system": _system_label(system, assembly, secondary),
                "system_key": system,
                "assembly": assembly,
                "secondary": secondary,
                "temperature": temperature,
                "stage": stage,
                "film_state": _film_state(stage),
                "coverage": coverage.get("total") if coverage else None,
                "coverage_sem": coverage.get("sem") if coverage else None,
                "uncovered": coverage.get("uncovered") if coverage else None,
                "primary_coverage": (
                    coverage.get("primary_mean") if coverage else None
                ),
                "secondary_coverage": (
                    coverage.get("secondary_mean") if coverage else None
                ),
                "overlap": coverage.get("overlap") if coverage else None,
                "void_patch_count": (
                    coverage.get("void_patch_count") if coverage else None
                ),
                "void_largest_patch": (
                    coverage.get("void_largest_patch") if coverage else None
                ),
                "void_mean_patch": (
                    coverage.get("void_mean_patch") if coverage else None
                ),
                "empty_sites": interface.get("empty") if interface else None,
                "primary_bound": (
                    interface.get("primary_bound") if interface else None
                ),
                "primary_bound_sem": (
                    primary_component.get("bound_sem")
                    if primary_component
                    else None
                ),
                "primary_unbound_anchor_density": (
                    interface.get("primary_unbound_anchor_density")
                    if interface
                    else None
                ),
                "secondary_bound": (
                    interface.get("secondary_bound") if interface else None
                ),
                "secondary_bound_sem": (
                    secondary_component.get("bound_sem")
                    if secondary_component
                    else None
                ),
                "secondary_unbound_anchor_density": (
                    interface.get("secondary_unbound_anchor_density")
                    if interface
                    else None
                ),
                "secondary_persistent": (
                    interface.get("secondary_persistent")
                    if interface
                    else None
                ),
                "secondary_p_height": (
                    interface.get("secondary_p_height") if interface else None
                ),
                "secondary_tilt": (
                    interface.get("secondary_tilt") if interface else None
                ),
                "dcz_0": terminal_populations[0],
                "dcz_1": terminal_populations[1],
                "dcz_2": terminal_populations[2],
                "packmol_seed": (
                    interface.get("packmol_seed") if interface else None
                ),
                "velocity_seed": (
                    interface.get("velocity_seed") if interface else None
                ),
                "status": status,
                "run_directory": (
                    interface.get("run_directory") if interface else system
                ),
            }
        )

    paired: dict[
        tuple[str, int | None, str], dict[str, dict[str, Any]]
    ] = {}
    for row in summary_rows:
        if row["assembly"] not in ("CoSAM", "Sequential"):
            continue
        paired.setdefault(
            (row["secondary"], row["temperature"], row["film_state"]),
            {},
        )[row["assembly"]] = row

    comparison_rows: list[dict[str, Any]] = []
    for (secondary, temperature, film_state), pair in sorted(
        paired.items(),
        key=lambda item: (
            item[0][0],
            math.inf if item[0][1] is None else item[0][1],
            0 if item[0][2] == "Compressed" else 1,
        ),
    ):
        if "CoSAM" not in pair or "Sequential" not in pair:
            continue
        cosam = pair["CoSAM"]
        sequential = pair["Sequential"]
        status = (
            "OK"
            if cosam["status"] == "OK" and sequential["status"] == "OK"
            else "CHECK"
        )
        comparison_rows.append(
            {
                "secondary": secondary,
                "temperature": temperature,
                "film_state": film_state,
                "cosam_coverage": cosam["coverage"],
                "sequential_coverage": sequential["coverage"],
                "delta_coverage": _delta(
                    cosam["coverage"], sequential["coverage"]
                ),
                "cosam_empty": cosam["empty_sites"],
                "sequential_empty": sequential["empty_sites"],
                "delta_empty": _delta(
                    cosam["empty_sites"], sequential["empty_sites"]
                ),
                "cosam_secondary_bound": cosam["secondary_bound"],
                "sequential_secondary_bound": sequential["secondary_bound"],
                "delta_secondary_bound": _delta(
                    cosam["secondary_bound"],
                    sequential["secondary_bound"],
                ),
                "cosam_secondary_persistent": cosam[
                    "secondary_persistent"
                ],
                "sequential_secondary_persistent": sequential[
                    "secondary_persistent"
                ],
                "delta_secondary_persistent": _delta(
                    cosam["secondary_persistent"],
                    sequential["secondary_persistent"],
                ),
                "cosam_secondary_p_height": cosam["secondary_p_height"],
                "sequential_secondary_p_height": sequential[
                    "secondary_p_height"
                ],
                "delta_secondary_p_height": _delta(
                    cosam["secondary_p_height"],
                    sequential["secondary_p_height"],
                ),
                "cosam_secondary_tilt": cosam["secondary_tilt"],
                "sequential_secondary_tilt": sequential["secondary_tilt"],
                "delta_secondary_tilt": _delta(
                    cosam["secondary_tilt"], sequential["secondary_tilt"]
                ),
                "cosam_secondary_unbound_anchor_density": cosam[
                    "secondary_unbound_anchor_density"
                ],
                "sequential_secondary_unbound_anchor_density": sequential[
                    "secondary_unbound_anchor_density"
                ],
                "delta_secondary_unbound_anchor_density": _delta(
                    cosam["secondary_unbound_anchor_density"],
                    sequential["secondary_unbound_anchor_density"],
                ),
                "cosam_void_largest_patch": cosam["void_largest_patch"],
                "sequential_void_largest_patch": sequential[
                    "void_largest_patch"
                ],
                "delta_void_largest_patch": _delta(
                    cosam["void_largest_patch"],
                    sequential["void_largest_patch"],
                ),
                "status": status,
            }
        )
    return summary_rows, comparison_rows


def collect_publication_results(
    prepared_root: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Collect both detailed analyses and return publication-facing rows."""
    coverage_results, _ = collect_coverage_results(prepared_root)
    (
        interface_results,
        interface_components,
        terminal_rows,
        _,
        _,
        _,
    ) = collect_interfacial_results(prepared_root)
    return build_publication_rows(
        coverage_results,
        interface_results,
        interface_components,
        terminal_rows,
    )


def build_correlation_rows(
    summary_rows: list[dict[str, Any]],
    aggregated_experimental: dict[tuple[str, str], dict[str, Any]],
    *,
    temperature: int = 300,
    film_state: str = "Compressed",
) -> list[dict[str, Any]]:
    """Join MD structural metrics against aggregated experimental JV data.

    Matches at the 300 K compressed-hold condition by default -- the
    as-deposited film state closest to how these devices are actually
    fabricated and measured (room temperature, no post-deposition thermal
    relaxation step).
    """
    md_by_key = {
        (row["secondary"], row["assembly"]): row
        for row in summary_rows
        if row["assembly"] in ("CoSAM", "Sequential")
        and row["temperature"] == temperature
        and row["film_state"] == film_state
    }
    rows = []
    for key in sorted(set(md_by_key) & set(aggregated_experimental)):
        secondary, assembly = key
        md = md_by_key[key]
        device = aggregated_experimental[key]
        rows.append(
            {
                "secondary": secondary,
                "assembly": assembly,
                "coverage": md["coverage"],
                "secondary_bound": md["secondary_bound"],
                "secondary_persistent": md["secondary_persistent"],
                "secondary_tilt": md["secondary_tilt"],
                "secondary_unbound_anchor_density": md[
                    "secondary_unbound_anchor_density"
                ],
                "void_largest_patch": md["void_largest_patch"],
                "void_patch_count": md["void_patch_count"],
                "voc_mean": device["voc_v"]["mean"],
                "voc_std": device["voc_v"]["std"],
                "jsc_mean": device["jsc_ma_cm2"]["mean"],
                "jsc_std": device["jsc_ma_cm2"]["std"],
                "ff_mean": device["ff_percent"]["mean"],
                "ff_std": device["ff_percent"]["std"],
                "pce_mean": device["pce_percent"]["mean"],
                "pce_std": device["pce_percent"]["std"],
                "n": device["n"],
            }
        )
    return rows


def create_publication_workbook(
    prepared_root: Path,
    output: Path,
    *,
    experimental_csv: Path | None = None,
) -> Path:
    """Write a compact draft table and direct CoSAM/sequential comparison.

    When ``experimental_csv`` is given, an additional "Structure-Property
    Correlation" sheet joins the MD metrics above against externally
    supplied device Voc/Jsc/FF/PCE data (see analysis.experimental). Omitting
    it leaves the workbook identical to today's shape.
    """
    api = _openpyxl()
    Workbook = api["Workbook"]
    FormulaRule = api["FormulaRule"]
    PatternFill = api["PatternFill"]
    Font = api["Font"]
    Alignment = api["Alignment"]

    summary_rows, comparison_rows = collect_publication_results(prepared_root)
    workbook = Workbook()
    workbook.properties.creator = "nio-md-prep"
    workbook.properties.title = "NiO SAM publication summary"
    summary_sheet = workbook.active
    summary_sheet.title = "Draft Summary"
    comparison_sheet = workbook.create_sheet("CoSAM vs Sequential")
    all_metrics_sheet = workbook.create_sheet("All Metrics")
    correlation_rows: list[dict[str, Any]] = []
    correlation_sheet = None
    if experimental_csv is not None:
        aggregated_experimental = aggregate_device_results(
            load_device_results(experimental_csv)
        )
        correlation_rows = build_correlation_rows(
            summary_rows, aggregated_experimental
        )
        correlation_sheet = workbook.create_sheet(
            "Structure-Property Correlation"
        )
    methods_sheet = workbook.create_sheet("Methods")

    summary_values = [
        [
            row["system"],
            row["temperature"],
            row["film_state"],
            _mean_sem(row["coverage"], row["coverage_sem"]),
            row["uncovered"],
            row["empty_sites"],
            _mean_sem(row["primary_bound"], row["primary_bound_sem"]),
            _mean_sem(row["secondary_bound"], row["secondary_bound_sem"]),
            row["secondary_persistent"],
            row["secondary_p_height"],
            row["secondary_tilt"],
            _terminal_populations(row),
            row["status"],
        ]
        for row in summary_rows
    ]
    _style_table_sheet(
        summary_sheet,
        SUMMARY_HEADERS,
        summary_values,
        "DraftSummaryTable",
        api,
    )
    summary_sheet.freeze_panes = "B2"
    for column_index in range(5, 12):
        for cell in summary_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cell:
                value_cell.number_format = "0.00"
    summary_widths = [
        32,
        15,
        16,
        31,
        18,
        23,
        31,
        33,
        27,
        22,
        20,
        31,
        18,
    ]
    for index, width in enumerate(summary_widths, 1):
        summary_sheet.column_dimensions[
            summary_sheet.cell(1, index).column_letter
        ].width = width
    for row in summary_sheet.iter_rows(
        min_row=2,
        min_col=2,
        max_col=len(DRAFT_HEADERS),
    ):
        for cell in row:
            cell.alignment = Alignment(horizontal="center", vertical="top")

    all_metrics_values = [
        [
            row["system"],
            row["assembly"],
            row["secondary"],
            row["temperature"],
            row["film_state"],
            row["coverage"],
            row["coverage_sem"],
            row["uncovered"],
            row["primary_coverage"],
            row["secondary_coverage"],
            row["overlap"],
            row["void_patch_count"],
            row["void_largest_patch"],
            row["void_mean_patch"],
            row["empty_sites"],
            row["primary_bound"],
            row["primary_bound_sem"],
            row["primary_unbound_anchor_density"],
            row["secondary_bound"],
            row["secondary_bound_sem"],
            row["secondary_unbound_anchor_density"],
            row["secondary_persistent"],
            row["secondary_p_height"],
            row["secondary_tilt"],
            row["dcz_0"],
            row["dcz_1"],
            row["dcz_2"],
            row["packmol_seed"],
            row["velocity_seed"],
            row["status"],
            row["run_directory"],
        ]
        for row in summary_rows
    ]
    _style_table_sheet(
        all_metrics_sheet,
        ALL_METRICS_HEADERS,
        all_metrics_values,
        "AllPublicationMetricsTable",
        api,
    )
    all_metrics_sheet.freeze_panes = "F2"
    for column_index in range(6, 28):
        for cells in all_metrics_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cells:
                value_cell.number_format = "0.00"
    for column_index in (4, 28, 29):
        for cells in all_metrics_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cells:
                value_cell.number_format = "#,##0"
    all_metrics_widths = [
        32, 19, 20, 15, 16, 22, 21, 18, 24, 26, 20, 18, 26, 26, 23, 23,
        24, 30, 25, 26, 30, 27, 22, 20, 22, 22, 22, 17, 17, 18, 34,
    ]
    for index, width in enumerate(all_metrics_widths, 1):
        all_metrics_sheet.column_dimensions[
            all_metrics_sheet.cell(1, index).column_letter
        ].width = width

    if correlation_sheet is not None:
        correlation_values = [
            [
                row["secondary"],
                row["assembly"],
                row["coverage"],
                row["secondary_bound"],
                row["secondary_persistent"],
                row["secondary_tilt"],
                row["secondary_unbound_anchor_density"],
                row["void_largest_patch"],
                row["void_patch_count"],
                _mean_sem(row["voc_mean"], row["voc_std"]),
                _mean_sem(row["jsc_mean"], row["jsc_std"]),
                _mean_sem(row["ff_mean"], row["ff_std"]),
                _mean_sem(row["pce_mean"], row["pce_std"]),
                row["n"],
            ]
            for row in correlation_rows
        ]
        _style_table_sheet(
            correlation_sheet,
            CORRELATION_HEADERS,
            correlation_values,
            "StructurePropertyCorrelationTable",
            api,
        )
        correlation_sheet.freeze_panes = "C2"
        for column_index in range(3, 9):
            for cells in correlation_sheet.iter_cols(
                min_col=column_index,
                max_col=column_index,
                min_row=2,
            ):
                for value_cell in cells:
                    value_cell.number_format = "0.00"
        correlation_widths = [20, 15, 22, 20, 26, 20, 32, 26, 16, 18, 22, 15, 18, 12]
        for index, width in enumerate(correlation_widths, 1):
            correlation_sheet.column_dimensions[
                correlation_sheet.cell(1, index).column_letter
            ].width = width
        correlation_sheet.sheet_properties.pageSetUpPr.fitToPage = True
        correlation_sheet.page_setup.orientation = "landscape"
        correlation_sheet.page_setup.fitToWidth = 1
        correlation_sheet.page_setup.fitToHeight = 0

        if correlation_rows:
            ScatterChart = api["ScatterChart"]
            Series = api["Series"]
            Reference = api["Reference"]
            last_row = len(correlation_rows) + 1
            pce_reference = Reference(
                correlation_sheet, min_col=13, min_row=2, max_row=last_row
            )
            for chart_title, x_col, x_title, anchor in (
                ("Coverage vs. PCE", 3, "Projected Coverage (%)", "P2"),
                (
                    "Unbound Anchor Density vs. PCE",
                    7,
                    "Secondary Unbound Anchor Density (per nm²)",
                    "P20",
                ),
            ):
                chart = ScatterChart()
                chart.title = chart_title
                chart.x_axis.title = x_title
                chart.y_axis.title = "PCE, Mean (%)"
                chart.style = 13
                x_reference = Reference(
                    correlation_sheet,
                    min_col=x_col,
                    min_row=2,
                    max_row=last_row,
                )
                series = Series(
                    pce_reference, x_reference, title_from_data=False
                )
                series.marker.symbol = "circle"
                series.graphicalProperties.line.noFill = True
                chart.series.append(series)
                correlation_sheet.add_chart(chart, anchor)

    comparison_values = [
        [
            row["secondary"],
            row["temperature"],
            row["film_state"],
            row["cosam_coverage"],
            row["sequential_coverage"],
            row["delta_coverage"],
            row["cosam_empty"],
            row["sequential_empty"],
            row["delta_empty"],
            row["cosam_secondary_bound"],
            row["sequential_secondary_bound"],
            row["delta_secondary_bound"],
            row["cosam_secondary_persistent"],
            row["sequential_secondary_persistent"],
            row["delta_secondary_persistent"],
            row["cosam_secondary_p_height"],
            row["sequential_secondary_p_height"],
            row["delta_secondary_p_height"],
            row["cosam_secondary_tilt"],
            row["sequential_secondary_tilt"],
            row["delta_secondary_tilt"],
            row["cosam_secondary_unbound_anchor_density"],
            row["sequential_secondary_unbound_anchor_density"],
            row["delta_secondary_unbound_anchor_density"],
            row["cosam_void_largest_patch"],
            row["sequential_void_largest_patch"],
            row["delta_void_largest_patch"],
            row["status"],
        ]
        for row in comparison_rows
    ]
    _style_table_sheet(
        comparison_sheet,
        COMPARISON_HEADERS,
        comparison_values,
        "CoSAMSequentialComparisonTable",
        api,
    )
    comparison_sheet.freeze_panes = "D2"
    for column_index in range(4, 28):
        for cell in comparison_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cell:
                value_cell.number_format = (
                    "+0.00;-0.00;0.00"
                    if column_index in (6, 9, 12, 15, 18, 21, 24, 27)
                    else "0.00"
                )
    comparison_sheet.column_dimensions["A"].width = 21
    comparison_sheet.column_dimensions["B"].width = 15
    comparison_sheet.column_dimensions["C"].width = 16
    for index in range(4, 28):
        comparison_sheet.column_dimensions[
            comparison_sheet.cell(1, index).column_letter
        ].width = 24
    comparison_sheet.column_dimensions["AB"].width = 18

    ok_fill = PatternFill("solid", fgColor="C6EFCE")
    ok_font = Font(color="006100")
    check_fill = PatternFill("solid", fgColor="FFC7CE")
    check_font = Font(color="9C0006")
    for sheet, column, row_count in (
        (summary_sheet, "M", len(summary_rows)),
        (all_metrics_sheet, "AD", len(summary_rows)),
        (comparison_sheet, "AB", len(comparison_rows)),
    ):
        if row_count:
            sheet.conditional_formatting.add(
                f"{column}2:{column}{row_count + 1}",
                FormulaRule(
                    formula=[f'${column}2="OK"'],
                    fill=ok_fill,
                    font=ok_font,
                ),
            )
            sheet.conditional_formatting.add(
                f"{column}2:{column}{row_count + 1}",
                FormulaRule(
                    formula=[f'${column}2<>"OK"'],
                    fill=check_fill,
                    font=check_font,
                ),
            )

    for sheet in (summary_sheet, comparison_sheet, all_metrics_sheet):
        sheet.sheet_properties.pageSetUpPr.fitToPage = True
        sheet.page_setup.orientation = "landscape"
        sheet.page_setup.fitToWidth = 1
        sheet.page_setup.fitToHeight = 0

    methods = [
        ["Field", "Value"],
        [
            "Workbook purpose",
            "Compact, draft-facing summary joined from the detailed projected-coverage and anchor-resolved interfacial analyses.",
        ],
        [
            "Preliminary status",
            "These rows describe a single assembly seed. Block SEM measures time-series sampling within one trajectory and is not independent-seed uncertainty. Revise inferential claims after the additional assembly seeds finish.",
        ],
        [
            "CoSAM vs Sequential deltas",
            "Every delta is Sequential minus CoSAM at matched secondary molecule, temperature, and compressed/relaxed film state.",
        ],
        [
            "Projected coverage",
            "Periodic union of projected elemental van der Waals disks divided by instantaneous x/y cell area.",
        ],
        [
            "Bound molecule",
            "At least one phosphonate terminal has an O atom within the fixed interfacial contact cutoff of a canonical exposed Ni site.",
        ],
        [
            "Empty exposed-Ni sites",
            "Fraction of the fixed canonical exposed-Ni set not contacted by any phosphonate terminal.",
        ],
        [
            "DCZ-4P terminal states",
            "Population of DCZ-4P molecules with zero, one, or two anchored phosphonate termini; blank for other molecules.",
        ],
        [
            "Film state",
            "Compressed is the fixed low-wall hold. Relaxed is the final hold after gradual upper-wall retraction.",
        ],
        [
            "Void patches",
            "Contiguous uncovered-cell regions on the same periodic projected-coverage grid (periodic in x and y); largest/mean size reported as a percent of the full cell area.",
        ],
        [
            "Unbound anchor density",
            "Areal density of phosphonate anchor groups (P + its 3 O neighbors) with no O atom within the contact cutoff of any exposed Ni site, per nm^2 of periodic x/y cell.",
        ],
        ["Source root", str(Path(prepared_root).resolve())],
        ["Detailed coverage workbook", "coverage_summary.xlsx"],
        ["Detailed interface workbook", "interface_structure_summary.xlsx"],
        ["Summary rows", len(summary_rows)],
        ["Matched CoSAM/sequential rows", len(comparison_rows)],
    ]
    if correlation_sheet is not None:
        methods.append(
            [
                "Structure-Property Correlation",
                f"Joins MD metrics at the 300 K compressed hold against "
                f"externally supplied device Voc/Jsc/FF/PCE data from "
                f"{Path(experimental_csv).resolve()}. Experimental values are "
                "pooled mean +/- sample SD across every supplied replicate/batch "
                "for a given secondary molecule and assembly. Not derived from "
                "this repository; supplied at run time.",
            ]
        )
        methods.append(
            ["Correlated rows", len(correlation_rows)]
        )
    for row in methods:
        methods_sheet.append(row)
    methods_sheet.sheet_view.showGridLines = False
    for cell in methods_sheet[1]:
        cell.fill = PatternFill("solid", fgColor="17365D")
        cell.font = Font(color="FFFFFF", bold=True)
    methods_sheet.column_dimensions["A"].width = 31
    methods_sheet.column_dimensions["B"].width = 115
    methods_sheet.freeze_panes = "A2"
    for row_index, row in enumerate(methods_sheet.iter_rows(), 1):
        row[1].alignment = Alignment(wrap_text=True, vertical="top")
        if row_index > 1:
            text = str(row[1].value or "")
            methods_sheet.row_dimensions[row_index].height = max(
                18,
                15 * max(1, (len(text) + 104) // 105),
            )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output)
    return output
