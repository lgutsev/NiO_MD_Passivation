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
    "Roughness RMS (Å)",
    "Empty Exposed-Ni Sites (%)",
    "Primary Bound, Mean ± Block SEM (%)",
    "Secondary Bound, Mean ± Block SEM (%)",
    "Secondary Persistent Bound (%)",
    "Secondary P Height (Å)",
    "Secondary Tilt (deg)",
    "DCZ-4P Terminals 0 / 1 / 2 (%)",
    "Approx. Potential Step (V)",
    "Deposition Exchange Rate (events/site/ns)",
    "QC Status",
]

ALL_METRICS_HEADERS = [
    "System",
    "Replicate Cohort",
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
    "Roughness RMS (Å)",
    "Roughness Mean Absolute (Å)",
    "Roughness Peak-to-Valley (Å)",
    "Empty Exposed-Ni Sites (%)",
    "Z Dipole Moment (e·Å)",
    "Approx. Potential Step (V)",
    "Deposition Exchange Events",
    "Deposition Exchange Rate (events/site/ns)",
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
    "Comparability Note",
]

# Public name retained for callers that imported the original constant.
SUMMARY_HEADERS = DRAFT_HEADERS

COMPARISON_HEADERS = [
    "Replicate Cohort",
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
    "CoSAM Roughness RMS (Å)",
    "Sequential Roughness RMS (Å)",
    "Δ Roughness RMS: Seq-CoSAM (Å)",
    "CoSAM Approx. Potential Step (V)",
    "Sequential Approx. Potential Step (V)",
    "Δ Potential Step: Seq-CoSAM (V)",
    "CoSAM Deposition Exchange Rate (events/site/ns)",
    "Sequential Deposition Exchange Rate (events/site/ns)",
    "Δ Exchange Rate: Seq-CoSAM (events/site/ns)",
    "QC Status",
]

CORRELATION_HEADERS = [
    "Secondary Molecule",
    "Assembly",
    "Temperature (K)",
    "Film State",
    "N MD Seeds",
    "Projected Coverage, MD Mean (%)",
    "Roughness RMS, MD Mean (Å)",
    "Secondary Bound, MD Mean (%)",
    "Secondary Persistent Bound, MD Mean (%)",
    "Secondary Tilt, MD Mean (deg)",
    "Secondary Unbound Anchor Density, MD Mean (per nm²)",
    "Largest Void Patch, MD Mean (% of cell)",
    "Void Patch Count, MD Mean",
    "Approx. Potential Step, MD Mean (V)",
    "Deposition Exchange Rate, MD Mean (events/site/ns)",
    "Voc Mean (V)",
    "Voc SD (V)",
    "Jsc Mean (mA/cm²)",
    "Jsc SD (mA/cm²)",
    "FF Mean (%)",
    "FF SD (%)",
    "PCE Mean (%)",
    "PCE SD (%)",
    "N Experimental Batches",
    "N Independent Units",
    "N Measurements",
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


def _run_key(row: dict[str, Any]) -> str:
    return str(row.get("run_directory") or row["system"])


def _replicate_cohort(run_directory: str) -> str:
    parent = str(Path(run_directory).parent)
    return "primary" if parent == "." else parent


def _mean_available(rows: list[dict[str, Any]], field: str) -> float | None:
    values = [
        float(row[field])
        for row in rows
        if row.get(field) is not None
    ]
    return sum(values) / len(values) if values else None


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


def _comparability_note(
    coverage: dict[str, Any], interface: dict[str, Any]
) -> str | None:
    """Explain why a coverage/interface pair should not be joined, or None if OK.

    Both analyzers can run against the same nominal system/temperature/stage
    with different --last-frames/--stride choices (or even different
    trajectories). Coverage and interface metrics computed over mismatched
    windows are not directly comparable, so this gate catches that before the
    numbers land in a "Δ Seq-CoSAM" column.
    """
    coverage_frames = coverage.get("frames")
    interface_frames = interface.get("frames")
    if (
        coverage_frames is not None
        and interface_frames is not None
        and coverage_frames != interface_frames
    ):
        return (
            f"frame window mismatch: coverage {coverage_frames} frames "
            f"(step {coverage.get('first_step')}-{coverage.get('last_step')}) "
            f"vs interface {interface_frames} frames "
            f"(step {interface.get('first_step')}-{interface.get('last_step')})"
        )
    coverage_stride = coverage.get("stride")
    interface_stride = interface.get("stride")
    if (
        coverage_stride is not None
        and interface_stride is not None
        and coverage_stride != interface_stride
    ):
        return (
            f"stride mismatch: coverage stride {coverage_stride} "
            f"vs interface stride {interface_stride}"
        )
    coverage_first = coverage.get("first_step")
    coverage_last = coverage.get("last_step")
    interface_first = interface.get("first_step")
    interface_last = interface.get("last_step")
    if None not in (coverage_first, coverage_last, interface_first, interface_last) and (
        coverage_first != interface_first or coverage_last != interface_last
    ):
        return (
            f"step-range mismatch: coverage step {coverage_first}-{coverage_last} "
            f"vs interface step {interface_first}-{interface_last}"
        )
    return None


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
        (_run_key(row), row["temperature"], row["stage"]): row
        for row in coverage_results
    }
    deposition_interface_by_run = {
        _run_key(row): row
        for row in interface_results
        if row["stage"] == "deposition"
    }
    interface_by_key = {
        (_run_key(row), row["temperature"], row["stage"]): row
        for row in interface_results
        if row["stage"] != "deposition"
    }
    component_by_key = {
        (
            _run_key(row),
            row["temperature"],
            row["stage"],
            row["name"],
        ): row
        for row in interface_components
        if row["stage"] != "deposition"
    }
    terminal_by_key = {
        (
            _run_key(row),
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
        run_directory, temperature, stage = key
        source_row = interface or coverage
        if source_row is None:
            continue
        system = source_row["system"]
        deposition_interface = deposition_interface_by_run.get(run_directory)
        secondary_raw = (
            interface.get("secondary_name")
            if interface is not None
            else coverage.get("secondary_name") if coverage is not None else None
        )
        secondary = _display_component(secondary_raw)
        assembly = _assembly(system)
        primary_name = interface.get("primary_name") if interface else None
        primary_component = component_by_key.get(
            (run_directory, temperature, stage, primary_name)
        )
        secondary_component = component_by_key.get(
            (run_directory, temperature, stage, secondary_raw)
        )

        terminal_populations = {0: None, 1: None, 2: None}
        if secondary_raw == "dcz-4p":
            for count in terminal_populations:
                terminal = terminal_by_key.get(
                    (
                        run_directory,
                        temperature,
                        stage,
                        secondary_raw,
                        count,
                    )
                )
                if terminal is not None:
                    terminal_populations[count] = terminal["population"]

        comparability_note = (
            _comparability_note(coverage, interface)
            if coverage is not None and interface is not None
            else None
        )
        if coverage is None:
            status = "MISSING COVERAGE"
        elif interface is None:
            status = "MISSING INTERFACE"
        elif comparability_note is not None:
            status = "INCOMPARABLE"
        elif coverage["status"] != "OK" or interface["status"] != "OK":
            status = "CHECK"
        else:
            status = "OK"

        summary_rows.append(
            {
                "system": _system_label(system, assembly, secondary),
                "system_key": system,
                "replicate_cohort": _replicate_cohort(run_directory),
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
                "roughness_rms": (
                    coverage.get("roughness_rms") if coverage else None
                ),
                "roughness_mean_absolute": (
                    coverage.get("roughness_mean_absolute")
                    if coverage
                    else None
                ),
                "roughness_peak_to_valley": (
                    coverage.get("roughness_peak_to_valley")
                    if coverage
                    else None
                ),
                "empty_sites": interface.get("empty") if interface else None,
                "z_dipole_moment": (
                    interface.get("z_dipole_moment")
                    if interface
                    else None
                ),
                "z_dipole_potential_step_volts": (
                    interface.get("z_dipole_potential_step_volts")
                    if interface
                    else None
                ),
                "deposition_exchange_events": (
                    deposition_interface.get(
                        "cross_component_exchange_events"
                    )
                    if deposition_interface
                    else None
                ),
                "deposition_exchange_rate": (
                    deposition_interface.get(
                        "exchange_rate_per_site_per_ns"
                    )
                    if deposition_interface
                    else None
                ),
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
                "comparability_note": comparability_note,
                "run_directory": (
                    run_directory
                ),
            }
        )

    paired: dict[
        tuple[str, str, int | None, str], dict[str, dict[str, Any]]
    ] = {}
    for row in summary_rows:
        if row["assembly"] not in ("CoSAM", "Sequential"):
            continue
        paired.setdefault(
            (
                row["replicate_cohort"],
                row["secondary"],
                row["temperature"],
                row["film_state"],
            ),
            {},
        )[row["assembly"]] = row

    comparison_rows: list[dict[str, Any]] = []
    for (cohort, secondary, temperature, film_state), pair in sorted(
        paired.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            math.inf if item[0][2] is None else item[0][2],
            0 if item[0][3] == "Compressed" else 1,
        ),
    ):
        if "CoSAM" not in pair or "Sequential" not in pair:
            continue
        cosam = pair["CoSAM"]
        sequential = pair["Sequential"]
        # An INCOMPARABLE side means its own metrics were computed over a
        # mismatched coverage/interface window; excluded outright rather
        # than risk a Seq-CoSAM delta built on non-comparable numbers.
        if "INCOMPARABLE" in (cosam["status"], sequential["status"]):
            continue
        status = (
            "OK"
            if cosam["status"] == "OK" and sequential["status"] == "OK"
            else "CHECK"
        )
        comparison_rows.append(
            {
                "replicate_cohort": cohort,
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
                "cosam_roughness_rms": cosam["roughness_rms"],
                "sequential_roughness_rms": sequential[
                    "roughness_rms"
                ],
                "delta_roughness_rms": _delta(
                    cosam["roughness_rms"],
                    sequential["roughness_rms"],
                ),
                "cosam_z_dipole_potential_step_volts": cosam[
                    "z_dipole_potential_step_volts"
                ],
                "sequential_z_dipole_potential_step_volts": sequential[
                    "z_dipole_potential_step_volts"
                ],
                "delta_z_dipole_potential_step_volts": _delta(
                    cosam["z_dipole_potential_step_volts"],
                    sequential["z_dipole_potential_step_volts"],
                ),
                "cosam_deposition_exchange_rate": cosam[
                    "deposition_exchange_rate"
                ],
                "sequential_deposition_exchange_rate": sequential[
                    "deposition_exchange_rate"
                ],
                "delta_deposition_exchange_rate": _delta(
                    cosam["deposition_exchange_rate"],
                    sequential["deposition_exchange_rate"],
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
    temperature: int | None = None,
    film_state: str | None = None,
) -> list[dict[str, Any]]:
    """Join seed-averaged MD states against batch-aware experimental JV data.

    By default every available compressed/relaxed 300/400 K state is retained
    so the artificial compressed wall and the experimental 100 C anneal are
    not silently mapped onto one preferred simulation state. Optional filters
    are retained for callers that need one explicit state.
    """
    md_groups: dict[
        tuple[str, str, int | None, str], list[dict[str, Any]]
    ] = {}
    for row in summary_rows:
        if row["assembly"] not in ("CoSAM", "Sequential"):
            continue
        if row["status"] != "OK":
            continue
        if temperature is not None and row["temperature"] != temperature:
            continue
        if film_state is not None and row["film_state"] != film_state:
            continue
        md_groups.setdefault(
            (
                row["secondary"],
                row["assembly"],
                row["temperature"],
                row["film_state"],
            ),
            [],
        ).append(row)

    rows = []
    for md_key, md_rows in sorted(
        md_groups.items(),
        key=lambda item: (
            math.inf if item[0][2] is None else item[0][2],
            0 if item[0][3] == "Compressed" else 1,
            item[0][0],
            item[0][1],
        ),
    ):
        secondary, assembly, md_temperature, md_film_state = md_key
        experimental_key = (secondary, assembly)
        if experimental_key not in aggregated_experimental:
            continue
        device = aggregated_experimental[experimental_key]
        rows.append(
            {
                "secondary": secondary,
                "assembly": assembly,
                "temperature": md_temperature,
                "film_state": md_film_state,
                "md_seed_count": len(md_rows),
                "coverage": _mean_available(md_rows, "coverage"),
                "roughness_rms": _mean_available(
                    md_rows, "roughness_rms"
                ),
                "secondary_bound": _mean_available(
                    md_rows, "secondary_bound"
                ),
                "secondary_persistent": _mean_available(
                    md_rows, "secondary_persistent"
                ),
                "secondary_tilt": _mean_available(
                    md_rows, "secondary_tilt"
                ),
                "secondary_unbound_anchor_density": _mean_available(
                    md_rows, "secondary_unbound_anchor_density"
                ),
                "void_largest_patch": _mean_available(
                    md_rows, "void_largest_patch"
                ),
                "void_patch_count": _mean_available(
                    md_rows, "void_patch_count"
                ),
                "z_dipole_potential_step_volts": _mean_available(
                    md_rows, "z_dipole_potential_step_volts"
                ),
                "deposition_exchange_rate": _mean_available(
                    md_rows, "deposition_exchange_rate"
                ),
                "voc_mean": device["voc_v"]["mean"],
                "voc_std": device["voc_v"]["std"],
                "jsc_mean": device["jsc_ma_cm2"]["mean"],
                "jsc_std": device["jsc_ma_cm2"]["std"],
                "ff_mean": device["ff_percent"]["mean"],
                "ff_std": device["ff_percent"]["std"],
                "pce_mean": device["pce_percent"]["mean"],
                "pce_std": device["pce_percent"]["std"],
                "n": device["n"],
                "n_batches": device["n_batches"],
                "n_independent_units": device[
                    "n_independent_units"
                ],
                "n_measurements": device["n_measurements"],
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
    supplied device Voc/Jsc/FF/PCE data (see analysis.experimental).
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
            row["roughness_rms"],
            row["empty_sites"],
            _mean_sem(row["primary_bound"], row["primary_bound_sem"]),
            _mean_sem(row["secondary_bound"], row["secondary_bound_sem"]),
            row["secondary_persistent"],
            row["secondary_p_height"],
            row["secondary_tilt"],
            _terminal_populations(row),
            row["z_dipole_potential_step_volts"],
            row["deposition_exchange_rate"],
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
    for column_index in (5, 6, 7, 10, 11, 12, 14, 15):
        for cell in summary_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cell:
                value_cell.number_format = "0.00"
    for value_cell in summary_sheet.iter_cols(
        min_col=2, max_col=2, min_row=2
    ):
        for cell in value_cell:
            cell.number_format = "#,##0"
    summary_widths = [
        32,
        15,
        16,
        31,
        18,
        21,
        23,
        31,
        33,
        27,
        22,
        20,
        31,
        26,
        32,
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
            row["replicate_cohort"],
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
            row["roughness_rms"],
            row["roughness_mean_absolute"],
            row["roughness_peak_to_valley"],
            row["empty_sites"],
            row["z_dipole_moment"],
            row["z_dipole_potential_step_volts"],
            row["deposition_exchange_events"],
            row["deposition_exchange_rate"],
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
            row["comparability_note"],
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
    all_metrics_sheet.freeze_panes = "G2"
    for column_index in range(7, 36):
        for cells in all_metrics_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cells:
                value_cell.number_format = "0.00"
    for column_index in (5, 13, 22, 36, 37):
        for cells in all_metrics_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cells:
                value_cell.number_format = "#,##0"
    all_metrics_widths = [
        32, 24, 19, 20, 15, 16, 22, 21, 18, 24, 26, 20, 18, 26, 26, 20,
        27, 28, 23, 23, 25, 25, 33, 23, 24, 30, 25, 26, 30, 27, 22, 20,
        22, 22, 22, 17, 17, 18, 34, 60,
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
                row["temperature"],
                row["film_state"],
                row["md_seed_count"],
                row["coverage"],
                row["roughness_rms"],
                row["secondary_bound"],
                row["secondary_persistent"],
                row["secondary_tilt"],
                row["secondary_unbound_anchor_density"],
                row["void_largest_patch"],
                row["void_patch_count"],
                row["z_dipole_potential_step_volts"],
                row["deposition_exchange_rate"],
                row["voc_mean"],
                row["voc_std"],
                row["jsc_mean"],
                row["jsc_std"],
                row["ff_mean"],
                row["ff_std"],
                row["pce_mean"],
                row["pce_std"],
                row["n_batches"],
                row["n_independent_units"],
                row["n_measurements"],
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
        correlation_sheet.freeze_panes = "F2"
        for column_index in range(6, 24):
            for cells in correlation_sheet.iter_cols(
                min_col=column_index,
                max_col=column_index,
                min_row=2,
            ):
                for value_cell in cells:
                    value_cell.number_format = "0.00"
        for column_index in (3, 5, 24, 25, 26):
            for cells in correlation_sheet.iter_cols(
                min_col=column_index,
                max_col=column_index,
                min_row=2,
            ):
                for value_cell in cells:
                    value_cell.number_format = "#,##0"
        correlation_widths = [
            20, 15, 15, 16, 13, 24, 24, 23, 28, 24, 34, 30, 23,
            31, 36, 17, 15, 23, 21, 16, 14, 17, 15, 21, 20, 17,
        ]
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
            pce_column = CORRELATION_HEADERS.index("PCE Mean (%)") + 1
            pce_reference = Reference(
                correlation_sheet,
                min_col=pce_column,
                min_row=2,
                max_row=last_row,
            )
            for chart_title, x_header, x_title, anchor in (
                (
                    "Coverage vs. PCE",
                    "Projected Coverage, MD Mean (%)",
                    "Projected Coverage (%)",
                    "AB2",
                ),
                (
                    "Unbound Anchor Density vs. PCE",
                    "Secondary Unbound Anchor Density, MD Mean (per nm²)",
                    "Secondary Unbound Anchor Density (per nm²)",
                    "AB20",
                ),
                (
                    "Roughness vs. PCE",
                    "Roughness RMS, MD Mean (Å)",
                    "Roughness RMS (Å)",
                    "AB38",
                ),
            ):
                chart = ScatterChart()
                chart.title = chart_title
                chart.x_axis.title = x_title
                chart.y_axis.title = "PCE, Mean (%)"
                chart.style = 13
                x_col = CORRELATION_HEADERS.index(x_header) + 1
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
                chart.legend = None
                correlation_sheet.add_chart(chart, anchor)

    comparison_values = [
        [
            row["replicate_cohort"],
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
            row["cosam_roughness_rms"],
            row["sequential_roughness_rms"],
            row["delta_roughness_rms"],
            row["cosam_z_dipole_potential_step_volts"],
            row["sequential_z_dipole_potential_step_volts"],
            row["delta_z_dipole_potential_step_volts"],
            row["cosam_deposition_exchange_rate"],
            row["sequential_deposition_exchange_rate"],
            row["delta_deposition_exchange_rate"],
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
    comparison_sheet.freeze_panes = "E2"
    delta_columns = {7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37}
    for column_index in range(5, len(COMPARISON_HEADERS)):
        for cell in comparison_sheet.iter_cols(
            min_col=column_index,
            max_col=column_index,
            min_row=2,
        ):
            for value_cell in cell:
                value_cell.number_format = (
                    "+0.00;-0.00;0.00"
                    if column_index in delta_columns
                    else "0.00"
                )
    comparison_sheet.column_dimensions["A"].width = 24
    comparison_sheet.column_dimensions["B"].width = 21
    comparison_sheet.column_dimensions["C"].width = 15
    comparison_sheet.column_dimensions["D"].width = 16
    for index in range(5, len(COMPARISON_HEADERS)):
        comparison_sheet.column_dimensions[
            comparison_sheet.cell(1, index).column_letter
        ].width = 24
    comparison_sheet.column_dimensions[
        comparison_sheet.cell(1, len(COMPARISON_HEADERS)).column_letter
    ].width = 18

    ok_fill = PatternFill("solid", fgColor="C6EFCE")
    ok_font = Font(color="006100")
    check_fill = PatternFill("solid", fgColor="FFC7CE")
    check_font = Font(color="9C0006")
    for sheet, headers, row_count in (
        (summary_sheet, DRAFT_HEADERS, len(summary_rows)),
        (all_metrics_sheet, ALL_METRICS_HEADERS, len(summary_rows)),
        (comparison_sheet, COMPARISON_HEADERS, len(comparison_rows)),
    ):
        if row_count:
            column = sheet.cell(
                1, headers.index("QC Status") + 1
            ).column_letter
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
            "Each detailed row retains its run directory and assembly seeds. Block SEM measures time-series sampling within one trajectory and is not independent-seed uncertainty. Correlation rows average every available MD seed for the same assembly, molecule, temperature, and film state.",
        ],
        [
            "CoSAM vs Sequential deltas",
            "Every delta is Sequential minus CoSAM at matched replicate cohort, secondary molecule, temperature, and compressed/relaxed film state.",
        ],
        [
            "Comparability gate",
            "A row's coverage and interface analyses were not necessarily run with the same --last-frames/--stride/trajectory. Status is INCOMPARABLE (with a Comparability Note on the All Metrics sheet explaining why) whenever the two sides' analyzed frame counts, step ranges, or strides disagree. INCOMPARABLE rows are excluded outright from the CoSAM vs Sequential delta table rather than risk a delta built on non-comparable numbers; they remain visible on Draft Summary/All Metrics so the underlying per-side numbers are not hidden.",
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
        [
            "Surface roughness",
            "RMS, mean-absolute, and peak-to-valley fluctuations of the sampled molecular surface-height field on the periodic projected-coverage grid.",
        ],
        [
            "Z dipole and potential step",
            "Charge-weighted z dipole per frame when trajectory charge data are available. The reported potential step is the corresponding slab-capacitor approximation and should be interpreted as a comparative descriptor, not a converged work-function calculation.",
        ],
        [
            "Deposition site exchange",
            "Confirmed handoffs between distinct molecular components at one canonical exposed-Ni site, normalized by analyzed deposition time and the number of canonical exposed sites. Empty or shared intermediate frames do not count as separate handoffs.",
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
                f"Joins every available MD temperature and compressed/relaxed "
                f"state against externally supplied device Voc/Jsc/FF/PCE data "
                f"from {Path(experimental_csv).resolve()}. MD values average "
                "available assembly seeds within each state. Forward/reverse "
                "scans are first collapsed within a device (or conservatively "
                "within a batch when device_id is absent), independent units "
                "are averaged within batch, and batch means receive equal "
                "weight. Sample SD is across batch means when multiple batches "
                "exist, otherwise across independent units. Experimental "
                "batch, unit, and measurement counts are retained. The "
                "experimental 100 °C anneal is not silently equated with the "
                "artificially compressed 300 K hold; all states remain visible. "
                "Device data are supplied at run time and are not derived from "
                "this repository.",
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
    methods_sheet.sheet_properties.pageSetUpPr.fitToPage = True
    methods_sheet.page_setup.orientation = "landscape"
    methods_sheet.page_setup.fitToWidth = 1
    methods_sheet.page_setup.fitToHeight = 0
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
