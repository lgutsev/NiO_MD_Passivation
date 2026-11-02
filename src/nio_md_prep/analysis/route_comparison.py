"""Replica-level comparison of molecular-assembly routes.

The unit of inference is one independently seeded run.  Frames from one
trajectory are never pooled as independent samples: each run contributes one
window mean, and uncertainty comes from the spread across independent
replicate cohorts.

Design facts encoded here (from the study configurations and run seeds):

* Within a molecular pair, co-assembly (CoSAM) and conventional sequential
  runs use the same inventory (180 Me-4PACz + N secondary), surface, cell,
  temperatures and hold/relaxation schedule, so CoSAM vs sequential is a
  matched-inventory route comparison.
* Each replicate cohort (``prepared``, ``prepared-rerun``,
  ``prepared-rerun2``) has its own Packmol/velocity seeds.  The three
  sequential runs of one cohort share that cohort's Me-4PACz stage-1 film, so
  route effects are estimated as within-cohort paired differences with the
  cohort as the replicate (n = 3).
* Gap-seeded I/II runs exist once each and were built from cohort 1's
  stage-1 film with an additional vertical deposition step; they are reported
  descriptively against the conventional-sequential replicate range, never
  as an independent replicate.
* 300 K/400 K branches and hold/relax stages of one run share their deposited
  structure, so they are repeated observations of one run, not replicates.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable
import json
import math
import re

COHORTS = {"prepared": 1, "prepared-rerun": 2, "prepared-rerun2": 3}
GAP_SEEDED = {"prepared-lego": "gap-seeded-I", "prepared-lego2": "gap-seeded-II"}
MODIFIERS = ("dcz-4p", "meo-2pacz", "meo-4padbc")
SECONDARY_COUNTS = {"dcz-4p": 59, "meo-2pacz": 107, "meo-4padbc": 77}
PRIMARY_COUNT = 180
DOSE_COUNT = 261

COVERAGE_METRICS = {
    "Total Coverage (%)": "projected_coverage_percent",
    "Near-Surface Coverage (%)": "near_surface_coverage_percent",
    "Anchor-Conditioned Coverage (%)": "anchor_conditioned_coverage_percent",
    "Largest Void Patch (% of cell)": "largest_void_percent",
    "Void Patch Count": "void_patch_count",
    "Mean Void Patch (% of cell)": "mean_void_percent",
    "Roughness RMS (Å)": "roughness_rms_angstrom",
}
INTERFACE_METRICS = {
    "Empty Sites (%)": "empty_ni_sites_percent",
    "Primary Bound (%)": "primary_bound_percent",
    "Secondary Bound (%)": "secondary_bound_percent",
    "Primary Tilt (deg)": "primary_tilt_deg",
    "Secondary Tilt (deg)": "secondary_tilt_deg",
    "Primary P Height (Å)": "primary_p_height_angstrom",
    "Secondary P Height (Å)": "secondary_p_height_angstrom",
    "Primary Core Height (Å)": "primary_core_height_angstrom",
    "Secondary Core Height (Å)": "secondary_core_height_angstrom",
}
SEM_COLUMN = "Block SEM (%)"
# "|"-separated paths into accessibility_summary.json
ACCESSIBILITY_METRICS = {
    "accessible_canopy_probe0_percent": "accessibility|canopy_probe0A|accessible_percent|mean",
    "accessible_canopy_probe1.4_percent": "accessibility|canopy_probe1.4A|accessible_percent|mean",
    "accessible_canopy_probe3_percent": "accessibility|canopy_probe3A|accessible_percent|mean",
    "accessible_h5_probe1.4_percent": "accessibility|h5A_probe1.4A|accessible_percent|mean",
    "accessible_h8_probe1.4_percent": "accessibility|h8A_probe1.4A|accessible_percent|mean",
    "accessible_h12_probe1.4_percent": "accessibility|h12A_probe1.4A|accessible_percent|mean",
    "gap_density_canopy_probe1.4_per_nm2": "accessibility|canopy_probe1.4A|gap_count_per_nm2|mean",
    "largest_gap_canopy_probe1.4_percent": "accessibility|canopy_probe1.4A|largest_gap_percent|mean",
    "persistent_exposed_canopy_probe1.4": "persistence|canopy_probe1.4A|persistently_exposed_fraction",
    "persistent_exposed_h8_probe1.4": "persistence|h8A_probe1.4A|persistently_exposed_fraction",
    "first_last_jaccard_canopy_probe1.4": "persistence|canopy_probe1.4A|first_last_jaccard",
    "thickness_mean_angstrom": "time_series|thickness_mean_angstrom|mean",
    "thickness_ge_10A_percent": "time_series|thickness_ge_10A_percent|mean",
    "bound_molecules_per_nm2": "time_series|bound_molecules_per_nm2|mean",
}


def classify_system(root: str, system: str) -> dict[str, Any]:
    """Map a prepared root and system slug onto route, modifier, cohort."""
    if root in GAP_SEEDED:
        route, cohort = GAP_SEEDED[root], 1
    elif root in COHORTS:
        cohort = COHORTS[root]
        if system == "me-4pacz-alone":
            route = "reference"
        elif system == "me-4pacz-high-dose":
            route = "dose-control"
        elif system.endswith("-cosam"):
            route = "cosam"
        elif "-then-" in system:
            route = "sequential"
        else:
            raise ValueError(f"unrecognized system {system!r} in {root!r}")
    else:
        raise ValueError(f"unrecognized prepared root {root!r}")
    modifier = next((m for m in MODIFIERS if m in system), "none")
    if route in {"reference"}:
        inventory = PRIMARY_COUNT
    elif route == "dose-control":
        inventory = DOSE_COUNT
    else:
        inventory = PRIMARY_COUNT + SECONDARY_COUNTS[modifier]
    return {
        "route": route,
        "modifier": modifier,
        "cohort": cohort,
        "independent_replicate": root in COHORTS,
        "inventory_molecules": inventory,
    }


def _stage_fields(stage: str) -> tuple[str, int]:
    match = re.fullmatch(r"(hold|relax)-(\d+)K", str(stage))
    if not match:
        raise ValueError(f"unrecognized stage {stage!r}")
    state = "compressed" if match.group(1) == "hold" else "relaxed"
    return state, int(match.group(2))


def load_workbook_records(lammps_directory: Path) -> list[dict[str, Any]]:
    """Long-format records from the consolidated coverage/interface workbooks."""
    import pandas as pd

    records = []
    lammps_directory = Path(lammps_directory)
    for root in [*COHORTS, *GAP_SEEDED]:
        coverage_path = lammps_directory / root / "coverage_summary.xlsx"
        interface_path = lammps_directory / root / "interface_structure_summary.xlsx"
        if not coverage_path.is_file() or not interface_path.is_file():
            continue
        for path, metrics in ((coverage_path, COVERAGE_METRICS), (interface_path, INTERFACE_METRICS)):
            frame = pd.read_excel(path, sheet_name="Results")
            for _, row in frame.iterrows():
                if not re.fullmatch(r"(hold|relax)-\d+K", str(row["Stage"])):
                    continue  # deposition rows describe assembly, not the held film
                info = classify_system(root, str(row["System"]))
                state, temperature = _stage_fields(row["Stage"])
                base = {
                    "prepared_root": root,
                    "system": str(row["System"]),
                    "stage": str(row["Stage"]),
                    "film_state": state,
                    "temperature_K": temperature,
                    "frames": int(row["Frames Analyzed"]),
                    "trajectory": str(row.get("Trajectory", "")),
                    "packmol_seed": row.get("Packmol Seed"),
                    "velocity_seed": row.get("Velocity Seed"),
                    "source": str(path.name),
                    **info,
                }
                for column, metric in metrics.items():
                    value = row.get(column)
                    if value is None or (isinstance(value, float) and math.isnan(value)):
                        continue
                    record = {**base, "metric": metric, "value": float(value)}
                    if metric == "projected_coverage_percent" and SEM_COLUMN in row:
                        record["within_run_block_sem"] = float(row[SEM_COLUMN])
                    records.append(record)
        records.extend(_derived_bound_records(records, root))
    return records


def _derived_bound_records(records, root) -> list[dict[str, Any]]:
    """Bound-molecule counts from bound percentages and the known inventory."""
    derived = []
    keyed: dict[tuple, dict[str, float]] = {}
    templates: dict[tuple, dict] = {}
    for record in records:
        if record["prepared_root"] != root:
            continue
        key = (record["system"], record["stage"])
        keyed.setdefault(key, {})[record["metric"]] = record["value"]
        templates.setdefault(key, record)
    for key, values in keyed.items():
        if "primary_bound_percent" not in values:
            continue
        template = templates[key]
        primary = DOSE_COUNT if template["route"] == "dose-control" else PRIMARY_COUNT
        count = values["primary_bound_percent"] / 100.0 * primary
        secondary = 0.0
        if template["modifier"] != "none" and "secondary_bound_percent" in values:
            secondary = values["secondary_bound_percent"] / 100.0 * SECONDARY_COUNTS[template["modifier"]]
        base = {k: v for k, v in template.items() if k not in {"metric", "value", "within_run_block_sem", "source"}}
        derived.append({**base, "metric": "bound_molecules", "value": count + secondary, "source": "derived"})
        derived.append({**base, "metric": "secondary_bound_molecules", "value": secondary, "source": "derived"})
        derived.append({
            **base,
            "metric": "bound_inventory_percent",
            "value": 100.0 * (count + secondary) / template["inventory_molecules"],
            "source": "derived",
        })
    return derived


def load_accessibility_records(manifest_csv: Path, metric_paths: dict[str, str]) -> list[dict[str, Any]]:
    """Long-format records from per-run ``accessibility_summary.json`` files.

    ``manifest_csv`` needs columns prepared_root, system, stage, summary_json.
    ``metric_paths`` maps output metric names to "|"-separated JSON paths, e.g.
    ``{"accessible_canopy_1.4": "accessibility|canopy_probe1.4A|accessible_percent|mean"}``.
    """
    import csv

    records = []
    with Path(manifest_csv).open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            summary = json.loads(Path(row["summary_json"]).read_text(encoding="utf-8"))
            info = classify_system(row["prepared_root"], row["system"])
            state, temperature = _stage_fields(row["stage"])
            for metric, dotted in metric_paths.items():
                value: Any = summary
                for part in _split_path(dotted):
                    value = value.get(part) if isinstance(value, dict) else None
                if value is None:
                    continue
                records.append({
                    "prepared_root": row["prepared_root"], "system": row["system"],
                    "stage": row["stage"], "film_state": state, "temperature_K": temperature,
                    "frames": summary.get("frames_analyzed"), "source": row["summary_json"],
                    "metric": metric, "value": float(value), **info,
                })
    return records


def _split_path(path: str) -> list[str]:
    # JSON keys such as "canopy_probe1.4A" contain dots, so paths use "|".
    return path.split("|")


def t_interval(values: Iterable[float], confidence: float = 0.95) -> dict[str, float]:
    from scipy import stats
    import numpy as np

    array = np.asarray([v for v in values if math.isfinite(v)], dtype=float)
    n = len(array)
    if n == 0:
        return {"n": 0, "mean": math.nan, "sd": math.nan, "ci_low": math.nan, "ci_high": math.nan}
    mean = float(array.mean())
    if n == 1:
        return {"n": 1, "mean": mean, "sd": math.nan, "ci_low": math.nan, "ci_high": math.nan,
                "min": mean, "max": mean}
    sd = float(array.std(ddof=1))
    half = float(stats.t.ppf(0.5 + confidence / 2.0, n - 1) * sd / math.sqrt(n))
    return {"n": n, "mean": mean, "sd": sd, "sem": sd / math.sqrt(n),
            "ci_low": mean - half, "ci_high": mean + half,
            "min": float(array.min()), "max": float(array.max())}


def replica_statistics(records) -> list[dict[str, Any]]:
    """Mean/SD/95% CI across independent replicate runs per condition."""
    groups: dict[tuple, list[dict]] = {}
    for record in records:
        key = (record["metric"], record["route"], record["modifier"], record["temperature_K"], record["film_state"])
        groups.setdefault(key, []).append(record)
    rows = []
    for (metric, route, modifier, temperature, state), members in sorted(groups.items(), key=str):
        values = [m["value"] for m in members if m["independent_replicate"]]
        single = [m["value"] for m in members if not m["independent_replicate"]]
        stats = t_interval(values if values else single)
        sems = [m["within_run_block_sem"] for m in members if "within_run_block_sem" in m]
        rows.append({
            "metric": metric, "route": route, "modifier": modifier,
            "temperature_K": temperature, "film_state": state,
            "independent_runs": len(values), "non_independent_runs": len(single),
            **stats,
            "mean_within_run_block_sem": (sum(sems) / len(sems)) if sems else math.nan,
        })
    return rows


def paired_route_differences(records, treatment: str = "sequential", control: str = "cosam"):
    """Within-cohort treatment-minus-control differences, cohort as replicate.

    Adds a modifier-averaged row ("all") whose replicate is the cohort mean of
    the three modifier differences, which respects the shared stage-1 film.
    """
    table: dict[tuple, dict[str, float]] = {}
    for record in records:
        if record["route"] not in {treatment, control} or not record["independent_replicate"]:
            continue
        key = (record["metric"], record["modifier"], record["temperature_K"], record["film_state"], record["cohort"])
        table.setdefault(key, {})[record["route"]] = record["value"]
    differences: dict[tuple, dict[int, float]] = {}
    for (metric, modifier, temperature, state, cohort), pair in table.items():
        if treatment in pair and control in pair:
            differences.setdefault((metric, modifier, temperature, state), {})[cohort] = pair[treatment] - pair[control]
    pooled: dict[tuple, dict[int, list[float]]] = {}
    for (metric, modifier, temperature, state), by_cohort in differences.items():
        for cohort, value in by_cohort.items():
            pooled.setdefault((metric, "all", temperature, state), {}).setdefault(cohort, []).append(value)
    for key, by_cohort in pooled.items():
        if all(len(v) == len(MODIFIERS) for v in by_cohort.values()):
            differences[key] = {c: sum(v) / len(v) for c, v in by_cohort.items()}
    rows = []
    for (metric, modifier, temperature, state), by_cohort in sorted(differences.items(), key=str):
        values = [by_cohort[c] for c in sorted(by_cohort)]
        stats = t_interval(values)
        positive = sum(v > 0 for v in values)
        negative = sum(v < 0 for v in values)
        if stats["n"] >= 2 and (stats["ci_low"] > 0 or stats["ci_high"] < 0):
            verdict = "supported"
        elif stats["n"] >= 2 and (positive == stats["n"] or negative == stats["n"]):
            verdict = "consistent sign, CI includes 0"
        else:
            verdict = "unresolved"
        rows.append({
            "metric": metric, "modifier": modifier, "temperature_K": temperature,
            "film_state": state, "contrast": f"{treatment} - {control}",
            **{f"cohort_{c}": by_cohort[c] for c in sorted(by_cohort)},
            **stats, "positive": positive, "negative": negative, "verdict": verdict,
        })
    return rows


def paired_vs_reference(records, reference: str = "reference") -> list[dict[str, Any]]:
    """Route-minus-Me-4PACz-reference differences within each cohort.

    In every cohort the sequential runs start from that cohort's reference
    film, so this contrast shows whether a CoSAM-vs-sequential gap comes
    from sequential gaining or from CoSAM losing relative to the primary film.
    """
    base: dict[tuple, float] = {}
    for record in records:
        if record["route"] == reference and record["independent_replicate"]:
            base[(record["metric"], record["temperature_K"], record["film_state"], record["cohort"])] = record["value"]
    shifted = []
    for record in records:
        if record["route"] in {"cosam", "sequential"} and record["independent_replicate"]:
            key = (record["metric"], record["temperature_K"], record["film_state"], record["cohort"])
            if key in base:
                shifted.append({**record, "value": record["value"] - base[key]})
    rows = []
    for stat in replica_statistics(shifted):
        stat = dict(stat)
        stat["contrast"] = f"{stat['route']} - {reference}"
        rows.append(stat)
    return rows


def gap_seeded_context(records) -> list[dict[str, Any]]:
    """Place each single gap-seeded run against the sequential replicate range."""
    sequential: dict[tuple, list[float]] = {}
    cosam: dict[tuple, list[float]] = {}
    cohort1_seq: dict[tuple, float] = {}
    for record in records:
        key = (record["metric"], record["modifier"], record["temperature_K"], record["film_state"])
        if record["route"] == "sequential" and record["independent_replicate"]:
            sequential.setdefault(key, []).append(record["value"])
            if record["cohort"] == 1:
                cohort1_seq[key] = record["value"]
        elif record["route"] == "cosam" and record["independent_replicate"]:
            cosam.setdefault(key, []).append(record["value"])
    rows = []
    for record in records:
        if not record["route"].startswith("gap-seeded"):
            continue
        key = (record["metric"], record["modifier"], record["temperature_K"], record["film_state"])
        seq = sequential.get(key, [])
        co = cosam.get(key, [])
        rows.append({
            "metric": record["metric"], "route": record["route"], "modifier": record["modifier"],
            "temperature_K": record["temperature_K"], "film_state": record["film_state"],
            "value": record["value"],
            "same_stage1_sequential_cohort1": cohort1_seq.get(key, math.nan),
            "sequential_min": min(seq) if seq else math.nan,
            "sequential_max": max(seq) if seq else math.nan,
            "cosam_min": min(co) if co else math.nan,
            "cosam_max": max(co) if co else math.nan,
            "inside_sequential_range": bool(seq and min(seq) <= record["value"] <= max(seq)),
        })
    return rows


def uncertainty_budget(records, metric: str = "projected_coverage_percent") -> list[dict[str, Any]]:
    """Within-run block SEM versus between-replica SD for one metric."""
    rows = []
    for stat in replica_statistics([r for r in records if r["metric"] == metric]):
        if stat["independent_runs"] >= 2 and math.isfinite(stat["mean_within_run_block_sem"]):
            rows.append({
                **{k: stat[k] for k in ("metric", "route", "modifier", "temperature_K", "film_state")},
                "between_replica_sd": stat["sd"],
                "mean_within_run_block_sem": stat["mean_within_run_block_sem"],
                "ratio": stat["sd"] / stat["mean_within_run_block_sem"]
                if stat["mean_within_run_block_sem"] > 0 else math.inf,
            })
    return rows
