"""Tables and figures for the replica-level route comparison."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from .route_comparison import (
    MODIFIERS,
    gap_seeded_context,
    paired_route_differences,
    paired_vs_reference,
    replica_statistics,
    t_interval,
    uncertainty_budget,
)

ROUTE_ORDER = ("cosam", "sequential", "gap-seeded-I", "gap-seeded-II")
ROUTE_LABELS = {
    "cosam": "CoSAM", "sequential": "Sequential",
    "gap-seeded-I": "Gap-seeded I", "gap-seeded-II": "Gap-seeded II",
    "reference": "Me-4PACz (180)", "dose-control": "Me-4PACz (261)",
}
# categorical slots 1-4 of the reference palette in fixed order; marker shape
# is a second encoding so identity never depends on color alone
ROUTE_STYLE = {
    "cosam": ("#2a78d6", "o"), "sequential": ("#eb6834", "s"),
    "gap-seeded-I": ("#1baf7a", "^"), "gap-seeded-II": ("#eda100", "D"),
    "reference": ("#7a7a74", "v"), "dose-control": ("#4a4a46", "P"),
}
MODIFIER_LABELS = {
    "dcz-4p": "DCZ-4P", "meo-2pacz": "MeO-2PACz", "meo-4padbc": "MeO-4PADBC", "all": "Mean of 3",
}
FIGURE_METRICS = {
    "projected_coverage_percent": "Projected coverage (%)",
    "near_surface_coverage_percent": "Near-surface (≤5 Å) coverage (%)",
    "largest_void_percent": "Largest projected void (% cell)",
    "bound_molecules": "Anchored molecules per cell",
}


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    import csv

    fields: list[str] = []
    for row in rows:
        fields.extend(key for key in row if key not in fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def trajectory_manifest(records) -> list[dict[str, Any]]:
    """One row per analyzed trajectory and source, with its design labels."""
    rows: dict[tuple, dict[str, Any]] = {}
    for record in records:
        if record.get("source") == "derived":
            continue
        key = tuple(record.get(f) for f in ("prepared_root", "system", "stage", "source"))
        row = rows.setdefault(key, {
            field: record.get(field) for field in (
                "prepared_root", "system", "stage", "route", "modifier", "cohort",
                "independent_replicate", "inventory_molecules", "temperature_K", "film_state",
                "frames", "trajectory", "packmol_seed", "velocity_seed", "source",
            )
        })
        for field in ("packmol_seed", "velocity_seed", "trajectory"):
            if row.get(field) in (None, "", "nan") and record.get(field) not in (None, ""):
                row[field] = record.get(field)
    return [rows[key] for key in sorted(rows, key=str)]


def write_route_report(records, output: Path, figure_metrics: dict[str, str] | None = None) -> Path:
    """Write replica, paired-difference, gap-seeded and budget tables + figures."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "records_long.csv", records)
    _write_csv(output / "analyzed_trajectories.csv", trajectory_manifest(records))
    _write_csv(output / "replica_statistics.csv", replica_statistics(records))
    _write_csv(output / "paired_sequential_minus_cosam.csv", paired_route_differences(records))
    _write_csv(output / "gap_seeded_context.csv", gap_seeded_context(records))
    _write_csv(output / "paired_route_minus_reference.csv", paired_vs_reference(records))
    budget = uncertainty_budget(records)
    if budget:
        _write_csv(output / "uncertainty_budget.csv", budget)
    metrics = figure_metrics or {
        key: label for key, label in FIGURE_METRICS.items()
        if any(r["metric"] == key for r in records)
    }
    if metrics:
        plot_route_levels(records, output / "fig_route_levels_300K_relaxed.png", metrics)
        plot_paired_differences(records, output / "fig_paired_differences.png", metrics)
    if budget:
        plot_uncertainty_budget(budget, output / "fig_uncertainty_budget.png")
    return output


def _pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _axes_style(axis) -> None:
    axis.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color("0.55")
    axis.tick_params(colors="0.3", labelsize=7.5)
    axis.grid(axis="y", color="0.9", lw=0.6)
    axis.set_axisbelow(True)


def plot_route_levels(records, path: Path, metrics: dict[str, str],
                      temperature: int = 300, state: str = "relaxed") -> Path:
    plt = _pyplot()
    count = len(metrics)
    columns = 2 if count > 1 else 1
    rows_n = (count + 1) // 2
    figure, axes = plt.subplots(rows_n, columns, figsize=(7.2, 2.8 * rows_n), squeeze=False)
    offsets = {route: (i - 1.5) * 0.18 for i, route in enumerate(ROUTE_ORDER)}
    for axis, (metric, label) in zip(axes.ravel(), metrics.items()):
        subset = [r for r in records if r["metric"] == metric
                  and r["temperature_K"] == temperature and r["film_state"] == state]
        for control, shade in (("reference", "0.86"), ("dose-control", "0.74")):
            values = [r["value"] for r in subset if r["route"] == control]
            if values:
                mean = sum(values) / len(values)
                axis.axhspan(min(values), max(values), color=shade, alpha=0.5, lw=0, zorder=0)
                axis.axhline(mean, color=shade, lw=1.0, zorder=0)
                # label at the band edge so the two controls never collide
                edge, align = (min(values), "top") if control == "reference" else (max(values), "bottom")
                axis.text(2.55, edge, ROUTE_LABELS[control], fontsize=6, color="0.35", va=align)
        for x, modifier in enumerate(MODIFIERS):
            for route in ROUTE_ORDER:
                color, marker = ROUTE_STYLE[route]
                values = [r["value"] for r in subset if r["route"] == route and r["modifier"] == modifier]
                if not values:
                    continue
                position = x + offsets[route]
                if len(values) > 1:
                    stats = t_interval(values)
                    axis.scatter([position] * len(values), values, s=16, facecolors="none",
                                 edgecolors=color, lw=1.0, marker=marker, zorder=3)
                    axis.plot([position, position], [stats["ci_low"], stats["ci_high"]], color=color, lw=1.5, zorder=2)
                    axis.scatter([position], [stats["mean"]], s=34, color=color, marker=marker,
                                 zorder=4, edgecolors="white", linewidths=0.8)
                else:
                    axis.scatter([position], values, s=34, color=color, marker=marker,
                                 zorder=4, edgecolors="white", linewidths=0.8)
        axis.set_xticks(range(len(MODIFIERS)), [MODIFIER_LABELS[m] for m in MODIFIERS])
        axis.set_xlim(-0.6, 3.15)
        axis.set_ylabel(label, fontsize=8, color="0.25")
        _axes_style(axis)
    for axis in axes.ravel()[count:]:
        axis.set_visible(False)
    handles = [plt.Line2D([], [], color=ROUTE_STYLE[r][0], marker=ROUTE_STYLE[r][1], ls="none",
                          label=ROUTE_LABELS[r]) for r in ROUTE_ORDER]
    figure.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=8)
    figure.text(
        0.5, 0.005,
        f"{temperature} K, {state} film. CoSAM/Sequential: open = cohort run, filled = mean, "
        "bar = 95% t-CI (n = 3 cohorts). Gap-seeded: one run each, built from the cohort-1 "
        "stage-1 film. Gray band: Me-4PACz controls (mean and range).",
        ha="center", fontsize=6, color="0.35", wrap=True,
    )
    figure.tight_layout(rect=(0, 0.05, 1, 0.94))
    figure.savefig(path, dpi=220)
    plt.close(figure)
    return path


def plot_paired_differences(records, path: Path, metrics: dict[str, str]) -> Path:
    plt = _pyplot()
    rows = paired_route_differences(records)
    states = [(300, "compressed"), (300, "relaxed"), (400, "compressed"), (400, "relaxed")]
    state_colors = {"compressed": "#7a7a74", "relaxed": "#eb6834"}
    figure, axes = plt.subplots(1, len(metrics), figsize=(2.5 * len(metrics) + 1.6, 4.4),
                                sharey=True, squeeze=False)
    axes = axes[0]
    labels: list[str] = []
    ticks: list[float] = []
    for index, (axis, (metric, label)) in enumerate(zip(axes, metrics.items())):
        y = 0.0
        for modifier in (*MODIFIERS, "all"):
            for temperature, state in states:
                row = next((r for r in rows if r["metric"] == metric and r["modifier"] == modifier
                            and r["temperature_K"] == temperature and r["film_state"] == state), None)
                if row is not None and row["n"] >= 2:
                    color = state_colors[state]
                    marker = "o" if temperature == 300 else "s"
                    axis.plot([row["ci_low"], row["ci_high"]], [y, y], color=color, lw=1.4)
                    cohort_values = [row[k] for k in row if k.startswith("cohort_")]
                    axis.scatter(cohort_values, [y] * len(cohort_values), color=color, s=7, alpha=0.55, zorder=2)
                    axis.scatter([row["mean"]], [y], color=color, marker=marker, s=26, zorder=3,
                                 edgecolors="white", linewidths=0.6)
                if index == 0:
                    ticks.append(y)
                    labels.append(f"{MODIFIER_LABELS[modifier]} · {temperature} K {state}")
                y -= 1.0
            y -= 0.6
        axis.axvline(0, color="0.4", lw=0.8)
        axis.set_xlabel(f"Δ {label}\n(Sequential − CoSAM)", fontsize=7.5, color="0.25")
        _axes_style(axis)
        axis.grid(axis="y", visible=False)
        axis.grid(axis="x", color="0.9", lw=0.6)
    axes[0].set_yticks(ticks, labels, fontsize=6.3)
    figure.text(
        0.5, 0.005,
        "Point = mean within-cohort difference; bar = 95% t-CI (n = 3, df = 2); small dots = cohorts. "
        "Gray = compressed hold, orange = relaxed; circle = 300 K, square = 400 K. "
        "States of one run share a deposited structure and are not independent.",
        ha="center", fontsize=6, color="0.35",
    )
    figure.tight_layout(rect=(0, 0.05, 1, 1))
    figure.savefig(path, dpi=220)
    plt.close(figure)
    return path


def plot_uncertainty_budget(budget, path: Path) -> Path:
    plt = _pyplot()
    figure, axis = plt.subplots(figsize=(4.2, 3.8))
    for route in ("reference", "dose-control", "cosam", "sequential"):
        subset = [r for r in budget if r["route"] == route]
        if not subset:
            continue
        color, marker = ROUTE_STYLE[route]
        axis.scatter([r["mean_within_run_block_sem"] for r in subset],
                     [r["between_replica_sd"] for r in subset],
                     color=color, marker=marker, s=24, label=ROUTE_LABELS[route],
                     edgecolors="white", linewidths=0.6)
    limits = (0.005, 30.0)
    axis.plot(limits, limits, color="0.6", lw=0.8, ls="--")
    axis.text(0.008, 0.012, "1:1", fontsize=7, color="0.5")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(*limits)
    axis.set_ylim(*limits)
    axis.set_xlabel("Within-run block SEM, projected coverage (%)", fontsize=8, color="0.25")
    axis.set_ylabel("Between-replica SD (%)", fontsize=8, color="0.25")
    _axes_style(axis)
    axis.legend(fontsize=7, frameon=False, loc="lower right")
    axis.text(0.006, 20.0, "points above 1:1: replica spread exceeds\nwithin-run SEM", fontsize=6.5, color="0.4", va="top")
    figure.tight_layout()
    figure.savefig(path, dpi=220)
    plt.close(figure)
    return path
