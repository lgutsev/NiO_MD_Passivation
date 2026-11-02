"""Coating continuity, probe accessibility, and anchoring persistence.

This module separates two questions that projected coverage alone mixes:

* collective organization -- which parts of the substrate a finite probe can
  still reach, how large those exposed regions are, and whether they persist;
* molecular attachment -- which phosphonate anchor groups contact the
  canonical exposed-Ni sites, with which denticity, and for how long.

Every height is measured from one consistent surface reference: the
lateral upper envelope of substrate (molecule-0) atom centres, obtained by a
periodic grey dilation of the per-cell maximum substrate height with a disk of
``reference_radius``.  The same reference is used for gating, molecular
heights, and film thickness.

Accessibility model (geometric only; no solvent is simulated):

* ``height_window=None`` (canopy): every ligand atom blocks a probe that
  approaches along the surface normal.  This is the conservative
  line-of-sight model and matches the projected-coverage convention.
* finite ``height_window`` H: only ligand atoms within H of the reference
  block, so overhanging material higher than H is ignored.  This models
  lateral access through an interfacial slab of thickness H.
* probe radius r: the accessible area is the periodic morphological opening
  of the exposed mask by a disk of radius r -- the set of surface points a
  disk of radius r can touch without overlapping any projected ligand van der
  Waals disk.  r = 0 reproduces the projected exposed area.

Frames within one trajectory are time-correlated.  Each scalar series is
reported with its statistical inefficiency and effective sample count; the
unit of inference for route comparisons is the run (independent seed), see
``route_comparison``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import json
import math

from .coverage import (
    VDW_RADII_ANGSTROM,
    _label_dependency,
    _periodic_patch_sizes,
    _resolve_trajectory,
    count_dump_frames,
    iter_dump_frames,
    provenance_fields,
    rasterize_periodic,
)
from .interfacial import (
    _index_by_atom_id,
    _minimum_image,
    _positions_for_ids,
    _surface_tree_for_frame,
    _wrap_xy,
    load_interfacial_topology,
)
from .model_scope import model_assessment


SCHEMA_VERSION = "1"
DEFAULT_PROBE_RADII = (0.0, 1.4, 2.0, 3.0)
DEFAULT_HEIGHT_WINDOWS: tuple[float | None, ...] = (None, 5.0, 8.0, 12.0)
DEFAULT_CONTACT_CUTOFFS = (3.0, 3.25, 3.5)
THICKNESS_THRESHOLDS = (5.0, 10.0, 15.0)


def _np():
    try:
        import numpy as np
        from scipy import ndimage
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise RuntimeError(
            "accessibility analysis requires NumPy and SciPy; "
            "install with: python -m pip install -e '.[analysis]'"
        ) from exc
    return np, ndimage


def window_label(height_window: float | None) -> str:
    return "canopy" if height_window is None else f"h{height_window:g}A"


def setting_key(height_window: float | None, probe_radius: float) -> str:
    return f"{window_label(height_window)}_probe{probe_radius:g}A"


def _cell_indices(np, x, y, bounds, shape):
    ny, nx = shape
    xlo, xhi = bounds[0]
    ylo, yhi = bounds[1]
    ix = np.floor(np.mod(x - xlo, xhi - xlo) / (xhi - xlo) * nx).astype(np.int64) % nx
    iy = np.floor(np.mod(y - ylo, yhi - ylo) / (yhi - ylo) * ny).astype(np.int64) % ny
    return iy, ix


def _disk(np, radius: float, dx: float, dy: float):
    rx = max(0, int(math.ceil(radius / dx)))
    ry = max(0, int(math.ceil(radius / dy)))
    xs = np.arange(-rx, rx + 1, dtype=float) * dx
    ys = np.arange(-ry, ry + 1, dtype=float) * dy
    return ys[:, None] ** 2 + xs[None, :] ** 2 <= radius**2 + 1.0e-12


def reference_surface(x, y, z, bounds, shape, radius: float):
    """Periodic upper envelope of substrate atom centres.

    Each grid cell receives the highest substrate z found within ``radius``
    laterally.  ``radius`` should exceed the largest lateral distance from any
    surface point to its nearest surface atom so that no cell is left
    undefined; 2.5 A does so for rock-salt NiO facets (Ni--O = 2.09 A).
    """
    np, ndimage = _np()
    ny, nx = shape
    iy, ix = _cell_indices(np, x, y, bounds, shape)
    heights = np.full(nx * ny, -np.inf)
    np.maximum.at(heights, iy * nx + ix, z)
    heights = heights.reshape(ny, nx)
    dx = (bounds[0][1] - bounds[0][0]) / nx
    dy = (bounds[1][1] - bounds[1][0]) / ny
    envelope = ndimage.grey_dilation(
        heights, footprint=_disk(np, radius, dx, dy), mode="wrap"
    )
    if not np.all(np.isfinite(envelope)):
        raise ValueError(
            "reference_radius leaves surface cells without a substrate atom; "
            "increase it"
        )
    return envelope


def periodic_opening(mask, radius: float, dx: float, dy: float):
    """Area of ``mask`` reachable by a disk of ``radius`` (periodic in x/y)."""
    np, ndimage = _np()
    if radius <= 0:
        return mask.copy()
    footprint = _disk(np, radius, dx, dy)
    pad_y, pad_x = footprint.shape[0], footprint.shape[1]
    padded = np.pad(mask, ((pad_y, pad_y), (pad_x, pad_x)), mode="wrap")
    centres = ndimage.binary_erosion(padded, structure=footprint, border_value=0)
    touched = ndimage.binary_dilation(centres, structure=footprint)
    return touched[pad_y:-pad_y, pad_x:-pad_x]


def statistical_inefficiency(values) -> float:
    """Statistical inefficiency g = 1 + 2 sum (1 - t/N) C(t).

    The sum stops at the first non-positive autocorrelation.  g = 1 for
    uncorrelated samples; the effective sample count is N / g.
    """
    np, _ = _np()
    series = np.asarray(values, dtype=float)
    series = series[np.isfinite(series)]
    count = len(series)
    if count < 3:
        return math.nan
    centred = series - series.mean()
    variance = float(np.dot(centred, centred) / count)
    if variance <= 0:
        return 1.0
    g = 1.0
    for lag in range(1, count - 1):
        correlation = float(np.dot(centred[:-lag], centred[lag:]) / (count - lag)) / variance
        if correlation <= 0:
            break
        g += 2.0 * correlation * (1.0 - lag / count)
    return max(1.0, g)


def series_statistics(values) -> dict[str, float | int | bool]:
    """Mean, frame SD, early/late halves, and correlation-aware sample size."""
    np, _ = _np()
    series = np.asarray(values, dtype=float)
    series = series[np.isfinite(series)]
    if len(series) == 0:
        return {"mean": math.nan, "frames": 0}
    g = statistical_inefficiency(series)
    n_eff = len(series) / g if math.isfinite(g) else math.nan
    half = len(series) // 2
    early = float(series[:half].mean()) if half else math.nan
    late = float(series[half:].mean()) if half else math.nan
    sd = float(series.std(ddof=1)) if len(series) > 1 else 0.0
    return {
        "mean": float(series.mean()),
        "frame_sd": sd,
        "frames": int(len(series)),
        "statistical_inefficiency": g,
        "effective_samples": n_eff,
        "correlated_sem": (
            sd / math.sqrt(n_eff) if math.isfinite(n_eff) and n_eff > 0 else math.nan
        ),
        "early_half_mean": early,
        "late_half_mean": late,
        "late_minus_early": late - early,
        # fewer than ~5 independent samples: within-run uncertainty is not
        # resolvable from this window and must not be quoted.
        "window_too_short": bool(not math.isfinite(n_eff) or n_eff < 5.0),
    }


def _jaccard(np, left, right) -> float:
    union = np.logical_or(left, right).sum()
    if union == 0:
        return math.nan
    return float(np.logical_and(left, right).sum() / union)


def _random_jaccard(a: float, b: float) -> float:
    """Expected Jaccard index of independent masks with area fractions a, b."""
    denominator = a + b - a * b
    return a * b / denominator if denominator > 0 else math.nan


def _run_lengths(np, column) -> list[int]:
    padded = np.concatenate(([0], column.astype(np.int8), [0]))
    changes = np.flatnonzero(np.diff(padded))
    return (changes[1::2] - changes[::2]).tolist()


def exposure_persistence(stack, lags_frames, persistent_fraction: float = 0.9):
    """Persistence of exposed regions from a (frames, ny, nx) boolean stack."""
    np, _ = _np()
    stack = np.asarray(stack, dtype=bool)
    frames = stack.shape[0]
    flat = stack.reshape(frames, -1)
    probability = flat.mean(axis=0)
    phi = float(flat.mean())
    result: dict[str, Any] = {
        "mean_exposed_fraction": phi,
        "persistently_exposed_fraction": float((probability >= persistent_fraction).mean()),
        "transiently_exposed_fraction": float(
            ((probability > 0) & (probability < persistent_fraction)).mean()
        ),
        "never_exposed_fraction": float((probability == 0).mean()),
        "persistence_threshold": persistent_fraction,
    }
    survival = []
    for lag in lags_frames:
        if lag <= 0 or lag >= frames:
            continue
        exposed_then = flat[:-lag]
        denominator = exposed_then.sum()
        conditional = (
            float((exposed_then & flat[lag:]).sum() / denominator)
            if denominator
            else math.nan
        )
        normalized = (
            (conditional - phi) / (1.0 - phi)
            if math.isfinite(conditional) and phi < 1.0
            else math.nan
        )
        survival.append(
            {"lag_frames": int(lag), "conditional_exposed": conditional, "normalized": normalized}
        )
    result["survival"] = survival
    decayed = [row["lag_frames"] for row in survival if row["normalized"] < math.exp(-1)]
    result["decorrelation_lag_frames"] = decayed[0] if decayed else None
    result["decorrelation_censored"] = not decayed
    half = frames // 2
    if half >= 1:
        early = flat[:half].mean(axis=0) >= persistent_fraction
        late = flat[half:].mean(axis=0) >= persistent_fraction
        result["early_late_persistent_jaccard"] = _jaccard(np, early, late)
        result["early_late_persistent_random_jaccard"] = _random_jaccard(
            float(early.mean()), float(late.mean())
        )
    first, last = flat[0], flat[-1]
    result["first_last_jaccard"] = _jaccard(np, first, last)
    result["first_last_random_jaccard"] = _random_jaccard(
        float(first.mean()), float(last.mean())
    )
    exposed_cells = np.flatnonzero(probability > 0)
    if len(exposed_cells):
        sample = exposed_cells[:: max(1, len(exposed_cells) // 20000)]
        lengths = []
        censored = 0
        for cell in sample:
            column = flat[:, cell]
            runs = _run_lengths(np, column)
            if column[0]:
                censored += 1
            if column[-1]:
                censored += 1
            lengths.extend(runs)
        result["exposure_run_mean_frames"] = float(np.mean(lengths)) if lengths else math.nan
        result["exposure_run_censored_fraction"] = (
            censored / len(lengths) if lengths else math.nan
        )
    return result, probability


def _cosine_switch(np, distances, r_on: float, r_off: float):
    weights = np.zeros_like(distances, dtype=float)
    weights[distances <= r_on] = 1.0
    middle = (distances > r_on) & (distances < r_off)
    weights[middle] = 0.5 * (1.0 + np.cos(math.pi * (distances[middle] - r_on) / (r_off - r_on)))
    return weights


def _histogram(np, values, edges):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    counts, _ = np.histogram(values, bins=edges)
    return counts


def _percentiles(np, values) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"mean": math.nan, "p10": math.nan, "p50": math.nan, "p90": math.nan}
    p10, p50, p90 = np.percentile(values, [10, 50, 90])
    return {"mean": float(values.mean()), "p10": float(p10), "p50": float(p50), "p90": float(p90)}


def analyze_accessibility(
    build_directory: Path,
    trajectory: Path | None = None,
    output: Path | None = None,
    *,
    grid_spacing: float = 0.25,
    probe_radii: tuple[float, ...] = DEFAULT_PROBE_RADII,
    height_windows: tuple[float | None, ...] = DEFAULT_HEIGHT_WINDOWS,
    persistence_settings: tuple[tuple[float | None, float], ...] = ((None, 1.4), (8.0, 1.4)),
    reference_radius: float = 2.5,
    thickness_grid_spacing: float = 1.0,
    contact_cutoffs: tuple[float, ...] = DEFAULT_CONTACT_CUTOFFS,
    primary_contact_cutoff: float = 3.25,
    smooth_contact: tuple[float, float] = (3.0, 3.5),
    persistent_fraction: float = 0.9,
    last_frames: int | None = 100,
    stride: int = 1,
    timestep_fs: float | None = None,
    surface_coordination_cutoff: float = 2.8,
    analysis_profile: str = "classical-ff",
    write_plots: bool = True,
) -> Path:
    """Analyze one trajectory window and return ``accessibility_summary.json``."""
    np, _ = _np()
    if grid_spacing <= 0 or thickness_grid_spacing <= 0 or reference_radius <= 0:
        raise ValueError("grid spacings and reference_radius must be positive")
    if any(radius < 0 for radius in probe_radii):
        raise ValueError("probe radii must be non-negative")
    if any(window is not None and window <= 0 for window in height_windows):
        raise ValueError("height windows must be positive or None")
    if primary_contact_cutoff not in contact_cutoffs:
        contact_cutoffs = tuple(sorted(set(contact_cutoffs) | {primary_contact_cutoff}))
    if not smooth_contact[0] < smooth_contact[1]:
        raise ValueError("smooth_contact must be (r_on, r_off) with r_on < r_off")
    if stride <= 0 or (last_frames is not None and last_frames <= 0):
        raise ValueError("stride and last_frames must be positive")
    for window, radius in persistence_settings:
        if window not in height_windows or radius not in probe_radii:
            raise ValueError(f"persistence setting {(window, radius)} is not analyzed")
    assessment = model_assessment(analysis_profile)

    build_directory = Path(build_directory)
    trajectory_path = _resolve_trajectory(build_directory, trajectory)
    topology = load_interfacial_topology(
        build_directory, surface_coordination_cutoff=surface_coordination_cutoff
    )
    output = Path(output) if output is not None else build_directory / "accessibility-analysis"
    output.mkdir(parents=True, exist_ok=True)

    frame_count = count_dump_frames(trajectory_path)
    if frame_count == 0:
        raise ValueError(f"{trajectory_path}: trajectory contains no frames")
    first_index = max(0, frame_count - last_frames) if last_frames is not None else 0

    elements = topology.atom_elements
    radius_of = {atom_id: VDW_RADII_ANGSTROM[element] for atom_id, element in elements.items()}
    molecules = topology.molecules
    component_names = [component.name for component in topology.components]
    component_of_molecule = np.asarray(
        [component_names.index(molecule.component.name) for molecule in molecules]
    )
    heavy_ids = np.concatenate([np.asarray(m.heavy_atom_ids) for m in molecules])
    heavy_segment = np.concatenate(
        [np.full(len(m.heavy_atom_ids), index) for index, m in enumerate(molecules)]
    )
    group_molecule = []
    oxygen_ids = []
    oxygen_group = []
    phosphorus_ids = []
    for index, molecule in enumerate(molecules):
        for phosphorus, oxygens in zip(molecule.phosphorus_ids, molecule.terminal_oxygen_ids):
            group = len(group_molecule)
            group_molecule.append(index)
            phosphorus_ids.append(phosphorus)
            oxygen_ids.extend(oxygens)
            oxygen_group.extend([group] * len(oxygens))
    group_molecule = np.asarray(group_molecule)
    oxygen_ids = np.asarray(oxygen_ids)
    oxygen_group = np.asarray(oxygen_group)
    oxygen_molecule = group_molecule[oxygen_group]
    first_phosphorus = np.asarray([m.phosphorus_ids[0] for m in molecules])
    core_ids = np.concatenate(
        [np.asarray(m.core_atom_ids or m.heavy_atom_ids) for m in molecules]
    )
    core_segment = np.concatenate(
        [np.full(len(m.core_atom_ids or m.heavy_atom_ids), i) for i, m in enumerate(molecules)]
    )
    core_counts = np.bincount(core_segment, minlength=len(molecules))

    label = _label_dependency()
    shape = thickness_shape = None
    settings = [(window, radius) for window in height_windows for radius in probe_radii]
    probability_sum: dict[str, Any] = {}
    stacks: dict[str, list] = {setting_key(*s): [] for s in persistence_settings}
    gap_areas: dict[str, list[float]] = {setting_key(*s): [] for s in settings}
    thickness_sum = None
    rows: list[dict[str, Any]] = []
    denticity_series = []
    molecule_records: dict[str, list] = {
        "p_height": [], "anchor_o_height": [], "top_height": [], "extension": [],
        "tilt": [], "bound": [], "component": [],
    }
    contact_distance_samples = []
    times_ps = []
    reference_map = None

    for frame in iter_dump_frames(trajectory_path):
        if frame.index < first_index or (frame.index - first_index) % stride:
            continue
        (xlo, xhi), (ylo, yhi) = frame.bounds[0], frame.bounds[1]
        lx, ly = xhi - xlo, yhi - ylo
        if shape is None:
            shape = (max(1, math.ceil(ly / grid_spacing)), max(1, math.ceil(lx / grid_spacing)))
            thickness_shape = (
                max(1, math.ceil(ly / thickness_grid_spacing)),
                max(1, math.ceil(lx / thickness_grid_spacing)),
            )
            for key in (setting_key(*s) for s in settings):
                probability_sum[key] = np.zeros(shape)
            thickness_sum = np.zeros(thickness_shape)
        dx, dy = lx / shape[1], ly / shape[0]
        lookup = _index_by_atom_id(np, frame)
        ligand = frame.molecule_ids > 0
        substrate = ~ligand
        reference_map = reference_surface(
            frame.x[substrate], frame.y[substrate], frame.z[substrate],
            frame.bounds, shape, reference_radius,
        )
        iy, ix = _cell_indices(np, frame.x, frame.y, frame.bounds, shape)
        height = frame.z - reference_map[iy, ix]
        radii = np.asarray([radius_of.get(int(i), 1.7) for i in frame.atom_ids])
        row: dict[str, Any] = {"frame_index": frame.index, "step": frame.step}
        if timestep_fs is not None:
            row["time_ps"] = frame.step * timestep_fs / 1000.0
            times_ps.append(row["time_ps"])

        for window in height_windows:
            blocking = ligand if window is None else ligand & (height <= window)
            covered = rasterize_periodic(
                frame.x[blocking], frame.y[blocking], radii[blocking], frame.bounds, shape
            )
            exposed = ~covered
            for radius in probe_radii:
                key = setting_key(window, radius)
                accessible = periodic_opening(exposed, radius, dx, dy)
                probability_sum[key] += accessible
                sizes = _periodic_patch_sizes(np, label, accessible)
                cell_area = dx * dy
                areas = [size * cell_area for size in sizes]
                gap_areas[key].extend(areas)
                row[f"{key}_accessible_percent"] = float(accessible.mean() * 100.0)
                row[f"{key}_gap_count_per_nm2"] = len(sizes) / (lx * ly / 100.0)
                row[f"{key}_largest_gap_percent"] = (
                    100.0 * max(sizes) / accessible.size if sizes else 0.0
                )
                if key in stacks:
                    stacks[key].append(np.packbits(accessible.ravel()))

        # film thickness: local top of ligand van der Waals spheres above reference
        tiy, tix = _cell_indices(np, frame.x[ligand], frame.y[ligand], frame.bounds, thickness_shape)
        top = np.zeros(thickness_shape[0] * thickness_shape[1])
        np.maximum.at(
            top, tiy * thickness_shape[1] + tix, np.maximum(height[ligand] + radii[ligand], 0.0)
        )
        top = top.reshape(thickness_shape)
        thickness_sum += top
        row["thickness_mean_angstrom"] = float(top.mean())
        covered_cells = top > 0
        row["thickness_covered_mean_angstrom"] = (
            float(top[covered_cells].mean()) if covered_cells.any() else 0.0
        )
        for threshold in THICKNESS_THRESHOLDS:
            row[f"thickness_ge_{threshold:g}A_percent"] = float((top >= threshold).mean() * 100.0)

        # molecular heights, extension, tilt
        heavy_index = lookup[heavy_ids]
        top_height = np.full(len(molecules), -np.inf)
        np.maximum.at(top_height, heavy_segment, height[heavy_index])
        oxygen_index = lookup[oxygen_ids]
        anchor_o = np.full(len(molecules), np.inf)
        np.minimum.at(anchor_o, oxygen_molecule, height[oxygen_index])
        p_height = height[lookup[first_phosphorus]]
        p_position = _positions_for_ids(np, frame, lookup, first_phosphorus)
        core_position = _positions_for_ids(np, frame, lookup, core_ids)
        delta = core_position - p_position[core_segment]
        delta[:, 0] = _minimum_image(delta[:, 0], lx)
        delta[:, 1] = _minimum_image(delta[:, 1], ly)
        vector = np.zeros((len(molecules), 3))
        np.add.at(vector, core_segment, delta)
        vector /= np.maximum(core_counts, 1)[:, None]
        norm = np.linalg.norm(vector, axis=1)
        tilt = np.degrees(np.arccos(np.clip(vector[:, 2] / np.where(norm > 0, norm, 1.0), -1, 1)))

        # head-group contacts with canonical exposed Ni sites
        _, tree3d, _, _, _ = _surface_tree_for_frame(np, frame, lookup, topology)
        oxygen_positions = _wrap_xy(
            np, _positions_for_ids(np, frame, lookup, oxygen_ids), frame.bounds
        )
        distances, _ = tree3d.query(oxygen_positions, k=1)
        contact_distance_samples.append(distances.astype(np.float32))
        for cutoff in contact_cutoffs:
            in_contact = distances <= cutoff
            denticity = np.bincount(oxygen_group, weights=in_contact, minlength=len(group_molecule))
            bound_groups = denticity >= 1
            row[f"bound_group_percent_{cutoff:g}A"] = float(bound_groups.mean() * 100.0)
            if cutoff == primary_contact_cutoff:
                denticity_series.append(denticity.astype(np.int8))
                molecule_bound = np.zeros(len(molecules), dtype=bool)
                np.logical_or.at(molecule_bound, group_molecule, bound_groups)
        smooth = _cosine_switch(np, distances, *smooth_contact)
        row["smooth_contacts_per_group"] = float(smooth.sum() / len(group_molecule))
        row["bound_molecule_percent"] = float(molecule_bound.mean() * 100.0)
        row["bound_molecules_per_nm2"] = float(molecule_bound.sum() / (lx * ly / 100.0))
        for name, values in (
            ("p_height", p_height), ("anchor_o_height", anchor_o), ("top_height", top_height),
            ("extension", top_height - anchor_o), ("tilt", tilt), ("bound", molecule_bound),
            ("component", component_of_molecule),
        ):
            molecule_records[name].append(np.asarray(values))
        rows.append(row)

    if not rows:
        raise ValueError("frame selection produced no trajectory frames")
    frames = len(rows)

    # ---------------- time series ----------------
    fieldnames = list(rows[0])
    with (output / "accessibility_timeseries.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    scalar_statistics = {
        name: series_statistics([row[name] for row in rows])
        for name in fieldnames
        if name not in {"frame_index", "step", "time_ps"}
    }

    frame_interval_ps = (
        float(np.median(np.diff(times_ps))) if len(times_ps) > 1 else None
    )
    lags = sorted({1, 2, 5, 10, 20, 50, 100, 200, 500} & set(range(1, frames)))

    # ---------------- accessibility & gaps ----------------
    maps: dict[str, Any] = {"reference_surface_angstrom": reference_map}
    accessibility: dict[str, Any] = {}
    for window, radius in settings:
        key = setting_key(window, radius)
        areas = np.asarray(gap_areas[key])
        total_area = float(areas.sum())
        maps[f"probability_{key}"] = probability_sum[key] / frames
        maps[f"gap_areas_{key}"] = areas.astype(np.float32)
        accessibility[key] = {
            "height_window_angstrom": window,
            "probe_radius_angstrom": radius,
            "accessible_percent": scalar_statistics[f"{key}_accessible_percent"],
            "gap_count_per_nm2": scalar_statistics[f"{key}_gap_count_per_nm2"],
            "largest_gap_percent": scalar_statistics[f"{key}_largest_gap_percent"],
            "gap_area_angstrom2": {
                **_percentiles(np, areas),
                "area_weighted_mean": (
                    float((areas**2).sum() / total_area) if total_area > 0 else 0.0
                ),
                "fraction_of_exposed_area_in_gaps_ge_100A2": (
                    float(areas[areas >= 100.0].sum() / total_area) if total_area > 0 else 0.0
                ),
            },
        }
    persistence: dict[str, Any] = {}
    for key, packed in stacks.items():
        stack = np.stack([np.unpackbits(p)[: shape[0] * shape[1]] for p in packed]).astype(bool)
        stack = stack.reshape(frames, *shape)
        result, _ = exposure_persistence(stack, lags, persistent_fraction)
        if frame_interval_ps:
            result["frame_interval_ps"] = frame_interval_ps
        persistence[key] = result
    maps["thickness_mean_angstrom"] = thickness_sum / frames

    # ---------------- molecular structure ----------------
    stacked = {name: np.stack(values) for name, values in molecule_records.items()}
    component_index = stacked["component"][0]
    height_edges = np.arange(-5.0, 60.01, 0.5)
    tilt_edges = np.arange(0.0, 180.01, 5.0)
    extension_edges = np.arange(0.0, 40.01, 0.5)
    structure: dict[str, Any] = {}
    for index, name in enumerate(component_names):
        selected = component_index == index
        if not selected.any():
            continue
        bound = stacked["bound"][:, selected]
        entry: dict[str, Any] = {"molecules": int(selected.sum())}
        for field in ("p_height", "anchor_o_height", "top_height", "extension", "tilt"):
            values = stacked[field][:, selected]
            entry[field] = {
                "all": _percentiles(np, values),
                "bound": _percentiles(np, values[bound]),
                "unbound": _percentiles(np, values[~bound]),
            }
            edges = tilt_edges if field == "tilt" else extension_edges if field == "extension" else height_edges
            maps[f"hist_{field}_{name}"] = _histogram(np, values, edges)
            maps[f"hist_{field}_{name}_bound"] = _histogram(np, values[bound], edges)
        structure[name] = entry
    maps["hist_edges_height"] = height_edges
    maps["hist_edges_tilt"] = tilt_edges
    maps["hist_edges_extension"] = extension_edges

    # ---------------- anchoring persistence ----------------
    denticity = np.stack(denticity_series)  # (frames, groups)
    bound = denticity >= 1
    occupancy = bound.mean(axis=0)
    group_component = component_index[group_molecule]
    anchoring: dict[str, Any] = {
        "primary_contact_cutoff_angstrom": primary_contact_cutoff,
        "smooth_contact_r_on_off_angstrom": list(smooth_contact),
        "persistent_fraction": persistent_fraction,
    }
    for index, name in enumerate(component_names):
        selected = group_component == index
        if not selected.any():
            continue
        groups_bound = bound[:, selected]
        groups_dent = denticity[:, selected]
        mode_population = {
            str(mode): float((groups_dent == mode).mean()) for mode in range(4)
        }
        changes = (np.diff(groups_dent, axis=0) != 0).sum()
        bound_frames = groups_bound[:-1].sum()
        same_mode = []
        for lag in lags:
            then = groups_dent[:-lag]
            mask = then >= 1
            same_mode.append(
                {
                    "lag_frames": lag,
                    "same_mode_given_bound": (
                        float((then[mask] == groups_dent[lag:][mask]).mean()) if mask.any() else math.nan
                    ),
                    "still_bound_given_bound": (
                        float(groups_bound[lag:][mask].mean()) if mask.any() else math.nan
                    ),
                }
            )
        half = frames // 2
        anchoring[name] = {
            "anchor_groups": int(selected.sum()),
            "bound_group_fraction": float(groups_bound.mean()),
            "ever_bound_fraction": float((occupancy[selected] > 0).mean()),
            "persistently_bound_fraction": float((occupancy[selected] >= persistent_fraction).mean()),
            "intermittent_fraction": float(
                ((occupancy[selected] > 0) & (occupancy[selected] < persistent_fraction)).mean()
            ),
            "denticity_population": mode_population,
            "mode_changes_per_bound_group_frame": (
                float(changes / bound_frames) if bound_frames else math.nan
            ),
            "early_half_bound_fraction": float(groups_bound[:half].mean()) if half else math.nan,
            "late_half_bound_fraction": float(groups_bound[half:].mean()) if half else math.nan,
            "mode_persistence": same_mode,
        }
    distances = np.concatenate(contact_distance_samples)
    contact_edges = np.arange(1.0, 8.001, 0.05)
    maps["hist_o_ni_nearest"] = _histogram(np, distances, contact_edges)
    maps["hist_edges_o_ni_nearest"] = contact_edges

    np.savez_compressed(output / "accessibility_maps.npz", **maps)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "analysis": "accessibility",
        "build_directory": str(build_directory),
        "trajectory": str(trajectory_path),
        "frames_analyzed": frames,
        "first_frame_index": rows[0]["frame_index"],
        "last_frame_index": rows[-1]["frame_index"],
        "first_step": rows[0]["step"],
        "last_step": rows[-1]["step"],
        "frame_interval_ps": frame_interval_ps,
        "model_assessment": assessment,
        "interpretation_scope": (
            "Geometric accessibility of the substrate to a rigid probe, and "
            "fixed-charge contact persistence. Not a measure of moisture "
            "resistance, electronic passivation, or degradation rate. Under a "
            "classical non-reactive force field, anchoring persistence is the "
            "persistence of P-O...Ni nonbonded contacts within the analyzed "
            "window, not a chemical bond lifetime."
        ),
        "settings": {
            "grid_spacing_angstrom": grid_spacing,
            "probe_radii_angstrom": list(probe_radii),
            "height_windows_angstrom": list(height_windows),
            "persistence_settings": [list(s) for s in persistence_settings],
            "reference": "periodic upper envelope of substrate atom centres",
            "reference_radius_angstrom": reference_radius,
            "thickness_grid_spacing_angstrom": thickness_grid_spacing,
            "contact_cutoffs_angstrom": list(contact_cutoffs),
            "surface_sites": topology.surface_site_source,
            "exposed_ni_sites": len(topology.surface_site_ids),
            "vdw_radii_angstrom": VDW_RADII_ANGSTROM,
            "last_frames": last_frames,
            "stride": stride,
            "timestep_fs": timestep_fs,
        },
        "provenance": provenance_fields(
            trajectory_path=trajectory_path,
            topology_path=build_directory / "topology_output.lmp",
            manifest_path=build_directory / "assembly_manifest.json",
            first_step=rows[0]["step"],
            last_step=rows[-1]["step"],
            first_frame_index=rows[0]["frame_index"],
            last_frame_index=rows[-1]["frame_index"],
            stride=stride,
            blocks=1,
        ),
        "accessibility": accessibility,
        "persistence": persistence,
        "structure": structure,
        "anchoring": anchoring,
        "time_series": scalar_statistics,
    }
    summary_path = output / "accessibility_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=_json_default), encoding="utf-8")
    if write_plots:
        try:
            plot_run(summary_path)
        except ImportError:  # pragma: no cover - matplotlib optional
            pass
    return summary_path


def _json_default(value):
    try:
        import numpy as np

        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
    except ImportError:  # pragma: no cover
        pass
    raise TypeError(f"not JSON serializable: {type(value)!r}")


def plot_run(summary_path: Path) -> Path:
    """Coverage-map, thickness and persistence figure for one analyzed run."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    np, _ = _np()
    summary_path = Path(summary_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    maps = np.load(summary_path.parent / "accessibility_maps.npz")
    keys = [setting_key(*s) for s in summary["settings"]["persistence_settings"]]
    figure, axes = plt.subplots(len(keys) + 2, 1, figsize=(7.0, 2.0 * (len(keys) + 2) + 1.5))
    for axis, key in zip(axes, keys):
        image = axis.imshow(
            maps[f"probability_{key}"], origin="lower", cmap="viridis", vmin=0, vmax=1,
            aspect="auto",
        )
        axis.set_title(f"exposure probability, {key}", fontsize=9)
        axis.set_xticks([])
        axis.set_yticks([])
        figure.colorbar(image, ax=axis, fraction=0.025)
    image = axes[len(keys)].imshow(
        maps["thickness_mean_angstrom"], origin="lower", cmap="magma", aspect="auto"
    )
    axes[len(keys)].set_title("mean film thickness above reference (A)", fontsize=9)
    axes[len(keys)].set_xticks([])
    axes[len(keys)].set_yticks([])
    figure.colorbar(image, ax=axes[len(keys)], fraction=0.025)
    axis = axes[-1]
    for key in keys:
        survival = summary["persistence"][key]["survival"]
        axis.plot(
            [row["lag_frames"] for row in survival],
            [row["normalized"] for row in survival],
            marker="o",
            label=key,
        )
    axis.axhline(math.exp(-1), color="0.6", lw=0.8, ls="--")
    axis.set_xscale("log")
    axis.set_xlabel("lag (frames)")
    axis.set_ylabel("exposure memory c(t)")
    axis.legend(fontsize=7)
    figure.tight_layout()
    path = summary_path.parent / "accessibility_overview.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path
