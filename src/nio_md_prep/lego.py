"""Coverage-guided initial placement for experimental stage-2 deposition."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import hashlib
import json

from .build import build
from .lammps import parse


def _numpy():
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "lego-style deposition requires NumPy; "
            "install with: python -m pip install -e '.[analysis]'"
        ) from exc
    return np


def _largest_periodic_true_run(mask) -> tuple[int, int]:
    """Return the start and length of the largest circular true run."""
    size = len(mask)
    if size == 0 or not bool(mask.any()):
        raise ValueError("coverage map contains no x interval below the threshold")
    if bool(mask.all()):
        raise ValueError(
            "coverage map marks the entire x direction as empty; "
            "a localized tunnel cannot be identified"
        )

    best_start = 0
    best_length = 0
    current_start = 0
    current_length = 0
    doubled = list(bool(value) for value in mask) * 2
    for index, is_empty in enumerate(doubled):
        if is_empty:
            if current_length == 0:
                current_start = index
            current_length = min(current_length + 1, size)
            if (
                current_length > best_length
                and current_start < size
            ):
                best_start = current_start
                best_length = current_length
        else:
            current_length = 0
    return best_start % size, best_length


def identify_x_tunnel(
    coverage_map: Path,
    primary_final: Path,
    *,
    occupancy_threshold: float = 0.20,
    packing_inset_angstrom: float = 4.0,
    minimum_gap_fraction: float = 0.12,
    minimum_tunnel_width_angstrom: float = 12.0,
    deposition_clearance_above_stage1_angstrom: float = 30.0,
    lateral_confinement: bool = False,
    wall_strength: float = 0.50,
    wall_cutoff_angstrom: float = 3.0,
) -> dict[str, Any]:
    """Identify and center the widest periodic low-occupancy x gap."""
    if not 0.0 <= occupancy_threshold < 1.0:
        raise ValueError("occupancy threshold must be in [0, 1)")
    if not 0.0 < minimum_gap_fraction < 1.0:
        raise ValueError("minimum gap fraction must be in (0, 1)")
    if packing_inset_angstrom <= 0.0:
        raise ValueError("packing inset must be positive")
    if lateral_confinement and packing_inset_angstrom <= wall_cutoff_angstrom:
        raise ValueError(
            "packing inset must exceed the tunnel-wall cutoff so Packmol "
            "starts every atom outside the wall-force zone"
        )
    if minimum_tunnel_width_angstrom <= 0.0:
        raise ValueError("minimum tunnel width must be positive")
    if deposition_clearance_above_stage1_angstrom <= 0.0:
        raise ValueError("deposition clearance above stage 1 must be positive")
    if lateral_confinement and (
        wall_strength <= 0.0 or wall_cutoff_angstrom <= 0.0
    ):
        raise ValueError("tunnel-wall strength and cutoff must be positive")

    np = _numpy()
    coverage_map = Path(coverage_map)
    primary_final = Path(primary_final)
    if not coverage_map.is_file():
        raise FileNotFoundError(coverage_map)
    if not primary_final.is_file():
        raise FileNotFoundError(primary_final)
    with np.load(coverage_map) as archive:
        if "probability" not in archive:
            raise ValueError(f"{coverage_map}: probability array is missing")
        probability = np.asarray(archive["probability"], dtype=float)
    if probability.ndim != 2 or min(probability.shape) < 2:
        raise ValueError(f"{coverage_map}: probability must be a 2D grid")
    if not np.isfinite(probability).all():
        raise ValueError(f"{coverage_map}: probability contains non-finite values")
    if probability.min() < -1.0e-9 or probability.max() > 1.0 + 1.0e-9:
        raise ValueError(f"{coverage_map}: probability lies outside [0, 1]")

    x_profile = probability.mean(axis=0)
    empty_columns = x_profile <= occupancy_threshold
    start_index, column_count = _largest_periodic_true_run(empty_columns)
    nx = probability.shape[1]
    gap_fraction = column_count / nx
    if gap_fraction < minimum_gap_fraction:
        raise ValueError(
            f"widest persistent x gap is only {gap_fraction:.4f} of the cell; "
            f"minimum requested fraction is {minimum_gap_fraction:.4f}"
        )

    primary = parse(primary_final)
    xlo, xhi = primary.bounds["x"]
    lx = xhi - xlo
    if lx <= 0.0:
        raise ValueError(f"{primary_final}: invalid x cell bounds")
    primary_max_z = max(
        float(atom.fields[6]) for atom in primary.sections["Atoms"]
    )
    deposition_endpoint = (
        primary_max_z + deposition_clearance_above_stage1_angstrom
    )
    gap_width = gap_fraction * lx
    packing_width = gap_width - 2.0 * packing_inset_angstrom
    if packing_width < minimum_tunnel_width_angstrom:
        raise ValueError(
            f"detected tunnel is {gap_width:.3f} A wide and leaves only "
            f"{packing_width:.3f} A after the packing inset; requested minimum "
            f"is {minimum_tunnel_width_angstrom:.3f} A"
        )

    original_start = start_index / nx
    original_end_unwrapped = (start_index + column_count) / nx
    original_center = (
        original_start + 0.5 * gap_fraction
    ) % 1.0
    shift_fraction = 0.5 - original_center
    if shift_fraction >= 0.5:
        shift_fraction -= 1.0
    elif shift_fraction < -0.5:
        shift_fraction += 1.0

    wall_lo_fraction = 0.5 - 0.5 * gap_fraction
    wall_hi_fraction = 0.5 + 0.5 * gap_fraction
    inset_fraction = packing_inset_angstrom / lx
    pack_lo_fraction = wall_lo_fraction + inset_fraction
    pack_hi_fraction = wall_hi_fraction - inset_fraction

    return {
        "method": (
            "periodic_x_gap_tunnel"
            if lateral_confinement
            else "periodic_x_gap_seeded"
        ),
        "coverage_map": str(coverage_map.resolve()),
        "coverage_map_sha256": hashlib.sha256(
            coverage_map.read_bytes()
        ).hexdigest(),
        "primary_final": str(primary_final.resolve()),
        "primary_final_sha256": hashlib.sha256(
            primary_final.read_bytes()
        ).hexdigest(),
        "probability_shape_yx": list(map(int, probability.shape)),
        "x_profile_statistic": "mean occupancy probability over fractional y",
        "occupancy_threshold": occupancy_threshold,
        "empty_column_count": int(column_count),
        "gap_start_index": int(start_index),
        "gap_wraps_periodic_x": bool(original_end_unwrapped > 1.0),
        "original_gap_start_fraction": original_start,
        "original_gap_end_fraction": original_end_unwrapped % 1.0,
        "original_gap_center_fraction": original_center,
        "gap_width_fraction": gap_fraction,
        "gap_width_angstrom": gap_width,
        "coordinate_shift_fraction_x": shift_fraction,
        "coordinate_shift_angstrom_x": shift_fraction * lx,
        "shifted_wall_x_fraction": [
            wall_lo_fraction,
            wall_hi_fraction,
        ],
        "shifted_packing_x_fraction": [
            pack_lo_fraction,
            pack_hi_fraction,
        ],
        "packing_inset_angstrom": packing_inset_angstrom,
        "packing_width_angstrom": packing_width,
        "minimum_gap_fraction": minimum_gap_fraction,
        "minimum_tunnel_width_angstrom": minimum_tunnel_width_angstrom,
        "stage1_max_z_angstrom": primary_max_z,
        "deposition_clearance_above_stage1_angstrom": (
            deposition_clearance_above_stage1_angstrom
        ),
        "deposition_wall_endpoint_angstrom": deposition_endpoint,
        "lateral_confinement": bool(lateral_confinement),
        "wall_style": "harmonic" if lateral_confinement else "none",
        "wall_strength_kcal_per_mol_angstrom2": (
            wall_strength if lateral_confinement else None
        ),
        "wall_cutoff_angstrom": (
            wall_cutoff_angstrom if lateral_confinement else None
        ),
        "wall_applies_to": (
            "stage2 atoms only" if lateral_confinement else "none"
        ),
        "wall_lifetime": (
            "minimization and moving-wall deposition only"
            if lateral_confinement
            else "no lateral wall is applied"
        ),
        "bias_warning": (
            "Stage-2 placement is deliberately coverage-guided, but all "
            "subsequent minimization and dynamics are laterally unconstrained. "
            "This remains a seeded control rather than an unbiased model of "
            "spontaneous gap finding."
        ),
    }


def _write_profile_csv(
    coverage_map: Path,
    plan: dict[str, Any],
    destination: Path,
) -> None:
    """Write the detected x profile and selected periodic gap for inspection."""
    np = _numpy()
    with np.load(coverage_map) as archive:
        probability = np.asarray(archive["probability"], dtype=float)
        if "fractional_x" in archive:
            fractional_x = np.asarray(archive["fractional_x"], dtype=float)
        else:
            fractional_x = (
                np.arange(probability.shape[1]) + 0.5
            ) / probability.shape[1]
    if fractional_x.ndim != 1 or len(fractional_x) != probability.shape[1]:
        raise ValueError(
            f"{coverage_map}: fractional_x does not match the probability grid"
        )
    selected = {
        (int(plan["gap_start_index"]) + offset) % probability.shape[1]
        for offset in range(int(plan["empty_column_count"]))
    }
    profile = probability.mean(axis=0)
    shift = float(plan["coordinate_shift_fraction_x"])
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "column_index",
                "original_fractional_x",
                "shifted_fractional_x",
                "mean_occupancy_probability",
                "selected_periodic_gap",
            ]
        )
        for index, (fraction, occupancy) in enumerate(
            zip(fractional_x, profile)
        ):
            writer.writerow(
                [
                    index,
                    f"{float(fraction):.12g}",
                    f"{(float(fraction) + shift) % 1.0:.12g}",
                    f"{float(occupancy):.12g}",
                    int(index in selected),
                ]
            )


def prepare_lego_stage2(
    config_path: Path,
    primary_final: Path,
    coverage_map: Path,
    output: Path,
    *,
    occupancy_threshold: float = 0.20,
    packing_inset_angstrom: float = 4.0,
    minimum_gap_fraction: float = 0.12,
    minimum_tunnel_width_angstrom: float = 12.0,
    deposition_clearance_above_stage1_angstrom: float = 30.0,
    lateral_confinement: bool = False,
    wall_strength: float = 0.50,
    wall_cutoff_angstrom: float = 3.0,
    packed_xyz: Path | None = None,
    packmol_seed: int | None = None,
    velocity_seed: int | None = None,
) -> Path:
    """Build sequential stage 2 with coverage-guided initial x placement."""
    plan = identify_x_tunnel(
        coverage_map,
        primary_final,
        occupancy_threshold=occupancy_threshold,
        packing_inset_angstrom=packing_inset_angstrom,
        minimum_gap_fraction=minimum_gap_fraction,
        minimum_tunnel_width_angstrom=minimum_tunnel_width_angstrom,
        deposition_clearance_above_stage1_angstrom=(
            deposition_clearance_above_stage1_angstrom
        ),
        lateral_confinement=lateral_confinement,
        wall_strength=wall_strength,
        wall_cutoff_angstrom=wall_cutoff_angstrom,
    )
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "lego_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    _write_profile_csv(
        Path(coverage_map),
        plan,
        output / "lego_tunnel_profile.csv",
    )
    result = build(
        config_path,
        output,
        primary_final=primary_final,
        packed_xyz=packed_xyz,
        packmol_seed=packmol_seed,
        velocity_seed=velocity_seed,
        deposition_guidance=plan,
    )
    if (output / "assembly_manifest.json").is_file():
        manifest = json.loads(
            (output / "assembly_manifest.json").read_text(encoding="utf-8")
        )
        manifest["deposition_guidance"] = plan
        (output / "assembly_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )
    return result
