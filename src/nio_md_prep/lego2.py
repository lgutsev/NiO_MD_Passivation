"""Two-dimensional coverage-void guidance for experimental stage-2 deposition."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import hashlib
import json
import math

from .build import build
from .lammps import parse


def _numpy():
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "lego2 deposition requires NumPy; "
            "install with: python -m pip install -e '.[analysis]'"
        ) from exc
    return np


def _load_probability(coverage_map: Path):
    np = _numpy()
    coverage_map = Path(coverage_map)
    if not coverage_map.is_file():
        raise FileNotFoundError(coverage_map)
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
    return np, probability


def _periodic_components(mask) -> list[list[tuple[int, int]]]:
    """Return four-connected true components on a periodic y/x grid."""
    np = _numpy()
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 2:
        raise ValueError("periodic component mask must be two-dimensional")
    ny, nx = mask.shape
    visited = np.zeros(mask.shape, dtype=bool)
    components: list[list[tuple[int, int]]] = []
    for row in range(ny):
        for column in range(nx):
            if not mask[row, column] or visited[row, column]:
                continue
            component: list[tuple[int, int]] = []
            stack = [(row, column)]
            visited[row, column] = True
            while stack:
                current_row, current_column = stack.pop()
                component.append((current_row, current_column))
                for neighbor_row, neighbor_column in (
                    ((current_row - 1) % ny, current_column),
                    ((current_row + 1) % ny, current_column),
                    (current_row, (current_column - 1) % nx),
                    (current_row, (current_column + 1) % nx),
                ):
                    if (
                        mask[neighbor_row, neighbor_column]
                        and not visited[neighbor_row, neighbor_column]
                    ):
                        visited[neighbor_row, neighbor_column] = True
                        stack.append((neighbor_row, neighbor_column))
            components.append(component)
    return components


def _center_shift_cells(indices, size: int) -> tuple[float, int]:
    """Return the circular center fraction and nearest-cell shift to 0.5."""
    np = _numpy()
    fractions = (np.asarray(indices, dtype=float) + 0.5) / size
    vector = np.exp(2.0j * math.pi * fractions).mean()
    if abs(vector) < 1.0e-10:
        raise ValueError(
            "largest void is distributed around the full periodic axis; "
            "a localized two-dimensional target cannot be centered"
        )
    center = float((np.angle(vector) / (2.0 * math.pi)) % 1.0)
    shift_cells = int(round((0.5 - center) * size))
    shift_cells = ((shift_cells + size // 2) % size) - size // 2
    return center, shift_cells


def _largest_true_rectangle(mask) -> tuple[int, int, int, int]:
    """Return top, bottom, left, right edges of the largest true rectangle."""
    np = _numpy()
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 2 or not bool(mask.any()):
        raise ValueError("void component contains no cells")
    ny, nx = mask.shape
    heights = np.zeros(nx, dtype=int)
    best_area = 0
    best = (0, 0, 0, 0)
    for row in range(ny):
        heights = np.where(mask[row], heights + 1, 0)
        stack: list[tuple[int, int]] = []
        for column in range(nx + 1):
            height = int(heights[column]) if column < nx else 0
            start = column
            while stack and stack[-1][1] > height:
                start_index, previous_height = stack.pop()
                area = previous_height * (column - start_index)
                if area > best_area:
                    best_area = area
                    best = (
                        row - previous_height + 1,
                        row + 1,
                        start_index,
                        column,
                    )
                start = start_index
            if not stack or stack[-1][1] < height:
                stack.append((start, height))
    if best_area == 0:
        raise ValueError("largest void contains no rectangular packing region")
    return best


def identify_2d_void(
    coverage_map: Path,
    primary_final: Path,
    *,
    occupancy_threshold: float = 0.20,
    packing_inset_angstrom: float = 4.0,
    minimum_void_area_fraction: float = 0.02,
    minimum_packing_width_angstrom: float = 12.0,
    minimum_packing_height_angstrom: float = 12.0,
    deposition_clearance_above_stage1_angstrom: float = 30.0,
    final_deposition_clearance_above_stage1_angstrom: float = 15.0,
    deposition_continuation_steps: int = 300000,
) -> dict[str, Any]:
    """Find, center, and inscribe a packing rectangle in the largest 2D void."""
    if not 0.0 <= occupancy_threshold < 1.0:
        raise ValueError("occupancy threshold must be in [0, 1)")
    if packing_inset_angstrom <= 0.0:
        raise ValueError("packing inset must be positive")
    if not 0.0 < minimum_void_area_fraction < 1.0:
        raise ValueError("minimum void area fraction must be in (0, 1)")
    if minimum_packing_width_angstrom <= 0.0:
        raise ValueError("minimum packing width must be positive")
    if minimum_packing_height_angstrom <= 0.0:
        raise ValueError("minimum packing height must be positive")
    if deposition_clearance_above_stage1_angstrom <= 0.0:
        raise ValueError("deposition clearance above stage 1 must be positive")
    if not (
        0.0
        < final_deposition_clearance_above_stage1_angstrom
        < deposition_clearance_above_stage1_angstrom
    ):
        raise ValueError(
            "final deposition clearance must be positive and below the "
            "initial deposition clearance"
        )
    if deposition_continuation_steps <= 0:
        raise ValueError("deposition continuation steps must be positive")

    np, probability = _load_probability(Path(coverage_map))
    empty = probability <= occupancy_threshold
    if not bool(empty.any()):
        raise ValueError("coverage map contains no 2D cells below the threshold")
    if bool(empty.all()):
        raise ValueError(
            "coverage map marks the entire surface as empty; "
            "a localized two-dimensional void cannot be identified"
        )
    components = _periodic_components(empty)
    component = max(components, key=len)
    void_fraction = len(component) / empty.size
    if void_fraction < minimum_void_area_fraction:
        raise ValueError(
            f"largest periodic 2D void is only {void_fraction:.4f} of the cell; "
            f"minimum requested fraction is {minimum_void_area_fraction:.4f}"
        )

    component_mask = np.zeros(empty.shape, dtype=bool)
    rows = np.asarray([row for row, _ in component], dtype=int)
    columns = np.asarray([column for _, column in component], dtype=int)
    component_mask[rows, columns] = True
    ny, nx = probability.shape
    center_y, shift_cells_y = _center_shift_cells(rows, ny)
    center_x, shift_cells_x = _center_shift_cells(columns, nx)
    shifted_component = np.roll(
        component_mask,
        shift=(shift_cells_y, shift_cells_x),
        axis=(0, 1),
    )
    shifted_probability = np.roll(
        probability,
        shift=(shift_cells_y, shift_cells_x),
        axis=(0, 1),
    )
    top, bottom, left, right = _largest_true_rectangle(shifted_component)

    primary_final = Path(primary_final)
    if not primary_final.is_file():
        raise FileNotFoundError(primary_final)
    primary = parse(primary_final)
    xlo, xhi = primary.bounds["x"]
    ylo, yhi = primary.bounds["y"]
    lx, ly = xhi - xlo, yhi - ylo
    if lx <= 0.0 or ly <= 0.0:
        raise ValueError(f"{primary_final}: invalid lateral cell bounds")

    rectangle_x_fraction = [left / nx, right / nx]
    rectangle_y_fraction = [top / ny, bottom / ny]
    rectangle_width = (right - left) / nx * lx
    rectangle_height = (bottom - top) / ny * ly
    packing_width = rectangle_width - 2.0 * packing_inset_angstrom
    packing_height = rectangle_height - 2.0 * packing_inset_angstrom
    if packing_width < minimum_packing_width_angstrom:
        raise ValueError(
            f"largest inscribed void rectangle is {rectangle_width:.3f} A "
            f"wide and leaves only {packing_width:.3f} A after the inset; "
            f"requested minimum is {minimum_packing_width_angstrom:.3f} A"
        )
    if packing_height < minimum_packing_height_angstrom:
        raise ValueError(
            f"largest inscribed void rectangle is {rectangle_height:.3f} A "
            f"high in y and leaves only {packing_height:.3f} A after the inset; "
            f"requested minimum is {minimum_packing_height_angstrom:.3f} A"
        )

    inset_x_fraction = packing_inset_angstrom / lx
    inset_y_fraction = packing_inset_angstrom / ly
    packing_x_fraction = [
        rectangle_x_fraction[0] + inset_x_fraction,
        rectangle_x_fraction[1] - inset_x_fraction,
    ]
    packing_y_fraction = [
        rectangle_y_fraction[0] + inset_y_fraction,
        rectangle_y_fraction[1] - inset_y_fraction,
    ]
    primary_max_z = max(
        float(atom.fields[6]) for atom in primary.sections["Atoms"]
    )
    deposition_endpoint = (
        primary_max_z + deposition_clearance_above_stage1_angstrom
    )
    final_deposition_endpoint = (
        primary_max_z + final_deposition_clearance_above_stage1_angstrom
    )
    rectangle_mask = np.zeros(empty.shape, dtype=bool)
    rectangle_mask[top:bottom, left:right] = True

    return {
        "method": "periodic_2d_void_seeded",
        "coverage_map": str(Path(coverage_map).resolve()),
        "coverage_map_sha256": hashlib.sha256(
            Path(coverage_map).read_bytes()
        ).hexdigest(),
        "primary_final": str(primary_final.resolve()),
        "primary_final_sha256": hashlib.sha256(
            primary_final.read_bytes()
        ).hexdigest(),
        "probability_shape_yx": list(map(int, probability.shape)),
        "occupancy_threshold": occupancy_threshold,
        "periodic_connectivity": "four-neighbor in x/y",
        "void_component_count": len(components),
        "largest_void_cell_count": len(component),
        "largest_void_area_fraction": void_fraction,
        "original_void_center_fraction_x": center_x,
        "original_void_center_fraction_y": center_y,
        "coordinate_shift_cells_x": shift_cells_x,
        "coordinate_shift_cells_y": shift_cells_y,
        "coordinate_shift_fraction_x": shift_cells_x / nx,
        "coordinate_shift_fraction_y": shift_cells_y / ny,
        "shifted_void_rectangle_x_fraction": rectangle_x_fraction,
        "shifted_void_rectangle_y_fraction": rectangle_y_fraction,
        "shifted_void_rectangle_cell_edges": {
            "top": top,
            "bottom": bottom,
            "left": left,
            "right": right,
        },
        "shifted_void_rectangle_area_fraction": (
            (right - left) * (bottom - top) / empty.size
        ),
        "shifted_packing_x_fraction": packing_x_fraction,
        "shifted_packing_y_fraction": packing_y_fraction,
        "packing_inset_angstrom": packing_inset_angstrom,
        "void_rectangle_width_angstrom": rectangle_width,
        "void_rectangle_height_angstrom": rectangle_height,
        "packing_width_angstrom": packing_width,
        "packing_height_angstrom": packing_height,
        "minimum_void_area_fraction": minimum_void_area_fraction,
        "minimum_packing_width_angstrom": minimum_packing_width_angstrom,
        "minimum_packing_height_angstrom": minimum_packing_height_angstrom,
        "stage1_max_z_angstrom": primary_max_z,
        "deposition_clearance_above_stage1_angstrom": (
            deposition_clearance_above_stage1_angstrom
        ),
        "deposition_wall_endpoint_angstrom": deposition_endpoint,
        "final_deposition_clearance_above_stage1_angstrom": (
            final_deposition_clearance_above_stage1_angstrom
        ),
        "final_deposition_wall_endpoint_angstrom": final_deposition_endpoint,
        "deposition_continuation_steps": int(deposition_continuation_steps),
        "lateral_confinement": False,
        "wall_style": "none",
        "wall_applies_to": "none",
        "wall_lifetime": "no lateral wall is applied",
        "shifted_probability": shifted_probability,
        "shifted_component_mask": shifted_component,
        "shifted_rectangle_mask": rectangle_mask,
        "bias_warning": (
            "Stage-2 placement is deliberately seeded inside the largest "
            "axis-aligned rectangle fully contained in the largest periodic "
            "2D low-occupancy component. All subsequent minimization and "
            "dynamics are laterally unconstrained. This is a geometric "
            "accessibility control, not unbiased spontaneous gap finding."
        ),
    }


def _write_void_outputs(
    output: Path,
    plan: dict[str, Any],
) -> None:
    np = _numpy()
    probability = np.asarray(plan.pop("shifted_probability"), dtype=float)
    component = np.asarray(plan.pop("shifted_component_mask"), dtype=bool)
    rectangle = np.asarray(plan.pop("shifted_rectangle_mask"), dtype=bool)
    np.savez_compressed(
        output / "lego2_void_map.npz",
        shifted_probability=probability,
        largest_periodic_void=component,
        packing_rectangle=rectangle,
        fractional_x=(np.arange(probability.shape[1]) + 0.5)
        / probability.shape[1],
        fractional_y=(np.arange(probability.shape[0]) + 0.5)
        / probability.shape[0],
    )
    with (output / "lego2_void_map.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "shifted_row_index",
                "shifted_column_index",
                "fractional_x",
                "fractional_y",
                "occupancy_probability",
                "largest_periodic_void",
                "packing_rectangle",
            ]
        )
        ny, nx = probability.shape
        for row in range(ny):
            for column in range(nx):
                writer.writerow(
                    [
                        row,
                        column,
                        f"{(column + 0.5) / nx:.12g}",
                        f"{(row + 0.5) / ny:.12g}",
                        f"{float(probability[row, column]):.12g}",
                        int(component[row, column]),
                        int(rectangle[row, column]),
                    ]
                )


def prepare_lego2_stage2(
    config_path: Path,
    primary_final: Path,
    coverage_map: Path,
    output: Path,
    *,
    occupancy_threshold: float = 0.20,
    packing_inset_angstrom: float = 4.0,
    minimum_void_area_fraction: float = 0.02,
    minimum_packing_width_angstrom: float = 12.0,
    minimum_packing_height_angstrom: float = 12.0,
    deposition_clearance_above_stage1_angstrom: float = 30.0,
    final_deposition_clearance_above_stage1_angstrom: float = 15.0,
    deposition_continuation_steps: int = 300000,
    packed_xyz: Path | None = None,
    packmol_seed: int | None = None,
    velocity_seed: int | None = None,
) -> Path:
    """Build sequential stage 2 inside a localized periodic 2D coverage void."""
    plan = identify_2d_void(
        coverage_map,
        primary_final,
        occupancy_threshold=occupancy_threshold,
        packing_inset_angstrom=packing_inset_angstrom,
        minimum_void_area_fraction=minimum_void_area_fraction,
        minimum_packing_width_angstrom=minimum_packing_width_angstrom,
        minimum_packing_height_angstrom=minimum_packing_height_angstrom,
        deposition_clearance_above_stage1_angstrom=(
            deposition_clearance_above_stage1_angstrom
        ),
        final_deposition_clearance_above_stage1_angstrom=(
            final_deposition_clearance_above_stage1_angstrom
        ),
        deposition_continuation_steps=deposition_continuation_steps,
    )
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    _write_void_outputs(output, plan)
    (output / "lego2_plan.json").write_text(
        json.dumps(plan, indent=2) + "\n",
        encoding="utf-8",
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
    manifest_path = output / "assembly_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["deposition_guidance"] = plan
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )
    return result
