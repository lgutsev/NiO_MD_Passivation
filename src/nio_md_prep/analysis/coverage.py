"""Coordinate-based projected coverage analysis for LAMMPS trajectories.

The analysis replaces image-color counting with the union of projected
van der Waals disks on a periodic fractional x/y grid.  It streams one dump
frame at a time and uses molecule-ID ranges from ``assembly_manifest.json``
to report total and per-component coverage.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator
import csv
import gzip
import json
import math
import re


ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.060,
    "Cl": 35.450,
    "Ni": 58.6934,
}

# Bondi-style radii where available, with conventional values for the
# remaining elements used by this project.
VDW_RADII_ANGSTROM = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Si": 2.10,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Ni": 1.63,
}

_DATA_SECTIONS = {
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
}

_TRAJECTORY_CANDIDATES = (
    "relax-300K.lammpstrj",
    "relax-400K.lammpstrj",
    "equilibration-300K.lammpstrj",
    "hold-300K.lammpstrj",
    "hold-400K.lammpstrj",
    "deposition.lammpstrj",
    "optimized.lammpstrj",
)


@dataclass(frozen=True)
class Component:
    name: str
    molecule_lo: int
    molecule_hi: int | None
    key: str


@dataclass
class DumpFrame:
    index: int
    step: int
    bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    atom_ids: Any
    molecule_ids: Any
    atom_types: Any
    x: Any
    y: Any
    z: Any
    charges: Any = None


def _dependencies() -> tuple[Any, Any]:
    try:
        import numpy as np
        from scipy.ndimage import grey_dilation
    except ImportError as exc:
        raise RuntimeError(
            "coverage analysis requires NumPy and SciPy; "
            "install with: python -m pip install -e '.[analysis]'"
        ) from exc
    return np, grey_dilation


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _required_line(handle, path: Path, expected: str | None = None) -> str:
    line = handle.readline()
    if not line:
        raise ValueError(f"{path}: unexpected end of LAMMPS dump")
    if expected is not None and not line.startswith(expected):
        raise ValueError(
            f"{path}: expected {expected!r}, found {line.rstrip()!r}"
        )
    return line


def _read_box(handle, path: Path, header: str):
    if any(tilt in header.split()[3:] for tilt in ("xy", "xz", "yz")):
        raise ValueError(
            f"{path}: triclinic dump boxes are not supported by projected coverage"
        )
    bounds = []
    for _ in range(3):
        fields = _required_line(handle, path).split()
        if len(fields) != 2:
            raise ValueError(f"{path}: malformed orthogonal BOX BOUNDS row")
        bounds.append((float(fields[0]), float(fields[1])))
    return tuple(bounds)


def _coordinate_column(columns: list[str], axis: str) -> tuple[int, bool]:
    for name in (axis, axis + "u"):
        if name in columns:
            return columns.index(name), False
    scaled = axis + "s"
    if scaled in columns:
        return columns.index(scaled), True
    raise ValueError(
        f"LAMMPS dump requires {axis}, {axis}u, or {axis}s coordinates"
    )


def iter_dump_frames(path: Path) -> Iterator[DumpFrame]:
    """Stream orthogonal LAMMPS custom-dump frames."""
    np, _ = _dependencies()
    with _open_text(path) as handle:
        frame_index = 0
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                raise ValueError(
                    f"{path}: expected 'ITEM: TIMESTEP', found {line.rstrip()!r}"
                )
            step = int(_required_line(handle, path).strip())
            _required_line(handle, path, "ITEM: NUMBER OF ATOMS")
            atom_count = int(_required_line(handle, path).strip())
            box_header = _required_line(handle, path, "ITEM: BOX BOUNDS")
            bounds = _read_box(handle, path, box_header)
            atom_header = _required_line(handle, path, "ITEM: ATOMS")
            columns = atom_header.split()[2:]
            for required in ("id", "mol", "type"):
                if required not in columns:
                    raise ValueError(
                        f"{path}: dump ATOMS columns are missing {required!r}"
                    )
            x_index, x_scaled = _coordinate_column(columns, "x")
            y_index, y_scaled = _coordinate_column(columns, "y")
            z_index, z_scaled = _coordinate_column(columns, "z")
            rows = [_required_line(handle, path) for _ in range(atom_count)]
            values = np.fromstring("".join(rows), sep=" ")
            expected = atom_count * len(columns)
            if values.size != expected:
                raise ValueError(
                    f"{path}: frame {frame_index} contains {values.size} "
                    f"numeric values; expected {expected}"
                )
            values = values.reshape(atom_count, len(columns))
            x = values[:, x_index]
            y = values[:, y_index]
            if x_scaled:
                xlo, xhi = bounds[0]
                x = xlo + x * (xhi - xlo)
            if y_scaled:
                ylo, yhi = bounds[1]
                y = ylo + y * (yhi - ylo)
            z = values[:, z_index]
            if z_scaled:
                zlo, zhi = bounds[2]
                z = zlo + z * (zhi - zlo)
            charges = values[:, columns.index("q")] if "q" in columns else None
            yield DumpFrame(
                index=frame_index,
                step=step,
                bounds=bounds,
                atom_ids=values[:, columns.index("id")].astype(np.int64),
                molecule_ids=values[:, columns.index("mol")].astype(np.int64),
                atom_types=values[:, columns.index("type")].astype(np.int64),
                x=x,
                y=y,
                z=z,
                charges=charges,
            )
            frame_index += 1


def count_dump_frames(path: Path) -> int:
    """Count frames without retaining atom arrays."""
    count = 0
    with _open_text(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                return count
            if not line.startswith("ITEM: TIMESTEP"):
                raise ValueError(
                    f"{path}: expected 'ITEM: TIMESTEP', found {line.rstrip()!r}"
                )
            _required_line(handle, path)
            _required_line(handle, path, "ITEM: NUMBER OF ATOMS")
            atom_count = int(_required_line(handle, path).strip())
            header = _required_line(handle, path, "ITEM: BOX BOUNDS")
            _read_box(handle, path, header)
            _required_line(handle, path, "ITEM: ATOMS")
            for _ in range(atom_count):
                _required_line(handle, path)
            count += 1


def _unique_key(name: str, used: set[str]) -> str:
    base = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_") or "component"
    key = base
    suffix = 2
    while key in used:
        key = f"{base}_{suffix}"
        suffix += 1
    used.add(key)
    return key


def load_components(manifest_path: Path) -> list[Component]:
    """Load non-surface molecule ranges from an assembly manifest."""
    if not manifest_path.is_file():
        return [Component("ligand", 1, None, "ligand")]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = manifest.get("components", [])
    components = []
    used: set[str] = set()
    for index, row in enumerate(rows, 1):
        molecule_ids = row.get("molecule_ids")
        if not molecule_ids or len(molecule_ids) != 2:
            continue
        lo = max(1, int(molecule_ids[0]))
        hi = int(molecule_ids[1])
        if hi < lo:
            continue
        name = str(row.get("component") or f"component-{index}")
        components.append(Component(name, lo, hi, _unique_key(name, used)))
    if not components:
        return [Component("ligand", 1, None, "ligand")]
    components.sort(key=lambda item: item.molecule_lo)
    for left, right in zip(components, components[1:]):
        if left.molecule_hi is None or left.molecule_hi >= right.molecule_lo:
            raise ValueError(
                f"{manifest_path}: overlapping positive molecule-ID ranges for "
                f"{left.name!r} and {right.name!r}"
            )
    return components


def _element_for_mass(atom_type: int, mass: float) -> str:
    matches = [
        element
        for element, expected in ATOMIC_MASSES.items()
        if abs(expected - mass) <= 0.02
    ]
    if len(matches) != 1:
        raise ValueError(
            f"ambiguous element for atom type {atom_type} with mass {mass}"
        )
    return matches[0]


def load_type_elements(topology_path: Path) -> dict[int, str]:
    """Read the LAMMPS Masses section and identify each atom type."""
    active = False
    masses: dict[int, float] = {}
    for raw in topology_path.read_text(encoding="utf-8").splitlines():
        bare = raw.split("#", 1)[0].strip()
        if bare == "Masses":
            active = True
            continue
        if active and bare in _DATA_SECTIONS:
            break
        if not active or not bare:
            continue
        fields = bare.split()
        if len(fields) >= 2 and fields[0].isdigit():
            masses[int(fields[0])] = float(fields[1])
    if not masses:
        raise ValueError(f"{topology_path}: Masses section is missing or empty")
    return {
        atom_type: _element_for_mass(atom_type, mass)
        for atom_type, mass in masses.items()
    }


def _component_selection(np, component: Component, molecule_ids):
    selected = molecule_ids >= component.molecule_lo
    if component.molecule_hi is not None:
        selected &= molecule_ids <= component.molecule_hi
    return selected


def _disk_footprint(np, radius: float, dx: float, dy: float):
    rx = max(1, int(math.ceil(radius / dx)))
    ry = max(1, int(math.ceil(radius / dy)))
    x = np.arange(-rx, rx + 1, dtype=float) * dx
    y = np.arange(-ry, ry + 1, dtype=float) * dy
    return y[:, None] ** 2 + x[None, :] ** 2 <= radius**2


def rasterize_periodic(
    x,
    y,
    radii,
    bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    shape: tuple[int, int],
):
    """Rasterize the union of periodic projected atomic disks."""
    np, grey_dilation = _dependencies()
    ny, nx = shape
    xlo, xhi = bounds[0]
    ylo, yhi = bounds[1]
    lx, ly = xhi - xlo, yhi - ylo
    if lx <= 0 or ly <= 0:
        raise ValueError("LAMMPS dump has non-positive lateral box dimensions")
    dx, dy = lx / nx, ly / ny
    mask = np.zeros((ny, nx), dtype=bool)
    if len(x) == 0:
        return mask
    wrapped_x = np.mod(x - xlo, lx)
    wrapped_y = np.mod(y - ylo, ly)
    ix = np.floor(wrapped_x / lx * nx).astype(np.int64) % nx
    iy = np.floor(wrapped_y / ly * ny).astype(np.int64) % ny
    for radius in np.unique(radii):
        selected = np.isclose(radii, radius, rtol=0.0, atol=1.0e-12)
        seeds = np.zeros((ny, nx), dtype=np.uint8)
        seeds[iy[selected], ix[selected]] = 1
        footprint = _disk_footprint(np, float(radius), dx, dy)
        mask |= grey_dilation(seeds, footprint=footprint, mode="wrap") > 0
    return mask


def _block_average(
    np, values, requested_blocks: int, suffix: str = "_percent"
) -> dict[str, float | int]:
    array = np.asarray(values, dtype=float)
    block_count = min(max(1, requested_blocks), len(array))
    blocks = [part for part in np.array_split(array, block_count) if len(part)]
    block_means = np.asarray([float(part.mean()) for part in blocks])
    block_sem = (
        float(block_means.std(ddof=1) / math.sqrt(len(block_means)))
        if len(block_means) > 1
        else 0.0
    )
    return {
        f"mean{suffix}": float(array.mean()),
        f"frame_std{suffix}": float(array.std(ddof=1)) if len(array) > 1 else 0.0,
        f"block_sem{suffix}": block_sem,
        "block_count": len(blocks),
    }


def _block_statistics(np, values, requested_blocks: int) -> dict[str, float | int]:
    return _block_average(np, values, requested_blocks, suffix="_percent")


def _label_dependency():
    try:
        from scipy.ndimage import label
    except ImportError as exc:
        raise RuntimeError(
            "void patch analysis requires SciPy; "
            "install with: python -m pip install -e '.[analysis]'"
        ) from exc
    return label


def _periodic_height_map(np, x, y, z, bounds, shape):
    """Max z per periodic x/y grid cell -- a top-of-film topography map.

    Uses every atom in the frame (ligands and substrate), so an uncovered
    patch reads as the bare, lower substrate height rather than a gap in
    the data; this is the height profile a nucleating perovskite layer
    would actually see.
    """
    ny, nx = shape
    xlo, xhi = bounds[0]
    ylo, yhi = bounds[1]
    lx, ly = xhi - xlo, yhi - ylo
    wrapped_x = np.mod(x - xlo, lx)
    wrapped_y = np.mod(y - ylo, ly)
    ix = np.floor(wrapped_x / lx * nx).astype(np.int64) % nx
    iy = np.floor(wrapped_y / ly * ny).astype(np.int64) % ny
    flat_index = iy * nx + ix
    heights = np.full(nx * ny, -np.inf)
    np.maximum.at(heights, flat_index, z)
    return heights.reshape(ny, nx)


def _roughness_statistics(np, height_map) -> dict[str, float]:
    valid = height_map[np.isfinite(height_map)]
    if valid.size == 0:
        return {
            "rms_angstrom": math.nan,
            "mean_absolute_angstrom": math.nan,
            "peak_to_valley_angstrom": math.nan,
        }
    mean_height = float(valid.mean())
    deviations = valid - mean_height
    return {
        "rms_angstrom": float(np.sqrt(np.mean(deviations**2))),
        "mean_absolute_angstrom": float(np.mean(np.abs(deviations))),
        "peak_to_valley_angstrom": float(valid.max() - valid.min()),
    }


def _periodic_patch_sizes(np, label, uncovered) -> list[int]:
    """Cell-count sizes of connected uncovered regions on a periodic x/y grid.

    scipy.ndimage.label does not wrap, so components are labeled normally
    and then merged across the periodic x/y seams with a small union-find.
    """
    structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    labels, num_labels = label(uncovered, structure=structure)
    if num_labels == 0:
        return []
    parent = list(range(num_labels + 1))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    ny, nx = uncovered.shape
    for i in range(ny):
        a, b = int(labels[i, 0]), int(labels[i, nx - 1])
        if a and b and a != b:
            union(a, b)
    for j in range(nx):
        a, b = int(labels[0, j]), int(labels[ny - 1, j])
        if a and b and a != b:
            union(a, b)
    roots = np.asarray([find(i) for i in range(num_labels + 1)])
    flat_labels = roots[labels.ravel()]
    flat_labels = flat_labels[flat_labels > 0]
    counts = np.bincount(flat_labels)
    return counts[counts > 0].tolist()


def _write_plots(
    output: Path,
    rows: list[dict[str, float | int]],
    probability,
    height_map,
    component_keys: list[str],
    timestep_fs: float | None,
) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    if timestep_fs is None:
        x = [int(row["step"]) for row in rows]
        xlabel = "LAMMPS step"
    else:
        x = [float(row["step"]) * timestep_fs / 1000.0 for row in rows]
        xlabel = "Time (ps)"

    figure, axis = plt.subplots(figsize=(7.2, 4.2))
    axis.plot(
        x,
        [row["coverage_total_percent"] for row in rows],
        color="black",
        linewidth=1.8,
        label="total",
    )
    for key in component_keys:
        axis.plot(
            x,
            [row[f"coverage_{key}_percent"] for row in rows],
            linewidth=1.1,
            label=key.replace("_", " "),
        )
    axis.set(xlabel=xlabel, ylabel="Projected coverage (%)", ylim=(0.0, 100.0))
    axis.legend(frameon=False, ncol=2)
    figure.tight_layout()
    timeseries_name = "coverage_vs_time.png"
    figure.savefig(output / timeseries_name, dpi=200)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7.2, 3.2))
    image = axis.imshow(
        probability,
        origin="lower",
        extent=(0.0, 1.0, 0.0, 1.0),
        aspect="auto",
        vmin=0.0,
        vmax=1.0,
        cmap="viridis",
    )
    axis.set(
        xlabel="Fractional x",
        ylabel="Fractional y",
        title="Ligand occupancy probability",
    )
    figure.colorbar(image, ax=axis, label="Occupied-frame fraction")
    figure.tight_layout()
    probability_name = "coverage_probability.png"
    figure.savefig(output / probability_name, dpi=200)
    plt.close(figure)

    import numpy as np

    figure, axis = plt.subplots(figsize=(7.2, 3.2))
    masked_height = np.ma.masked_invalid(height_map)
    image = axis.imshow(
        masked_height,
        origin="lower",
        extent=(0.0, 1.0, 0.0, 1.0),
        aspect="auto",
        cmap="magma",
    )
    axis.set(
        xlabel="Fractional x",
        ylabel="Fractional y",
        title="Top-of-film height map (last analyzed frame)",
    )
    figure.colorbar(image, ax=axis, label="Height (Å)")
    figure.tight_layout()
    height_map_name = "height_map.png"
    figure.savefig(output / height_map_name, dpi=200)
    plt.close(figure)
    return [timeseries_name, probability_name, height_map_name]


def _resolve_trajectory(build_directory: Path, trajectory: Path | None) -> Path:
    if trajectory is not None:
        resolved = trajectory if trajectory.is_absolute() else build_directory / trajectory
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        return resolved
    for name in _TRAJECTORY_CANDIDATES:
        candidate = build_directory / name
        if candidate.is_file():
            return candidate
    expected = ", ".join(_TRAJECTORY_CANDIDATES)
    raise FileNotFoundError(
        f"{build_directory}: no trajectory found; expected one of {expected}"
    )


def analyze_coverage(
    build_directory: Path,
    trajectory: Path | None = None,
    output: Path | None = None,
    *,
    grid_spacing: float = 0.20,
    radius_scale: float = 1.0,
    roughness_grid_spacing: float = 2.0,
    last_frames: int | None = 100,
    stride: int = 1,
    blocks: int = 5,
    exclude_hydrogen: bool = False,
    timestep_fs: float | None = None,
) -> Path:
    """Analyze projected ligand coverage and return the summary path."""
    np, _ = _dependencies()
    build_directory = Path(build_directory)
    if grid_spacing <= 0:
        raise ValueError("grid_spacing must be positive")
    if radius_scale <= 0:
        raise ValueError("radius_scale must be positive")
    if roughness_grid_spacing <= 0:
        raise ValueError("roughness_grid_spacing must be positive")
    if stride <= 0:
        raise ValueError("stride must be positive")
    if blocks <= 0:
        raise ValueError("blocks must be positive")
    if last_frames is not None and last_frames <= 0:
        raise ValueError("last_frames must be positive or None")
    if timestep_fs is not None and timestep_fs <= 0:
        raise ValueError("timestep_fs must be positive")

    trajectory_path = _resolve_trajectory(build_directory, trajectory)
    topology_path = build_directory / "topology_output.lmp"
    if not topology_path.is_file():
        raise FileNotFoundError(topology_path)
    components = load_components(build_directory / "assembly_manifest.json")
    type_elements = load_type_elements(topology_path)
    type_radii = {
        atom_type: VDW_RADII_ANGSTROM[element] * radius_scale
        for atom_type, element in type_elements.items()
    }

    output = Path(output) if output is not None else build_directory / "coverage-analysis"
    output.mkdir(parents=True, exist_ok=True)

    frame_count = count_dump_frames(trajectory_path)
    if frame_count == 0:
        raise ValueError(f"{trajectory_path}: trajectory contains no frames")
    first_index = (
        max(0, frame_count - last_frames) if last_frames is not None else 0
    )

    label = _label_dependency()
    shape = None
    roughness_shape = None
    probability_sum = None
    last_height_map = None
    rows: list[dict[str, float | int]] = []
    analyzed_indices: list[int] = []
    for frame in iter_dump_frames(trajectory_path):
        if frame.index < first_index or (frame.index - first_index) % stride:
            continue
        xlo, xhi = frame.bounds[0]
        ylo, yhi = frame.bounds[1]
        if shape is None:
            nx = max(1, int(math.ceil((xhi - xlo) / grid_spacing)))
            ny = max(1, int(math.ceil((yhi - ylo) / grid_spacing)))
            shape = (ny, nx)
            probability_sum = np.zeros(shape, dtype=np.float64)
            roughness_nx = max(1, int(math.ceil((xhi - xlo) / roughness_grid_spacing)))
            roughness_ny = max(1, int(math.ceil((yhi - ylo) / roughness_grid_spacing)))
            roughness_shape = (roughness_ny, roughness_nx)

        positive = frame.molecule_ids > 0
        mapped = np.zeros(len(frame.molecule_ids), dtype=bool)
        component_masks = []
        for component in components:
            selected = _component_selection(np, component, frame.molecule_ids)
            mapped |= selected
            if exclude_hydrogen:
                selected &= np.asarray(
                    [
                        type_elements.get(int(atom_type)) != "H"
                        for atom_type in frame.atom_types
                    ],
                    dtype=bool,
                )
            selected_types = frame.atom_types[selected]
            try:
                radii = np.asarray(
                    [type_radii[int(atom_type)] for atom_type in selected_types],
                    dtype=float,
                )
            except KeyError as exc:
                raise ValueError(
                    f"{topology_path}: no radius for atom type {exc.args[0]}"
                ) from exc
            component_masks.append(
                rasterize_periodic(
                    frame.x[selected],
                    frame.y[selected],
                    radii,
                    frame.bounds,
                    shape,
                )
            )
        if np.any(positive & ~mapped):
            missing = np.unique(frame.molecule_ids[positive & ~mapped])
            sample = ", ".join(map(str, missing[:10]))
            raise ValueError(
                "assembly manifest does not map all positive molecule IDs; "
                f"first missing IDs: {sample}"
            )

        stack = np.stack(component_masks, axis=0)
        occupancy_count = stack.sum(axis=0)
        total_mask = occupancy_count > 0
        overlap_mask = occupancy_count > 1
        probability_sum += total_mask
        patch_sizes = _periodic_patch_sizes(np, label, ~total_mask)
        total_cells = total_mask.size
        height_map = _periodic_height_map(
            np, frame.x, frame.y, frame.z, frame.bounds, roughness_shape
        )
        roughness = _roughness_statistics(np, height_map)
        last_height_map = height_map
        row: dict[str, float | int] = {
            "frame_index": frame.index,
            "step": frame.step,
            "box_area_angstrom2": (xhi - xlo) * (yhi - ylo),
            "coverage_total_percent": float(total_mask.mean() * 100.0),
            "coverage_uncovered_percent": float((~total_mask).mean() * 100.0),
            "coverage_overlap_percent": float(overlap_mask.mean() * 100.0),
            "void_patch_count": len(patch_sizes),
            "void_largest_patch_percent": (
                100.0 * max(patch_sizes) / total_cells if patch_sizes else 0.0
            ),
            "void_mean_patch_percent": (
                100.0 * (sum(patch_sizes) / len(patch_sizes)) / total_cells
                if patch_sizes
                else 0.0
            ),
            "roughness_rms_angstrom": roughness["rms_angstrom"],
            "roughness_mean_absolute_angstrom": roughness["mean_absolute_angstrom"],
            "roughness_peak_to_valley_angstrom": roughness["peak_to_valley_angstrom"],
        }
        for component, mask in zip(components, component_masks):
            row[f"coverage_{component.key}_percent"] = float(mask.mean() * 100.0)
        rows.append(row)
        analyzed_indices.append(frame.index)

    if not rows or shape is None or probability_sum is None:
        raise ValueError("frame selection produced no trajectory frames")

    fieldnames = list(rows[0])
    timeseries_path = output / "coverage_timeseries.csv"
    with timeseries_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    probability = probability_sum / len(rows)
    np.savez_compressed(
        output / "coverage_probability.npz",
        probability=probability,
        fractional_x=(np.arange(shape[1]) + 0.5) / shape[1],
        fractional_y=(np.arange(shape[0]) + 0.5) / shape[0],
    )
    np.savez_compressed(
        output / "height_map.npz",
        height_angstrom=last_height_map,
        fractional_x=(np.arange(roughness_shape[1]) + 0.5) / roughness_shape[1],
        fractional_y=(np.arange(roughness_shape[0]) + 0.5) / roughness_shape[0],
    )

    metrics = {
        "total": _block_statistics(
            np, [row["coverage_total_percent"] for row in rows], blocks
        ),
        "uncovered": _block_statistics(
            np, [row["coverage_uncovered_percent"] for row in rows], blocks
        ),
        "overlap": _block_statistics(
            np, [row["coverage_overlap_percent"] for row in rows], blocks
        ),
        "void_patch_count": _block_average(
            np, [row["void_patch_count"] for row in rows], blocks, suffix=""
        ),
        "void_largest_patch_percent": _block_statistics(
            np, [row["void_largest_patch_percent"] for row in rows], blocks
        ),
        "void_mean_patch_percent": _block_statistics(
            np, [row["void_mean_patch_percent"] for row in rows], blocks
        ),
        "roughness_rms": _block_average(
            np,
            [row["roughness_rms_angstrom"] for row in rows],
            blocks,
            suffix="_angstrom",
        ),
        "roughness_mean_absolute": _block_average(
            np,
            [row["roughness_mean_absolute_angstrom"] for row in rows],
            blocks,
            suffix="_angstrom",
        ),
        "roughness_peak_to_valley": _block_average(
            np,
            [row["roughness_peak_to_valley_angstrom"] for row in rows],
            blocks,
            suffix="_angstrom",
        ),
    }
    component_metrics = {
        component.name: {
            "key": component.key,
            "molecule_id_range": [
                component.molecule_lo,
                component.molecule_hi,
            ],
            **_block_statistics(
                np,
                [row[f"coverage_{component.key}_percent"] for row in rows],
                blocks,
            ),
        }
        for component in components
    }
    plot_files = _write_plots(
        output,
        rows,
        probability,
        last_height_map,
        [component.key for component in components],
        timestep_fs,
    )
    summary = {
        "method": "periodic_projected_vdw_union",
        "trajectory": str(trajectory_path),
        "topology": str(topology_path),
        "frames_in_trajectory": frame_count,
        "frames_analyzed": len(rows),
        "analyzed_frame_indices": analyzed_indices,
        "grid_target_spacing_angstrom": grid_spacing,
        "grid_shape_yx": list(shape),
        "roughness_grid_target_spacing_angstrom": roughness_grid_spacing,
        "roughness_grid_shape_yx": list(roughness_shape),
        "roughness_reference": (
            "max z per periodic x/y cell over every atom (ligands and "
            "substrate); RMS/mean-absolute/peak-to-valley computed over "
            "the deviation from the frame mean, then block-averaged"
        ),
        "radius_scale": radius_scale,
        "exclude_hydrogen": exclude_hydrogen,
        "vdw_radii_angstrom": VDW_RADII_ANGSTROM,
        "requested_blocks": blocks,
        "metrics": metrics,
        "components": component_metrics,
        "outputs": [
            timeseries_path.name,
            "coverage_probability.npz",
            "height_map.npz",
            *plot_files,
        ],
    }
    summary_path = output / "coverage_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary_path
