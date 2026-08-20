"""Anchor-resolved structural analysis for NiO/phosphonate trajectories.

This module complements projected coverage with observables that retain the
surface-normal coordinate and phosphonate connectivity.  It reports
component-resolved anchoring, exposed-Ni site ownership, molecular tilt and
orientational order, local-height density profiles, lateral two-dimensional
RDFs, contact persistence, and multi-terminal anchor populations.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator
import csv
import json
import math
import re

from ..config import ROOT
from ..lammps import DataFile, parse
from .coverage import (
    Component,
    count_dump_frames,
    iter_dump_frames,
    load_components,
    load_type_elements,
    provenance_fields,
    _resolve_trajectory,
)
from .model_scope import model_assessment


@dataclass(frozen=True)
class MoleculeDescriptor:
    molecule_id: int
    component: Component
    atom_ids: tuple[int, ...]
    heavy_atom_ids: tuple[int, ...]
    phosphorus_ids: tuple[int, ...]
    terminal_oxygen_ids: tuple[tuple[int, ...], ...]
    core_atom_ids: tuple[int, ...]


@dataclass(frozen=True)
class InterfacialTopology:
    components: tuple[Component, ...]
    molecules: tuple[MoleculeDescriptor, ...]
    surface_site_ids: tuple[int, ...]
    atom_elements: dict[int, str]
    atom_molecules: dict[int, int]
    atom_coordinates: dict[int, tuple[float, float, float]]
    lateral_bounds: tuple[tuple[float, float], tuple[float, float]]
    surface_coordination_cutoff: float
    surface_site_source: str


# Surface-dipole potential-step approximation (Neugebauer-Scheffler-style):
# a uniform parallel dipole sheet of areal density p [e/Angstrom] produces a
# potential step Delta V = p / epsilon0. Derived directly from SI constants
# (not a hardcoded literature number) so the conversion is auditable.
_ELEMENTARY_CHARGE_COULOMB = 1.602176634e-19
_VACUUM_PERMITTIVITY_F_PER_M = 8.8541878128e-12
_ANGSTROM_METERS = 1e-10
_DIPOLE_DENSITY_TO_VOLTS = _ELEMENTARY_CHARGE_COULOMB / (
    _VACUUM_PERMITTIVITY_F_PER_M * _ANGSTROM_METERS
)


def _z_dipole_moment(np, charges, z) -> float:
    """Net z-dipole moment (e*Angstrom) of a charge-neutral atom set.

    Sum(q_i * z_i) is independent of the z origin only when total charge is
    zero, which holds for a fully assembled system (validated elsewhere).
    """
    return float(np.sum(charges * z))


class SiteExchangeTracker:
    """Dwell/gap state machine for cross-component exposed-Ni site hand-offs.

    A site with zero or 2+ simultaneous owners in a frame is a pass-through
    (neither confirms nor ends a prior tenant's occupancy) by itself, but a
    long run of such frames should eventually forget who held the site --
    otherwise ``primary -> long vacancy -> secondary`` reads as one direct
    exchange with no time limit. Two knobs control this:

    - ``min_dwell_frames``: a candidate owner must be the site's sole owner
      for this many consecutive confirmed frames before a hand-off counts
      (default 1: single-frame confirmation, matching prior behavior).
    - ``max_vacancy_gap_frames``: once a site has gone unowned or shared for
      more than this many consecutive frames, its remembered prior owner is
      forgotten, so a later re-occupation starts fresh rather than counting
      as an exchange against a stale owner (default 3).
    """

    def __init__(
        self,
        site_count: int,
        *,
        min_dwell_frames: int = 1,
        max_vacancy_gap_frames: int = 3,
    ) -> None:
        if min_dwell_frames < 1:
            raise ValueError("min_dwell_frames must be >= 1")
        if max_vacancy_gap_frames < 0:
            raise ValueError("max_vacancy_gap_frames must be >= 0")
        self.min_dwell_frames = min_dwell_frames
        self.max_vacancy_gap_frames = max_vacancy_gap_frames
        self.last_distinct_owner: list[str | None] = [None] * site_count
        self._candidate_owner: list[str | None] = [None] * site_count
        self._candidate_streak: list[int] = [0] * site_count
        self._vacancy_streak: list[int] = [0] * site_count
        self.exchange_pair_counts: dict[str, int] = {}
        self.ever_contacted: set[int] = set()
        self._distinct_owners_seen: list[set[str]] = [
            set() for _ in range(site_count)
        ]

    def update(self, site_components: list[set[str]]) -> None:
        """Record one frame's site ownership."""
        for site_index, owners in enumerate(site_components):
            current_owner = next(iter(owners)) if len(owners) == 1 else None
            if current_owner is None:
                self._candidate_owner[site_index] = None
                self._candidate_streak[site_index] = 0
                self._vacancy_streak[site_index] += 1
                if self._vacancy_streak[site_index] > self.max_vacancy_gap_frames:
                    self.last_distinct_owner[site_index] = None
                continue
            self._vacancy_streak[site_index] = 0
            self.ever_contacted.add(site_index)
            self._distinct_owners_seen[site_index].add(current_owner)
            if self._candidate_owner[site_index] == current_owner:
                self._candidate_streak[site_index] += 1
            else:
                self._candidate_owner[site_index] = current_owner
                self._candidate_streak[site_index] = 1
            if self._candidate_streak[site_index] < self.min_dwell_frames:
                continue
            previous_owner = self.last_distinct_owner[site_index]
            if previous_owner is not None and previous_owner != current_owner:
                pair_key = f"{previous_owner}->{current_owner}"
                self.exchange_pair_counts[pair_key] = (
                    self.exchange_pair_counts.get(pair_key, 0) + 1
                )
            self.last_distinct_owner[site_index] = current_owner

    @property
    def contested_site_count(self) -> int:
        """Sites that have had confirmed sole ownership by >=2 distinct components."""
        return sum(
            1 for owners in self._distinct_owners_seen if len(owners) >= 2
        )


def _dependencies() -> tuple[Any, Any, Any]:
    try:
        import numpy as np
        from scipy.ndimage import gaussian_filter1d
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError(
            "interfacial analysis requires NumPy and SciPy; "
            "install with: python -m pip install -e '.[analysis]'"
        ) from exc
    return np, gaussian_filter1d, cKDTree


def _component_for_molecule(
    molecule_id: int, components: list[Component]
) -> Component | None:
    for component in components:
        if molecule_id < component.molecule_lo:
            continue
        if component.molecule_hi is None or molecule_id <= component.molecule_hi:
            return component
    return None


def _adjacency(data: DataFile) -> dict[int, set[int]]:
    graph = {int(row.fields[0]): set() for row in data.sections["Atoms"]}
    for row in data.sections.get("Bonds", []):
        left, right = map(int, row.fields[2:4])
        graph[left].add(right)
        graph[right].add(left)
    return graph


def _graph_distances(
    atom_ids: set[int], sources: tuple[int, ...], graph: dict[int, set[int]]
) -> dict[int, int]:
    distances = {source: 0 for source in sources}
    queue = deque(sources)
    while queue:
        atom_id = queue.popleft()
        for neighbor in graph.get(atom_id, ()):
            if neighbor not in atom_ids or neighbor in distances:
                continue
            distances[neighbor] = distances[atom_id] + 1
            queue.append(neighbor)
    return distances


def _core_atoms(
    atom_ids: set[int],
    phosphorus_ids: tuple[int, ...],
    graph: dict[int, set[int]],
    elements: dict[int, str],
) -> tuple[int, ...]:
    """Select the heavy molecular region farthest in bonds from phosphonate P."""
    candidates = [
        atom_id for atom_id in atom_ids if elements[atom_id] in {"C", "N"}
    ]
    if not candidates:
        candidates = [
            atom_id
            for atom_id in atom_ids
            if elements[atom_id] not in {"H", "O", "P"}
        ]
    if not candidates:
        return phosphorus_ids
    distances = _graph_distances(atom_ids, phosphorus_ids, graph)
    ranked = sorted(
        candidates,
        key=lambda atom_id: (distances.get(atom_id, -1), atom_id),
        reverse=True,
    )
    keep = max(3, int(math.ceil(0.30 * len(ranked))))
    return tuple(sorted(ranked[: min(keep, len(ranked))]))


def _periodic_tree(positions, lx: float, ly: float, *, dimensions: int):
    np, _, cKDTree = _dependencies()
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != dimensions:
        raise ValueError("periodic-tree coordinates have an unexpected shape")
    shifts = []
    for dx in (-lx, 0.0, lx):
        for dy in (-ly, 0.0, ly):
            shift = np.zeros(dimensions, dtype=float)
            shift[0] = dx
            shift[1] = dy
            shifts.append(shift)
    tiled = np.concatenate([positions + shift for shift in shifts], axis=0)
    originals = np.tile(np.arange(len(positions), dtype=np.int64), len(shifts))
    return cKDTree(tiled), originals


def _minimum_image(delta, length: float):
    np, _, _ = _dependencies()
    return delta - length * np.rint(delta / length)


def _static_exposed_ni_sites(
    data: DataFile,
    elements: dict[int, str],
    atom_molecules: dict[int, int],
    coordinates: dict[int, tuple[float, float, float]],
    cutoff: float,
) -> tuple[int, ...]:
    """Identify under-coordinated Ni atoms on the upper half of the slab."""
    np, _, _ = _dependencies()
    ni_ids = [
        atom_id
        for atom_id, element in elements.items()
        if element == "Ni" and atom_molecules[atom_id] == 0
    ]
    oxygen_ids = [
        atom_id
        for atom_id, element in elements.items()
        if element == "O" and atom_molecules[atom_id] == 0
    ]
    if not ni_ids or not oxygen_ids:
        raise ValueError("topology does not contain molecule-0 Ni/O surface atoms")
    ni_positions = np.asarray([coordinates[atom_id] for atom_id in ni_ids])
    oxygen_positions = np.asarray([coordinates[atom_id] for atom_id in oxygen_ids])
    xlo, xhi = data.bounds["x"]
    ylo, yhi = data.bounds["y"]
    tree, _ = _periodic_tree(
        oxygen_positions, xhi - xlo, yhi - ylo, dimensions=3
    )
    coordinations = np.asarray(
        [len(tree.query_ball_point(position, cutoff)) for position in ni_positions],
        dtype=int,
    )
    upper_midpoint = float(np.median(ni_positions[:, 2]))
    selected = [
        atom_id
        for atom_id, position, coordination in zip(
            ni_ids, ni_positions, coordinations
        )
        if position[2] >= upper_midpoint and coordination < 6
    ]
    if not selected:
        top = float(ni_positions[:, 2].max())
        selected = [
            atom_id
            for atom_id, position in zip(ni_ids, ni_positions)
            if position[2] >= top - 3.0
        ]
    if not selected:
        raise ValueError("no exposed upper-surface Ni sites could be identified")
    return tuple(sorted(selected))


def _canonical_exposed_ni_sites(
    data: DataFile,
    atom_elements: dict[int, str],
    atom_molecules: dict[int, int],
    cutoff: float,
) -> tuple[tuple[int, ...], str]:
    """Map one canonical pristine-slab site selection into an assembled system."""
    reference_path = (
        ROOT / "inputs" / "surfaces" / "corrugated-nio-110" / "surface.lmp"
    )
    if not reference_path.is_file():
        raise FileNotFoundError(
            f"canonical corrugated NiO reference is missing: {reference_path}"
        )
    reference = parse(reference_path)
    reference_type_elements = load_type_elements(reference_path)
    reference_rows = reference.sections["Atoms"]
    reference_elements = {
        int(row.fields[0]): reference_type_elements[int(row.fields[2])]
        for row in reference_rows
    }
    reference_molecules = {
        int(row.fields[0]): int(row.fields[1]) for row in reference_rows
    }
    reference_coordinates = {
        int(row.fields[0]): tuple(map(float, row.fields[4:7]))
        for row in reference_rows
    }
    selected_reference_ids = set(
        _static_exposed_ni_sites(
            reference,
            reference_elements,
            reference_molecules,
            reference_coordinates,
            cutoff,
        )
    )
    reference_surface_rows = [
        row for row in reference_rows if int(row.fields[1]) == 0
    ]
    assembled_surface_rows = [
        row for row in data.sections["Atoms"] if int(row.fields[1]) == 0
    ]
    if len(assembled_surface_rows) != len(reference_surface_rows):
        coordinates = {
            int(row.fields[0]): tuple(map(float, row.fields[4:7]))
            for row in data.sections["Atoms"]
        }
        return (
            _static_exposed_ni_sites(
                data,
                atom_elements,
                atom_molecules,
                coordinates,
                cutoff,
            ),
            "topology_fallback_for_noncanonical_surface",
        )
    selected = []
    for reference_row, assembled_row in zip(
        reference_surface_rows, assembled_surface_rows
    ):
        reference_id = int(reference_row.fields[0])
        assembled_id = int(assembled_row.fields[0])
        if reference_elements[reference_id] != atom_elements[assembled_id]:
            raise ValueError(
                "assembled surface atom order does not match the canonical "
                "corrugated NiO reference"
            )
        if atom_molecules[assembled_id] != 0:
            raise ValueError("canonical surface mapping reached a ligand atom")
        if reference_id in selected_reference_ids:
            selected.append(assembled_id)
    if len(selected) != len(selected_reference_ids):
        raise ValueError("canonical exposed-Ni mapping is incomplete")
    return (
        tuple(selected),
        str(reference_path.relative_to(ROOT)),
    )


def load_interfacial_topology(
    build_directory: Path,
    *,
    surface_coordination_cutoff: float = 2.8,
) -> InterfacialTopology:
    """Load molecular anchors, graph-derived cores, and exposed Ni sites."""
    if surface_coordination_cutoff <= 0:
        raise ValueError("surface_coordination_cutoff must be positive")
    build_directory = Path(build_directory)
    topology_path = build_directory / "topology_output.lmp"
    if not topology_path.is_file():
        raise FileNotFoundError(topology_path)
    data = parse(topology_path)
    type_elements = load_type_elements(topology_path)
    atom_rows = {int(row.fields[0]): row for row in data.sections["Atoms"]}
    atom_elements = {
        atom_id: type_elements[int(row.fields[2])]
        for atom_id, row in atom_rows.items()
    }
    atom_molecules = {
        atom_id: int(row.fields[1]) for atom_id, row in atom_rows.items()
    }
    atom_coordinates = {
        atom_id: tuple(map(float, row.fields[4:7]))
        for atom_id, row in atom_rows.items()
    }
    graph = _adjacency(data)
    components = load_components(build_directory / "assembly_manifest.json")
    molecule_atoms: dict[int, set[int]] = {}
    for atom_id, molecule_id in atom_molecules.items():
        if molecule_id > 0:
            molecule_atoms.setdefault(molecule_id, set()).add(atom_id)
    descriptors = []
    for molecule_id, atom_ids in sorted(molecule_atoms.items()):
        component = _component_for_molecule(molecule_id, components)
        if component is None:
            raise ValueError(
                f"assembly manifest does not map molecule ID {molecule_id}"
            )
        phosphorus = tuple(
            sorted(atom_id for atom_id in atom_ids if atom_elements[atom_id] == "P")
        )
        if not phosphorus:
            raise ValueError(
                f"molecule {molecule_id} ({component.name}) has no phosphorus anchor"
            )
        terminal_oxygens = []
        for phosphorus_id in phosphorus:
            oxygens = tuple(
                sorted(
                    atom_id
                    for atom_id in graph[phosphorus_id]
                    if atom_id in atom_ids and atom_elements[atom_id] == "O"
                )
            )
            if len(oxygens) != 3:
                raise ValueError(
                    f"molecule {molecule_id} P atom {phosphorus_id}: "
                    f"expected three phosphonate O neighbors, found {len(oxygens)}"
                )
            terminal_oxygens.append(oxygens)
        descriptors.append(
            MoleculeDescriptor(
                molecule_id=molecule_id,
                component=component,
                atom_ids=tuple(sorted(atom_ids)),
                heavy_atom_ids=tuple(
                    sorted(
                        atom_id
                        for atom_id in atom_ids
                        if atom_elements[atom_id] != "H"
                    )
                ),
                phosphorus_ids=phosphorus,
                terminal_oxygen_ids=tuple(terminal_oxygens),
                core_atom_ids=_core_atoms(
                    atom_ids, phosphorus, graph, atom_elements
                ),
            )
        )
    exposed_sites, surface_site_source = _canonical_exposed_ni_sites(
        data,
        atom_elements,
        atom_molecules,
        surface_coordination_cutoff,
    )
    return InterfacialTopology(
        components=tuple(components),
        molecules=tuple(descriptors),
        surface_site_ids=exposed_sites,
        atom_elements=atom_elements,
        atom_molecules=atom_molecules,
        atom_coordinates=atom_coordinates,
        lateral_bounds=(data.bounds["x"], data.bounds["y"]),
        surface_coordination_cutoff=surface_coordination_cutoff,
        surface_site_source=surface_site_source,
    )


def _selected_frames(
    path: Path, first_index: int, stride: int
) -> Iterator[Any]:
    for frame in iter_dump_frames(path):
        if frame.index >= first_index and (frame.index - first_index) % stride == 0:
            yield frame


def _index_by_atom_id(np, frame) -> Any:
    maximum = int(frame.atom_ids.max())
    lookup = np.full(maximum + 1, -1, dtype=np.int64)
    lookup[frame.atom_ids] = np.arange(len(frame.atom_ids), dtype=np.int64)
    return lookup


def _positions_for_ids(np, frame, lookup, atom_ids) -> Any:
    ids = np.asarray(atom_ids, dtype=np.int64)
    if len(ids) == 0 or int(ids.max()) >= len(lookup):
        raise ValueError("trajectory is missing topology atom IDs")
    indices = lookup[ids]
    if np.any(indices < 0):
        missing = ids[indices < 0]
        raise ValueError(
            "trajectory is missing topology atom IDs: "
            + ", ".join(map(str, missing[:10]))
        )
    return np.column_stack((frame.x[indices], frame.y[indices], frame.z[indices]))


def _wrap_xy(np, positions, bounds):
    result = np.asarray(positions, dtype=float).copy()
    for axis in (0, 1):
        lo, hi = bounds[axis]
        result[:, axis] = lo + np.mod(result[:, axis] - lo, hi - lo)
    return result


def _periodic_centroid(np, positions, lx: float, ly: float):
    positions = np.asarray(positions, dtype=float)
    reference = positions[0].copy()
    unwrapped = positions.copy()
    unwrapped[:, 0] = reference[0] + _minimum_image(
        positions[:, 0] - reference[0], lx
    )
    unwrapped[:, 1] = reference[1] + _minimum_image(
        positions[:, 1] - reference[1], ly
    )
    centroid = unwrapped.mean(axis=0)
    return centroid


def _unwrap_xy(np, positions, lx: float, ly: float):
    positions = np.asarray(positions, dtype=float)
    reference = positions[0].copy()
    unwrapped = positions.copy()
    unwrapped[:, 0] = reference[0] + _minimum_image(
        positions[:, 0] - reference[0], lx
    )
    unwrapped[:, 1] = reference[1] + _minimum_image(
        positions[:, 1] - reference[1], ly
    )
    return unwrapped


def _plane_orientation(np, positions, lx: float, ly: float):
    """Return plane-normal tilt and P2; 0 degrees means a flat molecular plane."""
    positions = _unwrap_xy(np, positions, lx, ly)
    if len(positions) < 3:
        return math.nan, math.nan
    centered = positions - positions.mean(axis=0)
    if np.linalg.matrix_rank(centered) < 2:
        return math.nan, math.nan
    _, _, vectors = np.linalg.svd(centered, full_matrices=False)
    normal = vectors[-1]
    cosine = min(1.0, max(0.0, float(abs(normal[2]))))
    tilt = math.degrees(math.acos(cosine))
    p2 = 0.5 * (3.0 * cosine * cosine - 1.0)
    return tilt, p2


def _anchor_axis_orientation(np, phosphorus_positions, lx: float, ly: float):
    """Return multi-anchor P--P axis tilt and absolute vertical separation."""
    if len(phosphorus_positions) < 2:
        return math.nan, math.nan
    positions = _unwrap_xy(np, phosphorus_positions, lx, ly)
    vector = positions[-1] - positions[0]
    horizontal = float(math.hypot(vector[0], vector[1]))
    vertical = float(abs(vector[2]))
    norm = math.hypot(horizontal, vertical)
    if norm <= 0:
        return math.nan, math.nan
    return math.degrees(math.atan2(horizontal, vertical)), vertical


def _surface_tree_for_frame(np, frame, lookup, topology):
    site_positions = _positions_for_ids(
        np, frame, lookup, topology.surface_site_ids
    )
    site_positions = _wrap_xy(np, site_positions, frame.bounds)
    lx = frame.bounds[0][1] - frame.bounds[0][0]
    ly = frame.bounds[1][1] - frame.bounds[1][0]
    tree3d, original3d = _periodic_tree(
        site_positions, lx, ly, dimensions=3
    )
    tree2d, original2d = _periodic_tree(
        site_positions[:, :2], lx, ly, dimensions=2
    )
    return site_positions, tree3d, original3d, tree2d, original2d


def _nearest_contact_distances(
    trajectory_path: Path,
    topology: InterfacialTopology,
    first_index: int,
    stride: int,
):
    np, _, _ = _dependencies()
    oxygen_ids = sorted(
        {
            atom_id
            for molecule in topology.molecules
            for terminal in molecule.terminal_oxygen_ids
            for atom_id in terminal
        }
    )
    values = []
    for frame in _selected_frames(trajectory_path, first_index, stride):
        lookup = _index_by_atom_id(np, frame)
        _, tree3d, _, _, _ = _surface_tree_for_frame(
            np, frame, lookup, topology
        )
        oxygen_positions = _positions_for_ids(np, frame, lookup, oxygen_ids)
        oxygen_positions = _wrap_xy(np, oxygen_positions, frame.bounds)
        distances, _ = tree3d.query(oxygen_positions, k=1)
        values.extend(map(float, distances))
    return np.asarray(values, dtype=float)


def _infer_contact_cutoff(distances):
    np, gaussian_filter1d, _ = _dependencies()
    edges = np.arange(1.0, 8.0001, 0.05)
    counts, edges = np.histogram(distances, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(counts.astype(float), 1.5)
    peak_region = (centers >= 1.5) & (centers <= 4.5)
    if not np.any(peak_region) or smooth[peak_region].max() <= 0:
        return 4.0, "fallback_no_contact_peak", centers, counts, smooth
    peak_indices = np.flatnonzero(peak_region)
    peak = int(peak_indices[np.argmax(smooth[peak_region])])
    search = np.flatnonzero(
        (centers >= centers[peak] + 0.25)
        & (centers <= min(6.0, centers[peak] + 2.5))
    )
    minimum = None
    for index in search[1:-1]:
        if smooth[index] <= smooth[index - 1] and smooth[index] <= smooth[index + 1]:
            minimum = int(index)
            break
    if minimum is None:
        cutoff = min(5.5, float(centers[peak] + 1.0))
        method = "fallback_peak_plus_1A"
    else:
        cutoff = float(centers[minimum])
        method = "first_smoothed_minimum_after_nearest_Ni_O_peak"
    return cutoff, method, centers, counts, smooth


def _statistics(np, values, blocks: int) -> dict[str, float | int]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) == 0:
        return {
            "mean": math.nan,
            "frame_std": math.nan,
            "block_sem": math.nan,
            "block_count": 0,
        }
    block_count = min(max(1, blocks), len(array))
    partitions = [part for part in np.array_split(array, block_count) if len(part)]
    block_means = np.asarray([part.mean() for part in partitions], dtype=float)
    sem = (
        float(block_means.std(ddof=1) / math.sqrt(len(block_means)))
        if len(block_means) > 1
        else 0.0
    )
    return {
        "mean": float(array.mean()),
        "frame_std": float(array.std(ddof=1)) if len(array) > 1 else 0.0,
        "block_sem": sem,
        "block_count": len(partitions),
    }


def _mean_or_nan(np, values) -> float:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    return float(array.mean()) if len(array) else math.nan


def _nearest_neighbor_mean(np, centers, lx: float, ly: float) -> float:
    centers = np.asarray(centers, dtype=float)
    if len(centers) < 2:
        return math.nan
    delta = centers[:, None, :] - centers[None, :, :]
    delta[:, :, 0] = _minimum_image(delta[:, :, 0], lx)
    delta[:, :, 1] = _minimum_image(delta[:, :, 1], ly)
    distances = np.sqrt((delta**2).sum(axis=2))
    np.fill_diagonal(distances, np.inf)
    return float(distances.min(axis=1).mean())


def _pair_distances_2d(np, left, right, lx, ly, *, same: bool):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if same:
        if len(left) < 2:
            return np.asarray([], dtype=float)
        i, j = np.triu_indices(len(left), k=1)
        delta = left[i] - left[j]
    else:
        if len(left) == 0 or len(right) == 0:
            return np.asarray([], dtype=float)
        delta = (left[:, None, :] - right[None, :, :]).reshape(-1, 2)
    delta[:, 0] = _minimum_image(delta[:, 0], lx)
    delta[:, 1] = _minimum_image(delta[:, 1], ly)
    return np.sqrt((delta**2).sum(axis=1))


def _residence_lengths(sequence: list[bool]) -> list[int]:
    lengths = []
    current = 0
    for value in sequence:
        if value:
            current += 1
        elif current:
            lengths.append(current)
            current = 0
    if current:
        lengths.append(current)
    return lengths


def analyze_interfacial_structure(
    build_directory: Path,
    trajectory: Path | None = None,
    output: Path | None = None,
    *,
    last_frames: int | None = 100,
    stride: int = 1,
    blocks: int = 5,
    contact_cutoff: float | None = 3.25,
    sensitivity_cutoffs: tuple[float, ...] = (3.0, 3.25, 3.5),
    surface_coordination_cutoff: float = 2.8,
    persistence_threshold: float = 0.50,
    z_min: float = -5.0,
    z_max: float = 60.0,
    z_bin_width: float = 0.50,
    rdf_bin_width: float = 0.25,
    rdf_rmax: float = 20.0,
    timestep_fs: float | None = None,
    exchange_min_dwell_frames: int = 1,
    exchange_max_vacancy_gap_frames: int = 3,
    analysis_profile: str = "classical-ff",
) -> Path:
    """Analyze interfacial structure and return ``interface_summary.json``."""
    np, _, _ = _dependencies()
    build_directory = Path(build_directory)
    if stride <= 0 or blocks <= 0:
        raise ValueError("stride and blocks must be positive")
    if last_frames is not None and last_frames <= 0:
        raise ValueError("last_frames must be positive or None")
    if contact_cutoff is not None and contact_cutoff <= 0:
        raise ValueError("contact_cutoff must be positive or None")
    if not sensitivity_cutoffs or any(value <= 0 for value in sensitivity_cutoffs):
        raise ValueError("sensitivity cutoffs must contain positive values")
    if not 0.0 < persistence_threshold <= 1.0:
        raise ValueError("persistence_threshold must be in (0, 1]")
    if not z_min < z_max or z_bin_width <= 0:
        raise ValueError("invalid z-profile range or bin width")
    if rdf_bin_width <= 0 or rdf_rmax <= 0:
        raise ValueError("RDF bin width and maximum must be positive")
    if timestep_fs is not None and timestep_fs <= 0:
        raise ValueError("timestep_fs must be positive")
    if exchange_min_dwell_frames < 1:
        raise ValueError("exchange_min_dwell_frames must be >= 1")
    if exchange_max_vacancy_gap_frames < 0:
        raise ValueError("exchange_max_vacancy_gap_frames must be >= 0")
    assessment = model_assessment(analysis_profile)

    trajectory_path = _resolve_trajectory(build_directory, trajectory)
    topology = load_interfacial_topology(
        build_directory,
        surface_coordination_cutoff=surface_coordination_cutoff,
    )
    output = (
        Path(output)
        if output is not None
        else build_directory / "interface-analysis"
    )
    output.mkdir(parents=True, exist_ok=True)

    frame_count = count_dump_frames(trajectory_path)
    if frame_count == 0:
        raise ValueError(f"{trajectory_path}: trajectory contains no frames")
    first_index = (
        max(0, frame_count - last_frames) if last_frames is not None else 0
    )
    nearest_distances = _nearest_contact_distances(
        trajectory_path, topology, first_index, stride
    )
    if len(nearest_distances) == 0:
        raise ValueError("no phosphonate oxygen distances were collected")
    inferred, cutoff_method, distance_centers, distance_counts, distance_smooth = (
        _infer_contact_cutoff(nearest_distances)
    )
    if contact_cutoff is None:
        contact_cutoff = inferred
    else:
        cutoff_method = "fixed_common_cutoff"
    sensitivity_cutoffs = tuple(
        sorted({float(contact_cutoff), *(float(value) for value in sensitivity_cutoffs)})
    )

    with (output / "contact_distance_histogram.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["distance_bin_center_angstrom", "count", "smoothed_count"]
        )
        writer.writerows(
            zip(distance_centers, distance_counts, distance_smooth)
        )

    molecules_by_component = {
        component.key: [
            molecule
            for molecule in topology.molecules
            if molecule.component.key == component.key
        ]
        for component in topology.components
    }
    bound_sequences = {
        molecule.molecule_id: [] for molecule in topology.molecules
    }
    component_series: dict[str, dict[str, list[float]]] = {
        component.key: {
            "bound_fraction_percent": [],
            "unbound_anchor_density_per_nm2": [],
            "mean_anchored_terminals": [],
            "tilt_degrees": [],
            "orientational_P2": [],
            "plane_normal_tilt_degrees": [],
            "plane_alignment_P2": [],
            "anchor_axis_tilt_degrees": [],
            "anchor_vertical_separation_angstrom": [],
            "phosphorus_height_angstrom": [],
            "core_height_angstrom": [],
            "core_minus_phosphorus_height_angstrom": [],
            "bound_phosphorus_height_angstrom": [],
            "unbound_phosphorus_height_angstrom": [],
            "bound_core_height_angstrom": [],
            "unbound_core_height_angstrom": [],
            "bound_tilt_degrees": [],
            "unbound_tilt_degrees": [],
            "nearest_neighbor_distance_angstrom": [],
        }
        for component in topology.components
    }
    sensitivity_series = {
        component.key: {
            cutoff: {
                "bound_fraction_percent": [],
                "mean_anchored_terminals": [],
            }
            for cutoff in sensitivity_cutoffs
        }
        for component in topology.components
    }
    terminal_series: dict[str, dict[int, list[float]]] = {
        component.key: {} for component in topology.components
    }
    site_series: dict[str, list[float]] = {
        "empty_percent": [],
        "shared_percent": [],
    }
    for component in topology.components:
        site_series[f"{component.key}_only_percent"] = []

    z_edges = np.arange(z_min, z_max + z_bin_width * 0.5, z_bin_width)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    z_counts = {
        component.key: {
            "heavy": np.zeros(len(z_centers), dtype=float),
            "phosphorus": np.zeros(len(z_centers), dtype=float),
        }
        for component in topology.components
    }
    z_volume_sum = np.zeros(len(z_centers), dtype=float)

    rdf_edges = np.arange(0.0, rdf_rmax + rdf_bin_width * 0.5, rdf_bin_width)
    rdf_centers = 0.5 * (rdf_edges[:-1] + rdf_edges[1:])
    rdf_pairs = []
    for left_index, left in enumerate(topology.components):
        for right in topology.components[left_index:]:
            rdf_pairs.append((left, right))
    rdf_actual = {
        f"{left.key}--{right.key}": np.zeros(len(rdf_centers), dtype=float)
        for left, right in rdf_pairs
    }
    rdf_ideal = {
        key: np.zeros(len(rdf_centers), dtype=float) for key in rdf_actual
    }

    rows: list[dict[str, float | int]] = []
    analyzed_indices = []
    analyzed_steps = []
    dipole_moment_series: list[float] = []
    dipole_density_series: list[float] = []
    potential_step_series: list[float] = []
    exchange_tracker = SiteExchangeTracker(
        len(topology.surface_site_ids),
        min_dwell_frames=exchange_min_dwell_frames,
        max_vacancy_gap_frames=exchange_max_vacancy_gap_frames,
    )
    for frame in _selected_frames(trajectory_path, first_index, stride):
        lookup = _index_by_atom_id(np, frame)
        (
            site_positions,
            tree3d,
            tree3d_originals,
            tree2d,
            tree2d_originals,
        ) = _surface_tree_for_frame(np, frame, lookup, topology)
        lx = frame.bounds[0][1] - frame.bounds[0][0]
        ly = frame.bounds[1][1] - frame.bounds[1][0]
        area = lx * ly
        if rdf_rmax > min(lx, ly) / 2.0:
            raise ValueError(
                f"rdf_rmax={rdf_rmax:g} exceeds half the shortest lateral box"
            )

        site_components = [set() for _ in topology.surface_site_ids]
        component_centers: dict[str, list[tuple[float, float]]] = {
            component.key: [] for component in topology.components
        }
        row: dict[str, float | int] = {
            "frame_index": frame.index,
            "step": frame.step,
            "surface_site_count": len(topology.surface_site_ids),
        }
        if frame.charges is not None:
            dipole_moment = _z_dipole_moment(np, frame.charges, frame.z)
            dipole_density = dipole_moment / area if area > 0 else math.nan
            potential_step = dipole_density * _DIPOLE_DENSITY_TO_VOLTS
            dipole_moment_series.append(dipole_moment)
            dipole_density_series.append(dipole_density)
            potential_step_series.append(potential_step)
            row["z_dipole_moment_e_angstrom"] = dipole_moment
            row["z_dipole_areal_density_e_per_angstrom"] = dipole_density
            row["approximate_potential_step_volts"] = potential_step
        for component in topology.components:
            molecule_bound = []
            terminal_counts = []
            sensitivity_terminal_counts = {
                cutoff: [] for cutoff in sensitivity_cutoffs
            }
            tilts = []
            p2_values = []
            plane_tilts = []
            plane_p2_values = []
            anchor_axis_tilts = []
            anchor_vertical_separations = []
            p_heights = []
            core_heights = []
            relative_heights = []
            unbound_anchor_count = 0
            molecules = molecules_by_component[component.key]
            for molecule in molecules:
                terminal_contacts = 0
                terminal_minimum_distances = []
                molecule_sites: set[int] = set()
                for terminal in molecule.terminal_oxygen_ids:
                    positions = _positions_for_ids(np, frame, lookup, terminal)
                    positions = _wrap_xy(np, positions, frame.bounds)
                    nearest_distances, _ = tree3d.query(positions, k=1)
                    terminal_minimum_distances.append(
                        float(np.min(nearest_distances))
                    )
                    neighbors = tree3d.query_ball_point(
                        positions, float(contact_cutoff)
                    )
                    sites = {
                        int(tree3d_originals[index])
                        for atom_neighbors in neighbors
                        for index in atom_neighbors
                    }
                    if sites:
                        terminal_contacts += 1
                        molecule_sites.update(sites)
                    else:
                        unbound_anchor_count += 1
                for site_index in molecule_sites:
                    site_components[site_index].add(component.key)
                bound = terminal_contacts > 0
                bound_sequences[molecule.molecule_id].append(bound)
                molecule_bound.append(bound)
                terminal_counts.append(terminal_contacts)
                for cutoff in sensitivity_cutoffs:
                    sensitivity_terminal_counts[cutoff].append(
                        sum(
                            distance <= cutoff
                            for distance in terminal_minimum_distances
                        )
                    )

                p_positions = _positions_for_ids(
                    np, frame, lookup, molecule.phosphorus_ids
                )
                core_positions = _positions_for_ids(
                    np, frame, lookup, molecule.core_atom_ids
                )
                p_centroid = _periodic_centroid(np, p_positions, lx, ly)
                core_centroid = _periodic_centroid(np, core_positions, lx, ly)
                plane_tilt, plane_p2 = _plane_orientation(
                    np, core_positions, lx, ly
                )
                anchor_axis_tilt, anchor_vertical_separation = (
                    _anchor_axis_orientation(np, p_positions, lx, ly)
                )
                vector = core_centroid - p_centroid
                vector[0] = _minimum_image(vector[0], lx)
                vector[1] = _minimum_image(vector[1], ly)
                horizontal = float(math.hypot(vector[0], vector[1]))
                vertical = float(abs(vector[2]))
                norm = math.hypot(horizontal, vertical)
                tilt = (
                    math.degrees(math.atan2(horizontal, vertical))
                    if norm > 0
                    else math.nan
                )
                cosine = vertical / norm if norm > 0 else math.nan
                p2 = (
                    0.5 * (3.0 * cosine * cosine - 1.0)
                    if math.isfinite(cosine)
                    else math.nan
                )
                p_xy = np.asarray([p_centroid[:2]])
                core_xy = np.asarray([core_centroid[:2]])
                _, nearest_p_tiled = tree2d.query(p_xy, k=1)
                _, nearest_core_tiled = tree2d.query(core_xy, k=1)
                nearest_p_site = int(
                    tree2d_originals[int(nearest_p_tiled[0])]
                )
                nearest_core_site = int(
                    tree2d_originals[int(nearest_core_tiled[0])]
                )
                p_surface_z = float(site_positions[nearest_p_site, 2])
                core_surface_z = float(site_positions[nearest_core_site, 2])
                tilts.append(tilt)
                p2_values.append(p2)
                plane_tilts.append(plane_tilt)
                plane_p2_values.append(plane_p2)
                anchor_axis_tilts.append(anchor_axis_tilt)
                anchor_vertical_separations.append(anchor_vertical_separation)
                p_height = float(p_centroid[2] - p_surface_z)
                core_height = float(core_centroid[2] - core_surface_z)
                p_heights.append(p_height)
                core_heights.append(core_height)
                relative_heights.append(float(core_centroid[2] - p_centroid[2]))
                component_centers[component.key].append(
                    (
                        float(
                            frame.bounds[0][0]
                            + np.mod(
                                core_centroid[0] - frame.bounds[0][0], lx
                            )
                        ),
                        float(
                            frame.bounds[1][0]
                            + np.mod(
                                core_centroid[1] - frame.bounds[1][0], ly
                            )
                        ),
                    )
                )

            bound_percent = (
                100.0 * sum(molecule_bound) / len(molecule_bound)
                if molecule_bound
                else math.nan
            )
            unbound_anchor_density_per_nm2 = (
                100.0 * unbound_anchor_count / area if area > 0 else math.nan
            )
            component_values = {
                "bound_fraction_percent": bound_percent,
                "unbound_anchor_density_per_nm2": unbound_anchor_density_per_nm2,
                "mean_anchored_terminals": _mean_or_nan(np, terminal_counts),
                "tilt_degrees": _mean_or_nan(np, tilts),
                "orientational_P2": _mean_or_nan(np, p2_values),
                "plane_normal_tilt_degrees": _mean_or_nan(np, plane_tilts),
                "plane_alignment_P2": _mean_or_nan(np, plane_p2_values),
                "anchor_axis_tilt_degrees": _mean_or_nan(
                    np, anchor_axis_tilts
                ),
                "anchor_vertical_separation_angstrom": _mean_or_nan(
                    np, anchor_vertical_separations
                ),
                "phosphorus_height_angstrom": _mean_or_nan(np, p_heights),
                "core_height_angstrom": _mean_or_nan(np, core_heights),
                "core_minus_phosphorus_height_angstrom": _mean_or_nan(
                    np, relative_heights
                ),
                "bound_phosphorus_height_angstrom": _mean_or_nan(
                    np,
                    [
                        value
                        for value, bound in zip(p_heights, molecule_bound)
                        if bound
                    ],
                ),
                "unbound_phosphorus_height_angstrom": _mean_or_nan(
                    np,
                    [
                        value
                        for value, bound in zip(p_heights, molecule_bound)
                        if not bound
                    ],
                ),
                "bound_core_height_angstrom": _mean_or_nan(
                    np,
                    [
                        value
                        for value, bound in zip(core_heights, molecule_bound)
                        if bound
                    ],
                ),
                "unbound_core_height_angstrom": _mean_or_nan(
                    np,
                    [
                        value
                        for value, bound in zip(core_heights, molecule_bound)
                        if not bound
                    ],
                ),
                "bound_tilt_degrees": _mean_or_nan(
                    np,
                    [
                        value
                        for value, bound in zip(tilts, molecule_bound)
                        if bound
                    ],
                ),
                "unbound_tilt_degrees": _mean_or_nan(
                    np,
                    [
                        value
                        for value, bound in zip(tilts, molecule_bound)
                        if not bound
                    ],
                ),
                "nearest_neighbor_distance_angstrom": _nearest_neighbor_mean(
                    np, component_centers[component.key], lx, ly
                ),
            }
            for metric, value in component_values.items():
                component_series[component.key][metric].append(value)
                row[f"{component.key}_{metric}"] = value
            for cutoff, counts in sensitivity_terminal_counts.items():
                sensitivity_series[component.key][cutoff][
                    "bound_fraction_percent"
                ].append(
                    100.0 * sum(value > 0 for value in counts) / len(counts)
                    if counts
                    else math.nan
                )
                sensitivity_series[component.key][cutoff][
                    "mean_anchored_terminals"
                ].append(_mean_or_nan(np, counts))
            maximum_terminals = max(
                (len(molecule.terminal_oxygen_ids) for molecule in molecules),
                default=0,
            )
            for count in range(maximum_terminals + 1):
                population = (
                    100.0
                    * sum(value == count for value in terminal_counts)
                    / len(terminal_counts)
                    if terminal_counts
                    else math.nan
                )
                terminal_series[component.key].setdefault(count, []).append(
                    population
                )
                row[f"{component.key}_terminal_{count}_percent"] = population

            heavy_ids = [
                atom_id
                for molecule in molecules
                for atom_id in molecule.heavy_atom_ids
            ]
            phosphorus_ids = [
                atom_id
                for molecule in molecules
                for atom_id in molecule.phosphorus_ids
            ]
            for group, atom_ids in (
                ("heavy", heavy_ids),
                ("phosphorus", phosphorus_ids),
            ):
                positions = _positions_for_ids(np, frame, lookup, atom_ids)
                positions = _wrap_xy(np, positions, frame.bounds)
                _, nearest_tiled = tree2d.query(positions[:, :2], k=1)
                nearest_sites = tree2d_originals[nearest_tiled]
                heights = positions[:, 2] - site_positions[nearest_sites, 2]
                z_counts[component.key][group] += np.histogram(
                    heights, bins=z_edges
                )[0]

        empty = sum(not owners for owners in site_components)
        shared = sum(len(owners) > 1 for owners in site_components)
        site_series["empty_percent"].append(
            100.0 * empty / len(site_components)
        )
        site_series["shared_percent"].append(
            100.0 * shared / len(site_components)
        )
        row["site_empty_percent"] = site_series["empty_percent"][-1]
        row["site_shared_percent"] = site_series["shared_percent"][-1]
        for component in topology.components:
            only = sum(owners == {component.key} for owners in site_components)
            value = 100.0 * only / len(site_components)
            site_series[f"{component.key}_only_percent"].append(value)
            row[f"site_{component.key}_only_percent"] = value

        exchange_tracker.update(site_components)
        row["site_competition_cumulative_exchange_events"] = sum(
            exchange_tracker.exchange_pair_counts.values()
        )

        shell_areas = math.pi * (rdf_edges[1:] ** 2 - rdf_edges[:-1] ** 2)
        for left, right in rdf_pairs:
            key = f"{left.key}--{right.key}"
            left_centers = component_centers[left.key]
            right_centers = component_centers[right.key]
            same = left.key == right.key
            distances = _pair_distances_2d(
                np, left_centers, right_centers, lx, ly, same=same
            )
            rdf_actual[key] += np.histogram(distances, bins=rdf_edges)[0]
            if same:
                pair_count = len(left_centers) * (len(left_centers) - 1) / 2
            else:
                pair_count = len(left_centers) * len(right_centers)
            rdf_ideal[key] += pair_count * shell_areas / area

        z_volume_sum += area * z_bin_width
        rows.append(row)
        analyzed_indices.append(frame.index)
        analyzed_steps.append(frame.step)

    if not rows:
        raise ValueError("frame selection produced no trajectory frames")

    with (output / "interface_timeseries.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    z_rows = []
    for component in topology.components:
        heavy_density = np.divide(
            z_counts[component.key]["heavy"],
            z_volume_sum,
            out=np.zeros_like(z_volume_sum),
            where=z_volume_sum > 0,
        )
        p_density = np.divide(
            z_counts[component.key]["phosphorus"],
            z_volume_sum,
            out=np.zeros_like(z_volume_sum),
            where=z_volume_sum > 0,
        )
        for height, heavy, phosphorus in zip(
            z_centers, heavy_density, p_density
        ):
            z_rows.append(
                {
                    "component": component.name,
                    "component_key": component.key,
                    "local_height_angstrom": float(height),
                    "heavy_atom_number_density_angstrom-3": float(heavy),
                    "phosphorus_number_density_angstrom-3": float(phosphorus),
                }
            )
    with (output / "z_density_profiles.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(z_rows[0]))
        writer.writeheader()
        writer.writerows(z_rows)

    rdf_rows = []
    for key in rdf_actual:
        g = np.divide(
            rdf_actual[key],
            rdf_ideal[key],
            out=np.full_like(rdf_actual[key], np.nan),
            where=rdf_ideal[key] > 0,
        )
        for radius, actual, ideal, value in zip(
            rdf_centers, rdf_actual[key], rdf_ideal[key], g
        ):
            rdf_rows.append(
                {
                    "component_pair": key,
                    "radius_angstrom": float(radius),
                    "g_2d": float(value),
                    "observed_pairs": float(actual),
                    "ideal_pairs": float(ideal),
                }
            )
    with (output / "lateral_rdf.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rdf_rows[0]))
        writer.writeheader()
        writer.writerows(rdf_rows)

    step_spacing = (
        float(np.median(np.diff(analyzed_steps))) if len(analyzed_steps) > 1 else 0.0
    )
    frame_spacing_ps = (
        step_spacing * timestep_fs / 1000.0
        if timestep_fs is not None
        else None
    )
    analyzed_step_span = (
        int(analyzed_steps[-1] - analyzed_steps[0])
        if len(analyzed_steps) > 1
        else 0
    )
    analyzed_time_ps = (
        analyzed_step_span * timestep_fs / 1000.0
        if timestep_fs is not None
        else None
    )
    exchange_event_count = sum(exchange_tracker.exchange_pair_counts.values())
    ever_contacted_site_count = len(exchange_tracker.ever_contacted)
    contested_site_count = exchange_tracker.contested_site_count
    analyzed_ns = (
        analyzed_time_ps / 1000.0
        if analyzed_time_ps is not None and analyzed_time_ps > 0.0
        else None
    )

    def _rate(denominator: int | None) -> float | None:
        if analyzed_ns is None or not denominator:
            return None
        return exchange_event_count / denominator / analyzed_ns

    exchange_rate_per_site_per_ns = _rate(len(topology.surface_site_ids))
    exchange_rate_per_ever_contacted_site_per_ns = _rate(ever_contacted_site_count)
    exchange_rate_per_contested_site_per_ns = _rate(contested_site_count)
    component_summaries = {}
    for component in topology.components:
        molecules = molecules_by_component[component.key]
        occupancies = [
            sum(bound_sequences[molecule.molecule_id])
            / len(bound_sequences[molecule.molecule_id])
            for molecule in molecules
        ]
        residence_frames = [
            length
            for molecule in molecules
            for length in _residence_lengths(
                bound_sequences[molecule.molecule_id]
            )
        ]
        terminal_counts = sorted(
            {len(molecule.terminal_oxygen_ids) for molecule in molecules}
        )
        component_summaries[component.name] = {
            "key": component.key,
            "molecule_id_range": [
                component.molecule_lo,
                component.molecule_hi,
            ],
            "molecule_count": len(molecules),
            "anchor_terminals_per_molecule": terminal_counts,
            "persistent_bound_molecules_percent": 100.0
            * sum(value >= persistence_threshold for value in occupancies)
            / len(occupancies),
            "mean_molecule_contact_occupancy_percent": 100.0
            * float(np.mean(occupancies)),
            "residence_event_count": len(residence_frames),
            "mean_residence_frames": (
                float(np.mean(residence_frames)) if residence_frames else 0.0
            ),
            "max_residence_frames": max(residence_frames, default=0),
            "mean_residence_ps": (
                float(np.mean(residence_frames)) * frame_spacing_ps
                if residence_frames and frame_spacing_ps is not None
                else None
            ),
            "max_residence_ps": (
                max(residence_frames, default=0) * frame_spacing_ps
                if frame_spacing_ps is not None
                else None
            ),
            "metrics": {
                metric: _statistics(np, values, blocks)
                for metric, values in component_series[component.key].items()
            },
            "contact_cutoff_sensitivity": {
                f"{cutoff:g}": {
                    metric: _statistics(np, values, blocks)
                    for metric, values in metrics.items()
                }
                for cutoff, metrics in sensitivity_series[
                    component.key
                ].items()
            },
            "terminal_state_populations_percent": {
                str(count): _statistics(np, values, blocks)
                for count, values in terminal_series[component.key].items()
            },
        }

    manifest_path = build_directory / "assembly_manifest.json"
    manifest = (
        json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest_path.is_file()
        else {}
    )
    summary = {
        "method": "anchor_resolved_interfacial_structure",
        "model_assessment": assessment,
        "trajectory": str(trajectory_path),
        "topology": str(build_directory / "topology_output.lmp"),
        "frames_in_trajectory": frame_count,
        "frames_analyzed": len(rows),
        "analyzed_frame_indices": analyzed_indices,
        "requested_blocks": blocks,
        "contact_cutoff_angstrom": float(contact_cutoff),
        "contact_cutoff_method": cutoff_method,
        "auto_inferred_contact_cutoff_angstrom": float(inferred),
        "contact_sensitivity_cutoffs_angstrom": list(sensitivity_cutoffs),
        "surface_coordination_cutoff_angstrom": surface_coordination_cutoff,
        "surface_site_definition": (
            "canonical exposed-Ni identities selected once from the "
            "authoritative pristine corrugated NiO slab and mapped by "
            "molecule-0 atom order into every assembled topology"
        ),
        "surface_site_source": topology.surface_site_source,
        "surface_site_count": len(topology.surface_site_ids),
        "random_seeds": {
            "packmol_seed": manifest.get("packmol_seed"),
            "velocity_seed": manifest.get("velocity_seed"),
        },
        "persistence_threshold": persistence_threshold,
        "timestep_fs": timestep_fs,
        "frame_spacing_ps": frame_spacing_ps,
        "site_ownership_metrics": {
            name: _statistics(np, values, blocks)
            for name, values in site_series.items()
        },
        "z_dipole": {
            "available": bool(dipole_moment_series),
            "workbook_reportable": assessment[
                "electronic_proxy_reportable"
            ],
            "reference": (
                "Sum(q_i * z_i) over every atom in the frame (ligands and "
                "substrate); translation-invariant because the assembled "
                "system is charge-neutral. Requires the trajectory dump to "
                "include a q (charge) column."
            ),
            "moment_e_angstrom": (
                _statistics(np, dipole_moment_series, blocks)
                if dipole_moment_series
                else None
            ),
            "areal_density_e_per_angstrom": (
                _statistics(np, dipole_density_series, blocks)
                if dipole_density_series
                else None
            ),
            "approximate_potential_step_volts": (
                _statistics(np, potential_step_series, blocks)
                if potential_step_series
                else None
            ),
            "potential_step_method": (
                "APPROXIMATE PROXY ONLY: assumes an idealized uniform "
                "parallel dipole sheet (Delta V = areal dipole density / "
                "epsilon0), ignoring lateral inhomogeneity and the real "
                "corrugated charge distribution the PPPM/slab-corrected "
                "force field actually uses. Treat as a qualitative, "
                "relative indicator across systems -- not an absolute "
                "work-function or UPS-comparable prediction."
            ),
        },
        "site_competition": {
            "workbook_reportable": assessment["kinetics_reportable"],
            "method": (
                "Counts transitions where an exposed-Ni site's sole "
                "occupying component changes to a DIFFERENT component "
                "between analyzed frames, confirmed by "
                "exchange_min_dwell_frames consecutive sole-owner frames "
                "and forgetting the prior occupant after more than "
                "exchange_max_vacancy_gap_frames consecutive empty/shared "
                "frames. Most informative when run over the deposition "
                "trajectory (molecules still competing for sites), not an "
                "already-settled compressed hold, where genuine exchange "
                "should be rare."
            ),
            "exchange_min_dwell_frames": exchange_min_dwell_frames,
            "exchange_max_vacancy_gap_frames": exchange_max_vacancy_gap_frames,
            "cross_component_exchange_event_count": exchange_event_count,
            "analyzed_step_span": analyzed_step_span,
            "analyzed_time_ps": analyzed_time_ps,
            "surface_site_count": len(topology.surface_site_ids),
            "ever_contacted_site_count": ever_contacted_site_count,
            "contested_site_count": contested_site_count,
            "exchange_rate_per_site_per_ns": (
                exchange_rate_per_site_per_ns
            ),
            "exchange_rate_per_ever_contacted_site_per_ns": (
                exchange_rate_per_ever_contacted_site_per_ns
            ),
            "exchange_rate_per_contested_site_per_ns": (
                exchange_rate_per_contested_site_per_ns
            ),
            "rate_normalization": (
                "cross-component exchange events divided by one of three "
                "denominators (the fixed canonical exposed-Ni site count, "
                "the count ever confirmed-occupied during this analyzed "
                "window, or the count confirmed-occupied by >=2 distinct "
                "components) and analyzed trajectory duration in ns; "
                "available only when --timestep-fs is set and at least two "
                "distinct timesteps are analyzed. The per-site rate is the "
                "most conservative; the per-contested-site rate answers "
                "'how fast do sites that ARE being competed for turn over'."
            ),
            "exchange_events_by_component_pair": dict(
                sorted(exchange_tracker.exchange_pair_counts.items())
            ),
        },
        "components": component_summaries,
        "z_profile": {
            "reference": "height above nearest exposed Ni site in periodic x/y",
            "minimum_angstrom": z_min,
            "maximum_angstrom": z_max,
            "bin_width_angstrom": z_bin_width,
        },
        "lateral_rdf": {
            "center_definition": (
                "centroid of graph-distant C/N core atoms for each molecule"
            ),
            "bin_width_angstrom": rdf_bin_width,
            "maximum_angstrom": rdf_rmax,
        },
        "outputs": [
            "interface_timeseries.csv",
            "contact_distance_histogram.csv",
            "z_density_profiles.csv",
            "lateral_rdf.csv",
        ],
        "provenance": provenance_fields(
            trajectory_path=trajectory_path,
            topology_path=build_directory / "topology_output.lmp",
            manifest_path=manifest_path,
            first_step=analyzed_steps[0] if analyzed_steps else None,
            last_step=analyzed_steps[-1] if analyzed_steps else None,
            first_frame_index=analyzed_indices[0] if analyzed_indices else None,
            last_frame_index=analyzed_indices[-1] if analyzed_indices else None,
            stride=stride,
            blocks=blocks,
        ),
    }
    summary_path = output / "interface_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary_path
