"""Tier-0 electronic decision gate: facet patches and molecular-dipole proxies.

This module answers one internal question before any DFT is proposed: do two
assembly routes differ, beyond their between-run spread, in quantities that
would control an interface work-function change?  It never produces a
work-function prediction.

For every analyzed frame it

1. classifies the canonical exposed Ni sites of the corrugated slab by their
   local surface normal (missing-O-neighbour direction), separating the 45
   degree {110} facets from (001) terraces;
2. enumerates DFT-sized facet windows (default 5 rows x 4 Ni = the 200-atom
   NiO(110) slab footprint) as fixed sets of site identities;
3. assigns each molecule to the window lying under its phosphonate P centroid
   along the facet normal, and records anchoring, exposed-site, tilt and
   lifted-material descriptors;
4. evaluates a Helmholtz dipole-layer proxy
   ``dphi = -(e / eps0) * sum(mu . n) / (eps_eff * A)`` from per-molecule
   dipoles, using either fixed force-field charges or a reference dipole
   (e.g. GFN2-xTB/DFT) rigidly rotated onto each conformer.

The proxy contains no bond dipole, charge transfer or substrate polarization.
Outputs are an internal gate under the ``classical-ff`` evidence scope and are
not intended for manuscript tables.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
import csv
import json
import math
import random

from ..config import ROOT
from ..lammps import parse
from .coverage import (
    count_dump_frames,
    load_type_elements,
    provenance_fields,
    _resolve_trajectory,
)
from .interfacial import (
    InterfacialTopology,
    _DIPOLE_DENSITY_TO_VOLTS,
    _adjacency,
    _dependencies,
    _index_by_atom_id,
    _minimum_image,
    _periodic_tree,
    _positions_for_ids,
    _selected_frames,
    _statistics,
    _wrap_xy,
    load_interfacial_topology,
)
from .model_scope import model_assessment


SUMMARY_NAME = "electronic_proxy_summary.json"
PATCH_CSV = "patch_candidates.csv"
FRAME_CSV = "frame_proxies.csv"
DEBYE_PER_E_ANGSTROM = 4.80320471
REPORTING_SCOPE = (
    "Internal Tier-0 decision gate. Molecular-dipole layer proxy only: no "
    "bond dipole, charge transfer, substrate polarization or depolarization "
    "beyond a scalar eps_eff. Not a work-function prediction and not for "
    "manuscript tables under the classical-ff evidence scope."
)

# Ideal facet normals that the missing-neighbour direction is snapped to.
_SQRT_HALF = math.sqrt(0.5)
_FACET_NORMALS = {
    "001": (0.0, 0.0, 1.0),
    "110+x": (_SQRT_HALF, 0.0, _SQRT_HALF),
    "110-x": (-_SQRT_HALF, 0.0, _SQRT_HALF),
    "110+y": (0.0, _SQRT_HALF, _SQRT_HALF),
    "110-y": (0.0, -_SQRT_HALF, _SQRT_HALF),
}


@dataclass(frozen=True)
class FacetFamily:
    key: str
    normal: tuple[float, float, float]
    t1_axis: int
    t2: tuple[float, float, float]
    site_ids: tuple[int, ...]
    row_spacing: float
    column_spacing: float


@dataclass(frozen=True)
class PatchWindow:
    window_id: str
    facet: str
    row_start: int
    column_start: int
    site_ids: tuple[int, ...]
    reference_area: float


@dataclass(frozen=True)
class ReferenceDipole:
    species: str
    dipole: tuple[float, float, float]
    coordinates: Any
    source: str


# --------------------------------------------------------------------------
# Surface facets and windows
# --------------------------------------------------------------------------


def classify_exposed_sites(
    topology: InterfacialTopology, *, coordination_cutoff: float | None = None
) -> dict[int, str]:
    """Label each exposed Ni site by its snapped missing-O-neighbour normal."""
    np, _, _ = _dependencies()
    cutoff = coordination_cutoff or topology.surface_coordination_cutoff
    oxygen_ids = [
        atom_id
        for atom_id, element in topology.atom_elements.items()
        if element == "O" and topology.atom_molecules[atom_id] == 0
    ]
    oxygen = np.asarray([topology.atom_coordinates[i] for i in oxygen_ids])
    (xlo, xhi), (ylo, yhi) = topology.lateral_bounds
    tree, _ = _periodic_tree(oxygen, xhi - xlo, yhi - ylo, dimensions=3)
    keys = list(_FACET_NORMALS)
    normals = np.asarray([_FACET_NORMALS[key] for key in keys])
    labels = {}
    for site_id in topology.surface_site_ids:
        position = np.asarray(topology.atom_coordinates[site_id])
        neighbours = tree.query_ball_point(position, cutoff)
        vector = -(tree.data[neighbours] - position).sum(axis=0)
        norm = float(np.linalg.norm(vector))
        if norm < 1e-6:
            labels[site_id] = "other"
            continue
        dots = normals @ (vector / norm)
        best = int(np.argmax(dots))
        labels[site_id] = keys[best] if dots[best] >= 0.95 else "other"
    return labels


def _circular_unwrap(np, values, lo: float, length: float):
    angles = 2.0 * math.pi * (values - lo) / length
    centre = math.atan2(float(np.sin(angles).mean()), float(np.cos(angles).mean()))
    centre_value = lo + (centre / (2.0 * math.pi)) * length
    return centre_value + _minimum_image(values - centre_value, length)


def build_facet_families(
    topology: InterfacialTopology,
    labels: dict[int, str],
    *,
    prefixes: tuple[str, ...] = ("110",),
) -> list[FacetFamily]:
    """Return facet families (default {110}) with their in-plane frame and lattice."""
    np, _, _ = _dependencies()
    families = []
    for key in sorted({label for label in labels.values() if label.startswith(prefixes)}):
        normal = np.asarray(_FACET_NORMALS[key])
        t1_axis = 0 if abs(normal[0]) < 0.1 else 1
        t1 = np.zeros(3)
        t1[t1_axis] = 1.0
        t2 = np.cross(normal, t1)
        t2 /= np.linalg.norm(t2)
        site_ids = tuple(sorted(s for s, label in labels.items() if label == key))
        u, v = _canonical_inplane(np, topology, site_ids, t1_axis, t2)
        rows = _cluster(np, v)
        row_values = np.asarray([v[rows == r].mean() for r in range(rows.max() + 1)])
        row_spacing = float(np.median(np.diff(row_values))) if len(row_values) > 1 else math.nan
        # In-row spacing from the sites themselves: a padded or
        # non-commensurate box must not distort the lattice estimate.
        in_row = np.sort(u[rows == 0])
        column_spacing = float(np.median(np.diff(in_row))) if len(in_row) > 1 else math.nan
        families.append(
            FacetFamily(
                key=key,
                normal=tuple(map(float, normal)),
                t1_axis=t1_axis,
                t2=tuple(map(float, t2)),
                site_ids=site_ids,
                row_spacing=row_spacing,
                column_spacing=column_spacing,
            )
        )
    return families


def _canonical_inplane(np, topology, site_ids, t1_axis, t2):
    positions = np.asarray([topology.atom_coordinates[s] for s in site_ids], dtype=float)
    other_axis = 1 - t1_axis
    lo, hi = topology.lateral_bounds[other_axis]
    positions[:, other_axis] = _circular_unwrap(np, positions[:, other_axis], lo, hi - lo)
    lo, hi = topology.lateral_bounds[t1_axis]
    u = lo + np.mod(positions[:, t1_axis] - lo, hi - lo)
    v = positions @ np.asarray(t2)
    return u, v


def _cluster(np, values, tolerance: float = 0.5):
    order = np.argsort(values)
    labels = np.empty(len(values), dtype=np.int64)
    current = 0
    previous = None
    for index in order:
        if previous is not None and values[index] - previous > tolerance:
            current += 1
        labels[index] = current
        previous = values[index]
    return labels


def build_patch_windows(
    topology: InterfacialTopology,
    families: list[FacetFamily],
    *,
    rows: int = 5,
    columns: int = 4,
    column_stride: int = 2,
) -> list[PatchWindow]:
    """Enumerate complete rows x columns site windows on each facet family."""
    np, _, _ = _dependencies()
    if rows <= 0 or columns <= 0 or column_stride <= 0:
        raise ValueError("window rows, columns and column stride must be positive")
    windows = []
    for family in families:
        u, v = _canonical_inplane(
            np, topology, family.site_ids, family.t1_axis, np.asarray(family.t2)
        )
        row_labels = _cluster(np, v)
        row_count = int(row_labels.max()) + 1
        row_values = np.asarray([v[row_labels == r].mean() for r in range(row_count)])
        lo, hi = topology.lateral_bounds[family.t1_axis]
        length = hi - lo
        span = columns * family.column_spacing
        site_array = np.asarray(family.site_ids)
        for row_start in range(0, row_count - rows + 1):
            gaps = np.diff(row_values[row_start : row_start + rows])
            if len(gaps) and gaps.max() > 1.5 * family.row_spacing:
                continue
            first_row = np.flatnonzero(row_labels == row_start)
            first_row = first_row[np.argsort(u[first_row])]
            for column_start in range(0, len(first_row), column_stride):
                u0 = u[first_row[column_start]] - 0.25 * family.column_spacing
                members = []
                complete = True
                for row in range(row_start, row_start + rows):
                    in_row = np.flatnonzero(row_labels == row)
                    offset = np.mod(u[in_row] - u0, length)
                    selected = in_row[offset < span - 1e-6]
                    if len(selected) != columns:
                        complete = False
                        break
                    members.extend(int(s) for s in site_array[selected])
                if not complete:
                    continue
                windows.append(
                    PatchWindow(
                        window_id=f"{family.key}:r{row_start}:c{column_start}",
                        facet=family.key,
                        row_start=row_start,
                        column_start=column_start,
                        site_ids=tuple(sorted(members)),
                        reference_area=rows
                        * columns
                        * family.row_spacing
                        * family.column_spacing,
                    )
                )
    return windows


# --------------------------------------------------------------------------
# Reference dipoles
# --------------------------------------------------------------------------


def _read_xyz_coordinates(np, path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0].split()[0])
    return np.asarray(
        [[float(v) for v in line.split()[1:4]] for line in lines[2 : 2 + count]]
    )


def _ligpargen_coordinates(np, species: str):
    path = ROOT / "inputs" / "molecules" / species / "ligpargen.lmp"
    data = parse(path)
    rows = sorted(data.sections["Atoms"], key=lambda row: int(row.fields[0]))
    return (
        np.asarray([[float(v) for v in row.fields[4:7]] for row in rows]),
        str(path.relative_to(ROOT)),
    )


def load_reference_dipoles(path: Path | None) -> dict[str, ReferenceDipole]:
    """Load ``{"species": {slug: {"dipole_e_angstrom"|"dipole_debye", "geometry"?}}}``.

    The dipole must be expressed in the frame of the reference geometry, whose
    atom order must equal the molecule's topology order.  The geometry
    defaults to the species' LigParGen coordinates.
    """
    if path is None:
        return {}
    np, _, _ = _dependencies()
    path = Path(path)
    document = json.loads(path.read_text(encoding="utf-8"))
    species_rows = document.get("species")
    if not isinstance(species_rows, dict) or not species_rows:
        raise ValueError(f"{path}: expected a non-empty 'species' mapping")
    references = {}
    for species, row in species_rows.items():
        if "dipole_e_angstrom" in row:
            dipole = [float(v) for v in row["dipole_e_angstrom"]]
        elif "dipole_debye" in row:
            dipole = [float(v) / DEBYE_PER_E_ANGSTROM for v in row["dipole_debye"]]
        else:
            raise ValueError(f"{path}: {species} needs dipole_e_angstrom or dipole_debye")
        if len(dipole) != 3:
            raise ValueError(f"{path}: {species} dipole must have three components")
        geometry = row.get("geometry")
        if geometry:
            geometry_path = Path(geometry)
            if not geometry_path.is_absolute():
                geometry_path = path.parent / geometry_path
            coordinates = _read_xyz_coordinates(np, geometry_path)
            geometry_source = str(geometry_path)
        else:
            coordinates, geometry_source = _ligpargen_coordinates(np, species)
        references[species] = ReferenceDipole(
            species=species,
            dipole=tuple(dipole),
            coordinates=coordinates - coordinates.mean(axis=0),
            source=f"{row.get('method', 'unspecified method')}; geometry {geometry_source}",
        )
    return references


def _kabsch(np, reference, target):
    """Rotation R minimizing |R reference - target| for centred point sets."""
    covariance = reference.T @ target
    u, _, vt = np.linalg.svd(covariance)
    sign = 1.0 if np.linalg.det(vt.T @ u.T) >= 0 else -1.0
    correction = np.diag([1.0, 1.0, sign])
    rotation = vt.T @ correction @ u.T
    rmsd = float(np.sqrt(((reference @ rotation.T - target) ** 2).sum(axis=1).mean()))
    return rotation, rmsd


# --------------------------------------------------------------------------
# Per-run preparation
# --------------------------------------------------------------------------


@dataclass
class _Prepared:
    topology: InterfacialTopology
    components: tuple[str, ...]
    ligand_ids: Any
    atom_molecule_index: Any
    heavy_mask: Any
    topology_charges: Any
    unwrap_levels: list[tuple[Any, Any]]
    molecule_keys: list[str]
    molecule_ids: list[int]
    phosphorus_index: list[Any]
    core_index: list[Any]
    molecule_slices: list[Any]
    terminal_molecule: Any
    terminal_oxygen_index: Any
    references: dict[int, ReferenceDipole]
    families: list[FacetFamily]
    competing_families: list[FacetFamily]
    windows: list[PatchWindow]
    site_index: dict[int, int]


def _alias(key: str, aliases: dict[str, str]) -> str:
    return aliases.get(key, key)


def _prepare(
    build_directory: Path,
    *,
    surface_coordination_cutoff: float,
    window_rows: int,
    window_columns: int,
    window_column_stride: int,
    aliases: dict[str, str],
    references: dict[str, ReferenceDipole],
    warnings: list[str],
) -> _Prepared:
    np, _, _ = _dependencies()
    topology = load_interfacial_topology(
        build_directory, surface_coordination_cutoff=surface_coordination_cutoff
    )
    data = parse(build_directory / "topology_output.lmp")
    graph = _adjacency(data)
    charges = {int(row.fields[0]): float(row.fields[3]) for row in data.sections["Atoms"]}

    ligand_ids = np.asarray(
        [a for m in topology.molecules for a in m.atom_ids], dtype=np.int64
    )
    position = {int(a): i for i, a in enumerate(ligand_ids)}
    molecule_index = np.empty(len(ligand_ids), dtype=np.int64)
    slices, phosphorus, cores = [], [], []
    terminal_molecule, terminal_oxygens = [], []
    unwrap_children: dict[int, list[tuple[int, int]]] = {}
    per_molecule_reference: dict[int, ReferenceDipole] = {}
    missing_reference: set[str] = set()
    for index, molecule in enumerate(topology.molecules):
        members = [position[a] for a in molecule.atom_ids]
        molecule_index[members] = index
        slices.append(np.asarray(members, dtype=np.int64))
        phosphorus.append(np.asarray([position[a] for a in molecule.phosphorus_ids]))
        cores.append(np.asarray([position[a] for a in molecule.core_atom_ids]))
        for terminal in molecule.terminal_oxygen_ids:
            terminal_molecule.append(index)
            terminal_oxygens.append([position[a] for a in terminal])
        # Breadth-first bond tree for periodic unwrapping, rooted at the first P.
        atom_set = set(molecule.atom_ids)
        root = molecule.phosphorus_ids[0]
        depth = {root: 0}
        queue = deque([root])
        while queue:
            atom = queue.popleft()
            for neighbour in graph.get(atom, ()):
                if neighbour in atom_set and neighbour not in depth:
                    depth[neighbour] = depth[atom] + 1
                    unwrap_children.setdefault(depth[neighbour], []).append(
                        (position[neighbour], position[atom])
                    )
                    queue.append(neighbour)
        if len(depth) != len(atom_set):
            raise ValueError(
                f"molecule {molecule.molecule_id} is not bond-connected; "
                "cannot unwrap it for dipole evaluation"
            )
        species = _alias(molecule.component.name, aliases)
        reference = references.get(species)
        if reference is None:
            missing_reference.add(species)
        elif len(reference.coordinates) != len(molecule.atom_ids):
            raise ValueError(
                f"reference geometry for {species} has {len(reference.coordinates)} "
                f"atoms; molecule {molecule.molecule_id} has {len(molecule.atom_ids)}"
            )
        else:
            per_molecule_reference[index] = reference
    if references and missing_reference:
        warnings.append(
            "no reference dipole for: " + ", ".join(sorted(missing_reference))
            + "; reference-dipole proxies cover only the listed species"
        )
    levels = [
        (
            np.asarray([c for c, _ in pairs], dtype=np.int64),
            np.asarray([p for _, p in pairs], dtype=np.int64),
        )
        for _, pairs in sorted(unwrap_children.items())
    ]
    labels = classify_exposed_sites(topology)
    families = build_facet_families(topology, labels)
    if not families:
        warnings.append("no {110} facet sites were identified; no patch windows")
    windows = build_patch_windows(
        topology,
        families,
        rows=window_rows,
        columns=window_columns,
        column_stride=window_column_stride,
    )
    keys = [_alias(m.component.name, aliases) for m in topology.molecules]
    components = tuple(dict.fromkeys(keys))
    return _Prepared(
        topology=topology,
        components=components,
        ligand_ids=ligand_ids,
        atom_molecule_index=molecule_index,
        heavy_mask=np.asarray(
            [topology.atom_elements[int(a)] != "H" for a in ligand_ids], dtype=bool
        ),
        topology_charges=np.asarray([charges[int(a)] for a in ligand_ids]),
        unwrap_levels=levels,
        molecule_keys=keys,
        molecule_ids=[m.molecule_id for m in topology.molecules],
        phosphorus_index=phosphorus,
        core_index=cores,
        molecule_slices=slices,
        terminal_molecule=np.asarray(terminal_molecule, dtype=np.int64),
        terminal_oxygen_index=np.asarray(terminal_oxygens, dtype=np.int64),
        references=per_molecule_reference,
        families=families,
        competing_families=build_facet_families(topology, labels, prefixes=("110", "001")),
        windows=windows,
        site_index={s: i for i, s in enumerate(topology.surface_site_ids)},
    )


# --------------------------------------------------------------------------
# Per-frame evaluation
# --------------------------------------------------------------------------


def _frame_state(np, prepared: _Prepared, frame, *, contact_cutoff, shield_cutoff):
    _, _, cKDTree = _dependencies()
    lookup = _index_by_atom_id(np, frame)
    lx = frame.bounds[0][1] - frame.bounds[0][0]
    ly = frame.bounds[1][1] - frame.bounds[1][0]
    lengths = np.asarray([lx, ly])

    raw = _positions_for_ids(np, frame, lookup, prepared.ligand_ids)
    positions = raw.copy()
    for children, parents in prepared.unwrap_levels:
        delta = raw[children, :2] - positions[parents, :2]
        positions[children, :2] = positions[parents, :2] + (
            delta - lengths * np.rint(delta / lengths)
        )
        positions[children, 2] = raw[children, 2]

    if frame.charges is not None:
        charges = np.asarray(frame.charges, dtype=float)[lookup[prepared.ligand_ids]]
        charge_source = "trajectory q column"
    else:
        charges = prepared.topology_charges
        charge_source = "topology_output.lmp charges"

    count = len(prepared.molecule_slices)
    p_centroid = np.zeros((count, 3))
    tilt_vector = np.zeros((count, 3))
    dipole_charges = np.zeros((count, 3))
    dipole_reference = np.full((count, 3), np.nan)
    net_charge = np.zeros(count)
    rmsd = []
    for index, members in enumerate(prepared.molecule_slices):
        coordinates = positions[members]
        centroid = coordinates.mean(axis=0)
        p_centroid[index] = positions[prepared.phosphorus_index[index]].mean(axis=0)
        tilt_vector[index] = (
            positions[prepared.core_index[index]].mean(axis=0) - p_centroid[index]
        )
        q = charges[members]
        net_charge[index] = float(q.sum())
        dipole_charges[index] = (q[:, None] * (coordinates - centroid)).sum(axis=0)
        reference = prepared.references.get(index)
        if reference is not None:
            rotation, deviation = _kabsch(np, reference.coordinates, coordinates - centroid)
            dipole_reference[index] = rotation @ np.asarray(reference.dipole)
            rmsd.append(deviation)

    sites = _wrap_xy(
        np,
        _positions_for_ids(np, frame, lookup, prepared.topology.surface_site_ids),
        frame.bounds,
    )
    site_tree, _ = _periodic_tree(sites, lx, ly, dimensions=3)
    oxygen = _wrap_xy(
        np, raw[prepared.terminal_oxygen_index.reshape(-1)], frame.bounds
    )
    distances, _ = site_tree.query(oxygen, k=1)
    terminal_anchored = (
        distances.reshape(-1, 3) <= contact_cutoff
    ).any(axis=1)
    anchored_terminals = np.bincount(
        prepared.terminal_molecule, weights=terminal_anchored, minlength=count
    ).astype(int)

    oxygen_tree, _ = _periodic_tree(oxygen, lx, ly, dimensions=3)
    site_anchored = np.asarray(
        [bool(oxygen_tree.query_ball_point(s, contact_cutoff)) for s in sites]
    )
    heavy = _wrap_xy(np, raw[prepared.heavy_mask], frame.bounds)
    heavy_tree, _ = _periodic_tree(heavy, lx, ly, dimensions=3)
    site_shielded = np.asarray(
        [bool(heavy_tree.query_ball_point(s, shield_cutoff)) for s in sites]
    )
    return {
        "lx": lx,
        "ly": ly,
        "positions": positions,
        "charges": charges,
        "charge_source": charge_source,
        "p_centroid": p_centroid,
        "tilt_vector": tilt_vector,
        "dipole_charges": dipole_charges,
        "dipole_reference": dipole_reference,
        "net_charge": net_charge,
        "reference_rmsd": rmsd,
        "anchored_terminals": anchored_terminals,
        "sites": sites,
        "site_anchored": site_anchored,
        "site_shielded": site_shielded,
    }


def _assign_to_family(np, prepared, state, family, *, height_below, height_above):
    """Index (into the exposed-site array) of the family site under each P centroid, or -1."""
    normal = np.asarray(family.normal)
    t1 = np.zeros(3)
    t1[family.t1_axis] = 1.0
    t2 = np.asarray(family.t2)
    family_index = np.asarray([prepared.site_index[s] for s in family.site_ids])
    sites = state["sites"][family_index]
    delta = state["p_centroid"][:, None, :] - sites[None, :, :]
    delta[:, :, 0] = _minimum_image(delta[:, :, 0], state["lx"])
    delta[:, :, 1] = _minimum_image(delta[:, :, 1], state["ly"])
    height = delta @ normal
    du = delta @ t1
    dv = delta @ t2
    valid = (
        (height >= -height_below)
        & (height <= height_above)
        & (np.abs(du) <= family.column_spacing)
        & (np.abs(dv) <= family.row_spacing)
    )
    inplane = np.where(valid, du * du + dv * dv, np.inf)
    best = np.argmin(inplane, axis=1)
    rows = np.arange(len(best))
    assigned = np.where(np.isfinite(inplane[rows, best]), family_index[best], -1)
    heights = np.where(assigned >= 0, height[rows, best], np.nan)
    return assigned, heights


def _home_assignments(np, prepared, state, *, height_below, height_above):
    """Assign each molecule to at most one facet: the valid one it lies closest to.

    Projection along one facet's normal can pass through a V-groove and reach
    material sitting on the opposite facet, so every family (including the
    (001) terraces, which receive no windows) competes and the smallest
    normal height wins.
    """
    results = {
        f.key: _assign_to_family(
            np, prepared, state, f, height_below=height_below, height_above=height_above
        )
        for f in prepared.competing_families
    }
    keys = list(results)
    heights = np.vstack([np.where(results[k][0] >= 0, results[k][1], np.inf) for k in keys])
    home = np.argmin(heights, axis=0)
    has_home = np.isfinite(heights.min(axis=0))
    assignments = {}
    for family in prepared.families:
        index = keys.index(family.key)
        assigned, family_heights = results[family.key]
        keep = has_home & (home == index)
        assignments[family.key] = (
            np.where(keep, assigned, -1),
            np.where(keep, family_heights, np.nan),
        )
    return assignments


def _molecule_max_height(np, prepared, state, molecule, site_position, normal):
    coordinates = state["positions"][prepared.molecule_slices[molecule]]
    anchor = state["p_centroid"][molecule]
    image = site_position.copy()
    image[0] = anchor[0] + _minimum_image(site_position[0] - anchor[0], state["lx"])
    image[1] = anchor[1] + _minimum_image(site_position[1] - anchor[1], state["ly"])
    return float(((coordinates - image) @ normal).max())


def _tilt_degrees(np, vector, axis) -> float:
    norm = float(np.linalg.norm(vector))
    if norm <= 0:
        return math.nan
    cosine = max(-1.0, min(1.0, float(vector @ axis) / norm))
    return math.degrees(math.acos(cosine))


def _dphi(dipole_normal: float, area: float, epsilon_eff: float) -> float:
    if not area > 0 or not math.isfinite(dipole_normal):
        return math.nan
    return -_DIPOLE_DENSITY_TO_VOLTS * dipole_normal / (epsilon_eff * area)


def _window_area_scale(np, prepared, family, state) -> float:
    (xlo, xhi), (ylo, yhi) = prepared.topology.lateral_bounds
    scale = np.asarray([state["lx"] / (xhi - xlo), state["ly"] / (yhi - ylo), 1.0])
    t2 = np.asarray(family.t2) * scale
    return float(scale[family.t1_axis] * np.linalg.norm(t2))


def _patch_rows(
    np,
    prepared,
    state,
    frame,
    *,
    epsilon_eff,
    height_below,
    height_above,
    lift_height,
):
    components = prepared.components
    has_reference = bool(prepared.references)
    family_lookup = {f.key: f for f in prepared.families}
    assignments = _home_assignments(
        np, prepared, state, height_below=height_below, height_above=height_above
    )
    keys = np.asarray(prepared.molecule_keys)
    rows = []
    for window in prepared.windows:
        family = family_lookup[window.facet]
        normal = np.asarray(family.normal)
        assigned, heights = assignments[window.facet]
        site_indices = np.asarray([prepared.site_index[s] for s in window.site_ids])
        members = np.flatnonzero(np.isin(assigned, site_indices))
        area = window.reference_area * _window_area_scale(np, prepared, family, state)
        row: dict[str, Any] = {
            "frame_index": frame.index,
            "step": frame.step,
            "window_id": window.window_id,
            "facet": window.facet,
            "row_start": window.row_start,
            "column_start": window.column_start,
            "n_sites": len(window.site_ids),
            "area_A2": area,
            "molecules_total": len(members),
        }
        anchored_total = 0
        lifted = 0
        tilts_all = []
        for component in components:
            selected = members[keys[members] == component]
            anchored = state["anchored_terminals"][selected]
            anchored_total += int(anchored.sum())
            tilts = [
                _tilt_degrees(np, state["tilt_vector"][m], normal) for m in selected
            ]
            tilts_all.extend(tilts)
            row[f"molecules_{component}"] = len(selected)
            row[f"anchored_terminals_{component}"] = int(anchored.sum())
            row[f"anchored_terminal_density_{component}_per_nm2"] = (
                100.0 * anchored.sum() / area if area > 0 else math.nan
            )
            row[f"multi_anchored_molecules_{component}"] = int((anchored >= 2).sum())
            row[f"unanchored_molecules_{component}"] = int((anchored == 0).sum())
            row[f"mean_tilt_deg_{component}"] = (
                float(np.nanmean(tilts)) if len(tilts) else math.nan
            )
        for molecule in members:
            if (
                state["anchored_terminals"][molecule] == 0
                and heights[molecule] > lift_height
            ):
                lifted += 1
        row["anchored_terminal_density_total_per_nm2"] = (
            100.0 * anchored_total / area if area > 0 else math.nan
        )
        row["lifted_molecules"] = lifted
        row["anchored_site_fraction"] = float(state["site_anchored"][site_indices].mean())
        row["exposed_site_fraction"] = float(
            1.0 - state["site_shielded"][site_indices].mean()
        )
        row["mean_tilt_deg"] = float(np.nanmean(tilts_all)) if tilts_all else math.nan
        row["mean_p_height_A"] = (
            float(np.nanmean(heights[members])) if len(members) else math.nan
        )
        row["max_molecule_height_A"] = (
            max(
                _molecule_max_height(
                    np,
                    prepared,
                    state,
                    int(m),
                    state["sites"][int(assigned[m])],
                    normal,
                )
                for m in members
            )
            if len(members)
            else math.nan
        )
        mu_charges = float((state["dipole_charges"][members] @ normal).sum())
        row["dipole_normal_charges_e_A"] = mu_charges
        row["dphi_layer_charges_eV"] = _dphi(mu_charges, area, epsilon_eff)
        if has_reference:
            reference = state["dipole_reference"][members]
            covered = np.isfinite(reference[:, 0])
            mu_reference = float((reference[covered] @ normal).sum())
            row["reference_coverage_fraction"] = (
                float(covered.mean()) if len(members) else math.nan
            )
            row["dipole_normal_reference_e_A"] = mu_reference
            row["dphi_layer_reference_eV"] = _dphi(mu_reference, area, epsilon_eff)
        row["molecule_ids"] = ";".join(str(prepared.molecule_ids[m]) for m in members)
        rows.append(row)
    return rows


def _frame_row(np, prepared, state, frame, *, epsilon_eff):
    area = state["lx"] * state["ly"]
    keys = np.asarray(prepared.molecule_keys)
    row: dict[str, Any] = {
        "frame_index": frame.index,
        "step": frame.step,
        "lx_A": state["lx"],
        "ly_A": state["ly"],
        "projected_area_A2": area,
        "molecules_total": len(keys),
    }
    anchored_total = 0
    for component in prepared.components:
        selected = keys == component
        anchored = state["anchored_terminals"][selected]
        anchored_total += int(anchored.sum())
        row[f"anchored_terminals_{component}"] = int(anchored.sum())
        row[f"anchored_terminal_density_{component}_per_nm2"] = 100.0 * anchored.sum() / area
        row[f"dipole_z_charges_{component}_e_A"] = float(
            state["dipole_charges"][selected, 2].sum()
        )
    row["anchored_terminal_density_total_per_nm2"] = 100.0 * anchored_total / area
    mu_charges = float(state["dipole_charges"][:, 2].sum())
    row["dipole_z_charges_e_A"] = mu_charges
    row["dphi_layer_charges_eV"] = _dphi(mu_charges, area, epsilon_eff)
    if prepared.references:
        reference = state["dipole_reference"]
        covered = np.isfinite(reference[:, 0])
        mu_reference = float(reference[covered, 2].sum())
        row["reference_coverage_fraction"] = float(covered.mean())
        row["dipole_z_reference_e_A"] = mu_reference
        row["dphi_layer_reference_eV"] = _dphi(mu_reference, area, epsilon_eff)
        row["reference_rmsd_mean_A"] = (
            float(np.mean(state["reference_rmsd"])) if state["reference_rmsd"] else math.nan
        )
    row["max_abs_molecule_net_charge_e"] = float(np.abs(state["net_charge"]).max())
    return row


# --------------------------------------------------------------------------
# Public analysis entry point
# --------------------------------------------------------------------------


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _percentiles(np, rows, column):
    values = np.asarray(
        [float(r[column]) for r in rows if r.get(column) not in (None, "")], dtype=float
    )
    values = values[np.isfinite(values)]
    if not len(values):
        return {"count": 0, "p10": math.nan, "p50": math.nan, "p90": math.nan}
    p10, p50, p90 = np.percentile(values, [10, 50, 90])
    return {"count": int(len(values)), "p10": float(p10), "p50": float(p50), "p90": float(p90)}


def parse_aliases(values: Iterable[str]) -> dict[str, str]:
    aliases = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"component alias {value!r} must look like OLD=NEW")
        old, new = (part.strip() for part in value.split("=", 1))
        if not old or not new:
            raise ValueError(f"component alias {value!r} must look like OLD=NEW")
        aliases[old] = new
    return aliases


def analyze_electronic_proxy(
    build_directory: Path,
    trajectory: Path | None = None,
    output: Path | None = None,
    *,
    last_frames: int | None = 100,
    stride: int = 1,
    blocks: int = 5,
    contact_cutoff: float = 3.25,
    shield_cutoff: float = 4.0,
    surface_coordination_cutoff: float = 2.8,
    window_rows: int = 5,
    window_columns: int = 4,
    window_column_stride: int = 2,
    height_below: float = 2.0,
    height_above: float = 30.0,
    lift_height: float = 6.0,
    epsilon_eff: float = 1.0,
    reference_dipoles: Path | None = None,
    component_aliases: dict[str, str] | None = None,
    export_conformers: int = 0,
    seed: int = 20261102,
) -> Path:
    """Run the Tier-0 proxy analysis and return ``electronic_proxy_summary.json``."""
    np, _, _ = _dependencies()
    if stride <= 0 or blocks <= 0:
        raise ValueError("stride and blocks must be positive")
    if last_frames is not None and last_frames <= 0:
        raise ValueError("last_frames must be positive or None")
    for name, value in (
        ("contact_cutoff", contact_cutoff),
        ("shield_cutoff", shield_cutoff),
        ("epsilon_eff", epsilon_eff),
        ("height_above", height_above),
    ):
        if not value > 0:
            raise ValueError(f"{name} must be positive")
    if height_below < 0 or export_conformers < 0:
        raise ValueError("height_below and export_conformers must be non-negative")
    build_directory = Path(build_directory)
    aliases = dict(component_aliases or {})
    warnings: list[str] = []
    references = load_reference_dipoles(reference_dipoles)
    trajectory_path = _resolve_trajectory(build_directory, trajectory)
    prepared = _prepare(
        build_directory,
        surface_coordination_cutoff=surface_coordination_cutoff,
        window_rows=window_rows,
        window_columns=window_columns,
        window_column_stride=window_column_stride,
        aliases=aliases,
        references=references,
        warnings=warnings,
    )
    output = Path(output) if output is not None else build_directory / "electronic-proxy"
    output.mkdir(parents=True, exist_ok=True)

    frame_count = count_dump_frames(trajectory_path)
    if frame_count == 0:
        raise ValueError(f"{trajectory_path}: trajectory contains no frames")
    first_index = max(0, frame_count - last_frames) if last_frames is not None else 0
    selected_indices = list(range(first_index, frame_count, stride))

    conformer_plan: dict[int, list[int]] = {}
    if export_conformers:
        generator = random.Random(seed)
        for component in prepared.components:
            molecules = [
                i for i, key in enumerate(prepared.molecule_keys) if key == component
            ]
            pairs = [(f, m) for f in selected_indices for m in molecules]
            for frame_index, molecule in generator.sample(
                pairs, min(export_conformers, len(pairs))
            ):
                conformer_plan.setdefault(frame_index, []).append(molecule)
    conformer_blocks: dict[str, list[str]] = {}

    patch_rows: list[dict[str, Any]] = []
    frame_rows: list[dict[str, Any]] = []
    first_step = last_step = None
    charge_sources: set[str] = set()
    for frame in _selected_frames(trajectory_path, first_index, stride):
        state = _frame_state(
            np, prepared, frame, contact_cutoff=contact_cutoff, shield_cutoff=shield_cutoff
        )
        charge_sources.add(state["charge_source"])
        frame_rows.append(_frame_row(np, prepared, state, frame, epsilon_eff=epsilon_eff))
        patch_rows.extend(
            _patch_rows(
                np,
                prepared,
                state,
                frame,
                epsilon_eff=epsilon_eff,
                height_below=height_below,
                height_above=height_above,
                lift_height=lift_height,
            )
        )
        for molecule in conformer_plan.get(frame.index, []):
            key = prepared.molecule_keys[molecule]
            conformer_blocks.setdefault(key, []).append(
                _conformer_block(np, prepared, state, frame, molecule)
            )
        first_step = frame.step if first_step is None else first_step
        last_step = frame.step
    if not frame_rows:
        raise ValueError("no frames were analyzed")

    _write_csv(output / FRAME_CSV, frame_rows)
    _write_csv(output / PATCH_CSV, patch_rows)
    conformer_files = {}
    if conformer_blocks:
        directory = output / "conformers"
        directory.mkdir(exist_ok=True)
        for key, blocks_text in conformer_blocks.items():
            path = directory / f"{key}.extxyz"
            path.write_text("".join(blocks_text), encoding="utf-8")
            conformer_files[key] = str(path.relative_to(output))

    max_charge = max(r["max_abs_molecule_net_charge_e"] for r in frame_rows)
    if max_charge > 1e-3:
        warnings.append(
            f"largest molecular net charge is {max_charge:.4f} e; molecular "
            "dipoles are evaluated about each molecule's geometric centroid"
        )
    whole_metrics = ["anchored_terminal_density_total_per_nm2", "dipole_z_charges_e_A", "dphi_layer_charges_eV"]
    whole_metrics += [f"anchored_terminal_density_{c}_per_nm2" for c in prepared.components]
    if references:
        whole_metrics += ["dipole_z_reference_e_A", "dphi_layer_reference_eV", "reference_rmsd_mean_A"]
    patch_metrics = [
        "anchored_terminal_density_total_per_nm2",
        "anchored_site_fraction",
        "exposed_site_fraction",
        "mean_tilt_deg",
        "lifted_molecules",
        "max_molecule_height_A",
        "dphi_layer_charges_eV",
    ] + [f"anchored_terminal_density_{c}_per_nm2" for c in prepared.components]
    if references:
        patch_metrics.append("dphi_layer_reference_eV")

    summary = {
        "analysis": "tier0-electronic-proxy",
        "reporting_scope": REPORTING_SCOPE,
        "model_assessment": model_assessment("classical-ff"),
        "build_directory": str(build_directory),
        "trajectory": str(trajectory_path),
        "provenance": provenance_fields(
            trajectory_path=trajectory_path,
            topology_path=build_directory / "topology_output.lmp",
            manifest_path=build_directory / "assembly_manifest.json",
            first_step=first_step,
            last_step=last_step,
            first_frame_index=frame_rows[0]["frame_index"],
            last_frame_index=frame_rows[-1]["frame_index"],
            stride=stride,
            blocks=blocks,
        ),
        "parameters": {
            "last_frames": last_frames,
            "stride": stride,
            "blocks": blocks,
            "contact_cutoff_A": contact_cutoff,
            "shield_cutoff_A": shield_cutoff,
            "surface_coordination_cutoff_A": surface_coordination_cutoff,
            "window_rows": window_rows,
            "window_columns": window_columns,
            "window_column_stride": window_column_stride,
            "height_below_A": height_below,
            "height_above_A": height_above,
            "lift_height_A": lift_height,
            "epsilon_eff": epsilon_eff,
            "component_aliases": aliases,
            "dphi_convention": (
                "dphi = -(e/eps0) * sum(mu . n) / (eps_eff * area); mu points "
                "from negative to positive charge; positive dphi raises the "
                "work function"
            ),
            "charge_source": sorted(charge_sources),
            "reference_dipoles": (
                {
                    "file": str(reference_dipoles),
                    "species": {k: v.source for k, v in references.items()},
                    "dipoles_e_angstrom": {k: list(v.dipole) for k, v in references.items()},
                    "approximation": "rigid rotation by heavy- and light-atom Kabsch fit",
                }
                if references
                else None
            ),
        },
        "components": list(prepared.components),
        "surface": {
            "site_source": prepared.topology.surface_site_source,
            "exposed_sites": len(prepared.topology.surface_site_ids),
            "facets": [
                {
                    "key": f.key,
                    "sites": len(f.site_ids),
                    "normal": list(f.normal),
                    "row_spacing_A": f.row_spacing,
                    "column_spacing_A": f.column_spacing,
                }
                for f in prepared.families
            ],
            "windows": len(prepared.windows),
            "window_note": (
                "windows overlap when column_stride < window_columns; patch "
                "distributions are descriptive, not independent samples"
            ),
        },
        "frames_analyzed": len(frame_rows),
        "whole_interface": {
            metric: _statistics(np, [r.get(metric, math.nan) for r in frame_rows], blocks)
            for metric in whole_metrics
        },
        "patches": {
            "candidates": len(patch_rows),
            "percentiles": {m: _percentiles(np, patch_rows, m) for m in patch_metrics},
        },
        "outputs": {
            "frame_csv": FRAME_CSV,
            "patch_csv": PATCH_CSV,
            "conformers": conformer_files,
        },
        "warnings": warnings,
    }
    path = output / SUMMARY_NAME
    path.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    return path


def _conformer_block(np, prepared, state, frame, molecule) -> str:
    members = prepared.molecule_slices[molecule]
    coordinates = state["positions"][members]
    charges = state["charges"][members]
    elements = [
        prepared.topology.atom_elements[int(prepared.ligand_ids[i])] for i in members
    ]
    mu_q = state["dipole_charges"][molecule]
    mu_ref = state["dipole_reference"][molecule]
    info = (
        f'frame_index={frame.index} step={frame.step} '
        f'molecule_id={prepared.molecule_ids[molecule]} '
        f'component={prepared.molecule_keys[molecule]} '
        f'anchored_terminals={int(state["anchored_terminals"][molecule])} '
        f'dipole_charges_e_A="{mu_q[0]:.6f} {mu_q[1]:.6f} {mu_q[2]:.6f}" '
        f'dipole_reference_e_A="{mu_ref[0]:.6f} {mu_ref[1]:.6f} {mu_ref[2]:.6f}" '
        'pbc="F F F" Properties=species:S:1:pos:R:3:charge:R:1'
    )
    lines = [str(len(members)), info]
    lines += [
        f"{element} {x:.6f} {y:.6f} {z:.6f} {q:.6f}"
        for element, (x, y, z), q in zip(elements, coordinates, charges)
    ]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Run-level gate and slot selection
# --------------------------------------------------------------------------


def _load_summary(path: Path) -> tuple[Path, dict[str, Any]]:
    path = Path(path)
    if path.is_dir():
        path = path / SUMMARY_NAME
    return path, json.loads(path.read_text(encoding="utf-8"))


def parse_labelled(values: Iterable[str]) -> list[tuple[str, Path]]:
    entries = []
    for value in values:
        if "=" not in value:
            raise ValueError(f"run {value!r} must look like LABEL=SUMMARY_JSON_OR_DIR")
        label, path = value.split("=", 1)
        if not label.strip():
            raise ValueError(f"run {value!r} has an empty label")
        entries.append((label.strip(), Path(path.strip())))
    return entries


def _sample_sd(values: list[float]) -> float:
    if len(values) < 2:
        return math.nan
    mean = sum(values) / len(values)
    return math.sqrt(sum((v - mean) ** 2 for v in values) / (len(values) - 1))


def summarize_electronic_gate(
    runs: list[tuple[str, Path]],
    output: Path,
    *,
    compare: tuple[str, str],
    threshold_sd: float = 2.0,
) -> Path:
    """Compare two labelled groups of run summaries against the Tier-0 gate."""
    if threshold_sd <= 0:
        raise ValueError("threshold_sd must be positive")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    loaded = [(label, *_load_summary(path)) for label, path in runs]
    metrics: dict[str, str] = {}
    for _, _, summary in loaded:
        for metric in summary["whole_interface"]:
            metrics[f"whole:{metric}"] = "whole_interface"
        for metric in summary["patches"]["percentiles"]:
            metrics[f"patch_median:{metric}"] = "patches"
    run_rows = []
    for label, path, summary in loaded:
        row: dict[str, Any] = {"label": label, "summary": str(path), "build_directory": summary["build_directory"]}
        for name, section in metrics.items():
            metric = name.split(":", 1)[1]
            if section == "whole_interface":
                value = summary["whole_interface"].get(metric, {}).get("mean", math.nan)
            else:
                value = summary["patches"]["percentiles"].get(metric, {}).get("p50", math.nan)
            row[name] = value
        run_rows.append(row)

    group_a, group_b = compare
    gate_rows = []
    for name in metrics:
        a = [float(r[name]) for r in run_rows if r["label"] == group_a and math.isfinite(float(r[name]))]
        b = [float(r[name]) for r in run_rows if r["label"] == group_b and math.isfinite(float(r[name]))]
        mean_a = sum(a) / len(a) if a else math.nan
        mean_b = sum(b) / len(b) if b else math.nan
        sd_a, sd_b = _sample_sd(a), _sample_sd(b)
        pooled = math.sqrt((sd_a**2 + sd_b**2) / 2.0) if math.isfinite(sd_a) and math.isfinite(sd_b) else math.nan
        difference = mean_b - mean_a
        ratio = abs(difference) / pooled if pooled and math.isfinite(pooled) and pooled > 0 else math.nan
        evaluable = len(a) >= 2 and len(b) >= 2 and math.isfinite(ratio)
        gate_rows.append(
            {
                "metric": name,
                f"n_{group_a}": len(a),
                f"n_{group_b}": len(b),
                f"mean_{group_a}": mean_a,
                f"mean_{group_b}": mean_b,
                f"sd_{group_a}": sd_a,
                f"sd_{group_b}": sd_b,
                "difference_b_minus_a": difference,
                "pooled_between_run_sd": pooled,
                "difference_over_pooled_sd": ratio,
                "evaluable": evaluable,
                "passes": bool(evaluable and ratio >= threshold_sd),
            }
        )
    primary = [
        "whole:dphi_layer_charges_eV",
        "whole:dphi_layer_reference_eV",
        "patch_median:dphi_layer_charges_eV",
        "patch_median:dphi_layer_reference_eV",
    ] + [m for m in metrics if m.startswith("whole:anchored_terminal_density_")]
    primary_rows = [r for r in gate_rows if r["metric"] in primary]
    decision = {
        "compare": {"a": group_a, "b": group_b},
        "threshold": f"|mean_b - mean_a| >= {threshold_sd} x pooled between-run SD, n >= 2 per group",
        "primary_metrics": [r["metric"] for r in primary_rows],
        "passing_primary_metrics": [r["metric"] for r in primary_rows if r["passes"]],
        "unevaluable_primary_metrics": [r["metric"] for r in primary_rows if not r["evaluable"]],
        "proceed_to_tier1": any(r["passes"] for r in primary_rows),
        "note": (
            "A pass justifies only the Tier-1 DFT timing benchmark. With three "
            "runs per group the SD estimate is itself uncertain; treat marginal "
            "passes as inconclusive. " + REPORTING_SCOPE
        ),
    }
    _write_csv(output / "electronic_gate_runs.csv", run_rows)
    _write_csv(output / "electronic_gate.csv", gate_rows)
    path = output / "electronic_gate.json"
    path.write_text(
        json.dumps({"decision": decision, "metrics": gate_rows, "runs": run_rows}, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    return path


def _read_candidates(summary_path: Path, summary: dict[str, Any]) -> list[dict[str, Any]]:
    path = summary_path.parent / summary["outputs"]["patch_csv"]
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key, value in list(row.items()):
            if key in {"window_id", "facet", "molecule_ids"}:
                continue
            try:
                row[key] = float(value) if value not in ("", None) else math.nan
            except ValueError:
                pass
    return rows


def _finite(row, key) -> float:
    value = row.get(key, math.nan)
    return float(value) if isinstance(value, (int, float)) and math.isfinite(value) else math.nan


def _median(values):
    ordered = sorted(v for v in values if math.isfinite(v))
    if not ordered:
        return math.nan
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else 0.5 * (ordered[middle - 1] + ordered[middle])


def _spread(values):
    ordered = sorted(v for v in values if math.isfinite(v))
    if len(ordered) < 2:
        return 1.0
    low = ordered[int(0.1 * (len(ordered) - 1))]
    high = ordered[int(0.9 * (len(ordered) - 1))]
    return high - low or 1.0


def select_patch_slots(
    runs: list[tuple[str, Path]],
    output: Path,
    *,
    compare: tuple[str, str],
    matched_pairs: int = 2,
    matched_tolerance: float = 0.10,
    reference_slots: int = 2,
) -> Path:
    """Choose Tier-1/2 patch slots: per-run median and extreme, plus matched pairs."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    density = "anchored_terminal_density_total_per_nm2"
    arrangement = ("mean_tilt_deg", "exposed_site_fraction")
    slots: list[dict[str, Any]] = []
    by_label: dict[str, list[dict[str, Any]]] = {}

    def add(label, summary_path, summary, row, criterion):
        slots.append(
            {
                "slot": len(slots) + 1,
                "label": label,
                "criterion": criterion,
                "build_directory": summary["build_directory"],
                "trajectory": summary["trajectory"],
                "summary": str(summary_path),
                "frame_index": int(row["frame_index"]),
                "step": int(row["step"]),
                "window_id": row["window_id"],
                "facet": row["facet"],
                density: _finite(row, density),
                "mean_tilt_deg": _finite(row, "mean_tilt_deg"),
                "exposed_site_fraction": _finite(row, "exposed_site_fraction"),
                "lifted_molecules": _finite(row, "lifted_molecules"),
                "dphi_layer_charges_eV": _finite(row, "dphi_layer_charges_eV"),
                "molecule_ids": row.get("molecule_ids", ""),
            }
        )

    for label, path in runs:
        summary_path, summary = _load_summary(path)
        candidates = [r for r in _read_candidates(summary_path, summary) if _finite(r, "molecules_total") > 0]
        if not candidates:
            continue
        for row in candidates:
            row["_run"] = (label, summary_path, summary)
        by_label.setdefault(label, []).extend(candidates)
        target = _median([_finite(r, density) for r in candidates])
        median_row = min(candidates, key=lambda r: (abs(_finite(r, density) - target), r["frame_index"]))
        add(label, summary_path, summary, median_row, "run-median anchored density")
        if label not in compare and reference_slots < 2:
            # Reference-only runs (e.g. single gap-seeded runs) may request
            # just the median slot.
            continue
        centres = {k: _median([_finite(r, k) for r in candidates]) for k in arrangement}
        spreads = {k: _spread([_finite(r, k) for r in candidates]) for k in arrangement}

        def extremeness(row):
            distance = sum(
                ((_finite(row, k) - centres[k]) / spreads[k]) ** 2
                for k in arrangement
                if math.isfinite(_finite(row, k))
            )
            multi = any(
                key.startswith("multi_anchored_molecules_") and _finite(row, key) > 0
                for key in row
            )
            return (multi, distance)

        add(label, summary_path, summary, max(candidates, key=extremeness), "arrangement extreme (multi-anchor preferred)")

    group_a, group_b = compare
    pairs = []
    for left in by_label.get(group_a, []):
        density_left = _finite(left, density)
        if not math.isfinite(density_left) or density_left <= 0:
            continue
        for right in by_label.get(group_b, []):
            density_right = _finite(right, density)
            if not math.isfinite(density_right):
                continue
            mean = 0.5 * (density_left + density_right)
            if abs(density_left - density_right) > matched_tolerance * mean:
                continue
            contrast = abs(_finite(left, "mean_tilt_deg") - _finite(right, "mean_tilt_deg")) / 10.0 + abs(
                _finite(left, "exposed_site_fraction") - _finite(right, "exposed_site_fraction")
            ) * 10.0
            if math.isfinite(contrast):
                pairs.append((contrast, left, right))
    pairs.sort(key=lambda item: -item[0])
    used_runs: set[tuple[str, str]] = set()
    chosen = 0
    for contrast, left, right in pairs:
        if chosen >= matched_pairs:
            break
        key = (str(left["_run"][1]), str(right["_run"][1]))
        if key in used_runs:
            continue
        used_runs.add(key)
        chosen += 1
        for row in (left, right):
            label, summary_path, summary = row["_run"]
            add(label, summary_path, summary, row, f"matched-density pair {chosen} (contrast {contrast:.2f})")
    path = output / "patch_slots.csv"
    _write_csv(path, slots)
    if chosen < matched_pairs:
        (output / "patch_slots_warnings.txt").write_text(
            f"only {chosen} of {matched_pairs} matched-density pairs were found within "
            f"+/-{100 * matched_tolerance:.0f}% anchored density\n",
            encoding="utf-8",
        )
    return path
