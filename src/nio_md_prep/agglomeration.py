"""Generate reproducible phosphonate agglomeration structures for DFT labeling."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shlex
import shutil
import subprocess
import tomllib
from dataclasses import dataclass
from pathlib import Path

from . import __version__
from .chemistry import phosphonate_roles
from .config import ROOT, missing_ligpargen, molecule_manifest
from .geometry import elements
from .lammps import DataFile, atom_coordinates, charge, parse


@dataclass(frozen=True)
class MoleculeTemplate:
    slug: str
    count: int
    source: Path
    symbols: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    masses: tuple[float, ...]
    phosphonate_p_indices: tuple[int, ...]
    phosphonate_c_indices: tuple[int, ...]
    net_charge: float


@dataclass(frozen=True)
class MoleculeInstance:
    slug: str
    molecule_id: int
    start: int
    stop: int
    masses: tuple[float, ...]
    phosphonate_p_indices: tuple[int, ...] = ()
    phosphonate_c_indices: tuple[int, ...] = ()


@dataclass(frozen=True)
class AgglomerateSpec:
    name: str
    templates: tuple[MoleculeTemplate, ...]
    replicas: int
    base_seed: int
    radius: float
    tolerance: float
    vacuum: float
    minimum_distance: float
    compact_to_distance: float | None
    scales: tuple[float, ...]


_VASP_OUTPUTS = {
    "POSCAR", "CONTCAR", "OUTCAR", "XDATCAR", "vasprun.xml", "OSZICAR",
    "WAVECAR", "CHG", "CHGCAR", "EIGENVAL", "DOSCAR", "PROCAR",
}
_REQUIRED_VASP_INPUTS = ("INCAR", "KPOINTS", "POTCAR")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atom_masses(data: DataFile) -> tuple[float, ...]:
    by_type = {
        int(row.fields[0]): float(row.fields[1])
        for row in data.sections["Masses"]
    }
    return tuple(by_type[int(atom.fields[2])] for atom in data.sections["Atoms"])


def _load_templates(config: dict) -> list[MoleculeTemplate]:
    specifications = config.get("molecules", [])
    if not specifications:
        raise ValueError("agglomeration config requires at least one [[molecules]] entry")
    templates: list[MoleculeTemplate] = []
    total = 0
    seen: set[str] = set()
    for specification in specifications:
        slug = str(specification.get("slug", "")).strip()
        count = int(specification.get("count", 0))
        if not slug or count < 1:
            raise ValueError("every agglomeration molecule requires a slug and positive count")
        if slug in seen:
            raise ValueError(f"duplicate agglomeration molecule slug: {slug}")
        seen.add(slug)
        folder, manifest = molecule_manifest(slug)
        source = folder / manifest.get("files", {}).get("ligpargen", "ligpargen.lmp")
        if not source.exists():
            raise missing_ligpargen(source, "nio-md-prep prepare-agglomeration ...")
        data = parse(source)
        try:
            expected_charge = float(
                manifest.get("molecule", {}).get("expected_net_charge")
            )
        except (TypeError, ValueError):
            raise ValueError(
                f"{folder / 'molecule.toml'}: expected_net_charge must be reviewed and numeric"
            ) from None
        actual_charge = charge(data)
        if abs(actual_charge - expected_charge) > 1.0e-6:
            raise ValueError(
                f"{source}: LigParGen charge {actual_charge:.8f} does not match "
                f"manifest expected_net_charge {expected_charge}"
            )
        template_symbols = tuple(elements(data))
        atom_position = {
            int(atom.fields[0]): index
            for index, atom in enumerate(data.sections["Atoms"])
        }
        roles = phosphonate_roles(data)
        adjacency = {atom_id: set() for atom_id in atom_position}
        for bond in data.sections.get("Bonds", []):
            left, right = map(int, bond.fields[2:4])
            adjacency[left].add(right)
            adjacency[right].add(left)
        p_atom_ids = sorted(
            (atom_id for atom_id, role in roles.items() if role == "P"),
            key=atom_position.__getitem__,
        )
        p_indices = tuple(atom_position[atom_id] for atom_id in p_atom_ids)
        c_indices = tuple(
            atom_position[
                next(
                    neighbor
                    for neighbor in adjacency[atom_id]
                    if template_symbols[atom_position[neighbor]] == "C"
                )
            ]
            for atom_id in p_atom_ids
        )
        templates.append(
            MoleculeTemplate(
                slug=slug,
                count=count,
                source=source,
                symbols=template_symbols,
                coordinates=tuple(atom_coordinates(data)),
                masses=_atom_masses(data),
                phosphonate_p_indices=p_indices,
                phosphonate_c_indices=c_indices,
                net_charge=actual_charge,
            )
        )
        total += count
    if total < 2:
        raise ValueError("agglomeration generation requires at least two molecules")
    return templates


def _agglomerate_specs(config: dict) -> tuple[list[AgglomerateSpec], bool]:
    defaults = config.get("agglomeration", {})
    explicit = config.get("agglomerates", [])
    raw_specs = explicit or [
        {
            "name": "default",
            "molecules": config.get("molecules", []),
        }
    ]
    specs: list[AgglomerateSpec] = []
    names: set[str] = set()
    seeds: set[int] = set()
    for index, raw in enumerate(raw_specs):
        name = str(raw.get("name", f"set-{index:02d}")).strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
            raise ValueError(
                f"invalid agglomerate name {name!r}; use letters, numbers, '_' or '-'"
            )
        if name in names:
            raise ValueError(f"duplicate agglomerate name: {name}")
        names.add(name)
        templates = tuple(_load_templates({"molecules": raw.get("molecules", [])}))
        replicas = int(raw.get("replicas", defaults.get("replicas", 4)))
        base_seed = int(raw.get("base_seed", defaults.get("base_seed", 2405367)))
        radius = float(raw.get("radius_angstrom", defaults.get("radius_angstrom", 14.0)))
        tolerance = float(
            raw.get(
                "packmol_tolerance_angstrom",
                defaults.get("packmol_tolerance_angstrom", 3.0),
            )
        )
        vacuum = float(raw.get("vacuum_angstrom", defaults.get("vacuum_angstrom", 8.0)))
        minimum_distance = float(
            raw.get(
                "minimum_distance_angstrom",
                defaults.get("minimum_distance_angstrom", 2.5),
            )
        )
        compact_value = raw.get(
            "compact_to_distance_angstrom",
            defaults.get("compact_to_distance_angstrom"),
        )
        compact_to_distance = (
            None if compact_value is None else float(compact_value)
        )
        scales = tuple(
            float(value)
            for value in raw.get(
                "center_scales", defaults.get("center_scales", [1.0])
            )
        )
        if replicas < 1 or not 0 < base_seed < 900000000:
            raise ValueError(
                f"{name}: replicas must be positive and base_seed must be between "
                "1 and 899999999"
            )
        if radius <= 0 or tolerance <= 0 or vacuum <= 0 or minimum_distance <= 0:
            raise ValueError(
                f"{name}: radius, tolerance, vacuum, and minimum distance must be positive"
            )
        if compact_to_distance is not None and compact_to_distance < minimum_distance:
            raise ValueError(
                f"{name}: compact_to_distance_angstrom must be at least "
                "minimum_distance_angstrom"
            )
        if (
            not scales
            or any(scale <= 0 for scale in scales)
            or len(set(scales)) != len(scales)
        ):
            raise ValueError(
                f"{name}: center_scales must contain distinct positive values"
            )
        resolved_seeds = set(range(base_seed, base_seed + replicas))
        if max(resolved_seeds) >= 900000000:
            raise ValueError(f"{name}: resolved Packmol seed exceeds 899999999")
        overlap = seeds & resolved_seeds
        if overlap:
            raise ValueError(
                f"{name}: Packmol seeds overlap another agglomerate family: "
                + ", ".join(str(seed) for seed in sorted(overlap))
            )
        seeds.update(resolved_seeds)
        specs.append(
            AgglomerateSpec(
                name=name,
                templates=templates,
                replicas=replicas,
                base_seed=base_seed,
                radius=radius,
                tolerance=tolerance,
                vacuum=vacuum,
                minimum_distance=minimum_distance,
                compact_to_distance=compact_to_distance,
                scales=scales,
            )
        )
    return specs, bool(explicit)


def _write_template_xyz(template: MoleculeTemplate, path: Path) -> None:
    rows = [
        f"{symbol} {x:.10f} {y:.10f} {z:.10f}"
        for symbol, (x, y, z) in zip(template.symbols, template.coordinates)
    ]
    path.write_text(
        f"{len(rows)}\n{template.slug} reconstructed from {template.source.name}\n"
        + "\n".join(rows) + "\n",
        encoding="utf-8",
    )


def _read_xyz(path: Path) -> tuple[list[str], list[tuple[float, float, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path}: incomplete XYZ file")
    count = int(lines[0])
    if len(lines[2:]) != count:
        raise ValueError(f"{path}: XYZ atom count mismatch")
    symbols: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for line_number, line in enumerate(lines[2:], 3):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"{path}:{line_number}: malformed XYZ row")
        symbols.append(fields[0])
        coordinates.append(tuple(float(value) for value in fields[1:4]))
    return symbols, coordinates


def _expected_layout(
    templates: list[MoleculeTemplate],
) -> tuple[list[str], list[MoleculeInstance]]:
    symbols: list[str] = []
    instances: list[MoleculeInstance] = []
    molecule_id = 0
    for template in templates:
        for _ in range(template.count):
            molecule_id += 1
            start = len(symbols)
            symbols.extend(template.symbols)
            instances.append(
                MoleculeInstance(
                    slug=template.slug,
                    molecule_id=molecule_id,
                    start=start,
                    stop=len(symbols),
                    masses=template.masses,
                    phosphonate_p_indices=tuple(
                        start + index for index in template.phosphonate_p_indices
                    ),
                    phosphonate_c_indices=tuple(
                        start + index for index in template.phosphonate_c_indices
                    ),
                )
            )
    return symbols, instances


def _nearest_phosphonate_pairs(
    coordinates: list[tuple[float, float, float]],
    instances: list[MoleculeInstance],
) -> list[dict]:
    """Greedily pair nearby P heads on different molecules, once per head."""
    if any(
        len(item.phosphonate_p_indices) != len(item.phosphonate_c_indices)
        for item in instances
    ):
        raise ValueError("every phosphonate P head requires one bonded carbon axis atom")
    heads = [
        (p_index, c_index, item.molecule_id)
        for item in instances
        for p_index, c_index in zip(
            item.phosphonate_p_indices, item.phosphonate_c_indices
        )
    ]
    candidates = sorted(
        (
            _distance(coordinates[left_atom], coordinates[right_atom]),
            left_atom,
            right_atom,
            left_carbon,
            right_carbon,
            left_molecule,
            right_molecule,
        )
        for position, (left_atom, left_carbon, left_molecule) in enumerate(heads)
        for right_atom, right_carbon, right_molecule in heads[position + 1:]
        if left_molecule != right_molecule
    )
    used: set[int] = set()
    pairs: list[dict] = []
    for (
        distance,
        left_atom,
        right_atom,
        left_carbon,
        right_carbon,
        left_molecule,
        right_molecule,
    ) in candidates:
        if left_atom in used or right_atom in used:
            continue
        used.update((left_atom, right_atom))
        pairs.append(
            {
                "left_atom_index": left_atom + 1,
                "right_atom_index": right_atom + 1,
                "left_carbon_atom_index": left_carbon + 1,
                "right_carbon_atom_index": right_carbon + 1,
                "left_molecule_id": left_molecule,
                "right_molecule_id": right_molecule,
                "initial_distance_angstrom": distance,
            }
        )
    return pairs


def _packmol_input(
    templates: list[MoleculeTemplate],
    template_paths: dict[str, Path],
    packed_name: str,
    seed: int,
    tolerance: float,
    radius: float,
) -> str:
    lines = [
        f"tolerance {tolerance:.6f}",
        "filetype xyz",
        f"output {packed_name}",
        f"seed {seed}",
        "",
    ]
    for template in templates:
        lines.extend(
            (
                f"structure {template_paths[template.slug].name}",
                f"  number {template.count}",
                f"  inside sphere 0.0 0.0 0.0 {radius:.6f}",
                "end structure",
                "",
            )
        )
    return "\n".join(lines)


def _center_of_mass(
    coordinates: list[tuple[float, float, float]], masses: tuple[float, ...]
) -> tuple[float, float, float]:
    total = sum(masses)
    return tuple(
        sum(mass * point[axis] for mass, point in zip(masses, coordinates)) / total
        for axis in range(3)
    )


def _scale_molecule_centers(
    coordinates: list[tuple[float, float, float]],
    instances: list[MoleculeInstance],
    scale: float,
) -> list[tuple[float, float, float]]:
    centers = [
        _center_of_mass(coordinates[item.start:item.stop], item.masses)
        for item in instances
    ]
    molecule_masses = [sum(item.masses) for item in instances]
    total_mass = sum(molecule_masses)
    cluster_center = tuple(
        sum(mass * center[axis] for mass, center in zip(molecule_masses, centers))
        / total_mass
        for axis in range(3)
    )
    result = list(coordinates)
    for item, center in zip(instances, centers):
        translation = tuple(
            (scale - 1.0) * (center[axis] - cluster_center[axis])
            for axis in range(3)
        )
        for atom_index in range(item.start, item.stop):
            point = coordinates[atom_index]
            result[atom_index] = tuple(
                point[axis] + translation[axis] for axis in range(3)
            )
    return result


def _compact_molecule_centers(
    coordinates: list[tuple[float, float, float]],
    symbols: list[str],
    instances: list[MoleculeInstance],
    target_distance: float | None,
) -> tuple[list[tuple[float, float, float]], float]:
    """Compact and connect a Packmol cluster without crossing its contact floor."""
    if target_distance is None:
        return list(coordinates), 1.0
    current = _closest_intermolecular_pair(coordinates, symbols, instances)[
        "minimum_distance_angstrom"
    ]
    contact_floor = min(current, target_distance)
    low, high = 0.0, 1.0
    if current > target_distance:
        # Find the smallest global center scale whose closest intermolecular
        # contact remains at the target. Intramolecular geometries are untouched.
        for _ in range(60):
            middle = (low + high) / 2.0
            candidate = _scale_molecule_centers(coordinates, instances, middle)
            nearest = _closest_intermolecular_pair(candidate, symbols, instances)[
                "minimum_distance_angstrom"
            ]
            if nearest >= contact_floor:
                high = middle
            else:
                low = middle
        result = _scale_molecule_centers(coordinates, instances, high)
    else:
        result = list(coordinates)

    def instance_distance(left: int, right: int) -> tuple[float, int, int]:
        best = (math.inf, -1, -1)
        first = instances[left]
        second = instances[right]
        for left_atom in range(first.start, first.stop):
            for right_atom in range(second.start, second.stop):
                distance = _distance(result[left_atom], result[right_atom])
                if distance < best[0]:
                    best = (distance, left_atom, right_atom)
        return best

    # A global scale can leave one molecule isolated when another pair reaches
    # the contact floor first. Grow a connected cluster by pulling the nearest
    # outsider along its closest atom-pair vector. Every move is bounded by the
    # first safe contact, so this cannot recreate short O--O overlaps.
    def components() -> list[set[int]]:
        remaining = set(range(len(instances)))
        groups: list[set[int]] = []
        while remaining:
            group = {remaining.pop()}
            changed = True
            while changed:
                changed = False
                neighbors = {
                    candidate
                    for candidate in remaining
                    if any(
                        instance_distance(candidate, member)[0]
                        <= contact_floor + 1.0e-5
                        for member in group
                    )
                }
                if neighbors:
                    group.update(neighbors)
                    remaining.difference_update(neighbors)
                    changed = True
            groups.append(group)
        return groups

    groups = components()
    while len(groups) > 1:
        distance, left_atom, right_atom, moving_group = min(
            (*instance_distance(left, right), right_group)
            for left_group in groups
            for right_group in groups
            if left_group is not right_group
            for left in left_group
            for right in right_group
        )
        delta = tuple(
            (result[right_atom][axis] - result[left_atom][axis])
            * (contact_floor - distance)
            / distance
            for axis in range(3)
        )
        moving_atoms = [
            atom
            for molecule in moving_group
            for atom in range(instances[molecule].start, instances[molecule].stop)
        ]
        stationary_atoms = [
            atom
            for molecule in range(len(instances))
            if molecule not in moving_group
            for atom in range(instances[molecule].start, instances[molecule].stop)
        ]

        def moved_clearance(fraction: float) -> float:
            minimum = math.inf
            for atom in moving_atoms:
                point = tuple(
                    result[atom][axis] + fraction * delta[axis]
                    for axis in range(3)
                )
                for other_atom in stationary_atoms:
                    minimum = min(minimum, _distance(point, result[other_atom]))
            return minimum

        move_low, move_high = 0.0, 1.0
        if moved_clearance(1.0) < contact_floor:
            for _ in range(60):
                middle = (move_low + move_high) / 2.0
                if moved_clearance(middle) >= contact_floor:
                    move_low = middle
                else:
                    move_high = middle
            fraction = move_low
        else:
            fraction = 1.0
        for atom in moving_atoms:
            result[atom] = tuple(
                result[atom][axis] + fraction * delta[axis]
                for axis in range(3)
            )
        new_groups = components()
        if len(new_groups) >= len(groups):
            raise RuntimeError("adaptive compaction could not connect molecule groups")
        groups = new_groups
    return result, high


def _closest_intermolecular_pair(
    coordinates: list[tuple[float, float, float]],
    symbols: list[str],
    instances: list[MoleculeInstance],
) -> dict:
    molecule_for_atom = [0] * len(coordinates)
    for item in instances:
        molecule_for_atom[item.start:item.stop] = [item.molecule_id] * (item.stop - item.start)
    minimum = math.inf
    pair: tuple[int, int] | None = None
    minimum_oo = math.inf
    oo_pair: tuple[int, int] | None = None
    for left, first in enumerate(coordinates):
        for right in range(left + 1, len(coordinates)):
            if molecule_for_atom[left] == molecule_for_atom[right]:
                continue
            second = coordinates[right]
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))
            if distance < minimum:
                minimum = distance
                pair = (left, right)
            if symbols[left] == symbols[right] == "O" and distance < minimum_oo:
                minimum_oo = distance
                oo_pair = (left, right)

    def describe(indices: tuple[int, int] | None, distance: float) -> dict | None:
        if indices is None:
            return None
        left, right = indices
        return {
            "distance_angstrom": distance,
            "left": {
                "source_atom_index": left + 1,
                "element": symbols[left],
                "molecule_id": molecule_for_atom[left],
            },
            "right": {
                "source_atom_index": right + 1,
                "element": symbols[right],
                "molecule_id": molecule_for_atom[right],
            },
        }

    return {
        "minimum_distance_angstrom": minimum,
        "closest_pair": describe(pair, minimum),
        "minimum_oxygen_oxygen_distance_angstrom": (
            minimum_oo if oo_pair is not None else None
        ),
        "closest_oxygen_oxygen_pair": describe(oo_pair, minimum_oo),
    }


def _maximum_nearest_molecule_distance(
    coordinates: list[tuple[float, float, float]],
    instances: list[MoleculeInstance],
) -> float:
    nearest = [math.inf] * len(instances)
    for left, first in enumerate(instances):
        for right in range(left + 1, len(instances)):
            second = instances[right]
            distance = min(
                _distance(coordinates[a], coordinates[b])
                for a in range(first.start, first.stop)
                for b in range(second.start, second.stop)
            )
            nearest[left] = min(nearest[left], distance)
            nearest[right] = min(nearest[right], distance)
    return max(nearest)


def _fixed_cell(
    variants: dict[float, list[tuple[float, float, float]]], vacuum: float
) -> float:
    maximum_span = 0.0
    for coordinates in variants.values():
        for axis in range(3):
            values = [point[axis] for point in coordinates]
            maximum_span = max(maximum_span, max(values) - min(values))
    return maximum_span + 2.0 * vacuum


def _center_in_cell(
    coordinates: list[tuple[float, float, float]], cell: float
) -> list[tuple[float, float, float]]:
    bounds = [
        (min(point[axis] for point in coordinates), max(point[axis] for point in coordinates))
        for axis in range(3)
    ]
    translations = [cell / 2.0 - (low + high) / 2.0 for low, high in bounds]
    return [
        tuple(point[axis] + translations[axis] for axis in range(3))
        for point in coordinates
    ]


def _distance(
    first: tuple[float, float, float], second: tuple[float, float, float]
) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def _maximum_intramolecular_distance_deviation(
    coordinates: list[tuple[float, float, float]],
    instances: list[MoleculeInstance],
    templates: list[MoleculeTemplate],
) -> float:
    by_slug = {template.slug: template for template in templates}
    maximum = 0.0
    for instance in instances:
        actual = coordinates[instance.start:instance.stop]
        reference = by_slug[instance.slug].coordinates
        for left in range(len(actual)):
            for right in range(left + 1, len(actual)):
                maximum = max(
                    maximum,
                    abs(
                        _distance(actual[left], actual[right])
                        - _distance(reference[left], reference[right])
                    ),
                )
    return maximum


def _potcar_elements(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return re.findall(r"^\s*TITEL\s*=\s*\S+\s+([A-Z][a-z]?)", text, re.MULTILINE)


def _reference_sanity(
    reference: Path | None,
    expected_elements: list[str],
    *,
    structures_only: bool,
) -> dict:
    if reference is None:
        if not structures_only:
            raise ValueError(
                "a complete VASP --reference-dir is required; use "
                "--structures-only only when calculation inputs are intentionally omitted"
            )
        return {
            "status": "STRUCTURES_ONLY",
            "potcar_status": "NOT_SUPPLIED",
            "expected_elements": expected_elements,
            "required_inputs": list(_REQUIRED_VASP_INPUTS),
        }
    if not reference.is_dir():
        raise FileNotFoundError(f"VASP reference directory not found: {reference}")
    missing = [name for name in _REQUIRED_VASP_INPUTS if not (reference / name).is_file()]
    if missing:
        raise ValueError(
            f"{reference}: incomplete VASP reference directory; missing "
            + ", ".join(missing)
        )
    potcar = reference / "POTCAR"
    actual_elements = _potcar_elements(potcar)
    if actual_elements != expected_elements:
        raise ValueError(
            f"{potcar}: POTCAR element order {actual_elements} does not match "
            f"generated POSCAR order {expected_elements}"
        )
    return {
        "status": "PASS",
        "potcar_status": "PASS",
        "expected_elements": expected_elements,
        "potcar_elements": actual_elements,
        "potcar_sha256": _sha256(potcar),
        "required_inputs": list(_REQUIRED_VASP_INPUTS),
    }


def _structure_sanity(
    coordinates: list[tuple[float, float, float]],
    symbols: list[str],
    instances: list[MoleculeInstance],
    templates: list[MoleculeTemplate],
    *,
    cell: float,
    requested_vacuum: float,
    minimum_distance: float,
    reference_sanity: dict,
    expected_elements: list[str],
    written_element_order: list[str],
    minimum_oxygen_oxygen_distance: float | None = None,
    maximum_nearest_molecule_distance: float | None = None,
    require_rigid_molecules: bool = True,
) -> dict:
    tolerance = 1.0e-7
    bounds = [
        (
            min(point[axis] for point in coordinates),
            max(point[axis] for point in coordinates),
        )
        for axis in range(3)
    ]
    spans = [upper - lower for lower, upper in bounds]
    lower_clearance = [lower for lower, _ in bounds]
    upper_clearance = [cell - upper for _, upper in bounds]
    minimum_face_clearance = min(lower_clearance + upper_clearance)
    periodic_image_lower_bound = min(cell - span for span in spans)
    contact = _closest_intermolecular_pair(coordinates, symbols, instances)
    nearest = contact["minimum_distance_angstrom"]
    farthest_nearest = _maximum_nearest_molecule_distance(coordinates, instances)
    rigidity_error = _maximum_intramolecular_distance_deviation(
        coordinates, instances, templates
    )
    expected_atom_count = sum(
        template.count * len(template.symbols) for template in templates
    )
    checks = {
        "finite_coordinates": all(
            math.isfinite(value) for point in coordinates for value in point
        ),
        "positive_cell": math.isfinite(cell) and cell > 0.0,
        "atom_count_matches": len(coordinates) == len(symbols) == expected_atom_count,
        "all_atoms_inside_cell": all(
            -tolerance <= value <= cell + tolerance
            for point in coordinates
            for value in point
        ),
        "requested_face_vacuum_present": (
            minimum_face_clearance + tolerance >= requested_vacuum
        ),
        "periodic_image_clearance_present": (
            periodic_image_lower_bound + tolerance >= 2.0 * requested_vacuum
        ),
        "intermolecular_distance_safe": nearest + tolerance >= minimum_distance,
        "oxygen_oxygen_distance_safe": (
            minimum_oxygen_oxygen_distance is None
            or contact["minimum_oxygen_oxygen_distance_angstrom"] is None
            or contact["minimum_oxygen_oxygen_distance_angstrom"] + tolerance
            >= minimum_oxygen_oxygen_distance
        ),
        "all_molecules_associated": (
            maximum_nearest_molecule_distance is None
            or farthest_nearest <= maximum_nearest_molecule_distance + tolerance
        ),
        "intramolecular_geometry_preserved": (
            not require_rigid_molecules or rigidity_error <= 1.0e-5
        ),
        "poscar_element_order_correct": written_element_order == expected_elements,
        "potcar_order_compatible": reference_sanity["potcar_status"] in {
            "PASS",
            "NOT_SUPPLIED",
        },
    }
    return {
        "schema_version": 1,
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "atom_count": len(coordinates),
        "element_order": written_element_order,
        "expected_element_order": expected_elements,
        "cell_angstrom": [cell, cell, cell],
        "coordinate_bounds_angstrom": {
            axis: [bounds[index][0], bounds[index][1]]
            for index, axis in enumerate("xyz")
        },
        "span_angstrom": {axis: spans[index] for index, axis in enumerate("xyz")},
        "lower_face_clearance_angstrom": {
            axis: lower_clearance[index] for index, axis in enumerate("xyz")
        },
        "upper_face_clearance_angstrom": {
            axis: upper_clearance[index] for index, axis in enumerate("xyz")
        },
        "requested_vacuum_angstrom": requested_vacuum,
        "minimum_face_clearance_angstrom": minimum_face_clearance,
        "periodic_image_separation_lower_bound_angstrom": periodic_image_lower_bound,
        "minimum_intermolecular_distance_angstrom": nearest,
        "closest_intermolecular_pair": contact["closest_pair"],
        "minimum_intermolecular_oxygen_oxygen_distance_angstrom": contact[
            "minimum_oxygen_oxygen_distance_angstrom"
        ],
        "closest_intermolecular_oxygen_oxygen_pair": contact[
            "closest_oxygen_oxygen_pair"
        ],
        "required_minimum_intermolecular_distance_angstrom": minimum_distance,
        "required_minimum_intermolecular_oxygen_oxygen_distance_angstrom": (
            minimum_oxygen_oxygen_distance
        ),
        "maximum_nearest_molecule_distance_angstrom": farthest_nearest,
        "required_maximum_nearest_molecule_distance_angstrom": (
            maximum_nearest_molecule_distance
        ),
        "maximum_intramolecular_distance_deviation_angstrom": rigidity_error,
        "vasp_reference": reference_sanity,
    }


def _write_sanity_report(case: Path, report: dict) -> None:
    (case / "sanity_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        report["status"],
        f"atoms={report['atom_count']}",
        f"elements={' '.join(report['element_order'])}",
        f"cell_angstrom={report['cell_angstrom'][0]:.8f}",
        (
            "minimum_face_clearance_angstrom="
            f"{report['minimum_face_clearance_angstrom']:.8f}"
        ),
        (
            "periodic_image_separation_lower_bound_angstrom="
            f"{report['periodic_image_separation_lower_bound_angstrom']:.8f}"
        ),
        (
            "minimum_intermolecular_distance_angstrom="
            f"{report['minimum_intermolecular_distance_angstrom']:.8f}"
        ),
        (
            "minimum_intermolecular_oxygen_oxygen_distance_angstrom="
            + (
                f"{report['minimum_intermolecular_oxygen_oxygen_distance_angstrom']:.8f}"
                if report["minimum_intermolecular_oxygen_oxygen_distance_angstrom"]
                is not None
                else "NOT_PRESENT"
            )
        ),
        (
            "maximum_nearest_molecule_distance_angstrom="
            f"{report['maximum_nearest_molecule_distance_angstrom']:.8f}"
        ),
        (
            "maximum_intramolecular_distance_deviation_angstrom="
            f"{report['maximum_intramolecular_distance_deviation_angstrom']:.10e}"
        ),
        f"potcar_status={report['vasp_reference']['potcar_status']}",
    ]
    lines.extend(
        f"check.{name}={'PASS' if passed else 'FAIL'}"
        for name, passed in report["checks"].items()
    )
    (case / "sanity_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _ordered_atoms(
    symbols: list[str],
    coordinates: list[tuple[float, float, float]],
    instances: list[MoleculeInstance],
) -> list[dict]:
    identity: list[tuple[str, int, int]] = [("", 0, 0)] * len(symbols)
    for item in instances:
        for local_index, atom_index in enumerate(range(item.start, item.stop), 1):
            identity[atom_index] = (item.slug, item.molecule_id, local_index)
    order = sorted(range(len(symbols)), key=lambda index: (symbols[index], index))
    return [
        {
            "element": symbols[index],
            "coordinate": coordinates[index],
            "source_atom_index": index + 1,
            "molecule_slug": identity[index][0],
            "molecule_id": identity[index][1],
            "local_atom_index": identity[index][2],
        }
        for index in order
    ]


def _write_poscar(atoms: list[dict], cell: float, path: Path, title: str) -> None:
    element_names = sorted({atom["element"] for atom in atoms})
    counts = [sum(atom["element"] == element for atom in atoms) for element in element_names]
    lines = [
        title,
        "1.0",
        f"{cell:.10f} 0.0000000000 0.0000000000",
        f"0.0000000000 {cell:.10f} 0.0000000000",
        f"0.0000000000 0.0000000000 {cell:.10f}",
        " ".join(element_names),
        " ".join(str(count) for count in counts),
        "Cartesian",
    ]
    lines.extend(
        f"{x:.10f} {y:.10f} {z:.10f}"
        for atom in atoms
        for x, y, z in (atom["coordinate"],)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_xyz(atoms: list[dict], path: Path, comment: str) -> None:
    rows = [
        f"{atom['element']} {atom['coordinate'][0]:.10f} "
        f"{atom['coordinate'][1]:.10f} {atom['coordinate'][2]:.10f}"
        for atom in atoms
    ]
    path.write_text(f"{len(rows)}\n{comment}\n" + "\n".join(rows) + "\n", encoding="utf-8")


def _write_raw_xyz(
    symbols: list[str],
    coordinates: list[tuple[float, float, float]],
    path: Path,
    comment: str,
) -> None:
    rows = [
        f"{symbol} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}"
        for symbol, point in zip(symbols, coordinates)
    ]
    path.write_text(
        f"{len(rows)}\n{comment}\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )


def _replace_incar_tags(path: Path, updates: dict[str, int | float]) -> None:
    text = path.read_text(encoding="utf-8")
    keys = {key.upper() for key in updates}
    kept = [
        line
        for line in text.splitlines()
        if not (
            "=" in line
            and line.split("=", 1)[0].strip().upper() in keys
        )
    ]
    kept.extend(["", "# Generated agglomeration MD schedule"])
    kept.extend(f"{key} = {value}" for key, value in updates.items())
    path.write_text("\n".join(kept).rstrip() + "\n", encoding="utf-8")


def _set_slurm_job_name(path: Path, job_name: str) -> None:
    if not path.is_file():
        return
    text = path.read_text(encoding="utf-8")
    replacement = f'#SBATCH -J "{job_name}"'
    if re.search(r"^#SBATCH\s+(?:-J|--job-name(?:=|\s))", text, re.MULTILINE):
        text = re.sub(
            r"^#SBATCH\s+(?:-J\s+.*|--job-name(?:=|\s).*)$",
            replacement,
            text,
            count=1,
            flags=re.MULTILINE,
        )
    else:
        lines = text.splitlines()
        lines.insert(1, replacement)
        text = "\n".join(lines) + "\n"
    path.write_text(text, encoding="utf-8")


def _prepare_vasp_schedule(
    case: Path,
    reference: Path,
    poscar: Path,
    *,
    job_prefix: str,
    hold_steps: int,
    heating_steps: int,
) -> list[dict]:
    runs = [
        (case / "300K", 300, 300, hold_steps, f"{job_prefix}_300K"),
        (
            case / "400K/01_heat_300_to_400K",
            300,
            400,
            heating_steps,
            f"{job_prefix}_heat",
        ),
        (
            case / "400K/02_hold_400K",
            400,
            400,
            hold_steps,
            f"{job_prefix}_400K",
        ),
    ]
    records: list[dict] = []
    for directory, start, end, steps, job_name in runs:
        directory.mkdir(parents=True, exist_ok=True)
        copied = _copy_reference(reference, directory)
        shutil.copy2(poscar, directory / "POSCAR")
        _replace_incar_tags(
            directory / "INCAR",
            {"IBRION": 0, "NSW": steps, "TEBEG": start, "TEEND": end},
        )
        _set_slurm_job_name(directory / "runvasp.sh", job_name[:64])
        records.append(
            {
                "path": str(directory.relative_to(case)),
                "temperature_start_K": start,
                "temperature_end_K": end,
                "steps": steps,
                "reference_files": copied,
            }
        )
    submit = case / "submit_temperature_jobs.sh"
    submit.write_text(
        """#!/bin/bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
job300=$(cd "$root/300K" && sbatch --parsable runvasp.sh)
jobheat=$(cd "$root/400K/01_heat_300_to_400K" && sbatch --parsable runvasp.sh)
job400=$(cd "$root/400K/02_hold_400K" && sbatch --parsable "--dependency=afterok:${jobheat}" runvasp.sh)
printf 'submitted 300K=%s heat=%s 400K-hold=%s\n' "$job300" "$jobheat" "$job400"
""",
        encoding="utf-8",
    )
    submit.chmod(0o755)
    hold_launcher = case / "400K/02_hold_400K/runvasp.sh"
    if hold_launcher.is_file():
        text = hold_launcher.read_text(encoding="utf-8")
        marker = "set -e\ncp ../01_heat_300_to_400K/CONTCAR POSCAR\n"
        lines = text.splitlines(keepends=True)
        insert_at = 1
        for index, line in enumerate(lines[1:], 1):
            if line.startswith("#SBATCH") or not line.strip():
                insert_at = index + 1
                continue
            break
        lines.insert(insert_at, marker)
        hold_launcher.write_text("".join(lines), encoding="utf-8")
    return records


def _write_xtb_batch(output: Path, cases: list[dict], settings: dict) -> None:
    case_list = output / "xtb_cases.tsv"
    case_list.write_text(
        "".join(
            f"{item['path']}\t{item['charge']}\t{item['uhf']}\n" for item in cases
        ),
        encoding="utf-8",
    )
    account = str(settings.get("account", "loni_perovsk27"))
    partition = str(settings.get("partition", "single"))
    cpus = int(settings.get("cpus", 1))
    walltime = str(settings.get("walltime", "12:00:00"))
    gfn = int(settings.get("gfn", 2))
    opt_level = str(settings.get("opt_level", "loose"))
    max_cycles = int(settings.get("max_cycles", 100))
    concurrency = int(settings.get("array_concurrency", 4))
    md_enabled = bool(settings.get("md_enabled", True))
    campaign_root = shlex.quote(str(output.resolve()))
    extractor = output / "extract_last_xtb_frame.py"
    extractor.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path

source = Path(sys.argv[1])
destination = Path(sys.argv[2])
lines = source.read_text(encoding="utf-8").splitlines()
frames = []
index = 0
while index < len(lines):
    while index < len(lines) and not lines[index].strip():
        index += 1
    if index >= len(lines):
        break
    try:
        atoms = int(lines[index].strip())
    except ValueError as error:
        raise SystemExit(f"{source}:{index + 1}: invalid XYZ atom count") from error
    stop = index + atoms + 2
    if stop > len(lines):
        raise SystemExit(f"{source}: incomplete final XYZ frame")
    frame = lines[index:stop]
    if any(len(row.split()) < 4 for row in frame[2:]):
        raise SystemExit(f"{source}: malformed XYZ coordinate row")
    frames.append(frame)
    index = stop
if not frames:
    raise SystemExit(f"{source}: no complete XYZ frames")
destination.write_text("\\n".join(frames[-1]) + "\\n", encoding="utf-8")
""",
        encoding="utf-8",
    )
    extractor.chmod(0o755)
    steering_writer = output / "write_xtb_steering_input.py"
    steering_writer.write_text(
        """#!/usr/bin/env python3
import json
import math
import sys
from pathlib import Path

source = Path(sys.argv[1])
plan_path = Path(sys.argv[2])
input_path = Path(sys.argv[3])
record_path = Path(sys.argv[4])
lines = source.read_text(encoding="utf-8").splitlines()
atoms = int(lines[0])
if len(lines) != atoms + 2:
    raise SystemExit(f"{source}: invalid single-frame XYZ")
coordinates = [tuple(map(float, row.split()[1:4])) for row in lines[2:]]
plan = json.loads(plan_path.read_text(encoding="utf-8"))

def angle(first, vertex, third):
    left = tuple(a - b for a, b in zip(first, vertex))
    right = tuple(a - b for a, b in zip(third, vertex))
    denominator = math.sqrt(sum(value * value for value in left)) * math.sqrt(
        sum(value * value for value in right)
    )
    cosine = sum(a * b for a, b in zip(left, right)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

realized = []
constraints = []
for pair in plan["pairs"]:
    left = pair["left_atom_index"] - 1
    right = pair["right_atom_index"] - 1
    left_carbon = pair["left_carbon_atom_index"] - 1
    right_carbon = pair["right_carbon_atom_index"] - 1
    current = math.dist(coordinates[left], coordinates[right])
    target = min(
        current,
        max(
            plan["minimum_target_pp_distance_angstrom"],
            current - plan["maximum_distance_reduction_angstrom"],
        ),
    )
    realized_pair = pair | {
        "post_optimization_distance_angstrom": current,
        "restraint_distance_angstrom": target,
        "left_axis_angle_degrees": angle(
            coordinates[left_carbon], coordinates[left], coordinates[right]
        ),
        "right_axis_angle_degrees": angle(
            coordinates[right_carbon], coordinates[right], coordinates[left]
        ),
        "restraint_axis_angle_degrees": plan["target_axis_angle_degrees"],
    }
    realized.append(realized_pair)
    constraints.append(
        f"  distance: {left + 1},{right + 1},{target:.8f}"
    )
    constraints.append(
        f"  angle: {left_carbon + 1},{left + 1},{right + 1},"
        f"{plan['target_axis_angle_degrees']:g}"
    )
    constraints.append(
        f"  angle: {right_carbon + 1},{right + 1},{left + 1},"
        f"{plan['target_axis_angle_degrees']:g}"
    )
input_path.write_text(
    f'''$seed {plan["seed"]}
$constrain
  force constant={plan["force_constant_atomic_units"]:g}
{chr(10).join(constraints)}
$end
$md
  temp={plan["temperature_K"]:g}
  time={plan["time_ps"]:g}
  dump={plan["dump_fs"]:g}
  step={plan["step_fs"]:g}
  nvt=1
  hmass=0
  shake=0
  sccacc={plan["sccacc"]:g}
$end
''',
    encoding="utf-8",
)
record_path.write_text(
    json.dumps({"pairs": realized}, indent=2) + "\\n", encoding="utf-8"
)
""",
        encoding="utf-8",
    )
    steering_writer.chmod(0o755)
    if md_enabled:
        workflow = f'''if [[ -s xtbopt.xyz && ! -s xtbopt_initial.xyz ]]; then
  mv xtbopt.xyz xtbopt_initial.xyz
  [[ ! -s xtbopt.log ]] || mv xtbopt.log xtbopt_initial.log
fi
if [[ ! -s xtbopt_initial.xyz ]]; then
  rm -f xtbopt.xyz
  "$XTB_EXE" input.xyz --gfn {gfn} --chrg "$charge" --uhf "$uhf" --opt {opt_level} --cycles {max_cycles} > xtb_initial_opt.out 2> xtb_initial_opt.err
  test -s xtbopt.xyz
  mv xtbopt.xyz xtbopt_initial.xyz
  [[ ! -s xtbopt.log ]] || mv xtbopt.log xtbopt_initial.log
fi
protocol_changed=false
if [[ -s xtb_stages.protocol ]]; then
  if ! cmp -s xtb_protocol.sha256 xtb_stages.protocol; then
    protocol_changed=true
  fi
elif [[ -s xtbsteered_last.xyz || -s xtbmd_unbiased_last.xyz || -s xtbfinal.xyz ]]; then
  protocol_changed=true
fi
if [[ "$protocol_changed" == true ]]; then
  old_tag=legacy
  if [[ -s xtb_stages.protocol ]]; then
    old_tag=$(head -c 12 xtb_stages.protocol)
  fi
  for artifact in md_steered.inp steering_restraints.json xtb_steered.trj xtbsteered_last.xyz xtb_unbiased.trj xtbmd_unbiased_last.xyz xtbfinal.xyz xtbfinal.log; do
    if [[ -e "$artifact" ]]; then
      archived="${{artifact}}.protocol_${{old_tag}}"
      if [[ -e "$archived" ]]; then
        rm -f "$artifact"
      else
        mv "$artifact" "$archived"
      fi
    fi
  done
  rm -f xtb_protocol.complete
fi
cp xtb_protocol.sha256 xtb_stages.protocol
md_start=xtbopt_initial.xyz
if [[ -s steering_plan.json ]]; then
  if [[ ! -s xtbsteered_last.xyz ]]; then
    "$PYTHON_EXE" "$root/write_xtb_steering_input.py" xtbopt_initial.xyz steering_plan.json md_steered.inp steering_restraints.json
    rm -f xtb.trj mdrestart xtbrestart xtbmdok
    "$XTB_EXE" xtbopt_initial.xyz --gfn {gfn} --chrg "$charge" --uhf "$uhf" --md --input md_steered.inp > xtb_steered_md.out 2> xtb_steered_md.err
    test -s xtb.trj
    "$PYTHON_EXE" "$root/extract_last_xtb_frame.py" xtb.trj xtbsteered_last.xyz
    mv xtb.trj xtb_steered.trj
    test -s xtbsteered_last.xyz
  fi
  md_start=xtbsteered_last.xyz
fi
if [[ ! -s xtbmd_unbiased_last.xyz ]]; then
  rm -f xtb.trj mdrestart xtbrestart xtbmdok
  "$XTB_EXE" "$md_start" --gfn {gfn} --chrg "$charge" --uhf "$uhf" --md --input md.inp > xtb_unbiased_md.out 2> xtb_unbiased_md.err
  test -s xtb.trj
  "$PYTHON_EXE" "$root/extract_last_xtb_frame.py" xtb.trj xtbmd_unbiased_last.xyz
  mv xtb.trj xtb_unbiased.trj
  test -s xtbmd_unbiased_last.xyz
fi
if [[ ! -s xtb_protocol.complete ]] || ! cmp -s xtb_protocol.sha256 xtb_protocol.complete; then
  if [[ -s xtbfinal.xyz && ! -e xtbfinal.pre_staged_md.xyz ]]; then
    mv xtbfinal.xyz xtbfinal.pre_staged_md.xyz
  else
    rm -f xtbfinal.xyz
  fi
  rm -f xtbopt.xyz
  "$XTB_EXE" xtbmd_unbiased_last.xyz --gfn {gfn} --chrg "$charge" --uhf "$uhf" --opt {opt_level} --cycles {max_cycles} > xtb_final_opt.out 2> xtb_final_opt.err
  test -s xtbopt.xyz
  mv xtbopt.xyz xtbfinal.xyz
  [[ ! -s xtbopt.log ]] || mv xtbopt.log xtbfinal.log
  cp xtb_protocol.sha256 xtb_protocol.complete
fi
test -s xtbfinal.xyz
cmp -s xtb_protocol.sha256 xtb_protocol.complete'''
        completion_test = (
            '[[ -s xtbfinal.xyz && -s xtb_protocol.complete ]] '
            '&& cmp -s xtb_protocol.sha256 xtb_protocol.complete'
        )
    else:
        workflow = f'''"$XTB_EXE" input.xyz --gfn {gfn} --chrg "$charge" --uhf "$uhf" --opt {opt_level} --cycles {max_cycles} > xtb.out 2> xtb.err
test -s xtbopt.xyz'''
        completion_test = "[[ -s xtbopt.xyz ]]"
    script = output / "run_xtb_array.sbatch"
    script.write_text(
        f"""#!/bin/bash
#SBATCH -p {partition}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c {cpus}
#SBATCH -t {walltime}
#SBATCH -A {account}
#SBATCH -J "me4pacz_xtb"
#SBATCH -o xtb.%A_%a.out
#SBATCH --array=0-{len(cases) - 1}%{concurrency}

set -euo pipefail
root={campaign_root}
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$root/xtb_cases.tsv")
IFS=$'\\t' read -r relative charge uhf <<< "$line"
work="$root/$relative"
cd "$work"
if {completion_test}; then
  echo "already complete: $relative"
  exit 0
fi

: "${{XTB_ENV:=/project/lgutsev/env/xtb_env}}"
if [[ -z "${{XTB_EXE:-}}" ]]; then
  if [[ -x "$XTB_ENV/bin/xtb" ]]; then
    XTB_EXE="$XTB_ENV/bin/xtb"
  else
    XTB_EXE=xtb
  fi
fi
if ! command -v "$XTB_EXE" >/dev/null 2>&1; then
  : "${{CONDA_SH:=$HOME/miniconda3/etc/profile.d/conda.sh}}"
  if [[ -f "$CONDA_SH" ]]; then
    source "$CONDA_SH"
    conda activate "$XTB_ENV"
  fi
fi
if ! command -v "$XTB_EXE" >/dev/null 2>&1; then
  echo "xTB not found; set XTB_EXE or provide a valid XTB_ENV" >&2
  exit 127
fi
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-{cpus}}}
export MKL_NUM_THREADS=1
ulimit -s unlimited
if [[ -x "$XTB_ENV/bin/python" ]]; then
  PYTHON_EXE="$XTB_ENV/bin/python"
else
  PYTHON_EXE=python3
fi
{workflow}
""",
        encoding="utf-8",
    )
    script.chmod(0o755)


def _md_seed(seed: int, offset: int) -> int:
    resolved = (seed + offset) % 899999999
    return resolved or 1


def _head_bias_for_replica(replica: int, replicas: int, settings: dict) -> bool:
    if not bool(settings.get("head_bias_enabled", True)):
        return False
    biased = int(settings.get("head_bias_replicas_per_family", 1))
    if biased < 0:
        raise ValueError("head_bias_replicas_per_family cannot be negative")
    return replica < min(biased, replicas)


def _write_xtb_md_inputs(
    directory: Path,
    settings: dict,
    seed: int,
    p_pairs: list[dict],
    steered: bool,
) -> dict:
    temperature = float(settings.get("md_temperature_K", 400.0))
    total_time_ps = float(settings.get("md_time_ps", 10.0))
    steering_time_ps = (
        float(settings.get("head_bias_time_ps", 4.0)) if steered else 0.0
    )
    unbiased_time_ps = total_time_ps - steering_time_ps
    dump_fs = float(settings.get("md_dump_fs", 50.0))
    step_fs = float(settings.get("md_step_fs", 1.0))
    sccacc = float(settings.get("md_sccacc", 1.0))
    target = float(settings.get("head_bias_target_pp_distance_angstrom", 4.5))
    maximum_reduction = float(
        settings.get("head_bias_max_distance_reduction_angstrom", 1.5)
    )
    target_axis_angle = float(
        settings.get("head_bias_target_axis_angle_degrees", 170.0)
    )
    force_constant = float(settings.get("head_bias_force_constant", 0.005))
    if min(temperature, total_time_ps, dump_fs, step_fs, sccacc) <= 0:
        raise ValueError("xTB MD temperature, time, dump, step, and sccacc must be positive")
    if steered and (not p_pairs or steering_time_ps <= 0 or unbiased_time_ps <= 0):
        raise ValueError(
            "steered xTB MD requires P pairs and 0 < head_bias_time_ps < md_time_ps"
        )
    if steered and (
        min(target, maximum_reduction, force_constant) <= 0
        or not 0 < target_axis_angle < 180
    ):
        raise ValueError(
            "xTB steering distance parameters and force constant must be positive, "
            "and the target axis angle must lie between 0 and 180 degrees"
        )

    def md_text(md_seed: int, time_ps: float) -> str:
        return f"""$seed {md_seed}
$md
  temp={temperature:g}
  time={time_ps:g}
  dump={dump_fs:g}
  step={step_fs:g}
  nvt=1
  hmass=0
  shake=0
  sccacc={sccacc:g}
$end
"""

    (directory / "md.inp").write_text(
        md_text(_md_seed(seed, 104729 if steered else 0), unbiased_time_ps),
        encoding="utf-8",
    )
    if steered:
        steering_plan = {
            "seed": _md_seed(seed, 0),
            "temperature_K": temperature,
            "time_ps": steering_time_ps,
            "dump_fs": dump_fs,
            "step_fs": step_fs,
            "sccacc": sccacc,
            "minimum_target_pp_distance_angstrom": target,
            "maximum_distance_reduction_angstrom": maximum_reduction,
            "target_axis_angle_degrees": target_axis_angle,
            "force_constant_atomic_units": force_constant,
            "pairs": p_pairs,
        }
        (directory / "steering_plan.json").write_text(
            json.dumps(steering_plan, indent=2) + "\n", encoding="utf-8"
        )
    else:
        (directory / "steering_plan.json").unlink(missing_ok=True)
        (directory / "md_steered.inp").unlink(missing_ok=True)
        (directory / "steering_restraints.json").unlink(missing_ok=True)

    protocol = {
        "schema_version": 3,
        "total_md_time_ps": total_time_ps,
        "temperature_K": temperature,
        "step_fs": step_fs,
        "dump_fs": dump_fs,
        "sccacc": sccacc,
        "head_bias": {
            "enabled": steered,
            "time_ps": steering_time_ps,
            "unbiased_time_ps": unbiased_time_ps,
            "minimum_target_pp_distance_angstrom": target if steered else None,
            "maximum_distance_reduction_angstrom": (
                maximum_reduction if steered else None
            ),
            "target_axis_angle_degrees": target_axis_angle if steered else None,
            "force_constant_atomic_units": force_constant if steered else None,
            "pairs": p_pairs if steered else [],
        },
        "final_stage": "unbiased GFN-xTB optimization",
    }
    protocol_path = directory / "xtb_protocol.json"
    protocol_path.write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    (directory / "xtb_protocol.sha256").write_text(
        _sha256(protocol_path) + "\n", encoding="utf-8"
    )
    return protocol


def _copy_reference(reference: Path | None, destination: Path) -> list[dict]:
    if reference is None:
        return []
    if not reference.is_dir():
        raise FileNotFoundError(f"VASP reference directory not found: {reference}")
    copied: list[dict] = []
    for source in sorted(reference.iterdir()):
        if not source.is_file() or source.name in _VASP_OUTPUTS:
            continue
        target = destination / source.name
        shutil.copy2(source, target)
        copied.append({"name": source.name, "sha256": _sha256(target)})
    return copied


def prepare_agglomeration(
    config_path: Path,
    output: Path,
    *,
    reference_dir: Path | None = None,
    packed_xyz: Path | None = None,
    structures_only: bool = False,
) -> Path:
    """Prepare size-stratified molecular clusters as VASP run packages."""
    if reference_dir is not None and structures_only:
        raise ValueError("--reference-dir and --structures-only are mutually exclusive")
    manifest_path = output / "agglomeration_manifest.json"
    resume = False
    if output.exists() and any(output.iterdir()):
        if manifest_path.is_file():
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
            resume = (
                previous.get("status") in {"PACKMOL_REQUIRED", "XTB_REQUIRED"}
                and previous.get("config_sha256") == _sha256(config_path)
            )
        if not resume:
            raise FileExistsError(f"refusing to replace non-empty output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    specs, grouped = _agglomerate_specs(config)
    xtb_settings = config.get("xtb", {})
    xtb_enabled = bool(xtb_settings.get("enabled", False))
    xtb_md_enabled = bool(xtb_settings.get("md_enabled", True))
    head_bias_enabled = bool(xtb_settings.get("head_bias_enabled", True))
    head_bias_replicas = int(
        xtb_settings.get("head_bias_replicas_per_family", 1)
    )
    if head_bias_replicas < 0:
        raise ValueError("head_bias_replicas_per_family cannot be negative")
    md_settings = config.get("vasp_md", {})
    hold_steps = int(md_settings.get("hold_steps", 3000))
    heating_steps = int(md_settings.get("heating_steps", 1000))
    if hold_steps < 1 or heating_steps < 1:
        raise ValueError("vasp_md hold_steps and heating_steps must be positive")
    total_replicas = sum(spec.replicas for spec in specs)
    if packed_xyz is not None and total_replicas != 1:
        raise ValueError(
            "--packed-xyz can only be used when the complete campaign has one replica"
        )

    template_dir = output / "templates"
    template_dir.mkdir(exist_ok=True)
    template_paths: dict[str, Path] = {}
    unique_templates: dict[str, MoleculeTemplate] = {}
    for spec in specs:
        for template in spec.templates:
            unique_templates.setdefault(template.slug, template)
    for template in unique_templates.values():
        path = template_dir / f"{template.slug}.xyz"
        _write_template_xyz(template, path)
        template_paths[template.slug] = path
    expected_elements = sorted(
        {
            symbol
            for template in unique_templates.values()
            for symbol in template.symbols
        }
    )
    reference_validation = _reference_sanity(
        reference_dir, expected_elements, structures_only=structures_only
    )
    case_root_name = "structures" if structures_only else "vasp_runs"
    executable = shutil.which("packmol")
    index_rows: list[dict] = []
    sanity_rows: list[dict] = []
    replica_records: list[dict] = []
    xtb_cases: list[dict] = []
    status = "COMPLETE"

    for spec in specs:
        expected_symbols, instances = _expected_layout(list(spec.templates))
        family_packmol = output / "packmol"
        family_cases = output / case_root_name
        if grouped:
            family_packmol /= spec.name
            family_cases /= spec.name
        for replica in range(spec.replicas):
            seed = spec.base_seed + replica
            work = family_packmol / f"replica_{replica:03d}"
            work.mkdir(parents=True, exist_ok=True)
            local_templates: dict[str, Path] = {}
            for template in spec.templates:
                source = template_paths[template.slug]
                destination = work / source.name
                shutil.copy2(source, destination)
                local_templates[template.slug] = destination
            packed = work / "packed.xyz"
            packmol_input = _packmol_input(
                list(spec.templates),
                local_templates,
                packed.name,
                seed,
                spec.tolerance,
                spec.radius,
            )
            (work / "packmol.inp").write_text(packmol_input, encoding="utf-8")
            if packed_xyz is not None:
                shutil.copy2(packed_xyz, packed)
            elif packed.exists():
                pass
            elif executable:
                with (work / "packmol.inp").open("rb") as source:
                    run = subprocess.run(
                        [executable],
                        cwd=work,
                        stdin=source,
                        capture_output=True,
                        text=False,
                        check=False,
                    )
                (work / "packmol.stdout").write_bytes(run.stdout)
                (work / "packmol.stderr").write_bytes(run.stderr)
                if run.returncode:
                    raise RuntimeError(
                        f"Packmol failed for {spec.name} replica {replica}: "
                        + run.stderr.decode(errors="replace")
                    )
            else:
                status = "PACKMOL_REQUIRED"
                replica_records.append(
                    {
                        "agglomerate": spec.name,
                        "replica": replica,
                        "seed": seed,
                        "status": "PACKMOL_REQUIRED",
                        "command": "packmol < packmol.inp",
                        "working_directory": str(work.relative_to(output)),
                    }
                )
                continue

            symbols, coordinates = _read_xyz(packed)
            if symbols != expected_symbols:
                raise ValueError(
                    f"{packed}: atom ordering differs from the configured "
                    "LigParGen templates"
                )
            coordinates, compact_scale = _compact_molecule_centers(
                coordinates,
                symbols,
                instances,
                spec.compact_to_distance,
            )
            uncentered_variants = {
                scale: _scale_molecule_centers(coordinates, instances, scale)
                for scale in spec.scales
            }
            fixed_unrelaxed_cell = _fixed_cell(uncentered_variants, spec.vacuum)
            structures: list[dict] = []
            for variant_index, scale in enumerate(spec.scales):
                case_name = (
                    f"r{replica:03d}_s{variant_index:02d}_{scale:.3f}"
                    .replace(".", "p")
                )
                raw_variant = uncentered_variants[scale]
                xtb_relative: Path | None = None
                xtb_protocol: dict | None = None
                if xtb_enabled:
                    xtb_work = output / "xtb"
                    if grouped:
                        xtb_work /= spec.name
                    xtb_work /= case_name
                    xtb_work.mkdir(parents=True, exist_ok=True)
                    xtb_relative = xtb_work.relative_to(output)
                    xtb_input = xtb_work / "input.xyz"
                    _write_raw_xyz(
                        symbols,
                        raw_variant,
                        xtb_input,
                        f"Packmol seed {seed}; adaptive center scale {compact_scale:.8f}",
                    )
                    if xtb_md_enabled:
                        p_pairs = _nearest_phosphonate_pairs(
                            raw_variant, instances
                        )
                        steered = _head_bias_for_replica(
                            replica, spec.replicas, xtb_settings
                        )
                        xtb_protocol = _write_xtb_md_inputs(
                            xtb_work,
                            xtb_settings,
                            seed + variant_index,
                            p_pairs,
                            steered,
                        )
                        realized_restraints = xtb_work / "steering_restraints.json"
                        if steered and realized_restraints.is_file():
                            xtb_protocol["head_bias"]["realized_pairs"] = json.loads(
                                realized_restraints.read_text(encoding="utf-8")
                            )["pairs"]
                    total_charge = sum(
                        template.count * template.net_charge
                        for template in spec.templates
                    )
                    rounded_charge = round(total_charge)
                    if abs(total_charge - rounded_charge) > 1.0e-2:
                        raise ValueError(
                            f"{spec.name}: xTB total charge {total_charge} is not integral"
                        )
                    xtb_cases.append(
                        {
                            "path": str(xtb_relative),
                            "charge": rounded_charge,
                            "uhf": int(xtb_settings.get("uhf", 0)),
                        }
                    )
                    optimized = xtb_work / (
                        "xtbfinal.xyz" if xtb_md_enabled else "xtbopt.xyz"
                    )
                    protocol_ready = True
                    if xtb_md_enabled:
                        expected_protocol = xtb_work / "xtb_protocol.sha256"
                        completed_protocol = xtb_work / "xtb_protocol.complete"
                        protocol_ready = (
                            completed_protocol.is_file()
                            and completed_protocol.read_text(encoding="utf-8")
                            == expected_protocol.read_text(encoding="utf-8")
                        )
                    if not optimized.is_file() or not protocol_ready:
                        status = "XTB_REQUIRED"
                        structures.append(
                            {
                                "case": case_name,
                                "status": "XTB_REQUIRED",
                                "xtb_working_directory": str(xtb_relative),
                            }
                        )
                        continue
                    optimized_symbols, raw_variant = _read_xyz(optimized)
                    if optimized_symbols != expected_symbols:
                        raise ValueError(
                            f"{optimized}: xTB changed atom ordering or element identities"
                        )
                cell = (
                    _fixed_cell({1.0: raw_variant}, spec.vacuum)
                    if xtb_enabled
                    else fixed_unrelaxed_cell
                )
                variant = _center_in_cell(raw_variant, cell)
                case = family_cases / case_name
                case.mkdir(parents=True, exist_ok=resume)
                ordered = _ordered_atoms(symbols, variant, instances)
                title = (
                    f"phosphonate agglomeration family={spec.name} "
                    f"replica={replica} center_scale={scale:g}"
                )
                _write_poscar(ordered, cell, case / "POSCAR", title)
                _write_xyz(ordered, case / "structure.xyz", title)
                written_element_order = (case / "POSCAR").read_text(
                    encoding="utf-8"
                ).splitlines()[5].split()
                atom_map = [
                    {
                        "poscar_index": index,
                        "element": atom["element"],
                        "source_atom_index": atom["source_atom_index"],
                        "molecule_slug": atom["molecule_slug"],
                        "molecule_id": atom["molecule_id"],
                        "local_atom_index": atom["local_atom_index"],
                    }
                    for index, atom in enumerate(ordered, 1)
                ]
                copied = []
                sanity = _structure_sanity(
                    variant,
                    symbols,
                    instances,
                    list(spec.templates),
                    cell=cell,
                    requested_vacuum=spec.vacuum,
                    minimum_distance=(
                        float(xtb_settings.get("minimum_distance_angstrom", 0.70))
                        if xtb_enabled
                        else spec.minimum_distance
                    ),
                    reference_sanity=reference_validation,
                    expected_elements=expected_elements,
                    written_element_order=written_element_order,
                    minimum_oxygen_oxygen_distance=(
                        float(
                            xtb_settings.get(
                                "minimum_oxygen_oxygen_distance_angstrom", 2.20
                            )
                        )
                        if xtb_enabled
                        else None
                    ),
                    maximum_nearest_molecule_distance=(
                        float(
                            xtb_settings.get(
                                "maximum_nearest_molecule_distance_angstrom", 6.0
                            )
                        )
                        if xtb_enabled
                        else None
                    ),
                    require_rigid_molecules=not xtb_enabled,
                )
                _write_sanity_report(case, sanity)
                if sanity["status"] != "PASS":
                    failed = [
                        name for name, passed in sanity["checks"].items() if not passed
                    ]
                    raise ValueError(
                        f"{case}: structural sanity check failed: "
                        f"{', '.join(failed)}; see sanity_report.json"
                    )
                nearest = sanity["minimum_intermolecular_distance_angstrom"]
                composition = [
                    {"slug": template.slug, "count": template.count}
                    for template in spec.templates
                ]
                schedules = []
                if not structures_only and xtb_enabled:
                    schedules = _prepare_vasp_schedule(
                        case,
                        reference_dir,
                        case / "POSCAR",
                        job_prefix=f"{spec.name}_r{replica:03d}",
                        hold_steps=hold_steps,
                        heating_steps=heating_steps,
                    )
                elif not structures_only:
                    copied = _copy_reference(reference_dir, case)
                case_manifest = {
                    "schema_version": 4,
                    "generator": "nio-md-prep prepare-agglomeration",
                    "tool_version": __version__,
                    "agglomerate": spec.name,
                    "replica": replica,
                    "packmol_seed": seed,
                    "center_scale": scale,
                    "adaptive_compaction_scale": compact_scale,
                    "xtb": (
                        {
                            "enabled": True,
                            "md_enabled": xtb_md_enabled,
                            "protocol": xtb_protocol,
                            "working_directory": str(xtb_relative),
                            "vasp_source_xyz": optimized.name,
                            "optimized_xyz_sha256": _sha256(optimized),
                        }
                        if xtb_enabled
                        else {"enabled": False}
                    ),
                    "vasp_md_schedule": schedules,
                    "cell_angstrom": [cell, cell, cell],
                    "minimum_intermolecular_distance_angstrom": nearest,
                    "sanity_status": sanity["status"],
                    "sanity_report": "sanity_report.json",
                    "composition": composition,
                    "reference_files": copied,
                    "vasp_ready": not structures_only,
                    "atom_map": atom_map,
                }
                (case / "agglomeration_manifest.json").write_text(
                    json.dumps(case_manifest, indent=2) + "\n", encoding="utf-8"
                )
                relative = case.relative_to(output)
                common = {
                    "agglomerate": spec.name,
                    "case": case_name,
                    "path": str(relative),
                    "replica": replica,
                    "packmol_seed": seed,
                    "center_scale": scale,
                    "atoms": len(symbols),
                    "molecules": len(instances),
                    "cell_angstrom": cell,
                    "minimum_intermolecular_distance_angstrom": nearest,
                }
                index_rows.append(
                    common
                    | {
                        "minimum_face_clearance_angstrom": sanity[
                            "minimum_face_clearance_angstrom"
                        ],
                        "periodic_image_separation_lower_bound_angstrom": sanity[
                            "periodic_image_separation_lower_bound_angstrom"
                        ],
                        "maximum_intramolecular_distance_deviation_angstrom": sanity[
                            "maximum_intramolecular_distance_deviation_angstrom"
                        ],
                        "minimum_intermolecular_oxygen_oxygen_distance_angstrom": sanity[
                            "minimum_intermolecular_oxygen_oxygen_distance_angstrom"
                        ],
                        "maximum_nearest_molecule_distance_angstrom": sanity[
                            "maximum_nearest_molecule_distance_angstrom"
                        ],
                        "potcar_status": sanity["vasp_reference"]["potcar_status"],
                        "sanity_status": sanity["status"],
                    }
                )
                sanity_rows.append(
                    {
                        "agglomerate": spec.name,
                        "case": case_name,
                        "path": str(relative),
                        "status": sanity["status"],
                        "minimum_face_clearance_angstrom": sanity[
                            "minimum_face_clearance_angstrom"
                        ],
                        "requested_vacuum_angstrom": spec.vacuum,
                        "periodic_image_separation_lower_bound_angstrom": sanity[
                            "periodic_image_separation_lower_bound_angstrom"
                        ],
                        "minimum_intermolecular_distance_angstrom": nearest,
                        "required_minimum_intermolecular_distance_angstrom": (
                            spec.minimum_distance
                        ),
                        "maximum_intramolecular_distance_deviation_angstrom": sanity[
                            "maximum_intramolecular_distance_deviation_angstrom"
                        ],
                        "minimum_intermolecular_oxygen_oxygen_distance_angstrom": sanity[
                            "minimum_intermolecular_oxygen_oxygen_distance_angstrom"
                        ],
                        "maximum_nearest_molecule_distance_angstrom": sanity[
                            "maximum_nearest_molecule_distance_angstrom"
                        ],
                        "potcar_status": sanity["vasp_reference"]["potcar_status"],
                    }
                )
                structures.append(
                    {
                        "case": case_name,
                        "path": str(relative),
                        "center_scale": scale,
                        "minimum_intermolecular_distance_angstrom": nearest,
                        "sanity_status": sanity["status"],
                    }
                )
            replica_records.append(
                {
                    "agglomerate": spec.name,
                    "replica": replica,
                    "seed": seed,
                    "status": (
                        "XTB_REQUIRED"
                        if any(item.get("status") == "XTB_REQUIRED" for item in structures)
                        else "COMPLETE"
                    ),
                    "packed_xyz_sha256": _sha256(packed),
                    "adaptive_compaction_scale": compact_scale,
                    "structures": structures,
                }
            )

    if xtb_enabled and xtb_cases:
        _write_xtb_batch(output, xtb_cases, xtb_settings)

    if index_rows:
        with (output / "agglomeration_index.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(index_rows[0]))
            writer.writeheader()
            writer.writerows(index_rows)
    if sanity_rows:
        with (output / "agglomeration_sanity.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(sanity_rows[0]))
            writer.writeheader()
            writer.writerows(sanity_rows)
    root_manifest = {
        "schema_version": 4,
        "generator": "nio-md-prep prepare-agglomeration",
        "tool_version": __version__,
        "status": status,
        "sanity_status": (
            "PASS"
            if status == "COMPLETE"
            and sanity_rows
            and all(row["status"] == "PASS" for row in sanity_rows)
            else "PENDING"
        ),
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "settings": (
            {
                "replicas": specs[0].replicas,
                "base_seed": specs[0].base_seed,
                "radius_angstrom": specs[0].radius,
                "packmol_tolerance_angstrom": specs[0].tolerance,
                "vacuum_angstrom": specs[0].vacuum,
                "minimum_distance_angstrom": specs[0].minimum_distance,
                "compact_to_distance_angstrom": specs[0].compact_to_distance,
                "center_scales": specs[0].scales,
            }
            if not grouped
            else config.get("agglomeration", {})
        ),
        "vasp_reference_sanity": reference_validation,
        "xtb": xtb_settings
        | {
            "enabled": xtb_enabled,
            "md_enabled": xtb_md_enabled,
            "head_bias_enabled": head_bias_enabled,
            "head_bias_replicas_per_family": head_bias_replicas,
            "md_temperature_K": float(
                xtb_settings.get("md_temperature_K", 400.0)
            ),
            "md_time_ps": float(xtb_settings.get("md_time_ps", 10.0)),
            "md_dump_fs": float(xtb_settings.get("md_dump_fs", 50.0)),
            "md_step_fs": float(xtb_settings.get("md_step_fs", 1.0)),
            "md_sccacc": float(xtb_settings.get("md_sccacc", 1.0)),
            "head_bias_time_ps": float(
                xtb_settings.get("head_bias_time_ps", 4.0)
            ),
            "head_bias_target_pp_distance_angstrom": float(
                xtb_settings.get("head_bias_target_pp_distance_angstrom", 4.5)
            ),
            "head_bias_force_constant": float(
                xtb_settings.get("head_bias_force_constant", 0.005)
            ),
            "head_bias_max_distance_reduction_angstrom": float(
                xtb_settings.get("head_bias_max_distance_reduction_angstrom", 1.5)
            ),
            "head_bias_target_axis_angle_degrees": float(
                xtb_settings.get("head_bias_target_axis_angle_degrees", 170.0)
            ),
        },
        "vasp_md": {
            "hold_steps": hold_steps,
            "heating_steps": heating_steps,
            "schedules": ["300K", "300K_to_400K", "400K"],
        },
        "output_mode": "structures_only" if structures_only else "vasp_runs",
        "case_root": case_root_name,
        "templates": [
            {
                "slug": template.slug,
                "atoms_per_molecule": len(template.symbols),
                "source": str(template.source.relative_to(ROOT)),
                "source_sha256": _sha256(template.source),
                "ligpargen_net_charge": template.net_charge,
            }
            for template in unique_templates.values()
        ],
        "molecules": (
            [
                {
                    "slug": template.slug,
                    "count": template.count,
                    "atoms_per_molecule": len(template.symbols),
                    "source": str(template.source.relative_to(ROOT)),
                    "source_sha256": _sha256(template.source),
                    "ligpargen_net_charge": template.net_charge,
                }
                for template in specs[0].templates
            ]
            if not grouped
            else []
        ),
        "agglomerates": [
            {
                "name": spec.name,
                "replicas": spec.replicas,
                "base_seed": spec.base_seed,
                "radius_angstrom": spec.radius,
                "packmol_tolerance_angstrom": spec.tolerance,
                "vacuum_angstrom": spec.vacuum,
                "minimum_distance_angstrom": spec.minimum_distance,
                "compact_to_distance_angstrom": spec.compact_to_distance,
                "center_scales": spec.scales,
                "composition": [
                    {"slug": template.slug, "count": template.count}
                    for template in spec.templates
                ],
            }
            for spec in specs
        ],
        "replicas": replica_records,
    }
    manifest_path.write_text(json.dumps(root_manifest, indent=2) + "\n", encoding="utf-8")
    return manifest_path
