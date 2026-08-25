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
from collections import Counter
from dataclasses import dataclass, replace
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
_REQUIRED_VASP_INPUTS = ("POSCAR", "INCAR", "KPOINTS", "POTCAR")
_VASP_RUNTIME_OUTPUTS = _VASP_OUTPUTS - {"POSCAR"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _existing_vasp_runtime_outputs(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file()
        and (
            path.name in _VASP_RUNTIME_OUTPUTS
            or path.name == "vasp.out"
            or (path.name.startswith("slurm-") and path.suffix == ".out")
        )
    )


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


def _read_poscar(
    path: Path,
) -> tuple[list[str], list[tuple[float, float, float]], dict]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 8:
        raise ValueError(f"{path}: incomplete POSCAR")
    scale_fields = lines[1].split()
    if len(scale_fields) != 1:
        raise ValueError(f"{path}: only a single positive POSCAR scale is supported")
    scale = float(scale_fields[0])
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError(f"{path}: POSCAR scale must be finite and positive")
    lattice = [
        tuple(scale * float(value) for value in lines[index].split()[:3])
        for index in range(2, 5)
    ]
    if any(len(lines[index].split()) < 3 for index in range(2, 5)):
        raise ValueError(f"{path}: malformed POSCAR lattice vector")
    element_names = lines[5].split()
    if not element_names or any(
        not re.fullmatch(r"[A-Z][a-z]?", item) for item in element_names
    ):
        raise ValueError(f"{path}: a VASP 5 element-symbol line is required")
    try:
        counts = [int(value) for value in lines[6].split()]
    except ValueError as error:
        raise ValueError(f"{path}: malformed POSCAR element counts") from error
    if len(counts) != len(element_names) or any(value < 0 for value in counts):
        raise ValueError(f"{path}: POSCAR element names/counts do not match")
    coordinate_line = 7
    if lines[coordinate_line].strip().lower().startswith("s"):
        coordinate_line += 1
    if coordinate_line >= len(lines):
        raise ValueError(f"{path}: missing POSCAR coordinate mode")
    mode = lines[coordinate_line].strip().lower()
    direct = mode.startswith("d")
    cartesian = mode.startswith(("c", "k"))
    if not direct and not cartesian:
        raise ValueError(f"{path}: coordinate mode must be Direct or Cartesian")
    atom_count = sum(counts)
    rows = lines[coordinate_line + 1:coordinate_line + 1 + atom_count]
    if len(rows) != atom_count:
        raise ValueError(f"{path}: POSCAR coordinate count mismatch")
    raw = []
    for offset, row in enumerate(rows, coordinate_line + 2):
        fields = row.split()
        if len(fields) < 3:
            raise ValueError(f"{path}:{offset}: malformed POSCAR coordinate")
        point = tuple(float(value) for value in fields[:3])
        if not all(math.isfinite(value) for value in point):
            raise ValueError(f"{path}:{offset}: non-finite POSCAR coordinate")
        if direct:
            point = tuple(
                sum(point[component] * lattice[component][axis] for component in range(3))
                for axis in range(3)
            )
        else:
            point = tuple(scale * value for value in point)
        raw.append(point)
    symbols = [
        element
        for element, count in zip(element_names, counts)
        for _ in range(count)
    ]
    return symbols, raw, {
        "title": lines[0],
        "element_order": element_names,
        "counts": counts,
        "coordinate_mode": "Direct" if direct else "Cartesian",
        "lattice_angstrom": lattice,
    }


_COVALENT_RADII_ANGSTROM = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "P": 1.07,
}


def _infer_molecular_adjacency(
    symbols: list[str],
    coordinates: list[tuple[float, float, float]],
    tolerance: float,
) -> list[set[int]]:
    adjacency = [set() for _ in symbols]
    for left, left_symbol in enumerate(symbols):
        if left_symbol == "H":
            candidates = []
            for right, right_symbol in enumerate(symbols):
                if right_symbol == "H":
                    continue
                cutoff = (
                    _COVALENT_RADII_ANGSTROM[left_symbol]
                    + _COVALENT_RADII_ANGSTROM[right_symbol]
                    + tolerance
                )
                distance = _distance(coordinates[left], coordinates[right])
                if 0.55 <= distance <= cutoff:
                    candidates.append((distance, right))
            if candidates:
                _, right = min(candidates)
                adjacency[left].add(right)
                adjacency[right].add(left)
            continue
        for right in range(left + 1, len(symbols)):
            right_symbol = symbols[right]
            if right_symbol == "H":
                continue
            try:
                cutoff = (
                    _COVALENT_RADII_ANGSTROM[left_symbol]
                    + _COVALENT_RADII_ANGSTROM[right_symbol]
                    + tolerance
                )
            except KeyError as error:
                raise ValueError(
                    f"no covalent radius is available for POSCAR element {error.args[0]}"
                ) from None
            distance = _distance(coordinates[left], coordinates[right])
            if 0.55 <= distance <= cutoff:
                adjacency[left].add(right)
                adjacency[right].add(left)
    return adjacency


def _labeled_graph_mapping(
    template_symbols: tuple[str, ...],
    template_adjacency: list[set[int]],
    target_symbols: list[str],
    target_adjacency: list[set[int]],
) -> list[int] | None:
    def signature(index: int, symbols, adjacency) -> tuple:
        return (
            symbols[index],
            len(adjacency[index]),
            tuple(sorted(Counter(symbols[item] for item in adjacency[index]).items())),
        )

    target_by_signature: dict[tuple, list[int]] = {}
    for index in range(len(target_symbols)):
        target_by_signature.setdefault(
            signature(index, target_symbols, target_adjacency), []
        ).append(index)
    candidates = {
        index: target_by_signature.get(
            signature(index, template_symbols, template_adjacency), []
        )
        for index in range(len(template_symbols))
    }
    if any(not items for items in candidates.values()):
        return None
    order = sorted(
        range(len(template_symbols)),
        key=lambda index: (len(candidates[index]), -len(template_adjacency[index])),
    )
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def search(position: int) -> bool:
        if position == len(order):
            return True
        template_index = order[position]
        for target_index in candidates[template_index]:
            if target_index in used:
                continue
            if any(
                ((other in template_adjacency[template_index])
                 != (mapping[other] in target_adjacency[target_index]))
                for other in mapping
            ):
                continue
            mapping[template_index] = target_index
            used.add(target_index)
            if search(position + 1):
                return True
            used.remove(target_index)
            del mapping[template_index]
        return False

    if not search(0):
        return None
    return [mapping[index] for index in range(len(template_symbols))]


def _restore_topology_order(
    template: MoleculeTemplate,
    poscar_symbols: list[str],
    poscar_coordinates: list[tuple[float, float, float]],
) -> tuple[list[tuple[float, float, float]], float]:
    data = parse(template.source)
    atom_position = {
        int(atom.fields[0]): index for index, atom in enumerate(data.sections["Atoms"])
    }
    template_adjacency = [set() for _ in template.symbols]
    for bond in data.sections.get("Bonds", []):
        left_id, right_id = map(int, bond.fields[2:4])
        left, right = atom_position[left_id], atom_position[right_id]
        template_adjacency[left].add(right)
        template_adjacency[right].add(left)
    for tolerance in (0.20, 0.30, 0.40, 0.50):
        target_adjacency = _infer_molecular_adjacency(
            poscar_symbols, poscar_coordinates, tolerance
        )
        mapping = _labeled_graph_mapping(
            template.symbols,
            template_adjacency,
            poscar_symbols,
            target_adjacency,
        )
        if mapping is not None:
            return [poscar_coordinates[index] for index in mapping], tolerance
    raise ValueError(
        "could not match the POSCAR bond graph to the LigParGen topology; "
        "check that POSCAR contains one intact molecule"
    )


def _template_from_reference_poscar(
    template: MoleculeTemplate,
    poscar: Path,
) -> tuple[MoleculeTemplate, dict]:
    poscar_symbols, poscar_coordinates, metadata = _read_poscar(poscar)
    expected_counts = {
        symbol: template.symbols.count(symbol) for symbol in sorted(set(template.symbols))
    }
    actual_counts = {
        symbol: poscar_symbols.count(symbol) for symbol in sorted(set(poscar_symbols))
    }
    if actual_counts != expected_counts:
        raise ValueError(
            f"{poscar}: composition {actual_counts} does not match one "
            f"{template.slug} molecule {expected_counts}"
        )
    try:
        reordered, graph_tolerance = _restore_topology_order(
            template, poscar_symbols, poscar_coordinates
        )
    except ValueError as error:
        raise ValueError(f"{poscar}: {error}") from None

    data = parse(template.source)
    atom_position = {
        int(atom.fields[0]): index
        for index, atom in enumerate(data.sections["Atoms"])
    }
    bond_deviations = []
    bond_lengths = []
    for bond in data.sections.get("Bonds", []):
        left_id, right_id = map(int, bond.fields[2:4])
        left, right = atom_position[left_id], atom_position[right_id]
        length = _distance(reordered[left], reordered[right])
        bond_lengths.append(length)
        reference_length = _distance(
            template.coordinates[left], template.coordinates[right]
        )
        bond_deviations.append(abs(length - reference_length))
    maximum_bond_deviation = max(bond_deviations, default=0.0)
    if (
        not bond_lengths
        or min(bond_lengths) < 0.55
        or max(bond_lengths) > 2.30
        or maximum_bond_deviation > 0.75
    ):
        raise ValueError(
            f"{poscar}: topology-based restoration to LigParGen order failed "
            "bond-geometry validation; check that POSCAR contains one intact molecule"
        )
    center = _center_of_mass(reordered, template.masses)
    centered = tuple(
        tuple(point[axis] - center[axis] for axis in range(3))
        for point in reordered
    )
    return replace(template, coordinates=centered), metadata | {
        "path": str(poscar),
        "sha256": _sha256(poscar),
        "slug": template.slug,
        "atom_count": len(centered),
        "mapping": "bond-graph isomorphism restored LigParGen atom order",
        "bond_inference_tolerance_angstrom": graph_tolerance,
        "translation": "center of mass shifted to origin",
        "minimum_topology_bond_angstrom": min(bond_lengths),
        "maximum_topology_bond_angstrom": max(bond_lengths),
        "maximum_bond_length_deviation_from_ligpargen_angstrom": (
            maximum_bond_deviation
        ),
    }


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


def _potcar_blocks(path: Path) -> dict[str, bytes]:
    """Return complete POTCAR datasets keyed by their TITEL element."""
    content = path.read_bytes()
    matches = list(
        re.finditer(
            rb"(?m)^\s*TITEL\s*=\s*\S+\s+([A-Z][a-z]?)",
            content,
        )
    )
    if not matches:
        raise ValueError(f"{path}: no POTCAR TITEL records found")
    blocks: dict[str, bytes] = {}
    for index, match in enumerate(matches):
        start = 0
        if index:
            marker = content.rfind(
                b"End of Dataset", matches[index - 1].end(), match.start()
            )
            newline = content.find(b"\n", marker) if marker >= 0 else -1
            start = newline + 1 if newline >= 0 else match.start()
        next_title = (
            matches[index + 1].start() if index + 1 < len(matches) else len(content)
        )
        marker = content.find(b"End of Dataset", match.end(), next_title)
        if marker >= 0:
            newline = content.find(b"\n", marker)
            stop = len(content) if newline < 0 else newline + 1
        else:
            stop = next_title
        element = match.group(1).decode("ascii")
        if element in blocks:
            raise ValueError(f"{path}: duplicate POTCAR dataset for {element}")
        blocks[element] = content[start:stop]
    return blocks


def _prepare_mixed_vasp_reference(
    references: dict[str, Path],
    template_slug: str,
    expected_elements: list[str],
    destination: Path,
) -> tuple[Path, dict]:
    """Build one calculation template while retaining per-species geometries."""
    if template_slug not in references:
        raise ValueError(
            f"VASP template slug {template_slug!r} has no matching --reference-dir"
        )
    for slug, reference in references.items():
        if not reference.is_dir():
            raise FileNotFoundError(f"VASP reference directory not found: {reference}")
        missing = [
            name for name in _REQUIRED_VASP_INPUTS if not (reference / name).is_file()
        ]
        if missing:
            raise ValueError(
                f"{reference} ({slug}): incomplete VASP reference directory; missing "
                + ", ".join(missing)
            )

    blocks_by_element: dict[str, tuple[bytes, str, Path]] = {}
    potcar_sources: dict[str, dict] = {}
    for slug, reference in references.items():
        potcar = reference / "POTCAR"
        potcar_sources[slug] = {
            "path": str(potcar),
            "sha256": _sha256(potcar),
            "elements": _potcar_elements(potcar),
        }
        for element, block in _potcar_blocks(potcar).items():
            previous = blocks_by_element.get(element)
            if previous is not None and previous[0] != block:
                raise ValueError(
                    f"conflicting POTCAR datasets for {element}: "
                    f"{previous[2]} ({previous[1]}) and {potcar} ({slug})"
                )
            if previous is None:
                blocks_by_element[element] = (block, slug, potcar)
    missing_elements = [
        element for element in expected_elements if element not in blocks_by_element
    ]
    if missing_elements:
        raise ValueError(
            "mixed VASP references do not provide POTCAR datasets for "
            + ", ".join(missing_elements)
        )

    destination.mkdir(parents=True, exist_ok=True)
    for stale in destination.iterdir():
        if stale.is_file():
            stale.unlink()
    primary = references[template_slug]
    copied: list[dict] = []
    for source in sorted(primary.iterdir()):
        if (
            not source.is_file()
            or source.name in _VASP_OUTPUTS
            or source.name == "POTCAR"
        ):
            continue
        target = destination / source.name
        shutil.copy2(source, target)
        copied.append({"name": source.name, "sha256": _sha256(target)})
    shutil.copy2(primary / "POSCAR", destination / "POSCAR")
    combined = b"".join(
        blocks_by_element[element][0] for element in expected_elements
    )
    if combined and not combined.endswith(b"\n"):
        combined += b"\n"
    (destination / "POTCAR").write_bytes(combined)
    return destination, {
        "mode": "mixed_species",
        "calculation_template_slug": template_slug,
        "calculation_template_directory": str(primary),
        "shared_files": copied,
        "potcar_sources": potcar_sources,
        "combined_potcar_elements": expected_elements,
        "combined_potcar_sha256": _sha256(destination / "POTCAR"),
        "combined_potcar_element_sources": {
            element: blocks_by_element[element][1] for element in expected_elements
        },
    }


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
if [[ "${1:-}" == "--dry-run" ]]; then
  printf '(cd %q && sbatch --parsable runvasp.sh)\n' "$root/300K"
  printf '(cd %q && sbatch --parsable runvasp.sh)\n' "$root/400K/01_heat_300_to_400K"
  printf '(cd %q && sbatch --parsable --dependency=afterok:<HEAT_JOB_ID> runvasp.sh)\n' "$root/400K/02_hold_400K"
  exit 0
elif [[ $# -ne 0 ]]; then
  echo "usage: $0 [--dry-run]" >&2
  exit 2
fi
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


def _write_vasp_campaign_launcher(output: Path, case_root_name: str) -> Path:
    launcher = output / "launch_vasp_runs.sh"
    launcher.write_text(
        f'''#!/bin/bash
set -euo pipefail
dry_run=false
assume_yes=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) dry_run=true ;;
    --yes|-y) assume_yes=true ;;
    -h|--help)
      echo "usage: $0 [--dry-run] [--yes]"
      exit 0
      ;;
    *)
      echo "unknown option: $1" >&2
      echo "usage: $0 [--dry-run] [--yes]" >&2
      exit 2
      ;;
  esac
  shift
done

root=$(cd "$(dirname "$0")" && pwd)
case_root="$root/{case_root_name}"
cases=()
while IFS= read -r -d '' submit; do
  cases+=("$(dirname "$submit")")
done < <(find "$case_root" -type f -name submit_temperature_jobs.sh -print0 | sort -z)

if [[ ${{#cases[@]}} -eq 0 ]]; then
  echo "no submit_temperature_jobs.sh files found under $case_root" >&2
  exit 1
fi

case_count=${{#cases[@]}}
stage_count=$((case_count * 3))
printf 'Found %d structure case(s) and %d VASP stage job(s) under %s\n' \
  "$case_count" "$stage_count" "$case_root"

if [[ "$dry_run" == true ]]; then
  echo "DRY RUN: no jobs will be submitted"
  for case_dir in "${{cases[@]}}"; do
    bash "$case_dir/submit_temperature_jobs.sh" --dry-run
  done
  exit 0
fi

if [[ "$assume_yes" != true ]]; then
  reply=""
  if ! read -r -p "Submit all $stage_count VASP stage jobs? [y/N] " reply; then
    reply=""
  fi
  case "$reply" in
    y|Y|yes|YES|Yes) ;;
    *)
      echo "Submission cancelled; no jobs were launched."
      exit 0
      ;;
  esac
fi

for case_dir in "${{cases[@]}}"; do
  printf 'Submitting %s\n' "$case_dir"
  bash "$case_dir/submit_temperature_jobs.sh"
done
echo "All VASP stage jobs were submitted."
''',
        encoding="utf-8",
    )
    launcher.chmod(0o755)
    return launcher


def _write_xtb_batch(output: Path, cases: list[dict], settings: dict) -> None:
    case_list = output / "xtb_cases.tsv"
    case_list.write_text(
        "".join(
            f"{item['path']}\t{item['charge']}\t{item['uhf']}\n" for item in cases
        ),
        encoding="utf-8",
    )
    account = str(settings.get("account", "loni_perovsk27"))
    partition = str(settings.get("partition", "workq"))
    cpus = int(settings.get("cpus", 8))
    walltime = str(settings.get("walltime", "24:00:00"))
    gfn = int(settings.get("gfn", 2))
    opt_level = str(settings.get("opt_level", "loose"))
    max_cycles = int(settings.get("max_cycles", 100))
    concurrency = int(settings.get("array_concurrency", 4))
    omp_stacksize = str(settings.get("omp_stacksize", "1G")).upper()
    if cpus < 1:
        raise ValueError("xtb.cpus must be positive")
    if concurrency < 1:
        raise ValueError("xtb.array_concurrency must be positive")
    if not re.fullmatch(r"[1-9][0-9]*[KMGT]", omp_stacksize):
        raise ValueError(
            "xtb.omp_stacksize must be a positive size such as '512M', '1G', or '4G'"
        )
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
        workflow = f'''input_changed=false
if [[ -s xtb_input.complete ]]; then
  if ! cmp -s xtb_input.sha256 xtb_input.complete; then
    input_changed=true
  fi
elif [[ -s xtbopt_initial.xyz || -s xtbsteered_last.xyz || -s xtbmd_unbiased_last.xyz || -s xtbfinal.xyz ]]; then
  input_changed=true
fi
if [[ "$input_changed" == true ]]; then
  old_tag=legacy
  if [[ -s xtb_input.complete ]]; then
    old_tag=$(head -c 12 xtb_input.complete)
  fi
  for artifact in xtbopt.xyz xtbopt.log xtbopt_initial.xyz xtbopt_initial.log xtb_initial_opt.out xtb_initial_opt.err md_steered.inp steering_restraints.json xtb_steered.trj xtbsteered_last.xyz xtb_steered_md.out xtb_steered_md.err xtb_unbiased.trj xtbmd_unbiased_last.xyz xtb_unbiased_md.out xtb_unbiased_md.err xtbfinal.xyz xtbfinal.log xtb_final_opt.out xtb_final_opt.err; do
    if [[ -e "$artifact" ]]; then
      archived="${{artifact}}.input_${{old_tag}}"
      if [[ -e "$archived" ]]; then
        rm -f "$artifact"
      else
        mv "$artifact" "$archived"
      fi
    fi
  done
  rm -f xtb_protocol.complete xtb_stages.protocol
fi
cp xtb_input.sha256 xtb_input.complete
if [[ -s xtbopt.xyz && ! -s xtbopt_initial.xyz ]]; then
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
#SBATCH --cpus-per-task={cpus}
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
threads=${{SLURM_CPUS_PER_TASK:-{cpus}}}
export OMP_NUM_THREADS="${{threads}},1"
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=dynamic
export OMP_STACKSIZE=${{XTB_OMP_STACKSIZE:-{omp_stacksize}}}
export MKL_NUM_THREADS="$threads"
export MKL_DYNAMIC=FALSE
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

    audit = output / "audit_xtb_runs.py"
    audit.write_text(
        f'''#!/usr/bin/env python3
"""Audit resumable xTB cases without changing calculation outputs."""

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MD_ENABLED = {md_enabled!r}
CPUS_PER_CASE = {cpus}
ARRAY_CONCURRENCY = {concurrency}
COMPLETE_STATES = {{"COMPLETE"}}


def same_hash(left: Path, right: Path) -> bool:
    return (
        left.is_file()
        and right.is_file()
        and left.read_bytes() == right.read_bytes()
    )


def classify(work: Path) -> tuple[str, str]:
    if not work.is_dir():
        return "MISSING_DIRECTORY", "case directory is absent"
    if not MD_ENABLED:
        if (work / "xtbopt.xyz").is_file():
            return "COMPLETE", "xtbopt.xyz present"
        if (work / "xtb.err").is_file() or (work / "xtb.out").is_file():
            return "OPTIMIZATION_STARTED", "xTB output exists; no xtbopt.xyz"
        return "PENDING", "optimization has not produced xtbopt.xyz"

    final = work / "xtbfinal.xyz"
    expected = work / "xtb_protocol.sha256"
    completed = work / "xtb_protocol.complete"
    if final.is_file() and same_hash(expected, completed):
        return "COMPLETE", "final geometry and current completion marker present"
    if final.is_file() and completed.is_file():
        return (
            "HASH_MISMATCH",
            "final geometry present but completion hash differs; rerun preparation "
            "with --regenerate-vasp to validate/migrate a legacy result",
        )
    if final.is_file():
        return "FINAL_UNMARKED", "xtbfinal.xyz present without completion marker"
    if completed.is_file():
        return "MARKER_WITHOUT_FINAL", "completion marker present without xtbfinal.xyz"
    if (work / "xtbmd_unbiased_last.xyz").is_file():
        return "UNBIASED_MD_COMPLETE", "final optimization remains"
    if (work / "xtbsteered_last.xyz").is_file():
        return "STEERING_COMPLETE", "unbiased MD and final optimization remain"
    if (work / "xtbopt_initial.xyz").is_file():
        return "INITIAL_OPT_COMPLETE", "MD and final optimization remain"
    if (work / "xtb_initial_opt.out").is_file():
        return "INITIAL_OPT_STARTED", "initial optimization output exists"
    return "PENDING", "initial optimization has not completed"


def rows() -> list[dict[str, object]]:
    result = []
    lines = (ROOT / "xtb_cases.tsv").read_text(encoding="utf-8").splitlines()
    for array_index, raw in enumerate(lines):
        if not raw.strip():
            continue
        relative, charge, uhf = raw.split("\\t")
        status, detail = classify(ROOT / relative)
        result.append(
            {{
                "array_index": array_index,
                "case": relative,
                "status": status,
                "charge": charge,
                "uhf": uhf,
                "detail": detail,
            }}
        )
    return result


def compact_indices(indices: list[int]) -> str:
    if not indices:
        return ""
    ranges = []
    start = previous = indices[0]
    for index in indices[1:]:
        if index == previous + 1:
            previous = index
            continue
        ranges.append(str(start) if start == previous else f"{{start}}-{{previous}}")
        start = previous = index
    ranges.append(str(start) if start == previous else f"{{start}}-{{previous}}")
    return ",".join(ranges)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--count-pending",
        action="store_true",
        help="print only the number of cases that are not safely complete",
    )
    parser.add_argument(
        "--pending-indices",
        action="store_true",
        help="print only the compact Slurm indices requiring submission",
    )
    parser.add_argument(
        "--summary-only", action="store_true", help="omit per-case status lines"
    )
    args = parser.parse_args()
    audited = rows()
    pending = sum(row["status"] not in COMPLETE_STATES for row in audited)
    if args.count_pending:
        print(pending)
        return 0
    if args.pending_indices:
        print(
            compact_indices(
                [
                    int(row["array_index"])
                    for row in audited
                    if row["status"] not in COMPLETE_STATES
                ]
            )
        )
        return 0

    counts = Counter(row["status"] for row in audited)
    completed = len(audited) - pending
    timestamp = datetime.now(timezone.utc).isoformat()
    csv_path = ROOT / "xtb_progress.csv"
    json_path = ROOT / "xtb_progress.json"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "array_index",
                "case",
                "status",
                "charge",
                "uhf",
                "detail",
            ],
        )
        writer.writeheader()
        writer.writerows(audited)
    json_path.write_text(
        json.dumps(
            {{
                "audited_at_utc": timestamp,
                "total": len(audited),
                "completed": completed,
                "pending": pending,
                "counts_by_status": dict(sorted(counts.items())),
                "cases": audited,
            }},
            indent=2,
        )
        + "\\n",
        encoding="utf-8",
    )

    percent = 100.0 if not audited else 100.0 * completed / len(audited)
    print(
        f"xTB progress: {{completed}}/{{len(audited)}} safely complete "
        f"({{percent:.1f}}%); {{pending}} pending or requiring attention."
    )
    print(
        f"Resources: {{CPUS_PER_CASE}} OpenMP core(s) per case; "
        f"array concurrency {{ARRAY_CONCURRENCY}}; at most "
        f"{{CPUS_PER_CASE * ARRAY_CONCURRENCY}} allocated cores."
    )
    print(
        "Stages: "
        + (", ".join(f"{{name}}={{count}}" for name, count in sorted(counts.items())) or "none")
    )
    if not args.summary_only:
        for row in audited:
            print(f"[{{row['status']}}] {{row['case']}} - {{row['detail']}}")
    print(f"Reports written: {{csv_path}} and {{json_path}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
''',
        encoding="utf-8",
    )
    audit.chmod(0o755)

    launcher = output / "launch_xtb.sh"
    launcher.write_text(
        """#!/bin/bash
set -euo pipefail

root=$(cd "$(dirname "$0")" && pwd)
audit="$root/audit_xtb_runs.py"
array="$root/run_xtb_array.sbatch"
dry_run=false
assume_yes=false

usage() {
  cat <<'EOF'
Usage: launch_xtb.sh [--dry-run] [--yes]

Audit xTB progress and submit the resumable Slurm array. The array itself
checks every case and skips cases carrying the current completion marker.
EOF
}

while (($#)); do
  case "$1" in
    --dry-run) dry_run=true ;;
    --yes|-y) assume_yes=true ;;
    --help|-h) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

"$audit"
pending=$("$audit" --count-pending)
pending_indices=$("$audit" --pending-indices)
if ((pending == 0)); then
  echo "All xTB cases are safely complete; no array was submitted."
  exit 0
fi

echo "$pending xTB case(s) remain pending or require attention."
array_selection="${pending_indices}%CONCURRENCY_PLACEHOLDER"
echo "Selected Slurm array indices: $array_selection"
if [[ "$dry_run" == true ]]; then
  echo "DRY RUN: no jobs will be submitted."
  printf 'Would run: cd %q && sbatch --array=%q %q\n' \
    "$root" "$array_selection" "$(basename "$array")"
  exit 0
fi

if [[ "$assume_yes" != true ]]; then
  read -r -p "Submit the resumable xTB array now? [y/N] " reply
  case "$reply" in
    y|Y|yes|YES|Yes) ;;
    *) echo "Submission cancelled; no jobs were launched."; exit 0 ;;
  esac
fi

(cd "$root" && sbatch --array="$array_selection" "$(basename "$array")")
echo "xTB array submitted. Re-run $audit at any time to check progress."
""".replace("CONCURRENCY_PLACEHOLDER", str(concurrency)),
        encoding="utf-8",
    )
    launcher.chmod(0o755)


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
    temperature = float(settings.get("md_temperature_K", 450.0))
    total_time_ps = float(settings.get("md_time_ps", 20.0))
    steering_time_ps = (
        float(settings.get("head_bias_time_ps", 6.0)) if steered else 0.0
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
        "schema_version": 4,
        "input_xyz_sha256": _sha256(directory / "input.xyz"),
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
    (directory / "xtb_input.sha256").write_text(
        protocol["input_xyz_sha256"] + "\n", encoding="utf-8"
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
    reference_dirs: dict[str, Path] | None = None,
    vasp_template_slug: str | None = None,
    packed_xyz: Path | None = None,
    structures_only: bool = False,
    regenerate_vasp: bool = False,
) -> Path:
    """Prepare size-stratified molecular clusters as VASP run packages."""
    if reference_dir is not None and reference_dirs is not None:
        raise ValueError("reference_dir and reference_dirs are mutually exclusive")
    if vasp_template_slug is not None and reference_dirs is None:
        raise ValueError(
            "--vasp-template-slug requires qualified --reference-dir values"
        )
    if (reference_dir is not None or reference_dirs is not None) and structures_only:
        raise ValueError("--reference-dir and --structures-only are mutually exclusive")
    if regenerate_vasp and structures_only:
        raise ValueError("--regenerate-vasp cannot be used with --structures-only")
    if reference_dir is not None:
        reference_dir = Path(reference_dir)
    if reference_dirs is not None:
        if not reference_dirs:
            raise ValueError("reference_dirs cannot be empty")
        reference_dirs = {
            str(slug): Path(directory) for slug, directory in reference_dirs.items()
        }
    manifest_path = output / "agglomeration_manifest.json"
    resume = False
    if output.exists() and any(output.iterdir()):
        if manifest_path.is_file():
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
            resumable_statuses = {"PACKMOL_REQUIRED", "XTB_REQUIRED"}
            if regenerate_vasp:
                resumable_statuses.add("COMPLETE")
            resume = (
                previous.get("status") in resumable_statuses
                and (
                    regenerate_vasp
                    or previous.get("config_sha256") == _sha256(config_path)
                )
            )
        if not resume:
            raise FileExistsError(f"refusing to replace non-empty output directory: {output}")
    if regenerate_vasp:
        runtime_outputs = _existing_vasp_runtime_outputs(output / "vasp_runs")
        if runtime_outputs:
            preview = ", ".join(
                str(path.relative_to(output)) for path in runtime_outputs[:5]
            )
            if len(runtime_outputs) > 5:
                preview += f", and {len(runtime_outputs) - 5} more"
            raise ValueError(
                "--regenerate-vasp refuses to modify a campaign containing VASP "
                f"runtime outputs: {preview}"
            )
    output.mkdir(parents=True, exist_ok=True)
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    specs, grouped = _agglomerate_specs(config)
    geometry_metadata: dict
    if reference_dir is not None:
        slugs = sorted(
            {template.slug for spec in specs for template in spec.templates}
        )
        if len(slugs) != 1:
            raise ValueError(
                "a single --reference-dir/POSCAR can define only one molecular "
                "species; this campaign contains " + ", ".join(slugs)
            )
        poscar = reference_dir / "POSCAR"
        if not poscar.is_file():
            raise ValueError(
                f"{reference_dir}: incomplete VASP reference directory; missing POSCAR"
            )
        representative = next(
            template
            for spec in specs
            for template in spec.templates
            if template.slug == slugs[0]
        )
        reference_template, geometry_metadata = _template_from_reference_poscar(
            representative, poscar
        )
        specs = [
            replace(
                spec,
                templates=tuple(
                    replace(template, coordinates=reference_template.coordinates)
                    for template in spec.templates
                ),
            )
            for spec in specs
        ]
        geometry_metadata = {"source": "reference_poscar"} | geometry_metadata
    elif reference_dirs is not None:
        required_slugs = {
            template.slug for spec in specs for template in spec.templates
        }
        supplied_slugs = set(reference_dirs)
        if supplied_slugs != required_slugs:
            missing = sorted(required_slugs - supplied_slugs)
            extra = sorted(supplied_slugs - required_slugs)
            details = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if extra:
                details.append("unknown " + ", ".join(extra))
            raise ValueError(
                "qualified VASP references must match campaign molecule slugs: "
                + "; ".join(details)
            )
        template_slug = vasp_template_slug or next(iter(reference_dirs))
        reference_templates: dict[str, MoleculeTemplate] = {}
        geometry_by_slug: dict[str, dict] = {}
        for slug, directory in reference_dirs.items():
            representative = next(
                template
                for spec in specs
                for template in spec.templates
                if template.slug == slug
            )
            poscar = directory / "POSCAR"
            if not poscar.is_file():
                raise ValueError(
                    f"{directory}: incomplete VASP reference directory; missing POSCAR"
                )
            reference_template, metadata = _template_from_reference_poscar(
                representative, poscar
            )
            reference_templates[slug] = reference_template
            geometry_by_slug[slug] = metadata
        specs = [
            replace(
                spec,
                templates=tuple(
                    replace(
                        template,
                        coordinates=reference_templates[template.slug].coordinates,
                    )
                    for template in spec.templates
                ),
            )
            for spec in specs
        ]
        geometry_metadata = {
            "source": "reference_poscars",
            "calculation_template_slug": template_slug,
            "by_slug": geometry_by_slug,
        }
    else:
        geometry_metadata = {
            "source": "ligpargen_coordinates",
            "note": "coordinate fallback used only for --structures-only",
        }
    xtb_settings = config.get("xtb", {})
    xtb_enabled = bool(xtb_settings.get("enabled", False))
    if regenerate_vasp and not xtb_enabled:
        raise ValueError("--regenerate-vasp requires an xTB-enabled campaign")
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
    calculation_reference = reference_dir
    mixed_reference_metadata = None
    if reference_dirs is not None:
        template_slug = vasp_template_slug or next(iter(reference_dirs))
        calculation_reference, mixed_reference_metadata = _prepare_mixed_vasp_reference(
            reference_dirs,
            template_slug,
            expected_elements,
            output / "vasp_reference",
        )
    reference_validation = _reference_sanity(
        calculation_reference, expected_elements, structures_only=structures_only
    )
    if mixed_reference_metadata is not None:
        reference_validation |= mixed_reference_metadata
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
            geometry_fingerprint = _sha256_text(
                packmol_input
                + "".join(
                    f"{slug}:{_sha256(path)}\n"
                    for slug, path in sorted(local_templates.items())
                )
            )
            fingerprint_path = work / "packmol_geometry.sha256"
            previous_fingerprint = (
                fingerprint_path.read_text(encoding="utf-8").strip()
                if fingerprint_path.is_file()
                else None
            )
            if packed.exists() and previous_fingerprint != geometry_fingerprint:
                old_tag = (previous_fingerprint or "legacy")[:12]
                archived = work / f"packed.xyz.geometry_{old_tag}"
                suffix = 1
                while archived.exists():
                    archived = work / f"packed.xyz.geometry_{old_tag}_{suffix}"
                    suffix += 1
                shutil.move(packed, archived)
            fingerprint_path.write_text(
                geometry_fingerprint + "\n", encoding="utf-8"
            )
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
                    optimized = xtb_work / (
                        "xtbfinal.xyz" if xtb_md_enabled else "xtbopt.xyz"
                    )
                    completed_protocol = xtb_work / "xtb_protocol.complete"
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
                    reuse_completed_xtb = (
                        regenerate_vasp
                        and optimized.is_file()
                        and (not xtb_md_enabled or completed_protocol.is_file())
                    )
                    if reuse_completed_xtb:
                        optimized_symbols, raw_variant = _read_xyz(optimized)
                        if optimized_symbols != expected_symbols:
                            raise ValueError(
                                f"{optimized}: xTB changed atom ordering or element "
                                "identities"
                            )
                        protocol_path = xtb_work / "xtb_protocol.json"
                        if protocol_path.is_file():
                            try:
                                xtb_protocol = json.loads(
                                    protocol_path.read_text(encoding="utf-8")
                                )
                            except json.JSONDecodeError:
                                xtb_protocol = {"legacy_protocol_readable": False}
                        else:
                            xtb_protocol = {"legacy_protocol_present": False}
                        xtb_protocol["regeneration_reused_completed_xtb"] = True
                        xtb_protocol["legacy_completion_marker_present"] = (
                            completed_protocol.is_file()
                        )
                        expected_protocol = xtb_work / "xtb_protocol.sha256"
                        if (
                            completed_protocol.is_file()
                            and expected_protocol.is_file()
                            and completed_protocol.read_text(encoding="utf-8")
                            != expected_protocol.read_text(encoding="utf-8")
                        ):
                            shutil.copy2(expected_protocol, completed_protocol)
                            xtb_protocol[
                                "legacy_completion_marker_migrated"
                            ] = True
                    else:
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
                        protocol_ready = True
                        if xtb_md_enabled:
                            expected_protocol = xtb_work / "xtb_protocol.sha256"
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
                        calculation_reference,
                        case / "POSCAR",
                        job_prefix=f"{spec.name}_r{replica:03d}",
                        hold_steps=hold_steps,
                        heating_steps=heating_steps,
                    )
                elif not structures_only:
                    copied = _copy_reference(calculation_reference, case)
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
    vasp_launcher = output / "launch_vasp_runs.sh"
    if status == "COMPLETE" and xtb_enabled and not structures_only and index_rows:
        _write_vasp_campaign_launcher(output, case_root_name)
    else:
        vasp_launcher.unlink(missing_ok=True)
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
        "molecular_geometry": geometry_metadata,
        "xtb": xtb_settings
        | {
            "enabled": xtb_enabled,
            "md_enabled": xtb_md_enabled,
            "head_bias_enabled": head_bias_enabled,
            "head_bias_replicas_per_family": head_bias_replicas,
            "md_temperature_K": float(
                xtb_settings.get("md_temperature_K", 450.0)
            ),
            "md_time_ps": float(xtb_settings.get("md_time_ps", 20.0)),
            "md_dump_fs": float(xtb_settings.get("md_dump_fs", 50.0)),
            "md_step_fs": float(xtb_settings.get("md_step_fs", 1.0)),
            "md_sccacc": float(xtb_settings.get("md_sccacc", 1.0)),
            "head_bias_time_ps": float(
                xtb_settings.get("head_bias_time_ps", 6.0)
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
        "regenerated_vasp": regenerate_vasp,
        "vasp_launcher": (
            vasp_launcher.name if vasp_launcher.is_file() else None
        ),
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
