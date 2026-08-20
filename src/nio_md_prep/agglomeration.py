"""Generate reproducible phosphonate agglomeration structures for DFT labeling."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import tomllib
from dataclasses import dataclass
from pathlib import Path

from . import __version__
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
    net_charge: float


@dataclass(frozen=True)
class MoleculeInstance:
    slug: str
    molecule_id: int
    start: int
    stop: int
    masses: tuple[float, ...]


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
        templates.append(
            MoleculeTemplate(
                slug=slug,
                count=count,
                source=source,
                symbols=tuple(elements(data)),
                coordinates=tuple(atom_coordinates(data)),
                masses=_atom_masses(data),
                net_charge=actual_charge,
            )
        )
        total += count
    if total < 2:
        raise ValueError("agglomeration generation requires at least two molecules")
    return templates


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
                )
            )
    return symbols, instances


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


def _minimum_intermolecular_distance(
    coordinates: list[tuple[float, float, float]],
    instances: list[MoleculeInstance],
) -> float:
    molecule_for_atom = [0] * len(coordinates)
    for item in instances:
        molecule_for_atom[item.start:item.stop] = [item.molecule_id] * (item.stop - item.start)
    minimum = math.inf
    for left, first in enumerate(coordinates):
        for right in range(left + 1, len(coordinates)):
            if molecule_for_atom[left] == molecule_for_atom[right]:
                continue
            second = coordinates[right]
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))
            minimum = min(minimum, distance)
    return minimum


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
    nearest = _minimum_intermolecular_distance(coordinates, instances)
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
        "intramolecular_geometry_preserved": rigidity_error <= 1.0e-5,
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
        "required_minimum_intermolecular_distance_angstrom": minimum_distance,
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
    """Prepare compact and separated molecular clusters as VASP run packages."""
    if reference_dir is not None and structures_only:
        raise ValueError("--reference-dir and --structures-only are mutually exclusive")
    manifest_path = output / "agglomeration_manifest.json"
    resume = False
    if output.exists() and any(output.iterdir()):
        if manifest_path.is_file():
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
            resume = (
                previous.get("status") == "PACKMOL_REQUIRED"
                and previous.get("config_sha256") == _sha256(config_path)
            )
        if not resume:
            raise FileExistsError(f"refusing to replace non-empty output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    settings = config.get("agglomeration", {})
    templates = _load_templates(config)
    replicas = int(settings.get("replicas", 4))
    base_seed = int(settings.get("base_seed", 2405367))
    radius = float(settings.get("radius_angstrom", 14.0))
    tolerance = float(settings.get("packmol_tolerance_angstrom", 3.0))
    vacuum = float(settings.get("vacuum_angstrom", 8.0))
    minimum_distance = float(settings.get("minimum_distance_angstrom", 2.5))
    scales = tuple(float(value) for value in settings.get("center_scales", [1.0]))
    if replicas < 1 or not 0 < base_seed < 900000000:
        raise ValueError("replicas must be positive and base_seed must be between 1 and 899999999")
    if radius <= 0 or tolerance <= 0 or vacuum <= 0 or minimum_distance <= 0:
        raise ValueError("radius, tolerance, vacuum, and minimum distance must be positive")
    if not scales or any(scale <= 0 for scale in scales) or len(set(scales)) != len(scales):
        raise ValueError("center_scales must contain distinct positive values")
    if packed_xyz is not None and replicas != 1:
        raise ValueError("--packed-xyz can only be used when replicas = 1")

    template_dir = output / "templates"
    template_dir.mkdir(exist_ok=True)
    template_paths: dict[str, Path] = {}
    for template in templates:
        path = template_dir / f"{template.slug}.xyz"
        _write_template_xyz(template, path)
        template_paths[template.slug] = path
    expected_symbols, instances = _expected_layout(templates)
    expected_elements = sorted(set(expected_symbols))
    reference_validation = _reference_sanity(
        reference_dir, expected_elements, structures_only=structures_only
    )
    case_root_name = "structures" if structures_only else "vasp_runs"
    executable = shutil.which("packmol")
    index_rows: list[dict] = []
    sanity_rows: list[dict] = []
    replica_records: list[dict] = []
    status = "COMPLETE"

    for replica in range(replicas):
        seed = base_seed + replica
        if seed >= 900000000:
            raise ValueError("resolved Packmol seed exceeds 899999999")
        work = output / "packmol" / f"replica_{replica:03d}"
        work.mkdir(parents=True, exist_ok=True)
        local_templates: dict[str, Path] = {}
        for slug, source in template_paths.items():
            destination = work / source.name
            shutil.copy2(source, destination)
            local_templates[slug] = destination
        packed = work / "packed.xyz"
        packmol_input = _packmol_input(
            templates, local_templates, packed.name, seed, tolerance, radius
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
                    f"Packmol failed for replica {replica}: "
                    + run.stderr.decode(errors="replace")
                )
        else:
            status = "PACKMOL_REQUIRED"
            replica_records.append(
                {
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
                f"{packed}: atom ordering differs from the configured LigParGen templates"
            )
        uncentered_variants = {
            scale: _scale_molecule_centers(coordinates, instances, scale)
            for scale in scales
        }
        cell = _fixed_cell(uncentered_variants, vacuum)
        structures: list[dict] = []
        for variant_index, scale in enumerate(scales):
            variant = _center_in_cell(uncentered_variants[scale], cell)
            case_name = f"r{replica:03d}_s{variant_index:02d}_{scale:.3f}".replace(".", "p")
            case = output / case_root_name / case_name
            case.mkdir(parents=True, exist_ok=resume)
            ordered = _ordered_atoms(symbols, variant, instances)
            title = f"phosphonate agglomeration replica={replica} center_scale={scale:g}"
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
            copied = _copy_reference(reference_dir, case)
            sanity = _structure_sanity(
                variant,
                symbols,
                instances,
                templates,
                cell=cell,
                requested_vacuum=vacuum,
                minimum_distance=minimum_distance,
                reference_sanity=reference_validation,
                expected_elements=expected_elements,
                written_element_order=written_element_order,
            )
            _write_sanity_report(case, sanity)
            if sanity["status"] != "PASS":
                failed = [
                    name for name, passed in sanity["checks"].items() if not passed
                ]
                raise ValueError(
                    f"{case}: structural sanity check failed: {', '.join(failed)}; "
                    "see sanity_report.json"
                )
            nearest = sanity["minimum_intermolecular_distance_angstrom"]
            case_manifest = {
                "schema_version": 2,
                "generator": "nio-md-prep prepare-agglomeration",
                "tool_version": __version__,
                "replica": replica,
                "packmol_seed": seed,
                "center_scale": scale,
                "cell_angstrom": [cell, cell, cell],
                "minimum_intermolecular_distance_angstrom": nearest,
                "sanity_status": sanity["status"],
                "sanity_report": "sanity_report.json",
                "composition": [
                    {"slug": template.slug, "count": template.count}
                    for template in templates
                ],
                "reference_files": copied,
                "vasp_ready": not structures_only,
                "atom_map": atom_map,
            }
            (case / "agglomeration_manifest.json").write_text(
                json.dumps(case_manifest, indent=2) + "\n", encoding="utf-8"
            )
            relative = case.relative_to(output)
            index_rows.append(
                {
                    "case": case_name,
                    "path": str(relative),
                    "replica": replica,
                    "packmol_seed": seed,
                    "center_scale": scale,
                    "atoms": len(symbols),
                    "molecules": len(instances),
                    "cell_angstrom": cell,
                    "minimum_intermolecular_distance_angstrom": nearest,
                    "minimum_face_clearance_angstrom": sanity[
                        "minimum_face_clearance_angstrom"
                    ],
                    "periodic_image_separation_lower_bound_angstrom": sanity[
                        "periodic_image_separation_lower_bound_angstrom"
                    ],
                    "maximum_intramolecular_distance_deviation_angstrom": sanity[
                        "maximum_intramolecular_distance_deviation_angstrom"
                    ],
                    "potcar_status": sanity["vasp_reference"]["potcar_status"],
                    "sanity_status": sanity["status"],
                }
            )
            sanity_rows.append(
                {
                    "case": case_name,
                    "path": str(relative),
                    "status": sanity["status"],
                    "minimum_face_clearance_angstrom": sanity[
                        "minimum_face_clearance_angstrom"
                    ],
                    "requested_vacuum_angstrom": vacuum,
                    "periodic_image_separation_lower_bound_angstrom": sanity[
                        "periodic_image_separation_lower_bound_angstrom"
                    ],
                    "minimum_intermolecular_distance_angstrom": nearest,
                    "required_minimum_intermolecular_distance_angstrom": minimum_distance,
                    "maximum_intramolecular_distance_deviation_angstrom": sanity[
                        "maximum_intramolecular_distance_deviation_angstrom"
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
                "replica": replica,
                "seed": seed,
                "status": "COMPLETE",
                "packed_xyz_sha256": _sha256(packed),
                "cell_angstrom": cell,
                "structures": structures,
            }
        )

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
        "schema_version": 2,
        "generator": "nio-md-prep prepare-agglomeration",
        "tool_version": __version__,
        "status": status,
        "sanity_status": (
            "PASS"
            if sanity_rows and all(row["status"] == "PASS" for row in sanity_rows)
            else "PENDING"
        ),
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "settings": {
            "replicas": replicas,
            "base_seed": base_seed,
            "radius_angstrom": radius,
            "packmol_tolerance_angstrom": tolerance,
            "vacuum_angstrom": vacuum,
            "minimum_distance_angstrom": minimum_distance,
            "center_scales": scales,
        },
        "vasp_reference_sanity": reference_validation,
        "output_mode": "structures_only" if structures_only else "vasp_runs",
        "case_root": case_root_name,
        "molecules": [
            {
                "slug": template.slug,
                "count": template.count,
                "atoms_per_molecule": len(template.symbols),
                "source": str(template.source.relative_to(ROOT)),
                "source_sha256": _sha256(template.source),
                "ligpargen_net_charge": template.net_charge,
            }
            for template in templates
        ],
        "replicas": replica_records,
    }
    manifest_path.write_text(json.dumps(root_manifest, indent=2) + "\n", encoding="utf-8")
    return manifest_path
