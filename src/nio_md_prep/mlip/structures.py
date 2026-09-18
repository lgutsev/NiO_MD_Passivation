"""ASE ``Atoms`` as the interchange structure for the MLIP layer.

For a full-system MLIP, symbols, coordinates, cell and PBC are everything a
potential needs, and ASE is the most direct route into MACE. So this layer
speaks ``Atoms`` -- and only this layer. The repository's topology-rich LAMMPS
representation (:mod:`nio_md_prep.lammps`) is untouched: bonds, angles,
dihedrals, impropers, molecule IDs and per-atom charges stay where the
classical workflows need them. :func:`from_lammps_data` is a *converter*, one
way, into the reduced view an MLIP consumes.

``ase`` is imported lazily inside the functions that need it, so importing
this module costs nothing.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from pathlib import Path

from .errors import ConfigError, MissingDependencyError

#: ``structure.format`` values this layer understands. ``auto`` defers to
#: ASE's own format detection; ``lammps-data-nio`` uses this repository's
#: parser plus its mass-to-element table, which is what makes a prepared
#: classical build readable by the MLIP layer without a round trip through ASE.
FORMATS = ("auto", "xyz", "extxyz", "cif", "vasp", "lammps-data", "lammps-data-nio")


def require_ase():
    """Import ASE, or explain what to install."""
    try:
        import ase  # noqa: F401
        from ase import Atoms  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without ase
        raise MissingDependencyError(
            "the MLIP structure layer",
            ("ase",),
            hint="Install the MLIP extra: pip install 'nio-md-prep[mlip]'",
        ) from exc
    import ase

    return ase


def load_structure(spec_or_path, *, format: str | None = None, index: int = -1):
    """Read a geometry into an ASE ``Atoms`` object.

    Accepts a :class:`~nio_md_prep.mlip.specs.StructureSpec` or a bare path.
    """
    path = getattr(spec_or_path, "path", spec_or_path)
    fmt = format if format is not None else getattr(spec_or_path, "format", "auto")
    idx = index if format is not None else getattr(spec_or_path, "index", index)
    path = Path(path)
    if fmt not in FORMATS:
        raise ConfigError(
            f"structure.format must be one of {', '.join(FORMATS)}; got {fmt!r}"
        )
    if not path.exists():
        raise FileNotFoundError(f"structure not found: {path}")
    if fmt == "lammps-data-nio":
        return from_lammps_data(path)
    require_ase()
    from ase.io import read

    kwargs = {} if fmt == "auto" else {"format": fmt}
    atoms = read(str(path), index=idx, **kwargs)
    if isinstance(atoms, list):  # a bare index=':' read
        atoms = atoms[-1]
    return atoms


def from_lammps_data(path: Path, *, type_map: Mapping[int, str] | None = None):
    """Convert a LigParGen/prepared LAMMPS data file into an ASE ``Atoms``.

    Uses this repository's own parser, so the cell bounds and coordinates are
    read exactly as the classical workflows read them. Elements come from
    ``type_map`` when the caller has one (a LAMMPS-native MLIP always does),
    and otherwise from :data:`nio_md_prep.geometry.ELEMENTS`, which infers
    them from the Masses section.

    Topology is deliberately dropped: a full-system MLIP does not consume
    bonds, and silently carrying them would suggest this conversion is
    reversible. It is not.
    """
    require_ase()
    from ase import Atoms

    from ..lammps import atom_coordinates, parse

    data = parse(Path(path))
    if type_map:
        symbols = [type_map[int(r.fields[2])] for r in data.sections["Atoms"]]
    else:
        from ..geometry import elements as infer_elements

        symbols = infer_elements(data)
    positions = atom_coordinates(data)
    cell = [data.bounds.get(axis, (0.0, 0.0))[1] - data.bounds.get(axis, (0.0, 0.0))[0]
            for axis in "xyz"]
    periodic = all(length > 0 for length in cell)
    return Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell if periodic else None,
        pbc=periodic,
    )


def to_lammps_type_order(atoms, type_map: Mapping[int, str]) -> list[int]:
    """Map each atom onto its LAMMPS type under ``type_map``.

    Raises before anything runs when the structure carries an element the
    type map has no slot for -- the LAMMPS analogue of element-coverage
    validation.
    """
    reverse: dict[str, int] = {}
    for type_id, symbol in sorted(type_map.items()):
        reverse.setdefault(symbol, int(type_id))
    missing = sorted({s for s in atoms.get_chemical_symbols() if s not in reverse})
    if missing:
        raise ConfigError(
            f"structure contains element(s) {', '.join(missing)} with no LAMMPS type in "
            f"potential.type_map (maps: {', '.join(sorted(reverse))})"
        )
    return [reverse[s] for s in atoms.get_chemical_symbols()]


def structure_elements(atoms) -> tuple[str, ...]:
    return tuple(sorted(set(atoms.get_chemical_symbols())))


def is_periodic(atoms) -> bool:
    return bool(getattr(atoms, "pbc", None) is not None and all(atoms.pbc))


def structure_digest(atoms) -> str:
    """A stable SHA256 over the interchange view of a structure.

    Hashes exactly what the MLIP layer consumes -- symbols, coordinates, cell,
    PBC -- rounded to 1e-8 Angstrom so that a round trip through a text format
    does not change the digest. Recorded in every manifest, so a result can
    always be tied back to the geometry that produced it.
    """
    payload = {
        "symbols": list(atoms.get_chemical_symbols()),
        "positions": [[round(float(c), 8) for c in row] for row in atoms.get_positions()],
        "cell": [[round(float(c), 8) for c in row] for row in atoms.get_cell()],
        "pbc": [bool(v) for v in atoms.pbc],
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def describe_structure(atoms) -> dict:
    """A compact, JSON-serialisable description for reports and manifests."""
    symbols = list(atoms.get_chemical_symbols())
    counts: dict[str, int] = {}
    for symbol in symbols:
        counts[symbol] = counts.get(symbol, 0) + 1
    cell = [[round(float(c), 8) for c in row] for row in atoms.get_cell()]
    volume = float(atoms.get_volume()) if is_periodic(atoms) else None
    return {
        "n_atoms": len(symbols),
        "elements": sorted(counts),
        "composition": dict(sorted(counts.items())),
        "periodic": is_periodic(atoms),
        "pbc": [bool(v) for v in atoms.pbc],
        "cell_angstrom": cell,
        "volume_angstrom3": volume,
        "sha256": structure_digest(atoms),
    }


def check_element_coverage(
    symbols: Iterable[str], supported: Iterable[str], *, label: str
) -> None:
    """Thin wrapper kept here so structure code reads in one place."""
    from .capabilities import CapabilitySet, check_elements

    check_elements(CapabilitySet(elements=frozenset(supported)), symbols, label=label)


__all__ = [
    "FORMATS",
    "require_ase",
    "load_structure",
    "from_lammps_data",
    "to_lammps_type_order",
    "structure_elements",
    "is_periodic",
    "structure_digest",
    "describe_structure",
    "check_element_coverage",
]
