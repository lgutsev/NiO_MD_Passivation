"""ASE ``Atoms`` as the interchange structure for the MLIP layer.

For a full-system MLIP, symbols, coordinates, cell and PBC are everything a
potential needs, and ASE is the most direct route into MACE. So this layer
speaks ``Atoms`` -- and only this layer. The repository's topology-rich LAMMPS
representation (:mod:`nio_md_prep.lammps`) is untouched: bonds, angles,
dihedrals, impropers, molecule IDs and per-atom charges stay where the
classical workflows need them. :func:`from_lammps_data` is a *converter*, one
way, into the reduced view an MLIP consumes.

Periodicity is **per cell axis** everywhere in this layer. ``atoms.pbc`` is a
three-tuple (one entry per cell vector a, b, c), and nothing here collapses
it to a single "periodic" boolean: a slab periodic along a and b only is a
different system from a cluster, and an engine that cannot tell them apart
must refuse the slab. :func:`periodic_axes` is the accessor;
:func:`is_periodic` survives with its old meaning, *fully* periodic.

Constraints are treated the same way: ASE ``FixAtoms`` is the one constraint
type the MLIP engines honour (:func:`fixed_atom_indices`); any other type is
refused by name rather than silently dropped by an engine that never looks.

``ase`` and ``numpy`` are imported lazily inside the functions that need
them, so importing this module costs nothing.
"""
from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from .errors import ConfigError, MissingDependencyError
from .specs import DEFAULT_VACUUM_GAP_THRESHOLD_ANGSTROM

#: ``structure.format`` values this layer understands. ``auto`` defers to
#: ASE's own format detection; ``lammps-data-nio`` uses this repository's
#: parser plus its mass-to-element table, which is what makes a prepared
#: classical build readable by the MLIP layer without a round trip through ASE.
FORMATS = ("auto", "xyz", "extxyz", "cif", "vasp", "lammps-data", "lammps-data-nio")
#: Formats whose files do not say which cell axes are periodic. A LAMMPS data
#: file has box bounds, but periodicity lives in the input deck's ``boundary``
#: command (and ASE's reader reports every such file as fully periodic), so
#: the caller must state ``pbc`` explicitly.
FORMATS_WITHOUT_PBC = ("lammps-data", "lammps-data-nio")

#: Rounding applied before hashing, so a text round trip keeps the digest.
_DIGEST_DECIMALS = 8


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


def load_structure(
    spec_or_path,
    *,
    format: str | None = None,
    index: int = -1,
    pbc: Sequence[bool] | None = None,
):
    """Read a geometry into an ASE ``Atoms`` object.

    Accepts a :class:`~nio_md_prep.mlip.specs.StructureSpec` or a bare path.
    ``pbc`` (or the spec's ``pbc``) replaces the periodicity the file
    declared; it is required for :data:`FORMATS_WITHOUT_PBC`. A cell axis
    declared periodic must have a non-zero cell vector.
    """
    path = getattr(spec_or_path, "path", spec_or_path)
    fmt = format if format is not None else getattr(spec_or_path, "format", "auto")
    idx = index if format is not None else getattr(spec_or_path, "index", index)
    if pbc is None:
        pbc = getattr(spec_or_path, "pbc", None)
    path = Path(path)
    if fmt not in FORMATS:
        raise ConfigError(
            f"structure.format must be one of {', '.join(FORMATS)}; got {fmt!r}"
        )
    if not path.exists():
        raise FileNotFoundError(f"structure not found: {path}")
    if fmt in FORMATS_WITHOUT_PBC and pbc is None:
        raise ConfigError(
            f"reading {path.name} as {fmt!r} needs an explicit pbc (structure.pbc): a "
            "LAMMPS data file does not record which box faces are periodic"
        )
    if fmt == "lammps-data-nio":
        return from_lammps_data(path, pbc=pbc)
    require_ase()
    from ase.io import read

    kwargs = {} if fmt == "auto" else {"format": fmt}
    atoms = read(str(path), index=idx, **kwargs)
    if isinstance(atoms, list):  # a bare index=':' read
        atoms = atoms[-1]
    if pbc is not None:
        atoms.pbc = _pbc_tuple(pbc)
    _check_periodic_cell(atoms)
    return atoms


def from_lammps_data(
    path: Path,
    *,
    pbc: Sequence[bool] | None = None,
    type_map: Mapping[int, str] | None = None,
):
    """Convert a LigParGen/prepared LAMMPS data file into an ASE ``Atoms``.

    Uses this repository's own parser, so coordinates are read exactly as
    the classical workflows read them. The box is taken completely:
    ``xlo xhi``/``ylo yhi``/``zlo zhi`` give the edge lengths, an ``xy xz yz``
    line gives the tilts of a (restricted) triclinic box, and the ``lo``
    corner is the origin -- ASE cells start at zero, so it is subtracted from
    every position. The resulting cell is the LAMMPS restricted form
    ``[[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]]``.

    ``pbc`` is required: which faces are periodic is stated by the input
    deck's ``boundary`` command, not by the data file (the classical slab
    builds use ``boundary p p f``, i.e. ``pbc=(True, True, False)``).

    Elements come from ``type_map`` when the caller has one (a LAMMPS-native
    MLIP always does), and otherwise from :data:`nio_md_prep.geometry.ELEMENTS`,
    which infers them from the Masses section. Only ``atom_style full``
    records (``id mol type q x y z ...``) are understood; any other declared
    style is refused rather than misread.

    Topology, partial charges, molecule IDs and image flags are deliberately
    dropped: a full-system MLIP does not consume them, and silently carrying
    them would suggest this conversion is reversible. It is not.
    """
    if pbc is None:
        raise ConfigError(
            f"converting {Path(path).name} needs an explicit pbc: a LAMMPS data file does "
            "not record which box faces are periodic (the classical slab builds use "
            "'boundary p p f', i.e. pbc = [true, true, false])"
        )
    pbc = _pbc_tuple(pbc)
    require_ase()
    import numpy as np
    from ase import Atoms

    from ..lammps import atom_coordinates, parse

    path = Path(path)
    header = _lammps_data_header(path)
    if header["atom_style"] not in (None, "full"):
        raise ConfigError(
            f"{path.name} declares 'Atoms # {header['atom_style']}'; only atom_style full "
            "data files (id mol type q x y z) can be converted by this reader"
        )
    data = parse(path)
    short = [r for r in data.sections["Atoms"] if len(r.fields) < 7]
    if short:
        raise ConfigError(
            f"{path.name}: Atoms record {short[0].render()!r} has fewer than the 7 fields "
            "of atom_style full (id mol type q x y z)"
        )
    missing = [axis for axis in "xyz" if axis not in data.bounds]
    if missing:
        raise ConfigError(
            f"{path.name} has no {' / '.join(f'{a}lo {a}hi' for a in missing)} line; the "
            "box must be complete to build a cell"
        )
    if type_map:
        try:
            symbols = [type_map[int(r.fields[2])] for r in data.sections["Atoms"]]
        except KeyError as exc:
            raise ConfigError(
                f"{path.name} uses LAMMPS atom type {exc.args[0]}, which the type map "
                f"({', '.join(f'{k}={v}' for k, v in sorted(type_map.items()))}) does not "
                "name"
            ) from None
    else:
        from ..geometry import elements as infer_elements

        symbols = infer_elements(data)
    lo = np.array([data.bounds[axis][0] for axis in "xyz"])
    lengths = np.array([data.bounds[axis][1] - data.bounds[axis][0] for axis in "xyz"])
    if np.any(lengths <= 0):
        raise ConfigError(f"{path.name}: box lengths {lengths.tolist()} must all be positive")
    xy, xz, yz = header["tilt"]
    cell = [
        [lengths[0], 0.0, 0.0],
        [xy, lengths[1], 0.0],
        [xz, yz, lengths[2]],
    ]
    positions = np.array(atom_coordinates(data), dtype=float) - lo
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)
    return atoms


_TILT = re.compile(
    r"^\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+xy\s+xz\s+yz\b"
)
_GENERAL_TRICLINIC = re.compile(r"^\s*[-+\d.eE]+\s+[-+\d.eE]+\s+[-+\d.eE]+\s+(avec|bvec|cvec)\b")
_ATOMS_HEADER = re.compile(r"^\s*Atoms\s*(?:#\s*(\S+))?\s*$")


def _lammps_data_header(path: Path) -> dict:
    """The box tilt and declared atom style, which the repository parser skips."""
    tilt = (0.0, 0.0, 0.0)
    atom_style = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        handle.readline()  # title line; free text, never parsed
        for line in handle:
            match = _ATOMS_HEADER.match(line)
            if match:
                atom_style = match.group(1)
                break
            if _GENERAL_TRICLINIC.match(line):
                raise ConfigError(
                    f"{path.name} uses the general-triclinic avec/bvec/cvec header, which "
                    "this reader does not convert; write the box in restricted form "
                    "(xlo xhi / ylo yhi / zlo zhi plus xy xz yz)"
                )
            match = _TILT.match(line)
            if match:
                tilt = tuple(float(match.group(i)) for i in (1, 2, 3))
    return {"tilt": tilt, "atom_style": atom_style}


def to_lammps_type_order(atoms, type_map: Mapping[int, str]) -> list[int]:
    """Map each atom onto its LAMMPS type under ``type_map``.

    Raises before anything runs when the structure carries an element the
    type map has no slot for -- the LAMMPS analogue of element-coverage
    validation. A type map maps each element to exactly one type (enforced
    by :class:`~nio_md_prep.mlip.specs.LammpsMlipPotentialSpec`), so the
    mapping is unambiguous.
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


# ---------------------------------------------------------------------------
# Periodicity
# ---------------------------------------------------------------------------


def periodic_axes(atoms) -> tuple[bool, bool, bool]:
    """Per-cell-axis periodicity, ``(pbc_a, pbc_b, pbc_c)``."""
    pbc = getattr(atoms, "pbc", None)
    if pbc is None:
        return (False, False, False)
    return _pbc_tuple(pbc)


def is_periodic(atoms) -> bool:
    """True only when *every* cell axis is periodic (bulk-like cells).

    Kept under its original name and meaning. Code that must distinguish a
    slab from a cluster uses :func:`periodic_axes` instead.
    """
    return all(periodic_axes(atoms))


def is_partially_periodic(atoms) -> bool:
    """Periodic along some cell axes but not all (a slab or a wire)."""
    axes = periodic_axes(atoms)
    return any(axes) and not all(axes)


def _pbc_tuple(pbc) -> tuple[bool, bool, bool]:
    if isinstance(pbc, bool):
        return (pbc, pbc, pbc)
    values = tuple(bool(v) for v in pbc)
    if len(values) != 3:
        raise ConfigError(f"pbc must have one entry per cell axis; got {pbc!r}")
    return values  # type: ignore[return-value]


def _check_periodic_cell(atoms) -> None:
    import numpy as np

    lengths = np.linalg.norm(np.array(atoms.get_cell()), axis=1)
    empty = [name for name, periodic, length in zip("abc", periodic_axes(atoms), lengths)
             if periodic and length == 0.0]
    if empty:
        raise ConfigError(
            f"structure is declared periodic along cell axis {', '.join(empty)} but that "
            "cell vector has zero length; give the cell or set pbc false for that axis"
        )


def largest_empty_gaps(atoms) -> tuple[float | None, float | None, float | None]:
    """Widest atom-free slab along each periodic lattice direction, in Angstrom.

    For cell axis *i*, the atoms' fractional coordinates along *i* are
    wrapped into [0, 1) and sorted; the largest spacing between neighbours
    (including the wrap from the last back to the first) is the widest empty
    slab bounded by lattice planes of constant fractional coordinate. It is
    converted to a distance with the interplanar spacing ``1 / |b_i|``
    (``b_i`` the reciprocal vector without the 2 pi), i.e. measured normal to
    the plane spanned by the other two cell vectors. ``None`` for
    non-periodic axes, where the cell edge is arbitrary.

    This is a geometric screen, not a physical definition of vacuum: a
    well-packed bulk crystal has gaps of about one interplanar spacing, a
    slab with 15 Angstrom of vacuum has one of about 15 Angstrom.
    """
    import numpy as np

    axes = periodic_axes(atoms)
    if not any(axes) or len(atoms) == 0:
        return (None, None, None)
    cell = atoms.cell.complete()
    fractional = np.linalg.solve(np.array(cell).T, np.array(atoms.get_positions()).T).T
    reciprocal = np.linalg.inv(np.array(cell)).T  # rows b_i with a_i . b_j = delta_ij
    gaps: list[float | None] = []
    for axis in range(3):
        if not axes[axis]:
            gaps.append(None)
            continue
        coordinate = np.sort(np.mod(fractional[:, axis], 1.0))
        spacing = np.diff(coordinate)
        wrap = coordinate[0] + 1.0 - coordinate[-1]
        widest = max(float(spacing.max()) if spacing.size else 0.0, float(wrap))
        gaps.append(widest / float(np.linalg.norm(reciprocal[axis])))
    return tuple(gaps)  # type: ignore[return-value]


def vacuum_gaps(
    atoms, *, threshold_angstrom: float = DEFAULT_VACUUM_GAP_THRESHOLD_ANGSTROM
) -> dict[int, float]:
    """``axis -> width`` for every periodic axis whose widest gap exceeds the threshold.

    An empty mapping means the cell looks bulk-like along every periodic
    axis. Used to refuse isotropic/anisotropic NPT on a vacuum slab that
    was written as a fully periodic cell (every POSCAR slab is).
    """
    return {
        axis: width
        for axis, width in enumerate(largest_empty_gaps(atoms))
        if width is not None and width > threshold_angstrom
    }


def outside_nonperiodic_faces(atoms) -> dict[int, tuple[int, ...]]:
    """Atoms outside the cell along a non-periodic axis, ``axis -> indices``.

    An engine that honours per-axis periodicity with fixed (not
    shrink-wrapped) faces -- LAMMPS ``boundary f`` -- loses or rejects an
    atom that lies outside ``[0, 1)`` in fractional coordinates along a
    non-periodic axis. This helper finds them so the engine can refuse the
    structure with the offending indices instead of wrapping, shifting or
    shrinking anything. Periodic axes are never reported (wrapping there is
    harmless). Needs a rank-3 cell: no box is synthesised for a cluster that
    has none.
    """
    import numpy as np

    axes = periodic_axes(atoms)
    if all(axes):
        return {}
    cell = np.array(atoms.get_cell())
    if atoms.cell.rank < 3 or abs(np.linalg.det(cell)) == 0.0:
        raise ConfigError(
            "this structure has no rank-3 cell, so atoms cannot be placed along its "
            "non-periodic axes; give it an explicit cell (box) -- none is synthesised"
        )
    fractional = np.linalg.solve(cell.T, np.array(atoms.get_positions()).T).T
    offenders: dict[int, tuple[int, ...]] = {}
    for axis in range(3):
        if axes[axis]:
            continue
        s = fractional[:, axis]
        bad = np.nonzero((s < 0.0) | (s >= 1.0))[0]
        if bad.size:
            offenders[axis] = tuple(int(i) for i in bad)
    return offenders


# ---------------------------------------------------------------------------
# Constraints
# ---------------------------------------------------------------------------


def fixed_atom_indices(atoms) -> tuple[int, ...]:
    """Sorted indices frozen by ASE ``FixAtoms``; refuse any other constraint.

    ``FixAtoms`` is the one constraint type the MLIP engines are expected to
    honour (frozen atoms in velocity initialisation, integration and the
    temperature degrees of freedom). Any other ASE constraint --
    ``FixBondLengths``, ``FixCom``, ``Hookean``, ``FixCartesian`` ... -- is
    refused by type name: an engine that ignores it would silently integrate
    different physics.
    """
    fixed: set[int] = set()
    unsupported: list[str] = []
    for constraint in getattr(atoms, "constraints", None) or ():
        if type(constraint).__name__ == "FixAtoms":
            fixed.update(int(i) for i in _fixatoms_indices(constraint))
        else:
            unsupported.append(type(constraint).__name__)
    if unsupported:
        raise ConfigError(
            f"structure carries ASE constraint(s) {', '.join(sorted(set(unsupported)))}; "
            "only FixAtoms (frozen atoms) is honoured by the MLIP engines. Remove the "
            "other constraints before running this job."
        )
    out_of_range = sorted(i for i in fixed if not 0 <= i < len(atoms))
    if out_of_range:
        raise ConfigError(
            f"FixAtoms names atom indices {out_of_range[:10]} outside 0..{len(atoms) - 1}"
        )
    return tuple(sorted(fixed))


def constraint_summary(atoms) -> dict:
    """A JSON-friendly record of the structure's constraints. Never raises."""
    fixed: set[int] = set()
    other: list[str] = []
    for constraint in getattr(atoms, "constraints", None) or ():
        if type(constraint).__name__ == "FixAtoms":
            fixed.update(int(i) for i in _fixatoms_indices(constraint))
        else:
            other.append(type(constraint).__name__)
    return {
        "fixed_atoms": sorted(fixed),
        "n_fixed": len(fixed),
        "unsupported": sorted(set(other)),
    }


def _fixatoms_indices(constraint) -> list[int]:
    getter = getattr(constraint, "get_indices", None)
    indices = getter() if callable(getter) else constraint.index
    return [int(i) for i in indices]


# ---------------------------------------------------------------------------
# Identity and description
# ---------------------------------------------------------------------------


def structure_digest(atoms) -> str:
    """A stable SHA256 over the interchange view of a structure.

    Hashes what the MLIP layer consumes, each value rounded to 1e-8 so that a
    round trip through a text format does not change the digest:

    - chemical symbols, positions, cell and per-axis PBC (always);
    - constraints, when any are attached: the sorted ``FixAtoms`` indices,
      and the type name plus ``todict()`` of any other constraint;
    - ``initial_magmoms``, ``initial_charges`` and ``tags``, each only when
      present *and* not all zero. ASE treats an absent array exactly like an
      all-zero one, so the two must hash alike; a ferromagnetic and an
      antiferromagnetic initialisation of the same NiO geometry must not.

    Still order-sensitive and image-sensitive: it identifies this exact
    ``Atoms`` object, not an equivalence class of structures. Recorded in
    every manifest, so a result can be tied back to the geometry, magnetic
    initialisation and constraints that produced it.
    """
    import numpy as np

    payload: dict = {
        "symbols": list(atoms.get_chemical_symbols()),
        "positions": _rounded(atoms.get_positions()),
        "cell": _rounded(atoms.get_cell()),
        "pbc": list(periodic_axes(atoms)),
    }
    constraints = _constraint_records(atoms)
    if constraints:
        payload["constraints"] = constraints
    for name in ("initial_magmoms", "initial_charges"):
        if atoms.has(name):
            values = np.asarray(atoms.arrays[name], dtype=float)
            if np.any(values != 0.0):
                payload[name] = _rounded(values)
    if atoms.has("tags"):
        tags = np.asarray(atoms.arrays["tags"])
        if np.any(tags != 0):
            payload["tags"] = [int(t) for t in tags]
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def _rounded(values):
    import numpy as np

    array = np.asarray(values, dtype=float)
    return np.round(array, _DIGEST_DECIMALS).tolist()


def _constraint_records(atoms) -> list:
    records = []
    for constraint in getattr(atoms, "constraints", None) or ():
        name = type(constraint).__name__
        if name == "FixAtoms":
            records.append({"type": name, "indices": sorted(_fixatoms_indices(constraint))})
            continue
        try:
            spec = _jsonable(constraint.todict())
        except Exception:  # pragma: no cover - a constraint without todict
            spec = repr(constraint)
        records.append({"type": name, "spec": spec})
    return sorted(records, key=lambda r: json.dumps(r, sort_keys=True))


def _jsonable(value):
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if hasattr(value, "tolist"):
        return _jsonable(value.tolist())
    if isinstance(value, float):
        return round(value, _DIGEST_DECIMALS)
    return value


def describe_structure(atoms) -> dict:
    """A compact, JSON-serialisable description for reports and manifests.

    ``periodic`` keeps its original meaning (fully periodic); ``pbc`` is the
    per-axis truth. ``constraints`` lists frozen atoms and names any
    unsupported constraint types (which a job will refuse), and
    ``largest_empty_gap_angstrom`` is :func:`largest_empty_gaps`, recorded
    so the NPT vacuum decision can be audited.
    """
    symbols = list(atoms.get_chemical_symbols())
    counts: dict[str, int] = {}
    for symbol in symbols:
        counts[symbol] = counts.get(symbol, 0) + 1
    cell = [[round(float(c), 8) for c in row] for row in atoms.get_cell()]
    volume = float(atoms.get_volume()) if is_periodic(atoms) else None
    gaps = largest_empty_gaps(atoms)
    return {
        "n_atoms": len(symbols),
        "elements": sorted(counts),
        "composition": dict(sorted(counts.items())),
        "periodic": is_periodic(atoms),
        "pbc": list(periodic_axes(atoms)),
        "cell_angstrom": cell,
        "volume_angstrom3": volume,
        "largest_empty_gap_angstrom": [None if g is None else round(g, 6) for g in gaps],
        "constraints": constraint_summary(atoms),
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
    "FORMATS_WITHOUT_PBC",
    "DEFAULT_VACUUM_GAP_THRESHOLD_ANGSTROM",
    "require_ase",
    "load_structure",
    "from_lammps_data",
    "to_lammps_type_order",
    "structure_elements",
    "periodic_axes",
    "is_periodic",
    "is_partially_periodic",
    "largest_empty_gaps",
    "vacuum_gaps",
    "outside_nonperiodic_faces",
    "fixed_atom_indices",
    "constraint_summary",
    "structure_digest",
    "describe_structure",
    "check_element_coverage",
]
