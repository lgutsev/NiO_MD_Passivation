"""Streaming, fail-closed reader for VASP ``vasprun.xml`` (plain, .gz, .bz2, .xz).

The reader transcribes the file; it never repairs, infers or pairs data across
blocks. Facts it relies on (see ``design/understand/vasp-format-research.md``
S1/S2/S3/S8, measured on real VASP 5.2-6.5 output):

* Top level: ``generator``, ``incar``, [``primitive_cell``], ``kpoints``,
  ``parameters``, ``atominfo``, ``structure name="initialpos"``, then one block
  per ionic step, then ``structure name="finalpos"`` and ``</modeling>``.
* A DFT ionic step is a ``<calculation>`` holding ``<scstep>`` x n, the
  ``<structure>`` at which its forces/stress/energy were evaluated, ``forces``,
  [``stress``], ``<energy>`` and ``<time>`` (plus DOS/eigenvalue blocks that are
  skipped here). Labels and geometry therefore always come from the same block.
* VASP-MLFF force-field-only steps are FLAT depth-1 sequences
  ``structure`` (no name) -> ``varray forces`` -> [``varray stress``] ->
  ``energy`` -> ``time`` interleaved with ``<calculation>`` blocks. They are
  yielded in file order with ``label_source="mlff"``; ``IonicStep.index``
  counts ALL steps (DFT and MLFF).
* ``finalpos`` is the geometry after the last update, not the last
  force-evaluated one; it is exposed on the trailer and never paired with labels.
* Interrupted runs end inside a block. The reader stops at the parse error,
  keeps every completed step, yields the in-progress block with
  ``complete=False`` and never synthesizes a frame from ``initialpos``.

Derived values (all pure functions of parsed blocks):

* :func:`derived_energies` -- F = calc ``e_fr_energy`` - PSTRESS*V (the
  force-consistent label) and E0/E_wo reconstructed from the last ``<scstep>``,
  version-gated against the VASP <= 6.0.8 calc-level mislabel, with
  ``energy_source``/``energy_rule``/parser provenance for every frame;
* :func:`scf_summary` -- electronic iterations and the last |dE|;
* :func:`resolve_selective` -- selective-dynamics flags (DIRECT basis) with the
  precedence initialpos > finalpos > POSCAR > CONTCAR and conflict detection;
* :func:`mlff_settings` -- VASP-MLFF state from ``<incar>``.

Memory stays bounded: ``ET.iterparse`` with start/end events, a depth-tracked
element stack, processed depth-1 children removed from the root, processed
children of a ``<calculation>`` removed from it, and uninteresting subtrees
(eigenvalues, DOS, projections, dielectric data) pruned element by element as
they are read.

Every problem string on a header, step or trailer starts with a
:data:`~nio_md_prep.dataset.model.REASONS` code (``non_finite``, ``bad_shape``,
``species_mismatch``, ``missing_energy`` ...) or ``missing``/``info``, then
``": "`` and a free-text detail, so the policy layer can map them mechanically.
"""

from __future__ import annotations

import bz2
import gzip
import hashlib
import lzma
import math
import re
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path
from typing import Any, BinaryIO, Callable, Iterator

from . import PARSER_NAME, PARSER_VERSION
from .errors import DatasetError, DependencyMissingError
from .model import AtomType, IonicStep, StepStructure, VasprunHeader, VasprunTrailer

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - exercised only without numpy
    raise DependencyMissingError(
        "numpy is required to read vasprun.xml; install nio-md-prep[dataset]"
    ) from exc


#: VASP's own electron-volt -> joule constant (VASP ``constant.F``:
#: ``EVTOJ=1.60217733E-19``). VASP adds ``PSTRESS*V`` to the calculation-level
#: ``e_fr_energy`` with this constant, and ASE 3.29 removes it with the same
#: constant (``ase/io/vasp.py``: ``pressure *= 1e-22 / EVTOJ``). Using it
#: reproduces ASE's free energy for ``vasprun_pstress.xml`` bit for bit
#: (-20.246959373169393 eV). The CODATA factor 1/1602.1766208 used for stress
#: labels differs by 4.4e-7 relative (3.1e-9 eV for that 11.4 A^3 cell).
VASP_EVTOJ = 1.60217733e-19
#: eV of a PSTRESS*V term per (kB * A^3), exactly as VASP computes it.
PV_EV_PER_KBAR_A3 = 1e-22 / VASP_EVTOJ

_COMPRESSION_MAGIC = ((b"\x1f\x8b", "gz"), (b"BZh", "bz2"), (b"\xfd7zXZ\x00", "xz"))
_VERSION = re.compile(r"\s*(\d+)\.(\d+)\.(\d+)")
_READ_ERRORS = (ET.ParseError, EOFError, OSError, zlib.error, lzma.LZMAError)

#: children of <calculation> that are transcribed; all other subtrees are pruned while streaming
_CALC_KEPT = frozenset({"scstep", "structure", "varray", "energy", "time"})
#: depth-1 elements whose subtree is kept until their end event
_TOP_KEPT = frozenset(
    {"generator", "incar", "kpoints", "parameters", "atominfo", "structure", "varray", "energy", "time"}
)
_HEADER_STRUCTURES = frozenset({"initialpos", "primitive_cell"})


class VasprunParseError(DatasetError):
    """The vasprun.xml header (generator ... atominfo) could not be read."""


# --------------------------------------------------------------------------
# Opening (compression by magic bytes) with a single-pass sha256
# --------------------------------------------------------------------------

def detect_compression(path: Path) -> str | None:
    """``"gz"``, ``"bz2"``, ``"xz"`` or ``None`` from the file's magic bytes (not its suffix)."""
    with Path(path).open("rb") as handle:
        head = handle.read(6)
    for magic, name in _COMPRESSION_MAGIC:
        if head.startswith(magic):
            return name
    return None


def open_maybe_compressed(path: Path) -> BinaryIO:
    """Open a plain/.gz/.bz2/.xz file for binary reading, choosing by magic bytes.

    Bytes are returned undecoded: vasprun.xml declares ISO-8859-1 in its XML
    header and expat must see the raw bytes; text readers decode themselves.
    """
    path = Path(path)
    compression = detect_compression(path)
    if compression == "gz":
        return gzip.open(path, "rb")
    if compression == "bz2":
        return bz2.open(path, "rb")
    if compression == "xz":
        return lzma.open(path, "rb")
    return path.open("rb")


class _HashingReader:
    """File-like wrapper hashing the raw (possibly compressed) bytes as they are read."""

    def __init__(self, handle: BinaryIO):
        self._handle = handle
        self._digest = hashlib.sha256()
        self.bytes = 0

    def readable(self) -> bool:
        return True

    def read(self, size: int = -1) -> bytes:
        data = self._handle.read(size)
        self._digest.update(data)
        self.bytes += len(data)
        return data

    def readinto(self, buffer) -> int:
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)

    def drain(self) -> None:
        while self.read(1 << 16):
            pass

    def hexdigest(self) -> str:
        return self._digest.hexdigest()

    def close(self) -> None:
        self._handle.close()


class HashedSource:
    """Decompressed binary stream plus the sha256 of the stored bytes (read once).

    Call :meth:`finish` after reading to hash any bytes the consumer did not
    read (e.g. after a parse error); ``sha256``/``nbytes`` are then final.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self.compression = detect_compression(self.path)
        self._raw = self.path.open("rb")
        self._hasher = _HashingReader(self._raw)
        self._decompressor: Any = None
        if self.compression == "gz":
            self._decompressor = gzip.GzipFile(fileobj=self._hasher, mode="rb")
        elif self.compression == "bz2":
            self._decompressor = bz2.BZ2File(self._hasher, mode="rb")
        elif self.compression == "xz":
            self._decompressor = lzma.LZMAFile(self._hasher, mode="rb")
        # read1(): at most one raw read per call, so every byte decompressed before
        # a truncated/corrupt tail reaches the parser before the error is raised
        # (a buffered read(n) would drop the partial chunk with the exception).
        self.stream: Any = _Read1(self._decompressor) if self._decompressor is not None else self._hasher
        self.sha256: str | None = None
        self.nbytes: int | None = None

    def finish(self) -> None:
        if self.sha256 is None and not self._raw.closed:
            self._hasher.drain()
            self.sha256 = self._hasher.hexdigest()
            self.nbytes = self._hasher.bytes

    def close(self) -> None:
        try:
            if self._decompressor is not None:
                self._decompressor.close()
        except _READ_ERRORS:  # a corrupt compressed tail may fail again on close
            pass
        self._raw.close()


class _Read1:
    """``read(n)`` -> ``read1(n)`` adapter for the decompressing readers."""

    def __init__(self, stream: Any):
        self._stream = stream

    def readable(self) -> bool:
        return True

    def read(self, size: int = -1) -> bytes:
        return self._stream.read1(size if size and size > 0 else 1 << 16)


# --------------------------------------------------------------------------
# Small typed-value helpers
# --------------------------------------------------------------------------

def parse_version(text: str | None) -> tuple[int, int, int] | None:
    """``"6.3.0"`` / ``"5.4.4.18Apr17-6-g9f103f2a35"`` -> ``(6, 3, 0)`` / ``(5, 4, 4)``."""
    if not text:
        return None
    match = _VERSION.match(text)
    return tuple(int(part) for part in match.groups()) if match else None  # type: ignore[return-value]


def _is_overflow(token: str) -> bool:
    return bool(token) and set(token) <= {"*"}


def _logical(token: str) -> bool | None:
    value = token.strip().strip(".").upper()
    if value in {"T", "TRUE"}:
        return True
    if value in {"F", "FALSE"}:
        return False
    return None


def _float(token: str, problems: list[str], what: str) -> float:
    """Float or NaN; records ``non_finite`` for overflow, junk, NaN or inf."""
    try:
        value = float(token)
    except ValueError:
        problems.append(f"non_finite: {what} is not a number ({token.strip()!r})")
        return math.nan
    if not math.isfinite(value):
        problems.append(f"non_finite: {what} is {token.strip()!r}")
    return value


def _typed_scalar(text: str, kind: str | None, problems: list[str], what: str) -> Any:
    """Typed ``<i>``/``<v>`` token; VASP overflow ``*****`` -> None (recorded)."""
    token = text.strip()
    if kind == "string":
        return token
    if _is_overflow(token):
        problems.append(f"non_finite: {what} printed as {token!r} (stored as None)")
        return None
    if kind == "logical":
        value = _logical(token)
        if value is None:
            problems.append(f"bad_value: {what} is not a logical ({token!r})")
        return value
    if kind == "int":
        try:
            return int(token)
        except ValueError:
            problems.append(f"bad_value: {what} is not an integer ({token!r})")
            return None
    try:  # no type attribute: VASP writes floats
        return float(token)
    except ValueError:
        problems.append(f"bad_value: {what} is not a number ({token!r})")
        return None


def _typed_element(elem: ET.Element, problems: list[str], what: str) -> Any:
    kind = elem.get("type")
    text = elem.text or ""
    if elem.tag == "v":
        if kind == "string":
            return text.strip()
        return [_typed_scalar(token, kind, problems, what) for token in text.split()]
    return _typed_scalar(text, kind, problems, what)


def _varray(elem: ET.Element, problems: list[str], what: str) -> Any:
    """``<varray>`` -> (rows, width) float64 array, or None when rows are ragged."""
    rows = [(v.text or "").split() for v in elem.findall("v")]
    widths = sorted({len(row) for row in rows})
    if len(widths) > 1:
        problems.append(f"bad_shape: {what} rows have {widths} values")
        return None
    width = widths[0] if widths else 0
    flat = [token for row in rows for token in row]
    try:
        array = np.array(flat, dtype=np.float64)
        if not np.all(np.isfinite(array)):
            problems.append(f"non_finite: {what} contains NaN or inf")
    except ValueError:
        local: list[str] = []
        array = np.array([_float(token, local, what) for token in flat], dtype=np.float64)
        problems.append(local[0] if local else f"non_finite: {what} contains non-numeric values")
    return array.reshape(len(rows), width)


def _logical_varray(elem: ET.Element, problems: list[str], what: str) -> Any:
    rows = [(v.text or "").split() for v in elem.findall("v")]
    flags: list[list[bool]] = []
    for number, row in enumerate(rows, 1):
        values = [_logical(token) for token in row]
        if len(values) != 3 or any(value is None for value in values):
            problems.append(f"bad_shape: {what} row {number} is not three logicals ({' '.join(row)!r})")
            return None
        flags.append([bool(value) for value in values])
    return np.array(flags, dtype=bool).reshape(len(flags), 3)


def _energy_block(elem: ET.Element, problems: list[str], what: str) -> dict[str, float]:
    values: dict[str, float] = {}
    for item in elem.findall("i"):
        name = (item.get("name") or "").strip()
        values[name] = _float(item.text or "", problems, f"{what} {name}")
    return values


def _times(elem: ET.Element) -> list[float]:
    values = []
    for token in (elem.text or "").split():
        try:  # timings are never labels: overflow is not a frame problem
            values.append(float(token))
        except ValueError:
            values.append(math.nan)
    return values


def _structure(elem: ET.Element, problems: list[str], what: str) -> StepStructure:
    """``<structure>``: basis rows (A), fractional positions, optional selective flags.

    Velocities (initialpos/finalpos of MD runs) and ``<nose>`` are ignored.
    """
    cell = fractional = positions = selective = None
    volume: float | None = None
    basis = elem.find("crystal/varray[@name='basis']")
    if basis is None:
        problems.append(f"bad_shape: {what} has no crystal/basis")
    else:
        cell = _varray(basis, problems, f"{what} basis")
        if cell is not None and cell.shape != (3, 3):
            problems.append(f"bad_shape: {what} basis has shape {cell.shape}")
    volume_elem = elem.find("crystal/i[@name='volume']")
    if volume_elem is not None:
        volume = _float(volume_elem.text or "", problems, f"{what} volume")
    frac_elem = elem.find("varray[@name='positions']")
    if frac_elem is None:
        problems.append(f"bad_shape: {what} has no positions")
    else:
        fractional = _varray(frac_elem, problems, f"{what} positions")
        if fractional is not None and fractional.ndim == 2 and fractional.shape[1] != 3:
            problems.append(f"bad_shape: {what} positions have {fractional.shape[1]} columns")
    if (
        cell is not None and cell.shape == (3, 3)
        and fractional is not None and fractional.ndim == 2 and fractional.shape[1] == 3
    ):
        positions = fractional @ cell
    sel_elem = elem.find("varray[@name='selective']")
    if sel_elem is not None:
        selective = _logical_varray(sel_elem, problems, f"{what} selective")
    return StepStructure(cell=cell, fractional=fractional, positions=positions, selective=selective, volume=volume)


# --------------------------------------------------------------------------
# Header sections
# --------------------------------------------------------------------------

def _generator(elem: ET.Element) -> dict[str, str]:
    return {(item.get("name") or "").strip(): (item.text or "").strip() for item in elem.findall("i")}


def _flat_items(elem: ET.Element, problems: list[str], section: str) -> dict[str, Any]:
    return {
        (item.get("name") or "").strip(): _typed_element(item, problems, f"{section}.{item.get('name')}")
        for item in elem
        if item.tag in {"i", "v"}
    }


def parse_parameters(
    elem: ET.Element, problems: list[str] | None = None, duplicates: list[str] | None = None
) -> dict[str, Any]:
    """Flatten ``<parameters>`` (nested ``<separator>`` blocks) into name -> typed value.

    Types: ``int`` -> int, ``logical`` (" T "/" F ") -> bool, ``string`` ->
    stripped str, no type -> float; ``<v>`` -> list. VASP's ``*****`` overflow
    -> None. Names repeated in different separators (e.g. OMEGAMAX, CSHIFT)
    keep their first value and are listed in ``duplicates``.
    """
    problems = [] if problems is None else problems
    duplicates = [] if duplicates is None else duplicates
    flat: dict[str, Any] = {}

    def walk(node: ET.Element) -> None:
        for child in node:
            if child.tag == "separator":
                walk(child)
            elif child.tag in {"i", "v"}:
                name = (child.get("name") or "").strip()
                value = _typed_element(child, problems, f"parameters.{name}")
                if name in flat:
                    if name not in duplicates:
                        duplicates.append(name)
                else:
                    flat[name] = value

    walk(elem)
    return flat


def _kpoints(elem: ET.Element, problems: list[str]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "scheme": None, "divisions": None, "usershift": None, "shift": None, "nkpts": None,
    }
    generation = elem.find("generation")
    if generation is not None:
        result["scheme"] = generation.get("param")
        for vector in generation.findall("v"):
            name = vector.get("name")
            if name in {"divisions", "usershift", "shift"}:
                result[name] = _typed_element(vector, problems, f"kpoints.{name}")
    kpointlist = elem.find("varray[@name='kpointlist']")
    if kpointlist is not None:
        result["nkpts"] = len(kpointlist.findall("v"))
    return result


def _array_rows(array: ET.Element) -> tuple[list[str], list[list[str]]]:
    fields = [(field.text or "").strip() for field in array.findall("field")]
    rows = [[(cell.text or "") for cell in rc.findall("c")] for rc in array.findall("set/rc")]
    return fields, rows


def _atominfo(elem: ET.Element, problems: list[str]) -> tuple[list[str], list[AtomType]]:
    arrays = {array.get("name"): array for array in elem.findall("array")}
    species: list[str] = []
    atom_types: list[AtomType] = []
    if "atoms" not in arrays or "atomtypes" not in arrays:
        raise VasprunParseError("atominfo lacks the 'atoms' or 'atomtypes' array")
    fields, rows = _array_rows(arrays["atoms"])
    try:
        element_col = fields.index("element")
        type_col = fields.index("atomtype")
    except ValueError as exc:
        raise VasprunParseError(f"atominfo 'atoms' array has unexpected fields {fields}") from exc
    atom_type_index = []
    for row in rows:
        species.append(row[element_col].strip())
        try:
            atom_type_index.append(int(row[type_col]))
        except ValueError:
            atom_type_index.append(-1)
    fields, rows = _array_rows(arrays["atomtypes"])
    column = {name: position for position, name in enumerate(fields)}
    for row in rows:
        def cell(name: str) -> str | None:
            return row[column[name]] if name in column and column[name] < len(row) else None

        def number(name: str) -> float | None:
            text = cell(name)
            try:
                return float(text) if text is not None else None
            except ValueError:
                problems.append(f"bad_value: atomtypes {name} is not a number ({text!r})")
                return None

        try:
            count = int(cell("atomspertype") or "")
        except ValueError:
            problems.append(f"bad_value: atomtypes atomspertype is not an integer ({cell('atomspertype')!r})")
            count = -1
        atom_types.append(AtomType(
            element=(cell("element") or "").strip(),
            count=count,
            mass=number("mass"),
            valence=number("valence"),
            pseudopotential=(cell("pseudopotential") or "").strip(),
        ))
    declared = elem.findtext("atoms")
    try:
        n_declared = int(declared) if declared is not None else None
    except ValueError:
        n_declared = None
    if n_declared is not None and n_declared != len(species):
        problems.append(f"species_mismatch: atominfo declares {n_declared} atoms but lists {len(species)}")
    expanded = [atom_type.element for atom_type in atom_types for _ in range(max(atom_type.count, 0))]
    if expanded != species:
        problems.append("species_mismatch: atominfo 'atoms' rows disagree with 'atomtypes' counts/order")
    expected_index = [number for number, atom_type in enumerate(atom_types, 1) for _ in range(max(atom_type.count, 0))]
    if atom_type_index != expected_index and expanded == species:
        problems.append("species_mismatch: atominfo atomtype indices are not contiguous per type")
    return species, atom_types


# --------------------------------------------------------------------------
# The streaming reader
# --------------------------------------------------------------------------

class _StepBuilder:
    def __init__(self, index: int, label_source: str):
        self.index = index
        self.label_source = label_source
        self.structure: StepStructure | None = None
        self.forces: Any = None
        self.forces_seen = False
        self.stress: Any = None
        self.energies: dict[str, float] = {}
        self.energy_seen = False
        self.scf: list[dict[str, float]] = []
        self.times: dict[str, list[float]] = {}
        self.problems: list[str] = []

    @property
    def what(self) -> str:
        return f"step {self.index}"

    def take(self, elem: ET.Element) -> None:
        """Transcribe one finished child block (a child of <calculation>, or a flat ML block)."""
        tag = elem.tag
        if tag == "scstep":
            energy = elem.find("energy")
            self.scf.append(
                _energy_block(energy, self.problems, f"{self.what} scstep {len(self.scf) + 1}")
                if energy is not None else {}
            )
        elif tag == "structure":
            if self.structure is not None:
                self.problems.append(f"bad_shape: {self.what} has more than one <structure>")
            self.structure = _structure(elem, self.problems, f"{self.what} structure")
        elif tag == "varray":
            name = elem.get("name")
            if name == "forces":
                self.forces_seen = True
                self.forces = _varray(elem, self.problems, f"{self.what} forces")
            elif name == "stress":
                self.stress = _varray(elem, self.problems, f"{self.what} stress")
        elif tag == "energy":
            self.energies = _energy_block(elem, self.problems, f"{self.what} energy")
            self.energy_seen = True
        elif tag == "time":
            self.times[elem.get("name") or ""] = _times(elem)

    def build(self, complete: bool) -> IonicStep:
        if complete:
            if not self.forces_seen:
                self.problems.append(f"missing_forces: {self.what} has no forces varray")
            if not self.energy_seen:
                self.problems.append(f"missing_energy: {self.what} has no <energy> block")
            if self.structure is None:
                self.problems.append(f"bad_shape: {self.what} has no <structure>")
        return IonicStep(
            index=self.index,
            complete=complete,
            structure=self.structure,
            forces=self.forces,
            stress_kbar_vasp=self.stress,
            energies=self.energies,
            scf_energies=self.scf,
            problems=self.problems,
            label_source=self.label_source,
            times=self.times,
        )


class VasprunReader:
    """Stream one vasprun.xml: the header is read on construction, steps lazily.

    >>> with VasprunReader(path) as reader:          # doctest: +SKIP
    ...     for step in reader.steps():
    ...         ...
    ...     trailer = reader.trailer()

    Raises :class:`VasprunParseError` when the header (through ``atominfo``)
    cannot be read. File-level truncation is never raised: it is reported on
    the trailer and on the incomplete step.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self._source = HashedSource(self.path)
        self._events = ET.iterparse(self._source.stream, events=("start", "end"))
        self._stack: list[ET.Element] = []
        self._root: ET.Element | None = None
        self._error: BaseException | None = None
        self._pending_start: ET.Element | None = None
        self._steps_started = False
        self._finished = False
        self._trailer = VasprunTrailer(
            closed=False, final_structure_present=False, truncated=False, error=None, steps_seen=0,
            compression=self._source.compression,
        )
        try:
            self.header = self._read_header()
        except BaseException:
            self.close()
            raise

    # -- context manager -------------------------------------------------
    def __enter__(self) -> "VasprunReader":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def close(self) -> None:
        self._source.close()

    # -- events ------------------------------------------------------------
    def _next(self) -> tuple[str, ET.Element] | None:
        if self._error is not None:
            return None
        try:
            return next(self._events)
        except StopIteration:
            return None
        except _READ_ERRORS as exc:
            self._error = exc
            return None

    def _depth(self) -> int:
        return len(self._stack) - 1

    def _detach(self, elem: ET.Element) -> None:
        """Remove a finished element from its (still open) parent."""
        if self._stack:
            parent = self._stack[-1]
            if len(parent) and parent[-1] is elem:
                del parent[-1]
            else:  # pragma: no cover - defensive; a finished element is always the last child
                parent.remove(elem)

    # -- header ------------------------------------------------------------
    def _read_header(self) -> VasprunHeader:
        problems: list[str] = []
        duplicates: list[str] = []
        sections: dict[str, Any] = {}
        while True:
            item = self._next()
            if item is None:
                break
            event, elem = item
            if event == "start":
                self._stack.append(elem)
                depth = self._depth()
                if depth == 0:
                    self._root = elem
                    if elem.tag != "modeling":
                        raise VasprunParseError(f"{self.path}: root element is <{elem.tag}>, not <modeling>")
                elif depth == 1 and self._starts_steps(elem):
                    self._pending_start = elem
                    break
                continue
            self._stack.pop()
            depth = len(self._stack)
            if depth == 0:
                self._trailer.closed = True
                break
            if depth != 1:
                continue
            tag = elem.tag
            if tag == "generator":
                sections["generator"] = _generator(elem)
            elif tag == "incar":
                sections["incar"] = _flat_items(elem, problems, "incar")
            elif tag == "kpoints":
                sections["kpoints"] = _kpoints(elem, problems)
            elif tag == "parameters":
                sections["parameters"] = parse_parameters(elem, problems, duplicates)
            elif tag == "atominfo":
                sections["atominfo"] = _atominfo(elem, problems)
            elif tag == "structure" and elem.get("name") == "initialpos":
                sections["initialpos"] = _structure(elem, problems, "initialpos")
            self._detach(elem)
        if "atominfo" not in sections:
            detail = f" ({self._error})" if self._error is not None else ""
            raise VasprunParseError(f"{self.path}: vasprun.xml header is incomplete: no <atominfo>{detail}")
        for required in ("generator", "incar", "kpoints", "parameters", "initialpos"):
            if required not in sections:
                problems.append(f"missing: <{required}> not found before the first ionic step")
        generator = sections.get("generator", {})
        species, atom_types = sections["atominfo"]
        initial = sections.get("initialpos")
        if initial is not None and initial.fractional is not None and len(initial.fractional) != len(species):
            problems.append(
                f"bad_shape: initialpos has {len(initial.fractional)} positions for {len(species)} atoms"
            )
        return VasprunHeader(
            generator=generator,
            incar=sections.get("incar", {}),
            parameters=sections.get("parameters", {}),
            kpoints=sections.get("kpoints", {}),
            species=species,
            atom_types=atom_types,
            initial_structure=initial,
            problems=problems,
            version_tuple=parse_version(generator.get("version")),
            parameter_duplicates=duplicates,
        )

    @staticmethod
    def _starts_steps(elem: ET.Element) -> bool:
        if elem.tag == "calculation":
            return True
        if elem.tag == "structure":
            return elem.get("name") not in _HEADER_STRUCTURES
        return elem.tag in {"varray", "energy"}

    # -- steps -------------------------------------------------------------
    def steps(self) -> Iterator[IonicStep]:
        """Yield every ionic step in file order (single use)."""
        if self._steps_started:
            raise RuntimeError("VasprunReader.steps() can only be iterated once")
        self._steps_started = True
        return self._iter_steps()

    def _iter_steps(self) -> Iterator[IonicStep]:
        trailer = self._trailer
        index = 0
        dft: _StepBuilder | None = None  # open <calculation>
        ml: _StepBuilder | None = None  # open flat MLFF step
        pending = self._pending_start
        self._pending_start = None
        try:
            while True:
                if pending is not None:  # start event consumed by the header reader
                    item: tuple[str, ET.Element] | None = ("start", pending)
                    pending = None
                    already_on_stack = True
                else:
                    item = self._next()
                    already_on_stack = False
                if item is None:
                    break
                event, elem = item
                if event == "start":
                    if not already_on_stack:
                        self._stack.append(elem)
                    if self._depth() != 1:
                        continue
                    tag = elem.tag
                    # A flat ML step is structure -> varray(s) -> energy, closed at its
                    # <energy>; only its <time> may follow. Anything else ends it.
                    if ml is not None:
                        continues = tag == "time" if ml.energy_seen else tag in {"varray", "energy"}
                        if not continues:
                            if not ml.energy_seen:
                                ml.problems.append(
                                    f"truncated_frame: flat MLFF {ml.what} ended without <energy>"
                                )
                            yield self._emit(ml, complete=ml.energy_seen)
                            ml = None
                    if tag == "calculation":
                        dft = _StepBuilder(index, "dft")
                        index += 1
                    elif tag == "structure" and elem.get("name") not in _HEADER_STRUCTURES | {"finalpos"}:
                        ml = _StepBuilder(index, "mlff")
                        index += 1
                    elif tag in {"varray", "energy"} and ml is None:
                        ml = _StepBuilder(index, "mlff")
                        ml.problems.append(f"bad_shape: flat MLFF {ml.what} has no <structure>")
                        index += 1
                    continue
                # end event
                self._stack.pop()
                depth = len(self._stack)
                if depth == 0:
                    trailer.closed = True
                    continue
                if depth == 1:
                    tag = elem.tag
                    if tag == "calculation" and dft is not None:
                        yield self._emit(dft, complete=True)
                        dft = None
                    elif tag == "structure" and elem.get("name") == "finalpos":
                        trailer.final_structure = _structure(elem, trailer.problems, "finalpos")
                        trailer.final_structure_present = True
                    elif ml is not None and tag in {"structure", "varray", "energy", "time"}:
                        ml.take(elem)
                        if tag == "time":
                            yield self._emit(ml, complete=True)
                            ml = None
                    self._detach(elem)
                    continue
                if depth == 2 and dft is not None and self._stack[1].tag == "calculation":
                    if elem.tag in _CALC_KEPT:
                        dft.take(elem)
                    self._detach(elem)
                    continue
                # depth >= 2: keep subtrees still needed, prune everything else
                top = self._stack[1].tag
                if top == "calculation":
                    if self._stack[2].tag not in _CALC_KEPT:
                        self._detach(elem)
                elif top not in _TOP_KEPT:
                    self._detach(elem)
            # iteration stopped: end of file or read error
            if self._error is not None:
                trailer.error = f"{type(self._error).__name__}: {self._error}"
                if trailer.closed:
                    trailer.problems.append(f"info: data after </modeling> ignored ({trailer.error})")
                else:
                    trailer.truncated = True
                    if len(self._stack) > 1:
                        trailer.partial_tail_tag = self._stack[1].tag
            if dft is not None:
                dft.problems.append(
                    f"truncated_frame: file ends inside <calculation> ({trailer.error or 'unexpected end'})"
                )
                yield self._emit(dft, complete=False)
            if ml is not None:
                if not ml.energy_seen:
                    ml.problems.append(
                        f"truncated_frame: file ends inside flat MLFF {ml.what} ({trailer.error or 'unexpected end'})"
                    )
                yield self._emit(ml, complete=ml.energy_seen)
        finally:
            self._finish()

    def _emit(self, builder: _StepBuilder, *, complete: bool) -> IonicStep:
        step = builder.build(complete)
        self._trailer.steps_seen += 1
        if complete:
            self._trailer.n_complete_steps += 1
        return step

    def _finish(self) -> None:
        if self._finished:
            return
        self._finished = True
        try:
            self._source.finish()
        except _READ_ERRORS:  # pragma: no cover - the raw file vanished mid-read
            pass
        self._trailer.source_sha256 = self._source.sha256
        self._trailer.source_bytes = self._source.nbytes
        self.close()

    def trailer(self) -> VasprunTrailer:
        """The trailer; valid only after :meth:`steps` has been exhausted."""
        if not self._finished:
            raise RuntimeError("read every ionic step before asking for the vasprun trailer")
        return self._trailer


def iter_vasprun(
    path: Path,
) -> tuple[VasprunHeader, Iterator[IonicStep], Callable[[], VasprunTrailer]]:
    """``(header, steps, trailer)`` for one vasprun.xml[.gz|.bz2|.xz].

    ``steps`` yields :class:`IonicStep` in file order (DFT and MLFF);
    ``trailer()`` may be called once ``steps`` is exhausted. The file is closed
    when iteration ends (or when the iterator is garbage collected).
    """
    reader = VasprunReader(path)
    return reader.header, reader.steps(), reader.trailer


def read_vasprun(path: Path) -> tuple[VasprunHeader, list[IonicStep], VasprunTrailer]:
    """Eager convenience wrapper (tests, small files)."""
    header, steps, trailer = iter_vasprun(path)
    collected = list(steps)
    return header, collected, trailer()


# --------------------------------------------------------------------------
# Derived quantities (amendments item 2 / research S1.7, S3.2, S2.1)
# --------------------------------------------------------------------------

#: Last VASP release in which the calculation-level ``e_wo_entrp``/``e_0_energy`` were
#: observed mislabelled (research S1.7: 5.2.2, 5.2.11, 5.2.12, 5.3.5, 5.4.1, 5.4.4, 6.0.8 --
#: ``e_wo_entrp`` holds E(sigma->0) and ``e_0_energy`` holds F - E_wo).
LAST_MISLABELLED_VERSION = (6, 0, 8)
#: First release verified consistent on real output (6.1.1, 6.1.2, 6.2.1, 6.3.0, 6.3.2, 6.4.1, 6.4.2).
#: Releases in between (6.0.9 .. 6.1.0) have not been observed and are treated as unverified.
FIRST_VERIFIED_CONSISTENT_VERSION = (6, 1, 1)
#: Maximum |direct calc-level - reconstructed| energy (eV) accepted for VASP >= 6.1.1. Both
#: sides are printed with 8 decimals, so honest agreement is ~2e-8 eV.
CALC_LEVEL_CROSSCHECK_TOL_EV = 1e-6

#: ``energy_version_gate`` values
GATE_MISLABELLED = "vasp<=6.0.8"
GATE_VERIFIED = "vasp>=6.1.1"
GATE_UNVERIFIED = "version_unverified"
GATE_UNKNOWN = "version_unknown"

#: ``energy_rule`` values (the rule used for the frame's energy set; F itself is always the
#: calculation-level ``e_fr_energy`` minus PSTRESS*V, which no VASP version mislabels)
RULE_CALC_LEVEL_DIRECT = "calc_level_direct"
RULE_RECONSTRUCTED_PREFIX = "reconstructed_last_scstep:"
RULE_RECONSTRUCTED_MISMATCH = RULE_RECONSTRUCTED_PREFIX + "calc_level_mismatch"
RULE_RECONSTRUCTED_MISSING = RULE_RECONSTRUCTED_PREFIX + "calc_level_missing"
RULE_FREE_ENERGY_ONLY = "free_energy_only:no_scstep_energies"
RULE_MLFF = "mlff_flat_step:not_a_dft_label"
RULE_UNAVAILABLE = "unavailable"

SOURCE_FREE_ENERGY = "vasprun:calculation.e_fr_energy-PSTRESS*V"
SOURCE_MLFF_FREE_ENERGY = "vasprun:mlff_flat_step.energy.e_fr_energy-PSTRESS*V"
_SOURCE_RECONSTRUCTED = {
    "energy_sigma0": "vasprun:(calculation.e_fr_energy-PSTRESS*V)+(scstep[-1].e_0_energy-scstep[-1].e_fr_energy)",
    "energy_no_entropy": "vasprun:(calculation.e_fr_energy-PSTRESS*V)+(scstep[-1].e_wo_entrp-scstep[-1].e_fr_energy)",
}
_SOURCE_DIRECT = {
    "energy_sigma0": "vasprun:calculation.e_0_energy-PSTRESS*V",
    "energy_no_entropy": "vasprun:calculation.e_wo_entrp-PSTRESS*V",
}
_CALC_NAME = {"energy_sigma0": "e_0_energy", "energy_no_entropy": "e_wo_entrp"}
ENERGY_QUANTITIES = ("free_energy", "energy_sigma0", "energy_no_entropy")


def energy_version_gate(version_tuple: tuple[int, int, int] | None) -> str:
    """Which calc-level energy rule a VASP version falls under (see :data:`LAST_MISLABELLED_VERSION`)."""
    if version_tuple is None:
        return GATE_UNKNOWN
    version = tuple(version_tuple)[:3]
    if version <= LAST_MISLABELLED_VERSION:
        return GATE_MISLABELLED
    if version >= FIRST_VERIFIED_CONSISTENT_VERSION:
        return GATE_VERIFIED
    return GATE_UNVERIFIED


def _cell_volume(step: IonicStep) -> float | None:
    structure = step.structure
    if structure is not None and structure.cell is not None and np.shape(structure.cell) == (3, 3):
        volume = abs(float(np.linalg.det(np.asarray(structure.cell, dtype=np.float64))))
        if math.isfinite(volume):
            return volume
    if structure is not None and structure.volume is not None and math.isfinite(structure.volume):
        return float(structure.volume)
    return None


def _close(a: float | None, b: float | None, tol: float) -> bool:
    return a is not None and b is not None and math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol


def derived_energies(
    step: IonicStep,
    pstress_kbar: float | None,
    version: tuple[int, int, int] | None = None,
) -> dict[str, Any]:
    """Force-consistent energies of one step plus how each was obtained.

    Values (eV; None when an input is absent, NaN inputs propagate as NaN):

    * ``free_energy`` F = calculation-level ``e_fr_energy`` - PSTRESS*V. This
      is the force-consistent label; VASP adds PV to the calculation-level
      value (see :data:`VASP_EVTOJ`) and no version mislabels ``e_fr_energy``.
    * ``energy_sigma0`` E0 and ``energy_no_entropy`` E_wo. Reconstructed from
      the final SCF step: E0 = F + (scstep[-1].e_0_energy - scstep[-1].e_fr_energy),
      E_wo = F + (scstep[-1].e_wo_entrp - scstep[-1].e_fr_energy) (the scstep
      differences exclude PV and Edisp, which cancel). The calculation-level
      ``e_0_energy``/``e_wo_entrp`` are used directly ONLY for versions verified
      unaffected by the VASP <= 6.0.8 mislabel (>= 6.1.1) and only when they
      agree with the reconstruction within :data:`CALC_LEVEL_CROSSCHECK_TOL_EV`;
      otherwise the reconstruction is used and ``calc_level_energy_mismatch``
      is flagged. Legacy (<= 6.0.8), unverified (6.0.9-6.1.0) and unknown
      versions always use the reconstruction.
    * ``additive_correction`` = calc e_fr - PV - scstep[-1].e_fr (Edisp for IVDW
      runs; ~0 otherwise), ``pv_term`` = PSTRESS*V, ``smearing_entropy_term`` =
      F - E_wo (= -TS), ``volume`` = |det(cell)| of the step's own structure.

    Provenance (recorded for every frame):

    * ``energy_source`` -- source expression of F (the default label);
      ``energy_sources``/``energy_rules`` -- the same per quantity;
    * ``energy_rule`` -- the rule applied to the frame's energy set:
      ``calc_level_direct`` | ``reconstructed_last_scstep:vasp<=6.0.8`` |
      ``reconstructed_last_scstep:version_unknown`` |
      ``reconstructed_last_scstep:version_unverified`` |
      ``reconstructed_last_scstep:calc_level_mismatch`` |
      ``reconstructed_last_scstep:calc_level_missing`` |
      ``free_energy_only:no_scstep_energies`` | ``mlff_flat_step:not_a_dft_label`` |
      ``unavailable``;
    * ``energy_version_gate``, ``vasp_version`` (as given);
    * ``calc_level_deltas`` -- {quantity: direct - reconstructed} where computable;
    * ``calc_level_pattern`` -- ``consistent`` | ``mislabelled`` | ``ambiguous`` |
      ``other`` | None: which layout the calc-level values actually follow;
    * ``energy_flags`` -- e.g. ``calc_level_energy_mismatch``,
      ``calc_level_mislabel_pattern``, ``pstress_unknown``, ``volume_unknown``,
      ``missing_e_fr_energy``, ``scstep_energies_missing``, ``non_finite_energy``;
    * ``parser``/``parser_version``.
    """
    gate = energy_version_gate(version)
    is_mlff = step.label_source != "dft"
    result: dict[str, Any] = {
        "free_energy": None, "energy_sigma0": None, "energy_no_entropy": None,
        "additive_correction": None, "pv_term": None, "smearing_entropy_term": None, "volume": None,
        "energy_source": SOURCE_MLFF_FREE_ENERGY if is_mlff else SOURCE_FREE_ENERGY,
        "energy_rule": RULE_UNAVAILABLE,
        "energy_sources": {name: None for name in ENERGY_QUANTITIES},
        "energy_rules": {name: RULE_UNAVAILABLE for name in ENERGY_QUANTITIES},
        "energy_version_gate": gate,
        "vasp_version": ".".join(str(part) for part in version) if version is not None else None,
        "calc_level_deltas": {},
        "calc_level_pattern": None,
        "energy_flags": [],
        "parser": PARSER_NAME,
        "parser_version": PARSER_VERSION,
    }
    flags: list[str] = result["energy_flags"]
    volume = _cell_volume(step)
    result["volume"] = volume
    e_fr = step.energies.get("e_fr_energy")
    if e_fr is None:
        flags.append("missing_e_fr_energy")
    if pstress_kbar is None:
        flags.append("pstress_unknown")
        return result
    pstress = float(pstress_kbar)
    if pstress == 0.0:
        pv: float | None = 0.0
    elif volume is not None:
        pv = pstress * PV_EV_PER_KBAR_A3 * volume
    else:
        pv = None
        flags.append("volume_unknown")
    result["pv_term"] = pv
    if e_fr is None or pv is None:
        return result
    free = e_fr - pv
    result["free_energy"] = free
    result["energy_sources"]["free_energy"] = result["energy_source"]
    if not math.isfinite(free):
        flags.append("non_finite_energy")
    if is_mlff:
        # VASP-MLFF prediction: never a DFT label; e_fr == e_wo == e_0 by construction.
        result["energy_rule"] = RULE_MLFF
        result["energy_rules"] = {name: RULE_MLFF for name in ENERGY_QUANTITIES}
        return result
    result["energy_rules"]["free_energy"] = RULE_CALC_LEVEL_DIRECT

    last = step.scf_energies[-1] if step.scf_energies else {}
    last_fr = last.get("e_fr_energy")
    reconstructed: dict[str, float | None] = {"energy_sigma0": None, "energy_no_entropy": None}
    if last_fr is not None:
        result["additive_correction"] = e_fr - pv - last_fr
        if "e_0_energy" in last:
            reconstructed["energy_sigma0"] = free + (last["e_0_energy"] - last_fr)
        if "e_wo_entrp" in last:
            reconstructed["energy_no_entropy"] = free + (last["e_wo_entrp"] - last_fr)
    if all(value is None for value in reconstructed.values()):
        flags.append("scstep_energies_missing")
        result["energy_rule"] = RULE_FREE_ENERGY_ONLY
        for name in ("energy_sigma0", "energy_no_entropy"):
            result["energy_rules"][name] = RULE_FREE_ENERGY_ONLY
        return result

    direct: dict[str, float | None] = {}
    for name, calc_name in _CALC_NAME.items():
        value = step.energies.get(calc_name)
        direct[name] = value - pv if value is not None else None
        if direct[name] is not None and reconstructed[name] is not None:
            result["calc_level_deltas"][name] = direct[name] - reconstructed[name]
    result["calc_level_pattern"] = _calc_level_pattern(step.energies, pv, free, reconstructed)

    if gate == GATE_VERIFIED:
        available = [name for name in _CALC_NAME if reconstructed[name] is not None]
        missing_direct = [name for name in available if direct[name] is None]
        agree = not missing_direct and all(
            _close(direct[name], reconstructed[name], CALC_LEVEL_CROSSCHECK_TOL_EV) for name in available
        )
        if missing_direct:
            rule = RULE_RECONSTRUCTED_MISSING
            flags.append("calc_level_energy_missing")
        elif agree:
            rule = RULE_CALC_LEVEL_DIRECT
        else:
            rule = RULE_RECONSTRUCTED_MISMATCH
            flags.append("calc_level_energy_mismatch")
        if result["calc_level_pattern"] == "mislabelled":
            flags.append("calc_level_mislabel_pattern")
    else:
        agree = False
        rule = RULE_RECONSTRUCTED_PREFIX + gate
        if result["calc_level_pattern"] == "mislabelled":
            flags.append("calc_level_mislabel_pattern")  # informational: the known legacy layout
    result["energy_rule"] = rule
    for name in _CALC_NAME:
        if reconstructed[name] is None:
            result["energy_rules"][name] = RULE_UNAVAILABLE
            continue
        use_direct = agree and direct[name] is not None
        result[name] = direct[name] if use_direct else reconstructed[name]
        result["energy_sources"][name] = _SOURCE_DIRECT[name] if use_direct else _SOURCE_RECONSTRUCTED[name]
        result["energy_rules"][name] = rule
    if result["energy_no_entropy"] is not None:
        result["smearing_entropy_term"] = free - result["energy_no_entropy"]
    for name in ("energy_sigma0", "energy_no_entropy"):
        value = result[name]
        if value is not None and not math.isfinite(value) and "non_finite_energy" not in flags:
            flags.append("non_finite_energy")
    return result


def _calc_level_pattern(
    calc: dict[str, float], pv: float, free: float, reconstructed: dict[str, float | None],
) -> str | None:
    """``consistent`` / ``mislabelled`` / ``ambiguous`` / ``other`` for the calc-level e_wo/e_0 layout.

    consistent: e_wo - PV == E_wo and e_0 - PV == E0; mislabelled (VASP <= 6.0.8):
    e_wo - PV == E0 and e_0 == F - E_wo. None when the inputs are incomplete.
    """
    e_wo, e_0 = calc.get("e_wo_entrp"), calc.get("e_0_energy")
    e0_r, ewo_r = reconstructed.get("energy_sigma0"), reconstructed.get("energy_no_entropy")
    if e_wo is None or e_0 is None or e0_r is None or ewo_r is None:
        return None
    tol = CALC_LEVEL_CROSSCHECK_TOL_EV
    consistent = _close(e_wo - pv, ewo_r, tol) and _close(e_0 - pv, e0_r, tol)
    mislabelled = _close(e_wo - pv, e0_r, tol) and _close(e_0, free - ewo_r, tol)
    if consistent and mislabelled:
        return "ambiguous"
    if consistent:
        return "consistent"
    if mislabelled:
        return "mislabelled"
    return "other"


def label_energy(derived: dict[str, Any], quantity: str = "free_energy") -> tuple[float | None, str | None, str]:
    """``(value, energy_source, energy_rule)`` of one quantity from :func:`derived_energies`."""
    if quantity not in ENERGY_QUANTITIES:
        raise ValueError(f"unknown energy quantity {quantity!r}; expected one of {ENERGY_QUANTITIES}")
    return derived.get(quantity), derived["energy_sources"].get(quantity), derived["energy_rules"].get(quantity)


def scf_summary(step: IonicStep) -> dict[str, Any]:
    """``{n_steps, last_dE}``: electronic iterations and |e_fr[n-1] - e_fr[n-2]| (None if n < 2)."""
    n_steps = len(step.scf_energies)
    last_dE = None
    if n_steps >= 2:
        a = step.scf_energies[-1].get("e_fr_energy")
        b = step.scf_energies[-2].get("e_fr_energy")
        if a is not None and b is not None:
            last_dE = abs(a - b)
    return {"n_steps": n_steps, "last_dE": last_dE}


#: Selective-dynamics source precedence (amendments item 5; research S1.5). VASP 6.x writes the
#: flags into both initialpos and finalpos, VASP 5.2.2 only into finalpos; POSCAR/CONTCAR are
#: fallbacks for flags only (never for geometry).
SELECTIVE_PRECEDENCE = ("initialpos", "finalpos", "POSCAR", "CONTCAR")


def _flags_of(source: Any) -> tuple[bool, Any]:
    """``(provided, flags)`` for a PoscarEvidence/StepStructure-like object, a raw array, or None."""
    if source is None:
        return False, None
    if hasattr(source, "selective"):
        return True, source.selective
    return True, np.asarray(source, dtype=bool)


def resolve_selective(
    header: VasprunHeader,
    trailer: VasprunTrailer | None = None,
    *,
    poscar: Any = None,
    contcar: Any = None,
) -> dict[str, Any]:
    """Selective-dynamics flags with precedence initialpos > finalpos > POSCAR > CONTCAR.

    ``poscar``/``contcar`` may be :class:`~nio_md_prep.dataset.model.PoscarEvidence`
    objects (``.selective`` None = file has no ``Selective dynamics`` line), raw
    (N, 3) bool arrays, or None (file not available). Flags are returned
    exactly as VASP stores them: DIRECT (fractional) basis, True = the
    coordinate may move; they are never converted to a Cartesian mask and never
    become ASE constraints (raw forces stay untouched).

    Only sources that carry flags are compared: a structure without a
    ``selective`` varray is silent (VASP 5.2.2 omits it from initialpos), but it
    is listed in ``sources_without_flags`` so a reviewer can see it. Any two
    flag-carrying sources that differ (values or shape), or flags whose row count
    differs from the vasprun atom count, are listed in ``conflicts`` -- the
    policy layer turns a conflict into ``quarantined/evidence_mismatch``.

    Returns ``{flags, source, basis, sources_with_flags, sources_without_flags,
    conflicts, n_fixed_atoms, n_partially_fixed_atoms, n_fixed_components}``.
    """
    candidates: list[tuple[str, bool, Any]] = [
        ("initialpos", header.initial_structure is not None,
         header.initial_structure.selective if header.initial_structure is not None else None),
        ("finalpos", trailer is not None and trailer.final_structure is not None,
         trailer.final_structure.selective if trailer is not None and trailer.final_structure is not None else None),
    ]
    for name, obj in (("POSCAR", poscar), ("CONTCAR", contcar)):
        provided, flags = _flags_of(obj)
        candidates.append((name, provided, flags))
    with_flags = [(name, np.asarray(flags, dtype=bool)) for name, provided, flags in candidates if provided and flags is not None]
    without = [name for name, provided, flags in candidates if provided and flags is None]
    conflicts: list[dict[str, Any]] = []
    n_atoms = len(header.species)
    for name, flags in with_flags:
        if flags.ndim != 2 or flags.shape != (n_atoms, 3):
            conflicts.append({"a": name, "b": "atominfo", "kind": "shape_differs",
                              "detail": f"{name} flags have shape {flags.shape} for {n_atoms} atoms"})
    for position, (name_a, flags_a) in enumerate(with_flags):
        for name_b, flags_b in with_flags[position + 1:]:
            if flags_a.shape != flags_b.shape:
                conflicts.append({"a": name_a, "b": name_b, "kind": "shape_differs",
                                  "detail": f"{flags_a.shape} vs {flags_b.shape}"})
            elif not np.array_equal(flags_a, flags_b):
                rows = sorted({int(i) for i in np.nonzero(np.any(flags_a != flags_b, axis=1))[0]})
                conflicts.append({"a": name_a, "b": name_b, "kind": "flags_differ",
                                  "detail": f"{len(rows)} atom(s) differ, first 0-based indices {rows[:10]}"})
    chosen_name, chosen = (with_flags[0] if with_flags else (None, None))
    summary = {"n_fixed_atoms": 0, "n_partially_fixed_atoms": 0, "n_fixed_components": 0}
    if chosen is not None and chosen.ndim == 2 and chosen.shape[1:] == (3,):
        frozen = ~chosen
        per_atom = frozen.sum(axis=1)
        summary = {
            "n_fixed_atoms": int(np.count_nonzero(per_atom == 3)),
            "n_partially_fixed_atoms": int(np.count_nonzero((per_atom > 0) & (per_atom < 3))),
            "n_fixed_components": int(frozen.sum()),
        }
    return {
        "flags": chosen,
        "source": chosen_name,
        "basis": "direct",
        "sources_with_flags": [name for name, _ in with_flags],
        "sources_without_flags": without,
        "conflicts": conflicts,
        **summary,
    }


def selective_flags(header: VasprunHeader, trailer: VasprunTrailer | None) -> tuple[Any, str | None, bool]:
    """Backwards-compatible vasprun-only view of :func:`resolve_selective`: ``(flags, source, conflict)``."""
    resolved = resolve_selective(header, trailer)
    return resolved["flags"], resolved["source"], bool(resolved["conflicts"])


# --------------------------------------------------------------------------
# VASP-MLFF settings (amendments item 1; research S8)
# --------------------------------------------------------------------------

#: ML_MODE values under which VASP performs no new ab initio calculations.
MLFF_MODES_WITHOUT_DFT = frozenset({"run", "select", "refit", "refitbayesian"})
#: Legacy ML_ISTART -> ML_MODE (research S8: 0/1 train, 3 select, 4 refit; 2 = prediction only
#: per the VASP wiki ML_ISTART page [INFERRED: not seen in the research sample]).
_ML_ISTART_MODE = {0: "train", 1: "train", 2: "run", 3: "select", 4: "refit"}


def mlff_settings(header: VasprunHeader) -> dict[str, Any]:
    """VASP-MLFF state from the vasprun ``<incar>`` (ML_* tags are not in ``<parameters>``).

    ``active`` is True/False, or None when the ``<incar>`` section was not read
    (then MLFF use cannot be excluded). ``mode`` is the effective ML_MODE
    (legacy ML_ISTART mapped). ``dft_labels``: ``"all"`` (MLFF off), ``"per_step"``
    (train: DFT ``<calculation>`` blocks interleaved with force-field flat steps),
    ``"none"`` (run/select/refit*), ``"not_pure_dft"`` (delta), ``"unknown"``.
    """
    incar = {str(key).upper(): value for key, value in (header.incar or {}).items()}
    incar_read = not any(problem.startswith("missing: <incar>") for problem in header.problems)
    raw_flag = incar.get("ML_LMLFF")
    active = raw_flag if isinstance(raw_flag, bool) else (_logical(str(raw_flag)) if raw_flag is not None else None)
    if active is None and raw_flag is None:
        active = False if incar_read else None
    mode = incar.get("ML_MODE")
    mode = str(mode).strip().lower() if mode is not None else None
    istart = incar.get("ML_ISTART")
    if mode is None and istart is not None:
        try:
            mode = _ML_ISTART_MODE.get(int(float(istart)))
        except (TypeError, ValueError):
            mode = None
    if active is False:
        dft_labels = "all"
    elif active is None:
        dft_labels = "unknown"
    elif mode in MLFF_MODES_WITHOUT_DFT:
        dft_labels = "none"
    elif mode == "delta":
        dft_labels = "not_pure_dft"
    elif mode in {None, "train"}:
        dft_labels = "per_step"
    else:
        dft_labels = "unknown"
    return {
        "active": active, "mode": mode, "ml_lmlff": raw_flag, "ml_mode": incar.get("ML_MODE"),
        "ml_istart": istart, "dft_labels": dft_labels, "source": "vasprun_incar",
    }
