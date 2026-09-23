"""Lossless, ASE-compatible extended-XYZ writing, verification and slicing.

Why this module writes the text itself
--------------------------------------
ASE 3.29 ``ase.io.write(format="extxyz")`` formats every per-atom float
column with ``%16.8f`` (``ase/io/extxyz.py`` ``output_column_format``).
VASP prints forces with 8 decimals, but Cartesian positions are computed
here as ``fractional @ cell`` and carry ~16 significant digits, so ASE's
writer would perturb them by up to 5e-9 A. ASE's *reader* parses every
column with ``float()``, so a file whose numbers are written with Python's
shortest round-trip ``repr`` is read back bit-for-bit. This module therefore
serialises :class:`ase.Atoms` to extended XYZ with ``repr`` precision,
following the exact comment-line grammar of ASE's reader
(``key_val_str_to_dict``), and proves the round trip by re-reading with
``ase.io.iread`` (:func:`read_back_and_verify`, and ``verify=True`` in
:func:`write_extxyz`).

Other reader behaviours handled explicitly (all measured on ASE 3.29):

* info/array names in ``ase.calculators.calculator.all_properties``
  (``energy``, ``forces``, ``stress``, ...) are moved into a
  ``SinglePointCalculator`` and ``move_mask`` becomes constraints ->
  label keys must avoid :data:`RESERVED_KEYS`;
* string values that look numeric or boolean (``"123"``, ``"3e12"``,
  ``"T"``, ``""``) are converted by the reader, and ``\\`` is an escape
  character -> such strings are written as ``"_JSON <json string>"``, which
  the reader decodes back to the identical ``str``;
* a value ``None`` would be read back as ``True`` -> ``None`` info values are
  omitted (absence == unknown; ``frames.jsonl`` keeps the explicit null);
* files are opened by ASE with the platform default encoding -> output is
  pure ASCII (non-ASCII text is ``\\uXXXX``-escaped inside ``_JSON``).

Selective dynamics are exported as the raw VASP flags in the DIRECT basis
(array ``vasp_selective_dynamics``, True = coordinate may move, plus info
``selective_dynamics_basis="direct"``); no ASE constraint is ever attached,
because ASE drops ``FixScaled`` and interprets ``move_mask`` itself.
Forces are always the raw VASP forces (never zeroed on fixed atoms).

Frame contract (fail closed, :class:`FrameContract`): every frame carries
the metadata in :data:`REQUIRED_INFO_KEYS` (``None`` is missing, so unknowns
must be explicit), ``stress_source`` when stress is exported or
``stress_reason`` when it is not (stress is omitted, never zero-filled), a
module-written ``label_set`` (``energy_forces``; ``energy_only`` frames have
no forces column and are allowed only with ``allow_energy_only``), and, for
clean data, ``label_source="dft"`` and ``scf_status="converged"``. The
``split`` key is reserved: it is inserted only into train/valid/test files
by :func:`slice_extxyz_by_frame_ids` (``add_info``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
import os
from pathlib import Path
import re
from typing import Any, Iterable, Iterator, Mapping, Sequence

from .errors import DatasetError, DependencyMissingError
from .fsio import AtomicFile, sha256_file

# --------------------------------------------------------------------------
# Names ASE's extended-XYZ reader treats specially
# --------------------------------------------------------------------------

#: ``ase.calculators.calculator.all_properties`` (ASE 3.29): the extxyz reader
#: moves info/array entries with these names into a SinglePointCalculator.
CALCULATOR_PROPERTY_KEYS = frozenset({
    "energy", "free_energy", "energies", "forces", "stress", "stresses",
    "dipole", "charges", "magmom", "magmoms", "polarization",
    "dielectric_tensor", "born_effective_charges",
})
#: Remaining ``ase.outputs.all_outputs`` names (not converted by extxyz today,
#: reserved so a label can never be mistaken for a calculator output).
OUTPUT_PROPERTY_KEYS = frozenset({
    "natoms", "nbands", "nkpts", "nspins", "fermi_level", "kpoint_weights",
    "ibz_kpoints", "eigenvalues", "occupations",
})
#: Names with structural meaning in the extxyz grammar or in ase.Atoms:
#: ``Lattice``/``Properties``/``pbc`` (comment line), ``virial``/``stress``
#: (SPECIAL_3_3_KEYS), ``uid`` (UNPROCESSED_KEYS), ``move_mask``
#: (constraints), the PROPERTY_NAME_MAP column names and ase.Atoms arrays.
STRUCTURAL_KEYS = frozenset({
    "lattice", "properties", "pbc", "virial", "uid", "comment", "move_mask",
    "species", "symbols", "pos", "positions", "z", "numbers", "charge",
    "initial_charges", "initial_magmoms", "masses", "momenta", "velocities", "tags",
})
#: Case-insensitive set of names that may not be used as label keys or as
#: exported info/array names.
RESERVED_KEYS = frozenset(
    name.lower() for name in CALCULATOR_PROPERTY_KEYS | OUTPUT_PROPERTY_KEYS | STRUCTURAL_KEYS
)

#: Per-atom arrays and info written by this module (not payload-configurable).
SELECTIVE_ARRAY = "vasp_selective_dynamics"
MAGMOM_INITIAL_ARRAY = "vasp_magmom_initial"
MAGMOM_FINAL_ARRAY = "vasp_magmom_final"
SELECTIVE_BASIS_KEY = "selective_dynamics_basis"
FRAME_ID_KEY = "frame_id"
#: Written only into train/valid/test files by :func:`slice_extxyz_by_frame_ids`.
SPLIT_KEY = "split"
#: ``energy_forces`` (default) or ``energy_only`` (no forces array; needs
#: ``FrameContract(allow_energy_only=True)``, i.e. ``--allow-energy-only``).
LABEL_SET_KEY = "label_set"
LABEL_SET_ENERGY_FORCES = "energy_forces"
LABEL_SET_ENERGY_ONLY = "energy_only"
MODULE_ARRAYS = (SELECTIVE_ARRAY, MAGMOM_INITIAL_ARRAY, MAGMOM_FINAL_ARRAY)
MODULE_INFO_KEYS = (FRAME_ID_KEY, SELECTIVE_BASIS_KEY, SPLIT_KEY, LABEL_SET_KEY)

#: Per-frame metadata every exported frame must carry (fail closed; ``None``
#: counts as missing, so unknown values must be spelled out, e.g. "unknown").
#: Together with the module keys (``frame_id``, ``label_set``, ``split`` in
#: split files) and the label keys this is the user-required metadata list:
#: source path and file type, ionic index, structure hash, lineage group,
#: energy derivation, stress availability, SCF status, DFT-vs-MLFF label
#: source, VASP version, magnetic class and policy result, campaign/family,
#: parser name/version and repository commit.
REQUIRED_INFO_KEYS = (
    "source",  # "<alias>:<relpath>/<file>"
    "source_file_type",  # vasprun | vasprun.gz | vasprun.bz2 | vasprun.xz
    "ionic_step",  # 0-based index among ALL ionic steps of the file (DFT and MLFF)
    "structure_key",  # exact, order-dependent structure hash
    "lineage_group",  # split unit: frames of one group never land in different splits
    "energy_source",  # e.g. "vasprun:calculation.e_fr_energy-PSTRESS*V"
    "energy_rule",  # calc_level_direct | reconstructed_last_scstep:vasp<=6.0.8 | ...
    "stress_available",  # bool; stress_source (True) or stress_reason (False) required too
    "scf_status",  # "converged" for clean data
    "label_source",  # "dft" for clean data (MLFF steps are never exported as DFT)
    "vasp_version",
    "magnetic_class",  # model.MAGNETIC_CLASSES
    "magnetic_policy",  # policy result, e.g. "accepted" or "override:<reason>"
    "campaign",
    "family",
    "parser",
    "parser_version",
    "repo_commit",
)
_REQUIRED_INT_KEYS = frozenset({"ionic_step"})
_REQUIRED_BOOL_KEYS = frozenset({"stress_available"})

#: Further provenance info keys the exporter may supply (spec "Output files");
#: informative, not enforced. ``None`` values are omitted.
PROVENANCE_INFO_KEYS = (
    "vasp_free_energy", "vasp_energy_no_entropy", "vasp_energy_sigma0", "energy_quantity",
    "run_id", "source_sha256", "time_fs", "calc_type", "pool_id",
    "group_id", "lineage_id", "config_type", "n_scf_steps", "total_magnetization",
    "stress_source", "stress_reason", "forces_reason", "vacuum_axes", "permutation_key",
    "scf_evidence", "magnetic_state_id", "magnetic_segment",
)


@dataclass(frozen=True)
class FrameContract:
    """What a frame must satisfy to be written.

    ``required_info``: info keys that must be present and not ``None``
    (default :data:`REQUIRED_INFO_KEYS`); ``allow_energy_only``: permit
    frames without forces (``label_set="energy_only"``, no forces column);
    ``clean``: refuse frames whose ``label_source`` is not ``"dft"`` or whose
    ``scf_status`` is not ``"converged"`` (MLFF predictions and failed or
    uncertain SCF are never exported as clean training data).
    """

    required_info: tuple[str, ...] = REQUIRED_INFO_KEYS
    allow_energy_only: bool = False
    clean: bool = True

    def as_dict(self) -> dict[str, Any]:
        return {"required_info": list(self.required_info), "allow_energy_only": self.allow_energy_only,
                "clean": self.clean}


#: Used when re-reading/comparing: checks shapes and finiteness only.
PERMISSIVE_CONTRACT = FrameContract(required_info=(), allow_energy_only=True, clean=False)

_KEY_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_SAFE_PLAIN_RE = re.compile(r"^[A-Za-z0-9_.:/#@+-]+$")
_BOOL_WORDS = frozenset({"T", "F", "true", "false", "True", "False", "TRUE", "FALSE"})
_JSON_PREFIX = "_JSON "


def _require_ase():
    try:
        import ase  # noqa: F401
        import ase.io  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without ASE
        raise DependencyMissingError(
            "writing/reading extended XYZ needs ASE: pip install 'nio-md-prep[dataset]'"
        ) from exc
    return ase


def _np():
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise DependencyMissingError("numpy is required for dataset export") from exc
    return np


# --------------------------------------------------------------------------
# Label keys
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class LabelKeys:
    energy: str = "REF_energy"
    forces: str = "REF_forces"
    stress: str = "REF_stress"

    def as_dict(self) -> dict[str, str]:
        return {"energy": self.energy, "forces": self.forces, "stress": self.stress}


def validate_label_keys(
    energy_key: str = "REF_energy", forces_key: str = "REF_forces", stress_key: str = "REF_stress"
) -> LabelKeys:
    """Check user label keys and return them as :class:`LabelKeys`.

    Rejected: non-identifiers (the extxyz grammar and ``Properties`` string
    need ``[A-Za-z_][A-Za-z0-9_]*``), any name in :data:`RESERVED_KEYS`
    (case-insensitive; ASE would turn it into a calculator result, constraint
    or structural field), names used by this module, and duplicates.
    """
    keys = {"energy": energy_key, "forces": forces_key, "stress": stress_key}
    problems = []
    for role, key in keys.items():
        if not isinstance(key, str) or not _KEY_RE.match(key):
            problems.append(f"{role} key {key!r} is not a valid extxyz name ([A-Za-z_][A-Za-z0-9_]*)")
        elif key.lower() in RESERVED_KEYS:
            problems.append(
                f"{role} key {key!r} is reserved: ASE's extxyz reader would convert it into a "
                "calculator result, constraint or structural field; use e.g. REF_" + role
            )
        elif key.lower() in {name.lower() for name in MODULE_ARRAYS + MODULE_INFO_KEYS + REQUIRED_INFO_KEYS}:
            problems.append(f"{role} key {key!r} collides with a field written by the exporter")
    lowered = [key.lower() for key in keys.values() if isinstance(key, str)]
    if len(set(lowered)) != len(lowered):
        problems.append(f"label keys must be distinct (case-insensitive): {keys}")
    if problems:
        raise DatasetError("invalid label keys: " + "; ".join(problems))
    return LabelKeys(energy_key, forces_key, stress_key)


# --------------------------------------------------------------------------
# Frame payload
# --------------------------------------------------------------------------

@dataclass
class FramePayload:
    """Everything written for one exported frame.

    Units: Angstrom, eV, eV/A, eV/A^3 (``stress`` already in the ASE sign
    convention, e.g. via :func:`model.vasp_stress_to_ase`). Geometry and labels
    must come from the same ionic step. ``info`` holds provenance scalars
    (see :data:`PROVENANCE_INFO_KEYS`); ``None`` values are omitted from the
    file. ``selective_dynamics`` is the raw VASP (N, 3) flag array in the
    DIRECT basis, True = may move. ``magmom_final`` belongs only on the step
    whose OUTCAR table it is. ``forces=None`` marks an energy-only frame; it
    is written only under ``FrameContract(allow_energy_only=True)`` and then
    carries ``label_set="energy_only"`` and no forces column.
    """

    frame_id: str
    species: Sequence[str]
    cell: Any
    positions: Any
    energy: float
    forces: Any | None
    stress: Any | None = None
    info: dict[str, Any] = field(default_factory=dict)
    selective_dynamics: Any | None = None
    magmom_initial: Any | None = None
    magmom_final: Any | None = None

    @classmethod
    def from_mapping(cls, value: "FramePayload | Mapping[str, Any]") -> "FramePayload":
        if isinstance(value, FramePayload):
            return value
        if not isinstance(value, Mapping):
            raise DatasetError(f"frame payload must be a FramePayload or mapping, got {type(value).__name__}")
        allowed = set(cls.__dataclass_fields__)
        unknown = sorted(set(value) - allowed)
        if unknown:
            raise DatasetError(f"unknown frame payload fields {unknown}; allowed: {sorted(allowed)}")
        missing = [name for name in ("frame_id", "species", "cell", "positions", "energy", "forces") if name not in value]
        if missing:
            raise DatasetError(f"frame payload lacks required fields {missing}")
        return cls(**dict(value))


def _finite_array(np, value, shape, what: str, frame_id: str):
    array = np.array(value, dtype=np.float64, copy=True)
    if array.shape != shape:
        raise DatasetError(f"frame {frame_id}: {what} has shape {array.shape}, expected {shape}")
    if not np.all(np.isfinite(array)):
        raise DatasetError(f"frame {frame_id}: {what} contains NaN/inf")
    return array


def _check_payload(payload: FramePayload, keys: LabelKeys, contract: "FrameContract | None" = None):
    """Validate and normalise a payload; returns plain numpy/python values."""
    np = _np()
    contract = contract or FrameContract()
    frame_id = payload.frame_id
    if not isinstance(frame_id, str) or not frame_id or not frame_id.isascii() or not frame_id.isprintable():
        raise DatasetError(f"frame_id must be a non-empty printable ASCII string, got {frame_id!r}")
    species = list(payload.species)
    if not species or not all(isinstance(item, str) and item for item in species):
        raise DatasetError(f"frame {frame_id}: species must be a non-empty list of element symbols")
    n_atoms = len(species)
    cell = _finite_array(np, payload.cell, (3, 3), "cell", frame_id)
    if abs(float(np.linalg.det(cell))) <= 0.0:
        raise DatasetError(f"frame {frame_id}: cell is singular")
    positions = _finite_array(np, payload.positions, (n_atoms, 3), "positions", frame_id)
    if payload.forces is None:
        if not contract.allow_energy_only:
            raise DatasetError(
                f"frame {frame_id}: no forces; a frame without forces never enters a force-training dataset "
                "(energy-only export needs allow_energy_only / --allow-energy-only)"
            )
        forces = None
    else:
        forces = _finite_array(np, payload.forces, (n_atoms, 3), "forces", frame_id)
    try:
        energy = float(payload.energy)
    except (TypeError, ValueError):
        raise DatasetError(f"frame {frame_id}: energy {payload.energy!r} is not a number") from None
    if isinstance(payload.energy, bool) or not math.isfinite(energy):
        raise DatasetError(f"frame {frame_id}: energy {payload.energy!r} is not a finite number")
    stress = None
    if payload.stress is not None:
        stress = _finite_array(np, payload.stress, (3, 3), "stress", frame_id)
    arrays: dict[str, Any] = {}
    if payload.selective_dynamics is not None:
        flags = np.array(payload.selective_dynamics, copy=True)
        if flags.dtype.kind != "b" or flags.shape != (n_atoms, 3):
            raise DatasetError(
                f"frame {frame_id}: selective_dynamics must be a ({n_atoms}, 3) bool array "
                f"(True = may move), got dtype {flags.dtype} shape {flags.shape}"
            )
        arrays[SELECTIVE_ARRAY] = flags
    if payload.magmom_initial is not None:
        arrays[MAGMOM_INITIAL_ARRAY] = _finite_array(np, payload.magmom_initial, (n_atoms,), "magmom_initial", frame_id)
    if payload.magmom_final is not None:
        arrays[MAGMOM_FINAL_ARRAY] = _finite_array(np, payload.magmom_final, (n_atoms,), "magmom_final", frame_id)

    label_names = {keys.energy.lower(), keys.forces.lower(), keys.stress.lower()}
    taken = label_names | {name.lower() for name in MODULE_ARRAYS + MODULE_INFO_KEYS}
    info: dict[str, Any] = {}
    for key, value in payload.info.items():
        if not isinstance(key, str) or not _KEY_RE.match(key):
            raise DatasetError(f"frame {frame_id}: info key {key!r} is not a valid extxyz name")
        if key.lower() in RESERVED_KEYS:
            raise DatasetError(f"frame {frame_id}: info key {key!r} is reserved by ASE's extxyz reader")
        if key.lower() in taken:
            raise DatasetError(f"frame {frame_id}: info key {key!r} collides with a label key or exporter field")
        if value is None:
            continue
        info[key] = _check_info_value(value, f"frame {frame_id}: info {key!r}")
    _check_contract(info, contract, frame_id, has_forces=forces is not None, has_stress=stress is not None)
    info[FRAME_ID_KEY] = frame_id
    info[LABEL_SET_KEY] = LABEL_SET_ENERGY_FORCES if forces is not None else LABEL_SET_ENERGY_ONLY
    if SELECTIVE_ARRAY in arrays:
        info[SELECTIVE_BASIS_KEY] = "direct"
    return species, cell, positions, energy, forces, stress, info, arrays


def _check_contract(info: Mapping[str, Any], contract: "FrameContract", frame_id: str, *,
                    has_forces: bool, has_stress: bool) -> None:
    """Required metadata, its types, stress/forces bookkeeping and the clean-data guard."""
    problems: list[str] = []
    required = list(contract.required_info)
    if required:
        required.append("stress_source" if has_stress else "stress_reason")
        if not has_forces:
            required.append("forces_reason")
    missing = [key for key in required if info.get(key) is None]
    if missing:
        problems.append(f"required metadata missing (None counts as missing; write 'unknown' explicitly): {missing}")
    for key in required:
        value = info.get(key)
        if value is None:
            continue
        if key in _REQUIRED_INT_KEYS:
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                problems.append(f"{key} must be a non-negative int, got {value!r}")
        elif key in _REQUIRED_BOOL_KEYS:
            if not isinstance(value, bool):
                problems.append(f"{key} must be a bool, got {value!r}")
        elif not isinstance(value, str) or not value.strip():
            problems.append(f"{key} must be a non-empty string, got {value!r}")
    available = info.get("stress_available")
    if isinstance(available, bool) and available != has_stress:
        problems.append(
            f"stress_available={available} but the stress label is {'present' if has_stress else 'absent'} "
            "(absent stress is omitted, never zero-filled)"
        )
    if contract.clean:
        if info.get("label_source") not in (None, "dft"):
            problems.append(f"label_source={info.get('label_source')!r}: only DFT labels are clean training data "
                            "(VASP-MLFF steps are never exported as DFT)")
        if info.get("scf_status") not in (None, "converged"):
            problems.append(f"scf_status={info.get('scf_status')!r}: failed or uncertain SCF is never clean data")
    if problems:
        raise DatasetError(f"frame {frame_id}: " + "; ".join(problems))


def _check_info_value(value, where: str):
    """Normalise an info value to str/bool/int/float or a JSON-able container."""
    np = _np()
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if not math.isfinite(number):
            raise DatasetError(f"{where} is not finite ({value!r})")
        return number
    if isinstance(value, str):
        return value
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return _check_homogeneous_list(list(value), where)
    if isinstance(value, Mapping):
        try:
            text = json.dumps(value, sort_keys=True, allow_nan=False)
        except (TypeError, ValueError) as exc:
            raise DatasetError(f"{where} is not JSON-serialisable: {exc}") from None
        return json.loads(text)
    raise DatasetError(f"{where} has unsupported type {type(value).__name__}")


def _check_homogeneous_list(items: list, where: str):
    """Lists must be homogeneous so ASE's ``np.array(json)`` cannot change types."""
    np = _np()
    flat: list = []

    def walk(node, depth):
        if isinstance(node, (list, tuple)):
            return [walk(child, depth + 1) for child in node]
        flat.append(node)
        return node

    nested = walk(items, 0)
    if not flat:
        return nested
    if all(isinstance(item, (bool, np.bool_)) for item in flat):
        kind = "bool"
    elif all(isinstance(item, str) for item in flat):
        kind = "str"
    elif all(isinstance(item, (int, np.integer)) and not isinstance(item, (bool, np.bool_)) for item in flat):
        kind = "int"
    elif all(isinstance(item, (int, float, np.integer, np.floating)) and not isinstance(item, (bool, np.bool_)) for item in flat):
        kind = "float"
    else:
        raise DatasetError(f"{where}: list elements must all be bool, all str, or all numbers")
    if kind == "float" and not all(math.isfinite(float(item)) for item in flat):
        raise DatasetError(f"{where}: list contains NaN/inf")
    convert = {"bool": bool, "str": str, "int": int, "float": float}[kind]

    def cast(node):
        if isinstance(node, list):
            return [cast(child) for child in node]
        return convert(node)

    result = cast(nested)
    if kind != "str":
        array = np.array(result)
        if array.dtype == object:
            raise DatasetError(f"{where}: numeric lists must be rectangular")
    return result


# --------------------------------------------------------------------------
# Payload -> ase.Atoms -> extxyz text
# --------------------------------------------------------------------------

def frame_to_atoms(frame_payload, *, keys: LabelKeys | None = None, contract: FrameContract | None = None):
    """Build an :class:`ase.Atoms` for one frame, without calculator or constraints.

    ``info``: label energy under ``keys.energy``, stress (3x3) under
    ``keys.stress`` when present, ``frame_id``, ``label_set``, every non-None
    payload info entry, and ``selective_dynamics_basis="direct"`` when flags
    exist. ``arrays``: raw forces under ``keys.forces`` (absent for
    ``label_set="energy_only"``) plus the optional ``vasp_selective_dynamics``
    (N,3 bool), ``vasp_magmom_initial`` and ``vasp_magmom_final`` (N,).
    ``pbc`` is always (True, True, True): VASP cells are 3-D periodic. The
    payload must satisfy ``contract`` (default :class:`FrameContract`).
    """
    _require_ase()
    from ase import Atoms

    keys = keys or LabelKeys()
    payload = FramePayload.from_mapping(frame_payload)
    species, cell, positions, energy, forces, stress, info, arrays = _check_payload(payload, keys, contract)
    try:
        atoms = Atoms(symbols=species, positions=positions, cell=cell, pbc=(True, True, True))
    except (KeyError, ValueError) as exc:
        raise DatasetError(f"frame {payload.frame_id}: cannot build atoms: {exc}") from None
    atoms.info[keys.energy] = energy
    if stress is not None:
        atoms.info[keys.stress] = stress
    atoms.info.update(info)
    if forces is not None:
        atoms.arrays[keys.forces] = forces
    for name, array in arrays.items():
        atoms.arrays[name] = array
    return atoms


def _repr_float(value) -> str:
    text = repr(float(value))
    if text in {"nan", "inf", "-inf"}:
        raise DatasetError(f"non-finite value {text} cannot be exported")
    return text


def _quote(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _plain_string_survives(text: str) -> bool:
    """True if ASE's reader returns ``text`` unchanged for ``key="text"``."""
    if not text or not _SAFE_PLAIN_RE.match(text):
        return False
    if text in _BOOL_WORDS or text.startswith("_JSON"):
        return False
    np = _np()
    tokens = re.findall(r"[^\s,]+", text)
    for dtype in (int, float):
        try:
            np.array(tokens, dtype=dtype)
        except (ValueError, OverflowError):
            continue
        return False
    return True


def _encode_info_value(value) -> str:
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return _repr_float(value)
    if isinstance(value, str):
        if _plain_string_survives(value):
            return _quote(value)
        return _quote(_JSON_PREFIX + json.dumps(value, ensure_ascii=True))
    np = _np()
    if isinstance(value, np.ndarray):
        value = value.tolist()
    return _quote(_JSON_PREFIX + json.dumps(value, ensure_ascii=True, allow_nan=False))


def _column_spec(name: str, array) -> tuple[str, int]:
    kind = array.dtype.kind
    if kind == "b":
        code = "L"
    elif kind == "f":
        code = "R"
    elif kind in "iu":
        code = "I"
    else:
        raise DatasetError(f"array {name!r} has unsupported dtype {array.dtype}")
    ncol = 1 if array.ndim == 1 else array.shape[1]
    return code, ncol


def atoms_to_extxyz_block(atoms, *, keys: LabelKeys | None = None) -> str:
    """Serialise one Atoms (as built by :func:`frame_to_atoms`) to extxyz text.

    Deterministic: columns ``species pos <forces> <other arrays sorted>``;
    comment ``Lattice Properties <energy> [<stress>] <info sorted> pbc``.
    Floats use Python ``repr`` (shortest string that round-trips exactly).
    """
    np = _np()
    keys = keys or LabelKeys()
    if atoms.calc is not None or atoms.constraints:
        raise DatasetError("exported Atoms must carry neither a calculator nor constraints")
    n_atoms = len(atoms)
    symbols = atoms.get_chemical_symbols()
    positions = np.asarray(atoms.positions, dtype=np.float64)
    extra = [name for name in atoms.arrays if name not in ("numbers", "positions")]
    ordered = ([keys.forces] if keys.forces in extra else []) + sorted(name for name in extra if name != keys.forces)
    columns = []
    props = ["species:S:1", "pos:R:3"]
    for name in ordered:
        if name.lower() in RESERVED_KEYS:
            raise DatasetError(f"array name {name!r} is reserved by ASE's extxyz reader")
        array = np.asarray(atoms.arrays[name])
        if array.shape[0] != n_atoms or array.ndim not in (1, 2):
            raise DatasetError(f"array {name!r} has shape {array.shape}; expected ({n_atoms},) or ({n_atoms}, k)")
        code, ncol = _column_spec(name, array)
        props.append(f"{name}:{code}:{ncol}")
        columns.append((code, array.reshape(n_atoms, ncol)))
    cell = np.asarray(atoms.cell.array, dtype=np.float64)
    parts = ['Lattice="' + " ".join(_repr_float(x) for x in cell.reshape(9)) + '"', "Properties=" + ":".join(props)]
    info = dict(atoms.info)
    for key in ("pbc", "Lattice", "Properties"):
        if key in info:
            raise DatasetError(f"info key {key!r} is structural and cannot be set")
    label_first = [key for key in (keys.energy, keys.stress) if key in info]
    for key in label_first + sorted(key for key in info if key not in label_first):
        if not _KEY_RE.match(key):
            raise DatasetError(f"info key {key!r} is not a valid extxyz name")
        parts.append(f"{key}={_encode_info_value(_check_info_value(info[key], f'info {key!r}'))}")
    parts.append('pbc="' + " ".join("T" if flag else "F" for flag in atoms.pbc) + '"')
    lines = [str(n_atoms), " ".join(parts)]
    for index in range(n_atoms):
        row = [f"{symbols[index]:<2}"]
        row.extend(f"{_repr_float(x):>22}" for x in positions[index])
        for code, array in columns:
            values = array[index]
            if code == "R":
                row.extend(f"{_repr_float(x):>22}" for x in values)
            elif code == "I":
                row.extend(f"{int(x):>8d}" for x in values)
            else:
                row.extend("T" if bool(x) else "F" for x in values)
        lines.append(" ".join(row))
    text = "\n".join(lines) + "\n"
    if not text.isascii():
        raise DatasetError("extxyz output must be ASCII (ASE reads with the platform encoding)")
    return text


# --------------------------------------------------------------------------
# Normalised frame view (shared by write-time digests and read-back checks)
# --------------------------------------------------------------------------

def _normalise_value(value):
    """What ASE's reader returns for an encoded info value, as plain Python."""
    np = _np()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, (list, tuple)):
        converted = json.loads(json.dumps(list(value), default=_json_default))
        array = np.array(converted)
        return array.tolist() if array.dtype.kind in "ifb" else converted
    return value


def _json_default(value):
    np = _np()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"not JSON serialisable: {type(value).__name__}")


def _typed(value):
    """Tag scalars with their type so int 1, float 1.0 and True never compare equal."""
    if isinstance(value, bool):
        return ["bool", value]
    if isinstance(value, int):
        return ["int", value]
    if isinstance(value, float):
        return ["float", repr(value)]
    if isinstance(value, str):
        return ["str", value]
    return ["json", json.dumps(value, sort_keys=True, default=_json_default)]


@dataclass
class _FrameView:
    frame_id: Any
    species: list[str]
    cell: Any
    pbc: tuple
    positions: Any
    energy: Any
    forces: Any
    stress: Any
    info: dict[str, Any]
    arrays: dict[str, Any]
    has_calc: bool = False
    constraints: int = 0

    def digest(self) -> str:
        np = _np()
        h = hashlib.sha256()

        def put(tag: str, payload: bytes):
            h.update(tag.encode() + b"\0" + len(payload).to_bytes(8, "little") + payload)

        def arr(value):
            return np.ascontiguousarray(value, dtype="<f8").tobytes() if value is not None else b"<none>"

        put("frame_id", json.dumps(_typed(self.frame_id)).encode())
        put("species", json.dumps(self.species).encode())
        put("cell", arr(self.cell))
        put("pbc", json.dumps([bool(x) for x in self.pbc]).encode())
        put("positions", arr(self.positions))
        put("energy", json.dumps(_typed(self.energy)).encode())
        put("forces", arr(self.forces))
        put("stress", arr(self.stress))
        put("info", json.dumps({k: _typed(v) for k, v in sorted(self.info.items())}).encode())
        for name in sorted(self.arrays):
            value = np.asarray(self.arrays[name])
            put("array:" + name + ":" + value.dtype.kind + str(value.shape), np.ascontiguousarray(value).tobytes())
        put("calc", json.dumps([self.has_calc, self.constraints]).encode())
        return h.hexdigest()


def _view_from_payload(payload, keys: LabelKeys) -> _FrameView:
    np = _np()
    payload = FramePayload.from_mapping(payload)
    species, cell, positions, energy, forces, stress, info, arrays = _check_payload(payload, keys, PERMISSIVE_CONTRACT)
    norm_info = {key: _normalise_value(value) for key, value in info.items()}
    norm_arrays = {}
    for name, value in arrays.items():
        norm_arrays[name] = value.astype(bool) if value.dtype.kind == "b" else value.astype(np.float64)
    return _FrameView(info[FRAME_ID_KEY], species, cell, (True, True, True), positions, float(energy),
                      forces, stress, {k: v for k, v in norm_info.items() if k != FRAME_ID_KEY}, norm_arrays)


def _view_from_atoms(atoms, keys: LabelKeys) -> _FrameView:
    np = _np()
    info = dict(atoms.info)
    frame_id = info.pop(FRAME_ID_KEY, None)
    energy = info.pop(keys.energy, None)
    stress = info.pop(keys.stress, None)
    arrays = {name: np.asarray(value) for name, value in atoms.arrays.items() if name not in ("numbers", "positions")}
    forces = arrays.pop(keys.forces, None)
    energy_value = _normalise_value(energy) if energy is not None else None
    stress_value = None if stress is None else np.asarray(stress, dtype=np.float64)
    norm_arrays = {name: (value.astype(bool) if value.dtype.kind == "b" else value.astype(np.float64))
                   for name, value in arrays.items()}
    return _FrameView(
        _normalise_value(frame_id) if frame_id is not None else None,
        list(atoms.get_chemical_symbols()),
        np.asarray(atoms.cell.array, dtype=np.float64),
        tuple(bool(x) for x in atoms.pbc),
        np.asarray(atoms.positions, dtype=np.float64),
        energy_value,
        None if forces is None else np.asarray(forces, dtype=np.float64),
        stress_value,
        {key: _normalise_value(value) for key, value in info.items()},
        norm_arrays,
        has_calc=atoms.calc is not None,
        constraints=len(atoms.constraints),
    )


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------

@dataclass
class ExtxyzWriteResult:
    path: Path
    frames: int
    sha256: str
    frame_ids: list[str]
    verified: bool
    label_sets: dict[str, int] = field(default_factory=dict)  # label_set -> frames
    stress_frames: int = 0

    def as_dict(self) -> dict[str, Any]:
        return {"path": self.path.as_posix(), "frames": self.frames, "sha256": self.sha256, "verified": self.verified,
                "label_sets": dict(sorted(self.label_sets.items())), "stress_frames": self.stress_frames}


def write_extxyz(
    path: Path,
    payloads: Iterable,
    *,
    keys: LabelKeys | None = None,
    verify: bool = True,
    contract: FrameContract | None = None,
) -> ExtxyzWriteResult:
    """Stream payloads to ``path`` atomically, losslessly and (by default) verified.

    Frames are written in the order given; ``frame_id`` values must be
    unique and every frame must satisfy ``contract`` (default
    :class:`FrameContract`: all :data:`REQUIRED_INFO_KEYS`, forces present,
    DFT labels with converged SCF). With ``verify=True`` the file is first
    written to a hidden staging sibling, re-read with ``ase.io.iread`` and
    compared frame by frame (digest of species, cell, pbc, positions, labels,
    info and arrays, all exact; no calculator, no constraints); only a
    verified file is moved to ``path``. Nothing is published on any failure.
    Output bytes depend only on the payloads (info order is irrelevant).
    Memory use is one frame plus one digest per frame.
    """
    _require_ase()
    keys = keys or LabelKeys()
    validate_label_keys(keys.energy, keys.forces, keys.stress)
    path = Path(path)
    staging = path.with_name(f".{path.name}.unverified-{os.getpid()}") if verify else path
    frame_ids: list[str] = []
    seen: set[str] = set()
    digests: list[str] = []
    label_sets: dict[str, int] = {}
    stress_frames = 0
    try:
        with AtomicFile(staging) as handle:
            for payload in payloads:
                payload = FramePayload.from_mapping(payload)
                if payload.frame_id in seen:
                    raise DatasetError(f"duplicate frame_id {payload.frame_id!r} in extxyz output")
                atoms = frame_to_atoms(payload, keys=keys, contract=contract)
                handle.write(atoms_to_extxyz_block(atoms, keys=keys))
                seen.add(payload.frame_id)
                frame_ids.append(payload.frame_id)
                label_set = atoms.info[LABEL_SET_KEY]
                label_sets[label_set] = label_sets.get(label_set, 0) + 1
                stress_frames += keys.stress in atoms.info
                if verify:
                    # Digest of the in-memory Atoms, normalised exactly as the
                    # ASE-read Atoms will be: equal digests <=> exact round trip.
                    digests.append(_view_from_atoms(atoms, keys).digest())
        if verify:
            _verify_digests(staging, digests, frame_ids, keys)
            os.replace(staging, path)
    finally:
        if verify and staging.exists():
            staging.unlink()
    return ExtxyzWriteResult(path, len(frame_ids), sha256_file(path), frame_ids, verify, label_sets, stress_frames)


def _verify_digests(path: Path, digests: list[str], frame_ids: list[str], keys: LabelKeys) -> None:
    import ase.io

    count = 0
    for index, atoms in enumerate(ase.io.iread(str(path), index=":", format="extxyz")):
        if index >= len(digests):
            raise DatasetError(f"round-trip check of {path.name}: file holds more frames than were written")
        view = _view_from_atoms(atoms, keys)
        if view.digest() != digests[index]:
            raise DatasetError(
                f"round-trip check of {path.name} failed at frame {index} ({frame_ids[index]!r}): "
                "the ASE-read frame differs from the payload; run read_back_and_verify for details"
            )
        count += 1
    if count != len(digests):
        raise DatasetError(f"round-trip check of {path.name}: read {count} frames, wrote {len(digests)}")


# --------------------------------------------------------------------------
# Read-back verification with field-level detail
# --------------------------------------------------------------------------

@dataclass
class RoundTripReport:
    path: str
    frames_expected: int = 0
    frames_read: int = 0
    mismatches: list[dict[str, Any]] = field(default_factory=list)
    mismatch_count: int = 0
    max_abs_diff: dict[str, float] = field(default_factory=dict)
    rtol: float = 0.0
    info_keys_checked: int = 0

    @property
    def ok(self) -> bool:
        return self.mismatch_count == 0 and self.frames_expected == self.frames_read

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path, "ok": self.ok, "frames_expected": self.frames_expected,
            "frames_read": self.frames_read, "mismatch_count": self.mismatch_count,
            "mismatches": self.mismatches, "max_abs_diff": dict(sorted(self.max_abs_diff.items())),
            "rtol": self.rtol, "info_keys_checked": self.info_keys_checked,
        }

    def raise_if_failed(self) -> "RoundTripReport":
        if not self.ok:
            first = self.mismatches[0] if self.mismatches else {}
            raise DatasetError(
                f"extxyz round trip failed for {self.path}: {self.mismatch_count} mismatches, "
                f"{self.frames_read}/{self.frames_expected} frames read; first: {first}"
            )
        return self


_MAX_REPORTED = 100


def read_back_and_verify(
    path: Path,
    expected_payloads: Iterable,
    *,
    keys: LabelKeys | None = None,
    rtol: float = 0.0,
) -> RoundTripReport:
    """Re-read ``path`` with ASE and compare it with the payloads, field by field.

    Checks, per frame and in order: ``frame_id``; species; cell; pbc (T T T);
    positions; energy, forces and stress under the label keys (exact equality
    by default; ``rtol`` > 0 allows ``|a-b| <= rtol*max(|a|,|b|)``); the
    exact set of info keys and their values/types; the set of per-atom arrays
    and their values; and that ASE attached no calculator and no constraints
    (i.e. no label was captured as a calculator result). Streams both sides.
    """
    _require_ase()
    import ase.io

    np = _np()
    keys = keys or LabelKeys()
    report = RoundTripReport(path=Path(path).as_posix(), rtol=rtol)

    def mismatch(index, frame_id, field_name, detail):
        report.mismatch_count += 1
        if len(report.mismatches) < _MAX_REPORTED:
            report.mismatches.append({"index": index, "frame_id": frame_id, "field": field_name, "detail": detail})

    def compare_array(index, frame_id, name, want, got):
        if want is None and got is None:
            return
        if want is None or got is None:
            mismatch(index, frame_id, name, "present in only one of file/payload")
            return
        want = np.asarray(want)
        got = np.asarray(got)
        if want.shape != got.shape:
            mismatch(index, frame_id, name, f"shape {got.shape} != {want.shape}")
            return
        if want.dtype.kind == "b" or got.dtype.kind == "b":
            if want.dtype.kind != got.dtype.kind or not np.array_equal(want, got):
                mismatch(index, frame_id, name, "boolean array differs")
            return
        diff = np.abs(got.astype(np.float64) - want.astype(np.float64))
        largest = float(diff.max()) if diff.size else 0.0
        report.max_abs_diff[name] = max(report.max_abs_diff.get(name, 0.0), largest)
        scale = np.maximum(np.abs(want), np.abs(got))
        if not np.all(diff <= rtol * scale):
            mismatch(index, frame_id, name, f"max |diff| {largest:.3e} exceeds rtol {rtol:g}")

    expected_iter = iter(expected_payloads)
    read_iter = ase.io.iread(str(path), index=":", format="extxyz")
    index = 0
    while True:
        want_payload = next(expected_iter, _END)
        got_atoms = next(read_iter, _END)
        if want_payload is _END and got_atoms is _END:
            break
        if want_payload is not _END:
            report.frames_expected += 1
        if got_atoms is not _END:
            report.frames_read += 1
        if want_payload is _END or got_atoms is _END:
            mismatch(index, None, "frame_count", "file and payload stream differ in length")
            index += 1
            continue
        want = _view_from_payload(want_payload, keys)
        got = _view_from_atoms(got_atoms, keys)
        fid = want.frame_id
        if got.frame_id != want.frame_id:
            mismatch(index, fid, "frame_id", f"file has {got.frame_id!r}")
        if got.species != want.species:
            mismatch(index, fid, "species", "species sequence differs")
        if tuple(got.pbc) != (True, True, True):
            mismatch(index, fid, "pbc", f"file pbc {got.pbc}")
        if got.has_calc:
            mismatch(index, fid, "calculator", "ASE attached a calculator (a key was read as a calculator result)")
        if got.constraints:
            mismatch(index, fid, "constraints", "ASE attached constraints")
        compare_array(index, fid, "cell", want.cell, got.cell)
        if got.positions.shape == want.positions.shape:
            compare_array(index, fid, "positions", want.positions, got.positions)
        else:
            mismatch(index, fid, "positions", f"shape {got.positions.shape} != {want.positions.shape}")
        if not isinstance(got.energy, float):
            mismatch(index, fid, keys.energy, f"missing or not a float in file ({got.energy!r})")
        else:
            compare_array(index, fid, keys.energy, np.array([want.energy]), np.array([got.energy]))
        compare_array(index, fid, keys.forces, want.forces, got.forces)
        compare_array(index, fid, keys.stress, want.stress, got.stress)
        if set(got.info) != set(want.info):
            mismatch(index, fid, "info_keys",
                     f"missing {sorted(set(want.info) - set(got.info))}, extra {sorted(set(got.info) - set(want.info))}")
        for key in sorted(set(got.info) & set(want.info)):
            report.info_keys_checked += 1
            if _typed(got.info[key]) != _typed(want.info[key]):
                mismatch(index, fid, f"info:{key}", f"file {got.info[key]!r} != payload {want.info[key]!r}")
        if set(got.arrays) != set(want.arrays):
            mismatch(index, fid, "arrays",
                     f"missing {sorted(set(want.arrays) - set(got.arrays))}, extra {sorted(set(got.arrays) - set(want.arrays))}")
        for name in sorted(set(got.arrays) & set(want.arrays)):
            compare_array(index, fid, name, want.arrays[name], got.arrays[name])
        index += 1
    return report


_END = object()


# --------------------------------------------------------------------------
# Label digest shared by export and audit
# --------------------------------------------------------------------------

def label_sha256(energy: float, forces, stress=None) -> str:
    """Digest of the exported labels (float64 little-endian bytes).

    ``E`` + energy, ``F`` + forces (N,3 C order), and ``S`` + stress (3x3)
    only when stress is exported. Used for ``FrameRecord.label_sha256`` and by
    the audit to compare ``frames.jsonl`` with the ASE-read file.
    """
    np = _np()
    h = hashlib.sha256()
    h.update(b"E" + np.array([float(energy)], dtype="<f8").tobytes())
    forces = np.ascontiguousarray(forces, dtype="<f8")
    h.update(b"F" + str(forces.shape).encode() + forces.tobytes())
    if stress is not None:
        stress = np.ascontiguousarray(stress, dtype="<f8")
        if stress.shape != (3, 3):
            raise DatasetError(f"stress must be 3x3 for label_sha256, got {stress.shape}")
        h.update(b"S" + stress.tobytes())
    return h.hexdigest()


# --------------------------------------------------------------------------
# Text-level access (stdlib only): frame blocks, frame ids, slicing
# --------------------------------------------------------------------------

def parse_comment_line(line: str) -> dict[str, str]:
    """Split an extxyz comment line into raw ``key -> value`` strings.

    Mirrors the tokenisation of ASE's ``key_val_str_to_dict`` (quotes,
    brackets, backslash escapes, bare keys = "T") without its type
    conversion. ``_JSON`` values are returned undecoded.
    """
    delimiters = {"'": "'", '"': '"', "{": "}", "[": "]"}
    pairs: list[list[list[str]]] = [[[]]]
    current = None
    escaped = False
    for char in line.strip():
        if escaped:
            pairs[-1][-1].append(char)
            escaped = False
        elif char == "\\":
            escaped = True
        elif current:
            if char == current:
                current = None
            else:
                pairs[-1][-1].append(char)
        elif char in delimiters:
            current = delimiters[char]
        elif char.isspace():
            if pairs == [[[]]] or pairs[-1][-1] == []:
                continue
            pairs.append([[]])
        elif char == "=":
            if pairs[-1] == [[]]:
                del pairs[-1]
            pairs[-1].append([])
        else:
            pairs[-1][-1].append(char)
    result: dict[str, str] = {}
    for pair in pairs:
        if not pair:
            continue
        if len(pair) == 1:
            key, value = "".join(pair[0]), "T"
        else:
            key, value = "".join(pair[0]), "=".join("".join(part) for part in pair[1:])
        if key:
            result[key] = value
    return result


def _decode_string_value(raw: str) -> str:
    if raw.startswith(_JSON_PREFIX):
        value = json.loads(raw[len(_JSON_PREFIX):])
        if not isinstance(value, str):
            raise DatasetError(f"expected a string value, got JSON {type(value).__name__}")
        return value
    return raw


_INT_RE = re.compile(r"^[+-]?\d+$")
_TRUE_WORDS = frozenset({"T", "True", "true", "TRUE"})
_FALSE_WORDS = frozenset({"F", "False", "false", "FALSE"})


def decode_comment_value(raw: str) -> Any:
    """Decode one raw comment-line value as written by this module (stdlib only).

    ``_JSON`` values -> their JSON value (lists stay lists); ``T``/``F`` ->
    bool; integers -> int; other numbers -> float; anything else -> str.
    Strings that would be ambiguous are always ``_JSON``-encoded by
    :func:`atoms_to_extxyz_block`, so the decoding is exact for our files.
    Whitespace-separated numbers (e.g. ``Lattice``) -> list of floats.
    """
    if raw.startswith(_JSON_PREFIX):
        return json.loads(raw[len(_JSON_PREFIX):])
    if raw in _TRUE_WORDS:
        return True
    if raw in _FALSE_WORDS:
        return False
    if _INT_RE.match(raw):
        return int(raw)
    try:
        return float(raw)
    except ValueError:
        pass
    tokens = raw.split()
    if len(tokens) > 1:
        try:
            return [float(token) for token in tokens]
        except ValueError:
            pass
    return raw


#: Comment-line fields that are not frame info.
_STRUCTURAL_FIELDS = ("Lattice", "Properties", "pbc")


@dataclass
class ExtxyzBlock:
    index: int
    frame_id: str | None
    n_atoms: int
    text: str  # the frame's exact lines, including the count and comment lines
    fields: dict[str, str] = field(default_factory=dict)  # raw comment-line key -> value strings

    def info(self) -> dict[str, Any]:
        """Decoded info of this frame (label energy/stress included; Lattice/Properties/pbc excluded)."""
        return {key: decode_comment_value(value) for key, value in self.fields.items()
                if key not in _STRUCTURAL_FIELDS}


def iter_extxyz_blocks(path: Path) -> Iterator[ExtxyzBlock]:
    """Yield frames of an extxyz file as exact text blocks (no ASE, no numpy)."""
    path = Path(path)
    with path.open("rb") as handle:
        index = 0
        while True:
            header = handle.readline()
            if not header:
                return
            if not header.strip():
                rest = handle.read()
                if rest.strip():
                    raise DatasetError(f"{path.name}: blank line inside the file before frame {index}")
                return
            try:
                n_atoms = int(header.decode("ascii").strip())
            except (UnicodeDecodeError, ValueError):
                raise DatasetError(f"{path.name}: frame {index} does not start with an atom count") from None
            if n_atoms < 0:
                raise DatasetError(f"{path.name}: frame {index} has a negative atom count")
            comment = handle.readline()
            if not comment.endswith(b"\n"):
                raise DatasetError(f"{path.name}: frame {index} is truncated")
            try:
                comment_text = comment.decode("ascii")
            except UnicodeDecodeError:
                raise DatasetError(f"{path.name}: frame {index} comment line is not ASCII") from None
            fields = parse_comment_line(comment_text)
            n_columns = _property_columns(fields.get("Properties", "species:S:1:pos:R:3"), path, index)
            lines = [header, comment]
            for _ in range(n_atoms):
                row = handle.readline()
                if not row.endswith(b"\n") or len(row.split()) != n_columns:
                    raise DatasetError(f"{path.name}: frame {index} is truncated or has a malformed atom line")
                lines.append(row)
            frame_id = _decode_string_value(fields[FRAME_ID_KEY]) if FRAME_ID_KEY in fields else None
            yield ExtxyzBlock(index, frame_id, n_atoms, b"".join(lines).decode("ascii"), fields)
            index += 1


def _property_columns(properties: str, path: Path, index: int) -> int:
    parts = properties.split(":")
    if len(parts) % 3 or not parts:
        raise DatasetError(f"{path.name}: frame {index} has a malformed Properties string")
    try:
        return sum(int(count) for count in parts[2::3])
    except ValueError:
        raise DatasetError(f"{path.name}: frame {index} has a malformed Properties string") from None


def read_frame_ids(path: Path) -> list[str | None]:
    """Frame ids of ``path`` in file order (text level)."""
    return [block.frame_id for block in iter_extxyz_blocks(path)]


_PBC_TAIL_RE = re.compile(r' pbc="[TF] [TF] [TF]"$')


def _add_info_to_block(text: str, add_info: Mapping[str, Any], fields: Mapping[str, str], where: str) -> str:
    """Insert ``key=value`` pairs into a block's comment line (before the trailing ``pbc``)."""
    count, comment, rest = text.split("\n", 2)
    pairs = []
    for key in sorted(add_info):
        if key in fields:
            raise DatasetError(f"{where}: frame already has an info key {key!r}")
        pairs.append(f"{key}={_encode_info_value(add_info[key])}")
    insertion = " " + " ".join(pairs)
    match = _PBC_TAIL_RE.search(comment)
    if match is not None:
        comment = comment[:match.start()] + insertion + comment[match.start():]
    else:
        comment = comment + insertion
    return "\n".join((count, comment, rest))


def _check_add_info(add_info: Mapping[str, Any] | None) -> dict[str, Any]:
    if not add_info:
        return {}
    checked = {}
    for key, value in add_info.items():
        if not isinstance(key, str) or not _KEY_RE.match(key) or key.lower() in RESERVED_KEYS:
            raise DatasetError(f"cannot add info key {key!r} (invalid or reserved)")
        if key in _STRUCTURAL_FIELDS or key == FRAME_ID_KEY:
            raise DatasetError(f"cannot add structural info key {key!r}")
        if value is None or isinstance(value, (list, tuple, dict, Mapping)):
            raise DatasetError(f"added info {key!r} must be a scalar str/int/float/bool, got {value!r}")
        checked[key] = _check_info_value(value, f"added info {key!r}")
    return checked


def slice_extxyz_by_frame_ids(
    src: Path,
    frame_ids: Iterable[str],
    dst: Path,
    *,
    add_info: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Copy the frames of ``src`` whose ``frame_id`` is requested into ``dst``.

    Frames are copied byte-for-byte in ``src`` order (the requested ids may be
    given in any order), except that ``add_info`` scalars (e.g.
    ``{"split": "train"}``) are inserted into each copied comment line just
    before the trailing ``pbc`` field; a frame that already has such a key is
    an error. Fails, publishing nothing, when a frame of ``src`` lacks a
    ``frame_id``, ``src`` repeats a frame_id, a requested id is repeated, or a
    requested id is absent from ``src``. ``dst`` is written atomically.
    Returns ``{path, frames, frame_ids (written order), sha256}``.
    """
    extra = _check_add_info(add_info)
    requested: set[str] = set()
    for frame_id in frame_ids:
        if frame_id in requested:
            raise DatasetError(f"frame_id {frame_id!r} requested twice for {Path(dst).name}")
        requested.add(frame_id)
    written: list[str] = []
    seen: set[str] = set()
    with AtomicFile(Path(dst)) as handle:
        for block in iter_extxyz_blocks(Path(src)):
            if block.frame_id is None:
                raise DatasetError(f"{Path(src).name}: frame {block.index} has no frame_id")
            if block.frame_id in seen:
                raise DatasetError(f"{Path(src).name}: frame_id {block.frame_id!r} occurs more than once")
            seen.add(block.frame_id)
            if block.frame_id in requested:
                text = block.text
                if extra:
                    text = _add_info_to_block(text, extra, block.fields, f"{Path(src).name} frame {block.index}")
                handle.write(text)
                written.append(block.frame_id)
        missing = requested - seen
        if missing:
            raise DatasetError(
                f"{len(missing)} requested frame_id(s) are not in {Path(src).name}: {sorted(missing)[:5]}"
            )
    return {"path": Path(dst).as_posix(), "frames": len(written), "frame_ids": written, "sha256": sha256_file(Path(dst))}
