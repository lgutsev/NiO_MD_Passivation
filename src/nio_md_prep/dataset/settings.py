"""Reference settings of a VASP run, method fingerprints and element-aware pools.

A dataset may only mix frames computed with the same electronic-structure
method. :func:`extract_settings` turns one run's evidence into three groups:

* ``method`` -- the *hard* fields that define the reference method (functional,
  cutoff, precision, spin/SOC, smearing, net charge, dipole correction, vdW,
  per-element POTCAR identity and per-element DFT+U). Runs whose method maps
  conflict never share a pool unless a reviewed override says so.
* ``sampling`` -- recorded, never gated (k-points, EDIFF, NELM, ionic/MD
  controls, parallelisation, the full ``<parameters>`` table ...).
* ``unknown`` / ``provenance`` -- which method fields could not be determined
  and where every value came from (``parameters`` > ``vasprun_incar`` >
  ``outcar`` > ``incar_file``), plus conflicts between sources.

Rules (design/phase1-dataset-design.md "Reference-settings pools" and
design/phase1-amendments.md items 3, 4, 11, 12):

* ENCUT is ``ENMAX`` in ``<parameters>``; ``*****`` overflow is None; MAGMOM
  is ignored here (it is part of the magnetic record, not of the method).
* IVDW and its parameters are absent from ``<parameters>``: they come from the
  vasprun ``<incar>``, then the OUTCAR echo, then the INCAR file. An INCAR file
  that is present but does not set IVDW means IVDW=0 (VASP default). Otherwise
  IVDW is unknown -- never guessed.
* ``GGA = --`` means "the functional named by the POTCAR LEXCH": resolved from
  the POTCAR file's LEXCH, else from the TITEL family prefix
  (PAW_PBE -> PE, PAW_GGA -> 91, PAW -> CA), else unknown.
* DFT+U per element ``{el: [L, U, J]}`` with elements that carry no +U (L < 0
  or U = J = 0, or LDAU off) omitted; ``LDAUTYPE`` recorded next to it and only
  compared between runs that both have +U elements.
* POTCAR identity per element: TITEL plus, when a POTCAR file is present and
  agrees with vasprun, the dataset-bytes sha256.
* Unknown values are the string :data:`UNKNOWN`. An unknown value equals the
  same field's unknown value in another run (same state of knowledge) but never
  a known value; a missing POTCAR dataset sha behaves the same way.
* Added beyond the spec list because they change energies and forces:
  ``LDIPOL``/``IDIPOL`` (dipole/monopole corrections) and ``LUSE_VDW`` (vdW-DF).

Overrides (TOML, :func:`load_overrides`)::

    [[equivalence]]
    field = "ISPIN"
    values = ["1", "2"]
    reason = "closed-shell molecules; spin polarisation verified to vanish"
    reviewed_by = "L. Gutsev"

map the listed values of one field onto one canonical token before
comparison; they are recorded verbatim in the manifest.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

from .errors import DatasetError
from .fsio import canonical_json, sha256_file

SETTINGS_SCHEMA = 1

#: Marker for a method value that could not be determined from any source.
UNKNOWN = "<unknown>"

#: Method fields compared as a whole between runs.
GLOBAL_METHOD_FIELDS: tuple[str, ...] = (
    "GGA", "METAGGA", "LHFCALC", "AEXX", "HFSCREEN", "IVDW", "IVDW_params", "LUSE_VDW",
    "ENCUT", "PREC", "LASPH", "ISPIN", "LSORBIT", "LNONCOLLINEAR", "NUPDOWN",
    "ISMEAR", "SIGMA", "net_charge", "LDIPOL", "IDIPOL",
)
#: Per-element method maps (compared only for elements two runs share).
ELEMENT_METHOD_FIELDS: tuple[str, ...] = ("potcar", "hubbard")
#: Compared only when both runs have at least one +U element.
CONDITIONAL_METHOD_FIELDS: tuple[str, ...] = ("LDAUTYPE",)

#: Field names an ``[[equivalence]]`` override may name.
OVERRIDABLE_FIELDS: tuple[str, ...] = GLOBAL_METHOD_FIELDS + (
    "LDAUTYPE", "potcar.titel", "potcar.sha256", "hubbard",
)

#: Sources in authority order: executed values first, the on-disk INCAR last
#: (it may have been edited after the run).
SOURCE_ORDER: tuple[str, ...] = ("parameters", "vasprun_incar", "outcar", "incar_file")

_SAMPLING_TAGS: tuple[tuple[str, str], ...] = (
    ("EDIFF", "float"), ("NELM", "int"), ("NELMIN", "int"), ("NELMDL", "int"), ("EDIFFG", "float"),
    ("IBRION", "int"), ("NSW", "int"), ("POTIM", "float"), ("ISIF", "int"), ("TEBEG", "float"),
    ("TEEND", "float"), ("MDALGO", "int"), ("SMASS", "float"), ("ALGO", "str"), ("IALGO", "int"),
    ("LREAL", "str"), ("ISYM", "int"), ("NCORE", "int"), ("NPAR", "int"), ("KPAR", "int"),
    ("ADDGRID", "bool"), ("LMAXMIX", "int"), ("LORBIT", "int"), ("NWRITE", "int"), ("ISTART", "int"),
    ("ICHARG", "int"), ("PSTRESS", "float"), ("NBLOCK", "int"), ("NELECT", "float"),
    ("ML_LMLFF", "bool"), ("ML_MODE", "str"), ("ML_ISTART", "int"),
    ("LEPSILON", "bool"), ("LCALCEPS", "bool"), ("LCHIMAG", "bool"), ("LOPTICS", "bool"),
    ("GGA_COMPAT", "bool"), ("LMETAGGA", "bool"),
)

_PREC = {"l": "low", "m": "medium", "h": "high", "n": "normal", "a": "accurate", "s": "single"}
_TITEL_FAMILY_LEXCH = {"PAW_PBE": "PE", "PAW_GGA": "91", "PAW": "CA"}
_VDW_TAG = re.compile(r"^(VDW_|LVDW)")
_REPEAT = re.compile(r"^(\d+)\*(.+)$")
_NUMBER = re.compile(r"^\s*[+-]?(\d*)(?:\.(\d*))?(?:[EeDd]([+-]?\d+))?\s*$")


# --------------------------------------------------------------------------
# Value coercion (shared with acceptance.py)
# --------------------------------------------------------------------------

def _first_token(value: str) -> str:
    text = str(value).strip().rstrip(";").strip()
    return text.split()[0] if text.split() else ""


def as_bool(value: Any) -> bool | None:
    """VASP logical: bool, 'T'/'F', '.TRUE.'/'.FALSE.', 'TRUE'/'FALSE' -> bool; else None."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        token = _first_token(value).strip(".").upper()
        if token in {"T", "TRUE"}:
            return True
        if token in {"F", "FALSE"}:
            return False
    return None


def as_float(value: Any) -> float | None:
    """Finite float from a number or a VASP/Fortran token ('1E-6', '0.1D-03'); else None."""
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        number = float(value)
    elif isinstance(value, (list, tuple)):
        return as_float(value[0]) if len(value) == 1 else None
    else:
        token = _first_token(value).replace("D", "E").replace("d", "e")
        try:
            number = float(token)
        except ValueError:
            return None
    return number if math.isfinite(number) else None


def as_int(value: Any) -> int | None:
    """Integer from an int, an integral float or a token ('2', '2.0'); else None."""
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, int):
        return value
    number = as_float(value)
    if number is None or number != int(number):
        return None
    return int(number)


def as_str(value: Any) -> str | None:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (list, tuple)):
        return " ".join(str(item) for item in value)
    text = str(value).strip()
    return text or None


def expand_list(value: Any) -> list[str] | None:
    """Tokens of a VASP list value with ``N*value`` repeats expanded."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        tokens = [str(item) for item in value]
    else:
        tokens = str(value).replace(";", " ").split()
    expanded: list[str] = []
    for token in tokens:
        repeat = _REPEAT.match(token)
        if repeat:
            expanded.extend([repeat.group(2)] * int(repeat.group(1)))
        else:
            expanded.append(token)
    return expanded


def as_float_list(value: Any) -> list[float] | None:
    tokens = expand_list(value)
    if tokens is None:
        return None
    numbers = [as_float(token) for token in tokens]
    return None if any(number is None for number in numbers) else [float(n) for n in numbers]  # type: ignore[arg-type]


def as_int_list(value: Any) -> list[int] | None:
    tokens = expand_list(value)
    if tokens is None:
        return None
    numbers = [as_int(token) for token in tokens]
    return None if any(number is None for number in numbers) else [int(n) for n in numbers]  # type: ignore[arg-type]


def canonical_float(value: float) -> float:
    """Float rounded to 10 significant digits, with -0.0 folded to 0.0."""
    rounded = float(f"{float(value):.10g}")
    return rounded + 0.0


_COERCE = {
    "bool": as_bool, "int": as_int, "float": as_float, "str": as_str,
    "float_list": as_float_list, "int_list": as_int_list,
}


def printed_tolerance(text: Any) -> float:
    """Half a unit in the last printed digit of a number token (0 for non-strings)."""
    if not isinstance(text, str):
        return 0.0
    match = _NUMBER.match(_first_token(text))
    if not match:
        return 0.0
    decimals = len(match.group(2) or "")
    exponent = int(match.group(3) or 0)
    return 0.5 * 10.0 ** (exponent - decimals)


def values_agree(a: Any, b: Any, *, raw_a: Any = None, raw_b: Any = None) -> bool:
    """Tolerant equality of two coerced values; numbers use their printed precision."""
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None:
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        tolerance = max(printed_tolerance(raw_a), printed_tolerance(raw_b))
        return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12 + tolerance * (1 + 1e-9))
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return False
        raw_a_list = expand_list(raw_a) if isinstance(raw_a, str) else None
        raw_b_list = expand_list(raw_b) if isinstance(raw_b, str) else None
        return all(
            values_agree(
                x, y,
                raw_a=raw_a_list[i] if raw_a_list and len(raw_a_list) == len(a) else None,
                raw_b=raw_b_list[i] if raw_b_list and len(raw_b_list) == len(b) else None,
            )
            for i, (x, y) in enumerate(zip(a, b))
        )
    if isinstance(a, str) and isinstance(b, str):
        return " ".join(a.upper().split()) == " ".join(b.upper().split())
    return a == b


def json_safe(value: Any) -> Any:
    """Plain JSON value: numpy -> python, tuples -> lists, non-finite floats -> None."""
    try:
        import numpy as np
    except ImportError:  # pragma: no cover - numpy is a dataset dependency
        np = None  # type: ignore[assignment]
    if np is not None:
        if isinstance(value, np.ndarray):
            return json_safe(value.tolist())
        if isinstance(value, np.generic):
            return json_safe(value.item())
    if isinstance(value, Mapping):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, (bool, int, str)):
        return value
    return str(value)


def titel_element(titel: str | None) -> str | None:
    """Element symbol of a POTCAR TITEL/pseudopotential ('PAW_PBE Ni_pv 06Sep2000' -> 'Ni')."""
    if not titel:
        return None
    tokens = titel.split()
    token = tokens[1] if len(tokens) >= 2 else tokens[0]
    token = token.split("/")[0]
    match = re.match(r"[A-Z][a-z]?", token)
    return match.group(0) if match else None


def normalize_titel(titel: Any) -> str | None:
    text = as_str(titel)
    return " ".join(text.split()) if text else None


# --------------------------------------------------------------------------
# Source lookup
# --------------------------------------------------------------------------

def _tags_of(evidence: Any) -> dict[str, Any]:
    """Upper-case tag -> value from an IncarEvidence-like object or a mapping."""
    if evidence is None:
        return {}
    tags = getattr(evidence, "tags", evidence)
    if not isinstance(tags, Mapping):
        return {}
    return {str(key).upper(): value for key, value in tags.items()}


class _Sources:
    def __init__(self, header, incar_file, outcar):
        self.tables: dict[str, dict[str, Any]] = {
            "parameters": {str(k).upper(): v for k, v in (getattr(header, "parameters", None) or {}).items()},
            "vasprun_incar": {str(k).upper(): v for k, v in (getattr(header, "incar", None) or {}).items()},
            "outcar": {str(k).upper(): v for k, v in (getattr(outcar, "executed_tags", None) or {}).items()},
            "incar_file": _tags_of(incar_file),
        }
        self.present = {name for name, table in self.tables.items() if table}
        if incar_file is not None:
            self.present.add("incar_file")
        if outcar is not None:
            self.present.add("outcar")
        self.used: set[str] = set()
        self.conflicts: list[dict[str, Any]] = []
        self.problems: list[str] = []

    def lookup(self, field_name: str, tags: Sequence[tuple[str, str]], kind: str,
               order: Sequence[str] = SOURCE_ORDER, normalize=None) -> tuple[Any, str | None]:
        """First coercible value; later sources that disagree are recorded as conflicts.

        ``tags`` maps source -> tag name (``("parameters", "ENMAX")``...); a
        source may appear once. Values that cannot be coerced are recorded as
        problems and skipped. ``normalize`` (e.g. PREC 'Acc' -> 'accurate') is
        applied before values are returned and compared.
        """
        coerce = _COERCE[kind]
        tag_for = dict(tags)
        found: list[tuple[str, Any, Any]] = []
        for source in order:
            tag = tag_for.get(source)
            if tag is None or tag not in self.tables[source]:
                continue
            raw = self.tables[source][tag]
            value = coerce(raw)
            if value is not None and normalize is not None:
                value = normalize(value)
            if value is None:
                if raw is not None:
                    self.problems.append(f"{field_name}: {source} {tag}={raw!r} could not be read as {kind}")
                continue
            found.append((source, value, raw))
        if not found:
            return None, None
        source, value, raw = found[0]
        self.used.add(source)
        disagree = {
            other: json_safe(other_value)
            for other, other_value, other_raw in found[1:]
            if not values_agree(value, other_value, raw_a=raw, raw_b=other_raw)
        }
        if disagree:
            self.conflicts.append({"field": field_name, "value": json_safe(value), "source": source,
                                   "others": disagree})
        return value, source


def _tag_everywhere(tag: str, *, parameters_tag: str | None = None) -> list[tuple[str, str]]:
    return [("parameters", parameters_tag or tag), ("vasprun_incar", tag), ("outcar", tag), ("incar_file", tag)]


# --------------------------------------------------------------------------
# Extraction
# --------------------------------------------------------------------------

def _potcar_lexch(potcar) -> str | None:
    datasets = list(getattr(potcar, "datasets", None) or [])
    values = {as_str(d.get("lexch")) for d in datasets if isinstance(d, Mapping)}
    values.discard(None)
    if len(values) == 1:
        return str(values.pop()).upper()
    return None


def _titel_family_lexch(titles: Sequence[str | None]) -> str | None:
    families = set()
    for titel in titles:
        if not titel:
            return None
        families.add(titel.split()[0].upper())
    if len(families) != 1:
        return None
    return _TITEL_FAMILY_LEXCH.get(families.pop())


def _prec(value: str | None) -> str | None:
    if not value:
        return None
    text = value.strip().lower()
    if text.startswith("singlen"):
        return "singlen"
    return _PREC.get(text[:1], text)


def _metagga(value: str | None) -> str | None:
    if value is None:
        return None
    text = value.strip().upper()
    return "--" if text in {"--", "NONE", "NO", ""} else text


def _dataset_sha(dataset: Mapping[str, Any]) -> str | None:
    return as_str(dataset.get("sha256_dataset_bytes") or dataset.get("dataset_sha256"))


def _potcar_entries(header, potcar, outcar, provenance: dict, hints: dict) -> dict[str, Any]:
    """Per-element {titel, sha256}; a list of entries if an element has several species blocks."""
    atom_types = list(getattr(header, "atom_types", None) or [])
    outcar_titles = [normalize_titel(t) for t in (getattr(outcar, "potcar_titles", None) or [])]
    datasets = [d for d in (getattr(potcar, "datasets", None) or []) if isinstance(d, Mapping)]
    usable_datasets = potcar is not None and len(datasets) == len(atom_types)
    if potcar is not None and not usable_datasets:
        hints["potcar_file"] = (
            f"POTCAR has {len(datasets)} datasets for {len(atom_types)} species; dataset hashes not attached"
        )
    per_element: dict[str, list[dict[str, Any]]] = {}
    for index, atom_type in enumerate(atom_types):
        element = as_str(atom_type.element) or UNKNOWN
        titel = normalize_titel(atom_type.pseudopotential)
        titel_source = "atominfo"
        if not titel and index < len(outcar_titles) and len(outcar_titles) == len(atom_types):
            titel, titel_source = outcar_titles[index], "outcar"
        if not titel and usable_datasets:
            titel, titel_source = normalize_titel(datasets[index].get("titel")), "potcar_file"
        sha = None
        if usable_datasets:
            dataset = datasets[index]
            same_element = as_str(dataset.get("element")) in {None, element}
            same_titel = titel is None or normalize_titel(dataset.get("titel")) in {None, titel}
            if same_element and same_titel:
                sha = _dataset_sha(dataset)
            else:
                hints["potcar_file"] = (
                    f"POTCAR dataset {index} ({dataset.get('titel')!r}) does not match species "
                    f"{element} ({titel!r}); dataset hashes not attached"
                )
        entry = {"titel": titel or UNKNOWN, "sha256": sha}
        provenance[f"potcar[{element}]"] = titel_source + ("+potcar_file" if sha else "")
        if entry not in per_element.setdefault(element, []):
            per_element[element].append(entry)
    if hints.get("potcar_file"):
        for entries in per_element.values():
            for entry in entries:
                entry["sha256"] = None
    return {
        element: entries[0] if len(entries) == 1 else sorted(entries, key=canonical_json)
        for element, entries in sorted(per_element.items())
    }


def _hubbard(header, sources: _Sources, provenance: dict) -> tuple[Any, Any]:
    """(per-element {el: [L, U, J]} or UNKNOWN, LDAUTYPE or None/UNKNOWN)."""
    atom_types = list(getattr(header, "atom_types", None) or [])
    ldau, ldau_source = sources.lookup("LDAU", _tag_everywhere("LDAU"), "bool")
    provenance["LDAU"] = ldau_source or "default(not set)"
    if ldau is None:
        # LDAU is always written to <parameters>; absence everywhere means VASP's default (off)
        # only when <parameters> itself was read.
        if not sources.tables["parameters"]:
            return UNKNOWN, UNKNOWN
        ldau = False
    if not ldau:
        return {}, None
    ldautype, type_source = sources.lookup("LDAUTYPE", _tag_everywhere("LDAUTYPE"), "int")
    arrays = {}
    for tag, kind in (("LDAUL", "int_list"), ("LDAUU", "float_list"), ("LDAUJ", "float_list")):
        arrays[tag], arrays[tag + "_source"] = sources.lookup(tag, _tag_everywhere(tag), kind)
    provenance["LDAUTYPE"] = type_source
    provenance["hubbard"] = ",".join(sorted({str(arrays[t + "_source"]) for t in ("LDAUL", "LDAUU", "LDAUJ")}))
    if any(arrays[tag] is None or len(arrays[tag]) != len(atom_types) for tag in ("LDAUL", "LDAUU", "LDAUJ")):
        return UNKNOWN, ldautype if ldautype is not None else UNKNOWN
    per_element: dict[str, list[list[float]]] = {}
    for index, atom_type in enumerate(atom_types):
        l_value = int(arrays["LDAUL"][index])
        u_value = canonical_float(arrays["LDAUU"][index])
        j_value = canonical_float(arrays["LDAUJ"][index])
        if l_value < 0 or (u_value == 0.0 and j_value == 0.0):
            continue
        triple = [l_value, u_value, j_value]
        entries = per_element.setdefault(as_str(atom_type.element) or UNKNOWN, [])
        if triple not in entries:
            entries.append(triple)
    hubbard = {el: (v[0] if len(v) == 1 else sorted(v)) for el, v in sorted(per_element.items())}
    if not hubbard:
        return {}, None
    return hubbard, ldautype if ldautype is not None else UNKNOWN


def extract_settings(header, *, incar_file=None, outcar=None, potcar=None) -> dict[str, Any]:
    """Method, sampling, unknown and provenance records of one run (JSON-safe).

    ``header`` is a :class:`~nio_md_prep.dataset.model.VasprunHeader`;
    ``incar_file`` an :class:`~nio_md_prep.dataset.model.IncarEvidence` (or a
    tag -> value mapping); ``outcar`` an OutcarEvidence; ``potcar`` a
    PotcarEvidence. Missing evidence is simply not used.
    """
    sources = _Sources(header, incar_file, outcar)
    provenance: dict[str, Any] = {}
    hints: dict[str, str] = {}
    method: dict[str, Any] = {}
    atom_types = list(getattr(header, "atom_types", None) or [])
    titles = [normalize_titel(t.pseudopotential) for t in atom_types]

    def put(name: str, value: Any, source: str | None) -> None:
        method[name] = UNKNOWN if value is None else value
        provenance[name] = source

    # functional ----------------------------------------------------------
    gga, gga_source = sources.lookup("GGA", _tag_everywhere("GGA"), "str", normalize=str.upper)
    if gga in {None, "--"}:
        resolved = _potcar_lexch(potcar)
        if resolved:
            gga, gga_source = resolved, "potcar_lexch"
        else:
            resolved = _titel_family_lexch(titles) if titles else None
            if resolved:
                gga, gga_source = resolved, "potcar_titel_family"
            else:
                if gga == "--":
                    hints["GGA"] = "GGA=-- (POTCAR default) but the POTCAR LEXCH could not be determined"
                gga, gga_source = None, None
    put("GGA", gga, gga_source)

    metagga, metagga_source = sources.lookup("METAGGA", _tag_everywhere("METAGGA"), "str", normalize=_metagga)
    if metagga is None:
        # VASP 6.3 <parameters> has LMETAGGA but no METAGGA; <incar> lists METAGGA when it was set.
        lmetagga = as_bool(sources.tables["parameters"].get("LMETAGGA"))
        if lmetagga is False:
            metagga, metagga_source = "--", "parameters:LMETAGGA"
        elif lmetagga is None and sources.tables["parameters"]:
            metagga, metagga_source = "--", "default(not in <incar>)"
    put("METAGGA", metagga, metagga_source)

    lhfcalc, source = sources.lookup("LHFCALC", _tag_everywhere("LHFCALC"), "bool")
    put("LHFCALC", lhfcalc, source)
    for name in ("AEXX", "HFSCREEN"):
        value, source = sources.lookup(name, _tag_everywhere(name), "float")
        if lhfcalc is False:
            method[name], provenance[name] = None, "not_applicable(LHFCALC=F)"
        else:
            put(name, canonical_float(value) if value is not None else None, source)

    # dispersion ----------------------------------------------------------
    ivdw, ivdw_source = sources.lookup(
        "IVDW", [("parameters", "IVDW"), ("vasprun_incar", "IVDW"), ("outcar", "IVDW"), ("incar_file", "IVDW")],
        "int",
    )
    if ivdw is None:
        # Legacy DFT-D2 switch: LVDW=.TRUE. is IVDW=1.
        lvdw, lvdw_source = sources.lookup("LVDW", _tag_everywhere("LVDW"), "bool")
        if lvdw is True:
            ivdw, ivdw_source = 1, f"{lvdw_source}:LVDW"
    if ivdw is None and "incar_file" in sources.present and "IVDW" not in sources.tables["incar_file"]:
        ivdw, ivdw_source = 0, "incar_file(not set)"
    if ivdw is None:
        version = getattr(header, "version_tuple", None)
        if version and version[0] >= 6 and "IVDW" not in sources.tables["vasprun_incar"]:
            hints["IVDW"] = "absent from <incar> (VASP >= 6 normally writes set tags there); still unknown"
        else:
            hints["IVDW"] = "IVDW is not in <parameters>; no <incar>, OUTCAR or INCAR value found"
    put("IVDW", ivdw, ivdw_source)
    vdw_params: dict[str, Any] = {}
    if isinstance(ivdw, int) and ivdw != 0:
        for source_name in ("incar_file", "outcar", "vasprun_incar"):  # later overrides earlier
            for tag, raw in sources.tables[source_name].items():
                if _VDW_TAG.match(tag):
                    number = as_float(raw)
                    flag = as_bool(raw)
                    vdw_params[tag] = (
                        flag if flag is not None else canonical_float(number) if number is not None else as_str(raw)
                    )
    method["IVDW_params"] = dict(sorted(vdw_params.items()))
    luse_vdw, source = sources.lookup("LUSE_VDW", _tag_everywhere("LUSE_VDW"), "bool")
    if luse_vdw is None and sources.tables["parameters"]:
        luse_vdw, source = False, "default(not set)"
    put("LUSE_VDW", luse_vdw, source)

    # basis, spin, smearing ------------------------------------------------
    encut, source = sources.lookup(
        "ENCUT", [("parameters", "ENMAX"), ("vasprun_incar", "ENCUT"), ("outcar", "ENCUT"), ("incar_file", "ENCUT")],
        "float",
    )
    put("ENCUT", canonical_float(encut) if encut is not None else None, source)
    prec, source = sources.lookup("PREC", _tag_everywhere("PREC"), "str", normalize=_prec)
    put("PREC", prec, source)
    for name, kind in (("LASPH", "bool"), ("ISPIN", "int"), ("LSORBIT", "bool"), ("LNONCOLLINEAR", "bool"),
                       ("ISMEAR", "int"), ("LDIPOL", "bool"), ("IDIPOL", "int")):
        value, source = sources.lookup(name, _tag_everywhere(name), kind)
        put(name, value, source)
    sigma, source = sources.lookup("SIGMA", _tag_everywhere("SIGMA"), "float")
    put("SIGMA", canonical_float(sigma) if sigma is not None else None, source)
    nupdown, source = sources.lookup("NUPDOWN", _tag_everywhere("NUPDOWN"), "float")
    ispin = method.get("ISPIN")
    if ispin == 1:
        method["NUPDOWN"], provenance["NUPDOWN"] = None, "not_applicable(ISPIN=1)"
    elif nupdown is not None and nupdown < 0:
        method["NUPDOWN"], provenance["NUPDOWN"] = None, f"{source}(not fixed)"
    else:
        put("NUPDOWN", canonical_float(nupdown) if nupdown is not None else None, source)

    # net charge ------------------------------------------------------------
    nelect, nelect_source = sources.lookup("NELECT", _tag_everywhere("NELECT"), "float")
    datasets = [d for d in (getattr(potcar, "datasets", None) or []) if isinstance(d, Mapping)]
    zvals: list[float | None] = []
    for index, atom_type in enumerate(atom_types):
        zval = as_float(atom_type.valence)
        if zval is None and len(datasets) == len(atom_types):
            zval = as_float(datasets[index].get("zval"))
        zvals.append(zval)
    if nelect is not None and atom_types and all(z is not None for z in zvals):
        total = sum(float(z) * int(t.count) for z, t in zip(zvals, atom_types))  # type: ignore[arg-type]
        put("net_charge", round(total - nelect, 6) + 0.0, f"atominfo+{nelect_source}")
    else:
        put("net_charge", None, None)

    # per element -----------------------------------------------------------
    method["potcar"] = _potcar_entries(header, potcar, outcar, provenance, hints)
    hubbard, ldautype = _hubbard(header, sources, provenance)
    method["hubbard"] = hubbard
    method["LDAUTYPE"] = ldautype

    unknown = sorted(
        [name for name in GLOBAL_METHOD_FIELDS if method.get(name) == UNKNOWN]
        + (["hubbard"] if hubbard == UNKNOWN else [])
        + (["LDAUTYPE"] if ldautype == UNKNOWN else [])
        + [f"potcar[{el}]" for el, entry in method["potcar"].items()
           if isinstance(entry, dict) and entry.get("titel") == UNKNOWN]
    )

    # sampling ----------------------------------------------------------------
    sampling: dict[str, Any] = {"KPOINTS": json_safe(getattr(header, "kpoints", None) or {})}
    sampling_sources: dict[str, str | None] = {}
    for tag, kind in _SAMPLING_TAGS:
        value, source = sources.lookup(tag, _tag_everywhere(tag), kind)
        if value is not None:
            sampling[tag] = json_safe(value)
            sampling_sources[tag] = source
    # vasprun prints EDIFF with 8 decimals: EDIFF < 5e-9 appears as 0.
    if sampling.get("EDIFF") == 0.0 and sampling_sources.get("EDIFF") in {"parameters", "vasprun_incar"}:
        for source_name in ("outcar", "incar_file"):
            precise = as_float(sources.tables[source_name].get("EDIFF"))
            if precise:
                sampling["EDIFF"], sampling_sources["EDIFF"] = precise, source_name
                hints["EDIFF"] = "vasprun prints EDIFF with 8 decimals (0.00000000); value taken from " + source_name
                break
        else:
            hints["EDIFF"] = "EDIFF printed as 0 in vasprun: either disabled or below 5e-9"
    sampling["parameters"] = json_safe(getattr(header, "parameters", None) or {})

    incar_files_conflicts = [c for c in sources.conflicts if "incar_file" in c["others"] or c["source"] == "incar_file"]
    return {
        "schema": SETTINGS_SCHEMA,
        "method": json_safe(method),
        "sampling": sampling,
        "unknown": unknown,
        "provenance": {
            "fields": json_safe(provenance),
            "sampling_fields": sampling_sources,
            "sources": sorted(sources.used | {s for s in ("potcar_file",) if potcar is not None}),
            "conflicts": sorted(sources.conflicts, key=canonical_json),
            "incar_file_differs": bool(incar_files_conflicts),
            "problems": sorted(set(sources.problems)),
            "hints": dict(sorted(hints.items())),
            "potcar": _potcar_detail(potcar),
        },
    }


def _potcar_detail(potcar) -> dict[str, Any] | None:
    """Recorded (not compared) POTCAR metadata: never any content."""
    if potcar is None:
        return None
    keep = ("symbol", "element", "titel", "vrhfin", "lexch", "zval", "pomass", "enmax", "enmin",
            "sha256_header", "sha256_header_verified", "sha256_dataset_bytes")
    return {
        "sha256_file": getattr(potcar, "sha256", None),
        "datasets": [
            {key: json_safe(d.get(key)) for key in keep if key in d}
            for d in (getattr(potcar, "datasets", None) or []) if isinstance(d, Mapping)
        ],
    }


#: Numeric tags whose OUTCAR echo must agree with vasprun (anything else is only recorded).
OUTCAR_CHECKED_FIELDS: tuple[str, ...] = (
    "ENCUT", "ISPIN", "ISMEAR", "SIGMA", "NUPDOWN", "NELECT", "IVDW", "LDAUL", "LDAUU", "LDAUJ", "LDAUTYPE",
    "NELM", "EDIFF", "IBRION", "NSW", "ISIF", "POTIM", "PSTRESS",
)


def outcar_conflicts(settings: Mapping[str, Any], fields: Sequence[str] = OUTCAR_CHECKED_FIELDS) -> list[dict[str, Any]]:
    """Conflicts between vasprun values and the OUTCAR echo for numeric ``fields`` (a different-run signal)."""
    wanted = set(fields)
    result = []
    for conflict in settings.get("provenance", {}).get("conflicts", []):
        if (
            conflict.get("field") in wanted
            and conflict.get("source") in {"parameters", "vasprun_incar"}
            and "outcar" in conflict.get("others", {})
        ):
            result.append(conflict)
    return result


# --------------------------------------------------------------------------
# Overrides
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Equivalence:
    field: str
    values: tuple[str, ...]
    reason: str
    reviewed_by: str

    @property
    def tokens(self) -> frozenset[str]:
        return frozenset(_normalize_token(self.field, value) for value in self.values)

    @property
    def representative(self) -> str:
        return "equiv(" + "|".join(sorted(self.tokens)) + ")"

    def as_dict(self) -> dict[str, Any]:
        return {"field": self.field, "values": list(self.values), "reason": self.reason,
                "reviewed_by": self.reviewed_by}


@dataclass(frozen=True)
class Overrides:
    entries: tuple[Equivalence, ...] = ()
    path: str | None = None
    sha256: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return {"path": self.path, "sha256": self.sha256, "equivalence": [e.as_dict() for e in self.entries]}

    def classes(self, field_name: str) -> list[Equivalence]:
        return [entry for entry in self.entries if entry.field == field_name]


NO_OVERRIDES = Overrides()


def _normalize_token(field_name: str, value: Any) -> str:
    """Comparable token for an override value or a method value."""
    if value is None:
        return "none"
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, (int, float)):
        return format(canonical_float(value), ".10g")
    if isinstance(value, (list, tuple)):
        return ",".join(_normalize_token(field_name, item) for item in value)
    if isinstance(value, Mapping):
        return canonical_json(json_safe(value))
    text = " ".join(str(value).split())
    if text == UNKNOWN:
        return UNKNOWN
    if field_name == "hubbard":
        return ",".join(_normalize_token("", part.strip()) for part in text.split(","))
    logical = as_bool(text) if text.upper().strip(".") in {"T", "F", "TRUE", "FALSE"} else None
    if logical is not None:
        return "T" if logical else "F"
    if text.lower() in {"none", "null"}:
        return "none"
    number = as_float(text) if _NUMBER.match(text) else None
    if number is not None:
        return format(canonical_float(number), ".10g")
    return text.upper()


def load_overrides(path: Path | None) -> Overrides:
    """Read an overrides TOML file (``[[equivalence]]`` tables); None -> no overrides."""
    if path is None:
        return NO_OVERRIDES
    import tomllib

    path = Path(path)
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise DatasetError(f"cannot read overrides {path.as_posix()}: {exc}") from exc
    unexpected = sorted(set(data) - {"equivalence"})
    if unexpected:
        raise DatasetError(f"overrides {path.as_posix()}: unknown top-level keys {unexpected}")
    raw_entries = data.get("equivalence", [])
    if not isinstance(raw_entries, list):
        raise DatasetError(f"overrides {path.as_posix()}: 'equivalence' must be an array of tables")
    entries: list[Equivalence] = []
    required = {"field", "values", "reason", "reviewed_by"}
    for number, raw in enumerate(raw_entries, 1):
        where = f"overrides {path.as_posix()} [[equivalence]] #{number}"
        if not isinstance(raw, dict):
            raise DatasetError(f"{where}: not a table")
        missing = sorted(required - set(raw))
        extra = sorted(set(raw) - required)
        if missing or extra:
            raise DatasetError(f"{where}: missing keys {missing}, unknown keys {extra}")
        field_name = raw["field"]
        if field_name not in OVERRIDABLE_FIELDS:
            raise DatasetError(f"{where}: field {field_name!r} is not one of {list(OVERRIDABLE_FIELDS)}")
        values = raw["values"]
        if not isinstance(values, list) or len(values) < 2 or not all(
            isinstance(v, (str, int, float, bool)) for v in values
        ):
            raise DatasetError(f"{where}: 'values' must list at least two scalar values")
        for key in ("reason", "reviewed_by"):
            if not isinstance(raw[key], str) or not raw[key].strip():
                raise DatasetError(f"{where}: {key!r} must be a non-empty string")
        entry = Equivalence(
            field=field_name,
            values=tuple(("T" if v else "F") if isinstance(v, bool) else str(v) for v in values),
            reason=raw["reason"].strip(),
            reviewed_by=raw["reviewed_by"].strip(),
        )
        if len(entry.tokens) < 2:
            raise DatasetError(f"{where}: values {list(entry.values)} name only one distinct value")
        for previous in entries:
            if previous.field == field_name and previous.tokens & entry.tokens:
                raise DatasetError(
                    f"{where}: values overlap an earlier equivalence for {field_name}: "
                    f"{sorted(previous.tokens & entry.tokens)}"
                )
        entries.append(entry)
    return Overrides(entries=tuple(entries), path=path.as_posix(), sha256=sha256_file(path))


def _apply(field_name: str, value: Any, overrides: Overrides) -> Any:
    token = _normalize_token(field_name, value)
    for entry in overrides.classes(field_name):
        if token in entry.tokens:
            return entry.representative
    return value


def canonical_method(settings_or_method: Mapping[str, Any], overrides: Overrides | None = None) -> dict[str, Any]:
    """The method map after override canonicalisation (a new dict; input untouched)."""
    method = settings_or_method.get("method", settings_or_method)
    method = copy.deepcopy(dict(method))
    overrides = overrides or NO_OVERRIDES
    if not overrides.entries:
        return method
    for name in GLOBAL_METHOD_FIELDS + CONDITIONAL_METHOD_FIELDS:
        if name in method:
            method[name] = _apply(name, method[name], overrides)
    potcar = method.get("potcar") or {}
    for element, entry in potcar.items():
        entries = entry if isinstance(entry, list) else [entry]
        for item in entries:
            if isinstance(item, dict):
                item["titel"] = _apply("potcar.titel", item.get("titel"), overrides)
                sha = item.get("sha256")
                replaced = _apply("potcar.sha256", UNKNOWN if sha is None else sha, overrides)
                item["sha256"] = None if replaced == UNKNOWN and sha is None else replaced
    hubbard = method.get("hubbard")
    if isinstance(hubbard, dict) and overrides.classes("hubbard"):
        for element in sorted(potcar):
            current = hubbard.get(element)
            replaced = _apply("hubbard", "none" if current is None else current, overrides)
            if replaced != ("none" if current is None else current):
                hubbard[element] = replaced
    return method


def method_fingerprint(settings: Mapping[str, Any], overrides: Overrides | None = None) -> str:
    """First 12 hex digits of sha256 over the canonical JSON of the (overridden) method map."""
    return hashlib.sha256(canonical_json(canonical_method(settings, overrides)).encode()).hexdigest()[:12]


# --------------------------------------------------------------------------
# Pools
# --------------------------------------------------------------------------

@dataclass
class Pool:
    """Runs whose (overridden) method maps are mutually compatible."""

    pool_id: str
    method: dict[str, Any]  # merged map: global fields, 'elements', 'potcar', 'hubbard', 'LDAUTYPE'
    run_ids: list[str] = field(default_factory=list)
    fingerprints: list[str] = field(default_factory=list)

    @property
    def unknown_fields(self) -> list[str]:
        return pool_unknown_fields(self.method)

    def as_dict(self) -> dict[str, Any]:
        return {"pool_id": self.pool_id, "method": self.method, "run_ids": list(self.run_ids),
                "fingerprints": list(self.fingerprints), "unknown_fields": self.unknown_fields}


def _elements(method: Mapping[str, Any]) -> list[str]:
    return sorted(method.get("elements") or (method.get("potcar") or {}).keys())


def _has_u(method: Mapping[str, Any]) -> bool:
    hubbard = method.get("hubbard")
    return isinstance(hubbard, dict) and bool(hubbard)


def _differences(a: Mapping[str, Any], b: Mapping[str, Any]) -> list[dict[str, Any]]:
    diffs: list[dict[str, Any]] = []
    for name in GLOBAL_METHOD_FIELDS:
        if a.get(name) != b.get(name):
            diffs.append({"field": name, "a": a.get(name), "b": b.get(name), "kind": "conflict"})
    if _has_u(a) and _has_u(b) and a.get("LDAUTYPE") != b.get("LDAUTYPE"):
        diffs.append({"field": "LDAUTYPE", "a": a.get("LDAUTYPE"), "b": b.get("LDAUTYPE"), "kind": "conflict"})
    elements_a, elements_b = set(_elements(a)), set(_elements(b))
    shared = sorted(elements_a & elements_b)
    potcar_a, potcar_b = a.get("potcar") or {}, b.get("potcar") or {}
    for element in shared:
        if potcar_a.get(element) != potcar_b.get(element):
            diffs.append({"field": f"potcar[{element}]", "a": potcar_a.get(element), "b": potcar_b.get(element),
                          "kind": "conflict"})
    hubbard_a, hubbard_b = a.get("hubbard"), b.get("hubbard")
    if (hubbard_a == UNKNOWN) != (hubbard_b == UNKNOWN):
        if shared:
            diffs.append({"field": "hubbard", "a": hubbard_a, "b": hubbard_b, "kind": "conflict"})
    elif isinstance(hubbard_a, dict) and isinstance(hubbard_b, dict):
        for element in shared:
            if hubbard_a.get(element) != hubbard_b.get(element):
                diffs.append({"field": f"hubbard[{element}]", "a": hubbard_a.get(element),
                              "b": hubbard_b.get(element), "kind": "conflict"})
    if elements_a != elements_b:
        diffs.append({"field": "elements", "a": sorted(elements_a - elements_b), "b": sorted(elements_b - elements_a),
                      "kind": "coverage"})
    return sorted(diffs, key=lambda d: (d["kind"] != "conflict", d["field"]))


def _merge(pool: dict[str, Any], method: Mapping[str, Any]) -> None:
    pool["elements"] = sorted(set(_elements(pool)) | set(_elements(method)))
    pool.setdefault("potcar", {}).update({k: v for k, v in (method.get("potcar") or {}).items()})
    pool["potcar"] = dict(sorted(pool["potcar"].items()))
    if isinstance(pool.get("hubbard"), dict) and isinstance(method.get("hubbard"), dict):
        pool["hubbard"] = dict(sorted({**pool["hubbard"], **method["hubbard"]}.items()))
    if pool.get("LDAUTYPE") is None and method.get("LDAUTYPE") is not None:
        pool["LDAUTYPE"] = method.get("LDAUTYPE")


def _new_pool_map(method: Mapping[str, Any]) -> dict[str, Any]:
    pool = copy.deepcopy(dict(method))
    pool["elements"] = _elements(method)
    return pool


def _run_items(runs: Iterable[Any]) -> list[tuple[str, Mapping[str, Any]]]:
    items = []
    for run in runs:
        if isinstance(run, tuple):
            run_id, settings = run
        else:
            run_id, settings = run.run_id, run.settings
        items.append((str(run_id), settings))
    return sorted(items, key=lambda item: item[0])


def build_pools(runs: Iterable[Any], *, overrides: Overrides | None = None) -> list[Pool]:
    """Element-aware greedy pooling of runs (sorted by run_id before processing).

    ``runs`` holds ``(run_id, settings)`` pairs or objects with ``run_id`` and
    ``settings`` (e.g. :class:`~nio_md_prep.dataset.model.RunRecord`). A run
    joins the first pool whose global fields equal its own and whose merged
    per-element POTCAR/+U maps agree on every shared element; otherwise it
    starts a new pool. ``pool_id`` = sha256[:12] of the merged map.
    """
    pools: list[dict[str, Any]] = []
    for run_id, settings in _run_items(runs):
        method = canonical_method(settings, overrides)
        fingerprint = method_fingerprint(settings, overrides)
        for pool in pools:
            if not [d for d in _differences(pool["method"], method) if d["kind"] == "conflict"]:
                _merge(pool["method"], method)
                pool["run_ids"].append(run_id)
                pool["fingerprints"].add(fingerprint)
                break
        else:
            pools.append({"method": _new_pool_map(method), "run_ids": [run_id], "fingerprints": {fingerprint}})
    result = []
    for pool in pools:
        pool_method = json_safe(pool["method"])
        pool_id = hashlib.sha256(canonical_json(pool_method).encode()).hexdigest()[:12]
        result.append(Pool(pool_id=pool_id, method=pool_method, run_ids=list(pool["run_ids"]),
                           fingerprints=sorted(pool["fingerprints"])))
    return result


def assign_pools(runs: Iterable[Any], *, overrides: Overrides | None = None) -> dict[str, list[str]]:
    """``{pool_id: [run_id, ...]}`` (see :func:`build_pools`)."""
    return {pool.pool_id: list(pool.run_ids) for pool in build_pools(runs, overrides=overrides)}


def describe_pool_differences(pool_a: Pool | Mapping[str, Any], pool_b: Pool | Mapping[str, Any]) -> list[dict[str, Any]]:
    """Field-level differences between two pools (or method maps / settings dicts).

    Each entry is ``{"field", "a", "b", "kind"}`` with kind ``"conflict"`` for
    values that prevent pooling and ``"coverage"`` for elements present in only
    one side (never a conflict by itself). Sorted: conflicts first, by field.
    """
    def as_map(value):
        if isinstance(value, Pool):
            return value.method
        return value.get("method", value) if isinstance(value, Mapping) else value

    return _differences(as_map(pool_a), as_map(pool_b))


def pool_unknown_fields(method: Mapping[str, Any]) -> list[str]:
    """Method fields of a pool (or run) whose value is :data:`UNKNOWN` (or a missing POTCAR hash)."""
    unknown = [name for name in GLOBAL_METHOD_FIELDS + CONDITIONAL_METHOD_FIELDS if method.get(name) == UNKNOWN]
    if method.get("hubbard") == UNKNOWN:
        unknown.append("hubbard")
    for element, entry in sorted((method.get("potcar") or {}).items()):
        entries = entry if isinstance(entry, list) else [entry]
        for item in entries:
            if isinstance(item, Mapping) and item.get("titel") == UNKNOWN:
                unknown.append(f"potcar[{element}].titel")
            if isinstance(item, Mapping) and item.get("sha256") is None:
                unknown.append(f"potcar[{element}].sha256")
    return sorted(set(unknown))
