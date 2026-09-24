"""``dataset scan`` and ``dataset export``: discovery -> classification -> accounting -> extxyz.

Pipeline (nothing is dropped silently; every candidate ends up in the accounting files):

1. **discovery** (:func:`discovery.discover_runs`) under aliased roots with the default exclusion
   rules; every run candidate inside a pruned tree becomes an excluded run record;
2. **inventory** selection (optional TOML; unlisted/excluded runs are recorded, not parsed);
3. **one streaming pass per run**: the evidence files (OUTCAR, OSZICAR, POTCAR, INCAR, POSCAR) are
   parsed once, then vasprun.xml is streamed once; settings (:func:`settings.extract_settings`),
   the run decision (:func:`acceptance.classify_run`, with agglomeration reference hashes), one
   decision per ionic step (:func:`acceptance.classify_steps`) and structure keys are computed
   from that single pass. Arrays of frames accepted by the label rules are spilled to a private
   work directory (``.work-<pid>`` inside the output directory), so memory holds one run;
4. **reference-settings pools** over runs with accepted frames (:func:`settings.build_pools`);
   ``export`` needs exactly one pool or ``--pool``; other pools -> ``reference_pool_not_selected``;
5. **exact duplicates** within a pool (:func:`duplicates.find_exact_duplicates`) plus
   shared-structure links over every keyed frame (:func:`duplicates.structure_links`);
6. **lineage** (:func:`lineage.resolve_lineage`) with those links (and near-duplicate links when
   enabled); frames of runs whose lineage is unresolved are ``quarantined/lineage_unresolved``
   (they could not be split safely);
7. declared **subsampling** (stride, max frames per run) -> ``excluded/subsampled``;
8. **accounting** for every candidate: ``runs.jsonl``, ``frames.jsonl``, ``exclusions.csv``,
   ``audit_report.json`` + ``audit.md`` (summaries by every audit dimension);
9. ``export`` only: ``dataset.extxyz`` written from the spill and verified by an ASE round trip
   before it is published (:func:`extxyz.write_extxyz`), then ``dataset_manifest.json`` last with
   file hashes, statistics, options, policy, tool/parser version, git commit and a
   ``content_sha256`` that excludes the ``invocation`` block (argv, created_at).

All outputs are staged in the work directory and moved into place only when complete (manifest
last); the work directory is always removed. Output bytes depend only on the data, the options
and the code version: never on root order, directory enumeration order or wall-clock time.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import datetime
import io
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np

from . import PARSER_NAME, PARSER_VERSION, SCHEMA_VERSION
from . import acceptance as acc
from . import discovery as disc
from . import duplicates as dup
from . import extxyz as xyz
from . import lineage as lin
from . import settings as st
from . import vaspfiles
from . import vasprun as vr
from .errors import DatasetError, LeakageError
from .fsio import atomic_write_json, atomic_write_text, canonical_json, git_provenance
from .fsio import prepare_output_dir, sha256_bytes, sha256_file, write_jsonl
from .model import ACCEPTED, Outcome, outcome_for

EXPORT_SCHEMA = "nio-md-prep.dataset-export"
SCAN_SCHEMA = "nio-md-prep.dataset-scan"
MANIFEST_SCHEMA_VERSION = 1
DATASET_FILE = "dataset.extxyz"
RUNS_FILE = "runs.jsonl"
FRAMES_FILE = "frames.jsonl"
EXCLUSIONS_FILE = "exclusions.csv"
AUDIT_JSON = "audit_report.json"
AUDIT_MD = "audit.md"
EXPORT_MANIFEST = "dataset_manifest.json"
SCAN_MANIFEST = "scan_manifest.json"
#: Manifest keys outside ``content_sha256`` (they change with every invocation).
NON_CONTENT_KEYS = ("content_sha256", "invocation")
EXCLUSION_COLUMNS = ("level", "id", "run_id", "status", "reason", "detail")
UNKNOWN_TEXT = "unknown"
UNITS = {
    "energy": "eV", "forces": "eV/Angstrom", "stress": "eV/Angstrom^3 (ASE sign: positive = tensile)",
    "positions": "Angstrom", "cell": "Angstrom", "time": "fs", "magnetization": "muB",
}
_MIN_POOL_PREFIX = 6


def _tool_version() -> str:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover
        return UNKNOWN_TEXT
    try:
        return version("nio-md-prep")
    except PackageNotFoundError:
        return UNKNOWN_TEXT


def _composition(species: Sequence[str]) -> tuple[str, str]:
    counts: dict[str, int] = {}
    for symbol in species:
        counts[symbol] = counts.get(symbol, 0) + 1
    return "".join(f"{s}{counts[s]}" for s in sorted(counts)), "-".join(sorted(counts))


def _species_blocks(species: Sequence[str]) -> list[list[Any]]:
    blocks: list[list[Any]] = []
    for symbol in species:
        if blocks and blocks[-1][0] == symbol:
            blocks[-1][1] += 1
        else:
            blocks.append([symbol, 1])
    return blocks


def _count(values: Iterable[Any]) -> dict[str, int]:
    result: dict[str, int] = {}
    for value in values:
        key = str(value)
        result[key] = result.get(key, 0) + 1
    return dict(sorted(result.items()))


# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------

@dataclass
class ExportOptions:
    """Everything that decides the content of a scan/export (recorded in the manifest)."""

    inventory: Path | None = None
    settings_overrides: Path | None = None
    lineage_policy: str | None = None
    pool: str | None = None
    policy: acc.Policy = field(default_factory=acc.Policy)
    label_keys: xyz.LabelKeys = field(default_factory=xyz.LabelKeys)
    include_parts: Sequence[str] = ()
    exclude_globs: Sequence[str] = ()
    include_globs: Sequence[str] = ()
    follow_symlinks: bool = False
    stride: int = 1
    max_frames_per_run: int | None = None

    def __post_init__(self):
        if isinstance(self.stride, bool) or not isinstance(self.stride, int) or self.stride < 1:
            raise DatasetError(f"--stride must be a positive integer, got {self.stride!r}")
        if self.max_frames_per_run is not None and (
                isinstance(self.max_frames_per_run, bool) or not isinstance(self.max_frames_per_run, int)
                or self.max_frames_per_run < 1):
            raise DatasetError(f"--max-frames-per-run must be a positive integer, got {self.max_frames_per_run!r}")
        lin.parse_lineage_policy(self.lineage_policy)
        xyz.validate_label_keys(self.label_keys.energy, self.label_keys.forces, self.label_keys.stress)
        self.include_parts = disc.check_include_parts(self.include_parts)
        self.exclude_globs = tuple(self.exclude_globs)
        self.include_globs = tuple(self.include_globs)
        if self.policy.magnetic_overrides and not self.policy.magnetic_override_reason:
            raise DatasetError(
                f"magnetic overrides {self.policy.magnetic_overrides} need a recorded reason "
                "(--magnetic-override-reason TEXT)"
            )

    def as_dict(self) -> dict[str, Any]:
        return {
            "inventory": self.inventory.as_posix() if self.inventory else None,
            "settings_overrides": self.settings_overrides.as_posix() if self.settings_overrides else None,
            "lineage_policy": self.lineage_policy, "pool": self.pool,
            "label_keys": self.label_keys.as_dict(), "include_parts": list(self.include_parts),
            "exclude_globs": list(self.exclude_globs), "include_globs": list(self.include_globs),
            "follow_symlinks": self.follow_symlinks, "stride": self.stride,
            "max_frames_per_run": self.max_frames_per_run,
        }


# --------------------------------------------------------------------------
# Per-run and per-frame state
# --------------------------------------------------------------------------

@dataclass
class _Run:
    ordinal: int
    run_id: str
    root_alias: str
    relpath: str
    discovered: disc.DiscoveredRun | None  # None: candidate inside a pruned directory
    record: dict[str, Any]
    outcome: Outcome
    classification: Outcome | None = None
    settings: dict[str, Any] | None = None
    species: list[str] = field(default_factory=list)
    selective: Any = None
    magmom_initial: list[float] | None = None
    spill: Path | None = None
    source: str | None = None
    source_file_type: str | None = None
    source_sha256: str | None = None
    vasp_version: str | None = None
    calc_type: str | None = None
    facts: dict[str, Any] = field(default_factory=dict)
    inventory_entry: disc.InventoryEntry | None = None
    lineage_input: lin.RunLineageInput | None = None
    lineage: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    pool_id: str | None = None
    frames: list["_Frame"] = field(default_factory=list)


@dataclass
class _Frame:
    frame_id: str
    run: _Run
    index: int
    record: dict[str, Any]
    classification: Outcome
    outcome: Outcome
    order_key: str | None = None
    permutation_key: str | None = None
    spill_row: int | None = None
    label_energy: float | None = None
    label_sha256: str | None = None
    duplicate_of: str | None = None
    global_flags: list[str] = field(default_factory=list)

    @property
    def accepted(self) -> bool:
        return self.outcome.status == ACCEPTED


class _Spill:
    """Arrays of label-accepted frames, one ``.npz`` per run (no pickles), read back in run order."""

    def __init__(self, directory: Path):
        self.directory = directory
        self.directory.mkdir(parents=True, exist_ok=True)
        self._path: Path | None = None
        self._data: dict[str, np.ndarray] | None = None

    def write(self, ordinal: int, rows: list[dict[str, Any]], n_atoms: int) -> Path:
        def stack(name, shape, dtype=np.float64):
            present = np.array([rows[i][name] is not None for i in range(len(rows))], dtype=bool)
            values = np.full((len(rows),) + shape, np.nan if dtype == np.float64 else 0, dtype=dtype)
            for i, row in enumerate(rows):
                if row[name] is not None:
                    values[i] = np.asarray(row[name], dtype=dtype)
            return values, present

        arrays: dict[str, np.ndarray] = {"index": np.array([row["index"] for row in rows], dtype=np.int64)}
        arrays["cell"], _ = stack("cell", (3, 3))
        arrays["positions"], _ = stack("positions", (n_atoms, 3))
        arrays["forces"], arrays["has_forces"] = stack("forces", (n_atoms, 3))
        arrays["stress"], arrays["has_stress"] = stack("stress", (3, 3))
        arrays["magmom_final"], arrays["has_magmom_final"] = stack("magmom_final", (n_atoms,))
        arrays["permutation"], _ = stack("permutation", (n_atoms,), np.int64)
        path = self.directory / f"run{ordinal:06d}.npz"
        with path.open("wb") as handle:
            np.savez(handle, **arrays)
        return path

    def load(self, path: Path) -> dict[str, np.ndarray]:
        if self._path != path:
            with np.load(path, allow_pickle=False) as data:
                self._data = {name: data[name] for name in data.files}
            self._path = path
        assert self._data is not None
        return self._data

    def row(self, frame: _Frame) -> dict[str, Any]:
        data = self.load(frame.run.spill)
        i = frame.spill_row
        if int(data["index"][i]) != frame.index:
            raise DatasetError(f"internal error: spilled arrays of {frame.frame_id} are out of order")
        return {
            "cell": data["cell"][i], "positions": data["positions"][i],
            "forces": data["forces"][i] if data["has_forces"][i] else None,
            "stress": data["stress"][i] if data["has_stress"][i] else None,
            "magmom_final": data["magmom_final"][i] if data["has_magmom_final"][i] else None,
            "permutation": data["permutation"][i],
        }


# --------------------------------------------------------------------------
# One run: evidence + one streaming vasprun pass + classification
# --------------------------------------------------------------------------

_EVIDENCE_PARSERS = (
    ("OUTCAR", vaspfiles.parse_outcar), ("OSZICAR", vaspfiles.parse_oszicar),
    ("POTCAR", vaspfiles.parse_potcar), ("INCAR", vaspfiles.parse_incar), ("POSCAR", vaspfiles.parse_poscar),
)


def _portable(text: str, found: disc.DiscoveredRun, alias: str) -> str:
    """Replace absolute paths under the run's root by ``alias:relpath`` (bytes never depend on where a root lives)."""
    root = Path(found.root_path)
    for form in sorted({str(root), root.as_posix()}, key=len, reverse=True):
        text = text.replace(form + os.sep, alias + ":").replace(form + "/", alias + ":").replace(form, alias + ":")
    return text.replace("\\", "/") if alias + ":" in text else text


def _parse_evidence(found: disc.DiscoveredRun, alias: str) -> tuple[dict[str, Any], dict[str, Any]]:
    parsed: dict[str, Any] = {}
    files: dict[str, Any] = {}
    for kind, parser in _EVIDENCE_PARSERS:
        path = found.evidence.get(kind)
        if path is None:
            continue
        try:
            evidence = parser(path)
        except OSError as exc:
            files[kind] = {"file": path.name, "used": False,
                           "error": _portable(f"{type(exc).__name__}: {exc}", found, alias)}
            continue
        parsed[kind] = evidence
        problems = list(getattr(evidence, "problems", []) or [])
        files[kind] = {"file": path.name, "sha256": evidence.sha256, "bytes": path.stat().st_size, "used": True,
                       "problems": problems[:20], "n_problems": len(problems)}
    for kind in sorted(set(found.evidence) - {k for k, _ in _EVIDENCE_PARSERS}):
        path = found.evidence[kind]
        try:
            files[kind] = {"file": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size,
                           "used": False}
        except OSError as exc:
            files[kind] = {"file": path.name, "used": False,
                           "error": _portable(f"{type(exc).__name__}: {exc}", found, alias)}
    return parsed, files


def _order_key(species: Sequence[str], structure: Any) -> dup.StructureKey | None:
    if structure is None or not species:
        return None
    return dup.step_structure_key(species, structure)


def _mag_frame_fields(magnetic: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "magnetic_class": magnetic.get("magnetic_class"),
        "magnetic_policy": magnetic.get("policy_result"),
        "total_magnetization": magnetic.get("total"),
        "magnetization_source": magnetic.get("total_source"),
        "magnetic_state_id": magnetic.get("magnetic_state_id"),
        "magnetic_segment": magnetic.get("segment"),
        "mag_d_total_first": magnetic.get("d_total_first"),
        "mag_d_total_previous": magnetic.get("d_total_previous"),
    }


def _process_run(run: _Run, *, policy: acc.Policy, overrides: st.Overrides, spill: _Spill | None,
                 cache: lin.ManifestCache) -> None:
    """The single streaming pass over one run (fills ``run.record`` and ``run.frames``)."""
    found = run.discovered
    assert found is not None and found.label_file is not None
    evidence, files = _parse_evidence(found, run.root_alias)
    run.record["evidence_files"] = files
    label_path = found.label_file
    run.source = f"{run.root_alias}:{run.relpath}/{label_path.name}" if run.relpath != "." else \
        f"{run.root_alias}:{label_path.name}"
    run.source_file_type = found.label_kind
    run.record.update({"source": run.source, "source_file": label_path.name, "source_file_type": found.label_kind})
    reference_hashes = lin.agglomeration_reference_hashes(found.run_dir, found.root_path, cache)
    common = dict(outcar=evidence.get("OUTCAR"), oszicar=evidence.get("OSZICAR"), potcar=evidence.get("POTCAR"),
                  poscar=evidence.get("POSCAR"), incar=evidence.get("INCAR"), policy=policy,
                  reference_hashes=reference_hashes)
    try:
        reader = vr.VasprunReader(label_path)
    except (vr.VasprunParseError, OSError) as exc:
        assessment = acc.classify_run(run_id=run.run_id, header=None,
                                      header_error=_portable(str(exc), found, run.root_alias), **common)
        try:
            run.source_sha256 = sha256_file(label_path)
        except OSError:
            run.source_sha256 = None
        run.record.update({"source_sha256": run.source_sha256, "source_bytes": _size(label_path),
                           "assessment": assessment.as_dict(), "frames_total": 0})
        run.classification = run.outcome = assessment.outcome
        return
    header = reader.header
    steps = list(reader.steps())  # the one pass; the reader closes itself at the end
    trailer = reader.trailer()
    settings = st.extract_settings(header, incar_file=evidence.get("INCAR"), outcar=evidence.get("OUTCAR"),
                                   potcar=evidence.get("POTCAR"))
    assessment = acc.classify_run(run_id=run.run_id, header=header, trailer=trailer, steps=steps,
                                  settings=settings, **common)
    step_results = acc.classify_steps(steps, assessment)
    species = list(assessment.species)
    run.species = species
    run.settings = settings
    run.selective = assessment.selective
    run.facts = dict(assessment.facts)
    run.calc_type = assessment.calc_type
    run.vasp_version = assessment.facts.get("version") or (header.generator or {}).get("version") or None
    run.source_sha256 = trailer.source_sha256
    run.classification = run.outcome = assessment.outcome
    n_atoms = len(species)
    magmom = assessment.facts.get("magmom_parameters")
    if assessment.facts.get("ispin") == 2 and magmom is not None and len(magmom) == n_atoms:
        run.magmom_initial = [float(v) for v in magmom]
    composition, elements = _composition(species) if species else (None, None)
    initial = _order_key(species, header.initial_structure)
    final = _order_key(species, trailer.final_structure)
    frame_keys: list[str] = []
    spill_rows: list[dict[str, Any]] = []
    for step, result in zip(steps, step_results):
        key = _order_key(species, step.structure)
        if key is not None:
            frame_keys.append(key.order_key)
        frame_id = f"{run.run_id}#{step.index:05d}"
        value, energy_source, energy_rule = vr.label_energy(result.energy, policy.energy_quantity)
        data = result.as_dict()
        stress_arr = result.stress.stress if result.stress.available else None
        record = {
            "frame_id": frame_id, "run_id": run.run_id, "index": result.index,
            "structure_key": key.order_key if key else None,
            "permutation_key": key.permutation_key if key else None,
            "classification": result.outcome.as_dict(),
            "label_quantity": policy.energy_quantity,
            "label_energy_source": energy_source, "label_energy_rule": energy_rule,
            "composition": composition, "elements": elements,
            **_mag_frame_fields(result.magnetic),
            **{k: v for k, v in data.items() if k not in ("outcome", "index")},
        }
        frame = _Frame(frame_id=frame_id, run=run, index=result.index, record=record,
                       classification=result.outcome, outcome=result.outcome,
                       order_key=key.order_key if key else None,
                       permutation_key=key.permutation_key if key else None)
        if result.accepted:
            forces = None if step.forces is None else np.asarray(step.forces, dtype=np.float64)
            frame.label_energy = result.label_energy
            frame.label_sha256 = xyz.label_sha256(result.label_energy, forces, stress_arr)
            record["label_sha256"] = frame.label_sha256
            if spill is not None:
                frame.spill_row = len(spill_rows)
                site = result.magnetic.get("site_moments")
                spill_rows.append({
                    "index": result.index,
                    "cell": np.asarray(step.structure.cell, dtype=np.float64),
                    "positions": np.asarray(step.structure.positions, dtype=np.float64),
                    "forces": forces, "stress": stress_arr,
                    "magmom_final": site if site is not None and len(site) == n_atoms else None,
                    "permutation": np.asarray(key.permutation, dtype=np.int64) if key else None,
                })
        run.frames.append(frame)
    if spill is not None and spill_rows:
        if any(row["permutation"] is None for row in spill_rows):
            raise DatasetError(f"internal error: an accepted frame of {run.run_id} has no structure key")
        run.spill = spill.write(run.ordinal, spill_rows, n_atoms)
    header_record = {
        "generator": dict(header.generator or {}), "vasp_version": run.vasp_version,
        "calc_type": assessment.calc_type, "calc_detail": assessment.calc_detail,
        "n_atoms": n_atoms or None, "species_blocks": _species_blocks(species),
        "composition": composition, "elements": elements,
        "source_sha256": trailer.source_sha256, "source_bytes": trailer.source_bytes,
        "compression": trailer.compression,
        "trailer": {"closed": trailer.closed, "truncated": trailer.truncated, "error": trailer.error,
                    "steps_seen": trailer.steps_seen, "n_complete_steps": trailer.n_complete_steps,
                    "final_structure_present": trailer.final_structure_present,
                    "partial_tail_tag": trailer.partial_tail_tag, "problems": list(trailer.problems)[:20]},
        "header_problems": list(header.problems)[:20],
        "settings": settings,
        "method_fingerprint": st.method_fingerprint(settings, overrides),
        "assessment": assessment.as_dict(),
        "magmom_initial": run.magmom_initial,
        "frames_total": len(run.frames),
        "initial_structure_key": initial.order_key if initial else None,
        "final_structure_key": final.order_key if final else None,
    }
    run.record.update(header_record)
    run.lineage_input = lin.RunLineageInput(
        run_id=run.run_id, root_alias=run.root_alias, relpath=run.relpath, root_path=found.root_path,
        run_dir=found.run_dir, nested_in=found.nested_in,
        inventory_lineage=run.inventory_entry.lineage if run.inventory_entry else None,
        inventory_metadata=dict(run.inventory_entry.metadata) if run.inventory_entry else {},
        initial_key=initial.order_key if initial else None, frame_keys=tuple(frame_keys),
        final_key=final.order_key if final else None,
    )


def _first_number(text: str | None) -> float | None:
    if text is None:
        return None
    token = str(text).replace(",", " ").split()
    if not token:
        return None
    try:
        return float(token[0].replace("d", "e").replace("D", "e"))
    except ValueError:
        return None


def _labelless_magnetic_evidence(run: _Run, policy: acc.Policy) -> dict[str, Any] | None:
    """Magnetic class of a run WITHOUT a label file, from INCAR (+ POSCAR/CONTCAR species) alone.

    Evidence only: the run keeps its missing-labels outcome; the class shows what the run would
    be if its labels were recovered (OutPackLite trees: ISPIN=2 without MAGMOM -> uncontrolled).
    """
    found = run.discovered
    if found is None or found.evidence.get("INCAR") is None:
        return None
    try:
        incar = vaspfiles.parse_incar(found.evidence["INCAR"])
    except OSError:
        return None
    tags = {str(k).upper(): v for k, v in (incar.tags or {}).items()}
    species: list[str] = []
    species_source = None
    for kind in ("POSCAR", "CONTCAR"):
        path = found.evidence.get(kind)
        if path is None:
            continue
        try:
            poscar = vaspfiles.parse_poscar(path)
        except (OSError, ValueError, DatasetError):
            continue
        if poscar.species and poscar.counts and len(poscar.species) == len(poscar.counts):
            species = [symbol for symbol, count in zip(poscar.species, poscar.counts) for _ in range(int(count))]
            species_source = kind
            break
    ispin_value = _first_number(tags.get("ISPIN"))
    ispin = int(ispin_value) if ispin_value is not None else 1  # VASP default ISPIN=1
    nupdown = _first_number(tags.get("NUPDOWN"))
    noncollinear = str(tags.get("LNONCOLLINEAR", "")).strip().upper().lstrip(".").startswith("T") or None
    magnetic = acc.classify_magnetism(
        species=species, ispin=ispin, nupdown=nupdown, noncollinear=noncollinear, magmom_initial=None,
        magmom_explicit="MAGMOM" in tags, dft_indices=[], totals={}, sites={}, policy=policy,
        magmom_source="incar_file",
    ).run
    would_be = None
    if magnetic.get("magnetic_species"):
        if magnetic["magnetic_class"] == "uncontrolled" and not policy.accept_uncontrolled:
            would_be = "magnetic_uncontrolled"
        elif magnetic["magnetic_class"] == "unknown" and not policy.accept_unknown:
            would_be = "magnetic_unknown"
    return st.json_safe({
        "evidence_only": True, "source": "INCAR" + (f"+{species_source}" if species_source else ""),
        "magnetic_class": magnetic.get("magnetic_class"), "class_reason": magnetic.get("class_reason"),
        "ispin": ispin, "ispin_source": "INCAR" if ispin_value is not None else "VASP default",
        "magmom_explicit": "MAGMOM" in tags, "nupdown": nupdown,
        "magnetic_species": magnetic.get("magnetic_species"), "n_atoms": len(species) or None,
        "run_reason_if_labelled": would_be,
    })


def _size(path: Path) -> int | None:
    try:
        return path.stat().st_size
    except OSError:
        return None


# --------------------------------------------------------------------------
# Analysis (shared by scan and export)
# --------------------------------------------------------------------------

@dataclass
class Analysis:
    roots: list[disc.RootSpec]
    options: ExportOptions
    runs: list[_Run]
    ignored: list[dict[str, Any]]
    discovery: dict[str, Any]
    pools: list[st.Pool]
    pool_frames: dict[str, int]
    selected_pool: st.Pool | None
    pool_error: str | None
    duplicates: dup.DuplicateReport
    near_duplicates: list[dict[str, Any]]
    lineage: lin.LineageResolution
    inventory: disc.Inventory | None
    inventory_unmatched: list[dict[str, Any]]
    overrides: st.Overrides
    git: dict[str, Any]
    spill: _Spill | None

    @property
    def frames(self) -> Iterator[_Frame]:
        for run in self.runs:
            yield from run.frames

    def exported(self) -> list[_Frame]:
        return [frame for frame in self.frames if frame.accepted]


def _pool_error(pools: list[st.Pool], pool_frames: Mapping[str, int]) -> str:
    ordered = sorted(pools, key=lambda p: (-pool_frames.get(p.pool_id, 0), p.pool_id))
    lines = [f"{len(pools)} reference-settings pools among accepted runs; export one with --pool ID "
             "(or reconcile settings with --settings-overrides):"]
    reference = ordered[0]
    for pool in ordered:
        lines.append(f"  {pool.pool_id}: {len(pool.run_ids)} run(s), {pool_frames.get(pool.pool_id, 0)} accepted "
                     f"frame(s), e.g. {pool.run_ids[0]}")
        if pool is not reference:
            diffs = [d for d in st.describe_pool_differences(reference, pool) if d["kind"] == "conflict"]
            text = "; ".join(f"{d['field']}: {d['a']!r} vs {d['b']!r}" for d in diffs[:6])
            lines.append(f"    differs from {reference.pool_id}: {text or 'element coverage only'}")
    return "\n".join(lines)


def _select_pool(pools: list[st.Pool], wanted: str | None) -> st.Pool | None:
    if wanted is None:
        return pools[0] if len(pools) == 1 else None
    wanted = wanted.strip().lower()
    if len(wanted) < _MIN_POOL_PREFIX:
        raise DatasetError(f"--pool needs at least {_MIN_POOL_PREFIX} characters of a pool id, got {wanted!r}")
    matches = [pool for pool in pools if pool.pool_id.startswith(wanted)]
    if len(matches) != 1:
        known = ", ".join(pool.pool_id for pool in pools) or "none"
        raise DatasetError(f"--pool {wanted!r} matches {len(matches)} of the pools among accepted runs ({known})")
    return matches[0]


def _metadata(run: _Run) -> dict[str, Any]:
    """campaign/family/config_type/temperature for records and extxyz info (inventory first)."""
    entry = run.inventory_entry
    meta = dict(run.lineage.get("metadata") or {})
    inv_meta = dict(entry.metadata) if entry else {}

    def first(*values):
        for value in values:
            if value is not None and str(value).strip():
                return value
        return None

    family = first(entry.family if entry else None, inv_meta.get("family"), meta.get("family"),
                   meta.get("structure_family"))
    campaign = first(entry.campaign if entry else None, inv_meta.get("campaign"), meta.get("campaign"))
    config_type = first(entry.config_type if entry else None, inv_meta.get("config_type"), family)
    temperature = first(inv_meta.get("temperature_K"), meta.get("temperature_K"), meta.get("temperature_k"))
    temperature_source = "inventory" if inv_meta.get("temperature_K") is not None else (
        "lineage" if temperature is not None else None)
    if temperature is None and run.calc_type == "md":
        tebeg, teend = run.facts.get("tebeg"), run.facts.get("teend")
        if tebeg is not None and (teend is None or teend < 0 or teend == tebeg):
            temperature, temperature_source = float(tebeg), "TEBEG"
        elif tebeg is not None and teend is not None:
            temperature, temperature_source = f"{tebeg:g}->{teend:g}", "TEBEG->TEEND"
    if isinstance(temperature, (int, float)) and not isinstance(temperature, bool):
        temperature = float(temperature)
    result = {
        "campaign": str(campaign) if campaign is not None else UNKNOWN_TEXT,
        "family": str(family) if family is not None else UNKNOWN_TEXT,
        "config_type": str(config_type) if config_type is not None else UNKNOWN_TEXT,
        "temperature_K": temperature, "temperature_source": temperature_source,
    }
    for key in ("campaign_stage", "replica", "stage", "structure_id", "family_source", "fam_coverage",
                "fam_arrangement", "fam_motif", "fam_ligand", "fam_anchor", "fam_parent"):
        if meta.get(key) is not None:
            result[key] = meta[key]
    return result


def analyze(roots: Sequence[Any], options: ExportOptions, *, output_dir: Path | None,
            spill_dir: Path | None) -> Analysis:
    """Discovery -> per-run pass -> pools -> duplicates -> lineage -> subsampling (no files written)."""
    specs = disc.parse_roots(roots)
    policy = options.policy
    overrides = st.load_overrides(options.settings_overrides)
    inventory = disc.load_inventory(options.inventory) if options.inventory is not None else None
    found, ignored = disc.discover_runs(
        specs, exclude_globs=options.exclude_globs, include_globs=options.include_globs,
        include_parts=options.include_parts, output_dir=output_dir, follow_symlinks=options.follow_symlinks,
    )
    decisions: dict[str, disc.InventoryDecision] = {}
    unmatched: list[dict[str, Any]] = []
    if inventory is not None:
        resolution = inventory.resolve(found, ignored=ignored, strict=True)
        decisions, unmatched = resolution.decisions, resolution.unmatched
    spill = _Spill(spill_dir) if spill_dir is not None else None
    cache = lin.ManifestCache()

    runs: list[_Run] = []
    for item in found:
        record: dict[str, Any] = {
            "run_id": item.run_id, "root_alias": item.root_alias, "relpath": item.relpath,
            "discovery": {k: v for k, v in item.as_dict().items() if k not in ("run_id", "root_alias", "relpath")},
            "parsed": False,
        }
        run = _Run(ordinal=0, run_id=item.run_id, root_alias=item.root_alias, relpath=item.relpath,
                   discovered=item, record=record, outcome=Outcome(ACCEPTED))
        decision = decisions.get(item.run_id)
        if decision is not None:
            run.inventory_entry = decision.entry
            record["inventory"] = {"selected": decision.selected, "reason": decision.reason,
                                   "entry": decision.entry.as_dict() if decision.entry else None}
        runs.append(run)
    for candidate in disc.pruned_run_candidates(ignored):
        record = {"run_id": candidate["run_id"], "root_alias": candidate["root_alias"],
                  "relpath": candidate["relpath"], "parsed": False,
                  "discovery": {"pruned": True, "discovery_reason": candidate["discovery_reason"],
                                "rule": candidate["rule"]}}
        outcome = Outcome(**candidate["outcome"])
        runs.append(_Run(ordinal=0, run_id=candidate["run_id"], root_alias=candidate["root_alias"],
                         relpath=candidate["relpath"], discovered=None, record=record, outcome=outcome,
                         classification=outcome))
    runs.sort(key=lambda r: (r.root_alias, r.relpath, r.discovered is None))
    seen: set[str] = set()
    for ordinal, run in enumerate(runs):
        if run.run_id in seen:
            raise DatasetError(f"internal error: run id {run.run_id!r} discovered twice")
        seen.add(run.run_id)
        run.ordinal = ordinal

    # 1. one pass per run -----------------------------------------------------------------
    for run in runs:
        found_run = run.discovered
        if found_run is None:
            continue
        decision = decisions.get(run.run_id)
        if decision is not None and not decision.selected:
            run.classification = run.outcome = outcome_for(
                decision.reason, f"inventory {inventory.path.name}: " + (
                    "entry has include = false" if decision.reason == "inventory_excluded" else "run not listed"))
        elif found_run.label_file is None:
            run.classification = run.outcome = found_run.missing_label_outcome()
            evidence = _labelless_magnetic_evidence(run, policy)
            if evidence is not None:
                run.record["magnetic_evidence"] = evidence
        else:
            run.record["parsed"] = True
            _process_run(run, policy=policy, overrides=overrides, spill=spill, cache=cache)
        if run.lineage_input is None:
            run.lineage_input = lin.RunLineageInput(
                run_id=run.run_id, root_alias=run.root_alias, relpath=run.relpath,
                root_path=found_run.root_path, run_dir=found_run.run_dir, nested_in=found_run.nested_in,
                inventory_lineage=run.inventory_entry.lineage if run.inventory_entry else None,
                inventory_metadata=dict(run.inventory_entry.metadata) if run.inventory_entry else {},
            )

    # 2. reference-settings pools over runs with label-accepted frames -----------------------
    candidates = [run for run in runs if run.settings is not None and any(f.accepted for f in run.frames)]
    pools = st.build_pools([(run.run_id, run.settings) for run in candidates], overrides=overrides)
    pool_of = {run_id: pool.pool_id for pool in pools for run_id in pool.run_ids}
    pool_frames: dict[str, int] = {}
    for run in runs:
        run.pool_id = pool_of.get(run.run_id)
        run.record["pool_id"] = run.pool_id
        if run.pool_id is not None:
            pool_frames[run.pool_id] = pool_frames.get(run.pool_id, 0) + sum(f.accepted for f in run.frames)
    selected = _select_pool(pools, options.pool)
    pool_error = _pool_error(pools, pool_frames) if selected is None and len(pools) > 1 else None

    # 3. exact duplicates (within a pool) + shared-structure links -----------------------------
    by_key: dict[str, list[_Frame]] = {}
    for run in runs:
        for frame in run.frames:
            if frame.accepted and frame.permutation_key and run.pool_id is not None:
                by_key.setdefault(frame.permutation_key, []).append(frame)
    dup_candidates = []
    tolerances = {"energy": policy.duplicate_energy_tolerance, "force": policy.duplicate_force_tolerance,
                  "magnetization": policy.duplicate_magnetization_tolerance,
                  "site_sign_threshold": policy.magmom_sign_threshold}
    for key in sorted(by_key):
        members = by_key[key]
        if len(members) < 2:
            continue
        for frame in members:
            arrays = spill.row(frame) if spill is not None else None
            dup_candidates.append(dup.DuplicateCandidate(
                frame_id=frame.frame_id, run_id=frame.run.run_id, pool_id=frame.run.pool_id,
                permutation_key=frame.permutation_key, order_key=frame.order_key,
                permutation=[int(i) for i in arrays["permutation"]] if arrays is not None else None,
                energy=frame.label_energy,
                forces=arrays["forces"] if arrays is not None else None,
                total_magnetization=frame.record.get("total_magnetization"),
                site_magmoms=(frame.record.get("magnetic") or {}).get("site_moments"),
                magnetic_state_id=frame.record.get("magnetic_state_id"),
            ))
    if spill is None and dup_candidates:
        raise DatasetError("internal error: duplicate search needs the spilled arrays")
    report = dup.find_exact_duplicates(dup_candidates, tolerances)
    shared = dup.structure_links({"permutation_key": f.permutation_key, "run_id": f.run.run_id}
                                 for run in runs for f in run.frames if f.permutation_key)
    links = list(report.links) + shared

    # 4. frame outcomes decided by pool selection and duplicates --------------------------------
    for run in runs:
        not_selected = selected is not None and run.pool_id is not None and run.pool_id != selected.pool_id
        if not_selected and run.outcome.status == ACCEPTED:
            run.outcome = outcome_for("reference_pool_not_selected",
                                      f"pool {run.pool_id}; exported pool {selected.pool_id}")
        for frame in run.frames:
            if not frame.accepted:
                continue
            if not_selected:
                frame.outcome = outcome_for("reference_pool_not_selected",
                                            f"pool {run.pool_id}; exported pool {selected.pool_id}")
            elif frame.frame_id in report.outcomes:
                frame.outcome = report.outcomes[frame.frame_id]
                frame.duplicate_of = report.duplicate_of.get(frame.frame_id)
            frame.global_flags.extend(report.flags.get(frame.frame_id, []))

    # 5. lineage (+ near-duplicates) -------------------------------------------------------------
    inputs = [run.lineage_input for run in runs if run.lineage_input is not None]
    resolution = lin.resolve_lineage(inputs, policy=options.lineage_policy, extra_links=links, cache=cache)
    near: list[dict[str, Any]] = []
    if policy.near_duplicate_tolerance is not None:
        near_candidates = []
        for run in runs:
            for frame in run.frames:
                if frame.accepted:
                    arrays = spill.row(frame) if spill is not None else None
                    if arrays is None:
                        raise DatasetError("internal error: near-duplicate search needs the spilled arrays")
                    near_candidates.append({
                        "frame_id": frame.frame_id, "group_id": resolution.group_of(run.run_id),
                        "species": run.species, "cell": arrays["cell"], "positions": arrays["positions"],
                        "run_id": run.run_id})
        near = dup.find_near_duplicates(near_candidates, policy.near_duplicate_tolerance)
        if near:
            near_links = [{"a": p["run_a"], "b": p["run_b"], "kind": "near_duplicate"} for p in near]
            resolution = lin.resolve_lineage(inputs, policy=options.lineage_policy, extra_links=links + near_links,
                                             cache=cache)
    for run in runs:
        record = resolution.runs.get(run.run_id)
        if record is None:
            run.lineage = {"lineage_id": None, "source": "not_applicable", "group_id": None, "lineage_group": None}
        else:
            run.lineage = dict(record)
        run.metadata = _metadata(run)
        unresolved = record is not None and record["source"] == "unresolved"
        if unresolved and run.outcome.status == ACCEPTED and any(f.accepted for f in run.frames):
            run.outcome = outcome_for("lineage_unresolved", "; ".join(
                [str(record["evidence"].get("why") or "")] + list(record["evidence"].get("notes") or [])).strip("; ")
                + " -- declare --lineage-policy run|parent|depth:N or an inventory lineage")
        for frame in run.frames:
            if frame.accepted and unresolved:
                frame.outcome = outcome_for("lineage_unresolved", f"run {run.run_id} has no lineage source")

    # 6. subsampling (declared, deterministic) -----------------------------------------------------
    if options.stride > 1 or options.max_frames_per_run is not None:
        for run in runs:
            kept = [frame for frame in run.frames if frame.accepted]
            if options.stride > 1:
                for position, frame in enumerate(kept):
                    if position % options.stride:
                        frame.outcome = outcome_for("subsampled", f"stride {options.stride}")
                kept = kept[::options.stride]
            limit = options.max_frames_per_run
            if limit is not None and len(kept) > limit:
                chosen = {int(round(x)) for x in np.linspace(0, len(kept) - 1, limit)}
                for position, frame in enumerate(kept):
                    if position not in chosen:
                        frame.outcome = outcome_for("subsampled", f"max_frames_per_run {limit}")

    return Analysis(
        roots=specs, options=options, runs=runs, ignored=ignored,
        discovery=disc.discovery_summary(found, ignored), pools=pools, pool_frames=pool_frames,
        selected_pool=selected, pool_error=pool_error, duplicates=report, near_duplicates=near,
        lineage=resolution, inventory=inventory, inventory_unmatched=unmatched, overrides=overrides,
        git=git_provenance(), spill=spill,
    )


# --------------------------------------------------------------------------
# Records
# --------------------------------------------------------------------------

def _frame_record(frame: _Frame) -> dict[str, Any]:
    run = frame.run
    record = dict(frame.record)
    record.update({
        "outcome": frame.outcome.as_dict(),
        "exported": frame.accepted,
        "pool_id": run.pool_id,
        "lineage_id": run.lineage.get("lineage_id"),
        "lineage_source": run.lineage.get("source"),
        "lineage_group": run.lineage.get("lineage_group"),
        "duplicate_of": frame.duplicate_of,
        "campaign": run.metadata.get("campaign"),
        "family": run.metadata.get("family"),
        "config_type": run.metadata.get("config_type"),
        "temperature_K": run.metadata.get("temperature_K"),
        "vasp_version": run.vasp_version,
        "calc_type": run.calc_type,
        "source": run.source,
        "source_file_type": run.source_file_type,
        "source_sha256": run.source_sha256,
    })
    if frame.global_flags:
        record["flags"] = sorted(set(list(record.get("flags") or []) + frame.global_flags))
    return st.json_safe(record)


def _run_record(run: _Run) -> dict[str, Any]:
    record = dict(run.record)
    statuses = _count(frame.outcome.status for frame in run.frames)
    record.update({
        "outcome": run.outcome.as_dict(),
        "classification": run.classification.as_dict() if run.classification is not None else None,
        "pool_id": run.pool_id,
        "lineage": {k: run.lineage.get(k) for k in ("lineage_id", "source", "evidence", "links")}
        if run.lineage else None,
        "lineage_group": run.lineage.get("lineage_group") if run.lineage else None,
        "metadata": run.metadata or None,
        "frames_total": len(run.frames) if run.record.get("parsed") else run.record.get("frames_total"),
        "frames_by_status": statuses,
        "frames_exported": sum(1 for frame in run.frames if frame.accepted),
    })
    return st.json_safe(record)


def _exclusion_rows(analysis: Analysis) -> list[dict[str, str]]:
    rows = []
    for run in analysis.runs:
        if run.outcome.status != ACCEPTED:
            rows.append({"level": "run", "id": run.run_id, "run_id": run.run_id, "status": run.outcome.status,
                         "reason": run.outcome.reason or "", "detail": run.outcome.detail})
        for frame in run.frames:
            if not frame.accepted:
                rows.append({"level": "frame", "id": frame.frame_id, "run_id": run.run_id,
                             "status": frame.outcome.status, "reason": frame.outcome.reason or "",
                             "detail": frame.outcome.detail})
    for item in analysis.ignored:
        if item.get("kind") == "run_dir":
            continue  # already a run record
        rows.append({"level": f"discovery:{item.get('kind')}", "id": f"{item['root_alias']}:{item['relpath']}",
                     "run_id": "", "status": "ignored", "reason": item.get("reason") or "",
                     "detail": " ".join(str(x) for x in (item.get("rule"), item.get("detail")) if x)})
    return rows


def _csv_text(rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=EXCLUSION_COLUMNS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: str(row.get(key, "")).replace("\r", " ").replace("\n", " ") for key in EXCLUSION_COLUMNS})
    return buffer.getvalue()


# --------------------------------------------------------------------------
# extxyz payloads
# --------------------------------------------------------------------------

def _repo_commit(git: Mapping[str, Any]) -> str:
    commit = git.get("commit")
    if not commit:
        return UNKNOWN_TEXT
    return commit + ("-dirty" if git.get("working_tree_dirty") else "")


def _payloads(analysis: Analysis, frames: Sequence[_Frame]) -> Iterator[xyz.FramePayload]:
    policy = analysis.options.policy
    commit = _repo_commit(analysis.git)
    assert analysis.spill is not None
    for frame in frames:
        run = frame.run
        rec = frame.record
        arrays = analysis.spill.row(frame)
        energies = rec.get("energy") or {}
        scf = rec.get("scf") or {}
        vacuum = rec.get("vacuum_axes") or []
        info: dict[str, Any] = {
            "source": run.source, "source_file_type": run.source_file_type, "ionic_step": frame.index,
            "structure_key": frame.order_key, "lineage_group": run.lineage.get("lineage_group"),
            "energy_source": rec.get("label_energy_source"), "energy_rule": rec.get("label_energy_rule"),
            "stress_available": arrays["stress"] is not None, "scf_status": scf.get("status"),
            "label_source": rec.get("label_source"), "vasp_version": run.vasp_version or UNKNOWN_TEXT,
            "magnetic_class": rec.get("magnetic_class") or UNKNOWN_TEXT,
            "magnetic_policy": rec.get("magnetic_policy") or UNKNOWN_TEXT,
            "campaign": run.metadata["campaign"], "family": run.metadata["family"],
            "parser": PARSER_NAME, "parser_version": PARSER_VERSION, "repo_commit": commit,
            # provenance (informative)
            "run_id": run.run_id, "source_sha256": run.source_sha256, "calc_type": run.calc_type,
            "pool_id": run.pool_id, "lineage_id": run.lineage.get("lineage_id"),
            "lineage_source": run.lineage.get("source"), "config_type": run.metadata["config_type"],
            "temperature_K": run.metadata.get("temperature_K"), "permutation_key": frame.permutation_key,
            "energy_quantity": policy.energy_quantity, "energy_version_gate": energies.get("energy_version_gate"),
            "vasp_free_energy": energies.get("free_energy"),
            "vasp_energy_sigma0": energies.get("energy_sigma0"),
            "vasp_energy_no_entropy": energies.get("energy_no_entropy"),
            "time_fs": rec.get("time_fs"), "n_scf_steps": scf.get("n_steps"), "scf_evidence": scf.get("evidence"),
            "total_magnetization": rec.get("total_magnetization"),
            "magnetization_source": rec.get("magnetization_source"),
            "magnetic_state_id": rec.get("magnetic_state_id"), "magnetic_segment": rec.get("magnetic_segment"),
            "mag_d_total_first": rec.get("mag_d_total_first"),
            "mag_d_total_previous": rec.get("mag_d_total_previous"),
            "vacuum_axes": ",".join(str(a) for a in vacuum) if vacuum else "none",
            "campaign_stage": run.metadata.get("campaign_stage"),
        }
        if arrays["stress"] is not None:
            info["stress_source"] = rec.get("stress_source")
        else:
            info["stress_reason"] = rec.get("stress_reason") or UNKNOWN_TEXT
        if arrays["forces"] is None:
            info["forces_reason"] = rec.get("forces_reason") or UNKNOWN_TEXT
        md_temperature = (rec.get("md") or {}).get("md_temperature_K")
        if md_temperature is not None:
            info["md_temperature_K"] = md_temperature
        flags = rec.get("flags") or []
        if flags:
            info["frame_flags"] = ";".join(sorted(set(flags)))
        yield xyz.FramePayload(
            frame_id=frame.frame_id, species=run.species, cell=arrays["cell"], positions=arrays["positions"],
            energy=frame.label_energy, forces=arrays["forces"], stress=arrays["stress"], info=info,
            selective_dynamics=run.selective, magmom_initial=run.magmom_initial,
            magmom_final=arrays["magmom_final"],
        )


# --------------------------------------------------------------------------
# Leakage self-check and statistics
# --------------------------------------------------------------------------

def check_group_consistency(frames: Iterable[_Frame]) -> dict[str, int]:
    """Every structure hash among exported frames maps to exactly one lineage group."""
    owners: dict[tuple[str, str], str] = {}
    groups: dict[str, int] = {}
    for frame in frames:
        group = frame.run.lineage.get("lineage_group")
        if not group or str(group).startswith("unresolved:"):
            raise LeakageError(f"exported frame {frame.frame_id} has no resolved lineage group")
        groups[group] = groups.get(group, 0) + 1
        for kind, key in (("structure_key", frame.order_key), ("permutation_key", frame.permutation_key)):
            if not key:
                raise LeakageError(f"exported frame {frame.frame_id} has no {kind}")
            other = owners.setdefault((kind, key), group)
            if other != group:
                raise LeakageError(f"{kind} {key[:12]} occurs in lineage groups {other!r} and {group!r}")
    return dict(sorted(groups.items()))


def _pool_summary(analysis: Analysis) -> list[dict[str, Any]]:
    result = []
    for pool in analysis.pools:
        result.append({
            "pool_id": pool.pool_id, "runs": len(pool.run_ids), "run_ids": list(pool.run_ids),
            "accepted_frames": analysis.pool_frames.get(pool.pool_id, 0),
            "selected": analysis.selected_pool is not None and pool.pool_id == analysis.selected_pool.pool_id,
            "fingerprints": list(pool.fingerprints), "unknown_fields": st.pool_unknown_fields(pool.method),
            "method": pool.method,
        })
    return result


def _counts(analysis: Analysis) -> dict[str, Any]:
    frames = list(analysis.frames)
    return {
        "runs": len(analysis.runs),
        "runs_parsed": sum(1 for run in analysis.runs if run.record.get("parsed")),
        "runs_by_status": _count(run.outcome.status for run in analysis.runs),
        "runs_by_reason": _count(f"{run.outcome.status}/{run.outcome.reason or 'accepted'}" for run in analysis.runs),
        "frames": len(frames),
        "frames_by_status": _count(frame.outcome.status for frame in frames),
        "frames_by_reason": _count(f"{frame.outcome.status}/{frame.outcome.reason or 'accepted'}" for frame in frames),
        "frames_exported": sum(1 for frame in frames if frame.accepted),
        "frames_label_accepted": sum(1 for frame in frames if frame.classification.status == ACCEPTED),
    }


def _stress_policy(policy: acc.Policy) -> dict[str, Any]:
    return {"include_stress": policy.include_stress, "allow_pstress": policy.allow_pstress,
            "allow_vacuum_stress": policy.allow_vacuum_stress, "vacuum_gap_threshold_A": policy.vacuum_gap_threshold,
            "source": acc.STRESS_SOURCE, "missing": "omitted with stress_available=False and stress_reason (never zeros)",
            "reasons": dict(acc.STRESS_REASONS)}


def _manifest_base(analysis: Analysis, kind: str) -> dict[str, Any]:
    options = analysis.options
    policy = options.policy
    groups = {}
    for frame in analysis.frames:
        if frame.accepted:
            group = frame.run.lineage.get("lineage_group")
            groups[group] = groups.get(group, 0) + 1
    links_by_kind = _count(link.get("kind") for link in analysis.lineage.links if not link.get("ignored"))
    lineage_sources = _count(run.lineage.get("source") for run in analysis.runs if run.lineage)
    return {
        "schema": EXPORT_SCHEMA if kind == "export" else SCAN_SCHEMA,
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "dataset_schema_version": SCHEMA_VERSION,
        "kind": kind,
        "tool": {"name": "nio-md-prep", "version": _tool_version(), "git": analysis.git},
        "parser": {"name": PARSER_NAME, "version": PARSER_VERSION},
        "options": options.as_dict(),
        "policy": policy.as_dict(),
        "roots": [spec.as_dict() for spec in sorted(analysis.roots, key=lambda s: s.alias)],
        "discovery_rules": disc.exclusion_rules_manifest(options.exclude_globs, options.include_globs,
                                                         options.include_parts),
        "discovery": analysis.discovery,
        "inventory": dict(analysis.inventory.as_dict(), unmatched=analysis.inventory_unmatched)
        if analysis.inventory is not None else None,
        "settings_overrides": analysis.overrides.as_dict() if analysis.overrides.entries or
        analysis.overrides.path else None,
        "label_keys": options.label_keys.as_dict(),
        "units": UNITS,
        "energy_quantity": policy.energy_quantity,
        "energy_convention": "F = vasprun <calculation> e_fr_energy - PSTRESS*V (force-consistent); "
                             "E0/E_wo reconstructed from the last <scstep> for VASP <= 6.0.8/unknown/unverified",
        "stress_policy": _stress_policy(policy),
        "frame_contract": xyz.FrameContract(allow_energy_only=policy.allow_energy_only).as_dict(),
        "pool": next((p for p in _pool_summary(analysis) if p["selected"]), None),
        "pools": _pool_summary(analysis),
        "pool_selection_required": analysis.pool_error is not None,
        "counts": _counts(analysis),
        "groups": {"exported_groups": len(groups), "frames_by_group": dict(sorted(groups.items()))},
        "lineage": {"policy": options.lineage_policy, "sources": lineage_sources,
                    "unresolved_runs": analysis.lineage.unresolved, "links_by_kind": links_by_kind,
                    "manifest_errors": analysis.lineage.as_dict()["manifest_errors"]},
        "duplicates": {**analysis.duplicates.as_dict()["counts"], "near_duplicate_pairs": len(analysis.near_duplicates),
                       "near_duplicates": analysis.near_duplicates[:1000]},
        "limitations": LIMITATIONS,
    }


LIMITATIONS = (
    "OUTCAR-only runs (no vasprun.xml) are recorded as missing_labels/no_vasprun; OUTCAR label extraction is not implemented",
    "InterfaceForge repair-prefix recovery (--include-repair-prefix) is not implemented; archived crash states stay excluded",
    "real NiO force validation is outstanding: local OutPackLite trees contain no vasprun.xml/OUTCAR labels",
    "legacy local NiO runs without explicit MAGMOM are quarantined as magnetically uncontrolled unless overridden",
)


def _content_sha256(manifest: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_json({k: v for k, v in manifest.items() if k not in NON_CONTENT_KEYS}).encode())


def _invocation(argv: Sequence[str] | None) -> dict[str, Any]:
    return {"argv": list(argv) if argv is not None else None,
            "created_at": datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()}


# --------------------------------------------------------------------------
# Public entry points
# --------------------------------------------------------------------------

@dataclass
class ExportResult:
    output: Path | None
    manifest: dict[str, Any]
    analysis: Analysis
    files: dict[str, Any]

    @property
    def counts(self) -> dict[str, Any]:
        return self.manifest["counts"]


def _write_accounting(analysis: Analysis, staging: Path) -> dict[str, dict[str, Any]]:
    from . import audit

    run_records = [_run_record(run) for run in analysis.runs]
    frame_records = [_frame_record(frame) for frame in analysis.frames]
    write_jsonl(staging / RUNS_FILE, run_records)
    write_jsonl(staging / FRAMES_FILE, frame_records)
    atomic_write_text(staging / EXCLUSIONS_FILE, _csv_text(_exclusion_rows(analysis)))
    summary = audit.summarize(run_records, frame_records, pools=_pool_summary(analysis),
                              lineage=analysis.lineage.as_dict()["counts"], ignored=analysis.ignored)
    atomic_write_json(staging / AUDIT_JSON, summary)
    atomic_write_text(staging / AUDIT_MD, audit.render_markdown(summary, title="Dataset accounting"))
    return {"summary": summary}


def _file_entry(path: Path, **extra: Any) -> dict[str, Any]:
    return {"sha256": sha256_file(path), "bytes": path.stat().st_size, **extra}


def _publish(staging: Path, output: Path, names: Sequence[str], manifest_name: str) -> None:
    old_manifest = output / manifest_name
    if old_manifest.exists():
        old_manifest.unlink()  # never leave a stale manifest describing half-replaced files
    for name in names:
        os.replace(staging / name, output / name)
    os.replace(staging / manifest_name, output / manifest_name)


def _work_dir(output: Path | None) -> tuple[Path, bool]:
    if output is None:
        return Path(tempfile.mkdtemp(prefix="nio-dataset-")), True
    work = output / f".work-{os.getpid()}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    return work, True


def scan(roots: Sequence[Any], *, output: Path | None = None, force: bool = False,
         options: ExportOptions | None = None, argv: Sequence[str] | None = None) -> ExportResult:
    """Accounting only (no extxyz): every run and frame with its outcome, pools, lineage, duplicates.

    With ``output`` the accounting files and ``scan_manifest.json`` are written there (same refusal
    rules as ``export``); without it nothing is written and the result carries the summary.
    """
    options = options or ExportOptions()
    if output is not None:
        output = prepare_output_dir(Path(output), force=force)
    work, _ = _work_dir(output)
    try:
        analysis = analyze(roots, options, output_dir=output, spill_dir=work / "spill")
        manifest = _manifest_base(analysis, "scan")
        files: dict[str, Any] = {}
        if output is not None:
            staging = work / "out"
            staging.mkdir()
            _write_accounting(analysis, staging)
            names = [RUNS_FILE, FRAMES_FILE, EXCLUSIONS_FILE, AUDIT_JSON, AUDIT_MD]
            files = {name: _file_entry(staging / name) for name in names}
            manifest["files"] = files
            manifest["content_sha256"] = _content_sha256(manifest)
            manifest["invocation"] = _invocation(argv)
            atomic_write_json(staging / SCAN_MANIFEST, manifest)
            _publish(staging, output, names, SCAN_MANIFEST)
        else:
            manifest["content_sha256"] = _content_sha256(manifest)
        return ExportResult(output, manifest, analysis, files)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def export(roots: Sequence[Any], output: Path, *, force: bool = False, options: ExportOptions | None = None,
           argv: Sequence[str] | None = None) -> ExportResult:
    """Full export: accounting files + verified ``dataset.extxyz`` + ``dataset_manifest.json``."""
    options = options or ExportOptions()
    xyz._require_ase()  # fail before any work when ASE is absent
    output = prepare_output_dir(Path(output), force=force)
    work, _ = _work_dir(output)
    try:
        analysis = analyze(roots, options, output_dir=output, spill_dir=work / "spill")
        if analysis.pool_error is not None:
            raise DatasetError(analysis.pool_error)
        exported = analysis.exported()
        if not exported:
            counts = _counts(analysis)
            top = sorted(counts["frames_by_reason"].items(), key=lambda kv: (-kv[1], kv[0]))[:8]
            raise DatasetError(
                f"no frame is exportable ({counts['frames']} frames in {counts['runs']} runs; "
                f"{counts['runs_by_reason']}); frames by outcome: {dict(top)}. "
                "Run 'nio-md-prep dataset scan --output DIR' for the full account."
            )
        group_counts = check_group_consistency(exported)
        staging = work / "out"
        staging.mkdir()
        _write_accounting(analysis, staging)
        contract = xyz.FrameContract(allow_energy_only=options.policy.allow_energy_only)
        result = xyz.write_extxyz(staging / DATASET_FILE, _payloads(analysis, exported),
                                  keys=options.label_keys, verify=True, contract=contract)
        if result.frame_ids != [frame.frame_id for frame in exported]:
            raise DatasetError("internal error: dataset.extxyz frame order differs from the accounting")
        names = [DATASET_FILE, RUNS_FILE, FRAMES_FILE, EXCLUSIONS_FILE, AUDIT_JSON, AUDIT_MD]
        files = {name: _file_entry(staging / name) for name in names}
        files[DATASET_FILE].update({"frames": result.frames, "label_sets": dict(sorted(result.label_sets.items())),
                                    "stress_frames": result.stress_frames, "round_trip_verified": result.verified})
        manifest = _manifest_base(analysis, "export")
        manifest["groups"]["frames_by_group"] = group_counts
        manifest["files"] = files
        manifest["content_sha256"] = _content_sha256(manifest)
        manifest["invocation"] = _invocation(argv)
        atomic_write_json(staging / EXPORT_MANIFEST, manifest)
        _publish(staging, output, names, EXPORT_MANIFEST)
        return ExportResult(output, manifest, analysis, files)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def load_manifest(export_dir: Path, *, name: str = EXPORT_MANIFEST) -> dict[str, Any]:
    """Read a dataset/scan manifest and check its ``content_sha256``."""
    import json

    path = Path(export_dir) / name
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise DatasetError(f"{path} does not exist; run 'nio-md-prep dataset export' first") from None
    except json.JSONDecodeError as exc:
        raise DatasetError(f"{path} is not valid JSON: {exc}") from None
    if manifest.get("schema") not in (EXPORT_SCHEMA, SCAN_SCHEMA):
        raise DatasetError(f"{path.name} is not a nio-md-prep dataset manifest")
    if manifest.get("content_sha256") != _content_sha256(manifest):
        raise DatasetError(f"{path.name}: content_sha256 does not match its content (edited or corrupted)")
    return manifest


__all__ = [
    "Analysis", "ExportOptions", "ExportResult", "analyze", "check_group_consistency", "export",
    "load_manifest", "scan", "DATASET_FILE", "EXPORT_MANIFEST", "SCAN_MANIFEST", "RUNS_FILE", "FRAMES_FILE",
    "EXCLUSIONS_FILE", "AUDIT_JSON", "AUDIT_MD", "LIMITATIONS",
]
