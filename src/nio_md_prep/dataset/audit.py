"""Re-verification and summaries of an exported dataset (``nio-md-prep dataset audit``).

Two layers:

* :func:`summarize` / :func:`render_markdown` -- deterministic counts of every run and frame by
  outcome state and reason, campaign, family, temperature, composition, lineage group, VASP
  version, force and stress availability, magnetic class (and a few more dimensions), plus the
  scalar audit figures (stress coverage, max |F - E0| per atom, lineage sources, pools, magnetic
  flags). The exporter writes these as ``audit_report.json`` / ``audit.md`` next to the data.
* :func:`audit_export` -- independent re-verification of an export directory: the manifest's
  ``content_sha256`` and every file hash; the accounting (every frame and run has exactly one
  valid outcome, counts agree with the manifest, the exported frames are exactly the accepted
  ones, in order); an ASE round trip of ``dataset.extxyz`` compared with the ``frames.jsonl`` label
  digests, lineage groups and structure keys (no calculator, no constraints, pbc T T T); structure
  hashes never shared across lineage groups; and, for each split directory, the split manifest,
  its file hashes and zero leakage between train/valid/test by lineage group, order-dependent
  structure key and permutation-invariant key. Optional ``verify_sources`` re-hashes the source
  vasprun files under relocated roots. Results go to ``audit_report.json`` + ``audit.md`` in the
  audit output directory (never into the hashed export files).

Nothing here is a scientific validation of the labels: the audit proves that the dataset is the
faithful, leak-free, fully accounted transcription of the VASP files it was built from.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .errors import DatasetError
from .fsio import atomic_write_json, atomic_write_text, prepare_output_dir, read_jsonl, sha256_file
from .model import ACCEPTED, OUTCOMES, REASON_STATUS, REASONS

AUDIT_SCHEMA = "nio-md-prep.dataset-audit"
AUDIT_SCHEMA_VERSION = 1
SUMMARY_SCHEMA = "nio-md-prep.dataset-summary"
MISSING = "(missing)"
#: Frame dimensions counted in every summary (record field -> label).
FRAME_DIMENSIONS = (
    "status", "reason", "campaign", "family", "temperature_K", "composition", "lineage_group",
    "vasp_version", "forces_available", "stress_available", "magnetic_class", "magnetic_policy", "label_source",
    "label_set", "calc_type", "source_file_type", "pool_id", "lineage_source", "energy_rule", "scf_status",
)
RUN_DIMENSIONS = (
    "status", "reason", "campaign", "family", "temperature_K", "composition", "vasp_version", "magnetic_class",
    "calc_type", "lineage_source", "pool_id", "source_file_type",
)
_MARKDOWN_ROWS = 40


def _label(value: Any) -> str:
    if value is None:
        return MISSING
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and value.is_integer():
        return str(int(value)) if abs(value) < 1e15 else repr(value)
    return str(value)


def _frame_value(record: Mapping[str, Any], dimension: str) -> Any:
    outcome = record.get("outcome") or {}
    if dimension == "status":
        return outcome.get("status")
    if dimension == "reason":
        return f"{outcome.get('status')}/{outcome.get('reason') or 'accepted'}"
    if dimension == "energy_rule":
        return record.get("label_energy_rule")
    if dimension == "scf_status":
        return (record.get("scf") or {}).get("status")
    return record.get(dimension)


def _run_value(record: Mapping[str, Any], dimension: str) -> Any:
    outcome = record.get("outcome") or {}
    metadata = record.get("metadata") or {}
    if dimension == "status":
        return outcome.get("status")
    if dimension == "reason":
        return f"{outcome.get('status')}/{outcome.get('reason') or 'accepted'}"
    if dimension in ("campaign", "family", "temperature_K"):
        return metadata.get(dimension)
    if dimension == "magnetic_class":
        magnetic = ((record.get("assessment") or {}).get("magnetic") or {}).get("magnetic_class")
        if magnetic is None and (record.get("magnetic_evidence") or {}).get("magnetic_class"):
            return f"{record['magnetic_evidence']['magnetic_class']} (evidence only, no labels)"
        return magnetic
    if dimension == "lineage_source":
        return (record.get("lineage") or {}).get("source")
    return record.get(dimension)


def _table(records: Iterable[Mapping[str, Any]], dimensions: Sequence[str], value_of, exported_of) -> dict[str, Any]:
    tables: dict[str, dict[str, dict[str, Any]]] = {name: {} for name in dimensions}
    for record in records:
        status = (record.get("outcome") or {}).get("status") or MISSING
        exported = exported_of(record)
        for name in dimensions:
            key = _label(value_of(record, name))
            row = tables[name].setdefault(key, {"total": 0, "exported": 0, "by_status": {}})
            row["total"] += 1
            row["exported"] += int(bool(exported))
            row["by_status"][status] = row["by_status"].get(status, 0) + 1
    return {name: {key: {**row, "by_status": dict(sorted(row["by_status"].items()))}
                   for key, row in sorted(rows.items())} for name, rows in tables.items()}


def summarize(runs: Sequence[Mapping[str, Any]], frames: Sequence[Mapping[str, Any]], *,
              pools: Sequence[Mapping[str, Any]] = (), lineage: Mapping[str, Any] | None = None,
              ignored: Sequence[Mapping[str, Any]] = ()) -> dict[str, Any]:
    """Counts by every audit dimension (deterministic; no timestamps, no absolute paths)."""
    exported = [f for f in frames if (f.get("outcome") or {}).get("status") == ACCEPTED]
    groups: dict[str, int] = {}
    for frame in exported:
        group = _label(frame.get("lineage_group"))
        groups[group] = groups.get(group, 0) + 1
    worst = None
    for frame in exported:
        energy = frame.get("energy") or {}
        f_value, e0 = energy.get("free_energy"), energy.get("energy_sigma0")
        n_atoms = frame.get("n_atoms")
        if f_value is not None and e0 is not None and n_atoms:
            value = abs(f_value - e0) / n_atoms
            if worst is None or value > worst["value"]:
                worst = {"value": value, "frame_id": frame["frame_id"]}
    with_stress = sum(1 for f in exported if f.get("stress_available"))
    energy_only = sum(1 for f in exported if f.get("label_set") == "energy_only")
    run_flags: dict[str, int] = {}
    for run in runs:
        for flag in (run.get("assessment") or {}).get("flags") or []:
            run_flags[flag] = run_flags.get(flag, 0) + 1
    magnetic_flags = {k: v for k, v in sorted(run_flags.items()) if "magnet" in k or k.startswith("mag")}
    transitions = sum(len((((run.get("assessment") or {}).get("magnetic")) or {}).get("transitions") or [])
                      for run in runs)
    ignored_by_reason: dict[str, int] = {}
    for item in ignored:
        key = f"{item.get('kind')}/{item.get('reason')}"
        ignored_by_reason[key] = ignored_by_reason.get(key, 0) + 1
    largest = sorted(groups.items(), key=lambda kv: (-kv[1], kv[0]))[:10]
    return {
        "schema": SUMMARY_SCHEMA,
        "totals": {
            "runs": len(runs), "frames": len(frames), "exported_frames": len(exported),
            "exported_groups": len(groups),
            "exported_runs": sum(1 for run in runs if (run.get("frames_exported") or 0) > 0),
            "parsed_runs": sum(1 for run in runs if run.get("parsed")),
        },
        "frames": _table(frames, FRAME_DIMENSIONS, _frame_value,
                         lambda r: (r.get("outcome") or {}).get("status") == ACCEPTED),
        "runs": _table(runs, RUN_DIMENSIONS, _run_value, lambda r: (r.get("frames_exported") or 0) > 0),
        "exported": {
            "stress_frames": with_stress,
            "stress_coverage": (with_stress / len(exported)) if exported else None,
            "energy_only_frames": energy_only,
            "max_abs_F_minus_E0_per_atom_eV": worst,
            "largest_groups": [{"lineage_group": k, "frames": v} for k, v in largest],
        },
        "magnetic": {"run_flags": magnetic_flags, "transitions": transitions},
        "run_flags": dict(sorted(run_flags.items())),
        "lineage": dict(lineage or {}),
        "pools": [{k: pool.get(k) for k in ("pool_id", "runs", "accepted_frames", "selected", "unknown_fields")}
                  for pool in pools],
        "discovery_ignored": dict(sorted(ignored_by_reason.items())),
        "notice": "Counts describe the VASP files found under the scanned roots. They do not validate "
                  "forces or energies against an independent reference calculation.",
    }


def _md_table(title: str, table: Mapping[str, Mapping[str, Any]], *, limit: int = _MARKDOWN_ROWS) -> list[str]:
    lines = [f"### {title}", "", "| value | total | exported | by state |", "|---|---:|---:|---|"]
    rows = sorted(table.items(), key=lambda kv: (-kv[1]["total"], kv[0]))
    for key, row in rows[:limit]:
        states = ", ".join(f"{s}: {n}" for s, n in row["by_status"].items())
        lines.append(f"| `{key}` | {row['total']} | {row['exported']} | {states} |")
    if len(rows) > limit:
        lines.append(f"| ... {len(rows) - limit} more | | | |")
    lines.append("")
    return lines


def render_markdown(summary: Mapping[str, Any], *, title: str = "Dataset audit",
                    checks: Sequence[Mapping[str, Any]] = (), extra: Sequence[str] = ()) -> str:
    totals = summary["totals"]
    exported = summary["exported"]
    lines = [f"# {title}", ""]
    lines += [f"- runs: {totals['runs']} ({totals['parsed_runs']} parsed, {totals['exported_runs']} with exported frames)",
              f"- frames: {totals['frames']} discovered, {totals['exported_frames']} exported "
              f"in {totals['exported_groups']} lineage groups",
              f"- stress on exported frames: {exported['stress_frames']}"
              + (f" ({exported['stress_coverage']:.1%})" if exported["stress_coverage"] is not None else ""),
              f"- energy-only exported frames: {exported['energy_only_frames']}"]
    worst = exported.get("max_abs_F_minus_E0_per_atom_eV")
    if worst:
        lines.append(f"- max |F - E0| per atom: {worst['value']:.3e} eV (`{worst['frame_id']}`)")
    lines += ["", f"> {summary['notice']}", ""]
    if checks:
        lines += ["## Verification", "", "| check | result | detail |", "|---|---|---|"]
        for check in checks:
            detail = str(check.get("detail") or "").replace("|", "/").replace("\n", " ")
            lines.append(f"| {check['name']} | {'ok' if check['ok'] else 'FAILED'} | {detail[:300]} |")
        lines.append("")
    lines += list(extra)
    lines += ["## Frames", ""]
    for name, table in summary["frames"].items():
        lines += _md_table(f"frames by {name}", table)
    lines += ["## Runs", ""]
    for name, table in summary["runs"].items():
        lines += _md_table(f"runs by {name}", table)
    if summary.get("pools"):
        lines += ["## Reference-settings pools", "", "| pool | runs | accepted frames | selected | unknown fields |",
                  "|---|---:|---:|---|---|"]
        for pool in summary["pools"]:
            lines.append(f"| `{pool['pool_id']}` | {pool['runs']} | {pool['accepted_frames']} | "
                         f"{pool['selected']} | {', '.join(pool.get('unknown_fields') or []) or '-'} |")
        lines.append("")
    if summary.get("magnetic", {}).get("run_flags"):
        lines += ["## Magnetic flags (runs)", ""]
        lines += [f"- `{k}`: {v}" for k, v in summary["magnetic"]["run_flags"].items()]
        lines.append("")
    if summary.get("discovery_ignored"):
        lines += ["## Discovery: ignored paths", ""]
        lines += [f"- `{k}`: {v}" for k, v in summary["discovery_ignored"].items()]
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


# --------------------------------------------------------------------------
# Re-verification of an export directory
# --------------------------------------------------------------------------

class _Checks:
    def __init__(self):
        self.items: list[dict[str, Any]] = []

    def add(self, name: str, ok: bool, detail: Any = "") -> bool:
        self.items.append({"name": name, "ok": bool(ok), "detail": detail if isinstance(detail, str) else
                           json.dumps(detail, sort_keys=True, default=str)[:2000]})
        return bool(ok)

    @property
    def ok(self) -> bool:
        return all(item["ok"] for item in self.items)


def _check_outcome(outcome: Any) -> str | None:
    if not isinstance(outcome, Mapping):
        return "missing outcome"
    status, reason = outcome.get("status"), outcome.get("reason")
    if status not in OUTCOMES:
        return f"unknown status {status!r}"
    if status == ACCEPTED:
        return None if reason is None else f"accepted with reason {reason!r}"
    if reason not in REASONS:
        return f"unknown reason {reason!r}"
    if REASON_STATUS[reason] != status:
        return f"reason {reason!r} reported as {status!r} (expected {REASON_STATUS[reason]!r})"
    return None


def _verify_files(export_dir: Path, manifest: Mapping[str, Any], checks: _Checks) -> None:
    files = manifest.get("files") or {}
    bad = []
    for name, entry in sorted(files.items()):
        path = export_dir / name
        if not path.is_file():
            bad.append(f"{name}: missing")
            continue
        if path.stat().st_size != entry.get("bytes"):
            bad.append(f"{name}: {path.stat().st_size} bytes, manifest {entry.get('bytes')}")
        if sha256_file(path) != entry.get("sha256"):
            bad.append(f"{name}: sha256 differs from the manifest")
    checks.add("file_hashes", not bad and bool(files), "; ".join(bad) or f"{len(files)} files verified")


def _verify_records(runs, frames, manifest, checks: _Checks) -> dict[str, Any]:
    problems: list[str] = []
    run_ids = [r.get("run_id") for r in runs]
    if len(set(run_ids)) != len(run_ids):
        problems.append("runs.jsonl repeats a run_id")
    known = set(run_ids)
    for record in runs:
        why = _check_outcome(record.get("outcome"))
        if why:
            problems.append(f"run {record.get('run_id')}: {why}")
    frame_ids = [f.get("frame_id") for f in frames]
    if len(set(frame_ids)) != len(frame_ids):
        problems.append("frames.jsonl repeats a frame_id")
    for record in frames:
        why = _check_outcome(record.get("outcome"))
        if why:
            problems.append(f"frame {record.get('frame_id')}: {why}")
        if record.get("run_id") not in known:
            problems.append(f"frame {record.get('frame_id')}: run {record.get('run_id')!r} not in runs.jsonl")
    counts = manifest.get("counts") or {}
    by_reason: dict[str, int] = {}
    for record in frames:
        outcome = record.get("outcome") or {}
        key = f"{outcome.get('status')}/{outcome.get('reason') or 'accepted'}"
        by_reason[key] = by_reason.get(key, 0) + 1
    if counts and dict(sorted(by_reason.items())) != counts.get("frames_by_reason"):
        problems.append("frame counts by reason differ from the manifest")
    if counts and counts.get("runs") != len(runs):
        problems.append(f"manifest counts {counts.get('runs')} runs, runs.jsonl has {len(runs)}")
    checks.add("accounting", not problems, "; ".join(problems[:20]) or
               f"{len(runs)} runs and {len(frames)} frames, each with exactly one valid outcome")
    return {"problems": problems}


def _round_trip(dataset: Path, exported: Sequence[Mapping[str, Any]], keys: Mapping[str, str],
                checks: _Checks) -> dict[str, Any]:
    from . import extxyz as xyz

    xyz._require_ase()
    import ase.io
    import numpy as np

    label_keys = xyz.LabelKeys(keys.get("energy", "REF_energy"), keys.get("forces", "REF_forces"),
                               keys.get("stress", "REF_stress"))
    expected = {f["frame_id"]: f for f in exported}
    order = [f["frame_id"] for f in exported]
    problems: list[str] = []
    read_ids: list[str] = []
    for atoms in ase.io.iread(str(dataset), index=":", format="extxyz"):
        view = xyz._view_from_atoms(atoms, label_keys)
        frame_id = view.frame_id
        read_ids.append(frame_id)
        record = expected.get(frame_id)
        if record is None:
            problems.append(f"{frame_id}: in dataset.extxyz but not an accepted frame of frames.jsonl")
            continue
        stress = None if view.stress is None else np.asarray(view.stress, dtype=np.float64).reshape(3, 3)
        if not isinstance(view.energy, float):
            problems.append(f"{frame_id}: energy missing")
            continue
        digest = xyz.label_sha256(view.energy, view.forces, stress)
        if digest != record.get("label_sha256"):
            problems.append(f"{frame_id}: label digest differs from frames.jsonl")
        info = view.info
        for key in ("lineage_group", "structure_key"):
            if info.get(key) != record.get(key):
                problems.append(f"{frame_id}: {key} {info.get(key)!r} != frames.jsonl {record.get(key)!r}")
        if bool(view.forces is not None) != bool(record.get("forces_available")):
            problems.append(f"{frame_id}: forces presence differs from frames.jsonl")
        if bool(stress is not None) != bool(record.get("stress_available")):
            problems.append(f"{frame_id}: stress presence differs from frames.jsonl")
        if tuple(view.pbc) != (True, True, True):
            problems.append(f"{frame_id}: pbc {view.pbc}")
        if view.has_calc or view.constraints:
            problems.append(f"{frame_id}: ASE attached a calculator or constraints")
        if len(view.species) != record.get("n_atoms"):
            problems.append(f"{frame_id}: {len(view.species)} atoms, frames.jsonl {record.get('n_atoms')}")
    if read_ids != order:
        problems.append(f"dataset.extxyz holds {len(read_ids)} frames; frames.jsonl accepts {len(order)} "
                        "(or the order differs)")
    checks.add("ase_round_trip", not problems, "; ".join(problems[:20]) or
               f"{len(read_ids)} frames re-read with ASE; label digests, lineage groups and structure keys match")
    return {"frames_read": len(read_ids), "problems": problems[:100]}


def _structure_leakage(frames: Iterable[Mapping[str, Any]], unit_of) -> list[str]:
    owners: dict[tuple[str, str], str] = {}
    problems = []
    for frame in frames:
        unit = unit_of(frame)
        for kind in ("lineage_group", "structure_key", "permutation_key"):
            key = frame.get(kind)
            if key is None:
                problems.append(f"{frame.get('frame_id')}: no {kind}")
                continue
            other = owners.setdefault((kind, str(key)), unit)
            if other != unit:
                problems.append(f"{kind} {str(key)[:40]} in both {other!r} and {unit!r}")
    return problems


def _verify_split(split_dir: Path, export_dir: Path, by_id: Mapping[str, Mapping[str, Any]],
                  exported_ids: set[str], checks: _Checks) -> dict[str, Any]:
    from .extxyz import SPLIT_KEY, iter_extxyz_blocks
    from .split import SPLITS, load_split_manifest

    name = f"split:{split_dir.name}"
    try:
        manifest = load_split_manifest(split_dir / "split_manifest.json")
    except DatasetError as exc:
        checks.add(name, False, str(exc))
        return {"dir": split_dir.as_posix(), "ok": False}
    problems: list[str] = []
    export_manifest = export_dir / "dataset_manifest.json"
    if manifest.get("export_manifest_sha256") != sha256_file(export_manifest):
        problems.append("split was made from a different export manifest")
    assigned: dict[str, str] = {}
    for split in SPLITS:
        entry = (manifest.get("files") or {}).get(split) or {}
        path = split_dir / f"{split}.extxyz"
        if not path.is_file() or sha256_file(path) != entry.get("sha256"):
            problems.append(f"{split}.extxyz missing or its sha256 differs from split_manifest.json")
            continue
        ids = []
        for block in iter_extxyz_blocks(path):
            ids.append(block.frame_id)
            if block.info().get(SPLIT_KEY) != split:
                problems.append(f"{split}.extxyz frame {block.frame_id} has split={block.info().get(SPLIT_KEY)!r}")
                break
        listed = list(manifest.get("frame_ids", {}).get(split, []))
        if sorted(ids) != sorted(listed):
            problems.append(f"{split}.extxyz frames differ from split_manifest.json")
        for frame_id in ids:
            if frame_id in assigned:
                problems.append(f"{frame_id} in {assigned[frame_id]} and {split}")
            assigned[frame_id] = split
    missing = exported_ids - set(assigned)
    extra = set(assigned) - exported_ids
    if missing or extra:
        problems.append(f"{len(missing)} exported frames in no split, {len(extra)} split frames not exported")
    leakage = _structure_leakage((by_id[f] for f in sorted(assigned) if f in by_id), lambda f: assigned[f["frame_id"]])
    # a lineage group / structure key shared by two splits is leakage
    problems += [f"leakage: {p}" for p in leakage]
    checks.add(name, not problems, "; ".join(problems[:20]) or
               f"{len(assigned)} frames; zero lineage-group/structure-key/permutation-key leakage between splits")
    return {"dir": split_dir.as_posix(), "ok": not problems, "frames": len(assigned),
            "leakage": leakage[:50], "problems": problems[:50]}


def _verify_sources(runs: Sequence[Mapping[str, Any]], roots: Mapping[str, Path], checks: _Checks) -> dict[str, Any]:
    checked = 0
    problems = []
    for run in runs:
        if not (run.get("frames_exported") or 0):
            continue
        alias = run.get("root_alias")
        if alias not in roots:
            continue
        relpath = run.get("relpath")
        directory = Path(roots[alias]) if relpath == "." else Path(roots[alias]) / relpath
        path = directory / str(run.get("source_file"))
        checked += 1
        if not path.is_file():
            problems.append(f"{run['run_id']}: {path.as_posix()} missing")
        elif sha256_file(path) != run.get("source_sha256"):
            problems.append(f"{run['run_id']}: source sha256 differs")
    checks.add("source_hashes", not problems and checked > 0,
               "; ".join(problems[:20]) or f"{checked} source files re-hashed")
    return {"checked": checked, "problems": problems[:100]}


def _parse_source_roots(values: Mapping[str, Any] | Sequence[str] | None) -> dict[str, Path]:
    if values is None:
        return {}
    if isinstance(values, Mapping):
        return {str(k): Path(v) for k, v in values.items()}
    roots = {}
    for value in values:
        if "=" not in str(value):
            raise DatasetError(f"--verify-sources expects alias=path, got {value!r}")
        alias, path = str(value).split("=", 1)
        roots[alias] = Path(path)
    return roots


def audit_export(export_dir: Path, *, split_dirs: Sequence[Path] = (), verify_sources: Any = None,
                 output: Path | None = None, force: bool = False) -> dict[str, Any]:
    """Re-verify an export (and splits); write ``audit_report.json`` + ``audit.md``; return the report.

    ``output`` defaults to ``<export_dir>/audit`` (owned by the audit and rewritten on every run);
    an explicit ``output`` must be new or empty unless ``force``. ``report["ok"]`` is False when
    any check fails (the CLI then exits with status 1).
    """
    from .export import AUDIT_JSON, AUDIT_MD, EXPORT_MANIFEST, FRAMES_FILE, RUNS_FILE, _content_sha256

    export_dir = Path(export_dir)
    manifest_path = export_dir / EXPORT_MANIFEST
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise DatasetError(f"{manifest_path} does not exist; audit an export directory") from None
    except json.JSONDecodeError as exc:
        raise DatasetError(f"{manifest_path} is not valid JSON: {exc}") from None
    checks = _Checks()
    checks.add("manifest_content_sha256", manifest.get("content_sha256") == _content_sha256(manifest),
               "content_sha256 matches" if manifest.get("content_sha256") == _content_sha256(manifest)
               else "content_sha256 does not match the manifest content (edited or corrupted)")
    _verify_files(export_dir, manifest, checks)
    runs = read_jsonl(export_dir / RUNS_FILE)
    frames = read_jsonl(export_dir / FRAMES_FILE)
    _verify_records(runs, frames, manifest, checks)
    exported = [f for f in frames if (f.get("outcome") or {}).get("status") == ACCEPTED]
    round_trip = _round_trip(export_dir / "dataset.extxyz", exported, manifest.get("label_keys") or {}, checks)
    leakage = _structure_leakage(exported, lambda f: f.get("lineage_group"))
    checks.add("group_structure_consistency", not leakage, "; ".join(leakage[:20]) or
               "every structure hash of the export belongs to exactly one lineage group")
    by_id = {f["frame_id"]: f for f in exported}
    splits = [_verify_split(Path(d), export_dir, by_id, set(by_id), checks) for d in split_dirs]
    sources = None
    roots = _parse_source_roots(verify_sources)
    if roots:
        sources = _verify_sources(runs, roots, checks)
    summary = summarize(runs, frames, pools=manifest.get("pools") or [],
                        lineage=(manifest.get("lineage") or {}).get("sources") or {})
    report = {
        "schema": AUDIT_SCHEMA, "schema_version": AUDIT_SCHEMA_VERSION,
        "export_manifest_sha256": sha256_file(manifest_path),
        "export_content_sha256": manifest.get("content_sha256"),
        "ok": checks.ok, "checks": checks.items, "round_trip": round_trip, "splits": splits,
        "sources": sources, "summary": summary, "limitations": list(manifest.get("limitations") or []),
    }
    if output is None:
        out = export_dir / "audit"
        out.mkdir(exist_ok=True)
    else:
        out = prepare_output_dir(Path(output), force=force)
    atomic_write_json(out / AUDIT_JSON, report)
    extra = ["## Limitations", ""] + [f"- {item}" for item in report["limitations"]] + [""]
    atomic_write_text(out / AUDIT_MD, render_markdown(summary, title="Dataset audit", checks=checks.items,
                                                      extra=extra))
    report["output"] = out.as_posix()
    return report


__all__ = ["audit_export", "render_markdown", "summarize", "FRAME_DIMENSIONS", "RUN_DIMENSIONS"]
