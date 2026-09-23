"""Seeded, lineage-group-level train/valid/test splits with a frozen test set on append.

The unit of assignment is the *lineage group* (``lineage_group``; the
exporter's union-find component of lineage, continuation, nested-run and
duplicate links, also accepted under its older name ``group_id``). Frames are
never split individually: all frames of a group always land in one split, so
correlated frames (one AIMD trajectory, one Packmol replica, a heating/hold
continuation, a restart chain, an InterfaceForge OPT/Step1/Step2 case,
exact duplicates) can never leak between train, valid and test.

Before assigning anything the splitter re-checks the lineage it is given and
raises :class:`LeakageError` when correlated frames carry different groups:
frames of one run (``run_id``) in two groups, or an identical structure
(``structure_key`` / ``permutation_key``) in two groups.

Algorithm (deterministic, independent of input order):

1. groups are sorted by id; a numpy ``PCG64(seed)`` permutation of that list
   gives each group a tie-break rank;
2. groups are taken largest first (frames), ties by rank, and each goes to the
   active split with the largest remaining *frame* deficit
   ``target - assigned`` (ties: larger requested fraction, then
   train < valid < test);
3. if a requested (non-zero) split is still empty and enough groups exist,
   the smallest group of the split with the largest surplus is moved there
   (whole groups only).

With ``stratify_by`` the same procedure runs independently inside each
stratum (a metadata value that must be constant within every group). With a
``previous`` split manifest, groups keep their previous split, the test split
is frozen (new groups go to train/valid with renormalised fractions unless
``grow_test``), a frozen test group that gains frames is an error, and a group
that now joins frames previously assigned to different splits raises
:class:`LeakageError`.

The manifest records per-split statistics (frames and groups by campaign,
family, temperature, composition, elements, VASP version, magnetic class,
stress availability, label set, ...) and coverage warnings for values that
could have been represented in every split but are not.
"""

# Step 2-3 adapted from lgutsev/InterfaceForge@4501e34 src/interfaceforge/leaf_collect.py
# assign_heritage_groups (MIT): largest-group-first greedy toward per-split targets with a
# deterministic seed, and the "move a whole lineage into an empty split" pass. Changed here:
# numpy PCG64 tie-break ranks over id-sorted groups (instead of random.Random shuffle),
# frame-count targets, strata, frozen previous assignments and leakage checks.

from __future__ import annotations

import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

from .errors import DatasetError, LeakageError
from .extxyz import SPLIT_KEY, iter_extxyz_blocks, read_frame_ids, slice_extxyz_by_frame_ids
from .fsio import atomic_write_json, atomic_write_text, canonical_json, prepare_output_dir, sha256_bytes, sha256_file

SPLITS = ("train", "valid", "test")
SPLIT_SCHEMA = "nio-md-prep.dataset-split"
SPLIT_SCHEMA_VERSION = 1
ALGORITHM = "group-greedy-largest-frame-deficit/v1"
RNG = "numpy.random.PCG64"
DEFAULT_SEED = 11
DEFAULT_FRACTIONS = (0.8, 0.1, 0.1)
#: Frame field holding the split unit; ``group_id`` is accepted as an alias.
GROUP_KEY = "lineage_group"
#: Achieved frame fractions further than this from the request produce a warning.
FRACTION_WARNING_TOLERANCE = 0.05
#: Keys added by :func:`write_split` that are not covered by ``content_sha256``.
NON_CONTENT_KEYS = ("content_sha256", "files", "summary_file", "dataset_extxyz_sha256", "export_manifest_file_sha256")
#: Metadata keys counted per split in ``statistics`` (missing -> "(missing)").
DEFAULT_STATS_KEYS = (
    "campaign", "family", "config_type", "temperature_K", "composition", "elements", "calc_type",
    "vasp_version", "magnetic_class", "stress_available", "label_set", "source_file_type",
)
#: Keys whose values should appear in every non-empty split when they have enough groups.
DEFAULT_COVERAGE_KEYS = ("campaign", "family", "temperature_K", "elements")
MISSING = "(missing)"
#: Structure-hash fields checked for cross-group leakage when present.
STRUCTURE_KEYS = ("structure_key", "permutation_key")


def _get(record: Any, key: str, default: Any = None) -> Any:
    if isinstance(record, Mapping):
        return record.get(key, default)
    return getattr(record, key, default)


def _outcome_status(record: Any) -> str | None:
    outcome = _get(record, "outcome")
    if outcome is None:
        return None
    if isinstance(outcome, Mapping):
        return outcome.get("status")
    return getattr(outcome, "status", None)


def _lineage_source(record: Any) -> str | None:
    lineage = _get(record, "lineage")
    if isinstance(lineage, Mapping) and lineage.get("source"):
        return lineage.get("source")
    return _get(record, "lineage_source")


def _field(record: Any, key: str) -> Any:
    """``metadata[key]`` first, then the top-level field."""
    metadata = _get(record, "metadata")
    if isinstance(metadata, Mapping) and key in metadata:
        return metadata[key]
    return _get(record, key)


def _group_of(record: Any, frame_id: str) -> Any:
    lineage_group = _field(record, GROUP_KEY)
    group_id = _field(record, "group_id")
    if lineage_group is not None and group_id is not None and lineage_group != group_id:
        raise DatasetError(
            f"frame {frame_id!r}: lineage_group {lineage_group!r} and group_id {group_id!r} disagree"
        )
    return lineage_group if lineage_group is not None else group_id


def _label(value: Any) -> str:
    if value is None:
        return MISSING
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    if isinstance(value, (list, tuple, dict)):
        return canonical_json(value)
    return str(value)


def _check_fractions(fractions: Sequence[float]) -> dict[str, float]:
    values = list(fractions)
    if len(values) != len(SPLITS):
        raise DatasetError(f"fractions must give {len(SPLITS)} values (train, valid, test), got {values}")
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise DatasetError(f"fractions must be finite and non-negative, got {values}")
    if abs(sum(values) - 1.0) > 1e-6:
        raise DatasetError(f"fractions must sum to 1, got {values} (sum {sum(values)})")
    if values[0] <= 0:
        raise DatasetError("the train fraction must be positive")
    return {name: float(value) for name, value in zip(SPLITS, values)}


def _check_seed(seed: Any) -> int:
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise DatasetError(f"seed must be a non-negative integer, got {seed!r}")
    return seed


def _ranks(group_ids: Sequence[str], seed: int) -> dict[str, int]:
    import numpy as np

    ordered = sorted(group_ids)
    permutation = np.random.Generator(np.random.PCG64(seed)).permutation(len(ordered))
    return {group_id: int(rank) for group_id, rank in zip(ordered, permutation)}


def _assign_groups(
    sizes: Mapping[str, int], fractions: Mapping[str, float], ranks: Mapping[str, int]
) -> dict[str, str]:
    """Greedy largest-first assignment of whole groups (see module docstring)."""
    active = [name for name in SPLITS if fractions[name] > 0]
    if not active:
        raise DatasetError("no split has a positive fraction")
    total = sum(sizes.values())
    target = {name: fractions[name] * total for name in SPLITS}
    assigned = {name: 0 for name in SPLITS}
    counts = {name: 0 for name in SPLITS}
    result: dict[str, str] = {}
    order = sorted(sizes, key=lambda group_id: (-sizes[group_id], ranks[group_id], group_id))
    for group_id in order:
        split = max(active, key=lambda name: (target[name] - assigned[name], fractions[name], -SPLITS.index(name)))
        result[group_id] = split
        assigned[split] += sizes[group_id]
        counts[split] += 1
    if len(order) >= len(active):
        for empty in [name for name in active if counts[name] == 0]:
            donors = [name for name in active if counts[name] > 1]
            if not donors:
                continue
            donor = max(donors, key=lambda name: (assigned[name] - target[name], -SPLITS.index(name)))
            members = [group_id for group_id in order if result[group_id] == donor]
            moved = min(members, key=lambda group_id: (sizes[group_id], ranks[group_id], group_id))
            result[moved] = empty
            assigned[donor] -= sizes[moved]
            assigned[empty] += sizes[moved]
            counts[donor] -= 1
            counts[empty] += 1
    return result


def _previous_index(previous: Mapping[str, Any]) -> tuple[dict[str, str], dict[str, str], str]:
    if previous.get("schema") != SPLIT_SCHEMA:
        raise DatasetError(f"previous split manifest has schema {previous.get('schema')!r}, expected {SPLIT_SCHEMA!r}")
    if previous.get("schema_version") != SPLIT_SCHEMA_VERSION:
        raise DatasetError(f"unsupported previous split schema_version {previous.get('schema_version')!r}")
    check_split_disjoint(previous)
    frame_split: dict[str, str] = {}
    for split in SPLITS:
        for frame_id in previous.get("frame_ids", {}).get(split, []):
            frame_split[frame_id] = split
    group_split = dict(previous.get("group_to_split", {}))
    for split in group_split.values():
        if split not in SPLITS:
            raise DatasetError(f"previous split manifest names unknown split {split!r}")
    digest = sha256_bytes(canonical_json(_content(previous)).encode())
    recorded = previous.get("content_sha256")
    if recorded is not None and recorded != digest:
        raise DatasetError("previous split manifest: content_sha256 does not match its content (edited or corrupted)")
    return frame_split, group_split, digest


def _content(manifest: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in manifest.items() if key not in NON_CONTENT_KEYS}


def check_split_disjoint(manifest: Mapping[str, Any]) -> None:
    """Raise :class:`LeakageError` unless every frame and group is in exactly one split."""
    owner: dict[str, str] = {}
    for split in SPLITS:
        for frame_id in manifest.get("frame_ids", {}).get(split, []):
            if frame_id in owner:
                raise LeakageError(f"frame {frame_id!r} is in both {owner[frame_id]} and {split}")
            owner[frame_id] = split
    group_split = manifest.get("group_to_split", {})
    for group_id, info in manifest.get("groups", {}).items():
        if group_split.get(group_id, info.get("split")) != info.get("split"):
            raise LeakageError(f"group {group_id!r} has inconsistent split records")


def _check_lineage_consistency(records: Mapping[str, Mapping[str, Any]]) -> dict[str, int]:
    """Correlated frames must share a group: one run -> one group, one structure -> one group."""
    checks = {"runs": 0, **{key: 0 for key in STRUCTURE_KEYS}}
    by_run: dict[str, set[str]] = {}
    by_structure: dict[tuple[str, str], set[str]] = {}
    for frame_id in sorted(records):
        record = records[frame_id]
        by_run.setdefault(record["run_id"], set()).add(record["group_id"])
        for key in STRUCTURE_KEYS:
            value = record[key]
            if value is not None:
                by_structure.setdefault((key, str(value)), set()).add(record["group_id"])
    checks["runs"] = len(by_run)
    split_runs = {run: sorted(groups) for run, groups in sorted(by_run.items()) if len(groups) > 1}
    if split_runs:
        run, groups = next(iter(split_runs.items()))
        raise LeakageError(
            f"{len(split_runs)} run(s) have frames in different lineage groups (e.g. run {run!r} in {groups}); "
            "all frames of one trajectory/run must share one lineage_group"
        )
    for (key, value), groups in sorted(by_structure.items()):
        checks[key] += 1
        if len(groups) > 1:
            raise LeakageError(
                f"identical structure ({key}={value[:16]}...) occurs in lineage groups {sorted(groups)}; exact "
                "duplicates must be unioned into one group (or deduplicated) before splitting"
            )
    return checks


def make_split(
    frames: Iterable[Any],
    *,
    seed: int = DEFAULT_SEED,
    fractions: Sequence[float] = DEFAULT_FRACTIONS,
    stratify_by: str | None = None,
    previous: Mapping[str, Any] | None = None,
    grow_test: bool = False,
    allow_underpopulated_strata: bool = False,
    export_manifest_sha256: str | None = None,
    stats_keys: Sequence[str] = DEFAULT_STATS_KEYS,
    coverage_keys: Sequence[str] = DEFAULT_COVERAGE_KEYS,
) -> dict[str, Any]:
    """Assign accepted frames to train/valid/test by lineage group; return the split manifest.

    ``frames``: accepted frame records (mappings or objects) with
    ``frame_id``, ``lineage_group`` (or ``group_id``), optional ``n_atoms``,
    ``run_id`` (default: the ``frame_id`` part before ``#``),
    ``structure_key``/``permutation_key``, ``metadata`` (dict; metadata keys
    are looked up there first, then as top-level fields), ``outcome`` and
    ``lineage``/``lineage_source``. Fails closed on: a non-accepted outcome, a
    missing or unresolved lineage group, duplicate frame ids, correlated
    frames in different groups (:class:`LeakageError`), fewer than 3 groups, a
    stratum value that varies inside a group or is missing, an
    under-populated stratum (unless allowed), a frozen test group that gained
    frames (unless ``grow_test``), and frames of one group previously in
    different splits (:class:`LeakageError`).
    """
    seed = _check_seed(seed)
    requested = _check_fractions(fractions)
    records: dict[str, dict[str, Any]] = {}
    unresolved: set[str] = set()
    not_accepted: list[str] = []
    for record in frames:
        frame_id = _get(record, "frame_id")
        if not isinstance(frame_id, str) or not frame_id:
            raise DatasetError(f"frame record without a frame_id: {record!r}"[:300])
        if frame_id in records:
            raise DatasetError(f"duplicate frame_id {frame_id!r} in split input")
        status = _outcome_status(record)
        if status is not None and status != "accepted":
            not_accepted.append(frame_id)
        group_id = _group_of(record, frame_id)
        run_id = _get(record, "run_id") or _field(record, "run_id") or frame_id.split("#", 1)[0]
        if not isinstance(group_id, str) or not group_id or _lineage_source(record) == "unresolved":
            unresolved.add(run_id)
        n_atoms = _get(record, "n_atoms")
        if n_atoms is not None and (isinstance(n_atoms, bool) or not isinstance(n_atoms, int) or n_atoms <= 0):
            raise DatasetError(f"frame {frame_id!r}: n_atoms must be a positive integer, got {n_atoms!r}")
        stratum = _field(record, stratify_by) if stratify_by else None
        records[frame_id] = {
            "group_id": group_id, "n_atoms": n_atoms, "stratum": stratum, "run_id": str(run_id),
            **{key: _field(record, key) for key in STRUCTURE_KEYS},
            "stats": {key: _label(_field(record, key)) for key in stats_keys},
        }
    if not_accepted:
        raise DatasetError(
            f"split input contains {len(not_accepted)} non-accepted frames (e.g. {not_accepted[:3]}); "
            "pass accepted frames only"
        )
    if unresolved:
        runs = sorted(unresolved)
        raise DatasetError(
            f"lineage is unresolved for {len(runs)} run(s) (e.g. {runs[:5]}); frames cannot be split safely. "
            "Provide an inventory with lineage ids or export with an explicit --lineage-policy "
            "(run|parent|depth:N)"
        )
    if not records:
        raise DatasetError("no frames to split")
    lineage_checks = _check_lineage_consistency(records)

    groups: dict[str, dict[str, Any]] = {}
    for frame_id in sorted(records):
        record = records[frame_id]
        group = groups.setdefault(record["group_id"], {"frame_ids": [], "atoms": 0, "strata": set()})
        group["frame_ids"].append(frame_id)
        group["atoms"] = None if group["atoms"] is None or record["n_atoms"] is None else group["atoms"] + record["n_atoms"]
        if stratify_by:
            group["strata"].add(canonical_json(record["stratum"]))
    if len(groups) < 3:
        raise DatasetError(
            f"only {len(groups)} independent lineage group(s); a leakage-free train/valid/test split needs at least 3"
        )
    strata_labels: dict[str, Any] = {}
    if stratify_by:
        missing = sorted(gid for gid, group in groups.items() if "null" in group["strata"])
        if missing:
            raise DatasetError(f"stratify_by={stratify_by!r}: value missing for groups {missing[:5]}")
        varying = sorted(gid for gid, group in groups.items() if len(group["strata"]) > 1)
        if varying:
            raise DatasetError(
                f"stratify_by={stratify_by!r} is not constant within group(s) {varying[:5]}; "
                "a group cannot be split, so it cannot belong to two strata"
            )
        for gid, group in groups.items():
            encoded = next(iter(group["strata"]))
            value = json.loads(encoded)
            label = _label(value)
            if label in strata_labels and canonical_json(strata_labels[label]) != encoded:
                raise DatasetError(f"stratum values {strata_labels[label]!r} and {value!r} have the same label")
            strata_labels[label] = value
            group["stratum"] = label
    else:
        for group in groups.values():
            group["stratum"] = None

    warnings: list[str] = []
    assignment: dict[str, str] = {}
    origin: dict[str, str] = {}
    previous_record = None
    new_fractions = dict(requested)
    if previous is not None:
        prev_frame_split, prev_group_split, prev_digest = _previous_index(previous)
        frozen = [] if grow_test else ["test"]
        for gid in sorted(groups):
            frame_ids = groups[gid]["frame_ids"]
            old_splits = {prev_frame_split[f] for f in frame_ids if f in prev_frame_split}
            if gid in prev_group_split:
                old_splits.add(prev_group_split[gid])
            if len(old_splits) > 1:
                examples = {split: [f for f in frame_ids if prev_frame_split.get(f) == split][:2] for split in sorted(old_splits)}
                raise LeakageError(
                    f"group {gid!r} now joins frames that the previous split placed in {sorted(old_splits)} "
                    f"(e.g. {examples}); previously separate groups were merged by new lineage/duplicate links. "
                    "Re-split from scratch (a new campaign) or remove the linking data"
                )
            if old_splits:
                split = old_splits.pop()
                gained = [f for f in frame_ids if f not in prev_frame_split]
                if split == "test" and gained and not grow_test:
                    raise DatasetError(
                        f"frozen test group {gid!r} gained {len(gained)} new frame(s) (e.g. {gained[:3]}); "
                        "pass grow_test=True (--grow-test) to extend the test set"
                    )
                if gained:
                    warnings.append(f"group {gid} kept in {split} and gained {len(gained)} frame(s)")
                assignment[gid] = split
                origin[gid] = "previous"
        current_frames = set(records)
        missing_by_split = {split: sorted(f for f, s in prev_frame_split.items() if s == split and f not in current_frames)
                            for split in SPLITS}
        for split in SPLITS:
            if missing_by_split[split]:
                warnings.append(
                    f"{len(missing_by_split[split])} frame(s) of the previous {split} split are absent now "
                    f"(e.g. {missing_by_split[split][:3]})"
                )
        if frozen:
            new_fractions = dict(requested)
            new_fractions["test"] = 0.0
            remaining = new_fractions["train"] + new_fractions["valid"]
            new_fractions = {name: value / remaining for name, value in new_fractions.items()}
        previous_fractions = previous.get("requested_fractions")
        if previous_fractions is not None and previous_fractions != requested:
            warnings.append(f"requested fractions {requested} differ from the previous split's {previous_fractions}")
        previous_record = {
            "content_sha256": prev_digest,
            "seed": previous.get("seed"),
            "frozen_splits": frozen,
            "groups_kept": sum(1 for value in origin.values() if value == "previous"),
            "groups_new": len(groups) - sum(1 for value in origin.values() if value == "previous"),
            # Previous group ids no longer present (their frames are absent, or they were renamed by merging).
            "previous_group_ids_absent": sorted(gid for gid in prev_group_split if gid not in groups),
            "frames_absent": {split: len(missing_by_split[split]) for split in SPLITS},
        }

    new_groups = sorted(gid for gid in groups if gid not in assignment)
    nonzero = [name for name in SPLITS if requested[name] > 0]
    strata_keys = sorted({groups[gid]["stratum"] for gid in groups}, key=lambda s: (s is None, s))
    underpopulated = []
    for stratum in strata_keys:
        members = [gid for gid in groups if groups[gid]["stratum"] == stratum]
        if stratify_by and len(members) < len(nonzero):
            underpopulated.append(stratum)
    if underpopulated:
        message = (f"stratify_by={stratify_by!r}: stratum/strata {underpopulated} have fewer groups than the "
                   f"{len(nonzero)} requested splits")
        if not allow_underpopulated_strata:
            raise DatasetError(message + "; pass allow_underpopulated_strata=True (--allow-underpopulated-strata) "
                               "to accept splits missing those strata")
        warnings.append(message + " (allowed)")
    for stratum in strata_keys:
        members = [gid for gid in new_groups if groups[gid]["stratum"] == stratum]
        if not members:
            continue
        sizes = {gid: len(groups[gid]["frame_ids"]) for gid in members}
        chosen = _assign_groups(sizes, new_fractions, _ranks(members, seed))
        for gid, split in chosen.items():
            assignment[gid] = split
            origin[gid] = "new"

    statistics, coverage = _statistics(records, groups, assignment, stats_keys, coverage_keys, nonzero)
    warnings.extend(coverage)
    manifest = _build_manifest(
        groups, assignment, origin, requested, new_fractions, seed, stratify_by, strata_labels,
        allow_underpopulated_strata, grow_test, previous_record, export_manifest_sha256, warnings, underpopulated,
        statistics, lineage_checks, stats_keys, coverage_keys,
    )
    check_split_disjoint(manifest)
    return manifest


def _statistics(records, groups, assignment, stats_keys, coverage_keys, nonzero):
    """Per key -> value -> split: frames and groups; plus coverage warnings."""
    group_of_frame = {frame_id: gid for gid, group in groups.items() for frame_id in group["frame_ids"]}
    table: dict[str, dict[str, dict[str, dict[str, Any]]]] = {}
    for key in stats_keys:
        rows: dict[str, dict[str, Any]] = {}
        for frame_id in sorted(records):
            value = records[frame_id]["stats"][key]
            gid = group_of_frame[frame_id]
            row = rows.setdefault(value, {split: {"frames": 0, "groups": set()} for split in SPLITS})
            cell = row[assignment[gid]]
            cell["frames"] += 1
            cell["groups"].add(gid)
        table[key] = {
            value: {split: {"frames": row[split]["frames"], "groups": len(row[split]["groups"])} for split in SPLITS}
            for value, row in sorted(rows.items())
        }
    warnings: list[str] = []
    split_groups = {split: sum(1 for gid in groups if assignment[gid] == split) for split in SPLITS}
    for key in coverage_keys:
        if key not in table:
            continue
        values = [value for value in table[key] if value != MISSING]
        for value in values:
            row = table[key][value]
            total_groups = sum(row[split]["groups"] for split in SPLITS)
            # A split with fewer groups than there are values cannot hold every value: not a warning.
            absent = [split for split in nonzero if row[split]["frames"] == 0 and split_groups[split] >= len(values)]
            if absent and total_groups >= len(nonzero):
                warnings.append(
                    f"coverage: {key}={value} ({total_groups} groups) is absent from {absent} "
                    "(consider --stratify-by " + key + ")"
                )
    return table, warnings


def _build_manifest(groups, assignment, origin, requested, new_fractions, seed, stratify_by, strata_labels,
                    allow_underpopulated, grow_test, previous_record, export_sha, warnings, underpopulated,
                    statistics, lineage_checks, stats_keys, coverage_keys):
    frame_ids = {split: [] for split in SPLITS}
    achieved = {split: {"groups": 0, "frames": 0, "atoms": 0} for split in SPLITS}
    atoms_known = all(group["atoms"] is not None for group in groups.values())
    for gid in sorted(groups):
        split = assignment[gid]
        frame_ids[split].extend(groups[gid]["frame_ids"])
        achieved[split]["groups"] += 1
        achieved[split]["frames"] += len(groups[gid]["frame_ids"])
        if atoms_known:
            achieved[split]["atoms"] += groups[gid]["atoms"]
    total_frames = sum(len(group["frame_ids"]) for group in groups.values())
    total_atoms = sum(group["atoms"] for group in groups.values()) if atoms_known else None
    for split in SPLITS:
        frame_ids[split].sort()
        row = achieved[split]
        row["group_fraction"] = row["groups"] / len(groups)
        row["frame_fraction"] = row["frames"] / total_frames
        if atoms_known:
            row["atom_fraction"] = row["atoms"] / total_atoms if total_atoms else 0.0
        else:
            row["atoms"] = None
            row["atom_fraction"] = None
        if abs(row["frame_fraction"] - requested[split]) > FRACTION_WARNING_TOLERANCE:
            warnings.append(
                f"{split}: achieved frame fraction {row['frame_fraction']:.4f} vs requested {requested[split]:.4f} "
                "(whole groups are indivisible)"
            )
        if requested[split] > 0 and row["groups"] == 0:
            warnings.append(f"{split}: requested fraction {requested[split]} but no group could be assigned")
    strata: dict[str, Any] = {}
    if stratify_by:
        for label in sorted(strata_labels):
            members = [gid for gid in groups if groups[gid]["stratum"] == label]
            per_split = {split: {"groups": 0, "frames": 0} for split in SPLITS}
            for gid in members:
                per_split[assignment[gid]]["groups"] += 1
                per_split[assignment[gid]]["frames"] += len(groups[gid]["frame_ids"])
            strata[label] = {
                "value": strata_labels[label],
                "groups": len(members),
                "frames": sum(len(groups[gid]["frame_ids"]) for gid in members),
                "per_split": per_split,
                "underpopulated": label in underpopulated,
            }
    manifest: dict[str, Any] = {
        "schema": SPLIT_SCHEMA,
        "schema_version": SPLIT_SCHEMA_VERSION,
        "algorithm": ALGORITHM,
        "rng": RNG,
        "seed": seed,
        "splits": list(SPLITS),
        "unit": GROUP_KEY,
        "requested_fractions": requested,
        "new_group_fractions": new_fractions,
        "stratify_by": stratify_by,
        "allow_underpopulated_strata": bool(allow_underpopulated),
        "grow_test": bool(grow_test),
        "export_manifest_sha256": export_sha,
        "previous": previous_record,
        "lineage_checks": lineage_checks,
        "totals": {"groups": len(groups), "frames": total_frames, "atoms": total_atoms},
        "achieved": achieved,
        "strata": strata,
        "stats_keys": list(stats_keys),
        "coverage_keys": list(coverage_keys),
        "statistics": statistics,
        "group_to_split": {gid: assignment[gid] for gid in sorted(groups)},
        "groups": {
            gid: {"split": assignment[gid], "frames": len(groups[gid]["frame_ids"]), "atoms": groups[gid]["atoms"],
                  "stratum": groups[gid]["stratum"], "origin": origin.get(gid, "new")}
            for gid in sorted(groups)
        },
        "frame_ids": frame_ids,
        "warnings": list(warnings),
    }
    manifest["content_sha256"] = sha256_bytes(canonical_json(_content(manifest)).encode())
    return manifest


def load_split_manifest(path: Path) -> dict[str, Any]:
    """Read a ``split_manifest.json``; verify schema, disjointness and ``content_sha256``."""
    path = Path(path)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise DatasetError(f"split manifest {path} does not exist") from None
    except json.JSONDecodeError as exc:
        raise DatasetError(f"split manifest {path} is not valid JSON: {exc}") from None
    if manifest.get("schema") != SPLIT_SCHEMA or manifest.get("schema_version") != SPLIT_SCHEMA_VERSION:
        raise DatasetError(f"{path.name} is not a {SPLIT_SCHEMA} v{SPLIT_SCHEMA_VERSION} manifest")
    expected = sha256_bytes(canonical_json(_content(manifest)).encode())
    if manifest.get("content_sha256") != expected:
        raise DatasetError(f"{path.name}: content_sha256 does not match its content (edited or corrupted)")
    check_split_disjoint(manifest)
    return manifest


# --------------------------------------------------------------------------
# Split records straight from dataset.extxyz (stdlib only)
# --------------------------------------------------------------------------

_SYMBOL_RE = re.compile(r"^[A-Z][a-z]?$")


def _composition(symbols: Sequence[str]) -> tuple[str, str]:
    counts: dict[str, int] = {}
    for symbol in symbols:
        counts[symbol] = counts.get(symbol, 0) + 1
    formula = "".join(f"{symbol}{counts[symbol]}" for symbol in sorted(counts))
    return formula, "-".join(sorted(counts))


def records_from_extxyz(path: Path) -> list[dict[str, Any]]:
    """Split records for every frame of an exported ``dataset.extxyz``.

    Reads the comment lines at text level (no ASE): ``frame_id``,
    ``lineage_group`` (and ``group_id`` if written), ``run_id``,
    ``structure_key``/``permutation_key`` and all other info as
    ``metadata``, plus ``n_atoms``, ``composition`` (e.g. ``Ni32O32``) and
    ``elements`` (``Ni-O``) from the atom lines. The export writes accepted
    frames only, so the records carry no outcome.
    """
    records = []
    for block in iter_extxyz_blocks(Path(path)):
        if block.frame_id is None:
            raise DatasetError(f"{Path(path).name}: frame {block.index} has no frame_id")
        info = block.info()
        if SPLIT_KEY in info:
            raise DatasetError(f"{Path(path).name}: frame {block.frame_id!r} already carries a split key; "
                               "split the export's dataset.extxyz, not a split file")
        symbols = [line.split(None, 1)[0] for line in block.text.split("\n")[2:2 + block.n_atoms]]
        if not all(_SYMBOL_RE.match(symbol) for symbol in symbols):
            raise DatasetError(f"{Path(path).name}: frame {block.frame_id!r} has malformed species")
        composition, elements = _composition(symbols)
        metadata = dict(info)
        metadata.setdefault("composition", composition)
        metadata.setdefault("elements", elements)
        record = {"frame_id": block.frame_id, "n_atoms": block.n_atoms, "metadata": metadata}
        for key in (GROUP_KEY, "group_id", "run_id", *STRUCTURE_KEYS):
            if info.get(key) is not None:
                record[key] = info[key]
        records.append(record)
    return records


# --------------------------------------------------------------------------
# Writing split files
# --------------------------------------------------------------------------

def _verify_split_file(path: Path, split: str, expected_ids: Sequence[str]) -> None:
    """Re-read a split file with ASE: same frames, split key, no calculator/constraints."""
    try:
        import ase.io
    except ImportError:  # pragma: no cover - ASE is part of the dataset extra
        return
    if not expected_ids:
        if path.read_bytes():
            raise DatasetError(f"{path.name}: expected an empty file")
        return
    got = []
    for atoms in ase.io.iread(str(path), index=":", format="extxyz"):
        if atoms.info.get(SPLIT_KEY) != split:
            raise DatasetError(f"{path.name}: frame {atoms.info.get('frame_id')!r} has split {atoms.info.get(SPLIT_KEY)!r}")
        if atoms.calc is not None or atoms.constraints:
            raise DatasetError(f"{path.name}: ASE attached a calculator or constraints to {atoms.info.get('frame_id')!r}")
        got.append(atoms.info.get("frame_id"))
    if got != list(expected_ids):
        raise DatasetError(f"{path.name}: ASE read {len(got)} frames, expected {len(expected_ids)} in dataset order")


def _summary_markdown(manifest: Mapping[str, Any]) -> str:
    lines = [
        "# Dataset split summary",
        "",
        f"- schema: `{manifest['schema']}` v{manifest['schema_version']}, algorithm `{manifest['algorithm']}`,"
        f" rng `{manifest['rng']}`, seed {manifest['seed']}",
        f"- unit: `{manifest['unit']}` (all frames of a lineage group share one split)",
        f"- stratify_by: `{manifest['stratify_by']}`; grow_test: {manifest['grow_test']}",
        f"- content_sha256: `{manifest['content_sha256']}`",
        "",
        "| split | requested | groups | frames | frame fraction |",
        "|---|---|---|---|---|",
    ]
    for split in SPLITS:
        row = manifest["achieved"][split]
        lines.append(f"| {split} | {manifest['requested_fractions'][split]:.4f} | {row['groups']} | {row['frames']} |"
                     f" {row['frame_fraction']:.4f} |")
    for key, values in manifest.get("statistics", {}).items():
        if list(values) == [MISSING]:
            continue
        lines += ["", f"## {key}", "", "| value | " + " | ".join(f"{s} frames (groups)" for s in SPLITS) + " |",
                  "|---|" + "---|" * len(SPLITS)]
        for value, row in values.items():
            cells = " | ".join(f"{row[s]['frames']} ({row[s]['groups']})" for s in SPLITS)
            lines.append(f"| {value} | {cells} |")
    files = manifest.get("files", {})
    if files:
        lines += ["", "## Files", "", "| file | frames | sha256 |", "|---|---|---|"]
        for split in SPLITS:
            lines.append(f"| {files[split]['path']} | {files[split]['frames']} | `{files[split]['sha256']}` |")
    lines += ["", "## Warnings", ""]
    lines += [f"- {warning}" for warning in manifest.get("warnings", [])] or ["- none"]
    return "\n".join(lines) + "\n"


def write_split(
    export_dir: Path,
    manifest: Mapping[str, Any],
    output_dir: Path,
    *,
    force: bool = False,
    dataset_name: str = "dataset.extxyz",
    verify: bool = True,
) -> dict[str, Any]:
    """Write ``train/valid/test.extxyz``, ``split_manifest.json`` and ``split_summary.md``.

    The split files are text slices of ``export_dir/dataset.extxyz`` (frames
    in dataset order, bytes unchanged except the inserted ``split="<name>"``
    info key). The split must cover exactly the exported frames. Output goes
    to a new (or ``force``-d) directory other than the export directory; with
    ``verify`` every split file is re-read with ASE (frame ids, split key, no
    calculator/constraints). ``split_manifest.json`` is written last, with
    each split file's frame count and sha256, the dataset and export-manifest
    hashes (outside ``content_sha256``). Repeated calls give identical bytes.
    """
    export_dir = Path(export_dir)
    output_dir = Path(output_dir)
    dataset = export_dir / dataset_name
    if not dataset.is_file():
        raise DatasetError(f"{dataset} does not exist; run 'dataset export' first")
    if output_dir.resolve() == export_dir.resolve():
        raise DatasetError("write the split to its own directory, not into the export directory")
    check_split_disjoint(manifest)
    requested = {frame_id for split in SPLITS for frame_id in manifest.get("frame_ids", {}).get(split, [])}
    present = read_frame_ids(dataset)
    if any(frame_id is None for frame_id in present):
        raise DatasetError(f"{dataset.name} contains frames without frame_id")
    present_set = set(present)
    if len(present_set) != len(present):
        raise DatasetError(f"{dataset.name} repeats frame ids")
    missing = sorted(requested - present_set)
    unassigned = sorted(present_set - requested)
    if missing or unassigned:
        raise DatasetError(
            f"split manifest does not match {dataset.name}: {len(missing)} split frame(s) absent from the export "
            f"(e.g. {missing[:3]}), {len(unassigned)} exported frame(s) in no split (e.g. {unassigned[:3]})"
        )
    export_manifest = export_dir / "dataset_manifest.json"
    export_sha = sha256_file(export_manifest) if export_manifest.is_file() else None
    recorded = manifest.get("export_manifest_sha256")
    if recorded is not None and export_sha is not None and recorded != export_sha:
        raise DatasetError(
            f"split manifest was made for export manifest {recorded[:12]}..., but {export_dir} has {export_sha[:12]}..."
        )
    prepare_output_dir(output_dir, force=force)
    files: dict[str, Any] = {}
    order = {frame_id: index for index, frame_id in enumerate(present)}
    for split in SPLITS:
        target = output_dir / f"{split}.extxyz"
        result = slice_extxyz_by_frame_ids(dataset, manifest["frame_ids"][split], target, add_info={SPLIT_KEY: split})
        if verify:
            _verify_split_file(target, split, sorted(manifest["frame_ids"][split], key=order.__getitem__))
        files[split] = {"path": f"{split}.extxyz", "frames": result["frames"], "sha256": result["sha256"]}
    written = dict(manifest)
    written["files"] = files
    written["dataset_extxyz_sha256"] = sha256_file(dataset)
    written["export_manifest_file_sha256"] = export_sha
    summary = output_dir / "split_summary.md"
    atomic_write_text(summary, _summary_markdown(written))
    written["summary_file"] = {"path": summary.name, "sha256": sha256_file(summary)}
    atomic_write_json(output_dir / "split_manifest.json", written)
    return written


def split_export(
    export_dir: Path,
    output_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    fractions: Sequence[float] = DEFAULT_FRACTIONS,
    stratify_by: str | None = None,
    previous: Path | Mapping[str, Any] | None = None,
    grow_test: bool = False,
    allow_underpopulated_strata: bool = False,
    force: bool = False,
    dataset_name: str = "dataset.extxyz",
    stats_keys: Sequence[str] = DEFAULT_STATS_KEYS,
    coverage_keys: Sequence[str] = DEFAULT_COVERAGE_KEYS,
    verify: bool = True,
) -> dict[str, Any]:
    """``dataset split``: records from ``dataset.extxyz`` -> :func:`make_split` -> :func:`write_split`.

    ``previous`` is a ``split_manifest.json`` path (verified with
    :func:`load_split_manifest`) or an already loaded manifest.
    """
    export_dir = Path(export_dir)
    dataset = export_dir / dataset_name
    if not dataset.is_file():
        raise DatasetError(f"{dataset} does not exist; run 'dataset export' first")
    if previous is not None and not isinstance(previous, Mapping):
        previous = load_split_manifest(Path(previous))
    export_manifest = export_dir / "dataset_manifest.json"
    manifest = make_split(
        records_from_extxyz(dataset),
        seed=seed, fractions=fractions, stratify_by=stratify_by, previous=previous, grow_test=grow_test,
        allow_underpopulated_strata=allow_underpopulated_strata,
        export_manifest_sha256=sha256_file(export_manifest) if export_manifest.is_file() else None,
        stats_keys=stats_keys, coverage_keys=coverage_keys,
    )
    return write_split(export_dir, manifest, output_dir, force=force, dataset_name=dataset_name, verify=verify)
