"""Discover VASP run directories under explicit, aliased calculation roots.

A *run directory* is a directory that directly contains a run-defining file:
``vasprun.xml`` (or ``vasprun.xml.gz|.bz2|.xz``), ``OUTCAR`` (or a compressed
OUTCAR), ``OSZICAR`` or ``vaspout.h5``. The label source is the vasprun
(plain preferred over compressed; the other copy is reported as
``duplicate_label_file``). Directories with outputs but no vasprun are still
returned, with ``label_file=None`` and the flag ``no_label_source`` (plus
``outpacklite_like`` for OSZICAR/XDATCAR-only packages), so the exporter
records them as ``missing_labels/no_vasprun``
(:meth:`DiscoveredRun.missing_label_outcome`) instead of dropping them.

Nothing is dropped silently: every pruned directory, every run candidate
inside a pruned directory, every non-canonical near-name (``vasprun.xml.1``,
``OUTCAR.bak``), every duplicate label/evidence file, every symlink not
followed, the output directory and every unreadable directory is returned in
the ``ignored`` list with a reason code (:data:`IGNORE_REASONS`).

Traversal is deterministic (entries sorted by name; results sorted by
``(alias, relpath)``), independent of root order and filesystem enumeration
order, never follows directory symlinks unless asked, keeps a visited set of
real paths (so junctions/symlink loops cannot recurse), and never descends
into the output directory.

The default exclusion list merges the Phase 1 spec with InterfaceForge tree
conventions (lgutsev/InterfaceForge@4501e34, docs; see
``design/understand/if-phase1-map.md`` section 2.8): ``.interfaceforge``
repair/rewind archives, ``precondition`` NSW=0 preconditioners, ``*archive*``,
``*backup*`` and upper-case ``X*`` (OutPackLite) parts, plus the agglomeration
stale-artifact names. The rules are data (:data:`DEFAULT_EXCLUSION_RULES`),
are recorded in manifests via :func:`exclusion_rules_manifest`, and a rule can
be disabled deliberately (``include_parts``; every run found that way carries
``included_excluded_part:<pattern>``). Run candidates inside excluded trees
become excluded run records via :func:`pruned_run_candidates`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import fnmatch
import os
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

from . import model
from .errors import DatasetError
from .fsio import sha256_file

# --------------------------------------------------------------------------
# File-name vocabulary
# --------------------------------------------------------------------------

#: Label sources in order of preference; value = ``label_kind``.
LABEL_FILES: tuple[tuple[str, str], ...] = (
    ("vasprun.xml", "vasprun"),
    ("vasprun.xml.gz", "vasprun.gz"),
    ("vasprun.xml.bz2", "vasprun.bz2"),
    ("vasprun.xml.xz", "vasprun.xz"),
)
_COMPRESSION_SUFFIXES = ("", ".gz", ".bz2", ".xz")
#: Evidence kinds -> accepted file names in order of preference.
EVIDENCE_FILES: dict[str, tuple[str, ...]] = {
    "OUTCAR": tuple("OUTCAR" + suffix for suffix in _COMPRESSION_SUFFIXES),
    "OSZICAR": tuple("OSZICAR" + suffix for suffix in _COMPRESSION_SUFFIXES),
    "XDATCAR": tuple("XDATCAR" + suffix for suffix in _COMPRESSION_SUFFIXES),
    "INCAR": ("INCAR",),
    "POTCAR": ("POTCAR",),
    "POSCAR": ("POSCAR",),
    "CONTCAR": ("CONTCAR",),
    "KPOINTS": ("KPOINTS",),
    "vaspout.h5": ("vaspout.h5",),
}
#: Evidence kinds whose presence makes a directory a run directory.
RUN_DEFINING_EVIDENCE = ("OUTCAR", "OSZICAR", "vaspout.h5")
#: Provenance/lineage files recorded when present in a run directory.
METADATA_FILES = (
    "agglomeration_manifest.json", "training_frame.json", "sanity_report.json",
    "provenance.json", "step1_repair.json", "opt_manifest.json", "step1_manifest.json",
    "step2_manifest.json",
)
_KNOWN_NAMES = {name for name, _ in LABEL_FILES} | {n for names in EVIDENCE_FILES.values() for n in names}
#: Prefixes (case-insensitive) of names that look like label/evidence files.
_NEAR_NAME_PREFIXES = ("vasprun", "outcar", "oszicar")

#: Discovery-level reason codes for ``ignored`` entries. Run-level codes that
#: already exist in :data:`model.REASONS` (``archived_path``,
#: ``preconditioning_run``) are reused for pruned directories.
IGNORE_REASONS: dict[str, str] = {
    "archived_path": "path matches an archived/stale-artifact or non-VASP-stage pattern",
    "preconditioning_run": "static wavefunction-preconditioning calculation directory",
    "user_excluded": "matched a user --exclude-glob",
    "output_dir": "the command's own output directory (never scanned)",
    "duplicate_label_file": "another copy of the label source exists in the same directory (plain vasprun preferred)",
    "duplicate_evidence_file": "another copy of this evidence file exists in the same directory (plain file preferred)",
    "not_a_label_file": "name resembles a VASP output but is not a canonical file name (e.g. vasprun.xml.1, OUTCAR.bak)",
    "symlink_not_followed": "directory symlink/junction not followed (follow_symlinks=False)",
    "symlink_loop": "directory already visited through another path (symlink loop or alias)",
    "broken_symlink": "symlink target does not exist",
    "unreadable_dir": "directory could not be listed",
    "no_vasp_outputs": "VASP input files but no vasprun.xml/OUTCAR/OSZICAR/vaspout.h5 (never run, or outputs removed)",
}


# --------------------------------------------------------------------------
# Exclusion rules
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class ExclusionRule:
    """A default rule applied to each path component below a root."""

    pattern: str  # fnmatch pattern for ONE path component
    reason: str  # archived_path | preconditioning_run
    note: str
    case_sensitive: bool = False
    descend: bool = True  # False: do not even look inside (VCS/cache internals)

    @property
    def rule_id(self) -> str:
        return f"part:{self.pattern}" + (":cs" if self.case_sensitive else "")

    def matches(self, name: str) -> bool:
        if self.case_sensitive:
            return fnmatch.fnmatchcase(name, self.pattern)
        return fnmatch.fnmatchcase(name.lower(), self.pattern.lower())

    def as_dict(self) -> dict[str, Any]:
        return {"rule": self.rule_id, "pattern": self.pattern, "reason": self.reason, "note": self.note,
                "case_sensitive": self.case_sensitive, "descend": self.descend}


DEFAULT_EXCLUSION_RULES: tuple[ExclusionRule, ...] = (
    ExclusionRule("precondition", "preconditioning_run",
                  "InterfaceForge NSW=0 preconditioner (duplicates the MD frame 0) or empty Step2 artefact"),
    ExclusionRule(".git", "archived_path", "version-control internals", descend=False),
    ExclusionRule("__pycache__", "archived_path", "Python cache", descend=False),
    ExclusionRule(".interfaceforge", "archived_path", "InterfaceForge repair/rewind archives (crash frames)"),
    ExclusionRule("xtb", "archived_path", "xTB stage files; xTB labels are never VASP labels"),
    ExclusionRule("rescue_xtb", "archived_path", "xTB rescue stage files"),
    ExclusionRule("packmol", "archived_path", "Packmol inputs/outputs, not calculations"),
    ExclusionRule("templates", "archived_path", "molecule templates, not calculations"),
    ExclusionRule("structures", "archived_path", "structures-only case tree (vasp_ready=false)"),
    ExclusionRule("vasp_reference", "archived_path", "copied reference/template VASP files, not a run"),
    ExclusionRule("*archive*", "archived_path",
                  "archive, archived*, *_archived*, restart_/refit_/stability_archive_* (IF rewinds)"),
    ExclusionRule("*backup*", "archived_path", "backup copies"),
    ExclusionRule("old", "archived_path", "stale copies"),
    ExclusionRule("old_*", "archived_path", "stale copies"),
    ExclusionRule("trash", "archived_path", "deleted/stale data"),
    ExclusionRule("*.bak", "archived_path", "backup copies"),
    ExclusionRule("*.protocol_*", "archived_path", "agglomeration stale xTB protocol artefact"),
    ExclusionRule("*.input_*", "archived_path", "agglomeration stale xTB input artefact"),
    ExclusionRule("*.geometry_*", "archived_path", "agglomeration stale Packmol geometry artefact"),
    ExclusionRule("X*", "archived_path",
                  "OutPackLite / packaged legacy trees (IF convention; upper-case X only)", case_sensitive=True),
)


def check_include_parts(include_parts: Iterable[str]) -> tuple[str, ...]:
    """Validate deliberately included default-excluded parts; returns canonical rule patterns.

    Each value must name one :data:`DEFAULT_EXCLUSION_RULES` pattern (case-
    insensitive), e.g. ``precondition``, ``.interfaceforge``, ``X*`` or
    ``*archive*``. Unknown names are an error (a typo must not silently keep
    a tree excluded or included).
    """
    by_name = {rule.pattern.lower(): rule.pattern for rule in DEFAULT_EXCLUSION_RULES}
    chosen: list[str] = []
    for value in include_parts:
        pattern = by_name.get(str(value).strip().lower())
        if pattern is None:
            raise DatasetError(
                f"--include-excluded-part {value!r} names no default exclusion rule; "
                f"choose from {[rule.pattern for rule in DEFAULT_EXCLUSION_RULES]}"
            )
        if pattern not in chosen:
            chosen.append(pattern)
    return tuple(sorted(chosen))


def exclusion_rules_manifest(
    exclude_globs: Sequence[str] = (),
    include_globs: Sequence[str] = (),
    include_parts: Sequence[str] = (),
) -> dict[str, Any]:
    """The effective discovery rules, for recording verbatim in manifests."""
    parts = check_include_parts(include_parts)
    return {
        "default_rules": [dict(rule.as_dict(), deliberately_included=rule.pattern in parts)
                          for rule in DEFAULT_EXCLUSION_RULES],
        "exclude_globs": list(exclude_globs),
        "include_globs": list(include_globs),
        "include_parts": list(parts),
        "precedence": "user exclude-glob > user include-glob > deliberately included part > default rules",
        "glob_semantics": "fnmatchcase against the posix relpath below the root and against '<alias>:<relpath>'; '*' crosses '/'",
        "label_files": [name for name, _ in LABEL_FILES],
        "run_defining_evidence": list(RUN_DEFINING_EVIDENCE),
    }


# --------------------------------------------------------------------------
# Roots
# --------------------------------------------------------------------------

_ALIAS_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


@dataclass(frozen=True)
class RootSpec:
    alias: str
    path: Path  # resolved absolute directory

    def as_dict(self) -> dict[str, str]:
        return {"alias": self.alias, "path": self.path.as_posix()}


def _check_alias(alias: str, source: str) -> str:
    if not _ALIAS_RE.match(alias):
        raise DatasetError(
            f"root alias {alias!r} (from {source!r}) must match [A-Za-z0-9][A-Za-z0-9_.-]*; "
            "pass an explicit alias as alias=path"
        )
    return alias


def parse_roots(values: Iterable[str | Path | RootSpec]) -> list[RootSpec]:
    """Parse ``alias=path`` or ``path`` root arguments.

    The default alias is the directory's basename. Aliases must be unique and
    safe for ids (no ``:``, ``/``, ``=``, whitespace); paths must exist, be
    directories, be distinct and must not contain one another (a nested root
    would be scanned twice). Order is preserved; discovery output does not
    depend on it.
    """
    specs: list[RootSpec] = []
    for value in values:
        if isinstance(value, RootSpec):
            alias, raw = value.alias, str(value.path)
        else:
            text = str(value)
            alias = None
            raw = text
            if isinstance(value, str) and "=" in text:
                candidate, rest = text.split("=", 1)
                if _ALIAS_RE.match(candidate) and rest:
                    alias, raw = candidate, rest
            if alias is None and isinstance(value, str) and text.startswith("="):
                raise DatasetError(f"root {text!r} has an empty alias")
        path = Path(raw).expanduser()
        if not path.exists():
            raise DatasetError(f"calculation root {raw!r} does not exist")
        if not path.is_dir():
            raise DatasetError(f"calculation root {raw!r} is not a directory")
        path = path.resolve()
        if alias is None:
            if not path.name:
                raise DatasetError(f"root {raw!r} has no basename; pass alias=path")
            alias = path.name
        specs.append(RootSpec(_check_alias(alias, str(value)), path))
    if not specs:
        raise DatasetError("at least one calculation root is required")
    aliases = [spec.alias for spec in specs]
    duplicates = sorted({alias for alias in aliases if aliases.count(alias) > 1})
    if duplicates:
        raise DatasetError(f"duplicate root aliases {duplicates}; pass explicit unique aliases as alias=path")
    for index, first in enumerate(specs):
        for second in specs[index + 1:]:
            if first.path == second.path:
                raise DatasetError(f"roots {first.alias!r} and {second.alias!r} are the same directory {first.path}")
            if _is_relative_to(first.path, second.path) or _is_relative_to(second.path, first.path):
                raise DatasetError(
                    f"roots {first.alias!r} ({first.path}) and {second.alias!r} ({second.path}) overlap; "
                    "a run would be discovered twice"
                )
    return specs


def _is_relative_to(path: Path, other: Path) -> bool:
    try:
        path.relative_to(other)
    except ValueError:
        return False
    return True


# --------------------------------------------------------------------------
# Discovered runs
# --------------------------------------------------------------------------

@dataclass
class DiscoveredRun:
    root_alias: str
    root_path: Path
    run_dir: Path
    relpath: str  # posix, relative to root_path; "." for the root itself
    label_file: Path | None
    label_kind: str | None  # vasprun | vasprun.gz | vasprun.bz2 | vasprun.xz | None
    evidence: dict[str, Path] = field(default_factory=dict)  # kind -> chosen file
    metadata_files: dict[str, Path] = field(default_factory=dict)
    ignored_candidates: list[dict[str, Any]] = field(default_factory=list)
    nested_in: str | None = None  # relpath of the enclosing run directory, if any
    flags: list[str] = field(default_factory=list)  # nested_run, included_by_glob, symlinked_label_file, ...

    @property
    def run_id(self) -> str:
        return f"{self.root_alias}:{self.relpath}"

    @property
    def has_label_source(self) -> bool:
        return self.label_file is not None

    def missing_label_outcome(self) -> "model.Outcome | None":
        """``missing_labels/no_vasprun`` for a run without vasprun (``None`` when one exists).

        OUTCAR-only, OSZICAR-only (e.g. OutPackLite packages: OSZICAR,
        XDATCAR, INCAR, CONTCAR) and vaspout.h5-only directories are recorded,
        never turned into labels: OSZICAR energies have ~8 significant digits
        and no forces, and OUTCAR-only extraction is not implemented.
        """
        if self.label_file is not None:
            return None
        present = sorted(path.name for path in self.evidence.values())
        detail = "outputs present: " + (", ".join(present) if present else "none")
        if "OSZICAR" in self.evidence and "OUTCAR" not in self.evidence:
            detail += "; OSZICAR/XDATCAR carry no forces and ~8-digit energies: never used as labels"
        return model.outcome_for("no_vasprun", detail)

    def as_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "root_alias": self.root_alias,
            "relpath": self.relpath,
            "label_file": self.label_file.name if self.label_file else None,
            "label_kind": self.label_kind,
            "evidence": {kind: path.name for kind, path in sorted(self.evidence.items())},
            "metadata_files": sorted(self.metadata_files),
            "nested_in": self.nested_in,
            "flags": list(self.flags),
            "ignored_candidates": [dict(item) for item in self.ignored_candidates],
        }


def _join(relpath: str, name: str) -> str:
    return name if relpath == "." else f"{relpath}/{name}"


def _ignored(alias: str, relpath: str, kind: str, reason: str, detail: str = "", rule: str | None = None) -> dict[str, Any]:
    if reason not in IGNORE_REASONS:
        raise AssertionError(f"unknown ignore reason {reason}")
    return {"root_alias": alias, "relpath": relpath, "path": f"{alias}:{relpath}", "kind": kind,
            "reason": reason, "rule": rule, "detail": detail}


def _glob_match(patterns: Sequence[str], alias: str, relpath: str) -> str | None:
    for pattern in patterns:
        if fnmatch.fnmatchcase(relpath, pattern) or fnmatch.fnmatchcase(f"{alias}:{relpath}", pattern):
            return pattern
    return None


def _default_rule(name: str, disabled: Sequence[str] = ()) -> ExclusionRule | None:
    """First default rule matching ``name`` that was not deliberately included."""
    for rule in DEFAULT_EXCLUSION_RULES:
        if rule.pattern not in disabled and rule.matches(name):
            return rule
    return None


def _included_parts(name: str, disabled: Sequence[str]) -> list[str]:
    """Deliberately included default rules that ``name`` matches."""
    return [rule.pattern for rule in DEFAULT_EXCLUSION_RULES if rule.pattern in disabled and rule.matches(name)]


def _is_link(entry: os.DirEntry) -> bool:
    if entry.is_symlink():
        return True
    is_junction = getattr(entry, "is_junction", None)
    return bool(is_junction and is_junction())


@dataclass
class _Walk:
    root: RootSpec
    exclude_globs: Sequence[str]
    include_globs: Sequence[str]
    include_parts: Sequence[str]
    output_dir: Path | None
    follow_symlinks: bool
    visited: set[str]
    runs: list[DiscoveredRun]
    ignored: list[dict[str, Any]]


def _classify_files(walk: _Walk, directory: Path, relpath: str, files: list[os.DirEntry]):
    """Pick label/evidence/metadata files; return (label, kind, evidence, metadata, local ignored, flags)."""
    alias = walk.root.alias
    names = {entry.name: entry for entry in files}
    local: list[dict[str, Any]] = []
    flags: list[str] = []
    label_path = None
    label_kind = None
    for name, kind in LABEL_FILES:
        if name not in names:
            continue
        if label_path is None:
            label_path, label_kind = directory / name, kind
            if _is_link(names[name]):
                flags.append("symlinked_label_file")
        else:
            local.append(_ignored(alias, _join(relpath, name), "file", "duplicate_label_file",
                                  f"label source is {label_path.name}"))
    evidence: dict[str, Path] = {}
    for kind, candidates in EVIDENCE_FILES.items():
        present = [name for name in candidates if name in names]
        if not present:
            continue
        evidence[kind] = directory / present[0]
        for extra in present[1:]:
            local.append(_ignored(alias, _join(relpath, extra), "file", "duplicate_evidence_file",
                                  f"{kind} evidence is {present[0]}"))
    metadata = {name: directory / name for name in METADATA_FILES if name in names}
    for name in sorted(names):
        if name in _KNOWN_NAMES:
            continue
        if name.lower().startswith(_NEAR_NAME_PREFIXES):
            local.append(_ignored(alias, _join(relpath, name), "file", "not_a_label_file",
                                  "only exact vasprun.xml[.gz|.bz2|.xz] / OUTCAR[...] / OSZICAR[...] names are read"))
    if label_path is None and _is_run_dir(label_path, evidence):
        flags.append("no_label_source")
        if "OSZICAR" in evidence and "OUTCAR" not in evidence:
            # OutPackLite / IF package_outputs layout: OSZICAR (+XDATCAR, INCAR, CONTCAR) only.
            flags.append("outpacklite_like")
    return label_path, label_kind, evidence, metadata, local, flags


def _is_run_dir(label_path, evidence) -> bool:
    return label_path is not None or any(kind in evidence for kind in RUN_DEFINING_EVIDENCE)


def _scan_dir(walk: _Walk, directory: Path, relpath: str):
    try:
        with os.scandir(directory) as iterator:
            entries = sorted(iterator, key=lambda entry: entry.name)
    except OSError as exc:
        walk.ignored.append(_ignored(walk.root.alias, relpath, "dir", "unreadable_dir", f"{type(exc).__name__}: {exc}"))
        return None
    files: list[os.DirEntry] = []
    dirs: list[os.DirEntry] = []
    for entry in entries:
        try:
            if entry.is_dir(follow_symlinks=True):
                dirs.append(entry)
            elif entry.is_file(follow_symlinks=True):
                files.append(entry)
            elif _is_link(entry):
                walk.ignored.append(_ignored(walk.root.alias, _join(relpath, entry.name), "symlink", "broken_symlink"))
        except OSError as exc:
            walk.ignored.append(_ignored(walk.root.alias, _join(relpath, entry.name), "dir", "unreadable_dir", str(exc)))
    return files, dirs


def _walk(
    walk: _Walk,
    directory: Path,
    relpath: str,
    parent_run: str | None,
    pruned: dict[str, Any] | None,
    inclusion: tuple[str, ...],
) -> None:
    """Depth-first traversal.

    ``pruned`` is the ignored-entry of an excluded ancestor (run candidates
    below it are reported, not returned); ``inclusion`` holds the run flags
    that record why this subtree is scanned although a default rule would
    have excluded it (``included_by_glob:<glob>``,
    ``included_excluded_part:<pattern>``).
    """
    real = os.path.realpath(directory)
    if real in walk.visited:
        walk.ignored.append(_ignored(walk.root.alias, relpath, "dir", "symlink_loop", f"real path {Path(real).as_posix()}"))
        return
    walk.visited.add(real)
    scanned = _scan_dir(walk, directory, relpath)
    if scanned is None:
        return
    files, dirs = scanned
    alias = walk.root.alias
    label_path, label_kind, evidence, metadata, local, flags = _classify_files(walk, directory, relpath, files)
    this_run = parent_run
    if _is_run_dir(label_path, evidence):
        if pruned is not None:
            outputs = sorted({path.name for path in evidence.values()} | ({label_path.name} if label_path else set()))
            walk.ignored.append(_ignored(alias, relpath, "run_dir", pruned["reason"],
                                         f"run candidate ({', '.join(outputs)}) inside excluded {pruned['path']}",
                                         pruned["rule"]))
        else:
            run = DiscoveredRun(alias, walk.root.path, directory, relpath, label_path, label_kind, evidence,
                                metadata, local, nested_in=parent_run, flags=list(flags))
            if parent_run is not None:
                run.flags.append("nested_run")
            run.flags.extend(flag for flag in inclusion if flag not in run.flags)
            walk.runs.append(run)
            walk.ignored.extend(local)
            this_run = relpath
    elif pruned is None:
        walk.ignored.extend(local)
        if "INCAR" in evidence:
            walk.ignored.append(_ignored(alias, relpath, "dir", "no_vasp_outputs",
                                         "inputs present: " + ", ".join(sorted(p.name for p in evidence.values()))))

    for entry in dirs:
        child_rel = _join(relpath, entry.name)
        child = directory / entry.name
        if _is_link(entry) and not walk.follow_symlinks:
            if pruned is None:
                target = Path(os.path.realpath(child)).as_posix()
                walk.ignored.append(_ignored(alias, child_rel, "symlink", "symlink_not_followed", f"target {target}"))
            continue
        if walk.output_dir is not None and Path(os.path.realpath(child)) == walk.output_dir:
            walk.ignored.append(_ignored(alias, child_rel, "dir", "output_dir"))
            continue
        user_ex = _glob_match(walk.exclude_globs, alias, child_rel)
        user_in = _glob_match(walk.include_globs, alias, child_rel)
        rule = _default_rule(entry.name, walk.include_parts)
        parts = tuple(f"included_excluded_part:{pattern}" for pattern in _included_parts(entry.name, walk.include_parts))
        if user_ex is not None:
            # A user exclude always wins, also over include-globs further down.
            record = pruned
            if pruned is None:
                record = _ignored(alias, child_rel, "dir", "user_excluded", f"--exclude-glob {user_ex}",
                                  f"exclude-glob:{user_ex}")
                walk.ignored.append(record)
            _walk(walk, child, child_rel, this_run, record, ())
        elif pruned is not None:
            if user_in is not None and pruned["reason"] != "user_excluded":
                _walk(walk, child, child_rel, this_run, None, parts + (f"included_by_glob:{user_in}",))
            elif rule is not None and not rule.descend:
                continue  # e.g. .git inside an already reported archive
            else:
                _walk(walk, child, child_rel, this_run, pruned, ())
        elif rule is not None and user_in is None:
            record = _ignored(alias, child_rel, "dir", rule.reason, rule.note, rule.rule_id)
            walk.ignored.append(record)
            if rule.descend:
                _walk(walk, child, child_rel, this_run, record, ())
        else:
            extra = parts + ((f"included_by_glob:{user_in}",) if rule is not None else ())
            _walk(walk, child, child_rel, this_run, None, inclusion + tuple(f for f in extra if f not in inclusion))


def discover_runs(
    roots: Sequence[RootSpec],
    *,
    exclude_globs: Sequence[str] = (),
    include_globs: Sequence[str] = (),
    include_parts: Sequence[str] = (),
    output_dir: Path | None = None,
    follow_symlinks: bool = False,
) -> tuple[list[DiscoveredRun], list[dict[str, Any]]]:
    """Find run directories under ``roots``.

    Returns ``(runs, ignored)``: runs sorted by ``(root_alias, relpath)``;
    ``ignored`` entries ``{root_alias, relpath, path, kind, reason, rule,
    detail}`` sorted by ``(root_alias, relpath, kind, reason)`` (``kind`` in
    dir|run_dir|file|symlink). ``run_dir`` entries are run candidates inside a
    pruned directory; :func:`pruned_run_candidates` turns them into excluded
    run records so no candidate calculation disappears from the audit.

    Precedence: a user exclude-glob always wins; a user include-glob
    overrides the default rules (also for a directory deep inside a
    default-excluded tree); ``include_parts`` deliberately disables named
    default rules (e.g. ``precondition``, ``.interfaceforge``, ``X*``) and
    tags every run found through them ``included_excluded_part:<pattern>``;
    default rules apply to path components below the root only (the root
    itself is always scanned).
    """
    if not roots:
        raise DatasetError("at least one calculation root is required")
    aliases = [root.alias for root in roots]
    if len(set(aliases)) != len(aliases):
        raise DatasetError(f"duplicate root aliases {aliases}")
    disabled = check_include_parts(include_parts)
    out = Path(os.path.realpath(Path(output_dir))) if output_dir is not None else None
    runs: list[DiscoveredRun] = []
    ignored: list[dict[str, Any]] = []
    for root in sorted(roots, key=lambda spec: spec.alias):
        if out is not None and (out == root.path or _is_relative_to(root.path, out)):
            raise DatasetError(f"root {root.alias!r} ({root.path}) lies inside the output directory {out}")
        walk = _Walk(root, tuple(exclude_globs), tuple(include_globs), disabled, out, follow_symlinks, set(), [], [])
        _walk(walk, root.path, ".", None, None, ())
        runs.extend(walk.runs)
        ignored.extend(walk.ignored)
    runs.sort(key=lambda run: (run.root_alias, run.relpath))
    ignored.sort(key=lambda item: (item["root_alias"], item["relpath"], item["kind"], item["reason"]))
    return runs, ignored


def pruned_run_candidates(ignored: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Run records for every run candidate inside an excluded directory.

    Each ``ignored`` entry of kind ``run_dir`` becomes ``{run_id, root_alias,
    relpath, discovery_reason, rule, outcome}`` with ``outcome`` built by
    :func:`model.outcome_for` (``excluded/archived_path`` or
    ``excluded/preconditioning_run``; a user ``--exclude-glob`` is reported as
    ``archived_path`` with the glob in ``detail``). The exporter writes them to
    ``runs.jsonl`` so an excluded calculation is visible, never silently gone.
    """
    records = []
    for item in ignored:
        if item.get("kind") != "run_dir":
            continue
        reason = item["reason"] if item["reason"] in model.REASONS else "archived_path"
        detail = item.get("detail") or ""
        if item["reason"] != reason:
            detail = f"{item['reason']} ({item.get('rule')}): {detail}"
        records.append({
            "run_id": f"{item['root_alias']}:{item['relpath']}",
            "root_alias": item["root_alias"],
            "relpath": item["relpath"],
            "discovery_reason": item["reason"],
            "rule": item.get("rule"),
            "outcome": model.outcome_for(reason, detail).as_dict(),
        })
    records.sort(key=lambda record: (record["root_alias"], record["relpath"]))
    return records


def discovery_summary(runs: Sequence[DiscoveredRun], ignored: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Deterministic counts for manifests and ``audit.md`` (nothing here is a decision)."""

    def count(values: Iterable[Any]) -> dict[str, int]:
        result: dict[str, int] = {}
        for value in values:
            key = str(value)
            result[key] = result.get(key, 0) + 1
        return dict(sorted(result.items()))

    return {
        "runs": len(runs),
        "runs_with_label_source": sum(1 for run in runs if run.label_file is not None),
        "runs_without_label_source": sum(1 for run in runs if run.label_file is None),
        "runs_by_root": count(run.root_alias for run in runs),
        "runs_by_label_kind": count(run.label_kind or "none" for run in runs),
        "run_flags": count(flag for run in runs for flag in run.flags),
        "ignored": len(ignored),
        "ignored_by_kind": count(item["kind"] for item in ignored),
        "ignored_by_reason": count(item["reason"] for item in ignored),
        "ignored_by_rule": count(item.get("rule") or "-" for item in ignored),
        "pruned_run_candidates_by_reason": count(item["reason"] for item in ignored if item["kind"] == "run_dir"),
    }


# --------------------------------------------------------------------------
# Inventory / selection manifest
# --------------------------------------------------------------------------

INVENTORY_SELECT = ("listed-only", "all")
_ENTRY_STRING_FIELDS = ("lineage", "family", "campaign", "config_type", "notes")
_ENTRY_FIELDS = ("path", "include", "metadata") + _ENTRY_STRING_FIELDS
_TOP_LEVEL_KEYS = ("schema_version", "select", "defaults", "run")


@dataclass(frozen=True)
class InventoryEntry:
    root_alias: str | None  # None: relative to the only root (or unique across roots)
    relpath: str
    include: bool
    lineage: str | None = None
    family: str | None = None
    campaign: str | None = None
    config_type: str | None = None
    notes: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)
    position: int = 0  # 0-based index of the [[run]] table in the file

    @property
    def label(self) -> str:
        return f"{self.root_alias}:{self.relpath}" if self.root_alias else self.relpath

    def as_dict(self) -> dict[str, Any]:
        return {"path": self.label, "include": self.include, "lineage": self.lineage, "family": self.family,
                "campaign": self.campaign, "config_type": self.config_type, "notes": self.notes,
                "metadata": dict(self.metadata), "position": self.position}


@dataclass
class InventoryDecision:
    selected: bool
    reason: str | None  # None | not_in_inventory | inventory_excluded (model.REASONS codes)
    entry: InventoryEntry | None


@dataclass
class InventoryResolution:
    decisions: dict[str, InventoryDecision]  # run_id -> decision
    unmatched: list[dict[str, Any]]  # entries that matched no discovered run (with an explanation)


@dataclass
class Inventory:
    path: Path
    sha256: str
    select: str
    defaults: dict[str, Any]
    entries: list[InventoryEntry]

    def as_dict(self) -> dict[str, Any]:
        return {"path": self.path.as_posix(), "sha256": self.sha256, "select": self.select,
                "defaults": dict(self.defaults), "entries": len(self.entries)}

    def resolve(
        self,
        runs: Sequence[DiscoveredRun],
        *,
        ignored: Sequence[Mapping[str, Any]] = (),
        strict: bool = True,
    ) -> InventoryResolution:
        """Decide for every discovered run whether the inventory selects it.

        ``select="listed-only"``: unlisted runs -> ``not_in_inventory``;
        ``include=false`` entries -> ``inventory_excluded``. Entries without
        an alias match a relpath under any root; matching more than one root is
        an error. Entries matching no run are returned in ``unmatched`` (with
        the discovery reason when the path was pruned); ``strict`` raises
        :class:`DatasetError` for them instead.
        """
        by_key = {(run.root_alias, run.relpath): run for run in runs}
        by_relpath: dict[str, list[DiscoveredRun]] = {}
        for run in runs:
            by_relpath.setdefault(run.relpath, []).append(run)
        assigned: dict[str, InventoryEntry] = {}
        unmatched: list[dict[str, Any]] = []
        for entry in self.entries:
            if entry.root_alias is not None:
                run = by_key.get((entry.root_alias, entry.relpath))
                matches = [run] if run is not None else []
            else:
                matches = by_relpath.get(entry.relpath, [])
            if len(matches) > 1:
                raise DatasetError(
                    f"inventory entry {entry.label!r} matches runs under several roots "
                    f"({[run.root_alias for run in matches]}); write it as alias:relpath"
                )
            if not matches:
                why = [item for item in ignored
                       if item["relpath"] == entry.relpath and entry.root_alias in (None, item["root_alias"])]
                unmatched.append({"entry": entry.as_dict(),
                                  "discovery": [dict(item) for item in why]})
                continue
            run = matches[0]
            if run.run_id in assigned:
                raise DatasetError(f"inventory lists run {run.run_id!r} twice (entries "
                                   f"{assigned[run.run_id].position} and {entry.position})")
            assigned[run.run_id] = entry
        if strict and unmatched:
            details = "; ".join(
                entry["entry"]["path"] + (f" (discovery: {entry['discovery'][0]['reason']})" if entry["discovery"] else "")
                for entry in unmatched[:10]
            )
            raise DatasetError(f"{len(unmatched)} inventory entries match no discovered run: {details}")
        decisions: dict[str, InventoryDecision] = {}
        for run in runs:
            entry = assigned.get(run.run_id)
            if entry is None:
                if self.select == "listed-only":
                    decisions[run.run_id] = InventoryDecision(False, "not_in_inventory", None)
                else:
                    decisions[run.run_id] = InventoryDecision(True, None, None)
            elif not entry.include:
                decisions[run.run_id] = InventoryDecision(False, "inventory_excluded", entry)
            else:
                decisions[run.run_id] = InventoryDecision(True, None, entry)
        return InventoryResolution(decisions, unmatched)


def _normalise_entry_path(raw: Any, where: str) -> tuple[str | None, str]:
    if not isinstance(raw, str) or not raw.strip():
        raise DatasetError(f"{where}: 'path' must be a non-empty string")
    text = raw.strip().replace("\\", "/")
    alias = None
    match = re.match(r"^([A-Za-z0-9][A-Za-z0-9_.-]*):(.*)$", text)
    is_drive = match is not None and len(match.group(1)) == 1 and match.group(2).startswith("/")
    if match is not None and not is_drive:
        alias, text = match.group(1), match.group(2)
    if text.startswith("/") or re.match(r"^[A-Za-z]:/", text):
        raise DatasetError(f"{where}: path {raw!r} must be relative to a root (use alias:relpath)")
    parts = [part for part in text.split("/") if part not in ("", ".")]
    if any(part == ".." for part in parts):
        raise DatasetError(f"{where}: path {raw!r} must not contain '..'")
    return alias, "/".join(parts) if parts else "."


def _check_fields(table: Mapping[str, Any], allowed: Sequence[str], where: str) -> None:
    unknown = sorted(set(table) - set(allowed))
    if unknown:
        raise DatasetError(f"{where}: unknown keys {unknown}; allowed: {sorted(allowed)}")
    if "include" in table and not isinstance(table["include"], bool):
        raise DatasetError(f"{where}: 'include' must be true or false")
    for name in _ENTRY_STRING_FIELDS:
        if name in table and not isinstance(table[name], str):
            raise DatasetError(f"{where}: {name!r} must be a string")
    if "metadata" in table:
        metadata = table["metadata"]
        if not isinstance(metadata, Mapping) or not all(
            isinstance(value, (str, int, float, bool)) for value in metadata.values()
        ):
            raise DatasetError(f"{where}: 'metadata' must be a table of scalar values")


def load_inventory(path: Path) -> Inventory:
    """Load an inventory/selection TOML.

    ::

        schema_version = 1            # optional
        select = "listed-only"        # or "all" (default "listed-only")
        [defaults]                    # optional; any entry field except path
        campaign = "agglo_2026"
        [[run]]
        path = "vasp_runs/n02/r000_s00_1p000/300K"   # relative to a root, or "alias:relpath"
        include = true
        lineage = "agglo-n02-r000"
        family = "n02"
        config_type = "n02_md"
        notes = "..."
        metadata = { temperature_K = 300 }

    Unknown keys are errors (typo protection); entries inherit ``[defaults]``.
    """
    import tomllib

    path = Path(path)
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise DatasetError(f"inventory {path} does not exist") from None
    except (tomllib.TOMLDecodeError, UnicodeDecodeError) as exc:
        raise DatasetError(f"inventory {path} is not valid TOML: {exc}") from None
    _check_top = sorted(set(data) - set(_TOP_LEVEL_KEYS))
    if _check_top:
        raise DatasetError(f"inventory {path.name}: unknown top-level keys {_check_top}; allowed: {list(_TOP_LEVEL_KEYS)}")
    if data.get("schema_version", 1) != 1:
        raise DatasetError(f"inventory {path.name}: unsupported schema_version {data.get('schema_version')!r}")
    select = data.get("select", "listed-only")
    if select not in INVENTORY_SELECT:
        raise DatasetError(f"inventory {path.name}: select must be one of {INVENTORY_SELECT}, got {select!r}")
    defaults = data.get("defaults", {})
    if not isinstance(defaults, Mapping):
        raise DatasetError(f"inventory {path.name}: [defaults] must be a table")
    _check_fields(defaults, [name for name in _ENTRY_FIELDS if name != "path"], f"inventory {path.name} [defaults]")
    rows = data.get("run", [])
    if not isinstance(rows, list):
        raise DatasetError(f"inventory {path.name}: use [[run]] tables")
    entries: list[InventoryEntry] = []
    seen: dict[tuple[str | None, str], int] = {}
    for position, row in enumerate(rows):
        where = f"inventory {path.name} [[run]] #{position + 1}"
        if not isinstance(row, Mapping):
            raise DatasetError(f"{where}: must be a table")
        _check_fields(row, _ENTRY_FIELDS, where)
        if "path" not in row:
            raise DatasetError(f"{where}: 'path' is required")
        alias, relpath = _normalise_entry_path(row["path"], where)
        key = (alias, relpath)
        if key in seen:
            raise DatasetError(f"{where}: duplicate path {row['path']!r} (also entry #{seen[key] + 1})")
        seen[key] = position
        merged = {**defaults, **row}
        metadata = {**dict(defaults.get("metadata", {})), **dict(row.get("metadata", {}))}
        entries.append(InventoryEntry(
            root_alias=alias, relpath=relpath, include=merged.get("include", True),
            lineage=merged.get("lineage"), family=merged.get("family"), campaign=merged.get("campaign"),
            config_type=merged.get("config_type"), notes=merged.get("notes"),
            metadata=metadata, position=position,
        ))
    return Inventory(path.resolve(), sha256_file(path), select, dict(defaults), entries)
