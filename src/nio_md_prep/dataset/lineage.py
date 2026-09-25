"""Lineage ids and split groups for discovered runs.

Every run gets one lineage record ``{lineage_id, source, evidence, metadata}``
from the first source that applies (priority order):

1. **inventory** -- an explicit ``lineage`` in the inventory TOML entry;
2. **agglomeration** -- the nio-md-prep agglomeration layout: walk up to the
   case ``agglomeration_manifest.json`` (has ``agglomerate`` and ``replica``, no
   ``replicas``) and the campaign root manifest (has ``replicas``);
   ``lineage_id = agglo:<alias>/<campaign relpath>/<agglomerate>/r<replica:03d>``
   (one Packmol replica = one lineage: all its scale variants, the 300 K run,
   the heating and hold continuations and every single point);
3. **interfaceforge** -- InterfaceForge campaign trees recognised by
   ``opt_manifest.json`` / ``step1_manifest.json`` / ``step2_manifest.json``
   (``format`` = ``interfaceforge-opt-manifest`` / ``-step1-series`` /
   ``-step2-series``). The run must be listed in its tree manifest; the join key
   is the manifest ``relative_path`` (shared by OPT, Step1 and every
   ``Step2_<T>K`` tree), never a directory-name coincidence and never the
   absolute paths stored in the manifests:
   ``lineage_id = iface:<alias>/<campaign parent relpath>/<relative_path>``;
4. **declared policy** -- ``run`` | ``parent`` | ``depth:N`` (explicit CLI choice);
5. otherwise ``source="unresolved"`` (the splitter refuses such runs).

Links (union-find over runs) then merge correlated runs into one split group:
same lineage id; nested run directories; continuation/restart chains (a run
whose initial structure equals another run's final structure or any of its
frames, or a launcher that copies another run's CONTCAR); runs that start from
the same structure; exact duplicates and shared structures across runs; and
near-duplicates when enabled. ``group_id`` is the smallest resolved lineage id
of the component (``unresolved:<run_id>`` only when nothing in it is resolved).
A run whose own lineage is unresolved stays ``unresolved`` even when linked:
fail closed.

The InterfaceForge tree layout and manifest formats recognised here are those
of lgutsev/InterfaceForge@4501e34 (src/interfaceforge/vasp.py step1/step2
prepare, step1_repair.py; MIT), as mapped in design/understand/if-phase1-map.md
sections 2 and 4.4. Only the conventions are used; no InterfaceForge code is
copied.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any, Iterable, Mapping, Sequence

from .errors import DatasetError
from .fsio import sha256_file

LINEAGE_SCHEMA = 1
SOURCES = ("inventory", "agglomeration", "interfaceforge", "policy", "unresolved")

AGGLO_MANIFEST = "agglomeration_manifest.json"
IF_MANIFESTS = {
    "opt_manifest.json": ("interfaceforge-opt-manifest", "opt"),
    "step1_manifest.json": ("interfaceforge-step1-series", "step1"),
    "step2_manifest.json": ("interfaceforge-step2-series", "step2"),
}
IF_REPAIR = "step1_repair.json"
IF_REPAIR_FORMAT = "interfaceforge-step1-repair"
_STEP2_TREE = re.compile(r"^Step2_(\d+(?:p\d+)?)K$")
_OH_PART = re.compile(r"^OH(\d+)$")
_LAUNCHERS = ("runvasp.sh", "run.slurm", "submit.sh", "job.sh")
_CP_CONTCAR = re.compile(r"(?:^|[;&|\s])cp\s+(?:-[a-zA-Z]+\s+)*['\"]?([^\s'\";&|]+?)/CONTCAR['\"]?\s+['\"]?(?:\./)?POSCAR['\"]?(?:\s|$|;)")
_POLICY = re.compile(r"^(run|parent|depth:(\d+))$")


# --------------------------------------------------------------------------
# Union-find
# --------------------------------------------------------------------------

class UnionFind:
    """Deterministic union-find over string ids (the representative is the smallest id)."""

    def __init__(self, items: Iterable[str] = ()):
        self._parent: dict[str, str] = {}
        for item in items:
            self.add(item)

    def add(self, item: str) -> None:
        self._parent.setdefault(item, item)

    def find(self, item: str) -> str:
        self.add(item)
        root = item
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[item] != root:  # path compression
            self._parent[item], item = root, self._parent[item]
        return root

    def union(self, a: str, b: str) -> str:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        low, high = sorted((ra, rb))
        self._parent[high] = low
        return low

    def components(self) -> dict[str, list[str]]:
        groups: dict[str, list[str]] = {}
        for item in sorted(self._parent):
            groups.setdefault(self.find(item), []).append(item)
        return groups


# --------------------------------------------------------------------------
# Records and inputs
# --------------------------------------------------------------------------

@dataclass
class LineageRecord:
    lineage_id: str | None
    source: str  # inventory | agglomeration | interfaceforge | policy:<name> | unresolved
    evidence: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {"lineage_id": self.lineage_id, "source": self.source, "evidence": dict(self.evidence),
                "metadata": dict(self.metadata)}


@dataclass
class RunLineageInput:
    """What lineage resolution needs to know about one run."""

    run_id: str
    root_alias: str
    relpath: str  # posix, relative to the root ("." for the root itself)
    root_path: Path | None = None
    run_dir: Path | None = None
    nested_in: str | None = None  # relpath of the enclosing run directory
    inventory_lineage: str | None = None
    inventory_metadata: Mapping[str, Any] = field(default_factory=dict)
    initial_key: str | None = None  # order key of <structure name="initialpos">
    frame_keys: Sequence[str] = ()  # order keys of every step's structure (DFT and MLFF)
    final_key: str | None = None  # order key of <structure name="finalpos">
    metadata: Mapping[str, Any] = field(default_factory=dict)  # e.g. temperature from TEBEG/TEEND

    @property
    def directory(self) -> Path | None:
        if self.run_dir is not None:
            return Path(self.run_dir)
        if self.root_path is not None:
            return Path(self.root_path) if self.relpath == "." else Path(self.root_path) / self.relpath
        return None


class ManifestCache:
    """Reads each JSON manifest once; unreadable files are recorded, never fatal."""

    def __init__(self):
        self._data: dict[str, Any] = {}
        self.errors: dict[str, str] = {}

    def load(self, path: Path) -> dict[str, Any] | None:
        key = os.path.normcase(os.path.abspath(path))
        if key not in self._data:
            value = None
            if path.is_file():
                try:
                    value = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
                    self.errors[Path(path).as_posix()] = f"{type(exc).__name__}: {exc}"
                    value = None
                if value is not None and not isinstance(value, dict):
                    self.errors[Path(path).as_posix()] = "manifest is not a JSON object"
                    value = None
            self._data[key] = value
        return self._data[key]


def _posix_rel(path: Path, root: Path) -> str:
    rel = Path(os.path.relpath(path, root)).as_posix()
    return "." if rel in ("", ".") else rel


def _join_id(prefix: str, *parts: str) -> str:
    kept = [p.strip("/") for p in parts if p not in (None, "", ".")]
    return prefix + "/".join(kept)


def _ancestors(start: Path, root: Path) -> list[Path]:
    """``start`` and its parents up to and including ``root`` (empty if start is outside root)."""
    start, root = Path(os.path.abspath(start)), Path(os.path.abspath(root))
    try:
        start.relative_to(root)
    except ValueError:
        return []
    chain = [start]
    while chain[-1] != root and chain[-1].parent != chain[-1]:
        chain.append(chain[-1].parent)
    return chain


# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------

def lineage_from_inventory(lineage: str | None, *, metadata: Mapping[str, Any] | None = None) -> LineageRecord | None:
    """An explicit inventory ``lineage`` string (``inv:<value>``)."""
    if lineage is None or not str(lineage).strip():
        return None
    value = str(lineage).strip()
    return LineageRecord(f"inv:{value}", "inventory", {"inventory_lineage": value}, dict(metadata or {}))


def _composition(value: Any) -> str | None:
    if isinstance(value, list):
        parts = []
        for item in value:
            if isinstance(item, Mapping) and "slug" in item:
                parts.append(f"{item['slug']}:{item.get('count')}")
        return ",".join(sorted(parts)) if parts else None
    return None


def _reference_hashes(*manifests: Mapping[str, Any] | None) -> list[str]:
    hashes: set[str] = set()

    def take(entries):
        for item in entries or []:
            if isinstance(item, Mapping) and isinstance(item.get("sha256"), str):
                hashes.add(item["sha256"].lower())

    for manifest in manifests:
        if not manifest:
            continue
        take(manifest.get("reference_files"))
        for job in manifest.get("vasp_jobs") or []:
            if isinstance(job, Mapping):
                take(job.get("reference_files"))
        sanity = manifest.get("vasp_reference_sanity")
        if isinstance(sanity, Mapping):
            take(sanity.get("shared_files"))
    return sorted(hashes)


def lineage_from_agglomeration(run_dir: Path, root: Path, *, alias: str,
                               cache: ManifestCache | None = None) -> LineageRecord | None:
    """Agglomeration replica lineage (see the module docstring), or None outside such a tree."""
    cache = cache or ManifestCache()
    case_dir = case = campaign_dir = campaign = None
    for directory in _ancestors(run_dir, root):
        data = cache.load(directory / AGGLO_MANIFEST)
        if data is None:
            continue
        if "replicas" in data:
            campaign_dir, campaign = directory, data
            break
        if case is None and "agglomerate" in data and "replica" in data:
            case_dir, case = directory, data
    if case is None:
        return None
    try:
        replica = int(case["replica"])
    except (TypeError, ValueError):
        return None
    agglomerate = str(case["agglomerate"])
    evidence: dict[str, Any] = {"case_manifest": _posix_rel(case_dir / AGGLO_MANIFEST, root)}
    if campaign_dir is not None:
        campaign_rel = _posix_rel(campaign_dir, root)
        evidence["campaign_manifest"] = _posix_rel(campaign_dir / AGGLO_MANIFEST, root)
    else:
        # partial copy without the root manifest: the directory above vasp_runs/structures, else the case parent
        parts = PurePosixPath(_posix_rel(case_dir, root)).parts
        anchor = next((i for i, part in enumerate(parts) if part in ("vasp_runs", "structures")), None)
        campaign_rel = "/".join(parts[:anchor]) if anchor is not None else "/".join(parts[:-1])
        campaign_rel = campaign_rel or "."
        evidence["campaign_manifest"] = None
        evidence["note"] = "campaign root manifest not found; campaign identity taken from the directory layout"
    lineage_id = _join_id(f"agglo:{alias}/", campaign_rel, agglomerate, f"r{replica:03d}")
    stage = _posix_rel(run_dir, case_dir)
    metadata = {
        "family": agglomerate,
        "campaign": campaign_rel,
        "structure_family": agglomerate,
        "replica": replica,
        "case": case_dir.name,
        "stage": stage,
        "composition": _composition(case.get("composition")),
        "packmol_seed": case.get("packmol_seed"),
        "center_scale": case.get("center_scale"),
        "vasp_training_mode": case.get("vasp_training_mode"),
        "config_sha256": (campaign or {}).get("config_sha256"),
    }
    evidence["reference_hashes"] = _reference_hashes(case, campaign)
    return LineageRecord(lineage_id, "agglomeration", evidence, {k: v for k, v in metadata.items() if v is not None})


def agglomeration_reference_hashes(run_dir: Path, root: Path, cache: ManifestCache | None = None) -> list[str]:
    """sha256 of files copied from reference/template dirs into this agglomeration run (for acceptance)."""
    record = lineage_from_agglomeration(run_dir, root, alias="_", cache=cache)
    return list(record.evidence.get("reference_hashes", [])) if record else []


def _if_tree(directory: Path, cache: ManifestCache) -> tuple[str, dict[str, Any], str] | None:
    """``(stage, manifest, filename)`` if ``directory`` is an InterfaceForge tree root."""
    for filename, (fmt, stage) in IF_MANIFESTS.items():
        data = cache.load(directory / filename)
        if data is not None and data.get("format") == fmt and isinstance(data.get("runs"), list):
            return stage, data, filename
    return None


def _entry_path(entry: Mapping[str, Any]) -> str | None:
    value = entry.get("relative_path")
    if not isinstance(value, str) or not value.strip():
        return None
    parts = [p for p in value.replace("\\", "/").split("/") if p not in ("", ".")]
    if any(p == ".." for p in parts):
        return None
    return "/".join(parts) or "."


def _temperature_label(text: str) -> float | None:
    try:
        return float(text.replace("p", "."))
    except ValueError:
        return None


def _normalise_family(provenance: Mapping[str, Any]) -> dict[str, Any]:
    """Both OPT provenance schemas -> fam_coverage/arrangement/motif/ligand/anchor/parent."""
    motif_map = {"capped": "terminal_oh", "terminal_hydroxyl": "terminal_oh", "terminal_oh": "terminal_oh",
                 "dissoc": "dissociated", "dissociated_water": "dissociated", "dissociated": "dissociated",
                 "bare": "bare"}
    family: dict[str, Any] = {}
    coverage = provenance.get("coverage_fraction", provenance.get("coverage"))
    if isinstance(coverage, (int, float)) and not isinstance(coverage, bool):
        family["fam_coverage"] = float(coverage)
    arrangement = provenance.get("arrangement", provenance.get("pattern_id"))
    if arrangement is not None:
        family["fam_arrangement"] = str(arrangement)
    motif = provenance.get("motif", provenance.get("hydroxylation_motif"))
    if motif is not None:
        family["fam_motif"] = motif_map.get(str(motif).lower(), str(motif))
    ligand = provenance.get("ligand")
    if ligand is not None:
        family["fam_ligand"] = str(ligand)
    docking = provenance.get("docking")
    anchor = provenance.get("anchor_position")
    if anchor is None and isinstance(docking, Mapping):
        anchor = docking.get("mode")
    if anchor is not None:
        family["fam_anchor"] = str(anchor)
    parent = provenance.get("parent_state", provenance.get("base_reference"))
    if parent is not None:
        family["fam_parent"] = str(parent)
    return family


def _family_from_name(relative_path: str) -> dict[str, Any]:
    family: dict[str, Any] = {}
    parts = PurePosixPath(relative_path).parts
    for part in parts:
        match = _OH_PART.match(part)
        if match:
            family["fam_coverage"] = int(match.group(1)) / 100.0
    leaf = parts[-1] if parts else ""
    for token, name in (("clustered", "fam_arrangement"), ("scattered", "fam_arrangement"),
                        ("full", "fam_arrangement")):
        if f"_{token}" in leaf:
            family[name] = token
    for token, motif in (("capped", "terminal_oh"), ("terminal_hydroxyl", "terminal_oh"),
                         ("dissoc", "dissociated"), ("dissociated_water", "dissociated")):
        if f"_{token}" in leaf:
            family["fam_motif"] = motif
            break
    for token in ("boundary", "bare", "hbond", "direct"):
        if leaf.endswith(f"_{token}"):
            family["fam_anchor"] = token
    return family


def _campaign_opt_root(tree_root: Path, stage: str, manifest: Mapping[str, Any], root: Path,
                       cache: ManifestCache, relative_path: str, depth: int = 0) -> tuple[Path | None, list[str]]:
    """The OPT tree this Step1/Step2 tree descends from (via ``source_root``), or None."""
    notes: list[str] = []
    if stage == "opt":
        return tree_root, notes
    if depth > 3:
        return None, ["source_root chain too deep"]
    source = manifest.get("source_root")
    candidate = None
    if isinstance(source, str) and source.strip():
        path = Path(source)
        if path.is_absolute() and path.is_dir() and _if_tree(path, cache) is not None:
            candidate = path
        else:
            sibling = tree_root.parent / PurePosixPath(source.replace("\\", "/")).name
            if sibling.is_dir() and _if_tree(sibling, cache) is not None:
                candidate = sibling
                notes.append("source_root resolved by name among sibling trees (absolute path not present)")
    if candidate is None:
        # fall back to a unique sibling OPT tree that lists this relative_path
        siblings = []
        for child in sorted(tree_root.parent.iterdir()) if tree_root.parent.is_dir() else []:
            tree = _if_tree(child, cache) if child.is_dir() else None
            if tree is not None and tree[0] == "opt" and any(
                    _entry_path(e) == relative_path for e in tree[1]["runs"] if isinstance(e, Mapping)):
                siblings.append(child)
        if len(siblings) == 1:
            notes.append("OPT tree found as the unique sibling tree listing this relative_path")
            return siblings[0], notes
        notes.append("source tree not resolved" if not siblings else "several sibling OPT trees list this relative_path")
        return None, notes
    tree = _if_tree(candidate, cache)
    if tree is None:
        return None, notes + ["source_root is not an InterfaceForge tree"]
    found, more = _campaign_opt_root(candidate, tree[0], tree[1], root, cache, relative_path, depth + 1)
    return found, notes + more


def lineage_from_interfaceforge(run_dir: Path, root: Path, *, alias: str,
                                cache: ManifestCache | None = None) -> tuple[LineageRecord | None, str | None]:
    """InterfaceForge campaign lineage, or ``(None, why)`` (``why`` None when not inside an IF tree)."""
    cache = cache or ManifestCache()
    run_dir = Path(os.path.abspath(run_dir))
    for directory in _ancestors(run_dir, root)[1:]:
        tree = _if_tree(directory, cache)
        if tree is None:
            continue
        stage, manifest, filename = tree
        rel = _posix_rel(run_dir, directory)
        entries = [e for e in manifest["runs"] if isinstance(e, Mapping)]
        exact = [e for e in entries if _entry_path(e) == rel]
        nested = [e for e in entries if _entry_path(e) not in (None, ".") and rel.startswith(_entry_path(e) + "/")]
        if exact:
            entry, relation = exact[0], "listed"
        elif nested:
            entry, relation = max(nested, key=lambda e: len(_entry_path(e) or "")), "nested_under_listed_run"
        else:
            return None, f"inside InterfaceForge {stage} tree {_posix_rel(directory, root)} but not listed in {filename}"
        relative_path = _entry_path(entry)
        opt_root, notes = _campaign_opt_root(directory, stage, manifest, root, cache, relative_path)
        campaign_dir = (opt_root or directory).parent
        try:
            Path(os.path.abspath(campaign_dir)).relative_to(Path(os.path.abspath(root)))
            campaign_rel = _posix_rel(campaign_dir, root)
        except ValueError:
            campaign_rel = _posix_rel(directory.parent, root)
            notes.append("OPT tree outside the scanned root; campaign taken from this tree")
        lineage_id = _join_id(f"iface:{alias}/", campaign_rel, relative_path)
        evidence: dict[str, Any] = {
            "manifest": _posix_rel(directory / filename, root), "relation": relation,
            "campaign_stage": stage, "tree": _posix_rel(directory, root),
            "opt_tree": _posix_rel(opt_root, root) if opt_root is not None else None, "notes": notes,
        }
        metadata: dict[str, Any] = {"campaign": campaign_rel, "campaign_stage": stage, "structure_id": relative_path}
        if relation != "listed":
            metadata["nested_path"] = rel
        match = _STEP2_TREE.match(directory.name)
        if stage == "step2":
            temperature = entry.get("temperature_k")
            if temperature is None and match:
                temperature = _temperature_label(match.group(1))
            if temperature is not None:
                metadata["temperature_k"] = float(temperature)
        elif stage == "step1" and manifest.get("temperature_k") is not None:
            metadata["temperature_k"] = manifest.get("temperature_k")
        for key in ("protocol", "profile"):
            if isinstance(manifest.get(key), str):
                metadata[f"if_{key}"] = manifest[key]
        family: dict[str, Any] = {}
        family_source = None
        if opt_root is not None:
            provenance = cache.load(opt_root / relative_path / "provenance.json")
            if provenance is not None:
                family, family_source = _normalise_family(provenance), "provenance"
        if not family:
            family, family_source = _family_from_name(relative_path), "name"
        metadata.update(family)
        metadata["family_source"] = family_source
        metadata["structure_family"] = PurePosixPath(relative_path).parts[0] if relative_path != "." else "."
        # sha edges recorded in the manifest (proof that the tree edge is real)
        verified = {}
        for sha_key, filename_on_disk in (("step1_poscar_sha256", "POSCAR"), ("step2_poscar_sha256", "POSCAR")):
            expected = entry.get(sha_key)
            path = run_dir / filename_on_disk
            if isinstance(expected, str) and path.is_file():
                verified[sha_key] = sha256_file(path).lower() == expected.lower()
        if verified:
            evidence["sha256_checks"] = verified
            evidence["lineage_verified"] = all(verified.values())
        repair = cache.load(run_dir / IF_REPAIR)
        if repair is not None and repair.get("format") == IF_REPAIR_FORMAT:
            archive_root = run_dir / ".interfaceforge" / "archive"
            archives = sorted(p.name for p in archive_root.glob("step1_repair_*")) if archive_root.is_dir() else []
            metadata["repair"] = {
                "status": repair.get("status"), "safe_prefix_steps": repair.get("safe_prefix_steps"),
                "rewind_frame": repair.get("rewind_frame"), "source": repair.get("source"),
                "first_bad_step": (repair.get("diagnostic") or {}).get("first_bad_step"),
                "scf_unreliable": (repair.get("diagnostic") or {}).get("scf_unreliable"),
                "archives": archives, "segment": len(archives),
                "limitation": "archived repair prefixes are excluded (--include-repair-prefix not implemented)",
            }
        return LineageRecord(lineage_id, "interfaceforge", evidence, metadata), None
    return None, None


def parse_lineage_policy(text: str | None) -> tuple[str, int | None] | None:
    """``run`` | ``parent`` | ``depth:N`` -> ``(name, N)``; None -> None."""
    if text is None:
        return None
    match = _POLICY.match(str(text).strip())
    if not match:
        raise DatasetError(f"lineage policy must be 'run', 'parent' or 'depth:N', got {text!r}")
    if match.group(2) is not None:
        depth = int(match.group(2))
        if depth < 1:
            raise DatasetError("lineage policy depth:N needs N >= 1")
        return "depth", depth
    return match.group(1), None


def lineage_from_policy(run: RunLineageInput, policy: str | None) -> LineageRecord | None:
    parsed = parse_lineage_policy(policy)
    if parsed is None:
        return None
    name, depth = parsed
    parts = list(PurePosixPath(run.relpath).parts) if run.relpath != "." else []
    if name == "run":
        lineage_id = f"run:{run.run_id}"
    elif name == "parent":
        lineage_id = f"dir:{run.root_alias}:" + ("/".join(parts[:-1]) or ".")
    else:
        lineage_id = f"dir:{run.root_alias}:" + ("/".join(parts[:depth]) or ".")
    label = name if depth is None else f"depth:{depth}"
    return LineageRecord(lineage_id, f"policy:{label}", {"policy": label})


def launcher_continuations(run_dir: Path) -> list[Path]:
    """Directories whose CONTCAR a launcher in ``run_dir`` copies to POSCAR (restart/continuation)."""
    found: list[Path] = []
    for name in _LAUNCHERS:
        path = Path(run_dir) / name
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            stripped = line.split("#", 1)[0]
            for match in _CP_CONTCAR.finditer(" " + stripped):
                target = Path(os.path.normpath(Path(run_dir) / match.group(1)))
                if target not in found:
                    found.append(target)
    return found


# --------------------------------------------------------------------------
# Resolution
# --------------------------------------------------------------------------

@dataclass
class LineageResolution:
    runs: dict[str, dict[str, Any]]  # run_id -> {lineage_id, source, evidence, metadata, group_id, links}
    groups: dict[str, list[str]]  # group_id -> run_ids
    links: list[dict[str, Any]]
    manifest_errors: dict[str, str] = field(default_factory=dict)
    policy: str | None = None

    def group_of(self, run_id: str) -> str:
        return self.runs[run_id]["group_id"]

    @property
    def unresolved(self) -> list[str]:
        return sorted(run_id for run_id, rec in self.runs.items() if rec["source"] == "unresolved")

    def as_dict(self) -> dict[str, Any]:
        sources: dict[str, int] = {}
        for record in self.runs.values():
            sources[record["source"]] = sources.get(record["source"], 0) + 1
        return {
            "schema": LINEAGE_SCHEMA, "policy": self.policy,
            "runs": {k: self.runs[k] for k in sorted(self.runs)},
            "groups": {k: self.groups[k] for k in sorted(self.groups)},
            "links": self.links, "unresolved": self.unresolved,
            "counts": {"runs": len(self.runs), "groups": len(self.groups), "links": len(self.links),
                       "by_source": dict(sorted(sources.items()))},
            "manifest_errors": dict(sorted(self.manifest_errors.items())),
        }


def _input(run: Any) -> RunLineageInput:
    if isinstance(run, RunLineageInput):
        return run
    get = (lambda key, default=None: run.get(key, default)) if isinstance(run, Mapping) else \
        (lambda key, default=None: getattr(run, key, default))
    return RunLineageInput(
        run_id=get("run_id"), root_alias=get("root_alias"), relpath=get("relpath"),
        root_path=get("root_path"), run_dir=get("run_dir"), nested_in=get("nested_in"),
        inventory_lineage=get("inventory_lineage"), inventory_metadata=get("inventory_metadata") or {},
        initial_key=get("initial_key"), frame_keys=get("frame_keys") or (), final_key=get("final_key"),
        metadata=get("metadata") or {},
    )


def lineage_for_run(run: RunLineageInput, *, policy: str | None = None,
                    cache: ManifestCache | None = None) -> LineageRecord:
    """The run's own lineage record, by source priority (links are applied by :func:`resolve_lineage`)."""
    cache = cache or ManifestCache()
    notes: list[str] = []
    record = lineage_from_inventory(run.inventory_lineage, metadata=run.inventory_metadata)
    directory = run.directory
    root = Path(run.root_path) if run.root_path is not None else None
    if record is None and directory is not None and root is not None:
        record = lineage_from_agglomeration(directory, root, alias=run.root_alias, cache=cache)
        if record is None:
            record, why = lineage_from_interfaceforge(directory, root, alias=run.root_alias, cache=cache)
            if why:
                notes.append(why)
    if record is None:
        record = lineage_from_policy(run, policy)
    if record is None:
        record = LineageRecord(None, "unresolved", {"why": "no inventory lineage, no agglomeration or "
                                                           "InterfaceForge manifest, no declared --lineage-policy"})
    if notes:
        record.evidence.setdefault("notes", [])
        record.evidence["notes"] = list(record.evidence["notes"]) + notes
    record.metadata = {**dict(run.metadata), **record.metadata}
    return record


def resolve_lineage(
    runs: Iterable[Any],
    *,
    policy: str | None = None,
    extra_links: Iterable[Mapping[str, Any]] = (),
    cache: ManifestCache | None = None,
) -> LineageResolution:
    """Lineage record + split group for every run (independent of input order).

    ``extra_links``: ``{a, b, kind}`` run-id pairs from the duplicate finder
    (exact duplicates, shared structures, near-duplicates). Links to unknown
    run ids are ignored (recorded with ``ignored=True``).
    """
    parse_lineage_policy(policy)  # validate early
    cache = cache or ManifestCache()
    inputs = sorted((_input(r) for r in runs), key=lambda r: r.run_id)
    by_id: dict[str, RunLineageInput] = {}
    for run in inputs:
        if run.run_id in by_id:
            raise DatasetError(f"duplicate run_id {run.run_id!r} in lineage resolution")
        by_id[run.run_id] = run
    records = {run.run_id: lineage_for_run(run, policy=policy, cache=cache) for run in inputs}

    uf = UnionFind(by_id)
    links: dict[tuple[str, str, str], dict[str, Any]] = {}

    def link(a: str, b: str, kind: str, **extra: Any) -> None:
        if a == b:
            return
        x, y = sorted((a, b))
        key = (x, y, kind)
        if key not in links:
            links[key] = {"a": x, "b": y, "kind": kind, **extra}
        uf.union(x, y)

    # same lineage id
    by_lineage: dict[str, list[str]] = {}
    for run_id, record in records.items():
        if record.lineage_id is not None:
            by_lineage.setdefault(record.lineage_id, []).append(run_id)
    for lineage_id, members in sorted(by_lineage.items()):
        for other in members[1:]:
            uf.union(members[0], other)
    # nested run directories
    for run in inputs:
        if run.nested_in is not None:
            parent = f"{run.root_alias}:{run.nested_in}"
            if parent in by_id:
                link(run.run_id, parent, "nested_run")
    # continuation / restart / common start by structure keys
    final_index: dict[str, list[str]] = {}
    frame_index: dict[str, list[str]] = {}
    initial_index: dict[str, list[str]] = {}
    for run in inputs:
        if run.final_key:
            final_index.setdefault(run.final_key, []).append(run.run_id)
        for key in set(run.frame_keys or ()):
            frame_index.setdefault(key, []).append(run.run_id)
        if run.initial_key:
            initial_index.setdefault(run.initial_key, []).append(run.run_id)
    for run in inputs:
        key = run.initial_key
        if not key:
            continue
        for other in final_index.get(key, []):
            if other != run.run_id:
                link(run.run_id, other, "continuation", detail=f"{run.run_id} starts from the final structure of {other}")
        for other in frame_index.get(key, []):
            if other != run.run_id:
                link(run.run_id, other, "restart_from_frame", detail=f"{run.run_id} starts from a frame of {other}")
        for other in initial_index.get(key, []):
            if other != run.run_id:
                link(run.run_id, other, "same_initial_structure")
    # launcher continuations (cp ../other/CONTCAR POSCAR)
    by_dir = {os.path.normcase(os.path.abspath(r.directory)): r.run_id for r in inputs if r.directory is not None}
    for run in inputs:
        if run.directory is None:
            continue
        for target in launcher_continuations(run.directory):
            other = by_dir.get(os.path.normcase(os.path.abspath(target)))
            if other is not None and other != run.run_id:
                link(run.run_id, other, "launcher_continuation")
    # duplicate finder links
    ignored: list[dict[str, Any]] = []
    for item in extra_links:
        a, b, kind = item.get("a") or item.get("run_a"), item.get("b") or item.get("run_b"), item.get("kind", "link")
        if a in by_id and b in by_id:
            link(a, b, str(kind))
        else:
            ignored.append({"a": a, "b": b, "kind": kind, "ignored": True})

    components = uf.components()
    groups: dict[str, list[str]] = {}
    group_of: dict[str, str] = {}
    for members in components.values():
        resolved = sorted(records[m].lineage_id for m in members
                          if records[m].source != "unresolved" and records[m].lineage_id is not None)
        group_id = resolved[0] if resolved else "unresolved:" + min(members)
        groups[group_id] = sorted(members)
        for member in members:
            group_of[member] = group_id
    link_list = [links[k] for k in sorted(links)]
    touching: dict[str, list[str]] = {run_id: [] for run_id in records}
    for item in link_list:
        touching[item["a"]].append(f"{item['kind']}:{item['b']}")
        touching[item["b"]].append(f"{item['kind']}:{item['a']}")
    per_run: dict[str, dict[str, Any]] = {}
    for run_id in sorted(records):
        record = records[run_id]
        per_run[run_id] = {**record.as_dict(), "group_id": group_of[run_id], "lineage_group": group_of[run_id],
                           "links": sorted(touching[run_id])}
    return LineageResolution(runs=per_run, groups=groups, links=link_list + ignored,
                             manifest_errors=dict(cache.errors), policy=policy)
