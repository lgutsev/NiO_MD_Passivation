"""Tests for dataset/discovery.py.

Every directory tree here is SYNTHETIC: files are tiny placeholders whose
content is the line ``SYNTHETIC TEST FIXTURE`` (discovery only looks at file
names, never parses them). Layouts imitate the agglomeration campaign tree
(main-reuse-map section 2), the InterfaceForge OPT/Step1/Step2 trees with
``precondition/`` and ``.interfaceforge/archive`` (if-phase1-map section 2)
and OutPackLite packages (OSZICAR/XDATCAR/INCAR, no OUTCAR/vasprun).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import random

import pytest

from nio_md_prep.dataset import discovery as D
from nio_md_prep.dataset import model
from nio_md_prep.dataset.errors import DatasetError

PLACEHOLDER = "SYNTHETIC TEST FIXTURE\n"


def _touch(root: Path, *relpaths: str) -> None:
    for relpath in relpaths:
        path = root / relpath
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(PLACEHOLDER, encoding="utf-8")


def _agglomeration_tree(root: Path) -> None:
    _touch(
        root,
        "agglomeration_manifest.json",
        "templates/mol.xyz",
        "packmol/n02/replica_000/packed.xyz",
        "xtb/n02/r000_s00_1p000/OUTCAR",  # an OUTCAR-named file under xtb/: never a VASP run
        "structures/n02/r000/POSCAR",
        "vasp_reference/INCAR", "vasp_reference/OUTCAR",
        "vasp_runs/n02/r000_s00_1p000/agglomeration_manifest.json",
        "vasp_runs/n02/r000_s00_1p000/POSCAR",
        "vasp_runs/n02/r000_s00_1p000/300K/INCAR",
        "vasp_runs/n02/r000_s00_1p000/300K/vasprun.xml",
        "vasp_runs/n02/r000_s00_1p000/300K/OUTCAR",
        "vasp_runs/n02/r000_s00_1p000/300K/OSZICAR",
        "vasp_runs/n02/r000_s00_1p000/300K/vasprun.xml.1",  # rerun leftover: not a label file
        "vasp_runs/n02/r000_s00_1p000/300K/OUTCAR.bak",
        "vasp_runs/n02/r000_s00_1p000/400K/01_heat_300_to_400K/vasprun.xml.gz",
        "vasp_runs/n02/r000_s00_1p000/400K/01_heat_300_to_400K/vasprun.xml",  # plain preferred over .gz
        "vasp_runs/n02/r000_s00_1p000/400K/02_hold_400K/vasprun.xml.xz",
        "vasp_runs/n02/r000_s00_1p000/400K/02_hold_400K/OUTCAR.gz",
        "vasp_runs/n02/r001_s00_1p000/300K/INCAR",  # prepared, never run
        "vasp_runs/n02/r001_s00_1p000/300K/POSCAR",
        "vasp_runs/n02/r002_s00_1p000/300K/OUTCAR",  # OUTCAR only
        "vasp_runs/n02/old_r003/300K/vasprun.xml",
        "vasp_runs/n02/r004_archived/300K/vasprun.xml",
        "vasp_runs/n02/backup_r005/300K/vasprun.xml",
        "vasp_runs/n02/r006.bak/300K/vasprun.xml",
        "vasp_runs/n02/packed.xyz.geometry_0123456789ab/OUTCAR",
    )


def _interfaceforge_tree(root: Path) -> None:
    rel = "OH25/NiO_m110_Big_U46_OH25_clustered_capped"
    _touch(
        root,
        "OPT/opt_manifest.json",
        f"OPT/{rel}/INCAR", f"OPT/{rel}/vasprun.xml", f"OPT/{rel}/OUTCAR", f"OPT/{rel}/provenance.json",
        "Step1/step1_manifest.json",
        f"Step1/{rel}/INCAR", f"Step1/{rel}/vasprun.xml", f"Step1/{rel}/OUTCAR", f"Step1/{rel}/step1_repair.json",
        f"Step1/{rel}/precondition/INCAR", f"Step1/{rel}/precondition/vasprun.xml", f"Step1/{rel}/precondition/OUTCAR",
        f"Step1/{rel}/.interfaceforge/archive/step1_repair_20260101T000000Z/vasprun.xml",
        f"Step1/{rel}/.interfaceforge/archive/step1_repair_20260101T000000Z/OUTCAR",
        "Step2_300K/step2_manifest.json",
        f"Step2_300K/{rel}/INCAR", f"Step2_300K/{rel}/vasprun.xml", f"Step2_300K/{rel}/OUTCAR",
        f"Step2_300K/{rel}/precondition/POSCAR",  # the Step2 artefact: inputs only, VASP never ran
        f"Step2_300K/{rel}/precondition/KPOINTS",
        "Step2_450K/restart_archive_1/vasprun.xml",
        "Step2_450K/refit_archive_2/OUTCAR",
        "Step2_450K/Backup_old/OUTCAR",
        "X_OutPackLite_MDNiO_m100_331_MD_Jul22_2024/NiO_331_Med_U46/OSZICAR",
        "X_OutPackLite_MDNiO_m100_331_MD_Jul22_2024/NiO_331_Med_U46/XDATCAR_FINAL",
        "X_OutPackLite_MDNiO_m100_331_MD_Jul22_2024/NiO_331_Med_U46/INCAR",
        "OutPackLite_unprefixed/NiO_331_ORich_U46/OSZICAR",  # same package without the X prefix
        "OutPackLite_unprefixed/NiO_331_ORich_U46/XDATCAR",
        "OutPackLite_unprefixed/NiO_331_ORich_U46/INCAR",
    )


def _all_candidates(root: Path) -> set[str]:
    """Every directory holding a run-defining output (ground truth for 'nothing silently dropped')."""
    defining = {name for name, _ in D.LABEL_FILES}
    for kind in D.RUN_DEFINING_EVIDENCE:
        defining |= set(D.EVIDENCE_FILES[kind])
    found = set()
    for dirpath, dirnames, filenames in os.walk(root):
        if defining & set(filenames):
            found.add(Path(dirpath).relative_to(root).as_posix() or ".")
    return found


def _by_rel(items):
    return {item["relpath"]: item for item in items}


# --------------------------------------------------------------------------
# agglomeration tree
# --------------------------------------------------------------------------

def test_agglomeration_tree_runs_labels_and_reported_exclusions(tmp_path):
    root = tmp_path / "agglo"
    _agglomeration_tree(root)
    runs, ignored = D.discover_runs(D.parse_roots([f"agg={root}"]))
    by_rel = {run.relpath: run for run in runs}
    assert sorted(by_rel) == [
        "vasp_runs/n02/r000_s00_1p000/300K",
        "vasp_runs/n02/r000_s00_1p000/400K/01_heat_300_to_400K",
        "vasp_runs/n02/r000_s00_1p000/400K/02_hold_400K",
        "vasp_runs/n02/r002_s00_1p000/300K",
    ]
    assert by_rel["vasp_runs/n02/r000_s00_1p000/300K"].label_kind == "vasprun"
    assert by_rel["vasp_runs/n02/r000_s00_1p000/400K/01_heat_300_to_400K"].label_kind == "vasprun"
    hold = by_rel["vasp_runs/n02/r000_s00_1p000/400K/02_hold_400K"]
    assert hold.label_kind == "vasprun.xz" and hold.evidence["OUTCAR"].name == "OUTCAR.gz"
    assert hold.run_id == "agg:vasp_runs/n02/r000_s00_1p000/400K/02_hold_400K"

    outcar_only = by_rel["vasp_runs/n02/r002_s00_1p000/300K"]
    assert outcar_only.label_file is None and "no_label_source" in outcar_only.flags
    outcome = outcar_only.missing_label_outcome()
    assert outcome.status == model.MISSING_LABELS and outcome.reason == "no_vasprun" and "OUTCAR" in outcome.detail
    assert by_rel["vasp_runs/n02/r000_s00_1p000/300K"].missing_label_outcome() is None

    ig = {(item["relpath"], item["reason"]) for item in ignored}
    for part in ("templates", "packmol", "xtb", "structures", "vasp_reference", "vasp_runs/n02/old_r003",
                 "vasp_runs/n02/r004_archived", "vasp_runs/n02/backup_r005", "vasp_runs/n02/r006.bak",
                 "vasp_runs/n02/packed.xyz.geometry_0123456789ab"):
        assert (part, "archived_path") in ig, part
    assert ("vasp_runs/n02/r000_s00_1p000/300K/vasprun.xml.1", "not_a_label_file") in ig
    assert ("vasp_runs/n02/r000_s00_1p000/300K/OUTCAR.bak", "not_a_label_file") in ig
    assert ("vasp_runs/n02/r000_s00_1p000/400K/01_heat_300_to_400K/vasprun.xml.gz", "duplicate_label_file") in ig
    assert ("vasp_runs/n02/r001_s00_1p000/300K", "no_vasp_outputs") in ig
    # run candidates inside excluded trees are reported one by one
    pruned = _by_rel(D.pruned_run_candidates(ignored))
    assert {"xtb/n02/r000_s00_1p000", "vasp_reference", "vasp_runs/n02/old_r003/300K"} <= set(pruned)
    assert pruned["vasp_runs/n02/old_r003/300K"]["outcome"] == {
        "status": "excluded", "reason": "archived_path",
        "detail": pruned["vasp_runs/n02/old_r003/300K"]["outcome"]["detail"]}
    assert _all_candidates(root) == set(by_rel) | set(pruned)  # nothing silently dropped


# --------------------------------------------------------------------------
# InterfaceForge + OutPackLite
# --------------------------------------------------------------------------

def test_interfaceforge_parts_and_outpacklite_are_excluded_and_reported(tmp_path):
    root = tmp_path / "iface"
    _interfaceforge_tree(root)
    rel = "OH25/NiO_m110_Big_U46_OH25_clustered_capped"
    runs, ignored = D.discover_runs(D.parse_roots([root]))
    by_rel = {run.relpath: run for run in runs}
    assert set(by_rel) == {f"OPT/{rel}", f"Step1/{rel}", f"Step2_300K/{rel}", "OutPackLite_unprefixed/NiO_331_ORich_U46"}
    assert by_rel[f"Step1/{rel}"].metadata_files.keys() == {"step1_repair.json"}
    assert by_rel[f"OPT/{rel}"].metadata_files.keys() == {"provenance.json"}

    package = by_rel["OutPackLite_unprefixed/NiO_331_ORich_U46"]
    assert package.label_file is None and {"no_label_source", "outpacklite_like"} <= set(package.flags)
    outcome = package.missing_label_outcome()
    assert outcome.status == "missing_labels" and "never used as labels" in outcome.detail

    pruned = _by_rel(D.pruned_run_candidates(ignored))
    assert pruned[f"Step1/{rel}/precondition"]["outcome"]["reason"] == "preconditioning_run"
    assert pruned[f"Step1/{rel}/precondition"]["outcome"]["status"] == "excluded"
    archive = f"Step1/{rel}/.interfaceforge/archive/step1_repair_20260101T000000Z"
    assert pruned[archive]["outcome"]["reason"] == "archived_path" and pruned[archive]["rule"] == "part:.interfaceforge"
    for part, rule in [("Step2_450K/restart_archive_1", "part:*archive*"), ("Step2_450K/refit_archive_2", "part:*archive*"),
                       ("Step2_450K/Backup_old", "part:*backup*"),
                       ("X_OutPackLite_MDNiO_m100_331_MD_Jul22_2024/NiO_331_Med_U46", "part:X*:cs")]:
        assert pruned[part]["rule"] == rule, part
    ig = {(item["relpath"], item["kind"]): item for item in ignored}
    step2_pre = ig[(f"Step2_300K/{rel}/precondition", "dir")]  # the empty Step2 artefact is reported, not a run
    assert step2_pre["reason"] == "preconditioning_run"
    assert _all_candidates(root) == set(by_rel) | set(pruned)

    summary = D.discovery_summary(runs, ignored)
    assert summary["runs"] == 4 and summary["runs_without_label_source"] == 1
    assert summary["run_flags"]["outpacklite_like"] == 1
    assert summary["pruned_run_candidates_by_reason"] == {"archived_path": 5, "preconditioning_run": 1}


def test_deliberately_included_parts_are_scanned_and_tagged(tmp_path):
    root = tmp_path / "iface"
    _interfaceforge_tree(root)
    rel = "OH25/NiO_m110_Big_U46_OH25_clustered_capped"
    runs, ignored = D.discover_runs(D.parse_roots([root]), include_parts=["precondition", "X*"])
    by_rel = {run.relpath: run for run in runs}
    pre = by_rel[f"Step1/{rel}/precondition"]
    assert "included_excluded_part:precondition" in pre.flags
    assert "nested_run" in pre.flags and pre.nested_in == f"Step1/{rel}"  # linked to the MD run it preconditions
    lite = by_rel["X_OutPackLite_MDNiO_m100_331_MD_Jul22_2024/NiO_331_Med_U46"]
    assert "included_excluded_part:X*" in lite.flags and lite.missing_label_outcome().status == "missing_labels"
    # the .interfaceforge archive stays excluded; the artefact precondition/ dir is now a plain no-output dir
    pruned = _by_rel(D.pruned_run_candidates(ignored))
    assert f"Step1/{rel}/.interfaceforge/archive/step1_repair_20260101T000000Z" in pruned
    assert f"Step2_300K/{rel}/precondition" not in {item["relpath"] for item in ignored}
    manifest = D.exclusion_rules_manifest(include_parts=["precondition", "X*"])
    assert manifest["include_parts"] == ["X*", "precondition"]
    flagged = {rule["pattern"] for rule in manifest["default_rules"] if rule["deliberately_included"]}
    assert flagged == {"X*", "precondition"}
    with pytest.raises(DatasetError, match="names no default exclusion rule"):
        D.discover_runs(D.parse_roots([root]), include_parts=["preconditon"])  # typo is an error


def test_including_interfaceforge_dir_still_excludes_archives_inside(tmp_path):
    root = tmp_path / "iface"
    _interfaceforge_tree(root)
    rel = "OH25/NiO_m110_Big_U46_OH25_clustered_capped"
    runs, ignored = D.discover_runs(D.parse_roots([root]), include_parts=[".interfaceforge"])
    pruned = _by_rel(D.pruned_run_candidates(ignored))
    archive = f"Step1/{rel}/.interfaceforge/archive/step1_repair_20260101T000000Z"
    assert pruned[archive]["rule"] == "part:*archive*"
    runs, ignored = D.discover_runs(D.parse_roots([root]), include_parts=[".interfaceforge", "*archive*"])
    run = {run.relpath: run for run in runs}[archive]
    assert {"included_excluded_part:.interfaceforge", "included_excluded_part:*archive*", "nested_run"} <= set(run.flags)


# --------------------------------------------------------------------------
# user globs
# --------------------------------------------------------------------------

def test_user_globs_include_exclude_precedence(tmp_path):
    root = tmp_path / "agglo"
    _agglomeration_tree(root)
    runs, ignored = D.discover_runs(D.parse_roots([f"agg={root}"]), include_globs=["vasp_runs/n02/old_r003"],
                                    exclude_globs=["vasp_runs/n02/r000_s00_1p000/400K"])
    by_rel = {run.relpath: run for run in runs}
    assert "included_by_glob:vasp_runs/n02/old_r003" in by_rel["vasp_runs/n02/old_r003/300K"].flags
    assert not any(rel.startswith("vasp_runs/n02/r000_s00_1p000/400K") for rel in by_rel)
    pruned = _by_rel(D.pruned_run_candidates(ignored))
    heat = pruned["vasp_runs/n02/r000_s00_1p000/400K/01_heat_300_to_400K"]
    assert heat["discovery_reason"] == "user_excluded" and heat["outcome"]["reason"] == "archived_path"
    assert "exclude-glob:vasp_runs/n02/r000_s00_1p000/400K" in heat["outcome"]["detail"]
    # exclude wins over include
    runs, _ = D.discover_runs(D.parse_roots([f"agg={root}"]), include_globs=["vasp_runs/n02/old_r003"],
                              exclude_globs=["vasp_runs/n02/old_r003"])
    assert "vasp_runs/n02/old_r003/300K" not in {run.relpath for run in runs}
    # alias-qualified globs
    runs, _ = D.discover_runs(D.parse_roots([f"agg={root}"]), exclude_globs=["agg:vasp_runs/n02/r002*"])
    assert "vasp_runs/n02/r002_s00_1p000/300K" not in {run.relpath for run in runs}


# --------------------------------------------------------------------------
# nesting, output dir, symlinks, determinism
# --------------------------------------------------------------------------

def test_nested_runs_are_flagged_and_output_dir_is_skipped(tmp_path):
    root = tmp_path / "root"
    _touch(root, "parent/vasprun.xml", "parent/child/vasprun.xml", "out/dataset.extxyz", "out/sub/vasprun.xml")
    runs, ignored = D.discover_runs(D.parse_roots([root]), output_dir=root / "out")
    by_rel = {run.relpath: run for run in runs}
    assert set(by_rel) == {"parent", "parent/child"}
    assert by_rel["parent/child"].nested_in == "parent" and "nested_run" in by_rel["parent/child"].flags
    assert {"relpath": "out", "reason": "output_dir"}.items() <= _by_rel(ignored)["out"].items()
    with pytest.raises(DatasetError, match="inside the output directory"):
        D.discover_runs(D.parse_roots([root / "parent"]), output_dir=root)


def test_symlink_loops_are_not_followed(tmp_path):
    root = tmp_path / "root"
    _touch(root, "a/vasprun.xml")
    try:
        os.symlink(root, root / "a" / "loop", target_is_directory=True)
    except (OSError, NotImplementedError):
        pytest.skip("directory symlinks not permitted on this system")
    runs, ignored = D.discover_runs(D.parse_roots([root]))
    assert [run.relpath for run in runs] == ["a"]
    assert _by_rel(ignored)["a/loop"]["reason"] == "symlink_not_followed"
    runs, ignored = D.discover_runs(D.parse_roots([root]), follow_symlinks=True)
    assert [run.relpath for run in runs] == ["a"]
    assert _by_rel(ignored)["a/loop"]["reason"] == "symlink_loop"


def test_discovery_is_independent_of_root_order_and_creation_order(tmp_path):
    def build(base: Path, order_seed: int):
        files = []
        for name in ("agglo", "iface"):
            scratch = tmp_path / f"list-{name}"
            (_agglomeration_tree if name == "agglo" else _interfaceforge_tree)(scratch)
            files += [(name, path.relative_to(scratch).as_posix()) for path in sorted(scratch.rglob("*")) if path.is_file()]
        random.Random(order_seed).shuffle(files)
        for name, relpath in files:
            _touch(base / name, relpath)
        return base

    one = build(tmp_path / "one", 1)
    two = build(tmp_path / "two", 2)

    def snapshot(base, reverse):
        roots = [f"agg={base / 'agglo'}", f"if={base / 'iface'}"]
        runs, ignored = D.discover_runs(D.parse_roots(roots[::-1] if reverse else roots))
        return json.dumps([[run.as_dict() for run in runs], ignored, D.pruned_run_candidates(ignored)],
                          sort_keys=True, default=str)

    assert snapshot(one, False) == snapshot(one, True) == snapshot(two, False)


# --------------------------------------------------------------------------
# roots
# --------------------------------------------------------------------------

def test_parse_roots_aliases_and_overlap(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b" / "inner").mkdir(parents=True)
    specs = D.parse_roots([f"first={tmp_path / 'a'}", tmp_path / "b"])
    assert [spec.alias for spec in specs] == ["first", "b"] and specs[1].path == (tmp_path / "b").resolve()
    assert specs[0].as_dict()["path"] == (tmp_path / "a").resolve().as_posix()
    with pytest.raises(DatasetError, match="duplicate root aliases"):
        D.parse_roots([f"x={tmp_path / 'a'}", f"x={tmp_path / 'b'}"])
    with pytest.raises(DatasetError, match="overlap"):
        D.parse_roots([tmp_path / "b", tmp_path / "b" / "inner"])
    with pytest.raises(DatasetError, match="same directory"):
        D.parse_roots([f"p={tmp_path / 'a'}", f"q={tmp_path / 'a'}"])
    with pytest.raises(DatasetError, match="does not exist"):
        D.parse_roots([tmp_path / "missing"])
    with pytest.raises(DatasetError, match="at least one"):
        D.parse_roots([])
    (tmp_path / "bad name").mkdir()
    with pytest.raises(DatasetError, match="alias"):
        D.parse_roots([tmp_path / "bad name"])


# --------------------------------------------------------------------------
# inventory
# --------------------------------------------------------------------------

def test_inventory_selects_excludes_and_reports_unmatched(tmp_path):
    root = tmp_path / "agglo"
    _agglomeration_tree(root)
    runs, ignored = D.discover_runs(D.parse_roots([f"agg={root}"]))
    inventory_path = tmp_path / "inventory.toml"
    inventory_path.write_text(
        "# SYNTHETIC TEST FIXTURE\n"
        "select = \"listed-only\"\n"
        "[defaults]\ncampaign = \"agglo_synthetic\"\n"
        "[[run]]\npath = \"vasp_runs/n02/r000_s00_1p000/300K\"\nlineage = \"agglo-n02-r000\"\nfamily = \"n02\"\n"
        "metadata = { temperature_K = 300 }\n"
        "[[run]]\npath = \"agg:vasp_runs/n02/r002_s00_1p000/300K\"\ninclude = false\nnotes = \"OUTCAR only\"\n"
        "[[run]]\npath = \"vasp_runs\\\\n02\\\\old_r003\\\\300K\"\n",
        encoding="utf-8",
    )
    inventory = D.load_inventory(inventory_path)
    assert inventory.select == "listed-only" and len(inventory.entries) == 3
    assert inventory.entries[2].relpath == "vasp_runs/n02/old_r003/300K"  # backslashes normalised
    with pytest.raises(DatasetError, match="match no discovered run.*archived_path"):
        inventory.resolve(runs, ignored=ignored)
    resolution = inventory.resolve(runs, ignored=ignored, strict=False)
    decisions = resolution.decisions
    selected = decisions["agg:vasp_runs/n02/r000_s00_1p000/300K"]
    assert selected.selected and selected.entry.campaign == "agglo_synthetic" and selected.entry.metadata == {"temperature_K": 300}
    assert decisions["agg:vasp_runs/n02/r002_s00_1p000/300K"].reason == "inventory_excluded"
    assert decisions["agg:vasp_runs/n02/r000_s00_1p000/400K/02_hold_400K"].reason == "not_in_inventory"
    assert all(decision.reason in (None, *model.REASONS) for decision in decisions.values())
    assert resolution.unmatched[0]["discovery"][0]["reason"] == "archived_path"


@pytest.mark.parametrize("text, message", [
    ("bogus = 1\n", "unknown top-level keys"),
    ("[[run]]\npath = \"a\"\nlinage = \"x\"\n", "unknown keys"),
    ("[[run]]\npath = \"../a\"\n", "must not contain"),
    ("[[run]]\npath = \"/abs/a\"\n", "relative"),
    ("[[run]]\npath = \"a\"\n[[run]]\npath = \"a/\"\n", "duplicate path"),
    ("[[run]]\npath = \"a\"\ninclude = \"yes\"\n", "true or false"),
    ("select = \"some\"\n", "select must be"),
    ("schema_version = 2\n", "schema_version"),
    ("[[run]]\ninclude = true\n", "'path' is required"),
    ("not toml ===", "not valid TOML"),
])
def test_inventory_validation_errors(tmp_path, text, message):
    path = tmp_path / "inventory.toml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(DatasetError, match=message):
        D.load_inventory(path)
