"""Lineage ids and split groups (inventory, agglomeration, InterfaceForge, policy, links).

Directory trees and manifests are SYNTHETIC, built in ``tmp_path`` with the
layouts written by ``nio-md-prep prepare-agglomeration`` and by InterfaceForge
(OPT / Step1 / Step2_<T>K trees joined on ``relative_path``). No VASP output is
needed: lineage only reads manifests, launchers and structure keys.
"""

from __future__ import annotations

import hashlib
import json
import os
import random
from pathlib import Path

import pytest

from nio_md_prep.dataset import lineage as lin
from nio_md_prep.dataset.errors import DatasetError

ALIAS = "loc"


def write_json(path: Path, data) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def run_input(root: Path, rel: str, **kw) -> lin.RunLineageInput:
    (root / rel).mkdir(parents=True, exist_ok=True)
    return lin.RunLineageInput(run_id=f"{ALIAS}:{rel}", root_alias=ALIAS, relpath=rel, root_path=root, **kw)


# --------------------------------------------------------------------------
# union-find
# --------------------------------------------------------------------------

def test_union_find_representative_is_the_smallest_id():
    uf = lin.UnionFind(["d", "b", "c", "a"])
    uf.union("d", "c")
    uf.union("c", "b")
    assert uf.find("d") == "b"
    assert uf.components() == {"a": ["a"], "b": ["b", "c", "d"]}


# --------------------------------------------------------------------------
# inventory
# --------------------------------------------------------------------------

def test_inventory_lineage_has_priority(tmp_path):
    agglo = _agglomeration_tree(tmp_path)
    run = run_input(tmp_path, f"{agglo}/r000_s00_1p000/300K", inventory_lineage="nio-bulk-seed7",
                    inventory_metadata={"family": "bulk"})
    record = lin.lineage_for_run(run)
    assert (record.lineage_id, record.source) == ("inv:nio-bulk-seed7", "inventory")
    assert record.metadata["family"] == "bulk"
    assert lin.lineage_from_inventory("  ") is None and lin.lineage_from_inventory(None) is None


# --------------------------------------------------------------------------
# agglomeration replicas
# --------------------------------------------------------------------------

def _agglomeration_tree(root: Path, *, campaign_manifest=True) -> str:
    campaign = "campaigns/pa_agglo"
    if campaign_manifest:
        write_json(root / campaign / lin.AGGLO_MANIFEST, {"replicas": [{"replica": 0}, {"replica": 1}],
                                                          "config_sha256": "c0ffee"})
    for case, replica, seed in (("r000_s00_1p000", 0, 11), ("r000_s01_1p100", 0, 11), ("r001_s00_1p000", 1, 12)):
        write_json(root / campaign / "vasp_runs" / case / lin.AGGLO_MANIFEST, {
            "agglomerate": "pa_dimer", "replica": replica, "packmol_seed": seed, "center_scale": 1.0,
            "vasp_training_mode": "vasp_md", "composition": [{"slug": "pa", "count": 2}],
            "reference_files": [{"path": "INCAR", "sha256": "AB" * 32}],
        })
    return f"{campaign}/vasp_runs"


def test_agglomeration_replica_is_one_lineage_across_scales_and_stages(tmp_path):
    base = _agglomeration_tree(tmp_path)
    rels = [f"{base}/r000_s00_1p000/300K", f"{base}/r000_s00_1p000/400K/01_heat_300_to_400K",
            f"{base}/r000_s00_1p000/400K/02_hold_400K", f"{base}/r000_s01_1p100/300K",
            f"{base}/r001_s00_1p000/300K"]
    runs = [run_input(tmp_path, rel) for rel in rels]
    result = lin.resolve_lineage(runs)
    r0 = f"agglo:{ALIAS}/campaigns/pa_agglo/pa_dimer/r000"
    r1 = f"agglo:{ALIAS}/campaigns/pa_agglo/pa_dimer/r001"
    assert [result.runs[r.run_id]["lineage_id"] for r in runs] == [r0, r0, r0, r0, r1]
    assert result.groups == {r0: sorted(r.run_id for r in runs[:4]), r1: [runs[4].run_id]}
    record = result.runs[runs[2].run_id]
    assert record["source"] == "agglomeration" and record["lineage_group"] == r0
    meta = record["metadata"]
    assert meta["stage"] == "400K/02_hold_400K" and meta["replica"] == 0 and meta["family"] == "pa_dimer"
    assert meta["composition"] == "pa:2" and meta["config_sha256"] == "c0ffee" and meta["packmol_seed"] == 11
    assert record["evidence"]["campaign_manifest"] == "campaigns/pa_agglo/agglomeration_manifest.json"
    assert result.unresolved == []


def test_agglomeration_reference_hashes_and_partial_copy(tmp_path):
    base = _agglomeration_tree(tmp_path, campaign_manifest=False)
    run_dir = tmp_path / base / "r000_s00_1p000" / "300K"
    run_dir.mkdir(parents=True)
    assert lin.agglomeration_reference_hashes(run_dir, tmp_path) == ["ab" * 32]
    record = lin.lineage_from_agglomeration(run_dir, tmp_path, alias=ALIAS)
    assert record.lineage_id == f"agglo:{ALIAS}/campaigns/pa_agglo/pa_dimer/r000"  # campaign from the layout
    assert record.evidence["campaign_manifest"] is None and "layout" in record.evidence["note"]
    assert lin.lineage_from_agglomeration(tmp_path / "elsewhere", tmp_path, alias=ALIAS) is None


def test_hold_launcher_copying_the_heat_contcar_is_a_continuation(tmp_path):
    base = _agglomeration_tree(tmp_path)
    heat = run_input(tmp_path, f"{base}/r000_s00_1p000/400K/01_heat_300_to_400K")
    hold = run_input(tmp_path, f"{base}/r000_s00_1p000/400K/02_hold_400K")
    (tmp_path / hold.relpath / "runvasp.sh").write_text(
        "#!/bin/bash\nset -e\ncp ../01_heat_300_to_400K/CONTCAR POSCAR\nsrun vasp_std\n", encoding="utf-8")
    targets = lin.launcher_continuations(tmp_path / hold.relpath)
    assert [os.path.normcase(os.path.abspath(p)) for p in targets] == [
        os.path.normcase(os.path.abspath(tmp_path / heat.relpath))]
    result = lin.resolve_lineage([heat, hold])
    assert [link["kind"] for link in result.links] == ["launcher_continuation"]
    assert result.runs[hold.run_id]["links"] == [f"launcher_continuation:{heat.run_id}"]


# --------------------------------------------------------------------------
# InterfaceForge campaign lineage via relative_path across OPT / Step1 / Step2
# --------------------------------------------------------------------------

CAPPED = "OH25/NiO_m110_Big_U46_OH25_clustered_capped"
BARE = "OH0/NiO_m110_Big_U46"


def _interfaceforge_campaign(root: Path, *, poscar_sha=None) -> dict[str, str]:
    camp = root / "if_campaign"
    write_json(camp / "OPT" / "opt_manifest.json", {"format": "interfaceforge-opt-manifest", "runs": [
        {"relative_path": CAPPED, "source": "/abs/cluster/path/OPT/" + CAPPED}, {"relative_path": BARE}]})
    write_json(camp / "OPT" / CAPPED / "provenance.json", {
        "schema": "interfaceforge.reactive-surface/v1", "coverage": 0.25, "arrangement": "clustered",
        "motif": "terminal_hydroxyl", "docking": {"mode": "direct"}, "parent_state": "bare"})
    (camp / "OPT" / BARE).mkdir(parents=True)
    step1_poscar = camp / "Step1" / CAPPED / "POSCAR"
    step1_poscar.parent.mkdir(parents=True)
    step1_poscar.write_text("SYNTHETIC TEST FIXTURE POSCAR\n", encoding="utf-8")
    sha = poscar_sha or hashlib.sha256(step1_poscar.read_bytes()).hexdigest()
    write_json(camp / "Step1" / "step1_manifest.json", {
        "format": "interfaceforge-step1-series", "source_root": "/work/cluster/if_campaign/OPT",
        "protocol": "training", "profile": "nio", "temperature_k": 300,
        "runs": [{"relative_path": CAPPED, "step1_poscar_sha256": sha}]})
    for label, temperature in (("300", 300.0), ("450p5", 450.5)):
        tree = camp / f"Step2_{label}K"
        write_json(tree / "step2_manifest.json", {
            "format": "interfaceforge-step2-series", "source_root": "/work/cluster/if_campaign/Step1",
            "runs": [{"relative_path": CAPPED, "temperature_k": temperature}]})
        (tree / CAPPED / "precondition").mkdir(parents=True)
    return {"opt": f"if_campaign/OPT/{CAPPED}", "opt_bare": f"if_campaign/OPT/{BARE}",
            "step1": f"if_campaign/Step1/{CAPPED}", "step2_300": f"if_campaign/Step2_300K/{CAPPED}",
            "step2_450": f"if_campaign/Step2_450p5K/{CAPPED}",
            "precondition": f"if_campaign/Step2_300K/{CAPPED}/precondition"}


def test_interfaceforge_relative_path_joins_opt_step1_and_step2(tmp_path):
    paths = _interfaceforge_campaign(tmp_path)
    runs = {name: run_input(tmp_path, rel) for name, rel in paths.items()}
    result = lin.resolve_lineage(runs.values())
    capped = f"iface:{ALIAS}/if_campaign/{CAPPED}"
    bare = f"iface:{ALIAS}/if_campaign/{BARE}"
    for name in ("opt", "step1", "step2_300", "step2_450", "precondition"):
        assert result.runs[runs[name].run_id]["lineage_id"] == capped, name
        assert result.runs[runs[name].run_id]["source"] == "interfaceforge"
    assert result.runs[runs["opt_bare"].run_id]["lineage_id"] == bare
    assert sorted(result.groups) == sorted([capped, bare]) and len(result.groups[capped]) == 5
    step2 = result.runs[runs["step2_450"].run_id]
    assert step2["metadata"]["campaign_stage"] == "step2" and step2["metadata"]["temperature_k"] == 450.5
    assert step2["metadata"]["fam_coverage"] == 0.25 and step2["metadata"]["fam_motif"] == "terminal_oh"
    assert step2["metadata"]["fam_anchor"] == "direct" and step2["metadata"]["family_source"] == "provenance"
    assert step2["evidence"]["opt_tree"] == "if_campaign/OPT"
    assert any("sibling" in note for note in step2["evidence"]["notes"])  # absolute source_root not present
    step1 = result.runs[runs["step1"].run_id]
    assert step1["metadata"]["temperature_k"] == 300 and step1["metadata"]["if_protocol"] == "training"
    assert step1["evidence"]["sha256_checks"] == {"step1_poscar_sha256": True}
    assert step1["evidence"]["lineage_verified"] is True
    nested = result.runs[runs["precondition"].run_id]
    assert nested["evidence"]["relation"] == "nested_under_listed_run"
    bare_record = result.runs[runs["opt_bare"].run_id]
    assert bare_record["metadata"]["family_source"] == "name" and bare_record["metadata"]["structure_family"] == "OH0"


def test_interfaceforge_sha_edge_mismatch_is_recorded(tmp_path):
    paths = _interfaceforge_campaign(tmp_path, poscar_sha="0" * 64)
    record, why = lin.lineage_from_interfaceforge(tmp_path / paths["step1"], tmp_path, alias=ALIAS)
    assert why is None and record.evidence["lineage_verified"] is False


def test_unlisted_directory_in_an_interfaceforge_tree_fails_closed(tmp_path):
    _interfaceforge_campaign(tmp_path)
    stray = run_input(tmp_path, "if_campaign/Step1/OH50/not_in_manifest")
    result = lin.resolve_lineage([stray])
    record = result.runs[stray.run_id]
    assert record["source"] == "unresolved" and record["lineage_id"] is None
    assert any("not listed in step1_manifest.json" in note for note in record["evidence"]["notes"])
    assert result.unresolved == [stray.run_id] and record["group_id"] == f"unresolved:{stray.run_id}"
    # a declared policy resolves it explicitly (and records that it did)
    record = lin.resolve_lineage([stray], policy="parent").runs[stray.run_id]
    assert record["source"] == "policy:parent" and record["lineage_id"] == f"dir:{ALIAS}:if_campaign/Step1/OH50"


def test_manifest_errors_are_recorded_not_fatal(tmp_path):
    bad = tmp_path / "broken" / "step1_manifest.json"
    bad.parent.mkdir(parents=True)
    bad.write_text("{not json", encoding="utf-8")
    run = run_input(tmp_path, "broken/case")
    result = lin.resolve_lineage([run])
    assert result.unresolved == [run.run_id]
    assert list(result.manifest_errors) == [bad.as_posix()]


# --------------------------------------------------------------------------
# declared policy
# --------------------------------------------------------------------------

@pytest.mark.parametrize("text, parsed", [(None, None), ("run", ("run", None)), ("parent", ("parent", None)),
                                          (" depth:2 ", ("depth", 2))])
def test_parse_lineage_policy(text, parsed):
    assert lin.parse_lineage_policy(text) == parsed


@pytest.mark.parametrize("text", ["dir", "depth:0", "depth:-1", "depth:x", ""])
def test_bad_lineage_policy_is_an_error(text):
    with pytest.raises(DatasetError, match="lineage policy"):
        lin.parse_lineage_policy(text)


def test_policy_lineage_ids():
    run = lin.RunLineageInput(run_id="loc:a/b/c", root_alias="loc", relpath="a/b/c")
    assert lin.lineage_from_policy(run, "run").lineage_id == "run:loc:a/b/c"
    assert lin.lineage_from_policy(run, "parent").lineage_id == "dir:loc:a/b"
    assert lin.lineage_from_policy(run, "depth:1").lineage_id == "dir:loc:a"
    assert lin.lineage_from_policy(run, "depth:1").source == "policy:depth:1"
    top = lin.RunLineageInput(run_id="loc:.", root_alias="loc", relpath=".")
    assert lin.lineage_from_policy(top, "parent").lineage_id == "dir:loc:."
    with pytest.raises(DatasetError):
        lin.resolve_lineage([run], policy="sideways")


# --------------------------------------------------------------------------
# continuation links, duplicates, unresolved fail-closed
# --------------------------------------------------------------------------

def _keyed(run_id, lineage=None, **kw):
    return lin.RunLineageInput(run_id=run_id, root_alias="loc", relpath=run_id.split(":", 1)[1],
                               inventory_lineage=lineage, **kw)


def test_structure_key_links_merge_restarts_and_common_starts():
    runs = [
        _keyed("loc:md1", "traj-1", initial_key="s0", frame_keys=["s0", "s1", "s2"], final_key="s2f"),
        _keyed("loc:md1_restart", "traj-2", initial_key="s2f", frame_keys=["s2f", "s3"]),  # CONTCAR restart
        _keyed("loc:md1_rewind", "traj-3", initial_key="s1", frame_keys=["s1", "s4"]),  # rewind to a frame
        _keyed("loc:other_T", "traj-4", initial_key="s0", frame_keys=["s0", "s9"]),  # same start, other T
        _keyed("loc:independent", "traj-5", initial_key="z0", frame_keys=["z0"]),
    ]
    result = lin.resolve_lineage(runs)
    kinds = {(link["a"], link["b"]): link["kind"] for link in result.links}
    assert kinds[("loc:md1", "loc:md1_restart")] == "continuation"
    assert kinds[("loc:md1", "loc:md1_rewind")] == "restart_from_frame"
    assert kinds[("loc:md1", "loc:other_T")] == "same_initial_structure"
    assert result.groups == {"inv:traj-1": ["loc:md1", "loc:md1_restart", "loc:md1_rewind", "loc:other_T"],
                             "inv:traj-5": ["loc:independent"]}
    # each run keeps its own lineage id; the split group is the component's smallest resolved id
    assert result.runs["loc:md1_restart"]["lineage_id"] == "inv:traj-2"
    assert result.group_of("loc:md1_restart") == "inv:traj-1"


def test_nested_runs_and_duplicate_links_join_groups():
    runs = [_keyed("loc:a", "A"), _keyed("loc:a/sp", "B", nested_in="a"), _keyed("loc:c", "C"), _keyed("loc:d", "D")]
    extra = [{"a": "loc:c", "b": "loc:d", "kind": "exact_duplicate"},
             {"run_a": "loc:c", "run_b": "loc:gone", "kind": "near_duplicate"}]
    result = lin.resolve_lineage(runs, extra_links=extra)
    assert result.groups == {"inv:A": ["loc:a", "loc:a/sp"], "inv:C": ["loc:c", "loc:d"]}
    ignored = [link for link in result.links if link.get("ignored")]
    assert ignored == [{"a": "loc:c", "b": "loc:gone", "kind": "near_duplicate", "ignored": True}]
    assert {link["kind"] for link in result.links if not link.get("ignored")} == {"nested_run", "exact_duplicate"}


def test_unresolved_runs_stay_unresolved_even_when_linked():
    runs = [_keyed("loc:known", "K", final_key="f"), _keyed("loc:mystery", initial_key="f"),
            _keyed("loc:alone")]
    result = lin.resolve_lineage(runs)
    assert result.unresolved == ["loc:alone", "loc:mystery"]
    assert result.runs["loc:mystery"]["source"] == "unresolved"
    assert result.group_of("loc:mystery") == "inv:K"  # grouped with its parent, still refused by the splitter
    assert result.group_of("loc:alone") == "unresolved:loc:alone"
    counts = result.as_dict()["counts"]
    assert counts["by_source"] == {"inventory": 1, "unresolved": 2} and counts["groups"] == 2


def test_resolution_is_independent_of_input_order():
    runs = [_keyed(f"loc:r{i}", f"L{i % 3}", initial_key=f"k{i % 4}", final_key=f"f{i}") for i in range(12)]
    first = lin.resolve_lineage(runs).as_dict()
    shuffled = runs[:]
    random.Random(7).shuffle(shuffled)
    assert lin.resolve_lineage(shuffled).as_dict() == first


def test_duplicate_run_ids_are_an_error():
    with pytest.raises(DatasetError, match="duplicate run_id"):
        lin.resolve_lineage([_keyed("loc:a", "A"), _keyed("loc:a", "B")])


def test_mapping_inputs_are_accepted():
    result = lin.resolve_lineage([{"run_id": "loc:a", "root_alias": "loc", "relpath": "a", "inventory_lineage": "A",
                                   "metadata": {"temperature_K": 300.0}}])
    assert result.runs["loc:a"]["metadata"] == {"temperature_K": 300.0}
    assert result.runs["loc:a"]["lineage_group"] == "inv:A"
