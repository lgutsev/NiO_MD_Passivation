"""End-to-end acceptance test: ``nio-md-prep dataset scan|export|split|audit`` on a SYNTHETIC tree.

Every VASP file here is SYNTHETIC (``tests/dataset_fixtures.py``; labelled in-file). The tree
mimics the layouts the exporter must handle:

* an agglomeration campaign with campaign/case manifests (three Packmol replicas; replica 0 has a
  300 K run, a heating and a hold continuation whose launcher copies the heating CONTCAR), a stale
  ``*.protocol_*`` artefact and a ``vasp_reference`` copy;
* generic directories covering every required failure mode (static/relaxation/MD, truncated MD,
  POTCAR species order, unconverged SCF, missing stress, exact and contradictory duplicates, two
  reference pools, compressed vasprun, restart chain, energy-only, VASP-MLFF steps, selective
  dynamics, magnetic transition, uncontrolled magnetism, OUTCAR-only, unreadable vasprun,
  archived/backup reruns, a non-canonical ``vasprun.xml.1``);
* an InterfaceForge-like campaign (OPT / Step1 / Step2_300K / Step2_450p5K joined on
  ``relative_path``) with a ``precondition/`` run and a ``.interfaceforge/archive`` rewind copy;
* an OutPackLite-like ``X*`` package (OSZICAR/XDATCAR/INCAR/CONTCAR, no labels).

Like real NiO data, every run is spin-polarized (ISPIN is a hard reference-settings field, so
mixing ISPIN=1 and ISPIN=2 runs would split the pool): runs carry an explicit AFM MAGMOM and
AFM final site moments unless a test case sets its own magnetism.

Nothing here validates real NiO forces.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

import numpy as np
import pytest

from dataset_fixtures import make_frames, rocksalt_nio, write_vasp_run

from nio_md_prep import cli as top_cli
from nio_md_prep.dataset import audit as au
from nio_md_prep.dataset import export as ex
from nio_md_prep.dataset import split as sp
from nio_md_prep.dataset.acceptance import Policy
from nio_md_prep.dataset.errors import DatasetError, LeakageError, OverwriteRefusedError
from nio_md_prep.dataset.extxyz import read_frame_ids
from nio_md_prep.dataset.fsio import read_jsonl, sha256_file

pytest.importorskip("ase")
import ase.io  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
MD = {"IBRION": 0, "POTIM": 1.0}
AFM_MAGMOM = "2 -2 2 -2 0 0 0 0"
AFM_SITES = [1.7, -1.7, 1.7, -1.7, 0.0, 0.0, 0.0, 0.0]
FM_SITES = [1.7, 1.7, 1.7, 1.7, 0.0, 0.0, 0.0, 0.0]
CAPPED = "OH25/NiO_m110_Big_U46_OH25_clustered_capped"
BARE = "OH0/NiO_m110_Big_U46"
AGGLO = "agglo/pa_campaign"
CASES = {0: "r000_s00_1p000", 1: "r001_s00_1p000", 2: "r002_s00_1p000"}
SYNTHETIC_TEXT = "SYNTHETIC TEST FIXTURE - not real VASP output\n"


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1, sort_keys=True), encoding="utf-8")


def structure(k: int):
    """A distinct Ni4O4 cell per run (no accidental cross-run structure links)."""
    return rocksalt_nio(a=4.10 + 0.013 * k)


def spin_incar(incar=None) -> dict:
    """Default NiO magnetism: ISPIN=2 with an explicit AFM MAGMOM (unless the case sets ISPIN)."""
    incar = dict(incar or {})
    if "ISPIN" not in incar:
        incar.update(ISPIN=2, MAGMOM=AFM_MAGMOM)
    return incar


def frames_for(k: int, n: int, **fields):
    species, cell, positions = structure(k)
    fields.setdefault("mag_total", 0.0)
    fields.setdefault("site_moments", AFM_SITES)
    return species, cell, make_frames(species, cell, positions, n, seed=1000 + k, base_energy=-50.0 - k, **fields)


# --------------------------------------------------------------------------
# The synthetic tree
# --------------------------------------------------------------------------

def _run_specs() -> list[tuple[str, str, dict]]:
    """(root, relpath, write_vasp_run kwargs) for every SYNTHETIC run (pure data, no I/O)."""
    specs: list[tuple[str, str, dict]] = []

    def add(root, rel, k, n, *, incar=None, frames=None, cell=None, species=None, **kw):
        if frames is None:
            species, cell, frames = frames_for(k, n)
        specs.append((root, rel, dict(species=species, cell=cell, frames=frames, incar=spin_incar(incar), **kw)))

    # agglomeration campaign (lineage from manifests)
    base = f"{AGGLO}/vasp_runs"
    add("calcs", f"{base}/{CASES[0]}/300K", 1, 4, incar=dict(MD, NSW=4, TEBEG=300))
    add("calcs", f"{base}/{CASES[0]}/400K/01_heat_300_to_400K", 2, 3, incar=dict(MD, NSW=3, TEBEG=300, TEEND=400))
    add("calcs", f"{base}/{CASES[0]}/400K/02_hold_400K", 3, 3, incar=dict(MD, NSW=3, TEBEG=400))
    add("calcs", f"{base}/{CASES[1]}/300K", 4, 4, incar=dict(MD, NSW=4, TEBEG=300))
    add("calcs", f"{base}/{CASES[2]}/300K", 5, 4, incar=dict(MD, NSW=4, TEBEG=300))
    add("calcs", f"{base}/{CASES[0]}.protocol_old/300K", 6, 2, incar=dict(MD, NSW=2))  # stale artefact
    add("calcs", f"{AGGLO}/vasp_reference", 7, 1)  # copied template, never a run
    # generic directories (need a declared lineage policy)
    add("calcs", "generic/static_a", 10, 1)
    add("calcs", "generic/relax_b", 11, 4, incar={"IBRION": 2, "NSW": 4})
    add("calcs", "generic/relax_b/backup_1", 12, 2, incar={"IBRION": 2, "NSW": 2})
    add("calcs", "generic/old/rerun", 13, 2, incar={"IBRION": 2, "NSW": 2})
    add("calcs", "generic/md_truncated", 14, 5, incar=dict(MD, NSW=5), truncate_after_steps=3,
        truncate_inside_step=True)
    add("calcs", "generic/species_mismatch", 15, 3, incar=dict(MD, NSW=3), potcar_file_order=["O", "Ni"])
    species, cell, frames = frames_for(16, 4)
    frames[1]["scf_marker"] = "not_reached"
    frames[2].update(n_scf=60, last_dE=1e-3, scf_marker=None)
    add("calcs", "generic/scf_fail", 16, 4, incar=dict(MD, NSW=4), species=species, cell=cell, frames=frames)
    add("calcs", "generic/isif0_md", 17, 4, incar=dict(MD, NSW=4, ISIF=0))
    species, cell, frames = frames_for(18, 3)
    add("calcs", "generic/dup_a", 18, 3, incar=dict(MD, NSW=3), species=species, cell=cell, frames=frames)
    add("calcs", "generic/dup_b", 18, 1, species=species, cell=cell, frames=[dict(frames[2])])
    species, cell, frames = frames_for(19, 3)
    add("calcs", "generic/conflict_a", 19, 3, incar=dict(MD, NSW=3), species=species, cell=cell, frames=frames)
    add("calcs", "generic/conflict_b", 19, 1, species=species, cell=cell,
        frames=[dict(frames[1], free_energy=frames[1]["free_energy"] + 0.5)])
    add("calcs", "generic/compressed_md", 20, 3, incar=dict(MD, NSW=3), compress="gz")
    species, cell, frames = frames_for(21, 3)
    add("calcs", "generic/restart_1", 21, 3, incar=dict(MD, NSW=3), species=species, cell=cell, frames=frames)
    _, _, later = frames_for(22, 2)
    later = [dict(f, positions=f["positions"] - structure(22)[2] + structure(21)[2]) for f in later]
    add("calcs", "generic/restart_2", 21, 3, incar=dict(MD, NSW=3), species=species, cell=cell,
        frames=[dict(frames[2])] + later)
    add("calcs", "generic/energy_only", 23, 3, incar=dict(MD, NSW=3),
        frames=frames_for(23, 3, omit={"forces"})[2], species=structure(23)[0], cell=structure(23)[1])
    add("calcs", "generic/mlff_md", 24, 5, incar=dict(MD, NSW=5, ML_LMLFF=True, ML_MODE="train"), mlff_steps=[1, 3])
    selective = np.ones((8, 3), dtype=bool)
    selective[:4] = False
    add("calcs", "generic/selective_relax", 25, 3, incar={"IBRION": 2, "NSW": 3}, selective=selective)
    species, cell, frames = frames_for(26, 4)
    for i, frame in enumerate(frames):
        frame.update(mag_total=6.8 if i >= 2 else 0.0, site_moments=FM_SITES if i >= 2 else AFM_SITES)
    add("calcs", "generic/afm_transition", 26, 4, incar=dict(MD, NSW=4, ISPIN=2, MAGMOM=AFM_MAGMOM),
        species=species, cell=cell, frames=frames)
    add("calcs", "generic/uncontrolled_nio", 27, 3, incar=dict(MD, NSW=3, ISPIN=2),
        frames=frames_for(27, 3, mag_total=4.0)[2], species=structure(27)[0], cell=structure(27)[1])
    add("calcs", "generic/other_pool", 28, 3, incar=dict(MD, NSW=3, ENCUT=520.0))
    add("calcs", "generic/outcar_only", 29, 2, incar=dict(MD, NSW=2))
    # InterfaceForge-like campaign (lineage via relative_path)
    add("ifc", f"if_campaign/OPT/{CAPPED}", 40, 3, incar={"IBRION": 2, "NSW": 3})
    add("ifc", f"if_campaign/OPT/{BARE}", 41, 3, incar={"IBRION": 2, "NSW": 3})
    species, cell, step1 = frames_for(42, 4)
    add("ifc", f"if_campaign/Step1/{CAPPED}", 42, 4, incar=dict(MD, NSW=4, TEBEG=300), species=species, cell=cell,
        frames=step1)
    add("ifc", f"if_campaign/Step1/{CAPPED}/.interfaceforge/archive/step1_repair_000", 43, 2,
        incar=dict(MD, NSW=2))
    for label, temperature, k in (("300", 300.0, 44), ("450p5", 450.5, 45)):
        _, _, extra = frames_for(k, 2)
        extra = [dict(f, positions=f["positions"] - structure(k)[2] + structure(42)[2]) for f in extra]
        add("ifc", f"if_campaign/Step2_{label}K/{CAPPED}", k, 3, incar=dict(MD, NSW=3, TEBEG=temperature),
            species=species, cell=cell, frames=[dict(step1[-1])] + extra)
    add("ifc", f"if_campaign/Step2_300K/{CAPPED}/precondition", 46, 1, species=species, cell=cell,
        frames=[dict(step1[-1])])
    return specs


def build_tree(base: Path, *, order_seed: int | None = None) -> dict[str, Path]:
    """Write the SYNTHETIC tree under ``base/calcs`` and ``base/ifc`` (creation order shuffled by seed)."""
    roots = {"calcs": base / "calcs", "ifc": base / "ifc"}
    specs = _run_specs()
    if order_seed is not None:
        random.Random(order_seed).shuffle(specs)
    for root, rel, kwargs in specs:
        write_vasp_run(roots[root] / rel, **kwargs)
    calcs, ifc = roots["calcs"], roots["ifc"]
    # agglomeration manifests + launcher continuation
    write_json(calcs / AGGLO / "agglomeration_manifest.json",
               {"replicas": [{"replica": r} for r in CASES], "config_sha256": "c0ffee"})
    for replica, case in CASES.items():
        write_json(calcs / AGGLO / "vasp_runs" / case / "agglomeration_manifest.json", {
            "agglomerate": "pa_dimer", "replica": replica, "packmol_seed": 10 + replica, "center_scale": 1.0,
            "vasp_training_mode": "vasp_md", "composition": [{"slug": "pa", "count": 2}],
            "reference_files": [{"path": "INCAR", "sha256": "ab" * 32}],
        })
    hold = calcs / AGGLO / "vasp_runs" / CASES[0] / "400K" / "02_hold_400K"
    (hold / "runvasp.sh").write_text("#!/bin/bash\n# SYNTHETIC\ncp ../01_heat_300_to_400K/CONTCAR POSCAR\n",
                                     encoding="utf-8", newline="\n")
    # non-canonical near-name, OUTCAR-only, unreadable vasprun
    relax = calcs / "generic" / "relax_b"
    shutil.copyfile(relax / "vasprun.xml", relax / "vasprun.xml.1")
    (calcs / "generic" / "outcar_only" / "vasprun.xml").unlink()
    broken = calcs / "generic" / "broken"
    broken.mkdir(parents=True)
    (broken / "vasprun.xml").write_text("<!-- SYNTHETIC TEST FIXTURE -->\nthis is not xml\n", encoding="utf-8")
    # OutPackLite-like package: outputs without labels
    package = calcs / "XNiO_OutPackLite" / "bulk_md"
    package.mkdir(parents=True)
    for name in ("OSZICAR", "XDATCAR", "INCAR", "CONTCAR"):
        (package / name).write_text(SYNTHETIC_TEXT, encoding="utf-8")
    # InterfaceForge manifests (relative_path join; POSCAR sha edge)
    camp = ifc / "if_campaign"
    write_json(camp / "OPT" / "opt_manifest.json", {"format": "interfaceforge-opt-manifest", "runs": [
        {"relative_path": CAPPED}, {"relative_path": BARE}]})
    write_json(camp / "OPT" / CAPPED / "provenance.json", {
        "schema": "interfaceforge.reactive-surface/v1", "coverage": 0.25, "arrangement": "clustered",
        "motif": "terminal_hydroxyl", "docking": {"mode": "direct"}, "parent_state": "bare"})
    sha = hashlib.sha256((camp / "Step1" / CAPPED / "POSCAR").read_bytes()).hexdigest()
    write_json(camp / "Step1" / "step1_manifest.json", {
        "format": "interfaceforge-step1-series", "source_root": "/work/cluster/if_campaign/OPT",
        "protocol": "training", "temperature_k": 300, "runs": [{"relative_path": CAPPED, "step1_poscar_sha256": sha}]})
    for label, temperature in (("300", 300.0), ("450p5", 450.5)):
        write_json(camp / f"Step2_{label}K" / "step2_manifest.json", {
            "format": "interfaceforge-step2-series", "source_root": "/work/cluster/if_campaign/Step1",
            "runs": [{"relative_path": CAPPED, "temperature_k": temperature}]})
    return roots


def rid(root: str, rel: str) -> str:
    return f"{root}:{rel}"


AGGLO_RUNS = [f"{AGGLO}/vasp_runs/{CASES[0]}/300K", f"{AGGLO}/vasp_runs/{CASES[0]}/400K/01_heat_300_to_400K",
              f"{AGGLO}/vasp_runs/{CASES[0]}/400K/02_hold_400K", f"{AGGLO}/vasp_runs/{CASES[1]}/300K",
              f"{AGGLO}/vasp_runs/{CASES[2]}/300K"]
IF_RUNS = {"opt": f"if_campaign/OPT/{CAPPED}", "bare": f"if_campaign/OPT/{BARE}",
           "step1": f"if_campaign/Step1/{CAPPED}", "step2_300": f"if_campaign/Step2_300K/{CAPPED}",
           "step2_450": f"if_campaign/Step2_450p5K/{CAPPED}"}


def main_options(pool: str | None = None, **policy) -> ex.ExportOptions:
    return ex.ExportOptions(lineage_policy="run", pool=pool, policy=Policy(include_stress=True, **policy))


def root_args(roots: dict[str, Path], order=("calcs", "ifc")) -> list[str]:
    return [f"{name}={roots[name]}" for name in order]


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    base = tmp_path_factory.mktemp("e2e")
    roots = build_tree(base / "tree")
    scan = ex.scan(root_args(roots), options=main_options())
    main_pool = next(p for p in scan.manifest["pools"] if "generic/other_pool" not in " ".join(p["run_ids"]))
    out = base / "export"
    result = ex.export(root_args(roots), out, options=main_options(main_pool["pool_id"]))
    split_dir = base / "split"
    split_manifest = sp.split_export(out, split_dir)
    return {"base": base, "roots": roots, "scan": scan, "pool": main_pool["pool_id"], "out": out,
            "result": result, "split": split_dir, "split_manifest": split_manifest}


def by_id(records):
    return {r.get("frame_id") or r.get("run_id"): r for r in records}


def outcome(record):
    return record["outcome"]["status"], record["outcome"]["reason"]


# --------------------------------------------------------------------------
# accounting: every candidate appears with exactly one outcome
# --------------------------------------------------------------------------

def test_scan_reports_two_pools_and_export_refuses_without_pool(world, tmp_path):
    manifest = world["scan"].manifest
    assert manifest["pool_selection_required"] is True and len(manifest["pools"]) == 2
    other = next(p for p in manifest["pools"] if p["pool_id"] != world["pool"])
    assert other["run_ids"] == [rid("calcs", "generic/other_pool")]
    with pytest.raises(DatasetError, match="reference-settings pools") as info:
        ex.export(root_args(world["roots"]), tmp_path / "out", options=main_options())
    assert "ENCUT" in str(info.value) and other["pool_id"] in str(info.value)
    assert not (tmp_path / "out" / "dataset_manifest.json").exists()
    assert not any(p.name.startswith(".work-") for p in (tmp_path / "out").iterdir())


def test_every_candidate_is_accounted_for(world):
    out = world["out"]
    runs = by_id(read_jsonl(out / "runs.jsonl"))
    expected_runs = {rid("calcs", r) for r in AGGLO_RUNS} | {rid("ifc", r) for r in IF_RUNS.values()} | {
        rid("calcs", f"generic/{name}") for name in (
            "static_a", "relax_b", "md_truncated", "species_mismatch", "scf_fail", "isif0_md", "dup_a", "dup_b",
            "conflict_a", "conflict_b", "compressed_md", "restart_1", "restart_2", "energy_only", "mlff_md",
            "selective_relax", "afm_transition", "uncontrolled_nio", "other_pool", "outcar_only", "broken")}
    pruned = {
        rid("calcs", f"{AGGLO}/vasp_runs/{CASES[0]}.protocol_old/300K"): ("excluded", "archived_path"),
        rid("calcs", f"{AGGLO}/vasp_reference"): ("excluded", "archived_path"),
        rid("calcs", "generic/relax_b/backup_1"): ("excluded", "archived_path"),
        rid("calcs", "generic/old/rerun"): ("excluded", "archived_path"),
        rid("calcs", "XNiO_OutPackLite/bulk_md"): ("excluded", "archived_path"),
        rid("ifc", f"if_campaign/Step1/{CAPPED}/.interfaceforge/archive/step1_repair_000"): ("excluded", "archived_path"),
        rid("ifc", f"if_campaign/Step2_300K/{CAPPED}/precondition"): ("excluded", "preconditioning_run"),
    }
    assert set(runs) == expected_runs | set(pruned)
    for run_id, want in pruned.items():
        assert outcome(runs[run_id]) == want and runs[run_id]["discovery"]["pruned"] is True
    g = lambda name: runs[rid("calcs", f"generic/{name}")]  # noqa: E731
    assert outcome(g("md_truncated")) == ("excluded", "run_incomplete")
    assert outcome(g("species_mismatch")) == ("quarantined", "species_mismatch")
    assert outcome(g("uncontrolled_nio")) == ("quarantined", "magnetic_uncontrolled")
    assert outcome(g("other_pool")) == ("excluded", "reference_pool_not_selected")
    assert outcome(g("outcar_only")) == ("missing_labels", "no_vasprun")
    assert outcome(g("broken")) == ("parse_error", "vasprun_unreadable")
    assert g("compressed_md")["source_file_type"] == "vasprun.gz" and outcome(g("compressed_md"))[0] == "accepted"
    assert g("afm_transition")["assessment"]["magnetic"]["magnetic_class"] == "controlled_transition"
    # every frame: one valid outcome; counts in the manifest agree
    frames = read_jsonl(out / "frames.jsonl")
    manifest = json.loads((out / "dataset_manifest.json").read_text(encoding="utf-8"))
    assert len(frames) == manifest["counts"]["frames"] == sum(r["frames_total"] or 0 for r in runs.values())
    reasons = {}
    for frame in frames:
        key = f"{frame['outcome']['status']}/{frame['outcome']['reason'] or 'accepted'}"
        reasons[key] = reasons.get(key, 0) + 1
    assert reasons == manifest["counts"]["frames_by_reason"]
    # exclusions.csv lists every non-accepted run and frame plus ignored discovery entries
    with (out / "exclusions.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    run_rows = {r["id"] for r in rows if r["level"] == "run"}
    frame_rows = {r["id"] for r in rows if r["level"] == "frame"}
    assert run_rows == {k for k, r in runs.items() if r["outcome"]["status"] != "accepted"}
    assert frame_rows == {f["frame_id"] for f in frames if f["outcome"]["status"] != "accepted"}
    near = [r for r in rows if r["level"] == "discovery:file" and r["reason"] == "not_a_label_file"]
    assert [r["id"] for r in near] == [rid("calcs", "generic/relax_b/vasprun.xml.1")]


def test_frame_outcomes_cover_the_required_failure_modes(world):
    frames = by_id(read_jsonl(world["out"] / "frames.jsonl"))
    f = lambda rel, i, root="calcs": frames[f"{root}:{rel}#{i:05d}"]  # noqa: E731
    # truncated MD (default exclude): complete frames excluded with the run, the tail rejected
    assert [outcome(f("generic/md_truncated", i)) for i in range(4)] == \
        [("excluded", "run_incomplete")] * 3 + [("rejected", "truncated_frame")]
    # unconverged SCF: OUTCAR 'not reached' marker and the NELM ceiling with |dE| >= EDIFF
    assert outcome(f("generic/scf_fail", 1)) == ("rejected", "scf_not_converged")
    assert outcome(f("generic/scf_fail", 2)) == ("rejected", "scf_not_converged")
    assert outcome(f("generic/scf_fail", 0))[0] == outcome(f("generic/scf_fail", 3))[0] == "accepted"
    # missing stress: omitted with a reason, never zeros; valid stress exported when requested
    assert f("generic/isif0_md", 0)["stress_available"] is False
    assert f("generic/isif0_md", 0)["stress_reason"] == "not_computed"
    assert f("generic/static_a", 0)["stress_available"] is True
    # exact duplicate (consistent labels) vs contradictory labels for an identical structure
    assert outcome(f("generic/dup_b", 0)) == ("excluded", "exact_duplicate")
    assert f("generic/dup_b", 0)["duplicate_of"] == rid("calcs", "generic/dup_a") + "#00002"
    assert outcome(f("generic/conflict_a", 1)) == outcome(f("generic/conflict_b", 0)) == \
        ("quarantined", "contradictory_labels")
    assert outcome(f("generic/conflict_a", 0))[0] == "accepted"
    # restart chain: the repeated first frame is a duplicate, and both runs share one lineage group
    assert outcome(f("generic/restart_2", 0)) == ("excluded", "exact_duplicate")
    assert f("generic/restart_1", 0)["lineage_group"] == f("generic/restart_2", 1)["lineage_group"]
    # forces are mandatory by default
    assert outcome(f("generic/energy_only", 0)) == ("missing_labels", "missing_forces")
    # VASP-MLFF steps are counted and indexed, never exported
    assert [outcome(f("generic/mlff_md", i))[0] for i in range(5)] == \
        ["accepted", "non_dft_mlff_step", "accepted", "non_dft_mlff_step", "accepted"]
    assert f("generic/mlff_md", 1)["label_source"] == "mlff"
    # magnetic transition: frames from the first transition on are quarantined
    assert [outcome(f("generic/afm_transition", i)) for i in range(4)] == \
        [("accepted", None)] * 2 + [("quarantined", "magnetic_order_changed")] * 2
    assert f("generic/afm_transition", 3)["mag_d_total_first"] == pytest.approx(6.8)
    assert {outcome(f("generic/uncontrolled_nio", i)) for i in range(3)} == {("quarantined", "magnetic_uncontrolled")}
    # other pool, species mismatch
    assert outcome(f("generic/other_pool", 0)) == ("excluded", "reference_pool_not_selected")
    assert outcome(f("generic/species_mismatch", 0)) == ("quarantined", "species_mismatch")
    # InterfaceForge: Step2 frame 0 == Step1 last frame -> duplicate of the Step1 frame (kept once)
    step1_last = rid("ifc", IF_RUNS["step1"]) + "#00003"
    for name in ("step2_300", "step2_450"):
        assert outcome(f(IF_RUNS[name], 0, "ifc")) == ("excluded", "exact_duplicate")
        assert f(IF_RUNS[name], 0, "ifc")["duplicate_of"] == step1_last
    assert outcome(frames[step1_last])[0] == "accepted"


def test_lineage_groups_follow_manifests_and_links(world):
    frames = [f for f in read_jsonl(world["out"] / "frames.jsonl") if f["outcome"]["status"] == "accepted"]
    group = {}
    for frame in frames:
        group.setdefault(frame["run_id"], set()).add(frame["lineage_group"])
    assert all(len(g) == 1 for g in group.values())
    group = {run_id: next(iter(g)) for run_id, g in group.items()}
    agglo = [group[rid("calcs", r)] for r in AGGLO_RUNS]
    assert agglo[:3] == ["agglo:calcs/agglo/pa_campaign/pa_dimer/r000"] * 3  # 300K + heating + hold
    assert agglo[3:] == ["agglo:calcs/agglo/pa_campaign/pa_dimer/r001", "agglo:calcs/agglo/pa_campaign/pa_dimer/r002"]
    capped = f"iface:ifc/if_campaign/{CAPPED}"
    assert {group[rid("ifc", IF_RUNS[n])] for n in ("opt", "step1", "step2_300", "step2_450")} == {capped}
    assert group[rid("ifc", IF_RUNS["bare"])] == f"iface:ifc/if_campaign/{BARE}"
    assert group[rid("calcs", "generic/static_a")] == "run:calcs:generic/static_a"
    runs = by_id(read_jsonl(world["out"] / "runs.jsonl"))
    step2 = runs[rid("ifc", IF_RUNS["step2_450"])]
    assert step2["metadata"]["temperature_K"] == 450.5 and step2["metadata"]["campaign_stage"] == "step2"
    assert step2["metadata"]["family"] == "OH25" and step2["lineage"]["source"] == "interfaceforge"
    heat = runs[rid("calcs", AGGLO_RUNS[1])]
    assert heat["metadata"]["family"] == "pa_dimer" and heat["metadata"]["campaign"] == AGGLO
    hold_links = runs[rid("calcs", AGGLO_RUNS[2])]["lineage"]["links"]
    assert f"launcher_continuation:{rid('calcs', AGGLO_RUNS[1])}" in hold_links


# --------------------------------------------------------------------------
# dataset.extxyz: exactly the accepted frames, exact labels, metadata contract
# --------------------------------------------------------------------------

def test_extxyz_holds_exactly_the_accepted_frames_with_exact_labels(world):
    out = world["out"]
    frames = read_jsonl(out / "frames.jsonl")
    accepted = [f["frame_id"] for f in frames if f["outcome"]["status"] == "accepted"]
    assert read_frame_ids(out / "dataset.extxyz") == accepted
    assert not any("mlff" in f["label_source"] for f in frames if f["frame_id"] in set(accepted))
    atoms_list = ase.io.read(out / "dataset.extxyz", index=":", format="extxyz")
    by_frame = {a.info["frame_id"]: a for a in atoms_list}
    for atoms in atoms_list:
        assert atoms.calc is None and not atoms.constraints and tuple(atoms.pbc) == (True, True, True)
        for key in ("source", "source_file_type", "ionic_step", "structure_key", "lineage_group", "energy_source",
                    "energy_rule", "stress_available", "scf_status", "label_source", "vasp_version",
                    "magnetic_class", "magnetic_policy", "campaign", "family", "parser", "parser_version",
                    "repo_commit", "label_set"):
            assert key in atoms.info, key
        assert atoms.info["scf_status"] == "converged" and atoms.info["label_source"] == "dft"
        assert ("REF_stress" in atoms.info) == atoms.info["stress_available"]
        assert ("stress_source" in atoms.info) != ("stress_reason" in atoms.info)
        assert "split" not in atoms.info
    # exact labels against the SYNTHETIC writer's parsed-back expectations
    species, cell, frames_in = frames_for(11, 4)
    check = write_vasp_run(world["base"] / "check_relax_b", species=species, cell=cell, frames=frames_in,
                           incar=spin_incar({"IBRION": 2, "NSW": 4}))
    for i, step in enumerate(check["expected"]["steps"]):
        atoms = by_frame[rid("calcs", "generic/relax_b") + f"#{i:05d}"]
        assert atoms.info["REF_energy"] == step["free_energy"]
        assert np.array_equal(atoms.arrays["REF_forces"], step["forces"])
        assert np.array_equal(atoms.positions, step["positions"])
        assert atoms.info["source"] == "calcs:generic/relax_b/vasprun.xml" and atoms.info["ionic_step"] == i
        assert atoms.info["energy_source"] == "vasprun:calculation.e_fr_energy-PSTRESS*V"
    # selective dynamics: raw direct-basis flags as an array, raw (non-zeroed) forces, no constraints
    sel = by_frame[rid("calcs", "generic/selective_relax") + "#00001"]
    flags = sel.arrays["vasp_selective_dynamics"]
    assert flags.dtype == bool and not flags[:4].any() and flags[4:].all()
    assert sel.info["selective_dynamics_basis"] == "direct" and np.abs(sel.arrays["REF_forces"][:4]).sum() > 0
    # magnetic metadata on a spin-polarized frame
    afm = by_frame[rid("calcs", "generic/afm_transition") + "#00001"]
    assert afm.info["magnetic_class"] == "controlled_transition" and afm.info["magnetic_policy"] == "accepted"
    assert np.allclose(afm.arrays["vasp_magmom_initial"], [2, -2, 2, -2, 0, 0, 0, 0])
    assert np.allclose(afm.arrays["vasp_magmom_final"], AFM_SITES)
    heat = by_frame[rid("calcs", AGGLO_RUNS[1]) + "#00000"]
    assert heat.info["temperature_K"] == "300->400" and heat.info["family"] == "pa_dimer"


def test_manifest_hashes_and_content_sha(world):
    out = world["out"]
    manifest = ex.load_manifest(out)
    for name, entry in manifest["files"].items():
        assert sha256_file(out / name) == entry["sha256"] and (out / name).stat().st_size == entry["bytes"]
    assert set(manifest["files"]) == {"dataset.extxyz", "runs.jsonl", "frames.jsonl", "exclusions.csv",
                                      "audit_report.json", "audit.md"}
    assert manifest["files"]["dataset.extxyz"]["round_trip_verified"] is True
    assert manifest["policy"]["non_default"] == {"include_stress": True}
    assert manifest["pool"]["pool_id"] == world["pool"] and manifest["options"]["lineage_policy"] == "run"
    assert manifest["parser"]["name"] == "nio-md-prep.vasprun-stream" and "git" in manifest["tool"]
    assert manifest["label_keys"] == {"energy": "REF_energy", "forces": "REF_forces", "stress": "REF_stress"}
    assert any("real NiO force validation is outstanding" in item for item in manifest["limitations"])
    assert "created_at" in manifest["invocation"]
    edited = dict(manifest, energy_quantity="energy_sigma0")
    assert ex._content_sha256(edited) != manifest["content_sha256"]
    summary = json.loads((out / "audit_report.json").read_text(encoding="utf-8"))
    for dimension in ("status", "reason", "campaign", "family", "temperature_K", "composition", "lineage_group",
                      "vasp_version", "forces_available", "stress_available", "magnetic_class"):
        assert dimension in summary["frames"], dimension
    assert "magnetic_class" in summary["runs"] and "## Frames" in (out / "audit.md").read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# split + audit + reload: zero leakage
# --------------------------------------------------------------------------

def test_split_audit_and_reload_have_zero_leakage(world, tmp_path):
    out, split_dir = world["out"], world["split"]
    frames = by_id(r for r in read_jsonl(out / "frames.jsonl") if r["outcome"]["status"] == "accepted")
    seen: dict[str, str] = {}
    owners: dict[tuple[str, str], str] = {}
    for name in ("train", "valid", "test"):
        for atoms in ase.io.iread(str(split_dir / f"{name}.extxyz"), index=":", format="extxyz"):
            frame_id = atoms.info["frame_id"]
            assert atoms.info["split"] == name and frame_id not in seen
            seen[frame_id] = name
            for kind in ("lineage_group", "structure_key", "permutation_key"):
                assert owners.setdefault((kind, atoms.info[kind]), name) == name, (kind, frame_id)
            assert frames[frame_id]["label_sha256"]
    assert set(seen) == set(frames)
    assert len({seen[f] for f in frames if frames[f]["run_id"] in {rid("calcs", r) for r in AGGLO_RUNS[:3]}}) == 1
    report = au.audit_export(out, split_dirs=[split_dir], verify_sources={k: str(v) for k, v in world["roots"].items()},
                             output=tmp_path / "audit")
    failed = [c for c in report["checks"] if not c["ok"]]
    assert report["ok"] and not failed, failed
    assert {c["name"] for c in report["checks"]} >= {"manifest_content_sha256", "file_hashes", "accounting",
                                                     "ase_round_trip", "group_structure_consistency",
                                                     "split:split", "source_hashes"}
    assert report["round_trip"]["frames_read"] == len(frames)
    # split reproducibility: same inputs -> identical bytes
    again = sp.split_export(out, tmp_path / "split_again")
    for name in ("train.extxyz", "valid.extxyz", "test.extxyz", "split_manifest.json", "split_summary.md"):
        assert (tmp_path / "split_again" / name).read_bytes() == (split_dir / name).read_bytes(), name
    assert again["content_sha256"] == world["split_manifest"]["content_sha256"]


def test_correlated_frames_cannot_land_in_different_splits(world, tmp_path):
    """One AIMD trajectory, one agglomeration replica, one InterfaceForge structure, one restart chain:
    each is one split, for every seed (and the split is not trivially 'everything in train')."""
    correlated = {
        "agglomeration replica 0 (300 K run + heating + hold)": [rid("calcs", r) for r in AGGLO_RUNS[:3]],
        "InterfaceForge capped structure (OPT, Step1, Step2 300 K, Step2 450.5 K)":
            [rid("ifc", IF_RUNS[n]) for n in ("opt", "step1", "step2_300", "step2_450")],
        "restart chain": [rid("calcs", "generic/restart_1"), rid("calcs", "generic/restart_2")],
        "one MD trajectory": [rid("calcs", "generic/compressed_md")],
    }
    accepted = [f for f in read_jsonl(world["out"] / "frames.jsonl") if f["outcome"]["status"] == "accepted"]
    runs_of = {}
    for frame in accepted:
        runs_of.setdefault(frame["run_id"], []).append(frame["frame_id"])
    for members in correlated.values():
        assert all(runs_of.get(run_id) for run_id in members), members  # every member contributes frames
    used_splits = set()
    for seed in (0, 1, 2, 3, 11, 12345):
        manifest = sp.split_export(world["out"], tmp_path / f"seed{seed}", seed=seed)
        where = {frame_id: name for name, ids in manifest["frame_ids"].items() for frame_id in ids}
        assert set(where) == {f["frame_id"] for f in accepted}
        for label, members in correlated.items():
            splits = {where[frame_id] for run_id in members for frame_id in runs_of[run_id]}
            assert len(splits) == 1, (seed, label, splits)
            used_splits |= splits
        assert sum(1 for ids in manifest["frame_ids"].values() if ids) >= 2, seed
    assert len(used_splits) >= 2  # the correlated sets really move between splits across seeds


def test_audit_detects_tampering(world, tmp_path):
    copy = tmp_path / "export"
    shutil.copytree(world["out"], copy)
    text = (copy / "frames.jsonl").read_text(encoding="utf-8").replace('"exported":true', '"exported":false', 1)
    (copy / "frames.jsonl").write_text(text, encoding="utf-8", newline="\n")
    report = au.audit_export(copy)
    assert not report["ok"] and not next(c for c in report["checks"] if c["name"] == "file_hashes")["ok"]
    assert Path(report["output"]) == copy / "audit" and (copy / "audit" / "audit.md").is_file()


def test_split_refuses_unresolved_lineage_and_impossible_strata(world, tmp_path):
    with pytest.raises(DatasetError, match="stratum|strata"):
        sp.split_export(world["out"], tmp_path / "strata", stratify_by="campaign")
    copy = tmp_path / "export"
    shutil.copytree(world["out"], copy)
    text = (copy / "dataset.extxyz").read_text(encoding="ascii")
    text = text.replace('lineage_source="policy:run"', 'lineage_source="unresolved"', 1)
    (copy / "dataset.extxyz").write_text(text, encoding="ascii", newline="\n")
    with pytest.raises(DatasetError, match="unresolved"):
        sp.split_export(copy, tmp_path / "split")


def test_unresolved_lineage_is_quarantined_not_exported(world):
    scan = ex.scan(root_args(world["roots"]), options=ex.ExportOptions(pool=world["pool"]))
    frames = by_id(ex._frame_record(f) for f in scan.analysis.frames)
    runs = {r.run_id: r for r in scan.analysis.runs}
    static = runs[rid("calcs", "generic/static_a")]
    assert (static.outcome.status, static.outcome.reason) == ("quarantined", "lineage_unresolved")
    assert "--lineage-policy" in static.outcome.detail
    assert outcome(frames[rid("calcs", "generic/static_a") + "#00000"]) == ("quarantined", "lineage_unresolved")
    assert outcome(frames[rid("ifc", IF_RUNS["opt"]) + "#00000"])[0] == "accepted"  # manifests resolve lineage
    assert outcome(frames[rid("calcs", AGGLO_RUNS[0]) + "#00000"])[0] == "accepted"
    assert rid("calcs", "generic/static_a") in scan.manifest["lineage"]["unresolved_runs"]


def test_frozen_test_split_on_append(world, tmp_path):
    first = tmp_path / "v1"
    ex.export([f"calcs={world['roots']['calcs']}"], first, options=main_options(world["pool"]))
    v1 = sp.split_export(first, tmp_path / "v1_split")
    v2 = sp.split_export(world["out"], tmp_path / "v2_split", previous=tmp_path / "v1_split" / "split_manifest.json")
    assert set(v1["frame_ids"]["test"]) == set(v2["frame_ids"]["test"])  # frozen: no new frames in test
    new = set(v2["frame_ids"]["train"] + v2["frame_ids"]["valid"]) - set(v1["frame_ids"]["train"] + v1["frame_ids"]["valid"])
    assert new and all(frame_id.startswith("ifc:") for frame_id in new)


# --------------------------------------------------------------------------
# determinism, overwrite refusal
# --------------------------------------------------------------------------

def _strip(manifest):
    return {k: v for k, v in manifest.items() if k not in ("invocation", "content_sha256", "roots")}


def test_deterministic_bytes_across_root_order_creation_order_and_repeats(world, tmp_path):
    options = main_options(world["pool"])
    repeat = ex.export(root_args(world["roots"]), tmp_path / "repeat", options=options)
    shuffled_roots = build_tree(tmp_path / "tree2", order_seed=7)
    other = ex.export(root_args(shuffled_roots, order=("ifc", "calcs")), tmp_path / "other", options=options)
    names = ["dataset.extxyz", "runs.jsonl", "frames.jsonl", "exclusions.csv", "audit_report.json", "audit.md"]
    for name in names:
        reference = (world["out"] / name).read_bytes()
        assert (tmp_path / "repeat" / name).read_bytes() == reference, name
        assert (tmp_path / "other" / name).read_bytes() == reference, name
    assert repeat.manifest["content_sha256"] == world["result"].manifest["content_sha256"]
    assert _strip(other.manifest) == _strip(world["result"].manifest)


def test_overwrite_is_refused_and_force_replaces_in_place(world, tmp_path):
    out = tmp_path / "out"
    options = main_options(world["pool"])
    ex.export(root_args(world["roots"]), out, options=options)
    (out / "notes.txt").write_text("user file\n", encoding="utf-8")
    with pytest.raises(OverwriteRefusedError):
        ex.export(root_args(world["roots"]), out, options=options)
    again = ex.export(root_args(world["roots"]), out, options=options, force=True)
    assert (out / "notes.txt").read_text(encoding="utf-8") == "user file\n"  # unrelated files untouched
    assert (out / "dataset.extxyz").read_bytes() == (world["out"] / "dataset.extxyz").read_bytes()
    assert again.manifest["content_sha256"] == world["result"].manifest["content_sha256"]
    assert not any(p.name.startswith(".") for p in out.iterdir())


# --------------------------------------------------------------------------
# deliberate policy choices
# --------------------------------------------------------------------------

def test_deliberate_inclusion_recovery_overrides_and_energy_only(world, tmp_path):
    roots = world["roots"]
    options = ex.ExportOptions(
        lineage_policy="run", pool=world["pool"], include_parts=("precondition", "X*"),
        policy=Policy(interrupted_run_policy="recover", accept_uncontrolled=True, allow_energy_only=True,
                      magnetic_override_reason="SYNTHETIC test: reviewed ferromagnetic start",
                      near_duplicate_tolerance=1e-4),
    )
    result = ex.export(root_args(roots), tmp_path / "out", options=options)
    runs = {r.run_id: r for r in result.analysis.runs}
    frames = by_id(read_jsonl(tmp_path / "out" / "frames.jsonl"))
    precondition = rid("ifc", f"if_campaign/Step2_300K/{CAPPED}/precondition")
    assert "included_excluded_part:precondition" in runs[precondition].record["discovery"]["flags"]
    assert outcome(frames[precondition + "#00000"]) == ("excluded", "exact_duplicate")
    package = runs[rid("calcs", "XNiO_OutPackLite/bulk_md")]
    assert (package.outcome.status, package.outcome.reason) == ("missing_labels", "no_vasprun")
    assert [outcome(frames[rid("calcs", "generic/md_truncated") + f"#{i:05d}"])[0] for i in range(4)] == \
        ["accepted"] * 3 + ["rejected"]
    unc = frames[rid("calcs", "generic/uncontrolled_nio") + "#00000"]
    assert unc["outcome"]["status"] == "accepted"
    assert unc["magnetic_policy"] == "accepted_by_override:accept_uncontrolled"
    energy_only = frames[rid("calcs", "generic/energy_only") + "#00001"]
    assert energy_only["outcome"]["status"] == "accepted" and energy_only["label_set"] == "energy_only"
    atoms = {a.info["frame_id"]: a for a in ase.io.read(tmp_path / "out" / "dataset.extxyz", index=":")}
    eo = atoms[rid("calcs", "generic/energy_only") + "#00001"]
    assert "REF_forces" not in eo.arrays and eo.info["label_set"] == "energy_only" and eo.info["forces_reason"]
    manifest = result.manifest
    assert manifest["policy"]["non_default"]["accept_uncontrolled"] is True
    assert manifest["policy"]["magnetic_override_reason"].startswith("SYNTHETIC")
    assert manifest["options"]["include_parts"] == ["X*", "precondition"]
    report = au.audit_export(tmp_path / "out")
    assert report["ok"], [c for c in report["checks"] if not c["ok"]]
    with pytest.raises(DatasetError, match="magnetic-override-reason"):
        ex.ExportOptions(policy=Policy(accept_uncontrolled=True))


def test_inventory_selection_and_lineage(world, tmp_path):
    inventory = tmp_path / "inventory.toml"
    inventory.write_text(
        'select = "listed-only"\n'
        '[[run]]\npath = "calcs:generic/static_a"\nlineage = "nio-static"\nfamily = "bulk"\ncampaign = "inv"\n'
        '[[run]]\npath = "calcs:generic/relax_b"\ninclude = false\n', encoding="utf-8")
    scan = ex.scan(root_args(world["roots"]), options=ex.ExportOptions(inventory=inventory))
    runs = {r.run_id: r for r in scan.analysis.runs}
    static = runs[rid("calcs", "generic/static_a")]
    assert static.outcome.status == "accepted" and static.lineage["lineage_id"] == "inv:nio-static"
    assert static.metadata["family"] == "bulk" and static.metadata["campaign"] == "inv"
    relax = runs[rid("calcs", "generic/relax_b")]
    assert (relax.outcome.status, relax.outcome.reason) == ("excluded", "inventory_excluded")
    other = runs[rid("calcs", "generic/dup_a")]
    assert (other.outcome.status, other.outcome.reason) == ("excluded", "not_in_inventory")
    assert scan.manifest["inventory"]["sha256"] == sha256_file(inventory)


# --------------------------------------------------------------------------
# CLI and import hygiene
# --------------------------------------------------------------------------

def test_cli_exit_codes(world, tmp_path, capsys):
    roots = world["roots"]
    common = ["--root", f"calcs={roots['calcs']}", "--root", f"ifc={roots['ifc']}", "--lineage-policy", "run",
              "--include-stress"]
    assert top_cli.main(["dataset", "scan", *common, "--output", str(tmp_path / "scan")]) == 0
    assert (tmp_path / "scan" / "scan_manifest.json").is_file() and (tmp_path / "scan" / "frames.jsonl").is_file()
    assert "reference-settings pools: 2" in capsys.readouterr().out
    with pytest.raises(SystemExit) as info:
        top_cli.main(["dataset", "export", *common, "--output", str(tmp_path / "bad")])
    assert info.value.code == 2 and "--pool" in capsys.readouterr().err
    export = ["dataset", "export", *common, "--pool", world["pool"], "--output", str(tmp_path / "out")]
    assert top_cli.main(export) == 0
    assert (tmp_path / "out" / "dataset.extxyz").read_bytes() == (world["out"] / "dataset.extxyz").read_bytes()
    manifest = ex.load_manifest(tmp_path / "out")
    assert manifest["invocation"]["argv"] == export
    with pytest.raises(SystemExit) as info:
        top_cli.main(export)
    assert info.value.code == 2 and "--force" in capsys.readouterr().err
    assert top_cli.main(["dataset", "split", str(tmp_path / "out"), "--output", str(tmp_path / "split")]) == 0
    assert top_cli.main(["dataset", "audit", str(tmp_path / "out"), "--split", str(tmp_path / "split"),
                         "--verify-sources", f"calcs={roots['calcs']}", "--verify-sources", f"ifc={roots['ifc']}"]) == 0
    (tmp_path / "out" / "exclusions.csv").write_text("tampered\n", encoding="utf-8")
    assert top_cli.main(["dataset", "audit", str(tmp_path / "out")]) == 1
    with pytest.raises(SystemExit) as info:
        top_cli.main(["dataset", "scan", *common, "--accept-uncontrolled-magnetism"])
    assert info.value.code == 2 and "magnetic-override-reason" in capsys.readouterr().err
    with pytest.raises(SystemExit) as info:
        top_cli.main(["dataset", "export", "--root", f"calcs={roots['calcs']}", "--output", str(tmp_path / "x"),
                      "--energy-key", "energy"])
    assert info.value.code == 2 and "reserved" in capsys.readouterr().err


def test_spotcheck_script_agrees_with_ase_vasp_reader(world, tmp_path):
    result = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "dataset_spotcheck.py"), str(world["out"]), "--n", "12",
         "--json", str(tmp_path / "spot.json")],
        capture_output=True, text=True, cwd=str(tmp_path), timeout=300,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    summary = json.loads((tmp_path / "spot.json").read_text(encoding="utf-8"))
    assert summary["ok"] and summary["checked"] == 12
    assert summary["max_d_free_energy"] <= 1e-6 and summary["max_d_forces"] <= 1e-6


def test_import_and_help_do_not_pull_heavy_dependencies(tmp_path):
    code = (
        "import sys, contextlib, io\n"
        "import nio_md_prep.dataset\n"
        "from nio_md_prep.dataset import cli as dcli\n"
        "from nio_md_prep import cli\n"
        "buf = io.StringIO()\n"
        "with contextlib.redirect_stdout(buf):\n"
        "    try:\n"
        "        cli.main(['dataset', '--help'])\n"
        "    except SystemExit as exc:\n"
        "        assert exc.code == 0, exc.code\n"
        "assert 'scan' in buf.getvalue() and 'audit' in buf.getvalue()\n"
        "heavy = sorted(m for m in sys.modules if m.split('.')[0] in {'ase', 'torch', 'mace', 'openmm', 'lammps'})\n"
        "print('HEAVY=' + ','.join(heavy))\n"
    )
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, cwd=str(tmp_path),
                            timeout=120)
    assert result.returncode == 0, result.stderr
    assert "HEAVY=\n" in result.stdout or result.stdout.strip().endswith("HEAVY="), result.stdout


def test_leakage_self_check_refuses_shared_structures(world):
    frames = [f for f in world["result"].analysis.frames if f.accepted][:2]
    a, b = frames
    original = b.order_key
    try:
        b.order_key = a.order_key
        if a.run.lineage.get("lineage_group") != b.run.lineage.get("lineage_group"):
            with pytest.raises(LeakageError):
                ex.check_group_consistency([a, b])
    finally:
        b.order_key = original
    assert ex.check_group_consistency([a])
