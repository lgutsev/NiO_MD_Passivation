"""Tests for dataset/split.py (lineage-group splits).

All frame records and extxyz files here are SYNTHETIC test data generated
in-file (random geometry and labels, invented run/group ids); they are not
VASP output and carry no physical meaning.
"""

from __future__ import annotations

import copy
import json
import random

import numpy as np
import pytest

from nio_md_prep.dataset import split as S
from nio_md_prep.dataset.errors import DatasetError, LeakageError, OverwriteRefusedError

SYNTHETIC = "SYNTHETIC-TEST-FIXTURE"


def _records(spec, *, seed=0, metadata=None):
    """SYNTHETIC frame records: spec = {group: {run: n_frames}}."""
    rng = random.Random(seed)
    records = []
    for group, runs in spec.items():
        for run, n_frames in runs.items():
            for index in range(n_frames):
                meta = {"generator": SYNTHETIC, "campaign": "synthetic"}
                if metadata:
                    meta.update(metadata(group, run, index))
                records.append({
                    "frame_id": f"synth:{run}#{index:05d}",
                    "run_id": f"synth:{run}",
                    "lineage_group": group,
                    "n_atoms": 4,
                    "structure_key": f"{group}/{run}/{index}/{rng.random()}",
                    "metadata": meta,
                })
    return records


def _random_spec(rng, n_groups):
    """Random correlated structure: groups hold 1-3 runs (restart chains, replicas), runs 1-25 frames."""
    spec = {}
    for g in range(n_groups):
        spec[f"g{g:03d}"] = {f"g{g:03d}/run{r}": rng.randint(1, 25) for r in range(rng.randint(1, 3))}
    return spec


def _split_of_frames(manifest):
    return {frame_id: split for split in S.SPLITS for frame_id in manifest["frame_ids"][split]}


# --------------------------------------------------------------------------
# the core guarantee: correlated frames never cross splits
# --------------------------------------------------------------------------

def test_correlated_frames_cannot_land_in_different_splits_over_many_seeds():
    """Property-style: for 300 random datasets and seeds, every lineage group and every run is in ONE split."""
    for trial in range(300):
        rng = random.Random(trial)
        spec = _random_spec(rng, rng.randint(3, 40))
        records = _records(spec, seed=trial)
        rng.shuffle(records)
        fractions = rng.choice([(0.8, 0.1, 0.1), (0.7, 0.15, 0.15), (0.6, 0.2, 0.2), (0.9, 0.1, 0.0)])
        manifest = S.make_split(records, seed=trial, fractions=fractions)
        owner = _split_of_frames(manifest)
        assert set(owner) == {record["frame_id"] for record in records}  # every frame exactly once
        by_group, by_run = {}, {}
        for record in records:
            by_group.setdefault(record["lineage_group"], set()).add(owner[record["frame_id"]])
            by_run.setdefault(record["run_id"], set()).add(owner[record["frame_id"]])
        assert all(len(splits) == 1 for splits in by_group.values()), trial
        assert all(len(splits) == 1 for splits in by_run.values()), trial
        assert all(manifest["group_to_split"][g] == next(iter(s)) for g, s in by_group.items())
        for split, fraction in zip(S.SPLITS, fractions):
            if fraction > 0 and len(spec) >= 3:
                assert manifest["achieved"][split]["groups"] >= 1, (trial, split)
            if fraction == 0:
                assert manifest["achieved"][split]["frames"] == 0
        S.check_split_disjoint(manifest)


def test_frames_of_one_run_in_two_groups_raise_leakage():
    records = _records({"gA": {"runA": 3}, "gB": {"runB": 3}, "gC": {"runC": 3}})
    records[0]["lineage_group"] = "gZ"  # one frame of runA claims another group
    with pytest.raises(LeakageError, match="different lineage groups"):
        S.make_split(records)


def test_identical_structure_in_two_groups_raises_leakage():
    records = _records({"gA": {"runA": 3}, "gB": {"runB": 3}, "gC": {"runC": 3}})
    records[4]["structure_key"] = records[0]["structure_key"]  # exact duplicate across groups
    with pytest.raises(LeakageError, match="identical structure"):
        S.make_split(records)
    records = _records({"gA": {"runA": 3}, "gB": {"runB": 3}, "gC": {"runC": 3}})
    for record in records:
        record["permutation_key"] = "same-for-all" if record["run_id"] != "synth:runC" else record["frame_id"]
    with pytest.raises(LeakageError, match="permutation_key"):
        S.make_split(records)


def test_duplicates_inside_one_group_are_fine():
    records = _records({"gA": {"runA": 3, "runA2": 2}, "gB": {"runB": 3}, "gC": {"runC": 3}})
    records[3]["structure_key"] = records[0]["structure_key"]  # runA2 frame 0 == runA frame 0 (restart)
    manifest = S.make_split(records)
    assert manifest["lineage_checks"]["runs"] == 4


# --------------------------------------------------------------------------
# reproducibility
# --------------------------------------------------------------------------

def test_split_is_reproducible_and_independent_of_input_order():
    spec = _random_spec(random.Random(5), 25)
    records = _records(spec, seed=1)
    first = S.make_split(records, seed=7)
    for shuffle_seed in range(5):
        shuffled = copy.deepcopy(records)
        random.Random(shuffle_seed).shuffle(shuffled)
        again = S.make_split(shuffled, seed=7)
        assert again == first
        assert json.dumps(again, sort_keys=True) == json.dumps(first, sort_keys=True)
    assignments = {json.dumps(S.make_split(records, seed=seed)["group_to_split"], sort_keys=True) for seed in range(20)}
    assert len(assignments) > 1  # the seed matters (tie-breaks), the input order does not


def test_default_seed_and_fraction_contract():
    records = _records({f"g{i}": {f"r{i}": 10} for i in range(10)})
    manifest = S.make_split(records)
    assert manifest["seed"] == 11 and manifest["requested_fractions"] == {"train": 0.8, "valid": 0.1, "test": 0.1}
    assert manifest["unit"] == "lineage_group" and manifest["rng"] == "numpy.random.PCG64"
    assert {k: v["frames"] for k, v in manifest["achieved"].items()} == {"train": 80, "valid": 10, "test": 10}
    for bad in [(0.5, 0.5), (0.8, 0.1, 0.2), (0.0, 0.5, 0.5), (1.1, -0.05, -0.05)]:
        with pytest.raises(DatasetError, match="fraction"):
            S.make_split(records, fractions=bad)
    for bad in (-1, 1.5, True, "3"):
        with pytest.raises(DatasetError, match="seed"):
            S.make_split(records, seed=bad)


# --------------------------------------------------------------------------
# fail-closed inputs
# --------------------------------------------------------------------------

def test_unresolved_lineage_is_refused_with_advice():
    records = _records({"gA": {"runA": 2}, "gB": {"runB": 2}, "gC": {"runC": 2}})
    records[0]["lineage_group"] = None
    with pytest.raises(DatasetError, match="lineage-policy"):
        S.make_split(records)
    records = _records({"gA": {"runA": 2}, "gB": {"runB": 2}, "gC": {"runC": 2}})
    records[2]["lineage"] = {"source": "unresolved"}
    with pytest.raises(DatasetError, match="unresolved"):
        S.make_split(records)


def test_non_accepted_duplicate_ids_and_too_few_groups_are_refused():
    records = _records({"gA": {"runA": 2}, "gB": {"runB": 2}, "gC": {"runC": 2}})
    records[1]["outcome"] = {"status": "quarantined", "reason": "magnetic_uncontrolled"}
    with pytest.raises(DatasetError, match="non-accepted"):
        S.make_split(records)
    records = _records({"gA": {"runA": 2}, "gB": {"runB": 2}, "gC": {"runC": 2}})
    with pytest.raises(DatasetError, match="duplicate frame_id"):
        S.make_split(records + [records[0]])
    with pytest.raises(DatasetError, match="at least 3"):
        S.make_split(_records({"gA": {"runA": 50}, "gB": {"runB": 50}}))
    with pytest.raises(DatasetError, match="no frames"):
        S.make_split([])


def test_group_id_is_an_alias_but_must_agree_with_lineage_group():
    records = _records({"gA": {"runA": 2}, "gB": {"runB": 2}, "gC": {"runC": 2}})
    for record in records:
        record["group_id"] = record.pop("lineage_group")
    assert S.make_split(records)["totals"]["groups"] == 3
    records[0]["lineage_group"] = "other"
    with pytest.raises(DatasetError, match="disagree"):
        S.make_split(records)


# --------------------------------------------------------------------------
# stratification and statistics
# --------------------------------------------------------------------------

def _temperature(group, run, index):
    return {"temperature_K": 300 if int(group[1:]) % 2 else 600}


def test_stratified_split_places_every_stratum_in_every_split():
    records = _records({f"g{i}": {f"r{i}": 5 + i % 4} for i in range(12)}, metadata=_temperature)
    manifest = S.make_split(records, stratify_by="temperature_K")
    assert set(manifest["strata"]) == {"300", "600"}
    for stratum in manifest["strata"].values():
        assert all(stratum["per_split"][split]["groups"] >= 1 for split in S.SPLITS)
    stats = manifest["statistics"]["temperature_K"]
    assert sum(stats["300"][s]["frames"] + stats["600"][s]["frames"] for s in S.SPLITS) == len(records)


def test_stratification_errors_are_explicit():
    records = _records({f"g{i}": {f"r{i}": 3} for i in range(6)}, metadata=_temperature)
    records[0]["metadata"]["temperature_K"] = 999  # varies inside group g0
    with pytest.raises(DatasetError, match="not constant"):
        S.make_split(records, stratify_by="temperature_K")
    records = _records({f"g{i}": {f"r{i}": 3} for i in range(6)}, metadata=_temperature)
    del records[0]["metadata"]["temperature_K"]
    for record in records:
        if record["lineage_group"] == "g0":
            record["metadata"].pop("temperature_K", None)
    with pytest.raises(DatasetError, match="missing"):
        S.make_split(records, stratify_by="temperature_K")
    records = _records({f"g{i}": {f"r{i}": 3} for i in range(5)},
                       metadata=lambda g, r, i: {"temperature_K": 450 if g == "g0" else 300})
    with pytest.raises(DatasetError, match="fewer groups"):
        S.make_split(records, stratify_by="temperature_K")
    manifest = S.make_split(records, stratify_by="temperature_K", allow_underpopulated_strata=True)
    assert manifest["strata"]["450"]["underpopulated"] and manifest["allow_underpopulated_strata"]
    assert any("allowed" in warning for warning in manifest["warnings"])


def test_statistics_and_coverage_warnings():
    records = _records({f"g{i}": {f"r{i}": 10} for i in range(10)},
                       metadata=lambda g, r, i: {"family": "n02" if g in ("g0", "g1", "g2") else "n04",
                                                 "magnetic_class": "controlled_consistent"})
    manifest = S.make_split(records, seed=1)
    family = manifest["statistics"]["family"]
    assert set(family) == {"n02", "n04"}
    assert sum(family[v][s]["frames"] for v in family for s in S.SPLITS) == 100
    assert manifest["statistics"]["magnetic_class"]["controlled_consistent"]["train"]["frames"] == 80
    assert manifest["statistics"]["vasp_version"] == {
        S.MISSING: {s: {"frames": manifest["achieved"][s]["frames"], "groups": manifest["achieved"][s]["groups"]}
                    for s in S.SPLITS}}
    assert all(w.startswith("coverage:") or "fraction" in w for w in manifest["warnings"])


# --------------------------------------------------------------------------
# append: frozen test, merged groups
# --------------------------------------------------------------------------

def test_frozen_test_on_append_and_new_groups_only_to_train_valid():
    spec = {f"g{i:02d}": {f"r{i:02d}": 5} for i in range(10)}
    first = S.make_split(_records(spec), seed=3)
    test_before = set(first["frame_ids"]["test"])
    assert test_before
    grown = dict(spec)
    grown.update({f"n{i:02d}": {f"new{i:02d}": 5} for i in range(10)})
    second = S.make_split(_records(grown), seed=3, previous=first)
    assert set(second["frame_ids"]["test"]) == test_before  # frozen: nothing added, nothing removed
    for gid, split in first["group_to_split"].items():
        assert second["group_to_split"][gid] == split  # previous assignments kept
    new_splits = {second["group_to_split"][gid] for gid in grown if gid.startswith("n")}
    assert new_splits <= {"train", "valid"} and second["new_group_fractions"]["test"] == 0.0
    assert second["previous"]["content_sha256"] == first["content_sha256"]
    assert second["previous"]["frozen_splits"] == ["test"]

    third = S.make_split(_records(grown), seed=3, previous=first, grow_test=True)
    assert set(third["frame_ids"]["test"]) >= test_before


def test_frozen_test_group_that_gains_frames_is_refused_unless_grow_test():
    spec = {f"g{i:02d}": {f"r{i:02d}": 5} for i in range(10)}
    first = S.make_split(_records(spec), seed=3)
    test_group = next(gid for gid, split in first["group_to_split"].items() if split == "test")
    spec[test_group] = {**spec[test_group], f"{test_group}-restart": 4}  # a continuation run appears later
    with pytest.raises(DatasetError, match="frozen test group"):
        S.make_split(_records(spec), seed=3, previous=first)
    grown = S.make_split(_records(spec), seed=3, previous=first, grow_test=True)
    assert grown["group_to_split"][test_group] == "test"


def test_groups_merged_across_previous_splits_raise_leakage():
    spec = {f"g{i:02d}": {f"r{i:02d}": 5} for i in range(10)}
    first = S.make_split(_records(spec), seed=3)
    train_group = next(g for g, s in first["group_to_split"].items() if s == "train")
    test_group = next(g for g, s in first["group_to_split"].items() if s == "test")
    records = _records(spec)
    for record in records:
        if record["lineage_group"] == test_group:
            record["lineage_group"] = train_group  # a new lineage link merged the two groups
    with pytest.raises(LeakageError, match="merged"):
        S.make_split(records, seed=3, previous=first)


def test_edited_previous_manifest_is_refused(tmp_path):
    spec = {f"g{i:02d}": {f"r{i:02d}": 5} for i in range(6)}
    first = S.make_split(_records(spec), seed=3)
    tampered = copy.deepcopy(first)
    moved = tampered["frame_ids"]["test"].pop()
    tampered["frame_ids"]["train"].append(moved)
    with pytest.raises(DatasetError, match="content_sha256"):
        S.make_split(_records(spec), previous=tampered)
    path = tmp_path / "split_manifest.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    with pytest.raises(DatasetError, match="content_sha256"):
        S.load_split_manifest(path)
    doubled = copy.deepcopy(first)
    doubled["frame_ids"]["valid"].append(doubled["frame_ids"]["train"][0])
    with pytest.raises(LeakageError):
        S.check_split_disjoint(doubled)


# --------------------------------------------------------------------------
# writing split files (needs ASE for the extxyz export)
# --------------------------------------------------------------------------

def _export_dir(tmp_path, spec, *, name="export"):
    """A SYNTHETIC export directory: dataset.extxyz (+ a dataset_manifest.json stand-in)."""
    pytest.importorskip("ase")
    from nio_md_prep.dataset import extxyz as X

    export = tmp_path / name
    export.mkdir()
    payloads = []
    rng = np.random.default_rng(0)
    cell = np.diag([6.0, 6.5, 7.0])
    for group, runs in spec.items():
        for run, n_frames in runs.items():
            for index in range(n_frames):
                species = ["Ni", "O", "Ni", "O"] if group != "g05" else ["Ni", "O", "O", "H"]
                payloads.append(X.FramePayload(
                    frame_id=f"synth:{run}#{index:05d}", species=species, cell=cell,
                    positions=rng.random((4, 3)) @ cell, energy=float(-20.0 - rng.random()),
                    forces=rng.normal(size=(4, 3)), stress=None,
                    info={
                        "generator": SYNTHETIC, "source": f"synth:{run}/vasprun.xml", "source_file_type": "vasprun",
                        "ionic_step": index, "structure_key": f"{run}#{index}", "lineage_group": group,
                        "run_id": f"synth:{run}", "energy_source": "vasprun:calculation.e_fr_energy-PSTRESS*V",
                        "energy_rule": "calc_level_direct", "stress_available": False, "stress_reason": "not_computed",
                        "scf_status": "converged", "label_source": "dft", "vasp_version": "6.4.2",
                        "magnetic_class": "controlled_consistent", "magnetic_policy": "accepted",
                        "campaign": "synthetic", "family": "n02" if int(group[1:]) % 2 else "n04",
                        "temperature_K": 300, "parser": "nio-md-prep.vasprun-stream", "parser_version": "1.0",
                        "repo_commit": "0" * 40,
                    },
                ))
    X.write_extxyz(export / "dataset.extxyz", payloads)
    (export / "dataset_manifest.json").write_text('{"synthetic": true}\n', encoding="utf-8")
    return export, payloads


def test_write_split_files_are_deterministic_carry_the_split_key_and_hashes(tmp_path):
    import ase.io

    from nio_md_prep.dataset.fsio import sha256_file

    spec = {f"g{i:02d}": {f"r{i:02d}a": 3, f"r{i:02d}b": 2} for i in range(8)}
    export, payloads = _export_dir(tmp_path, spec)
    first = S.split_export(export, tmp_path / "split1", seed=5)
    second = S.split_export(export, tmp_path / "split2", seed=5)
    for name in ("train.extxyz", "valid.extxyz", "test.extxyz", "split_manifest.json", "split_summary.md"):
        assert (tmp_path / "split1" / name).read_bytes() == (tmp_path / "split2" / name).read_bytes(), name
    assert first == second
    loaded = S.load_split_manifest(tmp_path / "split1" / "split_manifest.json")
    assert loaded["content_sha256"] == first["content_sha256"]
    assert loaded["dataset_extxyz_sha256"] == sha256_file(export / "dataset.extxyz")
    assert loaded["export_manifest_sha256"] == sha256_file(export / "dataset_manifest.json")
    total = 0
    by_id = {payload.frame_id: payload for payload in payloads}
    for split in S.SPLITS:
        path = tmp_path / "split1" / f"{split}.extxyz"
        assert loaded["files"][split]["sha256"] == sha256_file(path)
        frames = ase.io.read(path, index=":", format="extxyz") if path.stat().st_size else []
        assert len(frames) == loaded["files"][split]["frames"] == len(loaded["frame_ids"][split])
        total += len(frames)
        for atoms in frames:
            assert atoms.info["split"] == split and atoms.calc is None and atoms.constraints == []
            assert loaded["group_to_split"][atoms.info["lineage_group"]] == split
            payload = by_id[atoms.info["frame_id"]]
            assert np.array_equal(atoms.arrays["REF_forces"], payload.forces)
            assert atoms.info["REF_energy"] == payload.energy
    assert total == len(payloads)
    assert loaded["statistics"]["elements"].keys() == {"H-Ni-O", "Ni-O"}
    assert loaded["statistics"]["composition"].keys() == {"H1Ni1O2", "Ni2O2"}
    summary = (tmp_path / "split1" / "split_summary.md").read_text(encoding="utf-8")
    assert "lineage_group" in summary and first["content_sha256"] in summary
    assert loaded["summary_file"]["sha256"] == sha256_file(tmp_path / "split1" / "split_summary.md")


def test_write_split_refuses_overwrite_mismatch_and_self_output(tmp_path):
    spec = {f"g{i:02d}": {f"r{i:02d}": 3} for i in range(5)}
    export, _ = _export_dir(tmp_path, spec)
    S.split_export(export, tmp_path / "out")
    with pytest.raises(OverwriteRefusedError):
        S.split_export(export, tmp_path / "out")
    S.split_export(export, tmp_path / "out", force=True)
    with pytest.raises(DatasetError, match="own directory"):
        S.split_export(export, export, force=True)
    manifest = S.load_split_manifest(tmp_path / "out" / "split_manifest.json")
    partial = copy.deepcopy(manifest)
    partial["frame_ids"]["train"] = partial["frame_ids"]["train"][1:]
    with pytest.raises(DatasetError, match="in no split"):
        S.write_split(export, partial, tmp_path / "out2")
    other = copy.deepcopy(manifest)
    other["export_manifest_sha256"] = "f" * 64
    with pytest.raises(DatasetError, match="export manifest"):
        S.write_split(export, other, tmp_path / "out3")
    with pytest.raises(DatasetError, match="already carries a split key"):
        S.records_from_extxyz(tmp_path / "out" / "train.extxyz")


def test_split_export_append_keeps_test_frozen(tmp_path):
    spec = {f"g{i:02d}": {f"r{i:02d}": 3} for i in range(8)}
    export, _ = _export_dir(tmp_path, spec, name="v1")
    first = S.split_export(export, tmp_path / "s1", seed=2)
    grown = dict(spec, **{f"g{i:02d}": {f"r{i:02d}": 3} for i in range(8, 14)})
    export2, _ = _export_dir(tmp_path, grown, name="v2")
    second = S.split_export(export2, tmp_path / "s2", seed=2, previous=tmp_path / "s1" / "split_manifest.json")
    assert second["frame_ids"]["test"] == first["frame_ids"]["test"]
    assert (tmp_path / "s1" / "test.extxyz").read_bytes() == (tmp_path / "s2" / "test.extxyz").read_bytes()
