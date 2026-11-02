from pathlib import Path
import json
import math

import numpy as np
import pytest

from nio_md_prep.analysis.accessibility import (
    analyze_accessibility,
    exposure_persistence,
    periodic_opening,
    reference_surface,
    setting_key,
    statistical_inefficiency,
)
from nio_md_prep.analysis.route_comparison import (
    classify_system,
    gap_seeded_context,
    paired_route_differences,
    replica_statistics,
)


# ---------------------------------------------------------------- geometry


def test_opening_removes_slits_narrower_than_the_probe():
    mask = np.zeros((40, 80), dtype=bool)
    mask[:, 10:12] = True  # 0.5 A wide slit at 0.25 A spacing
    mask[:, 40:60] = True  # 5 A wide channel
    opened = periodic_opening(mask, 1.4, 0.25, 0.25)
    assert not opened[:, 10:12].any()
    assert opened[:, 40:60].all()
    assert periodic_opening(mask, 0.0, 0.25, 0.25).sum() == mask.sum()


def test_opening_is_periodic_across_the_cell_seam():
    mask = np.zeros((40, 80), dtype=bool)
    mask[:, :8] = True
    mask[:, -8:] = True  # one 4 A channel split by the x seam
    opened = periodic_opening(mask, 1.4, 0.25, 0.25)
    assert opened[:, :8].all() and opened[:, -8:].all()


def test_reference_surface_is_flat_for_a_flat_layer():
    xs, ys = np.meshgrid(np.arange(0, 10, 2.0), np.arange(0, 10, 2.0))
    bounds = ((0.0, 10.0), (0.0, 10.0), (-5.0, 20.0))
    reference = reference_surface(
        xs.ravel(), ys.ravel(), np.full(xs.size, 1.5), bounds, (40, 40), 2.5
    )
    assert np.allclose(reference, 1.5)


def test_reference_surface_rejects_a_radius_that_leaves_holes():
    bounds = ((0.0, 10.0), (0.0, 10.0), (-5.0, 20.0))
    with pytest.raises(ValueError, match="reference_radius"):
        reference_surface(np.array([5.0]), np.array([5.0]), np.array([0.0]), bounds, (40, 40), 1.0)


# ---------------------------------------------------------------- statistics


def test_statistical_inefficiency_detects_correlation():
    rng = np.random.default_rng(1)
    white = rng.normal(size=4000)
    assert statistical_inefficiency(white) == pytest.approx(1.0, abs=0.3)
    ar = np.zeros(4000)
    for i in range(1, len(ar)):
        ar[i] = 0.9 * ar[i - 1] + rng.normal()
    # analytic g = (1 + phi) / (1 - phi) = 19
    assert 12.0 < statistical_inefficiency(ar) < 26.0


def test_frozen_exposure_is_fully_persistent_and_censored():
    stack = np.zeros((20, 10, 10), dtype=bool)
    stack[:, :, :3] = True
    result, probability = exposure_persistence(stack, [1, 5, 10])
    assert result["persistently_exposed_fraction"] == pytest.approx(0.3)
    assert result["decorrelation_censored"]
    assert result["first_last_jaccard"] == pytest.approx(1.0)
    assert all(row["normalized"] == pytest.approx(1.0) for row in result["survival"])


def test_uncorrelated_exposure_has_no_memory():
    rng = np.random.default_rng(2)
    stack = rng.random((200, 30, 30)) < 0.3
    result, _ = exposure_persistence(stack, [1, 5])
    assert abs(result["survival"][0]["normalized"]) < 0.05
    assert result["persistently_exposed_fraction"] == 0.0
    assert result["decorrelation_lag_frames"] == 1


# ---------------------------------------------------------------- end to end

BOX = 12.0


def _substrate_rows():
    rows = []
    atom_id = 1
    for layer_z in (-2.0, 0.0):
        for i in range(6):
            for j in range(6):
                element_type = 4 if (i + j + (layer_z == 0.0)) % 2 == 0 else 2
                rows.append((atom_id, 0, element_type, 2.0 * i, 2.0 * j, layer_z))
                atom_id += 1
    return rows


def _molecule_rows(first_id, molecule, x, y, z0):
    # P, three O, three-carbon chain pointing up; O1 points down toward Ni
    return [
        (first_id, molecule, 3, x, y, z0 + 1.5),
        (first_id + 1, molecule, 2, x, y, z0),
        (first_id + 2, molecule, 2, x + 1.4, y, z0 + 2.0),
        (first_id + 3, molecule, 2, x - 1.4, y, z0 + 2.0),
        (first_id + 4, molecule, 1, x, y, z0 + 3.0),
        (first_id + 5, molecule, 1, x, y, z0 + 4.5),
        (first_id + 6, molecule, 1, x, y, z0 + 6.0),
    ]


def _bonds(first_id):
    p = first_id
    return [(p, p + 1), (p, p + 2), (p, p + 3), (p, p + 4), (p + 4, p + 5), (p + 5, p + 6)]


def _system(floating_height):
    substrate = _substrate_rows()
    first = len(substrate) + 1
    # Ni sits at even (i + j) on the top layer: (4, 4, 0) is Ni
    anchored = _molecule_rows(first, 1, 4.0, 4.0, 2.0)
    floating = _molecule_rows(first + 7, 2, 10.0, 10.0, floating_height)
    bonds = _bonds(first) + _bonds(first + 7)
    return substrate + anchored + floating, bonds


def _write_build(tmp_path: Path, heights) -> Path:
    build = tmp_path / "build"
    build.mkdir()
    rows, bonds = _system(heights[0])
    atoms = "\n".join(f"{a} {m} {t} 0.0 {x} {y} {z}" for a, m, t, x, y, z in rows)
    bond_lines = "\n".join(f"{i} 1 {l} {r}" for i, (l, r) in enumerate(bonds, 1))
    (build / "topology_output.lmp").write_text(
        f"""LAMMPS data

{len(rows)} atoms
{len(bonds)} bonds
4 atom types
1 bond types

0.0 {BOX} xlo xhi
0.0 {BOX} ylo yhi
-5.0 60.0 zlo zhi

Masses

1 12.011
2 15.999
3 30.974
4 58.6934

Atoms # full

{atoms}

Bonds

{bond_lines}
""",
        encoding="utf-8",
    )
    (build / "assembly_manifest.json").write_text(
        json.dumps({"components": [
            {"component": "primary", "molecule_ids": [1, 1]},
            {"component": "secondary", "molecule_ids": [2, 2]},
        ]}),
        encoding="utf-8",
    )
    frames = []
    for step, height in enumerate(heights):
        rows, _ = _system(height)
        body = "\n".join(f"{a} {m} {t} {x} {y} {z}" for a, m, t, x, y, z in rows)
        frames.append(
            f"ITEM: TIMESTEP\n{step * 1000}\nITEM: NUMBER OF ATOMS\n{len(rows)}\n"
            f"ITEM: BOX BOUNDS pp pp ff\n0 {BOX}\n0 {BOX}\n-5 60\n"
            f"ITEM: ATOMS id mol type x y z\n{body}\n"
        )
    (build / "hold-300K.lammpstrj").write_text("".join(frames), encoding="utf-8")
    return build


def test_canopy_and_height_gated_accessibility_separate_attachment_from_overlayer(tmp_path):
    build = _write_build(tmp_path, [30.0, 30.0, 30.0, 30.0])
    summary_path = analyze_accessibility(
        build, trajectory=build / "hold-300K.lammpstrj", last_frames=None,
        timestep_fs=0.5, write_plots=False,
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    access = summary["accessibility"]
    canopy = access[setting_key(None, 0.0)]["accessible_percent"]["mean"]
    gated = access[setting_key(5.0, 0.0)]["accessible_percent"]["mean"]
    # the floating molecule shadows the surface only in the canopy model
    assert gated > canopy
    # larger probes can never reach more area
    assert (
        access[setting_key(None, 3.0)]["accessible_percent"]["mean"]
        <= access[setting_key(None, 1.4)]["accessible_percent"]["mean"]
        <= canopy
    )
    anchoring = summary["anchoring"]
    assert anchoring["primary"]["bound_group_fraction"] == pytest.approx(1.0)
    assert anchoring["primary"]["denticity_population"]["1"] == pytest.approx(1.0)
    assert anchoring["secondary"]["bound_group_fraction"] == 0.0
    structure = summary["structure"]
    assert structure["primary"]["p_height"]["all"]["mean"] == pytest.approx(3.5, abs=0.01)
    assert structure["primary"]["extension"]["all"]["mean"] == pytest.approx(6.0, abs=0.01)
    assert structure["primary"]["tilt"]["all"]["mean"] < 5.0
    assert structure["secondary"]["p_height"]["all"]["mean"] == pytest.approx(31.5, abs=0.01)
    persistence = summary["persistence"][setting_key(None, 1.4)]
    assert persistence["decorrelation_censored"]  # static fixture never decorrelates
    assert summary["frame_interval_ps"] == pytest.approx(0.5)
    assert "Not a measure of moisture" in summary["interpretation_scope"]
    assert (summary_path.parent / "accessibility_maps.npz").is_file()


def test_detaching_molecule_changes_binding_but_not_its_identity(tmp_path):
    build = _write_build(tmp_path, [30.0, 2.0, 2.0, 2.0])
    summary = json.loads(
        analyze_accessibility(
            build, trajectory=build / "hold-300K.lammpstrj", last_frames=None, write_plots=False
        ).read_text(encoding="utf-8")
    )
    # (10, 10, 0) is Ni on the top layer, so the descended molecule anchors
    secondary = summary["anchoring"]["secondary"]
    assert secondary["bound_group_fraction"] == pytest.approx(0.75)
    assert secondary["ever_bound_fraction"] == 1.0
    assert secondary["persistently_bound_fraction"] == 0.0
    assert secondary["early_half_bound_fraction"] < secondary["late_half_bound_fraction"]


# ---------------------------------------------------------------- route comparison


def test_classify_system_routes_and_inventory():
    assert classify_system("prepared-rerun", "me-4pacz-then-dcz-4p")["route"] == "sequential"
    cosam = classify_system("prepared", "me-4pacz-meo-2pacz-cosam")
    assert cosam["inventory_molecules"] == 287 and cosam["cohort"] == 1
    lego = classify_system("prepared-lego2", "me-4pacz-then-dcz-4p-lego2-seeded")
    assert lego["route"] == "gap-seeded-II" and not lego["independent_replicate"]


def _records():
    base = {"metric": "m", "modifier": "dcz-4p", "temperature_K": 300, "film_state": "relaxed"}
    rows = []
    for cohort, (seq, co) in enumerate([(80.0, 70.0), (85.0, 72.0), (78.0, 75.0)], 1):
        rows.append({**base, "route": "sequential", "cohort": cohort, "value": seq, "independent_replicate": True})
        rows.append({**base, "route": "cosam", "cohort": cohort, "value": co, "independent_replicate": True})
    rows.append({**base, "route": "gap-seeded-I", "cohort": 1, "value": 90.0, "independent_replicate": False})
    return rows


def test_paired_differences_use_cohorts_as_replicates():
    rows = [r for r in paired_route_differences(_records()) if r["modifier"] == "dcz-4p"]
    assert len(rows) == 1
    row = rows[0]
    assert row["n"] == 3
    assert row["mean"] == pytest.approx((10 + 13 + 3) / 3)
    assert row["positive"] == 3
    assert row["verdict"] in {"supported", "consistent sign, CI includes 0"}


def test_gap_seeded_runs_are_contextualized_not_pooled():
    stats = replica_statistics(_records())
    sequential = next(s for s in stats if s["route"] == "sequential")
    assert sequential["independent_runs"] == 3
    context = gap_seeded_context(_records())
    assert context[0]["inside_sequential_range"] is False
    assert context[0]["same_stage1_sequential_cohort1"] == 80.0
    assert math.isfinite(context[0]["sequential_max"])


def test_route_report_writes_tables_and_figures(tmp_path):
    from nio_md_prep.analysis.route_report import write_route_report

    records = []
    for row in _records():
        records.append({**row, "metric": "projected_coverage_percent", "within_run_block_sem": 0.05})
    write_route_report(records, tmp_path / "report")
    for name in ("replica_statistics.csv", "paired_sequential_minus_cosam.csv",
                 "gap_seeded_context.csv", "fig_route_levels_300K_relaxed.png",
                 "fig_paired_differences.png", "fig_uncertainty_budget.png"):
        assert (tmp_path / "report" / name).is_file(), name


def test_accessibility_manifest_records_use_pipe_paths(tmp_path):
    from nio_md_prep.analysis.route_comparison import ACCESSIBILITY_METRICS, load_accessibility_records

    summary = tmp_path / "s.json"
    summary.write_text(json.dumps({
        "frames_analyzed": 500,
        "accessibility": {"canopy_probe1.4A": {"accessible_percent": {"mean": 7.5}}},
    }), encoding="utf-8")
    manifest = tmp_path / "m.csv"
    manifest.write_text(
        "prepared_root,system,stage,summary_json\n"
        f"prepared-rerun,me-4pacz-then-dcz-4p,relax-300K,{summary}\n",
        encoding="utf-8",
    )
    records = load_accessibility_records(manifest, ACCESSIBILITY_METRICS)
    assert len(records) == 1
    assert records[0]["metric"] == "accessible_canopy_probe1.4_percent"
    assert records[0]["value"] == 7.5 and records[0]["cohort"] == 2
