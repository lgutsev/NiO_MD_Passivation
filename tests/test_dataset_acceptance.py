"""Run- and frame-level acceptance policy (labels, SCF gate, MLFF, stress, magnetism).

Every VASP file used here is SYNTHETIC: written by ``tests/dataset_fixtures.py``
(labelled "SYNTHETIC TEST FIXTURE" in-file) and read back with the real
parsers, or built directly as model objects. Nothing here is real NiO data and
nothing here validates real NiO forces.
"""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest

from dataset_fixtures import make_frames, rocksalt_nio, write_vasp_run
from nio_md_prep.dataset import acceptance as ac
from nio_md_prep.dataset import vaspfiles as vf
from nio_md_prep.dataset.errors import DatasetError
from nio_md_prep.dataset.model import (
    ACCEPTED,
    EXCLUDED,
    FrameRecord,
    IncarEvidence,
    MAG_CONTROLLED_CONSISTENT,
    MAG_CONTROLLED_TRANSITION,
    MAG_NON_SPIN_POLARIZED,
    MAG_UNCONTROLLED,
    MAG_UNKNOWN,
    MISSING_LABELS,
    NON_DFT_MLFF_STEP,
    OUTCOMES,
    PARSE_ERROR,
    QUARANTINED,
    REJECTED,
    EV_PER_A3_IN_KBAR,
)
from nio_md_prep.dataset.vasprun import read_vasprun

MD = {"IBRION": 0, "NSW": 4, "POTIM": 1.0}
AFM_MAGMOM = "2 -2 2 -2 0 0 0 0"
AFM_SITES = [1.7, -1.7, 1.7, -1.7, 0.0, 0.0, 0.0, 0.0]
FM_SITES = [1.7, 1.7, 1.7, 1.7, 0.0, 0.0, 0.0, 0.0]


def nio(n=4, **fields):
    species, cell, positions = rocksalt_nio()
    return species, cell, make_frames(species, cell, positions, n, **fields)


def parse_run(out, *, drop=()):
    """Read a written SYNTHETIC run with the real parsers: (header, steps, trailer, evidence)."""
    header, steps, trailer = read_vasprun(out["vasprun"])
    evidence = {
        "outcar": vf.parse_outcar(out["outcar"]) if out["outcar"] else None,
        "oszicar": vf.parse_oszicar(out["oszicar"]) if out["oszicar"] else None,
        "potcar": vf.parse_potcar(out["potcar"]) if out["potcar"] else None,
        "poscar": vf.parse_poscar(out["poscar"]),
        "incar": vf.parse_incar(out["incar"]),
    }
    for name in drop:
        evidence[name] = None
    return header, steps, trailer, evidence


def assess(tmp_path, name="run", *, species=None, cell=None, frames=None, policy=None, drop=(),
           reference_hashes=(), mutate=None, **writer):
    """Write -> parse -> classify_run -> classify_steps (one SYNTHETIC run)."""
    if frames is None:
        species, cell, frames = nio()
    out = write_vasp_run(tmp_path / name, species=species, cell=cell, frames=frames, **writer)
    header, steps, trailer, evidence = parse_run(out, drop=drop)
    if mutate is not None:
        mutate(header, steps, trailer, evidence)
    run = ac.classify_run(run_id=f"t:{name}", header=header, trailer=trailer, steps=steps, policy=policy,
                          reference_hashes=reference_hashes, **evidence)
    return run, ac.classify_steps(steps, run), out


def statuses(steps):
    return [(s.outcome.status, s.outcome.reason) for s in steps]


# --------------------------------------------------------------------------
# accepted baseline, energy provenance, accounting
# --------------------------------------------------------------------------

def test_clean_md_run_is_accepted_with_energy_provenance(tmp_path):
    run, steps, out = assess(tmp_path, incar=MD)
    assert run.outcome.status == ACCEPTED and run.calc_type == "md" and run.findings == []
    assert statuses(steps) == [(ACCEPTED, None)] * 4
    expected = out["expected"]["steps"]
    for step, exp in zip(steps, expected):
        assert step.label_set == ac.LABEL_SET_ENERGY_FORCES and step.label_source == "dft"
        assert step.label_energy == pytest.approx(exp["free_energy"], abs=1e-12)
        assert step.energy["energy_source"].startswith("vasprun:calculation.e_fr_energy")
        assert step.energy["energy_rule"] == "calc_level_direct"
        assert step.energy["parser"] and step.energy["parser_version"]
        assert step.scf["status"] == "converged" and step.scf["evidence"] == "scstep_count+outcar_marker"
        assert step.time_fs == pytest.approx((step.index + 1) * 1.0)
    record = steps[0].record_fields()
    assert record["label_energy"] == steps[0].label_energy and record["outcome"].status == ACCEPTED
    assert set(record["energies"]) == set(ac.ENERGY_QUANTITIES)
    frame = FrameRecord(frame_id=f"{run.run_id}#00000", run_id=run.run_id, **record)  # fields map 1:1
    assert frame.stress_available is False and frame.stress_reason == "not_requested"
    for audit in (run.as_dict(), steps[0].as_dict()):  # audit records are JSON-safe
        assert json.loads(json.dumps(audit)) == audit


def test_legacy_vasp_energies_are_reconstructed_and_labelled(tmp_path):
    run, steps, out = assess(tmp_path, incar=MD, version="5.4.4")
    exp = out["expected"]["steps"][1]
    step = steps[1]
    assert step.outcome.status == ACCEPTED
    assert step.energy["energy_rule"] == "reconstructed_last_scstep:vasp<=6.0.8"
    assert step.energy["energy_version_gate"] == "vasp<=6.0.8"
    assert step.energy["energy_sigma0"] == pytest.approx(exp["energy_sigma0"], abs=1e-12)
    assert step.energy["energy_no_entropy"] == pytest.approx(exp["energy_no_entropy"], abs=1e-12)
    assert "calc_level_mislabel_pattern" in step.flags  # the known legacy layout, recorded not hidden


def test_every_step_gets_exactly_one_outcome(tmp_path):
    species, cell, frames = nio(5)
    frames[1]["omit"] = {"forces"}
    frames[2].update(scf_marker="not_reached")
    frames[3].update(n_scf=20, last_dE=1e-3, scf_marker=None)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=5),
                           mlff_steps=[4])
    assert len(steps) == run.n_steps == 5
    counts = ac.outcome_counts(steps)
    assert sum(n for bucket in counts.values() for n in bucket.values()) == 5
    assert set(counts) <= set(OUTCOMES)
    assert counts[NON_DFT_MLFF_STEP] == {"mlff_step": 1}


# --------------------------------------------------------------------------
# missing labels / energy-only
# --------------------------------------------------------------------------

def test_missing_forces_is_missing_labels_by_default(tmp_path):
    species, cell, frames = nio(3)
    frames[1]["omit"] = {"forces"}
    _, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=3))
    assert statuses(steps)[1] == (MISSING_LABELS, "missing_forces")
    assert steps[1].label_set is None
    assert [s.outcome.status for s in steps[::2]] == [ACCEPTED, ACCEPTED]


def test_energy_only_needs_the_explicit_flag_and_is_tagged(tmp_path):
    species, cell, frames = nio(3)
    frames[1]["omit"] = {"forces"}
    policy = ac.Policy(allow_energy_only=True)
    _, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=3), policy=policy)
    assert steps[1].outcome.status == ACCEPTED
    assert steps[1].label_set == ac.LABEL_SET_ENERGY_ONLY == "energy_only"
    assert "energy_only" in steps[1].flags and steps[1].max_force is None
    assert steps[1].as_dict()["forces_available"] is False and "forces" in steps[1].as_dict()["forces_reason"]
    assert steps[0].as_dict()["forces_available"] is True and steps[0].as_dict()["forces_reason"] is None
    assert steps[0].label_set == ac.LABEL_SET_ENERGY_FORCES
    assert policy.non_default() == {"allow_energy_only": True}


def test_missing_energy_is_missing_labels(tmp_path):
    species, cell, frames = nio(3)
    frames[2]["omit"] = {"energy"}
    _, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=3))
    assert statuses(steps)[2] == (MISSING_LABELS, "missing_energy")


def test_non_finite_forces_are_rejected(tmp_path):
    species, cell, frames = nio(3)
    frames[1]["overflow"] = {"forces"}
    _, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=3))
    assert statuses(steps)[1] == (REJECTED, "non_finite")


# --------------------------------------------------------------------------
# SCF gate: failed -> rejected, uncertain -> quarantined
# --------------------------------------------------------------------------

def test_scf_explicit_failures_are_rejected(tmp_path):
    species, cell, frames = nio(4)
    frames[1].update(scf_marker="not_reached")
    frames[2].update(n_scf=60, last_dE=1e-3, scf_marker=None)  # NELM ceiling, |dE| >= EDIFF
    frames[3].update(scf_marker="hard_stop")
    _, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=MD)
    assert statuses(steps) == [(ACCEPTED, None)] + [(REJECTED, "scf_not_converged")] * 3
    assert steps[1].scf["evidence"] == "outcar_marker" and steps[1].scf["converged"] is False
    assert steps[2].scf["evidence"] == "scstep_count" and "NELM" in steps[2].outcome.detail


def test_scf_uncertain_cases_are_quarantined(tmp_path):
    species, cell, frames = nio(4)
    frames[1].update(n_scf=20, last_dE=1e-3, scf_marker=None)  # left the loop early without meeting EDIFF
    frames[2].update(scf_marker=None)  # VASP 6 OUTCAR without the 'EDIFF is reached' marker
    frames[3].update(n_scf=60, last_dE=1e-6, scf_marker=None)  # at the NELM ceiling: ambiguous
    _, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=MD)
    assert statuses(steps) == [(ACCEPTED, None)] + [(QUARANTINED, "scf_convergence_unknown")] * 3
    assert all(s.scf["converged"] is None for s in steps[1:])


def test_scf_nelm_ceiling_is_accepted_only_with_a_vasp6_marker(tmp_path):
    species, cell, frames = nio(2)
    frames[0].update(n_scf=60, last_dE=1e-6)  # marker 'reached'
    _, steps6, _ = assess(tmp_path, "v6", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2))
    _, steps5, _ = assess(tmp_path, "v5", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                          version="5.4.4")
    assert steps6[0].outcome.status == ACCEPTED and steps6[0].scf["marker_is_proof"] is True
    assert statuses(steps5)[0] == (QUARANTINED, "scf_convergence_unknown")  # VASP 5 prints it even at NELM
    assert steps5[1].outcome.status == ACCEPTED


def test_scf_single_scstep_needs_the_marker(tmp_path):
    species, cell, frames = nio(2)
    frames[0].update(n_scf=1)
    _, with_outcar, _ = assess(tmp_path, "a", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2))
    _, without, _ = assess(tmp_path, "b", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                           write_outcar=False)
    assert with_outcar[0].outcome.status == ACCEPTED and with_outcar[0].scf["evidence"] == "outcar_marker"
    assert statuses(without)[0] == (QUARANTINED, "scf_convergence_unknown")
    assert without[1].outcome.status == ACCEPTED and without[1].scf["evidence"] == "scstep_count"


def test_scf_nwrite1_marker_exemption_is_flagged(tmp_path):
    species, cell, frames = nio(3)
    for frame in frames[1:]:
        frame["scf_marker"] = None
    _, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=3, NWRITE=1))
    assert [s.outcome.status for s in steps] == [ACCEPTED] * 3
    assert "scf_marker_not_printed" in steps[1].flags and "scf_marker_not_printed" not in steps[0].flags
    policy = ac.Policy(nwrite_marker_exemption=False)
    _, strict, _ = assess(tmp_path, "strict", species=species, cell=cell, frames=frames,
                          incar=dict(MD, NSW=3, NWRITE=1), policy=policy)
    assert statuses(strict)[1:] == [(QUARANTINED, "scf_convergence_unknown")] * 2


def test_scf_ediff_zero_is_quarantined_unless_allowed(tmp_path):
    species, cell, frames = nio(2)
    run, steps, _ = assess(tmp_path, "a", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                           ediff=0.0)
    assert "scf_criterion_disabled" in run.flags
    assert statuses(steps) == [(QUARANTINED, "scf_convergence_unknown")] * 2
    _, allowed, _ = assess(tmp_path, "b", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                           ediff=0.0, policy=ac.Policy(allow_ediff_zero=True))
    assert [s.outcome.status for s in allowed] == [ACCEPTED] * 2


def test_scf_unmappable_vasp6_outcar_is_quarantined(tmp_path):
    def drop_last_outcar_step(header, steps, trailer, evidence):
        outcar = evidence["outcar"]
        evidence["outcar"] = dataclasses.replace(
            outcar, free_energies=outcar.free_energies[:-1], scf_marker_kinds=outcar.scf_marker_kinds[:-1])

    run, steps, _ = assess(tmp_path, incar=MD, mutate=drop_last_outcar_step)
    assert run.outcar_map is None and "evidence_count_mismatch:outcar" in run.flags
    assert {s.outcome.reason for s in steps} == {"scf_convergence_unknown"}


# --------------------------------------------------------------------------
# VASP-MLFF
# --------------------------------------------------------------------------

def test_mlff_flat_steps_are_counted_but_never_labels(tmp_path):
    species, cell, frames = nio(3)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, NSW=3, ML_LMLFF=True, ML_MODE="train"), mlff_steps=[1])
    assert run.outcome.status == ACCEPTED and run.mlff["mode"] == "train" and run.mlff["n_mlff_steps"] == 1
    assert [s.index for s in steps] == [0, 1, 2]  # indexing keeps the MLFF step
    assert statuses(steps) == [(ACCEPTED, None), (NON_DFT_MLFF_STEP, "mlff_step"), (ACCEPTED, None)]
    mlff = steps[1]
    assert mlff.label_set is None and mlff.label_energy is None and mlff.label_source == "mlff"
    assert mlff.stress.available is False and mlff.stress.reason == "not_a_dft_label"


def test_mlff_prediction_runs_have_no_dft_frames(tmp_path):
    species, cell, frames = nio(3)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, NSW=3, ML_LMLFF=True, ML_MODE="run"), mlff_steps=[0, 1, 2])
    assert (run.outcome.status, run.outcome.reason) == (NON_DFT_MLFF_STEP, "mlff_steps")
    assert {s.outcome.status for s in steps} == {NON_DFT_MLFF_STEP}


@pytest.mark.parametrize("facts, expected", [
    ({}, (False, None)),
    ({"ml_mode": "none"}, (False, None)),
    ({"ml_lmlff": True, "ml_mode": "train"}, (True, "train")),
    ({"ml_lmlff": True, "ml_istart": 2}, (True, "run")),
    ({"ml_lmlff": True, "ml_istart": 0}, (True, "train")),
    ({"ml_lmlff": True}, (True, "unknown")),
    ({"ml_mode": "weird"}, (True, "unknown")),
])
def test_mlff_mode(facts, expected):
    assert ac.mlff_mode(facts) == expected


def test_mlff_unknown_mode_rejects_the_whole_run(tmp_path):
    species, cell, frames = nio(2)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, NSW=2, ML_LMLFF=True))
    assert run.outcome.reason == "mlff_steps" and "undeterminable" in run.outcome.detail
    assert {s.outcome.reason for s in steps} == {"mlff_steps"}


# --------------------------------------------------------------------------
# stress: omitted with a reason, never zeros
# --------------------------------------------------------------------------

def test_stress_is_not_exported_unless_requested(tmp_path):
    _, steps, _ = assess(tmp_path, incar=MD)
    stress = steps[0].stress
    assert (stress.available, stress.reason, stress.stress) == (False, "not_requested", None)
    assert stress.raw_present is True and stress.source == ac.STRESS_SOURCE


def test_requested_stress_is_converted_to_ase_sign(tmp_path):
    run, steps, out = assess(tmp_path, incar=MD, policy=ac.Policy(include_stress=True))
    raw = out["expected"]["steps"][0]["stress_kbar"]
    stress = steps[0].stress
    assert stress.available and stress.reason is None
    np.testing.assert_allclose(stress.stress, -0.5 * (raw + raw.T) / EV_PER_A3_IN_KBAR, rtol=0, atol=1e-15)
    assert steps[0].as_dict()["stress_available"] is True


def test_isif0_stress_is_absent_not_zero(tmp_path):
    _, steps, _ = assess(tmp_path, incar=dict(MD, ISIF=0), policy=ac.Policy(include_stress=True))
    for step in steps:
        assert (step.stress.available, step.stress.reason, step.stress.stress) == (False, "not_computed", None)
        assert step.stress.raw_present is False
        assert step.record_fields()["stress_available"] is False


@pytest.mark.parametrize("incar, policy, reason", [
    ({"ISIF": 1}, {}, "isif1_trace_only"),
    ({}, {"allow_pstress": False}, "pstress_nonzero"),
    ({}, {"allow_pstress": True}, None),
])
def test_stress_reasons(tmp_path, incar, policy, reason):
    _, steps, _ = assess(tmp_path, incar=dict(MD, **incar), pstress=5.0 if "allow_pstress" in policy else 0.0,
                         policy=ac.Policy(include_stress=True, **policy))
    assert steps[0].stress.reason == reason
    assert steps[0].stress.available is (reason is None)
    if reason is not None:
        assert steps[0].stress.stress is None


def test_vacuum_slab_stress_is_omitted_unless_allowed(tmp_path):
    species, _, positions = rocksalt_nio()
    cell = np.diag([4.17, 4.17, 20.0])
    frames = make_frames(species, cell, positions, 2)
    run, steps, _ = assess(tmp_path, "slab", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                           policy=ac.Policy(include_stress=True))
    assert steps[0].vacuum_axes == [2] and steps[0].vacuum_gaps[2] > 15.0
    assert steps[0].stress.reason == "vacuum_slab" and steps[0].stress.stress is None
    _, allowed, _ = assess(tmp_path, "slab2", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                           policy=ac.Policy(include_stress=True, allow_vacuum_stress=True))
    assert allowed[0].stress.available


def test_vacuum_axes_are_measured_along_plane_normals():
    cell = np.array([[4.0, 0.0, 0.0], [2.0, 3.4641016, 0.0], [0.0, 0.0, 30.0]])  # hexagonal, 60 degrees
    frac = np.array([[0.0, 0.0, 0.1], [0.5, 0.5, 0.2]])
    axes, gaps = ac.vacuum_axes(cell, frac, threshold=5.0)
    assert axes == [2] and gaps[2] == pytest.approx(0.9 * 30.0)
    assert gaps[0] == pytest.approx(0.5 * ac.interplanar_spacings(cell)[0])


# --------------------------------------------------------------------------
# magnetism: classes, default NiO quarantine, explicit overrides
# --------------------------------------------------------------------------

def test_non_spin_polarized_run(tmp_path):
    run, steps, _ = assess(tmp_path, incar=MD)
    assert run.magnetic.run["magnetic_class"] == MAG_NON_SPIN_POLARIZED
    assert run.magnetic.run["magnetic_species"] == []  # ISPIN=1: no magnetic degree of freedom
    assert all(s.magnetic["policy_result"] == "accepted" for s in steps)


def test_controlled_consistent_afm_run(tmp_path):
    species, cell, frames = nio(3, mag_total=0.0, site_moments=AFM_SITES)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, NSW=3, ISPIN=2, MAGMOM=AFM_MAGMOM))
    mag = run.magnetic.run
    assert mag["magnetic_class"] == MAG_CONTROLLED_CONSISTENT and mag["magmom_explicit"] is True
    assert mag["magmom_initial"] == [2.0, -2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    assert mag["magmom_initial_pattern"] == [1, -1, 1, -1] and mag["magnetic_species"] == ["Ni"]
    assert run.outcome.status == ACCEPTED and mag["policy_result"] == "accepted"
    frame = steps[1].magnetic
    assert frame["total"] == pytest.approx(0.0) and frame["total_source"] == "outcar"
    assert frame["site_moments"] == pytest.approx(AFM_SITES) and frame["site_source"] == "outcar_table_step"
    assert frame["d_total_first"] == pytest.approx(0.0) and frame["d_total_previous"] == pytest.approx(0.0)
    assert frame["pattern_matches_initial"] is True and frame["segment"] == 0
    assert len({s.magnetic["magnetic_state_id"] for s in steps}) == 1
    assert [s.outcome.status for s in steps] == [ACCEPTED] * 3
    record = steps[1].record_fields()
    assert record["total_magnetization"] == pytest.approx(0.0) and record["magnetic_segment"] == 0


def _afm_to_fm(n=4, flip_at=2):
    species, cell, frames = nio(n)
    for k, frame in enumerate(frames):
        fm = k >= flip_at
        frame.update(mag_total=6.8 if fm else 0.0, site_moments=FM_SITES if fm else AFM_SITES)
    return species, cell, frames


def test_controlled_transition_quarantines_frames_from_the_first_transition(tmp_path):
    species, cell, frames = _afm_to_fm()
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, ISPIN=2, MAGMOM=AFM_MAGMOM))
    mag = run.magnetic.run
    assert mag["magnetic_class"] == MAG_CONTROLLED_TRANSITION and mag["first_transition_index"] == 2
    assert mag["policy_result"] == "partially_quarantined" and run.outcome.status == ACCEPTED
    assert statuses(steps) == [(ACCEPTED, None)] * 2 + [(QUARANTINED, "magnetic_order_changed")] * 2
    assert [s.magnetic["segment"] for s in steps] == [0, 0, 1, 1]
    assert set(steps[2].magnetic["transition"]) == {"total_moment_jump", "site_pattern_change"}
    assert steps[2].magnetic["d_total_previous"] == pytest.approx(6.8)
    assert steps[3].magnetic["d_total_first"] == pytest.approx(6.8)
    assert steps[2].magnetic["pattern_matches_initial"] is False


def test_controlled_transition_override_is_explicit_and_recorded(tmp_path):
    species, cell, frames = _afm_to_fm()
    policy = ac.Policy(accept_transitions=True, magnetic_override_reason="SYNTHETIC test: reviewed FM branch")
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, ISPIN=2, MAGMOM=AFM_MAGMOM), policy=policy)
    assert [s.outcome.status for s in steps] == [ACCEPTED] * 4
    assert steps[2].magnetic["policy_result"] == "accepted_by_override:accept_transitions"
    assert steps[0].magnetic["policy_result"] == "accepted"
    assert run.magnetic.run["overrides"] == ["accept_transitions"]
    assert run.magnetic.run["override_reason"].startswith("SYNTHETIC")
    recorded = policy.as_dict()
    assert recorded["non_default"]["accept_transitions"] is True and recorded["schema"] == ac.ACCEPTANCE_SCHEMA


def test_total_moment_jump_without_site_table_is_a_state_change(tmp_path):
    species, cell, frames = nio(3)
    for frame, total in zip(frames, (0.0, 0.0, 3.0)):
        frame["mag_total"] = total
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, NSW=3, ISPIN=2, MAGMOM=AFM_MAGMOM))
    assert "magnetic_state_unresolved" in run.flags  # totals only, no LORBIT table
    assert statuses(steps)[2] == (QUARANTINED, "magnetic_state_change")
    assert steps[2].magnetic["transition"] == ["total_moment_jump"]


def test_initial_order_mismatch_is_a_transition_at_frame_zero(tmp_path):
    species, cell, frames = nio(2, mag_total=6.8, site_moments=FM_SITES)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames,
                           incar=dict(MD, NSW=2, ISPIN=2, MAGMOM=AFM_MAGMOM))
    assert run.magnetic.run["magnetic_class"] == MAG_CONTROLLED_TRANSITION
    assert "initial_order_mismatch" in steps[0].magnetic["transition"]
    assert "afm_total_moment_mismatch" in steps[0].magnetic["transition"]
    assert {s.outcome.reason for s in steps} == {"magnetic_order_changed"}


def test_uncontrolled_nio_run_is_quarantined_by_default(tmp_path):
    species, cell, frames = nio(3, mag_total=4.0)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=3, ISPIN=2))
    assert run.magnetic.run["magnetic_class"] == MAG_UNCONTROLLED and run.magnetic.run["magmom_explicit"] is False
    assert (run.outcome.status, run.outcome.reason) == (QUARANTINED, "magnetic_uncontrolled")
    assert statuses(steps) == [(QUARANTINED, "magnetic_uncontrolled")] * 3


def test_uncontrolled_override(tmp_path):
    species, cell, frames = nio(3, mag_total=4.0)
    policy = ac.Policy(accept_uncontrolled=True, magnetic_override_reason="SYNTHETIC: FM start is intended")
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=3, ISPIN=2),
                           policy=policy)
    assert run.outcome.status == ACCEPTED
    assert run.magnetic.run["policy_result"] == "accepted_by_override:accept_uncontrolled"
    assert [s.outcome.status for s in steps] == [ACCEPTED] * 3
    assert steps[0].magnetic["policy_result"] == "accepted_by_override:accept_uncontrolled"


def test_uncontrolled_run_without_magnetic_species_is_not_quarantined(tmp_path):
    species = ["O", "H", "H"]
    cell = np.eye(3) * 8.0
    positions = np.array([[4.0, 4.0, 4.0], [4.76, 4.59, 4.0], [3.24, 4.59, 4.0]])
    frames = make_frames(species, cell, positions, 2, mag_total=0.0, stress=False)
    run, steps, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2, ISPIN=2))
    assert run.magnetic.run["magnetic_class"] == MAG_UNCONTROLLED and run.magnetic.run["magnetic_species"] == []
    assert run.outcome.status == ACCEPTED and [s.outcome.status for s in steps] == [ACCEPTED] * 2


def test_unknown_magnetic_state_is_quarantined_unless_overridden(tmp_path):
    species, cell, frames = nio(2, mag_total=0.0)
    kwargs = dict(species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2, ISPIN=2, MAGMOM=AFM_MAGMOM),
                  write_outcar=False, write_oszicar=False)
    run, steps, _ = assess(tmp_path, "a", **kwargs)
    assert run.magnetic.run["magnetic_class"] == MAG_UNKNOWN
    assert (run.outcome.status, run.outcome.reason) == (QUARANTINED, "magnetic_unknown")
    assert statuses(steps) == [(QUARANTINED, "magnetic_unknown")] * 2
    run, steps, _ = assess(tmp_path, "b", policy=ac.Policy(accept_unknown=True), **kwargs)
    assert run.outcome.status == ACCEPTED and [s.outcome.status for s in steps] == [ACCEPTED] * 2


def _classify(**overrides):
    args = dict(species=["Ni", "Ni", "O", "O"], ispin=2, nupdown=None, noncollinear=False,
                magmom_initial=[2.0, -2.0, 0.0, 0.0], magmom_explicit=True, dft_indices=[0, 1],
                totals={0: (0.0, "outcar"), 1: (0.0, "outcar")}, sites={})
    args.update(overrides)
    return ac.classify_magnetism(**args)


def test_magnetic_class_table():
    assert _classify(ispin=1).run["magnetic_class"] == MAG_NON_SPIN_POLARIZED
    assert _classify(ispin=None).run["magnetic_class"] == MAG_UNKNOWN
    assert _classify().run["magnetic_class"] == MAG_CONTROLLED_CONSISTENT
    assert _classify(magmom_explicit=False).run["magnetic_class"] == MAG_UNCONTROLLED
    assert _classify(magmom_explicit=False, nupdown=0.0).run["magnetic_class"] == MAG_CONTROLLED_CONSISTENT
    unknown = _classify(magmom_explicit=None)  # no <incar> and no INCAR: cannot tell whether MAGMOM was set
    assert unknown.run["magnetic_class"] == MAG_UNKNOWN and "unknown whether MAGMOM" in unknown.run["class_reason"]
    assert unknown.run["policy_result"] == "quarantined:magnetic_unknown"
    no_evidence = _classify(totals={})
    assert no_evidence.run["magnetic_class"] == MAG_UNKNOWN
    jump = _classify(totals={0: (0.0, "outcar"), 1: (2.0, "outcar")})
    assert jump.run["magnetic_class"] == MAG_CONTROLLED_TRANSITION
    assert jump.frames[1]["policy_reason"] == "magnetic_state_change" and jump.frames[0]["policy_reason"] is None


def test_frame_without_evidence_in_a_controlled_run_is_quarantined():
    assessment = _classify(dft_indices=[0, 1, 2])  # frame 2 has no total and no site moments
    assert assessment.run["magnetic_class"] == MAG_CONTROLLED_CONSISTENT
    assert assessment.frames[2]["policy_reason"] == "magnetic_unknown"
    assert assessment.run["policy_result"] == "partially_quarantined"


def test_magnetic_species_and_sign_patterns():
    policy = ac.Policy()
    assert ac.magnetic_species_of(["Ni", "O", "H"], None, None, policy) == ["Ni"]
    assert ac.magnetic_species_of(["O", "H"], [1.0, 0.0], True, policy) == ["O"]  # explicit moment >= 0.5
    assert ac.magnetic_species_of(["O", "H"], [1.0, 0.0], False, policy) == []
    assert ac.sign_pattern([1.7, -1.7, 0.1], [0, 1, 2], 0.5) == (1, -1, 0)
    assert ac.sign_pattern([-1.7, 1.7, 0.1], [0, 1, 2], 0.5) == (1, -1, 0)  # global flip is the same state


def test_magmom_explicit_sources():
    assert ac._magmom_explicit({"incar": {"MAGMOM": [1.0]}}, None)[:2] == (True, "vasprun_incar")
    assert ac._magmom_explicit({"incar": {"ISPIN": 2}}, None)[:2] == (False, "vasprun_incar")
    incar = IncarEvidence("0" * 64, {"MAGMOM": "2*1"})
    assert ac._magmom_explicit({"incar": {}}, incar)[:2] == (True, "incar_file")
    assert ac._magmom_explicit({"incar": {}}, None)[:2] == (None, None)
    explicit, _, flags = ac._magmom_explicit({"incar": {"ISPIN": 2}}, incar)
    assert explicit is False and flags == ["magmom_incar_file_differs"]


# --------------------------------------------------------------------------
# run-level rules
# --------------------------------------------------------------------------

def test_unreadable_header_is_a_parse_error():
    run = ac.classify_run(run_id="t:x", header=None, header_error="ParseError: junk")
    assert (run.outcome.status, run.outcome.reason) == (PARSE_ERROR, "vasprun_unreadable")
    assert "junk" in run.outcome.detail


def test_not_vasp_and_unsupported_calculations_are_rejected(tmp_path):
    def not_vasp(header, *_):
        header.generator["program"] = "notvasp"

    run, steps, _ = assess(tmp_path, "a", incar=MD, mutate=not_vasp)
    assert (run.outcome.status, run.outcome.reason) == (REJECTED, "not_vasp")
    assert {(s.outcome.status, s.outcome.reason) for s in steps} == {(REJECTED, "not_vasp")}
    run, steps, _ = assess(tmp_path, "b", incar={"IBRION": 6, "NSW": 4})
    assert run.calc_type == "unsupported" and run.outcome.reason == "unsupported_calculation"
    assert {s.outcome.reason for s in steps} == {"unsupported_calculation"}
    run, _, _ = assess(tmp_path, "c", incar={"IBRION": 2, "NSW": 4}, policy=ac.Policy(allowed_calc_types=("md",)))
    assert run.calc_type == "relaxation" and run.outcome.reason == "unsupported_calculation"


@pytest.mark.parametrize("facts, calc_type", [
    ({"ibrion": -1, "nsw": 0}, "static"), ({"ibrion": 2, "nsw": 0}, "static"), ({"ibrion": 1, "nsw": 5}, "relaxation"),
    ({"ibrion": 0, "nsw": 5}, "md"), ({"ibrion": 44, "nsw": 5}, "other"), ({"ibrion": None, "nsw": 5}, "other"),
    ({"ibrion": 8, "nsw": 1}, "unsupported"), ({"ibrion": -1, "nsw": 0, "lepsilon": True}, "unsupported"),
])
def test_calc_type_of(facts, calc_type):
    assert ac.calc_type_of(facts)[0] == calc_type


def test_copied_reference_output_is_excluded_and_copied_evidence_ignored(tmp_path):
    species, cell, frames = nio(2)
    out = write_vasp_run(tmp_path / "ref", species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2))
    header, steps, trailer, evidence = parse_run(out)
    run = ac.classify_run(run_id="t:ref", header=header, trailer=trailer, steps=steps,
                          reference_hashes=[trailer.source_sha256.upper()], **evidence)
    assert (run.outcome.status, run.outcome.reason) == (EXCLUDED, "copied_reference_output")
    assert {s.outcome.reason for s in ac.classify_steps(steps, run)} == {"copied_reference_output"}
    run = ac.classify_run(run_id="t:ref", header=header, trailer=trailer, steps=steps,
                          reference_hashes=[evidence["outcar"].sha256], **evidence)
    assert run.outcome.status == ACCEPTED and "copied_reference_evidence_ignored:OUTCAR" in run.flags
    assert run.outcar is None and run.evidence["ignored"] == ["OUTCAR"]


def test_species_mismatch_is_quarantined(tmp_path):
    run, steps, _ = assess(tmp_path, incar=MD, potcar_file_order=["O", "Ni"])
    assert (run.outcome.status, run.outcome.reason) == (QUARANTINED, "species_mismatch")
    assert "POTCAR" in run.outcome.detail
    assert {s.outcome.reason for s in steps} == {"species_mismatch"}


def test_evidence_mismatch_is_quarantined(tmp_path):
    def shift_toten(header, steps, trailer, evidence):
        outcar = evidence["outcar"]
        energies = list(outcar.free_energies)
        energies[1] += 0.01
        evidence["outcar"] = dataclasses.replace(outcar, free_energies=energies)

    run, steps, _ = assess(tmp_path, incar=MD, mutate=shift_toten)
    assert (run.outcome.status, run.outcome.reason) == (QUARANTINED, "evidence_mismatch")
    assert "OUTCAR TOTEN" in run.outcome.detail
    assert {s.outcome.reason for s in steps} == {"evidence_mismatch"}


def test_incomplete_run_exclude_and_recover(tmp_path):
    species, cell, frames = nio(4)
    kwargs = dict(species=species, cell=cell, frames=frames, incar=MD, truncate_after_steps=2,
                  truncate_inside_step=True)
    run, steps, _ = assess(tmp_path, "a", **kwargs)
    assert (run.outcome.status, run.outcome.reason) == (EXCLUDED, "run_incomplete")
    assert statuses(steps) == [(EXCLUDED, "run_incomplete")] * 2 + [(REJECTED, "truncated_frame")]
    run, steps, _ = assess(tmp_path, "b", policy=ac.Policy(interrupted_run_policy="recover"), **kwargs)
    assert run.outcome.status == ACCEPTED and "recovered_incomplete_run" in run.flags
    assert statuses(steps) == [(ACCEPTED, None)] * 2 + [(REJECTED, "truncated_frame")]


def test_force_outlier_is_opt_in(tmp_path):
    _, steps, _ = assess(tmp_path, "a", incar=MD)
    assert all(s.outcome.status == ACCEPTED for s in steps)
    threshold = min(s.max_force for s in steps) / 2
    _, steps, _ = assess(tmp_path, "b", incar=MD, policy=ac.Policy(force_outlier_threshold=threshold))
    assert {s.outcome.reason for s in steps} == {"force_outlier"}


# --------------------------------------------------------------------------
# selective dynamics: explicit direct-basis flags, raw forces
# --------------------------------------------------------------------------

def test_selective_dynamics_is_an_explicit_direct_basis_array(tmp_path):
    species, cell, frames = nio(2)
    flags = np.ones((8, 3), dtype=bool)
    flags[:2] = False
    run, steps, out = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                             selective=flags)
    assert run.outcome.status == ACCEPTED
    np.testing.assert_array_equal(steps[0].selective, flags)
    info = run.as_dict()["selective"]
    assert info == {"present": True, "source": "vasprun:initialpos", "basis": "direct", "n_fixed_components": 6}
    assert steps[0].as_dict()["selective_dynamics"] is True
    # raw forces: fixed atoms keep their DFT forces (nothing is zeroed by the policy)
    assert steps[0].max_force == pytest.approx(
        float(np.max(np.linalg.norm(out["expected"]["steps"][0]["forces"], axis=1))))


def test_selective_flags_disagreeing_with_poscar_are_an_evidence_mismatch(tmp_path):
    species, cell, frames = nio(2)
    flags = np.ones((8, 3), dtype=bool)
    other = flags.copy()
    other[0] = False
    run, _, _ = assess(tmp_path, species=species, cell=cell, frames=frames, incar=dict(MD, NSW=2),
                       selective=flags, poscar_selective=other)
    assert (run.outcome.status, run.outcome.reason) == (QUARANTINED, "evidence_mismatch")
    assert "POSCAR" in run.outcome.detail


# --------------------------------------------------------------------------
# policy object
# --------------------------------------------------------------------------

def test_policy_validation_and_round_trip():
    policy = ac.Policy.from_mapping({"include_stress": True, "magnetic_species": ["Ni", "Co"],
                                     "force_outlier_threshold": 20})
    assert policy.magnetic_species == ("Co", "Ni") and policy.force_outlier_threshold == 20.0
    assert ac.Policy.from_mapping(policy.as_dict()) == policy
    with pytest.raises(DatasetError, match="unknown policy keys"):
        ac.Policy.from_mapping({"include_stres": True})
    with pytest.raises(DatasetError, match="interrupted_run_policy"):
        ac.Policy(interrupted_run_policy="maybe")
    with pytest.raises(DatasetError, match="true or false"):
        ac.Policy(accept_unknown="yes")
    with pytest.raises(DatasetError, match="positive number"):
        ac.Policy(mag_jump_threshold=0)
    with pytest.raises(DatasetError, match="element symbols"):
        ac.Policy(magnetic_species=("nickel",))
    assert ac.Policy().magnetic_overrides == [] and ac.DEFAULT_POLICY.non_default() == {}
