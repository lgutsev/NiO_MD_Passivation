"""Structure keys, exact duplicates, label consistency and near-duplicates.

All structures and labels are SYNTHETIC arrays built in this file.
"""

from __future__ import annotations

import numpy as np
import pytest

from dataset_fixtures import rocksalt_nio
from nio_md_prep.dataset import duplicates as dup
from nio_md_prep.dataset.errors import DatasetError
from nio_md_prep.dataset.model import EXCLUDED, QUARANTINED, StepStructure

SPECIES, CELL, POSITIONS = rocksalt_nio()
PERM = [5, 2, 7, 0, 3, 6, 1, 4]  # a relabelling that keeps nothing in place
FORCES = 0.3 * np.cos(np.arange(24, dtype=float).reshape(8, 3))
AFM = np.array([1.7, -1.7, 1.7, -1.7, 0.0, 0.0, 0.0, 0.0])


def key(species=SPECIES, cell=CELL, positions=POSITIONS, pbc=(True, True, True), **kw):
    return dup.structure_fingerprint(species, cell, pbc, positions, **kw)


# --------------------------------------------------------------------------
# keys
# --------------------------------------------------------------------------

def test_permutation_key_is_invariant_to_atom_order_and_order_key_is_not():
    base = key()
    permuted = key([SPECIES[i] for i in PERM], CELL, POSITIONS[PERM])
    assert permuted.permutation_key == base.permutation_key
    assert permuted.order_key != base.order_key
    # the canonical permutation aligns per-atom arrays of both copies
    np.testing.assert_array_equal(POSITIONS[list(base.permutation)], POSITIONS[PERM][list(permuted.permutation)])


def test_keys_are_invariant_to_periodic_wrapping():
    base = key()
    shifted = POSITIONS.copy()
    shifted[3] += CELL[0] - 2 * CELL[2]  # the same atom, two periodic images away
    assert key(positions=shifted).order_key == base.order_key
    frac = POSITIONS @ np.linalg.inv(CELL)
    near_one = frac.copy()
    near_one[0] = [1.0 - 1e-10, 1.0 - 1e-12, -1e-13]  # within the 1e-9 wrap snap of the atom at the origin
    assert key(fractional=near_one).order_key == key(fractional=frac).order_key
    negative_zero = frac.copy()
    negative_zero[0] = [-0.0, -0.0, -0.0]
    assert key(fractional=negative_zero).order_key == key(fractional=frac).order_key


def test_keys_resolve_1e_6_angstrom_and_see_species_cell_pbc():
    base = key()
    noisy = POSITIONS + 1e-9
    assert key(positions=noisy).order_key == base.order_key
    moved = POSITIONS.copy()
    moved[2, 1] += 1e-4
    assert key(positions=moved).permutation_key != base.permutation_key
    assert key(species=["Co"] + SPECIES[1:]).permutation_key != base.permutation_key
    assert key(cell=CELL * 1.001, positions=POSITIONS * 1.001).permutation_key != base.permutation_key
    assert key(pbc=(True, True, False)).permutation_key != base.permutation_key


def test_non_periodic_axes_are_not_wrapped():
    positions = POSITIONS.copy()
    positions[0, 2] -= CELL[2, 2]
    assert key(positions=positions, pbc=(True, True, False)).order_key != key(pbc=(True, True, False)).order_key
    assert key(positions=positions).order_key == key().order_key


def test_step_structure_key_is_none_for_unkeyable_structures():
    frac = POSITIONS @ np.linalg.inv(CELL)
    step = StepStructure(cell=CELL, fractional=frac, positions=POSITIONS)
    assert dup.step_structure_key(SPECIES, step).order_key == key().order_key
    assert dup.step_structure_key(SPECIES, None) is None
    assert dup.step_structure_key(SPECIES, StepStructure(cell=np.full((3, 3), np.nan), fractional=frac,
                                                         positions=POSITIONS)) is None
    assert dup.step_structure_key(SPECIES[:3], step) is None
    with pytest.raises(DatasetError, match="finite 3x3"):
        dup.structure_keys(SPECIES, np.eye(2), positions=POSITIONS)


# --------------------------------------------------------------------------
# exact duplicates
# --------------------------------------------------------------------------

def candidate(frame_id, *, run_id="r:a", pool="p1", perm=None, energy=-50.0, forces=FORCES, mag=None, sites=None,
              state=None):
    species, positions = SPECIES, POSITIONS
    forces = None if forces is None else np.asarray(forces, dtype=float)
    if perm is not None:
        species, positions = [SPECIES[i] for i in perm], POSITIONS[perm]
        forces = forces[perm] if forces is not None else None
        sites = np.asarray(sites)[perm] if sites is not None else None
    k = key(species, CELL, positions)
    return dup.DuplicateCandidate(frame_id=frame_id, run_id=run_id, pool_id=pool, permutation_key=k.permutation_key,
                                  order_key=k.order_key, permutation=k.permutation, energy=energy, forces=forces,
                                  total_magnetization=mag, site_magmoms=None if sites is None else list(sites),
                                  magnetic_state_id=state)


def test_consistent_duplicates_keep_the_lowest_frame_id_in_any_input_order():
    frames = [candidate("r:b#00003", run_id="r:b", energy=-50.0004), candidate("r:a#00007"),
              candidate("r:a#00002", forces=FORCES + 0.004)]
    for order in (frames, frames[::-1], frames[1:] + frames[:1]):
        report = dup.find_exact_duplicates(order)
        assert report.duplicate_of == {"r:a#00007": "r:a#00002", "r:b#00003": "r:a#00002"}
        assert {k: (o.status, o.reason) for k, o in report.outcomes.items()} == {
            "r:a#00007": (EXCLUDED, "exact_duplicate"), "r:b#00003": (EXCLUDED, "exact_duplicate")}
        assert report.clusters[0]["kept"] == "r:a#00002" and report.clusters[0]["consistent"] is True
    assert report.links == [{"a": "r:a", "b": "r:b", "kind": "exact_duplicate",
                             "key": frames[0].permutation_key}]
    counts = report.as_dict()["counts"]
    assert counts == {"clusters": 1, "exact_duplicate": 2, "contradictory_labels": 0, "links": 1}


def test_permuted_copy_with_permuted_forces_is_a_consistent_duplicate():
    report = dup.find_exact_duplicates([candidate("a#1"), candidate("b#1", run_id="r:b", perm=PERM)])
    assert report.duplicate_of == {"b#1": "a#1"}
    assert report.clusters[0]["comparisons"][0]["max_dF"] == pytest.approx(0.0)


def test_permuted_copy_whose_forces_were_not_permuted_is_contradictory():
    bad = candidate("b#1", perm=PERM)
    bad.forces = FORCES.copy()  # labels in the original atom order: misaligned
    report = dup.find_exact_duplicates([candidate("a#1"), bad])
    assert {o.reason for o in report.outcomes.values()} == {"contradictory_labels"}


@pytest.mark.parametrize("change, violation", [
    ({"energy": -50.01}, "|dE|"),
    ({"forces": FORCES + 0.05}, "max |dF|"),
    ({"forces": None}, "one copy has forces"),
    ({"mag": 2.0}, "mag_total"),
])
def test_contradictory_copies_are_all_quarantined(change, violation):
    first = candidate("a#1", mag=0.0)
    second = candidate("a#2", **{"mag": 0.0, **change})
    report = dup.find_exact_duplicates([second, first])
    assert {k: (o.status, o.reason) for k, o in report.outcomes.items()} == {
        "a#1": (QUARANTINED, "contradictory_labels"), "a#2": (QUARANTINED, "contradictory_labels")}
    assert report.duplicate_of == {} and report.clusters[0]["kept"] is None
    assert violation in report.outcomes["a#1"].detail


def test_one_bad_copy_quarantines_the_whole_cluster():
    report = dup.find_exact_duplicates([candidate("a#1"), candidate("a#2"), candidate("a#3", energy=-49.0)])
    assert len(report.outcomes) == 3 and {o.reason for o in report.outcomes.values()} == {"contradictory_labels"}


def test_different_magnetic_states_are_never_merged():
    same_order = dup.find_exact_duplicates([candidate("a#1", mag=0.0, state="afm1"),
                                            candidate("a#2", mag=0.0, state="afm2")])
    assert {o.reason for o in same_order.outcomes.values()} == {"contradictory_labels"}
    assert "magnetic_state_id" in same_order.outcomes["a#1"].detail


def test_permuted_copies_compare_aligned_site_sign_patterns():
    # the state id depends on each run's atom order; the aligned site moments decide
    same = dup.find_exact_duplicates([candidate("a#1", mag=0.0, sites=AFM, state="id-in-order-a"),
                                      candidate("b#1", run_id="r:b", perm=PERM, mag=0.0, sites=AFM,
                                                state="id-in-order-b")])
    assert same.duplicate_of == {"b#1": "a#1"}
    flipped = dup.find_exact_duplicates([candidate("a#1", mag=0.0, sites=AFM, state="x"),
                                         candidate("b#1", run_id="r:b", perm=PERM, mag=0.0, sites=-AFM,
                                                   state="y")])
    assert flipped.duplicate_of == {"b#1": "a#1"}  # a global spin flip is the same state
    other = np.array([1.7, 1.7, -1.7, -1.7, 0.0, 0.0, 0.0, 0.0])  # another AFM ordering, same total
    different = dup.find_exact_duplicates([candidate("a#1", mag=0.0, sites=AFM, state="x"),
                                           candidate("b#1", run_id="r:b", perm=PERM, mag=0.0, sites=other,
                                                     state="x")])
    assert {o.reason for o in different.outcomes.values()} == {"contradictory_labels"}
    unaligned = dup.find_exact_duplicates([candidate("a#1", mag=0.0, state="x"),
                                           candidate("b#1", run_id="r:b", perm=PERM, mag=0.0, state="x")])
    assert {o.reason for o in unaligned.outcomes.values()} == {"contradictory_labels"}  # fail closed


def test_missing_magnetization_on_one_copy_is_flagged_unverified():
    report = dup.find_exact_duplicates([candidate("a#1", mag=0.0), candidate("a#2")])
    assert report.duplicate_of == {"a#2": "a#1"}
    assert report.flags == {"a#1": ["duplicate_total_magnetization_unverified"],
                            "a#2": ["duplicate_total_magnetization_unverified"]}


def test_copies_in_different_pools_are_not_duplicates_but_are_linked():
    report = dup.find_exact_duplicates([candidate("a#1", pool="p1"), candidate("b#1", run_id="r:b", pool="p2",
                                                                               energy=-40.0)])
    assert report.outcomes == {} and report.clusters == []
    assert [link["kind"] for link in report.links] == ["same_structure_other_pool"]


def test_duplicate_frame_ids_are_an_error():
    with pytest.raises(DatasetError, match="duplicate frame_id"):
        dup.find_exact_duplicates([candidate("a#1"), candidate("a#1")])


def test_mapping_inputs_and_structure_links():
    k = key()
    frames = [{"frame_id": "a#1", "run_id": "r:a", "pool_id": "p", "permutation_key": k.permutation_key,
               "energy": -1.0},
              {"frame_id": "c#1", "run_id": "r:c", "pool_id": "p", "permutation_key": k.permutation_key,
               "energy": -1.0},
              {"frame_id": "b#1", "run_id": "r:b", "pool_id": "q", "permutation_key": "other"}]
    report = dup.find_exact_duplicates(frames)
    assert report.duplicate_of == {"c#1": "a#1"}
    assert dup.structure_links(frames) == [{"a": "r:a", "b": "r:c", "kind": "shared_structure",
                                            "key": k.permutation_key}]


def test_compare_labels_tolerances_are_configurable():
    a, b = candidate("a#1"), candidate("a#2", energy=-50.002)
    assert dup.compare_labels(a, b)["consistent"] is False
    assert dup.compare_labels(a, b, {"energy": 0.01})["consistent"] is True


# --------------------------------------------------------------------------
# near duplicates
# --------------------------------------------------------------------------

def near(frame_id, group, positions, cell=CELL):
    return dup.NearDuplicateCandidate(frame_id=frame_id, group_id=group, species=SPECIES, cell=cell,
                                      positions=positions, run_id=f"run-{frame_id}")


def test_near_duplicates_only_across_groups():
    pairs = dup.find_near_duplicates([near("a", "g1", POSITIONS), near("b", "g2", POSITIONS + 1e-3),
                                      near("c", "g1", POSITIONS + 2e-3), near("d", "g3", POSITIONS + 0.5)],
                                     tolerance=0.01)
    assert [(p["a"], p["b"]) for p in pairs] == [("a", "b"), ("b", "c")]
    assert pairs[0]["max_displacement"] == pytest.approx(np.sqrt(3) * 1e-3)
    assert pairs[0]["group_a"] == "g1" and pairs[0]["run_b"] == "run-b"


def test_near_duplicates_across_the_periodic_boundary():
    # atom 0 sits at -1e-17 A (fractional mod 1 rounds to exactly 1.0) in "b" and at +1e-4 A in "a"
    tiny = POSITIONS.copy()
    tiny[0] = [-4.17e-17, 0.0, 0.0]
    shifted = POSITIONS.copy()
    shifted[0] = [1e-4, 0.0, 0.0]
    pairs = dup.find_near_duplicates([near("b", "g1", tiny), near("a", "g2", shifted)], tolerance=0.01)
    assert [(p["a"], p["b"]) for p in pairs] == [("a", "b")]
    wrapped = POSITIONS.copy()
    wrapped[0] = [4.17 - 1e-4, 0.0, 0.0]  # the other side of the cell: 2e-4 A away under minimum image
    pairs = dup.find_near_duplicates([near("a", "g1", shifted), near("b", "g2", wrapped)], tolerance=0.01)
    assert len(pairs) == 1 and pairs[0]["max_displacement"] == pytest.approx(2e-4)


def test_near_duplicate_arguments_are_checked():
    with pytest.raises(DatasetError, match="positive"):
        dup.find_near_duplicates([], tolerance=0)
    frames = [near(f"f{i}", f"g{i}", POSITIONS) for i in range(4)]
    with pytest.raises(DatasetError, match="exceeded"):
        dup.find_near_duplicates(frames, tolerance=0.01, max_pairs=2)
    assert dup.find_near_duplicates([near("a", "g1", POSITIONS), near("b", "g2", POSITIONS, cell=CELL * 2)],
                                    tolerance=0.01) == []
