"""Structure keys, exact duplicates, label consistency and near-duplicates.

Keys (spec "Duplicates"): species sequence, pbc, cell rounded to 1e-6 A, and
positions wrapped along periodic axes (fractional into [0, 1), values within
1e-9 of 1 mapped to 0) then Cartesian rounded to 1e-6 A. ``order_key`` keeps
the atom order; ``permutation_key`` sorts atoms by (species, rounded
coordinates), so a relabelled copy of the same structure has the same
permutation key. Rounding is done on integers (``rint(x * 1e6)``), so -0.0 and
+0.0 cannot give different keys.

Exact duplicates are equal permutation keys within one reference pool. Their
labels are compared with numerical-noise tolerances (defaults 1e-3 eV,
1e-2 eV/A, 0.1 muB -- bounds for EDIFF 1e-5..1e-6 runs, configurable):
consistent copies keep the lowest ``frame_id`` and the others become
``excluded/exact_duplicate`` (``duplicate_of``); inconsistent copies (energy,
forces, total magnetization or site sign pattern) are all
``quarantined/contradictory_labels``. Every duplicate cluster, and every
structure shared across pools, is reported as a lineage link so the splitter
keeps the copies together.

Near-duplicates (opt-in) compare frames with the same species order, cell and
pbc across different groups only: max minimum-image displacement < tolerance.
They only union groups and are reported; they never remove a frame and never
merge labels of different magnetic states.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .errors import DatasetError
from .model import Outcome, outcome_for

KEY_SCHEMA = "nio-md-prep.structure-key/v1"
ROUND_A = 1e-6  # Angstrom
WRAP_SNAP = 1e-9  # fractional values within this of 1 are mapped to 0

#: energy eV, force eV/A (per component), magnetization muB (total moment); site_sign_threshold muB
#: (|m| below it counts as "no moment" when permuted copies' site sign patterns are compared).
DEFAULT_TOLERANCES = {"energy": 1e-3, "force": 1e-2, "magnetization": 0.1, "site_sign_threshold": 0.5}


# --------------------------------------------------------------------------
# Keys
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class StructureKey:
    order_key: str
    permutation_key: str
    #: ``permutation[k]`` = original atom index of the k-th atom in canonical (sorted) order;
    #: aligns per-atom arrays (forces, moments) of permuted copies.
    permutation: tuple[int, ...]

    def as_dict(self) -> dict[str, Any]:
        return {"order_key": self.order_key, "permutation_key": self.permutation_key}


def _wrapped_fractional(cell: np.ndarray, pbc: Sequence[bool], positions=None, fractional=None) -> np.ndarray:
    if fractional is None:
        if positions is None:
            raise DatasetError("structure key needs positions or fractional coordinates")
        fractional = np.linalg.solve(cell.T, np.asarray(positions, dtype=np.float64).T).T
    frac = np.array(fractional, dtype=np.float64, copy=True)
    if frac.ndim != 2 or frac.shape[1] != 3:
        raise DatasetError(f"fractional coordinates must be (N, 3), got {frac.shape}")
    for axis in range(3):
        if pbc[axis]:
            column = np.mod(frac[:, axis], 1.0)
            column[np.abs(column - 1.0) <= WRAP_SNAP] = 0.0
            frac[:, axis] = column
    return frac


def _digest(kind: str, species: Sequence[str], pbc: Sequence[bool], cell_int: np.ndarray, pos_int: np.ndarray) -> str:
    header = json.dumps({"schema": KEY_SCHEMA, "kind": kind, "species": list(species),
                         "pbc": [bool(p) for p in pbc]}, sort_keys=True).encode()
    hasher = hashlib.sha256(header)
    hasher.update(np.ascontiguousarray(cell_int, dtype="<i8").tobytes())
    hasher.update(np.ascontiguousarray(pos_int, dtype="<i8").tobytes())
    return hasher.hexdigest()


def structure_fingerprint(
    species: Sequence[str],
    cell: Any,
    pbc: Sequence[bool] = (True, True, True),
    positions: Any = None,
    *,
    fractional: Any = None,
) -> StructureKey:
    """Order/permutation keys (and the canonical atom permutation) of one structure.

    Pass ``fractional`` when available (vasprun stores fractional coordinates)
    so no matrix inversion enters the key.
    """
    species = [str(s).strip() for s in species]
    cell = np.asarray(cell, dtype=np.float64)
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
        raise DatasetError("structure key needs a finite 3x3 cell")
    pbc = [bool(p) for p in pbc]
    frac = _wrapped_fractional(cell, pbc, positions, fractional)
    if frac.shape[0] != len(species):
        raise DatasetError(f"{frac.shape[0]} positions for {len(species)} species")
    if not np.all(np.isfinite(frac)):
        raise DatasetError("structure key needs finite coordinates")
    cartesian = frac @ cell
    pos_int = np.rint(cartesian / ROUND_A).astype(np.int64)
    cell_int = np.rint(cell / ROUND_A).astype(np.int64)
    order_key = _digest("order", species, pbc, cell_int, pos_int)
    rank = {element: i for i, element in enumerate(sorted(set(species)))}
    codes = np.array([rank[s] for s in species], dtype=np.int64)
    # lexsort: last key is primary -> (species, x, y, z)
    permutation = np.lexsort((pos_int[:, 2], pos_int[:, 1], pos_int[:, 0], codes)) if len(species) else \
        np.zeros(0, dtype=np.int64)
    sorted_species = [species[i] for i in permutation]
    permutation_key = _digest("permutation", sorted_species, pbc, cell_int, pos_int[permutation])
    return StructureKey(order_key, permutation_key, tuple(int(i) for i in permutation))


def structure_keys(species: Sequence[str], cell: Any, pbc: Sequence[bool] = (True, True, True),
                   positions: Any = None, *, fractional: Any = None) -> tuple[str, str]:
    """``(order_key, permutation_key)`` (see :func:`structure_fingerprint`)."""
    key = structure_fingerprint(species, cell, pbc, positions, fractional=fractional)
    return key.order_key, key.permutation_key


def step_structure_key(species: Sequence[str], structure: Any, pbc: Sequence[bool] = (True, True, True)) -> StructureKey | None:
    """Key of a :class:`~nio_md_prep.dataset.model.StepStructure` (None if it cannot be keyed)."""
    if structure is None:
        return None
    try:
        return structure_fingerprint(species, structure.cell, pbc, getattr(structure, "positions", None),
                                     fractional=getattr(structure, "fractional", None))
    except (DatasetError, np.linalg.LinAlgError, TypeError, ValueError):
        return None


# --------------------------------------------------------------------------
# Label comparison
# --------------------------------------------------------------------------

def _get(record: Any, key: str, default: Any = None) -> Any:
    if isinstance(record, Mapping):
        return record.get(key, default)
    return getattr(record, key, default)


@dataclass
class DuplicateCandidate:
    """One frame as seen by the duplicate finder (build from FrameRecord + arrays)."""

    frame_id: str
    run_id: str
    pool_id: str | None
    permutation_key: str
    order_key: str | None = None
    permutation: Sequence[int] | None = None  # StructureKey.permutation
    energy: float | None = None  # the label energy (same quantity for the whole pool)
    forces: Any = None  # (N, 3) raw forces
    total_magnetization: float | None = None
    site_magmoms: Sequence[float] | None = None
    magnetic_state_id: str | None = None
    group_id: str | None = None


def _candidate(item: Any) -> DuplicateCandidate:
    if isinstance(item, DuplicateCandidate):
        return item
    return DuplicateCandidate(**{name: _get(item, name) for name in DuplicateCandidate.__dataclass_fields__
                                 if _get(item, name) is not None or name in {"frame_id", "run_id", "pool_id",
                                                                             "permutation_key"}})


def _aligned(values: Any, permutation: Sequence[int] | None) -> np.ndarray | None:
    if values is None:
        return None
    array = np.asarray(values, dtype=np.float64)
    if permutation is None:
        return array
    return array[np.asarray(permutation, dtype=np.int64)]


def _aligned_sign_pattern(moments: Any, permutation: Sequence[int] | None, threshold: float) -> tuple[int, ...] | None:
    """Site sign pattern in canonical atom order (+1/-1, 0 below ``threshold``), invariant under a global flip."""
    if moments is None or permutation is None:
        return None
    values = _aligned(moments, permutation)
    if values is None or values.ndim != 1:
        return None
    signs = [0 if abs(float(v)) < threshold else (1 if v > 0 else -1) for v in values]
    first = next((s for s in signs if s != 0), 1)
    return tuple(s * first for s in signs)


def compare_labels(a: Any, b: Any, tolerances: Mapping[str, float] | None = None) -> dict[str, Any]:
    """Compare two structural duplicates' labels; ``consistent`` is False on any violated bound.

    Per-atom arrays are aligned through the canonical permutation when both
    candidates carry one (or when the atom orders are identical). Forces that
    cannot be aligned are "unverifiable", which counts as inconsistent (fail
    closed). Magnetization is compared when both are known: total moment
    within the tolerance and identical site sign patterns (magnetic_state_id).
    """
    tol = dict(DEFAULT_TOLERANCES)
    tol.update(tolerances or {})
    a, b = _candidate(a), _candidate(b)
    result: dict[str, Any] = {"dE": None, "max_dF": None, "d_mag": None, "violations": [], "unverified": []}
    if a.energy is not None and b.energy is not None:
        result["dE"] = abs(float(a.energy) - float(b.energy))
        if not result["dE"] <= tol["energy"]:
            result["violations"].append(f"|dE| {result['dE']:.3g} > {tol['energy']:g} eV")
    else:
        result["unverified"].append("energy")
    same_order = a.order_key is not None and a.order_key == b.order_key
    if a.forces is not None and b.forces is not None:
        if same_order:
            fa, fb = _aligned(a.forces, None), _aligned(b.forces, None)
        elif a.permutation is not None and b.permutation is not None:
            fa, fb = _aligned(a.forces, a.permutation), _aligned(b.forces, b.permutation)
        else:
            fa = fb = None
        if fa is None or fb is None or fa.shape != fb.shape:
            result["violations"].append("forces cannot be aligned (atom permutation unknown)")
        else:
            result["max_dF"] = float(np.max(np.abs(fa - fb))) if fa.size else 0.0
            if not result["max_dF"] <= tol["force"]:
                result["violations"].append(f"max |dF| {result['max_dF']:.3g} > {tol['force']:g} eV/A")
    elif (a.forces is None) != (b.forces is None):
        result["violations"].append("one copy has forces, the other has none")
    if a.total_magnetization is not None and b.total_magnetization is not None:
        result["d_mag"] = abs(float(a.total_magnetization) - float(b.total_magnetization))
        if not result["d_mag"] <= tol["magnetization"]:
            result["violations"].append(f"|d mag_total| {result['d_mag']:.3g} > {tol['magnetization']:g} muB")
    elif (a.total_magnetization is None) != (b.total_magnetization is None):
        result["unverified"].append("total_magnetization")
    if a.magnetic_state_id is not None and b.magnetic_state_id is not None:
        if same_order:
            if a.magnetic_state_id != b.magnetic_state_id:
                result["violations"].append("different site sign patterns (magnetic_state_id)")
        else:
            # magnetic_state_id is built from the sign pattern in each run's own atom order, so permuted
            # copies are compared on their site moments aligned through the canonical permutation.
            pa = _aligned_sign_pattern(a.site_magmoms, a.permutation, tol["site_sign_threshold"])
            pb = _aligned_sign_pattern(b.site_magmoms, b.permutation, tol["site_sign_threshold"])
            if pa is None or pb is None or len(pa) != len(pb):
                result["violations"].append("site sign patterns cannot be aligned (atom permutation or moments unknown)")
            elif pa != pb:
                result["violations"].append("different site sign patterns (aligned site moments)")
    elif (a.magnetic_state_id is None) != (b.magnetic_state_id is None):
        result["unverified"].append("magnetic_state_id")
    result["consistent"] = not result["violations"]
    return result


# --------------------------------------------------------------------------
# Exact duplicates
# --------------------------------------------------------------------------

@dataclass
class DuplicateReport:
    outcomes: dict[str, Outcome] = field(default_factory=dict)  # frame_id -> exact_duplicate | contradictory_labels
    duplicate_of: dict[str, str] = field(default_factory=dict)  # removed frame_id -> kept frame_id
    clusters: list[dict[str, Any]] = field(default_factory=list)
    links: list[dict[str, Any]] = field(default_factory=list)  # {a, b, kind, key} between run_ids (lineage union)
    flags: dict[str, list[str]] = field(default_factory=dict)  # frame_id -> flags (e.g. magnetization_unverified)

    def as_dict(self) -> dict[str, Any]:
        return {
            "outcomes": {k: v.as_dict() for k, v in sorted(self.outcomes.items())},
            "duplicate_of": dict(sorted(self.duplicate_of.items())),
            "clusters": self.clusters,
            "links": self.links,
            "flags": {k: sorted(v) for k, v in sorted(self.flags.items())},
            "counts": {
                "clusters": len(self.clusters),
                "exact_duplicate": sum(1 for o in self.outcomes.values() if o.reason == "exact_duplicate"),
                "contradictory_labels": sum(1 for o in self.outcomes.values() if o.reason == "contradictory_labels"),
                "links": len(self.links),
            },
        }


def _link(a: str, b: str, kind: str, key: str) -> dict[str, Any]:
    x, y = sorted((a, b))
    return {"a": x, "b": y, "kind": kind, "key": key}


def find_exact_duplicates(frames: Iterable[Any], tolerances: Mapping[str, float] | None = None) -> DuplicateReport:
    """Exact duplicates (equal permutation key within one pool) and their label consistency.

    ``frames``: :class:`DuplicateCandidate` objects or mappings/records with
    the same fields -- normally the frames still accepted after the frame-level
    rules. Output is independent of input order.
    """
    candidates = sorted((_candidate(f) for f in frames), key=lambda c: c.frame_id)
    seen: set[str] = set()
    for candidate in candidates:
        if candidate.frame_id in seen:
            raise DatasetError(f"duplicate frame_id {candidate.frame_id!r} in duplicate search")
        seen.add(candidate.frame_id)
    report = DuplicateReport()
    by_key: dict[str, list[DuplicateCandidate]] = {}
    for candidate in candidates:
        if candidate.permutation_key:
            by_key.setdefault(candidate.permutation_key, []).append(candidate)
    link_set: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for key in sorted(by_key):
        members = by_key[key]
        if len(members) < 2:
            continue
        pools: dict[Any, list[DuplicateCandidate]] = {}
        for member in members:
            pools.setdefault(member.pool_id, []).append(member)
        runs = sorted({m.run_id for m in members})
        for i, run_a in enumerate(runs):  # a structure shared across runs links them, whatever the pool
            for run_b in runs[i + 1:]:
                kind = "exact_duplicate" if any(
                    len({m.run_id for m in group} & {run_a, run_b}) == 2 for group in pools.values()
                ) else "same_structure_other_pool"
                entry = _link(run_a, run_b, kind, key)
                link_set[(entry["a"], entry["b"], kind, key)] = entry
        for pool_id in sorted(pools, key=lambda p: (p is None, str(p))):
            group = pools[pool_id]
            if len(group) < 2:
                continue
            kept = group[0]  # lowest frame_id (candidates are sorted)
            comparisons = []
            consistent = True
            for i, first in enumerate(group):
                for second in group[i + 1:]:
                    comparison = compare_labels(first, second, tolerances)
                    comparisons.append({"a": first.frame_id, "b": second.frame_id,
                                        **{k: comparison[k] for k in ("dE", "max_dF", "d_mag", "violations",
                                                                      "unverified")}})
                    consistent = consistent and comparison["consistent"]
                    for name in comparison["unverified"]:
                        for frame_id in (first.frame_id, second.frame_id):
                            report.flags.setdefault(frame_id, [])
                            flag = f"duplicate_{name}_unverified"
                            if flag not in report.flags[frame_id]:
                                report.flags[frame_id].append(flag)
            member_ids = [m.frame_id for m in group]
            if consistent:
                for member in group[1:]:
                    report.outcomes[member.frame_id] = outcome_for("exact_duplicate", f"duplicate of {kept.frame_id}")
                    report.duplicate_of[member.frame_id] = kept.frame_id
            else:
                for member in group:
                    report.outcomes[member.frame_id] = outcome_for(
                        "contradictory_labels", f"{len(group)} copies of one structure disagree: "
                        + "; ".join(sorted({v for c in comparisons for v in c["violations"]}))[:300])
            report.clusters.append({
                "permutation_key": key, "pool_id": pool_id, "frames": member_ids,
                "kept": kept.frame_id if consistent else None, "consistent": consistent,
                "comparisons": comparisons,
            })
    report.links = [link_set[k] for k in sorted(link_set)]
    return report


def structure_links(frames: Iterable[Any]) -> list[dict[str, Any]]:
    """Run pairs that share any structure (permutation key), regardless of pool or outcome.

    Used for lineage union: a structure that appears in two runs correlates them
    even when one copy is rejected or belongs to another pool.
    """
    runs_by_key: dict[str, set[str]] = {}
    for frame in frames:
        key = _get(frame, "permutation_key")
        run_id = _get(frame, "run_id")
        if key and run_id:
            runs_by_key.setdefault(key, set()).add(run_id)
    links: dict[tuple[str, str], dict[str, Any]] = {}
    for key in sorted(runs_by_key):
        runs = sorted(runs_by_key[key])
        for i, run_a in enumerate(runs):
            for run_b in runs[i + 1:]:
                links.setdefault((run_a, run_b), _link(run_a, run_b, "shared_structure", key))
    return [links[k] for k in sorted(links)]


# --------------------------------------------------------------------------
# Near duplicates
# --------------------------------------------------------------------------

@dataclass
class NearDuplicateCandidate:
    frame_id: str
    group_id: str
    species: Sequence[str]
    cell: Any
    positions: Any
    pbc: Sequence[bool] = (True, True, True)
    run_id: str | None = None


def _near(item: Any) -> NearDuplicateCandidate:
    if isinstance(item, NearDuplicateCandidate):
        return item
    return NearDuplicateCandidate(
        frame_id=_get(item, "frame_id"), group_id=_get(item, "group_id"), species=_get(item, "species"),
        cell=_get(item, "cell"), positions=_get(item, "positions"),
        pbc=_get(item, "pbc", (True, True, True)) or (True, True, True), run_id=_get(item, "run_id"),
    )


def find_near_duplicates(frames: Iterable[Any], tolerance: float, *, max_pairs: int = 5_000_000) -> list[dict[str, Any]]:
    """Cross-group pairs with the same species order/cell/pbc and max min-image displacement < ``tolerance``.

    A spatial hash on atom 0's fractional position (bin width >= tolerance
    along each plane normal, periodic neighbours) limits the comparisons to
    plausible pairs; ``max_pairs`` bounds the work (DatasetError beyond it).
    Returns ``[{a, b, group_a, group_b, max_displacement}]`` sorted by frame id.
    """
    if not (isinstance(tolerance, (int, float)) and math.isfinite(tolerance) and tolerance > 0):
        raise DatasetError(f"near-duplicate tolerance must be a positive number, got {tolerance!r}")
    candidates = sorted((_near(f) for f in frames), key=lambda c: c.frame_id)
    buckets: dict[tuple, list[NearDuplicateCandidate]] = {}
    for candidate in candidates:
        cell = np.asarray(candidate.cell, dtype=np.float64)
        key = (tuple(str(s) for s in candidate.species),
               tuple(int(v) for v in np.rint(cell / ROUND_A).astype(np.int64).ravel()),
               tuple(bool(p) for p in candidate.pbc))
        buckets.setdefault(key, []).append(candidate)
    pairs: list[dict[str, Any]] = []
    compared = 0
    for key in sorted(buckets, key=lambda k: json.dumps([list(k[0]), list(k[1]), list(k[2])])):
        members = buckets[key]
        if len(members) < 2 or len({m.group_id for m in members}) < 2 or not key[0]:
            continue
        cell = np.asarray(members[0].cell, dtype=np.float64)
        inverse = np.linalg.inv(cell)
        pbc = np.array(key[2], dtype=bool)
        spacing = 1.0 / np.linalg.norm(inverse.T, axis=1)
        nbins = np.maximum(1, np.floor(spacing / tolerance).astype(np.int64))
        positions = [np.asarray(m.positions, dtype=np.float64) for m in members]
        grid: dict[tuple[int, int, int], list[int]] = {}
        cells_of = []
        for i, pos in enumerate(positions):
            frac0 = pos[0] @ inverse
            frac0 = np.where(pbc, np.mod(frac0, 1.0), frac0)
            raw = np.floor(frac0 * nbins).astype(np.int64)
            raw = np.where(pbc, np.mod(raw, nbins), raw)  # mod(-tiny, 1.0) can round to exactly 1.0
            b = tuple(int(v) for v in raw)
            cells_of.append(b)
            grid.setdefault(b, []).append(i)
        for i, pos_i in enumerate(positions):
            neighbours: set[int] = set()
            base = cells_of[i]
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        cell_index = []
                        for axis, delta in enumerate((dx, dy, dz)):
                            value = base[axis] + delta
                            if pbc[axis]:
                                value %= int(nbins[axis])
                            cell_index.append(value)
                        neighbours.update(grid.get(tuple(cell_index), ()))
            for j in sorted(neighbours):
                if j <= i or members[i].group_id == members[j].group_id:
                    continue
                compared += 1
                if compared > max_pairs:
                    raise DatasetError(f"near-duplicate search exceeded {max_pairs} comparisons; "
                                       "use a smaller tolerance or restrict the frames")
                delta = (positions[j] - pos_i) @ inverse
                delta[:, pbc] -= np.rint(delta[:, pbc])
                displacement = float(np.max(np.linalg.norm(delta @ cell, axis=1)))
                if displacement < tolerance:
                    a, b = members[i], members[j]
                    pairs.append({"a": a.frame_id, "b": b.frame_id, "group_a": a.group_id, "group_b": b.group_id,
                                  "run_a": a.run_id, "run_b": b.run_id, "max_displacement": displacement})
    return sorted(pairs, key=lambda p: (p["a"], p["b"]))
