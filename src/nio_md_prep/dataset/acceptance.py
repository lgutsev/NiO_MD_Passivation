"""Run- and frame-level acceptance policy for VASP training labels.

Every discovered run and every ionic frame gets exactly one
:class:`~nio_md_prep.dataset.model.Outcome`; non-accepted outcomes are always
built with :func:`~nio_md_prep.dataset.model.outcome_for`, so a reason code is
reported under exactly one status. The policy is conservative ("strict"):
anything that cannot be shown to be a converged, force-consistent DFT label is
not accepted.

Run level (:func:`classify_run`), in precedence order:

1. header unreadable -> ``parse_error/vasprun_unreadable``; program not VASP ->
   ``rejected/not_vasp``;
2. calculation type from IBRION/NSW; DFPT/linear response (IBRION 5-8,
   LEPSILON, LCALCEPS, LCHIMAG) and non-allowed ``other`` types ->
   ``rejected/unsupported_calculation``;
3. label file byte-identical to a reference/template file ->
   ``excluded/copied_reference_output`` (copied *evidence* files are ignored);
4. VASP-MLFF: ``ML_MODE`` run/select/refit*/delta or an undeterminable mode ->
   ``non_dft_mlff_step/mlff_steps``; ``train`` -> per step (flat MLFF steps are
   ``non_dft_mlff_step/mlff_step``, never DFT labels, but still counted);
5. species/order disagreement (vasprun, POTCAR, POSCAR, OUTCAR) ->
   ``quarantined/species_mismatch``;
6. OUTCAR/OSZICAR/POSCAR evidence contradicting vasprun (NIONS, executed
   settings, per-step free energies and forces, selective-dynamics flags) ->
   ``quarantined/evidence_mismatch``;
7. incomplete run (vasprun not closed or OUTCAR without the timing block) with
   policy ``exclude`` -> ``excluded/run_incomplete``; ``recover`` judges frames
   individually (the truncated tail is always ``rejected/truncated_frame``);
8. magnetic policy: a run with magnetic species whose class is ``uncontrolled``
   or ``unknown`` -> ``quarantined/magnetic_uncontrolled`` /
   ``quarantined/magnetic_unknown`` unless the matching override is set.

Frame level (:func:`classify_step`), first failing rule wins:
``mlff_step`` -> ``truncated_frame`` -> ``missing_forces``/``missing_energy`` ->
``bad_shape`` -> ``non_finite`` -> ``degenerate_cell`` -> SCF gate -> magnetic
policy -> ``force_outlier`` (opt-in) -> the run outcome (a non-accepted run
passes its outcome on to otherwise acceptable frames). Duplicates, reference
pools and subsampling are decided later by the exporter.

SCF gate (research S2.1, amendments item 7) per DFT step, from vasprun
(``<scstep>`` count n vs NELM, |dE| of the last two scsteps vs EDIFF) combined
with OUTCAR ``aborting loop`` markers, which are mapped to steps only when the
OUTCAR DFT-step count matches the vasprun DFT-step count. ``EDIFF is reached``
is proof only for VASP >= 6 (VASP 5 prints it even after hitting NELM).
Explicit failure -> ``rejected/scf_not_converged``; missing or contradictory
evidence -> ``quarantined/scf_convergence_unknown``.

Magnetism is an explicit audit dimension (not "jump = bad"): per run a class
(:data:`~nio_md_prep.dataset.model.MAGNETIC_CLASSES`) and per frame the total
moment, site moments, changes vs the first and the previous frame, the sign
pattern and a policy result. Default NiO policy: runs containing magnetic
species that are ``uncontrolled`` (ISPIN=2 without MAGMOM/NUPDOWN) or
``unknown`` are quarantined; ``controlled_transition`` runs have the frames from
the first transition on quarantined (a different branch). Overrides are
explicit :class:`Policy` switches recorded in the manifest.

Energies come from :func:`nio_md_prep.dataset.vasprun.derived_energies`, which
owns the VASP <= 6.0.8 calc-level mislabel rule and records ``energy_source``,
``energy_rule`` and the parser name/version for every frame.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
import hashlib
import json
import math
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .errors import DatasetError
from .model import (
    ACCEPTED,
    MAG_CONTROLLED_CONSISTENT,
    MAG_CONTROLLED_TRANSITION,
    MAG_NON_SPIN_POLARIZED,
    MAG_UNCONTROLLED,
    MAG_UNKNOWN,
    Outcome,
    outcome_for,
    vasp_stress_to_ase,
)
from .settings import UNKNOWN, as_bool, as_float, as_float_list, as_int, as_str, extract_settings, json_safe
from .settings import outcar_conflicts, titel_element
from . import vasprun as _vasprun

ACCEPTANCE_SCHEMA = 1

CALC_TYPES = ("static", "relaxation", "md", "other")
ENERGY_QUANTITIES = ("free_energy", "energy_sigma0", "energy_no_entropy")
INTERRUPTED_POLICIES = ("exclude", "recover")
#: Transition-metal default for NiO work (user requirement); elements whose explicit
#: |MAGMOM| >= Policy.magmom_species_threshold are added per run.
DEFAULT_MAGNETIC_SPECIES = ("Co", "Cr", "Cu", "Fe", "Mn", "Ni", "V")

LABEL_SET_ENERGY_FORCES = "energy_forces"
LABEL_SET_ENERGY_ONLY = "energy_only"
#: forces_reason of an accepted energy-only frame (allow_energy_only): VASP wrote no forces for the step.
FORCES_REASON_ENERGY_ONLY = "no <varray name=\"forces\"> in this <calculation> (accepted with allow_energy_only)"

#: stress_reason values (stress_available=False); None when the stress is exported.
STRESS_REASONS = {
    "not_computed": "no <varray name=\"stress\"> for this step (ISIF=0, the default for IBRION=0/LHFCALC)",
    "isif_unknown": "ISIF could not be determined, so the tensor's validity is unknown",
    "isif1_trace_only": "ISIF=1: only the trace (pressure) is valid",
    "pstress_nonzero": "PSTRESS != 0: whether the printed tensor includes it is undocumented (--allow-pstress)",
    "vacuum_slab": "vacuum gap detected; stress of a slab/cluster is diluted by the vacuum (--allow-vacuum-stress)",
    "not_requested": "stress labels were not requested (--include-stress)",
    "not_a_dft_label": "VASP-MLFF prediction",
}
STRESS_SOURCE = "vasprun:calculation.varray[stress] (kB, VASP sign) -> ASE sign eV/A^3, symmetrized"

#: SCF evidence summary values
SCF_EVIDENCE = ("scstep_count+outcar_marker", "scstep_count", "outcar_marker", "none")

_ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$")
_ML_ISTART_MODE = {0: "train", 1: "train", 2: "run", 3: "select", 4: "refit"}
_ML_NO_DFT = frozenset({"run", "select", "refit", "refitbayesian"})

#: Run-level reasons in precedence order (the first finding decides the run outcome).
RUN_REASON_PRECEDENCE = (
    "vasprun_unreadable", "not_vasp", "unsupported_calculation", "copied_reference_output", "mlff_steps",
    "species_mismatch", "evidence_mismatch", "run_incomplete", "magnetic_uncontrolled", "magnetic_unknown",
)
#: Run outcomes whose frames are not DFT training labels at all: frames inherit them directly
#: (before the frame-level label checks, which would be meaningless there).
NON_LABEL_RUN_REASONS = frozenset({
    "vasprun_unreadable", "not_vasp", "unsupported_calculation", "copied_reference_output", "mlff_steps",
})


# --------------------------------------------------------------------------
# Policy
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Policy:
    """Every switch and threshold of the acceptance policy (recorded verbatim in the manifest).

    Defaults are the strict NiO policy. Magnetic overrides
    (``accept_uncontrolled``, ``accept_unknown``, ``accept_transitions``) are
    deliberate, recorded decisions; ``magnetic_override_reason`` documents why.
    """

    strict: bool = True
    interrupted_run_policy: str = "exclude"  # exclude | recover
    allowed_calc_types: tuple[str, ...] = ("md", "relaxation", "static")
    energy_quantity: str = "free_energy"
    allow_energy_only: bool = False
    # SCF gate
    allow_ediff_zero: bool = False
    require_outcar_marker_vasp6: bool = True
    #: VASP with NWRITE <= 1 may print the per-step 'aborting loop' line only for the first ionic step;
    #: then (and only then) a missing marker falls back to the vasprun evidence (flag scf_marker_not_printed).
    nwrite_marker_exemption: bool = True
    # evidence cross-checks
    outcar_energy_tolerance: float = 5e-6  # eV, OUTCAR TOTEN vs F (covers 6-decimal VASP 5.2 printing)
    outcar_force_tolerance: float = 1e-6  # eV/A per component (compared as sum |F| over 3N components)
    oszicar_relative_tolerance: float = 1e-7  # OSZICAR prints 8 significant digits
    additive_correction_tolerance: float = 1e-6  # eV, unexplained calc e_fr - PV - scstep e_fr
    # stress
    include_stress: bool = False
    allow_vacuum_stress: bool = False
    allow_pstress: bool = False
    vacuum_gap_threshold: float = 5.0  # A, empty slab thickness along a lattice-plane normal
    # magnetism
    magnetic_species: tuple[str, ...] = DEFAULT_MAGNETIC_SPECIES
    magmom_species_threshold: float = 0.5  # muB; explicit |MAGMOM| >= this makes an element magnetic
    magmom_sign_threshold: float = 0.5  # muB; |m| below this is "no moment" in sign patterns
    mag_jump_threshold: float = 0.5  # muB, |total_k - total_(k-1)|
    mag_drift_threshold: float = 1.0  # muB, |total_k - total_first|
    mag_initial_total_tolerance: float = 1.0  # muB, |total_first - sum(MAGMOM)| for balanced (AFM) MAGMOM
    accept_uncontrolled: bool = False
    accept_unknown: bool = False
    accept_transitions: bool = False
    magnetic_override_reason: str | None = None
    # outliers
    force_outlier_threshold: float | None = None  # eV/A; None = off
    # duplicates (used by duplicates.py via the exporter)
    duplicate_energy_tolerance: float = 1e-3  # eV
    duplicate_force_tolerance: float = 1e-2  # eV/A
    duplicate_magnetization_tolerance: float = 0.1  # muB
    near_duplicate_tolerance: float | None = None  # A; None = off

    def __post_init__(self):
        def put(name, value):
            object.__setattr__(self, name, value)

        if self.interrupted_run_policy not in INTERRUPTED_POLICIES:
            raise DatasetError(f"interrupted_run_policy must be one of {INTERRUPTED_POLICIES}, "
                               f"got {self.interrupted_run_policy!r}")
        if self.energy_quantity not in ENERGY_QUANTITIES:
            raise DatasetError(f"energy_quantity must be one of {ENERGY_QUANTITIES}, got {self.energy_quantity!r}")
        calc_types = tuple(sorted(set(self.allowed_calc_types)))
        unknown = sorted(set(calc_types) - set(CALC_TYPES))
        if unknown or not calc_types:
            raise DatasetError(f"allowed_calc_types must be a non-empty subset of {CALC_TYPES}, got {calc_types}")
        put("allowed_calc_types", calc_types)
        species = tuple(sorted(set(str(s).strip() for s in self.magnetic_species)))
        bad = [s for s in species if not _ELEMENT_RE.match(s)]
        if bad:
            raise DatasetError(f"magnetic_species must be element symbols, got {bad}")
        put("magnetic_species", species)
        for spec in fields(self):
            value = getattr(self, spec.name)
            if spec.name in {"strict", "allow_energy_only", "allow_ediff_zero", "require_outcar_marker_vasp6",
                             "nwrite_marker_exemption", "include_stress", "allow_vacuum_stress", "allow_pstress",
                             "accept_uncontrolled", "accept_unknown", "accept_transitions"}:
                if not isinstance(value, bool):
                    raise DatasetError(f"policy {spec.name} must be true or false, got {value!r}")
            elif spec.name.endswith(("_tolerance", "_threshold")):
                optional = spec.name in {"force_outlier_threshold", "near_duplicate_tolerance"}
                if value is None and optional:
                    continue
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) \
                        or value <= 0:
                    raise DatasetError(f"policy {spec.name} must be a positive number, got {value!r}")
                put(spec.name, float(value))
        if self.magnetic_override_reason is not None and not str(self.magnetic_override_reason).strip():
            put("magnetic_override_reason", None)

    @property
    def magnetic_overrides(self) -> list[str]:
        return [name for name in ("accept_uncontrolled", "accept_unknown", "accept_transitions") if getattr(self, name)]

    def as_dict(self) -> dict[str, Any]:
        data = {spec.name: json_safe(getattr(self, spec.name)) for spec in fields(self)}
        data["schema"] = ACCEPTANCE_SCHEMA
        data["non_default"] = self.non_default()
        return data

    def non_default(self) -> dict[str, Any]:
        """Switches that differ from the strict defaults (what a reviewer must look at)."""
        default = Policy()
        return {spec.name: json_safe(getattr(self, spec.name)) for spec in fields(self)
                if getattr(self, spec.name) != getattr(default, spec.name)}

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "Policy":
        """Build from a TOML/JSON table; unknown keys are errors (typo protection)."""
        names = {spec.name for spec in fields(cls)}
        extra = sorted(set(data) - names - {"schema", "non_default"})
        if extra:
            raise DatasetError(f"unknown policy keys {extra}; allowed: {sorted(names)}")
        values = {key: (tuple(value) if isinstance(value, list) else value)
                  for key, value in data.items() if key in names}
        return cls(**values)


DEFAULT_POLICY = Policy()


# --------------------------------------------------------------------------
# Facts of a run (typed values from the settings record)
# --------------------------------------------------------------------------

def _get(record: Any, key: str, default: Any = None) -> Any:
    if isinstance(record, Mapping):
        return record.get(key, default)
    return getattr(record, key, default)


def _known(value: Any) -> Any:
    return None if value == UNKNOWN else value


def run_facts(header: Any, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The typed INCAR/parameter facts the policy uses (None = unknown; no VASP defaults are guessed)."""
    sampling = settings.get("sampling", {}) if settings else {}
    method = settings.get("method", {}) if settings else {}
    parameters = {str(k).upper(): v for k, v in (_get(header, "parameters") or {}).items()}

    def samp(name, coerce):
        value = sampling.get(name)
        return coerce(value) if value is not None else None

    ml_mode = samp("ML_MODE", as_str)
    facts = {
        "program": (as_str((_get(header, "generator") or {}).get("program")) or "").strip().lower() or None,
        "version": ".".join(map(str, header.version_tuple)) if _get(header, "version_tuple") else None,
        "version_tuple": tuple(_get(header, "version_tuple")) if _get(header, "version_tuple") else None,
        "ibrion": samp("IBRION", as_int),
        "nsw": samp("NSW", as_int),
        "isif": samp("ISIF", as_int),
        "potim": samp("POTIM", as_float),
        "nelm": samp("NELM", as_int),
        "nelmin": samp("NELMIN", as_int),
        "ediff": samp("EDIFF", as_float),
        "pstress": samp("PSTRESS", as_float),
        "lorbit": samp("LORBIT", as_int),
        "nwrite": samp("NWRITE", as_int),
        "istart": samp("ISTART", as_int),
        "icharg": samp("ICHARG", as_int),
        "tebeg": samp("TEBEG", as_float),
        "teend": samp("TEEND", as_float),
        "ml_lmlff": samp("ML_LMLFF", as_bool),
        "ml_mode": ml_mode.strip().lower() if ml_mode else None,
        "ml_istart": samp("ML_ISTART", as_int),
        "lepsilon": samp("LEPSILON", as_bool),
        "lcalceps": samp("LCALCEPS", as_bool),
        "lchimag": samp("LCHIMAG", as_bool),
        "algo": samp("ALGO", as_str),
        "ispin": as_int(_known(method.get("ISPIN"))),
        "nupdown": as_float(_known(method.get("NUPDOWN"))),
        "ismear": as_int(_known(method.get("ISMEAR"))),
        "lnoncollinear": as_bool(_known(method.get("LNONCOLLINEAR"))),
        "lsorbit": as_bool(_known(method.get("LSORBIT"))),
        "lhfcalc": as_bool(_known(method.get("LHFCALC"))),
        "ivdw": as_int(_known(method.get("IVDW"))),
        "luse_vdw": as_bool(_known(method.get("LUSE_VDW"))),
        "magmom_parameters": as_float_list(parameters.get("MAGMOM")),
    }
    return facts


def calc_type_of(facts: Mapping[str, Any]) -> tuple[str, str | None]:
    """``(calc_type, detail)``; calc_type ``unsupported`` marks response/DFPT runs (never labels)."""
    ibrion, nsw = facts.get("ibrion"), facts.get("nsw")
    response = [tag for tag in ("lepsilon", "lcalceps", "lchimag") if facts.get(tag) is True]
    if ibrion in {5, 6, 7, 8} or response:
        what = [f"IBRION={ibrion}"] if ibrion in {5, 6, 7, 8} else []
        what += [f"{tag.upper()}=T" for tag in response]
        return "unsupported", "linear-response/finite-difference calculation (" + ", ".join(what) + \
            "): trailing scsteps are response iterations, frames are not energy/force labels"
    if ibrion == -1 or nsw == 0:
        return "static", None
    if ibrion in {1, 2, 3}:
        return "relaxation", None
    if ibrion == 0:
        return "md", None
    if ibrion is None:
        return "other", f"IBRION unknown (NSW={nsw})"
    return "other", f"IBRION={ibrion}"


def mlff_mode(facts: Mapping[str, Any]) -> tuple[bool, str | None]:
    """``(active, mode)``: mode in train/run/select/refit/refitbayesian/delta/unknown (None when inactive)."""
    mode = facts.get("ml_mode")
    lmlff = facts.get("ml_lmlff")
    if mode in {"none", ""}:
        mode = None
    active = lmlff is True or mode is not None
    if not active:
        return False, None
    if mode is None:
        mode = _ML_ISTART_MODE.get(facts.get("ml_istart"), "unknown") if facts.get("ml_istart") is not None \
            else "unknown"
    if mode not in {"train", "run", "select", "refit", "refitbayesian", "delta"}:
        mode = "unknown"
    return True, mode


# --------------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------------

def interplanar_spacings(cell: Any) -> np.ndarray:
    """Distance between successive lattice planes (s_i = const) per axis: V/|a_j x a_k| = 1/|b_i|."""
    cell = np.asarray(cell, dtype=np.float64)
    reciprocal = np.linalg.inv(cell).T  # rows b_i with a_i . b_j = delta_ij (no 2 pi)
    return 1.0 / np.linalg.norm(reciprocal, axis=1)


def vacuum_axes(cell: Any, fractional: Any, *, threshold: float = 5.0,
                pbc: Sequence[bool] = (True, True, True)) -> tuple[list[int], list[float]]:
    """Lattice axes with an atom-free slab of at least ``threshold`` A, and the gap per axis.

    For axis i the atoms' wrapped fractional coordinates s_i are sorted on a
    circle; the largest gap (fraction) times the interplanar spacing d_i is the
    thickness of the empty slab bounded by lattice planes (a_j, a_k). This is
    exact for skewed cells (unlike a gap along the lattice vector itself).
    Non-periodic axes are reported as vacuum axes by definition.
    """
    frac = np.asarray(fractional, dtype=np.float64)
    spacing = interplanar_spacings(cell)
    axes: list[int] = []
    gaps: list[float] = []
    for axis in range(3):
        if not pbc[axis]:
            axes.append(axis)
            gaps.append(math.inf)
            continue
        if frac.shape[0] == 0:
            gaps.append(float(spacing[axis]))
            axes.append(axis)
            continue
        s = np.sort(np.mod(frac[:, axis], 1.0))
        circular = np.append(np.diff(s), 1.0 - s[-1] + s[0])
        gap = float(circular.max() * spacing[axis])
        gaps.append(gap)
        if gap >= threshold:
            axes.append(axis)
    return axes, gaps


def _finite(value: Any) -> bool:
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError):
        return False
    return bool(np.all(np.isfinite(array)))


# --------------------------------------------------------------------------
# Step descriptors (memory-light view used by the run-level checks)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class StepDescriptor:
    """What the run-level checks need from one ionic step (no arrays)."""

    index: int
    label_source: str
    complete: bool
    n_scf: int
    free_energy: float | None
    additive_correction: float | None
    force_abs_sum: float | None
    has_forces: bool


def describe_step(step: Any, *, pstress: float | None, version: tuple[int, int, int] | None = None) -> StepDescriptor:
    energies = _energies(step, pstress, version)
    forces = _get(step, "forces")
    force_sum = None
    if forces is not None:
        array = np.asarray(forces, dtype=np.float64)
        force_sum = float(np.abs(array).sum()) if array.size else 0.0
    return StepDescriptor(
        index=int(step.index), label_source=step.label_source, complete=bool(step.complete),
        n_scf=len(step.scf_energies or []), free_energy=energies.get("free_energy"),
        additive_correction=energies.get("additive_correction"), force_abs_sum=force_sum,
        has_forces=forces is not None,
    )


def _energies(step: Any, pstress: float | None, version: tuple[int, int, int] | None) -> dict[str, Any]:
    """vasprun.derived_energies (the parser owns the energy rule and its provenance)."""
    return dict(_vasprun.derived_energies(step, pstress, version))


def _as_descriptor(step: Any, pstress, version) -> StepDescriptor:
    return step if isinstance(step, StepDescriptor) else describe_step(step, pstress=pstress, version=version)


# --------------------------------------------------------------------------
# Species consistency
# --------------------------------------------------------------------------

def species_findings(header: Any, *, potcar: Any = None, poscar: Any = None, outcar: Any = None) -> list[str]:
    """Human-readable disagreements between vasprun species and the POTCAR/POSCAR/OUTCAR evidence."""
    problems: list[str] = []
    atom_types = list(_get(header, "atom_types") or [])
    types = [(str(t.element).strip(), int(t.count)) for t in atom_types]
    species = [str(s).strip() for s in (_get(header, "species") or [])]
    expanded = [element for element, count in types for _ in range(count)]
    if expanded and species and expanded != species:
        problems.append("vasprun <atominfo> atoms do not follow the atomtypes blocks")
    elements = [element for element, _ in types]
    counts = [count for _, count in types]
    titel_elements = [titel_element(t.pseudopotential) for t in atom_types]
    for index, (element, from_titel) in enumerate(zip(elements, titel_elements)):
        if from_titel is not None and from_titel != element:
            problems.append(f"vasprun type {index}: element {element} but pseudopotential {atom_types[index].pseudopotential!r}")
    for problem in _get(header, "problems") or []:
        if str(problem).startswith("species_mismatch"):
            problems.append(f"vasprun: {problem}")
    if potcar is not None:
        datasets = [d for d in (_get(potcar, "datasets") or []) if isinstance(d, Mapping)]
        potcar_elements = [as_str(d.get("element")) or titel_element(as_str(d.get("titel"))) for d in datasets]
        if potcar_elements and elements and potcar_elements != elements:
            problems.append(f"POTCAR file element order {potcar_elements} != vasprun {elements}")
    if poscar is not None:
        poscar_species = _get(poscar, "species")
        poscar_counts = list(_get(poscar, "counts") or [])
        if poscar_species is not None and elements and list(poscar_species) != elements:
            problems.append(f"POSCAR species {list(poscar_species)} != vasprun {elements}")
        if poscar_counts and counts and poscar_counts != counts:
            problems.append(f"POSCAR counts {poscar_counts} != vasprun {counts}")
    if outcar is not None:
        ions = _get(outcar, "ions_per_type")
        if ions is not None and counts and list(ions) != counts:
            problems.append(f"OUTCAR ions per type {list(ions)} != vasprun {counts}")
        titles = list(_get(outcar, "potcar_titles") or [])
        outcar_elements = [titel_element(t) for t in titles]
        if outcar_elements and elements and outcar_elements != elements:
            problems.append(f"OUTCAR POTCAR order {outcar_elements} != vasprun {elements}")
        for problem in _get(outcar, "problems") or []:
            if str(problem).startswith("species_mismatch"):
                problems.append(f"OUTCAR: {problem}")
    return problems


# --------------------------------------------------------------------------
# Magnetism
# --------------------------------------------------------------------------

def magnetic_species_of(species: Sequence[str], magmom: Sequence[float] | None, explicit: bool | None,
                        policy: Policy) -> list[str]:
    """Configured magnetic species present in the run, plus elements with explicit |MAGMOM| >= threshold."""
    present = set(species)
    chosen = present & set(policy.magnetic_species)
    if explicit and magmom is not None and len(magmom) == len(species):
        for element, moment in zip(species, magmom):
            if abs(float(moment)) >= policy.magmom_species_threshold:
                chosen.add(element)
    return sorted(chosen)


def sign_pattern(moments: Sequence[float], sites: Sequence[int], threshold: float) -> tuple[int, ...]:
    """Signs (+1/-1, 0 = |m| < threshold) on ``sites``, canonical under a global spin flip."""
    signs = [0 if abs(float(moments[i])) < threshold else (1 if moments[i] > 0 else -1) for i in sites]
    first = next((s for s in signs if s != 0), 1)
    return tuple(s * first for s in signs)


def _state_id(pattern: tuple[int, ...], ispin: int | None, nupdown: float | None) -> str:
    payload = json.dumps({"pattern": list(pattern), "ispin": ispin, "nupdown": nupdown}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


def magnetization_evidence(
    dft_indices: Sequence[int],
    all_indices: Sequence[int],
    *,
    outcar: Any = None,
    outcar_map: Mapping[int, int] | None = None,
    oszicar: Any = None,
    oszicar_map: Mapping[int, int] | None = None,
    n_atoms: int | None = None,
    complete: bool = True,
) -> tuple[dict[int, tuple[float, str]], dict[int, tuple[list[float], str]], list[str]]:
    """Per DFT step total moment ``{index: (value, source)}`` and site moments ``{index: (list, source)}``.

    Total: the OUTCAR per-step value (last ``number of electron ...
    magnetization`` of the electronic loop) where mapped, else OSZICAR
    ``mag=``. Sites: the per-step OUTCAR ``magnetization (x)`` table (LORBIT >=
    10), else the table VASP prints after the final step, attached to the last
    DFT step of a complete run only.
    """
    totals: dict[int, tuple[float, str]] = {}
    sites: dict[int, tuple[list[float], str]] = {}
    flags: list[str] = []
    outcar_totals = list(_get(outcar, "total_magnetizations") or []) if outcar is not None else []
    tables = dict(_get(outcar, "magnetization_tables") or {}) if outcar is not None else {}
    ionic = list(_get(oszicar, "ionic_lines") or []) if oszicar is not None else []
    for index in dft_indices:
        value, source = None, None
        if outcar_map is not None and index in outcar_map and outcar_map[index] < len(outcar_totals):
            candidate = outcar_totals[outcar_map[index]]
            if candidate is not None and math.isfinite(candidate):
                value, source = float(candidate), "outcar"
        if value is None and oszicar_map is not None and index in oszicar_map and oszicar_map[index] < len(ionic):
            candidate = ionic[oszicar_map[index]].get("mag")
            if candidate is not None and math.isfinite(candidate):
                value, source = float(candidate), "oszicar"
        if value is not None:
            totals[index] = (value, source)
        if outcar_map is not None and index in outcar_map and outcar_map[index] in tables:
            table = [float(v) for v in tables[outcar_map[index]]]
            if n_atoms is None or len(table) == n_atoms:
                sites[index] = (table, "outcar_table_step")
            else:
                flags.append("site_moment_table_wrong_length")
    final_table = _get(outcar, "final_magnetization_table") if outcar is not None else None
    if not sites and final_table is not None and dft_indices and complete and outcar_map is not None:
        table = [float(v) for v in final_table]
        if n_atoms is None or len(table) == n_atoms:
            sites[dft_indices[-1]] = (table, "outcar_table_final_only")
        else:
            flags.append("site_moment_table_wrong_length")
    return totals, sites, sorted(set(flags))


@dataclass
class MagneticAssessment:
    """Run-level magnetic class plus one record per DFT step (``frames[index]``)."""

    run: dict[str, Any]
    frames: dict[int, dict[str, Any]]

    def frame(self, index: int) -> dict[str, Any]:
        return self.frames.get(index, {"magnetic_class": self.run.get("magnetic_class"),
                                       "policy_result": "not_applicable"})

    def as_dict(self) -> dict[str, Any]:
        return {"run": json_safe(self.run), "frames": {str(k): json_safe(v) for k, v in sorted(self.frames.items())}}


def classify_magnetism(
    *,
    species: Sequence[str],
    ispin: int | None,
    nupdown: float | None,
    noncollinear: bool | None,
    magmom_initial: Sequence[float] | None,
    magmom_explicit: bool | None,
    dft_indices: Sequence[int],
    totals: Mapping[int, tuple[float, str]],
    sites: Mapping[int, tuple[Sequence[float], str]],
    policy: Policy = DEFAULT_POLICY,
    istart: int | None = None,
    magmom_source: str | None = None,
) -> MagneticAssessment:
    """Magnetic class of a run and the per-frame magnetic record + policy result.

    Transitions (each starts a new ``magnetic_segment``): total moment jump
    vs the previous DFT frame > ``mag_jump_threshold``; first crossing of
    |total - total_first| > ``mag_drift_threshold``; a change of the canonical
    site sign pattern vs the previous frame; the first frame whose pattern does
    not match the explicit MAGMOM pattern (``initial_order_mismatch``); for a
    balanced (antiferromagnetic, sum 0) explicit MAGMOM, a first total moment
    further than ``mag_initial_total_tolerance`` from 0.
    """
    species = [str(s) for s in species]
    n_atoms = len(species)
    fixed_nupdown = nupdown is not None and nupdown >= 0
    magmom = [float(m) for m in magmom_initial] if magmom_initial is not None else None
    if magmom is not None and len(magmom) != n_atoms:
        magmom_usable = None
    else:
        magmom_usable = magmom
    mag_species = magnetic_species_of(species, magmom_usable, magmom_explicit, policy) if ispin == 2 else []
    mag_sites = [i for i, element in enumerate(species) if element in mag_species]
    thr = policy.magmom_sign_threshold
    initial_pattern = None
    if magmom_explicit and magmom_usable is not None and mag_sites and ispin == 2:
        initial_pattern = sign_pattern(magmom_usable, mag_sites, thr)
    balanced = bool(
        magmom_explicit and magmom_usable is not None and any(abs(m) >= thr for m in magmom_usable)
        and abs(sum(magmom_usable)) < 1e-6
    )
    run: dict[str, Any] = {
        "ispin": ispin,
        "nupdown": nupdown,
        "nupdown_fixed": fixed_nupdown,
        "noncollinear": noncollinear,
        "magmom_explicit": magmom_explicit,
        "magmom_source": magmom_source,
        "magmom_initial": magmom if ispin == 2 else None,
        "magmom_initial_sum": float(sum(magmom)) if magmom is not None and ispin == 2 else None,
        "magmom_initial_pattern": list(initial_pattern) if initial_pattern is not None else None,
        "magnetic_species": mag_species,
        "istart": istart,
        "thresholds": {
            "magmom_sign_threshold": thr, "mag_jump_threshold": policy.mag_jump_threshold,
            "mag_drift_threshold": policy.mag_drift_threshold,
            "mag_initial_total_tolerance": policy.mag_initial_total_tolerance,
            "magmom_species_threshold": policy.magmom_species_threshold,
        },
        "flags": [],
    }
    if ispin == 2 and istart is not None and istart >= 1:
        run["flags"].append("magnetization_initialized_from_wavecar")

    frames: dict[int, dict[str, Any]] = {}
    transitions: list[dict[str, Any]] = []
    first_total = None
    prev_total = None
    first_sites = None
    prev_sites = None
    prev_pattern = None
    drift_crossed = False
    segment = 0
    n_total = n_site = 0
    for position, index in enumerate(dft_indices):
        total, total_source = totals.get(index, (None, None))
        site_values, site_source = sites.get(index, (None, None))
        site_list = [float(v) for v in site_values] if site_values is not None else None
        record: dict[str, Any] = {
            "total": total, "total_source": total_source,
            "site_moments": site_list, "site_source": site_source,
            "d_total_first": None, "d_total_previous": None,
            "max_d_site_first": None, "max_d_site_previous": None,
            "sign_pattern": None, "pattern_matches_initial": None, "magnetic_state_id": None,
            "transition": [],
        }
        kinds: list[str] = []
        if total is not None:
            n_total += 1
            if first_total is None:
                first_total = total
                if balanced and abs(total) > policy.mag_initial_total_tolerance:
                    kinds.append("afm_total_moment_mismatch")
            record["d_total_first"] = total - first_total
            if prev_total is not None:
                record["d_total_previous"] = total - prev_total
                if abs(total - prev_total) > policy.mag_jump_threshold:
                    kinds.append("total_moment_jump")
            if not drift_crossed and abs(total - first_total) > policy.mag_drift_threshold:
                drift_crossed = True
                if "total_moment_jump" not in kinds:
                    kinds.append("total_moment_drift")
            prev_total = total
        if site_list is not None:
            n_site += 1
            array = np.asarray(site_list)
            if first_sites is None:
                first_sites = array
            record["max_d_site_first"] = float(np.max(np.abs(array - first_sites))) if array.size else 0.0
            if prev_sites is not None and prev_sites.shape == array.shape:
                record["max_d_site_previous"] = float(np.max(np.abs(array - prev_sites))) if array.size else 0.0
            prev_sites = array
            if mag_sites:
                pattern = sign_pattern(site_list, mag_sites, thr)
                record["sign_pattern"] = list(pattern)
                record["magnetic_state_id"] = _state_id(pattern, ispin, nupdown)
                if initial_pattern is not None:
                    record["pattern_matches_initial"] = pattern == initial_pattern
                if prev_pattern is not None and pattern != prev_pattern:
                    kinds.append("site_pattern_change")
                elif prev_pattern is None and initial_pattern is not None and pattern != initial_pattern:
                    kinds.append("initial_order_mismatch")
                prev_pattern = pattern
        if kinds and ispin == 2:
            segment += 1  # segment 0 = the initialized state; every transition opens a new branch
            transitions.append({"index": index, "kinds": kinds})
        record["transition"] = kinds if ispin == 2 else []
        record["segment"] = segment if ispin == 2 else 0
        frames[index] = record

    # class ------------------------------------------------------------------------------
    if ispin is None:
        magnetic_class, why = MAG_UNKNOWN, "ISPIN unknown"
    elif ispin == 1:
        magnetic_class, why = MAG_NON_SPIN_POLARIZED, "ISPIN=1"
    elif not (magmom_explicit or fixed_nupdown):
        if magmom_explicit is None:
            magnetic_class, why = MAG_UNKNOWN, "ISPIN=2 and it is unknown whether MAGMOM was set (no <incar>/INCAR)"
        else:
            magnetic_class, why = MAG_UNCONTROLLED, "ISPIN=2 without explicit MAGMOM or fixed NUPDOWN (VASP default FM start)"
    elif n_total == 0 and n_site == 0:
        magnetic_class, why = MAG_UNKNOWN, "ISPIN=2 but no per-frame magnetization evidence (OUTCAR/OSZICAR)"
    elif transitions:
        magnetic_class, why = MAG_CONTROLLED_TRANSITION, f"{len(transitions)} magnetic transition(s)"
    else:
        magnetic_class, why = MAG_CONTROLLED_CONSISTENT, "explicit initialization, no state change observed"
    if ispin == 2 and noncollinear:
        run["flags"].append("noncollinear_moments_not_resolved")
    if ispin == 2 and n_site == 0 and n_total > 0:
        run["flags"].append("magnetic_state_unresolved")  # total moment only (no LORBIT table)
    first_transition = transitions[0]["index"] if transitions else None
    run.update({
        "magnetic_class": magnetic_class,
        "class_reason": why,
        "n_frames_with_total": n_total,
        "n_frames_with_site_moments": n_site,
        "transitions": transitions,
        "first_transition_index": first_transition,
        "overrides": policy.magnetic_overrides,
        "override_reason": policy.magnetic_override_reason,
    })

    # policy ------------------------------------------------------------------------------
    has_magnetic_species = bool(mag_species)
    run_reason = None
    if magnetic_class == MAG_UNCONTROLLED and has_magnetic_species:
        run_reason = None if policy.accept_uncontrolled else "magnetic_uncontrolled"
    elif magnetic_class == MAG_UNKNOWN and has_magnetic_species:
        run_reason = None if policy.accept_unknown else "magnetic_unknown"
    run_overrides = []
    if magnetic_class == MAG_UNCONTROLLED and has_magnetic_species and policy.accept_uncontrolled:
        run_overrides.append("accept_uncontrolled")
    if magnetic_class == MAG_UNKNOWN and has_magnetic_species and policy.accept_unknown:
        run_overrides.append("accept_unknown")
    first_kinds = transitions[0]["kinds"] if transitions else []
    transition_reason = (
        "magnetic_order_changed"
        if any(k in {"site_pattern_change", "initial_order_mismatch"} for k in first_kinds)
        else "magnetic_state_change"
    )
    counts: dict[str, int] = {}
    first_position = list(dft_indices).index(first_transition) if first_transition is not None else None
    for position, index in enumerate(dft_indices):
        record = frames[index]
        overrides = list(run_overrides)
        reason = run_reason
        if reason is None and ispin == 2 and first_position is not None and first_position <= position:
            if policy.accept_transitions:
                overrides.append("accept_transitions")
            else:
                reason = transition_reason
        if reason is None and ispin == 2 and has_magnetic_species and record["total"] is None \
                and record["site_moments"] is None and magnetic_class != MAG_UNKNOWN:
            if policy.accept_unknown:
                overrides.append("accept_unknown")
            else:
                reason = "magnetic_unknown"
        if reason is not None:
            result = f"quarantined:{reason}"
        elif overrides:
            result = "accepted_by_override:" + "+".join(sorted(set(overrides)))
        else:
            result = "accepted"
        record["magnetic_class"] = magnetic_class
        record["policy_result"] = result
        record["policy_reason"] = reason
        counts[result] = counts.get(result, 0) + 1
    run["policy_result"] = (
        f"quarantined:{run_reason}" if run_reason
        else ("accepted_by_override:" + "+".join(run_overrides)) if run_overrides
        else "partially_quarantined" if any(r.startswith("quarantined") for r in counts)
        else "accepted"
    )
    run["frame_policy_counts"] = dict(sorted(counts.items()))
    run["flags"] = sorted(set(run["flags"]))
    return MagneticAssessment(run=run, frames=frames)


def _magmom_explicit(header: Any, incar: Any) -> tuple[bool | None, str | None, list[str]]:
    """Whether MAGMOM was set by the user: vasprun <incar> (executed) first, then the INCAR file."""
    flags: list[str] = []
    header_incar = {str(k).upper() for k in (_get(header, "incar") or {})}
    file_tags = {str(k).upper() for k in ((_get(incar, "tags") or {}) if incar is not None else {})}
    if "MAGMOM" in header_incar:
        return True, "vasprun_incar", flags
    if header_incar:
        if "MAGMOM" in file_tags:
            flags.append("magmom_incar_file_differs")
        return False, "vasprun_incar", flags
    if incar is not None:
        return ("MAGMOM" in file_tags), "incar_file", flags
    return None, None, flags


# --------------------------------------------------------------------------
# Run assessment
# --------------------------------------------------------------------------

@dataclass
class RunAssessment:
    run_id: str
    outcome: Outcome
    policy: Policy
    calc_type: str | None = None
    calc_detail: str | None = None
    facts: dict[str, Any] = field(default_factory=dict)
    settings: dict[str, Any] = field(default_factory=dict)
    species: list[str] = field(default_factory=list)
    n_atoms: int | None = None
    findings: list[dict[str, str]] = field(default_factory=list)  # every run-level reason found (precedence order)
    flags: list[str] = field(default_factory=list)
    completion: dict[str, Any] = field(default_factory=dict)
    mlff: dict[str, Any] = field(default_factory=dict)
    evidence: dict[str, Any] = field(default_factory=dict)
    outcar_map: dict[int, int] | None = None  # vasprun step index -> OUTCAR DFT-step ordinal
    oszicar_map: dict[int, int] | None = None  # vasprun step index -> OSZICAR ionic line
    selective: Any = None  # (N, 3) bool, DIRECT basis, True = may move
    selective_source: str | None = None
    magnetic: MagneticAssessment | None = None
    outcar: Any = None
    oszicar: Any = None
    n_steps: int = 0
    n_dft_steps: int = 0

    @property
    def accepted(self) -> bool:
        return self.outcome.status == ACCEPTED

    def as_dict(self) -> dict[str, Any]:
        return json_safe({
            "run_id": self.run_id,
            "outcome": self.outcome.as_dict(),
            "calc_type": self.calc_type,
            "calc_detail": self.calc_detail,
            "facts": self.facts,
            "n_atoms": self.n_atoms,
            "n_steps": self.n_steps,
            "n_dft_steps": self.n_dft_steps,
            "findings": self.findings,
            "flags": self.flags,
            "completion": self.completion,
            "mlff": self.mlff,
            "evidence": self.evidence,
            "selective": {
                "present": self.selective is not None,
                "source": self.selective_source,
                "basis": "direct" if self.selective is not None else None,
                "n_fixed_components": int((~np.asarray(self.selective)).sum()) if self.selective is not None else 0,
            },
            "magnetic": self.magnetic.run if self.magnetic is not None else None,
        })


def _map_counts(n_evidence: int, indices_all: Sequence[int], indices_complete: Sequence[int]) -> dict[int, int] | None:
    if n_evidence == len(indices_all):
        return {index: k for k, index in enumerate(indices_all)}
    if n_evidence == len(indices_complete):
        return {index: k for k, index in enumerate(indices_complete)}
    return None


def classify_run(
    *,
    run_id: str,
    header: Any,
    trailer: Any = None,
    steps: Sequence[Any] = (),
    settings: Mapping[str, Any] | None = None,
    outcar: Any = None,
    oszicar: Any = None,
    potcar: Any = None,
    poscar: Any = None,
    incar: Any = None,
    policy: Policy | None = None,
    reference_hashes: Iterable[str] = (),
    header_error: str | None = None,
) -> RunAssessment:
    """Classify one run after its vasprun steps have been read.

    ``steps``: the run's :class:`~nio_md_prep.dataset.model.IonicStep` objects
    or :class:`StepDescriptor` summaries (file order, all steps incl. MLFF and
    an incomplete tail). ``reference_hashes``: sha256 digests of files copied
    from reference/template directories (agglomeration ``reference_files`` /
    ``shared_files``). Evidence objects are the parsers' OutcarEvidence,
    OszicarEvidence, PotcarEvidence, PoscarEvidence and IncarEvidence.
    """
    policy = policy or DEFAULT_POLICY
    if header is None:
        return RunAssessment(run_id=run_id, policy=policy,
                             outcome=outcome_for("vasprun_unreadable", header_error or "vasprun.xml header unreadable"),
                             findings=[{"reason": "vasprun_unreadable", "detail": header_error or ""}])
    if settings is None:
        settings = extract_settings(header, incar_file=incar, outcar=outcar, potcar=potcar)
    settings = dict(settings)
    facts = run_facts(header, settings)
    species = [str(s).strip() for s in (_get(header, "species") or [])]
    n_atoms = len(species) or None
    assessment = RunAssessment(run_id=run_id, outcome=Outcome(ACCEPTED), policy=policy, facts=facts,
                               settings=settings, species=species, n_atoms=n_atoms)
    findings: list[dict[str, str]] = []
    flags: list[str] = []

    def find(reason: str, detail: str) -> None:
        findings.append({"reason": reason, "detail": detail})

    # 1. program -------------------------------------------------------------------------------
    if facts["program"] is not None and not facts["program"].startswith("vasp"):
        find("not_vasp", f"generator program {facts['program']!r}")

    # 2. calculation type --------------------------------------------------------------------
    calc_type, detail = calc_type_of(facts)
    assessment.calc_type, assessment.calc_detail = calc_type, detail
    if calc_type == "unsupported":
        find("unsupported_calculation", detail or "")
    elif calc_type not in policy.allowed_calc_types:
        find("unsupported_calculation",
             f"calc_type {calc_type!r} ({detail or 'see facts'}) is not in allowed_calc_types {list(policy.allowed_calc_types)}")

    # 3. copied reference outputs --------------------------------------------------------------
    references = {str(h).lower() for h in reference_hashes if h}
    label_sha = (_get(trailer, "source_sha256") or "").lower() if trailer is not None else ""
    if label_sha and label_sha in references:
        find("copied_reference_output", "vasprun.xml is byte-identical to a reference/template file")
    ignored_evidence = []
    if outcar is not None and str(_get(outcar, "sha256") or "").lower() in references:
        ignored_evidence.append("OUTCAR")
        outcar = None
    if oszicar is not None and str(_get(oszicar, "sha256") or "").lower() in references:
        ignored_evidence.append("OSZICAR")
        oszicar = None
    for name in ignored_evidence:
        flags.append(f"copied_reference_evidence_ignored:{name}")
    assessment.outcar, assessment.oszicar = outcar, oszicar

    # step bookkeeping -----------------------------------------------------------------------------
    version = facts["version_tuple"]
    descriptors = [_as_descriptor(step, facts["pstress"], version) for step in steps]
    all_indices = [d.index for d in descriptors]
    complete_indices = [d.index for d in descriptors if d.complete]
    dft = [d for d in descriptors if d.label_source == "dft"]
    dft_indices = [d.index for d in dft]
    dft_complete = [d.index for d in dft if d.complete]
    ml_indices = [d.index for d in descriptors if d.label_source != "dft"]
    assessment.n_steps, assessment.n_dft_steps = len(descriptors), len(dft)

    # 4. MLFF ---------------------------------------------------------------------------------------
    active, mode = mlff_mode(facts)
    mlff = {"active": active, "mode": mode, "n_mlff_steps": len(ml_indices), "n_dft_steps": len(dft),
            "outcar_ml_blocks": _get(outcar, "ml_steps") if outcar is not None else None}
    if ml_indices and not active:
        flags.append("mlff_steps_without_ml_settings")
    if active:
        if mode in _ML_NO_DFT:
            find("mlff_steps", f"ML_MODE={mode}: no ab initio labels are computed")
        elif mode == "delta":
            find("mlff_steps", "ML_MODE=delta: labels are DFT + MLFF sums, not pure DFT")
        elif mode == "unknown":
            find("mlff_steps", "VASP-MLFF active but ML_MODE/ML_ISTART undeterminable: DFT steps cannot be trusted per step")
    assessment.mlff = mlff

    # 5. species ------------------------------------------------------------------------------------
    for problem in species_findings(header, potcar=potcar, poscar=poscar, outcar=outcar):
        find("species_mismatch", problem)

    # 6. evidence --------------------------------------------------------------------------------
    evidence: dict[str, Any] = {"ignored": ignored_evidence, "mismatches": []}
    mismatches: list[str] = evidence["mismatches"]
    if outcar is not None:
        nions = _get(outcar, "nions")
        if nions is not None and n_atoms is not None and nions != n_atoms:
            mismatches.append(f"OUTCAR NIONS={nions} != vasprun {n_atoms} atoms")
        for conflict in outcar_conflicts(settings):
            mismatches.append(
                f"OUTCAR {conflict['field']}={conflict['others'].get('outcar')!r} != vasprun {conflict['value']!r}"
            )
        n_out = len(_get(outcar, "free_energies") or [])
        assessment.outcar_map = _map_counts(n_out, dft_indices, dft_complete)
        evidence["outcar_dft_steps"] = n_out
        if assessment.outcar_map is None:
            flags.append("evidence_count_mismatch:outcar")
        else:
            totens = list(_get(outcar, "free_energies") or [])
            force_sums = list(_get(outcar, "force_abs_sums") or [])
            by_index = {d.index: d for d in dft}
            for index, k in sorted(assessment.outcar_map.items()):
                desc = by_index[index]
                if not desc.complete:
                    continue
                toten = totens[k] if k < len(totens) else None
                if toten is not None and desc.free_energy is not None and math.isfinite(toten) \
                        and math.isfinite(desc.free_energy) \
                        and abs(toten - desc.free_energy) > policy.outcar_energy_tolerance:
                    mismatches.append(f"step {index}: OUTCAR TOTEN {toten!r} != vasprun F {desc.free_energy!r}")
                fsum = force_sums[k] if k < len(force_sums) else None
                if fsum is not None and desc.force_abs_sum is not None and math.isfinite(fsum) and n_atoms:
                    tolerance = 3 * n_atoms * policy.outcar_force_tolerance + 1e-9
                    if abs(fsum - desc.force_abs_sum) > tolerance:
                        mismatches.append(f"step {index}: OUTCAR sum|F| {fsum:.6f} != vasprun {desc.force_abs_sum:.6f}")
    if oszicar is not None:
        ionic = list(_get(oszicar, "ionic_lines") or [])
        assessment.oszicar_map = _map_counts(len(ionic), all_indices, complete_indices)
        evidence["oszicar_ionic_steps"] = len(ionic)
        if assessment.oszicar_map is None:
            flags.append("evidence_count_mismatch:oszicar")
        else:
            by_index = {d.index: d for d in descriptors}
            for index, k in sorted(assessment.oszicar_map.items()):
                desc = by_index[index]
                value = ionic[k].get("F")
                if value is None or desc.free_energy is None or not desc.complete:
                    continue
                if not (math.isfinite(value) and math.isfinite(desc.free_energy)):
                    continue
                tolerance = policy.oszicar_relative_tolerance * max(abs(desc.free_energy), 1.0) + 1e-8
                if abs(value - desc.free_energy) > tolerance:
                    mismatches.append(f"step {index}: OSZICAR F {value!r} != vasprun F {desc.free_energy!r}")
    # selective dynamics: vasprun initialpos, else finalpos, else POSCAR; any disagreement is a mismatch
    flags_vasprun, source, conflict = _vasprun.selective_flags(header, trailer)
    poscar_flags = _get(poscar, "selective") if poscar is not None else None
    if conflict:
        mismatches.append("selective-dynamics flags differ between initialpos and finalpos")
    if flags_vasprun is not None and poscar_flags is not None:
        a, b = np.asarray(flags_vasprun, dtype=bool), np.asarray(poscar_flags, dtype=bool)
        if a.shape != b.shape or not np.array_equal(a, b):
            mismatches.append(f"selective-dynamics flags of vasprun {source} differ from POSCAR")
    if flags_vasprun is not None:
        assessment.selective, assessment.selective_source = np.asarray(flags_vasprun, dtype=bool), f"vasprun:{source}"
    elif poscar_flags is not None:
        assessment.selective, assessment.selective_source = np.asarray(poscar_flags, dtype=bool), "POSCAR"
    if assessment.selective is not None and n_atoms is not None and assessment.selective.shape != (n_atoms, 3):
        mismatches.append(f"selective-dynamics flags have shape {assessment.selective.shape} for {n_atoms} atoms")
    for problem in mismatches:
        find("evidence_mismatch", problem)
    assessment.evidence = evidence

    # 7. completion -----------------------------------------------------------------------------------
    closed = bool(_get(trailer, "closed")) if trailer is not None else False
    truncated = bool(_get(trailer, "truncated")) if trailer is not None else True
    outcar_completed = _get(outcar, "completed") if outcar is not None else None
    complete = closed and not truncated and (outcar is None or bool(outcar_completed))
    assessment.completion = {
        "complete": complete, "vasprun_closed": closed, "vasprun_truncated": truncated,
        "outcar_present": outcar is not None, "outcar_completed": outcar_completed,
        "ionic_converged": _get(outcar, "ionic_converged") if outcar is not None else None,
        "n_complete_steps": len(complete_indices), "n_steps": len(descriptors),
        "policy": policy.interrupted_run_policy,
        "trailer_error": _get(trailer, "error") if trailer is not None else None,
    }
    if not complete:
        why = []
        if not closed:
            why.append("vasprun.xml not closed")
        if truncated:
            why.append("vasprun.xml truncated")
        if outcar is not None and not outcar_completed:
            why.append("OUTCAR has no 'General timing' block")
        if policy.interrupted_run_policy == "exclude":
            find("run_incomplete", "; ".join(why))
        else:
            flags.append("recovered_incomplete_run")

    # other run flags ---------------------------------------------------------------------------------
    if facts.get("ismear") == -5:
        flags.append("forces_not_variational")
    if facts.get("ediff") == 0.0:
        flags.append("scf_criterion_disabled")
    corrections = [d.additive_correction for d in dft if d.additive_correction is not None
                   and math.isfinite(d.additive_correction)]
    vdw_known = (facts.get("ivdw") not in (None, 0)) or facts.get("luse_vdw") is True
    if corrections and max(abs(c) for c in corrections) > policy.additive_correction_tolerance and not vdw_known:
        flags.append("additive_energy_correction_unexplained")
    evidence["max_abs_additive_correction"] = max((abs(c) for c in corrections), default=None)

    # magnetism -----------------------------------------------------------------------------------------
    explicit, magmom_source, magmom_flags = _magmom_explicit(header, incar)
    flags.extend(magmom_flags)
    totals, sites, mag_flags = magnetization_evidence(
        dft_indices, all_indices, outcar=outcar, outcar_map=assessment.outcar_map, oszicar=oszicar,
        oszicar_map=assessment.oszicar_map, n_atoms=n_atoms, complete=complete,
    )
    flags.extend(mag_flags)
    assessment.magnetic = classify_magnetism(
        species=species, ispin=facts.get("ispin"), nupdown=facts.get("nupdown"),
        noncollinear=bool(facts.get("lnoncollinear") or facts.get("lsorbit")),
        magmom_initial=facts.get("magmom_parameters"), magmom_explicit=explicit,
        dft_indices=dft_indices, totals=totals, sites=sites, policy=policy,
        istart=facts.get("istart"), magmom_source=magmom_source,
    )
    flags.extend(assessment.magnetic.run.get("flags", []))
    # default NiO policy: an uncontrolled/unknown run with magnetic species is quarantined as a whole
    # (not only its frames), unless the matching override is set (recorded in magnetic.run["overrides"]).
    magnetic_policy = str(assessment.magnetic.run.get("policy_result") or "")
    if magnetic_policy.startswith("quarantined:"):
        find(magnetic_policy.split(":", 1)[1],
             f"magnetic class {assessment.magnetic.run.get('magnetic_class')}: "
             f"{assessment.magnetic.run.get('class_reason')} (magnetic species {assessment.magnetic.run.get('magnetic_species')})")

    # outcome ------------------------------------------------------------------------------------------
    ordered = sorted(findings, key=lambda f: (RUN_REASON_PRECEDENCE.index(f["reason"]), f["detail"]))
    assessment.findings = ordered
    assessment.flags = sorted(set(flags))
    if ordered:
        first = ordered[0]
        same = [f["detail"] for f in ordered if f["reason"] == first["reason"]]
        detail = "; ".join(same[:5]) + (f" (+{len(same) - 5} more)" if len(same) > 5 else "")
        assessment.outcome = outcome_for(first["reason"], detail)
    return assessment


# --------------------------------------------------------------------------
# Frame-level rules
# --------------------------------------------------------------------------

def _problem_codes(step: Any) -> dict[str, list[str]]:
    codes: dict[str, list[str]] = {}
    for problem in _get(step, "problems") or []:
        code, _, detail = str(problem).partition(":")
        codes.setdefault(code.strip(), []).append(detail.strip())
    return codes


def scf_gate(step: Any, run: RunAssessment, policy: Policy | None = None) -> tuple[Outcome | None, dict[str, Any]]:
    """Conservative SCF-convergence decision for one DFT step: ``(outcome or None, scf record)``.

    Returns None as outcome when the step is proven converged.
    """
    policy = policy or run.policy
    facts = run.facts
    nelm, nelmin, ediff = facts.get("nelm"), facts.get("nelmin"), facts.get("ediff")
    version = facts.get("version_tuple")
    vasp6 = version is not None and version[0] >= 6
    summary = _vasprun.scf_summary(step)
    n, last_de = summary["n_steps"], summary["last_dE"]
    record: dict[str, Any] = {
        "n_steps": n, "nelm": nelm, "nelmin": nelmin, "ediff": ediff, "last_dE": last_de,
        "converged": None, "status": "unknown", "evidence": "none", "evidence_detail": [],
        "outcar_mapped": run.outcar_map is not None, "outcar_marker": None,
        "outcar_last_dE": None, "outcar_last_deps": None, "oszicar_scf_steps": None,
        "marker_is_proof": bool(vasp6 and policy.require_outcar_marker_vasp6),
    }
    detail: list[str] = record["evidence_detail"]
    marker = None
    if run.outcar is not None and run.outcar_map is not None and step.index in run.outcar_map:
        k = run.outcar_map[step.index]
        kinds = list(_get(run.outcar, "scf_marker_kinds") or [])
        if not kinds:
            kinds = [{True: "reached", False: "not_reached"}.get(m) for m in (_get(run.outcar, "scf_converged_markers") or [])]
        marker = kinds[k] if k < len(kinds) else None
        changes = list(_get(run.outcar, "last_energy_changes") or [])
        if k < len(changes) and changes[k] is not None:
            record["outcar_last_dE"], record["outcar_last_deps"] = changes[k]
        record["outcar_marker"] = marker
    if run.oszicar is not None and run.oszicar_map is not None and step.index in run.oszicar_map:
        scf_steps = list(_get(run.oszicar, "scf_steps") or [])
        k = run.oszicar_map[step.index]
        if k < len(scf_steps):
            record["oszicar_scf_steps"] = scf_steps[k]
            if scf_steps[k] != n:
                detail.append(f"OSZICAR lists {scf_steps[k]} electronic steps, vasprun {n}")

    def done(outcome: Outcome | None, status: str, evidence: str, why: str) -> tuple[Outcome | None, dict]:
        record["status"] = status
        record["converged"] = {"converged": True, "not_converged": False}.get(status)
        record["evidence"] = evidence
        detail.append(why)
        return outcome, record

    # explicit failure markers are honoured for every version
    if marker in {"not_reached", "hard_stop"}:
        return done(outcome_for("scf_not_converged", f"OUTCAR: aborting loop ({marker})"),
                    "not_converged", "outcar_marker", f"OUTCAR marker {marker}")
    if ediff == 0.0:
        if policy.allow_ediff_zero:
            return done(None, "converged", "scstep_count", "EDIFF=0 (criterion disabled) accepted by allow_ediff_zero")
        return done(outcome_for("scf_convergence_unknown", "EDIFF=0 disables the convergence criterion"),
                    "unknown", "none", "EDIFF=0")
    if nelm is None or ediff is None:
        return done(outcome_for("scf_convergence_unknown", "NELM or EDIFF unknown"), "unknown", "none",
                    "NELM/EDIFF unknown")
    if n == 0:
        return done(outcome_for("scf_convergence_unknown", "no <scstep> in this <calculation>"), "unknown", "none",
                    "no scsteps")
    if nelmin is not None and nelmin > 1 and n < nelmin:
        detail.append(f"n_scf {n} < NELMIN {nelmin} (informational)")
    marker_ok = marker == "reached"
    marker_proof = marker_ok and record["marker_is_proof"]
    marker_needed = record["marker_is_proof"] and run.outcar is not None
    marker_exempt = False
    if marker_needed and marker is None and policy.nwrite_marker_exemption:
        nwrite = facts.get("nwrite")
        first_dft = min(run.outcar_map) if run.outcar_map else None
        if nwrite is not None and nwrite <= 1 and step.index != first_dft:
            marker_exempt = True
    if marker_needed and run.outcar_map is None:
        # OUTCAR present but its per-step markers cannot be mapped to vasprun steps
        marker_needed_unmappable = True
    else:
        marker_needed_unmappable = False

    if n >= nelm:
        if marker_proof and last_de is not None and last_de < ediff:
            return done(None, "converged", "scstep_count+outcar_marker",
                        f"n_scf {n} == NELM but |dE| {last_de:.3g} < EDIFF and VASP>=6 marker 'EDIFF is reached'")
        if last_de is not None and last_de < ediff:
            return done(outcome_for("scf_convergence_unknown",
                                    f"n_scf {n} reached NELM={nelm}; |dE| {last_de:.3g} < EDIFF but no VASP>=6 marker"),
                        "unknown", "scstep_count", "ambiguous at the NELM ceiling")
        return done(outcome_for("scf_not_converged", f"n_scf {n} >= NELM={nelm}"
                                + (f", |dE| {last_de:.3g} >= EDIFF {ediff:g}" if last_de is not None else "")),
                    "not_converged", "scstep_count", "NELM ceiling")
    # n < NELM
    if n < 2 or last_de is None:
        if marker_proof:
            return done(None, "converged", "outcar_marker", f"n_scf {n} < 2; VASP>=6 marker 'EDIFF is reached'")
        return done(outcome_for("scf_convergence_unknown", f"n_scf {n}: |dE| undefined and no VASP>=6 marker"),
                    "unknown", "none", "dE undefined")
    if last_de >= ediff:
        return done(outcome_for("scf_convergence_unknown",
                                f"loop left before NELM ({n} < {nelm}) but |dE| {last_de:.3g} >= EDIFF {ediff:g}"),
                    "unknown", "scstep_count", "early exit without meeting EDIFF")
    # vasprun says converged
    if not marker_needed:
        why = "n_scf < NELM and |dE| < EDIFF"
        if run.outcar is not None and marker_ok:
            why += " (VASP 5/unknown-version marker 'EDIFF is reached' recorded, not proof)"
        return done(None, "converged", "scstep_count", why)
    if marker_ok:
        return done(None, "converged", "scstep_count+outcar_marker", "n_scf < NELM, |dE| < EDIFF, VASP>=6 marker")
    if marker_exempt:
        record["flags"] = ["scf_marker_not_printed"]
        return done(None, "converged", "scstep_count",
                    "n_scf < NELM and |dE| < EDIFF; OUTCAR marker not printed for this step (NWRITE<=1)")
    if marker_needed_unmappable:
        return done(outcome_for("scf_convergence_unknown",
                                "VASP>=6 OUTCAR present but its step count does not match vasprun (markers unusable)"),
                    "unknown", "scstep_count", "OUTCAR markers unmappable")
    return done(outcome_for("scf_convergence_unknown",
                            f"VASP>=6 OUTCAR shows no 'EDIFF is reached' marker for this step (marker={marker!r})"),
                "unknown", "scstep_count", "VASP>=6 marker missing")


@dataclass
class StressAssessment:
    available: bool
    reason: str | None
    source: str | None
    raw_present: bool
    stress: Any = None  # (3, 3) float64, ASE sign, eV/A^3, symmetrized; only when available
    max_asymmetry_kbar: float | None = None
    volume: float | None = None

    def as_dict(self) -> dict[str, Any]:
        return json_safe({"stress_available": self.available, "stress_reason": self.reason,
                          "stress_source": self.source, "stress_raw_present": self.raw_present,
                          "stress_max_asymmetry_kbar": self.max_asymmetry_kbar, "volume": self.volume})


def stress_record(step: Any, facts: Mapping[str, Any], vacuum: Sequence[int], policy: Policy) -> StressAssessment:
    """Whether this step's stress is exported (never zeros: absent stays absent) and why not."""
    raw = _get(step, "stress_kbar_vasp")
    volume = None
    structure = _get(step, "structure")
    if structure is not None and _get(structure, "cell") is not None:
        volume = abs(float(np.linalg.det(np.asarray(structure.cell, dtype=np.float64))))
    if _get(step, "label_source", "dft") != "dft":
        return StressAssessment(False, "not_a_dft_label", None, raw is not None, volume=volume)
    if raw is None:
        return StressAssessment(False, "not_computed", None, False, volume=volume)
    array = np.asarray(raw, dtype=np.float64)
    asymmetry = float(np.max(np.abs(array - array.T)) / 2.0) if array.shape == (3, 3) else None
    isif, pstress = facts.get("isif"), facts.get("pstress")
    reason = None
    if isif is None:
        reason = "isif_unknown"
    elif isif == 1:
        reason = "isif1_trace_only"
    elif pstress is None or (pstress != 0 and not policy.allow_pstress):
        reason = "pstress_nonzero"
    elif vacuum and not policy.allow_vacuum_stress:
        reason = "vacuum_slab"
    elif not policy.include_stress:
        reason = "not_requested"
    if reason is not None:
        return StressAssessment(False, reason, STRESS_SOURCE, True, max_asymmetry_kbar=asymmetry, volume=volume)
    symmetric = 0.5 * (array + array.T)
    return StressAssessment(True, None, STRESS_SOURCE, True, stress=vasp_stress_to_ase(symmetric),
                            max_asymmetry_kbar=asymmetry, volume=volume)


@dataclass
class StepAssessment:
    """Decision and audit record for one ionic step (DFT or MLFF)."""

    index: int
    label_source: str
    complete: bool
    outcome: Outcome
    label_set: str | None  # energy_forces | energy_only | None (not a label)
    energy_quantity: str
    label_energy: float | None
    energy: dict[str, Any]  # vasprun.derived_energies result (values + energy_source/energy_rule/parser)
    scf: dict[str, Any]
    stress: StressAssessment
    magnetic: dict[str, Any]
    vacuum_axes: list[int]
    vacuum_gaps: list[float]
    max_force: float | None
    time_fs: float | None
    md: dict[str, float]
    n_atoms: int | None
    flags: list[str] = field(default_factory=list)
    selective: Any = None  # (N, 3) bool, DIRECT basis (the run's flags), never an ASE constraint

    @property
    def accepted(self) -> bool:
        return self.outcome.status == ACCEPTED

    def as_dict(self) -> dict[str, Any]:
        return json_safe({
            "index": self.index, "label_source": self.label_source, "complete": self.complete,
            "outcome": self.outcome.as_dict(), "label_set": self.label_set,
            "forces_available": self.label_set == LABEL_SET_ENERGY_FORCES,
            "forces_reason": FORCES_REASON_ENERGY_ONLY if self.label_set == LABEL_SET_ENERGY_ONLY else None,
            "energy_quantity": self.energy_quantity, "label_energy": self.label_energy,
            "energy": self.energy, "scf": self.scf, **self.stress.as_dict(), "magnetic": self.magnetic,
            "vacuum_axes": self.vacuum_axes, "vacuum_gaps_A": self.vacuum_gaps, "max_force": self.max_force,
            "time_fs": self.time_fs, "md": self.md, "n_atoms": self.n_atoms, "flags": self.flags,
            "selective_dynamics": self.selective is not None,
        })

    def record_fields(self) -> dict[str, Any]:
        """Values for the existing :class:`~nio_md_prep.dataset.model.FrameRecord` fields."""
        energies = {name: self.energy.get(name) for name in ENERGY_QUANTITIES if self.energy.get(name) is not None}
        return {
            "index": self.index, "n_atoms": self.n_atoms, "time_fs": self.time_fs, "energies": energies,
            "label_energy": self.label_energy, "max_force": self.max_force, "scf": dict(self.scf),
            "stress_available": self.stress.available, "stress_reason": self.stress.reason,
            "total_magnetization": self.magnetic.get("total"),
            "magnetization_source": self.magnetic.get("total_source"),
            "site_magmoms": self.magnetic.get("site_moments"),
            "magnetic_state_id": self.magnetic.get("magnetic_state_id"),
            "magnetic_segment": self.magnetic.get("segment"),
            "label_source": self.label_source, "md": dict(self.md), "flags": list(self.flags),
            "outcome": self.outcome,
        }


_MD_KEYS = {"kinetic": "md_kinetic_energy", "total": "md_total_energy", "lattice kinetic": "md_lattice_kinetic",
            "nosepot": "md_nose_potential", "nosekinetic": "md_nose_kinetic"}


def classify_step(step: Any, run: RunAssessment, policy: Policy | None = None) -> StepAssessment:
    """Frame-level decision (see the module docstring for the rule order)."""
    policy = policy or run.policy
    facts = run.facts
    n_atoms = run.n_atoms
    index = int(step.index)
    flags: list[str] = []
    codes = _problem_codes(step)
    is_dft = step.label_source == "dft"
    energy = _energies(step, facts.get("pstress"), facts.get("version_tuple"))
    flags.extend(energy.get("energy_flags") or [])
    label_energy = energy.get(policy.energy_quantity) if is_dft else None
    structure = _get(step, "structure")
    forces = _get(step, "forces")

    # geometry / vacuum (recorded for every step with a structure)
    axes: list[int] = []
    gaps: list[float] = []
    cell_ok = structure is not None and _get(structure, "cell") is not None \
        and np.shape(structure.cell) == (3, 3) and _finite(structure.cell)
    if cell_ok and _get(structure, "fractional") is not None and _finite(structure.fractional):
        try:
            axes, gaps = vacuum_axes(structure.cell, structure.fractional, threshold=policy.vacuum_gap_threshold)
        except np.linalg.LinAlgError:
            axes, gaps = [], []
    max_force = None
    if forces is not None and np.shape(forces) == (n_atoms, 3) and _finite(forces):
        max_force = float(np.max(np.linalg.norm(np.asarray(forces, dtype=np.float64), axis=1))) if n_atoms else 0.0
    time_fs = None
    if run.calc_type == "md" and facts.get("potim") is not None:
        time_fs = (index + 1) * float(facts["potim"])
    md = {}
    for key, name in _MD_KEYS.items():
        value = (_get(step, "energies") or {}).get(key)
        if value is not None:
            md[name] = value
    if run.oszicar is not None and run.oszicar_map is not None and index in run.oszicar_map:
        line = list(_get(run.oszicar, "ionic_lines") or [])[run.oszicar_map[index]]
        if "T" in line:
            md["md_temperature_K"] = line["T"]
    stress = stress_record(step, facts, axes, policy)
    magnetic = run.magnetic.frame(index) if run.magnetic is not None and is_dft else \
        {"magnetic_class": run.magnetic.run.get("magnetic_class") if run.magnetic else None,
         "policy_result": "not_applicable"}
    if facts.get("ismear") == -5:
        flags.append("forces_not_variational")
    add_corr = energy.get("additive_correction")
    if add_corr is not None and "additive_energy_correction_unexplained" in run.flags and \
            abs(add_corr) > policy.additive_correction_tolerance:
        flags.append("additive_energy_correction_unexplained")

    label_set = LABEL_SET_ENERGY_FORCES if forces is not None else LABEL_SET_ENERGY_ONLY
    scf: dict[str, Any] = {"n_steps": len(step.scf_energies or []), "status": "not_applicable"}

    def result(outcome: Outcome, label: str | None) -> StepAssessment:
        return StepAssessment(
            index=index, label_source=step.label_source, complete=bool(step.complete), outcome=outcome,
            label_set=label, energy_quantity=policy.energy_quantity, label_energy=label_energy, energy=energy,
            scf=scf, stress=stress, magnetic=magnetic, vacuum_axes=axes, vacuum_gaps=gaps, max_force=max_force,
            time_fs=time_fs, md=md, n_atoms=n_atoms, flags=sorted(set(flags + list(scf.get("flags", [])))),
            selective=run.selective,
        )

    # 1. MLFF: counted for indexing/lineage, never a DFT label
    if not is_dft:
        return result(outcome_for("mlff_step", f"VASP-MLFF {step.label_source} step"), None)
    # 1b. runs whose frames are not DFT labels at all (frame checks would be meaningless there)
    if not run.accepted and run.outcome.reason in NON_LABEL_RUN_REASONS:
        return result(outcome_for(run.outcome.reason, f"run-level: {run.outcome.detail}"), None)
    # 2. truncated
    if not step.complete:
        return result(outcome_for("truncated_frame", "file ends inside this ionic step"), None)
    # 3. missing labels
    if label_energy is None:
        why = "; ".join(energy.get("energy_flags") or []) or f"no {policy.energy_quantity}"
        return result(outcome_for("missing_energy", why), None)
    if forces is None and not policy.allow_energy_only:
        return result(outcome_for("missing_forces", "no <varray name=\"forces\"> (energy-only frames need allow_energy_only)"),
                      None)
    # 4. shapes
    shape_problems = list(codes.get("bad_shape", []))
    if structure is None or _get(structure, "positions") is None:
        shape_problems.append("no structure in this step")
    else:
        if np.shape(structure.positions) != (n_atoms, 3):
            shape_problems.append(f"positions shape {np.shape(structure.positions)} for {n_atoms} atoms")
        if np.shape(structure.cell) != (3, 3):
            shape_problems.append(f"cell shape {np.shape(structure.cell)}")
    if forces is not None and np.shape(forces) != (n_atoms, 3):
        shape_problems.append(f"forces shape {np.shape(forces)} for {n_atoms} atoms")
    raw_stress = _get(step, "stress_kbar_vasp")
    if raw_stress is not None and np.shape(raw_stress) != (3, 3):
        shape_problems.append(f"stress shape {np.shape(raw_stress)}")
    if shape_problems:
        return result(outcome_for("bad_shape", "; ".join(shape_problems[:3])), None)
    # 5. non-finite
    bad = list(codes.get("non_finite", []))
    for name, value in (("positions", structure.positions), ("cell", structure.cell), ("forces", forces),
                        ("stress", raw_stress)):
        if value is not None and not _finite(value):
            bad.append(f"{name} not finite")
    for name in ENERGY_QUANTITIES:
        value = energy.get(name)
        if value is not None and not math.isfinite(value):
            bad.append(f"{name} not finite")
    if bad:
        return result(outcome_for("non_finite", "; ".join(bad[:3])), None)
    # 6. cell
    det = float(np.linalg.det(np.asarray(structure.cell, dtype=np.float64)))
    if not det > 1e-6:
        return result(outcome_for("degenerate_cell", f"det(cell) = {det:.6g}"), None)
    # 7. SCF
    scf_outcome, scf = scf_gate(step, run, policy)
    if scf_outcome is not None:
        return result(scf_outcome, None)
    # 8. magnetic policy
    reason = magnetic.get("policy_reason")
    if reason:
        mrun = run.magnetic.run if run.magnetic is not None else {}
        return result(outcome_for(reason, f"magnetic class {magnetic.get('magnetic_class')}: "
                                  f"{mrun.get('class_reason', '')}; segment {magnetic.get('segment')}"), None)
    # 9. force outlier (opt-in)
    if policy.force_outlier_threshold is not None and max_force is not None and \
            max_force > policy.force_outlier_threshold:
        return result(outcome_for("force_outlier", f"max |F| {max_force:.4g} > {policy.force_outlier_threshold:g} eV/A"),
                      None)
    # 10. the run's own outcome
    if not run.accepted:
        return result(outcome_for(run.outcome.reason, f"run-level: {run.outcome.detail}"), None)
    if forces is None:
        flags.append("energy_only")
    return result(Outcome(ACCEPTED), label_set)


def classify_steps(steps: Iterable[Any], run: RunAssessment, policy: Policy | None = None) -> list[StepAssessment]:
    return [classify_step(step, run, policy) for step in steps]


def outcome_counts(assessments: Iterable[Any]) -> dict[str, dict[str, int]]:
    """``{status: {reason or 'accepted': n}}`` for runs or steps (anything with ``.outcome``)."""
    counts: dict[str, dict[str, int]] = {}
    for item in assessments:
        outcome = item.outcome
        bucket = counts.setdefault(outcome.status, {})
        key = outcome.reason or "accepted"
        bucket[key] = bucket.get(key, 0) + 1
    return {status: dict(sorted(bucket.items())) for status, bucket in sorted(counts.items())}
