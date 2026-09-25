"""Shared records for dataset discovery, parsing, classification and export.

Two layers live here:

* **Parser output** (:class:`VasprunHeader`, :class:`StepStructure`,
  :class:`IonicStep`, :class:`VasprunTrailer`) is a faithful, unit-labelled
  transcription of one ``vasprun.xml``. Nothing in it is inferred or repaired.
* **Accounting records** (:class:`RunRecord`, :class:`FrameRecord`) carry the
  classification decided by the policy modules, with an :class:`Outcome` for
  every discovered run and every discovered ionic frame.

Units in parser output are VASP's own except where a field name says
otherwise: positions are stored both fractional and Cartesian (Angstrom);
energies are eV; forces are eV/Angstrom; ``stress_kbar_vasp`` is the 3x3
tensor exactly as VASP prints it (kB, VASP sign: positive = compressive).
Conversion to the canonical ASE convention happens in exactly one place,
:func:`vasp_stress_to_ase`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# --------------------------------------------------------------------------
# Outcomes and reason codes
# --------------------------------------------------------------------------

ACCEPTED = "accepted"
REJECTED = "rejected"  # label validity failed (SCF not converged, non-finite, truncated, unsupported calculation)
QUARANTINED = "quarantined"  # plausible data that needs human review (contradictory labels, magnetic branch, evidence mismatch)
PARSE_ERROR = "parse_error"  # the source could not be parsed; nothing is silently dropped
MISSING_LABELS = "missing_labels"  # geometry present but energy and/or forces absent (no vasprun, no forces block, ...)
NON_DFT_MLFF_STEP = "non_dft_mlff_step"  # VASP-MLFF prediction; counted for indexing/lineage, never a DFT label
EXCLUDED = "excluded"  # valid data deliberately not used (duplicate, other reference pool, subsampled, archived path)

OUTCOMES = (ACCEPTED, REJECTED, QUARANTINED, PARSE_ERROR, MISSING_LABELS, NON_DFT_MLFF_STEP, EXCLUDED)

# Backwards-compatible alias used by early module drafts.
UNREADABLE = PARSE_ERROR

#: Magnetic-state policy classes (per run; frames inherit, transitions are per frame).
MAG_NON_SPIN_POLARIZED = "non_spin_polarized"  # ISPIN=1: no magnetic degree of freedom by construction
MAG_CONTROLLED_CONSISTENT = "controlled_consistent"  # explicit MAGMOM (or NUPDOWN) and no state change observed
MAG_CONTROLLED_TRANSITION = "controlled_transition"  # explicit initialization, but the state changed during the run
MAG_UNCONTROLLED = "uncontrolled"  # ISPIN=2 without explicit MAGMOM: VASP's default ferromagnetic start
MAG_UNKNOWN = "unknown"  # ISPIN=2 but no usable magnetization evidence
MAGNETIC_CLASSES = (
    MAG_NON_SPIN_POLARIZED,
    MAG_CONTROLLED_CONSISTENT,
    MAG_CONTROLLED_TRANSITION,
    MAG_UNCONTROLLED,
    MAG_UNKNOWN,
)

#: Machine-readable reason codes. Every non-accepted outcome uses one of these
#: (policy modules may add a free-text ``detail``). Keep this table the single
#: source of truth; docs/dataset-export.md lists the same codes.
REASONS: dict[str, str] = {
    # run-level
    "no_vasprun": "no vasprun.xml(.gz/.bz2/.xz) in the run directory; OUTCAR-only label extraction is not implemented",
    "vasprun_unreadable": "vasprun.xml header could not be parsed",
    "not_vasp": "generator program is not VASP",
    "unsupported_calculation": "IBRION/other settings describe a calculation whose frames are not energy/force labels (e.g. DFPT)",
    "mlff_steps": "VASP machine-learned force field active; ionic steps are not all DFT",
    "run_incomplete": "run did not finish (vasprun not closed / OUTCAR timing block absent) and the interrupted-run policy is 'exclude'",
    "species_mismatch": "species/order disagree between vasprun, POTCAR, POSCAR or OUTCAR",
    "evidence_mismatch": "OUTCAR/OSZICAR evidence disagrees with vasprun (different run, copied file, or corrupted output)",
    "copied_reference_output": "output file is byte-identical to a file copied from a reference/template directory",
    "archived_path": "path matches an archived/stale-artifact pattern",
    "not_in_inventory": "inventory manifest selects listed runs only and this run is not listed",
    "inventory_excluded": "inventory manifest explicitly excludes this run",
    "reference_pool_not_selected": "run belongs to a different reference-settings pool than the one exported",
    "preconditioning_run": "static wavefunction-preconditioning calculation (duplicates the following run's first frame)",
    # frame-level
    "truncated_frame": "ionic step block is incomplete (file ends inside it)",
    "missing_forces": "no forces for this ionic step",
    "missing_energy": "no free energy for this ionic step",
    "bad_shape": "array shape does not match the atom count",
    "non_finite": "NaN/inf or VASP overflow ('****') in positions, cell, energy, forces or stress",
    "degenerate_cell": "cell matrix is singular or left-handed",
    "scf_not_converged": "electronic loop hit NELM without meeting EDIFF (explicit evidence)",
    "scf_convergence_unknown": "no usable evidence that the electronic loop converged for this step",
    "mlff_step": "this ionic step was predicted by the VASP MLFF, not computed by DFT",
    "subsampled": "dropped by the declared stride/max-frames subsampling",
    "exact_duplicate": "exact structural duplicate of an earlier frame with consistent labels (kept once)",
    "contradictory_labels": "exact structural duplicate with inconsistent energy/forces/magnetization under the same reference settings",
    "magnetic_state_change": "total magnetization jumped between ionic steps beyond the declared threshold",
    "magnetic_order_changed": "final per-atom moment signs differ from the initial MAGMOM pattern (beyond a global flip)",
    "magnetic_uncontrolled": "spin-polarized run without explicit MAGMOM initialization (policy: quarantined unless overridden)",
    "magnetic_unknown": "spin-polarized run containing magnetic species but no magnetization evidence (policy: quarantined unless overridden)",
    "force_outlier": "max |F| exceeds the declared quarantine threshold",
    "run_not_accepted": "the parent run was not accepted",
    "lineage_unresolved": "no lineage source and no declared grouping policy (split refuses such frames)",
}

#: Each reason code maps to exactly one outcome status, so the same failure is
#: never reported under two different states.
REASON_STATUS: dict[str, str] = {
    "no_vasprun": MISSING_LABELS,
    "missing_forces": MISSING_LABELS,
    "missing_energy": MISSING_LABELS,
    "vasprun_unreadable": PARSE_ERROR,
    "not_vasp": REJECTED,
    "unsupported_calculation": REJECTED,
    "mlff_steps": NON_DFT_MLFF_STEP,
    "mlff_step": NON_DFT_MLFF_STEP,
    "run_incomplete": EXCLUDED,
    "species_mismatch": QUARANTINED,
    "evidence_mismatch": QUARANTINED,
    "copied_reference_output": EXCLUDED,
    "archived_path": EXCLUDED,
    "not_in_inventory": EXCLUDED,
    "inventory_excluded": EXCLUDED,
    "reference_pool_not_selected": EXCLUDED,
    "preconditioning_run": EXCLUDED,
    "truncated_frame": REJECTED,
    "bad_shape": PARSE_ERROR,
    "non_finite": REJECTED,
    "degenerate_cell": REJECTED,
    "scf_not_converged": REJECTED,
    "scf_convergence_unknown": QUARANTINED,
    "subsampled": EXCLUDED,
    "exact_duplicate": EXCLUDED,
    "contradictory_labels": QUARANTINED,
    "magnetic_state_change": QUARANTINED,
    "magnetic_order_changed": QUARANTINED,
    "magnetic_uncontrolled": QUARANTINED,
    "magnetic_unknown": QUARANTINED,
    "force_outlier": QUARANTINED,
    "run_not_accepted": EXCLUDED,
    "lineage_unresolved": QUARANTINED,
}
assert set(REASON_STATUS) == set(REASONS), "every reason needs exactly one status"


def outcome_for(reason: str, detail: str = "") -> "Outcome":
    """The one way policy code builds a non-accepted outcome: status follows the reason."""
    return Outcome(REASON_STATUS[reason], reason, detail)


@dataclass(frozen=True)
class Outcome:
    status: str
    reason: str | None = None
    detail: str = ""

    def __post_init__(self):
        if self.status not in OUTCOMES:
            raise ValueError(f"unknown outcome status {self.status!r}")
        if self.status != ACCEPTED and self.reason not in REASONS:
            raise ValueError(f"non-accepted outcome needs a known reason code; got {self.reason!r}")
        if self.status != ACCEPTED and REASON_STATUS[self.reason] != self.status:
            raise ValueError(
                f"reason {self.reason!r} is reported as {REASON_STATUS[self.reason]!r}, not {self.status!r}"
            )

    def as_dict(self) -> dict[str, Any]:
        return {"status": self.status, "reason": self.reason, "detail": self.detail}


# --------------------------------------------------------------------------
# Parser output (vasprun.xml transcription)
# --------------------------------------------------------------------------

@dataclass
class AtomType:
    """One ``<atominfo><array name="atomtypes">`` row."""

    element: str
    count: int
    mass: float | None
    valence: float | None  # ZVAL, used for the net-charge check
    pseudopotential: str  # e.g. "PAW_PBE Ni_pv 06Sep2000" as printed


@dataclass
class VasprunHeader:
    """Everything in vasprun.xml before the first ``<calculation>``."""

    generator: dict[str, str]  # program, version, subversion, platform, date, time
    incar: dict[str, Any]  # <incar> as written (user-given tags only)
    parameters: dict[str, Any]  # <parameters> flattened name -> typed value (full defaults)
    kpoints: dict[str, Any]  # generation scheme/divisions/shift + number of k-points
    species: list[str]  # per atom, in file order
    atom_types: list[AtomType]
    initial_structure: "StepStructure | None" = None
    #: Problems found while reading the header. Each entry starts with a REASONS
    #: code or a short tag followed by ": " (e.g. "non_finite: parameters.X ...").
    problems: list[str] = field(default_factory=list)
    version_tuple: tuple[int, int, int] | None = None  # parsed from generator["version"]
    #: <parameters> names seen in more than one separator; the first value is kept.
    parameter_duplicates: list[str] = field(default_factory=list)


@dataclass
class StepStructure:
    cell: Any  # (3, 3) float64, Angstrom, rows are lattice vectors
    fractional: Any  # (N, 3) float64
    positions: Any  # (N, 3) float64 Cartesian Angstrom = fractional @ cell
    selective: Any | None = None  # (N, 3) bool, True = may move (VASP 'T'), None if not selective dynamics
    volume: float | None = None  # as printed by VASP (A^3)


@dataclass
class IonicStep:
    """One ``<calculation>`` block: geometry and the labels computed for it."""

    index: int  # 0-based position of the block in the file
    complete: bool  # closing </calculation> was seen
    structure: StepStructure | None
    forces: Any | None  # (N, 3) eV/A
    stress_kbar_vasp: Any | None  # (3, 3) kB, VASP sign convention
    energies: dict[str, float]  # raw calc-level <energy> <i> values as printed (e_fr_energy, ...; MD: kinetic, total, ...)
    scf_energies: list[dict[str, float]]  # one dict per <scstep> (raw values)
    problems: list[str] = field(default_factory=list)  # parse problems local to this block (e.g. non-numeric value)
    #: "dft" for a <calculation> block; "mlff" for a VASP-MLFF force-field-only flat step
    #: (top-level structure/forces/energy sequence without <scstep>). Never exported as DFT.
    label_source: str = "dft"
    #: <time name="..."> values of the block (e.g. "totalsc": [cpu, wall] seconds); never labels.
    times: dict[str, list[float]] = field(default_factory=dict)


@dataclass
class VasprunTrailer:
    """What is known once iteration over the file has stopped."""

    closed: bool  # </modeling> seen
    final_structure_present: bool  # <structure name="finalpos"> seen
    truncated: bool  # XML ended early (ParseError at EOF)
    error: str | None  # parse error text, if any
    steps_seen: int
    final_structure: StepStructure | None = None  # <structure name="finalpos"> (NEVER paired with labels)
    n_complete_steps: int = 0  # steps yielded with complete=True
    partial_tail_tag: str | None = None  # depth-1 element open when the file ended early
    source_sha256: str | None = None  # sha256 of the file bytes as stored (compressed bytes for .gz etc.)
    source_bytes: int | None = None
    compression: str | None = None  # None | "gz" | "bz2" | "xz" (detected from magic bytes)
    problems: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------
# Evidence from auxiliary files (OUTCAR / OSZICAR / INCAR / POTCAR / POSCAR)
# --------------------------------------------------------------------------

@dataclass
class OutcarEvidence:
    sha256: str
    version: str | None
    nions: int | None
    ions_per_type: list[int] | None
    potcar_titles: list[str]
    scf_converged_markers: list[bool | None]  # per ionic step: True 'EDIFF is reached', False 'not reached', None absent
    free_energies: list[float]  # 'free  energy   TOTEN' per ionic step
    energies_no_entropy: list[float]
    energies_sigma0: list[float]
    magnetization_tables: dict[int, list[float]]  # ionic step index -> per-atom total moments (only where printed)
    total_magnetizations: list[float | None]  # per ionic step, last 'magnetization' of the electronic loop
    ionic_converged: bool  # 'reached required accuracy'
    completed: bool  # 'General timing and accounting informations for this job'
    problems: list[str] = field(default_factory=list)
    version_tuple: tuple[int, int, int] | None = None
    executed_tags: dict[str, str] = field(default_factory=dict)  # INCAR/parameter echo (IVDW, LDAUU, NELM, EDIFF, ...)
    #: Per DFT ionic step, memory-light cross-check values instead of full arrays:
    #: sum over atoms of |Fx|+|Fy|+|Fz| from TOTAL-FORCE (eV/A), and the last
    #: "total energy-change (2. order)" pair (dE, d_eps) of the electronic loop.
    force_abs_sums: list[float] = field(default_factory=list)
    last_energy_changes: list[tuple[float, float] | None] = field(default_factory=list)
    dispersion_energies: list[float | None] = field(default_factory=list)  # 'Edisp (eV)' per step when printed
    ml_steps: int = 0  # number of 'ML' prediction blocks seen (VASP-MLFF runs)
    #: Per DFT step: "reached" | "not_reached" | "hard_stop" | "other" | None (no 'aborting loop' line).
    scf_marker_kinds: list[str | None] = field(default_factory=list)
    scf_iteration_counts: list[int] = field(default_factory=list)  # 'Iteration k(n)' lines per DFT step
    stress_kbar: list[list[float] | None] = field(default_factory=list)  # 'in kB' XX YY ZZ XY YZ ZX per DFT step
    #: POTCAR headers echoed after the second half of the 'POTCAR:' lines (titel, vrhfin, lexch; absent with NWRITE=0).
    potcar_headers: list[dict[str, str]] = field(default_factory=list)
    #: 'magnetization (x)' table printed after the last step's summary (VASP repeats the final table).
    final_magnetization_table: list[float] | None = None
    forces: list[Any] | None = None  # full (N, 3) TOTAL-FORCE arrays per DFT step, only with keep_forces=True


@dataclass
class OszicarEvidence:
    sha256: str
    ionic_lines: list[dict[str, float]]  # per ionic step: F, E0, (T, E, EK, SP, SK), mag, dE
    scf_steps: list[int]  # electronic iterations per ionic step
    problems: list[str] = field(default_factory=list)
    line_kinds: list[str] = field(default_factory=list)  # per ionic line: "md" | "relax"


@dataclass
class PotcarEvidence:
    sha256: str  # whole-file digest (never the content)
    #: Per dataset (amendments item 11): symbol, element, titel, vrhfin, lexch, zval, pomass, enmax, enmin,
    #: sha256_header (SHA256 line or None), sha256_header_verified (bool or None when no header),
    #: sha256_verify_variant, sha256_dataset_bytes. Never any POTCAR content.
    datasets: list[dict[str, Any]]
    problems: list[str] = field(default_factory=list)


@dataclass
class IncarEvidence:
    """An INCAR file as VASP would read it: upper-case tag -> raw value string (last assignment wins)."""

    sha256: str
    tags: dict[str, str]
    duplicates: list[str] = field(default_factory=list)  # tags assigned more than once
    problems: list[str] = field(default_factory=list)


@dataclass
class PoscarEvidence:
    """A POSCAR/CONTCAR file (VASP 5 layout; VASP 4 files have species=None)."""

    sha256: str
    comment: str
    cell: Any  # (3, 3) Angstrom after applying the scale factor
    species_tokens: list[str] | None  # raw species line tokens (e.g. "Na_pv/6a2f546d")
    species: list[str] | None  # per type, element symbol only
    counts: list[int]
    coordinate_mode: str  # "direct" | "cartesian"
    fractional: Any  # (N, 3)
    positions: Any  # (N, 3) Cartesian Angstrom
    selective: Any | None = None  # (N, 3) bool, True = may move, DIRECT basis (VASP convention); None if absent
    problems: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------
# Accounting records
# --------------------------------------------------------------------------

@dataclass
class RunRecord:
    run_id: str  # "<alias>:<posix relpath of run dir>"
    root_alias: str
    relpath: str
    source_file: str | None  # file name of the label source (vasprun.xml[.gz])
    source_sha256: str | None
    source_bytes: int | None
    calc_type: str | None  # static | relaxation | md | other
    generator: dict[str, str] = field(default_factory=dict)
    settings: dict[str, Any] = field(default_factory=dict)  # see settings.extract_settings
    method_fingerprint: str | None = None
    pool_id: str | None = None
    species: list[str] = field(default_factory=list)
    n_atoms: int | None = None
    selective: bool = False
    magmom_initial: list[float] | None = None
    completion: dict[str, Any] = field(default_factory=dict)
    evidence_files: dict[str, Any] = field(default_factory=dict)  # name -> {sha256, bytes, used}
    metadata: dict[str, Any] = field(default_factory=dict)  # family, campaign, composition, temperature, config_type, ...
    lineage: dict[str, Any] = field(default_factory=dict)  # {lineage_id, source, evidence}
    group_id: str | None = None
    frames_total: int = 0
    frames_accepted: int = 0
    flags: list[str] = field(default_factory=list)  # non-fatal observations (e.g. magnetic_state_unresolved)
    outcome: Outcome = field(default_factory=lambda: Outcome(ACCEPTED))


@dataclass
class FrameRecord:
    frame_id: str  # "<run_id>#<index:05d>"
    run_id: str
    index: int
    n_atoms: int | None
    time_fs: float | None = None
    structure_key: str | None = None  # exact, order-dependent structure hash
    permutation_key: str | None = None  # exact, permutation-invariant structure hash
    label_sha256: str | None = None  # digest of the exported label arrays (energy, forces, [stress])
    energies: dict[str, float] = field(default_factory=dict)
    label_energy: float | None = None
    max_force: float | None = None
    scf: dict[str, Any] = field(default_factory=dict)  # {n_steps, nelm, last_dE, ediff, converged, evidence}
    stress_available: bool = False
    stress_reason: str | None = None
    total_magnetization: float | None = None
    magnetization_source: str | None = None  # outcar | oszicar
    site_magmoms: list[float] | None = None  # only on steps where VASP printed a per-ion table
    magnetic_state_id: str | None = None  # canonical sign pattern (global flip removed) + ispin + nupdown
    magnetic_segment: int | None = None
    label_source: str = "dft"
    md: dict[str, float] = field(default_factory=dict)  # kinetic, total, temperature_K when IBRION=0
    group_id: str | None = None
    duplicate_of: str | None = None
    flags: list[str] = field(default_factory=list)
    outcome: Outcome = field(default_factory=lambda: Outcome(ACCEPTED))


# --------------------------------------------------------------------------
# The one unit conversion
# --------------------------------------------------------------------------

#: 1 eV/A^3 = 160.21766208 GPa (CODATA e / 1e-30 m^3 * 1e-9); 1 kB = 0.1 GPa.
EV_PER_A3_IN_KBAR = 1602.1766208


def vasp_stress_to_ase(stress_kbar_vasp):
    """VASP stress (kB, positive = compressive) -> ASE stress (eV/A^3, positive = tensile).

    ASE's convention is sigma = (1/V) dE/d(strain); VASP prints the negative
    of that (a pressure-like quantity) in kB. The tensor is returned 3x3, in
    the same Cartesian basis as the cell rows.
    """
    import numpy as np

    stress = np.asarray(stress_kbar_vasp, dtype=np.float64)
    if stress.shape != (3, 3):
        raise ValueError(f"expected a 3x3 stress tensor, got shape {stress.shape}")
    return -stress / EV_PER_A3_IN_KBAR
