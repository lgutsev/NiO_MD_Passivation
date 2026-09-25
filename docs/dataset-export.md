# Audited VASP → MLIP dataset export

`nio-md-prep dataset scan|export|split|audit` turns completed VASP calculations into a canonical, MACE/ASE-readable
extended-XYZ dataset with an explicit account of every discovered calculation and ionic frame. Its purpose is to bridge

```
VASP calculations -> audited canonical dataset -> lineage-safe train/valid/test -> MACE training
```

without silently accepting questionable labels. Discovery, parsing, classification and splitting need only the standard
library and numpy; ASE is imported lazily to write and re-read extended XYZ. Nothing here imports torch, MACE, OpenMM or
LAMMPS. Install the optional dependencies with `pip install -e '.[dataset]'`.

## Quick start (CPU only)

```bash
nio-md-prep dataset scan   --root agglo=/path/to/Agglo_DFT --output scan/          # accounting only
nio-md-prep dataset export --root agglo=/path/to/Agglo_DFT --output export/        # dataset + accounting + manifest
nio-md-prep dataset split  export/ --output split/ --seed 11 --fractions 0.8,0.1,0.1
nio-md-prep dataset audit  export/ --split split/ --verify-sources agglo=/path/to/Agglo_DFT
python scripts/dataset_spotcheck.py export/ --n 50        # compare exported frames with ASE's own VASP reader
```

Roots are `ALIAS=PATH` (or a bare path, alias = directory name). The alias is part of every run id, so keep it stable:
manifests remain valid when the data are moved and `--verify-sources ALIAS=NEWPATH` re-hashes them in place.
Output directories must be new or empty; `--force` writes into a non-empty directory file by file (atomic replacement,
nothing unrelated is deleted). Outputs are byte-identical for the same inputs, independent of root order and of the
order in which the filesystem lists directories (only `created_at` in the manifest changes, and it is excluded from
`content_sha256`).

## What is read, and from where

- **Label source:** `vasprun.xml` (plain or `.gz/.bz2/.xz`; plain wins when both exist and the other is reported).
  Geometry, forces, stress and energies of a frame always come from the same `<calculation>` block. `finalpos`, CONTCAR
  and XDATCAR are never paired with labels. Directories with an OUTCAR but no vasprun are reported as
  `missing_labels/no_vasprun`: OUTCAR-only label extraction is not implemented.
- **Evidence:** OUTCAR (per-step SCF markers, total magnetization, site-moment tables, echoed settings such as IVDW and
  the +U arrays, TOTEN cross-check), OSZICAR (`mag=`, temperature, SCF iteration counts), POTCAR (TITEL/VRHFIN/LEXCH/
  SHA256 fingerprints only — the file content is never stored), INCAR and POSCAR/CONTCAR (selective-dynamics flags,
  explicit MAGMOM). Evidence is mapped onto steps only when step counts agree; disagreement is reported, never guessed.
- **Streaming:** vasprun is parsed incrementally, so memory stays flat for long AIMD trajectories.
- **xTB never:** only VASP outputs define runs; agglomeration `xtb/`, `packmol/`, `templates/`, `structures/` and
  `vasp_reference/` trees are never ingested.

## Label rules (strict defaults; every non-default choice is recorded in the manifest)

**Energy.** The label is the force-consistent free energy F (VASP forces are derivatives of F):
`F = calculation e_fr_energy − PSTRESS·V`. E(σ→0) and E(no entropy) are kept for audit. VASP ≤ 6.0.8 mislabels the
calculation-level `e_wo_entrp`/`e_0_energy`; for those versions, and whenever the version is unknown or unverified, they
are reconstructed from the final SCF step. For VASP ≥ 6.1.1 the calculation-level values are used only if they agree with
the reconstruction within 1e-6 eV (otherwise the reconstruction is used and the frame flagged). Each frame records
`energy_source`, `energy_rule` (`calc_level_direct`, `reconstructed_last_scstep:<reason>`, …), the version gate, the
parser name and version. `--energy-quantity` changes the label only when asked; it is never switched silently.

**Forces are mandatory.** A frame without forces is `missing_labels/missing_forces` and never enters a force-training
dataset. Energy-only frames are exported only with `--allow-energy-only`, tagged `label_set=energy_only`, with no forces
key and a recorded `forces_reason`.

**Stress may be absent.** It is exported only with `--include-stress` and never as zeros: frames without a trustworthy
tensor carry no stress key, `stress_available=False` and a `stress_reason` (`not_requested`, `not_computed` (ISIF=0),
`isif1_trace_only`, `pstress_nonzero`, `vacuum_slab`, …). Stress is converted to the ASE convention (eV/Å³, positive =
tensile). Stress of slabs/clusters with vacuum needs `--allow-vacuum-stress` (it scales with the vacuum-inclusive volume).

**SCF convergence is gated conservatively.** A frame is converged only when vasprun shows fewer SCF steps than NELM and a
final energy change below EDIFF, and, when an OUTCAR from VASP ≥ 6 is present, that step's loop-exit marker says EDIFF was
reached. Explicit failure (NELM ceiling without convergence, `EDIFF was not reached`, hard stop) → `rejected`; missing or
ambiguous evidence (e.g. a VASP 5 marker at the NELM ceiling, unmappable OUTCAR) → `quarantined/scf_convergence_unknown`.
Single-marker text is never the only evidence. EDIFF = 0 runs need `--allow-ediff-zero`.

**VASP MLFF steps are never DFT labels.** Force-field-only steps appear outside `<calculation>` blocks; they are counted for
step indexing and lineage and reported as `non_dft_mlff_step/mlff_step`. Runs whose ML mode performs no DFT
(`run`, `select`, `refit*`, `delta`) are `non_dft_mlff_step/mlff_steps`.

**Selective dynamics.** Flags are read from vasprun (`initialpos`, else `finalpos`), else POSCAR/CONTCAR; conflicting
sources are quarantined. They are exported as the explicit per-atom array `vasp_selective_dynamics` (N×3 bool, True = may
move) with `selective_dynamics_basis="direct"` — VASP applies them to fractional coordinates. No ASE constraint is
attached to exported frames, and constrained atoms keep their raw forces (never zeroed).

**Other rejections.** Truncated steps, NaN/inf or VASP `****` overflow, degenerate or left-handed cells, DFPT/linear-response
runs, non-VASP generators. Interrupted runs are excluded unless `--interrupted-runs recover`, which judges complete steps
individually (the truncated final step is always rejected).

## Magnetic state (an audit dimension, not a blanket rule)

Per frame the export keeps the total magnetization (and its source), site moments when VASP printed them (LORBIT ≥ 10),
the initial MAGMOM, the change relative to the first and previous frame, the run's magnetic class and the policy result.
Classes: `non_spin_polarized` (ISPIN=1), `controlled_consistent` (explicit MAGMOM/NUPDOWN, no state change),
`controlled_transition` (explicit initialization, total moment jumped by more than 0.5 μB or the site sign pattern changed
beyond a global flip), `uncontrolled` (ISPIN=2 without explicit MAGMOM — VASP's ferromagnetic default start), `unknown`
(ISPIN=2 without usable magnetization evidence).

Default policy for runs containing magnetic species (Co, Cr, Cu, Fe, Mn, Ni, V, plus any element with |MAGMOM| ≥ 0.5):
`uncontrolled` and `unknown` runs are quarantined, and frames after a controlled transition are quarantined, so distinct
electronic/magnetic branches are never mixed into one ordinary dataset. Overrides are explicit and recorded:
`--accept-uncontrolled-magnetism`, `--accept-unknown-magnetism`, `--accept-magnetic-transitions`, each requiring
`--magnetic-override-reason`. Exact duplicate geometries with different magnetization are `contradictory_labels`.

The legacy local NiO/Me-4PACz MD packages (OutPackLite, ISPIN=2 with no MAGMOM) are an example of `uncontrolled` data: they
also lack OUTCAR/vasprun and therefore force labels, so they are reported, never exported.

## Reference-settings pools

Frames computed with incompatible DFT settings are kept in separate pools; `export` writes exactly one pool (the only one,
or `--pool ID`). Method fields: functional (GGA/METAGGA/hybrid tags), IVDW, ENCUT, PREC, LASPH, ISPIN, spin–orbit/
non-collinear, fixed NUPDOWN, ISMEAR/SIGMA, net charge, per-element POTCAR identity and per-element effective DFT+U. Pools
are element-aware (a Ni-free cluster can share a pool with a NiO slab when all shared elements agree). A field whose value
cannot be determined is `<unknown>` and never equals a known value. k-points, EDIFF, NELM, ALGO, POTIM and similar sampling
settings are recorded but do not split pools. Reviewed equivalences go in a `--settings-overrides` TOML of
`[[equivalence]]` tables (`field`, `values`, `reason`, `reviewed_by`, all required) and are recorded in the manifest.

## Lineage and leakage-safe splits

Every frame carries `lineage_group`; `split` assigns whole groups, never individual frames. Lineage sources, in priority:
inventory TOML (`--inventory`), agglomeration manifests (one Packmol replica = one group, including its 300 K hold, heating
ramp, 400 K hold and single points), InterfaceForge campaign manifests (`opt_manifest.json`, `step1_manifest.json`,
`step2_manifest.json`, joined on `relative_path`, so OPT, Step1 and every Step2 temperature of a structure share a group),
then a declared `--lineage-policy run|parent|depth:N`. Runs are additionally linked when one starts from a frame of another
(restart/continuation), when they contain exact duplicate structures, and when nested. Runs without any lineage are
quarantined `lineage_unresolved` and `split` refuses them. `.interfaceforge/`, `precondition/` and `X*` (OutPackLite) path
parts are excluded by default and reported; `--include-excluded-part PATTERN` includes them deliberately.

`split` is seeded and deterministic, reports requested and achieved proportions in groups and frames, supports
`--stratify-by KEY` (refusing impossible stratifications unless `--allow-underpopulated-strata`), and with
`--previous split_manifest.json` keeps earlier assignments and freezes the test set (`--grow-test` to add to it). Groups that
merged across different previous splits raise a leakage error.

## Duplicates

Structures are keyed with wrap-aware, permutation-invariant hashes (1e-6 Å resolution). Exact duplicates in one pool with
consistent labels (|ΔE| ≤ 1e-3 eV, max |ΔF| ≤ 1e-2 eV/Å, |Δmag| ≤ 0.1 μB, aligned site-moment signs) keep the lowest frame
id; the others are `excluded/exact_duplicate`. Inconsistent copies are all `quarantined/contradictory_labels`. Duplicates
always join lineage groups. `--near-duplicate-tolerance Å` links (never removes) near-duplicates across groups.

## Outputs

`export/` contains `dataset.extxyz`, `runs.jsonl` and `frames.jsonl` (every discovered run and frame with its state),
`exclusions.csv`, `audit_report.json`, `audit.md` and `dataset_manifest.json` (inputs, aliases, options and full policy,
tool/parser version, git commit, label keys, units, pools, counts, per-file sha256, `content_sha256`). `split/` contains
`train.extxyz`, `valid.extxyz`, `test.extxyz`, `split_manifest.json` and `split_summary.md`. The audit reports counts by
campaign, family, temperature, composition, lineage group, VASP version, state/reason, force and stress availability and
magnetic class.

Each exported frame carries `REF_energy` and `REF_forces` (keys configurable; names that ASE would turn into calculator
results are refused), optional `REF_stress` (3×3), the energy variants `vasp_free_energy`, `vasp_energy_no_entropy`,
`vasp_energy_sigma0`, and the metadata `source`, `source_file_type`, `ionic_step`, `structure_key`, `lineage_group`,
`energy_source`, `energy_rule`, `stress_available`, `scf_status`, `label_source`, `vasp_version`, `magnetic_class`,
`magnetic_policy`, `campaign`, `family`, `parser`, `parser_version`, `repo_commit` (plus `split` in split files) and the
arrays `vasp_selective_dynamics`, `vasp_magmom_initial`, `vasp_magmom_final` where applicable. MACE reads the labels with
`--energy_key REF_energy --forces_key REF_forces [--stress_key REF_stress]`; a frame without the stress key is masked by MACE
(weight 0), which is why unavailable stress is omitted rather than zeroed.

## States and reasons

Every discovered candidate has exactly one state; each reason code belongs to one state.

| state | reason | meaning |
|---|---|---|
| `rejected` | `degenerate_cell` | cell matrix is singular or left-handed |
| `rejected` | `non_finite` | NaN/inf or VASP overflow ('****') in positions, cell, energy, forces or stress |
| `rejected` | `not_vasp` | generator program is not VASP |
| `rejected` | `scf_not_converged` | electronic loop hit NELM without meeting EDIFF (explicit evidence) |
| `rejected` | `truncated_frame` | ionic step block is incomplete (file ends inside it) |
| `rejected` | `unsupported_calculation` | IBRION/other settings describe a calculation whose frames are not energy/force labels (e.g. DFPT) |
| `quarantined` | `contradictory_labels` | exact structural duplicate with inconsistent energy/forces/magnetization under the same reference settings |
| `quarantined` | `evidence_mismatch` | OUTCAR/OSZICAR evidence disagrees with vasprun (different run, copied file, or corrupted output) |
| `quarantined` | `force_outlier` | max \|F\| exceeds the declared quarantine threshold |
| `quarantined` | `lineage_unresolved` | no lineage source and no declared grouping policy (split refuses such frames) |
| `quarantined` | `magnetic_order_changed` | final per-atom moment signs differ from the initial MAGMOM pattern (beyond a global flip) |
| `quarantined` | `magnetic_state_change` | total magnetization jumped between ionic steps beyond the declared threshold |
| `quarantined` | `magnetic_uncontrolled` | spin-polarized run without explicit MAGMOM initialization (policy: quarantined unless overridden) |
| `quarantined` | `magnetic_unknown` | spin-polarized run containing magnetic species but no magnetization evidence (policy: quarantined unless overridden) |
| `quarantined` | `scf_convergence_unknown` | no usable evidence that the electronic loop converged for this step |
| `quarantined` | `species_mismatch` | species/order disagree between vasprun, POTCAR, POSCAR or OUTCAR |
| `parse_error` | `bad_shape` | array shape does not match the atom count |
| `parse_error` | `vasprun_unreadable` | vasprun.xml header could not be parsed |
| `missing_labels` | `missing_energy` | no free energy for this ionic step |
| `missing_labels` | `missing_forces` | no forces for this ionic step |
| `missing_labels` | `no_vasprun` | no vasprun.xml(.gz/.bz2/.xz) in the run directory; OUTCAR-only label extraction is not implemented |
| `non_dft_mlff_step` | `mlff_step` | this ionic step was predicted by the VASP MLFF, not computed by DFT |
| `non_dft_mlff_step` | `mlff_steps` | VASP machine-learned force field active; ionic steps are not all DFT |
| `excluded` | `archived_path` | path matches an archived/stale-artifact pattern |
| `excluded` | `copied_reference_output` | output file is byte-identical to a file copied from a reference/template directory |
| `excluded` | `exact_duplicate` | exact structural duplicate of an earlier frame with consistent labels (kept once) |
| `excluded` | `inventory_excluded` | inventory manifest explicitly excludes this run |
| `excluded` | `not_in_inventory` | inventory manifest selects listed runs only and this run is not listed |
| `excluded` | `preconditioning_run` | static wavefunction-preconditioning calculation (duplicates the following run's first frame) |
| `excluded` | `reference_pool_not_selected` | run belongs to a different reference-settings pool than the one exported |
| `excluded` | `run_incomplete` | run did not finish (vasprun not closed / OUTCAR timing block absent) and the interrupted-run policy is 'exclude' |
| `excluded` | `run_not_accepted` | the parent run was not accepted |
| `excluded` | `subsampled` | dropped by the declared stride/max-frames subsampling |

## Validation on LONI (pull-only)

Real NiO force export has **not** been validated: the local copies of the NiO calculations are OutPackLite packages without
OUTCAR/vasprun.xml. Parser tests are synthetic (clearly labelled) plus ASE's shipped real VASP files read in place. On LONI,
from an up-to-date pull-only checkout with `.venv` holding `pip install -e '.[dataset]'`:

```bash
sbatch -A <allocation> --export=ALL,DATASET_ROOT="agglo=/path/to/Agglo_DFT iface=/path/to/NiO_campaign" \
    scripts/dataset_export_validation.sbatch
```

The job runs scan → export → split → audit (with source re-hashing) → `scripts/dataset_spotcheck.py`, which re-reads a
sample of exported frames from their original vasprun.xml with ASE's independent reader and reports the largest free-energy
and raw-force differences. Before trusting a dataset, inspect `scan/audit.md` (pools, states, magnetic classes) and compare a
few frames by hand against the original VASP energies and forces.

## Limitations

- OUTCAR-only runs are reported but not exported (no vasprun label source).
- InterfaceForge `step1-repair` archives stay excluded; recovering their validated prefix is not implemented.
- Real NiO force-export validation is outstanding (see above); synthetic tests do not establish it.
- Pool assignment is greedy in sorted run-id order: deterministic, but a run compatible with several pools joins the first.
- Stress training on vacuum slabs is a scientific choice; the defaults do not make it.
