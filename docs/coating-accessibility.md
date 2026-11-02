# Coating continuity, probe accessibility, and route comparison

Two additive post-processing commands. Neither changes the production MD
workflow, force fields, or existing `analyze-coverage` / `analyze-interface`
outputs.

| Command | Unit of work | Output |
|---|---|---|
| `nio-md-prep analyze-accessibility RUN_DIR --trajectory STAGE.lammpstrj` | one trajectory window | `accessibility-analysis-STAGE/` |
| `nio-md-prep compare-routes --workbooks DIR [--accessibility-manifest CSV]` | all runs | route tables + figures |

## What is measured

**Surface reference.** All heights are measured from one reference: the
periodic upper envelope of substrate atom centers (each grid cell takes the
highest substrate z within `--reference-radius`, default 2.5 Å; this exceeds
the largest point-to-atom lateral distance on rock-salt NiO facets, Ni–O =
2.09 Å). The command fails if any cell is left without a substrate atom.
`analyze-coverage` instead uses a per-cell substrate maximum with a global
fallback. On the corrugated slab at its 0.2 Å grid, 99.1 % of cells contain
no substrate atom and take the slab maximum (39.6 Å), while the local surface
envelope spans 27.1 to 39.6 Å. The fallback reference therefore sits a median
2.1 Å (90th percentile 10.4 Å) above the local surface. Its "near-surface
≤ 5 Å" and anchor-conditioned gates are effectively a flat plane 5 Å above
the slab top, which admits groove material up to about 17 Å above the local
floor. Use the new reference for any height-gated quantity.

**Collective organization (substrate accessibility).**

- *Blocking model.* `canopy`: every ligand atom blocks a probe approaching
  along the surface normal (line of sight; equals projected coverage at probe
  radius 0). `hNA`: only ligand atoms within N Å of the reference block, so
  overhangs above N are ignored. This models lateral access through an
  interfacial slab. Defaults: N = 5, 8, 12 Å.
- *Finite probe.* The accessible area is the periodic morphological opening
  of the exposed mask by a disk of radius r. It is the set of surface points
  a disk of radius r can touch without overlapping a projected ligand van der
  Waals disk. Defaults: r = 0 (projected), 1.4 Å (water-sized), 2.0, and 3.0 Å.
  Geometric only: no solvent, no substrate shape in the probe test.
- *Gaps.* Periodic 4-connected components of the accessible mask, giving the
  area distribution, number density per nm², largest gap, and the fraction of
  accessible area in gaps ≥ 1 nm².
- *Persistence.* From the full stack of accessible masks:
  - per-pixel exposure probability (the coverage map);
  - the persistently exposed fraction (exposed in ≥ 90 % of frames);
  - exposure memory c(τ) = (P(exposed at t+τ | exposed at t) − φ)/(1 − φ),
    with the first lag where c < 1/e. It is reported as *censored* when the
    window never decorrelates;
  - first/last-frame and early/late-half Jaccard overlap against the overlap
    expected for independent masks of the same areas, which tests memory of
    the starting configuration;
  - pixel exposure run lengths.

**Molecular attachment and layer structure.**

- P, anchor-O, and top-heavy-atom heights per molecule. *Vertical extension*
  is top heavy-atom height minus lowest anchor-O height. This is adapted from
  ⟨z_C,max⟩ − ⟨z_O⟩ in Wang et al. (arXiv:2610.32905), but applied to
  aromatic and carbazole cores. Their ≈ 8 Å water-protection onset is specific
  to alkyl chains on Li₆PS₅Cl and is *not* used as a threshold here.
- Tilt: angle between the P→core-centroid vector and the surface normal
  (0° upright, 90° flat, > 90° pointing toward the surface). Core atoms are
  the same graph-distant set used by `analyze-interface`.
- Film thickness map: local top of ligand van der Waals spheres above the
  reference on a 1 Å grid, with the fraction of the cell at ≥ 5, 10, and 15 Å.
- Head-group contacts: each phosphonate anchor group (one per P; DCZ-4P has
  two) counts terminal O atoms within the cutoff of a canonical exposed-Ni
  site (the same site set as `analyze-interface`). Denticity is 0–3.
  - Cutoffs 3.0/3.25/3.5 Å, matching the existing interface workbook. A
    cosine-switched smooth contact (r_on 3.0, r_off 3.5 Å) and the O–Ni
    nearest-distance histogram are written so the choice can be checked
    against the first minimum.
  - Persistence outputs: ever-bound, persistently bound (≥ 90 % of frames),
    and intermittent fractions; mode populations; mode changes per bound
    group-frame; P(same mode | bound) and P(still bound | bound) vs lag; and
    early/late halves.

**Model scope.** Under `classical-ff` the anchoring persistence is the
persistence of nonbonded P–O···Ni contacts in the analyzed window, not a
chemical bond lifetime. No output measures moisture resistance, electronic
passivation, or degradation.

## Uncertainty

Every per-frame series is reported with its statistical inefficiency g (the
autocorrelation sum, truncated at the first non-positive value), an effective
sample count N/g, a correlation-corrected SEM, and early/late half means.
`window_too_short` is set when N/g < 5, and then no within-run uncertainty
should be quoted.

`compare-routes` never pools frames. Each run contributes one window mean:

- Replicate statistics: mean, SD, and 95 % t-interval across the three
  independently seeded cohorts (`prepared`, `prepared-rerun`,
  `prepared-rerun2`).
- Route effect: within-cohort Sequential − CoSAM differences (n = 3, df = 2).
  A modifier-averaged row is computed per cohort, because the three
  sequential runs of one cohort share that cohort's Me-4PACz stage-1 film.
- Route − Me-4PACz-reference contrasts, to tell "sequential gains" apart from
  "CoSAM loses".
- Gap-seeded I/II (one run each, built from the cohort-1 stage-1 film with an
  extra vertical deposition step) are placed against the replicate ranges and
  never pooled.
- 300/400 K and compressed/relaxed states of one run share their deposited
  structure. They test robustness to perturbation, not independence.

Verdicts: `supported` when the 95 % CI excludes 0; `consistent sign, CI
includes 0` when all cohorts agree in sign; otherwise `unresolved`.

## Running on LONI

```bash
find prepared prepared-rerun prepared-rerun2 prepared-lego prepared-lego2 \
  -mindepth 1 -maxdepth 1 -type d -exec test -s {}/relax-300K.lammpstrj \; -print \
  | sort > accessibility_runs.txt
sbatch --array=1-$(wc -l < accessibility_runs.txt) \
  --export=ALL,RUN_LIST=accessibility_runs.txt scripts/analyze_accessibility_array.sbatch
# afterwards, on a login node
python scripts/collect_accessibility_manifest.py --output accessibility_manifest.csv
nio-md-prep compare-routes --workbooks . \
  --accessibility-manifest accessibility_manifest.csv --output route-comparison
```

By default every frame of each 500 ps hold is analyzed (dump interval 1 ps).
Measured cost is about 0.35 s per frame for the grid operations on the
production cell; the estimate is about 1 s per frame with parsing and contact
search, or about 10 min per trajectory.
