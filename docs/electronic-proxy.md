# Tier-0 electronic decision gate

`nio-md-prep analyze-electronic-proxy` decides whether a DFT work-function study of assembly routes is worth proposing. It launches no DFT. It asks one question: do two assembly routes (by default DCZ-4P CoSAM versus conventional sequential) differ by more than their between-run spread in quantities that would control an interface work-function change?

The outputs are an **internal decision gate** under the `classical-ff` evidence scope. They are not work-function predictions and do not belong in manuscript tables. The classical trajectories contain neutral, fully protonated `PO3H2` anchors with no Ni–O–P bonding or proton transfer, so the bond-dipole contribution to the work function is absent by construction.

## What is computed

### Facet classification

Each canonical exposed Ni site (see [`interfacial-analysis.md`](interfacial-analysis.md#exposed-ni-site)) is labelled by the direction of its missing O neighbours on the pristine slab. That direction is snapped to (001) or to one of the 45° {110} facets.

The supplied corrugated surface yields:

| Facet | Sites | Arrangement |
|---|---:|---|
| `110+y` | 210 | 7 rows × 30 Ni |
| `110-y` | 210 | 7 rows × 30 Ni |
| `001` | 180 | ridge top and groove bottom |

The {110} facets have the same termination as the 200-atom NiO(110) DFT reference slab.

### Patch windows

A window is a fixed set of facet-site identities: `--window-rows` facet rows (default 5) × `--window-columns` Ni per row (default 4). The default is 16.7 × 14.7 Å, the footprint of the 200-atom DFT slab. On the supplied surface, `--window-column-stride 2` gives 90 overlapping windows.

Because windows overlap, patch distributions are descriptive and are not independent samples.

### Molecule assignment

Each molecule's phosphonate-P centroid is projected along the facet normal onto the nearest site of that facet. A molecule is accepted when it lies between `--height-below` and `--height-above` along the normal and within one lattice spacing in-plane.

All facets, including (001), compete for each molecule, and the facet with the smallest normal height wins. Without this rule, projection through the V-groove would place material from the opposite facet into a patch.

Overlayer molecules are kept, because they contribute to a work-function change.

### Patch descriptors (`patch_candidates.csv`, one row per frame × window)

- molecules, anchored terminals, and anchored-terminal density per component;
- molecules with two or more anchored terminals, and unanchored molecules;
- lifted molecules: unanchored, with P more than `--lift-height` above the facet;
- fraction of window sites anchored, and fraction exposed (no ligand heavy atom within `--shield-cutoff`);
- tilt of the P→core vector from the facet normal;
- maximum molecular height.

### Dipole-layer proxy

```text
dphi = -(e / eps0) * sum_i (mu_i . n) / (eps_eff * A)
```

- μ points from negative to positive charge. A positive `dphi` raises the work function.
- For a patch, `n` is the facet normal and `A` is the window area, scaled for NPT.
- For the whole interface (`frame_proxies.csv`), `n = z` and `A = Lx·Ly`. A laterally periodic dipole layer sets the macroscopic potential step through the total z-dipole per projected area, so this quantity is well defined for the corrugated film.

Each molecular dipole μ_i comes from one of two sources:

- **Charges (always available):** `Σ q (r − centroid)` over the molecule, unwrapped through its bond graph, using the trajectory `q` column or the topology charges.
- **Reference (optional, `--reference-dipoles`):** a quantum-chemical dipole for each species (e.g. GFN2-xTB or DFT on the LigParGen geometry), rotated onto every conformer by a Kabsch fit. This rigid-rotation approximation ignores conformational changes of the dipole. The mean fit RMSD is reported so it can be judged.

The proxy omits the bond dipole, charge transfer and substrate polarization. Depolarization enters only through the scalar `--epsilon-eff`, so it rescales every value equally and cannot change the gate.

## Reference dipole file

```json
{
  "species": {
    "me-4pacz": {
      "dipole_debye": [0.0, 0.0, 0.0],
      "method": "GFN2-xTB gas phase, LigParGen geometry",
      "geometry": "optional/xyz/in/topology/atom/order.xyz"
    },
    "dcz-4p": {"dipole_e_angstrom": [0.0, 0.0, 0.0], "method": "..."}
  }
}
```

- The dipole must be expressed in the frame of the reference geometry. If no geometry is given, the default is `inputs/molecules/<slug>/ligpargen.lmp`.
- The geometry's atom order must match the topology order.
- Species are matched by manifest component name. Sequential and gap-seeded builds call the preformed layer `stage1_primary`, so pass `--component-alias stage1_primary=me-4pacz`. The alias is recorded in the summary.
- **The values above are placeholders.** Fill them only from a calculation you have run.

### Validating the rigid-rotation approximation

`--export-conformers N` writes N sampled MD conformers per component to `conformers/<component>.extxyz`. Each entry records its charge dipole and its rotated reference dipole. Run the same quantum-chemical method on a subset and compare against the rotated values before trusting the reference proxy.

## Running the gate on QBD

```bash
sbatch scripts/analyze_electronic_tier0.sbatch
sbatch --export=ALL,REFERENCE_DIPOLES=$PWD/reference_dipoles.json scripts/analyze_electronic_tier0.sbatch
```

The launcher analyzes the last 100 frames of `relax-300K.lammpstrj` in these runs:

- `prepared{,-rerun,-rerun2}/me-4pacz-dcz-4p-cosam`;
- `prepared{,-rerun,-rerun2}/me-4pacz-then-dcz-4p`;
- `prepared{,-rerun,-rerun2}/me-4pacz-alone`;
- the two gap-seeded DCZ-4P runs.

Missing trajectories are reported and skipped. It then runs:

- `summarize-electronic-proxy`, writing `electronic_gate.json`, `electronic_gate.csv` and `electronic_gate_runs.csv`;
- `select-patch-slots`, writing `patch_slots.csv`.

Both write to `electronic-tier0/dcz-4p/`. Per-run results go to `electronic-tier0/dcz-4p/runs/<root>/<system>/electronic-proxy-<stage>/` (`RUN_OUTPUT_ROOT`); the calculation trees are only read.

Overrides:

| Variable | Purpose |
|---|---|
| `SECONDARY` | secondary modifier slug |
| `STAGE` | trajectory stage |
| `LAST_FRAMES` | frames analyzed |
| `EXPORT_CONFORMERS` | conformers exported per component |
| `CONVENTIONAL_ROOTS` | conventional run roots |
| `INCLUDE_ME4PACZ_ALONE` | include the Me-4PACz-alone runs |
| `OUTPUT_ROOT` | output directory |
| `RUN_OUTPUT_ROOT` | per-run output directory (default `OUTPUT_ROOT/runs`) |

Before acting on the gate, resolve the comparability checks in the feasibility report: protocol equivalence across `prepared*` roots, stage-1 independence of the sequential runs, and side-view inspection of lifted DCZ-4P.

## Gate rule

For each metric:

- The run-level value is the trajectory mean (whole interface) or the median over patch candidates.
- The group difference is `mean_B − mean_A`, compared with the pooled between-run SD `sqrt((s_A² + s_B²)/2)`.
- A metric passes when the difference is at least `--threshold-sd` (default 2) times the pooled SD, with at least 2 runs per group.

`proceed_to_tier1` is true when any primary metric passes. The primary metrics are:

- whole-interface and patch-median `dphi` (charges and reference);
- whole-interface anchored-terminal densities.

With three runs per group the SD is itself uncertain, so treat marginal passes as inconclusive. A pass justifies only the Tier-1 DFT timing benchmark.

## Slot selection

For every run, `select-patch-slots` picks:

- the patch closest to the run-median anchored density;
- for runs in the compared groups (and for others when `--reference-slots 2`), the patch farthest from the run median in tilt and exposed fraction, preferring patches that contain a multi-anchored molecule.

It also picks `--matched-pairs` cross-group pairs. Each pair has anchored density within `--matched-tolerance` (default ±10%) and the largest arrangement contrast, and no two pairs come from the same run combination.

Each slot records its build directory, trajectory, frame index, step and window id, so the patch can be rebuilt for Tier 1.
