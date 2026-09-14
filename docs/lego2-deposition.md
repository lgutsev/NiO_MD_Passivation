# Gap-seeded II: localized periodic 2D-void accessibility control

## Scientific purpose

**Gap-seeded I** reduces the access problem to one dimension by averaging the Me-4PACz occupancy map over `y`, identifying the widest persistent low-occupancy interval in `x`, and initially delivering the secondary molecules into that stripe. This is useful for testing access through a broad channel, but a localized hole can be diluted by the `y` average or represented as more available area than is actually present.

**Gap-seeded II** therefore asks a more spatially resolved counterfactual question:

> If the secondary molecules are initially delivered inside the largest persistent localized two-dimensional void in the Me-4PACz film, can they occupy that region once all lateral constraints are removed?

This is a deliberately biased geometric accessibility control. It tests whether an identified persistent void is physically occupiable when the lateral search problem is reduced; it does not represent spontaneous self-assembly kinetics.

The implementation retains the historical internal identifier `lego2` in script names, directories, variables, and provenance files for compatibility with existing production data.

## Method

The builder consumes the numerical
`coverage-analysis-hold-300K/coverage_probability.npz`.

1. Mark grid cells with occupancy probability at or below `0.20`.
2. Find four-connected empty components with periodic connectivity in both `x` and `y`.
3. Select the largest periodic component by area.
4. Translate the completed stage-1 structure by whole coverage-grid cells so that the selected component is centered in both periodic directions. LAMMPS image flags are updated so unwrapped coordinates remain consistent.
5. Find the largest axis-aligned rectangle lying completely inside the shifted component.
6. Apply a 4 Å inset on all four sides and use the remaining rectangle as the Packmol region for the secondary molecules.
7. Apply no lateral wall during minimization, moving-wall deposition, continuation, or holds.
8. Use the same two-step vertical protocol as Gap-seeded I: the first wall ends 30 Å above the completed stage-1 maximum `z`, then the restart-safe continuation lowers it by another 15 Å without reinitializing velocities.

The inscribed rectangle prevents Packmol from treating covered cells inside the component's bounding box as empty. The resulting system is still a controlled initial-condition test rather than a prediction of unbiased assembly kinetics.

## Build the DCZ-4P pilot

Prerequisites:

- `prepared/me-4pacz-alone/held-300K.data`
- `prepared/me-4pacz-alone/coverage-analysis-hold-300K/coverage_probability.npz`
- `python -m pip install -e '.[analysis]'`

Submit:

```bash
sbatch scripts/build_lego2_systems.sbatch
```

The default array contains only index 0 and writes:

```text
prepared-lego2/me-4pacz-then-dcz-4p-lego2-seeded/
```

Inspect these files before launching dynamics:

- `lego2_plan.json`: selected periodic component, two-axis translation, inscribed rectangle, Packmol bounds, source hashes, and wall endpoints;
- `lego2_void_map.csv`: shifted occupancy, largest-component membership, and packing-rectangle membership for every grid cell;
- `lego2_void_map.npz`: the same arrays for plotting or numerical checks;
- `lego2-stage1-shifted.data`: translated stage-1 reference;
- `packmol.inp`, `deposition.in`, and `continue-deposition.in`;
- `validation_report.txt`.

Launch deposition and the required final compression:

```bash
sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego2,SYSTEM_SUFFIX=-lego2-seeded \
  scripts/run_sequential_deposition_array.sbatch

sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego2,SYSTEM_SUFFIX=-lego2-seeded \
  scripts/run_lego_deposition_continuation_array.sbatch
```

After the continuation promotes `deposited-continued.data` to `deposited.data`, launch the independent holds:

```bash
sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego2,SYSTEM_SUFFIX=-lego2-seeded \
  scripts/run_sequential_hold_array.sbatch

sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego2,SYSTEM_SUFFIX=-lego2-seeded \
  scripts/run_sequential_hold_400K_array.sbatch
```

## Extend to all secondary molecules

The array order is unchanged:

| Index | Secondary molecule |
|---:|---|
| 0 | DCZ-4P |
| 1 | MeO-2PACz |
| 2 | MeO-4PADBC |

After inspecting the DCZ-4P pilot:

```bash
sbatch --array=0-2 scripts/build_lego2_systems.sbatch
```

Use `--array=0-2` with the same `PREPARED_ROOT` and `SYSTEM_SUFFIX` for later stages.

## Controls

The builder accepts:

- `OCCUPANCY_THRESHOLD` (default `0.20`);
- `PACKING_INSET` (default `4.0` Å);
- `MINIMUM_VOID_AREA_FRACTION` (default `0.02`);
- `MINIMUM_PACKING_WIDTH` and `MINIMUM_PACKING_HEIGHT` (default `12.0` Å);
- `PACKMOL_SEED` and `VELOCITY_SEED`;
- `SOURCE_ROOT`, `COVERAGE_MAP`, and `PREPARED_ROOT`;
- the same vertical-clearance and continuation overrides as Gap-seeded I.

If the largest inscribed rectangle is too small for the requested minimum, the builder stops rather than silently placing molecules over covered cells. Even when the geometry passes, Packmol can reject an unrealistically dense request. That outcome indicates that the localized void cannot accommodate the full secondary inventory under the selected threshold and inset.

## Interpretation

Compare three matched calculations:

| Calculation | Initial secondary delivery | Question |
|---|---|---|
| Standard sequential | Full lateral region | Does the molecule find and occupy gaps without guidance? |
| Gap-seeded I | Widest persistent low-occupancy `x` stripe | Can it occupy a broad access channel when initially delivered there? |
| Gap-seeded II | Largest persistent localized periodic 2D void | Can it occupy the actual largest persistent hole when initially delivered there? |

Improved occupation in Gap-seeded II supports limited access/search as a plausible contributor to incomplete sequential coverage. A weak result despite targeted placement suggests that access alone is insufficient and that steric compatibility, orientation, or the effective interaction model may instead limit insertion.

Neither outcome establishes adsorption thermodynamics, proton transfer, bond exchange, electronic coupling, or unbiased kinetic rates.
