# Localized 2D-void “lego2” deposition

## Purpose

The original lego control averages the Me-4PACz occupancy map over `y`, finds
the widest persistent low-occupancy interval in `x`, and initially delivers
the secondary molecules into that stripe. That remains a useful 1D access
control, but a localized hole can be diluted by the `y` average or represented
as much more empty area than actually exists.

Lego2 is an additive experimental control that instead asks:

> If the secondary molecules are delivered inside the largest persistent,
> localized two-dimensional Me-4PACz void, can they occupy it after all
> lateral constraints are removed?

It writes a separate `prepared-lego2/` tree and does not modify the existing
`prepared-lego/` workflow or results.

## Method

The builder consumes the numerical
`coverage-analysis-hold-300K/coverage_probability.npz`.

1. Mark grid cells with occupancy probability at or below `0.20`.
2. Find four-connected empty components with periodic connectivity in both
   `x` and `y`.
3. Select the largest periodic component by area.
4. Translate the completed stage-1 structure by whole coverage-grid cells so
   the component is centered in both periodic directions. LAMMPS image flags
   are updated so unwrapped coordinates remain consistent.
5. Find the largest axis-aligned rectangle lying completely inside the
   shifted component.
6. Apply a 4 Å inset on all four sides and use the remaining rectangle as the
   Packmol region for the secondary molecules.
7. Apply no lateral wall during minimization, moving-wall deposition,
   continuation, or holds.
8. Use the same two-step vertical protocol as lego: the first wall ends 30 Å
   above the completed stage-1 maximum z, then the restart-safe continuation
   lowers it another 15 Å without reinitializing velocities.

The inscribed rectangle prevents Packmol from treating covered cells inside a
component's bounding box as empty. It is still a deliberately biased initial
condition, not spontaneous self-assembly kinetics.

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

- `lego2_plan.json`: selected periodic component, two-axis translation,
  inscribed rectangle, Packmol bounds, source hashes, and wall endpoints;
- `lego2_void_map.csv`: shifted occupancy, largest-component membership, and
  packing-rectangle membership for every grid cell;
- `lego2_void_map.npz`: the same arrays for plotting or numerical checks;
- `lego2-stage1-shifted.data`: translated stage-1 reference;
- `packmol.inp`, `deposition.in`, `continue-deposition.in`; and
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

After the continuation promotes `deposited-continued.data` to
`deposited.data`, launch the independent holds:

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

Use `--array=0-2` with the same `PREPARED_ROOT` and `SYSTEM_SUFFIX` for later
stages.

## Controls

The builder accepts:

- `OCCUPANCY_THRESHOLD` (default `0.20`);
- `PACKING_INSET` (default `4.0` Å);
- `MINIMUM_VOID_AREA_FRACTION` (default `0.02`);
- `MINIMUM_PACKING_WIDTH` and `MINIMUM_PACKING_HEIGHT` (default `12.0` Å);
- `PACKMOL_SEED` and `VELOCITY_SEED`;
- `SOURCE_ROOT`, `COVERAGE_MAP`, and `PREPARED_ROOT`; and
- the same vertical-clearance and continuation overrides as lego.

If the largest inscribed rectangle is too small for the requested minimum,
the builder stops rather than silently placing molecules over covered cells.
Even when the geometry passes, Packmol can reject an unrealistically dense
request. That outcome means the localized void cannot accommodate the full
secondary inventory under the chosen threshold and inset.

## Interpretation

Compare three matched calculations:

| Calculation | Initial secondary delivery | Question |
|---|---|---|
| Standard sequential | Full lateral region | Does the molecule find and occupy gaps without guidance? |
| Lego | Widest low-occupancy x stripe | Can it occupy a broad access channel? |
| Lego2 | Largest localized periodic 2D void | Can it occupy the actual largest persistent hole? |

An improved lego2 result supports access/search as a plausible limitation. It
does not establish unbiased kinetics, adsorption thermodynamics, proton
transfer, or electronic coupling.
