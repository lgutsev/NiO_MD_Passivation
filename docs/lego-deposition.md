# Coverage-guided, laterally free “lego” deposition

## Scientific purpose

The standard sequential calculation asks what the secondary SAM does without
lateral guidance after an equilibrated Me-4PACz layer has formed. The
lego-style calculation asks a narrower counterfactual question:

> If the secondary molecules are delivered directly over a persistent
> Me-4PACz coverage gap, can they enter and stabilize that gap?

This is therefore a deliberately biased geometric-control experiment. It can
distinguish failure to find the gap from failure to occupy the gap, but it
cannot be presented as spontaneous self-assembly kinetics.

## Method

The builder consumes the numerical
`coverage-analysis-hold-300K/coverage_probability.npz` produced by
`nio-md-prep analyze-coverage`. The rendered PNG is not used.

1. Average occupancy probability over fractional `y` to obtain an `x`
   occupancy profile.
2. Mark columns with mean occupancy at or below 0.20 as candidate empty
   columns.
3. Find the widest circular run of candidate columns. The circular search is
   essential because a gap split between the left and right plot edges is one
   physical gap under periodic boundary conditions.
4. Shift all stage-1 coordinates by a periodic `x` translation so that the
   gap is centered at fractional `x = 0.5`. Image flags are updated so that
   unwrapped coordinates remain consistent.
5. Pack only the secondary molecules inside the centered channel, with a
   4 Å inset from each side.
6. Apply no lateral wall during minimization or dynamics. The coverage map
   biases only the initial Packmol placement; stage 2 can spread freely in
   periodic x/y throughout deposition and all later stages.
7. Move the upper z-wall first to 30 Å above the maximum z of the completed
   Me-4PACz film. The initial wall position is measured after minimization and
   frozen at least 3.5 Å above the highest atom, avoiding atoms on or inside
   the wall at step 0.
8. In a separate restart-safe stage, lower the wall another 15 Å over 300,000
   steps (150 ps) without reinitializing velocities. The final compressed
   endpoint is therefore 15 Å above the stage-1 maximum z.

The lower recoil wall remains active exactly as in the standard sequential
workflow. The former harmonic tunnel remains reproducible as an explicit
confined control using `--lateral-wall`, but it is no longer the default.

## Build the first DCZ-4P control

Prerequisites:

- `prepared/me-4pacz-alone/held-300K.data`
- `prepared/me-4pacz-alone/coverage-analysis-hold-300K/coverage_probability.npz`
- the analysis dependencies installed with `python -m pip install -e
  '.[analysis]'`

Submit from the repository root:

```bash
sbatch scripts/build_lego_systems.sbatch
```

The script defaults to array index 0, DCZ-4P, and writes:

```text
prepared-lego/me-4pacz-then-dcz-4p-lego-seeded/
```

No simulation is submitted automatically. Before deposition, inspect:

- `lego_plan.json`: detected gap, periodic shift, packing bounds, stage-1
  maximum z, resolved deposition endpoint, and source-map hash;
- `lego_tunnel_profile.csv`: mean occupancy for every `x` column and the
  columns selected as the periodic gap;
- `packmol.inp`: the restricted stage-2 packing box;
- `lego-stage1-shifted.data`: the periodically translated stage-1 reference;
- `deposition.in`: laterally free stage-2 deposition with a measured safe
  initial wall position;
- `continue-deposition.in`: final 15 Å wall descent without velocity
  reinitialization;
- `validation_report.txt`: normal assembly validation.

Launch only the DCZ-4P deposition:

```bash
sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego-seeded \
  scripts/run_sequential_deposition_array.sbatch
```

Then complete the final gentle compression:

```bash
sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego-seeded \
  scripts/run_lego_deposition_continuation_array.sbatch
```

LAMMPS writes `deposited-continued.data`; only after successful completion
does the runner preserve the old state as
`deposited.pre-continuation.JOBID.data` and promote the new state to
`deposited.data`. The ordinary hold scripts therefore remain unchanged.

After the continuation completes, the unbiased 300 K and 400 K branches use
the same suffix mechanism:

```bash
sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego-seeded \
  scripts/run_sequential_hold_array.sbatch

sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego-seeded \
  scripts/run_sequential_hold_400K_array.sbatch
```

## Controls and interpretation

Compare at least these two systems using the same stage-1 source, molecular
inventory, Packmol seed, LAMMPS velocity seed, and downstream protocol:

| Calculation | Lateral delivery | Question answered |
|---|---|---|
| Standard sequential DCZ-4P | Unbiased full `x/y` region | Does DCZ-4P find and occupy gaps spontaneously in this model? |
| Lego-style sequential DCZ-4P | Coverage-guided initial `x` placement, laterally free dynamics | Can DCZ-4P occupy the gap when initially delivered there? |

Interpret the paired outcomes conservatively:

| Standard | Lego | Supported interpretation |
|---|---|---|
| weak | improved gap filling | access/search is a plausible limitation |
| weak | weak | delivery alone is insufficient; affinity, orientation, or steric compatibility may limit insertion |
| strong | strong | no special guidance is needed in the sampled trajectory |
| strong | weak | inspect tunnel width, forces, and initial packing; the bias may be obstructive |

The comparison should include secondary Ni-contact/anchor occupancy,
secondary height and orientation, persistent contacts, projected uncovered
area, and whether DCZ-4P remains above Me-4PACz rather than reaching NiO.

## Optional extension to the other secondary molecules

The system order matches the ordinary sequential scripts:

| Array index | Molecule |
|---:|---|
| 0 | DCZ-4P |
| 1 | MeO-2PACz |
| 2 | MeO-4PADBC |

After the DCZ-4P pilot is inspected, all three can be built with:

```bash
sbatch --array=0-2 scripts/build_lego_systems.sbatch
```

Use the same `--array=0-2`, `PREPARED_ROOT=prepared-lego`, and
`SYSTEM_SUFFIX=-lego-seeded` arguments for each later sequential runner.

## Parameter overrides

The builder accepts environment overrides without changing the study TOML:

```bash
sbatch --array=0 \
  --export=ALL,OCCUPANCY_THRESHOLD=0.20,PACKING_INSET=4.0,LEGO_DEPOSITION_CLEARANCE=30.0,LEGO_FINAL_DEPOSITION_CLEARANCE=15.0,LEGO_CONTINUATION_STEPS=300000 \
  scripts/build_lego_systems.sbatch
```

Other available controls are `COVERAGE_MAP`, `SOURCE_ROOT`, `PREPARED_ROOT`,
`LEGO_SUFFIX`, `MINIMUM_GAP_FRACTION`, `MINIMUM_TUNNEL_WIDTH`,
`PACKMOL_SEED`, `VELOCITY_SEED`, and `LEGO_ADDITIONAL_DROP` for upgrading
existing lego folders. Set `LEGO_LATERAL_WALL=1` only to
reproduce the former confined control; `LEGO_WALL_STRENGTH` and
`LEGO_WALL_CUTOFF` then set that optional wall. Record any non-default values
in the methods and treat threshold, clearance, or wall-strength scans as
sensitivity tests rather than independent replicates.

## Upgrading existing lego outputs

The continuation preparation job rewrites only input and provenance files; it
does not alter `topology_output.lmp`, `deposited.data`, trajectories, or
restarts:

```bash
sbatch --array=0-2 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego-seeded \
  scripts/prepare_lego_continuation_inputs.sbatch
```

This also repairs a failed initial deposition by regenerating `deposition.in`
with the measured safe wall start. Such a case can be resubmitted with the
ordinary sequential deposition runner before its continuation is launched.
