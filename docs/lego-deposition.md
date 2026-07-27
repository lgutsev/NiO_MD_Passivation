# Coverage-guided “lego” deposition

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
6. During the usual stage-2 SD/CG minimization and moving-upper-wall
   deposition, apply a harmonic `fix wall/region` side wall to the `stage2`
   group only. NiO and Me-4PACz receive no lateral guiding force.
7. Remove the lateral wall after deposition and before writing
   `deposited.data`. Compressed holds, decompression, and relaxed holds use the
   ordinary unbiased inputs.

The default harmonic wall is 0.50 kcal mol⁻¹ Å⁻² with a 3 Å cutoff. Packmol's
4 Å inset ensures that the initial layer-2 coordinates start outside the
wall-force zone. The lower recoil wall remains active exactly as in the
standard sequential workflow.

The generated syntax follows the LAMMPS
[`fix wall/region`](https://docs.lammps.org/fix_wall_region.html) and
[`region`](https://docs.lammps.org/region.html) definitions. The `y` and `z`
block faces are open, leaving only the two lateral `x` faces active.

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
prepared-lego/me-4pacz-then-dcz-4p-lego/
```

No simulation is submitted automatically. Before deposition, inspect:

- `lego_plan.json`: detected gap, periodic shift, tunnel bounds, wall
  parameters, and source-map hash;
- `lego_tunnel_profile.csv`: mean occupancy for every `x` column and the
  columns selected as the periodic gap;
- `packmol.inp`: the restricted stage-2 packing box;
- `lego-stage1-shifted.data`: the periodically translated stage-1 reference;
- `deposition.in`: the stage2-only tunnel and its removal after deposition;
- `validation_report.txt`: normal assembly validation.

Launch only the DCZ-4P deposition:

```bash
sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego \
  scripts/run_sequential_deposition_array.sbatch
```

After `deposited.data` is complete, the unbiased 300 K and 400 K branches use
the same suffix mechanism:

```bash
sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego \
  scripts/run_sequential_hold_array.sbatch

sbatch --array=0 \
  --export=ALL,PREPARED_ROOT=prepared-lego,SYSTEM_SUFFIX=-lego \
  scripts/run_sequential_hold_400K_array.sbatch
```

## Controls and interpretation

Compare at least these two systems using the same stage-1 source, molecular
inventory, Packmol seed, LAMMPS velocity seed, and downstream protocol:

| Calculation | Lateral delivery | Question answered |
|---|---|---|
| Standard sequential DCZ-4P | Unbiased full `x/y` region | Does DCZ-4P find and occupy gaps spontaneously in this model? |
| Lego-style sequential DCZ-4P | Coverage-guided `x` channel | Can DCZ-4P occupy the gap when it is delivered there? |

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
`SYSTEM_SUFFIX=-lego` arguments for each later sequential runner.

## Parameter overrides

The builder accepts environment overrides without changing the study TOML:

```bash
sbatch --array=0 \
  --export=ALL,OCCUPANCY_THRESHOLD=0.20,PACKING_INSET=4.0,LEGO_WALL_STRENGTH=0.50,LEGO_WALL_CUTOFF=3.0 \
  scripts/build_lego_systems.sbatch
```

Other available controls are `COVERAGE_MAP`, `SOURCE_ROOT`, `PREPARED_ROOT`,
`MINIMUM_GAP_FRACTION`, `MINIMUM_TUNNEL_WIDTH`, `PACKMOL_SEED`, and
`VELOCITY_SEED`. Record any non-default values in the methods and treat
threshold or wall-strength scans as sensitivity tests rather than independent
replicates.
