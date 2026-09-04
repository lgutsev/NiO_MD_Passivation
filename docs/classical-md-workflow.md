# Classical LAMMPS workflow

This page contains the production workflow that was previously embedded in the root README.

## Installation

Python 3.11 or newer and Packmol are required for building systems. LAMMPS is required for validation smoke tests and simulation.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Optional analysis and test dependencies:

```bash
python -m pip install -e '.[analysis]'
python -m pip install -e '.[test]'
pytest
```

On QBD, the supplied launchers load the current Python and LAMMPS modules. Build jobs use the `single` partition. Molecular-dynamics jobs use one complete 64-core `workq` node with one MPI rank per core.

## Building a primary or CoSAM system

Each study is defined in `studies/*.toml`. The configuration records molecule identities and counts, deposition and hold lengths, wall policy, and composition basis.

Build one system manually:

```bash
nio-md-prep build studies/me-4pacz-alone.toml \
    --output prepared/me-4pacz-alone
nio-md-prep validate prepared/me-4pacz-alone
```

On QBD, prepare and validate the five primary-deposition systems:

```bash
sbatch scripts/build_real_systems.sbatch
```

The array order is fixed across later primary/CoSAM stages:

| Array index | Prepared system |
|---:|---|
| 0 | `me-4pacz-alone` |
| 1 | `me-4pacz-meo-2pacz-cosam` |
| 2 | `me-4pacz-meo-4padbc-cosam` |
| 3 | `me-4pacz-dcz-4p-cosam` |
| 4 | `me-4pacz-high-dose` |

The workflow intentionally uses manual stages rather than an automatic dependency chain. After inspecting the build, submit deposition:

```bash
sbatch scripts/run_deposition_array.sbatch
```

After every required `deposited.data` exists, the temperature branches are independent and may run concurrently:

```bash
sbatch scripts/run_hold_array.sbatch
sbatch scripts/run_hold_400K_array.sbatch
```

The 400 K branch does not consume or modify `held-300K.data`; both branches read the same `deposited.data`.

## Sequential second-layer deposition

Sequential preparation is a genuine two-stage calculation, not a differently labeled CoSAM.

First, Me-4PACz is deposited by itself and its compressed 300 K hold is completed. The exact

```text
prepared/me-4pacz-alone/held-300K.data
```

state becomes the substrate for all three sequential systems. The compressed state is used rather than the later expanded-box continuation so the established film is preserved while sufficient vacuum remains above it for packing layer 2.

During stage-2 preparation:

1. every stage-1 coordinate, topology record, velocity, and box bound is preserved;
2. Packmol places only the secondary molecules above the existing film;
3. the completed Me-4PACz layer is fixed during a conservative steepest-descent/conjugate-gradient relaxation;
4. only newly added molecules receive fresh 300 K velocities;
5. the lock is removed and the complete system evolves together during moving-wall deposition.

Build after the Me-4PACz 300 K hold is available:

```bash
sbatch scripts/build_sequential_systems.sbatch
```

Sequential array order:

| Array index | Prepared system |
|---:|---|
| 0 | `me-4pacz-then-dcz-4p` |
| 1 | `me-4pacz-then-meo-2pacz` |
| 2 | `me-4pacz-then-meo-4padbc` |

Run deposition and then either or both temperature branches:

```bash
sbatch scripts/run_sequential_deposition_array.sbatch
sbatch scripts/run_sequential_hold_array.sbatch
sbatch scripts/run_sequential_hold_400K_array.sbatch
```

No builder submits a simulation job and no deposition job submits a hold job. Each stage stops after its own work so its output can be inspected before the next scientific step.

## Experimental seeded-accessibility controls

The optional LEGO controls deliberately place secondary material into an identified gap. They test geometric accessibility and should not be interpreted as evidence that the secondary species would spontaneously find that gap in an unbiased simulation.

- [`lego-deposition.md`](lego-deposition.md): 1D coverage-guided channel seeding.
- [`lego2-deposition.md`](lego2-deposition.md): periodic 2D void-component seeding.

Default builders:

```bash
sbatch scripts/build_lego_systems.sbatch
sbatch scripts/build_lego2_systems.sbatch
```

Review the generated planning JSON/CSV files before launching the corresponding sequential deposition runners.

## Simulation stages and wall policy

Every completed build writes `deposition.in`, `hold-300K.in`, `hold-400K.in`, `decompress-300K.in`, `decompress-400K.in`, `equilibrate-300K.in`, and `anneal-400K.in`.

| Stage | Source | Default duration | Main output |
|---|---|---:|---|
| Moving-wall deposition | `topology_output.lmp` | 600,000 × 0.5 fs = 300 ps | `deposited.data` |
| Compressed 300 K hold | `deposited.data` | 1,000,000 × 0.5 fs = 500 ps | `held-300K.data` |
| Gradual heating | `deposited.data` | 400,000 × 0.5 fs = 200 ps | `heated-400K.data` |
| Compressed 400 K hold | heated 400 K state | 1,000,000 × 0.5 fs = 500 ps | `held-400K.data` |
| 300 K wall retraction | `held-300K.data` | 400,000 × 0.5 fs = 200 ps | `decompressed-300K.data` |
| Relaxed 300 K hold | decompressed 300 K state | 1,000,000 × 0.5 fs = 500 ps | `relaxed-300K.data` |
| 400 K wall retraction | `held-400K.data` | 400,000 × 0.5 fs = 200 ps | `decompressed-400K.data` |
| Relaxed 400 K hold | decompressed 400 K state | 1,000,000 × 0.5 fs = 500 ps | `relaxed-400K.data` |
| Optional long 300 K continuation | `held-300K.data` | 5,000,000 steps | `equilibrated-300K.data` |
| Optional long 400 K continuation | `equilibrated-300K.data` | 3,000,000 steps | `annealed-400K.data` |

Primary and CoSAM depositions use the established displacement-limited CG relaxation. Sequential depositions use the more conservative SD/CG sequence while the completed first layer is fixed. Both write `optimized.data`, `optimized.restart`, `optimized.lammpstrj`, and `optimization-summary.txt` before dynamics.

The moving upper wall ends 30 Å above the NiO reference surface: `69.615 Å` for the supplied model. Its instantaneous location is written as `v_zwall` in deposition thermo output. Compressed holds keep the upper wall fixed at that endpoint. Lower recoil walls use LAMMPS `EDGE`.

The optional long continuations resolve a less restrictive upper wall from the preceding structure:

```text
ceil(max(120 Å, max_atom_z + 30 Å) / 10 Å) × 10 Å
```

Only the nonperiodic cell ceiling is expanded, retaining a 5 Å margin. A numeric `production_upper_wall` in the study file overrides the automatic policy.

## Relaxing the compressed films

The decompression branches preserve the compressed results and operate independently at 300 and 400 K. Each branch retracts the upper wall to the resolved production height over 200 ps and then keeps the distant wall fixed for a 500 ps relaxed-film hold. Only the periodic `x/y` dimensions are barostatted.

Regenerate only the decompression-stage inputs for all eight prepared systems:

```bash
sbatch scripts/build_decompression_inputs.sbatch
```

This does not rebuild topology or modify data, restart, trajectory, or log files. Then submit:

```bash
sbatch scripts/run_decompress_300K_array.sbatch
sbatch scripts/run_decompress_400K_array.sbatch
```

Common decompression array order:

| Array index | Prepared system |
|---:|---|
| 0 | `me-4pacz-alone` |
| 1 | `me-4pacz-meo-2pacz-cosam` |
| 2 | `me-4pacz-meo-4padbc-cosam` |
| 3 | `me-4pacz-dcz-4p-cosam` |
| 4 | `me-4pacz-high-dose` |
| 5 | `me-4pacz-then-dcz-4p` |
| 6 | `me-4pacz-then-meo-2pacz` |
| 7 | `me-4pacz-then-meo-4padbc` |

Final trajectories are `relax-300K.lammpstrj` and `relax-400K.lammpstrj`.

## DCZ-4P upper-wall controls

DCZ-4P can form a cohesive aggregate that follows the attractive part of the ordinary 12-6 upper wall during retraction. Two isolated controls test the wall contribution without replacing production results:

- `nonsticky`: same wall trajectory, but the 12-6 potential is truncated at `2^(1/6) sigma`, making it WCA-like and purely repulsive;
- `nowall`: upper wall is removed immediately while the lower recoil wall and observation windows are retained.

Prepare controls for sequential DCZ-4P in the supported prepared trees:

```bash
sbatch scripts/build_dcz_decompression_controls.sbatch
sbatch scripts/run_dcz_decompression_controls_array.sbatch
```

Control outputs use `relaxed-control-{nonsticky,nowall}-{300K,400K}.data` and matching trajectory/log labels. Existing production `decompressed-*` and `relaxed-*` files are never modified.

## Independent replicate trees

Builders and runners accept `PREPARED_ROOT`. New random seeds should be supplied only to build stages. Example seed-02:

```bash
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02,PACKMOL_SEED=202607242,VELOCITY_SEED=31415927 scripts/build_real_systems.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_deposition_array.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_hold_array.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_hold_400K_array.sbatch

sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02,PACKMOL_SEED=202608242,VELOCITY_SEED=31416927 scripts/build_sequential_systems.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_sequential_deposition_array.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_sequential_hold_array.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_sequential_hold_400K_array.sbatch
```

The 300 and 400 K branches are independent from each replicate's `deposited.data` and may run concurrently after deposition.

## Analysis after production

Use [`coverage-analysis.md`](coverage-analysis.md) for projected coverage, void topology, and roughness, and [`interfacial-analysis.md`](interfacial-analysis.md) for anchor-resolved observables, site ownership, orientation, and publication summary generation.
