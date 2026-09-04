# Operations, validation, and maintenance

This page collects repository-wide operational procedures that do not belong in the project overview.

## Full campaign inventory

Generate one workbook covering every run under every `prepared*` tree:

```bash
python scripts/inventory_runs.py
```

The command writes `run_inventory.xlsx` in the repository root. It inventories actual structure, trajectory, restart, log, marker, and analysis-report files; it does not infer completion from a Slurm exit code alone.

The workbook contains a campaign dashboard, one row per run, one row per stage, a prioritized action queue, analysis-report coverage, and a snapshot of the current user's Slurm jobs.

A completed hold loop with a restart but no `held-300K.data` or `held-400K.data` is marked `FINALIZE_NEEDED`, while missing DCZ control inputs are marked `INPUT_MISSING`. This makes failed or skipped array elements visible without manually opening every directory.

Select trees or choose the output filename explicitly:

```bash
python scripts/inventory_runs.py \
  --prepared-root prepared \
  --prepared-root prepared-rerun \
  --output status.xlsx
```

## Validation and safe input refresh

`nio-md-prep` checks atom and topology counts, type remapping, charge, coordinate bounds, molecule identities, force-field coverage, surface provenance, wall placement, stage ordering, and—when LAMMPS is available—zero-step executable inputs.

Build provenance and hashes are recorded in `assembly_manifest.json` and `input_hashes.json`.

Validation can be rerun as predecessor files become available:

```bash
nio-md-prep validate prepared/SYSTEM_NAME
```

If topology is already prepared and only generated LAMMPS inputs need an update, use:

```bash
nio-md-prep refresh-inputs studies/STUDY_NAME.toml \
    --output prepared/SYSTEM_NAME
```

`refresh-inputs` regenerates stage inputs and protocol notes only. It does not run Packmol, rebuild topology, or modify simulation data, restarts, logs, or trajectories.

## Archiving completed and restartable runs

Generated `prepared/` trees are not stored in Git. To retain every LAMMPS `.data` state and the small inputs needed to rerun it, create a ZIP archive:

```bash
nio-md-prep archive-runs prepared prepared-rerun prepared-rerun2 \
    --output archives/nio-md-rerun-files.zip
```

On QBD, the single-core archive launcher automatically includes every top-level directory matching `prepared*`:

```bash
sbatch scripts/archive_prepared_runs.sbatch
```

Select roots or the output filename explicitly:

```bash
sbatch --export=ALL,\
ARCHIVE_ROOTS="prepared-rerun prepared-rerun2",\
ARCHIVE_OUTPUT=archives/cosam-reruns.zip \
scripts/archive_prepared_runs.sbatch
```

The ZIP preserves the prepared-root and system-directory layout. It contains all `.data`, `.in`, `.lmp`, `.inp`, `.toml`, `.json`, `.txt`, and `.xyz` files, but excludes trajectories, logs, Slurm output, restart binaries, and analysis workbooks.

`ARCHIVE_MANIFEST.json` records the repository commit and SHA-256 checksum of every collected file. `README_RERUN.txt` provides minimal extraction and LAMMPS rerun instructions. Existing archives are not overwritten unless the direct CLI is given `--force`.

## Adding another passivant

1. Generate and download the LigParGen LAMMPS file independently.
2. Run `nio-md-prep init-molecule NAME`.
3. Place the downloaded file at the displayed path as `ligpargen.lmp`.
4. Review `molecule.toml`, especially expected charge, role, and number of phosphonic-acid anchors.
5. Run `nio-md-prep inspect-molecule NAME`.
6. Add a study TOML with explicit molecule counts.
7. Build and validate the study.

If Packmol is unavailable during a local build, the package retains `packmol.inp`. Run

```bash
packmol < packmol.inp
```

and then repeat the build with `--packed-xyz packed.xyz`.

## Why Excel is no longer part of system preparation

The historical workflow used a spreadsheet to paste Packmol coordinates into replicated LigParGen atom records, shift atom/molecule/topology IDs, append the NiO surface, and repair the LAMMPS header manually. This was a data-joining step rather than part of the scientific method.

The current assembler performs that operation programmatically:

- Packmol supplies coordinates and element order;
- LigParGen supplies charges, atom types, masses, coefficients, and bonded topology;
- the packaged NiO file supplies the unchanged surface;
- `nio-md-prep` verifies element order, remaps every namespace, assigns unique molecule IDs, and calculates counts, bounds, and total charge.

The legacy notebooks and historical examples remain as scientific provenance and a record of how the original calculation was constructed. They are not required for new builds.

## Analysis workbooks

The main analysis launchers build consolidated workbooks after trajectory processing:

```bash
sbatch scripts/analyze_coverage_holds.sbatch
sbatch scripts/analyze_interface_holds.sbatch
sbatch scripts/analyze_interface_depositions.sbatch
```

Typical outputs include:

```text
prepared/coverage_summary.xlsx
prepared/interface_structure_summary.xlsx
prepared/publication_summary.xlsx
```

The publication summary can optionally join private experimental device data. Keep unpublished device CSV files outside the repository working tree. The `.gitignore` safeguards are secondary; the primary rule is not to place private experimental data under version control.

See [`coverage-analysis.md`](coverage-analysis.md) and [`interfacial-analysis.md`](interfacial-analysis.md) for definitions, interpretation limits, and workbook contents.

## Documentation maintenance

The root README should remain a project overview. When a workflow becomes longer than a short quick-start example, add it to the relevant page in this directory and link it from [`docs/README.md`](README.md). This keeps user-facing entry points stable while allowing detailed procedures to evolve independently.
