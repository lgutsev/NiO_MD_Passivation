# Coordinate-based SAM coverage analysis

`nio-md-prep analyze-coverage` replaces OVITO screenshots and color-pixel
counting with a reproducible calculation from the LAMMPS trajectory.

For each selected frame, every ligand atom is projected onto the periodic
`x/y` cell as a disk with its elemental van der Waals radius. The union of
those disks is divided by the instantaneous lateral cell area. Molecule-ID
ranges in `assembly_manifest.json` provide separate primary and secondary
component masks:

```text
total coverage = primary coverage + secondary coverage - overlap
```

This quantity intentionally matches the meaning of the historical top-view
image measurement. It is a projected shielding footprint; it is not yet a
count of chemically bonded phosphonate anchors.

## Installation

Install the optional analysis dependencies once:

```bash
python -m pip install -e '.[analysis]'
```

## Usage

Analyze the final 100 frames of the preferred trajectory found in a prepared
directory:

```bash
nio-md-prep analyze-coverage prepared/me-4pacz-alone
```

Choose a trajectory and provide its timestep when a time axis in ps is wanted:

```bash
nio-md-prep analyze-coverage \
  prepared/me-4pacz-alone \
  --trajectory equilibration-300K.lammpstrj \
  --last-frames 100 \
  --timestep-fs 1.0
```

For the compressed deposition and hold trajectories in the current protocol,
use `--timestep-fs 0.5`. Use `--last-frames 0` to analyze the complete
trajectory. `--stride 5` analyzes every fifth selected frame.

The default target grid spacing is 0.20 Å. For a numerical convergence check,
repeat the calculation with `--grid-spacing 0.10` and compare the total mean:

```bash
nio-md-prep analyze-coverage prepared/me-4pacz-alone \
  --grid-spacing 0.10 \
  --output prepared/me-4pacz-alone/coverage-analysis-0.10A
```

## Outputs

The default output directory is `BUILD_DIRECTORY/coverage-analysis`:

- `coverage_timeseries.csv`: total, uncovered, overlap, and per-component
  coverage for every analyzed frame.
- `coverage_summary.json`: mean coverage, frame standard deviation, and
  block-based standard error.
- `coverage_probability.npz`: reusable fractional-coordinate occupancy map.
- `coverage_vs_time.png`: coverage time series.
- `coverage_probability.png`: time-averaged spatial occupancy map.

The trajectory reader streams atom arrays one frame at a time. When
`--last-frames` is used, it first performs a lightweight frame-counting pass,
then reads only the requested final window into the analysis.

## Interpretation and reproducibility

The default includes hydrogen atoms and uses fixed elemental radii. All radii
can be scaled with `--radius-scale`; `--exclude-hydrogen` provides a
heavy-atom-only sensitivity check. The chosen grid spacing, radius scale,
atom-radius table, frame indices, and statistical block count are recorded in
`coverage_summary.json`.

NPT changes in `Lx` and `Ly` are handled frame by frame. The occupancy map
uses a fixed fractional grid, while physical disk radii remain in Å.
Orthogonal periodic `x/y` boxes are supported. Triclinic boxes are rejected
explicitly rather than analyzed incorrectly.

## Analyze every completed hold on a compute node

Do not analyze many trajectories on the login/head node. Submit the supplied
single-CPU job from the repository root:

```bash
sbatch scripts/analyze_coverage_holds.sbatch
```

The job searches every immediate subdirectory of `prepared/` for non-empty
`hold-300K.lammpstrj` and `hold-400K.lammpstrj` files and processes them
sequentially. Results are kept separate:

```text
prepared/SYSTEM/coverage-analysis-hold-300K/
prepared/SYSTEM/coverage-analysis-hold-400K/
```

After processing the trajectories, the job also writes a consolidated workbook:

```text
prepared/coverage_summary.xlsx
```

Its `Results` sheet contains one row per system and hold temperature, including
total coverage, uncertainty, uncovered area, component overlap, and the primary
and secondary component coverages. The `Components` sheet provides one
normalized row per SAM component, and `Methods` records the coverage definition
and quality-control checks. A bar chart on `Results` gives a quick comparison
across all completed holds.

By default, the final 100 frames of each trajectory are analyzed. To analyze
all frames:

```bash
sbatch --export=ALL,LAST_FRAMES=0 scripts/analyze_coverage_holds.sbatch
```

The grid and statistical block count can also be overridden:

```bash
sbatch --export=ALL,GRID_SPACING=0.10,BLOCKS=5 \
  scripts/analyze_coverage_holds.sbatch
```

The job uses the `single` partition with one task and one CPU, and constrains
NumPy/SciPy math libraries to one thread.

The workbook can also be rebuilt from existing result directories without
rerunning the trajectory analysis:

```bash
nio-md-prep summarize-coverage prepared \
  --output prepared/coverage_summary.xlsx
```
