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
