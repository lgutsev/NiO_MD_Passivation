# Anchor-resolved interfacial structure analysis

`nio-md-prep analyze-interface` complements projected coverage with
surface-normal and connectivity-aware observables. It is intended to
distinguish three physically different outcomes that a top-view union cannot:

1. a secondary molecule displaces a primary phosphonate from NiO;
2. it fills an exposed Ni site without disturbing the primary film;
3. it remains above the primary film as an overlapping second layer.

The analysis is separate from `analyze-coverage` and writes a separate
workbook. Projected coverage remains useful as a shielding footprint, but it
is not reused as a chemical-anchoring criterion.

## Installation

The existing analysis extra contains every required package:

```bash
python -m pip install -e '.[analysis]'
```

## Single-trajectory usage

Analyze the final 100 frames:

```bash
nio-md-prep analyze-interface prepared/me-4pacz-then-dcz-4p \
    --trajectory hold-300K.lammpstrj \
    --last-frames 100 \
    --timestep-fs 0.5
```

The default contact cutoff is the common value `3.25 Å`. The workbook also
reports bound fractions at `3.0`, `3.25`, and `3.5 Å`, so qualitative
conclusions can be checked against the contact definition without rerunning
the trajectory analysis. The analyzer still records the nearest-distance
distribution and an automatically inferred minimum as diagnostics.

The main cutoff and sensitivity grid can be overridden:

```bash
nio-md-prep analyze-interface prepared/me-4pacz-then-dcz-4p \
    --trajectory hold-300K.lammpstrj \
    --contact-cutoff 3.25 \
    --sensitivity-cutoffs 3.0,3.25,3.5
```

Important controls include:

- `--last-frames 0`: analyze every trajectory frame.
- `--stride N`: analyze every Nth selected frame.
- `--blocks N`: block count for within-trajectory SEM.
- `--persistence-threshold 0.50`: fraction of frames required for a molecule
  to be classified as persistently bound.
- `--surface-coordination-cutoff 2.8`: Ni-O neighbor cutoff used once on the
  canonical pristine surface to define the fixed exposed-Ni identities.
- `--rdf-rmax 15.0` and `--rdf-bin-width 0.25`: lateral RDF range and binning.
- `--z-min`, `--z-max`, and `--z-bin-width`: local-height density grid.

`rdf-rmax` must not exceed half the shortest lateral box length.

## Definitions

### Exposed Ni site

The authoritative pristine corrugated-NiO topology is used once to select
molecule-0 Ni atoms in the upper half of the slab with fewer than six
molecule-0 O neighbors. Those canonical atom identities are mapped by
molecule-0 atom order into every CoSAM, sequential, compressed, relaxed, and
independent-seed topology. Distances use the actual coordinates in each
trajectory, but the chemical site identities and denominator remain fixed.

This avoids reclassifying thermally distorted bulk Ni atoms as new surface
sites in sequential systems.

### Phosphonate terminal and bound molecule

Every P atom and its three directly bonded O atoms are identified from the
LAMMPS bond graph. A terminal is anchored in a frame when at least one of
those O atoms lies within the recorded contact cutoff of an exposed Ni site.
A molecule is bound when one or more of its terminals are anchored.

For DCZ-4P, the `Terminal States` output therefore reports the populations
with zero, one, or two anchored phosphonate termini. This distinguishes an
unbound overlayer, a one-ended interfacial modifier, and a two-ended bridge.

### Ni-site ownership

Every exposed Ni site is classified per frame as:

- empty;
- contacted only by one named component; or
- shared by two or more components.

The workbook quality check verifies that empty, shared, and component-only
site fractions sum to 100%.

### Tilt and orientational order

The molecular core is selected reproducibly from the C/N atoms farthest in
the bond graph from the molecule's phosphonate P atom(s). The tilt angle is
formed between the phosphonate-P-to-core vector and the surface normal:

```text
0 degrees = upright
90 degrees = flat
```

The corresponding order parameter is:

```text
P2 = (3 cos(theta)^2 - 1) / 2
```

`P2 = 1` is upright and `P2 = -0.5` is flat. DCZ-4P has two terminal P atoms;
their centroid is used as the anchor-side reference, so its tilt should be
interpreted together with the zero/one/two-terminal populations.

The component table additionally reports:

- the normal of a least-squares plane through the graph-distant molecular
  core (`0°` means that plane is parallel to the surface);
- the P--P anchor-axis tilt for multi-phosphonate molecules;
- the absolute vertical separation of the two anchor P atoms; and
- bound-only and unbound-only height and tilt statistics.

These descriptors are preferable to the single P-to-core vector when
diagnosing whether bifunctional DCZ-4P lies above the primary SAM.

### Heights and density profiles

Phosphorus, graph-derived molecular-core, and atom-density heights are
referenced to the nearest exposed Ni site in periodic x/y. This local
reference follows the corrugation and is preferable to subtracting one global
surface z coordinate.

The z profile reports component-resolved heavy-atom and phosphorus number
density in inverse cubic angstroms.

### Lateral order

Molecular core centers are used for the two-dimensional same-component and
cross-component radial distribution functions. Pair counts are normalized by
the instantaneous lateral area. Mean nearest-neighbor distances and the
orientational `P2` value provide compact workbook-level ordering descriptors.

### Persistence and residence

For each molecule the analyzer retains the bound/unbound sequence across the
selected frames. The summary contains:

- mean contact occupancy;
- fraction of persistently bound molecules;
- number of continuous bound episodes;
- mean and maximum episode lengths in frames; and
- durations in ps when `--timestep-fs` is supplied.

These are trajectory-sampling measures. Block SEM and residence statistics
from one trajectory are not substitutes for independent assembly seeds.

## Per-trajectory outputs

The output directory contains:

- `interface_summary.json`: definitions, parameters, component summaries, and
  quality-control metadata;
- `interface_timeseries.csv`: frame-resolved site ownership, binding, tilt,
  order, height, nearest-neighbor, and terminal-state values;
- `contact_distance_histogram.csv`: nearest Ni-phosphonate-O distribution used
  to audit the automatic cutoff;
- `z_density_profiles.csv`: component-resolved local-height densities;
- `lateral_rdf.csv`: same-component and cross-component two-dimensional RDFs.

## Batch analysis and separate workbook

Submit the single-CPU batch job from the repository root:

```bash
sbatch scripts/analyze_interface_holds.sbatch
```

It analyzes every available compressed `hold-300K/400K` and decompressed
`relax-300K/400K` trajectory, then writes:

```text
prepared/interface_structure_summary.xlsx
```

The workbook contains:

- `Results`: one compact row per system/stage;
- `Components`: full component-resolved anchoring, persistence, orientation,
  height, nearest-neighbor, and residence statistics;
- `Terminal States`: zero/one/two-anchor populations;
- `Hold vs Relax`: decompressed-minus-compressed changes at matched
  temperature;
- `Z Profiles`: local-height number-density data;
- `Lateral RDF`: normalized pair-distribution data;
- `Cutoff Sensitivity`: component bound fractions at the common 3.0, 3.25,
  and 3.5 Å contact definitions;
- `Methods`: definitions and interpretation limits.

Rebuild the workbook from completed result directories without rereading the
trajectories:

```bash
nio-md-prep summarize-interface prepared \
    --output prepared/interface_structure_summary.xlsx
```

Submission-time controls use environment variables:

```bash
sbatch --export=ALL,LAST_FRAMES=0,CONTACT_CUTOFF=3.25,RDF_RMAX=15.0 \
    scripts/analyze_interface_holds.sbatch
```

The batch search is recursive, so trajectories under
`prepared/replicates/seed-*/` enter the same workbook with their Packmol and
LAMMPS velocity seeds recorded.

## Independent deposition seeds

The original completed calculations are seed 01. The existing builders accept
`PREPARED_ROOT`, `PACKMOL_SEED`, and `VELOCITY_SEED`, so no parallel set of
replicate scripts is required. For seed 02, submit each stage only after
checking the previous stage:

```bash
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02,PACKMOL_SEED=202607242,VELOCITY_SEED=31415927 scripts/build_real_systems.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_deposition_array.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_hold_array.sbatch

sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02,PACKMOL_SEED=202608242,VELOCITY_SEED=31416927 scripts/build_sequential_systems.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_sequential_deposition_array.sbatch
sbatch --export=ALL,PREPARED_ROOT=prepared/replicates/seed-02 scripts/run_sequential_hold_array.sbatch
```

For seed 03, replace `seed-02` with `seed-03`, use primary seeds
`202607243`/`27182818`, and use sequential-stage seeds
`202608243`/`27183818`. Each sequential replicate automatically consumes the
independently deposited and held Me-4PACz substrate under the same
`PREPARED_ROOT`. The existing scripts still refuse to overwrite completed
outputs.

## Interpretation boundary

The module detects geometric contact and organization under the supplied
fixed-topology force field. It does not prove proton transfer, formation of a
new covalent/coordinate Ni-O-P bond, adsorption free energy, or electronic
coupling. Those claims require targeted electronic-structure or reactive
calculations.
