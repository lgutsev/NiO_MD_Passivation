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

The default contact cutoff is `auto`. The analyzer builds the distribution of
the nearest exposed-Ni distance for every phosphonate O atom and chooses the
first smoothed minimum after the contact peak. Both the distribution and the
selected cutoff are saved. A fixed sensitivity-test cutoff can be requested:

```bash
nio-md-prep analyze-interface prepared/me-4pacz-then-dcz-4p \
    --trajectory hold-300K.lammpstrj \
    --contact-cutoff 3.75
```

Important controls include:

- `--last-frames 0`: analyze every trajectory frame.
- `--stride N`: analyze every Nth selected frame.
- `--blocks N`: block count for within-trajectory SEM.
- `--persistence-threshold 0.50`: fraction of frames required for a molecule
  to be classified as persistently bound.
- `--surface-coordination-cutoff 2.8`: Ni-O neighbor cutoff used only to
  identify under-coordinated upper-surface Ni sites.
- `--rdf-rmax 15.0` and `--rdf-bin-width 0.25`: lateral RDF range and binning.
- `--z-min`, `--z-max`, and `--z-bin-width`: local-height density grid.

`rdf-rmax` must not exceed half the shortest lateral box length.

## Definitions

### Exposed Ni site

The static topology is used to select molecule-0 Ni atoms in the upper half
of the slab with fewer than six molecule-0 O neighbors within the surface
coordination cutoff. This coordination-based definition follows the
corrugated upper surface rather than selecting only atoms near one global
maximum-z plane.

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
- `Methods`: definitions and interpretation limits.

Rebuild the workbook from completed result directories without rereading the
trajectories:

```bash
nio-md-prep summarize-interface prepared \
    --output prepared/interface_structure_summary.xlsx
```

Submission-time controls use environment variables:

```bash
sbatch --export=ALL,LAST_FRAMES=0,CONTACT_CUTOFF=auto,RDF_RMAX=15.0 \
    scripts/analyze_interface_holds.sbatch
```

## Interpretation boundary

The module detects geometric contact and organization under the supplied
fixed-topology force field. It does not prove proton transfer, formation of a
new covalent/coordinate Ni-O-P bond, adsorption free energy, or electronic
coupling. Those claims require targeted electronic-structure or reactive
calculations.
