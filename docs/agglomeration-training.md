# Phosphonate agglomeration structures for MLIP labeling

`nio-md-prep prepare-agglomeration` generates finite molecular clusters for
VASP energy/force labeling. In a single-species workflow, the single molecule
in `--reference-dir/POSCAR` is the authoritative geometry. A mixed workflow
accepts one slug-qualified reference directory per species and uses each
directory's single-molecule POSCAR for that species. The reviewed LigParGen
files under `inputs/molecules/` supply atom identity, topology, masses, and the
canonical atom order used by Packmol and xTB. No separate XYZ file is required.

The generator is deliberately independent of the classical NiO/LAMMPS force
field. Its output is reference data for a later MLIP dataset, not a classical
MD starting state.

## Generated ensemble

For each replica, the command asks Packmol for an independently randomized
spherical cluster. Packmol only supplies collision-free placement; it does not
guarantee that every ligand is associated. The generator therefore translates
whole molecules inward with a bisection search until the nearest intermolecular
contact reaches `compact_to_distance_angstrom`. This adaptive operation leaves
all bonds and angles unchanged and cannot cross the configured contact floor.
It replaces the old fixed compression factor that could make O--O contacts too
short in one seed while leaving another seed too dispersed.

The supplied Me-4PACz campaign varies cluster size as well
as position and orientation: three replicas each contain 2, 3, or 4 molecules;
two replicas contain 6 molecules; and one contains 8. This concentrates the
sampling budget on inexpensive small clusters instead of spending equally on
large, costly systems whose trajectories already contain many local contacts.

Every family uses only `center_scale = 1.0`. After adaptive compaction, each
candidate follows a reproducible xTB preconditioning protocol before a POSCAR
is written: loose GFN2-xTB optimization, 20 ps total NVT dynamics at 450 K with
a 1 fs step, then a final loose GFN2-xTB quench. Replica 0 of each size family
uses a 6 ps head-facing segment, followed by 14 ps of completely unbiased
dynamics. All remaining replicas run 20 ps without restraints. For every
greedily paired head, the generator restrains P--P separation and the two
`C--P...P` angles, where C is the carbon covalently bonded to that P. Driving
both angles toward 170 degrees makes the two molecular `C->P` head axes face one
another; P--P distance alone cannot distinguish sideways or back-to-back
approaches. Each P head is used at most once and heads from the same molecule
are never paired. No oxygen atom is restrained.

The distance target is calculated after the initial optimization: it is never
more than 1.5 A inward from that optimized P--P distance and never pushes an
already-close pair apart. The finite 170-degree angular target avoids a
singular exactly linear coordinate. The MD seeds are derived from the Packmol
seed.

The steering is only a starting-configuration sampling aid. The restraints are
absent from the long segment, the final quench, and every VASP calculation. xTB
provides attractive intermolecular relaxation that Packmol cannot, while
Packmol still provides the different starting topologies. xTB energies and
forces are not MLIP labels; only the subsequent VASP calculations provide
training labels. xTB's distance restraints follow its documented `$constrain`
input syntax: <https://github.com/grimme-lab/xtb/blob/main/man/xcontrol.7.adoc>.

Optional scale families remain supported. The generator translates whole
molecules so their centers of mass are multiplied by every configured
`center_scale` about the cluster center. Intramolecular coordinates are never
scaled. Values below one create controlled compression, one preserves the
Packmol cluster, and values above one provide separated controls.

Every scale belonging to one replica uses the same cubic cell. Consequently,
energy differences within that replica do not contain a changing-cell artifact.
The generator rejects any compressed variant whose nearest intermolecular atom
pair is below `minimum_distance_angstrom`.

Every completed case is also required to pass an automatic structural sanity
check. The check confirms that coordinates and cell dimensions are finite, all
atoms lie inside the cell, the requested vacuum is present at every cell face,
the cluster has at least twice that vacuum to its nearest periodic image,
intermolecular distances remain above the configured safety threshold, and
Packmol plus the rigid center scaling preserved every intramolecular distance.
The written POSCAR header is checked against the alphabetical element order.

Normal operation requires a complete, working VASP reference directory and
creates launch-ready calculation packages:

The supplied pilot contains 12 independent loose-start calculations across
five cluster sizes. Its 8 A face vacuum is a cost-conscious starting point;
compare a small subset against 10 A before treating absolute cluster energies
as converged isolated-cluster values.

```bash
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz.toml \
  --reference-dir /path/to/working/vasp-reference \
  --output agglomeration/me-4pacz
```

The reference directory must contain `POSCAR`, `INCAR`, `KPOINTS`, and `POTCAR`.
Its `POSCAR` must contain exactly one molecule with the same elemental
composition as the configured LigParGen molecule. POSCAR element grouping is
handled automatically, and atoms may appear in a different order within each
element block. The command infers the molecular bond graph, maps it to the
LigParGen topology, restores the canonical order, and checks every topology bond
before packing. It fails closed if the graphs are incompatible. Every other
reusable top-level file, including Slurm launchers, is also copied. The
reference `POSCAR` is used as the molecular template but is not copied verbatim:
the generated agglomerate replaces it. Existing VASP outputs are also excluded.
Completed calculations are under
`vasp_runs/`, one fully independent run directory per structure. With xTB
enabled the first invocation normally stops with `XTB_REQUIRED`:

```bash
cd agglomeration/me-4pacz
./audit_xtb_runs.py
./launch_xtb.sh --dry-run
./launch_xtb.sh
# after every array task has produced xtbfinal.xyz:
cd -
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz.toml \
  --reference-dir /path/to/working/vasp-reference \
  --output agglomeration/me-4pacz
```

`audit_xtb_runs.py` is safe to run while jobs are active. It reports the
completed fraction and classifies every case as pending, initial optimization,
steering, unbiased MD, final-but-unmarked, hash mismatch, or safely complete.
It also refreshes `xtb_progress.csv` and `xtb_progress.json` for automated
monitoring. `launch_xtb.sh` runs that audit automatically, reports how many
cases remain, and asks for one campaign-wide `y/N` confirmation before calling
Slurm. Use `--yes` for intentional noninteractive submission. The underlying
launcher compacts the unfinished TSV row numbers into a Slurm selection such
as `--array=2,5,7-10%4`, so completed folders do not receive array jobs. Each
selected task still checks its completion marker again immediately before xTB
as protection against a calculation finishing between the audit and submission.

The rerun converts `xtbfinal.xyz` automatically; no separate XYZ-to-POSCAR
command is needed. It validates the xTB atom order and rejects a quenched
structure with
an intermolecular O--O distance below 2.2 A or any molecule whose nearest
molecular neighbor is farther than 6.0 A. It then creates three VASP stages:
a 300 K hold, a 300-to-400 K heating ramp, and a 400 K hold. The latter is
submitted with an `afterok` dependency and copies the heating `CONTCAR` to its
own `POSCAR`. Run `submit_temperature_jobs.sh` inside a completed case to submit
all three jobs. A completed campaign also contains a root launcher covering
every structure case. Preview every `sbatch` command without submitting:

```bash
agglomeration/me-4pacz/launch_vasp_runs.sh --dry-run
```

Launch interactively with a single campaign-wide `y/N` confirmation:

```bash
agglomeration/me-4pacz/launch_vasp_runs.sh
```

For an intentional noninteractive submission, use `--yes`. The 400 K hold is
still submitted with an `afterok` dependency on its corresponding heating job.

### Regenerating VASP folders from an existing completed xTB campaign

Older campaigns whose root manifest is already `COMPLETE` can be upgraded or
rebuilt in place without rerunning Packmol or xTB. Repeat the original command
with `--regenerate-vasp`:

```bash
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz.toml \
  --reference-dir /path/to/working/vasp-reference \
  --output agglomeration/me-4pacz \
  --regenerate-vasp
```

The command reuses the existing completed `xtbfinal.xyz` files, rebuilds the
temperature-stage VASP inputs, and creates or refreshes
`launch_vasp_runs.sh`. For a mixed campaign, repeat both slug-qualified
`--reference-dir` arguments and the original `--vasp-template-slug`, then add
the same `--regenerate-vasp` flag. This mode tolerates a changed TOML checksum
for older campaigns; Packmol geometry and xTB protocol fingerprints still
invalidate incompatible intermediates instead of silently reusing them.
In explicit regeneration mode, an existing `xtbfinal.xyz` accompanied by
`xtb_protocol.complete` is treated as the operator-confirmed completed xTB
result even when the completion marker uses an older hash schema. A final XYZ
without the completion marker remains pending. The final XYZ atom identities,
ordering, distances, molecular association, and VASP compatibility are still
validated before any VASP launcher is created.
Regeneration fails closed if it finds VASP runtime outputs such as `OUTCAR`,
`CONTCAR`, `vasprun.xml`, or Slurm output, preventing an already-started
calculation from being overwritten.

## Mixed agglomerates from two VASP folders

Use a configuration whose composition contains both slugs; the ready-made
example is `examples/agglomeration/me-4pacz-dcz-4p.toml`. Supply one qualified
reference per slug:

```bash
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz-dcz-4p.toml \
  --reference-dir me-4pacz=/path/to/Me4PACz \
  --reference-dir dcz-4p=/path/to/DCZ4P \
  --vasp-template-slug me-4pacz \
  --output agglomeration/me-4pacz-dcz-4p
```

Every qualified directory must contain `POSCAR`, `INCAR`, `KPOINTS`, and
`POTCAR`. Each POSCAR must contain exactly one intact molecule matching its
slug's LigParGen topology. `--vasp-template-slug` selects the directory whose
`INCAR`, `KPOINTS`, launcher, and other reusable top-level files are copied to
the mixed calculations. If it is omitted, the first qualified
`--reference-dir` on the command line is used.

The tool extracts the required elemental datasets from the supplied POTCARs
and writes a combined POTCAR in the generated POSCAR's alphabetical element
order. When the same element occurs in both source POTCARs, its dataset must be
byte-identical; conflicting pseudopotentials fail closed instead of silently
mixing VASP setups. The root manifest records both geometry mappings, both
source POTCAR hashes, the chosen calculation-template slug, and the source slug
used for every combined POTCAR dataset.

With xTB enabled, continue exactly as for a pure campaign:

```bash
cd agglomeration/me-4pacz-dcz-4p
./audit_xtb_runs.py
./launch_xtb.sh
cd -

nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz-dcz-4p.toml \
  --reference-dir me-4pacz=/path/to/Me4PACz \
  --reference-dir dcz-4p=/path/to/DCZ4P \
  --vasp-template-slug me-4pacz \
  --output agglomeration/me-4pacz-dcz-4p
```

The second invocation validates the two-species xTB atom order and creates the
same 300 K, heating, and 400 K VASP stages used by single-species campaigns.
The xTB protocol fingerprints the generated mixed-species `input.xyz`. If a
reference POSCAR or packed geometry changes during a resumable campaign, the
launcher archives the stale initial optimization and downstream xTB products
before recalculating them; an old `xtbfinal.xyz` cannot be accepted for the new
geometry.

## Installing xTB on LONI

xTB does not need to be a LONI module. Install it once in project or work
storage with Conda/Miniforge, following the [official xTB installation
instructions](https://xtb-docs.readthedocs.io/en/latest/setup.html):

```bash
# after installing Miniforge, if conda is not already available
conda create -y -p /path/to/project/envs/xtb -c conda-forge xtb
conda activate /path/to/project/envs/xtb
xtb --version
```

Before submitting the generated array, either activate that environment or
export the executable explicitly. Slurm exports the submission environment by
default:

```bash
export XTB_EXE=/path/to/project/envs/xtb/bin/xtb
cd agglomeration/me-4pacz
./launch_xtb.sh
```

No exports are necessary for the current LONI installation: the generated
script first tries `/project/lgutsev/env/xtb_env/bin/xtb`. Nonstandard locations
can be selected without editing the generated script:

```bash
export XTB_ENV=/path/to/another/xtb/environment
cd agglomeration/me-4pacz
./launch_xtb.sh
```

The generated array embeds the campaign root's absolute path, so Slurm's spool
copy of the script cannot redirect it away from `xtb_cases.tsv`. If the campaign
directory is moved, regenerate the launcher or update its `root=` line. The
array uses one task and one CPU per aggregate and sets `ulimit -s unlimited`.
Each `xtb/...` directory records the complete staged provenance:

```text
input.xyz                 # compacted Packmol start
steering_plan.json        # paired P/C indices and bounded steering parameters
md_steered.inp            # generated after initial opt for optional 6 ps steering
steering_restraints.json  # realized P--P distances, angles, and targets
md.inp                    # 14 or 20 ps unbiased 450 K segment
xtb_protocol.json         # exact timings, seeds, steering rule, and P pair map
xtb_protocol.sha256       # expected protocol fingerprint
xtb_protocol.complete     # written only after the matching final quench
xtbopt_initial.xyz        # initial loose GFN2-xTB minimum
xtb_steered.trj           # optional restrained trajectory
xtbsteered_last.xyz       # optional steering endpoint
xtb_unbiased.trj          # unbiased trajectory used for handoff
xtbmd_unbiased_last.xyz   # last complete unbiased frame
xtbfinal.xyz              # final loose quench; source used for VASP
xtb_initial_opt.out/.err
xtb_steered_md.out/.err
xtb_unbiased_md.out/.err
xtb_final_opt.out/.err
```

The campaign root also contains `write_xtb_steering_input.py`, which calculates
the bounded targets from `xtbopt_initial.xyz` immediately before the steered
segment begins.

Completed stages are skipped on resubmission, so an interrupted task resumes
from its last valid file. Results from a different protocol fingerprint are not
accepted. The launcher preserves their stage files with a
`.protocol_<old-hash>` suffix, retains the reusable initial optimization, and
reruns steering, unbiased MD, and the final quench. Test one array index first with
`sbatch --array=0 run_xtb_array.sbatch` before releasing the complete array.

Geometry-only generation remains available, but must be requested explicitly:

```bash
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz.toml \
  --structures-only \
  --output agglomeration/me-4pacz
```

Regular top-level reference files are copied into every case, except existing
VASP structures and calculation outputs such as `POSCAR`, `CONTCAR`, `OUTCAR`,
`WAVECAR`, and `CHGCAR`. This permits local `POTCAR`, `INCAR`, `KPOINTS`, and
launcher reuse without placing licensed or generated files in the repository.
When `POTCAR` is present, its `TITEL` sequence must exactly match the generated
POSCAR header; a mismatch stops generation before calculations can be launched.
For the supplied Me-4PACz example, the required order is `C H N O P`. In
`--structures-only` mode, POTCAR is reported as `NOT_SUPPLIED`; the generated
folders are not represented as launch-ready calculations.

If Packmol is unavailable, the command writes one `packmol.inp` directory per
replica and records `PACKMOL_REQUIRED` in the root manifest. Run
`packmol < packmol.inp` inside every incomplete replica directory and rerun the
same command with the same output path. The second invocation validates the
produced atom order and completes the structures. `--packed-xyz` is available
for a one-replica offline or externally packed build. A geometry fingerprint is
stored beside every Packmol input. If the reference POSCAR changes, stale
`packed.xyz` files are archived with a `geometry_*` suffix and must be packed
again; the corresponding xTB protocol is then invalidated as well. In
`--structures-only` mode, where no VASP reference exists, coordinates fall back
to the LigParGen file and this provenance is recorded in the root manifest.

## Configuration

```toml
[agglomeration]
packmol_tolerance_angstrom = 3.0
vacuum_angstrom = 8.0
minimum_distance_angstrom = 2.5
compact_to_distance_angstrom = 3.0
center_scales = [1.00]

[xtb]
enabled = true
gfn = 2
opt_level = "loose"
max_cycles = 100
md_enabled = true
md_temperature_K = 450
md_time_ps = 20.0
md_step_fs = 1.0
md_dump_fs = 50.0
md_sccacc = 1.0
head_bias_enabled = true
head_bias_replicas_per_family = 1
head_bias_time_ps = 6.0
head_bias_target_pp_distance_angstrom = 4.5
head_bias_max_distance_reduction_angstrom = 1.5
head_bias_target_axis_angle_degrees = 170.0
head_bias_force_constant = 0.005
minimum_oxygen_oxygen_distance_angstrom = 2.20
maximum_nearest_molecule_distance_angstrom = 6.0

[vasp_md]
heating_steps = 1000
hold_steps = 3000

[[agglomerates]]
name = "n02"
replicas = 3
base_seed = 2405202
radius_angstrom = 11.5
molecules = [{ slug = "me-4pacz", count = 2 }]

# n03 and n04 also use 3 replicas; n06 uses 2; n08 uses 1.
```

Each `[[agglomerates]]` table can define a different composition. Mixtures use
multiple inline molecule entries, for example:

```toml
[[agglomerates]]
name = "mixed-n04"
replicas = 2
base_seed = 2405900
radius_angstrom = 14.0
molecules = [
  { slug = "me-4pacz", count = 3 },
  { slug = "dcz-4p", count = 1 },
]
```

Counts are explicit. The tool does not infer experimental solution ratios.
For launch-ready mixed VASP cases, pair this composition with the repeated
slug-qualified `--reference-dir` form shown above. A single unqualified
`--reference-dir /path` remains the backward-compatible single-species form.

## Outputs and provenance

```text
agglomeration/me-4pacz/
  agglomeration_manifest.json
  agglomeration_index.csv
  agglomeration_sanity.csv
  templates/
  packmol/n02/replica_000/
  packmol/n03/replica_000/
  packmol/n04/replica_000/
  packmol/n06/replica_000/
  packmol/n08/replica_000/
  xtb/n02/r000_s00_1p000/
    input.xyz
    xtbopt.xyz
    xtb.out
  vasp_runs/n02/r000_s00_1p000/
    POSCAR
    structure.xyz
    agglomeration_manifest.json
    sanity_report.json
    sanity_report.txt
    submit_temperature_jobs.sh
    300K/
      POSCAR
      INCAR                    # TEBEG=TEEND=300
      KPOINTS
      POTCAR
      runvasp.sh
    400K/01_heat_300_to_400K/  # TEBEG=300, TEEND=400
    400K/02_hold_400K/         # starts from the heating CONTCAR
```

`POSCAR` uses an alphabetical element header and groups coordinates in the same
order. The per-case manifest maps every POSCAR index back to its molecule slug,
molecule ID, local atom index, and original packed-XYZ index. The root manifest
records molecule/configuration hashes, seeds, scales, cell dimensions, and
minimum intermolecular distances.

`sanity_report.json` contains the individual boolean checks and exact cell,
bounding-box, face-clearance, periodic-image, intermolecular-distance, and
intramolecular-rigidity values. It explicitly identifies the closest
intermolecular atom pair and the closest intermolecular O--O pair, including
their molecule IDs, so a phosphonate-head collision is directly auditable.
`agglomeration_sanity.csv` provides one compact row per generated calculation.
A case is written as scientifically usable only with `sanity_status = PASS`.

All structures derived from the same Packmol replica are correlated and must
remain in the same training, validation, or test split. A low random-frame
error obtained after separating those variants across splits would constitute
data leakage.

## Scientific scope

These clusters cover ligand-ligand association and controlled separation. They
do not replace NiO-supported aggregates, isolated conformer sampling, proton-
transfer configurations, or short-range repulsive examples. Those should enter
the final DFT dataset as separate configuration families.
