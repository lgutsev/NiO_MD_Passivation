# Phosphonate agglomeration structures for MLIP labeling

`nio-md-prep prepare-agglomeration` generates finite molecular clusters for
VASP energy/force labeling. It reconstructs each molecule from the reviewed
LigParGen file already stored under `inputs/molecules/`; separate XYZ files are
not required.

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
candidate is relaxed with GFN2-xTB before a POSCAR is written. xTB provides the
attractive intermolecular relaxation that Packmol cannot, while Packmol still
provides the different starting topologies.

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

The reference directory must contain `INCAR`, `KPOINTS`, and `POTCAR`. Every
other reusable top-level file, including Slurm launchers, is also copied. The
reference `POSCAR` and existing VASP outputs are excluded because the generated
agglomerate replaces the reference structure. Completed calculations are under
`vasp_runs/`, one fully independent run directory per structure. With xTB
enabled the first invocation normally stops with `XTB_REQUIRED`:

```bash
cd agglomeration/me-4pacz
sbatch run_xtb_array.sbatch
# after every array task has produced xtbopt.xyz:
cd -
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz.toml \
  --reference-dir /path/to/working/vasp-reference \
  --output agglomeration/me-4pacz
```

The rerun validates the xTB atom order and rejects an optimized structure with
an intermolecular O--O distance below 2.2 A or any molecule whose nearest
molecular neighbor is farther than 6.0 A. It then creates three VASP stages:
a 300 K hold, a 300-to-400 K heating ramp, and a 400 K hold. The latter is
submitted with an `afterok` dependency and copies the heating `CONTCAR` to its
own `POSCAR`. Run `submit_temperature_jobs.sh` inside a completed case to submit
all three jobs.

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
sbatch agglomeration/me-4pacz/run_xtb_array.sbatch
```

For the older LONI setup used by this project, no exports are necessary when
`$HOME/miniconda3/etc/profile.d/conda.sh` exists and the environment is named
`xtb`: the generated script activates it automatically. Nonstandard locations
can be selected without editing the generated script:

```bash
export CONDA_SH=/path/to/miniconda3/etc/profile.d/conda.sh
export XTB_ENV=xtb
sbatch agglomeration/me-4pacz/run_xtb_array.sbatch
```

The array uses one task and eight OpenMP cores per aggregate, sets
`ulimit -s unlimited`, and writes `xtb.out`, `xtb.err`, and `xtbopt.xyz` in each
`xtb/...` directory. Test one array index first with
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
for a one-replica offline or externally packed build.

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
