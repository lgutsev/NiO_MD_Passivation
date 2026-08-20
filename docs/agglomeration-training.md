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
spherical cluster. It then translates whole molecules so that their centers of
mass are multiplied by every configured `center_scale` about the cluster center.
Intramolecular coordinates are never scaled. Values below one create controlled
compression, one preserves the Packmol cluster, and values above one provide
partially and strongly separated controls.

Every scale belonging to one replica uses the same cubic cell. Consequently,
energy differences within that replica do not contain a changing-cell artifact.
The generator rejects any compressed variant whose nearest intermolecular atom
pair is below `minimum_distance_angstrom`.

The supplied pilot is:

```bash
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz.toml \
  --output agglomeration/me-4pacz
```

To make each generated directory immediately runnable, point to a functioning
VASP reference directory:

```bash
nio-md-prep prepare-agglomeration \
  examples/agglomeration/me-4pacz.toml \
  --reference-dir /path/to/working/vasp-reference \
  --output agglomeration/me-4pacz
```

Regular top-level reference files are copied into every case, except existing
VASP structures and calculation outputs such as `POSCAR`, `CONTCAR`, `OUTCAR`,
`WAVECAR`, and `CHGCAR`. This permits local `POTCAR`, `INCAR`, `KPOINTS`, and
launcher reuse without placing licensed or generated files in the repository.

If Packmol is unavailable, the command writes one `packmol.inp` directory per
replica and records `PACKMOL_REQUIRED` in the root manifest. Run
`packmol < packmol.inp` inside every incomplete replica directory and rerun the
same command with the same output path. The second invocation validates the
produced atom order and completes the structures. `--packed-xyz` is available
for a one-replica offline or externally packed build.

## Configuration

```toml
[agglomeration]
replicas = 6
base_seed = 2405367
radius_angstrom = 12.0
packmol_tolerance_angstrom = 2.0
vacuum_angstrom = 12.0
minimum_distance_angstrom = 1.10
center_scales = [0.90, 1.00, 1.30, 1.80]

[[molecules]]
slug = "me-4pacz"
count = 4
```

Mixtures use multiple molecule entries, for example:

```toml
[[molecules]]
slug = "me-4pacz"
count = 3

[[molecules]]
slug = "dcz-4p"
count = 1
```

Counts are explicit. The tool does not infer experimental solution ratios.

## Outputs and provenance

```text
agglomeration/me-4pacz/
  agglomeration_manifest.json
  agglomeration_index.csv
  templates/
  packmol/replica_000/
  structures/r000_s00_0p900/
    POSCAR
    structure.xyz
    agglomeration_manifest.json
```

`POSCAR` uses an alphabetical element header and groups coordinates in the same
order. The per-case manifest maps every POSCAR index back to its molecule slug,
molecule ID, local atom index, and original packed-XYZ index. The root manifest
records molecule/configuration hashes, seeds, scales, cell dimensions, and
minimum intermolecular distances.

All structures derived from the same Packmol replica are correlated and must
remain in the same training, validation, or test split. A low random-frame
error obtained after separating those variants across splits would constitute
data leakage.

## Scientific scope

These clusters cover ligand-ligand association and controlled separation. They
do not replace NiO-supported aggregates, isolated conformer sampling, proton-
transfer configurations, or short-range repulsive examples. Those should enter
the final DFT dataset as separate configuration families.
