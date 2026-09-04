# Scientific model and project design

This page collects the scientific assumptions that previously occupied a large fraction of the root README.

## Scientific picture

The corrugated NiO surface was constructed from the NiO(110) model with Atomic Simulation Environment tools and is stored as an authoritative LAMMPS data file. Ligands are placed relatively loosely in the vacuum above it with Packmol. The workflow deliberately does not ask Packmol to create an already dense monolayer: Me-4PACz aggregates experimentally, and the simulation is intended to let the film form dynamically rather than assume its final packing.

After a restrained relaxation, a repulsive 12-6 Lennard-Jones upper wall moves toward the surface while the lower wall prevents atoms from leaving the bottom of the nonperiodic cell. The wall collects the initially dispersed ligands and creates the compressed film. The deposited state then becomes the common starting point for independent 300 and 400 K branches.

Three related questions are separated by construction:

- **Pure SAM:** what film does Me-4PACz form by itself?
- **CoSAM:** what happens when Me-4PACz and a secondary ligand compete for the surface during the same deposition?
- **Sequential deposition:** what happens when the secondary ligand is added only after an Me-4PACz film has already formed?

The high-dose pure Me-4PACz system controls for the possibility that an apparent coverage improvement is caused simply by supplying more molecules.

## Systems supplied with the repository

The base system contains 180 Me-4PACz molecules on the 21,060-atom NiO surface. Secondary counts are derived from the stated equal-stock-volume assumption using

```text
Nsecondary / Nprimary =
    (Csecondary / MWsecondary) / (Cprimary / MWprimary)
```

with 0.5 mg/mL Me-4PACz and 0.3 mg/mL secondary stock concentrations. Molecular weights are calculated from the imported LigParGen mass inventory, and the unrounded ratio and rounding decision are written to `ratio_report.json`.

| System | Ligand inventory | Total ligands | Final atoms | Purpose |
|---|---:|---:|---:|---|
| Me-4PACz | 180 Me-4PACz | 180 | 29,160 | Experimental baseline |
| Me-4PACz/MeO-2PACz | 180 + 107 | 287 | 33,547 | CoSAM or inventory-matched sequential pair |
| Me-4PACz/MeO-4PADBC | 180 + 77 | 257 | 33,703 | CoSAM or inventory-matched sequential pair |
| Me-4PACz/DCZ-4P | 180 + 59 | 239 | 33,644 | CoSAM or inventory-matched sequential pair |
| High-dose Me-4PACz | 261 Me-4PACz | 261 | 32,805 | Mean-molecule-count loading control |

The CoSAM and corresponding completed sequential system contain exactly the same molecular inventory. Their difference is deposition history, not total loading.

The 261-molecule control is the arithmetic mean of the three mixed totals (239, 257, and 287). It is a molecule-count control, not an exact anchor-count control, because DCZ-4P has two phosphonic-acid anchors.

An experimental 1:1 solution-mixing statement does not by itself establish equal molecular counts. Production counts therefore remain explicit in the study TOML files. A molar 1:1 alternative can be represented by setting both ligand counts to 180 and using

```toml
composition.basis = "molar_1_to_1"
```

## Force field

Passivant charges and bonded OPLS-AA parameters are generated independently with the [LigParGen server](https://traken.chem.yale.edu/ligpargen/) and stored in LAMMPS `.lmp` format. The package does not contact LigParGen, run BOSS, or invent missing parameters.

The classical model preserves the hybrid treatment developed for the original work:

- ligand-ligand nonbonded interactions use Lennard-Jones plus long-range Coulomb terms;
- NiO self-interactions use the OPENKIM Buckingham parameters from the reviewed surface model;
- ligand/NiO cross interactions use pseudo-LJ Ni/O values fitted against experiment and DFT calculations for an eight-layer NiO(110) slab, followed by geometric mixing;
- phosphonate LJ terms and selected torsions are corrected because the unmodified LigParGen phosphonate groups are too adhesive toward one another;
- electrostatics use PPPM with slab correction.

The generated force-field header is equivalent to:

```lammps
bond_style      harmonic
angle_style     harmonic
dihedral_style  hybrid opls charmm
improper_style  cvff
pair_style      hybrid lj/cut/coul/long 10.0 8.0 buck/coul/long 10.0 8.0
pair_modify     pair lj/cut/coul/long mix geometric
kspace_style    pppm 1e-4
kspace_modify   slab 3.0
special_bonds   amber
```

Force-field coefficients remain in `force_field_settings_lammps_with_header.lmp`. Generated state files use `write_data ... nocoeff`, preventing stale coefficients from being copied into intermediate structures.

## What the classical model can and cannot support

The fixed-charge force field is primarily a structural and packing model. Coverage, void topology, roughness, adsorption height, broad orientation, and qualitative anchoring are the safest observables.

Do not treat the classical trajectories as direct predictions of chemical reaction, charge transfer, bond replacement, electronic polarization, or work-function change. Kinetic quantities such as residence times and surface-site exchange are retained for auditability but are deliberately withheld from publication-facing outputs under the default `classical-ff` analysis profile.

A standard energy/force MLIP can improve the potential-energy surface and make exchange kinetics more meaningful only after validation against relevant DFT configurations; it still does not automatically provide electronic polarization or work-function changes. See [`interfacial-analysis.md`](interfacial-analysis.md) for the profile policy.

## Historical provenance

The original workflow used notebooks and spreadsheet assembly to combine Packmol coordinates, LigParGen records, the NiO surface, and manually repaired LAMMPS headers. The current package performs this data-joining step programmatically and records provenance in manifests and hashes. Historical notebooks and examples remain in the repository as a record of how the published calculations were constructed, not as the recommended production path.
