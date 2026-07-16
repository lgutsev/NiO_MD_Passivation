## LAMMPS simulations of a  passivated NiO corrugated surface
This code is present as-is with the goal of replicating the data presented in [doi/10.1002/aenm.202405367](https://onlinelibrary.wiley.com/doi/10.1002/aenm.202405367).  It can be modified to change the surface or the passivants used to suit the user’s task. 


For this project, the goal of the LAMMPs portion of the paper was to create a corrugated NiO-surface and passivate it with either just Me-4PACz or copassivate it with Me-4PACz and MPTMS-OH. After  applying the passivants to the surface via a moving wall potential we  ran it as a LAMMPS production run at 300K and 400K .  The features of the LAMMPs calculation are discussed more thoroughly in the draft. We will focus on file preparation here.

## Overall Workflow
The overall work process is shown in a diagram below:
![Plan](images/workflow.jpg)

## Corrugated Surface Preparation

In this module we will be focused on preparing the files to prepare the calculation. Some notes and details: 

The corrugated surface, with (011) & (011) cuts,  was prepared with Atomic Simulation Environment (ASE) tools. The passivant was added above the corrugated surface via PackMol. We didn’t attempt to use PackMol to pack the passivant tightly since there is experimental data indicating that  Me4PAcZ agglomerates. As such, we used a large vacuum above the surface and placed passivant relatively loosely. Afterwards, it was collected and applied to the surface during the LAMMPs simulation via a 12-6 LJ moving wall potential which ultimately replicated the agglomerating behavior seen in experiment.

## Forcefield Parameter Preparation
The OPLS-AA parameters and charges of the passivants  were generated with the LigParGen server:  <https://traken.chem.yale.edu/ligpargen/>  we want the .lmp format which we use as ligand.lmp

Now about the preparation of the files:

force\_field\_settings.lmp : which defines the pair, bond, angle, dihedral, improper force field coefficients.  The OPLS-AA parameters are generated as described above, note that they are LJ-style. For NiO we used OPENKim Buckingham parameters; however, the geometric mixing rule for generating mixed pair coefficients requires a 3-parameter forcefield (LJ). Thus, we had to generate our own NiO LJ parameters by fitting to experimental data & DFT calculations (110 slab, 8 layers deep). So the mixed parameters were LJ while the NiO parameters were Buckingham. Modifications to the OPLS-AA parameters also needed to be further modified in our case, the default phosphanate groups are too adhesive towards one another.  

### Header Preparation
We add this header manually:
\# Include the force field settings here

bond\_style      harmonic

angle\_style     harmonic

dihedral\_style  opls  

improper\_style  harmonic

pair\_style      lj/cut/coul/long 10.0 8.0 

kspace\_style    pppm 1e-4

kspace\_modify   slab 3.0

This is done because we need to manually modify the phosphanate groups anyways. 

## Topology and input file. 

topology\_output.lmp given the corrugated NiO.xyz file and the generated .lmp parameter file this is straightforward to generate. The masses list has to be modified afterwards since the surface atoms are not included, but they are included in the atoms list so it’s a small modification.

lammps.in prepared via the typical manner described in the manual of the program.

## Libraries Used

Atomic Simulation Environment (ASE) , along with its dependencies, and pandas.

## Conclusion

If there’s any problems or errors please report them to me @ <lgutsev@outlook.com>

Hope this code is helpful. 

## Reproducible preparation package

The notebook workflow is retained as a legacy reference. New builds use the offline, configuration-driven `nio-md-prep` package. Install it from this checkout with `python -m pip install -e .`. LigParGen and Codex are not expected to run on LONI; prepare and validate locally, then copy the completed build to LONI manually.

### Add a new molecule

1. Run LigParGen independently and download its LAMMPS `.lmp` output. This package does not contact LigParGen, invoke BOSS, or invent parameters.
2. Run `nio-md-prep init-molecule NAME`.
3. Put the downloaded file at the exact displayed path, named `ligpargen.lmp`.
4. Review the generated `molecule.toml`, including expected charge, role, anchors, and any reviewed override.
5. Run `nio-md-prep inspect-molecule NAME`.
6. Reference the slug in a study TOML and supply every molecule count explicitly.
7. Run `nio-md-prep build studies/your-study.toml --output builds/your-study`, then `nio-md-prep validate builds/your-study`.
8. Copy the completed directory to LONI manually.

The experimental 1:1 solution-mixing ratio does not by itself establish equal molecular counts in a simulation. The supplied study templates therefore contain no invented production counts. DCZ-4P is recorded as a two-phosphonic-acid-anchor molecule.

For CoSAM, both species are packed together and undergo one moving-wall protocol. Sequential preparation is genuinely two-stage: first produce and equilibrate the primary build, then run `nio-md-prep prepare-sequential-stage2 STUDY_CONFIG --primary-final primary_final.data --output stage2_directory`. The equilibrated stage-1 data is used as the existing component; the two species are never approximated as a common initial packing.

If Packmol is on `PATH`, a build runs it. Otherwise `packmol.inp` is retained and validation reports the exact offline command `packmol < packmol.inp`. The corrected examples are authoritative: hybrid LJ/Buckingham, geometric ligand/surface mixing, PPPM `1e-4`, slab correction `3.0`, `special_bonds amber`, `units real`, `atom_style full`, and `boundary p p f` are preserved.

## Why Excel is no longer required

The historical spreadsheet step pasted Packmol coordinates into replicated LigParGen atom records, shifted atom, molecule, and topology IDs, appended the shifted NiO atoms, and manually corrected LAMMPS header counts and box bounds. That was a manual data join, not part of the scientific method.

`nio-md-prep` now performs that join programmatically. Packmol supplies only coordinates and element labels; LigParGen supplies atom order, charges, types, masses, coefficients, and bonded topology; the packaged NiO file supplies the unchanged surface. The assembler checks the element sequence of every packed molecule, gives every physical ligand a unique positive molecule ID, remaps all topology and type namespaces from one authoritative mapping, assigns surface molecule ID 0, and calculates counts, bounds, and charge from the assembled records. Excel and manual pasting are no longer supported or required.

The component order for a CoSAM is study-manifest species order, then copies in Packmol order, then NiO. Sequential stage 2 preserves the stage-1 records and appends only the secondary ligand; it does not repack the established primary layer.

The user workflow is:

1. Obtain each LigParGen `.lmp` file manually.
2. Place it at `inputs/molecules/<molecule-name>/ligpargen.lmp`.
3. Run `nio-md-prep inspect-molecule <molecule-name>`.
4. Build the selected study.
5. If Packmol was unavailable, run `packmol < packmol.inp`, then rerun the build with `--packed-xyz packed.xyz`.
6. Run `nio-md-prep validate <build-directory>`.
7. Transfer the completed calculation files to LONI manually.
8. Run LAMMPS on LONI; neither Codex nor LigParGen is installed or run there.

LigParGen remains a manual external step. The legacy notebooks, historical examples, and paper Figure S27 are retained only as scientific provenance.



