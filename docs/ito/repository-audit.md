# Repository audit: reusing the NiO classical-MD workflow for ITO

Audited at `origin/main` 4c82afa (2026-09-26). Paths are relative to `src/nio_md_prep/` unless stated.

## Decision

The ITO extension is a separate subpackage, `nio_md_prep.ito`. It reuses the generic building blocks directly and does **not** route through `build.build()` or `validate.validate()`. Those two functions encode the published NiO model: fixed constants, a fixed surface hash and a fixed atom count. Changing them would put the reproducibility of the NiO paper outputs at risk.

Shared modules received additive changes only: In, Sn and I masses in `geometry.ELEMENTS`, and In, Sn and I masses and Bondi radii in `analysis/coverage.py`. No NiO code path changes behaviour.

## Reused as-is

| Stage | Function(s) | Used by ITO |
|---|---|---|
| Molecule preparation | `config.molecule_manifest`, `lammps.parse`, `geometry.write_xyz`, `chemistry.molecular_weight` | `ito.assemble`, `ito.checks` |
| Phosphonate FF corrections | `chemistry.correction_lines`, `phosphonate_roles` (Cao LJ and torsion fixes) | same, unchanged |
| Replication / type remapping | `lammps.replicate`, `lammps.write`, `build._coeff_lines`, `build._surface` | same (pair lines de-hybridised by `_strip_pair_style`) |
| Packing | `build._packmol`, `_fit_packed_molecules` (seeded, atom-order checked, rigid ≤0.1 Å boundary fixes) | same |
| Coverage | `analysis.coverage.analyze_coverage` (projected vdW disks on a 0.2 Å periodic grid; near-surface gate against the substrate (mol ≤ 0) height map; void patches; roughness; block statistics) | called unchanged by `ito.analysis` |
| Geometry helpers | `interfacial._plane_orientation`, `_periodic_tree`, `_nearest_neighbor_mean`, `SiteExchangeTracker` | available; `ito.analysis` currently has its own simpler tilt and cluster code |
| Reporting policy | `analysis.model_scope.model_assessment("classical-ff")` | ITO outputs carry the same scope statement |
| Coverage-guided placement | `lego.identify_x_tunnel`, `lego2.identify_2d_void` (work from `coverage_probability.npz`) | reusable for sequential ITO studies later |

## NiO-specific assumptions that must be replaced

| # | Location | Assumption | Why it breaks for ITO | ITO replacement |
|---|---|---|---|---|
| 1 | `build.py:12-22` | Hard-coded `HEADER`: hybrid `lj/cut/coul/long` + `buck/coul/long` | ITO has no validated Buckingham set in the pipeline | `ito.forcefield.HEADER`: single `lj/cut/coul/long`, geometric mixing, `diff ad` |
| 2 | `build.py:151-169` | Surface types identified by Ni/O mass | In, Sn, OH-O and H are not Ni/O | Types come from `surface_manifest.json` `type_ids` |
| 3 | `build.py:400-414` | Ni–O Buckingham and Cao pseudo-LJ Ni (0.1, 3.0) / O (0.21, 3.05) cross terms | Fitted to NiO(110) DFT and experiment | Named, source-tagged slab LJ sets (`ito.forcefield.PARAMETER_SETS`) |
| 4 | `build.py:381,483`; `validate.py:107-108` | Surface atoms = abs(q) == 2.0 | ITO charges are 1.575 / 2.1 / −1.05 / −0.95 / 0.425 | Surface = molecule ID 0 (`group slab molecule 0`) |
| 5 | `validate.py:109-114` | 21060 surface atoms, fixed SHA256, exactly two surface types | Different slab | Per-model manifest hash checked at assembly |
| 6 | `build.py:297` | Default packing region `2 2 45 123.1 39.7 145` (NiO box) | Different box | Region derived from slab bounds and top-atom height |
| 7 | `build.py:421-423`; `docs/project-design.md` | Me-4PACz is primary; 0.5 / 0.3 mg/mL stocks | ITO mixtures use other ratios and stocks | Explicit counts in `studies/ito/*.toml`, with the literature basis recorded |
| 8 | `build.py` stage inputs | `fix npt ... x y couple xy` for all atoms | Barostatting a rigid crystalline slab rescales its lattice | NVT on mobile atoms; slab frozen; lateral box = crystal |
| 9 | `chemistry.phosphonate_roles`; `interfacial.py:427-447`; `agglomeration.py:145` | Every molecule has P(=O)(OH)2 | PyCA-3F (carboxylic acid) has no P | Not needed for this pilot (Me-4PACz, MeO-2PACz); a carboxylic anchor path is required before PyCA-3F |
| 10 | `interfacial.py:255-375` | "Exposed Ni" sites: element Ni, CN < 6, hard-coded NiO reference path | No Ni | `ito.analysis`: In/Sn cation contacts plus H-bonds to slab O / OH |
| 11 | `interfacial.py:737`; `analyze_interface_*.sbatch` | Ni–O(P) contact cutoff 3.25 Å (sensitivity 3.0–3.5) | Calibrated for Ni | Same default, labelled a proximity criterion; needs recalibration against DFT (see `vasp-smoke.md`) |
| 12 | `geometry.py:5`; `coverage.py:28-54` | Mass and radius tables lack In, Sn, I | Element lookup raises | Fixed additively |
| 13 | `interfacial_report.py`, `publication_report.py`, `report.py` | NiO / exposed-Ni wording; Me-4PACz-primary labels; DCZ-4P special cases | Cosmetic, but misleading for ITO | ITO outputs use their own JSON summaries; the workbooks are not reused yet |
| 14 | `scripts/*.sbatch` | Hard-coded NiO study lists, array sizes, `me-4pacz-alone/held-300K.data` | Different studies | New `scripts/ito/` launchers (VASP smokes); the MD launcher is still to be written |
| 15 | Protocol | 0.5 fs, no SHAKE, 600k-step deposition on a 285 Å box | Cost | ITO pilot: SHAKE X–H at 1 fs, zhi 125 Å, PPPM `diff ad` (measured ~2x cheaper than `ik`) |

## Comparability with the NiO campaign

These settings are kept identical to the NiO campaign:
- LigParGen OPLS-AA ligand terms
- Cao phosphonate corrections
- `special_bonds amber`
- cutoffs 10 / 8 Å
- PPPM slab correction
- 300 K
- moving-wall deposition with a 30 Å endpoint clearance
- areal dose (3.45 molecules/nm²)

These settings differ, so an NiO-vs-ITO comparison must be read with them in mind:
- rigid vs flexible substrate
- NVT vs NPT(xy)
- SHAKE at 1 fs vs 0.5 fs
- a shorter protocol
