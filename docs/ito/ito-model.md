# Proposed ITO substrate model for classical SAM simulations

Status: design as implemented on branch `feat/ito-sam-extension`. Numbers and sources are in [`parameter-sources.md`](parameter-sources.md). The literature basis is in [`literature-matrix.md`](literature-matrix.md). DFT cross-checks are prepared in [`vasp-smoke.md`](vasp-smoke.md).

## 1. What the substrate has to represent

In a nonreactive, fixed-charge packing model the substrate sets three things:
1. the lateral lattice and site pattern the anchors see;
2. the electrostatic field above the surface: termination, hydroxyls, charges;
3. the dispersion well that holds physisorbed molecules.

Everything chemical is outside the model: condensation to P–O–In, deprotonation, charge transfer and work-function change.

## 2. Proxy In2O3 versus explicit Sn-doped ITO

| Aspect | In2O3(111) proxy (`in2o3-111-oh`) | Explicit ITO (`ito-111-oh`) |
|---|---|---|
| Composition | In1536O2256 + 48 OH pairs | 142 Sn (Sn/(In+Sn) = 0.092) + 71 O_i |
| Charge neutrality | exact | exact: each 2 Sn_In (+0.525 each) is compensated by one O_i (−1.05) |
| Physical meaning | Sputtered ITO is bixbyite with Sn substituting In, so the lattice and surface sites are In2O3-like | Adds the ionic defect-cluster picture (Frank–Köstlin); **cannot** include the free electrons that make ITO conductive |
| Surface Sn | none | 5 of 192 top In5c become Sn via clusters just below the surface; no segregation model |
| Use | pilot default | sensitivity: the single-molecule scan compares it with the proxy under identical settings |

Recommendation: use the proxy unless the single-molecule scan (and later DFT) shows Sn changes the adsorption statistics beyond placement-to-placement scatter.

**Reason:** in a fixed-charge model, explicit Sn changes the surface field only through the compensating O_i. That field is an artefact of the missing electron gas, not a physical Sn effect.

## 3. Termination

The bixbyite (111) direction stacks neutral, composition-symmetric O–In–O trilayers. Cutting between them gives a stoichiometric slab with identical faces and zero dipole: a Tasker type-II surface, checked numerically to about 1e-16 e/Å. This is the lowest-energy In2O3 facet, and ITO films are typically (222)-textured. The top face carries 12 In5c and 12 O3c per primitive cell.

Not modelled:
- (001) and (110) facets
- the reconstructed (001) surfaces
- step edges and polycrystalline grain boundaries

Real UV-ozone- or plasma-cleaned ITO has a few-nm RMS roughness, which exceeds anything a 57 × 50 Å slab can represent.

## 3b. Corrugated surface (ITO analogue of the NiO groove)

The NiO campaign uses a corrugated NiO(110) slab (`inputs/surfaces/corrugated-nio-110`):
- 125.1 × 41.7 Å box
- one symmetric V-groove per cell, running along x
- 45° walls built from seven monatomic 2.085 Å steps, 14.6 Å deep in total
- opening about 27 Å, plateau about 14.6 Å
- 25 Å of NiO under the floor

The working hypothesis is that SAM molecules accumulate in such geometries. `in2o3-111-oh-groove` reproduces the geometry on In2O3(111):

| | corrugated NiO(110) | `in2o3-111-oh-groove` |
|---|---|---|
| Box (x along groove) | 125.1 × 41.7 Å (52.2 nm²) | 123.9 × 42.9 Å (53.2 nm²) |
| Step | 2.085 Å monatomic | 2.92 Å neutral O–In–O trilayer |
| Depth | 7 steps = 14.60 Å | 5 steps = 14.60 Å |
| Wall angle | 45° | 45° (step run = step height) |
| Opening / plateau | ~27 / ~14.6 Å | ~29 / ~13.7 Å |
| Under the floor | 25.0 Å | 5 trilayers = 14.6 Å (rigid slab, so thickness only affects the electrostatic and dispersion tail) |
| Surface chemistry | bare NiO, formal ±2 | 3.38 OH/nm² per projected area, as on the flat pilot slab |

**Construction** (`substrate.carve_groove`, then `rotate_quarter`): the top trilayers are removed on a V profile, with round((depth − |u − u0|)/d) trilayers removed at profile coordinate u. Lateral cuts leave a rim, which is cleaned up in two passes:
- dangling atoms (O with CN ≤ 1, cations with CN ≤ 2) are removed;
- exact neutrality is restored by removing the lowest-coordinated exposed rim atoms.

All 25 removals are recorded in the manifest. Hydroxylation then treats groove floors and walls as exposed, including fourfold step-edge In and twofold step-edge O.

**Limits:**
- The steps are ideal trilayer staircases, not relaxed step structures.
- Real ITO roughness comes from grains and facets of other orientations.
- This groove is a controlled geometric test of accumulation, as the NiO one is. It is not a model of a specific ITO morphology.

## 3c. Hydroxylation series on the corrugated slabs (ITO and NiO)

Both corrugated substrates come at five coverages. The levels follow the InterfaceForge NiO-MLIP convention:
- coverage is the fraction of the exposed-cation inventory, meaning exposed cations with CN below bulk;
- the motif is dissociated water (OH on the cation plus H on a nearby lattice O);
- the arrangement is scattered, by seeded farthest-point selection.

Placements that would clash are skipped, and each manifest records the shortfall.

| Requested | ITO groove: `in2o3-111-groove-ohXXX` | NiO groove: `corrugated-nio-110-rigid-ohXXX` |
|---|---|---|
| 0% | 0 OH/nm² | 0 OH/nm² |
| 25% | 25% achieved, 3.8 OH/nm² | 25% achieved, 5.8 OH/nm² |
| 50% | 48%, 7.3 OH/nm² | 49%, 11.3 OH/nm² |
| 75% | 64%, 9.8 OH/nm², **saturated** | 70%, 16.1 OH/nm² |
| 100% | 64%, identical geometry to 75% | 90%, 20.7 OH/nm² |

Densities are OH groups (terminal plus protonated lattice O) per projected nm². Figures: `data/ito_hydroxylation_series.png` and `data/nio_hydroxylation_series.png`, drawn with `scripts/ito/render_surfaces.py`.

**What experiments suggest:**
- Cleaved NiO(100) single crystals show no OH on regular terrace sites; OH forms only at defects and steps.
- Polar NiO(111) dissociates water but saturates at about 0.25 ML (Zhao et al., ACS Catal. 2016, 10.1021/acscatal.6b01997).
- Device NiOx does carry OH/NiOOH in XPS, more on ALD or solution-processed films than on sputtered ones.
- In2O3(111) saturates at 3 dissociated H2O per 1×1 cell under UHV, i.e. 3.38 OH/nm² on the flat slab. That is close to the 25% level here.

Realistic surfaces are therefore expected at the **low end: 0–25% for NiO, about 25% for ITO**. The upper levels bracket heavily hydroxylated, plasma- or UV-ozone-treated surfaces.

**NiO variants.** The coordinates are the authoritative campaign slab (`inputs/surfaces/corrugated-nio-110`), which is not modified. The variants use formal ±2 charges, CLAYFF-style hydroxyls, and the campaign's Ni/O pseudo-LJ values on a rigid slab (`cao-nio-rigid`). They are **not** the flexible Buckingham production model. The 0% variant is the rigid-slab control needed to separate the rigid-slab effect from the hydroxylation effect.

## 4. Hydroxylation

Cleaned ITO is hydroxylated. Water dissociates on In2O3(111) into terminal In–OH and surface O–H. The model places dissociated water pairs geometrically:
- each OH sits along the missing octahedral bond of a top In5c;
- each H sits on the nearest O3c;
- sites are chosen by seeded farthest-point sampling.

Default density: 3 pairs per primitive cell (3.38 OH/nm²); see the literature matrix for the basis. Placement is **unrelaxed**, and the smallest new-atom distance is 2.07 Å (an H-bond-like H···O). A DFT relaxation of the hydroxylated slab is the first VASP smoke case.

The bare slab (`in2o3-111-bare`) is kept as an extreme reference with exposed In5c, not as a realistic surface.

## 5. Oxygen vacancies

Vacancies are not generated. A neutral V_O in a fixed-charge model needs either:
- a compensating cation-charge reduction, which is an arbitrary choice of which In become "In+", or
- a net-charged slab, which is unphysical under PPPM.

The electronic question (vacancy-mediated binding, band bending) needs DFT or a charge-aware MLIP.

## 6. Charges

Charges are formal charges × 0.525 plus H(O) = +0.425. This is the CLAYFF octahedral-oxide pattern: lattice O −1.05, hydroxyl O −0.95, H +0.425. It is the only scaling under which dissociative hydroxylation and (2 Sn + O_i) doping stay exactly neutral with CLAYFF hydroxyls.

Why not formal charges: the NiO workflow uses ±2 because its Buckingham potential was parameterised with them. With a rigid slab and formal +3/−2 charges, the field above the surface is strongly exaggerated relative to DFT-derived partial charges (In typically ~+1.5–1.9 in Bader-type analyses; to be checked on the VASP slabs).

## 7. Slab dynamics

The slab is **rigid**: slab–slab pairs are excluded, slab forces are zeroed and slab velocities are zero. A flexible In2O3 needs a validated interatomic potential (shell-model Buckingham sets exist in the literature; see matrix §C). That would mean the CORESHELL machinery and a separate validation campaign.

For a packing and coverage model at 300–400 K on a crystalline oxide, rigidity costs little. It does remove surface-hydroxyl librations and any substrate response to anchoring. Consequences:
- the lateral box is the crystal cell;
- the ensemble is NVT on the ligands (the NiO NPT(xy) protocol would rescale the lattice).

## 8. SAM–surface interaction

- The ligands keep their NiO-workflow OPLS-AA/LigParGen terms with the Cao phosphonate corrections.
- The phosphonic acids stay **protonated and neutral**.
- SAM–slab cross terms are geometric means of ligand and slab self terms, with two bracketing slab sets:
  - UFF In/Sn with a strong dispersion well;
  - CLAYFF-analogue bare cations with essentially no cation dispersion.
- Anchoring can therefore only appear as:
  - H-bonds: P–OH···O(slab) and P=O···H–O(slab);
  - Coulomb/dispersion proximity of phosphonate O to In/Sn.

The analysis reports these as contacts (default 3.25 Å O···cation, 2.5 Å H-bond), never as bonds.

## 9. What the nonreactive model cannot establish

- P–O–In bond formation, condensation with surface OH, or water release
- deprotonation state, and mono-, bi- or tridentate binding mode, even though experiments (RAIRS/XPS on ITO) indicate deprotonated, surface-bound phosphonates
- binding energies comparable to experiment or DFT: the classical E_ads is a physisorption number
- charge transfer, polarisation, interface dipoles or work-function shifts
- anchoring kinetics, exchange, or desorption barriers (withheld by `model_scope`)

What it can support, as in the NiO work: coverage, void topology, clustering and aggregation, broad tilt distributions, the adsorption-height distribution, and competition between the two components of a mixture for surface area.

## 10. Reference: the 2PACz/PyCA-3F classical-MD model, and how ours differs

Li et al. (*Nat. Commun.* 2024, 10.1038/s41467-024-51760-5; extract in `literature-matrix.md` §B3) used:
- LAMMPS
- OPLS-AA/LigParGen ligands with Meltzer phosphonic-group terms
- a **fixed In2O3 slab** with a Buckingham potential, with no Sn and no hydroxyls; facet, box size, molecule counts and parameter source are not stated
- nonbonded-only molecule–surface coupling (UFF LJ + Coulomb, no P–O–In bonds)
- NVT at 300 K for 1 ns

They cite Park et al. (*Nature* 624, 289 (2023)) as the recipe, but that recipe was run on SnO2(110), not In2O3. A better-documented In2O3(111) SAM-MD recipe uses the Walsh et al. 2009 In2O3 potential (10.1021/cm902280z) refit to LJ with geometric mixing (10.1038/s41467-025-58111-y, matrix B4-3).

| Choice | Li et al. 2024 | This model | Reason |
|---|---|---|---|
| Slab | fixed In2O3, Buckingham | fixed In2O3(111) from COD 2310009, explicit trilayer termination | same rigidity; our termination and provenance are explicit |
| Hydroxyls | none | 3.38 OH/nm² (3 dissociated H2O per 1×1 cell, verified saturation value, matrix §C) | cleaned ITO is hydroxylated; bare slab kept as a reference |
| Sn | none | optional explicit (2 Sn_In + O_i), used as a sensitivity arm | shows whether Sn matters within the model |
| Surface charges | not stated | CLAYFF-scaled formal (0.525) | neutral under hydroxylation and doping |
| Cross terms | UFF + Coulomb | two brackets: UFF-cation / CLAYFF-analogue cation, CLAYFF O | the unknown In dispersion is bracketed rather than assumed |
| Anchor | protonated, nonbonded | protonated, nonbonded | identical limitation (§9) |
| Length | 1 ns NVT | 0.9 ns pilot (deposition → hold → release → relax), 3 seeds | seed replicas; longer runs after review |

**What their nonreactive force field could not establish, and ours cannot either:** binding mode, deprotonation, P–O–In bonds, interface dipole or work function. Their structural claims about mixing and orientation are of the same kind this model supports. Their chemical-anchoring statements rest on experiment and DFT, not on the MD.

## 11. Upgrade path

1. DFT (VASP smokes, `vasp-smoke.md`): relax the hydroxylated slab. Compare the classical E_ads against DFT for the physisorbed pose, and quantify the chemisorption (bidentate) energy the force field cannot represent.
2. Flexible slab: Walsh 2009 In2O3 potential (matrix C7) with the NiO-style `buck/coul/long` hybrid. This needs its own validation and is not required for coverage statistics.
3. MLIP route: the same NiO MLIP machinery (InterfaceForge VASP → MACE) applies once DFT trees exist.
