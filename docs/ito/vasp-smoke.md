# ITO/SAM VASP reference inputs (smoke tests)

## Progress notes

- 2026-09-26: I read the ITO substrate, assemble, lammps and geometry modules, and the VASP conventions of
  InterfaceForge@4501e340. That covers the NiO(110)+phosphonate template INCAR, `profiles/potcar_pbe_54.yaml`,
  the `profiles/loni.yaml` VASP jobs, `launch_scripts/runvasp.sh`, and `vasp.py` (`resolve_potcar_root`,
  `assemble_potcar`).
- Implemented `src/nio_md_prep/ito/vasp.py`, `scripts/ito/prepare_vasp_smokes.py` and `tests/test_ito_vasp.py`.
- Regenerated `inputs/ito/vasp_smoke/` after the bixbyite coordinates were corrected to COD 2310009
  (In 24d x = 0.4663; O 48e 0.3912, 0.1558, 0.3796). The manifests record the new lattice block.
- 16/16 tests pass. The `run.sbatch` POTCAR, ZVAL and NELECT logic was dry-run locally against a fake PAW tree
  with `module` and `srun` stubbed. This included a negative ZVAL test, which exits with code 3. **No VASP job
  was run anywhere.**
- 2026-09-27 (coordinator): hydroxylation filters moved into `substrate.hydroxylate`; all cases regenerated
  (`--check` 74/74 identical); committed. Still nothing submitted.

## Purpose

The classical ITO/SAM model is nonreactive: a rigid slab, a fixed-charge slab with LJ mixing, and protonated
phosphonic acids. These cases give it a small set of DFT reference points. All of them use the orthorhombic 1x1
In2O3(111) cell (14.3076 x 24.7815 A, two O-In-O trilayers for smoke). They all share one common cell height
(30.0 A), so the same plane-wave basis and dipole treatment enter every energy difference.

| Case | Question it answers for the classical model |
|---|---|
| `slab-bare` | Bare (dry) In2O3(111) reference. Shows how far the rigid bulk-terminated slab is from a relaxed one: the first-steps forces and displacements of the free trilayer. |
| `slab-oh` | E_slab term of every E_ads. Also checks whether the geometric hydroxylation used by the FF slab is stable. |
| `slab-oh-sn` | Sn effect on the substrate: one neutral 2Sn_In + O_i cluster in the lower trilayer. Compare surface relaxation, charges and gap with `slab-oh`. It is also the E_slab reference for a later adsorption-on-Sn case, which is **not** generated yet. |
| `mol-<slug>` | E_mol term: the isolated protonated acid in the same cell. It is a rigid copy of the physisorbed start molecule, shifted only along z. |
| `ads-<slug>-phys` | Upright physisorbed acid, H-bonded to the surface OH/O. This is the configuration class the classical model represents: E_ads(phys), DFT vs classical. |
| `ads-<slug>-bidentate` | Bridging-bidentate P-O-In chemisorption, with both acidic H moved to surface O3c. It has exactly the same atoms as `phys`. E(bid) − E(phys) is the chemisorption energy the FF cannot capture. |

`<slug>` is `me-4pacz` or `meo-2pacz`.

## Generated cases (smoke, 2 trilayers)

| Case | Formula | Atoms | Fixed | NELECT | Min nonbonded (A) | Min mol–slab (A) | Vacuum (A) | Resources | Est. wall h |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| slab-bare | In64O96 | 160 | 80 | 1408 | 2.80 | – | 25.0 | workq × 64 | 0.3 |
| slab-oh | H12In64O102 | 178 | 80 | 1456 | 2.25 | – | 23.2 | workq × 64 | 0.3 |
| slab-oh-sn | H12In62O103Sn2 | 179 | 70 | 1464 | 2.22 | – | 23.2 | workq × 64 | 0.3 |
| mol-me-4pacz | C18H22NO3P | 45 | 0 | 122 | 1.74 | – | 18.7 | single × 16 | 0.05 |
| ads-me-4pacz-phys | C18H34In64NO105P | 223 | 80 | 1578 | 1.74 | 1.97 | 12.3 | workq × 64 | 0.4 |
| ads-me-4pacz-bidentate | C18H34In64NO105P | 223 | 80 | 1578 | 1.74 | 1.93 | 12.2 | workq × 64 | 0.4 |
| mol-meo-2pacz | C16H18NO5P | 41 | 0 | 122 | 1.75 | – | 19.9 | single × 16 | 0.05 |
| ads-meo-2pacz-phys | C16H30In64NO107P | 219 | 80 | 1578 | 1.75 | 1.98 | 13.6 | workq × 64 | 0.4 |
| ads-meo-2pacz-bidentate | C16H30In64NO107P | 219 | 80 | 1578 | 1.75 | 1.94 | 13.6 | workq × 64 | 0.4 |

Notes on the table:
- **Min nonbonded:** the minimum over all pairs that are not in the explicit bond list. The bond list covers
  LigParGen bonds, slab cation–anion pairs under 2.6 A, the constructed O–H bonds, and the bidentate In–O and
  transferred H–O bonds. The molecular minima (1.74–1.75 A) are geminal H···H of the LigParGen geometry.
- **Est. wall h:** an **estimate** from an uncalibrated scaling model, not a measurement. The model assumes
  about 30 s per electronic step for 1000 bands and 2×10⁵ Γ plane waves on 64 cores, and 65 electronic steps
  for NSW = 3.
- **Memory:** about 10 GB per slab case (estimate).
- **Measured values:** `inputs/ito/vasp_smoke/README.md` and each `case_manifest.json` carry the exact numbers.

**Electron counts** come from the ZVAL table in `POTCAR.spec`: In_d 13, Sn_d 14, O 6, P 5, N 5, C 4, H 1.
`run.sbatch` re-reads ZVAL from the assembled POTCAR and refuses to run if the valences or the NELECT total differ.

### Construction details and checks

- **Slab:**
  - `substrate.build_slab(1, 1, trilayers)` from the corrected COD 2310009 coordinates
  - In–O 2.13–2.25 A; 24 In5c and 24 O3c on each face
  - the bottom trilayer (In32O48) is `F F F`
- **Hydroxylation (`substrate.hydroxylate`, also called via `vasp.hydroxylate_checked`):**
  - 3 dissociated H2O per primitive cell, which is 6 pairs per orthorhombic cell, seed 20260926
  - clash filters, now in `substrate.hydroxylate` itself (2026-09-27): no terminal O within 2.5 A of a lattice
    O; every H at least 1.5 A from all atoms and at least 2.5 A from any cation (no acute In–O–H angles)
  - the first draft lacked these filters and placed one In5c class's terminal O 2.07–2.18 A from a lattice O;
    that is fixed, and the pilot slabs `inputs/ito/surfaces/*-oh` were rebuilt with the same routine
  - closest new contact on the pilot slab is 2.25 A (H···H between neighbouring OH), so pilot and DFT slabs now
    share one algorithm
- **Sn cluster:**
  - O_i sits on the empty 16c site of the lower trilayer that faces the slab interior (z = 2.39 A in builder
    coordinates)
  - the two nearest lower-trilayer In become Sn (Sn–O_i 2.40 A; unrelaxed O_i–O 2.22 A)
  - it is geometrically possible. The only catch is that this site sits in the fixed trilayer, so O_i, both Sn
    and their first shell (11 atoms) are freed from selective dynamics
  - Sn/(In+Sn) is 3.1 % in this small cell
- **Physisorption:**
  - the rigid, upright LigParGen conformer (P → tail axis along +z)
  - a grid over tilt (0 or 15°), z-rotation, lateral offset and P height (3.5–4.0 A above the top lattice-O
    plane); the selected P height is 3.5 A
  - it is placed above the same In pair used for the bidentate case, then each acidic H is turned about its
    P–O bond toward the nearest surface O
  - **Result:** one P-OH···O(terminal OH) contact (O···O 2.46/2.53 A, H···O 2.11/2.29 A for Me-4PACz/MeO-2PACz).
    The second P-OH O sits 3.5–3.7 A from the nearest surface O. It is a start geometry for relaxation, not an
    optimised H-bond network.
- **Bidentate:**
  - In5c pair 2/13, In–In 3.83 A; its two open octahedral sites are 2.35 A apart, close to the phosphonate
    O···O distance
  - the two former P-OH oxygens sit on those sites: In–O 2.10–2.15 A
  - the head is rotated about the O···O axis, and a greedy torsion drive runs over the acyclic P→ring bonds. The
    result is P–C tilted 33–37° and the tail axis within 10° of the normal.
  - both protons go to distinct nearby free O3c (acceptor–donor O···O 2.73–3.04 A)
- **Reference geometry:** `mol-<slug>` is exactly the phys molecule translated along z (tested). The first ionic
  step of `mol`, `phys` and `slab-oh` is therefore a single point on precisely the geometries in
  `geometry.extxyz`.

## Energy definitions (eV; all from the same cell, ENCUT and settings)

Let E_x be the final (relaxed; production) total energy of case x, and E1_x the energy of its first ionic step.
E1_x is a single point on `geometry.extxyz`: `OSZICAR` line 1, or the first "energy without entropy" in `OUTCAR`.

- **Physisorption energy** (the value the classical model must reproduce):
  E_ads(phys) = E(ads-phys) − E(slab-oh) − E(mol)
- **Dissociative bidentate adsorption energy** (same reference, so directly comparable):
  E_ads(bid) = E(ads-bid) − E(slab-oh) − E(mol)
- **Chemisorption energy the FF cannot represent** (no reference correction, same atoms):
  ΔE_chem = E(ads-bid) − E(ads-phys) = E_ads(bid) − E_ads(phys)
- **Rigid interaction energy** (smoke runs already give this): ΔE_int(rigid) = E1(ads-phys) − E1(slab-oh) − E1(mol).
  Compare it with the classical E_int = E_complex − E_slab − E_mol evaluated on the same three `geometry.extxyz`
  files. It isolates the interaction terms from relaxation and deformation energy.
- **Deformation terms** (production), for decomposing E_ads:
  - E_def,mol = E(mol, adsorbed geometry) − E(mol, relaxed)
  - E_def,slab = E(slab-oh, adsorbed geometry) − E(slab-oh, relaxed)
  - these need single points on geometries cut from the relaxed ads-phys case, which are not generated here
- **Sn effect** (future): E_ads(phys, Sn) − E_ads(phys). This needs an `ads-<slug>-phys-sn` case above the
  Sn-doped slab. Only the `slab-oh-sn` reference exists now.

## INCAR choices (and deviations from InterfaceForge)

The base is the InterfaceForge NiO(110)+phosphonate relaxation template
(`notebooks/nio_m110_hydroxylation/inputs/vasp_template/INCAR`). That template is PBE, ENCUT 520, PREC Accurate,
ALGO Normal, NELM 150, LASPH, ISMEAR 0/SIGMA 0.05, ISYM 0, IBRION 2/ISIF 2/POTIM 0.3, LDIPOL/IDIPOL 3 with a
structure-specific DIPOL, and NCORE 4.

- **ENCUT = 520 eV:** the InterfaceForge value for all its production relaxations and references. The hardest
  PAW here (O, N, C; ENMAX 400 eV, as recalled for potpaw_PBE.54) gives ENCUT/ENMAX = 1.3. The In_d, Sn_d, P and H
  ENMAX are all at most 255 eV. Matching the IF value keeps these energies consistent with the NiO/phosphonate
  data. The cell is fixed (ISIF = 2), so Pulay stress does not matter; 520 eV is conservative for adsorption
  energy differences. `run.sbatch` prints the real ENMAX lines into `potcar_titles.txt`.
- **Dispersion:** IVDW = 12, D3(BJ), as requested. **This is a deviation:** IF's NiO data use IVDW = 11 (D3 zero
  damping), so energies are not interchangeable with IF's NiO set.
- **EDIFF 1E-5** (IF relax template uses 1E-6). NELMIN = 6 is from IF's reactive-surface recovery settings.
- **LREAL = Auto** (IF relax template: .FALSE.; IF MD template: Auto). This is appropriate for 160–223 atoms. It is
  also used for the molecules, so that all terms of every E_ads share it.
- **ISPIN = 1:** In2O3 is nonmagnetic, 2Sn_In + O_i is ionically compensated (no free carriers), and the acids are
  closed shell. There is no MAGMOM.
- **No DFT+U:**
  - In 4d is a filled semicore shell, treated explicitly by In_d. It is not a correlated open shell like Ni 3d.
  - PBE underestimates the In2O3 gap (about 1 eV vs about 2.9 eV). This is not important for closed-shell
    adsorption geometries or for binding energies of neutral, non-charge-transferring acids at this level.
  - Band alignment or work-function conclusions would need a hybrid functional, or at least a check of the
    carbazole HOMO against the In2O3 VBM.
- **LDIPOL = .TRUE., IDIPOL = 3:**
  - used for every case, including the isolated molecules. The acid is polar and sits in a slab-shaped
    periodic cell, and E_ads needs identical treatment of all three terms.
  - DIPOL is set to the midpoint of the occupied z range, so the correction plane (DIPOL + c/2) falls at the
    middle of the vacuum.
  - AMIN = 0.01 follows IF's slab-alignment settings for LDIPOL charge sloshing.
- **Gamma-only KPOINTS:** the in-plane cell is at least 14.3 A. The run uses `vasp_gam`.
- **Modes:**
  - *smoke:* NSW 3, IBRION 2, EDIFFG −0.03, LWAVE/LCHARG .FALSE., 12 h walltime (4 h for molecules)
  - *production* (`--mode production`): NSW 300, EDIFFG −0.03, WAVECAR/CHGCAR kept, LORBIT 11, 72 h walltime
  - 3 trilayers is `--trilayers 3`, which gives In96O150H12 for slab-oh with 80 atoms still fixed
- **Omitted IF tags:** ADDGRID (not needed for energies; extra cost), and the NiO-specific
  LDAU/MAGMOM/ISPIN = 2 block.

**PAW datasets** (IF `profiles/potcar_pbe_54.yaml`, confirmed identical): In → In_d (4d semicore), Sn → Sn_d,
O → O, P → P, N → N, C → C, H → H.

## How to run on LONI

```bash
cd /path/to/NiO_MD_Passivation && git fetch && git checkout feat/ito-sam-extension && git pull
export IFACE_POTCAR_ROOT=/path/to/potpaw_PBE        # or VASP_PP_PATH; fallback ~/pot/potpaw_PBE
cd inputs/ito/vasp_smoke/slab-oh && sbatch run.sbatch   # one case
# or, after review, all nine:  bash inputs/ito/vasp_smoke/submit_all.sh
```

Each `run.sbatch` does the following, in order:
1. `sha256sum -c inputs.sha256` confirms the inputs are byte-identical to what was generated.
2. It builds POTCAR from `POTCAR.spec`, using the InterfaceForge search order
   (`IFACE_POTCAR_ROOT` → `$VASP_PP_PATH/potpaw_PBE` → `$VASP_PP_PATH` → `~/pot/potpaw_PBE`). An existing
   non-empty POTCAR is never overwritten, and it stops if any dataset is missing.
3. It checks the ZVAL order and the NELECT total against the spec.
4. It loads `vasp6/6.5.1-cpu` (override with `IFACE_VASP_MODULE`) and runs `srun vasp_gam` with
   OMP_NUM_THREADS = 1.
5. It writes `run_provenance.txt`: job, host, module, repo commit, input and POTCAR hashes, and time taken.

Resources:
- slab and adsorbate cases: `workq`, 1 node × 64 tasks
- molecules: `single`, 16 tasks
- account: `loni_perovsk27`

The generated POTCAR must not be committed. Add `POTCAR` (and `WAVECAR`/`CHGCAR`) to `.gitignore` before
committing results.

Regenerate or verify locally with `python scripts/ito/prepare_vasp_smokes.py [--check]`.

## What to bring back

Per case, bring back:
- `OSZICAR`, `OUTCAR` (gzip), `CONTCAR`, `vasprun.xml` (gzip)
- `vasp.cpu.<jobid>.out`, `run_provenance.txt`, `potcar_titles.txt`

Do not bring back POTCAR, WAVECAR or CHGCAR.

From the smoke runs, record:
- the first-step and final energies
- the number of SCF steps per ionic step
- `LOOP`/`LOOP+` timings and the "total amount of memory used", to calibrate the production cost estimate
- the largest forces on the free trilayer and on the adsorbate
- any SCF non-convergence or ZBRENT warnings

## Limitations and unverified items

- **No calculations were run.** All energies, timings and memory figures above are estimates or targets.
- **ZVAL and ENMAX** values are standard potpaw_PBE.54 values recalled here, not read from a licensed POTCAR. The
  sbatch verifies ZVAL at runtime. The LONI `single` partition limit for 16 tasks is not verified.
- **Rigid LigParGen conformers.** The Me-4PACz conformer is compact (maximum intramolecular distance 11.3 A), so
  the upright molecules are only about 10 A tall rather than about 15 A. The common cell height is therefore
  30 A, which still leaves at least 12 A of vacuum above every molecule and at least 23 A above the slabs.
- **Molecule reference.** `mol-*` is the physisorbed conformer, not a searched gas-phase minimum. In production
  it only relaxes locally.
- **Physisorbed start** has one good P-OH···O contact. The second P-OH O is 3.5–3.7 A from the surface. The
  smoke NSW = 3 will not optimise this, so smoke energies apart from the step-1 single points are not adsorption
  energies.
- **Geometric starting structures.** The bidentate start and the O_i are geometric constructions, not
  DFT-relaxed. The O_i starts with an unrelaxed O–O distance of 2.22 A.
- **Thin slab.** Two trilayers (5.8 A) with one fixed is too thin for converged surface energetics; use
  `--trilayers 3` for production. The 1x1 cell is the low-coverage limit: 1 molecule per 3.55 nm², that is
  0.28 nm⁻², compared with the pilot dose of 3.45 nm⁻².
- **Sn content.** Only one Sn cluster: 3.1 % Sn, below the 9.3 % ITO target. There is no free-carrier physics
  (by design, as in the classical model).
- **Line endings.** The manifests store SHA256 of LF files. This repo has `core.autocrlf = true` on Windows, so
  consider adding `inputs/ito/vasp_smoke/** -text` to `.gitattributes`, as was done for the ITO surfaces.
