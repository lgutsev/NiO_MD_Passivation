# ITO/SAM extension: status note (work in progress, 2026-09-26)

Branch `feat/ito-sam-extension` (from `origin/main` 4c82afa). Background task; the NiO campaign stays the
priority. **No NiO code path, input, or result is changed.** The only edits outside the new ITO files are additive
mass/radius entries (In, Sn, I) in `geometry.ELEMENTS` and `analysis/coverage.py`.

This is an interim commit made because the session is running out of usage. The first deliverable is only partly done. See "Remaining" below.

## What exists

| Path | Content |
|---|---|
| `src/nio_md_prep/ito/substrate.py` | Deterministic bixbyite In2O3(111) slab builder with checks, dissociative hydroxylation, and Sn doping as (2 Sn_In + O_i) |
| `src/nio_md_prep/ito/forcefield.py` | Two named, source-tagged slab LJ sets that bracket the unknown In/Sn dispersion |
| `src/nio_md_prep/ito/assemble.py` | Pilot assembly and LAMMPS stage inputs, with a zero-step validation |
| `src/nio_md_prep/ito/__main__.py` | `python -m nio_md_prep.ito build-substrate / build-pilot` |
| `inputs/ito/surfaces/{in2o3-111-bare,in2o3-111-oh,ito-111-oh}/` | `model.toml` plus the built `surface.lmp`/`.xyz` and `surface_manifest.json` (hashes, counts, dipole, sites) |
| `studies/ito/pilot-me-4pacz.toml` | Pilot 1: single-component SAM |
| `docs/ito/literature-matrix.md` | **Partial.** The literature agent hit the usage limit; only a few rows are in |

The `adsorption-scan` and `analyze` subcommands are wired in the CLI, but `checks.py` and `analysis.py` are **not written yet**.

## Repository audit (summary)

**Reusable unchanged:**
- `lammps.parse/replicate/write`
- `chemistry.correction_lines` (Cao phosphonate LJ and torsion corrections)
- `build._packmol` / `_surface` / `_coeff_lines`
- LigParGen molecule manifests
- coverage rasterisation, height maps, void patches and roughness (`analysis/coverage.py`, once the In/Sn masses are present)
- the periodic geometry, tilt, RDF and site-exchange helpers in `analysis/interfacial.py`
- `model_scope` reporting policy

**NiO-specific assumptions that must be replaced for ITO:**
- **Force field in `build.py`:**
  - hard-coded `HEADER` with hybrid Buckingham styles (`build.py:12`)
  - Ni/O identified by mass (`:151`)
  - Ni–O Buckingham and Cao pseudo-LJ Ni/O cross terms (`:400-414`)
- **Surface detection:** surface atoms are found by `|q|==2` (`build.py:381,483`; `validate.py:107`).
- **Validation (`validate.py:109-114`):** it requires exactly 21060 surface atoms, a fixed surface SHA, and two surface types.
- **Defaults and ratio report in `build.py`:**
  - packing region fixed to the 125.1 x 41.7 A box (`:297`)
  - ratio report assumes Me-4PACz is primary and 0.5/0.3 mg/mL stocks (`:421`)
  - NPT with x/y coupling, which is incompatible with a rigid crystalline slab
- **Anchor detection:** every molecule must carry a phosphonate (`chemistry.phosphonate_roles`), so PyCA-3F (a carboxylic acid) fails.
- **Interfacial analysis:**
  - "exposed Ni" sites (`interfacial.py:255-375`)
  - hard-coded NiO reference surface path (`:314`)
  - P-only anchors (`:427-447`)
  - Ni–O(P) contact cutoff 3.25 A (`:737`)
  - NiO wording in the report and publication modules
  - Me-4PACz-primary labels (`publication_report.py:157-190`)
- **Scripts:** the sbatch scripts hard-code the NiO study lists and array sizes.

These are handled by **not** routing ITO through `build.build()`. The ITO subpackage reuses the generic pieces instead.

## Proposed ITO model (implemented)

- **Structure:** bixbyite Ia-3, a = 10.117 A (Marezio 1966). The slab uses an orthorhombic (111) setting, [1-10] x [11-2] x [111], 14.31 x 24.78 A per cell. It is cut between neutral, symmetric O–In–O trilayers:
  - 6 trilayers, 17.5 A thick, In1536O2304
  - 4x2 cells, 28.37 nm²
  - stoichiometric, dipole-free by construction (checked, about 1e-16 e/A)
  - top face: 6.77 In5c and 6.77 O3c per nm²
- **Hydroxylation (pilot default):**
  - 3 dissociated H2O per primitive (111) cell (1.77 nm²), giving 3.38 OH/nm²
  - each pair is an OH on an In5c plus an H on the nearest O3c
  - placement is geometric, not DFT-relaxed
  - the density follows the ordered-hydroxyl picture reported for In2O3(111) by the Diebold group; the exact reference is still **to be verified** in the literature matrix
- **Explicit ITO:** Sn/(In+Sn) = 0.0924, which targets 90:10 wt% In2O3:SnO2. Sn is placed as neutral Frank–Köstlin-type (2 Sn_In + O_i) clusters, with O_i on empty 16c sites kept away from both faces.
- **What cannot be represented:**
  - free carriers
  - oxygen vacancies, which would need electronic compensation
- **Charges:** formal charge x 0.525. With hydroxyl H = 0.425 this reproduces the CLAYFF pattern exactly (O −1.05, OH-O −0.95, H +0.425, In +1.575, Sn +2.1). It keeps hydroxylation and Sn doping exactly neutral.
- **Slab dynamics:** the slab is rigid and its internal interactions are excluded. There is no validated flexible In2O3 potential in the pipeline. A rigid crystal is the conservative choice for a packing model.
- **SAM–slab interaction:** geometric mixing of OPLS ligand terms with the slab self terms. Two bracketing sets:
  - `uff-cation/clayff-anion`: UFF In/Sn (Rappé 1992, verified against the Open Babel `UFF.prm` table) with CLAYFF O
  - `clayff-cation/clayff-anion`: CLAYFF octahedral-Al analogue, so cation dispersion is effectively off
- **Anchoring:** the phosphonic acids stay protonated and neutral (LigParGen charges). Binding is physisorption and H-bonding only.
- **What this nonreactive model cannot establish:**
  - P–O–In bond formation or condensation
  - deprotonation state
  - mono-, bi- or tridentate binding mode
  - charge transfer, polarisation, or work-function shifts

## Pilot status

- Pilot 1 is Me-4PACz on `in2o3-111-oh`: 98 molecules, the NiO areal dose of 3.45 nm⁻². It **builds and passes** the zero-step LAMMPS check. The ligand–slab minimum distance is at least 3 A and the total charge is −0.0098 (LigParGen rounding).
- The build directory is outside git: `C:\Users\lguts\nio-wt\ito-runs\pilot-me-4pacz\s11`.
- **Local benchmark** (Windows, 4 OpenMP threads): 4.5 steps/s for 8394 atoms. The full pilot protocol (700k steps at 0.5 fs, 350 ps) would take about 43 h per placement locally.
- **Before any campaign:**
  - shrink zhi from 200 A to about 120 A, which cuts the PPPM grid
  - consider SHAKE with a 1 fs timestep
  - benchmark on one 64-core QBD `workq` node

## Remaining (not done)

1. Finish `literature-matrix.md`:
   - 2PACz, MeO-2PACz, Me-4PACz
   - the mixtures MeO-2PACz/I-2PACz, MeO-2PACz/Me-4PACz, 2PACz/PyCA-3F, Me-4PACz/PPA
   - classify ITO/SAM vs ITO/NiOx/SAM
   - record preparation and ratios
   - extract the 2PACz/PyCA-3F MD methodology
2. Choose the mixture from the literature ratio. The planned pilot 2 is MeO-2PACz/Me-4PACz, since both molecules are already parameterised. It needs a literature-justified ratio.
3. `ito/checks.py`: single-molecule adsorption scans (several placements and tilts, frozen slab, E_int = E_complex − E_slab − E_mol) for bare, OH and Sn-doped slabs under both LJ sets.
4. `ito/analysis.py`: contacts to In/Sn and H-bonds to surface OH, tilt, P-head clustering, and coverage via the existing coverage module, across 3 or more Packmol seeds.
5. LigParGen inputs the user must generate (the package never calls LigParGen): 2PACz, I-2PACz, PyCA-3F (no phosphonate, so the anchor logic needs a carboxylic-acid path), and phenylphosphonic acid derivatives.
6. Tests for the substrate builder: stoichiometry, neutrality, dipole, determinism.

**Stop point:** do not launch any adsorption or coverage campaign before review.
