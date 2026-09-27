# ITO/SAM extension: first deliverable (for review)

Branch `feat/ito-sam-extension`, from `origin/main` 4c82afa. This is background work; the NiO campaign stays the priority.

**NiO safety:**
- No NiO code path, input or result is changed.
- Shared modules got additive element entries only: In, Sn and I in `geometry.ELEMENTS` and `analysis/coverage.py`.
- The affected NiO tests pass: `test_coverage`, `test_build`, `test_chemistry`, `test_interfacial` and `test_lammps`, 48 passed.

**Stop point:** no adsorption or coverage campaign has been launched. The classical pilots and the VASP smokes are prepared for QBD and wait for review.

## Deliverable map

| Item requested | Where |
|---|---|
| Literature matrix (2PACz, MeO-2PACz, Me-4PACz; the four mixtures; ITO/SAM vs ITO/NiOx/SAM; preparation and ratios) | [`literature-matrix.md`](literature-matrix.md) |
| Repository audit: reusable pieces and every NiO-specific assumption | [`repository-audit.md`](repository-audit.md) |
| Proposed ITO model: proxy vs explicit Sn, termination, hydroxylation, vacancies, charges, SAM–surface terms, the PyCA-3F MD reference and what a nonreactive FF cannot establish | [`ito-model.md`](ito-model.md) |
| Parameter sources | [`parameter-sources.md`](parameter-sources.md) |
| Pilot plan and compute estimate | [`pilot-plan.md`](pilot-plan.md) |
| Single-molecule geometry and energy checks | [`adsorption-scan.md`](adsorption-scan.md) (local, classical) |
| DFT reference inputs (prepared, not run) | [`vasp-smoke.md`](vasp-smoke.md), `inputs/ito/vasp_smoke/` |

## Code and inputs

| Path | Content |
|---|---|
| `src/nio_md_prep/ito/substrate.py` | Bixbyite In2O3(111) slabs from COD 2310009 (checked stoichiometric, neutral and dipole-free); clash-filtered dissociative hydroxylation; (2 Sn_In + O_i) doping |
| `src/nio_md_prep/ito/forcefield.py` | Two source-tagged slab LJ sets that bracket In/Sn dispersion; CLAYFF-scaled charges |
| `src/nio_md_prep/ito/assemble.py` | Pilot assembly reusing the NiO LigParGen, Packmol and Cao-correction code. Rigid slab, NVT, SHAKE at 1 fs, PPPM `diff ad`, zero-step LAMMPS validation, full hash provenance |
| `src/nio_md_prep/ito/checks.py` | Single-molecule placement scan with exact fixed-mesh PPPM energy decomposition |
| `src/nio_md_prep/ito/analysis.py` | Contacts, H-bonds, tilt, P-head clustering, and per-component residence; coverage through the unchanged NiO coverage module |
| `src/nio_md_prep/ito/vasp.py`, `scripts/ito/prepare_vasp_smokes.py` | Nine VASP reference cases; POTCAR assembled only on the cluster |
| `inputs/ito/surfaces/{in2o3-111-bare,in2o3-111-oh,ito-111-oh}/` | `model.toml` → `surface.lmp`/`.xyz` plus `surface_manifest.json` |
| `studies/ito/` | `smoke-me-4pacz`, `pilot-me-4pacz`, `pilot-meo-2pacz-me-4pacz-1to1`, `pilot-me-4pacz-monolayer`, `adsorption-scan` |
| `scripts/ito/{build_pilots,run_pilot_array}.sbatch` | QBD launchers; smoke is the default target |
| `tests/test_ito_substrate.py`, `tests/test_ito_vasp.py` | 29 tests |

CLI: `python -m nio_md_prep.ito {build-substrate, build-pilot, adsorption-scan, analyze}`.

## Corrections made during this work

- **Bixbyite coordinates:** the first draft used O(48e) values transcribed from memory, off by up to 0.04 Å. They were replaced with the COD 2310009 (Marezio 1966) values and every slab was rebuilt.
- **Hydroxylation:** the first draft placed one In5c class's terminal O 2.07–2.18 Å from a lattice O, and some protons at acute In–O–H angles. Clash filters were added, the slabs rebuilt, and the scan rerun.

## Decisions needed at review

1. Approve the QBD pipeline smoke, then the 3 pilots × 3 seeds (`pilot-plan.md`, about 5–14 node-hours, an estimate).
2. Approve the 9 VASP smoke cases (about 3 node-hours, an estimate). The D3 variant: IVDW=12 was used here; InterfaceForge's NiO data use IVDW=11.
3. Keep or change the default slab parameter set. The scan shows how strongly this choice matters (`adsorption-scan.md`).
4. Future molecules need LigParGen files: 2PACz, I-2PACz, PyCA-3F (plus a carboxylic-acid anchor path), and phenylphosphonic acids.

## Checkpoint 2026-09-28 (resume here)

**Done since the first deliverable:**
- corrugated ITO slab matching the NiO groove (`in2o3-111-oh-groove`)
- hydroxylation series 0/25/50/75/100 % for corrugated ITO (`in2o3-111-groove-ohXXX`; 75 = 100, saturated at 64 %) and rigid corrugated NiO (`inputs/surfaces/corrugated-nio-110-rigid-ohXXX`)
- ball-and-stick renderer (`scripts/ito/render_surfaces.py`, figures in `data/`)
- 18 hydroxylation-series pilot studies (`studies/ito/hydroxylation-series.txt`) plus `scripts/ito/submit_hydroxylation_series.sh`
- ITO groove-vs-hydroxylation energy scan (`adsorption-scan.md`, last section): hydroxylation kills plateau binding, but the groove floor stays strongly binding
- the CLAYFF-cation set is invalid at exposed rim cations; a bug that left Ni out of the cation contact sets is fixed

**In flight / next:**
1. **NiO groove-vs-hydroxylation scan.** It was relaunched locally on 2026-09-28 (output `C:\Users\lguts\nio-wt\ito-runs\scan-nio-groove-oh`, log `scan-nio-groove-oh.log`). If it is not finished, rerun:
   `cd src && python -m nio_md_prep.ito adsorption-scan ../studies/ito/adsorption-scan-nio-groove-hydroxylation.toml --output <dir> --workers 6`
   Then add its table to `adsorption-scan.md` next to the ITO one.
2. **HPC, awaiting user approval:** smoke first, then the pilots and the VASP smokes (see `pilot-plan.md`).
3. **Not started:**
   - the carboxylic-acid anchor path (PyCA-3F)
   - LigParGen files for 2PACz, I-2PACz, PyCA-3F and PPA derivatives (user-supplied)
