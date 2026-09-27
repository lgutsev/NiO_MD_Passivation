# ITO pilot plan and compute estimate

**Nothing below has been launched.** The HPC pipeline smoke (step 1) is the only run proposed before review. The full adsorption and coverage campaign waits for review.

## Systems

All systems use the rigid, hydroxylated In2O3(111) slab `in2o3-111-oh`:
- 57.23 × 49.56 Å, 28.37 nm²
- 3984 atoms
- 3.38 OH/nm²
- parameter set `uff-cation/clayff-anion`

| Study (`studies/ito/`) | Molecules | Dose (nm⁻²) | Atoms | Why |
|---|---|---|---|---|
| `smoke-me-4pacz` | 98 Me-4PACz | 3.45 | 8394 | Pipeline smoke only (70 ps): build → MD → analysis on QBD, and a benchmark |
| `pilot-me-4pacz` | 98 Me-4PACz | 3.45 | 8394 | Single component; dose matched to the NiO baseline (180 on 52.17 nm²) for a direct substrate comparison |
| `pilot-meo-2pacz-me-4pacz-1to1` | 49 + 49 | 3.45 | 8198 | Justified mixture: 1:1 is the most frequent verified ratio (B2-4, B2-5); by mass ≈ by moles (MW within 1%) |
| `pilot-me-4pacz-monolayer` | 54 Me-4PACz | 1.90 | 6414 | Monolayer-scale dose inside the measured X-2PACz range 1.3–2.3 nm⁻². The 3.45 dose is an excess, like the NiO campaign |
| `smoke-me-4pacz-groove` | 184 Me-4PACz | 3.45 | ~18.5k | Pipeline smoke on the corrugated slab (70 ps) |
| `pilot-me-4pacz-groove` | 184 Me-4PACz | 3.45 | ~18.5k | Corrugated counterpart of the NiO baseline: tests groove accumulation |
| `pilot-meo-2pacz-me-4pacz-1to1-groove` | 92 + 92 | 3.45 | 18137 | Corrugated mixture: per-component groove vs plateau occupancy |

Each study uses three Packmol / velocity seeds (11/111, 12/112, 13/113), so the initial placements are independent.

Follow-ups, only after review:
- MeO-2PACz:Me-4PACz 1:2 (B2-2, the verified direct ITO/SAM stack)
- `ito-111-oh` (explicit Sn) and the `clayff-cation` parameter set as sensitivity arms
- a higher-OH slab for plasma-cleaned ITO
- 2PACz, I-2PACz and PyCA-3F once LigParGen files exist; PyCA-3F also needs a carboxylic-acid anchor path

## Protocol (per seed)

All stages use a rigid slab, NVT at 300 K (Nosé–Hoover, 100 fs) on the ligands, SHAKE on ligand X–H bonds, and a 1 fs timestep.

| Stage | Steps | Time | Walls |
|---|---|---|---|
| Minimize | SD 5k + CG 20k | — | — |
| Deposition | 300k | 300 ps | upper LJ wall ramps from 4 Å above the packed film to 30 Å above the slab top |
| Compressed hold | 200k | 200 ps | wall fixed at 30 Å |
| Release | 100k | 100 ps | wall retracts by 40 Å |
| Relaxed hold | 300k | 300 ps | analysis window |

Total: 0.9 ns. The literature MD for 2PACz/PyCA-3F used 1 ns, and Park et al. used 1 + 10 ns (matrix §B3). An aggregation claim therefore needs the relaxed hold extended to at least 5–10 ns once the pilot shows the protocol is sound.

## Observables (`python -m nio_md_prep.ito analyze`)

| Observable | Definition |
|---|---|
| Coverage | `analysis.coverage` unchanged: total projected coverage, near-surface (≤5 Å) coverage, P-anchored coverage, void patches, roughness |
| Surface residence | P within 6 Å of the slab top atom; surface density (nm⁻²) vs the experimental 1.3–2.3 nm⁻² |
| Contacts | phosphonate O within 3.25 Å of In/Sn; anchor H-bonds < 2.5 Å (P–OH···O_slab, P=O···H–O_slab). These are proximity measures, not bonds |
| Orientation | tilt of the P→core vector from the normal (surface-resident molecules) vs NEXAFS 61–65° carbazole-plane angle (matrix A3; the angle definitions differ, so compare with care) |
| Clustering | periodic single-linkage of surface P heads (7 Å) plus the fraction of molecules stranded above the first layer |
| Mixture | per-component surface-resident and anchored fractions, i.e. the adsorbed composition vs the 1:1 solution |
| Groove accumulation (corrugated) | surface-resident molecules per projected nm² on groove floor, walls and plateau; groove/plateau density ratio; per-component groove fraction. P height is measured from the local slab surface |

The mean over 3 seeds, and the scatter between them, is the unit of evidence. The NiO campaign reports block SEMs within runs, so seed-to-seed scatter is new information here.

## Compute estimate

**Local measurements** (Windows, 8394 atoms, per 100 MD steps):

| Setting | Time per 100 steps | Rate |
|---|---|---|
| PPPM `ik`, 1 thread | 8.2 s | ~12 steps/s |
| PPPM `diff ad`, 1 thread | 4.4 s | ~23 steps/s |
| PPPM `diff ad`, zhi 125 Å, 1 thread | 3.9 s | ~26 steps/s |

- 4 OpenMP threads helped little, because PPPM is 75–90% of the time.
- One pilot run is 0.9 M steps. Serially that would be about 10 h.

**QBD** (64-rank `workq` node): not measured. Assuming 10–25× the serial rate, a PPPM-bound 8k-atom system with a thin slab typically scales poorly past about 32 ranks.

| Item | Estimate |
|---|---|
| One pilot run | ~0.5–1.5 h |
| 9 flat pilot runs (3 studies × 3 seeds) | ~5–14 node-hours |
| 6 corrugated pilot runs (2 studies × 3 seeds, ~2.2× the atoms) | ~7–20 node-hours |
| Smoke (3 × 70 ps) | < 1 node-hour; it replaces these guesses with a measured `Performance:` line |

If scaling is poor, the fallback is to run 2–4 jobs per node with `-n 16`, using packed sub-node arrays.

Single-molecule energy scan (`studies/ito/adsorption-scan.toml`), measured locally: 144 placements, about 5 min each on one core with a 1 ps quench, so about 12 CPU-hours (about 2 h on 6 workers).

VASP smokes: see [`vasp-smoke.md`](vasp-smoke.md) for per-case electron counts and wall-time estimates.

## HPC commands (after `git pull` on QBD; nothing is submitted automatically)

```bash
sbatch scripts/ito/build_pilots.sbatch                    # builds smoke + 3 pilots x 3 seeds under prepared/ito/
sbatch scripts/ito/run_pilot_array.sbatch                 # smoke only (default STUDY=smoke-me-4pacz)
# after review of the smoke's Performance lines and analysis summaries:
STUDY=pilot-me-4pacz sbatch scripts/ito/run_pilot_array.sbatch
STUDY=pilot-meo-2pacz-me-4pacz-1to1 sbatch scripts/ito/run_pilot_array.sbatch
STUDY=pilot-me-4pacz-monolayer sbatch scripts/ito/run_pilot_array.sbatch
```

The QBD `.venv` needs `pip install -e '.[analysis]'` (numpy, scipy) for the analysis step. ASE is needed only to rebuild slabs; the built slabs are committed.

## Stop criteria for the smoke (all must hold before the pilots run)

- every stage finishes, and the ligand temperature stays at 300 ± 15 K after the first 10 ps
- there are no lost atoms and no SHAKE failures
- during the compressed hold, at least some molecules reach the surface (surface-resident fraction > 0)
- analysis JSON is produced for both the hold and the relax trajectories
- the measured `Performance:` line is recorded here, replacing the estimates above
