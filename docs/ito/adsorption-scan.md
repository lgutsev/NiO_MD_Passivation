# Single-molecule placement and energy checks (classical FF, local)

These are geometry and energy checks, not a campaign. They answer three questions before any film simulation:
- Does the substrate model matter?
- Does the slab parameter set matter?
- Does the force field support surface anchoring at all?

- **Run:** 2026-09-27, `python -m nio_md_prep.ito adsorption-scan studies/ito/adsorption-scan.toml`, on the corrected COD / clash-free slabs.
- **Raw rows:** [`data/adsorption-scan-results.json`](data/adsorption-scan-results.json). Those rows predate a key fix: `tilt_deg` there is the **final** tilt, and `initial_tilt_deg` was reconstructed from the placement index.

## Method

Each of the 108 placements is one neutral, protonated molecule:
- slabs: `in2o3-111-bare`, `in2o3-111-oh`, `ito-111-oh`
- molecules: Me-4PACz, MeO-2PACz
- both slab parameter sets
- 3 lateral sites × 3 initial tilts (0°, 45°, 80°), anchor down, lowest atom 2.5 Å above the top atom

Each placement is relaxed with the slab frozen: CG (≤1500 iterations), then a 0.5 ps NVT quench at 300 K, then CG again. The energies are:

```
E_int = E(complex) − E(slab) − E(molecule at the complex geometry)
E_ads = E(complex) − E(slab) − E(lowest isolated-molecule energy found)
```

- All terms use one box and a fixed PPPM mesh and g_ewald, so slab–slab terms cancel exactly. The regression test puts E_int = 0.005 kcal/mol at 35 Å separation.
- **E_int is the reliable number.** E_ads carries up to 48 kcal/mol of unconverged conformational deformation, because CG is capped and the quench is short. E_ads is reported, but it should not be interpreted.

## Results (kcal/mol; median over 9 placements, minimum in parentheses)

| Substrate | Molecule | UFF-cation E_int | CLAYFF-cation E_int | Anchor O–In/Sn < 3.25 Å | Anchor H-bond < 2.5 Å | Median P height above top atom (Å) | Median final tilt (°) |
|---|---|---|---|---|---|---|---|
| bare In2O3(111) | Me-4PACz | −48.1 | −33.3 | 22–44% of placements | 0–11% | 2.8–2.9 | 46 |
| bare In2O3(111) | MeO-2PACz | −47.8 | −29.8 | same range | 0% | 3.2–3.4 | 50–54 |
| In2O3(111)-OH | Me-4PACz | −6.2 | −3.2 | 0% | 0% | 4.2–4.7 | 54 |
| In2O3(111)-OH | MeO-2PACz | −4.8 | −3.1 | 0% | 0% | 5.0–5.1 | 65 |
| ITO(111)-OH (explicit Sn) | Me-4PACz | −4.2 | −3.1 | 0% | 0% | 4.8–4.9 | 56 |
| ITO(111)-OH (explicit Sn) | MeO-2PACz | −5.5 | −3.8 | 0% | 0% | 5.0–5.1 | 64 |

Pooled over both parameter sets and both molecules, E_int by initial tilt:

| Substrate | 0° (upright) | 45° | 80° (flat) |
|---|---|---|---|
| bare | −64.5 (−87.8) | −42.2 (−70.6) | −14.3 (−55.9) |
| In2O3-OH | −7.4 (−22.5) | −3.2 (−8.6) | −1.2 (−11.0) |
| ITO-OH | −4.7 (−18.3) | −3.6 (−9.3) | −3.2 (−7.8) |

**H-bonded start pose.** The classical energy was also evaluated at the physisorbed geometry built for DFT (`scripts/ito/classical_on_vasp_pose.py`). That pose has one P–OH···O H-bond on the 1×1 hydroxylated slab.

| Case | Parameter set | E_int at the DFT start geometry | E_int after minimizing the molecule |
|---|---|---|---|
| Me-4PACz | UFF | +4.4 | −11.4 |
| Me-4PACz | CLAYFF | +10.8 | −9.2 |
| MeO-2PACz | UFF | +4.2 | −13.4 |
| MeO-2PACz | CLAYFF | +9.4 | −8.0 |

## What this shows

1. **The substrate model dominates.** The bare slab binds the protonated acids 5–10× more strongly than the hydroxylated slabs, and it brings phosphonate O within 3.25 Å of In. The hydroxylated surfaces give weak physisorption, about −3 to −13 kcal/mol, with P 4–5 Å above the top atom and no local minimum with an anchor–surface contact. So the choice between bare and hydroxylated is not a detail. It decides whether the classical model has any anchoring at all.
2. **The parameter set matters mainly on the bare slab.** The UFF-cation set is 15–18 kcal/mol deeper there. On the hydroxylated slabs the two sets differ by only about 1–3 kcal/mol, because the OH layer screens the cations.
3. **Explicit Sn is not resolved at this level.** ITO-OH and In2O3-OH differ by less than the placement-to-placement scatter. Keeping the proxy as the pilot default is consistent with this.
4. **Anchoring on OH surfaces is weak and not H-bond-driven in local optimization.** No scan placement ended with an anchor H-bond; the closest acid-H to slab-O distance is 3.8 Å. The DFT-style H-bonded start relaxes to −8 to −13 kcal/mol, so the force field allows weak H-bonded physisorption. Local minimization from upright starts simply does not find it.

## Consequences for the pilots (review items)

- **Expectation:** on `in2o3-111-oh` the moving wall, not surface affinity, will drive film formation. Expect loosely bound first layers. The pilot therefore tests packing under confinement, which is also what the NiO protocol does, but anchoring statistics will be near zero unless contacts form under compression.
- **Decision 1:** should a **bare-slab arm** (the Li et al. 2024 choice, no OH) run as the "strong-anchoring" bracket next to the hydroxylated default?
- **Decision 2:** DFT calibration of the two numbers the classical model depends on most:
  - the `ads-*-phys` VASP case gives the DFT physisorption energy for the H-bonded pose, to compare against −8 to −13 kcal/mol here;
  - `ads-*-bidentate` gives the chemisorption energy that the protonated force field cannot represent at all.
  If DFT physisorption is much stronger than the force field's, the slab LJ or charges need refitting before any coverage claim.
- **E_ads caveat:** E_ads should not be used until the minimization caps are raised. That costs roughly 3–5× the CPU.

## Corrugated slab, site-resolved (first pass, 2026-09-27)

**Setup:** `studies/ito/adsorption-scan-groove.toml` on `in2o3-111-oh-groove`, the first build, in which CN 3 step-edge In were not yet eligible for hydroxylation. Sites are groove floor, wall midpoint and plateau centre, with 3 tilts each.

| Molecule | Set | Floor E_int (median) | Wall | Plateau |
|---|---|---|---|---|
| Me-4PACz | UFF-cation | −103 | −113 | −6.8 |
| MeO-2PACz | UFF-cation | −44 | −50 | −2.3 |
| Me-4PACz | CLAYFF-cation | −63 | −101 (min −216) | −3.8 |
| MeO-2PACz | CLAYFF-cation | −32 | −202 | −1.1 |

- **The CLAYFF-cation set is not usable on corrugated slabs.** In the strongest placements, phosphonate O sits 1.71–1.84 Å from CN 3–4 rim In, far below any real In–O bond (about 2.1–2.2 Å). The CLAYFF "ao" cation has essentially no core repulsion (ε ≈ 1e-6 kcal/mol), and on flat O/OH-terminated slabs the anion layer hides this. The set is now flagged in `forcefield.py`.
- **The UFF-cation set keeps all molecule–slab contacts ≥ 2.4 Å.** It still binds groove sites 15–50× more strongly than the plateau. Exposed low-CN step-edge In provides both Coulomb attraction and multi-sided dispersion.
- **The groove preference therefore depends on how the rim is terminated.** The follow-up scans on the hydroxylation series answer this; CN 3 In are now hydroxylatable:
  - `adsorption-scan-groove-hydroxylation.toml` (ITO 0/25/100%)
  - `adsorption-scan-nio-groove-hydroxylation.toml` (rigid NiO 0/25/100%)

## Corrugated ITO vs hydroxylation (UFF-cation set, 2026-09-28)

**Setup:** `studies/ito/adsorption-scan-groove-hydroxylation.toml` on `in2o3-111-groove-oh000/025/100`. The 100% surface is the saturated one, 64% achieved. Sites are groove floor, wall midpoint and plateau centre, with 3 tilts each. Heights are measured from the local surface. Raw rows: `data/adsorption-scan-ito-groove-hydroxylation.json`.

Median E_int (kcal/mol) over 3 placements; the two numbers are Me-4PACz / MeO-2PACz:

| Site | 0% OH | 25% OH | 100% (sat.) OH |
|---|---|---|---|
| Groove floor | −92 / −63 | −62 / −49 | −88 / −51 |
| Groove wall | −38 / −35 | −17 / −16 | −17 / −19 |
| Plateau | −54 / −43 | −8 / −4 | +5 / −2 |

- **Hydroxylation removes plateau binding** (from about −50 to about 0). On the bare slab, the plateau binding came from anchor O contacting exposed In (2–3 of 3 placements).
- **The groove floor stays strongly binding at every level.** At 25% it does so with no anchor–cation contacts at all. The concave floor encloses the molecule, so dispersion and electrostatics act from several sides. This is a geometric effect, not an anchoring one.
- **Within this model, the groove/plateau preference therefore grows with hydroxylation.** That is consistent with the hypothesis that SAMs accumulate in corrugations.
- **Caveats:**
  - only 3 local minimizations per cell;
  - the floor sites sit on ideal, unrelaxed trilayer steps;
  - these are classical physisorption energies, not binding energies.

  The MD pilots on the same surfaces, with groove/plateau densities, are the real test.
