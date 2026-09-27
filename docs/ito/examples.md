# ITO extension: worked examples

Every command runs from the repository root, except where noted. Scans and renders are local CPU jobs; MD pilots and VASP cases are for QBD only.

## 1. Build and look at a surface

```bash
cd src && python -m nio_md_prep.ito build-substrate ../inputs/ito/surfaces/in2o3-111-groove-oh025/model.toml
cd .. && python scripts/ito/render_surfaces.py out.png inputs/ito/surfaces/in2o3-111-groove-oh{000,025,050,075} --labels "0%" "25%" "50%" "75%"
```

The build rewrites `surface.lmp`, `surface.xyz` and `surface_manifest.json`. The manifest holds the hashes, counts, charges, dipole, hydroxylation shortfall and corrugation record. The render shows a ball-and-stick cross-section with bonds plus a top view of the hydroxyl positions. Examples: `data/ito_hydroxylation_series.png`, `data/nio_hydroxylation_series.png`.

## 2. Single-molecule energy scans (`studies/ito/*.toml`, `studies/ito/examples/*.toml`)

| Config | What it asks |
|---|---|
| `adsorption-scan.toml` | Flat slabs: bare vs OH vs Sn-doped, both slab LJ sets, random sites |
| `adsorption-scan-groove.toml` | Groove floor / wall / plateau on the pilot groove slab |
| `adsorption-scan-groove-hydroxylation.toml` | The same three sites on ITO at 0 / 25 / 100 % OH |
| `adsorption-scan-nio-groove-hydroxylation.toml` | The same on rigid NiO at 0 / 25 / 100 % OH |
| `examples/hydroxyl-h-core-s040.toml`, `-s107.toml` | Parameter-override check (section 3) |

```bash
cd src && python -m nio_md_prep.ito adsorption-scan ../studies/ito/adsorption-scan-groove-hydroxylation.toml --output ../scan-out --workers 6
```

Output: `scan-out/scan_results.json`, which holds one row per placement plus a summary, and one minimized `.xyz` per placement.

**Reading a row:**
- `e_int_kcal_mol` is the reliable number: it is computed at fixed geometry with exact fixed-mesh PPPM decomposition.
- `e_ads_kcal_mol` includes unconverged conformational deformation, so do not interpret it.
- Contact fields let you screen for artefacts. Real In/Ni–O bonds are about 2.0–2.2 Å, and hydrogen-bond H···O distances are about 1.6–2.0 Å. Anything shorter means the force field let atoms collapse:
  - `anchor_o_min_cation_angstrom`
  - `anchor_o_min_surface_h_angstrom`
  - `acid_h_min_surface_o_angstrom`

**Scan config keys:**
- `lateral_sites`: explicit fractional sites
- `floor_reference = "local"`: place relative to the local surface, not the top atom
- `quench_steps`, `min_steps`, `reference_min_steps`: cost and convergence
- `hydroxyl_h_lj = [eps, sigma]`: override the slab hydroxyl-H LJ

## 3. Example: diagnosing a force-field artefact with a parameter override

**Symptom.** In `adsorption-scan-nio-groove-hydroxylation.toml`, the NiO 100 % OH groove floor gave E_int ≈ −145 kcal/mol, with phosphonate O 1.36–1.50 Å from slab hydroxyl H. That is shorter than any hydrogen bond.

**Cause.** The slab hydroxyl H (`Hh`) was a bare +0.425 charge with no LJ core (CLAYFF `ho` convention). On a dense OH carpet, ligand O collapses onto it. The same thing happened earlier with the CLAYFF-cation set at under-coordinated rim In (see `adsorption-scan.md`).

**Test.** Rerun only the three collapsed starts, with the same seed and site, under two candidate cores:

```bash
cd src
python -m nio_md_prep.ito adsorption-scan ../studies/ito/examples/hydroxyl-h-core-s040.toml --output ../hh-s040 --workers 1
python -m nio_md_prep.ito adsorption-scan ../studies/ito/examples/hydroxyl-h-core-s107.toml --output ../hh-s107 --workers 1
```

- `s040`: ε 0.046 kcal/mol, σ 0.40 Å. These are the CHARMM-style polar-H values that the NiO workflow already uses for the acidic P–O–H hydrogen (`chemistry.CORRECTIONS`).
- `s107`: ε 0.046, σ 1.07 Å. With geometric mixing this puts the O···H repulsion onset near the typical 1.8 Å H-bond distance.

**Result.** Pending; the table is added when the runs finish.

**Decision rule.** Adopt the softest core that keeps every anchor-O···Hh contact at 1.6 Å or more, then rerun the 50–100 % OH scans.

## 4. Classical energy on a DFT geometry

```bash
python scripts/ito/classical_on_vasp_pose.py inputs/ito/vasp_smoke/ads-me-4pacz-phys --parameter-set uff-cation/clayff-anion
```

This prints the classical E_int at the VASP start geometry and after minimizing the molecule. Compare it with DFT E_ads once the `ads-*`, `slab-oh` and `mol-*` cases have run on QBD. Bidentate cases are refused, because the classical model has only the protonated acid.

## 5. Pilot assembly and analysis

```bash
cd src && python -m nio_md_prep.ito build-pilot ../studies/ito/pilot-me-4pacz-ito-groove-oh025.toml --output ../prepared/ito/test/s11 --packmol-seed 11 --velocity-seed 111
python -m nio_md_prep.ito analyze ../prepared/ito/test/s11 --trajectory ../prepared/ito/test/s11/relax.lammpstrj --last-frames 20
```

- The build writes `validation_report.txt`, which must say PASS, and `assembly_manifest.json` with full provenance.
- On corrugated slabs the analysis reports `groove_floor/groove_wall/plateau` densities and `groove_to_plateau_density_ratio`.
- The MD itself runs on QBD: `scripts/ito/run_pilot_array.sbatch` and `submit_hydroxylation_series.sh`.
