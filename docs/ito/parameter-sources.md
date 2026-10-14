# ITO extension: parameter sources

This file lists every number the ITO model uses and where it comes from.

Verification legend:
- **V**: the value was checked against the cited source or an authoritative tabulation during this work.
- **C**: the value is transcribed from a standard reference, but the table itself was not re-read here.
- **D**: design choice made here, and why.

## Crystal structure

| Quantity | Value | Source | Status |
|---|---|---|---|
| In2O3 bixbyite, Ia-3 (206), a | 10.117 Å | Marezio, *Acta Cryst.* 20, 723 (1966), doi:10.1107/S0365110X66001749, as deposited in COD 2310009 | V (COD CIF read 2026-09-26) |
| In(8b) | (1/4, 1/4, 1/4) | same | V |
| In(24d) x | 0.4663 | same | V |
| O(48e) | (0.3912, 0.1558, 0.3796) | same | V. An earlier draft used (0.3905, 0.1529, 0.3832) transcribed from memory; this was corrected and all slabs rebuilt |
| Internal check | In32O48 per cell; every In sixfold with In–O 2.133–2.247 Å; 16 empty 16c anion sites | computed (`substrate._bulk_with_vacancies`) | V |
| (111) setting | a[1-10], a[11-2], a[111] → 14.3076 × 24.7815 × 17.5232 Å | computed | V |
| Trilayer | In32O48, d(222) = a/(2√3) = 2.9205 Å, O6-O18-In32-O18-O6 | computed | V |
| Top face (bare) | 12 In5c and 12 O3c per primitive (111) cell (6.77 nm⁻² each) | computed; consistent with the In2O3(111) surface literature | V (count) / literature cross-check in the matrix |

## Surface chemistry

| Quantity | Value | Basis | Status |
|---|---|---|---|
| Hydroxylation | 3 dissociated H2O per primitive (111) cell → 3.38 OH/nm² (counting both O_wH and O_sH) | room-temperature saturation on In2O3(111): 10.1021/acsnano.7b06387 (abstract), 10.1021/acsnano.2c09115 (full text) | V (value). Placement is geometric, not the DFT site pattern. This is a UHV single-crystal value; plasma-treated ITO likely carries more OH |
| In–OH length | 2.10 Å | typical In–O | D (geometric; unrelaxed) |
| O–H length | 0.97 Å | typical hydroxyl | D |
| Sn content | Sn/(In+Sn) = 0.0924 (target 0.0928 = 90:10 wt% In2O3:SnO2) | standard ITO sputter-target composition | D |
| Sn compensation | 2 Sn_In + O_i on empty 16c sites, at least one trilayer below either face | Frank & Köstlin, *Appl. Phys. A* 27, 197 (1982) (defect-cluster model) | C |
| Oxygen vacancies | none | a fixed-charge model cannot hold the electrons a vacancy donates | D |

## Charges (rigid slab)

| Label | q (e) | Rule | Source |
|---|---|---|---|
| In | +1.575 | 0.525 × (+3) | CLAYFF octahedral-cation pattern: Cygan, Liang, Kalinichev, *J. Phys. Chem. B* 108, 1255 (2004) |
| Sn | +2.100 | 0.525 × (+4) | same scaling |
| O (lattice) | −1.050 | 0.525 × (−2) | CLAYFF ob |
| O (hydroxyl) | −0.950 | −0.525 − 0.425 | CLAYFF oh |
| H (hydroxyl) | +0.425 | — | CLAYFF ho |

The scale 0.525 is the one value for which the CLAYFF bulk-O and hydroxyl charges are simultaneously reproduced. It keeps dissociative hydroxylation and (2 Sn + O_i) doping exactly neutral. Status: **D**. It is not fitted to In2O3; DFT charge partitioning (Bader/DDEC) on the VASP smoke slabs would test it.

## Lennard-Jones (slab self terms; cross terms by geometric mixing)

| Set | In | Sn | O / O(H) | H(O) | Source |
|---|---|---|---|---|---|
| `uff-cation/clayff-anion` | ε 0.599 kcal/mol, σ 3.976 Å | ε 0.567, σ 3.913 | ε 0.1554, σ 3.1655 | 0 | UFF (Rappé et al., *JACS* 114, 10024 (1992)); x1 and D1 verified in Open Babel `data/UFF.prm` (In3+3: 4.463 / 0.599; Sn3: 4.392 / 0.567); σ = x1/2^(1/6). O and H from CLAYFF. **V** |
| `clayff-cation/clayff-anion` | ε 1.3298e-6, σ 4.2713 | same | ε 0.1554, σ 3.1655 | 0 | CLAYFF ao (octahedral Al) used as a bare-cation **analogue**, not an In parameterisation. **D** |

The two sets bracket the unknown In/Sn dispersion. Neither set is fitted to ITO adsorption data.

## Ligands (unchanged from the NiO workflow)

| Item | Source |
|---|---|
| OPLS-AA charges and bonded terms | LigParGen (Dodda et al., *Nucleic Acids Res.* 45, W331 (2017)), files in `inputs/molecules/*/ligpargen.lmp` |
| Phosphonate LJ and torsion corrections | Cao-corrected values in `chemistry.CORRECTIONS` (see `docs/project-design.md`) |
| Styles | harmonic / harmonic / opls+charmm / cvff; `special_bonds amber`; LJ 10 Å, Coulomb 8 Å; PPPM 1e-4 with slab 3.0 |

## Analysis radii

| Element | vdW radius (Å) | Source |
|---|---|---|
| In | 1.93 | Bondi (1964) / Mantina et al., *J. Phys. Chem. A* 113, 5806 (2009) |
| Sn | 2.17 | same |
| I | 1.98 | Bondi (1964) |
