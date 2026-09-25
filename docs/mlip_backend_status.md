# MLIP backend status and evidence (feat/mlip-hardening)

What has actually been exercised, on what, and what is still blocked. "Toy model" means
`tiny-nio.model`: a real MACE model trained with the real `mace_run_train` for 3 epochs on
SYNTHETIC Morse labels (sha256 `00020ad89117a4d3a9c2b9ead09bb051f40a590521f04fb1ba74bddd2921af95`, kept outside git).
It proves wiring and numerical interoperability only; **it is not a NiO potential, and nothing here
validates any model for NiO chemistry.**

Tested versions (Windows 11, 2026-09): Python 3.14.4 + ASE 3.29.0 + pip LAMMPS 22 Jul 2025 update 4
(ML-IAP, ML-SNAP; no PYTHON/KOKKOS packages); mace-torch 0.3.16 with torch 2.14.0+cpu (Python 3.14) and
torch 2.12.1 cpu (Python 3.12); OpenMM 8.5.2 + openmmml 1.7 (+ openmm-torch 1.5); platforms Reference, CPU,
OpenCL (NVIDIA RTX 5070 laptop GPU).

## Status table

| capability | implemented | synthetic/unit-tested | real backend (toy model or LJ) | real NiO model | notes / blocker |
|---|---|---|---|---|---|
| MACE via ASE (single point, smoke MD) | yes | yes | **yes** (toy model, real NiO(110) slab) | no | needs a trained NiO model |
| MACE via OpenMM-ML (single point, smoke MD) | yes | yes | **yes** (parity below) | no | NPT refused by policy |
| MACE via LAMMPS ML-IAP | yes (deck, export, launch plan) | yes (plan, export naming with real exporter) | **no** | no | needs LAMMPS with ML-IAP + PYTHON (+KOKKOS for GPU and for multi-layer MACE) |
| MACE via LAMMPS pair_style mace (legacy) | yes | yes | no | no | needs a LAMMPS build with ML-MACE; MACE v1.0 drops this artifact |
| LAMMPS-native pair styles (python and executable routes) | yes | yes | **yes** (lj/cut parity below) | no | needs a trained LAMMPS-native MLIP |
| LAMMPS via ASE LAMMPSlib | yes | yes | **yes** (lj/cut) | no | |
| GPU execution | plan + launch options | yes | OpenMM OpenCL only | no | CUDA torch / KOKKOS GPU LAMMPS not available here |
| Per-axis PBC, triclinic rotation back | yes | yes | **yes** (8 geometries, both LAMMPS routes; OpenMM non-reduced cell) | no | |
| Frozen atoms (FixAtoms) | yes (all engines) | yes | **yes** | no | ASE Nose-Hoover + FixAtoms refused (measured freeze) |
| Thermostats/barostats as requested | yes (tables below) | yes | **yes** (LAMMPS all maps; ASE; OpenMM Langevin/NH) | no | unsupported names refused at validate |
| 500-step smoke cap | yes | yes | yes | — | no production MD driver exists |

## Measured evidence

**Minimum MACE acceptance (ASE bridge, toy model, real structure).** `NiO_110_AFM_compromise.POSCAR`
(InterfaceForge NiO(110) slab, 200 atoms, 40 frozen), `pbc = (T, T, F)`, cpu, float64:
E = −1558.73271094 eV, max |F| = 0.969 eV/Å (both finite); per-atom energies sum to E within 4e-12 eV; manifest
status `completed`, declared model sha256 = observed sha256; resolved runtime (device, dtype, head, E0s) recorded.
20-step NVE smoke: 20 steps completed, 5 frames, drift −1.99e-4 eV/atom/ps (conservation test, NVE only).

**MACE ASE vs OpenMM-ML** (toy model, rattled bulk NiO 2×2×2 in the fcc-primitive cell, which is *not* in
OpenMM's reduced form and is rotated in and out; total energy convention; stress not compared: OpenMM has no virial):

| OpenMM platform (read back from the Context) | ΔE/atom (eV) | force RMSE (eV/Å) | max force-vector error (eV/Å) | max component error (eV/Å) |
|---|---|---|---|---|
| Reference | 0.0 | 2.0e-14 | 6.7e-14 | 5.2e-14 |
| CPU | 0.0 | 2.0e-14 | 6.7e-14 | 5.2e-14 |
| OpenCL, Precision=double | 0.0 | 1.4e-13 | 3.5e-13 | 3.1e-13 |
| OpenCL, Precision=single (OpenMM's default on this GPU machine) | 1.1e-7 | 5.0e-6 | 1.8e-5 | 1.5e-5 |

The last row is what an unpinned run gets: OpenMM picks OpenCL and single-precision accumulation. The
provenance records the platform and its properties as read back from the Context, so this is visible; set
`engine.platform` and `engine.platform_precision = "double"` for parity work. These are measurements, not tolerances.

**LAMMPS vs independent LJ reference** (`lj/cut`, both runtimes, 8 geometries: cubic, triclinic, rotated
triclinic, slabs normal to x/y/z with per-axis PBC, tilted slab, cluster): energy and forces agree to < 1e-10,
stress to < 1e-8 relative; LAMMPSlib vs native route: ΔE 1.6e-14 eV, stress 1.6e-10, per-atom 3e-16.
**Metal vs real units** (same physical LJ system): single point ΔE 5e-15 eV, ΔF 4e-16 eV/Å, stress 3.3e-8 relative;
20-step NVT (Langevin, Nose-Hoover) and NPT (MTK, Berendsen+CSVR) trajectories agree to < 1e-8 Å in positions and
< 6e-10 Å in cell against cell changes of 5e-3 Å — real-units NPT and stress are therefore supported. (A deliberate
bar-read-as-atm slip moves the Berendsen cell by 1.9e-6 Å, so the test would catch it.)
**Thermostats (LAMMPS)**: NVE, Nose-Hoover, Langevin, Berendsen and CSVR all rendered and run; conserved-energy
excursion ≈ 0.004 eV/atom vs total-energy change 0.07–0.09 eV/atom under a thermostat (reported as descriptive, not
as conservation); NVE excursion 4.3e-6 / 1.07e-6 / 3.2e-7 eV/atom at 1 / 0.5 / 0.25 fs (second order).

## Effective integrators per engine

| request | ASE | LAMMPS | OpenMM |
|---|---|---|---|
| NVE | VelocityVerlet (thermostat refused) | fix nve | VerletIntegrator |
| NVT langevin (default) | Langevin (fixcm=False; FixCom when no frozen atoms) | fix nve + fix langevin (zero yes, tally yes) | LangevinMiddleIntegrator |
| NVT nose-hoover | NoseHooverChainNVT (refused with frozen atoms) | fix nvt | NoseHooverIntegrator |
| NVT berendsen | NVTBerendsen | fix nve + fix temp/berendsen | refused |
| NVT csvr | Bussi | fix nve + fix temp/csvr | refused |
| NPT mtk | IsotropicMTKNPT / MaskedMTKNPT | fix npt (geometry-aware coupling; in-plane for slabs) | refused (policy) |
| NPT berendsen | NPTBerendsen / Inhomogeneous_NPTBerendsen | fix press/berendsen (not triclinic) | refused |
| NPT parrinello-rahman | refused | refused | refused |

NPT with frozen atoms is refused on every engine; isotropic/anisotropic NPT on a cell with a vacuum gap is refused.
The manifest records the integrator that actually ran with its resolved parameters.

## Blocked items and how to resume

1. **MACE → LAMMPS/ML-IAP execution and ASE↔ML-IAP parity** — needs a LAMMPS Python module built with ML-IAP +
   PYTHON (+ KOKKOS; multi-layer MACE requires KOKKOS even on CPU). On a LONI GPU node with such an environment:
   ```bash
   sbatch -A <allocation> -p <gpu partition> --export=ALL,PARITY_ENV=<conda env>,MACE_MODEL=<model.model>,\
   STRUCTURE=<POSCAR>,ELEMENTS=Ni,O,C,H,N,P,PBC=T,T,T scripts/mlip_mace_ase_lammps_parity.sbatch
   ```
   It exports the model (`mace_create_lammps_model --format mliap`), evaluates both routes and writes
   `parity_report.json` and a CANDIDATE tolerance file for human review (`src/nio_md_prep/mlip/tolerances.py`).
2. **GPU tests** — need CUDA torch (and a KOKKOS GPU LAMMPS for the ML-IAP route); run
   `NIO_MD_TEST_MACE_DEVICE=cuda pytest -m "mace or gpu" tests/test_mlip_backends.py`.
3. **Real NiO model evaluation** — needs a trained, reviewed NiO/phosphonate model; nothing in this branch
   establishes accuracy.
4. **Cross-engine tolerance reference** — `tests/data/mlip_cross_engine_tolerances.json` is not committed; a
   reviewer promotes a measured candidate by setting `status: accepted` and `reviewed_by`.
