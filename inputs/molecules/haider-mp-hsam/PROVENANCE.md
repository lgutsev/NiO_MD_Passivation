# haider-mp-hsam

MP-HSAM (4PABCz): (4-(9′-phenyl-9H,9′H-[3,3′-bicarbazol]-9-yl)butyl)phosphonic acid, C34H29N2O3P, 69 atoms,
neutral, fully protonated, one phosphonic-acid anchor.

- **Submitted:** `haider-mp-hsam.sdf` from `lgutsev/PVSK_Collab_Desk`,
  `projects/P23_haider_hsam/ligands/ligpargen_handoff/` (RDKit ETKDGv3 + MMFF94, seed 20261019).
- **LigParGen settings:** LigParGen 2.1, 1.14*CM1A-LBCC, net charge 0, 0 optimization rounds (2026-10-20).
- **`ligpargen.lmp`:** byte copy of the returned `LIG.lammps.lmp`, sha256
  `a653f23d21aa7334113b81878c9ebb53ea066aa898915c151886aa7c56cfc234`.
- **Validation:** done by the desk's `scripts/import_ligpargen.py`; see `ligands/accepted/haider-mp-hsam/import_report.json`.
  - graph isomorphic to the submitted structure, with atom order unchanged (RMSD 0.0005 Å);
  - charge sum 0 e;
  - phosphonate roles and 5 corrected torsions, mapped by `chemistry.py`;
  - mass 544.591 g/mol.
- **Charge scheme:** consistent with `dcz-4p` and `me-4pacz`. On the shared N-butyl-phosphonic-acid fragment,
  mean |Δq| is 0.013 e against `dcz-4p`.
- **Overrides:** none. The phosphonate corrections come from connectivity, as for `dcz-4p`.
