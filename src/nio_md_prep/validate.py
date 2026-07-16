from __future__ import annotations
from pathlib import Path
import shutil, subprocess
from .lammps import TOPOLOGY, charge, parse

def validate(folder: Path, packmol_ran: bool|None=None, primary_final: Path|None=None) -> bool:
    errors=[]; warnings=[]; data=parse(folder/"topology_output.lmp")
    atoms=data.sections["Atoms"]; ids=[int(r.fields[0]) for r in atoms]
    if len(ids)!=len(set(ids)) or set(ids)!=set(range(1,len(ids)+1)): errors.append("atom IDs are not unique and contiguous")
    atom_set=set(ids)
    for sec,nrefs in TOPOLOGY.items():
        rows=data.sections.get(sec,[]); rowids=[int(r.fields[0]) for r in rows]
        if len(rowids)!=len(set(rowids)): errors.append(f"duplicate {sec} IDs")
        for r in rows:
            if any(int(x) not in atom_set for x in r.fields[2:2+nrefs]): errors.append(f"invalid atom reference in {sec}")
    mass_types={int(r.fields[0]) for r in data.sections["Masses"]}; used={int(r.fields[2]) for r in atoms}
    if used-mass_types: errors.append(f"used atom types lack masses: {sorted(used-mass_types)}")
    for axis,i in zip("xyz",range(3)):
        lo,hi=data.bounds[axis]
        if any(not(lo<=float(r.fields[4+i])<=hi) for r in atoms): errors.append(f"coordinates outside {axis} bounds")
    # Surface is recognized by fixed charges from the validated asset.
    nio=[r for r in atoms if abs(abs(float(r.fields[3]))-2.0)<1e-8]
    if not nio or sum(float(r.fields[3]) for r in nio)!=0: errors.append("Ni/O fixed-charge population missing or unbalanced")
    if packmol_ran is False: warnings.append("Packmol was not found or not run; packed.xyz ordering/contact validation skipped. Run: packmol < packmol.inp")
    if primary_final and not primary_final.exists(): errors.append("sequential stage 2 did not receive the stage-1 final data")
    lmp=shutil.which("lmp") or shutil.which("lammps")
    if not lmp: warnings.append("LAMMPS executable not found; zero-step dry run skipped")
    report=["VALID" if not errors else "INVALID",f"atoms: {len(atoms)}",f"total charge: {charge(data):.8f}"]+[f"ERROR: {x}" for x in errors]+[f"WARNING: {x}" for x in warnings]
    (folder/"validation_report.txt").write_text("\n".join(report)+"\n",encoding="utf-8")
    if errors: raise ValueError("; ".join(errors))
    return True
