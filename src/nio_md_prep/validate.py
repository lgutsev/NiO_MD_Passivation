from __future__ import annotations
from pathlib import Path
import json, math, re, shutil
from .lammps import TOPOLOGY, charge, parse

def _records(data, section):
    return [r.fields for r in data.sections.get(section,[])]

def validate(folder: Path, packmol_ran: bool|None=None, primary_final: Path|None=None) -> bool:
    errors=[]; warnings=[]
    topology=folder/"topology_output.lmp"
    if not topology.exists():
        raise ValueError("packing is incomplete; topology_output.lmp was not produced. Run: packmol < packmol.inp")
    data=parse(topology); manifest=json.loads((folder/"assembly_manifest.json").read_text())
    atoms=data.sections["Atoms"]; ids=[int(r.fields[0]) for r in atoms]; by_id={int(r.fields[0]):r for r in atoms}
    for sec in ("Atoms","Bonds","Angles","Dihedrals","Impropers"):
        declared=data.declared_counts.get(sec)
        if declared is not None and declared!=data.count(sec): errors.append(f"header declares {declared} {sec}, found {data.count(sec)}")
    if len(ids)!=len(set(ids)): errors.append("atom IDs are not unique")
    atom_set=set(ids)
    positive_molecules={int(r.fields[1]) for r in atoms if int(r.fields[1])>0}
    expected=sum(c["count"] for c in manifest.get("components",[]))
    if len(positive_molecules)<expected: errors.append("physical ligand copies do not have unique positive molecule IDs")
    for sec,nrefs in TOPOLOGY.items():
        rows=data.sections.get(sec,[]); rowids=[int(r.fields[0]) for r in rows]
        if len(rowids)!=len(set(rowids)): errors.append(f"duplicate {sec} IDs")
        for r in rows:
            refs=[int(x) for x in r.fields[2:2+nrefs]]
            if any(x not in atom_set for x in refs): errors.append(f"invalid atom reference in {sec}"); continue
            memberships={int(by_id[x].fields[1]) for x in refs}
            if len(memberships)!=1: errors.append(f"{sec} term connects different molecules")
            if 0 in memberships: errors.append(f"surface atom appears in bonded {sec} term")
    mass_types={int(r.fields[0]) for r in data.sections["Masses"]}; used={int(r.fields[2]) for r in atoms}
    if used-mass_types: errors.append(f"used atom types lack masses: {sorted(used-mass_types)}")
    coefficient_sections={"Bonds":"Bond Coeffs","Angles":"Angle Coeffs","Dihedrals":"Dihedral Coeffs","Impropers":"Improper Coeffs"}
    for topology_section,coefficient_section in coefficient_sections.items():
        used_bonded={int(r.fields[1]) for r in data.sections.get(topology_section,[])}
        defined={int(r.fields[0]) for r in data.sections.get(coefficient_section,[])}
        if used_bonded-defined: errors.append(f"used {topology_section} types lack coefficients: {sorted(used_bonded-defined)}")
    for r in atoms:
        values=[float(x) for x in r.fields[3:7]]
        if not all(math.isfinite(x) for x in values): errors.append(f"non-finite atom value at ID {r.fields[0]}")
    for axis,i in zip("xyz",range(3)):
        lo,hi=data.bounds[axis]
        if any(not(lo<=float(r.fields[4+i])<=hi) for r in atoms): errors.append(f"coordinates outside {axis} bounds")
    nio=[r for r in atoms if abs(abs(float(r.fields[3]))-2.0)<1e-8]
    if not nio or abs(sum(float(r.fields[3]) for r in nio))>1e-8: errors.append("Ni/O fixed-charge population missing or unbalanced")
    ff=(folder/"force_field_settings_lammps_with_header.lmp").read_text()
    nio_types=sorted({int(r.fields[2]) for r in nio})
    ligand_types=sorted(used-set(nio_types))
    if len(nio_types)!=2: errors.append("expected exactly two Ni/O atom types")
    else:
        ni,o=nio_types
        for a,b in ((ni,ni),(ni,o),(o,o)):
            if not re.search(rf"^pair_coeff\s+{a}\s+{b}\s+buck/coul/long",ff,re.M): errors.append(f"missing Buckingham surface pair {a}-{b}")
        for surface in nio_types:
            for ligand in ligand_types:
                if not re.search(rf"^pair_coeff\s+{surface}\s+{ligand}\s+lj/cut/coul/long",ff,re.M): errors.append(f"missing LJ ligand-surface pair {surface}-{ligand}")
        if any(re.search(rf"^pair_coeff\s+{a}\s+{b}\s+buck/coul/long",ff,re.M) for a in ligand_types for b in ligand_types): errors.append("ligand-ligand pair accidentally uses Buckingham")
    mapped_types=set()
    for component in manifest.get("components",[]):
        for bounds in component.get("types",{}).values(): mapped_types.update(range(bounds[0],bounds[1]+1))
    if set(ligand_types)-mapped_types: errors.append("topology atom types disagree with assembly mapping")
    if abs(charge(data)-float(manifest["total_charge"]))>1e-6: errors.append("final charge disagrees with component/assembly manifest")
    if primary_final:
        old=parse(primary_final)
        for sec in ("Masses","Pair Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs","Improper Coeffs","Atoms","Bonds","Angles","Dihedrals","Impropers"):
            before=_records(old,sec); after=_records(data,sec)[:len(before)]
            if before!=after: errors.append(f"sequential stage-1 {sec} changed")
    if packmol_ran is False: errors.append("Packmol packing was not completed")
    if not (shutil.which("lmp") or shutil.which("lammps")): warnings.append("LAMMPS executable not found; zero-step dry run skipped")
    report=["VALID" if not errors else "INVALID",f"atoms: {len(atoms)}",f"molecules: {max(int(r.fields[1]) for r in atoms)}",f"total charge: {charge(data):.8f}"]
    if not errors: report.append("Spreadsheet assembly: not required; topology and coordinates were merged programmatically.")
    report += [f"ERROR: {x}" for x in errors]+[f"WARNING: {x}" for x in warnings]
    (folder/"validation_report.txt").write_text("\n".join(report)+"\n",encoding="utf-8")
    if errors: raise ValueError("; ".join(dict.fromkeys(errors)))
    return True
