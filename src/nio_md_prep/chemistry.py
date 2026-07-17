"""Chemistry-derived inventory and Cao phosphonate corrections."""
from __future__ import annotations
from .lammps import DataFile

CORRECTIONS = {
    "P": (0.200, 3.742),
    "P=O": (0.210, 2.960),
    "P-OH": (0.211, 3.066),
    "P-O-H": (0.046, 0.400),
}

def molecular_weight(data: DataFile) -> float:
    masses = {int(r.fields[0]): float(r.fields[1]) for r in data.sections["Masses"]}
    return sum(masses[int(a.fields[2])] for a in data.sections["Atoms"])

def _element(data: DataFile, atom_id: int) -> str:
    atom = next(a for a in data.sections["Atoms"] if int(a.fields[0]) == atom_id)
    mass = next(m for m in data.sections["Masses"] if int(m.fields[0]) == int(atom.fields[2]))
    token = mass.comment.split()[-1] if mass.comment else ""
    if token in {"H","C","N","O","P"}: return token
    value=float(mass.fields[1])
    for symbol,reference in (("H",1.008),("C",12.011),("N",14.007),("O",15.999),("P",30.974)):
        if abs(value-reference)<0.2: return symbol
    raise ValueError(f"cannot infer element for atom {atom_id}, mass {value}")

def phosphonate_roles(data: DataFile) -> dict[int,str]:
    adjacency={int(a.fields[0]):set() for a in data.sections["Atoms"]}
    for b in data.sections.get("Bonds",[]):
        a,c=map(int,b.fields[2:4]); adjacency[a].add(c); adjacency[c].add(a)
    roles={}
    for p in (i for i in adjacency if _element(data,i)=="P"):
        oxygens=[i for i in adjacency[p] if _element(data,i)=="O"]
        carbons=[i for i in adjacency[p] if _element(data,i)=="C"]
        if len(oxygens)!=3 or len(carbons)!=1:
            raise ValueError(f"ambiguous phosphonate at P atom {p}: expected three O and one C neighbors")
        roles[p]="P"
        for o in oxygens:
            hydrogens=[i for i in adjacency[o] if _element(data,i)=="H"]
            if len(hydrogens)==1: roles[o]="P-OH"; roles[hydrogens[0]]="P-O-H"
            elif len(adjacency[o])==1: roles[o]="P=O"
            else: raise ValueError(f"ambiguous phosphonate oxygen atom {o}")
    if not roles: raise ValueError("no phosphonate group identified from element/connectivity")
    return roles

def correction_lines(data: DataFile, atom_type_offset: int) -> tuple[list[str],int]:
    roles=phosphonate_roles(data); by_type={}
    atoms={int(a.fields[0]):a for a in data.sections["Atoms"]}
    for atom_id,role in roles.items():
        typ=int(atoms[atom_id].fields[2])+atom_type_offset
        if typ in by_type and by_type[typ]!=role: raise ValueError(f"atom type {typ} has conflicting phosphonate roles")
        by_type[typ]=role
    lines=[]
    for typ,role in sorted(by_type.items()):
        eps,sigma=CORRECTIONS[role]
        lines.append(f"pair_coeff {typ} {typ} lj/cut/coul/long {eps:.3f} {sigma:.3f} # Cao corrected {role}")
    return lines, sum(1 for role in roles.values() if role=="P")
