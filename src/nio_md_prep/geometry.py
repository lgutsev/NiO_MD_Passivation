from __future__ import annotations
from pathlib import Path
from .lammps import DataFile, atom_coordinates

ELEMENTS = {1.008:"H", 12.011:"C", 14.007:"N", 15.999:"O", 18.998:"F", 28.085:"Si", 30.974:"P", 32.060:"S", 35.450:"Cl", 58.6934:"Ni"}

def elements(data: DataFile) -> list[str]:
    masses = {int(r.fields[0]): float(r.fields[1]) for r in data.sections["Masses"]}
    by_type = {}
    for tid, mass in masses.items():
        hits = [symbol for expected, symbol in ELEMENTS.items() if abs(expected-mass) <= 0.02]
        if len(hits) != 1:
            raise ValueError(f"ambiguous element for atom type {tid} with mass {mass}; provide geometry.xyz")
        by_type[tid] = hits[0]
    return [by_type[int(r.fields[2])] for r in data.sections["Atoms"]]

def write_xyz(data: DataFile, path: Path, comment: str) -> None:
    symbols = elements(data); xyz = atom_coordinates(data)
    path.write_text(f"{len(xyz)}\n{comment}\n" + "".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s,(x,y,z) in zip(symbols,xyz)), encoding="utf-8")

def read_xyz(path: Path) -> tuple[list[str], list[tuple[float,float,float]]]:
    lines = path.read_text(encoding="utf-8").splitlines(); n = int(lines[0]); rows = lines[2:]
    if len(rows) != n: raise ValueError(f"{path}: XYZ count mismatch")
    symbols=[]; xyz=[]
    for row in rows:
        p=row.split(); symbols.append(p[0]); xyz.append(tuple(map(float,p[1:4])))
    return symbols, xyz
