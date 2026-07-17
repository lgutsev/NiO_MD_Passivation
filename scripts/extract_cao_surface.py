"""Deterministically extract molecule-ID-0 NiO from the Cao2025 mixed-SAM deck."""
from pathlib import Path
import hashlib
from nio_md_prep.lammps import DataFile, Record, parse, write

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"examples/Sigma_Big_NiO_110_Multi_Me4PACz_4x_180-360_swall/topology_output.lmp"
DEST=ROOT/"inputs/surfaces/corrugated-nio-110"

source=parse(SOURCE)
selected=[r for r in source.sections["Atoms"] if int(r.fields[1])==0]
old_types=sorted({int(r.fields[2]) for r in selected})
if len(selected)!=21060 or len(old_types)!=2:
    raise SystemExit(f"unexpected authoritative surface inventory: {len(selected)} atoms, types {old_types}")
type_map={old_types[0]:1,old_types[1]:2}
atoms=[]
for new_id,row in enumerate(selected,1):
    fields=row.fields.copy(); fields[0]=str(new_id); fields[1]="0"; fields[2]=str(type_map[int(fields[2])])
    atoms.append(Record(fields,row.comment))
data=DataFile("Cao2025 corrugated NiO surface extracted from successful mixed-SAM topology",{
    "Masses":[Record(["1","58.6934"],"Ni"),Record(["2","15.999"],"O")],
    "Atoms":atoms,
},{"x":source.bounds["x"],"y":source.bounds["y"],"z":source.bounds["z"]})
write(data,DEST/"surface.lmp")
with (DEST/"surface.xyz").open("w",encoding="utf-8") as out:
    out.write(f"{len(atoms)}\nCao2025 corrugated NiO; source SHA256 {hashlib.sha256(SOURCE.read_bytes()).hexdigest()}\n")
    for row in atoms:
        symbol="Ni" if int(row.fields[2])==1 else "O"
        out.write(f"{symbol} {row.fields[4]} {row.fields[5]} {row.fields[6]}\n")
# These are provenance-hashed data assets; force canonical LF bytes on Windows.
for path in (DEST/"surface.lmp",DEST/"surface.xyz"):
    path.write_bytes(path.read_bytes().replace(b"\r\n",b"\n"))
print(hashlib.sha256((DEST/"surface.lmp").read_bytes()).hexdigest())
