"""Lossless-enough parser and deterministic remapper for LigParGen data files."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

SECTIONS = (
    "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "Atoms", "Bonds", "Angles",
    "Dihedrals", "Impropers", "Velocities",
)
TOPOLOGY = {"Bonds": 2, "Angles": 3, "Dihedrals": 4, "Impropers": 4}
COUNT_NAMES = {"Atoms": "atoms", "Bonds": "bonds", "Angles": "angles", "Dihedrals": "dihedrals", "Impropers": "impropers"}
TYPE_NAMES = {"Masses": "atom", "Pair Coeffs": "atom", "Bond Coeffs": "bond", "Angle Coeffs": "angle", "Dihedral Coeffs": "dihedral", "Improper Coeffs": "improper", "Atoms": "atom", "Bonds": "bond", "Angles": "angle", "Dihedrals": "dihedral", "Impropers": "improper"}

@dataclass
class Record:
    fields: list[str]
    comment: str = ""

    def render(self) -> str:
        return " ".join(self.fields) + ((" # " + self.comment) if self.comment else "")

@dataclass
class DataFile:
    title: str = "LAMMPS data file"
    sections: dict[str, list[Record]] = field(default_factory=dict)
    bounds: dict[str, tuple[float, float]] = field(default_factory=dict)
    declared_types: dict[str, int] = field(default_factory=dict)
    declared_counts: dict[str, int] = field(default_factory=dict)

    def count(self, section: str) -> int:
        return len(self.sections.get(section, []))

    def type_count(self, category: str) -> int:
        relevant = "Masses" if category == "atom" else category.title() + " Coeffs"
        ids = [int(r.fields[0]) for r in self.sections.get(relevant, [])]
        if category == "atom":
            ids += [int(r.fields[2]) for r in self.sections.get("Atoms", [])]
        else:
            topology = category.title() + "s"
            ids += [int(r.fields[1]) for r in self.sections.get(topology, [])]
        return max(max(ids, default=0), self.declared_types.get(category, 0))

def _section_name(line: str) -> str | None:
    bare = line.split("#", 1)[0].strip()
    return bare if bare in SECTIONS else None

def parse(path: Path) -> DataFile:
    if not path.exists():
        raise FileNotFoundError(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    data = DataFile(title=lines[0].strip() if lines else "LAMMPS data file")
    current = None
    for line in lines[1:]:
        sec = _section_name(line)
        if sec:
            current = sec
            data.sections.setdefault(sec, [])
            continue
        stripped = line.strip()
        if not stripped:
            continue
        m = re.match(r"^\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([xyz])lo\s+\3hi", line)
        if m:
            data.bounds[m.group(3)] = (float(m.group(1)), float(m.group(2)))
            continue
        tm = re.match(r"^\s*(\d+)\s+(atom|bond|angle|dihedral|improper)\s+types\s*$", line)
        if tm:
            data.declared_types[tm.group(2)] = int(tm.group(1))
            continue
        cm = re.match(r"^\s*(\d+)\s+(atoms|bonds|angles|dihedrals|impropers)\s*$", line)
        if cm:
            data.declared_counts[cm.group(2).title()] = int(cm.group(1))
            continue
        if current and re.match(r"^\s*\d", line):
            payload, _, comment = line.partition("#")
            data.sections[current].append(Record(payload.split(), comment.strip()))
    if not data.sections.get("Masses") or not data.sections.get("Atoms"):
        raise ValueError(f"{path}: required Masses or Atoms section is missing")
    for name, nrefs in TOPOLOGY.items():
        for r in data.sections.get(name, []):
            if len(r.fields) < 2 + nrefs:
                raise ValueError(f"{path}: malformed {name} record: {r.render()}")
    return data

def write(data: DataFile, path: Path) -> None:
    order = ["Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers"]
    out = [data.title, ""]
    for sec, label in COUNT_NAMES.items():
        if data.count(sec) or sec == "Atoms": out.append(f"{data.count(sec)} {label}")
    out.append("")
    for cat in ("atom", "bond", "angle", "dihedral", "improper"):
        n = data.type_count(cat)
        if n: out.append(f"{n} {cat} types")
    out.append("")
    for axis in "xyz":
        lo, hi = data.bounds.get(axis, (0.0, 1.0))
        out.append(f"{lo:.8f} {hi:.8f} {axis}lo {axis}hi")
    for sec in order:
        records = data.sections.get(sec, [])
        if records:
            suffix = " # full" if sec == "Atoms" else ""
            out += ["", sec + suffix, ""] + [r.render() for r in records]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")

def atom_coordinates(data: DataFile) -> list[tuple[float, float, float]]:
    return [(float(r.fields[4]), float(r.fields[5]), float(r.fields[6])) for r in data.sections["Atoms"]]

def charge(data: DataFile) -> float:
    return sum(float(r.fields[3]) for r in data.sections["Atoms"])

def replicate(template: DataFile, count: int, type_offsets: dict[str, int], id_offsets: dict[str, int], molecule_offset: int, coordinates: list[tuple[float,float,float]]) -> tuple[dict[str,list[Record]], dict[str,int]]:
    """Replicate a complete template, remapping all IDs and topology references."""
    if len(coordinates) != count * template.count("Atoms"):
        raise ValueError("packed coordinate count does not match molecule replicas")
    result: dict[str, list[Record]] = {s: [] for s in SECTIONS}
    # Coefficients and masses are copied once with per-category offsets.
    for sec in ("Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs"):
        cat = TYPE_NAMES[sec]
        for r in template.sections.get(sec, []):
            f = r.fields.copy(); f[0] = str(int(f[0]) + type_offsets.get(cat, 0))
            result[sec].append(Record(f, r.comment))
    atom_n = template.count("Atoms")
    local_ids = [int(r.fields[0]) for r in template.sections["Atoms"]]
    for replica in range(count):
        atom_map = {old: id_offsets["atom"] + replica * atom_n + i + 1 for i, old in enumerate(local_ids)}
        for i, r in enumerate(template.sections["Atoms"]):
            f = r.fields.copy()
            f[0] = str(atom_map[int(f[0])]); f[1] = str(molecule_offset + replica + 1)
            f[2] = str(int(f[2]) + type_offsets["atom"])
            xyz = coordinates[replica * atom_n + i]; f[4:7] = [f"{v:.8f}" for v in xyz]
            result["Atoms"].append(Record(f, r.comment))
        for sec, nrefs in TOPOLOGY.items():
            base = id_offsets[sec.lower()[:-1]] + replica * template.count(sec)
            cat = TYPE_NAMES[sec]
            for i, r in enumerate(template.sections.get(sec, [])):
                f = r.fields.copy(); f[0] = str(base + i + 1); f[1] = str(int(f[1]) + type_offsets[cat])
                f[2:2+nrefs] = [str(atom_map[int(x)]) for x in f[2:2+nrefs]]
                result[sec].append(Record(f, r.comment))
    increments = {"atom": count*atom_n, "molecule": count}
    for sec in TOPOLOGY: increments[sec.lower()[:-1]] = count*template.count(sec)
    return result, increments
