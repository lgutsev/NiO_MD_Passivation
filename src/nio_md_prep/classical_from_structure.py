"""Type a given NiO slab + one phosphonic-acid SAM structure (e.g. a VASP POSCAR or AIMD frame) with the project's
classical force field, so the classical model can be compared with DFT and MLIPs on identical coordinates.

The force field is exactly the one build() writes: LigParGen OPLS-AA for the ligand with the phosphonate corrections
of chemistry.correction_lines, Buckingham Ni/O (+2/-2) for the oxide, and the pseudo-LJ ligand/surface cross terms of
build.surface_force_field_lines.

The classical model has a fixed, fully protonated phosphonic acid and no surface hydroxyls. A DFT structure whose
anchor has given protons to the oxide is therefore re-protonated when `reprotonate` is true: each missing acid H is
taken from a surface hydroxyl and placed 0.97 A from a P-bound O that carries no H (the O with the longest P-O bond
first, never the shortest P=O), pointing away from P and the slab. Structures with surface H left over after that
(hydroxylated surfaces) are refused; the classical model has no parameters for them. Every change is reported.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

from .agglomeration import (
    MoleculeTemplate, _infer_molecular_adjacency, _labeled_graph_mapping, _COVALENT_RADII_ANGSTROM,
)
from .build import HEADER, _coeff_lines, _merge, _surface, surface_force_field_lines
from .chemistry import correction_lines
from .config import molecule_dir, molecule_manifest
from .lammps import DataFile, Record, charge, parse, replicate, write

MASS = {"Ni": 58.6934, "O": 15.999}
OH_LENGTH = 0.97


def _mic(d, cell):
    return [d[k] - cell[k] * round(d[k] / cell[k]) if k < 2 else d[k] for k in range(3)]


def _dist(a, b, cell):
    return math.sqrt(sum(x * x for x in _mic([a[i] - b[i] for i in range(3)], cell)))


def _template(slug: str):
    _, manifest = molecule_manifest(slug)
    source = molecule_dir(slug) / manifest["files"]["ligpargen"]
    data = parse(source)
    from .geometry import elements  # mass -> element, as the rest of the package
    rows = sorted(data.sections["Atoms"], key=lambda r: int(r.fields[0]))
    if [int(r.fields[0]) for r in rows] != [int(r.fields[0]) for r in data.sections["Atoms"]]:
        raise ValueError(f"{slug}: ligpargen.lmp Atoms are not in ID order")
    symbols = tuple(elements(data))
    pos = {int(r.fields[0]): i for i, r in enumerate(rows)}
    adjacency = [set() for _ in symbols]
    for bond in data.sections.get("Bonds", []):
        a, b = pos[int(bond.fields[2])], pos[int(bond.fields[3])]
        adjacency[a].add(b)
        adjacency[b].add(a)
    return data, manifest, symbols, adjacency


def _split(symbols, coords, cell):
    """Molecule = atoms bonded (covalently, minimum image) to a C or P network; slab = Ni and the remaining O;
    surface H = H bonded to slab O."""
    n = len(symbols)
    heavy_org = [i for i, s in enumerate(symbols) if s in ("C", "P", "N")]
    mol = set(heavy_org)
    changed = True
    while changed:
        changed = False
        for i in range(n):
            if i in mol or symbols[i] == "Ni":
                continue
            for j in list(mol):
                cut = _COVALENT_RADII_ANGSTROM[symbols[i]] + _COVALENT_RADII_ANGSTROM[symbols[j]] + 0.3
                if symbols[i] == "H" and symbols[j] not in ("C", "N", "O"):
                    continue
                if symbols[i] == "O" and symbols[j] not in ("C", "P"):
                    continue
                if _dist(coords[i], coords[j], cell) <= cut:
                    mol.add(i)
                    changed = True
                    break
    slab = [i for i in range(n) if i not in mol and symbols[i] in ("Ni", "O")]
    surf_h = [i for i in range(n) if i not in mol and symbols[i] == "H"]
    return sorted(mol), slab, surf_h


def _unwrap(idx, coords, cell, symbols):
    """Unwrap molecule atoms through their bonds (minimum image in x, y)."""
    out = {idx[0]: list(coords[idx[0]])}
    todo = [idx[0]]
    rest = set(idx[1:])
    while todo:
        a = todo.pop()
        for b in list(rest):
            cut = _COVALENT_RADII_ANGSTROM[symbols[a]] + _COVALENT_RADII_ANGSTROM[symbols[b]] + 0.4
            d = _mic([coords[b][k] - coords[a][k] for k in range(3)], cell)
            if math.sqrt(sum(x * x for x in d)) <= cut:
                out[b] = [out[a][k] + d[k] for k in range(3)]
                rest.discard(b)
                todo.append(b)
    if rest:
        raise ValueError(f"molecule is not one bonded fragment ({len(rest)} atoms unattached)")
    return [out[i] for i in idx]


def build_from_structure(slug: str, symbols: list[str], coords: list[tuple[float, float, float]],
                         cell: tuple[float, float, float], output: Path, reprotonate: bool = True,
                         source: str = "") -> dict:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    cell = tuple(float(c) for c in cell)
    data, manifest, t_symbols, t_adj = _template(slug)
    mol, slab, surf_h = _split(symbols, coords, cell)
    m_sym = [symbols[i] for i in mol]
    m_xyz = _unwrap(mol, coords, cell, symbols)
    report = {"slug": slug, "source": source, "molecule_atoms": len(mol), "slab_atoms": len(slab),
              "surface_h_found": len(surf_h), "moved_h": []}
    need = Counter(t_symbols)["H"] - Counter(m_sym)["H"]
    if need < 0:
        raise ValueError(f"molecule has {-need} more H than the LigParGen template")
    if need:
        if not reprotonate:
            raise ValueError(f"anchor is missing {need} acid H; pass reprotonate=True or use another structure")
        # P-bound O without H, longest P-O first
        p_idx = [k for k, s in enumerate(m_sym) if s == "P"]
        cand = []
        for k, s in enumerate(m_sym):
            if s != "O":
                continue
            dp = [(math.dist(m_xyz[k], m_xyz[p]), p) for p in p_idx]
            d, p = min(dp)
            if d > 1.75:
                continue
            has_h = any(m_sym[h] == "H" and math.dist(m_xyz[k], m_xyz[h]) < 1.2 for h in range(len(m_sym)))
            if not has_h:
                cand.append((-d, k, p))
        # never the shortest P-O of each P (P=O)
        shortest = {}
        for k, s in enumerate(m_sym):
            if s == "O":
                d, p = min((math.dist(m_xyz[k], m_xyz[q]), q) for q in p_idx)
                if d <= 1.75 and (p not in shortest or d < shortest[p][0]):
                    shortest[p] = (d, k)
        cand = [c for c in sorted(cand) if c[1] not in {v[1] for v in shortest.values()}]
        if len(cand) < need or len(surf_h) < need:
            raise ValueError(f"cannot re-protonate: need {need} H, {len(cand)} acid O sites, {len(surf_h)} surface H")
        for _, k, p in cand[:need]:
            o, pp = m_xyz[k], m_xyz[p]
            v = [o[i] - pp[i] for i in range(3)]
            v[2] += 1.0  # tilt away from the slab
            norm = math.sqrt(sum(x * x for x in v))
            h_pos = [o[i] + OH_LENGTH * v[i] / norm for i in range(3)]
            # take the surface H nearest to this O
            hsrc = min(surf_h, key=lambda h: _dist(coords[h], o, cell))
            surf_h.remove(hsrc)
            report["moved_h"].append({"from_surface_index": hsrc, "to_acid_o_of_P": True,
                                      "displacement_A": round(_dist(coords[hsrc], h_pos, cell), 3)})
            m_sym.append("H")
            m_xyz.append(h_pos)
    if surf_h:
        raise ValueError(f"{len(surf_h)} surface H remain (hydroxylated surface): no classical parameters for them")
    # map to LigParGen order
    mapping = None
    for tol in (0.20, 0.30, 0.40, 0.50):
        adj = _infer_molecular_adjacency(m_sym, [tuple(x) for x in m_xyz], tol)
        mapping = _labeled_graph_mapping(t_symbols, t_adj, m_sym, adj)
        if mapping is not None:
            report["graph_tolerance_A"] = tol
            break
    if mapping is None:
        raise ValueError("molecule graph does not match the LigParGen topology")
    ordered = [tuple(m_xyz[i]) for i in mapping]
    result = DataFile("nio-md-prep classical_from_structure: " + (source or slug), {}, {})
    ids = {"atom": 0, "bond": 0, "angle": 0, "dihedral": 0, "improper": 0}
    types = {k: 0 for k in ids}
    pieces, inc = replicate(data, 1, types, ids, 0, ordered)
    _merge(result, pieces)
    lines, anchors, corrected = correction_lines(data, types["atom"], types["dihedral"])
    expected = int(manifest.get("molecule", {}).get("phosphonic_acid_anchors", 0))
    if anchors != expected:
        raise ValueError(f"{slug}: {anchors} phosphonate anchors found, manifest expects {expected}")
    for k in ids:
        ids[k] += inc.get(k, 0)
        types[k] += data.type_count(k)
    surface = DataFile("slab", {"Masses": [Record(["1", f"{MASS['Ni']}"], "Ni"), Record(["2", f"{MASS['O']}"], "O")],
                                "Atoms": []}, {})
    for n, i in enumerate(slab):
        q = "2.000" if symbols[i] == "Ni" else "-2.000"
        x, y, z = coords[i]
        surface.sections["Atoms"].append(Record([str(n + 1), "0", "1" if symbols[i] == "Ni" else "2", q,
                                                 f"{x % cell[0]:.6f}", f"{y % cell[1]:.6f}", f"{z:.6f}"]))
    rows, _ = _surface(surface, types, ids, 0)
    _merge(result, rows)
    surface_types = {types["atom"] + 1, types["atom"] + 2}
    zs = [float(r.fields[6]) for r in result.sections["Atoms"]]
    result.bounds = {"x": (0.0, cell[0]), "y": (0.0, cell[1]), "z": (min(zs) - 1.0, max(min(zs) - 1.0 + cell[2], max(zs) + 1.0))}
    write(result, output / "structure.data")
    ff = HEADER + "\n" + "\n".join(_coeff_lines(result) + lines) + "\n"
    ff += surface_force_field_lines(result, surface_types, lines)
    (output / "force_field.lmp").write_text(ff, encoding="utf-8")
    report.update({"atoms_total": len(result.sections["Atoms"]), "total_charge": round(charge(result), 6),
                   "phosphonate_anchors": anchors, "bounds": result.bounds,
                   "molecule_atom_ids": [1, inc["atom"]], "surface_types": sorted(surface_types)})
    (output / "classical_build.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return report
