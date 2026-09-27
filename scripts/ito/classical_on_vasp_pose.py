"""Classical E_int of a VASP adsorption case geometry (rigid slab, same decomposition as checks.py).

    python scripts/ito/classical_on_vasp_pose.py inputs/ito/vasp_smoke/ads-me-4pacz-phys [--parameter-set NAME]

Reads geometry.extxyz (slab labels + LigParGen ids), rebuilds the classical system with the
ITO charges/LJ, and prints E_int at the DFT start geometry and after minimizing the molecule
with the slab frozen.  The same numbers from DFT come from the ads/slab/mol VASP cases.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from nio_md_prep.config import molecule_manifest  # noqa: E402
from nio_md_prep.ito import checks  # noqa: E402
from nio_md_prep.ito.forcefield import PARAMETER_SETS, surface_pair_lines  # noqa: E402
from nio_md_prep.ito.substrate import Slab, charges, write_lammps  # noqa: E402
from nio_md_prep.lammps import parse, write  # noqa: E402


def read_extxyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    n = int(lines[0]); lat = [float(x) for x in re.search(r'Lattice="([^"]+)"', lines[1]).group(1).split()]
    rows = [l.split() for l in lines[2:2 + n]]
    return (lat[0], lat[4], lat[8]), rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", type=Path)
    ap.add_argument("--parameter-set", default="uff-cation/clayff-anion", choices=sorted(PARAMETER_SETS))
    a = ap.parse_args(argv)
    man = json.loads((a.case / "case_manifest.json").read_text(encoding="utf-8"))
    if "bidentate" in man["case"]:
        raise SystemExit("bidentate cases carry a deprotonated anchor; the classical model has only the protonated acid")
    (lx, ly, lz), rows = read_extxyz(a.case / "geometry.extxyz")
    slug = next(s for s in ("me-4pacz", "meo-2pacz") if s in man["case"])
    slab_rows = [r for r in rows if r[4] != "mol"]
    mol_rows = sorted((r for r in rows if r[4] == "mol"), key=lambda r: int(r[5]))
    slab = Slab([r[4] for r in slab_rows], np.array([[float(x) for x in r[1:4]] for r in slab_rows]), (lx, ly))
    q = charges(slab, 0.525, 0.425)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        zlo, zhi = -3.0, lz + 3.0
        type_ids = write_lammps(slab, q, tmp / "slab.lmp", (zlo, zhi), "vasp pose slab")
        slab_data = parse(tmp / "slab.lmp")
        folder, mm = molecule_manifest(slug)
        mol = parse(folder / mm["files"]["ligpargen"])
        xyz = np.array([[float(x) for x in r[1:4]] for r in mol_rows])
        system, o2n, lines = checks._system(slab_data, mol, xyz, zhi)
        ff = lines + surface_pair_lines({lab: o2n[old] for lab, old in type_ids.items()}, a.parameter_set)
        path = tmp / "pose.data"; write(system, path)
        n, ns = mol.count("Atoms"), slab_data.count("Atoms")
        mesh = [max(8, int(round(lx / 1.2))), max(8, int(round(ly / 1.2))), max(16, int(round(3 * (zhi - zlo) / 1.4)))]
        ks = checks.KSPACE.format(g=0.30, mx=mesh[0], my=mesh[1], mz=mesh[2])
        e_slab = checks._energy_run(path, ff, ks, n, ns, False, delete="mol")[0]
        e0 = checks._energy_run(path, ff, ks, n, ns, False)[0]
        m0 = checks._energy_run(path, ff, ks, n, ns, False, delete="slab")[0]
        e1, xmin = checks._energy_run(path, ff, ks, n, ns, True, min_steps=3000)
        p1 = checks._with_coords(path, tmp / "min.data", xmin)
        m1 = checks._energy_run(p1, ff, ks, n, ns, False, delete="slab")[0]
    out = {"case": man["case"], "parameter_set": a.parameter_set,
           "e_int_at_dft_start_kcal_mol": e0 - e_slab - m0, "e_int_after_molecule_minimization_kcal_mol": e1 - e_slab - m1,
           "rmsd_molecule_after_minimization_angstrom": float(np.sqrt(np.mean(np.sum((xmin - xyz) ** 2, axis=1))))}
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
