import math

import pytest

from nio_md_prep.classical_from_structure import build_from_structure
from nio_md_prep.config import molecule_dir
from nio_md_prep.geometry import elements
from nio_md_prep.lammps import parse


def _slab():
    a = 4.186
    b = a / math.sqrt(2)
    d = a / (2 * math.sqrt(2))
    sym, xyz = [], []
    for n in range(5):
        sx, sy = (n % 2) * a / 2, (n % 2) * b / 2
        for i in range(4):
            for j in range(5):
                sym += ["Ni", "O"]
                xyz += [(i * a + sx, j * b + sy, 5 + n * d), (i * a + a / 2 + sx, j * b + sy, 5 + n * d)]
    return sym, xyz, (4 * a, 5 * b, 48.0)


def _molecule(slug, z0):
    data = parse(molecule_dir(slug) / "ligpargen.lmp")
    sym = elements(data)
    xyz = [tuple(float(v) for v in r.fields[4:7]) for r in data.sections["Atoms"]]
    zmin = min(p[2] for p in xyz)
    cx = sum(p[0] for p in xyz) / len(xyz)
    cy = sum(p[1] for p in xyz) / len(xyz)
    return sym, [(p[0] - cx + 8.0, p[1] - cy + 7.0, p[2] - zmin + z0) for p in xyz]


def test_protonated_molecule_on_slab(tmp_path):
    ssym, sxyz, cell = _slab()
    msym, mxyz = _molecule("me-4pacz", 14.0)
    report = build_from_structure("me-4pacz", ssym + msym, sxyz + mxyz, cell, tmp_path / "out", reprotonate=False)
    assert report["atoms_total"] == 245
    assert report["slab_atoms"] == 200 and report["moved_h"] == []
    assert abs(report["total_charge"] + 0.0001) < 1e-3
    assert report["phosphonate_anchors"] == 1
    ff = (tmp_path / "out" / "force_field.lmp").read_text()
    assert "buck/coul/long 17397.5 0.35" in ff and "# Ni-" in ff


def test_reprotonation_and_hydroxyl_refusal(tmp_path):
    ssym, sxyz, cell = _slab()
    msym, mxyz = _molecule("me-4pacz", 14.0)
    # move one acid H onto a top-layer lattice O
    acid_h = [k for k, s in enumerate(msym) if s == "H"
              and any(msym[o] == "O" and math.dist(mxyz[k], mxyz[o]) < 1.2 for o in range(len(msym)))]
    k = acid_h[0]
    top_o = max((i for i, s in enumerate(ssym) if s == "O"), key=lambda i: sxyz[i][2])
    mxyz[k] = (sxyz[top_o][0], sxyz[top_o][1], sxyz[top_o][2] + 0.97)
    report = build_from_structure("me-4pacz", ssym + msym, sxyz + mxyz, cell, tmp_path / "re", reprotonate=True)
    assert len(report["moved_h"]) == 1 and report["atoms_total"] == 245
    with pytest.raises(ValueError):
        build_from_structure("me-4pacz", ssym + msym, sxyz + mxyz, cell, tmp_path / "no", reprotonate=False)
    # an extra surface hydroxyl H cannot be represented
    extra = [(sxyz[top_o][0] + 2.0, sxyz[top_o][1], sxyz[top_o][2] + 0.97)]
    with pytest.raises(ValueError):
        build_from_structure("me-4pacz", ssym + msym + ["H"], sxyz + mxyz + extra, cell, tmp_path / "oh")
