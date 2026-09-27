"""ITO extension: slab construction invariants, force-field text, and analysis helpers."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

pytest.importorskip("ase")

from nio_md_prep.ito import substrate as sub
from nio_md_prep.ito.assemble import _stage_inputs, _strip_pair_style
from nio_md_prep.ito.forcefield import PARAMETER_SETS, surface_pair_lines


def _slab(nx=1, ny=1, trilayers=4):
    return sub.build_slab(nx, ny, trilayers)


def test_bixbyite_cell_is_in32o48_with_sixfold_indium():
    symbols, pos, lengths = sub.oriented_cell()
    assert symbols.count("In") == 192 and symbols.count("O") == 288 and symbols.count("X") == 96
    assert np.allclose(lengths, [10.117 * math.sqrt(2), 10.117 * math.sqrt(6), 10.117 * math.sqrt(3)], atol=1e-6)


@pytest.mark.parametrize("trilayers", [2, 3, 6])
def test_slab_is_stoichiometric_neutral_and_dipole_free(trilayers):
    slab, _ = _slab(1, 1, trilayers)
    assert slab.count("In") == 32 * trilayers and slab.count("O") == 48 * trilayers
    q = sub.charges(slab, 0.525, 0.425)
    assert abs(q.sum()) < 1e-9
    assert abs(sub._dipole_per_area(slab, q)) < 1e-9
    if trilayers < 3:
        return  # no bulk-like interior to check
    cn, _ = sub.coordination(slab)
    z = slab.positions[:, 2]
    interior = (z > z.min() + 3.0) & (z < z.max() - 3.0)
    lab = np.array(slab.labels)
    assert set(cn[interior & (lab == "In")]) == {6}
    assert set(cn[interior & (lab == "O")]) == {4}


def test_top_face_has_12_in5c_and_12_o3c_per_primitive_cell():
    slab, _ = _slab(1, 1, 3)
    cn, _ = sub.coordination(slab)
    z = slab.positions[:, 2]; top = z > z.mean(); lab = np.array(slab.labels)
    # Orthorhombic cell = two primitive (111) surface cells.
    assert int(((lab == "In") & (cn == 5) & top).sum()) == 24
    assert int(((lab == "O") & (cn == 3) & top).sum()) == 24


def test_hydroxylation_and_doping_stay_neutral_and_deterministic():
    runs = []
    for _ in range(2):
        slab, sites = _slab(2, 1, 4)
        doping = sub.dope_sn(slab, sites, 0.09, seed=5)
        hyd = sub.hydroxylate(slab, 6, seed=3)
        q = sub.charges(slab, 0.525, 0.425)
        assert abs(q.sum()) < 1e-9
        assert slab.count("Hh") == 12 and slab.count("Oh") == 12
        assert slab.count("Sn") == 2 * doping["clusters"]
        assert hyd["min_new_atom_nonbonded_angstrom"] >= 1.5
        runs.append((list(slab.labels), slab.positions.copy()))
    assert runs[0][0] == runs[1][0] and np.array_equal(runs[0][1], runs[1][1])


def test_clayff_charge_pattern():
    slab, _ = _slab(1, 1, 2)
    sub.hydroxylate(slab, 2, seed=1)
    q = sub.charges(slab, 0.525, 0.425)
    by = {lab: round(float(q[slab.labels.index(lab)]), 6) for lab in set(slab.labels)}
    assert by == {"In": 1.575, "O": -1.05, "Oh": -0.95, "Hh": 0.425}


def test_committed_models_match_their_manifests(tmp_path):
    from nio_md_prep.config import ROOT
    for model in ("in2o3-111-bare", "in2o3-111-oh", "ito-111-oh", "in2o3-111-oh-groove"):
        folder = ROOT / "inputs" / "ito" / "surfaces" / model
        man = sub.build_from_model(folder / "model.toml", tmp_path / model)
        committed = json.loads((folder / "surface_manifest.json").read_text(encoding="utf-8"))
        assert man["surface_lmp_sha256"] == committed["surface_lmp_sha256"]
        assert (tmp_path / model / "surface.lmp").read_bytes() == (folder / "surface.lmp").read_bytes()


def test_surface_pair_lines_cover_every_label_and_set():
    for name, params in PARAMETER_SETS.items():
        cation = [c for c in ("In", "Sn", "Ni") if c in params]
        ids = {lab: 10 + k for k, lab in enumerate(cation + ["O", "Oh", "Hh"])}
        lines = surface_pair_lines(ids, name)
        assert [int(l.split()[1]) for l in lines] == sorted(ids.values())
        assert all("lj/cut/coul/long" not in l for l in lines)


def test_strip_pair_style_only_touches_pair_coeff():
    assert _strip_pair_style("pair_coeff 3 3 lj/cut/coul/long 0.2 3.7 # x") == "pair_coeff 3 3 0.2 3.7 # x"
    assert _strip_pair_style("dihedral_coeff 3 charmm 0.5 0 3 0.0") == "dihedral_coeff 3 charmm 0.5 0 3 0.0"


def test_stage_inputs_rigid_slab_nvt_and_shake():
    p = {"temperature": 300.0, "timestep": 1.0, "deposition_steps": 10, "hold_steps": 10, "release_steps": 10,
         "relax_steps": 10, "wall_clearance": 30.0, "release_height": 40.0, "dump_every": 5, "tdamp": 100.0, "shake": True}
    st = _stage_inputs(p, 19.0, 125.0, 7, 80.0, [2, 5])
    for name in ("deposition", "hold", "relax"):
        text = st[name]
        assert "neigh_modify exclude group slab slab" in text and "fix freeze slab setforce" in text
        assert "npt" not in text and "fix ensemble mobile nvt" in text
        assert "fix constrain mobile shake 0.0001 20 0 b 2 5" in text
    dep = st["deposition"]
    assert dep.index("minimize") < dep.index("fix constrain")
    with pytest.raises(ValueError):
        _stage_inputs(p | {"release_height": 200.0}, 19.0, 125.0, 7, 80.0, [2])


def test_periodic_clusters_merge_across_the_boundary():
    from nio_md_prep.ito.analysis import _periodic_clusters
    xy = np.array([[0.5, 5.0], [9.6, 5.0], [5.0, 5.0]])
    assert _periodic_clusters(xy, (10.0, 10.0), 1.5) == [2, 1]


def test_energy_decomposition_vanishes_at_large_separation(tmp_path):
    """Fixed-mesh PPPM decomposition: E_int -> 0 far from the slab, so slab-slab terms cancel exactly."""
    pytest.importorskip("lammps")
    from nio_md_prep.config import molecule_manifest
    from nio_md_prep.ito import checks
    from nio_md_prep.lammps import parse, write
    model = tmp_path / "model.toml"
    model.write_text('[model]\nid = "t"\n[slab]\nnx = 2\nny = 1\ntrilayers = 2\nzlo = -5.0\nzhi = 80.0\n'
                     '[charges]\nscale = 0.525\nhydroxyl_h = 0.425\n[hydroxylation]\npairs_per_primitive_cell = 3\nseed = 1\n')
    man = sub.build_from_model(model, tmp_path)
    slab = parse(tmp_path / "surface.lmp")
    folder, mm = molecule_manifest("me-4pacz")
    mol = parse(folder / mm["files"]["ligpargen"])
    frame = checks.molecule_frame(mol)
    ks = checks.KSPACE.format(g=0.30, mx=32, my=24, mz=96)
    xyz = checks.place(frame, (14.0, 12.0), man["z_top_atom_angstrom"] + 35.0, 0.0, 0.3, 0.2)
    system, o2n, lines = checks._system(slab, mol, xyz, 80.0)
    ff = lines + surface_pair_lines({lab: o2n[o] for lab, o in man["type_ids"].items()}, "uff-cation/clayff-anion")
    path = tmp_path / "far.data"; write(system, path)
    n, ns = mol.count("Atoms"), slab.count("Atoms")
    e_cx = checks._energy_run(path, ff, ks, n, ns, False)[0]
    e_slab = checks._energy_run(path, ff, ks, n, ns, False, delete="mol")[0]
    e_mol = checks._energy_run(path, ff, ks, n, ns, False, delete="slab")[0]
    assert abs(e_cx - e_slab - e_mol) < 0.05


def test_groove_slab_is_neutral_with_nio_like_profile(tmp_path):
    model = tmp_path / "model.toml"
    model.write_text("\n".join([
        '[model]', 'id = "g"',
        '[slab]', 'nx = 3', 'ny = 1', 'trilayers = 8', 'zlo = -5.0', 'zhi = 80.0',
        '[corrugation]', 'depth_trilayers = 5', 'wall_slope = 1.0', 'profile_axis = "x"',
        '[charges]', 'scale = 0.525', 'hydroxyl_h = 0.425',
        '[hydroxylation]', 'pairs_per_primitive_cell = 3', 'seed = 1', '']))
    man = sub.build_from_model(model, tmp_path)
    assert abs(man["total_charge"]) < 1e-9
    c = man["corrugation"]
    assert c["profile_axis"] == "y" and abs(c["depth_angstrom"] - 14.6026) < 1e-3
    assert man["box_angstrom"]["y"][1] == pytest.approx(3 * 10.117 * math.sqrt(2), abs=1e-6)
    assert man["hydroxylation"]["min_new_atom_nonbonded_angstrom"] >= 1.5
    # Lattice-only top heights: plateau minus floor ~ the 5-trilayer depth, minimum at the groove centre.
    from nio_md_prep.lammps import parse
    data = parse(tmp_path / "surface.lmp")
    lattice = {man["type_ids"]["In"], man["type_ids"]["O"]}
    xyz = np.array([[float(a.fields[4]), float(a.fields[5]), float(a.fields[6])] for a in data.sections["Atoms"] if int(a.fields[2]) in lattice])
    ly = man["box_angstrom"]["y"][1]
    tops = np.array([xyz[(xyz[:, 1] >= a) & (xyz[:, 1] < a + 2.0), 2].max() for a in np.arange(0, ly - 2.0, 2.0)])
    centre = c["groove_center_y_angstrom"]
    assert 13.0 < tops.max() - tops.min() < 16.5
    assert abs(np.arange(0, ly - 2.0, 2.0)[np.argmin(tops)] + 1.0 - centre) < 3.0


def test_groove_region_bands():
    from nio_md_prep.ito.analysis import _groove_regions, _region
    corr = {"depth_angstrom": 14.6, "wall_slope": 1.0, "step_run_angstrom": 2.92, "groove_center_y_angstrom": 21.0}
    reg = _groove_regions(corr, (100.0, 42.0))
    assert _region(21.5, reg, 42.0) == "groove_floor"
    assert _region(30.0, reg, 42.0) == "groove_wall"
    assert _region(40.0, reg, 42.0) == "plateau" and _region(1.0, reg, 42.0) == "plateau"
    assert sum(reg["area_nm2"].values()) == pytest.approx(42.0)


def test_local_surface_height_follows_a_step():
    from nio_md_prep.ito.analysis import _local_surface_height
    xs, ys = np.meshgrid(np.arange(0, 20, 1.0), np.arange(0, 20, 1.0))
    z = np.where(ys < 10, 10.0, 5.0)
    h = _local_surface_height(np.column_stack([xs.ravel(), ys.ravel(), z.ravel()]), (20.0, 20.0))
    assert h(5.0, 4.0) == 10.0 and h(5.0, 15.0) == 5.0
