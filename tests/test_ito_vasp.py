"""Checks for the ITO/SAM VASP reference-input generator (nio_md_prep.ito.vasp)."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pytest

from nio_md_prep.ito import substrate as substrate_mod
from nio_md_prep.ito import vasp
from nio_md_prep.ito.substrate import trilayer_spacing

ROOT = Path(__file__).parents[1]
COMMITTED = ROOT / "inputs" / "ito" / "vasp_smoke"
FORMAL = {"In": 3, "Sn": 4, "O": -2, "H": 1}


@pytest.fixture(scope="module")
def built():
    cases = vasp.build_cases(vasp.Settings())
    files = {}
    for c in cases:
        for name, text in vasp.case_files(c, vasp.Settings()).items():
            files[f"{c.name}/{name}"] = text
    return {c.name: c for c in cases}, files


def _manifest(files, case):
    return json.loads(files[f"{case}/case_manifest.json"])


def _read(text: str, fmt: str):
    from ase.io import read
    return read(io.StringIO(text), format=fmt)


def test_case_list(built):
    cases, _ = built
    assert list(cases) == ["slab-bare", "slab-oh", "slab-oh-sn",
                           "mol-me-4pacz", "ads-me-4pacz-phys", "ads-me-4pacz-bidentate",
                           "mol-meo-2pacz", "ads-meo-2pacz-phys", "ads-meo-2pacz-bidentate"]


def test_render_is_deterministic(built):
    _, files = built
    again = vasp.render(vasp.Settings())
    for rel, text in files.items():
        assert again[rel] == text, rel


def test_substrate_stoichiometry_and_ionic_neutrality(built):
    cases, _ = built
    assert cases["slab-bare"].counts() == {"In": 64, "O": 96}
    assert cases["slab-oh"].counts() == {"In": 64, "O": 102, "H": 12}          # 6 dissociated H2O
    assert cases["slab-oh-sn"].counts() == {"In": 62, "Sn": 2, "O": 103, "H": 12}  # + 2 Sn_In + O_i
    for name in ("slab-bare", "slab-oh", "slab-oh-sn"):
        assert sum(FORMAL[e] * n for e, n in cases[name].counts().items()) == 0, name


def test_molecule_and_adsorbate_stoichiometry(built):
    cases, _ = built
    assert vasp._formula(cases["mol-me-4pacz"].counts()) == "C18H22NO3P"
    assert vasp._formula(cases["mol-meo-2pacz"].counts()) == "C16H18NO5P"
    oh = cases["slab-oh"].counts()
    for slug in vasp.MOLECULES:
        mol = cases[f"mol-{slug}"].counts()
        want = {e: oh.get(e, 0) + mol.get(e, 0) for e in vasp.SPECIES_ORDER if oh.get(e, 0) + mol.get(e, 0)}
        assert cases[f"ads-{slug}-phys"].counts() == want
        # Bidentate conserves stoichiometry: E(bid) - E(phys) needs no reference correction.
        assert cases[f"ads-{slug}-bidentate"].counts() == cases[f"ads-{slug}-phys"].counts()


def test_electron_counts(built):
    cases, files = built
    expected = {"slab-bare": 1408, "slab-oh": 1456, "slab-oh-sn": 1464, "mol-me-4pacz": 122, "mol-meo-2pacz": 122,
                "ads-me-4pacz-phys": 1578, "ads-me-4pacz-bidentate": 1578,
                "ads-meo-2pacz-phys": 1578, "ads-meo-2pacz-bidentate": 1578}
    for name, ne in expected.items():
        assert vasp.nelect(cases[name].counts()) == ne
        assert _manifest(files, name)["nelect_expected"] == ne
        spec = files[f"{name}/POTCAR.spec"]
        assert f"# expected NELECT = {ne}" in spec
        rows = [line.split() for line in spec.splitlines() if line and not line.startswith("#")]
        assert sum(int(r[2]) * int(r[3]) for r in rows) == ne
        poscar = files[f"{name}/POSCAR"].splitlines()
        assert poscar[5].split() == [r[0] for r in rows]          # POTCAR order == POSCAR species order
        assert [int(x) for x in poscar[6].split()] == [int(r[2]) for r in rows]
        assert all(r[1] == vasp.POTCAR_DATASETS[r[0]] for r in rows)


def test_selective_dynamics(built):
    cases, files = built
    d = trilayer_spacing()
    for name, c in cases.items():
        lines = files[f"{name}/POSCAR"].splitlines()
        if c.kind == "molecule":
            assert "Selective dynamics" not in lines
            continue
        assert lines[7] == "Selective dynamics"
        flags = [line.split()[3:] for line in lines[9:]]
        n_fixed = sum(1 for f in flags if f == ["F", "F", "F"])
        assert n_fixed + sum(1 for f in flags if f == ["T", "T", "T"]) == len(c.symbols)
        assert n_fixed == _manifest(files, name)["constraints"]["fixed_atoms"]
        assert n_fixed == (70 if name == "slab-oh-sn" else 80)     # one trilayer = In32O48 (minus freed Sn shell)
        z = c.positions[:, 2]; slab = np.array([l != "mol" for l in c.labels])
        assert np.all(z[c.fixed] < z[slab].min() + d)             # only the bottom trilayer is fixed
        assert not np.any(c.fixed[~slab])


def test_poscar_and_extxyz_agree_and_parse(built):
    cases, files = built
    for name, c in cases.items():
        a = _read(files[f"{name}/POSCAR"], "vasp")
        b = _read(files[f"{name}/geometry.extxyz"], "extxyz")
        assert a.get_chemical_symbols() == b.get_chemical_symbols()
        assert np.allclose(a.positions, b.positions, atol=1e-7)
        assert np.allclose(a.cell.lengths(), [14.3076, 24.7815, c.notes["cell_c"]], atol=1e-3)
        fixed = set(a.constraints[0].index) if a.constraints else set()
        assert fixed == set(np.where(b.arrays["sd_fixed"])[0])


def test_minimum_distances(built):
    cases, files = built
    for name, c in cases.items():
        d = _manifest(files, name)["distances_angstrom"]
        assert d["min_nonbonded"] > 1.5, name
        a = _read(files[f"{name}/POSCAR"], "vasp")
        dist = a.get_all_distances(mic=True)[np.triu_indices(len(a), 1)]
        assert dist.min() > 0.95, name                           # shortest bond is O-H 0.965
        lo, hi = d["bonded_range"]
        assert 0.95 < lo and hi < 2.6, name
        if "min_molecule_slab_nonbonded" in d:
            assert d["min_molecule_slab_nonbonded"] > 1.8, name


def test_hydroxylation_is_clash_free(built):
    cases, _ = built
    h = cases["slab-oh"].notes["hydroxylation"]
    assert h["pairs"] == 6 and len(h["sites"]) == 6
    assert all(s["terminal_o_to_protonated_o"] >= substrate_mod.MIN_OO_TERMINAL for s in h["sites"])


def test_sn_cluster_in_lower_trilayer(built):
    cases, _ = built
    cl = cases["slab-oh-sn"].notes["sn_cluster"]
    assert cl["trilayer"] == 0
    assert 0.0 <= cl["o_i_site_angstrom"][2] < trilayer_spacing()
    assert max(cl["sn_to_oi_angstrom"]) < 2.6


def test_bidentate_geometry(built):
    cases, _ = built
    for slug in vasp.MOLECULES:
        b = cases[f"ads-{slug}-bidentate"].notes["bidentate"]
        assert len(b["in_o_bonds"]) == 2 and len(b["proton_transfers"]) == 2
        assert all(2.0 <= x["distance"] <= 2.3 for x in b["in_o_bonds"])
        assert len({x["acceptor_o3c_index"] for x in b["proton_transfers"]}) == 2


def test_molecule_reference_is_rigid_copy_of_physisorbed_molecule(built):
    cases, _ = built
    for slug in vasp.MOLECULES:
        m, p = cases[f"mol-{slug}"], cases[f"ads-{slug}-phys"]
        pm = {lid: x for lid, x, l in zip(p.lig_id, p.positions, p.labels) if l == "mol"}
        shift = np.array([pm[lid] - x for lid, x in zip(m.lig_id, m.positions)])
        shift[:, :2] -= np.array(m.lengths) * np.round(shift[:, :2] / np.array(m.lengths))
        assert np.allclose(shift, shift[0], atol=1e-6)
        assert abs(shift[0, 0]) < 1e-6 and abs(shift[0, 1]) < 1e-6


def test_incar_modes(built):
    cases, files = built
    incar = files["ads-me-4pacz-phys/INCAR"]
    for tag in ("IVDW = 12", "ISPIN = 1", "LDIPOL = .TRUE.", "IDIPOL = 3", "ENCUT   = 520", "PREC    = Accurate",
                "ISMEAR = 0", "SIGMA  = 0.05", "EDIFF   = 1E-5", "LREAL   = Auto", "NSW    = 3", "IBRION = 2"):
        assert tag in incar
    assert "LDAU" not in incar.replace("No DFT+U", "")
    assert files["slab-oh/KPOINTS"].splitlines()[2:4] == ["G", " 1 1 1"]
    prod = vasp.incar_text(cases["slab-oh"], vasp.Settings(mode="production"))
    assert "NSW    = 300" in prod and "EDIFFG = -0.03" in prod


def test_sbatch_and_hashes(built):
    _, files = built
    for rel, text in files.items():
        if not rel.endswith("inputs.sha256"):
            continue
        case = rel.split("/")[0]
        for line in text.splitlines():
            digest, name = line.split("  ")
            assert hashlib.sha256(files[f"{case}/{name}"].encode()).hexdigest() == digest
        sb = files[f"{case}/run.sbatch"]
        assert "#SBATCH -A loni_perovsk27" in sb and "vasp_gam" in sb and "Adapted from lgutsev/InterfaceForge@" in sb
        assert "POTCAR.spec" in sb and "sha256sum -c inputs.sha256" in sb
    assert not any(rel.endswith("/POTCAR") for rel in files)      # licensed data is never generated


def test_production_three_trilayers_is_one_flag_away():
    slab = vasp.build_substrate(vasp.Settings(trilayers=3), True, False)
    assert slab.counts() == {"In": 96, "O": 150, "H": 12}
    assert int(slab.fixed.sum()) == 80


@pytest.mark.skipif(not COMMITTED.exists(), reason="generated inputs not present")
def test_committed_inputs_match_generator(built):
    _, files = built
    for rel, text in files.items():
        if rel.endswith("case_manifest.json"):
            continue                                              # carries generator source hashes
        path = COMMITTED / rel
        assert path.exists(), rel
        assert path.read_bytes().replace(b"\r\n", b"\n") == text.encode("utf-8"), rel
