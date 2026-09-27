"""Single-molecule placement and interaction-energy checks on rigid ITO slabs.

For every (substrate, parameter set, molecule, placement) one neutral
molecule is placed above the slab, minimized with the slab frozen, and
decomposed as

    E_int  = E(complex) - E(slab) - E(molecule at the complex geometry)
    E_ads  = E(complex) - E(slab) - E(relaxed isolated molecule)

All three evaluations use one LAMMPS box and identical, explicitly fixed PPPM
mesh/g_ewald settings.  PPPM (and the slab dipole correction) is a quadratic
form of the assigned charge density, so the slab-slab terms cancel exactly
rather than to the PPPM accuracy; real-space terms are pairwise.  The results
are classical-FF numbers for model comparison only, not binding energies.
"""
from __future__ import annotations

import json
import math
import tomllib
from pathlib import Path

import numpy as np

from ..build import _coeff_lines, _surface
from ..chemistry import correction_lines, phosphonate_roles
from ..config import ROOT, molecule_manifest
from ..geometry import elements
from ..lammps import DataFile, parse, replicate, write
from .assemble import _strip_pair_style
from .forcefield import PARAMETER_SETS, surface_pair_lines

KSPACE = "kspace_style pppm 1e-5\nkspace_modify slab 3.0 diff ad gewald {g} mesh {mx} {my} {mz}\n"
STYLES = """units real
atom_style full
boundary p p f
bond_style harmonic
angle_style harmonic
dihedral_style hybrid opls charmm
improper_style cvff
pair_style lj/cut/coul/long 10.0 8.0
pair_modify mix geometric
special_bonds amber
"""


def _rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * k @ k


def _align(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    u = u / np.linalg.norm(u); v = v / np.linalg.norm(v)
    c = float(np.dot(u, v))
    if c > 1 - 1e-12: return np.eye(3)
    if c < -1 + 1e-12:
        perp = np.cross(u, [1, 0, 0]) if abs(u[0]) < 0.9 else np.cross(u, [0, 1, 0])
        return _rotation(perp, math.pi)
    return _rotation(np.cross(u, v), math.acos(c))


def molecule_frame(data: DataFile) -> dict:
    """Anchor (P) index, heavy-atom long axis, and element list of a LigParGen molecule."""
    xyz = np.array([[float(a.fields[4]), float(a.fields[5]), float(a.fields[6])] for a in data.sections["Atoms"]])
    el = elements(data)
    roles = phosphonate_roles(data)
    ids = [int(a.fields[0]) for a in data.sections["Atoms"]]
    p_idx = [ids.index(i) for i, r in roles.items() if r == "P"]
    heavy = [k for k, e in enumerate(el) if e not in ("H",)]
    p = xyz[p_idx].mean(axis=0)
    far = max(heavy, key=lambda k: np.linalg.norm(xyz[k] - p))
    # Core = heavy C/N atoms in the far 30 % (as in interfacial._core_atoms, by distance here).
    dist = np.linalg.norm(xyz - p, axis=1)
    cn = [k for k, e in enumerate(el) if e in ("C", "N")]
    core = sorted(cn, key=lambda k: -dist[k])[: max(3, int(round(0.3 * len(cn))))]
    return {"xyz": xyz, "elements": el, "p": p_idx, "axis": xyz[far] - p, "core": core,
            "anchor_o": [ids.index(i) for i, r in roles.items() if r in ("P=O", "P-OH")],
            "acid_h": [ids.index(i) for i, r in roles.items() if r == "P-O-H"]}


def place(frame: dict, xy: tuple[float, float], z_floor: float, tilt_deg: float, azimuth: float, spin: float) -> np.ndarray:
    """Anchor-down placement: P->far axis at ``tilt`` from +z, lowest atom at ``z_floor``."""
    xyz = frame["xyz"] - frame["xyz"][frame["p"]].mean(axis=0)
    xyz = xyz @ _rotation(frame["axis"], spin).T
    t = math.radians(tilt_deg)
    target = np.array([math.sin(t) * math.cos(azimuth), math.sin(t) * math.sin(azimuth), math.cos(t)])
    xyz = xyz @ _align(frame["axis"], target).T
    xyz[:, 0] += xy[0]; xyz[:, 1] += xy[1]
    xyz[:, 2] += z_floor - xyz[:, 2].min()
    return xyz


def _system(slab: DataFile, mol: DataFile, coords: np.ndarray, zhi: float) -> tuple[DataFile, dict, list[str]]:
    ids = {"atom": 0, "bond": 0, "angle": 0, "dihedral": 0, "improper": 0}; types = {k: 0 for k in ids}
    result = DataFile("ITO single-molecule check", {}, {})
    pieces, inc = replicate(mol, 1, types, ids, 0, [tuple(map(float, c)) for c in coords])
    for sec, rows in pieces.items(): result.sections.setdefault(sec, []).extend(rows)
    chem, _, _ = correction_lines(mol, 0, 0)
    for k in ids: ids[k] += inc.get(k, 0); types[k] += mol.type_count(k)
    srows, _ = _surface(slab, types, ids, 0)
    for sec, rows in srows.items(): result.sections.setdefault(sec, []).extend(rows)
    result.bounds = {"x": slab.bounds["x"], "y": slab.bounds["y"], "z": (slab.bounds["z"][0], zhi)}
    old_to_new = {int(r.fields[0]): types["atom"] + i + 1 for i, r in enumerate(slab.sections["Masses"])}
    lines = [_strip_pair_style(x) for x in _coeff_lines(result)] + [_strip_pair_style(x) for x in chem]
    for sec in ("Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs"):
        result.sections.pop(sec, None)
    return result, old_to_new, lines


def _contacts(xyz_mol: np.ndarray, frame: dict, slab_xyz: np.ndarray, slab_labels: list[str], lengths, top: float) -> dict:
    def pdist(a, b):
        d = b[None, :, :] - a[:, None, :]
        for ax in (0, 1): d[..., ax] -= lengths[ax] * np.round(d[..., ax] / lengths[ax])
        return np.linalg.norm(d, axis=2)
    lab = np.array(slab_labels)
    cat = slab_xyz[np.isin(lab, ["In", "Sn"])]
    ox = slab_xyz[np.isin(lab, ["O", "Oh"])]
    hh = slab_xyz[lab == "Hh"]
    ao = xyz_mol[frame["anchor_o"]]; ah = xyz_mol[frame["acid_h"]]
    out = {"anchor_o_min_cation_angstrom": float(pdist(ao, cat).min()),
           "acid_h_min_surface_o_angstrom": float(pdist(ah, ox).min()) if len(ah) else None,
           "anchor_o_min_surface_h_angstrom": float(pdist(ao, hh).min()) if len(hh) else None,
           "p_height_above_top_atom_angstrom": float(xyz_mol[frame["p"], 2].mean() - top)}
    hb = int((pdist(ah, ox) < 2.5).any(axis=1).sum()) if len(ah) else 0
    hb += int((pdist(ao, hh) < 2.5).any(axis=1).sum()) if len(hh) else 0
    out["anchor_hbonds_lt_2p5"] = hb
    out["anchor_o_cation_contacts_lt_3p25"] = int((pdist(ao, cat) < 3.25).any(axis=1).sum())
    v = xyz_mol[frame["core"]].mean(axis=0) - xyz_mol[frame["p"]].mean(axis=0)
    out["tilt_deg"] = float(math.degrees(math.acos(max(-1.0, min(1.0, v[2] / np.linalg.norm(v))))))
    out["min_molecule_z_above_top_angstrom"] = float(xyz_mol[:, 2].min() - top)
    return out


def _lammps():
    from lammps import lammps
    return lammps(cmdargs=["-screen", "none", "-log", "none", "-nocite"])


def _energy_run(data_path: Path, ff_lines: list[str], kspace: str, n_mol: int, n_slab: int, minimize: bool,
                delete: str | None = None, min_steps: int = 3000, quench_steps: int = 0, seed: int = 1,
                temperature: float = 300.0) -> tuple[float, np.ndarray]:
    """Energy of one evaluation; optional minimize -> NVT quench of the molecule -> minimize."""
    lmp = _lammps()
    try:
        lmp.commands_string(STYLES + f"read_data {data_path.as_posix()}\n" + "\n".join(ff_lines) + "\n" + kspace)
        lmp.commands_string(f"group mol id 1:{n_mol}\ngroup slab id {n_mol+1}:{n_mol+n_slab}\n")
        # Slab-slab pairs are excluded in every evaluation so they cancel in the differences.
        lmp.command("neigh_modify exclude group slab slab")
        if delete:
            lmp.command(f"delete_atoms group {delete} compress no bond yes")
        lmp.commands_string("neighbor 2.0 bin\nneigh_modify every 1 delay 0 check yes\n")
        if delete != "slab":
            lmp.command("fix freeze slab setforce 0.0 0.0 0.0")
        minimize_cmd = f"min_style cg\nmin_modify dmax 0.1 line quadratic\nminimize 0.0 1.0e-3 {min_steps} {10*min_steps}\n"
        if minimize:
            lmp.commands_string(minimize_cmd)
            if quench_steps > 0:
                lmp.commands_string(
                    f"velocity mol create {temperature} {seed} mom yes rot yes dist gaussian\n"
                    f"fix q mol nvt temp {temperature} 10.0 50.0\n"
                    "fix wlo mol wall/lj93 zlo EDGE 1.0 1.0 2.5 units box\n"
                    "fix whi mol wall/lj93 zhi EDGE 1.0 1.0 2.5 units box\n"
                    f"timestep 0.5\nrun {quench_steps}\nunfix q\nunfix wlo\nunfix whi\n")
                lmp.commands_string(minimize_cmd)
        lmp.command("run 0 post no")
        e = float(lmp.get_thermo("pe"))
        n = lmp.get_natoms()
        x = np.array(lmp.numpy.extract_atom("x"))[:n].copy()
        ids = np.array(lmp.numpy.extract_atom("id"))[:n].copy()
        order = np.argsort(ids)
        mol_xyz = x[order][ids[order] <= n_mol] if delete != "mol" else np.empty((0, 3))
        return e, mol_xyz
    finally:
        lmp.close()


def _placements(sc: dict) -> list[dict]:
    rng = np.random.default_rng(int(sc.get("seed", 1)))
    sites = [(float(rng.uniform()), float(rng.uniform())) for _ in range(int(sc["lateral_positions"]))]
    out = []
    for k, (fx, fy) in enumerate(sites):  # one lateral site shared by all tilts of that index
        for tilt in sc["tilts_deg"]:
            out.append({"index": len(out), "lateral": k, "tilt_deg": float(tilt), "fx": fx, "fy": fy,
                        "azimuth": float(rng.uniform(0, 2 * math.pi)), "spin": float(rng.uniform(0, 2 * math.pi))})
    return out


def _with_coords(path: Path, out: Path, xyz: np.ndarray) -> Path:
    data = parse(path)
    for a, x in zip(data.sections["Atoms"][: len(xyz)], xyz):
        a.fields[4:7] = [f"{v:.8f}" for v in x]
    write(data, out)
    return out


def _scan_tag(job: tuple) -> list[dict]:
    sc, sub, slug, pset, output = job
    sdir = ROOT / "inputs" / "ito" / "surfaces" / sub
    sman = json.loads((sdir / "surface_manifest.json").read_text(encoding="utf-8"))
    slab = parse(sdir / "surface.lmp")
    lengths = (slab.bounds["x"][1], slab.bounds["y"][1])
    top = float(sman["z_top_atom_angstrom"]); zhi = top + float(sc.get("box_height_above_top", 60.0))
    inv = {v: k for k, v in sman["type_ids"].items()}
    slab_labels = [inv[int(a.fields[2])] for a in slab.sections["Atoms"]]
    slab_xyz = np.array([[float(a.fields[4]), float(a.fields[5]), float(a.fields[6])] for a in slab.sections["Atoms"]])
    mesh = sc.get("mesh", [48, 42, 160])
    kspace = KSPACE.format(g=float(sc.get("gewald", 0.30)), mx=mesh[0], my=mesh[1], mz=mesh[2])
    folder, man = molecule_manifest(slug)
    mol = parse(folder / man["files"]["ligpargen"])
    frame = molecule_frame(mol)
    n_mol = mol.count("Atoms"); n_slab = slab.count("Atoms")
    quench = int(sc.get("quench_steps", 0))
    min_steps = int(sc.get("min_steps", 3000)); ref_steps = int(sc.get("reference_min_steps", 10000))
    tag = f"{sub}__{slug}__{pset.replace('/', '-')}"
    work = Path(output) / tag; work.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []; e_slab = None; refs: list[float] = []
    for pl in _placements(sc):
        coords = place(frame, (pl["fx"] * lengths[0], pl["fy"] * lengths[1]), top + float(sc.get("floor_gap", 2.5)),
                       pl["tilt_deg"], pl["azimuth"], pl["spin"])
        system, old_to_new, lines = _system(slab, mol, coords, zhi)
        ff = lines + surface_pair_lines({lab: old_to_new[old] for lab, old in sman["type_ids"].items()}, pset)
        path = work / f"placement-{pl['index']:02d}.data"; write(system, path)
        if e_slab is None:
            e_slab, _ = _energy_run(path, ff, kspace, n_mol, n_slab, False, delete="mol")
            refs.append(_energy_run(path, ff, kspace, n_mol, n_slab, True, delete="slab", min_steps=ref_steps,
                                    quench_steps=quench, seed=7)[0])
        e0, _ = _energy_run(path, ff, kspace, n_mol, n_slab, False)
        e_cx, xyz_min = _energy_run(path, ff, kspace, n_mol, n_slab, True, min_steps=min_steps, quench_steps=quench, seed=100 + pl["index"])
        fpath = _with_coords(path, work / f"placement-{pl['index']:02d}-min.data", xyz_min)
        e_mfix, _ = _energy_run(fpath, ff, kspace, n_mol, n_slab, False, delete="slab")
        refs.append(_energy_run(fpath, ff, kspace, n_mol, n_slab, True, delete="slab", min_steps=ref_steps)[0])
        (work / f"placement-{pl['index']:02d}-min.xyz").write_text(
            f"{n_mol}\n{tag} placement {pl['index']} minimized\n"
            + "".join(f"{e} {x:.5f} {y:.5f} {z:.5f}\n" for e, (x, y, z) in zip(frame["elements"], xyz_min)), encoding="utf-8")
        fpath.unlink(missing_ok=True)
        rows.append({"substrate": sub, "molecule": slug, "parameter_set": pset,
                     "index": pl["index"], "lateral": pl["lateral"], "initial_tilt_deg": pl["tilt_deg"],
                     "initial_xy_fraction": [round(pl["fx"], 4), round(pl["fy"], 4)], "quench_steps": quench,
                     "e_complex_initial": e0, "e_complex_min": e_cx, "e_slab": e_slab, "e_mol_at_complex_geometry": e_mfix,
                     "e_int_kcal_mol": e_cx - e_slab - e_mfix,
                     **_contacts(xyz_min, frame, slab_xyz, slab_labels, lengths, top)})
        r = rows[-1]
        print(f"{tag} p{pl['index']:02d} tilt0={pl['tilt_deg']:>3.0f} E_int={r['e_int_kcal_mol']:8.2f} "
              f"tilt={r['tilt_deg']:5.1f} hb={r['anchor_hbonds_lt_2p5']} cat={r['anchor_o_cation_contacts_lt_3p25']}", flush=True)
    # The lowest isolated-molecule energy reached from any start is the gas-phase reference.
    e_ref = min(refs)
    for r in rows:
        r["e_mol_reference"] = e_ref
        r["e_ads_kcal_mol"] = r["e_complex_min"] - r["e_slab"] - e_ref
        r["e_deformation_kcal_mol"] = r["e_mol_at_complex_geometry"] - e_ref
    (work / "rows.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return rows


def summarize_scan(rows: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        groups.setdefault((r["substrate"], r["molecule"], r["parameter_set"]), []).append(r)
    out = []
    for (sub, slug, pset), rs in sorted(groups.items()):
        e = np.array([r["e_ads_kcal_mol"] for r in rs]); ei = np.array([r["e_int_kcal_mol"] for r in rs])
        out.append({"substrate": sub, "molecule": slug, "parameter_set": pset, "n": len(rs),
                    "e_ads_min": float(e.min()), "e_ads_median": float(np.median(e)), "e_ads_max": float(e.max()),
                    "e_int_min": float(ei.min()), "e_int_median": float(np.median(ei)),
                    "fraction_with_anchor_hbond": float(np.mean([r["anchor_hbonds_lt_2p5"] > 0 for r in rs])),
                    "fraction_with_anchor_cation_contact": float(np.mean([r["anchor_o_cation_contacts_lt_3p25"] > 0 for r in rs])),
                    "tilt_median_deg": float(np.median([r["tilt_deg"] for r in rs])),
                    "p_height_median_angstrom": float(np.median([r["p_height_above_top_atom_angstrom"] for r in rs]))})
    return out


def adsorption_scan(config_path: Path, output: Path, workers: int | None = None) -> Path:
    import os
    from concurrent.futures import ProcessPoolExecutor
    with config_path.open("rb") as f:
        cfg = tomllib.load(f)
    output.mkdir(parents=True, exist_ok=True)
    sc = cfg["scan"]
    jobs = [(sc, sub, slug, pset, str(output)) for sub in cfg["substrates"] for slug in cfg["molecules"] for pset in cfg["parameter_sets"]]
    workers = workers or min(len(jobs), max(1, (os.cpu_count() or 2) // 2))
    if workers == 1:
        results = [_scan_tag(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(_scan_tag, jobs))
    rows = [r for chunk in results for r in chunk]
    (output / "scan_results.json").write_text(
        json.dumps({"config": cfg, "rows": rows, "summary": summarize_scan(rows)}, indent=2) + "\n", encoding="utf-8")
    return output / "scan_results.json"
