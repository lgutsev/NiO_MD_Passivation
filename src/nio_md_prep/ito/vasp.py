"""VASP reference inputs (smoke and production) for the classical ITO/SAM model.

Generates a small, deterministic set of DFT reference cases on the orthorhombic
1x1 bixbyite In2O3(111) cell (14.31 x 24.78 A) built by :mod:`.substrate`:

* ``slab-bare``, ``slab-oh``, ``slab-oh-sn``: substrate references;
* ``mol-<slug>``: the isolated (gas-phase) protonated phosphonic acid in the
  same cell, rigidly translated from the physisorbed start geometry;
* ``ads-<slug>-phys``: the upright, H-bonded physisorbed acid on ``slab-oh``
  (the configuration class the nonreactive classical model represents);
* ``ads-<slug>-bidentate``: a chemisorbed bridging-bidentate phosphonate on two
  adjacent In5c with both acidic protons moved to surface O3c (same atoms as the
  physisorbed case, so the two energies subtract directly).

Every case directory receives POSCAR, INCAR, KPOINTS, POTCAR.spec (never a
licensed POTCAR), geometry.extxyz, inputs.sha256, run.sbatch and
case_manifest.json.  See ``docs/ito/vasp-smoke.md``.

Adapted from lgutsev/InterfaceForge@4501e34081a63e3735336e693bea1eb40ac35cb1
notebooks/nio_m110_hydroxylation/inputs/vasp_template/INCAR (INCAR style),
profiles/potcar_pbe_54.yaml (PAW dataset choices), src/interfaceforge/vasp.py
(resolve_potcar_root / assemble_potcar search order and refusal rules) and
launch_scripts/runvasp.sh (LONI module and launch line) (MIT).
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..config import ROOT
from ..geometry import elements
from ..lammps import atom_coordinates, parse
from ..chemistry import phosphonate_roles
from .substrate import (BIXBYITE, IN_O_BOND_CUTOFF, IN_OH_LENGTH, O_H_LENGTH, _farthest_point, build_slab, coordination,
                        trilayer_spacing)

INTERFACEFORGE_SHA = "4501e34081a63e3735336e693bea1eb40ac35cb1"
ADAPTED_FROM = (
    f"Adapted from lgutsev/InterfaceForge@{INTERFACEFORGE_SHA} launch_scripts/runvasp.sh and "
    "src/interfaceforge/vasp.py (resolve_potcar_root, assemble_potcar) (MIT)"
)

# PBE PAW choices from InterfaceForge profiles/potcar_pbe_54.yaml.
POTCAR_DATASETS = {"In": "In_d", "Sn": "Sn_d", "O": "O", "P": "P", "N": "N", "C": "C", "H": "H"}
# Valence (ZVAL) of those datasets in potpaw_PBE.54.  Standard values, but not
# read from a POTCAR here (licensed); run.sbatch re-checks them against the
# assembled POTCAR and refuses to launch on any mismatch.
ZVAL = {"In_d": 13, "Sn_d": 14, "O": 6, "P": 5, "N": 5, "C": 4, "H": 1}
# ENMAX (eV) as recalled for potpaw_PBE.54; used only to justify ENCUT.  The
# sbatch prints the real ENMAX values from the assembled POTCAR.
ENMAX_EXPECTED = {"In_d": 239.2, "Sn_d": 241.1, "O": 400.0, "P": 255.0, "N": 400.0, "C": 400.0, "H": 250.0}
SPECIES_ORDER = ("In", "Sn", "O", "P", "N", "C", "H")
CATIONS = ("In", "Sn")
MOLECULES = ("me-4pacz", "meo-2pacz")
HYDROXYL_SEED = 20260926  # same seed as inputs/ito/surfaces/in2o3-111-oh/model.toml
MIN_NONBONDED = 1.5


@dataclass(frozen=True)
class Settings:
    mode: str = "smoke"
    trilayers: int = 2
    hydroxyl_pairs_per_primitive: int = 3
    hydroxyl_seed: int = HYDROXYL_SEED
    sn_trilayer: int = 0
    phys_height: float = 3.75          # P above the top lattice-O plane (A)
    bidentate_in_o: float = 2.15       # target In-O(P) bond length (A)
    vacuum_slab: float = 15.0
    vacuum_molecule: float = 12.0
    common_cell: bool = True
    encut: float = 520.0
    molecules: tuple[str, ...] = MOLECULES

    def __post_init__(self):
        if self.mode not in ("smoke", "production"):
            raise ValueError("mode must be 'smoke' or 'production'")
        if self.trilayers < 2:
            raise ValueError("at least two trilayers are required (one fixed, one relaxed)")
        if not 0 <= self.sn_trilayer < self.trilayers - 1:
            raise ValueError("the Sn cluster must sit below the top (adsorption) trilayer")


@dataclass
class Case:
    name: str
    kind: str                      # slab | molecule | adsorbate
    purpose: str
    symbols: list[str]
    positions: np.ndarray
    labels: list[str]              # substrate labels (In/Sn/O/Oh/Hh) or "mol"
    lig_id: list[int]              # LigParGen atom id for molecule atoms, 0 otherwise
    fixed: np.ndarray              # selective-dynamics F F F
    bonds: set[tuple[int, int]]
    lengths: tuple[float, float]
    notes: dict = field(default_factory=dict)

    def counts(self) -> dict[str, int]:
        return {s: self.symbols.count(s) for s in SPECIES_ORDER if s in self.symbols}


# --------------------------------------------------------------------------- geometry helpers

def _mic(delta: np.ndarray, lengths) -> np.ndarray:
    out = np.array(delta, dtype=float, copy=True)
    for axis, length in enumerate(lengths):
        if length:
            out[..., axis] -= length * np.round(out[..., axis] / length)
    return out


def _rotation(axis, angle: float) -> np.ndarray:
    axis = np.asarray(axis, float); axis = axis / np.linalg.norm(axis)
    x, y, z = axis; c, s = math.cos(angle), math.sin(angle); C = 1 - c
    return np.array([[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
                     [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
                     [z * x * C - y * s, z * y * C + x * s, c + z * z * C]])


def _align(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector a onto unit vector b."""
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); c = float(np.dot(a, b))
    if s < 1e-12:
        if c > 0:
            return np.eye(3)
        perp = np.cross(a, [1.0, 0.0, 0.0] if abs(a[0]) < 0.9 else [0.0, 1.0, 0.0])
        return _rotation(perp, math.pi)
    return _rotation(v, math.atan2(s, c))


def _pair_distances(a: np.ndarray, b: np.ndarray, lengths) -> np.ndarray:
    return np.linalg.norm(_mic(a[:, None, :] - b[None, :, :], lengths), axis=2)


# --------------------------------------------------------------------------- substrate

def _sn_cluster(slab, sites: np.ndarray, trilayer: int, d: float) -> dict:
    """One neutral 2 Sn_In + O_i cluster with O_i on an empty 16c site of ``trilayer``.

    O_i is taken on the anion sublayer of that trilayer nearest the slab centre
    (away from the bottom face); the two Sn replace the two nearest In of the
    same trilayer, both of which must be bonded (< 2.6 A) to O_i.
    """
    lo, hi = trilayer * d, (trilayer + 1) * d
    inside = sites[(sites[:, 2] >= lo) & (sites[:, 2] < hi)]
    if not len(inside):
        raise ValueError(f"no empty 16c anion site inside trilayer {trilayer}")
    levels = np.unique(np.round(inside[:, 2], 3))
    centre = 0.5 * (slab.positions[:, 2].min() + slab.positions[:, 2].max())
    level = levels[np.argmin(np.abs(levels - centre))]
    cands = inside[np.abs(inside[:, 2] - level) < 1e-2]
    order = np.lexsort((np.round(cands[:, 1], 3), np.round(cands[:, 0], 3)))
    pos = slab.positions
    reasons = []
    for k in order:
        site = cands[k]
        in_idx = np.array([i for i, s in enumerate(slab.labels) if s == "In" and lo <= pos[i, 2] < hi])
        dist = np.linalg.norm(_mic(pos[in_idx] - site, slab.lengths), axis=1)
        near = in_idx[np.argsort(dist, kind="stable")[:2]]
        dn = np.sort(dist)[:2]
        if dn[1] > IN_O_BOND_CUTOFF:
            reasons.append(f"site {np.round(site, 3).tolist()}: second In at {dn[1]:.2f} A")
            continue
        o_idx = [i for i, s in enumerate(slab.labels) if s in ("O", "Oh")]
        do = np.linalg.norm(_mic(pos[o_idx] - site, slab.lengths), axis=1)
        for i in near:
            slab.labels[int(i)] = "Sn"
        slab.labels.append("O")
        slab.positions = np.vstack([pos, site])
        return {
            "trilayer": trilayer,
            "o_i_site_angstrom": [round(float(x), 6) for x in site],
            "o_i_index_in_builder": len(slab.labels) - 1,
            "sn_indices_in_builder": [int(i) for i in near],
            "sn_to_oi_angstrom": [round(float(x), 4) for x in dn],
            "oi_to_nearest_o_angstrom": round(float(do.min()), 4),
            "compensation": "neutral Frank-Koestlin-type 2 Sn_In + O_i on an empty 16c anion site (unrelaxed)",
        }
    raise ValueError("no 16c site in the requested trilayer has two same-trilayer In neighbours: " + "; ".join(reasons))


def hydroxylate_checked(slab, pairs: int, seed: int) -> dict:
    """Kept for callers; the clash filters now live in :func:`substrate.hydroxylate`."""
    from .substrate import hydroxylate
    return hydroxylate(slab, pairs, seed)


def build_substrate(settings: Settings, hydroxylated: bool, sn: bool) -> Case:
    T = settings.trilayers
    slab, sites = build_slab(1, 1, T)
    d = trilayer_spacing()
    n0 = len(slab.labels)
    origin = ["lattice"] * n0
    bonds: set[tuple[int, int]] = set()
    notes: dict = {"trilayers": T, "slab_builder_notes": slab.notes}
    if hydroxylated:
        pairs = settings.hydroxyl_pairs_per_primitive * 2  # the orthorhombic 1x1 cell = 2 primitive (111) cells
        rec = hydroxylate_checked(slab, pairs, settings.hydroxyl_seed)
        for k, site in enumerate(rec["sites"]):
            o_t, h_t, h_j = n0 + 3 * k, n0 + 3 * k + 1, n0 + 3 * k + 2
            bonds |= {(o_t, h_t), tuple(sorted((site["protonated_o_index"], h_j)))}
            origin += ["terminal-OH-O", "terminal-OH-H", "bridging-OH-H"]
        notes["hydroxylation"] = rec
    if sn:
        notes["sn_cluster"] = _sn_cluster(slab, sites, settings.sn_trilayer, d)
        origin.append("interstitial-O")
    pos = slab.positions
    trilayer = np.clip(np.floor(pos[:, 2] / d).astype(int), 0, T - 1)
    trilayer[[i for i, o in enumerate(origin) if o.startswith(("terminal", "bridging"))]] = T - 1
    fixed = trilayer == 0
    # Cation-anion bonds by the substrate coordination criterion (2.6 A).
    _, nbr = coordination(slab)
    for i, s in enumerate(slab.labels):
        if s in CATIONS:
            bonds |= {tuple(sorted((i, j))) for j in nbr[i]}
    if sn:
        cl = notes["sn_cluster"]
        core = [cl["o_i_index_in_builder"], *cl["sn_indices_in_builder"]]
        freed = sorted({*core, *(j for i in core for j in nbr[i])})
        freed_fixed = [i for i in freed if fixed[i]]
        fixed[freed] = False
        cl["freed_from_selective_dynamics"] = len(freed_fixed)
        cl["freed_rule"] = "O_i, both Sn, and their first coordination shell are relaxed although they sit in the fixed trilayer"
    symbols = ["O" if x == "Oh" else "H" if x == "Hh" else x for x in slab.labels]
    notes["origin"] = origin
    notes["trilayer_index"] = trilayer.tolist()
    return Case("", "slab", "", symbols, pos.copy(), list(slab.labels), [0] * len(symbols), fixed, bonds,
                slab.lengths, notes)


# --------------------------------------------------------------------------- molecule

def load_molecule(slug: str) -> dict:
    folder = ROOT / "inputs" / "molecules" / slug
    data = parse(folder / "ligpargen.lmp")
    symbols = elements(data)
    xyz = np.array(atom_coordinates(data), float)
    ids = [int(r.fields[0]) for r in data.sections["Atoms"]]
    index = {a: k for k, a in enumerate(ids)}
    bonds = {tuple(sorted((index[int(b.fields[2])], index[int(b.fields[3])]))) for b in data.sections.get("Bonds", [])}
    roles = phosphonate_roles(data)
    p = [index[a] for a, r in roles.items() if r == "P"]
    if len(p) != 1:
        raise ValueError(f"{slug}: expected exactly one phosphonic acid anchor")
    p = p[0]
    acid = []
    for a, r in roles.items():
        if r == "P-OH":
            o = index[a]
            h = [j for (i, j) in bonds if i == o and symbols[j] == "H"] + [i for (i, j) in bonds if j == o and symbols[i] == "H"]
            acid.append((o, h[0]))
    oxo = [index[a] for a, r in roles.items() if r == "P=O"][0]
    return {"slug": slug, "symbols": symbols, "positions": xyz, "ids": ids, "bonds": bonds, "p": p,
            "oxo": oxo, "acid": sorted(acid), "source": str((folder / "ligpargen.lmp").relative_to(ROOT)).replace("\\", "/"),
            "source_sha256": hashlib.sha256((folder / "ligpargen.lmp").read_bytes()).hexdigest()}


def orient_upright(mol: dict) -> dict:
    """P at the origin, P -> tail-centroid axis along +z, acidic O-H pointing down.

    The tail centroid is the mean position of all heavy atoms outside the
    PO3H2 group.  Each acidic H is rotated about its P-O axis (5 degree grid)
    to the lowest z that keeps every intramolecular nonbonded contact >= 1.8 A.
    """
    x = mol["positions"] - mol["positions"][mol["p"]]
    s = mol["symbols"]
    head = {mol["p"], mol["oxo"], *(o for o, _ in mol["acid"]), *(h for _, h in mol["acid"])}
    tail = [i for i, e in enumerate(s) if e != "H" and i not in head]
    axis = x[tail].mean(axis=0)
    x = x @ _align(axis, np.array([0.0, 0.0, 1.0])).T
    heavy = [i for i, e in enumerate(s) if e != "H"]
    xy = x[heavy, :2] - x[heavy, :2].mean(axis=0)
    w, v = np.linalg.eigh(xy.T @ xy)
    main = v[:, np.argmax(w)]
    n_idx = [i for i, e in enumerate(s) if e == "N"]
    ref = x[n_idx[0], :2] if n_idx else x[tail[-1], :2]
    if float(np.dot(main, ref)) < 0:
        main = -main
    angle = math.atan2(main[0], main[1])  # rotate main in-plane axis onto +y
    x = x @ _rotation([0, 0, 1], angle).T
    for o, h in mol["acid"]:
        ax = x[o] - x[mol["p"]]
        best = None
        others = [k for k in range(len(s)) if k not in (o, h)]
        for step in range(72):
            R = _rotation(ax, math.radians(5 * step))
            cand = x[o] + R @ (x[h] - x[o])
            dmin = float(np.min(np.linalg.norm(x[others] - cand, axis=1)))
            if dmin >= 1.8 and (best is None or cand[2] < best[0] - 1e-9):
                best = (cand[2], cand, step)
        if best is None:
            raise ValueError(f"{mol['slug']}: no clash-free orientation for acidic H {h}")
        x[h] = best[1]
    out = dict(mol); out["positions"] = x
    out["orientation"] = {"axis": "P -> centroid of non-PO3H2 heavy atoms along +z",
                          "in_plane": "principal heavy-atom axis along +y (N-side positive)",
                          "acid_h": "rotated about P-O to minimum z, intramolecular nonbonded >= 1.8 A"}
    return out


# --------------------------------------------------------------------------- adsorption geometry

def _top_sites(case: Case) -> dict:
    slab_idx = [i for i, l in enumerate(case.labels) if l != "mol"]
    pos = case.positions
    z = pos[slab_idx, 2]; z_mid = 0.5 * (z.min() + z.max())
    origin = case.notes["origin"]
    cn = np.zeros(len(case.symbols), int)
    for i, j in case.bonds:
        if case.symbols[i] in CATIONS or case.symbols[j] in CATIONS:
            cn[i] += 1; cn[j] += 1
    in5c = [i for i in slab_idx if case.labels[i] in CATIONS and cn[i] == 5 and pos[i, 2] > z_mid]
    o3c = [i for i in slab_idx if case.labels[i] == "O" and origin[i] == "lattice" and cn[i] == 3 and pos[i, 2] > z_mid]
    lattice_o = [i for i in slab_idx if origin[i] == "lattice" and case.symbols[i] == "O"]
    z_ref = float(max(pos[i, 2] for i in lattice_o))
    return {"in5c": in5c, "o3c": o3c, "z_top_lattice_o": z_ref, "terminal_o": [i for i in slab_idx if origin[i] == "terminal-OH-O"]}


def _open_site(case: Case, i: int, length: float) -> np.ndarray:
    """Position of the missing sixth O of a surface In5c at ``length`` from it."""
    nb = [b for a, b in case.bonds if a == i] + [a for a, b in case.bonds if b == i]
    nb = [k for k in nb if case.symbols[k] == "O"]
    v = _mic(case.positions[nb] - case.positions[i], case.lengths[:2] + (0.0,))
    vac = -np.sum(v / np.linalg.norm(v, axis=1)[:, None], axis=0)
    return case.positions[i] + length * vac / np.linalg.norm(vac)


def choose_in_pair(case: Case, bond: float = 2.15) -> dict:
    """Adjacent free In5c pair for the bidentate anchor (and the physisorption site above it).

    Candidates are free top In5c whose open octahedral site (at ``bond``) is
    >= 2.5 A from every slab O and >= 2.0 A from every H; among pairs 3.0-4.0 A
    apart, the one whose two open sites are closest to the phosphonate O...O
    distance (2.5 A) wins, then the one farthest from terminal OH.
    """
    t = _top_sites(case)
    pos = case.positions; L = case.lengths; L3 = L[:2] + (0.0,)
    o_all = [i for i, s in enumerate(case.symbols) if s == "O"]
    h_all = [i for i, s in enumerate(case.symbols) if s == "H"]
    sites, rejected = {}, []
    for i in t["in5c"]:
        o = _open_site(case, i, bond)
        d_o = float(np.min(np.linalg.norm(_mic(pos[o_all] - o, L3), axis=1)))
        d_h = float(np.min(np.linalg.norm(_mic(pos[h_all] - o, L3), axis=1))) if h_all else 99.0
        if d_o >= 2.5 and d_h >= 2.0:
            sites[i] = o
        else:
            rejected.append(int(i))
    best = None
    keys = sorted(sites)
    for a_k, a in enumerate(keys):
        for b in keys[a_k + 1:]:
            dv = _mic(pos[b] - pos[a], L3); dab = float(np.linalg.norm(dv))
            if not 3.0 <= dab <= 4.0:
                continue
            sep = float(np.linalg.norm(_mic(sites[b] - sites[a], L3)))
            mid = pos[a] + 0.5 * dv
            room = float(np.min(np.linalg.norm(_mic(pos[t["terminal_o"], :2] - mid[:2], L[:2]), axis=1))) if t["terminal_o"] else 99.0
            key = (round(abs(sep - 2.5), 3), -round(room, 3), a, b)
            if best is None or key < best[0]:
                best = (key, {"in_indices": [int(a), int(b)], "in_in_distance": round(dab, 4),
                              "open_site_separation": round(sep, 4),
                              "open_sites": [[round(float(x), 6) for x in sites[a]], [round(float(x), 6) for x in sites[b]]],
                              "midpoint_xy": [round(float(mid[0]), 6), round(float(mid[1]), 6)],
                              "midpoint_z": round(float(mid[2]), 6),
                              "lateral_distance_to_nearest_terminal_oh": round(room, 4)})
    if best is None:
        raise ValueError("no adjacent free In5c pair (3.0-4.0 A) with clash-free open sites on the top face")
    return best[1] | {"z_top_lattice_o": t["z_top_lattice_o"], "in5c_rejected_crowded_open_site": rejected,
                      "selection": "free In5c with open site >= 2.5 A from O and >= 2.0 A from H; pair 3.0-4.0 A; "
                                   "open-site separation closest to 2.5 A"}


def _image_min(x: np.ndarray, L) -> float:
    out = 99.0
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i or j:
                shift = np.array([i * L[0], j * L[1], 0.0])
                out = min(out, float(np.min(np.linalg.norm(x[:, None, :] - (x + shift)[None, :, :], axis=2))))
    return out


def _clash_thresholds(m_sym, s_sym) -> np.ndarray:
    """Soft floors (A) for molecule-slab contacts used by the placement scores."""
    m_h = np.array([s == "H" for s in m_sym]); m_o = np.array([s == "O" for s in m_sym])
    s_h = np.array([s == "H" for s in s_sym]); s_cat = np.array([s in CATIONS for s in s_sym])
    thr = np.full((len(m_sym), len(s_sym)), 2.6)
    thr[m_h, :] = 1.9; thr[:, s_h] = 1.9
    thr[np.ix_(m_o, s_cat)] = 3.0
    return thr


def _intra_pairs(n: int, bonds) -> np.ndarray:
    """Mask of intramolecular pairs separated by more than two bonds."""
    adj = [set() for _ in range(n)]
    for i, j in bonds:
        adj[i].add(j); adj[j].add(i)
    mask = np.ones((n, n), bool)
    np.fill_diagonal(mask, False)
    for i in range(n):
        for j in adj[i]:
            mask[i, j] = False
            for k in adj[j]:
                mask[i, k] = False
    return mask


def _orient_acid_h(x: np.ndarray, mol: dict, spos: np.ndarray, s_sym, L3) -> list[dict]:
    """Turn each acidic H about its P-O bond toward the nearest surface O (in place).

    Floors: 1.8 A intramolecular, 1.6 A to slab heavy atoms, 2.0 A to slab H.
    """
    s_o = np.array([s == "O" for s in s_sym]); s_h = np.array([s == "H" for s in s_sym])
    out = []
    for o, h in mol["acid"]:
        ax = x[o] - x[mol["p"]]
        others = [k for k in range(len(x)) if k not in (o, h)]
        best = None
        for step in range(72):
            cand = x[o] + _rotation(ax, math.radians(5 * step)) @ (x[h] - x[o])
            d_mol = float(np.min(np.linalg.norm(x[others] - cand, axis=1)))
            d_s = np.linalg.norm(_mic(spos - cand, L3), axis=1)
            d_acc = float(np.min(d_s[s_o])); d_any = float(np.min(d_s[~s_h]))
            d_hh = float(np.min(d_s[s_h])) if s_h.any() else 99.0
            ok = d_mol >= 1.8 and d_any >= 1.6 and d_hh >= 2.0
            key = (not ok, round(d_acc, 4), step)
            if best is None or key < best[0]:
                best = (key, cand, d_acc, ok)
        x[h] = best[1]
        out.append({"lig_id": mol["ids"][h], "h_to_nearest_surface_o": round(best[2], 4), "clash_free": best[3]})
    return out


def place_phys(slab: Case, mol: dict, pair: dict, height: float) -> tuple[np.ndarray, dict]:
    """Rigid scan of the upright acid for an H-bonded physisorption start.

    Grid: axis tilt (0, or 15 deg toward four azimuths), z-rotation (10 deg),
    lateral P offset around the In-pair midpoint (+-2.5 A, 0.5 A) and P
    height ``height`` +- 0.25 A above the top lattice-O plane; only slab atoms
    within 2.5 A below that plane are scored.  Score (lower is better): each P-OH oxygen should have a
    surface O within 2.7 A (quadratic beyond; the H-bond O...O distance), plus
    quadratic clash penalties (acidic H excluded; P-OH O...surface O floor
    2.5 A).  The acidic H are then turned about their P-O bonds toward the
    nearest surface O acceptor.
    """
    L = slab.lengths; L3 = L[:2] + (0.0,)
    spos_all = slab.positions
    near = np.where(spos_all[:, 2] >= pair["z_top_lattice_o"] - 2.5)[0]
    spos = spos_all[near]
    s_o = np.array([slab.symbols[i] == "O" for i in near])
    acid_h = [h for _, h in mol["acid"]]; acid_o = [o for o, _ in mol["acid"]]
    use = [i for i in range(len(mol["symbols"])) if i not in acid_h]
    base = mol["positions"]
    thr = _clash_thresholds([mol["symbols"][i] for i in use], [slab.symbols[i] for i in near])
    loc = {i: k for k, i in enumerate(use)}
    ao = [loc[o] for o in acid_o]
    thr[np.ix_(ao, np.where(s_o)[0])] = 2.5
    grid = np.arange(-2.5, 2.501, 0.5)
    off = np.array([(dx, dy) for dx in grid for dy in grid])
    heights = (height - 0.25, height, height + 0.25)
    tilts = [(0, 0)] + [(15, az) for az in (0, 90, 180, 270)]
    best = None
    for tilt, az in tilts:
        r_tilt = _rotation([math.cos(math.radians(az)), math.sin(math.radians(az)), 0.0], math.radians(tilt))
        for theta in range(0, 360, 10):
            rot = base @ (_rotation([0, 0, 1], math.radians(theta)) @ r_tilt).T
            img = _image_min(rot, L)
            if img < 3.0:
                continue
            for hz in heights:
                x0 = rot + np.array([pair["midpoint_xy"][0], pair["midpoint_xy"][1], pair["z_top_lattice_o"] + hz])
                delta = np.repeat((x0[use][:, None, :] - spos[None, :, :])[None], len(off), axis=0)
                for axis in (0, 1):
                    delta[..., axis] += off[:, None, None, axis]
                    delta[..., axis] -= L[axis] * np.round(delta[..., axis] / L[axis])
                d = np.sqrt(np.einsum("omsk,omsk->oms", delta, delta))        # (offsets, mol, slab)
                clash = np.sum(np.clip(thr[None] - d, 0.0, None) ** 2, axis=(1, 2))
                doo = np.min(np.where(s_o[None, None, :], d[:, ao, :], np.inf), axis=2)
                hb = np.sum(np.clip(doo - 2.7, 0.0, None) ** 2, axis=1)
                score = np.round(hb + 10.0 * clash, 6)
                k = int(np.lexsort((np.abs(off[:, 0]) + np.abs(off[:, 1]), score))[0])
                key = (float(score[k]), tilt, round(float(abs(off[k, 0]) + abs(off[k, 1])), 3), abs(hz - height), theta, az)
                if best is None or key < best[0]:
                    x = x0 + np.array([off[k, 0], off[k, 1], 0.0])
                    best = (key, x, {"tilt_deg": tilt, "tilt_azimuth_deg": az, "theta_deg": theta,
                                     "offset_xy": [round(float(off[k, 0]), 3), round(float(off[k, 1]), 3)],
                                     "p_height_above_top_lattice_o": hz, "score": float(score[k]),
                                     "p_oh_o_to_nearest_surface_o": [round(float(v), 4) for v in doo[k]],
                                     "min_mol_image": round(img, 4)})
    if best is None:
        raise ValueError("no rotation keeps the molecule 3 A from its lateral images")
    x, info = best[1].copy(), best[2]
    info["acid_h"] = _orient_acid_h(x, mol, spos_all, slab.symbols, L3)
    info["min_mol_slab"] = round(float(np.min(_pair_distances(x, spos_all, L3))), 4)
    info["hbond_contacts"] = _hbond_contacts(slab, mol, x)
    info["method"] = ("rigid upright LigParGen conformer; tilt/z-rotation/offset/height grid scored on P-OH O...O(surface) <= 2.7 A "
                      "and clash floors; acidic H then turned toward the nearest surface O")
    return x, info


def _hbond_contacts(slab: Case, mol: dict, x: np.ndarray) -> list[dict]:
    L = slab.lengths
    o_slab = [i for i, s in enumerate(slab.symbols) if s == "O"]
    out = []
    head_o = [mol["oxo"], *(o for o, _ in mol["acid"])]
    for o in head_o:
        d = np.linalg.norm(_mic(slab.positions[o_slab] - x[o], L[:2] + (0.0,)), axis=1)
        k = int(np.argmin(d))
        role = "P=O" if o == mol["oxo"] else "P-OH"
        out.append({"lig_id": mol["ids"][o], "role": role, "nearest_surface_o_index": int(o_slab[k]),
                    "surface_o_origin": slab.notes["origin"][o_slab[k]], "o_o_distance": round(float(d[k]), 4)})
    for o, h in mol["acid"]:
        d = np.linalg.norm(_mic(slab.positions[o_slab] - x[h], L[:2] + (0.0,)), axis=1)
        k = int(np.argmin(d))
        out.append({"lig_id": mol["ids"][h], "role": "P-O-H", "nearest_surface_o_index": int(o_slab[k]),
                    "surface_o_origin": slab.notes["origin"][o_slab[k]], "h_o_distance": round(float(d[k]), 4)})
    return out


def _rotatable_path(n: int, bonds, start: int, symbols) -> list[tuple[int, int, list[int]]]:
    """Acyclic bonds on the shortest path from ``start`` (P) to the first ring atom, with their distal sides."""
    adj = [set() for _ in range(n)]
    for i, j in bonds:
        adj[i].add(j); adj[j].add(i)

    def side(a, b):
        seen, stack = {b}, [b]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in seen and not (u == b and v == a):
                    seen.add(v); stack.append(v)
        return seen

    def in_ring(a, b):
        return a in side(a, b)

    prev, queue = {start: None}, [start]
    while queue:
        u = queue.pop(0)
        for v in sorted(adj[u]):
            if v not in prev:
                prev[v] = u; queue.append(v)
    ring_atoms = sorted({a for a, b in bonds if in_ring(a, b)} | {b for a, b in bonds if in_ring(a, b)})
    target = min(ring_atoms, key=lambda r: (len(_path(prev, r)), r))
    path = _path(prev, target)
    out = []
    for a, b in zip(path, path[1:]):
        if not in_ring(a, b) and symbols[a] != "H" and symbols[b] != "H":
            out.append((a, b, sorted(side(a, b))))
    return out


def _path(prev: dict, node: int) -> list[int]:
    out = []
    while node is not None:
        out.append(node); node = prev[node]
    return out[::-1]


def place_bidentate(slab: Case, mol: dict, pair: dict, bond: float) -> tuple[np.ndarray, list[int], dict]:
    """Bridging-bidentate start for the doubly deprotonated acid.

    1. The two former P-OH oxygens are put on the open octahedral sites of the
       chosen In5c pair (their midpoint on the sites' midpoint, O...O along the
       site-site vector), and the rigid molecule is turned about that axis
       (5 deg grid) to minimise head-group clashes and keep P-C pointing up.
    2. The tail is then relaxed by a greedy torsion drive (10 deg grid, two
       passes) about the acyclic bonds from P to the first ring atom, scoring
       molecule-slab clashes, intramolecular clashes (>2 bonds apart, floor
       2.0 A), lateral-image contacts (< 3 A) and a reward for tail height.
    Bond lengths and angles are those of the LigParGen geometry.
    """
    L = slab.lengths; L3 = L[:2] + (0.0,)
    spos = slab.positions
    acid_h = {h for _, h in mol["acid"]}
    keep = [i for i in range(len(mol["symbols"])) if i not in acid_h]
    kl = {i: k for k, i in enumerate(keep)}
    x0 = mol["positions"][keep]
    sym = [mol["symbols"][i] for i in keep]
    bonds = [(kl[i], kl[j]) for i, j in mol["bonds"] if i in kl and j in kl]
    thr = _clash_thresholds(sym, slab.symbols)
    ina, inb = pair["in_indices"]
    for o, _ in mol["acid"]:
        thr[kl[o], [ina, inb]] = 0.0  # bonded partners handled by the site targets
    intra = _intra_pairs(len(sym), bonds)
    p = kl[mol["p"]]
    c_alpha = [j for i, j in bonds if i == p and sym[j] == "C"] + [i for i, j in bonds if j == p and sym[i] == "C"]
    c_alpha = c_alpha[0]
    head = [p, kl[mol["oxo"]], c_alpha]
    tail = [k for k, s in enumerate(sym) if s != "H" and k not in head and k not in [kl[o] for o, _ in mol["acid"]]]
    s1, s2 = (np.array(s) for s in pair["open_sites"])
    w = _mic(s2 - s1, L3); ms = s1 + 0.5 * w; axis = w / np.linalg.norm(w)
    rot_bonds = _rotatable_path(len(sym), bonds, p, sym)

    def score(x):
        d = _pair_distances(x, spos, L3)
        f = 10.0 * float(np.sum(np.clip(thr - d, 0.0, None) ** 2))
        di = np.linalg.norm(x[:, None] - x[None], axis=2)
        f += 10.0 * float(np.sum(np.clip(2.0 - di[intra], 0.0, None) ** 2))
        f += 10.0 * max(0.0, 3.0 - _image_min(x, L)) ** 2
        return f - 0.2 * float(np.mean(x[tail, 2]) - ms[2])

    best = None
    for perm in ((kl[mol["acid"][0][0]], kl[mol["acid"][1][0]]), (kl[mol["acid"][1][0]], kl[mol["acid"][0][0]])):
        o1, o2 = perm
        u = x0[o2] - x0[o1]
        y = (x0 - 0.5 * (x0[o1] + x0[o2])) @ _align(u, w).T + ms
        for step in range(72):
            z = (y - ms) @ _rotation(axis, math.radians(5 * step)).T + ms
            pc = z[c_alpha] - z[p]
            if z[p, 2] <= ms[2]:
                continue
            d = _pair_distances(z[head], spos, L3)
            f = 10.0 * float(np.sum(np.clip(thr[head] - d, 0.0, None) ** 2)) + (1.0 - pc[2] / np.linalg.norm(pc))
            key = (round(f, 6), perm, step)
            if best is None or key < best[0]:
                best = (key, z, perm, step)
    _, x, perm, phi = best
    x = x.copy()
    torsions = []
    for _pass in range(2):
        for a, b, dist in rot_bonds:
            ax = x[b] - x[a]
            cur = score(x); pick = (round(cur, 6), 0)
            trial = None
            for step in range(1, 36):
                y = x.copy()
                y[dist] = (x[dist] - x[b]) @ _rotation(ax, math.radians(10 * step)).T + x[b]
                sc = round(score(y), 6)
                if (sc, step) < pick:
                    pick = (sc, step); trial = y
            if trial is not None:
                x = trial
            torsions.append({"bond_lig_ids": [mol["ids"][keep[a]], mol["ids"][keep[b]]], "pass": _pass + 1, "rotation_deg": 10 * pick[1]})
    d = _pair_distances(x, spos, L3)
    tail_axis = x[tail].mean(axis=0) - x[p]
    info = {
        "score": round(float(score(x)), 6),
        "in_o_bonds": [{"lig_id": mol["ids"][keep[o]], "in_index": int(c), "distance": round(float(d[o, c]), 4)}
                       for o, c in zip(perm, (ina, inb))],
        "head_rotation_deg": 5 * phi,
        "p_c_tilt_deg": round(math.degrees(math.acos(float((x[c_alpha] - x[p])[2] / np.linalg.norm(x[c_alpha] - x[p])))), 2),
        "tail_axis_tilt_deg": round(math.degrees(math.acos(float(tail_axis[2] / np.linalg.norm(tail_axis)))), 2),
        "torsion_drive": torsions,
        "min_mol_image": round(_image_min(x, L), 4),
        "method": "P-O(former OH) oxygens on the In5c open octahedral sites, head rotation about the O...O axis, "
                  "greedy torsion drive of the P-to-ring acyclic bonds; LigParGen bond lengths/angles kept",
    }
    return x, keep, info | {"_perm": list(perm)}


def _o3c_proton(case: Case, j: int) -> np.ndarray:
    pos = case.positions; L = case.lengths
    nb = [b for a, b in case.bonds if a == j] + [a for a, b in case.bonds if b == j]
    nb = [k for k in nb if case.symbols[k] in CATIONS]
    v = _mic(pos[nb] - pos[j], L[:2] + (0.0,))
    out = -np.sum(v / np.linalg.norm(v, axis=1)[:, None], axis=0)
    z_hat = np.array([0.0, 0.0, 1.0])
    out = out / np.linalg.norm(out) if np.linalg.norm(out) > 1e-6 else z_hat
    if out[2] < 0.2:
        out = out + z_hat; out /= np.linalg.norm(out)
    return pos[j] + O_H_LENGTH * out


# --------------------------------------------------------------------------- case assembly

def _combine(slab: Case, mol: dict, x: np.ndarray, keep: list[int] | None = None) -> Case:
    keep = list(range(len(mol["symbols"]))) if keep is None else keep
    n = len(slab.symbols)
    loc = {i: n + k for k, i in enumerate(keep)}
    bonds = set(slab.bonds) | {tuple(sorted((loc[i], loc[j]))) for i, j in mol["bonds"] if i in loc and j in loc}
    notes = dict(slab.notes)
    notes["origin"] = slab.notes["origin"] + ["molecule"] * len(keep)
    notes["trilayer_index"] = slab.notes["trilayer_index"] + [-1] * len(keep)
    return Case("", "adsorbate", "", slab.symbols + [mol["symbols"][i] for i in keep], np.vstack([slab.positions, x]),
                slab.labels + ["mol"] * len(keep), slab.lig_id + [mol["ids"][i] for i in keep],
                np.concatenate([slab.fixed, np.zeros(len(keep), bool)]), bonds, slab.lengths, notes)


def bidentate_case(slab: Case, mol: dict, pair: dict, bond: float) -> Case:
    x, keep, info = place_bidentate(slab, mol, pair, bond)
    perm = info.pop("_perm")
    case = _combine(slab, mol, x, keep)
    n = len(slab.symbols)
    kloc = {i: k for k, i in enumerate(keep)}
    inv = {kloc[o]: o for o, _ in mol["acid"]}
    for o_k, cat in zip(perm, pair["in_indices"]):
        case.bonds.add(tuple(sorted((n + o_k, cat))))
    # Transfer each acidic H to the nearest free O3c (distinct, clash-free).
    t = _top_sites(slab)
    used = set()
    transfers = []
    for o_k in perm:
        o_orig = inv[o_k]
        h_orig = dict(mol["acid"])[o_orig]
        cands = [j for j in t["o3c"] if j not in used]
        d = np.linalg.norm(_mic(slab.positions[cands] - case.positions[n + o_k], slab.lengths[:2] + (0.0,)), axis=1)
        placed = None
        for k in np.argsort(d, kind="stable"):
            j = cands[int(k)]
            h = _o3c_proton(case, j)
            dd = np.linalg.norm(_mic(case.positions - h, slab.lengths[:2] + (0.0,)), axis=1)
            dd[j] = 99.0
            if dd.min() >= 1.6:
                placed = (j, h, float(d[k]), float(dd.min()))
                break
        if placed is None:
            raise ValueError("no clash-free O3c found for proton transfer")
        j, h, d_oo, d_h = placed
        used.add(j)
        case.symbols.append("H"); case.labels.append("Hh"); case.lig_id.append(0)
        case.positions = np.vstack([case.positions, h]); case.fixed = np.append(case.fixed, False)
        case.bonds.add((j, len(case.symbols) - 1))
        case.labels[j] = "Oh"
        case.notes["origin"] = case.notes["origin"] + ["transferred-H"]
        case.notes["trilayer_index"] = case.notes["trilayer_index"] + [slab.notes["trilayer_index"][j]]
        transfers.append({"from_lig_h_id": mol["ids"][h_orig], "donor_lig_o_id": mol["ids"][o_orig],
                          "acceptor_o3c_index": int(j), "acceptor_to_donor_o_distance": round(d_oo, 4),
                          "min_distance_new_h": round(d_h, 4)})
    info["proton_transfers"] = transfers
    case.notes["bidentate"] = info
    return case


def molecule_case(mol: dict, x: np.ndarray, lengths) -> Case:
    n = len(mol["symbols"])
    return Case("", "molecule", "", list(mol["symbols"]), x.copy(), ["mol"] * n, list(mol["ids"]), np.zeros(n, bool),
                set(mol["bonds"]), lengths, {"origin": ["molecule"] * n, "trilayer_index": [-1] * n})


PURPOSES = {
    "slab-bare": "Stoichiometric In2O3(111) reference; surface relaxation scale and dry-surface term for any bare-slab classical check.",
    "slab-oh": "Hydroxylated In2O3(111) reference (6 dissociated H2O / 1x1 cell); the E_slab term of every E_ads.",
    "slab-oh-sn": "slab-oh with one neutral 2Sn_In + O_i cluster; Sn effect on the surface (compare with slab-oh, and later with Sn-doped adsorption).",
    "mol": "Isolated protonated acid in the same cell (rigid copy of the physisorbed start geometry); the E_mol term of every E_ads.",
    "phys": "Upright physisorbed acid H-bonded to slab-oh: the configuration class the classical model represents; E_ads(phys) DFT vs classical.",
    "bidentate": "Bridging-bidentate P-O-In chemisorption with both acidic H moved to surface O3c (same atoms as phys): E(bid)-E(phys) is the chemisorption energy the classical model cannot capture.",
}


def build_cases(settings: Settings | None = None) -> list[Case]:
    settings = settings or Settings()
    bare = build_substrate(settings, False, False)
    oh = build_substrate(settings, True, False)
    ohsn = build_substrate(settings, True, True)
    bare.name, bare.purpose = "slab-bare", PURPOSES["slab-bare"]
    oh.name, oh.purpose = "slab-oh", PURPOSES["slab-oh"]
    ohsn.name, ohsn.purpose = "slab-oh-sn", PURPOSES["slab-oh-sn"]
    cases = [bare, oh, ohsn]
    pair = choose_in_pair(oh, settings.bidentate_in_o)
    for slug in settings.molecules:
        mol = orient_upright(load_molecule(slug))
        x, info = place_phys(oh, mol, pair, settings.phys_height)
        m = molecule_case(mol, x, oh.lengths)
        m.name, m.purpose = f"mol-{slug}", PURPOSES["mol"]
        m.notes["molecule"] = {"slug": slug, "source": mol["source"], "source_sha256": mol["source_sha256"],
                               "orientation": mol["orientation"],
                               "geometry": f"rigid z-translation of the ads-{slug}-phys molecule"}
        phys = _combine(oh, mol, x)
        phys.name, phys.purpose = f"ads-{slug}-phys", PURPOSES["phys"]
        phys.notes["physisorption"] = info | {"in_pair": pair}
        phys.notes["molecule"] = m.notes["molecule"] | {"geometry": "upright rigid conformer; acidic H turned toward acceptors"}
        bid = bidentate_case(oh, mol, pair, settings.bidentate_in_o)
        bid.name, bid.purpose = f"ads-{slug}-bidentate", PURPOSES["bidentate"]
        bid.notes["bidentate"]["in_pair"] = pair
        bid.notes["molecule"] = phys.notes["molecule"] | {"geometry": "doubly deprotonated; head on In pair, tail torsion drive"}
        cases += [m, phys, bid]
    _finalise_cell(cases, settings)
    for c in cases:
        c.notes["distances"] = distance_report(c)
        if c.notes["distances"]["min_nonbonded"] <= MIN_NONBONDED:
            raise ValueError(f"{c.name}: nonbonded contact {c.notes['distances']['min_nonbonded']:.3f} A <= {MIN_NONBONDED} A")
    return cases


def _finalise_cell(cases: list[Case], s: Settings) -> None:
    slab_bottom = min(float(c.positions[:, 2].min()) for c in cases if c.kind != "molecule")
    need = []
    for c in cases:
        span = float(np.ptp(c.positions[:, 2]))
        need.append(span + (s.vacuum_slab if c.kind == "slab" else s.vacuum_molecule))
    def roundup(v): return math.ceil(v * 2.0) / 2.0
    common = roundup(max(need))
    for c, n in zip(cases, need):
        cz = common if s.common_cell else roundup(n)
        if c.kind == "molecule":
            z = c.positions[:, 2]
            c.positions[:, 2] += 0.5 * cz - 0.5 * (z.min() + z.max())
        else:
            c.positions[:, 2] += 1.0 - slab_bottom
        c.positions[:, 0] %= c.lengths[0]; c.positions[:, 1] %= c.lengths[1]
        z = c.positions[:, 2]
        c.notes["cell_c"] = cz
        c.notes["vacuum_gap"] = round(cz - float(np.ptp(z)), 4)
        c.notes["dipol_centre_frac"] = round(0.5 * (float(z.min()) + float(z.max())) / cz, 6)


# --------------------------------------------------------------------------- checks

def distance_report(c: Case) -> dict:
    L = (c.lengths[0], c.lengths[1], c.notes["cell_c"])
    d = _pair_distances(c.positions, c.positions, L)
    n = len(c.symbols)
    np.fill_diagonal(d, np.inf)
    bonded = np.zeros((n, n), bool)
    for i, j in c.bonds:
        bonded[i, j] = bonded[j, i] = True
    nb = np.where(bonded, np.inf, d)
    i, j = np.unravel_index(np.argmin(nb), nb.shape)
    out = {"min_nonbonded": round(float(nb[i, j]), 4),
           "min_nonbonded_pair": [c.symbols[i] + str(i + 1), c.symbols[j] + str(j + 1)],
           "bonded_range": [round(float(d[bonded].min()), 4), round(float(d[bonded].max()), 4)] if bonded.any() else None,
           "bonded_pairs": int(bonded.sum() // 2)}
    mol = np.array([l == "mol" for l in c.labels])
    if mol.any() and (~mol).any():
        sub = nb[np.ix_(mol, ~mol)]
        out["min_molecule_slab_nonbonded"] = round(float(sub.min()), 4)
    sl = ~mol
    if sl.sum() > 1:
        out["min_slab_nonbonded"] = round(float(nb[np.ix_(sl, sl)].min()), 4)
    h = np.array([s == "H" for s in c.symbols])
    if h.any():
        out["min_h_nonbonded"] = round(float(nb[h].min()), 4)
    return out


def nelect(counts: dict[str, int]) -> int:
    return int(sum(n * ZVAL[POTCAR_DATASETS[e]] for e, n in counts.items()))


# --------------------------------------------------------------------------- writers

def _order(c: Case) -> list[int]:
    return sorted(range(len(c.symbols)), key=lambda i: (SPECIES_ORDER.index(c.symbols[i]), i))


def poscar_text(c: Case) -> str:
    idx = _order(c)
    species = [s for s in SPECIES_ORDER if s in c.symbols]
    lx, ly = c.lengths; lz = c.notes["cell_c"]
    lines = [f"{c.name} " + " ".join(species), "1.0",
             f"  {lx:.10f}  0.0000000000  0.0000000000", f"  0.0000000000  {ly:.10f}  0.0000000000",
             f"  0.0000000000  0.0000000000  {lz:.10f}",
             "  " + "  ".join(species), "  " + "  ".join(str(c.symbols.count(s)) for s in species)]
    sd = c.kind != "molecule"
    if sd:
        lines.append("Selective dynamics")
    lines.append("Cartesian")
    for i in idx:
        p = c.positions[i]
        row = f"  {p[0]:.8f}  {p[1]:.8f}  {p[2]:.8f}"
        if sd:
            row += "   F   F   F" if c.fixed[i] else "   T   T   T"
        lines.append(row)
    return "\n".join(lines) + "\n"


def extxyz_text(c: Case) -> str:
    idx = _order(c)
    lx, ly = c.lengths; lz = c.notes["cell_c"]
    head = (f'Lattice="{lx:.10f} 0.0 0.0 0.0 {ly:.10f} 0.0 0.0 0.0 {lz:.10f}" '
            'Properties=species:S:1:pos:R:3:label:S:1:lig_id:I:1:sd_fixed:L:1 '
            f'case={c.name} pbc="T T T"')
    rows = [f"{c.symbols[i]} {c.positions[i,0]:.8f} {c.positions[i,1]:.8f} {c.positions[i,2]:.8f} "
            f"{c.labels[i]} {c.lig_id[i]} {'T' if c.fixed[i] else 'F'}" for i in idx]
    return f"{len(idx)}\n{head}\n" + "\n".join(rows) + "\n"


def incar_text(c: Case, s: Settings) -> str:
    prod = s.mode == "production"
    nsw = 300 if prod else 3
    lines = [
        f"# ITO/SAM DFT reference point: {c.name} ({s.mode} mode)",
        "# Generated by nio_md_prep.ito.vasp; tag set follows the InterfaceForge NiO(110)+phosphonate",
        f"# relaxation template (InterfaceForge@{INTERFACEFORGE_SHA[:12]}); deviations are listed in docs/ito/vasp-smoke.md",
        f"SYSTEM = {c.name}",
        "",
        "# Electronic structure",
        "GGA     = PE",
        f"ENCUT   = {s.encut:.0f}",
        "PREC    = Accurate",
        "EDIFF   = 1E-5",
        "NELM    = 150",
        "NELMIN  = 6",
        "ALGO    = Normal",
        "LREAL   = Auto",
        "LASPH   = .TRUE.",
        "",
        "# Occupations",
        "ISMEAR = 0",
        "SIGMA  = 0.05",
        "",
        "# Closed shell: In2O3, ionically compensated 2Sn_In+O_i, and the neutral acids are nonmagnetic",
        "ISPIN = 1",
        "ISYM  = 0",
        "",
        "# No DFT+U: In 4d is a filled semicore shell (In_d PAW), not a correlated open shell",
        "# Grimme D3 with Becke-Johnson damping",
        "IVDW = 12",
        "",
        "# Fixed-cell ionic relaxation" + ("" if prod else " (smoke: 3 ionic steps; step 1 is a single point on geometry.extxyz)"),
        "IBRION = 2",
        "ISIF   = 2",
        f"NSW    = {nsw}",
        "EDIFFG = -0.03",
        "POTIM  = 0.3",
        "",
        "# Dipole correction along the normal; DIPOL puts the correction plane at the vacuum midpoint",
        "LDIPOL = .TRUE.",
        "IDIPOL = 3",
        f"DIPOL  = 0.5 0.5 {c.notes['dipol_centre_frac']:.6f}",
        "AMIN   = 0.01",
        "",
        "# Output",
        f"LWAVE  = {'.TRUE.' if prod else '.FALSE.'}",
        f"LCHARG = {'.TRUE.' if prod else '.FALSE.'}",
        "NWRITE = 1",
    ]
    if prod:
        lines.append("LORBIT = 11")
    lines += ["", "# Parallelization", "NCORE = 4"]
    return "\n".join(lines) + "\n"


KPOINTS_TEXT = "Automatic mesh   !  Gamma only (in-plane cell >= 14 A)\n0\nG\n 1 1 1\n"


def potcar_spec_text(c: Case) -> str:
    counts = c.counts()
    lines = ["# POTCAR.spec: PBE PAW datasets in POSCAR species order (potpaw_PBE.54 names as in",
             f"# InterfaceForge@{INTERFACEFORGE_SHA[:12]} profiles/potcar_pbe_54.yaml). The licensed POTCAR is",
             "# assembled on the cluster by run.sbatch; it is never stored in git.",
             "# element  dataset  count  zval_expected"]
    lines += [f"{e:<3} {POTCAR_DATASETS[e]:<5} {n:>4} {ZVAL[POTCAR_DATASETS[e]]:>3}" for e, n in counts.items()]
    lines.append(f"# expected NELECT = {nelect(counts)}")
    return "\n".join(lines) + "\n"


def resources(c: Case, s: Settings) -> dict:
    if c.kind == "molecule":
        return {"partition": "single", "ntasks": 16, "time": "72:00:00" if s.mode == "production" else "04:00:00"}
    return {"partition": "workq", "ntasks": 64, "time": "72:00:00" if s.mode == "production" else "12:00:00"}


def estimate(c: Case, s: Settings) -> dict:
    """Order-of-magnitude cost model; uncalibrated, replace with the smoke-run timings."""
    counts = c.counts(); ne = nelect(counts); nions = len(c.symbols)
    res = resources(c, s)
    groups = res["ntasks"] // 4
    nb = max(math.ceil((ne + 2) / 2) + max(nions // 2, 3), int(0.6 * ne))
    nb = int(math.ceil(nb / groups) * groups)
    bohr3 = (c.lengths[0] * c.lengths[1] * c.notes["cell_c"]) / 0.529177 ** 3
    npw = 0.5 * bohr3 * (s.encut / 13.605693) ** 1.5 / (6 * math.pi ** 2)  # Gamma-only half sphere
    # Calibration assumption: NB=1000, NPW=2e5 on 64 cores ~ 30 s per electronic step.
    t_scf = (6.7e-9 * nb * nb * npw + 1.7e-7 * nb * npw * math.log2(max(npw, 2))) / res["ntasks"]
    steps = 300 if s.mode == "production" else 3
    scf = 35 + 15 * (steps - 1) if s.mode == "smoke" else 35 + 12 * min(steps, 150)
    mem = 4 * nb * npw * 16 / 1e9 + 1.0
    return {"label": "ESTIMATE (uncalibrated scaling model; not a measurement)",
            "nelect": ne, "nions": nions, "nbands_vasp_default": nb, "plane_waves_gamma": int(round(npw, -3)),
            "cores": res["ntasks"], "seconds_per_electronic_step": round(t_scf, 1),
            "electronic_steps_assumed": scf, "wall_hours": round(t_scf * scf / 3600.0, 2),
            "memory_gb_total": round(mem, 1)}


def sbatch_text(c: Case, s: Settings) -> str:
    r = resources(c, s)
    ne = nelect(c.counts())
    return f"""#!/bin/bash
# {ADAPTED_FROM}
# ITO/SAM DFT reference point: {c.name} ({s.mode} mode). Submit from this directory: sbatch run.sbatch
#SBATCH -p {r['partition']}
#SBATCH -N 1
#SBATCH -n {r['ntasks']}
#SBATCH -c 1
#SBATCH -t {r['time']}
#SBATCH -A loni_perovsk27
#SBATCH -J ito_{c.name}
#SBATCH -o vasp.cpu.%j.out

set -euo pipefail
cd "${{SLURM_SUBMIT_DIR:-$PWD}}"

# 1. Inputs must match the committed hashes (catches CRLF conversion or local edits).
sha256sum -c inputs.sha256

# 2. Assemble POTCAR from the licensed PAW tree (InterfaceForge search order); never overwrite one.
if [ ! -s POTCAR ]; then
    root=""
    for cand in "${{IFACE_POTCAR_ROOT:-}}" "${{VASP_PP_PATH:+$VASP_PP_PATH/potpaw_PBE}}" "${{VASP_PP_PATH:-}}" "$HOME/pot/potpaw_PBE"; do
        if [ -n "$cand" ] && [ -d "$cand" ]; then root="$cand"; break; fi
    done
    if [ -z "$root" ]; then
        echo "No PBE PAW tree: set IFACE_POTCAR_ROOT or VASP_PP_PATH" >&2; exit 2
    fi
    rm -f POTCAR.tmp
    for ds in $(awk '!/^#/ && NF >= 2 {{print $2}}' POTCAR.spec); do
        src="$root/$ds/POTCAR"
        if [ ! -s "$src" ]; then echo "Missing licensed POTCAR source: $src" >&2; rm -f POTCAR.tmp; exit 2; fi
        cat "$src" >> POTCAR.tmp
    done
    mv POTCAR.tmp POTCAR
    echo "POTCAR assembled from $root"
fi

# 3. Refuse to run if the POTCAR valences disagree with POTCAR.spec.
got=$(grep -E 'POMASS.*ZVAL' POTCAR | sed -E 's/.*ZVAL *= *([0-9.]+).*/\\1/' | awk '{{printf "%d ", $1+0.5}}')
want=$(awk '!/^#/ && NF >= 4 {{printf "%d ", $4}}' POTCAR.spec)
if [ "$got" != "$want" ]; then echo "ZVAL mismatch: POTCAR [$got] vs POTCAR.spec [$want]" >&2; exit 3; fi
nel=$(paste <(awk '!/^#/ && NF >= 4 {{print $3}}' POTCAR.spec) <(echo $got | tr ' ' '\\n') | awk '{{s += $1*$2}} END {{print s}}')
if [ "$nel" != "{ne}" ]; then echo "NELECT $nel differs from expected {ne}" >&2; exit 3; fi
grep -E 'TITEL|ENMAX' POTCAR | tee potcar_titles.txt

# 4. Run (LONI VASP module convention from InterfaceForge runvasp.sh).
module purge
VASP_MODULE="${{IFACE_VASP_MODULE:-vasp6/6.5.1-cpu}}"
module load "$VASP_MODULE"
export OMP_NUM_THREADS=1
export SINGULARITYENV_OMP_NUM_THREADS=1
{{
    echo "case={c.name}"; echo "mode={s.mode}"; echo "job=${{SLURM_JOB_ID:-none}}"; echo "host=$(hostname)"
    echo "module=$VASP_MODULE"; echo "start=$(date -Is)"
    echo "repo_commit=$(git -C "$PWD" rev-parse HEAD 2>/dev/null || echo unknown)"
    sha256sum POSCAR INCAR KPOINTS POTCAR.spec POTCAR
}} > run_provenance.txt
SECONDS=0
srun -n "${{SLURM_NTASKS:-{r['ntasks']}}}" vasp_gam
echo "end=$(date -Is)" >> run_provenance.txt
echo "took $SECONDS sec." | tee -a run_provenance.txt
"""


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _source_hashes() -> dict[str, str]:
    here = Path(__file__).resolve().parent
    out = {}
    for p in (here / "vasp.py", here / "substrate.py"):
        out[str(p.relative_to(ROOT)).replace("\\", "/")] = hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    return out


def _formula(counts: dict[str, int]) -> str:
    """Hill-order formula (C, H first when C is present, otherwise alphabetical)."""
    keys = sorted(counts)
    if "C" in counts:
        keys = ["C"] + (["H"] if "H" in counts else []) + [k for k in keys if k not in ("C", "H")]
    return "".join(f"{k}{counts[k] if counts[k] > 1 else ''}" for k in keys)


def case_files(c: Case, s: Settings) -> dict[str, str]:
    files = {"POSCAR": poscar_text(c), "INCAR": incar_text(c, s), "KPOINTS": KPOINTS_TEXT,
             "POTCAR.spec": potcar_spec_text(c), "geometry.extxyz": extxyz_text(c)}
    files["inputs.sha256"] = "".join(f"{_sha(files[k])}  {k}\n" for k in ("POSCAR", "INCAR", "KPOINTS", "POTCAR.spec", "geometry.extxyz"))
    files["run.sbatch"] = sbatch_text(c, s)
    counts = c.counts()
    notes = {k: v for k, v in c.notes.items() if k not in ("origin", "trilayer_index")}
    manifest = {
        "case": c.name,
        "kind": c.kind,
        "mode": s.mode,
        "purpose": c.purpose,
        "formula": _formula(counts),
        "counts": counts,
        "atom_count": len(c.symbols),
        "nelect_expected": nelect(counts),
        "potcar": [{"element": e, "dataset": POTCAR_DATASETS[e], "zval_expected": ZVAL[POTCAR_DATASETS[e]]} for e in counts],
        "cell_angstrom": [round(c.lengths[0], 6), round(c.lengths[1], 6), c.notes["cell_c"]],
        "vacuum_gap_angstrom": c.notes["vacuum_gap"],
        "constraints": {
            "selective_dynamics": c.kind != "molecule",
            "fixed_atoms": int(c.fixed.sum()),
            "free_atoms": int((~c.fixed).sum()),
            "rule": "bottom O-In-O trilayer fixed" + ("; Sn cluster shell freed" if "sn_cluster" in c.notes else "") if c.kind != "molecule" else "none (isolated molecule)",
        },
        "distances_angstrom": c.notes["distances"],
        "estimate": estimate(c, s),
        "resources": resources(c, s),
        "settings": {"trilayers": s.trilayers, "hydroxyl_pairs_per_primitive_cell": s.hydroxyl_pairs_per_primitive,
                     "hydroxyl_seed": s.hydroxyl_seed, "phys_height": s.phys_height, "bidentate_in_o": s.bidentate_in_o,
                     "vacuum_slab": s.vacuum_slab, "vacuum_molecule": s.vacuum_molecule, "common_cell": s.common_cell,
                     "encut": s.encut},
        "construction": notes,
        "lattice": BIXBYITE,
        "provenance": {"generator": "nio_md_prep.ito.vasp", "source_sha256": _source_hashes(),
                       "interfaceforge_conventions": INTERFACEFORGE_SHA},
        "sha256": {k: _sha(v) for k, v in files.items()},
        "file_order_note": "POSCAR and geometry.extxyz share one atom order (species-sorted, stable)",
    }
    files["case_manifest.json"] = json.dumps(manifest, indent=2, default=_json_default) + "\n"
    return files


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def readme_text(cases: list[Case], s: Settings) -> str:
    rows = ["| case | kind | formula | atoms | fixed | NELECT | min nonbonded (A) | min mol-slab (A) | vacuum (A) | partition x cores | est. wall h |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|"]
    for c in cases:
        e = estimate(c, s); r = resources(c, s); d = c.notes["distances"]
        rows.append(f"| `{c.name}` | {c.kind} | {_formula(c.counts())} | {len(c.symbols)} | {int(c.fixed.sum())} | "
                    f"{nelect(c.counts())} | {d['min_nonbonded']:.3f} | {_fmt(d.get('min_molecule_slab_nonbonded'))} | "
                    f"{c.notes['vacuum_gap']:.2f} | {r['partition']} x {r['ntasks']} | {e['wall_hours']:.2f} |")
    purposes = [f"- `{c.name}`: {c.purpose}" for c in cases]
    return "\n".join([
        f"# ITO/SAM VASP reference inputs ({s.mode} mode, {s.trilayers} trilayers)",
        "",
        "Generated by `python scripts/ito/prepare_vasp_smokes.py`; do not edit by hand. Design, formulas and",
        "run instructions: [`docs/ito/vasp-smoke.md`](../../../docs/ito/vasp-smoke.md).",
        "",
        f"Common cell {cases[0].lengths[0]:.4f} x {cases[0].lengths[1]:.4f} x {cases[0].notes['cell_c']:.1f} A, Gamma only, "
        f"PBE+D3(BJ), ENCUT {s.encut:.0f} eV, NSW {300 if s.mode == 'production' else 3}.",
        "Wall-time column is an ESTIMATE from an uncalibrated scaling model, not a measurement.",
        "",
        *rows,
        "",
        *purposes,
        "",
        "Each directory: `POSCAR` (selective dynamics for slabs), `INCAR`, `KPOINTS`, `POTCAR.spec`, `geometry.extxyz`",
        "(same atom order as POSCAR, with substrate labels and LigParGen ids for classical cross-evaluation),",
        "`inputs.sha256`, `run.sbatch`, `case_manifest.json`. `submit_all.sh` submits every case (not run by the generator).",
        "",
    ])


def _fmt(v) -> str:
    return "-" if v is None else f"{v:.3f}"


def submit_all_text(cases: list[Case]) -> str:
    body = "\n".join(f'(cd "$here/{c.name}" && sbatch run.sbatch)' for c in cases)
    return ("#!/bin/bash\n# Submit every ITO VASP reference case. Review docs/ito/vasp-smoke.md first.\n"
            "set -euo pipefail\nhere=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n" + body + "\n")


def render(settings: Settings | None = None) -> dict[str, str]:
    """All output files, keyed by relative path; pure function of the inputs."""
    settings = settings or Settings()
    cases = build_cases(settings)
    out: dict[str, str] = {}
    for c in cases:
        for name, text in case_files(c, settings).items():
            out[f"{c.name}/{name}"] = text
    out["README.md"] = readme_text(cases, settings)
    out["submit_all.sh"] = submit_all_text(cases)
    return out


def write_all(output: Path, settings: Settings | None = None) -> dict[str, str]:
    files = render(settings)
    for rel, text in files.items():
        path = output / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    return files
