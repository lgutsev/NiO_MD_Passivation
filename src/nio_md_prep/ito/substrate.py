"""Deterministic In2O3(111) / Sn-doped ITO slab construction.

The slab is cut from the bixbyite structure (space group Ia-3, No. 206) in an
orthorhombic (111) setting with in-plane vectors a[1-10] and a[11-2] and
normal a[111].  Along [111] the crystal is a stack of neutral,
composition-symmetric O-In-O trilayers (In32O48 per orthorhombic cell,
d(222) = a/(2*sqrt(3))).  Cutting between trilayers therefore gives a
stoichiometric, dipole-free Tasker type-II slab with identical terminations,
which is checked rather than assumed.

Optional modifications, applied only to the top (adsorption) surface or the
slab interior:

* dissociative hydroxylation: OH- on a fivefold In (In5c) plus H+ on a
  nearby threefold surface O (O3c), placed geometrically (not DFT-relaxed);
* Sn doping as neutral Frank-Koestlin-type (2 Sn_In + O_i) clusters, O_i on
  the empty 16c fluorite anion sites of bixbyite.

Fixed-charge models cannot carry free electrons, so Sn donors are always
ionically compensated here and oxygen vacancies are not generated.
"""
from __future__ import annotations

import hashlib
import json
import math
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Marezio, Acta Cryst. 20, 723 (1966), as deposited in COD 2310009 (checked
# against the COD CIF on 2026-09-26); see docs/ito/parameter-sources.md.
BIXBYITE = {
    "a": 10.117,
    "in_8b": (0.25, 0.25, 0.25),
    "in_24d_x": 0.4663,
    "o_48e": (0.3912, 0.1558, 0.3796),
    "source": "Marezio, Acta Cryst. 20, 723 (1966), doi:10.1107/S0365110X66001749; coordinates from COD 2310009",
}
ORIENTATION = ((1, -1, 0), (1, 1, -2), (1, 1, 1))
MASSES = {"In": 114.818, "Sn": 118.710, "O": 15.999, "Oh": 15.999, "Hh": 1.008}
FORMAL = {"In": 3.0, "Sn": 4.0, "O": -2.0}
TYPE_ORDER = ("In", "Sn", "O", "Oh", "Hh")
IN_O_BOND_CUTOFF = 2.6
IN_OH_LENGTH = 2.10
O_H_LENGTH = 0.97


@dataclass
class Slab:
    labels: list[str]
    positions: np.ndarray
    lengths: tuple[float, float]
    notes: dict = field(default_factory=dict)

    def count(self, label: str) -> int:
        return sum(1 for x in self.labels if x == label)


def _bulk_with_vacancies():
    """Conventional bixbyite cell plus its 16 empty fluorite anion sites ("X")."""
    from ase import Atoms
    from ase.spacegroup import crystal

    a = BIXBYITE["a"]
    bulk = crystal(
        ["In", "In", "O"],
        basis=[BIXBYITE["in_8b"], (BIXBYITE["in_24d_x"], 0.0, 0.25), BIXBYITE["o_48e"]],
        spacegroup=206,
        cellpar=[a, a, a, 90, 90, 90],
    )
    if bulk.get_chemical_formula() != "In32O48":
        raise ValueError(f"bixbyite cell has formula {bulk.get_chemical_formula()}")
    ideal = np.array(
        [((i + 0.5) / 4, (j + 0.5) / 4, (k + 0.5) / 4) for i in range(4) for j in range(4) for k in range(4)]
    ) * a
    oxygen = bulk.positions[np.array(bulk.get_chemical_symbols()) == "O"]
    vacant = []
    for site in ideal:
        d = oxygen - site
        d -= a * np.round(d / a)
        if np.min(np.linalg.norm(d, axis=1)) > 1.0:
            vacant.append(site)
    if len(vacant) != 16:
        raise ValueError(f"expected 16 structural anion vacancies, found {len(vacant)}")
    full = bulk + Atoms("X16", positions=vacant, cell=bulk.cell, pbc=True)
    return full


def oriented_cell() -> tuple[list[str], np.ndarray, np.ndarray]:
    """Symbols, Cartesian positions, and lengths of the orthorhombic (111) cell."""
    from ase.build import make_supercell

    cell = make_supercell(_bulk_with_vacancies(), np.array(ORIENTATION))
    lengths = cell.cell.lengths()
    frac = cell.get_scaled_positions(wrap=True)
    return cell.get_chemical_symbols(), frac * lengths, lengths


def trilayer_spacing() -> float:
    return BIXBYITE["a"] / (2.0 * math.sqrt(3.0))


def build_slab(nx: int, ny: int, trilayers: int) -> tuple[Slab, np.ndarray]:
    """Stoichiometric (111) slab; returns the slab and interior O_i candidate sites."""
    if nx < 1 or ny < 1 or trilayers < 2:
        raise ValueError("slab needs nx, ny >= 1 and at least two trilayers")
    symbols, pos, lengths = oriented_cell()
    symbols = np.array(symbols)
    d = trilayer_spacing()
    in_z = pos[symbols == "In", 2]
    # In planes sit at trilayer centres; cut half a spacing away from one.
    cut = (np.min(in_z) + 0.5 * d) % d
    periods = math.ceil(trilayers * d / lengths[2]) + 1
    labels, coords = [], []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(periods):
                shift = np.array([ix * lengths[0], iy * lengths[1], iz * lengths[2] - cut])
                for s, p in zip(symbols, pos + shift):
                    if 0.0 <= p[2] < trilayers * d:
                        labels.append(str(s)); coords.append(p)
    coords = np.array(coords)
    lx, ly = nx * lengths[0], ny * lengths[1]
    is_x = np.array([s == "X" for s in labels])
    sites = coords[is_x]
    labels = [s for s in labels if s != "X"]
    coords = coords[~is_x]
    n_in, n_o = labels.count("In"), labels.count("O")
    if 3 * n_in != 2 * n_o:
        raise ValueError(f"slab is not stoichiometric: In{n_in} O{n_o}")
    slab = Slab(labels, coords, (lx, ly))
    slab.notes.update(
        {
            "orthorhombic_cell_angstrom": [round(float(x), 6) for x in lengths],
            "trilayer_spacing_angstrom": round(d, 6),
            "cut_offset_angstrom": round(float(cut), 6),
        }
    )
    return slab, sites


def _periodic_delta(delta: np.ndarray, lengths: tuple[float, float]) -> np.ndarray:
    out = delta.copy()
    for axis in (0, 1):
        out[..., axis] -= lengths[axis] * np.round(out[..., axis] / lengths[axis])
    return out


def coordination(slab: Slab) -> tuple[np.ndarray, list[list[int]]]:
    """Cation-anion coordination numbers with periodic x/y (cutoff 2.6 A)."""
    pos = slab.positions
    cation = np.array([x in ("In", "Sn") for x in slab.labels])
    anion = np.array([x in ("O", "Oh") for x in slab.labels])
    neighbors: list[list[int]] = [[] for _ in slab.labels]
    ci, ai = np.where(cation)[0], np.where(anion)[0]
    for i in ci:
        d = np.linalg.norm(_periodic_delta(pos[ai] - pos[i], slab.lengths), axis=1)
        for j in ai[d < IN_O_BOND_CUTOFF]:
            neighbors[i].append(int(j)); neighbors[int(j)].append(int(i))
    return np.array([len(n) for n in neighbors]), neighbors


def _farthest_point(points: np.ndarray, count: int, lengths, rng) -> list[int]:
    chosen = [int(rng.integers(len(points)))]
    dist = np.linalg.norm(_periodic_delta(points[:, :2] - points[chosen[0], :2], lengths), axis=1)
    while len(chosen) < count:
        nxt = int(np.argmax(dist))
        chosen.append(nxt)
        dist = np.minimum(dist, np.linalg.norm(_periodic_delta(points[:, :2] - points[nxt, :2], lengths), axis=1))
    return chosen


MIN_OO_TERMINAL = 2.5   # terminal-OH O to any lattice O
MIN_H_CONTACT = 1.5     # new H to any atom other than its own O
MIN_H_CATION = 2.5      # new H to any In/Sn (rejects acute M-O-H angles)


def hydroxylate(slab: Slab, pairs: int, seed: int, exposed_z_min: float | None = None,
                cation_cn: tuple[int, ...] = (5,), anion_cn: tuple[int, ...] = (3,)) -> dict:
    """Dissociate ``pairs`` waters on the top surface: In5c-OH + O3c-H.

    The terminal O goes along the missing-octahedron direction of an In5c.
    For one of the three top-face In5c classes of the bulk-terminated (111)
    face that position lies ~2.1-2.2 A from a lattice O, so such In5c are
    excluded before the seeded farthest-point selection.  The proton goes on
    the nearest O3c that is >= 2.5 A from the terminal O and whose H is
    >= 1.5 A from every other atom (new atoms included) and >= 2.5 A from every
    cation (no acute In-O-H angles).
    """
    if pairs <= 0:
        return {"pairs": 0}
    rng = np.random.default_rng(seed)
    cn, nbr = coordination(slab)
    pos = slab.positions
    L = slab.lengths
    # Flat slabs: the upper half.  Corrugated slabs pass a lower bound that keeps
    # groove floors and step walls while excluding the bottom face.
    z_mid = 0.5 * (pos[:, 2].min() + pos[:, 2].max()) if exposed_z_min is None else exposed_z_min
    top_in = [i for i, s in enumerate(slab.labels) if s == "In" and cn[i] in cation_cn and pos[i, 2] > z_mid]
    top_o = [i for i, s in enumerate(slab.labels) if s == "O" and cn[i] in anion_cn and pos[i, 2] > z_mid]
    lattice_o = [i for i, s in enumerate(slab.labels) if s == "O"]
    z_hat = np.array([0.0, 0.0, 1.0])
    terminal, excluded = {}, []
    cations = pos[[i for i, s in enumerate(slab.labels) if s in ("In", "Sn")]]
    for i in top_in:
        bonds = _periodic_delta(pos[nbr[i]] - pos[i], L)
        vac = -np.sum(bonds / np.linalg.norm(bonds, axis=1)[:, None], axis=0); vac /= np.linalg.norm(vac)
        o_t = pos[i] + IN_OH_LENGTH * vac
        h_t = o_t + O_H_LENGTH * (vac + z_hat) / np.linalg.norm(vac + z_hat)
        doo = float(np.min(np.linalg.norm(_periodic_delta(pos[lattice_o] - o_t, L), axis=1)))
        dhm = float(np.min(np.linalg.norm(_periodic_delta(cations - h_t, L), axis=1)))
        if doo < MIN_OO_TERMINAL or dhm < MIN_H_CATION:
            excluded.append({"in_index": int(i), "terminal_o_to_lattice_o": round(doo, 4), "terminal_h_to_cation": round(dhm, 4)})
        else:
            terminal[i] = (o_t, vac)
    allowed = sorted(terminal)
    if pairs > min(len(allowed), len(top_o)):
        raise ValueError(f"requested {pairs} hydroxyl pairs; only {len(allowed)} clash-free In5c / {len(top_o)} O3c on top")
    chosen = [allowed[k] for k in _farthest_point(pos[allowed], pairs, L, rng)]
    new_labels, new_pos, records, used = [], [], [], set()
    for i in chosen:
        o_t, vac = terminal[i]
        h_dir = vac + z_hat; h_dir /= np.linalg.norm(h_dir)
        h_t = o_t + O_H_LENGTH * h_dir
        cands = [j for j in top_o if j not in used]
        dist = np.linalg.norm(_periodic_delta(pos[cands] - o_t, L), axis=1)
        placed = None
        for k in np.argsort(dist, kind="stable"):
            j = cands[int(k)]
            if dist[k] < MIN_OO_TERMINAL:
                continue
            ob = _periodic_delta(pos[nbr[j]] - pos[j], L)
            out = -np.sum(ob / np.linalg.norm(ob, axis=1)[:, None], axis=0)
            out = out / np.linalg.norm(out) if np.linalg.norm(out) > 1e-6 else z_hat
            if out[2] < 0.2:
                out = out + z_hat; out /= np.linalg.norm(out)
            h_j = pos[j] + O_H_LENGTH * out
            others = np.vstack([np.delete(pos, j, axis=0), np.array(new_pos + [o_t, h_t]).reshape(-1, 3)])
            if (float(np.min(np.linalg.norm(_periodic_delta(others - h_j, L), axis=1))) >= MIN_H_CONTACT
                    and float(np.min(np.linalg.norm(_periodic_delta(cations - h_j, L), axis=1))) >= MIN_H_CATION):
                placed = (j, h_j, float(dist[k]))
                break
        if placed is None:
            raise ValueError(f"no acceptable O3c for the proton of In {i}")
        j, h_j, d_oo = placed
        used.add(j)
        slab.labels[j] = "Oh"
        new_labels += ["Oh", "Hh", "Hh"]; new_pos += [o_t, h_t, h_j]
        records.append({"in_index": int(i), "protonated_o_index": int(j), "vacancy_direction_z": round(float(vac[2]), 4),
                        "terminal_o_to_protonated_o": round(d_oo, 4)})
    n0 = len(pos)
    slab.labels += new_labels
    slab.positions = np.vstack([pos, np.array(new_pos)])
    # Closest non-bonded contact of any new atom against all atoms, new ones included.
    bonded = {(n0 + 3 * k, n0 + 3 * k + 1) for k in range(len(chosen))} | {(r["protonated_o_index"], n0 + 3 * k + 2) for k, r in enumerate(records)}         | {(r["in_index"], n0 + 3 * k) for k, r in enumerate(records)}
    min_sep = np.inf
    for a in range(n0, len(slab.labels)):
        d = np.linalg.norm(_periodic_delta(slab.positions - slab.positions[a], L), axis=1)
        d[a] = np.inf
        for b in range(len(d)):
            if (a, b) in bonded or (b, a) in bonded:
                d[b] = np.inf
        min_sep = min(min_sep, float(d.min()))
    return {"pairs": pairs, "top_in5c_available": len(top_in), "top_in5c_excluded_clash": excluded,
            "top_o3c_available": len(top_o), "seed": seed,
            "site_selection": "seeded periodic farthest-point over clash-free top In5c; O3c >= 2.5 A from terminal O; every H >= 1.5 A from all atoms and >= 2.5 A from cations",
            "min_new_atom_nonbonded_angstrom": round(min_sep, 4), "sites": records}


def groove_surface_height(coord: np.ndarray, period: float, center: float, top: float, depth: float, slope: float) -> np.ndarray:
    """Height of an ideal symmetric V-groove surface (periodic along the profile axis)."""
    d = coord - center
    d -= period * np.round(d / period)
    return top - np.clip(depth - slope * np.abs(d), 0.0, depth)


def carve_groove(slab: Slab, depth_trilayers: int, slope: float = 1.0, axis: int = 0,
                 center_fraction: float = 0.5, sites: np.ndarray | None = None) -> tuple[dict, np.ndarray | None]:
    """Cut a V-groove running perpendicular to ``axis`` out of a flat (111) slab.

    The walls are staircases of whole O-In-O trilayers: at profile coordinate u
    the number of removed top trilayers is round((depth - slope*|u-u0|)/d), so
    each step is one neutral trilayer (2.92 A) high, mirroring the monatomic
    2.085 A steps of the corrugated NiO(110) slab.  Lateral cuts through a
    trilayer at the step edges can leave a non-stoichiometric rim, so dangling
    atoms (O with CN<=1, cations with CN<=2) are removed and neutrality is then
    restored by removing the lowest-coordinated exposed rim atoms; every removal
    is recorded.
    """
    d = trilayer_spacing()
    period = slab.lengths[axis]
    center = center_fraction * period
    z = slab.positions[:, 2]
    n_tri = int(round((z.max() - z.min()) / d + 0.5))
    if depth_trilayers >= n_tri - 2:
        raise ValueError("groove must leave at least two full trilayers under its floor")
    depth = depth_trilayers * d
    u = slab.positions[:, axis]
    du = u - center; du -= period * np.round(du / period)
    removed_layers = np.rint(np.clip(depth - slope * np.abs(du), 0.0, depth) / d).astype(int)
    layer = np.floor(z / d + 1e-9).astype(int)
    keep = layer < (n_tri - removed_layers)
    carved = int((~keep).sum())
    slab.labels = [s for s, k in zip(slab.labels, keep) if k]
    slab.positions = slab.positions[keep]
    floor_z = (n_tri - depth_trilayers) * d
    exposed_z_min = 2.0 * d
    # Dangling atoms, then charge-neutralizing rim removals.
    removals = []
    formal = {"In": 3, "Sn": 4, "O": -2}
    for _ in range(10000):
        cn, nbr = coordination(slab)
        lab = slab.labels; zz = slab.positions[:, 2]
        dangling = [i for i, s in enumerate(lab) if zz[i] > exposed_z_min and ((s == "O" and cn[i] <= 1) or (s in ("In", "Sn") and cn[i] <= 2))]
        q = sum(formal[s] for s in lab)
        if dangling:
            pick = dangling[0]; why = "dangling"
        elif q != 0:
            want = "O" if q < 0 else "In"
            cands = [i for i, s in enumerate(lab) if s == want and zz[i] > exposed_z_min]
            pick = min(cands, key=lambda i: (cn[i], -zz[i], i)); why = "neutralize"
        else:
            break
        removals.append({"label": lab[pick], "cn": int(cn[pick]), "reason": why,
                         "position": [round(float(x), 4) for x in slab.positions[pick]]})
        del slab.labels[pick]
        slab.positions = np.delete(slab.positions, pick, axis=0)
    else:
        raise RuntimeError("rim neutralization did not converge")
    kept_sites = None
    if sites is not None and len(sites):
        # Interstitial sites must stay at least one trilayer below the local carved surface.
        su = sites[:, axis] - center; su -= period * np.round(su / period)
        s_removed = np.rint(np.clip(depth - slope * np.abs(su), 0.0, depth) / d).astype(int)
        s_top = (n_tri - s_removed) * d
        kept_sites = sites[sites[:, 2] < s_top - d]
    record = {
        "profile_axis_before_rotation": "xy"[axis], "period_angstrom": round(period, 6),
        "center_fraction": center_fraction, "depth_trilayers": depth_trilayers,
        "depth_angstrom": round(depth, 6), "wall_slope": slope,
        "step_height_angstrom": round(d, 6), "step_run_angstrom": round(d / slope, 6),
        "slab_trilayers": n_tri, "trilayers_under_floor": n_tri - depth_trilayers,
        "floor_top_z_nominal": round(floor_z, 6), "plateau_top_z_nominal": round(n_tri * d, 6),
        "opening_width_nominal_angstrom": round(2 * depth / slope, 6),
        "atoms_carved": carved, "rim_removals": removals, "exposed_z_min": round(exposed_z_min, 6),
    }
    return record, kept_sites


def rotate_quarter(slab: Slab, sites: np.ndarray | None = None) -> np.ndarray | None:
    """Rotate by +90 deg about z: (x, y) -> (y, Lx - x); keeps handedness."""
    lx, ly = slab.lengths
    p = slab.positions.copy()
    slab.positions = np.column_stack([p[:, 1] % ly, (lx - p[:, 0]) % lx, p[:, 2]])
    slab.lengths = (ly, lx)
    if sites is None:
        return None
    return np.column_stack([sites[:, 1] % ly, (lx - sites[:, 0]) % lx, sites[:, 2]])


def dope_sn(slab: Slab, sites: np.ndarray, target_fraction: float, seed: int, min_separation: float = 5.0) -> dict:
    """Neutral (2 Sn_In + O_i) clusters, O_i on empty 16c sites away from the top surface."""
    if target_fraction <= 0:
        return {"clusters": 0, "sn_fraction": 0.0}
    rng = np.random.default_rng(seed)
    n_in = slab.count("In")
    clusters = int(round(target_fraction * n_in / 2.0))
    z = slab.positions[:, 2]
    top, bottom = z.max(), z.min()
    # Keep interstitials at least one trilayer inside both faces.
    inner = sites[(sites[:, 2] > bottom + trilayer_spacing()) & (sites[:, 2] < top - trilayer_spacing())]
    order = rng.permutation(len(inner))
    chosen: list[np.ndarray] = []
    for k in order:
        p = inner[k]
        if all(np.linalg.norm(_periodic_delta(p - q, slab.lengths)) >= min_separation for q in chosen):
            chosen.append(p)
        if len(chosen) == clusters:
            break
    if len(chosen) < clusters:
        raise ValueError(f"could place only {len(chosen)} of {clusters} O_i sites at {min_separation} A separation")
    in_idx = np.array([i for i, s in enumerate(slab.labels) if s == "In"])
    substituted = []
    for p in chosen:
        free = np.array([i for i in in_idx if slab.labels[i] == "In"])
        d = np.linalg.norm(_periodic_delta(slab.positions[free] - p, slab.lengths), axis=1)
        for i in free[np.argsort(d)[:2]]:
            slab.labels[int(i)] = "Sn"; substituted.append(float(np.linalg.norm(_periodic_delta(slab.positions[int(i)] - p, slab.lengths))))
    slab.labels += ["O"] * len(chosen)
    slab.positions = np.vstack([slab.positions, np.array(chosen)])
    n_sn = slab.count("Sn")
    return {
        "clusters": len(chosen),
        "compensation": "Frank-Koestlin-type 2 Sn_In + O_i on empty 16c anion sites (unrelaxed)",
        "sn_fraction": round(n_sn / (n_sn + slab.count("In")), 6),
        "target_sn_fraction": target_fraction,
        "seed": seed,
        "sn_to_oi_distance_angstrom": [round(min(substituted), 4), round(max(substituted), 4)],
    }


def charges(slab: Slab, scale: float, hydroxyl_h: float) -> np.ndarray:
    """Formal charges times ``scale``; OH groups carry -scale in total.

    With scale=0.525 and hydroxyl_h=0.425 this reproduces the CLAYFF charge
    pattern (bulk O -1.05, hydroxyl O -0.95, H +0.425), so dissociative
    hydroxylation and (2 Sn + O_i) doping both remain exactly neutral.
    """
    out = np.empty(len(slab.labels))
    for i, s in enumerate(slab.labels):
        if s in FORMAL: out[i] = scale * FORMAL[s]
        elif s == "Oh": out[i] = -scale - hydroxyl_h
        elif s == "Hh": out[i] = hydroxyl_h
        else: raise ValueError(f"unknown slab label {s}")
    return out


def _dipole_per_area(slab: Slab, q: np.ndarray) -> float:
    z = slab.positions[:, 2]
    return float(np.sum(q * (z - z.mean())) / (slab.lengths[0] * slab.lengths[1]))


def _git_commit(root: Path) -> str | None:
    import subprocess
    try:
        run = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=True)
        return run.stdout.strip()
    except Exception:
        return None


def write_lammps(slab: Slab, q: np.ndarray, path: Path, z_bounds: tuple[float, float], title: str) -> dict[str, int]:
    present = [t for t in TYPE_ORDER if t in slab.labels]
    type_id = {t: k + 1 for k, t in enumerate(present)}
    lines = [title, "", f"{len(slab.labels)} atoms", "", f"{len(present)} atom types", ""]
    lines += [f"0.00000000 {slab.lengths[0]:.8f} xlo xhi", f"0.00000000 {slab.lengths[1]:.8f} ylo yhi",
              f"{z_bounds[0]:.8f} {z_bounds[1]:.8f} zlo zhi", "", "Masses", ""]
    lines += [f"{type_id[t]} {MASSES[t]} # {t}" for t in present]
    lines += ["", "Atoms # full", ""]
    for i, (s, p, c) in enumerate(zip(slab.labels, slab.positions, q), 1):
        x = p[0] % slab.lengths[0]; y = p[1] % slab.lengths[1]
        lines.append(f"{i} 0 {type_id[s]} {c:.6f} {x:.6f} {y:.6f} {p[2]:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return type_id


def build_from_model(model_path: Path, output: Path) -> dict:
    """Build a slab from ``model.toml`` and write surface.lmp/.xyz + manifest."""
    with model_path.open("rb") as f:
        model = tomllib.load(f)
    s = model["slab"]; h = model.get("hydroxylation", {}); dop = model.get("doping", {}); ch = model["charges"]
    slab, sites = build_slab(int(s["nx"]), int(s["ny"]), int(s["trilayers"]))
    # Primitive (111) surface cell = half the orthorhombic cell.
    primitive_cells = 2 * int(s["nx"]) * int(s["ny"])
    cn0, _ = coordination(slab)
    z0 = slab.positions[:, 2]
    bare_in5c = int(sum(1 for i, x in enumerate(slab.labels) if x == "In" and cn0[i] == 5 and z0[i] > z0.mean()))
    bare_o3c = int(sum(1 for i, x in enumerate(slab.labels) if x == "O" and cn0[i] == 3 and z0[i] > z0.mean()))
    corr = model.get("corrugation")
    corrugation = None
    hyd_kwargs: dict = {}
    if corr:
        corrugation, sites = carve_groove(slab, int(corr["depth_trilayers"]), float(corr.get("wall_slope", 1.0)),
                                          {"x": 0, "y": 1}[corr.get("profile_axis", "x")],
                                          float(corr.get("center_fraction", 0.5)), sites)
        # Groove floors and walls count as exposed; the bottom face does not.  Step-edge
        # cations can be fourfold, step-edge anions twofold.
        hyd_kwargs = {"exposed_z_min": corrugation["exposed_z_min"], "cation_cn": (4, 5), "anion_cn": (2, 3)}
    doping = dope_sn(slab, sites, float(dop.get("sn_fraction", 0.0)), int(dop.get("seed", 1)))
    pairs = int(round(float(h.get("pairs_per_primitive_cell", 0.0)) * primitive_cells))
    hydroxyl = hydroxylate(slab, pairs, int(h.get("seed", 1)), **hyd_kwargs)
    if corr and corr.get("rotate_groove_along_x", True):
        # Match the NiO slab: groove along x, profile along y.
        rotate_quarter(slab)
        corrugation["profile_axis"] = "y"
        corrugation["groove_center_y_angstrom"] = round((1.0 - corrugation["center_fraction"]) * slab.lengths[1], 6)
    q = charges(slab, float(ch["scale"]), float(ch.get("hydroxyl_h", 0.425)))
    total = float(q.sum())
    if abs(total) > 1e-6:
        raise ValueError(f"slab charge {total:.3e} is not neutral")
    cn, _ = coordination(slab)
    z = slab.positions[:, 2]
    top_cation = float(max(z[i] for i, x in enumerate(slab.labels) if x in ("In", "Sn")))
    output.mkdir(parents=True, exist_ok=True)
    zlo, zhi = float(s.get("zlo", -5.0)), float(s.get("zhi", 200.0))
    title = f"ITO extension slab {model['model']['id']} generated by nio_md_prep.ito.substrate"
    type_id = write_lammps(slab, q, output / "surface.lmp", (zlo, zhi), title)
    xyz = [f"{len(slab.labels)}", f"{model['model']['id']} Lattice=\"{slab.lengths[0]} 0 0 0 {slab.lengths[1]} 0 0 0 {zhi-zlo}\""]
    xyz += [f"{('O' if l == 'Oh' else 'H' if l == 'Hh' else l)} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}" for l, p in zip(slab.labels, slab.positions)]
    (output / "surface.xyz").write_text("\n".join(xyz) + "\n", encoding="utf-8")
    root = Path(__file__).resolve().parents[3]
    area_nm2 = slab.lengths[0] * slab.lengths[1] / 100.0
    top_in5c = sum(1 for i, x in enumerate(slab.labels) if x in ("In", "Sn") and cn[i] == 5 and z[i] > z.mean())
    manifest = {
        "model_id": model["model"]["id"],
        "description": model["model"].get("description", ""),
        "model_toml_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        "surface_lmp_sha256": hashlib.sha256((output / "surface.lmp").read_bytes()).hexdigest(),
        "code_commit": _git_commit(root),
        "lattice": BIXBYITE,
        "orientation_vectors_cubic": ORIENTATION,
        "box_angstrom": {"x": [0.0, slab.lengths[0]], "y": [0.0, slab.lengths[1]], "z": [zlo, zhi]},
        "area_nm2": round(area_nm2, 6),
        "trilayers": int(s["trilayers"]),
        "counts": {t: slab.count(t) for t in TYPE_ORDER},
        "atom_count": len(slab.labels),
        "type_ids": type_id,
        "charges": {"scheme": ch.get("scheme", "scaled-formal"), "scale": float(ch["scale"]), "hydroxyl_h": float(ch.get("hydroxyl_h", 0.425)),
                    "per_label": {t: round(float(q[slab.labels.index(t)]), 6) for t in TYPE_ORDER if t in slab.labels}},
        "total_charge": total,
        "dipole_z_e_per_angstrom": _dipole_per_area(slab, q),
        "z_top_cation_angstrom": round(top_cation, 6),
        "z_top_atom_angstrom": round(float(z.max()), 6),
        "z_bottom_atom_angstrom": round(float(z.min()), 6),
        "top_surface": {"bare_in5c": bare_in5c, "bare_o3c": bare_o3c, "bare_in5c_per_nm2": round(bare_in5c / area_nm2, 4),
                        "cation_5c_after_modification": top_in5c},
        "hydroxylation": hydroxyl | {"oh_groups_per_nm2": round(2 * hydroxyl["pairs"] / area_nm2, 4)},
        "doping": doping,
        "corrugation": corrugation,
        "limitations": [
            "geometric (unrelaxed) hydroxyl and interstitial placement; DFT relaxation not performed",
            "fixed point charges; no free carriers, polarization, or charge transfer",
            "no oxygen vacancies (would require electronic compensation)",
        ],
    }
    (output / "surface_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
