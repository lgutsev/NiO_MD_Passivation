"""Ball-and-stick pictures of rigid slab models (cross-sections and perspective views).

    python scripts/ito/render_surfaces.py OUT.png DIR [DIR ...] [--labels A B ...] [--ncols 5]

Each DIR holds surface.lmp + surface_manifest.json.  For every slab a thin slice
along the groove axis (x) is drawn as a y-z cross-section with bonds
(cation-O < 2.6 A, O-H < 1.15 A), atoms depth-sorted and shaded by x.  With
--perspective a second row shows a 3D ball-and-stick view of a groove segment.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from nio_md_prep.lammps import parse  # noqa: E402

COLORS = {"In": "#8e6fd1", "Sn": "#2a9d8f", "Ni": "#4c9a62", "O": "#d9453b", "Oh": "#f39c34", "Hh": "#f4f4f4"}
RADII = {"In": 0.52, "Sn": 0.50, "Ni": 0.45, "O": 0.36, "Oh": 0.42, "Hh": 0.24}
NAMES = {"In": "In", "Sn": "Sn", "Ni": "Ni", "O": "O (lattice)", "Oh": "O (hydroxyl)", "Hh": "H"}
CATIONS = {"In", "Sn", "Ni"}
MOL_COLORS = {"C": "#5a5a5a", "N": "#3b62c7", "O": "#d9453b", "P": "#f28c28", "H": "#ffffff", "F": "#6cc24a", "I": "#7b3f9e"}
MOL_RADII = {"C": 0.34, "N": 0.33, "O": 0.33, "P": 0.42, "H": 0.20, "F": 0.30, "I": 0.50}


def draw_molecule(ax, path: Path, lengths, slice_center_x: float | None = None):
    """Overlay a molecule (.xyz) on a y-z cross-section: bonds by distance, atoms depth-sorted."""
    rows = [l.split() for l in path.read_text(encoding="utf-8").splitlines()[2:] if l.strip()]
    el = [r[0] for r in rows]; xyz = np.array([[float(v) for v in r[1:4]] for r in rows])
    ref = xyz[[i for i, e in enumerate(el) if e == "P"][0]] if "P" in el else xyz.mean(axis=0)
    xyz[:, 1] -= lengths[1] * np.round((xyz[:, 1] - ref[1]) / lengths[1])
    for i in range(len(el)):
        for j in range(i + 1, len(el)):
            r = np.linalg.norm(xyz[i] - xyz[j]); cut = 1.25 if "H" in (el[i], el[j]) else 1.9
            if r < cut:
                ax.plot([xyz[i, 1], xyz[j, 1]], [xyz[i, 2], xyz[j, 2]], color="#222222", lw=1.6, zorder=20)
    import matplotlib.patches as mp
    for i in np.argsort(-xyz[:, 0]):
        ax.add_patch(mp.Circle((xyz[i, 1], xyz[i, 2]), MOL_RADII.get(el[i], 0.3), facecolor=MOL_COLORS.get(el[i], "#999999"),
                               edgecolor="#111111", lw=0.5, zorder=21))
    return xyz


def load(folder: Path):
    man = json.loads((folder / "surface_manifest.json").read_text(encoding="utf-8"))
    inv = {v: k for k, v in man["type_ids"].items()}
    data = parse(folder / "surface.lmp")
    lab = np.array([inv[int(a.fields[2])] for a in data.sections["Atoms"]])
    xyz = np.array([[float(a.fields[k]) for k in (4, 5, 6)] for a in data.sections["Atoms"]])
    return man, lab, xyz, (data.bounds["x"][1], data.bounds["y"][1])


def bonds(lab, xyz, lengths):
    from scipy.spatial import cKDTree
    p = xyz.copy(); p[:, 0] %= lengths[0]; p[:, 1] %= lengths[1]; p[:, 2] += 1e4
    tree = cKDTree(p, boxsize=[lengths[0], lengths[1], 1e6])
    out = []
    for i, j in tree.query_pairs(2.6):
        a, b = lab[i], lab[j]
        d = xyz[j] - xyz[i]
        d[0] -= lengths[0] * round(d[0] / lengths[0]); d[1] -= lengths[1] * round(d[1] / lengths[1])
        r = np.linalg.norm(d)
        pair = {a, b}
        metal_oxygen = bool(pair & CATIONS) and bool(pair & {"O", "Oh"}) and r < 2.6
        hydroxyl = pair == {"Oh", "Hh"} and r < 1.15
        if metal_oxygen or hydroxyl:
            out.append((i, j, d))
    return out


def cross_section(ax, lab, xyz, lengths, slab_thickness, title, zmin):
    import matplotlib.colors as mc
    sel = set(np.where((xyz[:, 0] % lengths[0]) < slab_thickness)[0].tolist())
    # Draw hydroxyls whole: pull in the O or H partner of any hydroxyl atom inside the slice.
    oh_idx = np.where(np.isin(lab, ["Oh", "Hh"]))[0]
    if len(oh_idx):
        from scipy.spatial import cKDTree
        p = xyz[oh_idx].copy(); p[:, 0] %= lengths[0]; p[:, 1] %= lengths[1]; p[:, 2] += 1e4
        tree = cKDTree(p, boxsize=[lengths[0], lengths[1], 1e6])
        for a, b in tree.query_pairs(1.15):
            i, j = oh_idx[a], oh_idx[b]
            if {lab[i], lab[j]} == {"Oh", "Hh"} and (i in sel or j in sel):
                sel.update((int(i), int(j)))
    sel = np.array(sorted(sel))
    s = set(sel.tolist())
    order = sel[np.argsort(-xyz[sel, 0])]  # far (large x) first
    depth = np.clip((xyz[:, 0] % lengths[0]) / slab_thickness, 0.0, 1.0)
    for i, j, d in bonds(lab[sel], xyz[sel], lengths):
        gi, gj = sel[i], sel[j]
        y0, z0 = xyz[gi, 1], xyz[gi, 2]; y1, z1 = y0 + d[1], z0 + d[2]
        shade = 0.35 + 0.45 * (depth[gi] + depth[gj]) / 2
        ax.plot([y0, y1], [z0, z1], color=str(min(0.75, 0.25 + 0.5 * (depth[gi] + depth[gj]) / 2)), lw=1.8, zorder=1, solid_capstyle="round")
    for i in order:
        c = np.array(mc.to_rgb(COLORS[lab[i]])); c = c * (1 - 0.45 * depth[i]) + 0.45 * depth[i] * np.ones(3)
        ax.add_patch(__import__("matplotlib.patches", fromlist=["Circle"]).Circle(
            (xyz[i, 1], xyz[i, 2]), RADII[lab[i]], facecolor=c, edgecolor="#222222", lw=0.4, zorder=2 + (1 - depth[i])))
    ax.set_xlim(0, lengths[1]); ax.set_ylim(zmin, xyz[:, 2].max() + 2.0)
    ax.set_aspect("equal"); ax.set_title(title, fontsize=9)
    ax.set_xlabel("y (Å)", fontsize=8); ax.set_ylabel("z (Å)", fontsize=8); ax.tick_params(labelsize=7)


def top_view(ax, lab, xyz, lengths, title):
    """Surface height map (grey) with hydroxyl O (orange) and their H (white)."""
    nx, ny = int(lengths[0] / 0.5), int(lengths[1] / 0.5)
    hm = np.full((ny, nx), -np.inf)
    lattice = np.isin(lab, ["In", "Sn", "Ni", "O"])
    ix = (xyz[lattice, 0] % lengths[0] / lengths[0] * nx).astype(int) % nx
    iy = (xyz[lattice, 1] % lengths[1] / lengths[1] * ny).astype(int) % ny
    np.maximum.at(hm, (iy, ix), xyz[lattice, 2])
    from scipy.ndimage import maximum_filter
    hm = maximum_filter(hm, size=5, mode="wrap")
    im = ax.imshow(hm, origin="lower", extent=[0, lengths[0], 0, lengths[1]], cmap="Greys_r", aspect="auto")
    oh = lab == "Oh"
    ax.scatter(xyz[oh, 0] % lengths[0], xyz[oh, 1] % lengths[1], s=9, c=COLORS["Oh"], edgecolors="#222222", linewidths=0.3)
    ax.set_title(title, fontsize=9); ax.set_xlabel("x (Å), groove axis", fontsize=8); ax.set_ylabel("y (Å)", fontsize=8)
    ax.tick_params(labelsize=7)
    return im


def perspective(ax, lab, xyz, lengths, xmax, zmin, title):
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    m = ((xyz[:, 0] % lengths[0]) < xmax) & (xyz[:, 2] > zmin)
    idx = np.where(m)[0]
    segs = [[xyz[idx[i]], xyz[idx[i]] + d] for i, j, d in bonds(lab[idx], xyz[idx], lengths)]
    ax.add_collection3d(Line3DCollection(segs, colors="#7a7a7a", linewidths=0.6))
    for L in ("In", "Sn", "Ni", "O", "Oh", "Hh"):
        s = idx[lab[idx] == L]
        if len(s):
            ax.scatter(xyz[s, 0], xyz[s, 1], xyz[s, 2], s=(RADII[L] * 22) ** 2 / 10, c=COLORS[L], edgecolors="#222222",
                       linewidths=0.3, depthshade=True)
    ax.set_box_aspect((xmax, lengths[1], xyz[idx, 2].max() - zmin))
    ax.view_init(elev=24, azim=-62); ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])


def main(argv=None) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path); ap.add_argument("dirs", type=Path, nargs="+")
    ap.add_argument("--labels", nargs="*"); ap.add_argument("--ncols", type=int, default=5)
    ap.add_argument("--perspective", action="store_true"); ap.add_argument("--title", default="")
    ap.add_argument("--molecules", nargs="*", type=Path, help="one .xyz per DIR to overlay on its cross-section (use - for none)")
    ap.add_argument("--no-top-view", action="store_true")
    a = ap.parse_args(argv)
    n = len(a.dirs); ncols = min(a.ncols, n); blocks = int(np.ceil(n / ncols)); rows_per = 1 if a.no_top_view else 2
    fig = plt.figure(figsize=(4.4 * ncols, (3.6 if a.no_top_view else 5.6) * blocks), dpi=150)
    gs = fig.add_gridspec(rows_per * blocks, ncols, height_ratios=([1.0] if a.no_top_view else [1.35, 0.75]) * blocks, hspace=0.5, wspace=0.22)
    present = set()
    for k, d in enumerate(a.dirs):
        man, lab, xyz, lengths = load(d); present |= set(lab.tolist())
        h = man.get("hydroxylation", {})
        oh = h.get("oh_groups_per_nm2", 0.0); frac = h.get("achieved_fraction_of_exposed")
        label = (a.labels[k] if a.labels else man["model_id"])
        sub = f"{label}\n{h.get('pairs', 0)} H2O dissociated · {oh:.2f} OH/nm²" + (f" · {frac:.0%} of exposed" if frac is not None else "")
        thick = 2.2 if "Ni" in set(lab.tolist()) else 3.5
        r0 = (k // ncols) * rows_per
        ax = fig.add_subplot(gs[r0, k % ncols])
        zmin = xyz[:, 2].max() - 22.0
        mol = a.molecules[k] if a.molecules and k < len(a.molecules) and str(a.molecules[k]) != "-" else None
        if mol is not None:
            mxyz = np.loadtxt(mol, skiprows=2, usecols=(1, 2, 3))
            thick_center = float(np.mean(mxyz[:, 0]) % lengths[0])
            shifted = xyz.copy(); shifted[:, 0] = (xyz[:, 0] - (thick_center - thick / 2)) % lengths[0]
            cross_section(ax, lab, shifted, lengths, thick, sub, zmin)
            top = draw_molecule(ax, mol, lengths)
            ax.set_ylim(zmin, max(xyz[:, 2].max(), top[:, 2].max()) + 2.0)
        else:
            cross_section(ax, lab, xyz, lengths, thick, sub, zmin)
        if not a.no_top_view:
            top_view(fig.add_subplot(gs[r0 + 1, k % ncols]), lab, xyz, lengths, "top view: height (grey) + hydroxyl O (orange)")
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor=COLORS[L], markeredgecolor="#222222", markersize=8, label=NAMES[L])
               for L in ("In", "Sn", "Ni", "O", "Oh", "Hh") if L in present]
    if a.molecules and any(str(m) != "-" for m in a.molecules):
        handles += [Line2D([], [], marker="o", ls="", markerfacecolor=MOL_COLORS[e], markeredgecolor="#111111", markersize=7, label=f"{e} (SAM)")
                    for e in ("C", "N", "P")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=len(handles), frameon=False, fontsize=9)
    if a.title: fig.suptitle(a.title, fontsize=12)
    fig.savefig(a.out, bbox_inches="tight")
    print(a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
