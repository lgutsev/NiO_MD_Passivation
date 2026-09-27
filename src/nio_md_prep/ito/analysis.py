"""Contacts, orientation, clustering and coverage for ITO pilot trajectories.

Structural, classical-FF observables only (see ``analysis/model_scope.py``):

* anchoring: phosphonate O within ``cation_cutoff`` of a slab In/Sn, and
  anchor H-bonds (acidic H to slab O, or phosphonate O to slab hydroxyl H,
  shorter than ``hbond_cutoff``).  "Contact" means proximity in a
  nonreactive model, never a P-O-In bond;
* P height above the *local* rigid-slab surface (max slab z within 2.5 A
  laterally; on flat slabs this equals the height above the top atom), tilt of the P -> core vector
  from the surface normal (core = farthest 30 % of C/N atoms, as in
  ``checks.molecule_frame``);
* clustering: periodic single-linkage of P heads of surface-resident
  molecules in x/y (``cluster_cutoff``), plus the fraction of molecules
  stranded above the first layer;
* coverage: delegated unchanged to ``analysis.coverage.analyze_coverage``;
* corrugated slabs: surface-resident molecules are binned into groove floor,
  groove walls and plateau (from the slab manifest's ``corrugation`` block),
  giving projected-area densities and the groove/plateau enrichment ratio.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ..analysis.coverage import analyze_coverage, iter_dump_frames
from ..config import ROOT, molecule_manifest
from ..lammps import parse
from .checks import molecule_frame


def _periodic_clusters(xy: np.ndarray, lengths: tuple[float, float], cutoff: float) -> list[int]:
    n = len(xy)
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for i in range(n):
        d = xy[i + 1:] - xy[i]
        for ax in (0, 1):
            d[:, ax] -= lengths[ax] * np.round(d[:, ax] / lengths[ax])
        for j in np.where(np.hypot(d[:, 0], d[:, 1]) < cutoff)[0] + i + 1:
            ri, rj = find(i), find(int(j))
            if ri != rj: parent[ri] = rj
    sizes: dict[int, int] = {}
    for i in range(n):
        r = find(i); sizes[r] = sizes.get(r, 0) + 1
    return sorted(sizes.values(), reverse=True)


def _min_periodic(a: np.ndarray, b: np.ndarray, lengths) -> np.ndarray:
    """Per-row minimum distance from each point in ``a`` to the set ``b``."""
    from scipy.spatial import cKDTree
    if len(b) == 0 or len(a) == 0:
        return np.full(len(a), np.inf)
    zshift = 1e5
    bb = b.copy(); aa = a.copy()
    for arr in (aa, bb):
        arr[:, 0] %= lengths[0]; arr[:, 1] %= lengths[1]; arr[:, 2] += zshift
    d, _ = cKDTree(bb, boxsize=[lengths[0], lengths[1], 1e7]).query(aa)
    return d


def _local_surface_height(slab_xyz: np.ndarray, lengths, radius: float = 2.5, spacing: float = 0.5):
    """Callable giving the max slab z within ``radius`` (periodic x/y) of a point."""
    nx = max(1, int(round(lengths[0] / spacing))); ny = max(1, int(round(lengths[1] / spacing)))
    grid = np.full((nx, ny), -np.inf)
    ix = np.floor((slab_xyz[:, 0] % lengths[0]) / lengths[0] * nx).astype(int) % nx
    iy = np.floor((slab_xyz[:, 1] % lengths[1]) / lengths[1] * ny).astype(int) % ny
    np.maximum.at(grid, (ix, iy), slab_xyz[:, 2])
    rx = int(np.ceil(radius / (lengths[0] / nx))); ry = int(np.ceil(radius / (lengths[1] / ny)))
    local = np.full_like(grid, -np.inf)
    for dx in range(-rx, rx + 1):
        for dy in range(-ry, ry + 1):
            if (dx * lengths[0] / nx) ** 2 + (dy * lengths[1] / ny) ** 2 <= radius ** 2:
                local = np.maximum(local, np.roll(np.roll(grid, dx, axis=0), dy, axis=1))

    def height(x: float, y: float) -> float:
        return float(local[int(np.floor((x % lengths[0]) / lengths[0] * nx)) % nx, int(np.floor((y % lengths[1]) / lengths[1] * ny)) % ny])
    return height


def _groove_regions(corr: dict, lengths) -> dict:
    """Projected groove floor / wall / plateau bands (groove runs along x, profile along y)."""
    half = corr["depth_angstrom"] / corr["wall_slope"]; floor = corr["step_run_angstrom"] / 2.0
    ly = lengths[1]
    return {"center_y": corr["groove_center_y_angstrom"], "floor_half_width": floor, "opening_half_width": half,
            "area_nm2": {"groove_floor": 2 * floor * lengths[0] / 100.0, "groove_wall": 2 * (half - floor) * lengths[0] / 100.0,
                         "plateau": (ly - 2 * half) * lengths[0] / 100.0}}


def _region(y: float, regions: dict, ly: float) -> str:
    dy = y - regions["center_y"]; dy -= ly * round(dy / ly); dy = abs(dy)
    if dy < regions["floor_half_width"]: return "groove_floor"
    if dy < regions["opening_half_width"]: return "groove_wall"
    return "plateau"


def _find_summary_values(node, keys: set[str], out: dict, prefix: str = "") -> None:
    if isinstance(node, dict):
        for k, v in node.items():
            if k in keys and isinstance(v, dict):
                out[prefix + k] = {kk: vv for kk, vv in v.items() if isinstance(vv, (int, float))}
            else:
                _find_summary_values(v, keys, out, prefix)


def analyze(build_directory: Path, trajectory: Path, output: Path | None = None, *, last_frames: int = 20,
            cation_cutoff: float = 3.25, hbond_cutoff: float = 2.5, surface_layer_height: float = 6.0,
            cluster_cutoff: float = 7.0, run_coverage: bool = True) -> Path:
    build_directory = Path(build_directory); trajectory = Path(trajectory)
    output = Path(output) if output else build_directory / f"ito-analysis-{trajectory.stem}"
    output.mkdir(parents=True, exist_ok=True)
    man = json.loads((build_directory / "assembly_manifest.json").read_text(encoding="utf-8"))
    top = float(man["substrate"]["z_top_atom_angstrom"])
    label_of_type = {int(t): lab for lab, t in man["surface"]["label_to_type"].items()}
    # Per-component local role indices from the LigParGen template.
    comps = []
    for c in man["components"]:
        folder, mm = molecule_manifest(c["component"])
        fr = molecule_frame(parse(folder / mm["files"]["ligpargen"]))
        comps.append({"name": c["component"], "atom_lo": c["atom_ids"][0], "n": c["atoms_per_molecule"],
                      "count": c["count"], "mol_lo": c["molecule_ids"][0], "frame": fr})
    frames = list(iter_dump_frames(trajectory))
    if not frames:
        raise ValueError(f"{trajectory}: no frames")
    frames = frames[-last_frames:] if last_frames else frames
    smanifest_path = ROOT / man["substrate"]["path"] / "surface_manifest.json"
    corrugation = json.loads(smanifest_path.read_text(encoding="utf-8")).get("corrugation") if smanifest_path.is_file() else None
    f0 = frames[0]; m0 = f0.molecule_ids <= 0
    lengths0 = (f0.bounds[0][1] - f0.bounds[0][0], f0.bounds[1][1] - f0.bounds[1][0])
    # The slab is rigid, so one height map from the first analyzed frame serves every frame.
    local_height = _local_surface_height(np.column_stack([f0.x[m0], f0.y[m0], f0.z[m0]]), lengths0)
    regions = _groove_regions(corrugation, lengths0) if corrugation else None
    per_frame = []
    last_molecules = []
    for f in frames:
        order = np.argsort(f.atom_ids)
        ids = f.atom_ids[order]; types = f.atom_types[order]
        xyz = np.column_stack([f.x[order], f.y[order], f.z[order]])
        lengths = (f.bounds[0][1] - f.bounds[0][0], f.bounds[1][1] - f.bounds[1][0])
        index = {int(i): k for k, i in enumerate(ids)}
        slab_mask = f.molecule_ids[order] <= 0
        labels = np.array([label_of_type.get(int(t), "?") for t in types])
        cations = xyz[slab_mask & np.isin(labels, ["In", "Sn", "Ni"])]
        slab_o = xyz[slab_mask & np.isin(labels, ["O", "Oh"])]
        slab_h = xyz[slab_mask & (labels == "Hh")]
        mols = []
        for c in comps:
            fr = c["frame"]
            for m in range(c["count"]):
                base = c["atom_lo"] + m * c["n"]
                sel = np.array([index[base + k] for k in range(c["n"])])
                mx = xyz[sel].copy()
                # Unwrap the molecule around its P atom.
                ref = mx[fr["p"][0]]
                for ax in (0, 1):
                    mx[:, ax] -= lengths[ax] * np.round((mx[:, ax] - ref[ax]) / lengths[ax])
                p = mx[fr["p"]].mean(axis=0)
                core = mx[fr["core"]].mean(axis=0)
                v = core - p
                ao = mx[fr["anchor_o"]]; ah = mx[fr["acid_h"]] if fr["acid_h"] else np.empty((0, 3))
                d_cat = _min_periodic(ao, cations, lengths)
                hb = int((_min_periodic(ah, slab_o, lengths) < hbond_cutoff).sum()) + int((_min_periodic(ao, slab_h, lengths) < hbond_cutoff).sum())
                mols.append({"component": c["name"], "molecule": c["mol_lo"] + m, "p_x": float(p[0] % lengths[0]), "p_y": float(p[1] % lengths[1]),
                             "p_height": float(p[2] - local_height(p[0], p[1])), "p_height_above_top_atom": float(p[2] - top),
                             "region": _region(p[1] % lengths[1], regions, lengths[1]) if regions else "flat", "tilt_deg": float(math.degrees(math.acos(max(-1, min(1, v[2] / np.linalg.norm(v)))))),
                             "min_anchor_o_cation": float(d_cat.min()), "cation_contacts": int((d_cat < cation_cutoff).sum()), "anchor_hbonds": hb})
        surf = [m for m in mols if m["p_height"] <= surface_layer_height]
        sizes = _periodic_clusters(np.array([[m["p_x"], m["p_y"]] for m in surf]), lengths, cluster_cutoff) if surf else []
        row = {"step": f.step, "n_molecules": len(mols), "surface_resident_fraction": len(surf) / len(mols),
               "cation_contact_fraction": float(np.mean([m["cation_contacts"] > 0 for m in mols])),
               "anchor_hbond_fraction": float(np.mean([m["anchor_hbonds"] > 0 for m in mols])),
               "anchored_fraction": float(np.mean([(m["cation_contacts"] > 0) or (m["anchor_hbonds"] > 0) for m in mols])),
               "p_height_median": float(np.median([m["p_height"] for m in mols])),
               "tilt_median_surface_resident": float(np.median([m["tilt_deg"] for m in surf])) if surf else None,
               "cluster_count": len(sizes), "largest_cluster_fraction_of_surface": (sizes[0] / len(surf)) if surf else None,
               "surface_density_nm2": len(surf) / (lengths[0] * lengths[1] / 100.0)}
        if regions:
            dens = {}
            for reg, area in regions["area_nm2"].items():
                n_reg = sum(1 for m in surf if m["region"] == reg)
                row[f"{reg}_surface_resident"] = n_reg
                dens[reg] = n_reg / area if area > 0 else float("nan")
                row[f"{reg}_density_nm2"] = dens[reg]
            groove_area = regions["area_nm2"]["groove_floor"] + regions["area_nm2"]["groove_wall"]
            groove_density = (row["groove_floor_surface_resident"] + row["groove_wall_surface_resident"]) / groove_area
            row["groove_density_nm2"] = groove_density
            row["groove_to_plateau_density_ratio"] = groove_density / dens["plateau"] if dens["plateau"] > 0 else None
            for c in comps:
                cs = [m for m in surf if m["component"] == c["name"]]
                row[f"{c['name']}_groove_fraction_of_surface_resident"] = (
                    sum(1 for m in cs if m["region"] != "plateau") / len(cs)) if cs else None
        for c in comps:
            cm = [m for m in mols if m["component"] == c["name"]]
            row[f"{c['name']}_surface_resident_fraction"] = float(np.mean([m["p_height"] <= surface_layer_height for m in cm]))
            row[f"{c['name']}_anchored_fraction"] = float(np.mean([(m["cation_contacts"] > 0) or (m["anchor_hbonds"] > 0) for m in cm]))
        per_frame.append(row)
        last_molecules = mols
    keys = sorted({k for r in per_frame for k, v in r.items() if isinstance(v, (int, float)) and not isinstance(v, bool) and k != "step"})
    summary = {"build_directory": str(build_directory), "trajectory": str(trajectory), "frames": len(per_frame),
               "corrugation_regions": regions,
               "steps": [per_frame[0]["step"], per_frame[-1]["step"]],
               "parameters": {"cation_cutoff": cation_cutoff, "hbond_cutoff": hbond_cutoff,
                              "surface_layer_height": surface_layer_height, "cluster_cutoff": cluster_cutoff},
               "model_scope": man.get("model_scope"),
               "mean": {k: float(np.mean([r[k] for r in per_frame if r.get(k) is not None])) for k in keys},
               "tilt_histogram_surface_resident_last_frame": np.histogram(
                   [m["tilt_deg"] for m in last_molecules if m["p_height"] <= surface_layer_height], bins=range(0, 181, 15))[0].tolist()}
    if run_coverage:
        cov_path = analyze_coverage(build_directory, trajectory, output / "coverage", last_frames=last_frames)
        cov = json.loads(Path(cov_path).read_text(encoding="utf-8"))
        picked: dict = {}
        _find_summary_values(cov, {"total", "near_surface", "anchor_conditioned", "void_largest_patch_percent", "roughness_rms"}, picked)
        summary["coverage"] = picked
    (output / "per_frame.json").write_text(json.dumps(per_frame, indent=2) + "\n", encoding="utf-8")
    (output / "last_frame_molecules.json").write_text(json.dumps(last_molecules, indent=2) + "\n", encoding="utf-8")
    path = output / "ito_summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return path
