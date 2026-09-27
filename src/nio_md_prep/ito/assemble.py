"""Assemble ITO/SAM pilot systems and their LAMMPS stage inputs.

Reuses the NiO workflow's LigParGen parser/replicator, Packmol wrapper and
Cao phosphonate corrections unchanged; replaces only the substrate, the
slab force field, and the thermostat/barostat protocol (rigid slab, NVT).
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import subprocess
import tempfile
from copy import deepcopy
from pathlib import Path

from ..build import _coeff_lines, _packmol, _surface
from ..chemistry import correction_lines, molecular_weight
from ..config import ROOT, load, molecule_manifest, missing_ligpargen
from ..geometry import write_xyz
from ..lammps import DataFile, TOPOLOGY, charge, parse, replicate, write
from .forcefield import HEADER, PARAMETER_SETS, surface_pair_lines

DEFAULTS = {
    "temperature": 300.0,
    "timestep": 0.5,
    "deposition_steps": 200000,
    "hold_steps": 200000,
    "release_steps": 100000,
    "relax_steps": 200000,
    "wall_clearance": 30.0,
    "release_height": 60.0,
    "dump_every": 1000,
    "tdamp": 100.0,
    # SHAKE on ligand X-H bonds; applied to dynamics only (never to minimization).
    "shake": False,
}


def _strip_pair_style(line: str) -> str:
    """Hybrid-style ``pair_coeff i j lj/cut/coul/long e s`` -> plain form."""
    return re.sub(r"^(pair_coeff\s+\d+\s+\d+\s+)lj/cut/coul/long\s+", r"\1", line)


def _init(source: str, shake: str = "") -> str:
    return f"""boundary p p f
processors * * 1
units real
atom_style full
read_data {source}
include force_field.lmp
group slab molecule 0
group mobile subtract all slab
neigh_modify exclude group slab slab
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
compute tmobile mobile temp
thermo 1000
thermo_style custom step c_tmobile pe ke etotal evdwl ecoul elong fnorm fmax
thermo_modify temp tmobile
fix freeze slab setforce 0.0 0.0 0.0
velocity slab set 0.0 0.0 0.0
{shake}"""


def _stage_inputs(p: dict, surface_top: float, zhi: float, velocity_seed: int, deposit_zstart: float,
                  shake_bond_types: list[int] | None = None) -> dict[str, str]:
    t = float(p["temperature"]); dt = float(p["timestep"])
    zend = surface_top + float(p["wall_clearance"])
    zrel = zend + float(p["release_height"])
    if zrel > zhi - 5.0:
        raise ValueError("release wall would leave the box; raise zhi")
    ns, nh, nr, nx = (int(p[k]) for k in ("deposition_steps", "hold_steps", "release_steps", "relax_steps"))
    dump = int(p["dump_every"]); tdamp = float(p["tdamp"])
    shake = ""
    if p.get("shake"):
        if not shake_bond_types:
            raise ValueError("shake requested but no ligand X-H bond types were found")
        shake = f"fix constrain mobile shake 0.0001 20 0 b {' '.join(map(str, shake_bond_types))}\n"
    deposit = f"""# ITO pilot: minimization + moving-wall deposition (rigid slab, NVT ligands)
{_init('topology_output.lmp')}fix walllo mobile wall/lj93 zlo EDGE 1.0 1.0 2.5 units box
min_style sd
min_modify dmax 0.005 line backtrack
minimize 0.0 10.0 5000 50000
min_style cg
min_modify dmax 0.01 line backtrack
minimize 0.0 1.0 20000 200000
print "optimization pe=$(pe) fnorm=$(fnorm) fmax=$(fmax)" file optimization-summary.txt screen yes
write_data optimized.data nocoeff
reset_timestep 0
{shake}variable zstart equal {deposit_zstart:.4f}
variable zend equal {zend:.4f}
variable zwall equal "v_zstart-(v_zstart-v_zend)*(step/{ns}.0)"
velocity mobile create {t} {velocity_seed} mom yes rot yes dist gaussian
fix wall mobile wall/lj126 zhi v_zwall 1.0 1.0 2.5 units box
fix ensemble mobile nvt temp {t} {t} {tdamp}
thermo_style custom step c_tmobile pe ke etotal evdwl ecoul elong v_zwall fnorm fmax
dump trajectory all custom {dump} deposition.lammpstrj id mol type q x y z
dump_modify trajectory sort id
timestep {dt}
run {ns} start 0 stop {ns}
undump trajectory
write_data deposited.data nocoeff
"""
    hold = f"""# ITO pilot: compressed-film hold at the deposition endpoint
{_init('deposited.data', shake)}fix walllo mobile wall/lj93 zlo EDGE 1.0 1.0 2.5 units box
fix wall mobile wall/lj126 zhi {zend:.4f} 1.0 1.0 2.5 units box
fix ensemble mobile nvt temp {t} {t} {tdamp}
dump trajectory all custom {dump} hold.lammpstrj id mol type q x y z
dump_modify trajectory sort id
timestep {dt}
run {nh}
undump trajectory
write_data held.data nocoeff
"""
    relax = f"""# ITO pilot: gradual wall release then relaxed-film hold
{_init('held.data', shake)}fix walllo mobile wall/lj93 zlo EDGE 1.0 1.0 2.5 units box
reset_timestep 0
variable zwall equal "{zend:.4f}+({zrel:.4f}-{zend:.4f})*(step/{nr}.0)"
fix wall mobile wall/lj126 zhi v_zwall 1.0 1.0 2.5 units box
fix ensemble mobile nvt temp {t} {t} {tdamp}
dump trajectory all custom {dump} release.lammpstrj id mol type q x y z
dump_modify trajectory sort id
timestep {dt}
run {nr} start 0 stop {nr}
undump trajectory
unfix wall
fix wall mobile wall/lj126 zhi {zrel:.4f} 1.0 1.0 2.5 units box
dump trajectory all custom {dump} relax.lammpstrj id mol type q x y z
dump_modify trajectory sort id
run {nx}
undump trajectory
write_data relaxed.data nocoeff
"""
    return {"deposition": deposit, "hold": hold, "relax": relax, "_zend": zend, "_zrel": zrel}


def _minimum_ligand_slab_distance(data: DataFile, lx: float, ly: float, cutoff: float = 3.0) -> float:
    import numpy as np
    from scipy.spatial import cKDTree
    atoms = data.sections["Atoms"]
    xyz = np.array([[float(a.fields[4]), float(a.fields[5]), float(a.fields[6])] for a in atoms])
    mol = np.array([int(a.fields[1]) for a in atoms])
    box = [lx, ly, 1e6]
    slab = xyz[mol == 0].copy(); lig = xyz[mol > 0].copy()
    for arr in (slab, lig):
        arr[:, 0] %= lx; arr[:, 1] %= ly; arr[:, 2] += 1e5
    d, _ = cKDTree(slab, boxsize=box).query(lig, distance_upper_bound=cutoff)
    return float(min(d.min(), cutoff))


def build_pilot(study_path: Path, output: Path, packmol_seed: int | None = None, velocity_seed: int | None = None,
                parameter_set: str | None = None, run_lammps: bool = True) -> Path:
    cfg = load(study_path)
    output.mkdir(parents=True, exist_ok=True)
    rnd = cfg.get("random", {})
    pseed = int(rnd.get("packmol_seed", 1) if packmol_seed is None else packmol_seed)
    vseed = int(rnd.get("velocity_seed", 1) if velocity_seed is None else velocity_seed)
    cfg["_resolved_random"] = {"packmol_seed": pseed, "velocity_seed": vseed}
    pset = parameter_set or cfg["forcefield"]["surface_parameter_set"]
    if pset not in PARAMETER_SETS:
        raise ValueError(f"unknown surface parameter set {pset!r}")
    substrate_dir = ROOT / cfg["substrate"]["path"]
    smanifest = json.loads((substrate_dir / "surface_manifest.json").read_text(encoding="utf-8"))
    surface_path = substrate_dir / "surface.lmp"
    if hashlib.sha256(surface_path.read_bytes()).hexdigest() != smanifest["surface_lmp_sha256"]:
        raise ValueError(f"{surface_path} does not match its surface_manifest.json hash")
    surface = parse(surface_path)
    lx = surface.bounds["x"][1] - surface.bounds["x"][0]; ly = surface.bounds["y"][1] - surface.bounds["y"][0]
    zlo, zhi = surface.bounds["z"]
    zhi = float(cfg.get("box", {}).get("zhi", zhi))
    top = float(smanifest["z_top_atom_angstrom"])
    pk = cfg.get("packing", {})
    z0 = top + float(pk.get("gap", 8.0)); z1 = z0 + float(pk.get("height", 60.0))
    inset = float(pk.get("lateral_inset", 1.0))
    region = f"{inset:.4f} {inset:.4f} {z0:.4f} {lx-inset:.4f} {ly-inset:.4f} {z1:.4f}"
    templates = []
    for spec in cfg["molecules"]:
        folder, manifest = molecule_manifest(spec["slug"])
        lmp = folder / manifest.get("files", {}).get("ligpargen", "ligpargen.lmp")
        if not lmp.exists():
            raise missing_ligpargen(lmp, f"python -m nio_md_prep.ito build-pilot {study_path} --output {output}")
        data = parse(lmp)
        expected = manifest["molecule"]["expected_net_charge"]
        if abs(charge(data) - float(expected)) > 1e-6:
            raise ValueError(f"{lmp}: charge {charge(data):.8f} != manifest expected_net_charge {expected}")
        xyz = output / f"{spec['slug']}.xyz"; write_xyz(data, xyz, spec["slug"])
        templates.append({"slug": spec["slug"], "data": data, "xyz": xyz, "count": int(spec["count"]), "atoms": data.count("Atoms"),
                          "mw": molecular_weight(data), "region": region, "folder": folder, "manifest": manifest})
    packed, coords, ok = _packmol(cfg, templates, output)
    if coords is None:
        raise RuntimeError("Packmol did not produce packed.xyz; install packmol or supply coordinates")
    result = DataFile("LAMMPS data file generated by nio_md_prep.ito", {}, {})
    ids = {"atom": 0, "bond": 0, "angle": 0, "dihedral": 0, "improper": 0}; types = {k: 0 for k in ids}; mol = 0
    manifest_out: dict = {"components": [], "packmol_seed": pseed, "velocity_seed": vseed}
    ff_lines: list[str] = []
    cursor = 0
    for t in templates:
        n = t["atoms"] * t["count"]; local = coords[cursor:cursor + n]; cursor += n
        pieces, inc = replicate(t["data"], t["count"], types, ids, mol, local)
        for sec, rows in pieces.items():
            result.sections.setdefault(sec, []).extend(rows)
        chem, anchors, corrected = correction_lines(t["data"], types["atom"], types["dihedral"])
        if anchors != int(t["manifest"]["molecule"].get("phosphonic_acid_anchors", 0)):
            raise ValueError(f"{t['slug']}: found {anchors} phosphonate anchors")
        ff_lines += [_strip_pair_style(x) for x in chem]
        manifest_out["components"].append({
            "component": t["slug"], "count": t["count"], "atoms_per_molecule": t["atoms"], "molecular_weight_g_mol": round(t["mw"], 6),
            "charge": charge(t["data"]) * t["count"], "atom_ids": [ids["atom"] + 1, ids["atom"] + inc["atom"]],
            "molecule_ids": [mol + 1, mol + inc["molecule"]],
            "types": {c: [types[c] + 1, types[c] + t["data"].type_count(c)] for c in types if t["data"].type_count(c)},
            "corrected_charmm_dihedral_types": sorted(corrected)})
        for k in ids:
            ids[k] += inc.get(k, 0); types[k] += t["data"].type_count(k)
        mol += inc["molecule"]
    srows, sinc = _surface(surface, types, ids, 0)
    for sec, rows in srows.items():
        result.sections.setdefault(sec, []).extend(rows)
    old_to_new = {int(r.fields[0]): types["atom"] + i + 1 for i, r in enumerate(surface.sections["Masses"])}
    label_to_type = {label: old_to_new[old] for label, old in smanifest["type_ids"].items()}
    manifest_out["surface"] = {"component": smanifest["model_id"], "atom_ids": [ids["atom"] + 1, ids["atom"] + sinc["atom"]],
                               "molecule_ids": [0, 0], "label_to_type": label_to_type, "charge": charge(surface)}
    result.bounds = deepcopy(surface.bounds)
    result.bounds["z"] = (zlo, zhi)
    if max(float(a.fields[6]) for a in result.sections["Atoms"]) > zhi - 5.0:
        raise ValueError(f"packed atoms reach within 5 A of box zhi={zhi}")
    masses = {int(r.fields[0]): float(r.fields[1]) for r in result.sections["Masses"]}
    atom_type = {int(a.fields[0]): int(a.fields[2]) for a in result.sections["Atoms"]}
    shake_types = sorted({int(b.fields[1]) for b in result.sections.get("Bonds", [])
                          if any(abs(masses[atom_type[int(i)]] - 1.008) < 0.01 for i in b.fields[2:4])})
    topo = deepcopy(result)
    for sec in ("Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs"):
        topo.sections.pop(sec, None)
    write(topo, output / "topology_output.lmp")
    ligand_coeffs = [_strip_pair_style(x) for x in _coeff_lines(result)]
    ff = HEADER + "\n" + "\n".join(ligand_coeffs + ff_lines + surface_pair_lines(label_to_type, pset)) + "\n"
    (output / "force_field.lmp").write_text(ff, encoding="utf-8")
    p = DEFAULTS | cfg.get("protocol", {})
    ligand_zmax = max(float(a.fields[6]) for a in result.sections["Atoms"] if int(a.fields[1]) > 0)
    zstart = ligand_zmax + 4.0
    stages = _stage_inputs(p, top, zhi, vseed, zstart, shake_types)
    for name in ("deposition", "hold", "relax"):
        (output / f"{name}.in").write_text(stages[name], encoding="utf-8")
    min_sep = _minimum_ligand_slab_distance(result, lx, ly)
    area = lx * ly / 100.0
    manifest_out.update({
        "status": "assembled", "study": str(study_path.resolve().relative_to(ROOT)) if study_path.resolve().is_relative_to(ROOT) else str(study_path),
        "substrate": {"path": cfg["substrate"]["path"], "model_id": smanifest["model_id"], "surface_lmp_sha256": smanifest["surface_lmp_sha256"],
                      "area_nm2": smanifest["area_nm2"], "z_top_atom_angstrom": top},
        "surface_parameter_set": pset, "surface_parameter_sources": PARAMETER_SETS[pset]["sources"],
        "packing_region": region, "areal_dose_molecules_per_nm2": round(sum(t["count"] for t in templates) / area, 4),
        "protocol": {k: p[k] for k in DEFAULTS}, "shake_bond_types": shake_types if p.get("shake") else [], "deposition_wall": {"zstart": zstart, "zend": stages["_zend"], "release": stages["_zrel"]},
        "total_charge": charge(result), "counts": {s.lower(): result.count(s) for s in ("Atoms", "Bonds", "Angles", "Dihedrals", "Impropers")},
        "box": result.bounds, "minimum_ligand_slab_distance_angstrom": round(min_sep, 4),
        "model_scope": "classical-ff: rigid slab, fixed charges, protonated neutral phosphonic acids; physisorption/H-bonding only",
        "source_hashes": {str(x.relative_to(ROOT)): hashlib.sha256(x.read_bytes()).hexdigest()
                          for x in [surface_path, substrate_dir / "model.toml"] + [t["folder"] / "molecule.toml" for t in templates]
                          + [t["folder"] / t["manifest"]["files"]["ligpargen"] for t in templates]},
    })
    errors = []
    expected_total = sum(charge(t["data"]) * t["count"] for t in templates) + charge(surface)
    if abs(manifest_out["total_charge"] - expected_total) > 1e-6 or abs(manifest_out["total_charge"]) > 0.05:
        errors.append(f"total charge {manifest_out['total_charge']:.6f} disagrees with components or is not ~0")
    if min_sep < 2.0:
        errors.append(f"ligand-slab minimum distance {min_sep:.3f} A < 2.0 A")
    report = [f"atoms: {result.count('Atoms')}", f"molecules: {mol}", f"total charge: {manifest_out['total_charge']:.6f}",
              f"ligand-slab minimum distance: {min_sep:.3f} A", f"surface parameter set: {pset}"]
    if run_lammps:
        exe = shutil.which("lmp") or shutil.which("lmp_serial") or shutil.which("lmp_mpi")
        if exe is None:
            report.append("WARNING: LAMMPS not found; zero-step checks skipped")
        else:
            with tempfile.TemporaryDirectory(prefix="ito-validate-") as tmp:
                for f in ("topology_output.lmp", "force_field.lmp"):
                    shutil.copy2(output / f, Path(tmp) / f)
                # Every run/minimize becomes a zero-step force evaluation of the real stage.
                smoke = re.sub(r"(?m)^(run|minimize)\s+.*$", "run 0 post no", stages["deposition"])
                smoke = re.sub(r"(?m)^(dump|dump_modify|undump|write_data|print).*$", "", smoke)
                (Path(tmp) / "smoke.in").write_text(smoke, encoding="utf-8")
                run = subprocess.run([exe, "-in", "smoke.in", "-log", "smoke.log"], cwd=tmp, capture_output=True, text=True)
                log = (Path(tmp) / "smoke.log").read_text(errors="replace") if (Path(tmp) / "smoke.log").exists() else run.stdout
                if run.returncode:
                    errors.append("LAMMPS zero-step deposition check failed: " + (run.stdout + run.stderr)[-800:])
                else:
                    m = re.findall(r"^\s*0\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", log, re.M)
                    report.append("LAMMPS zero-step deposition check: passed")
                    if m:
                        manifest_out["zero_step_thermo"] = dict(zip(["temp", "pe", "ke", "etotal", "evdwl", "ecoul", "elong", "fnorm", "fmax"], map(float, m[-1])))
                        if not all(math.isfinite(v) for v in manifest_out["zero_step_thermo"].values()):
                            errors.append("non-finite zero-step energy")
    status = "PASS" if not errors else "FAIL"
    (output / "validation_report.txt").write_text("\n".join([status] + report + [f"ERROR: {e}" for e in errors]) + "\n", encoding="utf-8")
    manifest_out["validation"] = status
    manifest_out["output_hashes"] = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(output.iterdir())
                                     if f.is_file() and f.suffix in (".lmp", ".in", ".xyz")}
    (output / "assembly_manifest.json").write_text(json.dumps(manifest_out, indent=2) + "\n", encoding="utf-8")
    if errors:
        raise ValueError("; ".join(errors))
    return output
