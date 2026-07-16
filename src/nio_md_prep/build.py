from __future__ import annotations
from pathlib import Path
import math
import hashlib, json, shutil, subprocess
from .config import ROOT, load, molecule_manifest, missing_ligpargen
from .geometry import read_xyz, write_xyz
from .lammps import DataFile, Record, TOPOLOGY, charge, parse, replicate, write

HEADER = """# Validated NiO/passivant force-field styles
bond_style      harmonic
angle_style     harmonic
dihedral_style  hybrid opls charmm
improper_style  cvff
pair_style      hybrid lj/cut/coul/long 10.0 8.0 buck/coul/long 10.0 8.0
pair_modify     pair lj/cut/coul/long mix geometric
kspace_style    pppm 1e-4
kspace_modify   slab 3.0
special_bonds   amber
"""

def _merge(target: DataFile, pieces: dict[str,list[Record]]) -> None:
    for sec, rows in pieces.items(): target.sections.setdefault(sec, []).extend(rows)

def _override_lines(path: Path, maps: dict[str,dict[int,int]]) -> list[str]:
    if not path.exists(): return []
    result=[]
    categories={"pair_coeff":"atom", "bond_coeff":"bond", "angle_coeff":"angle", "dihedral_coeff":"dihedral", "improper_coeff":"improper"}
    for raw in path.read_text(encoding="utf-8").splitlines():
        stripped=raw.strip()
        if not stripped or stripped.startswith("#"): result.append(raw); continue
        parts=stripped.split(); cat=categories.get(parts[0])
        if cat:
            indices=(1,2) if parts[0]=="pair_coeff" else (1,)
            for i in indices: parts[i]=str(maps[cat][int(parts[i])])
            result.append(" ".join(parts))
        else: result.append(raw)
    return result

def _coeff_lines(data: DataFile) -> list[str]:
    names={"Pair Coeffs":"pair_coeff", "Bond Coeffs":"bond_coeff", "Angle Coeffs":"angle_coeff", "Dihedral Coeffs":"dihedral_coeff", "Improper Coeffs":"improper_coeff"}
    out=[]
    for sec, command in names.items():
        for r in data.sections.get(sec, []):
            payload=r.fields.copy()
            if sec=="Dihedral Coeffs": payload.insert(1,"opls")
            out.append(f"{command} {' '.join(payload)}" + (f" # {r.comment}" if r.comment else ""))
    return out

def _packmol(config: dict, templates: list[dict], output: Path) -> tuple[Path,list[tuple[float,float,float]]|None]:
    lines=["tolerance 2.0", "filetype xyz", "output packed.xyz", "seed 202405367", ""]
    for item in templates:
        lines += [f"structure {item['xyz'].name}", f"  number {item['count']}", f"  inside box {item['region']}", "end structure", ""]
        destination=output/item["xyz"].name
        if item["xyz"].resolve()!=destination.resolve(): shutil.copy2(item["xyz"],destination)
    inp=output/"packmol.inp"; inp.write_text("\n".join(lines),encoding="utf-8")
    exe=shutil.which("packmol")
    if exe:
        with inp.open("rb") as source:
            run=subprocess.run([exe],cwd=output,stdin=source,capture_output=True,text=False)
        if run.returncode: raise RuntimeError(run.stderr.decode(errors="replace"))
    packed=output/"packed.xyz"
    if not packed.exists(): return packed,None
    _, coords=read_xyz(packed)
    expected=sum(i["atoms"]*i["count"] for i in templates)
    if len(coords)!=expected: raise ValueError(f"Packmol output has {len(coords)} atoms; expected {expected}")
    return packed,coords

def _surface(data: DataFile, offsets: dict[str,int], ids: dict[str,int], molecule: int) -> tuple[dict[str,list[Record]],dict[str,int]]:
    rows={s:[] for s in data.sections}
    type_map={int(r.fields[0]): offsets["atom"]+i+1 for i,r in enumerate(data.sections["Masses"])}
    for r in data.sections["Masses"]:
        f=r.fields.copy(); f[0]=str(type_map[int(f[0])]); rows["Masses"].append(Record(f,r.comment))
    for i,r in enumerate(data.sections["Atoms"]):
        f=r.fields.copy(); f[0]=str(ids["atom"]+i+1); f[1]=str(molecule); f[2]=str(type_map[int(f[2])]); rows["Atoms"].append(Record(f,r.comment))
    return rows,{"atom":len(data.sections["Atoms"]),"molecule":0}

def build(config_path: Path, output: Path, primary_final: Path|None=None) -> Path:
    cfg=load(config_path); output.mkdir(parents=True,exist_ok=True)
    protocol=cfg.get("study",{}).get("protocol","cosam")
    species=cfg.get("molecules",[])
    if not species or any("count" not in x for x in species): raise ValueError("every study molecule requires an explicit count")
    templates=[]; manifests=[]
    for spec in species:
        folder,manifest=molecule_manifest(spec["slug"]); lmp=folder/manifest.get("files",{}).get("ligpargen","ligpargen.lmp")
        if not lmp.exists(): raise missing_ligpargen(lmp,f"nio-md-prep build {config_path} --output {output}")
        data=parse(lmp)
        expected=manifest.get("molecule",{}).get("expected_net_charge")
        if not isinstance(expected,(int,float)):
            raise ValueError(f"{folder/'molecule.toml'}: expected_net_charge must be reviewed and numeric")
        if abs(charge(data)-float(expected))>1e-6:
            raise ValueError(f"{lmp}: LigParGen charge {charge(data):.8f} does not match manifest expected_net_charge {expected}")
        geometry_name=manifest.get("files",{}).get("geometry","")
        xyz=folder/geometry_name if geometry_name else output/f"{spec['slug']}.xyz"
        if geometry_name: xyz=folder/geometry_name
        else: write_xyz(data,xyz,spec["slug"])
        templates.append({"slug":spec["slug"],"data":data,"xyz":xyz,"count":int(spec["count"]),"atoms":data.count("Atoms"),"region":spec.get("region","2.0 2.0 70.0 123.1 39.7 260.0"),"folder":folder,"manifest":manifest})
        manifests.append(manifest)
    packed,coords=_packmol(cfg,templates,output)
    if coords is None:
        # A reproducible preview remains inspectable, but is deliberately not called usable.
        coords=[]
        for t in templates:
            _,one=read_xyz(t["xyz"]); coords.extend(one*t["count"])
    result=DataFile("LAMMPS data file generated by nio-md-prep",{}, {})
    ids={"atom":0,"bond":0,"angle":0,"dihedral":0,"improper":0}; types={k:0 for k in ids}; mol=0; cursor=0; type_map={"molecules":{},"surface":{}}
    overrides=[]
    for t in templates:
        n=t["atoms"]*t["count"]; local=coords[cursor:cursor+n]; cursor+=n
        pieces,inc=replicate(t["data"],t["count"],types,ids,mol,local); _merge(result,pieces)
        maps={cat:{i:i+types[cat] for i in range(1,t["data"].type_count(cat)+1)} for cat in types}
        type_map["molecules"][t["slug"]]=maps
        override=t["manifest"].get("files",{}).get("override","")
        if override: overrides += _override_lines(t["folder"]/override,maps)
        for k in ids: ids[k]+=inc.get(k,0); types[k]+=t["data"].type_count(k)
        mol+=inc["molecule"]
    surface_path=ROOT/cfg.get("surface",{}).get("data","inputs/surfaces/corrugated-nio-110/surface.lmp")
    surface=parse(primary_final if primary_final else surface_path)
    surface_rows,inc=_surface(surface,types,ids,mol); _merge(result,surface_rows)
    old_types=sorted({int(r.fields[2]) for r in surface.sections["Atoms"]}); new_types=sorted({int(r.fields[2]) for r in surface_rows["Atoms"]})
    type_map["surface"]={str(a):b for a,b in zip(old_types,new_types)}
    all_xyz=[(float(r.fields[4]),float(r.fields[5]),float(r.fields[6])) for r in result.sections["Atoms"]]
    margin=cfg.get("box",{}).get("margin",5.0)
    result.bounds={axis:(min(p[i] for p in all_xyz)-margin,max(p[i] for p in all_xyz)+margin) for i,axis in enumerate("xyz")}
    write(result,output/"topology_output.lmp")
    ff=HEADER+"\n"+"\n".join(_coeff_lines(result)+overrides)+"\n"
    ni,o=new_types[-2:]
    ff += f"pair_coeff {ni} {ni} buck/coul/long 0.000000e+00 1.000000e+00 0.000000e+00\n"
    ff += f"pair_coeff {ni} {o} buck/coul/long 17397.5 0.35 0.000000e+00\n"
    ff += f"pair_coeff {o} {o} buck/coul/long 524948.8 0.1490 643.23\n"
    # Validated pseudo-LJ surface values are used only to form ligand/surface
    # cross interactions; Ni/O self interactions remain Buckingham.
    ligand_pairs={int(r.fields[0]):(float(r.fields[1]),float(r.fields[2])) for r in result.sections.get("Pair Coeffs",[]) if int(r.fields[0]) not in (ni,o)}
    for surface_type,(seps,ssig),label in ((ni,(0.1,3.0),"Ni"),(o,(0.21,3.05),"O")):
        for ligand_type,(eps,sig) in sorted(ligand_pairs.items()):
            ff += f"pair_coeff {surface_type} {ligand_type} lj/cut/coul/long {math.sqrt(seps*eps):.3f} {math.sqrt(ssig*sig):.3f} # {label}-{ligand_type}\n"
    (output/"force_field_settings_lammps_with_header.lmp").write_text(ff,encoding="utf-8")
    _write_input(output,cfg)
    (output/"resolved_config.toml").write_text(config_path.read_text(encoding="utf-8"),encoding="utf-8")
    (output/"type_map.json").write_text(json.dumps(type_map,indent=2)+"\n",encoding="utf-8")
    hashes={}
    for t in templates:
        for p in (t["folder"]/"molecule.toml",t["folder"]/t["manifest"]["files"]["ligpargen"]): hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    hashes[str(surface_path.relative_to(ROOT)) if surface_path.is_relative_to(ROOT) else str(surface_path)]=hashlib.sha256(surface_path.read_bytes()).hexdigest()
    (output/"input_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n",encoding="utf-8")
    from .validate import validate
    validate(output, packmol_ran=packed.exists(), primary_final=primary_final)
    return output

def _write_input(output: Path,cfg:dict)->None:
    p=cfg.get("protocol",{}); steps=int(p.get("deposition_steps",300000)); temp=float(p.get("temperature",300.0)); zstart=float(p.get("wall_start",270)); zend=float(p.get("wall_end",70)); lower=float(p.get("lower_wall",-5))
    text=f"""# Generated from corrected successful examples
boundary p p f
units real
atom_style full
read_data topology_output.lmp
include force_field_settings_lammps_with_header.lmp
timestep 1.0
variable pressure equal 1.0
variable pressureDamp equal 500.0
variable zstart equal {zstart}
variable zend equal {zend}
variable zwall equal \"v_zstart - (v_zstart-v_zend)*(step/{steps}.0)\"
variable epsilon equal 1.0
variable sigma equal 1.0
variable cutoff equal 2.5
fix wall all wall/lj126 zhi v_zwall ${{epsilon}} ${{sigma}} ${{cutoff}}
fix walllo all wall/lj93 zlo {lower} ${{epsilon}} ${{sigma}} ${{cutoff}} units box
fix deposit all npt temp 5.0 {temp} 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy
run {steps}
unfix deposit
unfix wall
unfix walllo
write_data final.data
"""
    (output/"lammps.in").write_text(text,encoding="utf-8")
