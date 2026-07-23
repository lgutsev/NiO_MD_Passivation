from __future__ import annotations
from pathlib import Path
import math
from copy import deepcopy
import hashlib, json, re, shutil, subprocess
from .config import ROOT, load, molecule_dir, molecule_manifest, missing_ligpargen
from .geometry import read_xyz, write_xyz
from .lammps import DataFile, Record, TOPOLOGY, charge, parse, replicate, write
from .chemistry import molecular_weight, correction_lines

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
            if sec=="Pair Coeffs": payload=[payload[0],payload[0],"lj/cut/coul/long",*payload[1:]]
            elif sec=="Dihedral Coeffs": payload.insert(1,"opls")
            out.append(f"{command} {' '.join(payload)}" + (f" # {r.comment}" if r.comment else ""))
    return out

def _expected_symbols(templates: list[dict]) -> list[str]:
    expected=[]
    for item in templates:
        symbols,_=read_xyz(item["xyz"])
        expected.extend(symbols*item["count"])
    return expected

def _fit_packed_molecules(
    templates: list[dict],
    coords: list[tuple[float,float,float]],
    maximum_shift: float = 0.1,
) -> tuple[list[tuple[float,float,float]],list[dict]]:
    """Correct small Packmol boundary overshoots by rigid translation only."""
    corrected=list(coords); adjustments=[]; cursor=0; safety=1.0e-6
    for item in templates:
        bounds=list(map(float,item["region"].split()))
        for copy_index in range(item["count"]):
            start=cursor+copy_index*item["atoms"]; stop=start+item["atoms"]
            molecule=corrected[start:stop]; shifts=[]
            for axis in range(3):
                values=[xyz[axis] for xyz in molecule]
                lower,upper=bounds[axis],bounds[axis+3]
                if min(values)>=lower and max(values)<=upper:
                    shifts.append(0.0); continue
                allowed_low=lower+safety-min(values)
                allowed_high=upper-safety-max(values)
                if allowed_low>allowed_high:
                    raise ValueError(f"{item['slug']} molecule {copy_index+1} cannot fit inside its validated region")
                shift=max(allowed_low,min(0.0,allowed_high))
                if abs(shift)>maximum_shift:
                    raise ValueError(
                        f"{item['slug']} molecule {copy_index+1} requires a {abs(shift):.6f} A "
                        f"boundary correction; maximum allowed rigid shift is {maximum_shift:.6f} A"
                    )
                shifts.append(shift)
            if any(shifts):
                molecule=[tuple(value+shifts[axis] for axis,value in enumerate(xyz)) for xyz in molecule]
                corrected[start:stop]=molecule
                adjustments.append({"component":item["slug"],"copy":copy_index+1,"translation_angstrom":shifts})
        cursor+=item["atoms"]*item["count"]
    return corrected,adjustments

def _packmol(config: dict, templates: list[dict], output: Path, supplied: Path|None=None) -> tuple[Path,list[tuple[float,float,float]]|None,bool]:
    lines=["tolerance 2.0", "filetype xyz", "output packed.xyz", "seed 202405367", ""]
    for item in templates:
        lines += [f"structure {item['xyz'].name}", f"  number {item['count']}", f"  inside box {item['region']}", "end structure", ""]
        destination=output/item["xyz"].name
        if item["xyz"].resolve()!=destination.resolve(): shutil.copy2(item["xyz"],destination)
    inp=output/"packmol.inp"; inp.write_text("\n".join(lines),encoding="utf-8")
    exe=shutil.which("packmol") if supplied is None else None
    if supplied is not None:
        destination=output/"packed.xyz"
        if supplied.resolve()!=destination.resolve(): shutil.copy2(supplied,destination)
    elif exe:
        with inp.open("rb") as source:
            run=subprocess.run([exe],cwd=output,stdin=source,capture_output=True,text=False)
        if run.returncode: raise RuntimeError(run.stderr.decode(errors="replace"))
    packed=output/"packed.xyz"
    if not packed.exists(): return packed,None,False
    symbols, coords=read_xyz(packed)
    expected=sum(i["atoms"]*i["count"] for i in templates)
    if len(coords)!=expected: raise ValueError(f"Packmol output has {len(coords)} atoms; expected {expected}")
    expected_symbols=_expected_symbols(templates)
    for index,(actual,wanted) in enumerate(zip(symbols,expected_symbols),1):
        if actual!=wanted:
            raise ValueError(f"Packmol atom order mismatch at packed atom {index}: expected {wanted}, found {actual}; atoms must remain in LigParGen order")
    coords,adjustments=_fit_packed_molecules(templates,coords)
    if adjustments:
        packed.write_text(
            f"{len(coords)}\nPackmol output with validated rigid whole-molecule boundary corrections\n"
            + "".join(f"{symbol} {x:.8f} {y:.8f} {z:.8f}\n" for symbol,(x,y,z) in zip(symbols,coords)),
            encoding="utf-8",
        )
    (output/"packmol_adjustments.json").write_text(
        json.dumps({"method":"rigid_whole_molecule_translation","maximum_shift_angstrom":0.1,"adjustments":adjustments},indent=2)+"\n",
        encoding="utf-8",
    )
    cursor=0
    for item in templates:
        n=item["atoms"]*item["count"]; bounds=list(map(float,item["region"].split()))
        for x,y,z in coords[cursor:cursor+n]:
            if not(bounds[0]<=x<=bounds[3] and bounds[1]<=y<=bounds[4] and bounds[2]<=z<=bounds[5]):
                raise ValueError(f"{item['slug']} packed coordinate lies outside its validated region")
        cursor+=n
    return packed,coords,True

def _surface(data: DataFile, offsets: dict[str,int], ids: dict[str,int], molecule: int) -> tuple[dict[str,list[Record]],dict[str,int]]:
    rows={s:[] for s in data.sections}
    type_map={int(r.fields[0]): offsets["atom"]+i+1 for i,r in enumerate(data.sections["Masses"])}
    for r in data.sections["Masses"]:
        f=r.fields.copy(); f[0]=str(type_map[int(f[0])]); rows["Masses"].append(Record(f,r.comment))
    for i,r in enumerate(data.sections["Atoms"]):
        f=r.fields.copy(); f[0]=str(ids["atom"]+i+1); f[1]=str(molecule); f[2]=str(type_map[int(f[2])]); rows["Atoms"].append(Record(f,r.comment))
    return rows,{"atom":len(data.sections["Atoms"]),"molecule":0}

def build(config_path: Path, output: Path, primary_final: Path|None=None, packed_xyz: Path|None=None) -> Path:
    cfg=load(config_path); output.mkdir(parents=True,exist_ok=True)
    protocol=cfg.get("study",{}).get("protocol","cosam")
    species=cfg.get("molecules",[])
    if not species or any("count" not in x for x in species): raise ValueError("every study molecule requires an explicit resolved count")
    templates=[]; manifests=[]
    sequential_bounds=None; existing_max_z=None
    if primary_final:
        if not primary_final.exists(): raise FileNotFoundError(primary_final)
        existing=parse(primary_final); sequential_bounds=existing.bounds.copy()
        existing_max_z=max(float(a.fields[6]) for a in existing.sections["Atoms"])
        cfg["_sequential_primary_atom_count"]=existing.count("Atoms")
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
        region=spec.get("region","2.0 2.0 45.0 123.1 39.7 145.0")
        if primary_final:
            xlo,xhi=sequential_bounds["x"]; ylo,yhi=sequential_bounds["y"]; _,zhi=sequential_bounds["z"]
            zlo=existing_max_z+2.5
            if zlo>=zhi-5.0: raise ValueError("stage-1 box has insufficient space above the film for stage-2 packing")
            region=f"{xlo+2:.6f} {ylo+2:.6f} {zlo:.6f} {xhi-2:.6f} {yhi-2:.6f} {zhi-5:.6f}"
        templates.append({"slug":spec["slug"],"data":data,"xyz":xyz,"count":int(spec["count"]),"atoms":data.count("Atoms"),"mw":molecular_weight(data),"region":region,"folder":folder,"manifest":manifest})
        manifests.append(manifest)
    packed,coords,packed_ok=_packmol(cfg,templates,output,packed_xyz)
    if coords is None:
        plan={"status":"packing_incomplete","packmol_seed":202405367,"command":"packmol < packmol.inp","component_order":[t["slug"] for t in templates]+([] if primary_final else ["corrugated-nio-110"]),"molecules":[{"slug":t["slug"],"count":t["count"],"atoms_per_molecule":t["atoms"]} for t in templates]}
        (output/"assembly_manifest.json").write_text(json.dumps(plan,indent=2)+"\n",encoding="utf-8")
        (output/"validation_report.txt").write_text("INCOMPLETE\nPackmol has not run; no scientifically usable topology was produced.\nRun: packmol < packmol.inp\n",encoding="utf-8")
        return output
    if primary_final:
        result=parse(primary_final)
        primary_velocity_count=result.count("Velocities")
        if primary_velocity_count not in (0,result.count("Atoms")):
            raise ValueError("stage-1 data contains an incomplete Velocities section")
        preserve_primary_velocities=primary_velocity_count==result.count("Atoms")
        ids={"atom":max(int(r.fields[0]) for r in result.sections["Atoms"])}
        for sec in TOPOLOGY: ids[sec.lower()[:-1]]=max((int(r.fields[0]) for r in result.sections.get(sec,[])),default=0)
        types={k:result.type_count(k) for k in ids}; mol=max((int(r.fields[1]) for r in result.sections["Atoms"]),default=0)
    else:
        result=DataFile("LAMMPS data file generated by nio-md-prep",{}, {})
        preserve_primary_velocities=False
        ids={"atom":0,"bond":0,"angle":0,"dihedral":0,"improper":0}; types={k:0 for k in ids}; mol=0
    cursor=0; type_map={"components":[],"surface":{},"packmol_seed":202405367}
    if primary_final:
        existing_atoms=result.sections["Atoms"]
        type_map["components"].append({"component":"stage1_primary","count":len({int(r.fields[1]) for r in existing_atoms if int(r.fields[1])>0}),"atoms_per_molecule":None,"charge":charge(result),"atom_ids":[min(int(r.fields[0]) for r in existing_atoms),max(int(r.fields[0]) for r in existing_atoms)],"molecule_ids":[0,max(int(r.fields[1]) for r in existing_atoms)],"types":{cat:[1,result.type_count(cat)] for cat in types if result.type_count(cat)}})
    overrides=[]
    for t in templates:
        n=t["atoms"]*t["count"]; local=coords[cursor:cursor+n]; cursor+=n
        pieces,inc=replicate(t["data"],t["count"],types,ids,mol,local); _merge(result,pieces)
        if preserve_primary_velocities:
            result.sections["Velocities"].extend(Record([atom.fields[0],"0.0","0.0","0.0"]) for atom in pieces["Atoms"])
        maps={cat:{i:i+types[cat] for i in range(1,t["data"].type_count(cat)+1)} for cat in types}
        ranges={"component":t["slug"],"count":t["count"],"atoms_per_molecule":t["atoms"],"charge":charge(t["data"])*t["count"],"atom_ids":[ids["atom"]+1,ids["atom"]+inc["atom"]],"molecule_ids":[mol+1,mol+inc["molecule"]],"types":{cat:[types[cat]+1,types[cat]+t["data"].type_count(cat)] for cat in types if t["data"].type_count(cat)}}
        type_map["components"].append(ranges); type_map.setdefault("type_maps",{})[t["slug"]]=maps
        chemical_lines, anchors, corrected_dihedrals=correction_lines(t["data"],types["atom"],types["dihedral"])
        expected_anchors=int(t["manifest"].get("molecule",{}).get("phosphonic_acid_anchors",0))
        if anchors != expected_anchors: raise ValueError(f"{t['slug']}: found {anchors} phosphonate anchors; manifest expects {expected_anchors}")
        overrides += chemical_lines
        ranges["corrected_charmm_dihedral_types"]=sorted(corrected_dihedrals)
        for k in ids: ids[k]+=inc.get(k,0); types[k]+=t["data"].type_count(k)
        mol+=inc["molecule"]
    surface_path=ROOT/cfg.get("surface",{}).get("data","inputs/surfaces/corrugated-nio-110/surface.lmp")
    if not primary_final:
        surface=parse(surface_path); surface_rows,inc=_surface(surface,types,ids,0); _merge(result,surface_rows)
        old_types=sorted({int(r.fields[2]) for r in surface.sections["Atoms"]}); new_types=sorted({int(r.fields[2]) for r in surface_rows["Atoms"]})
        type_map["surface"]={"component":"corrugated-nio-110","atom_ids":[ids["atom"]+1,ids["atom"]+inc["atom"]],"molecule_ids":[0,0],"atom_types":[min(new_types),max(new_types)],"original_to_final":{str(a):b for a,b in zip(old_types,new_types)},"charge":charge(surface)}
    else:
        new_types=sorted({int(r.fields[2]) for r in result.sections["Atoms"] if abs(abs(float(r.fields[3]))-2.0)<1e-8})
    all_xyz=[(float(r.fields[4]),float(r.fields[5]),float(r.fields[6])) for r in result.sections["Atoms"]]
    margin=cfg.get("box",{}).get("margin",5.0)
    if primary_final: result.bounds=sequential_bounds
    else:
        # Preserve the authoritative surface x/y periodic cell; only z contains vacuum.
        result.bounds={"x":surface.bounds["x"],"y":surface.bounds["y"],"z":surface.bounds["z"]}
    topology_data=deepcopy(result)
    for section in ("Pair Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs","Improper Coeffs"):
        topology_data.sections.pop(section,None)
    write(topology_data,output/"topology_output.lmp")
    previous_coeffs=[]
    if primary_final:
        previous_force_field=primary_final.parent/"force_field_settings.lmp"
        if not previous_force_field.exists():
            previous_force_field=primary_final.parent/"force_field_settings_lammps_with_header.lmp"
        if previous_force_field.exists():
            previous_coeffs=[line for line in previous_force_field.read_text(encoding="utf-8").splitlines() if line.strip().startswith(("pair_coeff","bond_coeff","angle_coeff","dihedral_coeff","improper_coeff"))]
    ff=HEADER+"\n"+"\n".join(previous_coeffs+_coeff_lines(result)+overrides)+"\n"
    ni,o=new_types[-2:]
    ff += f"pair_coeff {ni} {ni} buck/coul/long 0.000000e+00 1.000000e+00 0.000000e+00\n"
    ff += f"pair_coeff {ni} {o} buck/coul/long 17397.5 0.35 0.000000e+00\n"
    ff += f"pair_coeff {o} {o} buck/coul/long 524948.8 0.1490 643.23\n"
    # Validated pseudo-LJ surface values are used only to form ligand/surface
    # cross interactions; Ni/O self interactions remain Buckingham.
    ligand_pairs={int(r.fields[0]):(float(r.fields[1]),float(r.fields[2])) for r in result.sections.get("Pair Coeffs",[]) if int(r.fields[0]) not in (ni,o)}
    # Cross terms must use the corrected phosphonate self parameters too.
    for line in overrides:
        fields=line.split()
        if fields[:1]==["pair_coeff"] and fields[1]==fields[2]:
            ligand_pairs[int(fields[1])]=(float(fields[4]),float(fields[5]))
    for surface_type,(seps,ssig),label in ((ni,(0.1,3.0),"Ni"),(o,(0.21,3.05),"O")):
        for ligand_type,(eps,sig) in sorted(ligand_pairs.items()):
            ff += f"pair_coeff {surface_type} {ligand_type} lj/cut/coul/long {math.sqrt(seps*eps):.3f} {math.sqrt(ssig*sig):.3f} # {label}-{ligand_type}\n"
    (output/"force_field_settings.lmp").write_text(ff,encoding="utf-8")
    (output/"force_field_settings_lammps_with_header.lmp").write_text(ff,encoding="utf-8")
    _write_input(output,cfg,result)
    ratio=cfg.get("composition",{})
    report={"basis":ratio.get("basis","explicit"),"provisional":ratio.get("basis")=="equal_stock_volume","stock_concentration_mg_ml":ratio.get("stock_concentration_mg_ml",{}),"stock_volume_assumption":ratio.get("stock_volume_assumption","n/a"),"rounding_rule":ratio.get("rounding_rule","explicit counts"),"molecules":[{"slug":t["slug"],"molecular_weight_g_mol":round(t["mw"],6),"count":t["count"]} for t in templates]}
    if ratio.get("basis") in {"equal_stock_volume","molar_1_to_1"}:
        secondary=templates[-1]; primary_count=int(ratio.get("primary_count",templates[0]["count"])); _,pm=molecule_manifest("me-4pacz"); pd=parse(molecule_dir("me-4pacz")/pm["files"]["ligpargen"]); pmw=molecular_weight(pd)
        exact=1.0 if ratio.get("basis")=="molar_1_to_1" else (0.3/secondary["mw"])/(0.5/pmw)
        report.update({"unrounded_concentration_mw_ratio":exact,"rounded_count_ratio":secondary["count"]/primary_count,"primary_reference_count":primary_count})
    (output/"ratio_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    p=cfg.get("protocol",{}); manual=p.get("production_upper_wall")
    policy=f'''\n[resolved_wall_policy]\nmode = "{"manual" if manual is not None else p.get("production_upper_wall_mode","auto")}"\nminimum = {float(p.get("production_upper_wall_min",120.0))}\nclearance = {float(p.get("production_wall_clearance",30.0))}\nrounding = {float(p.get("production_wall_rounding",10.0))}\nbox_margin = {float(p.get("production_box_margin",5.0))}\nmanual_override = {"true" if manual is not None else "false"}\n'''
    if manual is not None: policy+=f"manual_value = {float(manual)}\n"
    policy+='''\n[resolved_wall_results.equilibrate_300K]\nstatus = "DEFERRED"\nsource = "held-300K.data"\n\n[resolved_wall_results.anneal_400K]\nstatus = "DEFERRED"\nsource = "equilibrated-300K.data"\n'''
    (output/"resolved_config.toml").write_text(config_path.read_text(encoding="utf-8")+policy,encoding="utf-8")
    (output/"type_map.json").write_text(json.dumps(type_map,indent=2)+"\n",encoding="utf-8")
    hashes={}
    for t in templates:
        for p in (t["folder"]/"molecule.toml",t["folder"]/t["manifest"]["files"]["ligpargen"]): hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    provenance_surface=primary_final if primary_final else surface_path
    hashes[str(provenance_surface.relative_to(ROOT)) if provenance_surface.is_relative_to(ROOT) else str(provenance_surface)]=hashlib.sha256(provenance_surface.read_bytes()).hexdigest()
    (output/"input_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n",encoding="utf-8")
    type_map.update({"status":"assembled","source_hashes":hashes,"output_hashes":{},"total_charge":charge(result),"molecule_count":max(int(r.fields[1]) for r in result.sections["Atoms"]),"counts":{sec.lower():result.count(sec) for sec in ("Atoms","Bonds","Angles","Dihedrals","Impropers")},"box":result.bounds})
    (output/"assembly_manifest.json").write_text(json.dumps(type_map,indent=2)+"\n",encoding="utf-8")
    from .validate import validate
    validate(
        output,
        packmol_ran=packed_ok,
        primary_final=primary_final,
        run_lammps=primary_final is None,
    )
    return output

def _write_input(output: Path,cfg:dict,data:DataFile)->None:
    p=cfg.get("protocol",{}); steps=int(p.get("deposition_steps",600000)); temp=float(p.get("temperature",300.0))
    deposition_timestep=float(p.get("deposition_timestep",p.get("deposition_slow_timestep",0.5))); regular_timestep=float(p.get("timestep",1.0))
    hold_steps=int(p.get("hold_300K_steps",1000000)); hold_timestep=float(p.get("hold_300K_timestep",deposition_timestep))
    hold_400_steps=int(p.get("hold_400K_steps",hold_steps)); hold_400_timestep=float(p.get("hold_400K_timestep",hold_timestep))
    heat_400_steps=int(p.get("heat_400K_steps",400000)); heat_400_timestep=float(p.get("heat_400K_timestep",hold_400_timestep))
    if steps <= 0: raise ValueError("deposition_steps must be positive")
    if deposition_timestep <= 0: raise ValueError("deposition_timestep must be positive")
    if hold_steps <= 0: raise ValueError("hold_300K_steps must be positive")
    if hold_timestep <= 0: raise ValueError("hold_300K_timestep must be positive")
    if hold_400_steps <= 0: raise ValueError("hold_400K_steps must be positive")
    if hold_400_timestep <= 0: raise ValueError("hold_400K_timestep must be positive")
    if heat_400_steps <= 0: raise ValueError("heat_400K_steps must be positive")
    if heat_400_timestep <= 0: raise ValueError("heat_400K_timestep must be positive")
    def lower_wall(value: object) -> str:
        if isinstance(value,str):
            if value.upper()!="EDGE": raise ValueError("lower wall coordinate must be numeric or EDGE")
            return "EDGE"
        return str(float(value))
    deposition_lower=lower_wall(p.get("deposition_lower_wall",p.get("lower_wall","EDGE")))
    lower_300=lower_wall(p.get("production_lower_wall_300","EDGE")); lower_400=lower_wall(p.get("production_lower_wall_400","EDGE"))
    manual_upper=p.get("production_upper_wall"); wall_min=float(p.get("production_upper_wall_min",120.0)); wall_clearance=float(p.get("production_wall_clearance",30.0)); wall_rounding=float(p.get("production_wall_rounding",10.0)); box_margin=float(p.get("production_box_margin",5.0))
    surface_top=max(float(a.fields[6]) for a in data.sections["Atoms"] if abs(abs(float(a.fields[3]))-2.0)<1e-8)
    configured_endpoint=p.get("deposition_wall_endpoint")
    zend=float(configured_endpoint) if configured_endpoint is not None else surface_top+float(p.get("wall_clearance",30.0))
    zstart=max(zend+100.0,data.bounds["z"][1]-5.0)
    primary_atom_count=cfg.get("_sequential_primary_atom_count")
    sequential=primary_atom_count is not None
    if sequential:
        primary_atom_count=int(primary_atom_count)
        if not 0 < primary_atom_count < data.count("Atoms"):
            raise ValueError("sequential primary atom boundary is invalid")
        stage2_groups=f"group stage2 id {primary_atom_count+1}:{data.count('Atoms')}\ngroup stage1 subtract all stage2\n"
        minimization_lock="fix stage1_lock stage1 setforce 0.0 0.0 0.0\n"
        minimization_unlock="unfix stage1_lock\n"
        velocity_initialization=f"velocity stage2 create {temp} 214587 mom yes rot yes dist gaussian"
        deposition_temperature_start=temp
        deposition_label="Sequential stage-2 deposition onto equilibrated Me-4PACz"
        deposition_thermal_note=f"Stage 1 is fixed only during FIRE minimization; only stage-2 atoms receive new {temp:g} K velocities, and the full system then deposits at constant {temp:g} K."
    else:
        stage2_groups=""
        minimization_lock=""
        minimization_unlock=""
        velocity_initialization="velocity all create 5.0 214587 mom yes rot yes dist gaussian"
        deposition_temperature_start=5.0
        deposition_label="Cao deposition"
        deposition_thermal_note=f"The full system starts at 5 K and heats continuously to {temp:g} K during deposition."
    init=lambda source: f"""boundary p p f
processors * * 1
units real
atom_style full
read_data {source}
include force_field_settings_lammps_with_header.lmp
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
timestep {regular_timestep}
thermo 1000
thermo_style custom step temp pe ke etotal press pxx pyy lx ly lz fnorm fmax
variable pressure equal 1.0
variable pressureDamp equal 500.0
variable epsilon equal 1.0
variable sigma equal 1.0
variable cutoff equal 2.5
"""
    def resolution(source:str)->str:
        if manual_upper is not None:
            frozen=f"variable resolved_wall_hi equal {float(manual_upper)}\n"
        else:
            frozen=f"""compute stage_zmax all reduce max z
run 0 post no
variable measured_max_z equal c_stage_zmax
variable wall_candidate equal 0.5*({wall_min}+v_measured_max_z+{wall_clearance}+abs({wall_min}-(v_measured_max_z+{wall_clearance})))
variable rounded_wall equal ceil(v_wall_candidate/{wall_rounding})*{wall_rounding}
variable resolved_wall_hi equal ${{rounded_wall}}
"""
        return frozen+f"""variable required_box_zhi equal v_resolved_wall_hi+{box_margin}
variable frozen_box_zhi equal ${{required_box_zhi}}
if "$(zhi) < ${{frozen_box_zhi}}" then "change_box all z final $(zlo) ${{frozen_box_zhi}} units box"
print "Resolved production wall: ${{resolved_wall_hi}}"
print "Final box zhi: $(zhi)"
"""
    text=f"""# {deposition_label}; safely rerun only when deposited.data is absent
{init('topology_output.lmp')}{stage2_groups}fix walllo all wall/lj93 zlo {deposition_lower} ${{epsilon}} ${{sigma}} ${{cutoff}} units box
fix_modify walllo energy yes
timestep {deposition_timestep}
{minimization_lock}min_style fire
min_modify dmax 0.05
minimize 0.0 1.0 20000 200000
{minimization_unlock}variable optimization_step equal step
variable optimization_pe equal pe
variable optimization_fnorm equal fnorm
variable optimization_fmax equal fmax
variable optimization_force_target_met equal "fnorm <= 1.0"
print "step=${{optimization_step}} potential_energy_kcal_per_mol=${{optimization_pe}} force_two_norm_kcal_per_mol_A=${{optimization_fnorm}} force_max_kcal_per_mol_A=${{optimization_fmax}} force_target_kcal_per_mol_A=1.0 force_target_met=${{optimization_force_target_met}}" file optimization-summary.txt screen yes
write_data optimized.data nocoeff
write_restart optimized.restart
write_dump all custom optimized.lammpstrj id mol type q x y z fx fy fz modify sort id
reset_timestep 0
{velocity_initialization}
variable zstart equal {zstart}
variable zend equal {zend}
variable zwall equal \"v_zstart - (v_zstart-v_zend)*(step/{steps}.0)\"
thermo_style custom step temp pe ke etotal press pxx pyy lx ly lz v_zwall fnorm fmax
variable epsilon equal 1.0
variable sigma equal 1.0
variable cutoff equal 2.5
fix wall all wall/lj126 zhi v_zwall ${{epsilon}} ${{sigma}} ${{cutoff}}
fix deposit all npt temp {deposition_temperature_start} {temp} 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy
dump trajectory all custom 1000 deposition.lammpstrj id mol type q x y z
dump_modify trajectory sort id
restart 100000 deposition.restart.1 deposition.restart.2
timestep {deposition_timestep}
run {steps} start 0 stop {steps}
unfix deposit
unfix wall
unfix walllo
undump trajectory
write_data deposited.data nocoeff
write_restart deposited.restart
"""
    (output/"deposition.in").write_text(text,encoding="utf-8")
    (output/"lammps.in").write_text(text,encoding="utf-8")
    hold=f"""# Compressed-film hold: run only after deposition.in
{init('deposited.data')}variable hold_wall_hi equal {zend}
thermo_style custom step temp pe ke etotal press pxx pyy lx ly lz v_hold_wall_hi fnorm fmax
fix lo all wall/lj93 zlo {lower_300} ${{epsilon}} ${{sigma}} ${{cutoff}} units box
fix hi all wall/lj126 zhi ${{hold_wall_hi}} ${{epsilon}} ${{sigma}} ${{cutoff}} units box
fix_modify lo energy yes
fix_modify hi energy yes
fix ensemble all npt temp {temp} {temp} 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy
dump trajectory all custom 2000 hold-300K.lammpstrj id mol type q x y z
dump_modify trajectory sort id
restart 200000 hold-300K.restart.1 hold-300K.restart.2
timestep {hold_timestep}
run {hold_steps}
write_data held-300K.data nocoeff
write_restart held-300K.restart
"""
    hold_400=f"""# Independent gradual heating and compressed-film 400 K hold: run only after deposition.in
{init('deposited.data')}variable hold_wall_hi equal {zend}
thermo_style custom step temp pe ke etotal press pxx pyy lx ly lz v_hold_wall_hi fnorm fmax
fix lo all wall/lj93 zlo {lower_400} ${{epsilon}} ${{sigma}} ${{cutoff}} units box
fix hi all wall/lj126 zhi ${{hold_wall_hi}} ${{epsilon}} ${{sigma}} ${{cutoff}} units box
fix_modify lo energy yes
fix_modify hi energy yes
fix heating all npt temp {temp} 400.0 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy
dump trajectory all custom 2000 heat-400K.lammpstrj id mol type q x y z
dump_modify trajectory sort id
restart 100000 heat-400K.restart.1 heat-400K.restart.2
timestep {heat_400_timestep}
run {heat_400_steps}
unfix heating
undump trajectory
write_data heated-400K.data nocoeff
write_restart heated-400K.restart
fix hold all npt temp 400.0 400.0 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy
dump trajectory all custom 2000 hold-400K.lammpstrj id mol type q x y z
dump_modify trajectory sort id
restart 200000 hold-400K.restart.1 hold-400K.restart.2
timestep {hold_400_timestep}
run {hold_400_steps}
write_data held-400K.data nocoeff
write_restart held-400K.restart
"""
    eq=f"# Cao 300 K continuation: run only after hold-300K.in\n{init('held-300K.data')}{resolution('held-300K.data')}fix lo all wall/lj93 zlo {lower_300} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix hi all wall/lj93 zhi ${{resolved_wall_hi}} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix ensemble all npt temp {temp} {temp} 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy\ndump trajectory all custom 10000 equilibration-300K.lammpstrj id mol type q x y z\ndump_modify trajectory sort id\nrestart 500000 equilibration-300K.restart.1 equilibration-300K.restart.2\nrun 5000000\nwrite_data equilibrated-300K.data nocoeff\nwrite_restart equilibrated-300K.restart\n"
    anneal=f"# Cao 400 K continuation: run only after equilibrate-300K.in\n{init('equilibrated-300K.data')}{resolution('equilibrated-300K.data')}fix lo all wall/lj93 zlo {lower_400} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix hi all wall/lj93 zhi ${{resolved_wall_hi}} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix ensemble all npt temp 400.0 400.0 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy\ndump trajectory all custom 10000 anneal-400K.lammpstrj id mol type q x y z\ndump_modify trajectory sort id\nrestart 500000 anneal-400K.restart.1 anneal-400K.restart.2\nrun 3000000\nwrite_data annealed-400K.data nocoeff\nwrite_restart annealed-400K.restart\n"
    (output/"hold-300K.in").write_text(hold,encoding="utf-8")
    (output/"hold-400K.in").write_text(hold_400,encoding="utf-8")
    (output/"equilibrate-300K.in").write_text(eq,encoding="utf-8")
    (output/"anneal-400K.in").write_text(anneal,encoding="utf-8")
    for name,source in (("deposition",text),("hold-300K",hold),("hold-400K",hold_400),("equilibrate-300K",eq),("anneal-400K",anneal)):
        check=re.sub(r"(?m)^run\s+.*$","run 0",source)
        check=re.sub(r"(?m)^minimize\s+.*$","run 0 post no # minimization replaced by force evaluation for smoke validation",check)
        (output/f"validate-{name}.in").write_text("# Temporary-directory smoke form of the real stage\n"+check,encoding="utf-8")
    mode=f"manual override {float(manual_upper):g} A" if manual_upper is not None else f"automatic: ceil(max({wall_min:g}, max_atom_z+{wall_clearance:g})/{wall_rounding:g})*{wall_rounding:g} A"
    duration_ps=steps*deposition_timestep/1000.0; hold_duration_ps=hold_steps*hold_timestep/1000.0; heat_400_duration_ps=heat_400_steps*heat_400_timestep/1000.0; hold_400_duration_ps=hold_400_steps*hold_400_timestep/1000.0
    (output/"protocol_notes.txt").write_text(f"Deposition endpoint: {zend:.3f} A. Deposition uses {steps} steps at a continuous {deposition_timestep:g} fs timestep ({duration_ps:g} ps total), with a continuous upper-wall ramp and no timestep transition. {deposition_thermal_note} The minimized structure and force summary are written before deposition dynamics. Independent compressed-film temperature branches both start from deposited.data and keep zhi={zend:.3f} A. The 300 K branch holds for {hold_steps} steps at {hold_timestep:g} fs ({hold_duration_ps:g} ps). The 400 K branch first heats from {temp:g} to 400 K for {heat_400_steps} steps at {heat_400_timestep:g} fs ({heat_400_duration_ps:g} ps), writes heated-400K.data, then holds at 400 K for {hold_400_steps} steps at {hold_400_timestep:g} fs ({hold_400_duration_ps:g} ps). Production wall policy after the 300 K hold: {mode}; box margin {box_margin:g} A. Long 300 K source: held-300K.data (DEFERRED until present), zlo={lower_300}. 400 K source: equilibrated-300K.data (DEFERRED until present), zlo={lower_400}. Each resolved production wall is frozen before dynamics and zhi alone is expanded if required.\n",encoding="utf-8")

def refresh_inputs(config_path: Path, output: Path) -> Path:
    """Regenerate stage inputs without rebuilding or modifying simulation data."""
    topology=output/"topology_output.lmp"
    if not topology.is_file():
        raise FileNotFoundError(f"existing topology is missing: {topology}")
    cfg=load(config_path)
    manifest_path=output/"assembly_manifest.json"
    if manifest_path.is_file():
        manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
        primary=next((component for component in manifest.get("components",[]) if component.get("component")=="stage1_primary"),None)
        if primary: cfg["_sequential_primary_atom_count"]=int(primary["atom_ids"][1])
    _write_input(output,cfg,parse(topology))
    return output
