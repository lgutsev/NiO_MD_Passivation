from __future__ import annotations
from pathlib import Path
import math
from copy import deepcopy
import hashlib, json, shutil, subprocess
from .config import ROOT, load, molecule_manifest, missing_ligpargen
from .geometry import read_xyz, write_xyz
from .lammps import DataFile, Record, TOPOLOGY, charge, parse, replicate, write
from .chemistry import molecular_weight, correction_lines

HEADER = """# Validated NiO/passivant force-field styles
bond_style      harmonic
angle_style     harmonic
dihedral_style  hybrid opls
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

def _packmol(config: dict, templates: list[dict], output: Path, supplied: Path|None=None) -> tuple[Path,list[tuple[float,float,float]]|None,bool]:
    lines=["tolerance 2.0", "filetype xyz", "output packed.xyz", "seed 202405367", ""]
    for item in templates:
        lines += [f"structure {item['xyz'].name}", f"  number {item['count']}", f"  inside box {item['region']}", "end structure", ""]
        destination=output/item["xyz"].name
        if item["xyz"].resolve()!=destination.resolve(): shutil.copy2(item["xyz"],destination)
    inp=output/"packmol.inp"; inp.write_text("\n".join(lines),encoding="utf-8")
    exe=shutil.which("packmol") if supplied is None else None
    if supplied is not None:
        shutil.copy2(supplied,output/"packed.xyz")
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
        templates.append({"slug":spec["slug"],"data":data,"xyz":xyz,"count":int(spec["count"]),"atoms":data.count("Atoms"),"mw":molecular_weight(data),"region":spec.get("region","2.0 2.0 45.0 123.1 39.7 145.0"),"folder":folder,"manifest":manifest})
        manifests.append(manifest)
    packed,coords,packed_ok=_packmol(cfg,templates,output,packed_xyz)
    if coords is None:
        plan={"status":"packing_incomplete","packmol_seed":202405367,"command":"packmol < packmol.inp","component_order":[t["slug"] for t in templates]+([] if primary_final else ["corrugated-nio-110"]),"molecules":[{"slug":t["slug"],"count":t["count"],"atoms_per_molecule":t["atoms"]} for t in templates]}
        (output/"assembly_manifest.json").write_text(json.dumps(plan,indent=2)+"\n",encoding="utf-8")
        (output/"validation_report.txt").write_text("INCOMPLETE\nPackmol has not run; no scientifically usable topology was produced.\nRun: packmol < packmol.inp\n",encoding="utf-8")
        return output
    if primary_final:
        if not primary_final.exists(): raise FileNotFoundError(primary_final)
        result=parse(primary_final)
        ids={"atom":max(int(r.fields[0]) for r in result.sections["Atoms"])}
        for sec in TOPOLOGY: ids[sec.lower()[:-1]]=max((int(r.fields[0]) for r in result.sections.get(sec,[])),default=0)
        types={k:result.type_count(k) for k in ids}; mol=max((int(r.fields[1]) for r in result.sections["Atoms"]),default=0)
    else:
        result=DataFile("LAMMPS data file generated by nio-md-prep",{}, {})
        ids={"atom":0,"bond":0,"angle":0,"dihedral":0,"improper":0}; types={k:0 for k in ids}; mol=0
    cursor=0; type_map={"components":[],"surface":{},"packmol_seed":202405367}
    if primary_final:
        existing_atoms=result.sections["Atoms"]
        type_map["components"].append({"component":"stage1_primary","count":len({int(r.fields[1]) for r in existing_atoms if int(r.fields[1])>0}),"atoms_per_molecule":None,"charge":charge(result),"atom_ids":[min(int(r.fields[0]) for r in existing_atoms),max(int(r.fields[0]) for r in existing_atoms)],"molecule_ids":[0,max(int(r.fields[1]) for r in existing_atoms)],"types":{cat:[1,result.type_count(cat)] for cat in types if result.type_count(cat)}})
    overrides=[]
    for t in templates:
        n=t["atoms"]*t["count"]; local=coords[cursor:cursor+n]; cursor+=n
        pieces,inc=replicate(t["data"],t["count"],types,ids,mol,local); _merge(result,pieces)
        maps={cat:{i:i+types[cat] for i in range(1,t["data"].type_count(cat)+1)} for cat in types}
        ranges={"component":t["slug"],"count":t["count"],"atoms_per_molecule":t["atoms"],"charge":charge(t["data"])*t["count"],"atom_ids":[ids["atom"]+1,ids["atom"]+inc["atom"]],"molecule_ids":[mol+1,mol+inc["molecule"]],"types":{cat:[types[cat]+1,types[cat]+t["data"].type_count(cat)] for cat in types if t["data"].type_count(cat)}}
        type_map["components"].append(ranges); type_map.setdefault("type_maps",{})[t["slug"]]=maps
        chemical_lines, anchors=correction_lines(t["data"],types["atom"])
        expected_anchors=int(t["manifest"].get("molecule",{}).get("phosphonic_acid_anchors",0))
        if anchors != expected_anchors: raise ValueError(f"{t['slug']}: found {anchors} phosphonate anchors; manifest expects {expected_anchors}")
        overrides += chemical_lines
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
    result.bounds={axis:(min(p[i] for p in all_xyz)-margin,max(p[i] for p in all_xyz)+margin) for i,axis in enumerate("xyz")}
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
    report={"basis":ratio.get("basis","explicit"),"stock_concentration_mg_ml":ratio.get("stock_concentration_mg_ml",{}),"stock_volume_assumption":ratio.get("stock_volume_assumption","n/a"),"rounding_rule":ratio.get("rounding_rule","explicit counts"),"molecules":[{"slug":t["slug"],"molecular_weight_g_mol":round(t["mw"],6),"count":t["count"]} for t in templates]}
    if len(templates)==2:
        report["exact_secondary_to_primary_ratio"]=(templates[1]["count"]/templates[0]["count"])
    (output/"ratio_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    (output/"resolved_config.toml").write_text(config_path.read_text(encoding="utf-8"),encoding="utf-8")
    (output/"type_map.json").write_text(json.dumps(type_map,indent=2)+"\n",encoding="utf-8")
    hashes={}
    for t in templates:
        for p in (t["folder"]/"molecule.toml",t["folder"]/t["manifest"]["files"]["ligpargen"]): hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    provenance_surface=primary_final if primary_final else surface_path
    hashes[str(provenance_surface.relative_to(ROOT)) if provenance_surface.is_relative_to(ROOT) else str(provenance_surface)]=hashlib.sha256(provenance_surface.read_bytes()).hexdigest()
    (output/"input_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n",encoding="utf-8")
    outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file() and p.name!="assembly_manifest.json"}
    type_map.update({"status":"assembled","source_hashes":hashes,"output_hashes":outputs,"total_charge":charge(result),"molecule_count":max(int(r.fields[1]) for r in result.sections["Atoms"]),"counts":{sec.lower():result.count(sec) for sec in ("Atoms","Bonds","Angles","Dihedrals","Impropers")},"box":result.bounds})
    (output/"assembly_manifest.json").write_text(json.dumps(type_map,indent=2)+"\n",encoding="utf-8")
    from .validate import validate
    validate(output, packmol_ran=packed_ok, primary_final=primary_final)
    return output

def _write_input(output: Path,cfg:dict,data:DataFile)->None:
    p=cfg.get("protocol",{}); steps=int(p.get("deposition_steps",300000)); temp=float(p.get("temperature",300.0)); lower=float(p.get("lower_wall",-5))
    surface_top=max(float(a.fields[6]) for a in data.sections["Atoms"] if abs(abs(float(a.fields[3]))-2.0)<1e-8)
    zend=surface_top+float(p.get("wall_clearance",30.0)); zstart=max(zend+100.0,data.bounds["z"][1]-5.0)
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
write_data deposited.data
"""
    (output/"deposition.in").write_text(text,encoding="utf-8")
    (output/"lammps.in").write_text(text,encoding="utf-8")
    common="""boundary p p f\nunits real\natom_style full\ninclude force_field_settings_lammps_with_header.lmp\ntimestep 1.0\nvariable pressure equal 1.0\nvariable pressureDamp equal 500.0\nvariable epsilon equal 1.0\nvariable sigma equal 1.0\nvariable cutoff equal 2.5\n"""
    eq=f"# Cao 300 K equilibration: 5 ns\nread_data deposited.data\n{common}fix lo all wall/lj93 zlo {lower} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix hi all wall/lj93 zhi {zend} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix ensemble all npt temp {temp} {temp} 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy\nrun 5000000\nwrite_data equilibrated-300K.data\nwrite_restart equilibrated-300K.restart\n"
    anneal=f"# Cao 400 K anneal: 3 ns\nread_data equilibrated-300K.data\n{common}fix lo all wall/lj93 zlo {lower} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix hi all wall/lj93 zhi {zend} ${{epsilon}} ${{sigma}} ${{cutoff}} units box\nfix ensemble all npt temp 400.0 400.0 100.0 x ${{pressure}} ${{pressure}} ${{pressureDamp}} y ${{pressure}} ${{pressure}} ${{pressureDamp}} couple xy\nrun 3000000\nwrite_data annealed-400K.data\nwrite_restart annealed-400K.restart\n"
    (output/"equilibrate-300K.in").write_text(eq,encoding="utf-8")
    (output/"anneal-400K.in").write_text(anneal,encoding="utf-8")
    (output/"protocol_notes.txt").write_text("Cao2025 target: p p f, real/full, hybrid LJ/Buckingham, PPPM 1e-4, 1 fs, temperature damping 100 fs, moving wall, 5 ns at 300 K and 3 ns at 400 K. Established corrected executable decks use pressure damping 500 fs (retained here), despite the 1000 fs prose value. Legacy decks also contain 10 ns single-temperature production runs; the paper-target 5 ns + 3 ns split is generated here.\n",encoding="utf-8")
