from __future__ import annotations
from pathlib import Path
import hashlib, json, math, re, shutil, subprocess, tempfile, tomllib
from .lammps import TOPOLOGY, charge, parse

def _records(data, section):
    return [r.fields for r in data.sections.get(section,[])]

def _minimum_cross_separation(old_atoms,new_atoms,cutoff=2.0):
    cells={}
    for row in old_atoms:
        xyz=tuple(float(x) for x in row.fields[4:7]); key=tuple(math.floor(x/cutoff) for x in xyz); cells.setdefault(key,[]).append(xyz)
    minimum=float("inf")
    for row in new_atoms:
        xyz=tuple(float(x) for x in row.fields[4:7]); key=tuple(math.floor(x/cutoff) for x in xyz)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for other in cells.get((key[0]+dx,key[1]+dy,key[2]+dz),[]):
                        minimum=min(minimum,math.dist(xyz,other))
    # No candidates in neighboring cutoff cells proves every cross distance is
    # at least cutoff, without an O(N*M) all-pairs calculation.
    return cutoff if math.isinf(minimum) else minimum

def _wall_result(source: Path, policy: dict) -> dict:
    data=parse(source); max_z=max(float(a.fields[6]) for a in data.sections["Atoms"])
    manual=bool(policy.get("manual_override",False))
    if manual: wall=float(policy["manual_value"])
    else:
        candidate=max(float(policy["minimum"]),max_z+float(policy["clearance"]))
        quantum=float(policy["rounding"]); wall=math.ceil((candidate-1e-12)/quantum)*quantum
    final_zhi=max(data.bounds["z"][1],wall+float(policy["box_margin"]))
    return {"status":"RESOLVED","source":source.name,"measured_max_atom_z":max_z,"minimum_wall":float(policy["minimum"]),"clearance":float(policy["clearance"]),"rounding":float(policy["rounding"]),"box_margin":float(policy["box_margin"]),"resolved_wall":wall,"final_zlo":data.bounds["z"][0],"final_zhi":final_zhi,"mode":"manual" if manual else "auto"}

def _record_wall_results(folder: Path, results: dict[str,dict]) -> None:
    path=folder/"resolved_config.toml"; text=path.read_text()
    for stage,result in results.items():
        body="\n".join(f'{key} = "{value}"' if isinstance(value,str) else f"{key} = {str(value).lower() if isinstance(value,bool) else value}" for key,value in result.items())
        pattern=rf"(?ms)^\[resolved_wall_results\.{stage}\]\n.*?(?=^\[|\Z)"
        text=re.sub(pattern,f"[resolved_wall_results.{stage}]\n{body}\n\n",text)
    path.write_text(text,encoding="utf-8")
    base=(folder/"protocol_notes.txt").read_text().split("\nWALL RESOLUTION RESULTS\n",1)[0]
    lines=[base.rstrip(),"","WALL RESOLUTION RESULTS"]
    for stage,result in results.items(): lines.append(f"{stage}: "+", ".join(f"{k}={v}" for k,v in result.items()))
    (folder/"protocol_notes.txt").write_text("\n".join(lines)+"\n",encoding="utf-8")

def _refresh_output_hashes(folder: Path) -> None:
    path=folder/"assembly_manifest.json"; manifest=json.loads(path.read_text())
    manifest["output_hashes"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir()) if p.is_file() and p.name!="assembly_manifest.json"}
    path.write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")

def validate(folder: Path, packmol_ran: bool|None=None, primary_final: Path|None=None) -> bool:
    errors=[]; warnings=[]; deferred=[]; wall_results={}
    topology=folder/"topology_output.lmp"
    if not topology.exists():
        raise ValueError("packing is incomplete; topology_output.lmp was not produced. Run: packmol < packmol.inp")
    data=parse(topology); manifest=json.loads((folder/"assembly_manifest.json").read_text())
    atoms=data.sections["Atoms"]; ids=[int(r.fields[0]) for r in atoms]; by_id={int(r.fields[0]):r for r in atoms}
    for sec in ("Atoms","Bonds","Angles","Dihedrals","Impropers"):
        declared=data.declared_counts.get(sec)
        if declared is not None and declared!=data.count(sec): errors.append(f"header declares {declared} {sec}, found {data.count(sec)}")
    if len(ids)!=len(set(ids)): errors.append("atom IDs are not unique")
    atom_set=set(ids)
    positive_molecules={int(r.fields[1]) for r in atoms if int(r.fields[1])>0}
    expected=sum(c["count"] for c in manifest.get("components",[]))
    if len(positive_molecules)<expected: errors.append("physical ligand copies do not have unique positive molecule IDs")
    for sec,nrefs in TOPOLOGY.items():
        rows=data.sections.get(sec,[]); rowids=[int(r.fields[0]) for r in rows]
        if len(rowids)!=len(set(rowids)): errors.append(f"duplicate {sec} IDs")
        for r in rows:
            refs=[int(x) for x in r.fields[2:2+nrefs]]
            if any(x not in atom_set for x in refs): errors.append(f"invalid atom reference in {sec}"); continue
            memberships={int(by_id[x].fields[1]) for x in refs}
            if len(memberships)!=1: errors.append(f"{sec} term connects different molecules")
            if 0 in memberships: errors.append(f"surface atom appears in bonded {sec} term")
    mass_types={int(r.fields[0]) for r in data.sections["Masses"]}; used={int(r.fields[2]) for r in atoms}
    if used-mass_types: errors.append(f"used atom types lack masses: {sorted(used-mass_types)}")
    coefficient_commands={"Bonds":"bond_coeff","Angles":"angle_coeff","Dihedrals":"dihedral_coeff","Impropers":"improper_coeff"}
    ff=(folder/"force_field_settings_lammps_with_header.lmp").read_text()
    for topology_section,command in coefficient_commands.items():
        used_bonded={int(r.fields[1]) for r in data.sections.get(topology_section,[])}
        defined={int(match.group(1)) for match in re.finditer(rf"^{command}\s+(\d+)\s+",ff,re.M)}
        if used_bonded-defined: errors.append(f"used {topology_section} types lack coefficients: {sorted(used_bonded-defined)}")
    for r in atoms:
        values=[float(x) for x in r.fields[3:7]]
        if not all(math.isfinite(x) for x in values): errors.append(f"non-finite atom value at ID {r.fields[0]}")
    for axis,i in zip("xyz",range(3)):
        lo,hi=data.bounds[axis]
        if any(not(lo<=float(r.fields[4+i])<=hi) for r in atoms): errors.append(f"coordinates outside {axis} bounds")
    nio=[r for r in atoms if abs(abs(float(r.fields[3]))-2.0)<1e-8]
    if not nio or abs(sum(float(r.fields[3]) for r in nio))>1e-8: errors.append("Ni/O fixed-charge population missing or unbalanced")
    if len(nio)!=21060: errors.append(f"surface atom count is {len(nio)}; expected 21060")
    surface_path=Path(__file__).resolve().parents[2]/"inputs/surfaces/corrugated-nio-110/surface.lmp"
    if hashlib.sha256(surface_path.read_bytes()).hexdigest()!="09be07e828f487ab3599f2cf533c969ce6fad31391044055f4572bf81caff556": errors.append("authoritative surface SHA256 mismatch")
    nio_types=sorted({int(r.fields[2]) for r in nio})
    ligand_types=sorted(used-set(nio_types))
    if len(nio_types)!=2: errors.append("expected exactly two Ni/O atom types")
    else:
        ni,o=nio_types
        for a,b in ((ni,ni),(ni,o),(o,o)):
            if not re.search(rf"^pair_coeff\s+{a}\s+{b}\s+buck/coul/long",ff,re.M): errors.append(f"missing Buckingham surface pair {a}-{b}")
        for surface in nio_types:
            for ligand in ligand_types:
                if not re.search(rf"^pair_coeff\s+{surface}\s+{ligand}\s+lj/cut/coul/long",ff,re.M): errors.append(f"missing LJ ligand-surface pair {surface}-{ligand}")
        if any(re.search(rf"^pair_coeff\s+{a}\s+{b}\s+buck/coul/long",ff,re.M) for a in ligand_types for b in ligand_types): errors.append("ligand-ligand pair accidentally uses Buckingham")
    mapped_types=set()
    for component in manifest.get("components",[]):
        for bounds in component.get("types",{}).values(): mapped_types.update(range(bounds[0],bounds[1]+1))
    if set(ligand_types)-mapped_types: errors.append("topology atom types disagree with assembly mapping")
    if abs(charge(data)-float(manifest["total_charge"]))>1e-6: errors.append("final charge disagrees with component/assembly manifest")
    if primary_final:
        old=parse(primary_final)
        if old.bounds!=data.bounds: errors.append("sequential stage changed stage-1 box bounds")
        for sec in ("Masses","Atoms","Bonds","Angles","Dihedrals","Impropers"):
            before=_records(old,sec); after=_records(data,sec)[:len(before)]
            if before!=after: errors.append(f"sequential stage-1 {sec} changed")
        separation=_minimum_cross_separation(old.sections["Atoms"],data.sections["Atoms"][old.count("Atoms"):])
        if separation<2.0: errors.append(f"stage-1/stage-2 minimum separation {separation:.6f} is below 2.0 A")
        manifest["stage1_stage2_minimum_separation_lower_bound_angstrom"]=separation
        (folder/"assembly_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    required_order=("boundary p p f","units real","atom_style full","read_data ","include ")
    for input_name in ("deposition.in","hold-300K.in","hold-400K.in","equilibrate-300K.in","anneal-400K.in"):
        text=(folder/input_name).read_text(); positions=[text.find(token) for token in required_order]
        if any(x<0 for x in positions) or positions!=sorted(positions): errors.append(f"{input_name}: invalid initialization command order")
        for directive in ("thermo_style","dump trajectory","restart ","write_data","write_restart"):
            if directive not in text: errors.append(f"{input_name}: missing {directive.strip()} directive")
        for line in re.findall(r"(?m)^write_data\s+.*$",text):
            if not line.endswith(" nocoeff"): errors.append(f"{input_name}: write_data is not coefficient-free")
    deposition=(folder/"deposition.in").read_text(); hold=(folder/"hold-300K.in").read_text(); hold_400=(folder/"hold-400K.in").read_text(); eq=(folder/"equilibrate-300K.in").read_text(); anneal=(folder/"anneal-400K.in").read_text()
    if not re.search(r"variable zend equal 69\.615(?:0+)?",deposition): errors.append("deposition wall endpoint is not 69.615 A")
    if "read_data deposited.data" not in hold or "read_data held-300K.data" not in eq or "read_data equilibrated-300K.data" not in anneal: errors.append("continuation stage-local source selection is invalid")
    if not re.search(r"variable hold_wall_hi equal 69\.615(?:0+)?",hold): errors.append("300 K hold wall is not fixed at the deposition endpoint")
    if "read_data deposited.data" not in hold_400: errors.append("400 K hold does not branch from deposited.data")
    if not re.search(r"variable hold_wall_hi equal 69\.615(?:0+)?",hold_400): errors.append("400 K hold wall is not fixed at the deposition endpoint")
    if "fix ensemble all npt temp 400.0 400.0" not in hold_400: errors.append("400 K hold thermostat is invalid")
    resolved=tomllib.loads((folder/"resolved_config.toml").read_text()); policy=resolved["resolved_wall_policy"]
    for stage,source_name in (("equilibrate_300K","held-300K.data"),("anneal_400K","equilibrated-300K.data")):
        source=folder/source_name
        if not source.exists():
            result={"status":"DEFERRED","source":source_name,"measured_max_atom_z":"DEFERRED","minimum_wall":float(policy["minimum"]),"clearance":float(policy["clearance"]),"rounding":float(policy["rounding"]),"box_margin":float(policy["box_margin"]),"resolved_wall":"DEFERRED","final_zlo":"DEFERRED","final_zhi":"DEFERRED","mode":"manual" if policy.get("manual_override") else "auto"}; deferred.append(f"{stage}: waiting for real {source_name}")
        else:
            result=_wall_result(source,policy)
            if not result["measured_max_atom_z"] < result["resolved_wall"]: errors.append(f"{stage}: atom is not strictly below resolved wall")
            if result["final_zhi"] < result["resolved_wall"]+float(policy["box_margin"]): errors.append(f"{stage}: resolved wall lacks box margin")
            if not result["measured_max_atom_z"] < result["final_zhi"]: errors.append(f"{stage}: atom is not strictly below final box zhi")
        wall_results[stage]=result
    _record_wall_results(folder,wall_results)
    if packmol_ran is False: errors.append("Packmol packing was not completed")
    executable=shutil.which("lmp") or shutil.which("lammps")
    if not executable: warnings.append("LAMMPS executable not found; zero-step checks skipped")
    else:
        with tempfile.TemporaryDirectory(prefix="nio-md-validate-") as temp:
            target=Path(temp); shutil.copy2(topology,target/topology.name)
            shutil.copy2(folder/"force_field_settings_lammps_with_header.lmp",target/"force_field_settings_lammps_with_header.lmp")
            smoke_stages=[("deposition",None),("hold-300K","deposited.data"),("hold-400K","deposited.data"),("equilibrate-300K","held-300K.data"),("anneal-400K","equilibrated-300K.data")]
            for name,predecessor in smoke_stages:
                if predecessor and not (folder/predecessor).exists(): continue
                if predecessor: shutil.copy2(folder/predecessor,target/predecessor)
                source=folder/f"validate-{name}.in"; shutil.copy2(source,target/source.name)
                run=subprocess.run([executable,"-in",source.name],cwd=target,capture_output=True,text=True)
                if run.returncode: errors.append(f"LAMMPS zero-step {name} failed: {(run.stderr or run.stdout)[-500:]}")
                else: warnings.append(f"LAMMPS zero-step {name}: passed")
    status="INVALID" if errors else ("VALID_WITH_DEFERRED_STAGES" if deferred else "VALID")
    report=[status,f"atoms: {len(atoms)}",f"molecules: {max(int(r.fields[1]) for r in atoms)}",f"total charge: {charge(data):.8f}"]
    if not errors: report.append("Spreadsheet assembly: not required; topology and coordinates were merged programmatically.")
    for stage,result in wall_results.items(): report.append(f"WALL {stage}: "+", ".join(f"{k}={v}" for k,v in result.items()))
    report += [f"DEFERRED: {x}" for x in deferred]+[f"ERROR: {x}" for x in errors]+[f"WARNING: {x}" for x in warnings]
    (folder/"validation_report.txt").write_text("\n".join(report)+"\n",encoding="utf-8")
    _refresh_output_hashes(folder)
    if errors: raise ValueError("; ".join(dict.fromkeys(errors)))
    return True
