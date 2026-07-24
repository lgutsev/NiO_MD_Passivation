from pathlib import Path
import hashlib, json, shutil, subprocess
import pytest
from nio_md_prep.build import build, refresh_inputs
from nio_md_prep.geometry import elements
from nio_md_prep.lammps import Record, parse, write
from nio_md_prep.validate import _wall_result

ROOT=Path(__file__).parents[1]

def packed_fixture(tmp_path, shift=(15.0,10.0,62.0), swap=False):
    tmp_path.mkdir(parents=True,exist_ok=True)
    template=parse(ROOT/"inputs/molecules/me-4pacz/ligpargen.lmp")
    symbols=elements(template)
    if swap:
        i=next(i for i,x in enumerate(symbols) if x=="C"); j=next(i for i,x in enumerate(symbols) if x=="H"); symbols[i],symbols[j]=symbols[j],symbols[i]
    shifts=(shift,shift,shift) if isinstance(shift,(int,float)) else shift
    rows=[]
    for symbol,atom in zip(symbols,template.sections["Atoms"]): rows.append(f"{symbol} {float(atom.fields[4])+shifts[0]} {float(atom.fields[5])+shifts[1]} {float(atom.fields[6])+shifts[2]}")
    path=tmp_path/"packed.xyz"; path.write_text(f"{len(rows)}\nmock Packmol\n"+"\n".join(rows)+"\n")
    return path,template

def test_small_offline_build(tmp_path):
    packed,template=packed_fixture(tmp_path)
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"build",packed_xyz=packed)
    data=parse(out/"topology_output.lmp")
    assert data.count("Atoms")==21060+45
    assert (out/"type_map.json").exists()
    assert (out/"validation_report.txt").read_text().startswith("VALID_WITH_DEFERRED_STAGES")
    force_field=(out/"force_field_settings_lammps_with_header.lmp").read_text()
    assert "special_bonds   amber" in force_field
    assert "dihedral_style  hybrid opls charmm" in force_field
    for typ in (12,14,15,17,56):
        assert f"dihedral_coeff {typ} charmm 0.5333 0 3 0.000" in force_field
    first_pair=template.sections["Pair Coeffs"][0].fields
    assert f"pair_coeff 1 1 lj/cut/coul/long {first_pair[1]} {first_pair[2]}" in force_field
    first=data.sections["Atoms"][0]
    assert float(first.fields[4])==float(template.sections["Atoms"][0].fields[4])+15
    assert first.fields[3]==template.sections["Atoms"][0].fields[3]
    assert "Spreadsheet assembly: not required" in (out/"validation_report.txt").read_text()

def test_packmol_element_order_failure(tmp_path):
    packed,_=packed_fixture(tmp_path,swap=True)
    with pytest.raises(ValueError,match="Packmol atom order mismatch"):
        build(ROOT/"tests/data/small-study.toml",tmp_path/"bad",packed_xyz=packed)

def test_small_packmol_boundary_overshoot_is_corrected_rigidly(tmp_path):
    packed,template=packed_fixture(tmp_path/"fixture",shift=(18.59801,10.0,62.0))
    before=[tuple(map(float,row.split()[1:4])) for row in packed.read_text().splitlines()[2:]]
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"build",packed_xyz=packed)
    after=parse(out/"topology_output.lmp").sections["Atoms"][:template.count("Atoms")]
    after_xyz=[tuple(map(float,row.fields[4:7])) for row in after]
    shifts={tuple(round(b-a,8) for a,b in zip(old,new)) for old,new in zip(before,after_xyz)}
    assert len(shifts)==1
    assert max(x for x,_,_ in after_xyz)<20.0
    report=json.loads((out/"packmol_adjustments.json").read_text())
    assert report["method"]=="rigid_whole_molecule_translation"
    assert len(report["adjustments"])==1

def test_missing_packmol_is_incomplete(tmp_path,monkeypatch):
    monkeypatch.setattr("nio_md_prep.build.shutil.which",lambda _:None)
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"incomplete")
    assert not (out/"topology_output.lmp").exists()
    assert "packmol < packmol.inp" in (out/"validation_report.txt").read_text()

def test_sequential_stage_preserves_existing_records(tmp_path,monkeypatch):
    packed,template=packed_fixture(tmp_path/"one")
    stage1=build(ROOT/"tests/data/small-study.toml",tmp_path/"stage1",packed_xyz=packed)
    before=parse(stage1/"topology_output.lmp")
    before.sections["Velocities"]=[Record([atom.fields[0],"1.0","2.0","3.0"]) for atom in before.sections["Atoms"]]
    for atom in before.sections["Atoms"]:
        if int(atom.fields[1])==0:
            atom.fields[6]=f"{float(atom.fields[6])+0.25:.8f}"
        atom.fields.extend(("0","0","0"))
    primary=stage1/"held-300K.data"; write(before,primary)
    primary.write_text(
        primary.read_text().replace(
            f"{before.bounds['z'][0]:.8f} {before.bounds['z'][1]:.8f} zlo zhi",
            f"{before.bounds['z'][0]:.12f} {before.bounds['z'][1]+1e-10:.12f} zlo zhi",
        )
    )
    before=parse(primary)
    packed2,_=packed_fixture(tmp_path/"two",shift=(20.0,20.0,100.0))
    monkeypatch.setattr("nio_md_prep.validate.shutil.which",lambda _: "/unintended/user-local/lmp")
    stage2=build(ROOT/"tests/data/small-study.toml",tmp_path/"stage2",primary_final=primary,packed_xyz=packed2)
    after=parse(stage2/"topology_output.lmp")
    for section in ("Atoms","Bonds","Angles","Dihedrals","Impropers"):
        assert [r.fields for r in after.sections[section]][:before.count(section)]==[r.fields for r in before.sections[section]]
    assert after.count("Atoms")==before.count("Atoms")+45
    assert {len(atom.fields) for atom in after.sections["Atoms"]}=={10}
    assert all(atom.fields[-3:]==["0","0","0"] for atom in after.sections["Atoms"][before.count("Atoms"):])
    assert after.count("Velocities")==after.count("Atoms")
    assert [r.fields for r in after.sections["Velocities"]][:before.count("Atoms")]==[r.fields for r in before.sections["Velocities"]]
    assert all(r.fields[1:]==["0.0","0.0","0.0"] for r in after.sections["Velocities"][before.count("Atoms"):])
    for axis in before.bounds:
        assert after.bounds[axis]==pytest.approx(before.bounds[axis],abs=1e-7)
    assert "inside box 2.000000 2.000000 67.509910 123.100000 39.700000 265.000000" in (stage2/"packmol.inp").read_text()
    deposition=(stage2/"deposition.in").read_text()
    assert f"group stage2 id {before.count('Atoms')+1}:{after.count('Atoms')}" in deposition
    assert "group stage1 subtract all stage2" in deposition
    assert deposition.index("fix stage1_lock stage1 setforce 0.0 0.0 0.0") < deposition.index("minimize ") < deposition.index("unfix stage1_lock")
    assert "min_style fire" not in deposition
    assert deposition.index("min_style sd") < deposition.index("minimize 0.0 10.0 5000 50000")
    assert deposition.index("minimize 0.0 10.0 5000 50000") < deposition.index("min_style cg")
    assert deposition.index("min_style cg") < deposition.index("minimize 0.0 1.0 20000 200000")
    assert "min_modify dmax 0.005 line backtrack" in deposition
    assert "min_modify dmax 0.01 line backtrack" in deposition
    assert deposition.index("unfix stage1_lock") < deposition.index("velocity stage2 create 300.0")
    assert "fix deposit all npt temp 300.0 300.0" in deposition
    assert "variable zend equal 69.615" in deposition
    assert "velocity all create" not in deposition
    assert "LAMMPS zero-step checks skipped during build-only sequential preparation" in (stage2/"validation_report.txt").read_text()
    refresh_inputs(ROOT/"tests/data/small-study.toml",stage2)
    assert "velocity stage2 create 300.0" in (stage2/"deposition.in").read_text()

def test_generated_stage_inputs_are_ordered_and_restartable(tmp_path):
    packed,_=packed_fixture(tmp_path)
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"ordered",packed_xyz=packed)
    for name in ("deposition.in","hold-300K.in","hold-400K.in","decompress-300K.in","decompress-400K.in","equilibrate-300K.in","anneal-400K.in"):
        text=(out/name).read_text()
        positions=[text.index(x) for x in ("boundary p p f","units real","atom_style full","read_data ","include ")]
        assert positions==sorted(positions)
        assert "dump trajectory" in text and "restart " in text and "write_restart" in text
        assert all(line.endswith(" nocoeff") for line in text.splitlines() if line.startswith("write_data "))
    deposition=(out/"deposition.in").read_text(); hold=(out/"hold-300K.in").read_text(); hold_400=(out/"hold-400K.in").read_text(); decompress_300=(out/"decompress-300K.in").read_text(); decompress_400=(out/"decompress-400K.in").read_text(); eq=(out/"equilibrate-300K.in").read_text(); anneal=(out/"anneal-400K.in").read_text()
    assert "variable zend equal 69.615" in deposition
    assert "processors * * 1" in deposition
    assert "comm_modify cutoff 20.0" in deposition
    assert deposition.index("fix walllo all wall/lj93") < deposition.index("minimize ")
    assert deposition.count("fix walllo all wall/lj93")==1
    assert "fix walllo all wall/lj93 zlo EDGE" in deposition
    assert "fix_modify walllo energy yes" in deposition
    assert "min_style fire" in deposition
    assert "minimize 0.0 1.0 20000 200000" in deposition
    assert deposition.index("minimize ") < deposition.index("file optimization-summary.txt screen yes") < deposition.index("write_data optimized.data nocoeff") < deposition.index("reset_timestep 0")
    assert "write_restart optimized.restart" in deposition
    assert "write_dump all custom optimized.lammpstrj" in deposition
    assert "file optimization-summary.txt screen yes" in deposition
    assert "force_target_met=${optimization_force_target_met}" in deposition
    assert "thermo_style custom step temp pe ke etotal press pxx pyy lx ly lz fnorm fmax" in deposition
    assert "reset_timestep 0" in deposition
    assert "thermo_style custom step temp pe ke etotal press pxx pyy lx ly lz v_zwall fnorm fmax" in deposition
    assert "timestep 0.5\nrun 10 start 0 stop 10" in deposition
    assert deposition.count("run 10 start 0 stop 10")==1
    assert "read_data deposited.data" in hold
    assert "variable hold_wall_hi equal 69.615" in hold
    assert "fix hi all wall/lj126 zhi ${hold_wall_hi}" in hold
    assert "fix ensemble all npt temp 300.0 300.0" in hold
    assert "timestep 0.5\nrun 1000000" in hold
    assert "write_data held-300K.data nocoeff" in hold
    assert "read_data deposited.data" in hold_400
    assert "variable hold_wall_hi equal 69.615" in hold_400
    assert "fix hi all wall/lj126 zhi ${hold_wall_hi}" in hold_400
    assert "fix heating all npt temp 300.0 400.0" in hold_400
    assert "timestep 0.5\nrun 400000" in hold_400
    assert "write_data heated-400K.data nocoeff" in hold_400
    assert "write_restart heated-400K.restart" in hold_400
    assert hold_400.index("write_data heated-400K.data nocoeff") < hold_400.index("fix hold all npt temp 400.0 400.0")
    assert "timestep 0.5\nrun 1000000" in hold_400
    assert "write_data held-400K.data nocoeff" in hold_400
    for text,source,suffix,temperature in (
        (decompress_300,"held-300K.data","300K","300.0"),
        (decompress_400,"held-400K.data","400K","400.0"),
    ):
        assert f"read_data {source}" in text
        assert "variable compressed_wall_hi equal 69.615" in text
        assert "variable decompression_wall_hi equal" in text
        assert "(step/400000.0)" in text
        assert "fix hi all wall/lj126 zhi v_decompression_wall_hi" in text
        assert f"fix decompress all npt temp {temperature} {temperature}" in text
        assert "timestep 0.5\nrun 400000 start 0 stop 400000" in text
        assert f"write_data decompressed-{suffix}.data nocoeff" in text
        assert "variable relaxed_wall_hi equal ${resolved_wall_hi}" in text
        assert "fix hi all wall/lj126 zhi ${relaxed_wall_hi}" in text
        assert f"fix relax all npt temp {temperature} {temperature}" in text
        assert "timestep 0.5\nrun 1000000" in text
        assert f"write_data relaxed-{suffix}.data nocoeff" in text
        assert text.index(f"write_data decompressed-{suffix}.data nocoeff") < text.index(f"write_data relaxed-{suffix}.data nocoeff")
    assert "read_data held-300K.data" in eq and "fix lo all wall/lj93 zlo EDGE" in eq
    assert "read_data equilibrated-300K.data" in anneal and "fix lo all wall/lj93 zlo EDGE" in anneal
    for text in (eq,anneal):
        assert "compute stage_zmax all reduce max z" in text
        assert "variable resolved_wall_hi equal ${rounded_wall}" in text
        assert "change_box all z final" in text
        assert "fix hi all wall/lj93 zhi ${resolved_wall_hi}" in text
    for name,source in (("deposition","topology_output.lmp"),("hold-300K","deposited.data"),("hold-400K","deposited.data"),("decompress-300K","held-300K.data"),("decompress-400K","held-400K.data"),("equilibrate-300K","held-300K.data"),("anneal-400K","equilibrated-300K.data")):
        smoke=(out/f"validate-{name}.in").read_text()
        assert f"read_data {source}" in smoke and "run 0" in smoke
    assert "run 0 post no # minimization replaced by force evaluation" in (out/"validate-deposition.in").read_text()

def test_rebuild_reuses_own_packed_xyz_and_refreshes_hashes(tmp_path):
    packed,_=packed_fixture(tmp_path/"fixture")
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"reuse",packed_xyz=packed)
    topology_hash=hashlib.sha256((out/"topology_output.lmp").read_bytes()).hexdigest()
    build(ROOT/"tests/data/small-study.toml",out,packed_xyz=out/"packed.xyz")
    assert hashlib.sha256((out/"topology_output.lmp").read_bytes()).hexdigest()==topology_hash
    from nio_md_prep.validate import validate
    validate(out)
    manifest=json.loads((out/"assembly_manifest.json").read_text())
    actual={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!="assembly_manifest.json"}
    assert manifest["output_hashes"]==actual

def test_refresh_inputs_preserves_existing_results(tmp_path):
    packed,_=packed_fixture(tmp_path/"fixture")
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"refresh",packed_xyz=packed)
    deposited=out/"deposited.data"; deposited.write_text("completed deposition sentinel\n")
    held_300=out/"held-300K.data"; held_300.write_text("completed 300 K sentinel\n")
    refresh_inputs(ROOT/"tests/data/small-study.toml",out)
    assert deposited.read_text()=="completed deposition sentinel\n"
    assert held_300.read_text()=="completed 300 K sentinel\n"
    assert (out/"hold-400K.in").exists()
    assert (out/"decompress-300K.in").exists()
    assert (out/"decompress-400K.in").exists()

def test_geometry_aware_wall_resolution_and_manual_override(tmp_path):
    data=parse(ROOT/"inputs/surfaces/corrugated-nio-110/surface.lmp")
    data.sections["Atoms"][0].fields[6]="265.0"; source=tmp_path/"deposited.data"; write(data,source)
    auto={"minimum":120.0,"clearance":30.0,"rounding":10.0,"box_margin":5.0,"manual_override":False}
    first=_wall_result(source,auto); second=_wall_result(source,auto)
    assert first==second and first["resolved_wall"]==300.0
    assert first["final_zhi"]==305.0 and first["measured_max_atom_z"]==265.0
    manual={**auto,"manual_override":True,"manual_value":320.0}
    forced=_wall_result(source,manual)
    assert forced["resolved_wall"]==320.0 and forced["final_zhi"]==325.0 and forced["mode"]=="manual"

@pytest.mark.parametrize("manual,wall,zhi",[(None,300.0,305.0),(320.0,320.0,325.0)])
def test_lammps_real_predecessor_wall_resolution(tmp_path,manual,wall,zhi):
    executable=shutil.which("lmp_mpi") or shutil.which("lmp") or shutil.which("lammps")
    if not executable: pytest.skip("LAMMPS executable unavailable")
    packed,_=packed_fixture(tmp_path/"fixture")
    config=ROOT/"tests/data/small-study.toml"
    if manual is not None:
        text=config.read_text().replace('production_upper_wall_mode = "auto"',f'production_upper_wall_mode = "auto"\nproduction_upper_wall = {manual}')
        config=tmp_path/"manual.toml"; config.write_text(text)
    out=build(config,tmp_path/"build",packed_xyz=packed)
    predecessor=parse(out/"topology_output.lmp")
    ligand=[a for a in predecessor.sections["Atoms"] if int(a.fields[1])>0]
    shift=265.0-max(float(a.fields[6]) for a in ligand)
    for atom in ligand: atom.fields[6]=f"{float(atom.fields[6])+shift:.8f}"
    predecessor.bounds["z"]=(predecessor.bounds["z"][0],270.0); write(predecessor,out/"held-300K.data")
    smoke=tmp_path/"smoke"; smoke.mkdir()
    for name in ("held-300K.data","force_field_settings_lammps_with_header.lmp","validate-equilibrate-300K.in"):
        shutil.copy2(out/name,smoke/name)
    before=parse(smoke/"held-300K.data"); before_xyz={int(a.fields[0]):tuple(map(float,a.fields[4:7])) for a in before.sections["Atoms"]}
    run=subprocess.run([executable,"-in","validate-equilibrate-300K.in"],cwd=smoke,capture_output=True,text=True)
    assert run.returncode==0,run.stderr+run.stdout[-2000:]
    assert f"Resolved production wall: {int(wall)}" in run.stdout
    result=parse(smoke/"equilibrated-300K.data")
    assert result.bounds["z"]==(before.bounds["z"][0],zhi)
    after_xyz={int(a.fields[0]):tuple(map(float,a.fields[4:7])) for a in result.sections["Atoms"]}
    assert after_xyz.keys()==before_xyz.keys()
    for atom_id,after in after_xyz.items(): assert after==pytest.approx(before_xyz[atom_id])
