from pathlib import Path
import pytest
from nio_md_prep.build import build
from nio_md_prep.geometry import elements
from nio_md_prep.lammps import parse

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
    assert (out/"validation_report.txt").read_text().startswith("VALID")
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

def test_missing_packmol_is_incomplete(tmp_path,monkeypatch):
    monkeypatch.setattr("nio_md_prep.build.shutil.which",lambda _:None)
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"incomplete")
    assert not (out/"topology_output.lmp").exists()
    assert "packmol < packmol.inp" in (out/"validation_report.txt").read_text()

def test_sequential_stage_preserves_existing_records(tmp_path):
    packed,template=packed_fixture(tmp_path/"one")
    stage1=build(ROOT/"tests/data/small-study.toml",tmp_path/"stage1",packed_xyz=packed)
    before=parse(stage1/"topology_output.lmp")
    packed2,_=packed_fixture(tmp_path/"two",shift=(20.0,20.0,100.0))
    stage2=build(ROOT/"tests/data/small-study.toml",tmp_path/"stage2",primary_final=stage1/"topology_output.lmp",packed_xyz=packed2)
    after=parse(stage2/"topology_output.lmp")
    for section in ("Atoms","Bonds","Angles","Dihedrals","Impropers"):
        assert [r.fields for r in after.sections[section]][:before.count(section)]==[r.fields for r in before.sections[section]]
    assert after.count("Atoms")==before.count("Atoms")+45
    assert after.bounds==before.bounds
    assert "inside box 2.000000 2.000000 67.509910 123.100000 39.700000 265.000000" in (stage2/"packmol.inp").read_text()

def test_generated_stage_inputs_are_ordered_and_restartable(tmp_path):
    packed,_=packed_fixture(tmp_path)
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"ordered",packed_xyz=packed)
    for name in ("deposition.in","equilibrate-300K.in","anneal-400K.in"):
        text=(out/name).read_text()
        positions=[text.index(x) for x in ("boundary p p f","units real","atom_style full","read_data ","include ")]
        assert positions==sorted(positions)
        assert "dump trajectory" in text and "restart " in text and "write_restart" in text
