from pathlib import Path
import pytest
from nio_md_prep.geometry import elements
from nio_md_prep.lammps import parse, replicate

ROOT=Path(__file__).parents[1]

def test_golden_ligpargen_counts():
    a=parse(ROOT/"scripts/mnt/data/ligand_oplsaa.lmp")
    b=parse(ROOT/"scripts/mnt/data/ligand2_oplsaa.lmp")
    assert (a.count("Atoms"),b.count("Atoms"))==(45,18)
    assert (a.type_count("atom"),b.type_count("atom"))==(45,18)

def test_independent_offsets_and_connectivity():
    d=parse(ROOT/"scripts/mnt/data/ligand2_oplsaa.lmp")
    coords=[(float(r.fields[4]),float(r.fields[5]),float(r.fields[6])) for r in d.sections["Atoms"]]*2
    rows,inc=replicate(d,2,{"atom":45,"bond":47,"angle":83,"dihedral":120,"improper":27},{"atom":90,"bond":94,"angle":166,"dihedral":240,"improper":54},2,coords)
    assert len(rows["Atoms"])==36 and inc["bond"]==34
    assert min(int(r.fields[2]) for r in rows["Atoms"])==46
    assert {int(x) for r in rows["Bonds"] for x in r.fields[2:4]} <= {int(r.fields[0]) for r in rows["Atoms"]}

def test_example_180_360_counts():
    d=parse(ROOT/"examples/Sigma_Big_NiO_110_Multi_Me4PACz_4x_180-360_swall/topology_output.lmp")
    assert [d.count(x) for x in ("Atoms","Bonds","Angles","Dihedrals","Impropers")]==[35640,14580,25020,35640,7740]
    assert [d.type_count(x) for x in ("atom","bond","angle","dihedral","improper")]==[65,64,111,159,35]

def test_historical_unique_molecule_ids():
    primary=parse(ROOT/"scripts/mnt/data/ligand_oplsaa.lmp")
    secondary=parse(ROOT/"scripts/mnt/data/ligand2_oplsaa.lmp")
    pxyz=[tuple(map(float,r.fields[4:7])) for r in primary.sections["Atoms"]]*180
    prows,pinc=replicate(primary,180,{k:0 for k in ("atom","bond","angle","dihedral","improper")},{k:0 for k in ("atom","bond","angle","dihedral","improper")},0,pxyz)
    offsets={k:primary.type_count(k) for k in ("atom","bond","angle","dihedral","improper")}
    ids={"atom":pinc["atom"],"bond":pinc["bond"],"angle":pinc["angle"],"dihedral":pinc["dihedral"],"improper":pinc["improper"]}
    sxyz=[tuple(map(float,r.fields[4:7])) for r in secondary.sections["Atoms"]]*360
    srows,_=replicate(secondary,360,offsets,ids,180,sxyz)
    mids={int(r.fields[1]) for r in prows["Atoms"]+srows["Atoms"]}
    assert mids==set(range(1,541))
    assert len(prows["Atoms"])==8100 and len(srows["Atoms"])==6480

def test_missing_and_malformed(tmp_path):
    with pytest.raises(FileNotFoundError): parse(tmp_path/"missing.lmp")
    p=tmp_path/"bad.lmp"; p.write_text("LAMMPS\n\nMasses\n\n1 12.011\n")
    with pytest.raises(ValueError): parse(p)

def test_ambiguous_mass(tmp_path):
    p=tmp_path/"ambiguous.lmp"; p.write_text("LAMMPS\n\nMasses\n\n1 13.0\n\nAtoms # full\n\n1 1 1 0 0 0 0\n")
    with pytest.raises(ValueError,match="geometry.xyz"): elements(parse(p))
