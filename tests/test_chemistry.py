from pathlib import Path
import pytest
from nio_md_prep.chemistry import molecular_weight, phosphonate_roles, corrected_dihedral_types
from nio_md_prep.lammps import parse
from nio_md_prep.config import load

ROOT=Path(__file__).parents[1]

def test_ligpargen_molecular_weights_and_anchor_counts():
    expected={
        "me-4pacz":(331.352,1),
        "meo-2pacz":(335.296,1),
        "meo-4padbc":(463.470,1),
        "dcz-4p":(604.580,2),
    }
    for slug,(mw,anchors) in expected.items():
        data=parse(ROOT/"inputs"/"molecules"/slug/"ligpargen.lmp")
        assert molecular_weight(data)==pytest.approx(mw)
        assert list(phosphonate_roles(data).values()).count("P")==anchors
        corrected=corrected_dihedral_types(data)
        assert sum(int(row.fields[1]) in corrected for row in data.sections["Dihedrals"])==5*anchors

def test_equal_volume_stock_counts():
    primary_mw=331.352
    for secondary_mw,wanted in ((335.296,107),(463.470,77),(604.580,59)):
        exact=180*(0.3/secondary_mw)/(0.5/primary_mw)
        assert round(exact)==wanted

def test_all_study_total_atom_counts():
    molecule_atoms={"me-4pacz":45,"meo-2pacz":41,"meo-4padbc":59,"dcz-4p":76}
    expected={
        "me-4pacz-alone.toml":29160,
        "me-4pacz-high-dose.toml":32805,
        "me-4pacz-meo-2pacz-cosam.toml":33547,
        "me-4pacz-meo-4padbc-cosam.toml":33703,
        "me-4pacz-dcz-4p-cosam.toml":33644,
        "me-4pacz-then-meo-2pacz.toml":33547,
        "me-4pacz-then-meo-4padbc.toml":33703,
        "me-4pacz-then-dcz-4p.toml":33644,
    }
    for name,total in expected.items():
        cfg=load(ROOT/"studies"/name); ligand=sum(molecule_atoms[x["slug"]]*x["count"] for x in cfg["molecules"])
        ligand += 45*int(cfg.get("composition",{}).get("primary_count",0))
        assert 21060+ligand==total

def test_authoritative_charmm_phosphonate_parameter():
    text=(ROOT/"examples/CORRECTED_FORCE_FIELDS/force_field_settings_lammps_Me4PACz_MPTMS-OH.lmp").read_text()
    corrected=[line for line in text.splitlines() if line.strip().startswith("dihedral_coeff") and " charmm " in line]
    assert len(corrected)==5
    assert all("charmm 0.5333 0 3 0.000" in line for line in corrected)
