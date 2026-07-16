from pathlib import Path
from nio_md_prep.build import build
from nio_md_prep.lammps import parse

ROOT=Path(__file__).parents[1]

def test_small_offline_build(tmp_path):
    out=build(ROOT/"tests/data/small-study.toml",tmp_path/"build")
    data=parse(out/"topology_output.lmp")
    assert data.count("Atoms")==12960+45
    assert (out/"type_map.json").exists()
    assert (out/"validation_report.txt").read_text().startswith("VALID")
    assert "special_bonds   amber" in (out/"force_field_settings_lammps_with_header.lmp").read_text()
