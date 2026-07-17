from pathlib import Path
import pytest
from nio_md_prep.chemistry import molecular_weight, phosphonate_roles
from nio_md_prep.lammps import parse

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

def test_equal_volume_stock_counts():
    primary_mw=331.352
    for secondary_mw,wanted in ((335.296,107),(463.470,77),(604.580,59)):
        exact=180*(0.3/secondary_mw)/(0.5/primary_mw)
        assert round(exact)==wanted
