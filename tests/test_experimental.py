from pathlib import Path
import pytest

from nio_md_prep.analysis.experimental import (
    aggregate_device_results,
    load_device_results,
)


def test_load_device_results_parses_rows(tmp_path):
    path = tmp_path / "jv.csv"
    path.write_text(
        "secondary,assembly,batch,scan_direction,voc_v,jsc_ma_cm2,ff_percent,pce_percent\n"
        "Test-SAM,CoSAM,synthetic-a,R,1.100,10.000,70.0,7.700\n"
        "Test-SAM,CoSAM,synthetic-a,F,1.110,10.100,71.0,7.960\n",
        encoding="utf-8",
    )
    rows = load_device_results(path)
    assert len(rows) == 2
    assert rows[0]["secondary"] == "Test-SAM"
    assert rows[0]["assembly"] == "CoSAM"
    assert rows[0]["voc_v"] == pytest.approx(1.100)
    assert rows[1]["pce_percent"] == pytest.approx(7.960)


def test_load_device_results_requires_columns(tmp_path):
    path = tmp_path / "jv.csv"
    path.write_text("secondary,assembly\nTest-SAM,CoSAM\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing required column"):
        load_device_results(path)


def test_load_device_results_reports_bad_numeric_value(tmp_path):
    path = tmp_path / "jv.csv"
    path.write_text(
        "secondary,assembly,voc_v,jsc_ma_cm2,ff_percent,pce_percent\n"
        "Test-SAM,CoSAM,not-a-number,10.000,70.0,7.700\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="line 2"):
        load_device_results(path)


def test_aggregate_device_results_pools_replicates_and_batches():
    rows = [
        {
            "secondary": "Test-SAM",
            "assembly": "CoSAM",
            "batch": "a",
            "scan_direction": "R",
            "voc_v": 1.100,
            "jsc_ma_cm2": 10.000,
            "ff_percent": 70.0,
            "pce_percent": 7.700,
        },
        {
            "secondary": "Test-SAM",
            "assembly": "CoSAM",
            "batch": "b",
            "scan_direction": "",
            "voc_v": 1.050,
            "jsc_ma_cm2": 9.000,
            "ff_percent": 68.0,
            "pce_percent": 6.400,
        },
        {
            "secondary": "Test-SAM",
            "assembly": "Sequential",
            "batch": "a",
            "scan_direction": "R",
            "voc_v": 1.200,
            "jsc_ma_cm2": 12.000,
            "ff_percent": 80.0,
            "pce_percent": 11.500,
        },
    ]
    aggregated = aggregate_device_results(rows)
    assert set(aggregated) == {("Test-SAM", "CoSAM"), ("Test-SAM", "Sequential")}
    cosam = aggregated[("Test-SAM", "CoSAM")]
    assert cosam["n"] == 2
    assert cosam["pce_percent"]["mean"] == pytest.approx((7.700 + 6.400) / 2)
    assert cosam["pce_percent"]["std"] > 0.0
    sequential = aggregated[("Test-SAM", "Sequential")]
    assert sequential["n"] == 1
    assert sequential["pce_percent"]["std"] == 0.0
