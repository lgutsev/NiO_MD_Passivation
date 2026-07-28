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
    assert rows[0]["device_id"] == ""
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


def test_aggregate_device_results_weights_batches_equally():
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


def test_aggregate_device_results_collapses_forward_reverse_scans():
    rows = [
        {
            "secondary": "Test-SAM",
            "assembly": "CoSAM",
            "batch": "a",
            "device_id": "",
            "scan_direction": "R",
            "voc_v": 1.0,
            "jsc_ma_cm2": 10.0,
            "ff_percent": 70.0,
            "pce_percent": 7.0,
        },
        {
            "secondary": "Test-SAM",
            "assembly": "CoSAM",
            "batch": "a",
            "device_id": "",
            "scan_direction": "F",
            "voc_v": 1.2,
            "jsc_ma_cm2": 12.0,
            "ff_percent": 80.0,
            "pce_percent": 9.0,
        },
    ]

    result = aggregate_device_results(rows)[("Test-SAM", "CoSAM")]
    assert result["n_measurements"] == 2
    assert result["n_independent_units"] == 1
    assert result["n_batches"] == 1
    assert result["inferred_scan_pair_count"] == 1
    assert result["pce_percent"]["mean"] == pytest.approx(8.0)
    assert result["pce_percent"]["std"] == pytest.approx(0.0)


def test_aggregate_device_results_averages_devices_then_batches():
    def row(batch, device, value):
        return {
            "secondary": "Test-SAM",
            "assembly": "Sequential",
            "batch": batch,
            "device_id": device,
            "scan_direction": "",
            "voc_v": value,
            "jsc_ma_cm2": value,
            "ff_percent": value,
            "pce_percent": value,
        }

    # Batch a has three devices near 1; batch b has one device at 9.
    result = aggregate_device_results(
        [row("a", "a1", 1.0), row("a", "a2", 1.0), row("a", "a3", 1.0),
         row("b", "b1", 9.0)]
    )[("Test-SAM", "Sequential")]

    # Equal batch weighting gives (1 + 9) / 2, not the row-pooled value 3.
    assert result["pce_percent"]["mean"] == pytest.approx(5.0)
    assert result["n_independent_units"] == 4
    assert result["n_batches"] == 2
    assert result["pce_percent"]["uncertainty_basis"] == "between_batch"
