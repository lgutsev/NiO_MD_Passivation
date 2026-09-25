"""Cross-engine tolerance records: candidates are measurements, references are reviewed.

No backend is needed; the comparison report is a hand-built dict with the same
keys ``jobs.compare_engines`` produces.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from nio_md_prep.mlip import tolerances

ROOT = Path(__file__).resolve().parents[1]


def _report(delta_e=-2.0e-6, rmse=3.0e-5, vector=9.0e-5, component=7.0e-5, stress=None):
    return {
        "structure": {"sha256": "abc123"},
        "energy_convention": "total",
        "comparisons": {
            "ase->lammps": {
                "delta_energy_per_atom_eV": delta_e,
                "force_rmse_eV_per_A": rmse,
                "force_max_abs_error_eV_per_A": vector,
                "force_max_component_error_eV_per_A": component,
                "stress_max_abs_error_eV_per_A3": stress,
                "notes": [],
            }
        },
    }


def test_measured_metrics_takes_the_absolute_energy_difference_and_keeps_missing_stress_missing():
    measured = tolerances.measured_metrics(_report())
    row = measured["ase->lammps"]
    assert row["delta_energy_per_atom_eV"] == pytest.approx(2.0e-6)
    assert row["force_max_abs_error_eV_per_A"] == pytest.approx(9.0e-5)
    assert row["force_max_component_error_eV_per_A"] == pytest.approx(7.0e-5)
    assert row["stress_max_abs_error_eV_per_A3"] is None  # absent, never an invented zero


def test_a_candidate_is_never_accepted_as_a_reference(tmp_path):
    record = tolerances.candidate_record(_report(), measured_on={"host": "test"})
    assert record["status"] == tolerances.CANDIDATE_STATUS
    assert record["reviewed_by"] is None
    path = tmp_path / "tolerances.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="only a reviewed file"):
        tolerances.load_reference(path)


def test_a_reviewed_reference_loads_and_flags_exceedances(tmp_path):
    record = tolerances.candidate_record(_report(), measured_on={"host": "test"})
    record["status"] = tolerances.ACCEPTED_STATUS
    record["reviewed_by"] = "reviewer"
    path = tmp_path / "tolerances.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    reference = tolerances.load_reference(path)

    within = tolerances.measured_metrics(_report())
    assert tolerances.exceedances(within, reference) == []

    worse = tolerances.measured_metrics(_report(vector=1.0e-3))
    problems = tolerances.exceedances(worse, reference)
    assert len(problems) == 1 and "force_max_abs_error_eV_per_A" in problems[0]

    other_pair = {"ase->openmm": within["ase->lammps"]}
    assert "no reviewed tolerance" in tolerances.exceedances(other_pair, reference)[0]


@pytest.mark.parametrize(
    "mutation, message",
    [
        (lambda r: r.pop("reviewed_by"), "reviewer"),
        (lambda r: r.pop("measured_on"), "measured_on"),
        (lambda r: r.update(tolerances={}), "non-empty"),
        (lambda r: r.update(tolerances={"ase-lammps": {}}), "reference->candidate"),
        (lambda r: r["tolerances"]["ase->lammps"].update(bogus=1.0), "unknown tolerance metric"),
        (lambda r: r["tolerances"]["ase->lammps"].update(force_rmse_eV_per_A=-1.0), ">= 0"),
    ],
)
def test_malformed_references_fail_loudly(mutation, message):
    record = tolerances.candidate_record(_report(), measured_on={"host": "test"})
    record["status"] = tolerances.ACCEPTED_STATUS
    record["reviewed_by"] = "reviewer"
    mutation(record)
    with pytest.raises(ValueError, match=message):
        tolerances.validate_reference(record)


def test_absent_reference_is_none(tmp_path):
    assert tolerances.load_reference(tmp_path / "missing.json") is None


def _parity_script():
    path = ROOT / "scripts" / "mlip_mace_ase_lammps_parity.py"
    spec = importlib.util.spec_from_file_location("mlip_parity_script", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parity_script_parses_per_axis_pbc_and_uses_maces_export_name():
    script = _parity_script()
    assert script.parse_bools("T,T,F") == (True, True, False)
    assert script.parse_bools("true, false, 1") == (True, False, True)
    with pytest.raises(Exception):
        script.parse_bools("T,T")
    assert script.MLIAP_SUFFIX == "-mliap_lammps.pt"


def test_parity_script_refuses_a_non_empty_output_directory(tmp_path):
    script = _parity_script()
    out = tmp_path / "out"
    out.mkdir()
    (out / "old.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit):
        script.main(["--model", str(tmp_path / "m.model"), "--elements", "Ni,O",
                     "--structure", str(tmp_path / "s.xyz"), "--output", str(out)])
