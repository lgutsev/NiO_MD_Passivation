"""Result validation and trajectory diagnostics.

A result object is the last line of defence against an engine that returned
ghost-atom rows, a NaN energy or a truncated trajectory. These tests pin that
such output fails with a typed error, and that the energy diagnostics mean
what their names say: a drift is only reported where energy is conserved.
"""
import math

import pytest

from nio_md_prep.mlip.diagnostics import (
    DESCRIPTIVE_LABEL,
    descriptive_energy_change,
    ensemble_diagnostics,
    nve_drift,
    series_diagnostics,
    temperature_ndof,
)
from nio_md_prep.mlip.errors import ResultError
from nio_md_prep.mlip.results import PotentialResult, TrajectoryResult


def result(**overrides) -> PotentialResult:
    base = dict(
        energy_eV=-10.0,
        forces_eV_per_A=((0.1, 0.0, -0.1), (-0.1, 0.0, 0.1)),
        symbols=("Ni", "O"),
        energy_convention="total",
        engine="ase",
        potential="mock",
        implementation="mock-ase",
    )
    base.update(overrides)
    return PotentialResult(**base)


# --- PotentialResult ---------------------------------------------------------


@pytest.mark.parametrize(
    "overrides,message",
    [
        ({"energy_eV": math.nan}, "energy is not finite"),
        ({"energy_eV": math.inf}, "energy is not finite"),
        ({"forces_eV_per_A": ((0.1, 0.0), (0.0, 0.0, 0.0))}, "has 2 components"),
        ({"forces_eV_per_A": ((0.1, 0.0, math.inf), (0.0, 0.0, 0.0))}, "not finite"),
        # 3 rows for 2 atoms: what an nlocal+nghost array looks like.
        ({"forces_eV_per_A": ((0, 0, 0),) * 3}, "3 force rows for 2 atoms"),
        ({"symbols": (), "forces_eV_per_A": ()}, "at least one atom"),
        ({"stress_eV_per_A3": (0.0,) * 5}, "6-component Voigt"),
        ({"stress_eV_per_A3": (0.0,) * 5 + (math.nan,)}, "stress is not finite"),
        ({"per_atom_energy_eV": (-5.0, -3.0, -2.0)}, "3 per-atom energies for 2 atoms"),
        ({"per_atom_energy_eV": (-5.0, math.nan)}, "per-atom energy is not finite"),
    ],
)
def test_invalid_engine_output_is_a_typed_error(overrides, message):
    with pytest.raises(ResultError, match=message):
        result(**overrides)


def test_result_errors_are_still_value_errors():
    """Existing callers that caught ValueError keep working."""
    with pytest.raises(ValueError):
        result(energy_eV=math.nan)


def test_a_valid_result_is_normalised_to_floats_and_tuples():
    r = result(forces_eV_per_A=[[1, 0, 0], [0, 0, -1]], per_atom_energy_eV=[-6, -4])
    assert r.forces_eV_per_A == ((1.0, 0.0, 0.0), (0.0, 0.0, -1.0))
    assert r.per_atom_energy_eV == (-6.0, -4.0)
    assert r.max_force_eV_per_A == pytest.approx(1.0)


# --- TrajectoryResult --------------------------------------------------------


def trajectory(**overrides) -> TrajectoryResult:
    base = dict(
        steps=10,
        timestep_fs=1.0,
        ensemble="nve",
        trajectory_path=None,
        log_path=None,
        initial=result(),
        final=result(),
    )
    base.update(overrides)
    return TrajectoryResult(**base)


def test_trajectory_records_the_final_geometry_in_the_source_basis():
    t = trajectory(
        steps_completed=10,
        frames_written=3,
        final_positions_angstrom=[[0, 0, 0], [1, 1, 1]],
        final_cell_angstrom=[[4, 0, 0], [0, 4, 0], [0, 0, 4]],
        final_pbc=[True, True, False],
        integrator_resolved={"class": "VelocityVerlet"},
    )
    assert t.frames == 3  # backwards-compatible alias
    assert t.final_pbc == (True, True, False)
    record = t.as_dict()
    assert record["steps_completed"] == 10 and record["frames_written"] == 3
    assert record["final_cell_angstrom"][0] == [4.0, 0.0, 0.0]
    assert "final_positions_angstrom" not in record
    assert t.as_dict(include_arrays=True)["final_positions_angstrom"][1] == [1.0, 1.0, 1.0]


@pytest.mark.parametrize(
    "overrides,message",
    [
        ({"final": result(symbols=("O", "Ni"))}, "same order"),
        ({"final_positions_angstrom": [[0, 0, 0]]}, "2 rows of 3"),
        ({"final_positions_angstrom": [[0, 0, math.nan], [0, 0, 0]]}, "not finite"),
        ({"final_cell_angstrom": [[1, 0, 0], [0, 1, 0]]}, "3 rows of 3"),
        ({"final_pbc": [True, True]}, "three entries"),
        ({"frames_written": -1}, "non-negative integer"),
        ({"steps_completed": 2.5}, "non-negative integer"),
    ],
)
def test_inconsistent_trajectory_output_is_refused(overrides, message):
    with pytest.raises(ResultError, match=message):
        trajectory(**overrides)


# --- diagnostics -------------------------------------------------------------


def test_nve_drift_is_a_per_time_linear_fit():
    """E/N rising 2 meV per ps: drift 0.002 eV/atom/ps whatever the sampling."""
    n_atoms = 10
    times_fs = [0.0, 250.0, 500.0, 750.0, 1000.0]
    energies = [n_atoms * (-5.0 + 0.002 * t / 1000.0) for t in times_fs]
    report = nve_drift(times_fs, energies, n_atoms)
    assert report["energy_drift_eV_per_atom_per_ps"] == pytest.approx(0.002, rel=1e-12)
    assert report["max_abs_energy_excursion_eV_per_atom"] == pytest.approx(0.002, rel=1e-12)
    assert report["samples"] == 5 and report["duration_fs"] == 1000.0


def test_a_bounded_oscillation_is_not_a_drift():
    """A symplectic integrator's energy oscillates without trend.

    Ending a quarter period in (at a crest), the endpoint difference would
    read as ~1e-4 eV/atom/ps; the linear fit is an order of magnitude smaller.
    """
    times_fs = [float(t) for t in range(0, 1026, 5)]
    energies = [-50.0 + 1e-3 * math.sin(2 * math.pi * t / 100.0) for t in times_fs]
    report = nve_drift(times_fs, energies, 10)
    endpoint_rate = (energies[-1] - energies[0]) / 10 / (times_fs[-1] / 1000.0)
    assert endpoint_rate == pytest.approx(9.76e-5, rel=1e-3)
    assert abs(report["energy_drift_eV_per_atom_per_ps"]) < 0.1 * endpoint_rate
    assert report["max_abs_energy_excursion_eV_per_atom"] == pytest.approx(1e-4, rel=1e-3)


def test_thermostatted_runs_report_a_labelled_change_not_a_drift():
    times_fs = [0.0, 10.0, 20.0]
    energies = [-10.0, -9.0, -8.0]
    report = ensemble_diagnostics("nvt", times_fs, energies, 2)
    assert report["conservation_test"] is False
    assert report["label"] == DESCRIPTIVE_LABEL
    assert report["total_energy_change_eV_per_atom"] == pytest.approx(1.0)
    assert "energy_drift_eV_per_atom_per_ps" not in report
    assert descriptive_energy_change(times_fs, energies, 2)["label"] == DESCRIPTIVE_LABEL


def test_a_conserved_quantity_restores_a_conservation_test():
    """LAMMPS econserve / ASE NHC conserved energy: fitted like an NVE energy."""
    report = ensemble_diagnostics(
        "nvt",
        [0.0, 500.0, 1000.0],
        [-10.0, -9.0, -8.0],
        2,
        conserved_energy_eV=[-12.0, -12.0, -12.0],
        conserved_quantity="econserve",
    )
    assert report["conservation_test"] is True
    assert report["conserved_quantity"] == "econserve"
    assert report["conserved_quantity_drift"]["energy_drift_eV_per_atom_per_ps"] == 0.0
    with pytest.raises(ResultError, match="name of the quantity"):
        ensemble_diagnostics("nvt", [0.0, 1.0], [0.0, 0.0], 1, conserved_energy_eV=[0.0, 0.0])


def test_nve_reports_a_conservation_test():
    report = ensemble_diagnostics("nve", [0.0, 1.0], [-1.0, -1.0], 1)
    assert report["conservation_test"] is True and report["ensemble"] == "nve"


@pytest.mark.parametrize(
    "times,energies,n_atoms,message",
    [
        ([0.0], [1.0], 1, "at least two samples"),
        ([0.0, 1.0], [1.0], 1, "must align"),
        ([0.0, 0.0, 1.0], [1.0, 1.0, 1.0], 1, "strictly increasing"),
        ([0.0, 1.0], [1.0, math.nan], 1, "non-finite"),
        ([0.0, 1.0], [1.0, 1.0], 0, "positive integer"),
    ],
)
def test_malformed_series_are_refused(times, energies, n_atoms, message):
    with pytest.raises(ResultError, match=message):
        nve_drift(times, energies, n_atoms)


def test_series_payload_shape():
    report = series_diagnostics("nve", {"time_fs": [0, 1], "total_energy_eV": [0, 0]}, 1)
    assert report["energy_drift_eV_per_atom_per_ps"] == 0.0
    with pytest.raises(ResultError, match="time_fs"):
        series_diagnostics("nve", {"total_energy_eV": [0, 0]}, 1)


def test_temperature_degrees_of_freedom_convention():
    assert temperature_ndof(16) == 48
    assert temperature_ndof(16, com_removed=True) == 45
    # With frozen atoms the centre of mass is not free; nothing more is removed.
    assert temperature_ndof(16, n_fixed=4, com_removed=True) == 36
    with pytest.raises(ResultError, match="no kinetic degrees"):
        temperature_ndof(1, com_removed=True)
    with pytest.raises(ResultError, match="n_fixed"):
        temperature_ndof(2, n_fixed=3)
