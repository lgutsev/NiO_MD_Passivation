"""The whole pipeline, end to end, on the mock potential.

This is what makes the subsystem testable in ordinary CI: no torch, no
LAMMPS, no OpenMM, but every layer exercised -- configuration parsing, bridge
resolution, element-coverage validation, capability negotiation, execution,
canonical units, energy conventions, and a written manifest.
"""
import json
import textwrap
from pathlib import Path

import pytest

from nio_md_prep.mlip import jobs
from nio_md_prep.mlip.config import parse_job
from nio_md_prep.mlip.errors import (
    CapabilityError,
    ConfigError,
    ElementCoverageError,
)
from nio_md_prep.mlip.provenance import MANIFEST_NAME, read_manifest
from nio_md_prep.mlip.results import compare_results
from nio_md_prep.mlip.specs import StructureSpec

pytest.importorskip("ase")

MOCK_CONFIG = """
[potential]
kind = "mock"
elements = ["Ni", "O"]
epsilon_eV = 0.05
sigma_angstrom = 2.2
cutoff_angstrom = 6.0

[engine]
kind = "ase"

[simulation]
task = "singlepoint"
compute_stress = true
compute_per_atom_energy = true
"""

MOCK_MD_CONFIG = """
[potential]
kind = "mock"
elements = ["Ni", "O"]
epsilon_eV = 0.05
sigma_angstrom = 2.2

[engine]
kind = "ase"

[simulation]
task = "md"
ensemble = "nve"
timestep_fs = 0.5
steps = 20
seed = 12345
trajectory_interval = 5
log_interval = 5
"""


@pytest.fixture
def structure_file(tmp_path, rattled_nio_structure):
    from ase.io import write

    path = tmp_path / "nio.xyz"
    write(str(path), rattled_nio_structure, format="extxyz")
    return path


def write_config(tmp_path: Path, text: str, structure: Path | None = None) -> Path:
    body = textwrap.dedent(text)
    if structure is not None:
        body += f'\n[structure]\npath = "{structure.as_posix()}"\n'
    path = tmp_path / "job.toml"
    path.write_text(body, encoding="utf-8")
    return path


# --- validate -------------------------------------------------------------


def test_validate_executes_nothing_and_reports_the_route(tmp_path, structure_file):
    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    report = jobs.validate_job(job)
    assert report["ok"]
    assert report["route"]["implementation"] == "mock-ase"
    assert report["capabilities"]["energy"] and report["capabilities"]["forces"]
    assert report["structure"]["n_atoms"] == 16
    # Nothing was written: validation is a read-only operation.
    assert not list(tmp_path.glob("*.traj"))
    assert not (tmp_path / MANIFEST_NAME).exists()


def test_validate_rejects_an_uncovered_element(tmp_path, rattled_nio_structure):
    """A Ni/O model handed a fluorine must stop before anything runs."""
    from ase import Atom
    from ase.io import write

    contaminated = rattled_nio_structure.copy()
    contaminated.append(Atom("F", (0.5, 0.5, 0.5)))
    path = tmp_path / "contaminated.xyz"
    write(str(path), contaminated, format="extxyz")

    job = parse_job(write_config(tmp_path, MOCK_CONFIG, path))
    with pytest.raises(ElementCoverageError, match="F"):
        jobs.validate_job(job)


def test_validate_rejects_npt_on_openmm_before_anything_is_generated(
    tmp_path, structure_file
):
    """The real route, rejected on capabilities alone -- with no OpenMM installed.

    ``validate`` resolves and negotiates without importing a backend, so this
    runs in ordinary CI and proves the rejection happens before a batch job
    could exist.
    """
    config = textwrap.dedent(
        """
        [potential]
        kind = "mace"
        model = "nio.model"
        elements = ["Ni", "O"]

        [engine]
        kind = "openmm"

        [simulation]
        task = "md"
        ensemble = "npt"
        temperature_K = 300
        pressure_bar = 1.0
        timestep_fs = 0.5
        steps = 10
        """
    )
    job = parse_job(write_config(tmp_path, config, structure_file))
    with pytest.raises(CapabilityError) as excinfo:
        jobs.validate_job(job)
    message = str(excinfo.value)
    assert "stress" in message
    assert "integrates the simulation cell" in message


def test_validate_accepts_npt_on_a_route_that_reports_a_virial(tmp_path, structure_file):
    config = MOCK_CONFIG.replace(
        'task = "singlepoint"\ncompute_stress = true\ncompute_per_atom_energy = true',
        'task = "md"\nensemble = "npt"\ntemperature_K = 300\npressure_bar = 1.0\n'
        "timestep_fs = 0.5\nsteps = 10",
    )
    job = parse_job(write_config(tmp_path, config, structure_file))
    report = jobs.validate_job(job)
    assert report["ok"] and report["requirements"]["stress"] is True


# --- singlepoint ----------------------------------------------------------


def test_singlepoint_writes_canonical_results_and_a_manifest(tmp_path, structure_file):
    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    output = tmp_path / "run"
    report = jobs.run_singlepoint(job, output_dir=output)

    result = report["result"]
    assert result.n_atoms == 16
    assert result.native_units == "ase"
    assert result.energy_convention == "total"
    assert result.stress_eV_per_A3 is not None and len(result.stress_eV_per_A3) == 6
    assert result.per_atom_energy_eV is not None
    # A rattled cell has real forces; the site energies must sum to the total.
    assert result.max_force_eV_per_A > 0
    assert sum(result.per_atom_energy_eV) == pytest.approx(result.energy_eV, rel=1e-10)

    manifest = read_manifest(output)
    assert manifest["bridge"]["implementation"] == "mock-ase"
    assert manifest["structure"]["sha256"] == report["structure"]["sha256"]
    assert manifest["units"]["canonical_energy"] == "eV"
    assert manifest["results"]["singlepoint"]["energy_eV"] == pytest.approx(
        result.energy_eV
    )
    json.dumps(manifest)


def test_singlepoint_without_a_structure_says_so(tmp_path):
    job = parse_job(write_config(tmp_path, MOCK_CONFIG))
    with pytest.raises(ConfigError, match="needs a structure"):
        jobs.run_singlepoint(job, output_dir=tmp_path / "run")


def test_a_structure_may_be_supplied_separately(tmp_path, structure_file):
    job = parse_job(write_config(tmp_path, MOCK_CONFIG))
    report = jobs.run_singlepoint(
        job, output_dir=tmp_path / "run", structure=StructureSpec(path=structure_file)
    )
    assert report["result"].n_atoms == 16


# --- smoke-md -------------------------------------------------------------


def test_smoke_md_runs_a_short_trajectory_and_reports_drift(tmp_path, structure_file):
    job = parse_job(write_config(tmp_path, MOCK_MD_CONFIG, structure_file))
    output = tmp_path / "md"
    report = jobs.run_smoke_md(job, output_dir=output)

    trajectory = report["trajectory"]
    assert trajectory.steps == 20
    assert trajectory.frames >= 2
    assert Path(trajectory.trajectory_path).exists()
    assert Path(trajectory.log_path).exists()
    # NVE on a conservative analytic potential: drift is the diagnostic.
    assert trajectory.total_energy_drift_eV_per_atom is not None
    assert abs(trajectory.total_energy_drift_eV_per_atom) < 1e-2
    assert read_manifest(output)["results"]["trajectory"]["steps"] == 20


def test_smoke_md_refuses_to_become_a_production_run(tmp_path, structure_file):
    """The step ceiling is the point: this command is a diagnostic."""
    config = MOCK_MD_CONFIG.replace("steps = 20", "steps = 100000")
    job = parse_job(write_config(tmp_path, config, structure_file))
    with pytest.raises(ConfigError, match="not a replacement"):
        jobs.run_smoke_md(job, output_dir=tmp_path / "md")


def test_smoke_md_requires_an_md_task(tmp_path, structure_file):
    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    with pytest.raises(ConfigError, match="task = 'md'"):
        jobs.run_smoke_md(job, output_dir=tmp_path / "md")


# --- inspect --------------------------------------------------------------


def test_inspect_reports_the_matrix_without_importing_a_backend():
    import sys

    report = jobs.inspect_environment()
    assert "mace" in report["matrix"] and "lammps" in report["matrix"]
    unsupported = [r for r in report["routes"] if r["status"] != "supported"]
    assert unsupported and unsupported[0]["potential"] == "lammps"
    assert unsupported[0]["engine"] == "openmm"
    assert "torch" not in sys.modules


def test_inspect_with_a_config_reports_the_selected_route(tmp_path, structure_file):
    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    report = jobs.inspect_environment(job)
    assert report["selected_route"]["implementation"] == "mock-ase"
    assert report["route_availability"]["available"] is True
    assert report["potential"]["kind"] == "mock"


# --- comparison machinery -------------------------------------------------


def test_comparison_refuses_to_mix_conventions_without_reference_energies(
    tmp_path, structure_file
):
    """The ASE/OpenMM trap, caught by the comparison layer itself."""
    from dataclasses import replace

    from nio_md_prep.mlip.errors import EnergyConventionError

    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "run")["result"]
    as_interaction = replace(result, energy_convention="interaction", engine="openmm")
    with pytest.raises(EnergyConventionError, match="refusing to compare"):
        compare_results(result, as_interaction)


def test_comparison_harmonises_conventions_when_e0s_are_known(tmp_path, structure_file):
    from dataclasses import replace

    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "run")["result"]
    e0 = {"Ni": -5.78, "O": -4.95}
    offset = sum(e0[s] for s in result.symbols)
    shifted = replace(
        result,
        energy_eV=result.energy_eV - offset,
        energy_convention="interaction",
        engine="openmm",
        per_atom_energy_eV=None,
    )
    comparison = compare_results(result, shifted, atomic_reference_energies=e0)
    # After harmonisation the two are the same number.
    assert comparison.delta_energy_eV == pytest.approx(0.0, abs=1e-8)
    assert comparison.force_rmse_eV_per_A == pytest.approx(0.0, abs=1e-12)
    assert any("converted to" in note for note in comparison.notes)


def test_comparison_of_different_structures_is_refused(tmp_path, structure_file):
    from dataclasses import replace

    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "run")["result"]
    truncated = replace(
        result,
        symbols=result.symbols[:-1],
        forces_eV_per_A=result.forces_eV_per_A[:-1],
        per_atom_energy_eV=None,
    )
    with pytest.raises(ValueError, match="different structures"):
        compare_results(result, truncated)


def test_comparison_notes_a_missing_virial_instead_of_inventing_one(
    tmp_path, structure_file
):
    from dataclasses import replace

    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "run")["result"]
    stressless = replace(result, stress_eV_per_A3=None, engine="openmm")
    comparison = compare_results(result, stressless)
    assert comparison.stress_max_abs_error_eV_per_A3 is None
    assert any("did not report a virial" in note for note in comparison.notes)


def test_a_cross_engine_comparison_needs_at_least_two_engines(tmp_path, structure_file):
    job = parse_job(write_config(tmp_path, MOCK_CONFIG, structure_file))
    with pytest.raises(ConfigError, match="at least two engines"):
        jobs.compare_engines(job, ["ase"], output_dir=tmp_path / "cmp")


# --- CLI ------------------------------------------------------------------


def test_the_cli_dispatches_the_mlip_group(tmp_path, structure_file, capsys):
    from nio_md_prep.cli import main

    config = write_config(tmp_path, MOCK_CONFIG, structure_file)
    assert main(["mlip", "validate", str(config)]) == 0
    assert "mock-ase" in capsys.readouterr().out

    assert main(["mlip", "singlepoint", str(config), "--output", str(tmp_path / "cli")]) == 0
    assert "Single point via ase" in capsys.readouterr().out
    assert (tmp_path / "cli" / MANIFEST_NAME).exists()


def test_the_cli_reports_an_impossible_combination_as_an_error(tmp_path, capsys):
    """potential = lammps, engine = openmm must not produce a traceback."""
    from nio_md_prep.cli import main

    config = tmp_path / "impossible.toml"
    config.write_text(
        textwrap.dedent(
            """
            [potential]
            kind = "lammps"
            pair_style = "deepmd model.pb"
            pair_coeff = ["* *"]

            [potential.type_map]
            1 = "Ni"

            [engine]
            kind = "openmm"

            [simulation]
            task = "singlepoint"
            """
        ),
        encoding="utf-8",
    )
    with pytest.raises(SystemExit) as excinfo:
        main(["mlip", "validate", str(config)])
    assert excinfo.value.code == 2
    message = capsys.readouterr().err
    assert "error:" in message
    assert "compiled C++ inside LAMMPS" in message
    assert "Traceback" not in message
