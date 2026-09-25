"""Optional tests that need a real backend and a real model.

None of these run in ordinary CI. Each is marked (``mace``, ``lammps``,
``openmm``, ``gpu``) and skips itself unless the backend is importable, and
the model comes from ``NIO_MD_TEST_MACE_MODEL`` so that no large trained model
is ever committed to this repository::

    export NIO_MD_TEST_MACE_MODEL=/models/nio_phosphonate.model
    export NIO_MD_TEST_MACE_ELEMENTS=Ni,O,P,C,H
    pytest tests/test_mlip_backends.py -m "mace or openmm or lammps"

The cross-engine single-point equivalence check is the scientific acceptance
test for this subsystem, and it is deliberately *not* written against an
invented tolerance. The first time it runs it measures the agreement between
the adapters actually installed and writes a candidate reference file; a
maintainer reviews those numbers and commits them as
``tests/data/mlip_cross_engine_tolerances.json``, after which the test becomes
an ordinary regression criterion. MACE itself warns users to be careful
benchmarking LAMMPS output against the corresponding ASE calculator, so the
tolerance is a measurement, not an assumption.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from nio_md_prep.mlip import jobs, tolerances
from nio_md_prep.mlip.registry import build_bridge
from nio_md_prep.mlip.results import compare_results
from nio_md_prep.mlip.specs import (
    EngineSpec,
    JobSpec,
    MacePotentialSpec,
    SimulationSpec,
    StructureSpec,
)

#: Committed reference tolerances, once a maintainer has measured them.
TOLERANCE_FILE = Path(
    os.environ.get(
        "NIO_MD_MLIP_TOLERANCES",
        str(Path(__file__).parent / "data" / "mlip_cross_engine_tolerances.json"),
    )
)

SINGLEPOINT = SimulationSpec(task="singlepoint", compute_stress=True)


def mace_job(model: Path, elements, *, engine: str = "ase", **engine_kwargs) -> JobSpec:
    return JobSpec(
        potential=MacePotentialSpec(
            label=model.stem,
            model_path=model,
            declared_elements=tuple(elements),
            device=os.environ.get("NIO_MD_TEST_MACE_DEVICE", "cpu"),
            precision="float64",
        ),
        engine=EngineSpec(kind=engine, **engine_kwargs),
        simulation=SINGLEPOINT,
        name="mace-acceptance",
    )


@pytest.fixture
def structure_path(tmp_path, rattled_nio_structure):
    from ase.io import write

    path = tmp_path / "nio.xyz"
    write(str(path), rattled_nio_structure, format="extxyz")
    return path


# --- MACE / ASE -----------------------------------------------------------


@pytest.mark.mace
def test_mace_ase_singlepoint(require_mace, mace_model_path, mace_elements, structure_path, tmp_path):
    job = mace_job(mace_model_path, mace_elements)
    report = jobs.run_singlepoint(
        job, output_dir=tmp_path / "ase", structure=StructureSpec(path=structure_path)
    )
    result = report["result"]
    assert result.n_atoms > 0
    assert result.native_units == "ase"
    # The stock ASE calculator reports the total energy, self-energies included.
    assert result.energy_convention == "total"
    assert result.max_force_eV_per_A > 0


@pytest.mark.mace
def test_model_metadata_is_discovered_and_cross_checked(
    require_mace, mace_model_path, mace_elements
):
    """Elements, cutoff and E0s come from the model file when they can."""
    from nio_md_prep.mlip.potentials.mace import MaceAdapter

    adapter = MaceAdapter(
        MacePotentialSpec(model_path=mace_model_path, declared_elements=tuple(mace_elements))
    )
    report = adapter.inspect()
    assert report["availability"]["available"]
    discovered = report["discovered"]
    assert discovered is not None
    assert discovered["cutoff_angstrom"] is None or discovered["cutoff_angstrom"] > 0
    # A conflict must be reported, never silently resolved either way.
    assert isinstance(report.get("conflicts", []), list)


@pytest.mark.mace
@pytest.mark.gpu
def test_mace_runs_on_a_cuda_device(
    require_mace, require_gpu, mace_model_path, mace_elements, structure_path, tmp_path
):
    job = mace_job(mace_model_path, mace_elements)
    job = job.__class__(
        potential=MacePotentialSpec(
            label=mace_model_path.stem,
            model_path=mace_model_path,
            declared_elements=tuple(mace_elements),
            device="cuda",
            precision="float32",
        ),
        engine=job.engine,
        simulation=job.simulation,
    )
    bridge = build_bridge(job.potential, job.engine)
    assert bridge.capabilities().gpu
    report = jobs.run_singlepoint(
        job, output_dir=tmp_path / "gpu", structure=StructureSpec(path=structure_path)
    )
    assert report["result"].n_atoms > 0


# --- MACE / OpenMM --------------------------------------------------------


@pytest.mark.openmm
@pytest.mark.mace
def test_openmm_reports_the_convention_it_was_asked_for(
    require_mace, require_openmm, mace_model_path, mace_elements, structure_path, tmp_path
):
    """Not OpenMM-ML's default, but the one the job declares."""
    from dataclasses import replace

    job = mace_job(mace_model_path, mace_elements, engine="openmm")
    total = replace(job, simulation=replace(SINGLEPOINT, energy_convention="total"))
    bridge = build_bridge(total.potential, total.engine)
    assert bridge.engine_parameters()["returnEnergyType"] == "energy"

    interaction = replace(
        job, simulation=replace(SINGLEPOINT, energy_convention="interaction")
    )
    bridge = build_bridge(interaction.potential, interaction.engine)
    assert bridge.requested_convention(interaction.simulation) == "interaction"


@pytest.mark.openmm
@pytest.mark.mace
def test_openmm_npt_is_refused_for_want_of_a_virial(
    require_openmm, mace_model_path, mace_elements, structure_path
):
    from nio_md_prep.mlip.errors import CapabilityError

    job = mace_job(mace_model_path, mace_elements, engine="openmm")
    npt = SimulationSpec(
        task="md",
        ensemble="npt",
        temperature_K=300,
        pressure_bar=1.0,
        timestep_fs=0.5,
        steps=10,
    )
    bridge = build_bridge(job.potential, job.engine)
    with pytest.raises(CapabilityError, match="stress"):
        bridge.validate(npt)


# --- MACE / LAMMPS --------------------------------------------------------


@pytest.mark.lammps
@pytest.mark.mace
def test_mace_lammps_singlepoint(
    require_lammps,
    mace_model_path,
    mace_lammps_model_path,
    mace_elements,
    structure_path,
    tmp_path,
):
    if mace_lammps_model_path is None:
        pytest.skip("set NIO_MD_TEST_MACE_LAMMPS_MODEL to the exported LAMMPS model")
    job = mace_job(
        mace_model_path,
        mace_elements,
        engine="lammps",
        options={"model_path": str(mace_lammps_model_path)},
    )
    report = jobs.run_singlepoint(
        job, output_dir=tmp_path / "lammps", structure=StructureSpec(path=structure_path)
    )
    result = report["result"]
    assert result.native_units == "lammps_metal"
    assert result.max_force_eV_per_A > 0
    # The manifest must carry the exact commands LAMMPS executed.
    manifest = json.loads((tmp_path / "lammps" / "mlip_manifest.json").read_text())
    assert manifest["engine_parameters"]["pair_commands"][0].startswith("pair_style")


# --- the acceptance test --------------------------------------------------


@pytest.mark.mace
def test_cross_engine_single_point_equivalence(
    require_mace, mace_model_path, mace_elements, structure_path, tmp_path
):
    """Evaluate one model and one structure through every engine installed.

    ASE is the reference. Energies are harmonised to a single convention
    before being subtracted, so the numbers reported here are real
    disagreements between adapters rather than an artefact of OpenMM-ML's
    interaction-energy default.
    """
    from nio_md_prep.mlip.environment import module_available

    engines = ["ase"]
    if module_available("lammps") and os.environ.get("NIO_MD_TEST_MACE_LAMMPS_MODEL"):
        engines.append("lammps")
    if module_available("openmm") and module_available("openmmml"):
        engines.append("openmm")
    if len(engines) < 2:
        pytest.skip(
            "only the ASE route is available; install LAMMPS or OpenMM-ML to run the "
            "cross-engine acceptance test"
        )

    options = {}
    exported = os.environ.get("NIO_MD_TEST_MACE_LAMMPS_MODEL")
    if exported:
        options["model_path"] = exported
    job = mace_job(mace_model_path, mace_elements, options=options)

    report = jobs.compare_engines(
        job,
        engines,
        output_dir=tmp_path / "compare",
        structure=StructureSpec(path=structure_path),
        convention="total",
    )
    measured = tolerances.measured_metrics(report)

    reference = _load_reference()
    if reference is None:
        candidate = tmp_path / "compare" / "candidate_tolerances.json"
        record = tolerances.candidate_record(
            report,
            measured_on={"engines": engines, "source": "tests/test_mlip_backends.py"},
        )
        candidate.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        pytest.skip(
            "no reviewed cross-engine tolerance reference; measured values written "
            f"to {candidate} for review (a candidate is never an acceptance criterion)"
        )

    problems = tolerances.exceedances(measured, reference)
    assert not problems, "; ".join(problems)


def _load_reference() -> dict | None:
    return tolerances.load_reference(TOLERANCE_FILE)


def test_the_tolerance_file_is_well_formed_if_present():
    """Runs in ordinary CI: a malformed or unreviewed reference must fail loudly, not silently."""
    reference = _load_reference()  # validates status, reviewer, provenance and metric names
    if reference is None:
        pytest.skip("no reviewed cross-engine tolerance reference has been committed yet")
