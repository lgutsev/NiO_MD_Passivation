"""The MACE/ASE route against a real MACE calculator.

Two kinds of test, both marked ``mace``:

* **Self-contained** (need only mace-torch + torch): tiny *random* MACE
  models are built in the test (float64, float32 and a two-head model; a
  few kB each, well under a second to build). Random weights make them
  useless as potentials and ideal as wiring probes: every expectation here
  is an identity MACE itself defines -- ``total = sum(E0) + interaction``,
  the head's own E0 row, the dtype the parameters end up in, a periodic
  axis switched off changes the energy exactly as a vacuum wider than
  ``r_max`` does -- not a number copied from a previous run.
* **Trained model** (``NIO_MD_TEST_MACE_MODEL``, skipped without it): a
  single point and a 20-step NVE run on a NiO(110) slab
  (``NIO_MD_TEST_NIO_SLAB`` names a structure file; otherwise one is built
  with ASE), with finite energies/forces and a completed manifest.
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import replace
from pathlib import Path

import pytest

from nio_md_prep.mlip import jobs
from nio_md_prep.mlip.errors import CapabilityError, ConfigError, MissingDependencyError
from nio_md_prep.mlip.specs import (
    EngineSpec,
    JobSpec,
    MacePotentialSpec,
    SimulationSpec,
    StructureSpec,
)

pytestmark = pytest.mark.mace

#: Arbitrary per-element offsets baked into the random models (eV).
E0 = {"O": -4.95, "Ni": -5.78}
#: The second head's offsets in the two-head model.
E0_SECOND_HEAD = {"O": -1.0, "Ni": -2.0}
R_MAX = 4.5


def _build_model(torch, dtype, heads=None):
    import numpy as np
    from e3nn import o3
    from mace import modules, tools

    previous = torch.get_default_dtype()
    torch.set_default_dtype(dtype)
    try:
        torch.manual_seed(0)
        n_heads = len(heads) if heads else 1
        energies = (
            np.array([[E0["O"], E0["Ni"]], [E0_SECOND_HEAD["O"], E0_SECOND_HEAD["Ni"]]])
            if heads
            else np.array([E0["O"], E0["Ni"]])
        )
        config = dict(
            r_max=R_MAX, num_bessel=4, num_polynomial_cutoff=5, max_ell=1,
            interaction_cls=modules.interaction_classes["RealAgnosticResidualInteractionBlock"],
            interaction_cls_first=modules.interaction_classes[
                "RealAgnosticResidualInteractionBlock"],
            num_interactions=1, num_elements=2, hidden_irreps=o3.Irreps("8x0e"),
            MLP_irreps=o3.Irreps("4x0e"), gate=torch.nn.functional.silu,
            atomic_energies=energies, avg_num_neighbors=8.0,
            atomic_numbers=tools.AtomicNumberTable([8, 28]).zs, correlation=2,
            radial_type="bessel", atomic_inter_scale=[1.0] * n_heads,
            atomic_inter_shift=[0.0] * n_heads,
        )
        if heads:
            config["heads"] = heads
        return modules.ScaleShiftMACE(**config)
    finally:
        torch.set_default_dtype(previous)


@pytest.fixture(scope="module")
def tiny_models(require_mace, tmp_path_factory):
    """Random ScaleShiftMACE models on disk: ``float64``, ``float32``, ``two_head``."""
    import mace  # noqa: F401 - sets TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD
    import torch

    root = tmp_path_factory.mktemp("tiny-mace")
    paths = {}
    for name, dtype, heads in (
        ("float64", torch.float64, None),
        ("float32", torch.float32, None),
        ("two_head", torch.float64, ["pt", "nio"]),
    ):
        paths[name] = root / f"tiny_{name}.model"
        torch.save(_build_model(torch, dtype, heads), paths[name])
    return paths


@pytest.fixture(autouse=True)
def _fresh_discovery_cache():
    from nio_md_prep.mlip.potentials.mace import clear_discovery_cache

    clear_discovery_cache()
    yield
    clear_discovery_cache()


def nio_file(tmp_path, atoms, name="nio.xyz") -> Path:
    from ase.io import write

    path = tmp_path / name
    write(str(path), atoms, format="extxyz")
    return path


def reread(path: Path):
    """The structure exactly as the job reads it (extxyz keeps 8 decimals)."""
    from ase.io import read

    return read(str(path))


def job_for(model: Path, structure: Path, *, simulation=None, engine=None, **potential):
    return JobSpec(
        potential=MacePotentialSpec(
            label=model.stem, model_path=model, declared_elements=("Ni", "O"), **potential
        ),
        engine=engine or EngineSpec(kind="ase"),
        simulation=simulation or SimulationSpec(task="singlepoint", compute_stress=True,
                                                compute_per_atom_energy=True),
        structure=StructureSpec(path=structure),
    )


def direct_mace(model: Path, atoms, **kwargs):
    """The same geometry through a bare MACECalculator, for reference."""
    from mace.calculators import MACECalculator

    atoms = atoms.copy()
    atoms.calc = MACECalculator(model_paths=[str(model)], device="cpu",
                                default_dtype="float64", **kwargs)
    return atoms.get_potential_energy(), atoms.get_forces(apply_constraint=False), atoms.calc


# --- energy convention --------------------------------------------------------


def test_the_route_reports_macecalculators_total_energy(tmp_path, tiny_models,
                                                        rattled_nio_structure):
    np = pytest.importorskip("numpy")
    model = tiny_models["float64"]
    structure = nio_file(tmp_path, rattled_nio_structure)
    report = jobs.run_singlepoint(job_for(model, structure), output_dir=tmp_path / "sp")
    result = report["result"]
    energy, forces, _ = direct_mace(model, reread(structure))
    assert result.energy_convention == "total"
    assert result.energy_eV == pytest.approx(energy, abs=1e-10)
    np.testing.assert_allclose(result.forces_eV_per_A, forces, atol=1e-10)
    # MACE's per-atom energies include E0 and sum to the total.
    assert sum(result.per_atom_energy_eV) == pytest.approx(result.energy_eV, abs=1e-9)
    runtime = result.extras["runtime"]
    assert runtime["head"] == "Default" and runtime["device_effective"] == "cpu"
    assert runtime["dtype_effective"] == runtime["dtype_requested"] == "float64"
    assert runtime["kwargs_passed"]["default_dtype"] == "float64"
    manifest = json.loads((tmp_path / "sp" / "mlip_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["energy_convention"]["reported"] == "total"
    parameters = manifest["engine_parameters"]
    assert parameters["native_energy_convention"] == "total"
    assert parameters["atomic_reference_energies_eV"] == pytest.approx(E0)
    assert parameters["model_sha256"] == manifest["potential"]["model_sha256_observed"][str(model)]


@pytest.mark.parametrize("where", ["potential", "simulation"])
def test_an_interaction_request_subtracts_the_models_own_e0s(tmp_path, tiny_models,
                                                             rattled_nio_structure, where):
    model = tiny_models["float64"]
    structure = nio_file(tmp_path, rattled_nio_structure)
    total = jobs.run_singlepoint(job_for(model, structure), output_dir=tmp_path / "t")["result"]
    if where == "potential":
        job = job_for(model, structure, energy_convention="interaction")
    else:
        job = job_for(model, structure, simulation=SimulationSpec(
            task="singlepoint", compute_per_atom_energy=True, energy_convention="interaction"))
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "i")["result"]
    e0_sum = sum(E0[s] for s in rattled_nio_structure.get_chemical_symbols())
    assert result.energy_convention == "interaction"
    assert result.energy_eV == pytest.approx(total.energy_eV - e0_sum, abs=1e-9)
    assert result.extras["native_energy_convention"] == "total"
    assert result.extras["native_energy_eV"] == pytest.approx(total.energy_eV, abs=1e-12)
    assert result.extras["atomic_reference_energies_source"].startswith("model file")
    # The per-atom energies lose their own E0, and still sum to the reported energy.
    assert sum(result.per_atom_energy_eV) == pytest.approx(result.energy_eV, abs=1e-9)
    # MACE's own interaction split (node_energy) agrees.
    _, _, calculator = direct_mace(model, reread(structure))
    assert float(calculator.results["node_energy"].sum()) == pytest.approx(result.energy_eV,
                                                                           abs=1e-9)


def test_declared_e0s_are_cross_checked_against_the_model(tmp_path, tiny_models, nio_structure):
    structure = nio_file(tmp_path, nio_structure)
    model = tiny_models["float64"]
    agreeing = job_for(model, structure, atomic_reference_energies=dict(E0))
    assert jobs.validate_job(agreeing)["ok"]
    off = dict(E0, Ni=E0["Ni"] + 1e-6)
    with pytest.raises(ConfigError, match="model's own E0"):
        jobs.validate_job(job_for(model, structure, atomic_reference_energies=off))


def test_interaction_without_readable_model_e0s_is_refused(tmp_path, nio_structure):
    """Declared E0s do not stand in for the model's: no model file, no interaction."""
    job = job_for(tmp_path / "absent.model", nio_file(tmp_path, nio_structure),
                  energy_convention="interaction", atomic_reference_energies=dict(E0))
    with pytest.raises(CapabilityError, match="model's own E0s"):
        jobs.validate_job(job)


def test_declared_elements_the_model_lacks_are_refused(tmp_path, tiny_models, nio_structure):
    job = job_for(tiny_models["float64"], nio_file(tmp_path, nio_structure))
    job = replace(job, potential=replace(job.potential, declared_elements=("Ni", "O", "P")))
    with pytest.raises(ConfigError, match="declares P"):
        jobs.validate_job(job)


# --- heads, dtype, device, compile mode ----------------------------------------


def test_a_multi_head_model_needs_a_valid_head(tmp_path, tiny_models, rattled_nio_structure):
    model = tiny_models["two_head"]
    structure = nio_file(tmp_path, rattled_nio_structure)
    with pytest.raises(ConfigError, match="multi-head.*pt, nio"):
        jobs.validate_job(job_for(model, structure))
    # MACE itself would warn and evaluate the last head.
    with pytest.raises(ConfigError, match="not a head of this model"):
        jobs.validate_job(job_for(model, structure, head="nope"))

    results = {}
    for head in ("pt", "nio"):
        job = job_for(model, structure, head=head, simulation=SimulationSpec(
            task="singlepoint", energy_convention="interaction"))
        results[head] = jobs.run_singlepoint(job, output_dir=tmp_path / head)["result"]
        assert results[head].extras["runtime"]["head"] == head
    symbols = rattled_nio_structure.get_chemical_symbols()
    # Each head's interaction energy uses that head's own E0 row.
    for head, e0s in (("pt", E0), ("nio", E0_SECOND_HEAD)):
        total, _, calculator = direct_mace(model, reread(structure), head=head)
        assert calculator.head == head
        assert results[head].energy_eV == pytest.approx(
            total - sum(e0s[s] for s in symbols), abs=1e-9)
        assert results[head].extras["native_energy_eV"] == pytest.approx(total, abs=1e-10)


def test_a_float32_model_is_evaluated_in_the_requested_float64_and_says_so(
    tmp_path, tiny_models, rattled_nio_structure
):
    model = tiny_models["float32"]
    structure = nio_file(tmp_path, rattled_nio_structure)
    job = job_for(model, structure)
    plan = jobs.build_bridge(job.potential, job.engine).execution_plan(job.simulation)
    assert plan["precision"]["model_native"] == "float32"
    assert plan["precision"]["effective"] == "float64" and plan["precision"]["guaranteed"]
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "sp")["result"]
    runtime = result.extras["runtime"]
    assert runtime["dtype_model_native"] == "float32"
    assert runtime["dtype_effective"] == runtime["calculator_default_dtype"] == "float64"
    assert "converted" in runtime["dtype_note"]
    energy, _, _ = direct_mace(model, reread(structure))
    assert result.energy_eV == pytest.approx(energy, abs=1e-10)


def test_float32_evaluation_is_applied_when_requested(tmp_path, tiny_models, nio_structure):
    job = job_for(tiny_models["float64"], nio_file(tmp_path, nio_structure),
                  precision="float32")
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "sp")["result"]
    assert result.extras["runtime"]["dtype_effective"] == "float32"
    assert result.extras["runtime"]["dtype_model_native"] == "float64"


def test_an_unusable_cuda_device_is_caught_before_torch_deserialises(
    tmp_path, tiny_models, nio_structure
):
    import torch

    if torch.cuda.is_available():
        pytest.skip("this host has a CUDA device; the refusal path needs one without")
    job = job_for(tiny_models["float64"], nio_file(tmp_path, nio_structure), device="cuda")
    report = jobs.validate_job(job)
    assert report["availability"]["available"] is False
    assert "torch.cuda.is_available() is False" in report["availability"]["detail"]
    with pytest.raises(MissingDependencyError, match="no CUDA device"):
        jobs.run_singlepoint(job, output_dir=tmp_path / "gpu")
    assert not (tmp_path / "gpu" / "mlip_manifest.json").exists()


def test_compile_mode_is_validated(tmp_path, tiny_models, nio_structure):
    job = job_for(tiny_models["float64"], nio_file(tmp_path, nio_structure),
                  compile_mode="fastest-please")
    with pytest.raises(ConfigError, match="not a torch.compile mode"):
        jobs.validate_job(job)
    plan_job = replace(job, potential=replace(job.potential, compile_mode="default"))
    bridge = jobs.build_bridge(plan_job.potential, plan_job.engine)
    assert bridge.calculator_kwargs()["compile_mode"] == "default"


def test_engine_threads_set_torchs_thread_count(tmp_path, tiny_models, nio_structure):
    import torch

    before = torch.get_num_threads()
    try:
        job = job_for(tiny_models["float64"], nio_file(tmp_path, nio_structure),
                      engine=EngineSpec(kind="ase", threads=1))
        runtime = jobs.run_singlepoint(job, output_dir=tmp_path / "sp")["result"].extras["runtime"]
        assert runtime["torch_threads"] == 1
        assert runtime["torch_threads_source"].startswith("engine.threads")
    finally:
        torch.set_num_threads(before)


# --- one model load, periodicity, dynamics ------------------------------------------


def test_the_model_file_is_loaded_once_for_validation(tmp_path, tiny_models, nio_structure,
                                                      monkeypatch):
    import torch

    loads = []
    original = torch.load

    def counting(*args, **kwargs):
        loads.append(args[0] if args else kwargs.get("f"))
        return original(*args, **kwargs)

    monkeypatch.setattr(torch, "load", counting)
    job = job_for(tiny_models["float64"], nio_file(tmp_path, nio_structure))
    jobs.validate_job(job)
    assert len(loads) == 1  # discovery, cached by SHA256 for everything after
    jobs.run_singlepoint(job, output_dir=tmp_path / "sp")
    assert len(loads) == 2  # plus the calculator's own load
    md = replace(job, simulation=SimulationSpec(task="md", ensemble="nve", timestep_fs=1.0,
                                                steps=4, seed=3, temperature_K=100.0))
    jobs.run_smoke_md(md, output_dir=tmp_path / "md")
    assert len(loads) == 3


def test_a_non_periodic_axis_is_honoured(tmp_path, tiny_models):
    """pbc (T, T, F) on a gap-free cell equals pbc (T, T, T) with vacuum wider than r_max."""
    from ase.build import bulk

    model = tiny_models["float64"]
    dense = bulk("NiO", "rocksalt", a=4.17, cubic=True).repeat((2, 2, 1))
    dense.rattle(0.03, seed=5)
    slab = dense.copy()
    slab.pbc = (True, True, False)
    padded = dense.copy()
    padded.center(vacuum=2 * R_MAX, axis=2)  # pbc stays (T, T, T)
    energies = {}
    for name, atoms in (("dense", dense), ("slab", slab), ("padded", padded)):
        job = job_for(model, nio_file(tmp_path, atoms, f"{name}.xyz"),
                      simulation=SimulationSpec(task="singlepoint"))
        energies[name] = jobs.run_singlepoint(job, output_dir=tmp_path / name)["result"].energy_eV
    assert energies["slab"] == pytest.approx(energies["padded"], abs=1e-9)
    assert abs(energies["slab"] - energies["dense"]) > 1e-3


def test_a_short_nve_run_conserves_energy(tmp_path, tiny_models, rattled_nio_structure):
    job = job_for(tiny_models["float64"], nio_file(tmp_path, rattled_nio_structure),
                  simulation=SimulationSpec(task="md", ensemble="nve", timestep_fs=0.5,
                                            steps=20, seed=11, temperature_K=300.0,
                                            trajectory_interval=5, log_interval=1))
    report = jobs.run_smoke_md(job, output_dir=tmp_path / "md")
    trajectory = report["trajectory"]
    assert trajectory.steps_completed == 20 and trajectory.frames_written == 5
    assert trajectory.integrator_resolved["integrator"] == "VelocityVerlet"
    diagnostics = trajectory.diagnostics
    assert diagnostics["conservation_test"] is True
    assert math.isfinite(diagnostics["energy_drift_eV_per_atom_per_ps"])
    # 10 fs of velocity Verlet on a smooth random model: the excursion is tiny.
    assert diagnostics["max_abs_energy_excursion_eV_per_atom"] < 1e-3
    assert trajectory.extras["runtime"]["dtype_effective"] == "float64"
    manifest = json.loads((tmp_path / "md" / "mlip_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["results"]["trajectory"]["integrator_resolved"]["class"] == (
        "ase.md.verlet.VelocityVerlet")


# --- the execution plan -------------------------------------------------------------

PLAN_KEYS = {
    "potential_kind", "engine", "bridge", "implementation", "model_checkpoint",
    "exported_model", "elements", "energy_convention", "units", "device", "precision",
    "dynamics", "lammps", "openmm", "availability", "unmet_capabilities", "model_hashes",
}


def test_the_execution_plan_states_what_will_run(tmp_path, tiny_models, rattled_nio_structure):
    from ase import units

    model = tiny_models["float64"]
    job = job_for(model, nio_file(tmp_path, rattled_nio_structure), simulation=SimulationSpec(
        task="md", ensemble="nvt", thermostat="csvr", temperature_K=300.0, timestep_fs=2.0,
        steps=10, thermostat_damping_fs=50.0))
    bridge = jobs.build_bridge(job.potential, job.engine)
    plan = bridge.execution_plan(job.simulation, rattled_nio_structure)
    assert PLAN_KEYS <= set(plan)
    from nio_md_prep.mlip.potentials.mace import model_sha256

    sha = model_sha256(model)
    assert plan["model_checkpoint"] == {"path": str(model), "sha256": sha}
    assert plan["model_hashes"]["observed"] == {str(model): sha}
    assert plan["exported_model"] is None and plan["lammps"] is None and plan["openmm"] is None
    assert plan["implementation"] == "mace-ase-calculator" and plan["engine"] == "ase"
    assert plan["elements"] == ["Ni", "O"] and plan["energy_convention"] == "total"
    assert plan["units"]["native"] == "ase" and plan["units"]["pressure_unit"] == "eV/Angstrom^3"
    assert plan["device"] == {"requested": "cpu", "effective": "cpu", "guaranteed": True,
                              "note": plan["device"]["note"]}
    dynamics = plan["dynamics"]
    assert dynamics["integrator"] == "ase.md.bussi.Bussi" and dynamics["thermostat"] == "csvr"
    assert dynamics["timestep_native"] == pytest.approx(2.0 * units.fs, rel=1e-15)
    assert dynamics["thermostat_damping_fs"] == 50.0
    assert dynamics["thermostat_damping_native"] == pytest.approx(50.0 * units.fs, rel=1e-15)
    assert dynamics["temperature_ndof"] == 3 * len(rattled_nio_structure) - 3
    assert plan["availability"]["available"] and plan["unmet_capabilities"] == []
    assert plan["calculator"]["kwargs"]["head"] == "Default"


def test_the_execution_plan_lists_refusals_instead_of_raising(tmp_path, tiny_models,
                                                              nio_structure):
    job = job_for(tiny_models["two_head"], nio_file(tmp_path, nio_structure),
                  simulation=SimulationSpec(task="singlepoint", compute_per_atom_energy=True))
    plan = jobs.build_bridge(job.potential, job.engine).execution_plan(job.simulation,
                                                                       nio_structure)
    assert any("multi-head" in item for item in plan["unmet_capabilities"])
    assert plan["head"] is None


# --- a trained model on a NiO(110) slab (NIO_MD_TEST_MACE_MODEL) ------------------


@pytest.fixture
def nio_110_slab_path(tmp_path):
    """``NIO_MD_TEST_NIO_SLAB`` if set, else a 4-layer (2x2) NiO(110) slab built with ASE."""
    raw = os.environ.get("NIO_MD_TEST_NIO_SLAB")
    if raw:
        path = Path(raw)
        if not path.exists():
            pytest.skip(f"NIO_MD_TEST_NIO_SLAB points at {path}, which does not exist")
        return path
    from ase.build import bulk, surface

    slab = surface(bulk("NiO", "rocksalt", a=4.17), (1, 1, 0), layers=4, vacuum=8.0)
    slab = slab.repeat((2, 2, 1))
    return nio_file(tmp_path, slab, "nio110.xyz")


def test_a_trained_model_evaluates_a_nio_110_slab(require_mace, mace_model_path,
                                                  nio_110_slab_path, tmp_path):
    np = pytest.importorskip("numpy")
    from ase.io import read

    atoms = read(str(nio_110_slab_path))
    elements = tuple(sorted(set(atoms.get_chemical_symbols())))
    job = JobSpec(
        potential=MacePotentialSpec(label=mace_model_path.stem, model_path=mace_model_path,
                                    declared_elements=elements, precision="float64"),
        engine=EngineSpec(kind="ase"),
        # A stress tensor needs a cell periodic along all three axes: an ASE-built
        # vacuum slab keeps pbc (T, T, F); a POSCAR reads back fully periodic.
        simulation=SimulationSpec(task="singlepoint", compute_stress=bool(atoms.pbc.all())),
        structure=StructureSpec(path=nio_110_slab_path),
    )
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "sp")["result"]
    assert result.n_atoms == len(atoms) and math.isfinite(result.energy_eV)
    assert np.isfinite(np.array(result.forces_eV_per_A)).all()
    assert result.energy_convention == "total"
    manifest = json.loads((tmp_path / "sp" / "mlip_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["results"]["singlepoint"]["n_atoms"] == len(atoms)
    md = replace(job, simulation=SimulationSpec(task="md", ensemble="nve", timestep_fs=1.0,
                                                steps=20, seed=5, temperature_K=100.0,
                                                trajectory_interval=5, log_interval=1))
    trajectory = jobs.run_smoke_md(md, output_dir=tmp_path / "md")["trajectory"]
    assert trajectory.steps_completed == 20 and trajectory.frames_written == 5
    assert math.isfinite(trajectory.diagnostics["energy_drift_eV_per_atom_per_ps"])
