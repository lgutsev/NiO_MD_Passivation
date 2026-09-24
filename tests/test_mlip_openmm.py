"""The OpenMM route (MACE through OpenMM-ML), without and with OpenMM.

The first half runs in ordinary CI: box reduction math, the platform-property
builder, the precision/device/dynamics plans, DCD frame counting, the static
model read and ``execution_plan`` -- none of it imports torch or OpenMM.

The second half is marked ``openmm`` and ``mace`` and needs both installed. It
builds a tiny random ScaleShiftMACE model in-test (no model file is committed)
and checks the route against MACE's own ASE calculator. With the isolated
environment used for this project::

    set MKL_THREADING_LAYER=SEQUENTIAL
    set TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1
    env-openmm/python.exe -m pytest tests/test_mlip_openmm.py -m "openmm or not openmm"

(``MKL_THREADING_LAYER`` must be set before numpy loads MKL in that
environment; the weights-only switch is what lets torch >= 2.6 load a pickled
MACE model, which OpenMM-ML does with a plain ``torch.load``.)
"""
from __future__ import annotations

import os
import struct
import zipfile
from pathlib import Path

import pytest

from nio_md_prep.mlip.errors import CapabilityError, ConfigError, ResultError
from nio_md_prep.mlip.specs import EngineSpec, MacePotentialSpec, SimulationSpec

np = pytest.importorskip("numpy")
pytest.importorskip("ase")

from nio_md_prep.mlip.engines import openmm_engine  # noqa: E402

#: The keys every bridge's execution_plan() returns (the shared contract).
PLAN_KEYS = {
    "potential_kind", "engine", "bridge", "implementation", "model_checkpoint",
    "exported_model", "elements", "energy_convention", "units", "device", "precision",
    "dynamics", "lammps", "openmm", "availability", "unmet_capabilities", "model_hashes",
}
DYNAMICS_KEYS = {
    "ensemble", "integrator", "thermostat", "barostat", "barostat_coupling", "timestep_fs",
    "timestep_native", "thermostat_damping_fs", "thermostat_damping_native",
    "barostat_damping_fs", "barostat_damping_native",
}
GUARANTEE_KEYS = {"requested", "effective", "guaranteed", "note"}


def nvt(**overrides) -> SimulationSpec:
    settings = dict(
        task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=1.0, steps=20, seed=1234,
        trajectory_interval=5, log_interval=5,
    )
    settings.update(overrides)
    return SimulationSpec(**settings)


def openmm_engine_spec(**kwargs) -> EngineSpec:
    return EngineSpec(kind="openmm", **kwargs)


# ---------------------------------------------------------------------------
# Box reduction (pure)
# ---------------------------------------------------------------------------


def _check_reduced(box, cell):
    a, b, c = box.vectors_nm
    assert openmm_engine.box_violations(a, b, c) == []
    rotation = np.asarray(box.rotation)
    assert np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-12)
    change = np.asarray(box.lattice_change)
    assert abs(round(np.linalg.det(change))) == 1
    # Same lattice: the reduced vectors are an integer basis change of the
    # rotated source cell, and the cell volume is unchanged.
    rotated = np.asarray(cell) @ rotation.T
    assert np.allclose(np.asarray(box.vectors_nm) * 10.0, change @ rotated, atol=1e-9)
    assert abs(np.linalg.det(np.asarray(box.vectors_nm) * 10.0)) == pytest.approx(
        abs(np.linalg.det(cell)), rel=1e-12
    )


def test_the_fcc_primitive_cell_is_rejected_raw_and_reduced_exactly(nio_structure):
    """``bulk('NiO', 'rocksalt')`` is the cell PR #19 handed OpenMM as-is."""
    cell = nio_structure.get_cell().array
    raw_nm = cell / 10.0
    assert openmm_engine.box_violations(*raw_nm) != []
    box = openmm_engine.reduced_box(cell)
    _check_reduced(box, cell)
    assert not np.allclose(box.rotation, np.eye(3))  # a real rotation is exercised
    assert not box.basis_flipped


@pytest.mark.parametrize("seed", range(12))
def test_random_triclinic_cells_reduce_to_openmm_form(seed):
    rng = np.random.default_rng(seed)
    cell = rng.normal(size=(3, 3)) * 4.0 + np.eye(3) * 6.0
    if seed % 3 == 0:
        cell[2] *= -1.0  # left-handed basis
    box = openmm_engine.reduced_box(cell)
    _check_reduced(box, cell)
    assert box.basis_flipped == (np.linalg.det(cell) < 0)


def test_rotation_round_trip_preserves_vectors_and_physics(rattled_nio_structure):
    box = openmm_engine.reduced_box(rattled_nio_structure.get_cell())
    positions = rattled_nio_structure.get_positions()
    engine = box.to_engine(positions)
    assert np.allclose(box.to_source(engine), positions, atol=1e-12)
    # Distances are invariant under the rotation.
    assert np.allclose(
        np.linalg.norm(engine[1:] - engine[0], axis=1),
        np.linalg.norm(positions[1:] - positions[0], axis=1),
        atol=1e-12,
    )
    # The rotated cell is lower triangular (what OpenMM's first checks demand).
    rotated_cell = box.to_engine(rattled_nio_structure.get_cell().array)
    assert abs(rotated_cell[0, 1]) < 1e-12 and abs(rotated_cell[0, 2]) < 1e-12
    assert abs(rotated_cell[1, 2]) < 1e-12


def test_a_degenerate_cell_is_refused():
    with pytest.raises(ConfigError, match="degenerate"):
        openmm_engine.reduced_box([[4.0, 0, 0], [8.0, 0, 0], [0, 0, 4.0]])


# ---------------------------------------------------------------------------
# Platform properties and plans (pure)
# ---------------------------------------------------------------------------


def test_precision_property_is_refused_on_the_cpu_platform():
    with pytest.raises(ConfigError, match="'Precision' property"):
        openmm_engine.platform_properties("CPU", platform_precision="double")
    with pytest.raises(ConfigError, match="only the CUDA"):
        EngineSpec(kind="openmm", platform="CPU", platform_precision="double")


def test_threads_map_to_the_cpu_platform_only():
    assert openmm_engine.platform_properties("CPU", threads=4) == {"Threads": "4"}
    with pytest.raises(ConfigError, match="Threads"):
        openmm_engine.platform_properties("OpenCL", threads=4)
    with pytest.raises(ConfigError, match="named platform"):
        openmm_engine.platform_properties(None, threads=4)
    assert openmm_engine.platform_properties("Reference") == {}


def test_the_torch_device_maps_to_device_index_only_where_ordinals_agree():
    props = openmm_engine.platform_properties
    assert props("CUDA", device="cuda:1") == {"DeviceIndex": "1"}
    assert props("CUDA", device="cuda") == {"DeviceIndex": "0"}  # torch's current device
    assert props("HIP", device="cuda:2") == {"DeviceIndex": "2"}
    assert props("OpenCL", device="cuda:1") == {}  # OpenCL numbers devices differently
    assert props("CUDA", device="cpu") == {}
    assert props("CUDA", platform_precision="mixed", device="cuda:0") == {
        "Precision": "mixed",
        "DeviceIndex": "0",
    }


@pytest.mark.parametrize(
    "version, stored, guaranteed, refused",
    [
        ("1.7", "float64", True, False),
        ("1.7", "float32", False, True),
        ("1.7", None, False, False),
        ("1.6.1", "float32", False, True),
        ("1.8", "float32", True, False),
        ("1.8.0", None, True, False),
        (None, "float64", False, False),
        ("1.5", "float64", False, False),
    ],
)
def test_precision_is_guaranteed_only_when_the_installed_openmmml_delivers_it(
    version, stored, guaranteed, refused
):
    """1.6/1.7 cast inputs but not the model; 1.8 converts the model."""
    plan = openmm_engine.precision_plan(
        "float64", openmmml=version, model_dtype=stored, model_dtype_source="test"
    )
    assert GUARANTEE_KEYS <= set(plan)
    assert plan["guaranteed"] is guaranteed
    assert plan["effective"] == ("float64" if guaranteed else None)
    assert ("refusal" in plan) is refused
    assert plan["note"]
    assert plan["passed_as"] == {"createSystem": {"precision": "double"}}


def test_device_plan_is_guaranteed_by_the_device_keyword_and_explains_the_platform():
    plan = openmm_engine.device_plan("cuda:1", platform="CUDA", openmmml="1.7")
    assert plan["guaranteed"] and plan["effective"] == "cuda:1"
    assert plan["openmm_platform_device"] == {"DeviceIndex": "1"}
    opencl = openmm_engine.device_plan("cuda:1", platform="OpenCL", openmmml="1.7")
    assert opencl["openmm_platform_device"] is None
    assert "OpenCL device numbers are not CUDA ordinals" in opencl["note"]
    mixed = openmm_engine.device_plan("cpu", platform="CUDA", openmmml="1.7")
    assert "cross between host and device" in mixed["note"]
    unknown = openmm_engine.device_plan("cpu", platform=None, openmmml=None)
    assert not unknown["guaranteed"] and unknown["effective"] is None


def test_platform_precision_plan_never_claims_what_openmm_decides_later():
    reference = openmm_engine.platform_precision_plan("Reference", None)
    assert reference["effective"] == "double" and reference["guaranteed"]
    cuda = openmm_engine.platform_precision_plan("CUDA", "mixed")
    assert cuda["effective"] == "mixed" and cuda["guaranteed"]
    default = openmm_engine.platform_precision_plan("OpenCL", None)
    assert default["effective"] == "single" and not default["guaranteed"]
    cpu = openmm_engine.platform_precision_plan("CPU", None)
    assert cpu["effective"] is None and not cpu["guaranteed"]
    chosen = openmm_engine.platform_precision_plan(None, None)
    assert chosen["effective"] is None and not chosen["guaranteed"]


def test_dynamics_plan_resolves_thermostat_damping_in_both_unit_systems():
    langevin = openmm_engine.dynamics_plan(nvt(timestep_fs=0.5))
    assert DYNAMICS_KEYS <= set(langevin)
    assert langevin["integrator"] == "LangevinMiddleIntegrator"
    assert langevin["thermostat"] == "langevin" and langevin["thermostat_defaulted"]
    assert langevin["thermostat_damping_fs"] == 100.0
    assert langevin["thermostat_damping_native"]["value"] == pytest.approx(10.0)  # 1/ps
    assert langevin["thermostat_damping_native"]["parameter"] == "frictionCoeff"
    assert langevin["timestep_native"] == {"value": 0.0005, "unit": "ps"}
    hoover = openmm_engine.dynamics_plan(nvt(thermostat="nose-hoover", thermostat_damping_fs=50.0))
    assert hoover["integrator"] == "NoseHooverIntegrator" and not hoover["thermostat_defaulted"]
    assert hoover["thermostat_damping_native"]["value"] == pytest.approx(20.0)
    assert hoover["thermostat_damping_native"]["parameter"] == "collisionFrequency"
    nve = openmm_engine.dynamics_plan(
        SimulationSpec(task="md", ensemble="nve", timestep_fs=1.0, steps=5)
    )
    assert nve["integrator"] == "VerletIntegrator" and nve["thermostat_damping_native"] is None
    assert openmm_engine.dynamics_plan(SimulationSpec()) is None


@pytest.mark.parametrize("thermostat", ["berendsen", "csvr"])
def test_unimplemented_thermostats_are_refused_not_substituted(thermostat):
    with pytest.raises(ConfigError, match="never substituted"):
        openmm_engine.check_md_request(nvt(thermostat=thermostat))
    assert openmm_engine.dynamics_plan(nvt(thermostat=thermostat))["integrator"] is None


def test_npt_is_refused_by_policy_with_an_accurate_rationale():
    npt = SimulationSpec(
        task="md", ensemble="npt", temperature_K=300, pressure_bar=1.0, timestep_fs=1.0, steps=5
    )
    with pytest.raises(CapabilityError, match="MonteCarloBarostat needs only energies"):
        openmm_engine.check_md_request(npt)


def test_partially_periodic_structures_are_refused(nio_structure):
    slab = nio_structure.copy()
    slab.pbc = (True, True, False)
    with pytest.raises(CapabilityError, match="all three axes or none"):
        openmm_engine.check_structure(slab, md=False)


def test_child_seeds_are_reproducible_distinct_and_never_zero():
    first = openmm_engine.derive_seeds(1234)
    assert first == openmm_engine.derive_seeds(1234)
    assert first["velocity_seed"] != first["integrator_seed"]
    for seed in (1, 2, 2**31 - 1):
        for value in openmm_engine.derive_seeds(seed).values():
            assert 1 <= value <= 2**31 - 1


def _fake_dcd(path: Path, n_atoms: int, frames: int, *, box: bool, truncate: int = 0) -> None:
    """Bytes laid out exactly as ``openmm.app.DCDFile`` writes them."""
    header = struct.pack("<i4c9if", 84, b"C", b"O", b"R", b"D", frames, 0, 1, 0, 0, 0, 0, 0, 0, 0.02)
    header += struct.pack("<13i", 1 if box else 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 84, 164, 2)
    header += struct.pack("<80s", b"Created by OpenMM") + struct.pack("<80s", b"Created now")
    header += struct.pack("<4i", 164, 4, n_atoms, 4)
    frame = b""
    if box:
        frame += struct.pack("<i6di", 48, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 48)
    for _ in range(3):
        frame += struct.pack("<i", 4 * n_atoms) + b"\0" * (4 * n_atoms) + struct.pack("<i", 4 * n_atoms)
    data = header + frame * frames
    path.write_bytes(data[: len(data) - truncate] if truncate else data)


@pytest.mark.parametrize("box", [True, False])
def test_dcd_frames_are_counted_from_the_file(tmp_path, box):
    path = tmp_path / "t.dcd"
    _fake_dcd(path, 5, 3, box=box)
    assert openmm_engine.dcd_frame_count(path, 5) == 3
    with pytest.raises(ResultError, match="describes 5 atoms"):
        openmm_engine.dcd_frame_count(path, 6)
    _fake_dcd(path, 5, 3, box=box, truncate=7)
    with pytest.raises(ResultError, match="truncated or corrupt"):
        openmm_engine.dcd_frame_count(path, 5)
    with pytest.raises(ResultError, match="was not written"):
        openmm_engine.dcd_frame_count(tmp_path / "missing.dcd", 5)


# ---------------------------------------------------------------------------
# The bridge without OpenMM: static model read, gating, execution_plan
# ---------------------------------------------------------------------------


def _fake_model(path: Path, *globals_: str) -> Path:
    """A zip laid out like ``torch.save`` output, naming the given pickle globals."""
    body = b"\x80\x02"
    for name in globals_:
        module, attribute = name.split(" ")
        body += b"c" + module.encode() + b"\n" + attribute.encode() + b"\n0"
    body += b"N."
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{path.stem}/data.pkl", body)
    return path


def _bridge(model_path, *, engine=None, **potential):
    from nio_md_prep.mlip.bridges.mace_openmm import MaceOpenMMBridge

    spec = MacePotentialSpec(
        model_path=model_path, declared_elements=("Ni", "O"), **potential
    )
    return MaceOpenMMBridge(spec, engine or openmm_engine_spec(platform="CPU"))


def test_the_stored_dtype_and_class_are_read_without_torch(tmp_path, backend_import_guard):
    from nio_md_prep.mlip.bridges.mace_openmm import _peek_model

    double = _fake_model(
        tmp_path / "m64.model", "mace.modules.models ScaleShiftMACE", "torch DoubleStorage",
        "torch LongStorage",
    )
    assert _peek_model(double)["dtype"] == "float64"
    assert _peek_model(double)["model_class"] == "ScaleShiftMACE"
    mixed = _fake_model(tmp_path / "mix.model", "torch FloatStorage", "torch DoubleStorage")
    assert _peek_model(mixed)["dtype"] is None and "several dtypes" in _peek_model(mixed)["error"]
    junk = tmp_path / "junk.model"
    junk.write_bytes(b"not a zip")
    assert _peek_model(junk)["dtype"] is None and "error" in _peek_model(junk)
    backend_import_guard()


def test_protocol_4_stack_globals_are_read_the_same_way(tmp_path):
    """Protocol >= 4 names globals with STACK_GLOBAL instead of GLOBAL."""
    import pickletools

    from nio_md_prep.mlip.bridges.mace_openmm import _peek_model

    body = (
        b"\x80\x04"
        + b"\x8c\x05torch\x94"
        + b"\x8c\x0dDoubleStorage\x94"
        + b"\x93\x94"
        + b"0N."
    )
    assert [op.name for op, _, _ in pickletools.genops(body)].count("STACK_GLOBAL") == 1
    path = tmp_path / "p4.model"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("p4/data.pkl", body)
    assert _peek_model(path)["dtype"] == "float64"


def test_validate_refuses_head_compile_mode_thermostat_and_slabs(tmp_path, nio_structure):
    model = _fake_model(tmp_path / "m.model", "mace.modules.models ScaleShiftMACE", "torch DoubleStorage")
    with pytest.raises(CapabilityError, match="head"):
        _bridge(model, head="Default").validate(SimulationSpec())
    with pytest.raises(ConfigError, match="compile_mode"):
        _bridge(model, compile_mode="default").validate(SimulationSpec())
    with pytest.raises(ConfigError, match="berendsen"):
        _bridge(model).validate(nvt(thermostat="berendsen"))
    slab = nio_structure.copy()
    slab.pbc = (True, True, False)
    with pytest.raises((CapabilityError, ConfigError)):
        _bridge(model).validate(SimulationSpec(), slab)


def test_validate_refuses_a_precision_this_openmmml_would_not_deliver(tmp_path, monkeypatch):
    model = _fake_model(tmp_path / "m32.model", "mace.modules.models ScaleShiftMACE", "torch FloatStorage")
    monkeypatch.setattr(openmm_engine, "openmmml_version", lambda: "1.7")
    with pytest.raises(CapabilityError, match="does not convert it"):
        _bridge(model, precision="float64").validate(SimulationSpec())
    _bridge(model, precision="float32").validate(SimulationSpec())  # matches: accepted
    monkeypatch.setattr(openmm_engine, "openmmml_version", lambda: "1.8")
    _bridge(model, precision="float64").validate(SimulationSpec())  # 1.8 converts


def test_interaction_energy_is_refused_for_a_plain_mace_model(tmp_path):
    plain = _fake_model(tmp_path / "plain.model", "mace.modules.models MACE", "torch DoubleStorage")
    with pytest.raises(CapabilityError, match="interaction energy"):
        _bridge(plain).validate(SimulationSpec(energy_convention="interaction"))
    _bridge(plain).validate(SimulationSpec(energy_convention="total"))


def test_capabilities_match_what_the_route_executes(tmp_path):
    model = _fake_model(tmp_path / "m.model", "torch DoubleStorage")
    cpu = _bridge(model).capabilities()
    assert not cpu.stress and not cpu.per_atom_energy  # nothing but E and F comes back
    assert cpu.periodic and not cpu.partial_periodic and cpu.fixed_atoms
    assert not cpu.gpu
    # The model runs on the torch device whatever the OpenMM platform is.
    assert _bridge(model, device="cuda").capabilities().gpu


def test_execution_plan_reports_everything_without_importing_a_backend(
    tmp_path, rattled_nio_structure, backend_import_guard
):
    from ase.constraints import FixAtoms

    model = _fake_model(
        tmp_path / "m.model", "mace.modules.models ScaleShiftMACE", "torch DoubleStorage"
    )
    atoms = rattled_nio_structure.copy()
    atoms.set_constraint(FixAtoms(indices=[0, 1]))
    bridge = _bridge(model, engine=openmm_engine_spec(platform="CPU", threads=2))
    plan = bridge.execution_plan(nvt(thermostat="nose-hoover", thermostat_damping_fs=50.0), atoms)
    backend_import_guard()

    assert set(plan) >= PLAN_KEYS
    assert plan["potential_kind"] == "mace" and plan["engine"] == "openmm"
    assert plan["bridge"] == "MaceOpenMMBridge" and plan["implementation"] == "openmm-ml"
    assert plan["model_checkpoint"]["sha256"] and plan["model_checkpoint"]["stored_dtype"] == "float64"
    assert plan["exported_model"] is None and plan["lammps"] is None
    assert plan["elements"] == ["Ni", "O"] and plan["energy_convention"] == "total"
    assert plan["units"]["native"] == "openmm" and plan["units"]["pressure_unit"] is None
    assert GUARANTEE_KEYS <= set(plan["device"]) and GUARANTEE_KEYS <= set(plan["precision"])
    assert plan["precision"]["requested"] == "float64"
    assert plan["precision"]["platform_precision"]["passed_as"] is None
    dynamics = plan["dynamics"]
    assert DYNAMICS_KEYS <= set(dynamics)
    assert dynamics["integrator"] == "NoseHooverIntegrator"
    assert dynamics["removeCMMotion"] is False  # FixAtoms present
    assert dynamics["temperature_ndof"] == 3 * (len(atoms) - 2)
    assert dynamics["trajectory"]["expected_frames"] == 5
    assert plan["openmm"]["properties"] == {"Threads": "2"}
    kwargs = plan["openmm"]["create_system_kwargs"]
    assert kwargs["returnEnergyType"] == "energy" and kwargs["precision"] == "double"
    assert kwargs["device"] == "cpu" and kwargs["removeCMMotion"] is False
    assert plan["openmm"]["zero_mass_atoms"] == [0, 1]
    assert plan["openmm"]["box"]["lattice_change"]
    assert set(plan["availability"]) >= {"available", "missing"}
    assert plan["unmet_capabilities"] == []
    observed = plan["model_hashes"]["observed"]
    assert list(observed.values()) == [plan["model_checkpoint"]["sha256"]]
    assert plan["model_hashes"]["declared"] == {} and plan["model_hashes"]["match"] is None


def test_execution_plan_lists_refusals_instead_of_raising(tmp_path, nio_structure):
    model = _fake_model(tmp_path / "m.model", "torch DoubleStorage")
    slab = nio_structure.copy()
    slab.pbc = (True, True, False)
    plan = _bridge(model).execution_plan(nvt(thermostat="csvr"), slab)
    unmet = " | ".join(plan["unmet_capabilities"])
    assert "csvr" in unmet
    assert "partially periodic" in unmet or "partial_periodic" in unmet
    npt = SimulationSpec(
        task="md", ensemble="npt", temperature_K=300, pressure_bar=1.0, timestep_fs=1.0, steps=5
    )
    plan = _bridge(model).execution_plan(npt)
    assert any("stress" in item for item in plan["unmet_capabilities"])
    assert plan["dynamics"]["integrator"] is None


def test_execution_plan_checks_declared_hashes(tmp_path):
    model = _fake_model(tmp_path / "m.model", "torch DoubleStorage")
    wrong = _bridge(model, model_sha256="0" * 64).execution_plan(SimulationSpec())
    assert wrong["model_hashes"]["match"] is False
    right = _bridge(model).execution_plan(SimulationSpec())
    digest = right["model_checkpoint"]["sha256"]
    assert _bridge(model, model_sha256=digest).execution_plan(SimulationSpec())["model_hashes"][
        "match"
    ] is True


# ---------------------------------------------------------------------------
# Real OpenMM + OpenMM-ML + MACE (tiny random ScaleShiftMACE model)
# ---------------------------------------------------------------------------

E0 = {"O": -4.95, "Ni": -5.78}

#: Measured on this route (openmm 8.5.2, openmmml 1.7, mace 0.3.16, torch
#: 2.12.1 CPU, float64): |dE| <= 1.8e-14 eV and max|dF| <= 8.5e-16 eV/A against
#: MACE's ASE calculator. The limits leave ~4 orders of magnitude of room for
#: platform/thread-count differences while still catching any unit, rotation,
#: convention or E0 error (the smallest of which, the -3.3e-7 relative bias of
#: an exact eV<->kJ/mol constant, is ~3e-5 eV here).
ENERGY_TOLERANCE_EV = 1e-9
FORCE_TOLERANCE_EV_PER_A = 1e-10


def _make_tiny_model(dtype_name: str, path: Path) -> Path:
    import torch
    from e3nn import o3
    from mace import modules, tools

    dtype = getattr(torch, dtype_name)
    previous = torch.get_default_dtype()
    torch.set_default_dtype(dtype)
    try:
        torch.manual_seed(0)
        table = tools.AtomicNumberTable([8, 28])
        model = modules.ScaleShiftMACE(
            r_max=4.0,
            num_bessel=8,
            num_polynomial_cutoff=6,
            max_ell=2,
            interaction_cls=modules.interaction_classes["RealAgnosticResidualInteractionBlock"],
            interaction_cls_first=modules.interaction_classes[
                "RealAgnosticResidualInteractionBlock"
            ],
            num_interactions=2,
            num_elements=2,
            hidden_irreps=o3.Irreps("8x0e + 8x1o"),
            MLP_irreps=o3.Irreps("8x0e"),
            gate=torch.nn.functional.silu,
            atomic_energies=np.array([E0["O"], E0["Ni"]]),
            avg_num_neighbors=8.0,
            atomic_numbers=table.zs,
            correlation=2,
            radial_type="bessel",
            atomic_inter_scale=1.0,
            atomic_inter_shift=0.0,
        )
        torch.save(model, path)
    finally:
        torch.set_default_dtype(previous)
    return path


@pytest.fixture(scope="session")
def tiny_models(require_mace, require_openmm, tmp_path_factory):
    os.environ.setdefault("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", "1")
    root = tmp_path_factory.mktemp("tiny_mace")
    return {
        "float64": _make_tiny_model("float64", root / "tiny64.model"),
        "float32": _make_tiny_model("float32", root / "tiny32.model"),
    }


def _structures():
    from ase.build import bulk

    ortho = bulk("NiO", "rocksalt", a=4.17, cubic=True).repeat((2, 1, 1))
    ortho.rattle(0.05, seed=7)
    fcc = bulk("NiO", "rocksalt", a=4.17).repeat((2, 2, 2))
    fcc.rattle(0.05, seed=20250917)
    cluster = bulk("NiO", "rocksalt", a=4.17, cubic=True)
    cluster.rattle(0.05, seed=3)
    cluster.pbc = False
    return {"ortho": ortho, "fcc-primitive": fcc, "cluster": cluster}


def _mace_ase(model: Path, atoms, dtype: str = "float64"):
    from mace.calculators import MACECalculator

    reference = atoms.copy()
    reference.calc = MACECalculator(model_paths=[str(model)], device="cpu", default_dtype=dtype)
    return reference.get_potential_energy(), reference.get_forces()


def _read_dcd(path: Path, n_atoms: int):
    data = Path(path).read_bytes()
    box = struct.unpack("<i", data[48:52])[0]
    offset, frames = 276, []
    while offset < len(data):
        if box:
            offset += 56
        coordinates = []
        for _ in range(3):
            offset += 4
            coordinates.append(np.frombuffer(data, dtype="<f4", count=n_atoms, offset=offset))
            offset += 4 * n_atoms + 4
        frames.append(np.stack(coordinates, axis=1))
    return frames


@pytest.mark.openmm
@pytest.mark.mace
@pytest.mark.parametrize("platform", ["Reference", "CPU"])
@pytest.mark.parametrize("structure", ["ortho", "fcc-primitive", "cluster"])
@pytest.mark.parametrize("convention", ["total", "interaction"])
def test_singlepoint_matches_mace_ase(tiny_models, platform, structure, convention):
    atoms = _structures()[structure]
    model = tiny_models["float64"]
    energy, forces = _mace_ase(model, atoms)
    if convention == "interaction":
        energy -= sum(E0[s] for s in atoms.get_chemical_symbols())
    bridge = _bridge(model, engine=openmm_engine_spec(platform=platform))
    result = bridge.singlepoint(atoms, SimulationSpec(energy_convention=convention))
    assert result.energy_convention == convention
    assert abs(result.energy_eV - energy) <= ENERGY_TOLERANCE_EV
    assert np.abs(np.asarray(result.forces_eV_per_A) - forces).max() <= FORCE_TOLERANCE_EV_PER_A
    extras = result.extras
    assert extras["platform"]["name"] == platform
    assert extras["returnEnergyType"] == openmm_engine.ENERGY_TYPE[convention]
    assert extras["model_evaluation_dtype"] == "float64" and extras["model_precision"]["guaranteed"]
    if structure == "fcc-primitive":
        assert not np.allclose(extras["box"]["rotation"], np.eye(3))  # forces were rotated back
    if structure == "cluster":
        assert extras["box"] is None
        assert not extras["system"]["readback"]["uses_periodic_boundary_conditions"]


@pytest.mark.openmm
@pytest.mark.mace
def test_openmm_itself_rejects_precision_on_cpu(require_openmm):
    """Why the builder refuses it: OpenMM raises instead of ignoring the key."""
    import openmm

    system = openmm.System()
    system.addParticle(1.0)
    with pytest.raises(openmm.OpenMMException, match="Illegal property name"):
        openmm.Context(
            system,
            openmm.VerletIntegrator(0.001),
            openmm.Platform.getPlatformByName("CPU"),
            {"Precision": "double"},
        )


@pytest.mark.openmm
@pytest.mark.mace
def test_threads_reach_the_cpu_platform(tiny_models):
    atoms = _structures()["ortho"]
    bridge = _bridge(tiny_models["float64"], engine=openmm_engine_spec(platform="CPU", threads=2))
    platform = bridge.singlepoint(atoms, SimulationSpec()).extras["platform"]
    assert platform["properties"]["Threads"] == "2"  # read back from the Context
    assert platform["requested_properties"] == {"Threads": "2"}
    assert platform["requested_properties_verified"]


@pytest.mark.openmm
@pytest.mark.mace
def test_an_unavailable_platform_is_refused_with_the_available_ones(tiny_models):
    available = openmm_engine.available_platforms()
    missing = [name for name in ("CUDA", "HIP", "OpenCL") if name not in available]
    if not missing:
        pytest.skip("every OpenMM platform is installed here")
    bridge = _bridge(tiny_models["float64"], engine=openmm_engine_spec(platform=missing[0]))
    with pytest.raises(ConfigError, match="available platforms: " + ", ".join(available)):
        bridge.singlepoint(_structures()["ortho"], SimulationSpec())


@pytest.mark.openmm
@pytest.mark.mace
def test_a_precision_the_model_is_not_stored_in_is_refused_on_openmmml_1_7(tiny_models):
    version = openmm_engine.openmmml_version()
    if openmm_engine.version_tuple(version) >= openmm_engine.OPENMMML_CONVERTS_MODEL_DTYPE:
        pytest.skip(f"openmmml {version} converts the model to the requested dtype")
    atoms = _structures()["ortho"]
    bridge = _bridge(tiny_models["float32"], precision="float64")
    with pytest.raises(CapabilityError, match="does not convert it"):
        bridge.validate(SimulationSpec(), atoms)
    plan = bridge.execution_plan(SimulationSpec(), atoms)
    assert not plan["precision"]["guaranteed"] and plan["precision"]["effective"] is None
    assert any("does not convert it" in item for item in plan["unmet_capabilities"])


@pytest.mark.openmm
@pytest.mark.mace
def test_the_upstream_dtype_behaviour_the_refusal_relies_on(tiny_models):
    """openmmml 1.7 casts the inputs only: a float32 model asked for 'double' fails."""
    version = openmm_engine.openmmml_version()
    if openmm_engine.version_tuple(version) >= openmm_engine.OPENMMML_CONVERTS_MODEL_DTYPE:
        pytest.skip(f"openmmml {version} converts the model to the requested dtype")
    import openmm
    from openmmml import MLPotential

    atoms = _structures()["ortho"]
    box = openmm_engine.reduced_box(atoms.get_cell())
    topology = openmm_engine.build_topology(atoms, box)
    system = MLPotential("mace", modelPath=str(tiny_models["float32"])).createSystem(
        topology, precision="double", device="cpu", returnEnergyType="energy"
    )
    context = openmm.Context(
        system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference")
    )
    context.setPositions(openmm_engine.positions_nm(atoms, box))
    with pytest.raises(openmm.OpenMMException):
        context.getState(getEnergy=True)


@pytest.mark.openmm
@pytest.mark.mace
def test_float32_evaluation_matches_mace_ase_in_float32(tiny_models):
    atoms = _structures()["fcc-primitive"]
    energy, forces = _mace_ase(tiny_models["float32"], atoms, dtype="float32")
    bridge = _bridge(tiny_models["float32"], precision="float32")
    result = bridge.singlepoint(atoms, SimulationSpec())
    assert result.extras["model_evaluation_dtype"] == "float32"
    # float32 arithmetic on ~83 eV: a few ulps of the total energy.
    assert abs(result.energy_eV - energy) <= 1e-4
    assert np.abs(np.asarray(result.forces_eV_per_A) - forces).max() <= 1e-5


@pytest.mark.openmm
@pytest.mark.mace
def test_fixed_atoms_stay_fixed_in_langevin_dynamics(tiny_models, tmp_path):
    from ase.constraints import FixAtoms

    atoms = _structures()["ortho"]
    atoms.set_constraint(FixAtoms(indices=[0, 1, 2, 3]))
    bridge = _bridge(tiny_models["float64"], engine=openmm_engine_spec(platform="Reference"))
    trajectory = bridge.run_md(atoms, nvt(), workdir=tmp_path)
    frames = _read_dcd(trajectory.trajectory_path, len(atoms))
    assert len(frames) == trajectory.frames_written == 5
    for frame in frames[1:]:
        assert np.array_equal(frame[:4], frames[0][:4])  # massless: never moved
    assert max(np.abs(frame[4:] - frames[0][4:]).max() for frame in frames) > 1e-4
    final = np.asarray(trajectory.final_positions_angstrom)
    assert np.array_equal(final[:4], atoms.get_positions()[:4])
    assert trajectory.temperature_ndof == 3 * (len(atoms) - 4)
    resolved = trajectory.integrator_resolved
    assert resolved["name"] == "LangevinMiddleIntegrator" and resolved["removeCMMotion"] is False
    assert resolved["system"]["readback"]["zero_mass_particles"] == [0, 1, 2, 3]
    assert trajectory.constraints["fixed_atoms"] == [0, 1, 2, 3]


@pytest.mark.openmm
@pytest.mark.mace
@pytest.mark.parametrize("steps, interval, expected", [(20, 5, 5), (7, 5, 2), (6, 3, 3), (3, 10, 1)])
def test_frame_count_is_exact(tiny_models, tmp_path, steps, interval, expected):
    atoms = _structures()["ortho"]
    bridge = _bridge(tiny_models["float64"], engine=openmm_engine_spec(platform="Reference"))
    trajectory = bridge.run_md(
        atoms, nvt(steps=steps, trajectory_interval=interval), workdir=tmp_path
    )
    assert trajectory.steps_completed == steps
    assert trajectory.frames_written == expected
    assert openmm_engine.dcd_frame_count(trajectory.trajectory_path, len(atoms)) == expected
    assert len(_read_dcd(trajectory.trajectory_path, len(atoms))) == expected


@pytest.mark.openmm
@pytest.mark.mace
def test_a_non_periodic_cluster_runs_and_writes_a_boxless_dcd(tiny_models, tmp_path):
    atoms = _structures()["cluster"]
    bridge = _bridge(tiny_models["float64"], engine=openmm_engine_spec(platform="Reference"))
    trajectory = bridge.run_md(atoms, nvt(steps=10), workdir=tmp_path)
    assert trajectory.frames_written == 3
    assert struct.unpack("<i", Path(trajectory.trajectory_path).read_bytes()[48:52])[0] == 0


@pytest.mark.openmm
@pytest.mark.mace
def test_nve_conserves_energy_and_nose_hoover_reports_its_conserved_quantity(
    tiny_models, tmp_path
):
    atoms = _structures()["fcc-primitive"]
    bridge = _bridge(tiny_models["float64"], engine=openmm_engine_spec(platform="Reference"))
    nve = bridge.run_md(
        atoms,
        SimulationSpec(
            task="md", ensemble="nve", temperature_K=300.0, timestep_fs=0.5, steps=40, seed=11,
            trajectory_interval=10, log_interval=2,
        ),
        workdir=tmp_path / "nve",
    )
    assert nve.diagnostics["conservation_test"]
    # Measured 7.6e-7 eV/atom at 0.5 fs (second order in dt); 1e-5 is a
    # loose physical bound that a wrong force/unit/rotation would exceed.
    assert nve.diagnostics["max_abs_energy_excursion_eV_per_atom"] < 1e-5
    assert nve.integrator_resolved["name"] == "VerletIntegrator"
    assert nve.temperature_ndof == 3 * len(atoms) - 3  # CMMotionRemover present

    excursions = []
    for dt in (0.5, 0.25):
        hoover = bridge.run_md(
            atoms,
            nvt(thermostat="nose-hoover", timestep_fs=dt, steps=int(20 / dt), log_interval=1),
            workdir=tmp_path / f"nh{dt}",
        )
        resolved = hoover.integrator_resolved
        assert resolved["name"] == "NoseHooverIntegrator"
        assert resolved["openmm_thermostat_ndof"] == hoover.temperature_ndof == 3 * len(atoms) - 3
        assert resolved["collision_frequency_per_ps"] == pytest.approx(10.0)
        assert hoover.diagnostics["conserved_quantity"] == openmm_engine.NOSE_HOOVER_CONSERVED_QUANTITY
        excursions.append(
            hoover.diagnostics["conserved_quantity_drift"]["max_abs_energy_excursion_eV_per_atom"]
        )
    # OpenMM's NH kinetic energy is not synchronised with the positions, so
    # the heat-bath-corrected energy converges only at first order in dt
    # (measured 5.2e-5 -> 2.6e-5 eV/atom); it must at least converge.
    assert excursions[0] < 1e-4
    assert excursions[1] <= 0.6 * excursions[0]


@pytest.mark.openmm
@pytest.mark.mace
def test_single_points_and_md_record_identical_execution_settings(tiny_models, tmp_path):
    """Through jobs.run_smoke_md, so the manifest itself is what is checked."""
    import json

    from ase.io import write

    from nio_md_prep.mlip import jobs
    from nio_md_prep.mlip.specs import JobSpec, StructureSpec

    atoms = _structures()["fcc-primitive"]
    structure = tmp_path / "nio.xyz"
    write(str(structure), atoms, format="extxyz")
    job = JobSpec(
        potential=MacePotentialSpec(
            model_path=tiny_models["float64"], declared_elements=("Ni", "O"), precision="float64"
        ),
        engine=openmm_engine_spec(platform="CPU", threads=2),
        simulation=nvt(steps=10),
        name="openmm-settings",
    )
    report = jobs.run_smoke_md(
        job, output_dir=tmp_path / "run", structure=StructureSpec(path=structure)
    )
    manifest = json.loads(Path(report["manifest_path"]).read_text(encoding="utf-8"))
    trajectory = manifest["results"]["trajectory"]
    md_settings = trajectory["integrator_resolved"]["execution_settings"]
    assert trajectory["initial"]["extras"]["execution_settings"] == md_settings
    assert trajectory["final"]["extras"]["execution_settings"] == md_settings
    assert md_settings["platform"] == "CPU"
    assert md_settings["platform_properties"]["Threads"] == "2"
    assert md_settings["createSystem"]["precision"] == "double"
    assert md_settings["createSystem"]["device"] == "cpu"
    assert md_settings["model_dtype"] == "float64" and md_settings["model_device"] == "cpu"
    assert manifest["engine_parameters"]["model_precision"]["guaranteed"]
    assert trajectory["integrator_resolved"]["integrator_seed"] >= 1
