"""The shared execution-plan contract and the capability audit, route by route.

Every bridge's ``execution_plan`` carries the keys in
:data:`~nio_md_prep.mlip.bridges.base.EXECUTION_PLAN_KEYS`, computed without
running an engine or importing a backend: LAMMPS plans render their pair
commands and launch line with no LAMMPS installed, MACE plans report a missing
torch as unavailability rather than failing. The capability tests pin claims
to what each route actually executes (energy convention, GPU, precision,
per-atom energies).
"""
from __future__ import annotations

import json
import textwrap

import pytest

from nio_md_prep.mlip import environment, jobs
from nio_md_prep.mlip.bridges.base import (
    EXECUTION_PLAN_KEYS,
    PLAN_DYNAMICS_KEYS,
    PLAN_GUARANTEE_KEYS,
    Bridge,
    complete_plan,
    plan_problems,
)
from nio_md_prep.mlip.capabilities import CapabilitySet, RequirementSet
from nio_md_prep.mlip.config import parse_job
from nio_md_prep.mlip.engines import lammps_engine
from nio_md_prep.mlip.environment import Availability
from nio_md_prep.mlip.errors import CapabilityError
from nio_md_prep.mlip.registry import build_bridge
from nio_md_prep.mlip.results import PotentialResult
from nio_md_prep.mlip.specs import (
    EngineSpec,
    LammpsMlipPotentialSpec,
    MacePotentialSpec,
    MockPotentialSpec,
    SimulationSpec,
)

pytest.importorskip("ase")

SINGLEPOINT = SimulationSpec(task="singlepoint")
NVT = SimulationSpec(
    task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=1.0, steps=10, seed=3,
    thermostat_damping_fs=50.0,
)
LJ = LammpsMlipPotentialSpec(
    pair_style="lj/cut 6.0", pair_coeff=("* * 0.01 2.5",), type_map={1: "Ni", 2: "O"}
)


def fake_model(tmp_path, name="nio.model"):
    """A model file that exists but that nothing can load: plans must not try."""
    path = tmp_path / name
    path.write_bytes(b"not a real MACE checkpoint")
    return path


def mace(tmp_path, **kwargs) -> MacePotentialSpec:
    return MacePotentialSpec(
        label="nio", model_path=fake_model(tmp_path), declared_elements=("Ni", "O"), **kwargs
    )


def route_bridges(tmp_path):
    """One bridge per registered route (and the GPU variants), no backend needed."""
    return {
        "mock-ase": build_bridge(MockPotentialSpec(), EngineSpec(kind="ase")),
        "lammps-native": build_bridge(LJ, EngineSpec(kind="lammps")),
        "lammps-native-executable": build_bridge(LJ, EngineSpec(kind="lammps", runtime="executable")),
        "lammps-ase": build_bridge(LJ, EngineSpec(kind="ase")),
        "mace-ase": build_bridge(mace(tmp_path), EngineSpec(kind="ase")),
        "mace-openmm": build_bridge(mace(tmp_path), EngineSpec(kind="openmm")),
        "mace-mliap": build_bridge(mace(tmp_path), EngineSpec(kind="lammps")),
        "mace-mliap-cuda": build_bridge(
            mace(tmp_path, device="cuda"), EngineSpec(kind="lammps", runtime="executable")
        ),
        "mace-pair": build_bridge(
            mace(tmp_path), EngineSpec(kind="lammps", runtime="executable"),
            implementation="pair-mace",
        ),
    }


ROUTES = (
    "mock-ase", "lammps-native", "lammps-native-executable", "lammps-ase", "mace-ase",
    "mace-openmm", "mace-mliap", "mace-mliap-cuda", "mace-pair",
)


# --- the contract --------------------------------------------------------------


@pytest.mark.parametrize("route", ROUTES)
@pytest.mark.parametrize("simulation", [SINGLEPOINT, NVT], ids=["singlepoint", "nvt"])
def test_every_route_returns_the_full_plan_contract(tmp_path, nio_structure, route, simulation):
    bridge = route_bridges(tmp_path)[route]
    plan = bridge.execution_plan(simulation, nio_structure)
    assert plan_problems(plan) == []
    assert set(EXECUTION_PLAN_KEYS) <= set(plan)
    assert plan["bridge"] == type(bridge).__name__
    assert plan["implementation"] == bridge.implementation
    assert plan["label"] == bridge.label
    for key in ("device", "precision"):
        assert set(PLAN_GUARANTEE_KEYS) <= set(plan[key])
        assert isinstance(plan[key]["guaranteed"], bool)
    if simulation.task == "md":
        assert set(PLAN_DYNAMICS_KEYS) <= set(plan["dynamics"])
        assert plan["dynamics"]["timestep_fs"] == 1.0
        assert plan["dynamics"]["timestep_native"] is not None
    else:
        assert plan["dynamics"] is None
    assert (plan["lammps"] is None) is (bridge.engine_kind != "lammps" and route != "lammps-ase")
    assert (plan["openmm"] is None) is (bridge.engine_kind != "openmm")
    json.dumps(plan, default=str)  # the --json view must serialise


@pytest.mark.parametrize("route", ROUTES)
def test_plan_refusals_are_listed_not_raised(tmp_path, nio_structure, route):
    """A structure with an uncovered element: listed as unmet, and validate still refuses."""
    from ase import Atom

    contaminated = nio_structure.copy()
    contaminated.append(Atom("F", (0.3, 0.3, 0.3)))
    bridge = route_bridges(tmp_path)[route]
    plan = bridge.execution_plan(SINGLEPOINT, contaminated)
    if bridge.capabilities().elements is None:
        pytest.skip("route does not restrict elements")
    assert any("F" in item for item in plan["unmet_capabilities"])
    with pytest.raises(CapabilityError):
        bridge.validate(SINGLEPOINT, contaminated)


def test_lammps_plans_render_pair_commands_and_launch_without_lammps(tmp_path, nio_structure):
    plans = {
        name: route_bridges(tmp_path)[name].execution_plan(NVT, nio_structure)
        for name in ("lammps-native-executable", "mace-mliap-cuda", "mace-pair")
    }
    native = plans["lammps-native-executable"]["lammps"]
    assert native["pair_style"] == "pair_style lj/cut 6.0"
    assert native["pair_coeff"] == ["pair_coeff * * 0.01 2.5"]
    assert native["launch_command"].startswith("lmp")
    assert native["kokkos_args"] == [] and native["env"] == {}

    mliap = plans["mace-mliap-cuda"]["lammps"]
    assert mliap["pair_style"].startswith("pair_style mliap unified nio.model-mliap_lammps.pt")
    assert mliap["pair_coeff"] == ["pair_coeff * * Ni O"]
    assert mliap["kokkos_args"] == list(lammps_engine.KOKKOS_GPU_ARGS)
    assert " ".join(lammps_engine.KOKKOS_GPU_ARGS) in mliap["launch_command"]

    pair = plans["mace-pair"]["lammps"]
    assert pair["pair_style"] == "pair_style mace no_domain_decomposition"
    assert pair["env"] == {"CUDA_VISIBLE_DEVICES": ""}
    assert pair["launch_command"].startswith("CUDA_VISIBLE_DEVICES")
    # The MD section is in LAMMPS metal units: 1 fs = 0.001 ps.
    assert plans["mace-pair"]["dynamics"]["timestep_native"] == 0.001
    assert plans["mace-pair"]["dynamics"]["thermostat_damping_native"] == pytest.approx(0.05)


def hide_modules(monkeypatch, *names):
    """Make ``find_spec`` report ``names`` as absent, as on a machine without them."""
    real = environment.find_spec

    def find_spec(name, *args, **kwargs):
        if name.split(".", 1)[0] in names:
            return None
        return real(name, *args, **kwargs)

    monkeypatch.setattr(environment, "find_spec", find_spec)


@pytest.mark.parametrize("engine", ["ase", "openmm"])
def test_mace_plans_without_torch_report_unavailability_honestly(
    tmp_path, nio_structure, monkeypatch, engine
):
    hide_modules(monkeypatch, "torch", "mace", "openmm", "openmmml")
    bridge = build_bridge(mace(tmp_path), EngineSpec(kind=engine))
    plan = bridge.execution_plan(NVT, nio_structure)
    assert plan_problems(plan) == []
    assert plan["availability"]["available"] is False
    assert "torch" in plan["availability"]["missing"]
    # Nothing is promised about a runtime that is not there.
    assert plan["device"]["guaranteed"] is False
    assert plan["precision"]["guaranteed"] is False
    assert plan["model_checkpoint"]["sha256"] is not None


def test_mace_lammps_plan_without_torch_does_not_claim_the_export_dtype(
    tmp_path, nio_structure, monkeypatch
):
    hide_modules(monkeypatch, "torch", "mace")
    spec = mace(tmp_path)
    (tmp_path / "nio.model-mliap_lammps.pt").write_bytes(b"export")
    plan = build_bridge(spec, EngineSpec(kind="lammps")).execution_plan(SINGLEPOINT, nio_structure)
    assert plan_problems(plan) == []
    assert plan["exported_model"]["sha256"] is not None
    assert plan["precision"]["effective"] is None
    assert plan["precision"]["guaranteed"] is False


class _MinimalBridge(Bridge):
    potential_kind = "mock"
    engine_kind = "ase"
    implementation = "minimal"

    def capabilities(self) -> CapabilitySet:
        return CapabilitySet(energy=True, forces=True, elements=frozenset({"Ni", "O"}))

    def availability(self) -> Availability:
        return Availability(False, ("nothing",), "a test route")

    def singlepoint(self, atoms, simulation):  # pragma: no cover - never run
        raise AssertionError("a plan must not execute")


def test_the_default_plan_satisfies_the_contract(nio_structure):
    bridge = _MinimalBridge(MockPotentialSpec(), EngineSpec(kind="ase"))
    plan = bridge.execution_plan(NVT, nio_structure)
    assert plan_problems(plan) == []
    assert plan["bridge"] == "_MinimalBridge" and plan["implementation"] == "minimal"
    assert plan["device"]["guaranteed"] is False and plan["precision"]["guaranteed"] is False
    assert plan["dynamics"]["ensemble"] == "nvt" and plan["dynamics"]["integrator"] is None
    assert plan["availability"] == {"available": False, "missing": ["nothing"], "detail": "a test route"}
    assert plan["lammps"] is None and plan["openmm"] is None


def test_complete_plan_fills_gaps_and_never_invents_a_guarantee():
    plan = complete_plan({"device": {"requested": "cuda"}, "extra": 1})
    assert plan_problems(plan) == []
    assert plan["device"] == {"requested": "cuda", "effective": None, "guaranteed": False, "note": None}
    assert plan["precision"]["guaranteed"] is False
    assert plan["availability"]["available"] is False
    assert plan["unmet_capabilities"] == [] and plan["extra"] == 1
    assert plan["exported_model"] is None and plan["dynamics"] is None
    assert plan_problems({"device": None}) != []


# --- validate_job and the CLI --------------------------------------------------

MOCK_MD = """
[potential]
kind = "mock"
elements = ["Ni", "O"]

[engine]
kind = "ase"

[simulation]
task = "md"
ensemble = "{ensemble}"
temperature_K = 300
pressure_bar = 1.0
timestep_fs = 1.0
steps = 5

[structure]
path = "{structure}"
"""


def write_job(tmp_path, atoms, ensemble="nvt"):
    from ase.io import write

    structure = tmp_path / "structure.xyz"
    write(str(structure), atoms, format="extxyz")
    config = tmp_path / "job.toml"
    config.write_text(
        textwrap.dedent(MOCK_MD).format(ensemble=ensemble, structure=structure.as_posix()),
        encoding="utf-8",
    )
    return config


def test_validate_job_returns_the_execution_plan(tmp_path, rattled_nio_structure):
    report = jobs.validate_job(parse_job(write_job(tmp_path, rattled_nio_structure)))
    assert report["ok"]
    plan = report["execution_plan"]
    assert plan_problems(plan) == []
    assert plan["dynamics"]["integrator"] and plan["unmet_capabilities"] == []


def test_non_strict_validation_shows_the_plan_of_a_refused_job(tmp_path, rattled_nio_structure):
    slab = rattled_nio_structure.copy()
    slab.pbc = (True, True, False)
    job = parse_job(write_job(tmp_path, slab, ensemble="npt"))
    with pytest.raises(Exception):
        jobs.validate_job(job)
    report = jobs.validate_job(job, strict=False)
    assert report["ok"] is False
    assert report["execution_plan"]["unmet_capabilities"]
    assert report["unmet"] == report["execution_plan"]["unmet_capabilities"]


def test_the_cli_prints_a_readable_plan_and_json(tmp_path, rattled_nio_structure, capsys):
    from nio_md_prep.cli import main

    config = write_job(tmp_path, rattled_nio_structure)
    assert main(["mlip", "validate", str(config)]) == 0
    text = capsys.readouterr().out
    for label in ("Execution plan", "route:", "checkpoint:", "exported:", "elements:", "energy:",
                  "units:", "device:", "precision:", "dynamics:", "timestep", "available:",
                  "unmet:", "model hashes:"):
        assert label in text, label
    assert main(["mlip", "validate", str(config), "--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert plan_problems(report["execution_plan"]) == []


def test_the_cli_shows_the_plan_then_refuses(tmp_path, rattled_nio_structure, capsys):
    from nio_md_prep.cli import main

    slab = rattled_nio_structure.copy()
    slab.pbc = (True, True, False)
    config = write_job(tmp_path, slab, ensemble="npt")
    with pytest.raises(SystemExit) as excinfo:
        main(["mlip", "validate", str(config)])
    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "Execution plan" in captured.out and "INVALID" in captured.out
    assert "error:" in captured.err and "Traceback" not in captured.err


# --- capability audit: claims equal execution ----------------------------------


@pytest.mark.parametrize("implementation", ["mliap", "pair-mace"])
def test_mace_on_lammps_reports_the_total_energy_and_refuses_interaction(tmp_path, implementation):
    """Both exports sum MACE node energies, which include the E0s: total, never relabelled."""
    spec = mace(tmp_path, energy_convention="interaction", atomic_reference_energies={"Ni": -1.0, "O": -2.0})
    bridge = build_bridge(spec, EngineSpec(kind="lammps"), implementation=implementation)
    capabilities = bridge.capabilities()
    assert capabilities.native_energy_convention == "total"
    assert capabilities.energy_conventions == frozenset({"total"})
    assert bridge.as_lammps_spec().energy_convention == "total"
    assert any("energy convention 'interaction'" in p for p in bridge.plan_refusals(SINGLEPOINT))
    asked = SimulationSpec(task="singlepoint", energy_convention="interaction")
    plain = build_bridge(mace(tmp_path), EngineSpec(kind="lammps"), implementation=implementation)
    with pytest.raises(CapabilityError, match="interaction"):
        plain.validate(asked)
    assert plain.execution_plan(SINGLEPOINT, None)["energy_convention"] == "total"


def test_pair_mace_claims_a_gpu_only_where_the_device_is_checked(tmp_path):
    executable = build_bridge(
        mace(tmp_path, device="cuda"), EngineSpec(kind="lammps", runtime="executable"),
        implementation="pair-mace",
    )
    python_route = build_bridge(
        mace(tmp_path, device="cuda"), EngineSpec(kind="lammps", runtime="python"),
        implementation="pair-mace",
    )
    assert executable.capabilities().gpu
    assert not python_route.capabilities().gpu
    assert python_route.execution_plan(SINGLEPOINT, None)["device"]["guaranteed"] is False


def test_mliap_claims_a_gpu_only_with_a_kokkos_gpu_launch(tmp_path):
    cpu = build_bridge(mace(tmp_path), EngineSpec(kind="lammps"))
    cuda = build_bridge(mace(tmp_path, device="cuda"), EngineSpec(kind="lammps"))
    assert not cpu.capabilities().gpu and cpu.launch().accelerator["kokkos_gpus"] == 0
    assert cuda.capabilities().gpu and cuda.launch().accelerator["kokkos_gpus"] >= 1
    for bridge in (cpu, cuda):
        assert not bridge.capabilities().per_atom_energy


@pytest.mark.parametrize("engine", ["lammps", "ase"])
def test_lammps_pair_styles_guarantee_no_precision(engine):
    bridge = build_bridge(LJ, EngineSpec(kind=engine))
    assert bridge.capabilities().precisions == frozenset()
    problems = RequirementSet(precision="float64").unmet(bridge.capabilities())
    assert problems and "guarantees no precision" in problems[0]
    assert RequirementSet().unmet(bridge.capabilities()) == []


def test_a_lammps_native_route_delivers_the_convention_it_claims():
    """Declared E0s make 'interaction' reachable; the result must then be in it."""
    spec = LammpsMlipPotentialSpec(
        pair_style="lj/cut 6.0", pair_coeff=("* * 0.01 2.5",), type_map={1: "Ni", 2: "O"},
        atomic_reference_energies={"Ni": -1.5, "O": -0.5},
    )
    bridge = build_bridge(spec, EngineSpec(kind="lammps"))
    assert "interaction" in bridge.capabilities().energy_conventions
    raw = PotentialResult(
        energy_eV=-10.0, forces_eV_per_A=[[0.0, 0.0, 0.0]] * 2, symbols=["Ni", "O"],
        energy_convention="total", engine="lammps", potential="lammps",
        implementation="native", native_units="lammps_metal",
    )
    asked = SimulationSpec(task="singlepoint", energy_convention="interaction")
    reported = bridge.in_requested_convention(raw, asked)
    assert reported.energy_convention == "interaction"
    assert reported.energy_eV == pytest.approx(-10.0 - (-1.5 - 0.5))
    assert reported.extras["native_energy_eV"] == -10.0
    assert bridge.in_requested_convention(raw, SINGLEPOINT) is raw
    assert bridge.execution_plan(asked, None)["energy_convention"] == "interaction"


def test_gpu_claims_match_the_plans_device(tmp_path, nio_structure):
    """A route that claims a GPU plans a GPU; one that plans the CPU claims none."""
    for name, bridge in route_bridges(tmp_path).items():
        plan = bridge.execution_plan(SINGLEPOINT, nio_structure)
        effective = str(plan["device"]["effective"] or "")
        if bridge.capabilities().gpu:
            assert effective.startswith(("cuda", "gpu")), name
        if effective == "cpu" or effective.startswith("cpu "):
            assert not bridge.capabilities().gpu, name
