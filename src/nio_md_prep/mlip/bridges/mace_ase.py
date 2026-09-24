"""MACE on ASE: the direct route, and the reference for cross-engine checks.

MACE ships its own ASE calculator, so this bridge is thin -- which is exactly
why it is the reference. Its units are already canonical and it introduces no
export step that could change a number. The other two MACE routes are
compared against this one.

**Energy convention.** ``MACECalculator``'s ``energy`` is always the *total*
energy: ``ScaleShiftMACE`` returns ``sum(E0) + interaction``
(``mace/modules/models.py``), and the per-atom ``energies`` include E0 too.
So the native convention of this route is ``total``, whatever the
configuration says. ``potential.energy_convention`` (or, taking precedence,
``simulation.energy_convention``) is the convention to *report*: an
``interaction`` request is answered by subtracting the model's own E0s for
the evaluated head -- read from the model file, never the declared ones,
which are only cross-checked (1e-8 eV) -- and refused when the model's E0s
cannot be read.

**What is passed to MACE, and read back.** The calculator is built with
``model_paths``, ``device``, ``default_dtype`` (``potential.precision``),
``head`` (resolved and validated against the model's heads: MACE itself
falls back to the *last* head for an unknown one, with only a warning) and
``compile_mode`` (validated against torch.compile's modes). The device is
pre-checked (CUDA/MPS present, ordinal in range) before torch deserialises
anything onto it, as part of this route's availability.
``engine.threads`` sets torch's intra-op thread count. After the build, the
evaluated head, the parameters' device and dtype, the calculator's
``default_dtype`` and ``use_compile`` and torch's thread count are read back
from the live calculator and recorded with every result, next to the model
file's native dtype -- a float32-trained model evaluated in float64 is
up-cast, and says so.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..errors import CapabilityError, ConfigError, MlipError
from ..potentials.mace import MaceAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import ASE, INTERACTION, TOTAL
from ..engines import ase_engine
from ..engines.ase_engine import AseEngine
from .base import Bridge

CALCULATOR = "mace.calculators.MACECalculator"


class MaceAseBridge(Bridge):
    potential_kind = "mace"
    engine_kind = "ase"
    implementation = "mace-ase-calculator"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = MaceAdapter(potential)
        self.runtime = AseEngine(engine)
        self._calculator = None
        self._runtime_report: dict | None = None

    # -- description ------------------------------------------------------

    def capabilities(self) -> CapabilitySet:
        return self.adapter.capabilities().intersect(
            self._route_capabilities(
                self.runtime.capabilities(),
                energy_convention=TOTAL,
                native_units=ASE.name,
                notes=(
                    "MACECalculator reports the total energy (the model's E0s included); "
                    "the interaction energy is total minus the model's own E0s",
                ),
            )
        )

    def availability(self) -> Availability:
        model = self.adapter.availability()
        if not model:
            return model
        engine = self.runtime.availability()
        if not engine:
            return engine
        problem = self._device_problem()
        if problem:
            return Availability(available=False, missing=(), detail=problem)
        return engine

    def atomic_reference_energies(self):
        """The model's own E0s for the evaluated head: what MACE adds to its total."""
        return self.adapter.model_atomic_reference_energies()

    def reported_convention(self, simulation: SimulationSpec | None = None) -> str:
        if simulation is not None and simulation.energy_convention:
            return simulation.energy_convention
        return self.potential.energy_convention

    def _e0_unavailable_reason(self) -> str:
        availability = self.adapter.availability()
        if not availability:
            return availability.detail or f"missing {', '.join(availability.missing)}"
        if self.adapter.discovery_error:
            return f"the model file could not be read ({self.adapter.discovery_error})"
        discovered = self.adapter._discovered_if_available() or {}
        if discovered.get("head_error"):
            return discovered["head_error"]
        return "the model's atomic-energies layout is not recognised"

    # -- gating -----------------------------------------------------------

    def requirements(self, simulation: SimulationSpec, atoms=None):
        """The shared requirements, plus the convention this route must report.

        ``potential.energy_convention`` is a reporting request on this route
        (the calculator's own convention is always ``total``), so it becomes
        the requirement when the simulation does not name one. An
        ``interaction`` request with no readable model E0s is refused here,
        with the reason, before negotiation reduces it to "not available".
        """
        requirements = super().requirements(simulation, atoms)
        target = self.reported_convention(simulation)
        if requirements.energy_convention is None:
            reasons = dict(requirements.reasons)
            reasons["energy_convention"] = (
                f"potential.energy_convention = {target!r} is the convention to report"
            )
            requirements = replace(requirements, energy_convention=target, reasons=reasons)
        if requirements.energy_convention == INTERACTION and not self.atomic_reference_energies():
            raise CapabilityError(
                f"{self.label} cannot report the interaction energy:",
                [
                    "MACECalculator returns the total energy (the model's E0s included); "
                    "the interaction energy is that minus the model's own E0s for the "
                    f"evaluated head, which are not known here: {self._e0_unavailable_reason()}. "
                    "Declared potential.atomic_reference_energies cannot stand in for them "
                    "(they are only cross-checked). Report energy_convention = 'total'."
                ],
            )
        return requirements

    def check_simulation(self, simulation: SimulationSpec, atoms=None) -> None:
        """Model/config conflicts, then the ASE integrator map."""
        self._raise_model_problems()
        ase_engine.check_simulation(simulation, atoms, options=self.engine.options)

    def _raise_model_problems(self) -> None:
        problems = self.adapter.model_problems()
        if problems:
            raise ConfigError(
                f"{self.potential.label}: the configuration disagrees with the model file "
                f"{self.potential.model_path}:\n" + "\n".join(f"  - {p}" for p in problems)
            )

    def _device_problem(self) -> str | None:
        """Why ``potential.device`` cannot be used here, or ``None``.

        Asked only once the model route is otherwise available (torch is
        installed and the file exists), so it never imports torch on a
        machine that cannot run this route anyway.
        """
        device = self.potential.device
        kind, _, ordinal = device.partition(":")
        if kind == "cpu":
            return None
        import torch

        if kind == "cuda":
            if not torch.cuda.is_available():
                build = (
                    "a CPU-only build"
                    if getattr(torch.version, "cuda", None) is None
                    else f"built for CUDA {torch.version.cuda}"
                )
                return (
                    f"potential.device = {device!r}, but torch {torch.__version__} ({build}) "
                    "sees no CUDA device here (torch.cuda.is_available() is False)"
                )
            count = torch.cuda.device_count()
            index = int(ordinal or 0)
            if index >= count:
                return (
                    f"potential.device = {device!r}, but torch sees only {count} CUDA "
                    f"device(s) here (ordinals 0..{count - 1})"
                )
        elif kind == "mps":
            if not torch.backends.mps.is_available():
                return f"potential.device = {device!r}, but torch reports MPS unavailable"
        return None

    # -- the calculator ---------------------------------------------------

    def calculator_kwargs(self) -> dict:
        """The keyword arguments ``MACECalculator`` is constructed with."""
        try:
            head = self.adapter.resolved_head()
        except ConfigError:
            head = self.potential.head
        kwargs = {
            "model_paths": [str(self.potential.model_path)],
            "device": self.potential.device,
            "default_dtype": self.potential.precision,
            "enable_cueq": False,
        }
        if head is not None:
            kwargs["head"] = head
        if self.potential.compile_mode is not None:
            kwargs["compile_mode"] = self.potential.compile_mode
        return kwargs

    def calculator(self):
        """Build MACE's ASE calculator. The first torch import happens here.

        The device is pre-checked (through :meth:`availability`) before
        torch deserialises the model onto it.
        """
        if self._calculator is not None:
            return self._calculator
        self.require_available()
        self._raise_model_problems()
        import torch
        from mace.calculators import MACECalculator

        if self.engine.threads is not None:
            torch.set_num_threads(int(self.engine.threads))
        kwargs = self.calculator_kwargs()
        calculator = MACECalculator(**kwargs)
        self._runtime_report = self._read_back(calculator, kwargs)
        self._calculator = calculator
        return calculator

    def _read_back(self, calculator, kwargs: dict) -> dict:
        """What the live calculator actually evaluates with, read from it."""
        import mace
        import torch

        parameter = next(iter(calculator.models[0].parameters()), None)
        expected_device = torch.device(self.potential.device)
        actual_device = parameter.device if parameter is not None else None
        report = {
            "calculator": CALCULATOR,
            "kwargs_passed": {k: (list(v) if isinstance(v, list) else v) for k, v in kwargs.items()},
            "model_sha256": self.adapter.model_sha256(),
            "head": calculator.head,
            "available_heads": list(calculator.available_heads),
            "device_requested": self.potential.device,
            "device_effective": str(actual_device) if actual_device is not None else None,
            "dtype_requested": self.potential.precision,
            "dtype_model_native": self.adapter.model_dtype(),
            "dtype_effective": (
                str(parameter.dtype).replace("torch.", "") if parameter is not None else None
            ),
            "calculator_default_dtype": calculator.default_dtype,
            "compile_mode": self.potential.compile_mode,
            "use_compile": bool(getattr(calculator, "use_compile", False)),
            "enable_cueq": False,
            "r_max_angstrom": float(calculator.r_max),
            "torch_threads": int(torch.get_num_threads()),
            "torch_threads_source": (
                "engine.threads (torch.set_num_threads)"
                if self.engine.threads is not None
                else "torch default"
            ),
            "versions": {"torch": torch.__version__, "mace": getattr(mace, "__version__", None)},
        }
        mismatches = []
        if report["dtype_effective"] != self.potential.precision:
            mismatches.append(
                f"model parameters are {report['dtype_effective']}, not the requested "
                f"{self.potential.precision}"
            )
        if actual_device is not None and (
            actual_device.type != expected_device.type
            or (expected_device.index is not None and actual_device.index != expected_device.index)
        ):
            mismatches.append(
                f"model parameters are on {actual_device}, not the requested "
                f"{self.potential.device}"
            )
        if kwargs.get("head") is not None and calculator.head != kwargs["head"]:
            mismatches.append(f"MACE evaluates head {calculator.head!r}, not {kwargs['head']!r}")
        if kwargs.get("compile_mode") is not None and not report["use_compile"]:
            mismatches.append(f"compile_mode {kwargs['compile_mode']!r} was not applied")
        if mismatches:
            raise MlipError(
                "MACECalculator did not apply the requested settings: " + "; ".join(mismatches)
            )
        native = report["dtype_model_native"]
        if native and native != report["dtype_effective"]:
            report["dtype_note"] = (
                f"the model was saved in {native} and is evaluated in "
                f"{report['dtype_effective']} (MACECalculator converted its parameters)"
            )
        return report

    # -- records ----------------------------------------------------------

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        try:
            head = self.adapter.resolved_head()
        except ConfigError as exc:
            head = {"refused": str(exc)}
        e0s = self.atomic_reference_energies()
        parameters = {
            "implementation": self.implementation,
            "calculator": CALCULATOR,
            "calculator_kwargs": self.calculator_kwargs(),
            "model_paths": [str(self.potential.model_path)],
            "model_sha256": self.adapter.model_sha256(),
            "device": self.potential.device,
            "default_dtype": self.potential.precision,
            "model_native_dtype": self.adapter.model_dtype(),
            "head": head,
            "compile_mode": self.potential.compile_mode,
            "enable_cueq": False,
            "threads": self.engine.threads,
            "native_units": ASE.name,
            "native_energy_convention": TOTAL,
            "energy_convention": self.reported_convention(simulation),
            "atomic_reference_energies_eV": e0s,
            "atomic_reference_energies_source": (
                self.adapter.atomic_reference_energies_source() if e0s else None
            ),
            "forces": "raw calculator forces (constraints are not applied to them)",
            "runtime": dict(self._runtime_report) if self._runtime_report else None,
        }
        if simulation is not None and simulation.task == "md":
            parameters["integrator"] = ase_engine.integrator_parameters(
                simulation, self.engine.options
            )
        return parameters

    def execution_plan(self, simulation: SimulationSpec, atoms=None) -> dict:
        """What would run, resolved without building the calculator."""
        plan = ase_engine.base_execution_plan(self, simulation, atoms)
        sha = self.adapter.model_sha256()
        available = plan["availability"]["available"]
        target = self.reported_convention(simulation)
        e0s = self.atomic_reference_energies()
        native_dtype = self.adapter.model_dtype()
        precision = self.potential.precision
        try:
            head = self.adapter.resolved_head()
        except ConfigError:
            head = None  # the refusal is already in unmet_capabilities
        plan.update(
            model_checkpoint={"path": str(self.potential.model_path), "sha256": sha},
            energy_convention=target,
            energy_convention_native=TOTAL,
            energy_conversion=(
                None
                if target == TOTAL
                else "interaction = MACECalculator total - sum of the model's E0s"
                + ("" if e0s else " (refused: model E0s unknown)")
            ),
            atomic_reference_energies_eV=e0s,
            head=head,
            calculator={"class": CALCULATOR, "kwargs": self.calculator_kwargs()},
            device={
                "requested": self.potential.device,
                "effective": self.potential.device if available else None,
                "guaranteed": bool(available),
                "note": (
                    "MACECalculator moves the model to this device (model.to(device)); "
                    "the parameters' device is read back after the build and a mismatch "
                    "is an error"
                    if available
                    else f"not runnable here: {plan['availability']['detail']}"
                ),
            },
            precision={
                "requested": precision,
                "effective": precision if available else None,
                "guaranteed": bool(available),
                "model_native": native_dtype,
                "note": (
                    f"MACECalculator(default_dtype={precision!r}) converts the parameters "
                    "and builds its input tensors in this dtype; read back after the build"
                    + (
                        f". The model was saved in {native_dtype}, so it is converted"
                        if native_dtype and native_dtype != precision
                        else ""
                    )
                ),
            },
            threads={
                "requested": self.engine.threads,
                "applied_with": "torch.set_num_threads" if self.engine.threads else None,
            },
            model_hashes={
                "declared": dict(self.potential.declared_hashes()),
                "observed": {str(self.potential.model_path): sha} if sha else {},
            },
        )
        return plan

    # -- execution --------------------------------------------------------

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        self.validate(simulation, atoms)
        calculator = self.calculator()
        payload = ase_engine.singlepoint(
            atoms,
            calculator,
            compute_stress=simulation.compute_stress or simulation.ensemble == "npt",
            compute_per_atom_energy=simulation.compute_per_atom_energy,
        )
        result = self._result(
            payload,
            energy_convention=TOTAL,
            native_units=ASE.name,
            **payload["native"],
            runtime=dict(self._runtime_report or {}),
        )
        return ase_engine.restate(
            result,
            self.reported_convention(simulation),
            atomic_reference_energies=self.atomic_reference_energies(),
            source=self.adapter.atomic_reference_energies_source() or "model file",
        )

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        self.validate(simulation, atoms)
        endpoint = self.as_singlepoint(simulation)
        initial = self.singlepoint(atoms, endpoint)
        payload = ase_engine.run_md(
            atoms,
            self.calculator(),
            simulation,
            workdir=Path(workdir),
            options=self.engine.options,
        )
        final = self.singlepoint(payload["atoms"], endpoint)
        trajectory = self._trajectory(payload, simulation, initial, final)
        extras = {**trajectory.extras, "runtime": dict(self._runtime_report or {})}
        return replace(trajectory, extras=extras)


__all__ = ["MaceAseBridge", "CALCULATOR"]
