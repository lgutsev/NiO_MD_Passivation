"""What a bridge is: the one object that knows both halves.

A bridge owns exactly the knowledge that neither a potential nor an engine
has on its own -- how *this* model is handed to *that* engine, which units and
energy convention come back out, and what engine-specific parameters were
generated along the way. Everything a bridge does not need to know (what a
thermostat is, what a MACE cutoff is) stays in the engine and potential
modules.

:meth:`Bridge.validate` is the gate. It runs element-coverage validation and
capability negotiation before any engine object is constructed, which is what
makes ``mlip validate`` meaningful and what keeps a job from producing a
plausible number for a structure the model was never trained on.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import ClassVar

from ..capabilities import CapabilitySet, RequirementSet, check_elements, negotiate
from ..environment import Availability
from ..errors import MissingDependencyError, NotImplementedYetError
from ..results import PotentialResult, TrajectoryResult
from ..specs import EngineSpec, PotentialSpec, SimulationSpec, as_singlepoint
from ..structures import is_periodic


class Bridge(ABC):
    """Base class for every cell of the compatibility matrix."""

    potential_kind: ClassVar[str] = "abstract"
    engine_kind: ClassVar[str] = "abstract"
    implementation: ClassVar[str] = "abstract"

    def __init__(
        self,
        potential: PotentialSpec,
        engine: EngineSpec,
        *,
        registration=None,
    ) -> None:
        self.potential = potential
        self.engine = engine
        self.registration = registration

    # -- description ------------------------------------------------------

    @property
    def label(self) -> str:
        return f"{self.potential_kind} -> {self.engine_kind} ({self.implementation})"

    @abstractmethod
    def capabilities(self) -> CapabilitySet:
        """What this route can actually deliver, potential and engine combined."""

    @abstractmethod
    def availability(self) -> Availability:
        """Whether this route can run on this machine, probed without importing."""

    def engine_parameters(self) -> dict:
        """Engine-specific parameters, recorded verbatim in the manifest."""
        return {
            "implementation": self.implementation,
            "native_units": self.capabilities().native_units,
            "energy_convention": self.capabilities().native_energy_convention,
        }

    def atomic_reference_energies(self) -> dict[str, float] | None:
        """E0s for converting between energy conventions, if any are known."""
        declared = getattr(self.potential, "atomic_reference_energies", None)
        return dict(declared) if declared else None

    def _route_capabilities(
        self,
        engine_capabilities: CapabilitySet,
        *,
        energy_convention: str,
        native_units: str,
        notes: tuple[str, ...] = (),
    ) -> CapabilitySet:
        """Restate an engine's capabilities with *this route's* conventions.

        The engine module cannot know them: whether MACE-on-OpenMM reports an
        interaction energy is a fact about the pairing, not about OpenMM. A
        route can convert between conventions only when atomic reference
        energies are available, so that is derived here rather than assumed.
        """
        from dataclasses import replace

        convertible = (
            frozenset({"total", "interaction"})
            if self.atomic_reference_energies()
            else frozenset()
        )
        return replace(
            engine_capabilities,
            native_energy_convention=energy_convention,
            convertible_energy_conventions=convertible,
            native_units=native_units,
            notes=engine_capabilities.notes + tuple(notes),
        )

    # -- gating -----------------------------------------------------------

    #: The single-point form of an MD spec, used for a run's endpoints.
    as_singlepoint = staticmethod(as_singlepoint)

    def requirements(self, simulation: SimulationSpec, atoms=None) -> RequirementSet:
        periodic = bool(atoms is not None and is_periodic(atoms))
        return simulation.required_capabilities(periodic=periodic)

    def validate(self, simulation: SimulationSpec, atoms=None) -> CapabilitySet:
        """Refuse an impossible job before constructing anything.

        Element coverage first, because "this model has never seen fluorine"
        is a more useful message than "this route cannot report a stress".
        """
        capabilities = self.capabilities()
        if atoms is not None:
            check_elements(
                capabilities, atoms.get_chemical_symbols(), label=self.potential.label
            )
        negotiate(capabilities, self.requirements(simulation, atoms), label=self.label)
        return capabilities

    def require_available(self) -> None:
        availability = self.availability()
        if not availability:
            raise MissingDependencyError(
                self.label, availability.missing, hint=availability.detail
            )

    # -- execution --------------------------------------------------------

    @abstractmethod
    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        """Evaluate one geometry, returning canonical units and a stated convention."""

    def run_md(
        self, atoms, simulation: SimulationSpec, *, workdir: Path
    ) -> TrajectoryResult:
        """Run a short diagnostic trajectory. Not production dynamics."""
        raise NotImplementedYetError(
            f"{self.label} does not implement molecular dynamics in this release"
        )

    # -- helpers for subclasses -------------------------------------------

    def _result(self, payload: dict, *, energy_convention: str, native_units: str, **extras):
        return PotentialResult(
            energy_eV=payload["energy_eV"],
            forces_eV_per_A=payload["forces_eV_per_A"],
            symbols=payload["symbols"],
            stress_eV_per_A3=payload.get("stress_eV_per_A3"),
            per_atom_energy_eV=payload.get("per_atom_energy_eV"),
            wall_time_s=payload.get("wall_time_s"),
            energy_convention=energy_convention,
            engine=self.engine_kind,
            potential=self.potential_kind,
            implementation=self.implementation,
            native_units=native_units,
            extras=extras or payload.get("native", {}),
        )

    def _trajectory(self, payload: dict, simulation: SimulationSpec, initial, final):
        return TrajectoryResult(
            steps=simulation.steps,
            timestep_fs=simulation.timestep_fs,
            ensemble=simulation.ensemble,
            frames=payload.get("frames", 0),
            trajectory_path=payload.get("trajectory_path"),
            log_path=payload.get("log_path"),
            initial=initial,
            final=final,
            temperature_start_K=payload.get("temperature_start_K"),
            temperature_end_K=payload.get("temperature_end_K"),
            max_temperature_K=payload.get("max_temperature_K"),
            total_energy_drift_eV_per_atom=payload.get("total_energy_drift_eV_per_atom"),
            wall_time_s=payload.get("wall_time_s"),
            extras={"integrator": payload.get("integrator")},
        )


__all__ = ["Bridge"]
