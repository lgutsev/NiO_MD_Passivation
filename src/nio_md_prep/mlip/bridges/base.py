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
plausible number for a structure the model was never trained on. Geometry
facts (per-axis PBC, frozen atoms, vacuum gaps) enter through
:meth:`Bridge.requirements`; engine-specific refusals (an unsupported
thermostat, a platform property the platform lacks) belong in
:meth:`Bridge.check_simulation`, which both ``validate`` paths call.
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

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        """Engine-specific parameters, recorded verbatim in the manifest.

        ``simulation`` is the job's simulation spec when one exists (the
        manifest writer and ``validate`` always pass it), so a route can
        record what it will actually request -- e.g. the energy type implied
        by ``simulation.energy_convention`` -- rather than a default.
        """
        capabilities = self.capabilities()
        return {
            "implementation": self.implementation,
            "native_units": capabilities.native_units,
            "energy_convention": capabilities.native_energy_convention,
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
        """What this job needs from this route, including the structure's geometry.

        With a structure, the per-axis PBC, the ``FixAtoms`` indices (any
        other ASE constraint is refused here, by name) and -- for a fully
        periodic NPT job -- the vacuum gaps along each lattice direction are
        passed to :meth:`SimulationSpec.required_capabilities`, which refuses
        geometry the request cannot support. ``engine.precision``, when set,
        becomes a precision requirement, and the engine kind an engine
        requirement, so both are negotiated rather than ignored.
        """
        from dataclasses import replace

        if atoms is None:
            requirements = simulation.required_capabilities()
        else:
            from ..structures import fixed_atom_indices, periodic_axes, vacuum_gaps

            pbc = periodic_axes(atoms)
            gaps = {}
            if simulation.ensemble == "npt" and all(pbc):
                gaps = vacuum_gaps(
                    atoms,
                    threshold_angstrom=simulation.resolved_vacuum_gap_threshold_angstrom,
                )
            requirements = simulation.required_capabilities(
                pbc=pbc,
                fixed_atoms=len(fixed_atom_indices(atoms)),
                vacuum_gaps=gaps,
            )
        reasons = dict(requirements.reasons)
        precision = requirements.precision
        if self.engine.precision is not None:
            precision = self.engine.precision
            reasons["precision"] = (
                "engine.precision is a cross-check and must match what the route "
                "evaluates in (potential.precision)"
            )
        return replace(
            requirements, engine=self.engine_kind, precision=precision, reasons=reasons
        )

    def check_simulation(self, simulation: SimulationSpec, atoms=None) -> None:
        """Engine-specific refusals, before anything is constructed. Default: none.

        Overridden by bridges whose engine supports only part of what a
        :class:`SimulationSpec` can express -- a thermostat or barostat
        coupling the engine does not implement, a platform property the
        platform lacks, ``engine.threads`` on a route that cannot apply it,
        a structure outside the engine's box conventions. Raise
        :class:`~nio_md_prep.mlip.errors.ConfigError` (a bad request) or
        :class:`~nio_md_prep.mlip.errors.CapabilityError` (a route limit).
        Called by :meth:`validate` and by ``jobs.validate_job``; must not
        execute the engine or import a heavy backend.
        """

    def validate(self, simulation: SimulationSpec, atoms=None) -> CapabilitySet:
        """Refuse an impossible job before constructing anything.

        Element coverage first, because "this model has never seen fluorine"
        is a more useful message than "this route cannot report a stress";
        then capability negotiation; then :meth:`check_simulation`.
        """
        capabilities = self.capabilities()
        if atoms is not None:
            check_elements(
                capabilities, atoms.get_chemical_symbols(), label=self.potential.label
            )
        negotiate(capabilities, self.requirements(simulation, atoms), label=self.label)
        self.check_simulation(simulation, atoms)
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
        """Build a validated :class:`PotentialResult` from an engine payload.

        Single points never carry a trajectory: any ``trajectory_path`` in a
        single-point payload is ignored here, and engines should set it to
        ``None``.
        """
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
        """Build a :class:`TrajectoryResult` from an engine's MD payload.

        The payload contract (every key optional unless stated; absent means
        "not reported", recorded as ``None``, never guessed):

        ``steps_completed``
            The engine's own step counter after the run.
        ``frames_written``
            Frames counted by reading the trajectory back. The legacy
            ``frames`` key is *not* trusted as verified; it is kept in
            ``extras["frames_reported_unverified"]``.
        ``trajectory_path`` / ``log_path``
        ``final_positions`` / ``final_cell`` / ``final_pbc``
            Final geometry in the source basis (rotated back from any
            engine-internal frame).
        ``integrator_resolved``
            Mapping of what actually integrated; ``integrator`` (a name) is
            the legacy spelling and fills ``{"name": ...}`` when the mapping is
            absent.
        ``energy_series``
            ``{"time_fs": [...], "total_energy_eV": [...]}`` plus optional
            ``"conserved_energy_eV"`` and ``"conserved_quantity"``, sampled at
            unique, increasing times. Turned into ``diagnostics`` by
            :func:`nio_md_prep.mlip.diagnostics.series_diagnostics`, so every
            engine's drift means the same thing.
        ``temperature_start_K`` / ``temperature_end_K`` / ``max_temperature_K``
        / ``temperature_ndof`` / ``constraints`` / ``wall_time_s``
        """
        from ..diagnostics import series_diagnostics

        series = payload.get("energy_series")
        if series is not None:
            diagnostics = series_diagnostics(simulation.ensemble, series, initial.n_atoms)
        else:
            diagnostics = {
                "available": False,
                "reason": "the engine reported no energy time series for this run",
            }
        integrator = payload.get("integrator_resolved")
        if integrator is None:
            integrator = (
                {"name": payload["integrator"]} if payload.get("integrator") else {}
            )
        extras: dict = {"integrator": payload.get("integrator")}
        if "frames" in payload and "frames_written" not in payload:
            extras["frames_reported_unverified"] = payload["frames"]
        return TrajectoryResult(
            steps=simulation.steps,
            timestep_fs=simulation.timestep_fs,
            ensemble=simulation.ensemble,
            trajectory_path=payload.get("trajectory_path"),
            log_path=payload.get("log_path"),
            initial=initial,
            final=final,
            steps_completed=payload.get("steps_completed"),
            frames_written=payload.get("frames_written"),
            final_positions_angstrom=payload.get("final_positions"),
            final_cell_angstrom=payload.get("final_cell"),
            final_pbc=payload.get("final_pbc"),
            integrator_resolved=dict(integrator),
            diagnostics=diagnostics,
            temperature_start_K=payload.get("temperature_start_K"),
            temperature_end_K=payload.get("temperature_end_K"),
            max_temperature_K=payload.get("max_temperature_K"),
            temperature_ndof=payload.get("temperature_ndof"),
            constraints=payload.get("constraints"),
            wall_time_s=payload.get("wall_time_s"),
            extras=extras,
        )


__all__ = ["Bridge"]
