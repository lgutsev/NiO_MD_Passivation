"""MACE on OpenMM, through OpenMM-ML.

The energy convention is the whole story here. OpenMM-ML's MACE
implementation distinguishes the interaction energy from the energy including
atomic self-energies, and *interaction* is its default. A comparison against
the ASE calculator that ignores this looks catastrophically wrong -- thousands
of eV apart on a few thousand atoms -- while the forces agree to machine
precision, because the difference is a per-element constant with zero
gradient.

So this bridge:

* asks OpenMM-ML for whichever convention the job declares, rather than
  taking the default silently;
* records the convention that ran in the manifest and on every result;
* refuses to advertise a convention it cannot reach, which for a model with
  no known atomic reference energies means it can only report the one
  OpenMM-ML gives it.

OpenMM reports no stress tensor through this route, so ``stress`` is false and
a constant-pressure job is rejected during capability negotiation.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..errors import EnergyConventionError
from ..potentials.mace import MaceAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import INTERACTION, OPENMM
from ..engines import openmm_engine
from ..engines.openmm_engine import ENERGY_TYPE, OpenMMEngine
from .base import Bridge

#: What OpenMM-ML gives you if nobody says otherwise.
OPENMM_ML_DEFAULT_CONVENTION = INTERACTION


class MaceOpenMMBridge(Bridge):
    potential_kind = "mace"
    engine_kind = "openmm"
    implementation = "openmm-ml"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = MaceAdapter(potential)
        self.runtime = OpenMMEngine(engine)
        self._systems: dict[str, object] = {}

    # -- conventions ------------------------------------------------------

    def requested_convention(self, simulation: SimulationSpec | None = None) -> str:
        """Which convention this run will ask OpenMM-ML for.

        The job's explicit request wins; then the potential's own convention,
        if this route can produce it; otherwise OpenMM-ML's default, which is
        the interaction energy.
        """
        if simulation is not None and simulation.energy_convention:
            return simulation.energy_convention
        if self.potential.energy_convention in ENERGY_TYPE:
            return self.potential.energy_convention
        return OPENMM_ML_DEFAULT_CONVENTION

    def native_convention(self) -> str:
        return self.requested_convention()

    def capabilities(self) -> CapabilitySet:
        engine_capabilities = replace(
            self.runtime.capabilities(),
            gpu=self.runtime.capabilities().gpu and self.potential.wants_gpu,
        )
        route = self._route_capabilities(
            engine_capabilities,
            energy_convention=self.native_convention(),
            native_units=OPENMM.name,
            notes=(
                f"OpenMM-ML returnEnergyType={ENERGY_TYPE[self.native_convention()]!r}",
                "OpenMM-ML's own default is the interaction energy; this route states "
                "the convention explicitly instead of inheriting it",
            ),
        )
        # OpenMM-ML can produce either convention directly, whether or not this
        # subsystem knows the E0s -- so both are reachable here even when a
        # post-hoc conversion would not be.
        route = replace(
            route, convertible_energy_conventions=frozenset(ENERGY_TYPE)
        )
        return self.adapter.capabilities().intersect(route)

    def availability(self) -> Availability:
        model = self.adapter.availability()
        if not model:
            return model
        return self.runtime.availability()

    def atomic_reference_energies(self):
        return self.adapter.atomic_reference_energies()

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        convention = self.native_convention()
        return {
            "implementation": self.implementation,
            "potential": "openmmml.MLPotential('mace')",
            "model_path": str(self.potential.model_path),
            "returnEnergyType": ENERGY_TYPE[convention],
            "energy_convention": convention,
            "platform": self.engine.platform,
            "precision": self.engine.precision,
            "native_units": OPENMM.name,
            "native_energy_unit": OPENMM.energy,
            "native_length_unit": OPENMM.length,
            "atomic_reference_energies_known": bool(self.atomic_reference_energies()),
        }

    # -- execution --------------------------------------------------------

    def system(self, atoms, simulation: SimulationSpec | None = None):
        """Build (and cache) the OpenMM System for one energy convention.

        Keyed by convention, because asking OpenMM-ML for the interaction
        energy and asking it for the total energy produce different Systems --
        a single cache would silently hand back the wrong one.
        """
        convention = self.requested_convention(simulation)
        if convention not in self._systems:
            self.require_available()
            system, _ = openmm_engine.build_system(
                atoms, self.potential.model_path, energy_convention=convention
            )
            self._systems[convention] = system
        return self._systems[convention]

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        self.validate(simulation, atoms)
        convention = self.requested_convention(simulation)
        if convention not in ENERGY_TYPE:
            raise EnergyConventionError(
                f"OpenMM-ML cannot report the {convention!r} energy convention"
            )
        payload = openmm_engine.singlepoint(
            atoms, self.system(atoms, simulation), self.engine
        )
        return self._result(
            payload, energy_convention=convention, native_units=OPENMM.name
        )

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        self.validate(simulation, atoms)
        endpoint = self.as_singlepoint(simulation)
        initial = self.singlepoint(atoms, endpoint)
        payload = openmm_engine.run_md(
            atoms,
            self.system(atoms, simulation),
            self.engine,
            simulation,
            workdir=Path(workdir),
        )
        final = self.singlepoint(payload["atoms"], endpoint)
        return self._trajectory(payload, simulation, initial, final)


__all__ = ["MaceOpenMMBridge", "OPENMM_ML_DEFAULT_CONVENTION"]
