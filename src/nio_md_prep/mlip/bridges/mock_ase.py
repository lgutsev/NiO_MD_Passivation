"""The analytic test potential on ASE. Diagnostic infrastructure, not science.

Exists so the whole pipeline -- configuration parsing, bridge resolution,
capability negotiation, element-coverage validation, unit and convention
handling, provenance and both execution commands -- is exercised end to end in
ordinary CI, on a machine with no torch, no LAMMPS and no OpenMM.
"""
from __future__ import annotations

from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..potentials.mock import MOCK_WARNING, MockAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import ASE
from ..engines import ase_engine
from ..engines.ase_engine import AseEngine
from .base import Bridge


class MockAseBridge(Bridge):
    potential_kind = "mock"
    engine_kind = "ase"
    implementation = "mock-ase"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = MockAdapter(potential)
        self.runtime = AseEngine(engine)
        self._calculator = None

    def capabilities(self) -> CapabilitySet:
        return self.adapter.capabilities().intersect(
            self._route_capabilities(
                self.runtime.capabilities(),
                energy_convention=self.potential.energy_convention,
                native_units=ASE.name,
                notes=(MOCK_WARNING,),
            )
        )

    def availability(self) -> Availability:
        model = self.adapter.availability()
        return model if not model else self.runtime.availability()

    def calculator(self):
        if self._calculator is None:
            self.require_available()
            self._calculator = self.adapter.calculator()
        return self._calculator

    def engine_parameters(self) -> dict:
        return {
            "implementation": self.implementation,
            "calculator": "nio_md_prep.mlip._mock_calculator.MockLennardJones",
            "epsilon_eV": self.potential.epsilon_eV,
            "sigma_angstrom": self.potential.sigma_angstrom,
            "cutoff_angstrom": self.potential.cutoff_angstrom,
            "native_units": ASE.name,
            "energy_convention": self.potential.energy_convention,
            "warning": MOCK_WARNING,
        }

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        payload = ase_engine.singlepoint(
            atoms,
            self.calculator(),
            compute_stress=simulation.compute_stress or simulation.ensemble == "npt",
            compute_per_atom_energy=simulation.compute_per_atom_energy,
        )
        return self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=ASE.name,
        )

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        self.validate(simulation, atoms)
        endpoint = self.as_singlepoint(simulation)
        initial = self.singlepoint(atoms, endpoint)
        payload = ase_engine.run_md(
            atoms, self.calculator(), simulation, workdir=Path(workdir)
        )
        final = self.singlepoint(payload["atoms"], endpoint)
        return self._trajectory(payload, simulation, initial, final)


__all__ = ["MockAseBridge"]
