"""The analytic test potential on ASE. Diagnostic infrastructure, not science.

Exists so the whole pipeline -- configuration parsing, bridge resolution,
capability negotiation, element-coverage validation, unit and convention
handling, provenance and both execution commands -- is exercised end to end in
ordinary CI, on a machine with no torch, no LAMMPS and no OpenMM.

The mock reports the energy convention its configuration declares (there are
no real self-energies in a Lennard-Jones form); a job asking for the other
convention gets it converted with the declared ``atomic_reference_energies``,
and is refused during negotiation when none are declared.
"""
from __future__ import annotations

from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..errors import ConfigError
from ..potentials.mock import MOCK_WARNING, MockAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import ASE
from ..engines import ase_engine
from ..engines.ase_engine import AseEngine
from .base import Bridge, complete_plan


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

    def check_simulation(self, simulation: SimulationSpec, atoms=None) -> None:
        """The ASE integrator map, plus: the analytic calculator has no threads knob."""
        if self.engine.threads is not None:
            raise ConfigError(
                "engine.threads cannot be applied on the mock/ASE route: the analytic "
                "calculator is single-threaded numpy. Remove engine.threads."
            )
        ase_engine.check_simulation(simulation, atoms, options=self.engine.options)

    def calculator(self):
        if self._calculator is None:
            self.require_available()
            self._calculator = self.adapter.calculator()
        return self._calculator

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        parameters = {
            "implementation": self.implementation,
            "calculator": "nio_md_prep.mlip._mock_calculator.MockLennardJones",
            "epsilon_eV": self.potential.epsilon_eV,
            "sigma_angstrom": self.potential.sigma_angstrom,
            "cutoff_angstrom": self.potential.cutoff_angstrom,
            "native_units": ASE.name,
            "native_energy_convention": self.potential.energy_convention,
            "energy_convention": self._reported_convention(simulation),
            "forces": "raw calculator forces (constraints are not applied to them)",
            "warning": MOCK_WARNING,
        }
        if simulation is not None and simulation.task == "md":
            parameters["integrator"] = ase_engine.integrator_parameters(
                simulation, self.engine.options
            )
        return parameters

    def execution_plan(self, simulation: SimulationSpec, atoms=None) -> dict:
        """What would run, resolved without building the calculator."""
        plan = ase_engine.base_execution_plan(self, simulation, atoms)
        target = self._reported_convention(simulation)
        native = self.potential.energy_convention
        plan.update(
            energy_convention=target,
            energy_convention_native=native,
            energy_conversion=(
                None
                if target == native
                else f"{target} = {native} {'-' if target == 'interaction' else '+'} sum of "
                "the declared potential.atomic_reference_energies"
            ),
            calculator={
                "class": "nio_md_prep.mlip._mock_calculator.MockLennardJones",
                "kwargs": {
                    "epsilon_eV": self.potential.epsilon_eV,
                    "sigma_angstrom": self.potential.sigma_angstrom,
                    "cutoff_angstrom": self.potential.cutoff_angstrom,
                    "supported_elements": list(self.potential.elements),
                },
            },
            device={
                "requested": None,
                "effective": "cpu",
                "guaranteed": True,
                "note": "the analytic calculator is numpy on the CPU; there is no device knob",
            },
            precision={
                "requested": self.engine.precision,
                "effective": "float64",
                "guaranteed": True,
                "note": "numpy float64 throughout",
            },
            warning=MOCK_WARNING,
        )
        return complete_plan(plan)

    def _reported_convention(self, simulation: SimulationSpec | None) -> str:
        if simulation is not None and simulation.energy_convention:
            return simulation.energy_convention
        return self.potential.energy_convention

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        payload = ase_engine.singlepoint(
            atoms,
            self.calculator(),
            compute_stress=simulation.compute_stress or simulation.ensemble == "npt",
            compute_per_atom_energy=simulation.compute_per_atom_energy,
        )
        result = self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=ASE.name,
        )
        return ase_engine.restate(
            result,
            self._reported_convention(simulation),
            atomic_reference_energies=self.atomic_reference_energies(),
            source="potential.atomic_reference_energies (declared)",
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
        return self._trajectory(payload, simulation, initial, final)


__all__ = ["MockAseBridge"]
