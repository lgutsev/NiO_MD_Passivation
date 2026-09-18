"""MACE on ASE: the direct route, and the reference for cross-engine checks.

MACE ships its own ASE calculator, so this bridge is thin -- which is exactly
why it is the reference. Its energy convention is the model's own (total
energy, including atomic self-energies), its units are already canonical, and
it introduces no export step that could change a number. The other two MACE
routes are compared against this one.
"""
from __future__ import annotations

from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..potentials.mace import MaceAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import ASE
from ..engines import ase_engine
from ..engines.ase_engine import AseEngine
from .base import Bridge


class MaceAseBridge(Bridge):
    potential_kind = "mace"
    engine_kind = "ase"
    implementation = "mace-ase-calculator"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = MaceAdapter(potential)
        self.runtime = AseEngine(engine)
        self._calculator = None

    def capabilities(self) -> CapabilitySet:
        return self.adapter.capabilities().intersect(
            self._route_capabilities(
                self.runtime.capabilities(),
                energy_convention=self.potential.energy_convention,
                native_units=ASE.name,
            )
        )

    def availability(self) -> Availability:
        model = self.adapter.availability()
        if not model:
            return model
        return self.runtime.availability()

    def atomic_reference_energies(self):
        return self.adapter.atomic_reference_energies()

    def calculator(self):
        """Build MACE's ASE calculator. The first torch import happens here."""
        if self._calculator is not None:
            return self._calculator
        self.require_available()
        from mace.calculators import MACECalculator

        self._calculator = MACECalculator(
            model_paths=[str(self.potential.model_path)],
            device=self.potential.device,
            default_dtype=self.potential.precision,
        )
        return self._calculator

    def engine_parameters(self) -> dict:
        return {
            "implementation": self.implementation,
            "calculator": "mace.calculators.MACECalculator",
            "model_paths": [str(self.potential.model_path)],
            "device": self.potential.device,
            "default_dtype": self.potential.precision,
            "native_units": ASE.name,
            "energy_convention": self.potential.energy_convention,
            "atomic_reference_energies_known": bool(self.atomic_reference_energies()),
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
        initial = self.singlepoint(atoms, simulation)
        payload = ase_engine.run_md(
            atoms, self.calculator(), simulation, workdir=Path(workdir)
        )
        final = self.singlepoint(payload["atoms"], simulation)
        return self._trajectory(payload, simulation, initial, final)


__all__ = ["MaceAseBridge"]
