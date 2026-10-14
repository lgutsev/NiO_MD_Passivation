"""A deterministic analytic potential, for testing this machinery only.

Not science, and it says so in every report it produces. It exists because
the valuable part of this subsystem -- configuration parsing, bridge
resolution, capability negotiation, unit and convention handling, provenance,
and the ``singlepoint`` and ``smoke-md`` code paths -- must be testable in
ordinary CI on a machine with no torch, no LAMMPS and no OpenMM.

The functional form is a shifted Lennard-Jones pair potential evaluated with
ASE's own neighbour list, which makes it periodic-aware, conservative (so
energy drift in an NVE smoke test is meaningful), and analytically
differentiable for a per-atom energy decomposition.
"""
from __future__ import annotations

from ..capabilities import CapabilitySet
from ..environment import Availability, probe
from ..errors import ConfigError
from ..specs import MockPotentialSpec
from .base import PotentialAdapter

MOCK_WARNING = (
    "the 'mock' potential is a diagnostic Lennard-Jones stand-in, not a trained "
    "model; never use it for production numbers"
)


class MockAdapter(PotentialAdapter):
    kind = "mock"
    requires = ("ase", "numpy")

    spec: MockPotentialSpec

    def capabilities(self) -> CapabilitySet:
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=True,
            per_atom_energy=True,
            periodic=True,
            gpu=False,
            elements=frozenset(self.spec.elements),
            precisions=frozenset({"float64"}),
            engines=frozenset({"ase"}),
            native_energy_convention=self.spec.energy_convention,
            convertible_energy_conventions=(
                frozenset({"total", "interaction"})
                if self.spec.atomic_reference_energies
                else frozenset()
            ),
            notes=(MOCK_WARNING,),
        )

    def availability(self) -> Availability:
        return probe(self.requires, detail="analytic test potential")

    def calculator(self):
        """An ASE calculator implementing the analytic form."""
        from .._mock_calculator import MockLennardJones

        return MockLennardJones(
            epsilon_eV=self.spec.epsilon_eV,
            sigma_angstrom=self.spec.sigma_angstrom,
            cutoff_angstrom=self.spec.cutoff_angstrom,
            supported_elements=self.spec.elements,
        )

    def inspect(self) -> dict:
        report = super().inspect()
        report.update(
            {
                "epsilon_eV": self.spec.epsilon_eV,
                "sigma_angstrom": self.spec.sigma_angstrom,
                "cutoff_angstrom": self.spec.cutoff_angstrom,
                "warning": MOCK_WARNING,
                "availability": self.availability().as_dict(),
                "discovered": {"functional_form": "shifted Lennard-Jones", "trained": False},
            }
        )
        return report


def build_adapter(spec: MockPotentialSpec) -> MockAdapter:
    if not isinstance(spec, MockPotentialSpec):
        raise ConfigError(f"expected a mock potential spec; got {type(spec).__name__}")
    return MockAdapter(spec)


__all__ = ["MOCK_WARNING", "MockAdapter", "build_adapter"]
