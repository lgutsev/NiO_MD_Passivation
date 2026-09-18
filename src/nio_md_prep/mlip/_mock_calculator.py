"""The ASE calculator behind the ``mock`` potential. Testing infrastructure.

Kept out of ``potentials/`` because it imports ASE and numpy at module level;
``potentials.mock`` imports it lazily, so the spec and capability layers stay
import-free.

The form is a Lennard-Jones pair potential shifted to zero at the cutoff:

    E = sum_{i<j} 4 eps [ (sigma/r)^12 - (sigma/r)^6 ] - E_shift ,  r < r_cut

Conservative and analytically differentiable, which is what makes an NVE
energy-drift smoke test meaningful. The virial is accumulated pairwise, so
the stress tensor is exact for this form rather than a finite difference.
"""
from __future__ import annotations

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.neighborlist import NeighborList


class MockLennardJones(Calculator):
    """A deterministic, periodic-aware analytic calculator. Not a trained model."""

    implemented_properties = ["energy", "energies", "free_energy", "forces", "stress"]

    def __init__(
        self,
        *,
        epsilon_eV: float = 0.01,
        sigma_angstrom: float = 2.5,
        cutoff_angstrom: float = 6.0,
        supported_elements=(),
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.epsilon = float(epsilon_eV)
        self.sigma = float(sigma_angstrom)
        self.cutoff = float(cutoff_angstrom)
        self.supported_elements = tuple(supported_elements)
        self._shift = self._pair_energy(np.array([self.cutoff]))[0]

    def _pair_energy(self, r):
        c6 = (self.sigma / r) ** 6
        return 4.0 * self.epsilon * (c6 * c6 - c6)

    def _pair_energy_and_derivative(self, r):
        c6 = (self.sigma / r) ** 6
        energy = 4.0 * self.epsilon * (c6 * c6 - c6) - self._shift
        # dE/dr
        derivative = 4.0 * self.epsilon * (-12.0 * c6 * c6 + 6.0 * c6) / r
        return energy, derivative

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        atoms = self.atoms
        n_atoms = len(atoms)
        site_energies = np.zeros(n_atoms)
        forces = np.zeros((n_atoms, 3))
        virial = np.zeros((3, 3))

        neighbors = NeighborList(
            [self.cutoff / 2.0] * n_atoms, skin=0.0, self_interaction=False, bothways=False
        )
        neighbors.update(atoms)
        positions = atoms.get_positions()
        cell = atoms.get_cell()

        for i in range(n_atoms):
            indices, offsets = neighbors.get_neighbors(i)
            if len(indices) == 0:
                continue
            vectors = positions[indices] + offsets @ cell - positions[i]
            distances = np.linalg.norm(vectors, axis=1)
            inside = (distances < self.cutoff) & (distances > 1e-8)
            if not np.any(inside):
                continue
            vectors = vectors[inside]
            distances = distances[inside]
            targets = indices[inside]
            energy, derivative = self._pair_energy_and_derivative(distances)
            # Half the pair energy to each partner, so site energies sum to E.
            np.add.at(site_energies, targets, 0.5 * energy)
            site_energies[i] += 0.5 * np.sum(energy)
            pair_forces = (derivative / distances)[:, None] * vectors
            forces[i] += np.sum(pair_forces, axis=0)
            np.add.at(forces, targets, -pair_forces)
            virial += vectors.T @ pair_forces

        self.results["energy"] = float(np.sum(site_energies))
        self.results["free_energy"] = self.results["energy"]
        self.results["energies"] = site_energies
        self.results["forces"] = forces
        if atoms.cell.rank == 3 and atoms.get_volume() > 0:
            stress = virial / atoms.get_volume()
            stress = 0.5 * (stress + stress.T)
            self.results["stress"] = np.array(
                [
                    stress[0, 0],
                    stress[1, 1],
                    stress[2, 2],
                    stress[1, 2],
                    stress[0, 2],
                    stress[0, 1],
                ]
            )


__all__ = ["MockLennardJones"]
