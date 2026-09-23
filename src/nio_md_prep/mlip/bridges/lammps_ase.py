"""A LAMMPS-native MLIP driven from ASE.

LAMMPS still evaluates the potential -- this is not a reimplementation. ASE's
``LAMMPSlib`` calculator starts a LAMMPS instance in-process, feeds it the
pair commands, and drives it, which is what makes ASE's optimisers,
integrators and analysis available to a pair style that only exists inside
LAMMPS.

The ``type_map`` is what makes the bridge possible: it maps an ASE ``Atoms``
object's elements onto the LAMMPS types the ``pair_coeff`` line expects.
Without it the two would agree only by accident.

Note what this bridge does *not* imply. Driving a LAMMPS pair style from ASE
works because ASE can call into LAMMPS. There is no equivalent route into
OpenMM, and none is faked: see the explicitly-unsupported ``lammps -> openmm``
cell in :mod:`nio_md_prep.mlip.registry`.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, module_available
from ..potentials.lammps_mlip import LammpsMlipAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import lammps_unit_system
from ..engines import ase_engine
from ..engines.ase_engine import AseEngine
from .base import Bridge


class LammpsAseBridge(Bridge):
    potential_kind = "lammps"
    engine_kind = "ase"
    implementation = "ase-lammps"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = LammpsMlipAdapter(potential)
        self.runtime = AseEngine(engine)
        self._calculator = None

    def capabilities(self) -> CapabilitySet:
        """ASE's LAMMPS calculator carries energy, forces and stress.

        Not a per-atom energy decomposition: ``LAMMPSlib`` does not expose
        ``pe/atom`` as an ASE property, so a job asking for site energies is
        pointed at the native LAMMPS route instead of receiving nothing.
        """
        engine_capabilities = replace(
            self.runtime.capabilities(),
            per_atom_energy=False,
            notes=self.runtime.capabilities().notes
            + (
                "ASE's LAMMPSlib calculator does not expose per-atom energies; "
                "use engine = 'lammps' for a site-energy decomposition",
            ),
        )
        return self.adapter.capabilities().intersect(
            self._route_capabilities(
                engine_capabilities,
                energy_convention=self.potential.energy_convention,
                # LAMMPSlib returns ASE-convention eV/Angstrom regardless of the
                # LAMMPS unit style, so the canonical system is what comes out.
                native_units="ase",
            )
        )

    def availability(self) -> Availability:
        missing = self.adapter.missing_model_files()
        if missing:
            return Availability(
                False, (), f"model file(s) not found: {', '.join(missing)}"
            )
        ase_available = self.runtime.availability()
        if not ase_available:
            return ase_available
        if not module_available("lammps"):
            return Availability(
                False,
                ("lammps",),
                "ASE's LAMMPSlib calculator needs the LAMMPS python module",
            )
        return Availability(True, (), "ASE driving an in-process LAMMPS")

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        return {
            "implementation": self.implementation,
            "calculator": "ase.calculators.lammpslib.LAMMPSlib",
            "native_units": "ase",
            "energy_convention": self.potential.energy_convention,
            "lammps_units": self.potential.units,
            "lammps_unit_system": lammps_unit_system(self.potential.units).name,
            "atom_types": {v: k for k, v in sorted(self.potential.type_map.items())},
            "framework": self.adapter.framework,
            "required_packages": list(self.adapter.required_packages()),
            # The exact strings handed to the in-process LAMMPS, verbatim.
            "pair_commands": list(self.potential.render_pair_commands()),
        }

    def calculator(self):
        """Build ASE's LAMMPSlib calculator around this pair style."""
        if self._calculator is not None:
            return self._calculator
        self.require_available()
        from ase.calculators.lammpslib import LAMMPSlib

        commands = list(self.potential.render_pair_commands())
        self._calculator = LAMMPSlib(
            lmpcmds=commands,
            atom_types={
                symbol: type_id
                for type_id, symbol in sorted(self.potential.type_map.items())
            },
            lammps_header=[
                f"units {self.potential.units}",
                f"atom_style {self.potential.atom_style}",
                "atom_modify map array sort 0 0",
            ],
            keep_alive=True,
            log_file=str(self.engine.options.get("log_file", "lammpslib.log")),
        )
        return self._calculator

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        payload = ase_engine.singlepoint(
            atoms,
            self.calculator(),
            compute_stress=simulation.compute_stress or simulation.ensemble == "npt",
            compute_per_atom_energy=False,
        )
        return self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units="ase",
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


__all__ = ["LammpsAseBridge"]
