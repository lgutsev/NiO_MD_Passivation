"""A LAMMPS-native MLIP running in LAMMPS, as intended.

The simplest cell of the matrix and the one with the least to go wrong: the
pair style already lives inside LAMMPS, so this bridge renders the deck,
writes the structure in the potential's own type order, runs it, and converts
the results out of LAMMPS's units.

Nothing here knows which framework the pair style belongs to. DeepMD, ML-IAP,
PACE and a MACE pair style all arrive as the same
:class:`~nio_md_prep.mlip.specs.LammpsMlipPotentialSpec`.
"""
from __future__ import annotations

from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..potentials.lammps_mlip import LammpsMlipAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import lammps_unit_system
from ..engines import lammps_engine
from ..engines.lammps_engine import LammpsEngine
from .base import Bridge


class LammpsNativeBridge(Bridge):
    potential_kind = "lammps"
    engine_kind = "lammps"
    implementation = "native"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = LammpsMlipAdapter(potential)
        self.runtime = LammpsEngine(engine)

    def capabilities(self) -> CapabilitySet:
        return self.adapter.capabilities().intersect(
            self._route_capabilities(
                self.runtime.capabilities(),
                energy_convention=self.potential.energy_convention,
                native_units=self.potential.unit_system_name,
            )
        )

    def availability(self) -> Availability:
        missing = self.adapter.missing_model_files()
        if missing:
            return Availability(
                available=False,
                missing=(),
                detail=f"model file(s) not found: {', '.join(missing)}",
            )
        return self.runtime.availability()

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        return {
            "implementation": self.implementation,
            "native_units": self.potential.unit_system_name,
            "energy_convention": self.potential.energy_convention,
            "lammps_units": self.potential.units,
            "atom_style": self.potential.atom_style,
            "framework": self.adapter.framework,
            "required_packages": list(self.adapter.required_packages()),
            "type_map": {str(k): v for k, v in sorted(self.potential.type_map.items())},
            # The exact strings LAMMPS executes, preserved verbatim.
            "pair_commands": list(self.potential.render_pair_commands()),
        }

    def deck(self, atoms, simulation: SimulationSpec) -> tuple[str, ...]:
        """Render the input deck without running it. Used by ``mlip validate``."""
        from ..structures import is_periodic

        return lammps_engine.render_deck(
            self.potential,
            simulation,
            n_types=len(self.potential.type_map),
            periodic=is_periodic(atoms),
        )

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        workdir = Path(self.engine.options.get("workdir", ".")) / "lammps_singlepoint"
        payload = self._run(atoms, self.as_singlepoint(simulation), workdir)
        return self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=self.potential.unit_system_name,
        )

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        capabilities = self.validate(simulation, atoms)
        endpoint = self.as_singlepoint(simulation)
        initial = self.singlepoint(atoms, endpoint)
        payload = self._run(atoms, simulation, Path(workdir))
        final = self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=self.potential.unit_system_name,
        )
        interval = max(1, simulation.trajectory_interval)
        payload["frames"] = simulation.steps // interval + 1
        return self._trajectory(payload, simulation, initial, final)

    def _run(self, atoms, simulation: SimulationSpec, workdir: Path) -> dict:
        from ..structures import is_periodic

        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        lammps_engine.write_data_file(
            atoms, self.potential, workdir / lammps_engine.DATA_FILE
        )
        commands = lammps_engine.render_deck(
            self.potential,
            simulation,
            n_types=len(self.potential.type_map),
            periodic=is_periodic(atoms),
        )
        raw = lammps_engine.run_deck(
            commands,
            workdir=workdir,
            engine_spec=self.engine,
            n_atoms=len(atoms),
            want_per_atom=simulation.compute_per_atom_energy,
        )
        return _to_canonical(raw, self.potential, atoms, commands)


def _to_canonical(raw: dict, potential, atoms, commands) -> dict:
    """Convert LAMMPS output into canonical units, explicitly and once."""
    units = potential.units
    thermo = raw["thermo"]
    per_atom = raw.get("per_atom_energy")
    if per_atom is not None:
        per_atom = [lammps_engine.energy_to_eV(v, units) for v in per_atom]
    final_atoms = atoms.copy()
    if raw.get("positions"):
        final_atoms.set_positions(raw["positions"])
    return {
        "energy_eV": lammps_engine.energy_to_eV(thermo["pe"], units),
        "forces_eV_per_A": lammps_engine.forces_to_eV_per_A(raw["forces"], units),
        "stress_eV_per_A3": lammps_engine.stress_from_pressure(thermo, units),
        "per_atom_energy_eV": per_atom,
        "symbols": list(atoms.get_chemical_symbols()),
        "wall_time_s": raw.get("wall_time_s"),
        "atoms": final_atoms,
        "frames": 0,
        "trajectory_path": str(
            Path(raw["deck_path"]).parent / lammps_engine.TRAJECTORY_FILE
        ),
        "log_path": raw.get("log_path"),
        "temperature_start_K": None,
        "temperature_end_K": thermo.get("temp"),
        "max_temperature_K": None,
        "total_energy_drift_eV_per_atom": None,
        "integrator": "lammps fix mlip_integrate",
        "native": {
            "route": raw.get("route"),
            "thermo": thermo,
            "units": lammps_unit_system(units).name,
            "deck_path": raw.get("deck_path"),
            "pair_commands": list(potential.render_pair_commands()),
            "deck": list(commands),
        },
    }


__all__ = ["LammpsNativeBridge"]
