"""A LAMMPS-native MLIP running in LAMMPS, as intended.

The simplest cell of the matrix and the one with the least to go wrong: the
pair style already lives inside LAMMPS, so this bridge stages the model
files, writes the structure in the potential's own type order and frame,
renders the deck, runs it, and converts the results out of LAMMPS's units
and back into the source basis -- all through
:func:`nio_md_prep.mlip.engines.lammps_engine.run_job`.

Nothing here knows which framework the pair style belongs to. DeepMD, ML-IAP,
PACE and a MACE pair style all arrive as the same
:class:`~nio_md_prep.mlip.specs.LammpsMlipPotentialSpec`.
"""
from __future__ import annotations

from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..errors import ConfigError
from ..potentials.lammps_mlip import LammpsMlipAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..engines import lammps_engine
from ..engines.lammps_engine import LammpsEngine
from .base import Bridge, complete_plan


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

    def launch(self) -> lammps_engine.LammpsLaunch:
        return lammps_engine.resolve_launch(self.engine)

    def check_simulation(self, simulation: SimulationSpec, atoms=None) -> None:
        """Launch options, thermostat/barostat mapping and the geometry preflight.

        ``engine.precision`` is refused: a LAMMPS-native pair style evaluates
        in whatever precision its package and model were built with, which
        nothing here can read or change, so the cross-check could only ever
        be recorded, never honoured.
        """
        if self.engine.precision is not None:
            raise ConfigError(
                f"engine.precision = {self.engine.precision!r} cannot be applied or verified "
                f"for the LAMMPS-native pair style {self.potential.pair_style.split()[0]!r}: "
                "its precision is fixed by the LAMMPS package and model it was built with. "
                "Remove engine.precision."
            )
        lammps_engine.check_request(self.potential, self.engine, simulation, atoms)

    def execution_plan(self, simulation: SimulationSpec, atoms=None) -> dict:
        """What a run would execute, resolved without running (``mlip validate``)."""
        launch = None
        launch_error = None
        try:
            launch = self.launch()
        except ConfigError as exc:
            launch_error = str(exc)
        accelerator = launch.accelerator if launch is not None else {"gpu": False}
        files = lammps_engine.model_file_plan(
            self.potential.model_paths, self.potential.model_hashes
        )
        first = self.potential.model_paths[0] if self.potential.model_paths else None
        plan = {
            "potential_kind": self.potential_kind,
            "engine": self.engine_kind,
            "implementation": self.implementation,
            "model_checkpoint": {
                "path": str(first) if first is not None else None,
                "sha256": files["observed"].get(str(first)) if first is not None else None,
                "model_files": [
                    {"path": key, "sha256": value, "exists": value is not None}
                    for key, value in files["observed"].items()
                ],
            },
            "exported_model": None,
            "elements": list(self.potential.elements),
            "type_map": {str(k): v for k, v in sorted(self.potential.type_map.items())},
            "energy_convention": self.reported_energy_convention(simulation),
            "energy_convention_native": self.potential.energy_convention,
            "units": lammps_engine.units_plan(self.potential.units),
            "device": {
                "requested": None,
                "effective": (
                    f"gpu ({accelerator['gpu_via']}, selected by engine.lammps_args)"
                    if accelerator["gpu"]
                    else "cpu"
                ),
                "guaranteed": not accelerator["gpu"],
                "note": (
                    "a LAMMPS-native pair style has no device setting here; engine.lammps_args "
                    "select an accelerator, and whether the build has that backend is known "
                    "only when LAMMPS starts"
                    if accelerator["gpu"]
                    else "no accelerator switches: the pair style runs on the host CPU"
                ),
            },
            "precision": {
                "requested": self.engine.precision,
                "effective": None,
                "guaranteed": False,
                "note": (
                    "set by the pair style's package and model build; not readable before a "
                    "run (engine.precision is refused for this route)"
                ),
            },
            "dynamics": lammps_engine.dynamics_plan(
                simulation, units=self.potential.units, atoms=atoms, options=self.engine.options
            ),
            "lammps": lammps_engine.lammps_plan(self.potential, launch, atoms=atoms),
            "openmm": None,
            "model_hashes": files,
            "deck": lammps_engine.deck_plan(
                self.potential, simulation, atoms, options=self.engine.options
            ),
            **lammps_engine.bridge_plan_basics(self, simulation, atoms),
        }
        if launch_error is not None:
            plan["lammps"]["launch_refused"] = launch_error
        return complete_plan(plan)

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
            # The exact strings LAMMPS executes (model files renamed to their
            # staged names in the job directory), preserved verbatim.
            **lammps_engine.describe_request(self.potential, self.engine, simulation, self.launch()),
        }

    def deck(self, atoms, simulation: SimulationSpec) -> tuple[str, ...]:
        """Render the input deck for ``atoms`` without running it."""
        from ..structures import fixed_atom_indices

        geometry = lammps_engine.prepare_geometry(atoms)
        md = simulation.task == "md"
        fixed = fixed_atom_indices(atoms) if md else ()
        return lammps_engine.render_deck(
            lammps_engine.rewrite_model_tokens(self.potential),
            simulation,
            n_types=len(self.potential.type_map),
            pbc=geometry.pbc,
            fixed_ids=[i + 1 for i in fixed],
            n_atoms=len(atoms),
            vacuum=(
                lammps_engine.vacuum_axes(atoms, simulation)
                if md and simulation.ensemble == "npt"
                else ()
            ),
            options=self.engine.options,
            triclinic=geometry.triclinic,
        )

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        workdir = Path(self.engine.options.get("workdir", ".")) / "lammps_singlepoint"
        payload = self._run(atoms, self.as_singlepoint(simulation), workdir)
        result = self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=self.potential.unit_system_name,
        )
        return self.in_requested_convention(result, simulation)

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
        final = self.in_requested_convention(final, simulation)
        return self._trajectory(payload, simulation, initial, final)

    def _run(self, atoms, simulation: SimulationSpec, workdir: Path) -> dict:
        return lammps_engine.run_job(
            self.potential,
            atoms,
            simulation,
            engine_spec=self.engine,
            workdir=Path(workdir),
            launch=self.launch(),
        )


__all__ = ["LammpsNativeBridge"]
