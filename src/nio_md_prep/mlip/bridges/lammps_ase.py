"""A LAMMPS-native MLIP driven from ASE.

LAMMPS still evaluates the potential -- this is not a reimplementation. ASE's
``LAMMPSlib`` calculator starts a LAMMPS instance in-process, feeds it the
pair commands, and drives it, which is what makes ASE's optimisers,
integrators and analysis available to a pair style that only exists inside
LAMMPS.

The ``type_map`` is what makes the bridge possible: it maps an ASE ``Atoms``
object's elements onto the LAMMPS types the ``pair_coeff`` line expects.
Without it the two would agree only by accident.

What ``LAMMPSlib`` (ASE 3.29, ``ase/calculators/lammpslib.py``) actually does,
and therefore what this route claims:

* energy and forces, converted from the LAMMPS unit style by ASE;
* stress from thermo ``pxx..pxy`` through ``Prism.tensor2_to_ase``. Its
  single points never hand LAMMPS any velocities (``calculate`` runs
  ``propagate`` with ``n_steps = 0``), so that is the virial-only stress;
* per-atom energies (``energies``, from ``compute pe/atom``) -- claimed only
  when the pair style is known to tally them, and checked to sum to the
  energy;
* per-axis boundaries: ``p`` for a periodic axis, ``f`` for a non-periodic
  one with a finite cell vector. Atoms outside such a fixed face would be
  lost, so the same preflight as the native route refuses them;
* one MPI rank only (it raises for ``world_size != 1``).

Model files named by the pair commands are staged into the job directory,
hashed, and passed to LAMMPS by absolute path, because the in-process LAMMPS
resolves relative names against this process's working directory. The
LAMMPS log goes to ``lammpslib.log`` in the job directory.

Note what this bridge does *not* imply. Driving a LAMMPS pair style from ASE
works because ASE can call into LAMMPS. There is no equivalent route into
OpenMM, and none is faked: see the explicitly-unsupported ``lammps -> openmm``
cell in :mod:`nio_md_prep.mlip.registry`.
"""
from __future__ import annotations

import inspect
from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, module_available
from ..errors import ConfigError, ResultError
from ..potentials.lammps_mlip import LammpsMlipAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import lammps_unit_system
from ..engines import ase_engine, lammps_engine
from ..engines.ase_engine import AseEngine
from .base import Bridge

LOG_FILE = "lammpslib.log"
#: Job-directory name for single points when no MD workdir is given.
SINGLEPOINT_DIR = "lammps_ase_singlepoint"


class LammpsAseBridge(Bridge):
    potential_kind = "lammps"
    engine_kind = "ase"
    implementation = "ase-lammps"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = LammpsMlipAdapter(potential)
        self.runtime = AseEngine(engine)
        self._calculators: dict[str, tuple[object, list[dict]]] = {}

    def capabilities(self) -> CapabilitySet:
        """What ASE's LAMMPSlib calculator delivers for this pair style.

        Energy, forces, virial-only stress, per-axis periodicity and (through
        the ASE engine) FixAtoms. Per-atom energies only when the pair style's
        framework is known to tally them; no GPU (no accelerator switches are
        passed to the in-process LAMMPS).
        """
        engine_capabilities = replace(
            self.runtime.capabilities(),
            gpu=False,
            notes=self.runtime.capabilities().notes
            + (
                "ASE's LAMMPSlib evaluates the pair style in-process on one MPI rank; its "
                "stress is the virial-only thermo pressure rotated to the source basis",
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
                "ASE's LAMMPSlib calculator needs the LAMMPS python module (an lmp "
                "executable is not enough)",
            )
        return Availability(True, (), "ASE driving an in-process LAMMPS (python module)")

    # -- validation -------------------------------------------------------

    def check_simulation(self, simulation: SimulationSpec, atoms=None) -> None:
        """Refuse what LAMMPSlib or the ASE engine would not do as asked.

        ``engine.threads`` and ``engine.precision`` cannot be applied to the
        in-process LAMMPS; the structure gets the native route's geometry
        preflight (LAMMPSlib uses the same fixed faces); the MD request is
        checked against the ASE engine's integrator map.
        """
        if self.engine.threads is not None:
            raise ConfigError(
                "engine.threads cannot be applied on the lammps -> ase route: LAMMPSlib "
                "starts LAMMPS without accelerator switches and this route passes none. "
                "Remove engine.threads, or use engine = 'lammps' (OpenMP threads)."
            )
        if self.engine.precision is not None:
            raise ConfigError(
                f"engine.precision = {self.engine.precision!r} cannot be applied or verified "
                "for a LAMMPS-native pair style: its precision is fixed by the LAMMPS package "
                "and model. Remove engine.precision."
            )
        lammps_engine.check_units(self.potential.units, simulation)
        lammps_engine.staged_model_names(self.potential)
        if atoms is not None:
            lammps_engine.prepare_geometry(atoms)
        check = getattr(ase_engine, "check_simulation", None)
        if check is not None:
            check(simulation, atoms, options=self.engine.options)
        elif simulation.task == "md" and (
            simulation.ensemble == "npt" or simulation.thermostat == "csvr"
        ):
            # An ASE engine without an explicit integrator map runs every NPT
            # request with ase.md.npt.NPT and has no CSVR thermostat.
            raise ConfigError(
                "this ASE engine has no explicit integrator map for "
                f"{'npt' if simulation.ensemble == 'npt' else 'thermostat csvr'}; refused "
                "rather than run with a different integrator"
            )

    # -- the calculator ---------------------------------------------------

    def _workdir(self, workdir: Path | None) -> Path:
        if workdir is not None:
            return Path(workdir)
        return Path(self.engine.options.get("workdir", ".")) / SINGLEPOINT_DIR

    def log_path(self, workdir: Path | None = None) -> Path:
        """Where LAMMPSlib's log goes: ``engine.options.log_file`` or the job directory."""
        override = self.engine.options.get("log_file")
        if override:
            return Path(str(override))
        return self._workdir(workdir) / LOG_FILE

    def calculator(self, workdir: Path | None = None):
        """ASE's LAMMPSlib around this pair style, with staged, hashed model files."""
        workdir = self._workdir(workdir)
        key = str(workdir.resolve())
        if key in self._calculators:
            return self._calculators[key][0]
        self.require_available()
        from ase.calculators.lammpslib import LAMMPSlib

        staged, records = lammps_engine.stage_model_files(self.potential, workdir, absolute=True)
        log_path = self.log_path(workdir)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        calculator = LAMMPSlib(
            lmpcmds=list(staged.render_pair_commands()),
            atom_types={
                symbol: type_id for type_id, symbol in sorted(self.potential.type_map.items())
            },
            lammps_header=self.lammps_header(),
            keep_alive=True,
            log_file=str(log_path),
        )
        self._calculators[key] = (calculator, records)
        return calculator

    def lammps_header(self) -> list[str]:
        header = [
            f"units {self.potential.units}",
            f"atom_style {self.potential.atom_style}",
            "atom_modify map array sort 0 0",
        ]
        if self.potential.newton:
            header.insert(0, f"newton {self.potential.newton}")
        return header

    def _staged_records(self, workdir: Path | None) -> list[dict]:
        entry = self._calculators.get(str(self._workdir(workdir).resolve()))
        return list(entry[1]) if entry else []

    # -- description ------------------------------------------------------

    def _integrator(self, simulation: SimulationSpec | None) -> dict | None:
        if simulation is None or simulation.task != "md":
            return None
        describe = getattr(ase_engine, "integrator_parameters", None)
        if describe is None:
            return {"resolved_by": "the ASE engine at run time (no integrator map available)"}
        return describe(simulation, self.engine.options)

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        parameters = {
            "implementation": self.implementation,
            "calculator": "ase.calculators.lammpslib.LAMMPSlib",
            "native_units": "ase",
            "energy_convention": self.potential.energy_convention,
            "lammps_units": self.potential.units,
            "lammps_unit_system": lammps_unit_system(self.potential.units).name,
            "lammps_header": self.lammps_header(),
            "atom_types": {v: k for k, v in sorted(self.potential.type_map.items())},
            "framework": self.adapter.framework,
            "required_packages": list(self.adapter.required_packages()),
            # The strings handed to the in-process LAMMPS: model files are
            # replaced by the absolute path of their staged copy at run time.
            "pair_commands": list(self.potential.render_pair_commands()),
            "staged_model_files": [
                {"source": str(p), "staged_name": name}
                for p, name in zip(
                    self.potential.model_paths,
                    lammps_engine.staged_model_names(self.potential).values(),
                )
            ],
            "log_file": str(self.log_path()),
            "stress": "LAMMPSlib thermo pxx..pxy with zero velocities (virial only)",
        }
        integrator = self._integrator(simulation)
        if integrator is not None:
            parameters["integrator"] = integrator
        return parameters

    def execution_plan(self, simulation: SimulationSpec, atoms=None) -> dict:
        """What a run would execute, resolved without running (``mlip validate``)."""
        from ase import units as ase_units

        files = lammps_engine.model_file_plan(
            self.potential.model_paths, self.potential.model_hashes
        )
        first = self.potential.model_paths[0] if self.potential.model_paths else None
        integrator = self._integrator(simulation)
        dynamics = None
        if integrator is not None:
            timestep = simulation.timestep_fs
            t_damp = (
                simulation.resolved_thermostat_damping_fs if simulation.ensemble != "nve" else None
            )
            p_damp = (
                simulation.resolved_barostat_damping_fs if simulation.ensemble == "npt" else None
            )
            dynamics = {
                "ensemble": simulation.ensemble,
                "integrator": integrator.get("class") or integrator.get("integrator"),
                "thermostat": integrator.get("thermostat"),
                "barostat": integrator.get("barostat"),
                "barostat_coupling": integrator.get("barostat_coupling"),
                "timestep_fs": timestep,
                "timestep_native": timestep * ase_units.fs,
                "native_time_unit": (
                    f"ASE time unit (Angstrom*sqrt(amu/eV); 1 fs = {ase_units.fs!r})"
                ),
                "thermostat_damping_fs": t_damp,
                "thermostat_damping_native": None if t_damp is None else t_damp * ase_units.fs,
                "barostat_damping_fs": p_damp,
                "barostat_damping_native": None if p_damp is None else p_damp * ase_units.fs,
                "resolved": integrator,
            }
            if "refused" in integrator:
                dynamics["refused"] = integrator["refused"]
        pair_coeff = [
            c if c.startswith("pair_coeff") else f"pair_coeff {c}"
            for c in self.potential.pair_coeff
        ]
        boundary = "per axis from structure.pbc: p (periodic) or f (fixed face)"
        if atoms is not None:
            try:
                boundary = lammps_engine.boundary_string(
                    lammps_engine.prepare_geometry(atoms).pbc
                )
            except ConfigError as exc:
                boundary = f"refused: {exc}"
        return {
            "potential_kind": self.potential_kind,
            "engine": self.engine_kind,
            "bridge": self.label,
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
            "energy_convention": self.potential.energy_convention,
            "units": {
                **lammps_engine.units_plan(self.potential.units),
                "native": "ase",
                "conversion": (
                    "ASE LAMMPSlib converts LAMMPS units to eV/Angstrom "
                    "(ase.calculators.lammps.convert)"
                ),
            },
            "device": {
                "requested": None,
                "effective": "cpu",
                "guaranteed": True,
                "note": "LAMMPSlib starts LAMMPS without accelerator switches",
            },
            "precision": {
                "requested": self.engine.precision,
                "effective": None,
                "guaranteed": False,
                "note": "set by the pair style's package and model build",
            },
            "dynamics": dynamics,
            "lammps": {
                "pair_style": f"pair_style {self.potential.pair_style}",
                "pair_coeff": pair_coeff,
                "pair_commands": list(self.potential.render_pair_commands()),
                "model_file_tokens": "replaced by the absolute path of the staged copy",
                "lammps_header": self.lammps_header(),
                "boundary": boundary,
                "launch_command": (
                    "python: ase.calculators.lammpslib.LAMMPSlib -> lammps.lammps(cmdargs="
                    f"['-echo', 'log', '-log', {str(self.log_path())!r}, '-screen', 'none', "
                    "'-nocite'])"
                ),
                "launch_argv": None,
                "kokkos_args": [],
                "env": {},
                "route": "python (ASE LAMMPSlib, in-process, one MPI rank)",
            },
            "openmm": None,
            "model_hashes": files,
            **lammps_engine.bridge_plan_basics(self, simulation, atoms),
        }

    # -- execution --------------------------------------------------------

    def _evaluate(self, atoms, simulation: SimulationSpec, workdir: Path | None) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        calculator = self.calculator(workdir)
        want_per_atom = bool(simulation.compute_per_atom_energy)
        payload = ase_engine.singlepoint(
            atoms,
            calculator,
            compute_stress=simulation.compute_stress or simulation.ensemble == "npt",
            compute_per_atom_energy=want_per_atom,
        )
        if want_per_atom:
            total = float(sum(payload["per_atom_energy_eV"]))
            energy = float(payload["energy_eV"])
            if abs(total - energy) > 1e-6 * max(1.0, abs(energy)):
                raise ResultError(
                    f"LAMMPSlib per-atom energies sum to {total!r} eV but its energy is "
                    f"{energy!r} eV; the pair style does not tally a complete per-atom energy"
                )
        payload["native"] = {
            "calculator": "ase.calculators.lammpslib.LAMMPSlib",
            "log_path": str(self.log_path(workdir)),
            "staged_model_files": self._staged_records(workdir),
            "lammps_units": self.potential.units,
        }
        return self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units="ase",
        )

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        return self._evaluate(atoms, simulation, None)

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        self.validate(simulation, atoms)
        workdir = Path(workdir)
        endpoint = self.as_singlepoint(simulation)
        initial = self._evaluate(atoms, endpoint, workdir)
        kwargs = {}
        if "options" in inspect.signature(ase_engine.run_md).parameters:
            kwargs["options"] = self.engine.options
        payload = ase_engine.run_md(
            atoms, self.calculator(workdir), simulation, workdir=workdir, **kwargs
        )
        payload.setdefault("log_path", str(self.log_path(workdir)))
        final = self._evaluate(payload["atoms"], endpoint, workdir)
        return self._trajectory(payload, simulation, initial, final)


__all__ = ["LammpsAseBridge"]
