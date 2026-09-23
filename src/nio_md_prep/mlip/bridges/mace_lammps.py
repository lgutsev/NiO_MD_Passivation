"""MACE on LAMMPS, through ML-IAP or through a MACE pair style.

Two registered implementations, because clusters differ:

``mliap``
    The newer MACE/LAMMPS interface, via ``pair_style mliap unified``. This is
    the route that brings GPU acceleration, multi-GPU inference and atomic
    virials, so it is preferred by default and is the only one that advertises
    a per-atom energy decomposition here.

``pair-mace``
    A LAMMPS built with a MACE pair style. Kept because not every site has an
    ML-IAP build, and declared more conservatively: per-atom virials are not
    assumed.

Both routes run the *exported* model, not the checkpoint the ASE calculator
loads. That export is a real step where numbers can move, which is why MACE
advises care when benchmarking LAMMPS output against the ASE calculator, and
why cross-engine equivalence is an explicit acceptance test in this subsystem
rather than an assumption.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..potentials.mace import MaceAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import LammpsMlipPotentialSpec, SimulationSpec
from ..units import LAMMPS_METAL
from ..engines import lammps_engine
from ..engines.lammps_engine import LammpsEngine
from .base import Bridge

#: Model file suffix conventionally produced by each export route.
EXPORT_SUFFIX = {"mliap": "-mliap.pt", "pair-mace": "-lammps.pt"}


class _MaceLammpsBridge(Bridge):
    """Shared machinery: render a pair style from a MACE spec, then run LAMMPS."""

    potential_kind = "mace"
    engine_kind = "lammps"
    #: Whether this route is trusted to report a per-atom energy decomposition.
    per_atom_energy = False
    #: LAMMPS units this route runs in. MACE is an eV/Angstrom model, so metal.
    units = "metal"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = MaceAdapter(potential)
        self.runtime = LammpsEngine(engine)

    # -- the pair style this route renders --------------------------------

    def exported_model_path(self) -> Path:
        """Where the LAMMPS-ready export of this model is expected to live.

        ``engine.options.model_path`` overrides it, because the export is
        produced outside this package and may be named anything.
        """
        override = self.engine.options.get("model_path")
        if override:
            return Path(override)
        suffix = EXPORT_SUFFIX[self.implementation]
        return self.potential.model_path.with_name(self.potential.model_path.stem + suffix)

    def pair_style(self) -> str:
        raise NotImplementedError

    def pair_coeff(self) -> tuple[str, ...]:
        raise NotImplementedError

    def required_packages(self) -> tuple[str, ...]:
        return tuple(self.registration.packages) if self.registration else ()

    def as_lammps_spec(self) -> LammpsMlipPotentialSpec:
        """Project the MACE spec onto the LAMMPS-native representation.

        This is the point of having one LAMMPS MLIP representation: a MACE
        model routed through LAMMPS becomes an ordinary LAMMPS MLIP spec, and
        the deck renderer, the type map and the manifest handle it exactly as
        they handle DeepMD or PACE.
        """
        elements = self.potential.elements
        return LammpsMlipPotentialSpec(
            label=self.potential.label,
            pair_style=self.pair_style(),
            pair_coeff=self.pair_coeff(),
            type_map={index + 1: symbol for index, symbol in enumerate(elements)},
            required_packages=self.required_packages(),
            units=self.units,
            atom_style="atomic",
            model_paths=(self.exported_model_path(),),
            framework="mace",
            energy_convention=self.potential.energy_convention,
            atomic_reference_energies=self.potential.atomic_reference_energies,
        )

    # -- capabilities -----------------------------------------------------

    def capabilities(self) -> CapabilitySet:
        engine_capabilities = replace(
            self.runtime.capabilities(),
            per_atom_energy=self.per_atom_energy,
            gpu=self.potential.wants_gpu,
        )
        return self.adapter.capabilities().intersect(
            self._route_capabilities(
                engine_capabilities,
                energy_convention=self.potential.energy_convention,
                native_units=LAMMPS_METAL.name,
                notes=(
                    "runs the exported LAMMPS model, not the checkpoint the ASE "
                    "calculator loads; verify equivalence before production use",
                ),
            )
        )

    def availability(self) -> Availability:
        exported = self.exported_model_path()
        if not exported.exists():
            return Availability(
                available=False,
                missing=(),
                detail=(
                    f"the LAMMPS export of this model was not found at {exported}. "
                    "Export it with MACE's own tooling, or point "
                    "engine.options.model_path at it."
                ),
            )
        return self.runtime.availability()

    def atomic_reference_energies(self):
        return self.adapter.atomic_reference_energies()

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        spec = self.as_lammps_spec()
        return {
            "implementation": self.implementation,
            "native_units": LAMMPS_METAL.name,
            "energy_convention": self.potential.energy_convention,
            "lammps_units": self.units,
            "atom_style": spec.atom_style,
            "type_map": {str(k): v for k, v in sorted(spec.type_map.items())},
            "required_packages": list(spec.required_packages),
            "exported_model_path": str(self.exported_model_path()),
            # The exact strings LAMMPS will execute, preserved verbatim.
            "pair_commands": list(spec.render_pair_commands()),
        }

    # -- execution --------------------------------------------------------

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        workdir = Path(self.engine.options.get("workdir", ".")) / "lammps_singlepoint"
        payload = _run_lammps(
            self.as_lammps_spec(),
            atoms,
            self.as_singlepoint(simulation),
            engine_spec=self.engine,
            workdir=workdir,
            want_per_atom=simulation.compute_per_atom_energy,
        )
        return self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=LAMMPS_METAL.name,
        )

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        capabilities = self.validate(simulation, atoms)
        spec = self.as_lammps_spec()
        initial = self.singlepoint(atoms, self.as_singlepoint(simulation))
        payload = _run_lammps(
            spec,
            atoms,
            simulation,
            engine_spec=self.engine,
            workdir=Path(workdir),
            want_per_atom=simulation.compute_per_atom_energy,
        )
        final = self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=LAMMPS_METAL.name,
        )
        return self._trajectory(payload, simulation, initial, final)


class MaceLammpsMliapBridge(_MaceLammpsBridge):
    implementation = "mliap"
    per_atom_energy = True

    def pair_style(self) -> str:
        return f"mliap unified {self.exported_model_path().name} 0"

    def pair_coeff(self) -> tuple[str, ...]:
        return (f"* * {' '.join(self.potential.elements)}",)


class MaceLammpsPairStyleBridge(_MaceLammpsBridge):
    implementation = "pair-mace"
    per_atom_energy = False

    def pair_style(self) -> str:
        return "mace no_domain_decomposition"

    def pair_coeff(self) -> tuple[str, ...]:
        model = self.exported_model_path().name
        return (f"* * {model} {' '.join(self.potential.elements)}",)


def _run_lammps(
    spec: LammpsMlipPotentialSpec,
    atoms,
    simulation: SimulationSpec,
    *,
    engine_spec,
    workdir: Path,
    want_per_atom: bool,
) -> dict:
    """Render, write, run and convert. Shared by both MACE/LAMMPS routes."""
    from ..structures import is_periodic

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    lammps_engine.write_data_file(atoms, spec, workdir / lammps_engine.DATA_FILE)
    commands = lammps_engine.render_deck(
        spec,
        simulation,
        n_types=len(spec.type_map),
        periodic=is_periodic(atoms),
    )
    raw = lammps_engine.run_deck(
        commands,
        workdir=workdir,
        engine_spec=engine_spec,
        n_atoms=len(atoms),
        want_per_atom=want_per_atom,
    )
    payload = _to_canonical(raw, spec, atoms, commands)
    if simulation.task == "md":
        interval = max(1, simulation.trajectory_interval)
        payload["frames"] = simulation.steps // interval + 1
    return payload


def _to_canonical(raw: dict, spec: LammpsMlipPotentialSpec, atoms, commands) -> dict:
    thermo = raw["thermo"]
    energy = lammps_engine.energy_to_eV(thermo["pe"], spec.units)
    forces = lammps_engine.forces_to_eV_per_A(raw["forces"], spec.units)
    per_atom = raw.get("per_atom_energy")
    if per_atom is not None:
        per_atom = [lammps_engine.energy_to_eV(v, spec.units) for v in per_atom]
    final_atoms = atoms.copy()
    if raw.get("positions"):
        final_atoms.set_positions(raw["positions"])
    return {
        "energy_eV": energy,
        "forces_eV_per_A": forces,
        "stress_eV_per_A3": lammps_engine.stress_from_pressure(thermo, spec.units),
        "per_atom_energy_eV": per_atom,
        "symbols": list(atoms.get_chemical_symbols()),
        "wall_time_s": raw.get("wall_time_s"),
        "atoms": final_atoms,
        "trajectory_path": str(Path(raw["deck_path"]).parent / lammps_engine.TRAJECTORY_FILE),
        "log_path": raw.get("log_path"),
        "frames": 0,
        "temperature_start_K": None,
        "temperature_end_K": thermo.get("temp"),
        "max_temperature_K": None,
        "total_energy_drift_eV_per_atom": None,
        "integrator": "lammps fix mlip_integrate",
        "native": {
            "route": raw.get("route"),
            "thermo": thermo,
            "units": spec.unit_system_name,
            "deck_path": raw.get("deck_path"),
            "pair_commands": list(spec.render_pair_commands()),
            "deck": list(commands),
        },
    }


__all__ = ["MaceLammpsMliapBridge", "MaceLammpsPairStyleBridge", "EXPORT_SUFFIX"]
