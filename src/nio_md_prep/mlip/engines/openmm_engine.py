"""The OpenMM engine, via OpenMM-ML.

Three things about this engine are load-bearing and easy to get wrong:

1. **Units.** OpenMM speaks kJ/mol and nm. Every number crossing this
   boundary goes through :mod:`nio_md_prep.mlip.units`; nothing is compared
   in OpenMM's units.
2. **Energy convention.** OpenMM-ML's MACE implementation distinguishes the
   interaction energy from the energy including atomic self-energies, and
   *interaction* is its default. This module asks for whichever convention
   the job declares and records which one it got, so an ASE/OpenMM comparison
   never silently compares the two.
3. **No stress.** OpenMM reports no virial here. The capability set says so,
   and a constant-pressure job is refused during negotiation rather than run
   with an absent stress tensor.

Only full-system MLIP is implemented. OpenMM-ML can build mixed ML/MM systems,
but that machinery -- and the still-active questions around bonds crossing an
ML/MM boundary -- is deliberately out of scope here.
"""
from __future__ import annotations

import time
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, probe
from ..errors import ConfigError, MissingDependencyError
from ..specs import SimulationSpec
from ..units import INTERACTION, OPENMM, TOTAL, to_canonical
from .base import EngineRuntime

#: This subsystem's convention name -> OpenMM-ML's ``returnEnergyType`` value.
ENERGY_TYPE = {
    INTERACTION: "interaction_energy",
    TOTAL: "energy",
}

GPU_PLATFORMS = ("CUDA", "OpenCL", "HIP")
DEFAULT_FRICTION_PER_PS = 1.0


class OpenMMEngine(EngineRuntime):
    kind = "openmm"
    requires = ("openmm", "openmmml")
    native_units = OPENMM.name

    def capabilities(self) -> CapabilitySet:
        """What OpenMM can carry -- notably not a stress tensor.

        ``stress=False`` and ``per_atom_energy=False`` are the honest answers
        for this route, and they are what make an NPT request fail during
        negotiation instead of producing a trajectory with no virial behind it.
        """
        platform = self.spec.platform
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=False,
            per_atom_energy=False,
            periodic=True,
            gpu=platform is None or platform in GPU_PLATFORMS,
            elements=None,
            precisions=None,
            engines=frozenset({"openmm"}),
            native_units=self.native_units,
            notes=(
                "OpenMM reports no stress tensor through OpenMM-ML, so constant-"
                "pressure jobs are rejected rather than run without a virial",
                "no per-atom energy decomposition is exposed by this route",
            ),
        )

    def availability(self) -> Availability:
        return probe(self.requires, detail="OpenMM-ML machine-learned potentials")


def require_openmm():
    try:
        import openmm  # noqa: F401
        import openmmml  # noqa: F401
    except ImportError as exc:
        raise MissingDependencyError(
            "the OpenMM engine",
            ("openmm", "openmm-ml"),
            hint="Install the OpenMM extra: pip install 'nio-md-prep[openmm]'",
        ) from exc


def build_topology(atoms):
    """Build an OpenMM ``Topology`` from an ASE ``Atoms``.

    One chain, one residue, one atom per site. That is sufficient for a
    full-system MLIP -- the model reads elements and coordinates, not
    residues -- and it avoids inventing a topology this layer does not have.
    """
    require_openmm()
    from openmm import Vec3
    from openmm.app import Element, Topology
    from openmm.unit import nanometer

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MLIP", chain)
    for symbol in atoms.get_chemical_symbols():
        topology.addAtom(symbol, Element.getBySymbol(symbol), residue)
    if all(atoms.pbc):
        cell = atoms.get_cell()
        vectors = [Vec3(*(float(c) / 10.0 for c in row)) for row in cell]
        topology.setPeriodicBoxVectors(vectors * nanometer)
    return topology


def positions_nm(atoms):
    from openmm import Vec3
    from openmm.unit import nanometer

    return [Vec3(*(float(c) / 10.0 for c in row)) for row in atoms.get_positions()] * nanometer


def build_system(atoms, model_path, *, energy_convention: str, potential_name: str = "mace"):
    """Create the OpenMM ``System`` for a model, under an explicit convention.

    ``energy_convention`` is passed through to OpenMM-ML's ``returnEnergyType``
    rather than left at its default, so the convention in the manifest is the
    convention that ran.
    """
    require_openmm()
    from openmmml import MLPotential

    if energy_convention not in ENERGY_TYPE:
        raise ConfigError(
            f"energy convention {energy_convention!r} has no OpenMM-ML equivalent; "
            f"expected one of {', '.join(ENERGY_TYPE)}"
        )
    topology = build_topology(atoms)
    potential = MLPotential(potential_name, modelPath=str(model_path))
    system = potential.createSystem(
        topology, returnEnergyType=ENERGY_TYPE[energy_convention]
    )
    return system, topology


def make_context(system, atoms, engine_spec, *, timestep_fs: float = 1.0):
    require_openmm()
    import openmm
    from openmm import unit

    integrator = openmm.VerletIntegrator(timestep_fs * 0.001 * unit.picosecond)
    platform = None
    properties = {}
    if engine_spec.platform:
        platform = openmm.Platform.getPlatformByName(engine_spec.platform)
        if engine_spec.precision and engine_spec.platform in GPU_PLATFORMS:
            properties["Precision"] = {"float32": "single", "float64": "double"}[
                engine_spec.precision
            ]
    context = (
        openmm.Context(system, integrator, platform, properties)
        if platform is not None
        else openmm.Context(system, integrator)
    )
    context.setPositions(positions_nm(atoms))
    if all(atoms.pbc):
        from openmm import Vec3

        cell = atoms.get_cell()
        context.setPeriodicBoxVectors(
            *[Vec3(*(float(c) / 10.0 for c in row)) for row in cell]
        )
    return context, integrator


def singlepoint(atoms, system, engine_spec) -> dict:
    """Evaluate one geometry and convert straight out of OpenMM's units."""
    require_openmm()
    started = time.perf_counter()
    context, _ = make_context(system, atoms, engine_spec)
    state = context.getState(getEnergy=True, getForces=True)
    from openmm import unit

    energy_kj_per_mol = state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    forces_kj_per_mol_nm = [
        [float(c) for c in row.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)]
        for row in state.getForces()
    ]
    return {
        "energy_eV": to_canonical(energy_kj_per_mol, "energy", OPENMM),
        "forces_eV_per_A": to_canonical(forces_kj_per_mol_nm, "force", OPENMM),
        "stress_eV_per_A3": None,
        "per_atom_energy_eV": None,
        "symbols": list(atoms.get_chemical_symbols()),
        "wall_time_s": time.perf_counter() - started,
        "native": {
            "energy_kJ_per_mol": energy_kj_per_mol,
            "units": OPENMM.name,
        },
    }


def run_md(atoms, system, engine_spec, simulation: SimulationSpec, *, workdir: Path) -> dict:
    """Run a short diagnostic trajectory through OpenMM."""
    require_openmm()
    import openmm
    from openmm import app, unit

    if simulation.ensemble == "npt":
        raise ConfigError(
            "constant-pressure dynamics are not offered by this OpenMM route: no "
            "stress tensor is available to integrate the cell against"
        )
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    timestep = simulation.timestep_fs * 0.001 * unit.picosecond
    temperature = (simulation.temperature_K or 0.0) * unit.kelvin
    if simulation.ensemble == "nvt":
        friction = (
            1000.0 / simulation.thermostat_damping_fs
            if simulation.thermostat_damping_fs
            else DEFAULT_FRICTION_PER_PS
        ) / unit.picosecond
        integrator = openmm.LangevinMiddleIntegrator(temperature, friction, timestep)
    else:
        integrator = openmm.VerletIntegrator(timestep)
    if simulation.seed is not None and hasattr(integrator, "setRandomNumberSeed"):
        integrator.setRandomNumberSeed(int(simulation.seed))

    topology = build_topology(atoms)
    platform = (
        openmm.Platform.getPlatformByName(engine_spec.platform)
        if engine_spec.platform
        else None
    )
    simulation_object = (
        app.Simulation(topology, system, integrator, platform)
        if platform is not None
        else app.Simulation(topology, system, integrator)
    )
    simulation_object.context.setPositions(positions_nm(atoms))
    if simulation.temperature_K:
        simulation_object.context.setVelocitiesToTemperature(
            temperature, int(simulation.seed) if simulation.seed is not None else 0
        )

    trajectory_path = workdir / "smoke_md.dcd"
    log_path = workdir / "smoke_md.log"
    simulation_object.reporters.append(
        app.DCDReporter(str(trajectory_path), max(1, simulation.trajectory_interval))
    )
    simulation_object.reporters.append(
        app.StateDataReporter(
            str(log_path),
            max(1, simulation.log_interval),
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
        )
    )

    started = time.perf_counter()
    initial = simulation_object.context.getState(getEnergy=True)
    initial_total = (
        initial.getPotentialEnergy() + initial.getKineticEnergy()
    ).value_in_unit(unit.kilojoule_per_mole)
    simulation_object.step(simulation.steps)
    final_state = simulation_object.context.getState(getEnergy=True, getPositions=True)
    final_total = (
        final_state.getPotentialEnergy() + final_state.getKineticEnergy()
    ).value_in_unit(unit.kilojoule_per_mole)
    wall_time = time.perf_counter() - started

    final_atoms = atoms.copy()
    final_atoms.set_positions(
        [
            [float(c) * 10.0 for c in row.value_in_unit(unit.nanometer)]
            for row in final_state.getPositions()
        ]
    )
    drift_kj = (final_total - initial_total) / len(atoms)
    return {
        "atoms": final_atoms,
        "frames": max(1, simulation.steps // max(1, simulation.trajectory_interval)),
        "trajectory_path": str(trajectory_path),
        "log_path": str(log_path),
        "total_energy_drift_eV_per_atom": to_canonical(drift_kj, "energy", OPENMM),
        "wall_time_s": wall_time,
        "integrator": type(integrator).__name__,
    }


__all__ = [
    "ENERGY_TYPE",
    "GPU_PLATFORMS",
    "OpenMMEngine",
    "require_openmm",
    "build_topology",
    "build_system",
    "make_context",
    "singlepoint",
    "run_md",
]
