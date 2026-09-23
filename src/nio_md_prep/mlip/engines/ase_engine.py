"""The ASE engine: single points, optimisation and short trajectories.

ASE's native units are already the canonical ones (eV, eV/Angstrom,
eV/Angstrom^3), so no conversion happens here -- but the unit system is still
recorded on every result, because "no conversion was needed" and "no
conversion was done" must be distinguishable in a manifest.

Everything ASE-specific is imported inside the functions, so this module can
be imported to ask a capability question on a machine without ASE.
"""
from __future__ import annotations

import time
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, probe
from ..errors import ConfigError, MissingDependencyError
from ..results import PotentialResult
from ..specs import DEFAULT_BAROSTAT_DAMPING_FS, DEFAULT_THERMOSTAT_DAMPING_FS, SimulationSpec
from ..units import ASE
from .base import EngineRuntime

#: Thermostat -> the ASE integrator that implements it. ``nose-hoover`` maps to
#: ASE's Nose-Hoover chain NVT, which is only present in newer ASE releases;
#: the import failure is turned into a clear message rather than a traceback.
NVT_INTEGRATORS = ("langevin", "nose-hoover", "berendsen")
DEFAULT_THERMOSTAT = "langevin"


class AseEngine(EngineRuntime):
    kind = "ase"
    requires = ("ase", "numpy")
    native_units = ASE.name

    def capabilities(self) -> CapabilitySet:
        """ASE can carry anything a calculator chooses to report.

        So every flag is True here, and the route's real limits come from
        intersecting with the potential's capability set.
        """
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=True,
            per_atom_energy=True,
            periodic=True,
            # ASE keeps atoms.pbc per axis and applies FixAtoms to momenta in
            # its integrators and to the degrees of freedom in get_temperature.
            partial_periodic=True,
            fixed_atoms=True,
            gpu=True,
            elements=None,
            precisions=None,
            engines=frozenset({"ase"}),
            native_units=self.native_units,
            notes=("ASE reports eV and eV/Angstrom natively; no unit conversion applies",),
        )

    def availability(self) -> Availability:
        return probe(self.requires, detail="ASE optimisers and integrators")


def require_ase():
    try:
        import ase  # noqa: F401
    except ImportError as exc:
        raise MissingDependencyError(
            "the ASE engine",
            ("ase",),
            hint="Install the MLIP extra: pip install 'nio-md-prep[mlip]'",
        ) from exc


def singlepoint(
    atoms,
    calculator,
    *,
    compute_stress: bool = False,
    compute_per_atom_energy: bool = False,
) -> dict:
    """Evaluate one geometry and return canonical quantities as plain data.

    Returns a mapping rather than a :class:`PotentialResult` so the calling
    bridge can attach the provenance fields it alone knows (which
    implementation, which convention, which native units).
    """
    require_ase()
    started = time.perf_counter()
    atoms = atoms.copy()
    atoms.calc = calculator
    energy = float(atoms.get_potential_energy())
    forces = [[float(c) for c in row] for row in atoms.get_forces()]
    stress = None
    if compute_stress:
        stress = [float(v) for v in atoms.get_stress(voigt=True)]
    per_atom = None
    if compute_per_atom_energy:
        per_atom = [float(v) for v in atoms.get_potential_energies()]
    return {
        "energy_eV": energy,
        "forces_eV_per_A": forces,
        "stress_eV_per_A3": stress,
        "per_atom_energy_eV": per_atom,
        "symbols": list(atoms.get_chemical_symbols()),
        "wall_time_s": time.perf_counter() - started,
    }


def optimize(atoms, calculator, simulation: SimulationSpec, *, workdir: Path | None = None):
    """Relax a geometry with BFGS. Returns the relaxed ``Atoms`` and a summary."""
    require_ase()
    from ase.optimize import BFGS

    atoms = atoms.copy()
    atoms.calc = calculator
    logfile = str(Path(workdir) / "optimize.log") if workdir else None
    started = time.perf_counter()
    optimizer = BFGS(atoms, logfile=logfile)
    converged = optimizer.run(fmax=simulation.fmax_eV_per_A, steps=simulation.max_optimizer_steps)
    return atoms, {
        "converged": bool(converged),
        "steps": int(optimizer.get_number_of_steps()),
        "fmax_eV_per_A": simulation.fmax_eV_per_A,
        "wall_time_s": time.perf_counter() - started,
    }


def run_md(
    atoms,
    calculator,
    simulation: SimulationSpec,
    *,
    workdir: Path,
    trajectory_name: str = "smoke_md.traj",
) -> dict:
    """Run a short trajectory and report what a diagnostic run should report.

    The point of this function is not production dynamics; it is to prove that
    a route integrates stably. So it returns the energy drift per atom, the
    temperature at both ends, and the peak temperature -- the three numbers
    that reveal a broken timestep, a broken neighbour list or a broken
    unit conversion within a few hundred steps.
    """
    require_ase()
    import numpy as np
    from ase import units as ase_units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    atoms = atoms.copy()
    atoms.calc = calculator

    temperature = simulation.temperature_K or 0.0
    if temperature > 0:
        rng = np.random.default_rng(simulation.seed)
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature, rng=rng)
        Stationary(atoms)

    timestep = simulation.timestep_fs * ase_units.fs
    dynamics = _make_dynamics(atoms, simulation, timestep)

    trajectory_path = workdir / trajectory_name
    log_path = workdir / "smoke_md.log"
    from ase.io.trajectory import Trajectory

    trajectory = Trajectory(str(trajectory_path), "w", atoms)
    dynamics.attach(trajectory.write, interval=max(1, simulation.trajectory_interval))

    # (MD step, potential, kinetic, temperature), one row per step: ASE calls
    # observers at step 0 and the explicit calls below would repeat a step.
    history: list[tuple[int, float, float, float]] = []

    def record():
        step = int(dynamics.nsteps)
        if history and history[-1][0] == step:
            return
        potential = atoms.get_potential_energy()
        kinetic = atoms.get_kinetic_energy()
        history.append((step, potential, kinetic, atoms.get_temperature()))

    dynamics.attach(record, interval=max(1, simulation.log_interval))

    started = time.perf_counter()
    record()
    dynamics.run(simulation.steps)
    record()
    trajectory.close()
    wall_time = time.perf_counter() - started

    from ..diagnostics import temperature_ndof
    from ..structures import constraint_summary

    temperatures = [t for _, _, _, t in history]
    log_path.write_text(
        "# step time_fs potential_eV kinetic_eV total_eV temperature_K\n"
        + "".join(
            f"{step} {step * simulation.timestep_fs:.6f} {potential:.10f} {kinetic:.10f} "
            f"{potential + kinetic:.10f} {temp:.6f}\n"
            for step, potential, kinetic, temp in history
        ),
        encoding="utf-8",
    )
    steps_completed = int(dynamics.nsteps)
    interval = max(1, simulation.trajectory_interval)
    frames_written = _count_frames(trajectory_path)
    expected_frames = steps_completed // interval + 1
    if frames_written != expected_frames:
        from ..errors import ResultError

        raise ResultError(
            f"{trajectory_path} holds {frames_written} frames; {steps_completed} steps "
            f"written every {interval} steps (including step 0) should give "
            f"{expected_frames}"
        )
    constraints = constraint_summary(atoms)
    return {
        "atoms": atoms,
        "steps_completed": steps_completed,
        "frames_written": frames_written,
        "trajectory_path": str(trajectory_path),
        "log_path": str(log_path),
        "final_positions": atoms.get_positions().tolist(),
        "final_cell": atoms.get_cell().tolist(),
        "final_pbc": [bool(v) for v in atoms.pbc],
        "energy_series": {
            "time_fs": [step * simulation.timestep_fs for step, _, _, _ in history],
            "total_energy_eV": [potential + kinetic for _, potential, kinetic, _ in history],
        },
        "temperature_start_K": temperatures[0] if temperatures else None,
        "temperature_end_K": temperatures[-1] if temperatures else None,
        "max_temperature_K": max(temperatures) if temperatures else None,
        # ASE's get_temperature: 3N minus the DOF FixAtoms removes; the
        # Stationary COM removal is not subtracted.
        "temperature_ndof": temperature_ndof(len(atoms), n_fixed=constraints["n_fixed"]),
        "constraints": (
            {
                "fixed_atoms": constraints["fixed_atoms"],
                "n_fixed": constraints["n_fixed"],
                "method": "ASE FixAtoms: momenta of frozen atoms zeroed by the integrator",
            }
            if constraints["n_fixed"]
            else None
        ),
        "wall_time_s": wall_time,
        "integrator": type(dynamics).__name__,
    }


def _make_dynamics(atoms, simulation: SimulationSpec, timestep: float):
    """Pick the ASE integrator for the requested ensemble."""
    from ase import units as ase_units

    ensemble = simulation.ensemble
    if ensemble == "nve":
        from ase.md.verlet import VelocityVerlet

        return VelocityVerlet(atoms, timestep)

    damping_fs = simulation.thermostat_damping_fs or DEFAULT_THERMOSTAT_DAMPING_FS
    thermostat = simulation.thermostat or DEFAULT_THERMOSTAT

    if ensemble == "nvt":
        if thermostat == "langevin":
            from ase.md.langevin import Langevin

            # ASE's Langevin friction is a rate; 1/tau with tau the damping time.
            return Langevin(
                atoms,
                timestep,
                temperature_K=simulation.temperature_K,
                friction=1.0 / (damping_fs * ase_units.fs),
                rng=_rng(simulation.seed),
            )
        if thermostat == "berendsen":
            from ase.md.nvtberendsen import NVTBerendsen

            return NVTBerendsen(
                atoms,
                timestep,
                temperature_K=simulation.temperature_K,
                taut=damping_fs * ase_units.fs,
            )
        if thermostat == "nose-hoover":
            try:
                from ase.md.nose_hoover_chain import NoseHooverChainNVT
            except ImportError as exc:
                raise ConfigError(
                    "this ASE installation has no Nose-Hoover chain NVT integrator; "
                    "use thermostat = 'langevin' or 'berendsen', or upgrade ASE"
                ) from exc
            return NoseHooverChainNVT(
                atoms,
                timestep,
                temperature_K=simulation.temperature_K,
                tdamp=damping_fs * ase_units.fs,
            )
        raise ConfigError(
            f"thermostat {thermostat!r} is not implemented for the ASE engine; "
            f"available: {', '.join(NVT_INTEGRATORS)}"
        )

    if ensemble == "npt":
        return _make_npt(atoms, simulation, timestep, damping_fs)

    raise ConfigError(f"ensemble {ensemble!r} is not implemented for the ASE engine")


#: The bulk-modulus guess ASE's NPT needs to turn a barostat time into its
#: ``pfactor``. Only a diagnostic default -- NPT here is for proving that a
#: route integrates, not for production equations of state.
DEFAULT_BULK_MODULUS_GPA = 100.0


def _make_npt(atoms, simulation: SimulationSpec, timestep: float, damping_fs: float):
    import numpy as np
    from ase import units as ase_units
    from ase.md.npt import NPT

    from ..units import EV_PER_ANGSTROM3_IN_BAR

    cell = np.array(atoms.get_cell())
    if np.abs(np.tril(cell, -1)).max() > 1e-8:
        raise ConfigError(
            "ASE's NPT integrator requires an upper-triangular cell, and this "
            "structure's cell is not. Rotate the cell first, or run the constant-"
            "pressure diagnostic through the LAMMPS engine instead."
        )
    pressure_eV_per_A3 = simulation.pressure_bar / EV_PER_ANGSTROM3_IN_BAR
    barostat_fs = simulation.barostat_damping_fs or DEFAULT_BAROSTAT_DAMPING_FS
    bulk_modulus = DEFAULT_BULK_MODULUS_GPA
    return NPT(
        atoms,
        timestep,
        temperature_K=simulation.temperature_K,
        externalstress=pressure_eV_per_A3,
        ttime=damping_fs * ase_units.fs,
        pfactor=(barostat_fs * ase_units.fs) ** 2 * bulk_modulus * ase_units.GPa,
    )


def _count_frames(path: Path) -> int:
    """Frames in an ASE trajectory, read back from disk. Unreadable is an error."""
    from ase.io.trajectory import Trajectory

    from ..errors import ResultError

    if not Path(path).exists():
        raise ResultError(f"the trajectory {path} was not written")
    try:
        with Trajectory(str(path), "r") as handle:
            return len(handle)
    except Exception as exc:
        raise ResultError(f"the trajectory {path} cannot be read: {exc}") from exc


def _rng(seed: int | None):
    import numpy as np

    return np.random.default_rng(seed)


def build_result(payload: dict, **provenance) -> PotentialResult:
    """Wrap :func:`singlepoint` output with the bridge's provenance fields."""
    return PotentialResult(
        energy_eV=payload["energy_eV"],
        forces_eV_per_A=payload["forces_eV_per_A"],
        symbols=payload["symbols"],
        stress_eV_per_A3=payload["stress_eV_per_A3"],
        per_atom_energy_eV=payload["per_atom_energy_eV"],
        wall_time_s=payload["wall_time_s"],
        **provenance,
    )


__all__ = [
    "AseEngine",
    "NVT_INTEGRATORS",
    "DEFAULT_THERMOSTAT",
    "require_ase",
    "singlepoint",
    "optimize",
    "run_md",
    "build_result",
]
