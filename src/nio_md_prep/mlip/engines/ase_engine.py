"""The ASE engine: single points, optimisation and short trajectories.

ASE's native units are already the canonical ones (eV, eV/Angstrom,
eV/Angstrom^3), so no conversion happens here -- but the unit system is still
recorded on every result, because "no conversion was needed" and "no
conversion was done" must be distinguishable in a manifest.

Everything ASE-specific is imported inside the functions, so this module can
be imported to ask a capability question on a machine without ASE.

**What integrates.** :func:`resolve_integrator` maps a request onto one ASE
integrator, and :func:`check_simulation` refuses -- at ``validate`` time --
anything this engine would otherwise have to substitute:

==========  ==============================  ====================================
ensemble    request                         ASE integrator
==========  ==============================  ====================================
nve         (a thermostat is refused)       ``VelocityVerlet``
nvt         ``langevin`` (default)          ``Langevin(fixcm=False)``
nvt         ``berendsen``                   ``NVTBerendsen(fixcm=False)``
nvt         ``nose-hoover``                 ``NoseHooverChainNVT`` (refused with
                                            frozen atoms)
nvt         ``csvr``                        ``Bussi``
npt         ``mtk`` (default), isotropic    ``IsotropicMTKNPT``
npt         ``mtk``, anisotropic            ``MaskedMTKNPT(mask=(1, 1, 1))``
npt         ``berendsen``, isotropic        ``NPTBerendsen``
npt         ``berendsen``, anisotropic      ``Inhomogeneous_NPTBerendsen``
                                            (axis-aligned orthorhombic cells)
npt         ``parrinello-rahman``;          refused
            ``barostat_coupling='in-plane'``
==========  ==============================  ====================================

The MTK integrators thermostat with their own Nose-Hoover chain and the
Berendsen barostat with a Berendsen thermostat, so an NPT ``thermostat`` must
match the barostat (an unset barostat follows the thermostat). None of these
needs a triangular cell: the isotropic integrators scale the whole cell
matrix, so ASE's lower-triangular standard cells and upper-triangular ones
are equally accepted.

**Frozen atoms and the centre of mass.** ``FixAtoms`` is honoured through
ASE's own constraint machinery: velocities are drawn with the constraint
applied (frozen atoms start at rest), the constraint-aware integrators keep
frozen atoms exactly in place, and the temperature degrees of freedom exclude
them. ASE 3.29's ``NoseHooverChainNVT`` keeps private momenta that ignore
``FixAtoms`` and froze the free atoms of a half-frozen slab to 0 K
(measured), so that pairing is refused. Without frozen atoms the
centre-of-mass momentum is removed at initialisation, and the integrators
that apply ASE constraints (``VelocityVerlet``, ``Langevin``,
``NVTBerendsen``, ``Bussi``) run with an ``ase.constraints.FixCom`` attached
-- ASE's replacement for the deprecated ``Langevin(fixcm=True)`` -- so their
degree-of-freedom count, and therefore the Bussi/Berendsen targets and the
reported temperature, is ``3N - 3``. The MTK integrators refuse ASE
constraints, and they and ``NoseHooverChainNVT`` thermostat against ``3N``;
they run without ``FixCom`` and report ``3N``. ``FixCom`` belongs to the run
only: it is not written to the trajectory frames and is removed from the
returned atoms.
"""
from __future__ import annotations

import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, probe
from ..errors import ConfigError, MissingDependencyError, ResultError
from ..results import PotentialResult
# The shared damping defaults, re-exported: every engine resolves unset
# damping times through SimulationSpec.resolved_*_damping_fs.
from ..specs import (  # noqa: F401
    DEFAULT_BAROSTAT_DAMPING_FS,
    DEFAULT_THERMOSTAT_DAMPING_FS,
    SimulationSpec,
)
from ..units import ASE
from .base import EngineRuntime

#: NVT thermostat -> (module, class) of the ASE integrator that implements it.
NVT_INTEGRATORS: dict[str, tuple[str, str]] = {
    "langevin": ("ase.md.langevin", "Langevin"),
    "berendsen": ("ase.md.nvtberendsen", "NVTBerendsen"),
    "nose-hoover": ("ase.md.nose_hoover_chain", "NoseHooverChainNVT"),
    "csvr": ("ase.md.bussi", "Bussi"),
}
#: (barostat, coupling) -> (module, class) of the ASE NPT integrator.
NPT_INTEGRATORS: dict[tuple[str, str], tuple[str, str]] = {
    ("mtk", "isotropic"): ("ase.md.nose_hoover_chain", "IsotropicMTKNPT"),
    ("mtk", "anisotropic"): ("ase.md.nose_hoover_chain", "MaskedMTKNPT"),
    ("berendsen", "isotropic"): ("ase.md.nptberendsen", "NPTBerendsen"),
    ("berendsen", "anisotropic"): ("ase.md.nptberendsen", "Inhomogeneous_NPTBerendsen"),
}
#: The thermostat each ASE barostat integrator carries with it.
NPT_THERMOSTATS = {"mtk": "nose-hoover", "berendsen": "berendsen"}
NVE_INTEGRATOR = ("ase.md.verlet", "VelocityVerlet")
DEFAULT_THERMOSTAT = "langevin"
DEFAULT_BAROSTAT = "mtk"
DEFAULT_BAROSTAT_COUPLING = "isotropic"

#: Integrators that apply ASE constraints through ``set_positions`` /
#: ``set_momenta`` and count degrees of freedom with
#: ``get_number_of_degrees_of_freedom``; they run with ``FixCom`` when no atom
#: is frozen. The others keep private momenta (NHC, MTK) or remove the
#: centre-of-mass momentum themselves (NPTBerendsen, ``fixcm=True``).
CONSTRAINT_AWARE_INTEGRATORS = frozenset({"VelocityVerlet", "Langevin", "NVTBerendsen", "Bussi"})

#: The bulk-modulus guess the Berendsen barostat needs for its
#: compressibility. Only a diagnostic default -- NPT here proves that a route
#: integrates, it is not an equation of state -- and it is recorded on every
#: run. ``engine.options["bulk_modulus_GPa"]`` overrides it (NiO is ~200 GPa).
DEFAULT_BULK_MODULUS_GPA = 100.0
#: The Nose-Hoover chain / MTK settings used (ASE's own defaults), recorded.
NHC_CHAIN = {"tchain": 3, "tloop": 1}
MTK_CHAIN = {"tchain": 3, "pchain": 3, "tloop": 1, "ploop": 1}


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


# ---------------------------------------------------------------------------
# Request -> integrator
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IntegratorPlan:
    """The ASE integrator a request resolves to, before anything is built.

    ``defaults`` names the fields this engine filled in (an unset thermostat,
    barostat or coupling), so a manifest reader can tell a choice from a
    default.
    """

    ensemble: str
    module: str
    integrator: str
    thermostat: str | None = None
    barostat: str | None = None
    barostat_coupling: str | None = None
    defaults: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def constraint_aware(self) -> bool:
        return self.integrator in CONSTRAINT_AWARE_INTEGRATORS

    def as_dict(self) -> dict:
        return {
            "ensemble": self.ensemble,
            "integrator": self.integrator,
            "module": self.module,
            "thermostat": self.thermostat,
            "barostat": self.barostat,
            "barostat_coupling": self.barostat_coupling,
            "defaults_applied": list(self.defaults),
            "notes": list(self.notes),
        }


def resolve_integrator(simulation: SimulationSpec) -> IntegratorPlan:
    """Map an MD request onto one ASE integrator, or raise :class:`ConfigError`.

    Pure: imports nothing and reads no structure. Whether the installed ASE
    has the class, and the structure-dependent refusals, are
    :func:`check_simulation`'s.
    """
    if simulation.task != "md":
        raise ConfigError(f"no integrator applies to task {simulation.task!r}")
    ensemble = simulation.ensemble
    if ensemble == "nve":
        return IntegratorPlan("nve", *NVE_INTEGRATOR)
    if ensemble == "nvt":
        thermostat = simulation.thermostat or DEFAULT_THERMOSTAT
        if thermostat not in NVT_INTEGRATORS:
            raise ConfigError(
                f"thermostat {thermostat!r} is not implemented for the ASE engine; "
                f"available: {', '.join(NVT_INTEGRATORS)}"
            )
        return IntegratorPlan(
            "nvt",
            *NVT_INTEGRATORS[thermostat],
            thermostat=thermostat,
            defaults=() if simulation.thermostat else ("thermostat",),
        )
    if ensemble == "npt":
        return _resolve_npt(simulation)
    raise ConfigError(f"ensemble {ensemble!r} is not implemented for the ASE engine")


def _resolve_npt(simulation: SimulationSpec) -> IntegratorPlan:
    defaults: list[str] = []
    barostat = simulation.barostat
    if barostat == "parrinello-rahman":
        raise ConfigError(
            "barostat = 'parrinello-rahman' is not offered on the ASE engine: ASE's only "
            "Parrinello-Rahman-type integrator is the deprecated MelchionnaNPT "
            "(ase.md.npt.NPT), which the ASE documentation discourages for its "
            "oscillations. Use barostat = 'mtk' (Martyna-Tobias-Klein) or 'berendsen' on "
            "ASE, or run Parrinello-Rahman-style control through LAMMPS."
        )
    if barostat is None:
        # An unset barostat follows the thermostat it must be paired with.
        barostat = "berendsen" if simulation.thermostat == "berendsen" else DEFAULT_BAROSTAT
        defaults.append("barostat")
    if barostat not in NPT_THERMOSTATS:
        raise ConfigError(
            f"barostat {barostat!r} is not implemented for the ASE engine; available: "
            f"{', '.join(NPT_THERMOSTATS)}"
        )
    implied = NPT_THERMOSTATS[barostat]
    if simulation.thermostat is not None and simulation.thermostat != implied:
        raise ConfigError(
            f"on the ASE engine barostat = {barostat!r} integrates with its own "
            f"{implied} thermostat, so thermostat = {simulation.thermostat!r} cannot be "
            f"honoured. Leave the thermostat unset (or set it to {implied!r}), or choose "
            + (
                "barostat = 'berendsen' with thermostat = 'berendsen'."
                if barostat == "mtk"
                else "barostat = 'mtk' with thermostat = 'nose-hoover'."
            )
        )
    if simulation.thermostat is None:
        defaults.append("thermostat")
    coupling = simulation.barostat_coupling
    if coupling == "in-plane":
        raise ConfigError(
            "barostat_coupling = 'in-plane' (barostatting only the periodic in-plane "
            "directions of a slab) is implemented only on the LAMMPS engine. The ASE "
            "engine barostats fully periodic bulk cells only: use coupling 'isotropic' or "
            "'anisotropic', or run the slab through LAMMPS."
        )
    if coupling is None:
        coupling = DEFAULT_BAROSTAT_COUPLING
        defaults.append("barostat_coupling")
    if (barostat, coupling) not in NPT_INTEGRATORS:
        raise ConfigError(
            f"barostat {barostat!r} with coupling {coupling!r} is not implemented for the "
            "ASE engine"
        )
    notes: tuple[str, ...] = ()
    if coupling == "anisotropic" and barostat == "mtk":
        notes = (
            "MaskedMTKNPT(mask=(1, 1, 1)): the a, b and c cell lengths fluctuate "
            "independently; unlike MTKNPT, cell shear is not integrated",
        )
    return IntegratorPlan(
        "npt",
        *NPT_INTEGRATORS[(barostat, coupling)],
        thermostat=implied,
        barostat=barostat,
        barostat_coupling=coupling,
        defaults=tuple(defaults),
        notes=notes,
    )


def bulk_modulus_GPa(options: Mapping | None) -> float:
    """``engine.options["bulk_modulus_GPa"]``, or the diagnostic default."""
    value = (options or {}).get("bulk_modulus_GPa", DEFAULT_BULK_MODULUS_GPA)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
        raise ConfigError(
            f"engine.options.bulk_modulus_GPa must be a positive number; got {value!r}"
        )
    return float(value)


def check_simulation(
    simulation: SimulationSpec, atoms=None, *, options: Mapping | None = None
) -> IntegratorPlan | None:
    """Refuse what the ASE engine cannot integrate as asked; return the plan.

    Called by the ASE bridges' ``check_simulation`` (so ``mlip validate``
    refuses before a manifest exists) and again at the start of
    :func:`run_md`, for any bridge that did not. ``None`` for a task that
    integrates nothing.
    """
    if simulation.task != "md":
        return None
    plan = resolve_integrator(simulation)
    if plan.barostat == "berendsen":
        bulk_modulus_GPa(options)
    _check_ase_has(plan)
    if atoms is None:
        return plan
    from ..structures import fixed_atom_indices

    fixed = fixed_atom_indices(atoms)
    if fixed and len(fixed) == len(atoms):
        raise ConfigError(
            f"every one of the {len(atoms)} atoms is frozen by FixAtoms; there is nothing "
            "to integrate"
        )
    if fixed and plan.integrator == "NoseHooverChainNVT":
        raise ConfigError(
            f"thermostat = 'nose-hoover' with {len(fixed)} frozen atom(s) is refused on the "
            "ASE engine: ASE's NoseHooverChainNVT keeps private momenta that ignore "
            "FixAtoms and thermostats against 3N degrees of freedom, and it drove the free "
            "atoms of a half-frozen Cu slab to 0 K (measured, ASE 3.29). Use thermostat = "
            "'langevin' or 'csvr' (Bussi, which counts the constrained degrees of freedom)."
        )
    if plan.integrator == "Inhomogeneous_NPTBerendsen" and not _is_axis_aligned(atoms):
        raise ConfigError(
            "barostat = 'berendsen' with barostat_coupling = 'anisotropic' uses ASE's "
            "Inhomogeneous_NPTBerendsen, which scales each cell vector by the Cartesian "
            "diagonal stress along x, y, z; that is the stress along the cell vector only "
            "for an axis-aligned orthorhombic cell, and this cell is not one. Use "
            "barostat = 'mtk' (MaskedMTKNPT) or an isotropic coupling."
        )
    return plan


def _check_ase_has(plan: IntegratorPlan) -> None:
    """Feature-detect the integrator in the installed ASE, when ASE is installed."""
    from importlib import import_module

    try:
        module = import_module(plan.module)
    except ImportError:
        try:
            import ase  # noqa: F401
        except ImportError:
            return  # no ASE here at all: availability() reports that
        module = None
    if module is None or not hasattr(module, plan.integrator):
        import ase

        raise ConfigError(
            f"this ASE ({ase.__version__}) has no {plan.module}.{plan.integrator}, which the "
            f"{plan.ensemble} request resolves to; upgrade ASE (the MTK, Bussi and "
            "Nose-Hoover chain integrators need ASE >= 3.26, MaskedMTKNPT a newer one) or "
            "choose another thermostat/barostat"
        )


def _is_axis_aligned(atoms, tolerance: float = 1e-10) -> bool:
    cell = [[float(c) for c in row] for row in atoms.get_cell()]
    scale = max(1.0, max(abs(c) for row in cell for c in row))
    return all(
        abs(cell[i][j]) <= tolerance * scale for i in range(3) for j in range(3) if i != j
    )


# ---------------------------------------------------------------------------
# Single points and optimisation
# ---------------------------------------------------------------------------


def singlepoint(
    atoms,
    calculator,
    *,
    compute_stress: bool = False,
    compute_per_atom_energy: bool = False,
) -> dict:
    """Evaluate one geometry and return canonical quantities as plain data.

    Forces are the calculator's **raw** forces
    (``get_forces(apply_constraint=False)``): ASE's default would zero the
    rows of ``FixAtoms`` atoms, which makes a validation or cross-engine
    comparison disagree with an engine that ignores constraints and biases
    ``max_force``. The constraints are recorded under ``native`` instead,
    with the largest force on a free atom.

    Returns a mapping rather than a :class:`PotentialResult` so the calling
    bridge can attach the provenance fields it alone knows (which
    implementation, which convention, which native units).
    """
    require_ase()
    from ..structures import constraint_summary

    started = time.perf_counter()
    atoms = atoms.copy()
    atoms.calc = calculator
    energy = float(atoms.get_potential_energy(apply_constraint=False))
    forces = [[float(c) for c in row] for row in atoms.get_forces(apply_constraint=False)]
    stress = None
    if compute_stress:
        stress = [float(v) for v in atoms.get_stress(voigt=True, apply_constraint=False)]
    per_atom = None
    if compute_per_atom_energy:
        per_atom = [float(v) for v in atoms.get_potential_energies()]
    summary = constraint_summary(atoms)
    native: dict = {
        "forces": "raw calculator forces: get_forces(apply_constraint=False)",
        "constraints": None,
    }
    if summary["n_fixed"] or summary["unsupported"]:
        frozen = set(summary["fixed_atoms"])
        native["constraints"] = {
            "fixed_atoms": summary["fixed_atoms"],
            "n_fixed": summary["n_fixed"],
            "other": summary["unsupported"],
            "applied_to_forces": False,
        }
        native["max_force_free_atoms_eV_per_A"] = max(
            (
                sum(c * c for c in row) ** 0.5
                for index, row in enumerate(forces)
                if index not in frozen
            ),
            default=0.0,
        )
    return {
        "energy_eV": energy,
        "forces_eV_per_A": forces,
        "stress_eV_per_A3": stress,
        "per_atom_energy_eV": per_atom,
        "symbols": list(atoms.get_chemical_symbols()),
        "wall_time_s": time.perf_counter() - started,
        "trajectory_path": None,
        "native": native,
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


# ---------------------------------------------------------------------------
# Molecular dynamics
# ---------------------------------------------------------------------------


def run_md(
    atoms,
    calculator,
    simulation: SimulationSpec,
    *,
    workdir: Path,
    trajectory_name: str = "smoke_md.traj",
    options: Mapping | None = None,
) -> dict:
    """Run a short trajectory and report what a diagnostic run should report.

    The point of this function is not production dynamics; it is to prove that
    a route integrates stably. Everything it reports is read from the run
    itself: the integrator's step counter, the frames counted back from the
    written trajectory (step 0 included; a missing, unreadable or short file
    raises :class:`~nio_md_prep.mlip.errors.ResultError`), one energy sample
    per logged MD step (taken only by an observer, so no step is sampled
    twice), the conserved quantity where the integrator has one, the
    temperature degrees of freedom, and the integrator with its resolved
    parameters and seed (drawn and recorded when ``simulation.seed`` is
    ``None``).

    ``options`` is the engine's ``options`` mapping (``bulk_modulus_GPa``).
    """
    require_ase()
    import numpy as np
    from ase import units as ase_units
    from ase.constraints import FixCom
    from ase.io.trajectory import Trajectory

    from ..specs import resolve_seed
    from ..structures import fixed_atom_indices

    plan = check_simulation(simulation, atoms, options=options)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    atoms = atoms.copy()
    atoms.calc = calculator
    input_constraints = list(atoms.constraints)
    fixed = fixed_atom_indices(atoms)

    # One seed, two independent streams: velocities and the integrator's noise.
    seed = resolve_seed(simulation.seed)
    velocity_stream, dynamics_stream = np.random.SeedSequence(seed).spawn(2)
    rng_velocities = np.random.default_rng(velocity_stream)
    rng_dynamics = np.random.default_rng(dynamics_stream)

    policy = centre_of_mass_policy(plan, len(atoms), len(fixed))
    use_fixcom = policy["fixcom"]
    if use_fixcom:
        atoms.set_constraint(input_constraints + [FixCom()])
    velocities = _initialise_velocities(
        atoms, simulation.temperature_K, rng_velocities, stationary=policy["stationary"]
    )
    com_removal = policy["description"]
    ndof = policy["temperature_ndof"]
    ase_ndof = int(atoms.get_number_of_degrees_of_freedom())
    if ase_ndof != ndof:
        raise ResultError(
            f"ASE counts {ase_ndof} degrees of freedom but the shared convention gives "
            f"{ndof}; refusing to report temperatures that mean different things"
        )

    timestep = simulation.timestep_fs * ase_units.fs
    dynamics, resolved = _make_dynamics(
        atoms, simulation, plan, timestep, rng_dynamics, options=options
    )
    conserved_name, conserved = _conserved_energy(dynamics, atoms)

    trajectory_path = workdir / trajectory_name
    log_path = workdir / "smoke_md.log"
    interval = max(1, simulation.trajectory_interval)
    trajectory = Trajectory(str(trajectory_path), "w")

    def write_frame():
        # Frames carry the structure's own constraints (FixAtoms), not the
        # run's FixCom, so every frame reloads as a valid input structure.
        run_constraints = atoms.constraints
        atoms.set_constraint(list(input_constraints))
        try:
            trajectory.write(atoms)
        finally:
            atoms.set_constraint(run_constraints)

    # One row per logged MD step, taken only by this observer: ASE calls it
    # once at nsteps == 0 and then after every step whose number is a
    # multiple of the interval, so step numbers are unique by construction.
    # The final step is sampled once more after the run when it is not a
    # multiple, so the energy diagnostics always span the whole run.
    history: list[tuple[int, float, float, float, float | None]] = []

    def sample():
        potential = float(atoms.get_potential_energy())
        kinetic = float(atoms.get_kinetic_energy())
        history.append(
            (
                int(dynamics.nsteps),
                potential,
                kinetic,
                float(atoms.get_temperature()),
                float(conserved()) if conserved else None,
            )
        )

    dynamics.attach(write_frame, interval=interval)
    dynamics.attach(sample, interval=max(1, simulation.log_interval))

    started = time.perf_counter()
    try:
        dynamics.run(simulation.steps)
    finally:
        trajectory.close()
    wall_time = time.perf_counter() - started
    steps_completed = int(dynamics.nsteps)
    if not history or history[-1][0] != steps_completed:
        sample()
    final_temperature = float(atoms.get_temperature())
    atoms.set_constraint(input_constraints)

    steps = [row[0] for row in history]
    if any(b <= a for a, b in zip(steps, steps[1:])):
        raise ResultError(f"the energy log repeats or reorders MD steps: {steps}")
    frames_written = verify_frames(trajectory_path, steps_completed, interval)

    has_conserved = conserved is not None
    log_path.write_text(
        "# step time_fs potential_eV kinetic_eV total_eV temperature_K"
        + (" conserved_eV" if has_conserved else "")
        + f"   (temperatures over {ndof} degrees of freedom)\n"
        + "".join(
            f"{step} {step * simulation.timestep_fs:.6f} {potential:.10f} {kinetic:.10f} "
            f"{potential + kinetic:.10f} {temp:.6f}"
            + (f" {value:.10f}" if has_conserved else "")
            + "\n"
            for step, potential, kinetic, temp, value in history
        ),
        encoding="utf-8",
    )
    series: dict = {
        "time_fs": [step * simulation.timestep_fs for step in steps],
        "total_energy_eV": [potential + kinetic for _, potential, kinetic, _, _ in history],
    }
    if has_conserved:
        series["conserved_energy_eV"] = [row[4] for row in history]
        series["conserved_quantity"] = conserved_name

    temperatures = [row[3] for row in history]
    resolved.update(
        {
            "seed": seed,
            "seed_source": "simulation.seed" if simulation.seed is not None else "drawn",
            "rng": (
                "numpy SeedSequence(seed).spawn(2): stream 0 draws the initial "
                "velocities, stream 1 the integrator's noise"
            ),
            "velocity_initialisation": velocities,
            "centre_of_mass": com_removal,
            "temperature_ndof": ndof,
            "conserved_quantity": conserved_name,
            "ase_todict": _jsonable(dynamics.todict()),
        }
    )
    constraints = None
    if fixed:
        constraints = {
            "fixed_atoms": list(fixed),
            "n_fixed": len(fixed),
            "method": (
                "ASE FixAtoms: velocities drawn with the constraint applied, positions "
                "and momenta of frozen atoms held by the integrator, frozen atoms "
                "excluded from the temperature degrees of freedom"
            ),
        }
    return {
        "atoms": atoms,
        "steps_completed": steps_completed,
        "frames_written": frames_written,
        "trajectory_path": str(trajectory_path),
        "log_path": str(log_path),
        "final_positions": atoms.get_positions().tolist(),
        "final_cell": atoms.get_cell().tolist(),
        "final_pbc": [bool(v) for v in atoms.pbc],
        "energy_series": series,
        "temperature_start_K": temperatures[0] if temperatures else None,
        "temperature_end_K": final_temperature,
        "max_temperature_K": max(temperatures + [final_temperature]),
        "temperature_ndof": ndof,
        "constraints": constraints,
        "wall_time_s": wall_time,
        "integrator": plan.integrator,
        "integrator_resolved": resolved,
    }


def _initialise_velocities(atoms, temperature_K, rng, *, stationary: bool) -> str:
    """Draw Maxwell-Boltzmann momenta with the constraints applied; describe it."""
    from ase.md import velocitydistribution

    if not temperature_K:
        if not atoms.has("momenta"):
            return "none: atoms start at rest"
        # Input momenta, projected through the constraints like drawn ones.
        atoms.set_momenta(atoms.get_momenta(), apply_constraint=True)
        return "momenta carried by the input structure, constraints applied"
    thermalize = getattr(velocitydistribution, "thermalize_momenta", None)
    if thermalize is not None:  # ASE >= 3.29, where MaxwellBoltzmannDistribution is deprecated
        thermalize(atoms, temperature_K=temperature_K, rng=rng)
        how = "ase thermalize_momenta"
    else:
        velocitydistribution.MaxwellBoltzmannDistribution(
            atoms, temperature_K=temperature_K, rng=rng
        )
        how = "ase MaxwellBoltzmannDistribution"
    if stationary:
        velocitydistribution.Stationary(atoms)
        how += " + Stationary (centre-of-mass momentum removed, temperature preserved)"
    return f"{how} at {temperature_K:g} K, constraints applied"


def integrator_settings(
    simulation: SimulationSpec, plan: IntegratorPlan, *, options: Mapping | None = None
) -> tuple[dict, dict]:
    """What the planned integrator is constructed with, and the record of it.

    Returns ``(kwargs, resolved)``: ``kwargs`` are the keyword arguments
    passed to the ASE class after ``(atoms, timestep)``, in ASE's native
    units (time in ``ase.units.fs`` multiples, pressure and compressibility in
    eV/Angstrom^3 and its inverse); ``resolved`` records the same values in
    fs / bar / GPa for a manifest. The integrator's random generator, when
    it takes one, is added by :func:`_make_dynamics`. Pure apart from the
    lazy ``ase.units`` import, so :func:`execution_plan` shows exactly the
    numbers :func:`run_md` will pass.
    """
    from ase import units as ase_units

    from ..units import EV_PER_ANGSTROM3_IN_BAR

    resolved: dict = {
        **plan.as_dict(),
        "class": f"{plan.module}.{plan.integrator}",
        "timestep_fs": float(simulation.timestep_fs),
        "temperature_K": simulation.temperature_K,
    }
    if plan.ensemble == "nve":
        return {}, resolved

    damping_fs = simulation.resolved_thermostat_damping_fs
    tau = damping_fs * ase_units.fs
    temperature = simulation.temperature_K
    resolved["thermostat_damping_fs"] = damping_fs
    if plan.integrator == "Langevin":
        # ASE's Langevin friction is a rate: 1/tau, in inverse ASE time.
        resolved.update(friction_per_fs=1.0 / damping_fs, fixcm=False)
        return dict(temperature_K=temperature, friction=1.0 / tau, fixcm=False), resolved
    if plan.integrator == "NVTBerendsen":
        resolved.update(taut_fs=damping_fs, fixcm=False)
        return dict(temperature_K=temperature, taut=tau, fixcm=False), resolved
    if plan.integrator == "NoseHooverChainNVT":
        resolved.update(tdamp_fs=damping_fs, **NHC_CHAIN)
        return dict(temperature_K=temperature, tdamp=tau, **NHC_CHAIN), resolved
    if plan.integrator == "Bussi":
        resolved.update(taut_fs=damping_fs)
        return dict(temperature_K=temperature, taut=tau), resolved

    pressure_au = simulation.pressure_bar / EV_PER_ANGSTROM3_IN_BAR
    barostat_fs = simulation.resolved_barostat_damping_fs
    resolved.update(
        pressure_bar=float(simulation.pressure_bar),
        pressure_eV_per_A3=pressure_au,
        barostat_damping_fs=barostat_fs,
    )
    if plan.integrator in ("IsotropicMTKNPT", "MaskedMTKNPT"):
        kwargs = dict(
            temperature_K=temperature,
            pressure_au=pressure_au,
            tdamp=tau,
            pdamp=barostat_fs * ase_units.fs,
            **MTK_CHAIN,
        )
        resolved.update(tdamp_fs=damping_fs, pdamp_fs=barostat_fs, **MTK_CHAIN)
        if plan.integrator == "MaskedMTKNPT":
            kwargs["mask"] = (True, True, True)
            resolved["mask"] = [True, True, True]
        return kwargs, resolved

    # Berendsen barostat: needs a compressibility, i.e. a bulk-modulus guess.
    bulk_modulus = bulk_modulus_GPa(options)
    compressibility_au = 1.0 / (bulk_modulus * ase_units.GPa)
    resolved.update(
        taut_fs=damping_fs,
        taup_fs=barostat_fs,
        bulk_modulus_GPa=bulk_modulus,
        bulk_modulus_source=(
            "engine.options.bulk_modulus_GPa"
            if "bulk_modulus_GPa" in (options or {})
            else "DEFAULT_BULK_MODULUS_GPA (a diagnostic default)"
        ),
        compressibility_per_GPa=1.0 / bulk_modulus,
        compressibility_A3_per_eV=compressibility_au,
        fixcm=True,
    )
    kwargs = dict(
        temperature_K=temperature,
        pressure_au=pressure_au,
        taut=tau,
        taup=barostat_fs * ase_units.fs,
        compressibility_au=compressibility_au,
    )
    if plan.integrator == "Inhomogeneous_NPTBerendsen":
        kwargs["mask"] = (1, 1, 1)
        resolved["mask"] = [1, 1, 1]
    return kwargs, resolved


#: Integrators whose constructor takes the run's random generator.
STOCHASTIC_INTEGRATORS = frozenset({"Langevin", "Bussi"})


def _make_dynamics(
    atoms,
    simulation: SimulationSpec,
    plan: IntegratorPlan,
    timestep: float,
    rng,
    *,
    options: Mapping | None = None,
):
    """Build the planned ASE integrator; return it and its resolved parameters."""
    from importlib import import_module

    cls = getattr(import_module(plan.module), plan.integrator)
    kwargs, resolved = integrator_settings(simulation, plan, options=options)
    if plan.integrator in STOCHASTIC_INTEGRATORS:
        kwargs["rng"] = rng
    return cls(atoms, timestep, **kwargs), resolved


def centre_of_mass_policy(plan: IntegratorPlan, n_atoms: int, n_fixed: int) -> dict:
    """How the centre of mass is treated for ``plan``, and the resulting DOF.

    ``fixcom``: an ``ase.constraints.FixCom`` is attached for the run (no
    frozen atoms, constraint-aware integrator); ``stationary``: the
    centre-of-mass momentum is removed once at initialisation.
    """
    from ..diagnostics import temperature_ndof

    fixcom = not n_fixed and plan.constraint_aware
    if n_fixed:
        description = "none: frozen atoms pin the system, the centre of mass is not free"
    elif fixcom:
        description = (
            "ase.constraints.FixCom during the run (centre-of-mass position and momentum "
            "held; not written to trajectory frames)"
        )
    else:
        description = (
            f"Stationary at initialisation only; {plan.integrator} "
            + (
                "removes the centre-of-mass momentum itself (fixcm=True)"
                if plan.integrator == "NPTBerendsen"
                else "conserves the zero centre-of-mass momentum"
            )
            + " and counts 3N degrees of freedom"
        )
    return {
        "fixcom": fixcom,
        "stationary": not n_fixed and not fixcom,
        "temperature_ndof": temperature_ndof(n_atoms, n_fixed=n_fixed, com_removed=fixcom),
        "description": description,
    }


def _conserved_energy(dynamics, atoms):
    """The integrator's conserved quantity as ``(name, callable)``, or ``(None, None)``.

    The NHC and MTK integrators expose ``get_conserved_energy()`` (thermostat
    and, for NPT, barostat and ``P V`` terms included). For Bussi the energy
    the thermostat has put in is ``transferred_energy``, so ``E_total -
    transferred_energy`` is conserved. Langevin and Berendsen have none.
    """
    name = type(dynamics).__name__
    if name == "Bussi" and hasattr(dynamics, "transferred_energy"):
        return (
            "Bussi: total energy minus transferred_energy",
            lambda: atoms.get_potential_energy()
            + atoms.get_kinetic_energy()
            - dynamics.transferred_energy,
        )
    getter = getattr(dynamics, "get_conserved_energy", None)
    if callable(getter):
        return f"{name}.get_conserved_energy", getter
    return None, None


def verify_frames(path: Path, steps_completed: int, interval: int) -> int:
    """Count the trajectory's frames and require exactly the expected number.

    ``steps_completed // interval + 1`` (step 0 included). A missing,
    unreadable or short trajectory raises
    :class:`~nio_md_prep.mlip.errors.ResultError` -- it is a failed run, not
    a successful zero-frame one.
    """
    frames = _count_frames(path)
    expected = steps_completed // interval + 1
    if frames != expected:
        raise ResultError(
            f"{path} holds {frames} frames; {steps_completed} steps written every "
            f"{interval} steps (including step 0) should give {expected}"
        )
    return frames


def _count_frames(path: Path) -> int:
    """Frames in an ASE trajectory, read back from disk. Unreadable is an error."""
    from ase.io.trajectory import Trajectory

    if not Path(path).exists():
        raise ResultError(f"the trajectory {path} was not written")
    try:
        with Trajectory(str(path), "r") as handle:
            count = len(handle)
            if count:
                handle[count - 1]  # the last frame must decode, not only the index
            return count
    except Exception as exc:
        raise ResultError(f"the trajectory {path} cannot be read: {exc}") from exc


def _jsonable(value):
    """``dyn.todict()`` with numpy scalars/arrays turned into plain Python."""
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        return _jsonable(tolist())
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return repr(value)


def integrator_parameters(
    simulation: SimulationSpec, options: Mapping | None = None
) -> dict | None:
    """What :func:`run_md` will build for ``simulation``, for a manifest.

    The same :func:`integrator_settings` record the run returns (minus what
    only the run knows: the seed actually used and ``todict()``). Refusals
    are reported, not raised: ``engine_parameters`` must describe a job even
    when validation would refuse it.
    """
    if simulation.task != "md":
        return None
    try:
        plan = resolve_integrator(simulation)
        _, report = integrator_settings(simulation, plan, options=options)
    except ConfigError as exc:
        return {"refused": str(exc)}
    except ImportError:  # no ASE here: describe the plan without unit values
        plan = resolve_integrator(simulation)
        report = {**plan.as_dict(), "class": f"{plan.module}.{plan.integrator}"}
    return report


#: The native time unit of every ASE integrator argument.
ASE_TIME_UNIT = "ASE time unit (Angstrom*sqrt(amu/eV), 1 fs = ase.units.fs = 0.0982269 of it)"
ASE_PRESSURE_UNIT = "eV/Angstrom^3"


def dynamics_plan(
    simulation: SimulationSpec, atoms=None, *, options: Mapping | None = None
) -> dict | None:
    """The ``dynamics`` block of a bridge's execution plan, or ``None`` (no MD).

    Every time is given in fs and in the ASE-native value the integrator is
    constructed with (``timestep_native``, ``thermostat_damping_native``,
    ``barostat_damping_native``; ``friction_native`` for Langevin). With a
    structure, the centre-of-mass treatment and the temperature degrees of
    freedom are resolved too. Raises :class:`ConfigError` for a refused
    request, exactly as :func:`check_simulation` would.
    """
    if simulation.task != "md":
        return None
    from ase import units as ase_units

    plan = check_simulation(simulation, atoms, options=options)
    kwargs, resolved = integrator_settings(simulation, plan, options=options)
    fs = ase_units.fs
    dynamics: dict = {
        "ensemble": plan.ensemble,
        "integrator": f"{plan.module}.{plan.integrator}",
        "thermostat": plan.thermostat,
        "barostat": plan.barostat,
        "barostat_coupling": plan.barostat_coupling,
        "timestep_fs": float(simulation.timestep_fs),
        "timestep_native": float(simulation.timestep_fs) * fs,
        "thermostat_damping_fs": None,
        "thermostat_damping_native": None,
        "barostat_damping_fs": None,
        "barostat_damping_native": None,
        "native_time_unit": ASE_TIME_UNIT,
        "defaults_applied": list(plan.defaults),
        "notes": list(plan.notes),
        "temperature_K": simulation.temperature_K,
        "seed": simulation.seed,
        "seed_source": (
            "simulation.seed"
            if simulation.seed is not None
            else "drawn at run time (secrets) and recorded in integrator_resolved"
        ),
        "rng": (
            "numpy SeedSequence(seed).spawn(2): stream 0 draws the initial velocities, "
            "stream 1 the integrator's noise"
        ),
        "constructor_kwargs_native": _jsonable(kwargs),
        "resolved": resolved,
    }
    if plan.ensemble != "nve":
        damping = simulation.resolved_thermostat_damping_fs
        dynamics.update(thermostat_damping_fs=damping, thermostat_damping_native=damping * fs)
        if plan.integrator == "Langevin":
            dynamics["friction_native"] = 1.0 / (damping * fs)
    if plan.ensemble == "npt":
        damping = simulation.resolved_barostat_damping_fs
        dynamics.update(
            barostat_damping_fs=damping,
            barostat_damping_native=damping * fs,
            pressure_bar=float(simulation.pressure_bar),
            pressure_native=resolved["pressure_eV_per_A3"],
            pressure_native_unit=ASE_PRESSURE_UNIT,
        )
    if atoms is not None:
        from ..structures import fixed_atom_indices

        policy = centre_of_mass_policy(plan, len(atoms), len(fixed_atom_indices(atoms)))
        dynamics.update(
            centre_of_mass=policy["description"],
            temperature_ndof=policy["temperature_ndof"],
        )
    return dynamics


def refused_dynamics(simulation: SimulationSpec, reason: str) -> dict:
    """The ``dynamics`` block for an MD request no ASE integrator will run.

    Same keys as :func:`dynamics_plan`: the request as written, no
    integrator and no native values (nothing would be constructed), and the
    reason, which is also in the plan's ``unmet_capabilities``.
    """
    return {
        "ensemble": simulation.ensemble,
        "integrator": None,
        "thermostat": simulation.thermostat,
        "barostat": simulation.barostat,
        "barostat_coupling": simulation.barostat_coupling,
        "timestep_fs": simulation.timestep_fs,
        "timestep_native": None,
        "thermostat_damping_fs": simulation.thermostat_damping_fs,
        "thermostat_damping_native": None,
        "barostat_damping_fs": simulation.barostat_damping_fs,
        "barostat_damping_native": None,
        "native_time_unit": ASE_TIME_UNIT,
        "refused": reason,
    }


def base_execution_plan(bridge, simulation: SimulationSpec, atoms=None) -> dict:
    """The route-independent part of an ASE bridge's ``execution_plan``.

    Runs every check ``validate`` runs -- element coverage, capability
    negotiation, the bridge's ``check_simulation`` -- but collects the
    refusals into ``unmet_capabilities`` instead of raising, so the plan can
    be shown for a job that would be refused. Nothing is executed and no
    calculator is built. The bridge fills in ``model_checkpoint``,
    ``energy_convention``, ``device``, ``precision`` and ``model_hashes``.
    """
    from ..capabilities import check_elements
    from ..errors import MlipError

    unmet: list[str] = []
    capabilities = bridge.capabilities()
    if atoms is not None:
        try:
            check_elements(
                capabilities, atoms.get_chemical_symbols(), label=bridge.potential.label
            )
        except MlipError as exc:
            unmet.append(str(exc))
    try:
        unmet += bridge.requirements(simulation, atoms).unmet(capabilities)
    except MlipError as exc:
        unmet.append(str(exc))
    try:
        bridge.check_simulation(simulation, atoms)
    except MlipError as exc:
        unmet.append(str(exc))
    dynamics = None
    if simulation.task == "md":
        try:
            dynamics = dynamics_plan(simulation, atoms, options=bridge.engine.options)
        except MlipError as exc:
            dynamics = refused_dynamics(simulation, str(exc))
    availability = bridge.availability()
    elements = (
        sorted(capabilities.elements)
        if capabilities.elements is not None
        else sorted(bridge.potential.elements)
    )
    return {
        "potential_kind": bridge.potential_kind,
        "engine": bridge.engine_kind,
        "bridge": f"{type(bridge).__module__}.{type(bridge).__qualname__}",
        "implementation": bridge.implementation,
        "model_checkpoint": None,
        "exported_model": None,
        "elements": elements,
        "energy_convention": None,
        "units": {
            "native": ASE.name,
            "pressure_unit": ASE_PRESSURE_UNIT,
            "energy": "eV",
            "length": "Angstrom",
            "force": "eV/Angstrom",
            "time": ASE_TIME_UNIT,
            "conversion": "none: ASE's native units are the canonical ones",
        },
        "device": None,
        "precision": None,
        "dynamics": dynamics,
        "lammps": None,
        "openmm": None,
        "availability": {
            "available": bool(availability),
            "missing": list(availability.missing),
            "detail": availability.detail,
        },
        "unmet_capabilities": unmet,
        "model_hashes": {"declared": dict(bridge.potential.declared_hashes()), "observed": {}},
        "task": simulation.task,
        "forces": "raw calculator forces: get_forces(apply_constraint=False)",
    }


def restate(
    result: PotentialResult,
    target: str | None,
    *,
    atomic_reference_energies,
    source: str,
) -> PotentialResult:
    """Report ``result`` under ``target``, recording the conversion in ``extras``.

    The ASE routes compute in one convention (the calculator's) and report in
    the one the job asked for; the calculator's own number, the E0s' origin
    and the direction of the conversion stay in ``extras`` so the two can be
    told apart in a manifest.
    """
    if target is None or target == result.energy_convention:
        return result
    from dataclasses import replace

    converted = result.in_convention(target, atomic_reference_energies=atomic_reference_energies)
    extras = {
        **result.extras,
        "native_energy_convention": result.energy_convention,
        "native_energy_eV": result.energy_eV,
        "energy_conversion": (
            f"{target} = {result.energy_convention} "
            f"{'-' if target == 'interaction' else '+'} sum of atomic reference energies (E0)"
        ),
        "atomic_reference_energies_source": source,
    }
    return replace(converted, extras=extras)


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
    "IntegratorPlan",
    "NVT_INTEGRATORS",
    "NPT_INTEGRATORS",
    "NPT_THERMOSTATS",
    "CONSTRAINT_AWARE_INTEGRATORS",
    "DEFAULT_THERMOSTAT",
    "DEFAULT_BAROSTAT",
    "DEFAULT_BAROSTAT_COUPLING",
    "DEFAULT_BULK_MODULUS_GPA",
    "require_ase",
    "resolve_integrator",
    "check_simulation",
    "integrator_parameters",
    "integrator_settings",
    "base_execution_plan",
    "dynamics_plan",
    "refused_dynamics",
    "centre_of_mass_policy",
    "singlepoint",
    "optimize",
    "run_md",
    "verify_frames",
    "restate",
    "build_result",
]
