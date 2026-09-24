"""MACE on OpenMM, through OpenMM-ML.

The energy convention is the whole story here. OpenMM-ML's MACE
implementation distinguishes the interaction energy from the energy including
atomic self-energies, and *interaction* is its default. A comparison against
the ASE calculator that ignores this looks catastrophically wrong -- thousands
of eV apart on a few thousand atoms -- while the forces agree to machine
precision, because the difference is a per-element constant with zero
gradient.

So this bridge:

* asks OpenMM-ML for whichever convention the job declares, rather than
  taking the default silently, and records the convention it asked for
  (``engine_parameters(simulation)`` uses the job's simulation, not the
  potential's default);
* uses the **model's own** atomic reference energies (E0) for any
  conversion, because OpenMM-ML's ``interaction_energy`` subtracts exactly
  those; E0s declared in the configuration that disagree are reported;
* reads the model once per file hash, at execution time only -- validation
  never loads torch.

Everything the engine needs that the potential spec owns is passed
explicitly: ``potential.device`` and ``potential.precision`` go to
OpenMM-ML's ``createSystem``, ``atoms.info['charge']`` / ``['spin']`` go to
MACE's ``total_charge`` / ``total_spin`` exactly as ASE's ``MACECalculator``
reads them. OpenMM-ML cannot select a head of a multi-head model, and has no
``compile_mode``; both are refused.

No stress tensor is reported through this route, so ``stress`` is false and a
constant-pressure job is rejected during capability negotiation (and, by
policy, in :meth:`MaceOpenMMBridge.check_simulation`).

Precision is the subtle one. OpenMM-ML 1.6/1.7 cast the *inputs* to the
requested dtype but never the model's parameters, so ``potential.precision``
is only achievable when the model file is stored in that dtype (1.8 converts
the model). Validation reads the stored dtype statically from the pickled
model -- no torch import -- and refuses a known mismatch; execution checks it
again against the loaded model. :meth:`MaceOpenMMBridge.execution_plan`
reports ``guaranteed = False`` whenever it cannot establish the outcome.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..errors import CapabilityError, ConfigError, EnergyConventionError, MlipError, ResultError
from ..potentials.mace import MaceAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import SimulationSpec
from ..units import INTERACTION, OPENMM
from ..engines import openmm_engine
from ..engines.openmm_engine import ENERGY_TYPE, OpenMMEngine
from .base import Bridge

#: What OpenMM-ML gives you if nobody says otherwise.
OPENMM_ML_DEFAULT_CONVENTION = INTERACTION

#: Declared and model E0s closer than this (eV) are the same reference.
E0_CONFLICT_TOLERANCE_EV = 1e-8

#: MACE model classes whose output has no ``interaction_energy`` key (the
#: plain ``mace.modules.models.MACE``; ScaleShiftMACE adds it). Refused at
#: validate for the interaction convention when the stored class is known.
NO_INTERACTION_ENERGY_CLASSES = frozenset({"MACE"})

#: torch storage class -> dtype name, for the static model read.
_STORAGE_DTYPES = {
    "FloatStorage": "float32",
    "DoubleStorage": "float64",
    "HalfStorage": "float16",
    "BFloat16Storage": "bfloat16",
}


class MaceOpenMMBridge(Bridge):
    potential_kind = "mace"
    engine_kind = "openmm"
    implementation = "openmm-ml"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = MaceAdapter(potential)
        self.runtime = OpenMMEngine(engine)
        self._systems: dict[tuple, tuple] = {}
        self._metadata: dict[str, dict] = {}
        self._peeks: dict[tuple, dict] = {}

    # -- conventions ------------------------------------------------------

    def requested_convention(self, simulation: SimulationSpec | None = None) -> str:
        """Which convention this run will ask OpenMM-ML for.

        The job's explicit request wins; then the potential's own convention,
        if this route can produce it; otherwise OpenMM-ML's default, which is
        the interaction energy.
        """
        if simulation is not None and simulation.energy_convention:
            return simulation.energy_convention
        if self.potential.energy_convention in ENERGY_TYPE:
            return self.potential.energy_convention
        return OPENMM_ML_DEFAULT_CONVENTION

    def native_convention(self) -> str:
        return self.requested_convention()

    def capabilities(self) -> CapabilitySet:
        """The route's capabilities, derived without loading the model.

        OpenMM-ML produces either convention directly (``returnEnergyType``),
        so both are reachable whether or not E0s are known.
        """
        engine_capabilities = self.runtime.capabilities()
        route = replace(
            engine_capabilities,
            gpu=engine_capabilities.gpu and self.potential.wants_gpu,
            native_energy_convention=self.native_convention(),
            convertible_energy_conventions=frozenset(ENERGY_TYPE),
            native_units=OPENMM.name,
            notes=engine_capabilities.notes
            + (
                f"OpenMM-ML returnEnergyType={ENERGY_TYPE[self.native_convention()]!r}",
                "OpenMM-ML's own default is the interaction energy; this route states "
                "the convention explicitly instead of inheriting it",
            ),
        )
        return self.adapter.capabilities().intersect(route)

    def availability(self) -> Availability:
        model = self.adapter.availability()
        if not model:
            return model
        return self.runtime.availability()

    # -- model metadata (execution time only) ------------------------------

    def model_metadata(self) -> dict:
        """Elements, cutoff, E0s, class, heads and dtype, read once per file hash."""
        from ..provenance import sha256_file

        digest = sha256_file(self.potential.model_path)
        if digest not in self._metadata:
            self._metadata[digest] = {"sha256": digest, **_read_model(self.potential.model_path)}
        return self._metadata[digest]

    def stored_model(self) -> dict:
        """Class and stored float dtype of the model file, read without torch.

        Cached per file identity (path, size, mtime). ``{}`` when the file
        does not exist.
        """
        path = self.potential.model_path
        try:
            stat = path.stat()
        except OSError:
            return {}
        key = (str(path), stat.st_size, stat.st_mtime_ns)
        if key not in self._peeks:
            self._peeks[key] = _peek_model(path)
        return self._peeks[key]

    def precision_plan(self, *, loaded: bool = False) -> dict:
        """:func:`openmm_engine.precision_plan` for this model and the installed OpenMM-ML.

        ``loaded=True`` uses the dtype of the model as torch loads it
        (execution time); otherwise the statically read stored dtype.
        """
        if loaded:
            dtype = self.model_metadata().get("dtype")
            source = "the model's parameters, loaded with torch"
        else:
            stored = self.stored_model()
            dtype = stored.get("dtype")
            source = stored.get("source")
        return openmm_engine.precision_plan(
            self.potential.precision,
            openmmml=openmm_engine.openmmml_version(),
            model_dtype=dtype,
            model_dtype_source=source,
        )

    def device_plan(self) -> dict:
        return openmm_engine.device_plan(
            self.potential.device,
            platform=self.engine.platform,
            openmmml=openmm_engine.openmmml_version(),
        )

    def atomic_reference_energies(self):
        """The model's own E0s: what OpenMM-ML's interaction energy subtracts.

        ``None`` when the model file cannot be read here (then this route
        cannot run either) or its E0 layout is not recognised.
        """
        if not self.availability():
            return None
        e0s = self.model_metadata().get("atomic_reference_energies_eV")
        return dict(e0s) if e0s else None

    def _e0_conflicts(self, model_e0s) -> list[str]:
        declared = self.potential.atomic_reference_energies
        if not declared or not model_e0s:
            return []
        conflicts = []
        for symbol, value in sorted(dict(declared).items()):
            if symbol in model_e0s and abs(model_e0s[symbol] - value) > E0_CONFLICT_TOLERANCE_EV:
                conflicts.append(
                    f"{symbol}: declared {value!r} eV, model {model_e0s[symbol]!r} eV"
                )
        return conflicts

    def _check_model(self, convention: str) -> dict:
        """Refuse model features OpenMM-ML cannot evaluate as asked."""
        metadata = self.model_metadata()
        heads = metadata.get("heads") or []
        if len(heads) > 1:
            raise CapabilityError(
                f"{self.potential.model_path.name} is a multi-head MACE model (heads: "
                f"{', '.join(heads)}); OpenMM-ML evaluates the first head and cannot "
                "select one. Use the ASE engine with potential.head, or export a "
                "single-head model."
            )
        if convention == INTERACTION and not metadata.get("scale_shift"):
            raise CapabilityError(
                f"{metadata.get('model_class')} does not output an interaction energy "
                "(only ScaleShiftMACE does), so OpenMM-ML cannot report the interaction "
                "convention for it; request energy_convention = 'total'"
            )
        precision = self.precision_plan(loaded=True)
        if not precision["guaranteed"]:
            raise CapabilityError(
                precision.get("refusal")
                or f"potential.precision = {self.potential.precision!r} cannot be "
                f"guaranteed on this route: {precision['note']}"
            )
        device = self.device_plan()
        if not device["guaranteed"]:  # pragma: no cover - require_openmm refuses first
            raise CapabilityError(f"potential.device cannot be guaranteed: {device['note']}")
        return {**metadata, "precision_plan": precision, "device_plan": device}

    # -- gating -----------------------------------------------------------

    def check_simulation(self, simulation: SimulationSpec, atoms=None) -> None:
        """Refuse what this route cannot honour, without importing OpenMM."""
        if self.potential.head is not None:
            raise CapabilityError(
                f"potential.head = {self.potential.head!r}: OpenMM-ML's MACE route has no "
                "head selection (it evaluates the model's first head). Use the ASE engine."
            )
        if self.potential.compile_mode is not None:
            raise ConfigError(
                f"potential.compile_mode = {self.potential.compile_mode!r} is not supported "
                "by OpenMM-ML's MACE route, which would ignore it; remove it"
            )
        openmm_engine.check_request(self.engine, device=self.potential.device)
        openmm_engine.check_md_request(simulation)
        if atoms is not None:
            openmm_engine.check_structure(atoms, md=simulation.task == "md")
        # Statically read (no torch): refuses a precision this OpenMM-ML
        # release is known not to deliver for this model file.
        refusal = self.precision_plan().get("refusal")
        if refusal:
            raise CapabilityError(refusal)
        stored_class = self.stored_model().get("model_class")
        if (
            self.requested_convention(simulation) == INTERACTION
            and stored_class in NO_INTERACTION_ENERGY_CLASSES
        ):
            raise CapabilityError(
                f"{stored_class} does not output an interaction energy (only "
                "ScaleShiftMACE does), so OpenMM-ML cannot report the interaction "
                "convention for it; request energy_convention = 'total'"
            )

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        """What this job will request, recorded without loading the model.

        ``model_precision`` and ``device`` are the static plans (with their
        ``guaranteed`` flags); what actually ran is in the results'
        ``execution_settings``.
        """
        convention = self.requested_convention(simulation)
        try:
            properties = openmm_engine.check_request(self.engine, device=self.potential.device)
        except ConfigError as exc:
            properties = {"refused": str(exc)}
        parameters = {
            "implementation": self.implementation,
            "potential": "openmmml.MLPotential('mace')",
            "model_path": str(self.potential.model_path),
            "returnEnergyType": ENERGY_TYPE[convention],
            "energy_convention": convention,
            "createSystem": {
                "returnEnergyType": ENERGY_TYPE[convention],
                "precision": openmm_engine.MODEL_PRECISION[self.potential.precision],
                "device": self.potential.device,
                "removeCMMotion": (
                    "False for single points; for MD True unless FixAtoms atoms are present"
                ),
                "charge": "atoms.info['charge'] (default 0)",
                "multiplicity": "atoms.info['spin'] (default 1), as ASE's MACECalculator",
            },
            "platform": self.engine.platform,
            "platform_precision": self.engine.platform_precision,
            "threads": self.engine.threads,
            "precision": self.engine.precision,
            "platform_properties_requested": properties,
            "model_precision": self.precision_plan(),
            "device": self.device_plan(),
            "native_units": OPENMM.name,
            "native_energy_unit": OPENMM.energy,
            "native_length_unit": OPENMM.length,
            "energy_scale_kJ_per_mol_per_eV": openmm_engine.OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV,
            "atomic_reference_energies_source": (
                "the model file (read at execution); OpenMM-ML's interaction energy "
                "subtracts the model's own E0s"
            ),
        }
        if simulation is not None and simulation.task == "md":
            dynamics = openmm_engine.dynamics_plan(simulation)
            parameters["integrator"] = dynamics["integrator"]
            parameters["thermostat"] = dynamics["thermostat"]
            parameters["dynamics"] = dynamics
            if dynamics["thermostat"]:
                parameters["thermostat_damping_fs"] = dynamics["thermostat_damping_fs"]
        return parameters

    def execution_plan(self, simulation: SimulationSpec, atoms=None) -> dict:
        """The complete plan for this job on this route, without executing anything.

        Imports neither torch nor OpenMM: the model is hashed and its class
        and stored dtype are read statically, the OpenMM-ML version comes from
        package metadata. Refusals are listed in ``unmet_capabilities`` (this
        method reports; :meth:`validate` refuses). Values OpenMM decides only
        when a Context exists (an unnamed platform, its property defaults)
        are ``None`` with a note, and ``guaranteed`` is true only where the
        installed OpenMM-ML is known to apply the request.
        """
        from ..diagnostics import temperature_ndof
        from ..provenance import sha256_file
        from ..specs import path_key
        from ..structures import fixed_atom_indices, periodic_axes

        unmet: list[str] = []

        def attempt(check, *args, **kwargs):
            try:
                return check(*args, **kwargs)
            except MlipError as exc:
                message = str(exc)
                if message not in unmet:
                    unmet.append(message)
                return None

        convention = self.requested_convention(simulation)
        capabilities = self.capabilities()
        model_path = Path(self.potential.model_path)
        exists = model_path.exists()
        digest = sha256_file(model_path) if exists else None
        stored = self.stored_model() if exists else {}
        md = simulation.task == "md"

        if atoms is not None:
            missing = capabilities.supports_elements(atoms.get_chemical_symbols())
            if missing:
                unmet.append(
                    f"elements {', '.join(sorted(missing))} are not covered by "
                    f"{self.potential.label}"
                )
        requirements = attempt(self.requirements, simulation, atoms)
        if requirements is not None:
            unmet.extend(p for p in requirements.unmet(capabilities) if p not in unmet)
        attempt(self.check_simulation, simulation, atoms)

        box = None
        fixed: list[int] = []
        charge = multiplicity = None
        sources: tuple = (None, None)
        if atoms is not None:
            fixed = list(attempt(fixed_atom_indices, atoms) or ())
            if all(periodic_axes(atoms)):
                box = attempt(openmm_engine.reduced_box, atoms.get_cell())
            charge, multiplicity, sources = _charge_and_spin(atoms)
        remove_cm_motion = bool(md and not fixed)
        properties = attempt(
            openmm_engine.check_request, self.engine, device=self.potential.device
        )

        dynamics = openmm_engine.dynamics_plan(simulation)
        if dynamics is not None:
            seed_policy = (
                "velocity and integrator seeds are two SeedSequence children of the "
                "recorded seed"
            )
            if simulation.seed is None:
                seed_policy += "; the seed is drawn at execution and recorded"
            velocities = "zero"
            if simulation.temperature_K:
                velocities = "setVelocitiesToTemperature(T, velocity_seed)"
                if remove_cm_motion:
                    velocities += ", then the centre-of-mass velocity is removed"
            dynamics = {
                **dynamics,
                "removeCMMotion": remove_cm_motion,
                "temperature_ndof": (
                    temperature_ndof(
                        len(atoms), n_fixed=len(fixed), com_removed=remove_cm_motion
                    )
                    if atoms is not None
                    else None
                ),
                "velocities": velocities,
                "seed": simulation.seed,
                "seed_policy": seed_policy,
                "trajectory": {
                    "format": "DCD",
                    "interval_steps": simulation.trajectory_interval,
                    "expected_frames": simulation.steps // simulation.trajectory_interval + 1,
                    "includes_step_0": True,
                },
            }

        platform = self.engine.platform
        if exists:
            precision = self.precision_plan()
        else:
            precision = openmm_engine.precision_plan(
                self.potential.precision,
                openmmml=openmm_engine.openmmml_version(),
                model_dtype=None,
                model_dtype_source="the model file does not exist",
            )
        precision["platform_precision"] = openmm_engine.platform_precision_plan(
            platform, self.engine.platform_precision
        )
        declared = self.potential.declared_hashes()
        observed = {str(model_path): digest} if digest else {}
        observed_by_key = {path_key(p): sha for p, sha in observed.items()}
        if platform is None:
            platform_note = (
                "OpenMM picks the fastest available platform when the Context is "
                "created; the choice is read back and recorded"
            )
        else:
            platform_note = (
                "named platform; whether it is installed is checked when OpenMM is "
                "imported at execution (an unavailable one is refused with the list of "
                "available platforms)"
            )
        return {
            "potential_kind": self.potential_kind,
            "engine": self.engine_kind,
            "bridge": type(self).__name__,
            "implementation": self.implementation,
            "label": self.label,
            "model_checkpoint": {
                "path": str(model_path),
                "exists": exists,
                "sha256": digest,
                "model_class": stored.get("model_class"),
                "stored_dtype": stored.get("dtype"),
                "stored_dtype_source": stored.get("source"),
            },
            "exported_model": None,
            "exported_model_note": (
                "none: OpenMM-ML loads the MACE checkpoint itself (torch.load in "
                "MACEPotentialImpl.addForces); there is no export step"
            ),
            "elements": sorted(self.potential.elements),
            "structure_elements": (
                sorted(set(atoms.get_chemical_symbols())) if atoms is not None else None
            ),
            "energy_convention": convention,
            "energy_convention_detail": {
                "requested": convention,
                "returnEnergyType": ENERGY_TYPE[convention],
                "openmmml_default": ENERGY_TYPE[OPENMM_ML_DEFAULT_CONVENTION],
                "atomic_reference_energies_source": "the model's own E0s (read at execution)",
            },
            "units": {
                "native": OPENMM.name,
                "energy": OPENMM.energy,
                "length": OPENMM.length,
                "time": "ps",
                "pressure_unit": None,
                "pressure_note": "no pressure or stress is produced on this route (NPT refused)",
                "energy_scale_kJ_per_mol_per_eV": (
                    openmm_engine.OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV
                ),
                "length_scale_A_per_nm": openmm_engine.OPENMMML_LENGTH_SCALE_A_PER_NM,
            },
            "device": self.device_plan(),
            "precision": precision,
            "dynamics": dynamics,
            "lammps": None,
            "openmm": {
                "platform": {"requested": platform, "effective": platform, "note": platform_note},
                "properties": properties,
                "create_system_kwargs": {
                    "returnEnergyType": ENERGY_TYPE[convention],
                    "precision": openmm_engine.MODEL_PRECISION[self.potential.precision],
                    "device": self.potential.device,
                    "charge": charge,
                    "multiplicity": multiplicity,
                    "charge_source": sources[0],
                    "multiplicity_source": sources[1],
                    "removeCMMotion": remove_cm_motion,
                },
                "masses": "ASE masses (atoms.get_masses()); FixAtoms atoms zero",
                "zero_mass_atoms": fixed,
                "box": box.as_dict() if box is not None else None,
                "versions": {
                    "openmm": openmm_engine.distribution_version("openmm"),
                    "openmmml": openmm_engine.openmmml_version(),
                    "source": "package metadata (nothing imported)",
                },
            },
            "availability": self.availability().as_dict(),
            "unmet_capabilities": unmet,
            "model_hashes": {
                "declared": dict(declared),
                "observed": observed,
                "match": (
                    all(
                        observed_by_key.get(path_key(p), "").lower() == sha.lower()
                        for p, sha in declared.items()
                    )
                    if declared
                    else None
                ),
            },
        }

    # -- execution --------------------------------------------------------

    def system(self, atoms, simulation: SimulationSpec | None = None, *, md: bool = False):
        """Build (and cache) the OpenMM System for this structure and request.

        Keyed by everything OpenMM-ML bakes into the System at ``createSystem``
        time -- convention, element order, periodicity, cell (the box and
        one-hot table live inside the force) -- plus the masses (zero for
        FixAtoms atoms), the CMMotionRemover, dtype, device, charge and
        multiplicity. Returns ``(system, record, box)``.
        """
        from ..structures import fixed_atom_indices, periodic_axes

        convention = self.requested_convention(simulation)
        fixed = fixed_atom_indices(atoms)
        remove_cm_motion = bool(md and not fixed)
        charge, multiplicity, sources = _charge_and_spin(atoms)
        pbc = periodic_axes(atoms)
        key = (
            convention,
            tuple(atoms.get_chemical_symbols()),
            pbc,
            tuple(float(x) for x in atoms.get_cell().array.flat) if all(pbc) else None,
            tuple(float(m) for m in atoms.get_masses()),
            fixed,
            remove_cm_motion,
            self.potential.precision,
            self.potential.device,
            charge,
            multiplicity,
        )
        if key not in self._systems:
            self.require_available()
            box = openmm_engine.check_structure(atoms, md=md)
            system, _, record = openmm_engine.build_system(
                atoms,
                self.potential.model_path,
                energy_convention=convention,
                precision=self.potential.precision,
                device=self.potential.device,
                remove_cm_motion=remove_cm_motion,
                fixed_atoms=fixed,
                charge=charge,
                multiplicity=multiplicity,
                box=box,
            )
            record["charge_source"], record["multiplicity_source"] = sources
            self._systems[key] = (system, record, box)
        return self._systems[key]

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        self.validate(simulation, atoms)
        convention = self.requested_convention(simulation)
        if convention not in ENERGY_TYPE:
            raise EnergyConventionError(
                f"OpenMM-ML cannot report the {convention!r} energy convention"
            )
        self.require_available()
        metadata = self._check_model(convention)
        system, record, box = self.system(atoms, simulation)
        payload = openmm_engine.singlepoint(
            atoms, system, self.engine, box=box, device=self.potential.device
        )
        model_e0s = metadata.get("atomic_reference_energies_eV")
        extras = {
            **payload["native"],
            "returnEnergyType": ENERGY_TYPE[convention],
            "system": record,
            "execution_settings": _execution_settings(record, payload["native"]["platform"], metadata),
            "model_sha256": metadata["sha256"],
            "model_class": metadata.get("model_class"),
            "model_native_dtype": metadata.get("dtype"),
            "model_evaluation_dtype": metadata["precision_plan"]["effective"],
            "model_precision": metadata["precision_plan"],
            "device": metadata["device_plan"],
            "atomic_reference_energies_eV": dict(model_e0s) if model_e0s else None,
            "atomic_reference_energies_source": "model",
        }
        conflicts = self._e0_conflicts(model_e0s)
        if conflicts:
            extras["declared_e0_conflicts"] = conflicts
        return self._result(
            payload, energy_convention=convention, native_units=OPENMM.name, **extras
        )

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        """A smoke trajectory whose endpoints and dynamics share one execution setup.

        The endpoint single points and the dynamics use the same platform,
        platform properties, model device and dtype (``execution_settings``,
        read back from each Context); a difference is a
        :class:`~nio_md_prep.mlip.errors.ResultError`, not a footnote.
        """
        from ..structures import fixed_atom_indices

        self.validate(simulation, atoms)
        endpoint = self.as_singlepoint(simulation)
        initial = self.singlepoint(atoms, endpoint)
        metadata = self._check_model(self.requested_convention(endpoint))
        fixed = fixed_atom_indices(atoms)
        system, record, box = self.system(atoms, simulation, md=True)
        payload = openmm_engine.run_md(
            atoms,
            system,
            self.engine,
            simulation,
            workdir=Path(workdir),
            box=box,
            device=self.potential.device,
            fixed_atoms=fixed,
            remove_cm_motion=record["createSystem"]["removeCMMotion"],
        )
        resolved = payload["integrator_resolved"]
        resolved["system"] = record
        resolved["execution_settings"] = _execution_settings(
            record, resolved["platform"], metadata
        )
        final = self.singlepoint(payload["atoms"], endpoint)
        for label, point in (("initial", initial), ("final", final)):
            if point.extras["execution_settings"] != resolved["execution_settings"]:
                raise ResultError(
                    f"the {label} single point ran with different execution settings than "
                    f"the dynamics: {point.extras['execution_settings']} vs "
                    f"{resolved['execution_settings']}"
                )
        trajectory = self._trajectory(payload, simulation, initial, final)
        trajectory.extras["trajectory_frame"] = payload.get("trajectory_frame")
        trajectory.extras["final_velocities"] = payload.get("final_velocities")
        trajectory.extras["execution_settings_identical_for_endpoints_and_dynamics"] = True
        return trajectory


def _execution_settings(record: dict, platform: dict, metadata: dict) -> dict:
    """Everything that decides how the model and OpenMM evaluate, as it actually ran.

    ``createSystem`` arguments as passed (without ``removeCMMotion``, which
    differs by design between single points and MD), the platform name and
    every property value read back from the Context, and the model's
    device/dtype as the installed OpenMM-ML applies them.
    """
    create = {k: v for k, v in record["createSystem"].items() if k != "removeCMMotion"}
    return {
        "createSystem": create,
        "platform": platform["name"],
        "platform_properties": dict(platform["properties"]),
        "model_device": metadata["device_plan"]["effective"],
        "model_dtype": metadata["precision_plan"]["effective"],
        "openmmml": record["versions"].get("openmmml"),
        "openmm": record["versions"].get("openmm"),
    }


def _charge_and_spin(atoms) -> tuple[float, float, tuple[str, str]]:
    """MACE's total charge and spin inputs, read as ASE's MACECalculator reads them."""
    info = getattr(atoms, "info", {}) or {}
    charge = float(info["charge"]) if "charge" in info else 0.0
    spin = float(info["spin"]) if "spin" in info else 1.0
    return (
        charge,
        spin,
        (
            "atoms.info['charge']" if "charge" in info else "default 0",
            "atoms.info['spin']" if "spin" in info else "default 1",
        ),
    )


def _peek_model(path: Path) -> dict:
    """The model class and stored float dtype of a ``torch.save``d model, without torch.

    ``torch.save`` writes a zip archive whose ``data.pkl`` names every tensor
    storage class (``torch FloatStorage`` ...) and the model's class as pickle
    globals; :mod:`pickletools` lists them without executing anything. A
    model whose floating-point storages are not all one dtype, or a file
    that is not in that format, reports ``dtype = None`` (unknown), never a
    guess.
    """
    import pickletools
    import zipfile

    record: dict = {
        "model_class": None,
        "dtype": None,
        "storage_dtypes": [],
        "source": "static read of the pickled model (zip data.pkl globals; torch not imported)",
    }
    try:
        with zipfile.ZipFile(path) as archive:
            names = [n for n in archive.namelist() if n == "data.pkl" or n.endswith("/data.pkl")]
            if not names:
                record["error"] = "no data.pkl in the archive"
                return record
            data = archive.read(names[0])
        strings: list[str] = []
        dtypes: set[str] = set()
        for opcode, argument, _ in pickletools.genops(data):
            if opcode.name in ("SHORT_BINUNICODE", "BINUNICODE", "UNICODE", "BINUNICODE8"):
                strings.append(argument)
                continue
            if opcode.name == "GLOBAL":
                module, _, name = str(argument).partition(" ")
            elif opcode.name == "STACK_GLOBAL" and len(strings) >= 2:
                module, name = strings[-2], strings[-1]
            else:
                continue
            if module == "torch" and name in _STORAGE_DTYPES:
                dtypes.add(_STORAGE_DTYPES[name])
            elif (
                record["model_class"] is None
                and module.startswith("mace.modules")
                and name.endswith("MACE")
            ):
                record["model_class"] = name
    except (OSError, zipfile.BadZipFile, ValueError, EOFError) as exc:
        record["error"] = f"{type(exc).__name__}: {exc}"
        return record
    record["storage_dtypes"] = sorted(dtypes)
    if len(dtypes) == 1:
        record["dtype"] = next(iter(dtypes))
    elif dtypes:
        record["error"] = f"floating-point storages of several dtypes: {sorted(dtypes)}"
    return record


def _read_model(path: Path) -> dict:
    """Load the model on the CPU once and describe what OpenMM-ML will see."""
    import mace  # noqa: F401 - sets TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD for e3nn
    import torch
    from mace.modules import ScaleShiftMACE

    from ..potentials.mace import _model_atomic_energies, _model_cutoff, _model_elements

    model = torch.load(str(path), map_location="cpu", weights_only=False)
    elements = _model_elements(model)
    heads = getattr(model, "heads", None)
    parameters = next(iter(model.parameters()), None)
    return {
        "elements": elements,
        "cutoff_angstrom": _model_cutoff(model),
        "atomic_reference_energies_eV": _model_atomic_energies(model, elements),
        "model_class": type(model).__name__,
        "scale_shift": isinstance(model, ScaleShiftMACE),
        "heads": [str(h) for h in heads] if heads is not None else [],
        "dtype": str(parameters.dtype).replace("torch.", "") if parameters is not None else None,
    }


__all__ = ["MaceOpenMMBridge", "OPENMM_ML_DEFAULT_CONVENTION", "E0_CONFLICT_TOLERANCE_EV"]
