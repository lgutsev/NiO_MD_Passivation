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
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..errors import CapabilityError, ConfigError, EnergyConventionError
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
        return metadata

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

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        """What this job will request, recorded without loading the model."""
        convention = self.requested_convention(simulation)
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
            "platform_properties_requested": openmm_engine.check_request(
                self.engine, device=self.potential.device
            ),
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
            thermostat = openmm_engine.resolved_thermostat(simulation)
            parameters["integrator"] = (
                openmm_engine.SUPPORTED_THERMOSTATS.get(thermostat, thermostat)
                if thermostat
                else "VerletIntegrator"
            )
            parameters["thermostat"] = thermostat
            if thermostat:
                parameters["thermostat_damping_fs"] = simulation.resolved_thermostat_damping_fs
        return parameters

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
            "model_sha256": metadata["sha256"],
            "model_class": metadata.get("model_class"),
            "model_native_dtype": metadata.get("dtype"),
            "model_evaluation_dtype": self.potential.precision,
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
        from ..structures import fixed_atom_indices

        self.validate(simulation, atoms)
        endpoint = self.as_singlepoint(simulation)
        initial = self.singlepoint(atoms, endpoint)
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
        payload["integrator_resolved"]["system"] = record
        final = self.singlepoint(payload["atoms"], endpoint)
        trajectory = self._trajectory(payload, simulation, initial, final)
        trajectory.extras["trajectory_frame"] = payload.get("trajectory_frame")
        return trajectory


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
