"""MACE on LAMMPS, through ML-IAP or through a MACE pair style.

Both routes run the *exported* model, not the checkpoint the ASE calculator
loads. MACE's exporter (``mace_create_lammps_model``) appends a suffix to the
checkpoint's *full* file name (``nio.model`` -> ``nio.model-mliap_lammps.pt``
or ``nio.model-lammps.pt``; see :data:`EXPORT_SUFFIX`), and that is where the
export is looked for unless ``engine.options.exported_model_path`` names it.
The export is staged into the job directory and hashed, the pair commands
name the staged file, and -- where torch can read the file -- its dtype,
element table, head and MACE runtime flags are compared with the request
before anything runs.

``mliap`` (preferred)
    ``pair_style mliap unified <export> 0``. Needs a LAMMPS build with the
    ML-IAP *and* PYTHON packages. MACE models with more than one interaction
    layer exchange ghost-atom features between layers through
    ``forward_exchange``, which only the KOKKOS ML-IAP coupling provides, so
    the default coupling is KOKKOS for every device:

    * ``device = "cuda"``: ``-k on g 1 -sf kk -pk kokkos newton on neigh half``
      (MACE's documented arguments) on a KOKKOS build with a GPU back end,
      ``activate_mliappy_kokkos`` in the python route;
    * ``device = "cpu"``: ``-k on t N -sf kk -pk kokkos newton on neigh half``
      on a host-only KOKKOS build. MACE refuses CPU tensors under KOKKOS
      unless its ``MACE_ALLOW_CPU`` flag is true, and mace-torch reads that
      flag when the export is *created* and pickles it into the file (setting
      it when LAMMPS runs has no effect, measured with mace-torch 0.3.16), so
      the export's stored value is what is checked, and nothing is set in the
      LAMMPS environment for it.

    ``engine.options.mliap_coupling = "plain"`` selects the non-KOKKOS python
    coupling (``activate_mliappy``) on the CPU, for single-layer models only.

``pair-mace`` (legacy)
    ``pair_style mace`` with a TorchScript export. MACE's deprecation table
    announces that v1.0 drops that artifact; every run records the warning.
    ``no_domain_decomposition`` is rendered only for a single-rank run.
    ``pair_mace`` puts the model on CUDA whenever LibTorch sees a GPU, so a
    CPU request hides the GPUs from the LAMMPS process
    (``CUDA_VISIBLE_DEVICES=""``); the device it reports on stdout is checked
    after an executable-route run.

Neither route claims a per-atom energy decomposition: the MACE site energies
LAMMPS would tally have not been demonstrated against a reference here.
"""
from __future__ import annotations

import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, module_available
from ..errors import ConfigError, MlipError, ModelIntegrityError
from ..potentials.mace import MaceAdapter
from ..results import PotentialResult, TrajectoryResult
from ..specs import LammpsMlipPotentialSpec, SimulationSpec
from ..units import LAMMPS_METAL
from ..engines import lammps_engine
from ..engines.lammps_engine import LammpsEngine
from .base import Bridge

#: Suffix MACE's exporter appends to the checkpoint's *full* path. mace-torch
#: 0.3.16, ``mace/cli/create_lammps_model.py:107``:
#: ``torch.save(lammps_model, model_path + "-mliap_lammps.pt")`` (``--format
#: mliap``) and ``:110``: ``lammps_model_compiled.save(model_path +
#: "-lammps.pt")`` (``--format libtorch``, the default).
EXPORT_SUFFIX = {"mliap": "-mliap_lammps.pt", "pair-mace": "-lammps.pt"}
#: The exporter's ``--format`` value for each route.
EXPORT_FORMAT = {"mliap": "mliap", "pair-mace": "libtorch"}
#: ``create_lammps_model.py`` line that writes each route's file (mace-torch 0.3.16).
_EXPORT_SOURCE_LINE = {"mliap": 107, "pair-mace": 110}

#: ``engine.options`` keys read by these routes.
EXPORTED_MODEL_OPTION = "exported_model_path"
#: The older spelling of :data:`EXPORTED_MODEL_OPTION`, still accepted.
LEGACY_EXPORTED_MODEL_OPTION = "model_path"
EXPORTED_SHA256_OPTION = "exported_model_sha256"
MLIAP_COUPLING_OPTION = "mliap_coupling"
MLIAP_COUPLINGS = ("kokkos", "plain")

#: KOKKOS arguments after ``-k on ...`` for both ML-IAP KOKKOS routes (MACE docs).
KOKKOS_SUFFIX_ARGS = ("-sf", "kk", "-pk", "kokkos", "newton", "on", "neigh", "half")

PAIR_MACE_LEGACY_WARNING = (
    "pair-mace is a legacy route: MACE's deprecation table announces that v1.0 drops "
    "the LAMMPS_MACE TorchScript wrapper and the -lammps.pt artifact; the ML-IAP export "
    "is MACE's one supported LAMMPS artifact"
)
MLIAP_PLAIN_NOTE = (
    "the plain (non-KOKKOS) ML-IAP python coupling has no forward_exchange, which MACE "
    "needs for every interaction layer after the first (mace/tools/utils.py LAMMPS_MP): "
    "only single-layer MACE models run on it"
)
MLIAP_ALLOW_CPU_NOTE = (
    "MACE refuses CPU tensors under KOKKOS unless MACE_ALLOW_CPU is true; mace-torch reads "
    "MACE_ALLOW_CPU/MACE_FORCE_CPU in LAMMPS_MLIAP_MACE.__init__ (lammps_mliap_mace.py "
    "MACELammpsConfig) when the export is created and pickles them into the file, so they "
    "must be set when exporting; in the LAMMPS run environment they have no effect "
    "(measured, mace-torch 0.3.16) and are not set"
)
PER_ATOM_NOTE = (
    "per-atom (site) energies are not claimed for MACE on LAMMPS: pair_mliap and pair_mace "
    "tally them, but that has not been demonstrated against a reference for these exports"
)

_EXPORT_CACHE: dict[tuple[str, str, str], dict] = {}


def default_export_path(checkpoint: Path, implementation: str) -> Path:
    """Where ``mace_create_lammps_model <checkpoint>`` writes this route's export."""
    checkpoint = Path(checkpoint)
    return checkpoint.with_name(checkpoint.name + EXPORT_SUFFIX[implementation])


def inspect_export(path: Path, implementation: str) -> dict:
    """What the export file says about itself, when torch can read it here.

    ML-IAP exports are pickles of ``LAMMPS_MLIAP_MACE`` (``torch.load``; the
    mace package must be importable), pair-mace exports are TorchScript
    (``torch.jit.load``). Read: dtype (``r_max``), element table, number of
    interaction layers, exported head index, and for ML-IAP the MACE
    runtime flags stored at export time. Never raises: an unreadable file
    is reported with the reason, and nothing is marked verified. Cached by
    path, SHA256 and route.
    """
    from ..provenance import sha256_file

    path = Path(path)
    report: dict = {"readable": False, "path": str(path), "reader": None, "note": ""}
    if not path.is_file():
        report["note"] = "the export file does not exist"
        return report
    sha = sha256_file(path)
    key = (str(path.resolve()), sha, implementation)
    if key in _EXPORT_CACHE:
        return dict(_EXPORT_CACHE[key])
    needed = ("torch", "mace") if implementation == "mliap" else ("torch",)
    missing = [name for name in needed if not module_available(name)]
    if missing:
        report["note"] = (
            f"{', '.join(missing)} not importable in {sys.executable}; the export's dtype, "
            "elements and flags were not read"
        )
        return report
    try:
        import torch
        from ase.data import chemical_symbols

        if implementation == "mliap":
            report["reader"] = "torch.load(weights_only=False)"
            wrapper = torch.load(str(path), map_location="cpu", weights_only=False)
            inner = wrapper.model
            config = getattr(wrapper, "config", None)
            report.update(
                cls=type(wrapper).__name__,
                dtype=str(getattr(wrapper, "dtype", inner.r_max.dtype)).replace("torch.", ""),
                elements=list(getattr(wrapper, "element_types", [])),
                num_interactions=int(inner.num_interactions),
                head_index=int(inner.head.reshape(-1)[0]),
                heads=list(getattr(getattr(inner, "model", None), "heads", []) or []),
                allow_cpu=getattr(config, "allow_cpu", None),
                force_cpu=getattr(config, "force_cpu", None),
            )
        else:
            report["reader"] = "torch.jit.load"
            module = torch.jit.load(str(path), map_location="cpu")
            buffers = dict(module.named_buffers(recurse=False))
            try:
                heads = list(module.model.heads)
            except Exception:  # not every TorchScript export keeps the list
                heads = None
            report.update(
                cls=type(module).__name__,
                dtype=str(buffers["r_max"].dtype).replace("torch.", ""),
                elements=[chemical_symbols[int(z)] for z in buffers["atomic_numbers"].tolist()],
                num_interactions=int(buffers["num_interactions"].reshape(-1)[0]),
                head_index=int(buffers["head"].reshape(-1)[0]),
                heads=heads,
            )
        report["readable"] = True
        report["sha256"] = sha
    except Exception as exc:  # a corrupt or unexpected export
        report["readable"] = False
        report["note"] = f"could not read the export: {type(exc).__name__}: {exc}"
        return report
    _EXPORT_CACHE[key] = dict(report)
    return report


class _MaceLammpsBridge(Bridge):
    """Shared machinery: render a pair style from a MACE spec, then run LAMMPS."""

    potential_kind = "mace"
    engine_kind = "lammps"
    #: LAMMPS units this route runs in. MACE is an eV/Angstrom model, so metal.
    units = "metal"

    def __init__(self, potential, engine, *, registration=None) -> None:
        super().__init__(potential, engine, registration=registration)
        self.adapter = MaceAdapter(potential)
        self.runtime = LammpsEngine(engine)

    # -- the exported model -----------------------------------------------

    def _override(self) -> tuple[Path | None, str | None]:
        options = self.engine.options
        given = {
            name: options.get(name)
            for name in (EXPORTED_MODEL_OPTION, LEGACY_EXPORTED_MODEL_OPTION)
            if options.get(name)
        }
        if len(given) == 2 and Path(str(given[EXPORTED_MODEL_OPTION])) != Path(
            str(given[LEGACY_EXPORTED_MODEL_OPTION])
        ):
            raise ConfigError(
                f"engine.options.{EXPORTED_MODEL_OPTION} and engine.options."
                f"{LEGACY_EXPORTED_MODEL_OPTION} name different exported models; give one"
            )
        for name in (EXPORTED_MODEL_OPTION, LEGACY_EXPORTED_MODEL_OPTION):
            if name in given:
                return Path(str(given[name])), name
        return None, None

    def default_exported_model_path(self) -> Path:
        return default_export_path(self.potential.model_path, self.implementation)

    def exported_model_path(self) -> Path:
        """The LAMMPS-ready export: the override if given, else MACE's own file name."""
        override, _ = self._override()
        return override if override is not None else self.default_exported_model_path()

    def declared_export_sha256(self) -> str | None:
        value = self.engine.options.get(EXPORTED_SHA256_OPTION)
        if value is None:
            return None
        text = str(value).strip().lower()
        if len(text) != 64 or any(c not in "0123456789abcdef" for c in text):
            raise ConfigError(
                f"engine.options.{EXPORTED_SHA256_OPTION} must be a 64-character hex SHA256; "
                f"got {value!r}"
            )
        return text

    def export_record(self) -> dict:
        """Checkpoint and export, as ``mlip validate`` and the manifest report them."""
        from ..provenance import sha256_file

        exported = self.exported_model_path()
        _, option = self._override()
        checkpoint = self.potential.model_path
        suffix = EXPORT_SUFFIX[self.implementation]
        return {
            "checkpoint": {
                "path": str(checkpoint),
                "exists": checkpoint.is_file(),
                "sha256": sha256_file(checkpoint) if checkpoint.is_file() else None,
                "declared_sha256": self.potential.model_sha256,
            },
            "exported_model": {
                "path": str(exported),
                "exists": exported.is_file(),
                "sha256": sha256_file(exported) if exported.is_file() else None,
                "declared_sha256": self.declared_export_sha256(),
                "implementation": self.implementation,
                "export_format": EXPORT_FORMAT[self.implementation],
                "source": (
                    f"engine.options.{option}"
                    if option
                    else f"MACE's export name: checkpoint path + {suffix!r} (mace-torch 0.3.16 "
                    f"mace/cli/create_lammps_model.py:{_EXPORT_SOURCE_LINE[self.implementation]})"
                ),
                "default_path": str(self.default_exported_model_path()),
                "exporter": self.exporter_command(),
            },
        }

    def exporter_command(self) -> str:
        """The ``mace_create_lammps_model`` call that writes this route's export."""
        checkpoint = self.potential.model_path
        env = ""
        if self.implementation == "mliap" and _cuda_ordinal(self.potential.device) is None:
            # Stored in the export at creation time (see MLIAP_ALLOW_CPU_NOTE).
            env = "MACE_ALLOW_CPU=true "
        argv = [
            "mace_create_lammps_model", str(checkpoint),
            "--format", EXPORT_FORMAT[self.implementation],
            "--dtype", self.potential.precision,
        ]
        if self.potential.head:
            argv += ["--head", self.potential.head]
        return env + subprocess.list2cmdline(argv)

    # -- the pair style this route renders --------------------------------

    def staged_model_name(self) -> str:
        """The export's name inside the job directory (spaces cannot be LAMMPS tokens)."""
        return self.exported_model_path().name.replace(" ", "_")

    def pair_style(self) -> str:
        raise NotImplementedError

    def pair_coeff(self) -> tuple[str, ...]:
        raise NotImplementedError

    def launch(self) -> lammps_engine.LammpsLaunch:
        raise NotImplementedError

    def required_packages(self) -> tuple[str, ...]:
        return tuple(self.registration.packages) if self.registration else ()

    def as_lammps_spec(self) -> LammpsMlipPotentialSpec:
        """Project the MACE spec onto the LAMMPS-native representation.

        A MACE model routed through LAMMPS becomes an ordinary LAMMPS MLIP
        spec, so the deck renderer, the type map, model staging (with the
        declared export hash) and the manifest handle it exactly as they
        handle DeepMD or PACE. Both MACE pair styles require ``newton on``,
        which is rendered explicitly.
        """
        elements = self.potential.elements
        exported = self.exported_model_path()
        declared = self.declared_export_sha256()
        return LammpsMlipPotentialSpec(
            label=self.potential.label,
            pair_style=self.pair_style(),
            pair_coeff=self.pair_coeff(),
            type_map={index + 1: symbol for index, symbol in enumerate(elements)},
            required_packages=self.required_packages(),
            units=self.units,
            atom_style="atomic",
            model_paths=(exported,),
            model_hashes={str(exported): declared} if declared else {},
            framework="mace",
            energy_convention=self.potential.energy_convention,
            atomic_reference_energies=self.potential.atomic_reference_energies,
            newton="on",
        )

    # -- capabilities -----------------------------------------------------

    def gpu_active(self) -> bool:
        """Whether the generated runtime puts the model on a GPU."""
        raise NotImplementedError

    def capabilities(self) -> CapabilitySet:
        """Energy, forces and the global virial; no per-atom energies; GPU only if activated."""
        try:
            gpu = self.gpu_active()
        except ConfigError:
            gpu = False
        engine_capabilities = replace(
            self.runtime.capabilities(), per_atom_energy=False, gpu=gpu
        )
        return self.adapter.capabilities().intersect(
            self._route_capabilities(
                engine_capabilities,
                energy_convention=self.potential.energy_convention,
                native_units=LAMMPS_METAL.name,
                notes=(
                    "runs the exported LAMMPS model, not the checkpoint the ASE "
                    "calculator loads; verify equivalence before production use",
                    PER_ATOM_NOTE,
                ),
            )
        )

    def availability(self) -> Availability:
        exported = self.exported_model_path()
        if not exported.is_file():
            record = self.export_record()["exported_model"]
            return Availability(
                available=False,
                missing=(),
                detail=(
                    f"the {self.implementation} LAMMPS export of this model was not found at "
                    f"{exported} ({record['source']}). Create it with `{record['exporter']}`, "
                    f"or name it with engine.options.{EXPORTED_MODEL_OPTION}."
                ),
            )
        try:
            launch = self.launch()
        except ConfigError as exc:
            return Availability(False, (), str(exc))
        return lammps_engine.launch_availability(launch)

    def atomic_reference_energies(self):
        return self.adapter.atomic_reference_energies()

    # -- validation -------------------------------------------------------

    def check_export(self) -> dict:
        """Compare what the export says about itself with the request.

        Returns the export inspection record. Raises :class:`ConfigError` when
        a readable export disagrees: another dtype than
        ``potential.precision``, missing declared elements, another head, or
        a runtime flag this route needs. An unreadable export is not an
        error here (it is reported as unverified), except that
        ``potential.head`` cannot be honoured without reading it.
        """
        exported = self.exported_model_path()
        info = inspect_export(exported, self.implementation)
        if not info["readable"]:
            if self.potential.head is not None:
                raise ConfigError(
                    f"potential.head = {self.potential.head!r}: a LAMMPS export carries the "
                    "head chosen when it was created (mace_create_lammps_model --head), and "
                    f"this export could not be read to check it ({info['note']})"
                )
            return info
        if info["dtype"] != self.potential.precision:
            raise ConfigError(
                f"the LAMMPS export {exported.name} holds a {info['dtype']} model, but "
                f"potential.precision is {self.potential.precision}; re-export with "
                f"mace_create_lammps_model --dtype {self.potential.precision}, or set "
                f"potential.precision = {info['dtype']!r}"
            )
        missing = sorted(set(self.potential.elements) - set(info["elements"]))
        if missing:
            raise ConfigError(
                f"potential.elements declares {', '.join(missing)}, which the LAMMPS export "
                f"{exported.name} does not cover (export covers: {', '.join(info['elements'])})"
            )
        if self.potential.head is not None:
            heads = info.get("heads")
            if not heads:
                raise ConfigError(
                    f"potential.head = {self.potential.head!r}, but the export's head list "
                    "could not be read, so the exported head cannot be checked"
                )
            if self.potential.head not in heads:
                raise ConfigError(
                    f"potential.head = {self.potential.head!r} is not a head of this model "
                    f"(heads: {', '.join(map(str, heads))})"
                )
            if heads.index(self.potential.head) != info["head_index"]:
                raise ConfigError(
                    f"the LAMMPS export {exported.name} was created for head "
                    f"{heads[info['head_index']]!r}, not potential.head = "
                    f"{self.potential.head!r}; re-export with --head {self.potential.head}"
                )
        return info

    def check_simulation(self, simulation: SimulationSpec, atoms=None) -> None:
        """Launch/device compatibility, the export itself, then the shared LAMMPS refusals."""
        if self.potential.compile_mode is not None:
            raise ConfigError(
                f"potential.compile_mode = {self.potential.compile_mode!r} applies to MACE's "
                "ASE calculator; a LAMMPS export is fixed when it is created, so the setting "
                "would be ignored. Remove it for engine = 'lammps'."
            )
        launch = self.launch()
        self.declared_export_sha256()
        self.check_export()
        lammps_engine.check_request(
            self.as_lammps_spec(), self.engine, simulation, atoms, launch=launch
        )

    # -- description ------------------------------------------------------

    def device_plan(self, launch) -> dict:
        raise NotImplementedError

    def precision_plan(self) -> dict:
        info = inspect_export(self.exported_model_path(), self.implementation)
        observed = info.get("dtype") if info["readable"] else None
        return {
            "requested": self.potential.precision,
            "effective": observed,
            "guaranteed": observed == self.potential.precision,
            "note": (
                f"read from the export ({info['reader']}): LAMMPS evaluates the model in the "
                "dtype it was exported with"
                if info["readable"]
                else f"the export's dtype was not verified: {info['note']}. LAMMPS evaluates "
                "the dtype the export was created with (mace_create_lammps_model --dtype)"
            ),
        }

    def warnings(self) -> tuple[str, ...]:
        return ()

    def engine_parameters(self, simulation: SimulationSpec | None = None) -> dict:
        spec = self.as_lammps_spec()
        try:
            launch = self.launch()
            refused = None
        except ConfigError as exc:
            launch = None
            refused = str(exc)
        record = self.export_record()
        report = {
            "implementation": self.implementation,
            "native_units": LAMMPS_METAL.name,
            "energy_convention": self.potential.energy_convention,
            "lammps_units": self.units,
            "atom_style": spec.atom_style,
            "type_map": {str(k): v for k, v in sorted(spec.type_map.items())},
            "required_packages": list(spec.required_packages),
            "model_checkpoint": record["checkpoint"],
            "exported_model": {
                **record["exported_model"],
                "inspection": inspect_export(self.exported_model_path(), self.implementation),
            },
            "exported_model_path": record["exported_model"]["path"],
            "precision": self.precision_plan(),
            "warnings": list(self.warnings()),
        }
        if launch is not None:
            report["device"] = self.device_plan(launch)
            report.update(
                lammps_engine.describe_request(spec, self.engine, simulation, launch)
            )
        else:
            report["launch"] = {"refused": refused}
            report["pair_commands"] = list(
                lammps_engine.rewrite_model_tokens(spec).render_pair_commands()
            )
        return report

    def execution_plan(self, simulation: SimulationSpec, atoms=None) -> dict:
        """What a run would execute, resolved without running (``mlip validate``)."""
        spec = self.as_lammps_spec()
        launch = None
        launch_error = None
        try:
            launch = self.launch()
        except ConfigError as exc:
            launch_error = str(exc)
        record = self.export_record()
        exported = record["exported_model"]
        inspection = inspect_export(self.exported_model_path(), self.implementation)
        declared = {}
        if self.potential.model_sha256:
            declared[record["checkpoint"]["path"]] = self.potential.model_sha256
        if exported["declared_sha256"]:
            declared[exported["path"]] = exported["declared_sha256"]
        plan = {
            "potential_kind": self.potential_kind,
            "engine": self.engine_kind,
            "bridge": self.label,
            "implementation": self.implementation,
            "model_checkpoint": {
                "path": record["checkpoint"]["path"],
                "sha256": record["checkpoint"]["sha256"],
                "exists": record["checkpoint"]["exists"],
            },
            "exported_model": {
                "path": exported["path"],
                "exists": exported["exists"],
                "sha256": exported["sha256"],
                "implementation": self.implementation,
                "export_format": exported["export_format"],
                "source": exported["source"],
                "default_path": exported["default_path"],
                "exporter": exported["exporter"],
                "inspection": inspection,
            },
            "elements": list(self.potential.elements),
            "type_map": {str(k): v for k, v in sorted(spec.type_map.items())},
            "energy_convention": self.potential.energy_convention,
            "units": lammps_engine.units_plan(self.units),
            "device": (
                self.device_plan(launch)
                if launch is not None
                else {
                    "requested": self.potential.device,
                    "effective": None,
                    "guaranteed": False,
                    "note": f"refused: {launch_error}",
                }
            ),
            "precision": self.precision_plan(),
            "dynamics": lammps_engine.dynamics_plan(
                simulation, units=self.units, atoms=atoms, options=self.engine.options
            ),
            "lammps": lammps_engine.lammps_plan(spec, launch, atoms=atoms),
            "openmm": None,
            "model_hashes": {
                "declared": declared,
                "observed": {
                    record["checkpoint"]["path"]: record["checkpoint"]["sha256"],
                    exported["path"]: exported["sha256"],
                },
            },
            "warnings": list(self.warnings()),
            "deck": lammps_engine.deck_plan(spec, simulation, atoms, options=self.engine.options),
            **lammps_engine.bridge_plan_basics(self, simulation, atoms),
        }
        if launch_error is not None:
            plan["lammps"]["launch_refused"] = launch_error
        return plan

    # -- execution --------------------------------------------------------

    def _run(self, atoms, simulation: SimulationSpec, workdir: Path) -> dict:
        self.require_available()
        info = self.check_export()
        record = self.export_record()
        launch = self.launch()
        payload = lammps_engine.run_job(
            self.as_lammps_spec(),
            atoms,
            simulation,
            engine_spec=self.engine,
            workdir=Path(workdir),
            launch=launch,
        )
        native = payload["native"]
        staged = native["staged_model_files"][0]
        if record["exported_model"]["sha256"] != staged["sha256"]:
            raise ModelIntegrityError(
                Path(staged["staged"]), record["exported_model"]["sha256"], staged["sha256"]
            )
        native["model_checkpoint"] = record["checkpoint"]
        native["exported_model"] = {**record["exported_model"], "inspection": info}
        native["precision"] = self.precision_plan()
        native["warnings"] = list(self.warnings())
        native["device"] = {
            **self.device_plan(launch),
            "observed": self.observed_device(native, Path(workdir)),
        }
        return payload

    def observed_device(self, native: dict, workdir: Path) -> str | None:
        """The device the route reports having used. Default: not reported."""
        return None

    def singlepoint(self, atoms, simulation: SimulationSpec) -> PotentialResult:
        capabilities = self.validate(simulation, atoms)
        workdir = Path(self.engine.options.get("workdir", ".")) / "lammps_singlepoint"
        payload = self._run(atoms, self.as_singlepoint(simulation), workdir)
        return self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=LAMMPS_METAL.name,
        )

    def run_md(self, atoms, simulation: SimulationSpec, *, workdir: Path) -> TrajectoryResult:
        capabilities = self.validate(simulation, atoms)
        initial = self.singlepoint(atoms, self.as_singlepoint(simulation))
        payload = self._run(atoms, simulation, Path(workdir))
        final = self._result(
            payload,
            energy_convention=capabilities.native_energy_convention,
            native_units=LAMMPS_METAL.name,
        )
        return self._trajectory(payload, simulation, initial, final)


def _cuda_ordinal(device: str) -> int | None:
    head, _, ordinal = device.partition(":")
    if head != "cuda":
        return None
    return int(ordinal) if ordinal else 0


class MaceLammpsMliapBridge(_MaceLammpsBridge):
    implementation = "mliap"

    def pair_style(self) -> str:
        # "The 0 is mandatory" (MACE docs): ghostneigh_flag = 0.
        return f"mliap unified {self.staged_model_name()} 0"

    def pair_coeff(self) -> tuple[str, ...]:
        return (f"* * {' '.join(self.potential.elements)}",)

    def coupling(self) -> str:
        value = self.engine.options.get(MLIAP_COUPLING_OPTION, "kokkos")
        if value not in MLIAP_COUPLINGS:
            raise ConfigError(
                f"engine.options.{MLIAP_COUPLING_OPTION} must be one of "
                f"{', '.join(MLIAP_COUPLINGS)}; got {value!r}"
            )
        if value == "plain" and _cuda_ordinal(self.potential.device) is not None:
            raise ConfigError(
                f"engine.options.{MLIAP_COUPLING_OPTION} = 'plain' runs MACE on the CPU (the "
                "non-KOKKOS coupling never moves data to a GPU), but potential.device is "
                f"{self.potential.device!r}; use the default KOKKOS coupling for a GPU"
            )
        return value

    def launch(self) -> lammps_engine.LammpsLaunch:
        """ML-IAP + PYTHON always; KOKKOS with a GPU or a host-only back end by device.

        With ``engine.lammps_args`` given they are used verbatim instead of the
        defaults, but they must still select what the device needs: KOKKOS on
        with the ``kk`` suffix, and ``g`` >= 1 exactly when the device is CUDA.
        """
        device = self.potential.device
        if device.startswith("mps"):
            raise ConfigError(
                "potential.device = 'mps' has no LAMMPS equivalent; use device = 'cpu' or "
                "'cuda' for MACE on LAMMPS"
            )
        required = [("package", "ML-IAP"), ("package", "PYTHON")]
        coupling = self.coupling()
        if coupling == "plain":
            if self.engine.lammps_args and lammps_engine.parse_accelerator_args(
                self.engine.lammps_args
            )["kokkos"]:
                raise ConfigError(
                    f"engine.options.{MLIAP_COUPLING_OPTION} = 'plain' but engine.lammps_args "
                    "enable KOKKOS; drop one of them"
                )
            return lammps_engine.resolve_launch(
                self.engine, activate="mliappy", required=required, notes=(MLIAP_PLAIN_NOTE,)
            )
        ordinal = _cuda_ordinal(device)
        cuda = ordinal is not None
        required += [("package", "KOKKOS"), ("kokkos_backend", "gpu" if cuda else "host-only")]
        user_args = tuple(self.engine.lammps_args)
        threads = self.engine.threads
        if user_args:
            accelerator = lammps_engine.parse_accelerator_args(user_args)
            if not accelerator["kokkos"] or accelerator["suffix"] != "kk":
                raise ConfigError(
                    "the ML-IAP MACE route couples through KOKKOS: engine.lammps_args must "
                    "enable it ('-k', 'on', ..., '-sf', 'kk', ...), or leave lammps_args empty "
                    "for the defaults"
                )
            if cuda != (accelerator["kokkos_gpus"] >= 1):
                raise ConfigError(
                    f"potential.device = {device!r} but engine.lammps_args ask KOKKOS for "
                    f"{accelerator['kokkos_gpus']} GPU(s); they must agree"
                )
            extra: tuple[str, ...] = ()
        elif cuda:
            if ordinal != 0:
                raise ConfigError(
                    f"potential.device = {device!r}: KOKKOS assigns GPUs to MPI ranks itself, "
                    "starting at the first visible device; select the GPU with "
                    "CUDA_VISIBLE_DEVICES, or give the KOKKOS arguments explicitly in "
                    "engine.lammps_args"
                )
            kokkos = ("-k", "on", "g", "1") + (("t", str(threads)) if threads else ())
            extra = (*kokkos, *KOKKOS_SUFFIX_ARGS)
        else:
            extra = ("-k", "on", "t", str(threads or 1), *KOKKOS_SUFFIX_ARGS)
        return lammps_engine.resolve_launch(
            self.engine,
            extra_args=extra,
            activate="mliappy_kokkos",
            required=required,
            threads_in_extra_args=True,
            notes=() if cuda else (MLIAP_ALLOW_CPU_NOTE,),
        )

    def gpu_active(self) -> bool:
        return bool(self.launch().accelerator["gpu"])

    def check_export(self) -> dict:
        info = super().check_export()
        if not info["readable"]:
            return info
        coupling = self.coupling()
        if coupling == "plain" and info["num_interactions"] > 1:
            raise ConfigError(
                f"the export has {info['num_interactions']} interaction layers, which the "
                f"plain ML-IAP coupling cannot evaluate ({MLIAP_PLAIN_NOTE}); use the default "
                f"KOKKOS coupling (drop engine.options.{MLIAP_COUPLING_OPTION})"
            )
        if coupling == "kokkos" and _cuda_ordinal(self.potential.device) is not None:
            if info.get("force_cpu"):
                raise ConfigError(
                    f"potential.device = {self.potential.device!r}, but the export was created "
                    "with MACE_FORCE_CPU=true, which mace-torch stores in the export and which "
                    "keeps the model on the CPU under KOKKOS; re-export without it"
                )
        if coupling == "kokkos" and _cuda_ordinal(self.potential.device) is None:
            if not info.get("force_cpu") and not info.get("allow_cpu"):
                raise ConfigError(
                    "potential.device = 'cpu' on the KOKKOS ML-IAP route, but the export was "
                    "created without MACE_ALLOW_CPU, and mace-torch stores that flag in the "
                    "export: LAMMPS_MLIAP_MACE would refuse the CPU tensors ('GPU requested but "
                    "tensor is on CPU'). Re-export with MACE_ALLOW_CPU=true "
                    f"({self.export_record()['exported_model']['exporter']})"
                )
        return info

    def device_plan(self, launch) -> dict:
        accelerator = launch.accelerator
        if self.coupling() == "plain":
            return {
                "requested": self.potential.device,
                "effective": "cpu",
                "guaranteed": True,
                "note": (
                    "non-KOKKOS ML-IAP coupling: LAMMPS_MLIAP_MACE keeps the model on the CPU "
                    "(lammps_mliap_mace.py _initialize_device)"
                ),
            }
        if accelerator["gpu"]:
            return {
                "requested": self.potential.device,
                "effective": f"cuda (KOKKOS, {accelerator['kokkos_gpus']} GPU(s) per node)",
                "guaranteed": False,
                "note": (
                    "LAMMPS_MLIAP_MACE follows the device of the KOKKOS data; the command line "
                    f"({' '.join(accelerator['kokkos_args'])}) requests the GPU and the build "
                    "probe requires a GPU KOKKOS back end, but the device is only fixed when "
                    "LAMMPS starts"
                ),
            }
        info = inspect_export(self.exported_model_path(), self.implementation)
        return {
            "requested": self.potential.device,
            "effective": "cpu (KOKKOS host execution space)",
            "guaranteed": False,
            "note": (
                "a host-only KOKKOS build keeps the data on the CPU; MACE accepts it only with "
                "MACE_ALLOW_CPU stored in the export ("
                + (
                    f"export allow_cpu={info.get('allow_cpu')}, force_cpu={info.get('force_cpu')}"
                    if info["readable"]
                    else f"not verified: {info['note']}"
                )
                + ")"
            ),
        }


class MaceLammpsPairStyleBridge(_MaceLammpsBridge):
    implementation = "pair-mace"

    def single_rank(self) -> bool:
        return not self.engine.mpi_launcher

    def pair_style(self) -> str:
        # pair_mace.cpp: no_domain_decomposition hands the model every atom of
        # one domain ("TODO: add check against MPI rank") -- single rank only.
        return "mace no_domain_decomposition" if self.single_rank() else "mace"

    def pair_coeff(self) -> tuple[str, ...]:
        return (f"* * {self.staged_model_name()} {' '.join(self.potential.elements)}",)

    def launch(self) -> lammps_engine.LammpsLaunch:
        """``pair mace`` in the build; GPUs hidden for a CPU request.

        pair_mace.cpp chooses CUDA whenever ``torch::cuda::is_available()``,
        whatever was requested, so ``device = "cpu"`` runs LAMMPS with
        ``CUDA_VISIBLE_DEVICES=""``. A CUDA ordinal other than 0 is refused:
        pair_mace picks the GPU by MPI local rank.
        """
        device = self.potential.device
        if device.startswith("mps"):
            raise ConfigError("potential.device = 'mps' has no LAMMPS equivalent")
        ordinal = _cuda_ordinal(device)
        if ordinal not in (None, 0):
            raise ConfigError(
                f"potential.device = {device!r}: pair_mace selects the GPU by MPI local rank "
                "(cuda:<local rank>); select GPUs with CUDA_VISIBLE_DEVICES in the job "
                "environment and request device = 'cuda'"
            )
        env = {} if ordinal is not None else {"CUDA_VISIBLE_DEVICES": ""}
        return lammps_engine.resolve_launch(
            self.engine,
            required=[("pair", "mace")],
            env=env,
            notes=(
                "pair_mace selects CUDA whenever LibTorch sees a GPU (pair_mace.cpp "
                "coeff); a CPU request hides the GPUs with CUDA_VISIBLE_DEVICES=''",
            ),
        )

    def gpu_active(self) -> bool:
        self.launch()
        return _cuda_ordinal(self.potential.device) is not None

    def device_plan(self, launch) -> dict:
        cuda = _cuda_ordinal(self.potential.device) is not None
        executable = launch.route == "executable"
        if not cuda:
            return {
                "requested": self.potential.device,
                "effective": "cpu",
                "guaranteed": executable,
                "note": (
                    "LAMMPS runs with CUDA_VISIBLE_DEVICES='', so pair_mace finds no GPU"
                    + (
                        ""
                        if executable
                        else "; on the python route the variable is set in this process, "
                        "which hides GPUs only if CUDA was not initialised here before"
                    )
                ),
            }
        return {
            "requested": self.potential.device,
            "effective": "cuda if LibTorch sees a GPU, else cpu (pair_mace.cpp coeff)",
            "guaranteed": False,
            "note": (
                "the device pair_mace reports on stdout is checked after an executable-route "
                "run and a CPU fallback is an error; the python route cannot read it"
            ),
        }

    def warnings(self) -> tuple[str, ...]:
        return (PAIR_MACE_LEGACY_WARNING,)

    def observed_device(self, native: dict, workdir: Path) -> str | None:
        """pair_mace's device line from stdout (executable route); refuse a mismatch."""
        if native.get("route") != "executable":
            return None
        try:
            text = (Path(workdir) / lammps_engine.STDOUT_FILE).read_text(
                encoding="utf-8", errors="replace"
            )
        except OSError:
            return None
        observed = None
        if "setting device type to torch::kCUDA" in text:
            observed = "cuda"
        elif "setting device type to torch::kCPU" in text:
            observed = "cpu"
        requested = "cuda" if _cuda_ordinal(self.potential.device) is not None else "cpu"
        if observed is not None and observed != requested:
            raise MlipError(
                f"pair_mace ran on {observed} but potential.device is "
                f"{self.potential.device!r} (pair_mace chooses CUDA whenever LibTorch sees a "
                "GPU)"
            )
        return observed


__all__ = [
    "MaceLammpsMliapBridge",
    "MaceLammpsPairStyleBridge",
    "EXPORT_SUFFIX",
    "EXPORT_FORMAT",
    "EXPORTED_MODEL_OPTION",
    "LEGACY_EXPORTED_MODEL_OPTION",
    "EXPORTED_SHA256_OPTION",
    "MLIAP_COUPLING_OPTION",
    "PAIR_MACE_LEGACY_WARNING",
    "default_export_path",
    "inspect_export",
]
