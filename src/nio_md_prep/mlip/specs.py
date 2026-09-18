"""Declarative specifications: what to run, with what, under what conventions.

The central architectural rule of this subsystem lives here: **the potential
type and the dynamics engine are independent concepts**. A
:class:`PotentialSpec` says what the model is and nothing about how it is
integrated; an :class:`EngineSpec` says which dynamics engine drives it and
nothing about which model it is; a :class:`SimulationSpec` says what physics is
wanted and nothing about either. Joining them is the registry's job, not a
field on any of these objects.

Specs are pure data. Constructing one imports nothing heavier than
``pathlib``: no torch, no OpenMM, no LAMMPS. That is what lets
``nio-md-prep mlip validate`` check a configuration on a laptop for a job that
will only ever run on a GPU node.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

from .capabilities import RequirementSet
from .errors import ConfigError
from .units import ENERGY_CONVENTIONS, TOTAL, lammps_unit_system

PRECISIONS = ("float32", "float64")
DEVICES = ("cpu", "cuda", "mps")
ENGINE_KINDS = ("ase", "lammps", "openmm")
#: ``optimize`` is expressible and fully negotiated (it requires energy and
#: forces), but no CLI command runs one in this release; the classical
#: relaxation workflows remain the way production geometries are relaxed.
TASKS = ("singlepoint", "optimize", "md")
ENSEMBLES = ("nve", "nvt", "npt")

#: Only full-system MLIP is implemented. ``selection`` is reserved so that a
#: later ML/MM or fixed-ML-region capability can be expressed in the same
#: schema without a breaking change; it is rejected today rather than being
#: silently treated as ``all``.
REGIONS = ("all", "selection")
IMPLEMENTED_REGIONS = ("all",)

THERMOSTATS = ("langevin", "nose-hoover", "berendsen", "csvr")
BAROSTATS = ("mtk", "berendsen", "parrinello-rahman")


def _one_of(value, allowed: Sequence[str], field_name: str):
    if value not in allowed:
        raise ConfigError(
            f"{field_name} must be one of {', '.join(repr(a) for a in allowed)}; got {value!r}"
        )
    return value


def _positive(value, field_name: str):
    if value is None:
        return None
    if not isinstance(value, (int, float)) or value <= 0:
        raise ConfigError(f"{field_name} must be a positive number; got {value!r}")
    return value


def _device(value: str) -> str:
    """Accept ``cpu``, ``cuda``, ``cuda:1``, ``mps``."""
    head = value.split(":", 1)[0]
    if head not in DEVICES:
        raise ConfigError(
            f"potential.device must start with one of {', '.join(DEVICES)}; got {value!r}"
        )
    if ":" in value and not value.split(":", 1)[1].isdigit():
        raise ConfigError(f"potential.device ordinal must be an integer; got {value!r}")
    return value


def _elements(values: Iterable[str], field_name: str) -> tuple[str, ...]:
    symbols = tuple(dict.fromkeys(values))
    if not symbols:
        raise ConfigError(f"{field_name} must list at least one chemical element")
    for symbol in symbols:
        if not symbol or not symbol[0].isupper() or not symbol.isalpha():
            raise ConfigError(f"{field_name} contains an invalid element symbol: {symbol!r}")
    return symbols


def _sha256(value: str | None, field_name: str) -> str | None:
    if value is None:
        return None
    text = value.strip().lower()
    if len(text) != 64 or any(c not in "0123456789abcdef" for c in text):
        raise ConfigError(f"{field_name} must be a 64-character hex SHA256; got {value!r}")
    return text


# ---------------------------------------------------------------------------
# Structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StructureSpec:
    """Where the geometry comes from, and how to read it.

    The MLIP layer's interchange structure is an ASE ``Atoms`` object: for a
    full-system MLIP, symbols, coordinates, cell and PBC are all a potential
    needs, and ASE is also the most direct route into MACE. The repository's
    topology-rich LAMMPS representation is untouched; ``format =
    "lammps-data-nio"`` converts one into ``Atoms`` without the reverse being
    implied.
    """

    path: Path
    format: str = "auto"
    index: int = -1
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        if self.label is None:
            object.__setattr__(self, "label", self.path.name)


# ---------------------------------------------------------------------------
# Potentials
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PotentialSpec:
    """Base class: everything a potential declares independently of any engine."""

    kind: ClassVar[str] = "abstract"

    label: str = "potential"
    energy_convention: str = TOTAL
    atomic_reference_energies: Mapping[str, float] | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        _one_of(self.energy_convention, ENERGY_CONVENTIONS, "potential.energy_convention")
        if self.atomic_reference_energies is not None:
            object.__setattr__(
                self,
                "atomic_reference_energies",
                {str(k): float(v) for k, v in dict(self.atomic_reference_energies).items()},
            )

    @property
    def elements(self) -> tuple[str, ...]:
        raise NotImplementedError

    def model_files(self) -> tuple[Path, ...]:
        """Files whose SHA256 belongs in the manifest."""
        return ()

    def declared_hashes(self) -> dict[str, str]:
        """Declared ``path -> sha256`` pairs, verified before execution."""
        return {}

    def as_dict(self) -> dict:
        return {
            "kind": self.kind,
            "label": self.label,
            "elements": list(self.elements),
            "energy_convention": self.energy_convention,
            "atomic_reference_energies": (
                dict(self.atomic_reference_energies)
                if self.atomic_reference_energies
                else None
            ),
            "notes": self.notes,
        }


@dataclass(frozen=True)
class MacePotentialSpec(PotentialSpec):
    """A trained MACE model, described independently of ASE/LAMMPS/OpenMM.

    ``implementation`` is a *preference order*, not a hard selection: for
    LAMMPS, ``("mliap", "pair_mace")`` asks for the ML-IAP route first --
    the newer interface, which brings GPU acceleration, multi-GPU inference
    and atomic virials -- and falls back to a plain MACE pair style only if
    ML-IAP is unavailable. An empty tuple means "let the registry choose".

    ``cutoff_angstrom`` is optional because it is usually discoverable from
    the model file itself; when it is declared and the file disagrees,
    :mod:`nio_md_prep.mlip.potentials.mace` reports the conflict instead of
    trusting either silently.
    """

    kind: ClassVar[str] = "mace"

    model_path: Path = Path("model.pt")
    model_sha256: str | None = None
    declared_elements: tuple[str, ...] = ()
    device: str = "cpu"
    precision: str = "float64"
    cutoff_angstrom: float | None = None
    implementation: tuple[str, ...] = ()
    model_format: str = "mace-torch"
    compile_mode: str | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        object.__setattr__(self, "model_path", Path(self.model_path))
        object.__setattr__(self, "model_sha256", _sha256(self.model_sha256, "potential.sha256"))
        object.__setattr__(
            self, "declared_elements", _elements(self.declared_elements, "potential.elements")
        )
        object.__setattr__(self, "device", _device(self.device))
        _one_of(self.precision, PRECISIONS, "potential.precision")
        _positive(self.cutoff_angstrom, "potential.cutoff_angstrom")
        object.__setattr__(self, "implementation", tuple(self.implementation))

    @property
    def elements(self) -> tuple[str, ...]:
        return self.declared_elements

    @property
    def wants_gpu(self) -> bool:
        return self.device.startswith(("cuda", "mps"))

    def model_files(self) -> tuple[Path, ...]:
        return (self.model_path,)

    def declared_hashes(self) -> dict[str, str]:
        return {str(self.model_path): self.model_sha256} if self.model_sha256 else {}

    def as_dict(self) -> dict:
        return {
            **super().as_dict(),
            "model_path": str(self.model_path),
            "model_sha256": self.model_sha256,
            "device": self.device,
            "precision": self.precision,
            "cutoff_angstrom": self.cutoff_angstrom,
            "implementation": list(self.implementation),
            "model_format": self.model_format,
            "compile_mode": self.compile_mode,
        }


@dataclass(frozen=True)
class LammpsMlipPotentialSpec(PotentialSpec):
    """A LAMMPS-native machine-learned potential, described by its pair style.

    Deliberately *not* tied to one framework. ``pair_style`` is the primary
    datum, so ``deepmd``, ``mliap``, ``pace``, a MACE pair style, or a pair
    style that does not exist yet are all expressible without a code change;
    ``framework`` is a label for provenance and package checking, not a
    discriminator that the rest of the subsystem branches on.

    ``type_map`` is the LAMMPS-type-to-element mapping. It is what makes
    element-coverage validation possible for a model whose element list is
    otherwise buried in a binary file, and it is what the ASE bridge needs to
    map an ``Atoms`` object onto LAMMPS types.
    """

    kind: ClassVar[str] = "lammps"

    pair_style: str = ""
    pair_coeff: tuple[str, ...] = ()
    type_map: Mapping[int, str] = field(default_factory=dict)
    required_packages: tuple[str, ...] = ()
    units: str = "metal"
    atom_style: str = "atomic"
    model_paths: tuple[Path, ...] = ()
    model_hashes: Mapping[str, str] = field(default_factory=dict)
    framework: str | None = None
    extra_commands: tuple[str, ...] = ()
    newton: str | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if not self.pair_style.strip():
            raise ConfigError("potential.pair_style is required for a LAMMPS-native MLIP")
        if not self.pair_coeff:
            raise ConfigError("potential.pair_coeff must list at least one pair_coeff line")
        object.__setattr__(self, "pair_style", self.pair_style.strip())
        object.__setattr__(self, "pair_coeff", tuple(str(c).strip() for c in self.pair_coeff))
        type_map = {int(k): str(v) for k, v in dict(self.type_map).items()}
        if not type_map:
            raise ConfigError(
                "potential.type_map is required: without a LAMMPS-type-to-element "
                "mapping the element coverage of this potential cannot be validated"
            )
        if sorted(type_map) != list(range(1, len(type_map) + 1)):
            raise ConfigError(
                "potential.type_map must cover LAMMPS types 1..N with no gaps; got "
                f"types {sorted(type_map)}"
            )
        _elements(type_map.values(), "potential.type_map")
        object.__setattr__(self, "type_map", type_map)
        object.__setattr__(self, "required_packages", tuple(self.required_packages))
        # Rejects unit styles with no honest eV/Angstrom mapping (e.g. lj).
        lammps_unit_system(self.units)
        object.__setattr__(self, "model_paths", tuple(Path(p) for p in self.model_paths))
        object.__setattr__(
            self,
            "model_hashes",
            {str(k): _sha256(str(v), f"potential.model_hashes[{k}]") for k, v in
             dict(self.model_hashes).items()},
        )
        object.__setattr__(self, "extra_commands", tuple(self.extra_commands))
        if self.newton is not None:
            _one_of(self.newton, ("on", "off"), "potential.newton")

    @property
    def elements(self) -> tuple[str, ...]:
        return tuple(self.type_map[t] for t in sorted(self.type_map))

    @property
    def unit_system_name(self) -> str:
        return lammps_unit_system(self.units).name

    def model_files(self) -> tuple[Path, ...]:
        return self.model_paths

    def declared_hashes(self) -> dict[str, str]:
        return dict(self.model_hashes)

    def render_pair_commands(self) -> tuple[str, ...]:
        """The exact LAMMPS commands this potential contributes, in order.

        Preserved verbatim in the manifest: reproducing a LAMMPS run means
        reproducing these strings, not re-deriving them from a config.
        """
        lines = [f"pair_style {self.pair_style}"]
        lines += [f"pair_coeff {c}" if not c.startswith("pair_coeff") else c
                  for c in self.pair_coeff]
        lines += list(self.extra_commands)
        return tuple(lines)

    def as_dict(self) -> dict:
        return {
            **super().as_dict(),
            "pair_style": self.pair_style,
            "pair_coeff": list(self.pair_coeff),
            "type_map": {str(k): v for k, v in sorted(self.type_map.items())},
            "required_packages": list(self.required_packages),
            "units": self.units,
            "unit_system": self.unit_system_name,
            "atom_style": self.atom_style,
            "model_paths": [str(p) for p in self.model_paths],
            "model_hashes": dict(self.model_hashes),
            "framework": self.framework,
            "extra_commands": list(self.extra_commands),
            "newton": self.newton,
            "rendered_commands": list(self.render_pair_commands()),
        }


@dataclass(frozen=True)
class MockPotentialSpec(PotentialSpec):
    """A deterministic analytic potential used to exercise this machinery.

    Not science. It exists so that configuration parsing, bridge resolution,
    capability negotiation, unit handling, provenance and the ``singlepoint``
    and ``smoke-md`` code paths can be tested end to end in ordinary CI, on a
    machine with no torch, no LAMMPS and no OpenMM.
    """

    kind: ClassVar[str] = "mock"

    declared_elements: tuple[str, ...] = ("H",)
    epsilon_eV: float = 0.01
    sigma_angstrom: float = 2.5
    cutoff_angstrom: float = 6.0
    precision: str = "float64"

    def __post_init__(self) -> None:
        super().__post_init__()
        object.__setattr__(
            self, "declared_elements", _elements(self.declared_elements, "potential.elements")
        )
        _positive(self.epsilon_eV, "potential.epsilon_eV")
        _positive(self.sigma_angstrom, "potential.sigma_angstrom")
        _positive(self.cutoff_angstrom, "potential.cutoff_angstrom")
        _one_of(self.precision, PRECISIONS, "potential.precision")

    @property
    def elements(self) -> tuple[str, ...]:
        return self.declared_elements

    def as_dict(self) -> dict:
        return {
            **super().as_dict(),
            "epsilon_eV": self.epsilon_eV,
            "sigma_angstrom": self.sigma_angstrom,
            "cutoff_angstrom": self.cutoff_angstrom,
            "precision": self.precision,
        }


# ---------------------------------------------------------------------------
# Engines
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EngineSpec:
    """Which dynamics engine runs the model, and how it is configured.

    Carries nothing about the potential. ``platform`` is OpenMM's notion
    (``CPU``/``CUDA``/``OpenCL``/``Reference``); ``executable`` and
    ``command`` are LAMMPS's; ``threads`` applies to all three. Unknown
    engine-specific knobs go in ``options`` and are recorded in the manifest
    verbatim.
    """

    kind: str = "ase"
    platform: str | None = None
    precision: str | None = None
    threads: int | None = None
    executable: str | None = None
    options: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _one_of(self.kind, ENGINE_KINDS, "engine.kind")
        if self.precision is not None:
            _one_of(self.precision, PRECISIONS, "engine.precision")
        if self.threads is not None and (not isinstance(self.threads, int) or self.threads < 1):
            raise ConfigError(f"engine.threads must be a positive integer; got {self.threads!r}")
        object.__setattr__(self, "options", dict(self.options))

    def as_dict(self) -> dict:
        return {
            "kind": self.kind,
            "platform": self.platform,
            "precision": self.precision,
            "threads": self.threads,
            "executable": self.executable,
            "options": dict(self.options),
        }


# ---------------------------------------------------------------------------
# Simulations
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SimulationSpec:
    """The physics asked for, stated independently of model and engine.

    :meth:`required_capabilities` is the important method: it is what lets the
    resolver reject an impossible job before a batch script exists. NPT needs
    a trustworthy stress/virial route; an analysis job asking for per-atom
    energies needs a potential that decomposes them; a geometry optimisation
    needs forces.
    """

    task: str = "singlepoint"
    ensemble: str | None = None
    temperature_K: float | None = None
    pressure_bar: float | None = None
    timestep_fs: float | None = None
    steps: int = 0
    seed: int | None = None
    thermostat: str | None = None
    barostat: str | None = None
    thermostat_damping_fs: float | None = None
    barostat_damping_fs: float | None = None
    fmax_eV_per_A: float = 0.05
    max_optimizer_steps: int = 200
    compute_stress: bool = False
    compute_per_atom_energy: bool = False
    trajectory_interval: int = 10
    log_interval: int = 10
    region: str = "all"
    selection: str | None = None
    energy_convention: str | None = None

    def __post_init__(self) -> None:
        _one_of(self.task, TASKS, "simulation.task")
        _one_of(self.region, REGIONS, "simulation.region")
        if self.region not in IMPLEMENTED_REGIONS:
            raise ConfigError(
                f"simulation.region = {self.region!r} is reserved but not implemented. "
                "This subsystem currently applies the MLIP to the whole system only; "
                "ML/MM partitioning, fixed ML regions and electrostatic embedding are "
                'out of scope, so use region = "all".'
            )
        if self.selection is not None and self.region == "all":
            raise ConfigError(
                'simulation.selection is meaningless with region = "all"; remove it'
            )
        if self.energy_convention is not None:
            _one_of(self.energy_convention, ENERGY_CONVENTIONS, "simulation.energy_convention")
        if self.task == "md":
            self._check_md()
        elif self.ensemble is not None:
            raise ConfigError(
                f"simulation.ensemble is only meaningful for task = 'md'; got task "
                f"{self.task!r} with ensemble {self.ensemble!r}"
            )
        if self.task == "optimize":
            _positive(self.fmax_eV_per_A, "simulation.fmax_eV_per_A")
            if self.max_optimizer_steps < 1:
                raise ConfigError("simulation.max_optimizer_steps must be at least 1")

    def _check_md(self) -> None:
        if self.ensemble is None:
            raise ConfigError("simulation.ensemble is required for task = 'md'")
        _one_of(self.ensemble, ENSEMBLES, "simulation.ensemble")
        _positive(self.timestep_fs, "simulation.timestep_fs")
        if self.timestep_fs is None:
            raise ConfigError("simulation.timestep_fs is required for task = 'md'")
        if not isinstance(self.steps, int) or self.steps < 1:
            raise ConfigError("simulation.steps must be a positive integer for task = 'md'")
        if self.ensemble in ("nvt", "npt"):
            if self.temperature_K is None:
                raise ConfigError(
                    f"simulation.temperature_K is required for the {self.ensemble} ensemble"
                )
            _positive(self.temperature_K, "simulation.temperature_K")
        if self.ensemble == "npt" and self.pressure_bar is None:
            raise ConfigError("simulation.pressure_bar is required for the npt ensemble")
        if self.thermostat is not None:
            _one_of(self.thermostat, THERMOSTATS, "simulation.thermostat")
        if self.barostat is not None:
            _one_of(self.barostat, BAROSTATS, "simulation.barostat")
        if self.ensemble != "npt" and self.barostat is not None:
            raise ConfigError("simulation.barostat is only meaningful for the npt ensemble")

    @property
    def needs_velocities(self) -> bool:
        return self.task == "md"

    def required_capabilities(self, *, periodic: bool = False) -> RequirementSet:
        """Translate the requested physics into capabilities that must exist."""
        reasons = {"energy": f"every {self.task} job reports an energy"}
        requirements = RequirementSet(energy=True, reasons=reasons)

        if self.task in ("optimize", "md"):
            requirements = requirements.requiring(
                "forces",
                "geometry optimisation and molecular dynamics integrate forces"
                if self.task == "md"
                else "geometry optimisation minimises against forces",
            )
        else:
            requirements = requirements.requiring(
                "forces", "single-point evaluation reports forces alongside the energy"
            )

        if self.ensemble == "npt":
            requirements = requirements.requiring(
                "stress",
                "the npt ensemble integrates the simulation cell against the "
                "stress/virial, so an absent or untrustworthy virial route is fatal",
            )
        elif self.compute_stress:
            requirements = requirements.requiring(
                "stress", "simulation.compute_stress was requested"
            )

        if self.compute_per_atom_energy:
            requirements = requirements.requiring(
                "per_atom_energy",
                "this job requests a per-atom energy decomposition, which not every "
                "potential/engine route exposes",
            )

        if periodic:
            requirements = requirements.requiring(
                "periodic", "the input structure declares periodic boundary conditions"
            )

        if self.energy_convention is not None:
            from dataclasses import replace as _replace

            requirements = _replace(requirements, energy_convention=self.energy_convention)
        return requirements

    def as_dict(self) -> dict:
        return {
            "task": self.task,
            "ensemble": self.ensemble,
            "temperature_K": self.temperature_K,
            "pressure_bar": self.pressure_bar,
            "timestep_fs": self.timestep_fs,
            "steps": self.steps,
            "seed": self.seed,
            "thermostat": self.thermostat,
            "barostat": self.barostat,
            "thermostat_damping_fs": self.thermostat_damping_fs,
            "barostat_damping_fs": self.barostat_damping_fs,
            "fmax_eV_per_A": self.fmax_eV_per_A,
            "max_optimizer_steps": self.max_optimizer_steps,
            "compute_stress": self.compute_stress,
            "compute_per_atom_energy": self.compute_per_atom_energy,
            "trajectory_interval": self.trajectory_interval,
            "log_interval": self.log_interval,
            "region": self.region,
            "selection": self.selection,
            "energy_convention": self.energy_convention,
        }


def as_singlepoint(simulation: SimulationSpec) -> SimulationSpec:
    """Derive the single-point form of any simulation spec.

    A :class:`SimulationSpec` validates itself, so the ensemble and integrator
    settings must be dropped rather than carried along -- ``ensemble`` is
    meaningless for ``task = "singlepoint"``. The stress request survives,
    because the endpoints of an NPT run are exactly where a missing virial
    should still be visible.
    """
    if simulation.task == "singlepoint":
        return simulation
    from dataclasses import replace

    return replace(
        simulation,
        task="singlepoint",
        ensemble=None,
        steps=0,
        timestep_fs=None,
        thermostat=None,
        barostat=None,
        compute_stress=simulation.compute_stress or simulation.ensemble == "npt",
    )


@dataclass(frozen=True)
class JobSpec:
    """One validated MLIP job: structure + potential + engine + simulation."""

    potential: PotentialSpec
    engine: EngineSpec
    simulation: SimulationSpec
    structure: StructureSpec | None = None
    name: str = "mlip-job"
    source: Path | None = None

    def __post_init__(self) -> None:
        if self.source is not None:
            object.__setattr__(self, "source", Path(self.source))

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "source": str(self.source) if self.source else None,
            "potential": self.potential.as_dict(),
            "engine": self.engine.as_dict(),
            "simulation": self.simulation.as_dict(),
            "structure": (
                {
                    "path": str(self.structure.path),
                    "format": self.structure.format,
                    "index": self.structure.index,
                    "label": self.structure.label,
                }
                if self.structure
                else None
            ),
        }


#: ``potential.kind`` -> spec class. Adding a potential family means adding a
#: spec class here and registering its bridges; nothing else branches on kind.
POTENTIAL_SPECS: dict[str, type[PotentialSpec]] = {
    MacePotentialSpec.kind: MacePotentialSpec,
    LammpsMlipPotentialSpec.kind: LammpsMlipPotentialSpec,
    MockPotentialSpec.kind: MockPotentialSpec,
}


__all__ = [
    "PRECISIONS",
    "DEVICES",
    "ENGINE_KINDS",
    "TASKS",
    "ENSEMBLES",
    "REGIONS",
    "IMPLEMENTED_REGIONS",
    "THERMOSTATS",
    "BAROSTATS",
    "StructureSpec",
    "PotentialSpec",
    "MacePotentialSpec",
    "LammpsMlipPotentialSpec",
    "MockPotentialSpec",
    "EngineSpec",
    "SimulationSpec",
    "as_singlepoint",
    "JobSpec",
    "POTENTIAL_SPECS",
]
