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
#: How the barostat couples the cell degrees of freedom. ``in-plane`` barostats
#: only the two in-plane directions of a slab (LAMMPS ``x P P Pd y P P Pd
#: couple xy``); which engine supports which coupling is the engine's decision.
BAROSTAT_COUPLINGS = ("isotropic", "anisotropic", "in-plane")

#: The one default thermostat damping time, used by *every* engine when
#: ``simulation.thermostat_damping_fs`` is unset: Langevin friction 1/tau (ASE,
#: OpenMM), ``Tdamp`` (LAMMPS), ``taut``/``tdamp`` (ASE Berendsen/NHC/Bussi).
#: One constant, so the same configuration never means a 10x different
#: thermostat on two engines.
DEFAULT_THERMOSTAT_DAMPING_FS = 100.0
#: The one default barostat damping time (``Pdamp``/``pfactor`` time scale).
DEFAULT_BAROSTAT_DAMPING_FS = 1000.0
#: A fully periodic cell whose atoms leave an empty slab wider than this along
#: a lattice direction is treated as a vacuum slab, not a bulk crystal: an
#: isotropic/anisotropic barostat would scale the vacuum with the solid.
DEFAULT_VACUUM_GAP_THRESHOLD_ANGSTROM = 5.0

#: ``engine.runtime``: how an engine with more than one execution route is
#: driven. Only LAMMPS has more than one; see :class:`EngineSpec`.
RUNTIMES = ("auto", "python", "executable")
#: OpenMM's platform ``Precision`` property values.
PLATFORM_PRECISIONS = ("single", "mixed", "double")
#: OpenMM platforms that have a ``Precision`` property at all. The CPU platform
#: has only ``Threads``/``DeterministicForces``; Reference has none.
OPENMM_GPU_PLATFORMS = ("CUDA", "OpenCL", "HIP")
#: Largest seed accepted: OpenMM and LAMMPS both take a signed 32-bit integer.
MAX_SEED = 2**31 - 1

#: SimulationSpec fields that only an MD task reads; refused for other tasks.
_MD_ONLY_FIELDS = (
    "thermostat",
    "barostat",
    "thermostat_damping_fs",
    "barostat_damping_fs",
    "barostat_coupling",
)


def _one_of(value, allowed: Sequence[str], field_name: str):
    if value not in allowed:
        raise ConfigError(
            f"{field_name} must be one of {', '.join(repr(a) for a in allowed)}; got {value!r}"
        )
    return value


def _positive(value, field_name: str):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
        raise ConfigError(f"{field_name} must be a positive number; got {value!r}")
    return value


def _positive_int(value, field_name: str):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigError(f"{field_name} must be a positive integer; got {value!r}")
    return value


def _flag(value, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigError(f"{field_name} must be true or false; got {value!r}")
    return value


def _arguments(value, field_name: str) -> tuple[str, ...]:
    """A command-line argument vector: a list of non-empty strings, never one string.

    A single string would have to be split, and no splitting rule is right for
    both ``mpirun -np 4`` and ``C:/Program Files/MPI/mpiexec.exe``.
    """
    if isinstance(value, str):
        raise ConfigError(
            f"{field_name} must be a list of arguments (e.g. [\"mpirun\", \"-np\", \"4\"]), "
            f"not a single string; got {value!r}"
        )
    values = tuple(value)
    for item in values:
        if not isinstance(item, str) or not item:
            raise ConfigError(f"{field_name} entries must be non-empty strings; got {item!r}")
    return values


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


def path_key(path) -> str:
    """The comparison key for a file path: absolute, normalised, case-folded
    where the filesystem is case-insensitive.

    Declared hashes and hashed files are matched on this key, so ``nio.pb``
    and ``/abs/dir/nio.pb`` name the same file whenever they resolve to it.
    """
    import os

    return os.path.normcase(os.path.abspath(os.fspath(path)))


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

    ``pbc`` is an explicit per-axis periodicity, e.g. ``(True, True, False)``
    for a slab. It is *required* for ``lammps-data-nio`` (a LAMMPS data file
    does not say which box faces are periodic; the deck's ``boundary`` command
    does), and for any other format it replaces the periodicity the file
    declared -- the way to say that a POSCAR slab, which ASE always reads as
    fully periodic, is periodic in-plane only. The value used is recorded in
    the manifest's structure description.
    """

    path: Path
    format: str = "auto"
    index: int = -1
    label: str | None = None
    pbc: tuple[bool, bool, bool] | None = None

    def __post_init__(self) -> None:
        from .structures import FORMATS, FORMATS_WITHOUT_PBC

        object.__setattr__(self, "path", Path(self.path))
        if self.label is None:
            object.__setattr__(self, "label", self.path.name)
        _one_of(self.format, FORMATS, "structure.format")
        if isinstance(self.index, bool) or not isinstance(self.index, int):
            raise ConfigError(f"structure.index must be an integer; got {self.index!r}")
        if self.pbc is not None:
            object.__setattr__(self, "pbc", _pbc(self.pbc, "structure.pbc"))
        elif self.format in FORMATS_WITHOUT_PBC:
            raise ConfigError(
                f"structure.format = {self.format!r} needs an explicit structure.pbc: a "
                "LAMMPS data file does not record which box faces are periodic (the "
                "input deck's boundary command does). The classical slab builds use "
                "'boundary p p f', i.e. pbc = [true, true, false]."
            )


def _pbc(value, field_name: str) -> tuple[bool, bool, bool]:
    values = tuple(value) if not isinstance(value, (str, bool)) else ()
    if len(values) != 3 or not all(isinstance(v, bool) for v in values):
        raise ConfigError(
            f"{field_name} must list three booleans, one per cell axis "
            f"(e.g. [true, true, false] for a slab); got {value!r}"
        )
    return values  # type: ignore[return-value]


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

    ``implementation`` is a *preference order* over registered implementation
    names: for LAMMPS, ``("mliap", "pair-mace")`` selects the ML-IAP route --
    the newer interface, which brings GPU acceleration, multi-GPU inference
    and atomic virials -- whenever the registry has it for the engine, and
    the MACE pair style only when it does not. Availability on this machine
    is *not* consulted: a preferred route that cannot run fails with
    :class:`~nio_md_prep.mlip.errors.MissingDependencyError` rather than
    silently switching to another route. A name that no engine registers for
    MACE is refused by the registry. An empty tuple means "let the registry
    choose".

    ``head`` selects one head of a multi-head MACE model (MACE's own ``head``
    keyword). ``None`` leaves the choice to the model, which only works for a
    single-head model or one with a head named ``default``.

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
    head: str | None = None

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
        implementation = (
            (self.implementation,)
            if isinstance(self.implementation, str)
            else tuple(self.implementation)
        )
        for name in implementation:
            if not isinstance(name, str) or not name:
                raise ConfigError(
                    f"potential.implementation entries must be non-empty strings; got {name!r}"
                )
        object.__setattr__(self, "implementation", implementation)
        if self.head is not None and (not isinstance(self.head, str) or not self.head.strip()):
            raise ConfigError(f"potential.head must be a non-empty string; got {self.head!r}")

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
            "head": self.head,
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
        shared = {
            symbol: sorted(t for t, s in type_map.items() if s == symbol)
            for symbol in dict.fromkeys(type_map.values())
        }
        shared = {symbol: types for symbol, types in shared.items() if len(types) > 1}
        if shared:
            listing = "; ".join(f"{symbol}: types {types}" for symbol, types in shared.items())
            raise ConfigError(
                f"potential.type_map maps more than one LAMMPS type to one element ({listing}). "
                "An ASE structure carries elements, not LAMMPS types, so nothing says "
                "which atom gets which type: the native route would put every such atom "
                "on the first type and the ASE route on the last. Distinct magnetic "
                "sublattices (e.g. spin-up and spin-down Ni in AFM NiO) cannot be "
                "expressed through type_map yet; use exactly one type per element."
            )
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
        listed = {path_key(p) for p in self.model_paths}
        orphans = sorted(k for k in self.model_hashes if path_key(k) not in listed)
        if orphans:
            raise ConfigError(
                f"potential.model_hashes names file(s) {', '.join(orphans)} that are not "
                "in potential.model_files, so they would never be hashed; list each "
                "hashed file in model_files (both are resolved against the "
                "configuration's directory)"
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


#: OpenMM platform names this subsystem knows. A typo is refused at validate
#: time instead of surfacing as a raw ``OpenMMException`` at run time; whether
#: the named platform is actually installed is checked by the engine.
OPENMM_PLATFORMS = ("Reference", "CPU", "CUDA", "OpenCL", "HIP")

#: EngineSpec fields that only one engine interprets. Setting one for another
#: engine is refused rather than silently ignored.
_LAMMPS_ONLY_FIELDS = ("executable", "mpi_launcher", "lammps_args", "timeout_s")
_OPENMM_ONLY_FIELDS = ("platform", "platform_precision")


@dataclass(frozen=True)
class EngineSpec:
    """Which dynamics engine runs the model, and how it is configured.

    Carries nothing about the potential. Fields that only one engine
    interprets are refused for the others, never silently dropped:

    ``runtime`` (LAMMPS)
        ``"python"`` drives the LAMMPS python module in-process,
        ``"executable"`` runs ``executable`` (default ``lmp``) as a
        subprocess, and ``"auto"`` lets the engine choose -- but an explicitly
        given ``executable``, ``mpi_launcher`` or ``timeout_s`` is only
        honoured by the executable route, so ``auto`` must choose it then.
        ``runtime = "python"`` together with any of those three is refused.
        Other engines have one route and accept only ``"auto"``.
    ``mpi_launcher`` (LAMMPS)
        Argument vector prepended to the executable, e.g.
        ``["mpirun", "-np", "4"]``. A list, never one string.
    ``lammps_args`` (LAMMPS)
        Extra LAMMPS command-line arguments (e.g. the KOKKOS
        ``["-k", "on", "g", "1", "-sf", "kk"]``), passed as ``cmdargs`` on
        the python route and appended on the executable route.
    ``timeout_s`` (LAMMPS executable route)
        Wall-clock limit for the subprocess.
    ``platform`` / ``platform_precision`` (OpenMM)
        ``platform`` is one of :data:`OPENMM_PLATFORMS`.
        ``platform_precision`` is the platform's ``Precision`` property
        (``single``/``mixed``/``double``), which only CUDA, OpenCL and HIP
        have, so it requires one of those platforms to be named explicitly.
        It controls force accumulation and integration, not the dtype the
        MACE model is evaluated in (that is ``potential.precision``).
    ``threads`` (every engine; each engine applies it or refuses it)
        ASE with MACE: torch intra-op threads. LAMMPS: OpenMP threads per
        MPI rank. OpenMM: the CPU platform's ``Threads`` property (refused on
        other platforms). An engine or route that cannot honour it must raise
        from :meth:`~nio_md_prep.mlip.bridges.base.Bridge.check_simulation`.
    ``precision``
        A cross-check, not a second knob: when set it must be a precision the
        route offers (for MACE, ``potential.precision``), which capability
        negotiation enforces. ``potential.precision`` is what is used.

    Unknown engine-specific knobs go in ``options`` and are recorded in the
    manifest verbatim.
    """

    kind: str = "ase"
    platform: str | None = None
    precision: str | None = None
    threads: int | None = None
    executable: str | None = None
    runtime: str = "auto"
    mpi_launcher: tuple[str, ...] = ()
    lammps_args: tuple[str, ...] = ()
    platform_precision: str | None = None
    timeout_s: float | None = None
    options: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _one_of(self.kind, ENGINE_KINDS, "engine.kind")
        if self.precision is not None:
            _one_of(self.precision, PRECISIONS, "engine.precision")
        if self.threads is not None:
            _positive_int(self.threads, "engine.threads")
        _one_of(self.runtime, RUNTIMES, "engine.runtime")
        object.__setattr__(
            self, "mpi_launcher", _arguments(self.mpi_launcher, "engine.mpi_launcher")
        )
        object.__setattr__(
            self, "lammps_args", _arguments(self.lammps_args, "engine.lammps_args")
        )
        _positive(self.timeout_s, "engine.timeout_s")
        if self.executable is not None and (
            not isinstance(self.executable, str) or not self.executable.strip()
        ):
            raise ConfigError(
                f"engine.executable must be a non-empty string; got {self.executable!r}"
            )
        if self.platform_precision is not None:
            _one_of(self.platform_precision, PLATFORM_PRECISIONS, "engine.platform_precision")
        if self.platform is not None:
            _one_of(self.platform, OPENMM_PLATFORMS, "engine.platform")
        object.__setattr__(self, "options", dict(self.options))
        self._check_engine_specific_fields()

    def _check_engine_specific_fields(self) -> None:
        if self.kind != "lammps":
            given = [name for name in _LAMMPS_ONLY_FIELDS if getattr(self, name)]
            if self.runtime != "auto":
                given.insert(0, "runtime")
            if given:
                raise ConfigError(
                    f"engine.{', engine.'.join(given)} only apply to the LAMMPS engine; "
                    f"the {self.kind} engine has a single in-process route and would "
                    "ignore them"
                )
        elif self.runtime == "python":
            conflicting = [
                name for name in ("executable", "mpi_launcher", "timeout_s") if getattr(self, name)
            ]
            if conflicting:
                raise ConfigError(
                    f"engine.runtime = 'python' drives the LAMMPS python module in this "
                    f"process, so engine.{', engine.'.join(conflicting)} would be ignored; "
                    "use runtime = 'executable' (or 'auto') to run a LAMMPS executable"
                )
        if self.kind != "openmm":
            given = [name for name in _OPENMM_ONLY_FIELDS if getattr(self, name) is not None]
            if given:
                raise ConfigError(
                    f"engine.{', engine.'.join(given)} only apply to the OpenMM engine; "
                    f"the {self.kind} engine would ignore them"
                )
        elif self.platform_precision is not None and self.platform not in OPENMM_GPU_PLATFORMS:
            raise ConfigError(
                "engine.platform_precision sets OpenMM's platform 'Precision' property, "
                f"which only the {', '.join(OPENMM_GPU_PLATFORMS)} platforms have; name one "
                f"of them as engine.platform (got platform {self.platform!r})"
            )

    def retargeted(self, kind: str) -> "EngineSpec":
        """This spec moved to another engine kind, keeping what that engine reads.

        Used by the cross-engine comparison, which evaluates one configuration
        on several engines: an OpenMM ``platform`` means nothing to ASE, and
        carrying it over would be refused by validation.
        """
        from dataclasses import replace

        if kind == self.kind:
            return self
        reset: dict = {}
        if kind != "lammps":
            reset.update(
                runtime="auto", executable=None, mpi_launcher=(), lammps_args=(), timeout_s=None
            )
        if kind != "openmm":
            reset.update(platform=None, platform_precision=None)
        return replace(self, kind=kind, **reset)

    def as_dict(self) -> dict:
        return {
            "kind": self.kind,
            "platform": self.platform,
            "precision": self.precision,
            "threads": self.threads,
            "executable": self.executable,
            "runtime": self.runtime,
            "mpi_launcher": list(self.mpi_launcher),
            "lammps_args": list(self.lammps_args),
            "platform_precision": self.platform_precision,
            "timeout_s": self.timeout_s,
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

    Thermostat and barostat fields are only accepted where they mean
    something: a thermostat (or its damping time) is refused for NVE and for
    non-MD tasks, a barostat, its damping time and ``barostat_coupling`` for
    anything but NPT. Unset damping times resolve to
    :data:`DEFAULT_THERMOSTAT_DAMPING_FS` / :data:`DEFAULT_BAROSTAT_DAMPING_FS`
    on every engine (see :attr:`resolved_thermostat_damping_fs`).

    ``barostat_coupling`` says which cell degrees of freedom the barostat
    moves: ``isotropic`` (one scale factor), ``anisotropic`` (the three cell
    lengths independently) or ``in-plane`` (only the two periodic in-plane
    directions of a slab). ``None`` lets the engine pick its default, which
    must be recorded. Which engine supports which coupling is the engine's
    decision (:meth:`~nio_md_prep.mlip.bridges.base.Bridge.check_simulation`).

    ``vacuum_gap_threshold_angstrom`` (default
    :data:`DEFAULT_VACUUM_GAP_THRESHOLD_ANGSTROM`) is the widest empty slab
    along a lattice direction that a fully periodic cell may contain and
    still be treated as bulk for NPT.

    ``seed`` must lie in ``1..MAX_SEED``: LAMMPS requires a positive seed and
    OpenMM treats 0 as "choose one at random", so 0 cannot mean the same
    thing everywhere. ``None`` means "draw one"; engines draw it with
    :func:`resolve_seed` and record the value they used.
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
    barostat_coupling: str | None = None
    vacuum_gap_threshold_angstrom: float | None = None
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
        _flag(self.compute_stress, "simulation.compute_stress")
        _flag(self.compute_per_atom_energy, "simulation.compute_per_atom_energy")
        _positive_int(self.trajectory_interval, "simulation.trajectory_interval")
        _positive_int(self.log_interval, "simulation.log_interval")
        _positive(self.vacuum_gap_threshold_angstrom, "simulation.vacuum_gap_threshold_angstrom")
        if self.seed is not None:
            if isinstance(self.seed, bool) or not isinstance(self.seed, int) or not (
                1 <= self.seed <= MAX_SEED
            ):
                raise ConfigError(
                    f"simulation.seed must be an integer in 1..{MAX_SEED} (LAMMPS requires a "
                    "positive seed and OpenMM reads 0 as 'pick one at random'); "
                    f"got {self.seed!r}"
                )
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
        else:
            if self.ensemble is not None:
                raise ConfigError(
                    f"simulation.ensemble is only meaningful for task = 'md'; got task "
                    f"{self.task!r} with ensemble {self.ensemble!r}"
                )
            given = [name for name in _MD_ONLY_FIELDS if getattr(self, name) is not None]
            if given:
                raise ConfigError(
                    f"simulation.{', simulation.'.join(given)} only apply to task = 'md'; "
                    f"this is a {self.task!r} job, which would ignore them"
                )
        if self.task == "optimize":
            _positive(self.fmax_eV_per_A, "simulation.fmax_eV_per_A")
            _positive_int(self.max_optimizer_steps, "simulation.max_optimizer_steps")

    def _check_md(self) -> None:
        if self.ensemble is None:
            raise ConfigError("simulation.ensemble is required for task = 'md'")
        _one_of(self.ensemble, ENSEMBLES, "simulation.ensemble")
        _positive(self.timestep_fs, "simulation.timestep_fs")
        if self.timestep_fs is None:
            raise ConfigError("simulation.timestep_fs is required for task = 'md'")
        if isinstance(self.steps, bool) or not isinstance(self.steps, int) or self.steps < 1:
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
        if self.barostat_coupling is not None:
            _one_of(self.barostat_coupling, BAROSTAT_COUPLINGS, "simulation.barostat_coupling")
        _positive(self.thermostat_damping_fs, "simulation.thermostat_damping_fs")
        _positive(self.barostat_damping_fs, "simulation.barostat_damping_fs")
        if self.ensemble == "nve":
            given = [
                name
                for name in ("thermostat", "thermostat_damping_fs")
                if getattr(self, name) is not None
            ]
            if given:
                raise ConfigError(
                    f"simulation.{', simulation.'.join(given)} cannot apply to the nve "
                    "ensemble: NVE integrates without a thermostat, so the setting would "
                    "be silently ignored. Remove it, or use ensemble = 'nvt'."
                )
        if self.ensemble != "npt":
            given = [
                name
                for name in ("barostat", "barostat_damping_fs", "barostat_coupling")
                if getattr(self, name) is not None
            ]
            if given:
                raise ConfigError(
                    f"simulation.{', simulation.'.join(given)} "
                    f"{'is' if len(given) == 1 else 'are'} only meaningful for the npt "
                    "ensemble"
                )

    @property
    def needs_velocities(self) -> bool:
        return self.task == "md"

    @property
    def resolved_thermostat_damping_fs(self) -> float:
        """The thermostat damping time every engine uses (fs)."""
        return float(self.thermostat_damping_fs or DEFAULT_THERMOSTAT_DAMPING_FS)

    @property
    def resolved_barostat_damping_fs(self) -> float:
        """The barostat damping time every engine uses (fs)."""
        return float(self.barostat_damping_fs or DEFAULT_BAROSTAT_DAMPING_FS)

    @property
    def resolved_vacuum_gap_threshold_angstrom(self) -> float:
        return float(
            self.vacuum_gap_threshold_angstrom or DEFAULT_VACUUM_GAP_THRESHOLD_ANGSTROM
        )

    def required_capabilities(
        self,
        *,
        periodic: bool = False,
        pbc: Sequence[bool] | None = None,
        fixed_atoms: int = 0,
        vacuum_gaps: Mapping[int, float] | None = None,
    ) -> RequirementSet:
        """Translate the requested physics into capabilities that must exist.

        Geometry facts come from the structure (see
        :meth:`~nio_md_prep.mlip.bridges.base.Bridge.requirements`, which
        computes them with :mod:`nio_md_prep.mlip.structures`):

        ``pbc``
            Per-axis periodicity. Any periodic axis requires the ``periodic``
            capability; a mix of periodic and non-periodic axes also requires
            ``partial_periodic``. The tuple is carried as
            :attr:`RequirementSet.periodic_axes`. ``periodic=True`` without
            ``pbc`` is the older spelling of ``pbc=(True, True, True)``.
        ``fixed_atoms``
            Number of atoms frozen by an ASE ``FixAtoms`` constraint. MD or
            optimisation with frozen atoms requires the ``fixed_atoms``
            capability; NPT with frozen atoms is refused on every engine.
        ``vacuum_gaps``
            ``axis -> width`` of empty slabs wider than the threshold, for a
            fully periodic cell. NPT with an isotropic or anisotropic
            barostat is refused when any are present.

        Stress (``compute_stress``) and NPT need a cell periodic along all
        three axes. The one exception is NPT with ``barostat_coupling =
        "in-plane"`` on a structure with at least two periodic axes, which is
        left to the engine to accept or refuse. Geometry refusals raise
        :class:`~nio_md_prep.mlip.errors.ConfigError`: they are problems with
        the request, not with any route's capabilities.
        """
        reasons = {"energy": f"every {self.task} job reports an energy"}
        requirements = RequirementSet(energy=True, reasons=reasons)
        if pbc is None and periodic:
            pbc = (True, True, True)
        axes = None if pbc is None else tuple(bool(p) for p in pbc)
        if axes is not None and len(axes) != 3:
            raise ConfigError(f"pbc must have one entry per cell axis; got {pbc!r}")
        self._check_geometry(axes, fixed_atoms=fixed_atoms, vacuum_gaps=vacuum_gaps or {})

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

        from dataclasses import replace as _replace

        if axes is not None:
            requirements = _replace(requirements, periodic_axes=axes)
            if any(axes):
                requirements = requirements.requiring(
                    "periodic",
                    f"the input structure is periodic along {_axis_names(axes)}",
                )
            if any(axes) and not all(axes):
                requirements = requirements.requiring(
                    "partial_periodic",
                    f"the input structure is periodic along {_axis_names(axes)} only; "
                    "each axis must be honoured independently, not collapsed to "
                    "all-or-nothing periodicity",
                )

        if fixed_atoms and self.task in ("md", "optimize"):
            requirements = requirements.requiring(
                "fixed_atoms",
                f"the structure freezes {fixed_atoms} atom(s) with FixAtoms, which "
                "velocity initialisation, integration and the temperature degrees of "
                "freedom must respect",
            )

        if self.energy_convention is not None:
            requirements = _replace(requirements, energy_convention=self.energy_convention)
        return requirements

    def _check_geometry(
        self,
        axes: tuple[bool, bool, bool] | None,
        *,
        fixed_atoms: int,
        vacuum_gaps: Mapping[int, float],
    ) -> None:
        """Refuse stress/NPT requests the structure's geometry cannot support."""
        npt = self.ensemble == "npt"
        if fixed_atoms and npt:
            raise ConfigError(
                f"npt with {fixed_atoms} frozen atom(s) is refused on every engine: the "
                "barostat's pressure and kinetic terms would include atoms that cannot "
                "move. Remove the FixAtoms constraint or use nvt."
            )
        if axes is None:
            return
        in_plane = npt and self.barostat_coupling == "in-plane"
        if in_plane:
            if sum(axes) < 2:
                raise ConfigError(
                    "barostat_coupling = 'in-plane' needs at least two periodic axes to "
                    f"barostat; this structure is periodic along {_axis_names(axes)}"
                )
            return
        if (npt or self.compute_stress) and not all(axes):
            what = "the npt ensemble" if npt else "simulation.compute_stress"
            raise ConfigError(
                f"{what} needs a cell periodic along all three axes: a stress tensor "
                "is normalised by the cell volume and the virial is only defined for "
                f"periodic directions. This structure is periodic along "
                f"{_axis_names(axes)}."
                + (
                    " For a slab, barostat_coupling = 'in-plane' barostats only the "
                    "periodic in-plane directions on engines that support it."
                    if npt and sum(axes) == 2
                    else ""
                )
            )
        if npt and vacuum_gaps:
            gaps = ", ".join(
                f"{width:.2f} Angstrom along lattice vector {'abc'[axis]}"
                for axis, width in sorted(vacuum_gaps.items())
            )
            raise ConfigError(
                f"npt with barostat_coupling = {self.barostat_coupling or 'engine default'!r} "
                f"on a fully periodic cell that contains an empty gap ({gaps}; threshold "
                f"{self.resolved_vacuum_gap_threshold_angstrom:g} Angstrom): this looks "
                "like a vacuum slab, and an isotropic or anisotropic barostat would scale "
                "the vacuum together with the solid. Use barostat_coupling = 'in-plane' on "
                "an engine that supports it, remove the vacuum, or raise "
                "simulation.vacuum_gap_threshold_angstrom if the gap is intended."
            )

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
            "barostat_coupling": self.barostat_coupling,
            "vacuum_gap_threshold_angstrom": self.vacuum_gap_threshold_angstrom,
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
    should still be visible -- except for an in-plane NPT run, whose
    structure is not periodic along all three axes and so has no full stress
    tensor to report.
    """
    if simulation.task == "singlepoint":
        return simulation
    from dataclasses import replace

    npt_stress = simulation.ensemble == "npt" and simulation.barostat_coupling != "in-plane"
    return replace(
        simulation,
        task="singlepoint",
        ensemble=None,
        steps=0,
        timestep_fs=None,
        thermostat=None,
        barostat=None,
        thermostat_damping_fs=None,
        barostat_damping_fs=None,
        barostat_coupling=None,
        compute_stress=simulation.compute_stress or npt_stress,
    )


def resolve_seed(seed: int | None) -> int:
    """The seed an engine actually uses: ``seed``, or a freshly drawn one.

    Engines call this once per run and record the returned value, so a run
    with ``seed = None`` is still reproducible from its manifest.
    """
    if seed is not None:
        return int(seed)
    import secrets

    return secrets.randbelow(MAX_SEED) + 1


def _axis_names(axes: Sequence[bool]) -> str:
    """``(True, True, False)`` -> ``"cell axes a, b"``; pbc is per cell vector."""
    names = [name for name, periodic in zip("abc", axes) if periodic]
    if not names:
        return "no cell axis"
    return ("cell axis " if len(names) == 1 else "cell axes ") + ", ".join(names)


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
                    "pbc": list(self.structure.pbc) if self.structure.pbc else None,
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
    "BAROSTAT_COUPLINGS",
    "DEFAULT_THERMOSTAT_DAMPING_FS",
    "DEFAULT_BAROSTAT_DAMPING_FS",
    "DEFAULT_VACUUM_GAP_THRESHOLD_ANGSTROM",
    "RUNTIMES",
    "PLATFORM_PRECISIONS",
    "OPENMM_PLATFORMS",
    "OPENMM_GPU_PLATFORMS",
    "MAX_SEED",
    "resolve_seed",
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
