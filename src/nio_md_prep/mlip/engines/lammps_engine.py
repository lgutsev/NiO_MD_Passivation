"""The LAMMPS engine: geometry, deck rendering, execution and read-back.

Rendering is separated from running on purpose. The exact ``pair_style`` and
``pair_coeff`` strings -- and the whole input deck around them -- are what a
LAMMPS run actually *is*, so :func:`render_deck` is pure, testable, and
callable on a machine with no LAMMPS at all. ``mlip validate`` records the
pair commands (with model files renamed to their staged names), the resolved
execution route and the thermostat/barostat mapping; the complete deck of a
run is written to ``in.lammps`` in the job directory before LAMMPS starts and
stored verbatim in the result's ``native`` record.

Geometry is passed through without being changed:

* ``boundary`` is set per axis (``p`` or ``f``) from ``atoms.pbc``. A
  non-periodic axis uses a *fixed* face (``f``), so the box volume is the
  source cell volume, and an atom outside ``[0, 1)`` along such an axis is
  refused before LAMMPS runs (:func:`prepare_geometry`) rather than lost,
  wrapped or shrink-wrapped.
* Every cell is written in LAMMPS's restricted-triclinic form through ASE's
  ``Prism`` (a rotation, never a shear), and forces, the stress tensor,
  positions and the final cell are rotated back into the source basis
  before anything is returned.

Execution has two routes, chosen by ``engine.runtime`` (:func:`resolve_launch`):

* the LAMMPS python module, whose per-atom data are read with
  ``gather``/``gather_atoms`` -- ordered by atom ID, ghost atoms excluded,
  correct under MPI -- after checking the atom count and the IDs;
* an ``lmp`` executable (optionally behind ``engine.mpi_launcher``), whose
  results are read from dumps written with ``sort id`` and ``%.17g`` and from
  ``print`` lines, frame by frame, with the atom count, IDs and timesteps
  checked. An explicitly configured executable, launcher or timeout selects
  this route; it is never silently replaced by an installed module.

Units: LAMMPS ``metal`` already speaks eV and Angstrom, ``real`` speaks
kcal/mol, and LAMMPS reports *pressure* (bar for ``metal``, atm for ``real``)
where this subsystem's canonical stress is eV/Angstrom^3 with the opposite
sign. The stress is the *virial-only* pressure (``compute pressure NULL
virial``) -- the configurational stress an ASE calculator reports, without the
kinetic term -- converted by
:func:`~nio_md_prep.mlip.units.lammps_pressure_tensor_to_stress`. Barostat
targets are converted the other way with
:func:`~nio_md_prep.mlip.units.pressure_to_lammps`.
"""
from __future__ import annotations

import contextlib
import functools
import json
import os
import shutil
import subprocess
import sys
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, executable_available, module_available
from ..errors import (
    ConfigError,
    MissingDependencyError,
    MlipError,
    ModelIntegrityError,
    ResultError,
)
from ..specs import (
    LammpsMlipPotentialSpec,
    SimulationSpec,
    path_key,
)
from ..units import (
    LAMMPS_METAL,
    lammps_pressure_tensor_to_stress,
    lammps_pressure_unit,
    lammps_unit_system,
    pressure_to_lammps,
)
from .base import EngineRuntime

DATA_FILE = "structure.lmp"
DECK_FILE = "in.lammps"
DUMP_FILE = "forces.dump"
THERMO_FILE = "thermo_result.txt"
SERIES_FILE = "thermo_series.dat"
TRAJECTORY_FILE = "smoke_md.lammpstrj"
#: The trajectory converted to the source basis (element symbols, source cell).
TRAJECTORY_EXTXYZ = "smoke_md.extxyz"
LOG_FILE = "log.lammps"
#: stdout/stderr of the executable route, kept in the job directory.
STDOUT_FILE = "lammps.stdout"
STDERR_FILE = "lammps.stderr"

#: LAMMPS timestep units per unit style. ``metal`` counts picoseconds, ``real``
#: femtoseconds; getting this wrong is a factor of 1000 in the dynamics.
TIMESTEP_PER_FS = {"metal": 1e-3, "real": 1.0}
#: Thermostat/barostat damping is quoted in the same time unit as the timestep.
DAMPING_PER_FS = TIMESTEP_PER_FS
#: Femtoseconds per LAMMPS time unit. Times are *divided* by this, so e.g.
#: 1.5 fs renders as exactly ``0.0015`` rather than ``1.5 * 1e-3``.
FS_PER_TIME_UNIT = {"metal": 1000.0, "real": 1.0}

#: What has been verified for each LAMMPS unit style, and against what: ``metal``
#: against an independent numpy Lennard-Jones reference (energy, forces, virial
#: stress, every geometry class), ``real`` against the same system run in
#: ``metal`` units (energy, forces and stress equal to 1e-7 relative; NVE, NVT
#: with every thermostat and NPT with both barostats following the same
#: trajectory to 1e-8 Angstrom, which a bar-for-atm mix-up fails), on both
#: execution routes: ``tests/test_mlip_lammps_physics.py``. A request needing
#: anything not listed for its unit style is refused by :func:`check_units`.
VERIFIED_UNIT_STYLES = {
    "metal": frozenset({"energy", "forces", "stress", "nve", "nvt", "npt"}),
    "real": frozenset({"energy", "forces", "stress", "nve", "nvt", "npt"}),
}

#: Thermo keywords read back after every run (``step`` included, so the
#: executable route can report what LAMMPS actually completed).
THERMO_KEYS = (
    "step", "pe", "ke", "etotal", "temp", "press", "ecouple", "econserve",
    "lx", "ly", "lz", "xy", "xz", "yz",
)
#: Per-step series written during MD by ``fix print`` (column order).
SERIES_KEYS = ("step", "pe", "ke", "etotal", "ecouple", "econserve", "temp")
VIRIAL_COMPUTE = "mlip_virial"
#: ``compute pressure`` global vector order: xx, yy, zz, xy, xz, yz.
VIRIAL_COMPONENTS = ("xx", "yy", "zz", "xy", "xz", "yz")
PE_ATOM_COMPUTE = "mlip_pe_atom"
UNWRAPPED_COMPUTE = "mlip_unwrapped"
TEMP_COMPUTE = "mlip_temp"
FROZEN_GROUP = "mlip_frozen"
MOBILE_GROUP = "mlip_mobile"

#: ``simulation.thermostat`` -> how LAMMPS implements it. Every entry was run
#: against the pinned LAMMPS build; anything else is refused at validate time.
THERMOSTATS = {
    "nose-hoover": "fix nvt (Nose-Hoover chain, LAMMPS default tchain 3)",
    "langevin": "fix nve + fix langevin (zero yes, tally yes)",
    "berendsen": "fix nve + fix temp/berendsen",
    "csvr": "fix nve + fix temp/csvr (Bussi-Donadio-Parrinello)",
}
DEFAULT_THERMOSTAT = "nose-hoover"
#: ``simulation.barostat`` -> how LAMMPS implements it.
BAROSTATS = {
    "mtk": "fix npt / fix nph (Martyna-Tobias-Klein, Nose-Hoover chains)",
    "berendsen": "fix press/berendsen on top of the thermostat's integrator",
}
DEFAULT_BAROSTAT = "mtk"
#: The one engine option the Berendsen barostat needs: LAMMPS's own default
#: modulus (10 pressure units) is an LJ value, far too soft for a solid.
BERENDSEN_MODULUS_OPTION = "barostat_bulk_modulus_bar"

#: KOKKOS/GPU command-line arguments for MACE ML-IAP (MACE docs).
KOKKOS_GPU_ARGS = ("-k", "on", "g", "1", "-sf", "kk", "-pk", "kokkos", "newton", "on", "neigh", "half")
#: Command-line switches that already choose an accelerator/suffix. Combined
#: with ``engine.threads`` they would be ambiguous, so that is refused.
_ACCELERATOR_SWITCHES = ("-sf", "-suffix", "-pk", "-package", "-k", "-kokkos")

#: Tilts smaller than this (Angstrom) are rounding noise from the rotation and
#: are written as an orthogonal box.
_TILT_NOISE_ANGSTROM = 1e-10
#: Largest acceptable round-trip error of the positions LAMMPS read back.
_POSITION_ROUNDTRIP_ANGSTROM = 1e-6


class LammpsEngine(EngineRuntime):
    kind = "lammps"
    requires = ()
    native_units = LAMMPS_METAL.name

    def capabilities(self) -> CapabilitySet:
        """What the LAMMPS engine carries; whether a pair style supplies it is separate.

        ``per_atom_energy`` means only that the deck can ask for ``compute
        pe/atom`` (and checks its sum against ``pe``); a route claims it only
        when its pair style is known to tally per-atom energies. ``gpu`` is
        true only when the resolved command line actually puts the pair style
        on a GPU (KOKKOS ``-k on g N -sf kk`` or the GPU package), never
        merely because a LAMMPS build could.
        """
        try:
            gpu = resolve_launch(self.spec).accelerator["gpu"]
        except ConfigError:
            gpu = False
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=True,
            per_atom_energy=True,
            periodic=True,
            # Per-axis boundary (p/f), atoms outside a fixed face refused.
            partial_periodic=True,
            # FixAtoms -> mlip_frozen/mlip_mobile groups (see render_deck).
            fixed_atoms=True,
            gpu=gpu,
            elements=None,
            precisions=None,
            engines=frozenset({"lammps"}),
            native_units=self.native_units,
            notes=(
                "LAMMPS reports pressure (bar in metal units, atm in real units); this "
                "subsystem returns the virial-only stress in eV/Angstrom^3 with the sign "
                "flipped, rotated back to the source basis",
            ),
        )

    def availability(self) -> Availability:
        try:
            launch = resolve_launch(self.spec)
        except ConfigError as exc:
            return Availability(False, (), str(exc))
        return launch_availability(launch)


# ---------------------------------------------------------------------------
# Execution route
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LammpsLaunch:
    """How LAMMPS is started, resolved once from an :class:`EngineSpec`.

    ``route`` is ``"python"`` (the in-process module) or ``"executable"``.
    ``cmdargs`` are LAMMPS command-line arguments, used on both routes.
    ``activate`` names the ML-IAP python coupling to activate in the python
    route (``"mliappy"`` / ``"mliappy_kokkos"``). ``required`` lists build
    features, ``("package", "PYTHON")`` or ``("pair", "mace")``, checked
    against the LAMMPS build before anything runs.
    """

    route: str
    executable: str | None = None
    prefix: tuple[str, ...] = ()
    cmdargs: tuple[str, ...] = ()
    env: Mapping[str, str] = field(default_factory=dict)
    timeout_s: float | None = None
    threads: int | None = None
    activate: str | None = None
    required: tuple[tuple[str, str], ...] = ()
    reason: str = ""
    notes: tuple[str, ...] = ()

    def argv(self, deck: str = DECK_FILE, log: str = LOG_FILE) -> list[str] | None:
        """The executable route's exact argument vector (``None`` for the python route)."""
        if self.route != "executable":
            return None
        return [*self.prefix, str(self.executable), *self.cmdargs, "-in", deck, "-log", log]

    def python_cmdargs(self, log: str = LOG_FILE) -> list[str]:
        """The ``cmdargs`` the python route hands to ``lammps.lammps``."""
        return [*self.cmdargs, "-log", log, "-screen", "none"]

    def command_line(self) -> str:
        """One human-readable line saying how LAMMPS is started, as recorded."""
        env = " ".join(f"{k}={v if v else repr(v)}" for k, v in sorted(self.env.items()))
        prefix = f"{env} " if env else ""
        if self.route == "executable":
            return prefix + subprocess.list2cmdline(self.argv())
        steps = [f"lmp = lammps.lammps(cmdargs={self.python_cmdargs()!r})"]
        if self.activate:
            steps.append(f"lammps.mliap.activate_{self.activate}(lmp)")
        steps.append(f"lmp.commands_list(<{DECK_FILE}>)")
        return prefix + "python: " + "; ".join(steps)

    @property
    def accelerator(self) -> dict:
        """What the command-line switches select: KOKKOS (host or GPU), GPU package, OMP."""
        return parse_accelerator_args(self.cmdargs)

    def as_dict(self) -> dict:
        accelerator = self.accelerator
        return {
            "route": self.route,
            "executable": self.executable,
            "mpi_launcher": list(self.prefix),
            "cmdargs": list(self.cmdargs),
            "argv": self.argv(),
            "command_line": self.command_line(),
            "kokkos_args": list(accelerator["kokkos_args"]),
            "accelerator": accelerator,
            "env": dict(self.env),
            "timeout_s": self.timeout_s,
            "threads": self.threads,
            "activate": self.activate,
            "required_build_features": [f"{kind}:{name}" for kind, name in self.required],
            "reason": self.reason,
            "notes": list(self.notes),
        }


def parse_accelerator_args(cmdargs: Sequence[str]) -> dict:
    """Read the accelerator a LAMMPS command line selects, without guessing.

    ``-k on g 1 t 4`` enables KOKKOS (``g N`` GPUs per node, ``t N`` host
    threads), ``-sf kk|gpu|omp`` sets the style suffix and ``-pk <package>
    ...`` configures a package. ``gpu`` is true only when the switches
    actually put the pair style on a GPU: KOKKOS enabled with ``g`` >= 1 and
    the ``kk`` suffix, or the GPU package with the ``gpu`` suffix. Whether the
    LAMMPS build *has* a GPU backend is a separate, build-level question.
    """
    tokens = list(cmdargs)
    kokkos_on = False
    kokkos_opts: dict[str, str] = {}
    kokkos_args: list[str] = []
    suffix = None
    packages: dict[str, list[str]] = {}
    i = 0

    def trailing(start: int) -> list[str]:
        values = []
        while start < len(tokens) and not str(tokens[start]).startswith("-"):
            values.append(str(tokens[start]))
            start += 1
        return values

    while i < len(tokens):
        token = str(tokens[i])
        if token in ("-k", "-kokkos"):
            values = trailing(i + 1)
            kokkos_args += [token, *values]
            if values:
                kokkos_on = values[0] == "on"
                rest = values[1:]
                kokkos_opts.update(dict(zip(rest[0::2], rest[1::2])))
            i += 1 + len(values)
        elif token in ("-sf", "-suffix"):
            values = trailing(i + 1)[:1]
            suffix = values[0] if values else None
            if suffix == "kk":
                kokkos_args += [token, *values]
            i += 1 + len(values)
        elif token in ("-pk", "-package"):
            values = trailing(i + 1)
            if values:
                packages[values[0]] = values[1:]
                if values[0] == "kokkos":
                    kokkos_args += [token, *values]
            i += 1 + len(values)
        else:
            i += 1
    try:
        kokkos_gpus = int(kokkos_opts.get("g", "0"))
    except ValueError:
        kokkos_gpus = 0
    kokkos_gpu = kokkos_on and suffix == "kk" and kokkos_gpus >= 1
    gpu_package = suffix == "gpu" or "gpu" in packages
    return {
        "kokkos": kokkos_on,
        "kokkos_gpus": kokkos_gpus if kokkos_on else 0,
        "kokkos_threads": kokkos_opts.get("t") if kokkos_on else None,
        "kokkos_args": kokkos_args,
        "suffix": suffix,
        "packages": packages,
        "gpu": bool(kokkos_gpu or gpu_package),
        "gpu_via": "KOKKOS" if kokkos_gpu else ("GPU package" if gpu_package else None),
    }


def resolve_launch(
    engine_spec,
    *,
    extra_args: Sequence[str] = (),
    activate: str | None = None,
    required: Sequence[tuple[str, str]] = (),
    notes: Sequence[str] = (),
    env: Mapping[str, str] | None = None,
    threads_in_extra_args: bool = False,
) -> LammpsLaunch:
    """Resolve ``engine.runtime`` into one concrete, recorded way of starting LAMMPS.

    ``auto`` picks the executable route when ``executable``, ``mpi_launcher``
    or ``timeout_s`` is set (only that route can honour them), otherwise the
    python module when it is importable, otherwise ``lmp`` on ``PATH``.
    ``engine.threads`` becomes OpenMP threads (``-pk omp N -sf omp``, plus
    ``OMP_NUM_THREADS``); it is refused together with explicit accelerator
    switches in ``lammps_args``. A caller whose ``extra_args`` already carry
    the thread count (KOKKOS ``t N``) passes ``threads_in_extra_args``.
    ``extra_args`` are the route's defaults and are replaced, not extended, by
    a non-empty ``engine.lammps_args``. ``env`` is set for the LAMMPS process
    on both routes (for the python route, around the in-process run).
    """
    runtime = getattr(engine_spec, "runtime", "auto") if engine_spec else "auto"
    executable = getattr(engine_spec, "executable", None) if engine_spec else None
    prefix = tuple(getattr(engine_spec, "mpi_launcher", ()) or ()) if engine_spec else ()
    timeout = getattr(engine_spec, "timeout_s", None) if engine_spec else None
    user_args = tuple(getattr(engine_spec, "lammps_args", ()) or ()) if engine_spec else ()
    threads = getattr(engine_spec, "threads", None) if engine_spec else None

    if runtime == "auto":
        if executable or prefix or timeout:
            route, reason = "executable", (
                "runtime auto: engine.executable/mpi_launcher/timeout_s is set, which only "
                "the executable route honours"
            )
        elif module_available("lammps"):
            route, reason = "python", "runtime auto: the LAMMPS python module is importable"
        else:
            route, reason = "executable", (
                "runtime auto: the LAMMPS python module is not importable"
            )
    else:
        route, reason = runtime, f"runtime = {runtime!r} was requested explicitly"

    cmdargs = list(user_args)
    launch_env: dict[str, str] = {str(k): str(v) for k, v in dict(env or {}).items()}
    if threads and not (threads_in_extra_args and not user_args):
        clashing = [arg for arg in user_args if arg in _ACCELERATOR_SWITCHES]
        if not user_args:
            clashing += [arg for arg in extra_args if arg in _ACCELERATOR_SWITCHES]
        if clashing:
            raise ConfigError(
                f"engine.threads = {threads} would add '-pk omp {threads} -sf omp', but the "
                f"LAMMPS arguments already choose an accelerator ({', '.join(clashing)}); "
                "remove engine.threads or put the thread settings in engine.lammps_args"
            )
        cmdargs += ["-pk", "omp", str(threads), "-sf", "omp"]
    if threads:
        launch_env["OMP_NUM_THREADS"] = str(threads)
    if not user_args:
        cmdargs += list(extra_args)
    return LammpsLaunch(
        route=route,
        executable=(executable or "lmp") if route == "executable" else None,
        prefix=prefix if route == "executable" else (),
        cmdargs=tuple(cmdargs),
        env=launch_env,
        timeout_s=timeout if route == "executable" else None,
        threads=threads,
        activate=activate,
        required=tuple(required),
        reason=reason,
        notes=tuple(notes),
    )


def launch_availability(launch: LammpsLaunch) -> Availability:
    """Whether ``launch`` can start LAMMPS here, and -- when it names required
    build features -- whether the build has them (probed in a subprocess, so
    no LAMMPS is imported into this process)."""
    if launch.route == "python":
        if not module_available("lammps"):
            return Availability(
                False,
                ("lammps",),
                "engine.runtime = 'python' needs the LAMMPS python module, which is not "
                "importable",
            )
        where = "LAMMPS python module"
    else:
        if not _executable_found(launch.executable):
            return Availability(
                False,
                ("lammps",),
                f"the LAMMPS executable {launch.executable!r} was not found",
            )
        missing_launcher = launch.prefix and not _executable_found(launch.prefix[0])
        if missing_launcher:
            return Availability(
                False, (), f"the MPI launcher {launch.prefix[0]!r} was not found"
            )
        where = f"{launch.executable} executable"
    if not launch.required:
        return Availability(True, (), where)
    try:
        build = probe_build(launch.route, launch.executable)
    except MlipError as exc:
        return Availability(False, (), f"could not probe the LAMMPS build: {exc}")
    missing = missing_features(build, launch.required)
    if missing:
        return Availability(False, (), missing_features_message(missing, build))
    return Availability(True, (), f"{where}; build has {', '.join(n for _, n in launch.required)}")


def _executable_found(name: str | None) -> bool:
    return bool(name) and (executable_available(name) or Path(name).is_file())


#: KOKKOS back ends that put data on a GPU (``accelerator_config`` / ``lmp -h`` names).
KOKKOS_GPU_BACKENDS = frozenset({"cuda", "hip", "sycl"})


@functools.lru_cache(maxsize=8)
def probe_build(route: str, executable: str | None = None) -> dict:
    """Installed packages, pair styles and accelerator back ends of a LAMMPS build.

    Runs without importing LAMMPS here. Python route: a fresh interpreter
    creates a LAMMPS instance and reports ``installed_packages``,
    ``available_styles("pair")`` and ``accelerator_config``. Executable route:
    ``lmp -h`` is parsed ("Installed packages:", "* Pair styles:" and the
    "<PACKAGE> package API:" lines of "Accelerator configuration:").
    """
    if route == "python":
        code = (
            "import json, lammps\n"
            "l = lammps.lammps(cmdargs=['-screen', 'none', '-log', 'none'])\n"
            "acc = {k: {'api': list(v.get('api', [])), 'precision': list(v.get('precision', []))}\n"
            "       for k, v in dict(l.accelerator_config).items()}\n"
            "out = {'version': l.version(), 'packages': list(l.installed_packages),\n"
            "       'pair_styles': sorted(l.available_styles('pair')), 'accelerators': acc}\n"
            "l.close()\n"
            "print(json.dumps(out))\n"
        )
        try:
            completed = subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True, text=True, timeout=300, check=False,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise MlipError(f"probing the LAMMPS python module failed: {exc}") from exc
        if completed.returncode != 0:
            raise MlipError(completed.stderr.strip()[-500:] or "the probe exited non-zero")
        data = json.loads(completed.stdout.strip().splitlines()[-1])
        packages = sorted(data.get("packages", []))
        accelerators = {
            name: {
                "api": sorted(str(a).lower() for a in info.get("api", [])),
                "precision": list(info.get("precision", [])),
            }
            for name, info in dict(data.get("accelerators", {})).items()
            if name in packages
        }
        return {
            "route": route,
            "version": data.get("version"),
            "packages": packages,
            "pair_styles": sorted(data.get("pair_styles", [])),
            "accelerators": accelerators,
        }
    try:
        completed = subprocess.run(
            [executable or "lmp", "-h"], capture_output=True, text=True, timeout=300, check=False
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise MlipError(f"running '{executable} -h' failed: {exc}") from exc
    return {"route": route, **parse_help(completed.stdout)}


def parse_help(text: str) -> dict:
    """Packages, pair styles and accelerator back ends from ``lmp -h`` output."""
    lines = text.splitlines()
    version = next((line.strip() for line in lines if line.strip()), None)
    packages: list[str] = []
    pair_styles: list[str] = []
    accelerators: dict[str, dict] = {}
    section = None
    for line in lines:
        stripped = line.strip()
        if " package API:" in stripped:
            name, _, api = stripped.partition(" package API:")
            accelerators.setdefault(name.strip(), {})["api"] = sorted(
                token.lower() for token in api.split()
            )
            continue
        if " package precision:" in stripped:
            name, _, precision = stripped.partition(" package precision:")
            accelerators.setdefault(name.strip(), {})["precision"] = precision.split()
            continue
        if stripped.startswith("Installed packages:"):
            section = "packages"
            continue
        if stripped.startswith("List of individual style options"):
            section = None
            continue
        if stripped.startswith("* "):
            section = "pair" if stripped.startswith("* Pair styles:") else None
            continue
        if section == "packages":
            packages += stripped.split()
        elif section == "pair":
            pair_styles += stripped.split()
    return {
        "version": version,
        "packages": sorted(set(packages)),
        "pair_styles": sorted(set(pair_styles)),
        "accelerators": accelerators,
    }


def missing_features(build: Mapping, required: Sequence[tuple[str, str]]) -> list[tuple[str, str]]:
    """Required build features the probed build lacks.

    Kinds: ``("package", NAME)``, ``("pair", STYLE)`` and
    ``("kokkos_backend", "gpu" | "host-only")`` -- a KOKKOS build whose back
    ends include a GPU (CUDA/HIP/SYCL), or one with host back ends only (a
    GPU-enabled KOKKOS build does not run without a GPU).
    """
    packages = set(build.get("packages", ()))
    pair_styles = set(build.get("pair_styles", ()))
    kokkos_api = set((build.get("accelerators", {}).get("KOKKOS") or {}).get("api", ()))
    missing = []
    for kind, name in required:
        if kind == "package":
            present = name in packages
        elif kind == "pair":
            present = name in pair_styles
        elif kind == "kokkos_backend":
            if "KOKKOS" not in packages:
                present = True  # reported once, as the missing KOKKOS package
            elif name == "gpu":
                present = bool(kokkos_api & KOKKOS_GPU_BACKENDS)
            else:
                present = bool(kokkos_api) and not kokkos_api & KOKKOS_GPU_BACKENDS
        else:
            raise ValueError(f"unknown build feature kind {kind!r}")
        if not present:
            missing.append((kind, name))
    return missing


def missing_features_message(missing, build: Mapping) -> str:
    labels = {"package": "package", "pair": "pair style", "kokkos_backend": "KOKKOS back end"}
    listing = ", ".join(f"{labels.get(k, k)} {n}" for k, n in missing)
    kokkos_api = (build.get("accelerators", {}).get("KOKKOS") or {}).get("api")
    detail = f"; its KOKKOS API is {' '.join(kokkos_api)}" if kokkos_api else ""
    return (
        f"this LAMMPS build ({build.get('version') or 'unknown version'}, {build.get('route')} "
        f"route) lacks {listing}{detail}"
    )


# ---------------------------------------------------------------------------
# Geometry (numpy + ASE; no LAMMPS)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LammpsGeometry:
    """A structure as LAMMPS will see it, and the rotation back to the source basis.

    ``rotation`` is ASE ``Prism.rot_mat``: a row vector ``v`` in the LAMMPS
    frame is ``v @ rotation.T`` in the source frame; a tensor ``T`` is
    ``rotation @ T @ rotation.T``.
    """

    pbc: tuple[bool, bool, bool]
    boundary: str
    rotation: tuple[tuple[float, float, float], ...]
    lammps_cell: tuple[tuple[float, float, float], ...]
    triclinic: bool
    prism: object = field(repr=False, compare=False, default=None)

    def as_dict(self) -> dict:
        return {
            "pbc": list(self.pbc),
            "boundary": self.boundary,
            "triclinic": self.triclinic,
            "lammps_cell_angstrom": [list(row) for row in self.lammps_cell],
            "rotation_to_source": [list(row) for row in self.rotation],
            "convention": "restricted triclinic via ASE Prism (rotation only)",
        }


def boundary_string(pbc: Sequence[bool]) -> str:
    """``(True, True, False)`` -> ``"p p f"``: fixed, never shrink-wrapped, faces."""
    return " ".join("p" if p else "f" for p in pbc)


def prepare_geometry(atoms) -> LammpsGeometry:
    """Refuse geometry LAMMPS would change, then build the restricted-triclinic frame.

    A rank-3 cell is required (none is synthesised for a cluster). Along a
    non-periodic axis every atom must lie in ``[0, 1)`` fractional, because a
    fixed LAMMPS face would otherwise lose it; the offending indices are
    listed instead of the structure being wrapped or shifted.
    """
    import numpy as np
    from ase.calculators.lammps.coordinatetransform import Prism

    from ..structures import outside_nonperiodic_faces, periodic_axes

    pbc = periodic_axes(atoms)
    cell = np.array(atoms.get_cell(), dtype=float)
    if atoms.cell.rank < 3 or abs(np.linalg.det(cell)) == 0.0:
        raise ConfigError(
            "LAMMPS needs a rank-3 simulation box, and this structure has none "
            f"(cell {cell.tolist()}). Give it an explicit cell -- for a cluster, a box "
            "larger than the cluster plus the cutoff, with pbc false -- none is synthesised."
        )
    outside = outside_nonperiodic_faces(atoms)
    if outside:
        detail = "; ".join(
            f"cell axis {'abc'[axis]}: {len(indices)} atom(s), e.g. indices "
            f"{list(indices[:8])}"
            for axis, indices in sorted(outside.items())
        )
        raise ConfigError(
            "atoms lie outside the cell along a non-periodic axis, where LAMMPS uses a "
            f"fixed boundary ('f') and would lose them ({detail}). Centre the structure "
            "in its cell (e.g. atoms.center(axis=...)) or enlarge the cell along that "
            "axis; nothing is wrapped or shrink-wrapped automatically."
        )
    prism = Prism(cell, pbc)
    lammps_cell = np.array(prism.lammps_tilt, dtype=float)
    tilts = lammps_cell[(1, 2, 2), (0, 0, 1)]
    return LammpsGeometry(
        pbc=pbc,
        boundary=boundary_string(pbc),
        rotation=tuple(tuple(float(v) for v in row) for row in prism.rot_mat),
        lammps_cell=tuple(tuple(float(v) for v in row) for row in lammps_cell),
        triclinic=bool(np.abs(tilts).max() > _TILT_NOISE_ANGSTROM),
        prism=prism,
    )


def vacuum_axes(atoms, simulation: SimulationSpec) -> tuple[int, ...]:
    """Lattice axes of a fully periodic cell whose empty gap exceeds the threshold."""
    from ..structures import periodic_axes, vacuum_gaps

    if not all(periodic_axes(atoms)):
        return ()
    gaps = vacuum_gaps(atoms, threshold_angstrom=simulation.resolved_vacuum_gap_threshold_angstrom)
    return tuple(sorted(gaps))


def write_data_file(
    atoms, potential: LammpsMlipPotentialSpec, path: Path, *, geometry: LammpsGeometry | None = None
) -> Path:
    """Write the structure as a LAMMPS data file in the potential's type order.

    The per-atom ``type`` column is always written explicitly from
    ``potential.type_map`` -- an ``Atoms`` object read from a LAMMPS data file
    carries its own ``type`` array, which ASE's writer would otherwise prefer
    over ``specorder`` and so silently mislabel elements. Positions are
    written in the frame of ``geometry`` (restricted triclinic, rotated with
    ASE ``Prism``), unwrapped, at 17 significant digits; a cell with a tilt
    above rounding noise is always written as triclinic.
    """
    import numpy as np
    from ase.io.lammpsdata import write_lammps_data

    from ..structures import to_lammps_type_order

    types = to_lammps_type_order(atoms, potential.type_map)  # validates coverage first
    specorder = [potential.type_map[t] for t in sorted(potential.type_map)]
    masses = np.asarray(atoms.get_masses(), dtype=float)
    types_array = np.asarray(types, dtype=int)
    for type_id in sorted(set(types)):
        values = masses[types_array == type_id]
        if values.size and float(values.max() - values.min()) > 0.0:
            raise ConfigError(
                f"atoms of LAMMPS type {type_id} ({potential.type_map[type_id]}) carry "
                f"different masses ({sorted(set(values.tolist()))[:4]}); a LAMMPS type has "
                "one mass, so this structure cannot be written without changing it"
            )
    geometry = geometry or prepare_geometry(atoms)
    staged = atoms.copy()
    staged.set_constraint()
    staged.arrays["type"] = types_array
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        write_lammps_data(
            handle,
            staged,
            specorder=specorder,
            masses=True,
            units=potential.units,
            atom_style=potential.atom_style,
            prismobj=geometry.prism,
            force_skew=geometry.triclinic,
        )
    return path


# ---------------------------------------------------------------------------
# Model files: staging, hashing and pair-command rewriting
# ---------------------------------------------------------------------------


def staged_model_names(potential: LammpsMlipPotentialSpec) -> dict[str, str]:
    """``path_key(model file) -> file name inside the job directory``.

    The staged name is the file's own name with spaces replaced (LAMMPS
    splits arguments on whitespace). Two model files with the same name are
    refused: a pair command that names one of them would be ambiguous.
    """
    names: dict[str, str] = {}
    owner: dict[str, str] = {}
    for model in potential.model_paths:
        name = Path(model).name.replace(" ", "_")
        key = path_key(model)
        if name in owner and owner[name] != key:
            raise ConfigError(
                f"potential.model_files lists two files named {Path(model).name!r}; each "
                "model file is staged into the job directory under its own name, so the "
                "names must be distinct"
            )
        owner[name] = key
        names[key] = name
    return names


def _rewrite_text(text: str, by_name: Mapping[str, str]) -> str:
    tokens = text.split()
    out = []
    for token in tokens:
        bare = token.strip("\"'")
        staged = by_name.get(Path(bare).name) if bare else None
        out.append(staged if staged is not None else token)
    return " ".join(out)


def rewrite_model_tokens(
    potential: LammpsMlipPotentialSpec, *, targets: Mapping[str, str] | None = None
) -> LammpsMlipPotentialSpec:
    """The potential with every model-file argument renamed to its staged name.

    A ``pair_style``/``pair_coeff``/extra-command token names a model file
    when its file name equals the name of a listed ``model_files`` entry
    (``nio.pb``, ``models/nio.pb`` and ``/abs/nio.pb`` all do). LAMMPS runs
    in the job directory, so the staged name is what it must open -- never a
    file of the same name found through ``LAMMPS_POTENTIALS``. ``targets``
    (``path_key(model) -> token``) replaces the staged names, e.g. with
    absolute paths for a LAMMPS instance that does not run in the job
    directory.
    """
    if not potential.model_paths:
        return potential
    names = dict(targets) if targets is not None else staged_model_names(potential)
    by_name = {Path(model).name: names[path_key(model)] for model in potential.model_paths}
    return replace(
        potential,
        pair_style=_rewrite_text(potential.pair_style, by_name),
        pair_coeff=tuple(_rewrite_text(c, by_name) for c in potential.pair_coeff),
        extra_commands=tuple(_rewrite_text(c, by_name) for c in potential.extra_commands),
    )


def _command_tokens(commands: Sequence[str]) -> list[str]:
    """LAMMPS arguments of ``commands``, double quotes removed."""
    import shlex

    text = " ".join(commands)
    try:
        return shlex.split(text, posix=True)
    except ValueError:
        return text.split()


def lammps_path_token(path: Path) -> str:
    """An absolute path as one LAMMPS argument (forward slashes, quoted if needed)."""
    text = Path(path).resolve().as_posix()
    return f'"{text}"' if any(ch.isspace() for ch in text) else text


def stage_model_files(
    potential: LammpsMlipPotentialSpec, workdir: Path, *, absolute: bool = False
) -> tuple[LammpsMlipPotentialSpec, list[dict]]:
    """Copy (or hard-link) every model file into ``workdir`` and verify it.

    Returns the potential with rewritten pair commands and one record per
    file (source, staged path, both SHA256s, method). A declared
    ``model_hashes`` entry must match the staged file. With ``absolute``, the
    pair commands name the staged files by absolute path (for a LAMMPS that
    does not run inside ``workdir``, e.g. ASE's LAMMPSlib).
    """
    from ..provenance import sha256_file

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    names = staged_model_names(potential)
    declared = {path_key(k): v for k, v in potential.model_hashes.items()}
    targets = (
        {key: lammps_path_token(workdir / name) for key, name in names.items()}
        if absolute
        else None
    )
    rewritten = rewrite_model_tokens(potential, targets=targets)
    commands = _command_tokens(rewritten.render_pair_commands())
    records = []
    for model in potential.model_paths:
        source = Path(model)
        if not source.is_file():
            raise MissingDependencyError(
                f"the model file {source}", (), hint="potential.model_files names a file that does not exist"
            )
        target = workdir / names[path_key(source)]
        source_sha = sha256_file(source)
        method = "already staged"
        if path_key(target) != path_key(source):
            if target.exists():
                target.unlink()
            try:
                os.link(source, target)
                method = "hard link"
            except OSError:
                shutil.copy2(source, target)
                method = "copy"
        staged_sha = sha256_file(target)
        if staged_sha != source_sha:
            raise ModelIntegrityError(target, source_sha, staged_sha)
        expected = declared.get(path_key(source))
        if expected is not None and expected.lower() != staged_sha:
            raise ModelIntegrityError(source, expected, staged_sha)
        records.append(
            {
                "source": str(source),
                "staged": str(target),
                "staged_name": target.name,
                "sha256": staged_sha,
                "declared_sha256": expected,
                "method": method,
                "referenced_by_pair_commands": (
                    (targets[path_key(source)].strip('"') if absolute else target.name)
                    in commands
                ),
            }
        )
    return rewritten, records


# ---------------------------------------------------------------------------
# Thermostat / barostat mapping (pure)
# ---------------------------------------------------------------------------


def _num(value) -> str:
    """Shortest exact decimal for a number (``repr`` of the float)."""
    value = float(value)
    if value == int(value) and abs(value) < 1e15:
        return str(int(value))
    return repr(value)


def time_in_units(fs: float, units: str) -> float:
    return float(fs) / FS_PER_TIME_UNIT[units]


def id_ranges(ids: Sequence[int]) -> list[str]:
    """``[1, 2, 3, 7, 9, 10]`` -> ``["1:3", "7", "9:10"]`` (LAMMPS ``group id`` syntax)."""
    values = sorted(set(int(i) for i in ids))
    ranges: list[str] = []
    start = previous = None
    for value in values:
        if start is None:
            start = previous = value
        elif value == previous + 1:
            previous = value
        else:
            ranges.append(f"{start}:{previous}" if previous > start else str(start))
            start = previous = value
    if start is not None:
        ranges.append(f"{start}:{previous}" if previous > start else str(start))
    return ranges


def barostat_dimensions(simulation: SimulationSpec, pbc: Sequence[bool], vacuum: Sequence[int]):
    """``(style, dims, coupling_name)`` for the barostat on this geometry.

    ``style`` is ``"iso"``/``"aniso"`` (all three dimensions) or
    ``"couple"`` (two named dimensions coupled together, ``dims``). LAMMPS
    barostats only periodic dimensions, and an isotropic/anisotropic barostat
    on a cell with a vacuum gap would scale the vacuum; both are refused.
    """
    coupling = simulation.barostat_coupling
    periodic = [axis for axis in range(3) if pbc[axis]]
    if coupling == "in-plane":
        if len(periodic) == 2:
            dims = tuple(periodic)
        elif len(periodic) == 3 and len(vacuum) == 1:
            dims = tuple(axis for axis in range(3) if axis != vacuum[0])
        else:
            raise ConfigError(
                "barostat_coupling = 'in-plane' needs a plane to barostat: either exactly "
                "two periodic axes (a slab with pbc like [true, true, false]) or a fully "
                "periodic cell with exactly one vacuum gap wider than "
                f"{simulation.resolved_vacuum_gap_threshold_angstrom:g} Angstrom; this "
                f"structure has periodic axes {['abc'[a] for a in periodic]} and vacuum "
                f"along {['abc'[a] for a in vacuum] or 'none'}"
            )
        return "couple", dims, "in-plane"
    if len(periodic) != 3 or vacuum:
        raise ConfigError(
            f"barostat_coupling = {coupling or 'isotropic (default)'!r} barostats all three "
            "cell dimensions, which needs a fully periodic cell without a vacuum gap; use "
            "barostat_coupling = 'in-plane' for a slab"
        )
    if coupling == "anisotropic":
        return "aniso", (0, 1, 2), "anisotropic"
    return "iso", (0, 1, 2), "isotropic"


def plan_md(
    simulation: SimulationSpec,
    *,
    units: str,
    pbc: Sequence[bool] = (True, True, True),
    vacuum: Sequence[int] = (),
    n_fixed: int = 0,
    seed: int = 1,
    options: Mapping | None = None,
    triclinic: bool = False,
) -> dict:
    """The fixes implementing ``simulation`` in LAMMPS, or a :class:`ConfigError`.

    Returns ``{"commands": [...], "record": {...}}``. ``record`` is what
    :attr:`TrajectoryResult.integrator_resolved` reports: the fixes, the
    requested and resolved thermostat/barostat, damping times (fs and LAMMPS
    units), coupling, targets in both unit systems, the seed, and whether the
    run has a conserved quantity. Pure; used at validate time too.
    """
    options = dict(options or {})
    ensemble = simulation.ensemble
    group = MOBILE_GROUP if n_fixed else "all"
    t_damp = time_in_units(simulation.resolved_thermostat_damping_fs, units)
    record: dict = {
        "engine": "lammps",
        "ensemble": ensemble,
        "group": group,
        "temperature_compute": TEMP_COMPUTE,
        "timestep_fs": simulation.timestep_fs,
        "timestep_lammps": time_in_units(simulation.timestep_fs, units),
        "lammps_time_unit": "ps" if units == "metal" else "fs",
        "seed": seed,
    }
    commands: list[str] = []
    temperature = simulation.temperature_K
    if ensemble == "nve":
        commands.append(f"fix mlip_integrate {group} nve")
        record.update(
            fixes=["nve"],
            thermostat=None,
            description="fix nve (velocity Verlet)",
            conserved_quantity="etotal",
        )
        return {"commands": commands, "record": record}

    requested = simulation.thermostat
    thermostat = requested or DEFAULT_THERMOSTAT
    if thermostat not in THERMOSTATS:
        raise ConfigError(
            f"thermostat {thermostat!r} is not implemented for the LAMMPS engine; "
            f"available: {', '.join(THERMOSTATS)}"
        )
    T = _num(temperature)
    record.update(
        thermostat_requested=requested,
        thermostat=thermostat,
        thermostat_default_used=requested is None,
        thermostat_damping_fs=simulation.resolved_thermostat_damping_fs,
        thermostat_damping_lammps=t_damp,
        temperature_K=temperature,
        velocity_seed=seed,
    )
    thermostat_fixes: list[str] = []
    fix_modify: list[str] = []
    if thermostat == "langevin":
        thermostat_fixes.append(
            f"fix mlip_thermostat {group} langevin {T} {T} {_num(t_damp)} {seed} zero yes tally yes"
        )
        record["thermostat_seed"] = seed
    elif thermostat == "berendsen":
        thermostat_fixes.append(f"fix mlip_thermostat {group} temp/berendsen {T} {T} {_num(t_damp)}")
        fix_modify.append(f"fix_modify mlip_thermostat temp {TEMP_COMPUTE}")
    elif thermostat == "csvr":
        thermostat_fixes.append(f"fix mlip_thermostat {group} temp/csvr {T} {T} {_num(t_damp)} {seed}")
        fix_modify.append(f"fix_modify mlip_thermostat temp {TEMP_COMPUTE}")
        record["thermostat_seed"] = seed

    if ensemble == "nvt":
        if thermostat == "nose-hoover":
            commands.append(f"fix mlip_integrate {group} nvt temp {T} {T} {_num(t_damp)}")
            commands.append(f"fix_modify mlip_integrate temp {TEMP_COMPUTE}")
            fixes = ["nvt"]
        else:
            commands.append(f"fix mlip_integrate {group} nve")
            commands += thermostat_fixes + fix_modify
            fixes = ["nve", thermostat_fixes[0].split()[3]]
        record.update(
            fixes=fixes,
            description=THERMOSTATS[thermostat],
            conserved_quantity="econserve",
        )
        return {"commands": commands, "record": record}

    # --- npt -------------------------------------------------------------
    if n_fixed:
        raise ConfigError(
            "npt with frozen atoms is refused on the LAMMPS engine: fix npt's pressure and "
            "temperature computes act on group all"
        )
    requested_barostat = simulation.barostat
    barostat = requested_barostat or DEFAULT_BAROSTAT
    if barostat not in BAROSTATS:
        raise ConfigError(
            f"barostat {barostat!r} is not implemented for the LAMMPS engine (no exactly "
            f"equivalent LAMMPS fix is used here); available: {', '.join(BAROSTATS)}"
        )
    style, dims, coupling = barostat_dimensions(simulation, pbc, vacuum)
    pressure_native = pressure_to_lammps(simulation.pressure_bar, units)
    p_damp = time_in_units(simulation.resolved_barostat_damping_fs, units)
    P, Pd = _num(pressure_native), _num(p_damp)
    if style in ("iso", "aniso"):
        coupling_args = f"{style} {P} {P} {Pd}"
        barostatted = ["x", "y", "z"]
    else:
        names = ["xyz"[d] for d in dims]
        coupling_args = (
            f"{names[0]} {P} {P} {Pd} {names[1]} {P} {P} {Pd} couple {''.join(names)}"
        )
        barostatted = names
    pressure_unit = lammps_pressure_unit(units)[0]
    record.update(
        barostat_requested=requested_barostat,
        barostat=barostat,
        barostat_default_used=requested_barostat is None,
        barostat_coupling_requested=simulation.barostat_coupling,
        barostat_coupling=coupling,
        barostatted_lammps_dimensions=barostatted,
        barostat_damping_fs=simulation.resolved_barostat_damping_fs,
        barostat_damping_lammps=p_damp,
        pressure_bar=simulation.pressure_bar,
        pressure_lammps=pressure_native,
        pressure_lammps_unit=pressure_unit,
    )
    if barostat == "mtk":
        if thermostat == "nose-hoover":
            commands.append(
                f"fix mlip_integrate all npt temp {T} {T} {_num(t_damp)} {coupling_args}"
            )
            fixes = ["npt"]
        else:
            commands.append(f"fix mlip_integrate all nph {coupling_args}")
            commands += thermostat_fixes + fix_modify
            fixes = ["nph", thermostat_fixes[0].split()[3]]
        conserved = "econserve"
        description = f"{BAROSTATS['mtk']}; thermostat: {THERMOSTATS[thermostat]}"
    else:
        if triclinic:
            raise ConfigError(
                "barostat = 'berendsen' is refused for a triclinic cell on the LAMMPS engine: "
                "fix press/berendsen rejects triclinic boxes ('Cannot use fix press/berendsen "
                "with triclinic box'); use barostat = 'mtk' (fix npt/nph handles tilts)"
            )
        modulus_bar = options.get(BERENDSEN_MODULUS_OPTION)
        if (
            isinstance(modulus_bar, bool)
            or not isinstance(modulus_bar, (int, float))
            or not modulus_bar > 0
        ):
            raise ConfigError(
                "barostat = 'berendsen' on LAMMPS needs the system's bulk modulus: set "
                f"engine.options.{BERENDSEN_MODULUS_OPTION} (bar, e.g. 1.9e6 for NiO); LAMMPS's "
                "own default of 10 pressure units is a Lennard-Jones value"
            )
        modulus_native = pressure_to_lammps(modulus_bar, units)
        if thermostat == "nose-hoover":
            commands.append(f"fix mlip_integrate all nvt temp {T} {T} {_num(t_damp)}")
            fixes = ["nvt"]
        else:
            commands.append("fix mlip_integrate all nve")
            commands += thermostat_fixes + fix_modify
            fixes = ["nve", thermostat_fixes[0].split()[3]]
        commands.append(
            f"fix mlip_barostat all press/berendsen {coupling_args} modulus {_num(modulus_native)}"
        )
        fixes.append("press/berendsen")
        record.update(bulk_modulus_bar=float(modulus_bar), bulk_modulus_lammps=modulus_native)
        # press/berendsen's work on the box is not tallied into ecouple.
        conserved = None
        description = f"{BAROSTATS['berendsen']}; thermostat: {THERMOSTATS[thermostat]}"
    record.update(fixes=fixes, description=description, conserved_quantity=conserved)
    return {"commands": commands, "record": record}


# ---------------------------------------------------------------------------
# Deck rendering (pure; no LAMMPS required)
# ---------------------------------------------------------------------------


def wants_stress(simulation: SimulationSpec) -> bool:
    """Whether the job's single-point form reports a stress (see ``as_singlepoint``)."""
    from ..specs import as_singlepoint

    return bool(as_singlepoint(simulation).compute_stress)


def render_deck(
    potential: LammpsMlipPotentialSpec,
    simulation: SimulationSpec,
    *,
    n_types: int,
    pbc: Sequence[bool] = (True, True, True),
    fixed_ids: Sequence[int] = (),
    n_atoms: int | None = None,
    seed: int = 1,
    vacuum: Sequence[int] = (),
    options: Mapping | None = None,
    triclinic: bool = False,
    data_file: str = DATA_FILE,
    dump_file: str = DUMP_FILE,
    thermo_file: str = THERMO_FILE,
    series_file: str = SERIES_FILE,
    trajectory_file: str = TRAJECTORY_FILE,
) -> tuple[str, ...]:
    """Render the complete input deck as an ordered tuple of commands.

    Pure: no files are written and no LAMMPS is needed. The pair commands come
    from :meth:`LammpsMlipPotentialSpec.render_pair_commands`, so what is
    validated, what is stored in the manifest and what is executed are the
    same strings. ``fixed_ids`` are 1-based LAMMPS atom IDs frozen by
    ``FixAtoms``; they are excluded from velocity creation, integration, the
    thermostat and the temperature compute. The per-atom energy compute and
    the virial compute are emitted only when the job asks for them.
    """
    if n_types != len(potential.type_map):
        raise ValueError(
            f"deck asks for {n_types} atom types but the potential maps "
            f"{len(potential.type_map)}"
        )
    pbc = tuple(bool(p) for p in pbc)
    md = simulation.task == "md"
    per_atom = bool(simulation.compute_per_atom_energy)
    stress = wants_stress(simulation)
    lines = [f"# generated by nio-md-prep mlip ({potential.label})", "clear"]
    if potential.newton:
        lines.append(f"newton {potential.newton}")
    lines += [
        f"units {potential.units}",
        f"atom_style {potential.atom_style}",
        f"boundary {boundary_string(pbc)}",
        "atom_modify map yes",
        f"read_data {data_file}",
    ]
    lines.extend(potential.render_pair_commands())
    lines.append("")
    group = "all"
    if fixed_ids:
        if n_atoms is not None and len(set(fixed_ids)) >= n_atoms:
            raise ConfigError("every atom is frozen by FixAtoms; there is nothing to integrate")
        lines.append(f"group {FROZEN_GROUP} id {' '.join(id_ranges(fixed_ids))}")
        lines.append(f"group {MOBILE_GROUP} subtract all {FROZEN_GROUP}")
        group = MOBILE_GROUP
    lines.append(f"compute {TEMP_COMPUTE} {group} temp")
    if fixed_ids:
        # 3 (N - N_frozen) kinetic DOF: with frozen atoms the centre of mass is
        # not a free coordinate (diagnostics.temperature_ndof).
        lines.append(f"compute_modify {TEMP_COMPUTE} extra/dof 0")
    lines.append(f"compute {UNWRAPPED_COMPUTE} all property/atom xu yu zu")
    if per_atom:
        lines.append(f"compute {PE_ATOM_COMPUTE} all pe/atom")
    thermo = list(THERMO_KEYS)
    if stress:
        lines.append(f"compute {VIRIAL_COMPUTE} all pressure NULL virial")
        thermo += [f"c_{VIRIAL_COMPUTE}[{i}]" for i in range(1, 7)]
    lines.append("thermo_style custom " + " ".join(thermo))
    lines.append(f"thermo_modify temp {TEMP_COMPUTE} norm no")

    if md:
        lines.extend(
            _md_block(
                potential,
                simulation,
                group=group,
                frozen=bool(fixed_ids),
                pbc=pbc,
                vacuum=vacuum,
                n_fixed=len(set(fixed_ids)),
                seed=seed,
                options=options,
                triclinic=triclinic,
                series_file=series_file,
                trajectory_file=trajectory_file,
            )
        )
    lines.extend(_endpoint_block(dump_file, per_atom=per_atom, post=md))
    lines.extend(_thermo_print_block(thermo_file, stress=stress))
    return tuple(lines)


def _series_format() -> str:
    return " ".join(
        "$(step)" if key == "step" else f"$({key}:%.17g)" for key in SERIES_KEYS
    )


def _md_block(
    potential: LammpsMlipPotentialSpec,
    simulation: SimulationSpec,
    *,
    group: str,
    frozen: bool,
    pbc,
    vacuum,
    n_fixed: int,
    seed: int,
    options,
    triclinic: bool,
    series_file: str,
    trajectory_file: str,
) -> list[str]:
    units = potential.units
    plan = plan_md(
        simulation, units=units, pbc=pbc, vacuum=vacuum, n_fixed=n_fixed, seed=seed,
        options=options, triclinic=triclinic,
    )
    log_every = max(1, simulation.log_interval)
    traj_every = max(1, simulation.trajectory_interval)
    elements = " ".join(potential.type_map[t] for t in sorted(potential.type_map))
    lines = [f"timestep {_num(time_in_units(simulation.timestep_fs, units))}"]
    if simulation.temperature_K:
        # mom yes removes the centre-of-mass momentum when nothing is frozen;
        # rot no everywhere (periodic images have no rigid-body rotation).
        lines.append(
            f"velocity {group} create {_num(simulation.temperature_K)} {seed} "
            f"mom {'no' if frozen else 'yes'} rot no dist gaussian temp {TEMP_COMPUTE}"
        )
    lines.extend(plan["commands"])
    lines += [
        f"thermo {log_every}",
        f'fix mlip_series all print {log_every} "{_series_format()}" file {series_file} '
        f'screen no title "# {" ".join(SERIES_KEYS)}"',
        f"dump mlip_traj all custom {traj_every} {trajectory_file} id type xu yu zu",
        f"dump_modify mlip_traj sort id format float %.17g element {elements}",
        f"run {simulation.steps}",
        "undump mlip_traj",
    ]
    if simulation.steps % traj_every:
        # The last step is not a multiple of the dump interval: append it, so
        # the trajectory always ends at the final geometry.
        lines.append(
            f"write_dump all custom {trajectory_file} id type xu yu zu modify sort id "
            f"format float %.17g element {elements} append yes"
        )
    # Close fix print's file before appending the final step to it.
    lines.append("unfix mlip_series")
    if simulation.steps % log_every:
        lines.append(f'print "{_series_format()}" append {series_file} screen no')
    # The final single point must be the potential alone: thermostat forces
    # (Langevin friction/noise) are in f after a step, so every fix is removed
    # and the final geometry re-evaluated with run 0.
    for fix_id in _fix_ids(plan["commands"]):
        lines.append(f"unfix {fix_id}")
    return lines


def _fix_ids(commands: Sequence[str]) -> list[str]:
    ids = []
    for command in commands:
        parts = command.split()
        if parts and parts[0] == "fix" and parts[1] not in ids:
            ids.append(parts[1])
    return list(reversed(ids))


def _endpoint_block(dump_file: str, *, per_atom: bool, post: bool) -> list[str]:
    columns = "id type xu yu zu fx fy fz" + (f" c_{PE_ATOM_COMPUTE}" if per_atom else "")
    return [
        "thermo 1",
        f"dump mlip_forces all custom 1 {dump_file} {columns}",
        "dump_modify mlip_forces sort id format float %.17g",
        "run 0 post no" if post else "run 0",
        "undump mlip_forces",
    ]


def _thermo_print_block(thermo_file: str, *, stress: bool) -> list[str]:
    """Print the thermo quantities read back, one ``name value`` per line at %.17g.

    Used by the executable route, which has no way to reach into LAMMPS's
    memory. Harmless under the Python route.
    """
    entries = [(key, key) for key in THERMO_KEYS]
    if stress:
        entries += [
            (f"virial_{name}", f"c_{VIRIAL_COMPUTE}[{index}]")
            for index, name in enumerate(VIRIAL_COMPONENTS, start=1)
        ]
    lines = []
    for index, (label, expression) in enumerate(entries):
        redirect = "file" if index == 0 else "append"
        lines.append(f'print "{label} $({expression}:%.17g)" {redirect} {thermo_file} screen no')
    return lines


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def run_deck(
    commands,
    *,
    workdir: Path,
    launch: LammpsLaunch,
    n_atoms: int,
    want_per_atom: bool = False,
    want_stress: bool = False,
) -> dict:
    """Execute a rendered deck and read the final state back in the LAMMPS frame.

    Always writes the deck to ``workdir`` first, whichever route runs it, so a
    failed job leaves behind exactly the input that failed. Returns numpy
    arrays ordered by atom ID: ``forces``, ``positions`` (unwrapped), and
    ``per_atom_energy`` when requested, plus ``thermo``, ``virial`` (LAMMPS
    order xx, yy, zz, xy, xz, yz, when requested), ``box`` and
    ``steps_completed``.
    """
    # Absolute: the python route changes into it, and error messages read the
    # log through this path from inside it.
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    deck_path = workdir / DECK_FILE
    deck_path.write_text("\n".join(commands) + "\n", encoding="utf-8")
    for stale in (
        THERMO_FILE, DUMP_FILE, SERIES_FILE, TRAJECTORY_FILE, TRAJECTORY_EXTXYZ,
        STDOUT_FILE, STDERR_FILE,
    ):
        with contextlib.suppress(FileNotFoundError):
            (workdir / stale).unlink()

    availability = launch_availability(launch) if launch.route == "executable" else None
    if availability is not None and not availability:
        raise MissingDependencyError(
            "running a LAMMPS deck",
            availability.missing or ("lammps",),
            hint=f"{availability.detail}. The rendered deck was written to {deck_path}.",
        )
    if launch.route == "python":
        if not module_available("lammps"):
            raise MissingDependencyError(
                "running a LAMMPS deck through the python module",
                ("lammps",),
                hint=f"The rendered deck was written to {deck_path}; run it with an 'lmp' "
                "executable (engine.runtime = 'executable') or install the lammps module.",
            )
        raw = _run_with_module(commands, workdir, launch, n_atoms, want_per_atom, want_stress)
    else:
        raw = _run_with_executable(workdir, deck_path, launch, n_atoms, want_per_atom, want_stress)
    raw.update(
        deck_path=str(deck_path),
        log_path=str(workdir / LOG_FILE),
        launch=launch.as_dict(),
    )
    return raw


def _check_build_in_process(lmp, launch: LammpsLaunch) -> dict:
    packages = sorted(lmp.installed_packages)
    build = {
        "route": "python",
        "version": lmp.version(),
        "packages": packages,
        "pair_styles": sorted(
            name for kind, name in launch.required if kind == "pair" and lmp.has_style("pair", name)
        ),
        "accelerators": {
            name: {"api": sorted(str(a).lower() for a in info.get("api", []))}
            for name, info in dict(lmp.accelerator_config).items()
            if name in packages
        },
    }
    missing = missing_features(build, launch.required)
    if missing:
        raise MissingDependencyError(
            "this LAMMPS route",
            (),
            hint=missing_features_message(missing, build)
            + ". Use a LAMMPS build with these packages/styles.",
        )
    return build


def _run_with_module(commands, workdir, launch, n_atoms, want_per_atom, want_stress) -> dict:
    import numpy as np
    from lammps import LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, lammps

    started = time.perf_counter()
    # The deck names its data, dump and thermo files relatively, so the run has
    # to happen in the job directory (process-global chdir; not thread-safe).
    # `shell cd` would break on any path with a space in it. launch.env is set
    # around the run (process-global too) and restored afterwards.
    with _working_directory(workdir), _environment(launch.env):
        lmp = lammps(cmdargs=launch.python_cmdargs())
        try:
            build = _check_build_in_process(lmp, launch) if launch.required else {
                "route": "python", "version": lmp.version()
            }
            if launch.activate:
                _activate_mliappy(lmp, launch.activate)
            try:
                lmp.commands_list(list(commands))
            except Exception as exc:  # LAMMPS raises a bare Exception with its message
                raise MlipError(
                    f"LAMMPS failed in {workdir}: {str(exc).strip()}\n"
                    + _log_tail(workdir / LOG_FILE)
                ) from exc
            natoms = int(lmp.get_natoms())
            if natoms != n_atoms:
                raise ResultError(
                    f"LAMMPS holds {natoms} atoms after the run but the structure has "
                    f"{n_atoms} (atoms lost through a fixed boundary?)"
                )
            ids = np.ctypeslib.as_array(lmp.gather_atoms("id", 0, 1)).astype(int).copy()
            if not np.array_equal(ids, np.arange(1, n_atoms + 1)):
                raise ResultError("LAMMPS atom IDs are not the consecutive 1..N this deck wrote")
            forces = np.ctypeslib.as_array(lmp.gather("f", 1, 3)).reshape(n_atoms, 3).copy()
            positions = (
                np.ctypeslib.as_array(lmp.gather(f"c_{UNWRAPPED_COMPUTE}", 1, 3))
                .reshape(n_atoms, 3)
                .copy()
            )
            per_atom = None
            if want_per_atom:
                per_atom = (
                    np.ctypeslib.as_array(lmp.gather(f"c_{PE_ATOM_COMPUTE}", 1, 1))
                    .reshape(n_atoms)
                    .copy()
                )
            thermo = {key: float(lmp.get_thermo(key)) for key in THERMO_KEYS}
            virial = None
            if want_stress:
                virial = np.array(
                    lmp.numpy.extract_compute(VIRIAL_COMPUTE, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR),
                    dtype=float,
                )[:6].copy()
            boxlo, _boxhi, xy, yz, xz, _periodicity, _changed = lmp.extract_box()
            steps = int(lmp.extract_global("ntimestep"))
        finally:
            lmp.close()
    return {
        "route": "python-module",
        "build": build,
        "thermo": thermo,
        "forces": forces,
        "positions": positions,
        "per_atom_energy": per_atom,
        "virial": virial,
        "box": _box(boxlo, thermo, xy=xy, xz=xz, yz=yz),
        "steps_completed": steps,
        "wall_time_s": time.perf_counter() - started,
    }


def _activate_mliappy(lmp, which: str) -> None:
    try:
        if which == "mliappy_kokkos":
            from lammps.mliap import activate_mliappy_kokkos as activate
        else:
            from lammps.mliap import activate_mliappy as activate
        activate(lmp)
    except Exception as exc:
        raise MissingDependencyError(
            f"activating the ML-IAP python coupling ({which})",
            ("lammps.mliap",),
            hint=f"{type(exc).__name__}: {exc}",
        ) from exc


def _run_with_executable(workdir, deck_path, launch, n_atoms, want_per_atom, want_stress) -> dict:
    import numpy as np

    started = time.perf_counter()
    command = launch.argv(deck=deck_path.name, log=LOG_FILE)
    env = {**os.environ, **dict(launch.env)} if launch.env else None
    try:
        completed = subprocess.run(
            command,
            cwd=str(workdir),
            capture_output=True,
            text=True,
            check=False,
            timeout=launch.timeout_s,
            env=env,
        )
    except subprocess.TimeoutExpired as exc:
        raise MlipError(
            f"{' '.join(command)} did not finish within engine.timeout_s = {launch.timeout_s} s "
            f"in {workdir}\n" + _log_tail(workdir / LOG_FILE)
        ) from exc
    # Kept verbatim: some pair styles (pair_mace) report their device on the
    # process's stdout, not in the LAMMPS log.
    (workdir / STDOUT_FILE).write_text(completed.stdout or "", encoding="utf-8")
    (workdir / STDERR_FILE).write_text(completed.stderr or "", encoding="utf-8")
    if completed.returncode != 0:
        tail = (completed.stdout or completed.stderr or "").strip().splitlines()[-20:]
        raise MlipError(
            f"{launch.executable} exited with status {completed.returncode} for {deck_path}.\n"
            + "\n".join(tail)
        )
    thermo = read_thermo(workdir / THERMO_FILE)
    steps = int(round(thermo["step"]))
    frame = read_dump(workdir / DUMP_FILE, n_atoms=n_atoms)[-1]
    if frame.timestep != steps:
        raise ResultError(
            f"{workdir / DUMP_FILE} ends at timestep {frame.timestep}, but LAMMPS reported "
            f"step {steps}"
        )
    virial = None
    if want_stress:
        virial = np.array([thermo[f"virial_{name}"] for name in VIRIAL_COMPONENTS], dtype=float)
    per_atom = frame.column(f"c_{PE_ATOM_COMPUTE}") if want_per_atom else None
    return {
        "route": "executable",
        "build": {"route": "executable", "version": _log_version(workdir / LOG_FILE)},
        "thermo": thermo,
        "forces": frame.columns3("fx", "fy", "fz"),
        "positions": frame.columns3("xu", "yu", "zu"),
        "per_atom_energy": per_atom,
        "virial": virial,
        "box": _box(frame.box["lo"], thermo, xy=thermo["xy"], xz=thermo["xz"], yz=thermo["yz"]),
        "steps_completed": steps,
        "stdout_tail": (completed.stdout or "").strip().splitlines()[-5:],
        "stdout_path": str(workdir / STDOUT_FILE),
        "wall_time_s": time.perf_counter() - started,
    }


def _box(boxlo, thermo, *, xy, xz, yz) -> dict:
    return {
        "lo": [float(v) for v in boxlo],
        "lammps_cell": [
            [thermo["lx"], 0.0, 0.0],
            [float(xy), thermo["ly"], 0.0],
            [float(xz), float(yz), thermo["lz"]],
        ],
    }


@contextlib.contextmanager
def _working_directory(path: Path):
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


@contextlib.contextmanager
def _environment(values: Mapping[str, str]):
    """Set environment variables for the duration of an in-process run."""
    saved = {key: os.environ.get(key) for key in values}
    os.environ.update({key: str(value) for key, value in values.items()})
    try:
        yield
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _log_tail(path: Path, lines: int = 20) -> str:
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return f"(no LAMMPS log at {path})"
    return "\n".join(text.strip().splitlines()[-lines:])


def _log_version(path: Path) -> str | None:
    try:
        with Path(path).open(encoding="utf-8", errors="replace") as handle:
            first = handle.readline().strip()
    except OSError:
        return None
    return first or None


# ---------------------------------------------------------------------------
# Reading LAMMPS output
# ---------------------------------------------------------------------------


def read_thermo(path: Path) -> dict[str, float]:
    """``name value`` lines written by the deck's final print block."""
    path = Path(path)
    if not path.exists():
        raise ResultError(f"LAMMPS did not write {path}; the run did not reach its thermo block")
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            values[parts[0]] = float(parts[1])
        except ValueError:
            continue
    missing = [key for key in THERMO_KEYS if key not in values]
    if missing:
        raise ResultError(f"{path} lacks thermo value(s) {', '.join(missing)}")
    return values


@dataclass(frozen=True)
class DumpFrame:
    """One validated frame of a LAMMPS custom dump (rows sorted by atom ID)."""

    timestep: int
    n_atoms: int
    box: dict
    columns: tuple[str, ...]
    rows: tuple[tuple[str, ...], ...]

    def column(self, name: str):
        import numpy as np

        if name not in self.columns:
            raise ResultError(f"dump frame at timestep {self.timestep} has no column {name!r}")
        index = self.columns.index(name)
        return np.array([float(row[index]) for row in self.rows], dtype=float)

    def columns3(self, a: str, b: str, c: str):
        import numpy as np

        return np.stack([self.column(a), self.column(b), self.column(c)], axis=1)


def read_dump(path: Path, *, n_atoms: int, expected_timesteps: Sequence[int] | None = None) -> list[DumpFrame]:
    """Read every frame of a custom dump, validating each one.

    Each frame must have ``NUMBER OF ATOMS == n_atoms``, a three-line box, an
    ``ATOMS`` header with an ``id`` column, exactly ``n_atoms`` complete rows
    and IDs ``1..N`` in order (``dump_modify sort id``). ``expected_timesteps``,
    when given, must equal the frames' timesteps exactly. Anything else -- a
    missing file, a truncated frame, an extra or absent frame -- is a
    :class:`ResultError`, never a shorter result.
    """
    path = Path(path)
    if not path.exists():
        raise ResultError(f"LAMMPS did not write {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    frames: list[DumpFrame] = []
    i = 0

    def fail(message: str):
        raise ResultError(f"{path}: frame {len(frames) + 1}: {message}")

    def expect(prefix: str) -> str:
        nonlocal i
        if i >= len(lines) or not lines[i].startswith(prefix):
            fail(f"expected {prefix!r} at line {i + 1}, got "
                 f"{lines[i]!r}" if i < len(lines) else f"expected {prefix!r}; the file ends")
        header = lines[i]
        i += 1
        return header

    def value_line() -> str:
        nonlocal i
        if i >= len(lines):
            fail("the file ends inside a frame header (truncated)")
        line = lines[i]
        i += 1
        return line

    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        expect("ITEM: TIMESTEP")
        try:
            timestep = int(value_line().split()[0])
        except (ValueError, IndexError):
            fail("the timestep is not an integer")
        expect("ITEM: NUMBER OF ATOMS")
        try:
            count = int(value_line().split()[0])
        except (ValueError, IndexError):
            fail("the atom count is not an integer")
        if count != n_atoms:
            fail(f"holds {count} atoms; the structure has {n_atoms}")
        box_header = expect("ITEM: BOX BOUNDS")
        box_lines = []
        for _ in range(3):
            try:
                box_lines.append([float(v) for v in value_line().split()])
            except ValueError:
                fail("a box-bounds line is not numeric")
        columns = tuple(expect("ITEM: ATOMS").split()[2:])
        if "id" not in columns:
            fail("the ATOMS section has no id column")
        rows = []
        for _ in range(count):
            if i >= len(lines) or lines[i].startswith("ITEM:"):
                fail(f"only {len(rows)} of {count} atom rows before the frame ends (truncated)")
            parts = tuple(lines[i].split())
            i += 1
            if len(parts) != len(columns):
                fail(f"an atom row has {len(parts)} fields for {len(columns)} columns (truncated)")
            rows.append(parts)
        id_index = columns.index("id")
        try:
            ids = [int(row[id_index]) for row in rows]
        except ValueError:
            fail("an atom id is not an integer")
        if ids != list(range(1, count + 1)):
            fail("atom IDs are not 1..N in order (missing 'dump_modify sort id', or duplicates)")
        frames.append(
            DumpFrame(
                timestep=timestep,
                n_atoms=count,
                box=_dump_box(box_header, box_lines),
                columns=columns,
                rows=tuple(rows),
            )
        )
    if not frames:
        raise ResultError(f"{path} contains no frames")
    if expected_timesteps is not None:
        found = [frame.timestep for frame in frames]
        if found != list(expected_timesteps):
            raise ResultError(
                f"{path} holds frames at timesteps {found[:6]}{'...' if len(found) > 6 else ''} "
                f"({len(found)} frames); the run should have written "
                f"{list(expected_timesteps)[:6]}{'...' if len(expected_timesteps) > 6 else ''} "
                f"({len(expected_timesteps)} frames)"
            )
    return frames


def _dump_box(header: str, rows: list[list[float]]) -> dict:
    """Dump box bounds -> ``lo`` corner and restricted-triclinic cell (LAMMPS frame)."""
    triclinic = "xy xz yz" in header
    if triclinic:
        (xlo_b, xhi_b, xy), (ylo_b, yhi_b, xz), (zlo_b, zhi_b, yz) = (r[:3] for r in rows)
        xlo = xlo_b - min(0.0, xy, xz, xy + xz)
        xhi = xhi_b - max(0.0, xy, xz, xy + xz)
        ylo = ylo_b - min(0.0, yz)
        yhi = yhi_b - max(0.0, yz)
        zlo, zhi = zlo_b, zhi_b
    else:
        (xlo, xhi), (ylo, yhi), (zlo, zhi) = (r[:2] for r in rows)
        xy = xz = yz = 0.0
    return {
        "lo": [xlo, ylo, zlo],
        "lammps_cell": [[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]],
        "boundary": header.split()[3:] if not triclinic else header.split()[6:],
    }


def read_series(path: Path, *, steps: int, interval: int) -> list[list[float]]:
    """The MD thermo series written by ``fix print``, one row per logged step.

    Rows are de-duplicated by step (the final step can be printed twice) and
    must cover exactly ``0, interval, 2*interval, ..., steps``.
    """
    path = Path(path)
    if not path.exists():
        raise ResultError(f"LAMMPS did not write the thermo series {path}")
    rows: dict[int, list[float]] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) != len(SERIES_KEYS):
            raise ResultError(f"{path}:{number}: {len(parts)} values for {len(SERIES_KEYS)} columns")
        try:
            values = [float(v) for v in parts]
        except ValueError:
            raise ResultError(f"{path}:{number}: non-numeric value in {line!r}") from None
        rows[int(round(values[0]))] = values
    expected = expected_steps(steps, interval)
    if sorted(rows) != expected:
        raise ResultError(
            f"{path} logged steps {sorted(rows)[:6]}...; expected {expected[:6]}... "
            f"({len(expected)} rows)"
        )
    return [rows[step] for step in expected]


def expected_steps(steps: int, interval: int) -> list[int]:
    """Steps a dump/print every ``interval`` writes over a run of ``steps``, plus the last."""
    values = list(range(0, steps + 1, max(1, interval)))
    if values[-1] != steps:
        values.append(steps)
    return values


# ---------------------------------------------------------------------------
# One job, end to end
# ---------------------------------------------------------------------------


def check_units(units: str, simulation: SimulationSpec | None) -> None:
    """Refuse a request whose quantities were not verified for this unit style."""
    verified = VERIFIED_UNIT_STYLES.get(units)
    if verified is None:
        raise ConfigError(
            f"LAMMPS units {units!r} are not supported: nothing has been verified for them "
            f"(verified: {', '.join(VERIFIED_UNIT_STYLES)})"
        )
    if simulation is None:
        return
    needed = {"energy", "forces"}
    if wants_stress(simulation):
        needed.add("stress")
    if simulation.task == "md":
        needed.add(simulation.ensemble)
    unverified = sorted(needed - verified)
    if unverified:
        raise ConfigError(
            f"LAMMPS units {units!r}: {', '.join(unverified)} has not been verified against a "
            "reference for this unit style and is refused; use units = 'metal'"
        )


def check_request(
    potential: LammpsMlipPotentialSpec,
    engine_spec,
    simulation: SimulationSpec,
    atoms=None,
    *,
    launch: LammpsLaunch | None = None,
) -> None:
    """Every LAMMPS-specific refusal, without running (or importing) LAMMPS.

    Launch options, thermostat/barostat support, and -- with a structure --
    the geometry preflight, frozen atoms and the barostat coupling for this
    geometry. Raises :class:`ConfigError`.
    """
    if launch is None:
        resolve_launch(engine_spec)
    check_units(potential.units, simulation)
    staged_model_names(potential)
    options = dict(getattr(engine_spec, "options", {}) or {})
    if atoms is None:
        if simulation.task == "md":
            if simulation.ensemble == "npt" and simulation.barostat_coupling == "in-plane":
                # Which dimensions are barostatted depends on the structure.
                probe = replace(simulation, barostat_coupling=None)
                plan_md(probe, units=potential.units, options=options)
            else:
                plan_md(simulation, units=potential.units, options=options)
        return
    from ..structures import fixed_atom_indices

    geometry = prepare_geometry(atoms)
    fixed = fixed_atom_indices(atoms)
    if simulation.task == "md":
        if fixed and len(fixed) >= len(atoms):
            raise ConfigError("every atom is frozen by FixAtoms; there is nothing to integrate")
        plan_md(
            simulation,
            units=potential.units,
            pbc=geometry.pbc,
            vacuum=vacuum_axes(atoms, simulation) if simulation.ensemble == "npt" else (),
            n_fixed=len(fixed),
            options=options,
            triclinic=geometry.triclinic,
        )


# ---------------------------------------------------------------------------
# Execution plan (what `mlip validate` shows; pure, runs nothing)
# ---------------------------------------------------------------------------

#: Placeholder for a velocity seed that is drawn only when the job runs.
SEED_AT_RUN_TIME = "<drawn at run time>"


def units_plan(units: str) -> dict:
    """The LAMMPS unit style and the conversions this subsystem applies to it."""
    system = lammps_unit_system(units)
    pressure_name, pressure_in_bar = lammps_pressure_unit(units)
    return {
        "native": system.name,
        "lammps_units": units,
        "energy_unit": system.energy,
        "length_unit": system.length,
        "time_unit": "ps" if units == "metal" else "fs",
        "pressure_unit": pressure_name,
        "pressure_unit_in_bar": pressure_in_bar,
        "energy_to_eV": system.energy_in_eV,
        "verified": sorted(VERIFIED_UNIT_STYLES.get(units, ())),
        "stress": (
            f"virial-only pressure tensor ({pressure_name}) -> eV/Angstrom^3, sign flipped, "
            "rotated to the source basis"
        ),
    }


#: Keys of every ``dynamics`` section of an execution plan (see :func:`dynamics_plan`).
DYNAMICS_PLAN_KEYS = (
    "ensemble", "integrator", "fixes", "commands", "thermostat", "thermostat_requested",
    "barostat", "barostat_requested", "barostat_coupling", "barostatted_lammps_dimensions",
    "timestep_fs", "timestep_native", "native_time_unit", "thermostat_damping_fs",
    "thermostat_damping_native", "barostat_damping_fs", "barostat_damping_native",
    "temperature_K", "pressure_bar", "pressure_native", "pressure_native_unit",
    "conserved_quantity", "velocity_seed", "group",
)


def dynamics_plan(
    simulation: SimulationSpec,
    *,
    units: str,
    atoms=None,
    options: Mapping | None = None,
) -> dict | None:
    """The integrator LAMMPS will run for ``simulation``, in user and native units.

    ``None`` for a task that integrates nothing. A request the LAMMPS engine
    refuses yields ``{"refused": message}`` (``validate`` raises the same
    message); without a structure an in-plane coupling cannot be resolved to
    LAMMPS dimensions and says so.
    """
    if simulation is None or simulation.task != "md":
        return None
    unresolved_coupling = False
    try:
        if atoms is not None:
            from ..structures import fixed_atom_indices

            geometry = prepare_geometry(atoms)
            fixed = fixed_atom_indices(atoms)
            plan = plan_md(
                simulation,
                units=units,
                pbc=geometry.pbc,
                vacuum=vacuum_axes(atoms, simulation) if simulation.ensemble == "npt" else (),
                n_fixed=len(fixed),
                seed=simulation.seed or SEED_AT_RUN_TIME,
                options=options,
                triclinic=geometry.triclinic,
            )
        else:
            probe = simulation
            if simulation.ensemble == "npt" and simulation.barostat_coupling == "in-plane":
                probe = replace(simulation, barostat_coupling=None)
                unresolved_coupling = True
            plan = plan_md(
                probe, units=units, seed=simulation.seed or SEED_AT_RUN_TIME, options=options
            )
    except ConfigError as exc:
        # Same keys as an accepted plan, so a reader never has to guess.
        refused = dict.fromkeys(DYNAMICS_PLAN_KEYS)
        refused.update(
            ensemble=simulation.ensemble,
            thermostat_requested=simulation.thermostat,
            barostat_requested=simulation.barostat,
            barostat_coupling=simulation.barostat_coupling,
            timestep_fs=simulation.timestep_fs,
            refused=str(exc),
        )
        return refused
    record = plan["record"]
    dynamics = {
        "ensemble": simulation.ensemble,
        "integrator": record["description"],
        "fixes": list(record["fixes"]),
        "commands": list(plan["commands"]),
        "thermostat": record.get("thermostat"),
        "thermostat_requested": simulation.thermostat,
        "barostat": record.get("barostat"),
        "barostat_requested": simulation.barostat,
        "barostat_coupling": (
            "in-plane (LAMMPS dimensions resolved from the structure)"
            if unresolved_coupling
            else record.get("barostat_coupling")
        ),
        "barostatted_lammps_dimensions": (
            None if unresolved_coupling else record.get("barostatted_lammps_dimensions")
        ),
        "timestep_fs": simulation.timestep_fs,
        "timestep_native": record["timestep_lammps"],
        "native_time_unit": record["lammps_time_unit"],
        "thermostat_damping_fs": record.get("thermostat_damping_fs"),
        "thermostat_damping_native": record.get("thermostat_damping_lammps"),
        "barostat_damping_fs": record.get("barostat_damping_fs"),
        "barostat_damping_native": record.get("barostat_damping_lammps"),
        "temperature_K": record.get("temperature_K"),
        "pressure_bar": record.get("pressure_bar"),
        "pressure_native": record.get("pressure_lammps"),
        "pressure_native_unit": record.get("pressure_lammps_unit"),
        "conserved_quantity": record.get("conserved_quantity"),
        "velocity_seed": simulation.seed if simulation.seed is not None else SEED_AT_RUN_TIME,
        "group": record.get("group"),
    }
    if "bulk_modulus_bar" in record:
        dynamics["bulk_modulus_bar"] = record["bulk_modulus_bar"]
        dynamics["bulk_modulus_native"] = record["bulk_modulus_lammps"]
    return dynamics


def lammps_plan(
    potential: LammpsMlipPotentialSpec, launch: LammpsLaunch | None, *, atoms=None
) -> dict:
    """The ``lammps`` section of an execution plan: rendered pair commands and launch."""
    rewritten = rewrite_model_tokens(potential)
    pair_coeff = [
        c if c.startswith("pair_coeff") else f"pair_coeff {c}" for c in rewritten.pair_coeff
    ]
    section = {
        "pair_style": f"pair_style {rewritten.pair_style}",
        "pair_coeff": pair_coeff,
        "pair_commands": list(rewritten.render_pair_commands()),
        "newton": potential.newton,
        "units": potential.units,
        "atom_style": potential.atom_style,
        "boundary": "per axis from structure.pbc: p (periodic) or f (fixed face)",
    }
    if atoms is not None:
        try:
            section["boundary"] = boundary_string(prepare_geometry(atoms).pbc)
        except ConfigError as exc:
            section["boundary"] = f"refused: {exc}"
    if launch is None:
        section.update(launch_command=None, launch_argv=None, kokkos_args=[], env={}, route=None)
        return section
    accelerator = launch.accelerator
    section.update(
        route=launch.route,
        launch_command=launch.command_line(),
        launch_argv=launch.argv(),
        cmdargs=list(launch.cmdargs),
        kokkos_args=list(accelerator["kokkos_args"]),
        accelerator=accelerator,
        activate=launch.activate,
        env=dict(launch.env),
        required_build_features=[f"{kind}:{name}" for kind, name in launch.required],
        launch_reason=launch.reason,
        notes=list(launch.notes),
    )
    return section


def deck_plan(
    potential: LammpsMlipPotentialSpec,
    simulation: SimulationSpec,
    atoms,
    *,
    options: Mapping | None = None,
) -> list[str] | None:
    """The complete deck a run of ``atoms`` would execute (seed placeholder if drawn)."""
    if atoms is None or simulation is None:
        return None
    from ..structures import fixed_atom_indices

    try:
        geometry = prepare_geometry(atoms)
        md = simulation.task == "md"
        fixed = fixed_atom_indices(atoms) if md else ()
        return list(
            render_deck(
                rewrite_model_tokens(potential),
                simulation,
                n_types=len(potential.type_map),
                pbc=geometry.pbc,
                fixed_ids=[i + 1 for i in fixed],
                n_atoms=len(atoms),
                seed=simulation.seed or SEED_AT_RUN_TIME,
                vacuum=(
                    vacuum_axes(atoms, simulation) if md and simulation.ensemble == "npt" else ()
                ),
                options=options,
                triclinic=geometry.triclinic,
            )
        )
    except ConfigError:
        return None


def model_file_plan(paths: Sequence[Path], declared: Mapping[str, str]) -> dict:
    """``{declared, observed}`` SHA256s of model files, keyed by path as configured."""
    from ..provenance import sha256_file

    declared_by_key = {path_key(k): v for k, v in dict(declared).items() if v}
    observed: dict[str, str | None] = {}
    declared_out: dict[str, str] = {}
    for path in paths:
        key = str(path)
        observed[key] = sha256_file(path) if Path(path).is_file() else None
        if path_key(path) in declared_by_key:
            declared_out[key] = declared_by_key[path_key(path)]
    return {"declared": declared_out, "observed": observed}


def bridge_plan_basics(bridge, simulation: SimulationSpec | None, atoms=None) -> dict:
    """Availability and unmet capabilities of ``bridge`` for this job (no run)."""
    availability = bridge.availability()
    unmet: list[str] = []
    if simulation is not None:
        capabilities = bridge.capabilities()
        try:
            unmet = list(bridge.requirements(simulation, atoms).unmet(capabilities))
        except ConfigError as exc:
            unmet = [f"refused: {exc}"]
        if atoms is not None and capabilities.elements is not None:
            uncovered = sorted(set(atoms.get_chemical_symbols()) - set(capabilities.elements))
            unmet += [f"element {symbol} is not covered by the potential" for symbol in uncovered]
    return {
        "availability": {
            "available": bool(availability),
            "missing": list(availability.missing),
            "detail": availability.detail,
        },
        "unmet_capabilities": unmet,
    }


def describe_request(potential: LammpsMlipPotentialSpec, engine_spec, simulation, launch) -> dict:
    """What ``mlip validate`` and the manifest record for a LAMMPS route."""
    rewritten = rewrite_model_tokens(potential)
    names = staged_model_names(potential)
    report = {
        "launch": launch.as_dict(),
        "boundary": "per axis from structure.pbc: p (periodic) or f (fixed face)",
        "pair_commands": list(rewritten.render_pair_commands()),
        "pair_commands_as_configured": list(potential.render_pair_commands()),
        "staged_model_files": [
            {"source": str(p), "staged_name": names[path_key(p)]} for p in potential.model_paths
        ],
        "stress": "virial only (compute pressure NULL virial), rotated to the source basis",
        "pressure_unit": {
            "name": lammps_pressure_unit(potential.units)[0],
            "in_bar": lammps_pressure_unit(potential.units)[1],
        },
    }
    if simulation is not None and simulation.task == "md":
        try:
            probe = (
                replace(simulation, barostat_coupling=None)
                if simulation.ensemble == "npt" and simulation.barostat_coupling == "in-plane"
                else simulation
            )
            plan = plan_md(probe, units=potential.units, options=getattr(engine_spec, "options", {}))
            record = dict(plan["record"])
            record.pop("seed", None)
            record.pop("velocity_seed", None)
            record.pop("thermostat_seed", None)
            if simulation.ensemble == "npt":
                record["barostat_coupling"] = "resolved from the structure's geometry at run time"
            report["integrator"] = record
        except ConfigError as exc:
            report["integrator"] = {"refused": str(exc)}
    return report


def run_job(
    potential: LammpsMlipPotentialSpec,
    atoms,
    simulation: SimulationSpec,
    *,
    engine_spec,
    workdir: Path,
    launch: LammpsLaunch | None = None,
) -> dict:
    """Stage, render, run and convert one job; return a canonical payload.

    The payload carries energy/forces/stress/per-atom energies in canonical
    units and the source basis, the final geometry, and -- for MD -- the
    verified trajectory (converted to ``smoke_md.extxyz`` in the source
    basis), the logged energy series and the resolved integrator.
    """
    import numpy as np

    from ..diagnostics import temperature_ndof
    from ..specs import resolve_seed
    from ..structures import constraint_summary, fixed_atom_indices

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    launch = launch or resolve_launch(engine_spec)
    md = simulation.task == "md"
    geometry = prepare_geometry(atoms)
    fixed = fixed_atom_indices(atoms) if md else ()
    vacuum = vacuum_axes(atoms, simulation) if md and simulation.ensemble == "npt" else ()
    seed = resolve_seed(simulation.seed) if md else None
    options = dict(getattr(engine_spec, "options", {}) or {})

    staged, model_records = stage_model_files(potential, workdir)
    write_data_file(atoms, staged, workdir / DATA_FILE, geometry=geometry)
    commands = render_deck(
        staged,
        simulation,
        n_types=len(staged.type_map),
        pbc=geometry.pbc,
        fixed_ids=[i + 1 for i in fixed],
        n_atoms=len(atoms),
        seed=seed or 1,
        vacuum=vacuum,
        options=options,
        triclinic=geometry.triclinic,
    )
    want_per_atom = bool(simulation.compute_per_atom_energy)
    want_stress = wants_stress(simulation)
    raw = run_deck(
        commands,
        workdir=workdir,
        launch=launch,
        n_atoms=len(atoms),
        want_per_atom=want_per_atom,
        want_stress=want_stress,
    )

    units = potential.units
    system = lammps_unit_system(units)
    rotation = np.array(geometry.rotation)
    thermo = raw["thermo"]
    energy = float(thermo["pe"]) * system.energy_in_eV
    forces = (np.asarray(raw["forces"], dtype=float) @ rotation.T) * system.force_in_eV_per_angstrom
    lo = np.array(raw["box"]["lo"], dtype=float)
    lammps_cell = np.array(raw["box"]["lammps_cell"], dtype=float)
    positions = (np.asarray(raw["positions"], dtype=float) - lo) @ rotation.T
    final_cell = lammps_cell @ rotation.T

    per_atom = None
    if want_per_atom:
        per_atom_native = np.asarray(raw["per_atom_energy"], dtype=float)
        total = float(per_atom_native.sum())
        if abs(total - float(thermo["pe"])) > 1e-6 * max(1.0, abs(float(thermo["pe"]))):
            raise ResultError(
                f"the per-atom energies sum to {total!r} but LAMMPS reports pe = "
                f"{thermo['pe']!r}; the pair style does not tally a complete per-atom "
                "energy, so none is returned"
            )
        per_atom = (per_atom_native * system.energy_in_eV).tolist()

    stress = None
    if want_stress:
        xx, yy, zz, xy, xz, yz = (float(v) for v in raw["virial"])
        tensor = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        stress = list(lammps_pressure_tensor_to_stress(rotation @ tensor @ rotation.T, units))

    native = {
        "route": raw["route"],
        "launch": raw["launch"],
        "build": raw.get("build"),
        "units": system.name,
        "thermo": thermo,
        "deck_path": raw["deck_path"],
        "deck": list(commands),
        "pair_commands": list(staged.render_pair_commands()),
        "geometry": geometry.as_dict(),
        "staged_model_files": model_records,
        "steps_completed": raw["steps_completed"],
        "stress_convention": (
            "virial-only pressure tensor (compute pressure NULL virial), rotated to the "
            "source basis, stress = -pressure" if want_stress else None
        ),
    }
    payload = {
        "energy_eV": energy,
        "forces_eV_per_A": forces.tolist(),
        "stress_eV_per_A3": stress,
        "per_atom_energy_eV": per_atom,
        "symbols": list(atoms.get_chemical_symbols()),
        "wall_time_s": raw.get("wall_time_s"),
        "trajectory_path": None,
        "log_path": raw.get("log_path"),
        "native": native,
    }
    if not md:
        error = float(np.abs(positions - atoms.get_positions()).max())
        native["position_roundtrip_max_error_angstrom"] = error
        if error > _POSITION_ROUNDTRIP_ANGSTROM:
            raise ResultError(
                f"positions read back from LAMMPS differ from the input by {error:.3g} "
                "Angstrom after rotating back; the structure was not passed through intact"
            )
        return payload

    # --- MD: trajectory, series, final geometry, integrator --------------
    steps = raw["steps_completed"]
    if steps != simulation.steps:
        raise ResultError(
            f"LAMMPS completed {steps} steps; the job asked for {simulation.steps}"
        )
    traj_every = max(1, simulation.trajectory_interval)
    frames = read_dump(
        workdir / TRAJECTORY_FILE,
        n_atoms=len(atoms),
        expected_timesteps=expected_steps(steps, traj_every),
    )
    extxyz = write_source_frame_trajectory(
        frames, atoms, geometry, workdir / TRAJECTORY_EXTXYZ, units=units
    )
    series = read_series(workdir / SERIES_FILE, steps=steps, interval=max(1, simulation.log_interval))
    plan = plan_md(
        simulation, units=units, pbc=geometry.pbc, vacuum=vacuum, n_fixed=len(fixed),
        seed=seed, options=options, triclinic=geometry.triclinic,
    )
    record = dict(plan["record"])
    columns = {key: [row[index] for row in series] for index, key in enumerate(SERIES_KEYS)}
    to_eV = system.energy_in_eV
    energy_series = {
        "time_fs": [step * simulation.timestep_fs for step in columns["step"]],
        "total_energy_eV": [value * to_eV for value in columns["etotal"]],
    }
    conserved = record.get("conserved_quantity")
    if conserved == "econserve":
        energy_series["conserved_energy_eV"] = [value * to_eV for value in columns["econserve"]]
        energy_series["conserved_quantity"] = "LAMMPS econserve (pe + ke + ecouple)"
    n_fixed = len(fixed)
    ndof = temperature_ndof(len(atoms), n_fixed=n_fixed, com_removed=not n_fixed)
    record["temperature_ndof"] = ndof
    constraints = constraint_summary(atoms)
    payload.update(
        {
            "steps_completed": steps,
            "frames_written": len(frames),
            "trajectory_path": str(extxyz),
            "final_positions": positions.tolist(),
            "final_cell": final_cell.tolist(),
            "final_pbc": list(geometry.pbc),
            "integrator_resolved": record,
            "integrator": record["description"],
            "energy_series": energy_series,
            "temperature_start_K": columns["temp"][0],
            "temperature_end_K": columns["temp"][-1],
            "max_temperature_K": max(columns["temp"]),
            "temperature_ndof": ndof,
            "constraints": (
                {
                    "fixed_atoms": constraints["fixed_atoms"],
                    "n_fixed": constraints["n_fixed"],
                    "method": (
                        f"LAMMPS groups: {FROZEN_GROUP} (FixAtoms) excluded from velocity "
                        "creation, integration, thermostat and the temperature compute; "
                        "forces on frozen atoms are reported unmodified"
                    ),
                }
                if n_fixed
                else None
            ),
        }
    )
    native.update(
        raw_trajectory_path=str(workdir / TRAJECTORY_FILE),
        series_path=str(workdir / SERIES_FILE),
        seed=seed,
        seed_drawn=simulation.seed is None,
        final_box_lammps={"lo": lo.tolist(), "cell": lammps_cell.tolist()},
    )
    return payload


def write_source_frame_trajectory(frames, atoms, geometry: LammpsGeometry, path: Path, *, units: str):
    """Convert validated dump frames to extxyz in the source basis.

    Each frame's cell and unwrapped positions are rotated back with the
    run's fixed rotation (the cell may change under NPT) and positions are
    taken relative to the box corner, so fractional coordinates are preserved.
    """
    import numpy as np
    from ase import Atoms
    from ase.io import write

    rotation = np.array(geometry.rotation)
    length = lammps_unit_system(units).length_in_angstrom
    images = []
    for frame in frames:
        lo = np.array(frame.box["lo"], dtype=float)
        cell = np.array(frame.box["lammps_cell"], dtype=float) @ rotation.T * length
        positions = (frame.columns3("xu", "yu", "zu") - lo) @ rotation.T * length
        image = Atoms(
            symbols=atoms.get_chemical_symbols(), positions=positions, cell=cell, pbc=geometry.pbc
        )
        image.info["timestep"] = frame.timestep
        images.append(image)
    path = Path(path)
    write(str(path), images, format="extxyz")
    return path


__all__ = [
    "DATA_FILE",
    "DECK_FILE",
    "DUMP_FILE",
    "THERMO_FILE",
    "SERIES_FILE",
    "TRAJECTORY_FILE",
    "TRAJECTORY_EXTXYZ",
    "LOG_FILE",
    "TIMESTEP_PER_FS",
    "FS_PER_TIME_UNIT",
    "VERIFIED_UNIT_STYLES",
    "THERMO_KEYS",
    "THERMOSTATS",
    "BAROSTATS",
    "KOKKOS_GPU_ARGS",
    "LammpsEngine",
    "LammpsLaunch",
    "LammpsGeometry",
    "DumpFrame",
    "resolve_launch",
    "launch_availability",
    "probe_build",
    "parse_help",
    "boundary_string",
    "prepare_geometry",
    "write_data_file",
    "staged_model_names",
    "rewrite_model_tokens",
    "stage_model_files",
    "plan_md",
    "render_deck",
    "run_deck",
    "read_thermo",
    "read_dump",
    "read_series",
    "expected_steps",
    "check_units",
    "check_request",
    "describe_request",
    "run_job",
]
