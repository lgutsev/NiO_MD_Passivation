"""``mlip_manifest.json``: what was run, with what, and under what conventions.

Provenance comes before long MD, not after it. Once several MACE committee
members exist -- and later DeepMD or DPA models -- a directory full of
trajectories is worthless unless each one records exactly which model
produced it. So every MLIP job directory gets a machine-readable manifest
containing the model SHA256, the input structure hash, the git commit, the
potential implementation, the simulation engine, package versions,
CUDA/device information, the element mapping, precision, energy convention,
thermostat/integrator settings, the seed, and the generated engine
parameters.

For LAMMPS that last item is literal: the exact rendered ``pair_style`` and
``pair_coeff`` command strings are preserved, not re-derivable instructions
for producing them.

Hashing and version probing are done without importing any MLIP backend:
:func:`package_versions` reads distribution metadata, and
:func:`device_report` reads torch's CUDA state only when torch is already
imported in this process or the caller explicitly asks for the probe.

A manifest carries a ``status`` -- ``prepared`` (assembled and
hash-verified, nothing started), ``running`` (the engine has been invoked),
``completed`` or ``failed`` -- and an ``error`` record for a failure, and is
always replaced atomically (:func:`write_manifest`), so a crash mid-write can
never leave a truncated file behind.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from importlib import metadata
from importlib.util import find_spec
from pathlib import Path

from .errors import ModelIntegrityError, ProvenanceError

MANIFEST_NAME = "mlip_manifest.json"
MANIFEST_VERSION = 1

#: The lifecycle of a job directory, in order. ``failed`` can follow either
#: ``prepared`` or ``running``.
STATUSES = ("prepared", "running", "completed", "failed")

#: Distributions worth recording whenever they are installed. Absent ones are
#: reported as ``null`` rather than omitted, so a manifest shows what was
#: *not* present as clearly as what was.
TRACKED_DISTRIBUTIONS = (
    "nio-md-prep",
    "ase",
    "numpy",
    "torch",
    "mace-torch",
    "openmm",
    # The OpenMM-ML distribution is published as "openmmml" (setup.py
    # name='openmmml'); "openmm-ml" is not a distribution and always read None.
    "openmmml",
    "openmm-torch",
    "nnpops",
    "lammps",
    "deepmd-kit",
)


def sha256_file(path: Path, *, chunk: int = 1 << 20) -> str:
    """Stream a file's SHA256. Model files are large; do not read them whole."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"cannot hash a file that does not exist: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_model_files(paths: Iterable[Path]) -> dict[str, str]:
    return {str(Path(p)): sha256_file(p) for p in paths}


def verify_hashes(declared: Mapping[str, str], actual: Mapping[str, str]) -> None:
    """Refuse to run when a model file is not the one the configuration names.

    Paths are matched on :func:`~nio_md_prep.mlip.specs.path_key` (absolute,
    normalised), so ``nio.pb`` and ``/config/dir/nio.pb`` are the same file.
    """
    from .specs import path_key

    observed_by_key = {path_key(p): sha for p, sha in actual.items()}
    for path, expected in declared.items():
        observed = observed_by_key.get(path_key(path))
        if observed is None:
            raise ProvenanceError(
                f"configuration declares a sha256 for {path}, but that file was not "
                "hashed; the path in the configuration may be wrong"
            )
        if observed.lower() != expected.lower():
            raise ModelIntegrityError(path, expected, observed)


def git_commit(root: Path | None = None) -> dict:
    """The repository revision, if this is a checkout. Never fatal."""
    root = Path(root) if root else Path(__file__).resolve().parents[3]
    try:
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "dirty": None, "root": str(root)}
    if commit.returncode != 0:
        return {"commit": None, "dirty": None, "root": str(root)}
    return {
        "commit": commit.stdout.strip() or None,
        "dirty": bool(dirty.stdout.strip()) if dirty.returncode == 0 else None,
        "root": str(root),
    }


def package_versions(names: Iterable[str] = TRACKED_DISTRIBUTIONS) -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def device_report(requested_device: str | None = None, *, probe_torch: bool = False) -> dict:
    """Describe the compute device without importing torch behind the caller's back.

    torch's CUDA state is read only when torch is *already imported* in this
    process (a MACE route has loaded it, so reading it costs nothing), or
    when ``probe_torch=True`` explicitly asks for the import (``mlip inspect
    --probe-torch``; a GPU MACE job's manifest). Otherwise an installed torch
    is reported by version from distribution metadata with ``probed: False``,
    and an absent one as ``None``. Importing torch takes seconds and
    initialises CUDA, which a classical or ``validate`` path must never do.
    """
    report: dict = {
        "requested": requested_device,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
        "hostname": platform.node(),
        "cpu_count": os.cpu_count(),
        "torch": None,
    }
    torch = sys.modules.get("torch")
    if torch is None:
        if find_spec("torch") is None:
            return report
        if not probe_torch:
            try:
                version = metadata.version("torch")
            except metadata.PackageNotFoundError:  # pragma: no cover - odd install
                version = None
            report["torch"] = {
                "version": version,
                "probed": False,
                "note": "torch is installed but was not imported, so CUDA was not probed",
            }
            return report
        try:
            import torch
        except Exception as exc:  # pragma: no cover - a broken torch install
            report["torch"] = {"error": f"{type(exc).__name__}: {exc}"}
            return report
    info: dict = {
        "probed": True,
        "version": getattr(torch, "__version__", None),
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_version": getattr(getattr(torch, "version", None), "cuda", None),
        "device_count": 0,
        "devices": [],
    }
    if info["cuda_available"]:
        info["device_count"] = torch.cuda.device_count()
        for index in range(info["device_count"]):
            try:
                properties = torch.cuda.get_device_properties(index)
                info["devices"].append(
                    {
                        "index": index,
                        "name": properties.name,
                        "total_memory_bytes": properties.total_memory,
                        "capability": f"{properties.major}.{properties.minor}",
                    }
                )
            except Exception:  # pragma: no cover - driver quirks
                info["devices"].append({"index": index, "name": None})
    report["torch"] = info
    return report


def build_manifest(
    *,
    job,
    registration,
    capabilities=None,
    requirements=None,
    structure=None,
    engine_parameters: Mapping | None = None,
    results: Mapping | None = None,
    extra: Mapping | None = None,
    verify: bool = True,
    probe_torch: bool = False,
) -> dict:
    """Assemble the manifest for one MLIP job.

    ``engine_parameters`` is where an engine records what it actually
    generated: for LAMMPS, the rendered command block verbatim; for OpenMM,
    the potential name, energy type and platform; for ASE, the calculator and
    integrator settings.
    """
    potential = job.potential
    declared = potential.declared_hashes()
    actual: dict[str, str] = {}
    missing: list[str] = []
    for path in potential.model_files():
        if Path(path).exists():
            actual[str(path)] = sha256_file(path)
        else:
            missing.append(str(path))
    if verify and declared:
        verify_hashes(declared, actual)

    simulation = job.simulation
    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "kind": "mlip",
        "name": job.name,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config_source": str(job.source) if job.source else None,
        "git": git_commit(),
        "potential": {
            **potential.as_dict(),
            "model_sha256_observed": actual,
            "model_files_missing": missing,
        },
        "engine": job.engine.as_dict(),
        "bridge": {
            "potential": registration.potential,
            "engine": registration.engine,
            "implementation": registration.implementation,
            "summary": registration.summary,
            "status": registration.status,
            "requires": list(registration.requires),
            "packages": list(registration.packages),
            "notes": list(registration.notes),
        },
        "simulation": {
            **simulation.as_dict(),
            "integrator": _integrator_report(simulation),
        },
        "units": {
            "canonical_energy": "eV",
            "canonical_force": "eV/Angstrom",
            "canonical_length": "Angstrom",
            "canonical_stress": "eV/Angstrom^3",
            "engine_native": (
                engine_parameters.get("native_units") if engine_parameters else None
            ),
        },
        "energy_convention": {
            "potential_native": potential.energy_convention,
            "reported": (
                engine_parameters.get("energy_convention") if engine_parameters else None
            ),
            "requested": simulation.energy_convention,
        },
        "element_mapping": _element_mapping(potential),
        "environment": {
            "packages": package_versions(),
            # Refreshed by jobs after the run, when a MACE route has imported
            # torch and its CUDA state can be read without a new import.
            "device": device_report(getattr(potential, "device", None), probe_torch=probe_torch),
        },
        "engine_parameters": dict(engine_parameters or {}),
        "status": "prepared",
        "error": None,
    }
    if structure is not None:
        manifest["structure"] = dict(structure)
    if capabilities is not None:
        manifest["capabilities"] = capabilities.as_dict()
    if requirements is not None:
        manifest["requirements"] = requirements.as_dict()
    if results is not None:
        manifest["results"] = dict(results)
    if extra:
        manifest.update(extra)
    return manifest


def _integrator_report(simulation) -> dict:
    """Thermostat/barostat settings as *requested*, pulled out where a reader will look.

    The damping times are shown after default resolution, since that is what
    every engine applies. What actually integrated (class or fix, resolved
    parameters, drawn seed) is the engine's report, under
    ``results.trajectory.integrator_resolved``.
    """
    md = simulation.task == "md"
    thermostatted = md and simulation.ensemble in ("nvt", "npt")
    return {
        "task": simulation.task,
        "ensemble": simulation.ensemble,
        "thermostat": simulation.thermostat,
        "barostat": simulation.barostat,
        "barostat_coupling": simulation.barostat_coupling,
        "thermostat_damping_fs": simulation.thermostat_damping_fs,
        "barostat_damping_fs": simulation.barostat_damping_fs,
        "resolved_thermostat_damping_fs": (
            simulation.resolved_thermostat_damping_fs if thermostatted else None
        ),
        "resolved_barostat_damping_fs": (
            simulation.resolved_barostat_damping_fs
            if md and simulation.ensemble == "npt"
            else None
        ),
        "vacuum_gap_threshold_angstrom": (
            simulation.resolved_vacuum_gap_threshold_angstrom
            if md and simulation.ensemble == "npt"
            else None
        ),
        "timestep_fs": simulation.timestep_fs,
        "steps": simulation.steps,
        "temperature_K": simulation.temperature_K,
        "pressure_bar": simulation.pressure_bar,
        "seed": simulation.seed,
    }


def _element_mapping(potential) -> dict:
    mapping: dict = {"elements": list(potential.elements)}
    type_map = getattr(potential, "type_map", None)
    if type_map:
        mapping["lammps_type_map"] = {str(k): v for k, v in sorted(type_map.items())}
    if getattr(potential, "atomic_reference_energies", None):
        mapping["atomic_reference_energies_eV"] = dict(potential.atomic_reference_energies)
    return mapping


def write_manifest(directory: Path, manifest: Mapping) -> Path:
    """Atomically write ``mlip_manifest.json`` into a job directory; return its path.

    The JSON is written to a temporary file in the same directory, flushed to
    disk, and moved over the manifest with :func:`os.replace`, which is atomic
    on POSIX and Windows. A reader therefore sees either the previous complete
    manifest or the new one, never a partial file. A ``status`` outside
    :data:`STATUSES` is refused.
    """
    status = manifest.get("status")
    if status is not None and status not in STATUSES:
        raise ProvenanceError(
            f"manifest status must be one of {', '.join(STATUSES)}; got {status!r}"
        )
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / MANIFEST_NAME
    text = json.dumps(manifest, indent=2, sort_keys=False, default=_json_default) + "\n"
    temporary = directory / f".{MANIFEST_NAME}.{os.getpid()}.tmp"
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def mark_manifest(
    directory: Path, manifest: dict, status: str, *, error: BaseException | None = None
) -> Path:
    """Set ``status`` (and, for a failure, the ``error`` record) and rewrite the manifest.

    ``error`` is recorded as ``{"type", "message"}`` plus ``finished_utc``;
    ``completed`` and ``failed`` also stamp ``finished_utc``.
    """
    if status not in STATUSES:
        raise ProvenanceError(f"unknown manifest status {status!r}")
    manifest["status"] = status
    if status in ("completed", "failed"):
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    if error is not None:
        manifest["error"] = {"type": type(error).__name__, "message": str(error)}
    return write_manifest(directory, manifest)


def _json_default(value):
    """Serialise the few non-JSON types a manifest legitimately carries."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (set, frozenset)):
        return sorted(value, key=str)
    if hasattr(value, "tolist"):  # numpy arrays and scalars
        return value.tolist()
    if hasattr(value, "as_dict"):
        return value.as_dict()
    return str(value)


def read_manifest(directory: Path) -> dict:
    path = Path(directory)
    if path.is_dir():
        path = path / MANIFEST_NAME
    if not path.exists():
        raise FileNotFoundError(f"no MLIP manifest at {path}")
    return json.loads(path.read_text(encoding="utf-8"))


__all__ = [
    "MANIFEST_NAME",
    "MANIFEST_VERSION",
    "STATUSES",
    "TRACKED_DISTRIBUTIONS",
    "sha256_file",
    "hash_model_files",
    "verify_hashes",
    "git_commit",
    "package_versions",
    "device_report",
    "build_manifest",
    "write_manifest",
    "mark_manifest",
    "read_manifest",
]
