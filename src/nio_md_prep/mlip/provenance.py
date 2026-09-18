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
:func:`device_report` only touches torch if torch is already importable.
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
    "openmm-ml",
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
    """Refuse to run when a model file is not the one the configuration names."""
    for path, expected in declared.items():
        observed = actual.get(str(path))
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


def device_report(requested_device: str | None = None) -> dict:
    """Describe the compute device, touching torch only if it is already there.

    ``find_spec`` avoids importing torch on a machine that does not have it,
    which keeps ``mlip inspect`` fast and keeps a plain ``nio-md-prep``
    install free of a CUDA dependency.
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
    if find_spec("torch") is None:
        return report
    try:
        import torch
    except Exception as exc:  # pragma: no cover - a broken torch install
        report["torch"] = {"error": f"{type(exc).__name__}: {exc}"}
        return report
    info: dict = {
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
            "device": device_report(getattr(potential, "device", None)),
        },
        "engine_parameters": dict(engine_parameters or {}),
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
    """Thermostat/barostat settings, pulled out where a reader will look."""
    return {
        "task": simulation.task,
        "ensemble": simulation.ensemble,
        "thermostat": simulation.thermostat,
        "barostat": simulation.barostat,
        "thermostat_damping_fs": simulation.thermostat_damping_fs,
        "barostat_damping_fs": simulation.barostat_damping_fs,
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
    """Write ``mlip_manifest.json`` into a job directory and return its path."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / MANIFEST_NAME
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=False, default=str) + "\n",
        encoding="utf-8",
    )
    return path


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
    "TRACKED_DISTRIBUTIONS",
    "sha256_file",
    "hash_model_files",
    "verify_hashes",
    "git_commit",
    "package_versions",
    "device_report",
    "build_manifest",
    "write_manifest",
    "read_manifest",
]
