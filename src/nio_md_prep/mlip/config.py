"""Parsing and validating MLIP job configurations.

Deliberately independent of execution. Parsing a configuration imports
``tomllib`` and this package's spec classes and nothing else, so a
configuration destined for a GPU node can be checked on a laptop, in CI, or
in a pull-request review, long before anything is installed.

A minimal MACE-on-ASE configuration::

    [potential]
    kind = "mace"
    model = "/path/to/nio_phosphonate.model"
    device = "cuda"
    precision = "float32"

    [engine]
    kind = "ase"

    [simulation]
    task = "md"
    ensemble = "nvt"
    temperature_K = 400
    timestep_fs = 0.5
    steps = 10000
    seed = 12345

and the LAMMPS-native alternative::

    [potential]
    kind = "lammps"
    pair_style = "mliap unified mace_nio.pt 0"
    pair_coeff = ["* * Ni O P C H"]
    units = "metal"
    model_files = ["mace_nio.pt"]

    [potential.type_map]
    1 = "Ni"
    2 = "O"
    3 = "P"
    4 = "C"
    5 = "H"

    [engine]
    kind = "lammps"

Unknown keys are errors, not warnings: a typo in ``temperature_K`` must not
quietly produce a 0 K trajectory.
"""
from __future__ import annotations

import tomllib
from collections.abc import Mapping
from pathlib import Path

from .errors import ConfigError
from .specs import (
    POTENTIAL_SPECS,
    EngineSpec,
    JobSpec,
    LammpsMlipPotentialSpec,
    MacePotentialSpec,
    MockPotentialSpec,
    PotentialSpec,
    SimulationSpec,
    StructureSpec,
)

TOP_LEVEL_SECTIONS = ("potential", "engine", "simulation", "structure", "job")

#: Configuration key -> spec field, per potential kind. Keys absent from the
#: mapping are rejected, which is what makes a typo fail loudly.
_MACE_KEYS = {
    "kind": None,
    "label": "label",
    "model": "model_path",
    "model_path": "model_path",
    "sha256": "model_sha256",
    "model_sha256": "model_sha256",
    "elements": "declared_elements",
    "device": "device",
    "precision": "precision",
    "cutoff_angstrom": "cutoff_angstrom",
    "implementation": "implementation",
    "model_format": "model_format",
    "compile_mode": "compile_mode",
    "energy_convention": "energy_convention",
    "atomic_reference_energies": "atomic_reference_energies",
    "notes": "notes",
}

_LAMMPS_KEYS = {
    "kind": None,
    "label": "label",
    "pair_style": "pair_style",
    "pair_coeff": "pair_coeff",
    "type_map": "type_map",
    "required_packages": "required_packages",
    "units": "units",
    "atom_style": "atom_style",
    "model_files": "model_paths",
    "model_paths": "model_paths",
    "model_hashes": "model_hashes",
    "framework": "framework",
    "extra_commands": "extra_commands",
    "newton": "newton",
    "energy_convention": "energy_convention",
    "atomic_reference_energies": "atomic_reference_energies",
    "notes": "notes",
}

_MOCK_KEYS = {
    "kind": None,
    "label": "label",
    "elements": "declared_elements",
    "epsilon_eV": "epsilon_eV",
    "sigma_angstrom": "sigma_angstrom",
    "cutoff_angstrom": "cutoff_angstrom",
    "precision": "precision",
    "energy_convention": "energy_convention",
    "atomic_reference_energies": "atomic_reference_energies",
    "notes": "notes",
}

POTENTIAL_KEYS = {
    "mace": _MACE_KEYS,
    "lammps": _LAMMPS_KEYS,
    "mock": _MOCK_KEYS,
}

ENGINE_KEYS = {
    "kind": "kind",
    "platform": "platform",
    "precision": "precision",
    "threads": "threads",
    "executable": "executable",
    "options": "options",
}

SIMULATION_KEYS = {name: name for name in SimulationSpec.__dataclass_fields__}

STRUCTURE_KEYS = {"path": "path", "format": "format", "index": "index", "label": "label"}

#: Path-valued keys, resolved relative to the configuration file's directory
#: so a config can be moved with its model without becoming machine-specific.
_PATH_FIELDS = {"model_path", "path"}
_PATH_LIST_FIELDS = {"model_paths"}


def load_config(path: Path) -> dict:
    """Read a TOML configuration file. No validation, no imports beyond tomllib."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MLIP configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def parse_job(path: Path) -> JobSpec:
    """Parse and fully validate a configuration file into a :class:`JobSpec`."""
    path = Path(path)
    return parse_job_mapping(load_config(path), base_dir=path.parent, source=path)


def parse_job_mapping(
    payload: Mapping,
    *,
    base_dir: Path | None = None,
    source: Path | None = None,
) -> JobSpec:
    """Parse an already-loaded mapping. Used by tests and by programmatic callers."""
    _reject_unknown(payload, TOP_LEVEL_SECTIONS, "configuration")
    if "potential" not in payload:
        raise ConfigError("configuration is missing the [potential] section")
    if "engine" not in payload:
        raise ConfigError("configuration is missing the [engine] section")

    potential = parse_potential(payload["potential"], base_dir=base_dir)
    engine = parse_engine(payload["engine"])
    simulation = parse_simulation(payload.get("simulation", {}))
    structure = (
        parse_structure(payload["structure"], base_dir=base_dir)
        if "structure" in payload
        else None
    )
    job = payload.get("job", {})
    _reject_unknown(job, ("name",), "[job]")
    return JobSpec(
        potential=potential,
        engine=engine,
        simulation=simulation,
        structure=structure,
        name=str(job.get("name", source.stem if source else "mlip-job")),
        source=source,
    )


def parse_potential(payload: Mapping, *, base_dir: Path | None = None) -> PotentialSpec:
    """Build the right :class:`PotentialSpec` subclass for ``potential.kind``."""
    if "kind" not in payload:
        raise ConfigError(
            "[potential] is missing 'kind'; expected one of "
            f"{', '.join(sorted(POTENTIAL_SPECS))}"
        )
    kind = str(payload["kind"])
    if kind not in POTENTIAL_SPECS:
        from .errors import UnknownPotentialError

        raise UnknownPotentialError(
            f"unknown potential.kind {kind!r}; known kinds: "
            f"{', '.join(sorted(POTENTIAL_SPECS))}"
        )
    keys = POTENTIAL_KEYS[kind]
    _reject_unknown(payload, keys, f"[potential] (kind = {kind!r})")
    kwargs = _map_keys(payload, keys, base_dir=base_dir)
    if kind == "mace":
        kwargs.setdefault("label", Path(kwargs.get("model_path", "mace")).stem)
        if "declared_elements" not in kwargs:
            raise ConfigError(
                "[potential] for a MACE model must declare 'elements': the element "
                "coverage is validated against every structure before anything runs, "
                "and guessing it from the model file at parse time would require "
                "loading torch"
            )
        return MacePotentialSpec(**kwargs)
    if kind == "lammps":
        kwargs.setdefault("label", str(kwargs.get("framework") or "lammps-mlip"))
        return LammpsMlipPotentialSpec(**kwargs)
    kwargs.setdefault("label", "mock")
    return MockPotentialSpec(**kwargs)


def parse_engine(payload: Mapping) -> EngineSpec:
    _reject_unknown(payload, ENGINE_KEYS, "[engine]")
    if "kind" not in payload:
        raise ConfigError("[engine] is missing 'kind'")
    return EngineSpec(**_map_keys(payload, ENGINE_KEYS))


def parse_simulation(payload: Mapping) -> SimulationSpec:
    _reject_unknown(payload, SIMULATION_KEYS, "[simulation]")
    return SimulationSpec(**_map_keys(payload, SIMULATION_KEYS))


def parse_structure(payload: Mapping, *, base_dir: Path | None = None) -> StructureSpec:
    _reject_unknown(payload, STRUCTURE_KEYS, "[structure]")
    if "path" not in payload:
        raise ConfigError("[structure] is missing 'path'")
    return StructureSpec(**_map_keys(payload, STRUCTURE_KEYS, base_dir=base_dir))


def _reject_unknown(payload: Mapping, allowed, section: str) -> None:
    if not isinstance(payload, Mapping):
        raise ConfigError(f"{section} must be a table")
    unknown = sorted(set(payload) - set(allowed))
    if unknown:
        raise ConfigError(
            f"{section} has unknown key(s): {', '.join(unknown)}. "
            f"Known keys: {', '.join(sorted(allowed))}"
        )


def _map_keys(payload: Mapping, keys: Mapping, *, base_dir: Path | None = None) -> dict:
    kwargs: dict = {}
    for key, value in payload.items():
        field = keys[key]
        if field is None:
            continue
        if field in _PATH_FIELDS:
            value = _resolve(value, base_dir)
        elif field in _PATH_LIST_FIELDS:
            value = tuple(_resolve(item, base_dir) for item in value)
        elif field in ("implementation", "declared_elements", "pair_coeff",
                       "required_packages", "extra_commands"):
            value = tuple(value) if not isinstance(value, str) else (value,)
        elif field == "type_map":
            value = {int(k): str(v) for k, v in dict(value).items()}
        if field in kwargs and kwargs[field] != value:
            raise ConfigError(
                f"conflicting aliases for {field!r} in the same section; use one spelling"
            )
        kwargs[field] = value
    return kwargs


def _resolve(value, base_dir: Path | None) -> Path:
    path = Path(value)
    if base_dir is not None and not path.is_absolute():
        return (Path(base_dir) / path).resolve()
    return path


def describe_config(job: JobSpec) -> dict:
    """A JSON-serialisable echo of the parsed configuration, for reports."""
    return job.as_dict()


__all__ = [
    "TOP_LEVEL_SECTIONS",
    "POTENTIAL_KEYS",
    "ENGINE_KEYS",
    "SIMULATION_KEYS",
    "STRUCTURE_KEYS",
    "load_config",
    "parse_job",
    "parse_job_mapping",
    "parse_potential",
    "parse_engine",
    "parse_simulation",
    "parse_structure",
    "describe_config",
]
