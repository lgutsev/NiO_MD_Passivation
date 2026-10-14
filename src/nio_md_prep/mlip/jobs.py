"""The four MLIP operations, as functions the CLI is a thin wrapper around.

``inspect``
    What this machine can do: which cells of the compatibility matrix are
    reachable, what is installed, what the model file says about itself.
    Reads nothing heavier than distribution metadata unless a model is named.

``validate``
    Resolve the bridge, negotiate capabilities, check element coverage and
    render the engine parameters -- without executing anything. This is the
    command that must catch an impossible job before a batch script exists.

``singlepoint``
    Evaluate one geometry and write a manifest.

``smoke-md``
    Run a deliberately tiny trajectory to prove a route integrates. Not a
    replacement for the existing production deposition and relaxation
    workflows, and not intended to become one in this release.

Every operation that executes anything writes ``mlip_manifest.json`` first --
before the numbers exist -- so a crashed job still leaves behind a record of
what was attempted.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from pathlib import Path

from .capabilities import check_elements
from .environment import module_available
from .errors import ConfigError, MlipError
from .provenance import (
    build_manifest,
    device_report,
    package_versions,
    write_manifest,
)
from .registry import REGISTRY, build_bridge, compatibility_matrix, resolve_bridge
from .results import compare_results
from .specs import JobSpec, StructureSpec, as_singlepoint
from .structures import describe_structure, load_structure
from .units import describe_conventions

#: A smoke test is meant to take seconds. Anything longer is a production run
#: wearing a diagnostic's name, so it is refused rather than quietly accepted.
SMOKE_MD_MAX_STEPS = 500


def inspect_environment(job: JobSpec | None = None) -> dict:
    """Report what this machine can run, and what a named model contains."""
    routes = []
    for registration in REGISTRY.entries():
        row = registration.as_dict()
        if registration.supported:
            missing = [m for m in registration.requires if not module_available(m)]
            row["modules_available"] = not missing
            row["missing_modules"] = missing
        else:
            row["modules_available"] = False
            row["missing_modules"] = []
        routes.append(row)

    report = {
        "matrix": compatibility_matrix(),
        "routes": routes,
        "packages": package_versions(),
        "device": device_report(
            getattr(job.potential, "device", None) if job else None
        ),
        "energy_conventions": describe_conventions(),
        "canonical_units": {
            "energy": "eV",
            "force": "eV/Angstrom",
            "length": "Angstrom",
            "stress": "eV/Angstrom^3",
        },
    }
    if job is not None:
        from .potentials import build_adapter

        adapter = build_adapter(job.potential)
        report["potential"] = adapter.inspect()
        report["job"] = job.as_dict()
        registration = resolve_bridge(
            job.potential.kind,
            job.engine.kind,
            preferences=tuple(getattr(job.potential, "implementation", ()) or ()),
        )
        bridge = registration.instantiate(job.potential, job.engine)
        report["selected_route"] = registration.as_dict()
        report["route_availability"] = bridge.availability().as_dict()
        report["capabilities"] = bridge.capabilities().as_dict()
    return report


def validate_job(job: JobSpec, *, structure: StructureSpec | None = None) -> dict:
    """Resolve, negotiate and render -- but execute nothing.

    Deliberately ordered so the most informative failure comes first: an
    impossible combination, then a missing element, then an unmet capability.
    """
    registration = resolve_bridge(
        job.potential.kind,
        job.engine.kind,
        preferences=tuple(getattr(job.potential, "implementation", ()) or ()),
    )
    bridge = registration.instantiate(job.potential, job.engine)
    capabilities = bridge.capabilities()

    atoms = None
    structure_report = None
    structure_spec = structure or job.structure
    if structure_spec is not None:
        atoms = load_structure(structure_spec)
        structure_report = describe_structure(atoms)
        check_elements(
            capabilities, atoms.get_chemical_symbols(), label=job.potential.label
        )

    requirements = bridge.requirements(job.simulation, atoms)
    unmet = requirements.unmet(capabilities)
    availability = bridge.availability()

    report = {
        "ok": not unmet,
        "job": job.as_dict(),
        "route": registration.as_dict(),
        "capabilities": capabilities.as_dict(),
        "requirements": requirements.as_dict(),
        "unmet": unmet,
        "availability": availability.as_dict(),
        "structure": structure_report,
        "engine_parameters": bridge.engine_parameters(),
    }
    if unmet:
        from .errors import CapabilityError

        raise CapabilityError(
            f"{bridge.label} cannot satisfy this job:", unmet
        )
    return report


def run_singlepoint(
    job: JobSpec,
    *,
    output_dir: Path,
    structure: StructureSpec | None = None,
) -> dict:
    """Evaluate one geometry through the resolved route, with provenance."""
    structure_spec = structure or job.structure
    if structure_spec is None:
        raise ConfigError(
            "a single-point evaluation needs a structure: add a [structure] section "
            "to the configuration, or pass --structure"
        )
    simulation = as_singlepoint(job.simulation)
    job = replace(job, simulation=simulation)
    return _execute(job, structure_spec, output_dir=Path(output_dir), md=False)


def run_smoke_md(
    job: JobSpec,
    *,
    output_dir: Path,
    structure: StructureSpec | None = None,
    max_steps: int = SMOKE_MD_MAX_STEPS,
) -> dict:
    """Run a tiny diagnostic trajectory. Refuses to be a production run."""
    structure_spec = structure or job.structure
    if structure_spec is None:
        raise ConfigError(
            "a smoke test needs a structure: add a [structure] section to the "
            "configuration, or pass --structure"
        )
    if job.simulation.task != "md":
        raise ConfigError(
            f"smoke-md needs simulation.task = 'md'; this configuration says "
            f"{job.simulation.task!r}"
        )
    if job.simulation.steps > max_steps:
        raise ConfigError(
            f"smoke-md runs a diagnostic trajectory of at most {max_steps} steps; this "
            f"configuration asks for {job.simulation.steps}. This command is not a "
            "replacement for the existing production deposition and relaxation "
            "workflows -- run those through the classical pipeline, or lower "
            "simulation.steps for a diagnostic."
        )
    return _execute(job, structure_spec, output_dir=Path(output_dir), md=True)


def _execute(job: JobSpec, structure_spec: StructureSpec, *, output_dir: Path, md: bool) -> dict:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    registration = resolve_bridge(
        job.potential.kind,
        job.engine.kind,
        preferences=tuple(getattr(job.potential, "implementation", ()) or ()),
    )
    engine = replace(
        job.engine,
        options={**dict(job.engine.options), "workdir": str(output_dir)},
    )
    bridge = registration.instantiate(job.potential, engine)

    atoms = load_structure(structure_spec)
    structure_report = describe_structure(atoms)
    capabilities = bridge.validate(job.simulation, atoms)
    bridge.require_available()

    # Written before the run, so a crash still leaves a record of the attempt.
    manifest = build_manifest(
        job=job,
        registration=registration,
        capabilities=capabilities,
        requirements=bridge.requirements(job.simulation, atoms),
        structure={**structure_report, "path": str(structure_spec.path)},
        engine_parameters=bridge.engine_parameters(),
    )
    manifest_path = write_manifest(output_dir, manifest)

    if md:
        trajectory = bridge.run_md(atoms, job.simulation, workdir=output_dir)
        results = {"trajectory": trajectory.as_dict()}
        payload = {"trajectory": trajectory}
    else:
        result = bridge.singlepoint(atoms, job.simulation)
        results = {"singlepoint": result.as_dict(include_arrays=False)}
        payload = {"result": result}

    manifest["results"] = results
    write_manifest(output_dir, manifest)
    return {
        "manifest_path": str(manifest_path),
        "output_dir": str(output_dir),
        "route": registration.as_dict(),
        "structure": structure_report,
        "results": results,
        **payload,
    }


def compare_engines(
    job: JobSpec,
    engines: Iterable[str],
    *,
    output_dir: Path,
    structure: StructureSpec | None = None,
    convention: str | None = None,
) -> dict:
    """Evaluate one geometry through several engines and report the differences.

    The scientific acceptance test for this subsystem. The first engine listed
    is the reference (ASE, for MACE). Every result is harmonised to a single
    energy convention before the energies are subtracted -- comparing OpenMM's
    interaction energy with ASE's total energy would otherwise produce a
    difference of thousands of eV while the forces agree perfectly.

    No tolerance is applied here. The caller decides what counts as passing,
    against a documented reference established from the actual adapters.
    """
    structure_spec = structure or job.structure
    if structure_spec is None:
        raise ConfigError("a cross-engine comparison needs a structure")
    engines = list(engines)
    if len(engines) < 2:
        raise ConfigError("a cross-engine comparison needs at least two engines")

    output_dir = Path(output_dir)
    atoms = load_structure(structure_spec)
    simulation = as_singlepoint(job.simulation)

    results = {}
    parameters = {}
    reference_energies = None
    for kind in engines:
        engine_spec = replace(
            job.engine,
            kind=kind,
            options={**dict(job.engine.options), "workdir": str(output_dir / kind)},
        )
        bridge = build_bridge(job.potential, engine_spec)
        bridge.require_available()
        results[kind] = bridge.singlepoint(atoms, simulation)
        parameters[kind] = bridge.engine_parameters()
        if reference_energies is None:
            reference_energies = bridge.atomic_reference_energies()

    reference_kind = engines[0]
    target = convention or results[reference_kind].energy_convention
    comparisons = {}
    for kind in engines[1:]:
        comparison = compare_results(
            results[reference_kind],
            results[kind],
            atomic_reference_energies=reference_energies,
            convention=target,
        )
        comparisons[f"{reference_kind}->{kind}"] = comparison.as_dict()

    report = {
        "structure": describe_structure(atoms),
        "reference_engine": reference_kind,
        "energy_convention": target,
        "atomic_reference_energies_known": bool(reference_energies),
        "results": {k: v.as_dict(include_arrays=False) for k, v in results.items()},
        "engine_parameters": parameters,
        "comparisons": comparisons,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    import json

    (output_dir / "cross_engine_comparison.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8"
    )
    return report


__all__ = [
    "SMOKE_MD_MAX_STEPS",
    "inspect_environment",
    "validate_job",
    "run_singlepoint",
    "run_smoke_md",
    "compare_engines",
    "MlipError",
]
