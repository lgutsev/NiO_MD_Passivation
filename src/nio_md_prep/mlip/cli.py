"""``nio-md-prep mlip ...`` -- three levels of commitment, in order.

``inspect`` and ``validate`` are the production-ready commands in this
release: they answer "what can this machine run" and "is this job possible"
without executing anything. ``singlepoint`` and ``smoke-md`` execute, but
``smoke-md`` is a diagnostic with a hard step ceiling, not a production MD
driver -- the existing classical deposition and relaxation workflows remain
the way production trajectories are produced.

``compare`` drives the cross-engine single-point equivalence check that is
this subsystem's scientific acceptance test.

Everything here is a thin shell over :mod:`nio_md_prep.mlip.jobs`. Building
the top-level parser imports this module (to register the ``mlip`` group)
and the configuration/spec layer, but never ASE, torch, MACE, OpenMM or
LAMMPS; ``tests/test_mlip_isolation.py`` checks that in a fresh interpreter.
"""
from __future__ import annotations

import json
from pathlib import Path

from .config import parse_job
from .specs import StructureSpec

COMMANDS = ("inspect", "validate", "singlepoint", "smoke-md", "compare")


def add_parser(subparsers) -> None:
    """Attach the ``mlip`` command group to the top-level parser."""
    parser = subparsers.add_parser(
        "mlip",
        help="machine-learned interatomic potential jobs (preliminary subsystem)",
    )
    group = parser.add_subparsers(dest="mlip_command", required=True)

    inspect = group.add_parser(
        "inspect", help="report environment, installed routes and model capabilities"
    )
    inspect.add_argument("config", type=Path, nargs="?")
    inspect.add_argument("--json", action="store_true", help="emit the full report as JSON")
    inspect.add_argument(
        "--probe-torch",
        action="store_true",
        help="import torch to report CUDA devices (slow; off by default)",
    )

    validate = group.add_parser(
        "validate", help="check a configuration and its route without executing anything"
    )
    validate.add_argument("config", type=Path)
    validate.add_argument("--structure", type=Path)
    validate.add_argument("--json", action="store_true")

    singlepoint = group.add_parser(
        "singlepoint", help="evaluate one geometry and write an MLIP manifest"
    )
    singlepoint.add_argument("config", type=Path)
    singlepoint.add_argument("--output", type=Path, required=True)
    singlepoint.add_argument("--structure", type=Path)
    singlepoint.add_argument("--json", action="store_true")

    smoke = group.add_parser(
        "smoke-md", help="run a tiny diagnostic trajectory (not a production MD driver)"
    )
    smoke.add_argument("config", type=Path)
    smoke.add_argument("--output", type=Path, required=True)
    smoke.add_argument("--structure", type=Path)
    smoke.add_argument("--json", action="store_true")

    compare = group.add_parser(
        "compare",
        help="cross-engine single-point equivalence check for one model and structure",
    )
    compare.add_argument("config", type=Path)
    compare.add_argument("--output", type=Path, required=True)
    compare.add_argument(
        "--engines",
        default="ase,lammps,openmm",
        help="comma-separated engines to compare; the first is the reference",
    )
    compare.add_argument("--structure", type=Path)
    compare.add_argument("--energy-convention", choices=("total", "interaction"))
    compare.add_argument("--json", action="store_true")


def run(args) -> int:
    """Dispatch one ``mlip`` subcommand. Errors propagate to the top-level handler."""
    from . import jobs

    structure = (
        StructureSpec(path=args.structure) if getattr(args, "structure", None) else None
    )

    if args.mlip_command == "inspect":
        job = parse_job(args.config) if args.config else None
        report = jobs.inspect_environment(job, probe_torch=args.probe_torch)
        _emit(report, args.json, _format_inspect)
        return 0

    job = parse_job(args.config)

    if args.mlip_command == "validate":
        report = jobs.validate_job(job, structure=structure)
        _emit(report, args.json, _format_validate)
        return 0

    if args.mlip_command == "singlepoint":
        report = jobs.run_singlepoint(job, output_dir=args.output, structure=structure)
        _emit(report, args.json, _format_singlepoint)
        return 0

    if args.mlip_command == "smoke-md":
        report = jobs.run_smoke_md(job, output_dir=args.output, structure=structure)
        _emit(report, args.json, _format_smoke)
        return 0

    engines = [e.strip() for e in args.engines.split(",") if e.strip()]
    report = jobs.compare_engines(
        job,
        engines,
        output_dir=args.output,
        structure=structure,
        convention=args.energy_convention,
    )
    _emit(report, args.json, _format_compare)
    return 0


def _emit(report: dict, as_json: bool, formatter) -> None:
    if as_json:
        print(json.dumps(_serialisable(report), indent=2, default=str))
    else:
        print(formatter(report))


def _serialisable(report: dict) -> dict:
    """Drop the live objects the JSON view has no use for."""
    return {k: v for k, v in report.items() if k not in ("result", "trajectory")}


# ---------------------------------------------------------------------------
# Human-readable formatting
# ---------------------------------------------------------------------------


def _format_inspect(report: dict) -> str:
    lines = [
        "Compatibility matrix (potential -> engine):",
        "  (+ python modules present; the engine's own packages -- an ML-IAP build,",
        "   a MACE pair style -- are verified by the engine at run time)",
    ]
    for route in report["routes"]:
        if route["status"] != "supported":
            mark = "x"
            state = "unsupported"
        elif route["modules_available"]:
            mark = "+"
            state = "modules present"
        else:
            mark = "-"
            state = "missing: " + ", ".join(route["missing_modules"])
        lines.append(
            f"  [{mark}] {route['potential']:>6} -> {route['engine']:<7} "
            f"{route['implementation']:<22} {state}"
        )
        lines.append(f"        {route['summary']}")
        if route["status"] != "supported" and route["reason"]:
            lines.append(f"        reason: {route['reason']}")

    lines.append("")
    lines.append("Packages:")
    for name, version in report["packages"].items():
        lines.append(f"  {name:<14} {version or '(not installed)'}")

    device = report["device"]
    lines.append("")
    lines.append(f"Device: {device['platform']} ({device['cpu_count']} CPUs)")
    torch_info = device.get("torch")
    if torch_info and torch_info.get("probed") is False:
        lines.append(
            f"  torch {torch_info.get('version')} installed; CUDA not probed "
            "(use --probe-torch)"
        )
    elif torch_info and "error" not in torch_info:
        lines.append(
            f"  torch {torch_info['version']}, CUDA available: "
            f"{torch_info['cuda_available']} ({torch_info['device_count']} device(s))"
        )
        for entry in torch_info.get("devices", []):
            lines.append(f"    [{entry['index']}] {entry.get('name')}")
    elif torch_info is None:
        lines.append("  torch is not installed (no GPU routes available)")

    lines.append("")
    lines.append("Energy conventions:")
    for name, description in report["energy_conventions"].items():
        lines.append(f"  {name:<12} {description}")

    if "potential" in report:
        lines.append("")
        lines.append(_format_potential(report["potential"]))
        lines.append("")
        lines.append(f"Selected route: {report['selected_route']['implementation']}")
        availability = report["route_availability"]
        lines.append(
            "  available: "
            + ("yes" if availability["available"] else f"no -- {availability['detail']}")
        )
        lines.append(_format_capabilities(report["capabilities"]))
    return "\n".join(lines)


def _format_potential(potential: dict) -> str:
    lines = [f"Potential: {potential['label']} ({potential['kind']})"]
    lines.append(f"  declared elements: {', '.join(potential['declared_elements'])}")
    lines.append(f"  energy convention: {potential['energy_convention']}")
    if potential.get("model_path"):
        lines.append(f"  model: {potential['model_path']}")
        lines.append(f"  device/precision: {potential.get('device')} / {potential.get('precision')}")
    if potential.get("rendered_commands"):
        lines.append("  rendered LAMMPS commands:")
        lines += [f"    {command}" for command in potential["rendered_commands"]]
    discovered = potential.get("discovered")
    if discovered:
        lines.append("  discovered from the model file:")
        for key, value in discovered.items():
            if key == "atomic_reference_energies_eV" and value:
                value = ", ".join(f"{k}={v:.4f}" for k, v in value.items())
            lines.append(f"    {key}: {value}")
    elif potential.get("discovery_note"):
        lines.append(f"  {potential['discovery_note']}")
    for conflict in potential.get("conflicts", []):
        lines.append(f"  CONFLICT: {conflict}")
    return "\n".join(lines)


def _format_capabilities(capabilities: dict) -> str:
    flags = [
        name
        for name in ("energy", "forces", "stress", "per_atom_energy", "periodic", "gpu")
        if capabilities[name]
    ]
    lines = ["  capabilities: " + (", ".join(flags) or "none")]
    if capabilities["elements"]:
        lines.append("  elements: " + ", ".join(capabilities["elements"]))
    lines.append(
        f"  native units: {capabilities['native_units']}; energy convention: "
        f"{capabilities['native_energy_convention']}"
    )
    for note in capabilities["notes"]:
        lines.append(f"  note: {note}")
    return "\n".join(lines)


def _format_validate(report: dict) -> str:
    route = report["route"]
    lines = [
        f"Route: {route['potential']} -> {route['engine']} ({route['implementation']})",
        f"  {route['summary']}",
        _format_capabilities(report["capabilities"]),
    ]
    availability = report["availability"]
    lines.append(
        "  runnable here: "
        + ("yes" if availability["available"] else f"no -- {availability['detail']}")
    )
    if report["structure"]:
        structure = report["structure"]
        pbc = "".join("T" if p else "F" for p in structure["pbc"])
        lines.append(
            f"  structure: {structure['n_atoms']} atoms, elements "
            f"{', '.join(structure['elements'])}, pbc={pbc}, "
            f"fixed atoms={structure['constraints']['n_fixed']}, "
            f"sha256={structure['sha256'][:16]}..."
        )
    parameters = report["engine_parameters"]
    if parameters.get("pair_commands"):
        lines.append("  rendered LAMMPS commands:")
        lines += [f"    {command}" for command in parameters["pair_commands"]]
    if parameters.get("returnEnergyType"):
        lines.append(f"  OpenMM-ML returnEnergyType: {parameters['returnEnergyType']}")
    lines.append("")
    lines.append("Configuration is valid for this route." if report["ok"] else "INVALID")
    return "\n".join(lines)


def _format_singlepoint(report: dict) -> str:
    result = report["results"]["singlepoint"]
    lines = [
        f"Single point via {result['engine']} ({result['implementation']})",
        f"  energy: {result['energy_eV']:.8f} eV "
        f"({result['energy_per_atom_eV']:.8f} eV/atom, {result['energy_convention']})",
        f"  max force: {result['max_force_eV_per_A']:.8f} eV/Angstrom",
        f"  atoms: {result['n_atoms']}, native units: {result['native_units']}",
        f"  manifest: {report['manifest_path']}",
    ]
    return "\n".join(lines)


def _format_smoke(report: dict) -> str:
    trajectory = report["results"]["trajectory"]
    diagnostics = trajectory.get("diagnostics") or {}
    completed = trajectory.get("steps_completed")
    lines = [
        f"Smoke MD: {trajectory['steps']} steps of {trajectory['timestep_fs']} fs "
        f"({trajectory['ensemble']}); engine reports "
        f"{completed if completed is not None else 'no'} steps completed",
        f"  frames written: "
        f"{trajectory['frames_written'] if trajectory['frames_written'] is not None else 'unverified'}"
        f" -> {trajectory['trajectory_path']}",
    ]
    if trajectory["temperature_start_K"] is not None:
        lines.append(
            f"  temperature: {trajectory['temperature_start_K']:.1f} K -> "
            f"{trajectory['temperature_end_K']:.1f} K "
            f"(peak {trajectory['max_temperature_K']:.1f} K)"
        )
    if "energy_drift_eV_per_atom_per_ps" in diagnostics:
        lines.append(
            f"  NVE energy drift: {diagnostics['energy_drift_eV_per_atom_per_ps']:.3e} "
            f"eV/atom/ps (max excursion "
            f"{diagnostics['max_abs_energy_excursion_eV_per_atom']:.3e} eV/atom)"
        )
    elif "total_energy_change_eV_per_atom" in diagnostics:
        lines.append(
            f"  total-energy change: {diagnostics['total_energy_change_eV_per_atom']:.3e} "
            f"eV/atom ({diagnostics['label']})"
        )
        conserved = diagnostics.get("conserved_quantity_drift")
        if conserved:
            lines.append(
                f"  {diagnostics['conserved_quantity']} drift: "
                f"{conserved['energy_drift_eV_per_atom_per_ps']:.3e} eV/atom/ps"
            )
    elif diagnostics.get("available") is False:
        lines.append(f"  energy diagnostics: none ({diagnostics.get('reason')})")
    lines.append(f"  manifest: {report['manifest_path']}")
    return "\n".join(lines)


def _format_compare(report: dict) -> str:
    lines = [
        f"Cross-engine single-point comparison (reference: {report['reference_engine']})",
        f"  harmonised energy convention: {report['energy_convention']}",
        f"  atomic reference energies known: {report['atomic_reference_energies_known']}",
        f"  structure: {report['structure']['n_atoms']} atoms, "
        f"sha256={report['structure']['sha256'][:16]}...",
        "",
    ]
    for name, result in report["results"].items():
        lines.append(
            f"  {name:<8} E = {result['energy_eV']:.8f} eV "
            f"({result['energy_convention']}, native {result['native_units']})"
        )
    lines.append("")
    for pair, comparison in report["comparisons"].items():
        lines.append(f"  {pair}")
        lines.append(
            f"    dE/N            {comparison['delta_energy_per_atom_eV']:+.3e} eV/atom"
        )
        lines.append(
            f"    force RMSE      {comparison['force_rmse_eV_per_A']:.3e} eV/Angstrom"
        )
        lines.append(
            f"    max |dF|        {comparison['force_max_abs_error_eV_per_A']:.3e} eV/Angstrom"
        )
        if comparison["stress_max_abs_error_eV_per_A3"] is not None:
            lines.append(
                f"    max |dsigma|    "
                f"{comparison['stress_max_abs_error_eV_per_A3']:.3e} eV/Angstrom^3"
            )
        for note in comparison["notes"]:
            lines.append(f"    note: {note}")
    return "\n".join(lines)


__all__ = ["COMMANDS", "add_parser", "run"]
