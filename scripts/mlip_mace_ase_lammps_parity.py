#!/usr/bin/env python
"""Evaluate one MACE model on one structure through ASE and through LAMMPS/ML-IAP, and compare.

This is the preferred Phase-2 acceptance check, meant for a machine whose LAMMPS
build has ML-IAP with Python support (and KOKKOS for GPU / multi-layer MACE), such
as a LONI GPU node. It

1. optionally exports the checkpoint with ``mace_create_lammps_model --format mliap``
   (the export command, its return code and the export's sha256 are recorded);
2. resolves and saves the LAMMPS execution plan (what ``nio-md-prep mlip validate``
   shows) and refuses to continue if the route is not runnable;
3. runs the single point through the ASE bridge (reference) and the LAMMPS/ML-IAP
   bridge with ``jobs.compare_engines`` in the *total* energy convention;
4. writes the raw comparison report and a CANDIDATE tolerance file. The candidate
   is a measurement for human review, never an acceptance criterion.

Metrics: dE/atom, force RMSE, maximum per-atom force-vector error, maximum force
component error, and stress only when both routes report a virial.

Example (LONI GPU node, after activating an environment with mace-torch and a
LAMMPS Python module built with ML-IAP + PYTHON + KOKKOS):

    python scripts/mlip_mace_ase_lammps_parity.py \
        --model /path/to/model.model --elements Ni,O,C,H,N,P \
        --structure /path/to/POSCAR --pbc T,T,T --device cuda --export \
        --output parity/model-structure
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

_REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if _REPO_SRC.is_dir() and str(_REPO_SRC) not in sys.path:
    sys.path.insert(0, str(_REPO_SRC))

from nio_md_prep.mlip import jobs, tolerances  # noqa: E402
from nio_md_prep.mlip.provenance import git_commit, package_versions, sha256_file  # noqa: E402
from nio_md_prep.mlip.registry import build_bridge  # noqa: E402
from nio_md_prep.mlip.specs import (  # noqa: E402
    EngineSpec,
    JobSpec,
    MacePotentialSpec,
    SimulationSpec,
    StructureSpec,
    as_singlepoint,
)
from nio_md_prep.mlip.structures import load_structure  # noqa: E402

MLIAP_SUFFIX = "-mliap_lammps.pt"  # mace-torch 0.3.16 mace/cli/create_lammps_model.py:107


def parse_bools(text: str) -> tuple[bool, bool, bool]:
    values = [part.strip().upper() for part in text.split(",")]
    if len(values) != 3 or any(v not in {"T", "F", "TRUE", "FALSE", "1", "0"} for v in values):
        raise argparse.ArgumentTypeError("pbc must be three comma-separated booleans, e.g. T,T,F")
    return tuple(v in {"T", "TRUE", "1"} for v in values)  # type: ignore[return-value]


def export_mliap(model: Path, *, dtype: str, device: str, command: str) -> dict:
    """Run MACE's own exporter; record exactly what ran and what it produced."""
    argv = [command, str(model), "--format", "mliap", "--dtype", dtype]
    env = dict(os.environ)
    if device == "cpu":
        # MACE reads these when the export is created and pickles them into it;
        # setting them later, at LAMMPS run time, has no effect.
        env["MACE_ALLOW_CPU"] = "true"
    completed = subprocess.run(argv, capture_output=True, text=True, env=env)
    exported = model.with_name(model.name + MLIAP_SUFFIX)
    return {
        "argv": argv,
        "env": {"MACE_ALLOW_CPU": env.get("MACE_ALLOW_CPU")},
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-4000:],
        "exported_model": str(exported),
        "exported_sha256": sha256_file(exported) if exported.exists() else None,
    }


def build_job(args) -> JobSpec:
    options = {}
    if args.exported:
        options["exported_model_path"] = str(args.exported)
    return JobSpec(
        potential=MacePotentialSpec(
            label=args.model.stem,
            model_path=args.model,
            declared_elements=tuple(args.elements),
            device=args.device,
            precision=args.precision,
            implementation=("mliap",),
        ),
        engine=EngineSpec(kind="lammps", options=options, runtime=args.lammps_runtime,
                          executable=args.lammps_executable),
        simulation=SimulationSpec(task="singlepoint", compute_stress=args.stress),
        name="mace-ase-lammps-parity",
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--model", type=Path, required=True, help="MACE checkpoint (.model)")
    parser.add_argument("--elements", type=lambda s: [e.strip() for e in s.split(",") if e.strip()],
                        required=True, help="comma list, e.g. Ni,O")
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--pbc", type=parse_bools, default=None, help="override per-axis PBC, e.g. T,T,F")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--precision", default="float64", choices=("float32", "float64"))
    parser.add_argument("--export", action="store_true", help="run mace_create_lammps_model --format mliap first")
    parser.add_argument("--export-command", default="mace_create_lammps_model")
    parser.add_argument("--exported", type=Path, default=None, help="explicit exported-model path override")
    parser.add_argument("--lammps-runtime", default="auto", choices=("auto", "python", "executable"))
    parser.add_argument("--lammps-executable", default=None)
    parser.add_argument("--stress", action="store_true", help="also compare stress (fully periodic cells only)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.output.exists() and any(args.output.iterdir()):
        parser.error(f"{args.output} is not empty; choose a new output directory")
    args.output.mkdir(parents=True, exist_ok=True)
    record: dict = {
        "model": str(args.model),
        "model_sha256": sha256_file(args.model),
        "structure": str(args.structure),
        "structure_file_sha256": sha256_file(args.structure),
        "device": args.device,
        "precision": args.precision,
        "host": platform.node(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "packages": package_versions(),
        "git": git_commit(),
    }

    if args.export:
        record["export"] = export_mliap(args.model, dtype=args.precision, device=args.device,
                                        command=args.export_command)
        if record["export"]["returncode"] != 0:
            (args.output / "parity_report.json").write_text(json.dumps(record, indent=2, default=str) + "\n")
            print("export failed; see parity_report.json", file=sys.stderr)
            return 3

    job = build_job(args)
    structure = StructureSpec(path=args.structure, pbc=args.pbc) if args.pbc else StructureSpec(path=args.structure)
    atoms = load_structure(structure)
    lammps_bridge = build_bridge(job.potential, job.engine)
    plan = lammps_bridge.execution_plan(as_singlepoint(job.simulation), atoms)
    record["lammps_execution_plan"] = plan
    if not plan["availability"]["available"] or plan["unmet_capabilities"]:
        (args.output / "parity_report.json").write_text(json.dumps(record, indent=2, default=str) + "\n")
        availability = plan["availability"]
        print("the LAMMPS/ML-IAP route is not runnable here:",
              availability.get("detail") or availability.get("missing"),
              plan["unmet_capabilities"] or "", file=sys.stderr)
        return 4

    report = jobs.compare_engines(job, ["ase", "lammps"], output_dir=args.output / "compare",
                                  structure=structure, convention="total")
    record["comparison"] = report
    candidate = tolerances.candidate_record(report, measured_on={
        key: record[key] for key in ("host", "platform", "python", "packages", "git",
                                     "model_sha256", "device", "precision")
    } | {"exported_sha256": (record.get("export") or {}).get("exported_sha256")})
    (args.output / "parity_report.json").write_text(json.dumps(record, indent=2, default=str) + "\n")
    (args.output / "candidate_tolerances.json").write_text(json.dumps(candidate, indent=2) + "\n")

    for pair, metrics in tolerances.measured_metrics(report).items():
        print(pair)
        for metric, value in metrics.items():
            print(f"  {metric:38s} {'n/a' if value is None else f'{value:.3e}'}")
    print(f"raw report: {args.output / 'parity_report.json'}")
    print(f"CANDIDATE tolerances (not accepted; review before use): {args.output / 'candidate_tolerances.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
