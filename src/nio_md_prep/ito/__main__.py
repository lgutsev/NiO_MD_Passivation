"""``python -m nio_md_prep.ito`` command line for the ITO extension."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m nio_md_prep.ito")
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("build-substrate", help="build a slab from inputs/ito/surfaces/<model>/model.toml")
    s.add_argument("model", type=Path)
    s.add_argument("--output", type=Path, help="defaults to the model.toml directory")
    s = sub.add_parser("build-pilot", help="assemble a pilot study (studies/ito/*.toml)")
    s.add_argument("study", type=Path); s.add_argument("--output", type=Path, required=True)
    s.add_argument("--packmol-seed", type=int); s.add_argument("--velocity-seed", type=int)
    s.add_argument("--parameter-set"); s.add_argument("--no-lammps", action="store_true")
    s = sub.add_parser("adsorption-scan", help="single-molecule placement/energy checks on rigid slabs")
    s.add_argument("config", type=Path); s.add_argument("--output", type=Path, required=True)
    s.add_argument("--workers", type=int)
    s = sub.add_parser("analyze", help="contacts, orientation, clustering and coverage of an ITO pilot trajectory")
    s.add_argument("build_directory", type=Path); s.add_argument("--trajectory", type=Path, required=True)
    s.add_argument("--output", type=Path); s.add_argument("--last-frames", type=int, default=20)
    a = p.parse_args(argv)
    if a.command == "build-substrate":
        from .substrate import build_from_model
        m = build_from_model(a.model, a.output or a.model.parent)
        print(json.dumps({k: m[k] for k in ("model_id", "atom_count", "counts", "total_charge", "dipole_z_e_per_angstrom")}, indent=2))
    elif a.command == "build-pilot":
        from .assemble import build_pilot
        out = build_pilot(a.study, a.output, a.packmol_seed, a.velocity_seed, a.parameter_set, not a.no_lammps)
        print((out / "validation_report.txt").read_text())
    elif a.command == "adsorption-scan":
        from .checks import adsorption_scan
        print(adsorption_scan(a.config, a.output, a.workers))
    elif a.command == "analyze":
        from .analysis import analyze
        print(analyze(a.build_directory, a.trajectory, a.output, last_frames=a.last_frames))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
