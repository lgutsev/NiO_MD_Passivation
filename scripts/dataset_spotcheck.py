"""Independent spot check of an exported dataset against the original VASP files, using ASE only.

Reads N exported frames from ``dataset.extxyz`` (ASE extxyz reader) and the same ionic steps
straight from the source ``vasprun.xml`` with ASE's own VASP reader (``format="vasp-xml"``, not
the nio-md-prep parser), then compares the force-consistent free energy, the raw forces
(``apply_constraint=False``), species and positions. Prints per-frame and maximum absolute
differences; exit status 0 when every difference is within tolerance, 1 otherwise, 2 on usage
errors. Only the standard library and ASE are imported (no nio_md_prep code), so the comparison
is independent of the exporter's parser.

Usage::

    python scripts/dataset_spotcheck.py EXPORT_DIR [--n 20] [--seed 11]
        [--root ALIAS=PATH ...] [--energy-tol 1e-6] [--force-tol 1e-6] [--json OUT.json]

The exported label is F = e_fr_energy - PSTRESS*V. Current ASE (3.29 checked on ASE's own
``vasprun_pstress.xml``) also removes the PV term; older ASE versions may return the raw
``e_fr_energy``. With PSTRESS != 0 the comparison therefore accepts either convention and reports
which one matched per frame (``ase_energy_convention``); with PSTRESS = 0 both are identical.

``--root`` points an alias at a relocated calculation root; by default the root paths recorded
in ``dataset_manifest.json`` are used. The source step is located by counting the ``<calculation>``
(DFT) steps before the exported ``ionic_step`` in ``frames.jsonl`` (VASP-MLFF flat steps are not
``<calculation>`` blocks and ASE does not read them).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import random
import sys


def _load(export_dir: Path):
    manifest = json.loads((export_dir / "dataset_manifest.json").read_text(encoding="utf-8"))
    frames = []
    with (export_dir / "frames.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                frames.append(json.loads(line))
    runs = {}
    with (export_dir / "runs.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                record = json.loads(line)
                runs[record["run_id"]] = record
    return manifest, frames, runs


def _choose(exported: list[str], n: int, seed: int) -> list[str]:
    if n >= len(exported):
        return list(exported)
    chosen = set(random.Random(seed).sample(range(len(exported)), n))
    return [frame_id for i, frame_id in enumerate(exported) if i in chosen]


def _max_abs(a, b) -> float:
    import numpy as np

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return math.inf
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def _iread_until(iterator, last: int):
    """The first ``last + 1`` images of an ASE reader (stops before reading further steps)."""
    for k, atoms in enumerate(iterator):
        yield atoms
        if k >= last:
            return


def spotcheck(export_dir: Path, *, n: int = 20, seed: int = 11, roots: dict[str, Path] | None = None,
              energy_tol: float = 1e-6, force_tol: float = 1e-6, position_tol: float = 1e-6) -> dict:
    import ase.io

    export_dir = Path(export_dir)
    manifest, frames, runs = _load(export_dir)
    keys = manifest.get("label_keys") or {}
    energy_key = keys.get("energy", "REF_energy")
    forces_key = keys.get("forces", "REF_forces")
    quantity = manifest.get("energy_quantity", "free_energy")
    root_paths = {entry["alias"]: Path(entry["path"]) for entry in manifest.get("roots") or []}
    root_paths.update(roots or {})
    exported = [f["frame_id"] for f in frames if (f.get("outcome") or {}).get("status") == "accepted"]
    chosen = _choose(exported, n, seed)
    by_id = {f["frame_id"]: f for f in frames}
    # DFT ordinal of every exported step within its source file
    ordinal: dict[str, int] = {}
    counters: dict[str, int] = {}
    for record in sorted(frames, key=lambda f: (f["run_id"], f["index"])):
        if record.get("label_source") == "dft":
            ordinal[record["frame_id"]] = counters.get(record["run_id"], 0)
            counters[record["run_id"]] = counters.get(record["run_id"], 0) + 1
    # our side: the exported frames
    ours = {}
    wanted = set(chosen)
    for atoms in ase.io.iread(str(export_dir / "dataset.extxyz"), index=":", format="extxyz"):
        frame_id = atoms.info.get("frame_id")
        if frame_id in wanted:
            forces = atoms.arrays.get(forces_key)
            ours[frame_id] = {
                "energy": float(atoms.info[energy_key]),
                "free_energy": float(atoms.info.get("vasp_free_energy", atoms.info[energy_key])),
                "forces": None if forces is None else forces.copy(),
                "species": atoms.get_chemical_symbols(), "positions": atoms.positions.copy(),
            }
    # reference side: ASE's VASP reader on the original files, one pass per source file
    per_source: dict[Path, list[str]] = {}
    problems: list[str] = []
    for frame_id in chosen:
        record = by_id[frame_id]
        run = runs[record["run_id"]]
        alias, relpath = run["root_alias"], run["relpath"]
        if alias not in root_paths:
            problems.append(f"{frame_id}: no path for root alias {alias!r} (pass --root {alias}=PATH)")
            continue
        base = root_paths[alias] if relpath == "." else root_paths[alias] / relpath
        per_source.setdefault(base / run["source_file"], []).append(frame_id)
    results = []
    for source, ids in sorted(per_source.items()):
        if not source.is_file():
            problems.append(f"{source.as_posix()}: missing")
            continue
        targets = {ordinal[frame_id]: frame_id for frame_id in ids if frame_id in ordinal}
        last = max(targets) if targets else -1
        try:
            reference_steps = list(_iread_until(ase.io.iread(str(source), index=":", format="vasp-xml"), last))
        except Exception as exc:  # noqa: BLE001 - any reader failure is a spot-check problem
            problems.append(f"{source.as_posix()}: ASE could not read it ({type(exc).__name__}: {exc})")
            continue
        for k, atoms in enumerate(reference_steps):
            if k in targets:
                frame_id = targets[k]
                mine = ours.get(frame_id)
                if mine is None:
                    problems.append(f"{frame_id}: not found in dataset.extxyz")
                    continue
                reference_f = float(atoms.get_potential_energy(force_consistent=True))
                reference_forces = atoms.get_forces(apply_constraint=False)
                ours_f = mine["energy"] if quantity == "free_energy" else mine["free_energy"]
                pv_term = float((by_id[frame_id].get("energy") or {}).get("pv_term") or 0.0)
                d_minus_pv, d_raw = abs(ours_f - reference_f), abs(ours_f + pv_term - reference_f)
                convention = "F=e_fr-PV" if d_minus_pv <= d_raw else "raw e_fr"
                entry = {
                    "frame_id": frame_id, "source": source.name, "calculation_index": k,
                    "pv_term": pv_term, "ase_energy_convention": convention if pv_term else "PSTRESS=0",
                    "d_free_energy": min(d_minus_pv, d_raw),
                    "d_forces_max": None if mine["forces"] is None else _max_abs(mine["forces"], reference_forces),
                    "d_positions_max": _max_abs(mine["positions"], atoms.positions),
                    "species_equal": mine["species"] == atoms.get_chemical_symbols(),
                }
                results.append(entry)
        missing = sorted(set(targets.values()) - {r["frame_id"] for r in results})
        problems += [f"{frame_id}: calculation not found by ASE in {source.name}" for frame_id in missing]
    finite = results
    summary = {
        "export_dir": export_dir.as_posix(), "requested": n, "checked": len(results), "seed": seed,
        "energy_quantity": quantity,
        "max_d_free_energy": max((r["d_free_energy"] for r in finite), default=None),
        "max_d_forces": max((r["d_forces_max"] for r in finite if r["d_forces_max"] is not None), default=None),
        "max_d_positions": max((r["d_positions_max"] for r in finite), default=None),
        "tolerances": {"energy": energy_tol, "forces": force_tol, "positions": position_tol},
        "reader": "ase.io vasp-xml (independent of nio-md-prep's parser)",
        "frames": results, "problems": problems,
    }
    failures = [r["frame_id"] for r in results if r["d_free_energy"] > energy_tol
                or (r["d_forces_max"] is not None and r["d_forces_max"] > force_tol)
                or r["d_positions_max"] > position_tol or not r["species_equal"]]
    summary["failures"] = failures
    summary["ok"] = not failures and not problems and len(results) == len(chosen) and bool(results)
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("export_dir", type=Path)
    parser.add_argument("--n", type=int, default=20, help="number of exported frames to check (default 20)")
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--root", action="append", default=[], metavar="ALIAS=PATH")
    parser.add_argument("--energy-tol", type=float, default=1e-6, help="eV (default 1e-6)")
    parser.add_argument("--force-tol", type=float, default=1e-6, help="eV/A (default 1e-6)")
    parser.add_argument("--position-tol", type=float, default=1e-6, help="A (default 1e-6)")
    parser.add_argument("--json", type=Path, help="also write the full result as JSON")
    args = parser.parse_args(argv)
    roots = {}
    for value in args.root:
        if "=" not in value:
            parser.error(f"--root expects ALIAS=PATH, got {value!r}")
        alias, path = value.split("=", 1)
        roots[alias] = Path(path)
    import importlib.util

    if importlib.util.find_spec("ase") is None:
        print("error: ASE is required (pip install 'nio-md-prep[dataset]')", file=sys.stderr)
        return 2
    if not (args.export_dir / "dataset_manifest.json").is_file():
        print(f"error: {args.export_dir} is not an export directory", file=sys.stderr)
        return 2
    summary = spotcheck(args.export_dir, n=args.n, seed=args.seed, roots=roots, energy_tol=args.energy_tol,
                        force_tol=args.force_tol, position_tol=args.position_tol)
    for entry in summary["frames"]:
        forces = "n/a" if entry["d_forces_max"] is None else f"{entry['d_forces_max']:.3e}"
        print(f"{entry['frame_id']}: |dF_energy| {entry['d_free_energy']:.3e} eV, max |dF| {forces} eV/A, "
              f"max |dx| {entry['d_positions_max']:.3e} A")
    for problem in summary["problems"]:
        print(f"problem: {problem}")
    print(f"checked {summary['checked']} frame(s) (requested {args.n}); max |dE| {summary['max_d_free_energy']}, max |dF| {summary['max_d_forces']}, "
          f"max |dx| {summary['max_d_positions']}; {'OK' if summary['ok'] else 'FAILED'}")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if summary["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
