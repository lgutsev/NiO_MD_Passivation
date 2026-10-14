"""``nio-md-prep dataset scan|export|split|audit`` (argument parsing and dispatch).

Importing this module and building the parser import only the standard library, so
``nio-md-prep dataset --help`` works without numpy/ASE; the dataset modules are imported when a
subcommand runs. Exit status: 0 success; 1 an audit found problems; 2 a deliberate refusal
(``error: ...`` from the top-level CLI: bad arguments, overwrite refused, several reference pools,
leakage, missing dependency, ...).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Sequence

ENERGY_QUANTITIES = ("free_energy", "energy_sigma0", "energy_no_entropy")


def _add_selection(p: argparse.ArgumentParser) -> None:
    p.add_argument("roots", nargs="*", metavar="[ALIAS=]ROOT", help="calculation roots (also --root)")
    p.add_argument("--root", action="append", default=[], metavar="ALIAS=PATH",
                   help="calculation root with an explicit alias (repeatable); the alias is part of every run id")
    p.add_argument("--inventory", type=Path, help="selection/lineage TOML ([[run]] path/include/lineage/family/...)")
    p.add_argument("--settings-overrides", type=Path,
                   help="TOML of reviewed [[equivalence]] tables for reference-settings pooling")
    p.add_argument("--lineage-policy", metavar="run|parent|depth:N",
                   help="declared grouping for runs without inventory/agglomeration/InterfaceForge lineage")
    p.add_argument("--pool", metavar="ID", help="reference-settings pool to export (id or unique prefix, >= 6 chars)")
    p.add_argument("--include-excluded-part", action="append", default=[], metavar="PATTERN",
                   help="deliberately scan a default-excluded path part (e.g. precondition, .interfaceforge, 'X*')")
    p.add_argument("--exclude-glob", action="append", default=[], metavar="GLOB", help="posix relpath glob to skip")
    p.add_argument("--include-glob", action="append", default=[], metavar="GLOB",
                   help="posix relpath glob scanned even inside default-excluded trees")
    p.add_argument("--follow-symlinks", action="store_true", help="follow directory symlinks (loops are detected)")
    # policy
    g = p.add_argument_group("label policy (strict defaults; every change is recorded in the manifest)")
    g.add_argument("--policy-file", type=Path, help="TOML/JSON table of acceptance.Policy fields (flags override it)")
    g.add_argument("--energy-quantity", choices=ENERGY_QUANTITIES, help="label energy (default free_energy, "
                   "the force-consistent F)")
    g.add_argument("--allow-energy-only", action="store_true",
                   help="export frames without forces as label_set=energy_only (never the default)")
    g.add_argument("--include-stress", action="store_true", help="export valid stress tensors (default: omitted "
                   "with stress_reason=not_requested)")
    g.add_argument("--allow-pstress", action="store_true", help="accept stress from runs with PSTRESS != 0")
    g.add_argument("--allow-vacuum-stress", action="store_true", help="accept stress of slabs/clusters with vacuum")
    g.add_argument("--interrupted-runs", choices=("exclude", "recover"),
                   help="incomplete runs: exclude (default) or judge complete frames individually")
    g.add_argument("--allow-ediff-zero", action="store_true", help="accept EDIFF=0 runs (criterion disabled)")
    g.add_argument("--allowed-calc-types", metavar="LIST", help="comma list of md,relaxation,static,other")
    g.add_argument("--force-outlier-threshold", type=float, metavar="EV_PER_A",
                   help="quarantine frames with max |F| above this")
    g.add_argument("--near-duplicate-tolerance", type=float, metavar="A",
                   help="link (never remove) near-duplicate frames across groups (max displacement in A)")
    g.add_argument("--magnetic-species", metavar="LIST", help="comma list of magnetic elements (default Co,Cr,Cu,Fe,Mn,Ni,V)")
    g.add_argument("--accept-uncontrolled-magnetism", action="store_true",
                   help="OVERRIDE: accept ISPIN=2 runs with magnetic species but no explicit MAGMOM")
    g.add_argument("--accept-unknown-magnetism", action="store_true",
                   help="OVERRIDE: accept spin-polarized runs without magnetization evidence")
    g.add_argument("--accept-magnetic-transitions", action="store_true",
                   help="OVERRIDE: accept frames after a magnetic state change")
    g.add_argument("--magnetic-override-reason", metavar="TEXT",
                   help="required with any magnetic override; recorded in the manifest policy")
    # subsampling
    s = p.add_argument_group("declared subsampling (dropped frames are excluded/subsampled)")
    s.add_argument("--stride", type=int, default=1, help="keep every N-th accepted frame of each run")
    s.add_argument("--max-frames-per-run", type=int, help="keep at most N evenly spaced accepted frames per run")


def add_parser(sub) -> argparse.ArgumentParser:
    """Attach ``dataset`` to the top-level ``nio-md-prep`` subparsers."""
    dataset = sub.add_parser(
        "dataset", help="audited VASP -> MLIP extended-XYZ dataset (scan, export, split, audit)",
        description="Audited VASP -> MLIP dataset export. See docs/dataset-export.md.",
    )
    commands = dataset.add_subparsers(dest="dataset_command", required=True)

    scan = commands.add_parser("scan", help="classify every run and frame; write the accounting only")
    _add_selection(scan)
    scan.add_argument("--output", type=Path, help="write runs.jsonl/frames.jsonl/exclusions.csv/audit files here")
    scan.add_argument("--force", action="store_true", help="write into a non-empty output directory")
    scan.add_argument("--json", action="store_true", help="print the scan manifest as JSON")

    export = commands.add_parser("export", help="write dataset.extxyz + accounting + dataset_manifest.json")
    _add_selection(export)
    export.add_argument("--output", type=Path, required=True)
    export.add_argument("--force", action="store_true", help="write into a non-empty output directory")
    export.add_argument("--energy-key", default="REF_energy")
    export.add_argument("--forces-key", default="REF_forces")
    export.add_argument("--stress-key", default="REF_stress")

    split = commands.add_parser("split", help="lineage-group train/valid/test split of an export")
    split.add_argument("export_dir", type=Path)
    split.add_argument("--output", type=Path, required=True)
    split.add_argument("--seed", type=int, default=11)
    split.add_argument("--fractions", default="0.8,0.1,0.1", metavar="TRAIN,VALID,TEST")
    split.add_argument("--stratify-by", metavar="KEY")
    split.add_argument("--allow-underpopulated-strata", action="store_true")
    split.add_argument("--previous", type=Path, metavar="SPLIT_MANIFEST",
                       help="keep earlier assignments; the test split is frozen unless --grow-test")
    split.add_argument("--grow-test", action="store_true")
    split.add_argument("--force", action="store_true")

    audit = commands.add_parser("audit", help="re-verify an export (hashes, ASE round trip, split leakage)")
    audit.add_argument("export_dir", type=Path)
    audit.add_argument("--split", action="append", default=[], type=Path, metavar="SPLIT_DIR")
    audit.add_argument("--verify-sources", action="append", default=[], metavar="ALIAS=PATH",
                       help="re-hash the source vasprun files under (relocated) roots")
    audit.add_argument("--output", type=Path, help="report directory (default EXPORT_DIR/audit)")
    audit.add_argument("--force", action="store_true")
    return dataset


def _csv(text: str | None) -> list[str] | None:
    if text is None:
        return None
    return [item.strip() for item in text.split(",") if item.strip()]


def _policy(a) -> Any:
    from .acceptance import Policy
    from .errors import DatasetError

    base: dict[str, Any] = {}
    if a.policy_file is not None:
        text = Path(a.policy_file).read_text(encoding="utf-8")
        if a.policy_file.suffix.lower() == ".json":
            base = json.loads(text)
        else:
            import tomllib

            try:
                base = tomllib.loads(text)
            except tomllib.TOMLDecodeError as exc:
                raise DatasetError(f"policy file {a.policy_file}: {exc}") from None
        base = dict(base.get("policy", base))
    flags = {
        "energy_quantity": a.energy_quantity,
        "interrupted_run_policy": a.interrupted_runs,
        "force_outlier_threshold": a.force_outlier_threshold,
        "near_duplicate_tolerance": a.near_duplicate_tolerance,
        "magnetic_override_reason": a.magnetic_override_reason,
    }
    for name, value in flags.items():
        if value is not None:
            base[name] = value
    for name, attr in (("allow_energy_only", "allow_energy_only"), ("include_stress", "include_stress"),
                       ("allow_pstress", "allow_pstress"), ("allow_vacuum_stress", "allow_vacuum_stress"),
                       ("allow_ediff_zero", "allow_ediff_zero"),
                       ("accept_uncontrolled", "accept_uncontrolled_magnetism"),
                       ("accept_unknown", "accept_unknown_magnetism"),
                       ("accept_transitions", "accept_magnetic_transitions")):
        if getattr(a, attr):
            base[name] = True
    if a.allowed_calc_types is not None:
        base["allowed_calc_types"] = _csv(a.allowed_calc_types)
    if a.magnetic_species is not None:
        base["magnetic_species"] = _csv(a.magnetic_species)
    return Policy.from_mapping(base)


def _options(a, *, label_keys=None):
    from .export import ExportOptions
    from .extxyz import LabelKeys, validate_label_keys

    keys = label_keys or LabelKeys()
    if label_keys is None and hasattr(a, "energy_key"):
        keys = validate_label_keys(a.energy_key, a.forces_key, a.stress_key)
    return ExportOptions(
        inventory=a.inventory, settings_overrides=a.settings_overrides, lineage_policy=a.lineage_policy,
        pool=a.pool, policy=_policy(a), label_keys=keys, include_parts=tuple(a.include_excluded_part),
        exclude_globs=tuple(a.exclude_glob), include_globs=tuple(a.include_glob),
        follow_symlinks=a.follow_symlinks, stride=a.stride, max_frames_per_run=a.max_frames_per_run,
    )


def _roots(a) -> list[str]:
    from .errors import DatasetError

    roots = list(a.roots) + list(a.root)
    if not roots:
        raise DatasetError("give at least one calculation root (positional or --root ALIAS=PATH)")
    return roots


def _print_counts(manifest: dict[str, Any], out) -> None:
    counts = manifest["counts"]
    print(f"runs: {counts['runs']} ({counts['runs_parsed']} parsed)", file=out)
    for key, value in counts["runs_by_reason"].items():
        print(f"  run   {key}: {value}", file=out)
    print(f"frames: {counts['frames']} (label-accepted {counts['frames_label_accepted']}, "
          f"exportable {counts['frames_exported']})", file=out)
    for key, value in counts["frames_by_reason"].items():
        print(f"  frame {key}: {value}", file=out)
    pools = manifest.get("pools") or []
    if pools:
        print(f"reference-settings pools: {len(pools)}", file=out)
        for pool in pools:
            mark = " (selected)" if pool["selected"] else ""
            print(f"  {pool['pool_id']}: {pool['runs']} run(s), {pool['accepted_frames']} accepted frame(s){mark}",
                  file=out)
    unresolved = (manifest.get("lineage") or {}).get("unresolved_runs") or []
    if unresolved:
        print(f"lineage unresolved for {len(unresolved)} run(s) (quarantined); declare --lineage-policy "
              "or inventory lineage", file=out)


def run(a, argv: Sequence[str] | None = None) -> int:
    """Dispatch a parsed ``dataset`` command; returns the exit status.

    ``argv`` is the argument list the caller parsed (recorded in the manifest ``invocation``);
    ``None`` means the process command line.
    """
    command = a.dataset_command
    argv = [str(x) for x in argv] if argv is not None else (list(sys.argv[1:]) if sys.argv else None)
    if command == "scan":
        from .export import scan

        result = scan(_roots(a), output=a.output, force=a.force, options=_options(a), argv=argv)
        if a.json:
            print(json.dumps(result.manifest, indent=2, sort_keys=True))
        else:
            _print_counts(result.manifest, sys.stdout)
            if result.manifest.get("pool_selection_required"):
                print("export will need --pool: several reference-settings pools hold accepted frames")
            if a.output is not None:
                print(f"accounting written to {result.output}")
        return 0
    if command == "export":
        from .export import export

        result = export(_roots(a), a.output, force=a.force, options=_options(a), argv=argv)
        _print_counts(result.manifest, sys.stdout)
        files = result.manifest["files"]["dataset.extxyz"]
        print(f"dataset.extxyz: {files['frames']} frames in {result.manifest['groups']['exported_groups']} lineage "
              f"groups, ASE round trip verified; manifest {result.output / 'dataset_manifest.json'}")
        return 0
    if command == "split":
        from .errors import DatasetError
        from .split import split_export

        try:
            fractions = tuple(float(x) for x in a.fractions.split(","))
        except ValueError:
            raise DatasetError(f"--fractions must be three numbers like 0.8,0.1,0.1, got {a.fractions!r}") from None
        manifest = split_export(a.export_dir, a.output, seed=a.seed, fractions=fractions, stratify_by=a.stratify_by,
                                previous=a.previous, grow_test=a.grow_test,
                                allow_underpopulated_strata=a.allow_underpopulated_strata, force=a.force)
        files = manifest["files"]
        print("split written to " + str(a.output) + ": " + ", ".join(
            f"{name} {files[name]['frames']} frames" for name in ("train", "valid", "test")))
        for warning in manifest.get("warnings") or []:
            print(f"warning: {warning}")
        return 0
    if command == "audit":
        from .audit import audit_export

        report = audit_export(a.export_dir, split_dirs=a.split, verify_sources=a.verify_sources or None,
                              output=a.output, force=a.force)
        for check in report["checks"]:
            print(f"{'ok    ' if check['ok'] else 'FAILED'} {check['name']}: {check['detail'][:200]}")
        print(f"audit report: {report['output']}")
        return 0 if report["ok"] else 1
    raise AssertionError(command)  # pragma: no cover - argparse restricts the choices


def main(argv: Sequence[str] | None = None) -> int:
    """Stand-alone entry (``python -m nio_md_prep.dataset.cli ...``), same exit codes as the main CLI."""
    parser = argparse.ArgumentParser(prog="nio-md-prep")
    sub = parser.add_subparsers(dest="command", required=True)
    add_parser(sub)
    args = parser.parse_args(argv)
    try:
        return run(args, argv=argv)
    except (ValueError, FileNotFoundError, FileExistsError, RuntimeError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0  # pragma: no cover


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
