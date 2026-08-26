from __future__ import annotations

import argparse
import json
from pathlib import Path

from .build import build, refresh_inputs
from .config import missing_ligpargen, molecule_dir, molecule_manifest
from .geometry import elements
from .lammps import charge, parse
from .validate import validate
from .analysis.model_scope import MODEL_PROFILES

TEMPLATE='''[molecule]\ndisplay_name = "{name}"\nslug = "{slug}"\nrole = "secondary"\nexpected_net_charge = 0.0\nphosphonic_acid_anchors = 1\nanchor_atom_ids = []\n\n[files]\nligpargen = "ligpargen.lmp"\ngeometry = ""\noverride = ""\n'''


def _agglomeration_completion_message(
    manifest: dict,
    output: Path,
    *,
    structures_only: bool,
    regenerated_vasp: bool = False,
) -> str:
    case_root = output / manifest.get(
        "case_root", "structures" if structures_only else "vasp_runs"
    )
    completed_cases = sum(
        1
        for replica in manifest.get("replicas", [])
        for structure in replica.get("structures", [])
        if structure.get("path")
    )
    if structures_only:
        return (
            f"Agglomeration structures created under {case_root}: "
            f"{completed_cases} structure case(s)."
        )
    if manifest.get("xtb", {}).get("enabled", False):
        schedules = manifest.get("vasp_md", {}).get("schedules", [])
        stage_count = completed_cases * len(schedules)
        schedule_text = ", ".join(schedules)
        launcher = output / manifest.get("vasp_launcher", "launch_vasp_runs.sh")
        verb = "regenerated" if regenerated_vasp else "created"
        return (
            f"Agglomeration launch-ready VASP folders {verb} under {case_root}: "
            f"{completed_cases} structure case(s), {stage_count} temperature-stage "
            f"run folder(s) ({schedule_text}). Each structure case contains "
            f"submit_temperature_jobs.sh. Preview submissions with {launcher} "
            f"--dry-run; launch with confirmation using {launcher}."
        )
    return (
        f"Agglomeration launch-ready VASP folders created under {case_root}: "
        f"{completed_cases} structure case(s)."
    )


def _agglomeration_xtb_pending_message(manifest: dict, output: Path) -> str:
    structures = [
        structure
        for replica in manifest.get("replicas", [])
        for structure in replica.get("structures", [])
    ]
    completed = sum(bool(structure.get("path")) for structure in structures)
    pending = sum(
        structure.get("status") == "XTB_REQUIRED" for structure in structures
    )
    audit = output / "audit_xtb_runs.py"
    launcher = output / "launch_xtb.sh"
    return (
        f"Agglomeration xTB status for {output}: {completed} completed case(s), "
        f"{pending} pending case(s). Audit stages with {audit}; preview a "
        f"resumable submission with {launcher} --dry-run, or launch with "
        f"confirmation using {launcher}. Completed cases are detected and "
        "skipped. Rerun this preparation command with the same output path "
        "after the remaining array tasks finish."
    )


def main(argv=None)->int:
    p=argparse.ArgumentParser(prog="nio-md-prep"); sub=p.add_subparsers(dest="command",required=True)
    i=sub.add_parser("init-molecule"); i.add_argument("name")
    i=sub.add_parser("inspect-molecule"); i.add_argument("name")
    i=sub.add_parser("build"); i.add_argument("config",type=Path); i.add_argument("--output",type=Path,required=True); i.add_argument("--packed-xyz",type=Path); i.add_argument("--packmol-seed",type=int); i.add_argument("--velocity-seed",type=int)
    i=sub.add_parser("refresh-inputs"); i.add_argument("config",type=Path); i.add_argument("--output",type=Path,required=True)
    i=sub.add_parser("validate"); i.add_argument("output",type=Path); i.add_argument("--primary-final",type=Path)
    i=sub.add_parser("prepare-agglomeration")
    i.add_argument("config",type=Path)
    i.add_argument("--output",type=Path,required=True)
    agglomeration_mode=i.add_mutually_exclusive_group(required=True)
    agglomeration_mode.add_argument(
        "--reference-dir",
        action="append",
        metavar="[SLUG=]DIR",
        help=(
            "complete VASP template; repeat as SLUG=DIR for mixed agglomerates "
            "(the first qualified reference supplies shared calculation inputs)"
        ),
    )
    agglomeration_mode.add_argument("--structures-only",action="store_true",help="write POSCAR/XYZ structures without VASP calculation inputs")
    agglomeration_mode.add_argument(
        "--regenerate-xtb-launcher",
        action="store_true",
        help=(
            "rebuild only the xTB audit, workq carousel, and launcher files "
            "for an existing pure or mixed campaign"
        ),
    )
    i.add_argument(
        "--vasp-template-slug",
        help="qualified reference whose INCAR, KPOINTS, launchers, and other shared files are reused",
    )
    i.add_argument("--packed-xyz",type=Path,help="use a supplied ordered Packmol XYZ for a one-replica offline build")
    i.add_argument(
        "--regenerate-vasp",
        action="store_true",
        help=(
            "rebuild VASP folders and the campaign launcher from completed xTB "
            "results in an existing output directory"
        ),
    )
    i=sub.add_parser(
        "audit-agglomerations",
        help="audit xTB progress across every agglomeration campaign under a root",
    )
    i.add_argument("root",type=Path,nargs="?",default=Path("."))
    i.add_argument(
        "--output-prefix",
        type=Path,
        help="report path without extension (default: ROOT/agglomeration_xtb_audit)",
    )
    i=sub.add_parser("archive-runs")
    i.add_argument("roots",nargs="+",type=Path,help="prepared/work roots to archive")
    i.add_argument("--output",type=Path,required=True)
    i.add_argument("--force",action="store_true",help="replace an existing archive")
    i=sub.add_parser("prepare-sequential-stage2"); i.add_argument("config",type=Path); i.add_argument("--primary-final",type=Path,required=True); i.add_argument("--output",type=Path,required=True); i.add_argument("--packed-xyz",type=Path); i.add_argument("--packmol-seed",type=int); i.add_argument("--velocity-seed",type=int)
    i=sub.add_parser("prepare-lego-stage2")
    i.add_argument("config",type=Path)
    i.add_argument("--primary-final",type=Path,required=True)
    i.add_argument("--coverage-map",type=Path,required=True)
    i.add_argument("--output",type=Path,required=True)
    i.add_argument("--packed-xyz",type=Path)
    i.add_argument("--packmol-seed",type=int)
    i.add_argument("--velocity-seed",type=int)
    i.add_argument("--occupancy-threshold",type=float,default=0.20)
    i.add_argument("--packing-inset",type=float,default=4.0)
    i.add_argument("--minimum-gap-fraction",type=float,default=0.12)
    i.add_argument("--minimum-tunnel-width",type=float,default=12.0)
    i.add_argument("--deposition-clearance",type=float,default=30.0)
    i.add_argument("--final-deposition-clearance",type=float,default=15.0)
    i.add_argument("--continuation-steps",type=int,default=300000)
    i.add_argument("--lateral-wall",action="store_true")
    i.add_argument("--wall-strength",type=float,default=0.50)
    i.add_argument("--wall-cutoff",type=float,default=3.0)
    i=sub.add_parser("prepare-lego2-stage2")
    i.add_argument("config",type=Path)
    i.add_argument("--primary-final",type=Path,required=True)
    i.add_argument("--coverage-map",type=Path,required=True)
    i.add_argument("--output",type=Path,required=True)
    i.add_argument("--packed-xyz",type=Path)
    i.add_argument("--packmol-seed",type=int)
    i.add_argument("--velocity-seed",type=int)
    i.add_argument("--occupancy-threshold",type=float,default=0.20)
    i.add_argument("--packing-inset",type=float,default=4.0)
    i.add_argument("--minimum-void-area-fraction",type=float,default=0.02)
    i.add_argument("--minimum-packing-width",type=float,default=12.0)
    i.add_argument("--minimum-packing-height",type=float,default=12.0)
    i.add_argument("--deposition-clearance",type=float,default=30.0)
    i.add_argument("--final-deposition-clearance",type=float,default=15.0)
    i.add_argument("--continuation-steps",type=int,default=300000)
    i=sub.add_parser("prepare-lego-continuation")
    i.add_argument("config",type=Path)
    i.add_argument("--output",type=Path,required=True)
    i.add_argument("--additional-drop",type=float,default=15.0)
    i.add_argument("--steps",type=int,default=300000)
    i=sub.add_parser("analyze-coverage")
    i.add_argument("build_directory",type=Path)
    i.add_argument("--trajectory",type=Path)
    i.add_argument("--output",type=Path)
    i.add_argument("--grid-spacing",type=float,default=0.20)
    i.add_argument("--radius-scale",type=float,default=1.0)
    i.add_argument("--roughness-grid-spacing",type=float,default=2.0,help="grid resolution for the top-of-film height/roughness map")
    i.add_argument("--near-surface-height",type=float,default=5.0,help="height above the local substrate reference (A) gating the near-surface/anchor-conditioned coverage metrics")
    i.add_argument("--last-frames",type=int,default=100,help="analyze the final N frames; use 0 for all frames")
    i.add_argument("--stride",type=int,default=1)
    i.add_argument("--blocks",type=int,default=5)
    i.add_argument("--exclude-hydrogen",action="store_true")
    i.add_argument("--timestep-fs",type=float)
    i=sub.add_parser("summarize-coverage")
    i.add_argument("prepared_root",type=Path)
    i.add_argument("--output",type=Path,required=True)
    i=sub.add_parser("analyze-interface")
    i.add_argument("build_directory",type=Path)
    i.add_argument("--trajectory",type=Path)
    i.add_argument("--output",type=Path)
    i.add_argument("--last-frames",type=int,default=100,help="analyze the final N frames; use 0 for all frames")
    i.add_argument("--stride",type=int,default=1)
    i.add_argument("--blocks",type=int,default=5)
    i.add_argument("--contact-cutoff",default="3.25",help="Ni-phosphonate-O cutoff in A, or 'auto' (default: 3.25)")
    i.add_argument("--sensitivity-cutoffs",default="3.0,3.25,3.5",help="comma-separated contact cutoffs for bound-fraction sensitivity")
    i.add_argument("--surface-coordination-cutoff",type=float,default=2.8)
    i.add_argument("--persistence-threshold",type=float,default=0.50)
    i.add_argument("--z-min",type=float,default=-5.0)
    i.add_argument("--z-max",type=float,default=60.0)
    i.add_argument("--z-bin-width",type=float,default=0.50)
    i.add_argument("--rdf-bin-width",type=float,default=0.25)
    i.add_argument("--rdf-rmax",type=float,default=15.0)
    i.add_argument("--timestep-fs",type=float)
    i.add_argument("--exchange-min-dwell-frames",type=int,default=1,help="consecutive sole-owner frames required to confirm a site hand-off")
    i.add_argument("--exchange-max-vacancy-gap-frames",type=int,default=3,help="consecutive empty/shared frames after which a site's prior owner is forgotten")
    i.add_argument("--analysis-profile",choices=MODEL_PROFILES,default="classical-ff",help="trajectory model evidence scope; classical-ff withholds kinetic and dipole claims from workbooks")
    i=sub.add_parser("summarize-interface")
    i.add_argument("prepared_root",type=Path)
    i.add_argument("--output",type=Path,required=True)
    i=sub.add_parser("summarize-publication")
    i.add_argument("prepared_root",type=Path)
    i.add_argument("--output",type=Path,required=True)
    i.add_argument("--experimental",type=Path,help="CSV of experimental device Voc/Jsc/FF/PCE results to correlate against MD metrics")
    a=p.parse_args(argv)
    try:
        if a.command=="init-molecule":
            slug=a.name.lower().replace("_","-").replace(" ","-"); folder=molecule_dir(slug); folder.mkdir(parents=True,exist_ok=True)
            manifest=folder/"molecule.toml"
            if not manifest.exists(): manifest.write_text(TEMPLATE.format(name=a.name,slug=slug),encoding="utf-8")
            print(f"Created {manifest}\nPlace the independently generated LigParGen LAMMPS output at:\n{(folder/'ligpargen.lmp').resolve()}")
        elif a.command=="inspect-molecule":
            folder,m=molecule_manifest(a.name); path=folder/m["files"].get("ligpargen","ligpargen.lmp")
            if not path.exists(): raise missing_ligpargen(path,f"nio-md-prep inspect-molecule {a.name}")
            d=parse(path); symbols=elements(d); print(f"{m['molecule']['display_name']}: {d.count('Atoms')} atoms, charge {charge(d):.8f}, elements {' '.join(symbols)}")
        elif a.command=="build": build(a.config,a.output,packed_xyz=a.packed_xyz,packmol_seed=a.packmol_seed,velocity_seed=a.velocity_seed); print(f"Build plan written to {a.output}")
        elif a.command=="refresh-inputs": refresh_inputs(a.config,a.output); print(f"Stage inputs refreshed in {a.output}; simulation data were not modified")
        elif a.command=="validate": validate(a.output,primary_final=a.primary_final); print((a.output/"validation_report.txt").read_text(),end="")
        elif a.command=="prepare-agglomeration":
            from .agglomeration import prepare_agglomeration
            reference_dir = None
            reference_dirs = None
            if a.reference_dir:
                qualified = ["=" in value for value in a.reference_dir]
                if len(a.reference_dir) == 1 and not qualified[0]:
                    reference_dir = Path(a.reference_dir[0])
                else:
                    if not all(qualified):
                        raise ValueError(
                            "mixed references must all use --reference-dir SLUG=DIR"
                        )
                    reference_dirs = {}
                    for value in a.reference_dir:
                        slug, directory = value.split("=", 1)
                        slug, directory = slug.strip(), directory.strip()
                        if not slug or not directory:
                            raise ValueError(
                                "mixed references must use --reference-dir SLUG=DIR"
                            )
                        if slug in reference_dirs:
                            raise ValueError(f"duplicate VASP reference slug: {slug}")
                        reference_dirs[slug] = Path(directory)
            manifest=prepare_agglomeration(
                a.config,
                a.output,
                reference_dir=reference_dir,
                reference_dirs=reference_dirs,
                vasp_template_slug=a.vasp_template_slug,
                packed_xyz=a.packed_xyz,
                structures_only=a.structures_only,
                regenerate_vasp=a.regenerate_vasp,
                regenerate_xtb_launcher=a.regenerate_xtb_launcher,
            )
            if a.regenerate_xtb_launcher:
                print(
                    f"xTB launch files regenerated in {a.output}. Audit with "
                    f"{a.output / 'audit_xtb_runs.py'}; preview the pending-only "
                    f"workq carousel with {a.output / 'launch_xtb.sh'} --dry-run."
                )
                return 0
            manifest_data=json.loads(manifest.read_text(encoding="utf-8"))
            status=manifest_data["status"]
            if status=="PACKMOL_REQUIRED":
                print(f"Agglomeration Packmol inputs written to {a.output}; run packmol in each replica directory and rerun this command with the same output path")
            elif status=="XTB_REQUIRED":
                print(_agglomeration_xtb_pending_message(manifest_data, a.output))
            else:
                print(
                    _agglomeration_completion_message(
                        manifest_data,
                        a.output,
                        structures_only=a.structures_only,
                        regenerated_vasp=a.regenerate_vasp,
                    )
                )
        elif a.command=="audit-agglomerations":
            from .agglomeration_audit import audit_agglomerations
            audit_agglomerations(a.root, a.output_prefix)
        elif a.command=="archive-runs":
            from .archive import create_run_archive
            result=create_run_archive(a.roots,a.output,force=a.force)
            size_mib=result.uncompressed_bytes/(1024*1024)
            print(f"Rerun archive written to {result.path} ({result.file_count} files, {size_mib:.1f} MiB uncompressed)")
        elif a.command=="analyze-coverage":
            from .analysis.coverage import analyze_coverage
            summary=analyze_coverage(
                a.build_directory,
                trajectory=a.trajectory,
                output=a.output,
                grid_spacing=a.grid_spacing,
                radius_scale=a.radius_scale,
                roughness_grid_spacing=a.roughness_grid_spacing,
                near_surface_height=a.near_surface_height,
                last_frames=None if a.last_frames==0 else a.last_frames,
                stride=a.stride,
                blocks=a.blocks,
                exclude_hydrogen=a.exclude_hydrogen,
                timestep_fs=a.timestep_fs,
            )
            print(f"Coverage analysis written to {summary.parent}")
        elif a.command=="summarize-coverage":
            from .analysis.report import create_coverage_workbook
            workbook=create_coverage_workbook(a.prepared_root,a.output)
            print(f"Coverage workbook written to {workbook}")
        elif a.command=="analyze-interface":
            from .analysis.interfacial import analyze_interfacial_structure
            if str(a.contact_cutoff).lower()=="auto":
                contact_cutoff=None
            else:
                try:
                    contact_cutoff=float(a.contact_cutoff)
                except ValueError as exc:
                    raise ValueError(
                        "--contact-cutoff must be a positive number or 'auto'"
                    ) from exc
            summary=analyze_interfacial_structure(
                a.build_directory,
                trajectory=a.trajectory,
                output=a.output,
                last_frames=None if a.last_frames==0 else a.last_frames,
                stride=a.stride,
                blocks=a.blocks,
                contact_cutoff=contact_cutoff,
                sensitivity_cutoffs=tuple(
                    float(value)
                    for value in a.sensitivity_cutoffs.split(",")
                    if value.strip()
                ),
                surface_coordination_cutoff=a.surface_coordination_cutoff,
                persistence_threshold=a.persistence_threshold,
                z_min=a.z_min,
                z_max=a.z_max,
                z_bin_width=a.z_bin_width,
                rdf_bin_width=a.rdf_bin_width,
                rdf_rmax=a.rdf_rmax,
                timestep_fs=a.timestep_fs,
                exchange_min_dwell_frames=a.exchange_min_dwell_frames,
                exchange_max_vacancy_gap_frames=a.exchange_max_vacancy_gap_frames,
                analysis_profile=a.analysis_profile,
            )
            print(f"Interfacial analysis written to {summary.parent}")
        elif a.command=="summarize-interface":
            from .analysis.interfacial_report import create_interfacial_workbook
            workbook=create_interfacial_workbook(a.prepared_root,a.output)
            print(f"Interfacial workbook written to {workbook}")
        elif a.command=="summarize-publication":
            from .analysis.publication_report import create_publication_workbook
            workbook=create_publication_workbook(a.prepared_root,a.output,experimental_csv=a.experimental)
            print(f"Publication workbook written to {workbook}")
        elif a.command=="prepare-lego-stage2":
            from .lego import prepare_lego_stage2
            prepare_lego_stage2(
                a.config,
                a.primary_final,
                a.coverage_map,
                a.output,
                occupancy_threshold=a.occupancy_threshold,
                packing_inset_angstrom=a.packing_inset,
                minimum_gap_fraction=a.minimum_gap_fraction,
                minimum_tunnel_width_angstrom=a.minimum_tunnel_width,
                deposition_clearance_above_stage1_angstrom=(
                    a.deposition_clearance
                ),
                final_deposition_clearance_above_stage1_angstrom=(
                    a.final_deposition_clearance
                ),
                deposition_continuation_steps=a.continuation_steps,
                lateral_confinement=a.lateral_wall,
                wall_strength=a.wall_strength,
                wall_cutoff_angstrom=a.wall_cutoff,
                packed_xyz=a.packed_xyz,
                packmol_seed=a.packmol_seed,
                velocity_seed=a.velocity_seed,
            )
            print(f"Lego-style sequential stage 2 written to {a.output}")
        elif a.command=="prepare-lego2-stage2":
            from .lego2 import prepare_lego2_stage2
            prepare_lego2_stage2(
                a.config,
                a.primary_final,
                a.coverage_map,
                a.output,
                occupancy_threshold=a.occupancy_threshold,
                packing_inset_angstrom=a.packing_inset,
                minimum_void_area_fraction=a.minimum_void_area_fraction,
                minimum_packing_width_angstrom=a.minimum_packing_width,
                minimum_packing_height_angstrom=a.minimum_packing_height,
                deposition_clearance_above_stage1_angstrom=(
                    a.deposition_clearance
                ),
                final_deposition_clearance_above_stage1_angstrom=(
                    a.final_deposition_clearance
                ),
                deposition_continuation_steps=a.continuation_steps,
                packed_xyz=a.packed_xyz,
                packmol_seed=a.packmol_seed,
                velocity_seed=a.velocity_seed,
            )
            print(f"Lego2 2D-void sequential stage 2 written to {a.output}")
        elif a.command=="prepare-lego-continuation":
            from .lego import prepare_lego_continuation
            prepare_lego_continuation(
                a.config,
                a.output,
                additional_drop_angstrom=a.additional_drop,
                continuation_steps=a.steps,
            )
            print(
                "Lego deposition continuation and downstream inputs written "
                f"to {a.output}; simulation data were not modified"
            )
        else:
            build(a.config,a.output,primary_final=a.primary_final,packed_xyz=a.packed_xyz,packmol_seed=a.packmol_seed,velocity_seed=a.velocity_seed)
            print(f"Sequential stage 2 written to {a.output}")
    except (ValueError,FileNotFoundError,FileExistsError,RuntimeError) as e:
        p.exit(2,f"error: {e}\n")
    return 0
