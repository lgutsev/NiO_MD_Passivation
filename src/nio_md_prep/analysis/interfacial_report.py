"""Consolidated Excel reporting for anchor-resolved interface analyses."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import json
import math

from .report import _openpyxl, _style_table_sheet, _temperature


RESULT_HEADERS = [
    "Run",
    "System",
    "Temperature (K)",
    "Stage",
    "Frames Analyzed",
    "Contact Cutoff (Å)",
    "Cutoff Method",
    "Exposed Ni Sites",
    "Empty Sites (%)",
    "Shared Sites (%)",
    "Primary Component",
    "Primary-only Sites (%)",
    "Primary Bound (%)",
    "Primary Persistent (%)",
    "Primary Tilt (deg)",
    "Primary P2",
    "Primary P Height (Å)",
    "Primary Core Height (Å)",
    "Secondary Component",
    "Secondary-only Sites (%)",
    "Secondary Bound (%)",
    "Secondary Persistent (%)",
    "Secondary Tilt (deg)",
    "Secondary P2",
    "Secondary P Height (Å)",
    "Secondary Core Height (Å)",
    "Site Balance QC (%)",
    "Terminal Balance QC (%)",
    "Status",
    "Trajectory",
    "Summary JSON",
    "Run Directory",
    "Packmol Seed",
    "Velocity Seed",
]

COMPONENT_HEADERS = [
    "Run",
    "System",
    "Temperature (K)",
    "Stage",
    "Component",
    "Molecule ID Start",
    "Molecule ID End",
    "Molecule Count",
    "Anchor Terminals / Molecule",
    "Bound Fraction (%)",
    "Block SEM (%)",
    "Persistent Bound (%)",
    "Mean Contact Occupancy (%)",
    "Anchored Terminals / Molecule",
    "Tilt (deg)",
    "Orientational P2",
    "P Height (Å)",
    "Core Height (Å)",
    "Core-P Height (Å)",
    "Nearest Neighbor (Å)",
    "Residence Events",
    "Mean Residence (frames)",
    "Max Residence (frames)",
    "Mean Residence (ps)",
    "Max Residence (ps)",
    "Plane-normal Tilt (deg)",
    "Plane Alignment P2",
    "Anchor-axis Tilt (deg)",
    "Anchor Δz (Å)",
    "Bound P Height (Å)",
    "Unbound P Height (Å)",
    "Bound Core Height (Å)",
    "Unbound Core Height (Å)",
    "Bound Tilt (deg)",
    "Unbound Tilt (deg)",
    "Run Directory",
    "Packmol Seed",
    "Velocity Seed",
]

TERMINAL_HEADERS = [
    "Run",
    "System",
    "Temperature (K)",
    "Stage",
    "Component",
    "Anchored Terminal Count",
    "Population (%)",
    "Block SEM (%)",
    "Frame SD (%)",
]

COMPARISON_HEADERS = [
    "Run Directory",
    "Temperature (K)",
    "Compressed Stage",
    "Relaxed Stage",
    "Δ Empty Sites (%)",
    "Δ Primary Bound (%)",
    "Δ Secondary Bound (%)",
    "Δ Primary Tilt (deg)",
    "Δ Secondary Tilt (deg)",
    "Δ Primary P2",
    "Δ Secondary P2",
]

Z_PROFILE_HEADERS = [
    "Run",
    "System",
    "Temperature (K)",
    "Stage",
    "Component",
    "Local Height (Å)",
    "Heavy-atom Density (Å^-3)",
    "Phosphorus Density (Å^-3)",
]

RDF_HEADERS = [
    "Run",
    "System",
    "Temperature (K)",
    "Stage",
    "Component Pair",
    "Radius (Å)",
    "g2D(r)",
    "Observed Pairs",
    "Ideal Pairs",
]

CUTOFF_HEADERS = [
    "Run",
    "System",
    "Temperature (K)",
    "Stage",
    "Component",
    "Contact Cutoff (Å)",
    "Bound Fraction (%)",
    "Block SEM (%)",
    "Anchored Terminals / Molecule",
    "Run Directory",
]


def _metric(component: dict[str, Any] | None, name: str, field: str = "mean"):
    if component is None:
        return None
    return component.get("metrics", {}).get(name, {}).get(field)


def _component_rows(summary: dict[str, Any]) -> list[dict[str, Any]]:
    components = []
    for name, values in summary.get("components", {}).items():
        molecule_range = values.get("molecule_id_range", [None, None])
        components.append(
            {
                "name": name,
                "lo": molecule_range[0] if len(molecule_range) else None,
                "hi": molecule_range[1] if len(molecule_range) > 1 else None,
                "values": values,
            }
        )
    components.sort(
        key=lambda row: (
            math.inf if row["lo"] is None else int(row["lo"]),
            row["name"],
        )
    )
    return components


def collect_interfacial_results(
    prepared_root: Path,
) -> tuple[
    list[dict],
    list[dict],
    list[dict],
    list[dict],
    list[dict],
    list[dict],
]:
    """Collect normalized result, component, terminal, z-profile, and RDF rows."""
    prepared_root = Path(prepared_root)
    paths = sorted(
        set(
            prepared_root.rglob(
                "interface-analysis-hold-*/interface_summary.json"
            )
        )
        | set(
            prepared_root.rglob(
                "interface-analysis-relax-*/interface_summary.json"
            )
        )
    )
    if not paths:
        raise FileNotFoundError(
            f"{prepared_root}: no compressed-hold or relaxed-hold "
            "interfacial summaries were found"
        )

    results = []
    components_out = []
    terminals_out = []
    z_out = []
    rdf_out = []
    cutoff_out = []
    for summary_path in paths:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        run_directory_path = summary_path.parents[1]
        system = run_directory_path.name
        run_directory = str(run_directory_path.relative_to(prepared_root))
        stage = summary_path.parent.name.removeprefix("interface-analysis-")
        trajectory = str(summary.get("trajectory", ""))
        temperature = _temperature(stage, trajectory)
        run = (
            f"{run_directory} | {stage}"
            if temperature is None
            else f"{run_directory} | {temperature} K | {stage.split('-', 1)[0]}"
        )
        random_seeds = summary.get("random_seeds", {})
        component_rows = _component_rows(summary)
        primary = component_rows[0] if component_rows else None
        secondary = component_rows[1] if len(component_rows) > 1 else None
        site_metrics = summary.get("site_ownership_metrics", {})
        only_sum = 0.0
        for component in component_rows:
            key = component["values"]["key"]
            only_sum += float(
                site_metrics.get(f"{key}_only_percent", {}).get("mean", 0.0)
            )
        empty = float(site_metrics.get("empty_percent", {}).get("mean", 0.0))
        shared = float(site_metrics.get("shared_percent", {}).get("mean", 0.0))
        site_balance = empty + shared + only_sum - 100.0
        terminal_balance = 0.0
        for component in component_rows:
            populations = component["values"].get(
                "terminal_state_populations_percent", {}
            )
            if populations:
                terminal_balance = max(
                    terminal_balance,
                    abs(
                        sum(
                            float(values.get("mean", 0.0))
                            for values in populations.values()
                        )
                        - 100.0
                    ),
                )
        status = (
            "OK"
            if abs(site_balance) <= 0.05 and terminal_balance <= 0.05
            else "CHECK"
        )
        try:
            relative_summary = str(summary_path.relative_to(prepared_root.parent))
        except ValueError:
            relative_summary = str(summary_path)

        def comp_value(component, key):
            return component["values"].get(key) if component else None

        def comp_metric(component, key):
            return _metric(component["values"], key) if component else None

        result = {
            "run": run,
            "system": system,
            "temperature": temperature,
            "stage": stage,
            "frames": int(summary["frames_analyzed"]),
            "cutoff": float(summary["contact_cutoff_angstrom"]),
            "cutoff_method": str(summary["contact_cutoff_method"]),
            "surface_sites": int(summary["surface_site_count"]),
            "empty": empty,
            "shared": shared,
            "primary_name": primary["name"] if primary else None,
            "primary_only": (
                site_metrics.get(
                    f"{primary['values']['key']}_only_percent", {}
                ).get("mean")
                if primary
                else None
            ),
            "primary_bound": comp_metric(primary, "bound_fraction_percent"),
            "primary_unbound_anchor_density": comp_metric(
                primary, "unbound_anchor_density_per_nm2"
            ),
            "primary_persistent": comp_value(
                primary, "persistent_bound_molecules_percent"
            ),
            "primary_tilt": comp_metric(primary, "tilt_degrees"),
            "primary_p2": comp_metric(primary, "orientational_P2"),
            "primary_p_height": comp_metric(
                primary, "phosphorus_height_angstrom"
            ),
            "primary_core_height": comp_metric(
                primary, "core_height_angstrom"
            ),
            "secondary_name": secondary["name"] if secondary else None,
            "secondary_only": (
                site_metrics.get(
                    f"{secondary['values']['key']}_only_percent", {}
                ).get("mean")
                if secondary
                else None
            ),
            "secondary_bound": comp_metric(secondary, "bound_fraction_percent"),
            "secondary_unbound_anchor_density": comp_metric(
                secondary, "unbound_anchor_density_per_nm2"
            ),
            "secondary_persistent": comp_value(
                secondary, "persistent_bound_molecules_percent"
            ),
            "secondary_tilt": comp_metric(secondary, "tilt_degrees"),
            "secondary_p2": comp_metric(secondary, "orientational_P2"),
            "secondary_p_height": comp_metric(
                secondary, "phosphorus_height_angstrom"
            ),
            "secondary_core_height": comp_metric(
                secondary, "core_height_angstrom"
            ),
            "site_balance": site_balance,
            "terminal_balance": terminal_balance,
            "status": status,
            "trajectory": trajectory,
            "summary_path": relative_summary,
            "run_directory": run_directory,
            "packmol_seed": random_seeds.get("packmol_seed"),
            "velocity_seed": random_seeds.get("velocity_seed"),
        }
        results.append(result)

        for component in component_rows:
            values = component["values"]
            metrics = values.get("metrics", {})
            anchor_counts = values.get("anchor_terminals_per_molecule", [])
            components_out.append(
                {
                    "run": run,
                    "system": system,
                    "temperature": temperature,
                    "stage": stage,
                    "name": component["name"],
                    "lo": component["lo"],
                    "hi": component["hi"],
                    "count": values.get("molecule_count"),
                    "anchor_terminals": ", ".join(map(str, anchor_counts)),
                    "bound": metrics["bound_fraction_percent"]["mean"],
                    "bound_sem": metrics["bound_fraction_percent"]["block_sem"],
                    "persistent": values.get(
                        "persistent_bound_molecules_percent"
                    ),
                    "occupancy": values.get(
                        "mean_molecule_contact_occupancy_percent"
                    ),
                    "anchored_terminals": metrics[
                        "mean_anchored_terminals"
                    ]["mean"],
                    "tilt": metrics["tilt_degrees"]["mean"],
                    "p2": metrics["orientational_P2"]["mean"],
                    "p_height": metrics["phosphorus_height_angstrom"]["mean"],
                    "core_height": metrics["core_height_angstrom"]["mean"],
                    "relative_height": metrics[
                        "core_minus_phosphorus_height_angstrom"
                    ]["mean"],
                    "nearest_neighbor": metrics[
                        "nearest_neighbor_distance_angstrom"
                    ]["mean"],
                    "residence_events": values.get("residence_event_count"),
                    "mean_residence_frames": values.get(
                        "mean_residence_frames"
                    ),
                    "max_residence_frames": values.get("max_residence_frames"),
                    "mean_residence_ps": values.get("mean_residence_ps"),
                    "max_residence_ps": values.get("max_residence_ps"),
                    "plane_tilt": _metric(
                        values, "plane_normal_tilt_degrees"
                    ),
                    "plane_p2": _metric(values, "plane_alignment_P2"),
                    "anchor_axis_tilt": _metric(
                        values, "anchor_axis_tilt_degrees"
                    ),
                    "anchor_vertical_separation": _metric(
                        values, "anchor_vertical_separation_angstrom"
                    ),
                    "bound_p_height": _metric(
                        values, "bound_phosphorus_height_angstrom"
                    ),
                    "unbound_p_height": _metric(
                        values, "unbound_phosphorus_height_angstrom"
                    ),
                    "bound_core_height": _metric(
                        values, "bound_core_height_angstrom"
                    ),
                    "unbound_core_height": _metric(
                        values, "unbound_core_height_angstrom"
                    ),
                    "bound_tilt": _metric(values, "bound_tilt_degrees"),
                    "unbound_tilt": _metric(values, "unbound_tilt_degrees"),
                    "run_directory": run_directory,
                    "packmol_seed": random_seeds.get("packmol_seed"),
                    "velocity_seed": random_seeds.get("velocity_seed"),
                }
            )
            for cutoff, cutoff_values in sorted(
                values.get("contact_cutoff_sensitivity", {}).items(),
                key=lambda item: float(item[0]),
            ):
                bound_values = cutoff_values.get(
                    "bound_fraction_percent", {}
                )
                anchored_values = cutoff_values.get(
                    "mean_anchored_terminals", {}
                )
                cutoff_out.append(
                    {
                        "run": run,
                        "system": system,
                        "temperature": temperature,
                        "stage": stage,
                        "component": component["name"],
                        "cutoff": float(cutoff),
                        "bound": bound_values.get("mean"),
                        "sem": bound_values.get("block_sem"),
                        "anchored": anchored_values.get("mean"),
                        "run_directory": run_directory,
                    }
                )
            for count, population in sorted(
                values.get("terminal_state_populations_percent", {}).items(),
                key=lambda item: int(item[0]),
            ):
                terminals_out.append(
                    {
                        "run": run,
                        "system": system,
                        "temperature": temperature,
                        "stage": stage,
                        "component": component["name"],
                        "count": int(count),
                        "population": population["mean"],
                        "sem": population["block_sem"],
                        "std": population["frame_std"],
                    }
                )

        z_path = summary_path.parent / "z_density_profiles.csv"
        if z_path.is_file():
            with z_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    z_out.append(
                        {
                            "run": run,
                            "system": system,
                            "temperature": temperature,
                            "stage": stage,
                            "component": row["component"],
                            "height": float(row["local_height_angstrom"]),
                            "heavy": float(
                                row["heavy_atom_number_density_angstrom-3"]
                            ),
                            "phosphorus": float(
                                row["phosphorus_number_density_angstrom-3"]
                            ),
                        }
                    )
        rdf_path = summary_path.parent / "lateral_rdf.csv"
        if rdf_path.is_file():
            with rdf_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    rdf_out.append(
                        {
                            "run": run,
                            "system": system,
                            "temperature": temperature,
                            "stage": stage,
                            "pair": row["component_pair"],
                            "radius": float(row["radius_angstrom"]),
                            "g": float(row["g_2d"]),
                            "observed": float(row["observed_pairs"]),
                            "ideal": float(row["ideal_pairs"]),
                        }
                    )

    stage_order = {"hold": 0, "relax": 1}
    results.sort(
        key=lambda row: (
            row["system"],
            math.inf if row["temperature"] is None else row["temperature"],
            stage_order.get(row["stage"].split("-", 1)[0], 9),
            row["stage"],
        )
    )
    components_out.sort(
        key=lambda row: (
            row["system"],
            math.inf if row["temperature"] is None else row["temperature"],
            stage_order.get(row["stage"].split("-", 1)[0], 9),
            math.inf if row["lo"] is None else row["lo"],
        )
    )
    return results, components_out, terminals_out, z_out, rdf_out, cutoff_out


def _stage_comparisons(results: list[dict]) -> list[list[Any]]:
    grouped: dict[tuple[str, int | None], dict[str, dict]] = {}
    for row in results:
        kind = row["stage"].split("-", 1)[0]
        grouped.setdefault(
            (row["run_directory"], row["temperature"]), {}
        )[kind] = row
    comparisons = []
    for (run_directory, temperature), stages in sorted(grouped.items()):
        if "hold" not in stages or "relax" not in stages:
            continue
        compressed = stages["hold"]
        relaxed = stages["relax"]

        def delta(name):
            left = compressed.get(name)
            right = relaxed.get(name)
            return None if left is None or right is None else right - left

        comparisons.append(
            [
                run_directory,
                temperature,
                compressed["stage"],
                relaxed["stage"],
                delta("empty"),
                delta("primary_bound"),
                delta("secondary_bound"),
                delta("primary_tilt"),
                delta("secondary_tilt"),
                delta("primary_p2"),
                delta("secondary_p2"),
            ]
        )
    return comparisons


def create_interfacial_workbook(prepared_root: Path, output: Path) -> Path:
    """Create the separate anchor-resolved structural-analysis workbook."""
    api = _openpyxl()
    Workbook = api["Workbook"]
    BarChart = api["BarChart"]
    Reference = api["Reference"]
    FormulaRule = api["FormulaRule"]
    PatternFill = api["PatternFill"]
    Font = api["Font"]

    results, components, terminals, z_rows, rdf_rows, cutoff_rows = (
        collect_interfacial_results(prepared_root)
    )
    workbook = Workbook()
    workbook.properties.creator = "nio-md-prep"
    workbook.properties.title = "NiO phosphonate interfacial structure summary"
    results_sheet = workbook.active
    results_sheet.title = "Results"
    component_sheet = workbook.create_sheet("Components")
    terminal_sheet = workbook.create_sheet("Terminal States")
    comparison_sheet = workbook.create_sheet("Hold vs Relax")
    z_sheet = workbook.create_sheet("Z Profiles")
    rdf_sheet = workbook.create_sheet("Lateral RDF")
    cutoff_sheet = workbook.create_sheet("Cutoff Sensitivity")
    methods_sheet = workbook.create_sheet("Methods")

    result_rows = [
        [
            row["run"],
            row["system"],
            row["temperature"],
            row["stage"],
            row["frames"],
            row["cutoff"],
            row["cutoff_method"],
            row["surface_sites"],
            row["empty"],
            row["shared"],
            row["primary_name"],
            row["primary_only"],
            row["primary_bound"],
            row["primary_persistent"],
            row["primary_tilt"],
            row["primary_p2"],
            row["primary_p_height"],
            row["primary_core_height"],
            row["secondary_name"],
            row["secondary_only"],
            row["secondary_bound"],
            row["secondary_persistent"],
            row["secondary_tilt"],
            row["secondary_p2"],
            row["secondary_p_height"],
            row["secondary_core_height"],
            row["site_balance"],
            row["terminal_balance"],
            row["status"],
            row["trajectory"],
            row["summary_path"],
            row["run_directory"],
            row["packmol_seed"],
            row["velocity_seed"],
        ]
        for row in results
    ]
    _style_table_sheet(
        results_sheet,
        RESULT_HEADERS,
        result_rows,
        "InterfaceResultsTable",
        api,
    )
    for column in (
        "F",
        "I",
        "J",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "T",
        "U",
        "V",
        "W",
        "X",
        "Y",
        "Z",
        "AA",
        "AB",
    ):
        for cell in results_sheet[column][1:]:
            cell.number_format = "0.00"
    for column in ("C", "E", "H"):
        for cell in results_sheet[column][1:]:
            cell.number_format = "#,##0"
    for index, width in enumerate(
        [
            46,
            34,
            14,
            16,
            16,
            18,
            35,
            17,
            17,
            18,
            24,
            22,
            19,
            21,
            18,
            13,
            19,
            21,
            24,
            24,
            21,
            23,
            20,
            15,
            21,
            23,
            19,
            22,
            11,
            48,
            48,
            44,
            17,
            17,
        ],
        1,
    ):
        results_sheet.column_dimensions[
            results_sheet.cell(1, index).column_letter
        ].width = width
    if results:
        last = len(results) + 1
        results_sheet.conditional_formatting.add(
            f"AC2:AC{last}",
            FormulaRule(
                formula=['$AC2="OK"'],
                fill=PatternFill("solid", fgColor="C6EFCE"),
                font=Font(color="006100"),
            ),
        )
        results_sheet.conditional_formatting.add(
            f"AC2:AC{last}",
            FormulaRule(
                formula=['$AC2="CHECK"'],
                fill=PatternFill("solid", fgColor="FFC7CE"),
                font=Font(color="9C0006"),
            ),
        )
        chart = BarChart()
        chart.type = "bar"
        chart.style = 10
        chart.title = "Primary bound fraction by trajectory"
        chart.x_axis.title = "Bound molecules (%)"
        chart.y_axis.title = "System | stage"
        chart.height = max(7.0, min(16.0, 2.5 + 0.45 * len(results)))
        chart.width = 17.0
        chart.add_data(
            Reference(results_sheet, min_col=13, min_row=1, max_row=last),
            titles_from_data=True,
        )
        chart.set_categories(
            Reference(results_sheet, min_col=1, min_row=2, max_row=last)
        )
        chart.legend = None
        results_sheet.add_chart(chart, "AF2")

    component_rows = [
        [
            row["run"],
            row["system"],
            row["temperature"],
            row["stage"],
            row["name"],
            row["lo"],
            row["hi"],
            row["count"],
            row["anchor_terminals"],
            row["bound"],
            row["bound_sem"],
            row["persistent"],
            row["occupancy"],
            row["anchored_terminals"],
            row["tilt"],
            row["p2"],
            row["p_height"],
            row["core_height"],
            row["relative_height"],
            row["nearest_neighbor"],
            row["residence_events"],
            row["mean_residence_frames"],
            row["max_residence_frames"],
            row["mean_residence_ps"],
            row["max_residence_ps"],
            row["plane_tilt"],
            row["plane_p2"],
            row["anchor_axis_tilt"],
            row["anchor_vertical_separation"],
            row["bound_p_height"],
            row["unbound_p_height"],
            row["bound_core_height"],
            row["unbound_core_height"],
            row["bound_tilt"],
            row["unbound_tilt"],
            row["run_directory"],
            row["packmol_seed"],
            row["velocity_seed"],
        ]
        for row in components
    ]
    _style_table_sheet(
        component_sheet,
        COMPONENT_HEADERS,
        component_rows,
        "InterfaceComponentsTable",
        api,
    )
    for column in ("J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "V", "X", "Y"):
        for cell in component_sheet[column][1:]:
            cell.number_format = "0.00"
    for column in ("C", "F", "G", "H", "U", "W"):
        for cell in component_sheet[column][1:]:
            cell.number_format = "#,##0"
    for index in range(1, len(COMPONENT_HEADERS) + 1):
        component_sheet.column_dimensions[
            component_sheet.cell(1, index).column_letter
        ].width = 19 if index > 5 else 28
    component_sheet.column_dimensions["AJ"].width = 44

    terminal_rows = [
        [
            row["run"],
            row["system"],
            row["temperature"],
            row["stage"],
            row["component"],
            row["count"],
            row["population"],
            row["sem"],
            row["std"],
        ]
        for row in terminals
    ]
    _style_table_sheet(
        terminal_sheet,
        TERMINAL_HEADERS,
        terminal_rows,
        "TerminalStatesTable",
        api,
    )
    for column in ("G", "H", "I"):
        for cell in terminal_sheet[column][1:]:
            cell.number_format = "0.00"
    for column, width in zip(
        "ABCDEFGHI", [46, 34, 14, 16, 25, 24, 18, 16, 16]
    ):
        terminal_sheet.column_dimensions[column].width = width

    comparison_rows = _stage_comparisons(results)
    _style_table_sheet(
        comparison_sheet,
        COMPARISON_HEADERS,
        comparison_rows,
        "HoldRelaxComparisonTable",
        api,
    )
    for column in "EFGHIJK":
        for cell in comparison_sheet[column][1:]:
            cell.number_format = "0.00"
    for column, width in zip(
        "ABCDEFGHIJK", [34, 14, 18, 18, 18, 22, 24, 21, 23, 17, 19]
    ):
        comparison_sheet.column_dimensions[column].width = width

    z_values = [
        [
            row["run"],
            row["system"],
            row["temperature"],
            row["stage"],
            row["component"],
            row["height"],
            row["heavy"],
            row["phosphorus"],
        ]
        for row in z_rows
    ]
    _style_table_sheet(
        z_sheet,
        Z_PROFILE_HEADERS,
        z_values,
        "ZProfilesTable",
        api,
    )
    for column in ("F", "G", "H"):
        for cell in z_sheet[column][1:]:
            cell.number_format = "0.0000"
    for column, width in zip(
        "ABCDEFGH", [46, 34, 14, 16, 25, 18, 27, 27]
    ):
        z_sheet.column_dimensions[column].width = width

    rdf_values = [
        [
            row["run"],
            row["system"],
            row["temperature"],
            row["stage"],
            row["pair"],
            row["radius"],
            row["g"],
            row["observed"],
            row["ideal"],
        ]
        for row in rdf_rows
    ]
    _style_table_sheet(
        rdf_sheet,
        RDF_HEADERS,
        rdf_values,
        "LateralRDFTable",
        api,
    )
    for column in ("F", "G", "H", "I"):
        for cell in rdf_sheet[column][1:]:
            cell.number_format = "0.0000"
    for column, width in zip(
        "ABCDEFGHI", [46, 34, 14, 16, 30, 16, 16, 18, 18]
    ):
        rdf_sheet.column_dimensions[column].width = width

    cutoff_values = [
        [
            row["run"],
            row["system"],
            row["temperature"],
            row["stage"],
            row["component"],
            row["cutoff"],
            row["bound"],
            row["sem"],
            row["anchored"],
            row["run_directory"],
        ]
        for row in cutoff_rows
    ]
    _style_table_sheet(
        cutoff_sheet,
        CUTOFF_HEADERS,
        cutoff_values,
        "CutoffSensitivityTable",
        api,
    )
    for column in ("F", "G", "H", "I"):
        for cell in cutoff_sheet[column][1:]:
            cell.number_format = "0.00"
    for column, width in zip(
        "ABCDEFGHIJ", [46, 34, 14, 16, 25, 20, 21, 16, 29, 44]
    ):
        cutoff_sheet.column_dimensions[column].width = width

    methods_sheet.sheet_view.showGridLines = False
    methods = [
        ["Field", "Value"],
        [
            "Workbook purpose",
            "Anchor-resolved structure, ordering, and wall-relaxation comparisons; separate from projected coverage.",
        ],
        [
            "Bound molecule",
            "At least one phosphonate terminal has an O atom within the recorded cutoff of an exposed Ni site.",
        ],
        [
            "Contact cutoff",
            "A common 3.25 A cutoff is the default. The Cutoff Sensitivity sheet reports bound fractions at 3.0, 3.25, and 3.5 A; auto-inference remains available diagnostically.",
        ],
        [
            "Persistent molecule",
            "Bound in at least the recorded persistence-threshold fraction of analyzed frames.",
        ],
        [
            "Surface sites",
            "Canonical exposed-Ni atom identities selected once from the authoritative pristine corrugated NiO slab and mapped into every assembled topology.",
        ],
        [
            "Tilt",
            "Angle between the phosphonate-P-to-graph-distant-core vector and the surface normal; 0° upright, 90° flat.",
        ],
        [
            "Molecular plane",
            "Plane-normal angle from graph-distant core atoms; 0° means the molecular plane is parallel to the surface. The anchor-axis fields apply to multi-phosphonate molecules such as DCZ-4P.",
        ],
        [
            "Bound/unbound structure",
            "Phosphorus height, core height, and tilt are also reported separately for molecules classified as bound and unbound at the common cutoff.",
        ],
        [
            "Orientational P2",
            "(3 cos²(theta)-1)/2; 1 upright, -0.5 flat.",
        ],
        [
            "Heights",
            "Height relative to the nearest exposed Ni site in periodic x/y, accommodating the corrugated surface.",
        ],
        [
            "Lateral RDF",
            "Two-dimensional molecular-core pair distribution normalized by instantaneous lateral area.",
        ],
        [
            "Uncertainty",
            "Block SEM describes time-series sampling within a trajectory; it is not independent-seed uncertainty.",
        ],
        ["Source root", str(Path(prepared_root).resolve())],
        ["Runs included", len(results)],
    ]
    for row in methods:
        methods_sheet.append(row)
    for cell in methods_sheet[1]:
        cell.fill = PatternFill("solid", fgColor="17365D")
        cell.font = Font(color="FFFFFF", bold=True)
    methods_sheet.column_dimensions["A"].width = 25
    methods_sheet.column_dimensions["B"].width = 115
    methods_sheet.freeze_panes = "A2"
    for row_index, row in enumerate(methods_sheet.iter_rows(), start=1):
        row[1].alignment = api["Alignment"](wrap_text=True, vertical="top")
        if row_index > 1:
            value = str(row[1].value or "")
            wrapped_lines = max(1, (len(value) + 104) // 105)
            methods_sheet.row_dimensions[row_index].height = max(
                18, 15 * wrapped_lines
            )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output)
    return output
