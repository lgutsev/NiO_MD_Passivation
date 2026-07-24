from pathlib import Path
import csv
import json
import shutil

import pytest

from nio_md_prep.analysis.coverage import load_type_elements
from nio_md_prep.analysis.interfacial import (
    _canonical_exposed_ni_sites,
    analyze_interfacial_structure,
)
from nio_md_prep.analysis.interfacial_report import (
    COMPONENT_HEADERS,
    RESULT_HEADERS,
    create_interfacial_workbook,
)
from nio_md_prep.lammps import parse


ATOM_ROWS = [
    (1, 1, 3, 2.0, 2.0, 2.5),
    (2, 1, 2, 2.0, 2.0, 1.8),
    (3, 1, 2, 2.5, 2.0, 2.3),
    (4, 1, 2, 1.5, 2.0, 2.3),
    (5, 1, 1, 2.0, 2.0, 3.5),
    (6, 1, 1, 2.0, 2.0, 6.0),
    (7, 2, 3, 8.0, 8.0, 2.5),
    (8, 2, 2, 8.0, 8.0, 1.8),
    (9, 2, 2, 8.5, 8.0, 2.3),
    (10, 2, 2, 7.5, 8.0, 2.3),
    (11, 2, 1, 8.0, 8.0, 3.5),
    (12, 2, 1, 6.5, 8.0, 6.0),
    (13, 2, 3, 2.0, 8.0, 10.0),
    (14, 2, 2, 2.0, 8.0, 9.3),
    (15, 2, 2, 2.5, 8.0, 9.8),
    (16, 2, 2, 1.5, 8.0, 9.8),
    (17, 2, 1, 2.0, 8.0, 8.5),
    (18, 2, 1, 5.0, 8.0, 6.0),
    (19, 0, 4, 2.0, 2.0, 0.0),
    (20, 0, 4, 8.0, 8.0, 0.0),
    (21, 0, 2, 2.0, 2.0, -2.0),
    (22, 0, 2, 8.0, 8.0, -2.0),
]


def test_canonical_surface_sites_do_not_change_with_thermal_coordinates():
    root = Path(__file__).parents[1]
    path = root / "inputs/surfaces/corrugated-nio-110/surface.lmp"
    data = parse(path)
    type_elements = load_type_elements(path)
    atom_elements = {
        int(row.fields[0]): type_elements[int(row.fields[2])]
        for row in data.sections["Atoms"]
    }
    atom_molecules = {
        int(row.fields[0]): int(row.fields[1])
        for row in data.sections["Atoms"]
    }
    reference_ids, source = _canonical_exposed_ni_sites(
        data, atom_elements, atom_molecules, 2.8
    )
    for index, row in enumerate(data.sections["Atoms"]):
        row.fields[6] = f"{float(row.fields[6]) + 0.75 * ((index % 5) - 2):.8f}"
    distorted_ids, distorted_source = _canonical_exposed_ni_sites(
        data, atom_elements, atom_molecules, 2.8
    )
    assert len(reference_ids) == 600
    assert distorted_ids == reference_ids
    assert distorted_source == source

BONDS = [
    (1, 1, 1, 2),
    (2, 1, 1, 3),
    (3, 1, 1, 4),
    (4, 1, 1, 5),
    (5, 1, 5, 6),
    (6, 1, 7, 8),
    (7, 1, 7, 9),
    (8, 1, 7, 10),
    (9, 1, 7, 11),
    (10, 1, 11, 12),
    (11, 1, 12, 18),
    (12, 1, 18, 17),
    (13, 1, 17, 13),
    (14, 1, 13, 14),
    (15, 1, 13, 15),
    (16, 1, 13, 16),
]


def write_topology(path: Path) -> None:
    atoms = "\n".join(
        f"{atom_id} {molecule} {atom_type} 0.0 {x} {y} {z}"
        for atom_id, molecule, atom_type, x, y, z in ATOM_ROWS
    )
    bonds = "\n".join(
        f"{bond_id} {bond_type} {left} {right}"
        for bond_id, bond_type, left, right in BONDS
    )
    path.write_text(
        f"""LAMMPS data

22 atoms
16 bonds
4 atom types
1 bond types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
-5.0 20.0 zlo zhi

Masses

1 12.011
2 15.999
3 30.974
4 58.6934

Atoms # full

{atoms}

Bonds

{bonds}
""",
        encoding="utf-8",
    )


def dump_frame(step: int, *, second_terminal_bound: bool) -> str:
    rows = []
    for atom_id, molecule, atom_type, x, y, z in ATOM_ROWS:
        if second_terminal_bound and atom_id in {13, 14, 15, 16, 17}:
            y -= 6.0
            z -= 7.5
        rows.append(f"{atom_id} {molecule} {atom_type} {x} {y} {z}")
    return f"""ITEM: TIMESTEP
{step}
ITEM: NUMBER OF ATOMS
22
ITEM: BOX BOUNDS pp pp ff
0 10
0 10
-5 20
ITEM: ATOMS id mol type x y z
{chr(10).join(rows)}
"""


def build_fixture(tmp_path: Path) -> Path:
    build = tmp_path / "prepared" / "mixed-system"
    build.mkdir(parents=True)
    write_topology(build / "topology_output.lmp")
    (build / "assembly_manifest.json").write_text(
        json.dumps(
            {
                "components": [
                    {"component": "primary", "molecule_ids": [1, 1]},
                    {"component": "dcz-4p", "molecule_ids": [2, 2]},
                ]
            }
        ),
        encoding="utf-8",
    )
    (build / "hold-300K.lammpstrj").write_text(
        dump_frame(0, second_terminal_bound=False)
        + dump_frame(1000, second_terminal_bound=True),
        encoding="utf-8",
    )
    return build


def test_anchor_orientation_terminal_rdf_and_density_outputs(tmp_path):
    build = build_fixture(tmp_path)
    output = build / "interface-analysis-hold-300K"
    summary_path = analyze_interfacial_structure(
        build,
        trajectory=Path("hold-300K.lammpstrj"),
        output=output,
        last_frames=2,
        blocks=2,
        contact_cutoff=3.0,
        surface_coordination_cutoff=2.8,
        z_min=-3.0,
        z_max=15.0,
        z_bin_width=1.0,
        rdf_bin_width=0.5,
        rdf_rmax=4.0,
        timestep_fs=0.5,
    )

    summary = json.loads(summary_path.read_text())
    assert summary["method"] == "anchor_resolved_interfacial_structure"
    assert summary["surface_site_count"] == 2
    assert summary["contact_cutoff_method"] == "fixed_common_cutoff"
    assert summary["contact_sensitivity_cutoffs_angstrom"] == [3.0, 3.25, 3.5]
    primary = summary["components"]["primary"]
    dcz = summary["components"]["dcz-4p"]
    assert primary["metrics"]["bound_fraction_percent"]["mean"] == 100.0
    assert dcz["terminal_state_populations_percent"]["1"]["mean"] == 50.0
    assert dcz["terminal_state_populations_percent"]["2"]["mean"] == 50.0
    assert dcz["persistent_bound_molecules_percent"] == 100.0
    assert set(dcz["contact_cutoff_sensitivity"]) == {"3", "3.25", "3.5"}
    assert "plane_normal_tilt_degrees" in dcz["metrics"]
    assert "anchor_axis_tilt_degrees" in dcz["metrics"]
    assert "bound_phosphorus_height_angstrom" in dcz["metrics"]
    assert summary["site_ownership_metrics"]["shared_percent"]["mean"] == 25.0
    assert summary["site_ownership_metrics"]["empty_percent"]["mean"] == 0.0
    for name in (
        "interface_timeseries.csv",
        "contact_distance_histogram.csv",
        "z_density_profiles.csv",
        "lateral_rdf.csv",
    ):
        assert (output / name).is_file()
    with (output / "lateral_rdf.csv").open() as handle:
        pairs = {row["component_pair"] for row in csv.DictReader(handle)}
    assert {"primary--primary", "primary--dcz_4p", "dcz_4p--dcz_4p"} <= pairs


def test_separate_workbook_contains_normalized_sheets(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    build = build_fixture(tmp_path)
    analyze_interfacial_structure(
        build,
        trajectory=Path("hold-300K.lammpstrj"),
        output=build / "interface-analysis-hold-300K",
        contact_cutoff=3.0,
        surface_coordination_cutoff=2.8,
        z_min=-3.0,
        z_max=15.0,
        z_bin_width=1.0,
        rdf_bin_width=0.5,
        rdf_rmax=4.0,
        timestep_fs=0.5,
    )
    hold_output = build / "interface-analysis-hold-300K"
    relax_output = build / "interface-analysis-relax-300K"
    shutil.copytree(hold_output, relax_output)
    relax_summary_path = relax_output / "interface_summary.json"
    relax_summary = json.loads(relax_summary_path.read_text())
    relax_summary["trajectory"] = str(build / "relax-300K.lammpstrj")
    relax_summary_path.write_text(json.dumps(relax_summary), encoding="utf-8")
    output = build.parent / "interface_structure_summary.xlsx"
    assert create_interfacial_workbook(build.parent, output) == output

    workbook = openpyxl.load_workbook(output)
    assert workbook.sheetnames == [
        "Results",
        "Components",
        "Terminal States",
        "Hold vs Relax",
        "Z Profiles",
        "Lateral RDF",
        "Cutoff Sensitivity",
        "Methods",
    ]
    assert [cell.value for cell in workbook["Results"][1]] == RESULT_HEADERS
    assert [cell.value for cell in workbook["Components"][1]] == COMPONENT_HEADERS
    assert workbook["Results"]["M2"].value == pytest.approx(100.0)
    assert workbook["Results"]["AC2"].value == "OK"
    assert "InterfaceResultsTable" in workbook["Results"].tables
    assert "InterfaceComponentsTable" in workbook["Components"].tables
    assert "CutoffSensitivityTable" in workbook["Cutoff Sensitivity"].tables
    assert workbook["Hold vs Relax"].max_row == 2
    assert workbook["Hold vs Relax"]["A2"].value == "mixed-system"
