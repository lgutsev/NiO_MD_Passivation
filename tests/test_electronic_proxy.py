from pathlib import Path
import csv
import json
import math

import numpy as np
import pytest

from nio_md_prep.analysis.coverage import iter_dump_frames
from nio_md_prep.analysis.electronic_proxy import (
    PATCH_CSV,
    _frame_state,
    _prepare,
    analyze_electronic_proxy,
    build_facet_families,
    build_patch_windows,
    classify_exposed_sites,
    load_reference_dipoles,
    parse_aliases,
    select_patch_slots,
    summarize_electronic_gate,
)
from nio_md_prep.analysis.interfacial import (
    _DIPOLE_DENSITY_TO_VOLTS,
    load_interfacial_topology,
)
from nio_md_prep.cli import main


D = 2.085
NX, NY = 8, 20
LX, LY = NX * D, NY * D
S = math.sqrt(0.5)
N_MINUS = np.array([0.0, -S, S])
N_PLUS = np.array([0.0, S, S])
X = np.array([1.0, 0.0, 0.0])

# Ligand template in its own frame (e1, e2, e3); e3 is the molecular axis.
TEMPLATE = [
    ("P", (0.0, 0.0, 0.0), -0.5),
    ("O", (0.0, 0.0, -1.2), 0.0),
    ("O", (1.5, 0.0, 0.0), 0.0),
    ("O", (-1.5, 0.0, 0.0), 0.0),
    ("H", (2.5, 0.0, 0.0), 0.0),
    ("C", (0.0, 0.0, 1.8), 0.0),
    ("C", (0.0, 0.0, 3.3), 0.0),
    ("C", (0.0, 0.0, 4.8), 0.5),
]
TEMPLATE_BONDS = [(0, 1), (0, 2), (0, 3), (2, 4), (0, 5), (5, 6), (6, 7)]
TYPES = {"Ni": 1, "Os": 2, "C": 3, "O": 4, "P": 5, "H": 6}
MASSES = {1: 58.6934, 2: 15.999, 3: 12.011, 4: 15.999, 5: 30.974, 6: 1.008}


def top_plane(j: int) -> int:
    """1 groove row, 6 rising steps, 7 ridge rows, 6 falling steps."""
    if j == 0:
        return 3
    if j <= 6:
        return 3 + j
    if j <= 13:
        return 10
    return 3 + (20 - j)


def surface_atoms():
    atoms = []
    for i in range(NX):
        for j in range(NY):
            for k in range(top_plane(j) + 1):
                element = "Ni" if (i + j + k) % 2 == 0 else "Os"
                atoms.append((element, (i * D, j * D, k * D), 2.0 if element == "Ni" else -2.0))
    return atoms


def frame_vectors(n):
    t1 = X
    t2 = np.cross(n, t1)
    return t1, t2, n


def place(origin, e1, e2, e3):
    return [
        (element, tuple(origin + a * e1 + b * e2 + c * e3), q)
        for element, (a, b, c), q in TEMPLATE
    ]


def system():
    surface = surface_atoms()
    # 110-y middle row (j=4, k=7); 110+y middle row (j=16, k=6).
    site_a = np.array([3 * D, 4 * D, 7 * D])
    # B: lifted 8 A above row j=6 of the same facet, lying along x and
    # crossing the periodic x boundary.
    site_b = np.array([7 * D, 6 * D, 9 * D])
    t1, t2, n = frame_vectors(N_MINUS)
    molecule_a = place(site_a + 3.2 * n, t1, t2, n)
    molecule_b = place(site_b + 8.0 * n, t2, n, t1)
    m1, m2, m3 = frame_vectors(N_PLUS)
    site_c1 = np.array([2 * D, 16 * D, 6 * D])
    site_c2 = np.array([4 * D, 16 * D, 6 * D])
    p1, p2 = site_c1 + 3.2 * m3, site_c2 + 3.2 * m3
    molecule_c = [
        ("P", tuple(p1), 0.0),
        ("O", tuple(p1 - 1.2 * m3), 0.0),
        ("O", tuple(p1 + 1.4 * m2), 0.0),
        ("O", tuple(p1 - 1.4 * m2), 0.0),
        ("P", tuple(p2), 0.0),
        ("O", tuple(p2 - 1.2 * m3), 0.0),
        ("O", tuple(p2 + 1.4 * m2), 0.0),
        ("O", tuple(p2 - 1.4 * m2), 0.0),
        ("C", tuple(p1 + 1.5 * m3), 0.0),
        ("C", tuple(0.5 * (p1 + p2) + 2.5 * m3), 0.0),
        ("C", tuple(p2 + 1.5 * m3), 0.0),
    ]
    bonds_c = [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7), (0, 8), (8, 9), (9, 10), (10, 4)]
    return surface, [(molecule_a, TEMPLATE_BONDS), (molecule_b, TEMPLATE_BONDS), (molecule_c, bonds_c)]


def write_build(root: Path, *, primary_name: str = "stage1_primary") -> Path:
    build = root / "prepared" / "fixture"
    build.mkdir(parents=True)
    surface, molecules = system()
    rows, bonds = [], []
    atom_id = 0
    for element, xyz, q in surface:
        atom_id += 1
        rows.append((atom_id, 0, TYPES[element], q, xyz))
    for molecule_index, (atoms, molecule_bonds) in enumerate(molecules, 1):
        first = atom_id + 1
        for element, xyz, q in atoms:
            atom_id += 1
            rows.append((atom_id, molecule_index, TYPES[element], q, xyz))
        for left, right in molecule_bonds:
            bonds.append((first + left, first + right))
    atoms_text = "\n".join(
        f"{i} {m} {t} {q} {x:.6f} {y:.6f} {z:.6f}" for i, m, t, q, (x, y, z) in rows
    )
    bonds_text = "\n".join(f"{n} 1 {a} {b}" for n, (a, b) in enumerate(bonds, 1))
    masses = "\n".join(f"{t} {m}" for t, m in MASSES.items())
    (build / "topology_output.lmp").write_text(
        f"""LAMMPS data

{len(rows)} atoms
{len(bonds)} bonds
6 atom types
1 bond types

0.0 {LX} xlo xhi
0.0 {LY} ylo yhi
-5.0 60.0 zlo zhi

Masses

{masses}

Atoms # full

{atoms_text}

Bonds

{bonds_text}
""",
        encoding="utf-8",
    )
    (build / "assembly_manifest.json").write_text(
        json.dumps(
            {
                "components": [
                    {"component": primary_name, "molecule_ids": [1, 2]},
                    {"component": "dcz-4p", "molecule_ids": [3, 3]},
                ]
            }
        ),
        encoding="utf-8",
    )
    frames = []
    for step in (0, 1000):
        lines = []
        for i, m, t, q, (x, y, z) in rows:
            lines.append(f"{i} {m} {t} {q} {x % LX:.6f} {y % LY:.6f} {z:.6f}")
        frames.append(
            f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(rows)}\n"
            f"ITEM: BOX BOUNDS pp pp ff\n0 {LX}\n0 {LY}\n-5 60\n"
            f"ITEM: ATOMS id mol type q x y z\n" + "\n".join(lines) + "\n"
        )
    (build / "relax-300K.lammpstrj").write_text("".join(frames), encoding="utf-8")
    return build


def write_reference(root: Path) -> Path:
    geometry = root / "me-template.xyz"
    geometry.write_text(
        f"{len(TEMPLATE)}\nlocal frame\n"
        + "\n".join(f"{e} {a} {b} {c}" for e, (a, b, c), _ in TEMPLATE)
        + "\n",
        encoding="utf-8",
    )
    path = root / "reference_dipoles.json"
    path.write_text(
        json.dumps(
            {
                "species": {
                    "me-4pacz": {
                        "dipole_e_angstrom": [0.0, 0.0, 1.0],
                        "geometry": geometry.name,
                        "method": "synthetic test",
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    return path


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_facets_match_the_real_corrugated_surface_counts():
    root = Path(__file__).parents[1]
    from nio_md_prep.analysis.interfacial import InterfacialTopology, _canonical_exposed_ni_sites
    from nio_md_prep.analysis.coverage import load_type_elements
    from nio_md_prep.lammps import parse

    path = root / "inputs/surfaces/corrugated-nio-110/surface.lmp"
    data = parse(path)
    types = load_type_elements(path)
    elements = {int(r.fields[0]): types[int(r.fields[2])] for r in data.sections["Atoms"]}
    molecules = {int(r.fields[0]): int(r.fields[1]) for r in data.sections["Atoms"]}
    coordinates = {int(r.fields[0]): tuple(map(float, r.fields[4:7])) for r in data.sections["Atoms"]}
    sites, source = _canonical_exposed_ni_sites(data, elements, molecules, 2.8)
    topology = InterfacialTopology(
        components=(),
        molecules=(),
        surface_site_ids=sites,
        atom_elements=elements,
        atom_molecules=molecules,
        atom_coordinates=coordinates,
        lateral_bounds=(data.bounds["x"], data.bounds["y"]),
        surface_coordination_cutoff=2.8,
        surface_site_source=source,
    )
    labels = classify_exposed_sites(topology)
    counts = {key: list(labels.values()).count(key) for key in set(labels.values())}
    assert counts == {"110+y": 210, "110-y": 210, "001": 180}
    families = build_facet_families(topology, labels)
    for family in families:
        assert family.row_spacing == pytest.approx(D * math.sqrt(2), abs=0.01)
        assert family.column_spacing == pytest.approx(2 * D, abs=0.01)
    windows = build_patch_windows(topology, families, column_stride=1)
    # 7 facet rows -> 3 row offsets; 30 Ni per row -> 30 column offsets; 2 facets.
    assert len(windows) == 2 * 3 * 30
    assert all(len(w.site_ids) == 20 for w in windows)
    assert windows[0].reference_area == pytest.approx(20 * 2 * D * D * math.sqrt(2), rel=1e-3)


def test_synthetic_facets_and_windows(tmp_path):
    build = write_build(tmp_path)
    topology = load_interfacial_topology(build)
    labels = classify_exposed_sites(topology)
    families = build_facet_families(topology, labels)
    assert [f.key for f in families] == ["110+y", "110-y"]
    assert all(len(f.site_ids) == 28 for f in families)
    windows = build_patch_windows(topology, families, column_stride=4)
    assert len(windows) == 6


def test_unwrapping_dipoles_and_reference_rotation(tmp_path):
    build = write_build(tmp_path)
    references = load_reference_dipoles(write_reference(tmp_path))
    warnings: list[str] = []
    prepared = _prepare(
        build,
        surface_coordination_cutoff=2.8,
        window_rows=5,
        window_columns=4,
        window_column_stride=4,
        aliases={"stage1_primary": "me-4pacz"},
        references=references,
        warnings=warnings,
    )
    assert warnings == ["no reference dipole for: dcz-4p; reference-dipole proxies cover only the listed species"]
    frame = next(iter_dump_frames(build / "relax-300K.lammpstrj"))
    state = _frame_state(np, prepared, frame, contact_cutoff=3.25, shield_cutoff=4.0)
    t1, t2, n = frame_vectors(N_MINUS)
    # Molecule B crosses the periodic x boundary; unwrapping must keep it rigid.
    np.testing.assert_allclose(state["dipole_charges"][0], 2.4 * n, atol=1e-6)
    np.testing.assert_allclose(state["dipole_charges"][1], 2.4 * t1, atol=1e-6)
    np.testing.assert_allclose(state["dipole_reference"][0], n, atol=1e-6)
    np.testing.assert_allclose(state["dipole_reference"][1], t1, atol=1e-6)
    assert np.all(np.isnan(state["dipole_reference"][2]))
    assert max(state["reference_rmsd"]) < 1e-6
    assert state["anchored_terminals"].tolist() == [1, 0, 2]


def test_analysis_outputs_patch_descriptors_and_proxies(tmp_path):
    build = write_build(tmp_path)
    summary_path = analyze_electronic_proxy(
        build,
        last_frames=None,
        window_column_stride=4,
        reference_dipoles=write_reference(tmp_path),
        component_aliases={"stage1_primary": "me-4pacz"},
        export_conformers=2,
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["components"] == ["me-4pacz", "dcz-4p"]
    assert summary["frames_analyzed"] == 2
    assert summary["model_assessment"]["profile"] == "classical-ff"
    assert summary["surface"]["windows"] == 6
    rows = read_csv(summary_path.parent / PATCH_CSV)
    assert len(rows) == 12
    area = 20 * 2 * D * D * math.sqrt(2)
    minus = [r for r in rows if r["facet"] == "110-y"]
    plus = [r for r in rows if r["facet"] == "110+y"]
    assert sorted(int(r["row_start"]) for r in minus) == [0, 0, 1, 1, 2, 2]
    for row in minus:
        # A sits in facet row 3 (all windows); B in row 5 (windows r1, r2).
        with_b = int(row["row_start"]) >= 1
        assert int(row["molecules_total"]) == (2 if with_b else 1)
        assert int(row["anchored_terminals_me-4pacz"]) == 1
        assert int(row["lifted_molecules"]) == (1 if with_b else 0)
        assert float(row["area_A2"]) == pytest.approx(area, rel=1e-6)
        assert float(row["dipole_normal_charges_e_A"]) == pytest.approx(2.4, abs=1e-6)
        expected = -_DIPOLE_DENSITY_TO_VOLTS * 2.4 / area
        assert float(row["dphi_layer_charges_eV"]) == pytest.approx(expected, rel=1e-6)
        assert float(row["dphi_layer_reference_eV"]) == pytest.approx(
            -_DIPOLE_DENSITY_TO_VOLTS * 1.0 / area, rel=1e-6
        )
        assert float(row["mean_tilt_deg_me-4pacz"]) == pytest.approx(
            45.0 if with_b else 0.0, abs=1e-6
        )
        # A's top C and all of B sit 8.0 A above their facet sites.
        assert float(row["max_molecule_height_A"]) == pytest.approx(8.0, abs=1e-6)
    for row in plus:
        # The V-groove must not leak A or B into the opposite facet's patches.
        assert int(row["molecules_total"]) == 1
        assert int(row["molecules_dcz-4p"]) == 1
        assert int(row["multi_anchored_molecules_dcz-4p"]) == 1
        assert int(row["lifted_molecules"]) == 0
        assert float(row["mean_tilt_deg_dcz-4p"]) == pytest.approx(0.0, abs=1e-6)
    frames = read_csv(summary_path.parent / "frame_proxies.csv")
    # Whole-interface z dipole: A contributes 2.4 * n_z, B lies along x.
    assert float(frames[0]["dipole_z_charges_e_A"]) == pytest.approx(2.4 * S, abs=1e-6)
    assert float(frames[0]["anchored_terminals_dcz-4p"]) == 2
    conformers = summary["outputs"]["conformers"]
    assert set(conformers) == {"me-4pacz", "dcz-4p"}
    text = (summary_path.parent / conformers["me-4pacz"]).read_text(encoding="utf-8")
    assert text.count("Properties=species:S:1:pos:R:3:charge:R:1") == 2


def test_cli_runs_analysis(tmp_path):
    build = write_build(tmp_path, primary_name="me-4pacz")
    output = tmp_path / "out"
    assert main(
        [
            "analyze-electronic-proxy",
            str(build),
            "--output",
            str(output),
            "--last-frames",
            "0",
            "--window-column-stride",
            "4",
        ]
    ) == 0
    assert (output / "electronic_proxy_summary.json").is_file()


def fake_run(root: Path, name: str, dphi: float, density: float, tilt_shift: float) -> Path:
    directory = root / name
    directory.mkdir(parents=True)
    rows = []
    for frame in range(3):
        for window in range(2):
            rows.append(
                {
                    "frame_index": frame,
                    "step": 1000 * frame,
                    "window_id": f"110+y:r{window}:c0",
                    "facet": "110+y",
                    "molecules_total": 4,
                    "anchored_terminal_density_total_per_nm2": density + 0.1 * window,
                    "mean_tilt_deg": 30.0 + tilt_shift + 5.0 * frame,
                    "exposed_site_fraction": 0.2 + 0.05 * window,
                    "lifted_molecules": 0,
                    "multi_anchored_molecules_dcz-4p": int(frame == 2),
                    "dphi_layer_charges_eV": dphi,
                    "molecule_ids": "1;2",
                }
            )
    with (directory / PATCH_CSV).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "build_directory": f"prepared/{name}",
        "trajectory": f"prepared/{name}/relax-300K.lammpstrj",
        "whole_interface": {
            "dphi_layer_charges_eV": {"mean": dphi},
            "anchored_terminal_density_total_per_nm2": {"mean": density},
        },
        "patches": {"percentiles": {"dphi_layer_charges_eV": {"p50": dphi}}},
        "outputs": {"patch_csv": PATCH_CSV},
    }
    (directory / "electronic_proxy_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    return directory


def test_gate_passes_only_when_difference_exceeds_between_run_spread(tmp_path):
    runs = [
        ("cosam", fake_run(tmp_path, "c1", -0.50, 1.0, 0)),
        ("cosam", fake_run(tmp_path, "c2", -0.52, 1.1, 0)),
        ("cosam", fake_run(tmp_path, "c3", -0.48, 0.9, 0)),
        ("sequential", fake_run(tmp_path, "s1", -0.80, 1.05, 10)),
        ("sequential", fake_run(tmp_path, "s2", -0.82, 0.95, 10)),
        ("sequential", fake_run(tmp_path, "s3", -0.78, 1.0, 10)),
    ]
    path = summarize_electronic_gate(runs, tmp_path / "gate", compare=("cosam", "sequential"))
    result = json.loads(path.read_text(encoding="utf-8"))
    metrics = {row["metric"]: row for row in result["metrics"]}
    assert metrics["whole:dphi_layer_charges_eV"]["passes"] is True
    assert metrics["whole:anchored_terminal_density_total_per_nm2"]["passes"] is False
    assert result["decision"]["proceed_to_tier1"] is True

    single = summarize_electronic_gate(runs[:4], tmp_path / "gate1", compare=("cosam", "sequential"))
    decision = json.loads(single.read_text(encoding="utf-8"))["decision"]
    assert decision["proceed_to_tier1"] is False
    assert "whole:dphi_layer_charges_eV" in decision["unevaluable_primary_metrics"]


def test_slot_selection_per_run_and_matched_pairs(tmp_path):
    runs = [
        ("cosam", fake_run(tmp_path, "c1", -0.5, 1.0, 0)),
        ("sequential", fake_run(tmp_path, "s1", -0.8, 1.0, 10)),
        ("gap-seeded-i", fake_run(tmp_path, "g1", -0.9, 1.2, 20)),
    ]
    path = select_patch_slots(
        runs, tmp_path / "slots", compare=("cosam", "sequential"), matched_pairs=1, reference_slots=1
    )
    slots = read_csv(path)
    criteria = [(s["label"], s["criterion"].split(" (")[0]) for s in slots]
    assert criteria[:5] == [
        ("cosam", "run-median anchored density"),
        ("cosam", "arrangement extreme"),
        ("sequential", "run-median anchored density"),
        ("sequential", "arrangement extreme"),
        ("gap-seeded-i", "run-median anchored density"),
    ]
    extremes = [s for s in slots if s["criterion"].startswith("arrangement extreme")]
    assert all(int(s["frame_index"]) == 2 for s in extremes)
    pair = [s for s in slots if s["criterion"].startswith("matched-density pair 1")]
    assert [s["label"] for s in pair] == ["cosam", "sequential"]
    assert not (tmp_path / "slots" / "patch_slots_warnings.txt").exists()


def test_alias_parsing_rejects_malformed_values():
    assert parse_aliases(["stage1_primary=me-4pacz"]) == {"stage1_primary": "me-4pacz"}
    with pytest.raises(ValueError):
        parse_aliases(["stage1_primary"])
