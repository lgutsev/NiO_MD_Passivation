"""LAMMPS deck rendering, launch resolution and output parsing -- no LAMMPS needed.

The rendered deck is what a LAMMPS run *is*, and ``mlip validate`` must be
able to show it from a laptop, so everything here is pure: the pair commands,
units (metal ps vs real fs, bar vs atm), the thermostat/barostat fixes, frozen
atoms, the launch command line and what it activates, the dump/print readers,
and the MACE export naming. ``tests/test_mlip_lammps_physics.py`` runs the
same decks through a real LAMMPS.
"""
from pathlib import Path

import pytest

from nio_md_prep.mlip.bridges import mace_lammps
from nio_md_prep.mlip.config import parse_job_mapping
from nio_md_prep.mlip.engines import lammps_engine
from nio_md_prep.mlip.errors import ConfigError, ResultError
from nio_md_prep.mlip.potentials.lammps_mlip import LammpsMlipAdapter
from nio_md_prep.mlip.registry import build_bridge
from nio_md_prep.mlip.specs import (
    EngineSpec,
    LammpsMlipPotentialSpec,
    MacePotentialSpec,
    SimulationSpec,
)
from nio_md_prep.mlip.units import ATM_IN_BAR

DEEPMD = LammpsMlipPotentialSpec(
    label="nio-deepmd",
    pair_style="deepmd models/nio.pb",
    pair_coeff=("* *",),
    type_map={1: "Ni", 2: "O"},
    units="metal",
    framework="deepmd",
    model_paths=(Path("models/nio.pb"),),
)
LJ_REAL = LammpsMlipPotentialSpec(
    pair_style="lj/cut 6.0", pair_coeff=("* * 1.0 2.5",), type_map={1: "Ni"}, units="real"
)

SINGLEPOINT = SimulationSpec(task="singlepoint", compute_per_atom_energy=True)
NVT = SimulationSpec(
    task="md",
    ensemble="nvt",
    temperature_K=400,
    timestep_fs=0.5,
    steps=100,
    seed=12345,
    thermostat_damping_fs=100.0,
    trajectory_interval=10,
)


def md(**kwargs):
    base = dict(task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=1.0, steps=20, seed=7)
    base.update(kwargs)
    return SimulationSpec(**base)


def render(potential=DEEPMD, simulation=SINGLEPOINT, **kwargs):
    return lammps_engine.render_deck(
        lammps_engine.rewrite_model_tokens(potential),
        simulation,
        n_types=len(potential.type_map),
        **kwargs,
    )


def line_starting(deck, prefix):
    return next(line for line in deck if line.startswith(prefix))


# --- structure and pair commands ------------------------------------------


def test_the_deck_sets_units_and_atom_style_before_reading_data():
    deck = render()
    assert "units metal" in deck and "atom_style atomic" in deck
    assert deck.index("units metal") < deck.index(f"read_data {lammps_engine.DATA_FILE}")


def test_model_files_are_named_by_their_staged_name_after_the_structure():
    """LAMMPS runs in the job directory, where the model is staged under its own name."""
    deck = render()
    assert "pair_style deepmd nio.pb" in deck
    assert "pair_coeff * *" in deck
    assert deck.index(f"read_data {lammps_engine.DATA_FILE}") < deck.index("pair_style deepmd nio.pb")


def test_boundary_is_per_axis_and_never_shrink_wrapped():
    assert "boundary p p f" in render(pbc=(True, True, False))
    assert "boundary f p p" in render(pbc=(False, True, True))
    assert "boundary f f f" in render(pbc=(False, False, False))
    assert "boundary p p p" in render()


def test_the_single_point_dumps_sorted_full_precision_forces_and_optional_computes():
    deck = render()
    assert "run 0" in deck
    dump = line_starting(deck, "dump mlip_forces")
    assert "id type xu yu zu fx fy fz c_mlip_pe_atom" in dump
    assert "dump_modify mlip_forces sort id format float %.17g" in deck
    assert "compute mlip_pe_atom all pe/atom" in deck
    # No stress requested -> no virial compute.
    assert not any("pressure NULL virial" in line for line in deck)
    plain = render(simulation=SimulationSpec(task="singlepoint", compute_stress=True))
    assert "compute mlip_virial all pressure NULL virial" in plain
    assert not any("pe/atom" in line for line in plain)


def test_thermo_values_are_printed_at_full_precision_for_the_executable_route():
    deck = render(simulation=SimulationSpec(task="singlepoint", compute_stress=True))
    assert 'print "step $(step:%.17g)" file thermo_result.txt screen no' in deck
    assert 'print "pe $(pe:%.17g)" append thermo_result.txt screen no' in deck
    assert 'print "virial_xy $(c_mlip_virial[4]:%.17g)" append thermo_result.txt screen no' in deck


def test_write_data_file_uses_the_type_map_not_an_existing_type_array(tmp_path, rattled_nio_structure):
    """An Atoms read from a LAMMPS data file carries a ``type`` array; it must not win."""
    atoms = rattled_nio_structure.copy()
    atoms.arrays["type"] = [9] * len(atoms)  # nonsense types that ASE would otherwise prefer
    path = lammps_engine.write_data_file(atoms, DEEPMD, tmp_path / "structure.lmp")
    text = path.read_text()
    rows = text.split("Atoms")[1].strip().splitlines()[1:]
    types = [int(row.split()[1]) for row in rows if row.strip()]
    expected = [1 if s == "Ni" else 2 for s in atoms.get_chemical_symbols()]
    assert types == expected


def test_a_type_map_with_one_element_twice_is_refused():
    with pytest.raises(ConfigError, match="more than one LAMMPS type to one element"):
        LammpsMlipPotentialSpec(pair_style="x", pair_coeff=("* *",), type_map={1: "Ni", 2: "Ni"})


def test_two_model_files_with_one_name_are_refused():
    spec = LammpsMlipPotentialSpec(
        pair_style="hybrid a/nio.pb b/nio.pb", pair_coeff=("* *",), type_map={1: "Ni"},
        model_paths=(Path("a/nio.pb"), Path("b/nio.pb")),
    )
    with pytest.raises(ConfigError, match="two files named"):
        lammps_engine.staged_model_names(spec)


# --- time, damping and pressure units --------------------------------------


def test_metal_timestep_and_damping_are_picoseconds():
    """A factor of 1000 here would silently ruin every trajectory."""
    deck = render(simulation=NVT)
    assert "timestep 0.0005" in deck
    assert line_starting(deck, "fix mlip_integrate").endswith("nvt temp 400 400 0.1")
    assert lammps_engine.TIMESTEP_PER_FS == {"metal": 1e-3, "real": 1.0}


def test_real_timestep_and_damping_are_femtoseconds():
    deck = render(potential=LJ_REAL, simulation=NVT)
    assert "timestep 0.5" in deck
    assert line_starting(deck, "fix mlip_integrate").endswith("nvt temp 400 400 100")


def test_barostat_targets_are_bar_in_metal_and_atm_in_real():
    """Never bar-as-atm: 1.01325 bar is exactly 1 atm."""
    sim = md(ensemble="npt", pressure_bar=ATM_IN_BAR, barostat_damping_fs=1000.0)
    metal = line_starting(render(simulation=sim), "fix mlip_integrate")
    real = line_starting(render(potential=LJ_REAL, simulation=sim), "fix mlip_integrate")
    assert metal.endswith("iso 1.01325 1.01325 1")
    assert real.endswith("iso 1 1 1000")
    plan = lammps_engine.plan_md(sim, units="real")["record"]
    assert plan["pressure_lammps"] == pytest.approx(1.0) and plan["pressure_lammps_unit"] == "atm"


def test_the_units_plan_names_the_pressure_unit_and_what_was_verified():
    real = lammps_engine.units_plan("real")
    assert real["pressure_unit"] == "atm" and real["pressure_unit_in_bar"] == ATM_IN_BAR
    assert real["time_unit"] == "fs" and lammps_engine.units_plan("metal")["time_unit"] == "ps"
    assert {"stress", "npt"} <= set(real["verified"])


def test_an_unverified_unit_quantity_is_refused(monkeypatch):
    monkeypatch.setitem(
        lammps_engine.VERIFIED_UNIT_STYLES, "real", frozenset({"energy", "forces"})
    )
    with pytest.raises(ConfigError, match="stress has not been verified"):
        lammps_engine.check_units("real", SimulationSpec(task="singlepoint", compute_stress=True))
    with pytest.raises(ConfigError, match="npt"):
        lammps_engine.check_units("real", md(ensemble="npt", pressure_bar=1.0))
    lammps_engine.check_units("real", SimulationSpec(task="singlepoint"))


# --- thermostats and barostats: implemented or refused -------------------


@pytest.mark.parametrize(
    ("thermostat", "fixes"),
    [
        (None, ["fix mlip_integrate all nvt temp 300 300 0.1"]),
        ("nose-hoover", ["fix mlip_integrate all nvt temp 300 300 0.1"]),
        ("langevin", ["fix mlip_integrate all nve",
                      "fix mlip_thermostat all langevin 300 300 0.1 7 zero yes tally yes"]),
        ("berendsen", ["fix mlip_integrate all nve",
                       "fix mlip_thermostat all temp/berendsen 300 300 0.1"]),
        ("csvr", ["fix mlip_integrate all nve", "fix mlip_thermostat all temp/csvr 300 300 0.1 7"]),
    ],
)
def test_every_thermostat_renders_the_fix_that_implements_it(thermostat, fixes):
    deck = render(simulation=md(thermostat=thermostat), seed=7)
    rendered = [line for line in deck if line.startswith(("fix mlip_integrate", "fix mlip_thermostat"))]
    assert rendered == fixes
    record = lammps_engine.plan_md(md(thermostat=thermostat), units="metal", seed=7)["record"]
    assert record["thermostat"] == (thermostat or "nose-hoover")
    assert record["thermostat_default_used"] is (thermostat is None)
    assert record["thermostat_damping_fs"] == 100.0 and record["thermostat_damping_lammps"] == 0.1


def test_nve_integrates_without_a_thermostat_and_conserves_etotal():
    deck = render(simulation=md(ensemble="nve"))
    assert "fix mlip_integrate all nve" in deck
    assert not any("mlip_thermostat" in line for line in deck)
    assert lammps_engine.plan_md(md(ensemble="nve"), units="metal")["record"]["conserved_quantity"] == "etotal"


def test_every_fix_is_removed_before_the_final_single_point():
    """Langevin friction must not be in the forces the run returns."""
    deck = render(simulation=md(thermostat="langevin"))
    final = deck.index("run 0 post no")
    assert deck.index("unfix mlip_thermostat") < final and deck.index("unfix mlip_integrate") < final


@pytest.mark.parametrize(
    ("coupling", "pbc", "expected"),
    [
        (None, (True, True, True), "all npt temp 300 300 0.1 iso 1 1 1"),
        ("anisotropic", (True, True, True), "all npt temp 300 300 0.1 aniso 1 1 1"),
        ("in-plane", (True, True, False), "all npt temp 300 300 0.1 x 1 1 1 y 1 1 1 couple xy"),
        ("in-plane", (False, True, True), "all npt temp 300 300 0.1 y 1 1 1 z 1 1 1 couple yz"),
        ("in-plane", (True, False, True), "all npt temp 300 300 0.1 x 1 1 1 z 1 1 1 couple xz"),
    ],
)
def test_barostat_coupling_follows_the_geometry(coupling, pbc, expected):
    sim = md(ensemble="npt", pressure_bar=1.0, barostat_coupling=coupling)
    assert line_starting(render(simulation=sim, pbc=pbc), "fix mlip_integrate").endswith(expected)


def test_isotropic_npt_is_refused_for_a_slab_or_a_vacuum_gap():
    sim = md(ensemble="npt", pressure_bar=1.0)
    with pytest.raises(ConfigError, match="in-plane"):
        lammps_engine.plan_md(sim, units="metal", pbc=(True, True, False))
    with pytest.raises(ConfigError, match="vacuum"):
        lammps_engine.plan_md(sim, units="metal", vacuum=(2,))
    # A fully periodic cell with one vacuum gap barostats the other two axes.
    inplane = md(ensemble="npt", pressure_bar=1.0, barostat_coupling="in-plane")
    record = lammps_engine.plan_md(inplane, units="metal", vacuum=(2,))["record"]
    assert record["barostatted_lammps_dimensions"] == ["x", "y"]


def test_parrinello_rahman_is_refused_not_mapped():
    with pytest.raises(ConfigError, match="parrinello-rahman"):
        lammps_engine.plan_md(md(ensemble="npt", pressure_bar=1.0, barostat="parrinello-rahman"), units="metal")


def test_berendsen_barostat_needs_a_modulus_and_an_orthogonal_box():
    sim = md(ensemble="npt", pressure_bar=1.0, barostat="berendsen", thermostat="csvr")
    with pytest.raises(ConfigError, match="bulk modulus"):
        lammps_engine.plan_md(sim, units="metal")
    options = {"barostat_bulk_modulus_bar": 1.9e6}
    with pytest.raises(ConfigError, match="triclinic"):
        lammps_engine.plan_md(sim, units="metal", options=options, triclinic=True)
    plan = lammps_engine.plan_md(sim, units="real", options=options)
    assert plan["commands"][-1].endswith(f"modulus {1.9e6 / ATM_IN_BAR!r}")
    assert plan["record"]["fixes"] == ["nve", "temp/csvr", "press/berendsen"]
    assert plan["record"]["conserved_quantity"] is None


def test_frozen_atoms_are_excluded_from_velocity_integration_and_temperature():
    deck = render(simulation=md(thermostat="langevin"), fixed_ids=[1, 2, 3, 7], n_atoms=16, seed=7)
    assert "group mlip_frozen id 1:3 7" in deck
    assert "group mlip_mobile subtract all mlip_frozen" in deck
    assert "compute mlip_temp mlip_mobile temp" in deck
    assert "compute_modify mlip_temp extra/dof 0" in deck
    assert line_starting(deck, "velocity").startswith("velocity mlip_mobile create 300 7 mom no")
    assert "fix mlip_integrate mlip_mobile nve" in deck
    assert line_starting(deck, "fix mlip_thermostat").startswith("fix mlip_thermostat mlip_mobile langevin")


def test_npt_with_frozen_atoms_or_everything_frozen_is_refused():
    with pytest.raises(ConfigError, match="frozen"):
        lammps_engine.plan_md(md(ensemble="npt", pressure_bar=1.0), units="metal", n_fixed=2)
    with pytest.raises(ConfigError, match="every atom is frozen"):
        render(simulation=md(), fixed_ids=[1, 2], n_atoms=2)


def test_the_trajectory_always_ends_at_the_final_step():
    deck = render(simulation=md(steps=25, trajectory_interval=10, log_interval=10))
    assert any(line.startswith("write_dump all custom smoke_md.lammpstrj") for line in deck)
    assert lammps_engine.expected_steps(25, 10) == [0, 10, 20, 25]


# --- launch: route, command line, accelerators ----------------------------


def test_an_explicit_executable_is_never_replaced_by_the_python_module():
    launch = lammps_engine.resolve_launch(EngineSpec(kind="lammps", executable="lmp_mpi"))
    assert launch.route == "executable"
    assert launch.argv() == ["lmp_mpi", "-in", "in.lammps", "-log", "log.lammps"]
    mpi = lammps_engine.resolve_launch(
        EngineSpec(kind="lammps", runtime="executable", mpi_launcher=("mpirun", "-np", "4"))
    )
    assert mpi.argv()[:4] == ["mpirun", "-np", "4", "lmp"]


def test_threads_become_openmp_and_are_recorded_in_the_command_and_env():
    launch = lammps_engine.resolve_launch(EngineSpec(kind="lammps", runtime="executable", threads=4))
    assert launch.cmdargs == ("-pk", "omp", "4", "-sf", "omp")
    assert launch.env == {"OMP_NUM_THREADS": "4"}
    assert launch.command_line().startswith("OMP_NUM_THREADS=4 lmp -pk omp 4 -sf omp")
    with pytest.raises(ConfigError, match="already choose an accelerator"):
        lammps_engine.resolve_launch(EngineSpec(kind="lammps", threads=4, lammps_args=("-sf", "gpu")))


@pytest.mark.parametrize(
    ("args", "gpu", "kokkos"),
    [
        (("-k", "on", "g", "1", "-sf", "kk", "-pk", "kokkos", "newton", "on", "neigh", "half"), True, True),
        (("-k", "on", "t", "4", "-sf", "kk"), False, True),
        (("-k", "on", "g", "1"), False, True),  # no kk suffix: pair styles stay on the host
        (("-sf", "gpu"), True, False),
        (("-pk", "omp", "4", "-sf", "omp"), False, False),
        ((), False, False),
    ],
)
def test_gpu_is_claimed_only_when_the_switches_put_the_pair_style_on_a_gpu(args, gpu, kokkos):
    accelerator = lammps_engine.parse_accelerator_args(args)
    assert accelerator["gpu"] is gpu and accelerator["kokkos"] is kokkos
    engine = EngineSpec(kind="lammps", lammps_args=args)
    assert lammps_engine.LammpsEngine(engine).capabilities().gpu is gpu


def test_lmp_help_is_parsed_for_packages_pair_styles_and_kokkos_backends():
    text = (
        "Large-scale Atomic/Molecular Massively Parallel Simulator - 22 Jul 2025\n\n"
        "Accelerator configuration:\n\nKOKKOS package API: CUDA Serial\n"
        "KOKKOS package precision: double\n\n"
        "Installed packages:\n\nKOKKOS ML-IAP ML-SNAP PYTHON\n\n"
        "List of individual style options included in this LAMMPS executable\n\n"
        "* Pair styles:\n\nlj/cut          mliap           mliap/kk\n\n"
        "* Fix styles:\n\nnve\n"
    )
    build = lammps_engine.parse_help(text)
    assert build["packages"] == ["KOKKOS", "ML-IAP", "ML-SNAP", "PYTHON"]
    assert build["pair_styles"] == ["lj/cut", "mliap", "mliap/kk"]
    assert build["accelerators"]["KOKKOS"]["api"] == ["cuda", "serial"]
    assert lammps_engine.missing_features(build, [("kokkos_backend", "gpu")]) == []
    assert lammps_engine.missing_features(build, [("kokkos_backend", "host-only")]) == [
        ("kokkos_backend", "host-only")
    ]
    assert lammps_engine.missing_features({"packages": ["ML-IAP"]}, [("package", "PYTHON")]) == [
        ("package", "PYTHON")
    ]


# --- reading LAMMPS output -------------------------------------------------


def dump_frame(timestep, rows, *, n=None):
    header = [
        "ITEM: TIMESTEP", str(timestep), "ITEM: NUMBER OF ATOMS", str(n if n is not None else len(rows)),
        "ITEM: BOX BOUNDS pp pp pp", "0 10", "0 10", "0 10", "ITEM: ATOMS id type xu yu zu",
    ]
    return header + [f"{i} 1 {x} 0 0" for i, x in rows]


def test_a_multi_frame_dump_is_read_frame_by_frame(tmp_path):
    path = tmp_path / "t.lammpstrj"
    path.write_text("\n".join(dump_frame(0, [(1, 0.5), (2, 1.5)]) + dump_frame(10, [(1, 0.6), (2, 1.6)])))
    frames = lammps_engine.read_dump(path, n_atoms=2, expected_timesteps=[0, 10])
    assert [f.timestep for f in frames] == [0, 10]
    assert frames[1].columns3("xu", "yu", "zu")[:, 0].tolist() == [0.6, 1.6]
    assert frames[0].box["lammps_cell"] == [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]


def test_a_truncated_short_or_unsorted_dump_is_an_error_not_a_shorter_result(tmp_path):
    path = tmp_path / "t.lammpstrj"
    good = dump_frame(0, [(1, 0.5), (2, 1.5)])
    path.write_text("\n".join(good + dump_frame(10, [(1, 0.6)], n=2)))
    with pytest.raises(ResultError, match="truncated"):
        lammps_engine.read_dump(path, n_atoms=2)
    path.write_text("\n".join(good))
    with pytest.raises(ResultError, match="frames at timesteps"):
        lammps_engine.read_dump(path, n_atoms=2, expected_timesteps=[0, 10])
    path.write_text("\n".join(dump_frame(0, [(2, 0.5), (1, 1.5)])))
    with pytest.raises(ResultError, match="sort id"):
        lammps_engine.read_dump(path, n_atoms=2)
    with pytest.raises(ResultError, match="holds 2 atoms"):
        lammps_engine.read_dump(path, n_atoms=3)
    path.write_text("\n".join(good[:6]))
    with pytest.raises(ResultError, match="truncated"):
        lammps_engine.read_dump(path, n_atoms=2)


def test_the_thermo_series_must_cover_every_logged_step(tmp_path):
    path = tmp_path / "s.dat"
    path.write_text("# step pe ke etotal ecouple econserve temp\n0 1 2 3 0 3 300\n5 1 2 3 0 3 301\n")
    with pytest.raises(ResultError, match="expected"):
        lammps_engine.read_series(path, steps=10, interval=5)
    path.write_text(path.read_text() + "10 1 2 3 0 3 302\n10 1 2 3 0 3 302\n")
    assert len(lammps_engine.read_series(path, steps=10, interval=5)) == 3


# --- framework hints --------------------------------------------------------


def test_per_atom_energy_is_not_claimed_for_mace_through_lammps():
    unified = LammpsMlipPotentialSpec(pair_style="mliap unified m.pt 0", pair_coeff=("* * Ni",), type_map={1: "Ni"})
    snap = LammpsMlipPotentialSpec(pair_style="mliap model linear a descriptor sna b", pair_coeff=("* * Ni",), type_map={1: "Ni"})
    pair_mace = LammpsMlipPotentialSpec(pair_style="mace no_domain_decomposition", pair_coeff=("* * m.pt Ni",), type_map={1: "Ni"})
    assert not LammpsMlipAdapter(unified).capabilities().per_atom_energy
    assert not LammpsMlipAdapter(pair_mace).capabilities().per_atom_energy
    assert LammpsMlipAdapter(snap).capabilities().per_atom_energy


def test_engine_precision_is_refused_for_a_lammps_native_pair_style():
    bridge = build_bridge(DEEPMD, EngineSpec(kind="lammps", precision="float64"))
    with pytest.raises(ConfigError, match="engine.precision"):
        bridge.check_simulation(SimulationSpec(task="singlepoint"))


# --- MACE projected onto LAMMPS ---------------------------------------------


def mace_spec(tmp_path, **kwargs):
    return MacePotentialSpec(
        label="nio",
        model_path=tmp_path / "nio.model",
        declared_elements=("Ni", "O", "P", "C", "H"),
        **kwargs,
    )


def test_the_export_path_is_what_mace_create_lammps_model_writes(tmp_path):
    """create_lammps_model.py:107 / :110 append to the checkpoint's full path."""
    mliap = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps"))
    pair = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps"), implementation="pair-mace")
    assert mliap.exported_model_path() == tmp_path / "nio.model-mliap_lammps.pt"
    assert pair.exported_model_path() == tmp_path / "nio.model-lammps.pt"
    assert mace_lammps.EXPORT_SUFFIX == {"mliap": "-mliap_lammps.pt", "pair-mace": "-lammps.pt"}


def test_the_export_path_can_be_overridden_and_is_reported(tmp_path):
    export = tmp_path / "custom export.pt"
    export.write_bytes(b"not a real model")
    engine = EngineSpec(kind="lammps", options={"exported_model_path": str(export)})
    bridge = build_bridge(mace_spec(tmp_path), engine)
    record = bridge.export_record()
    assert bridge.exported_model_path() == export
    assert record["exported_model"]["exists"] is True
    assert record["exported_model"]["source"] == "engine.options.exported_model_path"
    assert len(record["exported_model"]["sha256"]) == 64
    assert record["checkpoint"]["exists"] is False and record["checkpoint"]["sha256"] is None
    assert bridge.pair_style() == "mliap unified custom_export.pt 0"
    legacy = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps", options={"model_path": str(export)}))
    assert legacy.exported_model_path() == export
    both = EngineSpec(kind="lammps", options={"model_path": "a.pt", "exported_model_path": "b.pt"})
    with pytest.raises(ConfigError, match="different exported models"):
        build_bridge(mace_spec(tmp_path), both).exported_model_path()


def test_a_config_file_resolves_the_exported_model_path_against_its_directory(tmp_path):
    job = parse_job_mapping(
        {
            "potential": {"kind": "mace", "model_path": "nio.model", "elements": ["Ni", "O"]},
            "engine": {"kind": "lammps", "options": {"exported_model_path": "exports/nio.pt"}},
            "simulation": {"task": "singlepoint"},
        },
        base_dir=tmp_path,
    )
    assert Path(job.engine.options["exported_model_path"]) == (tmp_path / "exports" / "nio.pt").resolve()


def test_a_missing_export_is_unavailable_with_the_exporter_command(tmp_path):
    bridge = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps"))
    availability = bridge.availability()
    assert not availability
    assert "nio.model-mliap_lammps.pt" in availability.detail
    assert "MACE_ALLOW_CPU=true mace_create_lammps_model" in availability.detail
    assert "--format mliap --dtype float64" in availability.detail


def test_mliap_cpu_uses_host_kokkos_and_records_the_launch(tmp_path):
    """Multi-layer MACE needs KOKKOS forward_exchange, so the CPU default is host KOKKOS."""
    bridge = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps", runtime="executable"))
    launch = bridge.launch()
    assert launch.cmdargs == ("-k", "on", "t", "1", "-sf", "kk", "-pk", "kokkos", "newton", "on", "neigh", "half")
    assert launch.activate == "mliappy_kokkos"
    assert ("kokkos_backend", "host-only") in launch.required and ("package", "PYTHON") in launch.required
    # MACE_ALLOW_CPU is pickled into the export when it is created; setting it
    # for the run would have no effect, so it is not set (or recorded).
    assert launch.env == {}
    assert not bridge.capabilities().gpu
    assert not bridge.capabilities().per_atom_energy


def test_mliap_cuda_activates_kokkos_on_the_gpu_and_only_then_claims_it(tmp_path):
    bridge = build_bridge(mace_spec(tmp_path, device="cuda"), EngineSpec(kind="lammps", runtime="executable"))
    launch = bridge.launch()
    assert launch.cmdargs == lammps_engine.KOKKOS_GPU_ARGS
    assert ("kokkos_backend", "gpu") in launch.required
    assert bridge.capabilities().gpu
    assert "-k on g 1 -sf kk -pk kokkos newton on neigh half" in launch.command_line()
    python_route = build_bridge(mace_spec(tmp_path, device="cuda"), EngineSpec(kind="lammps", runtime="python"))
    assert "activate_mliappy_kokkos(lmp)" in python_route.launch().command_line()


def test_mliap_launch_arguments_must_agree_with_the_device(tmp_path):
    host = EngineSpec(kind="lammps", lammps_args=("-k", "on", "t", "2", "-sf", "kk"))
    with pytest.raises(ConfigError, match="must agree"):
        build_bridge(mace_spec(tmp_path, device="cuda"), host).launch()
    no_kokkos = EngineSpec(kind="lammps", lammps_args=("-sf", "omp"))
    with pytest.raises(ConfigError, match="couples through KOKKOS"):
        build_bridge(mace_spec(tmp_path), no_kokkos).launch()
    with pytest.raises(ConfigError, match="CUDA_VISIBLE_DEVICES"):
        build_bridge(mace_spec(tmp_path, device="cuda:1"), EngineSpec(kind="lammps")).launch()
    with pytest.raises(ConfigError, match="mps"):
        build_bridge(mace_spec(tmp_path, device="mps"), EngineSpec(kind="lammps")).launch()


def test_the_plain_mliap_coupling_is_cpu_only(tmp_path):
    plain = EngineSpec(kind="lammps", options={"mliap_coupling": "plain"})
    launch = build_bridge(mace_spec(tmp_path), plain).launch()
    assert launch.cmdargs == () and launch.activate == "mliappy"
    assert "forward_exchange" in launch.notes[0]
    with pytest.raises(ConfigError, match="plain"):
        build_bridge(mace_spec(tmp_path, device="cuda"), plain).launch()


def test_mace_pair_commands_render_newton_on_and_staged_names(tmp_path):
    bridge = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps"))
    spec = bridge.as_lammps_spec()
    deck = lammps_engine.render_deck(spec, SINGLEPOINT, n_types=5)
    assert "newton on" in deck and deck.index("newton on") < deck.index("units metal")
    assert "pair_style mliap unified nio.model-mliap_lammps.pt 0" in deck
    assert "pair_coeff * * Ni O P C H" in deck
    assert spec.type_map == {1: "Ni", 2: "O", 3: "P", 4: "C", 5: "H"}


def test_pair_mace_is_legacy_and_single_rank_only_skips_domain_decomposition(tmp_path):
    single = build_bridge(mace_spec(tmp_path), EngineSpec(kind="lammps"), implementation="pair-mace")
    assert single.pair_style() == "mace no_domain_decomposition"
    assert single.pair_coeff() == ("* * nio.model-lammps.pt Ni O P C H",)
    assert mace_lammps.PAIR_MACE_LEGACY_WARNING in single.warnings()
    assert single.launch().env == {"CUDA_VISIBLE_DEVICES": ""}  # pair_mace picks CUDA whenever it can
    mpi = build_bridge(
        mace_spec(tmp_path),
        EngineSpec(kind="lammps", mpi_launcher=("mpirun", "-np", "2")),
        implementation="pair-mace",
    )
    assert mpi.pair_style() == "mace"


def test_compile_mode_is_refused_for_a_lammps_export(tmp_path):
    export = tmp_path / "nio.model-mliap_lammps.pt"
    export.write_bytes(b"x")
    bridge = build_bridge(mace_spec(tmp_path, compile_mode="default"), EngineSpec(kind="lammps"))
    with pytest.raises(ConfigError, match="compile_mode"):
        bridge.check_simulation(SimulationSpec(task="singlepoint"))


# --- execution plans ------------------------------------------------------

PLAN_KEYS = {
    "potential_kind", "engine", "bridge", "implementation", "model_checkpoint", "exported_model",
    "elements", "energy_convention", "units", "device", "precision", "dynamics", "lammps",
    "openmm", "availability", "unmet_capabilities", "model_hashes",
}
DYNAMICS_KEYS = {
    "ensemble", "integrator", "thermostat", "barostat", "barostat_coupling", "timestep_fs",
    "timestep_native", "thermostat_damping_fs", "thermostat_damping_native", "barostat_damping_fs",
    "barostat_damping_native",
}


def check_plan_schema(plan):
    assert PLAN_KEYS <= set(plan)
    assert {"path", "sha256"} <= set(plan["model_checkpoint"])
    if plan["exported_model"] is not None:
        assert {"path", "exists", "sha256", "implementation"} <= set(plan["exported_model"])
    assert {"native", "pressure_unit"} <= set(plan["units"])
    for key in ("device", "precision"):
        assert {"requested", "effective", "guaranteed", "note"} <= set(plan[key])
    if plan["dynamics"] is not None:
        assert DYNAMICS_KEYS <= set(plan["dynamics"])
    assert {"pair_style", "pair_coeff", "launch_command", "kokkos_args", "env"} <= set(plan["lammps"])
    assert plan["openmm"] is None
    assert {"available", "missing"} <= set(plan["availability"])
    assert isinstance(plan["unmet_capabilities"], list)
    assert {"declared", "observed"} <= set(plan["model_hashes"])


@pytest.mark.parametrize("engine_kind", ["lammps", "ase"])
def test_lammps_native_execution_plans_have_the_shared_schema(engine_kind, rattled_nio_structure):
    lj = LammpsMlipPotentialSpec(
        pair_style="lj/cut 6.0", pair_coeff=("* * 0.01 2.5",), type_map={1: "Ni", 2: "O"}
    )
    bridge = build_bridge(lj, EngineSpec(kind=engine_kind))
    for simulation in (SimulationSpec(task="singlepoint", compute_stress=True), md(thermostat="langevin")):
        plan = bridge.execution_plan(simulation, rattled_nio_structure)
        check_plan_schema(plan)
    assert plan["dynamics"]["thermostat"] in ("langevin", None)
    if engine_kind == "lammps":
        assert plan["dynamics"]["thermostat_damping_native"] == 0.1
        assert plan["dynamics"]["timestep_native"] == 0.001
        assert plan["lammps"]["boundary"] == "p p p"
        assert plan["device"] == {
            "requested": None, "effective": "cpu", "guaranteed": True,
            "note": "no accelerator switches: the pair style runs on the host CPU",
        }


def test_a_refused_dynamics_request_keeps_the_schema(rattled_nio_structure):
    lj = LammpsMlipPotentialSpec(pair_style="lj/cut 6.0", pair_coeff=("* * 0.01 2.5",), type_map={1: "Ni", 2: "O"})
    plan = build_bridge(lj, EngineSpec(kind="lammps")).execution_plan(
        md(ensemble="npt", pressure_bar=1.0, barostat="parrinello-rahman"), rattled_nio_structure
    )
    check_plan_schema(plan)
    assert "parrinello-rahman" in plan["dynamics"]["refused"]


@pytest.mark.parametrize("implementation", ["mliap", "pair-mace"])
@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_mace_lammps_execution_plans_report_checkpoint_export_and_launch(tmp_path, implementation, device):
    from nio_md_prep.mlip.provenance import sha256_file

    checkpoint = tmp_path / "nio.model"
    checkpoint.write_bytes(b"checkpoint")
    export = tmp_path / ("nio.model" + mace_lammps.EXPORT_SUFFIX[implementation])
    export.write_bytes(b"export")
    bridge = build_bridge(
        mace_spec(tmp_path, device=device),
        EngineSpec(kind="lammps", runtime="executable"),
        implementation=implementation,
    )
    plan = bridge.execution_plan(md(ensemble="nve"), None)
    check_plan_schema(plan)
    assert plan["model_checkpoint"]["sha256"] == sha256_file(checkpoint)
    assert plan["exported_model"]["sha256"] == sha256_file(export)
    assert plan["exported_model"]["path"] == str(export)
    assert plan["exported_model"]["implementation"] == implementation
    assert plan["model_hashes"]["observed"] == {
        str(checkpoint): sha256_file(checkpoint), str(export): sha256_file(export)
    }
    assert plan["lammps"]["launch_command"].startswith(("lmp", "CUDA_VISIBLE_DEVICES"))
    kokkos = plan["lammps"]["kokkos_args"]
    if implementation == "mliap" and device == "cuda":
        assert kokkos == list(lammps_engine.KOKKOS_GPU_ARGS)
        assert plan["device"]["effective"].startswith("cuda")
    elif implementation == "mliap":
        assert kokkos[:4] == ["-k", "on", "t", "1"]
    else:
        assert kokkos == []
        assert mace_lammps.PAIR_MACE_LEGACY_WARNING in plan["warnings"]
    # Nothing can promise the device before LAMMPS starts, except a CPU run of
    # an executable with every GPU hidden from it.
    assert plan["device"]["guaranteed"] is (implementation == "pair-mace" and device == "cpu")
    # torch cannot read these fake exports here, so the precision is unverified.
    assert plan["precision"]["guaranteed"] is False and plan["precision"]["effective"] is None
    assert plan["dynamics"]["timestep_native"] == 0.001
