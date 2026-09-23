"""LAMMPS routes against physics that can be checked independently.

Everything here runs a real LAMMPS (the python module, an ``lmp`` executable,
or ASE's LAMMPSlib) and compares with something that does not share its code:

* a Lennard-Jones mixture evaluated by a numpy reference (energy, forces,
  virial stress, per-atom energies) in every geometry class the engine
  handles -- cubic, primitive (triclinic), rotated, slabs open along x, y
  and z, a tilted slab and a cluster;
* the same system in LAMMPS ``real`` units, which must reproduce the
  ``metal`` results and trajectories (so kcal/mol, fs and atm are all
  converted, and a bar-for-atm slip would show);
* thermostats that must hold a temperature NVE cannot, with the thermostat
  energy tallied into the conserved quantity; barostats that move only the
  dimensions they should; frozen atoms that do not move.

The MACE export tests at the end run mace-torch's own exporter and are marked
``mace``; they run only where mace-torch is installed.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from importlib.util import find_spec
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("ase")

from ase.build import bulk  # noqa: E402
from ase.constraints import FixAtoms  # noqa: E402
from ase.geometry import cellpar_to_cell  # noqa: E402
from ase.neighborlist import neighbor_list  # noqa: E402

from nio_md_prep.mlip.engines import lammps_engine  # noqa: E402
from nio_md_prep.mlip.errors import MissingDependencyError, MlipError, ResultError  # noqa: E402
from nio_md_prep.mlip.potentials import lammps_mlip  # noqa: E402
from nio_md_prep.mlip.registry import build_bridge  # noqa: E402
from nio_md_prep.mlip.specs import (  # noqa: E402
    EngineSpec,
    LammpsMlipPotentialSpec,
    MacePotentialSpec,
    SimulationSpec,
)
from nio_md_prep.mlip.units import ATM_IN_BAR, KCAL_PER_MOL_IN_EV  # noqa: E402

pytestmark = pytest.mark.lammps

HAVE_MODULE = find_spec("lammps") is not None
HAVE_EXECUTABLE = shutil.which("lmp") is not None

CUTOFF = 6.0
#: (epsilon eV, sigma Angstrom) of an NiO-like LJ mixture; Ni-O at 2.085 A
#: sits near its minimum (2^(1/6) * 1.9 = 2.13 A), so the rocksalt cell is stable.
LJ = {("Ni", "Ni"): (0.10, 2.6), ("O", "O"): (0.05, 2.6), ("Ni", "O"): (0.07, 1.9)}


def lj_potential(units: str = "metal", *, scale_to_units: bool = True) -> LammpsMlipPotentialSpec:
    """The LJ mixture as a LAMMPS pair style; in ``real`` units epsilon is kcal/mol."""
    factor = 1.0 / KCAL_PER_MOL_IN_EV if units == "real" and scale_to_units else 1.0
    coeffs = []
    for (a, b), (eps, sig) in LJ.items():
        i, j = (1 if a == "Ni" else 2), (1 if b == "Ni" else 2)
        coeffs.append(f"{i} {j} {eps * factor!r} {sig!r}")
    return LammpsMlipPotentialSpec(
        label=f"lj-{units}",
        pair_style=f"lj/cut {CUTOFF}",
        pair_coeff=tuple(coeffs),
        type_map={1: "Ni", 2: "O"},
        units=units,
    )


def reference_lj(atoms):
    """Independent LJ: energy, forces, per-atom energies and virial stress (eV, A)."""
    i, j, D = neighbor_list("ijD", atoms, CUTOFF)
    symbols = atoms.get_chemical_symbols()
    energy = 0.0
    per_atom = np.zeros(len(atoms))
    forces = np.zeros((len(atoms), 3))
    virial = np.zeros((3, 3))
    for a, b, d in zip(i, j, D):
        key = (symbols[a], symbols[b]) if (symbols[a], symbols[b]) in LJ else (symbols[b], symbols[a])
        eps, sig = LJ[key]
        r = float(np.linalg.norm(d))
        sr6 = (sig / r) ** 6
        pair = 4.0 * eps * (sr6 * sr6 - sr6)
        energy += 0.5 * pair
        per_atom[a] += 0.5 * pair
        dphi = 4.0 * eps * (-12.0 * sr6 * sr6 + 6.0 * sr6) / r
        forces[a] += dphi * d / r
        virial += 0.5 * dphi * np.outer(d, d) / r
    stress = virial / atoms.get_volume()
    voigt = np.array([stress[0, 0], stress[1, 1], stress[2, 2], stress[1, 2], stress[0, 2], stress[0, 1]])
    return energy, forces, per_atom, voigt


# --- geometries -------------------------------------------------------------


def nio_cubic():
    atoms = bulk("NiO", "rocksalt", a=4.17, cubic=True).repeat((2, 2, 2))
    atoms.rattle(0.05, seed=3)
    return atoms


def nio_primitive():
    atoms = bulk("NiO", "rocksalt", a=4.17).repeat((2, 2, 2))
    atoms.rattle(0.05, seed=3)
    return atoms


def nio_rotated():
    atoms = nio_primitive()
    atoms.rotate(33, "x", rotate_cell=True)
    atoms.rotate(21, "z", rotate_cell=True)
    return atoms


def nio_slab(axis):
    atoms = nio_cubic()
    atoms.center(vacuum=8.0, axis=axis)
    pbc = [True, True, True]
    pbc[axis] = False
    atoms.pbc = pbc
    return atoms


def nio_tilted_slab():
    atoms = bulk("NiO", "rocksalt", a=4.17, cubic=True).repeat((2, 2, 2))
    fractional = atoms.get_scaled_positions() * [1, 1, 8.34 / 25.0] + [0, 0, 0.3]
    atoms.set_cell(cellpar_to_cell([8.34, 8.34, 25.0, 80.0, 75.0, 90.0]))
    atoms.set_scaled_positions(fractional)
    atoms.rattle(0.05, seed=5)
    atoms.pbc = (True, True, False)
    return atoms


def nio_cluster():
    atoms = bulk("NiO", "rocksalt", a=4.17, cubic=True)
    atoms.rattle(0.05, seed=3)
    atoms.center(vacuum=8.0)
    atoms.pbc = False
    return atoms


GEOMETRIES = {
    "cubic": (nio_cubic, "p p p", False),
    "primitive": (nio_primitive, "p p p", True),
    "rotated": (nio_rotated, "p p p", True),
    "slab_x": (lambda: nio_slab(0), "f p p", False),
    "slab_y": (lambda: nio_slab(1), "p f p", False),
    "slab_z": (lambda: nio_slab(2), "p p f", False),
    "tilted_slab": (nio_tilted_slab, "p p f", True),
    "cluster": (nio_cluster, "f f f", False),
}


# --- routes -----------------------------------------------------------------


@pytest.fixture(scope="module", autouse=True)
def forget_the_lammps_module():
    """Unload the LAMMPS python package after this module's in-process runs.

    The python route imports ``lammps`` into this test process by design. Left
    in ``sys.modules``, it would make every later in-process "imports no
    backend" check (``backend_import_guard``) skip itself; removing it keeps
    those checks meaningful. The shared library stays loaded, and a later
    import simply binds to it again.
    """
    preloaded = {name for name in sys.modules if name == "lammps" or name.startswith("lammps.")}
    yield
    for name in [n for n in sys.modules if n == "lammps" or n.startswith("lammps.")]:
        if name not in preloaded:
            del sys.modules[name]


@pytest.fixture(params=["python", "executable"])
def route(request):
    if request.param == "python" and not HAVE_MODULE:
        pytest.skip("the LAMMPS python module is not installed")
    if request.param == "executable" and not HAVE_EXECUTABLE:
        pytest.skip("no 'lmp' executable on PATH")
    return request.param


@pytest.fixture
def module_route():
    if not HAVE_MODULE:
        pytest.skip("the LAMMPS python module is not installed")
    return "python"


def native(potential, route, tmp_path, **engine):
    options = {"workdir": str(tmp_path), **engine.pop("options", {})}
    return build_bridge(potential, EngineSpec(kind="lammps", runtime=route, options=options, **engine))


def md_spec(**kwargs):
    base = dict(task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=1.0, steps=20, seed=11,
                trajectory_interval=5, log_interval=5)
    base.update(kwargs)
    return SimulationSpec(**base)


def series_temperatures(workdir: Path):
    rows = np.loadtxt(Path(workdir) / lammps_engine.SERIES_FILE)
    return rows[:, lammps_engine.SERIES_KEYS.index("temp")]


# --- single points: every geometry, both routes -----------------------------


@pytest.mark.parametrize("geometry", sorted(GEOMETRIES))
def test_lj_single_point_matches_an_independent_reference(geometry, route, tmp_path):
    build, boundary, triclinic = GEOMETRIES[geometry]
    atoms = build()
    periodic = all(atoms.pbc)
    result = native(lj_potential(), route, tmp_path).singlepoint(
        atoms, SimulationSpec(task="singlepoint", compute_stress=periodic)
    )
    energy, forces, _, stress = reference_lj(atoms)
    assert result.energy_eV == pytest.approx(energy, abs=1e-10)
    assert np.abs(np.array(result.forces_eV_per_A) - forces).max() < 1e-10
    if periodic:
        # LAMMPS's own bar-per-eV/A^3 constant (nktv2p) limits this to ~1e-9 relative.
        assert np.abs(np.array(result.stress_eV_per_A3) - stress).max() < 1e-8 * max(1.0, np.abs(stress).max()) + 1e-10
    geometry_record = result.extras["geometry"]
    assert geometry_record["boundary"] == boundary
    assert geometry_record["triclinic"] is triclinic
    assert result.extras["position_roundtrip_max_error_angstrom"] < 1e-9
    assert result.extras["route"] == ("python-module" if route == "python" else "executable")


def test_the_executable_route_keeps_its_command_and_output(tmp_path):
    if not HAVE_EXECUTABLE:
        pytest.skip("no 'lmp' executable on PATH")
    result = native(lj_potential(), "executable", tmp_path).singlepoint(
        nio_cubic(), SimulationSpec(task="singlepoint")
    )
    launch = result.extras["launch"]
    workdir = tmp_path / "lammps_singlepoint"
    assert launch["argv"] == ["lmp", "-in", "in.lammps", "-log", "log.lammps"]
    assert launch["command_line"] == "lmp -in in.lammps -log log.lammps"
    assert (workdir / lammps_engine.STDOUT_FILE).read_text().strip()
    assert (workdir / "in.lammps").read_text().splitlines() == result.extras["deck"]


def test_an_executable_timeout_is_honoured(tmp_path):
    if not HAVE_EXECUTABLE:
        pytest.skip("no 'lmp' executable on PATH")
    bridge = native(lj_potential(), "executable", tmp_path, timeout_s=1e-4)
    with pytest.raises(MlipError, match="did not finish within engine.timeout_s"):
        bridge.singlepoint(nio_cubic(), SimulationSpec(task="singlepoint"))


def test_openmp_threads_are_applied_and_do_not_change_the_answer(route, tmp_path):
    atoms = nio_cubic()
    serial = native(lj_potential(), route, tmp_path / "serial").singlepoint(atoms, SimulationSpec(task="singlepoint"))
    threaded = native(lj_potential(), route, tmp_path / "omp", threads=2).singlepoint(
        atoms, SimulationSpec(task="singlepoint")
    )
    assert threaded.extras["launch"]["cmdargs"] == ["-pk", "omp", "2", "-sf", "omp"]
    assert threaded.extras["launch"]["env"] == {"OMP_NUM_THREADS": "2"}
    assert threaded.energy_eV == pytest.approx(serial.energy_eV, abs=1e-10)
    assert np.abs(np.array(threaded.forces_eV_per_A) - np.array(serial.forces_eV_per_A)).max() < 1e-10


def test_a_pre_existing_type_array_does_not_relabel_elements(route, tmp_path):
    atoms = nio_cubic()
    reference = reference_lj(atoms)[0]
    atoms.arrays["type"] = np.array([2 if s == "Ni" else 1 for s in atoms.get_chemical_symbols()])
    result = native(lj_potential(), route, tmp_path).singlepoint(atoms, SimulationSpec(task="singlepoint"))
    assert result.energy_eV == pytest.approx(reference, abs=1e-10)


def test_per_atom_energies_are_returned_only_when_they_sum_to_the_total(route, tmp_path, monkeypatch):
    # lj/cut tallies eatom; the framework table only lists MLIP pair styles.
    monkeypatch.setitem(lammps_mlip.KNOWN_FRAMEWORKS, "lj/cut", {"per_atom_energy": True, "stress": True})
    atoms = nio_primitive()
    result = native(lj_potential(), route, tmp_path).singlepoint(
        atoms, SimulationSpec(task="singlepoint", compute_per_atom_energy=True)
    )
    _, _, per_atom, _ = reference_lj(atoms)
    assert np.abs(np.array(result.per_atom_energy_eV) - per_atom).max() < 1e-10
    assert sum(result.per_atom_energy_eV) == pytest.approx(result.energy_eV, abs=1e-9)

    real_run_deck = lammps_engine.run_deck

    def incomplete(*args, **kwargs):
        raw = real_run_deck(*args, **kwargs)
        raw["per_atom_energy"] = raw["per_atom_energy"] * 0.5  # a style that tallies half
        return raw

    monkeypatch.setattr(lammps_engine, "run_deck", incomplete)
    with pytest.raises(ResultError, match="does not tally a complete per-atom energy"):
        native(lj_potential(), route, tmp_path / "bad").singlepoint(
            atoms, SimulationSpec(task="singlepoint", compute_per_atom_energy=True)
        )


# --- metal vs real ----------------------------------------------------------


def test_real_units_are_kcal_per_mol_so_the_same_numbers_scale_by_the_conversion(route, tmp_path):
    """The same pair_coeff numbers mean kcal/mol in real units: E, F and stress scale by one factor."""
    atoms = nio_primitive()
    sim = SimulationSpec(task="singlepoint", compute_stress=True)
    metal = native(lj_potential("metal"), route, tmp_path / "m").singlepoint(atoms, sim)
    real = native(lj_potential("real", scale_to_units=False), route, tmp_path / "r").singlepoint(atoms, sim)
    assert real.energy_eV / metal.energy_eV == pytest.approx(KCAL_PER_MOL_IN_EV, rel=1e-12)
    assert np.array(real.forces_eV_per_A) / KCAL_PER_MOL_IN_EV == pytest.approx(np.array(metal.forces_eV_per_A), rel=1e-10, abs=1e-14)
    # Stress passes through LAMMPS's per-unit-style pressure constants (bar vs atm).
    assert np.array(real.stress_eV_per_A3) / KCAL_PER_MOL_IN_EV == pytest.approx(np.array(metal.stress_eV_per_A3), rel=1e-6)
    assert real.native_units == "lammps_real" and metal.native_units == "lammps_metal"


def test_real_units_reproduce_the_metal_single_point(route, tmp_path):
    atoms = nio_primitive()
    sim = SimulationSpec(task="singlepoint", compute_stress=True)
    metal = native(lj_potential("metal"), route, tmp_path / "m").singlepoint(atoms, sim)
    real = native(lj_potential("real"), route, tmp_path / "r").singlepoint(atoms, sim)
    assert real.energy_eV == pytest.approx(metal.energy_eV, abs=1e-10)
    assert np.abs(np.array(real.forces_eV_per_A) - np.array(metal.forces_eV_per_A)).max() < 1e-10
    stress = np.array(metal.stress_eV_per_A3)
    assert np.abs(np.array(real.stress_eV_per_A3) - stress).max() < 1e-6 * np.abs(stress).max()


CROSS_UNIT_RUNS = {
    "nve": dict(ensemble="nve"),
    "nvt-langevin": dict(thermostat="langevin"),
    "nvt-nose-hoover": dict(thermostat="nose-hoover"),
    "npt-mtk": dict(ensemble="npt", pressure_bar=5000.0),
    "npt-berendsen-csvr": dict(ensemble="npt", pressure_bar=5000.0, barostat="berendsen", thermostat="csvr"),
}
BERENDSEN = {"barostat_bulk_modulus_bar": 1.9e6}


def cross_unit_run(units, route, workdir, case, *, pressure_scale=1.0):
    kwargs = dict(CROSS_UNIT_RUNS[case])
    if "pressure_bar" in kwargs:
        kwargs["pressure_bar"] *= pressure_scale
    options = BERENDSEN if kwargs.get("barostat") == "berendsen" else {}
    # press/berendsen needs an orthogonal box, so every case uses the cubic cell.
    return native(lj_potential(units), route, workdir, options=options).run_md(
        nio_cubic(), md_spec(**kwargs), workdir=workdir
    )


@pytest.mark.parametrize("case", sorted(CROSS_UNIT_RUNS))
def test_real_units_follow_the_metal_trajectory(case, route, tmp_path):
    """Timestep (ps vs fs), damping, velocities, and barostat targets (bar vs atm) all agree."""
    metal = cross_unit_run("metal", route, tmp_path / "metal", case)
    real = cross_unit_run("real", route, tmp_path / "real", case)
    positions = np.array(metal.final_positions_angstrom)
    cell = np.array(metal.final_cell_angstrom)
    # Measured: 1e-8 A (positions) and 6e-10 A (cells) -- LAMMPS's per-unit-style
    # constants agree to ~1e-8, nothing else differs.
    assert np.abs(np.array(real.final_positions_angstrom) - positions).max() < 1e-7
    cell_difference = np.abs(np.array(real.final_cell_angstrom) - cell).max()
    assert cell_difference < 1e-8
    assert real.temperature_end_K == pytest.approx(metal.temperature_end_K, rel=1e-6)
    assert np.abs(positions - nio_cubic().positions).max() > 1e-3  # the atoms did move
    if case.startswith("npt"):
        assert np.abs(cell - nio_cubic().cell.array).max() > 1e-4  # the cell did move
        # A target read as atm when it is bar (1.3 % too high) moves the cell
        # measurably differently: the agreement above is not a coincidence.
        slipped = cross_unit_run("real", route, tmp_path / "slip", case, pressure_scale=ATM_IN_BAR)
        slip = np.abs(np.array(slipped.final_cell_angstrom) - cell).max()
        assert slip > 1e-6 and slip > 100 * cell_difference


# --- thermostats and barostats run as rendered ------------------------------


@pytest.mark.parametrize("thermostat", ["nose-hoover", "langevin", "berendsen", "csvr"])
def test_every_thermostat_holds_the_temperature_nve_cannot(thermostat, route, tmp_path):
    """Started at 600 K, a near-harmonic crystal equipartitions to ~300 K under NVE."""
    atoms = nio_cubic()
    common = dict(temperature_K=600.0, steps=300, trajectory_interval=50, log_interval=5, seed=1)
    bridge = native(lj_potential(), route, tmp_path)
    nve = bridge.run_md(atoms, md_spec(ensemble="nve", **common), workdir=tmp_path / "nve")
    nvt = bridge.run_md(
        atoms, md_spec(thermostat=thermostat, thermostat_damping_fs=10.0, **common), workdir=tmp_path / "nvt"
    )
    nve_T = series_temperatures(tmp_path / "nve")
    nvt_T = series_temperatures(tmp_path / "nvt")
    assert nve_T[len(nve_T) // 2:].mean() < 420.0
    assert 480.0 < nvt_T[len(nvt_T) // 2:].mean() < 720.0
    assert nvt.temperature_start_K == pytest.approx(600.0, rel=1e-9)
    record = nvt.integrator_resolved
    assert record["thermostat"] == thermostat
    assert record["thermostat_damping_fs"] == 10.0 and record["thermostat_damping_lammps"] == 0.01
    expected = {"nose-hoover": ["nvt"], "langevin": ["nve", "langevin"],
                "berendsen": ["nve", "temp/berendsen"], "csvr": ["nve", "temp/csvr"]}[thermostat]
    assert record["fixes"] == expected
    deck = (tmp_path / "nvt" / "in.lammps").read_text()
    assert f"fix mlip_integrate all {expected[0]}" in deck
    # The thermostat's energy exchange is tallied (ecouple), so econserve is flat
    # where etotal is not.
    diagnostics = nvt.diagnostics
    assert diagnostics["conserved_quantity"] == "LAMMPS econserve (pe + ke + ecouple)"
    excursion = diagnostics["conserved_quantity_drift"]["max_abs_energy_excursion_eV_per_atom"]
    assert excursion < 0.25 * abs(diagnostics["total_energy_change_eV_per_atom"])
    assert nvt.steps_completed == 300 and nvt.frames_written == 7


def test_nve_conserves_energy_with_second_order_error(route, tmp_path):
    """Velocity Verlet: the energy excursion shrinks ~4x when the timestep halves.

    The LJ is shifted to zero at the cutoff (``pair_modify shift yes``), so the
    energy has no jumps when pairs cross it and only integration error remains.
    """
    from dataclasses import replace

    shifted = replace(lj_potential(), extra_commands=("pair_modify shift yes",))
    excursions = []
    for dt in (1.0, 0.5):
        result = native(shifted, route, tmp_path / str(dt)).run_md(
            nio_cubic(), md_spec(ensemble="nve", steps=200, timestep_fs=dt), workdir=tmp_path / str(dt)
        )
        assert result.diagnostics["conservation_test"] is True
        assert result.integrator_resolved["fixes"] == ["nve"]
        excursions.append(result.diagnostics["max_abs_energy_excursion_eV_per_atom"])
    assert excursions[0] < 1e-4  # measured 4e-6 eV/atom at 1 fs
    assert 2.5 < excursions[0] / excursions[1] < 6.0  # measured 4.0


def test_md_endpoint_forces_are_the_potential_alone(route, tmp_path):
    """Langevin friction and noise are removed before the final single point."""
    result = native(lj_potential(), route, tmp_path).run_md(
        nio_cubic(), md_spec(thermostat="langevin", thermostat_damping_fs=5.0), workdir=tmp_path
    )
    final = nio_cubic()
    final.set_cell(result.final_cell_angstrom)
    final.positions = result.final_positions_angstrom
    _, forces, _, _ = reference_lj(final)
    assert np.abs(np.array(result.final.forces_eV_per_A) - forces).max() < 1e-9


def test_the_trajectory_is_written_in_the_source_basis(route, tmp_path):
    from ase.io import read

    atoms = nio_rotated()
    result = native(lj_potential(), route, tmp_path).run_md(
        atoms, md_spec(ensemble="nve", steps=12, trajectory_interval=5), workdir=tmp_path
    )
    frames = read(result.trajectory_path, index=":")
    assert [frame.info["timestep"] for frame in frames] == [0, 5, 10, 12]
    assert result.frames_written == 4
    # extxyz stores 8 decimals; the rotation back is exact to rounding.
    assert np.abs(frames[0].positions - atoms.positions).max() < 1e-7
    assert np.abs(frames[0].cell.array - atoms.cell.array).max() < 1e-7
    assert np.abs(frames[-1].positions - np.array(result.final_positions_angstrom)).max() < 1e-7
    assert frames[0].get_chemical_symbols() == atoms.get_chemical_symbols()


@pytest.mark.parametrize(
    ("builder", "coupling", "fixed_axis"),
    [
        (nio_cubic, None, None),
        (nio_cubic, "anisotropic", None),
        (lambda: nio_slab(2), "in-plane", 2),
        (lambda: nio_slab(0), "in-plane", 0),
    ],
)
def test_barostats_move_only_the_dimensions_they_should(builder, coupling, fixed_axis, route, tmp_path):
    atoms = builder()
    result = native(lj_potential(), route, tmp_path).run_md(
        atoms, md_spec(ensemble="npt", pressure_bar=20000.0, barostat_coupling=coupling, steps=100,
                       barostat_damping_fs=200.0),
        workdir=tmp_path,
    )
    before = atoms.cell.lengths()
    after = np.array(result.final_cell_angstrom)
    lengths = np.linalg.norm(after, axis=1)
    change = lengths - before
    periodic = [axis for axis in range(3) if axis != fixed_axis]
    assert all(change[axis] < -1e-3 for axis in periodic)  # compressed along every barostatted axis
    if fixed_axis is not None:
        assert abs(change[fixed_axis]) < 1e-9  # the open axis is never barostatted
        assert change[periodic[0]] == pytest.approx(change[periodic[1]], rel=1e-9)  # coupled
    if coupling is None:
        assert np.ptp(change) < 1e-9  # isotropic
    assert result.integrator_resolved["barostat"] == "mtk"


def test_frozen_atoms_do_not_move_and_are_excluded_from_the_temperature(route, tmp_path):
    atoms = nio_slab(2)
    bottom = [i for i, z in enumerate(atoms.positions[:, 2]) if z < atoms.positions[:, 2].min() + 1.0]
    atoms.set_constraint(FixAtoms(indices=bottom))
    result = native(lj_potential(), route, tmp_path).run_md(
        atoms, md_spec(thermostat="langevin", steps=40), workdir=tmp_path
    )
    final = np.array(result.final_positions_angstrom)
    mobile = [i for i in range(len(atoms)) if i not in bottom]
    assert np.abs(final[bottom] - atoms.positions[bottom]).max() < 1e-9
    assert np.abs(final[mobile] - atoms.positions[mobile]).max() > 1e-3
    assert result.temperature_ndof == 3 * len(mobile)
    assert result.temperature_start_K == pytest.approx(300.0, rel=1e-9)  # velocities on mobile atoms only
    assert result.constraints["fixed_atoms"] == bottom


# --- model files -------------------------------------------------------------


def write_lj_table(path: Path, n: int = 2000):
    """LJ as a LAMMPS pair table (energy and force columns), one section per pair."""
    lines = ["# LJ mixture as tables"]
    for (a, b), (eps, sig) in LJ.items():
        lines += ["", f"{a}{b}", f"N {n} R 1.0 {CUTOFF}", ""]
        for k, r in enumerate(np.linspace(1.0, CUTOFF, n).tolist(), start=1):
            sr6 = (sig / r) ** 6
            lines.append(f"{k} {r!r} {4 * eps * (sr6 * sr6 - sr6)!r} {4 * eps * (12 * sr6 * sr6 - 6 * sr6) / r!r}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def table_potential(table: Path) -> LammpsMlipPotentialSpec:
    return LammpsMlipPotentialSpec(
        label="lj-table",
        pair_style="table linear 2000",
        pair_coeff=tuple(
            # Quoted: LAMMPS splits arguments on whitespace, and this path has a space.
            f'{1 if a == "Ni" else 2} {1 if b == "Ni" else 2} "{table.as_posix()}" {a}{b} {CUTOFF}'
            for (a, b) in LJ
        ),
        type_map={1: "Ni", 2: "O"},
        model_paths=(table,),
    )


def test_model_files_are_staged_hashed_and_named_for_each_lammps(module_route, tmp_path):
    """Native LAMMPS opens the staged name in the job directory; LAMMPSlib the absolute path."""
    from nio_md_prep.mlip.provenance import sha256_file

    table = write_lj_table(tmp_path / "my models" / "lj.table")
    potential = table_potential(table)
    atoms = nio_cubic()
    sim = SimulationSpec(task="singlepoint", compute_stress=True)
    by_native = native(potential, "python", tmp_path / "native").singlepoint(atoms, sim)
    job_dir = tmp_path / "ase job"
    lammpslib = build_bridge(potential, EngineSpec(kind="ase", options={"workdir": str(job_dir)}))
    by_ase = lammpslib.singlepoint(atoms, sim)

    staged_native = by_native.extras["staged_model_files"][0]
    staged_ase = by_ase.extras["staged_model_files"][0]
    for record in (staged_native, staged_ase):
        assert record["sha256"] == sha256_file(table)
        assert record["referenced_by_pair_commands"] is True
        assert Path(record["staged"]).name == "lj.table"
    assert all(" lj.table " in c for c in by_native.extras["pair_commands"][1:])
    pair_commands = lammpslib.calculator(job_dir / lammps_engine_singlepoint_dir()).parameters.lmpcmds
    absolute = Path(staged_ase["staged"]).resolve().as_posix()
    assert any(f'"{absolute}"' in c for c in pair_commands)  # quoted: the path has a space
    assert (job_dir / lammps_engine_singlepoint_dir() / "lammpslib.log").is_file()
    # Both LAMMPS instances read the same table; the table is close to the analytic LJ.
    assert by_ase.energy_eV == pytest.approx(by_native.energy_eV, abs=1e-9)
    assert np.abs(np.array(by_ase.forces_eV_per_A) - np.array(by_native.forces_eV_per_A)).max() < 1e-9
    assert np.abs(np.array(by_ase.stress_eV_per_A3) - np.array(by_native.stress_eV_per_A3)).max() < 1e-9
    assert by_native.energy_eV == pytest.approx(reference_lj(atoms)[0], abs=1e-3)


def lammps_engine_singlepoint_dir() -> str:
    from nio_md_prep.mlip.bridges.lammps_ase import SINGLEPOINT_DIR

    return SINGLEPOINT_DIR


def test_a_declared_model_hash_that_does_not_match_is_refused(module_route, tmp_path):
    from dataclasses import replace

    from nio_md_prep.mlip.errors import ModelIntegrityError

    table = write_lj_table(tmp_path / "lj.table")
    potential = replace(table_potential(table), model_hashes={str(table): "0" * 64})
    with pytest.raises(ModelIntegrityError):
        native(potential, "python", tmp_path / "job").singlepoint(nio_cubic(), SimulationSpec(task="singlepoint"))


# --- LAMMPSlib (ASE) route ----------------------------------------------------


def test_lammpslib_matches_the_native_route_including_stress_and_per_atom(module_route, tmp_path, monkeypatch):
    monkeypatch.setitem(lammps_mlip.KNOWN_FRAMEWORKS, "lj/cut", {"per_atom_energy": True, "stress": True})
    atoms = nio_primitive()
    sim = SimulationSpec(task="singlepoint", compute_stress=True, compute_per_atom_energy=True)
    by_native = native(lj_potential(), "python", tmp_path / "n").singlepoint(atoms, sim)
    by_ase = build_bridge(lj_potential(), EngineSpec(kind="ase", options={"workdir": str(tmp_path / "a")})).singlepoint(atoms, sim)
    assert by_ase.energy_eV == pytest.approx(by_native.energy_eV, abs=1e-10)
    assert np.abs(np.array(by_ase.forces_eV_per_A) - np.array(by_native.forces_eV_per_A)).max() < 1e-10
    assert np.abs(np.array(by_ase.stress_eV_per_A3) - np.array(by_native.stress_eV_per_A3)).max() < 1e-9
    assert np.abs(np.array(by_ase.per_atom_energy_eV) - np.array(by_native.per_atom_energy_eV)).max() < 1e-10


def test_lammpslib_md_runs_through_the_ase_engine(module_route, tmp_path):
    bridge = build_bridge(lj_potential(), EngineSpec(kind="ase"))
    result = bridge.run_md(nio_cubic(), md_spec(thermostat="langevin"), workdir=tmp_path)
    assert result.steps_completed == 20
    assert result.frames_written == 5
    assert (tmp_path / "lammpslib.log").is_file()


# --- MACE routes on a LAMMPS without the packages they need ------------------


def fake_mace(tmp_path, implementation):
    from nio_md_prep.mlip.bridges.mace_lammps import EXPORT_SUFFIX

    checkpoint = tmp_path / "nio.model"
    checkpoint.write_bytes(b"checkpoint")
    (tmp_path / ("nio.model" + EXPORT_SUFFIX[implementation])).write_bytes(b"export")
    return MacePotentialSpec(label="nio", model_path=checkpoint, declared_elements=("Ni", "O"))


@pytest.mark.parametrize(
    ("implementation", "options", "missing"),
    [
        ("mliap", {}, ("package PYTHON", "package KOKKOS")),
        ("mliap", {"mliap_coupling": "plain"}, ("package PYTHON",)),
        ("pair-mace", {}, ("pair style mace",)),
    ],
)
def test_mace_on_a_lammps_without_its_packages_fails_clearly(implementation, options, missing, route, tmp_path):
    build = lammps_engine.probe_build(route, "lmp" if route == "executable" else None)
    if "PYTHON" in build["packages"] or "mace" in build["pair_styles"]:
        pytest.skip("this LAMMPS build has the PYTHON package or pair mace")
    bridge = build_bridge(
        fake_mace(tmp_path, implementation),
        EngineSpec(kind="lammps", runtime=route, options=options),
        implementation=implementation,
    )
    availability = bridge.availability()
    assert not availability
    for name in missing:
        assert name in availability.detail
    with pytest.raises(MissingDependencyError, match="lacks"):
        bridge.singlepoint(nio_cubic(), SimulationSpec(task="singlepoint"))
    plan = bridge.execution_plan(SimulationSpec(task="singlepoint"), nio_cubic())
    assert plan["availability"]["available"] is False


# --- MACE's own exporter (mace-torch) -------------------------------------------

TINY_MODEL = r"""
import sys, numpy as np, torch
from e3nn import o3
from mace import modules
torch.set_default_dtype(torch.float64)
model = modules.ScaleShiftMACE(
    r_max=4.0, num_bessel=4, num_polynomial_cutoff=5, max_ell=1,
    interaction_cls=modules.interaction_classes["RealAgnosticResidualInteractionBlock"],
    interaction_cls_first=modules.interaction_classes["RealAgnosticInteractionBlock"],
    num_interactions=2, num_elements=2, hidden_irreps=o3.Irreps("8x0e"),
    MLP_irreps=o3.Irreps("8x0e"), atomic_energies=np.array([-1.0, -2.0]),
    avg_num_neighbors=8.0, atomic_numbers=[8, 28], correlation=2, gate=torch.nn.functional.silu,
    atomic_inter_scale=1.0, atomic_inter_shift=0.0)
torch.save(model, sys.argv[1])
"""
#: mace_create_lammps_model's own main(); only the e3nn->cuequivariance conversion
#: (which needs the optional cuequivariance package) is replaced by the identity,
#: so the file name and the pickled flags are exactly what the exporter writes.
EXPORT_MLIAP = r"""
import sys
import mace.cli.create_lammps_model as exporter
exporter.run_e3nn_to_cueq = lambda model: model
sys.argv = ["mace_create_lammps_model", sys.argv[1], "--format", "mliap", "--dtype", sys.argv[2]]
exporter.main()
"""


def _mace_env(**extra):
    env = {**os.environ, "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD": "1", **extra}
    for name in ("MACE_ALLOW_CPU", "MACE_FORCE_CPU"):
        if name not in extra:
            env.pop(name, None)
    return env


def _run(args, **env):
    completed = subprocess.run(
        [sys.executable, *args], capture_output=True, text=True, timeout=900, env=_mace_env(**env), check=False
    )
    assert completed.returncode == 0, completed.stderr[-2000:]
    return completed


@pytest.fixture(scope="module")
def tiny_mace(tmp_path_factory):
    if find_spec("mace") is None or find_spec("torch") is None:
        pytest.skip("mace-torch and torch are not installed")
    directory = tmp_path_factory.mktemp("mace export")
    checkpoint = directory / "nio_tiny.model"
    _run(["-c", TINY_MODEL, str(checkpoint)])
    return checkpoint


@pytest.mark.mace
def test_the_libtorch_export_lands_where_the_pair_mace_route_looks(tiny_mace):
    before = set(tiny_mace.parent.iterdir())
    _run(["-m", "mace.cli.create_lammps_model", str(tiny_mace), "--dtype", "float32"])
    written = sorted(set(tiny_mace.parent.iterdir()) - before)
    spec = MacePotentialSpec(model_path=tiny_mace, declared_elements=("Ni", "O"), precision="float32")
    bridge = build_bridge(spec, EngineSpec(kind="lammps"), implementation="pair-mace")
    assert written == [bridge.exported_model_path()]
    assert written[0].name == "nio_tiny.model-lammps.pt"
    record = bridge.export_record()["exported_model"]
    assert record["exists"] and len(record["sha256"]) == 64
    info = bridge.check_export()
    assert info["readable"] and info["dtype"] == "float32"
    assert sorted(info["elements"]) == ["Ni", "O"] and info["num_interactions"] == 2
    double = build_bridge(
        MacePotentialSpec(model_path=tiny_mace, declared_elements=("Ni", "O"), precision="float64"),
        EngineSpec(kind="lammps"),
        implementation="pair-mace",
    )
    with pytest.raises(Exception, match="float32 model"):
        double.check_export()


@pytest.mark.mace
def test_the_mliap_export_lands_where_the_mliap_route_looks_and_carries_its_cpu_flag(tiny_mace, tmp_path):
    from nio_md_prep.mlip.errors import ConfigError

    checkpoint = tmp_path / "nio.model"
    shutil.copy2(tiny_mace, checkpoint)
    spec = MacePotentialSpec(model_path=checkpoint, declared_elements=("Ni", "O"))
    bridge = build_bridge(spec, EngineSpec(kind="lammps"))
    # Exported without MACE_ALLOW_CPU: the flag is pickled as False, and the
    # CPU KOKKOS route refuses it whatever the run environment would say.
    _run(["-c", EXPORT_MLIAP, str(checkpoint), "float64"])
    assert bridge.exported_model_path() == tmp_path / "nio.model-mliap_lammps.pt"
    assert bridge.exported_model_path().is_file()
    info = mace_inspect(bridge)
    assert info["readable"] and info["dtype"] == "float64" and info["allow_cpu"] is False
    with pytest.raises(ConfigError, match="MACE_ALLOW_CPU"):
        bridge.check_export()
    # Plain (single-layer-only) coupling refuses the two-layer model.
    plain = build_bridge(spec, EngineSpec(kind="lammps", options={"mliap_coupling": "plain"}))
    with pytest.raises(ConfigError, match="2 interaction layers"):
        plain.check_export()
    # Re-exported with MACE_ALLOW_CPU=true: stored, and accepted.
    _run(["-c", EXPORT_MLIAP, str(checkpoint), "float64"], MACE_ALLOW_CPU="true")
    info = mace_inspect(bridge)
    assert info["allow_cpu"] is True
    assert bridge.check_export()["num_interactions"] == 2
    plan = bridge.execution_plan(SimulationSpec(task="singlepoint"), None)
    assert plan["precision"] == {
        "requested": "float64", "effective": "float64", "guaranteed": True,
        "note": plan["precision"]["note"],
    }
    assert plan["exported_model"]["inspection"]["allow_cpu"] is True


def mace_inspect(bridge):
    from nio_md_prep.mlip.bridges import mace_lammps

    return mace_lammps.inspect_export(bridge.exported_model_path(), bridge.implementation)
