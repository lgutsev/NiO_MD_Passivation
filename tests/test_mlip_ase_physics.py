"""The ASE engine's physics, checked against expectations that mean something.

Runs in ordinary CI (ASE + numpy only): the mock Lennard-Jones potential and
ASE's EMT for copper stand in for a trained model. Each expectation is
physical rather than a snapshot:

* single-point forces are the calculator's raw forces, frozen atoms included,
  and match an independent brute-force numpy Lennard-Jones sum;
* Langevin dynamics keeps FixAtoms atoms exactly in place and holds the free
  atoms at the target temperature within a stated statistical window;
* NVE drift is a conservation test whose velocity-Verlet error scales as
  dt^2; a thermostatted run reports only a descriptive energy change, plus
  the integrator's own conserved quantity where it has one;
* the requested thermostat/barostat is what integrates, or the request is
  refused before anything runs;
* step numbers and frame counts come from the run and the file, and a missing
  or corrupt trajectory is a failure.
"""
from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest

from nio_md_prep.mlip import jobs
from nio_md_prep.mlip.errors import ConfigError, MlipError, ResultError
from nio_md_prep.mlip.specs import (
    MAX_SEED,
    EngineSpec,
    JobSpec,
    MockPotentialSpec,
    SimulationSpec,
    StructureSpec,
)

pytest.importorskip("ase")
np = pytest.importorskip("numpy")

from ase.constraints import FixAtoms  # noqa: E402

from nio_md_prep.mlip._mock_calculator import MockLennardJones  # noqa: E402
from nio_md_prep.mlip.engines import ase_engine  # noqa: E402

#: Near the LJ minimum for the NiO nearest-neighbour distance (2.085 A =
#: 2^(1/6) sigma), so a 300 K run stays a solid instead of exploding.
EPSILON, SIGMA, CUTOFF = 0.05, 1.86, 6.0


def mock_calculator():
    return MockLennardJones(epsilon_eV=EPSILON, sigma_angstrom=SIGMA, cutoff_angstrom=CUTOFF)


def mock_job(structure: Path, *, elements=("Ni", "O"), sigma=SIGMA, engine=None, **simulation):
    defaults = dict(task="md", timestep_fs=1.0, steps=40, seed=7, trajectory_interval=10,
                    log_interval=5)
    defaults.update(simulation)
    return JobSpec(
        potential=MockPotentialSpec(
            label="mock", declared_elements=tuple(elements), epsilon_eV=EPSILON,
            sigma_angstrom=sigma, cutoff_angstrom=CUTOFF,
        ),
        engine=engine or EngineSpec(kind="ase"),
        simulation=SimulationSpec(**defaults),
        structure=StructureSpec(path=structure),
    )


def write_structure(tmp_path: Path, atoms, name="structure.xyz") -> Path:
    from ase.io import write

    path = tmp_path / name
    write(str(path), atoms, format="extxyz")
    return path


@pytest.fixture
def pinned_nio(rattled_nio_structure):
    """The rattled NiO cell with its first four atoms frozen."""
    atoms = rattled_nio_structure.copy()
    atoms.set_constraint(FixAtoms(indices=[0, 1, 2, 3]))
    return atoms


# --- single points: raw forces ---------------------------------------------


def reference_lj(atoms, epsilon=EPSILON, sigma=SIGMA, cutoff=CUTOFF):
    """Shifted LJ energy and forces by brute force over periodic images, in numpy.

    Independent of the calculator under test: no neighbour list, every image
    in a +-2 cell block, each unordered pair counted once.
    """
    positions = atoms.get_positions()
    cell = np.array(atoms.get_cell())
    images = [np.array([i, j, k]) @ cell for i in range(-2, 3) for j in range(-2, 3)
              for k in range(-2, 3)]
    shift = 4 * epsilon * ((sigma / cutoff) ** 12 - (sigma / cutoff) ** 6)
    energy, forces = 0.0, np.zeros_like(positions)
    for i in range(len(atoms)):
        for j in range(len(atoms)):
            for offset in images:
                d = positions[j] + offset - positions[i]
                r = float(np.linalg.norm(d))
                if r < 1e-8 or r >= cutoff:
                    continue
                sr6 = (sigma / r) ** 6
                energy += 0.5 * (4 * epsilon * (sr6 * sr6 - sr6) - shift)
                dEdr = 4 * epsilon * (-12 * sr6 * sr6 + 6 * sr6) / r
                forces[i] += dEdr * d / r
    return energy, forces


def test_single_point_forces_are_raw_on_frozen_atoms_and_match_a_numpy_reference(pinned_nio):
    payload = ase_engine.singlepoint(pinned_nio, mock_calculator())
    energy, forces = reference_lj(pinned_nio)
    reported = np.array(payload["forces_eV_per_A"])
    assert payload["energy_eV"] == pytest.approx(energy, rel=1e-12)
    np.testing.assert_allclose(reported, forces, rtol=0, atol=1e-10)
    # The frozen rows carry real forces, not the zeros get_forces() would give.
    assert np.linalg.norm(reported[:4], axis=1).min() > 1e-3
    constrained = pinned_nio.copy()
    constrained.calc = mock_calculator()
    assert np.abs(constrained.get_forces()[:4]).max() == 0.0  # ASE's default, for contrast
    native = payload["native"]
    assert native["constraints"]["fixed_atoms"] == [0, 1, 2, 3]
    assert native["constraints"]["applied_to_forces"] is False
    assert native["max_force_free_atoms_eV_per_A"] == pytest.approx(
        np.linalg.norm(forces[4:], axis=1).max(), abs=1e-10
    )
    assert payload["trajectory_path"] is None


def test_the_bridge_records_the_constraints_of_a_single_point(tmp_path, pinned_nio):
    from ase.io import read

    path = write_structure(tmp_path, pinned_nio)
    job = mock_job(path, task="singlepoint", timestep_fs=None, steps=0, seed=None)
    result = jobs.run_singlepoint(job, output_dir=tmp_path / "sp")["result"]
    assert result.extras["constraints"]["n_fixed"] == 4
    assert result.max_force_eV_per_A > 0
    reloaded = read(str(path))  # extxyz keeps FixAtoms and 1e-8 A positions
    assert [type(c).__name__ for c in reloaded.constraints] == ["FixAtoms"]
    _, forces = reference_lj(reloaded)
    np.testing.assert_allclose(np.array(result.forces_eV_per_A), forces, atol=1e-10)


# --- frozen atoms under a thermostat -----------------------------------------


def test_langevin_holds_frozen_atoms_and_thermalises_the_free_ones(tmp_path):
    """Cu(111) EMT slab, bottom half frozen: 32 free atoms, 96 degrees of freedom.

    Window: the canonical per-sample spread of T is T*sqrt(2/96) ~ 43 K; the
    400-sample time average (after 0.4 ps of equilibration, damping 20 fs)
    scattered by ~3 K over five seeds, and a 4-5 fs Langevin step biases the
    kinetic temperature low by ~1-2 % (measured against 2 fs). +-5 % covers
    both with margin; a thermostat that ignored FixAtoms (ASE's NHC reads
    ~0 K here) or counted 3N (192 DOF, ~150 K) fails by far more.
    """
    from ase.build import fcc111
    from ase.calculators.emt import EMT
    from ase.io.trajectory import Trajectory
    from ase.units import kB

    slab = fcc111("Cu", size=(4, 4, 4), vacuum=6.0)
    frozen = np.flatnonzero(slab.get_tags() >= 3)  # the two bottom layers
    slab.set_constraint(FixAtoms(indices=frozen))
    simulation = SimulationSpec(
        task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=4.0, steps=500, seed=11,
        thermostat="langevin", thermostat_damping_fs=20.0, trajectory_interval=5, log_interval=1,
    )
    payload = ase_engine.run_md(slab, EMT(), simulation, workdir=tmp_path)

    free = np.setdiff1d(np.arange(len(slab)), frozen)
    assert payload["temperature_ndof"] == 3 * len(free) == 96
    assert payload["constraints"]["fixed_atoms"] == frozen.tolist()
    rows = np.loadtxt(payload["log_path"])
    assert rows[:, 0].tolist() == list(range(501))
    mean_temperature = rows[100:, 5].mean()
    assert 285.0 < mean_temperature < 315.0

    initial = slab.get_positions()[frozen]
    with Trajectory(payload["trajectory_path"]) as frames:
        assert len(frames) == 500 // 5 + 1 == payload["frames_written"]
        kinetic_free = []
        for frame in frames:
            # Exactly in place, not merely close: FixAtoms resets the positions.
            assert np.array_equal(frame.get_positions()[frozen], initial)
            momenta, masses = frame.get_momenta(), frame.get_masses()
            assert np.array_equal(momenta[frozen], np.zeros((len(frozen), 3)))
            kinetic_free.append(0.5 * np.sum(momenta[free] ** 2 / masses[free, None]))
            assert [type(c).__name__ for c in frame.constraints] == ["FixAtoms"]
    # The temperature column is the free atoms' kinetic temperature.
    t_free = 2 * np.array(kinetic_free) / (3 * len(free) * kB)
    np.testing.assert_allclose(t_free, rows[::5, 5], rtol=0, atol=1e-5)  # log keeps 6 decimals
    resolved = payload["integrator_resolved"]
    assert resolved["integrator"] == "Langevin" and resolved["fixcm"] is False
    assert resolved["friction_per_fs"] == pytest.approx(1 / 20.0)
    assert resolved["centre_of_mass"].startswith("none")


def test_nose_hoover_with_frozen_atoms_is_refused_before_anything_runs(tmp_path, pinned_nio):
    job = mock_job(write_structure(tmp_path, pinned_nio), ensemble="nvt", temperature_K=300.0,
                   thermostat="nose-hoover")
    with pytest.raises(ConfigError, match="nose-hoover.*frozen"):
        jobs.validate_job(job)
    with pytest.raises(ConfigError, match="nose-hoover"):
        jobs.run_smoke_md(job, output_dir=tmp_path / "md")
    assert not (tmp_path / "md" / "mlip_manifest.json").exists()
    # The engine refuses on its own too, for bridges that skip the hook.
    with pytest.raises(ConfigError, match="0 K"):
        ase_engine.run_md(pinned_nio, mock_calculator(), job.simulation, workdir=tmp_path / "e")


def test_langevin_and_csvr_with_frozen_atoms_are_accepted(tmp_path, pinned_nio):
    path = write_structure(tmp_path, pinned_nio)
    for thermostat in ("langevin", "csvr", "berendsen"):
        job = mock_job(path, ensemble="nvt", temperature_K=300.0, thermostat=thermostat)
        assert jobs.validate_job(job)["ok"]


# --- thermostats, conserved quantities, diagnostics -------------------------


@pytest.mark.parametrize(
    "thermostat, integrator, conserved",
    [
        (None, "Langevin", None),
        ("langevin", "Langevin", None),
        ("berendsen", "NVTBerendsen", None),
        ("csvr", "Bussi", "Bussi: total energy minus transferred_energy"),
        ("nose-hoover", "NoseHooverChainNVT", "NoseHooverChainNVT.get_conserved_energy"),
    ],
)
def test_the_requested_thermostat_is_what_integrates(
    tmp_path, rattled_nio_structure, thermostat, integrator, conserved
):
    job = mock_job(write_structure(tmp_path, rattled_nio_structure), ensemble="nvt",
                   temperature_K=300.0, thermostat=thermostat, steps=100)
    report = jobs.run_smoke_md(job, output_dir=tmp_path / "md")
    trajectory = report["trajectory"]
    resolved = trajectory.integrator_resolved
    assert resolved["integrator"] == integrator
    assert resolved["ase_todict"]["md-type"] == integrator
    assert resolved["thermostat_damping_fs"] == 100.0  # the shared default
    assert ("thermostat" in resolved["defaults_applied"]) is (thermostat is None)
    diagnostics = trajectory.diagnostics
    # A thermostat exchanges energy with a bath: the total-energy change is
    # descriptive, and never reported under a drift name.
    assert diagnostics["label"] == "descriptive; not conserved under a thermostat"
    assert "energy_drift_eV_per_atom_per_ps" not in diagnostics
    if conserved is None:
        assert diagnostics["conservation_test"] is False
        assert "conserved_quantity" not in diagnostics
    else:
        # The extended-system energy is conserved to the integrator's accuracy:
        # measured ~1e-5 eV/atom/ps against a descriptive change of ~1e-2 eV/atom.
        assert diagnostics["conservation_test"] is True
        assert diagnostics["conserved_quantity"] == conserved
        drift = diagnostics["conserved_quantity_drift"]
        assert abs(drift["energy_drift_eV_per_atom_per_ps"]) < 1e-3
        assert drift["max_abs_energy_excursion_eV_per_atom"] < 1e-4
    manifest = report["results"]["trajectory"]
    assert manifest["integrator_resolved"]["integrator"] == integrator
    # NHC thermostats 3N; the constraint-aware integrators run with FixCom.
    assert trajectory.temperature_ndof == (48 if integrator == "NoseHooverChainNVT" else 45)


def test_bussi_holds_the_target_temperature(tmp_path, rattled_nio_structure):
    """Bussi rescales towards the canonical kinetic energy over 3N-3 DOF; it does not heat."""
    payload = ase_engine.run_md(
        rattled_nio_structure, mock_calculator(),
        SimulationSpec(task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=1.0,
                       steps=500, seed=3, thermostat="csvr", thermostat_damping_fs=20.0,
                       trajectory_interval=50, log_interval=1),
        workdir=tmp_path,
    )
    rows = np.loadtxt(payload["log_path"])
    # 45 DOF: per-sample spread 300*sqrt(2/45) ~ 63 K; 0.4 ps of a 20 fs
    # thermostat averages it to a few K. NVE from the same start ends ~365 K.
    assert 270.0 < rows[100:, 5].mean() < 330.0


def test_nve_drift_is_a_per_ps_conservation_test_that_scales_as_dt_squared(
    tmp_path, rattled_nio_structure
):
    """Velocity Verlet's energy error is O(dt^2): halving dt quarters the excursion."""
    excursions = {}
    for dt in (1.0, 2.0):
        job = mock_job(write_structure(tmp_path, rattled_nio_structure), ensemble="nve",
                       temperature_K=300.0, timestep_fs=dt, steps=int(200 / dt), log_interval=2)
        diagnostics = jobs.run_smoke_md(job, output_dir=tmp_path / f"dt{dt}")[
            "trajectory"].diagnostics
        assert diagnostics["conservation_test"] is True and diagnostics["duration_fs"] == 200.0
        assert abs(diagnostics["energy_drift_eV_per_atom_per_ps"]) < 1e-3
        excursions[dt] = diagnostics["max_abs_energy_excursion_eV_per_atom"]
    # Measured 4.0 (4.2e-6 -> 1.7e-5 eV/atom).
    assert 3.0 < excursions[2.0] / excursions[1.0] < 5.0


def test_nve_refuses_a_thermostat():
    with pytest.raises(ConfigError, match="nve"):
        SimulationSpec(task="md", ensemble="nve", timestep_fs=1.0, steps=10, thermostat="langevin")


# --- steps, frames, trajectories --------------------------------------------


def test_step_numbers_and_frames_come_from_the_run(tmp_path, rattled_nio_structure):
    """23 steps, frames every 5 and log every 5: steps 0,5,10,15,20 and the final 23 -- once each."""
    payload = ase_engine.run_md(
        rattled_nio_structure, mock_calculator(),
        SimulationSpec(task="md", ensemble="nve", temperature_K=300.0, timestep_fs=1.0,
                       steps=23, seed=5, trajectory_interval=5, log_interval=5),
        workdir=tmp_path,
    )
    assert payload["steps_completed"] == 23
    assert payload["frames_written"] == 23 // 5 + 1 == 5
    steps = np.loadtxt(payload["log_path"])[:, 0].astype(int).tolist()
    assert steps == [0, 5, 10, 15, 20, 23]
    assert payload["energy_series"]["time_fs"] == [0.0, 5.0, 10.0, 15.0, 20.0, 23.0]
    # The end temperature is the final state's, not the last logged row's.
    final = payload["atoms"]
    assert payload["temperature_end_K"] == pytest.approx(
        2 * final.get_kinetic_energy() / (45 * 8.617333262e-5), rel=1e-6
    )
    np.testing.assert_array_equal(payload["final_positions"], final.get_positions())
    assert final.constraints == []  # the run's FixCom is not handed back


def test_a_missing_or_corrupt_trajectory_is_a_failure(tmp_path, rattled_nio_structure):
    payload = ase_engine.run_md(
        rattled_nio_structure, mock_calculator(),
        SimulationSpec(task="md", ensemble="nve", timestep_fs=1.0, steps=10,
                       trajectory_interval=5),
        workdir=tmp_path / "ok",
    )
    good = Path(payload["trajectory_path"])
    assert ase_engine.verify_frames(good, 10, 5) == 3
    with pytest.raises(ResultError, match="should give 4"):
        ase_engine.verify_frames(good, 15, 5)  # a short file is not a complete run
    with pytest.raises(ResultError, match="was not written"):
        ase_engine.verify_frames(tmp_path / "absent.traj", 10, 5)
    garbage = tmp_path / "garbage.traj"
    garbage.write_bytes(b"this is not an ASE trajectory" * 10)
    with pytest.raises(ResultError, match="cannot be read"):
        ase_engine.verify_frames(garbage, 10, 5)
    truncated = tmp_path / "truncated.traj"
    truncated.write_bytes(good.read_bytes()[: good.stat().st_size // 2])
    with pytest.raises(ResultError):
        ase_engine.verify_frames(truncated, 10, 5)
    assert issubclass(ResultError, MlipError)


def test_a_run_whose_trajectory_comes_up_short_fails_with_a_failed_manifest(
    tmp_path, rattled_nio_structure, monkeypatch
):
    from ase.io import trajectory as ase_trajectory

    original = ase_trajectory.TrajectoryWriter.write
    written = []

    def drop_the_last_frame(self, atoms=None, **kwargs):
        written.append(1)
        if len(written) < 3:
            original(self, atoms, **kwargs)

    monkeypatch.setattr(ase_trajectory.TrajectoryWriter, "write", drop_the_last_frame)
    job = mock_job(write_structure(tmp_path, rattled_nio_structure), ensemble="nve",
                   steps=10, trajectory_interval=5)
    with pytest.raises(ResultError, match="holds 2 frames"):
        jobs.run_smoke_md(job, output_dir=tmp_path / "md")
    from nio_md_prep.mlip.provenance import read_manifest

    assert read_manifest(tmp_path / "md")["status"] == "failed"


# --- seeds -------------------------------------------------------------------


def test_velocity_and_thermostat_noise_use_independent_streams(tmp_path, rattled_nio_structure):
    """The same seed must not feed the same normals to the velocities and the noise."""
    seed = 1234
    atoms = rattled_nio_structure
    simulation = SimulationSpec(task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=1.0,
                                steps=1, seed=seed, thermostat="langevin", log_interval=1)
    first = ase_engine.run_md(atoms, mock_calculator(), simulation, workdir=tmp_path / "a")
    velocity_stream, noise_stream = np.random.SeedSequence(seed).spawn(2)
    velocity_normals = np.random.default_rng(velocity_stream).standard_normal((len(atoms), 3))
    noise_normals = np.random.default_rng(noise_stream).standard_normal((len(atoms), 3))
    assert not np.allclose(velocity_normals, noise_normals)
    # And the recorded seed reproduces the run bit for bit.
    again = ase_engine.run_md(atoms, mock_calculator(), simulation, workdir=tmp_path / "b")
    np.testing.assert_array_equal(first["final_positions"], again["final_positions"])
    assert first["integrator_resolved"]["seed"] == seed


def test_an_unset_seed_is_drawn_recorded_and_reproducible(tmp_path, rattled_nio_structure):
    simulation = SimulationSpec(task="md", ensemble="nvt", temperature_K=300.0, timestep_fs=1.0,
                                steps=5, seed=None, log_interval=1)
    first = ase_engine.run_md(rattled_nio_structure, mock_calculator(), simulation,
                              workdir=tmp_path / "a")
    resolved = first["integrator_resolved"]
    assert resolved["seed_source"] == "drawn" and 1 <= resolved["seed"] <= MAX_SEED
    replay = ase_engine.run_md(rattled_nio_structure, mock_calculator(),
                               replace(simulation, seed=resolved["seed"]), workdir=tmp_path / "b")
    np.testing.assert_array_equal(first["final_positions"], replay["final_positions"])


# --- NPT ---------------------------------------------------------------------


def hcp_mg(upper_triangular: bool = False):
    """hcp Mg; ASE's standard cell is lower-triangular, the rotated one upper."""
    from ase.build import bulk

    atoms = bulk("Mg", "hcp", a=3.21, c=5.21).repeat((3, 3, 2))
    if upper_triangular:
        a, c = 3.21 * 3, 5.21 * 2
        atoms.set_cell([[a * math.sqrt(3) / 2, -a / 2, 0.0], [0.0, a, 0.0], [0.0, 0.0, c]],
                       scale_atoms=True)
    return atoms


@pytest.mark.parametrize("upper", [False, True])
def test_mtk_npt_accepts_lower_and_upper_triangular_cells(tmp_path, upper):
    atoms = hcp_mg(upper)
    cell = np.array(atoms.get_cell())
    assert np.all(np.triu(cell, 1) == 0) != upper and np.all(np.tril(cell, -1) == 0) == upper
    job = mock_job(write_structure(tmp_path, atoms), elements=("Mg",), sigma=2.86,
                   ensemble="npt", temperature_K=300.0, pressure_bar=1.0, steps=20)
    report = jobs.run_smoke_md(job, output_dir=tmp_path / "md")
    trajectory = report["trajectory"]
    resolved = trajectory.integrator_resolved
    assert resolved["integrator"] == "IsotropicMTKNPT"
    assert set(resolved["defaults_applied"]) == {"barostat", "thermostat", "barostat_coupling"}
    assert resolved["pdamp_fs"] == 1000.0 and resolved["tdamp_fs"] == 100.0
    assert resolved["pressure_eV_per_A3"] == pytest.approx(1.0 / 1.602176634e6, rel=1e-9)
    assert trajectory.diagnostics["conserved_quantity"] == "IsotropicMTKNPT.get_conserved_energy"
    # The cell evolved isotropically, and the returned cell is the evolved one.
    final = np.array(trajectory.final_cell_angstrom)
    ratio = final[np.abs(cell) > 1e-9] / cell[np.abs(cell) > 1e-9]
    assert not np.allclose(ratio, 1.0, atol=1e-12)
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-10)
    assert trajectory.temperature_ndof == 3 * len(atoms)  # MTK thermostats 3N


def test_rotating_the_cell_does_not_change_the_energy(tmp_path):
    lower = ase_engine.singlepoint(hcp_mg(False), MockLennardJones(
        epsilon_eV=EPSILON, sigma_angstrom=2.86, cutoff_angstrom=CUTOFF))
    upper = ase_engine.singlepoint(hcp_mg(True), MockLennardJones(
        epsilon_eV=EPSILON, sigma_angstrom=2.86, cutoff_angstrom=CUTOFF))
    assert upper["energy_eV"] == pytest.approx(lower["energy_eV"], rel=1e-10)


def test_anisotropic_mtk_uses_masked_mtk(tmp_path):
    job = mock_job(write_structure(tmp_path, hcp_mg()), elements=("Mg",), sigma=2.86,
                   ensemble="npt", temperature_K=300.0, pressure_bar=1.0, steps=10,
                   barostat="mtk", barostat_coupling="anisotropic")
    resolved = jobs.run_smoke_md(job, output_dir=tmp_path / "md")["trajectory"].integrator_resolved
    assert resolved["integrator"] == "MaskedMTKNPT" and resolved["mask"] == [True, True, True]


def test_berendsen_npt_records_its_bulk_modulus(tmp_path):
    path = write_structure(tmp_path, hcp_mg())
    job = mock_job(path, elements=("Mg",), sigma=2.86, ensemble="npt", temperature_K=300.0,
                   pressure_bar=1.0, steps=10, barostat="berendsen")
    resolved = jobs.run_smoke_md(job, output_dir=tmp_path / "a")["trajectory"].integrator_resolved
    assert resolved["integrator"] == "NPTBerendsen" and resolved["thermostat"] == "berendsen"
    assert resolved["bulk_modulus_GPa"] == 100.0 and resolved["compressibility_per_GPa"] == 0.01
    assert "diagnostic default" in resolved["bulk_modulus_source"]
    # 1 GPa = 1/160.21766 eV/A^3, so B = 200 GPa is a compressibility of 0.80109 A^3/eV.
    stiff = replace(job, engine=EngineSpec(kind="ase", options={"bulk_modulus_GPa": 200.0}))
    resolved = jobs.run_smoke_md(stiff, output_dir=tmp_path / "b")["trajectory"].integrator_resolved
    assert resolved["compressibility_A3_per_eV"] == pytest.approx(160.21766208 / 200.0, rel=1e-8)
    with pytest.raises(ConfigError, match="bulk_modulus_GPa"):
        jobs.validate_job(replace(job, engine=EngineSpec(kind="ase",
                                                         options={"bulk_modulus_GPa": -1})))


def test_anisotropic_berendsen_needs_an_axis_aligned_cell(tmp_path):
    from ase.build import bulk

    job_hcp = mock_job(write_structure(tmp_path, hcp_mg(), "hcp.xyz"), elements=("Mg",),
                       sigma=2.86, ensemble="npt", temperature_K=300.0, pressure_bar=1.0,
                       barostat="berendsen", barostat_coupling="anisotropic")
    with pytest.raises(ConfigError, match="orthorhombic"):
        jobs.validate_job(job_hcp)
    cubic = bulk("Cu", "fcc", a=3.61, cubic=True).repeat((2, 2, 2))
    job_cubic = replace(
        job_hcp, structure=StructureSpec(path=write_structure(tmp_path, cubic, "cu.xyz")),
        potential=replace(job_hcp.potential, declared_elements=("Cu",), sigma_angstrom=2.27),
    )
    resolved = jobs.run_smoke_md(job_cubic, output_dir=tmp_path / "md")[
        "trajectory"].integrator_resolved
    assert resolved["integrator"] == "Inhomogeneous_NPTBerendsen" and resolved["mask"] == [1, 1, 1]


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"barostat": "parrinello-rahman"}, "MelchionnaNPT"),
        ({"thermostat": "langevin"}, "own nose-hoover thermostat"),
        ({"barostat": "berendsen", "thermostat": "nose-hoover"}, "own berendsen thermostat"),
        ({"thermostat": "csvr"}, "cannot be honoured"),
    ],
)
def test_npt_requests_the_ase_engine_cannot_honour_are_refused(tmp_path, overrides, message):
    job = mock_job(write_structure(tmp_path, hcp_mg()), elements=("Mg",), sigma=2.86,
                   ensemble="npt", temperature_K=300.0, pressure_bar=1.0, **overrides)
    with pytest.raises(ConfigError, match=message):
        jobs.validate_job(job)


def test_in_plane_npt_on_a_slab_is_refused_on_ase(tmp_path):
    slab = hcp_mg()
    slab.pbc = (True, True, False)
    job = mock_job(write_structure(tmp_path, slab), elements=("Mg",), sigma=2.86,
                   ensemble="npt", temperature_K=300.0, pressure_bar=1.0,
                   barostat_coupling="in-plane")
    with pytest.raises(ConfigError, match="only on the LAMMPS engine"):
        jobs.validate_job(job)
    isotropic = replace(job, simulation=replace(job.simulation, barostat_coupling=None))
    with pytest.raises(ConfigError, match="periodic along all three axes"):
        jobs.validate_job(isotropic)


# --- route options --------------------------------------------------------------


def test_threads_are_refused_on_the_mock_route(tmp_path, rattled_nio_structure):
    job = mock_job(write_structure(tmp_path, rattled_nio_structure), ensemble="nve",
                   engine=EngineSpec(kind="ase", threads=4))
    with pytest.raises(ConfigError, match="engine.threads"):
        jobs.validate_job(job)


def test_the_manifest_records_the_planned_integrator(tmp_path, rattled_nio_structure):
    job = mock_job(write_structure(tmp_path, rattled_nio_structure), ensemble="nvt",
                   temperature_K=300.0, thermostat="csvr")
    parameters = jobs.validate_job(job)["engine_parameters"]
    assert parameters["integrator"]["class"] == "ase.md.bussi.Bussi"
    assert parameters["integrator"]["thermostat_damping_fs"] == 100.0


# --- the execution plan -------------------------------------------------------------

PLAN_KEYS = {
    "potential_kind", "engine", "bridge", "implementation", "model_checkpoint",
    "exported_model", "elements", "energy_convention", "units", "device", "precision",
    "dynamics", "lammps", "openmm", "availability", "unmet_capabilities", "model_hashes",
}
DYNAMICS_KEYS = {
    "ensemble", "integrator", "thermostat", "barostat", "barostat_coupling", "timestep_fs",
    "timestep_native", "thermostat_damping_fs", "thermostat_damping_native",
    "barostat_damping_fs", "barostat_damping_native",
}


def test_the_mock_execution_plan_is_what_the_run_constructs(tmp_path, pinned_nio):
    """The plan's native numbers are the ones ASE's integrator reports it was built with."""
    from ase import units

    structure = write_structure(tmp_path, pinned_nio)
    job = mock_job(structure, ensemble="nvt", temperature_K=300.0, thermostat="langevin",
                   thermostat_damping_fs=40.0, timestep_fs=2.0, steps=4)
    bridge = jobs.build_bridge(job.potential, job.engine)
    plan = bridge.execution_plan(job.simulation, pinned_nio)
    assert PLAN_KEYS <= set(plan) and DYNAMICS_KEYS <= set(plan["dynamics"])
    assert plan["potential_kind"] == "mock" and plan["engine"] == "ase"
    # The shared contract: a model-free route still has a model_checkpoint section.
    assert plan["model_checkpoint"] == {"path": None, "sha256": None}
    assert plan["exported_model"] is None
    assert plan["lammps"] is None and plan["openmm"] is None
    assert plan["units"] == {**plan["units"], "native": "ase", "pressure_unit": "eV/Angstrom^3"}
    assert plan["device"]["effective"] == "cpu" and plan["device"]["guaranteed"] is True
    assert plan["precision"]["effective"] == "float64" and plan["precision"]["guaranteed"] is True
    assert plan["availability"]["available"] and plan["unmet_capabilities"] == []
    dynamics = plan["dynamics"]
    assert dynamics["integrator"] == "ase.md.langevin.Langevin"
    assert dynamics["timestep_native"] == pytest.approx(2.0 * units.fs, rel=1e-15)
    assert dynamics["thermostat_damping_native"] == pytest.approx(40.0 * units.fs, rel=1e-15)
    assert dynamics["barostat"] is None and dynamics["barostat_damping_native"] is None
    assert dynamics["temperature_ndof"] == 3 * (len(pinned_nio) - 4)

    trajectory = jobs.run_smoke_md(job, output_dir=tmp_path / "md")["trajectory"]
    todict = trajectory.integrator_resolved["ase_todict"]
    assert todict["timestep"] == pytest.approx(dynamics["timestep_native"], rel=1e-15)
    assert todict["friction"] == pytest.approx(dynamics["friction_native"], rel=1e-15)
    assert trajectory.integrator_resolved["class"] == dynamics["integrator"]
    assert trajectory.temperature_ndof == dynamics["temperature_ndof"]


def test_the_mock_npt_plan_gives_the_barostat_in_native_units(tmp_path):
    from ase import units

    atoms = hcp_mg()
    job = mock_job(write_structure(tmp_path, atoms), elements=("Mg",), sigma=2.86,
                   ensemble="npt", temperature_K=300.0, pressure_bar=1000.0,
                   barostat="berendsen", barostat_damping_fs=500.0)
    plan = jobs.build_bridge(job.potential, job.engine).execution_plan(job.simulation, atoms)
    dynamics = plan["dynamics"]
    assert dynamics["integrator"] == "ase.md.nptberendsen.NPTBerendsen"
    assert dynamics["thermostat"] == "berendsen" and dynamics["barostat"] == "berendsen"
    assert dynamics["barostat_damping_native"] == pytest.approx(500.0 * units.fs, rel=1e-15)
    # 1000 bar in ASE's native eV/Angstrom^3 -- never atm -- exactly as passed to ASE.
    from nio_md_prep.mlip.units import EV_PER_ANGSTROM3_IN_BAR

    assert dynamics["pressure_native"] == 1000.0 / EV_PER_ANGSTROM3_IN_BAR
    assert dynamics["constructor_kwargs_native"]["pressure_au"] == dynamics["pressure_native"]
    # ase.units defaults to CODATA 2014 (e = 1.6021766208e-19 C, vs the exact SI 2019
    # value used here): 8e-9 relative, far inside any statistical NPT tolerance.
    assert dynamics["pressure_native"] == pytest.approx(1000.0 * units.bar, rel=1e-7)
    assert dynamics["pressure_native_unit"] == "eV/Angstrom^3"


def test_a_refused_md_request_is_listed_in_the_plan_not_raised(tmp_path):
    atoms = hcp_mg()
    job = mock_job(write_structure(tmp_path, atoms), elements=("Mg",), sigma=2.86,
                   ensemble="npt", temperature_K=300.0, pressure_bar=1.0,
                   barostat="parrinello-rahman")
    with pytest.raises(ConfigError, match="MelchionnaNPT"):
        jobs.validate_job(job)
    plan = jobs.build_bridge(job.potential, job.engine).execution_plan(job.simulation, atoms)
    assert any("MelchionnaNPT" in item for item in plan["unmet_capabilities"])
    assert plan["unmet_capabilities"]
    dynamics = plan["dynamics"]
    assert DYNAMICS_KEYS <= set(dynamics)
    assert dynamics["integrator"] is None and dynamics["timestep_native"] is None
    assert dynamics["barostat"] == "parrinello-rahman" and "MelchionnaNPT" in dynamics["refused"]


def test_ase_npt_is_refused_with_frozen_atoms_and_on_vacuum_slabs():
    """NPT must never scale frozen atoms or vacuum silently on the ASE engine."""
    from ase.build import bulk, fcc111
    from ase.constraints import FixAtoms

    from nio_md_prep.mlip.engines import ase_engine
    from nio_md_prep.mlip.errors import ConfigError
    from nio_md_prep.mlip.specs import SimulationSpec

    npt = SimulationSpec(task="md", ensemble="npt", steps=10, timestep_fs=1.0,
                         temperature_K=300.0, pressure_bar=1.0, barostat="mtk")

    frozen = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((2, 2, 2))
    frozen.set_constraint(FixAtoms(indices=[0, 1]))
    with pytest.raises(ConfigError, match="frozen atom"):
        ase_engine.check_simulation(npt, frozen)

    slab = fcc111("Cu", size=(2, 2, 3), vacuum=8.0)
    slab.pbc = (True, True, True)  # a POSCAR-style slab: fully periodic with vacuum
    with pytest.raises(ConfigError, match="vacuum gap"):
        ase_engine.check_simulation(npt, slab)

    bulk_cell = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((2, 2, 2))
    assert ase_engine.check_simulation(npt, bulk_cell).integrator == "IsotropicMTKNPT"
