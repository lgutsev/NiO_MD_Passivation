"""Tests for the VASP parsers: dataset/vasprun.py and dataset/vaspfiles.py.

Two kinds of input are used:

* SYNTHETIC runs written by ``tests/dataset_fixtures.py`` (labelled synthetic
  in every file; not real VASP output) for edge cases no local real file has:
  VASP <= 6.0.8 calc-level energy mislabel, PSTRESS, Edisp, missing blocks,
  unconverged SCF, MLFF flat steps, truncation, compression, selective
  dynamics, spin polarisation with/without MAGMOM, species-order mismatches.
* REAL VASP output shipped with ASE (``ase/test/testdata/vasp``), read in
  place and never copied (amendments item 15); those tests skip when ASE or
  the files are absent. They validate the parsers on real formats only -- no
  real NiO force export can be validated on this machine (no NiO OUTCAR or
  vasprun.xml exists locally).

Where ASE is available it is used as an INDEPENDENT reader: our energies,
forces and geometry must equal what ``ase.io.read`` returns for the same file.
"""

from __future__ import annotations

import bz2
import gzip
import hashlib
import math
import tracemalloc
from pathlib import Path

import numpy as np
import pytest

import dataset_fixtures as fx
from nio_md_prep.dataset import PARSER_NAME, PARSER_VERSION
from nio_md_prep.dataset import vaspfiles as vf
from nio_md_prep.dataset import vasprun as vr
from nio_md_prep.dataset.model import IonicStep, StepStructure, vasp_stress_to_ase


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _nio(n_frames=3, **frame_fields):
    species, cell, positions = fx.rocksalt_nio()
    frames = fx.make_frames(species, cell, positions, n_frames, **frame_fields)
    return species, cell, frames


def _run(tmp_path, name="run", n_frames=3, frame_fields=None, **kwargs):
    species, cell, frames = _nio(n_frames, **(frame_fields or {}))
    return fx.write_vasp_run(tmp_path / name, species=species, cell=cell, frames=frames, **kwargs)


def _read(run):
    return vr.read_vasprun(run["vasprun"])


def _problem_codes(problems):
    return [problem.split(":", 1)[0] for problem in problems]


def _ase_testdata() -> Path:
    ase = pytest.importorskip("ase")
    directory = Path(ase.__file__).parent / "test" / "testdata" / "vasp"
    if not directory.is_dir():
        pytest.skip("ASE test data directory not shipped with this ASE install")
    return directory


def _ase_file(name: str) -> Path:
    path = _ase_testdata() / name
    if not path.is_file():
        pytest.skip(f"ASE test data file {name} not present")
    return path


def _synthetic_step(energies, scf_energies, *, label_source="dft", cell=None):
    """A SYNTHETIC IonicStep (not VASP output) for unit tests of derived_energies."""
    cell = np.eye(3) * 3.0 if cell is None else np.asarray(cell, dtype=float)
    fractional = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    structure = StepStructure(cell=cell, fractional=fractional, positions=fractional @ cell)
    return IonicStep(index=0, complete=True, structure=structure, forces=np.zeros((2, 3)), stress_kbar_vasp=None,
                     energies=dict(energies), scf_energies=[dict(item) for item in scf_energies],
                     label_source=label_source)


# --------------------------------------------------------------------------
# synthetic fixtures are labelled synthetic in-file
# --------------------------------------------------------------------------

def test_fixture_files_are_labelled_synthetic(tmp_path):
    run = _run(tmp_path, n_frames=2)
    text = run["vasprun"].read_text(encoding="latin-1")
    assert fx.SYNTHETIC_COMMENT in text and fx.SYNTHETIC_SUBVERSION in text
    for key in ("outcar", "oszicar", "incar", "poscar", "kpoints", "potcar"):
        assert "SYNTHETIC" in run[key].read_text(encoding="latin-1"), key
    header, _steps, _trailer = _read(run)
    assert header.generator["subversion"] == fx.SYNTHETIC_SUBVERSION


# --------------------------------------------------------------------------
# vasprun: header, steps, geometry/label pairing
# --------------------------------------------------------------------------

def test_modern_md_run_header_and_steps(tmp_path):
    run = _run(tmp_path, n_frames=4, incar={"IBRION": 0, "POTIM": 2.0, "ISIF": 2})
    header, steps, trailer = _read(run)
    expected = run["expected"]
    assert header.generator["program"] == "vasp" and header.version_tuple == (6, 4, 2)
    assert header.species == expected["species"]
    assert [(t.element, t.count) for t in header.atom_types] == [("Ni", 4), ("O", 4)]
    assert [t.pseudopotential for t in header.atom_types] == expected["titles"]
    assert header.parameters["IBRION"] == 0 and isinstance(header.parameters["IBRION"], int)
    assert header.parameters["NELM"] == 60 and header.parameters["EDIFF"] == pytest.approx(1e-5)
    assert header.parameters["LNONCOLLINEAR"] is False
    assert header.parameters["MAGMOM"] == [1.0] * 8  # VASP writes the default even for ISPIN=1
    assert header.parameter_duplicates == ["OMEGAMAX"]
    assert header.problems == []
    assert header.kpoints["scheme"] == "Gamma" and header.kpoints["nkpts"] == 1
    assert trailer.closed and not trailer.truncated and trailer.final_structure_present
    assert trailer.steps_seen == trailer.n_complete_steps == 4
    assert [step.index for step in steps] == [0, 1, 2, 3]
    for step, exp in zip(steps, expected["steps"]):
        assert step.complete and step.label_source == "dft" and step.problems == []
        # geometry and labels come from the same <calculation>
        np.testing.assert_array_equal(step.structure.cell, exp["cell"])
        np.testing.assert_array_equal(step.structure.fractional, exp["fractional"])
        np.testing.assert_allclose(step.structure.positions, exp["positions"], rtol=0, atol=1e-12)
        np.testing.assert_array_equal(step.forces, exp["forces"])
        np.testing.assert_array_equal(step.stress_kbar_vasp, exp["stress_kbar"])
        assert step.energies["e_fr_energy"] == exp["calc_energies"]["e_fr_energy"]
        assert step.energies["kinetic"] == pytest.approx(0.1)  # MD block kept raw, never a label
        assert len(step.scf_energies) == exp["n_scf"] == 12
        assert step.structure.selective is None  # flags never live inside <calculation>
    # finalpos is exposed but is NOT the last force-evaluated geometry
    assert not np.array_equal(trailer.final_structure.fractional, steps[-1].structure.fractional)
    # sha256 of the stored bytes
    assert trailer.source_sha256 == hashlib.sha256(run["vasprun"].read_bytes()).hexdigest()
    assert trailer.source_bytes == run["vasprun"].stat().st_size


def test_iter_vasprun_is_streaming_and_single_use(tmp_path):
    run = _run(tmp_path, n_frames=2)
    header, steps, trailer = vr.iter_vasprun(run["vasprun"])
    with pytest.raises(RuntimeError):
        trailer()  # not before the steps are exhausted
    assert len(list(steps)) == 2
    assert trailer().closed
    with vr.VasprunReader(run["vasprun"]) as reader:
        list(reader.steps())
        with pytest.raises(RuntimeError):
            reader.steps()


def test_parsing_is_deterministic(tmp_path):
    run = _run(tmp_path, n_frames=3)
    first = _read(run)
    second = _read(run)
    assert first[2].source_sha256 == second[2].source_sha256
    for a, b in zip(first[1], second[1]):
        assert a.energies == b.energies and a.scf_energies == b.scf_energies
        np.testing.assert_array_equal(a.forces, b.forces)


def _streaming_peak(path):
    tracemalloc.start()
    try:
        _header, steps, trailer = vr.iter_vasprun(path)
        count = sum(1 for _step in steps)
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert trailer().closed
    return count, peak


def test_streaming_memory_does_not_grow_with_steps(tmp_path):
    small = _run(tmp_path, name="small", n_frames=50)
    large = _run(tmp_path, name="large", n_frames=800)
    n_small, peak_small = _streaming_peak(small["vasprun"])
    n_large, peak_large = _streaming_peak(large["vasprun"])
    size_large = large["vasprun"].stat().st_size
    assert (n_small, n_large) == (50, 800)
    # a whole-tree parse needs ~10x the file size; streaming stays flat in the step count
    assert peak_large < 1.5 * peak_small + 200_000, (peak_small, peak_large)
    assert peak_large < size_large / 5, (peak_large, size_large)


# --------------------------------------------------------------------------
# energies: reconstruction, version gating, provenance
# --------------------------------------------------------------------------

def test_modern_energies_are_calc_level_direct_after_cross_check(tmp_path):
    run = _run(tmp_path, n_frames=3, version="6.4.2", incar={"IBRION": 2})
    header, steps, _ = _read(run)
    for step, exp in zip(steps, run["expected"]["steps"]):
        d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
        assert d["free_energy"] == exp["free_energy"]
        assert d["energy_sigma0"] == pytest.approx(exp["energy_sigma0"], abs=1e-9)
        assert d["energy_no_entropy"] == pytest.approx(exp["energy_no_entropy"], abs=1e-9)
        assert d["energy_rule"] == vr.RULE_CALC_LEVEL_DIRECT == "calc_level_direct"
        assert d["energy_source"] == "vasprun:calculation.e_fr_energy-PSTRESS*V"
        assert d["energy_sources"]["energy_sigma0"] == "vasprun:calculation.e_0_energy-PSTRESS*V"
        assert d["energy_rules"] == {name: "calc_level_direct" for name in vr.ENERGY_QUANTITIES}
        assert d["energy_version_gate"] == "vasp>=6.1.1" and d["vasp_version"] == "6.4.2"
        assert d["calc_level_pattern"] == "consistent" and d["energy_flags"] == []
        assert all(abs(delta) <= 1e-8 for delta in d["calc_level_deltas"].values())
        assert d["parser"] == PARSER_NAME and d["parser_version"] == PARSER_VERSION
        assert d["smearing_entropy_term"] == pytest.approx(0.004, abs=1e-9)
        assert d["pv_term"] == 0.0 and abs(d["additive_correction"]) < 1e-9


@pytest.mark.parametrize("version", ["5.2.12", "5.4.4", "6.0.8"])
def test_legacy_mislabelled_calc_energies_are_reconstructed(tmp_path, version):
    run = _run(tmp_path, n_frames=2, version=version, incar={"IBRION": 2})
    assert run["expected"]["calc_energy_bug"]
    header, steps, _ = _read(run)
    for step, exp in zip(steps, run["expected"]["steps"]):
        d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
        # F is the calculation-level e_fr_energy (never mislabelled)
        assert d["free_energy"] == exp["free_energy"]
        assert d["energy_rules"]["free_energy"] == "calc_level_direct"
        # E0 / E_wo reconstructed from the last scstep, not read from the mislabelled calc level
        assert d["energy_sigma0"] == pytest.approx(exp["energy_sigma0"], abs=1e-9)
        assert d["energy_no_entropy"] == pytest.approx(exp["energy_no_entropy"], abs=1e-9)
        naive_e0 = step.energies["e_0_energy"]
        assert abs(naive_e0 - exp["energy_sigma0"]) > 10.0  # holds F - E_wo, not E0
        assert d["energy_rule"] == "reconstructed_last_scstep:vasp<=6.0.8"
        assert d["energy_rules"]["energy_sigma0"] == d["energy_rule"]
        assert d["energy_sources"]["energy_sigma0"].startswith("vasprun:(calculation.e_fr_energy-PSTRESS*V)+")
        assert d["calc_level_pattern"] == "mislabelled"
        assert "calc_level_mislabel_pattern" in d["energy_flags"]
        assert "calc_level_energy_mismatch" not in d["energy_flags"]


@pytest.mark.parametrize(
    ("version", "gate", "rule"),
    [
        ("6.0.8", "vasp<=6.0.8", "reconstructed_last_scstep:vasp<=6.0.8"),
        ("6.1.0", "version_unverified", "reconstructed_last_scstep:version_unverified"),
        ("6.1.1", "vasp>=6.1.1", "calc_level_direct"),
        ("6.5.0", "vasp>=6.1.1", "calc_level_direct"),
    ],
)
def test_energy_version_gate_boundaries(tmp_path, version, gate, rule):
    run = _run(tmp_path, n_frames=1, version=version, calc_energy_bug=False)
    header, steps, _ = _read(run)
    assert vr.energy_version_gate(header.version_tuple) == gate
    d = vr.derived_energies(steps[0], header.parameters["PSTRESS"], header.version_tuple)
    assert d["energy_rule"] == rule
    assert d["energy_sigma0"] == pytest.approx(run["expected"]["steps"][0]["energy_sigma0"], abs=1e-9)


def test_unknown_version_is_reconstructed(tmp_path):
    run = _run(tmp_path, n_frames=1, version="6.4.2")
    text = run["vasprun"].read_text(encoding="latin-1").replace(">6.4.2  <", ">unknown-build  <")
    run["vasprun"].write_text(text, encoding="latin-1")
    header, steps, _ = _read(run)
    assert header.version_tuple is None
    d = vr.derived_energies(steps[0], header.parameters["PSTRESS"], header.version_tuple)
    assert d["energy_version_gate"] == "version_unknown"
    assert d["energy_rule"] == "reconstructed_last_scstep:version_unknown"
    # the default (no version passed) is the fail-safe reconstruction too
    assert vr.derived_energies(steps[0], 0.0)["energy_rule"] == "reconstructed_last_scstep:version_unknown"


def test_modern_version_with_disagreeing_calc_level_uses_reconstruction_and_flags(tmp_path):
    run = _run(tmp_path, n_frames=1, version="6.4.2", calc_energy_bug=True)
    header, steps, _ = _read(run)
    d = vr.derived_energies(steps[0], header.parameters["PSTRESS"], header.version_tuple)
    assert d["energy_rule"] == "reconstructed_last_scstep:calc_level_mismatch"
    assert "calc_level_energy_mismatch" in d["energy_flags"]
    assert d["energy_sigma0"] == pytest.approx(run["expected"]["steps"][0]["energy_sigma0"], abs=1e-9)


def test_cross_check_threshold_is_one_micro_ev():
    """SYNTHETIC IonicSteps: calc-level e_0 off by 5e-7 eV is accepted, 5e-6 eV is not."""
    last = {"e_fr_energy": -10.0, "e_wo_entrp": -10.004, "e_0_energy": -10.002}
    for offset, rule in ((5e-7, "calc_level_direct"), (5e-6, "reconstructed_last_scstep:calc_level_mismatch")):
        calc = {"e_fr_energy": -10.0, "e_wo_entrp": -10.004, "e_0_energy": -10.002 + offset}
        d = vr.derived_energies(_synthetic_step(calc, [last, last]), 0.0, (6, 3, 0))
        assert d["energy_rule"] == rule
        assert d["calc_level_deltas"]["energy_sigma0"] == pytest.approx(offset, abs=1e-12)
        expected = -10.002 + offset if rule == "calc_level_direct" else -10.002
        assert d["energy_sigma0"] == pytest.approx(expected, abs=1e-12)


def test_energy_edge_cases_are_explicit():
    """SYNTHETIC IonicSteps for missing inputs; nothing is guessed."""
    last = {"e_fr_energy": -10.0, "e_wo_entrp": -10.004, "e_0_energy": -10.002}
    calc = {"e_fr_energy": -10.0, "e_wo_entrp": -10.004, "e_0_energy": -10.002}
    d = vr.derived_energies(_synthetic_step(calc, [last]), None, (6, 3, 0))
    assert d["free_energy"] is None and "pstress_unknown" in d["energy_flags"] and d["energy_rule"] == "unavailable"
    d = vr.derived_energies(_synthetic_step({}, [last]), 0.0, (6, 3, 0))
    assert d["free_energy"] is None and "missing_e_fr_energy" in d["energy_flags"]
    d = vr.derived_energies(_synthetic_step(calc, []), 0.0, (6, 3, 0))
    assert d["free_energy"] == -10.0 and d["energy_sigma0"] is None
    assert d["energy_rule"] == "free_energy_only:no_scstep_energies" and "scstep_energies_missing" in d["energy_flags"]
    d = vr.derived_energies(_synthetic_step({"e_fr_energy": -10.0}, [last]), 0.0, (6, 3, 0))
    assert d["energy_rule"] == "reconstructed_last_scstep:calc_level_missing"
    assert d["energy_sigma0"] == pytest.approx(-10.002)
    d = vr.derived_energies(_synthetic_step({"e_fr_energy": math.nan}, [last]), 0.0, (6, 3, 0))
    assert "non_finite_energy" in d["energy_flags"]
    value, source, rule = vr.label_energy(vr.derived_energies(_synthetic_step(calc, [last]), 0.0, (6, 3, 0)))
    assert (value, source, rule) == (-10.0, "vasprun:calculation.e_fr_energy-PSTRESS*V", "calc_level_direct")
    with pytest.raises(ValueError):
        vr.label_energy({}, "energy")


def test_pstress_pv_term_is_removed_with_the_vasp_constant(tmp_path):
    run = _run(tmp_path, n_frames=2, pstress=5.0, incar={"IBRION": 2})
    header, steps, _ = _read(run)
    assert header.parameters["PSTRESS"] == 5.0
    for step, exp in zip(steps, run["expected"]["steps"]):
        d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
        volume = abs(np.linalg.det(step.structure.cell))
        assert d["pv_term"] == 5.0 * vr.PV_EV_PER_KBAR_A3 * volume
        assert d["free_energy"] == step.energies["e_fr_energy"] - d["pv_term"] == exp["free_energy"]
        assert d["free_energy"] == pytest.approx(exp["free_energy_input"], abs=1e-8)
        assert abs(d["additive_correction"]) < 1e-7  # only the 8-decimal print rounding remains
        assert d["energy_rule"] == "calc_level_direct"


def test_dispersion_offset_is_an_additive_correction(tmp_path):
    run = _run(tmp_path, n_frames=2, frame_fields={"edisp": -0.58552}, incar={"IBRION": 2, "IVDW": 12})
    header, steps, _ = _read(run)
    assert "IVDW" not in header.parameters and header.incar["IVDW"] == 12  # as in real VASP
    for step, exp in zip(steps, run["expected"]["steps"]):
        d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
        assert d["free_energy"] == exp["free_energy"]
        assert d["additive_correction"] == pytest.approx(-0.58552, abs=2e-8)
        assert d["energy_sigma0"] == pytest.approx(exp["free_energy"] - 0.002, abs=2e-8)
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.dispersion_energies == pytest.approx([-0.58552, -0.58552])
    assert outcar.executed_tags["IVDW"] == "12"


def test_energies_match_ase_as_independent_reader(tmp_path):
    ase_io = pytest.importorskip("ase.io")
    for version, pstress in (("6.4.2", 0.0), ("5.4.4", 3.0)):
        run = _run(tmp_path, name=f"ase_{version}", n_frames=3, version=version, pstress=pstress, incar={"IBRION": 2})
        header, steps, _ = _read(run)
        images = ase_io.read(run["vasprun"], index=":", format="vasp-xml")
        assert len(images) == len(steps)
        for atoms, step in zip(images, steps):
            d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
            assert atoms.get_potential_energy(force_consistent=True) == d["free_energy"]
            assert atoms.get_potential_energy() == d["energy_sigma0"]
            np.testing.assert_array_equal(atoms.get_forces(apply_constraint=False), step.forces)
            np.testing.assert_allclose(atoms.get_positions(), step.structure.positions, rtol=0, atol=1e-12)
            np.testing.assert_array_equal(np.asarray(atoms.get_cell()), step.structure.cell)
            np.testing.assert_allclose(atoms.get_stress(voigt=False), vasp_stress_to_ase(step.stress_kbar_vasp),
                                       rtol=1e-6, atol=1e-12)


# --------------------------------------------------------------------------
# missing labels, stress, non-finite values
# --------------------------------------------------------------------------

def test_missing_forces_block_is_reported(tmp_path):
    species, cell, frames = _nio(3)
    frames[1]["omit"] = {"forces"}
    run = fx.write_vasp_run(tmp_path / "noforces", species=species, cell=cell, frames=frames, incar={"IBRION": 2})
    _header, steps, trailer = _read(run)
    assert trailer.closed and len(steps) == 3
    assert steps[1].complete and steps[1].forces is None
    assert "missing_forces" in _problem_codes(steps[1].problems)
    assert steps[0].problems == [] and steps[2].problems == []
    outcar = vf.parse_outcar(run["outcar"])
    assert math.isnan(outcar.force_abs_sums[1])
    assert any(p.startswith("missing_forces: OUTCAR DFT step 1") for p in outcar.problems)


def test_missing_energy_block_is_reported(tmp_path):
    species, cell, frames = _nio(2)
    frames[0]["omit"] = {"energy"}
    run = fx.write_vasp_run(tmp_path / "noenergy", species=species, cell=cell, frames=frames, incar={"IBRION": 2})
    header, steps, _ = _read(run)
    assert "missing_energy" in _problem_codes(steps[0].problems)
    d = vr.derived_energies(steps[0], header.parameters["PSTRESS"], header.version_tuple)
    assert d["free_energy"] is None and "missing_e_fr_energy" in d["energy_flags"]


def test_isif0_run_has_no_stress_anywhere(tmp_path):
    run = _run(tmp_path, n_frames=3, incar={"IBRION": 0, "ISIF": 0})
    header, steps, _ = _read(run)
    assert header.parameters["ISIF"] == 0
    for step in steps:
        assert step.stress_kbar_vasp is None and step.problems == []  # absent, never zeros
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.stress_kbar == [None, None, None]


def test_overflow_values_are_non_finite_problems(tmp_path):
    species, cell, frames = _nio(2)
    frames[1]["overflow"] = {"forces"}
    run = fx.write_vasp_run(tmp_path / "overflow", species=species, cell=cell, frames=frames, incar={"IBRION": 2},
                            raw_parameters={"ENAUG": (None, "   *********")})
    header, steps, _ = _read(run)
    assert header.parameters["ENAUG"] is None
    assert any(p.startswith("non_finite: parameters.ENAUG") for p in header.problems)
    assert "non_finite" in _problem_codes(steps[1].problems)
    assert np.isnan(steps[1].forces[0, 0]) and np.isfinite(steps[1].forces[1:]).all()
    assert steps[0].problems == []


def test_header_species_inconsistency_is_reported(tmp_path):
    run = _run(tmp_path, n_frames=1)
    text = run["vasprun"].read_text(encoding="latin-1").replace("<atoms>       8 </atoms>", "<atoms>       9 </atoms>")
    run["vasprun"].write_text(text, encoding="latin-1")
    header, _steps, _ = _read(run)
    assert any(p.startswith("species_mismatch: atominfo declares 9 atoms") for p in header.problems)


# --------------------------------------------------------------------------
# parse failures are never silent
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "content",
    [
        b"",
        b"this is not xml at all\n",
        b'<?xml version="1.0"?>\n<notvasp><generator/></notvasp>\n',
        b'<?xml version="1.0" encoding="ISO-8859-1"?>\n<!-- SYNTHETIC -->\n<modeling>\n <generator>\n  <i name="program"',
    ],
)
def test_unreadable_vasprun_raises_a_parse_error(tmp_path, content):
    path = tmp_path / "vasprun.xml"
    path.write_bytes(content)
    with pytest.raises(vr.VasprunParseError):
        vr.read_vasprun(path)


def test_header_without_steps_is_readable_but_truncated(tmp_path):
    run = _run(tmp_path, n_frames=2)
    text = run["vasprun"].read_text(encoding="latin-1")
    run["vasprun"].write_text(text.split(" <calculation>")[0], encoding="latin-1")
    header, steps, trailer = _read(run)
    assert header.species and steps == []
    assert trailer.truncated and not trailer.closed and trailer.steps_seen == 0


def test_garbage_evidence_files_never_raise(tmp_path):
    junk = tmp_path / "junk"
    junk.write_bytes(b"\x00\xff garbage \n more garbage\r\n")
    outcar = vf.parse_outcar(junk)
    assert outcar.free_energies == [] and outcar.problems
    oszicar = vf.parse_oszicar(junk)
    assert oszicar.ionic_lines == []
    poscar = vf.parse_poscar(junk)
    assert poscar.cell is None and poscar.problems
    potcar = vf.parse_potcar(junk)
    assert potcar.datasets == [] and potcar.problems
    incar = vf.parse_incar(junk)
    assert incar.problems


# --------------------------------------------------------------------------
# SCF convergence evidence
# --------------------------------------------------------------------------

def test_unconverged_scf_step_evidence_vasp6(tmp_path):
    species, cell, frames = _nio(3)
    frames[1].update(n_scf=60, last_dE=5e-3, scf_marker="not_reached")
    run = fx.write_vasp_run(tmp_path / "unconv", species=species, cell=cell, frames=frames, nelm=60, ediff=1e-5,
                            incar={"IBRION": 2})
    header, steps, _ = _read(run)
    summaries = [vr.scf_summary(step) for step in steps]
    assert [s["n_steps"] for s in summaries] == [12, 60, 12]
    assert summaries[1]["n_steps"] >= header.parameters["NELM"]
    assert summaries[1]["last_dE"] == pytest.approx(5e-3, rel=1e-6) and summaries[1]["last_dE"] > header.parameters["EDIFF"]
    assert summaries[0]["last_dE"] < header.parameters["EDIFF"]
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.scf_marker_kinds == ["reached", "not_reached", "reached"]
    assert outcar.scf_converged_markers == [True, False, True]
    assert outcar.scf_iteration_counts == [12, 60, 12]
    markers = vf.outcar_scf_markers(outcar, n_dft_steps=3)
    assert markers["mapped"] and markers["marker_is_proof"] and markers["kinds"][1] == "not_reached"
    oszicar = vf.parse_oszicar(run["oszicar"])
    assert oszicar.scf_steps == [12, 60, 12]


def test_vasp5_reached_marker_is_not_proof(tmp_path):
    species, cell, frames = _nio(2)
    frames[0].update(n_scf=60, last_dE=5e-3, scf_marker="reached")  # VASP 5 prints 'reached' even at NELM
    run = fx.write_vasp_run(tmp_path / "v5", species=species, cell=cell, frames=frames, version="5.4.4",
                            nelm=60, incar={"IBRION": 2})
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.version_tuple == (5, 4, 4)
    markers = vf.outcar_scf_markers(outcar, n_dft_steps=2)
    assert markers["mapped"] and not markers["marker_is_proof"]
    assert vr.scf_summary(_read(run)[1][0])["n_steps"] == 60


def test_hard_stop_marker_and_marker_count_gate(tmp_path):
    species, cell, frames = _nio(3)
    frames[0]["scf_marker"] = "hard_stop"
    frames[2]["scf_marker"] = None  # e.g. NWRITE<=1 prints the exit line for some steps only
    run = fx.write_vasp_run(tmp_path / "markers", species=species, cell=cell, frames=frames, incar={"IBRION": 2})
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.scf_marker_kinds == ["hard_stop", "reached", None]
    assert outcar.scf_converged_markers == [False, True, None]
    gated = vf.outcar_scf_markers(outcar, n_dft_steps=3)
    assert not gated["mapped"] and gated["kinds"] is None and gated["reason"].startswith("marker_count_mismatch")
    assert vf.outcar_scf_markers(outcar, n_dft_steps=4)["reason"].startswith("step_count_mismatch")


# --------------------------------------------------------------------------
# VASP-MLFF flat steps
# --------------------------------------------------------------------------

def test_mlff_flat_steps_are_yielded_but_never_dft(tmp_path):
    run = _run(tmp_path, n_frames=6, mlff_steps=(1, 2, 4),
               incar={"IBRION": 0, "POTIM": 2.0, "ISIF": 0, "ML_LMLFF": True, "ML_MODE": "train"})
    header, steps, trailer = _read(run)
    assert [s.index for s in steps] == list(range(6))  # the index counts ALL steps
    assert [s.label_source for s in steps] == ["dft", "mlff", "mlff", "dft", "mlff", "dft"]
    assert trailer.n_complete_steps == 6 and trailer.closed
    for step in steps:
        assert step.complete and step.problems == []
        if step.label_source == "mlff":
            assert step.scf_energies == []
            assert step.energies["e_fr_energy"] == step.energies["e_wo_entrp"] == step.energies["e_0_energy"]
            d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
            assert d["energy_rule"] == "mlff_flat_step:not_a_dft_label"
            assert d["energy_source"].startswith("vasprun:mlff_flat_step")
            assert d["energy_sigma0"] is None
        else:
            assert len(step.scf_energies) == 12
    # ML_* tags are only in <incar>, never in <parameters>
    assert "ML_LMLFF" not in header.parameters and header.incar["ML_LMLFF"] is True
    settings = vr.mlff_settings(header)
    assert settings["active"] is True and settings["mode"] == "train" and settings["dft_labels"] == "per_step"
    outcar = vf.parse_outcar(run["outcar"])
    assert len(outcar.free_energies) == 3  # DFT steps only
    assert outcar.ml_steps == 6  # an ML prediction block precedes every step in train mode
    oszicar = vf.parse_oszicar(run["oszicar"])
    assert len(oszicar.ionic_lines) == 6 and oszicar.scf_steps == [12, 0, 0, 12, 0, 12]
    ase_io = pytest.importorskip("ase.io")
    assert len(ase_io.read(run["vasprun"], index=":", format="vasp-xml")) == 3  # ASE drops ML steps


@pytest.mark.parametrize(
    ("incar", "dft_labels", "mode"),
    [
        ({}, "all", None),
        ({"ML_LMLFF": True, "ML_MODE": "run"}, "none", "run"),
        ({"ML_LMLFF": True, "ML_MODE": "refitbayesian"}, "none", "refitbayesian"),
        ({"ML_LMLFF": True, "ML_MODE": "delta"}, "not_pure_dft", "delta"),
        ({"ML_LMLFF": True, "ML_ISTART": 3}, "none", "select"),
        ({"ML_LMLFF": False, "ML_MODE": "train"}, "all", "train"),
    ],
)
def test_mlff_settings_modes(tmp_path, incar, dft_labels, mode):
    run = _run(tmp_path, n_frames=1, incar=dict(incar, IBRION=-1))
    header, _steps, _ = _read(run)
    settings = vr.mlff_settings(header)
    assert settings["dft_labels"] == dft_labels and settings["mode"] == mode


# --------------------------------------------------------------------------
# truncation and compression
# --------------------------------------------------------------------------

def test_truncated_between_steps(tmp_path):
    run = _run(tmp_path, n_frames=4, truncate_after_steps=2, incar={"IBRION": 0})
    _header, steps, trailer = _read(run)
    assert [s.complete for s in steps] == [True, True]
    assert trailer.truncated and not trailer.closed and not trailer.final_structure_present
    assert trailer.n_complete_steps == trailer.steps_seen == 2 and trailer.partial_tail_tag is None
    assert trailer.error and "ParseError" in trailer.error
    outcar = vf.parse_outcar(run["outcar"])
    assert not outcar.completed and len(outcar.free_energies) == 2


@pytest.mark.parametrize("mode", [True, "mid_forces"])
def test_truncated_inside_a_step(tmp_path, mode):
    run = _run(tmp_path, n_frames=4, truncate_after_steps=2, truncate_inside_step=mode, incar={"IBRION": 0})
    _header, steps, trailer = _read(run)
    assert len(steps) == 3 and [s.complete for s in steps] == [True, True, False]
    partial = steps[2]
    assert partial.index == 2 and "truncated_frame" in _problem_codes(partial.problems)
    assert partial.energies == {}  # never synthesized from initialpos or a neighbour
    if mode == "mid_forces":
        assert partial.forces is None
    else:
        np.testing.assert_array_equal(partial.forces, run["expected"]["steps"][2]["forces"])
    assert trailer.truncated and trailer.partial_tail_tag == "calculation"
    assert trailer.n_complete_steps == 2 and trailer.steps_seen == 3
    outcar = vf.parse_outcar(run["outcar"])
    assert len(outcar.free_energies) == 2
    assert any(p.startswith("truncated_frame: OUTCAR ends inside DFT ionic step 2") for p in outcar.problems)
    oszicar = vf.parse_oszicar(run["oszicar"])
    assert len(oszicar.ionic_lines) == 2 and "truncated_frame" in _problem_codes(oszicar.problems)


def test_truncated_inside_an_mlff_flat_step(tmp_path):
    run = _run(tmp_path, n_frames=4, mlff_steps=(1, 2, 3), truncate_after_steps=2, truncate_inside_step=True,
               incar={"IBRION": 0, "ISIF": 0, "ML_LMLFF": True})
    _header, steps, trailer = _read(run)
    assert [s.label_source for s in steps] == ["dft", "mlff", "mlff"]
    assert [s.complete for s in steps] == [True, True, False]
    assert "truncated_frame" in _problem_codes(steps[2].problems)
    assert trailer.truncated


@pytest.mark.parametrize("compress", ["gz", "bz2", "xz"])
def test_compressed_vasprun_parses_identically(tmp_path, compress):
    plain = _run(tmp_path, name="plain", n_frames=3)
    packed = _run(tmp_path, name="packed", n_frames=3, compress=compress)
    assert packed["vasprun"].name == f"vasprun.xml.{compress}"
    assert vr.detect_compression(packed["vasprun"]) == compress
    _, steps_a, trailer_a = _read(plain)
    _, steps_b, trailer_b = _read(packed)
    assert trailer_b.compression == compress and trailer_a.compression is None
    assert trailer_b.source_sha256 == hashlib.sha256(packed["vasprun"].read_bytes()).hexdigest()
    for a, b in zip(steps_a, steps_b):
        assert a.energies == b.energies
        np.testing.assert_array_equal(a.forces, b.forces)


def test_truncated_bz2_single_block_is_unreadable_not_silent(tmp_path):
    """bz2 decompresses whole 900 kB blocks: a cut inside the only block yields no bytes at all."""
    run = _run(tmp_path, n_frames=8, compress="bz2")
    data = run["vasprun"].read_bytes()
    run["vasprun"].write_bytes(data[: int(len(data) * 0.8)])
    with pytest.raises(vr.VasprunParseError, match="Compressed file ended"):
        _read(run)


@pytest.mark.parametrize("compress", ["gz", "xz"])
def test_truncated_compressed_stream_keeps_complete_steps(tmp_path, compress):
    run = _run(tmp_path, n_frames=8, compress=compress)
    data = run["vasprun"].read_bytes()
    run["vasprun"].write_bytes(data[: int(len(data) * 0.8)])
    _header, steps, trailer = _read(run)
    assert trailer.truncated and not trailer.closed and trailer.error
    complete = [s for s in steps if s.complete]
    assert 1 <= len(complete) < 8
    expected = run["expected"]["steps"]
    for step in complete:
        np.testing.assert_array_equal(step.forces, expected[step.index]["forces"])
    assert all("truncated_frame" in _problem_codes(s.problems) for s in steps if not s.complete)
    assert trailer.source_sha256 == hashlib.sha256(run["vasprun"].read_bytes()).hexdigest()


def test_data_after_modeling_is_not_truncation(tmp_path):
    run = _run(tmp_path, n_frames=2)
    with run["vasprun"].open("ab") as handle:
        handle.write(b"<junk/>\n")
    _header, steps, trailer = _read(run)
    assert trailer.closed and not trailer.truncated and len(steps) == 2
    assert any(p.startswith("info: data after </modeling>") for p in trailer.problems)


def test_compressed_text_evidence_is_read(tmp_path):
    run = _run(tmp_path, n_frames=2, incar={"IBRION": 2})
    outcar_gz = run["dir"] / "OUTCAR.gz"
    outcar_gz.write_bytes(gzip.compress(run["outcar"].read_bytes(), mtime=0))
    oszicar_bz2 = run["dir"] / "OSZICAR.bz2"
    oszicar_bz2.write_bytes(bz2.compress(run["oszicar"].read_bytes()))
    assert vf.parse_outcar(outcar_gz).free_energies == vf.parse_outcar(run["outcar"]).free_energies
    assert vf.parse_oszicar(oszicar_bz2).ionic_lines == vf.parse_oszicar(run["oszicar"]).ionic_lines


# --------------------------------------------------------------------------
# selective dynamics
# --------------------------------------------------------------------------

def _flags():
    flags = np.ones((8, 3), dtype=bool)
    flags[0] = flags[1] = False  # fully fixed
    flags[2] = [True, True, False]  # partially fixed (FixScaled-like)
    return flags


def test_selective_dynamics_from_initialpos_with_partial_flags(tmp_path):
    flags = _flags()
    run = _run(tmp_path, n_frames=2, selective=flags, incar={"IBRION": 2})
    header, steps, trailer = _read(run)
    np.testing.assert_array_equal(header.initial_structure.selective, flags)
    poscar = vf.parse_poscar(run["poscar"])
    np.testing.assert_array_equal(poscar.selective, flags)
    resolved = vr.resolve_selective(header, trailer, poscar=poscar)
    assert resolved["source"] == "initialpos" and resolved["basis"] == "direct"
    assert resolved["sources_with_flags"] == ["initialpos", "finalpos", "POSCAR"]
    assert resolved["conflicts"] == []
    np.testing.assert_array_equal(resolved["flags"], flags)
    assert (resolved["n_fixed_atoms"], resolved["n_partially_fixed_atoms"], resolved["n_fixed_components"]) == (2, 1, 7)
    # raw forces are untouched on fixed atoms
    np.testing.assert_array_equal(steps[0].forces, run["expected"]["steps"][0]["forces"])
    assert np.abs(steps[0].forces[:2]).sum() > 0
    assert vr.selective_flags(header, trailer)[1] == "initialpos"


def test_selective_dynamics_precedence_fallbacks(tmp_path):
    flags = _flags()
    run = _run(tmp_path, name="v522", n_frames=2, selective=flags, selective_in=("finalpos",))  # VASP 5.2.2 layout
    header, _steps, trailer = _read(run)
    resolved = vr.resolve_selective(header, trailer, poscar=vf.parse_poscar(run["poscar"]))
    assert resolved["source"] == "finalpos" and resolved["sources_without_flags"] == ["initialpos"]
    assert resolved["conflicts"] == []
    run = _run(tmp_path, name="poscar_only", n_frames=2, selective=flags, selective_in=())
    header, _steps, trailer = _read(run)
    resolved = vr.resolve_selective(header, trailer, poscar=vf.parse_poscar(run["poscar"]))
    assert resolved["source"] == "POSCAR"
    np.testing.assert_array_equal(resolved["flags"], flags)
    resolved = vr.resolve_selective(header, trailer, contcar=flags)
    assert resolved["source"] == "CONTCAR"
    none = vr.resolve_selective(header, trailer)
    assert none["flags"] is None and none["source"] is None


def test_selective_dynamics_disagreement_is_a_conflict(tmp_path):
    flags = _flags()
    other = flags.copy()
    other[5] = [False, False, False]
    run = _run(tmp_path, n_frames=2, selective=flags, poscar_selective=other)
    header, _steps, trailer = _read(run)
    resolved = vr.resolve_selective(header, trailer, poscar=vf.parse_poscar(run["poscar"]))
    assert resolved["source"] == "initialpos"
    kinds = {(c["a"], c["b"], c["kind"]) for c in resolved["conflicts"]}
    assert kinds == {("initialpos", "POSCAR", "flags_differ"), ("finalpos", "POSCAR", "flags_differ")}
    wrong_shape = vr.resolve_selective(header, trailer, contcar=np.ones((7, 3), dtype=bool))
    assert any(c["kind"] == "shape_differs" for c in wrong_shape["conflicts"])


def test_selective_dynamics_ase_constraints_do_not_touch_raw_forces(tmp_path):
    ase_io = pytest.importorskip("ase.io")
    flags = _flags()
    run = _run(tmp_path, n_frames=1, selective=flags, incar={"IBRION": -1})
    _header, steps, _ = _read(run)
    atoms = ase_io.read(run["vasprun"], index=0, format="vasp-xml")
    assert atoms.constraints  # ASE turns the flags into constraints on read ...
    np.testing.assert_array_equal(atoms.get_forces(apply_constraint=False), steps[0].forces)  # ... raw stays raw


# --------------------------------------------------------------------------
# magnetism: MAGMOM, total and site moments per step
# --------------------------------------------------------------------------

def test_spin_polarized_with_magmom_and_moment_jump(tmp_path):
    species, cell, frames = _nio(3)
    afm = [1.7, -1.7, 1.7, -1.7, 0.0, 0.0, 0.0, 0.0]
    flipped = [1.7, 1.7, 1.7, -1.7, 0.0, 0.0, 0.0, 0.0]
    for frame, mag, sites in zip(frames, (0.0, 0.0, 3.4), (afm, afm, flipped)):
        frame.update(mag_total=mag, site_moments=sites)
    run = fx.write_vasp_run(tmp_path / "afm", species=species, cell=cell, frames=frames,
                            incar={"IBRION": 0, "ISPIN": 2, "MAGMOM": "2 -2 2 -2 4*0", "LORBIT": 11})
    header, _steps, _ = _read(run)
    assert header.parameters["ISPIN"] == 2
    assert header.incar["MAGMOM"] == [2.0, -2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    assert header.parameters["MAGMOM"] == [2.0, -2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.total_magnetizations == [0.0, 0.0, 3.4]
    assert sorted(outcar.magnetization_tables) == [0, 1, 2]
    assert outcar.magnetization_tables[0] == afm and outcar.magnetization_tables[2] == flipped
    assert outcar.final_magnetization_table == flipped
    assert outcar.executed_tags["ISPIN"] == "2" and outcar.executed_tags["LORBIT"] == "11"
    oszicar = vf.parse_oszicar(run["oszicar"])
    values = vf.oszicar_step_values(oszicar, 3)
    assert values["mapped"] and values["mag"] == [0.0, 0.0, 3.4] and values["temperature_K"] == [300.0] * 3
    assert not vf.oszicar_step_values(oszicar, 4)["mapped"]
    incar = vf.parse_incar(run["incar"])
    assert vf.vasp_floats(incar.tags["MAGMOM"]) == [2.0, -2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]


def test_spin_polarized_without_magmom_uses_vasp_default(tmp_path):
    run = _run(tmp_path, n_frames=2, frame_fields={"mag_total": 8.0}, incar={"IBRION": 0, "ISPIN": 2})
    header, _steps, _ = _read(run)
    assert header.parameters["ISPIN"] == 2
    assert "MAGMOM" not in header.incar  # not set by the user ...
    assert header.parameters["MAGMOM"] == [1.0] * 8  # ... VASP's ferromagnetic default start
    assert "MAGMOM" not in vf.parse_incar(run["incar"]).tags
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.total_magnetizations == [8.0, 8.0] and outcar.magnetization_tables == {}


def test_non_spin_polarized_run_has_no_moment_values(tmp_path):
    run = _run(tmp_path, n_frames=2, incar={"IBRION": 2})
    outcar = vf.parse_outcar(run["outcar"])
    assert outcar.total_magnetizations == [None, None]
    assert all("mag" not in line for line in vf.parse_oszicar(run["oszicar"]).ionic_lines)


def test_site_moment_tables_are_capped(tmp_path):
    species, cell, frames = _nio(4, mag_total=0.0, site_moments=[1.7, -1.7, 1.7, -1.7, 0, 0, 0, 0])
    run = fx.write_vasp_run(tmp_path / "cap", species=species, cell=cell, frames=frames,
                            incar={"IBRION": 0, "ISPIN": 2, "LORBIT": 11})
    outcar = vf.parse_outcar(run["outcar"], max_site_moment_values=16)
    assert sorted(outcar.magnetization_tables) == [0, 1]
    assert outcar.final_magnetization_table is not None
    assert any("2 per-step magnetization (x) tables not retained" in p for p in outcar.problems)
    assert len(vf.parse_outcar(run["outcar"], max_site_moment_values=None).magnetization_tables) == 4


# --------------------------------------------------------------------------
# OUTCAR cross-check values
# --------------------------------------------------------------------------

def test_outcar_cross_check_values_match_vasprun(tmp_path):
    run = _run(tmp_path, n_frames=3, incar={"IBRION": 2})
    header, steps, _ = _read(run)
    outcar = vf.parse_outcar(run["outcar"], keep_forces=True)
    assert outcar.version == "6.4.2" and outcar.nions == 8 and outcar.ions_per_type == [4, 4]
    assert outcar.potcar_titles == run["expected"]["titles"]
    assert [h["titel"] for h in outcar.potcar_headers] == run["expected"]["titles"]
    assert outcar.completed and outcar.ionic_converged
    for k, (step, exp) in enumerate(zip(steps, run["expected"]["steps"])):
        d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
        assert abs(outcar.free_energies[k] - d["free_energy"]) <= 5e-9
        assert abs(outcar.energies_sigma0[k] - d["energy_sigma0"]) <= 5e-8
        assert abs(outcar.force_abs_sums[k] - np.abs(step.forces).sum()) <= 24 * 5e-7
        np.testing.assert_allclose(outcar.forces[k], step.forces, rtol=0, atol=5e-7)
        s = exp["stress_kbar"]
        assert outcar.stress_kbar[k] == pytest.approx([s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[2, 0]], abs=5e-6)
        assert outcar.last_energy_changes[k] is not None
    assert outcar.executed_tags["NELM"] == "60" and outcar.executed_tags["IBRION"] == "2"


def test_outcar_digit_minus_repair_and_fortran_floats():
    assert vf.fortran_float("0.2737684-111") == pytest.approx(0.2737684e-111)
    assert vf.fortran_float("1.5D+02") == 150.0
    with pytest.raises(ValueError):
        vf.fortran_float("*****")
    assert vf._repair_merged_numbers("-12353.08821-12353.08821").split() == ["-12353.08821", "-12353.08821"]


# --------------------------------------------------------------------------
# species / POTCAR / POSCAR
# --------------------------------------------------------------------------

def test_potcar_fingerprints_without_content(tmp_path):
    for mode, verified in ((True, True), (False, None), ("bad", False)):
        run = _run(tmp_path, name=f"pot_{mode}", n_frames=1, potcar_hash=mode)
        evidence = vf.parse_potcar(run["potcar"])
        assert evidence.sha256 == hashlib.sha256(run["potcar"].read_bytes()).hexdigest()
        assert [d["element"] for d in evidence.datasets] == ["Ni", "O"]
        assert [d["sha256_header_verified"] for d in evidence.datasets] == [verified, verified]
        for dataset, expected in zip(evidence.datasets, run["expected"]["potcar"]):
            assert dataset["titel"] == expected["titel"] and dataset["zval"] == expected["zval"]
            assert dataset["sha256_dataset_bytes"] == expected["sha256_dataset_bytes"]
            assert set(dataset) == {"titel", "vrhfin", "lexch", "pomass", "zval", "enmax", "enmin", "symbol",
                                    "element", "sha256_header", "sha256_dataset_bytes", "sha256_header_verified",
                                    "sha256_verify_variant"}
            assert not any("COPYR" in str(value) for value in dataset.values())
        assert evidence.problems == []


def test_species_order_mismatch_is_exposed(tmp_path):
    run = _run(tmp_path, n_frames=1, potcar_file_order=["O", "Ni"], poscar_species=["O", "Ni"])
    header, _steps, _ = _read(run)
    potcar = vf.parse_potcar(run["potcar"])
    poscar = vf.parse_poscar(run["poscar"])
    outcar = vf.parse_outcar(run["outcar"])
    vasprun_order = [t.element for t in header.atom_types]
    assert vasprun_order == ["Ni", "O"]
    assert [d["element"] for d in potcar.datasets] == ["O", "Ni"] != vasprun_order
    assert poscar.species == ["O", "Ni"] != vasprun_order
    assert [vf.normalize_species_token(t.split()[1]) for t in outcar.potcar_titles] == vasprun_order
    assert potcar.problems == []  # the POTCAR itself is internally consistent


def test_poscar_hashed_species_tokens_and_selective(tmp_path):
    run = _run(tmp_path, n_frames=1, poscar_species=["Ni_pv/6a2f546d", "O"], selective=_flags())
    poscar = vf.parse_poscar(run["poscar"])
    assert poscar.species_tokens == ["Ni_pv/6a2f546d", "O"] and poscar.species == ["Ni", "O"]
    assert poscar.counts == [4, 4] and poscar.coordinate_mode == "direct" and poscar.problems == []
    np.testing.assert_array_equal(poscar.selective, _flags())
    np.testing.assert_allclose(poscar.positions, run["expected"]["steps"][0]["positions"], atol=1e-12)


def test_poscar_cartesian_and_scale_factors(tmp_path):
    text = "\n".join([
        "SYNTHETIC TEST FIXTURE - not real VASP input",
        "  2.0",
        "  1.0 0.0 0.0", "  0.0 1.0 0.0", "  0.0 0.0 1.5",
        "  Ni O", "  1 1",
        "Selective dynamics",
        "Cartesian",
        "  0.0 0.0 0.0  F F F",
        "  1.0 1.0 1.5  T T F",
        "",
    ])
    path = tmp_path / "POSCAR"
    path.write_text(text, encoding="utf-8")
    poscar = vf.parse_poscar(path)
    np.testing.assert_allclose(poscar.cell, np.diag([2.0, 2.0, 3.0]))
    np.testing.assert_allclose(poscar.positions, [[0, 0, 0], [2.0, 2.0, 3.0]])
    np.testing.assert_allclose(poscar.fractional, [[0, 0, 0], [1.0, 1.0, 1.0]])
    assert poscar.selective.tolist() == [[False, False, False], [True, True, False]]
    path.write_text(text.replace("  2.0\n", "  -12.0\n", 1), encoding="utf-8")  # negative scale = volume
    assert abs(np.linalg.det(vf.parse_poscar(path).cell)) == pytest.approx(12.0)


def test_incar_parser_rules(tmp_path):
    path = tmp_path / "INCAR"
    path.write_text("\n".join([
        "# SYNTHETIC TEST FIXTURE - not real VASP input",
        "SYSTEM = NiO slab ! trailing comment",
        "ISPIN = 2 ; MAGMOM = 2*2.0 2*-2.0 4*0",
        "ENCUT = 520   # comment",
        "LDAU = .TRUE.",
        "LDAUU = 0 4.6 \\",
        "  0",
        "ispin = 2",
        "this line is not an assignment",
        "",
    ]), encoding="utf-8")
    incar = vf.parse_incar(path)
    assert incar.tags["SYSTEM"] == "NiO slab" and incar.tags["ENCUT"] == "520"
    assert vf.expand_vasp_list(incar.tags["MAGMOM"]) == ["2.0", "2.0", "-2.0", "-2.0", "0", "0", "0", "0"]
    assert incar.tags["LDAUU"].split() == ["0", "4.6", "0"]
    assert vf.vasp_bool(incar.tags["LDAU"]) is True and vf.vasp_bool("F") is False and vf.vasp_bool("x") is None
    assert vf.vasp_number(incar.tags["ENCUT"]) == 520.0
    assert incar.duplicates == ["ISPIN"]
    assert any(p.startswith("bad_value: INCAR line 9") for p in incar.problems)


# --------------------------------------------------------------------------
# OSZICAR line kinds
# --------------------------------------------------------------------------

def test_oszicar_md_and_relaxation_lines(tmp_path):
    md = _run(tmp_path, name="md", n_frames=3, frame_fields={"mag_total": 1.5}, incar={"IBRION": 0, "ISPIN": 2})
    relax = _run(tmp_path, name="relax", n_frames=3, incar={"IBRION": 2})
    md_osz = vf.parse_oszicar(md["oszicar"])
    assert md_osz.line_kinds == ["md"] * 3 and md_osz.problems == []
    assert [line["step"] for line in md_osz.ionic_lines] == [1, 2, 3]
    assert [line["T"] for line in md_osz.ionic_lines] == [300.0] * 3
    assert [line["mag"] for line in md_osz.ionic_lines] == [1.5] * 3
    for line, exp in zip(md_osz.ionic_lines, md["expected"]["steps"]):
        assert line["F"] == pytest.approx(exp["free_energy_input"], rel=1e-7)  # 8 significant digits only
        assert line["EK"] == pytest.approx(0.1)
    relax_osz = vf.parse_oszicar(relax["oszicar"])
    assert relax_osz.line_kinds == ["relax"] * 3
    assert set(relax_osz.ionic_lines[0]) == {"step", "F", "E0", "dE"}
    assert relax_osz.scf_steps == [12, 12, 12]


# --------------------------------------------------------------------------
# REAL VASP output shipped with ASE (read in place; amendments item 15)
# --------------------------------------------------------------------------

def test_real_vasprun_pstress():
    path = _ase_file("vasprun_pstress.xml")
    header, steps, trailer = vr.read_vasprun(path)
    assert header.generator["version"] == "6.3.0" and header.version_tuple == (6, 3, 0)
    assert header.species == ["C", "C"] and header.atom_types[0].pseudopotential == "PAW C 22Mar2012"
    assert header.parameters["PSTRESS"] == 1.0 and header.parameters["ISPIN"] == 1
    assert header.parameters["NELM"] == 60 and header.parameters["EDIFF"] == 1e-4
    assert header.parameters["MAGMOM"] == [1.0, 1.0] and header.parameters["LEPSILON"] is False
    assert header.problems == [] and sorted(header.parameter_duplicates) == ["CSHIFT", "OMEGAMAX"]
    assert trailer.closed and not trailer.truncated and trailer.final_structure_present
    assert len(steps) == 1 and len(steps[0].scf_energies) == 9
    step = steps[0]
    assert step.energies["e_fr_energy"] == -20.23985977
    d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
    assert d["free_energy"] == -20.246959373169393  # calc e_fr - PV, bit-identical to ASE
    assert d["energy_sigma0"] == pytest.approx(-20.247200953169394, abs=1e-12)
    assert d["pv_term"] == pytest.approx(0.0070996, abs=1e-7)
    assert d["energy_rule"] == "calc_level_direct" and d["calc_level_pattern"] == "consistent"
    assert vr.scf_summary(step)["last_dE"] == pytest.approx(8.0e-8, abs=1e-12)
    np.testing.assert_allclose(vasp_stress_to_ase(step.stress_kbar_vasp).diagonal(), [0.10360632132874024] * 3,
                               rtol=1e-12)
    assert np.all(step.forces == 0.0)
    ase_io = pytest.importorskip("ase.io")
    atoms = ase_io.read(path, index=-1, format="vasp-xml")
    assert atoms.get_potential_energy(force_consistent=True) == d["free_energy"]
    assert atoms.get_potential_energy() == d["energy_sigma0"]
    np.testing.assert_array_equal(atoms.get_forces(apply_constraint=False), step.forces)


def test_real_vasprun_dfpt_exposes_response_run_evidence():
    path = _ase_file("vasprun_dfpt.xml")
    header, steps, trailer = vr.read_vasprun(path)
    assert header.version_tuple == (6, 3, 2) and header.species == ["Na", "Cl"]
    assert header.parameters["LEPSILON"] is True and header.incar["LEPSILON"] is True  # -> non-training run
    assert trailer.closed and len(steps) == 1 and len(steps[0].scf_energies) == 80
    step = steps[0]
    assert step.energies["e_fr_energy"] == -6.74587304
    assert step.scf_energies[-1]["e_fr_energy"] == pytest.approx(-0.21489544)  # response iterations
    d = vr.derived_energies(step, header.parameters["PSTRESS"], header.version_tuple)
    assert d["free_energy"] == -6.74587304
    assert abs(d["additive_correction"]) > 1.0  # trailing scsteps are not the ground state
    np.testing.assert_allclose(vasp_stress_to_ase(step.stress_kbar_vasp).diagonal(), [-0.022525029251756264] * 3,
                               rtol=1e-12)


def test_real_outcar_example_1():
    path = _ase_file("OUTCAR_example_1")
    outcar = vf.parse_outcar(path, keep_forces=True)
    assert outcar.version.startswith("5.3.3") and outcar.version_tuple == (5, 3, 3)
    assert outcar.nions == 18 and outcar.ions_per_type == [18]
    assert outcar.potcar_titles == ["PAW_PBE Ni 02Aug2007"]
    assert outcar.potcar_headers == []  # NWRITE=0: no TITEL/VRHFIN echo
    assert outcar.free_energies == [-68.22868532]
    assert outcar.energies_no_entropy == [-68.23570214] and outcar.energies_sigma0 == [-68.23102426]
    assert outcar.scf_marker_kinds == ["reached"] and outcar.scf_iteration_counts == [63]
    assert outcar.total_magnetizations == [pytest.approx(17.0482806)]
    assert outcar.magnetization_tables == {} and outcar.final_magnetization_table is None  # LORBIT=0
    assert outcar.stress_kbar == [[-4.29429, -4.58894, -4.50342, 0.50047, -0.94049, 0.36481]]
    assert outcar.completed and not outcar.ionic_converged and outcar.ml_steps == 0
    assert outcar.problems == []
    forces = outcar.forces[0]
    np.testing.assert_array_equal(forces[0], [0.030415, -0.114705, -0.460114])
    np.testing.assert_array_equal(forces[17], [-0.00017, -0.427894, -0.52364])
    assert np.linalg.norm(forces, axis=1).max() == pytest.approx(2.395914408832043, rel=1e-12)
    tags = outcar.executed_tags
    assert (tags["NELM"], tags["EDIFF"], tags["ISPIN"], tags["LORBIT"], tags["NWRITE"]) == ("120", "0.1E-03", "2", "0", "0")
    markers = vf.outcar_scf_markers(outcar, n_dft_steps=1)
    assert markers["mapped"] and not markers["marker_is_proof"]  # VASP 5: informational only
    ase_io = pytest.importorskip("ase.io")
    atoms = ase_io.read(path, index=-1, format="vasp-out")
    assert atoms.get_potential_energy(force_consistent=True) == outcar.free_energies[0]
    assert atoms.get_potential_energy() == outcar.energies_sigma0[0]
    np.testing.assert_array_equal(atoms.get_forces(apply_constraint=False), forces)
    assert atoms.get_magnetic_moment() == pytest.approx(outcar.total_magnetizations[0])


@pytest.mark.parametrize(
    ("name", "kind"),
    [("convergence_OUTCAR_y_y", "reached"), ("convergence_OUTCAR_y_n", "not_reached"),
     ("convergence_OUTCAR_n_y", "reached")],
)
def test_real_convergence_fragments_are_open_steps(name, kind):
    """ASE's hand-trimmed OUTCAR fragments end inside a step: never counted, marker reported."""
    outcar = vf.parse_outcar(_ase_file(name))
    assert outcar.scf_marker_kinds == [] and outcar.free_energies == []
    assert any(p.startswith("truncated_frame: OUTCAR ends inside DFT ionic step 0") and repr(kind) in p
               for p in outcar.problems)
    assert outcar.executed_tags["EDIFF"] == "0.1E-03"
    marker_line = [line for line in _ase_file(name).read_text().splitlines() if "aborting loop" in line][0]
    assert vf.scf_marker_kind(marker_line) == kind


def test_real_poscar_example_1_is_vasp4_format():
    poscar = vf.parse_poscar(_ase_file("POSCAR_example_1"))
    assert poscar.species is None and poscar.counts == [18] and poscar.coordinate_mode == "cartesian"
    assert "missing: POSCAR has no species line (VASP 4 format)" in poscar.problems
    np.testing.assert_allclose(poscar.cell, np.eye(3) * 17.93435)
    assert poscar.selective is None
