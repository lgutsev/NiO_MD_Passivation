"""Tests for dataset/extxyz.py.

All frames here are SYNTHETIC test data generated in-file (random geometry and
labels); they are not VASP output and carry no physical meaning.
"""

from __future__ import annotations

import io

import numpy as np
import pytest

ase = pytest.importorskip("ase")
import ase.io  # noqa: E402
from ase.calculators.calculator import all_properties  # noqa: E402
from ase.io import extxyz as ase_extxyz  # noqa: E402

from nio_md_prep.dataset import extxyz as X  # noqa: E402
from nio_md_prep.dataset.errors import DatasetError  # noqa: E402

SYNTHETIC = "SYNTHETIC-TEST-FIXTURE"


def _required_info(index: int = 0, *, stress: bool = True) -> dict:
    """The metadata X.REQUIRED_INFO_KEYS demands (SYNTHETIC values)."""
    info = {
        "source": f"synth:run/{index // 3}/vasprun.xml",
        "source_file_type": "vasprun",
        "ionic_step": index,
        "structure_key": f"{index:064x}",
        "lineage_group": f"synth:group/{index // 3}",
        "energy_source": "vasprun:calculation.e_fr_energy-PSTRESS*V",
        "energy_rule": "calc_level_direct",
        "stress_available": stress,
        "scf_status": "converged",
        "label_source": "dft",
        "vasp_version": "6.4.2",
        "magnetic_class": "controlled_consistent",
        "magnetic_policy": "accepted",
        "campaign": "synthetic",
        "family": "unknown",
        "parser": "nio-md-prep.vasprun-stream",
        "parser_version": "1.0",
        "repo_commit": "0" * 40,
    }
    if stress:
        info["stress_source"] = "vasprun:calculation.stress"
    else:
        info["stress_reason"] = "not_computed"
    return info


def _payload(index: int = 0, *, n_atoms: int = 6, seed: int = 0, **overrides) -> X.FramePayload:
    """A SYNTHETIC frame: triclinic cell, positions = fractional @ cell (full float64 digits)."""
    rng = np.random.default_rng(seed + index)
    cell = np.array([[8.34, 0.0, 0.0], [0.4170000000000001, 8.34, 0.0], [0.1, 0.2, 21.7]])
    fractional = rng.random((n_atoms, 3)) * 1.2 - 0.1  # a few atoms outside [0, 1): never wrapped
    species = ["Ni", "O"] * (n_atoms // 2) + ["P"] * (n_atoms % 2)
    selective = np.ones((n_atoms, 3), dtype=bool)
    selective[0] = False  # fully fixed atom
    selective[1, 2] = False  # partially fixed atom (FixScaled-like)
    fields = dict(
        frame_id=f"synth:run/{index // 3}#{index:05d}",
        species=species,
        cell=cell,
        positions=fractional @ cell,
        energy=-100.0 - index - 1.0 / 3.0,
        forces=np.round(rng.normal(size=(n_atoms, 3)), 8),
        stress=rng.normal(size=(3, 3)) * 1e-3,
        info={
            **_required_info(index),
            "generator": SYNTHETIC,
            "vasp_free_energy": -100.0 - index - 1.0 / 3.0,
            "vasp_energy_sigma0": -100.0 - index - 0.3,
            "run_id": f"synth:run/{index // 3}",
            "time_fs": 0.5 * (index + 1),
            "pool_id": "123456789012",  # all digits: ASE would read an int
        },
        selective_dynamics=selective,
        magmom_initial=np.array([2.0, 0.0, -2.0, 0.0, 2.0, 0.0][:n_atoms]),
    )
    fields.update(overrides)
    return X.FramePayload(**fields)


# --------------------------------------------------------------------------
# reserved keys
# --------------------------------------------------------------------------

def test_reserved_keys_cover_everything_ase_extxyz_converts():
    converted = set(all_properties) | set(ase_extxyz.per_atom_properties) | set(ase_extxyz.per_config_properties)
    structural = {"move_mask"} | set(ase_extxyz.SPECIAL_3_3_KEYS) | set(ase_extxyz.UNPROCESSED_KEYS)
    structural |= set(ase_extxyz.PROPERTY_NAME_MAP) | set(ase_extxyz.PROPERTY_NAME_MAP.values())
    from ase.outputs import all_outputs

    wanted = {name.lower() for name in converted | structural | set(all_outputs)}
    assert wanted <= X.RESERVED_KEYS, sorted(wanted - X.RESERVED_KEYS)
    for name in ("energy", "free_energy", "forces", "stress", "stresses", "energies", "magmom", "magmoms",
                 "charges", "dipole", "move_mask", "virial", "pbc", "lattice", "properties"):
        assert name in X.RESERVED_KEYS


@pytest.mark.parametrize("bad", ["energy", "Energy", "FORCES", "free_energy", "stress", "magmoms", "charges",
                                 "move_mask", "virial", "pbc", "Lattice", "dipole", "energies"])
def test_validate_label_keys_rejects_names_ase_would_convert(bad):
    with pytest.raises(DatasetError, match="reserved"):
        X.validate_label_keys(energy_key=bad)
    with pytest.raises(DatasetError, match="reserved"):
        X.validate_label_keys(forces_key=bad)


def test_validate_label_keys_rejects_invalid_or_duplicate_names():
    for bad in ["REF energy", "1abc", "", "a:b", "frame_id", "vasp_selective_dynamics"]:
        with pytest.raises(DatasetError):
            X.validate_label_keys(energy_key=bad)
    with pytest.raises(DatasetError, match="distinct"):
        X.validate_label_keys("REF_x", "ref_X", "REF_stress")
    keys = X.validate_label_keys("E_dft", "F_dft", "S_dft")
    assert keys.as_dict() == {"energy": "E_dft", "forces": "F_dft", "stress": "S_dft"}


def test_reserved_label_key_really_becomes_a_calculator_result_in_ase():
    """Why the check exists: ASE moves info 'energy'/arrays 'forces' into a calculator."""
    atoms = ase.Atoms("NiO", positions=[[0, 0, 0], [2.0, 0, 0]], cell=np.eye(3) * 4, pbc=True)
    atoms.info["energy"] = -1.0
    atoms.arrays["forces"] = np.ones((2, 3))
    buffer = io.StringIO()
    ase.io.write(buffer, atoms, format="extxyz", write_results=False)
    back = ase.io.read(io.StringIO(buffer.getvalue()), format="extxyz")
    assert "energy" not in back.info and "forces" not in back.arrays
    assert back.calc is not None and set(back.calc.results) >= {"energy", "forces"}


# --------------------------------------------------------------------------
# frame_to_atoms
# --------------------------------------------------------------------------

def test_frame_to_atoms_has_no_calculator_or_constraints_and_keeps_raw_forces():
    payload = _payload()
    atoms = X.frame_to_atoms(payload)
    assert atoms.calc is None
    assert atoms.constraints == []
    assert atoms.pbc.tolist() == [True, True, True]
    assert atoms.get_chemical_symbols() == list(payload.species)
    assert np.array_equal(atoms.positions, payload.positions)
    # raw forces on the fully fixed atom 0 are kept, not zeroed
    assert np.array_equal(atoms.arrays["REF_forces"], payload.forces)
    assert np.any(atoms.arrays["REF_forces"][0] != 0)
    flags = atoms.arrays[X.SELECTIVE_ARRAY]
    assert flags.dtype == bool and flags.shape == (6, 3)
    assert flags[0].tolist() == [False, False, False] and flags[1].tolist() == [True, True, False]
    assert atoms.info[X.SELECTIVE_BASIS_KEY] == "direct"
    assert atoms.info["REF_energy"] == payload.energy
    assert atoms.info["REF_stress"].shape == (3, 3)
    assert atoms.info["frame_id"] == payload.frame_id
    assert "vasp_magmom_final" not in atoms.arrays


def test_frame_to_atoms_accepts_mapping_and_custom_keys():
    payload = _payload()
    mapping = {name: getattr(payload, name) for name in X.FramePayload.__dataclass_fields__}
    keys = X.validate_label_keys("E_dft", "F_dft", "S_dft")
    atoms = X.frame_to_atoms(mapping, keys=keys)
    assert "E_dft" in atoms.info and "F_dft" in atoms.arrays and "S_dft" in atoms.info
    with pytest.raises(DatasetError, match="unknown frame payload fields"):
        X.frame_to_atoms({**mapping, "bogus": 1})


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"forces": np.full((6, 3), np.nan)}, "NaN"),
        ({"forces": np.zeros((5, 3))}, "shape"),
        ({"energy": float("inf")}, "finite"),
        ({"energy": True}, "finite"),
        ({"stress": np.zeros(6)}, "shape"),
        ({"cell": np.zeros((3, 3))}, "singular"),
        ({"selective_dynamics": np.ones((6, 3), dtype=int)}, "bool"),
        ({"info": {"energy": 1.0}}, "reserved"),
        ({"info": {"REF_forces": 1.0}}, "collides"),
        ({"info": {"frame_id": "x"}}, "collides"),
        ({"info": {"bad key": 1.0}}, "valid extxyz name"),
        ({"info": {"x": float("nan")}}, "finite"),
        ({"info": {"x": [1, "a"]}}, "list elements"),
        ({"info": {"x": object()}}, "unsupported"),
        ({"frame_id": "caf\u00e9"}, "ASCII"),
    ],
)
def test_invalid_payloads_are_rejected(overrides, message):
    with pytest.raises(DatasetError, match=message):
        X.frame_to_atoms(_payload(**overrides))


# --------------------------------------------------------------------------
# precision: ASE's writer vs ours
# --------------------------------------------------------------------------

def test_ase_default_writer_loses_position_precision_but_ours_round_trips_exactly(tmp_path):
    payload = _payload()
    atoms = X.frame_to_atoms(payload)
    ase_path = tmp_path / "ase_default.extxyz"
    ase.io.write(ase_path, atoms, format="extxyz")
    via_ase = ase.io.read(ase_path, format="extxyz")
    ase_error = np.max(np.abs(via_ase.positions - payload.positions))
    assert ase_error > 0.0  # %16.8f per-atom columns: up to 5e-9 A lost
    assert ase_error < 1e-8

    ours = tmp_path / "ours.extxyz"
    X.write_extxyz(ours, [payload])
    back = ase.io.read(ours, format="extxyz")
    assert np.array_equal(back.positions, payload.positions)
    assert np.array_equal(back.cell.array, payload.cell)
    assert np.array_equal(back.arrays["REF_forces"], payload.forces)
    assert back.info["REF_energy"] == payload.energy
    assert np.array_equal(back.info["REF_stress"], payload.stress)
    assert back.calc is None and back.constraints == []


# --------------------------------------------------------------------------
# write + verify
# --------------------------------------------------------------------------

def _awkward_info():
    return {
        "generator": SYNTHETIC,
        "digits": "123456789012",
        "exponent_like": "3e1234567890",
        "bool_word": "T",
        "empty": "",
        "json_like": "_JSON [1]",
        "windows_path": "C:\\Users\\x y\\OUTCAR",
        "quoted": 'say "hi" = {x} [y]',
        "unicode": "Ni\u2013O caf\u00e9",
        "newline": "a\nb",
        "posix_source": "synth:run/0/vasprun.xml",
        "count": 12,
        "negative_zero": -0.0,
        "big": 1e16,
        "flag": False,
        "axes": [2],
        "no_axes": [],
        "names": ["a", "b"],
        "matrix": [[1.0, 2.5], [3.0, 4.0]],
        "record": {"a": 1, "b": [1, 2], "c": "x"},
        "unknown": None,
    }


def test_round_trip_is_exact_for_every_field_and_awkward_info_values(tmp_path):
    payloads = [_payload(i) for i in range(3)]
    payloads[1] = _payload(1, info={**_required_info(1), **_awkward_info()},
                           magmom_final=np.array([1.71, 0.01, -1.69, 0.0, 1.7, -0.02]))
    payloads[2] = _payload(2, stress=None, selective_dynamics=None, magmom_initial=None,
                           info=_required_info(2, stress=False))
    path = tmp_path / "dataset.extxyz"
    result = X.write_extxyz(path, payloads)
    assert result.verified and result.frames == 3 and result.frame_ids == [p.frame_id for p in payloads]
    text = path.read_bytes()
    assert text.isascii() and b"\r\n" not in text

    report = X.read_back_and_verify(path, payloads)
    assert report.ok, report.as_dict()
    assert report.frames_read == 3 and report.info_keys_checked > 20
    assert all(value == 0.0 for value in report.max_abs_diff.values())

    frames = ase.io.read(path, index=":", format="extxyz")
    info = frames[1].info
    assert info["digits"] == "123456789012" and isinstance(info["digits"], str)
    assert info["exponent_like"] == "3e1234567890"
    assert info["bool_word"] == "T" and info["empty"] == "" and info["json_like"] == "_JSON [1]"
    assert info["windows_path"] == "C:\\Users\\x y\\OUTCAR"
    assert info["quoted"] == 'say "hi" = {x} [y]'
    assert info["unicode"] == "Ni\u2013O caf\u00e9" and info["newline"] == "a\nb"
    assert info["count"] == 12 and info["flag"] is False and info["big"] == 1e16
    assert info["axes"].tolist() == [2] and info["no_axes"].tolist() == [] and info["names"] == ["a", "b"]
    assert info["record"] == {"a": 1, "b": [1, 2], "c": "x"}
    assert "unknown" not in info  # None is omitted, never written as a bare (True) key
    assert np.array_equal(frames[1].arrays["vasp_magmom_final"], payloads[1].magmom_final)
    assert "REF_stress" not in frames[2].info and "stress_available" in frames[2].info
    assert X.SELECTIVE_ARRAY not in frames[2].arrays and X.SELECTIVE_BASIS_KEY not in frames[2].info
    for atoms in frames:
        assert atoms.calc is None and atoms.constraints == []


def test_write_is_byte_deterministic_and_independent_of_info_order(tmp_path):
    first = [_payload(i) for i in range(4)]
    second = []
    for payload in first:
        reversed_info = dict(reversed(list(payload.info.items())))
        second.append(_payload(int(payload.info["ionic_step"]), info=reversed_info))
    assert [list(p.info) for p in first] != [list(p.info) for p in second]
    a = X.write_extxyz(tmp_path / "a.extxyz", first)
    b = X.write_extxyz(tmp_path / "b.extxyz", second)
    assert (tmp_path / "a.extxyz").read_bytes() == (tmp_path / "b.extxyz").read_bytes()
    assert a.sha256 == b.sha256


def test_write_rejects_duplicate_frame_ids_and_publishes_nothing(tmp_path):
    path = tmp_path / "out.extxyz"
    with pytest.raises(DatasetError, match="duplicate frame_id"):
        X.write_extxyz(path, [_payload(0), _payload(0)])
    assert list(tmp_path.iterdir()) == []


def test_bad_payload_mid_stream_publishes_nothing(tmp_path):
    path = tmp_path / "out.extxyz"
    payloads = [_payload(0), _payload(1, forces=np.full((6, 3), np.inf))]
    with pytest.raises(DatasetError, match="NaN/inf"):
        X.write_extxyz(path, payloads)
    assert list(tmp_path.iterdir()) == []


def test_verification_catches_precision_loss_before_publishing(tmp_path, monkeypatch):
    """If the writer ever lost digits, verify=True must refuse to publish."""
    monkeypatch.setattr(X, "_repr_float", lambda value: "%.6f" % float(value))
    path = tmp_path / "out.extxyz"
    with pytest.raises(DatasetError, match="round-trip check"):
        X.write_extxyz(path, [_payload(0)])
    assert list(tmp_path.iterdir()) == []


def test_read_back_reports_field_mismatches_order_and_count(tmp_path):
    payloads = [_payload(i) for i in range(3)]
    path = tmp_path / "d.extxyz"
    X.write_extxyz(path, payloads)

    changed = [_payload(i) for i in range(3)]
    changed[1].forces = changed[1].forces * (1 + 1e-13)
    report = X.read_back_and_verify(path, changed)
    assert not report.ok and {m["field"] for m in report.mismatches} == {"REF_forces"}
    assert X.read_back_and_verify(path, changed, rtol=1e-12).ok
    with pytest.raises(DatasetError, match="round trip failed"):
        report.raise_if_failed()

    swapped = [payloads[1], payloads[0], payloads[2]]
    fields = {m["field"] for m in X.read_back_and_verify(path, swapped).mismatches}
    assert "frame_id" in fields

    short = X.read_back_and_verify(path, payloads[:2])
    assert not short.ok and short.frames_read == 3 and short.frames_expected == 2

    extra_info = [_payload(0, info={**payloads[0].info, "added": 1}), payloads[1], payloads[2]]
    assert "info_keys" in {m["field"] for m in X.read_back_and_verify(path, extra_info).mismatches}


def test_read_back_flags_a_label_that_ase_captured_as_calculator_result(tmp_path):
    path = tmp_path / "raw.extxyz"
    atoms = X.frame_to_atoms(_payload(0))
    atoms.info["energy"] = -5.0  # bypass validation: simulate a foreign file with a reserved key
    path.write_text(X.atoms_to_extxyz_block(atoms), encoding="ascii")
    report = X.read_back_and_verify(path, [_payload(0)])
    assert "calculator" in {m["field"] for m in report.mismatches}


# --------------------------------------------------------------------------
# text level: comment parsing, frame ids, slicing
# --------------------------------------------------------------------------

def test_parse_comment_line_tokenises_like_ase(tmp_path):
    path = tmp_path / "d.extxyz"
    X.write_extxyz(path, [_payload(1, info={**_required_info(1), **_awkward_info()})])
    comment = path.read_text(encoding="ascii").splitlines()[1]
    ours = X.parse_comment_line(comment)
    theirs = ase_extxyz.key_val_str_to_dict(comment)
    assert set(ours) == set(theirs)
    assert ours["frame_id"] == theirs["frame_id"]


def test_numeric_looking_frame_ids_survive_text_and_ase_paths(tmp_path):
    payloads = [_payload(0, frame_id="00012"), _payload(1, frame_id="T"), _payload(2, frame_id="1e5")]
    path = tmp_path / "ids.extxyz"
    X.write_extxyz(path, payloads)
    assert X.read_frame_ids(path) == ["00012", "T", "1e5"]
    assert [a.info["frame_id"] for a in ase.io.read(path, index=":", format="extxyz")] == ["00012", "T", "1e5"]


def test_slice_copies_requested_frames_bytewise_in_source_order(tmp_path):
    payloads = [_payload(i) for i in range(5)]
    src = tmp_path / "dataset.extxyz"
    X.write_extxyz(src, payloads)
    blocks = {block.frame_id: block.text for block in X.iter_extxyz_blocks(src)}
    dst = tmp_path / "test.extxyz"
    wanted = [payloads[3].frame_id, payloads[1].frame_id]
    result = X.slice_extxyz_by_frame_ids(src, wanted, dst)
    assert result["frame_ids"] == [payloads[1].frame_id, payloads[3].frame_id]
    assert dst.read_text(encoding="ascii") == blocks[payloads[1].frame_id] + blocks[payloads[3].frame_id]
    assert X.read_back_and_verify(dst, [payloads[1], payloads[3]]).ok

    empty = X.slice_extxyz_by_frame_ids(src, [], tmp_path / "empty.extxyz")
    assert empty["frames"] == 0 and (tmp_path / "empty.extxyz").read_bytes() == b""


def test_slice_fails_closed(tmp_path):
    payloads = [_payload(i) for i in range(3)]
    src = tmp_path / "dataset.extxyz"
    X.write_extxyz(src, payloads)
    dst = tmp_path / "out.extxyz"
    with pytest.raises(DatasetError, match="not in"):
        X.slice_extxyz_by_frame_ids(src, [payloads[0].frame_id, "synth:missing#00000"], dst)
    assert not dst.exists()
    with pytest.raises(DatasetError, match="twice"):
        X.slice_extxyz_by_frame_ids(src, [payloads[0].frame_id, payloads[0].frame_id], dst)
    doubled = tmp_path / "doubled.extxyz"
    doubled.write_bytes(src.read_bytes() * 2)
    with pytest.raises(DatasetError, match="more than once"):
        X.slice_extxyz_by_frame_ids(doubled, [payloads[0].frame_id], dst)
    truncated = tmp_path / "truncated.extxyz"
    truncated.write_bytes(src.read_bytes()[:-40])
    with pytest.raises(DatasetError, match="truncated"):
        X.slice_extxyz_by_frame_ids(truncated, [payloads[0].frame_id], dst)
    assert not dst.exists()


def test_label_sha256_is_stable_and_label_sensitive():
    payload = _payload()
    digest = X.label_sha256(payload.energy, payload.forces, payload.stress)
    assert digest == X.label_sha256(payload.energy, payload.forces.copy(), payload.stress.copy())
    assert digest != X.label_sha256(payload.energy, payload.forces)
    assert digest != X.label_sha256(np.nextafter(payload.energy, 0.0), payload.forces, payload.stress)


# --------------------------------------------------------------------------
# frame contract: required metadata, forces, stress, clean data
# --------------------------------------------------------------------------

def test_required_metadata_list_is_complete_and_enforced(tmp_path):
    wanted = {"source", "source_file_type", "ionic_step", "structure_key", "lineage_group", "energy_source",
              "energy_rule", "stress_available", "scf_status", "label_source", "vasp_version", "magnetic_class",
              "magnetic_policy", "campaign", "family", "parser", "parser_version", "repo_commit"}
    assert wanted <= set(X.REQUIRED_INFO_KEYS)
    for key in X.REQUIRED_INFO_KEYS:
        info = _required_info(0)
        info[key] = None  # None is missing: unknowns must be written explicitly
        with pytest.raises(DatasetError, match="required metadata missing"):
            X.write_extxyz(tmp_path / "x.extxyz", [_payload(0, info=info)])
    assert list(tmp_path.iterdir()) == []
    for key, bad in [("ionic_step", True), ("ionic_step", -1), ("ionic_step", "3"), ("stress_available", 1),
                     ("vasp_version", ""), ("lineage_group", 5)]:
        with pytest.raises(DatasetError, match=key):
            X.frame_to_atoms(_payload(0, info={**_required_info(0), key: bad}))
    # every required key is really written and read back by ASE
    path = tmp_path / "ok.extxyz"
    X.write_extxyz(path, [_payload(0)])
    back = ase.io.read(path, format="extxyz")
    assert set(X.REQUIRED_INFO_KEYS) <= set(back.info)
    assert back.info["label_set"] == "energy_forces" and back.info["lineage_group"] == "synth:group/0"
    assert back.info["structure_key"] == f"{0:064x}" and isinstance(back.info["structure_key"], str)
    assert back.info["repo_commit"] == "0" * 40 and isinstance(back.info["repo_commit"], str)


def test_missing_forces_never_enter_a_force_dataset_unless_energy_only_is_explicit(tmp_path):
    info = {**_required_info(0), "forces_reason": "no_forces_block"}
    energy_only = _payload(0, forces=None, info=info)
    with pytest.raises(DatasetError, match="no forces"):
        X.write_extxyz(tmp_path / "default.extxyz", [energy_only])
    assert list(tmp_path.iterdir()) == []

    contract = X.FrameContract(allow_energy_only=True)
    with pytest.raises(DatasetError, match="forces_reason"):
        X.write_extxyz(tmp_path / "e.extxyz", [_payload(0, forces=None)], contract=contract)
    path = tmp_path / "energy_only.extxyz"
    result = X.write_extxyz(path, [energy_only, _payload(1)], contract=contract)
    assert result.label_sets == {"energy_only": 1, "energy_forces": 1}
    comment = path.read_text(encoding="ascii").splitlines()[1]
    assert "REF_forces" not in comment and 'label_set="energy_only"' in comment
    frames = ase.io.read(path, index=":", format="extxyz")
    assert "REF_forces" not in frames[0].arrays and frames[0].info["label_set"] == "energy_only"
    assert frames[0].calc is None and frames[0].info["REF_energy"] == energy_only.energy
    assert frames[1].info["label_set"] == "energy_forces" and "REF_forces" in frames[1].arrays
    assert X.read_back_and_verify(path, [energy_only, _payload(1)]).ok


def test_payload_cannot_set_module_keys_split_or_label_set():
    for key in ("split", "label_set", "frame_id", "selective_dynamics_basis"):
        with pytest.raises(DatasetError, match="collides"):
            X.frame_to_atoms(_payload(0, info={**_required_info(0), key: "train"}))
    for bad in ("source", "lineage_group", "split", "label_set"):
        with pytest.raises(DatasetError, match="collides"):
            X.validate_label_keys(energy_key=bad)


def test_missing_stress_is_omitted_never_zero_and_must_be_explained(tmp_path):
    with pytest.raises(DatasetError, match="stress_available=True but the stress label is absent"):
        X.frame_to_atoms(_payload(0, stress=None, info={**_required_info(0, stress=False), "stress_available": True}))
    info = _required_info(0, stress=False)
    del info["stress_reason"]
    with pytest.raises(DatasetError, match="stress_reason"):
        X.frame_to_atoms(_payload(0, stress=None, info=info))
    with pytest.raises(DatasetError, match="stress_available=False but the stress label is present"):
        X.frame_to_atoms(_payload(0, info={**_required_info(0), "stress_available": False}))
    path = tmp_path / "nostress.extxyz"
    X.write_extxyz(path, [_payload(0, stress=None, info=_required_info(0, stress=False))])
    back = ase.io.read(path, format="extxyz")
    assert "REF_stress" not in back.info and back.info["stress_available"] is False
    assert back.info["stress_reason"] == "not_computed" and back.calc is None


@pytest.mark.parametrize("key, value", [("label_source", "mlff"), ("scf_status", "not_converged"),
                                        ("scf_status", "unknown")])
def test_clean_contract_refuses_mlff_and_unconverged_frames(tmp_path, key, value):
    payload = _payload(0, info={**_required_info(0), key: value})
    with pytest.raises(DatasetError, match=key):
        X.write_extxyz(tmp_path / "clean.extxyz", [payload])
    assert list(tmp_path.iterdir()) == []
    review = tmp_path / "review.extxyz"  # e.g. a quarantine file for human review: explicitly not clean
    X.write_extxyz(review, [payload], contract=X.FrameContract(clean=False))
    assert ase.io.read(review, format="extxyz").info[key] == value


def test_selective_dynamics_exported_as_raw_direct_flags_without_ase_constraints(tmp_path):
    payload = _payload(0)
    path = tmp_path / "sd.extxyz"
    X.write_extxyz(path, [payload])
    comment = path.read_text(encoding="ascii").splitlines()[1]
    assert "move_mask" not in comment and "vasp_selective_dynamics:L:3" in comment
    back = ase.io.read(path, format="extxyz")
    assert back.constraints == [] and back.calc is None
    assert np.array_equal(back.arrays["vasp_selective_dynamics"], payload.selective_dynamics)
    assert back.info["selective_dynamics_basis"] == "direct"
    assert np.array_equal(back.arrays["REF_forces"], payload.forces)  # fixed atom 0 keeps its raw force


# --------------------------------------------------------------------------
# text level: decoding and the split key
# --------------------------------------------------------------------------

def test_block_info_decoding_matches_ase_for_scalars_and_strings(tmp_path):
    path = tmp_path / "d.extxyz"
    info = {**_required_info(1), **_awkward_info()}
    X.write_extxyz(path, [_payload(1, info=info)])
    block = next(X.iter_extxyz_blocks(path))
    ours = block.info()
    theirs = ase.io.read(path, format="extxyz").info
    for key, value in info.items():
        if value is None:
            assert key not in ours
            continue
        if isinstance(value, (list, dict)):
            assert ours[key] == value
        else:
            assert ours[key] == value and type(ours[key]) is type(value), key
            assert ours[key] == theirs[key]
    assert ours["frame_id"] == block.frame_id and ours["label_set"] == "energy_forces"
    assert X.decode_comment_value("T") is True and X.decode_comment_value("-12") == -12
    assert X.decode_comment_value("1.5") == 1.5 and X.decode_comment_value("abc") == "abc"


def test_slice_adds_the_split_key_and_otherwise_keeps_bytes(tmp_path):
    payloads = [_payload(i) for i in range(4)]
    src = tmp_path / "dataset.extxyz"
    X.write_extxyz(src, payloads)
    blocks = {block.frame_id: block.text for block in X.iter_extxyz_blocks(src)}
    dst = tmp_path / "valid.extxyz"
    X.slice_extxyz_by_frame_ids(src, [payloads[2].frame_id], dst, add_info={"split": "valid"})
    text = dst.read_text(encoding="ascii")
    original = blocks[payloads[2].frame_id]
    assert text.replace(' split="valid"', "", 1) == original
    assert text.splitlines()[1].endswith('split="valid" pbc="T T T"')
    back = ase.io.read(dst, format="extxyz")
    assert back.info["split"] == "valid" and back.calc is None and back.constraints == []
    assert np.array_equal(back.positions, payloads[2].positions)
    assert np.array_equal(back.arrays["REF_forces"], payloads[2].forces)
    again = tmp_path / "again.extxyz"
    with pytest.raises(DatasetError, match="already has an info key"):
        X.slice_extxyz_by_frame_ids(dst, [payloads[2].frame_id], again, add_info={"split": "test"})
    for bad in ({"energy": 1.0}, {"pbc": "x"}, {"frame_id": "x"}, {"split": None}, {"split": ["a"]}):
        with pytest.raises(DatasetError):
            X.slice_extxyz_by_frame_ids(src, [payloads[0].frame_id], again, add_info=bad)
    assert not again.exists()
