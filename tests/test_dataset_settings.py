"""Reference-settings extraction, method fingerprints, overrides and element-aware pools.

All inputs are SYNTHETIC model objects built in this file (no VASP output is
parsed), except the real-format check that reads ASE's shipped
``vasprun_pstress.xml`` in place (skipped when ASE's test data is absent).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from nio_md_prep.dataset import settings as st
from nio_md_prep.dataset.errors import DatasetError
from nio_md_prep.dataset.model import AtomType, IncarEvidence, OutcarEvidence, PotcarEvidence, VasprunHeader

VALENCE = {"Ni": 10.0, "O": 6.0, "H": 1.0, "C": 4.0, "P": 5.0}
TITEL = {"Ni": "PAW_PBE Ni 02Aug2007", "O": "PAW_PBE O 08Apr2002", "H": "PAW_PBE H 15Jun2001",
         "C": "PAW_PBE C 08Apr2002", "P": "PAW_PBE P 06Sep2000"}


def header(blocks=(("Ni", 2), ("O", 2)), *, params=None, incar=None, version=(6, 4, 2), titles=None):
    """A SYNTHETIC VasprunHeader (no file behind it)."""
    atom_types = [AtomType(el, n, 1.0, VALENCE[el], (titles or TITEL)[el]) for el, n in blocks]
    species = [el for el, n in blocks for _ in range(n)]
    nelect = sum(VALENCE[el] * n for el, n in blocks)
    parameters = {
        "ENMAX": 520.0, "PREC": "accurate", "GGA": "PE", "LASPH": True, "ISPIN": 1, "LSORBIT": False,
        "LNONCOLLINEAR": False, "NUPDOWN": -1.0, "ISMEAR": 0, "SIGMA": 0.05, "NELECT": nelect,
        "LDAU": False, "LHFCALC": False, "NELM": 60, "EDIFF": 1e-6, "IBRION": 2, "NSW": 10, "ISIF": 2,
        "PSTRESS": 0.0, "POTIM": 0.3, "MAGMOM": [1.0] * len(species), "LDIPOL": False, "IDIPOL": 0,
    }
    parameters.update(params or {})
    return VasprunHeader(
        generator={"program": "vasp", "version": ".".join(map(str, version))},
        incar=dict(incar if incar is not None else {"ENCUT": 520.0, "IVDW": 11}),
        parameters=parameters, kpoints={"divisions": [1, 1, 1]}, species=species, atom_types=atom_types,
        version_tuple=version,
    )


def potcar(elements, *, sha_prefix="d", lexch="PE"):
    return PotcarEvidence(sha256="f" * 64, datasets=[
        {"element": el, "titel": TITEL[el], "lexch": lexch, "zval": VALENCE[el],
         "sha256_dataset_bytes": sha_prefix * 63 + str(i)} for i, el in enumerate(elements)
    ])


def outcar(**tags):
    return OutcarEvidence(sha256="0" * 64, version="6.4.2", nions=4, ions_per_type=[2, 2], potcar_titles=[],
                          scf_converged_markers=[], free_energies=[], energies_no_entropy=[], energies_sigma0=[],
                          magnetization_tables={}, total_magnetizations=[], ionic_converged=False, completed=True,
                          executed_tags={k: str(v) for k, v in tags.items()})


# --------------------------------------------------------------------------
# extraction
# --------------------------------------------------------------------------

def test_method_fields_come_from_parameters_with_encut_as_enmax():
    s = st.extract_settings(header())
    method = s["method"]
    assert method["ENCUT"] == 520.0 and s["provenance"]["fields"]["ENCUT"] == "parameters"
    assert method["PREC"] == "accurate" and method["GGA"] == "PE"
    assert method["ISPIN"] == 1 and method["NUPDOWN"] is None  # not applicable for ISPIN=1
    assert method["net_charge"] == 0.0
    assert method["IVDW"] == 11 and s["provenance"]["fields"]["IVDW"] == "vasprun_incar"
    assert method["hubbard"] == {} and method["LDAUTYPE"] is None
    assert method["potcar"]["Ni"] == {"titel": "PAW_PBE Ni 02Aug2007", "sha256": None}
    assert "MAGMOM" not in method  # magnetism is a separate record
    assert s["sampling"]["EDIFF"] == 1e-6 and s["sampling"]["NELM"] == 60


def test_ivdw_is_never_guessed():
    # not in <parameters>, not in <incar>, no OUTCAR/INCAR -> unknown
    s = st.extract_settings(header(incar={"ENCUT": 520.0}))
    assert s["method"]["IVDW"] == st.UNKNOWN and "IVDW" in s["unknown"]
    # the OUTCAR echo resolves it
    s = st.extract_settings(header(incar={"ENCUT": 520.0}), outcar=outcar(IVDW=12))
    assert s["method"]["IVDW"] == 12 and s["provenance"]["fields"]["IVDW"] == "outcar"
    # an INCAR file that does not set IVDW means the VASP default 0
    s = st.extract_settings(header(incar={"ENCUT": 520.0}), incar_file=IncarEvidence("1" * 64, {"ENCUT": "520"}))
    assert s["method"]["IVDW"] == 0 and s["provenance"]["fields"]["IVDW"] == "incar_file(not set)"


def test_gga_default_resolved_from_potcar_lexch_or_titel_family():
    s = st.extract_settings(header(params={"GGA": "--"}), potcar=potcar(["Ni", "O"]))
    assert s["method"]["GGA"] == "PE" and s["provenance"]["fields"]["GGA"] == "potcar_lexch"
    s = st.extract_settings(header(params={"GGA": "--"}))
    assert s["method"]["GGA"] == "PE" and s["provenance"]["fields"]["GGA"] == "potcar_titel_family"
    lda = {"Ni": "PAW Ni 02Aug2007", "O": "PAW O 08Apr2002"}
    assert st.extract_settings(header(params={"GGA": "--"}, titles=lda))["method"]["GGA"] == "CA"


def test_hubbard_is_per_element_and_ignores_u_zero():
    params = {"LDAU": True, "LDAUTYPE": 2, "LDAUL": [2, -1], "LDAUU": [4.6, 0.0], "LDAUJ": [0.0, 0.0]}
    s = st.extract_settings(header(params=params))
    assert s["method"]["hubbard"] == {"Ni": [2, 4.6, 0.0]} and s["method"]["LDAUTYPE"] == 2
    # LDAU on but arrays of the wrong length -> unknown, never a guess
    bad = dict(params, LDAUU=[4.6])
    s = st.extract_settings(header(params=bad))
    assert s["method"]["hubbard"] == st.UNKNOWN and "hubbard" in s["unknown"]


def test_potcar_dataset_hash_attached_only_when_potcar_matches_species():
    s = st.extract_settings(header(), potcar=potcar(["Ni", "O"]))
    assert s["method"]["potcar"]["Ni"]["sha256"].startswith("ddd")
    s = st.extract_settings(header(), potcar=potcar(["O", "Ni"]))  # wrong order: never attach hashes
    assert all(entry["sha256"] is None for entry in s["method"]["potcar"].values())
    assert "potcar_file" in s["provenance"]["hints"]


def test_ediff_printed_as_zero_is_taken_from_a_precise_source():
    s = st.extract_settings(header(params={"EDIFF": 0.0}), outcar=outcar(EDIFF="0.1E-08"))
    assert s["sampling"]["EDIFF"] == pytest.approx(1e-9) and s["provenance"]["sampling_fields"]["EDIFF"] == "outcar"
    s = st.extract_settings(header(params={"EDIFF": 0.0}))
    assert s["sampling"]["EDIFF"] == 0.0 and "EDIFF" in s["provenance"]["hints"]


def test_outcar_disagreement_is_recorded_as_conflict():
    s = st.extract_settings(header(), outcar=outcar(ENCUT="400.0"))
    conflicts = st.outcar_conflicts(s)
    assert [c["field"] for c in conflicts] == ["ENCUT"]


# --------------------------------------------------------------------------
# fingerprints, overrides, pools
# --------------------------------------------------------------------------

def test_fingerprint_is_stable_and_sensitive_to_method_only():
    a = st.extract_settings(header())
    b = st.extract_settings(header(params={"NELM": 200, "EDIFF": 1e-5, "POTIM": 0.5}))  # sampling only
    c = st.extract_settings(header(params={"ENMAX": 400.0}))
    assert st.method_fingerprint(a) == st.method_fingerprint(b)
    assert st.method_fingerprint(a) != st.method_fingerprint(c)
    assert len(st.method_fingerprint(a)) == 12


def test_pools_are_element_aware_and_split_on_conflicts():
    nio = st.extract_settings(header())
    water = st.extract_settings(header(blocks=(("H", 2), ("O", 1))))
    low_cut = st.extract_settings(header(params={"ENMAX": 400.0}))
    other_o = st.extract_settings(header(titles=dict(TITEL, O="PAW_PBE O_s 07Sep2000")))
    pools = st.build_pools([("r:nio", nio), ("r:water", water), ("r:low", low_cut), ("r:os", other_o)])
    members = sorted(sorted(p.run_ids) for p in pools)
    assert members == [["r:low"], ["r:nio", "r:water"], ["r:os"]]
    merged = next(p for p in pools if "r:nio" in p.run_ids)
    assert merged.method["elements"] == ["H", "Ni", "O"]
    diffs = st.describe_pool_differences(nio, low_cut)
    assert diffs[0]["field"] == "ENCUT" and diffs[0]["kind"] == "conflict"
    coverage = st.describe_pool_differences(nio, water)
    assert [d["kind"] for d in coverage] == ["coverage"]
    # input order does not matter
    again = st.build_pools([("r:os", other_o), ("r:low", low_cut), ("r:water", water), ("r:nio", nio)])
    assert sorted(p.pool_id for p in pools) == sorted(p.pool_id for p in again)


def test_unknown_never_pools_with_known_unless_overridden(tmp_path):
    known = st.extract_settings(header())
    unknown = st.extract_settings(header(incar={"ENCUT": 520.0}))  # IVDW unknown
    assert len(st.build_pools([("a", known), ("b", unknown)])) == 2
    path = tmp_path / "overrides.toml"
    path.write_text(
        '[[equivalence]]\nfield = "IVDW"\nvalues = ["11", "<unknown>"]\n'
        'reason = "SYNTHETIC test: INCAR on LONI verified IVDW=11"\nreviewed_by = "tester"\n', encoding="utf-8")
    overrides = st.load_overrides(path)
    assert overrides.as_dict()["equivalence"][0]["reviewed_by"] == "tester" and overrides.sha256
    assert len(st.build_pools([("a", known), ("b", unknown)], overrides=overrides)) == 1


@pytest.mark.parametrize("body, message", [
    ('[[equivalence]]\nfield = "ISPIN"\nvalues = ["1", "2"]\n', "missing keys"),
    ('[[equivalence]]\nfield = "NELM"\nvalues = ["60", "100"]\nreason = "x"\nreviewed_by = "y"\n', "not one of"),
    ('[[equivalence]]\nfield = "ISPIN"\nvalues = ["1", "1.0"]\nreason = "x"\nreviewed_by = "y"\n', "only one distinct"),
    ('[other]\n', "unknown top-level"),
])
def test_override_file_errors_are_explicit(tmp_path, body, message):
    path = tmp_path / "overrides.toml"
    path.write_text(body, encoding="utf-8")
    with pytest.raises(DatasetError, match=message):
        st.load_overrides(path)


def test_pool_unknown_fields_lists_missing_potcar_hashes():
    s = st.extract_settings(header(incar={"ENCUT": 520.0}))
    fields = st.pool_unknown_fields(s["method"])
    assert "IVDW" in fields and "potcar[Ni].sha256" in fields


# --------------------------------------------------------------------------
# real format (ASE's shipped real VASP 6.3.0 file, read in place)
# --------------------------------------------------------------------------

def test_real_vasprun_pstress_settings():
    ase = pytest.importorskip("ase")
    path = Path(ase.__file__).parent / "test" / "testdata" / "vasp" / "vasprun_pstress.xml"
    if not path.is_file():
        pytest.skip("ASE test data not installed")
    from nio_md_prep.dataset.vasprun import read_vasprun

    head, _, _ = read_vasprun(path)
    s = st.extract_settings(head)
    assert s["method"]["ENCUT"] == 400.0 and s["method"]["ISPIN"] == 1
    assert s["method"]["GGA"] == "CA"  # GGA='--' with LDA 'PAW C' POTCAR
    assert s["sampling"]["PSTRESS"] == 1.0 and s["sampling"]["NELM"] == 60
    assert s["method"]["potcar"]["C"]["titel"] == "PAW C 22Mar2012"


def test_pools_split_on_hubbard_u_but_not_on_elements_without_u():
    u = {"LDAU": True, "LDAUTYPE": 2, "LDAUL": [2, -1], "LDAUU": [4.6, 0.0], "LDAUJ": [0.0, 0.0]}
    nio_u = st.extract_settings(header(params=u))
    nio_u_other = st.extract_settings(header(params=dict(u, LDAUU=[6.2, 0.0])))
    nio_plain = st.extract_settings(header())
    water = st.extract_settings(header(blocks=(("H", 2), ("O", 1))))  # no Ni: +U on Ni does not apply
    pools = st.build_pools([("r:u", nio_u), ("r:u62", nio_u_other), ("r:plain", nio_plain), ("r:water", water)])
    members = sorted(sorted(p.run_ids) for p in pools)
    # greedy in run_id order: the Ni-free water run joins the first compatible pool (r:plain)
    assert members == [["r:plain", "r:water"], ["r:u"], ["r:u62"]]
    diffs = st.describe_pool_differences(nio_u, nio_u_other)
    assert [d["field"] for d in diffs if d["kind"] == "conflict"] == ["hubbard[Ni]"]
    assert st.assign_pools([("r:u", nio_u), ("r:water", water)]) == {
        p.pool_id: p.run_ids for p in st.build_pools([("r:u", nio_u), ("r:water", water)])}
