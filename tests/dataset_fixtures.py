"""SYNTHETIC VASP run writer shared by the ``tests/test_dataset_*.py`` suites.

SYNTHETIC TEST FIXTURE - not real VASP output. Every file written here is
labelled synthetic in-file:

* ``vasprun.xml``: an XML comment ``SYNTHETIC TEST FIXTURE - not real VASP
  output`` and ``<generator>`` subversion ``SYNTHETIC-TEST-FIXTURE``;
* ``OUTCAR``/``OSZICAR``/``INCAR``/``POSCAR``/``KPOINTS``: a
  ``SYNTHETIC TEST FIXTURE`` line;
* ``POTCAR``: header lines only (TITEL/VRHFIN/LEXCH/POMASS/ZVAL/ENMAX/``End of
  Dataset``) with ``SYNTHETIC`` titles and a COPYR line saying it is fake. No
  licensed pseudopotential data is ever written.

The layouts follow the real-format facts in
``design/understand/vasp-format-research.md`` (S1 vasprun structure, S2 SCF
markers, S3 energies, S5 magnetization, S6 POTCAR headers, S7 OUTCAR/OSZICAR),
so the parsers can be exercised on edge cases no local real file contains
(PSTRESS, VASP <= 6.0.8 energy mislabel, truncation, MLFF flat steps, '*****').

Main entry point: :func:`write_vasp_run`. Everything is deterministic: no
randomness unless a ``seed`` is passed to :func:`make_frames`.

Energy bookkeeping (what the writer puts where, per DFT step; eV)::

    F      = frame["free_energy"]            force-consistent free energy (the label)
    E_wo   = F + frame["e_wo_offset"]         energy without entropy
    E0     = F + frame["e0_offset"]           energy(sigma->0)
    Edisp  = frame.get("edisp", 0.0)          additive dispersion term (IVDW runs)
    PV     = pstress[kB] * V[A^3] * 1e-22/1.60217733e-19   (VASP's constant)

    last <scstep>:        e_fr = F - Edisp, e_wo = e_fr + e_wo_offset, e_0 = e_fr + e0_offset
    <calculation><energy> e_fr = F + PV,    e_wo = E_wo + PV,           e_0 = E0 + PV
      (VASP <= 6.0.8 bug, reproduced when version <= 6.0.8 or calc_energy_bug=True:
       e_wo holds E0 + PV and e_0 holds F - E_wo)
    OUTCAR summary:       TOTEN = F, energy without entropy = E_wo, energy(sigma->0) = E0
    OSZICAR ionic line:   F= F, E0= E0 (8 significant digits)

The returned ``expected`` values are computed from the *formatted* numbers
that were written (8 decimals in vasprun), so a correct parser must reproduce
them to floating-point precision.
"""

from __future__ import annotations

import bz2
import gzip
import hashlib
import lzma
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

SYNTHETIC_COMMENT = "SYNTHETIC TEST FIXTURE - not real VASP output"
SYNTHETIC_SUBVERSION = "SYNTHETIC-TEST-FIXTURE"
#: VASP's eV->J constant; PV = PSTRESS[kB] * V[A^3] * PV_FACTOR
PV_FACTOR = 1e-22 / 1.60217733e-19

#: element -> (ZVAL, POMASS, ENMAX); plausible, clearly not a pseudopotential library
ELEMENT_DATA = {
    "H": (1.0, 1.000, 250.0), "C": (4.0, 12.011, 400.0), "N": (5.0, 14.001, 400.0),
    "O": (6.0, 16.000, 400.0), "Na": (1.0, 22.990, 102.0), "P": (5.0, 30.974, 255.0),
    "S": (6.0, 32.066, 259.0), "Cl": (7.0, 35.453, 262.0), "Ni": (10.0, 58.690, 270.0),
}

MARKER_REACHED = "------------------------ aborting loop because EDIFF is reached ----------------------------------------"
MARKER_NOT_REACHED = "------------------------ aborting loop EDIFF was not reached (unconverged)  ----------------------------"
MARKER_HARD_STOP = "------------------------ aborting loop because hard stop was set ---------------------------------------"
_DASHES = " " + "-" * 83


# --------------------------------------------------------------------------
# Structures and frames
# --------------------------------------------------------------------------

def rocksalt_nio(a: float = 4.17) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Conventional rocksalt Ni4O4 (species blocks Ni then O): ``(species, cell, cartesian positions)``."""
    fcc = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]], dtype=float)
    frac = np.vstack([fcc, fcc + [0.5, 0, 0]]) % 1.0
    cell = np.eye(3) * a
    return ["Ni"] * 4 + ["O"] * 4, cell, frac @ cell


def make_frames(
    species: Sequence[str],
    cell: Any,
    positions: Any,
    n_frames: int,
    *,
    seed: int | None = None,
    base_energy: float = -50.0,
    amplitude: float = 0.02,
    stress: bool = True,
    **frame_fields: Any,
) -> list[dict[str, Any]]:
    """Deterministic, plausible frames (small displacements, net-zero forces).

    With ``seed=None`` the displacements are a fixed trigonometric pattern;
    with a seed they come from ``numpy.random.default_rng(seed)``. Extra
    keyword arguments are copied into every frame (e.g. ``n_scf=12``).
    """
    n_atoms = len(species)
    base = np.asarray(positions, dtype=float)
    rng = np.random.default_rng(seed) if seed is not None else None
    frames = []
    for k in range(n_frames):
        if rng is not None:
            shift = rng.normal(scale=amplitude, size=(n_atoms, 3))
            forces = rng.normal(scale=0.3, size=(n_atoms, 3))
        else:
            grid = np.arange(n_atoms * 3, dtype=float).reshape(n_atoms, 3)
            shift = amplitude * np.sin(grid * 0.7 + k * 1.3)
            forces = 0.3 * np.cos(grid * 1.1 + k * 0.9)
        forces -= forces.mean(axis=0)  # VASP removes the net force
        frame = {
            "positions": base + shift,
            "forces": forces,
            "free_energy": base_energy - 0.01 * k + 0.001 * math.sin(k),
            "e_wo_offset": -0.004,
            "e0_offset": -0.002,
            "stress_kbar": (np.array([[-4.2, 0.5, 0.3], [0.5, -4.5, -0.9], [0.3, -0.9, -4.4]]) + 0.1 * k)
            if stress else None,
            "n_scf": 12,
            "mag_total": None,
            "site_moments": None,
            "scf_marker": "reached",
        }
        frame.update(frame_fields)
        frames.append(frame)
    return frames


def _species_blocks(species: Sequence[str]) -> list[tuple[str, int]]:
    blocks: list[tuple[str, int]] = []
    for element in species:
        if blocks and blocks[-1][0] == element:
            blocks[-1] = (element, blocks[-1][1] + 1)
        else:
            blocks.append((element, 1))
    return blocks


# --------------------------------------------------------------------------
# Number formatting
# --------------------------------------------------------------------------

def fortran_e(value: float, digits: int, *, leading_zero_negative: bool = False) -> str:
    """Fortran ``E`` style: ``0.1234567E-05``; negatives as ``-.1234567E-05`` (OSZICAR) or ``-0.12...``."""
    if value == 0 or not math.isfinite(value):
        return f"0.{'0' * digits}E+00" if value == 0 else str(value)
    exponent = math.floor(math.log10(abs(value))) + 1
    mantissa = round(abs(value) / 10 ** exponent, digits)
    if mantissa >= 1.0:
        mantissa /= 10
        exponent += 1
    text = f"{mantissa:.{digits}f}"
    if value < 0:
        text = "-" + (text if leading_zero_negative else text[1:])
    return f"{text}E{exponent:+03d}"


def _f8(value: float) -> str:
    return f"{value:16.8f}"


def _parse_back(text: str) -> float:
    return float(text)


def _row(values: Iterable[float], overflow: bool = False) -> str:
    cells = [_f8(value) for value in values]
    if overflow:
        cells[0] = "   *************"
    return " ".join(cells)


# --------------------------------------------------------------------------
# vasprun.xml pieces
# --------------------------------------------------------------------------

def _xml_i(name: str, value: Any, kind: str | None = None) -> str:
    if kind is None:
        if isinstance(value, bool):
            kind = "logical"
        elif isinstance(value, int):
            kind = "int"
        elif isinstance(value, str):
            kind = "string"
    if kind == "logical":
        text = " T  " if value else " F  "
    elif kind == "int":
        text = f"{int(value):6d}"
    elif kind == "string":
        text = str(value)
    else:
        text = _f8(float(value))
    type_attr = f' type="{kind}"' if kind else ""
    return f'<i{type_attr} name="{name}">{text}</i>'


def _xml_v(name: str, values: Sequence[Any], kind: str | None = None) -> str:
    if kind == "int":
        text = " ".join(f"{int(value):8d}" for value in values)
    elif kind == "logical":
        text = " ".join(" T " if value else " F " for value in values)
    else:
        text = " ".join(_f8(float(value)) for value in values)
    type_attr = f' type="{kind}"' if kind else ""
    return f'<v{type_attr} name="{name}">{text}</v>'


def _structure_xml(
    name: str | None, cell_text: list[list[str]], frac_text: list[list[str]], volume_text: str,
    selective: Any = None, velocities: bool = False, indent: str = " ",
) -> list[str]:
    attr = f' name="{name}" ' if name else ""
    out = [f"{indent}<structure{attr}>", f"{indent} <crystal>", f'{indent}  <varray name="basis" >']
    out += [f"{indent}   <v>{' '.join(row)} </v>" for row in cell_text]
    out += [f"{indent}  </varray>", f'{indent}  <i name="volume">{volume_text} </i>', f'{indent}  <varray name="rec_basis" >']
    cell = np.array([[float(x) for x in row] for row in cell_text])
    rec = np.linalg.inv(cell).T
    out += [f"{indent}   <v>{_row(r)} </v>" for r in rec]
    out += [f"{indent}  </varray>", f"{indent} </crystal>", f'{indent} <varray name="positions" >']
    out += [f"{indent}  <v>{' '.join(row)} </v>" for row in frac_text]
    out += [f"{indent} </varray>"]
    if selective is not None:
        out += [f'{indent} <varray type="logical" name="selective" >']
        out += [
            f'{indent}  <v type="logical" >{"".join(" T " if flag else " F " for flag in row)}</v>'
            for row in np.asarray(selective, dtype=bool)
        ]
        out += [f"{indent} </varray>"]
    if velocities:
        out += [f'{indent} <varray name="velocities" >']
        out += [f"{indent}  <v>{_row([0.001, -0.002, 0.003])} </v>" for _ in frac_text]
        out += [f"{indent} </varray>"]
    out += [f"{indent}</structure>"]
    return out


def _energy_items(values: dict[str, float], indent: str, overflow: Iterable[str] = ()) -> list[str]:
    lines = []
    for name, value in values.items():
        text = "   ************ " if name in overflow else f"{_f8(value)} "
        lines.append(f'{indent}<i name="{name}">{text}</i>')
    return lines


_SC_COMPONENTS = ("alphaZ", "ewald", "hartreedc", "XCdc", "pawpsdc", "pawaedc", "eentropy", "bandstr", "atom")


# --------------------------------------------------------------------------
# The writer
# --------------------------------------------------------------------------

def write_vasp_run(
    dirpath: Path,
    *,
    species: Sequence[str],
    cell: Any,
    frames: Sequence[dict[str, Any]],
    incar: dict[str, Any] | None = None,
    nelm: int = 60,
    ediff: float = 1e-5,
    version: str = "6.4.2",
    truncate_after_steps: int | None = None,
    truncate_inside_step: bool | str = False,
    write_outcar: bool = True,
    write_oszicar: bool = True,
    potcar_titles: Sequence[str] | None = None,
    selective: Any = None,
    mlff_steps: Iterable[int] = (),
    pstress: float = 0.0,
    compress: str | None = None,
    selective_in: Sequence[str] = ("initialpos", "finalpos"),
    calc_energy_bug: bool | None = None,
    potcar_hash: bool | str = True,
    outcar_potcar_headers: bool = True,
    poscar_species: Sequence[str] | None = None,
    raw_parameters: dict[str, tuple[str | None, str]] | None = None,
    ionic_converged: bool | None = None,
    write_potcar: bool = True,
    potcar_file_order: Sequence[str] | None = None,
    poscar_selective: Any = "same",
) -> dict[str, Any]:
    """Write a SYNTHETIC VASP run directory; return written paths and expected parse results.

    Parameters (``frames`` entries; all arrays Cartesian Angstrom / eV / eV/A):

    ``positions`` (N,3), ``forces`` (N,3), ``free_energy`` F, ``e_wo_offset``
    (E_wo - F), ``e0_offset`` (E0 - F), ``stress_kbar`` (3,3 VASP sign) or None
    (no stress varray, like ISIF=0), ``n_scf`` (number of ``<scstep>``),
    ``mag_total`` (OUTCAR/OSZICAR total moment; None -> not printed for
    ISPIN=1, 0.0 for ISPIN=2), ``site_moments`` (N,) OUTCAR ``magnetization
    (x)`` table or None, ``scf_marker`` ``'reached'``/``'not_reached'``/None
    (OUTCAR ``aborting loop`` line; also ``'hard_stop'``). Optional: ``cell`` (per-frame 3x3),
    ``edisp`` (additive dispersion energy), ``last_dE`` (|dE| between the last
    two scsteps; default EDIFF/10), ``d_eps``, ``kinetic`` (MD),
    ``overflow`` (set of ``"forces"``, ``"positions"``, ``"stress"``,
    ``"e_fr_energy"``: written as ``*****`` in vasprun), ``omit`` (set of
    ``"forces"``, ``"stress"``, ``"energy"``: the block is not written to
    vasprun; ``"forces"`` also drops the OUTCAR TOTAL-FORCE block; use
    ``n_scf=0`` for a calculation without ``<scstep>``).
    ``stress_kbar`` is ignored (no stress anywhere) when ``ISIF=0``, as in VASP.

    Run options: ``incar`` tags (written to INCAR, ``<incar>`` and, where VASP
    does so, ``<parameters>``; IVDW and ML_* are deliberately NOT put in
    ``<parameters>``, as in real VASP); ``nelm``/``ediff``; ``version`` (VASP
    <= 6.0.8 reproduces the calc-level e_wo/e_0 mislabel unless
    ``calc_energy_bug`` says otherwise); ``pstress`` (kB; adds PV to the
    calc-level energies); ``selective`` (N,3 bool, True = may move) written to
    POSCAR and to the vasprun structures named in ``selective_in``;
    ``mlff_steps`` (frame indices written as VASP-MLFF flat force-field steps;
    OUTCAR then gets an ML prediction block on every step);
    ``truncate_after_steps=k`` (only k steps complete, no finalpos, no
    ``</modeling>``, no OUTCAR timing) and ``truncate_inside_step`` (True: the
    file ends after step k's forces/stress; ``"mid_forces"``: inside its forces
    varray); ``compress`` ``None|'gz'|'bz2'|'xz'`` (vasprun only);
    ``raw_parameters`` ``{name: (type_attr, text)}`` appended verbatim to
    ``<parameters>`` (e.g. ``("int", "*****")``); ``potcar_hash`` True (valid
    SHA256 line), False (none) or ``"bad"``; ``poscar_species`` (override
    POSCAR species tokens, e.g. ``["Ni_pv/6a2f546d", "O"]``);
    ``potcar_file_order`` (element order of the POTCAR FILE only, e.g. ``["O",
    "Ni"]`` for a POTCAR that disagrees with the vasprun/OUTCAR species);
    ``poscar_selective`` (``"same"`` = ``selective``; None = no ``Selective
    dynamics`` in POSCAR; or an (N,3) bool array that differs from vasprun).

    Returns a dict with ``dir`` and the paths (``vasprun``, ``outcar``,
    ``oszicar``, ``potcar``, ``incar``, ``poscar``, ``kpoints``; None when not
    written) plus ``expected``: run-level facts and one ``steps`` entry per
    vasprun step written (see the module docstring).
    """
    dirpath = Path(dirpath)
    dirpath.mkdir(parents=True, exist_ok=True)
    species = list(species)
    n_atoms = len(species)
    blocks = _species_blocks(species)
    incar = dict(incar or {})
    mlff = set(mlff_steps)
    version_tuple = tuple(int(part) for part in version.split(".")[:3])
    bug = calc_energy_bug if calc_energy_bug is not None else version_tuple <= (6, 0, 8)
    ibrion = int(incar.get("IBRION", 0 if len(frames) > 1 else -1))
    nsw = int(incar.get("NSW", len(frames) if len(frames) > 1 else 0))
    isif = int(incar.get("ISIF", 2))
    potim = float(incar.get("POTIM", 1.0 if ibrion == 0 else 0.5))
    ispin = int(incar.get("ISPIN", 2 if any(f.get("mag_total") is not None for f in frames) else 1))
    lorbit = int(incar.get("LORBIT", 11 if any(f.get("site_moments") is not None for f in frames) else 0))
    ismear = int(incar.get("ISMEAR", 0))
    sigma = float(incar.get("SIGMA", 0.05))
    encut = float(incar.get("ENCUT", 400.0))
    is_md = ibrion == 0
    titles = list(potcar_titles) if potcar_titles is not None else [f"PAW_PBE {el} SYNTHETIC" for el, _ in blocks]
    if len(titles) != len(blocks):
        raise ValueError("potcar_titles must have one entry per species block")
    zvals = [ELEMENT_DATA.get(el, (1.0, 1.0, 300.0))[0] for el, _ in blocks]
    masses = [ELEMENT_DATA.get(el, (1.0, 1.0, 300.0))[1] for el, _ in blocks]
    nelect = float(incar.get("NELECT", sum(z * count for z, (_, count) in zip(zvals, blocks))))
    magmom = incar.get("MAGMOM")
    magmom_values = [float(x) for x in str(magmom).split()] if magmom is not None and "*" not in str(magmom) else None
    if magmom is not None and "*" in str(magmom):
        magmom_values = []
        for token in str(magmom).split():
            if "*" in token:
                count, value = token.split("*")
                magmom_values += [float(value)] * int(count)
            else:
                magmom_values.append(float(token))
    if magmom_values is None:
        magmom_values = [1.0] * n_atoms  # VASP default, present even for ISPIN=1
    n_written = len(frames) if truncate_after_steps is None else min(truncate_after_steps, len(frames))
    partial_index = n_written if (truncate_after_steps is not None and truncate_inside_step and n_written < len(frames)) else None
    truncated = truncate_after_steps is not None

    # ---- per-step numbers, formatted once and parsed back for expectations ----
    steps: list[dict[str, Any]] = []
    for k, frame in enumerate(frames):
        frame_cell = np.asarray(frame.get("cell", cell), dtype=float)
        cell_text = [[_f8(x) for x in row] for row in frame_cell]
        cell_back = np.array([[float(x) for x in row] for row in cell_text])
        positions = np.asarray(frame["positions"], dtype=float)
        frac = positions @ np.linalg.inv(frame_cell)
        frac = frac - np.floor(frac)
        frac[np.isclose(frac, 1.0, atol=5e-9)] = 0.0
        frac_text = [[_f8(x) for x in row] for row in frac]
        frac_back = np.array([[float(x) for x in row] for row in frac_text])
        volume = abs(float(np.linalg.det(cell_back)))
        pv = pstress * PV_FACTOR * volume if pstress else 0.0
        F = float(frame["free_energy"])
        e_wo_offset = float(frame.get("e_wo_offset", -0.004))
        e0_offset = float(frame.get("e0_offset", e_wo_offset / 2))
        edisp = float(frame.get("edisp", 0.0))
        label_source = "mlff" if k in mlff else "dft"
        n_scf = 0 if label_source == "mlff" else int(frame.get("n_scf", 12))
        last_dE = float(frame.get("last_dE", ediff / 10))
        forces = np.asarray(frame["forces"], dtype=float)
        forces_text = [[_f8(x) for x in row] for row in forces]
        stress = frame.get("stress_kbar") if isif != 0 and "stress" not in set(frame.get("omit", ())) else None
        stress_text = [[_f8(x) for x in row] for row in np.asarray(stress, dtype=float)] if stress is not None else None
        overflow = set(frame.get("overflow", ()))
        # scstep energies (formatted)
        last_fr = F - edisp
        scf: list[dict[str, str]] = []
        for i in range(n_scf):
            back = n_scf - 1 - i
            delta = 0.0 if back == 0 else min(last_dE * 3.0 ** (back - 1), 50.0)
            e_fr = last_fr + delta
            scf.append({
                "e_fr_energy": _f8(e_fr), "e_wo_entrp": _f8(e_fr + e_wo_offset), "e_0_energy": _f8(e_fr + e0_offset),
            })
        if label_source == "mlff":
            calc = {"e_fr_energy": F + pv, "e_wo_entrp": F + pv, "e_0_energy": F + pv}
        elif bug:
            calc = {"e_fr_energy": F + pv, "e_wo_entrp": F + e0_offset + pv, "e_0_energy": -e_wo_offset}
        else:
            calc = {"e_fr_energy": F + pv, "e_wo_entrp": F + e_wo_offset + pv, "e_0_energy": F + e0_offset + pv}
        if is_md:
            kinetic = float(frame.get("kinetic", 0.1))
            calc.update({"kinetic": kinetic, "lattice kinetic": 0.0, "nosepot": 0.0, "nosekinetic": 0.0,
                         "total": calc["e_fr_energy"] + kinetic})
        calc_text = {name: _f8(value) for name, value in calc.items()}
        calc_back = {name: float(text) for name, text in calc_text.items()}
        f_expected = calc_back["e_fr_energy"] - (pstress * PV_FACTOR * volume if pstress else 0.0)
        expected: dict[str, Any] = {
            "index": k, "label_source": label_source, "complete": k < n_written,
            "cell": cell_back, "fractional": frac_back, "positions": frac_back @ cell_back, "volume": volume,
            "forces": np.array([[float(x) for x in row] for row in forces_text]),
            "stress_kbar": np.array([[float(x) for x in row] for row in stress_text]) if stress_text else None,
            "calc_energies": calc_back, "n_scf": n_scf, "free_energy": f_expected, "pv_term": pv,
            "free_energy_input": F, "edisp": edisp,
        }
        if scf:
            last = {name: float(text) for name, text in scf[-1].items()}
            expected["energy_sigma0"] = f_expected + (last["e_0_energy"] - last["e_fr_energy"])
            expected["energy_no_entropy"] = f_expected + (last["e_wo_entrp"] - last["e_fr_energy"])
            expected["additive_correction"] = calc_back["e_fr_energy"] - pv - last["e_fr_energy"]
            expected["last_dE"] = abs(float(scf[-1]["e_fr_energy"]) - float(scf[-2]["e_fr_energy"])) if n_scf >= 2 else None
        else:
            expected.update(energy_sigma0=None, energy_no_entropy=None, additive_correction=None, last_dE=None)
        steps.append({
            "frame": frame, "cell_text": cell_text, "frac_text": frac_text, "volume": volume, "forces_text": forces_text,
            "stress_text": stress_text, "scf": scf, "calc_text": calc_text, "label_source": label_source,
            "overflow": overflow, "F": F, "E_wo": F + e_wo_offset, "E0": F + e0_offset, "n_scf": n_scf,
            "last_dE": last_dE, "expected": expected, "edisp": edisp,
        })

    # ---- vasprun.xml ----
    parameters = _parameters_xml(
        incar=incar, nelm=nelm, ediff=ediff, ispin=ispin, magmom=magmom_values, ismear=ismear, sigma=sigma,
        encut=encut, nelect=nelect, ibrion=ibrion, nsw=nsw, isif=isif, potim=potim, pstress=pstress,
        lorbit=lorbit, raw=raw_parameters or {},
    )
    xml: list[str] = ['<?xml version="1.0" encoding="ISO-8859-1"?>', f"<!-- {SYNTHETIC_COMMENT} -->", "<modeling>"]
    xml += [" <generator>", '  <i name="program" type="string">vasp </i>',
            f'  <i name="version" type="string">{version}  </i>',
            f'  <i name="subversion" type="string">{SYNTHETIC_SUBVERSION} </i>',
            '  <i name="platform" type="string">SYNTHETIC </i>', '  <i name="date" type="string">2026 01 01 </i>',
            '  <i name="time" type="string">00:00:00 </i>', " </generator>", " <incar>"]
    for name, value in incar.items():
        if name == "MAGMOM":
            xml.append("  " + _xml_v("MAGMOM", magmom_values))
        else:
            xml.append("  " + _xml_i(name, value, None if not isinstance(value, float) else ""))
    xml += [" </incar>", " <kpoints>", '  <generation param="Gamma">',
            '   ' + _xml_v("divisions", [1, 1, 1], "int"), '   ' + _xml_v("usershift", [0, 0, 0]),
            '   ' + _xml_v("genvec1", [1, 0, 0]), '   ' + _xml_v("genvec2", [0, 1, 0]),
            '   ' + _xml_v("genvec3", [0, 0, 1]), '   ' + _xml_v("shift", [0, 0, 0]), "  </generation>",
            '  <varray name="kpointlist" >', f"   <v>{_row([0, 0, 0])} </v>", "  </varray>",
            '  <varray name="weights" >', f"   <v>{_f8(1.0)} </v>", "  </varray>", " </kpoints>"]
    xml += parameters
    xml += [" <atominfo>", f"  <atoms>{n_atoms:8d} </atoms>", f"  <types>{len(blocks):8d} </types>",
            '  <array name="atoms" >', '   <dimension dim="1">ion</dimension>',
            '   <field type="string">element</field>', '   <field type="int">atomtype</field>', "   <set>"]
    type_of = [t for t, (_, count) in enumerate(blocks, 1) for _ in range(count)]
    xml += [f"    <rc><c>{el:2s}</c><c>{t:4d}</c></rc>" for el, t in zip(species, type_of)]
    xml += ["   </set>", "  </array>", '  <array name="atomtypes" >', '   <dimension dim="1">type</dimension>',
            '   <field type="int">atomspertype</field>', '   <field type="string">element</field>',
            "   <field>mass</field>", "   <field>valence</field>", '   <field type="string">pseudopotential</field>',
            "   <set>"]
    xml += [
        f"    <rc><c>{count:4d}</c><c>{el:2s}</c><c>{mass:16.8f}</c><c>{zval:16.8f}</c><c>  {title}                  </c></rc>"
        for (el, count), mass, zval, title in zip(blocks, masses, zvals, titles)
    ]
    xml += ["   </set>", "  </array>", " </atominfo>"]
    first = steps[0]
    xml += _structure_xml(
        "initialpos", first["cell_text"], first["frac_text"], _f8(first["volume"]),
        selective if "initialpos" in selective_in else None, velocities=is_md,
    )
    stop_text: str | None = None
    for k, step in enumerate(steps):
        if k >= n_written and k != partial_index:
            break
        partial = k == partial_index
        body = _step_xml(step, partial=partial, mode=truncate_inside_step)
        if partial:
            stop_text = "\n".join(body)
            break
        xml += body
    if not truncated:
        last = steps[-1]
        final_frac = [[_f8((float(x) + 0.0003) % 1.0) for x in row] for row in last["frac_text"]]
        xml += _structure_xml(
            "finalpos", last["cell_text"], final_frac, _f8(last["volume"]),
            selective if "finalpos" in selective_in else None,
        )
        xml.append("</modeling>")
    text = "\n".join(xml) + "\n"
    if stop_text is not None:
        text += stop_text
    vasprun = _write_maybe_compressed(dirpath / "vasprun.xml", text.encode("latin-1"), compress)

    # ---- inputs ----
    incar_path = dirpath / "INCAR"
    incar_lines = [f"# {SYNTHETIC_COMMENT.replace('output', 'input')}"]
    for name, value in incar.items():
        if isinstance(value, bool):
            value = ".TRUE." if value else ".FALSE."
        incar_lines.append(f"{name} = {value}")
    incar_path.write_text("\n".join(incar_lines) + "\n", encoding="utf-8", newline="\n")
    poscar_path = dirpath / "POSCAR"
    tokens = list(poscar_species) if poscar_species is not None else [el for el, _ in blocks]
    poscar = [SYNTHETIC_COMMENT.replace("output", "input"), "   1.00000000000000"]
    poscar += ["  " + " ".join(f"{float(x):20.14f}" for x in row) for row in first["cell_text"]]
    poscar += ["   " + "   ".join(tokens), "   " + "   ".join(str(count) for _, count in blocks)]
    poscar_flags = selective if isinstance(poscar_selective, str) and poscar_selective == "same" else poscar_selective
    if poscar_flags is not None:
        poscar.append("Selective dynamics")
    poscar.append("Direct")
    flags = np.asarray(poscar_flags, dtype=bool) if poscar_flags is not None else None
    for i, row in enumerate(first["frac_text"]):
        line = "  " + " ".join(f"{float(x):20.16f}" for x in row)
        if flags is not None:
            line += "   " + "   ".join("T" if flag else "F" for flag in flags[i])
        poscar.append(line)
    poscar_path.write_text("\n".join(poscar) + "\n", encoding="utf-8", newline="\n")
    kpoints_path = dirpath / "KPOINTS"
    kpoints_path.write_text(f"{SYNTHETIC_COMMENT.replace('output', 'input')}\n0\nGamma\n 1 1 1\n 0 0 0\n",
                            encoding="utf-8", newline="\n")
    potcar_path = None
    potcar_expected: list[dict[str, Any]] = []
    if write_potcar:
        if potcar_file_order is not None:
            file_elements = list(potcar_file_order)
            file_titles = [f"PAW_PBE {el} SYNTHETIC" for el in file_elements]
        else:
            file_elements, file_titles = [el for el, _ in blocks], titles
        potcar_path, potcar_expected = _write_potcar(dirpath / "POTCAR", file_elements, file_titles, potcar_hash)

    # ---- OUTCAR / OSZICAR ----
    outcar_path = oszicar_path = None
    context = {
        "version": version, "titles": titles, "blocks": blocks, "n_atoms": n_atoms, "nelm": nelm, "ediff": ediff,
        "ibrion": ibrion, "nsw": nsw, "isif": isif, "potim": potim, "ispin": ispin, "lorbit": lorbit,
        "ismear": ismear, "sigma": sigma, "encut": encut, "nelect": nelect, "pstress": pstress, "incar": incar,
        "zvals": zvals, "masses": masses, "mlff": bool(mlff), "outcar_potcar_headers": outcar_potcar_headers,
    }
    dft_written = sum(1 for s in steps[:n_written] if s["label_source"] == "dft")
    if write_outcar:
        outcar_path = dirpath / "OUTCAR"
        outcar_path.write_text(
            _outcar_text(context, steps, n_written, partial_index, truncated, ionic_converged),
            encoding="latin-1", newline="\n",
        )
    if write_oszicar:
        oszicar_path = dirpath / "OSZICAR"
        oszicar_path.write_text(_oszicar_text(context, steps, n_written, partial_index), encoding="latin-1", newline="\n")

    expected_steps = [s["expected"] for s in steps[:n_written]]
    if partial_index is not None:
        partial_expected = dict(steps[partial_index]["expected"], complete=False)
        expected_steps.append(partial_expected)
    return {
        "dir": dirpath, "vasprun": vasprun, "outcar": outcar_path, "oszicar": oszicar_path, "potcar": potcar_path,
        "incar": incar_path, "poscar": poscar_path, "kpoints": kpoints_path,
        "expected": {
            "species": species, "n_atoms": n_atoms, "nelm": nelm, "ediff": ediff, "pstress": pstress,
            "version": version, "ibrion": ibrion, "ispin": ispin, "nelect": nelect, "magmom": magmom_values,
            "titles": titles, "truncated": truncated, "closed": not truncated,
            "n_steps_seen": len(expected_steps), "n_complete": n_written,
            "outcar_dft_steps": dft_written, "ml_steps_outcar": (n_written + (1 if partial_index is not None else 0)) if mlff else 0,
            "steps": expected_steps, "potcar": potcar_expected, "calc_energy_bug": bug,
        },
    }


def _write_maybe_compressed(path: Path, data: bytes, compress: str | None) -> Path:
    if compress is None:
        path.write_bytes(data)
        return path
    target = path.with_name(path.name + "." + compress)
    if compress == "gz":
        with gzip.GzipFile(filename="", mode="wb", fileobj=target.open("wb"), mtime=0) as handle:
            handle.write(data)
    elif compress == "bz2":
        target.write_bytes(bz2.compress(data))
    elif compress == "xz":
        target.write_bytes(lzma.compress(data))
    else:
        raise ValueError(f"unknown compression {compress!r}")
    return target


def _parameters_xml(
    *, incar, nelm, ediff, ispin, magmom, ismear, sigma, encut, nelect, ibrion, nsw, isif, potim, pstress, lorbit, raw,
) -> list[str]:
    ldau = bool(incar.get("LDAU", False))
    out = [" <parameters>", '  <separator name="general" >', '   ' + _xml_i("SYSTEM", "SYNTHETIC", "string"),
           "  </separator>", '  <separator name="electronic" >',
           '   ' + _xml_i("PREC", str(incar.get("PREC", "normal")).lower(), "string"),
           '   ' + _xml_i("ENMAX", encut, ""), '   ' + _xml_i("EDIFF", ediff, ""),
           '   ' + _xml_i("NELECT", nelect, ""), '   ' + _xml_i("NBANDS", 32, "int"),
           '   <separator name="electronic smearing" >', '    ' + _xml_i("ISMEAR", ismear, "int"),
           '    ' + _xml_i("SIGMA", sigma, ""), "   </separator>", '   <separator name="electronic spin" >',
           '    ' + _xml_i("ISPIN", ispin, "int"), '    ' + _xml_i("LNONCOLLINEAR", False, "logical"),
           '    ' + _xml_v("MAGMOM", magmom), '    ' + _xml_i("NUPDOWN", float(incar.get("NUPDOWN", -1.0)), ""),
           '    ' + _xml_i("LSORBIT", False, "logical"), "   </separator>",
           '   <separator name="electronic exchange-correlation" >',
           '    ' + _xml_i("LASPH", bool(incar.get("LASPH", False)), "logical"), "   </separator>",
           '   <separator name="electronic convergence" >', '    ' + _xml_i("NELM", nelm, "int"),
           '    ' + _xml_i("NELMDL", -5, "int"), '    ' + _xml_i("NELMIN", int(incar.get("NELMIN", 2)), "int"),
           "   </separator>", "  </separator>", '  <separator name="grids" >',
           '   ' + _xml_i("GGA", str(incar.get("GGA", "--")), "string"), "  </separator>",
           '  <separator name="ionic" >', '   ' + _xml_i("EDIFFG", float(incar.get("EDIFFG", ediff * 10)), ""),
           '   ' + _xml_i("NSW", nsw, "int"), '   ' + _xml_i("IBRION", ibrion, "int"),
           '   ' + _xml_i("ISIF", isif, "int"), '   ' + _xml_i("POTIM", potim, ""),
           '   ' + _xml_i("PSTRESS", pstress, ""), '   ' + _xml_i("SMASS", float(incar.get("SMASS", -3.0)), ""),
           "  </separator>", '  <separator name="ionic md" >',
           '   ' + _xml_i("TEBEG", float(incar.get("TEBEG", 0.0)), ""),
           '   ' + _xml_i("TEEND", float(incar.get("TEEND", incar.get("TEBEG", 0.0))), ""),
           '   ' + _xml_i("NBLOCK", int(incar.get("NBLOCK", 1)), "int"), "  </separator>",
           '  <separator name="dos" >', '   ' + _xml_i("LORBIT", lorbit, "int"), "  </separator>",
           '  <separator name="linear response parameters" >',
           '   ' + _xml_i("LEPSILON", bool(incar.get("LEPSILON", False)), "logical"),
           '   ' + _xml_i("OMEGAMAX", -1.0, ""), "  </separator>", '  <separator name="response functions" >',
           '   ' + _xml_i("OMEGAMAX", -30.0, ""), "  </separator>"]
    out += ['  ' + _xml_i("LDAU", ldau, "logical")]
    if ldau:
        n_types = len(str(incar.get("LDAUL", "")).split())
        out += ['  ' + _xml_i("LDAUTYPE", int(incar.get("LDAUTYPE", 2)), "int"),
                '  ' + _xml_v("LDAUL", [int(x) for x in str(incar.get("LDAUL", "")).split()] or [0] * n_types, "int"),
                '  ' + _xml_v("LDAUU", [float(x) for x in str(incar.get("LDAUU", "")).split()]),
                '  ' + _xml_v("LDAUJ", [float(x) for x in str(incar.get("LDAUJ", "")).split()])]
    if raw:
        out.append('  <separator name="synthetic extras" >')
        for name, (kind, text) in raw.items():
            type_attr = f' type="{kind}"' if kind else ""
            out.append(f'   <i{type_attr} name="{name}">{text}</i>')
        out.append("  </separator>")
    out.append(" </parameters>")
    return out


def _step_xml(step: dict[str, Any], *, partial: bool, mode: bool | str) -> list[str]:
    frame = step["frame"]
    overflow = step["overflow"]
    frac_text = step["frac_text"]
    if "positions" in overflow:
        frac_text = [list(row) for row in frac_text]
        frac_text[0][0] = "   *************"
    forces_rows = [" ".join(row) for row in step["forces_text"]]
    if "forces" in overflow:
        forces_rows[0] = "   ************* " + " ".join(step["forces_text"][0][1:])
    omit = set(frame.get("omit", ()))
    if step["label_source"] == "mlff":
        out = _structure_xml(None, step["cell_text"], frac_text, _f8(step["volume"]))
        out += [' <varray name="forces" >'] + [f"  <v>{row} </v>" for row in forces_rows]
        if partial:
            return out[:-1] if mode == "mid_forces" else out + [" </varray>"]
        out.append(" </varray>")
        if step["stress_text"]:
            out += [' <varray name="stress" >'] + [f"  <v>{' '.join(r)} </v>" for r in step["stress_text"]] + [" </varray>"]
        out += [" <energy>"] + _energy_items({n: float(t) for n, t in step["calc_text"].items()}, "  ", overflow)
        out += [" </energy>", ' <time name="totalsc">    0.01    0.01</time>']
        return out
    out = [" <calculation>"]
    for i, energies in enumerate(step["scf"]):
        out += ["  <scstep>", '   <time name="dav">    0.01    0.01</time>', "   <energy>"]
        if i in (0, len(step["scf"]) - 1):
            out += [f'    <i name="{name}">{_f8(0.0)} </i>' for name in _SC_COMPONENTS]
        out += [f'    <i name="{name}">{text} </i>' for name, text in energies.items()]
        out += ["   </energy>", "  </scstep>"]
    out += _structure_xml(None, step["cell_text"], frac_text, _f8(step["volume"]), indent="  ")
    if partial and mode == "mid_forces":
        half = max(1, len(forces_rows) // 2)
        out += ['  <varray name="forces" >'] + [f"   <v>{row} </v>" for row in forces_rows[:half]]
        return out + [f"   <v>{forces_rows[half % len(forces_rows)][:20]}"]
    if "forces" not in omit:
        out += ['  <varray name="forces" >'] + [f"   <v>{row} </v>" for row in forces_rows] + ["  </varray>"]
    if step["stress_text"]:
        stress_rows = [" ".join(r) for r in step["stress_text"]]
        if "stress" in overflow:
            stress_rows[0] = "   ************* " + " ".join(step["stress_text"][0][1:])
        out += ['  <varray name="stress" >'] + [f"   <v>{row} </v>" for row in stress_rows] + ["  </varray>"]
    if partial:
        return out
    if "energy" not in omit:
        out += ["  <energy>"] + _energy_items({n: float(t) for n, t in step["calc_text"].items()}, "   ", overflow)
        out += ["  </energy>"]
    out += ['  <time name="totalsc">    0.11    0.12</time>',
            "  <eigenvalues>", "   <array>", '    <dimension dim="1">band</dimension>', "    <set>",
            '     <set comment="spin 1">', "      <r>   -5.0000    1.0000 </r>", "      <r>    1.0000    0.0000 </r>",
            "     </set>", "    </set>", "   </array>", "  </eigenvalues>",
            '  <separator name="orbital magnetization" >', '   <v name="MAGDIPOLE">      0.00000000      0.00000000      0.00000000</v>',
            "  </separator>", "  <dos>", '   <i name="efermi">      1.00000000 </i>', "  </dos>", " </calculation>"]
    return out


def _write_potcar(path: Path, elements, titles, potcar_hash) -> tuple[Path, list[dict[str, Any]]]:
    data = b""
    expected = []
    for element, title in zip(elements, titles):
        zval, pomass, enmax = ELEMENT_DATA.get(element, (1.0, 1.0, 300.0))
        lines = [
            f"  {title}",
            f"   {zval:.13f}",
            " parameters from PSCTR are:",
            "   COPYR  = (c) SYNTHETIC TEST FIXTURE - not a real POTCAR, no licensed content",
            f"   VRHFIN ={element}: SYNTHETIC",
            "   LEXCH  = PE",
            f"   TITEL  = {title}",
            f"   POMASS = {pomass:9.3f}; ZVAL   = {zval:9.3f}    mass and valenz",
            f"   ENMAX  = {enmax:9.3f}; ENMIN  = {enmax * 0.75:9.3f} eV",
            "   SYNTHETIC TEST FIXTURE - fake pseudopotential header, no data follows",
            " End of Dataset",
        ]
        filtered = "\n".join(line for line in lines if not line.strip().startswith(("SHA256", "COPYR"))) + "\n"
        digest = hashlib.sha256(filtered.encode()).hexdigest()
        header = None
        if potcar_hash:
            header = digest if potcar_hash is True else hashlib.sha256(b"tampered").hexdigest()
            lines.insert(3, f"   SHA256 = {header} {title.split()[1]}/POTCAR")
        text = ("\n".join(lines) + "\n").encode()
        expected.append({
            "titel": title, "element": element, "zval": zval, "pomass": pomass, "enmax": enmax,
            "sha256_header": header, "sha256_dataset_bytes": hashlib.sha256(text).hexdigest(),
            "sha256_header_verified": None if not potcar_hash else potcar_hash is True,
        })
        data += text
    path.write_bytes(data)
    return path, expected


def _outcar_text(ctx, steps, n_written, partial_index, truncated, ionic_converged) -> str:
    version = ctx["version"]
    blocks = ctx["blocks"]
    incar = ctx["incar"]
    lines = [
        f" vasp.{version} 01Jan26 {SYNTHETIC_SUBVERSION} complex",
        " executed on             SYNTHETIC date 2026.01.01  00:00:00",
        " running on    1 total cores",
        f" {SYNTHETIC_COMMENT}",
        "",
        " INCAR:",
    ]
    lines += [f" POTCAR:    {title}" for title in ctx["titles"]]
    for (element, _), title, zval, mass in zip(blocks, ctx["titles"], ctx["zvals"], ctx["masses"]):
        lines.append(f" POTCAR:    {title}")
        if ctx["outcar_potcar_headers"]:
            lines += [f"   VRHFIN ={element}: SYNTHETIC", "   LEXCH  = PE", f"   TITEL  = {title}",
                      f"   POMASS = {mass:8.3f}; ZVAL   = {zval:8.3f}    mass and valenz"]
    lines += [
        "", " exchange correlation table for  LEXCH =        8", "",
        " Dimension of arrays:",
        "   k-points           NKPTS =      1   k-points in BZ     NKDIM =      1   number of bands    NBANDS=     32",
        f"   number of dos      NEDOS =    301   number of ions     NIONS = {ctx['n_atoms']:6d}",
        "   ions per type =          " + "".join(f"{count:4d}" for _, count in blocks),
        "",
        " Startparameter for this run:",
        "   NWRITE =      2    write-flag & timer",
        f"   PREC   = {str(incar.get('PREC', 'normal')).lower()[:6]:6s}    normal or accurate (medium, high low for compatibility)",
        "   ISTART =      0    job   : 0-new  1-cont  2-samecut",
        f"   ISPIN  = {ctx['ispin']:6d}    spin polarized calculation?",
        "   LNONCOLLINEAR =      F non collinear calculations",
        "   LSORBIT =      F    spin-orbit coupling",
        "",
        " Electronic Relaxation 1",
        f"   ENCUT  = {ctx['encut']:6.1f} eV  29.40 Ry    5.42 a.u.",
        f"   NELM   = {ctx['nelm']:6d};   NELMIN= {int(incar.get('NELMIN', 2)):2d}; NELMDL= -5     # of ELM steps",
        f"   EDIFF  = {fortran_e(ctx['ediff'], 1)}   stopping-criterion for ELM",
        "",
        " Ionic relaxation",
        f"   EDIFFG = {fortran_e(float(incar.get('EDIFFG', ctx['ediff'] * 10)), 1)}   stopping-criterion for IOM",
        f"   NSW    = {ctx['nsw']:6d}    number of steps for IOM",
        f"   IBRION = {ctx['ibrion']:6d}    ionic relax: 0-MD 1-quasi-New 2-CG",
        f"   ISIF   = {ctx['isif']:6d}    stress and relaxation",
        f"   POTIM  = {ctx['potim']:6.4f}    time-step for ionic-motion",
        f"   TEBEG  = {float(incar.get('TEBEG', 0.0)):6.1f};   TEEND  = {float(incar.get('TEEND', incar.get('TEBEG', 0.0))):6.1f} temperature during run",
        f"   PSTRESS= {ctx['pstress']:6.1f} pullay stress",
        f"   NELECT = {ctx['nelect']:12.4f}    total number of electrons",
        f"   NUPDOWN= {float(incar.get('NUPDOWN', -1.0)):12.4f}    fix difference up-down",
        "",
        " DOS related values:",
        f"   ISMEAR = {ctx['ismear']:5d};   SIGMA  = {ctx['sigma']:6.2f}  broadening in eV -4-tet -1-fermi 0-gaus",
        "",
        " Write flags",
        f"   LORBIT = {ctx['lorbit']:6d}    0 simple, 1 ext, 2 COOP (PROOUT)",
        "",
        " Exchange correlation treatment:",
        f"   GGA     =    {str(incar.get('GGA', '--')):2s}    GGA type",
        "   LEXCH   =     8    internal setting for exchange type",
    ]
    if incar.get("LDAU"):
        lines += [f" LDA+U is selected, type is set to LDAUTYPE = {int(incar.get('LDAUTYPE', 2)):2d}",
                  f" angular momentum for each species LDAUL = {'    '.join(str(incar.get('LDAUL', '')).split())}",
                  f" U (eV)           for each species LDAUU = {'  '.join(str(incar.get('LDAUU', '')).split())}",
                  f" J (eV)           for each species LDAUJ = {'  '.join(str(incar.get('LDAUJ', '')).split())}"]
    if "IVDW" in incar:
        lines += ["", f"   IVDW         = {incar['IVDW']}"]
    ml_tags = [(name, value) for name, value in incar.items() if name.startswith("ML_")]
    if ml_tags:
        lines += ["", " Machine learning (SYNTHETIC echo):"]
        for name, value in ml_tags:
            value = "T" if value is True else "F" if value is False else value
            lines.append(f"   {name} = {value}")
    lines += ["", f" number of electron {ctx['nelect']:15.7f} magnetization {float(ctx['n_atoms']):15.7f}",
              " (initial moment above: SYNTHETIC header value, not a step value)", ""]

    for k, step in enumerate(steps):
        if k >= n_written and k != partial_index:
            break
        partial = k == partial_index
        lines += _outcar_step(ctx, step, k, partial)
        if partial:
            break
    if not truncated:
        last = steps[n_written - 1]
        if last["frame"].get("site_moments") is not None and last["label_source"] == "dft":
            lines += _mag_table(last["frame"]["site_moments"])  # VASP repeats the final table
        if ionic_converged if ionic_converged is not None else ctx["ibrion"] in (1, 2, 3):
            lines += ["", " reached required accuracy - stopping structural energy minimisation"]
        lines += ["", "", " General timing and accounting informations for this job:",
                  " ========================================================", "",
                  "                  Total CPU time used (sec):        1.000", "                         Elapsed time (sec):        1.000"]
    return "\n".join(lines) + "\n"


def _mag_table(moments) -> list[str]:
    out = ["", " magnetization (x)", "", "# of ion       s       p       d       tot", "-" * 42]
    for i, m in enumerate(moments, 1):
        out.append(f"{i:5d}        0.000   0.000 {m:7.3f} {m:8.3f}")
    out += ["-" * 42, f"tot          0.000   0.000 {sum(moments):7.3f} {sum(moments):8.3f}", ""]
    return out


def _outcar_step(ctx, step, k, partial) -> list[str]:
    frame = step["frame"]
    out: list[str] = []
    positions = np.array([[float(x) for x in row] for row in step["frac_text"]]) @ np.array(
        [[float(x) for x in row] for row in step["cell_text"]]
    )
    forces = np.array([[float(x) for x in row] for row in step["forces_text"]])
    if ctx["mlff"]:
        out += ["", "  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)",
                "  Direction    XX          YY          ZZ          XY          YZ          ZX",
                "  " + "-" * 85, "  Total        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000",
                "  in kB        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000", "",
                " POSITION                                       TOTAL-FORCE (eV/Angst) (ML)", _DASHES]
        out += [f" {p[0]:12.5f} {p[1]:12.5f} {p[2]:12.5f}   {f[0]:14.6f} {f[1]:13.6f} {f[2]:13.6f}"
                for p, f in zip(positions, forces * 0.9)]
        out += [_DASHES, "    total drift:                                0.000000      0.000000      0.000000", "",
                "  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)", "  ---------------------------------------------------",
                f"  free  energy ML TOTEN  = {step['F'] + 0.01:18.8f} eV", "",
                f"  ML energy  without entropy= {step['F'] + 0.01:17.8f}  ML energy(sigma->0) = {step['F'] + 0.01:17.8f}", ""]
    if step["label_source"] == "mlff":
        if not partial:
            out += _outcar_md_block(ctx, step)
        return out
    mag = frame.get("mag_total")
    mag_text = "" if ctx["ispin"] == 1 else f" {float(mag if mag is not None else 0.0):15.7f}"
    scf = step["scf"]
    n_iter = len(scf) if not partial else (max(1, len(scf) // 2) if scf else 0)
    previous = None
    for i in range(n_iter):
        e_fr = float(scf[i]["e_fr_energy"])
        d_e = e_fr if previous is None else e_fr - previous
        previous = e_fr
        d_eps = frame.get("d_eps", -d_e / 2) if i == len(scf) - 1 else -d_e / 2
        out += [
            "", f"----------------------------------------- Iteration {k + 1:4d}({i + 1:4d})  ---------------------------------------",
            "", "    POTLOK:  cpu time    0.0100: real time    0.0100", "",
            " eigenvalue-minimisations  :    32",
            f" total energy-change (2. order) :{fortran_e(d_e, 7, leading_zero_negative=True):>14s}  ({fortran_e(d_eps, 7, leading_zero_negative=True)})",
            f" number of electron {ctx['nelect']:15.7f} magnetization{mag_text}",
            f" augmentation part  {1.0:15.7f} magnetization{mag_text}",
            "", "  free energy    TOTEN  = " + f"{e_fr + step['edisp']:18.8f} eV", "",
            f"  energy without entropy = {e_fr + step['edisp'] + (step['E_wo'] - step['F']):17.8f}  energy(sigma->0) = {e_fr + step['edisp'] + (step['E0'] - step['F']):17.8f}",
        ]
        if step["edisp"] and i == n_iter - 1:
            out.append(f"  Edisp (eV)  {step['edisp']:12.5f}")
    if partial:
        return out
    out.append("")
    marker = frame.get("scf_marker", "reached")
    if marker == "reached":
        out.append(MARKER_REACHED)
    elif marker == "not_reached":
        out.append(MARKER_NOT_REACHED)
    elif marker == "hard_stop":
        out.append(MARKER_HARD_STOP)
    if frame.get("site_moments") is not None:
        out += _mag_table(frame["site_moments"])
    if step["stress_text"]:
        s = np.array([[float(x) for x in row] for row in step["stress_text"]])
        voigt = [s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[2, 0]]
        total = [value * step["volume"] / 1602.1766208 for value in voigt]
        out += ["", "  FORCE on cell =-STRESS in cart. coord.  units (eV):",
                "  Direction    XX          YY          ZZ          XY          YZ          ZX",
                "  " + "-" * 86, "  Total   " + "".join(f"{v:12.5f}" for v in total),
                "  in kB   " + "".join(f"{v:12.5f}" for v in voigt),
                f"  external pressure = {sum(voigt[:3]) / 3:10.2f} kB  Pullay stress = {ctx['pstress']:10.2f} kB"]
    cell = np.array([[float(x) for x in row] for row in step["cell_text"]])
    out += ["", " VOLUME and BASIS-vectors are now :", " " + "-" * 77,
            f"  volume of cell : {step['volume']:12.2f}", "      direct lattice vectors                 reciprocal lattice vectors"]
    rec = np.linalg.inv(cell).T
    out += [f"  {row[0]:13.9f}{row[1]:13.9f}{row[2]:13.9f}  {r[0]:13.9f}{r[1]:13.9f}{r[2]:13.9f}" for row, r in zip(cell, rec)]
    if "forces" not in set(frame.get("omit", ())):
        out += ["", " POSITION                                       TOTAL-FORCE (eV/Angst)", _DASHES]
        out += [f" {p[0]:12.5f} {p[1]:12.5f} {p[2]:12.5f}   {f[0]:14.6f} {f[1]:13.6f} {f[2]:13.6f}"
                for p, f in zip(positions, forces)]
        out += [_DASHES, "    total drift:                                0.000001     -0.000002      0.000000"]
    out += ["", " " + "-" * 104, "", "",
            "  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)", "  ---------------------------------------------------",
            f"  free  energy   TOTEN  = {step['F']:18.8f} eV", "",
            f"  energy  without entropy= {step['E_wo']:17.8f}  energy(sigma->0) = {step['E0']:17.8f}"]
    if ctx["pstress"]:
        pv = step["expected"]["pv_term"]
        out.append(f"  enthalpy is  TOTEN    = {step['F'] + pv:18.8f} eV   P V= {pv:14.8f}")
    out += _outcar_md_block(ctx, step)
    return out


def _outcar_md_block(ctx, step) -> list[str]:
    out = [""]
    if ctx["ibrion"] == 0:
        kinetic = float(step["frame"].get("kinetic", 0.1))
        out += ["  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)", "  ---------------------------------------------------",
                f"% ion-electron   TOTEN  = {step['F']:18.6f}  see above",
                f"  kinetic energy EKIN   = {kinetic:18.6f}", "  kin. lattice  EKIN_LAT=           0.000000  (temperature  300.00 K)",
                "  nose potential ES     =           0.000000", "  nose kinetic   EPS    =           0.000000",
                "  ---------------------------------------------------", f"  total energy   ETOTAL = {step['F'] + kinetic:18.6f} eV", ""]
    out += ["     LOOP+:  cpu time    0.1000: real time    0.1000", ""]
    return out


def _oszicar_text(ctx, steps, n_written, partial_index) -> str:
    lines = [f"   {SYNTHETIC_COMMENT}",
             "       N       E                     dE             d eps       ncg     rms          rms(c)"]
    for k, step in enumerate(steps):
        if k >= n_written and k != partial_index:
            break
        partial = k == partial_index
        scf = step["scf"]
        n_iter = len(scf) if not partial else (max(1, len(scf) // 2) if scf else 0)
        previous = None
        for i in range(n_iter):
            e_fr = float(scf[i]["e_fr_energy"])
            d_e = e_fr if previous is None else e_fr - previous
            previous = e_fr
            lines.append(f"DAV: {i + 1:3d}    {fortran_e(e_fr, 12, leading_zero_negative=True)}   "
                         f"{fortran_e(d_e, 5, leading_zero_negative=True)}   {fortran_e(-d_e / 2, 5, leading_zero_negative=True)}    32   0.100E+00")
        if partial:
            break
        frame = step["frame"]
        mag = frame.get("mag_total")
        mag_text = f"  mag={float(mag if mag is not None else 0.0):11.4f}" if ctx["ispin"] == 2 else ""
        if ctx["ibrion"] == 0:
            kinetic = float(frame.get("kinetic", 0.1))
            lines.append(
                f"{k + 1:7d} T= {300.0:6.0f}. E= {fortran_e(step['F'] + kinetic, 8)} F= {fortran_e(step['F'], 8)} "
                f"E0= {fortran_e(step['E0'], 8)}  EK= {fortran_e(kinetic, 5, leading_zero_negative=True)} "
                f"SP= 0.00E+00 SK= 0.00E+00{mag_text}"
            )
        else:
            d_e = step["F"] - (steps[k - 1]["F"] if k else 0.0)
            lines.append(f"{k + 1:4d} F= {fortran_e(step['F'], 8)} E0= {fortran_e(step['E0'], 8)}  d E ={fortran_e(d_e, 6)}{mag_text}")
    return "\n".join(lines) + "\n"
