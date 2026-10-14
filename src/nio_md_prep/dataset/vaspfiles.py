"""Evidence parsers for the files next to ``vasprun.xml``: OUTCAR, OSZICAR, POTCAR, INCAR, POSCAR/CONTCAR.

None of these files is a label source. They supply cross-checks (OUTCAR
energies/forces per ionic step, SCF exit markers, per-ion moments), step counts
and magnetization (OSZICAR), pseudopotential identity (POTCAR fingerprints,
never content), what the user asked for (INCAR) and selective-dynamics flags
(POSCAR). Formats follow ``design/understand/vasp-format-research.md`` S2, S5,
S6, S7 and amendments items 9-11, measured on real VASP 5.2-6.5 output.

All parsers stream or read small files, never raise on malformed content
(problems are returned on the evidence object; each problem string starts with
a ``model.REASONS`` code or ``missing``/``info``/``bad_value`` and ``": "``),
and raise only ``OSError`` when the file cannot be opened.
"""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path
from typing import Any, Iterator

from .errors import DependencyMissingError
from .fsio import sha256_bytes
from .model import IncarEvidence, OszicarEvidence, OutcarEvidence, PoscarEvidence, PotcarEvidence
from .vasprun import HashedSource, parse_version

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover
    raise DependencyMissingError("numpy is required; install nio-md-prep[dataset]") from exc

_LINE_CHUNK = 1 << 20
_FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def _iter_text_lines(path: Path, problems: list[str]) -> tuple[Iterator[str], HashedSource]:
    """Latin-1 decoded lines (no newline) of a plain or compressed text file, hashed in the same pass."""
    source = HashedSource(path)

    def generate() -> Iterator[str]:
        remainder = b""
        try:
            while True:
                try:
                    chunk = source.stream.read(_LINE_CHUNK)
                except (EOFError, OSError, ValueError) as exc:
                    problems.append(f"truncated_frame: {Path(path).name} could not be read to the end ({exc})")
                    break
                if not chunk:
                    break
                lines = (remainder + chunk).split(b"\n")
                remainder = lines.pop()
                for line in lines:
                    yield line.rstrip(b"\r").decode("latin-1")
            if remainder:
                yield remainder.rstrip(b"\r").decode("latin-1")
        finally:
            source.finish()
            source.close()

    return generate(), source


def fortran_float(token: str) -> float:
    """float() that also accepts Fortran output without 'E' (``0.2737684-111``) and 'D' exponents.

    Raises ValueError for anything else (e.g. ``*****`` overflow).
    """
    try:
        return float(token)
    except ValueError:
        text = token.strip().replace("D", "E").replace("d", "e")
        try:
            return float(text)
        except ValueError:
            match = re.fullmatch(r"([-+]?\d*\.\d*)([-+]\d{2,3})", text)
            if match:
                return float(f"{match.group(1)}e{match.group(2)}")
            raise


def _repair_merged_numbers(line: str) -> str:
    """ASE's fixed-width repair: ``-12353.08821-12353.08821`` -> ``-12353.08821 -12353.08821``."""
    return re.sub(r"(\d)-(\d)", r"\1 -\2", line)


def _floats(tokens: list[str], problems: list[str], what: str) -> list[float]:
    values = []
    for token in tokens:
        try:
            value = fortran_float(token)
        except ValueError:
            problems.append(f"non_finite: {what} contains {token!r}")
            value = math.nan
        values.append(value)
    return values


# --------------------------------------------------------------------------
# OUTCAR
# --------------------------------------------------------------------------

_RE_VERSION = re.compile(r"^\s*vasp\.(\S+)")
_RE_POTCAR = re.compile(r"^\s*POTCAR:\s*(.*?)\s*$")
_RE_TITEL = re.compile(r"^\s*TITEL\s*=\s*(.*?)\s*$")
_RE_VRHFIN = re.compile(r"^\s*VRHFIN\s*=\s*(.*?)\s*$")
_RE_POTCAR_LEXCH = re.compile(r"^\s*LEXCH\s*=\s*(\S+)\s*$")  # POTCAR form; the integer echo has trailing text
_RE_NIONS = re.compile(r"\bNIONS\s*=\s*(\d+)")
_RE_IONS_PER_TYPE = re.compile(r"^\s*ions per type\s*=\s*(.*?)\s*$")
_RE_ITERATION = re.compile(r"Iteration\s+(\d+)\s*\(\s*(\d+)\s*\)")
_RE_ABORT = re.compile(r"aborting loop", re.IGNORECASE)
_RE_DE = re.compile(r"total energy-change \(2\. order\)\s*:\s*(\S+)\s*\(\s*(\S+?)\s*\)")
_RE_MAGTOT = re.compile(r"^\s*number of electron\s+(\S+)\s+magnetization\s*(.*)$")
_RE_MAG_TABLE = re.compile(r"^\s*magnetization \(x\)\s*$")
_RE_DASHES = re.compile(r"^\s*-{5,}\s*$")
_RE_STRESS_KB = re.compile(r"^\s*in kB\s+(.*)$")
_RE_EDISP = re.compile(r"^\s*Edisp \(eV\)\s*(\S+)")
_RE_FORCES = re.compile(r"^\s*POSITION\s+TOTAL-FORCE \(eV/Angst\)\s*(\(ML\))?\s*$")
_RE_ML_CELL = re.compile(r"^\s*ML FORCE on cell")
_RE_ML_END = re.compile(r"^\s*ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM")
_RE_ML_E0 = re.compile(r"^\s*ML energy\s+without entropy")
_RE_STEP_END = re.compile(r"^\s*FREE ENERGIE OF THE ION-ELECTRON SYSTEM")
_RE_TOTEN = re.compile(r"^\s*free  energy   TOTEN\s*=\s*(\S+)\s+eV")
_RE_E0 = re.compile(r"^\s*energy  without entropy=\s*(\S+)\s+energy\(sigma->0\)\s*=\s*(\S+)")
_RE_TIMING = re.compile(r"General timing and accounting informations for this job")
_RE_LDAU_SELECTED = re.compile(r"LDA\+U is selected")
_RE_DFTD3 = re.compile(r"^\s*DFTD3\s+(V\S.*?)\s*$")

#: single-valued tags taken from the OUTCAR INCAR echo / "Startparameter" echo (header only; last wins)
ECHO_SCALAR_TAGS = (
    "NELM", "NELMIN", "NELMDL", "EDIFF", "EDIFFG", "IBRION", "NSW", "ISIF", "ISPIN", "LORBIT", "NUPDOWN",
    "IVDW", "ENCUT", "PREC", "GGA", "METAGGA", "ISMEAR", "SIGMA", "PSTRESS", "POTIM", "NWRITE", "LDAU",
    "LDAUTYPE", "LASPH", "LSORBIT", "LNONCOLLINEAR", "LHFCALC", "AEXX", "HFSCREEN", "NELECT", "TEBEG",
    "TEEND", "SMASS", "MDALGO", "NBLOCK", "ISYM", "LREAL", "ALGO", "IALGO", "ISTART", "ICHARG", "ADDGRID",
    "LMAXMIX", "LEPSILON", "LCALCEPS", "LCHIMAG", "NBANDS", "VDW_S6", "VDW_S8", "VDW_SR", "VDW_A1", "VDW_A2",
)
#: multi-valued tags: the rest of the line is kept (whitespace normalized)
ECHO_LIST_TAGS = ("LDAUL", "LDAUU", "LDAUJ", "MAGMOM", "DIPOL")
_RE_ECHO_SCALAR = re.compile(
    r"(?<![A-Za-z0-9_])(" + "|".join(ECHO_SCALAR_TAGS) + r"|ML_[A-Z0-9_]+)\s*=\s*([^\s;]+)"
)
_RE_ECHO_LIST = re.compile(r"(?<![A-Za-z0-9_])(" + "|".join(ECHO_LIST_TAGS) + r")\s*=\s*([^;!#]*)")


#: Default cap on retained per-step site-moment values (``magnetization (x)`` tables). A
#: 200-atom LORBIT>=10 MD of 10000 steps prints 2e6 values (~64 MB as Python floats); tables
#: beyond the cap are counted but not kept (the final table is always kept).
DEFAULT_MAX_SITE_MOMENT_VALUES = 2_000_000


class _OutcarStep:
    def __init__(self):
        self.iterations = 0
        self.marker: str | None = None
        self.energy_change: tuple[float, float] | None = None
        self.magnetization: float | None = None
        self.table: list[float] | None = None
        self.stress: list[float] | None = None
        self.edisp: float | None = None
        self.force_abs_sum: float | None = None
        self.forces: Any = None

    def has_content(self) -> bool:
        return bool(self.iterations or self.marker or self.energy_change is not None
                    or self.force_abs_sum is not None or self.stress is not None)


def scf_marker_kind(line: str) -> str:
    """Classify an OUTCAR ``aborting loop`` line: ``reached`` | ``not_reached`` | ``hard_stop`` | ``other``.

    Real wordings (research S2.1, S10.4): ``aborting loop because EDIFF is reached``
    (VASP 5 prints it even when NELM was hit -- not proof there), ``aborting loop
    EDIFF was not reached (unconverged)`` (VASP 6), ``aborting loop because hard
    stop was set``.
    """
    lowered = line.lower()
    if "because ediff is reached" in lowered:
        return "reached"
    if "ediff was not reached" in lowered:
        return "not_reached"
    if "hard stop" in lowered:
        return "hard_stop"
    return "other"


_marker_kind = scf_marker_kind  # draft-era name


def parse_outcar(
    path: Path, *, keep_forces: bool = False, max_site_moment_values: int | None = DEFAULT_MAX_SITE_MOMENT_VALUES,
) -> OutcarEvidence:
    """Stream an OUTCAR (plain or compressed) into :class:`OutcarEvidence`.

    Per DFT ionic step (a step closes at the ``FREE ENERGIE OF THE ION-ELECTRON
    SYSTEM`` line that is not prefixed by ``ML``; its ``free  energy   TOTEN``
    and ``energy  without entropy=`` lines follow it): SCF exit marker, number
    of ``Iteration`` blocks, last ``total energy-change (2. order)`` pair, last
    ``number of electron ... magnetization`` value, ``magnetization (x)`` table
    (last column; s/p/d[/f] layouts), ``in kB`` stress (6 values, VASP order XX
    YY ZZ XY YZ ZX), ``Edisp``, and sum |F| over the TOTAL-FORCE rows. VASP-MLFF
    prediction blocks (``ML FORCE on cell``, ``TOTAL-FORCE ... (ML)``, ``ML FREE
    ENERGIE``, ``free  energy ML TOTEN``) are counted in ``ml_steps`` and never
    mixed into DFT values. Header facts are read before the first ionic step
    (the header ends at the first ``Iteration``, ``total energy-change``,
    ``aborting loop``, ML or ``FREE ENERGIE`` line). Every per-step list has
    exactly one entry per closed DFT step (NaN/None where the value was not
    printed); a step still open at EOF is never counted and is reported as a
    ``truncated_frame`` problem naming its SCF marker, if any.

    Memory stays O(one step) plus the per-step scalars: force rows are reduced
    to sum |F| unless ``keep_forces``; per-step site-moment tables are retained
    until ``max_site_moment_values`` values are stored (None = no cap), beyond
    which they are counted in a problem but dropped (the final table is kept).
    """
    problems: list[str] = []
    lines, source = _iter_text_lines(Path(path), problems)

    version: str | None = None
    nions: int | None = None
    ions_per_type: list[int] | None = None
    potcar_lines: list[str] = []
    potcar_blocks: list[dict[str, str]] = []
    executed: dict[str, str] = {}
    header = True

    markers: list[bool | None] = []
    marker_kinds: list[str | None] = []
    iteration_counts: list[int] = []
    free_energies: list[float] = []
    energies_wo: list[float] = []
    energies_e0: list[float] = []
    tables: dict[int, list[float]] = {}
    magnetizations: list[float | None] = []
    force_sums: list[float] = []
    energy_changes: list[tuple[float, float] | None] = []
    dispersion: list[float | None] = []
    stresses: list[list[float] | None] = []
    all_forces: list[Any] | None = [] if keep_forces else None
    ml_steps = 0
    ionic_converged = False
    completed = False

    step = _OutcarStep()
    in_scf = False  # an 'Iteration' line was seen since the last step end
    ml_block = False
    awaiting_summary: int | None = None  # index of a closed step whose TOTEN/E0 lines are due
    post_table: list[float] | None = None  # table printed after a step summary (VASP repeats the final one)
    mode: str | None = None  # multi-line block being read
    mode_lines = 0
    block_rows: list[str] = []

    table_values_kept = 0
    tables_dropped = 0

    def close_block() -> None:
        nonlocal post_table, table_values_kept, tables_dropped
        if mode in {"table_rows"}:
            values = []
            for row in block_rows:
                tokens = row.split()
                if not tokens:
                    continue
                values.extend(_floats(tokens[-1:], problems, "OUTCAR magnetization (x) table"))
            if nions is not None and len(values) != nions:
                problems.append(
                    f"bad_shape: OUTCAR magnetization (x) table has {len(values)} rows for {nions} ions"
                )
            if in_scf:
                if max_site_moment_values is not None and table_values_kept + len(values) > max_site_moment_values:
                    tables_dropped += 1
                else:
                    table_values_kept += len(values)
                    step.table = values
            else:
                post_table = values
        elif mode == "forces_rows":
            rows = []
            for row in block_rows:
                tokens = _repair_merged_numbers(row).split()
                if len(tokens) != 6:
                    problems.append(
                        f"bad_shape: OUTCAR TOTAL-FORCE row of DFT step {len(free_energies)} has {len(tokens)} values"
                    )
                    rows.append([math.nan] * 3)
                    continue
                rows.append(_floats(tokens[3:], problems, f"OUTCAR TOTAL-FORCE of DFT step {len(free_energies)}"))
            array = np.array(rows, dtype=np.float64).reshape(len(rows), 3)
            if nions is not None and len(rows) != nions:
                problems.append(f"bad_shape: OUTCAR DFT step {len(free_energies)} has {len(rows)} force rows for {nions} ions")
            step.force_abs_sum = float(np.abs(array).sum())
            if keep_forces:
                step.forces = array

    def end_header() -> None:
        nonlocal header
        header = False

    unattributed_edisp = 0
    for line in lines:
        if mode is not None:
            mode_lines += 1
            if mode.endswith("_wait"):
                if _RE_DASHES.match(line):
                    mode = mode[: -len("_wait")] + "_rows"
                    block_rows = []
                    mode_lines = 0
                elif mode_lines > 6:
                    problems.append(f"bad_shape: OUTCAR {mode[:-5]} header without its dashed separator; block ignored")
                    mode = None
                continue
            if _RE_DASHES.match(line):
                if mode != "mlforces_rows":
                    close_block()
                mode = None
                block_rows = []
                continue
            block_rows.append(line)
            if nions is not None and len(block_rows) > nions + 2:
                problems.append(f"bad_shape: OUTCAR {mode[:-5]} block longer than NIONS={nions}; block ignored")
                mode = None
                block_rows = []
            continue

        if "Iteration" in line and _RE_ITERATION.search(line):
            if header:
                end_header()
            if not in_scf:
                in_scf = True
                ml_block = False
                post_table = None
                awaiting_summary = None
            step.iterations += 1
            continue

        if header:
            if version is None and (match := _RE_VERSION.match(line)):
                version = match.group(1)
                continue
            if match := _RE_POTCAR.match(line):
                potcar_lines.append(match.group(1))
                potcar_blocks.append({"potcar_line": match.group(1)})
                continue
            if potcar_blocks and (match := _RE_TITEL.match(line)):
                potcar_blocks[-1].setdefault("titel", match.group(1))
                continue
            if potcar_blocks and (match := _RE_VRHFIN.match(line)):
                potcar_blocks[-1].setdefault("vrhfin", match.group(1))
                continue
            if potcar_blocks and (match := _RE_POTCAR_LEXCH.match(line)):
                potcar_blocks[-1].setdefault("lexch", match.group(1))
                continue
            if match := _RE_NIONS.search(line):
                nions = int(match.group(1))
            if match := _RE_IONS_PER_TYPE.match(line):
                try:
                    ions_per_type = [int(token) for token in match.group(1).split()]
                except ValueError:
                    problems.append(f"bad_value: OUTCAR 'ions per type' line unreadable: {match.group(1)!r}")
                continue
            if _RE_LDAU_SELECTED.search(line):
                executed.setdefault("LDAU", "T")
            if match := _RE_DFTD3.match(line):
                executed["DFTD3"] = match.group(1)
            if "=" in line:
                for match in _RE_ECHO_LIST.finditer(line):
                    executed[match.group(1)] = " ".join(match.group(2).split())
                for match in _RE_ECHO_SCALAR.finditer(line):
                    executed[match.group(1)] = match.group(2)
            if not (
                _RE_ML_CELL.match(line) or _RE_ML_END.match(line) or _RE_STEP_END.match(line)
                or "aborting loop" in line or ("energy-change" in line and _RE_DE.search(line))
            ):
                continue
            end_header()

        # ---- VASP-MLFF prediction blocks: counted, never mixed into DFT values ----
        if "ML" in line:
            if _RE_ML_CELL.match(line):
                ml_block = True
                continue
            if _RE_ML_END.match(line):
                ml_steps += 1
                ml_block = True
                continue
            if _RE_ML_E0.match(line):
                ml_block = False
                continue
        if "POSITION" in line and (match := _RE_FORCES.match(line)):
            mode = "mlforces_wait" if match.group(1) else "forces_wait"
            mode_lines = 0
            continue
        if "aborting loop" in line:
            if not ml_block:
                step.marker = _marker_kind(line)
            continue
        if "energy-change" in line and (match := _RE_DE.search(line)):
            pair = _floats([match.group(1), match.group(2)], problems, "OUTCAR energy change")
            step.energy_change = (pair[0], pair[1])
            continue
        if "number of electron" in line and (match := _RE_MAGTOT.match(line)):
            if in_scf:
                tokens = match.group(2).split()
                if len(tokens) == 1:
                    step.magnetization = _floats(tokens, problems, "OUTCAR total magnetization")[0]
                else:  # ISPIN=1 prints no value; noncollinear prints a vector
                    step.magnetization = None
            continue
        if "magnetization (x)" in line and _RE_MAG_TABLE.match(line):
            mode = "table_wait"
            mode_lines = 0
            continue
        if "in kB" in line and (match := _RE_STRESS_KB.match(line)):
            if not ml_block:
                tokens = _repair_merged_numbers(match.group(1)).split()
                if len(tokens) == 6:
                    step.stress = _floats(tokens, problems, "OUTCAR 'in kB' stress")
                else:
                    problems.append(f"bad_shape: OUTCAR 'in kB' line has {len(tokens)} values")
            continue
        if "Edisp" in line and (match := _RE_EDISP.match(line)):
            # Attributed only inside a step's SCF chunk or its summary; where VASP prints
            # it relative to the Iteration blocks has not been verified on real IVDW output.
            value = _floats([match.group(1)], problems, "OUTCAR Edisp")[0]
            if in_scf:
                step.edisp = value
            elif awaiting_summary is not None and dispersion[awaiting_summary] is None:
                dispersion[awaiting_summary] = value
            else:
                unattributed_edisp += 1
            continue
        if "FREE ENERGIE" in line and _RE_STEP_END.match(line):
            index = len(free_energies)
            if not in_scf:
                problems.append(f"info: OUTCAR DFT step {index} summary without preceding SCF iterations")
            markers.append({"reached": True, "not_reached": False, "hard_stop": False}.get(step.marker or "", None))
            marker_kinds.append(step.marker)
            iteration_counts.append(step.iterations)
            free_energies.append(math.nan)
            energies_wo.append(math.nan)
            energies_e0.append(math.nan)
            if step.table is not None:
                tables[index] = step.table
            magnetizations.append(step.magnetization)
            if step.force_abs_sum is None:
                problems.append(f"missing_forces: OUTCAR DFT step {index} has no TOTAL-FORCE block")
            force_sums.append(step.force_abs_sum if step.force_abs_sum is not None else math.nan)
            energy_changes.append(step.energy_change)
            dispersion.append(step.edisp)
            stresses.append(step.stress)
            if all_forces is not None:
                all_forces.append(step.forces)
            step = _OutcarStep()
            in_scf = False
            ml_block = False
            awaiting_summary = index
            continue
        if "TOTEN" in line and (match := _RE_TOTEN.match(line)):
            if awaiting_summary is not None and math.isnan(free_energies[awaiting_summary]):
                free_energies[awaiting_summary] = _floats([match.group(1)], problems, "OUTCAR TOTEN")[0]
            else:
                problems.append("info: OUTCAR 'free  energy   TOTEN' line outside an ionic-step summary ignored")
            continue
        if "without entropy" in line and (match := _RE_E0.match(line)):
            if awaiting_summary is not None:
                values = _floats([match.group(1), match.group(2)], problems, "OUTCAR energy without entropy")
                energies_wo[awaiting_summary], energies_e0[awaiting_summary] = values
                awaiting_summary = None
            continue
        if "reached required accuracy" in line:
            ionic_converged = True
            continue
        if "General timing" in line and _RE_TIMING.search(line):
            completed = True
            continue

    if mode is not None:
        problems.append("truncated_frame: OUTCAR ends inside a table/force block")
    if in_scf or step.has_content():
        problems.append(
            f"truncated_frame: OUTCAR ends inside DFT ionic step {len(free_energies)} (not counted; "
            f"{step.iterations} Iteration lines, SCF marker {step.marker!r})"
        )
    if tables_dropped:
        problems.append(
            f"info: {tables_dropped} per-step magnetization (x) tables not retained "
            f"(max_site_moment_values={max_site_moment_values})"
        )
    for index, value in enumerate(free_energies):
        if math.isnan(value) or math.isnan(energies_e0[index]):
            problems.append(f"missing_energy: OUTCAR DFT step {index} summary energies not printed")
    if unattributed_edisp:
        problems.append(f"info: {unattributed_edisp} 'Edisp (eV)' lines outside ionic-step chunks were not attributed")

    titles, title_problem = _potcar_titles(potcar_lines)
    if title_problem:
        problems.append(title_problem)
    if ions_per_type is not None and titles and len(titles) != len(ions_per_type):
        problems.append(
            f"species_mismatch: OUTCAR lists {len(titles)} POTCAR titles but {len(ions_per_type)} ion types"
        )
    if ions_per_type is not None and nions is not None and sum(ions_per_type) != nions:
        problems.append(f"species_mismatch: OUTCAR ions per type sum {sum(ions_per_type)} != NIONS {nions}")
    headers = [
        {key: block[key] for key in ("titel", "vrhfin", "lexch") if key in block}
        for block in potcar_blocks
        if "titel" in block
    ]
    if headers and len(headers) != len(titles):
        problems.append(f"info: {len(headers)} POTCAR header echoes for {len(titles)} POTCAR titles")

    return OutcarEvidence(
        sha256=source.sha256 or "",
        version=version,
        nions=nions,
        ions_per_type=ions_per_type,
        potcar_titles=titles,
        scf_converged_markers=markers,
        free_energies=free_energies,
        energies_no_entropy=energies_wo,
        energies_sigma0=energies_e0,
        magnetization_tables=tables,
        total_magnetizations=magnetizations,
        ionic_converged=ionic_converged,
        completed=completed,
        problems=problems,
        version_tuple=parse_version(version),
        executed_tags=executed,
        force_abs_sums=force_sums,
        last_energy_changes=energy_changes,
        dispersion_energies=dispersion,
        ml_steps=ml_steps,
        scf_marker_kinds=marker_kinds,
        scf_iteration_counts=iteration_counts,
        stress_kbar=stresses,
        potcar_headers=headers,
        final_magnetization_table=post_table,
        forces=all_forces,
    )


#: VASP major version from which the OUTCAR exit line distinguishes converged from
#: unconverged loops (VASP 5 prints "EDIFF is reached" even when NELM was hit; research S2.1).
MARKER_PROOF_MIN_VERSION = (6, 0, 0)


def outcar_scf_markers(outcar: OutcarEvidence, n_dft_steps: int) -> dict[str, Any]:
    """Map OUTCAR ``aborting loop`` markers onto vasprun DFT steps -- only when the counts match.

    ``mapped`` is True only if the OUTCAR closed exactly ``n_dft_steps`` DFT
    steps AND every one of them printed a marker (NWRITE<=1 prints it for the
    first ionic step only; a truncated OUTCAR has fewer steps). Otherwise the
    markers are not attributed to steps at all (``kinds`` None, ``reason`` says
    why). ``marker_is_proof`` is True only for VASP >= 6, where the line
    distinguishes converged from unconverged loops; for VASP 5 and unknown
    versions a ``reached`` marker is informational, never proof.
    """
    kinds = list(outcar.scf_marker_kinds)
    n_markers = sum(kind is not None for kind in kinds)
    version = outcar.version_tuple
    proof = version is not None and tuple(version) >= MARKER_PROOF_MIN_VERSION
    if len(kinds) != n_dft_steps:
        reason = f"step_count_mismatch: OUTCAR closed {len(kinds)} DFT steps, vasprun has {n_dft_steps}"
    elif n_markers != n_dft_steps:
        reason = f"marker_count_mismatch: {n_markers} markers for {n_dft_steps} DFT steps"
    else:
        reason = None
    return {
        "mapped": reason is None,
        "reason": reason,
        "kinds": kinds if reason is None else None,
        "n_outcar_steps": len(kinds),
        "n_markers": n_markers,
        "marker_is_proof": proof,
        "outcar_version": outcar.version,
    }


def oszicar_step_values(oszicar: OszicarEvidence, n_steps: int) -> dict[str, Any]:
    """Per-step OSZICAR values (``mag``, ``T``, SCF counts) -- only when the ionic-line count matches.

    ``n_steps`` counts ALL ionic steps (DFT and MLFF; OSZICAR writes a line for
    force-field-only steps too). Values are never labels (8 significant digits).
    """
    lines = oszicar.ionic_lines
    if len(lines) != n_steps:
        return {"mapped": False, "reason": f"step_count_mismatch: OSZICAR has {len(lines)} ionic lines, "
                                           f"expected {n_steps}", "mag": None, "temperature_K": None,
                "scf_steps": None, "n_lines": len(lines)}
    mags = [line.get("mag") for line in lines]
    return {
        "mapped": True, "reason": None,
        "mag": mags if any(value is not None for value in mags) else None,
        "temperature_K": [line.get("T") for line in lines] if any("T" in line for line in lines) else None,
        "scf_steps": list(oszicar.scf_steps), "n_lines": len(lines),
    }


def _potcar_titles(lines: list[str]) -> tuple[list[str], str | None]:
    """VASP echoes ``POTCAR:`` 2 x NTYP times; the first half is the species order."""
    if not lines:
        return [], "missing: OUTCAR has no 'POTCAR:' lines"
    half, odd = divmod(len(lines), 2)
    if odd:
        return list(lines), f"info: OUTCAR has an odd number ({len(lines)}) of 'POTCAR:' lines; all kept"
    if lines[:half] != lines[half:]:
        return lines[:half], "info: the two halves of the OUTCAR 'POTCAR:' lines differ; first half kept"
    return lines[:half], None


# --------------------------------------------------------------------------
# OSZICAR
# --------------------------------------------------------------------------

_OSZ_RELAX = re.compile(r"^\s*(\d+)\s+F=\s*(\S+)\s+E0=\s*(\S+)\s+d E\s*=\s*(\S+)(?:\s+mag=\s*(.+))?")
_OSZ_MD = re.compile(
    r"^\s*(\d+)\s+T=\s*(\S+)\s+E=\s*(\S+)\s+F=\s*(\S+)\s+E0=\s*(\S+)\s+EK=\s*(\S+)\s+SP=\s*(\S+)"
    r"\s+SK=\s*(\S+)(?:\s+mag=\s*(.+))?"
)
_OSZ_SCF = re.compile(r"^\s*([A-Z][A-Za-z]{1,3})\s?:\s*(\d+)\s")


def parse_oszicar(path: Path) -> OszicarEvidence:
    """OSZICAR ionic lines (relaxation/static and MD forms) and SCF iterations per ionic step.

    Only for step counts, temperature, ``mag=`` and SCF counts: the ionic
    lines carry 8 significant digits and are never labels. ``scf_steps[k]`` is
    the number of electronic lines (DAV/RMM/CG/...) before ionic line k (0 for
    VASP-MLFF force-field-only steps).
    """
    problems: list[str] = []
    lines, source = _iter_text_lines(Path(path), problems)
    ionic: list[dict[str, float]] = []
    kinds: list[str] = []
    scf_steps: list[int] = []
    pending_scf = 0
    for line in lines:
        if _OSZ_SCF.match(line):
            pending_scf += 1
            continue
        match = _OSZ_MD.match(line)
        if match:
            keys = ("T", "E", "F", "E0", "EK", "SP", "SK")
            values = _floats(list(match.groups()[1:8]), problems, f"OSZICAR MD line {match.group(1)}")
            record: dict[str, Any] = {"step": int(match.group(1)), **dict(zip(keys, values))}
            mag = match.group(9)
            kinds.append("md")
        else:
            match = _OSZ_RELAX.match(line)
            if not match:
                continue
            keys = ("F", "E0", "dE")
            values = _floats(list(match.groups()[1:4]), problems, f"OSZICAR ionic line {match.group(1)}")
            record = {"step": int(match.group(1)), **dict(zip(keys, values))}
            mag = match.group(5)
            kinds.append("relax")
        if mag is not None:
            tokens = mag.split()
            mags = _floats(tokens, problems, f"OSZICAR mag= of step {record['step']}")
            if len(mags) == 1:
                record["mag"] = mags[0]
            else:
                for axis, value in zip("xyz", mags):
                    record[f"mag_{axis}"] = value
        ionic.append(record)
        scf_steps.append(pending_scf)
        pending_scf = 0
    if pending_scf:
        problems.append(f"truncated_frame: OSZICAR ends with {pending_scf} SCF lines after the last ionic line")
    numbers = [int(record["step"]) for record in ionic]
    if numbers != list(range(1, len(numbers) + 1)):
        problems.append("info: OSZICAR ionic step numbers are not 1..N in order")
    return OszicarEvidence(
        sha256=source.sha256 or "", ionic_lines=ionic, scf_steps=scf_steps, problems=problems, line_kinds=kinds,
    )


# --------------------------------------------------------------------------
# POTCAR (fingerprints only; the licensed content is never stored)
# --------------------------------------------------------------------------

_POT_KEYS = {
    "titel": re.compile(rb"^\s*TITEL\s*=\s*(.*?)\s*$"),
    "vrhfin": re.compile(rb"^\s*VRHFIN\s*=\s*(.*?)\s*$"),
    "lexch": re.compile(rb"^\s*LEXCH\s*=\s*(\S+)"),
    "pomass": re.compile(rb"POMASS\s*=\s*([^;\s]+)"),
    "zval": re.compile(rb"ZVAL\s*=\s*([^;\s]+)"),
    "enmax": re.compile(rb"ENMAX\s*=\s*([^;\s]+)"),
    "enmin": re.compile(rb"ENMIN\s*=\s*([^;\s]+)"),
    "sha256_header": re.compile(rb"^\s*SHA256\s*=\s*([0-9A-Fa-f]{64})"),
}
_POT_NUMERIC = ("pomass", "zval", "enmax", "enmin")
_END_OF_DATASET = b"End of Dataset"


def _hash_without_sha_copyr(lines: list[bytes]) -> str:
    kept = [line for line in lines if not line.strip().startswith((b"SHA256", b"COPYR"))]
    return hashlib.sha256(b"\n".join(kept)).hexdigest()


def parse_potcar(path: Path) -> PotcarEvidence:
    """Per-dataset POTCAR fingerprints (amendments item 11); the content is never kept.

    Datasets end at a line ``End of Dataset``. For each: ``symbol`` (TITEL's
    second token, e.g. ``Ni_pv``), ``element``, ``titel``, ``vrhfin``,
    ``lexch``, ``zval``, ``pomass``, ``enmax``, ``enmin``, ``sha256_header``
    (the ``SHA256 =`` line of hashed releases, else None),
    ``sha256_header_verified`` (None without a header line) and
    ``sha256_verify_variant`` (which reconstruction matched: ``"exact"`` = the
    dataset's own lines, or ``"pymatgen_strip"`` = pymatgen's
    ``f"{chunk.strip()}\\nEnd of Dataset\\n"``), and ``sha256_dataset_bytes``
    (the dataset's exact bytes). ``sha256`` is the whole file's digest.
    """
    path = Path(path)
    data = path.read_bytes()
    problems: list[str] = []
    datasets: list[dict[str, Any]] = []
    raw_lines = data.split(b"\n")
    start_line = 0
    offset = 0
    pieces = data.split(_END_OF_DATASET)  # pymatgen's splitting, for the alternative reconstruction
    for number, line in enumerate(raw_lines):
        line_end = offset + len(line) + 1
        if line.strip() == _END_OF_DATASET:
            chunk_lines = [raw.rstrip(b"\r") for raw in raw_lines[start_line : number + 1]]
            dataset_bytes = data[_line_offset(raw_lines, start_line) : min(line_end, len(data))]
            datasets.append(_potcar_dataset(chunk_lines, dataset_bytes, pieces, len(datasets), problems))
            start_line = number + 1
        offset = line_end
    trailing = b"\n".join(raw_lines[start_line:]).strip()
    if trailing:
        problems.append(
            "missing: POTCAR has content after the last 'End of Dataset'"
            if datasets else "missing: POTCAR has no 'End of Dataset' line"
        )
    if not datasets:
        problems.append("missing: no POTCAR datasets found")
    return PotcarEvidence(sha256=sha256_bytes(data), datasets=datasets, problems=problems)


def _line_offset(lines: list[bytes], index: int) -> int:
    return sum(len(line) + 1 for line in lines[:index])


def _potcar_dataset(
    lines: list[bytes], dataset_bytes: bytes, pieces: list[bytes], position: int, problems: list[str]
) -> dict[str, Any]:
    found: dict[str, Any] = {}
    for line in lines:
        for key, pattern in _POT_KEYS.items():
            if key not in found and (match := pattern.search(line)):
                found[key] = match.group(1).decode("latin-1").strip()
    record: dict[str, Any] = {key: found.get(key) for key in ("titel", "vrhfin", "lexch")}
    for key in _POT_NUMERIC:
        text = found.get(key)
        try:
            record[key] = float(text) if text is not None else None
        except ValueError:
            record[key] = None
            problems.append(f"bad_value: POTCAR dataset {position} {key.upper()} is {text!r}")
    titel = record["titel"]
    symbol = titel.split()[1] if titel and len(titel.split()) >= 2 else None
    record["symbol"] = symbol
    element = symbol.split("/")[0].split("_")[0] if symbol else None
    vrhfin_element = None
    if record["vrhfin"]:
        match = re.match(r"([A-Z][a-z]?)", record["vrhfin"])
        vrhfin_element = match.group(1) if match else None
    if element and vrhfin_element and element != vrhfin_element:
        problems.append(
            f"species_mismatch: POTCAR dataset {position} TITEL element {element} != VRHFIN element {vrhfin_element}"
        )
    record["element"] = element or vrhfin_element
    if titel is None:
        problems.append(f"missing: POTCAR dataset {position} has no TITEL")
    header = found.get("sha256_header")
    record["sha256_header"] = header.lower() if header else None
    record["sha256_dataset_bytes"] = hashlib.sha256(dataset_bytes).hexdigest()
    verified: bool | None = None
    variant: str | None = None
    if header:
        candidates = {"exact": _hash_without_sha_copyr(lines + ([b""] if dataset_bytes.endswith(b"\n") else []))}
        if position < len(pieces):
            text = pieces[position].strip()
            stripped = (text + b"\n" + _END_OF_DATASET + b"\n").replace(b"\r\n", b"\n")
            candidates["pymatgen_strip"] = _hash_without_sha_copyr(stripped.split(b"\n"))
        verified = False
        for name, digest in candidates.items():
            if digest == header.lower():
                verified, variant = True, name
                break
    record["sha256_header_verified"] = verified
    record["sha256_verify_variant"] = variant
    return record


# --------------------------------------------------------------------------
# INCAR
# --------------------------------------------------------------------------

# Adapted from lgutsev/InterfaceForge@f9aa4a2 src/interfaceforge/derivative_probe.py
# _incar_assignments (lines 157-171) (MIT). Changes: undecodable bytes are replaced and
# recorded instead of raising; unparseable assignments are recorded as problems instead of
# raising; duplicate tags are recorded (VASP uses the last one, as does this parser).
_INCAR_ASSIGNMENT = re.compile(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)\s*")


def parse_incar(path: Path) -> IncarEvidence:
    """INCAR -> upper-case tag -> raw value string.

    Handles ``;``-separated assignments, ``!``/``#`` comments and
    backslash-continued lines. Tags are case-insensitive (stored upper case);
    values are kept verbatim (stripped) for :func:`expand_vasp_list`,
    :func:`vasp_bool` and :func:`vasp_number`.
    """
    raw = Path(path).read_bytes()
    problems: list[str] = []
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("utf-8", errors="replace")
        problems.append("info: INCAR is not valid UTF-8; undecodable bytes were replaced")
    physical = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    logical: list[tuple[int, str]] = []  # (first physical line number, joined text)
    pending: tuple[int, str] | None = None
    for number, line in enumerate(physical, 1):
        if pending is not None:
            number, line = pending[0], pending[1] + " " + line
            pending = None
        if line.endswith("\\"):
            pending = (number, line[:-1])
            continue
        logical.append((number, line))
    if pending is not None:
        logical.append(pending)
    tags: dict[str, str] = {}
    duplicates: list[str] = []
    for number, line in logical:
        active = re.split(r"[!#]", line, maxsplit=1)[0]
        for assignment in active.split(";"):
            if not assignment.strip():
                continue
            match = _INCAR_ASSIGNMENT.fullmatch(assignment)
            if not match:
                problems.append(f"bad_value: INCAR line {number}: cannot parse {assignment.strip()!r}")
                continue
            tag = match.group(1).upper()
            if tag in tags and tag not in duplicates:
                duplicates.append(tag)
            tags[tag] = match.group(2)
    return IncarEvidence(sha256=sha256_bytes(raw), tags=tags, duplicates=duplicates, problems=problems)


# Adapted from lgutsev/InterfaceForge@4501e34 src/interfaceforge/vasp.py _vasp_list_length
# (lines 435-440) (MIT). Changed to return the expanded tokens instead of their count.
def expand_vasp_list(value: str) -> list[str]:
    """Expand VASP ``N*value`` shorthand: ``"2*1.0 -2"`` -> ``["1.0", "1.0", "-2"]``."""
    tokens: list[str] = []
    for token in str(value).split():
        repeat = re.fullmatch(r"(\d+)\*(.+)", token)
        if repeat:
            tokens.extend([repeat.group(2)] * int(repeat.group(1)))
        else:
            tokens.append(token)
    return tokens


def vasp_floats(value: str) -> list[float]:
    """``expand_vasp_list`` as floats (ValueError on a non-number)."""
    return [fortran_float(token) for token in expand_vasp_list(value)]


def vasp_bool(value: str | None) -> bool | None:
    """``.TRUE.``/``T``/``true`` -> True, ``.FALSE.``/``F`` -> False, else None."""
    if value is None:
        return None
    token = str(value).strip().split()[0].strip(".").upper() if str(value).strip() else ""
    if token in {"T", "TRUE"}:
        return True
    if token in {"F", "FALSE"}:
        return False
    return None


def vasp_number(value: str | None) -> float | None:
    """First token as a float, or None."""
    if value is None or not str(value).split():
        return None
    try:
        return fortran_float(str(value).split()[0])
    except ValueError:
        return None


# --------------------------------------------------------------------------
# POSCAR / CONTCAR
# --------------------------------------------------------------------------

_ELEMENT = re.compile(r"[A-Z][a-z]?")


def normalize_species_token(token: str) -> str:
    """``Na_pv/6a2f546d`` -> ``Na``; ``Fe_pv`` -> ``Fe`` (VASP >= 6.4 CONTCAR/XDATCAR tokens)."""
    return token.split("/")[0].split("_")[0]


def parse_poscar(path: Path) -> PoscarEvidence:
    """POSCAR/CONTCAR: species line (VASP 5; ``Na_pv/6a2f546d`` tokens), counts,
    optional ``Selective dynamics`` flags (DIRECT basis, True = may move) and
    Direct/Cartesian coordinates. Velocity blocks are ignored.

    A single positive scale multiplies the lattice (and Cartesian positions); a
    negative scale is the cell volume. Three scale values are recorded as a
    problem and leave cell/positions None (species, counts and flags are still read).
    """
    raw = Path(path).read_bytes()
    problems: list[str] = []
    lines = raw.decode("latin-1").replace("\r\n", "\n").replace("\r", "\n").split("\n")

    def fail(message: str) -> PoscarEvidence:
        problems.append(message)
        return PoscarEvidence(
            sha256=sha256_bytes(raw), comment=lines[0] if lines else "", cell=None, species_tokens=None,
            species=None, counts=[], coordinate_mode="direct", fractional=None, positions=None, problems=problems,
        )

    if len(lines) < 8:
        return fail("bad_shape: POSCAR is too short")
    comment = lines[0].strip()
    try:
        scale = [fortran_float(token) for token in lines[1].split()]
        lattice = np.array([[fortran_float(t) for t in lines[i].split()[:3]] for i in (2, 3, 4)], dtype=np.float64)
    except ValueError:
        return fail("bad_value: POSCAR scale or lattice lines are not numbers")
    if lattice.shape != (3, 3):
        return fail("bad_shape: POSCAR lattice is not 3 x 3")
    cursor = 5
    tokens = lines[cursor].split()
    species_tokens: list[str] | None = None
    if tokens and not all(token.isdigit() for token in tokens):
        species_tokens = tokens
        cursor += 1
    try:
        counts = [int(token) for token in lines[cursor].split()]
    except ValueError:
        return fail(f"bad_value: POSCAR ion-count line unreadable: {lines[cursor]!r}")
    if not counts:
        return fail("bad_value: POSCAR has no ion counts")
    cursor += 1
    species = None
    if species_tokens is not None:
        species = [normalize_species_token(token) for token in species_tokens]
        bad = [token for token in species if not _ELEMENT.fullmatch(token)]
        if bad:
            problems.append(f"bad_value: POSCAR species tokens are not element symbols: {bad}")
        if len(species) != len(counts):
            problems.append(f"species_mismatch: POSCAR has {len(species)} species but {len(counts)} counts")
    else:
        problems.append("missing: POSCAR has no species line (VASP 4 format)")
    selective_enabled = cursor < len(lines) and lines[cursor].strip()[:1].lower() == "s"
    if selective_enabled:
        cursor += 1
    if cursor >= len(lines):
        return fail("bad_shape: POSCAR ends before the coordinate-mode line")
    mode_char = lines[cursor].strip()[:1].lower()
    coordinate_mode = "cartesian" if mode_char in {"c", "k"} else "direct"
    cursor += 1
    nions = sum(counts)
    rows = lines[cursor : cursor + nions]
    if len(rows) != nions or any(not row.split() for row in rows):
        return fail(f"bad_shape: POSCAR has fewer than {nions} coordinate rows")
    coords = []
    flags = []
    for number, row in enumerate(rows, 1):
        fields = row.split()
        try:
            coords.append([fortran_float(token) for token in fields[:3]])
        except ValueError:
            coords.append([math.nan] * 3)
            problems.append(f"non_finite: POSCAR coordinate row {number} is not numeric")
        if selective_enabled:
            values = [vasp_bool(token) for token in fields[3:6]]
            if len(values) != 3 or any(value is None for value in values):
                problems.append(f"bad_shape: POSCAR selective-dynamics row {number} lacks T/F flags")
                values = [None, None, None]
            flags.append(values)
    coords_array = np.array(coords, dtype=np.float64).reshape(nions, 3)
    selective = None
    if selective_enabled and not any(value is None for row in flags for value in row):
        selective = np.array(flags, dtype=bool).reshape(nions, 3)
    cell = fractional = positions = None
    if len(scale) == 1:
        factor = scale[0]
        if factor < 0:
            volume = abs(float(np.linalg.det(lattice)))
            factor = (-factor / volume) ** (1.0 / 3.0) if volume > 0 else math.nan
        cell = lattice * factor
        if coordinate_mode == "cartesian":
            positions = coords_array * factor
            try:
                fractional = positions @ np.linalg.inv(cell)
            except np.linalg.LinAlgError:
                problems.append("degenerate_cell: POSCAR lattice is singular")
        else:
            fractional = coords_array
            positions = fractional @ cell
    else:
        problems.append(f"bad_value: POSCAR scale line has {len(scale)} values; cell/positions not derived")
    return PoscarEvidence(
        sha256=sha256_bytes(raw), comment=comment, cell=cell, species_tokens=species_tokens, species=species,
        counts=counts, coordinate_mode=coordinate_mode, fractional=fractional, positions=positions,
        selective=selective, problems=problems,
    )


def parse_poscar_selective(path: Path) -> Any:
    """Just the (N, 3) selective-dynamics flags of a POSCAR/CONTCAR (None if absent)."""
    return parse_poscar(path).selective
