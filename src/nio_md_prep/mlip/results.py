"""The common result object every bridge returns, and cross-engine comparison.

A :class:`PotentialResult` is always canonical: eV, eV/Angstrom,
eV/Angstrom^3, with the energy convention it was produced under recorded
explicitly. Bridges convert from their engine's native units on the way out
and record what those native units were, so a wrong number can be traced to a
conversion rather than guessed at.

:func:`compare_results` refuses to compare two results whose energy
conventions differ and which cannot be reconciled. That refusal is the whole
point: an ASE/OpenMM MACE comparison that ignores atomic self-energies looks
catastrophically wrong while the physics is fine.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from typing import Any

from .errors import EnergyConventionError, ResultError
from .units import CANONICAL, convert_energy_convention


@dataclass(frozen=True)
class PotentialResult:
    """One evaluated geometry, in canonical units, under a declared convention.

    Validated on construction, so an engine that returns a NaN energy, a
    force row that is not a 3-vector, ghost-atom rows, or a per-atom array of
    the wrong length fails with :class:`~nio_md_prep.mlip.errors.ResultError`
    instead of producing a plausible-looking result: at least one atom; every
    energy, force and stress component finite; exactly one 3-component force
    row per atom; a 6-component Voigt stress; one per-atom energy per atom.
    """

    energy_eV: float
    forces_eV_per_A: tuple[tuple[float, float, float], ...]
    symbols: tuple[str, ...]
    energy_convention: str
    engine: str
    potential: str
    implementation: str
    native_units: str = CANONICAL.name
    stress_eV_per_A3: tuple[float, ...] | None = None
    per_atom_energy_eV: tuple[float, ...] | None = None
    wall_time_s: float | None = None
    extras: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "symbols", tuple(str(s) for s in self.symbols))
        n_atoms = len(self.symbols)
        if n_atoms < 1:
            raise ResultError("a result must describe at least one atom; got none")
        energy = float(self.energy_eV)
        if not math.isfinite(energy):
            raise ResultError(f"energy is not finite: {self.energy_eV!r}")
        object.__setattr__(self, "energy_eV", energy)
        forces = tuple(tuple(float(c) for c in row) for row in self.forces_eV_per_A)
        if len(forces) != n_atoms:
            raise ResultError(
                f"{len(forces)} force rows for {n_atoms} atoms (ghost or missing atoms "
                "in the engine output?)"
            )
        for index, row in enumerate(forces):
            if len(row) != 3:
                raise ResultError(
                    f"force row {index} has {len(row)} components; expected 3"
                )
            if not all(math.isfinite(c) for c in row):
                raise ResultError(f"force on atom {index} is not finite: {row}")
        object.__setattr__(self, "forces_eV_per_A", forces)
        if self.stress_eV_per_A3 is not None:
            stress = tuple(float(v) for v in self.stress_eV_per_A3)
            if len(stress) != 6:
                raise ResultError(
                    "stress must be the 6-component Voigt vector "
                    f"(xx, yy, zz, yz, xz, xy); got {len(stress)} components"
                )
            if not all(math.isfinite(v) for v in stress):
                raise ResultError(f"stress is not finite: {stress}")
            object.__setattr__(self, "stress_eV_per_A3", stress)
        if self.per_atom_energy_eV is not None:
            per_atom = tuple(float(v) for v in self.per_atom_energy_eV)
            if len(per_atom) != n_atoms:
                raise ResultError(
                    f"{len(per_atom)} per-atom energies for {n_atoms} atoms"
                )
            if not all(math.isfinite(v) for v in per_atom):
                raise ResultError("a per-atom energy is not finite")
            object.__setattr__(self, "per_atom_energy_eV", per_atom)

    @property
    def n_atoms(self) -> int:
        return len(self.symbols)

    @property
    def energy_per_atom_eV(self) -> float:
        return self.energy_eV / self.n_atoms

    @property
    def max_force_eV_per_A(self) -> float:
        return max(
            (math.sqrt(sum(c * c for c in row)) for row in self.forces_eV_per_A), default=0.0
        )

    def in_convention(self, target: str, *, atomic_reference_energies=None) -> "PotentialResult":
        """Return this result restated under ``target``.

        Forces and stresses are returned unchanged: the two conventions differ
        by a composition-dependent constant, whose gradient is zero.
        """
        if target == self.energy_convention:
            return self
        energy = convert_energy_convention(
            self.energy_eV,
            self.symbols,
            source=self.energy_convention,
            target=target,
            atomic_reference_energies=atomic_reference_energies,
        )
        per_atom = self.per_atom_energy_eV
        if per_atom is not None and atomic_reference_energies:
            sign = 1.0 if target == "total" else -1.0
            per_atom = tuple(
                value + sign * atomic_reference_energies[symbol]
                for value, symbol in zip(per_atom, self.symbols)
            )
        return replace(
            self,
            energy_eV=energy,
            energy_convention=target,
            per_atom_energy_eV=per_atom,
        )

    def as_dict(self, *, include_arrays: bool = True) -> dict:
        payload: dict[str, Any] = {
            "engine": self.engine,
            "potential": self.potential,
            "implementation": self.implementation,
            "native_units": self.native_units,
            "canonical_units": {
                "energy": "eV",
                "forces": "eV/Angstrom",
                "stress": "eV/Angstrom^3",
                "length": "Angstrom",
            },
            "energy_convention": self.energy_convention,
            "n_atoms": self.n_atoms,
            "energy_eV": self.energy_eV,
            "energy_per_atom_eV": self.energy_per_atom_eV,
            "max_force_eV_per_A": self.max_force_eV_per_A,
            "wall_time_s": self.wall_time_s,
            "extras": dict(self.extras),
        }
        if include_arrays:
            payload["forces_eV_per_A"] = [list(row) for row in self.forces_eV_per_A]
            payload["stress_eV_per_A3"] = (
                list(self.stress_eV_per_A3) if self.stress_eV_per_A3 else None
            )
            payload["per_atom_energy_eV"] = (
                list(self.per_atom_energy_eV) if self.per_atom_energy_eV else None
            )
        return payload


@dataclass(frozen=True)
class TrajectoryResult:
    """A short diagnostic trajectory. Deliberately not a production MD result.

    What the engine *reported* is kept apart from what was *requested*:

    ``steps``
        The requested number of steps (``simulation.steps``).
    ``steps_completed``
        The step counter the engine reported after the run (ASE
        ``dyn.nsteps``, LAMMPS ``ntimestep``, OpenMM ``currentStep``);
        ``None`` when the engine did not report one. A run that stopped
        early must say so here, not look complete.
    ``frames_written``
        Frames counted by reading the written trajectory back, never
        estimated from ``steps // interval``; ``None`` when not verified.
    ``final_positions_angstrom`` / ``final_cell_angstrom`` / ``final_pbc``
        The final geometry in the *source* basis: the same Cartesian frame and
        cell-vector convention as the input structure, after rotating back
        from any engine-internal frame (LAMMPS restricted triclinic, OpenMM
        reduced box). The cell changes under NPT, so it must be returned.
    ``integrator_resolved``
        What actually integrated: class or fix names and every resolved
        parameter (damping times after defaults, friction, chain length,
        coupling, compressibility, seed ...).
    ``diagnostics``
        :func:`nio_md_prep.mlip.diagnostics.ensemble_diagnostics` output.
        ``conservation_test`` says whether any number in it is a
        conservation test; an NVT total-energy change is labelled
        descriptive. Empty with ``available: False`` when the engine
        supplied no energy series.
    ``temperature_ndof``
        The kinetic degrees of freedom the temperatures were computed with
        (see :func:`nio_md_prep.mlip.diagnostics.temperature_ndof`).
    ``constraints``
        What the engine applied, e.g. ``{"fixed_atoms": [...], "n_fixed": 4,
        "method": "zero particle mass"}``; ``None`` when there were none or
        the engine did not report them.
    """

    steps: int
    timestep_fs: float
    ensemble: str | None
    trajectory_path: str | None
    log_path: str | None
    initial: PotentialResult
    final: PotentialResult
    steps_completed: int | None = None
    frames_written: int | None = None
    final_positions_angstrom: tuple[tuple[float, float, float], ...] | None = None
    final_cell_angstrom: tuple[tuple[float, float, float], ...] | None = None
    final_pbc: tuple[bool, bool, bool] | None = None
    integrator_resolved: dict[str, Any] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    temperature_start_K: float | None = None
    temperature_end_K: float | None = None
    max_temperature_K: float | None = None
    temperature_ndof: int | None = None
    constraints: dict[str, Any] | None = None
    wall_time_s: float | None = None
    extras: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n_atoms = self.initial.n_atoms
        if self.final.symbols != self.initial.symbols:
            raise ResultError(
                "the final frame does not have the initial frame's atoms in the same order"
            )
        for name in ("steps_completed", "frames_written", "temperature_ndof"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or int(value) != value or value < 0):
                raise ResultError(f"{name} must be a non-negative integer; got {value!r}")
        if self.final_positions_angstrom is not None:
            positions = tuple(tuple(float(c) for c in row) for row in self.final_positions_angstrom)
            if len(positions) != n_atoms or any(len(row) != 3 for row in positions):
                raise ResultError(
                    f"final positions must be {n_atoms} rows of 3 coordinates"
                )
            if not all(math.isfinite(c) for row in positions for c in row):
                raise ResultError("a final position is not finite")
            object.__setattr__(self, "final_positions_angstrom", positions)
        if self.final_cell_angstrom is not None:
            cell = tuple(tuple(float(c) for c in row) for row in self.final_cell_angstrom)
            if len(cell) != 3 or any(len(row) != 3 for row in cell):
                raise ResultError("the final cell must be 3 rows of 3 components")
            if not all(math.isfinite(c) for row in cell for c in row):
                raise ResultError("a final cell component is not finite")
            object.__setattr__(self, "final_cell_angstrom", cell)
        if self.final_pbc is not None:
            pbc = tuple(bool(v) for v in self.final_pbc)
            if len(pbc) != 3:
                raise ResultError(f"final_pbc needs three entries; got {self.final_pbc!r}")
            object.__setattr__(self, "final_pbc", pbc)

    @property
    def frames(self) -> int | None:
        """Backwards-compatible alias of :attr:`frames_written`."""
        return self.frames_written

    def as_dict(self, *, include_arrays: bool = False) -> dict:
        payload = {
            "steps": self.steps,
            "steps_completed": self.steps_completed,
            "timestep_fs": self.timestep_fs,
            "ensemble": self.ensemble,
            "frames_written": self.frames_written,
            "trajectory_path": self.trajectory_path,
            "log_path": self.log_path,
            "temperature_start_K": self.temperature_start_K,
            "temperature_end_K": self.temperature_end_K,
            "max_temperature_K": self.max_temperature_K,
            "temperature_ndof": self.temperature_ndof,
            "integrator_resolved": dict(self.integrator_resolved),
            "diagnostics": dict(self.diagnostics),
            "constraints": dict(self.constraints) if self.constraints else None,
            "final_cell_angstrom": (
                [list(row) for row in self.final_cell_angstrom]
                if self.final_cell_angstrom is not None
                else None
            ),
            "final_pbc": list(self.final_pbc) if self.final_pbc is not None else None,
            "wall_time_s": self.wall_time_s,
            "initial": self.initial.as_dict(include_arrays=include_arrays),
            "final": self.final.as_dict(include_arrays=include_arrays),
            "extras": dict(self.extras),
        }
        if include_arrays:
            payload["final_positions_angstrom"] = (
                [list(row) for row in self.final_positions_angstrom]
                if self.final_positions_angstrom is not None
                else None
            )
        return payload


@dataclass(frozen=True)
class ComparisonResult:
    """The metrics of the cross-engine single-point equivalence check."""

    reference: str
    candidate: str
    n_atoms: int
    energy_convention: str
    delta_energy_eV: float
    delta_energy_per_atom_eV: float
    force_rmse_eV_per_A: float
    force_max_abs_error_eV_per_A: float
    force_max_component_error_eV_per_A: float
    stress_max_abs_error_eV_per_A3: float | None = None
    reference_max_force_eV_per_A: float | None = None
    notes: tuple[str, ...] = ()

    def as_dict(self) -> dict:
        return {
            "reference": self.reference,
            "candidate": self.candidate,
            "n_atoms": self.n_atoms,
            "energy_convention": self.energy_convention,
            "delta_energy_eV": self.delta_energy_eV,
            "delta_energy_per_atom_eV": self.delta_energy_per_atom_eV,
            "force_rmse_eV_per_A": self.force_rmse_eV_per_A,
            "force_max_abs_error_eV_per_A": self.force_max_abs_error_eV_per_A,
            "force_max_component_error_eV_per_A": self.force_max_component_error_eV_per_A,
            "stress_max_abs_error_eV_per_A3": self.stress_max_abs_error_eV_per_A3,
            "reference_max_force_eV_per_A": self.reference_max_force_eV_per_A,
            "notes": list(self.notes),
        }


def compare_results(
    reference: PotentialResult,
    candidate: PotentialResult,
    *,
    atomic_reference_energies=None,
    convention: str | None = None,
) -> ComparisonResult:
    """Compare two single-point results after harmonising their conventions.

    Raises :class:`EnergyConventionError` rather than comparing energies
    produced under different conventions. Force and stress comparison is
    unaffected by convention, so the caller still gets a useful message that
    names the shortfall instead of a meaningless energy difference.
    """
    if reference.symbols != candidate.symbols:
        raise ResultError(
            "cannot compare results for different structures: "
            f"{len(reference.symbols)} vs {len(candidate.symbols)} atoms, or a different "
            "element order"
        )
    target = convention or reference.energy_convention
    notes: list[str] = []
    if reference.energy_convention != candidate.energy_convention:
        notes.append(
            f"{reference.engine} reported the {reference.energy_convention} energy and "
            f"{candidate.engine} the {candidate.energy_convention} energy; both were "
            f"converted to {target} using the model's atomic reference energies (E0)"
        )
    try:
        left = reference.in_convention(
            target, atomic_reference_energies=atomic_reference_energies
        )
        right = candidate.in_convention(
            target, atomic_reference_energies=atomic_reference_energies
        )
    except EnergyConventionError as exc:
        raise EnergyConventionError(
            f"refusing to compare {reference.engine} ({reference.energy_convention}) with "
            f"{candidate.engine} ({candidate.energy_convention}): {exc}"
        ) from exc

    delta_energy = right.energy_eV - left.energy_eV
    squared = 0.0
    max_component = 0.0
    max_vector = 0.0
    for a, b in zip(left.forces_eV_per_A, right.forces_eV_per_A):
        diff = [bi - ai for ai, bi in zip(a, b)]
        squared += sum(d * d for d in diff)
        max_component = max(max_component, max(abs(d) for d in diff))
        max_vector = max(max_vector, math.sqrt(sum(d * d for d in diff)))
    n_atoms = left.n_atoms
    force_rmse = math.sqrt(squared / (3 * n_atoms)) if n_atoms else 0.0

    stress_error = None
    if left.stress_eV_per_A3 is not None and right.stress_eV_per_A3 is not None:
        stress_error = max(
            abs(b - a) for a, b in zip(left.stress_eV_per_A3, right.stress_eV_per_A3)
        )
    elif (left.stress_eV_per_A3 is None) != (right.stress_eV_per_A3 is None):
        missing = left.engine if left.stress_eV_per_A3 is None else right.engine
        notes.append(f"stress not compared: {missing} did not report a virial")

    return ComparisonResult(
        reference=f"{left.engine}:{left.implementation}",
        candidate=f"{right.engine}:{right.implementation}",
        n_atoms=n_atoms,
        energy_convention=target,
        delta_energy_eV=delta_energy,
        delta_energy_per_atom_eV=delta_energy / n_atoms if n_atoms else 0.0,
        force_rmse_eV_per_A=force_rmse,
        force_max_abs_error_eV_per_A=max_vector,
        force_max_component_error_eV_per_A=max_component,
        stress_max_abs_error_eV_per_A3=stress_error,
        reference_max_force_eV_per_A=left.max_force_eV_per_A,
        notes=tuple(notes),
    )


def rms(values: Sequence[float]) -> float:
    values = list(values)
    if not values:
        return 0.0
    return math.sqrt(sum(v * v for v in values) / len(values))


__all__ = [
    "PotentialResult",
    "TrajectoryResult",
    "ComparisonResult",
    "compare_results",
    "rms",
]
