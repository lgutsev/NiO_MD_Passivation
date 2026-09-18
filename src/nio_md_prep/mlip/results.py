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

from .errors import EnergyConventionError
from .units import CANONICAL, convert_energy_convention


@dataclass(frozen=True)
class PotentialResult:
    """One evaluated geometry, in canonical units, under a declared convention."""

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
        object.__setattr__(self, "symbols", tuple(self.symbols))
        object.__setattr__(
            self,
            "forces_eV_per_A",
            tuple(tuple(float(c) for c in row) for row in self.forces_eV_per_A),
        )
        if len(self.forces_eV_per_A) != len(self.symbols):
            raise ValueError(
                f"{len(self.forces_eV_per_A)} force rows for {len(self.symbols)} atoms"
            )
        if self.stress_eV_per_A3 is not None:
            stress = tuple(float(v) for v in self.stress_eV_per_A3)
            if len(stress) != 6:
                raise ValueError(
                    "stress must be the 6-component Voigt vector "
                    f"(xx, yy, zz, yz, xz, xy); got {len(stress)} components"
                )
            object.__setattr__(self, "stress_eV_per_A3", stress)
        if self.per_atom_energy_eV is not None:
            object.__setattr__(
                self, "per_atom_energy_eV", tuple(float(v) for v in self.per_atom_energy_eV)
            )

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
    """A short diagnostic trajectory. Deliberately not a production MD result."""

    steps: int
    timestep_fs: float
    ensemble: str | None
    frames: int
    trajectory_path: str | None
    log_path: str | None
    initial: PotentialResult
    final: PotentialResult
    temperature_start_K: float | None = None
    temperature_end_K: float | None = None
    total_energy_drift_eV_per_atom: float | None = None
    max_temperature_K: float | None = None
    wall_time_s: float | None = None
    extras: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "steps": self.steps,
            "timestep_fs": self.timestep_fs,
            "ensemble": self.ensemble,
            "frames": self.frames,
            "trajectory_path": self.trajectory_path,
            "log_path": self.log_path,
            "temperature_start_K": self.temperature_start_K,
            "temperature_end_K": self.temperature_end_K,
            "max_temperature_K": self.max_temperature_K,
            "total_energy_drift_eV_per_atom": self.total_energy_drift_eV_per_atom,
            "wall_time_s": self.wall_time_s,
            "initial": self.initial.as_dict(include_arrays=False),
            "final": self.final.as_dict(include_arrays=False),
            "extras": dict(self.extras),
        }


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
        raise ValueError(
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
