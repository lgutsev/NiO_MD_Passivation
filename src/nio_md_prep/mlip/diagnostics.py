"""Trajectory diagnostics with explicit physical meaning.

A short MD run is only worth its manifest entry if the number it reports
means what its name says. Two different quantities are easy to conflate:

**NVE energy drift.** Without a thermostat the total energy is a conserved
quantity of the exact dynamics, so its systematic change per unit time is a
test of the integrator, the timestep, the force/energy consistency of the
potential and every unit conversion in between. :func:`nve_drift` reports the
slope of a least-squares linear fit of the total energy per atom against time
(eV/atom/ps) and the largest excursion from the starting value -- the fit is
far less sensitive to the bounded oscillation of a symplectic integrator than
an endpoint difference is.

**Total-energy change under a thermostat.** A thermostat exchanges energy
with a bath by design, so ``E_total(t_end) - E_total(t_0)`` is not a
conservation test at all: a healthy Langevin run starting from a cold lattice
gains energy for as long as it heats up. :func:`descriptive_energy_change`
reports it under an explicit ``"descriptive; not conserved under a
thermostat"`` label. When an engine exposes the conserved quantity of its
extended system (LAMMPS ``econserve``, ASE ``NoseHooverChainNVT``
``get_conserved_energy()``, OpenMM Nose-Hoover heat-bath energy), that series
*is* a conservation test and :func:`ensemble_diagnostics` fits it the same
way as an NVE total energy.

:func:`ensemble_diagnostics` chooses between them by ensemble, which is the
only entry point engines should need. Everything here is pure Python on
plain sequences; nothing imports numpy.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .errors import ResultError

#: The label every thermostatted/barostatted total-energy change carries.
DESCRIPTIVE_LABEL = "descriptive; not conserved under a thermostat"
#: The label an NVE (or conserved-quantity) drift carries.
CONSERVATION_LABEL = "conservation test: linear-fit drift of a conserved energy"

FS_PER_PS = 1000.0


def nve_drift(
    times_fs: Sequence[float], total_energy_eV: Sequence[float], n_atoms: int
) -> dict:
    """Energy-conservation drift of a constant-energy trajectory.

    Returns ``energy_drift_eV_per_atom_per_ps`` (the least-squares slope of
    ``E(t)/N`` against ``t`` in ps), ``max_abs_energy_excursion_eV_per_atom``
    (``max |E(t) - E(t_0)| / N``), the number of samples and the sampled
    duration. Needs at least two samples at strictly increasing times.
    """
    times, energies = _series(times_fs, total_energy_eV, n_atoms)
    per_atom = [e / n_atoms for e in energies]
    times_ps = [t / FS_PER_PS for t in times]
    return {
        "energy_drift_eV_per_atom_per_ps": _slope(times_ps, per_atom),
        "max_abs_energy_excursion_eV_per_atom": max(abs(e - per_atom[0]) for e in per_atom),
        "samples": len(times),
        "duration_fs": times[-1] - times[0],
        "method": "least-squares linear fit of total energy per atom against time",
        "label": CONSERVATION_LABEL,
    }


def descriptive_energy_change(
    times_fs: Sequence[float], total_energy_eV: Sequence[float], n_atoms: int
) -> dict:
    """End-minus-start total energy per atom, labelled as *not* a conservation test."""
    times, energies = _series(times_fs, total_energy_eV, n_atoms)
    return {
        "total_energy_change_eV_per_atom": (energies[-1] - energies[0]) / n_atoms,
        "samples": len(times),
        "duration_fs": times[-1] - times[0],
        "label": DESCRIPTIVE_LABEL,
    }


def ensemble_diagnostics(
    ensemble: str | None,
    times_fs: Sequence[float],
    total_energy_eV: Sequence[float],
    n_atoms: int,
    *,
    conserved_energy_eV: Sequence[float] | None = None,
    conserved_quantity: str | None = None,
) -> dict:
    """The diagnostics appropriate to ``ensemble``, with explicit semantics.

    Always contains ``ensemble`` and ``conservation_test`` (bool). For
    ``nve`` it is :func:`nve_drift` of the total energy. For any other
    ensemble it is :func:`descriptive_energy_change` of the total energy,
    plus -- when the engine supplied ``conserved_energy_eV`` (the conserved
    quantity of the thermostatted/barostatted extended system, named by
    ``conserved_quantity``) -- a ``conserved_quantity_drift`` block computed
    exactly like an NVE drift. Only then is ``conservation_test`` true.
    """
    if ensemble == "nve":
        return {"ensemble": ensemble, "conservation_test": True,
                **nve_drift(times_fs, total_energy_eV, n_atoms)}
    report = {
        "ensemble": ensemble,
        "conservation_test": False,
        **descriptive_energy_change(times_fs, total_energy_eV, n_atoms),
    }
    if conserved_energy_eV is not None:
        if not conserved_quantity:
            raise ResultError(
                "a conserved-energy series needs the name of the quantity it is "
                "(e.g. 'econserve' or 'NoseHooverChainNVT.get_conserved_energy')"
            )
        drift = nve_drift(times_fs, conserved_energy_eV, n_atoms)
        report["conservation_test"] = True
        report["conserved_quantity"] = conserved_quantity
        report["conserved_quantity_drift"] = {
            "energy_drift_eV_per_atom_per_ps": drift["energy_drift_eV_per_atom_per_ps"],
            "max_abs_energy_excursion_eV_per_atom": drift[
                "max_abs_energy_excursion_eV_per_atom"
            ],
            "method": drift["method"],
            "label": CONSERVATION_LABEL,
        }
    return report


def series_diagnostics(ensemble: str | None, series: Mapping, n_atoms: int) -> dict:
    """:func:`ensemble_diagnostics` on an engine's ``energy_series`` payload.

    ``series`` is ``{"time_fs": [...], "total_energy_eV": [...]}`` plus,
    optionally, ``"conserved_energy_eV"`` and ``"conserved_quantity"``. This
    is the shape :meth:`~nio_md_prep.mlip.bridges.base.Bridge._trajectory`
    expects engines to return under ``payload["energy_series"]``.
    """
    try:
        times = series["time_fs"]
        totals = series["total_energy_eV"]
    except (KeyError, TypeError):
        raise ResultError(
            "energy_series must provide 'time_fs' and 'total_energy_eV' sequences"
        ) from None
    return ensemble_diagnostics(
        ensemble,
        times,
        totals,
        n_atoms,
        conserved_energy_eV=series.get("conserved_energy_eV"),
        conserved_quantity=series.get("conserved_quantity"),
    )


def temperature_ndof(n_atoms: int, *, n_fixed: int = 0, com_removed: bool = False) -> int:
    """Kinetic degrees of freedom used to turn kinetic energy into a temperature.

    ``3 N`` minus ``3`` per frozen atom, minus ``3`` more when the centre-of-
    mass momentum is removed and *no* atom is frozen (with frozen atoms the
    centre of mass is not a free coordinate, so no further DOF is removed).
    Engines report the value they used as ``temperature_ndof``; this function
    is the shared convention so two engines' temperatures can be compared.
    """
    if n_atoms < 1:
        raise ResultError(f"a temperature needs at least one atom; got {n_atoms}")
    if not 0 <= n_fixed <= n_atoms:
        raise ResultError(f"n_fixed must lie in 0..{n_atoms}; got {n_fixed}")
    ndof = 3 * (n_atoms - n_fixed)
    if com_removed and n_fixed == 0:
        ndof -= 3
    if ndof < 1:
        raise ResultError(
            f"no kinetic degrees of freedom remain ({n_atoms} atoms, {n_fixed} frozen, "
            f"centre of mass {'removed' if com_removed else 'kept'})"
        )
    return ndof


def _series(times_fs, energies_eV, n_atoms: int) -> tuple[list[float], list[float]]:
    if isinstance(n_atoms, bool) or not isinstance(n_atoms, int) or n_atoms < 1:
        raise ResultError(f"n_atoms must be a positive integer; got {n_atoms!r}")
    times = [float(t) for t in times_fs]
    energies = [float(e) for e in energies_eV]
    if len(times) != len(energies):
        raise ResultError(
            f"{len(times)} sample times for {len(energies)} energies; the series must align"
        )
    if len(times) < 2:
        raise ResultError(
            f"an energy diagnostic needs at least two samples; got {len(times)}"
        )
    if not all(math.isfinite(v) for v in times + energies):
        raise ResultError("the energy series contains a non-finite value")
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ResultError(
            "sample times must be strictly increasing (a repeated step, e.g. the first "
            "frame logged twice, must be removed before fitting)"
        )
    return times, energies


def _slope(x: Sequence[float], y: Sequence[float]) -> float:
    n = len(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sxx = sum((xi - mean_x) ** 2 for xi in x)
    sxy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    return sxy / sxx


__all__ = [
    "DESCRIPTIVE_LABEL",
    "CONSERVATION_LABEL",
    "nve_drift",
    "descriptive_energy_change",
    "ensemble_diagnostics",
    "series_diagnostics",
    "temperature_ndof",
]
