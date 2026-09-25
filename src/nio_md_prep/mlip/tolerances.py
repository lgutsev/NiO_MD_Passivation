"""Cross-engine tolerance records: measured candidates versus reviewed references.

A cross-engine comparison (``jobs.compare_engines``) reports how closely two
routes agree for one model and structure. Those numbers become a regression
criterion only after a person has reviewed them. This module keeps the two
kinds of record apart:

* :func:`candidate_record` turns measured metrics into a *candidate* file. It
  has the same layout as a reference but ``status`` says it is not accepted,
  and :func:`load_reference` refuses to use it.
* A *reference* is a candidate that a reviewer has edited: ``status`` set to
  :data:`ACCEPTED_STATUS` and ``reviewed_by`` filled in. Tests and scripts load
  it with :func:`load_reference`, which validates the layout.

No tolerance is ever invented here; every number comes from a measurement.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

CANDIDATE_STATUS = "candidate: measured, NOT accepted; requires human review"
ACCEPTED_STATUS = "accepted"

#: Metrics carried from ``ComparisonResult.as_dict()``. The energy difference is
#: stored as an absolute value; ``force_max_abs_error_eV_per_A`` is the largest
#: per-atom force-vector error, ``force_max_component_error_eV_per_A`` the largest
#: single Cartesian component error. Stress is ``None`` unless both routes
#: reported a virial.
METRICS = (
    "delta_energy_per_atom_eV",
    "force_rmse_eV_per_A",
    "force_max_abs_error_eV_per_A",
    "force_max_component_error_eV_per_A",
    "stress_max_abs_error_eV_per_A3",
)


def measured_metrics(report: Mapping[str, Any]) -> dict[str, dict[str, float | None]]:
    """``{"ase->lammps": {metric: value}}`` from a ``compare_engines`` report."""
    measured: dict[str, dict[str, float | None]] = {}
    for pair, values in report["comparisons"].items():
        row: dict[str, float | None] = {}
        for metric in METRICS:
            value = values.get(metric)
            if value is not None and metric == "delta_energy_per_atom_eV":
                value = abs(value)
            row[metric] = value
        measured[pair] = row
    return measured


def candidate_record(report: Mapping[str, Any], *, measured_on: Mapping[str, Any]) -> dict:
    """A reviewable candidate tolerance record built from one comparison report."""
    return {
        "status": CANDIDATE_STATUS,
        "reviewed_by": None,
        "note": (
            "Measured agreement between the routes named below, for one model and one "
            "structure. These are measurements, not acceptance criteria. To adopt them, "
            "review the numbers (and whether this structure and model are representative), "
            f"set status to '{ACCEPTED_STATUS}', fill in reviewed_by, and commit the file "
            "as tests/data/mlip_cross_engine_tolerances.json. Interoperability on a toy "
            "model does not validate a potential for NiO chemistry."
        ),
        "measured_on": dict(measured_on),
        "structure_sha256": report.get("structure", {}).get("sha256"),
        "energy_convention": report.get("energy_convention"),
        "tolerances": measured_metrics(report),
    }


def validate_reference(record: Mapping[str, Any]) -> None:
    """Raise ``ValueError`` unless ``record`` is a well-formed, reviewed reference."""
    if record.get("status") != ACCEPTED_STATUS:
        raise ValueError(
            f"tolerance file status is {record.get('status')!r}; only a reviewed file "
            f"with status {ACCEPTED_STATUS!r} can serve as a reference"
        )
    if not record.get("reviewed_by"):
        raise ValueError("a reference tolerance file must name its reviewer (reviewed_by)")
    if not record.get("measured_on"):
        raise ValueError("a reference tolerance file must record where it was measured (measured_on)")
    tolerances = record.get("tolerances")
    if not isinstance(tolerances, Mapping) or not tolerances:
        raise ValueError("a reference tolerance file needs a non-empty 'tolerances' table")
    for pair, limits in tolerances.items():
        if "->" not in pair:
            raise ValueError(f"tolerance key {pair!r} is not of the form 'reference->candidate'")
        for metric, value in limits.items():
            if metric not in METRICS:
                raise ValueError(f"unknown tolerance metric {pair}.{metric}")
            if value is not None and not value >= 0:
                raise ValueError(f"tolerance {pair}.{metric} must be >= 0 or null; got {value!r}")


def load_reference(path: Path) -> dict | None:
    """The reviewed reference at ``path``, ``None`` if absent; malformed or unreviewed files raise."""
    path = Path(path)
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    validate_reference(record)
    return record


def exceedances(measured: Mapping[str, Mapping[str, float | None]], reference: Mapping[str, Any]) -> list[str]:
    """Human-readable list of measured metrics above the reviewed limits (empty = within)."""
    problems = []
    for pair, values in measured.items():
        limits = reference["tolerances"].get(pair)
        if limits is None:
            problems.append(f"{pair}: no reviewed tolerance for this route pair")
            continue
        for metric, value in values.items():
            limit = limits.get(metric)
            if limit is None or value is None:
                continue
            if value > limit:
                problems.append(f"{pair} {metric} = {value:.3e} exceeds the reviewed limit {limit:.3e}")
    return problems


__all__ = [
    "ACCEPTED_STATUS",
    "CANDIDATE_STATUS",
    "METRICS",
    "candidate_record",
    "exceedances",
    "load_reference",
    "measured_metrics",
    "validate_reference",
]
