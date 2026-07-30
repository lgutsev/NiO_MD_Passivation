"""Loading and aggregation of externally supplied experimental device results.

This module never bundles real measurement data itself. Callers pass a CSV
path at run time (kept outside version control while the data is
unpublished); only the column schema and aggregation logic live here.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
import csv
import math

REQUIRED_COLUMNS = (
    "secondary",
    "assembly",
    "voc_v",
    "jsc_ma_cm2",
    "ff_percent",
    "pce_percent",
)

_METRIC_FIELDS = ("voc_v", "jsc_ma_cm2", "ff_percent", "pce_percent")


def load_device_results(path: Path) -> list[dict[str, Any]]:
    """Parse a CSV of experimental JV measurements.

    Required columns: secondary, assembly, voc_v, jsc_ma_cm2, ff_percent,
    pce_percent. ``batch``, ``device_id``, and ``scan_direction`` are
    optional free-text provenance columns carried through unchanged when
    present. Supplying ``device_id`` is strongly recommended whenever more
    than one device in a batch has forward/reverse scans.
    """
    path = Path(path)
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: CSV has no header row")
        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"{path}: missing required column(s): {', '.join(missing)}"
            )
        rows: list[dict[str, Any]] = []
        for line_number, record in enumerate(reader, start=2):
            try:
                rows.append(
                    {
                        "secondary": record["secondary"].strip(),
                        "assembly": record["assembly"].strip(),
                        "batch": (record.get("batch") or "").strip(),
                        "device_id": (
                            record.get("device_id") or ""
                        ).strip(),
                        "scan_direction": (
                            record.get("scan_direction") or ""
                        ).strip(),
                        "voc_v": float(record["voc_v"]),
                        "jsc_ma_cm2": float(record["jsc_ma_cm2"]),
                        "ff_percent": float(record["ff_percent"]),
                        "pce_percent": float(record["pce_percent"]),
                    }
                )
            except (KeyError, ValueError) as exc:
                raise ValueError(
                    f"{path}: line {line_number}: {exc}"
                ) from exc
    if not rows:
        raise ValueError(f"{path}: no data rows found")
    return rows


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _sample_std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = _mean(values)
    return math.sqrt(
        sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    )


def aggregate_device_results(
    rows: list[dict[str, Any]],
) -> dict[tuple[str, str], dict[str, Any]]:
    """Aggregate measurements without treating F/R scans as new devices.

    The calculation is deliberately hierarchical:

    1. average repeated scans for one explicit ``device_id``;
    2. average independent devices within each batch;
    3. give every batch equal weight in the final condition mean.

    When ``device_id`` is absent, rows with a non-empty ``scan_direction`` in
    the same condition and batch are conservatively treated as repeated scans
    of one device. Rows without a scan direction are treated as independent
    measurement units. If a batch contains F/R scans for more than one
    device, an explicit ``device_id`` is required to avoid ambiguity.

    ``std`` is the sample SD across batch means when two or more batches are
    present, otherwise across independent device/measurement-unit means.
    """
    if not rows:
        raise ValueError("device-result aggregation requires at least one row")

    measurement_groups: dict[
        tuple[str, str, str, str], list[dict[str, Any]]
    ] = {}
    for row_index, row in enumerate(rows):
        condition = (row["secondary"], row["assembly"])
        batch = str(row.get("batch") or "").strip()
        device_id = str(row.get("device_id") or "").strip()
        scan_direction = str(row.get("scan_direction") or "").strip()
        if device_id:
            unit_id = f"device:{device_id}"
        elif scan_direction:
            unit_id = "inferred-scan-pair"
        else:
            unit_id = f"row:{row_index}"
        measurement_groups.setdefault(
            (*condition, batch, unit_id), []
        ).append(row)

    independent_units: dict[
        tuple[str, str, str], list[dict[str, Any]]
    ] = {}
    for (secondary, assembly, batch, unit_id), measurements in (
        measurement_groups.items()
    ):
        unit = {
            field: _mean([float(row[field]) for row in measurements])
            for field in _METRIC_FIELDS
        }
        unit["measurement_count"] = len(measurements)
        unit["inferred_scan_pair"] = unit_id == "inferred-scan-pair"
        independent_units.setdefault(
            (secondary, assembly, batch), []
        ).append(unit)

    batch_results: dict[
        tuple[str, str], list[dict[str, Any]]
    ] = {}
    for (secondary, assembly, batch), units in independent_units.items():
        batch_result = {
            field: _mean([float(unit[field]) for unit in units])
            for field in _METRIC_FIELDS
        }
        batch_result.update(
            {
                "batch": batch,
                "independent_unit_count": len(units),
                "measurement_count": sum(
                    int(unit["measurement_count"]) for unit in units
                ),
                "inferred_scan_pair_count": sum(
                    int(bool(unit["inferred_scan_pair"])) for unit in units
                ),
                "unit_values": units,
            }
        )
        batch_results.setdefault((secondary, assembly), []).append(
            batch_result
        )

    aggregated: dict[tuple[str, str], dict[str, Any]] = {}
    for key, batches in batch_results.items():
        all_units = [
            unit
            for batch in batches
            for unit in batch["unit_values"]
        ]
        entry: dict[str, Any] = {
            # Backward-compatible alias: n now means independent units, not
            # raw CSV rows.
            "n": len(all_units),
            "n_independent_units": len(all_units),
            "n_measurements": sum(
                int(batch["measurement_count"]) for batch in batches
            ),
            "n_batches": len(batches),
            "batches": sorted(str(batch["batch"]) for batch in batches),
            "inferred_scan_pair_count": sum(
                int(batch["inferred_scan_pair_count"]) for batch in batches
            ),
        }
        for field in _METRIC_FIELDS:
            batch_values = [float(batch[field]) for batch in batches]
            mean = _mean(batch_values)
            if len(batch_values) > 1:
                uncertainty_values = batch_values
                uncertainty_basis = "between_batch"
            else:
                uncertainty_values = [
                    float(unit[field]) for unit in all_units
                ]
                uncertainty_basis = "between_independent_unit"
            entry[field] = {
                "mean": mean,
                "std": _sample_std(uncertainty_values),
                "uncertainty_basis": uncertainty_basis,
            }
        aggregated[key] = entry
    return aggregated
