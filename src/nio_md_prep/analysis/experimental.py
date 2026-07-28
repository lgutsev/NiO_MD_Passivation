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
    pce_percent. ``batch`` and ``scan_direction`` are optional free-text
    provenance columns carried through unchanged when present.
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


def aggregate_device_results(
    rows: list[dict[str, Any]],
) -> dict[tuple[str, str], dict[str, Any]]:
    """Group device rows by (secondary, assembly); report mean/std/n per metric.

    Pools every row for a given (secondary, assembly) pair regardless of
    ``batch``/``scan_direction`` -- a deliberate simplification when rows
    come from more than one measurement campaign; group by batch separately
    beforehand if that pooling isn't appropriate for a given dataset.
    """
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault((row["secondary"], row["assembly"]), []).append(row)

    aggregated: dict[tuple[str, str], dict[str, Any]] = {}
    for key, group_rows in groups.items():
        entry: dict[str, Any] = {"n": len(group_rows)}
        for field in _METRIC_FIELDS:
            values = [row[field] for row in group_rows]
            mean = sum(values) / len(values)
            if len(values) > 1:
                variance = sum((v - mean) ** 2 for v in values) / (
                    len(values) - 1
                )
                std = math.sqrt(variance)
            else:
                std = 0.0
            entry[field] = {"mean": mean, "std": std}
        aggregated[key] = entry
    return aggregated
