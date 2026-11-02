"""Collect accessibility_summary.json files into a compare-routes manifest.

Writes one row per analyzed trajectory with the analyzed step/frame window,
so the manifest doubles as the record of which trajectories and time windows
entered the route comparison.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import csv
import json

ROOTS = ("prepared", "prepared-rerun", "prepared-rerun2", "prepared-lego", "prepared-lego2")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for root in ROOTS:
        for summary_path in sorted((args.repository / root).glob("*/accessibility-analysis-*/accessibility_summary.json")):
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            stage = summary_path.parent.name.removeprefix("accessibility-analysis-")
            rows.append({
                "prepared_root": root,
                "system": summary_path.parent.parent.name,
                "stage": stage,
                "summary_json": str(summary_path),
                "trajectory": summary["trajectory"],
                "frames_analyzed": summary["frames_analyzed"],
                "first_step": summary["first_step"],
                "last_step": summary["last_step"],
                "frame_interval_ps": summary["frame_interval_ps"],
                "model_profile": summary["model_assessment"].get("profile"),
                "git_commit": summary["provenance"].get("git_commit"),
            })
    if not rows:
        raise SystemExit("no accessibility_summary.json files found")
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} analyzed trajectories -> {args.output}")


if __name__ == "__main__":
    main()
