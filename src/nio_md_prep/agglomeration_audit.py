from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ATTENTION_STATES = {
    "FINAL_UNMARKED",
    "HASH_MISMATCH",
    "MARKER_WITHOUT_FINAL",
    "MISSING_DIRECTORY",
}


def _campaign_manifests(root: Path) -> list[tuple[Path, dict]]:
    candidates = {
        path.parent for path in root.rglob("agglomeration_manifest.json")
    }
    candidates.update(path for path in root.rglob("*_agglo") if path.is_dir())
    if root.name.endswith("_agglo"):
        candidates.add(root)
    campaigns: list[tuple[Path, dict]] = []
    for campaign in sorted(candidates):
        manifest_path = campaign / "agglomeration_manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except OSError:
            campaigns.append((campaign, {"_manifest_error": "manifest is missing"}))
            continue
        except json.JSONDecodeError as exc:
            campaigns.append(
                (campaign, {"_manifest_error": f"manifest is invalid JSON: {exc}"})
            )
            continue
        if isinstance(manifest, dict) and isinstance(manifest.get("replicas"), list):
            campaigns.append((campaign, manifest))
        elif campaign.name.endswith("_agglo"):
            campaigns.append(
                (campaign, {"_manifest_error": "manifest lacks campaign replicas"})
            )
    return campaigns


def _campaign_label(root: Path, campaign: Path) -> str:
    try:
        relative = campaign.relative_to(root)
    except ValueError:
        return str(campaign)
    return "." if str(relative) == "." else str(relative)


def _audit_campaign(root: Path, campaign: Path, manifest: dict) -> dict[str, object]:
    row: dict[str, object] = {
        "campaign": _campaign_label(root, campaign),
        "path": str(campaign.resolve()),
        "campaign_status": "ATTENTION",
        "manifest_status": str(manifest.get("status", "UNKNOWN")),
        "total": 0,
        "completed": 0,
        "pending": 0,
        "attention": 0,
        "percent_complete": 0.0,
        "counts_by_status": {},
        "detail": "",
    }
    if manifest.get("_manifest_error"):
        row["detail"] = str(manifest["_manifest_error"])
        return row
    case_list = campaign / "xtb_cases.tsv"
    audit = campaign / "audit_xtb_runs.py"
    if not case_list.is_file():
        if manifest.get("status") == "PACKMOL_REQUIRED":
            row["campaign_status"] = "PACKMOL_REQUIRED"
            row["detail"] = "Packmol structures must be generated before xTB"
        elif not manifest.get("xtb", {}).get("enabled", False):
            row["campaign_status"] = "XTB_DISABLED"
            row["detail"] = "xTB is disabled for this campaign"
        else:
            row["detail"] = "xtb_cases.tsv is missing"
        return row
    if not audit.is_file():
        row["detail"] = (
            "audit_xtb_runs.py is missing; regenerate this campaign's xTB launcher"
        )
        return row

    result = subprocess.run(
        [sys.executable, str(audit), "--summary-only"],
        cwd=campaign,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        row["detail"] = (
            f"audit failed with exit code {result.returncode}: "
            f"{(result.stderr or result.stdout).strip()}"
        )
        return row
    progress_path = campaign / "xtb_progress.json"
    try:
        progress = json.loads(progress_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        row["detail"] = f"audit did not produce readable xtb_progress.json: {exc}"
        return row

    counts = progress.get("counts_by_status", {})
    if not isinstance(counts, dict):
        counts = {}
    total = int(progress.get("total", 0))
    completed = int(progress.get("completed", 0))
    pending = int(progress.get("pending", max(0, total - completed)))
    attention = sum(int(counts.get(status, 0)) for status in ATTENTION_STATES)
    row.update(
        total=total,
        completed=completed,
        pending=pending,
        attention=attention,
        percent_complete=round(100.0 if total == 0 else 100.0 * completed / total, 1),
        counts_by_status=dict(sorted(counts.items())),
    )
    if attention:
        row["campaign_status"] = "ATTENTION"
        row["detail"] = "one or more cases have inconsistent/missing completion data"
    elif pending:
        row["campaign_status"] = "IN_PROGRESS"
        row["detail"] = "xTB work remains"
    else:
        row["campaign_status"] = "COMPLETE"
        row["detail"] = "all xTB cases are safely complete"
    return row


def audit_agglomerations(root: Path, output_prefix: Path | None = None) -> dict:
    root = root.resolve()
    if not root.is_dir():
        raise ValueError(f"agglomeration audit root is not a directory: {root}")
    campaigns = _campaign_manifests(root)
    if not campaigns:
        raise ValueError(f"no agglomeration campaigns found under {root}")

    rows = [_audit_campaign(root, campaign, manifest) for campaign, manifest in campaigns]
    status_counts = Counter(str(row["campaign_status"]) for row in rows)
    totals = {
        "campaigns": len(rows),
        "cases": sum(int(row["total"]) for row in rows),
        "completed": sum(int(row["completed"]) for row in rows),
        "pending": sum(int(row["pending"]) for row in rows),
        "attention": sum(int(row["attention"]) for row in rows),
    }
    report = {
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        **totals,
        "campaign_counts_by_status": dict(sorted(status_counts.items())),
        "campaign_results": rows,
    }

    prefix = (output_prefix or (root / "agglomeration_xtb_audit")).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = prefix.with_suffix(".csv")
    json_path = prefix.with_suffix(".json")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "campaign", "path", "campaign_status", "manifest_status", "total",
            "completed", "pending", "attention", "percent_complete",
            "counts_by_status", "detail",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            serialized = dict(row)
            serialized["counts_by_status"] = json.dumps(
                serialized["counts_by_status"], sort_keys=True
            )
            writer.writerow(serialized)
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print("Campaign                         Status           Complete   Pending  Attention")
    print("-------------------------------- ---------------- ---------- --------- ---------")
    for row in rows:
        print(
            f"{str(row['campaign'])[:32]:32} "
            f"{str(row['campaign_status']):16} "
            f"{int(row['completed']):4}/{int(row['total']):<5} "
            f"{int(row['pending']):9} {int(row['attention']):9}"
        )
    for row in rows:
        if row["campaign_status"] == "ATTENTION":
            print(f"ATTENTION {row['campaign']}: {row['detail']}")
    print(
        f"Overall xTB progress: {totals['completed']}/{totals['cases']} safely complete; "
        f"{totals['pending']} pending or requiring attention across "
        f"{totals['campaigns']} campaign(s)."
    )
    print(f"Reports written: {csv_path} and {json_path}")
    return report
