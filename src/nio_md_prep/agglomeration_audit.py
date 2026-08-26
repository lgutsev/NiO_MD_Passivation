from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ATTENTION_STATES = {
    "FINAL_UNMARKED",
    "HASH_MISMATCH",
    "MARKER_WITHOUT_FINAL",
    "MISSING_DIRECTORY",
}

# Campaign-level statuses that mean "the audit itself could not run to
# completion" (missing manifest, missing auditor, crashed subprocess,
# unreadable/self-inconsistent output) are kept as "ERROR", distinct from
# "ATTENTION" which means the audit ran fine and found specific case(s) that
# need a human to look at them. Conflating the two was the original tool's
# main blind spot: a campaign could show up as "ATTENTION" with no way to
# tell whether that meant "one xTB case has a hash mismatch" or "the audit
# script for this campaign doesn't even exist".
ERROR_DETAIL_DEFAULT = "audit did not run to completion for this campaign"

DEFAULT_AUDIT_TIMEOUT_SECONDS = 300.0
MAX_DETAIL_CHARS = 4000

CAMPAIGN_FIELDS = [
    "campaign", "path", "campaign_status", "manifest_status", "total",
    "completed", "pending", "attention", "percent_complete",
    "counts_by_status", "data_consistent", "audit_seconds", "detail",
]
CASE_FIELDS = [
    "campaign", "campaign_path", "array_index", "case", "status", "charge",
    "uhf", "detail",
]


def _as_int(value: object, default: int = 0) -> int:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default


def _truncate(text: str, limit: int = MAX_DETAIL_CHARS) -> str:
    text = text.strip()
    if len(text) <= limit:
        return text
    return text[:limit] + f" ... [truncated, {len(text) - limit} more chars]"


def _discover_campaign_dirs(root: Path) -> tuple[set[Path], list[str]]:
    """Find candidate campaign directories under root.

    Uses os.walk with an error callback instead of Path.rglob so that one
    unreadable/locked subdirectory (permissions, a broken symlink, a network
    mount that dropped mid-scan) does not abort the whole audit -- it is
    recorded as a warning and the walk continues into every other subtree.
    """
    candidates: set[Path] = set()
    warnings: list[str] = []

    def _on_error(exc: OSError) -> None:
        warnings.append(f"could not scan {exc.filename}: {exc.strerror or exc}")

    for dirpath, _dirnames, filenames in os.walk(root, onerror=_on_error):
        current = Path(dirpath)
        if "agglomeration_manifest.json" in filenames:
            candidates.add(current)
        if current.name.endswith("_agglo"):
            candidates.add(current)
    if root.name.endswith("_agglo"):
        candidates.add(root)
    return candidates, warnings


def _campaign_manifests(root: Path) -> tuple[list[tuple[Path, dict]], list[str]]:
    candidates, warnings = _discover_campaign_dirs(root)
    campaigns: list[tuple[Path, dict]] = []
    for campaign in sorted(candidates):
        manifest_path = campaign / "agglomeration_manifest.json"
        try:
            manifest_text = manifest_path.read_text(encoding="utf-8")
        except OSError:
            campaigns.append((campaign, {"_manifest_error": "manifest is missing"}))
            continue
        except UnicodeDecodeError as exc:
            campaigns.append(
                (campaign, {"_manifest_error": f"manifest is not valid UTF-8: {exc}"})
            )
            continue
        try:
            manifest = json.loads(manifest_text)
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
    return campaigns, warnings


def _campaign_label(root: Path, campaign: Path) -> str:
    try:
        relative = campaign.relative_to(root)
    except ValueError:
        return str(campaign)
    return "." if str(relative) == "." else str(relative)


def _case_rows(progress: dict, campaign_label: str, campaign_path: str) -> list[dict]:
    cases = progress.get("cases")
    if not isinstance(cases, list):
        return []
    rows = []
    for case in cases:
        if not isinstance(case, dict):
            continue
        rows.append(
            {
                "campaign": campaign_label,
                "campaign_path": campaign_path,
                "array_index": case.get("array_index"),
                "case": case.get("case"),
                "status": str(case.get("status", "UNKNOWN")),
                "charge": case.get("charge"),
                "uhf": case.get("uhf"),
                "detail": str(case.get("detail", "")),
            }
        )
    return rows


def _cross_check_progress(progress: dict) -> str | None:
    """Sanity-check a campaign's self-reported totals against its own rows.

    audit_xtb_runs.py computes counts_by_status/total/completed/pending from
    the same "cases" list it writes into xtb_progress.json, so these must
    always agree. A mismatch means the file was hand-edited, written by an
    incompatible version of the script, or truncated by a crash mid-write --
    all worth surfacing rather than silently trusting the summary numbers.
    """
    cases = progress.get("cases")
    if not isinstance(cases, list) or not cases:
        return None
    actual_counts = Counter(
        str(case.get("status", "UNKNOWN")) for case in cases if isinstance(case, dict)
    )
    reported_counts = progress.get("counts_by_status", {})
    if not isinstance(reported_counts, dict):
        reported_counts = {}
    reported_counts_normalized = {str(k): _as_int(v) for k, v in reported_counts.items()}
    if dict(actual_counts) != reported_counts_normalized:
        return (
            "counts_by_status does not match the per-case rows in the same file "
            f"(reported {reported_counts_normalized}, computed {dict(actual_counts)})"
        )
    reported_total = progress.get("total")
    if reported_total is not None and _as_int(reported_total, -1) != len(cases):
        return (
            f"total ({reported_total}) does not match the number of case rows "
            f"({len(cases)})"
        )
    return None


def _audit_campaign(
    root: Path, campaign: Path, manifest: dict, timeout: float
) -> tuple[dict[str, object], list[dict]]:
    label = _campaign_label(root, campaign)
    row: dict[str, object] = {
        "campaign": label,
        "path": str(campaign.resolve()),
        "campaign_status": "ERROR",
        "manifest_status": str(manifest.get("status", "UNKNOWN")),
        "total": 0,
        "completed": 0,
        "pending": 0,
        "attention": 0,
        "percent_complete": 0.0,
        "counts_by_status": {},
        "data_consistent": True,
        "audit_seconds": None,
        "detail": ERROR_DETAIL_DEFAULT,
    }
    if manifest.get("_manifest_error"):
        row["detail"] = str(manifest["_manifest_error"])
        return row, []

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
        return row, []
    if not audit.is_file():
        row["detail"] = (
            "audit_xtb_runs.py is missing; regenerate this campaign's xTB "
            "launcher with --regenerate-xtb-launcher"
        )
        return row, []

    started = time.monotonic()
    try:
        result = subprocess.run(
            [sys.executable, str(audit), "--summary-only"],
            cwd=campaign,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        row["audit_seconds"] = round(time.monotonic() - started, 2)
        row["detail"] = f"audit_xtb_runs.py did not finish within {timeout:.0f}s"
        return row, []
    except OSError as exc:
        row["detail"] = f"could not launch audit_xtb_runs.py: {exc}"
        return row, []
    row["audit_seconds"] = round(time.monotonic() - started, 2)

    if result.returncode:
        row["detail"] = _truncate(
            f"audit failed with exit code {result.returncode}: "
            f"{(result.stderr or result.stdout)}"
        )
        return row, []
    progress_path = campaign / "xtb_progress.json"
    try:
        progress_text = progress_path.read_text(encoding="utf-8")
    except OSError as exc:
        row["detail"] = f"audit did not produce a readable xtb_progress.json: {exc}"
        return row, []
    try:
        progress = json.loads(progress_text)
    except json.JSONDecodeError as exc:
        row["detail"] = f"xtb_progress.json is not valid JSON: {exc}"
        return row, []
    if not isinstance(progress, dict):
        row["detail"] = "xtb_progress.json did not contain a JSON object"
        return row, []

    inconsistency = _cross_check_progress(progress)

    counts = progress.get("counts_by_status", {})
    if not isinstance(counts, dict):
        counts = {}
    total = _as_int(progress.get("total"))
    completed = _as_int(progress.get("completed"))
    pending = _as_int(progress.get("pending"), max(0, total - completed))
    attention = sum(_as_int(counts.get(status)) for status in ATTENTION_STATES)
    row.update(
        total=total,
        completed=completed,
        pending=pending,
        attention=attention,
        percent_complete=round(100.0 if total == 0 else 100.0 * completed / total, 1),
        counts_by_status=dict(sorted(counts.items())),
    )

    cases = _case_rows(progress, label, str(row["path"]))

    if inconsistency:
        row["campaign_status"] = "ATTENTION"
        row["data_consistent"] = False
        row["detail"] = f"inconsistent audit output: {inconsistency}"
    elif attention:
        row["campaign_status"] = "ATTENTION"
        row["detail"] = (
            f"{attention} case(s) have inconsistent/missing completion data; "
            "see the Attention rows/sheet for the exact case IDs"
        )
    elif pending:
        row["campaign_status"] = "IN_PROGRESS"
        row["detail"] = f"{pending} case(s) still pending"
    else:
        row["campaign_status"] = "COMPLETE"
        row["detail"] = "all xTB cases are safely complete"
    return row, cases


def _write_csv_reports(
    prefix: Path,
    rows: list[dict[str, object]],
    all_cases: list[dict[str, object]],
    attention_cases: list[dict[str, object]],
) -> dict[str, Path]:
    csv_path = prefix.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CAMPAIGN_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            serialized = dict(row)
            serialized["counts_by_status"] = json.dumps(
                serialized.get("counts_by_status", {}), sort_keys=True
            )
            writer.writerow(serialized)

    cases_path = prefix.parent / f"{prefix.name}_cases.csv"
    with cases_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CASE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_cases)

    attention_path = prefix.parent / f"{prefix.name}_attention.csv"
    with attention_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CASE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(attention_cases)

    return {"campaigns": csv_path, "cases": cases_path, "attention": attention_path}


def _write_xlsx_report(prefix: Path, report: dict) -> Path | None:
    """Write a multi-sheet workbook mirroring the CSV/JSON reports.

    Best-effort: openpyxl is only declared under the project's "analysis"
    extra, and this audit is a core CLI feature that must keep working
    (via CSV/JSON) even where that extra isn't installed.
    """
    try:
        from .analysis.report import _openpyxl, _style_table_sheet
    except ImportError:
        return None
    try:
        api = _openpyxl()
    except RuntimeError as exc:
        print(f"Note: skipping the multi-sheet workbook ({exc})")
        return None

    workbook = api["Workbook"]()
    workbook.properties.creator = "nio-md-prep"
    workbook.properties.title = "NiO agglomeration xTB audit"

    campaign_headers = [
        "Campaign", "Status", "Manifest Status", "Total", "Completed", "Pending",
        "Attention", "% Complete", "Data Consistent", "Audit Seconds", "Detail", "Path",
    ]
    campaign_rows = [
        [
            row["campaign"], row["campaign_status"], row["manifest_status"], row["total"],
            row["completed"], row["pending"], row["attention"], row["percent_complete"],
            row.get("data_consistent", True), row.get("audit_seconds"), row["detail"],
            row["path"],
        ]
        for row in report["campaign_results"]
    ]
    sheet = workbook.active
    sheet.title = "Campaigns"
    _style_table_sheet(sheet, campaign_headers, campaign_rows, "CampaignsTable", api)

    case_headers = ["Campaign", "Array Index", "Case", "Status", "Charge", "UHF", "Detail"]

    def _rows_for(cases: list[dict]) -> list[list]:
        return [
            [c["campaign"], c["array_index"], c["case"], c["status"], c["charge"], c["uhf"], c["detail"]]
            for c in cases
        ]

    for title, key, table_name in (
        ("All Cases", "case_rows", "AllCasesTable"),
        ("Attention", "attention_cases", "AttentionTable"),
        ("Pending", "pending_cases", "PendingTable"),
    ):
        ws = workbook.create_sheet(title)
        _style_table_sheet(ws, case_headers, _rows_for(report[key]), table_name, api)

    status_sheet = workbook.create_sheet("Status Counts")
    _style_table_sheet(
        status_sheet,
        ["Status", "Case Count", "Needs Attention"],
        [
            [status, count, status in ATTENTION_STATES]
            for status, count in sorted(report["case_status_counts"].items())
        ],
        "StatusCountsTable",
        api,
    )

    issues_sheet = workbook.create_sheet("Campaign Issues")
    _style_table_sheet(
        issues_sheet,
        ["Campaign", "Status", "Detail", "Path"],
        [
            [row["campaign"], row["campaign_status"], row["detail"], row["path"]]
            for row in report["campaign_issues"]
        ],
        "CampaignIssuesTable",
        api,
    )

    if report.get("scan_warnings"):
        warn_sheet = workbook.create_sheet("Scan Warnings")
        _style_table_sheet(
            warn_sheet,
            ["Warning"],
            [[warning] for warning in report["scan_warnings"]],
            "ScanWarningsTable",
            api,
        )

    xlsx_path = prefix.with_suffix(".xlsx")
    workbook.save(xlsx_path)
    return xlsx_path


def _print_case_sample(cases: list[dict], limit: int = 8) -> None:
    for case in cases[:limit]:
        print(f"    [{case['status']}] {case['case']} - {case['detail']}")
    remaining = len(cases) - limit
    if remaining > 0:
        print(f"    ... and {remaining} more (see the case-level CSV/workbook for the full list)")


def _print_console_report(report: dict, written_paths: list[Path]) -> None:
    for warning in report.get("scan_warnings", []):
        print(f"WARNING: {warning}")

    print("Campaign                         Status           Complete   Pending  Attention")
    print("-------------------------------- ---------------- ---------- --------- ---------")
    for row in report["campaign_results"]:
        print(
            f"{str(row['campaign'])[:32]:32} "
            f"{str(row['campaign_status']):16} "
            f"{int(row['completed']):4}/{int(row['total']):<5} "
            f"{int(row['pending']):9} {int(row['attention']):9}"
        )

    by_campaign: dict[str, list[dict]] = {}
    for case in report["attention_cases"] + report["pending_cases"]:
        by_campaign.setdefault(str(case["campaign"]), []).append(case)

    for row in report["campaign_results"]:
        campaign = str(row["campaign"])
        if row["campaign_status"] == "ERROR":
            print(f"ERROR {campaign} ({row['path']}): {row['detail']}")
        elif row["campaign_status"] in {"PACKMOL_REQUIRED", "XTB_DISABLED"}:
            print(f"BLOCKED {campaign}: {row['detail']}")
        elif row["campaign_status"] == "ATTENTION":
            print(f"ATTENTION {campaign}: {row['detail']}")
            _print_case_sample(
                [c for c in by_campaign.get(campaign, []) if c["status"] in ATTENTION_STATES]
            )
        elif row["campaign_status"] == "IN_PROGRESS":
            pending_here = [c for c in by_campaign.get(campaign, []) if c["status"] not in ATTENTION_STATES]
            if pending_here:
                print(f"IN_PROGRESS {campaign}: {row['detail']}")
                _print_case_sample(pending_here)

    print(
        f"Overall xTB progress: {report['completed']}/{report['cases']} safely complete; "
        f"{report['pending']} pending or requiring attention across "
        f"{report['campaigns']} campaign(s)."
    )
    error_campaigns = sum(
        1 for row in report["campaign_results"] if row["campaign_status"] == "ERROR"
    )
    if error_campaigns:
        print(
            f"{error_campaigns} campaign(s) could not be audited at all -- see the "
            "ERROR line(s) above and the Campaign Issues sheet/rows before trusting "
            "the totals."
        )
    print("Reports written: " + ", ".join(str(path) for path in written_paths))


def audit_agglomerations(
    root: Path,
    output_prefix: Path | None = None,
    *,
    timeout: float = DEFAULT_AUDIT_TIMEOUT_SECONDS,
) -> dict:
    root = root.resolve()
    if not root.is_dir():
        raise ValueError(f"agglomeration audit root is not a directory: {root}")
    campaigns, scan_warnings = _campaign_manifests(root)
    if not campaigns:
        raise ValueError(f"no agglomeration campaigns found under {root}")

    rows: list[dict[str, object]] = []
    all_cases: list[dict[str, object]] = []
    for campaign, manifest in campaigns:
        row, cases = _audit_campaign(root, campaign, manifest, timeout)
        rows.append(row)
        all_cases.extend(cases)

    attention_cases = [c for c in all_cases if c["status"] in ATTENTION_STATES]
    pending_cases = [
        c for c in all_cases
        if c["status"] != "COMPLETE" and c["status"] not in ATTENTION_STATES
    ]
    case_status_counts = Counter(str(c["status"]) for c in all_cases)
    campaign_status_counts = Counter(str(row["campaign_status"]) for row in rows)
    campaign_issues = [
        row for row in rows
        if row["campaign_status"] in {"ERROR", "PACKMOL_REQUIRED", "XTB_DISABLED"}
    ]

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
        "scan_warnings": scan_warnings,
        **totals,
        "campaign_counts_by_status": dict(sorted(campaign_status_counts.items())),
        "case_status_counts": dict(sorted(case_status_counts.items())),
        "campaign_results": rows,
        "case_rows": all_cases,
        "attention_cases": attention_cases,
        "pending_cases": pending_cases,
        "campaign_issues": campaign_issues,
    }

    prefix = (output_prefix or (root / "agglomeration_xtb_audit")).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_paths = _write_csv_reports(prefix, rows, all_cases, attention_cases)
    json_path = prefix.with_suffix(".json")
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    xlsx_path = _write_xlsx_report(prefix, report)

    written_paths = [csv_paths["campaigns"], csv_paths["cases"], csv_paths["attention"], json_path]
    if xlsx_path is not None:
        written_paths.append(xlsx_path)
    _print_console_report(report, written_paths)
    return report
