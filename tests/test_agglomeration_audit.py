import json
from pathlib import Path

from openpyxl import load_workbook

from nio_md_prep.agglomeration_audit import audit_agglomerations


def _campaign(root: Path, name: str, progress: dict | None) -> Path:
    campaign = root / name
    campaign.mkdir()
    (campaign / "agglomeration_manifest.json").write_text(
        json.dumps({"status": "XTB_REQUIRED", "xtb": {"enabled": True}, "replicas": []}),
        encoding="utf-8",
    )
    (campaign / "xtb_cases.tsv").write_text("xtb/case\t0\t0\n", encoding="utf-8")
    if progress is not None:
        (campaign / "audit_xtb_runs.py").write_text(
            "#!/usr/bin/env python3\n"
            "import json\n"
            f"data = {progress!r}\n"
            "open('xtb_progress.json', 'w').write(json.dumps(data))\n",
            encoding="utf-8",
        )
    return campaign


def test_audit_agglomerations_consolidates_campaigns(tmp_path, capsys):
    _campaign(
        tmp_path,
        "Me4PACz_agglo",
        {
            "total": 2,
            "completed": 2,
            "pending": 0,
            "counts_by_status": {"COMPLETE": 2},
            "cases": [
                {"array_index": 0, "case": "xtb/case0", "status": "COMPLETE", "charge": "0", "uhf": "0", "detail": "done"},
                {"array_index": 1, "case": "xtb/case1", "status": "COMPLETE", "charge": "0", "uhf": "0", "detail": "done"},
            ],
        },
    )
    _campaign(
        tmp_path,
        "mixed_agglo",
        {
            "total": 3,
            "completed": 1,
            "pending": 2,
            "counts_by_status": {"COMPLETE": 1, "PENDING": 2},
            "cases": [
                {"array_index": 0, "case": "xtb/case0", "status": "COMPLETE", "charge": "0", "uhf": "0", "detail": "done"},
                {"array_index": 1, "case": "xtb/case1", "status": "PENDING", "charge": "0", "uhf": "0", "detail": "not started"},
                {"array_index": 2, "case": "xtb/case2", "status": "PENDING", "charge": "0", "uhf": "0", "detail": "not started"},
            ],
        },
    )
    report = audit_agglomerations(tmp_path)
    assert report["campaigns"] == 2
    assert report["cases"] == 5
    assert report["completed"] == 3
    assert report["pending"] == 2
    assert report["campaign_counts_by_status"] == {"COMPLETE": 1, "IN_PROGRESS": 1}
    assert (tmp_path / "agglomeration_xtb_audit.csv").is_file()
    assert (tmp_path / "agglomeration_xtb_audit.json").is_file()
    assert (tmp_path / "agglomeration_xtb_audit_cases.csv").is_file()
    assert (tmp_path / "agglomeration_xtb_audit_incomplete.csv").is_file()
    workbook_path = tmp_path / "agglomeration_xtb_audit.xlsx"
    assert workbook_path.is_file()
    assert load_workbook(workbook_path, read_only=True).sheetnames == [
        "Overview",
        "Campaigns",
        "Incomplete Cases",
        "Attention",
        "Pending",
        "All Cases",
        "Case Integrity",
        "Artifacts",
        "Logs",
        "Protocols",
        "Resources",
        "Slurm Jobs",
        "Status Counts",
        "Campaign Issues",
        "Scan Warnings",
        "Audit Metadata",
        "Status Definitions",
    ]
    assert len(report["incomplete_cases"]) == 2
    assert "Overall xTB progress: 3/5 safely complete" in capsys.readouterr().out


def test_case_diagnostics_identify_exact_checkpoint_and_log(tmp_path):
    campaign = _campaign(
        tmp_path,
        "diagnostic_agglo",
        {
            "total": 1,
            "completed": 0,
            "pending": 1,
            "counts_by_status": {"STEERING_COMPLETE": 1},
            "cases": [
                {
                    "array_index": 4,
                    "case": "xtb/n08/r000_s00_1p000",
                    "status": "STEERING_COMPLETE",
                    "charge": "0",
                    "uhf": "0",
                    "detail": "unbiased MD and final optimization remain",
                }
            ],
        },
    )
    work = campaign / "xtb/n08/r000_s00_1p000"
    work.mkdir(parents=True)
    (work / "xtbopt_initial.xyz").write_text("initial\n", encoding="utf-8")
    (work / "xtbsteered_last.xyz").write_text("steered\n", encoding="utf-8")
    (work / "xtb_unbiased_md.err").write_text(
        "DUE TO TIME LIMIT\n", encoding="utf-8"
    )
    (campaign / "xtb-pool.12345.out").write_text(
        ">>> Starting xTB case index 4 at Tue Aug 26 12:00:00 UTC 2026\n",
        encoding="utf-8",
    )
    (campaign / "xtb-pool.12345.err").write_text(
        "<<< Finished xTB case index 4 with ERROR 1 at Tue Aug 26 12:10:00 UTC 2026\n",
        encoding="utf-8",
    )
    report = audit_agglomerations(tmp_path, include_slurm=False)
    case = report["incomplete_cases"][0]
    assert case["campaign"] == "diagnostic_agglo"
    assert case["array_index"] == 4
    assert case["case"] == "xtb/n08/r000_s00_1p000"
    assert case["initial_opt_complete"] is True
    assert case["steered_md_complete"] is True
    assert case["unbiased_md_complete"] is False
    assert case["likely_issue"] == "WALLTIME_EXCEEDED"
    assert case["latest_log"] == "xtb_unbiased_md.err"
    assert case["latest_job_id"] == "12345"
    assert case["worker_state"] == "FAILED"
    assert case["worker_exit_code"] == "1"


def test_audit_agglomerations_flags_missing_auditor(tmp_path):
    _campaign(tmp_path, "legacy_agglo", None)
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    assert row["campaign_status"] == "ERROR"
    assert "regenerate" in row["detail"]


def test_audit_agglomerations_flags_agglo_folder_without_manifest(tmp_path):
    (tmp_path / "broken_agglo").mkdir()
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    assert row["campaign_status"] == "ERROR"
    assert row["detail"] == "manifest is missing"


def test_audit_agglomerations_detects_self_inconsistent_progress(tmp_path):
    _campaign(
        tmp_path,
        "corrupt_agglo",
        {
            "total": 99,
            "completed": 99,
            "pending": 0,
            "counts_by_status": {"COMPLETE": 99},
            "cases": [
                {
                    "array_index": 0,
                    "case": "xtb/case0",
                    "status": "PENDING",
                    "charge": "0",
                    "uhf": "0",
                    "detail": "not started",
                }
            ],
        },
    )
    row = audit_agglomerations(tmp_path)["campaign_results"][0]
    assert row["campaign_status"] == "ATTENTION"
    assert row["data_consistent"] is False
    assert "inconsistent audit output" in row["detail"]


def test_audit_agglomerations_strict_exit_code_via_cli(tmp_path, monkeypatch, capsys):
    from nio_md_prep import cli

    _campaign(tmp_path, "legacy_agglo", None)
    monkeypatch.chdir(tmp_path)
    assert cli.main(["audit-agglomerations", "--strict", "--no-slurm"]) == 1
    capsys.readouterr()
