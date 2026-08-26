import json
from pathlib import Path

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
                {"array_index": 0, "case": "xtb/case0", "status": "COMPLETE", "charge": "0", "uhf": "0", "detail": "ok"},
                {"array_index": 1, "case": "xtb/case1", "status": "COMPLETE", "charge": "0", "uhf": "0", "detail": "ok"},
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
                {"array_index": 0, "case": "xtb/case0", "status": "COMPLETE", "charge": "0", "uhf": "0", "detail": "ok"},
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
    assert (tmp_path / "agglomeration_xtb_audit_attention.csv").is_file()

    # The per-case detail that was previously only inside each campaign's own
    # xtb_progress.json must now be surfaced in the consolidated report, so
    # you can tell exactly which case is still pending without hunting
    # through every campaign directory.
    assert len(report["case_rows"]) == 5
    pending_cases = {c["case"] for c in report["pending_cases"]}
    assert pending_cases == {"xtb/case1", "xtb/case2"}

    out = capsys.readouterr().out
    assert "Overall xTB progress: 3/5 safely complete" in out
    assert "xtb/case1" in out or "IN_PROGRESS mixed_agglo" in out


def test_audit_agglomerations_flags_missing_auditor_as_error(tmp_path):
    _campaign(tmp_path, "legacy_agglo", None)
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    # This campaign never even produced case-level data, which is a
    # different situation from "cases exist and need attention" -- it gets
    # its own ERROR status so the two are never confused.
    assert row["campaign_status"] == "ERROR"
    assert "regenerate" in row["detail"]
    assert row in report["campaign_issues"]


def test_audit_agglomerations_flags_agglo_folder_without_manifest(tmp_path):
    (tmp_path / "broken_agglo").mkdir()
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    assert row["campaign_status"] == "ERROR"
    assert row["detail"] == "manifest is missing"
    assert row in report["campaign_issues"]


def test_audit_agglomerations_flags_attention_cases_with_case_detail(tmp_path):
    _campaign(
        tmp_path,
        "attention_agglo",
        {
            "total": 2,
            "completed": 1,
            "pending": 0,
            "counts_by_status": {"COMPLETE": 1, "HASH_MISMATCH": 1},
            "cases": [
                {"array_index": 0, "case": "xtb/case0", "status": "COMPLETE", "charge": "0", "uhf": "0", "detail": "ok"},
                {
                    "array_index": 1,
                    "case": "xtb/case1",
                    "status": "HASH_MISMATCH",
                    "charge": "0",
                    "uhf": "0",
                    "detail": "final geometry present but completion hash differs",
                },
            ],
        },
    )
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    assert row["campaign_status"] == "ATTENTION"
    assert len(report["attention_cases"]) == 1
    assert report["attention_cases"][0]["case"] == "xtb/case1"
    assert report["attention_cases"][0]["status"] == "HASH_MISMATCH"


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
                {"array_index": 0, "case": "xtb/case0", "status": "PENDING", "charge": "0", "uhf": "0", "detail": "not started"},
            ],
        },
    )
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    assert row["campaign_status"] == "ATTENTION"
    assert row["data_consistent"] is False
    assert "inconsistent audit output" in row["detail"]


def test_audit_agglomerations_strict_exit_code_via_cli(tmp_path, monkeypatch, capsys):
    from nio_md_prep import cli

    _campaign(tmp_path, "legacy_agglo", None)
    monkeypatch.chdir(tmp_path)
    rc = cli.main(["audit-agglomerations", "--strict"])
    assert rc == 1
    capsys.readouterr()
