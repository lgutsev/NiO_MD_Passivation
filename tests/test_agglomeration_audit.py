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
        {"total": 2, "completed": 2, "pending": 0, "counts_by_status": {"COMPLETE": 2}},
    )
    _campaign(
        tmp_path,
        "mixed_agglo",
        {
            "total": 3,
            "completed": 1,
            "pending": 2,
            "counts_by_status": {"COMPLETE": 1, "PENDING": 2},
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
    assert "Overall xTB progress: 3/5 safely complete" in capsys.readouterr().out


def test_audit_agglomerations_flags_missing_auditor(tmp_path):
    _campaign(tmp_path, "legacy_agglo", None)
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    assert row["campaign_status"] == "ATTENTION"
    assert "regenerate" in row["detail"]


def test_audit_agglomerations_flags_agglo_folder_without_manifest(tmp_path):
    (tmp_path / "broken_agglo").mkdir()
    report = audit_agglomerations(tmp_path)
    row = report["campaign_results"][0]
    assert row["campaign_status"] == "ATTENTION"
    assert row["detail"] == "manifest is missing"
