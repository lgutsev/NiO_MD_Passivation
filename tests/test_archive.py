from __future__ import annotations

import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import pytest

from nio_md_prep.archive import ARCHIVED_SUFFIXES, create_run_archive


def test_create_run_archive_collects_rerun_files_and_excludes_bulk_outputs(tmp_path):
    root_a=tmp_path/"prepared-rerun"
    root_b=tmp_path/"prepared-rerun2"
    system_a=root_a/"me-4pacz-alone"
    system_b=root_b/"me-4pacz-dcz-4p-cosam"
    system_a.mkdir(parents=True)
    system_b.mkdir(parents=True)

    included={
        system_a/"deposited.data":b"LAMMPS state\n",
        system_a/"hold-300K.in":b"read_data deposited.data\n",
        system_a/"topology_output.lmp":b"LAMMPS topology\n",
        system_a/"packmol.inp":b"seed 202405367\n",
        system_a/"resolved_config.toml":b"[study]\n",
        system_a/"assembly_manifest.json":b"{}\n",
        system_a/"protocol_notes.txt":b"protocol\n",
        system_a/"packed.xyz":b"0\npacked\n",
        system_b/"held-400K.data":b"held state\n",
    }
    excluded={
        system_a/"hold-300K.lammpstrj":b"trajectory\n",
        system_a/"log.hold-300K.lammps":b"log\n",
        system_a/"held-300K.restart":b"restart\n",
        system_a/"coverage_summary.xlsx":b"workbook\n",
        system_a/"nio.hold300.out":b"slurm\n",
    }
    for path,content in (included|excluded).items():
        path.write_bytes(content)

    output=tmp_path/"reruns.zip"
    result=create_run_archive(
        [root_a,root_b],
        output,
        provenance_directory=tmp_path,
    )

    assert result.path==output.resolve()
    assert result.file_count==len(included)
    assert result.uncompressed_bytes==sum(len(content) for content in included.values())
    with ZipFile(output) as archive:
        names=set(archive.namelist())
        assert "README_RERUN.txt" in names
        assert "ARCHIVE_MANIFEST.json" in names
        expected_names={
            (Path(path.parents[1].name)/path.parent.name/path.name).as_posix()
            for path in included
        }
        assert expected_names <= names
        assert not {
            (Path(path.parents[1].name)/path.parent.name/path.name).as_posix()
            for path in excluded
        } & names
        manifest=json.loads(archive.read("ARCHIVE_MANIFEST.json"))

    assert manifest["format"]=="nio-md-prep-rerun-archive"
    assert manifest["source_roots"]==["prepared-rerun","prepared-rerun2"]
    assert set(manifest["selection"]["suffixes"])==ARCHIVED_SUFFIXES
    records={item["path"]:item for item in manifest["files"]}
    for path,content in included.items():
        archive_name=(Path(path.parents[1].name)/path.parent.name/path.name).as_posix()
        assert records[archive_name]["size_bytes"]==len(content)
        assert records[archive_name]["sha256"]==hashlib.sha256(content).hexdigest()


def test_create_run_archive_refuses_overwrite_without_force(tmp_path):
    root=tmp_path/"prepared"
    root.mkdir()
    (root/"deposited.data").write_text("state\n")
    output=tmp_path/"reruns.zip"
    output.write_text("existing")

    with pytest.raises(FileExistsError,match="--force"):
        create_run_archive([root],output,provenance_directory=tmp_path)

    result=create_run_archive(
        [root],
        output,
        force=True,
        provenance_directory=tmp_path,
    )
    assert result.file_count==1
    assert ZipFile(output).testzip() is None


def test_create_run_archive_does_not_follow_symlinked_directories(tmp_path):
    root=tmp_path/"prepared"
    system=root/"me-4pacz-alone"
    system.mkdir(parents=True)
    (system/"deposited.data").write_text("state\n")

    outside=tmp_path/"outside-secret"
    outside.mkdir()
    (outside/"secret.data").write_text("should not be archived\n")

    link=system/"linked"
    try:
        link.symlink_to(outside,target_is_directory=True)
    except OSError:
        pytest.skip("creating directory symlinks is not permitted in this environment")

    output=tmp_path/"reruns.zip"
    result=create_run_archive([root],output,provenance_directory=tmp_path)

    assert result.file_count==1
    with ZipFile(output) as archive:
        names=set(archive.namelist())
    assert not any("secret.data" in name for name in names)


def test_create_run_archive_requires_unique_root_names(tmp_path):
    first=tmp_path/"one"/"prepared"
    second=tmp_path/"two"/"prepared"
    first.mkdir(parents=True)
    second.mkdir(parents=True)
    (first/"one.data").write_text("one\n")
    (second/"two.data").write_text("two\n")

    with pytest.raises(ValueError,match="unique directory names"):
        create_run_archive(
            [first,second],
            tmp_path/"reruns.zip",
            provenance_directory=tmp_path,
        )
