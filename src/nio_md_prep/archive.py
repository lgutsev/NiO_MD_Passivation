from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Iterable, Sequence
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


ARCHIVED_SUFFIXES = frozenset(
    {
        ".data",
        ".in",
        ".inp",
        ".json",
        ".lmp",
        ".toml",
        ".txt",
        ".xyz",
    }
)

ARCHIVE_README = """NiO_MD_Passivation rerun archive
=====================================

This archive contains restart-critical LAMMPS state and input files collected
from the requested prepared/work roots. Directory structure is preserved.

Included by default:
  *.data  LAMMPS states written at each completed simulation stage
  *.in    LAMMPS stage inputs
  *.lmp   assembled topology and force-field inputs
  *.inp   Packmol inputs
  *.toml, *.json, *.txt, *.xyz  resolved configuration and build provenance

Excluded by design:
  trajectories, LAMMPS logs, Slurm output, restart binaries, and analysis files

After extraction, run a stage from its own system directory so its relative
read_data and include paths resolve correctly. For example:

  cd prepared-rerun/me-4pacz-alone
  srun -n 64 lmp_mpi -in hold-300K.in

ARCHIVE_MANIFEST.json lists every archived source file, its uncompressed size,
and its SHA-256 checksum. It also records the source repository commit when
that information was available at archive time.
"""


@dataclass(frozen=True)
class ArchiveResult:
    path: Path
    file_count: int
    uncompressed_bytes: int


def _git_provenance(directory: Path) -> dict[str, object]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=directory,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=directory,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return {"commit": None, "working_tree_dirty": None}
    return {"commit": commit, "working_tree_dirty": dirty}


def _selected_files(root: Path) -> Iterable[Path]:
    for path in sorted(root.rglob("*")):
        if (
            path.is_file()
            and not path.is_symlink()
            and path.suffix.lower() in ARCHIVED_SUFFIXES
        ):
            yield path


def _write_source_file(
    archive: ZipFile,
    source: Path,
    archive_name: str,
) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    info = ZipInfo.from_file(source, arcname=archive_name)
    info.compress_type = ZIP_DEFLATED
    with source.open("rb") as input_file, archive.open(
        info,
        mode="w",
        force_zip64=True,
    ) as output_file:
        while chunk := input_file.read(1024 * 1024):
            digest.update(chunk)
            size += len(chunk)
            output_file.write(chunk)
    return digest.hexdigest(), size


def create_run_archive(
    roots: Sequence[Path],
    output: Path,
    *,
    force: bool = False,
    provenance_directory: Path | None = None,
) -> ArchiveResult:
    """Archive restart-critical files from one or more prepared/work roots."""
    if not roots:
        raise ValueError("at least one prepared/work root is required")

    resolved_roots: list[tuple[str, Path]] = []
    labels: set[str] = set()
    for supplied_root in roots:
        root = supplied_root.expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"prepared/work root does not exist: {supplied_root}")
        label = root.name
        if label in labels:
            raise ValueError(
                f"archive roots must have unique directory names; duplicate: {label}"
            )
        labels.add(label)
        resolved_roots.append((label, root))

    output = output.expanduser().resolve()
    if output.suffix.lower() != ".zip":
        raise ValueError("archive output must end in .zip")
    if output.exists() and not force:
        raise FileExistsError(
            f"archive already exists: {output}; use --force to replace it"
        )

    selected: list[tuple[Path, str]] = []
    for label, root in resolved_roots:
        for source in _selected_files(root):
            selected.append((source, (Path(label) / source.relative_to(root)).as_posix()))
    selected.sort(key=lambda item: item[1])
    if not selected:
        suffixes = ", ".join(sorted(ARCHIVED_SUFFIXES))
        raise ValueError(f"no archiveable files found; selected suffixes are: {suffixes}")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    if temporary.exists():
        temporary.unlink()

    manifest_files: list[dict[str, object]] = []
    total_bytes = 0
    try:
        with ZipFile(
            temporary,
            mode="w",
            compression=ZIP_DEFLATED,
            compresslevel=6,
            allowZip64=True,
        ) as archive:
            archive.writestr("README_RERUN.txt", ARCHIVE_README)
            for source, archive_name in selected:
                checksum, size = _write_source_file(archive, source, archive_name)
                manifest_files.append(
                    {
                        "path": archive_name,
                        "size_bytes": size,
                        "sha256": checksum,
                    }
                )
                total_bytes += size

            provenance_root = (
                provenance_directory.resolve()
                if provenance_directory is not None
                else Path.cwd().resolve()
            )
            manifest = {
                "format": "nio-md-prep-rerun-archive",
                "format_version": 1,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "repository": _git_provenance(provenance_root),
                "source_roots": [label for label, _ in resolved_roots],
                "selection": {
                    "suffixes": sorted(ARCHIVED_SUFFIXES),
                    "symlinks": "excluded",
                },
                "file_count": len(manifest_files),
                "uncompressed_bytes": total_bytes,
                "files": manifest_files,
            }
            archive.writestr(
                "ARCHIVE_MANIFEST.json",
                json.dumps(manifest, indent=2) + "\n",
            )
        temporary.replace(output)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise

    return ArchiveResult(
        path=output,
        file_count=len(manifest_files),
        uncompressed_bytes=total_bytes,
    )
