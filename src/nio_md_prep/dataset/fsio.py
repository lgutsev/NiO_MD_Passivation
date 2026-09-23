"""Hashing, atomic writes and output-directory safety for dataset commands."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Any, Iterable

from .errors import OverwriteRefusedError

_CHUNK = 1 << 20


def sha256_file(path: Path) -> str:
    """SHA-256 of a file's bytes, streamed so multi-GB outputs never load at once."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    """Deterministic JSON: sorted keys, no whitespace variation, no NaN."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _temporary_sibling(path: Path) -> Path:
    return path.with_name(f".{path.name}.tmp-{os.getpid()}")


def atomic_write_bytes(path: Path, data: bytes) -> None:
    """Write ``data`` to ``path`` so readers never see a partial file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_sibling(path)
    try:
        with temporary.open("wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_text(path: Path, text: str) -> None:
    # newline="\n" semantics: encode explicitly so Windows never rewrites line endings.
    atomic_write_bytes(path, text.encode("utf-8"))


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


class AtomicFile:
    """Context manager: stream text to a hidden temporary, publish on success only.

    On an exception the temporary is removed and ``path`` is left untouched.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self.temporary = _temporary_sibling(self.path)
        self._handle = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.temporary.open("w", encoding="utf-8", newline="\n")
        return self._handle

    def __exit__(self, exc_type, exc, tb):
        assert self._handle is not None
        self._handle.flush()
        os.fsync(self._handle.fileno())
        self._handle.close()
        if exc_type is None:
            os.replace(self.temporary, self.path)
        elif self.temporary.exists():
            self.temporary.unlink()
        return False


def write_jsonl(path: Path, records: Iterable[dict]) -> int:
    """Atomically write one canonical-JSON record per line; returns the count."""
    count = 0
    with AtomicFile(path) as handle:
        for record in records:
            handle.write(canonical_json(record))
            handle.write("\n")
            count += 1
    return count


def read_jsonl(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def prepare_output_dir(path: Path, *, force: bool) -> Path:
    """Create ``path`` or refuse to write into an existing non-empty directory.

    ``force`` allows writing into a non-empty directory; files are then
    replaced one by one, atomically, and nothing unrelated is deleted.
    """
    path = Path(path)
    if path.exists() and not path.is_dir():
        raise OverwriteRefusedError(f"{path} exists and is not a directory")
    if path.is_dir() and any(path.iterdir()) and not force:
        raise OverwriteRefusedError(
            f"{path} already contains files; choose a new output directory or pass --force"
        )
    path.mkdir(parents=True, exist_ok=True)
    return path


def git_provenance(directory: Path | None = None) -> dict[str, Any]:
    """Commit of the code that produced a dataset.

    Untracked files are ignored (``--untracked-files=no``): generated data
    trees next to a checkout must not make every export look "dirty".
    """
    directory = Path(directory) if directory is not None else Path(__file__).resolve().parent
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=directory, check=True,
            capture_output=True, text=True, timeout=20,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain", "--untracked-files=no"], cwd=directory,
                check=True, capture_output=True, text=True, timeout=20,
            ).stdout.strip()
        )
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "working_tree_dirty": None}
    return {"commit": commit, "working_tree_dirty": dirty}
