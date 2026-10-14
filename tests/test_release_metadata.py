"""Pin that the release version is declared consistently across metadata files."""
import re
from pathlib import Path

import nio_md_prep

ROOT = Path(__file__).resolve().parents[1]


def _first_match(path, pattern):
    match = re.search(pattern, (ROOT / path).read_text(encoding="utf-8"), re.MULTILINE)
    assert match, f"no version found in {path}"
    return match.group(1)


def test_version_matches_pyproject():
    assert _first_match("pyproject.toml", r'^version\s*=\s*"([^"]+)"') == nio_md_prep.__version__


def test_version_matches_citation():
    assert _first_match("CITATION.cff", r"^version:\s*(\S+)") == nio_md_prep.__version__
