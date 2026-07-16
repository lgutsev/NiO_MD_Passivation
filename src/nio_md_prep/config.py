from __future__ import annotations
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[2]

def load(path: Path) -> dict:
    with path.open("rb") as f: return tomllib.load(f)

def molecule_dir(slug: str) -> Path:
    return ROOT / "inputs" / "molecules" / slug

def molecule_manifest(slug: str) -> tuple[Path, dict]:
    folder = molecule_dir(slug); path = folder / "molecule.toml"
    if not path.exists(): raise FileNotFoundError(f"molecule manifest not found: {path}")
    return folder, load(path)

def missing_ligpargen(path: Path, rerun: str) -> FileNotFoundError:
    return FileNotFoundError(
        f"LigParGen input is missing.\nExpected path: {path.resolve()}\n"
        "Run LigParGen independently (outside this package), download its LAMMPS .lmp output, "
        f"and save it as ligpargen.lmp at that path.\nRerun: {rerun}"
    )
