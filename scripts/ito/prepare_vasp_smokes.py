#!/usr/bin/env python
"""Generate the ITO/SAM VASP reference inputs (smoke or production).

Examples (from the repository root):

    python scripts/ito/prepare_vasp_smokes.py                       # smoke set -> inputs/ito/vasp_smoke
    python scripts/ito/prepare_vasp_smokes.py --check               # regenerate in memory, compare, write nothing
    python scripts/ito/prepare_vasp_smokes.py --mode production --trilayers 3 --output inputs/ito/vasp_production

Nothing is submitted: each case directory gets a run.sbatch, and the output
root gets submit_all.sh, for review and manual submission on LONI.
See docs/ito/vasp-smoke.md.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from nio_md_prep.ito.vasp import MOLECULES, Settings, render, write_all


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output", type=Path, default=None, help="output root (default inputs/ito/vasp_smoke or vasp_production)")
    p.add_argument("--mode", choices=("smoke", "production"), default="smoke",
                   help="smoke: NSW=3, 12 h / 4 h walltime; production: NSW=300, 72 h, WAVECAR/CHGCAR kept")
    p.add_argument("--trilayers", type=int, default=2, help="O-In-O trilayers (bottom one fixed); 3 for production")
    p.add_argument("--molecules", nargs="+", default=list(MOLECULES), help="molecule slugs under inputs/molecules")
    p.add_argument("--per-case-cell", action="store_true",
                   help="size c per case instead of one common cell for all cases (not recommended for E_ads)")
    p.add_argument("--check", action="store_true", help="regenerate in memory and compare with the output tree; write nothing")
    a = p.parse_args(argv)
    settings = Settings(mode=a.mode, trilayers=a.trilayers, molecules=tuple(a.molecules), common_cell=not a.per_case_cell)
    out = a.output or ROOT / "inputs" / "ito" / ("vasp_smoke" if a.mode == "smoke" else "vasp_production")
    if a.check:
        files = render(settings)
        bad = []
        for rel, text in files.items():
            path = out / rel
            have = path.read_bytes().replace(b"\r\n", b"\n") if path.exists() else None
            if have is None or hashlib.sha256(have).digest() != hashlib.sha256(text.encode("utf-8")).digest():
                bad.append(rel)
        print(f"{len(files) - len(bad)}/{len(files)} files match {out}")
        for rel in bad:
            print(f"DIFFERS: {rel}")
        return 1 if bad else 0
    files = write_all(out, settings)
    print(f"wrote {len(files)} files under {out}")
    print((out / "README.md").read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
