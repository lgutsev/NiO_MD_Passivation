"""The MLIP subsystem must not change what a classical command imports.

The top-level parser registers the ``mlip`` group on every invocation, so the
MLIP configuration layer is loaded even for ``nio-md-prep build --help``.
That is acceptable only while it stays free of every heavy backend. These
checks run in a *fresh interpreter* (``sys.executable -c ...``): an
in-process ``sys.modules`` check is order-dependent and meaningless once any
earlier test has imported torch.

The subprocess imports this worktree's ``src`` explicitly, because the
interpreter's site-packages may hold an editable install of another checkout.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"

#: Never imported by building the CLI parser or by a classical command's parser.
FORBIDDEN = ("torch", "mace", "openmm", "openmmml", "lammps", "ase")

PROBE = textwrap.dedent(
    """
    import json, sys
    sys.path.insert(0, {src!r})
    import nio_md_prep
    assert nio_md_prep.__file__.startswith({src!r}), nio_md_prep.__file__
    from nio_md_prep.cli import main
    argv = json.loads(sys.argv[1])
    try:
        code = main(argv)
    except SystemExit as exc:
        code = exc.code
    forbidden = set(json.loads(sys.argv[2]))
    leaked = sorted(m for m in sys.modules if m.split(".")[0] in forbidden)
    print(json.dumps({{"code": code, "leaked": leaked}}))
    """
)


def run_cli_in_fresh_interpreter(argv, forbidden=FORBIDDEN, cwd=None) -> dict:
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            PROBE.format(src=str(SRC)),
            json.dumps(list(argv)),
            json.dumps(list(forbidden)),
        ],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
        env=env,
        cwd=cwd,
    )
    assert completed.returncode == 0, completed.stderr
    return json.loads(completed.stdout.strip().splitlines()[-1])


@pytest.mark.parametrize(
    "argv",
    [
        ["--help"],
        ["validate", "--help"],
        ["build", "--help"],
        ["analyze-coverage", "--help"],
        ["mlip", "--help"],
        ["mlip", "validate", "--help"],
    ],
    ids=lambda argv: " ".join(argv),
)
def test_building_the_cli_parser_imports_no_backend(argv):
    outcome = run_cli_in_fresh_interpreter(argv)
    assert outcome["code"] == 0
    assert outcome["leaked"] == [], (
        f"'nio-md-prep {' '.join(argv)}' imported {outcome['leaked']}"
    )


def test_mlip_validate_without_a_structure_imports_no_heavy_backend(tmp_path):
    """``mlip validate`` must stay executable on a laptop with nothing installed.

    ASE is allowed here (the mock route is ASE-backed); torch, MACE, OpenMM and
    LAMMPS are not.
    """
    config = tmp_path / "mock.toml"
    config.write_text(
        textwrap.dedent(
            """
            [potential]
            kind = "mock"
            elements = ["Ni", "O"]

            [engine]
            kind = "ase"

            [simulation]
            task = "singlepoint"
            """
        ),
        encoding="utf-8",
    )
    outcome = run_cli_in_fresh_interpreter(
        ["mlip", "validate", str(config)],
        forbidden=("torch", "mace", "openmm", "openmmml", "lammps"),
    )
    assert outcome["code"] == 0
    assert outcome["leaked"] == []
