"""Shared fixtures for the MLIP test suite.

Ordinary CI runs the unmarked tests: configuration parsing, bridge
resolution, unit conversions, provenance, capability failures and the mock
potential. Nothing there needs torch, LAMMPS or OpenMM.

The expensive tests are marked ``mace``, ``lammps``, ``openmm`` or ``gpu`` and
skip themselves unless the backend is importable. A trained model is supplied
through ``NIO_MD_TEST_MACE_MODEL`` -- deliberately an environment variable, so
that no large model file is ever committed to this repository.

No fixture here imports torch (or any other backend) into the test process:
availability is probed with ``find_spec``, and the one question that needs a
real import -- is CUDA usable? -- is answered in a subprocess. An in-process
import would leak into ``sys.modules`` and make every later "imports no
backend" check meaningless.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

#: Top-level modules the MLIP layer must never import on its cheap paths.
BACKEND_MODULES = ("torch", "mace", "openmm", "openmmml", "lammps")

#: Environment variable naming a trained MACE model for the optional tests.
MACE_MODEL_ENV = "NIO_MD_TEST_MACE_MODEL"
#: Optional: a LAMMPS-ready export of the same model, for the LAMMPS routes.
MACE_LAMMPS_MODEL_ENV = "NIO_MD_TEST_MACE_LAMMPS_MODEL"
#: Optional: comma-separated elements the test model covers.
MACE_ELEMENTS_ENV = "NIO_MD_TEST_MACE_ELEMENTS"


def _module_available(name: str) -> bool:
    from importlib.util import find_spec

    try:
        return find_spec(name) is not None
    except (ImportError, ValueError):
        return False


@pytest.fixture(scope="session")
def mace_model_path() -> Path:
    """The trained MACE model named by the environment, or skip."""
    raw = os.environ.get(MACE_MODEL_ENV)
    if not raw:
        pytest.skip(
            f"set {MACE_MODEL_ENV} to a trained MACE model to run this test; no model "
            "is committed to this repository"
        )
    path = Path(raw)
    if not path.exists():
        pytest.skip(f"{MACE_MODEL_ENV} points at {path}, which does not exist")
    return path


@pytest.fixture(scope="session")
def mace_elements() -> tuple[str, ...]:
    raw = os.environ.get(MACE_ELEMENTS_ENV, "Ni,O,P,C,H")
    return tuple(part.strip() for part in raw.split(",") if part.strip())


@pytest.fixture(scope="session")
def mace_lammps_model_path() -> Path | None:
    raw = os.environ.get(MACE_LAMMPS_MODEL_ENV)
    return Path(raw) if raw else None


@pytest.fixture(scope="session")
def require_mace():
    if not (_module_available("mace") and _module_available("torch")):
        pytest.skip("mace-torch and torch are not installed")


@pytest.fixture(scope="session")
def require_openmm():
    if not (_module_available("openmm") and _module_available("openmmml")):
        pytest.skip("openmm and openmmml (OpenMM-ML) are not installed")


@pytest.fixture(scope="session")
def require_lammps():
    if not _module_available("lammps"):
        pytest.skip("the LAMMPS python module is not installed")


@pytest.fixture(scope="session")
def require_gpu():
    """Skip unless torch sees a CUDA device -- asked in a subprocess, not here."""
    if not _module_available("torch"):
        pytest.skip("torch is not installed")
    probe = subprocess.run(
        [sys.executable, "-c", "import sys, torch; sys.exit(0 if torch.cuda.is_available() else 3)"],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if probe.returncode == 3:
        pytest.skip("no CUDA device is available")
    if probe.returncode != 0:
        pytest.skip(f"torch could not be imported to probe CUDA: {probe.stderr.strip()[-300:]}")


@pytest.fixture
def backend_import_guard():
    """Check that the code under test imports no MLIP backend.

    Call the returned function after exercising the code. If a backend was
    already imported by an earlier test in this process, an in-process check
    cannot tell whether *this* code imported it, so the test is skipped with
    that reason; ``tests/test_mlip_isolation.py`` covers the same property
    in a fresh interpreter.
    """
    preloaded = sorted(name for name in BACKEND_MODULES if name in sys.modules)
    if preloaded:
        pytest.skip(
            f"{', '.join(preloaded)} already imported in this test process by an earlier "
            "test; the subprocess isolation tests cover this check"
        )

    def check() -> None:
        leaked = sorted(name for name in BACKEND_MODULES if name in sys.modules)
        assert not leaked, f"the code under test imported {', '.join(leaked)}"

    return check


@pytest.fixture
def nio_structure():
    """A small periodic NiO cell, used as the interchange structure."""
    ase_build = pytest.importorskip("ase.build")
    return ase_build.bulk("NiO", "rocksalt", a=4.17).repeat((2, 2, 2))


@pytest.fixture
def rattled_nio_structure(nio_structure):
    """The same cell, displaced so that the forces are not all zero by symmetry."""
    atoms = nio_structure.copy()
    atoms.rattle(stdev=0.05, seed=20250917)
    return atoms
