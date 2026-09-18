"""Probing what is installed, without importing it.

Every availability question in this subsystem goes through
:func:`module_available`, which uses ``importlib.util.find_spec``. That is the
mechanism that keeps the promise in this package's docstring: ``mlip inspect``
can report that a MACE/OpenMM route is unavailable on this machine without
importing torch, and installing plain ``nio-md-prep`` never drags in PyTorch,
CUDA, OpenMM or LAMMPS.
"""
from __future__ import annotations

import shutil
from collections.abc import Iterable
from dataclasses import dataclass
from importlib.util import find_spec


def module_available(name: str) -> bool:
    """True when ``name`` could be imported, without importing it.

    ``find_spec`` on a dotted name imports the parent packages, so only
    top-level names are probed here; that is all any caller needs.
    """
    try:
        return find_spec(name.split(".", 1)[0]) is not None
    except (ImportError, ValueError):  # pragma: no cover - broken namespace packages
        return False


def missing_modules(names: Iterable[str]) -> tuple[str, ...]:
    return tuple(name for name in names if not module_available(name))


def executable_available(name: str | None) -> bool:
    return bool(name) and shutil.which(name) is not None


@dataclass(frozen=True)
class Availability:
    """Whether a route can run here, and what is missing if it cannot."""

    available: bool
    missing: tuple[str, ...] = ()
    detail: str = ""

    def __bool__(self) -> bool:
        return self.available

    def as_dict(self) -> dict:
        return {
            "available": self.available,
            "missing": list(self.missing),
            "detail": self.detail,
        }


def probe(names: Iterable[str], *, detail: str = "") -> Availability:
    missing = missing_modules(names)
    return Availability(available=not missing, missing=missing, detail=detail)


__all__ = [
    "module_available",
    "missing_modules",
    "executable_available",
    "Availability",
    "probe",
]
