"""Potential families: what a model is, independently of how it is integrated.

Each module here describes one family -- MACE, LAMMPS-native MLIPs, and an
analytic test potential -- and none of them knows anything about timesteps or
thermostats. Adapters are looked up lazily by ``potential.kind`` so that
importing this package costs nothing.
"""
from __future__ import annotations

import importlib

from ..errors import UnknownPotentialError

#: ``potential.kind`` -> module providing ``build_adapter(spec)``.
ADAPTER_MODULES = {
    "mace": ".mace",
    "lammps": ".lammps_mlip",
    "mock": ".mock",
}


def build_adapter(spec):
    """Return the :class:`PotentialAdapter` for a potential spec."""
    try:
        module_name = ADAPTER_MODULES[spec.kind]
    except KeyError:
        raise UnknownPotentialError(
            f"no potential adapter for kind {spec.kind!r}; known kinds: "
            f"{', '.join(sorted(ADAPTER_MODULES))}"
        ) from None
    module = importlib.import_module(module_name, package=__name__)
    return module.build_adapter(spec)


__all__ = ["ADAPTER_MODULES", "build_adapter"]
