"""Dynamics engines: how a potential is integrated, independently of which one.

ASE, LAMMPS and OpenMM each get a module that knows about single points,
optimisation and short diagnostic trajectories, and nothing about MACE,
DeepMD or any other model family. Joining the two halves is the registry's
job; see :mod:`nio_md_prep.mlip.registry`.
"""
from __future__ import annotations

from .base import ENGINE_MODULES, EngineRuntime

__all__ = ["ENGINE_MODULES", "EngineRuntime"]
