"""Bridges: the cells of the compatibility matrix.

Each module joins one potential family to one engine. They are imported
lazily by :mod:`nio_md_prep.mlip.registry`, never eagerly from here, so that
nothing in this package pulls in torch, OpenMM or LAMMPS until a bridge is
actually built.
"""
from __future__ import annotations

from .base import Bridge

__all__ = ["Bridge"]
