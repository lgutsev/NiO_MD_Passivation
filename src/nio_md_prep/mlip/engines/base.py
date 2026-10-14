"""What an engine can carry, independently of which potential it carries.

An engine's capability set describes the *transport*, not the model: whether
this engine can report a stress tensor at all, whether it can decompose a
per-atom energy, whether it handles periodic cells, and which units it speaks
natively. Intersecting that with a potential's capability set is what produces
the honest answer for a route -- see
:meth:`nio_md_prep.mlip.capabilities.CapabilitySet.intersect`.

The one asymmetry worth stating plainly: OpenMM reports no stress tensor, so
``openmm`` advertises ``stress=False``. A constant-pressure job routed through
OpenMM is therefore rejected during capability negotiation rather than run
with a silently absent virial.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..specs import EngineSpec

#: ``engine.kind`` -> module implementing the runtime. Lazily imported.
ENGINE_MODULES = {
    "ase": ".ase_engine",
    "lammps": ".lammps_engine",
    "openmm": ".openmm_engine",
}


class EngineRuntime(ABC):
    """Base class for the engine-side half of a bridge."""

    kind: ClassVar[str] = "abstract"
    requires: ClassVar[tuple[str, ...]] = ()
    native_units: ClassVar[str] = "canonical"

    def __init__(self, spec: EngineSpec) -> None:
        self.spec = spec

    @abstractmethod
    def capabilities(self) -> CapabilitySet:
        """What this engine can carry, before a potential constrains it."""

    @abstractmethod
    def availability(self) -> Availability:
        """Whether this engine can run here."""

    def describe(self) -> dict:
        return {
            "kind": self.kind,
            "native_units": self.native_units,
            "requires": list(self.requires),
            "availability": self.availability().as_dict(),
            "spec": self.spec.as_dict(),
        }


__all__ = ["ENGINE_MODULES", "EngineRuntime"]
