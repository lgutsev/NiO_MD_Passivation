"""What a potential family knows about itself, independently of any engine.

A potential adapter answers three questions and no others:

* what can this model produce (:meth:`capabilities`),
* can it run here (:meth:`availability`),
* what does its model file actually contain (:meth:`inspect`).

It never integrates anything. Anything about timesteps, thermostats or
trajectories belongs to an engine, and anything about wiring the two together
belongs to a bridge.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar

from ..capabilities import CapabilitySet
from ..environment import Availability
from ..specs import PotentialSpec


class PotentialAdapter(ABC):
    """Base class for the model-side half of a bridge."""

    kind: ClassVar[str] = "abstract"
    #: Modules needed to evaluate this model at all, engine-independent.
    requires: ClassVar[tuple[str, ...]] = ()

    def __init__(self, spec: PotentialSpec) -> None:
        self.spec = spec

    @abstractmethod
    def capabilities(self) -> CapabilitySet:
        """What this model can produce, before an engine constrains it."""

    @abstractmethod
    def availability(self) -> Availability:
        """Whether this model can be evaluated on this machine."""

    def inspect(self) -> dict:
        """Report what can be discovered about the model file itself.

        Implementations must degrade gracefully: on a machine without the
        backend installed, this returns what the spec declares and says that
        nothing was read from the file, rather than raising.
        """
        return {
            "kind": self.kind,
            "label": self.spec.label,
            "declared_elements": list(self.spec.elements),
            "energy_convention": self.spec.energy_convention,
            "discovered": None,
        }


__all__ = ["PotentialAdapter"]
