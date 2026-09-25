"""Capability advertisement and negotiation.

Nothing here assumes a potential can do everything. A bridge advertises a
:class:`CapabilitySet`; a job declares a :class:`RequirementSet`; the resolver
refuses the job *before* a batch script is generated if the two do not meet.

The capability set that matters is always the **intersection** of what the
potential can produce and what the engine can carry. A MACE model can produce
a stress tensor, but a given LAMMPS route may not expose atomic virials; an
OpenMM route has no stress reporting at all, so a constant-pressure job is
rejected there rather than being run with a silently absent virial.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, replace

from .errors import CapabilityError, ElementCoverageError
from .units import TOTAL

#: Boolean capability names, in report order. ``elements``, ``precisions``,
#: ``engines`` and the convention fields are handled separately because they
#: are sets rather than flags.
FLAGS = (
    "energy",
    "forces",
    "stress",
    "per_atom_energy",
    "periodic",
    "partial_periodic",
    "fixed_atoms",
    "gpu",
)

#: Flags that describe what the *engine route* does with a structure rather
#: than what a model can compute: any full-system potential can be evaluated
#: on a slab or with frozen atoms, but whether the engine honours per-axis
#: periodicity or ASE ``FixAtoms`` is the engine's business. In
#: :meth:`CapabilitySet.intersect` they are taken from the engine side, so a
#: potential adapter never has to declare them. Both default to ``False``:
#: an engine that has not been taught to honour them refuses such jobs.
ENGINE_FLAGS = ("partial_periodic", "fixed_atoms")

FLAG_DESCRIPTIONS = {
    "energy": "total potential energy",
    "forces": "per-atom forces",
    "stress": "stress tensor / virial",
    "per_atom_energy": "per-atom (site) energy decomposition",
    "periodic": "periodic cells with the minimum-image convention",
    "partial_periodic": (
        "cells periodic along some axes only (slabs, wires), each axis honoured "
        "independently"
    ),
    "fixed_atoms": (
        "ASE FixAtoms honoured in velocity initialisation, integration and the "
        "temperature degrees of freedom"
    ),
    "gpu": "GPU-accelerated evaluation",
}


@dataclass(frozen=True)
class CapabilitySet:
    """What a potential, an engine, or a potential/engine bridge can deliver.

    ``elements`` and ``precisions`` use ``None`` to mean "unrestricted", which
    is the honest answer for an engine: LAMMPS does not care which elements a
    model covers. Intersecting an unrestricted set with a restricted one
    yields the restricted one. An *empty* ``precisions`` set means the route
    can guarantee no precision at all (a LAMMPS pair style's precision is
    fixed by its build), so any precision requirement is refused.
    """

    energy: bool = True
    forces: bool = False
    stress: bool = False
    per_atom_energy: bool = False
    periodic: bool = False
    partial_periodic: bool = False
    fixed_atoms: bool = False
    gpu: bool = False
    elements: frozenset[str] | None = None
    precisions: frozenset[str] | None = None
    engines: frozenset[str] | None = None
    native_energy_convention: str = TOTAL
    convertible_energy_conventions: frozenset[str] = frozenset()
    native_units: str = "canonical"
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("elements", "precisions", "engines", "convertible_energy_conventions"):
            value = getattr(self, name)
            if value is not None and not isinstance(value, frozenset):
                object.__setattr__(self, name, frozenset(value))
        # Notes accumulate as potential and engine sets are intersected; the same
        # warning arriving from both sides should be said once.
        object.__setattr__(self, "notes", tuple(dict.fromkeys(self.notes)))

    @property
    def energy_conventions(self) -> frozenset[str]:
        """Conventions this bridge can report, native plus convertible ones."""
        return frozenset({self.native_energy_convention}) | self.convertible_energy_conventions

    def intersect(self, other: "CapabilitySet") -> "CapabilitySet":
        """Combine a potential's capabilities with an engine's.

        Flags are ANDed: a capability survives only if both sides have it.
        The :data:`ENGINE_FLAGS` are the exception: they come from ``other``
        (the engine side) alone. Set-valued fields are intersected, with
        ``None`` acting as the identity. The energy convention and native
        units come from ``other`` too, because that is what actually produces
        numbers.
        """
        return CapabilitySet(
            **{
                name: (
                    getattr(other, name)
                    if name in ENGINE_FLAGS
                    else getattr(self, name) and getattr(other, name)
                )
                for name in FLAGS
            },
            elements=_intersect(self.elements, other.elements),
            precisions=_intersect(self.precisions, other.precisions),
            engines=_intersect(self.engines, other.engines),
            native_energy_convention=other.native_energy_convention,
            convertible_energy_conventions=other.convertible_energy_conventions,
            native_units=other.native_units,
            notes=self.notes + other.notes,
        )

    def with_notes(self, *notes: str) -> "CapabilitySet":
        return replace(self, notes=self.notes + tuple(notes))

    def supports_elements(self, symbols: Iterable[str]) -> frozenset[str]:
        """Return the symbols this capability set does *not* cover."""
        if self.elements is None:
            return frozenset()
        return frozenset(symbols) - self.elements

    def as_dict(self) -> dict:
        return {
            **{name: getattr(self, name) for name in FLAGS},
            "elements": None if self.elements is None else sorted(self.elements),
            "precisions": None if self.precisions is None else sorted(self.precisions),
            "engines": None if self.engines is None else sorted(self.engines),
            "native_energy_convention": self.native_energy_convention,
            "energy_conventions": sorted(self.energy_conventions),
            "native_units": self.native_units,
            "notes": list(self.notes),
        }


def _intersect(a: frozenset | None, b: frozenset | None) -> frozenset | None:
    if a is None:
        return b
    if b is None:
        return a
    return a & b


@dataclass(frozen=True)
class RequirementSet:
    """What a job needs, with a reason attached to every requirement.

    The reasons are not decoration: they are what turns "unsupported" into an
    actionable message such as ``stress/virial: the npt ensemble integrates
    the cell against the virial``.
    """

    energy: bool = True
    forces: bool = False
    stress: bool = False
    per_atom_energy: bool = False
    periodic: bool = False
    partial_periodic: bool = False
    fixed_atoms: bool = False
    gpu: bool = False
    elements: frozenset[str] = frozenset()
    precision: str | None = None
    engine: str | None = None
    energy_convention: str | None = None
    #: Per-axis periodicity of the structure the job runs on, when known.
    #: Carried so an engine can act on (and a manifest can record) exactly
    #: which cell axes are periodic, not just whether any are.
    periodic_axes: tuple[bool, bool, bool] | None = None
    reasons: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.elements, frozenset):
            object.__setattr__(self, "elements", frozenset(self.elements))
        if not isinstance(self.reasons, dict):
            object.__setattr__(self, "reasons", dict(self.reasons))
        if self.periodic_axes is not None:
            axes = tuple(bool(v) for v in self.periodic_axes)
            if len(axes) != 3:
                raise ValueError(f"periodic_axes needs three entries; got {self.periodic_axes!r}")
            object.__setattr__(self, "periodic_axes", axes)

    def merge(self, other: "RequirementSet") -> "RequirementSet":
        reasons = dict(self.reasons)
        reasons.update(other.reasons)
        return RequirementSet(
            **{name: getattr(self, name) or getattr(other, name) for name in FLAGS},
            elements=self.elements | other.elements,
            precision=other.precision or self.precision,
            engine=other.engine or self.engine,
            energy_convention=other.energy_convention or self.energy_convention,
            periodic_axes=(
                other.periodic_axes if other.periodic_axes is not None else self.periodic_axes
            ),
            reasons=reasons,
        )

    def requiring(self, capability: str, reason: str) -> "RequirementSet":
        reasons = dict(self.reasons)
        reasons[capability] = reason
        return replace(self, **{capability: True}, reasons=reasons)

    def unmet(self, capabilities: CapabilitySet) -> list[str]:
        """List every requirement ``capabilities`` fails to satisfy."""
        problems: list[str] = []
        for name in FLAGS:
            if getattr(self, name) and not getattr(capabilities, name):
                reason = self.reasons.get(name, "required by this job")
                problems.append(
                    f"{FLAG_DESCRIPTIONS[name]} ({name}) is not available: {reason}"
                )
        if self.precision is not None and capabilities.precisions is not None:
            if self.precision not in capabilities.precisions:
                reason = self.reasons.get("precision")
                offered = (
                    f"offered: {', '.join(sorted(capabilities.precisions))}"
                    if capabilities.precisions
                    else "this route guarantees no precision"
                )
                problems.append(
                    f"precision {self.precision!r} is not available ({offered})"
                    + (f": {reason}" if reason else "")
                )
        if self.engine is not None and capabilities.engines is not None:
            if self.engine not in capabilities.engines:
                problems.append(
                    f"engine {self.engine!r} is not supported by this potential "
                    f"(supported: {', '.join(sorted(capabilities.engines))})"
                )
        if self.energy_convention is not None:
            if self.energy_convention not in capabilities.energy_conventions:
                problems.append(
                    f"energy convention {self.energy_convention!r} cannot be reported "
                    f"(available: {', '.join(sorted(capabilities.energy_conventions))})"
                )
        return problems

    def as_dict(self) -> dict:
        return {
            **{name: getattr(self, name) for name in FLAGS},
            "elements": sorted(self.elements),
            "precision": self.precision,
            "engine": self.engine,
            "energy_convention": self.energy_convention,
            "periodic_axes": None if self.periodic_axes is None else list(self.periodic_axes),
            "reasons": dict(self.reasons),
        }


def check_elements(
    capabilities: CapabilitySet,
    symbols: Iterable[str],
    *,
    label: str,
) -> None:
    """Fail closed when a structure carries an element the model never saw.

    Called before any engine is constructed. A Ni/O/P/C/H model handed an
    extra atom type must stop here, not produce a plausible-looking number.
    """
    symbols = list(symbols)
    missing = capabilities.supports_elements(symbols)
    if missing:
        raise ElementCoverageError(missing, capabilities.elements or (), label)


def negotiate(
    capabilities: CapabilitySet,
    requirements: RequirementSet,
    *,
    label: str,
) -> None:
    """Raise :class:`CapabilityError` unless ``capabilities`` covers everything."""
    problems = requirements.unmet(capabilities)
    if problems:
        raise CapabilityError(f"{label} cannot satisfy this job:", problems)


__all__ = [
    "FLAGS",
    "ENGINE_FLAGS",
    "FLAG_DESCRIPTIONS",
    "CapabilitySet",
    "RequirementSet",
    "check_elements",
    "negotiate",
    "CapabilityError",
]
