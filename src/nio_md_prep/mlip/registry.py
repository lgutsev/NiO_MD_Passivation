"""The bridge registry: the one place that knows which pairings exist.

Potential type and dynamics engine are independent concepts. Their *joining*
is a third concept -- a bridge -- and it lives here, in a table, rather than
being scattered as ``if potential == "mace"`` tests through preparation code.

The table below is the compatibility matrix, in executable form:

===================  ======================  ============================  ============================
Potential            ASE                     LAMMPS                        OpenMM
===================  ======================  ============================  ============================
MACE                 direct ASE calculator   via MACE ML-IAP (or pair      via OpenMM-ML
                                             style)
LAMMPS-native MLIP   via the ASE -> LAMMPS   native                        unsupported
                     bridge
===================  ======================  ============================  ============================

The single unsupported cell is registered explicitly, with a reason, so that
asking for it produces an :class:`~nio_md_prep.mlip.errors.UnsupportedCombinationError`
rather than an import traceback from a half-written adapter. There is
deliberately no generic "LAMMPS potential -> OpenMM" adapter: a LAMMPS pair
style is compiled C++ inside LAMMPS, and pretending otherwise would fabricate
physics. A model that has a separate OpenMM (or ASE) implementation should be
registered as its own potential kind, the way MACE is.

Bridge classes are loaded lazily through an import string, so importing this
module -- and therefore ``nio-md-prep`` as a whole -- never pulls in torch,
CUDA, OpenMM or LAMMPS.
"""
from __future__ import annotations

import importlib
from collections.abc import Iterable
from dataclasses import dataclass

from .errors import (
    UnknownEngineError,
    UnknownPotentialError,
    UnsupportedCombinationError,
)
from .specs import ENGINE_KINDS, EngineSpec, PotentialSpec

SUPPORTED = "supported"
UNSUPPORTED = "unsupported"


@dataclass(frozen=True)
class BridgeRegistration:
    """One cell of the compatibility matrix.

    ``module``/``attribute`` name the bridge class without importing it.
    ``requires`` lists the importable modules the route needs, so
    ``mlip inspect`` and ``mlip validate`` can report availability without
    attempting the import chain.
    """

    potential: str
    engine: str
    implementation: str
    summary: str
    status: str = SUPPORTED
    module: str | None = None
    attribute: str | None = None
    requires: tuple[str, ...] = ()
    packages: tuple[str, ...] = ()
    priority: int = 100
    reason: str = ""
    notes: tuple[str, ...] = ()

    @property
    def key(self) -> tuple[str, str]:
        return (self.potential, self.engine)

    @property
    def supported(self) -> bool:
        return self.status == SUPPORTED

    def load(self):
        """Import and return the bridge class. The first heavy import happens here."""
        if not self.supported:
            raise UnsupportedCombinationError(self.potential, self.engine, self.reason)
        module = importlib.import_module(self.module, package=__package__)
        return getattr(module, self.attribute)

    def instantiate(self, potential: PotentialSpec, engine: EngineSpec):
        return self.load()(potential, engine, registration=self)

    def as_dict(self) -> dict:
        return {
            "potential": self.potential,
            "engine": self.engine,
            "implementation": self.implementation,
            "status": self.status,
            "summary": self.summary,
            "requires": list(self.requires),
            "packages": list(self.packages),
            "reason": self.reason,
            "notes": list(self.notes),
        }


class BridgeRegistry:
    """An ordered table of registrations, keyed by (potential kind, engine kind)."""

    def __init__(self) -> None:
        self._entries: list[BridgeRegistration] = []

    def register(self, registration: BridgeRegistration) -> BridgeRegistration:
        existing = [
            e
            for e in self._entries
            if e.key == registration.key and e.implementation == registration.implementation
        ]
        if existing:
            raise ValueError(
                f"bridge {registration.implementation!r} is already registered for "
                f"{registration.potential} -> {registration.engine}"
            )
        self._entries.append(registration)
        self._entries.sort(key=lambda e: (e.potential, e.engine, e.priority, e.implementation))
        return registration

    def entries(self) -> tuple[BridgeRegistration, ...]:
        return tuple(self._entries)

    def potentials(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(e.potential for e in self._entries))

    def engines(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(e.engine for e in self._entries))

    def for_pair(self, potential: str, engine: str) -> tuple[BridgeRegistration, ...]:
        return tuple(e for e in self._entries if e.key == (potential, engine))

    def resolve(
        self,
        potential: str,
        engine: str,
        *,
        implementation: str | None = None,
        preferences: Iterable[str] = (),
    ) -> BridgeRegistration:
        """Select the implementation for a (potential, engine) pair.

        Raises a specific, actionable error for each way this can fail: an
        unknown potential kind, an unknown engine kind, a combination that is
        registered as impossible, and a requested implementation that does not
        exist for an otherwise-valid pair.
        """
        if potential not in self.potentials():
            raise UnknownPotentialError(
                f"unknown potential kind {potential!r}; registered kinds: "
                f"{', '.join(self.potentials())}"
            )
        if engine not in ENGINE_KINDS:
            raise UnknownEngineError(
                f"unknown engine kind {engine!r}; known engines: {', '.join(ENGINE_KINDS)}"
            )
        candidates = self.for_pair(potential, engine)
        if not candidates:
            raise UnsupportedCombinationError(
                potential,
                engine,
                "no bridge is registered for this combination",
                alternatives=self._alternatives(potential),
            )
        blocked = [c for c in candidates if not c.supported]
        usable = [c for c in candidates if c.supported]
        if not usable:
            entry = blocked[0]
            raise UnsupportedCombinationError(
                potential,
                engine,
                entry.reason or entry.summary,
                alternatives=self._alternatives(potential),
            )
        if implementation is not None:
            for entry in usable:
                if entry.implementation == implementation:
                    return entry
            raise UnsupportedCombinationError(
                potential,
                engine,
                f"no implementation named {implementation!r} "
                f"(available: {', '.join(e.implementation for e in usable)})",
            )
        for preferred in preferences:
            for entry in usable:
                if entry.implementation == preferred:
                    return entry
        return usable[0]

    def _alternatives(self, potential: str) -> list[str]:
        return [
            f"{e.potential} -> {e.engine} ({e.implementation}: {e.summary})"
            for e in self._entries
            if e.potential == potential and e.supported
        ]

    def matrix(self) -> dict[str, dict[str, list[dict]]]:
        """The compatibility matrix as plain data, for reporting and tests."""
        table: dict[str, dict[str, list[dict]]] = {}
        for potential in self.potentials():
            row: dict[str, list[dict]] = {}
            for engine in ENGINE_KINDS:
                row[engine] = [e.as_dict() for e in self.for_pair(potential, engine)]
            table[potential] = row
        return table


REGISTRY = BridgeRegistry()


def register(**kwargs) -> BridgeRegistration:
    return REGISTRY.register(BridgeRegistration(**kwargs))


# ---------------------------------------------------------------------------
# The compatibility matrix
# ---------------------------------------------------------------------------

register(
    potential="mace",
    engine="ase",
    implementation="mace-ase-calculator",
    summary="MACE's own ASE calculator, driven by ASE's optimisers and integrators",
    module=".bridges.mace_ase",
    attribute="MaceAseBridge",
    requires=("mace", "torch", "ase"),
    packages=("mace-torch", "torch", "ase"),
    priority=10,
    notes=(
        "reports the total energy including atomic self-energies (E0)",
        "the reference route for cross-engine equivalence checks",
    ),
)

register(
    potential="mace",
    engine="lammps",
    implementation="mliap",
    summary="MACE through the LAMMPS ML-IAP interface (pair_style mliap unified)",
    module=".bridges.mace_lammps",
    attribute="MaceLammpsMliapBridge",
    requires=("lammps", "torch"),
    packages=("LAMMPS with the ML-IAP package", "mace-torch", "torch"),
    priority=10,
    notes=(
        "the newer MACE/LAMMPS route: GPU acceleration, multi-GPU inference "
        "and atomic virials",
        "needs a model exported to the ML-IAP format before use",
        "MACE advises care when benchmarking LAMMPS output against the ASE "
        "calculator; cross-engine equivalence is an explicit acceptance test here",
    ),
)

register(
    potential="mace",
    engine="lammps",
    implementation="pair-mace",
    summary="MACE through a compiled-in MACE pair style (pair_style mace)",
    module=".bridges.mace_lammps",
    attribute="MaceLammpsPairStyleBridge",
    requires=("lammps",),
    packages=("LAMMPS built with the MACE pair style",),
    priority=20,
    notes=(
        "the older route; kept because not every cluster build has ML-IAP",
        "per-atom virials are not assumed to be available",
    ),
)

register(
    potential="mace",
    engine="openmm",
    implementation="openmm-ml",
    summary="MACE through OpenMM-ML's MACEPotential",
    module=".bridges.mace_openmm",
    attribute="MaceOpenMMBridge",
    requires=("openmmml", "openmm", "mace", "torch"),
    packages=("openmm-ml", "openmm", "mace-torch", "torch"),
    priority=10,
    notes=(
        "accepts both pretrained (mace-off / mace-mp) and locally trained models",
        "defaults to the INTERACTION energy (atomic self-energies removed); this "
        "subsystem harmonises the convention before any comparison",
        "OpenMM reports no stress tensor here, so constant-pressure jobs are "
        "rejected by capability negotiation rather than run without a virial",
    ),
)

register(
    potential="lammps",
    engine="lammps",
    implementation="native",
    summary="a LAMMPS-native MLIP pair style running in LAMMPS, as intended",
    module=".bridges.lammps_native",
    attribute="LammpsNativeBridge",
    requires=("lammps",),
    packages=("LAMMPS with the pair style's package",),
    priority=10,
    notes=("the rendered pair_style and pair_coeff lines are preserved verbatim",),
)

register(
    potential="lammps",
    engine="ase",
    implementation="ase-lammps",
    summary="a LAMMPS-native MLIP driven from ASE through ASE's LAMMPS calculator",
    module=".bridges.lammps_ase",
    attribute="LammpsAseBridge",
    requires=("ase", "lammps"),
    packages=("ase", "LAMMPS (python module or executable)"),
    priority=10,
    notes=(
        "LAMMPS still evaluates the potential; ASE only drives it, which is what "
        "makes ASE optimisers and analysis available",
        "the type_map is what maps an ASE Atoms object onto LAMMPS types",
    ),
)

register(
    potential="lammps",
    engine="openmm",
    implementation="none",
    status=UNSUPPORTED,
    summary="not possible",
    reason=(
        "a LAMMPS-native MLIP is compiled C++ inside LAMMPS; OpenMM cannot call it, "
        "and no generic LAMMPS-pair-style-to-OpenMM adapter exists or is being "
        "faked here. If the underlying model has a separate OpenMM (or ASE) "
        "implementation, register that implementation as its own potential kind -- "
        "the way MACE is registered -- rather than routing it through LAMMPS"
    ),
    notes=("registered explicitly so the failure is a clear error, not a traceback",),
)

register(
    potential="mock",
    engine="ase",
    implementation="mock-ase",
    summary="an analytic test potential in ASE; exercises this machinery without torch",
    module=".bridges.mock_ase",
    attribute="MockAseBridge",
    requires=("ase",),
    packages=("ase",),
    priority=10,
    notes=("diagnostic only: not a scientific potential",),
)


def resolve_bridge(
    potential: str,
    engine: str,
    *,
    implementation: str | None = None,
    preferences: Iterable[str] = (),
) -> BridgeRegistration:
    """Select the implementation joining ``potential`` to ``engine``.

    Returns the registration, not a live bridge: resolution is cheap and
    import-free, which is what makes ``mlip validate`` usable on a machine
    that cannot run the job.
    """
    return REGISTRY.resolve(
        potential, engine, implementation=implementation, preferences=preferences
    )


def build_bridge(
    potential: PotentialSpec,
    engine: EngineSpec,
    *,
    implementation: str | None = None,
):
    """Resolve and instantiate the bridge for two specs."""
    preferences = tuple(getattr(potential, "implementation", ()) or ())
    registration = resolve_bridge(
        potential.kind, engine.kind, implementation=implementation, preferences=preferences
    )
    return registration.instantiate(potential, engine)


def compatibility_matrix() -> dict:
    return REGISTRY.matrix()


def registrations() -> tuple[BridgeRegistration, ...]:
    return REGISTRY.entries()


__all__ = [
    "SUPPORTED",
    "UNSUPPORTED",
    "BridgeRegistration",
    "BridgeRegistry",
    "REGISTRY",
    "register",
    "resolve_bridge",
    "build_bridge",
    "compatibility_matrix",
    "registrations",
]
