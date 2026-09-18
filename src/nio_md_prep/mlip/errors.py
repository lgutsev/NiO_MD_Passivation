"""Error taxonomy for the MLIP subsystem.

Every failure this subsystem raises deliberately is one of these types. The
point is that an impossible request (``potential = "lammps"`` with
``engine = "openmm"``), a missing optional dependency, and an unmet scientific
requirement are three different situations that a caller -- or a test -- must
be able to tell apart without reading a traceback.

``MlipError`` derives from ``RuntimeError`` and ``ConfigError`` additionally
from ``ValueError`` so that ``nio_md_prep.cli.main``'s existing top-level
``except (ValueError, FileNotFoundError, FileExistsError, RuntimeError)``
turns them into a clean ``error: ...`` exit rather than a traceback.
"""
from __future__ import annotations

from collections.abc import Iterable


class MlipError(RuntimeError):
    """Base class for every deliberate MLIP-subsystem failure."""


class ConfigError(MlipError, ValueError):
    """A configuration file or spec is malformed, incomplete or contradictory."""


class UnknownPotentialError(ConfigError):
    """``potential.kind`` names something no bridge is registered for."""


class UnknownEngineError(ConfigError):
    """``engine.kind`` names something no bridge is registered for."""


class UnsupportedCombinationError(MlipError):
    """This potential kind cannot drive this engine, and no adapter can fix that.

    This is a statement about the architecture, not about the current machine.
    Installing more packages will not make the combination work; a different
    engine, or a separate implementation of the model for that engine, is
    required. Contrast with :class:`MissingDependencyError`.
    """

    def __init__(
        self,
        potential: str,
        engine: str,
        reason: str,
        *,
        alternatives: Iterable[str] = (),
    ) -> None:
        self.potential = potential
        self.engine = engine
        self.reason = reason
        self.alternatives = tuple(alternatives)
        message = f"potential {potential!r} cannot run on engine {engine!r}: {reason}"
        if self.alternatives:
            message += "\nSupported instead: " + ", ".join(self.alternatives)
        super().__init__(message)


class MissingDependencyError(MlipError):
    """The combination is supported, but this machine lacks the packages for it.

    Unlike :class:`UnsupportedCombinationError` this is recoverable by
    installing something, so the message names what to install.
    """

    def __init__(self, what: str, modules: Iterable[str], *, hint: str = "") -> None:
        self.what = what
        self.modules = tuple(modules)
        self.hint = hint
        missing = ", ".join(self.modules) if self.modules else "an optional dependency"
        message = f"{what} requires {missing}, which is not importable here."
        if hint:
            message += f"\n{hint}"
        super().__init__(message)


class NotImplementedYetError(MlipError, NotImplementedError):
    """A route exists but does not implement this operation in this release.

    Derives from ``MlipError`` -- and therefore ``RuntimeError`` -- so the
    top-level CLI's existing exception handling turns it into a clean
    ``error: ...`` message without that handler needing to change.
    """


class CapabilityError(MlipError):
    """The resolved potential/engine pair cannot supply what the job requires."""

    def __init__(self, summary: str, unmet: Iterable[str] = ()) -> None:
        self.unmet = tuple(unmet)
        message = summary
        if self.unmet:
            message += "\n" + "\n".join(f"  - {item}" for item in self.unmet)
        super().__init__(message)


class ElementCoverageError(CapabilityError):
    """The structure contains elements the potential was never trained on.

    Validated before anything is executed: a Ni/O/P/C/H model must never
    silently receive an extra atom type and report a number anyway.
    """

    def __init__(self, missing: Iterable[str], supported: Iterable[str], label: str) -> None:
        self.missing = tuple(sorted(missing))
        self.supported = tuple(sorted(supported))
        super().__init__(
            f"structure contains element(s) {', '.join(self.missing)} that potential "
            f"{label!r} does not cover (covers: {', '.join(self.supported) or 'nothing declared'})"
        )


class EnergyConventionError(MlipError):
    """Two energies cannot be compared or converted under a common convention.

    Raised instead of returning a number that silently mixes an interaction
    energy with one that includes atomic self-energies.
    """


class UnitError(MlipError):
    """A quantity was about to cross a unit boundary without being converted."""


class ModelIntegrityError(MlipError):
    """A model file's SHA256 does not match the one the configuration declares."""

    def __init__(self, path, expected: str, actual: str) -> None:
        self.path = path
        self.expected = expected
        self.actual = actual
        super().__init__(
            f"model file {path} has sha256 {actual}, but the configuration declares "
            f"{expected}. Refusing to run: provenance would be wrong."
        )


class ProvenanceError(MlipError):
    """A job manifest cannot be assembled from the information available."""


__all__ = [
    "MlipError",
    "ConfigError",
    "UnknownPotentialError",
    "UnknownEngineError",
    "UnsupportedCombinationError",
    "MissingDependencyError",
    "NotImplementedYetError",
    "CapabilityError",
    "ElementCoverageError",
    "EnergyConventionError",
    "UnitError",
    "ModelIntegrityError",
    "ProvenanceError",
]
