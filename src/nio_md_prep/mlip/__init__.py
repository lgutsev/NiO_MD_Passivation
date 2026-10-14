"""Machine-learned interatomic potentials: specifications, routes, provenance.

A preliminary subsystem, isolated from the classical preparation workflows in
:mod:`nio_md_prep.build`, :mod:`nio_md_prep.agglomeration`,
:mod:`nio_md_prep.lammps` and :mod:`nio_md_prep.analysis`. Nothing here changes
how an existing classical study is built; the two live side by side.

The architectural rule
----------------------
**Potential type and dynamics engine are independent concepts.** A MACE model
does not know whether it will be integrated by ASE, LAMMPS or OpenMM; an
engine does not know which model family it is driving. Their joining is a
third thing -- a *bridge* -- and every bridge lives in one table, in
:mod:`nio_md_prep.mlip.registry`, rather than as scattered ``if mace:`` tests
through preparation code.

The compatibility matrix
------------------------
==================  =====================  ==========================  ====================
Potential           ASE                    LAMMPS                      OpenMM
==================  =====================  ==========================  ====================
MACE                supported directly     supported via MACE/ML-IAP   via OpenMM-ML
LAMMPS-native MLIP  via ASE -> LAMMPS      supported natively          **unsupported**
==================  =====================  ==========================  ====================

The unsupported cell is registered explicitly, so asking for it raises
:class:`~nio_md_prep.mlip.errors.UnsupportedCombinationError` rather than
failing somewhere inside an import. There is deliberately no generic
"LAMMPS pair style -> OpenMM" adapter: a model with a separate OpenMM
implementation should be registered as its own potential kind, the way MACE is.

Scope
-----
Full-system MLIP only. ML/MM partitioning, fixed ML regions, electrostatic
embedding and bonds crossing an ML/MM boundary are out of scope; ``region``
exists in the schema with ``"all"`` as its only implemented value so that a
later ``"selection"`` needs no schema change.

Lazy imports
------------
Importing this package -- and therefore ``nio-md-prep`` -- pulls in no torch,
CUDA, OpenMM or LAMMPS. Specs, capabilities, configuration parsing, bridge
resolution and provenance are all import-free; a backend is imported only when
a bridge is actually constructed. The names below resolve on first access.
"""
from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

__all__ = [
    # errors
    "MlipError",
    "ConfigError",
    "UnsupportedCombinationError",
    "MissingDependencyError",
    "CapabilityError",
    "ElementCoverageError",
    "EnergyConventionError",
    "ModelIntegrityError",
    # specs
    "PotentialSpec",
    "MacePotentialSpec",
    "LammpsMlipPotentialSpec",
    "MockPotentialSpec",
    "EngineSpec",
    "SimulationSpec",
    "StructureSpec",
    "JobSpec",
    # capabilities
    "CapabilitySet",
    "RequirementSet",
    # registry
    "resolve_bridge",
    "build_bridge",
    "compatibility_matrix",
    # results and provenance
    "PotentialResult",
    "TrajectoryResult",
    "ComparisonResult",
    "compare_results",
    "build_manifest",
    "write_manifest",
    "read_manifest",
    "MANIFEST_NAME",
    # configuration
    "parse_job",
    "parse_job_mapping",
    # jobs
    "inspect_environment",
    "validate_job",
    "run_singlepoint",
    "run_smoke_md",
    "compare_engines",
]

#: Exported name -> the submodule that defines it. Everything is resolved on
#: first attribute access, so ``from nio_md_prep import mlip`` stays cheap.
_EXPORTS = {
    "MlipError": "errors",
    "ConfigError": "errors",
    "UnsupportedCombinationError": "errors",
    "MissingDependencyError": "errors",
    "CapabilityError": "errors",
    "ElementCoverageError": "errors",
    "EnergyConventionError": "errors",
    "ModelIntegrityError": "errors",
    "PotentialSpec": "specs",
    "MacePotentialSpec": "specs",
    "LammpsMlipPotentialSpec": "specs",
    "MockPotentialSpec": "specs",
    "EngineSpec": "specs",
    "SimulationSpec": "specs",
    "StructureSpec": "specs",
    "JobSpec": "specs",
    "CapabilitySet": "capabilities",
    "RequirementSet": "capabilities",
    "resolve_bridge": "registry",
    "build_bridge": "registry",
    "compatibility_matrix": "registry",
    "PotentialResult": "results",
    "TrajectoryResult": "results",
    "ComparisonResult": "results",
    "compare_results": "results",
    "build_manifest": "provenance",
    "write_manifest": "provenance",
    "read_manifest": "provenance",
    "MANIFEST_NAME": "provenance",
    "parse_job": "config",
    "parse_job_mapping": "config",
    "inspect_environment": "jobs",
    "validate_job": "jobs",
    "run_singlepoint": "jobs",
    "run_smoke_md": "jobs",
    "compare_engines": "jobs",
}


def __getattr__(name: str):
    try:
        module_name = _EXPORTS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    module = importlib.import_module(f".{module_name}", __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(__all__)


if TYPE_CHECKING:  # pragma: no cover - import-time cost is the whole point
    from .capabilities import CapabilitySet, RequirementSet
    from .config import parse_job, parse_job_mapping
    from .errors import (
        CapabilityError,
        ConfigError,
        ElementCoverageError,
        EnergyConventionError,
        MissingDependencyError,
        MlipError,
        ModelIntegrityError,
        UnsupportedCombinationError,
    )
    from .jobs import (
        compare_engines,
        inspect_environment,
        run_singlepoint,
        run_smoke_md,
        validate_job,
    )
    from .provenance import MANIFEST_NAME, build_manifest, read_manifest, write_manifest
    from .registry import build_bridge, compatibility_matrix, resolve_bridge
    from .results import ComparisonResult, PotentialResult, TrajectoryResult, compare_results
    from .specs import (
        EngineSpec,
        JobSpec,
        LammpsMlipPotentialSpec,
        MacePotentialSpec,
        MockPotentialSpec,
        PotentialSpec,
        SimulationSpec,
        StructureSpec,
    )
