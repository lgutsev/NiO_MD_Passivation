"""The OpenMM engine, via OpenMM-ML.

Several things about this engine are load-bearing and easy to get wrong:

1. **Units.** OpenMM speaks kJ/mol and nm. OpenMM-ML's MACE force multiplies
   the model's eV by a hard-coded ``96.4853`` kJ/mol per eV (and eV/Angstrom
   by ``964.853``), so this module converts the potential energy and forces
   back with *that* constant (:data:`OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV`,
   recorded as ``energy_scale_kJ_per_mol_per_eV``), not the exact CODATA
   value: the exact value would leave a -3.3e-7 relative bias on every
   energy. Kinetic energies are OpenMM's own (amu nm^2/ps^2 = kJ/mol) and
   are converted with the exact constant.
2. **Energy convention.** OpenMM-ML's MACE implementation distinguishes the
   interaction energy from the energy including atomic self-energies, and
   *interaction* is its default. This module asks for whichever convention
   the job declares and records which one it got.
3. **Execution settings reach the Context.** Platform and platform
   properties (``Precision``, ``Threads``, ``DeviceIndex``) are applied to
   the single-point Context *and* to the MD ``Simulation``, and the platform
   OpenMM actually used is recorded with every property value it reports.
   The MACE torch device and dtype are OpenMM-ML ``createSystem`` arguments
   (``device=``, ``precision=``), which only OpenMM-ML >= 1.6 (the
   ``PythonForce`` implementation) honours; older releases are refused.
4. **Box vectors.** OpenMM requires a *reduced* lower-triangular box
   (``a`` along x, ``b`` in the xy plane, ``a_x >= 2|b_x|`` ...), checked
   with exact comparisons. :func:`reduced_box` rotates the cell into ASE's
   standard form, reduces it the way OpenMM does, and returns the rotation
   so positions go in rotated and forces and positions come back in the
   source basis.
5. **Geometry and constraints.** OpenMM-ML's MACE periodicity is
   all-or-nothing, so a partially periodic structure is refused. ASE
   ``FixAtoms`` atoms get zero particle mass (OpenMM never moves a massless
   particle); every other atom gets its ASE mass.
6. **No NPT, by policy.** No stress tensor is reported through this route,
   and constant-pressure dynamics are refused. A ``MonteCarloBarostat``
   would need only energies; the refusal is a policy until an OpenMM NPT
   path has been validated against another engine, not a claim that OpenMM
   cannot do it.

Only full-system MLIP is implemented. OpenMM-ML can build mixed ML/MM systems,
but that machinery -- and the still-active questions around bonds crossing an
ML/MM boundary -- is deliberately out of scope here.
"""
from __future__ import annotations

import math
import struct
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, probe
from ..errors import CapabilityError, ConfigError, MissingDependencyError, ResultError
from ..specs import (
    DEFAULT_THERMOSTAT_DAMPING_FS,
    MAX_SEED,
    OPENMM_GPU_PLATFORMS,
    SimulationSpec,
)
from ..units import INTERACTION, KJ_PER_MOL_IN_EV, OPENMM, TOTAL
from .base import EngineRuntime

#: This subsystem's convention name -> OpenMM-ML's ``returnEnergyType`` value.
ENERGY_TYPE = {
    INTERACTION: "interaction_energy",
    TOTAL: "energy",
}

#: Kept for backward compatibility; the canonical list lives in specs.
GPU_PLATFORMS = OPENMM_GPU_PLATFORMS

#: OpenMM-ML's MACE unit factors (``_computeMACE``: ``energyScale = 96.4853``,
#: ``lengthScale = 10.0``; unchanged from 1.2 through 1.8). Results are
#: converted back with exactly these, so an eV the model produced is the eV
#: reported.
OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV = 96.4853
OPENMMML_LENGTH_SCALE_A_PER_NM = 10.0

#: OpenMM-ML releases before 1.6 evaluated MACE through a TorchForce on the
#: Context's device and ignored ``device=``; 1.6 moved to ``PythonForce`` and
#: the ``device`` keyword. Only the newer semantics are supported.
MIN_OPENMMML_VERSION = (1, 6)

#: ``potential.precision`` -> OpenMM-ML ``createSystem(precision=...)``: the
#: dtype the MACE model is evaluated in (without it, the model's own dtype).
MODEL_PRECISION = {"float32": "single", "float64": "double"}

#: Requested thermostat -> the OpenMM integrator that implements it. ``None``
#: resolves to :data:`DEFAULT_THERMOSTAT`. Everything else is refused.
SUPPORTED_THERMOSTATS = {
    "langevin": "LangevinMiddleIntegrator",
    "nose-hoover": "NoseHooverIntegrator",
}
DEFAULT_THERMOSTAT = "langevin"

#: Langevin friction when no damping time is given: 1/tau with the shared
#: default tau (specs.DEFAULT_THERMOSTAT_DAMPING_FS), so OpenMM and ASE apply
#: the same thermostat to the same configuration.
DEFAULT_FRICTION_PER_PS = 1000.0 / DEFAULT_THERMOSTAT_DAMPING_FS

#: How far (relative) a reduced-form inequality may be violated by rounding
#: before it is treated as a real violation. Only exact ties (e.g. the fcc
#: lattice, where ``b_x = a_x / 2`` exactly) land inside it.
BOX_TIE_TOLERANCE = 1e-12

#: Molar gas constant in kJ/(mol K) (exact since the 2019 SI redefinition).
MOLAR_GAS_CONSTANT_KJ_PER_MOL_K = 1.380649e-23 * 6.02214076e23 / 1000.0

NPT_POLICY = (
    "constant-pressure dynamics are refused on the OpenMM route by policy: no "
    "stress tensor is reported through OpenMM-ML here, and an OpenMM NPT path "
    "(a MonteCarloBarostat needs only energies) has not been validated against "
    "another engine. Use the ASE or LAMMPS engine for NPT."
)


class OpenMMEngine(EngineRuntime):
    kind = "openmm"
    requires = ("openmm", "openmmml")
    native_units = OPENMM.name

    def capabilities(self) -> CapabilitySet:
        """What OpenMM can carry through OpenMM-ML -- notably not a stress tensor.

        ``stress=False`` and ``per_atom_energy=False`` are the honest answers
        for this route; ``stress=False`` is also what makes an NPT request
        fail during negotiation (see :data:`NPT_POLICY` for why NPT is refused
        on this route). Frozen atoms are honoured (zero particle mass);
        partially periodic cells are not, because OpenMM-ML's MACE
        periodicity is all-or-nothing.
        """
        platform = self.spec.platform
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=False,
            per_atom_energy=False,
            periodic=True,
            partial_periodic=False,
            fixed_atoms=True,
            gpu=platform is None or platform in OPENMM_GPU_PLATFORMS,
            elements=None,
            precisions=None,
            engines=frozenset({"openmm"}),
            native_units=self.native_units,
            notes=(
                "no stress tensor is reported through OpenMM-ML; constant-pressure "
                "jobs are refused by policy until an OpenMM NPT path is validated",
                "no per-atom energy decomposition is exposed by this route",
                "OpenMM-ML's MACE periodicity is all-or-nothing: slabs/wires are refused",
                "FixAtoms atoms are frozen by zero particle mass",
            ),
        )

    def availability(self) -> Availability:
        return probe(self.requires, detail="OpenMM-ML machine-learned potentials")


# ---------------------------------------------------------------------------
# Pure helpers (no OpenMM import): box reduction, platform properties, checks
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReducedBox:
    """A periodic cell in the reduced form OpenMM accepts, and the way back.

    ``vectors_nm`` are exactly the three box vectors handed to OpenMM (nm),
    satisfying :func:`box_violations` with OpenMM's exact comparisons.
    ``rotation`` is the orthogonal ``Q`` of ASE's ``Cell.standard_form``:
    a source-frame row vector ``r`` is ``r @ Q.T`` in the engine frame, and
    an engine-frame force or position ``f`` is ``f @ Q`` in the source frame.
    ``lattice_change`` is the integer matrix ``M`` with
    ``vectors = M @ (source cell @ Q.T)`` (in Angstrom): ``|det M| = 1``, so
    both bases generate the same lattice and the physics is unchanged.
    ``basis_flipped`` records that a left-handed input basis had its third
    vector negated (same lattice) so a positive-diagonal form exists.
    ``tie_adjustment_nm`` is the largest component moved onto an exact tie
    (``b_x = a_x/2`` ...) that rounding had broken by at most
    :data:`BOX_TIE_TOLERANCE` relative.
    """

    vectors_nm: tuple[tuple[float, float, float], ...]
    rotation: tuple[tuple[float, float, float], ...]
    lattice_change: tuple[tuple[int, int, int], ...]
    basis_flipped: bool
    tie_adjustment_nm: float

    def to_engine(self, vectors):
        """Source-frame row vectors (any length unit) -> engine frame."""
        import numpy as np

        return np.asarray(vectors, dtype=float) @ np.asarray(self.rotation).T

    def to_source(self, vectors):
        """Engine-frame row vectors (forces, positions) -> source frame."""
        import numpy as np

        return np.asarray(vectors, dtype=float) @ np.asarray(self.rotation)

    def as_dict(self) -> dict:
        return {
            "vectors_nm": [list(v) for v in self.vectors_nm],
            "rotation": [list(r) for r in self.rotation],
            "lattice_change": [list(r) for r in self.lattice_change],
            "basis_flipped": self.basis_flipped,
            "tie_adjustment_nm": self.tie_adjustment_nm,
            "convention": (
                "engine = source @ rotation.T; source = engine @ rotation "
                "(ASE Cell.standard_form, then OpenMM reducePeriodicBoxVectors)"
            ),
        }


def box_violations(a: Sequence[float], b: Sequence[float], c: Sequence[float]) -> list[str]:
    """OpenMM's box-vector checks, with its exact comparisons.

    Mirrors ``System::setDefaultPeriodicBoxVectors`` and
    ``ContextImpl::setPeriodicBoxVectors`` (OpenMM 8.x): an empty list means
    OpenMM accepts ``(a, b, c)``.
    """
    problems = []
    if a[1] != 0.0 or a[2] != 0.0:
        problems.append("First periodic box vector must be parallel to x.")
    if b[2] != 0.0:
        problems.append("Second periodic box vector must be in the x-y plane.")
    if (
        a[0] <= 0.0
        or b[1] <= 0.0
        or c[2] <= 0.0
        or a[0] < 2 * abs(b[0])
        or a[0] < 2 * abs(c[0])
        or b[1] < 2 * abs(c[1])
    ):
        problems.append("Periodic box vectors must be in reduced form.")
    return problems


def reduced_box(cell) -> ReducedBox:
    """The OpenMM box for a fully periodic ASE cell (rows = lattice vectors, Angstrom).

    1. A left-handed basis has its third vector negated (the lattice is
       unchanged; only the basis orientation is).
    2. ``Cell.standard_form()`` rotates the cell to lower-triangular form
       with a positive diagonal (``rcell @ Q = cell``); its upper triangle is
       exactly zero.
    3. The vectors (in nm) are reduced exactly as OpenMM's
       ``reducePeriodicBoxVectors`` does: ``c -= b*round(c_y/b_y)``,
       ``c -= a*round(c_x/a_x)``, ``b -= a*round(b_x/a_x)``.
    4. An inequality broken only by rounding at an exact tie is restored by
       setting that component to exactly ``+-a_x/2`` (or ``+-b_y/2``).
    5. OpenMM's checks (:func:`box_violations`) are re-run exactly; anything
       left is a :class:`~nio_md_prep.mlip.errors.ConfigError`.
    """
    import numpy as np
    from ase.cell import Cell

    matrix = np.array(cell, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ConfigError(f"a periodic OpenMM box needs a finite 3x3 cell; got {cell!r}")
    determinant = float(np.linalg.det(matrix))
    scale = float(np.prod(np.linalg.norm(matrix, axis=1)))
    if scale == 0.0 or abs(determinant) <= 1e-10 * scale:
        raise ConfigError(
            "the cell is degenerate (zero volume); OpenMM needs three linearly "
            "independent box vectors for a periodic system"
        )
    flipped = determinant < 0.0
    basis = matrix.copy()
    if flipped:
        basis[2] *= -1.0
    rcell, rotation = Cell(basis).standard_form()
    lower = np.array(rcell, dtype=float)
    lower[0, 1] = lower[0, 2] = lower[1, 2] = 0.0
    a, b, c = ([float(v) / OPENMMML_LENGTH_SCALE_A_PER_NM for v in row] for row in lower)
    if a[0] <= 0.0 or b[1] <= 0.0 or c[2] <= 0.0:  # pragma: no cover - standard_form guarantees
        raise ConfigError(f"standard form of the cell has a non-positive diagonal: {lower}")

    def minus(u, v, k):
        return [u[i] - k * v[i] for i in range(3)]

    c = minus(c, b, round(c[1] / b[1]))
    c = minus(c, a, round(c[0] / a[0]))
    b = minus(b, a, round(b[0] / a[0]))
    adjustment = 0.0
    for vector, component, reference in ((b, 0, a[0]), (c, 0, a[0]), (c, 1, b[1])):
        excess = 2 * abs(vector[component]) - reference
        if excess > 0.0:
            if excess > BOX_TIE_TOLERANCE * reference:
                raise ConfigError(
                    f"could not reduce the cell to OpenMM's form (component "
                    f"{component} of a box vector exceeds half the reference by "
                    f"{excess:.3e} nm): {[a, b, c]}"
                )
            tied = math.copysign(reference / 2.0, vector[component])
            adjustment = max(adjustment, abs(tied - vector[component]))
            vector[component] = tied
    problems = box_violations(a, b, c)
    if problems:  # pragma: no cover - the steps above establish every inequality
        raise ConfigError(f"OpenMM would reject the reduced box {[a, b, c]}: {problems}")
    reduced_angstrom = np.array([a, b, c]) * OPENMMML_LENGTH_SCALE_A_PER_NM
    change = reduced_angstrom @ np.linalg.inv(matrix @ rotation.T)
    lattice_change = np.rint(change)
    if np.abs(change - lattice_change).max() > 1e-6 or round(abs(np.linalg.det(lattice_change))) != 1:
        raise ConfigError(  # pragma: no cover - reduction is an integer basis change
            f"box reduction changed the lattice (basis change {change.tolist()})"
        )
    return ReducedBox(
        vectors_nm=tuple(tuple(float(x) for x in v) for v in (a, b, c)),
        rotation=tuple(tuple(float(x) for x in row) for row in rotation),
        lattice_change=tuple(tuple(int(x) for x in row) for row in lattice_change),
        basis_flipped=flipped,
        tie_adjustment_nm=adjustment,
    )


def device_ordinal(device: str | None) -> int | None:
    """``"cuda:1"`` -> ``1``; ``"cuda"``, ``"cpu"``, ``"mps"`` -> ``None``."""
    if not device or ":" not in device:
        return None
    return int(device.split(":", 1)[1])


def platform_properties(
    platform: str | None,
    *,
    platform_precision: str | None = None,
    threads: int | None = None,
    device: str | None = None,
) -> dict[str, str]:
    """The OpenMM platform properties a request maps to, or a ConfigError.

    - ``platform_precision`` -> ``Precision``: CUDA/OpenCL/HIP only (the CPU
      platform has only ``Threads``/``DeterministicForces``, Reference none;
      OpenMM raises "Illegal property name" otherwise).
    - ``threads`` -> ``Threads``: the CPU platform only.
    - ``device = "cuda:N"`` -> ``DeviceIndex = "N"`` on a GPU platform, so the
      OpenMM context and the MACE model use the same GPU. A bare ``"cuda"``
      leaves ``DeviceIndex`` to OpenMM (the value it picks is recorded).

    With ``platform=None`` OpenMM chooses the platform, so no property can be
    passed and any of these requests is refused.
    """
    properties: dict[str, str] = {}
    if platform is None:
        given = [
            name
            for name, value in (
                ("engine.platform_precision", platform_precision),
                ("engine.threads", threads),
            )
            if value is not None
        ]
        if given:
            raise ConfigError(
                f"{' and '.join(given)} map to OpenMM platform properties, which can "
                "only be set on a named platform; set engine.platform (e.g. 'CPU' for "
                "threads, 'CUDA' or 'OpenCL' for platform_precision)"
            )
        return properties
    if platform_precision is not None:
        if platform not in OPENMM_GPU_PLATFORMS:
            raise ConfigError(
                f"engine.platform_precision sets OpenMM's 'Precision' property, which the "
                f"{platform} platform does not have (only {', '.join(OPENMM_GPU_PLATFORMS)} "
                "do); remove it, or choose one of those platforms"
            )
        properties["Precision"] = platform_precision
    if threads is not None:
        if platform != "CPU":
            raise ConfigError(
                f"engine.threads maps to the CPU platform's 'Threads' property; the "
                f"{platform} platform has no such property. Remove engine.threads or "
                "use engine.platform = 'CPU'."
            )
        properties["Threads"] = str(int(threads))
    ordinal = device_ordinal(device)
    if ordinal is not None and platform in OPENMM_GPU_PLATFORMS:
        properties["DeviceIndex"] = str(ordinal)
    return properties


def check_request(engine_spec, *, device: str | None) -> dict[str, str]:
    """Validate-time check of the engine request; returns the properties it maps to."""
    return platform_properties(
        engine_spec.platform,
        platform_precision=engine_spec.platform_precision,
        threads=engine_spec.threads,
        device=device,
    )


def resolved_thermostat(simulation: SimulationSpec) -> str | None:
    """The thermostat this engine will run for ``simulation`` (``None`` for NVE)."""
    if simulation.ensemble != "nvt":
        return None
    return simulation.thermostat or DEFAULT_THERMOSTAT


def check_md_request(simulation: SimulationSpec) -> None:
    """Refuse MD settings this route cannot honour, before anything runs."""
    if simulation.task != "md":
        return
    if simulation.ensemble == "npt":
        raise CapabilityError(NPT_POLICY)
    thermostat = resolved_thermostat(simulation)
    if thermostat is not None and thermostat not in SUPPORTED_THERMOSTATS:
        raise ConfigError(
            f"thermostat {thermostat!r} is not implemented by the OpenMM engine; "
            f"supported: {', '.join(f'{k} ({v})' for k, v in SUPPORTED_THERMOSTATS.items())}"
        )


def check_structure(atoms, *, md: bool) -> ReducedBox | None:
    """Refuse geometry OpenMM-ML cannot represent; return the box it would use."""
    from ..structures import fixed_atom_indices, periodic_axes

    pbc = periodic_axes(atoms)
    if any(pbc) and not all(pbc):
        raise CapabilityError(
            f"the structure is periodic along cell axes {[a for a, p in zip('abc', pbc) if p]} "
            "only; OpenMM-ML's MACE force is periodic along all three axes or none, so a "
            "slab or wire would silently become a 3D crystal or a cluster. Use the ASE or "
            "LAMMPS engine, or make the structure fully periodic or fully non-periodic."
        )
    fixed = fixed_atom_indices(atoms)
    if md and len(fixed) == len(atoms):
        raise ConfigError(
            "every atom is frozen by FixAtoms; there is nothing for the dynamics to move"
        )
    return reduced_box(atoms.get_cell()) if all(pbc) else None


def derive_seeds(seed: int) -> dict[str, int]:
    """Independent child seeds for velocity initialisation and the integrator.

    OpenMM's ``setVelocitiesToTemperature`` and a Langevin integrator given
    the *same* seed may draw from identically seeded generators; spawning
    two children of one ``SeedSequence`` keeps them independent while the
    run stays reproducible from the one recorded ``seed``. Every child lies
    in ``1..MAX_SEED`` (OpenMM reads 0 as "pick one at random").
    """
    import numpy as np

    children = np.random.SeedSequence(int(seed)).spawn(2)
    velocity, integrator = (
        int(child.generate_state(1, dtype=np.uint32)[0]) % MAX_SEED + 1 for child in children
    )
    return {"velocity_seed": velocity, "integrator_seed": integrator}


def energy_to_eV(energy_kj_per_mol: float) -> float:
    """An OpenMM-ML MACE potential energy (kJ/mol) back to the model's eV."""
    return float(energy_kj_per_mol) / OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV


def forces_to_eV_per_A(forces_kj_per_mol_nm):
    """OpenMM-ML MACE forces (kJ/mol/nm) back to the model's eV/Angstrom."""
    import numpy as np

    return np.asarray(forces_kj_per_mol_nm, dtype=float) / (
        OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV * OPENMMML_LENGTH_SCALE_A_PER_NM
    )


def kinetic_to_eV(energy_kj_per_mol: float) -> float:
    """An OpenMM kinetic (or heat-bath) energy, which is exactly kJ/mol, in eV."""
    return float(energy_kj_per_mol) * KJ_PER_MOL_IN_EV


#: DCD layout written by ``openmm.app.DCDFile`` (little-endian Fortran records).
_DCD_HEADER_BYTES = 276
_DCD_BOX_RECORD_BYTES = 56


def dcd_frame_count(path, n_atoms: int) -> int:
    """Frames in a DCD file, read from its header and cross-checked by size.

    The header's ``NSET`` (int32 at byte 8) is what ``DCDFile`` updates on
    every frame; the file length must equal the header plus that many
    frames of ``n_atoms`` coordinates (with a unit-cell record when the
    header's box flag is set). A missing, truncated or inconsistent file is a
    :class:`~nio_md_prep.mlip.errors.ResultError`, never "zero frames".
    """
    path = Path(path)
    if not path.exists():
        raise ResultError(f"the trajectory {path} was not written")
    data = path.read_bytes()
    if len(data) < _DCD_HEADER_BYTES or data[4:8] != b"CORD" or struct.unpack("<i", data[:4])[0] != 84:
        raise ResultError(f"the trajectory {path} is not a readable DCD file")
    frames = struct.unpack("<i", data[8:12])[0]
    box_flag = struct.unpack("<i", data[44:48])[0]
    atoms_in_file = struct.unpack("<i", data[268:272])[0]
    if atoms_in_file != n_atoms:
        raise ResultError(
            f"the trajectory {path} describes {atoms_in_file} atoms; the run had {n_atoms}"
        )
    frame_bytes = (_DCD_BOX_RECORD_BYTES if box_flag else 0) + 3 * (8 + 4 * n_atoms)
    expected_size = _DCD_HEADER_BYTES + frames * frame_bytes
    if frames < 0 or len(data) != expected_size:
        raise ResultError(
            f"the trajectory {path} is truncated or corrupt: its header says {frames} "
            f"frames ({expected_size} bytes) but the file has {len(data)} bytes"
        )
    return frames


# ---------------------------------------------------------------------------
# OpenMM-side helpers (lazy imports)
# ---------------------------------------------------------------------------


def openmmml_version() -> str | None:
    from importlib import metadata

    try:
        return metadata.version("openmmml")
    except metadata.PackageNotFoundError:
        return None


def _version_tuple(version: str) -> tuple[int, ...]:
    parts = []
    for piece in version.split(".")[:3]:
        digits = "".join(ch for ch in piece if ch.isdigit())
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


def require_openmm() -> dict:
    """Import OpenMM and OpenMM-ML, refusing releases with other device semantics.

    Returns ``{"openmm": version, "openmmml": version}``.
    """
    try:
        import openmm
        import openmmml  # noqa: F401
    except ImportError as exc:
        raise MissingDependencyError(
            "the OpenMM engine",
            ("openmm", "openmmml"),
            hint="Install the OpenMM extra: pip install 'nio-md-prep[openmm]'",
        ) from exc
    version = openmmml_version()
    if version is not None:
        supported = _version_tuple(version) >= MIN_OPENMMML_VERSION
    else:
        # No distribution metadata (a source checkout): fall back to the
        # feature that defines the >=1.6 semantics.
        supported = hasattr(openmm, "PythonForce")
    if not supported:
        raise MissingDependencyError(
            "the OpenMM engine",
            ("openmmml>=1.6",),
            hint=(
                f"OpenMM-ML {version or 'unknown'} is installed. Releases before 1.6 "
                "evaluate MACE through a TorchForce on the Context's device and ignore "
                "the device= argument, so potential.device could not be honoured. "
                "Upgrade: pip install 'openmmml>=1.6'."
            ),
        )
    return {"openmm": getattr(openmm, "__version__", None), "openmmml": version}


def available_platforms() -> list[str]:
    require_openmm()
    import openmm

    return [
        openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())
    ]


def resolve_platform(engine_spec, *, device: str | None):
    """``(Platform or None, properties)`` for this request; unknown -> ConfigError."""
    require_openmm()
    import openmm

    properties = check_request(engine_spec, device=device)
    if engine_spec.platform is None:
        return None, properties
    available = available_platforms()
    if engine_spec.platform not in available:
        raise ConfigError(
            f"OpenMM platform {engine_spec.platform!r} is not available on this machine; "
            f"available platforms: {', '.join(available)}"
        )
    return openmm.Platform.getPlatformByName(engine_spec.platform), properties


def describe_platform(context, requested: dict[str, str]) -> dict:
    """The platform a Context actually runs on, with every property's value."""
    platform = context.getPlatform()
    return {
        "name": platform.getName(),
        "properties": {
            name: platform.getPropertyValue(context, name) for name in platform.getPropertyNames()
        },
        "requested_properties": dict(requested),
    }


def check_torch_device(device: str) -> None:
    """Refuse a MACE device torch cannot use, before OpenMM-ML loads the model."""
    if not device.startswith(("cuda", "mps")):
        return
    import torch

    if device.startswith("cuda"):
        count = torch.cuda.device_count() if torch.cuda.is_available() else 0
        ordinal = device_ordinal(device) or 0
        if ordinal >= count:
            raise ConfigError(
                f"potential.device = {device!r}, but torch sees {count} CUDA device(s) "
                "on this machine"
            )
    elif not (getattr(torch.backends, "mps", None) and torch.backends.mps.is_available()):
        raise ConfigError("potential.device = 'mps', but torch reports MPS unavailable")


def build_topology(atoms, box: ReducedBox | None = None):
    """Build an OpenMM ``Topology`` from an ASE ``Atoms``.

    One chain, one residue, one atom per site. That is sufficient for a
    full-system MLIP -- the model reads elements and coordinates, not
    residues -- and it avoids inventing a topology this layer does not have.
    A periodic structure gets the reduced box (:func:`reduced_box`), never
    the raw ASE cell rows; a partially periodic one is refused.
    """
    require_openmm()
    from openmm import Vec3
    from openmm.app import Element, Topology
    from openmm.unit import nanometer

    from ..structures import periodic_axes

    pbc = periodic_axes(atoms)
    if any(pbc) and not all(pbc):
        check_structure(atoms, md=False)
    if all(pbc) and box is None:
        box = reduced_box(atoms.get_cell())
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MLIP", chain)
    for symbol in atoms.get_chemical_symbols():
        topology.addAtom(symbol, Element.getBySymbol(symbol), residue)
    if box is not None:
        topology.setPeriodicBoxVectors([Vec3(*v) for v in box.vectors_nm] * nanometer)
    return topology


def positions_nm(atoms, box: ReducedBox | None = None):
    """Positions in nm, rotated into the engine frame when a box is used."""
    from openmm import Vec3
    from openmm.unit import nanometer

    positions = atoms.get_positions()
    if box is not None:
        positions = box.to_engine(positions)
    return [
        Vec3(*(float(c) / OPENMMML_LENGTH_SCALE_A_PER_NM for c in row)) for row in positions
    ] * nanometer


def build_system(
    atoms,
    model_path,
    *,
    energy_convention: str,
    potential_name: str = "mace",
    precision: str = "float64",
    device: str = "cpu",
    remove_cm_motion: bool = False,
    fixed_atoms: Sequence[int] = (),
    charge: float = 0.0,
    multiplicity: float = 1.0,
    box: ReducedBox | None = None,
):
    """Create the OpenMM ``System`` for a model, with every argument explicit.

    ``energy_convention`` -> ``returnEnergyType``, ``precision`` (the
    potential's) -> OpenMM-ML ``precision='single'|'double'``, ``device`` ->
    ``device=``, ``charge``/``multiplicity`` -> MACE ``total_charge`` /
    ``total_spin``. Particle masses are replaced by the ASE masses, and
    ``fixed_atoms`` get zero mass. ``remove_cm_motion`` adds OpenMM-ML's
    ``CMMotionRemover``. Returns ``(system, topology, record)`` where
    ``record`` says what was passed.
    """
    versions = require_openmm()
    from openmm import unit
    from openmmml import MLPotential

    if energy_convention not in ENERGY_TYPE:
        raise ConfigError(
            f"energy convention {energy_convention!r} has no OpenMM-ML equivalent; "
            f"expected one of {', '.join(ENERGY_TYPE)}"
        )
    if precision not in MODEL_PRECISION:
        raise ConfigError(f"precision {precision!r} has no OpenMM-ML equivalent")
    check_torch_device(device)
    topology = build_topology(atoms, box)
    arguments = {
        "returnEnergyType": ENERGY_TYPE[energy_convention],
        "precision": MODEL_PRECISION[precision],
        "device": device,
        "charge": float(charge),
        "multiplicity": float(multiplicity),
    }
    potential = MLPotential(potential_name, modelPath=str(model_path))
    system = potential.createSystem(topology, removeCMMotion=bool(remove_cm_motion), **arguments)
    if system.getNumParticles() != len(atoms):
        raise ResultError(
            f"OpenMM-ML built a System with {system.getNumParticles()} particles for "
            f"{len(atoms)} atoms"
        )
    fixed = sorted({int(i) for i in fixed_atoms})
    masses = atoms.get_masses()
    for index, mass in enumerate(masses):
        system.setParticleMass(index, 0.0 if index in fixed else float(mass) * unit.dalton)
    record = {
        "createSystem": {"removeCMMotion": bool(remove_cm_motion), **arguments},
        "masses": "ASE masses (atoms.get_masses()); FixAtoms atoms zero",
        "zero_mass_atoms": fixed,
        "periodic": box is not None,
        "box": box.as_dict() if box is not None else None,
        "versions": versions,
    }
    return system, topology, record


def make_context(
    system,
    atoms,
    engine_spec,
    *,
    box: ReducedBox | None = None,
    device: str | None = None,
    timestep_fs: float = 1.0,
):
    """A Context on the requested platform with the requested properties."""
    require_openmm()
    import openmm
    from openmm import Vec3, unit

    platform, properties = resolve_platform(engine_spec, device=device)
    integrator = openmm.VerletIntegrator(timestep_fs * 0.001 * unit.picosecond)
    context = (
        openmm.Context(system, integrator, platform, properties)
        if platform is not None
        else openmm.Context(system, integrator)
    )
    if box is not None:
        context.setPeriodicBoxVectors(*[Vec3(*v) for v in box.vectors_nm])
    context.setPositions(positions_nm(atoms, box))
    return context, integrator, describe_platform(context, properties)


def singlepoint(atoms, system, engine_spec, *, box: ReducedBox | None = None, device=None) -> dict:
    """Evaluate one geometry and convert straight out of OpenMM's units.

    Forces are rotated back into the source frame; the energy is invariant.
    """
    require_openmm()
    import numpy as np
    from openmm import unit

    started = time.perf_counter()
    context, integrator, platform_record = make_context(
        system, atoms, engine_spec, box=box, device=device
    )
    try:
        state = context.getState(getEnergy=True, getForces=True)
        energy_kj_per_mol = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces_native = np.asarray(
            state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer),
            dtype=float,
        )
    finally:
        del context, integrator
    if forces_native.shape != (len(atoms), 3):
        raise ResultError(
            f"OpenMM returned forces of shape {forces_native.shape} for {len(atoms)} atoms"
        )
    forces = forces_to_eV_per_A(forces_native)
    if box is not None:
        forces = box.to_source(forces)
    return {
        "energy_eV": energy_to_eV(energy_kj_per_mol),
        "forces_eV_per_A": forces.tolist(),
        "stress_eV_per_A3": None,
        "per_atom_energy_eV": None,
        "symbols": list(atoms.get_chemical_symbols()),
        "trajectory_path": None,
        "wall_time_s": time.perf_counter() - started,
        "native": {
            "energy_kJ_per_mol": energy_kj_per_mol,
            "units": OPENMM.name,
            "energy_scale_kJ_per_mol_per_eV": OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV,
            "length_scale_A_per_nm": OPENMMML_LENGTH_SCALE_A_PER_NM,
            "platform": platform_record,
            "box": box.as_dict() if box is not None else None,
            "forces_frame": "source (rotated back from the OpenMM box frame)"
            if box is not None
            else "source",
        },
    }


def make_integrator(simulation: SimulationSpec, seeds: dict[str, int]):
    """``(integrator, resolved)`` for the requested ensemble and thermostat."""
    require_openmm()
    import openmm
    from openmm import unit

    check_md_request(simulation)
    timestep = simulation.timestep_fs * 0.001 * unit.picosecond
    thermostat = resolved_thermostat(simulation)
    resolved: dict = {
        "ensemble": simulation.ensemble,
        "thermostat_requested": simulation.thermostat,
        "thermostat": thermostat,
        "timestep_fs": simulation.timestep_fs,
    }
    if thermostat is None:
        integrator = openmm.VerletIntegrator(timestep)
    else:
        damping_fs = simulation.resolved_thermostat_damping_fs
        rate_per_ps = 1000.0 / damping_fs
        temperature = simulation.temperature_K * unit.kelvin
        resolved.update(temperature_K=simulation.temperature_K, damping_fs=damping_fs)
        if thermostat == "langevin":
            integrator = openmm.LangevinMiddleIntegrator(
                temperature, rate_per_ps / unit.picosecond, timestep
            )
            integrator.setRandomNumberSeed(seeds["integrator_seed"])
            resolved.update(
                friction_per_ps=rate_per_ps, integrator_seed=seeds["integrator_seed"]
            )
        else:
            integrator = openmm.NoseHooverIntegrator(
                temperature, rate_per_ps / unit.picosecond, timestep
            )
            resolved.update(collision_frequency_per_ps=rate_per_ps)
    resolved["name"] = type(integrator).__name__
    return integrator, resolved


class _FrameReporter:
    """Writes DCD frames through ``openmm.app.DCDFile`` on a handle we own."""

    def __init__(self, handle, topology, timestep_fs: float, interval: int) -> None:
        from openmm import unit
        from openmm.app import DCDFile

        self.interval = interval
        self.count = 0
        self._dcd = DCDFile(
            handle, topology, timestep_fs * 0.001 * unit.picosecond, firstStep=0, interval=interval
        )

    def describeNextReport(self, simulation):  # noqa: N802 - OpenMM reporter API
        steps = self.interval - simulation.currentStep % self.interval
        return {"steps": steps, "periodic": False, "include": ["positions"]}

    def report(self, simulation, state):
        self._dcd.writeModel(
            state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors()
        )
        self.count += 1


class _SeriesReporter:
    """Samples step, energies and temperature; writes the log on a handle we own."""

    COLUMNS = "# step time_fs potential_eV kinetic_eV total_eV temperature_K conserved_eV\n"

    def __init__(self, handle, interval: int, *, timestep_fs: float, ndof: int, heat_bath: bool):
        self.handle = handle
        self.interval = interval
        self.timestep_fs = timestep_fs
        self.ndof = ndof
        self.heat_bath = heat_bath
        self.rows: list[tuple] = []
        handle.write(self.COLUMNS)

    def describeNextReport(self, simulation):  # noqa: N802 - OpenMM reporter API
        steps = self.interval - simulation.currentStep % self.interval
        return {"steps": steps, "periodic": None, "include": ["energy"]}

    def report(self, simulation, state):
        from openmm import unit

        step = int(simulation.currentStep)
        if self.rows and self.rows[-1][0] == step:
            return
        potential_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        kinetic_kj = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        potential = energy_to_eV(potential_kj)
        kinetic = kinetic_to_eV(kinetic_kj)
        temperature = 2.0 * kinetic_kj / (self.ndof * MOLAR_GAS_CONSTANT_KJ_PER_MOL_K)
        conserved = None
        if self.heat_bath:
            # The chain starts at rest; OpenMM 8.5 Reference/CPU crash when
            # computeHeatBathEnergy() is called before the first step.
            bath_kj = (
                0.0
                if step == 0
                else simulation.integrator.computeHeatBathEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                )
            )
            conserved = potential + kinetic + kinetic_to_eV(bath_kj)
        row = (step, step * self.timestep_fs, potential, kinetic, potential + kinetic,
               temperature, conserved)
        self.rows.append(row)
        self.handle.write(
            f"{step} {row[1]:.6f} {potential:.12f} {kinetic:.12f} {potential + kinetic:.12f} "
            f"{temperature:.6f} {'' if conserved is None else f'{conserved:.12f}'}\n"
        )


def _remove_com_velocity(context, system) -> None:
    """Subtract the mass-weighted mean velocity of the massive particles."""
    import numpy as np
    from openmm import unit

    masses = np.array(
        [system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles())]
    )
    velocities = np.asarray(
        context.getState(getVelocities=True)
        .getVelocities(asNumpy=True)
        .value_in_unit(unit.nanometer / unit.picosecond)
    )
    massive = masses > 0.0
    com = (masses[massive, None] * velocities[massive]).sum(axis=0) / masses[massive].sum()
    velocities[massive] -= com
    context.setVelocities(velocities * (unit.nanometer / unit.picosecond))


def run_md(
    atoms,
    system,
    engine_spec,
    simulation: SimulationSpec,
    *,
    workdir: Path,
    box: ReducedBox | None = None,
    device: str | None = None,
    fixed_atoms: Sequence[int] = (),
    remove_cm_motion: bool = False,
) -> dict:
    """Run a short diagnostic trajectory through OpenMM.

    ``system`` must have been built for ``atoms`` with the same ``box``,
    ``fixed_atoms`` (zero masses) and ``remove_cm_motion``. Returns the
    payload :meth:`~nio_md_prep.mlip.bridges.base.Bridge._trajectory`
    expects: ``steps_completed`` from ``simulation.currentStep``, frames
    counted from the DCD header (step 0 included, like the ASE engine, so
    ``steps // interval + 1``; any mismatch is a ResultError), final
    positions rotated back to the source frame, an energy series, and
    temperatures computed with explicit degrees of freedom
    (``3N - 3 n_fixed - 3`` when a CMMotionRemover is present).
    """
    versions = require_openmm()
    import numpy as np
    from openmm import Vec3, app, unit

    from ..diagnostics import temperature_ndof
    from ..specs import resolve_seed

    if simulation.ensemble == "npt":
        raise CapabilityError(NPT_POLICY)
    check_md_request(simulation)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    n_atoms = len(atoms)
    fixed = sorted({int(i) for i in fixed_atoms})
    ndof = temperature_ndof(n_atoms, n_fixed=len(fixed), com_removed=bool(remove_cm_motion))

    seed = resolve_seed(simulation.seed)
    seeds = derive_seeds(seed)
    integrator, resolved = make_integrator(simulation, seeds)
    platform, properties = resolve_platform(engine_spec, device=device)
    topology = build_topology(atoms, box)
    simulation_object = (
        app.Simulation(topology, system, integrator, platform, properties)
        if platform is not None
        else app.Simulation(topology, system, integrator)
    )
    context = simulation_object.context
    if box is not None:
        context.setPeriodicBoxVectors(*[Vec3(*v) for v in box.vectors_nm])
    context.setPositions(positions_nm(atoms, box))
    if simulation.temperature_K:
        context.setVelocitiesToTemperature(
            simulation.temperature_K * unit.kelvin, seeds["velocity_seed"]
        )
        velocities = f"setVelocitiesToTemperature({simulation.temperature_K} K, velocity_seed)"
        if remove_cm_motion:
            # setVelocitiesToTemperature leaves a random centre-of-mass
            # velocity that the CMMotionRemover would delete at the first step:
            # a one-off kinetic-energy drop that an NVE drift would misread.
            _remove_com_velocity(context, system)
            velocities += ", centre-of-mass velocity removed (CMMotionRemover present)"
    else:
        velocities = "zero"
    resolved.update(
        seed=seed,
        seed_drawn=simulation.seed is None,
        velocity_seed=seeds["velocity_seed"] if simulation.temperature_K else None,
        velocities=velocities,
        removeCMMotion=bool(remove_cm_motion),
        temperature_ndof=ndof,
        platform=describe_platform(context, properties),
        versions=versions,
    )
    if resolved["name"] == "NoseHooverIntegrator":
        chain = integrator.getThermostat(0)
        resolved.update(
            chain_length=chain.getChainLength(),
            num_multi_time_steps=chain.getNumMultiTimeSteps(),
            num_yoshida_suzuki=chain.getNumYoshidaSuzukiTimeSteps(),
            openmm_thermostat_ndof=chain.getNumDegreesOfFreedom(),
        )

    trajectory_path = workdir / "smoke_md.dcd"
    log_path = workdir / "smoke_md.log"
    frame_interval = max(1, simulation.trajectory_interval)
    started = time.perf_counter()
    with open(trajectory_path, "wb") as dcd_handle, open(log_path, "w", encoding="utf-8") as log:
        frames = _FrameReporter(dcd_handle, topology, simulation.timestep_fs, frame_interval)
        series = _SeriesReporter(
            log,
            max(1, simulation.log_interval),
            timestep_fs=simulation.timestep_fs,
            ndof=ndof,
            heat_bath=resolved["name"] == "NoseHooverIntegrator",
        )
        initial_state = context.getState(getPositions=True, getEnergy=True)
        frames.report(simulation_object, initial_state)
        series.report(simulation_object, initial_state)
        simulation_object.reporters.extend([frames, series])
        try:
            simulation_object.step(simulation.steps)
        finally:
            simulation_object.reporters.clear()
        steps_completed = int(simulation_object.currentStep)
        final_state = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
        series.report(simulation_object, final_state)
    wall_time = time.perf_counter() - started

    frames_written = dcd_frame_count(trajectory_path, n_atoms)
    expected_frames = steps_completed // frame_interval + 1
    if frames_written != frames.count or frames_written != expected_frames:
        raise ResultError(
            f"{trajectory_path} holds {frames_written} frames ({frames.count} written); "
            f"{steps_completed} steps written every {frame_interval} steps (including "
            f"step 0) should give {expected_frames}"
        )

    positions = np.asarray(final_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    velocities_nm_ps = np.asarray(
        final_state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
    )
    positions = positions * OPENMMML_LENGTH_SCALE_A_PER_NM
    velocities_A_fs = velocities_nm_ps * OPENMMML_LENGTH_SCALE_A_PER_NM / 1000.0
    if box is not None:
        positions = box.to_source(positions)
        velocities_A_fs = box.to_source(velocities_A_fs)
    source = atoms.get_positions()
    if fixed:
        # OpenMM never moves a massless particle; the rotation round trip
        # costs ~1e-15 Angstrom, so restore the exact input coordinates.
        moved = np.abs(positions[fixed] - source[fixed]).max()
        if moved > 1e-9:
            raise ResultError(
                f"a zero-mass (FixAtoms) atom moved by {moved:.3e} Angstrom in OpenMM"
            )
        positions[fixed] = source[fixed]
        velocities_A_fs[fixed] = 0.0
    from ase import units as ase_units

    final_atoms = atoms.copy()
    final_atoms.set_positions(positions, apply_constraint=False)
    final_atoms.set_momenta(
        atoms.get_masses()[:, None] * velocities_A_fs / ase_units.fs, apply_constraint=False
    )

    rows = series.rows
    temperatures = [row[5] for row in rows]
    energy_series = {
        "time_fs": [row[1] for row in rows],
        "total_energy_eV": [row[4] for row in rows],
    }
    if resolved["name"] == "NoseHooverIntegrator":
        energy_series["conserved_energy_eV"] = [row[6] for row in rows]
        energy_series["conserved_quantity"] = (
            "potential + kinetic + NoseHooverIntegrator.computeHeatBathEnergy()"
        )
    return {
        "atoms": final_atoms,
        "steps_completed": steps_completed,
        "frames_written": frames_written,
        "trajectory_path": str(trajectory_path),
        "trajectory_frame": (
            "DCD positions are unwrapped, in the OpenMM box frame (engine = source @ "
            "rotation.T); the DCD unit cell is the reduced box"
            if box is not None
            else "DCD positions are unwrapped, in the source frame"
        ),
        "log_path": str(log_path),
        "final_positions": positions.tolist(),
        "final_cell": atoms.get_cell().tolist(),
        "final_pbc": [bool(v) for v in atoms.pbc],
        "energy_series": energy_series,
        "temperature_start_K": temperatures[0],
        "temperature_end_K": temperatures[-1],
        "max_temperature_K": max(temperatures),
        "temperature_ndof": ndof,
        "constraints": (
            {
                "fixed_atoms": fixed,
                "n_fixed": len(fixed),
                "method": "zero particle mass (OpenMM does not move massless particles)",
            }
            if fixed
            else None
        ),
        "wall_time_s": wall_time,
        "integrator": resolved["name"],
        "integrator_resolved": resolved,
    }


__all__ = [
    "ENERGY_TYPE",
    "GPU_PLATFORMS",
    "OPENMMML_ENERGY_SCALE_KJ_PER_MOL_PER_EV",
    "OPENMMML_LENGTH_SCALE_A_PER_NM",
    "MIN_OPENMMML_VERSION",
    "MODEL_PRECISION",
    "SUPPORTED_THERMOSTATS",
    "DEFAULT_THERMOSTAT",
    "NPT_POLICY",
    "OpenMMEngine",
    "ReducedBox",
    "box_violations",
    "reduced_box",
    "device_ordinal",
    "platform_properties",
    "check_request",
    "check_md_request",
    "check_structure",
    "resolved_thermostat",
    "derive_seeds",
    "dcd_frame_count",
    "energy_to_eV",
    "forces_to_eV_per_A",
    "kinetic_to_eV",
    "require_openmm",
    "available_platforms",
    "resolve_platform",
    "describe_platform",
    "build_topology",
    "positions_nm",
    "build_system",
    "make_context",
    "make_integrator",
    "singlepoint",
    "run_md",
]
