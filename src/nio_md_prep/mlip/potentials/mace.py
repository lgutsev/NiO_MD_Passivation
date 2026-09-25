"""The MACE potential family, described independently of ASE/LAMMPS/OpenMM.

Everything torch-flavoured is behind a lazy import. :meth:`MaceAdapter.inspect`
is written so that a machine without torch still gets a useful report -- the
declared metadata, plus an explicit statement that nothing was read.

A trained MACE model carries its own element table, cutoff, heads and
per-element atomic reference energies (E0, one row per head). Those are
exactly the things this subsystem needs and cannot safely guess, so when they
*are* discoverable they are cross-checked against what the configuration
declares (:meth:`MaceAdapter.model_problems`), and a disagreement is refused
rather than silently resolved in either direction.

**One load per model file.** :meth:`MaceAdapter.discover` loads the model on
the CPU (so validating a CUDA job needs no GPU) and caches the metadata --
never the model object -- by the file's SHA256 at module level, so every
adapter, bridge and capability query in a process shares one ``torch.load``
per file. The SHA256 itself is cached by ``(path, size, mtime)``.
"""
from __future__ import annotations

from pathlib import Path

from ..capabilities import CapabilitySet
from ..environment import Availability, probe
from ..errors import ConfigError, MissingDependencyError
from ..specs import MacePotentialSpec
from .base import PotentialAdapter

#: A declared E0 further than this from the model's own (eV) is a conflict.
E0_CONFLICT_TOLERANCE_EV = 1e-8
#: The ``mode`` values ``torch.compile`` accepts, which MACE's ASE calculator
#: forwards unchanged (``mace/calculators/mace.py``, ``compile_mode``).
COMPILE_MODES = ("default", "reduce-overhead", "max-autotune", "max-autotune-no-cudagraphs")

_SHA256_CACHE: dict[tuple[str, int, int], str] = {}
_DISCOVERY_CACHE: dict[str, dict] = {}


def model_sha256(path) -> str:
    """SHA256 of a model file, cached by ``(resolved path, size, mtime_ns)``."""
    from ..provenance import sha256_file

    path = Path(path)
    stat = path.stat()
    key = (str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    digest = _SHA256_CACHE.get(key)
    if digest is None:
        digest = _SHA256_CACHE[key] = sha256_file(path)
    return digest


def clear_discovery_cache() -> None:
    """Forget every cached model hash and discovery (tests; a rewritten file)."""
    _SHA256_CACHE.clear()
    _DISCOVERY_CACHE.clear()


class MaceAdapter(PotentialAdapter):
    kind = "mace"
    requires = ("mace", "torch")

    spec: MacePotentialSpec

    def __init__(self, spec: MacePotentialSpec) -> None:
        super().__init__(spec)
        #: Why the last :meth:`discover` attempt failed, when it did.
        self.discovery_error: str | None = None

    def capabilities(self) -> CapabilitySet:
        """What a MACE model can produce, before any engine constrains it.

        MACE is a gradient-domain model: energy, forces, stress and a per-atom
        (site) energy decomposition all come out of the same forward pass, and
        periodic cells are supported natively. Whether a *route* exposes each
        of them is the bridge's claim, not this one. GPU support is a
        property of the requested device. When the model file can be read,
        the element set is the declared elements the model actually covers.
        """
        elements = frozenset(self.spec.elements)
        discovered = self._discovered_if_available()
        if discovered and discovered.get("elements"):
            elements &= frozenset(discovered["elements"])
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=True,
            per_atom_energy=True,
            periodic=True,
            gpu=self.spec.wants_gpu,
            elements=elements,
            precisions=frozenset({self.spec.precision}),
            engines=frozenset({"ase", "lammps", "openmm"}),
            native_energy_convention=self.spec.energy_convention,
            convertible_energy_conventions=(
                frozenset({"total", "interaction"})
                if self.atomic_reference_energies()
                else frozenset()
            ),
            notes=(
                f"MACE model {self.spec.model_path.name} on {self.spec.device} "
                f"in {self.spec.precision}",
            ),
        )

    def availability(self) -> Availability:
        if not self.spec.model_path.exists():
            return Availability(
                available=False,
                missing=(),
                detail=f"model file not found: {self.spec.model_path}",
            )
        return probe(self.requires, detail="MACE model evaluation")

    # -- the model file ---------------------------------------------------

    def model_sha256(self) -> str | None:
        """The model file's SHA256, or ``None`` when the file does not exist."""
        if not self.spec.model_path.exists():
            return None
        return model_sha256(self.spec.model_path)

    def load_model(self, map_location=None):
        """Load the torch model (on ``map_location``, default the spec's device)."""
        if not self.spec.model_path.exists():
            raise FileNotFoundError(f"MACE model not found: {self.spec.model_path}")
        try:
            import mace  # noqa: F401 - sets TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD for e3nn
            import torch
        except ImportError as exc:
            raise MissingDependencyError(
                f"loading the MACE model {self.spec.model_path.name}",
                ("mace", "torch"),
                hint="Install the MACE extra: pip install 'nio-md-prep[mace]'",
            ) from exc
        return torch.load(
            self.spec.model_path,
            map_location=map_location or self._torch_device(),
            weights_only=False,
        )

    def _torch_device(self):
        return self.spec.device

    def torch_dtype(self):
        import torch

        return {"float32": torch.float32, "float64": torch.float64}[self.spec.precision]

    def discover(self) -> dict:
        """Element table, cutoff, heads, E0s and dtype, read from the model file.

        Loaded once per file (on the CPU) and cached by SHA256; the returned
        mapping is a fresh copy with this spec's head resolved into
        ``head`` / ``atomic_reference_energies_eV`` (``None`` plus
        ``head_error`` when the head cannot be resolved).
        """
        digest = model_sha256(self.spec.model_path)
        facts = _DISCOVERY_CACHE.get(digest)
        if facts is None:
            model = self.load_model(map_location="cpu")
            facts = _describe_model(model)
            facts["sha256"] = digest
            _DISCOVERY_CACHE[digest] = facts
            del model
        report = {key: _copy(value) for key, value in facts.items()}
        try:
            head = resolve_head(report, self.spec.head)
        except ConfigError as exc:
            report.update(head=None, head_error=str(exc), atomic_reference_energies_eV=None)
        else:
            by_head = report.get("atomic_reference_energies_by_head") or {}
            report.update(head=head, atomic_reference_energies_eV=by_head.get(head))
        return report

    def _discovered_if_available(self) -> dict | None:
        """:meth:`discover`, or ``None`` when the model cannot be read here."""
        if not self.availability():
            return None
        try:
            discovered = self.discover()
        except Exception as exc:  # a corrupt or foreign checkpoint
            self.discovery_error = f"{type(exc).__name__}: {exc}"
            return None
        self.discovery_error = None
        return discovered

    def resolved_head(self) -> str | None:
        """The head MACE will evaluate (raises for an unknown or ambiguous head)."""
        discovered = self._discovered_if_available()
        if discovered is None:
            return self.spec.head
        return resolve_head(discovered, self.spec.head)

    def model_dtype(self) -> str | None:
        """The dtype the model's parameters were saved in, when readable."""
        discovered = self._discovered_if_available()
        return discovered.get("dtype") if discovered else None

    # -- cross-checks -----------------------------------------------------

    def inspect(self) -> dict:
        """Read elements, cutoff, heads and E0s from the model file when possible."""
        report = super().inspect()
        report.update(
            {
                "model_path": str(self.spec.model_path),
                "model_exists": self.spec.model_path.exists(),
                "device": self.spec.device,
                "precision": self.spec.precision,
                "head": self.spec.head,
                "declared_cutoff_angstrom": self.spec.cutoff_angstrom,
                "implementation_preference": list(self.spec.implementation),
            }
        )
        availability = self.availability()
        report["availability"] = availability.as_dict()
        if not availability:
            report["discovered"] = None
            report["discovery_note"] = (
                "model metadata not read: "
                + (availability.detail or f"missing {', '.join(availability.missing)}")
            )
            return report
        try:
            report["discovered"] = self.discover()
        except Exception as exc:  # a corrupt or unexpected checkpoint
            report["discovered"] = None
            report["discovery_note"] = f"could not read model metadata: {type(exc).__name__}: {exc}"
            return report
        report["conflicts"] = self._conflicts(report["discovered"])
        return report

    def _conflicts(self, discovered: dict | None) -> list[str]:
        """Every disagreement between the configuration and the model file."""
        if not discovered:
            return []
        conflicts: list[str] = []
        found = discovered.get("elements")
        if found:
            declared = set(self.spec.elements)
            extra = declared - set(found)
            if extra:
                conflicts.append(
                    f"potential.elements declares {', '.join(sorted(extra))}, which the "
                    f"model does not cover (model covers: {', '.join(found)})"
                )
            unlisted = set(found) - declared
            if unlisted:
                conflicts.append(
                    f"the model covers {', '.join(sorted(unlisted))}, which "
                    "potential.elements does not list; structures containing those "
                    "elements will be rejected by element-coverage validation"
                )
        cutoff = discovered.get("cutoff_angstrom")
        if cutoff and self.spec.cutoff_angstrom:
            if abs(cutoff - self.spec.cutoff_angstrom) > 1e-6:
                conflicts.append(
                    f"potential.cutoff_angstrom is {self.spec.cutoff_angstrom}, but the "
                    f"model reports {cutoff}"
                )
        if discovered.get("head_error"):
            conflicts.append(discovered["head_error"])
        conflicts += self._e0_conflicts(discovered.get("atomic_reference_energies_eV"))
        return conflicts

    def _e0_conflicts(self, model_e0s: dict | None) -> list[str]:
        declared = self.spec.atomic_reference_energies
        if not declared or not model_e0s:
            return []
        problems = []
        for symbol, value in sorted(declared.items()):
            if symbol not in model_e0s:
                continue
            if abs(model_e0s[symbol] - value) > E0_CONFLICT_TOLERANCE_EV:
                problems.append(
                    f"potential.atomic_reference_energies[{symbol!r}] = {value!r} eV, but the "
                    f"model's own E0 is {model_e0s[symbol]!r} eV (tolerance "
                    f"{E0_CONFLICT_TOLERANCE_EV:g} eV)"
                )
        return problems

    def model_problems(self) -> list[str]:
        """What makes this configuration unsafe to evaluate with this model file.

        A declared element the model does not cover, a cutoff, head or E0
        that disagrees with the model, and an unknown ``compile_mode``. The
        model-derived checks run only when the file can be read here (the
        route cannot run otherwise, and says so through availability). An
        element the model covers but the configuration does not list is
        harmless -- element coverage refuses such structures -- and is not
        a problem.
        """
        problems: list[str] = []
        if self.spec.compile_mode is not None and self.spec.compile_mode not in COMPILE_MODES:
            problems.append(
                f"potential.compile_mode = {self.spec.compile_mode!r} is not a torch.compile "
                f"mode; use one of {', '.join(COMPILE_MODES)} or leave it unset"
            )
        discovered = self._discovered_if_available()
        if discovered is None:
            if self.discovery_error:
                problems.append(
                    f"{self.spec.model_path.name} could not be read as a MACE model "
                    f"({self.discovery_error})"
                )
            return problems
        found = discovered.get("elements") or []
        if found:
            extra = set(self.spec.elements) - set(found)
            if extra:
                problems.append(
                    f"potential.elements declares {', '.join(sorted(extra))}, which "
                    f"{self.spec.model_path.name} does not cover (the model covers "
                    f"{', '.join(found)}); remove them from potential.elements"
                )
        cutoff = discovered.get("cutoff_angstrom")
        if cutoff and self.spec.cutoff_angstrom and abs(cutoff - self.spec.cutoff_angstrom) > 1e-6:
            problems.append(
                f"potential.cutoff_angstrom is {self.spec.cutoff_angstrom}, but "
                f"{self.spec.model_path.name} was trained with r_max = {cutoff}"
            )
        if discovered.get("head_error"):
            problems.append(discovered["head_error"])
        problems += self._e0_conflicts(discovered.get("atomic_reference_energies_eV"))
        return problems

    # -- atomic reference energies ----------------------------------------

    def model_atomic_reference_energies(self) -> dict[str, float] | None:
        """The model's own E0s for the resolved head, or ``None`` when unreadable.

        These -- not declared ones -- are what MACE adds into its total energy
        (``ScaleShiftMACE``: ``total = sum(E0) + interaction``), so they are
        the only ones a MACE total can be converted with.
        """
        discovered = self._discovered_if_available()
        if not discovered:
            return None
        e0s = discovered.get("atomic_reference_energies_eV")
        return dict(e0s) if e0s else None

    def atomic_reference_energies(self) -> dict[str, float] | None:
        """E0s: the model's own when readable, else the declared ones (unverified).

        :meth:`atomic_reference_energies_source` says which. A declared value
        that disagrees with the model is refused by :meth:`model_problems`.
        """
        model = self.model_atomic_reference_energies()
        if model:
            return model
        if self.spec.atomic_reference_energies:
            return dict(self.spec.atomic_reference_energies)
        return None

    def atomic_reference_energies_source(self) -> str | None:
        if self.model_atomic_reference_energies():
            return f"model file (head {self.resolved_head()!r})"
        if self.spec.atomic_reference_energies:
            return "potential.atomic_reference_energies (declared; not verified against the model)"
        return None


def resolve_head(discovered: dict, requested: str | None) -> str:
    """The head MACE's calculator will evaluate, or :class:`ConfigError`.

    Mirrors ``MACECalculator``'s rule (one head: that one; several: the one
    named ``default``, case-insensitively) -- except that an unknown requested
    head is refused, where MACE only logs a warning and silently evaluates
    the model's *last* head.
    """
    heads = list(discovered.get("heads") or ["Default"])
    if requested is not None:
        if requested not in heads:
            raise ConfigError(
                f"potential.head = {requested!r} is not a head of this model (heads: "
                f"{', '.join(heads)}); MACE would silently evaluate its last head, "
                f"{heads[-1]!r}, instead"
            )
        return requested
    if len(heads) == 1:
        return heads[0]
    default = [head for head in heads if head.lower() == "default"]
    if default:
        return default[0]
    raise ConfigError(
        f"this is a multi-head MACE model (heads: {', '.join(heads)}) with no head named "
        "'default'; set potential.head to the head to evaluate"
    )


def _describe_model(model) -> dict:
    """Plain-data facts about a loaded MACE model (nothing torch-typed)."""
    elements = _model_elements(model)
    heads = getattr(model, "heads", None)
    heads = [str(h) for h in heads] if heads is not None else ["Default"]
    parameters = next(iter(model.parameters()), None)
    try:
        from mace.modules import ScaleShiftMACE

        scale_shift = isinstance(model, ScaleShiftMACE)
    except ImportError:  # pragma: no cover - mace is what loaded the model
        scale_shift = None
    return {
        "elements": elements,
        "cutoff_angstrom": _model_cutoff(model),
        "heads": heads,
        "atomic_reference_energies_by_head": _model_atomic_energies_by_head(
            model, elements, heads
        ),
        "model_class": type(model).__name__,
        "scale_shift": scale_shift,
        "dtype": str(parameters.dtype).replace("torch.", "") if parameters is not None else None,
    }


def _copy(value):
    if isinstance(value, dict):
        return {k: _copy(v) for k, v in value.items()}
    if isinstance(value, list):
        return list(value)
    return value


def _model_elements(model) -> list[str]:
    """Best-effort element table from a loaded MACE model.

    MACE checkpoints have carried this as ``atomic_numbers`` on the model or
    inside ``z_table`` depending on version, so several shapes are tried and
    an unrecognised one yields an empty list rather than a wrong answer.
    """
    from ase.data import chemical_symbols

    numbers = None
    for attribute in ("atomic_numbers", "z_table"):
        candidate = getattr(model, attribute, None)
        if candidate is None:
            continue
        zs = getattr(candidate, "zs", candidate)
        try:
            numbers = [int(z) for z in zs]
            break
        except TypeError:
            continue
    if not numbers:
        return []
    return [chemical_symbols[z] for z in numbers]


def _model_cutoff(model) -> float | None:
    for attribute in ("r_max", "cutoff"):
        value = getattr(model, attribute, None)
        if value is None:
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            try:
                return float(value.item())
            except Exception:
                continue
    return None


def _e0_rows(model) -> list[list[float]] | None:
    """The ``AtomicEnergiesBlock`` buffer as rows (one per head), or ``None``.

    MACE stores ``[n_elements]`` for a single head and ``[n_heads,
    n_elements]`` for several (``mace/modules/blocks.py``,
    ``AtomicEnergiesBlock``: ``atleast_2d(atomic_energies).T`` is multiplied
    by the one-hot element matrix, giving one column per head).
    """
    block = getattr(model, "atomic_energies_fn", None)
    values = getattr(block, "atomic_energies", None) if block is not None else None
    if values is None:
        return None
    try:
        array = values.detach().cpu()
        rows = array.tolist() if array.dim() == 2 else [array.flatten().tolist()]
    except Exception:
        return None
    return [[float(v) for v in row] for row in rows]


def _model_atomic_energies_by_head(model, elements: list[str], heads: list[str]):
    """``{head: {symbol: E0}}``, or ``None`` for an unrecognised layout."""
    rows = _e0_rows(model)
    if rows is None or not elements or len(rows) != len(heads):
        return None
    if any(len(row) != len(elements) for row in rows):
        return None
    return {head: dict(zip(elements, row)) for head, row in zip(heads, rows)}


def _model_atomic_energies(model, elements: list[str]) -> dict[str, float] | None:
    """Per-element E0s of a single-head model, as a symbol-keyed mapping in eV.

    ``None`` for a multi-head model (use :func:`_model_atomic_energies_by_head`)
    or an unrecognised layout: guessing here would corrupt every
    energy-convention conversion downstream.
    """
    rows = _e0_rows(model)
    if rows is None or len(rows) != 1 or not elements or len(rows[0]) != len(elements):
        return None
    return dict(zip(elements, rows[0]))


def build_adapter(spec: MacePotentialSpec) -> MaceAdapter:
    if not isinstance(spec, MacePotentialSpec):
        raise ConfigError(f"expected a MACE potential spec; got {type(spec).__name__}")
    return MaceAdapter(spec)


__all__ = [
    "COMPILE_MODES",
    "E0_CONFLICT_TOLERANCE_EV",
    "MaceAdapter",
    "build_adapter",
    "clear_discovery_cache",
    "model_sha256",
    "resolve_head",
]
