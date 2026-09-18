"""The MACE potential family, described independently of ASE/LAMMPS/OpenMM.

Everything torch-flavoured is behind a lazy import. :meth:`MaceAdapter.inspect`
is the only method that opens the model file, and it is written so that a
machine without torch still gets a useful report -- the declared metadata,
plus an explicit statement that nothing was read.

A trained MACE model carries its own element table, cutoff and per-element
atomic reference energies (E0). Those are exactly the three things this
subsystem needs and cannot safely guess, so when they *are* discoverable they
are cross-checked against what the configuration declares, and a disagreement
is reported rather than silently resolved in either direction.
"""
from __future__ import annotations

from ..capabilities import CapabilitySet
from ..environment import Availability, probe
from ..errors import ConfigError, MissingDependencyError
from ..specs import MacePotentialSpec
from .base import PotentialAdapter


class MaceAdapter(PotentialAdapter):
    kind = "mace"
    requires = ("mace", "torch")

    spec: MacePotentialSpec

    def capabilities(self) -> CapabilitySet:
        """What a MACE model can produce, before any engine constrains it.

        MACE is a gradient-domain model: energy, forces, stress and a per-atom
        (site) energy decomposition all come out of the same forward pass, and
        periodic cells are supported natively. GPU support is a property of
        the model's device, which is why it is read from the spec rather than
        assumed.
        """
        return CapabilitySet(
            energy=True,
            forces=True,
            stress=True,
            per_atom_energy=True,
            periodic=True,
            gpu=self.spec.wants_gpu,
            elements=frozenset(self.spec.elements),
            precisions=frozenset({self.spec.precision}),
            engines=frozenset({"ase", "lammps", "openmm"}),
            native_energy_convention=self.spec.energy_convention,
            convertible_energy_conventions=(
                frozenset({"total", "interaction"})
                if self.spec.atomic_reference_energies
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

    def load_model(self):
        """Load the torch model. The first genuinely heavy call in this package."""
        if not self.spec.model_path.exists():
            raise FileNotFoundError(f"MACE model not found: {self.spec.model_path}")
        try:
            import torch
        except ImportError as exc:
            raise MissingDependencyError(
                f"loading the MACE model {self.spec.model_path.name}",
                ("torch",),
                hint="Install the MACE extra: pip install 'nio-md-prep[mace]'",
            ) from exc
        return torch.load(
            self.spec.model_path, map_location=self._torch_device(), weights_only=False
        )

    def _torch_device(self):
        return self.spec.device

    def torch_dtype(self):
        import torch

        return {"float32": torch.float32, "float64": torch.float64}[self.spec.precision]

    def inspect(self) -> dict:
        """Read elements, cutoff and E0s from the model file when possible."""
        report = super().inspect()
        report.update(
            {
                "model_path": str(self.spec.model_path),
                "model_exists": self.spec.model_path.exists(),
                "device": self.spec.device,
                "precision": self.spec.precision,
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

    def discover(self) -> dict:
        """Pull element table, cutoff and atomic reference energies from the model."""
        model = self.load_model()
        elements = _model_elements(model)
        cutoff = _model_cutoff(model)
        e0s = _model_atomic_energies(model, elements)
        return {
            "elements": elements,
            "cutoff_angstrom": cutoff,
            "atomic_reference_energies_eV": e0s,
            "model_class": type(model).__name__,
        }

    def _conflicts(self, discovered: dict | None) -> list[str]:
        """Report -- never silently resolve -- disagreements with the config."""
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
        return conflicts

    def atomic_reference_energies(self) -> dict[str, float] | None:
        """E0s from the configuration, falling back to the model file.

        These are what make an ASE/OpenMM energy comparison possible at all,
        so the configuration wins when it declares them: a committee member
        re-fitted against different references must be describable.
        """
        if self.spec.atomic_reference_energies:
            return dict(self.spec.atomic_reference_energies)
        if not self.availability():
            return None
        try:
            discovered = self.discover()
        except Exception:
            return None
        return discovered.get("atomic_reference_energies_eV") or None


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


def _model_atomic_energies(model, elements: list[str]) -> dict[str, float] | None:
    """Per-element E0s, as a symbol-keyed mapping in eV.

    MACE stores these in an ``AtomicEnergiesBlock``; the layout has moved
    between versions, so an unrecognised one returns ``None`` and the caller
    falls back to whatever the configuration declares. Guessing here would
    corrupt every energy-convention conversion downstream.
    """
    block = getattr(model, "atomic_energies_fn", None)
    values = getattr(block, "atomic_energies", None) if block is not None else None
    if values is None:
        return None
    try:
        flat = [float(v) for v in values.detach().cpu().flatten().tolist()]
    except Exception:
        try:
            flat = [float(v) for v in values]
        except Exception:
            return None
    if not elements or len(flat) != len(elements):
        return None
    return dict(zip(elements, flat))


def build_adapter(spec: MacePotentialSpec) -> MaceAdapter:
    if not isinstance(spec, MacePotentialSpec):
        raise ConfigError(f"expected a MACE potential spec; got {type(spec).__name__}")
    return MaceAdapter(spec)


__all__ = ["MaceAdapter", "build_adapter"]
