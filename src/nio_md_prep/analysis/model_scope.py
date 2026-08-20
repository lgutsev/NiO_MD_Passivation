"""Force-field-aware evidence scopes for trajectory-derived observables."""
from __future__ import annotations

from typing import Any


MODEL_PROFILES = (
    "classical-ff",
    "mlip",
    "charge-aware-mlip",
)


def model_assessment(profile: str) -> dict[str, Any]:
    """Return the reporting policy associated with a trajectory model.

    All observables remain in the detailed JSON audit trail.  The booleans
    below control whether model-sensitive values are promoted into workbook
    tables intended for interpretation or manuscript drafting.
    """
    if profile not in MODEL_PROFILES:
        choices = ", ".join(MODEL_PROFILES)
        raise ValueError(f"analysis profile must be one of: {choices}")

    kinetics_reportable = profile in {"mlip", "charge-aware-mlip"}
    electronic_proxy_reportable = profile == "charge-aware-mlip"
    withheld = []
    if not kinetics_reportable:
        withheld.append(
            "contact persistence, residence times, site-exchange counts, "
            "and exchange rates"
        )
    if not electronic_proxy_reportable:
        withheld.append("z-dipole and approximate potential-step proxy")

    if profile == "classical-ff":
        note = (
            "Classical fixed-charge force field: emphasize coverage, void "
            "topology, roughness, adsorption heights, broad orientation, "
            "and qualitative anchoring. Kinetic and electronic-adjacent "
            "observables are retained only in JSON and withheld from "
            "workbook-facing values."
        )
    elif profile == "mlip":
        note = (
            "MLIP trajectory: kinetic observables are eligible for reporting "
            "only after validation against adsorption geometries, relative "
            "binding energies, forces, and relevant barrier configurations. "
            "A standard energy/force MLIP does not validate dipoles, charge "
            "transfer, or work-function changes."
        )
    else:
        note = (
            "Charge-aware MLIP trajectory: kinetic and dipole-proxy values "
            "are eligible for reporting only after separate validation of "
            "the potential, predicted charges/dipoles, and relevant barrier "
            "configurations. The potential-step remains an idealized proxy, "
            "not a direct work-function prediction."
        )

    return {
        "profile": profile,
        "kinetics_reportable": kinetics_reportable,
        "electronic_proxy_reportable": electronic_proxy_reportable,
        "primary_structural_observables": [
            "near-surface coverage and uncovered fraction",
            "void topology",
            "top-of-film roughness and height distributions",
            "adsorption heights and broad tilt/orientational distributions",
        ],
        "qualitative_observables": [
            "exposed-Ni site ownership",
            "bound and multi-anchor fractions",
            "component RDFs and density profiles",
            "dangling-anchor density and decompression-induced changes",
        ],
        "withheld_from_workbooks": withheld,
        "interpretation_note": note,
    }


def assessment_from_summary(summary: dict[str, Any]) -> dict[str, Any]:
    """Read an assessment, treating legacy summaries conservatively."""
    assessment = summary.get("model_assessment")
    if (
        isinstance(assessment, dict)
        and assessment.get("profile") in MODEL_PROFILES
    ):
        return model_assessment(str(assessment["profile"]))
    return model_assessment("classical-ff")
