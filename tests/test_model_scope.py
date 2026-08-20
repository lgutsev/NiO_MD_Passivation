import pytest

from nio_md_prep.analysis.model_scope import (
    assessment_from_summary,
    model_assessment,
)


def test_classical_profile_withholds_model_sensitive_observables():
    assessment = model_assessment("classical-ff")
    assert assessment["kinetics_reportable"] is False
    assert assessment["electronic_proxy_reportable"] is False
    assert assessment["withheld_from_workbooks"]


def test_mlip_profile_enables_kinetics_but_not_dipole_proxy():
    assessment = model_assessment("mlip")
    assert assessment["kinetics_reportable"] is True
    assert assessment["electronic_proxy_reportable"] is False


def test_charge_aware_profile_enables_both_sensitive_families():
    assessment = model_assessment("charge-aware-mlip")
    assert assessment["kinetics_reportable"] is True
    assert assessment["electronic_proxy_reportable"] is True


def test_legacy_summary_is_treated_conservatively():
    assessment = assessment_from_summary({})
    assert assessment["profile"] == "classical-ff"


def test_unknown_profile_is_rejected():
    with pytest.raises(ValueError, match="analysis profile"):
        model_assessment("unvalidated-magic")
