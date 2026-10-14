"""The MLIP manifest: what every job directory must record.

Provenance is not optional decoration here. Once several MACE committee
members exist -- and later DeepMD or DPA models -- a trajectory is only
interpretable if its directory says exactly which model produced it, under
which conventions, with which rendered engine commands.
"""
import json
from pathlib import Path

import pytest

from nio_md_prep.mlip import provenance
from nio_md_prep.mlip.errors import ModelIntegrityError
from nio_md_prep.mlip.registry import resolve_bridge
from nio_md_prep.mlip.specs import (
    EngineSpec,
    JobSpec,
    LammpsMlipPotentialSpec,
    MacePotentialSpec,
    SimulationSpec,
)

#: Fields a reader must be able to find in any MLIP manifest.
REQUIRED_TOP_LEVEL = (
    "manifest_version",
    "created_utc",
    "git",
    "potential",
    "engine",
    "bridge",
    "simulation",
    "units",
    "energy_convention",
    "element_mapping",
    "environment",
    "engine_parameters",
)


def make_model(tmp_path: Path, name: str = "model.pt", body: bytes = b"weights") -> Path:
    path = tmp_path / name
    path.write_bytes(body)
    return path


def mace_job(tmp_path: Path, **overrides) -> JobSpec:
    model = make_model(tmp_path)
    potential = MacePotentialSpec(
        label="nio_phosphonate",
        model_path=model,
        model_sha256=provenance.sha256_file(model),
        declared_elements=("Ni", "O", "P", "C", "H"),
        device="cuda",
        precision="float32",
        atomic_reference_energies={"Ni": -5.78, "O": -4.95, "P": -5.41, "C": -8.02, "H": -1.12},
        **overrides,
    )
    return JobSpec(
        potential=potential,
        engine=EngineSpec(kind="ase"),
        simulation=SimulationSpec(
            task="md",
            ensemble="nvt",
            temperature_K=400,
            timestep_fs=0.5,
            steps=10000,
            seed=12345,
            thermostat="langevin",
        ),
        name="mace-nvt",
    )


def test_manifest_records_everything_a_rerun_needs(tmp_path):
    job = mace_job(tmp_path)
    registration = resolve_bridge("mace", "ase")
    manifest = provenance.build_manifest(
        job=job,
        registration=registration,
        engine_parameters={"native_units": "ase", "energy_convention": "total"},
    )
    for key in REQUIRED_TOP_LEVEL:
        assert key in manifest, key

    assert manifest["potential"]["model_sha256"] == job.potential.model_sha256
    assert manifest["potential"]["model_sha256_observed"]
    assert manifest["potential"]["precision"] == "float32"
    assert manifest["potential"]["device"] == "cuda"
    assert manifest["bridge"]["implementation"] == "mace-ase-calculator"
    assert manifest["simulation"]["integrator"]["thermostat"] == "langevin"
    assert manifest["simulation"]["integrator"]["seed"] == 12345
    assert manifest["element_mapping"]["elements"] == ["Ni", "O", "P", "C", "H"]
    assert manifest["element_mapping"]["atomic_reference_energies_eV"]["Ni"] == -5.78
    assert manifest["units"]["canonical_energy"] == "eV"
    assert manifest["energy_convention"]["potential_native"] == "total"
    assert "packages" in manifest["environment"]
    assert "device" in manifest["environment"]


def test_manifest_preserves_lammps_commands_verbatim(tmp_path):
    """Reproducing a LAMMPS run means reproducing these exact strings."""
    model = make_model(tmp_path, "nio.pb")
    potential = LammpsMlipPotentialSpec(
        label="nio-deepmd",
        pair_style="deepmd nio.pb",
        pair_coeff=("* *",),
        type_map={1: "Ni", 2: "O"},
        model_paths=(model,),
        framework="deepmd",
    )
    job = JobSpec(
        potential=potential,
        engine=EngineSpec(kind="lammps"),
        simulation=SimulationSpec(task="singlepoint"),
    )
    manifest = provenance.build_manifest(
        job=job,
        registration=resolve_bridge("lammps", "lammps"),
        engine_parameters={"pair_commands": list(potential.render_pair_commands())},
    )
    assert manifest["potential"]["rendered_commands"] == [
        "pair_style deepmd nio.pb",
        "pair_coeff * *",
    ]
    assert manifest["engine_parameters"]["pair_commands"] == [
        "pair_style deepmd nio.pb",
        "pair_coeff * *",
    ]
    assert manifest["element_mapping"]["lammps_type_map"] == {"1": "Ni", "2": "O"}


def test_a_changed_model_file_stops_the_job(tmp_path):
    """Provenance that lies is worse than none: refuse rather than record it."""
    job = mace_job(tmp_path)
    job.potential.model_path.write_bytes(b"different weights")
    with pytest.raises(ModelIntegrityError, match="Refusing to run"):
        provenance.build_manifest(
            job=job, registration=resolve_bridge("mace", "ase")
        )


def test_hash_verification_can_be_deferred_for_reporting(tmp_path):
    job = mace_job(tmp_path)
    job.potential.model_path.write_bytes(b"different weights")
    manifest = provenance.build_manifest(
        job=job, registration=resolve_bridge("mace", "ase"), verify=False
    )
    observed = manifest["potential"]["model_sha256_observed"]
    assert list(observed.values())[0] != job.potential.model_sha256


def test_a_missing_model_file_is_reported_not_hidden(tmp_path):
    potential = MacePotentialSpec(
        model_path=tmp_path / "absent.pt", declared_elements=("Ni",)
    )
    job = JobSpec(
        potential=potential,
        engine=EngineSpec(kind="ase"),
        simulation=SimulationSpec(),
    )
    manifest = provenance.build_manifest(job=job, registration=resolve_bridge("mace", "ase"))
    assert manifest["potential"]["model_files_missing"] == [str(tmp_path / "absent.pt")]


def test_manifest_round_trips_through_disk(tmp_path):
    job = mace_job(tmp_path)
    manifest = provenance.build_manifest(job=job, registration=resolve_bridge("mace", "ase"))
    path = provenance.write_manifest(tmp_path / "run", manifest)
    assert path.name == provenance.MANIFEST_NAME
    assert provenance.read_manifest(tmp_path / "run")["name"] == "mace-nvt"
    assert json.loads(path.read_text(encoding="utf-8"))["manifest_version"] == 1


def test_sha256_matches_hashlib(tmp_path):
    import hashlib

    path = make_model(tmp_path, "blob.bin", b"x" * 5000)
    assert provenance.sha256_file(path) == hashlib.sha256(b"x" * 5000).hexdigest()


def test_device_report_does_not_import_torch_when_absent():
    import sys

    report = provenance.device_report("cuda")
    assert report["requested"] == "cuda"
    assert "platform" in report and "cpu_count" in report
    if "torch" not in sys.modules:
        assert report["torch"] is None


def test_package_versions_report_absent_packages_as_null():
    versions = provenance.package_versions(("nio-md-prep", "definitely-not-installed"))
    assert versions["nio-md-prep"]
    assert versions["definitely-not-installed"] is None
