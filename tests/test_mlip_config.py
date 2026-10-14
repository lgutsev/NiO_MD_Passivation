"""Configuration parsing and spec validation, with no execution anywhere.

Parsing must be usable on a machine that cannot run the job, so these tests
also pin that importing the config layer pulls in no backend.
"""
import sys
import textwrap
from pathlib import Path

import pytest

from nio_md_prep.mlip.config import parse_job, parse_job_mapping
from nio_md_prep.mlip.errors import ConfigError, UnknownPotentialError
from nio_md_prep.mlip.specs import (
    EngineSpec,
    LammpsMlipPotentialSpec,
    MacePotentialSpec,
    SimulationSpec,
)

ROOT = Path(__file__).parents[1]

MACE_CONFIG = """
[potential]
kind = "mace"
model = "nio_phosphonate.model"
elements = ["Ni", "O", "P", "C", "H"]
device = "cuda"
precision = "float32"

[engine]
kind = "ase"

[simulation]
task = "md"
ensemble = "nvt"
temperature_K = 400
timestep_fs = 0.5
steps = 10000
seed = 12345
"""

LAMMPS_CONFIG = """
[potential]
kind = "lammps"
pair_style = "mliap unified mace_nio.pt 0"
pair_coeff = ["* * Ni O P C H"]
units = "metal"
model_files = ["mace_nio.pt"]

[potential.type_map]
1 = "Ni"
2 = "O"
3 = "P"
4 = "C"
5 = "H"

[engine]
kind = "lammps"

[simulation]
task = "singlepoint"
"""


def write(tmp_path: Path, text: str, name: str = "job.toml") -> Path:
    path = tmp_path / name
    path.write_text(textwrap.dedent(text), encoding="utf-8")
    return path


def test_parsing_imports_no_backend(tmp_path):
    parse_job(write(tmp_path, MACE_CONFIG))
    assert "torch" not in sys.modules
    assert "openmm" not in sys.modules


def test_mace_config_round_trips_the_documented_example(tmp_path):
    job = parse_job(write(tmp_path, MACE_CONFIG))
    assert isinstance(job.potential, MacePotentialSpec)
    assert job.potential.device == "cuda"
    assert job.potential.precision == "float32"
    assert job.potential.elements == ("Ni", "O", "P", "C", "H")
    assert job.engine == EngineSpec(kind="ase")
    assert job.simulation.ensemble == "nvt"
    assert job.simulation.temperature_K == 400
    assert job.simulation.steps == 10000
    assert job.simulation.seed == 12345


def test_model_paths_resolve_against_the_config_directory(tmp_path):
    job = parse_job(write(tmp_path, MACE_CONFIG))
    assert job.potential.model_path.is_absolute()
    assert job.potential.model_path.name == "nio_phosphonate.model"


def test_lammps_config_keeps_the_pair_style_verbatim(tmp_path):
    job = parse_job(write(tmp_path, LAMMPS_CONFIG))
    potential = job.potential
    assert isinstance(potential, LammpsMlipPotentialSpec)
    assert potential.pair_style == "mliap unified mace_nio.pt 0"
    assert potential.type_map == {1: "Ni", 2: "O", 3: "P", 4: "C", 5: "H"}
    assert potential.elements == ("Ni", "O", "P", "C", "H")
    assert potential.render_pair_commands() == (
        "pair_style mliap unified mace_nio.pt 0",
        "pair_coeff * * Ni O P C H",
    )


@pytest.mark.parametrize(
    "pair_style,framework",
    [
        ("deepmd model.pb", "deepmd"),
        ("pace", "pace"),
        ("mliap unified mace.pt 0", "mliap"),
        ("mace no_domain_decomposition", "mace"),
        ("some_future_style model.bin", "some_future_style"),
    ],
)
def test_lammps_mlip_is_not_synonymous_with_one_framework(pair_style, framework):
    """Every one of these must be expressible without a code change."""
    from nio_md_prep.mlip.potentials.lammps_mlip import LammpsMlipAdapter

    spec = LammpsMlipPotentialSpec(
        pair_style=pair_style, pair_coeff=("* *",), type_map={1: "Ni"}
    )
    assert LammpsMlipAdapter(spec).framework == framework


def test_unknown_keys_are_errors_not_warnings(tmp_path):
    """A typo must not silently become a 0 K trajectory."""
    config = MACE_CONFIG.replace("temperature_K = 400", "temperature = 400")
    with pytest.raises(ConfigError, match="unknown key"):
        parse_job(write(tmp_path, config))


def test_unknown_potential_kind_names_the_known_ones(tmp_path):
    config = MACE_CONFIG.replace('kind = "mace"', 'kind = "dpa3"')
    with pytest.raises(UnknownPotentialError, match="known kinds"):
        parse_job(write(tmp_path, config))


def test_mace_must_declare_its_elements(tmp_path):
    config = MACE_CONFIG.replace('elements = ["Ni", "O", "P", "C", "H"]\n', "")
    with pytest.raises(ConfigError, match="must declare 'elements'"):
        parse_job(write(tmp_path, config))


def test_lammps_mlip_must_declare_a_type_map():
    with pytest.raises(ConfigError, match="type_map is required"):
        LammpsMlipPotentialSpec(pair_style="deepmd m.pb", pair_coeff=("* *",))


def test_type_map_must_have_no_gaps():
    with pytest.raises(ConfigError, match="1..N with no gaps"):
        LammpsMlipPotentialSpec(
            pair_style="deepmd m.pb", pair_coeff=("* *",), type_map={1: "Ni", 3: "O"}
        )


def test_lj_units_are_refused_at_spec_construction():
    with pytest.raises(Exception, match="no canonical"):
        LammpsMlipPotentialSpec(
            pair_style="x", pair_coeff=("* *",), type_map={1: "Ni"}, units="lj"
        )


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("device", "tpu", "must start with one of"),
        ("precision", "float16", "must be one of"),
        ("model_sha256", "not-a-hash", "hex SHA256"),
    ],
)
def test_mace_spec_rejects_impossible_values(field, value, message):
    kwargs = {"model_path": Path("m.pt"), "declared_elements": ("Ni",), field: value}
    with pytest.raises(ConfigError, match=message):
        MacePotentialSpec(**kwargs)


def test_md_requires_an_ensemble_and_a_timestep():
    with pytest.raises(ConfigError, match="ensemble is required"):
        SimulationSpec(task="md")
    with pytest.raises(ConfigError, match="timestep_fs"):
        SimulationSpec(task="md", ensemble="nve", steps=10)


def test_nvt_requires_a_temperature():
    with pytest.raises(ConfigError, match="temperature_K is required"):
        SimulationSpec(task="md", ensemble="nvt", timestep_fs=0.5, steps=10)


def test_npt_requires_a_pressure():
    with pytest.raises(ConfigError, match="pressure_bar is required"):
        SimulationSpec(
            task="md", ensemble="npt", timestep_fs=0.5, steps=10, temperature_K=300
        )


def test_ensemble_is_meaningless_without_md():
    with pytest.raises(ConfigError, match="only meaningful for task"):
        SimulationSpec(task="singlepoint", ensemble="nvt")


def test_region_selection_is_reserved_but_refused():
    """The ML/MM hook exists in the schema and is explicitly not implemented."""
    with pytest.raises(ConfigError, match="reserved but not implemented"):
        SimulationSpec(region="selection", selection="index 1-10")


def test_region_all_is_the_implemented_meaning():
    assert SimulationSpec().region == "all"


def test_every_shipped_example_parses():
    examples = sorted((ROOT / "examples" / "mlip").glob("*.toml"))
    assert examples, "the documented examples must exist"
    for path in examples:
        job = parse_job(path)
        assert job.potential.kind and job.engine.kind


def test_a_mapping_can_be_parsed_without_a_file():
    job = parse_job_mapping(
        {
            "potential": {"kind": "mock", "elements": ["Ni"]},
            "engine": {"kind": "ase"},
            "simulation": {"task": "singlepoint"},
        }
    )
    assert job.potential.kind == "mock"
