"""Configuration parsing and spec validation, with no execution anywhere.

Parsing must be usable on a machine that cannot run the job, so these tests
also pin that importing the config layer pulls in no backend.
"""
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


def test_parsing_imports_no_backend(tmp_path, backend_import_guard):
    parse_job(write(tmp_path, MACE_CONFIG))
    backend_import_guard()


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


# --- fields added for the engine routes -----------------------------------


def test_engine_fields_parse_from_toml(tmp_path):
    config = LAMMPS_CONFIG.replace(
        '[engine]\nkind = "lammps"',
        '[engine]\nkind = "lammps"\nruntime = "executable"\nexecutable = "lmp_mpi"\n'
        'mpi_launcher = ["mpirun", "-np", "4"]\nlammps_args = ["-sf", "omp"]\n'
        "timeout_s = 600\nthreads = 2",
    )
    engine = parse_job(write(tmp_path, config)).engine
    assert engine.runtime == "executable"
    assert engine.mpi_launcher == ("mpirun", "-np", "4")
    assert engine.lammps_args == ("-sf", "omp")
    assert engine.timeout_s == 600 and engine.threads == 2
    assert engine.as_dict()["mpi_launcher"] == ["mpirun", "-np", "4"]


def test_an_argument_vector_must_not_be_one_string():
    with pytest.raises(ConfigError, match="list of arguments"):
        EngineSpec(kind="lammps", mpi_launcher="mpirun -np 4")


@pytest.mark.parametrize(
    "kwargs,message",
    [
        ({"kind": "ase", "executable": "lmp"}, "only apply to the LAMMPS engine"),
        ({"kind": "openmm", "runtime": "executable"}, "only apply to the LAMMPS engine"),
        ({"kind": "ase", "platform": "CUDA"}, "only apply to the OpenMM engine"),
        ({"kind": "lammps", "runtime": "python", "executable": "lmp"}, "would be ignored"),
        ({"kind": "lammps", "runtime": "python", "timeout_s": 10}, "would be ignored"),
        ({"kind": "lammps", "runtime": "gpu"}, "engine.runtime must be one of"),
        ({"kind": "openmm", "platform": "Cuda"}, "engine.platform must be one of"),
        ({"kind": "openmm", "platform": "CPU", "platform_precision": "mixed"}, "only the CUDA"),
        ({"kind": "openmm", "platform_precision": "double"}, "only the CUDA"),
        ({"kind": "openmm", "platform": "CUDA", "platform_precision": "half"}, "must be one of"),
        ({"kind": "ase", "threads": 0}, "positive integer"),
        ({"kind": "lammps", "timeout_s": -1}, "positive number"),
    ],
)
def test_engine_fields_are_refused_where_they_would_be_ignored(kwargs, message):
    with pytest.raises(ConfigError, match=message):
        EngineSpec(**kwargs)


def test_platform_precision_is_accepted_on_a_gpu_platform():
    engine = EngineSpec(kind="openmm", platform="CUDA", platform_precision="mixed")
    assert engine.as_dict()["platform_precision"] == "mixed"


def test_retargeting_keeps_only_what_the_target_engine_reads():
    openmm = EngineSpec(kind="openmm", platform="CUDA", platform_precision="double", threads=4)
    ase = openmm.retargeted("ase")
    assert ase.kind == "ase" and ase.platform is None and ase.platform_precision is None
    assert ase.threads == 4
    lammps = EngineSpec(kind="lammps", executable="lmp", mpi_launcher=("mpirun",))
    assert lammps.retargeted("openmm").executable is None


def test_mace_head_parses(tmp_path):
    config = MACE_CONFIG.replace('precision = "float32"', 'precision = "float32"\nhead = "pbe_u"')
    potential = parse_job(write(tmp_path, config)).potential
    assert potential.head == "pbe_u"
    assert potential.as_dict()["head"] == "pbe_u"


def test_structure_pbc_parses_and_is_required_for_lammps_data(tmp_path):
    base = LAMMPS_CONFIG + '\n[structure]\npath = "slab.lmp"\nformat = "lammps-data-nio"\n'
    with pytest.raises(ConfigError, match="needs an explicit structure.pbc"):
        parse_job(write(tmp_path, base))
    job = parse_job(write(tmp_path, base + "pbc = [true, true, false]\n"))
    assert job.structure.pbc == (True, True, False)
    assert job.as_dict()["structure"]["pbc"] == [True, True, False]
    with pytest.raises(ConfigError, match="three booleans"):
        parse_job(write(tmp_path, base + "pbc = [1, 1, 0]\n"))


def test_an_unknown_structure_format_is_refused_at_parse_time(tmp_path):
    config = MACE_CONFIG + '\n[structure]\npath = "x.gro"\nformat = "gromacs"\n'
    with pytest.raises(ConfigError, match="structure.format"):
        parse_job(write(tmp_path, config))


# --- simulation fields that only mean something in some ensembles ---------


def md(**kwargs):
    base = {"task": "md", "ensemble": "nvt", "temperature_K": 300.0, "timestep_fs": 0.5,
            "steps": 10}
    base.update(kwargs)
    return SimulationSpec(**base)


@pytest.mark.parametrize(
    "kwargs,message",
    [
        ({"ensemble": "nve", "thermostat": "langevin"}, "cannot apply to the nve"),
        ({"ensemble": "nve", "thermostat_damping_fs": 50.0}, "cannot apply to the nve"),
        ({"barostat_coupling": "isotropic"}, "only meaningful for the npt"),
        ({"barostat_damping_fs": 500.0}, "only meaningful for the npt"),
        ({"ensemble": "npt", "pressure_bar": 1.0, "barostat_coupling": "xy"}, "must be one of"),
        ({"thermostat_damping_fs": 0.0}, "positive number"),
        ({"seed": 0}, "1..2147483647"),
        ({"seed": 2**31}, "1..2147483647"),
        ({"seed": True}, "1..2147483647"),
        ({"trajectory_interval": 0}, "positive integer"),
        ({"log_interval": 2.5}, "positive integer"),
        ({"compute_stress": "yes"}, "true or false"),
        ({"vacuum_gap_threshold_angstrom": -1.0}, "positive number"),
    ],
)
def test_md_settings_are_refused_where_they_would_be_ignored(kwargs, message):
    with pytest.raises(ConfigError, match=message):
        md(**kwargs)


def test_thermostat_settings_are_refused_for_a_single_point():
    with pytest.raises(ConfigError, match="only apply to task = 'md'"):
        SimulationSpec(task="singlepoint", thermostat="langevin")


def test_the_single_point_form_drops_every_md_only_setting():
    from nio_md_prep.mlip.specs import as_singlepoint

    npt = md(ensemble="npt", pressure_bar=1.0, thermostat="langevin", barostat="mtk",
             thermostat_damping_fs=50.0, barostat_damping_fs=500.0, barostat_coupling="isotropic")
    endpoint = as_singlepoint(npt)
    assert endpoint.task == "singlepoint" and endpoint.compute_stress
    assert endpoint.barostat_coupling is None and endpoint.thermostat_damping_fs is None
    in_plane = as_singlepoint(md(ensemble="npt", pressure_bar=1.0, barostat_coupling="in-plane"))
    assert not in_plane.compute_stress  # a slab has no full stress tensor


def test_one_default_damping_time_for_every_engine():
    """Previously ASE used tau = 100 fs and OpenMM 1/ps (tau = 1000 fs) by default."""
    from nio_md_prep.mlip import specs
    from nio_md_prep.mlip.engines import ase_engine, openmm_engine

    assert specs.DEFAULT_THERMOSTAT_DAMPING_FS == 100.0
    assert md().resolved_thermostat_damping_fs == 100.0
    assert md(thermostat_damping_fs=25.0).resolved_thermostat_damping_fs == 25.0
    assert ase_engine.DEFAULT_THERMOSTAT_DAMPING_FS is specs.DEFAULT_THERMOSTAT_DAMPING_FS
    # OpenMM's Langevin friction is 1/tau: tau = 100 fs is 10 per ps.
    assert openmm_engine.DEFAULT_FRICTION_PER_PS == pytest.approx(10.0)
    # LAMMPS reads simulation.resolved_thermostat_damping_fs; its rendering in
    # metal/real time units is pinned in test_mlip_lammps_deck.py.


def test_a_drawn_seed_is_valid_everywhere():
    from nio_md_prep.mlip.specs import MAX_SEED, resolve_seed

    assert resolve_seed(12345) == 12345
    drawn = resolve_seed(None)
    assert 1 <= drawn <= MAX_SEED


# --- potential fields ------------------------------------------------------


def test_a_duplicate_element_type_map_is_refused_with_the_reason():
    """Two Ni types (AFM sublattices) cannot be told apart from an ASE structure."""
    with pytest.raises(ConfigError, match="magnetic") as excinfo:
        LammpsMlipPotentialSpec(
            pair_style="deepmd m.pb", pair_coeff=("* *",), type_map={1: "Ni", 2: "Ni", 3: "O"}
        )
    assert "Ni: types [1, 2]" in str(excinfo.value)


def test_a_non_integer_type_map_key_is_a_config_error(tmp_path):
    config = LAMMPS_CONFIG.replace('1 = "Ni"', 'Ni = "Ni"')
    with pytest.raises(ConfigError, match="LAMMPS type numbers"):
        parse_job(write(tmp_path, config))


def test_model_hash_keys_resolve_like_model_files(tmp_path):
    """The shipped lammps-native.toml pattern: a relative file with a correct hash."""
    from nio_md_prep.mlip import provenance
    from nio_md_prep.mlip.errors import ModelIntegrityError
    from nio_md_prep.mlip.registry import resolve_bridge

    model = tmp_path / "nio.pb"
    model.write_bytes(b"frozen graph")
    digest = provenance.sha256_file(model)
    config = LAMMPS_CONFIG.replace(
        'model_files = ["mace_nio.pt"]', 'model_files = ["nio.pb"]'
    ) + f'\n[potential.model_hashes]\n"nio.pb" = "{digest}"\n'
    job = parse_job(write(tmp_path, config))
    (key,) = job.potential.model_hashes
    assert Path(key).is_absolute() and Path(key).name == "nio.pb"
    # Before the fix this raised ProvenanceError although the hash was right.
    manifest = provenance.build_manifest(job=job, registration=resolve_bridge("lammps", "lammps"))
    assert list(manifest["potential"]["model_sha256_observed"].values()) == [digest]
    model.write_bytes(b"a different graph")
    with pytest.raises(ModelIntegrityError):
        provenance.build_manifest(job=job, registration=resolve_bridge("lammps", "lammps"))


def test_a_hash_for_an_unlisted_file_is_refused_at_parse_time():
    with pytest.raises(ConfigError, match="not in potential.model_files"):
        LammpsMlipPotentialSpec(
            pair_style="deepmd m.pb",
            pair_coeff=("* *",),
            type_map={1: "Ni"},
            model_paths=(Path("m.pb"),),
            model_hashes={"other.pb": "0" * 64},
        )


@pytest.mark.parametrize("name", ["pair_mace", ""])
def test_implementation_names_are_checked(name):
    """The registered name is 'pair-mace'; the old docstring spelled it 'pair_mace'."""
    from nio_md_prep.mlip.registry import build_bridge

    if not name:
        with pytest.raises(ConfigError, match="non-empty"):
            MacePotentialSpec(declared_elements=("Ni",), implementation=(name,))
        return
    spec = MacePotentialSpec(declared_elements=("Ni",), implementation=(name,))
    with pytest.raises(ConfigError, match="no engine registers"):
        build_bridge(spec, EngineSpec(kind="lammps"))
