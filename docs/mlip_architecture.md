# Machine-learned interatomic potentials: architecture

`nio_md_prep.mlip` is a new, isolated subsystem for running machine-learned
interatomic potentials. It does not change how a classical study is built.
`agglomeration.py`, `build.py`, `lammps.py`, `lego.py`, `lego2.py` and
`analysis/` are untouched, and no existing workflow gains an MLIP code path.
The two live side by side.

This document is the contract: the compatibility matrix, the units and energy
conventions, the dependencies, what a job directory must record, how to extend
the subsystem, and what is explicitly not supported.

## The architectural rule

**Potential type and dynamics engine are independent concepts.**

A MACE model does not know whether ASE, LAMMPS or OpenMM will integrate it. An
engine does not know which model family it is driving. Their joining is a third
thing — a *bridge* — and every bridge lives in one table, in
`nio_md_prep.mlip.registry`, rather than as scattered `if potential == "mace"`
tests through preparation code.

The practical consequence: adding DeepMD does not touch any engine, adding a
new engine does not touch any potential, and asking for a combination that
cannot exist produces a typed error from one place.

## Compatibility matrix

| Potential | ASE | LAMMPS | OpenMM |
|---|---|---|---|
| **MACE** | supported directly (MACE's own ASE calculator) | supported via MACE/ML-IAP, with a MACE pair style as a fallback | supported via OpenMM-ML |
| **LAMMPS-native MLIP** | supported through the ASE → LAMMPS bridge | supported natively | **unsupported** unless that model has a separate OpenMM/ASE implementation |

`mock` is a seventh cell: an analytic Lennard-Jones stand-in on ASE, used to
exercise this machinery in CI on a machine with no backend. It is diagnostic
infrastructure and says so in every report it produces.

The matrix is executable. `nio-md-prep mlip inspect` prints it with per-route
availability, `registry.compatibility_matrix()` returns it as data, and
`tests/test_mlip_registry.py` asserts every cell.

### Why `lammps → openmm` is unsupported

A LAMMPS-native MLIP is compiled C++ inside LAMMPS. OpenMM cannot call it, and
there is deliberately **no generic "LAMMPS pair style → OpenMM" adapter** here;
writing one would mean reimplementing someone else's model and calling the
result the same potential.

The cell is registered *explicitly*, with a reason, so the failure is a clean
`UnsupportedCombinationError` naming the routes that do work — not an import
traceback from a half-written adapter:

```text
$ nio-md-prep mlip validate lammps-on-openmm.toml
error: potential 'lammps' cannot run on engine 'openmm': a LAMMPS-native MLIP
is compiled C++ inside LAMMPS; OpenMM cannot call it, and no generic
LAMMPS-pair-style-to-OpenMM adapter exists or is being faked here. If the
underlying model has a separate OpenMM (or ASE) implementation, register that
implementation as its own potential kind -- the way MACE is registered --
rather than routing it through LAMMPS
Supported instead: lammps -> ase (...), lammps -> lammps (...)
```

If a model *does* have a separate OpenMM implementation, it becomes its own
potential kind with its own bridges, exactly as MACE is. That is the extension
path, and it is honest about what is actually being run.

### Why MACE → LAMMPS has two implementations

`mliap` is preferred and chosen by default: the newer MACE/LAMMPS interface,
which brings GPU acceleration, multi-GPU inference and atomic virials. It is
the only MACE/LAMMPS route here that advertises a per-atom energy
decomposition.

`pair-mace` is kept because not every cluster build has ML-IAP, and is declared
more conservatively — per-atom virials are not assumed. Select it with
`potential.implementation = ["pair-mace"]`, which is a *preference order*, not
a hard selection; an unrecognised preference falls back to the default rather
than failing.

Both routes run the **exported** LAMMPS model, not the checkpoint the ASE
calculator loads. That export is a real step where numbers can move, which is
why MACE advises care when benchmarking LAMMPS output against the
corresponding ASE calculator, and why cross-engine equivalence is an explicit
acceptance test below rather than an assumption.

## Scope

Implemented: **full-system MLIP**. The potential is applied to every atom.

Deliberately out of scope in this release, and not partially implemented:

- ML/MM partitioning and classical + MLIP overlays
- fixed ML regions
- electrostatic embedding
- bonds crossing an ML/MM boundary
- the NiO/phosphonate head-group hybrid idea

OpenMM-ML can build mixed systems, and there is still active work around bonds
crossing an ML/MM boundary; none of that is exposed here. The schema keeps the
hook: `simulation.region` accepts `"all"` and `"selection"`, but `"selection"`
is rejected with an explanation rather than silently treated as `"all"`, so a
later ML-region capability needs no breaking schema change.

`smoke-md` is **not** a production MD driver and refuses to become one: it caps
at 500 steps and says so. Production deposition and relaxation still run
through the existing classical workflows.

## Module map

```text
src/nio_md_prep/mlip/
  errors.py          the error taxonomy; every deliberate failure is one of these
  units.py           unit systems and the two energy conventions
  capabilities.py    CapabilitySet, RequirementSet, negotiation, element coverage
  specs.py           PotentialSpec, EngineSpec, SimulationSpec, StructureSpec, JobSpec
  config.py          TOML parsing; independent of execution
  structures.py      ASE Atoms as the interchange structure, plus converters
  registry.py        the bridge registry: the compatibility matrix, in code
  results.py         PotentialResult, TrajectoryResult, ComparisonResult
  provenance.py      mlip_manifest.json
  environment.py     availability probing that never imports what it probes
  jobs.py            inspect / validate / singlepoint / smoke-md / compare
  cli.py             the `nio-md-prep mlip ...` command group
  potentials/        model families: mace.py, lammps_mlip.py, mock.py
  engines/           ase_engine.py, lammps_engine.py, openmm_engine.py
  bridges/           one module per cell of the matrix
```

Everything above `potentials/`, `engines/` and `bridges/` is pure data and
logic. Constructing a spec, parsing a configuration, resolving a bridge or
building a manifest imports nothing heavier than `pathlib` and `tomllib`.

## Lazy imports

Installing plain `nio-md-prep` must not require PyTorch, CUDA, OpenMM or
LAMMPS, so:

- backends are imported **inside** the functions that need them, never at
  module scope;
- bridge classes are registered by import string and loaded only when a bridge
  is actually built;
- availability is probed with `importlib.util.find_spec`, which answers
  "could this be imported" without importing it;
- `nio_md_prep.mlip.__getattr__` resolves the public names on first access.

`tests/test_mlip_registry.py::test_resolution_imports_no_backend` and
`tests/test_mlip_config.py::test_parsing_imports_no_backend` assert this, so a
stray top-level `import torch` fails CI.

### Dependencies

| Extra | Pulls in | For |
|---|---|---|
| *(none)* | `openpyxl` | the classical workflows, unchanged |
| `mlip` | `ase`, `numpy` | specs, structures, the mock potential |
| `mace` | `mlip` + `mace-torch`, `torch` | MACE on ASE and on LAMMPS |
| `openmm` | `mace` + `openmm`, `openmm-ml` | MACE on OpenMM |

LAMMPS is not a pip dependency: the engine uses the LAMMPS Python module when
it is importable and otherwise drives an `lmp` executable, and if neither
exists it writes the rendered deck and says where it put it. Install the torch
build matching the target machine's CUDA separately; no extra pins one.

## Units and energy conventions

Two independent things can make two correct engines disagree about the same
model and the same geometry. Both are treated as first-class metadata.

### Units

| System | Energy | Length | Notes |
|---|---|---|---|
| canonical | eV | Å | the comparison system; stress in eV/Å³ |
| ASE | eV | Å | already canonical |
| LAMMPS `metal` | eV | Å | already canonical, but pressure is in **bars** |
| LAMMPS `real` | kcal/mol | Å | **not** canonical |
| OpenMM | kJ/mol | nm | **not** canonical |

Every `PotentialResult` is in canonical units and records which native units it
came from, so a wrong number can be traced to a conversion rather than guessed
at. Engines still run natively. An OpenMM kJ/mol/nm output is never compared
with an ASE eV/Å output without passing through `units.to_canonical`.

Two LAMMPS-specific conversions are explicit, never implicit: LAMMPS reports
*pressure* in bars where canonical stress is eV/Å³, and pressure is the
negative of stress under ASE's convention. `units = "lj"` is refused outright —
it has no absolute energy scale, so treating it as eV would be a fiction.

Conversion factors are built from the exact SI-defining constants (the
elementary charge and the Avogadro constant), not copied from a table.
`tests/test_mlip_units.py` cross-checks them against `ase.units`' CODATA 2018
set and documents the ninth-figure gap against ASE's CODATA 2014 default.

### Energy conventions

| Convention | Meaning | Who reports it |
|---|---|---|
| `total` | includes per-element atomic self-energies (the isolated-atom references, `E0`) | MACE's ASE calculator; MACE LAMMPS pair styles |
| `interaction` | atomic self-energies removed | **OpenMM-ML's MACE default** (`returnEnergyType="interaction_energy"`) |

This matters enormously. On a few thousand atoms the offset is thousands of
eV, so an ASE/OpenMM comparison can look catastrophically wrong while the
forces agree to machine precision — a per-element constant has zero gradient.

So:

- the OpenMM bridge asks for whichever convention the job declares, rather than
  inheriting OpenMM-ML's default silently, and records what it got;
- converting between conventions needs per-element `E0` values, taken from
  `potential.atomic_reference_energies` or discovered from the model file. When
  neither is available, `compare_results` **refuses** rather than returning a
  meaningless difference;
- forces and stresses are returned unchanged by a convention conversion,
  because the offset is a composition-dependent constant.

## Configuration

Parsing is independent of execution. A configuration destined for a GPU node
can be checked on a laptop with nothing installed. Unknown keys are errors, not
warnings — a typo in `temperature_K` must not quietly produce a 0 K trajectory.

```toml
[potential]
kind = "mace"
model = "/path/to/nio_phosphonate.model"
sha256 = "…"                       # verified before the job runs
elements = ["Ni", "O", "P", "C", "H"]
device = "cuda"                    # cpu | cuda | cuda:1 | mps
precision = "float32"              # float32 | float64
cutoff_angstrom = 6.0              # optional; cross-checked against the model
implementation = ["mliap"]         # preference order, not a hard selection
energy_convention = "total"

[potential.atomic_reference_energies]  # E0, in eV; needed to cross conventions
Ni = -5.7801

[engine]
kind = "ase"                       # ase | lammps | openmm
platform = "CUDA"                  # OpenMM only
precision = "float64"
threads = 8
executable = "lmp"                 # LAMMPS only

[simulation]
task = "md"                        # singlepoint | optimize | md
ensemble = "nvt"                   # nve | nvt | npt
temperature_K = 400
timestep_fs = 0.5
steps = 10000
seed = 12345
region = "all"                     # "selection" is reserved, not implemented

[structure]
path = "prepared/nio_slab.xyz"
format = "auto"                    # …or "lammps-data-nio" for a prepared build
```

The LAMMPS-native alternative. Nothing here is DeepMD-specific: swapping
`pair_style` for `pace …`, `mliap unified …` or a MACE pair style needs no
code change.

```toml
[potential]
kind = "lammps"
framework = "deepmd"               # a label for provenance, not a discriminator
pair_style = "deepmd nio_phosphonate.pb"
pair_coeff = ["* *"]
units = "metal"                    # metal | real
atom_style = "atomic"
required_packages = ["USER-DEEPMD / DEEPMD"]
model_files = ["nio_phosphonate.pb"]

[potential.type_map]               # mandatory: type -> element
1 = "Ni"
2 = "O"
3 = "P"
4 = "C"
5 = "H"

[potential.model_hashes]
"nio_phosphonate.pb" = "…"

[engine]
kind = "lammps"
```

`type_map` is mandatory because it is what makes element coverage checkable for
a model whose element list is otherwise buried in a binary file, and what maps
an ASE `Atoms` object onto LAMMPS types when the same potential is driven
through the ASE bridge.

Worked examples: `examples/mlip/`.

## Structures

The MLIP layer's interchange structure is an ASE `Atoms` object. For a
full-system MLIP, symbols, coordinates, cell and PBC are everything a potential
needs, and ASE is the most direct route into MACE.

The repository's topology-rich LAMMPS representation (`nio_md_prep.lammps`) is
**untouched**. Bonds, angles, dihedrals, impropers, molecule IDs and per-atom
charges stay where the classical workflows need them.
`structures.from_lammps_data` is a converter, one way, into the reduced view an
MLIP consumes — it drops topology deliberately, because carrying it would
suggest the conversion is reversible.

Every structure is hashed (`structures.structure_digest`) over exactly the
interchange view, rounded to 1e-8 Å so a text round trip does not change the
digest. That hash goes in the manifest.

### Element coverage

Validated before anything is constructed. A Ni/O/P/C/H model handed an extra
atom type raises `ElementCoverageError` naming the offending element and what
the model does cover — it never produces a plausible-looking number instead.

## Capability negotiation

Nothing assumes every model returns everything. A bridge advertises a
`CapabilitySet`; a job declares a `RequirementSet`; the resolver refuses before
a batch job is generated.

Capabilities: `energy`, `forces`, `stress`, `per_atom_energy`, `periodic`,
`gpu`, plus supported `elements`, `precisions`, `engines`, the native energy
convention and which conventions are reachable.

A route's capabilities are the **intersection** of the potential's and the
engine's: a MACE model can produce a stress, but OpenMM does not report one, so
`mace → openmm` cannot.

Requirements are derived from the simulation, each with a reason attached:

| Job | Requires |
|---|---|
| any | `energy` |
| `singlepoint` | `forces` |
| `optimize` | `energy` + `forces` |
| `md` (nve/nvt) | `energy` + `forces` |
| `md` (**npt**) | additionally `stress` — the cell is integrated against the virial |
| `compute_per_atom_energy` | `per_atom_energy` |
| a periodic input structure | `periodic` |

The reasons are what turn "unsupported" into something actionable:

```text
error: mace -> openmm (openmm-ml) cannot satisfy this job:
  - stress tensor / virial (stress) is not available: the npt ensemble
    integrates the simulation cell against the stress/virial, so an absent or
    untrustworthy virial route is fatal
```

## Provenance

Provenance comes before long MD, not after it. Once several MACE committee
members exist — and later DeepMD or DPA models — a directory full of
trajectories is worthless unless each one records which model produced it.

Every executing job writes `mlip_manifest.json` into its output directory
**before** the run, so a crash still leaves a record of the attempt, and
rewrites it with results afterwards. It contains:

- model SHA256 (declared *and* observed) and any missing model files
- input structure hash, composition, cell and PBC
- git commit and whether the working tree was dirty
- the potential implementation and the simulation engine, with the bridge's
  own notes
- package versions (`ase`, `numpy`, `torch`, `mace-torch`, `openmm`,
  `openmm-ml`, `lammps`, `deepmd-kit`, …), absent ones recorded as `null`
- CUDA/device information, probed without importing torch when torch is absent
- element mapping, including the LAMMPS type map and the `E0` values
- precision and energy convention (potential-native, requested, reported)
- thermostat/integrator settings and the seed
- the generated engine parameters

For LAMMPS the last item is literal: the exact rendered `pair_style` and
`pair_coeff` command strings are preserved, not instructions for re-deriving
them. For OpenMM it records `returnEnergyType`, the platform and the precision.

A model file whose SHA256 does not match the configuration stops the job with
`ModelIntegrityError`. Provenance that lies is worse than none.

## Commands

Three levels of commitment. Only the first two are production-ready in this
release.

```bash
# What can this machine run? What does this model contain?
nio-md-prep mlip inspect [CONFIG] [--json]

# Is this job possible? Resolves, negotiates and renders — executes nothing.
nio-md-prep mlip validate CONFIG [--structure PATH] [--json]

# Evaluate one geometry, write a manifest.
nio-md-prep mlip singlepoint CONFIG --output DIR [--structure PATH]

# A tiny diagnostic trajectory. Capped at 500 steps.
nio-md-prep mlip smoke-md CONFIG --output DIR [--structure PATH]

# The cross-engine acceptance test.
nio-md-prep mlip compare CONFIG --output DIR --engines ase,lammps,openmm
```

`mlip inspect` with no configuration reports the matrix, per-route
availability, installed packages, device information and the energy
conventions. With a configuration it adds the selected route and what could be
read from the model file.

There is deliberately no `mlip md`: this subsystem does not replace the
existing production deposition and relaxation workflows.

## The scientific acceptance test

Cross-engine single-point equivalence. For the same MACE model and the same
structure, evaluate energy and forces through ASE, LAMMPS and OpenMM wherever
they are installed, with OpenMM configured to use the same energy convention,
and record:

- ΔE/N (eV/atom)
- force RMSE (eV/Å)
- maximum force error (eV/Å)
- stress/virial differences where the route supports them

ASE is the reference. Energies are harmonised to one convention before being
subtracted, so what is reported is real disagreement between adapters rather
than OpenMM-ML's interaction-energy default.

**The tolerance is measured, not invented.** The first run on a machine with
two or more backends writes a candidate reference file with the agreement
actually achieved in double precision. A maintainer reviews those numbers and
commits them as `tests/data/mlip_cross_engine_tolerances.json`, after which the
test is an ordinary regression criterion. Until then it skips with an
explanation. MACE itself warns users to be careful benchmarking LAMMPS output
against the corresponding ASE calculator, so an assumed tolerance would be
exactly the wrong thing to encode.

## Testing

Ordinary CI runs the unmarked tests and needs no backend:

- `test_mlip_units.py` — unit conversions and energy conventions
- `test_mlip_registry.py` — the compatibility matrix, including the
  intentional `UnsupportedCombinationError`
- `test_mlip_config.py` — configuration parsing and spec validation
- `test_mlip_capabilities.py` — negotiation and element-coverage failures
- `test_mlip_provenance.py` — the manifest
- `test_mlip_structures.py` — the interchange structure and converters
- `test_mlip_lammps_deck.py` — deck rendering, which is pure
- `test_mlip_pipeline.py` — the whole pipeline, end to end, on the mock potential

Expensive tests are marked `mace`, `lammps`, `openmm` and `gpu`, and skip
themselves unless the backend is importable. No trained model is committed;
one is supplied through an environment variable:

```bash
export NIO_MD_TEST_MACE_MODEL=/models/nio_phosphonate.model
export NIO_MD_TEST_MACE_ELEMENTS=Ni,O,P,C,H
export NIO_MD_TEST_MACE_LAMMPS_MODEL=/models/nio_phosphonate-mliap.pt  # optional
pytest tests/test_mlip_backends.py -m "mace or lammps or openmm"
```

## Extending

### How to add DeepMD (or PACE, or any LAMMPS pair style)

Nothing. A DeepMD model is already expressible:

```toml
[potential]
kind = "lammps"
framework = "deepmd"
pair_style = "deepmd nio.pb"
pair_coeff = ["* *"]
```

`LammpsMlipPotentialSpec` is not tied to a framework — `pair_style` is the
datum and `framework` is a label. Both the native LAMMPS route and the
ASE → LAMMPS route work immediately.

Two optional refinements:

1. Add a row to `KNOWN_FRAMEWORKS` in `potentials/lammps_mlip.py` so the
   required LAMMPS package is inferred and the capability advertisement is
   accurate (does this style compute `eatom`? a virial?). Capabilities there
   are declared **conservatively**: a job asking for something the style does
   not compute must fail loudly rather than receive zeros.
2. If DeepMD gains a first-class non-LAMMPS route (a native ASE calculator, or
   an OpenMM implementation), give it its own potential kind — a
   `DeepMDPotentialSpec` in `specs.py`, an adapter in `potentials/`, and one
   `register(...)` call per new cell. Do **not** try to route it through the
   LAMMPS spec; that would claim LAMMPS is running when it is not.

### How to add a new engine

1. Add its kind to `ENGINE_KINDS` in `specs.py`.
2. Add `engines/<name>_engine.py` with an `EngineRuntime` subclass declaring
   what the engine can *carry* — and honestly: if it reports no stress, say so,
   the way `openmm_engine` does. Add the single-point and MD functions.
3. Add one bridge module per potential family that can reach it, and one
   `register(...)` call each in `registry.py`.
4. For any family that **cannot** reach it, add a `register(..., status=UNSUPPORTED,
   reason=...)` row. An explicit refusal with a reason is part of the matrix,
   not an omission from it.
5. Add its native units to `units.py` if they are not already canonical, and
   convert at the bridge boundary — never in the engine's caller.

### How to add a new potential family

1. A `PotentialSpec` subclass in `specs.py`, plus an entry in
   `POTENTIAL_SPECS` and a key map in `config.py`.
2. An adapter in `potentials/` answering three questions: what can this model
   produce, can it run here, what does its model file say about itself.
3. One bridge module per reachable engine, and one `register(...)` call each.

No existing engine, bridge or spec changes.

## Explicitly unsupported

| Request | Result |
|---|---|
| `potential = "lammps"`, `engine = "openmm"` | `UnsupportedCombinationError`, with the reason and the working alternatives |
| `simulation.region = "selection"` | `ConfigError`: reserved, not implemented — full-system MLIP only |
| NPT through OpenMM | `CapabilityError`: no stress tensor available to integrate the cell against |
| per-atom energies through the ASE → LAMMPS bridge | `CapabilityError`: `LAMMPSlib` does not expose them; use `engine = "lammps"` |
| per-atom energies through `pair-mace` | `CapabilityError`: not assumed for that route |
| `units = "lj"` | `UnitError`: no absolute energy scale, so no canonical mapping |
| comparing two conventions with no `E0` values | `EnergyConventionError`, rather than a meaningless number |
| a structure with an element the model does not cover | `ElementCoverageError`, before anything runs |
| a model file whose SHA256 does not match | `ModelIntegrityError`, before anything runs |
| `mlip smoke-md` with more than 500 steps | `ConfigError`: this is a diagnostic, not a production MD driver |
| ML/MM, fixed ML regions, electrostatic embedding, ML/MM bond crossing | not implemented, and not partially implemented |
