# Documentation

The root [`README.md`](../README.md) is intentionally a short project overview. This directory contains the working scientific and operational documentation.

## Start here

- [`project-design.md`](project-design.md) — scientific picture, supplied systems, composition logic, and classical force field.
- [`classical-md-workflow.md`](classical-md-workflow.md) — installation, primary/CoSAM construction, sequential deposition, temperature branches, wall retraction, and controls.
- [`operations.md`](operations.md) — validation, safe input refresh, campaign inventory, archiving, adding passivants, and provenance.

## Analysis

- [`coverage-analysis.md`](coverage-analysis.md) — coordinate-based projected coverage, uncertainty, void topology, and consolidated workbooks.
- [`interfacial-analysis.md`](interfacial-analysis.md) — Ni-site ownership, anchoring, orientations, RDFs, kinetics policy, dipole proxy limits, and structure-property summary generation.

## Geometric-control experiments

- [`lego-deposition.md`](lego-deposition.md) — coverage-guided 1D gap seeding.
- [`lego2-deposition.md`](lego2-deposition.md) — localized periodic 2D-void seeding.

These are seeded accessibility controls, not replacements for the unbiased sequential-deposition calculation.

## MLIP training-data construction

- [`agglomeration-training.md`](agglomeration-training.md) — reproducible phosphonate agglomerates, mixed agglomerates, xTB optimization/MD/quench sampling, VASP training-set preparation, validation, campaign regeneration, and audit/status tools.

## Suggested reading paths

### Reproduce the classical SAM study

1. [`project-design.md`](project-design.md)
2. [`classical-md-workflow.md`](classical-md-workflow.md)
3. [`coverage-analysis.md`](coverage-analysis.md)
4. [`interfacial-analysis.md`](interfacial-analysis.md)
5. [`operations.md`](operations.md)

### Extend the study with another passivant

1. [`project-design.md`](project-design.md)
2. [Adding another passivant](operations.md#adding-another-passivant)
3. [`classical-md-workflow.md`](classical-md-workflow.md)

### Generate phosphonate cluster data for MLIP development

1. [`agglomeration-training.md`](agglomeration-training.md)
2. [`operations.md`](operations.md) for repository-wide validation and archiving conventions when needed.

## Documentation policy

Keep the root README general. New details should normally go into the closest specialist page here and be linked from this index. Long command sequences, campaign-specific recovery procedures, method definitions, and interpretation caveats should not accumulate in the root README.
