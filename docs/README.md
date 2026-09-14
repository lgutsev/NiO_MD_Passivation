# Documentation

The root [`README.md`](../README.md) provides the project overview. This directory contains the scientific, analysis, and operational documentation needed to reproduce and extend the workflows.

## Start here

- [`project-design.md`](project-design.md) — scientific model, supplied systems, composition logic, and classical force field.
- [`classical-md-workflow.md`](classical-md-workflow.md) — installation, primary/CoSAM construction, sequential deposition, temperature branches, wall retraction, and production controls.
- [`operations.md`](operations.md) — validation, safe input refresh, campaign inventory, archiving, adding passivants, and provenance.

## Analysis

- [`coverage-analysis.md`](coverage-analysis.md) — coordinate-based projected coverage, uncertainty, periodic void topology, and consolidated workbooks.
- [`interfacial-analysis.md`](interfacial-analysis.md) — Ni-site ownership, anchoring, orientations, RDFs, kinetics policy, dipole-proxy limits, and structure-property summaries.

## Gap-seeded accessibility controls

- [`lego-deposition.md`](lego-deposition.md) — **Gap-seeded I:** coverage-guided placement into the widest persistent low-coverage stripe, followed by laterally unconstrained dynamics.
- [`lego2-deposition.md`](lego2-deposition.md) — **Gap-seeded II:** placement into the largest persistent periodic 2D void, followed by laterally unconstrained dynamics.

These are deliberately biased initial-condition controls. They test whether secondary molecules can occupy persistent gaps when the search/access problem is reduced; they are not replacements for the unbiased sequential-deposition calculation and should not be interpreted as spontaneous self-assembly kinetics.

The implementation retains historical internal identifiers such as `lego`, `lego2`, `prepared-lego`, and `prepared-lego2` so existing scripts and production data remain compatible. Publication-facing and scientific documentation uses the names **Gap-seeded I** and **Gap-seeded II**.

## MLIP training-data construction

- [`agglomeration-training.md`](agglomeration-training.md) — reproducible phosphonate agglomerates, mixed agglomerates, xTB optimization/MD/quench sampling, VASP training-set preparation, validation, campaign regeneration, and audit/status tools.

## Suggested reading paths

### Reproduce the classical SAM study

1. [`project-design.md`](project-design.md)
2. [`classical-md-workflow.md`](classical-md-workflow.md)
3. [`coverage-analysis.md`](coverage-analysis.md)
4. [`interfacial-analysis.md`](interfacial-analysis.md)
5. [`operations.md`](operations.md)

### Reproduce the gap-seeded controls

1. [`coverage-analysis.md`](coverage-analysis.md)
2. [`lego-deposition.md`](lego-deposition.md) for Gap-seeded I
3. [`lego2-deposition.md`](lego2-deposition.md) for Gap-seeded II
4. [`interfacial-analysis.md`](interfacial-analysis.md) for matched structural comparison

### Extend the study with another passivant

1. [`project-design.md`](project-design.md)
2. [Adding another passivant](operations.md#adding-another-passivant)
3. [`classical-md-workflow.md`](classical-md-workflow.md)

### Generate phosphonate cluster data for MLIP development

1. [`agglomeration-training.md`](agglomeration-training.md)
2. [`operations.md`](operations.md) for repository-wide validation and archiving conventions when needed.

## Documentation policy

Keep the root README general. New details should normally go into the closest specialist page here and be linked from this index. Long command sequences, campaign-specific recovery procedures, method definitions, and interpretation caveats should not accumulate in the root README.
