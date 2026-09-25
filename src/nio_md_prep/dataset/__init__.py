"""Audited VASP -> MLIP dataset export (``nio-md-prep dataset ...``).

This package turns completed VASP calculations into MACE/ASE-compatible
extended-XYZ datasets with a machine-readable account of every discovered run
and ionic frame. It depends only on the standard library and numpy for
discovery, parsing, classification and splitting; ASE is imported lazily, and
only to write and re-read extended XYZ. No module here may import torch, MACE,
OpenMM or LAMMPS.

The scientific contract is documented in ``docs/dataset-export.md``. In short:
labels and geometry always come from the same ``<calculation>`` block of
``vasprun.xml``; frames without positive electronic-convergence evidence are
never accepted in strict mode; incompatible reference settings are kept in
separate pools; and correlated frames share a split group.
"""

SCHEMA_VERSION = 1
PARSER_NAME = "nio-md-prep.vasprun-stream"
PARSER_VERSION = "1.0"

__all__ = ["PARSER_NAME", "PARSER_VERSION", "SCHEMA_VERSION"]
