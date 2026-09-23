"""Dataset errors.

``DatasetError`` subclasses ``ValueError`` so the top-level CLI's existing
``except (ValueError, ...)`` handler reports it as ``error: ...`` with exit
status 2, like every other deliberate failure in ``nio-md-prep``.
"""

from __future__ import annotations


class DatasetError(ValueError):
    """A deliberate, user-actionable failure of a dataset command."""


class OverwriteRefusedError(DatasetError, FileExistsError):
    """The output location already holds data and ``--force`` was not given."""


class LeakageError(DatasetError):
    """A split would place correlated or duplicate frames in different subsets."""


class DependencyMissingError(DatasetError):
    """An optional dependency (numpy, ASE) needed by this operation is absent."""
