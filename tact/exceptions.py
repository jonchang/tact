class TactError(Exception):
    """Base class for errors raised by TACT."""


class DisjointConstraintError(TactError):
    """Exception raised when a set of constraints lead to a disjoint implied age interval."""
