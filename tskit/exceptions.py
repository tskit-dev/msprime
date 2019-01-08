"""
Exceptions defined in tskit.
"""
from _tskit import TskitException
from _tskit import LibraryError
from _tskit import FileFormatError
from _tskit import VersionTooNewError
from _tskit import VersionTooOldError

# Some exceptions are defined in the low-level module. In particular, the
# superclass of all exceptions for tskit is defined here. We define the
# docstrings here to avoid difficulties with compiling C code on
# readthedocs.

# TODO finalise this when working out the docs structure for tskit on rtd.

try:
    TskitException.__doc__ = "Superclass of all exceptions defined in tskit."
    LibraryError.__doc__ = "Generic low-level error raised by the C library."
    FileFormatError.__doc__ = "An error was detected in the file format."
    VersionTooNewError.__doc__ = """
    The version of the file is too new and cannot be read by the library.
    """
    VersionTooOldError.__doc__ = """
    The version of the file is too old and cannot be read by the library.
    """
except AttributeError:
    # Python2 throws attribute error. Ignore.
    pass


class DuplicatePositionsError(TskitException):
    """
    Duplicate positions in the list of sites.
    """


class ProvenanceValidationError(TskitException):
    """
    A JSON document did non validate against the provenance schema.
    """
