"""
Exceptions defined in tskit.
"""

# Duplicate the low-level library error so that at least client code doesn't need to
# import _tskit. But, this means that we do have

from _tskit import LibraryError


class TskitException(Exception):
    """
    Superclass of all exceptions thrown.
    """


class DuplicatePositionsError(TskitException):
    """
    Duplicate positions in the list of sites.
    """


class FileFormatError(TskitException):
    """
    Some file format error was detected.
    """


class VersionTooNewError(FileFormatError):
    """
    The version of the file is too new and cannot be read by the library.
    """


class VersionTooOldError(FileFormatError):
    """
    The version of the file is too old and cannot be read by the library.
    """


class ProvenanceValidationError(TskitException):
    """
    A JSON document did non validate against the provenance schema.
    """
