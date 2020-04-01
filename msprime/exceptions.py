#
# Copyright (C) 2017 University of Oxford
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msprime is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msprime.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Exceptions defined in mprime.
"""


class MsprimeException(Exception):
    """
    Superclass of all exceptions thrown.
    """


class DuplicatePositionsError(MsprimeException):
    """
    Duplicate positions in the list of sites.
    """


class FileFormatError(MsprimeException):
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


class ProvenanceValidationError(MsprimeException):
    """
    A JSON document did non validate against the provenance schema.
    """
