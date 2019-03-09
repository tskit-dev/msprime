#
# Copyright (C) 2016-2017 University of Oxford
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
Common provenance methods used to determine the state and versions
of various dependencies and the OS.
"""
import tskit

import _msprime

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass


def get_provenance_dict(parameters=None):
    """
    Returns a dictionary encoding an execution of msprime conforming to the
    tskit provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {
            "name": "msprime",
            "version": __version__,
        },
        "parameters": parameters,
        "environment": get_environment()
    }
    return document


def _get_environment():
    gsl_version = ".".join(map(str, _msprime.get_gsl_version()))
    libraries = {"gsl": {"version": gsl_version}}
    return tskit.provenance.get_environment(extra_libs=libraries)


_environment = None


def get_environment():
    """
    Returns a dictionary describing the environment in which tskit
    is currently running.
    """
    # Everything here is fixed so we cache it
    global _environment
    if _environment is None:
        _environment = _get_environment()
    return _environment
