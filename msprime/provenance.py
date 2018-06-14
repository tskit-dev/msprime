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
from __future__ import print_function
from __future__ import division

import platform

import _msprime

__version__ = "undefined"
try:
    from . import _version
    __version__ = _version.version
except ImportError:
    pass


_gsl_version = _msprime.get_gsl_version()


def get_environment():
    """
    Returns a dictionary describing the environment in which msprime
    is currently running.
    """
    env = {
        "libraries": {
            "gsl": {
                "version": _gsl_version
            },
            "kastore": {
                # Hard coding this for now as there is no way to get version
                # information from the kastore C API. See
                # https://github.com/tskit-dev/kastore/issues/41
                # We could import the kastore module here and use its version,
                # but this is not the same as the C code we have compiler against.
                "version": (0, 1, 0)
            }
        },
        "os": {
            "system": platform.system(),
            "node": platform.node(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "python": {
            "implementation": platform.python_implementation(),
            "version": platform.python_version_tuple(),
        }
    }
    return env


def get_provenance_dict(command, parameters=None):
    """
    Returns a dictionary encoding an execution of msprime.

    Note: this format is incomplete and provisional.
    """
    document = {
        "software": "msprime",
        "version": __version__,
        "command": command,
        "parameters": parameters,
        "environment": get_environment()
    }
    return document
