#
# Copyright (C) 2015-2020 University of Oxford
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
Core functions and classes used throughout msprime.
"""
import os
import random

from msprime import _msprime

__version__ = "undefined"
try:
    from . import _version

    __version__ = _version.version
except ImportError:
    pass


# Make sure the GSL error handler is turned off so that we can be sure that
# we don't abort on errors. This can be reset by using the function
# _msprime.restore_gsl_error_handler(), which will set the error handler to
# the value that it had before this function was called.
_msprime.unset_gsl_error_handler()


# Some machinery here for generating default random seeds. We need a map
# indexed by process ID here because we cannot use a global variable
# to store the state across multiple processes. Copy-on-write semantics
# for child processes means that they inherit the state of the parent
# process, so if we just keep a global variable without indexing by
# PID, child processes will share the same random generator as the
# parent.

_seed_rng_map = {}


def get_seed_rng():
    return _seed_rng_map.get(os.getpid(), None)


def clear_seed_rng():
    _seed_rng_map.pop(os.getpid(), None)


def get_random_seed():
    global _seed_rng_map
    pid = os.getpid()
    if pid not in _seed_rng_map:
        # If we don't provide a seed to Random(), Python will seed either
        # from a system source of randomness (i.e., /dev/urandom) or the
        # current time if this is not available. Thus, our seed rng should
        # be unique, even across different processes.
        _seed_rng_map[pid] = random.Random()
    return _seed_rng_map[pid].randint(1, 2 ** 32 - 1)


def isinteger(value):
    """
    Returns True if the specified value can be converted losslessly to an
    integer.
    """
    try:
        int_val = int(value)
        float_val = float(value)
        return int_val == float_val
    except (ValueError, TypeError):
        return False
