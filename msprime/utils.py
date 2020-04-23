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
Various utilities shared across modules.
"""
import os
import random


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
    return _seed_rng_map[pid].randint(1, 2**32 - 1)


# Note, this is no longer needed as we have Python 3.6 minimum.
# https://github.com/tskit-dev/msprime/issues/960
def almost_equal(a, b, rel_tol=1e-9, abs_tol=0.0):
    """
    Returns true if the specified pair of numbers are equal to
    within the specified tolerances.

    The signature and implementation are taken from PEP 485,
    https://www.python.org/dev/peps/pep-0485/
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
