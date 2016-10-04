#
# Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
Module responsible for computing various statistics on tree sequences.
"""
from __future__ import division
from __future__ import print_function

import struct
import sys

try:
    import numpy as np
    _numpy_imported = True
except ImportError:
    _numpy_imported = False

import msprime
import _msprime


def check_numpy():
    if not _numpy_imported:
        raise RuntimeError("numpy is required for statistical calculations.")


class LdCalculator(object):
    """
    Class for calculating linkage disequilibrium coefficients between
    pairs of mutations in a :class:`.TreeSequence`.
    """

    def __init__(self, tree_sequence):
        check_numpy()
        self._tree_sequence = tree_sequence
        self._ll_ld_calculator = _msprime.LdCalculator(
            tree_sequence.get_ll_tree_sequence())
        item_size = struct.calcsize('d')
        self._buffer = bytearray(
            tree_sequence.get_num_mutations() * item_size)

    def get_r2(self, a, b):
        """
        Returns the value of the r^2 statistic between the pair of mutations
        at the specified indexes.
        """
        return self._ll_ld_calculator.get_r2(a, b)

    def get_r2_array(
            self, a, direction=msprime.FORWARD, max_mutations=None,
            max_distance=None):
        if max_mutations is None:
            max_mutations = -1
        if max_distance is None:
            max_distance = sys.float_info.max
        num_values = self._ll_ld_calculator.get_r2_array(
            self._buffer, a, direction=direction,
            max_mutations=max_mutations, max_distance=max_distance)
        return np.frombuffer(self._buffer, "d", num_values)

    def get_r2_matrix(self):
        m = self._tree_sequence.get_num_mutations()
        A = np.ones((m, m), dtype=float)
        for j in range(m - 1):
            a = self.get_r2_array(j)
            A[j, j + 1:] = a
            A[j + 1:, j] = a
        return A
