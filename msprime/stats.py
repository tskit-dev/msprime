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

import _msprime


def check_numpy():
    if not _numpy_imported:
        raise RuntimeError("numpy is required for statistical calculations.")


class LdCalculator(object):
    """
    Class for calculating `linkage disequilibrium
    <https://en.wikipedia.org/wiki/Linkage_disequilibrium>`_ coefficients
    between pairs of mutations in a :class:`.TreeSequence`. This class requires
    the `numpy <http://www.numpy.org/>`_ library.

    This class supports multithreaded access using the Python :mod:`threading`
    module. Separate instances of :class:`.LdCalculator` referencing the
    same tree sequence can operate in parallel in multiple threads.
    See the :ref:`sec-tutorial-threads` section in the :ref:`sec-tutorial`
    for an example of how use multiple threads to calculate LD values
    efficiently.

    :param TreeSequence tree_sequence: The tree sequence containing the
        mutations we are interested in.
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
        Returns the value of the :math:`r^2` statistic between the pair of
        mutations at the specified indexes. This method is *not* an efficient
        method for computing large numbers of pairwise; please use either
        :meth:`.get_r2_array` or :meth:`.get_r2_matrix` for this purpose.

        :param int a: The index of the first mutation.
        :param int b: The index of the second mutation.
        :return: The value of :math:`r^2` between the mutations at indexes
            ``a`` and ``b``.
        :rtype: float
        """
        return self._ll_ld_calculator.get_r2(a, b)

    def get_r2_array(
            self, a, direction=1, max_mutations=None, max_distance=None):
        """
        Returns the value of the :math:`r^2` statistic between the focal
        mutation at index :math:`a` and a set of other mutations. The method
        operates by starting at the focal mutation and iterating over adjacent
        mutations (in either the forward or backwards direction) until either a
        maximum number of other mutations have been considered (using the
        ``max_mutations`` parameter), a maximum distance in sequence
        coordinates has been reached (using the ``max_distance`` parameter) or
        the start/end of the sequence has been reached. For every mutation
        :math:`b` considered, we then insert the value of :math:`r^2` between
        :math:`a` and :math:`b` at the corresponding index in an array, and
        return the entire array. If the returned array is :math:`x` and
        ``direction`` is :const:`msprime.FORWARD` then :math:`x[0]` is the
        value of the statistic for :math:`a` and :math:`a + 1`, :math:`x[1]`
        the value for :math:`a` and :math:`a + 2`, etc. Similarly, if
        ``direction`` is :const:`msprime.REVERSE` then :math:`x[0]` is the
        value of the statistic for :math:`a` and :math:`a - 1`, :math:`x[1]`
        the value for :math:`a` and :math:`a - 2`, etc.

        :param int a: The index of the focal mutation.
        :param int direction: The direction in which to travel when
            examining other mutations. Must be either
            :const:`msprime.FORWARD` or :const:`msprime.REVERSE`. Defaults
            to :const:`msprime.FORWARD`.
        :param int max_mutations: The maximum number of mutations to return
            :math:`r^2` values for. Defaults to as many mutations as
            possible.
        :param float max_distance: The maximum absolute distance between
            the focal mutation and those for which :math:`r^2` values
            are returned.
        :return: An array of double precision floating point values
            representing the :math:`r^2` values for mutations in the
            specified direction.
        :rtype: numpy.ndarray
        :warning: For efficiency reasons, the underlying memory used to
            store the returned array is shared between calls. Therefore,
            if you wish to store the results of a single call to
            ``get_r2_array()`` for later processing you **must** take a
            copy of the array!
        """
        if max_mutations is None:
            max_mutations = -1
        if max_distance is None:
            max_distance = sys.float_info.max
        num_values = self._ll_ld_calculator.get_r2_array(
            self._buffer, a, direction=direction,
            max_mutations=max_mutations, max_distance=max_distance)
        return np.frombuffer(self._buffer, "d", num_values)

    def get_r2_matrix(self):
        """
        Returns the complete :math:`m \\times m` matrix of pairwise
        :math:`r^2` values in a tree sequence with :math:`m` mutations.

        :return: An 2 dimensional square array of double precision
            floating point values representing the :math:`r^2` values for
            all pairs of mutations.
        :rtype: numpy.ndarray
        """
        m = self._tree_sequence.get_num_mutations()
        A = np.ones((m, m), dtype=float)
        for j in range(m - 1):
            a = self.get_r2_array(j)
            A[j, j + 1:] = a
            A[j + 1:, j] = a
        return A
