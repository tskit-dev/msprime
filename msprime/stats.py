#
# Copyright (C) 2016 University of Oxford
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

import threading
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
    See the :ref:`sec_tutorial_threads` section in the :ref:`sec_tutorial`
    for an example of how use multiple threads to calculate LD values
    efficiently.

    .. note:: This class does not currently support sites that have more than one
        mutation. Using it on such a tree sequence will raise a LibraryError with
        an "Unsupported operation" message.

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
        # To protect low-level C code, only one method may execute on the
        # low-level objects at one time.
        self._instance_lock = threading.Lock()

    def get_r2(self, a, b):
        # Deprecated alias for r2(a, b)
        return self.r2(a, b)

    def r2(self, a, b):
        """
        Returns the value of the :math:`r^2` statistic between the pair of
        mutations at the specified indexes. This method is *not* an efficient
        method for computing large numbers of pairwise; please use either
        :meth:`.r2_array` or :meth:`.r2_matrix` for this purpose.

        :param int a: The index of the first mutation.
        :param int b: The index of the second mutation.
        :return: The value of :math:`r^2` between the mutations at indexes
            ``a`` and ``b``.
        :rtype: float
        """
        with self._instance_lock:
            return self._ll_ld_calculator.get_r2(a, b)

    def get_r2_array(self, a, direction=1, max_mutations=None, max_distance=None):
        # Deprecated alias for r2_array
        return self.r2_array(a, direction, max_mutations, max_distance)

    def r2_array(self, a, direction=1, max_mutations=None, max_distance=None):
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
        with self._instance_lock:
            num_values = self._ll_ld_calculator.get_r2_array(
                self._buffer, a, direction=direction,
                max_mutations=max_mutations, max_distance=max_distance)
        return np.frombuffer(self._buffer, "d", num_values)

    def get_r2_matrix(self):
        # Deprecated alias for r2_matrix
        return self.r2_matrix()

    def r2_matrix(self):
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


class GeneralStatCalculator(object):
    """
    A common class for BranchLengthStatCalculator and SiteStatCalculator -- those
    implemment different `tree_stat_vector()` methods, but given that
    general-purpose function, many statistics are computed in the same way.

    .. warning::
        This interface is still in beta, and may change in the future.
    """

    def __init__(self, tree_sequence):
        check_numpy()
        self.tree_sequence = tree_sequence

    def divergence(self, sample_sets, windows):
        """
        Finds the divergence between pairs of samples as described in
        mean_pairwise_tmrca_matrix (which uses this function).  Returns the
        upper triangle (including the diagonal) in row-major order, so if the
        output is `x`, then:

        >>> k=0
        >>> for w in range(len(windows)-1):
        >>>     for i in range(len(sample_sets)):
        >>>         for j in range(i,len(sample_sets)):
        >>>             trmca[i,j] = tmrca[j,i] = x[w][k]/2.0
        >>>             k += 1

        will fill out the matrix of mean TMRCAs in the `i`th window between (and
        within) each group of samples in `sample_sets` in the matrix `tmrca`.
        (This is because divergence is one-half TMRCA.) Alternatively, if
        `names` labels the sample_sets, the output labels are:

        >>> [".".join(names[i],names[j]) for i in range(len(names))
        >>>         for j in range(i,len(names))]

        :param list sample_sets: A list of sets of IDs of samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of the upper triangle of divergences in row-major
            order, including the diagonal.
        """
        ns = len(sample_sets)
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i]*(n[j]-x[j]))
                    for i in range(ns) for j in range(i, ns)]

        out = self.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        # corrects the diagonal for self comparisons
        # and note factor of two for tree length -> real time
        for w in range(len(windows)-1):
            k = 0
            for i in range(ns):
                for j in range(i, ns):
                    if i == j:
                        if n[i] == 1:
                            out[w][k] = np.nan
                        else:
                            out[w][k] /= float(n[i] * (n[i] - 1))
                    else:
                        out[w][k] /= float(n[i] * n[j])
                    k += 1

        return out

    def divergence_matrix(self, sample_sets, windows):
        """
        Finds the mean divergence  between pairs of samples from each set of
        samples and in each window. Returns a numpy array indexed by (window,
        sample_set, sample_set).  Diagonal entries are corrected so that the
        value gives the mean divergence for *distinct* samples, but it is not
        checked whether the sample_sets are disjoint (so offdiagonals are not
        corrected).  For this reason, if an element of `sample_sets` has only
        one element, the corresponding diagonal will be NaN.

        The mean divergence between two samples is defined to be the mean: (as
        a TreeStat) length of all edges separating them in the tree, or (as a
        SiteStat) density of segregating sites, at a uniformly chosen position
        on the genome.

        :param list sample_sets: A list of sets of IDs of samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of the upper triangle of mean TMRCA values in row-major
            order, including the diagonal.
        """
        x = self.divergence(sample_sets, windows)
        ns = len(sample_sets)
        nw = len(windows) - 1
        A = np.ones((nw, ns, ns), dtype=float)
        for w in range(nw):
            k = 0
            for i in range(ns):
                for j in range(i, ns):
                    A[w, i, j] = A[w, j, i] = x[w][k]
                    k += 1
        return A

    def Y3_vector(self, sample_sets, windows, indices):
        """
        Finds the 'Y' statistic between three sample_sets.  The sample_sets should
        be disjoint (the computation works fine, but if not the result depends
        on the amount of overlap).  If the sample_sets are A, B, and C, then the
        result gives the mean total length of any edge in the tree between a
        and the most recent common ancestor of b and c, where a, b, and c are
        random draws from A, B, and C respectively; or the density of mutations
        segregating a|bc.

        The result is, for each window, a vector whose k-th entry is
            Y(sample_sets[indices[k][0]], sample_sets[indices[k][1]],
              sample_sets[indices[k][2]]).

        :param list sample_sets: A list of *three* lists of IDs of samples: (A,B,C).
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of triples of indices of sample_sets.
        :return: A list of numeric vectors of length equal to the length of
            indices, computed separately on each window.
        """
        for u in indices:
            if not len(u) == 3:
                raise ValueError("All indices should be of length 3.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i] * (n[j] - x[j]) * (n[k] - x[k]))
                    for i, j, k in indices]

        out = self.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)

        # move this division outside of f(x) so it only has to happen once
        # corrects the diagonal for self comparisons
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                out[w][u] /= float(n[indices[u][0]] * n[indices[u][1]]
                                   * n[indices[u][2]])

        return out

    def Y2_vector(self, sample_sets, windows, indices):
        """
        Finds the 'Y' statistic for two groups of samples in sample_sets.
        The sample_sets should be disjoint (the computation works fine, but if
        not the result depends on the amount of overlap).
        If the sample_sets are A and B then the result gives the mean total length
        of any edge in the tree between a and the most recent common ancestor of
        b and c, where a, b, and c are random draws from A, B, and B
        respectively (without replacement).

        The result is, for each window, a vector whose k-th entry is
            Y2(sample_sets[indices[k][0]], sample_sets[indices[k][1]]).

        :param list sample_sets: A list of lists of IDs of leaves.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of pairs of indices of sample_sets.
        :return: A list of numeric vectors of length equal to the length of
            indices, computed separately on each window.
        """
        for u in indices:
            if not len(u) == 2:
                raise ValueError("All indices should be of length 2.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i] * (n[j] - x[j]) * (n[j] - x[j] - 1))
                    for i, j in indices]

        out = self.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                out[w][u] /= float(n[indices[u][0]] * n[indices[u][1]]
                                   * (n[indices[u][1]]-1))

        return out

    def Y1_vector(self, sample_sets, windows):
        """
        Finds the 'Y1' statistic within each set of samples in sample_sets. The
        sample_sets should be disjoint (the computation works fine, but if not
        the result depends on the amount of overlap).  For the sample set A, the
        result gives the mean total length of any edge in the tree between a
        and the most recent common ancestor of b and c, where a, b, and c are
        random draws from A, without replacement.

        The result is, for each window, a vector whose k-th entry is
            Y1(sample_sets[k]).

        :param list sample_sets: A list of sets of IDs of samples, each of length
            at least 3.
        :param iterable windows: The breakpoints of the windows (including
            start and end, so has one more entry than number of windows).
        :return: A list of numeric vectors of length equal to the length of
            sample_sets, computed separately on each window.
        """
        for x in sample_sets:
            if len(x) < 3:
                raise ValueError("All sample_sets should be of length at least 3.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(z * (m - z) * (m - z - 1)) for m, z in zip(n, x)]

        out = self.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)
        for w in range(len(windows)-1):
            for u in range(len(sample_sets)):
                out[w][u] /= float(n[u] * (n[u]-1) * (n[u]-2))

        return out

    def Y2(self, sample_sets, windows):
        return self.Y2_vector(sample_sets, windows, indices=[(0, 1)])

    def Y3(self, sample_sets, windows):
        """
        Finds the 'Y' statistic between the three groups of samples in
        sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).  If the
        sample_sets are A, B, and C, then the result gives the mean total
        length of any edge in the tree between a and the most recent common
        ancestor of b and c, where a, b, and c are random draws from A, B, and
        C respectively.

        :param list sample_sets: A list of *three* sets of IDs of samples: (A,B,C).
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of numeric values computed separately on each window.
        """
        return self.Y3_vector(sample_sets, windows, indices=[(0, 1, 2)])

    def f4_vector(self, sample_sets, windows, indices):
        """
        Finds the Patterson's f4 statistics between multiple subsets of four
        groups of sample_sets. The sample_sets should be disjoint (the computation
        works fine, but if not the result depends on the amount of overlap).

        :param list sample_sets: A list of four sets of IDs of samples: (A,B,C,D)
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of 4-tuples of indices of sample_sets.
        :return: A list of values of f4(A,B;C,D) of length equal to the length of
            indices, computed separately on each window.
        """
        for u in indices:
            if not len(u) == 4:
                raise ValueError("All tuples in indices should be of length 4.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i] * x[k] * (n[j] - x[j]) * (n[l] - x[l])
                          - x[i] * x[l] * (n[j] - x[j]) * (n[k] - x[k]))
                    for i, j, k, l in indices]

        out = self.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        # corrects the diagonal for self comparisons
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                out[w][u] /= float(n[indices[u][0]] * n[indices[u][1]]
                                   * n[indices[u][2]] * n[indices[u][3]])

        return out

    def f4(self, sample_sets, windows):
        """
        Finds the Patterson's f4 statistics between the four groups of samples
        in sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).

        :param list sample_sets: A list of four sets of IDs of samples: (A,B,C,D)
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of values of f4(A,B;C,D) computed separately on each window.
        """
        if not len(sample_sets) == 4:
            raise ValueError("sample_sets should be of length 4.")
        return self.f4_vector(sample_sets, windows, indices=[(0, 1, 2, 3)])

    def f3_vector(self, sample_sets, windows, indices):
        """
        Finds the Patterson's f3 statistics between multiple subsets of three
        groups of samples in sample_sets. The sample_sets should be disjoint (the
        computation works fine, but if not the result depends on the amount of
        overlap).

        f3(A;B,C) is f4(A,B;A,C) corrected to not include self comparisons.

        If A does not contain at least three samples, the result is NaN.

        :param list sample_sets: A list of sets of IDs of samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of triples of indices of sample_sets.
        :return: A list of values of f3(A,B,C) computed separately on each window.
        """
        for u in indices:
            if not len(u) == 3:
                raise ValueError("All tuples in indices should be of length 3.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i] * (x[i] - 1) * (n[j] - x[j]) * (n[k] - x[k])
                          - x[i] * (n[i] - x[i]) * (n[j] - x[j]) * x[k])
                    for i, j, k in indices]

        out = self.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                if n[indices[u][0]] == 1:
                    out[w][u] = np.nan
                else:
                    out[w][u] /= float(n[indices[u][0]] * (n[indices[u][0]]-1)
                                       * n[indices[u][1]] * n[indices[u][2]])

        return out

    def f3(self, sample_sets, windows):
        """
        Finds the Patterson's f3 statistics between the three groups of samples
        in sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).

        f3(A;B,C) is f4(A,B;A,C) corrected to not include self comparisons.

        :param list sample_sets: A list of *three* sets of IDs of samples: (A,B,C),
            with the first set having at least two samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of values of f3(A,B,C) computed separately on each window.
        """
        if not len(sample_sets) == 3:
            raise ValueError("sample_sets should be of length 3.")
        return self.f3_vector(sample_sets, windows, indices=[(0, 1, 2)])

    def f2_vector(self, sample_sets, windows, indices):
        """
        Finds the Patterson's f2 statistics between multiple subsets of pairs
        of samples in sample_sets. The sample_sets should be disjoint (the
        computation works fine, but if not the result depends on the amount of
        overlap).

        f2(A;B) is f4(A,B;A,B) corrected to not include self comparisons.

        :param list sample_sets: A list of sets of IDs of samples, each having at
            least two samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of pairs of indices of sample_sets.
        :return: A list of values of f2(A,C) computed separately on each window.
        """
        for u in indices:
            if not len(u) == 2:
                raise ValueError("All tuples in indices should be of length 2.")
        n = [len(x) for x in sample_sets]
        for xlen in n:
            if not xlen > 1:
                raise ValueError("All sample_sets must have at least two samples.")

        def f(x):
            return [float(x[i] * (x[i] - 1) * (n[j] - x[j]) * (n[j] - x[j] - 1)
                          - x[i] * (n[i] - x[i]) * (n[j] - x[j]) * x[j])
                    for i, j in indices]

        out = self.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                out[w][u] /= float(n[indices[u][0]] * (n[indices[u][0]]-1)
                                   * n[indices[u][1]] * (n[indices[u][1]] - 1))

        return out

    def f2(self, sample_sets, windows):
        """
        Finds the Patterson's f2 statistics between the three groups of samples
        in sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).

        f2(A;B) is f4(A,B;A,B) corrected to not include self comparisons.

        :param list sample_sets: A list of *two* sets of IDs of samples: (A,B),
            each having at least two samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of values of f2(A,B) computed separately on each window.
        """
        if not len(sample_sets) == 2:
            raise ValueError("sample_sets should be of length 2.")
        return self.f2_vector(sample_sets, windows, indices=[(0, 1)])

    def tree_stat(self, sample_sets, weight_fun):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a weight.  A branch in a tree is weighted by weight_fun(x),
        where x[i] is the number of samples in sample_sets[i] below that
        branch.  This finds the sum of all counted branches for each tree,
        and averages this across the tree sequence, weighted by genomic length.
        '''
        out = self.tree_stat_vector(sample_sets, lambda x: [weight_fun(x)])
        assert len(out) == 1 and len(out[0]) == 1
        return out[0][0]

    def tree_stat_windowed(self, sample_sets, weight_fun, windows=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a boolean.  A branch in a tree is weighted by weight_fun(x),
        where x[i] is the number of samples in sample_sets[i] below that
        branch.  This finds the sum of all counted branches for each tree,
        and averages this across the tree sequence, weighted by genomic length.
        '''
        out = self.tree_stat_vector(sample_sets, lambda x: [weight_fun(x)], windows)
        assert len(out[0]) == 1
        return [x[0] for x in out]


class BranchLengthStatCalculator(GeneralStatCalculator):
    """
    Class for calculating a broad class of tree statistics.  These are all
    calculated using :meth:``BranchLengthStatCalculator.tree_stat_vector`` as the
    underlying engine.  This class requires the `numpy
    <http://www.numpy.org/>`_ library.

    .. warning::
        This interface is still in beta, and may change in the future.

    :param TreeSequence tree_sequence: The tree sequence mutations we are
        interested in.
    """

    def tree_stat_vector(self, sample_sets, weight_fun, windows=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a
        function whose argument is a list of integers of the same length as
        sample_sets that returns a list of numbers.  A branch in a tree is
        weighted by weight_fun(x) + weight_fun(n-x), where x[i] is the number
        of samples in sample_sets[i] below that branch, and n[i]-x[i] is the
        number *not* below that branch.  This finds the sum of this weight for
        all branches in each tree, and averages this across the tree sequence,
        weighted by genomic length.

        It does this separately for each window [windows[i], windows[i+1]) and
        returns the values in a list.  Note that windows cannot be overlapping,
        but overlapping windows can be achieved by (a) computing staistics on a
        small window size and (b) averaging neighboring windows, by additivity
        of the statistics.
        '''
        if windows is None:
            windows = (0, self.tree_sequence.sequence_length)
        for U in sample_sets:
            if ((not isinstance(U, list)) or
               len(U) != len(set(U))):
                raise ValueError(
                    "elements of sample_sets must be lists without repeated elements.")
            if len(U) == 0:
                raise ValueError("elements of sample_sets cannot be empty.")
            for u in U:
                if not self.tree_sequence.node(u).is_sample():
                    raise ValueError("Not all elements of sample_sets are samples.")
        num_windows = len(windows) - 1
        if windows[0] != 0.0:
            raise ValueError(
                "Windows must start at the start of the sequence (at 0.0).")
        if windows[-1] != self.tree_sequence.sequence_length:
            raise ValueError("Windows must extend to the end of the sequence.")
        for k in range(num_windows):
            if windows[k + 1] <= windows[k]:
                raise ValueError("Windows must be increasing.")
        # below we actually just keep track of x, not (x,xbar), so here's the
        #   weighting function we actually use of just x:
        num_sample_sets = len(sample_sets)

        # this how we apply the weight function to both below the branch and
        # above it
        n = [len(x) for x in sample_sets]

        def wfn(x):
            ax = [nn - xx for nn, xx in zip(n, x)]
            return [a + b for a, b in zip(weight_fun(x), weight_fun(ax))]

        # initialize
        n_out = len(wfn([0 for a in range(num_sample_sets)]))

        S = [[0.0 for j in range(n_out)] for _ in range(num_windows)]
        L = [0.0 for j in range(n_out)]
        N = self.tree_sequence.num_nodes
        X = [[int(u in a) for a in sample_sets] for u in range(N)]
        # we will essentially construct the tree
        pi = [-1 for j in range(N)]
        node_time = [self.tree_sequence.node(u).time for u in range(N)]
        # keep track of where we are for the windows
        chrom_pos = 0.0
        # index of *left-hand* end of the current window
        window_num = 0
        for interval, records_out, records_in in self.tree_sequence.edge_diffs():
            length = interval[1] - interval[0]
            for sign, records in ((-1, records_out), (+1, records_in)):
                for edge in records:
                    dx = [0 for k in range(num_sample_sets)]
                    if sign == +1:
                        pi[edge.child] = edge.parent
                    for k in range(num_sample_sets):
                        dx[k] += sign * X[edge.child][k]
                    w = wfn(X[edge.child])
                    dt = (node_time[pi[edge.child]] - node_time[edge.child])
                    for j in range(n_out):
                        L[j] += sign * dt * w[j]
                    if sign == -1:
                        pi[edge.child] = -1
                    old_w = wfn(X[edge.parent])
                    for k in range(num_sample_sets):
                        X[edge.parent][k] += dx[k]
                    if pi[edge.parent] != -1:
                        w = wfn(X[edge.parent])
                        dt = (node_time[pi[edge.parent]] - node_time[edge.parent])
                        for j in range(n_out):
                            L[j] += dt * (w[j]-old_w[j])
                    # propagate change up the tree
                    u = pi[edge.parent]
                    if u != -1:
                        next_u = pi[u]
                        while u != -1:
                            old_w = wfn(X[u])
                            for k in range(num_sample_sets):
                                X[u][k] += dx[k]
                            # need to update X for the root,
                            # but the root does not have a branch length
                            if next_u != -1:
                                w = wfn(X[u])
                                dt = (node_time[pi[u]] - node_time[u])
                                for j in range(n_out):
                                    L[j] += dt*(w[j] - old_w[j])
                            u = next_u
                            next_u = pi[next_u]
            while chrom_pos + length >= windows[window_num + 1]:
                # wrap up the last window
                this_length = windows[window_num + 1] - chrom_pos
                window_length = windows[window_num + 1] - windows[window_num]
                for j in range(n_out):
                    S[window_num][j] += L[j] * this_length
                    S[window_num][j] /= window_length
                length -= this_length
                # start the next
                if window_num < num_windows - 1:
                    window_num += 1
                    chrom_pos = windows[window_num]
                else:
                    # skips the else statement below
                    break
            else:
                for j in range(n_out):
                    S[window_num][j] += L[j] * length
                chrom_pos += length
        return S

    def site_frequency_spectrum(self, sample_set, windows=None):
        '''
        Computes the expected *derived* (unfolded) site frequency spectrum,
        based on tree lengths, separately in each window.

        :param list sample_set: A list of IDs of samples of length n.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of lists of length n, one for each window, whose kth
            entry gives the total length of any branches in the marginal trees
            over that window that are ancestral to exactly k of the samples,
            divided by the length of the window.
        '''
        if windows is None:
            windows = (0, self.tree_sequence.sequence_length)
        if ((not isinstance(sample_set, list)) or
           len(sample_set) != len(set(sample_set))):
            raise ValueError(
                "elements of sample_sets must be lists without repeated elements.")
        if len(sample_set) == 0:
            raise ValueError("elements of sample_sets cannot be empty.")
        for u in sample_set:
            if not self.tree_sequence.node(u).is_sample():
                raise ValueError("Not all elements of sample_sets are samples.")
        num_windows = len(windows) - 1
        if windows[0] != 0.0:
            raise ValueError(
                "Windows must start at the start of the sequence (at 0.0).")
        if windows[-1] != self.tree_sequence.sequence_length:
            raise ValueError("Windows must extend to the end of the sequence.")
        for k in range(num_windows):
            if windows[k + 1] <= windows[k]:
                raise ValueError("Windows must be increasing.")
        n_out = len(sample_set)
        S = [[0.0 for j in range(n_out)] for _ in range(num_windows)]
        L = [0.0 for j in range(n_out)]
        N = self.tree_sequence.num_nodes
        X = [int(u in sample_set) for u in range(N)]
        # we will essentially construct the tree
        pi = [-1 for j in range(N)]
        node_time = [self.tree_sequence.node(u).time for u in range(N)]
        # keep track of where we are for the windows
        chrom_pos = 0.0
        # index of *left-hand* end of the current window
        window_num = 0
        for interval, records_out, records_in in self.tree_sequence.edge_diffs():
            length = interval[1] - interval[0]
            for sign, records in ((-1, records_out), (+1, records_in)):
                for edge in records:
                    dx = 0
                    if sign == +1:
                        pi[edge.child] = edge.parent
                    dx += sign * X[edge.child]
                    dt = (node_time[pi[edge.child]] - node_time[edge.child])
                    if X[edge.child] > 0:
                        L[X[edge.child] - 1] += sign * dt
                    if sign == -1:
                        pi[edge.child] = -1
                    old_X = X[edge.parent]
                    X[edge.parent] += dx
                    if pi[edge.parent] != -1:
                        dt = (node_time[pi[edge.parent]] - node_time[edge.parent])
                        if X[edge.parent] > 0:
                            L[X[edge.parent] - 1] += dt
                        if old_X > 0:
                            L[old_X - 1] -= dt
                    # propagate change up the tree
                    u = pi[edge.parent]
                    if u != -1:
                        next_u = pi[u]
                        while u != -1:
                            old_X = X[u]
                            X[u] += dx
                            # need to update X for the root,
                            # but the root does not have a branch length
                            if next_u != -1:
                                dt = (node_time[pi[u]] - node_time[u])
                                if X[u] > 0:
                                    L[X[u] - 1] += dt
                                if old_X > 0:
                                    L[old_X - 1] -= dt
                            u = next_u
                            next_u = pi[next_u]
            while chrom_pos + length >= windows[window_num + 1]:
                # wrap up the last window
                this_length = windows[window_num + 1] - chrom_pos
                window_length = windows[window_num + 1] - windows[window_num]
                for j in range(n_out):
                    S[window_num][j] += L[j] * this_length
                    S[window_num][j] /= window_length
                length -= this_length
                # start the next
                if window_num < num_windows - 1:
                    window_num += 1
                    chrom_pos = windows[window_num]
                else:
                    # skips the else statement below
                    break
            else:
                for j in range(n_out):
                    S[window_num][j] += L[j] * length
                chrom_pos += length
        return S


class SiteStatCalculator(GeneralStatCalculator):
    """
    Class for calculating a broad class of single-site statistics.  These are
    all calculated using :meth:``SiteStatCalculator.tree_stat_vector`` as the
    underlying engine.  This class requires the `numpy
    <http://www.numpy.org/>`_ library.

    .. warning::
        This interface is still in beta, and may change in the future.

    :param TreeSequence tree_sequence: The tree sequence mutations we are
        interested in.
    """

    def __init__(self, tree_sequence):
        check_numpy()
        self.tree_sequence = tree_sequence

    def tree_stat_vector(self, sample_sets, weight_fun, windows=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a
        function whose argument is a list of integers of the same length as
        sample_sets that returns a list of numbers.  Each allele is weighted by
        weight_fun(x), where x[i] is the number of samples in sample_sets[i]
        that inherit that allele.  This finds the sum of this weight for all
        polymorphic sites, and divides by the sequence length.

        It does this separately for each window [windows[i], windows[i+1]) and
        returns the values in a list.  Note that windows cannot be overlapping,
        but overlapping windows can be achieved by (a) computing staistics on a
        small window size and (b) averaging neighboring windows, by additivity
        of the statistics.
        '''
        if windows is None:
            windows = (0, self.tree_sequence.sequence_length)
        for U in sample_sets:
            if ((not isinstance(U, list)) or
               len(U) != len(set(U))):
                raise ValueError(
                    "elements of sample_sets must be lists without repeated elements.")
            if len(U) == 0:
                raise ValueError("elements of sample_sets cannot be empty.")
            for u in U:
                if not self.tree_sequence.node(u).is_sample():
                    raise ValueError("Not all elements of sample_sets are samples.")
        num_windows = len(windows) - 1
        if windows[0] != 0.0:
            raise ValueError(
                "Windows must start at the start of the sequence (at 0.0).")
        if windows[-1] != self.tree_sequence.sequence_length:
            raise ValueError("Windows must extend to the end of the sequence.")
        for k in range(num_windows):
            if windows[k + 1] <= windows[k]:
                raise ValueError("Windows must be increasing.")
        num_sample_sets = len(sample_sets)
        num_sites = self.tree_sequence.num_sites
        n = [len(x) for x in sample_sets]
        n_out = len(weight_fun([0 for a in range(num_sample_sets)]))
        # we store the final answers here
        S = [[0.0 for j in range(n_out)] for _ in range(num_windows)]
        if num_sites == 0:
            return S
        N = self.tree_sequence.num_nodes
        # initialize: with no tree, each node is either in a sample set or not
        X = [[int(u in a) for a in sample_sets] for u in range(N)]
        # we will construct the tree here
        pi = [-1 for j in range(N)]
        # keep track of which site we're looking at
        sites = self.tree_sequence.sites()
        ns = 0  # this will record number of sites seen so far
        s = next(sites)
        # index of *left-hand* end of the current window
        window_num = 0
        while s.position > windows[window_num + 1]:
            window_num += 1
        for interval, records_out, records_in in self.tree_sequence.edge_diffs():
            # if we've done all the sites then stop
            if ns == num_sites:
                break
            # update the tree
            for sign, records in ((-1, records_out), (+1, records_in)):
                for edge in records:
                    dx = [0 for k in range(num_sample_sets)]
                    if sign == +1:
                        pi[edge.child] = edge.parent
                    for k in range(num_sample_sets):
                        dx[k] += sign * X[edge.child][k]
                    if sign == -1:
                        pi[edge.child] = -1
                    for k in range(num_sample_sets):
                        X[edge.parent][k] += dx[k]
                    # propagate change up the tree
                    u = pi[edge.parent]
                    if u != -1:
                        next_u = pi[u]
                        while u != -1:
                            for k in range(num_sample_sets):
                                X[u][k] += dx[k]
                            u = next_u
                            next_u = pi[next_u]
            # loop over sites in this tree
            while s.position < interval[1]:
                if s.position > windows[window_num + 1]:
                    # finalize this window and move to the next
                    window_length = windows[window_num + 1] - windows[window_num]
                    for j in range(n_out):
                        S[window_num][j] /= window_length
                    # may need to advance through empty windows
                    while s.position > windows[window_num + 1]:
                        window_num += 1
                nm = len(s.mutations)
                if nm > 0:
                    U = {s.ancestral_state: list(n)}
                    for mut in s.mutations:
                        if mut.derived_state not in U:
                            U[mut.derived_state] = [0 for _ in range(num_sample_sets)]
                        for k in range(num_sample_sets):
                            U[mut.derived_state][k] += X[mut.node][k]
                        parent_state = get_derived_state(s, mut.parent)
                        if parent_state not in U:
                            U[parent_state] = [0 for _ in range(num_sample_sets)]
                        for k in range(num_sample_sets):
                            U[parent_state][k] -= X[mut.node][k]
                    for a in U:
                        w = weight_fun(U[a])
                        for j in range(n_out):
                            S[window_num][j] += w[j]
                ns += 1
                if ns == num_sites:
                    break
                s = next(sites)
        # wrap up the final window
        window_length = windows[window_num + 1] - windows[window_num]
        for j in range(n_out):
            S[window_num][j] /= window_length
        return S

    def site_frequency_spectrum(self, sample_set, windows=None):
        '''
        Computes the folded site frequency spectrum in sample_set,
        independently in windows.

        :param list sample_set: A list of IDs of samples of length n.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of lists of length n, one for each window, whose kth
            entry gives the number of mutations in that window at which a mutation
            is seen by exactly k of the samples, divided by the window length.
        '''
        if windows is None:
            windows = (0, self.tree_sequence.sequence_length)
        if ((not isinstance(sample_set, list)) or
           len(sample_set) != len(set(sample_set))):
            raise ValueError(
                "sample_set must not contain repeated elements.")
        if len(sample_set) == 0:
            raise ValueError("sample_set cannot be empty.")
        for u in sample_set:
            if not self.tree_sequence.node(u).is_sample():
                raise ValueError("Not all elements of sample_set are samples.")
        num_windows = len(windows) - 1
        if windows[0] != 0.0:
            raise ValueError(
                "Windows must start at the start of the sequence (at 0.0).")
        if windows[-1] != self.tree_sequence.sequence_length:
            raise ValueError("Windows must extend to the end of the sequence.")
        for k in range(num_windows):
            if windows[k + 1] <= windows[k]:
                raise ValueError("Windows must be increasing.")
        num_sites = self.tree_sequence.num_sites
        n = len(sample_set)
        n_out = n
        # we store the final answers here
        S = [[0.0 for j in range(n_out)] for _ in range(num_windows)]
        if num_sites == 0:
            return S
        N = self.tree_sequence.num_nodes
        # initialize: with no tree, each node is either in a sample set or not
        X = [int(u in sample_set) for u in range(N)]
        # we will construct the tree here
        pi = [-1 for j in range(N)]
        # keep track of which site we're looking at
        sites = self.tree_sequence.sites()
        ns = 0  # this will record number of sites seen so far
        s = next(sites)
        # index of *left-hand* end of the current window
        window_num = 0
        while s.position > windows[window_num + 1]:
            window_num += 1
        for interval, records_out, records_in in self.tree_sequence.edge_diffs():
            # if we've done all the sites then stop
            if ns == num_sites:
                break
            # update the tree
            for sign, records in ((-1, records_out), (+1, records_in)):
                for edge in records:
                    dx = 0
                    if sign == +1:
                        pi[edge.child] = edge.parent
                    dx += sign * X[edge.child]
                    if sign == -1:
                        pi[edge.child] = -1
                    X[edge.parent] += dx
                    # propagate change up the tree
                    u = pi[edge.parent]
                    if u != -1:
                        next_u = pi[u]
                        while u != -1:
                            X[u] += dx
                            u = next_u
                            next_u = pi[next_u]
            # loop over sites in this tree
            while s.position < interval[1]:
                if s.position > windows[window_num + 1]:
                    # finalize this window and move to the next
                    window_length = windows[window_num + 1] - windows[window_num]
                    for j in range(n_out):
                        S[window_num][j] /= window_length
                    # may need to advance through empty windows
                    while s.position > windows[window_num + 1]:
                        window_num += 1
                nm = len(s.mutations)
                if nm > 0:
                    U = {s.ancestral_state: n}
                    for mut in s.mutations:
                        if mut.derived_state not in U:
                            U[mut.derived_state] = 0
                        U[mut.derived_state] += X[mut.node]
                        parent_state = get_derived_state(s, mut.parent)
                        if parent_state not in U:
                            U[parent_state] = 0
                        U[parent_state] -= X[mut.node]
                    for a in U:
                        if U[a] > 0:
                            S[window_num][U[a] - 1] += 1.0
                ns += 1
                if ns == num_sites:
                    break
                s = next(sites)
        # wrap up the final window
        window_length = windows[window_num + 1] - windows[window_num]
        for j in range(n_out):
            S[window_num][j] /= window_length
        return S


def get_derived_state(site, mut_id):
    """
    Find the derived state of the mutation with id `mut_id` at site `site`.
    """
    if mut_id == -1:
        state = site.ancestral_state
    else:
        for m in site.mutations:
            if m.id == mut_id:
                state = m.derived_state
    return state
