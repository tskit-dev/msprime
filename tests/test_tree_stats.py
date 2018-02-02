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
Test cases for generalized statistic computation.
"""
from __future__ import print_function
from __future__ import division


import unittest
import random

import numpy as np
import numpy.testing as nt

import six

import msprime
import tests.tsutil as tsutil


def path_length(tr, x, y):
    L = 0
    mrca = tr.mrca(x, y)
    for u in x, y:
        while u != mrca:
            L += tr.branch_length(u)
            u = tr.parent(u)
    return L


class PythonBranchLengthStatCalculator(object):
    """
    Python implementations of various ("tree") branch-length statistics -
    inefficient but more clear what they are doing.
    """

    def __init__(self, tree_sequence):
        self.tree_sequence = tree_sequence

    def divergence(self, X, Y, begin=0.0, end=None):
        '''
        Computes average pairwise diversity between a random choice from x
        and a random choice from y over the window specified.
        '''
        if end is None:
            end = self.tree_sequence.sequence_length
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            SS = 0
            for x in X:
                for y in Y:
                    SS += path_length(tr, x, y) / 2.0
            S += SS*(min(end, tr.interval[1]) - max(begin, tr.interval[0]))
        return S/((end-begin)*len(X)*len(Y))

    def tree_length_diversity(self, X, Y, begin=0.0, end=None):
        '''
        Computes average pairwise diversity between a random choice from x
        and a random choice from y over the window specified.
        '''
        if end is None:
            end = self.tree_sequence.sequence_length
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            SS = 0
            for x in X:
                for y in Y:
                    SS += path_length(tr, x, y)
            S += SS*(min(end, tr.interval[1]) - max(begin, tr.interval[0]))
        return S/((end-begin)*len(X)*len(Y))

    def Y3(self, X, Y, Z, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            this_length = min(end, tr.interval[1]) - max(begin, tr.interval[0])
            for x in X:
                for y in Y:
                    for z in Z:
                        xy_mrca = tr.mrca(x, y)
                        xz_mrca = tr.mrca(x, z)
                        yz_mrca = tr.mrca(y, z)
                        if xy_mrca == xz_mrca:
                            #   /\
                            #  / /\
                            # x y  z
                            S += path_length(tr, x, yz_mrca) * this_length
                        elif xy_mrca == yz_mrca:
                            #   /\
                            #  / /\
                            # y x  z
                            S += path_length(tr, x, xz_mrca) * this_length
                        elif xz_mrca == yz_mrca:
                            #   /\
                            #  / /\
                            # z x  y
                            S += path_length(tr, x, xy_mrca) * this_length
        return S/((end - begin) * len(X) * len(Y) * len(Z))

    def Y2(self, X, Y, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            this_length = min(end, tr.interval[1]) - max(begin, tr.interval[0])
            for x in X:
                for y in Y:
                    for z in set(Y) - set([y]):
                        xy_mrca = tr.mrca(x, y)
                        xz_mrca = tr.mrca(x, z)
                        yz_mrca = tr.mrca(y, z)
                        if xy_mrca == xz_mrca:
                            #   /\
                            #  / /\
                            # x y  z
                            S += path_length(tr, x, yz_mrca) * this_length
                        elif xy_mrca == yz_mrca:
                            #   /\
                            #  / /\
                            # y x  z
                            S += path_length(tr, x, xz_mrca) * this_length
                        elif xz_mrca == yz_mrca:
                            #   /\
                            #  / /\
                            # z x  y
                            S += path_length(tr, x, xy_mrca) * this_length
        return S/((end - begin) * len(X) * len(Y) * (len(Y)-1))

    def Y1(self, X, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            this_length = min(end, tr.interval[1]) - max(begin, tr.interval[0])
            for x in X:
                for y in set(X) - set([x]):
                    for z in set(X) - set([x, y]):
                        xy_mrca = tr.mrca(x, y)
                        xz_mrca = tr.mrca(x, z)
                        yz_mrca = tr.mrca(y, z)
                        if xy_mrca == xz_mrca:
                            #   /\
                            #  / /\
                            # x y  z
                            S += path_length(tr, x, yz_mrca) * this_length
                        elif xy_mrca == yz_mrca:
                            #   /\
                            #  / /\
                            # y x  z
                            S += path_length(tr, x, xz_mrca) * this_length
                        elif xz_mrca == yz_mrca:
                            #   /\
                            #  / /\
                            # z x  y
                            S += path_length(tr, x, xy_mrca) * this_length
        return S/((end - begin) * len(X) * (len(X)-1) * (len(X)-2))

    def f4(self, A, B, C, D, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        for U in A, B, C, D:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("A,B,C, and D cannot contain repeated elements.")
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            this_length = min(end, tr.interval[1]) - max(begin, tr.interval[0])
            SS = 0
            for a in A:
                for b in B:
                    for c in C:
                        for d in D:
                            SS += path_length(tr, tr.mrca(a, c), tr.mrca(b, d))
                            SS -= path_length(tr, tr.mrca(a, d), tr.mrca(b, c))
            S += SS * this_length
        return S / ((end - begin) * len(A) * len(B) * len(C) * len(D))

    def f3(self, A, B, C, begin=0.0, end=None):
        # this is f4(A,B;A,C) but drawing distinct samples from A
        if end is None:
            end = self.tree_sequence.sequence_length
        assert(len(A) > 1)
        for U in A, B, C:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("A, B and C cannot contain repeated elements.")
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            this_length = min(end, tr.interval[1]) - max(begin, tr.interval[0])
            SS = 0
            for a in A:
                for b in B:
                    for c in set(A) - set([a]):
                        for d in C:
                            SS += path_length(tr, tr.mrca(a, c), tr.mrca(b, d))
                            SS -= path_length(tr, tr.mrca(a, d), tr.mrca(b, c))
            S += SS * this_length
        return S / ((end - begin) * len(A) * (len(A) - 1) * len(B) * len(C))

    def f2(self, A, B, begin=0.0, end=None):
        # this is f4(A,B;A,B) but drawing distinct samples from A and B
        if end is None:
            end = self.tree_sequence.sequence_length
        assert(len(A) > 1)
        for U in A, B:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("A and B cannot contain repeated elements.")
        S = 0
        for tr in self.tree_sequence.trees():
            if tr.interval[1] <= begin:
                continue
            if tr.interval[0] >= end:
                break
            this_length = min(end, tr.interval[1]) - max(begin, tr.interval[0])
            SS = 0
            for a in A:
                for b in B:
                    for c in set(A) - set([a]):
                        for d in set(B) - set([b]):
                            SS += path_length(tr, tr.mrca(a, c), tr.mrca(b, d))
                            SS -= path_length(tr, tr.mrca(a, d), tr.mrca(b, c))
            S += SS * this_length
        return S / ((end - begin) * len(A) * (len(A) - 1) * len(B) * (len(B) - 1))

    def tree_stat(self, sample_sets, weight_fun, begin=0.0, end=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a number.  Each branch in a tree is weighted by weight_fun(x),
        where x[i] is the number of samples in sample_sets[i] below that
        branch.  This finds the sum of all counted branches for each tree,
        and averages this across the tree sequence ts, weighted by genomic length.

        This version is inefficient as it iterates over all nodes in each tree.
        '''
        out = self.tree_stat_vector(sample_sets,
                                    lambda x: [weight_fun(x)],
                                    begin=begin, end=end)
        if len(out) > 1:
            raise ValueError("Expecting output of length 1.")
        return out[0]

    def tree_stat_vector(self, sample_sets, weight_fun, begin=0.0, end=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a list of numbers; there will be one output for each element.
        For each value, each branch in a tree is weighted by weight_fun(x),
        where x[i] is the number of samples in sample_sets[i] below that
        branch.  This finds the sum of all counted branches for each tree,
        and averages this across the tree sequence ts, weighted by genomic length.

        This version is inefficient as it iterates over all nodes in each tree.
        '''
        for U in sample_sets:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("elements of sample_sets",
                                 "cannot contain repeated elements.")
        if end is None:
            end = self.tree_sequence.sequence_length
        tr_its = [
            self.tree_sequence.trees(tracked_samples=x, sample_counts=True)
            for x in sample_sets]
        n = [len(U) for U in sample_sets]
        n_out = len(weight_fun([0 for a in sample_sets]))
        S = [0.0 for j in range(n_out)]
        for k in range(self.tree_sequence.num_trees):
            trs = [next(x) for x in tr_its]
            root = trs[0].root
            tr_len = min(end, trs[0].interval[1]) - max(begin, trs[0].interval[0])
            if tr_len > 0:
                for node in trs[0].nodes():
                    if node != root:
                        x = [tr.num_tracked_samples(node) for tr in trs]
                        nx = [a - b for a, b in zip(n, x)]
                        w = [a + b for a, b in zip(weight_fun(x), weight_fun(nx))]
                        for j in range(n_out):
                            S[j] += w[j] * trs[0].branch_length(node) * tr_len
        for j in range(n_out):
            S[j] /= (end-begin)
        return S

    def site_frequency_spectrum(self, sample_set, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        n_out = len(sample_set)
        S = [0.0 for j in range(n_out)]
        for t in self.tree_sequence.trees(tracked_samples=sample_set,
                                          sample_counts=True):
            root = t.root
            tr_len = min(end, t.interval[1]) - max(begin, t.interval[0])
            if tr_len > 0:
                for node in t.nodes():
                    if node != root:
                        x = t.num_tracked_samples(node)
                        if x > 0:
                            S[x - 1] += t.branch_length(node) * tr_len
        for j in range(n_out):
            S[j] /= (end-begin)
        return S


class PythonSiteStatCalculator(object):
    """
    Python implementations of various single-site statistics -
    inefficient but more clear what they are doing.
    """

    def __init__(self, tree_sequence):
        self.tree_sequence = tree_sequence

    def divergence(self, X, Y, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        haps = list(self.tree_sequence.haplotypes())
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = 0
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                for x in X:
                    for y in Y:
                        if (haps[x][k] != haps[y][k]):
                            # x|y
                            S += 1
        return S/((end - begin) * len(X) * len(Y))

    def Y3(self, X, Y, Z, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        haps = list(self.tree_sequence.haplotypes())
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = 0
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                for x in X:
                    for y in Y:
                        for z in Z:
                            if ((haps[x][k] != haps[y][k])
                               and (haps[x][k] != haps[z][k])):
                                # x|yz
                                S += 1
        return S/((end - begin) * len(X) * len(Y) * len(Z))

    def Y2(self, X, Y, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        haps = list(self.tree_sequence.haplotypes())
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = 0
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                for x in X:
                    for y in Y:
                        for z in set(Y) - set([y]):
                            if ((haps[x][k] != haps[y][k])
                               and (haps[x][k] != haps[z][k])):
                                # x|yz
                                S += 1
        return S/((end - begin) * len(X) * len(Y) * (len(Y) - 1))

    def Y1(self, X, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        haps = list(self.tree_sequence.haplotypes())
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = 0
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                for x in X:
                    for y in set(X) - set([x]):
                        for z in set(X) - set([x, y]):
                            if ((haps[x][k] != haps[y][k])
                               and (haps[x][k] != haps[z][k])):
                                # x|yz
                                S += 1
        return S/((end - begin) * len(X) * (len(X) - 1) * (len(X) - 2))

    def f4(self, A, B, C, D, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        for U in A, B, C, D:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("A,B,C, and D cannot contain repeated elements.")
        haps = list(self.tree_sequence.haplotypes())
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = 0
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                for a in A:
                    for b in B:
                        for c in C:
                            for d in D:
                                if ((haps[a][k] == haps[c][k])
                                   and (haps[a][k] != haps[d][k])
                                   and (haps[a][k] != haps[b][k])):
                                    # ac|bd
                                    S += 1
                                elif ((haps[a][k] == haps[d][k])
                                      and (haps[a][k] != haps[c][k])
                                      and (haps[a][k] != haps[b][k])):
                                    # ad|bc
                                    S -= 1
        return S / ((end - begin) * len(A) * len(B) * len(C) * len(D))

    def f3(self, A, B, C, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        for U in A, B, C:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("A,B,and C cannot contain repeated elements.")
        haps = list(self.tree_sequence.haplotypes())
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = 0
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                for a in A:
                    for b in B:
                        for c in set(A) - set([a]):
                            for d in C:
                                if ((haps[a][k] == haps[c][k])
                                   and (haps[a][k] != haps[d][k])
                                   and (haps[a][k] != haps[b][k])):
                                    # ac|bd
                                    S += 1
                                elif ((haps[a][k] == haps[d][k])
                                      and (haps[a][k] != haps[c][k])
                                      and (haps[a][k] != haps[b][k])):
                                    # ad|bc
                                    S -= 1
        return S / ((end - begin) * len(A) * len(B) * len(C) * (len(A) - 1))

    def f2(self, A, B, begin=0.0, end=None):
        if end is None:
            end = self.tree_sequence.sequence_length
        for U in A, B:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("A,and B cannot contain repeated elements.")
        haps = list(self.tree_sequence.haplotypes())
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = 0
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                for a in A:
                    for b in B:
                        for c in set(A) - set([a]):
                            for d in set(B) - set([b]):
                                if ((haps[a][k] == haps[c][k])
                                   and (haps[a][k] != haps[d][k])
                                   and (haps[a][k] != haps[b][k])):
                                    # ac|bd
                                    S += 1
                                elif ((haps[a][k] == haps[d][k])
                                      and (haps[a][k] != haps[c][k])
                                      and (haps[a][k] != haps[b][k])):
                                    # ad|bc
                                    S -= 1
        return S / ((end - begin) * len(A) * len(B)
                    * (len(A) - 1) * (len(B) - 1))

    def tree_stat_vector(self, sample_sets, weight_fun, begin=0.0, end=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a list of numbers; there will be one output for each element.
        For each value, each allele in a tree is weighted by weight_fun(x), where
        x[i] is the number of samples in sample_sets[i] that inherit that allele.
        This finds the sum of this value for all alleles at all polymorphic sites,
        and across the tree sequence ts, weighted by genomic length.

        This version is inefficient as it works directly with haplotypes.
        '''
        for U in sample_sets:
            if max([U.count(x) for x in set(U)]) > 1:
                raise ValueError("elements of sample_sets",
                                 "cannot contain repeated elements.")
        if end is None:
            end = self.tree_sequence.sequence_length
        haps = list(self.tree_sequence.haplotypes())
        n_out = len(weight_fun([0 for a in sample_sets]))
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = [0.0 for j in range(n_out)]
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                all_g = [haps[j][k] for j in range(self.tree_sequence.num_samples)]
                g = [[haps[j][k] for j in u] for u in sample_sets]
                for a in set(all_g):
                    x = [h.count(a) for h in g]
                    w = weight_fun(x)
                    for j in range(n_out):
                        S[j] += w[j]
        for j in range(n_out):
            S[j] /= (end - begin)
        return S

    def tree_stat(self, sample_sets, weight_fun, begin=0.0, end=None):
        '''
        This provides a non-vectorized interface to `tree_stat_vector()`.
        '''
        out = self.tree_stat_vector(sample_sets, lambda x: [weight_fun(x)],
                                    begin=begin, end=end)
        if len(out) > 1:
            raise ValueError("Expecting output of length 1.")
        return out[0]

    def site_frequency_spectrum(self, sample_set, begin=0.0, end=None):
        '''
        '''
        if end is None:
            end = self.tree_sequence.sequence_length
        haps = list(self.tree_sequence.haplotypes())
        n_out = len(sample_set)
        site_positions = [x.position for x in self.tree_sequence.sites()]
        S = [0.0 for j in range(n_out)]
        for k in range(self.tree_sequence.num_sites):
            if (site_positions[k] >= begin) and (site_positions[k] < end):
                all_g = [haps[j][k] for j in range(self.tree_sequence.num_samples)]
                g = [haps[j][k] for j in sample_set]
                for a in set(all_g):
                    x = g.count(a)
                    if x > 0:
                        S[x - 1] += 1.0
        for j in range(n_out):
            S[j] /= (end - begin)
        return S


def upper_tri_to_matrix(x):
    """
    Given x, a vector of entries of the upper triangle of a matrix
    in row-major order, including the diagonal, return the corresponding matrix.
    """
    # n^2 + n = 2 u => n = (-1 + sqrt(1 + 8*u))/2
    n = int((np.sqrt(1 + 8 * len(x)) - 1)/2.0)
    out = np.ones((n, n))
    k = 0
    for i in range(n):
        for j in range(i, n):
            out[i, j] = out[j, i] = x[k]
            k += 1
    return out


class TestStatsInterface(unittest.TestCase):
    """
    Tests basic stat calculator interface.
    """

    def test_interface(self):
        self.assertRaises(TypeError, msprime.GeneralStatCalculator)
        self.assertRaises(TypeError, msprime.SiteStatCalculator)
        self.assertRaises(TypeError, msprime.BranchLengthStatCalculator)


class GeneralStatsTestCase(unittest.TestCase):
    """
    Tests of statistic computation.  Derived classes should have attributes
    `stat_class` and `py_stat_class`.
    """
    random_seed = 123456

    def assertListAlmostEqual(self, x, y):
        self.assertEqual(len(x), len(y))
        for a, b in zip(x, y):
            self.assertAlmostEqual(a, b)

    def assertArrayEqual(self, x, y):
        nt.assert_equal(x, y)

    def assertArrayAlmostEqual(self, x, y):
        nt.assert_array_almost_equal(x, y)

    def compare_stats(self, ts, tree_fn, sample_sets, index_length,
                      tsc_fn=None, tsc_vector_fn=None):
        """
        Use to compare a tree sequence method tsc_vector_fn to a single-window-based
        implementation tree_fn that takes index_length leaf sets at once.  Pass
        index_length=0 to signal that tsc_fn does not take an 'indices' argument;
        otherwise, gives the length of each of the tuples.

        Here are the arguments these functions will get:
            tree_fn(sample_set[i], ... , sample_set[k], begin=left, end=right)
            tsc_vector_fn(sample_sets, windows, indices)
            ... or tsc_vector_fn(sample_sets, windows)
            tsc_fn(sample_sets, windows)
        """
        assert(len(sample_sets) >= index_length)
        nl = len(sample_sets)
        windows = [k * ts.sequence_length / 20 for k in
                   [0] + sorted(random.sample(range(1, 20), 4)) + [20]]
        indices = [random.sample(range(nl), max(1, index_length)) for _ in range(5)]
        leafset_args = [[sample_sets[i] for i in ii] for ii in indices]
        win_args = [{'begin': windows[i], 'end': windows[i+1]}
                    for i in range(len(windows)-1)]
        tree_vals = [[tree_fn(*a, **b) for a in leafset_args] for b in win_args]
        # flatten if necessary
        if isinstance(tree_vals[0][0], list):
            tree_vals = [[x for a in b for x in a] for b in tree_vals]

        if tsc_vector_fn is not None:
            if index_length > 0:
                tsc_vector_vals = tsc_vector_fn(sample_sets, windows, indices)
            else:
                tsc_vector_vals = tsc_vector_fn([sample_sets[i[0]] for i in indices],
                                                windows)
            # print("vector:")
            # print(tsc_vector_vals)
            # print(tree_vals)
            self.assertEqual(len(tree_vals), len(windows)-1)
            self.assertEqual(len(tsc_vector_vals), len(windows)-1)
            for i in range(len(windows)-1):
                self.assertListAlmostEqual(tsc_vector_vals[i], tree_vals[i])

        if tsc_fn is not None:
            tsc_vals_orig = [tsc_fn(*([ls] + [windows])) for ls in leafset_args]
            tsc_vals = [[x[k][0] for x in tsc_vals_orig] for k in range(len(windows)-1)]
            # print("not:")
            # print(tsc_vals)
            # print(tree_vals)
            self.assertEqual(len(tsc_vals), len(windows)-1)
            for i in range(len(windows)-1):
                self.assertListAlmostEqual(tsc_vals[i], tree_vals[i])

    def compare_sfs(self, ts, tree_fn, sample_sets, tsc_fn):
        """
        """
        for sample_set in sample_sets:
            windows = [k * ts.sequence_length / 20 for k in
                       [0] + sorted(random.sample(range(1, 20), 4)) + [20]]
            win_args = [{'begin': windows[i], 'end': windows[i+1]}
                        for i in range(len(windows)-1)]
            tree_vals = [tree_fn(sample_set, **b) for b in win_args]

            tsc_vals = tsc_fn(sample_set, windows)
            self.assertEqual(len(tsc_vals), len(windows) - 1)
            for i in range(len(windows) - 1):
                self.assertListAlmostEqual(tsc_vals[i], tree_vals[i])

    def check_tree_stat_interface(self, ts):
        samples = list(ts.samples())
        tsc = self.stat_class(ts)

        def wfn(x):
            return [1]

        # empty sample sets will raise an error
        self.assertRaises(ValueError, tsc.tree_stat_vector,
                          samples[0:2] + [], wfn)
        # sample_sets must be lists without repeated elements
        self.assertRaises(ValueError, tsc.tree_stat_vector,
                          samples[0:2], wfn)
        self.assertRaises(ValueError, tsc.tree_stat_vector,
                          [samples[0:2], [samples[2], samples[2]]], wfn)
        # and must all be samples
        self.assertRaises(ValueError, tsc.tree_stat_vector,
                          [samples[0:2], [max(samples)+1]], wfn)
        # windows must start at 0.0, be increasing, and extend to the end
        self.assertRaises(ValueError, tsc.tree_stat_vector,
                          [samples[0:2], samples[2:4]], wfn,
                          [0.1, ts.sequence_length])
        self.assertRaises(ValueError, tsc.tree_stat_vector,
                          [samples[0:2], samples[2:4]], wfn,
                          [0.0, 0.8*ts.sequence_length])
        self.assertRaises(ValueError, tsc.tree_stat_vector,
                          [samples[0:2], samples[2:4]], wfn,
                          [0.0, 0.8*ts.sequence_length, 0.4*ts.sequence_length,
                           ts.sequence_length])

    def check_sfs_interface(self, ts):
        samples = ts.samples()
        tsc = self.stat_class(ts)

        # empty sample sets will raise an error
        self.assertRaises(ValueError, tsc.site_frequency_spectrum, [])
        # sample_sets must be lists without repeated elements
        self.assertRaises(ValueError, tsc.site_frequency_spectrum,
                          [samples[2], samples[2]])
        # and must all be samples
        self.assertRaises(ValueError, tsc.site_frequency_spectrum,
                          [samples[0], max(samples)+1])
        # windows must start at 0.0, be increasing, and extend to the end
        self.assertRaises(ValueError, tsc.site_frequency_spectrum,
                          samples[0:2], [0.1, ts.sequence_length])
        self.assertRaises(ValueError, tsc.site_frequency_spectrum,
                          samples[0:2], [0.0, 0.8*ts.sequence_length])
        self.assertRaises(ValueError, tsc.site_frequency_spectrum,
                          samples[0:2],
                          [0.0, 0.8*ts.sequence_length, 0.4*ts.sequence_length,
                           ts.sequence_length])

    def check_tree_stat_vector(self, ts):
        # test the general tree_stat_vector() machinery
        self.check_tree_stat_interface(ts)
        samples = random.sample(list(ts.samples()), 12)
        A = [[samples[0], samples[1], samples[6]],
             [samples[2], samples[3], samples[7]],
             [samples[4], samples[5], samples[8]],
             [samples[9], samples[10], samples[11]]]
        tsc = self.stat_class(ts)
        py_tsc = self.py_stat_class(ts)

        # a made-up example
        def tsf(sample_sets, windows, indices):
            def f(x):
                return [x[i] + 2.0 * x[j] + 3.5 * x[k] for i, j, k in indices]
            return tsc.tree_stat_vector(sample_sets, weight_fun=f, windows=windows)

        def py_tsf(X, Y, Z, begin, end):
            def f(x):
                return x[0] + 2.0 * x[1] + 3.5 * x[2]
            return py_tsc.tree_stat([X, Y, Z], weight_fun=f,
                                    begin=begin, end=end)

        self.compare_stats(ts, py_tsf, A, 3, tsc_vector_fn=tsf)

    def check_sfs(self, ts):
        # check site frequency spectrum
        self.check_sfs_interface(ts)
        A = [random.sample(list(ts.samples()), 2),
             random.sample(list(ts.samples()), 4),
             random.sample(list(ts.samples()), 8),
             random.sample(list(ts.samples()), 10),
             random.sample(list(ts.samples()), 12)]
        tsc = self.stat_class(ts)
        py_tsc = self.py_stat_class(ts)

        self.compare_sfs(ts, py_tsc.site_frequency_spectrum, A,
                         tsc.site_frequency_spectrum)

    def check_f_interface(self, ts):
        tsc = self.stat_class(ts)
        # sample sets must have at least two samples
        self.assertRaises(ValueError, tsc.f2_vector,
                          [[0, 1], [3]], [0, ts.sequence_length], [(0, 1)])

    def check_f_stats(self, ts):
        self.check_f_interface(ts)
        samples = random.sample(list(ts.samples()), 12)
        A = [[samples[0], samples[1], samples[2]],
             [samples[3], samples[4]],
             [samples[5], samples[6]],
             [samples[7], samples[8]],
             [samples[9], samples[10], samples[11]]]
        tsc = self.stat_class(ts)
        py_tsc = self.py_stat_class(ts)
        self.compare_stats(ts, py_tsc.f2, A, 2,
                           tsc_fn=tsc.f2, tsc_vector_fn=tsc.f2_vector)
        self.compare_stats(ts, py_tsc.f3, A, 3,
                           tsc_fn=tsc.f3, tsc_vector_fn=tsc.f3_vector)
        self.compare_stats(ts, py_tsc.f4, A, 4,
                           tsc_fn=tsc.f4, tsc_vector_fn=tsc.f4_vector)

    def check_Y_stat(self, ts):
        samples = random.sample(list(ts.samples()), 12)
        A = [[samples[0], samples[1], samples[6]],
             [samples[2], samples[3], samples[7]],
             [samples[4], samples[5], samples[8]],
             [samples[9], samples[10], samples[11]]]
        tsc = self.stat_class(ts)
        py_tsc = self.py_stat_class(ts)
        self.compare_stats(ts, py_tsc.Y3, A, 3,
                           tsc_fn=tsc.Y3, tsc_vector_fn=tsc.Y3_vector)
        self.compare_stats(ts, py_tsc.Y2, A, 2,
                           tsc_fn=tsc.Y2, tsc_vector_fn=tsc.Y2_vector)
        self.compare_stats(ts, py_tsc.Y1, A, 0,
                           tsc_vector_fn=tsc.Y1_vector)


class SpecificTreesTestCase(GeneralStatsTestCase):
    seed = 21

    def test_case_1(self):
        # With mutations:
        #
        # 1.0          6
        # 0.7         / \                                    5
        #            /   X                                  / \
        # 0.5       X     4                4               /   4
        #          /     / \              / \             /   X X
        # 0.4     X     X   \            X   3           X   /   \
        #        /     /     X          /   / X         /   /     \
        # 0.0   0     1       2        1   0   2       0   1       2
        #          (0.0, 0.2),        (0.2, 0.8),       (0.8, 1.0)
        #
        branch_true_diversity_01 = 2*(1 * (0.2-0) + 0.5 * (0.8-0.2) + 0.7 * (1.0-0.8))
        branch_true_diversity_02 = 2*(1 * (0.2-0) + 0.4 * (0.8-0.2) + 0.7 * (1.0-0.8))
        branch_true_diversity_12 = 2*(0.5 * (0.2-0) + 0.5 * (0.8-0.2) + 0.5 * (1.0-0.8))
        branch_true_Y = 0.2*(1 + 0.5) + 0.6*(0.4) + 0.2*(0.7+0.2)
        site_true_Y = 3 + 0 + 1

        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           0.4
        4       0           0.5
        5       0           0.7
        6       0           1.0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.2     0.8     3       0,2
        0.0     0.2     4       1,2
        0.2     0.8     4       1,3
        0.8     1.0     4       1,2
        0.8     1.0     5       0,4
        0.0     0.2     6       0,4
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.05        0
        1   0.1         0
        2   0.11        0
        3   0.15        0
        4   0.151       0
        5   0.3         0
        6   0.6         0
        7   0.9         0
        8   0.95        0
        9   0.951       0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state
        0       4       1
        1       0       1
        2       2       1
        3       0       1
        4       1       1
        5       1       1
        6       2       1
        7       0       1
        8       1       1
        9       2       1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        branch_tsc = msprime.BranchLengthStatCalculator(ts)
        py_branch_tsc = PythonBranchLengthStatCalculator(ts)
        site_tsc = msprime.SiteStatCalculator(ts)
        py_site_tsc = PythonSiteStatCalculator(ts)

        # diversity between 0 and 1
        A = [[0], [1]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(2*n[0]*n[1])

        # tree lengths:
        self.assertAlmostEqual(py_branch_tsc.tree_length_diversity([0], [1]),
                               branch_true_diversity_01)
        self.assertAlmostEqual(branch_tsc.tree_stat(A, f),
                               branch_true_diversity_01)
        self.assertAlmostEqual(py_branch_tsc.tree_stat(A, f),
                               branch_true_diversity_01)

        # mean diversity between [0, 1] and [0, 2]:
        branch_true_mean_diversity = (0 + branch_true_diversity_02
                                      + branch_true_diversity_01
                                      + branch_true_diversity_12)/4
        A = [[0, 1], [0, 2]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/8.0

        # tree lengths:
        self.assertAlmostEqual(py_branch_tsc.tree_length_diversity(A[0], A[1]),
                               branch_true_mean_diversity)
        self.assertAlmostEqual(branch_tsc.tree_stat(A, f),
                               branch_true_mean_diversity)
        self.assertAlmostEqual(py_branch_tsc.tree_stat(A, f),
                               branch_true_mean_diversity)

        # Y-statistic for (0/12)
        A = [[0], [1, 2]]

        def f(x):
            return float(((x[0] == 1) and (x[1] == 0))
                         or ((x[0] == 0) and (x[1] == 2)))/2.0

        # tree lengths:
        branch_tsc_Y = branch_tsc.Y3([[0], [1], [2]], [0.0, 1.0])[0][0]
        py_branch_tsc_Y = py_branch_tsc.Y3([0], [1], [2], 0.0, 1.0)
        self.assertAlmostEqual(branch_tsc_Y, branch_true_Y)
        self.assertAlmostEqual(py_branch_tsc_Y, branch_true_Y)
        self.assertAlmostEqual(branch_tsc.tree_stat(A, f), branch_true_Y)
        self.assertAlmostEqual(py_branch_tsc.tree_stat(A, f), branch_true_Y)

        # sites, Y:
        site_tsc_Y = site_tsc.Y3([[0], [1], [2]], [0.0, 1.0])[0][0]
        py_site_tsc_Y = py_site_tsc.Y3([0], [1], [2], 0.0, 1.0)
        self.assertAlmostEqual(site_tsc_Y, site_true_Y)
        self.assertAlmostEqual(py_site_tsc_Y, site_true_Y)
        self.assertAlmostEqual(site_tsc.tree_stat(A, f), site_true_Y)
        self.assertAlmostEqual(py_site_tsc.tree_stat(A, f), site_true_Y)

    def test_case_odds_and_ends(self):
        # Tests having (a) the first site after the first window, and
        # (b) no samples having the ancestral state.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           0.5
        3       0           1.0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.0     0.5     2       0,1
        0.5     1.0     3       0,1
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.65        0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state   parent
        0       0       1               -1
        0       1       2               -1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        site_tsc = msprime.SiteStatCalculator(ts)
        py_site_tsc = PythonSiteStatCalculator(ts)

        # recall that divergence returns the upper triangle
        # with nans on the diag in this case
        py_div = [[np.nan, py_site_tsc.divergence([0], [1], 0.0, 0.5), np.nan],
                  [np.nan, py_site_tsc.divergence([0], [1], 0.5, 1.0), np.nan]]
        div = site_tsc.divergence([[0], [1]], [0.0, 0.5, 1.0])
        self.assertListEqual(py_div[0], div[0])
        self.assertListEqual(py_div[1], div[1])

    def test_case_recurrent_muts(self):
        # With mutations:
        #
        # 1.0          6
        # 0.7         / \                                    5
        #           (0)  \                                  /(6)
        # 0.5      (1)    4                4               /   4
        #          /     / \              / \             /  (7|8)
        # 0.4    (2)   (3)  \           (4)  3           /   /   \
        #        /     /     \          /   /(5)        /   /     \
        # 0.0   0     1       2        1   0   2       0   1       2
        #          (0.0, 0.2),        (0.2, 0.8),       (0.8, 1.0)
        # genotypes:
        #       0     2       0        1   0   1       0   2       3
        site_true_Y = 0 + 1 + 1

        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           0.4
        4       0           0.5
        5       0           0.7
        6       0           1.0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.2     0.8     3       0,2
        0.0     0.2     4       1,2
        0.2     0.8     4       1,3
        0.8     1.0     4       1,2
        0.8     1.0     5       0,4
        0.0     0.2     6       0,4
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.05        0
        1   0.3         0
        2   0.9         0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state   parent
        0       0       1               -1
        0       0       2               0
        0       0       0               1
        0       1       2               -1
        1       1       1               -1
        1       2       1               -1
        2       4       1               -1
        2       1       2               6
        2       2       3               6
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        site_tsc = msprime.SiteStatCalculator(ts)
        py_site_tsc = PythonSiteStatCalculator(ts)

        # Y3:
        site_tsc_Y = site_tsc.Y3([[0], [1], [2]], [0.0, 1.0])[0][0]
        py_site_tsc_Y = py_site_tsc.Y3([0], [1], [2], 0.0, 1.0)
        self.assertAlmostEqual(site_tsc_Y, site_true_Y)
        self.assertAlmostEqual(py_site_tsc_Y, site_true_Y)

    def test_case_2(self):
        # Here are the trees:
        # t                  |              |              |             |
        #
        # 0       --3--      |     --3--    |     --3--    |    --3--    |    --3--
        #        /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \
        # 1     4   |   5    |   4   |   5  |   4       5  |  4       5  |  4       5
        #       |\ / \ /|    |   |\   \     |   |\     /   |  |\     /   |  |\     /|
        # 2     | 6   7 |    |   | 6   7    |   | 6   7    |  | 6   7    |  | 6   7 |
        #       | |\ /| |    |   *  \  |    |   |  \  |    |  |  \       |  |  \    | ...
        # 3     | | 8 | |    |   |   8 *    |   |   8 |    |  |   8      |  |   8   |
        #       | |/ \| |    |   |  /  |    |   |  /  |    |  |  / \     |  |  / \  |
        # 4     | 9  10 |    |   * 9  10    |   | 9  10    |  | 9  10    |  | 9  10 |
        #       |/ \ / \|    |   |  \   \   |   |  \   \   |  |  \   \   |  |  \    |
        # 5     0   1   2    |   0   1   2  |   0   1   2  |  0   1   2  |  0   1   2
        #
        #                    |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5
        # ... continued:
        # t                  |             |             |             |
        #
        # 0         --3--    |    --3--    |    --3--    |    --3--    |    --3--
        #          /     \   |   /     \   |   /     \   |   /     \   |   /  |  \
        # 1       4       5  |  4       5  |  4       5  |  4       5  |  4   |   5
        #         |\     /|  |   \     /|  |   \     /|  |   \     /|  |     /   /|
        # 2       | 6   7 |  |    6   7 |  |    6   7 |  |    6   7 |  |    6   7 |
        #         |  *    *  |     \    |  |       *  |  |    |  /  |  |    |  /  |
        # 3  ...  |   8   |  |      8   |  |      8   |  |    | 8   |  |    | 8   |
        #         |  / \  |  |     / \  |  |     * \  |  |    |  \  |  |    |  \  |
        # 4       | 9  10 |  |    9  10 |  |    9  10 |  |    9  10 |  |    9  10 |
        #         |    /  |  |   /   /  |  |   /   /  |  |   /   /  |  |   /   /  |
        # 5       0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2
        #
        #         0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0
        #
        # Above, subsequent mutations are backmutations.

        # divergence betw 0 and 1
        branch_true_diversity_01 = 2*(0.6*4 + 0.2*2 + 0.2*5)
        # divergence betw 1 and 2
        branch_true_diversity_12 = 2*(0.2*5 + 0.2*2 + 0.3*5 + 0.3*4)
        # divergence betw 0 and 2
        branch_true_diversity_02 = 2*(0.2*5 + 0.2*4 + 0.3*5 + 0.1*4 + 0.2*5)
        # mean divergence between 0, 1 and 0, 2
        branch_true_mean_diversity = (
            0 + branch_true_diversity_02 + branch_true_diversity_01
            + branch_true_diversity_12) / 4
        # Y(0;1, 2)
        branch_true_Y = 0.2*4 + 0.2*(4+2) + 0.2*4 + 0.2*2 + 0.2*(5+1)

        # site stats
        # Y(0;1, 2)
        site_true_Y = 1

        nodes = six.StringIO("""\
        is_sample       time    population
        1       0.000000        0
        1       0.000000        0
        1       0.000000        0
        0       5.000000        0
        0       4.000000        0
        0       4.000000        0
        0       3.000000        0
        0       3.000000        0
        0       2.000000        0
        0       1.000000        0
        0       1.000000        0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.500000        1.000000        10      1
        0.000000        0.400000        10      2
        0.600000        1.000000        9       0
        0.000000        0.500000        9       1
        0.800000        1.000000        8       10
        0.200000        0.800000        8       9,10
        0.000000        0.200000        8       9
        0.700000        1.000000        7       8
        0.000000        0.200000        7       10
        0.800000        1.000000        6       9
        0.000000        0.700000        6       8
        0.400000        1.000000        5       2,7
        0.100000        0.400000        5       7
        0.600000        0.900000        4       6
        0.000000        0.600000        4       0,6
        0.900000        1.000000        3       4,5,6
        0.100000        0.900000        3       4,5
        0.000000        0.100000        3       4,5,7
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.0         0
        1   0.55        0
        2   0.75        0
        3   0.85        0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state   parent
        0       0       1               -1
        0       10      1               -1
        0       0       0               0
        1       8       1               -1
        1       2       1               -1
        2       8       1               -1
        2       9       0               5
        """)
        ts = msprime.load_text(nodes=nodes, edges=edges, sites=sites,
                               mutations=mutations, strict=False)
        branch_tsc = msprime.BranchLengthStatCalculator(ts)
        py_branch_tsc = PythonBranchLengthStatCalculator(ts)
        site_tsc = msprime.SiteStatCalculator(ts)
        py_site_tsc = PythonSiteStatCalculator(ts)

        # divergence between 0 and 1
        A = [[0], [1]]

        def f(x):
            return float((x[0] > 0) != (x[1] > 0))/2.0

        # tree lengths:
        self.assertAlmostEqual(py_branch_tsc.tree_length_diversity([0], [1]),
                               branch_true_diversity_01)
        self.assertAlmostEqual(branch_tsc.tree_stat(A, f),
                               branch_true_diversity_01)
        self.assertAlmostEqual(py_branch_tsc.tree_stat(A, f),
                               branch_true_diversity_01)

        # mean divergence between 0, 1 and 0, 2
        A = [[0, 1], [0, 2]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/8.0

        # tree lengths:
        self.assertAlmostEqual(py_branch_tsc.tree_length_diversity(A[0], A[1]),
                               branch_true_mean_diversity)
        self.assertAlmostEqual(branch_tsc.tree_stat(A, f),
                               branch_true_mean_diversity)
        self.assertAlmostEqual(py_branch_tsc.tree_stat(A, f),
                               branch_true_mean_diversity)

        # Y-statistic for (0/12)
        A = [[0], [1, 2]]

        def f(x):
            return float(((x[0] == 1) and (x[1] == 0))
                         or ((x[0] == 0) and (x[1] == 2)))/2.0

        # tree lengths:
        self.assertAlmostEqual(py_branch_tsc.Y3([0], [1], [2]), branch_true_Y)
        self.assertAlmostEqual(branch_tsc.tree_stat(A, f), branch_true_Y)
        self.assertAlmostEqual(py_branch_tsc.tree_stat(A, f), branch_true_Y)

        # sites:
        site_tsc_Y = site_tsc.Y3([[0], [1], [2]], [0.0, 1.0])[0][0]
        py_site_tsc_Y = py_site_tsc.Y3([0], [1], [2], 0.0, 1.0)
        self.assertAlmostEqual(site_tsc_Y, site_true_Y)
        self.assertAlmostEqual(py_site_tsc_Y, site_true_Y)
        self.assertAlmostEqual(site_tsc.tree_stat(A, f), site_true_Y)
        self.assertAlmostEqual(py_site_tsc.tree_stat(A, f), site_true_Y)

    def test_small_sim(self):
        orig_ts = msprime.simulate(4, random_seed=self.random_seed,
                                   mutation_rate=0.0,
                                   recombination_rate=3.0)
        ts = tsutil.jukes_cantor(orig_ts, num_sites=3, mu=3,
                                 multiple_per_node=True, seed=self.seed)
        branch_tsc = msprime.BranchLengthStatCalculator(ts)
        py_branch_tsc = PythonBranchLengthStatCalculator(ts)
        site_tsc = msprime.SiteStatCalculator(ts)
        py_site_tsc = PythonSiteStatCalculator(ts)

        A = [[0], [1], [2]]
        self.assertAlmostEqual(branch_tsc.Y3(A, [0.0, 1.0])[0][0],
                               py_branch_tsc.Y3(*A))
        self.assertAlmostEqual(site_tsc.Y3(A, [0.0, 1.0])[0][0],
                               py_site_tsc.Y3(*A))

        A = [[0], [1, 2]]
        self.assertAlmostEqual(branch_tsc.Y2(A, [0.0, 1.0])[0][0],
                               py_branch_tsc.Y2(*A))
        self.assertAlmostEqual(site_tsc.Y2(A, [0.0, 1.0])[0][0],
                               py_site_tsc.Y2(*A))


class BranchLengthStatsTestCase(GeneralStatsTestCase):
    """
    Tests of tree statistic computation.
    """
    stat_class = msprime.BranchLengthStatCalculator
    py_stat_class = PythonBranchLengthStatCalculator

    def get_ts(self):
        for N in [12, 15, 20]:
            yield msprime.simulate(N, random_seed=self.random_seed,
                                   recombination_rate=10)

    def check_pairwise_diversity(self, ts):
        samples = random.sample(list(ts.samples()), 2)
        tsc = msprime.BranchLengthStatCalculator(ts)
        py_tsc = PythonBranchLengthStatCalculator(ts)
        A_one = [[samples[0]], [samples[1]]]
        A_many = [random.sample(list(ts.samples()), 2),
                  random.sample(list(ts.samples()), 2)]
        for A in (A_one, A_many):
            n = [len(a) for a in A]

            def f(x):
                return float(x[0]*(n[1]-x[1]))/float(n[0]*n[1])

            self.assertAlmostEqual(
                py_tsc.tree_stat(A, f),
                py_tsc.tree_length_diversity(A[0], A[1]))
            self.assertAlmostEqual(
                tsc.tree_stat(A, f),
                py_tsc.tree_length_diversity(A[0], A[1]))

    def check_divergence_matrix(self, ts):
        # nonoverlapping samples
        samples = random.sample(list(ts.samples()), 6)
        tsc = msprime.BranchLengthStatCalculator(ts)
        py_tsc = PythonBranchLengthStatCalculator(ts)
        A = [samples[0:3], samples[3:5], samples[5:6]]
        windows = [0.0, ts.sequence_length/2, ts.sequence_length]
        ts_values = tsc.divergence(A, windows)
        ts_matrix_values = tsc.divergence_matrix(A, windows)
        self.assertListEqual([len(x) for x in ts_values], [len(samples), len(samples)])
        assert(len(A[2]) == 1)
        self.assertListEqual([x[5] for x in ts_values], [np.nan, np.nan])
        self.assertEqual(len(ts_values), len(ts_matrix_values))
        for w in range(len(ts_values)):
            self.assertArrayEqual(
                ts_matrix_values[w, :, :], upper_tri_to_matrix(ts_values[w]))
        here_values = np.array([[[py_tsc.tree_length_diversity(A[i], A[j],
                                                               begin=windows[k],
                                                               end=windows[k+1])
                                  for i in range(len(A))]
                                 for j in range(len(A))]
                                for k in range(len(windows)-1)])
        for k in range(len(windows)-1):
            for i in range(len(A)):
                for j in range(len(A)):
                    if i == j:
                        if len(A[i]) == 1:
                            here_values[k, i, i] = np.nan
                        else:
                            here_values[k, i, i] /= (len(A[i])-1)/len(A[i])
                    else:
                        here_values[k, j, i]
        for k in range(len(windows)-1):
            self.assertArrayAlmostEqual(here_values[k], ts_matrix_values[k])

    def test_errors(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=10)
        tsc = msprime.BranchLengthStatCalculator(ts)
        self.assertRaises(ValueError,
                          tsc.divergence, [[0], [11]], [0, ts.sequence_length])
        self.assertRaises(ValueError,
                          tsc.divergence, [[0], [1]], [0, ts.sequence_length/2])
        self.assertRaises(ValueError,
                          tsc.divergence, [[0], [1]], [ts.sequence_length/2,
                                                       ts.sequence_length])
        self.assertRaises(ValueError,
                          tsc.divergence, [[0], [1]], [0.0, 2.0, 1.0,
                                                       ts.sequence_length])
        # errors for not enough sample_sets
        self.assertRaises(ValueError,
                          tsc.f4, [[0, 1], [2], [3]], [0, ts.sequence_length])
        self.assertRaises(ValueError,
                          tsc.f3, [[0], [2]], [0, ts.sequence_length])
        self.assertRaises(ValueError,
                          tsc.f2, [[0], [1], [2]], [0, ts.sequence_length])
        # errors if indices aren't of the right length
        self.assertRaises(ValueError,
                          tsc.Y3_vector, [[0], [1], [2]], [0, ts.sequence_length],
                          [[0, 1]])
        self.assertRaises(ValueError,
                          tsc.f4_vector, [[0], [1], [2], [3]], [0, ts.sequence_length],
                          [[0, 1]])
        self.assertRaises(ValueError,
                          tsc.f3_vector, [[0], [1], [2], [3]], [0, ts.sequence_length],
                          [[0, 1]])
        self.assertRaises(ValueError,
                          tsc.f2_vector, [[0], [1], [2], [3]], [0, ts.sequence_length],
                          [[0, 1, 2]])

    def test_windowization(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        samples = random.sample(list(ts.samples()), 2)
        tsc = msprime.BranchLengthStatCalculator(ts)
        py_tsc = PythonBranchLengthStatCalculator(ts)
        A_one = [[samples[0]], [samples[1]]]
        A_many = [random.sample(list(ts.samples()), 2),
                  random.sample(list(ts.samples()), 2)]
        some_breaks = list(set([0.0, ts.sequence_length/2, ts.sequence_length] +
                               random.sample(list(ts.breakpoints()), 5)))
        some_breaks.sort()
        tiny_breaks = ([(k / 4) * list(ts.breakpoints())[1] for k in range(4)] +
                       [ts.sequence_length])
        wins = [[0.0, ts.sequence_length],
                [0.0, ts.sequence_length/2, ts.sequence_length],
                tiny_breaks,
                some_breaks]

        with self.assertRaises(ValueError):
            tsc.tree_stat_vector(A_one, lambda x: 1.0,
                                 windows=[0.0, 1.0, ts.sequence_length+1.1])

        for A in (A_one, A_many):
            for windows in wins:
                n = [len(a) for a in A]

                def f(x):
                    return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(2*n[0]*n[1])

                def g(x):
                    return [f(x)]

                tsdiv_v = tsc.tree_stat_vector(A, g, windows)
                tsdiv_vx = [x[0] for x in tsdiv_v]
                tsdiv = tsc.tree_stat_windowed(A, f, windows)
                pydiv = [py_tsc.tree_length_diversity(A[0], A[1], windows[k],
                                                      windows[k+1])
                         for k in range(len(windows)-1)]
                self.assertEqual(len(tsdiv), len(windows)-1)
                self.assertListAlmostEqual(tsdiv, pydiv)
                self.assertListEqual(tsdiv, tsdiv_vx)

    def test_tree_stat_vector_interface(self):
        ts = msprime.simulate(10)
        tsc = msprime.BranchLengthStatCalculator(ts)

        def f(x):
            return [1.0]

        # Duplicated samples raise an error
        self.assertRaises(ValueError, tsc.tree_stat_vector, [[1, 1]], f)
        self.assertRaises(ValueError, tsc.tree_stat_vector, [[1], [2, 2]], f)
        # Make sure the basic call doesn't throw an exception
        tsc.tree_stat_vector([[1, 2]], f)
        # Check for bad windows
        for bad_start in [-1, 1, 1e-7]:
            self.assertRaises(
                ValueError, tsc.tree_stat_vector, [[1, 2]], f,
                [bad_start, ts.sequence_length])
        for bad_end in [0, ts.sequence_length - 1, ts.sequence_length + 1]:
            self.assertRaises(
                ValueError, tsc.tree_stat_vector, [[1, 2]], f,
                [0, bad_end])
        # Windows must be increasing.
        self.assertRaises(
            ValueError, tsc.tree_stat_vector, [[1, 2]], f, [0, 1, 1])

    def test_sfs_interface(self):
        ts = msprime.simulate(10)
        tsc = msprime.BranchLengthStatCalculator(ts)

        # Duplicated samples raise an error
        self.assertRaises(ValueError, tsc.site_frequency_spectrum, [1, 1])
        self.assertRaises(ValueError, tsc.site_frequency_spectrum, [])
        self.assertRaises(ValueError, tsc.site_frequency_spectrum, [0, 11])
        # Check for bad windows
        for bad_start in [-1, 1, 1e-7]:
            self.assertRaises(
                ValueError, tsc.site_frequency_spectrum, [1, 2],
                [bad_start, ts.sequence_length])
        for bad_end in [0, ts.sequence_length - 1, ts.sequence_length + 1]:
            self.assertRaises(
                ValueError, tsc.site_frequency_spectrum, [1, 2],
                [0, bad_end])
        # Windows must be increasing.
        self.assertRaises(
            ValueError, tsc.site_frequency_spectrum, [1, 2], [0, 1, 1])

    def test_branch_general_stats(self):
        for ts in self.get_ts():
            self.check_tree_stat_vector(ts)

    def test_branch_f_stats(self):
        for ts in self.get_ts():
            self.check_f_stats(ts)

    def test_branch_Y_stats(self):
        for ts in self.get_ts():
            self.check_Y_stat(ts)

    def test_diversity(self):
        for ts in self.get_ts():
            self.check_pairwise_diversity(ts)
            self.check_divergence_matrix(ts)

    def test_branch_sfs(self):
        for ts in self.get_ts():
            self.check_sfs(ts)


class SiteStatsTestCase(GeneralStatsTestCase):
    """
    Tests of site statistic computation.
    """
    stat_class = msprime.SiteStatCalculator
    py_stat_class = PythonSiteStatCalculator
    seed = 23

    def get_ts(self):
        for mut in [0.0, 3.0]:
            yield msprime.simulate(20, random_seed=self.random_seed,
                                   mutation_rate=mut,
                                   recombination_rate=3.0)
        ts = msprime.simulate(20, random_seed=self.random_seed,
                              mutation_rate=0.0,
                              recombination_rate=3.0)
        for mpn in [False, True]:
            for num_sites in [10, 100]:
                mut_ts = tsutil.jukes_cantor(ts, num_sites=num_sites, mu=3,
                                             multiple_per_node=mpn, seed=self.seed)
                yield mut_ts

    def check_pairwise_diversity_mutations(self, ts):
        py_tsc = PythonSiteStatCalculator(ts)
        samples = random.sample(list(ts.samples()), 2)
        A = [[samples[0]], [samples[1]]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(2*n[0]*n[1])

        self.assertAlmostEqual(
            py_tsc.tree_stat(A, f), ts.pairwise_diversity(samples=samples))

    def test_pairwise_diversity(self):
        ts = msprime.simulate(20, random_seed=self.random_seed, recombination_rate=100)
        self.check_pairwise_diversity_mutations(ts)

    def test_site_general_stats(self):
        for ts in self.get_ts():
            self.check_tree_stat_vector(ts)

    def test_site_f_stats(self):
        for ts in self.get_ts():
            self.check_f_stats(ts)

    def test_site_Y_stats(self):
        for ts in self.get_ts():
            self.check_Y_stat(ts)

    def test_site_sfs(self):
        for ts in self.get_ts():
            self.check_sfs(ts)
