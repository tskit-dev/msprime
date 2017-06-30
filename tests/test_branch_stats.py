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
Test cases for branch length statistic computation.
"""
from __future__ import print_function
from __future__ import division


import unittest
import random

import six

import msprime

##
# This tests implementation of the algorithm described in branch-lengths-methods.md
##


def path_length(tr, x, y):
    L = 0
    mrca = tr.mrca(x, y)
    for u in x, y:
        while u != mrca:
            L += tr.branch_length(u)
            u = tr.parent(u)
    return L


def branch_length_diversity(ts, X, Y, begin=0.0, end=None):
    '''
    Computes average pairwise diversity between a random choice from x
    and a random choice from y over the window specified.
    '''
    if end is None:
        end = ts.sequence_length
    S = 0
    for tr in ts.trees():
        if tr.interval[1] < begin:
            continue
        if tr.interval[0] > end:
            break
        SS = 0
        for x in X:
            for y in Y:
                SS += path_length(tr, x, y)
        S += SS*(min(end, tr.interval[1]) - max(begin, tr.interval[0]))
    return S/((end-begin)*len(X)*len(Y))


def branch_length_diversity_window(ts, X, Y, windows):
    out = [branch_length_diversity(ts, X, Y, windows[k], windows[k+1])
           for k in range(len(windows)-1)]
    return out


def branch_length_Y(ts, x, y, z, begin=0.0, end=None):
    if end is None:
        end = ts.sequence_length
    S = 0
    for tr in ts.trees():
        if tr.interval[1] < begin:
            continue
        if tr.interval[0] > end:
            break
        xy_mrca = tr.mrca(x, y)
        xz_mrca = tr.mrca(x, z)
        yz_mrca = tr.mrca(y, z)
        if xy_mrca == xz_mrca:
            #   /\
            #  / /\
            # x y  z
            S += path_length(tr, x, yz_mrca)*tr.length
        elif xy_mrca == yz_mrca:
            #   /\
            #  / /\
            # y x  z
            S += path_length(tr, x, xz_mrca)*tr.length
        elif xz_mrca == yz_mrca:
            #   /\
            #  / /\
            # z x  y
            S += path_length(tr, x, xy_mrca)*tr.length
    return S/(end-begin)


def branch_length_f4(ts, A, B, C, D):
    for U in A, B, C, D:
        if max([U.count(x) for x in set(U)]) > 1:
            raise ValueError("A,B,C, and D cannot contain repeated elements.")
    S = 0
    for tr in ts.trees():
        SS = 0
        for a in A:
            for b in B:
                for c in C:
                    for d in D:
                        SS += path_length(tr, tr.mrca(a, c), tr.mrca(b, d))
                        SS -= path_length(tr, tr.mrca(a, d), tr.mrca(b, c))
        S += SS*tr.length
    return S/(ts.sequence_length*len(A)*len(B)*len(C)*len(D))


def branch_stats_node_iter(ts, leaf_sets, weight_fun, method='length'):
    '''
    Here leaf_sets is a list of lists of leaves, and weight_fun is a function
    whose argument is a list of integers of the same length as leaf_sets
    that returns a number.  Each branch in a tree is weighted by weight_fun(x),
    where x[i] is the number of leaves in leaf_sets[i] below that
    branch.  This finds the sum of all counted branches for each tree,
    and averages this across the tree sequence ts, weighted by genomic length.

    If method='mutations' instead, then branch lengths will be measured in
    numbers of mutations instead of time.

    This version is inefficient as it iterates over all nodes in each tree.
    '''
    out = branch_stats_vector_node_iter(ts, leaf_sets, lambda x: [weight_fun(x)], method)
    if len(out) > 1:
        raise ValueError("Expecting output of length 1.")
    return out[0]


def branch_stats_vector_node_iter(ts, leaf_sets, weight_fun, method='length'):
    '''
    Here leaf_sets is a list of lists of leaves, and weight_fun is a function
    whose argument is a list of integers of the same length as leaf_sets
    that returns a list of numbers; there will be one output for each element.
    For each value, each branch in a tree is weighted by weight_fun(x),
    where x[i] is the number of leaves in leaf_sets[i] below that
    branch.  This finds the sum of all counted branches for each tree,
    and averages this across the tree sequence ts, weighted by genomic length.

    If method='mutations' instead, then branch lengths will be measured in
    numbers of mutations instead of time.

    This version is inefficient as it iterates over all nodes in each tree.
    '''
    for U in leaf_sets:
        if max([U.count(x) for x in set(U)]) > 1:
            raise ValueError("elements of leaf_sets cannot contain repeated elements.")
    tr_its = [ts.trees(
        tracked_leaves=x,
        leaf_counts=True,
        leaf_lists=True) for x in leaf_sets]
    n_out = len(weight_fun([0 for a in leaf_sets]))
    S = [0.0 for j in range(n_out)]
    for k in range(ts.num_trees):
        trs = [next(x) for x in tr_its]
        root = trs[0].root
        tr_len = trs[0].length
        if method == 'length':
            for node in trs[0].nodes():
                if node != root:
                    x = [tr.num_tracked_leaves(node) for tr in trs]
                    w = weight_fun(x)
                    for j in range(n_out):
                        S[j] += w[j] * trs[0].branch_length(node) * tr_len
        elif method == 'mutations':
            count_nodes = dict(
                [(node, weight_fun([tr.num_tracked_leaves(node) for tr in trs]))
                    for node in trs[0].nodes() if node != root])
            # print(count_nodes)
            # TODO update this to use sites rather than mutations.
            for mut in trs[0].mutations():
                # print(mut)
                for j in range(n_out):
                    # TODO: this is the theoretical method
                    # that assumes we can distinguish recurrent mutations
                    S[j] += count_nodes[mut.node][j]
        else:
            raise(TypeError("Unknown method "+method))
    for j in range(n_out):
        S[j] /= ts.get_sequence_length()
    return S


class BranchStatsTestCase(unittest.TestCase):
    """
    Tests of branch statistic computation.
    """
    random_seed = 123456

    def assertListAlmostEqual(self, x, y):
        for a, b in zip(x, y):
            self.assertAlmostEqual(a, b)

    def check_vectorization(self, ts):
        samples = random.sample(ts.samples(), 3)
        A = [[samples[0]], [samples[1]], [samples[2]]]

        def f(x):
            return [float((x[0] > 0) != (x[1] > 0)),
                    float((x[0] > 0) != (x[2] > 0)),
                    float((x[1] > 0) != (x[2] > 0))]

        self.assertListAlmostEqual(
                branch_stats_vector_node_iter(ts, A, f, method='mutations'),
                [ts.pairwise_diversity(samples=[samples[0], samples[1]]),
                 ts.pairwise_diversity(samples=[samples[0], samples[2]]),
                 ts.pairwise_diversity(samples=[samples[1], samples[2]])])
        self.assertListAlmostEqual(
                branch_stats_vector_node_iter(ts, A, f, method='length'),
                [branch_length_diversity(ts, A[0], A[1]),
                 branch_length_diversity(ts, A[0], A[2]),
                 branch_length_diversity(ts, A[1], A[2])])
        self.assertListAlmostEqual(
                ts.branch_stats_vector(A, f)[0],
                [branch_length_diversity(ts, A[0], A[1]),
                 branch_length_diversity(ts, A[0], A[2]),
                 branch_length_diversity(ts, A[1], A[2])])

    def check_windowization(self, ts):
        samples = random.sample(ts.samples(), 2)
        A_one = [[samples[0]], [samples[1]]]
        A_many = [random.sample(ts.samples(), 2),
                  random.sample(ts.samples(), 2)]
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
            ts.branch_stats_vector(A_one, lambda x: 1.0,
                                   windows=[0.0, 1.0, ts.sequence_length+1.1])

        for A in (A_one, A_many):
            for windows in wins:
                n = [len(a) for a in A]

                def f(x):
                    return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(n[0]*n[1])

                tsdiv_v = ts.branch_stats_vector(A, lambda x: [f(x)], windows)
                tsdiv_vx = [x[0] for x in tsdiv_v]
                tsdiv = ts.branch_stats_windowed(A, f, windows)
                pydiv = branch_length_diversity_window(ts, A[0], A[1], windows)
                self.assertEqual(len(tsdiv), len(windows)-1)
                self.assertListAlmostEqual(tsdiv, pydiv)
                self.assertListEqual(tsdiv, tsdiv_vx)

    def check_pairwise_diversity(self, ts):
        samples = random.sample(ts.samples(), 2)
        A_one = [[samples[0]], [samples[1]]]
        A_many = [random.sample(ts.samples(), 2),
                  random.sample(ts.samples(), 2)]
        for A in (A_one, A_many):
            n = [len(a) for a in A]

            def f(x):
                return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(n[0]*n[1])

            self.assertAlmostEqual(
                    branch_stats_node_iter(ts, A, f, method='length'),
                    branch_length_diversity(ts, A[0], A[1]))
            self.assertAlmostEqual(
                    ts.branch_stats(A, f),
                    branch_length_diversity(ts, A[0], A[1]))

    def check_tmrca_matrix(self, ts):
        A = [random.sample(ts.samples(), 3),
             random.sample(ts.samples(), 2),
             random.sample(ts.samples(), 1)]
        windows = [0.0, ts.sequence_length/2, ts.sequence_length]
        ts_values = ts.mean_pairwise_tmrca(A, windows)
        self.assertListEqual([len(x) for x in ts_values], [6, 6])
        assert(len(A[2]) == 1)
        self.assertListEqual([x[5] for x in ts_values], [0.0, 0.0])
        here_values = [[branch_length_diversity(ts, A[i], A[j], begin=windows[k],
                                                end=windows[k+1])
                        for i in range(len(A)) for j in range(i, len(A))]
                       for k in range(len(windows)-1)]
        for k in range(len(windows)-1):
            self.assertListAlmostEqual(here_values[k], ts_values[k])

    def check_f3_matrix(self, ts):
        A = [random.sample(ts.samples(), 3),
             random.sample(ts.samples(), 2),
             random.sample(ts.samples(), 1)]
        windows = [0.0, ts.sequence_length/2, ts.sequence_length]
        ts_values = ts.mean_pairwise_f3(A, windows)
        self.assertListEqual([len(x) for x in ts_values], [6, 6])
        assert(len(A[2]) == 1)
        self.assertListEqual([x[5] for x in ts_values], [0.0, 0.0])
        here_values = [[branch_length_Y(ts, A[i], A[j], begin=windows[k],
                                                end=windows[k+1])
                        for i in range(len(A)) for j in range(i, len(A))]
                       for k in range(len(windows)-1)]
        for k in range(len(windows)-1):
            self.assertListAlmostEqual(here_values[k], ts_values[k])

    def check_pairwise_diversity_mutations(self, ts):
        samples = random.sample(ts.samples(), 2)
        A = [[samples[0]], [samples[1]]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/float(n[0]*n[1])

        self.assertAlmostEqual(
                branch_stats_node_iter(ts, A, f, method='mutations'),
                ts.pairwise_diversity(samples=samples))

    def check_Y_stat(self, ts):
        samples = random.sample(ts.samples(), 3)
        A = [[samples[0]], samples[1:3]]

        def f(x):
            return float(((x[0] == 1) and (x[1] == 0)) or ((x[0] == 0) and (x[1] == 2)))

        self.assertAlmostEqual(
                ts.branch_stats(A, f),
                branch_length_Y(ts, A[0][0], A[1][0], A[1][1]))
        self.assertAlmostEqual(
                branch_stats_node_iter(ts, A, f, method='length'),
                branch_length_Y(ts, A[0][0], A[1][0], A[1][1]))

    def check_f4_stat(self, ts):
        samples = random.sample(ts.samples(), 4)
        A_zero = [[samples[0]], [samples[0]], [samples[1]], [samples[1]]]
        A_f1 = [[samples[0]], [samples[1]], [samples[0]], [samples[1]]]
        A_one = [[samples[0]], [samples[1]], [samples[2]], [samples[3]]]
        A_many = [random.sample(ts.samples(), 3),
                  random.sample(ts.samples(), 3),
                  random.sample(ts.samples(), 3),
                  random.sample(ts.samples(), 3)]
        for A in (A_zero, A_f1, A_one, A_many):

            def f(x):
                return ((float(x[0])/len(A[0])-float(x[1])/len(A[1]))
                        * (float(x[2])/len(A[2])-float(x[3])/len(A[3])))

            self.assertAlmostEqual(
                    ts.branch_stats(A, f),
                    branch_length_f4(ts, A[0], A[1], A[2], A[3]))
            self.assertAlmostEqual(
                    branch_stats_node_iter(ts, A, f, method='length'),
                    branch_length_f4(ts, A[0], A[1], A[2], A[3]))

    def test_errors(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=10)
        self.assertRaises(AssertionError,
                          ts.mean_pairwise_tmrca, [[0], [1.12]], [0, ts.sequence_length])
        self.assertRaises(AssertionError,
                          ts.mean_pairwise_tmrca, [[0], [11]], [0, ts.sequence_length])
        self.assertRaises(ValueError,
                          ts.mean_pairwise_tmrca, [[0], [1]], [0, ts.sequence_length/2])
        self.assertRaises(ValueError,
                          ts.mean_pairwise_tmrca, [[0], [1]], [ts.sequence_length/2,
                                                          ts.sequence_length])
        self.assertRaises(ValueError,
                          ts.mean_pairwise_tmrca, [[0], [1]], [0.0, 2.0, 1.0,
                                                          ts.sequence_length])

    def test_pairwise_diversity(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_pairwise_diversity(ts)
        self.check_pairwise_diversity_mutations(ts)

    def test_Y_stat(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_Y_stat(ts)

    def test_f4(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_f4_stat(ts)

    def test_vectorization(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_vectorization(ts)

    def test_windowization(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_windowization(ts)

    def test_derived_functions(self):
        '''
        Test implementation of statistics using these functions.
        '''
        ts = msprime.simulate(10, random_seed=self.random_seed, recombination_rate=100)
        self.check_tmrca_matrix(ts)
        self.check_f3_matrix(ts)

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
        true_diversity_01 = 2*(1 * (0.2-0) + 0.5 * (0.8-0.2) + 0.7 * (1.0-0.8))
        true_diversity_02 = 2*(1 * (0.2-0) + 0.4 * (0.8-0.2) + 0.7 * (1.0-0.8))
        true_diversity_12 = 2*(0.5 * (0.2-0) + 0.5 * (0.8-0.2) + 0.5 * (1.0-0.8))
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
        edgesets = six.StringIO("""\
        left    right   parent  children
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
            nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)
        self.check_pairwise_diversity(ts)
        self.check_pairwise_diversity_mutations(ts)
        self.check_Y_stat(ts)
        self.check_vectorization(ts)

        # diversity between 0 and 1
        A = [[0], [1]]

        def f(x):
            return float((x[0] > 0) != (x[1] > 0))

        # branch lengths:
        self.assertAlmostEqual(branch_length_diversity(ts, [0], [1]),
                               true_diversity_01)
        self.assertAlmostEqual(ts.branch_stats(A, f),
                               true_diversity_01)
        self.assertAlmostEqual(branch_stats_node_iter(ts, A, f),
                               true_diversity_01)

        # mean diversity between [0, 1] and [0, 2]:
        true_mean_diversity = (0 + true_diversity_02
                               + true_diversity_01 + true_diversity_12)/4
        A = [[0, 1], [0, 2]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/4.0

        # branch lengths:
        self.assertAlmostEqual(branch_length_diversity(ts, A[0], A[1]),
                               true_mean_diversity)
        self.assertAlmostEqual(ts.branch_stats(A, f),
                               true_mean_diversity)
        self.assertAlmostEqual(branch_stats_node_iter(ts, A, f),
                               true_mean_diversity)

        # Y-statistic for (0/12)
        A = [[0], [1, 2]]

        def f(x):
            return ((x[0] == 1) and (x[1] == 0)) or ((x[0] == 0) and (x[1] == 2))

        # branch lengths:
        true_Y = 0.2*(1 + 0.5) + 0.6*(0.4) + 0.2*(0.7+0.2)
        self.assertAlmostEqual(branch_length_Y(ts, 0, 1, 2), true_Y)
        self.assertAlmostEqual(ts.branch_stats(A, f), true_Y)
        self.assertAlmostEqual(branch_stats_node_iter(ts, A, f), true_Y)

    def test_case_2(self):
        # Here are the trees:
        # t                  |              |              |             |
        #
        # 0       --3--      |     --3--    |     --3--    |    --3--    |    --3--
        #        /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \
        # 1     4   |   5    |   4   |   5  |   4       5  |  4       5  |  4       5
        #       |\ / \ /|    |   |\   \     |   |\     /   |  |\     /   |  |\     /|
        # 2     | 6   7 |    |   | 6   7    |   | 6   7    |  | 6   7    |  | 6   7 |
        #       | |\ /| |    |   |  \  |    |   |  \  |    |  |  \       |  |  \    | ...
        # 3     | | 8 | |    |   |   8 |    |   |   8 |    |  |   8      |  |   8   |
        #       | |/ \| |    |   |  /  |    |   |  /  |    |  |  / \     |  |  / \  |
        # 4     | 9  10 |    |   | 9  10    |   | 9  10    |  | 9  10    |  | 9  10 |
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
        #         |  \    |  |     \    |  |       /  |  |    |  /  |  |    |  /  |
        # 3  ...  |   8   |  |      8   |  |      8   |  |    | 8   |  |    | 8   |
        #         |  / \  |  |     / \  |  |     / \  |  |    |  \  |  |    |  \  |
        # 4       | 9  10 |  |    9  10 |  |    9  10 |  |    9  10 |  |    9  10 |
        #         |    /  |  |   /   /  |  |   /   /  |  |   /   /  |  |   /   /  |
        # 5       0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2
        #
        #         0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0

        # divergence betw 0 and 1
        true_diversity_01 = 2*(0.6*4 + 0.2*2 + 0.2*5)
        # divergence betw 1 and 2
        true_diversity_12 = 2*(0.2*5 + 0.2*2 + 0.3*5 + 0.3*4)
        # divergence betw 0 and 2
        true_diversity_02 = 2*(0.2*5 + 0.2*4 + 0.3*5 + 0.1*4 + 0.2*5)
        # mean divergence between 0, 1 and 0, 2
        true_mean_diversity = (
                0 + true_diversity_02 + true_diversity_01 + true_diversity_12) / 4
        # Y(0;1, 2)
        true_Y = 0.2*4 + 0.2*(4+2) + 0.2*4 + 0.2*2 + 0.2*(5+1)

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
        edgesets = six.StringIO("""\
        left    right   parent  children
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
        ts = msprime.load_text(nodes=nodes, edgesets=edgesets)

        self.check_pairwise_diversity(ts)
        self.check_pairwise_diversity_mutations(ts)
        self.check_Y_stat(ts)
        self.check_vectorization(ts)

        # divergence between 0 and 1
        A = [[0], [1]]

        def f(x):
            return (x[0] > 0) != (x[1] > 0)

        # branch lengths:
        self.assertAlmostEqual(branch_length_diversity(ts, [0], [1]),
                               true_diversity_01)
        self.assertAlmostEqual(ts.branch_stats(A, f),
                               true_diversity_01)
        self.assertAlmostEqual(branch_stats_node_iter(ts, A, f),
                               true_diversity_01)

        # mean divergence between 0, 1 and 0, 2
        A = [[0, 1], [0, 2]]
        n = [len(a) for a in A]

        def f(x):
            return float(x[0]*(n[1]-x[1]) + (n[0]-x[0])*x[1])/4.0

        # branch lengths:
        self.assertAlmostEqual(branch_length_diversity(ts, A[0], A[1]),
                               true_mean_diversity)
        self.assertAlmostEqual(ts.branch_stats(A, f),
                               true_mean_diversity)
        self.assertAlmostEqual(branch_stats_node_iter(ts, A, f),
                               true_mean_diversity)

        # Y-statistic for (0/12)
        A = [[0], [1, 2]]

        def f(x):
            return ((x[0] == 1) and (x[1] == 0)) or ((x[0] == 0) and (x[1] == 2))

        # branch lengths:
        self.assertAlmostEqual(branch_length_Y(ts, 0, 1, 2), true_Y)
        self.assertAlmostEqual(ts.branch_stats(A, f), true_Y)
        self.assertAlmostEqual(branch_stats_node_iter(ts, A, f), true_Y)

    def test_branch_stats_vector_interface(self):
        ts = msprime.simulate(10)

        def f(x):
            return [1.0]

        # Duplicated leaves raise an error
        self.assertRaises(ValueError, ts.branch_stats_vector, [[1, 1]], f)
        self.assertRaises(ValueError, ts.branch_stats_vector, [[1], [2, 2]], f)
        # Make sure the basic call doesn't throw an exception
        ts.branch_stats_vector([[1, 2]], f)
        # Check for bad windows
        for bad_start in [-1, 1, 1e-7]:
            self.assertRaises(
                ValueError, ts.branch_stats_vector, [[1, 2]], f,
                [bad_start, ts.sequence_length])
        for bad_end in [0, ts.sequence_length - 1, ts.sequence_length + 1]:
            self.assertRaises(
                ValueError, ts.branch_stats_vector, [[1, 2]], f,
                [0, bad_end])
        # Windows must be increasing.
        self.assertRaises(
            ValueError, ts.branch_stats_vector, [[1, 2]], f, [0, 1, 1])
