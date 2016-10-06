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
Test cases for stats calculations in msprime.
"""
from __future__ import print_function
from __future__ import division

import unittest
import random

import numpy as np

import msprime


def get_r2_matrix(ts):
    """
    Returns the matrix for the specified tree sequence. This is computed
    via a straightforward Python algorithm.
    """
    n = ts.get_sample_size()
    m = ts.get_num_mutations()
    A = np.zeros((m, m), dtype=float)
    for t1 in ts.trees():
        for mA in t1.mutations():
            A[mA.index, mA.index] = 1
            fA = t1.get_num_leaves(mA.node) / n
            leaves = list(t1.leaves(mA.node))
            for t2 in ts.trees(tracked_leaves=leaves):
                for mB in t2.mutations():
                    if mB.position > mA.position:
                        fB = t2.get_num_leaves(mB.node) / n
                        fAB = t2.get_num_tracked_leaves(mB.node) / n
                        D = fAB - fA * fB
                        r2 = D * D / (fA * fB * (1 - fA) * (1 - fB))
                        A[mA.index, mB.index] = r2
                        A[mB.index, mA.index] = r2
    return A


class TestLdCalculator(unittest.TestCase):
    """
    Tests for the LdCalculator class.
    """

    num_test_mutations = 50

    def verify_matrix(self, ts):
        m = ts.get_num_mutations()
        ldc = msprime.LdCalculator(ts)
        A = ldc.get_r2_matrix()
        self.assertEqual(A.shape, (m, m))
        B = get_r2_matrix(ts)
        self.assertTrue(np.allclose(A, B))

        # Now look at each row in turn, and verify it's the same
        # when we use get_r2 directly.
        for j in range(m):
            a = ldc.get_r2_array(j, direction=msprime.FORWARD)
            b = A[j, j + 1:]
            self.assertEqual(a.shape[0], m - j - 1)
            self.assertEqual(b.shape[0], m - j - 1)
            self.assertTrue(np.allclose(a, b))
            a = ldc.get_r2_array(j, direction=msprime.REVERSE)
            b = A[j, :j]
            self.assertEqual(a.shape[0], j)
            self.assertEqual(b.shape[0], j)
            self.assertTrue(np.allclose(a[::-1], b))

        # Now check every cell in the matrix in turn.
        for j in range(m):
            for k in range(m):
                self.assertAlmostEqual(ldc.get_r2(j, k), A[j, k])

    def verify_max_distance(self, ts):
        """
        Verifies that the max_distance parameter works as expected.
        """
        mutations = list(ts.mutations())
        ldc = msprime.LdCalculator(ts)
        A = ldc.get_r2_matrix()
        j = len(mutations) // 2
        for k in range(j):
            x = mutations[j + k].position - mutations[j].position
            a = ldc.get_r2_array(j, max_distance=x)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j + 1: j + 1 + k], a))
            x = mutations[j].position - mutations[j - k].position
            a = ldc.get_r2_array(j, max_distance=x, direction=msprime.REVERSE)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j - k: j], a[::-1]))
        L = ts.get_sequence_length()
        m = len(mutations)
        a = ldc.get_r2_array(0, max_distance=L)
        self.assertEqual(a.shape[0], m - 1)
        self.assertTrue(np.allclose(A[0, 1:], a))
        a = ldc.get_r2_array(m - 1, max_distance=L, direction=msprime.REVERSE)
        self.assertEqual(a.shape[0], m - 1)
        self.assertTrue(np.allclose(A[m - 1, :-1], a[::-1]))

    def verify_max_mutations(self, ts):
        """
        Verifies that the max mutations parameter works as expected.
        """
        mutations = list(ts.mutations())
        ldc = msprime.LdCalculator(ts)
        A = ldc.get_r2_matrix()
        j = len(mutations) // 2
        for k in range(j):
            a = ldc.get_r2_array(j, max_mutations=k)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j + 1: j + 1 + k], a))
            a = ldc.get_r2_array(j, max_mutations=k, direction=msprime.REVERSE)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j - k: j], a[::-1]))

    def test_single_tree_simulated_mutations(self):
        ts = msprime.simulate(20, mutation_rate=10)
        mutations = random.sample(
            list(ts.mutations()), self.num_test_mutations)
        ts.set_mutations(sorted(mutations))
        self.verify_matrix(ts)
        self.verify_max_distance(ts)

    def test_single_tree_regular_mutations(self):
        ts = msprime.simulate(
            self.num_test_mutations, length=self.num_test_mutations)
        mutations = [(j, j) for j in range(self.num_test_mutations)]
        ts.set_mutations(mutations)
        self.verify_matrix(ts)
        self.verify_max_distance(ts)

    def test_tree_sequence_regular_mutations(self):
        ts = msprime.simulate(
            self.num_test_mutations, recombination_rate=1,
            length=self.num_test_mutations)
        self.assertGreater(ts.get_num_trees(), 10)
        mutations = [(j, j) for j in range(self.num_test_mutations)]
        ts.set_mutations(mutations)
        self.verify_matrix(ts)
        self.verify_max_distance(ts)

    def test_tree_sequence_simulated_mutations(self):
        ts = msprime.simulate(20, mutation_rate=10, recombination_rate=10)
        self.assertGreater(ts.get_num_trees(), 10)
        mutations = random.sample(
            list(ts.mutations()), self.num_test_mutations)
        ts.set_mutations(sorted(mutations))
        self.verify_matrix(ts)
        self.verify_max_distance(ts)
        self.verify_max_mutations(ts)
