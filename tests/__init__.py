"""
Common code for the msprime test cases.
"""
from __future__ import print_function
from __future__ import division

import collections
import os
import random
import tempfile
import unittest

def setUp():
    # Make random tests reproducible.
    random.seed(210)


class MsprimeTestCase(unittest.TestCase):
    """
    Superclass of all tests msprime simulator test cases.
    """

    def verify_sparse_tree(self, n, pi, tau):
        """
        Verifies that the specified sparse_tree is a consistent coalescent history
        for a sample of size n.
        """
        self.assertEqual(set(pi.keys()), set(tau.keys()))
        self.assertEqual(len(pi), 2 * n - 1)
        self.assertEqual(len(tau), 2 * n - 1)
        # Zero should not be a node
        self.assertNotIn(0, pi)
        self.assertNotIn(0, tau)
        # 1 to n inclusive should always be nodes
        for j in range(1, n + 1):
            self.assertIn(j, pi)
            self.assertIn(j, tau)
        num_children = collections.defaultdict(int)
        roots = 0
        for j in pi.keys():
            num_children[pi[j]] += 1
            roots += pi[j] == 0
        if roots != 1:
            print(pi)
        self.assertEqual(roots, 1)
        # nodes 1 to n are leaves.
        for j in range(1, n + 1):
            self.assertNotEqual(pi[j], 0)
            self.assertEqual(tau[j], 0)
            self.assertEqual(num_children[j], 0)
        # All non-leaf nodes should be binary with non-zero times.
        taup = {}
        for j in pi.keys():
            if j > n:
                self.assertEqual(num_children[j], 2)
                self.assertGreater(tau[j], 0.0)
                taup[j] = tau[j]
        # times of non leaves should be distinct
        self.assertEqual(len(set(taup)), len(taup))
        # Times of leaves should be zero, and increasing up the tree
        for j in range(1, n + 1):
            self.assertEqual(tau[j], 0.0)
            last_time = -1
            k = j
            while k != 0:
                self.assertNotEqual(k, pi[k])
                self.assertGreater(tau[k], last_time)
                last_time = tau[k]
                k = pi[k]

    def verify_sparse_trees(self, n, m, trees):
        """
        Verifies that the specified set of sparse trees is consistent with the specified
        paramters.
        """
        s = 0
        for l, pi, tau in trees:
            self.verify_sparse_tree(n, pi, tau)
            self.assertTrue(l > 0)
            s += l
        self.assertEqual(s, m)

    def verify_haplotypes(self, n, haplotypes):
        """
        Verify that the specified set of haplotypes are consistent.
        """
        self.assertEqual(len(haplotypes), n)
        for h in haplotypes:
            self.assertEqual(len(h), len(haplotypes[0]))
        # Examine each column; we must have a mixture of 0s and 1s
        for k in range(len(haplotypes[0])):
            zeros = 0
            ones = 0
            for j in range(n):
                zeros += haplotypes[j][k] == '0'
                ones += haplotypes[j][k] == '1'
            self.assertEqual(zeros + ones, n)
