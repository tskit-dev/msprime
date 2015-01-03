"""
Common code for the msprime test cases.
"""
from __future__ import print_function
from __future__ import division

import os
import tempfile
import unittest
import random

def setUp():
    # Make random tests reproducible.
    random.seed(2)


class MsprimeTestCase(unittest.TestCase):
    """
    Superclass of all tests msprime simulator test cases.
    """

    def setUp(self):
        fd, tf = tempfile.mkstemp(prefix="msp_test_", suffix=".dat")
        os.close(fd)
        self._treefile = tf

    def tearDown(self):
        os.unlink(self._treefile)

    def verify_tree(self, n, pi, tau):
        """
        Verifies that the specified tree is a consistent coalescent history
        for a sample of size n.
        """
        self.assertEqual(len(pi), 2 * n)
        self.assertEqual(len(tau), 2 * n)
        # leading value must be zero.
        self.assertEqual(pi[0], -1)
        self.assertEqual(tau[0], -1)
        num_children = [0 for j in range(0, 2 * n)]
        for j in range(1, 2 * n):
            num_children[pi[j]] += 1
        self.assertEqual(num_children[0], 1)
        # nodes 1 to n are leaves.
        for j in range(1, n + 1):
            self.assertNotEqual(pi[j], 0)
            self.assertEqual(tau[j], 0)
            self.assertEqual(num_children[j], 0)
        # All non-leaf nodes should be binary with non-zero times.
        for j in range(n + 1, 2 * n):
            self.assertEqual(num_children[j], 2)
            self.assertGreater(tau[j], 0.0)
        # times of non leaves should be distinct and increasing.
        taup = [tau[j] for j in range(n + 1, 2 * n)]
        self.assertEqual(len(set(taup)), len(taup))
        self.assertEqual(taup, sorted(taup))

    def verify_trees(self, n, m, trees):
        """
        Verifies that the specified set of trees is consistent with the specified
        paramters.
        """
        s = 0
        for l, pi, tau in trees:
            self.verify_tree(n, pi, tau)
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
