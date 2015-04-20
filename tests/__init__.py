"""
Common code for the msprime test cases.
"""
from __future__ import print_function
from __future__ import unicode_literals
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
        # verify the root is equal for all leaves
        root = 1
        while root in pi:
            root = pi[root]
        for j in range(1, n + 1):
            k = j
            while k in pi:
                k = pi[k]
            self.assertEqual(k, root)
        self.assertIn(root, tau)
        self.assertEqual(set(list(pi.keys()) + [root]), set(tau.keys()))
        self.assertEqual(len(pi), 2 * n - 2)
        self.assertEqual(len(tau), 2 * n - 1)
        # Zero should not be a node
        self.assertNotIn(0, pi)
        self.assertNotIn(0, tau)
        # 1 to n inclusive should always be nodes
        for j in range(1, n + 1):
            self.assertIn(j, pi)
            self.assertIn(j, tau)
        num_children = collections.defaultdict(int)
        for j in pi.keys():
            num_children[pi[j]] += 1
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
            while k in pi:
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

    def verify_tree_sequence(self, sim, tree_sequence):
        """
        Verify that varoius ways of looking at the specified tree sequence
        are equivalent.
        """
        n = tree_sequence.get_sample_size()
        m = tree_sequence.get_num_loci()
        self.assertEqual(sim.get_sample_size(), n)
        self.assertEqual(sim.get_num_loci(), m)
        # TODO verify the rest of the metadata.

        sparse_trees = [
            (l, dict(pi), dict(tau)) for l, pi, tau in tree_sequence.sparse_trees()]
        self.verify_sparse_trees(n, m, sparse_trees)
        self.assertLess(len(sparse_trees), sim.get_num_breakpoints())
        # Verify the diffs
        diffs = tree_sequence.diffs()
        l, records_out, records_in = next(diffs)
        self.assertGreaterEqual(l, 1)
        self.assertEqual(len(records_out), 0)
        self.assertEqual(len(records_in), n - 1)
        pi = {}
        tau = {j:0 for j in range(1, n + 1)}
        for children, parent, t in records_in:
            pi[children[0]] = parent
            pi[children[1]] = parent
            tau[parent] = t
        l_p, pi_p, tau_p = sparse_trees[0]
        self.assertEqual(l_p, l)
        self.assertEqual(pi_p, pi)
        self.assertEqual(tau_p, tau)
        j = 1
        for l, records_out, records_in in diffs:
            self.assertGreaterEqual(l, 1)
            self.assertEqual(len(records_out), len(records_in))
            for children, parent, time in records_out:
                del pi[children[0]]
                del pi[children[1]]
                del tau[parent]
            for children, parent, time in records_in:
                pi[children[0]] = parent
                pi[children[1]] = parent
                tau[parent] = time
            l_p, pi_p, tau_p = sparse_trees[j]
            self.assertEqual(l_p, l)
            self.assertEqual(pi_p, pi)
            self.assertEqual(tau_p, tau)
            j += 1
