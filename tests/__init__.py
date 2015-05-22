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

class PythonTreeSequence(object):
    """
    A python implementation of the TreeDiffIterator algorithm.
    """
    def __init__(self, tree_sequence, breakpoints=None):
        self._tree_sequence = tree_sequence
        self._sample_size = tree_sequence.get_sample_size()
        self._breakpoints = breakpoints

    def records(self):
        for j in range(self._tree_sequence.get_num_records()):
            yield self._tree_sequence.get_record(j)

    def _diffs(self):
        n = self._sample_size
        left = 0
        used_records = collections.defaultdict(list)
        records_in = []
        for l, r, node, children, t in self.records():
            if l != left:
                yield l - left, used_records[left], records_in
                del used_records[left]
                records_in = []
                left = l
            used_records[r].append((node, children, t))
            records_in.append((node, children, t))
        yield r - left, used_records[left], records_in

    def _diffs_with_breaks(self):
        k = 1
        x = 0
        b = self._breakpoints
        for length, records_out, records_in in self._diffs():
            x += length
            yield b[k] - b[k - 1], records_out, records_in
            while self._breakpoints[k] != x:
                k += 1
                yield b[k] - b[k - 1], [], []
            k += 1

    def diffs(self, all_breaks=False):
        if all_breaks:
            return self._diffs_with_breaks()
        else:
            return self._diffs()

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
        # verify the root is equal for all leaves
        root = 1
        while pi[root] != 0:
            root = pi[root]
        for j in range(1, n + 1):
            k = j
            while pi[k] != 0:
                k = pi[k]
            if k != root:
                print("ERROR!!")
                print(pi)
                print(j)
            self.assertEqual(k, root)
        self.assertIn(root, tau)
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


