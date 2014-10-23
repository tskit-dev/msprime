"""
Test cases for the high level interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import os
import random
import unittest
import tempfile

import msprime

class MsprimeTestCase(unittest.TestCase):
    """
    Superclass of all tests msprime simulator test cases.
    """
    def verify_tree(self, n, pi, tau):
        """
        Verifies that the specified tree is a consistent coalescent history
        for a sample of size n.
        """
        self.assertEqual(len(pi), 2 * n)
        self.assertEqual(len(tau), 2 * n)
        # leading value must be zero.
        self.assertEqual(pi[0], 0)
        self.assertEqual(tau[0], 0)
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

class TestSingleLocusSimulation(MsprimeTestCase):
    """
    Tests on the single locus simulations.
    """
    def test_simple_cases(self):
        for n in range(2, 10):
            pi, tau = msprime.simulate_tree(n)
            self.verify_tree(n, pi, tau)
        for n in [11, 13, 19, 101]:
            pi, tau = msprime.simulate_tree(n)
            self.verify_tree(n, pi, tau)

    def test_error_cases(self):
        for n in [-100, -1, 0, 1]:
            self.assertRaises(msprime.InputError, msprime.simulate_tree, n)
        for n in ["", None, "2", 2.2, 1e5]:
            self.assertRaises(TypeError, msprime.simulate_tree, n)


class TestMultiLocusSimulation(MsprimeTestCase):
    """
    Tests on the single locus simulations.
    """

    def test_simple_cases(self):
        m = 1
        r = 0.1
        for n in range(2, 10):
            self.verify_trees(n, m, msprime.simulate_trees(n, m, r))
        n = 4
        for m in range(1, 10):
            self.verify_trees(n, m, msprime.simulate_trees(n, m, r))
        m = 100
        for r in [0.001, 0.01, 0.1, 1.0]:
            self.verify_trees(n, m, msprime.simulate_trees(n, m, r))

    def test_error_cases(self):
        def f(n, m, r):
            return [t for t in msprime.simulate_trees(n, m, r)]
        for n in [-100, -1, 0, 1]:
            self.assertRaises(msprime.InputError, f, n, 1, 1.0)
        for n in ["", None, "2", 2.2, 1e5]:
            self.assertRaises(TypeError, f, n, 1, 1.0)

class TestTreeSimulator(MsprimeTestCase):
    """
    Runs tests on the underlying TreeSimulator object.
    """
    def setUp(self):
        fd, tf = tempfile.mkstemp(prefix="msp_test_", suffix=".dat")
        os.close(fd)
        self._treefile = tf

    def tearDown(self):
        os.unlink(self._treefile)

    def test_full_simulation(self):
        n = 10
        m = 100
        r = 0.1
        # TODO add some different n and m values here.
        ts = msprime.TreeSimulator(n, self._treefile)
        # todo verify all the setters.
        # self.assertEqual(ts.get_sample_size(), n)
        ts.set_recombination_rate(r)
        ts.set_num_loci(m)
        self.assertTrue(ts.run())
        tf = msprime.TreeFile(self._treefile, 'u')
        self.assertFalse(tf.issorted())
        tf.sort()
        self.assertTrue(tf.issorted())
        tf.close()
        tf = msprime.TreeFile(self._treefile)
        l = [t for t in tf]
        self.verify_trees(n, m, l)
        # self.assertEqual(len(l), sim.get_num_trees())


