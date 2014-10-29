"""
Test cases for the high level interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import os
import random
import unittest

import msprime
import tests

class TestSingleLocusSimulation(tests.MsprimeTestCase):
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


class TestMultiLocusSimulation(tests.MsprimeTestCase):
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
        for r in [0.001, 0.01]:
            self.verify_trees(n, m, msprime.simulate_trees(n, m, r))

    def test_error_cases(self):
        def f(n, m, r):
            return [t for t in msprime.simulate_trees(n, m, r)]
        for n in [-100, -1, 0, 1]:
            self.assertRaises(msprime.InputError, f, n, 1, 1.0)
        for n in ["", None, "2", 2.2, 1e5]:
            self.assertRaises(TypeError, f, n, 1, 1.0)

class TestTreeSimulator(tests.MsprimeTestCase):
    """
    Runs tests on the underlying TreeSimulator object.
    """
    def verify_simulation(self, n, m, r):
        """
        Verifies a simulation for the specified parameters.
        """
        ts = msprime.TreeSimulator(n, self._treefile)
        # todo verify all the setters.
        # self.assertEqual(ts.get_sample_size(), n)
        ts.set_scaled_recombination_rate(r)
        ts.set_num_loci(m)
        self.assertTrue(ts.run())
        msprime.sort_tree_file(self._treefile)
        tf = msprime.TreeFile(self._treefile)
        l = [t for t in tf]
        self.verify_trees(n, m, l)
        self.assertLessEqual(len(l), ts.get_num_breakpoints())

    def test_random_parameters(self):
        num_random_sims = 1
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            self.verify_simulation(n, m, r)
