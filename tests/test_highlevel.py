"""
Test cases for the high level interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import os
import random
import re
import unittest

# Used to disable the newick and haplotype tests while we're developing
# the new tree format.
from nose.tools import nottest

import msprime

import tests


def sparse_tree_to_newick(pi, tau, precision):
    """
    Converts the specified oriented tree to an ms-compatible Newick tree.
    """
    c = {}
    branch_lengths = {}
    root = 1
    for child, parent in pi.items():
        if parent in c:
            c[parent].append(child)
            c[parent] = sorted(c[parent])
        else:
            c[parent] = [child]
        s = "{0:.{1}f}".format(tau[parent] - tau[child], precision)
        branch_lengths[child] = s
        if parent not in pi:
            root = parent
    return _build_newick(root, root, c, branch_lengths).encode()


def _build_newick(node, root, tree, branch_lengths):
    if node in tree:
        c1, c2 = tree[node]
        s1 = _build_newick(c1, root, tree, branch_lengths)
        s2 = _build_newick(c2, root, tree, branch_lengths)
        if node == root:
            # The root node is treated differently
            s = "({0},{1});".format(s1, s2)
        else:
            s = "({0},{1}):{2}".format(
                s1, s2, branch_lengths[node])
    else:
        # Leaf node
        s = "{0}:{1}".format(node, branch_lengths[node])
    return s


class TestSingleLocusSimulation(tests.MsprimeTestCase):
    """
    Tests on the single locus simulations.
    """
    def test_simple_cases(self):
        for n in range(2, 10):
            pi, tau = msprime.simulate_tree(n)
            self.verify_sparse_tree(n, pi, tau)
        for n in [11, 13, 19, 101]:
            pi, tau = msprime.simulate_tree(n)
            self.verify_sparse_tree(n, pi, tau)

    def test_error_cases(self):
        for n in [-100, -1, 0, 1]:
            self.assertRaises(ValueError, msprime.simulate_tree, n)
        for n in ["", None, "2", 2.2, 1e5]:
            self.assertRaises(TypeError, msprime.simulate_tree, n)

    def test_models(self):
        # Exponential growth of 0 and constant model should be identical.
        m1 = msprime.ExponentialPopulationModel(alpha=0.0, start_time=0.0)
        m2 = msprime.ConstantPopulationModel(size=1.0, start_time=0.0)
        for n in [2, 10, 100]:
            # TODO this _should_ be the same as running with no population
            # models, but it's not. need to investigate.
            pi1, tau1 = msprime.simulate_tree(n, random_seed=1,
                    population_models=[m1])
            pi2, tau2 = msprime.simulate_tree(n, random_seed=1,
                    population_models=[m2])
            self.assertEqual(pi1, pi2)
            self.assertEqual(tau1, tau2)
        # TODO add more tests!



class TestMultiLocusSimulation(tests.MsprimeTestCase):
    """
    Tests on the single locus simulations.
    """

    def test_simple_cases(self):
        m = 1
        r = 0.1
        for n in range(2, 10):
            self.verify_sparse_trees(n, m, msprime.simulate_trees(n, m, r))
        n = 4
        for m in range(1, 10):
            self.verify_sparse_trees(n, m, msprime.simulate_trees(n, m, r))
        m = 100
        for r in [0.001, 0.01]:
            self.verify_sparse_trees(n, m, msprime.simulate_trees(n, m, r))

    def test_error_cases(self):
        def f(n, m, r):
            return [t for t in msprime.simulate_trees(n, m, r)]
        for n in [-100, -1, 0, 1]:
            self.assertRaises(ValueError, f, n, 1, 1.0)
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
        sim = msprime.TreeSimulator(n)
        # TODO verify all the setters.
        self.assertEqual(sim.get_sample_size(), n)
        sim.set_scaled_recombination_rate(r)
        self.assertEqual(sim.get_scaled_recombination_rate(), r)
        sim.set_num_loci(m)
        self.assertEqual(sim.get_num_loci(), m)
        tree_sequence = sim.run()
        self.verify_tree_sequence(sim, tree_sequence)

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            self.verify_simulation(n, m, r)



class TestHaplotypeGenerator(tests.MsprimeTestCase):
    """
    Tests the haplotype generation code.
    """
    def verify_simulation(self, n, m, r, theta):
        """
        Verifies a simulation for the specified parameters.
        """
        ts = msprime.TreeSimulator(n, self._treefile)
        ts.set_scaled_recombination_rate(r)
        ts.set_num_loci(m)
        self.assertTrue(ts.run())
        msprime.sort_tree_file(self._treefile)
        hg = msprime.HaplotypeGenerator(self._treefile, theta)
        self.verify_haplotypes(n, hg.get_haplotypes())

    @nottest
    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            theta = random.uniform(0, 2)
            self.verify_simulation(n, m, r,theta)

class TestNewickConversion(tests.MsprimeTestCase):
    """
    Test the newick tree generation code.
    """
    def verify_trees(self, tree_sequence):
        """
        Verifies that the specified tree is converted to Newick correctly.
        """
        precision = 5
        old_trees = [(l, sparse_tree_to_newick(pi, tau, precision))
                for l, pi, tau in tree_sequence.sparse_trees()]
        new_trees = list(tree_sequence.newick_trees(precision))
        self.assertEqual(len(new_trees), len(old_trees))
        for (l1, t1), (l2, t2) in zip(new_trees, old_trees):
            self.assertEqual(l1, l2)
            self.assertEqual(t1, t2)

    def test_simple_cases(self):
        cases = [
            (2, 1, 0),
            (2, 10, 0.1),
            (4, 10, 0.1),
            (10, 10, 0.1),
            (20, 1, 0),
            (20, 10, 0.1),
            (10, 50, 1.0),
        ]
        for n, m, r in cases:
            ts = msprime.TreeSimulator(n)
            ts.set_random_seed(1)
            ts.set_scaled_recombination_rate(r)
            ts.set_num_loci(m)
            tree_sequence = ts.run()
            self.verify_tree_sequence(ts, tree_sequence)
            self.verify_trees(tree_sequence)

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            ts = msprime.TreeSimulator(n)
            ts.set_scaled_recombination_rate(r)
            ts.set_num_loci(m)
            tree_sequence = ts.run()
            self.verify_tree_sequence(ts, tree_sequence)
            self.verify_trees(tree_sequence)


