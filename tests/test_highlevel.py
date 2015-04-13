"""
Test cases for the high level interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import re
import os
import random
import unittest

# Used to disable the newick and haplotype tests while we're developing
# the new tree format.
from nose.tools import nottest

import msprime

import tests

def oriented_tree_to_newick(pi, tau, precision):
    """
    Converts the specified oriented tree to an ms-compatible Newick tree.
    """
    # Build a top-down linked tree using a dict. Each node has a list
    # of its children
    d = {}
    n = len(pi) // 2
    # We also want the branch lengths; time from a node back to its parent
    b = {2 * n - 1: None}
    for j in range(1, n + 1):
        u = j
        d[u] = None
        not_done = True
        while pi[u] != 0 and not_done:
            b[u] = "{0:.{1}f}".format(tau[pi[u]] - tau[u], precision)
            if pi[u] in d:
                # This is the second time we've seen this node
                not_done = False
            else:
                d[pi[u]] = []
            d[pi[u]].append(u)
            # make sure the children are sorted so we can compare trees easily
            if len(d[pi[u]]) == 2:
                d[pi[u]] = sorted(d[pi[u]])
            u = pi[u]
    return _build_newick(2 * n - 1, d, b)

def _build_newick(node, tree, branch_lengths):
    l = branch_lengths[node]
    if tree[node] is not None:
        c1, c2 = tree[node]
        s1 = _build_newick(c1, tree, branch_lengths)
        s2 = _build_newick(c2, tree, branch_lengths)
        if l is None:
            # The root node is treated differently
            s = "({0},{1});".format(s1, s2)
        else:
            s = "({0},{1}):{2}".format(s1, s2, l)
    else:
        s = "{0}:{1}".format(node, l)
    return s

class TestSingleLocusSimulation(tests.MsprimeTestCase):
    """
    Tests on the single locus simulations.
    """
    def test_simple_cases(self):
        for n in range(2, 10):
            pi, tau = msprime.simulate_tree(n)
            self.verify_dense_tree(n, pi, tau)
        for n in [11, 13, 19, 101]:
            pi, tau = msprime.simulate_tree(n)
            self.verify_dense_tree(n, pi, tau)

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
            self.verify_dense_trees(n, m, msprime.simulate_trees(n, m, r))
        n = 4
        for m in range(1, 10):
            self.verify_dense_trees(n, m, msprime.simulate_trees(n, m, r))
        m = 100
        for r in [0.001, 0.01]:
            self.verify_dense_trees(n, m, msprime.simulate_trees(n, m, r))

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
    def verify_simulation(self, n, m, r, squash_records):
        """
        Verifies a simulation for the specified parameters.
        """
        ts = msprime.TreeSimulator(n, self._treefile)
        # todo verify all the setters.
        # self.assertEqual(ts.get_sample_size(), n)
        ts.set_scaled_recombination_rate(r)
        ts.set_squash_records(squash_records)
        ts.set_num_loci(m)
        self.assertTrue(ts.run())
        self.assertEqual(squash_records, ts.get_squash_records())
        msprime.sort_tree_file(self._treefile)
        tf = msprime.TreeFile(self._treefile)
        sparse_trees = [(l, dict(pi), dict(tau)) for l, pi, tau in tf.sparse_trees()]
        self.verify_sparse_trees(n, m, sparse_trees)
        tf = msprime.TreeFile(self._treefile)
        dense_trees = [(l, list(pi), list(tau)) for l, pi, tau in tf.dense_trees()]
        self.verify_dense_trees(n, m, dense_trees)
        self.assertEqual(len(sparse_trees), len(dense_trees))
        for (l_sparse, sparse_pi, sparse_tau), (l_dense, dense_pi, dense_tau) in zip(
                sparse_trees, dense_trees):
            self.assertEqual(l_sparse, l_dense)
            self.assertTreesEqual(n, sparse_pi, sparse_tau, dense_pi, dense_tau)
        # With inline segment merging, the number of trees we get back is
        # not necessarily equal to the number of breakpoints.
        self.assertLessEqual(len(sparse_trees), ts.get_num_breakpoints())
        # TODO implement full tree access that uses the breakpoints so that
        # we can fully recover the set of trees including the adjacent
        # identical trees.
        # If record squashing is on, we won't necessarily get a distinct
        # tree for every recombination breakpoint.
        # if ts.get_squash_records():
        #     self.assertLessEqual(len(sparse_trees), ts.get_num_breakpoints())
        # else:
        #     self.assertEqual(len(sparse_trees), ts.get_num_breakpoints())

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            s = bool(random.randint(0, 1))
            self.verify_simulation(n, m, r, s)



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
    def verify_trees(self, treefile):
        """
        Verifies that the specified tree is converted to Newick correctly.
        """
        def strip_tree(newick):
            """
            Strips all time information out of the specified newick tree.
            """
            return re.sub(":\d", "", newick)
        # We set the precision to 0 here to avoid problems that occur when
        # Python and C using different rounding strategies. This allows us
        # to remove the times completely, so we're just comparing the
        # structure of the trees.
        precision = 0
        old_trees = [(l, oriented_tree_to_newick(pi, tau, precision))
                for l, pi, tau in msprime.TreeFile(treefile).trees()]
        new_trees = list(msprime.TreeFile(treefile).newick_trees(precision))
        self.assertEqual(len(new_trees), len(old_trees))
        for (l1, t1), (l2, t2) in zip(new_trees, old_trees):
            self.assertEqual(l1, l2)
            s1 = strip_tree(t1)
            s2 = strip_tree(t2)
            self.assertEqual(s1, s2)

    @nottest
    def test_simple_cases(self):
        cases = [
            (2, 1, 0),
            (2, 10, 0.1),
            (4, 10, 0.1),
            (10, 10, 0.1),
            (20, 1, 0),
            (20, 10, 0.1),
        ]
        for n, m, r in cases:
            ts = msprime.TreeSimulator(n, self._treefile)
            ts.set_random_seed(1)
            ts.set_scaled_recombination_rate(r)
            ts.set_num_loci(m)
            self.assertTrue(ts.run())
            msprime.sort_tree_file(self._treefile)
            self.verify_trees(self._treefile)

    @nottest
    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            squash_records = bool(random.randint(0, 1))
            r = random.random()
            ts = msprime.TreeSimulator(n, self._treefile)
            ts.set_scaled_recombination_rate(r)
            ts.set_num_loci(m)
            ts.set_squash_records(squash_records)
            self.assertTrue(ts.run())
            msprime.sort_tree_file(self._treefile)
            self.verify_trees(self._treefile)


