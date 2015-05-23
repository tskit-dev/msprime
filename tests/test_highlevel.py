"""
Test cases for the high level interface to msprime.
"""
from __future__ import print_function
from __future__ import division

try:
    # We use the zip as iterator functionality here.
    from future_builtins import zip
except ImportError:
    # This fails for Python 3.x, but that's fine.
    pass

import os
import random
import re
import itertools
import unittest

import msprime

import tests


def sparse_tree_to_newick(st, precision):
    """
    Converts the specified sparse tree to an ms-compatible Newick tree.
    """
    branch_lengths = {}
    stack = [st.root]
    while len(stack) > 0:
        node = stack.pop()
        if st.children[0][node] != 0:
            for c in range(2):
                child = st.children[c][node]
                stack.append(child)
                s = "{0:.{1}f}".format(
                        st.time[node] - st.time[child], precision)
                branch_lengths[child] = s
    return _build_newick(st.root, st.root, st, branch_lengths)

def _build_newick(node, root, tree, branch_lengths):
    c1 = tree.children[0][node]
    c2 = tree.children[1][node]
    if c1 != 0:
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


class HighLevelTestCase(tests.MsprimeTestCase):
    """
    Superclass of tests on the high level interface.
    """
    def verify_sparse_tree(self, st):
        used_nodes = set()
        for j in range(1, st.sample_size + 1):
            self.assertEqual(st.time[j], 0)
            # verify the path to root
            u = j
            times = []
            while st.parent[u] != 0:
                used_nodes.add(u)
                v = st.parent[u]
                times.append(st.time[v])
                self.assertGreaterEqual(st.time[v], 0.0)
                self.assertIn(u, [st.children[c][v] for c in range(2)])
                u = v
            self.assertEqual(u, st.root)
            self.assertEqual(times, sorted(times))
        used_nodes.add(st.root)
        self.assertEqual(len(used_nodes), 2 * st.sample_size - 1)
        # for every entry other than used_nodes we should have an empty row
        for j in range(st.num_nodes):
            if j not in used_nodes:
                self.assertEqual(st.parent[j], 0)
                self.assertEqual(st.time[j], 0)
                for c in range(2):
                    self.assertEqual(st.children[c][j], 0)
        # To a top-down traversal, and make sure we meet all the leaves.
        stack = [st.root]
        leaves = []
        while len(stack) > 0:
            u = stack.pop()
            self.assertNotEqual(u, 0)
            if st.children[0][u] == 0:
                leaves.append(u)
                self.assertEqual(st.children[1][u], 0)
            else:
                for c in range(2):
                    stack.append(st.children[c][u])
        self.assertEqual(sorted(leaves), list(range(1, st.sample_size + 1)))

    def verify_sparse_trees(self, ts):
        pts = tests.PythonTreeSequence(ts.get_ll_tree_sequence())
        iter1 = ts.sparse_trees()
        iter2 = pts.sparse_trees()
        length = 0
        for st1, st2 in zip(iter1, iter2):
            self.assertEqual(st1.sample_size, ts.get_sample_size())
            self.assertEqual(st1, st2)
            self.assertEqual(st1.left, length)
            self.assertGreaterEqual(st1.left, 0)
            self.assertGreater(st1.right, st1.left)
            self.assertLessEqual(st1.right, ts.get_num_loci())
            length += st1.right - st1.left
            self.verify_sparse_tree(st1)
        self.assertEqual(length, ts.get_num_loci())
        self.assertRaises(StopIteration, next, iter1)
        self.assertRaises(StopIteration, next, iter2)




class TestSingleLocusSimulation(HighLevelTestCase):
    """
    Tests on the single locus simulations.
    """
    def test_simple_cases(self):
        for n in range(2, 10):
            st = msprime.generate_tree(n)
            self.verify_sparse_tree(st)
        for n in [11, 13, 19, 101]:
            st = msprime.generate_tree(n)
            self.verify_sparse_tree(st)

    def test_error_cases(self):
        for n in [-100, -1, 0, 1]:
            self.assertRaises(ValueError, msprime.generate_tree, n)
        for n in ["", None, "2", 2.2, 1e5]:
            self.assertRaises(TypeError, msprime.generate_tree, n)

    def test_models(self):
        # Exponential growth of 0 and constant model should be identical.
        m1 = msprime.ExponentialPopulationModel(alpha=0.0, start_time=0.0)
        m2 = msprime.ConstantPopulationModel(size=1.0, start_time=0.0)
        for n in [2, 10, 100]:
            # TODO this _should_ be the same as running with no population
            # models, but it's not. need to investigate.
            st1 = msprime.generate_tree(
                n, random_seed=1, population_models=[m1])
            st2 = msprime.generate_tree(
                n, random_seed=1, population_models=[m2])
            self.assertEqual(st1, st2)
        # TODO add more tests!



class TestMultiLocusSimulation(HighLevelTestCase):
    """
    Tests on the single locus simulations.
    """

    def test_simple_cases(self):
        m = 1
        r = 0.1
        for n in range(2, 10):
            self.verify_sparse_trees(msprime.generate_tree_sequence(n, m, r))
        n = 4
        for m in range(1, 10):
            self.verify_sparse_trees(msprime.generate_tree_sequence(n, m, r))
        m = 100
        for r in [0.001, 0.01]:
            self.verify_sparse_trees(msprime.generate_tree_sequence(n, m, r))

    def test_error_cases(self):
        def f(n, m, r):
            return msprime.generate_tree_sequence(n, m, r)
        for n in [-100, -1, 0, 1]:
            self.assertRaises(ValueError, f, n, 1, 1.0)
        for n in ["", None, "2", 2.2, 1e5]:
            self.assertRaises(TypeError, f, n, 1, 1.0)

class TestTreeSimulator(HighLevelTestCase):
    """
    Runs tests on the underlying TreeSimulator object.
    """
    def verify_parameters(self, sim, tree_sequence):
        parameters = tree_sequence.get_parameters()
        self.assertIsInstance(parameters, dict)
        self.assertEqual(parameters["sample_size"], sim.get_sample_size())
        self.assertEqual(parameters["num_loci"], sim.get_num_loci())
        self.assertEqual(
            parameters["scaled_recombination_rate"],
            sim.get_scaled_recombination_rate())
        self.assertEqual(parameters["random_seed"], sim.get_random_seed())
        models = [m.get_ll_model() for m in sim.get_population_models()]
        self.assertEqual(parameters["population_models"], models)

    def verify_environment(self, tree_sequence):
        environment = tree_sequence.get_environment()
        self.assertIsInstance(environment, dict)
        self.assertGreater(len(environment), 0)

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
        sim.run()
        tree_sequence = sim.get_tree_sequence()
        # TODO verify the records are the same
        self.verify_sparse_trees(tree_sequence)
        # TODO reenable test parameters and environment
        # self.verify_parameters(sim, tree_sequence)
        # self.verify_environment(tree_sequence)
        # TODO save the tree_sequence to a file and verify equality
        # between the two.

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 100)
            r = random.random()
            self.verify_simulation(n, m, r)



class TestHaplotypeGenerator(HighLevelTestCase):
    """
    Tests the haplotype generation code.
    """

    def verify_haplotypes(self, n, haplotype_strings):
        """
        Verify that the specified set of haplotypes are consistent.
        """
        self.assertEqual(len(haplotype_strings), n)
        for h in haplotype_strings:
            self.assertEqual(len(h), len(haplotype_strings[0]))
        # Examine each column; we must have a mixture of 0s and 1s
        for k in range(len(haplotype_strings[0])):
            zeros = 0
            ones = 0
            for j in range(n):
                b = haplotype_strings[j][k]
                zeros += b == '0'
                ones += b == '1'
            self.assertGreater(zeros, 0)
            self.assertEqual(zeros + ones, n)

    def verify_simulation(self, n, m, r, theta):
        """
        Verifies a simulation for the specified parameters.
        """
        ts = msprime.TreeSimulator(n)
        ts.set_scaled_recombination_rate(r)
        ts.set_num_loci(m)
        ts.run()
        tree_sequence = ts.get_tree_sequence()
        tree_sequence.generate_mutations(theta)
        hg = msprime.HaplotypeGenerator(tree_sequence)
        haplotypes = list(hg.haplotypes())
        for h in haplotypes:
            self.assertEqual(len(h), tree_sequence.get_num_mutations())
        self.verify_haplotypes(n, haplotypes)

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            theta = random.uniform(0, 2)
            self.verify_simulation(n, m, r,theta)

class TestNewickConversion(HighLevelTestCase):
    """
    Test the newick tree generation code.
    """
    def verify_trees(self, tree_sequence, breakpoints):
        """
        Verifies that the specified tree is converted to Newick correctly.
        """
        def strip_tree(newick):
            """
            Strips all time information out of the specified newick tree.
            """
            s = newick.replace(":0", "")
            s = s.replace(":1", "")
            return s
        # We set the precision to 0 here to avoid problems that occur when
        # Python and C using different rounding strategies. This allows us
        # to remove the times completely, so we're just comparing the
        # structure of the trees.
        precision = 0
        old_trees = [
            (st.right - st.left, sparse_tree_to_newick(st, precision))
                for st in tree_sequence.sparse_trees()]
        new_trees = list(tree_sequence.newick_trees(precision))
        self.assertEqual(len(new_trees), len(old_trees))
        for (l1, t1), (l2, t2) in zip(new_trees, old_trees):
            self.assertEqual(l1, l2)
            self.assertEqual(strip_tree(t1), strip_tree(t2))
        # TODO test the form of the trees when we're using breakpoints.



    def verify_all_breakpoints(self, tree_sequence, breakpoints):
        """
        Verifies that we get the correct list of trees when we use
        the all_breakpoints option for newick generation.
        """
        trees = list(tree_sequence.newick_trees(2, breakpoints))
        self.assertEqual(len(trees), len(breakpoints) - 1)
        j = 0
        s = 0
        for length, _ in trees:
            self.assertGreater(length, 0)
            self.assertEqual(s, breakpoints[j])
            s += length
            j += 1
        self.assertEqual(s, tree_sequence.get_num_loci())
        pts = tests.PythonTreeSequence(
                tree_sequence.get_ll_tree_sequence(), breakpoints)
        diffs = list(pts.diffs(all_breaks=True))
        self.assertEqual(len(diffs), len(trees))
        for j in range(1, len(diffs)):
            if len(diffs[j][1]) == 0:
                # If the list of diffs is empty, we should have the
                # same tree as the last one.
                self.assertEqual(trees[j][1], trees[j - 1][1])

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
            ts.run()
            tree_sequence = ts.get_tree_sequence()
            breakpoints = ts.get_breakpoints()
            self.verify_trees(tree_sequence, breakpoints)
            self.verify_all_breakpoints(tree_sequence, breakpoints)

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            ts = msprime.TreeSimulator(n)
            ts.set_scaled_recombination_rate(r)
            ts.set_num_loci(m)
            ts.run()
            tree_sequence = ts.get_tree_sequence()
            breakpoints = ts.get_breakpoints()
            self.verify_trees(tree_sequence, breakpoints)
            self.verify_all_breakpoints(tree_sequence, breakpoints)


class TestTreeSequence(HighLevelTestCase):
    """
    Tests for the tree sequence object.
    """
    def get_example_tree_sequences(self):
        for n in [2, 3, 10, 100]:
            for m in [1, 2, 10, 100]:
                for rho in [0, 0.1, 10]:
                    ts = msprime.generate_tree_sequence(n, m, rho)
                    yield ts

    def test_sparse_trees(self):
        for ts in self.get_example_tree_sequences():
            self.verify_sparse_trees(ts)


    def verify_tree_diffs(self, ts):
        pts = tests.PythonTreeSequence(ts.get_ll_tree_sequence())
        iter1 = ts.diffs()
        iter2 = pts.diffs()
        for t1, t2 in zip(iter1, iter2):
            self.assertEqual(t1, t2)
        self.assertRaises(StopIteration, next, iter1)
        self.assertRaises(StopIteration, next, iter2)


    def test_tree_diffs(self):
        for ts in self.get_example_tree_sequences():
            self.verify_tree_diffs(ts)
