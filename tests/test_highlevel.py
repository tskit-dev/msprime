"""
Test cases for the high level interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import os
import random
import re
import unittest

import msprime

import tests

def sparse_tree_to_newick(pi, tau, precision):
    """
    Converts the specified oriented tree to an ms-compatible Newick tree.
    """
    c = {}
    branch_lengths = {}
    root = 1
    for child, node in pi.items():
        if node == 0:
            root = child
        else:
            if node in c:
                c[node].append(child)
                c[node] = sorted(c[node])
            else:
                c[node] = [child]
            s = "{0:.{1}f}".format(tau[node] - tau[child], precision)
            branch_lengths[child] = s
    return _build_newick(root, root, c, branch_lengths)


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


class HighLevelTestCase(tests.MsprimeTestCase):
    """
    Superclass of tests on the high level interface.
    """


    def verify_tree_sequence(self, sim, tree_sequence):
        """
        Verify that varoius ways of looking at the specified tree sequence
        are equivalent.
        """
        self.verify_breakpoints(tree_sequence)
        n = tree_sequence.get_sample_size()
        m = tree_sequence.get_num_loci()
        self.assertEqual(sim.get_sample_size(), n)
        self.assertEqual(sim.get_num_loci(), m)
        # TODO verify the rest of the metadata.

        # check the records are in the right form.
        last_l = 0
        last_t = 0
        for l, r, children, node, t in tree_sequence.records():
            self.assertGreaterEqual(l, last_l)
            if last_l != l:
                last_t = 0.0
                last_l = l
            self.assertGreaterEqual(t, last_t)
            last_t = t

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
        # Now we create the trees using the simplest algorithm where we
        # insert the root immediately after we build it and remove it
        # before we update.
        pi = {}
        tau = {j:0 for j in range(1, n + 1)}
        for node, children, t in records_in:
            pi[children[0]] = node
            pi[children[1]] = node
            tau[node] = t
        # insert the root
        v = 1
        while v in pi:
            v = pi[v]
        pi[v] = 0
        l_p, pi_p, tau_p = sparse_trees[0]
        self.assertEqual(l_p, l)
        self.assertEqual(pi_p, pi)
        self.assertEqual(tau_p, tau)
        del pi[v]
        j = 1
        s = 0
        for l, records_out, records_in in diffs:
            self.assertGreaterEqual(l, 1)
            self.assertEqual(len(records_out), len(records_in))
            # Update the sparse tree
            for node, children, time in records_out:
                for c in children:
                    del pi[c]
                del tau[node]
            for node, children, time in records_in:
                for c in children:
                    pi[c] = node
                tau[node] = time
            v = 1
            while v in pi:
                v = pi[v]
            pi[v] = 0
            l_p, pi_p, tau_p = sparse_trees[j]
            self.assertEqual(l_p, l)
            self.assertEqual(pi_p, pi)
            self.assertEqual(tau_p, tau)
            del pi[v]
            j += 1

    def verify_breakpoints(self, tree_sequence):
        n = tree_sequence.get_sample_size()
        m = tree_sequence.get_num_loci()
        all_breaks = tree_sequence.diffs(all_breaks=True)
        distinct_trees = tree_sequence.diffs(all_breaks=False)
        x_all = 0
        x_distinct = 0
        for l1, out1, in1 in distinct_trees:
            x_distinct += l1
            l2, out2, in2 = next(all_breaks)
            x_all += l2
            self.assertEqual(out1, out2)
            self.assertEqual(in1, in2)
            while x_distinct != x_all:
                l2, out2, in2 = next(all_breaks)
                self.assertEqual(0, len(out2))
                self.assertEqual(0, len(in2))
                x_all += l2
        total = 0
        num_breaks = 0
        breakpoints = tree_sequence.get_breakpoints()
        for l, records_out, records_in in tree_sequence.diffs(all_breaks=True):
            num_breaks += 1
            total += l
        self.assertEqual(num_breaks, len(breakpoints) - 1)
        self.assertEqual(num_breaks, len(breakpoints) - 1)
        self.assertEqual(total,  m)
        total = 0
        num_breaks = 0
        for l, records_out, records_in in tree_sequence.diffs(all_breaks=False):
            total += l
            num_breaks += 1
        self.assertEqual(total,  m)
        self.assertLessEqual(num_breaks, len(breakpoints) - 1)



class TestSingleLocusSimulation(HighLevelTestCase):
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



class TestMultiLocusSimulation(HighLevelTestCase):
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
        tree_sequence = sim.run()
        self.assertEqual(tree_sequence.get_sample_size(), n)
        self.assertEqual(tree_sequence.get_num_loci(), m)
        # TODO reenable test parameters and environment
        # self.verify_parameters(sim, tree_sequence)
        # self.verify_environment(tree_sequence)
        self.verify_tree_sequence(sim, tree_sequence)
        # TODO save the tree_sequence to a file and verify equality
        # between the two.

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
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
        tree_sequence = ts.run()
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
    def verify_trees(self, tree_sequence):
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
        old_trees = [(l, sparse_tree_to_newick(pi, tau, precision))
                for l, pi, tau in tree_sequence.sparse_trees()]
        new_trees = list(tree_sequence.newick_trees(precision))
        self.assertEqual(len(new_trees), len(old_trees))
        for (l1, t1), (l2, t2) in zip(new_trees, old_trees):
            self.assertEqual(l1, l2)
            self.assertEqual(strip_tree(t1), strip_tree(t2))

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


