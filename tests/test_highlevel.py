#
# Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msprime is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msprime.  If not, see <http://www.gnu.org/licenses/>.
#
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

import math
import os
import random
import tempfile
import unittest
import xml.etree

import msprime
import tests


def get_pairwise_diversity(haplotypes):
    """
    Returns the value of pi for the specified haplotypes.
    """
    # Very simplistic algorithm...
    n = len(haplotypes)
    pi = 0
    for k in range(n):
        for j in range(k):
            for u, v in zip(haplotypes[j], haplotypes[k]):
                pi += u != v
    return 2 * pi / (n * (n - 1))


def sparse_tree_to_newick(st, precision):
    """
    Converts the specified sparse tree to an ms-compatible Newick tree.
    """
    branch_lengths = {}
    root = st.get_root()
    stack = [root]
    while len(stack) > 0:
        node = stack.pop()
        if st.is_internal(node):
            for child in st.get_children(node):
                stack.append(child)
                s = "{0:.{1}f}".format(
                    st.get_time(node) - st.get_time(child), precision)
                branch_lengths[child] = s
    return _build_newick(root, root, st, branch_lengths)


def _build_newick(node, root, tree, branch_lengths):
    c1, c2 = tree.get_children(node)
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


class TestHarmonicNumber(unittest.TestCase):
    """
    Tests for the harmonic number calculation.
    """

    def test_harmonic_number(self):
        def H(n):
            return sum(1 / k for k in range(1, n + 1))
        for n in range(10, 1000, 100):
            self.assertAlmostEqual(msprime.harmonic_number(n), H(n), 1)


class TestMsCommandLine(tests.MsprimeTestCase):
    """
    Tests the output of the get_ms_command_line method.
    """
    def test_executable(self):
        sim = msprime.TreeSimulator(10)
        line = sim.get_ms_command_line()
        self.assertEqual(line[0], "ms")
        line = sim.get_ms_command_line("otherms")
        self.assertEqual(line[0], "otherms")

    def test_sample_size(self):
        for n in [2, 10, 100]:
            sim = msprime.TreeSimulator(n)
            self.assertEqual(sim.get_ms_command_line()[1], str(n))

    def test_recombination(self):
        for m in [10, 100, 1000]:
            for r in [0.125, 1.0, 10]:
                rho = r * (m - 1)
                sim = msprime.TreeSimulator(10)
                sim.set_num_loci(m)
                sim.set_scaled_recombination_rate(r)
                args = sim.get_ms_command_line()
                self.assertEqual(args[-3], "-r")
                self.assertEqual(args[-2], str(rho))
                self.assertEqual(args[-1], str(m))

    def test_mutation(self):
        for m in [1, 10, 100, 1000]:
            for u in [0.125, 1.0, 10]:
                mu = u * m
                sim = msprime.TreeSimulator(10)
                sim.set_scaled_recombination_rate(1.0)
                sim.set_num_loci(m)
                args = sim.get_ms_command_line(scaled_mutation_rate=u)
                self.assertEqual(args[-2], "-t")
                self.assertEqual(args[-1], str(mu))

    def test_trees(self):
        sim = msprime.TreeSimulator(10)
        self.assertIn("-T", sim.get_ms_command_line())
        self.assertIn("-T", sim.get_ms_command_line(output_trees=True))
        self.assertNotIn("-T", sim.get_ms_command_line(output_trees=False))
        self.assertIn("-T", sim.get_ms_command_line(scaled_mutation_rate=1.0))
        self.assertIn("-T", sim.get_ms_command_line(
            scaled_mutation_rate=1.0, output_trees=True))
        self.assertNotIn("-T", sim.get_ms_command_line(
            scaled_mutation_rate=1.0, output_trees=False))

    def test_num_replicates(self):
        for j in [1, 100, 1000]:
            sim = msprime.TreeSimulator(10)
            args = sim.get_ms_command_line(num_replicates=j)
            self.assertEqual(str(j), args[2])

    # TODO Test population models.


class HighLevelTestCase(tests.MsprimeTestCase):
    """
    Superclass of tests on the high level interface.
    """

    def verify_sparse_tree_mrcas(self, st):
        # Check the mrcas
        oriented_forest = [st.get_parent(j) for j in range(st.get_root() + 1)]
        mrca_calc = tests.MRCACalculator(oriented_forest)
        # We've done exhaustive tests elsewhere, no need to go
        # through the combinations.
        for j in range(1, st.get_root() + 1):
            mrca = st.get_mrca(1, j)
            self.assertEqual(mrca, mrca_calc.get_mrca(1, j))
            self.assertEqual(st.get_time(mrca), st.get_tmrca(1, j))

    def verify_sparse_tree_branch_lengths(self, st):
        for j in range(1, st.get_sample_size() + 1):
            u = j
            while st.get_parent(u) != 0:
                l = st.get_time(st.get_parent(u)) - st.get_time(u)
                self.assertGreater(l, 0.0)
                self.assertEqual(st.get_branch_length(u), l)
                u = st.get_parent(u)

    def verify_sparse_tree_structure(self, st):
        used_nodes = set()
        for j in range(1, st.get_sample_size() + 1):
            self.assertEqual(st.get_time(j), 0)
            # verify the path to root
            u = j
            times = []
            while st.get_parent(u) != 0:
                used_nodes.add(u)
                v = st.get_parent(u)
                times.append(st.get_time(v))
                self.assertGreaterEqual(st.get_time(v), 0.0)
                self.assertIn(u, st.get_children(v))
                u = v
            self.assertEqual(u, st.get_root())
            self.assertEqual(times, sorted(times))
        used_nodes.add(st.get_root())
        self.assertEqual(len(used_nodes), 2 * st.get_sample_size() - 1)
        # for every entry other than used_nodes we should have an empty row
        for j in range(st.get_root()):
            if j not in used_nodes:
                self.assertEqual(st.get_parent(j), 0)
                self.assertEqual(st.get_time(j), 0)
                for c in st.get_children(j):
                    self.assertEqual(c, 0)
        # To a top-down traversal, and make sure we meet all the leaves.
        stack = [st.get_root()]
        leaves = []
        while len(stack) > 0:
            u = stack.pop()
            self.assertNotEqual(u, 0)
            if st.is_leaf(u):
                leaves.append(u)
                self.assertEqual(st.get_children(u)[0], 0)
                self.assertEqual(st.get_children(u)[1], 0)
            else:
                for c in reversed(st.get_children(u)):
                    stack.append(c)
            # Check that we get the correct number of leaves at each
            # node.
            self.assertEqual(st.get_num_leaves(u), len(list(st.leaves(u))))
            self.assertEqual(st.get_num_tracked_leaves(u), 0)
        self.assertEqual(
            sorted(leaves), list(range(1, st.get_sample_size() + 1)))
        # Check the parent dict
        pi = st.get_parent_dict()
        self.assertEqual(len(pi), 2 * st.get_sample_size() - 1)
        self.assertEqual(pi[st.get_root()], 0)
        for k, v in pi.items():
            self.assertEqual(st.get_parent(k), v)

    def verify_sparse_tree(self, st):
        self.verify_sparse_tree_mrcas(st)
        self.verify_sparse_tree_branch_lengths(st)
        self.verify_sparse_tree_structure(st)

    def verify_sparse_trees(self, ts):
        pts = tests.PythonTreeSequence(ts.get_ll_tree_sequence())
        iter1 = ts.trees()
        iter2 = pts.trees()
        length = 0
        for st1, st2 in zip(iter1, iter2):
            self.assertEqual(st1.get_sample_size(), ts.get_sample_size())
            root = 1
            while st1.get_parent(root) != 0:
                root = st1.get_parent(root)
            self.assertEqual(root, st1.get_root())
            self.assertEqual(st1, st2)
            self.assertFalse(st1 != st2)
            l, r = st1.get_interval()
            self.assertEqual(l, length)
            self.assertGreaterEqual(l, 0)
            self.assertGreater(r, l)
            self.assertLessEqual(r, ts.get_num_loci())
            length += r - l
            self.verify_sparse_tree(st1)
        self.assertEqual(length, ts.get_num_loci())
        self.assertRaises(StopIteration, next, iter1)
        self.assertRaises(StopIteration, next, iter2)

    def verify_haplotype_statistics(self, ts):
        """
        Verifies the statistics calculated for the haplotypes
        in the specified tree sequence.
        """
        pi1 = ts.get_pairwise_diversity()
        pi2 = get_pairwise_diversity(list(ts.haplotypes()))
        self.assertAlmostEqual(pi1, pi2)
        self.assertGreaterEqual(pi1, 0.0)
        self.assertFalse(math.isnan(pi1))

    def verify_mutations(self, ts):
        """
        Verify the mutations on this tree sequence make sense.
        """
        self.verify_haplotype_statistics(ts)
        all_mutations = list(ts.mutations())
        # Mutations must be sorted by position
        self.assertEqual(
            all_mutations, sorted(all_mutations, key=lambda x: x[0]))
        self.assertEqual(len(all_mutations), ts.get_num_mutations())
        all_tree_mutations = []
        for st in ts.trees():
            tree_mutations = list(st.mutations())
            self.assertEqual(st.get_num_mutations(), len(tree_mutations))
            all_tree_mutations.extend(tree_mutations)
            for position, node in tree_mutations:
                left, right = st.get_interval()
                self.assertTrue(left <= position < right)
                self.assertNotEqual(st.get_parent(node), 0)
        self.assertEqual(all_tree_mutations, all_mutations)
        pts = tests.PythonTreeSequence(ts.get_ll_tree_sequence())
        iter1 = ts.trees()
        iter2 = pts.trees()
        for st1, st2 in zip(iter1, iter2):
            self.assertEqual(st1, st2)


class TestSingleLocusSimulation(HighLevelTestCase):
    """
    Tests on the single locus simulations.
    """
    def test_simple_cases(self):
        for n in range(2, 10):
            st = msprime.simulate_tree(n)
            self.verify_sparse_tree(st)
        for n in [11, 13, 19, 101]:
            st = msprime.simulate_tree(n)
            self.verify_sparse_tree(st)

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
            st1 = msprime.simulate_tree(
                n, random_seed=1, population_models=[m1])
            st2 = msprime.simulate_tree(
                n, random_seed=1, population_models=[m2])
            self.assertEqual(st1, st2)
        # TODO add more tests!

    def test_initial_growth_rate(self):
        alpha = 5.0
        random_seed = 1234
        n = 10
        ts1 = msprime.simulate(
            n, random_seed=random_seed, population_models=[
                msprime.ExponentialPopulationModel(0.0, alpha)])
        sim = msprime.TreeSimulator(n)
        sim.set_random_seed(random_seed)
        sim.set_population_configurations([
            msprime.PopulationConfiguration(n, growth_rate=alpha)])
        sim.run()
        ts2 = sim.get_tree_sequence()
        self.assertEqual(list(ts1.records()), list(ts2.records()))

    def test_initial_pop_size(self):
        size = 5.0
        random_seed = 1234
        n = 10
        ts1 = msprime.simulate(
            n, random_seed=random_seed, population_models=[
                msprime.ConstantPopulationModel(0.0, size)])
        sim = msprime.TreeSimulator(n)
        sim.set_random_seed(random_seed)
        sim.set_population_configurations([
            msprime.PopulationConfiguration(n, initial_size=size)])
        sim.run()
        ts2 = sim.get_tree_sequence()
        self.assertEqual(list(ts1.records()), list(ts2.records()))

    def test_mixed_events(self):
        random_seed = 1234
        n = 10
        models = [
            msprime.ConstantPopulationModel(0.0, 2.0),
            msprime.ConstantPopulationModel(0.01, 5.0),
            msprime.ExponentialPopulationModel(0.5, 3.0)]
        ts1 = msprime.simulate(
            n, random_seed=random_seed, population_models=models)
        sim = msprime.TreeSimulator(n)
        sim.set_random_seed(random_seed)
        sim.set_population_configurations([
            msprime.PopulationConfiguration(n, initial_size=2)])
        sim.set_demographic_events([
            msprime.SizeChangeEvent(0.01, 5.0),
            msprime.GrowthRateChangeEvent(0.5, 3.0)])
        sim.run()
        ts2 = sim.get_tree_sequence()
        self.assertEqual(list(ts1.records()), list(ts2.records()))


class TestMultiLocusSimulation(HighLevelTestCase):
    """
    Tests on the single locus simulations.
    """

    def test_simple_cases(self):
        m = 1
        r = 0.1
        for n in range(2, 10):
            self.verify_sparse_trees(msprime.simulate(n, m, r))
        n = 4
        for m in range(1, 10):
            self.verify_sparse_trees(msprime.simulate(n, m, r))
        m = 100
        for r in [0.001, 0.01]:
            self.verify_sparse_trees(msprime.simulate(n, m, r))

    def test_error_cases(self):
        def f(n, m, r):
            return msprime.simulate(n, m, r)
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
        config = sim.get_configuration()
        self.assertEqual(config, parameters)

    def verify_dump_load(self, tree_sequence):
        """
        Dump the tree sequence and verify we can load again from the same
        file.
        """
        with tempfile.NamedTemporaryFile() as f:
            tree_sequence.dump(f.name)
            other = msprime.load(f.name)
        records = list(tree_sequence.records())
        other_records = list(other.records())
        self.assertEqual(records, other_records)

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
        seed = 1
        sim.set_random_seed(seed)
        self.assertEqual(sim.get_random_seed(), seed)
        sim.run()
        self.assertEqual(sim.get_num_breakpoints(), len(sim.get_breakpoints()))
        self.assertGreater(sim.get_used_memory(), 0)
        self.assertGreater(sim.get_time(), 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
        self.assertGreater(sim.get_max_memory(), 0)
        tree_sequence = sim.get_tree_sequence()
        t = 0.0
        for _, _, _, _, time in tree_sequence.records():
            if time > t:
                t = time
        self.assertEqual(sim.get_time(), t)
        self.assertGreater(sim.get_num_common_ancestor_events(), 0)
        self.assertGreaterEqual(sim.get_num_recombination_events(), 0)
        self.assertGreaterEqual(sim.get_total_num_migration_events(), 0)
        self.assertGreaterEqual(sim.get_num_multiple_recombination_events(), 0)
        self.verify_sparse_trees(tree_sequence)
        self.verify_parameters(sim, tree_sequence)
        self.verify_dump_load(tree_sequence)

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 100)
            r = random.random()
            self.verify_simulation(n, m, r)

    def test_perf_parameters(self):
        sim = msprime.TreeSimulator(10)
        # Before we call run, all perf parameters should be None
        self.assertIsNone(sim.get_avl_node_block_size())
        self.assertIsNone(sim.get_segment_block_size())
        self.assertIsNone(sim.get_node_mapping_block_size())
        self.assertIsNone(sim.get_coalescence_record_block_size())
        sim.run()
        self.assertGreater(sim.get_avl_node_block_size(), 0)
        self.assertGreater(sim.get_segment_block_size(), 0)
        self.assertGreater(sim.get_node_mapping_block_size(), 0)
        self.assertGreater(sim.get_coalescence_record_block_size(), 0)
        sim.reset()
        sim.set_avl_node_block_size(1)
        sim.set_segment_block_size(1)
        sim.set_node_mapping_block_size(1)
        sim.set_coalescence_record_block_size(1)
        self.assertEqual(sim.get_avl_node_block_size(), 1)
        self.assertEqual(sim.get_segment_block_size(), 1)
        self.assertEqual(sim.get_node_mapping_block_size(), 1)
        self.assertEqual(sim.get_coalescence_record_block_size(), 1)

    def test_bad_inputs(self):
        sim = msprime.TreeSimulator(10)
        for bad_type in ["xd", None, [], 4.4]:
            self.assertRaises(TypeError, msprime.TreeSimulator, bad_type)
            self.assertRaises(TypeError, sim.set_num_loci, bad_type)
        for bad_value in [-1, 0, 2**32]:
            self.assertRaises(ValueError, msprime.TreeSimulator, bad_value)
            self.assertRaises(ValueError, sim.set_num_loci, bad_value)
        self.assertRaises(ValueError, msprime.TreeSimulator, 1)
        self.assertRaises(ValueError, sim.set_scaled_recombination_rate, -1)


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
            self.assertGreater(ones, 0)
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
        haplotypes = list(tree_sequence.haplotypes())
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
            self.verify_simulation(n, m, r, theta)


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
            (st.get_length(), sparse_tree_to_newick(st, precision))
            for st in tree_sequence.trees()]
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
        bp = [0] + breakpoints + [tree_sequence.get_num_loci()]
        self.assertEqual(len(trees), len(bp) - 1)
        j = 0
        s = 0
        for length, _ in trees:
            self.assertGreater(length, 0)
            self.assertEqual(s, bp[j])
            s += length
            j += 1
        self.assertEqual(s, tree_sequence.get_num_loci())
        pts = tests.PythonTreeSequence(
            tree_sequence.get_ll_tree_sequence(), bp)
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
                    yield msprime.simulate(n, m, rho)

    def test_sparse_trees(self):
        for ts in self.get_example_tree_sequences():
            self.verify_sparse_trees(ts)

    def test_mutations(self):
        all_zero = True
        for ts in self.get_example_tree_sequences():
            ts.generate_mutations(0)
            self.assertEqual(ts.get_num_mutations(), 0)
            for st in ts.trees():
                self.assertEqual(st.get_num_mutations(), 0)
            # choose a mutation rate that hopefully guarantees mutations,
            # but not too many.
            mu = 100 / ts.get_num_loci()
            ts.generate_mutations(mu)
            if ts.get_num_mutations() > 0:
                all_zero = False
                self.verify_mutations(ts)
            muts = [[], [(0, 1)], [(0, 1), (0, 2)]]
            for mutations in muts:
                ts.set_mutations(mutations)
                self.assertEqual(ts.get_num_mutations(), len(mutations))
                self.assertEqual(list(ts.mutations()), mutations)
                self.verify_mutations(ts)
        self.assertFalse(all_zero)

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

    def verify_tracked_leaves(self, ts):
        # Should be empty list by default.
        for tree in ts.trees():
            for u in tree.nodes():
                self.assertEqual(tree.get_num_tracked_leaves(u), 0)
        tracked_leaves = [1, 2]
        for tree in ts.trees(tracked_leaves):
            nu = [0 for j in range(ts.get_num_nodes() + 1)]
            for j in tracked_leaves:
                u = j
                while u != 0:
                    nu[u] += 1
                    u = tree.get_parent(u)
            for u, count in enumerate(nu):
                self.assertEqual(tree.get_num_tracked_leaves(u), count)

    def test_tracked_leaves(self):
        for ts in self.get_example_tree_sequences():
            self.verify_tracked_leaves(ts)


class TestSparseTree(HighLevelTestCase):
    """
    Some simple tests on the API for the sparse tree.
    """
    def get_tree(self):
        return msprime.simulate_tree(
            10, random_seed=1, scaled_mutation_rate=1)

    def test_str(self):
        t = self.get_tree()
        self.assertIsInstance(str(t), str)
        self.assertEqual(str(t), str(t.get_parent_dict()))

    def test_leaves(self):
        t = self.get_tree()
        n = t.get_sample_size()
        all_leaves = list(t.leaves(t.get_root()))
        self.assertEqual(sorted(all_leaves), list(range(1, n + 1)))
        for j in range(1, n + 1):
            self.assertEqual(list(t.leaves(j)), [j])

        def test_func(t, u):
            """
            Simple test definition of the traversal.
            """
            stack = [u]
            while len(stack) > 0:
                v = stack.pop()
                if t.is_internal(v):
                    for c in reversed(t.get_children(v)):
                        stack.append(c)
                else:
                    yield v
        for u in t.nodes():
            l1 = list(t.leaves(u))
            l2 = list(test_func(t, u))
            self.assertEqual(l1, l2)
            self.assertEqual(t.get_num_leaves(u), len(l1))

    def test_draw(self):
        t = self.get_tree()
        with tempfile.NamedTemporaryFile() as f:
            w = 123
            h = 456
            t.draw(f.name, w, h, show_times=True)
            self.assertGreater(os.path.getsize(f.name), 0)
            # Check some basic stuff about the SVG output.
            f.seek(0)
            root = xml.etree.ElementTree.fromstring(f.read())
            self.assertEqual(root.tag, "{http://www.w3.org/2000/svg}svg")
            width = int(root.attrib["width"])
            self.assertEqual(w, width)
            height = int(root.attrib["height"])
            self.assertEqual(h, height)

    def test_traversals(self):
        t1 = self.get_tree()
        t2 = tests.PythonSparseTree.from_sparse_tree(t1)
        self.assertEqual(list(t1.nodes()), list(t2.nodes()))
        self.assertEqual(list(t1.nodes()), list(t1.nodes(t1.get_root())))
        self.assertEqual(
            list(t1.nodes()),
            list(t1.nodes(t1.get_root(), "preorder")))
        for u in t1.nodes():
            self.assertEqual(list(t1.nodes(u)), list(t2.nodes(u)))
        self.assertRaises(ValueError, t1.nodes, None, "bad order")
