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
import json
import os
import random
import tempfile
import unittest
import xml.etree

import msprime
import tests


def simple_get_pairwise_diversity(haplotypes):
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


def get_pairwise_diversity(tree_sequence, samples=None):
    """
    This is the exact algorithm used by the low-level C code
    and should return identical results.
    """
    if samples is None:
        tracked_leaves = list(range(tree_sequence.get_sample_size()))
    else:
        tracked_leaves = list(samples)
    if len(tracked_leaves) < 2:
        raise ValueError("len(samples) must be >= 2")
    pi = 0
    k = len(tracked_leaves)
    denom = k * (k - 1) / 2
    for t in tree_sequence.trees(tracked_leaves=tracked_leaves):
        for _, node in t.mutations():
            j = t.get_num_tracked_leaves(node)
            pi += j * (k - j) / denom
    return pi


def sparse_tree_to_newick(st, precision, Ne):
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
                length = (st.get_time(node) - st.get_time(child)) / (4 * Ne)
                s = "{0:.{1}f}".format(length, precision)
                branch_lengths[child] = s
    return _build_newick(root, root, st, branch_lengths)


def _build_newick(node, root, tree, branch_lengths):
    c1, c2 = tree.get_children(node)
    if c1 != msprime.NULL_NODE:
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
        s = "{0}:{1}".format(node + 1, branch_lengths[node])
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
        L = 1
        recomb_map = msprime.RecombinationMap.uniform_map(
            length=L, rate=0, num_loci=L)
        sim = msprime.simulator_factory(10, recombination_map=recomb_map)
        line = sim.get_ms_command_line()
        self.assertEqual(line[0], "ms")
        line = sim.get_ms_command_line("otherms")
        self.assertEqual(line[0], "otherms")

    def test_sample_size(self):
        L = 1
        recomb_map = msprime.RecombinationMap.uniform_map(
            length=L, rate=0, num_loci=L)
        for n in [2, 10, 100]:
            sim = msprime.simulator_factory(n, recombination_map=recomb_map)
            self.assertEqual(sim.get_ms_command_line()[1], str(n))

    def test_recombination(self):
        for L in [10, 100, 1000]:
            for r in [0.125, 1.0, 10]:
                rho = r * L
                recomb_map = msprime.RecombinationMap.uniform_map(
                    length=L, rate=r, num_loci=L)
                sim = msprime.simulator_factory(
                    10, recombination_map=recomb_map)
                args = sim.get_ms_command_line()
                self.assertEqual(args[-3], "-r")
                self.assertEqual(float(args[-2]), float(str(rho)))
                self.assertEqual(float(args[-1]), L)

    def test_mutation(self):
        for L in [1, 10, 100, 1000]:
            for u in [0.125, 1.0, 10]:
                mu = u * L * 4
                recomb_map = msprime.RecombinationMap.uniform_map(
                    length=L, rate=0, num_loci=L)
                sim = msprime.simulator_factory(
                    10, recombination_map=recomb_map)
                args = sim.get_ms_command_line(mutation_rate=u)
                self.assertEqual(args[-2], "-t")
                self.assertEqual(float(args[-1]), mu)

    def test_trees(self):
        r = 0
        L = 1
        recomb_map = msprime.RecombinationMap.uniform_map(
            length=L, rate=r, num_loci=L)
        sim = msprime.simulator_factory(
            10, recombination_map=recomb_map)
        self.assertIn("-T", sim.get_ms_command_line())
        self.assertIn("-T", sim.get_ms_command_line(output_trees=True))
        self.assertNotIn("-T", sim.get_ms_command_line(output_trees=False))
        self.assertIn("-T", sim.get_ms_command_line(mutation_rate=1.0))
        self.assertIn("-T", sim.get_ms_command_line(
            mutation_rate=1.0, output_trees=True))
        self.assertNotIn("-T", sim.get_ms_command_line(
            mutation_rate=1.0, output_trees=False))

    def test_num_replicates(self):
        L = 1
        for j in [1, 100, 1000]:
            recomb_map = msprime.RecombinationMap.uniform_map(
                length=L, rate=0, num_loci=L)
            sim = msprime.simulator_factory(
                10, recombination_map=recomb_map)
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
        for j in range(st.get_root() + 1):
            mrca = st.get_mrca(0, j)
            self.assertEqual(mrca, mrca_calc.get_mrca(0, j))
            if mrca != msprime.NULL_NODE:
                self.assertEqual(st.get_time(mrca), st.get_tmrca(0, j))

    def verify_sparse_tree_branch_lengths(self, st):
        for j in range(st.get_sample_size()):
            u = j
            while st.get_parent(u) != msprime.NULL_NODE:
                l = st.get_time(st.get_parent(u)) - st.get_time(u)
                self.assertGreater(l, 0.0)
                self.assertEqual(st.get_branch_length(u), l)
                u = st.get_parent(u)

    def verify_sparse_tree_structure(self, st):
        used_nodes = set()
        for j in range(st.get_sample_size()):
            self.assertEqual(st.get_time(j), 0)
            # verify the path to root
            u = j
            times = []
            while st.get_parent(u) != msprime.NULL_NODE:
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
                self.assertEqual(st.get_parent(j), msprime.NULL_NODE)
                self.assertEqual(st.get_time(j), 0)
                for c in st.get_children(j):
                    self.assertEqual(c, msprime.NULL_NODE)
        # To a top-down traversal, and make sure we meet all the leaves.
        stack = [st.get_root()]
        leaves = []
        while len(stack) > 0:
            u = stack.pop()
            self.assertNotEqual(u, msprime.NULL_NODE)
            if st.is_leaf(u):
                leaves.append(u)
                self.assertEqual(st.get_children(u)[0], msprime.NULL_NODE)
                self.assertEqual(st.get_children(u)[1], msprime.NULL_NODE)
            else:
                for c in reversed(st.get_children(u)):
                    stack.append(c)
            # Check that we get the correct number of leaves at each
            # node.
            self.assertEqual(st.get_num_leaves(u), len(list(st.leaves(u))))
            self.assertEqual(st.get_num_tracked_leaves(u), 0)
        self.assertEqual(sorted(leaves), list(range(st.get_sample_size())))
        # Check the parent dict
        pi = st.get_parent_dict()
        self.assertEqual(len(pi), 2 * st.get_sample_size() - 1)
        self.assertEqual(pi[st.get_root()], msprime.NULL_NODE)
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
        num_trees = 0
        breakpoints = [0]
        for st1, st2 in zip(iter1, iter2):
            self.assertEqual(st1.get_sample_size(), ts.get_sample_size())
            root = 0
            while st1.get_parent(root) != msprime.NULL_NODE:
                root = st1.get_parent(root)
            self.assertEqual(root, st1.get_root())
            self.assertEqual(st1, st2)
            self.assertFalse(st1 != st2)
            l, r = st1.get_interval()
            breakpoints.append(r)
            self.assertEqual(l, length)
            self.assertGreaterEqual(l, 0)
            self.assertGreater(r, l)
            self.assertLessEqual(r, ts.get_sequence_length())
            length += r - l
            self.verify_sparse_tree(st1)
            num_trees += 1
        self.assertEqual(breakpoints, list(ts.breakpoints()))
        self.assertEqual(length, ts.get_sequence_length())
        self.assertEqual(ts.get_num_trees(), num_trees)
        self.assertRaises(StopIteration, next, iter1)
        self.assertRaises(StopIteration, next, iter2)

    def verify_haplotype_statistics(self, ts):
        """
        Verifies the statistics calculated for the haplotypes
        in the specified tree sequence.
        """
        haplotypes = list(ts.haplotypes())
        pi1 = ts.get_pairwise_diversity()
        pi2 = simple_get_pairwise_diversity(haplotypes)
        pi3 = get_pairwise_diversity(ts)
        self.assertAlmostEqual(pi1, pi2)
        self.assertEqual(pi1, pi3)
        self.assertGreaterEqual(pi1, 0.0)
        self.assertFalse(math.isnan(pi1))
        # Check for a subsample.
        samples = range(ts.get_sample_size() // 2 + 1)
        pi1 = ts.get_pairwise_diversity(samples)
        pi2 = simple_get_pairwise_diversity([haplotypes[j] for j in samples])
        pi3 = get_pairwise_diversity(ts, samples)
        self.assertAlmostEqual(pi1, pi2)
        self.assertEqual(pi1, pi3)
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
                self.assertNotEqual(st.get_parent(node), msprime.NULL_NODE)
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
            st = next(msprime.simulate(n).trees())
            self.verify_sparse_tree(st)
        for n in [11, 13, 19, 101]:
            st = next(msprime.simulate(n).trees())
            self.verify_sparse_tree(st)

    def test_models(self):
        # Exponential growth of 0 and constant model should be identical.
        for n in [2, 10, 100]:
            m1 = msprime.PopulationConfiguration(n, growth_rate=0)
            m2 = msprime.PopulationConfiguration(n, initial_size=1.0)
            st1 = next(msprime.simulate(
                n, random_seed=1, population_configurations=[m1]).trees())
            st2 = next(msprime.simulate(
                n, random_seed=1, population_configurations=[m2]).trees())
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
            self.verify_sparse_trees(msprime.simulate(n, m, r))
        n = 4
        for m in range(1, 10):
            self.verify_sparse_trees(msprime.simulate(n, m, r))
        m = 100
        for r in [0.001, 0.01]:
            self.verify_sparse_trees(msprime.simulate(n, m, r))

    def test_error_cases(self):
        def f(n, m, r):
            return msprime.simulate(
                sample_size=n, length=m, recombination_rate=r)
        for n in [-100, -1, 0, 1, None]:
            self.assertRaises(ValueError, f, n, 1, 1.0)
        for n in ["", "2", 2.2, 1e5]:
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
            sim.get_per_locus_scaled_recombination_rate())
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
        recomb_map = msprime.RecombinationMap.uniform_map(m, r, num_loci=m)
        rng = msprime.RandomGenerator(1)
        sim = msprime.simulator_factory(
            n, recombination_map=recomb_map, random_generator=rng)
        self.assertEqual(sim.get_random_generator(), rng)
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
        for record in tree_sequence.records():
            if record.time > t:
                t = record.time
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
        sim = msprime.simulator_factory(10)
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
        recomb_map = msprime.RecombinationMap.uniform_map(1, 0)
        for bad_type in ["xd", None, [], 4.4]:
            self.assertRaises(
                TypeError, msprime.TreeSimulator, bad_type, recomb_map)
            self.assertRaises(
                TypeError, msprime.TreeSimulator, 2, bad_type)
        for bad_value in [-1, 0, 2**32]:
            self.assertRaises(
                ValueError, msprime.TreeSimulator, bad_value, recomb_map)
        self.assertRaises(ValueError, msprime.TreeSimulator, 1, recomb_map)


class TestHaplotypeGenerator(HighLevelTestCase):
    """
    Tests the haplotype generation code.
    """

    def verify_haplotypes_variants(self, n, haplotypes, variants):
        """
        Verify that the specified set of haplotypes and variants are
        consistent.
        """
        self.assertEqual(len(haplotypes), n)
        m = len(haplotypes[0])
        for h in haplotypes:
            self.assertEqual(len(h), m)
        self.assertEqual(len(variants), m)
        # Examine each column in H; we must have a mixture of 0s and 1s
        for k in range(m):
            zeros = 0
            ones = 0
            col = ""
            for j in range(n):
                b = haplotypes[j][k]
                zeros += b == '0'
                ones += b == '1'
                col += b
            self.assertGreater(zeros, 0)
            self.assertGreater(ones, 0)
            self.assertEqual(zeros + ones, n)
            self.assertEqual(col, variants[k])

    def verify_simulation(self, n, m, r, theta):
        """
        Verifies a simulation for the specified parameters.
        """
        recomb_map = msprime.RecombinationMap.uniform_map(m, r, m)
        tree_sequence = msprime.simulate(
            n, recombination_map=recomb_map, mutation_rate=theta)
        haplotypes = list(tree_sequence.haplotypes())
        for h in haplotypes:
            self.assertEqual(len(h), tree_sequence.get_num_mutations())
        variants = [variant for _, variant in tree_sequence.variants()]
        for v in variants:
            self.assertEqual(len(v), tree_sequence.get_sample_size())
        self.verify_haplotypes_variants(n, haplotypes, variants)
        self.assertEqual(
            [pos for pos, _ in tree_sequence.variants()],
            [pos for pos, _ in tree_sequence.mutations()])

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
    def verify_trees(self, tree_sequence, breakpoints, Ne):
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
            (st.get_length(), sparse_tree_to_newick(st, precision, Ne))
            for st in tree_sequence.trees()]
        new_trees = list(tree_sequence.newick_trees(precision, Ne=Ne))
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
        bp = [0] + breakpoints + [tree_sequence.get_sequence_length()]
        self.assertEqual(len(trees), len(bp) - 1)
        j = 0
        s = 0
        for length, _ in trees:
            self.assertGreater(length, 0)
            self.assertEqual(s, bp[j])
            s += length
            j += 1
        self.assertEqual(s, tree_sequence.get_sequence_length())
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
            (2, 1, 0, 0.25),
            (2, 10, 0.1, 1),
            (4, 10, 0.1, 100),
            (10, 10, 0.1,  10),
            (20, 1, 0, 1025),
            (20, 10, 0.1, 100),
            (10, 50, 1.0, 1e6),
        ]
        for n, m, r, Ne in cases:
            recomb_map = msprime.RecombinationMap.uniform_map(m, r, m)
            ts = msprime.simulator_factory(
                n, Ne=Ne, recombination_map=recomb_map)
            ts.run()
            tree_sequence = ts.get_tree_sequence()
            breakpoints = ts.get_breakpoints()
            self.verify_trees(tree_sequence, breakpoints, Ne)
            self.verify_all_breakpoints(tree_sequence, breakpoints)

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 100)
            r = random.random()
            Ne = random.uniform(1, 20)
            recomb_map = msprime.RecombinationMap.uniform_map(m, r, m)
            ts = msprime.simulator_factory(
                    n, Ne=Ne, recombination_map=recomb_map)
            ts.run()
            tree_sequence = ts.get_tree_sequence()
            breakpoints = ts.get_breakpoints()
            self.verify_trees(tree_sequence, breakpoints, Ne)
            self.verify_all_breakpoints(tree_sequence, breakpoints)


class TestTreeSequence(HighLevelTestCase):
    """
    Tests for the tree sequence object.
    """
    def get_example_tree_sequences(self):
        for n in [2, 3, 10, 100]:
            for m in [1, 2, 32]:
                for rho in [0, 0.1, 0.5]:
                    recomb_map = msprime.RecombinationMap.uniform_map(
                        m, rho, num_loci=m)
                    ts = msprime.simulate(
                        n, recombination_map=recomb_map)
                    yield ts

    def test_sparse_trees(self):
        for ts in self.get_example_tree_sequences():
            self.verify_sparse_trees(ts)

    def test_mutations(self):
        rng = msprime.RandomGenerator(3)
        all_zero = True
        for ts in self.get_example_tree_sequences():
            ts.generate_mutations(0, rng)
            self.assertEqual(ts.get_num_mutations(), 0)
            for st in ts.trees():
                self.assertEqual(st.get_num_mutations(), 0)
            # choose a mutation rate that hopefully guarantees mutations,
            # but not too many.
            mu = 10 / ts.get_sequence_length()
            ts.generate_mutations(mu, rng)
            if ts.get_num_mutations() > 0:
                all_zero = False
                self.verify_mutations(ts)
            muts = [[], [(0, 0)], [(0, 0), (0, 1)]]
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
        tracked_leaves = [0, 1]
        for tree in ts.trees(tracked_leaves):
            nu = [0 for j in range(ts.get_num_nodes())]
            for j in tracked_leaves:
                u = j
                while u != msprime.NULL_NODE:
                    nu[u] += 1
                    u = tree.get_parent(u)
            for u, count in enumerate(nu):
                self.assertEqual(tree.get_num_tracked_leaves(u), count)

    def test_tracked_leaves(self):
        for ts in self.get_example_tree_sequences():
            self.verify_tracked_leaves(ts)

    def test_get_pairwise_diversity(self):
        for ts in self.get_example_tree_sequences():
            n = ts.get_sample_size()
            self.assertRaises(ValueError, ts.get_pairwise_diversity, [])
            self.assertRaises(ValueError, ts.get_pairwise_diversity, [1])
            self.assertRaises(ValueError, ts.get_pairwise_diversity, [1, n])
            self.assertEqual(
                ts.get_pairwise_diversity(),
                ts.get_pairwise_diversity(range(n)))
            self.assertEqual(
                ts.get_pairwise_diversity([0, 1]),
                ts.get_pairwise_diversity([1, 0]))

    def test_get_population(self):
        for ts in self.get_example_tree_sequences():
            n = ts.get_sample_size()
            self.assertRaises(ValueError, ts.get_population, -1)
            self.assertRaises(ValueError, ts.get_population, n)
            self.assertRaises(ValueError, ts.get_population, n + 1)
            self.assertEqual(ts.get_population(0), 0)
            self.assertEqual(ts.get_population(n - 1), 0)

    def test_get_samples(self):
        for ts in self.get_example_tree_sequences():
            n = ts.get_sample_size()
            samples = list(range(n))
            self.assertEqual(ts.get_samples(), samples)
            self.assertEqual(ts.get_samples(0), samples)
            self.assertEqual(ts.get_samples(msprime.NULL_POPULATION), [])
            self.assertEqual(ts.get_samples(1), [])


class TestSparseTree(HighLevelTestCase):
    """
    Some simple tests on the API for the sparse tree.
    """
    def get_tree(self):
        ts = msprime.simulate(10, random_seed=1, mutation_rate=1)
        return next(ts.trees())

    def test_str(self):
        t = self.get_tree()
        self.assertIsInstance(str(t), str)
        self.assertEqual(str(t), str(t.get_parent_dict()))

    def test_leaves(self):
        t = self.get_tree()
        n = t.get_sample_size()
        all_leaves = list(t.leaves(t.get_root()))
        self.assertEqual(sorted(all_leaves), list(range(n)))
        for j in range(n):
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

    def test_total_branch_length(self):
        t1 = self.get_tree()
        bl = 0
        root = t1.get_root()
        for node in t1.nodes():
            if node != root:
                bl += t1.get_branch_length(node)
        self.assertGreater(bl, 0)
        self.assertEqual(t1.get_total_branch_length(), bl)


class TestRecombinationMap(unittest.TestCase):
    """
    Tests the code for recombination map.
    """

    def verify_coordinate_conversion(self, positions, rates):
        """
        Verifies coordinate conversions by the specified RecombinationMap
        instance.
        """
        num_loci = 10
        rm = msprime.RecombinationMap(positions, rates, num_loci)
        other_rm = tests.PythonRecombinationMap(positions, rates, num_loci)
        self.assertEqual(
            rm.get_total_recombination_rate(),
            other_rm.get_total_recombination_rate())
        num_random_trials = 10
        num_systematic_trials = 10
        values = [random.random() for j in range(num_random_trials)]
        for j in range(num_systematic_trials):
            values.append(j * 1 / num_systematic_trials)
        values += positions
        for x in values:
            # x is a physical coordinate
            y = rm.physical_to_genetic(x)
            self.assertEqual(y, other_rm.physical_to_genetic(x))
            self.assertTrue(0 <= y <= num_loci)
            z = rm.genetic_to_physical(y)
            self.assertAlmostEqual(x, z)

            # Now x is a genetic coordinate
            y = rm.genetic_to_physical(x)
            self.assertTrue(0 <= y <= 1)
            self.assertAlmostEqual(y, other_rm.genetic_to_physical(x))
            z = rm.physical_to_genetic(y)
            self.assertAlmostEqual(x, z)

    def test_zero_rate_values(self):
        # When we have a zero rate in some interval we no longer have a
        # bijective function, since all the physical coordinates in this
        # interval map to a single genetic coordinate.
        positions = [0, 0.25, 0.5, 0.75, 1]
        rates = [1, 0, 1, 0, 0]
        num_loci = 100
        rm = msprime.RecombinationMap(positions, rates, num_loci)
        other_rm = tests.PythonRecombinationMap(positions, rates, num_loci)
        self.assertEqual(0.5, rm.get_total_recombination_rate())
        self.assertEqual(0.5, other_rm.get_total_recombination_rate())
        # Between 0 and 0.25 and 0.5 and 0.75 we should be able to map 1-1
        # in physical coordinates.
        for x in [0, 0.125, 0.25, 0.50001, 0.66, 0.75]:
            y = rm.physical_to_genetic(x)
            self.assertEqual(y, other_rm.physical_to_genetic(x))
            self.assertTrue(0 <= y <= num_loci)
            z = rm.genetic_to_physical(y)
            self.assertAlmostEqual(x, z)
        # All physical coordinates within the 0 region should map down to
        # the first point.
        for start, end in [(0.25, 0.5), (0.75, 1)]:
            for x in [start + delta for delta in [0, 0.01, 0.1]] + [end]:
                y = rm.physical_to_genetic(x)
                self.assertEqual(y, other_rm.physical_to_genetic(x))
                self.assertTrue(0 <= y <= num_loci)
                z = rm.genetic_to_physical(y)
                self.assertEqual(z, start)

    def test_one_rate(self):
        num_loci = 1024
        for rate in [0.1, 1.0, 10]:
            positions = [0, 1]
            rates = [rate, 0]
            rm = msprime.RecombinationMap(positions, rates, num_loci)
            self.assertEqual(rate, rm.get_total_recombination_rate())
            self.verify_coordinate_conversion(positions, rates)

    def test_simple_map(self):
        positions = [0, 0.25, 0.5, 0.75, 1]
        rates = [0.125, 0.25, 0.5, 0.75, 0]
        self.verify_coordinate_conversion(positions, rates)

    def test_random_map(self):
        for size in [2, 3, 4, 100]:
            positions = [0] + sorted(
                random.random() for _ in range(size - 2)) + [1]
            rates = [random.random() for _ in range(size - 1)] + [0]
            self.verify_coordinate_conversion(positions, rates)

    def test_zero_rate(self):
        positions = [0, 1]
        rates = [0, 0]
        for m in [1, 10]:
            rm = msprime.RecombinationMap(positions, rates, m)
            other_rm = tests.PythonRecombinationMap(positions, rates, m)
            self.assertEqual(0.0, rm.get_total_recombination_rate())
            self.assertEqual(0.0, other_rm.get_total_recombination_rate())
            # All values should map directly to themselves.
            for x in [0, 0.24, 0.33, 0.99, 1]:
                self.assertEqual(rm.genetic_to_physical(m * x), x)
                self.assertEqual(other_rm.genetic_to_physical(m * x), x)
                self.assertEqual(other_rm.physical_to_genetic(x), x * m)
                self.assertEqual(rm.physical_to_genetic(x), x * m)

    def test_simple_examples(self):
        rm = msprime.RecombinationMap([0, 0.9, 1], [2, 1, 0], 10)
        self.assertAlmostEqual(rm.get_total_recombination_rate(), 1.9)
        rm = msprime.RecombinationMap([0, 0.5, 0.6, 1], [2, 1, 2, 0], 100)
        self.assertAlmostEqual(rm.get_total_recombination_rate(), 1.9)

    def test_read_hapmap(self):
        with tempfile.NamedTemporaryFile("w+") as f:
            print("HEADER", file=f)
            print("chr1 0 1", file=f)
            print("chr1 1 5 x", file=f)
            print("s    2 0 x x x", file=f)
            f.flush()
            rm = msprime.RecombinationMap.read_hapmap(f.name)
            self.assertEqual(rm.get_positions(), [0, 1, 2])
            self.assertEqual(rm.get_rates(), [1e-8, 5e-8, 0])


def get_ll_demographic_events(ll_sim):
    """
    Utility function to get the low-level demographic events from the
    specified low-level simulator.
    """
    d = json.loads(ll_sim.get_configuration_json())
    return d["demographic_events"]


class TestSimulatorFactory(unittest.TestCase):
    """
    Tests that the simulator factory high-level function correctly
    creates simulators with the required parameter values.
    """
    def test_default_random_seed(self):
        sim = msprime.simulator_factory(10)
        rng = sim.get_random_generator()
        self.assertIsInstance(rng, msprime.RandomGenerator)
        self.assertNotEqual(rng.get_seed(), 0)

    def test_random_seed(self):
        seed = 12345
        rng = msprime.RandomGenerator(seed)
        sim = msprime.simulator_factory(10, random_generator=rng)
        rng = sim.get_random_generator()
        self.assertEqual(rng, sim.get_random_generator())
        self.assertEqual(rng.get_seed(), seed)

    def test_sample_size(self):
        self.assertRaises(ValueError, msprime.simulator_factory)
        self.assertRaises(ValueError, msprime.simulator_factory, 1)
        self.assertRaises(
            ValueError, msprime.simulator_factory, sample_size=1)
        for n in [2, 100, 1000]:
            sim = msprime.simulator_factory(n)
            self.assertEqual(sim.get_sample_size(), n)
            ll_sim = sim.create_ll_instance()
            self.assertEqual(ll_sim.get_sample_size(), n)

    def test_effective_population_size(self):
        def f(Ne):
            return msprime.simulator_factory(10, Ne=Ne)
        for bad_value in [-1, -1e16, 0]:
            self.assertRaises(ValueError, f, bad_value)
        for Ne in [1, 10, 1e5]:
            sim = f(Ne)
            self.assertEqual(sim.get_effective_population_size(), Ne)
        # Test the default.
        sim = msprime.simulator_factory(10)

    def test_population_configurations(self):
        def f(configs):
            return msprime.simulator_factory(
                population_configurations=configs)
        for bad_type in [10, ["sdf"], "sdfsd"]:
            self.assertRaises(TypeError, f, bad_type)
        # Just test the basic equalities here. The actual
        # configuration options are tested elewhere.
        for N in range(1, 10):
            pop_configs = [
                msprime.PopulationConfiguration(5) for _ in range(N)]
            sample_size = 5 * N
            sim = msprime.simulator_factory(
                population_configurations=pop_configs)
            self.assertEqual(
                sim.get_population_configurations(), pop_configs)
            self.assertEqual(
                sim.get_sample_size(), sample_size)
            ll_sim = sim.create_ll_instance()
            self.assertEqual(len(ll_sim.get_population_configuration()), N)
        # The default is a single population
        sim = msprime.simulator_factory(10)
        ll_sim = sim.create_ll_instance()
        self.assertEqual(len(ll_sim.get_population_configuration()), 1)

    def test_sample_size_population_configuration(self):
        for d in range(1, 5):
            # Zero sample size is always an error
            configs = [msprime.PopulationConfiguration(0) for _ in range(d)]
            self.assertRaises(
                ValueError, msprime.simulator_factory,
                population_configurations=configs)
            # Disagreements between sample sizes raises and error
            configs = [msprime.PopulationConfiguration(1) for _ in range(d)]
            self.assertRaises(
                ValueError,  msprime.simulator_factory,
                sample_size=d - 1,
                population_configurations=configs)

    def test_migration_matrix(self):
        # Cannot specify a migration matrix without population
        # configurations
        self.assertRaises(
            ValueError, msprime.simulator_factory, 10,
            migration_matrix=[])
        for N in range(1, 10):
            pop_configs = [
                msprime.PopulationConfiguration(5) for _ in range(N)]
            sim = msprime.simulator_factory(
                N * 5, population_configurations=pop_configs)
            ll_sim = sim.create_ll_instance()
            # If we don't specify a matrix, it's 0 everywhere.
            matrix = [0 for j in range(N * N)]
            self.assertEqual(ll_sim.get_migration_matrix(), matrix)

            def f(hl_matrix):
                return msprime.simulator_factory(
                    N * 5, population_configurations=pop_configs,
                    migration_matrix=hl_matrix)
            hl_matrix = [
                [(j + k) * int(j != k) for j in range(N)] for k in range(N)]
            sim = f(hl_matrix)
            self.assertEqual(sim.get_migration_matrix(), hl_matrix)
            ll_sim = sim.create_ll_instance()
            Ne = sim.get_effective_population_size()
            ll_matrix = [4 * Ne * v for row in hl_matrix for v in row]
            self.assertEqual(ll_sim.get_migration_matrix(), ll_matrix)
            for bad_type in ["", {}, 234]:
                self.assertRaises(TypeError, f, bad_type)
            # Now check for the structure of the matrix.
            hl_matrix[0][0] = "bad value"
            sim = f(hl_matrix)
            self.assertRaises(TypeError, sim.create_ll_instance)
            hl_matrix[0] = None
            self.assertRaises(TypeError, f, hl_matrix)
            hl_matrix[0] = []
            self.assertRaises(ValueError, f, hl_matrix)

    def test_default_migration_matrix(self):
        sim = msprime.simulator_factory(10)
        ll_sim = sim.create_ll_instance()
        self.assertEqual(ll_sim.get_migration_matrix(), [0.0])

    def test_demographic_events(self):
        for bad_type in ["sdf", 234, [12], [None]]:
            self.assertRaises(
                TypeError, msprime.simulator_factory, 2,
                demographic_events=bad_type)
        # TODO test for bad values.

    def test_default_demographic_events(self):
        sim = msprime.simulator_factory(10)
        ll_sim = sim.create_ll_instance()
        self.assertEqual(get_ll_demographic_events(ll_sim), [])

    def test_recombination_rate(self):
        def f(recomb_rate):
            return msprime.simulator_factory(
                10, recombination_rate=recomb_rate)
        for bad_type in ["", {}, []]:
            self.assertRaises(TypeError, f, bad_type)
        for bad_value in [-1, -1e15]:
            self.assertRaises(ValueError, f, bad_value)
        for rate in [0, 1e-3, 10]:
            sim = f(rate)
            recomb_map = sim.get_recombinatation_map()
            self.assertEqual(recomb_map.get_positions(), [0, 1], [rate, 0])
            self.assertEqual(
                recomb_map.get_num_loci(),
                msprime.RecombinationMap.DEFAULT_NUM_LOCI)

    def test_scaled_recombination_rate(self):
        values = [
            (10, 0.1, 0.1),
            (0.1, 1, 10),
            (1e-8, 10**4, 10**8),
            (1e-8, 10**5, 10**9),
        ]
        for rate, Ne, length in values:
            sim = msprime.simulator_factory(
                10, Ne=Ne, recombination_rate=rate, length=length)
            num_loci = msprime.RecombinationMap.DEFAULT_NUM_LOCI
            total_rate = 4 * Ne * length * rate
            per_locus_rate = total_rate / (num_loci - 1)
            # We expect all these rates to be positive.
            self.assertGreater(per_locus_rate, 0)
            ll_sim = sim.create_ll_instance()
            self.assertAlmostEqual(
                per_locus_rate, ll_sim.get_scaled_recombination_rate())

    def test_recombination_map(self):
        def f(recomb_map):
            return msprime.simulator_factory(
                10, recombination_map=recomb_map)
        self.assertRaises(TypeError, f, "wrong type")
        for n in range(2, 10):
            positions = list(range(n))
            rates = [0.1 * j for j in range(n - 1)] + [0.0]
            recomb_map = msprime.RecombinationMap(positions, rates)
            sim = msprime.simulator_factory(10, recombination_map=recomb_map)
            self.assertEqual(sim.get_recombinatation_map(), recomb_map)
            self.assertEqual(recomb_map.get_positions(), positions)
            self.assertEqual(recomb_map.get_rates(), rates)
            ll_sim = sim.create_ll_instance()
            self.assertEqual(ll_sim.get_num_loci(), recomb_map.get_num_loci())

    def test_combining_recomb_map_and_rate_length(self):
        recomb_map = msprime.RecombinationMap([0, 1], [1, 0])
        self.assertRaises(
            ValueError, msprime.simulator_factory, 10,
            recombination_map=recomb_map, length=1)
        self.assertRaises(
            ValueError, msprime.simulator_factory, 10,
            recombination_map=recomb_map, recombination_rate=100)
        self.assertRaises(
            ValueError, msprime.simulator_factory, 10,
            recombination_map=recomb_map, length=1,
            recombination_rate=1)


class TestSimulateInterface(unittest.TestCase):
    """
    Some simple test cases for the simulate() interface.
    """
    def test_defaults(self):
        n = 10
        ts = msprime.simulate(n)
        self.assertIsInstance(ts, msprime.TreeSequence)
        self.assertEqual(ts.get_sample_size(), n)
        self.assertEqual(ts.get_num_trees(), 1)
        self.assertEqual(ts.get_num_mutations(), 0)

    def test_replicates(self):
        n = 20
        num_replicates = 10
        count = 0
        for ts in msprime.simulate(n, num_replicates=num_replicates):
            count += 1
            self.assertIsInstance(ts, msprime.TreeSequence)
            self.assertEqual(ts.get_sample_size(), n)
            self.assertEqual(ts.get_num_trees(), 1)
        self.assertEqual(num_replicates, count)

    def test_mutations(self):
        n = 10
        ts = msprime.simulate(n, mutation_rate=10)
        self.assertIsInstance(ts, msprime.TreeSequence)
        self.assertEqual(ts.get_sample_size(), n)
        self.assertEqual(ts.get_num_trees(), 1)
        self.assertGreater(ts.get_num_mutations(), 0)

    def test_recombination(self):
        n = 10
        ts = msprime.simulate(n, recombination_rate=10)
        self.assertIsInstance(ts, msprime.TreeSequence)
        self.assertEqual(ts.get_sample_size(), n)
        self.assertGreater(ts.get_num_trees(), 1)
        self.assertEqual(ts.get_num_mutations(), 0)
