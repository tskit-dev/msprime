#
# Copyright (C) 2015-2017 University of Oxford
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

import collections
import gzip
import itertools
import json
import math
import os
import random
import shutil
import sys
import six
import tempfile
import unittest
import xml.etree

import numpy as np

import msprime
import _msprime
import tests


def get_example_tree_sequences():
    for n in [2, 3, 10, 100]:
        for m in [1, 2, 32]:
            for rho in [0, 0.1, 0.5]:
                recomb_map = msprime.RecombinationMap.uniform_map(m, rho, num_loci=m)
                ts = msprime.simulate(
                    n, recombination_map=recomb_map, mutation_rate=0.1)
                yield ts
    for ts in get_bottleneck_examples():
        yield ts
    ts = msprime.simulate(30, length=20, recombination_rate=1)
    assert ts.num_trees > 1
    yield make_alternating_back_mutations(ts)


def get_bottleneck_examples():
    """
    Returns an iterator of example tree sequences with nonbinary
    trees.
    """
    bottlenecks = [
        msprime.SimpleBottleneck(0.01, proportion=0.05),
        msprime.SimpleBottleneck(0.02, proportion=0.25),
        msprime.SimpleBottleneck(0.03, proportion=1)]
    for n in [3, 10, 100]:
        ts = msprime.simulate(
            n, length=100, recombination_rate=1,
            demographic_events=bottlenecks,
            random_seed=n)
        yield ts


def get_back_mutation_examples():
    """
    Returns an iterator of example tree sequences with nonbinary
    trees.
    """
    ts = msprime.simulate(10, random_seed=1)
    yield make_alternating_back_mutations(ts)
    for ts in get_bottleneck_examples():
        yield make_alternating_back_mutations(ts)


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
        for mutation in t.mutations():
            j = t.get_num_tracked_leaves(mutation.node)
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
    if tree.is_leaf(node):
        s = "{0}:{1}".format(node + 1, branch_lengths[node])
    else:
        c1, c2 = tree.get_children(node)
        s1 = _build_newick(c1, root, tree, branch_lengths)
        s2 = _build_newick(c2, root, tree, branch_lengths)
        if node == root:
            # The root node is treated differently
            s = "({0},{1});".format(s1, s2)
        else:
            s = "({0},{1}):{2}".format(
                s1, s2, branch_lengths[node])
    return s


def simplify_tree_sequence(ts, samples):
    """
    Simple tree-by-tree algorithm to get a simplify of a tree sequence.
    """
    if len(samples) < 2:
        raise ValueError("Must have at least two samples")
    # TODO this algorithm is partially refactored to use the tables API.
    # It should be updated to properly support non binary mutations and
    # to also remove references to the CoalescenceRecord object. The algorithm
    # can also be clarified, as it's quite muddled at the moment what a 'record'
    # is.

    # TODO remove
    for site in ts.sites():
        assert site.ancestral_state == "0"
        for mutation in site.mutations:
            assert mutation.derived_state == "1"

    num_nodes = ts.get_num_nodes()
    active_records = {}
    new_records = []
    site_records = collections.defaultdict(list)
    for tree in ts.trees(tracked_leaves=samples):
        parent = [msprime.NULL_NODE for j in range(num_nodes)]
        children = collections.defaultdict(list)
        for leaf in samples:
            u = leaf
            v = tree.get_parent(u)
            while v != msprime.NULL_NODE:
                is_parent = (
                    tree.get_num_tracked_leaves(v) > tree.get_num_tracked_leaves(u))
                if parent[u] == msprime.NULL_NODE and is_parent:
                    parent[u] = v
                    children[v].append(u)
                    children[v].sort()
                    u = v
                v = tree.get_parent(v)
        removed = []
        for u, record in active_records.items():
            if u not in children:
                active_records[u][1] = tree.get_interval()[0]
                new_records.append(msprime.CoalescenceRecord(*record))
                removed.append(u)
        for u in removed:
            del active_records[u]
        for u, c in children.items():
            if u not in active_records:
                active_records[u] = list(tree.get_interval()) + [
                    u, tuple(c), tree.get_time(u), tree.get_population(u)]
            elif active_records[u][3] != tuple(c):
                active_records[u][1] = tree.get_interval()[0]
                new_records.append(msprime.CoalescenceRecord(*active_records[u]))
                active_records[u][0] = tree.get_interval()[0]
                active_records[u][3] = tuple(c)
        subset_root = samples[0]
        while parent[subset_root] != msprime.NULL_NODE:
            subset_root = parent[subset_root]
        # Now find the new nodes for all mutations that can be mapped back in
        for site in tree.sites():
            new_nodes = []
            for mut in site.mutations:
                stack = [mut.node]
                while not len(stack) == 0:
                    u = stack.pop()
                    if parent[u] != msprime.NULL_NODE or u in children:
                        if u != subset_root:
                            new_nodes.append(u)
                        break
                    stack.extend(tree.get_children(u))
            if len(new_nodes) > 0:
                site_records[site.position].extend(new_nodes)
    for record in active_records.values():
        record[1] = ts.get_sequence_length()
        new_records.append(msprime.CoalescenceRecord(*record))
    new_records.sort(key=lambda r: (r.time, r.node, r.left))

    # Now compress the nodes.
    node_map = [msprime.NULL_NODE for _ in range(num_nodes)]
    new_nodes = msprime.NodeTable(num_nodes)
    new_edgesets = msprime.EdgesetTable(ts.num_edgesets, 2 * ts.num_edgesets)
    for j, u in enumerate(samples):
        node_map[u] = j
        node = ts.node(u)
        new_nodes.add_row(flags=node.flags, time=node.time, population=node.population)
    next_node = len(samples)
    for record in new_records:
        for u in list(record.children) + [record.node]:
            if node_map[u] == msprime.NULL_NODE:
                node = ts.node(u)
                new_nodes.add_row(
                    flags=node.flags, time=node.time, population=node.population)
                node_map[u] = next_node
                next_node += 1
        children = tuple(sorted(node_map[c] for c in record.children))
        new_edgesets.add_row(
            left=record.left, right=record.right, parent=node_map[record.node],
            children=children)
    new_sites = msprime.SiteTable()
    new_mutations = msprime.MutationTable()
    for j, position in enumerate(sorted(site_records.keys())):
        # We must sort the mapped nodes by nonincreasing time.
        nodes = sorted(
            [node_map[u] for u in site_records[position]],
            key=lambda v: -ts.node(v).time)
        # TODO get the ancestral_state and derived_state properly.
        new_sites.add_row(position=position, ancestral_state="0")
        for node in nodes:
            new_mutations.add_row(site=j, node=node, derived_state="1")

    return msprime.load_tables(
        nodes=new_nodes, edgesets=new_edgesets, sites=new_sites,
        mutations=new_mutations)


def make_alternating_back_mutations(ts):
    """
    Returns a copy of the specified tree sequence with a sequence of
    alternating mutations along each path in each tree.
    """
    nodes = msprime.NodeTable()
    edgesets = msprime.EdgesetTable()
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()

    site = 0
    for tree in ts.trees():
        sites.add_row(position=tree.interval[0], ancestral_state="0")
        state = {tree.root: 0}
        for u in tree.nodes():
            if u != tree.root:
                state[u] = (state[u] + 1) % 2
            for v in tree.children(u):
                state[v] = state[u]
        del state[tree.root]
        # Ensure we have some variation in our samples.
        s = sum(state[u] for u in tree.leaves(tree.root))
        if s == 0 or s == tree.sample_size:
            del state[next(tree.leaves(tree.root))]
        site_mutations = sorted([(-tree.time(u), u) for u in state.keys()])
        for _, u in sorted(site_mutations):
            mutations.add_row(site, u, str(state[u]))
        site += 1
    ts.dump_tables(nodes=nodes, edgesets=edgesets)
    return msprime.load_tables(
        nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)


class TestHarmonicNumber(unittest.TestCase):
    """
    Tests for the harmonic number calculation.
    """

    def test_harmonic_number(self):
        def H(n):
            return sum(1 / k for k in range(1, n + 1))
        for n in range(10, 1000, 100):
            self.assertAlmostEqual(msprime.harmonic_number(n), H(n), 1)


class TestAlmostEqual(unittest.TestCase):
    """
    Simple tests to ensure that the almost_equal() method is sensible.
    """

    def test_defaults(self):
        eps = sys.float_info.epsilon
        equal = [
            (1, 1), (0, 0), (1 + eps, 1), (1, 1 - eps),
            (10.000000000001, 10.0)]
        for a, b in equal:
            self.assertAlmostEqual(a, b)
            self.assertTrue(msprime.almost_equal(a, b))

    def test_near_zero(self):
        eps = sys.float_info.epsilon
        equal = [(0, 0), (eps, 0), (0, -eps), (-eps, eps)]
        for a, b in equal:
            self.assertAlmostEqual(a, b)
            self.assertTrue(
                msprime.almost_equal(a, b, abs_tol=1e-9))
        not_equal = [(0, 0.0000001), (-0.0000001, 0)]
        for a, b in not_equal:
            self.assertNotAlmostEqual(a, b)
            self.assertFalse(
                msprime.almost_equal(a, b, abs_tol=1e-9))


class HighLevelTestCase(tests.MsprimeTestCase):
    """
    Superclass of tests on the high level interface.
    """
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="msp_hl_testcase_")
        self.temp_file = os.path.join(self.temp_dir, "generic")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

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
        self.assertLessEqual(len(used_nodes), 2 * st.get_sample_size() - 1)
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
                self.assertEqual(len(st.get_children(u)), 0)
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
        self.assertLessEqual(len(pi), 2 * st.get_sample_size() - 1)
        self.assertNotIn(st.get_root(), pi)
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
            self.assertAlmostEqual(l, length)
            self.assertGreaterEqual(l, 0)
            self.assertGreater(r, l)
            self.assertLessEqual(r, ts.get_sequence_length())
            length += r - l
            self.verify_sparse_tree(st1)
            num_trees += 1
        self.assertEqual(breakpoints, list(ts.breakpoints()))
        self.assertAlmostEqual(length, ts.get_sequence_length())
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
        self.assertAlmostEqual(pi1, pi3)
        self.assertGreaterEqual(pi1, 0.0)
        self.assertFalse(math.isnan(pi1))
        # Check for a subsample.
        samples = range(ts.get_sample_size() // 2 + 1)
        pi1 = ts.get_pairwise_diversity(samples)
        pi2 = simple_get_pairwise_diversity([haplotypes[j] for j in samples])
        pi3 = get_pairwise_diversity(ts, samples)
        self.assertAlmostEqual(pi1, pi2)
        self.assertAlmostEqual(pi1, pi3)
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
        j = 0
        for st in ts.trees():
            tree_mutations = list(st.mutations())
            self.assertEqual(st.get_num_mutations(), len(tree_mutations))
            all_tree_mutations.extend(tree_mutations)
            for mutation in tree_mutations:
                left, right = st.get_interval()
                self.assertTrue(left <= mutation.position < right)
                self.assertNotEqual(
                    st.get_parent(mutation.node), msprime.NULL_NODE)
                self.assertEqual(mutation.index, j)
                j += 1
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
                random_seed=1, population_configurations=[m1]).trees())
            st2 = next(msprime.simulate(
                random_seed=1, population_configurations=[m2]).trees())
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

    def test_nonbinary_cases(self):
        for ts in get_bottleneck_examples():
            self.verify_sparse_trees(ts)

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

    def verify_dump_load(self, tree_sequence):
        """
        Dump the tree sequence and verify we can load again from the same
        file.
        """
        tree_sequence.dump(self.temp_file)
        other = msprime.load(self.temp_file)
        records = list(tree_sequence.edgesets())
        other_records = list(other.edgesets())
        self.assertEqual(records, other_records)
        haplotypes = list(tree_sequence.haplotypes())
        other_haplotypes = list(other.haplotypes())
        self.assertEqual(haplotypes, other_haplotypes)

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
        for record in tree_sequence.nodes():
            if record.time > t:
                t = record.time
        self.assertEqual(sim.get_time(), t)
        self.assertGreater(sim.get_num_common_ancestor_events(), 0)
        self.assertGreaterEqual(sim.get_num_recombination_events(), 0)
        self.assertGreaterEqual(sim.get_total_num_migration_events(), 0)
        self.assertGreaterEqual(sim.get_num_multiple_recombination_events(), 0)
        self.verify_sparse_trees(tree_sequence)
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
        for bad_type in ["xd", None, 4.4]:
            self.assertRaises(
                TypeError, msprime.TreeSimulator, [(0, 0), (0, 0)], bad_type)
        self.assertRaises(
            ValueError, msprime.TreeSimulator, [], recomb_map)
        self.assertRaises(
            ValueError, msprime.TreeSimulator, [(0, 0)], recomb_map)


class TestVariantGenerator(HighLevelTestCase):
    """
    Tests the variants() method to ensure the output is consistent.
    """
    def get_tree_sequence(self):
        ts = msprime.simulate(
            10, length=10, recombination_rate=1, mutation_rate=10)
        self.assertGreater(ts.get_num_mutations(), 10)
        return ts

    def test_as_bytes(self):
        ts = self.get_tree_sequence()
        n = ts.get_sample_size()
        m = ts.get_num_mutations()
        A = np.zeros((m, n), dtype='u1')
        B = np.zeros((m, n), dtype='u1')
        for variant in ts.variants():
            A[variant.index] = variant.genotypes
        for variant in ts.variants(as_bytes=True):
            self.assertIsInstance(variant.genotypes, bytes)
            B[variant.index] = np.fromstring(variant.genotypes, np.uint8) - ord('0')
        self.assertTrue(np.all(A == B))
        bytes_variants = list(ts.variants(as_bytes=True))
        for j, variant in enumerate(bytes_variants):
            self.assertEqual(j, variant.index)
            row = np.fromstring(variant.genotypes, np.uint8) - ord('0')
            self.assertTrue(np.all(A[j] == row))

    def test_site_information(self):
        ts = self.get_tree_sequence()
        for site, variant in zip(ts.sites(), ts.variants()):
            self.assertEqual(site.position, variant.position)
            self.assertEqual(site, variant.site)

    def test_no_mutations(self):
        ts = msprime.simulate(10)
        self.assertEqual(ts.get_num_mutations(), 0)
        variants = list(ts.variants())
        self.assertEqual(len(variants), 0)

    def test_recurrent_mutations_over_leaves(self):
        ts = self.get_tree_sequence()
        num_sites = 5
        sites = [
            msprime.Site(
                index=j, ancestral_state="0",
                position=j * ts.sequence_length / num_sites,
                mutations=[
                    msprime.Mutation(site=j, node=u, derived_state="1")
                    for u in range(ts.sample_size)])
            for j in range(num_sites)]
        ts = ts.copy(sites)
        sites = list(ts.sites())
        variants = list(ts.variants(as_bytes=True))
        self.assertEqual(len(variants), num_sites)
        for site, variant in zip(sites, variants):
            self.assertEqual(site.position, variant.position)
            self.assertEqual(site, variant.site)
            self.assertEqual(site.index, variant.index)
            self.assertEqual(variant.genotypes, b'1' * ts.sample_size)
        # Now try without as_bytes
        for variant in ts.variants():
            self.assertTrue(np.all(variant.genotypes == np.ones(ts.sample_size)))

    def test_recurrent_mutations_errors(self):
        ts = self.get_tree_sequence()
        tree = next(ts.trees())
        for u in tree.nodes():
            for leaf in tree.leaves(u):
                if leaf != u:
                    site = msprime.Site(
                        index=0, ancestral_state="0", position=0, mutations=[
                            msprime.Mutation(site=0, derived_state="1", node=u),
                            msprime.Mutation(site=0, derived_state="1", node=leaf)])
            ts_new = ts.copy(sites=[site])
            self.assertRaises(_msprime.LibraryError, list, ts_new.variants())


class TestHaplotypeGenerator(HighLevelTestCase):
    """
    Tests the haplotype generation code.
    """

    def verify_haplotypes(self, n, haplotypes):
        """
        Verify that the specified set of haplotypes is consistent.
        """
        self.assertEqual(len(haplotypes), n)
        m = len(haplotypes[0])
        for h in haplotypes:
            self.assertEqual(len(h), m)
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

    def verify_tree_sequence(self, tree_sequence):
        n = tree_sequence.sample_size
        m = tree_sequence.num_sites
        haplotypes = list(tree_sequence.haplotypes())
        A = np.zeros((n, m), dtype='u1')
        B = np.zeros((n, m), dtype='u1')
        for j, h in enumerate(haplotypes):
            self.assertEqual(len(h), m)
            A[j] = np.fromstring(h, np.uint8) - ord('0')
        for variant in tree_sequence.variants():
            B[:, variant.index] = variant.genotypes
        self.assertTrue(np.all(A == B))
        self.verify_haplotypes(n, haplotypes)

    def verify_simulation(self, n, m, r, theta):
        """
        Verifies a simulation for the specified parameters.
        """
        recomb_map = msprime.RecombinationMap.uniform_map(m, r, m)
        tree_sequence = msprime.simulate(
            n, recombination_map=recomb_map, mutation_rate=theta)
        self.verify_tree_sequence(tree_sequence)

    def test_random_parameters(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 1000)
            r = random.random()
            theta = random.uniform(0, 2)
            self.verify_simulation(n, m, r, theta)

    def test_nonbinary_trees(self):
        for ts in get_bottleneck_examples():
            self.verify_tree_sequence(ts)

    def test_recurrent_mutations_over_leaves(self):
        for ts in get_bottleneck_examples():
            num_sites = 5
            sites = [
                msprime.Site(
                    index=j, ancestral_state="0",
                    position=j * ts.sequence_length / num_sites,
                    mutations=[
                        msprime.Mutation(site=j, node=u, derived_state="1")
                        for u in range(ts.sample_size)])
                for j in range(num_sites)]
            ts_new = ts.copy(sites)
            ones = "1" * num_sites
            for h in ts_new.haplotypes():
                self.assertEqual(ones, h)

    def test_recurrent_mutations_errors(self):
        for ts in get_bottleneck_examples():
            tree = next(ts.trees())
            for u in tree.children(tree.root):
                sites = [
                    msprime.Site(
                        index=0, position=0, ancestral_state="0",
                        mutations=[
                            msprime.Mutation(site=0, node=tree.root, derived_state="1"),
                            msprime.Mutation(site=0, node=u, derived_state="1"),
                        ])]
                ts_new = ts.copy(sites)
                self.assertRaises(_msprime.LibraryError, ts_new.haplotypes)

    def test_back_mutations(self):
        for ts in get_back_mutation_examples():
            self.verify_tree_sequence(ts)


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
            sim = msprime.simulator_factory(
                n, Ne=Ne, recombination_map=recomb_map)
            sim.run()
            tree_sequence = sim.get_tree_sequence()
            breakpoints = sim.get_breakpoints()
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

    def test_sparse_trees(self):
        for ts in get_example_tree_sequences():
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
        for ts in get_example_tree_sequences():
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
        for ts in get_example_tree_sequences():
            self.verify_tracked_leaves(ts)

    def test_trees_interface(self):
        ts = list(get_example_tree_sequences())[0]
        # The defaults should make sense and count leaves.
        # get_num_tracked_leaves
        for t in ts.trees():
            self.assertEqual(t.get_num_leaves(0), 1)
            self.assertEqual(t.get_num_tracked_leaves(0), 0)
            self.assertEqual(list(t.leaves(0)), [0])

        for t in ts.trees(leaf_counts=False):
            self.assertEqual(t.get_num_leaves(0), 1)
            self.assertRaises(RuntimeError, t.get_num_tracked_leaves, 0)
            self.assertEqual(list(t.leaves(0)), [0])

        for t in ts.trees(leaf_counts=True):
            self.assertEqual(t.get_num_leaves(0), 1)
            self.assertEqual(t.get_num_tracked_leaves(0), 0)
            self.assertEqual(list(t.leaves(0)), [0])

        for t in ts.trees(leaf_counts=True, tracked_leaves=[0]):
            self.assertEqual(t.get_num_leaves(0), 1)
            self.assertEqual(t.get_num_tracked_leaves(0), 1)
            self.assertEqual(list(t.leaves(0)), [0])

        for t in ts.trees(leaf_lists=True, leaf_counts=True):
            self.assertEqual(t.get_num_leaves(0), 1)
            self.assertEqual(t.get_num_tracked_leaves(0), 0)
            self.assertEqual(list(t.leaves(0)), [0])

        for t in ts.trees(leaf_lists=True, leaf_counts=False):
            self.assertEqual(t.get_num_leaves(0), 1)
            self.assertRaises(RuntimeError, t.get_num_tracked_leaves, 0)
            self.assertEqual(list(t.leaves(0)), [0])

        # This is a bit weird as we don't seem to actually execute the
        # method until it is iterated.
        self.assertRaises(
            ValueError, list, ts.trees(leaf_counts=False, tracked_leaves=[0]))

    @unittest.skip("pi on recurrent mutations")
    def test_get_pairwise_diversity(self):
        for ts in get_example_tree_sequences():
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
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            self.assertRaises(ValueError, ts.get_population, -1)
            self.assertRaises(ValueError, ts.get_population, n)
            self.assertRaises(ValueError, ts.get_population, n + 1)
            self.assertEqual(ts.get_population(0), 0)
            self.assertEqual(ts.get_population(n - 1), 0)

    def test_get_time(self):
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            self.assertRaises(ValueError, ts.get_time, -1)
            self.assertRaises(ValueError, ts.get_time, n)
            self.assertRaises(ValueError, ts.get_time, n + 1)
            self.assertEqual(ts.get_time(0), 0)
            self.assertEqual(ts.get_time(n - 1), 0)

    def test_get_samples(self):
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            samples = list(range(n))
            self.assertEqual(ts.get_samples(), samples)
            self.assertEqual(ts.get_samples(0), samples)
            self.assertEqual(ts.get_samples(msprime.NULL_POPULATION), [])
            self.assertEqual(ts.get_samples(1), [])

    def test_write_vcf_interface(self):
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            for bad_ploidy in [-1, 0, n + 1]:
                self.assertRaises(ValueError, ts.write_vcf, self.temp_file, bad_ploidy)

    def verify_simplify_topology(self, ts, sample):
        new_ts = ts.simplify(sample)
        sample_map = {k: j for j, k in enumerate(sample)}
        old_trees = ts.trees()
        old_tree = next(old_trees)
        self.assertGreaterEqual(ts.get_num_trees(), new_ts.get_num_trees())
        for new_tree in new_ts.trees():
            new_left, new_right = new_tree.get_interval()
            old_left, old_right = old_tree.get_interval()
            # Skip ahead on the old tree until new_left is within its interval
            while old_right <= new_left:
                old_tree = next(old_trees)
                old_left, old_right = old_tree.get_interval()
            # If the TMRCA of all pairs of samples is the same, then we have the
            # same information. We limit this to at most 500 pairs
            pairs = itertools.islice(itertools.combinations(sample, 2), 500)
            for pair in pairs:
                mapped_pair = [sample_map[u] for u in pair]
                mrca1 = old_tree.get_mrca(*pair)
                mrca2 = new_tree.get_mrca(*mapped_pair)
                self.assertEqual(old_tree.get_time(mrca1), new_tree.get_time(mrca2))
                self.assertEqual(
                    old_tree.get_population(mrca1), new_tree.get_population(mrca2))

    def verify_simplify_mutations(self, ts, sample):
        # Get the allele counts within the subset.
        allele_counts = {mut.position: 0 for mut in ts.mutations()}
        sample_map = {k: j for j, k in enumerate(sample)}
        leaves = {mut.position: [] for mut in ts.mutations()}
        for tree in ts.trees(tracked_leaves=sample):
            for site in tree.sites():
                for mut in site.mutations:
                    allele_counts[site.position] += tree.get_num_tracked_leaves(mut.node)
                    for u in tree.leaves(mut.node):
                        if u in sample_map:
                            leaves[site.position].append(sample_map[u])
        new_ts = ts.simplify(sample)
        self.assertLessEqual(new_ts.get_num_sites(), ts.get_num_sites())
        self.assertLessEqual(new_ts.get_num_mutations(), ts.get_num_mutations())
        for tree in new_ts.trees():
            for site in tree.sites():
                leaf_count = 0
                new_leaves = []
                for mut in site.mutations:
                    leaf_count += tree.get_num_leaves(mut.node)
                    new_leaves.extend(tree.leaves(mut.node))
                self.assertEqual(sorted(leaves[site.position]), sorted(new_leaves))
                self.assertEqual(leaf_count, allele_counts[site.position])

    def verify_simplify_equality(self, ts, sample):
        s1 = ts.simplify(sample)
        s2 = simplify_tree_sequence(ts, sample)
        self.assertEqual(list(s1.edgesets()), list(s2.edgesets()))
        self.assertEqual(list(s1.nodes()), list(s2.nodes()))
        self.assertEqual(list(s1.sites()), list(s2.sites()))
        self.assertEqual(list(s1.haplotypes()), list(s2.haplotypes()))
        self.assertEqual(
            list(s1.variants(as_bytes=True)), list(s2.variants(as_bytes=True)))

    def verify_simplify_variants(self, ts, sample):
        subset = ts.simplify(sample)
        s = np.array(sample)
        full_genotypes = np.empty((ts.num_sites, ts.sample_size))
        full_positions = np.empty(ts.num_sites)
        for variant in ts.variants():
            full_positions[variant.index] = variant.position
            full_genotypes[variant.index] = variant.genotypes
        subset_genotypes = np.empty((subset.num_sites, subset.sample_size))
        subset_positions = np.empty(subset.num_sites)
        for variant in subset.variants():
            subset_positions[variant.index] = variant.position
            subset_genotypes[variant.index] = variant.genotypes
        j = 0
        for sg, sp in zip(subset_genotypes, subset_positions):
            while full_positions[j] < sp:
                unique = np.unique(full_genotypes[j][s])
                self.assertEqual(unique.shape[0], 1)
                self.assertIn(unique[0], [0, 1])
                j += 1
            self.assertEqual(full_positions[j], sp)
            self.assertTrue(np.all(sg == full_genotypes[j][s]))
            j += 1
        while j < ts.num_sites:
            unique = np.unique(full_genotypes[j][s])
            self.assertEqual(unique.shape[0], 1)
            self.assertIn(unique[0], [0, 1])
            j += 1

    @unittest.skip("simplify with lots of back mutations")
    def test_simplify(self):
        num_mutations = 0
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            num_mutations += ts.get_num_mutations()
            if n > 2:
                sample_sizes = set([2, max(2, n // 2), n - 1])
                for k in sample_sizes:
                    subset = random.sample(range(ts.get_sample_size()), k)
                    self.verify_simplify_topology(ts, subset)
                    self.verify_simplify_mutations(ts, subset)
                    self.verify_simplify_equality(ts, subset)
                    self.verify_simplify_variants(ts, subset)
        self.assertGreater(num_mutations, 0)

    def test_simplify_bugs(self):
        prefix = "tests/data/simplify-bugs/"
        j = 1
        while True:
            nodes_file = os.path.join(prefix, "{:02d}-nodes.txt".format(j))
            if not os.path.exists(nodes_file):
                break
            edgesets_file = os.path.join(prefix, "{:02d}-edgesets.txt".format(j))
            sites_file = os.path.join(prefix, "{:02d}-sites.txt".format(j))
            mutations_file = os.path.join(prefix, "{:02d}-mutations.txt".format(j))
            with open(nodes_file) as nodes, \
                    open(edgesets_file) as edgesets,\
                    open(sites_file) as sites,\
                    open(mutations_file) as mutations:
                ts = msprime.load_text(
                    nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)
            samples = list(range(ts.sample_size))
            self.verify_simplify_equality(ts, samples)
            self.verify_simplify_topology(ts, samples)
            self.verify_simplify_mutations(ts, samples)
            self.verify_simplify_variants(ts, samples)
            j += 1
        self.assertGreater(j, 1)

    @unittest.skip("pi on recurrent mutations")
    def test_apis(self):
        for ts in get_example_tree_sequences():
            self.assertEqual(ts.get_ll_tree_sequence(), ts.ll_tree_sequence)
            self.assertEqual(ts.get_provenance(), ts.provenance)
            self.assertEqual(ts.get_sample_size(), ts.sample_size)
            self.assertEqual(ts.get_sequence_length(), ts.sequence_length)
            self.assertEqual(ts.num_edgesets, ts.num_records)
            self.assertEqual(ts.get_num_trees(), ts.num_trees)
            self.assertEqual(ts.get_num_mutations(), ts.num_mutations)
            self.assertEqual(ts.get_num_nodes(), ts.num_nodes)
            self.assertEqual(
                ts.get_pairwise_diversity(), ts.pairwise_diversity())
            samples = range(ts.get_sample_size() // 2 + 1)
            self.assertEqual(
                ts.get_pairwise_diversity(samples), ts.pairwise_diversity(samples))
            for s in samples:
                self.assertEqual(ts.get_time(s), ts.time(s))
                p = ts.get_population(s)
                self.assertEqual(p, ts.population(s))
                self.assertEqual(ts.get_samples(p), ts.samples(p))
            self.assertEqual(ts.get_samples(), ts.samples())

    def test_copy(self):
        for ts1 in get_example_tree_sequences():
            ts2 = ts1.copy()
            self.assertNotEqual(id(ts1), id(ts2))
            self.assertEqual(list(ts1.edgesets()), list(ts2.edgesets()))
            self.assertEqual(list(ts1.nodes()), list(ts2.nodes()))
            self.assertEqual(list(ts1.mutations()), list(ts2.mutations()))
            site_lists = [[], list(ts1.sites())[:-1]]
            for sites in site_lists:
                ts2 = ts1.copy(sites=sites)
                self.assertNotEqual(id(ts1), id(ts2))
                self.assertEqual(list(ts1.edgesets()), list(ts2.edgesets()))
                self.assertEqual(list(ts1.nodes()), list(ts2.nodes()))
                self.assertEqual(sites, list(ts2.sites()))

    def test_generate_mutations_on_tree_sequence(self):
        some_mutations = False
        for ts in get_example_tree_sequences():
            nodes = msprime.NodeTable()
            edgesets = msprime.EdgesetTable()
            sites = msprime.SiteTable()
            mutations = msprime.MutationTable()
            ts.dump_tables(nodes=nodes, edgesets=edgesets)
            mutgen = msprime.MutationGenerator(msprime.RandomGenerator(1), 10)
            mutgen.generate(nodes, edgesets, sites, mutations)
            if mutations.num_rows > 0:
                some_mutations = True
            tsp = msprime.load_tables(
                nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)
            self.assertEqual(tsp.num_mutations, mutations.num_rows)
        self.assertTrue(some_mutations)

    def test_sites(self):
        some_sites = False
        for ts in get_example_tree_sequences():
            previous_pos = -1
            for index, site in enumerate(ts.sites()):
                self.assertTrue(0 <= site.position < ts.sequence_length)
                self.assertGreater(site.position, previous_pos)
                self.assertGreater(len(site.mutations), 0)
                self.assertEqual(site.index, index)
                self.assertEqual(site.ancestral_state, '0')
                previous_pos = site.position
                self.assertGreater(len(site.mutations), 0)
                for mutation in site.mutations:
                    self.assertEqual(mutation.site, site.index)
                    self.assertIn(mutation.derived_state, ['0', '1'])
                    self.assertTrue(0 <= mutation.node < ts.num_nodes)
                some_sites = True
        self.assertTrue(some_sites)

    def test_sites_mutations(self):
        # Check that the mutations iterator returns the correct values.
        for ts in get_example_tree_sequences():
            mutations = list(ts.mutations())
            other_mutations = []
            for site in ts.sites():
                for mut in site.mutations:
                    other_mutations.append(msprime.DeprecatedMutation(
                        position=site.position, node=mut.node, index=site.index))
            self.assertEqual(mutations, other_mutations)


class TestTreeSequenceTextIO(HighLevelTestCase):
    """
    Tests for the tree sequence text IO.
    """

    def verify_nodes_format(self, ts, nodes_file, precision):
        """
        Verifies that the nodes we output have the correct form.
        """
        def convert(v):
            return "{:.{}f}".format(v, precision)
        output_nodes = nodes_file.read().splitlines()
        self.assertEqual(len(output_nodes) - 1, ts.num_nodes)
        self.assertEqual(
            list(output_nodes[0].split()),
            ["is_sample", "time", "population"])
        for node, line in zip(ts.nodes(), output_nodes[1:]):
            splits = line.split("\t")
            self.assertEqual(str(node.is_sample()), splits[0])
            self.assertEqual(convert(node.time), splits[1])
            self.assertEqual(str(node.population), splits[2])

    def verify_edgesets_format(self, ts, edgesets_file, precision):
        """
        Verifies that the edgesets we output have the correct form.
        """
        def convert(v):
            return "{:.{}f}".format(v, precision)
        output_edgesets = edgesets_file.read().splitlines()
        self.assertEqual(len(output_edgesets) - 1, ts.num_edgesets)
        self.assertEqual(
            list(output_edgesets[0].split()),
            ["left", "right", "parent", "children"])
        for edgeset, line in zip(ts.edgesets(), output_edgesets[1:]):
            splits = line.split("\t")
            self.assertEqual(convert(edgeset.left), splits[0])
            self.assertEqual(convert(edgeset.right), splits[1])
            self.assertEqual(str(edgeset.parent), splits[2])
            self.assertEqual(",".join(map(str, edgeset.children)), splits[3])

    def verify_sites_format(self, ts, sites_file, precision):
        """
        Verifies that the sites we output have the correct form.
        """
        def convert(v):
            return "{:.{}f}".format(v, precision)
        output_sites = sites_file.read().splitlines()
        self.assertEqual(len(output_sites) - 1, ts.num_sites)
        self.assertEqual(
            list(output_sites[0].split()),
            ["position", "ancestral_state"])
        for site, line in zip(ts.sites(), output_sites[1:]):
            splits = line.split("\t")
            self.assertEqual(convert(site.position), splits[0])

    def verify_mutations_format(self, ts, mutations_file, precision):
        """
        Verifies that the mutationss we output have the correct form.
        """
        def convert(v):
            return "{:.{}f}".format(v, precision)
        output_mutations = mutations_file.read().splitlines()
        self.assertEqual(len(output_mutations) - 1, ts.num_mutations)
        self.assertEqual(
            list(output_mutations[0].split()),
            ["site", "node", "derived_state"])
        mutations = [mut for site in ts.sites() for mut in site.mutations]
        for mutation, line in zip(mutations, output_mutations[1:]):
            splits = line.split("\t")
            self.assertEqual(str(mutation.site), splits[0])
            self.assertEqual(str(mutation.node), splits[1])
            self.assertEqual(str(mutation.derived_state), splits[2])

    def test_output_format(self):
        for ts in get_example_tree_sequences():
            for precision in [2, 7]:
                nodes_file = six.StringIO()
                edgesets_file = six.StringIO()
                sites_file = six.StringIO()
                mutations_file = six.StringIO()
                ts.dump_text(
                    nodes=nodes_file, edgesets=edgesets_file, sites=sites_file,
                    mutations=mutations_file, precision=precision)
                nodes_file.seek(0)
                edgesets_file.seek(0)
                sites_file.seek(0)
                mutations_file.seek(0)
                self.verify_nodes_format(ts, nodes_file, precision)
                self.verify_edgesets_format(ts, edgesets_file, precision)
                self.verify_sites_format(ts, sites_file, precision)
                self.verify_mutations_format(ts, mutations_file, precision)

    def verify_approximate_equality(self, ts1, ts2):
        """
        Verifies that the specified tree sequences are approximately
        equal, taking into account the error incurred in exporting to text.
        """
        self.assertEqual(ts1.sample_size, ts2.sample_size)
        self.assertAlmostEqual(ts1.sequence_length, ts2.sequence_length)
        self.assertEqual(ts1.num_nodes, ts2.num_nodes)
        self.assertEqual(ts1.num_edgesets, ts2.num_edgesets)
        self.assertEqual(ts1.num_sites, ts2.num_sites)
        self.assertEqual(ts1.num_mutations, ts2.num_mutations)

        checked = 0
        for n1, n2 in zip(ts1.nodes(), ts2.nodes()):
            self.assertEqual(n1.population, n2.population)
            self.assertEqual(n1.name, n2.name)
            self.assertAlmostEqual(n1.time, n2.time)
            checked += 1
        self.assertEqual(checked, ts1.num_nodes)

        checked = 0
        for r1, r2 in zip(ts1.edgesets(), ts2.edgesets()):
            checked += 1
            self.assertAlmostEqual(r1.left, r2.left)
            self.assertAlmostEqual(r1.right, r2.right)
            self.assertEqual(r1.parent, r2.parent)
            self.assertEqual(r1.children, r2.children)
        self.assertEqual(ts1.num_edgesets, checked)

        checked = 0
        for s1, s2 in zip(ts1.sites(), ts2.sites()):
            checked += 1
            self.assertAlmostEqual(s1.position, s2.position)
            self.assertAlmostEqual(s1.ancestral_state, s2.ancestral_state)
            self.assertEqual(s1.mutations, s2.mutations)
        self.assertEqual(ts1.num_sites, checked)

        # Check the trees
        check = 0
        for t1, t2 in zip(ts1.trees(), ts2.trees()):
            self.assertEqual(list(t1.nodes()), list(t2.nodes()))
            check += 1
        self.assertEqual(check, ts1.get_num_trees())

    def test_text_record_round_trip(self):
        for ts1 in get_example_tree_sequences():
            nodes_file = six.StringIO()
            edgesets_file = six.StringIO()
            sites_file = six.StringIO()
            mutations_file = six.StringIO()
            ts1.dump_text(
                nodes=nodes_file, edgesets=edgesets_file, sites=sites_file,
                mutations=mutations_file, precision=9)
            nodes_file.seek(0)
            edgesets_file.seek(0)
            sites_file.seek(0)
            mutations_file.seek(0)
            ts2 = msprime.load_text(
                nodes=nodes_file, edgesets=edgesets_file, sites=sites_file,
                mutations=mutations_file)
            self.verify_approximate_equality(ts1, ts2)

    def test_empty_files(self):
        nodes_file = six.StringIO("is_sample\ttime\n")
        edgesets_file = six.StringIO("left\tright\tparent\tchildren\n")
        sites_file = six.StringIO("position\tancestral_state\n")
        mutations_file = six.StringIO("site\tnode\tderived_state\n")
        ts = msprime.load_text(
            nodes=nodes_file, edgesets=edgesets_file, sites=sites_file,
            mutations=mutations_file)
        self.assertEqual(ts.num_nodes, 0)
        self.assertEqual(ts.num_edgesets, 0)
        self.assertEqual(ts.num_sites, 0)
        self.assertEqual(ts.num_mutations, 0)


class TestSparseTree(HighLevelTestCase):
    """
    Some simple tests on the API for the sparse tree.
    """
    def get_tree(self, leaf_lists=False):
        ts = msprime.simulate(10, random_seed=1, mutation_rate=1)
        return next(ts.trees(leaf_lists=leaf_lists))

    def test_str(self):
        t = self.get_tree()
        self.assertIsInstance(str(t), str)
        self.assertEqual(str(t), str(t.get_parent_dict()))

    def test_leaves(self):
        for leaf_lists in [True, False]:
            t = self.get_tree(leaf_lists)
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
        w = 123
        h = 456
        t.draw(self.temp_file, w, h, show_times=True)
        self.assertGreater(os.path.getsize(self.temp_file), 0)
        with open(self.temp_file) as f:
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
        orders = ["inorder", "postorder", "levelorder", "breadthfirst"]
        for test_order in orders:
            self.assertEqual(
                list(t1.nodes(order=test_order)),
                list(t1.nodes(t1.get_root(), order=test_order)))
            self.assertEqual(
                list(t1.nodes(order=test_order)),
                list(t1.nodes(t1.get_root(), test_order)))
            self.assertEqual(
               list(t1.nodes(order=test_order)),
               list(t2.nodes(order=test_order)))
            for u in t1.nodes():
                self.assertEqual(
                    list(t1.nodes(u, test_order)),
                    list(t2.nodes(u, test_order)))
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

    def test_apis(self):
        # tree properties
        t1 = self.get_tree()
        self.assertEqual(t1.get_root(), t1.root)
        self.assertEqual(t1.get_index(), t1.index)
        self.assertEqual(t1.get_interval(), t1.interval)
        self.assertEqual(t1.get_length(), t1.length)
        self.assertEqual(t1.get_sample_size(), t1.sample_size)
        self.assertEqual(t1.get_num_mutations(), t1.num_mutations)
        self.assertEqual(t1.get_parent_dict(), t1.parent_dict)
        self.assertEqual(t1.get_time_dict(), t1.time_dict)
        self.assertEqual(t1.get_total_branch_length(), t1.total_branch_length)
        # node properties
        root = t1.get_root()
        for node in t1.nodes():
            if node != root:
                self.assertEqual(t1.get_time(node), t1.time(node))
                self.assertEqual(t1.get_parent(node), t1.parent(node))
                self.assertEqual(t1.get_children(node), t1.children(node))
                self.assertEqual(t1.get_population(node), t1.population(node))
                self.assertEqual(t1.get_num_leaves(node), t1.num_leaves(node))
                self.assertEqual(t1.get_branch_length(node),
                                 t1.branch_length(node))
                self.assertEqual(t1.get_num_tracked_leaves(node),
                                 t1.num_tracked_leaves(node))

        pairs = itertools.islice(itertools.combinations(t1.nodes(), 2), 50)
        for pair in pairs:
            self.assertEqual(t1.get_mrca(*pair), t1.mrca(*pair))
            self.assertEqual(t1.get_tmrca(*pair), t1.tmrca(*pair))


class TestRecombinationMap(HighLevelTestCase):
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

    def test_read_hapmap_simple(self):
        with open(self.temp_file, "w+") as f:
            print("HEADER", file=f)
            print("chr1 0 1", file=f)
            print("chr1 1 5 x", file=f)
            print("s    2 0 x x x", file=f)
        rm = msprime.RecombinationMap.read_hapmap(self.temp_file)
        self.assertEqual(rm.get_positions(), [0, 1, 2])
        self.assertEqual(rm.get_rates(), [1e-8, 5e-8, 0])

    def test_read_hapmap_nonzero_start(self):
        with open(self.temp_file, "w+") as f:
            print("HEADER", file=f)
            print("chr1 1 5 x", file=f)
            print("s    2 0 x x x", file=f)
        rm = msprime.RecombinationMap.read_hapmap(self.temp_file)
        self.assertEqual(rm.get_positions(), [0, 1, 2])
        self.assertEqual(rm.get_rates(), [0, 5e-8, 0])

    def test_read_hapmap_gzipped(self):
        try:
            filename = self.temp_file + ".gz"
            with gzip.open(filename, "w+") as f:
                f.write(b"HEADER\n")
                f.write(b"chr1 0 1\n")
                f.write(b"chr1 1 5.5\n")
                f.write(b"s    2 0\n")
            rm = msprime.RecombinationMap.read_hapmap(filename)
            self.assertEqual(rm.get_positions(), [0, 1, 2])
            self.assertEqual(rm.get_rates(), [1e-8, 5.5e-8, 0])
        finally:
            os.unlink(filename)


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
            samples = ll_sim.get_samples()
            self.assertEqual(len(samples), n)
            for sample in samples:
                self.assertEqual(sample[0], 0)
                self.assertEqual(sample[1], 0)

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
            configs = [msprime.PopulationConfiguration(2) for _ in range(d)]
            sim = msprime.simulator_factory(population_configurations=configs)
            self.assertEqual(sim.get_sample_size(), 2 * d)
            samples = []
            for j in range(d):
                samples += [
                    msprime.Sample(population=j, time=0) for _ in range(2)]
            self.assertEqual(sim.get_samples(), samples)
            ll_sim = sim.create_ll_instance()
            self.assertEqual(ll_sim.get_samples(), samples)

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
                population_configurations=pop_configs)
            ll_sim = sim.create_ll_instance()
            # If we don't specify a matrix, it's 0 everywhere.
            matrix = [0 for j in range(N * N)]
            self.assertEqual(ll_sim.get_migration_matrix(), matrix)

            def f(hl_matrix):
                return msprime.simulator_factory(
                    population_configurations=pop_configs,
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

    def test_sample_combination_errors(self):
        # Make sure that the various ways we can specify the samples
        # operate correctly.
        s = msprime.Sample(time=0.0, population=0)
        self.assertRaises(ValueError, msprime.simulator_factory)
        # Cannot provide sample_size with either population configurations
        # or samples
        self.assertRaises(
            ValueError, msprime.simulator_factory,
            sample_size=2, samples=[s, s])
        pop_configs = [
            msprime.PopulationConfiguration(sample_size=2)]
        self.assertRaises(
            ValueError, msprime.simulator_factory,
            sample_size=2, population_configurations=pop_configs)
        # If we provide samples and population_configurations we cannot
        # have a sample size for the config.
        pop_configs = [
            msprime.PopulationConfiguration(sample_size=2)]
        self.assertRaises(
            ValueError, msprime.simulator_factory,
            samples=[s, s], population_configurations=pop_configs)
        pop_configs = [
            msprime.PopulationConfiguration(sample_size=None),
            msprime.PopulationConfiguration(sample_size=2)]
        self.assertRaises(
            ValueError, msprime.simulator_factory,
            samples=[s, s], population_configurations=pop_configs)

    def test_samples(self):
        pop_configs = [
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration()]
        samples = [
            msprime.Sample(population=0, time=0),
            msprime.Sample(population=1, time=1),
            msprime.Sample(population=2, time=2)]
        # Ne = 1/4 to keep in coalescence units.
        sim = msprime.simulator_factory(
            Ne=1/4, samples=samples, population_configurations=pop_configs)
        self.assertEqual(sim.get_samples(), samples)
        ll_sim = sim.create_ll_instance()
        self.assertEqual(ll_sim.get_samples(), samples)


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
        self.assertEqual(ts.get_sequence_length(), 1)
        self.assertEqual(len(ts.provenance), 1)

    def test_provenance(self):
        ts = msprime.simulate(10)
        self.assertEqual(len(ts.provenance), 1)
        d = json.loads(ts.provenance[0].decode())
        # TODO check the form of the dictionary
        for ts in msprime.simulate(10, num_replicates=10):
            self.assertEqual(len(ts.provenance), 1)
            d = json.loads(ts.provenance[0].decode())
            self.assertGreater(len(d), 0)

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

    def test_mutation_generator(self):
        n = 10
        rng = msprime.RandomGenerator(1)
        mutgen = msprime.MutationGenerator(rng, 10)
        ts = msprime.simulate(n, mutation_generator=mutgen)
        self.assertIsInstance(ts, msprime.TreeSequence)
        self.assertEqual(ts.get_sample_size(), n)
        self.assertEqual(ts.get_num_trees(), 1)
        self.assertGreater(ts.get_num_mutations(), 0)

    def test_mutation_interface(self):
        for bad_type in ["x", [], {}]:
            self.assertRaises(
                TypeError, msprime.simulate, 10, mutation_rate=bad_type)
        mutgen = msprime.MutationGenerator(msprime.RandomGenerator(1), 1)
        self.assertRaises(
            ValueError, msprime.simulate, 10, mutation_generator=mutgen,
            mutation_rate=1)

    def test_recombination(self):
        n = 10
        ts = msprime.simulate(n, recombination_rate=10)
        self.assertIsInstance(ts, msprime.TreeSequence)
        self.assertEqual(ts.get_sample_size(), n)
        self.assertGreater(ts.get_num_trees(), 1)
        self.assertEqual(ts.get_num_mutations(), 0)


class TestNodeOrdering(HighLevelTestCase):
    """
    Verify that we can use any node ordering for internal nodes
    and get the same topologies.
    """
    num_random_permutations = 10

    def verify_tree_sequences_equal(self, ts1, ts2, approx=False):
        self.assertEqual(ts1.get_num_trees(), ts2.get_num_trees())
        self.assertEqual(ts1.get_sample_size(), ts2.get_sample_size())
        self.assertEqual(ts1.get_num_nodes(), ts2.get_num_nodes())
        j = 0
        for r1, r2 in zip(ts1.edgesets(), ts2.edgesets()):
            self.assertEqual(r1.parent, r2.parent)
            self.assertEqual(r1.children, r2.children)
            if approx:
                self.assertAlmostEqual(r1.left, r2.left)
                self.assertAlmostEqual(r1.right, r2.right)
            else:
                self.assertEqual(r1.left, r2.left)
                self.assertEqual(r1.right, r2.right)
            j += 1
        self.assertEqual(ts1.num_edgesets, j)
        j = 0
        for n1, n2 in zip(ts1.nodes(), ts2.nodes()):
            self.assertEqual(n1.name, n2.name)
            self.assertEqual(n1.population, n2.population)
            if approx:
                self.assertAlmostEqual(n1.time, n2.time)
            else:
                self.assertEqual(n1.time, n2.time)
            j += 1
        self.assertEqual(ts1.num_nodes, j)

    def verify_random_permutation(self, ts):
        n = ts.sample_size
        node_map = {}
        for j in range(n):
            node_map[j] = j
        internal_nodes = list(range(n, ts.num_nodes))
        random.shuffle(internal_nodes)
        for j, node in enumerate(internal_nodes):
            node_map[n + j] = node
        node_table = msprime.NodeTable()
        # Insert the new nodes into the table.
        inv_node_map = {v: k for k, v in node_map.items()}
        for j in range(ts.num_nodes):
            node = ts.node(inv_node_map[j])
            node_table.add_row(
                flags=node.flags, time=node.time, population=node.population)
        edgeset_table = msprime.EdgesetTable()
        for e in ts.edgesets():
            edgeset_table.add_row(
                left=e.left, right=e.right, parent=node_map[e.parent],
                children=tuple(
                    sorted([node_map[c] for c in e.children])))
        other_ts = msprime.load_tables(nodes=node_table, edgesets=edgeset_table)

        self.assertEqual(ts.get_num_trees(), other_ts.get_num_trees())
        self.assertEqual(ts.get_sample_size(), other_ts.get_sample_size())
        self.assertEqual(ts.get_num_nodes(), other_ts.get_num_nodes())
        j = 0
        for t1, t2 in zip(ts.trees(), other_ts.trees()):
            # Verify the topologies are identical. We do this by traversing
            # upwards to the root for every leaf and checking if we map to
            # the correct node and time.
            for u in range(n):
                v_orig = u
                v_map = u
                while v_orig != msprime.NULL_NODE:
                    self.assertEqual(node_map[v_orig], v_map)
                    self.assertEqual(
                        t1.get_time(v_orig),
                        t2.get_time(v_map))
                    v_orig = t1.get_parent(v_orig)
                    v_map = t2.get_parent(v_map)
                self.assertEqual(v_orig, msprime.NULL_NODE)
                self.assertEqual(v_map, msprime.NULL_NODE)
            j += 1
        self.assertEqual(j, ts.get_num_trees())
        # Verify we can dump this new tree sequence OK.
        other_ts.dump(self.temp_file)
        ts3 = msprime.load(self.temp_file)
        self.verify_tree_sequences_equal(other_ts, ts3)
        nodes_file = six.StringIO()
        edgesets_file = six.StringIO()
        # Also verify we can read the text version.
        other_ts.dump_text(nodes=nodes_file, edgesets=edgesets_file, precision=14)
        nodes_file.seek(0)
        edgesets_file.seek(0)
        ts3 = msprime.load_text(nodes_file, edgesets_file)
        self.verify_tree_sequences_equal(other_ts, ts3, True)

    def test_single_locus(self):
        ts = msprime.simulate(7)
        for _ in range(self.num_random_permutations):
            self.verify_random_permutation(ts)

    def test_multi_locus(self):
        ts = msprime.simulate(20, recombination_rate=10)
        for _ in range(self.num_random_permutations):
            self.verify_random_permutation(ts)

    def test_nonbinary(self):
        ts = msprime.simulate(
            sample_size=20, recombination_rate=10,
            demographic_events=[
                msprime.SimpleBottleneck(time=0.5, proportion=1)])
        # Make sure this really has some non-binary nodes
        found = False
        for r in ts.edgesets():
            if len(r.children) > 2:
                found = True
        self.assertTrue(found)
        for _ in range(self.num_random_permutations):
            self.verify_random_permutation(ts)
