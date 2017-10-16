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
import six
import sys
import tempfile
import unittest

import numpy as np

import msprime
import _msprime
import tests
import tests.tsutil as tsutil


def get_uniform_mutations(num_mutations, sequence_length, nodes):
    """
    Returns n evenly mutations over the specified list of nodes.
    """
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    for j in range(num_mutations):
        sites.add_row(
            position=j * (sequence_length / num_mutations), ancestral_state='0')
        mutations.add_row(site=j, derived_state='1', node=nodes[j % len(nodes)])
    return sites, mutations


def insert_gap(ts, position, length):
    """
    Inserts a gap of the specified size into the specified tree sequence.
    This involves: (1) breaking all edges that intersect with this point;
    and (2) shifting all coordinates greater than this value up by the
    gap length.
    """
    new_edges = []
    for e in ts.edges():
        if e.left < position < e.right:
            new_edges.append([e.left, position, e.parent, e.child])
            new_edges.append([position, e.right, e.parent, e.child])
        else:
            new_edges.append([e.left, e.right, e.parent, e.child])

    # Now shift up all coordinates.
    for e in new_edges:
        # Left coordinates == position get shifted
        if e[0] >= position:
            e[0] += length
        # Right coordinates == position do not get shifted
        if e[1] > position:
            e[1] += length
    tables = ts.dump_tables()
    edges = msprime.EdgeTable()
    for left, right, parent, child in new_edges:
        edges.add_row(left, right, parent, child)
    msprime.sort_tables(nodes=tables.nodes, edges=edges)
    # Throw in a bunch of mutations over the whole sequence on the samples.
    L = ts.sequence_length + length
    sites, mutations = get_uniform_mutations(100, L, list(ts.samples()))
    return msprime.load_tables(
            nodes=tables.nodes, edges=edges, sites=sites, mutations=mutations)


def get_gap_examples():
    """
    Returns example tree sequences that contain gaps within the list of
    edges.
    """
    ts = msprime.simulate(20, random_seed=56, recombination_rate=1)

    assert ts.num_trees > 1

    gap = 0.0125
    for x in [0, 0.1, 0.5, 0.75]:
        ts = insert_gap(ts, x, gap)
        found = False
        for t in ts.trees():
            if t.interval[0] == x:
                assert t.interval[1] == x + gap
                assert len(t.parent_dict) == 0
                found = True
        assert found
        yield ts
    # Give an example with a gap at the end.
    ts = msprime.simulate(10, random_seed=5, recombination_rate=1)
    t = ts.dump_tables()
    L = 2
    sites, mutations = get_uniform_mutations(100, L, list(ts.samples()))
    ts_new = ts.load_tables(
        nodes=t.nodes, edges=t.edges, sites=sites, mutations=mutations,
        sequence_length=L)
    yield ts_new


def get_internal_samples_examples():
    """
    Returns example tree sequences with internal samples.
    """
    n = 5
    ts = msprime.simulate(n, random_seed=10, mutation_rate=5)
    assert ts.num_mutations > 0
    tables = ts.dump_tables()
    nodes = tables.nodes
    flags = nodes.flags
    # Set all nodes to be samples.
    flags[:] = msprime.NODE_IS_SAMPLE
    nodes.set_columns(flags=flags, time=nodes.time, population=nodes.population)
    ts = msprime.load_tables(
        nodes=nodes, edges=tables.edges,
        sites=tables.sites, mutations=tables.mutations)
    yield ts
    # Set just internal nodes to be samples.
    flags[:] = 0
    flags[n:] = msprime.NODE_IS_SAMPLE
    nodes.set_columns(flags=flags, time=nodes.time, population=nodes.population)
    ts = msprime.load_tables(
        nodes=nodes, edges=tables.edges,
        sites=tables.sites, mutations=tables.mutations)
    yield ts
    # Set a mixture of internal and leaf samples.
    flags[:] = 0
    flags[n // 2: n + n // 2] = msprime.NODE_IS_SAMPLE
    nodes.set_columns(flags=flags, time=nodes.time, population=nodes.population)
    ts = msprime.load_tables(
        nodes=nodes, edges=tables.edges,
        sites=tables.sites, mutations=tables.mutations)
    yield ts


def get_decapitated_examples():
    """
    Returns example tree sequences in which the oldest edges have been removed.
    """
    ts = msprime.simulate(10, random_seed=1234)
    yield tsutil.decapitate(ts, ts.num_edges // 2)

    ts = msprime.simulate(20, recombination_rate=1, random_seed=1234)
    assert ts.num_trees > 2
    yield tsutil.decapitate(ts, ts.num_edges // 4)


def get_example_tree_sequences(back_mutations=True, gaps=True, internal_samples=True):
    if gaps:
        for ts in get_decapitated_examples():
            yield ts
        for ts in get_gap_examples():
            yield ts
    if internal_samples:
        for ts in get_internal_samples_examples():
            yield ts
    for n in [2, 3, 10, 100]:
        for m in [1, 2, 32]:
            for rho in [0, 0.1, 0.5]:
                recomb_map = msprime.RecombinationMap.uniform_map(m, rho, num_loci=m)
                ts = msprime.simulate(
                    n, recombination_map=recomb_map, mutation_rate=0.1)
                yield ts
    for ts in get_bottleneck_examples():
        yield ts
    ts = msprime.simulate(15, length=4, recombination_rate=1)
    assert ts.num_trees > 1
    if back_mutations:
        yield tsutil.insert_branch_mutations(ts)


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
    yield tsutil.insert_branch_mutations(ts)
    for ts in get_bottleneck_examples():
        yield tsutil.insert_branch_mutations(ts)


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
        tracked_samples = tree_sequence.get_samples()
    else:
        tracked_samples = list(samples)
    if len(tracked_samples) < 2:
        raise ValueError("len(samples) must be >= 2")
    pi = 0
    k = len(tracked_samples)
    denom = k * (k - 1) / 2
    for t in tree_sequence.trees(tracked_samples=tracked_samples):
        for mutation in t.mutations():
            j = t.get_num_tracked_samples(mutation.node)
            pi += j * (k - j) / denom
    return pi


def simplify_tree_sequence(ts, samples, filter_invariant_sites=True):
    """
    Simple tree-by-tree algorithm to get a simplify of a tree sequence.
    """
    s = tests.Simplifier(ts, samples, filter_invariant_sites=filter_invariant_sites)
    return s.simplify()


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
        oriented_forest = [st.get_parent(j) for j in range(st.num_nodes)]
        mrca_calc = tests.MRCACalculator(oriented_forest)
        # We've done exhaustive tests elsewhere, no need to go
        # through the combinations.
        for j in range(st.num_nodes):
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
        roots = set()
        for u in st.samples():
            # verify the path to root
            self.assertTrue(st.is_sample(u))
            times = []
            while st.get_parent(u) != msprime.NULL_NODE:
                v = st.get_parent(u)
                times.append(st.get_time(v))
                self.assertGreaterEqual(st.get_time(v), 0.0)
                self.assertIn(u, st.get_children(v))
                u = v
            roots.add(u)
            self.assertEqual(times, sorted(times))
        self.assertEqual(sorted(list(roots)), sorted(st.roots))
        self.assertEqual(len(st.roots), st.num_roots)
        u = st.left_root
        roots = []
        while u != msprime.NULL_NODE:
            roots.append(u)
            u = st.right_sib(u)
        self.assertEqual(roots, st.roots)
        # To a top-down traversal, and make sure we meet all the samples.
        samples = []
        for root in st.roots:
            stack = [root]
            while len(stack) > 0:
                u = stack.pop()
                self.assertNotEqual(u, msprime.NULL_NODE)
                if st.is_sample(u):
                    samples.append(u)
                if st.is_leaf(u):
                    self.assertEqual(len(st.get_children(u)), 0)
                else:
                    for c in reversed(st.get_children(u)):
                        stack.append(c)
                # Check that we get the correct number of samples at each
                # node.
                self.assertEqual(st.get_num_samples(u), len(list(st.samples(u))))
                self.assertEqual(st.get_num_tracked_samples(u), 0)
        self.assertEqual(sorted(samples), sorted(st.samples()))
        # Check the parent dict
        pi = st.get_parent_dict()
        for root in st.roots:
            self.assertNotIn(root, pi)
        for k, v in pi.items():
            self.assertEqual(st.get_parent(k), v)
        self.assertEqual(st.num_samples(), len(samples))
        self.assertEqual(sorted(st.samples()), sorted(samples))

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
            roots = set()
            for u in ts.samples():
                root = u
                while st1.get_parent(root) != msprime.NULL_NODE:
                    root = st1.get_parent(root)
                roots.add(root)
            self.assertEqual(sorted(list(roots)), sorted(st1.roots))
            if len(roots) > 1:
                with self.assertRaises(ValueError):
                    st1.root
            else:
                self.assertEqual(st1.root, list(roots)[0])
            self.assertEqual(st2, st1)
            self.assertFalse(st2 != st1)
            l, r = st1.get_interval()
            breakpoints.append(r)
            self.assertAlmostEqual(l, length)
            self.assertGreaterEqual(l, 0)
            self.assertGreater(r, l)
            self.assertLessEqual(r, ts.get_sequence_length())
            length += r - l
            self.verify_sparse_tree(st1)
            num_trees += 1
        self.assertRaises(StopIteration, next, iter1)
        self.assertRaises(StopIteration, next, iter2)
        self.assertEqual(ts.get_num_trees(), num_trees)
        self.assertEqual(breakpoints, list(ts.breakpoints()))
        self.assertAlmostEqual(length, ts.get_sequence_length())

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
        num_samples = ts.get_sample_size() // 2 + 1
        samples = list(ts.samples())[:num_samples]
        pi1 = ts.get_pairwise_diversity(samples)
        pi2 = simple_get_pairwise_diversity([haplotypes[j] for j in range(num_samples)])
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
            self.assertEqual(st1.parent_dict, st2.parent_dict)
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


class TestSimulator(HighLevelTestCase):
    """
    Runs tests on the underlying Simulator object.
    """

    def verify_dump_load(self, tree_sequence):
        """
        Dump the tree sequence and verify we can load again from the same
        file.
        """
        tree_sequence.dump(self.temp_file)
        other = msprime.load(self.temp_file)
        records = list(tree_sequence.edges())
        other_records = list(other.edges())
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
        self.assertEqual(sim.random_generator, rng)
        sim.run()
        self.assertEqual(sim.num_breakpoints, len(sim.breakpoints))
        self.assertGreater(sim.used_memory, 0)
        self.assertGreater(sim.time, 0)
        self.assertGreater(sim.num_avl_node_blocks, 0)
        self.assertGreater(sim.num_segment_blocks, 0)
        self.assertGreater(sim.num_node_mapping_blocks, 0)
        self.assertGreater(sim.num_node_blocks, 0)
        self.assertGreater(sim.num_edge_blocks, 0)
        self.assertGreater(sim.max_memory, 0)
        tree_sequence = sim.get_tree_sequence()
        t = 0.0
        for record in tree_sequence.nodes():
            if record.time > t:
                t = record.time
        self.assertEqual(sim.time, t)
        self.assertGreater(sim.num_common_ancestor_events, 0)
        self.assertGreaterEqual(sim.num_recombination_events, 0)
        self.assertGreaterEqual(sim.total_num_migration_events, 0)
        self.assertGreaterEqual(sim.num_multiple_recombination_events, 0)
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
        self.assertGreater(sim.avl_node_block_size, 0)
        self.assertGreater(sim.segment_block_size, 0)
        self.assertGreater(sim.node_mapping_block_size, 0)
        self.assertGreater(sim.node_block_size, 0)
        self.assertGreater(sim.edge_block_size, 0)
        sim.reset()
        sim.avl_node_block_size = 1
        sim.segment_block_size = 1
        sim.node_mapping_block_size = 1
        sim.node_block_size = 1
        sim.edge_block_size = 1
        self.assertEqual(sim.avl_node_block_size, 1)
        self.assertEqual(sim.segment_block_size, 1)
        self.assertEqual(sim.node_mapping_block_size, 1)
        self.assertEqual(sim.node_block_size, 1)
        self.assertEqual(sim.edge_block_size, 1)

    def test_bad_inputs(self):
        recomb_map = msprime.RecombinationMap.uniform_map(1, 0)
        for bad_type in ["xd", None, 4.4]:
            self.assertRaises(
                TypeError, msprime.Simulator, [(0, 0), (0, 0)], bad_type)
        self.assertRaises(ValueError, msprime.Simulator, [], recomb_map)
        self.assertRaises(ValueError, msprime.Simulator, [(0, 0)], recomb_map)


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

    @unittest.skip("Recurrent mutations")
    def test_recurrent_mutations_over_samples(self):
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

    @unittest.skip("Recurrent mutations")
    def test_recurrent_mutations_errors(self):
        ts = self.get_tree_sequence()
        tree = next(ts.trees())
        for u in tree.nodes():
            for sample in tree.samples(u):
                if sample != u:
                    site = msprime.Site(
                        index=0, ancestral_state="0", position=0, mutations=[
                            msprime.Mutation(site=0, derived_state="1", node=u),
                            msprime.Mutation(site=0, derived_state="1", node=sample)])
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

    @unittest.skip("Recurrent mutations")
    def test_recurrent_mutations_over_samples(self):
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

    @unittest.skip("Recurrent mutation error")
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
                self.assertRaises(_msprime.LibraryError, list, ts_new.haplotypes())

    def test_back_mutations(self):
        for ts in get_back_mutation_examples():
            self.verify_tree_sequence(ts)


class TestTreeSequence(HighLevelTestCase):
    """
    Tests for the tree sequence object.
    """

    def test_sparse_trees(self):
        for ts in get_example_tree_sequences():
            self.verify_sparse_trees(ts)

    def test_mutations(self):
        # TODO enable the back_mutations here once this has been implemented
        # for pi and variants.
        for ts in get_example_tree_sequences(back_mutations=False):
            self.verify_mutations(ts)

    def verify_edge_diffs(self, ts):
        pts = tests.PythonTreeSequence(ts.get_ll_tree_sequence())
        d1 = list(ts.edge_diffs())
        d2 = list(pts.edge_diffs())
        self.assertEqual(d1, d2)

        # check that we have the correct set of children at all nodes.
        children = collections.defaultdict(set)
        trees = iter(ts.trees())
        tree = next(trees)
        last_right = 0
        for (left, right), edges_out, edges_in in ts.edge_diffs():
            assert left == last_right
            last_right = right
            for edge in edges_out:
                children[edge.parent].remove(edge.child)
            for edge in edges_in:
                children[edge.parent].add(edge.child)
            while tree.interval[1] <= left:
                tree = next(trees)
            # print(left, right, tree.interval)
            self.assertTrue(left >= tree.interval[0])
            self.assertTrue(right <= tree.interval[1])
            for u in tree.nodes():
                if tree.is_internal(u):
                    self.assertIn(u, children)
                    self.assertEqual(children[u], set(tree.children(u)))

    def test_edge_diffs(self):
        for ts in get_example_tree_sequences():
            self.verify_edge_diffs(ts)

    def verify_edgesets(self, ts):
        """
        Verifies that the edgesets we return are equivalent to the original edges.
        """
        new_edges = []
        for edgeset in ts.edgesets():
            self.assertEqual(edgeset.children, sorted(edgeset.children))
            self.assertGreater(len(edgeset.children), 0)
            for child in edgeset.children:
                new_edges.append(msprime.Edge(
                    edgeset.left, edgeset.right, edgeset.parent, child))
        # squash the edges.
        t = ts.dump_tables().nodes.time
        new_edges.sort(key=lambda e: (t[e.parent], e.parent, e.child, e.left))

        squashed = []
        last_e = new_edges[0]
        for e in new_edges[1:]:
            condition = (
                e.parent != last_e.parent or
                e.child != last_e.child or
                e.left != last_e.right)
            if condition:
                squashed.append(last_e)
                last_e = e
            last_e.right = e.right
        squashed.append(last_e)
        edges = list(ts.edges())
        self.assertEqual(len(squashed), len(edges))
        self.assertEqual(edges, squashed)

    def test_edgesets(self):
        for ts in get_example_tree_sequences():
            self.verify_edgesets(ts)

    def verify_tracked_samples(self, ts):
        # Should be empty list by default.
        for tree in ts.trees():
            self.assertEqual(tree.get_num_tracked_samples(), 0)
            for u in tree.nodes():
                self.assertEqual(tree.get_num_tracked_samples(u), 0)
        samples = list(ts.samples())
        tracked_samples = samples[:2]
        for tree in ts.trees(tracked_samples):
            if len(tree.parent_dict) == 0:
                # This is a crude way of checking if we have multiple roots.
                # We'll need to fix this code up properly when we support multiple
                # roots and remove this check
                break
            nu = [0 for j in range(ts.get_num_nodes())]
            self.assertEqual(tree.get_num_tracked_samples(), len(tracked_samples))
            for j in tracked_samples:
                u = j
                while u != msprime.NULL_NODE:
                    nu[u] += 1
                    u = tree.get_parent(u)
            for u, count in enumerate(nu):
                self.assertEqual(tree.get_num_tracked_samples(u), count)

    def test_tracked_samples(self):
        for ts in get_example_tree_sequences():
            self.verify_tracked_samples(ts)

    def test_deprecated_sample_aliases(self):
        for ts in get_example_tree_sequences():
            # Ensure that we get the same results from the various combinations
            # of leaf_lists, sample_lists etc.
            samples = list(ts.samples())[:2]
            # tracked leaves/samples
            trees_new = ts.trees(tracked_samples=samples)
            trees_old = ts.trees(tracked_leaves=samples)
            for t_new, t_old in zip(trees_new, trees_old):
                for u in t_new.nodes():
                    self.assertEqual(
                        t_new.num_tracked_samples(u), t_old.get_num_tracked_leaves(u))
            for on in [True, False]:
                # sample/leaf counts
                trees_new = ts.trees(sample_counts=on)
                trees_old = ts.trees(leaf_counts=on)
                for t_new, t_old in zip(trees_new, trees_old):
                    for u in t_new.nodes():
                        self.assertEqual(t_new.num_samples(u), t_old.get_num_leaves(u))
                        self.assertEqual(
                            list(t_new.samples(u)), list(t_old.get_leaves(u)))
                trees_new = ts.trees(sample_lists=on)
                trees_old = ts.trees(leaf_lists=on)
                for t_new, t_old in zip(trees_new, trees_old):
                    for u in t_new.nodes():
                        self.assertEqual(t_new.num_samples(u), t_old.get_num_leaves(u))
                        self.assertEqual(
                            list(t_new.samples(u)), list(t_old.get_leaves(u)))

    def verify_samples(self, ts):
        # We should get the same list of samples if we use the low-level
        # sample lists or a simple traversal.
        samples1 = []
        for t in ts.trees(sample_lists=False):
            samples1.append(list(t.samples()))
        samples2 = []
        for t in ts.trees(sample_lists=True):
            samples2.append(list(t.samples()))
        self.assertEqual(samples1, samples2)

    def test_samples(self):
        for ts in get_example_tree_sequences():
            self.verify_samples(ts)

    def test_trees_interface(self):
        ts = list(get_example_tree_sequences())[0]
        # The defaults should make sense and count samples.
        # get_num_tracked_samples
        for t in ts.trees():
            self.assertEqual(t.get_num_samples(0), 1)
            self.assertEqual(t.get_num_tracked_samples(0), 0)
            self.assertEqual(list(t.samples(0)), [0])

        for t in ts.trees(sample_counts=False):
            self.assertEqual(t.get_num_samples(0), 1)
            self.assertRaises(RuntimeError, t.get_num_tracked_samples, 0)
            self.assertEqual(list(t.samples(0)), [0])

        for t in ts.trees(sample_counts=True):
            self.assertEqual(t.get_num_samples(0), 1)
            self.assertEqual(t.get_num_tracked_samples(0), 0)
            self.assertEqual(list(t.samples(0)), [0])

        for t in ts.trees(sample_counts=True, tracked_samples=[0]):
            self.assertEqual(t.get_num_samples(0), 1)
            self.assertEqual(t.get_num_tracked_samples(0), 1)
            self.assertEqual(list(t.samples(0)), [0])

        for t in ts.trees(sample_lists=True, sample_counts=True):
            self.assertEqual(t.get_num_samples(0), 1)
            self.assertEqual(t.get_num_tracked_samples(0), 0)
            self.assertEqual(list(t.samples(0)), [0])

        for t in ts.trees(sample_lists=True, sample_counts=False):
            self.assertEqual(t.get_num_samples(0), 1)
            self.assertRaises(RuntimeError, t.get_num_tracked_samples, 0)
            self.assertEqual(list(t.samples(0)), [0])

        # This is a bit weird as we don't seem to actually execute the
        # method until it is iterated.
        self.assertRaises(
            ValueError, list, ts.trees(sample_counts=False, tracked_samples=[0]))

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
            N = ts.get_num_nodes()
            self.assertRaises(ValueError, ts.get_population, -1)
            self.assertRaises(ValueError, ts.get_population, N)
            self.assertRaises(ValueError, ts.get_population, N + 1)
            self.assertEqual(ts.get_population(0), 0)
            self.assertEqual(ts.get_population(N - 1), 0)

    def test_get_time(self):
        for ts in get_example_tree_sequences():
            N = ts.get_num_nodes()
            self.assertRaises(ValueError, ts.get_time, -1)
            self.assertRaises(ValueError, ts.get_time, N)
            self.assertRaises(ValueError, ts.get_time, N + 1)
            for u in range(N):
                self.assertEqual(ts.get_time(u), ts.node(u).time)

    def test_get_samples(self):
        for ts in get_example_tree_sequences():
            samples = []
            for u in range(ts.num_nodes):
                if ts.node(u).is_sample():
                    samples.append(u)
            self.assertEqual(ts.samples(), samples)
            self.assertEqual(ts.samples(0), samples)
            self.assertEqual(ts.samples(msprime.NULL_POPULATION), [])
            self.assertEqual(ts.samples(1), [])
            for t in ts.trees():
                self.assertEqual(sorted(list(t.samples())), samples)

    def test_write_vcf_interface(self):
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            for bad_ploidy in [-1, 0, n + 1]:
                self.assertRaises(ValueError, ts.write_vcf, self.temp_file, bad_ploidy)

    def verify_simplify_topology(self, ts, sample):
        new_ts, node_map = ts.simplify(sample, map_nodes=True)
        if len(sample) == 0:
            self.assertEqual(new_ts.num_nodes, 0)
            self.assertEqual(new_ts.num_edges, 0)
            self.assertEqual(new_ts.num_sites, 0)
            self.assertEqual(new_ts.num_mutations, 0)
        elif len(sample) == 1:
            self.assertEqual(new_ts.num_nodes, 1)
            self.assertEqual(new_ts.num_edges, 0)
        for u in range(ts.num_nodes):
            old_node = ts.node(u)
            if node_map[u] != msprime.NULL_NODE:
                new_node = new_ts.node(node_map[u])
                self.assertEqual(old_node.time, new_node.time)
                self.assertEqual(old_node.population, new_node.population)
                self.assertEqual(old_node.name, new_node.name)
        for u in sample:
            old_node = ts.node(u)
            new_node = new_ts.node(node_map[u])
            self.assertEqual(old_node.flags, new_node.flags)
            self.assertEqual(old_node.time, new_node.time)
            self.assertEqual(old_node.population, new_node.population)
            self.assertEqual(old_node.name, new_node.name)
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
            # If the MRCA of all pairs of samples is the same, then we have the
            # same information. We limit this to at most 500 pairs
            pairs = itertools.islice(itertools.combinations(sample, 2), 500)
            for pair in pairs:
                mapped_pair = [node_map[u] for u in pair]
                mrca1 = old_tree.get_mrca(*pair)
                mrca2 = new_tree.get_mrca(*mapped_pair)
                if mrca1 == msprime.NULL_NODE:
                    self.assertEqual(mrca2, mrca1)
                else:
                    self.assertEqual(mrca2, node_map[mrca1])
                    self.assertEqual(old_tree.get_time(mrca1), new_tree.get_time(mrca2))
                    self.assertEqual(
                        old_tree.get_population(mrca1), new_tree.get_population(mrca2))

    def verify_simplify_haplotypes(self, ts, samples):
        sub_ts, node_map = ts.simplify(
                samples, map_nodes=True, filter_invariant_sites=False)
        # Sites tables should be equal
        self.assertEqual(ts.tables.sites, sub_ts.tables.sites)
        sub_haplotypes = dict(zip(sub_ts.samples(), sub_ts.haplotypes()))
        all_haplotypes = dict(zip(ts.samples(), ts.haplotypes()))
        mapped_ids = []
        for node_id, h in all_haplotypes.items():
            mapped_node_id = node_map[node_id]
            if mapped_node_id in sub_haplotypes:
                self.assertEqual(h, sub_haplotypes[mapped_node_id])
                mapped_ids.append(mapped_node_id)
        self.assertEqual(sorted(mapped_ids), sorted(sub_ts.samples()))

    def verify_simplify_equality(self, ts, sample):
        for filter_invariant_sites in [False, True]:
            s1, node_map1 = ts.simplify(
                sample, map_nodes=True, filter_invariant_sites=filter_invariant_sites)
            t1 = s1.dump_tables()
            s2, node_map2 = simplify_tree_sequence(
                ts, sample, filter_invariant_sites=filter_invariant_sites)
            t2 = s2.dump_tables()
            self.assertTrue(all(node_map1 == node_map2))
            self.assertEqual(t1.nodes, t2.nodes)
            self.assertEqual(t1.edges, t2.edges)
            self.assertEqual(t1.migrations, t2.migrations)
            self.assertEqual(t1.sites, t2.sites)
            self.assertEqual(t1.mutations, t2.mutations)

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
        if len(sample) > 0:
            while j < ts.num_sites:
                unique = np.unique(full_genotypes[j][s])
                self.assertEqual(unique.shape[0], 1)
                self.assertIn(unique[0], [0, 1])
                j += 1

    def test_simplify(self):
        num_mutations = 0
        print("\n\nFIXME\n\nenable internal samples once 1..n property is restored")
        for ts in get_example_tree_sequences(internal_samples=False):
            n = ts.get_sample_size()
            num_mutations += ts.get_num_mutations()
            sample_sizes = {0, 1}
            if n > 2:
                sample_sizes |= set([2, max(2, n // 2), n - 1])
            for k in sample_sizes:
                subset = random.sample(list(ts.samples()), k)
                self.verify_simplify_topology(ts, subset)
                self.verify_simplify_haplotypes(ts, subset)
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
            edges_file = os.path.join(prefix, "{:02d}-edges.txt".format(j))
            sites_file = os.path.join(prefix, "{:02d}-sites.txt".format(j))
            mutations_file = os.path.join(prefix, "{:02d}-mutations.txt".format(j))
            with open(nodes_file) as nodes, \
                    open(edges_file) as edges,\
                    open(sites_file) as sites,\
                    open(mutations_file) as mutations:
                ts = msprime.load_text(
                    nodes=nodes, edges=edges, sites=sites, mutations=mutations)
            samples = list(ts.samples())
            self.verify_simplify_equality(ts, samples)
            j += 1
        self.assertGreater(j, 1)

    def test_apis(self):
        ts = msprime.simulate(10, random_seed=1)
        self.assertEqual(ts.get_ll_tree_sequence(), ts.ll_tree_sequence)
        self.assertEqual(ts.get_provenance(), ts.provenance)
        self.assertEqual(ts.get_sample_size(), ts.sample_size)
        self.assertEqual(ts.get_sample_size(), ts.num_samples)
        self.assertEqual(ts.get_sequence_length(), ts.sequence_length)
        self.assertEqual(ts.num_edges, ts.num_records)
        self.assertEqual(ts.get_num_trees(), ts.num_trees)
        self.assertEqual(ts.get_num_mutations(), ts.num_mutations)
        self.assertEqual(ts.get_num_nodes(), ts.num_nodes)
        self.assertEqual(
            ts.get_pairwise_diversity(), ts.pairwise_diversity())
        samples = ts.samples()
        self.assertEqual(
            ts.get_pairwise_diversity(samples), ts.pairwise_diversity(samples))
        for s in samples:
            self.assertEqual(ts.get_time(s), ts.time(s))
            p = ts.get_population(s)
            self.assertEqual(p, ts.population(s))
            self.assertEqual(ts.get_samples(p), ts.samples(p))
        self.assertEqual(ts.get_samples(), ts.samples())

    @unittest.skip("Copy is broken for back mutations??")
    def test_copy(self):
        for ts1 in get_example_tree_sequences():
            ts2 = ts1.copy()
            self.assertNotEqual(id(ts1), id(ts2))
            self.assertEqual(ts1.sequence_length, ts2.sequence_length)
            self.assertEqual(list(ts1.edges()), list(ts2.edges()))
            self.assertEqual(list(ts1.nodes()), list(ts2.nodes()))
            self.assertEqual(list(ts1.mutations()), list(ts2.mutations()))
            site_lists = [[], list(ts1.sites())[:-1]]
            for sites in site_lists:
                ts2 = ts1.copy(sites=sites)
                self.assertNotEqual(id(ts1), id(ts2))
                self.assertEqual(list(ts1.edges()), list(ts2.edges()))
                self.assertEqual(list(ts1.nodes()), list(ts2.nodes()))
                self.assertEqual(sites, list(ts2.sites()))

    def test_generate_mutations_on_tree_sequence(self):
        some_mutations = False
        for ts in get_example_tree_sequences():
            nodes = msprime.NodeTable()
            edges = msprime.EdgeTable()
            sites = msprime.SiteTable()
            mutations = msprime.MutationTable()
            ts.dump_tables(nodes=nodes, edges=edges)
            mutgen = msprime.MutationGenerator(msprime.RandomGenerator(1), 10)
            mutgen.generate(nodes, edges, sites, mutations)
            if mutations.num_rows > 0:
                some_mutations = True
            tsp = msprime.load_tables(
                nodes=nodes, edges=edges, sites=sites, mutations=mutations)
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

    def test_removed_methods(self):
        ts = next(get_example_tree_sequences())
        self.assertRaises(NotImplementedError, ts.diffs)
        self.assertRaises(NotImplementedError, ts.newick_trees)


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

    def verify_samples_format(self, ts, samples_file, precision):
        """
        Verifies that the samples we output have the correct form:
        same as nodes, but with first column equal to 'id'.
        """
        def convert(v):
            return "{:.{}f}".format(v, precision)
        output_nodes = samples_file.read().splitlines()
        self.assertEqual(len(output_nodes) - 1, len(ts.samples()))
        self.assertEqual(
            list(output_nodes[0].split()),
            ["id", "is_sample", "time", "population"])
        sample_nodes = [ts.node(x) for x in ts.samples()]
        for node_id, node, line in zip(ts.samples(), sample_nodes,
                                       output_nodes[1:]):
            splits = line.split("\t")
            self.assertEqual(str(node_id), splits[0])
            self.assertEqual(str(node.is_sample()), splits[1])
            self.assertEqual(convert(node.time), splits[2])
            self.assertEqual(str(node.population), splits[3])

    def verify_edges_format(self, ts, edges_file, precision):
        """
        Verifies that the edges we output have the correct form.
        """
        def convert(v):
            return "{:.{}f}".format(v, precision)
        output_edges = edges_file.read().splitlines()
        self.assertEqual(len(output_edges) - 1, ts.num_edges)
        self.assertEqual(
            list(output_edges[0].split()),
            ["left", "right", "parent", "child"])
        for edge, line in zip(ts.edges(), output_edges[1:]):
            splits = line.split("\t")
            self.assertEqual(convert(edge.left), splits[0])
            self.assertEqual(convert(edge.right), splits[1])
            self.assertEqual(str(edge.parent), splits[2])
            self.assertEqual(str(edge.child), splits[3])

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
                edges_file = six.StringIO()
                sites_file = six.StringIO()
                mutations_file = six.StringIO()
                ts.dump_text(
                    nodes=nodes_file, edges=edges_file, sites=sites_file,
                    mutations=mutations_file, precision=precision)
                nodes_file.seek(0)
                edges_file.seek(0)
                sites_file.seek(0)
                mutations_file.seek(0)
                self.verify_nodes_format(ts, nodes_file, precision)
                self.verify_edges_format(ts, edges_file, precision)
                self.verify_sites_format(ts, sites_file, precision)
                self.verify_mutations_format(ts, mutations_file, precision)

    def test_dump_samples_text(self):
        for ts in get_example_tree_sequences():
            for precision in [2, 7]:
                samples_file = six.StringIO()
                ts.dump_samples_text(samples_file, precision=precision)
                samples_file.seek(0)
                self.verify_samples_format(ts, samples_file, precision)

    def verify_approximate_equality(self, ts1, ts2):
        """
        Verifies that the specified tree sequences are approximately
        equal, taking into account the error incurred in exporting to text.
        """
        self.assertEqual(ts1.sample_size, ts2.sample_size)
        self.assertAlmostEqual(ts1.sequence_length, ts2.sequence_length)
        self.assertEqual(ts1.num_nodes, ts2.num_nodes)
        self.assertEqual(ts1.num_edges, ts2.num_edges)
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
        for r1, r2 in zip(ts1.edges(), ts2.edges()):
            checked += 1
            self.assertAlmostEqual(r1.left, r2.left)
            self.assertAlmostEqual(r1.right, r2.right)
            self.assertEqual(r1.parent, r2.parent)
            self.assertEqual(r1.child, r2.child)
        self.assertEqual(ts1.num_edges, checked)

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

    @unittest.skip("Add parent column to mutation text format")
    def test_text_record_round_trip(self):
        for ts1 in get_example_tree_sequences():
            nodes_file = six.StringIO()
            edges_file = six.StringIO()
            sites_file = six.StringIO()
            mutations_file = six.StringIO()
            ts1.dump_text(
                nodes=nodes_file, edges=edges_file, sites=sites_file,
                mutations=mutations_file, precision=16)
            nodes_file.seek(0)
            edges_file.seek(0)
            sites_file.seek(0)
            mutations_file.seek(0)
            ts2 = msprime.load_text(
                nodes=nodes_file, edges=edges_file, sites=sites_file,
                mutations=mutations_file, sequence_length=ts1.sequence_length)
            self.verify_approximate_equality(ts1, ts2)

    def test_empty_files(self):
        nodes_file = six.StringIO("is_sample\ttime\n")
        edges_file = six.StringIO("left\tright\tparent\tchild\n")
        sites_file = six.StringIO("position\tancestral_state\n")
        mutations_file = six.StringIO("site\tnode\tderived_state\n")
        self.assertRaises(
            _msprime.LibraryError, msprime.load_text,
            nodes=nodes_file, edges=edges_file, sites=sites_file,
            mutations=mutations_file)

    def test_empty_files_sequence_length(self):
        nodes_file = six.StringIO("is_sample\ttime\n")
        edges_file = six.StringIO("left\tright\tparent\tchild\n")
        sites_file = six.StringIO("position\tancestral_state\n")
        mutations_file = six.StringIO("site\tnode\tderived_state\n")
        ts = msprime.load_text(
                nodes=nodes_file, edges=edges_file, sites=sites_file,
                mutations=mutations_file, sequence_length=100)
        self.assertEqual(ts.sequence_length, 100)
        self.assertEqual(ts.num_nodes, 0)
        self.assertEqual(ts.num_edges, 0)
        self.assertEqual(ts.num_sites, 0)
        self.assertEqual(ts.num_edges, 0)


class TestSparseTree(HighLevelTestCase):
    """
    Some simple tests on the API for the sparse tree.
    """
    def get_tree(self, sample_lists=False):
        ts = msprime.simulate(10, random_seed=1, mutation_rate=1)
        return next(ts.trees(sample_lists=sample_lists))

    def test_str(self):
        t = self.get_tree()
        self.assertIsInstance(str(t), str)
        self.assertEqual(str(t), str(t.get_parent_dict()))

    def test_samples(self):
        for sample_lists in [True, False]:
            t = self.get_tree(sample_lists)
            n = t.get_sample_size()
            all_samples = list(t.samples(t.get_root()))
            self.assertEqual(sorted(all_samples), list(range(n)))
            for j in range(n):
                self.assertEqual(list(t.samples(j)), [j])

            def test_func(t, u):
                """
                Simple test definition of the traversal.
                """
                stack = [u]
                while len(stack) > 0:
                    v = stack.pop()
                    if t.is_sample(v):
                        yield v
                    if t.is_internal(v):
                        for c in reversed(t.get_children(v)):
                            stack.append(c)
            for u in t.nodes():
                l1 = list(t.samples(u))
                l2 = list(test_func(t, u))
                self.assertEqual(l1, l2)
                self.assertEqual(t.get_num_samples(u), len(l1))

    def verify_newick(self, tree):
        """
        Verifies that we output the newick tree as expected.
        """
        if tree.num_roots == 1:
            py_tree = tests.PythonSparseTree.from_sparse_tree(tree)
            newick1 = tree.newick(precision=0, time_scale=0)
            newick2 = py_tree.newick(precision=0, time_scale=0)
            self.assertEqual(newick1, newick2)
        else:
            self.assertRaises(_msprime.LibraryError, tree.newick)

    def test_newick(self):
        for ts in get_example_tree_sequences():
            for tree in ts.trees():
                self.verify_newick(tree)

    def test_traversals(self):
        for ts in get_example_tree_sequences():
            tree = next(ts.trees())
            self.verify_traversals(tree)

    def verify_traversals(self, tree):
        t1 = tree
        t2 = tests.PythonSparseTree.from_sparse_tree(t1)
        self.assertEqual(list(t1.nodes()), list(t2.nodes()))
        orders = ["inorder", "postorder", "levelorder", "breadthfirst"]
        if tree.num_roots == 1:
            self.assertRaises(ValueError, list, t1.nodes(order="bad order"))
            self.assertEqual(list(t1.nodes()), list(t1.nodes(t1.get_root())))
            self.assertEqual(
                list(t1.nodes()),
                list(t1.nodes(t1.get_root(), "preorder")))
            for u in t1.nodes():
                self.assertEqual(list(t1.nodes(u)), list(t2.nodes(u)))
            for test_order in orders:
                self.assertEqual(
                    sorted(list(t1.nodes())),
                    sorted(list(t1.nodes(order=test_order))))
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
        else:
            for test_order in orders:
                all_nodes = []
                for root in t1.roots:
                    self.assertEqual(
                        list(t1.nodes(root, order=test_order)),
                        list(t2.nodes(root, order=test_order)))
                    all_nodes.extend(t1.nodes(root, order=test_order))
                self.assertEqual(all_nodes, list(t1.nodes(order=test_order)))

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
        self.assertEqual(t1.get_total_branch_length(), t1.total_branch_length)
        # node properties
        root = t1.get_root()
        for node in t1.nodes():
            if node != root:
                self.assertEqual(t1.get_time(node), t1.time(node))
                self.assertEqual(t1.get_parent(node), t1.parent(node))
                self.assertEqual(t1.get_children(node), t1.children(node))
                self.assertEqual(t1.get_population(node), t1.population(node))
                self.assertEqual(t1.get_num_samples(node), t1.num_samples(node))
                self.assertEqual(t1.get_branch_length(node),
                                 t1.branch_length(node))
                self.assertEqual(t1.get_num_tracked_samples(node),
                                 t1.num_tracked_samples(node))

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
        rng = sim.random_generator
        self.assertIsInstance(rng, msprime.RandomGenerator)
        self.assertNotEqual(rng.get_seed(), 0)

    def test_random_seed(self):
        seed = 12345
        rng = msprime.RandomGenerator(seed)
        sim = msprime.simulator_factory(10, random_generator=rng)
        self.assertEqual(rng, sim.random_generator)
        self.assertEqual(rng.get_seed(), seed)

    def test_sample_size(self):
        self.assertRaises(ValueError, msprime.simulator_factory)
        self.assertRaises(ValueError, msprime.simulator_factory, 1)
        self.assertRaises(
            ValueError, msprime.simulator_factory, sample_size=1)
        for n in [2, 100, 1000]:
            sim = msprime.simulator_factory(n)
            self.assertEqual(sim.sample_size, n)
            ll_sim = sim.create_ll_instance()
            self.assertEqual(ll_sim.get_num_samples(), n)
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
            self.assertEqual(sim.model.population_size, Ne)
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
            pop_configs = [msprime.PopulationConfiguration(5) for _ in range(N)]
            sample_size = 5 * N
            sim = msprime.simulator_factory(population_configurations=pop_configs)
            self.assertEqual(sim.population_configurations, pop_configs)
            self.assertEqual(sim.sample_size, sample_size)
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
                ValueError, msprime.simulator_factory, population_configurations=configs)
            configs = [msprime.PopulationConfiguration(2) for _ in range(d)]
            sim = msprime.simulator_factory(population_configurations=configs)
            self.assertEqual(sim.sample_size, 2 * d)
            samples = []
            for j in range(d):
                samples += [msprime.Sample(population=j, time=0) for _ in range(2)]
            self.assertEqual(sim.samples, samples)
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
            self.assertEqual(sim.migration_matrix, hl_matrix)
            ll_sim = sim.create_ll_instance()
            ll_matrix = [v for row in hl_matrix for v in row]
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
            recomb_map = sim.recombination_map
            self.assertEqual(recomb_map.get_positions(), [0, 1], [rate, 0])
            self.assertEqual(
                recomb_map.get_num_loci(),
                msprime.RecombinationMap.DEFAULT_NUM_LOCI)

    def test_recombination_rate_scaling(self):
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
            total_rate = length * rate
            per_locus_rate = total_rate / (num_loci - 1)
            # We expect all these rates to be positive.
            self.assertGreater(per_locus_rate, 0)
            ll_sim = sim.create_ll_instance()
            self.assertAlmostEqual(per_locus_rate, ll_sim.get_recombination_rate())

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
            self.assertEqual(sim.recombination_map, recomb_map)
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
        self.assertEqual(sim.samples, samples)
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
        for r1, r2 in zip(ts1.edges(), ts2.edges()):
            self.assertEqual(r1.parent, r2.parent)
            self.assertEqual(r1.child, r2.child)
            if approx:
                self.assertAlmostEqual(r1.left, r2.left)
                self.assertAlmostEqual(r1.right, r2.right)
            else:
                self.assertEqual(r1.left, r2.left)
                self.assertEqual(r1.right, r2.right)
            j += 1
        self.assertEqual(ts1.num_edges, j)
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
        edge_table = msprime.EdgeTable()
        for e in ts.edges():
            edge_table.add_row(
                left=e.left, right=e.right, parent=node_map[e.parent],
                child=node_map[e.child])
        msprime.sort_tables(nodes=node_table, edges=edge_table)
        other_ts = msprime.load_tables(nodes=node_table, edges=edge_table)

        self.assertEqual(ts.get_num_trees(), other_ts.get_num_trees())
        self.assertEqual(ts.get_sample_size(), other_ts.get_sample_size())
        self.assertEqual(ts.get_num_nodes(), other_ts.get_num_nodes())
        j = 0
        for t1, t2 in zip(ts.trees(), other_ts.trees()):
            # Verify the topologies are identical. We do this by traversing
            # upwards to the root for every sample and checking if we map to
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
        edges_file = six.StringIO()
        # Also verify we can read the text version.
        other_ts.dump_text(nodes=nodes_file, edges=edges_file, precision=14)
        nodes_file.seek(0)
        edges_file.seek(0)
        ts3 = msprime.load_text(nodes_file, edges_file)
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
        for t in ts.trees():
            for u in t.nodes():
                if len(t.children(u)) > 2:
                    found = True
                    break
            if found:
                break
        self.assertTrue(found)
        for _ in range(self.num_random_permutations):
            self.verify_random_permutation(ts)
