#
# Copyright (C) 2015-2018 University of Oxford
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
import datetime
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
import warnings
import multiprocessing
import uuid as _uuid

import numpy as np

import msprime
import _msprime
import tests
import tests.tsutil as tsutil
import tests.simplify as simplify


def insert_uniform_mutations(tables, num_mutations, nodes):
    """
    Returns n evenly mutations over the specified list of nodes.
    """
    for j in range(num_mutations):
        tables.sites.add_row(
            position=j * (tables.sequence_length / num_mutations), ancestral_state='0',
            metadata=json.dumps({"index": j}).encode())
        tables.mutations.add_row(
            site=j, derived_state='1', node=nodes[j % len(nodes)],
            metadata=json.dumps({"index": j}).encode())


def get_table_collection_copy(tables, sequence_length):
    """
    Returns a copy of the specified table collection with the specified
    sequence length.
    """
    ll_tables = _msprime.TableCollection(
        individuals=tables.individuals.ll_table,
        nodes=tables.nodes.ll_table,
        edges=tables.edges.ll_table,
        migrations=tables.migrations.ll_table,
        sites=tables.sites.ll_table,
        mutations=tables.mutations.ll_table,
        populations=tables.populations.ll_table,
        provenances=tables.provenances.ll_table,
        sequence_length=sequence_length)
    assert ll_tables.sequence_length == sequence_length
    tables = msprime.TableCollection(ll_tables=ll_tables)
    assert tables.sequence_length == sequence_length
    return tables


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
    L = ts.sequence_length + length
    tables = get_table_collection_copy(tables, L)
    tables.edges.clear()
    tables.sites.clear()
    tables.mutations.clear()
    for left, right, parent, child in new_edges:
        tables.edges.add_row(left, right, parent, child)
    tables.sort()
    # Throw in a bunch of mutations over the whole sequence on the samples.
    insert_uniform_mutations(tables, 100, list(ts.samples()))
    return tables.tree_sequence()


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
    tables = get_table_collection_copy(ts.dump_tables(), 2)
    tables.sites.clear()
    tables.mutations.clear()
    insert_uniform_mutations(tables, 100, list(ts.samples()))
    yield tables.tree_sequence()


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
    yield tables.tree_sequence()

    # Set just internal nodes to be samples.
    flags[:] = 0
    flags[n:] = msprime.NODE_IS_SAMPLE
    nodes.set_columns(flags=flags, time=nodes.time, population=nodes.population)
    yield tables.tree_sequence()

    # Set a mixture of internal and leaf samples.
    flags[:] = 0
    flags[n // 2: n + n // 2] = msprime.NODE_IS_SAMPLE
    nodes.set_columns(flags=flags, time=nodes.time, population=nodes.population)
    yield tables.tree_sequence()


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
    seed = 1
    for n in [2, 3, 10, 100]:
        for m in [1, 2, 32]:
            for rho in [0, 0.1, 0.5]:
                recomb_map = msprime.RecombinationMap.uniform_map(m, rho, num_loci=m)
                ts = msprime.simulate(
                    recombination_map=recomb_map, mutation_rate=0.1,
                    random_seed=seed,
                    population_configurations=[
                        msprime.PopulationConfiguration(n),
                        msprime.PopulationConfiguration(0)],
                    migration_matrix=[[0, 1], [1, 0]])
                ts = tsutil.insert_random_ploidy_individuals(ts, 4, seed=seed)
                yield tsutil.add_random_metadata(ts, seed=seed)
                seed += 1
    for ts in get_bottleneck_examples():
        yield msprime.mutate(
            ts, rate=0.1, random_seed=seed,
            model=msprime.InfiniteSites(msprime.NUCLEOTIDES))
    ts = msprime.simulate(15, length=4, recombination_rate=1)
    assert ts.num_trees > 1
    if back_mutations:
        yield tsutil.insert_branch_mutations(ts, mutations_per_branch=2)
    ts = tsutil.insert_multichar_mutations(ts)
    yield ts
    yield tsutil.add_random_metadata(ts)


def get_bottleneck_examples():
    """
    Returns an iterator of example tree sequences with nonbinary
    trees.
    """
    bottlenecks = [
        msprime.SimpleBottleneck(0.01, 0, proportion=0.05),
        msprime.SimpleBottleneck(0.02, 0, proportion=0.25),
        msprime.SimpleBottleneck(0.03, 0, proportion=1)]
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
    for j in [1, 2, 3]:
        yield tsutil.insert_branch_mutations(ts, mutations_per_branch=j)
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


def simplify_tree_sequence(ts, samples, filter_sites=True):
    """
    Simple tree-by-tree algorithm to get a simplify of a tree sequence.
    """
    s = simplify.Simplifier(
        ts, samples, filter_sites=filter_sites)
    return s.simplify()


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
                length = st.get_time(st.get_parent(u)) - st.get_time(u)
                self.assertGreater(length, 0.0)
                self.assertEqual(st.get_branch_length(u), length)
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
        self.assertIsNotNone(other.file_uuid)
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
        self.assertGreater(sim.time, 0)
        self.assertGreater(sim.num_avl_node_blocks, 0)
        self.assertGreater(sim.num_segment_blocks, 0)
        self.assertGreater(sim.num_node_mapping_blocks, 0)
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
        sim.reset()
        sim.avl_node_block_size = 1
        sim.segment_block_size = 1
        sim.node_mapping_block_size = 1
        self.assertEqual(sim.avl_node_block_size, 1)
        self.assertEqual(sim.segment_block_size, 1)
        self.assertEqual(sim.node_mapping_block_size, 1)

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
            10, length=10, recombination_rate=1, mutation_rate=10, random_seed=3)
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

    def test_as_bytes_fails(self):
        ts = tsutil.insert_multichar_mutations(self.get_tree_sequence())
        self.assertRaises(ValueError, list, ts.variants(as_bytes=True))

    def test_multichar_alleles(self):
        ts = tsutil.insert_multichar_mutations(self.get_tree_sequence())
        for var in ts.variants():
            self.assertEqual(len(var.alleles), 2)
            self.assertEqual(var.site.ancestral_state, var.alleles[0])
            self.assertEqual(var.site.mutations[0].derived_state, var.alleles[1])
            self.assertTrue(all(0 <= var.genotypes))
            self.assertTrue(all(var.genotypes <= 1))

    def test_many_alleles(self):
        ts = self.get_tree_sequence()
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        # This gives us a total of 360 permutations.
        alleles = list(map("".join, itertools.permutations('ABCDEF', 4)))
        self.assertGreater(len(alleles), 255)
        tables.sites.add_row(0, alleles[0])
        parent = -1
        num_alleles = 1
        for allele in alleles[1:]:
            ts = tables.tree_sequence()
            if num_alleles > 255:
                self.assertRaises(_msprime.LibraryError, next, ts.variants())
            else:
                var = next(ts.variants())
                self.assertEqual(len(var.alleles), num_alleles)
                self.assertEqual(list(var.alleles), alleles[:num_alleles])
                self.assertEqual(
                    var.alleles[var.genotypes[0]], alleles[num_alleles - 1])
                for u in ts.samples():
                    if u != 0:
                        self.assertEqual(var.alleles[var.genotypes[u]], alleles[0])
            tables.mutations.add_row(0, 0, allele, parent=parent)
            parent += 1
            num_alleles += 1

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

    def test_genotype_matrix(self):
        ts = self.get_tree_sequence()
        G = np.empty((ts.num_sites, ts.num_samples), dtype=np.uint8)
        for v in ts.variants():
            G[v.index, :] = v.genotypes
        self.assertTrue(np.array_equal(G, ts.genotype_matrix()))

    def test_recurrent_mutations_over_samples(self):
        ts = self.get_tree_sequence()
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        num_sites = 5
        for j in range(num_sites):
            tables.sites.add_row(
                position=j * ts.sequence_length / num_sites,
                ancestral_state="0")
            for u in range(ts.sample_size):
                tables.mutations.add_row(site=j, node=u, derived_state="1")
        ts = tables.tree_sequence()
        variants = list(ts.variants())
        self.assertEqual(len(variants), num_sites)
        for site, variant in zip(ts.sites(), variants):
            self.assertEqual(site.position, variant.position)
            self.assertEqual(site, variant.site)
            self.assertEqual(site.id, variant.index)
            self.assertEqual(variant.alleles, ("0", "1"))
            self.assertTrue(np.all(variant.genotypes == np.ones(ts.sample_size)))

    def test_recurrent_mutations_errors(self):
        ts = self.get_tree_sequence()
        tree = next(ts.trees())
        tables = ts.dump_tables()
        for u in tree.nodes():
            for sample in tree.samples(u):
                if sample != u:
                    tables.sites.clear()
                    tables.mutations.clear()
                    site = tables.sites.add_row(position=0, ancestral_state="0")
                    tables.mutations.add_row(site=site, node=u, derived_state="1")
                    tables.mutations.add_row(site=site, node=sample, derived_state="1")
                    ts_new = tables.tree_sequence()
                    self.assertRaises(_msprime.LibraryError, list, ts_new.variants())

    def test_zero_samples(self):
        ts = self.get_tree_sequence()
        for var1, var2 in zip(ts.variants(), ts.variants(samples=[])):
            self.assertEqual(var1.site, var2.site)
            self.assertEqual(var1.alleles, var2.alleles)
            self.assertEqual(var2.genotypes.shape[0], 0)

    def test_samples(self):
        n = 4
        ts = msprime.simulate(
            n, length=5, recombination_rate=1, mutation_rate=5, random_seed=2)
        self.assertGreater(ts.num_sites, 1)
        samples = list(range(n))
        # Generate all possible sample lists.
        for j in range(n + 1):
            for s in itertools.permutations(samples, j):
                s = np.array(s, dtype=np.int32)
                count = 0
                for var1, var2 in zip(ts.variants(), ts.variants(samples=s)):
                    self.assertEqual(var1.site, var2.site)
                    self.assertEqual(var1.alleles, var2.alleles)
                    self.assertEqual(var2.genotypes.shape, (len(s),))
                    self.assertTrue(np.array_equal(var1.genotypes[s], var2.genotypes))
                    count += 1
                self.assertEqual(count, ts.num_sites)

    def test_non_sample_samples(self):
        # We don't have to use sample nodes. This does make the terminology confusing
        # but it's probably still the best option.
        ts = msprime.simulate(
            10, length=5, recombination_rate=1, mutation_rate=5, random_seed=2)
        tables = ts.dump_tables()
        tables.nodes.set_columns(
            flags=np.zeros_like(tables.nodes.flags) + msprime.NODE_IS_SAMPLE,
            time=tables.nodes.time)
        all_samples_ts = tables.tree_sequence()
        self.assertEqual(all_samples_ts.num_samples, ts.num_nodes)

        count = 0
        samples = range(ts.num_nodes)
        for var1, var2 in zip(all_samples_ts.variants(), ts.variants(samples=samples)):
            self.assertEqual(var1.site, var2.site)
            self.assertEqual(var1.alleles, var2.alleles)
            self.assertEqual(var2.genotypes.shape, (len(samples),))
            self.assertTrue(np.array_equal(var1.genotypes, var2.genotypes))
            count += 1
        self.assertEqual(count, ts.num_sites)


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

    def test_acgt_mutations(self):
        ts = msprime.simulate(10, mutation_rate=10)
        self.assertGreater(ts.num_sites, 0)
        tables = ts.tables
        sites = tables.sites
        mutations = tables.mutations
        sites.set_columns(
            position=sites.position,
            ancestral_state=np.zeros(ts.num_sites, dtype=np.int8) + ord("A"),
            ancestral_state_offset=np.arange(ts.num_sites + 1, dtype=np.uint32))
        mutations.set_columns(
            site=mutations.site,
            node=mutations.node,
            derived_state=np.zeros(ts.num_sites, dtype=np.int8) + ord("T"),
            derived_state_offset=np.arange(ts.num_sites + 1, dtype=np.uint32))
        tsp = msprime.load_tables(**tables.asdict())
        H = [h.replace("0", "A").replace("1", "T") for h in ts.haplotypes()]
        self.assertEqual(H, list(tsp.haplotypes()))

    def test_multiletter_mutations(self):
        ts = msprime.simulate(10)
        tables = ts.tables
        sites = tables.sites
        sites.add_row(0, "ACTG")
        tsp = msprime.load_tables(**tables.asdict())
        self.assertRaises(_msprime.LibraryError, list, tsp.haplotypes())

    def test_recurrent_mutations_over_samples(self):
        for ts in get_bottleneck_examples():
            num_sites = 5
            tables = ts.dump_tables()
            for j in range(num_sites):
                tables.sites.add_row(
                    position=j * ts.sequence_length / num_sites,
                    ancestral_state="0")
                for u in range(ts.sample_size):
                    tables.mutations.add_row(site=j, node=u, derived_state="1")
            ts_new = tables.tree_sequence()
            ones = "1" * num_sites
            for h in ts_new.haplotypes():
                self.assertEqual(ones, h)

    def test_recurrent_mutations_errors(self):
        for ts in get_bottleneck_examples():
            tables = ts.dump_tables()
            tree = next(ts.trees())
            for u in tree.children(tree.root):
                tables.sites.clear()
                tables.mutations.clear()
                site = tables.sites.add_row(position=0, ancestral_state="0")
                tables.mutations.add_row(site=site, node=u, derived_state="1")
                tables.mutations.add_row(site=site, node=tree.root, derived_state="1")
                ts_new = tables.tree_sequence()
                self.assertRaises(_msprime.LibraryError, list, ts_new.haplotypes())
                ts_new.haplotypes()

    def test_back_mutations(self):
        for ts in get_back_mutation_examples():
            self.verify_tree_sequence(ts)


class TestNumpySamples(unittest.TestCase):
    """
    Tests that we correctly handle samples as numpy arrays when passed to
    various methods.
    """
    def get_tree_sequence(self, num_demes=4):
        n = 40
        return msprime.simulate(
            samples=[
                msprime.Sample(time=0, population=j % num_demes) for j in range(n)],
            population_configurations=[
                msprime.PopulationConfiguration() for _ in range(num_demes)],
            migration_matrix=[
                [int(j != k) for j in range(num_demes)] for k in range(num_demes)],
            random_seed=1,
            mutation_rate=10)

    def test_samples(self):
        d = 4
        ts = self.get_tree_sequence(d)
        self.assertTrue(np.array_equal(
            ts.samples(), np.arange(ts.num_samples, dtype=np.int32)))
        total = 0
        for pop in range(d):
            subsample = ts.samples(pop)
            total += subsample.shape[0]
            self.assertTrue(np.array_equal(subsample, ts.samples(population=pop)))
            self.assertEqual(
                list(subsample),
                [node.id for node in ts.nodes()
                    if node.population == pop and node.is_sample()])
        self.assertEqual(total, ts.num_samples)

    def test_genotype_matrix_indexing(self):
        num_demes = 4
        ts = self.get_tree_sequence(num_demes)
        G = ts.genotype_matrix()
        for d in range(num_demes):
            samples = ts.samples(population=d)
            total = 0
            for tree in ts.trees(tracked_samples=samples):
                for mutation in tree.mutations():
                    total += tree.num_tracked_samples(mutation.node)
            self.assertEqual(total, np.sum(G[:, samples]))

    def test_genotype_indexing(self):
        num_demes = 6
        ts = self.get_tree_sequence(num_demes)
        for d in range(num_demes):
            samples = ts.samples(population=d)
            total = 0
            for tree in ts.trees(tracked_samples=samples):
                for mutation in tree.mutations():
                    total += tree.num_tracked_samples(mutation.node)
            other_total = 0
            for variant in ts.variants():
                other_total += np.sum(variant.genotypes[samples])
            self.assertEqual(total, other_total)

    def test_pairwise_diversity(self):
        num_demes = 6
        ts = self.get_tree_sequence(num_demes)
        pi1 = ts.pairwise_diversity(ts.samples())
        pi2 = ts.pairwise_diversity()
        self.assertEqual(pi1, pi2)
        for d in range(num_demes):
            samples = ts.samples(population=d)
            pi1 = ts.pairwise_diversity(samples)
            pi2 = ts.pairwise_diversity(list(samples))
            self.assertEqual(pi1, pi2)

    def test_simplify(self):
        num_demes = 3
        ts = self.get_tree_sequence(num_demes)
        sts = ts.simplify(samples=ts.samples())
        self.assertEqual(ts.num_samples, sts.num_samples)
        for d in range(num_demes):
            samples = ts.samples(population=d)
            sts = ts.simplify(samples=samples)
            self.assertEqual(sts.num_samples, samples.shape[0])


class TestTreeSequence(HighLevelTestCase):
    """
    Tests for the tree sequence object.
    """

    def test_sparse_trees(self):
        for ts in get_example_tree_sequences():
            self.verify_sparse_trees(ts)

    def test_mutations(self):
        for ts in get_example_tree_sequences():
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

    def verify_coalescence_records(self, ts):
        """
        Checks that the coalescence records we output are correct.
        """
        edgesets = list(ts.edgesets())
        records = list(ts.records())
        self.assertEqual(len(edgesets), len(records))
        for edgeset, record in zip(edgesets, records):
            self.assertEqual(edgeset.left, record.left)
            self.assertEqual(edgeset.right, record.right)
            self.assertEqual(edgeset.parent, record.node)
            self.assertEqual(edgeset.children, record.children)
            parent = ts.node(edgeset.parent)
            self.assertEqual(parent.time, record.time)
            self.assertEqual(parent.population, record.population)

    def test_coalescence_records(self):
        for ts in get_example_tree_sequences():
            self.verify_coalescence_records(ts)

    def test_compute_mutation_parent(self):
        for ts in get_example_tree_sequences():
            tables = ts.dump_tables()
            before = tables.mutations.parent[:]
            tables.compute_mutation_parents()
            parent = ts.tables.mutations.parent
            self.assertTrue(np.array_equal(parent, before))

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
            pops = set(node.population for node in ts.nodes())
            for pop in pops:
                subsample = ts.samples(pop)
                self.assertTrue(np.array_equal(subsample, ts.samples(population=pop)))
                self.assertTrue(np.array_equal(subsample, ts.samples(population_id=pop)))
                self.assertEqual(
                    list(subsample),
                    [node.id for node in ts.nodes()
                        if node.population == pop and node.is_sample()])
            self.assertRaises(ValueError, ts.samples, population=0, population_id=0)

    def test_first(self):
        for ts in get_example_tree_sequences():
            t1 = ts.first()
            t2 = next(ts.trees())
            self.assertFalse(t1 is t2)
            self.assertEqual(t1.parent_dict, t2.parent_dict)

    def test_trees_interface(self):
        ts = list(get_example_tree_sequences())[0]
        # The defaults should make sense and count samples.
        # get_num_tracked_samples
        for t in ts.trees():
            self.assertEqual(t.get_num_samples(0), 1)
            self.assertEqual(t.get_num_tracked_samples(0), 0)
            self.assertEqual(list(t.samples(0)), [0])
            self.assertIs(t.tree_sequence, ts)

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

    def test_get_pairwise_diversity(self):
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            self.assertRaises(ValueError, ts.get_pairwise_diversity, [])
            self.assertRaises(ValueError, ts.get_pairwise_diversity, [1])
            self.assertRaises(ValueError, ts.get_pairwise_diversity, [1, n])
            samples = list(ts.samples())
            if any(len(site.mutations) > 1 for site in ts.sites()):
                # Multi-mutations are not currenty supported when computing pi.
                self.assertRaises(
                    _msprime.LibraryError, ts.get_pairwise_diversity)
            else:
                self.assertEqual(
                    ts.get_pairwise_diversity(),
                    ts.get_pairwise_diversity(samples))
                self.assertEqual(
                    ts.get_pairwise_diversity(samples[:2]),
                    ts.get_pairwise_diversity(reversed(samples[:2])))

    def test_populations(self):
        more_than_zero = False
        for ts in get_example_tree_sequences():
            N = ts.num_populations
            if N > 0:
                more_than_zero = True
            pops = list(ts.populations())
            self.assertEqual(len(pops), N)
            for j in range(N):
                self.assertEqual(pops[j], ts.population(j))
                self.assertEqual(pops[j].id, j)
                self.assertTrue(isinstance(pops[j].metadata, bytes))
        self.assertTrue(more_than_zero)

    def test_individuals(self):
        more_than_zero = False
        mapped_to_nodes = False
        for ts in get_example_tree_sequences():
            ind_node_map = collections.defaultdict(list)
            for node in ts.nodes():
                if node.individual != msprime.NULL_INDIVIDUAL:
                    ind_node_map[node.individual].append(node.id)
            if len(ind_node_map) > 0:
                mapped_to_nodes = True
            N = ts.num_individuals
            if N > 0:
                more_than_zero = True
            inds = list(ts.individuals())
            self.assertEqual(len(inds), N)
            for j in range(N):
                self.assertEqual(inds[j], ts.individual(j))
                self.assertEqual(inds[j].id, j)
                self.assertTrue(isinstance(inds[j].metadata, bytes))
                self.assertTrue(isinstance(inds[j].location, np.ndarray))
                self.assertTrue(isinstance(inds[j].nodes, np.ndarray))
                self.assertEqual(ind_node_map[j], list(inds[j].nodes))

        self.assertTrue(more_than_zero)
        self.assertTrue(mapped_to_nodes)

    def test_get_population(self):
        # Deprecated interface for ts.node(id).population
        for ts in get_example_tree_sequences():
            N = ts.get_num_nodes()
            self.assertRaises(ValueError, ts.get_population, -1)
            self.assertRaises(ValueError, ts.get_population, N)
            self.assertRaises(ValueError, ts.get_population, N + 1)
            for node in [0, N - 1]:
                self.assertEqual(ts.get_population(node), ts.node(node).population)

    def test_get_time(self):
        # Deprecated interface for ts.node(id).time
        for ts in get_example_tree_sequences():
            N = ts.get_num_nodes()
            self.assertRaises(ValueError, ts.get_time, -1)
            self.assertRaises(ValueError, ts.get_time, N)
            self.assertRaises(ValueError, ts.get_time, N + 1)
            for u in range(N):
                self.assertEqual(ts.get_time(u), ts.node(u).time)

    def test_write_vcf_interface(self):
        for ts in get_example_tree_sequences():
            n = ts.get_sample_size()
            for bad_ploidy in [-1, 0, n + 1]:
                self.assertRaises(ValueError, ts.write_vcf, self.temp_file, bad_ploidy)

    def verify_simplify_provenance(self, ts):
        new_ts = ts.simplify()
        self.assertEqual(new_ts.num_provenances, ts.num_provenances + 1)
        old = list(ts.provenances())
        new = list(new_ts.provenances())
        self.assertEqual(old, new[:-1])
        # TODO call verify_provenance on this.
        self.assertGreater(len(new[-1].timestamp), 0)
        self.assertGreater(len(new[-1].record), 0)

        new_ts = ts.simplify(record_provenance=False)
        self.assertEqual(new_ts.tables.provenances, ts.tables.provenances)

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
        # The output samples should be 0...n
        self.assertEqual(new_ts.num_samples, len(sample))
        self.assertEqual(list(range(len(sample))), list(new_ts.samples()))
        for j in range(new_ts.num_samples):
            self.assertEqual(node_map[sample[j]], j)
        for u in range(ts.num_nodes):
            old_node = ts.node(u)
            if node_map[u] != msprime.NULL_NODE:
                new_node = new_ts.node(node_map[u])
                self.assertEqual(old_node.time, new_node.time)
                self.assertEqual(old_node.population, new_node.population)
                self.assertEqual(old_node.metadata, new_node.metadata)
        for u in sample:
            old_node = ts.node(u)
            new_node = new_ts.node(node_map[u])
            self.assertEqual(old_node.flags, new_node.flags)
            self.assertEqual(old_node.time, new_node.time)
            self.assertEqual(old_node.population, new_node.population)
            self.assertEqual(old_node.metadata, new_node.metadata)
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

    def verify_simplify_equality(self, ts, sample):
        for filter_sites in [False, True]:
            s1, node_map1 = ts.simplify(
                sample, map_nodes=True, filter_sites=filter_sites)
            t1 = s1.dump_tables()
            s2, node_map2 = simplify_tree_sequence(ts, sample, filter_sites=filter_sites)
            t2 = s2.dump_tables()
            self.assertEqual(s1.num_samples,  len(sample))
            self.assertEqual(s2.num_samples,  len(sample))
            self.assertTrue(all(node_map1 == node_map2))
            self.assertEqual(t1.individuals, t2.individuals)
            self.assertEqual(t1.nodes, t2.nodes)
            self.assertEqual(t1.edges, t2.edges)
            self.assertEqual(t1.migrations, t2.migrations)
            self.assertEqual(t1.sites, t2.sites)
            self.assertEqual(t1.mutations, t2.mutations)
            self.assertEqual(t1.populations, t2.populations)

    def verify_simplify_variants(self, ts, sample):
        subset = ts.simplify(sample)
        sample_map = {u: j for j, u in enumerate(ts.samples())}
        # Need to map IDs back to their sample indexes
        s = np.array([sample_map[u] for u in sample])
        # Build a map of genotypes by position
        full_genotypes = {}
        for variant in ts.variants():
            alleles = [variant.alleles[g] for g in variant.genotypes]
            full_genotypes[variant.position] = alleles
        for variant in subset.variants():
            if variant.position in full_genotypes:
                a1 = [full_genotypes[variant.position][u] for u in s]
                a2 = [variant.alleles[g] for g in variant.genotypes]
                self.assertEqual(a1, a2)

    def test_simplify(self):
        num_mutations = 0
        for ts in get_example_tree_sequences():
            self.verify_simplify_provenance(ts)
            n = ts.get_sample_size()
            num_mutations += ts.get_num_mutations()
            sample_sizes = {0, 1}
            if n > 2:
                sample_sizes |= set([2, max(2, n // 2), n - 1])
            for k in sample_sizes:
                subset = random.sample(list(ts.samples()), k)
                self.verify_simplify_topology(ts, subset)
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
                    nodes=nodes, edges=edges, sites=sites, mutations=mutations,
                    strict=False)
            samples = list(ts.samples())
            self.verify_simplify_equality(ts, samples)
            j += 1
        self.assertGreater(j, 1)

    def test_simplify_migrations_fails(self):
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(10),
                msprime.PopulationConfiguration(10)],
            migration_matrix=[[0, 1], [1, 0]],
            random_seed=2,
            record_migrations=True)
        self.assertGreater(ts.num_migrations, 0)
        # We don't support simplify with migrations, so should fail.
        with self.assertRaises(_msprime.LibraryError):
            ts.simplify()

    def test_deprecated_apis(self):
        ts = msprime.simulate(10, random_seed=1)
        self.assertEqual(ts.get_ll_tree_sequence(), ts.ll_tree_sequence)
        self.assertEqual(ts.get_sample_size(), ts.sample_size)
        self.assertEqual(ts.get_sample_size(), ts.num_samples)
        self.assertEqual(ts.get_sequence_length(), ts.sequence_length)
        self.assertEqual(ts.get_num_trees(), ts.num_trees)
        self.assertEqual(ts.get_num_mutations(), ts.num_mutations)
        self.assertEqual(ts.get_num_nodes(), ts.num_nodes)
        self.assertEqual(ts.get_pairwise_diversity(), ts.pairwise_diversity())
        samples = ts.samples()
        self.assertEqual(
            ts.get_pairwise_diversity(samples), ts.pairwise_diversity(samples))
        self.assertTrue(np.array_equal(ts.get_samples(), ts.samples()))

    def test_generate_mutations_on_tree_sequence(self):
        some_mutations = False
        for ts in get_example_tree_sequences():
            tables = ts.dump_tables()
            mutgen = msprime.MutationGenerator(msprime.RandomGenerator(1), 10)
            tables.sites.clear()
            tables.mutations.clear()
            mutgen.generate(tables.ll_tables)
            if tables.mutations.num_rows > 0:
                some_mutations = True
            tsp = msprime.load_tables(**tables.asdict())
            self.assertEqual(tsp.num_mutations, tables.mutations.num_rows)
        self.assertTrue(some_mutations)

    def test_sites(self):
        some_sites = False
        for ts in get_example_tree_sequences():
            tables = ts.dump_tables()
            sites = tables.sites
            mutations = tables.mutations
            self.assertEqual(ts.num_sites, len(sites))
            self.assertEqual(ts.num_mutations, len(mutations))
            previous_pos = -1
            mutation_index = 0
            ancestral_state = msprime.unpack_strings(
                sites.ancestral_state, sites.ancestral_state_offset)
            derived_state = msprime.unpack_strings(
                mutations.derived_state, mutations.derived_state_offset)

            for index, site in enumerate(ts.sites()):
                s2 = ts.site(site.id)
                self.assertEqual(s2, site)
                self.assertEqual(site.position, sites.position[index])
                self.assertGreater(site.position, previous_pos)
                previous_pos = site.position
                self.assertEqual(ancestral_state[index], site.ancestral_state)
                self.assertEqual(site.id, index)
                for mutation in site.mutations:
                    m2 = ts.mutation(mutation.id)
                    self.assertEqual(m2, mutation)
                    self.assertEqual(mutation.site, site.id)
                    self.assertEqual(mutation.site, mutations.site[mutation_index])
                    self.assertEqual(mutation.node, mutations.node[mutation_index])
                    self.assertEqual(mutation.parent, mutations.parent[mutation_index])
                    self.assertEqual(mutation.id, mutation_index)
                    self.assertEqual(
                        derived_state[mutation_index], mutation.derived_state)
                    mutation_index += 1
                some_sites = True
            total_sites = 0
            for tree in ts.trees():
                self.assertEqual(len(list(tree.sites())), tree.num_sites)
                total_sites += tree.num_sites
            self.assertEqual(ts.num_sites, total_sites)
            self.assertEqual(mutation_index, len(mutations))
        self.assertTrue(some_sites)

    def verify_mutations(self, ts):
        other_mutations = []
        for site in ts.sites():
            for mutation in site.mutations:
                other_mutations.append(mutation)
        mutations = list(ts.mutations())
        self.assertEqual(ts.num_mutations, len(other_mutations))
        self.assertEqual(ts.num_mutations, len(mutations))
        for mut, other_mut in zip(mutations, other_mutations):
            # We cannot compare these directly as the mutations obtained
            # from the mutations iterator will have extra deprecated
            # attributes.
            self.assertEqual(mut.id, other_mut.id)
            self.assertEqual(mut.site, other_mut.site)
            self.assertEqual(mut.parent, other_mut.parent)
            self.assertEqual(mut.node, other_mut.node)
            self.assertEqual(mut.metadata, other_mut.metadata)
            # Check the deprecated attrs.
            self.assertEqual(mut.position, ts.site(mut.site).position)
            self.assertEqual(mut.index, mut.site)

    def test_sites_mutations(self):
        # Check that the mutations iterator returns the correct values.
        for ts in get_example_tree_sequences():
            self.verify_mutations(ts)

    def test_removed_methods(self):
        ts = next(get_example_tree_sequences())
        self.assertRaises(NotImplementedError, ts.get_num_records)
        self.assertRaises(NotImplementedError, ts.diffs)
        self.assertRaises(NotImplementedError, ts.newick_trees)

    def test_zlib_compression_warning(self):
        ts = msprime.simulate(5, random_seed=1)
        with warnings.catch_warnings(record=True) as w:
            ts.dump(self.temp_file, zlib_compression=True)
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[0].category, RuntimeWarning))
        with warnings.catch_warnings(record=True) as w:
            ts.dump(self.temp_file, zlib_compression=False)
            self.assertEqual(len(w), 0)

    def test_tables_sequence_length_round_trip(self):
        for sequence_length in [0.1, 1, 10, 100]:
            ts = msprime.simulate(5, length=sequence_length, random_seed=1)
            self.assertEqual(ts.sequence_length, sequence_length)
            tables = ts.tables
            self.assertEqual(tables.sequence_length, sequence_length)
            new_ts = tables.tree_sequence()
            self.assertEqual(new_ts.sequence_length, sequence_length)


class TestFileUuid(HighLevelTestCase):
    """
    Tests that the file UUID attribute is handled correctly.
    """
    def validate(self, ts):
        self.assertIsNone(ts.file_uuid)
        ts.dump(self.temp_file)
        other_ts = msprime.load(self.temp_file)
        self.assertIsNotNone(other_ts.file_uuid)
        self.assertTrue(len(other_ts.file_uuid), 36)
        uuid = other_ts.file_uuid
        other_ts = msprime.load(self.temp_file)
        self.assertEqual(other_ts.file_uuid, uuid)
        self.assertEqual(ts.tables, other_ts.tables)

        # Check that the UUID is well-formed.
        parsed = _uuid.UUID("{" + uuid + "}")
        self.assertEqual(str(parsed), uuid)

        # Save the same tree sequence to the file. We should get a different UUID.
        ts.dump(self.temp_file)
        other_ts = msprime.load(self.temp_file)
        self.assertIsNotNone(other_ts.file_uuid)
        self.assertNotEqual(other_ts.file_uuid, uuid)

        # Even saving a ts that has a UUID to another file changes the UUID
        old_uuid = other_ts.file_uuid
        other_ts.dump(self.temp_file)
        self.assertEqual(other_ts.file_uuid, old_uuid)
        other_ts = msprime.load(self.temp_file)
        self.assertIsNotNone(other_ts.file_uuid)
        self.assertNotEqual(other_ts.file_uuid, old_uuid)

        # Tables dumped from this ts are a deep copy, so they don't have
        # the file_uuid.
        tables = other_ts.dump_tables()
        self.assertIsNone(tables.file_uuid)

        # For now, ts.tables also returns a deep copy. This will hopefully
        # change in the future thoug.
        self.assertIsNone(ts.tables.file_uuid)

    def test_simple_simulation(self):
        ts = msprime.simulate(2, random_seed=1)
        self.validate(ts)

    def test_empty_tables(self):
        tables = msprime.TableCollection(1)
        self.validate(tables.tree_sequence())


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
            ["id", "is_sample", "time", "population", "individual", "metadata"])
        for node, line in zip(ts.nodes(), output_nodes[1:]):
            splits = line.split("\t")
            self.assertEqual(str(node.id), splits[0])
            self.assertEqual(str(node.is_sample()), splits[1])
            self.assertEqual(convert(node.time), splits[2])
            self.assertEqual(str(node.population), splits[3])
            self.assertEqual(str(node.individual), splits[4])
            self.assertEqual(tests.base64_encode(node.metadata), splits[5])

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
            ["position", "ancestral_state", "metadata"])
        for site, line in zip(ts.sites(), output_sites[1:]):
            splits = line.split("\t")
            self.assertEqual(convert(site.position), splits[0])
            self.assertEqual(site.ancestral_state, splits[1])
            self.assertEqual(tests.base64_encode(site.metadata), splits[2])

    def verify_mutations_format(self, ts, mutations_file, precision):
        """
        Verifies that the mutations we output have the correct form.
        """
        def convert(v):
            return "{:.{}f}".format(v, precision)
        output_mutations = mutations_file.read().splitlines()
        self.assertEqual(len(output_mutations) - 1, ts.num_mutations)
        self.assertEqual(
            list(output_mutations[0].split()),
            ["site", "node", "derived_state", "parent", "metadata"])
        mutations = [mut for site in ts.sites() for mut in site.mutations]
        for mutation, line in zip(mutations, output_mutations[1:]):
            splits = line.split("\t")
            self.assertEqual(str(mutation.site), splits[0])
            self.assertEqual(str(mutation.node), splits[1])
            self.assertEqual(str(mutation.derived_state), splits[2])
            self.assertEqual(str(mutation.parent), splits[3])
            self.assertEqual(tests.base64_encode(mutation.metadata), splits[4])

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
            self.assertEqual(n1.metadata, n2.metadata)
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
            self.assertEqual(s1.metadata, s2.metadata)
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
            edges_file = six.StringIO()
            sites_file = six.StringIO()
            mutations_file = six.StringIO()
            individuals_file = six.StringIO()
            populations_file = six.StringIO()
            tables = ts1.dump_tables()
            ts1 = msprime.load_tables(
                sequence_length=ts1.sequence_length, **tables.asdict())
            ts1.dump_text(
                nodes=nodes_file, edges=edges_file, sites=sites_file,
                mutations=mutations_file, individuals=individuals_file,
                populations=populations_file, precision=16)
            nodes_file.seek(0)
            edges_file.seek(0)
            sites_file.seek(0)
            mutations_file.seek(0)
            individuals_file.seek(0)
            populations_file.seek(0)
            ts2 = msprime.load_text(
                nodes=nodes_file, edges=edges_file, sites=sites_file,
                mutations=mutations_file, individuals=individuals_file,
                populations=populations_file,
                sequence_length=ts1.sequence_length,
                strict=True)
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

    def verify_mutations(self, tree):
        self.assertGreater(tree.num_mutations, 0)
        other_mutations = []
        for site in tree.sites():
            for mutation in site.mutations:
                other_mutations.append(mutation)
        mutations = list(tree.mutations())
        self.assertEqual(tree.num_mutations, len(other_mutations))
        self.assertEqual(tree.num_mutations, len(mutations))
        for mut, other_mut in zip(mutations, other_mutations):
            # We cannot compare these directly as the mutations obtained
            # from the mutations iterator will have extra deprecated
            # attributes.
            self.assertEqual(mut.id, other_mut.id)
            self.assertEqual(mut.site, other_mut.site)
            self.assertEqual(mut.parent, other_mut.parent)
            self.assertEqual(mut.node, other_mut.node)
            self.assertEqual(mut.metadata, other_mut.metadata)
            # Check the deprecated attrs.
            self.assertEqual(mut.position, tree.tree_sequence.site(mut.site).position)
            self.assertEqual(mut.index, mut.site)

    def test_simple_mutations(self):
        tree = self.get_tree()
        self.verify_mutations(tree)

    def test_complex_mutations(self):
        ts = tsutil.insert_branch_mutations(msprime.simulate(10, random_seed=1))
        self.verify_mutations(ts.first())

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
        # TODO to make this work we may need to clamp the precision of node
        # times because Python and C float printing algorithms work slightly
        # differently. Seems to work OK now, so leaving alone.
        if tree.num_roots == 1:
            py_tree = tests.PythonSparseTree.from_sparse_tree(tree)
            newick1 = tree.newick(precision=16)
            newick2 = py_tree.newick()
            self.assertEqual(newick1, newick2)

            # Make sure we get the same results for a leaf root.
            newick1 = tree.newick(root=0, precision=16)
            newick2 = py_tree.newick(root=0)
            self.assertEqual(newick1, newick2)

            # When we specify the node_labels we should get precisely the
            # same result as we are using Python code now.
            for precision in [0, 3, 19]:
                newick1 = tree.newick(precision=precision, node_labels={})
                newick2 = py_tree.newick(precision=precision, node_labels={})
                self.assertEqual(newick1, newick2)
        else:
            self.assertRaises(ValueError, tree.newick)
            for root in tree.roots:
                py_tree = tests.PythonSparseTree.from_sparse_tree(tree)
                newick1 = tree.newick(precision=16, root=root)
                newick2 = py_tree.newick(root=root)
                self.assertEqual(newick1, newick2)

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
            self.assertEqual(len(sim.samples), n)
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
            self.assertEqual(len(sim.samples), sample_size)
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
            self.assertEqual(len(sim.samples), 2 * d)
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
            # Try with equivalent numpy array.
            sim = f(np.array(hl_matrix))
            self.assertEqual(sim.migration_matrix, hl_matrix)
            ll_sim = sim.create_ll_instance()
            ll_matrix = [v for row in hl_matrix for v in row]
            self.assertEqual(ll_sim.get_migration_matrix(), ll_matrix)
            for bad_type in [234, 1.2]:
                self.assertRaises(TypeError, f, bad_type)
            # Iterables should raise a value error.
            for bad_type in [{}, ""]:
                self.assertRaises(ValueError, f, bad_type)
            # Now check for the structure of the matrix.
            hl_matrix[0][0] = "bad value"
            sim = f(hl_matrix)
            self.assertRaises(TypeError, sim.create_ll_instance)
            hl_matrix[0] = None
            self.assertRaises(TypeError, f, hl_matrix)
            hl_matrix[0] = []
            self.assertRaises(ValueError, f, hl_matrix)
            # Simple numpy array.
            hl_matrix = np.ones((N, N))
            np.fill_diagonal(hl_matrix, 0)
            sim = f(hl_matrix)
            self.assertTrue(np.array_equal(np.array(sim.migration_matrix), hl_matrix))

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
        self.assertEqual(len(list(ts.provenances())), 1)

    def test_numpy_random_seed(self):
        seed = np.array([12345], dtype=np.int64)[0]
        self.assertEqual(seed.dtype, np.int64)
        ts1 = msprime.simulate(10, random_seed=seed)
        ts2 = msprime.simulate(10, random_seed=seed)
        self.assertEqual(ts1.tables.nodes, ts2.tables.nodes)

    def verify_provenance(self, provenance):
        """
        Checks that the specified provenance object has the right sort of
        properties.
        """
        # Generate the ISO 8601 time for now, without the high precision suffix,
        # and compare the prefixes.
        today = datetime.date.today().isoformat()
        k = len(today)
        self.assertEqual(provenance.timestamp[:k], today)
        self.assertEqual(provenance.timestamp[k], "T")
        d = json.loads(provenance.record)
        self.assertGreater(len(d), 0)
        # TODO check the format of the dictionary.

    def test_provenance(self):
        ts = msprime.simulate(10)
        self.assertEqual(ts.num_provenances, 1)
        self.verify_provenance(ts.provenance(0))
        # TODO check the form of the dictionary
        for ts in msprime.simulate(10, num_replicates=10):
            self.assertEqual(ts.num_provenances, 1)
            self.verify_provenance(ts.provenance(0))

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

    def test_no_mutations_with_start_time(self):
        with self.assertRaises(ValueError):
            msprime.simulate(10, mutation_rate=10, start_time=3)
        # But fine if we set start_time = None
        ts = msprime.simulate(10, mutation_rate=10, start_time=None, random_seed=1)
        self.assertGreater(ts.num_sites, 0)

    def test_mutation_generator_unsupported(self):
        n = 10
        mutgen = msprime.MutationGenerator(msprime.RandomGenerator(1), 1)
        with self.assertRaises(ValueError):
            msprime.simulate(n, mutation_generator=mutgen)

    def test_mutation_interface(self):
        for bad_type in ["x", [], {}]:
            self.assertRaises(
                TypeError, msprime.simulate, 10, mutation_rate=bad_type)

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
            self.assertEqual(n1.metadata, n2.metadata)
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
        other_tables = msprime.TableCollection(ts.sequence_length)
        # Insert the new nodes into the table.
        inv_node_map = {v: k for k, v in node_map.items()}
        for j in range(ts.num_nodes):
            node = ts.node(inv_node_map[j])
            other_tables.nodes.add_row(
                flags=node.flags, time=node.time, population=node.population)
        for e in ts.edges():
            other_tables.edges.add_row(
                left=e.left, right=e.right, parent=node_map[e.parent],
                child=node_map[e.child])
        for _ in range(ts.num_populations):
            other_tables.populations.add_row()
        other_tables.sort()
        other_ts = other_tables.tree_sequence()

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
                msprime.SimpleBottleneck(time=0.5, population=0, proportion=1)])
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


# Convenience method for getting seeds in a subprocess.
def _get_seed(x):
    return msprime.simulations._get_random_seed()


class TestDefaultRandomSeeds(unittest.TestCase):
    """
    Tests for the default random seed generator.
    """

    def test_seed_generator_init(self):
        msprime.simulations._clear_seed_rng()
        seed = msprime.simulations._get_random_seed()
        self.assertGreater(seed, 0)
        self.assertIsNotNone(msprime.simulations._get_seed_rng())

    def test_unique(self):
        n = 100
        msprime.simulations._clear_seed_rng()
        seeds1 = [msprime.simulations._get_random_seed() for _ in range(n)]
        self.assertEqual(len(set(seeds1)), n)
        seeds2 = [msprime.simulations._get_random_seed() for _ in range(n)]
        self.assertEqual(len(set(seeds2)), n)
        self.assertEqual(len(set(seeds2)) + len(set(seeds2)), 2 * n)

    def test_unique_multiple_processes_no_init(self):
        n = 100
        msprime.simulations._clear_seed_rng()
        # Would use with block here, but not supported in Py < 3.3.
        pool = multiprocessing.Pool(5)
        seeds = pool.map(_get_seed, range(n))
        self.assertEqual(len(set(seeds)), n)
        pool.terminate()
        pool.join()

    def test_unique_multiple_processes_init(self):
        n = 100
        msprime.simulations._get_random_seed()
        self.assertIsNotNone(msprime.simulations._get_seed_rng())
        # Would use with block here, but not supported in Py < 3.3.
        pool = multiprocessing.Pool(5)
        seeds = pool.map(_get_seed, range(n))
        self.assertEqual(len(set(seeds)), n)
        pool.terminate()
        pool.join()
