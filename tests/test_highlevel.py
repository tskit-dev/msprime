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
import datetime
import json
import os
import random
import shutil
import sys
import tempfile
import unittest
import multiprocessing

import numpy as np

import msprime
import tests


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


# TODO Not sure why we need this method in msprime. Remove?
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


class TestFullArg(unittest.TestCase):
    """
    Tests for recording the full ARG.
    """
    def verify(self, sim, multiple_mergers=False):
        sim.run()
        tree_sequence = sim.get_tree_sequence()
        # Check if we have multiple merger somewhere.
        found = False
        for edgeset in tree_sequence.edgesets():
            if len(edgeset.children) > 2:
                found = True
                break
        self.assertEqual(multiple_mergers, found)

        flags = tree_sequence.tables.nodes.flags
        time = tree_sequence.tables.nodes.time
        # TODO add checks for migrations.
        re_nodes = np.where(flags == msprime.NODE_IS_RE_EVENT)[0]
        ca_nodes = np.where(flags == msprime.NODE_IS_CA_EVENT)[0]
        coal_nodes = np.where(flags == 0)[0]
        # There should be two recombination nodes for every event
        self.assertTrue(np.array_equal(
            time[re_nodes[::2]],  # Even indexes
            time[re_nodes[1::2]]))  # Odd indexes
        self.assertEqual(re_nodes.shape[0] / 2,  sim.num_recombination_events)
        if not multiple_mergers:
            self.assertEqual(
                ca_nodes.shape[0] + coal_nodes.shape[0], sim.num_common_ancestor_events)
        # After simplification, all the RE and CA nodes should be gone.
        ts_simplified = tree_sequence.simplify()
        new_flags = ts_simplified.tables.nodes.flags
        new_time = ts_simplified.tables.nodes.time
        self.assertEqual(np.sum(new_flags == msprime.NODE_IS_RE_EVENT), 0)
        self.assertEqual(np.sum(new_flags == msprime.NODE_IS_CA_EVENT), 0)
        # All coal nodes from the original should be identical to the originals
        self.assertTrue(np.array_equal(time[coal_nodes], new_time[new_flags == 0]))
        self.assertLessEqual(ts_simplified.num_nodes, tree_sequence.num_nodes)
        self.assertLessEqual(ts_simplified.num_edges, tree_sequence.num_edges)
        return tree_sequence

    def test_no_recombination(self):
        rng = msprime.RandomGenerator(1)
        sim = msprime.simulator_factory(10, random_generator=rng, record_full_arg=True)
        ts = self.verify(sim)
        ts_simplified = ts.simplify()
        t1 = ts.tables
        t2 = ts_simplified.tables
        self.assertEqual(t1.nodes, t2.nodes)
        self.assertEqual(t1.edges, t2.edges)

    def test_recombination_n25(self):
        rng = msprime.RandomGenerator(10)
        sim = msprime.simulator_factory(
            25, recombination_rate=1, record_full_arg=True, random_generator=rng)
        self.verify(sim)

    def test_recombination_n5(self):
        rng = msprime.RandomGenerator(10)
        sim = msprime.simulator_factory(
            5, recombination_rate=10, record_full_arg=True, random_generator=rng)
        self.verify(sim)

    def test_recombination_n50(self):
        rng = msprime.RandomGenerator(100)
        sim = msprime.simulator_factory(
            50, recombination_rate=2, record_full_arg=True, random_generator=rng)
        self.verify(sim)

    def test_recombination_n100(self):
        rng = msprime.RandomGenerator(100)
        sim = msprime.simulator_factory(
            100, recombination_rate=0.2, record_full_arg=True, random_generator=rng)
        self.verify(sim)

    def test_multimerger(self):
        rng = msprime.RandomGenerator(1234)
        sim = msprime.simulator_factory(
            100, recombination_rate=0.1, record_full_arg=True,
            random_generator=rng, demographic_events=[
                msprime.InstantaneousBottleneck(time=0.1, population=0, strength=5)])
        self.verify(sim, multiple_mergers=True)


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

    def test_length(self):
        for bad_length in [-1, 0, -1e-6]:
            with self.assertRaises(ValueError):
                msprime.simulator_factory(10, length=bad_length)

    def test_num_labels(self):
        for bad_value in [-1, 0, 0.1]:
            with self.assertRaises(ValueError):
                msprime.simulator_factory(10, num_labels=bad_value)

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
            self.assertEqual(sim.model.reference_size, Ne)
        # Test the default.
        sim = msprime.simulator_factory(10)
        self.assertEqual(sim.model.reference_size, 1)

    def test_population_configurations(self):
        def f(configs):
            return msprime.simulator_factory(population_configurations=configs)
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
            sim.run()
            events = np.array(sim.num_migration_events)
            self.assertEqual(events.shape, (N, N))
            self.assertTrue(np.all(events >= 0))

    def test_default_migration_matrix(self):
        sim = msprime.simulator_factory(10)
        ll_sim = sim.create_ll_instance()
        self.assertEqual(ll_sim.get_migration_matrix(), [0.0])

    def test_demographic_events(self):
        for bad_type in ["sdf", 234, [12], [None]]:
            self.assertRaises(
                TypeError, msprime.simulator_factory, 2, demographic_events=bad_type)
        # TODO test for bad values.

    def test_recombination_rate(self):
        def f(recomb_rate):
            return msprime.simulator_factory(10, recombination_rate=recomb_rate)
        for bad_type in ["", {}, []]:
            self.assertRaises(TypeError, f, bad_type)
        for bad_value in [-1, -1e15]:
            self.assertRaises(ValueError, f, bad_value)
        for rate in [0, 1e-3, 10]:
            sim = f(rate)
            recomb_map = sim.recombination_map
            self.assertEqual(recomb_map.get_positions(), [0, 1], [rate, 0])
            self.assertEqual(
                recomb_map.get_num_loci(), msprime.RecombinationMap.DEFAULT_NUM_LOCI)
            self.assertEqual(sim.num_loci, recomb_map.get_num_loci())

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
            self.assertAlmostEqual(
                sim.recombination_map.get_per_locus_recombination_rate(),
                per_locus_rate)

    def test_recombination_map(self):
        def f(recomb_map):
            return msprime.simulator_factory(10, recombination_map=recomb_map)
        self.assertRaises(TypeError, f, "wrong type")
        for n in range(2, 10):
            positions = list(range(n))
            rates = [0.1 * j for j in range(n - 1)] + [0.0]
            recomb_map = msprime.RecombinationMap(positions, rates)
            sim = msprime.simulator_factory(10, recombination_map=recomb_map)
            self.assertEqual(sim.recombination_map, recomb_map)
            self.assertEqual(recomb_map.get_positions(), positions)
            self.assertEqual(recomb_map.get_rates(), rates)
            self.assertEqual(sim.num_loci, recomb_map.get_num_loci())
            ll_sim = sim.create_ll_instance()
            self.assertEqual(ll_sim.get_num_loci(), recomb_map.get_num_loci())

    def test_zero_recombination_map(self):
        # test that beginning and trailing zero recombination regions in the
        # recomb map are included in the sequence
        for n in range(3, 10):
            positions = list(range(n))
            rates = [0.0, 0.2] + [0.0] * (n - 2)
            recomb_map = msprime.RecombinationMap(positions, rates)
            ts = msprime.simulate(10, recombination_map=recomb_map)
            self.assertEqual(ts.sequence_length, n - 1)
            self.assertEqual(min(ts.tables.edges.left), 0.0)
            self.assertEqual(max(ts.tables.edges.right), n - 1.0)

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

    def test_mean_recombination_rate(self):
        # Some quick sanity checks.
        recomb_map = msprime.RecombinationMap([0, 1], [1, 0])
        mean_rr = recomb_map.mean_recombination_rate
        self.assertEqual(mean_rr, 1.0)

        recomb_map = msprime.RecombinationMap([0, 1, 2], [1, 0, 0])
        mean_rr = recomb_map.mean_recombination_rate
        self.assertEqual(mean_rr, 0.5)

        recomb_map = msprime.RecombinationMap([0, 1, 2], [0, 0, 0])
        mean_rr = recomb_map.mean_recombination_rate
        self.assertEqual(mean_rr, 0.0)

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

    def test_num_labels(self):
        # Running simulations with different numbers of labels in the default
        # setting should have no effect.
        tables = [
            msprime.simulate(10, num_labels=num_labels, random_seed=1).tables
            for num_labels in range(1, 5)]
        for t in tables:
            t.provenances.clear()
        for t in tables:
            self.assertEqual(t, tables[0])


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
