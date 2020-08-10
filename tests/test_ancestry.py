#
# Copyright (C) 2015-2020 University of Oxford
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
Test cases for basic ancestry simulation operations.
"""
import datetime
import json
import logging
import os
import random
import tempfile
import unittest
import warnings

import numpy as np
import tskit

import _msprime
import msprime


def has_discrete_genome(ts):
    """
    Returns True if the specified tree sequence has discrete genome coordinates.
    """
    tables = ts.tables
    edges_left = np.all(tables.edges.left == np.floor(tables.edges.left))
    edges_right = np.all(tables.edges.right == np.floor(tables.edges.right))
    migrations_left = np.all(tables.migrations.left == np.floor(tables.migrations.left))
    migrations_right = np.all(
        tables.migrations.right == np.floor(tables.migrations.right)
    )
    sites = np.all(tables.sites.position == np.floor(tables.sites.position))
    return edges_left and edges_right and migrations_left and migrations_right and sites


def get_bottleneck_examples():
    """
    Returns an iterator of example tree sequences with nonbinary
    trees.
    """
    bottlenecks = [
        msprime.SimpleBottleneck(0.01, 0, proportion=0.05),
        msprime.SimpleBottleneck(0.02, 0, proportion=0.25),
        msprime.SimpleBottleneck(0.03, 0, proportion=1),
    ]
    for n in [3, 10, 100]:
        ts = msprime.simulate(
            n,
            length=100,
            recombination_rate=1,
            demographic_events=bottlenecks,
            random_seed=n,
        )
        yield ts


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
        self.assertTrue(
            np.array_equal(time[re_nodes[::2]], time[re_nodes[1::2]])  # Even indexes
        )  # Odd indexes
        self.assertEqual(re_nodes.shape[0] / 2, sim.num_recombination_events)
        if not multiple_mergers:
            self.assertEqual(
                ca_nodes.shape[0] + coal_nodes.shape[0], sim.num_common_ancestor_events
            )
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
        rng = _msprime.RandomGenerator(1)
        sim = msprime.simulator_factory(10, random_generator=rng, record_full_arg=True)
        ts = self.verify(sim)
        ts_simplified = ts.simplify()
        t1 = ts.tables
        t2 = ts_simplified.tables
        self.assertEqual(t1.nodes, t2.nodes)
        self.assertEqual(t1.edges, t2.edges)

    def test_recombination_n25(self):
        rng = _msprime.RandomGenerator(10)
        sim = msprime.simulator_factory(
            25, recombination_rate=1, record_full_arg=True, random_generator=rng
        )
        self.verify(sim)

    def test_recombination_n5(self):
        rng = _msprime.RandomGenerator(10)
        sim = msprime.simulator_factory(
            5, recombination_rate=10, record_full_arg=True, random_generator=rng
        )
        self.verify(sim)

    def test_recombination_n50(self):
        rng = _msprime.RandomGenerator(100)
        sim = msprime.simulator_factory(
            50, recombination_rate=2, record_full_arg=True, random_generator=rng
        )
        self.verify(sim)

    def test_recombination_n100(self):
        rng = _msprime.RandomGenerator(100)
        sim = msprime.simulator_factory(
            100, recombination_rate=0.2, record_full_arg=True, random_generator=rng
        )
        self.verify(sim)

    def test_multimerger(self):
        rng = _msprime.RandomGenerator(1234)
        sim = msprime.simulator_factory(
            100,
            recombination_rate=0.1,
            record_full_arg=True,
            random_generator=rng,
            demographic_events=[
                msprime.InstantaneousBottleneck(time=0.1, population=0, strength=5)
            ],
        )
        self.verify(sim, multiple_mergers=True)


class TestSimulator(unittest.TestCase):
    """
    Runs tests on the underlying Simulator object.
    """

    def verify_simulation(self, n, m, r):
        """
        Verifies a simulation for the specified parameters.
        """
        recomb_map = msprime.RecombinationMap.uniform_map(m, r, discrete=True)
        rng = _msprime.RandomGenerator(1)
        sim = msprime.simulator_factory(
            n, recombination_map=recomb_map, random_generator=rng
        )
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
        self.assertGreaterEqual(np.sum(sim.num_migration_events), 0)
        self.assertGreaterEqual(sim.num_multiple_recombination_events, 0)

    def test_random_parameters(self):
        num_random_sims = 10
        for _ in range(num_random_sims):
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

    def test_event_chunk(self):
        sim = msprime.simulator_factory(10)
        for bad_chunk in [-(2 ** 32), -1, 0]:
            with self.assertRaises(ValueError):
                sim.run(event_chunk=bad_chunk)
        sim.reset()
        sim.run(event_chunk=2 ** 32 + 1)
        sim.reset()
        sim.run(event_chunk=2 ** 64 + 1)

    def test_debug_func(self):
        sim = msprime.simulator_factory(10)
        count = 0

        def f(sim):
            nonlocal count
            count += 1

        sim.run(event_chunk=1, debug_func=f)
        self.assertGreater(count, 0)

    def test_info_logging(self):
        sim = msprime.simulator_factory(10)
        with self.assertLogs("msprime.ancestry", logging.INFO) as log:
            sim.run()
            self.assertEqual(len(log.output), 2)
            self.assertEqual(
                log.output[0],
                (
                    "INFO:msprime.ancestry:Running model {'name': 'hudson'} "
                    "until max time: inf"
                ),
            )
            self.assertTrue(
                log.output[1].startswith("INFO:msprime.ancestry:Completed at time")
            )

    def test_debug_logging(self):
        sim = msprime.simulator_factory(3)
        with self.assertLogs("msprime.ancestry", logging.DEBUG) as log:
            sim.run(event_chunk=1)
            self.assertEqual(len(log.output), 3)
            self.assertTrue(log.output[0].startswith("INFO"))
            self.assertTrue(log.output[-1].startswith("INFO"))
            self.assertTrue(log.output[1].startswith("DEBUG:msprime.ancestry:time="))

    def test_debug_logging_dtwf(self):
        sim = msprime.simulator_factory(3, Ne=10, model="dtwf")
        with self.assertLogs("msprime.ancestry", logging.DEBUG) as log:
            sim.run(event_chunk=1)
            self.assertGreaterEqual(len(log.output), 3)
            self.assertTrue(log.output[0].startswith("INFO"))
            self.assertTrue(log.output[-1].startswith("INFO"))
            self.assertTrue(log.output[1].startswith("DEBUG:msprime.ancestry:time="))


class TestDemographyFactory(unittest.TestCase):
    """
    Tests fo the demography_factory function.
    """

    def test_mixed_old_and_new_style(self):
        demography = msprime.Demography()

        def f(
            population_configurations=None,
            migration_matrix=None,
            demographic_events=None,
        ):
            msprime.demography_factory(
                Ne=1,
                demography=demography,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix,
                demographic_events=demographic_events,
            )

        with self.assertRaises(ValueError):
            f(population_configurations=[])
        with self.assertRaises(ValueError):
            f(migration_matrix=[[]])
        with self.assertRaises(ValueError):
            f(demographic_events=[])

    def test_input_demography_copied(self):
        d1 = msprime.Demography.island_model(2, 1, Ne=100)
        d2 = msprime.demography_factory(
            Ne=None,
            demography=d1,
            population_configurations=None,
            migration_matrix=None,
            demographic_events=None,
        )
        self.assertEqual(d1, d2)
        self.assertIsNot(d1, d2)
        self.assertIsNot(d1.populations[0], d2.populations[0])
        self.assertIsNot(d1.populations[1], d2.populations[1])
        self.assertIsNot(d1.migration_matrix, d2.migration_matrix)

    def test_Ne_does_not_override_demography(self):
        d1 = msprime.Demography.island_model(2, 1, Ne=100)
        self.assertEqual(d1.populations[0].initial_size, 100)
        self.assertEqual(d1.populations[1].initial_size, 100)
        d2 = msprime.demography_factory(
            Ne=1234,
            demography=d1,
            population_configurations=None,
            migration_matrix=None,
            demographic_events=None,
        )
        self.assertEqual(d2.populations[0].initial_size, 100)
        self.assertEqual(d2.populations[1].initial_size, 100)

    def test_Ne_overwrites_size_none(self):
        d1 = msprime.Demography.island_model(2, 1, Ne=None)
        self.assertEqual(d1.populations[0].initial_size, None)
        self.assertEqual(d1.populations[1].initial_size, None)
        d2 = msprime.demography_factory(
            Ne=1234,
            demography=d1,
            population_configurations=None,
            migration_matrix=None,
            demographic_events=None,
        )
        self.assertEqual(d2.populations[0].initial_size, 1234)
        self.assertEqual(d2.populations[1].initial_size, 1234)

        d1.populations[0].initial_size = 100
        d1.populations[1].initial_size = None
        d2 = msprime.demography_factory(
            Ne=1234,
            demography=d1,
            population_configurations=None,
            migration_matrix=None,
            demographic_events=None,
        )
        self.assertEqual(d2.populations[0].initial_size, 100)
        self.assertEqual(d2.populations[1].initial_size, 1234)


class TestSimulatorFactory(unittest.TestCase):
    """
    Tests that the simulator factory high-level function correctly
    creates simulators with the required parameter values.
    """

    def test_default_random_seed(self):
        sim = msprime.simulator_factory(10)
        rng = sim.random_generator
        self.assertIsInstance(rng, _msprime.RandomGenerator)
        self.assertNotEqual(rng.get_seed(), 0)

    def test_random_generator(self):
        seed = 12345
        rng = _msprime.RandomGenerator(seed)
        sim = msprime.simulator_factory(10, random_generator=rng)
        self.assertEqual(rng, sim.random_generator)
        self.assertEqual(rng.get_seed(), seed)

    def test_random_seed(self):
        seed = 12345
        sim = msprime.simulator_factory(10, random_seed=seed)
        self.assertEqual(sim.random_generator.get_seed(), seed)

        # It's an error to specify both seed and generator.
        with self.assertRaises(ValueError):
            msprime.simulator_factory(
                10, random_seed=1234, random_generator=_msprime.RandomGenerator(1234)
            )

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
        self.assertRaises(ValueError, msprime.simulator_factory, sample_size=1)
        for n in [2, 100, 1000]:
            sim = msprime.simulator_factory(n)
            self.assertEqual(sim.num_samples, n)
            self.assertEqual(len(sim.samples), n)
            for sample in sim.samples:
                self.assertEqual(sample[0], 0)
                self.assertEqual(sample[1], 0)

    def test_effective_population_size(self):
        def f(Ne):
            return msprime.simulator_factory(10, Ne=Ne)

        for bad_value in [-1, -1e16, 0]:
            self.assertRaises(ValueError, f, bad_value)
        for Ne in [1, 10, 1e5]:
            sim = f(Ne)
            self.assertEqual(sim.demography.populations[0].initial_size, Ne)
        # Test the default.
        sim = msprime.simulator_factory(10)
        self.assertEqual(sim.demography.populations[0].initial_size, 1)

    def test_discrete_genome_continuous_length(self):
        for bad_length in [0.1, 1.1, 1000.1]:
            with self.assertRaises(ValueError):
                msprime.simulator_factory(10, discrete_genome=True, length=bad_length)

    def test_population_configurations(self):
        def f(configs):
            return msprime.simulator_factory(population_configurations=configs)

        for bad_type in [10, ["sdf"], "sdfsd"]:
            self.assertRaises(TypeError, f, bad_type)
        # Just test the basic equalities here. The actual
        # configuration options are tested elewhere.
        for N in range(1, 10):
            pop_configs = [
                msprime.PopulationConfiguration(5, initial_size=5) for _ in range(N)
            ]
            sample_size = 5 * N
            sim = msprime.simulator_factory(population_configurations=pop_configs)
            self.assertEqual(len(sim.demography.populations), len(pop_configs))
            for pop, pop_config in zip(sim.demography.populations, pop_configs):
                self.assertEqual(pop.initial_size, pop_config.initial_size)
                self.assertEqual(pop.growth_rate, pop_config.growth_rate)
            self.assertEqual(len(sim.samples), sample_size)
            self.assertEqual(len(sim.population_configuration), N)
        # The default is a single population
        sim = msprime.simulator_factory(10)
        self.assertEqual(len(sim.population_configuration), 1)

    def test_sample_size_population_configuration(self):
        for d in range(1, 5):
            # Zero sample size is always an error
            configs = [msprime.PopulationConfiguration(0) for _ in range(d)]
            self.assertRaises(
                ValueError, msprime.simulator_factory, population_configurations=configs
            )
            configs = [msprime.PopulationConfiguration(2) for _ in range(d)]
            sim = msprime.simulator_factory(population_configurations=configs)
            self.assertEqual(len(sim.samples), 2 * d)
            samples = []
            for j in range(d):
                samples += [msprime.Sample(population=j, time=0) for _ in range(2)]
            self.assertEqual(sim.samples, samples)

    def test_migration_matrix(self):
        # Cannot specify a migration matrix without population
        # configurations
        self.assertRaises(
            ValueError, msprime.simulator_factory, 10, migration_matrix=[]
        )
        for N in range(1, 10):
            pop_configs = [msprime.PopulationConfiguration(5) for _ in range(N)]
            sim = msprime.simulator_factory(population_configurations=pop_configs)
            # If we don't specify a matrix, it's 0 everywhere.
            matrix = np.zeros((N, N))
            np.testing.assert_array_equal(sim.migration_matrix, matrix)

            def f(matrix):
                return msprime.simulator_factory(
                    population_configurations=pop_configs, migration_matrix=matrix
                )

            matrix = [[(j + k) * int(j != k) for j in range(N)] for k in range(N)]
            sim = f(matrix)
            np.testing.assert_array_equal(sim.demography.migration_matrix, matrix)
            # Try with equivalent numpy array.
            sim = f(np.array(matrix))
            np.testing.assert_array_equal(sim.demography.migration_matrix, matrix)
            np.testing.assert_array_equal(sim.migration_matrix, matrix)
            for bad_type in [{}, "", 234, 1.2]:
                self.assertRaises(ValueError, f, bad_type)
            # Now check for the structure of the matrix.
            matrix[0][0] = "bad value"
            self.assertRaises(ValueError, f, matrix)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                matrix[0] = None
                self.assertRaises(ValueError, f, matrix)
                matrix[0] = []
                self.assertRaises(ValueError, f, matrix)
            # Simple numpy array.
            matrix = np.ones((N, N))
            np.fill_diagonal(matrix, 0)
            sim = f(matrix)
            np.testing.assert_array_equal(
                np.array(sim.demography.migration_matrix), matrix
            )
            sim.run()
            events = np.array(sim.num_migration_events)
            self.assertEqual(events.shape, (N, N))
            self.assertTrue(np.all(events >= 0))

    def test_default_migration_matrix(self):
        sim = msprime.simulator_factory(10)
        self.assertEqual(sim.migration_matrix, [0.0])

    def test_demographic_events(self):
        for bad_type in ["sdf", 234, [12], [None]]:
            self.assertRaises(
                TypeError, msprime.simulator_factory, 2, demographic_events=bad_type
            )
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
            self.assertEqual(sim.sequence_length, recomb_map.get_sequence_length())

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
            self.assertEqual(sim.sequence_length, recomb_map.get_sequence_length())

    def test_combining_recomb_map_and_rate_length(self):
        recomb_map = msprime.RecombinationMap([0, 1], [1, 0])
        self.assertRaises(
            ValueError,
            msprime.simulator_factory,
            10,
            recombination_map=recomb_map,
            length=1,
        )
        self.assertRaises(
            ValueError,
            msprime.simulator_factory,
            10,
            recombination_map=recomb_map,
            recombination_rate=100,
        )
        self.assertRaises(
            ValueError,
            msprime.simulator_factory,
            10,
            recombination_map=recomb_map,
            length=1,
            recombination_rate=1,
        )

    def test_sample_combination_errors(self):
        # Make sure that the various ways we can specify the samples
        # operate correctly.
        s = msprime.Sample(time=0.0, population=0)
        self.assertRaises(ValueError, msprime.simulator_factory)
        # Cannot provide sample_size with either population configurations
        # or samples
        self.assertRaises(
            ValueError, msprime.simulator_factory, sample_size=2, samples=[s, s]
        )
        pop_configs = [msprime.PopulationConfiguration(sample_size=2)]
        self.assertRaises(
            ValueError,
            msprime.simulator_factory,
            sample_size=2,
            population_configurations=pop_configs,
        )
        # If we provide samples and population_configurations we cannot
        # have a sample size for the config.
        pop_configs = [msprime.PopulationConfiguration(sample_size=2)]
        self.assertRaises(
            ValueError,
            msprime.simulator_factory,
            samples=[s, s],
            population_configurations=pop_configs,
        )
        pop_configs = [
            msprime.PopulationConfiguration(sample_size=None),
            msprime.PopulationConfiguration(sample_size=2),
        ]
        self.assertRaises(
            ValueError,
            msprime.simulator_factory,
            samples=[s, s],
            population_configurations=pop_configs,
        )

    def test_samples(self):
        pop_configs = [
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
        ]
        samples = [
            msprime.Sample(population=0, time=0),
            msprime.Sample(population=1, time=1),
            msprime.Sample(population=2, time=2),
        ]
        sim = msprime.simulator_factory(
            samples=samples, population_configurations=pop_configs
        )
        self.assertEqual(sim.samples, samples)

    def test_new_old_style_model_changes_equal(self):
        models = [
            msprime.SweepGenicSelection(
                position=j, start_frequency=j, end_frequency=j, alpha=j, dt=j,
            )
            for j in range(1, 10)
        ]
        # Old style
        sim = msprime.simulator_factory(
            sample_size=2,
            Ne=10,
            demographic_events=[
                msprime.SimulationModelChange(None, model) for model in models
            ],
        )
        self.assertEqual(len(sim.model_change_events), len(models))
        for event, model in zip(sim.model_change_events, models):
            self.assertEqual(event.model, model)

        sim2 = msprime.simulator_factory(
            sample_size=2,
            Ne=10,
            model=[None]
            + [msprime.SimulationModelChange(None, model) for model in models],
        )
        self.assertEqual(sim.model_change_events, sim2.model_change_events)

    def test_model_change_old_style(self):
        main_model = msprime.SmcApproxCoalescent()
        sim = msprime.simulator_factory(
            Ne=100,
            sample_size=2,
            model=main_model,
            demographic_events=[
                msprime.SimulationModelChange(1, msprime.DiscreteTimeWrightFisher()),
                msprime.SimulationModelChange(2, None),
            ],
        )
        self.assertEqual(len(sim.model_change_events), 2)
        self.assertEqual(sim.model_change_events[0].time, 1)
        # When model=None we change to the standard coalescent
        self.assertEqual(sim.model_change_events[1].time, 2)
        self.assertEqual(sim.model_change_events[1].model.name, "hudson")

        # This should be the same in new notation
        sim = msprime.simulator_factory(
            Ne=100, sample_size=2, model=[main_model, (1, "dtwf"), (2, None)],
        )
        self.assertEqual(len(sim.model_change_events), 2)
        self.assertEqual(sim.model_change_events[0].time, 1)
        # When model=None we change to the standard coalescent
        self.assertEqual(sim.model_change_events[1].time, 2)
        self.assertEqual(sim.model_change_events[1].model.name, "hudson")

    def test_bad_sample_population_reference(self):
        # What happens when we reference a population that doesn't exist?
        with self.assertRaises(ValueError) as ve:
            msprime.simulate(
                samples=[
                    msprime.Sample(population=0, time=0),
                    msprime.Sample(population=1, time=0),
                ]
            )
        self.assertEqual(
            str(ve.exception), "Invalid population reference '1' in sample at index 1"
        )

        with self.assertRaises(ValueError) as ve:
            msprime.simulate(
                samples=[
                    msprime.Sample(population=0, time=0),
                    msprime.Sample(population=0, time=0),
                    msprime.Sample(population=-1, time=0),
                ]
            )
        self.assertEqual(
            str(ve.exception), "Negative population ID in sample at index 2"
        )


class TestSimulateInterface(unittest.TestCase):
    """
    Some simple test cases for the simulate() interface.
    """

    def test_defaults(self):
        n = 10
        ts = msprime.simulate(n)
        self.assertIsInstance(ts, tskit.TreeSequence)
        self.assertEqual(ts.get_sample_size(), n)
        self.assertEqual(ts.get_num_trees(), 1)
        self.assertEqual(ts.get_num_mutations(), 0)
        self.assertEqual(ts.get_sequence_length(), 1)
        self.assertEqual(len(list(ts.provenances())), 1)

    def test_positional_args_not_allowed(self):
        with self.assertRaises(TypeError):
            msprime.simulate(2, 100)

    def test_discrete_genome_recombination_map(self):
        # Cannot specify discrete_genome and recombination_map at once
        recomb_map = msprime.RecombinationMap.uniform_map(10, 0.1)
        with self.assertRaises(ValueError):
            msprime.simulate(10, discrete_genome=True, recombination_map=recomb_map)

    def test_discrete_genome_no_mutations(self):
        def run_sim(discrete_genome=None):
            return msprime.simulate(
                10,
                length=2,
                recombination_rate=1,
                discrete_genome=discrete_genome,
                random_seed=2134,
            )

        ts_discrete = run_sim(True)
        self.assertGreater(ts_discrete.num_trees, 1)
        self.assertTrue(has_discrete_genome(ts_discrete))

        ts_continuous = run_sim(False)
        self.assertGreater(ts_continuous.num_trees, 1)
        self.assertFalse(has_discrete_genome(ts_continuous))

        ts_default = run_sim()
        tables_default = ts_default.dump_tables()
        tables_continuous = ts_continuous.dump_tables()
        tables_continuous.provenances.clear()
        tables_default.provenances.clear()
        self.assertEqual(tables_default, tables_continuous)

    def test_discrete_genome_mutations(self):
        def run_sim(discrete_genome=None):
            return msprime.simulate(
                10,
                length=2,
                recombination_rate=1,
                mutation_rate=1,
                discrete_genome=discrete_genome,
                random_seed=2134,
            )

        ts_discrete = run_sim(True)
        self.assertGreater(ts_discrete.num_trees, 1)
        self.assertGreater(ts_discrete.num_sites, 1)
        self.assertTrue(has_discrete_genome(ts_discrete))

        ts_continuous = run_sim(False)
        self.assertGreater(ts_continuous.num_trees, 1)
        self.assertGreater(ts_discrete.num_sites, 1)
        self.assertFalse(has_discrete_genome(ts_continuous))

        ts_default = run_sim()
        tables_default = ts_default.dump_tables()
        tables_continuous = ts_continuous.dump_tables()
        tables_continuous.provenances.clear()
        tables_default.provenances.clear()
        self.assertEqual(tables_default, tables_continuous)

    def test_discrete_genome_migrations(self):
        def run_sim(discrete_genome=None):
            demography = msprime.Demography.stepping_stone_1d(2, 0.1)
            samples = demography.sample(5, 5)
            return msprime.simulate(
                samples=samples,
                demography=demography,
                length=5,
                recombination_rate=1,
                discrete_genome=discrete_genome,
                record_migrations=True,
                random_seed=2134,
            )

        ts_discrete = run_sim(True)
        self.assertGreater(ts_discrete.num_trees, 1)
        self.assertGreater(ts_discrete.num_migrations, 1)
        self.assertTrue(has_discrete_genome(ts_discrete))

        ts_continuous = run_sim(False)
        self.assertGreater(ts_continuous.num_trees, 1)
        self.assertGreater(ts_continuous.num_migrations, 1)
        self.assertFalse(has_discrete_genome(ts_continuous))

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

    def test_provenance(self):
        ts = msprime.simulate(10)
        self.assertEqual(ts.num_provenances, 1)
        self.verify_provenance(ts.provenance(0))
        for ts in msprime.simulate(10, num_replicates=10):
            self.assertEqual(ts.num_provenances, 1)
            self.verify_provenance(ts.provenance(0))

    def test_end_time(self):
        ts = msprime.simulate(15, recombination_rate=2, random_seed=8, end_time=0.1)
        for tree in ts.trees():
            for root in tree.roots:
                self.assertEqual(tree.time(root), 0.1)

    def test_replicates(self):
        n = 20
        num_replicates = 10
        count = 0
        for ts in msprime.simulate(n, num_replicates=num_replicates):
            count += 1
            self.assertIsInstance(ts, tskit.TreeSequence)
            self.assertEqual(ts.get_sample_size(), n)
            self.assertEqual(ts.get_num_trees(), 1)
        self.assertEqual(num_replicates, count)

    def test_mutations(self):
        n = 10
        ts = msprime.simulate(n, mutation_rate=10)
        self.assertIsInstance(ts, tskit.TreeSequence)
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
        mutgen = msprime.mutations._simple_mutation_generator(
            1, 1, _msprime.RandomGenerator(1)
        )
        with self.assertRaises(ValueError):
            msprime.simulate(n, mutation_generator=mutgen)

    def test_mutation_interface(self):
        for bad_type in [{}, self]:
            self.assertRaises(TypeError, msprime.simulate, 10, mutation_rate=bad_type)
        for bad_value in ["x", [], [[], []]]:
            self.assertRaises(ValueError, msprime.simulate, 10, mutation_rate=bad_value)

    def test_recombination(self):
        n = 10
        ts = msprime.simulate(n, recombination_rate=10)
        self.assertIsInstance(ts, tskit.TreeSequence)
        self.assertEqual(ts.sample_size, n)
        self.assertGreater(ts.num_trees, 1)
        self.assertEqual(ts.num_mutations, 0)

    def test_gene_conversion_simple_map(self):
        n = 10
        ts = msprime.simulate(
            n,
            gene_conversion_rate=1,
            gene_conversion_track_length=1,
            recombination_map=msprime.RecombinationMap.uniform_map(
                10, 1, discrete=True
            ),
        )
        self.assertIsInstance(ts, tskit.TreeSequence)
        self.assertEqual(ts.num_samples, n)
        self.assertGreater(ts.num_trees, 1)

    def test_gene_conversion_continuous(self):
        rm = msprime.RecombinationMap.uniform_map(10, 1, discrete=False)
        with self.assertRaises(ValueError):
            msprime.simulate(
                10,
                gene_conversion_rate=1,
                gene_conversion_track_length=1,
                recombination_map=rm,
            )

    @unittest.skip("Cannot use GC with default recomb map")
    def test_gene_conversion_default_map(self):
        n = 10
        # FIXME we have to be quite delicate with the GC code at the moment.
        # If we take the default where we have a very large number of loci,
        # we might be getting overflows. It's not clear what happening in any case.
        ts = msprime.simulate(n, gene_conversion_rate=1, gene_conversion_track_length=1)
        self.assertIsInstance(ts, tskit.TreeSequence)
        self.assertEqual(ts.num_samples, n)
        self.assertGreater(ts.num_trees, 1)

    def test_num_labels(self):
        # Running simulations with different numbers of labels in the default
        # setting should have no effect.
        tables = [
            msprime.simulate(10, num_labels=num_labels, random_seed=1).tables
            for num_labels in range(1, 5)
        ]
        for t in tables:
            t.provenances.clear()
        for t in tables:
            self.assertEqual(t, tables[0])

    def test_replicate_index(self):
        tables_1 = list(msprime.simulate(10, num_replicates=5, random_seed=1))[4].tables
        tables_2 = msprime.simulate(10, replicate_index=4, random_seed=1).tables
        tables_1.provenances.clear()
        tables_2.provenances.clear()
        self.assertEqual(tables_1, tables_2)

        with self.assertRaises(ValueError) as cm:
            msprime.simulate(5, replicate_index=5)
        self.assertEqual(
            "Cannot specify replicate_index without random_seed as this "
            "has the same effect as not specifying replicate_index i.e. a "
            "random tree sequence",
            str(cm.exception),
        )
        with self.assertRaises(ValueError) as cm:
            msprime.simulate(5, random_seed=1, replicate_index=5, num_replicates=26)
        self.assertEqual(
            "Cannot specify replicate_index with num_replicates as only "
            "the replicate_index specified will be returned.",
            str(cm.exception),
        )


class TestRecombinationMap(unittest.TestCase):
    """
    Tests for the RecombinationMap class.
    """

    # TODO these are incomplete.
    def test_discrete(self):
        for truthy in [True, False, {}, None, "ser"]:
            rm = msprime.RecombinationMap.uniform_map(1, 0, discrete=truthy)
            self.assertEqual(rm.discrete, bool(truthy))

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

        # Test mean_recombination_rate is correct after reading from
        # a hapmap file. RecombinationMap.read_hapmap() ignores the cM
        # field, so here we test against using the cM field directly.
        def hapmap_rr(hapmap_file):
            first_pos = 0
            with open(hapmap_file) as f:
                next(f)  # skip header
                for line in f:
                    pos, rate, cM = map(float, line.split()[1:4])
                    if cM == 0:
                        first_pos = pos
            return cM / 100 / (pos - first_pos)

        hapmap = """chr pos        rate                    cM
                    1   4283592    3.79115663174456        0
                    1   4361401    0.0664276817058413      0.294986106359414
                    1   7979763   10.9082897515584         0.535345505591925
                    1   8007051    0.0976780648822495      0.833010916332456
                    1   8762788    0.0899929572085616      0.906829844052373
                    1   9477943    0.0864382908650907      0.971188757364862
                    1   9696341    4.76495005895746        0.990066707213216
                    1   9752154    0.0864316558730679      1.25601286485381
                    1   9881751    0.0                     1.26721414815999"""
        with tempfile.TemporaryDirectory() as temp_dir:
            hapfile = os.path.join(temp_dir, "hapmap.txt")
            with open(hapfile, "w") as f:
                f.write(hapmap)
            recomb_map = msprime.RecombinationMap.read_hapmap(f.name)
            mean_rr = recomb_map.mean_recombination_rate
            mean_rr2 = hapmap_rr(hapfile)
        self.assertAlmostEqual(mean_rr, mean_rr2, places=15)


class TestReprRoundTrip(unittest.TestCase):
    """
    Tests that we can eval the repr of objects to round trip them.
    """

    def assert_repr_round_trip(self, obj_list):
        for obj in obj_list:
            obj_copy = eval(repr(obj), globals(), msprime.__dict__)
            self.assertEqual(obj_copy, obj)
            self.assertFalse(obj_copy is obj)

    def test_population(self):
        examples = [
            msprime.Population(),
            msprime.Population(initial_size=2),
            msprime.Population(growth_rate=5),
            msprime.Population(initial_size=234, growth_rate=10),
        ]
        self.assert_repr_round_trip(examples)

    def test_population_parameters_change(self):
        examples = [
            msprime.PopulationParametersChange(time=1, initial_size=1),
            msprime.PopulationParametersChange(time=1, growth_rate=2),
            msprime.PopulationParametersChange(time=1, growth_rate=1, population=2),
            msprime.PopulationParametersChange(
                time=3, initial_size=3, growth_rate=1, population=2
            ),
        ]
        self.assert_repr_round_trip(examples)

    def test_migration_rate_change(self):
        examples = [
            msprime.MigrationRateChange(time=1, rate=1),
            msprime.MigrationRateChange(time=1, rate=1, source=1, dest=2),
        ]
        self.assert_repr_round_trip(examples)

    def test_mass_migration(self):
        examples = [
            msprime.MassMigration(time=1, source=1, dest=2),
            msprime.MassMigration(time=1, source=1, dest=2, proportion=0.2),
        ]
        self.assert_repr_round_trip(examples)

    def test_simulation_model_change(self):
        examples = [
            msprime.SimulationModelChange(),
            msprime.SimulationModelChange(model="hudson"),
            msprime.SimulationModelChange(model=msprime.DiscreteTimeWrightFisher()),
            msprime.SimulationModelChange(
                model=msprime.BetaCoalescent(alpha=1, truncation_point=2)
            ),
        ]
        self.assert_repr_round_trip(examples)

    def test_simple_bottleneck(self):
        examples = [
            msprime.SimpleBottleneck(time=10, population=2),
            msprime.SimpleBottleneck(time=10, population=2, proportion=0.5),
        ]
        self.assert_repr_round_trip(examples)

    def test_instantaneous_bottleneck(self):
        examples = [
            msprime.InstantaneousBottleneck(time=10, population=1),
            msprime.InstantaneousBottleneck(time=10, population=1, strength=10),
        ]
        self.assert_repr_round_trip(examples)

    def test_census_event(self):
        examples = [
            msprime.CensusEvent(time=10),
        ]
        self.assert_repr_round_trip(examples)

    def test_simulation_models(self):
        examples = [
            msprime.StandardCoalescent(),
            msprime.SmcApproxCoalescent(),
            msprime.SmcPrimeApproxCoalescent(),
            msprime.DiscreteTimeWrightFisher(),
            msprime.WrightFisherPedigree(),
            msprime.BetaCoalescent(),
            msprime.BetaCoalescent(alpha=1, truncation_point=10),
            msprime.DiracCoalescent(),
            msprime.DiracCoalescent(psi=1234, c=56),
            msprime.SweepGenicSelection(
                position=1, start_frequency=0.5, end_frequency=0.9, alpha=1, dt=1e-4,
            ),
        ]
        self.assert_repr_round_trip(examples)
