#
# Copyright (C) 2016 University of Oxford
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
Test cases for demographic events in msprime.
"""
from __future__ import print_function
from __future__ import division

import itertools
import math
import tempfile
import unittest

import numpy as np

import msprime


class TestBadDemographicParameters(unittest.TestCase):
    """
    Tests for nonsensical demographic parameters.
    """

    def test_bad_population_size(self):
        for bad_size in [-1, 0, -1e300]:
            self.assertRaises(
                ValueError, msprime.PopulationParametersChange, 0, initial_size=bad_size)
            self.assertRaises(
                ValueError, msprime.PopulationConfiguration, initial_size=bad_size)

    def test_bad_sample_size(self):
        for bad_size in [-1, -1e300]:
            self.assertRaises(
                ValueError, msprime.PopulationConfiguration, sample_size=bad_size)


class TestDeprecatedParameters(unittest.TestCase):
    """
    Tests to check that aliased parameters are handled correctly.
    """
    def test_mass_migration_dest(self):
        self.assertRaises(
            ValueError, msprime.MassMigration, time=0, source=0, dest=0, destination=0)
        for j in range(10):
            e = msprime.MassMigration(time=0.1, source=0, dest=j, proportion=0.5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.source, 0)
            self.assertEqual(e.dest, j)
            self.assertEqual(e.proportion, 0.5)
            e = msprime.MassMigration(0.1, 0, j, 0.5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.source, 0)
            self.assertEqual(e.dest, j)
            self.assertEqual(e.proportion, 0.5)
            e = msprime.MassMigration(time=0.1, source=0, destination=j, proportion=0.5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.source, 0)
            self.assertEqual(e.dest, j)
            self.assertEqual(e.proportion, 0.5)

    def test_population_parameters_population_id(self):
        self.assertRaises(
            ValueError, msprime.PopulationParametersChange,
            time=0, initial_size=1.1, growth_rate=0.1, population=1,
            population_id=1)

        for j in range(10):
            e = msprime.PopulationParametersChange(
                time=0.1, initial_size=1.1, growth_rate=0.1, population=j)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.initial_size, 1.1)
            self.assertEqual(e.growth_rate, 0.1)
            self.assertEqual(e.population, j)
            e = msprime.PopulationParametersChange(0.1, 1.1, 0.1, j)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.initial_size, 1.1)
            self.assertEqual(e.growth_rate, 0.1)
            self.assertEqual(e.population, j)
            e = msprime.PopulationParametersChange(
                time=0.1, initial_size=1.1, growth_rate=0.1, population_id=j)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.initial_size, 1.1)
            self.assertEqual(e.growth_rate, 0.1)
            self.assertEqual(e.population, j)

    def test_simple_bottleneck_population_id(self):
        self.assertRaises(
            ValueError, msprime.SimpleBottleneck,
            time=0, population=1, proportion=5, population_id=1)

        for j in range(10):
            e = msprime.SimpleBottleneck(time=0.1, population=j, proportion=5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.population, j)
            self.assertEqual(e.proportion, 5)
            e = msprime.SimpleBottleneck(0.1, j, 5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.population, j)
            self.assertEqual(e.proportion, 5)
            e = msprime.SimpleBottleneck(time=0.1, population_id=j, proportion=5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.population, j)
            self.assertEqual(e.proportion, 5)

    def test_instantaneous_bottleneck_population_id(self):
        self.assertRaises(
            ValueError, msprime.InstantaneousBottleneck,
            time=0, population=1, strength=1, population_id=1)

        for j in range(10):
            e = msprime.InstantaneousBottleneck(time=0.1, population=j, strength=5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.population, j)
            self.assertEqual(e.strength, 5)
            e = msprime.InstantaneousBottleneck(0.1, j, 5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.population, j)
            self.assertEqual(e.strength, 5)
            e = msprime.InstantaneousBottleneck(time=0.1, population_id=j, strength=5)
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.population, j)
            self.assertEqual(e.strength, 5)


class TestRateConversions(unittest.TestCase):
    """
    Tests for the demographic events interface.
    """
    def test_size_change(self):
        g = 100
        new_size = 512
        event = msprime.PopulationParametersChange(time=g, initial_size=new_size)
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": -1,
            "initial_size": new_size
        }
        self.assertEqual(event.get_ll_representation(1), ll_event)

    def test_growth_rate_change(self):
        g = 512
        growth_rate = 1
        event = msprime.PopulationParametersChange(
            time=g, growth_rate=growth_rate, population=1)
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": 1,
            "growth_rate": growth_rate
        }
        self.assertEqual(event.get_ll_representation(1), ll_event)

    def test_growth_rate_and_size_change(self):
        g = 1024
        growth_rate = 2
        initial_size = 8192
        event = msprime.PopulationParametersChange(
            time=g, initial_size=initial_size,
            growth_rate=growth_rate, population=1)
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": 1,
            "initial_size": initial_size,
            "growth_rate": growth_rate
        }
        self.assertEqual(event.get_ll_representation(1), ll_event)

    def test_migration_rate_change(self):
        g = 1024
        migration_rate = 0.125
        d = 2
        event = msprime.MigrationRateChange(time=g, rate=migration_rate)
        ll_event = {
            "type": "migration_rate_change",
            "time": g,
            "matrix_index": -1,
            "migration_rate": migration_rate
        }
        self.assertEqual(event.get_ll_representation(d), ll_event)


class TestDemographyDebuggerOutput(unittest.TestCase):
    """
    Tests for the demography debug interface.
    """
    def verify_debug(
            self, population_configurations, migration_matrix, demographic_events):

        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        # Check the reprs
        s = repr(dd.epochs)
        self.assertGreater(len(s), 0)
        with tempfile.TemporaryFile("w+") as f:
            dd.print_history(f)
            f.seek(0)
            debug_output = f.read()
        # TODO when there is better output, write some tests to
        # verify its format.
        self.assertGreater(len(debug_output), 0)

    def test_zero_samples(self):
        population_configurations = [
            msprime.PopulationConfiguration(0)]
        self.verify_debug(population_configurations, [[0]], [])
        self.verify_debug(None, None, [])

    def test_one_population(self):
        population_configurations = [
            msprime.PopulationConfiguration(10)]
        migration_matrix = [[0]]
        demographic_events = [
            msprime.PopulationParametersChange(0.1, initial_size=2),
            msprime.PopulationParametersChange(0.1, growth_rate=10)]
        self.verify_debug(
            population_configurations, migration_matrix, demographic_events)

    def test_no_events(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10)]
        migration_matrix = [[0, 0], [0, 0]]
        self.verify_debug(population_configurations, migration_matrix, [])

    def test_demographic_events(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10)]
        migration_matrix = [[0, 0], [0, 0]]
        demographic_events = [
            msprime.PopulationParametersChange(0.1, initial_size=2),
            msprime.PopulationParametersChange(0.1, growth_rate=10),
            msprime.MassMigration(0.2, source=1, dest=0),
            msprime.MigrationRateChange(0.2, rate=0),
            msprime.MigrationRateChange(0.4, matrix_index=(0, 1), rate=1),
            msprime.MigrationRateChange(0.4, matrix_index=(1, 0), rate=1),
            msprime.InstantaneousBottleneck(0.5, population=0, strength=100)]
        self.verify_debug(
            population_configurations, migration_matrix, demographic_events)


class TestDemographyDebugger(unittest.TestCase):
    """
    Tests for the demography debugger. Ensure that we compute the correct
    population sizes etc.
    """

    def test_equal_after(self):
        population_configurations = [
            msprime.PopulationConfiguration(sample_size=10),
            msprime.PopulationConfiguration(sample_size=20)]
        msprime.DemographyDebugger(
            population_configurations=population_configurations)
        self.assertEqual(population_configurations[0].sample_size, 10)
        self.assertEqual(population_configurations[1].sample_size, 20)

    def test_one_pop_zero_events(self):
        dd = msprime.DemographyDebugger(
            population_configurations=[msprime.PopulationConfiguration()])
        self.assertEqual(len(dd.epochs), 1)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertTrue(math.isinf(e.end_time))
        self.assertEqual(len(e.demographic_events), 0)
        self.assertEqual(len(e.populations), 1)
        self.assertEqual(e.migration_matrix, [[0]])
        pop = e.populations[0]
        self.assertEqual(pop.growth_rate, 0)
        self.assertEqual(pop.start_size, 1)
        self.assertEqual(pop.end_size, 1)

    def test_two_pop_different_sizes(self):
        dd = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=10),
                msprime.PopulationConfiguration(initial_size=20)])
        self.assertEqual(len(dd.epochs), 1)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertTrue(math.isinf(e.end_time))
        self.assertEqual(len(e.demographic_events), 0)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 0], [0, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
        self.assertEqual(e.populations[0].start_size, 10)
        self.assertEqual(e.populations[1].start_size, 20)

    def test_two_pop_different_growth_rates(self):
        g1 = 0.1
        g2 = 0.5
        p0_end_size = 10 * math.exp(-g1 * 10)
        p1_end_size = 20 * math.exp(-g2 * 10)
        dd = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=10, growth_rate=g1),
                msprime.PopulationConfiguration(initial_size=20, growth_rate=g2)],
            demographic_events=[
                msprime.PopulationParametersChange(time=10, growth_rate=0)])
        # Make sure we're testing the __repr__ paths.
        s = repr(dd)
        self.assertGreater(len(s), 0)
        self.assertEqual(len(dd.epochs), 2)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertEqual(e.end_time, 10)
        self.assertEqual(len(e.demographic_events), 0)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 0], [0, 0]])
        self.assertEqual(e.populations[0].start_size, 10)
        self.assertEqual(e.populations[0].end_size, p0_end_size)
        self.assertEqual(e.populations[1].start_size, 20)
        self.assertEqual(e.populations[1].end_size, p1_end_size)

        e = dd.epochs[1]
        self.assertEqual(e.start_time, 10)
        self.assertTrue(math.isinf(e.end_time))
        self.assertEqual(len(e.demographic_events), 1)
        d = e.demographic_events[0]
        self.assertEqual(d.time, 10)
        self.assertEqual(d.growth_rate, 0)
        self.assertEqual(d.initial_size, None)
        self.assertEqual(d.population, -1)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 0], [0, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
        self.assertEqual(e.populations[0].start_size, p0_end_size)
        self.assertEqual(e.populations[1].start_size, p1_end_size)

    def test_two_pop_update_migration_rate(self):
        dd = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=100),
                msprime.PopulationConfiguration(initial_size=100)],
            demographic_events=[
                msprime.MigrationRateChange(time=20, rate=1),
                msprime.MigrationRateChange(time=22, rate=1.7, matrix_index=(0, 1))],
            migration_matrix=[[0, 0.25], [0, 0]])
        self.assertEqual(len(dd.epochs), 3)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertEqual(e.end_time, 20)
        self.assertEqual(e.migration_matrix, [[0, 0.25], [0, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
            self.assertEqual(pop.start_size, 100)

        e = dd.epochs[1]
        self.assertEqual(e.start_time, 20)
        self.assertTrue(e.end_time, 22)
        self.assertEqual(len(e.demographic_events), 1)
        d = e.demographic_events[0]
        self.assertEqual(d.matrix_index, None)
        self.assertEqual(d.time, 20)
        self.assertEqual(d.rate, 1)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 1], [1, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
            self.assertEqual(pop.start_size, 100)

        e = dd.epochs[2]
        self.assertEqual(e.start_time, 22)
        self.assertTrue(math.isinf(e.end_time))
        self.assertEqual(len(e.demographic_events), 1)
        d = e.demographic_events[0]
        self.assertEqual(d.matrix_index, (0, 1))
        self.assertEqual(d.time, 22)
        self.assertEqual(d.rate, 1.7)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 1.7], [1, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
            self.assertEqual(pop.start_size, 100)

    def test_two_pop_change_growth_rates(self):
        alpha = 0.33
        N0 = 1000
        N1 = 10
        t1 = 5
        t2 = 10
        t3 = 15
        dd = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=N0, growth_rate=alpha),
                msprime.PopulationConfiguration(initial_size=N1, growth_rate=0)],
            demographic_events=[
                # p1 changes growth rate to alpha
                msprime.PopulationParametersChange(
                    population=1, time=t1, growth_rate=alpha),
                # p0 changes growth rate to -alpha
                msprime.PopulationParametersChange(
                    population=0, time=t2, growth_rate=-alpha),
                # Both change growth_rate to 0 at t3.
                msprime.PopulationParametersChange(time=t3, growth_rate=0)])

        self.assertEqual(len(dd.epochs), 4)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertEqual(e.end_time, t1)
        self.assertEqual(len(e.demographic_events), 0)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 0], [0, 0]])
        self.assertEqual(e.populations[0].start_size, N0)
        n0 = N0 * math.exp(-alpha * t1)
        self.assertEqual(e.populations[0].end_size, n0)
        self.assertEqual(e.populations[1].start_size, N1)
        self.assertEqual(e.populations[1].end_size, N1)

        e = dd.epochs[1]
        self.assertEqual(e.start_time, t1)
        self.assertEqual(e.end_time, t2)
        self.assertEqual(len(e.demographic_events), 1)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 0], [0, 0]])
        self.assertEqual(e.populations[0].start_size, n0)
        n0 = N0 * math.exp(-alpha * t2)
        self.assertEqual(e.populations[0].end_size, n0)
        self.assertEqual(e.populations[1].start_size, N1)
        n1 = N1 * math.exp(-alpha * (t2 - t1))
        self.assertEqual(e.populations[1].end_size, n1)

        e = dd.epochs[2]
        self.assertEqual(e.start_time, t2)
        self.assertEqual(e.end_time, t3)
        self.assertEqual(len(e.demographic_events), 1)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 0], [0, 0]])
        self.assertEqual(e.populations[0].start_size, n0)
        n0 = n0 * math.exp(alpha * (t3 - t2))
        self.assertEqual(e.populations[0].end_size, n0)
        self.assertEqual(e.populations[1].start_size, n1)
        n1 = N1 * math.exp(-alpha * (t3 - t1))
        self.assertEqual(e.populations[1].end_size, n1)

        e = dd.epochs[3]
        self.assertEqual(e.start_time, t3)
        self.assertTrue(math.isinf(e.end_time))
        self.assertEqual(len(e.demographic_events), 1)
        self.assertEqual(len(e.populations), 2)
        self.assertEqual(e.migration_matrix, [[0, 0], [0, 0]])
        self.assertEqual(e.populations[0].start_size, n0)
        self.assertEqual(e.populations[0].end_size, n0)
        self.assertEqual(e.populations[1].start_size, n1)
        self.assertEqual(e.populations[1].end_size, n1)


class TestCoalescenceLocations(unittest.TestCase):
    """
    Tests that coalescences happen in demes that they are supposed to
    for simple models.
    """
    def test_two_pops_single_sample(self):
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(0),
        ]
        t = 5
        demographic_events = [
            msprime.MassMigration(time=t, source=0, dest=2),
            msprime.MassMigration(time=t, source=1, dest=2),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        self.assertEqual(tree.root, 2)
        self.assertGreater(tree.time(2), t)
        self.assertEqual(tree.population(0), 0)
        self.assertEqual(tree.population(1), 1)
        self.assertEqual(tree.population(2), 2)
        self.assertEqual(ts.node(0).population, 0)
        self.assertEqual(ts.node(1).population, 1)
        self.assertEqual(list(ts.samples()), [0, 1])
        self.assertEqual(list(ts.samples(0)), [0])
        self.assertEqual(list(ts.samples(1)), [1])

    def test_two_pops_multiple_samples(self):
        # Made absolutely sure that all samples have coalesced within
        # the source deme
        n = 10
        t = 100
        population_configurations = [
            msprime.PopulationConfiguration(n // 2),
            msprime.PopulationConfiguration(n // 2),
            msprime.PopulationConfiguration(0),
        ]
        demographic_events = [
            msprime.MassMigration(time=t, source=0, dest=2),
            msprime.MassMigration(time=t, source=1, dest=2),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        self.assertEqual(tree.root, 2 * n - 2)
        self.assertGreater(tree.time(tree.root), t)
        for j in range(n // 2):
            self.assertEqual(tree.population(j), 0)
            self.assertEqual(tree.population(n // 2 + j), 1)
            self.assertEqual(ts.get_population(j), 0)
            self.assertEqual(ts.get_population(n // 2 + j), 1)
        self.assertEqual(tree.population(tree.root), 2)

        self.assertTrue(np.array_equal(
            ts.samples(0), np.arange(n // 2, dtype=np.int32)))
        self.assertTrue(np.array_equal(
            ts.samples(1), np.arange(n // 2, n, dtype=np.int32)))
        self.assertTrue(np.array_equal(
            ts.samples(2), np.array([], dtype=np.int32)))

        # self.assertEqual(ts.samples(0), list(range(n // 2)))
        # self.assertEqual(ts.samples(1), list(range(n // 2, n)))
        # self.assertEqual(ts.samples(2), [])

    def test_three_pops_migration(self):
        n = 9
        t = 100
        population_configurations = [
            msprime.PopulationConfiguration(n // 3),
            msprime.PopulationConfiguration(n // 3),
            msprime.PopulationConfiguration(n // 3),
        ]
        # Start migrating everyone into 0 after t
        demographic_events = [
            msprime.MigrationRateChange(time=t, matrix_index=(1, 0), rate=1),
            msprime.MigrationRateChange(time=t, matrix_index=(2, 0), rate=1),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        self.assertEqual(tree.root, 2 * n - 2)
        self.assertGreater(tree.time(tree.root), t)
        for j in range(n // 3):
            self.assertEqual(tree.population(j), 0)
            self.assertEqual(tree.population(n // 3 + j), 1)
            self.assertEqual(tree.population(2 * (n // 3) + j), 2)
            self.assertEqual(ts.get_population(j), 0)
            self.assertEqual(ts.get_population(n // 3 + j), 1)
            self.assertEqual(ts.get_population(2 * (n // 3) + j), 2)
        # The MRCAs of 0, 1 and 3 must have occured in deme 0
        self.assertEqual(tree.population(tree.get_mrca(0, n // 3)), 0)
        self.assertEqual(
            tree.population(tree.get_mrca(0, 2 * (n // 3))), 0)
        # The MRCAs of all the samples within each deme must have
        # occured within that deme
        for k in range(3):
            deme_samples = range(k * (n // 3), (k + 1) * (n // 3))
            for u, v in itertools.combinations(deme_samples, 2):
                mrca_pop = tree.population(tree.get_mrca(u, v))
                self.assertEqual(k, mrca_pop)
        self.assertTrue(np.array_equal(
            ts.samples(0), np.arange(n // 3, dtype=np.int32)))
        self.assertTrue(np.array_equal(
            ts.samples(1), np.arange(n // 3, 2 * (n // 3), dtype=np.int32)))
        self.assertTrue(np.array_equal(
            ts.samples(2), np.arange(2 * (n // 3), n, dtype=np.int32)))

    def test_four_pops_three_mass_migrations(self):
        t1 = 1
        t2 = 100
        t3 = 200
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
        ]
        # We migrate the lineages to the next step by step.
        demographic_events = [
            msprime.MassMigration(time=t1, source=0, dest=1),
            msprime.MassMigration(time=t2, source=1, dest=2),
            msprime.MassMigration(time=t3, source=2, dest=3),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        # Check the leaves have the correct population.
        for j in range(4):
            self.assertEqual(tree.population(j), j)
            self.assertEqual(ts.get_population(j), j)
            self.assertEqual(ts.samples(j), [j])
        # The MRCA of 0 and 1 should happen in 1 at time > t1, and < t2
        u = tree.get_mrca(0, 1)
        self.assertEqual(u, 4)
        self.assertEqual(tree.population(u), 1)
        g = tree.time(u) * 4
        self.assertTrue(t1 < g < t2)
        # The MRCA of 0, 1 and 2 should happen in 2 at time > t2 and < t3
        u = tree.get_mrca(0, 2)
        self.assertEqual(u, 5)
        self.assertEqual(tree.population(u), 2)
        self.assertTrue(t2 < tree.time(u) < t3)
        # The MRCA of 0, 1, 2 and 3 should happen in 3 at time > t3
        u = tree.get_mrca(0, 3)
        self.assertEqual(u, 6)
        self.assertEqual(tree.population(u), 3)
        self.assertGreater(tree.time(u), t3)

    def test_empty_demes(self):
        t1 = 1
        t2 = 100
        t3 = 200
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(0),
            msprime.PopulationConfiguration(0),
            msprime.PopulationConfiguration(1),
        ]
        # We migrate the lineages to the next step by step.
        demographic_events = [
            msprime.MassMigration(time=t1, source=0, dest=1),
            msprime.MassMigration(time=t2, source=1, dest=2),
            msprime.MassMigration(time=t3, source=2, dest=3),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        # Check the leaves have the correct population.
        self.assertEqual(tree.population(0), 0)
        self.assertEqual(tree.population(1), 3)
        self.assertEqual(ts.node(0).population, 0)
        self.assertEqual(ts.node(1).population, 3)
        self.assertEqual(list(ts.samples(0)), [0])
        self.assertEqual(list(ts.samples(1)), [])
        self.assertEqual(list(ts.samples(2)), [])
        self.assertEqual(list(ts.samples(3)), [1])
        # The MRCA of 0, 1 in 3 at time > t3
        u = tree.get_mrca(0, 1)
        self.assertEqual(u, 2)
        self.assertEqual(tree.population(u), 3)
        g = tree.time(u) * 4
        self.assertGreater(g, t3)

    def test_migration_rate_directionality(self):
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(0),
        ]
        t = 5
        demographic_events = [
            msprime.MigrationRateChange(time=t, rate=1, matrix_index=(0, 2)),
            msprime.MigrationRateChange(time=t, rate=1, matrix_index=(1, 2)),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        self.assertEqual(tree.root, 2)
        self.assertGreater(tree.time(2), t / 4)
        self.assertEqual(tree.population(0), 0)
        self.assertEqual(tree.population(1), 1)
        self.assertEqual(tree.population(2), 2)
        self.assertEqual(ts.node(0).population, 0)
        self.assertEqual(ts.node(1).population, 1)
        self.assertEqual(list(ts.samples()), [0, 1])
        self.assertEqual(list(ts.samples(0)), [0])
        self.assertEqual(list(ts.samples(1)), [1])

    def test_many_demes(self):
        num_demes = 300
        population_configurations = [
            msprime.PopulationConfiguration(1)] + [
            msprime.PopulationConfiguration(0) for _ in range(num_demes - 2)
            ] + [msprime.PopulationConfiguration(1)]
        t = 5
        demographic_events = [
            msprime.MassMigration(time=t, source=0, dest=num_demes - 1),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        self.assertEqual(tree.root, 2)
        self.assertGreater(tree.time(2), t)
        self.assertEqual(tree.population(0), 0)
        self.assertEqual(tree.population(1), num_demes - 1)
        self.assertEqual(tree.population(2), num_demes - 1)
        self.assertEqual(ts.node(0).population, 0)
        self.assertEqual(ts.node(1).population, num_demes - 1)

    def test_instantaneous_bottleneck_locations(self):
        population_configurations = [
            msprime.PopulationConfiguration(100),
            msprime.PopulationConfiguration(100),
            msprime.PopulationConfiguration(100),
        ]
        strength = 1e6  # Make sure everyone coalescences.
        t1 = 0.0001
        t2 = 0.0002
        t3 = 0.0003
        t4 = 0.0004
        demographic_events = [
            msprime.InstantaneousBottleneck(time=t1, population=0, strength=strength),
            msprime.InstantaneousBottleneck(time=t2, population=1, strength=strength),
            msprime.InstantaneousBottleneck(time=t3, population=2, strength=strength),
            msprime.MassMigration(time=t4, source=2, dest=0),
            msprime.MassMigration(time=t4, source=1, dest=0)
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        tree = next(ts.trees())
        self.assertGreater(tree.time(tree.root), t4)
        self.assertEqual(tree.population(tree.root), 0)
        # The parent of all the samples from each deme should be in that deme.
        for pop in range(3):
            parents = [
                tree.get_parent(u) for u in ts.samples(population=pop)]
            for v in parents:
                self.assertEqual(tree.population(v), pop)


class TestMigrationRecords(unittest.TestCase):
    """
    Tests that migrations happen where they should for simple models.
    """
    def verify_migrations(self, ts):
        """
        Verifies that the migrations for the specified tree sequence
        have the required properties.
        """
        migrations = list(ts.migrations())
        self.assertEqual(ts.num_migrations, len(migrations))
        oldest_t = max(node.time for node in ts.nodes())
        for mig in migrations:
            self.assertGreaterEqual(mig.left, 0)
            self.assertLessEqual(mig.right, ts.sequence_length)
            self.assertNotEqual(mig.source, mig.dest)
            self.assertTrue(0 <= mig.node < ts.num_nodes)
            self.assertTrue(0 <= mig.time < oldest_t)
        # Migrations must be listed in non-decreasing time order.
        for j in range(1, len(migrations)):
            self.assertTrue(migrations[j - 1].time <= migrations[j].time)

    def test_two_pops_single_sample(self):
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(0),
        ]
        t = 5
        demographic_events = [
            msprime.MassMigration(time=t, source=0, dest=2),
            msprime.MassMigration(time=t, source=1, dest=2),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1, record_migrations=True)
        self.verify_migrations(ts)
        migrations = list(ts.migrations())
        self.assertEqual(len(migrations), 2)
        m0 = migrations[0]
        m1 = migrations[1]
        self.assertEqual(m0.left, 0)
        self.assertEqual(m0.right, 1)
        self.assertEqual(m0.source, 0)
        self.assertEqual(m0.dest, 2)
        self.assertEqual(m0.time, t)
        self.assertEqual(m1.left, 0)
        self.assertEqual(m1.right, 1)
        self.assertEqual(m1.source, 1)
        self.assertEqual(m1.dest, 2)
        self.assertEqual(m1.time, t)

    def test_two_pops_asymmetric_migrations(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            migration_matrix=[
                [0, 0],
                [1, 0]],
            # Can migrate from 1 to 0 but not vice-versa
            random_seed=1, record_migrations=True)
        self.verify_migrations(ts)
        migrations = list(ts.migrations())
        self.assertGreater(len(migrations), 0)
        for mig in migrations:
            self.assertGreater(mig.time, 0)
            self.assertEqual(ts.node(mig.node).population, 1)
            self.assertEqual(mig.source, 1)
            self.assertEqual(mig.dest, 0)
            self.assertEqual(mig.left, 0)
            self.assertEqual(mig.right, 1)

    def test_two_pops_asymmetric_migrations_recombination(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        ts = msprime.simulate(
            recombination_rate=1,
            population_configurations=population_configurations,
            migration_matrix=[
                [0, 0],
                [1, 0]],
            # Can migrate from 1 to 0 but not vice-versa
            random_seed=1, record_migrations=True)
        self.verify_migrations(ts)
        self.assertGreater(ts.num_trees, 1)
        migrations = list(ts.migrations())
        self.assertGreater(len(migrations), 0)
        for mig in migrations:
            self.assertGreater(mig.time, 0)
            self.assertEqual(ts.node(mig.node).population, 1)
            self.assertEqual(mig.source, 1)
            self.assertEqual(mig.dest, 0)
            self.assertGreaterEqual(mig.left, 0)
            self.assertLessEqual(mig.right, 1)

    def test_two_pops_mass_migration_recombination(self):
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
        ]
        ts = msprime.simulate(
            recombination_rate=10,
            population_configurations=population_configurations,
            demographic_events=[
                msprime.MassMigration(time=20, source=1, dest=0, proportion=1)
            ],
            random_seed=1, record_migrations=True)
        self.verify_migrations(ts)
        self.assertGreater(ts.num_trees, 10)
        migrations = list(ts.migrations())
        self.assertGreater(len(migrations), 0)
        for mig in migrations:
            self.assertGreater(mig.time, 0)
            self.assertEqual(mig.node, 1)
            self.assertEqual(mig.source, 1)
            self.assertEqual(mig.dest, 0)
            self.assertGreaterEqual(mig.left, 0)
            self.assertLessEqual(mig.right, 1)


class TestTimeUnits(unittest.TestCase):
    """
    Tests for time conversion between generations and coalescent
    units.
    """
    def test_coalescence_after_growth_rate_change(self):
        Ne = 10000
        # Migrations and bottleneck occured 100 generations ago.
        g = 100
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
        ]
        # At this time, we migrate the lineage in 1 to 0, and
        # have a very strong negative growth rate, resulting in almost instant
        # coalescence.
        demographic_events = [
            msprime.MassMigration(time=g, source=1, dest=0),
            msprime.PopulationParametersChange(time=g, growth_rate=1000),
        ]
        reps = msprime.simulate(
            Ne=Ne,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1, num_replicates=10, record_migrations=True)
        for ts in reps:
            migrations = list(ts.migrations())
            self.assertGreater(len(migrations), 0)
            for mr in migrations:
                self.assertAlmostEqual(g, mr.time)
            tree = next(ts.trees())
            u = tree.get_mrca(0, 1)
            self.assertEqual(u, 2)
            self.assertAlmostEqual(g, tree.time(u), places=1)

    def test_coalescence_after_size_change(self):
        Ne = 20000
        # Migrations and bottleneck occured 100 generations ago.
        g = 1000
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
        ]
        # At this time, we migrate the lineage in 1 to 0, and
        # have a very strong bottleneck, resulting in almost instant
        # coalescence.
        demographic_events = [
            msprime.MassMigration(time=g, source=1, dest=0),
            msprime.PopulationParametersChange(time=g, initial_size=1e-3),
        ]
        reps = msprime.simulate(
            Ne=Ne,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1, num_replicates=10)
        for ts in reps:
            tree = next(ts.trees())
            u = tree.get_mrca(0, 1)
            self.assertEqual(u, 2)
            self.assertAlmostEqual(g, tree.time(u), places=1)

    def test_instantaneous_bottleneck(self):
        Ne = 0.5
        # Bottleneck occured 0.1 coalescent units ago
        t = 0.1
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        # At this time, we migrate the lineages in 1 to 0, and
        # have a very strong bottleneck, resulting in instant
        # coalescence.
        demographic_events = [
            msprime.MassMigration(time=t, source=1, dest=0),
            msprime.InstantaneousBottleneck(time=t, population=0, strength=100)
        ]
        reps = msprime.simulate(
            Ne=Ne,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1, num_replicates=10)
        for ts in reps:
            tree = next(ts.trees())
            self.assertAlmostEqual(t, tree.time(tree.root), places=5)


class TestLowLevelConversions(unittest.TestCase):
    """
    Checks that we convert to the correct low-level values when we
    do the rescalings from generations.
    """
    def test_population_configuration_defaults(self):
        conf = msprime.PopulationConfiguration()
        self.assertIsNone(conf.sample_size)
        d = conf.get_ll_representation()
        dp = {
            "initial_size": None,
            "growth_rate": 0
        }
        self.assertEqual(d, dp)

    def test_population_configuration_initial_size(self):
        for initial_size in [1, 10, 1000]:
            conf = msprime.PopulationConfiguration(initial_size=initial_size)
            self.assertIsNone(conf.sample_size)
            d = conf.get_ll_representation()
            dp = {
                "initial_size": initial_size,
                "growth_rate": 0
            }
            self.assertEqual(d, dp)

    def test_population_configuration_growth_rate(self):
        sample_size = 8
        for growth_rate in [1, 10, -10]:
            conf = msprime.PopulationConfiguration(sample_size, growth_rate=growth_rate)
            self.assertEqual(conf.sample_size, sample_size)
            d = conf.get_ll_representation()
            dp = {
                "initial_size": None,
                "growth_rate": growth_rate
            }
            self.assertEqual(d, dp)

    def test_population_parameters_change_time(self):
        for Ne in [1, 10, 1000]:
            for g in [0.1, 1, 100, 1e6]:
                event = msprime.PopulationParametersChange(time=g, initial_size=Ne)
                d = event.get_ll_representation(1)
                dp = {
                    "time": g,
                    "population": -1,
                    "type": "population_parameters_change",
                    "initial_size": Ne}
                self.assertEqual(d, dp)

    def test_population_parameters_change_initial_size(self):
        g = 100
        for Ne in [1, 10, 1000]:
            for initial_size in [0.01, 1, 100, 1e6]:
                event = msprime.PopulationParametersChange(
                    time=g, initial_size=initial_size)
                d = event.get_ll_representation(1)
                dp = {
                    "time": g,
                    "population": -1,
                    "type": "population_parameters_change",
                    "initial_size": initial_size}
                self.assertEqual(d, dp)

    def test_population_parameters_change_growth_rate(self):
        g = 100
        for growth_rate in [0.01, 1, 100, 1e6]:
            event = msprime.PopulationParametersChange(time=g, growth_rate=growth_rate)
            d = event.get_ll_representation(1)
            dp = {
                "time": g,
                "population": -1,
                "type": "population_parameters_change",
                "growth_rate": growth_rate}
            self.assertEqual(d, dp)

    def test_population_parameters_change_population(self):
        g = 100
        Ne = 10
        for population in range(3):
            event = msprime.PopulationParametersChange(
                time=g, initial_size=Ne, population=population)
            d = event.get_ll_representation(1)
            dp = {
                "time": g,
                "population": population,
                "type": "population_parameters_change",
                "initial_size": Ne}
            self.assertEqual(d, dp)

    def test_migration_rate_change_time(self):
        for g in [0.1, 1, 100, 1e6]:
            event = msprime.MigrationRateChange(time=g, rate=0)
            d = event.get_ll_representation(1)
            dp = {
                "time": g,
                "type": "migration_rate_change",
                "migration_rate": 0,
                "matrix_index": -1}
            self.assertEqual(d, dp)

    def test_migration_rate_change_matrix_index(self):
        g = 51
        for N in range(1, 5):
            for index in itertools.permutations(range(N), 2):
                event = msprime.MigrationRateChange(time=g, rate=0, matrix_index=index)
                d = event.get_ll_representation(N)
                dp = {
                    "time": g,
                    "type": "migration_rate_change",
                    "migration_rate": 0,
                    "matrix_index": index[0] * N + index[1]}
                self.assertEqual(d, dp)

    def test_migration_rate_change_rate(self):
        g = 1234
        for rate in [0, 1e-6, 10, 1e6]:
            event = msprime.MigrationRateChange(time=g, rate=rate)
            d = event.get_ll_representation(1)
            dp = {
                "time": g,
                "type": "migration_rate_change",
                "migration_rate": rate,
                "matrix_index": -1}
            self.assertEqual(d, dp)

    def test_mass_migration_time(self):
        for g in [0.1, 1, 100, 1e6]:
            event = msprime.MassMigration(time=g, source=0, dest=1)
            d = event.get_ll_representation(1)
            dp = {
                "time": g,
                "type": "mass_migration",
                "source": 0,
                "dest": 1,
                "proportion": 1}
            self.assertEqual(d, dp)

    def test_mass_migration_source_dest(self):
        g = 51
        for source, dest in itertools.permutations(range(4), 2):
            event = msprime.MassMigration(time=g, source=source, dest=dest)
            d = event.get_ll_representation(1)
            dp = {
                "time": g,
                "type": "mass_migration",
                "source": source,
                "dest": dest,
                "proportion": 1}
            self.assertEqual(d, dp)

    def test_mass_migration_proportion(self):
        g = 51
        for p in [0, 1e-6, 0.4, 1]:
            event = msprime.MassMigration(time=g, source=0, dest=1, proportion=p)
            d = event.get_ll_representation(1)
            dp = {
                "time": g,
                "type": "mass_migration",
                "source": 0,
                "dest": 1,
                "proportion": p}
            self.assertEqual(d, dp)

    def test_migration_matrix(self):
        m = [
            [0, 1, 2],
            [3, 0, 4],
            [5, 6, 0]]
        sim = msprime.simulator_factory(
            population_configurations=[
                msprime.PopulationConfiguration(1),
                msprime.PopulationConfiguration(1),
                msprime.PopulationConfiguration(1)],
            migration_matrix=m)
        self.assertEqual(sim.migration_matrix, m)

    def test_instantaneous_bottleneck(self):
        g = 51
        for population in [0, 1, 5]:
            for strength in [0, 100, 1000, 1e9]:
                event = msprime.InstantaneousBottleneck(
                    time=g, population=population, strength=strength)
                d = event.get_ll_representation(1)
                dp = {
                    "time": g,
                    "type": "instantaneous_bottleneck",
                    "population": population,
                    "strength": strength}
                self.assertEqual(d, dp)


class TestHistoricalSampling(unittest.TestCase):
    """
    Tests to make sure historical sampling works correctly.
    """
    def test_two_samples(self):
        sampling_time = 1.01
        for recombination_rate in [0, 1]:
            ts = msprime.simulate(
                recombination_rate=recombination_rate,
                samples=[
                    msprime.Sample(0, 0), msprime.Sample(0, sampling_time)])
            for t in ts.trees():
                self.assertEqual(t.get_time(0), 0)
                self.assertEqual(t.get_time(1), sampling_time)
                self.assertEqual(t.get_parent(0), t.get_parent(1))
                self.assertEqual(t.get_parent(1), t.get_parent(0))
                self.assertGreater(t.get_time(t.get_parent(0)), sampling_time)

    def test_different_times(self):
        st1 = 1.01
        st2 = 2.01
        st3 = 3.01
        ts = msprime.simulate(
            samples=[
                msprime.Sample(0, 0),
                msprime.Sample(0, st1),
                msprime.Sample(0, st2),
                msprime.Sample(0, st3)])
        t = next(ts.trees())
        self.assertEqual(t.get_time(0), 0)
        self.assertEqual(t.get_time(1), st1)
        self.assertEqual(t.get_time(2), st2)
        self.assertEqual(t.get_time(3), st3)
        self.assertGreater(t.get_time(t.get_parent(1)), st1)
        self.assertGreater(t.get_time(t.get_parent(2)), st2)
        self.assertGreater(t.get_time(t.get_parent(3)), st3)
        self.assertGreater(t.get_time(t.get_root()), st3)

    def test_old_sampling_time(self):
        # This is an enormously long time in coalescent time, so we should
        # coalesce quickly after the samples are introduced.
        sampling_time = 1000.01
        n = 5
        samples = [
            msprime.Sample(0, sampling_time) for j in range(n - 1)] + [
            msprime.Sample(0, 0)]
        ts = msprime.simulate(Ne=1/4, samples=samples)
        t = next(ts.trees())
        for j in range(n - 1):
            self.assertEqual(t.get_time(j), sampling_time)
        self.assertEqual(t.get_time(n - 1), 0)
        # Allow it to be within 10 coalescent time units.
        self.assertLess(t.get_time(t.get_root()), sampling_time + 10)

    def test_two_samples_mass_migration(self):
        sampling_time = 2.01
        migration_time = 4.33
        ts = msprime.simulate(
            samples=[
                msprime.Sample(0, 0),
                msprime.Sample(1, sampling_time)],
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration()],
            demographic_events=[
                msprime.MassMigration(
                    time=migration_time, source=1, dest=0)])
        t = next(ts.trees())
        self.assertEqual(t.get_time(0), 0)
        self.assertEqual(t.get_time(1), sampling_time)
        self.assertGreater(t.get_time(2), migration_time)
        self.assertEqual(t.get_population(0), 0)
        self.assertEqual(t.get_population(1), 1)
        self.assertEqual(t.get_population(2), 0)

    def test_interleaved_migrations(self):
        t1 = 1.5
        t2 = 10.5
        t3 = 50.5
        ts = msprime.simulate(
            Ne=1/4,
            samples=[
                msprime.Sample(0, 0),
                msprime.Sample(1, t1),
                msprime.Sample(2, t2),
                msprime.Sample(3, t3)],
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration()],
            demographic_events=[
                msprime.MassMigration(time=t1, source=0, dest=1),
                msprime.MassMigration(time=t2, source=1, dest=2),
                msprime.MassMigration(time=t3, source=2, dest=3)])
        t = next(ts.trees())
        self.assertEqual(t.get_time(0), 0)
        self.assertEqual(t.get_time(1), t1)
        self.assertEqual(t.get_time(2), t2)
        self.assertEqual(t.get_time(3), t3)
        self.assertEqual(t.get_population(0), 0)
        self.assertEqual(t.get_population(1), 1)
        self.assertEqual(t.get_population(2), 2)
        self.assertEqual(t.get_population(3), 3)
        self.assertEqual(t.get_population(4), 1)
        self.assertEqual(t.get_population(5), 2)
        self.assertEqual(t.get_population(6), 3)
        self.assertTrue(t1 < t.get_time(4) < t2)
        self.assertTrue(t2 < t.get_time(5) < t3)
        self.assertGreater(t.get_time(6), t3)
