#
# Copyright (C) 2016-2020 University of Oxford
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
import io
import itertools
import json
import math
import random
import unittest
import warnings
from unittest import mock

import numpy as np
import scipy.linalg
import tskit

import msprime
from msprime import _msprime
from msprime import ancestry


class TestNePopulationSizeEquivalence(unittest.TestCase):
    """
    Test that setting Ne as a parameter of the population model and
    at individual populations is the same thing.
    """

    def assert_tree_sequences_equal(self, ts1, ts2):
        t1 = ts1.dump_tables()
        t2 = ts2.dump_tables()
        # We can't compare directly byte-for-byte because there'll be slight
        # differences in the computed times.
        self.assertEqual(len(t1.nodes), len(t2.nodes))
        np.testing.assert_array_equal(t1.nodes.flags, t2.nodes.flags)
        np.testing.assert_array_almost_equal(t1.nodes.time, t2.nodes.time)
        np.testing.assert_array_equal(t1.nodes.population, t2.nodes.population)
        np.testing.assert_array_equal(t1.nodes.individual, t2.nodes.individual)
        np.testing.assert_array_equal(t1.nodes.metadata, t2.nodes.metadata)
        np.testing.assert_array_equal(
            t1.nodes.metadata_offset, t2.nodes.metadata_offset
        )
        self.assertEqual(len(t1.edges), len(t2.edges))
        np.testing.assert_array_almost_equal(t1.edges.left, t2.edges.left)
        np.testing.assert_array_almost_equal(t1.edges.right, t2.edges.right)
        np.testing.assert_array_equal(t1.edges.parent, t2.edges.parent)
        np.testing.assert_array_equal(t1.edges.child, t2.edges.child)

    def test_one_population(self):
        random_seed = 1012
        for Ne in [1e-6, 0.5, 1, 1e6]:
            samples = [msprime.Sample(0, 0)] * 10
            ts1 = msprime.simulate(samples=samples, Ne=Ne, random_seed=random_seed)
            ts2 = msprime.simulate(
                samples=samples,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne)
                ],
                random_seed=random_seed,
            )
            self.assert_tree_sequences_equal(ts1, ts2)
            # Specifying Ne as well should be the same thing
            ts3 = msprime.simulate(
                samples=samples,
                Ne=Ne,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne)
                ],
                random_seed=random_seed,
            )
            self.assert_tree_sequences_equal(ts1, ts3)

    def test_many_populations_equal_size_no_recombination(self):
        random_seed = 1012
        num_pops = 4
        migration_matrix = np.ones((num_pops, num_pops)) * 1e-3
        np.fill_diagonal(migration_matrix, 0)

        for Ne in [0.5, 1, 1.5]:
            samples = [msprime.Sample(population=j, time=0) for j in range(num_pops)]
            ts1 = msprime.simulate(
                samples=samples,
                Ne=Ne,
                population_configurations=[
                    msprime.PopulationConfiguration() for _ in range(num_pops)
                ],
                migration_matrix=migration_matrix,
                random_seed=random_seed,
            )
            ts2 = msprime.simulate(
                samples=samples,
                migration_matrix=migration_matrix,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne)
                    for _ in range(num_pops)
                ],
                random_seed=random_seed,
            )
            self.assert_tree_sequences_equal(ts1, ts2)
            # Also setting initial_size=Ne is the same.
            ts3 = msprime.simulate(
                samples=samples,
                Ne=Ne,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne)
                    for _ in range(num_pops)
                ],
                migration_matrix=migration_matrix,
                random_seed=random_seed,
            )
            self.assert_tree_sequences_equal(ts1, ts3)

    def test_many_populations_equal_size_recombination(self):
        random_seed = 1012234
        num_pops = 3
        recombination_rate = 1
        migration_matrix = np.ones((num_pops, num_pops)) * 1e-3
        np.fill_diagonal(migration_matrix, 0)
        recombination_rate = 0.5

        for Ne in [0.5, 1, 1.5]:
            samples = [
                msprime.Sample(population=j, time=0) for j in range(num_pops)
            ] * 4
            ts1 = msprime.simulate(
                samples=samples,
                Ne=Ne,
                recombination_rate=recombination_rate,
                population_configurations=[
                    msprime.PopulationConfiguration() for _ in range(num_pops)
                ],
                migration_matrix=migration_matrix,
                random_seed=random_seed,
            )
            self.assertGreater(ts1.num_trees, 1)
            ts2 = msprime.simulate(
                samples=samples,
                recombination_rate=recombination_rate,
                migration_matrix=migration_matrix,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne)
                    for _ in range(num_pops)
                ],
                random_seed=random_seed,
            )
            self.assert_tree_sequences_equal(ts1, ts2)
            # Also setting initial_size=Ne is the same.
            ts3 = msprime.simulate(
                samples=samples,
                Ne=Ne,
                recombination_rate=recombination_rate,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne)
                    for _ in range(num_pops)
                ],
                migration_matrix=migration_matrix,
                random_seed=random_seed,
            )
            self.assert_tree_sequences_equal(ts1, ts3)

    def test_many_populations_different_sizes(self):
        random_seed = 1012
        num_pops = 5
        migration_matrix = np.ones((num_pops, num_pops)) * 1e-3
        np.fill_diagonal(migration_matrix, 0)

        for Ne in [0.5, 1, 1.5]:
            samples = [
                msprime.Sample(population=j, time=0) for j in range(num_pops)
            ] * 5
            ts1 = msprime.simulate(
                samples=samples,
                Ne=Ne,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne * (j + 1))
                    for j in reversed(range(num_pops))
                ],
                migration_matrix=migration_matrix,
                random_seed=random_seed,
            )
            ts2 = msprime.simulate(
                samples=samples,
                migration_matrix=migration_matrix,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=Ne * (j + 1))
                    for j in reversed(range(num_pops))
                ],
                random_seed=random_seed,
            )
            self.assert_tree_sequences_equal(ts1, ts2)


class TestIntrospectionInterface(unittest.TestCase):
    """
    Tests that we have meaningful repr and str functions for all the
    classes used in the demography hierarchy.
    """

    def test_population_parameters_change(self):
        event = msprime.PopulationParametersChange(1.0, population=1, initial_size=2.0)
        repr_s = (
            "PopulationParametersChange(time=1.0, initial_size=2.0, "
            "growth_rate=None, population=1)"
        )
        self.assertEqual(repr(event), repr_s)
        str_s = "Population parameter change for 1: initial_size -> 2.0"
        self.assertEqual(str(event), str_s)

        event = msprime.PopulationParametersChange(1.0, population=1, growth_rate=2.0)
        str_s = "Population parameter change for 1: growth_rate -> 2.0"
        self.assertEqual(str(event), str_s)

        event = msprime.PopulationParametersChange(
            1.0, population=1, growth_rate=2.0, initial_size=4.0
        )
        str_s = (
            "Population parameter change for 1: initial_size -> 4.0 growth_rate -> 2.0"
        )
        self.assertEqual(str(event), str_s)

    def test_migration_rate_change(self):
        event = msprime.MigrationRateChange(time=1, rate=2)
        repr_s = "MigrationRateChange(time=1, rate=2, source=-1, dest=-1)"
        str_s = "Migration rate change to 2 everywhere"
        self.assertEqual(repr(event), repr_s)
        self.assertEqual(str(event), str_s)

        event = msprime.MigrationRateChange(time=1, rate=2, source=0, dest=1)
        repr_s = "MigrationRateChange(time=1, rate=2, source=0, dest=1)"
        str_s = "Migration rate change for (0, 1) to 2"
        self.assertEqual(repr(event), repr_s)
        self.assertEqual(str(event), str_s)

    def test_mass_migration(self):
        event = msprime.MassMigration(time=1, proportion=0.5, source=0, dest=1)
        repr_s = "MassMigration(time=1, source=0, dest=1, proportion=0.5)"
        str_s = (
            "Mass migration: Lineages moved with probability 0.5 backwards in time "
            "with source 0 & dest 1 (equivalent to migration from 1 to 0 "
            "forwards in time)"
        )
        self.assertEqual(repr(event), repr_s)
        # Too much hassle to track the exact whitespace in the output string.
        self.assertEqual(str(event).split(), str_s.split())

    def test_simple_bottleneck(self):
        event = msprime.SimpleBottleneck(time=1, population=1, proportion=0.5)
        repr_s = "SimpleBottleneck(time=1, population=1, proportion=0.5)"
        str_s = "Simple bottleneck: lineages in population 1 coalesce probability 0.5"
        self.assertEqual(repr(event), repr_s)
        self.assertEqual(str(event), str_s)

    def test_instantaneous_bottleneck(self):
        event = msprime.InstantaneousBottleneck(time=1, population=1, strength=1.5)
        repr_s = "InstantaneousBottleneck(time=1, population=1, strength=1.5)"
        str_s = (
            "Instantaneous bottleneck in population 1: equivalent to 1.5 "
            "generations of the coalescent"
        )
        self.assertEqual(repr(event), repr_s)
        self.assertEqual(str(event), str_s)

    def test_census(self):
        event = msprime.CensusEvent(time=1)
        repr_s = "CensusEvent(time=1)"
        str_s = "Census event"
        self.assertEqual(repr(event), repr_s)
        self.assertEqual(str(event), str_s)


class TestDemographicEventsHaveExtraLLParameter(unittest.TestCase):
    """
    For legacy reasons the DemographicEvent.get_ll_representation took
    a "num_populations" argument. In versions of stdpopsim using
    msprime < 1.0, this parameter was specified when testing for
    model equality (even though it did nothing). To ensure we're not
    breaking older versions of stdpopsim, we keep this extra parameter
    and test it here to make sure it works.

    Since the functionality is only really used as a developer tool
    in stdpopsim, we can get rid of the extra parameters once
    stdpopsim 0.2 (which won't use this API) has been released.

    See https://github.com/tskit-dev/msprime/issues/1037
    """

    def test_demographic_events_have_param(self):
        events = [
            msprime.PopulationParametersChange(1.0, population=1, initial_size=2.0),
            msprime.MigrationRateChange(1.0, 1.0),
            msprime.MassMigration(1.0, 0),
            msprime.SimpleBottleneck(1.0, 0),
            msprime.InstantaneousBottleneck(1.0),
            msprime.CensusEvent(1.0),
        ]
        for event in events:
            ll_config1 = event.get_ll_representation()
            ll_config2 = event.get_ll_representation(None)
            self.assertEqual(ll_config1, ll_config2)


class TestTimeTravelErrors(unittest.TestCase):
    """
    It is possible to specify models in msprime that result in malformed
    tree sequences where the parent node has time equal to its child.
    We throw an error in this case and expect the user to fix their model.
    """

    def test_multiple_bottlenecks(self):
        with self.assertRaises(_msprime.LibraryError):
            for model in ["hudson", "smc"]:
                msprime.simulate(
                    model=model,
                    sample_size=100,
                    demographic_events=[
                        msprime.SimpleBottleneck(
                            time=0.1, population=0, proportion=0.75
                        ),
                        msprime.SimpleBottleneck(
                            time=0.1, population=0, proportion=1.0
                        ),
                    ],
                    random_seed=1,
                )

    def test_tiny_population_size(self):
        # Derived from bug report in #570.
        n = 3
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=10, initial_size=10000, growth_rate=0
            )
            for k in range(n)
        ]
        demographic_events = [
            msprime.PopulationParametersChange(
                time=0.00001, initial_size=1e-18, population_id=2, growth_rate=0
            ),
            msprime.MassMigration(time=0.02, source=1, destination=2, proportion=1.0),
            msprime.MigrationRateChange(time=0.02, rate=0),
        ]
        M = [[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
        with self.assertRaises(_msprime.LibraryError):
            msprime.simulate(
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                migration_matrix=M,
                recombination_rate=0.0,
                mutation_rate=0.0,
                random_seed=1,
            )


class TestBadDemographicParameters(unittest.TestCase):
    """
    Tests for nonsensical demographic parameters.
    """

    def test_bad_population_size(self):
        for bad_size in [-1, -1e300]:
            self.assertRaises(
                ValueError, msprime.PopulationParametersChange, 0, initial_size=bad_size
            )
            self.assertRaises(
                ValueError,
                msprime.PopulationConfiguration,
                initial_size=bad_size,
                sample_size=2,
            )

    def test_bad_sample_size(self):
        for bad_size in [-1, -1e300]:
            self.assertRaises(
                ValueError,
                msprime.PopulationConfiguration,
                initial_size=1,
                sample_size=bad_size,
            )

    def test_dtwf_bottleneck(self):
        with self.assertRaises(_msprime.LibraryError):
            msprime.simulate(
                sample_size=2,
                model="dtwf",
                demographic_events=[msprime.SimpleBottleneck(time=0.1, population=0)],
                random_seed=1,
            )


class TestBadDemographicEvents(unittest.TestCase):
    """
    Tests for input errors when creating demographic events.
    """

    def test_growth_rate_or_initial_size(self):
        self.assertRaises(ValueError, msprime.PopulationParametersChange, time=0)

    def test_bad_simulation_model(self):
        for model in [[], {}]:
            des = [msprime.SimulationModelChange(time=0, model=model)]
            with self.assertRaises(TypeError):
                msprime.simulate(10, demographic_events=des)


class TestZeroPopulationSize(unittest.TestCase):
    """
    Tests to check that a population with zero size is correctly handled.
    """

    def setUp(self):
        # Build a three-population demography, where pop 0 is a strict ancestor
        # of pop 1 and pop 2, such that pop 0 goes extinct when pop 1 and pop 2
        # come into existence.
        # Each population's size is set to zero outside its existence time.
        self.N = 1000
        self.T = 500
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=0),
            msprime.PopulationConfiguration(initial_size=self.N),
            msprime.PopulationConfiguration(initial_size=self.N),
        ]
        self.demographic_events = [
            msprime.MassMigration(self.T, source=1, dest=0),
            msprime.MassMigration(self.T, source=2, dest=0),
            msprime.PopulationParametersChange(
                self.T, initial_size=self.N, population_id=0
            ),
            msprime.PopulationParametersChange(self.T, initial_size=0, population_id=1),
            msprime.PopulationParametersChange(self.T, initial_size=0, population_id=2),
        ]

    def test_demography_with_zero_sizes(self):
        for samples in [
            [msprime.Sample(1, 0)],
            [msprime.Sample(1, 0), msprime.Sample(2, 0)],
            # Check ancient sampling.
            [msprime.Sample(0, self.T), msprime.Sample(1, 0)],
        ]:
            ts = msprime.simulate(
                population_configurations=self.population_configurations,
                demographic_events=self.demographic_events,
                samples=samples * 10,
            )
            self.assertTrue(
                all(tree.num_roots == 1 for tree in ts.trees()), msg=f"{samples}"
            )

    def test_cant_sample_a_zero_sized_population(self):
        for bad_samples, err in [
            ([msprime.Sample(0, 0)], _msprime.InputError),
            # Check ancient sampling.
            ([msprime.Sample(1, self.T), msprime.Sample(1, 0)], _msprime.LibraryError),
        ]:
            with self.assertRaises(err, msg=f"{bad_samples}"):
                msprime.simulate(
                    population_configurations=self.population_configurations,
                    demographic_events=self.demographic_events,
                    samples=bad_samples * 10,
                )

    def test_demography_debugger(self):
        msprime.DemographyDebugger(
            population_configurations=self.population_configurations,
            demographic_events=self.demographic_events,
        )

        with self.assertRaises(ValueError):
            msprime.DemographyDebugger(
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=0)
                ]
            )


class TestDeprecatedParameters(unittest.TestCase):
    """
    Tests to check that aliased parameters are handled correctly.
    """

    def test_mass_migration_dest(self):
        self.assertRaises(
            ValueError, msprime.MassMigration, time=0, source=0, dest=0, destination=0
        )
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
            ValueError,
            msprime.PopulationParametersChange,
            time=0,
            initial_size=1.1,
            growth_rate=0.1,
            population=1,
            population_id=1,
        )

        for j in range(10):
            e = msprime.PopulationParametersChange(
                time=0.1, initial_size=1.1, growth_rate=0.1, population=j
            )
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
                time=0.1, initial_size=1.1, growth_rate=0.1, population_id=j
            )
            self.assertEqual(e.time, 0.1)
            self.assertEqual(e.initial_size, 1.1)
            self.assertEqual(e.growth_rate, 0.1)
            self.assertEqual(e.population, j)


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
            "initial_size": new_size,
        }
        self.assertEqual(event.get_ll_representation(), ll_event)

    def test_growth_rate_change(self):
        g = 512
        growth_rate = 1
        event = msprime.PopulationParametersChange(
            time=g, growth_rate=growth_rate, population=1
        )
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": 1,
            "growth_rate": growth_rate,
        }
        self.assertEqual(event.get_ll_representation(), ll_event)

    def test_growth_rate_and_size_change(self):
        g = 1024
        growth_rate = 2
        initial_size = 8192
        event = msprime.PopulationParametersChange(
            time=g, initial_size=initial_size, growth_rate=growth_rate, population=1
        )
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": 1,
            "initial_size": initial_size,
            "growth_rate": growth_rate,
        }
        self.assertEqual(event.get_ll_representation(), ll_event)

    def test_migration_rate_change(self):
        g = 1024
        migration_rate = 0.125
        event = msprime.MigrationRateChange(time=g, rate=migration_rate)
        ll_event = {
            "type": "migration_rate_change",
            "time": g,
            "source": -1,
            "dest": -1,
            "migration_rate": migration_rate,
        }
        self.assertEqual(event.get_ll_representation(), ll_event)


class TestDemographyDebuggerOutput(unittest.TestCase):
    """
    Tests for the demography debug interface.
    """

    def verify_debug(
        self, population_configurations, migration_matrix, demographic_events
    ):

        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events,
        )
        # Check the reprs
        s = repr(dd.epochs)
        self.assertGreater(len(s), 0)
        buff = io.StringIO()
        dd.print_history(buff)
        debug_output = buff.getvalue()
        self.assertEqual(debug_output, str(dd))
        # TODO when there is better output, write some tests to
        # verify its format.
        self.assertGreater(len(debug_output), 0)

    def test_zero_samples(self):
        population_configurations = [msprime.PopulationConfiguration(0)]
        self.verify_debug(population_configurations, [[0]], [])
        self.verify_debug(None, None, [])

    def test_one_population(self):
        population_configurations = [msprime.PopulationConfiguration(10)]
        migration_matrix = [[0]]
        demographic_events = [
            msprime.PopulationParametersChange(0.1, initial_size=2),
            msprime.PopulationParametersChange(0.1, growth_rate=10),
        ]
        self.verify_debug(
            population_configurations, migration_matrix, demographic_events
        )

    def test_no_events(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        migration_matrix = [[0, 0], [0, 0]]
        self.verify_debug(population_configurations, migration_matrix, [])

    def test_demographic_events(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        migration_matrix = [[0, 0], [0, 0]]
        demographic_events = [
            msprime.PopulationParametersChange(0.1, initial_size=2),
            msprime.PopulationParametersChange(0.1, growth_rate=10),
            msprime.MigrationRateChange(0.2, matrix_index=(0, 1), rate=1),
            msprime.MigrationRateChange(0.2, matrix_index=(1, 0), rate=1),
            msprime.MassMigration(0.4, source=1, dest=0),
            msprime.MigrationRateChange(0.4, rate=0),
            msprime.InstantaneousBottleneck(0.5, population=0, strength=100),
            msprime.CensusEvent(0.55),
        ]
        self.verify_debug(
            population_configurations, migration_matrix, demographic_events
        )


class TestDemographyDebugger(unittest.TestCase):
    """
    Tests for the demography debugger. Ensure that we compute the correct
    population sizes etc.
    """

    def verify_arrays(self, dd):
        """
        Check that the array properties that we generate are computed correctly.
        """
        pop_size = dd.population_size_history
        times = dd.epoch_times
        self.assertEqual(dd.num_epochs, times.shape[0])
        self.assertEqual(dd.num_epochs, pop_size.shape[1])
        self.assertEqual(dd.num_populations, pop_size.shape[0])
        for j in range(dd.num_epochs):
            self.assertEqual(dd.epochs[j].start_time, times[j])
            for k in range(dd.num_populations):
                self.assertEqual(dd.epochs[j].populations[k].start_size, pop_size[k, j])

    def test_equal_after(self):
        population_configurations = [
            msprime.PopulationConfiguration(sample_size=10),
            msprime.PopulationConfiguration(sample_size=20),
        ]
        msprime.DemographyDebugger(population_configurations=population_configurations)
        self.assertEqual(population_configurations[0].sample_size, 10)
        self.assertEqual(population_configurations[1].sample_size, 20)

    def test_model_change_events(self):
        population_configurations = [msprime.PopulationConfiguration(sample_size=10)]
        demographic_events = [msprime.SimulationModelChange(1, "hudson")]
        with self.assertRaises(TypeError):
            msprime.DemographyDebugger(
                population_configurations=population_configurations,
                demographic_events=demographic_events,
            )

    def test_one_pop_zero_events(self):
        dd = msprime.DemographyDebugger(
            population_configurations=[msprime.PopulationConfiguration()]
        )
        self.assertEqual(len(dd.epochs), 1)
        self.verify_arrays(dd)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertEqual(dd.epoch_times[0], 0)
        self.assertEqual(dd.population_size_history.shape[0], 1)
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
                msprime.PopulationConfiguration(initial_size=20),
            ]
        )
        self.verify_arrays(dd)
        self.assertEqual(len(dd.epochs), 1)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertEqual(dd.population_size_history.shape[0], 2)
        self.assertEqual(dd.population_size_history[0][0], 10)
        self.assertEqual(dd.population_size_history[1][0], 20)
        self.assertTrue(math.isinf(e.end_time))
        self.assertEqual(len(e.demographic_events), 0)
        self.assertEqual(len(e.populations), 2)
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
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
                msprime.PopulationConfiguration(initial_size=20, growth_rate=g2),
            ],
            demographic_events=[
                msprime.PopulationParametersChange(time=10, growth_rate=0)
            ],
        )
        self.verify_arrays(dd)
        # Make sure we're testing the __repr__ paths.
        s = repr(dd)
        self.assertGreater(len(s), 0)
        self.assertEqual(len(dd.epochs), 2)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertEqual(e.end_time, 10)
        self.assertEqual(len(e.demographic_events), 0)
        self.assertEqual(len(e.populations), 2)
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
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
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
        self.assertEqual(e.populations[0].start_size, p0_end_size)
        self.assertEqual(e.populations[1].start_size, p1_end_size)

    def test_two_pop_update_migration_rate(self):
        dd = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=100),
                msprime.PopulationConfiguration(initial_size=100),
            ],
            demographic_events=[
                msprime.MigrationRateChange(time=20, rate=1),
                msprime.MigrationRateChange(time=22, rate=1.7, matrix_index=(0, 1)),
            ],
            migration_matrix=[[0, 0.25], [0, 0]],
        )
        self.verify_arrays(dd)
        self.assertEqual(len(dd.epochs), 3)
        e = dd.epochs[0]
        self.assertEqual(e.start_time, 0)
        self.assertEqual(e.end_time, 20)
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0.25], [0, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
            self.assertEqual(pop.start_size, 100)

        e = dd.epochs[1]
        self.assertEqual(e.start_time, 20)
        self.assertTrue(e.end_time, 22)
        self.assertEqual(len(e.demographic_events), 1)
        d = e.demographic_events[0]
        self.assertEqual(d.source, -1)
        self.assertEqual(d.dest, -1)
        self.assertEqual(d.time, 20)
        self.assertEqual(d.rate, 1)
        self.assertEqual(len(e.populations), 2)
        np.testing.assert_array_equal(e.migration_matrix, [[0, 1], [1, 0]])
        for pop in e.populations:
            self.assertEqual(pop.growth_rate, 0)
            self.assertEqual(pop.start_size, pop.end_size)
            self.assertEqual(pop.start_size, 100)

        e = dd.epochs[2]
        self.assertEqual(e.start_time, 22)
        self.assertTrue(math.isinf(e.end_time))
        self.assertEqual(len(e.demographic_events), 1)
        d = e.demographic_events[0]
        self.assertEqual(d.source, 0)
        self.assertEqual(d.dest, 1)
        self.assertEqual(d.time, 22)
        self.assertEqual(d.rate, 1.7)
        self.assertEqual(len(e.populations), 2)
        np.testing.assert_array_equal(e.migration_matrix, [[0, 1.7], [1, 0]])
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
        for model in ["dtwf", "hudson"]:
            dd = msprime.DemographyDebugger(
                model=model,
                population_configurations=[
                    msprime.PopulationConfiguration(initial_size=N0, growth_rate=alpha),
                    msprime.PopulationConfiguration(initial_size=N1, growth_rate=0),
                ],
                demographic_events=[
                    # p1 changes growth rate to alpha
                    msprime.PopulationParametersChange(
                        population=1, time=t1, growth_rate=alpha
                    ),
                    # p0 changes growth rate to -alpha
                    msprime.PopulationParametersChange(
                        population=0, time=t2, growth_rate=-alpha
                    ),
                    # Both change growth_rate to 0 at t3.
                    msprime.PopulationParametersChange(time=t3, growth_rate=0),
                ],
            )
            self.verify_arrays(dd)

            self.assertEqual(len(dd.epochs), 4)
            e = dd.epochs[0]
            self.assertEqual(e.start_time, 0)
            self.assertEqual(e.end_time, t1)
            self.assertEqual(dd.epoch_times[0], 0)
            self.assertEqual(dd.epoch_times[1], t1)
            self.assertEqual(len(e.demographic_events), 0)
            self.assertEqual(len(e.populations), 2)
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
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
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
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
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
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
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
            self.assertEqual(e.populations[0].start_size, n0)
            self.assertEqual(e.populations[0].end_size, n0)
            self.assertEqual(e.populations[1].start_size, n1)
            self.assertEqual(e.populations[1].end_size, n1)

    def check_model_misspecification_warning(self, misspecify):
        population_configurations = [
            msprime.PopulationConfiguration(initial_size=1000),
            msprime.PopulationConfiguration(initial_size=1000),
        ]
        migration_matrix = [[0, 1e-5], [0, 0]]
        demographic_events = [
            msprime.MassMigration(time=1000, dest=0, source=1, proportion=1),
        ]
        if not misspecify:
            demographic_events.append(msprime.MigrationRateChange(time=1000, rate=0))
        msprime.DemographyDebugger(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            migration_matrix=migration_matrix,
        )

    def test_misspecification_warning(self):
        with warnings.catch_warnings(record=True) as w:
            self.check_model_misspecification_warning(misspecify=True)
            self.assertEqual(len(w), 1)

    def test_no_misspecification_no_warning(self):
        with warnings.catch_warnings(record=True) as w:
            self.check_model_misspecification_warning(misspecify=False)
            self.assertEqual(len(w), 0)


class TestDemographyTrajectories(unittest.TestCase):
    """
    Tests that methods msprime.DemographyDebugger.population_size_trajectory
    and msprime.DemographyDebugger.coalescence_rate_trajectory are giving
    correct predictions for trivial demographic models
    """

    def subdivide(self, steps):
        inter = steps[:-1] + np.diff(steps) / 2
        double_steps = sorted(np.concatenate([steps, inter]))
        return double_steps

    def one_pop_example(self):
        # One population:
        # - population size today is 100
        # - size change to 200 at t=100 generations ago
        # - population 1 has a size change to 100 at t=300 generations ago
        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=100)
            ],
            demographic_events=[
                msprime.PopulationParametersChange(
                    time=100, initial_size=200, population_id=0
                ),
                msprime.PopulationParametersChange(time=300, initial_size=100),
            ],
        )
        return ddb

    def test_one_pop(self):
        # a) The coalescence rates through time of two samples are:
        #     (1/(2 * 100) until time 100),
        #     (1/(2 * 200) between times 100 and 300),
        #     (1/(2 * 100) past t=300)
        # b) The probabilities of not yet having coalesced for two samples:
        #     (exp(-t/200) until time 100)
        #     (exp(-100/200 - (t-100)/400) between times 100 and 300)
        #     (exp(-100/200 - (300-100)/400 - (t-300)/200) past t=300)
        ddb = self.one_pop_example()
        steps = np.linspace(0, 400, 101)
        rates, P = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2])
        for time, rA, pA in zip(steps, rates, P):
            if time >= 0 and time < 100:
                self.assertAlmostEqual(rA, 1 / (2 * 100))
                self.assertAlmostEqual(pA, np.exp(-time / 200))
            elif time >= 100 and time < 200:
                self.assertAlmostEqual(rA, 1 / (2 * 200))
                self.assertAlmostEqual(pA, np.exp(-100 / 200 - (time - 100) / 400))
            elif time >= 200 and time < 300:
                self.assertAlmostEqual(rA, 1 / (2 * 200))
                self.assertAlmostEqual(pA, np.exp(-100 / 200 - (time - 100) / 400))
            else:
                self.assertAlmostEqual(rA, 1 / (2 * 100))
                self.assertAlmostEqual(
                    pA, np.exp(-100 / 200 - (300 - 100) / 400 - (time - 300) / 200)
                )

    def test_nan_coalrate(self):
        # The returned rate will be nan at times where the probability
        # of not having coalesced is floating-point-zero
        ddb = self.one_pop_example()
        steps = np.array([1000000 * k for k in range(1, 4)])
        rates, P = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2])
        for rA, pA in zip(rates, P):
            self.assertTrue(np.isnan(rA))
            self.assertEqual(pA, 0)

    def test_mean_one_pop(self):
        # As above. The mean time to coalescence should be
        #  200 * (1 - exp(-100/200))
        #  + 400 * exp(-100/200) * (1 - exp(-(300-100)/400))
        #  + 200 * exp(-100/200 - (300-100)/400)
        ddb = self.one_pop_example()
        m = ddb.mean_coalescence_time(num_samples=[2])
        truth = (
            200 * (1 - np.exp(-100 / 200))
            + 400 * np.exp(-100 / 200) * (1 - np.exp(-(300 - 100) / 400))
            + 200 * np.exp(-100 / 200 - (300 - 100) / 400)
        )
        self.assertLess(abs(m - truth) / truth, 5e-3)
        # test passing steps
        m = ddb.mean_coalescence_time(num_samples=[2], steps=[0, 100, 1000])
        self.assertLess(abs(m - truth) / truth, 5e-3)

    def test_logging(self):
        ddb = self.one_pop_example()
        with mock.patch("msprime.demography.logger.debug") as mocked_debug:
            ddb.mean_coalescence_time([2])
        self.assertEqual(mocked_debug.call_count, 3)

    def test_mean_errors(self):
        ddb = self.one_pop_example()
        self.assertRaises(
            ValueError,
            ddb.mean_coalescence_time,
            num_samples=[2],
            steps=[0, 10],
            max_iter=1,
        )

    def two_pop_example(self):
        # Have a history with two populations, with no migration;
        # check that coalescence rate of lineages from
        # a single pop reflects the size history of that pop:
        # - population sizes today are [100, 1000]
        # - population 1 had a size change to 200 at t=100 generations ago
        # - there is a mass migration event that entirely merges population 2 (source)
        # into population 1 (dest) at t=200 generations ago
        # - population 1 has a size change to 100 at t=300 generations ago
        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=1e2),
                msprime.PopulationConfiguration(initial_size=1e3),
            ],
            demographic_events=[
                msprime.PopulationParametersChange(
                    time=100, initial_size=200, population_id=0
                ),
                msprime.MassMigration(time=200, source=1, dest=0),
                msprime.PopulationParametersChange(time=300, initial_size=100),
            ],
            migration_matrix=[[0, 0], [0, 0]],
        )
        return ddb

    def test_two_pop(self):
        # Then, we can test the following:
        # a) The coalescence rates through time of two samples from population 1 is:
        #     (1/(2 * 100) until time 100),
        #     (1/(2 * 200) between times 100 and 300),
        #     (1/(2 * 100) past t=300)
        # b) The coalescence rates through time of two samples from population 2 is:
        #     (1/(2 * 1000) until time 200),
        #     (1/(2 * 200) between times 200 and 300),
        #     (1/(2 * 100) past t=300)
        # c) The coalescence rates through time of one sample from population 1 and
        #     one sample from population 2 is:
        #     (0 until time 200),
        #     (1/(2 * 200) between times 200 and 300),
        #     (1/(2 * 100) past t=300)
        # d) The probabilities of not yet having coalesced for two samples
        #     from population 1 are:
        #     (exp(-t/200) until time 100)
        #     (exp(-100/200 - (t-100)/400) between times 100 and 300)
        #     (exp(-100/200 - (300-100)/400 - (t-300)/200) past t=300)
        # e) The probabilities of not yet having coalesced for two samples
        #     from population 2 are:
        #     (exp(-t/2000) until time 200)
        #     (exp(-200/2000 - (t-200)/400) between times 200 and 300)
        #     (exp(-200/2000 - (300-200)/400 - (t-300)/200) past t=300)
        # f) The probabilities of not yet having coalesced for one sample
        #     from population 1 and one from population 2 are;
        #     (1.0 until time 200)
        #     (exp(- (t-200)/400) between times 200 and 300)
        #     (exp(-(300-200)/400 - (t-300)/200) past t=300)
        # Also note that when comparing to P, the time intervals should include the
        # right endpoints but for comparing to rates, they should include the
        # left endpoints.
        ddb = self.two_pop_example()
        steps = np.linspace(0, 400, 1001)
        rates, P = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2, 0])
        rates1, P1 = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[0, 2])
        rates2, P2 = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[1, 1])
        for time_step in range(len(steps)):
            time = steps[time_step]
            rA = rates[time_step]
            rB = rates1[time_step]
            rC = rates2[time_step]
            pA = P[time_step]
            pB = P1[time_step]
            pC = P2[time_step]
            if time >= 0 and time < 100:
                self.assertAlmostEqual(rA, 1 / (2 * 100))
                self.assertAlmostEqual(rB, 1 / (2 * 1000))
                self.assertAlmostEqual(rC, 0)
                self.assertAlmostEqual(pA, np.exp(-time / 200))
                self.assertAlmostEqual(pB, np.exp(-time / 2000))
                self.assertAlmostEqual(pC, 1.0)
            elif time == 100:
                self.assertAlmostEqual(rA, 1 / (2 * 200))
                self.assertAlmostEqual(rB, 1 / (2 * 1000))
                self.assertAlmostEqual(rC, 0)
                self.assertAlmostEqual(pA, np.exp(-time / 200))
                self.assertAlmostEqual(pB, np.exp(-time / 2000))
                self.assertAlmostEqual(pC, 1.0)
            elif time > 100 and time < 200:
                self.assertAlmostEqual(rA, 1 / (2 * 200))
                self.assertAlmostEqual(rB, 1 / (2 * 1000))
                self.assertAlmostEqual(rC, 0)
                self.assertAlmostEqual(pA, np.exp(-100 / 200 - (time - 100) / 400))
                self.assertAlmostEqual(pB, np.exp(-time / 2000))
                self.assertAlmostEqual(pC, 1.0)
            elif time == 200:
                self.assertAlmostEqual(rA, 1 / (2 * 200))
                self.assertAlmostEqual(rB, 1 / (2 * 200))
                self.assertAlmostEqual(rC, 1 / (2 * 200))
                self.assertAlmostEqual(pA, np.exp(-100 / 200 - (time - 100) / 400))
                self.assertAlmostEqual(pB, np.exp(-time / 2000))
                self.assertAlmostEqual(pC, 1.0)
            elif time > 200 and time < 300:
                self.assertAlmostEqual(rA, 1 / (2 * 200))
                self.assertAlmostEqual(rB, 1 / (2 * 200))
                self.assertAlmostEqual(rC, 1 / (2 * 200))
                self.assertAlmostEqual(pA, np.exp(-100 / 200 - (time - 100) / 400))
                self.assertAlmostEqual(pB, np.exp(-200 / 2000 - (time - 200) / 400))
                self.assertAlmostEqual(pC, np.exp(-(time - 200) / 400))
            elif time == 300:
                self.assertAlmostEqual(rA, 1 / (2 * 100))
                self.assertAlmostEqual(rB, 1 / (2 * 100))
                self.assertAlmostEqual(rC, 1 / (2 * 100))
                self.assertAlmostEqual(pA, np.exp(-100 / 200 - (time - 100) / 400))
                self.assertAlmostEqual(pB, np.exp(-200 / 2000 - (time - 200) / 400))
                self.assertAlmostEqual(pC, np.exp(-(time - 200) / 400))
            else:
                self.assertAlmostEqual(rA, 1 / (2 * 100))
                self.assertAlmostEqual(rB, 1 / (2 * 100))
                self.assertAlmostEqual(rC, 1 / (2 * 100))
                self.assertAlmostEqual(
                    pA, np.exp(-100 / 200 - (300 - 100) / 400 - (time - 300) / 200)
                )
                self.assertAlmostEqual(
                    pB, np.exp(-200 / 2000 - (300 - 200) / 400 - (time - 300) / 200)
                )
                self.assertAlmostEqual(
                    pC, np.exp(-(300 - 200) / 400 - (time - 300) / 200)
                )

    def test_mean_two_pop(self):
        # As above. The mean times to coalescence should be:
        # pop 1:  200 * (1 - exp(-100/200))
        #        + 400 * exp(-100/200) * (1 - exp(-(300-100)/400))
        #        + 200 * exp(-100/200-(300-100)/400)
        # pop 2:  2000 * (1 - exp(-200/2000))
        #        + 400 * exp(-200/2000) * (1 - exp(-(300-200)/400))
        #        + 200 * exp(-200/2000-(300-200)/400)
        # pop 1-2: 200 + 400 * (1 - exp(-(300-200)/400))
        #        + 200 * exp(-(300-200)/400)
        ddb = self.two_pop_example()
        mA = ddb.mean_coalescence_time(num_samples=[2, 0])
        mB = ddb.mean_coalescence_time(num_samples=[0, 2])
        mC = ddb.mean_coalescence_time(num_samples=[1, 1])
        tA = (
            200 * (1 - np.exp(-100 / 200))
            + 400 * np.exp(-100 / 200) * (1 - np.exp(-(300 - 100) / 400))
            + 200 * np.exp(-100 / 200 - (300 - 100) / 400)
        )
        tB = (
            2000 * (1 - np.exp(-200 / 2000))
            + 400 * np.exp(-200 / 2000) * (1 - np.exp(-(300 - 200) / 400))
            + 200 * np.exp(-200 / 2000 - (300 - 200) / 400)
        )
        tC = (
            200
            + 400 * (1 - np.exp(-(300 - 200) / 400))
            + 200 * np.exp(-(300 - 200) / 400)
        )
        self.assertLess(abs(mA - tA) / tA, 5e-3)
        self.assertLess(abs(mB - tB) / tB, 5e-3)
        self.assertLess(abs(mC - tC) / tC, 5e-3)

    def test_constant_sizes(self):
        # With constant population sizes, results should not depend on the steps
        N_A = 1e2
        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=N_A)
            ],
            demographic_events=[
                msprime.PopulationParametersChange(time=100, initial_size=200),
                msprime.PopulationParametersChange(time=300, initial_size=100),
            ],
        )
        steps = np.linspace(0, 400, 21)
        rates, P = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2])
        steps2 = self.subdivide(steps)
        rates2, P2 = ddb.coalescence_rate_trajectory(steps=steps2, num_samples=[2])
        assert np.all(steps == steps2[::2])
        self.assertLess(max(np.abs(rates - rates2[::2])), 1e-6)
        self.assertLess(max(np.abs(P - P2[::2])), 1e-6)

    def test_convergence(self):
        # Have two populations with very high migration; check they act as a single
        # population:
        # - population sizes today are [100, 1000]
        # - migration rate from population 1 to 2 of 20, and from 2 to 1 of 10
        #     Then, we test:
        # a) The coalescence rate for two samples from population 1 is
        #     (1/3) * (1/(2100)) + (2/3) * (1/(21000))
        # b) Same for two samples from population 2 and for one sample from each
        #     population.
        N_A = 1e2
        N_B = 1e3
        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=N_A),
                msprime.PopulationConfiguration(initial_size=N_B),
            ],
            demographic_events=[],
            migration_matrix=[[0, 0.5], [0.25, 0]],
        )
        steps = np.linspace(0, 400, 401)
        rates, PP = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2, 0])
        rates1, PP = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[0, 2])
        rates2, PP = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[1, 1])
        for time_step in range(len(steps)):
            time = steps[time_step]
            rA = rates[time_step]
            rB = rates1[time_step]
            rC = rates2[time_step]
            if time > 100:
                r = ((1 / 3) ** 2) * (1 / (2 * 100)) + ((2 / 3) ** 2) * (1 / (2 * 1000))
                self.assertAlmostEqual(rA, r, places=4)
                self.assertAlmostEqual(rB, r, places=4)
                self.assertAlmostEqual(rC, r, places=4)

    def test_rate_size_equality(self):
        # This tests a trivial model with two populations and no migration.
        # The population size trajectories should match the
        # given 1/ 2*coalescent rates for the respective sampled population.
        N_A = 1e7
        N_B = 1e6
        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=N_A),
                msprime.PopulationConfiguration(initial_size=N_B),
            ],
            demographic_events=[
                msprime.PopulationParametersChange(
                    time=100, initial_size=1e5, population_id=0, growth_rate=3e-3
                ),
                msprime.PopulationParametersChange(
                    time=300, initial_size=1e6, population_id=0, growth_rate=7e-4
                ),
            ],
            migration_matrix=[[0, 0], [0, 0]],
        )
        steps = np.linspace(0, 400, 1001)
        rates, P = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2, 0])
        pop_sizes = ddb.population_size_trajectory(steps=steps)
        for r, p in zip(rates, pop_sizes[:, 0]):
            self.assertAlmostEqual(r, 1 / (2 * p))
        rates, P = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[0, 2])
        for r, p in zip(rates, pop_sizes[:, 1]):
            self.assertAlmostEqual(r, 1 / (2 * p))

    def test_rate_size_equality_with_population_merge(self):
        # Test equality between two population that split with continued
        # migration
        N_A = 1e7
        N_B = 1e6

        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=N_A),
                msprime.PopulationConfiguration(initial_size=N_B),
            ],
            demographic_events=[
                msprime.PopulationParametersChange(
                    time=100, initial_size=1e5, population_id=0, growth_rate=0
                ),
                msprime.PopulationParametersChange(
                    time=100, initial_size=1e5, population_id=1, growth_rate=7e-4
                ),
                msprime.MassMigration(time=100, source=1, dest=0),
                msprime.PopulationParametersChange(
                    time=200, initial_size=1e6, population_id=0, growth_rate=-1e-2
                ),
                msprime.PopulationParametersChange(
                    time=300, initial_size=1e6, population_id=0, growth_rate=0
                ),
            ],
            migration_matrix=[[0, 0], [0, 0]],
        )
        steps = np.linspace(0, 400, 401)
        rates, P = ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2, 0])
        pop_sizes = ddb.population_size_trajectory(steps=steps)
        for r, p in zip(rates[301:], pop_sizes[:, 0][301:]):
            self.assertAlmostEqual(1 / (2 * r), p)

    def test_double_step_validation(self):
        # Test that the double step validation throws a warning
        # with small step sizes
        N_A = 1e3
        N_B = 1e4
        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=N_A),
                msprime.PopulationConfiguration(initial_size=N_B),
            ],
            demographic_events=[
                msprime.PopulationParametersChange(
                    time=100, initial_size=1e2, population_id=0, growth_rate=3e-3
                ),
                msprime.PopulationParametersChange(
                    time=300, initial_size=1e3, population_id=0, growth_rate=7e-4
                ),
            ],
            migration_matrix=[[0, 0], [0, 0]],
        )
        steps = np.linspace(0, 400, 2)
        with self.assertWarns(UserWarning):
            ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2, 0])
        # Test coalescence rates without double step validation
        steps = np.linspace(0, 400, 401)
        rates, P = ddb.coalescence_rate_trajectory(
            steps=steps, num_samples=[2, 0], double_step_validation=False
        )

    def test_value_errors(self):
        # test all user input domains which should raise ValuErrors.
        ddb = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=1e2),
            ],
            demographic_events=[],
            migration_matrix=[[0]],
        )
        steps = np.linspace(0, 10, 11)
        # Test when num_pops != len(num_samples), we throw error
        with self.assertRaises(ValueError):
            ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2, 0])
        # Test that when steps are not strictly increasing values, we throw error.
        with self.assertRaises(ValueError):
            ddb.coalescence_rate_trajectory(
                steps=np.flip(steps, axis=0), num_samples=[2]
            )
        # Test that when steps are negative, we throw error
        with self.assertRaises(ValueError):
            ddb.coalescence_rate_trajectory(
                steps=np.linspace(-1, 10, 11), num_samples=[2]
            )

    def get_random_example(self):
        """
        A big annoying example.
        """
        random.seed(23)
        Ne = 100
        N = 4
        pop_sizes = [random.uniform(0.01, 10) * Ne for _ in range(N)]
        growth_rates = [random.uniform(-0.01, 0.01) for _ in range(N)]
        migration_matrix = [
            [random.random() * (i != j) for j in range(N)] for i in range(N)
        ]
        sample_sizes = [random.randint(2, 10) for _ in range(N)]
        population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=k, sample_size=n, growth_rate=r
            )
            for k, n, r in zip(pop_sizes, sample_sizes, growth_rates)
        ]
        demographic_events = []
        for i in [0, 1]:
            n = random.uniform(0.01, 10)
            r = 0
            demographic_events.append(
                msprime.PopulationParametersChange(
                    time=100, initial_size=n, growth_rate=r, population_id=i
                )
            )
        for ij in [(0, 1), (2, 3), (0, 3)]:
            demographic_events.append(
                msprime.MigrationRateChange(180, random.random(), matrix_index=ij)
            )
        demographic_events.append(
            msprime.MassMigration(time=200, source=3, dest=0, proportion=0.3)
        )
        for i in [1, 3]:
            n = random.uniform(0.01, 10)
            r = random.uniform(-0.01, 0.01)
            demographic_events.append(
                msprime.PopulationParametersChange(
                    time=210, initial_size=n, growth_rate=r, population_id=i
                )
            )

        ddb = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            migration_matrix=migration_matrix,
        )
        return ddb

    def test_random_example(self):
        ddb = self.get_random_example()
        num_samples = list(range(ddb.num_populations))
        rates, P = ddb.coalescence_rate_trajectory(
            steps=np.linspace(0, 200, 2001), num_samples=num_samples
        )
        self.assertTrue(np.all(rates >= 0))
        self.assertTrue(np.all(P >= 0))
        self.assertTrue(np.all(P <= 1))
        self.assertTrue(np.all(np.diff(P) <= 0))
        coaltime = ddb.mean_coalescence_time(num_samples=num_samples)
        self.assertGreater(coaltime, 0)
        coaltime2 = ddb.mean_coalescence_time(
            num_samples=num_samples, steps=np.linspace(0, 200, 501)
        )
        self.assertLess(abs(coaltime - coaltime2), 2)


class TestMatrixExponential(unittest.TestCase):
    """
    Test cases for the matrix exponential function.
    """

    def verify(self, A):
        E1 = scipy.linalg.expm(A)
        E2 = msprime.demography._matrix_exponential(A)
        self.assertEqual(E1.shape, E2.shape)
        self.assertTrue(np.allclose(E1, E2))

    def test_singleton(self):
        for j in range(10):
            A = np.array([[j]])
            self.verify(A)

    def test_zeros(self):
        for j in range(1, 10):
            A = np.zeros((j, j))
            self.verify(A)

    def test_ones_minus_diagonal(self):
        # If we got to larger values we start getting complex number results.
        # (k x k) matrices of ones, but with (-k) on the diagonal, for k >= 2.
        for j in range(2, 5):
            A = np.ones((j, j))
            A = A - (2 * np.eye(j))
            self.verify(A)

    def test_singleton_against_exp(self):
        # a 1 x 1 matrix consisting of just 0 (compared to exp(0) = 1)
        # a 1 x 1 matrix consisting of just -1 (compared to exp(-1))
        for t in [0, -1]:
            A = msprime.demography._matrix_exponential([[t]])
            B = np.exp(t)
            self.assertEqual(A, B)

    def test_identity_exp(self):
        # (-1) * np.eye(k), compared to exp(-1) * np.eye(k)
        for k in range(2, 5):
            A = msprime.demography._matrix_exponential((-1) * np.eye(k))
            B = np.exp(-1) * np.eye(k)
            self.assertTrue(np.allclose(A, B))


class TestEventTimes(unittest.TestCase):
    """
    Tests that demographic events occur when they should.
    """

    def test_event_at_start_time(self):
        for start_time in [0, 10, 20]:
            ts = msprime.simulate(
                population_configurations=[
                    msprime.PopulationConfiguration(2),
                    msprime.PopulationConfiguration(0),
                ],
                demographic_events=[
                    msprime.MassMigration(time=start_time, source=0, dest=1),
                ],
                random_seed=1,
                start_time=start_time,
                record_migrations=True,
            )
            migrations = list(ts.migrations())
            self.assertEqual(len(migrations), 2)
            self.assertEqual(migrations[0].time, start_time)
            self.assertEqual(migrations[1].time, start_time)
            nodes = list(ts.nodes())
            self.assertEqual(nodes[0].population, 0)
            self.assertEqual(nodes[1].population, 0)
            self.assertEqual(nodes[2].population, 1)

    def test_negative_times(self):
        with self.assertRaises(ValueError):
            msprime.simulate(
                sample_size=10,
                demographic_events=[
                    msprime.PopulationParametersChange(time=-1, initial_size=2)
                ],
            )

    def test_event_before_start_time(self):
        for start_time in [10, 20]:
            for time in [start_time - 1, start_time - 1e-6]:
                ts = msprime.simulate(
                    sample_size=10,
                    start_time=start_time,
                    demographic_events=[
                        msprime.PopulationParametersChange(time=time, initial_size=2)
                    ],
                )
                self.assertEqual(ts.first().num_roots, 1)


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
            random_seed=1,
        )
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
        self.assertEqual(ts.num_populations, 3)

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
            random_seed=1,
        )
        tree = next(ts.trees())
        self.assertEqual(tree.root, 2 * n - 2)
        self.assertGreater(tree.time(tree.root), t)
        for j in range(n // 2):
            self.assertEqual(tree.population(j), 0)
            self.assertEqual(tree.population(n // 2 + j), 1)
            self.assertEqual(ts.get_population(j), 0)
            self.assertEqual(ts.get_population(n // 2 + j), 1)
        self.assertEqual(tree.population(tree.root), 2)

        self.assertTrue(
            np.array_equal(ts.samples(0), np.arange(n // 2, dtype=np.int32))
        )
        self.assertTrue(
            np.array_equal(ts.samples(1), np.arange(n // 2, n, dtype=np.int32))
        )
        self.assertTrue(np.array_equal(ts.samples(2), np.array([], dtype=np.int32)))
        self.assertEqual(ts.num_populations, 3)

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
            random_seed=1,
        )
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
        self.assertEqual(tree.population(tree.get_mrca(0, 2 * (n // 3))), 0)
        # The MRCAs of all the samples within each deme must have
        # occured within that deme
        for k in range(3):
            deme_samples = range(k * (n // 3), (k + 1) * (n // 3))
            for u, v in itertools.combinations(deme_samples, 2):
                mrca_pop = tree.population(tree.get_mrca(u, v))
                self.assertEqual(k, mrca_pop)
        self.assertTrue(
            np.array_equal(ts.samples(0), np.arange(n // 3, dtype=np.int32))
        )
        self.assertTrue(
            np.array_equal(
                ts.samples(1), np.arange(n // 3, 2 * (n // 3), dtype=np.int32)
            )
        )
        self.assertTrue(
            np.array_equal(ts.samples(2), np.arange(2 * (n // 3), n, dtype=np.int32))
        )

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
            random_seed=1,
        )
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
            random_seed=1,
        )
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

    def test_empty_demes_model_changes(self):
        t1 = 1
        t2 = 100
        t3 = 200
        dt = 10
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(0),
            msprime.PopulationConfiguration(0),
            msprime.PopulationConfiguration(1),
        ]
        # We migrate the lineages to the next step by step and intersperse
        # them with model change events
        demographic_events = [
            msprime.MassMigration(time=t1, source=0, dest=1),
            msprime.SimulationModelChange(t1 + dt, "dtwf"),
            msprime.MassMigration(time=t2, source=1, dest=2),
            msprime.SimulationModelChange(t2 + dt, "hudson"),
            msprime.MassMigration(time=t3, source=2, dest=3),
            msprime.SimulationModelChange(t3 + dt, "dtwf"),
        ]
        ts = msprime.simulate(
            Ne=100,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
        )
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
        self.assertEqual(ts.num_populations, 4)

    def test_migration_rate_directionality(self):
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(0),
        ]
        t = 5
        demographic_events = [
            msprime.MigrationRateChange(time=t, rate=1, source=0, dest=2),
            msprime.MigrationRateChange(time=t, rate=1, source=1, dest=2),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
        )
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

    def test_migration_rate_directionality_from_ts(self):
        tables = tskit.TableCollection(1)
        for _ in range(3):
            tables.populations.add_row()
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=1)

        population_configurations = [
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
        ]
        t = 5
        demographic_events = [
            msprime.MigrationRateChange(time=t, rate=1, matrix_index=(0, 2)),
            msprime.MigrationRateChange(time=t, rate=1, matrix_index=(1, 2)),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            from_ts=tables.tree_sequence(),
            start_time=0,
            random_seed=1,
        )
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
        population_configurations = (
            [msprime.PopulationConfiguration(1)]
            + [msprime.PopulationConfiguration(0) for _ in range(num_demes - 2)]
            + [msprime.PopulationConfiguration(1)]
        )
        t = 5
        demographic_events = [
            msprime.MassMigration(time=t, source=0, dest=num_demes - 1),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
        )
        tree = next(ts.trees())
        self.assertEqual(tree.root, 2)
        self.assertGreater(tree.time(2), t)
        self.assertEqual(tree.population(0), 0)
        self.assertEqual(tree.population(1), num_demes - 1)
        self.assertEqual(tree.population(2), num_demes - 1)
        self.assertEqual(ts.node(0).population, 0)
        self.assertEqual(ts.node(1).population, num_demes - 1)

    def test_many_demes_from_ts(self):
        num_demes = 300
        tables = tskit.TableCollection(1)
        for _ in range(num_demes):
            tables.populations.add_row()
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, time=0, population=num_demes - 1
        )
        population_configurations = (
            [msprime.PopulationConfiguration()]
            + [msprime.PopulationConfiguration() for _ in range(num_demes - 2)]
            + [msprime.PopulationConfiguration()]
        )
        t = 5
        demographic_events = [
            msprime.MassMigration(time=t, source=0, dest=num_demes - 1),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            from_ts=tables.tree_sequence(),
            start_time=0,
            random_seed=1,
        )
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
            msprime.MassMigration(time=t4, source=1, dest=0),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
        )
        tree = next(ts.trees())
        self.assertGreater(tree.time(tree.root), t4)
        self.assertEqual(tree.population(tree.root), 0)
        # The parent of all the samples from each deme should be in that deme.
        for pop in range(3):
            parents = [tree.get_parent(u) for u in ts.samples(population=pop)]
            for v in parents:
                self.assertEqual(tree.population(v), pop)
        self.assertEqual(ts.num_populations, 3)


class MigrationRecordsMixin:
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

    def verify_two_pops_single_sample(self, ts, t):
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
            Ne=100,
            model=self.model,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
            record_migrations=True,
        )
        self.verify_two_pops_single_sample(ts, t)

    def test_two_pops_single_sample_from_ts(self):
        tables = tskit.TableCollection(1)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=1)
        for _ in range(3):
            tables.populations.add_row()
        population_configurations = [
            msprime.PopulationConfiguration() for _ in range(3)
        ]
        t = 5
        demographic_events = [
            msprime.MassMigration(time=t, source=0, dest=2),
            msprime.MassMigration(time=t, source=1, dest=2),
        ]
        ts = msprime.simulate(
            Ne=100,
            model=self.model,
            from_ts=tables.tree_sequence(),
            start_time=0,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            record_migrations=True,
            random_seed=1,
        )
        self.verify_two_pops_single_sample(ts, t)

    def verify_two_pops_asymmetric_migrations(self, ts):
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
        # All lineages should end up in population 1
        for tree in ts.trees():
            node = ts.node(tree.root)
            self.assertEqual(node.population, 0)

    def test_two_pops_asymmetric_migrations(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        ts = msprime.simulate(
            model=self.model,
            population_configurations=population_configurations,
            migration_matrix=[[0, 0], [1, 0]],
            # Can migrate from 1 to 0 but not vice-versa
            random_seed=1,
            record_migrations=True,
        )
        self.verify_two_pops_asymmetric_migrations(ts)

    def test_two_pops_asymmetric_migrations_from_ts(self):
        tables = tskit.TableCollection(1)
        for _ in range(10):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
        for _ in range(10):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=1)
        tables.populations.add_row()
        tables.populations.add_row()

        population_configurations = [
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
        ]
        ts = msprime.simulate(
            model=self.model,
            population_configurations=population_configurations,
            migration_matrix=[[0, 0], [1, 0]],
            # Can migrate from 1 to 0 but not vice-versa
            random_seed=1,
            record_migrations=True,
            from_ts=tables.tree_sequence(),
            start_time=0,
        )
        self.verify_two_pops_asymmetric_migrations(ts)

    def test_two_pops_asymmetric_migrations_recombination(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        ts = msprime.simulate(
            model=self.model,
            recombination_rate=1,
            population_configurations=population_configurations,
            migration_matrix=[[0, 0], [1, 0]],
            # Can migrate from 1 to 0 but not vice-versa
            random_seed=1,
            record_migrations=True,
        )
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
            model=self.model,
            recombination_rate=10,
            population_configurations=population_configurations,
            demographic_events=[
                msprime.MassMigration(time=20, source=1, dest=0, proportion=1)
            ],
            random_seed=1,
            record_migrations=True,
        )
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


class TestMigrationRecordsHudson(unittest.TestCase, MigrationRecordsMixin):
    model = "hudson"


class TestMigrationRecordsSmc(unittest.TestCase, MigrationRecordsMixin):
    model = "smc"


class TestMigrationRecordsSmcPrime(unittest.TestCase, MigrationRecordsMixin):
    model = "smc_prime"


class TestMigrationRecordsDtwf(unittest.TestCase, MigrationRecordsMixin):
    model = "dtwf"


class TestFullArgMigration(unittest.TestCase):
    """
    Tests for migration with the full ARG.
    """

    def verify_two_pops_full_arg(self, ts):
        migrations = ts.tables.migrations
        edges = ts.tables.edges
        nodes = ts.tables.nodes
        for mig in migrations:
            self.assertEqual(nodes[mig.node].flags, msprime.NODE_IS_MIG_EVENT)
            self.assertEqual(nodes[mig.node].time, mig.time)
            self.assertEqual(nodes[mig.node].population, mig.dest)
            e1 = np.where(edges.parent == mig.node)
            e2 = np.where(edges.left == mig.left)
            e = np.intersect1d(e1[0], e2[0])
            self.assertEqual(len(e), 1)
            e = np.asscalar(e)
            self.assertEqual(edges[e].right, mig.right)
            self.assertEqual(nodes[edges[e].child].population, mig.source)
        for edge in edges:
            if nodes[edge.parent].flags == msprime.NODE_IS_MIG_EVENT:
                m1 = np.where(migrations.node == edge.parent)
                m2 = np.where(migrations.left == edge.left)
                m = np.intersect1d(m1[0], m2[0])
                self.assertEqual(len(m), 1)
                m = np.asscalar(m)
                self.assertEqual(migrations[m].right, edge.right)
                self.assertEqual(migrations[m].time, nodes[edge.parent].time)
                self.assertEqual(migrations[m].source, nodes[edge.child].population)
                self.assertEqual(migrations[m].dest, nodes[edge.parent].population)

    def test_full_arg_migration(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            migration_matrix=[[0, 1], [1, 0]],
            random_seed=1,
            recombination_rate=0.1,
            record_migrations=True,
            record_full_arg=True,
        )
        self.verify_two_pops_full_arg(ts)

    def test_full_arg_migration_smc(self):
        for model in ["smc", "smc_prime"]:

            population_configurations = [
                msprime.PopulationConfiguration(10),
                msprime.PopulationConfiguration(10),
            ]
            ts = msprime.simulate(
                population_configurations=population_configurations,
                migration_matrix=[[0, 1], [1, 0]],
                random_seed=101,
                recombination_rate=0.1,
                model=model,
                record_migrations=True,
                record_full_arg=True,
            )
            self.verify_two_pops_full_arg(ts)


class TimeUnitsMixin:
    """
    Tests for time conversion between generations and coalescent
    units.
    """

    def test_coalescence_after_size_change(self):
        Ne = 2000
        # Migrations and bottleneck occured 100 generations ago.
        g = 100
        population_configurations = [
            msprime.PopulationConfiguration(1),
            msprime.PopulationConfiguration(1),
        ]
        # At this time, we migrate the lineage in 1 to 0, and
        # have a very strong bottleneck, resulting in almost instant
        # coalescence.
        demographic_events = [
            msprime.MassMigration(time=g, source=1, dest=0),
            msprime.PopulationParametersChange(time=g, initial_size=1),
        ]
        reps = msprime.simulate(
            Ne=Ne,
            model=self.model,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
            num_replicates=10,
        )
        for ts in reps:
            tree = ts.first()
            u = tree.get_mrca(0, 1)
            self.assertEqual(u, 2)
            self.assertLessEqual(g, tree.time(u))


class TestTimeUnitsHudson(unittest.TestCase, TimeUnitsMixin):
    model = "hudson"

    # We don't run this test in the DTWF case because extreme growth rates like
    # this are problematic.
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
            model=self.model,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
            num_replicates=10,
            record_migrations=True,
        )
        for ts in reps:
            migrations = list(ts.migrations())
            self.assertGreater(len(migrations), 0)
            for mr in migrations:
                self.assertAlmostEqual(g, mr.time)
            tree = next(ts.trees())
            u = tree.get_mrca(0, 1)
            self.assertEqual(u, 2)
            self.assertAlmostEqual(g, tree.time(u), places=1)

    # Bottlenecks are not supported in the DTWF.
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
            msprime.InstantaneousBottleneck(time=t, population=0, strength=100),
        ]
        reps = msprime.simulate(
            Ne=Ne,
            model=self.model,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1,
            num_replicates=10,
        )
        for ts in reps:
            tree = next(ts.trees())
            self.assertAlmostEqual(t, tree.time(tree.root), places=5)


class TestTimeUnitsWrightFisher(unittest.TestCase, TimeUnitsMixin):
    model = "dtwf"


class TestLowLevelConversions(unittest.TestCase):
    """
    Checks that we convert to the correct low-level values when we
    do the rescalings from generations.
    """

    def test_population_configuration_defaults(self):
        conf = msprime.PopulationConfiguration()
        self.assertIsNone(conf.sample_size)
        d = conf.asdict()
        dp = {
            "initial_size": None,
            "growth_rate": 0,
            "metadata": None,
            "sample_size": None,
        }
        self.assertEqual(d, dp)

    def test_population_configuration_initial_size(self):
        for initial_size in [1, 10, 1000]:
            conf = msprime.PopulationConfiguration(initial_size=initial_size)
            self.assertIsNone(conf.sample_size)
            d = conf.asdict()
            dp = {
                "initial_size": initial_size,
                "growth_rate": 0,
                "metadata": None,
                "sample_size": None,
            }
            self.assertEqual(d, dp)

    def test_population_configuration_growth_rate(self):
        sample_size = 8
        for growth_rate in [1, 10, -10]:
            conf = msprime.PopulationConfiguration(sample_size, growth_rate=growth_rate)
            self.assertEqual(conf.sample_size, sample_size)
            d = conf.asdict()
            dp = {
                "initial_size": None,
                "growth_rate": growth_rate,
                "metadata": None,
                "sample_size": sample_size,
            }
            self.assertEqual(d, dp)

    def test_population_parameters_change_time(self):
        for Ne in [1, 10, 1000]:
            for g in [0.1, 1, 100, 1e6]:
                event = msprime.PopulationParametersChange(time=g, initial_size=Ne)
                d = event.get_ll_representation()
                dp = {
                    "time": g,
                    "population": -1,
                    "type": "population_parameters_change",
                    "initial_size": Ne,
                }
                self.assertEqual(d, dp)

    def test_population_parameters_change_initial_size(self):
        g = 100
        for initial_size in [0.01, 1, 100, 1e6]:
            event = msprime.PopulationParametersChange(
                time=g, initial_size=initial_size
            )
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "population": -1,
                "type": "population_parameters_change",
                "initial_size": initial_size,
            }
            self.assertEqual(d, dp)

    def test_population_parameters_change_growth_rate(self):
        g = 100
        for growth_rate in [0.01, 1, 100, 1e6]:
            event = msprime.PopulationParametersChange(time=g, growth_rate=growth_rate)
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "population": -1,
                "type": "population_parameters_change",
                "growth_rate": growth_rate,
            }
            self.assertEqual(d, dp)

    def test_population_parameters_change_population(self):
        g = 100
        Ne = 10
        for population in range(3):
            event = msprime.PopulationParametersChange(
                time=g, initial_size=Ne, population=population
            )
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "population": population,
                "type": "population_parameters_change",
                "initial_size": Ne,
            }
            self.assertEqual(d, dp)

    def test_migration_rate_change_time(self):
        for g in [0.1, 1, 100, 1e6]:
            event = msprime.MigrationRateChange(time=g, rate=0)
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "type": "migration_rate_change",
                "migration_rate": 0,
                "source": -1,
                "dest": -1,
            }
            self.assertEqual(d, dp)

    def test_migration_rate_change_matrix_index(self):
        g = 51
        for N in range(1, 5):
            for index in itertools.permutations(range(N), 2):
                event = msprime.MigrationRateChange(
                    time=g, rate=0, source=index[0], dest=index[1]
                )
                d = event.get_ll_representation()
                dp = {
                    "time": g,
                    "type": "migration_rate_change",
                    "migration_rate": 0,
                    "source": index[0],
                    "dest": index[1],
                }
                self.assertEqual(d, dp)

                # Check the deprecated form
                event = msprime.MigrationRateChange(time=g, rate=0, matrix_index=index)
                d = event.get_ll_representation()
                self.assertEqual(d, dp)

    def test_migration_rate_change_rate(self):
        g = 1234
        for rate in [0, 1e-6, 10, 1e6]:
            event = msprime.MigrationRateChange(time=g, rate=rate)
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "type": "migration_rate_change",
                "migration_rate": rate,
                "source": -1,
                "dest": -1,
            }
            self.assertEqual(d, dp)

    def test_mass_migration_time(self):
        for g in [0.1, 1, 100, 1e6]:
            event = msprime.MassMigration(time=g, source=0, dest=1)
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "type": "mass_migration",
                "source": 0,
                "dest": 1,
                "proportion": 1,
            }
            self.assertEqual(d, dp)

    def test_mass_migration_source_dest(self):
        g = 51
        for source, dest in itertools.permutations(range(4), 2):
            event = msprime.MassMigration(time=g, source=source, dest=dest)
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "type": "mass_migration",
                "source": source,
                "dest": dest,
                "proportion": 1,
            }
            self.assertEqual(d, dp)

    def test_mass_migration_proportion(self):
        g = 51
        for p in [0, 1e-6, 0.4, 1]:
            event = msprime.MassMigration(time=g, source=0, dest=1, proportion=p)
            d = event.get_ll_representation()
            dp = {
                "time": g,
                "type": "mass_migration",
                "source": 0,
                "dest": 1,
                "proportion": p,
            }
            self.assertEqual(d, dp)

    def test_migration_matrix(self):
        m = [[0, 1, 2], [3, 0, 4], [5, 6, 0]]
        sim = ancestry._parse_simulate(
            population_configurations=[
                msprime.PopulationConfiguration(1),
                msprime.PopulationConfiguration(1),
                msprime.PopulationConfiguration(1),
            ],
            migration_matrix=m,
        )
        np.testing.assert_array_equal(sim.demography.migration_matrix, m)

    def test_instantaneous_bottleneck(self):
        g = 51
        for population in [0, 1, 5]:
            for strength in [0, 100, 1000, 1e9]:
                event = msprime.InstantaneousBottleneck(
                    time=g, population=population, strength=strength
                )
                d = event.get_ll_representation()
                dp = {
                    "time": g,
                    "type": "instantaneous_bottleneck",
                    "population": population,
                    "strength": strength,
                }
                self.assertEqual(d, dp)


class HistoricalSamplingMixin:
    """
    Tests to make sure historical sampling works correctly.
    """

    def test_two_samples(self):
        N = 100
        sampling_time = 1.01 * N
        for recombination_rate in [0, 1]:
            ts = msprime.simulate(
                Ne=N,
                model=self.model,
                recombination_map=msprime.RecombinationMap.uniform_map(
                    length=1, rate=recombination_rate,
                ),
                samples=[msprime.Sample(0, 0), msprime.Sample(0, sampling_time)],
                random_seed=3,
                discrete_genome=True,
            )
            for t in ts.trees():
                self.assertEqual(t.get_time(0), 0)
                self.assertEqual(t.get_time(1), sampling_time)
                self.assertEqual(t.get_parent(0), t.get_parent(1))
                self.assertEqual(t.get_parent(1), t.get_parent(0))
                self.assertGreater(t.get_time(t.get_parent(0)), sampling_time)

    def test_two_samples_start_time(self):
        N = 10
        sampling_time = 10.01 * N
        for start_time in [0, sampling_time / 2, sampling_time, 10000 * sampling_time]:
            ts = msprime.simulate(
                Ne=N,
                start_time=start_time,
                model=self.model,
                random_seed=3,
                samples=[msprime.Sample(0, 0), msprime.Sample(0, sampling_time)],
            )
            nodes = list(ts.nodes())
            self.assertEqual(ts.num_nodes, 3)
            self.assertEqual(nodes[0].time, 0)
            self.assertEqual(nodes[1].time, sampling_time)
            self.assertGreater(nodes[2].time, sampling_time)
            self.assertGreater(nodes[2].time, start_time)

    def test_different_times(self):
        N = 50
        st1 = 1.01 * N
        st2 = 2.01 * N
        st3 = 3.01 * N
        ts = msprime.simulate(
            model=self.model,
            Ne=N,
            samples=[
                msprime.Sample(0, 0),
                msprime.Sample(0, st1),
                msprime.Sample(0, st2),
                msprime.Sample(0, st3),
            ],
        )
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
        N = 100
        sampling_time = N * 100.01
        n = 5
        samples = [msprime.Sample(0, sampling_time) for j in range(n - 1)] + [
            msprime.Sample(0, 0)
        ]
        ts = msprime.simulate(Ne=N, samples=samples, model=self.model, random_seed=4)
        time = [node.time for node in ts.nodes()]
        for j in range(n - 1):
            self.assertEqual(time[j], sampling_time)
        self.assertEqual(time[n - 1], 0)
        # Allow it to be within 10 coalescent time units.
        self.assertLess(time[-1], sampling_time + 10 * N)

    def test_start_time_invariance(self):
        for N in [10, 100, 128]:
            offset = None
            # The difference between the start time and the coalescence
            # should be invariant.
            for start_time in [0, 10, 20, 50]:
                ts = msprime.simulate(
                    2, Ne=N, start_time=start_time, model=self.model, random_seed=2
                )
                time = [node.time for node in ts.nodes()]
                self.assertEqual(time[0], 0)
                self.assertEqual(time[1], 0)
                self.assertGreater(time[2], start_time)
                if offset is None:
                    offset = time[2] - start_time
                else:
                    self.assertAlmostEqual(offset, time[2] - start_time)

    def test_negative_start_time(self):
        ts = msprime.simulate(2, Ne=10, start_time=-1, model=self.model, random_seed=2)
        tables = ts.tables
        self.assertEqual(tables.nodes[0].time, 0)
        self.assertEqual(tables.nodes[1].time, 0)
        self.assertEqual(len(tables.edges), 2)

    def test_start_time_before_sample_time(self):
        # If all samples are > than the start time, it doesn't affect
        # the simulation.
        samples = [msprime.Sample(population=0, time=10)] * 2
        for start_time in [-100, 0, 9, 9.999]:
            ts = msprime.simulate(
                samples=samples,
                Ne=10,
                start_time=start_time,
                model=self.model,
                random_seed=2,
            )
            tables = ts.tables
            self.assertEqual(tables.nodes[0].time, 10)
            self.assertEqual(tables.nodes[1].time, 10)
            self.assertEqual(len(tables.edges), 2)

    def test_two_samples_mass_migration(self):
        N = 200
        sampling_time = 2.01 * N
        migration_time = 4.33 * N
        ts = msprime.simulate(
            model=self.model,
            random_seed=10,
            Ne=N,
            samples=[msprime.Sample(0, 0), msprime.Sample(1, sampling_time)],
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration(),
            ],
            demographic_events=[
                msprime.MassMigration(time=migration_time, source=1, dest=0)
            ],
        )
        t = next(ts.trees())
        self.assertEqual(t.get_time(0), 0)
        self.assertEqual(t.get_time(1), sampling_time)
        self.assertGreaterEqual(t.get_time(2), migration_time)
        self.assertEqual(t.get_population(0), 0)
        self.assertEqual(t.get_population(1), 1)
        self.assertEqual(t.get_population(2), 0)

    def test_events_before_sampling(self):
        # Demographic events that are scheduled for before sampling time
        # should go ahead and be applied.
        ts = msprime.simulate(
            model=self.model,
            samples=[msprime.Sample(0, time=100), msprime.Sample(1, time=100)],
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=10),
                msprime.PopulationConfiguration(initial_size=10),
            ],
            # Without this migration rate change, we would never coalesce
            demographic_events=[msprime.MigrationRateChange(99, rate=0.1)],
            record_migrations=True,
            random_seed=2,
        )
        tables = ts.tables
        self.assertEqual(len(tables.edges), 2)
        self.assertGreater(len(tables.migrations), 1)

    def test_interleaved_migrations(self):
        N = 100
        t1 = 1.5 * N
        t2 = 10.5 * N
        t3 = 50.5 * N
        ts = msprime.simulate(
            model=self.model,
            Ne=N,
            samples=[
                msprime.Sample(0, 0),
                msprime.Sample(1, t1),
                msprime.Sample(2, t2),
                msprime.Sample(3, t3),
            ],
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration(),
            ],
            demographic_events=[
                msprime.MassMigration(time=t1, source=0, dest=1),
                msprime.MassMigration(time=t2, source=1, dest=2),
                msprime.MassMigration(time=t3, source=2, dest=3),
            ],
            random_seed=2,
        )
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


class TestHistoricalSamplingHudson(unittest.TestCase, HistoricalSamplingMixin):
    model = "hudson"

    def test_sampling_time_invariance(self):
        for N in [10, 100, 128]:
            offset = None
            # The difference between the sampling time and the coalescence
            # should be invariant.
            for sampling_time in [0, 10, 20, 50]:
                samples = [msprime.Sample(0, sampling_time), msprime.Sample(0, 0)]
                ts = msprime.simulate(
                    Ne=N, samples=samples, model=self.model, random_seed=2
                )
                time = [node.time for node in ts.nodes()]
                self.assertEqual(time[0], sampling_time)
                self.assertEqual(time[1], 0)
                if offset is None:
                    offset = time[2] - sampling_time
                else:
                    self.assertAlmostEqual(offset, time[2] - sampling_time)


class TestHistoricalSamplingWrightFisher(unittest.TestCase, HistoricalSamplingMixin):
    model = "dtwf"

    def test_simultaneous_historical_samples(self):
        N = 10
        samples = [msprime.Sample(0, 0), msprime.Sample(0, 1.1), msprime.Sample(0, 1.2)]
        ts = msprime.simulate(Ne=N, samples=samples, model=self.model, random_seed=2)
        time = [node.time for node in ts.nodes()]
        self.assertEqual(time[0], 0)
        self.assertEqual(time[1], 1.1)
        self.assertEqual(time[2], 1.2)


class EndTimeMixin:
    """
    Tests for the max_time parameter.
    """

    def verify_empty_tree_sequence(self, n, ts):
        self.assertEqual(ts.num_edges, 0)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, n)
        self.assertEqual(ts.num_samples, n)
        tree = ts.first()
        self.assertEqual(tree.num_roots, n)

    def verify_incomplete_tree_sequence(self, n, max_time, ts):
        self.assertEqual(ts.num_samples, n)
        time = ts.tables.nodes.time
        for tree in ts.trees():
            # Every sample with time <= max_time will end on a path
            # with time == max_time
            for u in tree.samples():
                if time[u] <= max_time:
                    while tree.parent(u) != tskit.NULL:
                        u = tree.parent(u)
                    self.assertEqual(ts.node(u).time, max_time)
                else:
                    self.assertEqual(tree.parent(u), tskit.NULL)
        max_roots = max(tree.num_roots for tree in ts.trees())
        self.assertGreater(max_roots, 1)

    def test_zero_time(self):
        n = 10
        for n in [2, 10, 100]:
            ts = msprime.simulate(n, end_time=0, model=self.model)
            self.verify_empty_tree_sequence(n, ts)

    def test_negative(self):
        with self.assertRaises(ValueError):
            msprime.simulate(3, end_time=-1, model=self.model)

    def test_large_time(self):
        seed = 1
        ts1 = msprime.simulate(
            10, Ne=100, end_time=1e10, model=self.model, random_seed=seed
        )
        ts2 = msprime.simulate(10, Ne=100, model=self.model, random_seed=seed)
        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()
        tables1.provenances.clear()
        tables2.provenances.clear()
        self.assertEqual(tables1, tables2)

    def test_small_time(self):
        n = 100
        max_time = 100
        ts = msprime.simulate(
            n, Ne=1000, end_time=max_time, model=self.model, random_seed=10
        )
        self.verify_incomplete_tree_sequence(n, max_time, ts)

    def test_migrations(self):
        n = 10
        max_time = 100
        ts = msprime.simulate(
            Ne=1000,
            end_time=max_time,
            model=self.model,
            random_seed=10,
            record_migrations=True,
            migration_matrix=[[0, 1], [1, 0]],
            population_configurations=[
                msprime.PopulationConfiguration(n),
                msprime.PopulationConfiguration(n),
            ],
        )
        self.verify_incomplete_tree_sequence(2 * n, max_time, ts)
        self.assertGreater(ts.num_migrations, 0)
        self.assertTrue(np.all(ts.tables.migrations.time < max_time))

    def test_ancient_samples(self):
        n = 40
        samples = [msprime.Sample(time=j, population=0) for j in range(n)]
        max_time = 20
        ts = msprime.simulate(
            samples=samples, Ne=10, end_time=max_time, model=self.model, random_seed=100
        )
        self.verify_incomplete_tree_sequence(n, max_time, ts)
        nodes = ts.tables.nodes
        self.assertTrue(np.array_equal(nodes.time[:n], np.arange(n)))
        self.assertGreater(len(nodes), n)
        tree = ts.first()
        self.assertGreater(tree.num_roots, 1)

    def test_all_ancient_samples(self):
        n = 40
        samples = [msprime.Sample(time=j + 1, population=0) for j in range(n)]
        max_time = 20
        ts = msprime.simulate(
            samples=samples, Ne=10, end_time=max_time, model=self.model, random_seed=100
        )
        self.verify_incomplete_tree_sequence(n, max_time, ts)
        nodes = ts.tables.nodes
        self.assertTrue(np.array_equal(nodes.time[:n], np.arange(n) + 1))
        self.assertGreater(len(nodes), n)
        tree = ts.first()
        self.assertGreater(tree.num_roots, 1)

    def test_ancient_samples_equal_time(self):
        max_time = 10
        samples = [
            msprime.Sample(time=0, population=0),
            msprime.Sample(time=max_time, population=0),
        ]
        ts = msprime.simulate(
            samples=samples,
            Ne=10,
            end_time=max_time,
            model=self.model,
            random_seed=1000,
        )
        self.verify_incomplete_tree_sequence(2, max_time, ts)
        self.assertEqual(ts.num_nodes, 3)
        self.assertEqual(ts.num_edges, 1)

    def test_demographic_events(self):
        max_time = 100
        n = 20
        ts = msprime.simulate(
            n,
            Ne=1000,
            end_time=max_time,
            model=self.model,
            demographic_events=[
                msprime.SimpleBottleneck(time=max_time, population=0, proportion=1)
            ],
            random_seed=1000,
        )
        self.verify_incomplete_tree_sequence(n, max_time, ts)


class TestEndTimeHudson(unittest.TestCase, EndTimeMixin):
    model = "hudson"


class TestEndTimeWrightFisher(unittest.TestCase, EndTimeMixin):
    model = "dtwf"


class TestEventsBetweenGenerationsWrightFisher(unittest.TestCase):
    """
    Tests that events occuring between generations in the DTWF are
    handled correctly.
    """

    def test_4_populations(self):
        migration_matrix = np.zeros((4, 4))
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=10, initial_size=10, growth_rate=0
            ),
            msprime.PopulationConfiguration(
                sample_size=10, initial_size=10, growth_rate=0
            ),
            msprime.PopulationConfiguration(
                sample_size=0, initial_size=10, growth_rate=0
            ),
            msprime.PopulationConfiguration(
                sample_size=0, initial_size=10, growth_rate=0
            ),
        ]
        demographic_events = [
            msprime.PopulationParametersChange(population=1, time=0.1, initial_size=5),
            msprime.PopulationParametersChange(population=0, time=0.2, initial_size=5),
            msprime.MassMigration(time=1.1, source=0, dest=2),
            msprime.MassMigration(time=1.2, source=1, dest=3),
            msprime.MigrationRateChange(time=2.1, rate=0.3, source=2, dest=3),
            msprime.MigrationRateChange(time=2.2, rate=0.3, source=3, dest=2),
        ]
        ts = msprime.simulate(
            migration_matrix=migration_matrix,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=2,
            model="dtwf",
        )
        for node in ts.nodes():
            self.assertEqual(node.time, int(node.time))


class TestPopulationMetadata(unittest.TestCase):
    """
    Tests for the metadata behaviour on populations.
    """

    def test_simple_case(self):
        md = {"x": "y"}
        ts = msprime.simulate(
            population_configurations=[msprime.PopulationConfiguration(2, metadata=md)],
            random_seed=1,
        )
        self.assertEqual(ts.num_populations, 1)
        pop = ts.population(0)
        self.assertEqual(md, json.loads(pop.metadata.decode()))

    def test_default(self):
        ts = msprime.simulate(
            population_configurations=[msprime.PopulationConfiguration(2)],
            random_seed=1,
        )
        self.assertEqual(ts.num_populations, 1)
        pop = ts.population(0)
        self.assertEqual(b"", pop.metadata)

        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(2, metadata=None)
            ],
            random_seed=1,
        )
        self.assertEqual(ts.num_populations, 1)
        pop = ts.population(0)
        self.assertEqual(b"", pop.metadata)

    def test_errors(self):
        for bad_metadata in [b"asdf", Exception]:
            popconf = msprime.PopulationConfiguration(2, metadata=bad_metadata)
            pop = msprime.Population.from_old_style(popconf)
            with self.assertRaises(TypeError):
                pop.temporary_hack_for_encoding_old_style_metadata()

    def test_multi_population(self):
        for num_pops in range(1, 10):
            pop_configs = [
                msprime.PopulationConfiguration(2, metadata={"x": "x" * j})
                for j in range(num_pops)
            ]
            ts = msprime.simulate(
                population_configurations=pop_configs, random_seed=1, end_time=1
            )
            self.assertEqual(ts.num_populations, num_pops)
            for j in range(num_pops):
                pop = ts.population(j)
                self.assertEqual(
                    pop_configs[j].metadata, json.loads(pop.metadata.decode())
                )


class TestCensusEvent(unittest.TestCase):
    """
    Tests of the census demographic event.
    """

    def verify(self, ts, census_time):
        """
        Verifies that a census event has been added correctly.
        """
        census_ids = np.where(ts.tables.nodes.flags == msprime.NODE_IS_CEN_EVENT)[0]
        for u in census_ids:
            self.assertEqual(ts.tables.nodes.time[u], census_time)
        self.assertGreater(len(census_ids), 1)
        # Check that all samples have a census ancestor on each tree.
        for tree in ts.trees():
            leaves = []
            census_nodes = [u for u in census_ids if u in list(tree.nodes())]
            for node in census_nodes:
                self.assertEqual(len(tree.children(node)), 1)
                le = list(tree.leaves(node))
                leaves += le
            leaves.sort()
            self.assertEqual(leaves, [i for i in range(0, ts.num_samples)])

    def test_simple_case(self):
        census_time = 0.5
        ts = msprime.simulate(
            sample_size=5,
            random_seed=1,
            demographic_events=[msprime.CensusEvent(time=census_time)],
        )
        self.verify(ts, census_time)

    def test_multiple_trees(self):
        census_time = 0.5
        ts = msprime.simulate(
            sample_size=5,
            random_seed=1,
            recombination_rate=0.4,
            demographic_events=[msprime.CensusEvent(time=census_time)],
        )
        self.verify(ts, census_time)

    def test_population_IDs(self):
        census_time = 100
        pop = msprime.PopulationConfiguration(sample_size=8, initial_size=500)
        mig_rate_change = msprime.MigrationRateChange(time=200, rate=0.05)
        ts = msprime.simulate(
            population_configurations=[pop, pop],
            length=1000,
            demographic_events=[msprime.CensusEvent(time=census_time), mig_rate_change],
            recombination_rate=1e-5,
            random_seed=142,
        )
        self.verify(ts, census_time)
        # Since there is no migration between generations 0 - 200, the census nodes
        # should have the same population label as their children in the trees.
        census_ids = np.where(ts.tables.nodes.flags == msprime.NODE_IS_CEN_EVENT)[0]
        nodes = ts.tables.nodes
        for row in ts.tables.edges:
            if row.parent in census_ids:
                self.assertEqual(
                    nodes.population[row.parent], nodes.population[row.child]
                )

    def test_census_at_existing_node_time(self):
        with self.assertRaises(_msprime.LibraryError):
            msprime.simulate(
                sample_size=2,
                random_seed=3,
                demographic_events=[msprime.CensusEvent(time=0)],
            )

    def test_migration_time_equals_census_time(self):
        census_time = 100
        pop = msprime.PopulationConfiguration(sample_size=8, initial_size=500)
        mig_rate_change = msprime.MigrationRateChange(time=census_time, rate=0.05)
        # If the census is before the migration in the list, census nodes should have
        # the same population as their children.
        ts = msprime.simulate(
            population_configurations=[pop, pop],
            length=1000,
            demographic_events=[msprime.CensusEvent(time=census_time), mig_rate_change],
            recombination_rate=1e-5,
            random_seed=142,
        )
        self.verify(ts, census_time)
        census_ids = np.where(ts.tables.nodes.flags == msprime.NODE_IS_CEN_EVENT)[0]
        nodes = ts.tables.nodes
        for row in ts.tables.edges:
            if row.parent in census_ids:
                self.assertEqual(
                    nodes.population[row.parent], nodes.population[row.child]
                )
        # If the census is after the migration in the list, census nodes should have
        # the same population as their parents.
        divergence = msprime.MassMigration(
            time=census_time, source=1, dest=0, proportion=1
        )
        ts = msprime.simulate(
            population_configurations=[pop, pop],
            length=1000,
            demographic_events=[divergence, msprime.CensusEvent(time=census_time)],
            recombination_rate=1e-5,
            random_seed=12,
        )
        self.verify(ts, census_time)
        census_ids = np.where(ts.tables.nodes.flags == msprime.NODE_IS_CEN_EVENT)[0]
        nodes = ts.tables.nodes
        for row in ts.tables.edges:
            if row.child in census_ids:
                self.assertEqual(
                    nodes.population[row.parent], nodes.population[row.child]
                )

    def test_no_census_nodes_above_root_nodes(self):
        ts = msprime.simulate(sample_size=2, random_seed=525)
        self.assertTrue(all(ts.tables.nodes.flags) != msprime.NODE_IS_CEN_EVENT)
        tsc = msprime.simulate(
            sample_size=2,
            random_seed=525,
            demographic_events=[msprime.CensusEvent(time=2000)],
        )
        self.assertEqual(ts.tables.nodes, tsc.tables.nodes)


class TestPossibleLineages(unittest.TestCase):
    """
    Tests for checking where lineages are possible within the demography debugger.
    """

    def test_possible_lineages_no_migration(self):
        samples = [
            msprime.Sample(time=0, population=0),
            msprime.Sample(time=0, population=1),
            msprime.Sample(time=0, population=2),
            msprime.Sample(time=0, population=3),
        ]
        dem_events = [
            msprime.MassMigration(time=50, source=3, destination=2),
            msprime.MassMigration(time=100, source=2, destination=1),
            msprime.MassMigration(time=150, source=1, destination=0),
        ]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events, population_configurations=pop_config
        )
        lineages = dd.possible_lineage_locations(samples=samples)
        self.assertTrue(len(lineages) == 4)
        self.assertTrue(np.all(lineages[(0, 50)] == [1, 1, 1, 1]))
        self.assertTrue(np.all(lineages[(50, 100)] == [1, 1, 1, 0]))
        self.assertTrue(np.all(lineages[(100, 150)] == [1, 1, 0, 0]))
        self.assertTrue(np.all(lineages[(150, np.inf)] == [1, 0, 0, 0]))

    def test_possible_lineages_with_migration(self):
        # draw sample from first population, which has migrants from second
        samples = [msprime.Sample(time=0, population=0)]
        mig_mat = [[0, 0.1], [0.1, 0]]
        dem_events = [msprime.MassMigration(time=50, source=1, destination=0)]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        # suppress the migrations-after-merging warning
        with warnings.catch_warnings(record=True):
            dd = msprime.DemographyDebugger(
                demographic_events=dem_events,
                population_configurations=pop_config,
                migration_matrix=mig_mat,
            )
        lineages = dd.possible_lineage_locations(samples=samples)
        self.assertTrue(len(lineages) == 1)
        self.assertTrue(np.all(lineages[(0, np.inf)] == [1, 1]))

    def test_possible_lineages_ancient_samples(self):
        samples = [
            msprime.Sample(time=0, population=0),
            msprime.Sample(time=0, population=1),
            msprime.Sample(time=100, population=1),
        ]
        dem_events = [
            msprime.MassMigration(time=50, source=1, destination=0),
            msprime.MassMigration(time=150, source=1, destination=0),
        ]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events, population_configurations=pop_config
        )
        lineages = dd.possible_lineage_locations(samples=samples)
        self.assertTrue(len(lineages) == 4)
        self.assertTrue(np.all(lineages[(0, 50)] == [1, 1]))
        self.assertTrue(np.all(lineages[(50, 100)] == [1, 0]))
        self.assertTrue(np.all(lineages[(100, 150)] == [1, 1]))
        self.assertTrue(np.all(lineages[(150, np.inf)] == [1, 0]))

    def test_possible_lineages_complex_history(self):
        samples = [
            msprime.Sample(time=0, population=0),
            msprime.Sample(time=0, population=2),
            msprime.Sample(time=0, population=3),
            msprime.Sample(time=500, population=4),
        ]
        mig_mat = [
            [0, 0, 0, 0, 0],
            [0, 0, 1e-5, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ]
        dem_events = [
            msprime.MassMigration(time=100, source=0, destination=1, proportion=0.1),
            msprime.MassMigration(time=200, source=3, destination=2),
            msprime.MigrationRateChange(time=200, rate=0),
            msprime.MassMigration(time=300, source=1, destination=0),
            msprime.MassMigration(time=400, source=2, destination=4, proportion=0.1),
            msprime.MassMigration(time=600, source=2, destination=0),
            msprime.MassMigration(time=700, source=4, destination=0),
        ]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events,
            population_configurations=pop_config,
            migration_matrix=mig_mat,
        )
        lineages = dd.possible_lineage_locations(samples=samples)
        self.assertTrue(len(lineages) == 7)
        self.assertTrue(np.all(lineages[(0, 100)] == [1, 0, 1, 1, 0]))
        self.assertTrue(np.all(lineages[(100, 200)] == [1, 1, 1, 1, 0]))
        self.assertTrue(np.all(lineages[(200, 300)] == [1, 1, 1, 0, 0]))
        self.assertTrue(np.all(lineages[(300, 400)] == [1, 0, 1, 0, 0]))
        self.assertTrue(np.all(lineages[(400, 600)] == [1, 0, 1, 0, 1]))
        self.assertTrue(np.all(lineages[(600, 700)] == [1, 0, 0, 0, 1]))
        self.assertTrue(np.all(lineages[(700, np.inf)] == [1, 0, 0, 0, 0]))


class TestLineageProbabilities(unittest.TestCase):
    """
    Tests for checking where lineages are possible within the demography debugger.
    """

    def two_pop_example(self, a, b):
        # Just migration between two populations:
        #  the transition probability matrix for the lineage random walk is:
        def f(t):
            return np.array(
                [
                    [
                        np.exp(-(a + b) * t) + (1 - np.exp(-(a + b) * t)) * b / (a + b),
                        (1 - np.exp(-(a + b) * t)) * a / (a + b),
                    ],
                    [
                        (1 - np.exp(-(a + b) * t)) * b / (a + b),
                        np.exp(-(a + b) * t) + (1 - np.exp(-(a + b) * t)) * a / (a + b),
                    ],
                ]
            )

        dd = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            migration_matrix=[[0, a], [b, 0]],
        )
        return dd, f

    def verify_simulation(self, dd):
        samples = [
            msprime.Sample(pop, time)
            for pop in range(dd.num_populations)
            for time in dd.epoch_times
        ]
        reps = msprime.simulate(
            samples=samples,
            demography=dd.demography,
            end_time=max(dd.epoch_times + 1),
            num_replicates=100,
            random_seed=42,
        )
        for ts in reps:
            t = ts.first()
            for u in ts.samples():
                p = u
                lineage = []
                while p != tskit.NULL:
                    lineage.append(ts.node(p))
                    p = t.parent(p)
                probs = dd.lineage_probabilities(
                    [n.time for n in lineage], sample_time=lineage[0].time
                )
                locs = dd.possible_lineage_locations(
                    samples=[
                        msprime.Sample(time=n.time, population=n.population)
                        for n in lineage
                    ]
                )
                pop = lineage[0].population
                for j, n in enumerate(lineage):
                    self.assertGreater(probs[j, pop, n.population], 0.0)
                    for epoch in locs:
                        if n.time >= epoch[0] and n.time < epoch[1]:
                            self.assertTrue(locs[epoch][n.population])

    def test_two_pop(self):
        for _b in [2, 0]:
            dd, f = self.two_pop_example(1, 2)
            times = np.linspace(0, 10, 21)
            for st in [0.0, 2.5]:
                P = dd.lineage_probabilities(times + st, sample_time=st)
                for j, t in enumerate(times):
                    self.assertTrue(np.allclose(P[j, :, :], f(t)))
            self.verify_simulation(dd)

    def test_lineage_probabilities_tree(self):
        dem_events = [
            msprime.MassMigration(time=50, source=3, destination=2),
            msprime.MassMigration(time=100, source=2, destination=1),
            msprime.MassMigration(time=150, source=1, destination=0),
        ]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events, population_configurations=pop_config
        )
        P_out = dd.lineage_probabilities([10, 50, 60, 100, 101, 200])
        self.assertTrue(np.all([np.sum(P) == len(pop_config) for P in P_out]))
        self.assertTrue(np.all(np.diag(P_out[0]) == [1, 1, 1, 1]))
        self.assertTrue(np.all(np.diag(P_out[1]) == [1, 1, 1, 1]))
        self.assertTrue(np.all(probs == [1, 0, 0, 0] for probs in P_out[5]))
        self.verify_simulation(dd)

    def test_lineage_probabilities_pulse(self):
        f_pulse = 0.3
        dem_events = [
            msprime.MassMigration(time=1, source=1, destination=0, proportion=f_pulse),
        ]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events, population_configurations=pop_config
        )
        P_out = dd.lineage_probabilities([2])
        self.assertTrue(np.allclose(P_out[0], [[1, 0], [f_pulse, 1 - f_pulse]]))
        self.verify_simulation(dd)

    def test_lineage_probabilities_continuous_migration(self):
        mig_mat = [[0, 0.01], [0.01, 0]]
        dem_events = [
            msprime.MassMigration(time=100, source=1, destination=0, proportion=1),
            msprime.MigrationRateChange(time=100, rate=0),
        ]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events,
            population_configurations=pop_config,
            migration_matrix=mig_mat,
        )
        P_out = dd.lineage_probabilities([0, 50, 100, 150])
        self.assertTrue(np.all(P_out[0] == np.eye(len(pop_config))))
        self.assertTrue(np.all(P_out[1] > 0))
        self.assertTrue(np.all(P_out[1] > 0))
        # checking if close because of precision of _matrix_exponential function
        self.assertTrue(np.all(np.isclose(P_out[3], [[1, 0], [1, 0]])))
        self.verify_simulation(dd)

        mig_mat = [[0, 0.01], [0, 0]]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events,
            population_configurations=pop_config,
            migration_matrix=mig_mat,
        )
        P_out = dd.lineage_probabilities([0, 50, 100, 150])
        self.assertTrue(np.all(P_out[0] == np.eye(len(pop_config))))
        self.assertTrue(abs(P_out[1][1][0]) < np.finfo(float).eps)
        self.assertTrue(abs(P_out[2][1][0]) < np.finfo(float).eps)
        # machine precision instead of zero because of _matrix_exponential function
        self.assertTrue(np.all(np.isclose(P_out[3], [[1, 0], [1, 0]])))
        self.verify_simulation(dd)

    def test_sampling_time(self):
        mig_mat = [[0, 0.01], [0.02, 0]]
        dem_events = [
            msprime.MassMigration(time=100, source=1, destination=0, proportion=1),
            msprime.MigrationRateChange(time=100, rate=0),
        ]
        pop_config = [
            msprime.PopulationConfiguration(initial_size=100),
            msprime.PopulationConfiguration(initial_size=100),
        ]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events,
            population_configurations=pop_config,
            migration_matrix=mig_mat,
        )
        P_out = dd.lineage_probabilities([0], sample_time=1)
        self.assertTrue(np.all(P_out[0] == 0))
        P_out = dd.lineage_probabilities([0, 1, 2], sample_time=1)
        self.assertTrue(np.all(P_out[0] == 0))
        self.assertTrue(np.all(P_out[1] == [[1, 0], [0, 1]]))
        self.assertTrue(np.allclose(np.sum(P_out[2], axis=1), 1))
        P_out = dd.lineage_probabilities([99, 100, 101], sample_time=1)
        self.assertTrue(np.all(P_out[0] > 0))
        self.assertTrue(np.all(P_out[1] > 0))
        self.assertTrue(np.allclose(P_out[2], [[1, 0], [1, 0]]))
        self.verify_simulation(dd)


class TestPreCannedModels(unittest.TestCase):
    """
    Tests for the specialised models returned by static methods
    on the Demography.
    """

    def assertZeroDiagonal(self, A):
        self.assertTrue(np.all(np.diagonal(A)) == 0)


class TestIslandModel(TestPreCannedModels):
    def test_errors(self):
        for bad_N in [-1, 0, 0.1]:
            with self.assertRaises(ValueError):
                msprime.Demography.island_model(bad_N, 0.1)
        for bad_m in [-1, -1e5]:
            with self.assertRaises(ValueError):
                msprime.Demography.island_model(1, bad_m)

        for bad_Ne in [-1, 0]:
            with self.assertRaises(ValueError):
                msprime.Demography.island_model(1, 1, Ne=bad_Ne)

    def test_one_pop(self):
        model = msprime.Demography.island_model(1, 1)
        self.assertEqual(len(model.populations), 1)
        self.assertEqual(len(model.migration_matrix), 1)
        ts = msprime.simulate(samples=model.sample(2), demography=model, random_seed=1)
        self.assertEqual(ts.num_populations, 1)

    def test_migration(self):
        for N in [1, 2, 5]:
            model = msprime.Demography.island_model(N, 0.1)
            self.assertEqual(len(model.populations), N)
            self.assertEqual(model.migration_matrix.shape, (N, N))
            self.assertZeroDiagonal(model.migration_matrix)
            self.assertTrue(
                np.all(model.migration_matrix[~np.eye(N, dtype=bool)] == 0.1)
            )
            ts = msprime.simulate(
                samples=model.sample(*([2] * N)), demography=model, random_seed=1
            )
            self.assertEqual(ts.num_populations, N)
            self.assertEqual(ts.num_samples, 2 * N)

    def test_Ne(self):
        # By default, Ne is 1
        model = msprime.Demography.island_model(2, 0.1)
        self.assertEqual(model.populations[0].initial_size, 1)
        self.assertEqual(model.populations[1].initial_size, 1)
        for Ne in [0.1, 1, 10]:
            model = msprime.Demography.island_model(2, 0.1, Ne=Ne)
            self.assertEqual(model.populations[0].initial_size, Ne)
            self.assertEqual(model.populations[1].initial_size, Ne)


class TestSteppingStoneModel(TestPreCannedModels):
    def test_errors(self):
        for bad_N in [-1, 0, 0.1]:
            with self.assertRaises(ValueError):
                msprime.Demography.stepping_stone_1d(bad_N, 0.1)
        for bad_m in [-1, -1e5]:
            with self.assertRaises(ValueError):
                msprime.Demography.stepping_stone_1d(1, bad_m)

        for bad_Ne in [-1, 0]:
            with self.assertRaises(ValueError):
                msprime.Demography.stepping_stone_1d(1, 1, Ne=bad_Ne)

    def test_one_pop(self):
        for circular in [True, False]:
            model = msprime.Demography.stepping_stone_1d(1, 1, circular=circular)
            self.assertEqual(len(model.populations), 1)
            self.assertEqual(len(model.migration_matrix), 1)
            ts = msprime.simulate(
                samples=model.sample(2), demography=model, random_seed=1
            )
            self.assertEqual(ts.num_populations, 1)

    def test_migration_circular(self):
        m = 0.3
        for N in [2, 3, 5]:
            model = msprime.Demography.stepping_stone_1d(N, m)
            # Circular is the default
            self.assertEqual(
                model, msprime.Demography.stepping_stone_1d(N, m, circular=True)
            )
            self.assertEqual(len(model.populations), N)
            self.assertEqual(model.migration_matrix.shape, (N, N))
            self.assertZeroDiagonal(model.migration_matrix)
            for j in range(N):
                adjacent = [(j - 1) % N, (j + 1) % N]
                for k in range(N):
                    if k in adjacent:
                        self.assertEqual(model.migration_matrix[j, k], m)
                    else:
                        self.assertEqual(model.migration_matrix[j, k], 0)
            ts = msprime.simulate(
                samples=model.sample(*([2] * N)), demography=model, random_seed=1
            )
            self.assertEqual(ts.num_populations, N)
            self.assertEqual(ts.num_samples, 2 * N)

    def test_migration_line_two_pops(self):
        m = 1
        model = msprime.Demography.stepping_stone_1d(2, m, circular=False)
        self.assertEqual(len(model.populations), 2)
        self.assertEqual(model.migration_matrix.shape, (2, 2))
        self.assertTrue(np.all(model.migration_matrix == 0))

    def test_migration_line(self):
        m = 0.3
        for N in [3, 4, 5]:
            model = msprime.Demography.stepping_stone_1d(N, m, circular=False)
            self.assertEqual(len(model.populations), N)
            self.assertEqual(model.migration_matrix.shape, (N, N))
            self.assertZeroDiagonal(model.migration_matrix)
            for j in range(N):
                adjacent = []
                if j > 0:
                    adjacent.append(j - 1)
                if j < N - 1:
                    adjacent.append(j + 1)
                for k in range(N):
                    if k in adjacent:
                        self.assertEqual(model.migration_matrix[j, k], m)
                    else:
                        self.assertEqual(model.migration_matrix[j, k], 0)
            ts = msprime.simulate(
                samples=model.sample(*([2] * N)), demography=model, random_seed=1
            )
            self.assertEqual(ts.num_populations, N)
            self.assertEqual(ts.num_samples, 2 * N)

    def test_Ne(self):
        # By default, Ne is 1
        model = msprime.Demography.stepping_stone_1d(2, 0.1)
        self.assertEqual(model.populations[0].initial_size, 1)
        self.assertEqual(model.populations[1].initial_size, 1)
        for Ne in [0.1, 1, 10]:
            model = msprime.Demography.stepping_stone_1d(2, 0.1, Ne=Ne)
            self.assertEqual(model.populations[0].initial_size, Ne)
            self.assertEqual(model.populations[1].initial_size, Ne)


class TestDemographyObject(unittest.TestCase):
    """
    Basic tests for the demography object.
    """

    def test_equality(self):
        m1 = msprime.Demography.island_model(2, 1 / 3)
        m2 = msprime.Demography.island_model(2, 1 / 3)
        self.assertEqual(m1, m2)
        self.assertEqual(m2, m1)
        self.assertEqual(m1, m1)
        self.assertFalse(m1 != m2)
        self.assertFalse(m1 != m1)

        m3 = msprime.Demography.island_model(2, 1 / 3 + 0.001)
        self.assertNotEqual(m1, m3)
        self.assertTrue(m1 != m3)

        self.assertNotEqual(m1, None)
        self.assertNotEqual(m1, [])

    def test_debug(self):
        model = msprime.Demography.island_model(2, 1 / 3)
        dbg1 = model.debug()
        self.assertEqual(dbg1.demography, model)
        dbg2 = msprime.DemographyDebugger(demography=model)
        self.assertEqual(dbg1.demography, dbg2.demography)
        self.assertEqual(str(dbg1), str(dbg2))

    def test_positional_sampling_errors(self):
        model = msprime.Demography.island_model(2, 1)
        with self.assertRaises(ValueError):
            # Sampling from no populations is an error (this is almost
            # certainly a mistake by the user).
            model.sample()
        for bad_sample in [(1, -1), (-1,), (0, -10)]:
            with self.assertRaises(ValueError):
                model.sample(*bad_sample)
        with self.assertRaises(TypeError):
            model.sample(6.6)
        with self.assertRaises(ValueError):
            model.sample(0, 0, 1)

    def test_positional_samples_two_populations(self):
        model = msprime.Demography.island_model(2, 1)
        self.assertEqual(model.sample(1), [msprime.Sample(0, 0)])
        self.assertEqual(model.sample(0, 1), [msprime.Sample(1, 0)])
        self.assertEqual(
            model.sample(1, 1), [msprime.Sample(0, 0), msprime.Sample(1, 0)]
        )
        self.assertEqual(
            model.sample(2, 0), [msprime.Sample(0, 0), msprime.Sample(0, 0)]
        )
        # Drawing 0 samples is OK
        self.assertEqual(model.sample(0), [])
        self.assertEqual(model.sample(0, 0), [])
        self.assertEqual(
            model.sample(3, 1), [msprime.Sample(0, 0)] * 3 + [msprime.Sample(1, 0)]
        )

    def test_positional_samples_n_populations(self):
        for n in [1, 2, 3, 5]:
            model = msprime.Demography.island_model(n, 1)
            samples = model.sample(10)
            self.assertEqual(samples, [msprime.Sample(0, 0)] * 10)
            samples = model.sample(*np.ones(n, dtype=int))
            self.assertEqual(samples, [msprime.Sample(j, 0) for j in range(n)])
            samples = model.sample(*np.zeros(n, dtype=int))
            self.assertEqual(samples, [])
            samples = model.sample(*range(n))
            self.assertEqual(
                samples,
                list(itertools.chain(*[[msprime.Sample(j, 0)] * j for j in range(n)])),
            )

    def test_keyword_sampling_errors(self):
        model = msprime.Demography.island_model(2, 1)
        model.populations[0].name = "A"
        model.populations[1].name = "B"
        with self.assertRaises(ValueError):
            # Sampling from no populations is an error (this is almost
            # certainly a mistake by the user).
            model.sample(**{})
        for bad_sample in [{"A": 1, "B": -1}, {"A": -1}, {"A": 0, "B": -10}]:
            with self.assertRaises(ValueError):
                model.sample(**bad_sample)
        with self.assertRaises(TypeError):
            model.sample(A=6.6)
        with self.assertRaises(ValueError):
            model.sample(C=1)
        with self.assertRaises(ValueError):
            model.sample(**{"AC": 1})

    def test_keyword_samples_two_populations(self):
        model = msprime.Demography.island_model(2, 1)
        model.populations[0].name = "A"
        model.populations[1].name = "B"
        self.assertEqual(model.sample(A=1), [msprime.Sample(0, 0)])
        self.assertEqual(model.sample(B=1), [msprime.Sample(1, 0)])
        self.assertEqual(
            model.sample(A=1, B=1), [msprime.Sample(0, 0), msprime.Sample(1, 0)]
        )
        # Samples are returned **in the order specified**. This is guaranteed
        # since Python 3.6
        self.assertEqual(
            model.sample(B=1, A=1), [msprime.Sample(1, 0), msprime.Sample(0, 0)]
        )
        self.assertEqual(
            model.sample(A=2, B=0), [msprime.Sample(0, 0), msprime.Sample(0, 0)]
        )
        # Drawing 0 samples is OK
        self.assertEqual(model.sample(A=0), [])
        self.assertEqual(model.sample(B=0), [])
        self.assertEqual(model.sample(A=0, B=0), [])
        self.assertEqual(
            model.sample(A=3, B=1), [msprime.Sample(0, 0)] * 3 + [msprime.Sample(1, 0)]
        )

    def test_mixed_positional_and_keyword(self):
        model = msprime.Demography.island_model(2, 1)
        model.populations[0].name = "A"
        model.populations[1].name = "B"
        with self.assertRaises(ValueError):
            model.sample(0, A=1)


class TestDemographyFromOldStyle(unittest.TestCase):
    """
    Tests the method for creating a demography object from the old
    style population_configurations, migration_matrix and demographic_events
    parameters.
    """

    def test_defaults(self):
        demog = msprime.Demography.from_old_style()
        self.assertEqual(demog.num_populations, 1)
        self.assertEqual(list(demog.migration_matrix), [[0]])
        self.assertEqual(list(demog.events), [])

    def test_pop_configs_defaults(self):
        for n in range(1, 5):
            pop_configs = [msprime.PopulationConfiguration() for _ in range(n)]
            demog = msprime.Demography.from_old_style(
                population_configurations=pop_configs
            )
            self.assertEqual(demog.num_populations, n)
            np.testing.assert_array_equal(demog.migration_matrix, np.zeros((n, n)))
            self.assertEqual(list(demog.events), [])

    def test_migration_matrix(self):
        for n in range(1, 5):
            pop_configs = [msprime.PopulationConfiguration() for _ in range(n)]
            M = np.ones((n, n))
            demog = msprime.Demography.from_old_style(
                population_configurations=pop_configs, migration_matrix=M
            )
            self.assertEqual(demog.num_populations, n)
            np.testing.assert_array_equal(demog.migration_matrix, M)
            self.assertEqual(list(demog.events), [])

    def test_demographic_events(self):
        events = [
            msprime.PopulationParametersChange(time=j, initial_size=2)
            for j in range(10)
        ]
        demog = msprime.Demography.from_old_style(demographic_events=events)
        self.assertEqual(demog.num_populations, 1)
        self.assertEqual(list(demog.migration_matrix), [[0]])
        self.assertEqual(events, demog.events)


class TestPopulationFromOldStyle(unittest.TestCase):
    """
    Tests the method for creating a Population object from the old
    style PopulationConfiguration.
    """

    # TODO figure out what to do with metadata

    def test_defaults(self):
        pop_config = msprime.PopulationConfiguration()
        pop = msprime.Population.from_old_style(pop_config)
        self.assertEqual(pop_config.initial_size, pop.initial_size)
        self.assertEqual(pop_config.growth_rate, pop.growth_rate)

    def test_values(self):
        pop_config = msprime.PopulationConfiguration(
            initial_size=1234, growth_rate=5678
        )
        pop = msprime.Population.from_old_style(pop_config)
        self.assertEqual(pop_config.initial_size, pop.initial_size)
        self.assertEqual(pop_config.growth_rate, pop.growth_rate)
