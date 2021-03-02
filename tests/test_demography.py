#
# Copyright (C) 2016-2021 University of Oxford
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
import textwrap
import unittest
import warnings
import xml
from unittest import mock

import numpy as np
import pytest
import scipy.linalg
import stdpopsim
import tskit

import msprime
from msprime import _msprime


def all_events_example_demography(*, integer_ids=False):
    demography = msprime.Demography.isolated_model([10] * 10)
    demography.add_population_parameters_change(0.1, initial_size=2)
    demography.add_population_parameters_change(0.1, growth_rate=10)
    demography.add_migration_rate_change(0.4, rate=0)
    demography.add_census(0.55)
    if integer_ids:
        demography.add_migration_rate_change(0.2, source=0, dest=1, rate=1)
        demography.add_symmetric_migration_rate_change(
            0.3, populations=[0, 1], rate=0.5
        )
        demography.add_mass_migration(0.4, source=1, dest=0, proportion=0.5)
        demography.add_population_split(0.4, derived=[3, 4], ancestral=5)
        demography.add_admixture(
            0.45, derived=7, ancestral=[8, 9], proportions=[0.5, 0.5]
        )
        demography.add_instantaneous_bottleneck(0.5, population=0, strength=100)
        demography.add_simple_bottleneck(0.56, population=1, proportion=0.1)
    else:
        demography.add_migration_rate_change(0.2, source="pop_0", dest="pop_1", rate=1)
        demography.add_symmetric_migration_rate_change(
            0.3, populations=["pop_0", "pop_1"], rate=0.5
        )
        demography.add_mass_migration(0.4, source="pop_1", dest="pop_0", proportion=0.5)
        demography.add_population_split(
            0.4, derived=["pop_3", "pop_4"], ancestral="pop_5"
        )
        demography.add_admixture(
            0.45, derived="pop_7", ancestral=["pop_8", "pop_9"], proportions=[0.5, 0.5]
        )
        demography.add_instantaneous_bottleneck(0.5, population="pop_0", strength=100)
        demography.add_simple_bottleneck(0.56, population="pop_1", proportion=0.1)
    demography.sort_events()
    return demography


class TestNePopulationSizeEquivalence:
    """
    Test that setting Ne as a parameter of the population model and
    at individual populations is the same thing.
    """

    def assert_tree_sequences_equal(self, ts1, ts2):
        t1 = ts1.dump_tables()
        t2 = ts2.dump_tables()
        # We can't compare directly byte-for-byte because there'll be slight
        # differences in the computed times.
        assert len(t1.nodes) == len(t2.nodes)
        np.testing.assert_array_equal(t1.nodes.flags, t2.nodes.flags)
        np.testing.assert_array_almost_equal(t1.nodes.time, t2.nodes.time)
        np.testing.assert_array_equal(t1.nodes.population, t2.nodes.population)
        np.testing.assert_array_equal(t1.nodes.individual, t2.nodes.individual)
        np.testing.assert_array_equal(t1.nodes.metadata, t2.nodes.metadata)
        np.testing.assert_array_equal(
            t1.nodes.metadata_offset, t2.nodes.metadata_offset
        )
        assert len(t1.edges) == len(t2.edges)
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
            assert ts1.num_trees > 1
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


class TestIntrospectionInterface:
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
        assert repr(event) == repr_s
        assert str(event) == repr_s

    def test_migration_rate_change(self):
        event = msprime.MigrationRateChange(time=1, rate=2)
        repr_s = "MigrationRateChange(time=1, rate=2, source=-1, dest=-1)"
        assert repr(event) == repr_s
        assert str(event) == repr_s

    def test_symmetric_migration_rate_change(self):
        event = msprime.SymmetricMigrationRateChange(time=1, populations=[0, 1], rate=2)
        repr_s = "SymmetricMigrationRateChange(time=1, populations=[0, 1], rate=2)"
        assert repr(event) == repr_s
        assert str(event) == repr_s

    def test_mass_migration(self):
        event = msprime.MassMigration(time=1, proportion=0.5, source=0, dest=1)
        repr_s = "MassMigration(time=1, source=0, dest=1, proportion=0.5)"
        assert repr(event) == repr_s
        assert str(event) == repr_s

    def test_population_split(self):
        event = msprime.PopulationSplit(time=1, derived=[0], ancestral=1)
        repr_s = "PopulationSplit(time=1, derived=[0], ancestral=1)"
        assert repr(event) == repr_s
        assert str(event) == repr_s

    def test_simple_bottleneck(self):
        event = msprime.SimpleBottleneck(time=1, population=1, proportion=0.5)
        repr_s = "SimpleBottleneck(time=1, population=1, proportion=0.5)"
        assert repr(event) == repr_s
        assert str(event) == repr_s

    def test_instantaneous_bottleneck(self):
        event = msprime.InstantaneousBottleneck(time=1, population=1, strength=1.5)
        repr_s = "InstantaneousBottleneck(time=1, population=1, strength=1.5)"
        assert repr(event) == repr_s
        assert str(event) == repr_s

    def test_census(self):
        event = msprime.CensusEvent(time=1)
        repr_s = "CensusEvent(time=1)"
        assert repr(event) == repr_s
        assert str(event) == repr_s


class TestDemographicEventsHaveExtraLLParameter:
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
        demography = all_events_example_demography()
        for event in demography.events:
            ll_config1 = event.get_ll_representation()
            ll_config2 = event.get_ll_representation(None)
            assert ll_config1 == ll_config2


class TestTimeTravelErrors:
    """
    It is possible to specify models in msprime that result in malformed
    tree sequences where the parent node has time equal to its child.
    We throw an error in this case and expect the user to fix their model.
    """

    def test_multiple_bottlenecks(self):
        with pytest.raises(_msprime.LibraryError):
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
        with pytest.raises(_msprime.LibraryError):
            msprime.simulate(
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                migration_matrix=M,
                recombination_rate=0.0,
                mutation_rate=0.0,
                random_seed=1,
            )


class TestBadDemographicParameters:
    """
    Tests for nonsensical demographic parameters.
    """

    def test_bad_population_size(self):
        for bad_size in [-1, -1e300]:
            with pytest.raises(ValueError):
                msprime.PopulationParametersChange(0, initial_size=bad_size)
            with pytest.raises(ValueError):
                msprime.PopulationConfiguration(
                    initial_size=bad_size,
                    sample_size=2,
                )

    def test_bad_sample_size(self):
        for bad_size in [-1, -1e300]:
            with pytest.raises(ValueError):
                msprime.PopulationConfiguration(
                    initial_size=1,
                    sample_size=bad_size,
                )

    def test_dtwf_bottleneck(self):
        with pytest.raises(_msprime.LibraryError):
            msprime.simulate(
                sample_size=2,
                model="dtwf",
                demographic_events=[msprime.SimpleBottleneck(time=0.1, population=0)],
                random_seed=1,
            )


class TestBadDemographicEvents:
    """
    Tests for input errors when creating demographic events.
    """

    def test_growth_rate_or_initial_size(self):
        with pytest.raises(ValueError):
            msprime.PopulationParametersChange(time=0)

    def test_bad_simulation_model(self):
        for model in [[], {}]:
            des = [msprime.SimulationModelChange(time=0, model=model)]
            with pytest.raises(TypeError):
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
            assert all(tree.num_roots == 1 for tree in ts.trees()), f"{samples}"

    def test_cant_sample_a_zero_sized_population(self):
        for bad_samples, err in [
            ([msprime.Sample(0, 0)], _msprime.InputError),
            # Check ancient sampling.
            ([msprime.Sample(1, self.T), msprime.Sample(1, 0)], _msprime.LibraryError),
        ]:
            with pytest.raises(err, match="Bad sample"):
                msprime.simulate(
                    population_configurations=self.population_configurations,
                    demographic_events=self.demographic_events,
                    samples=bad_samples * 10,
                )


class TestDeprecatedInterfaces:
    """
    Tests to check that aliased parameters are handled correctly.
    """

    def test_mass_migration_dest(self):
        with pytest.raises(ValueError):
            msprime.MassMigration(time=0, source=0, dest=0, destination=0)
        for j in range(10):
            e = msprime.MassMigration(time=0.1, source=0, dest=j, proportion=0.5)
            assert e.time == 0.1
            assert e.source == 0
            assert e.dest == j
            assert e.proportion == 0.5
            e = msprime.MassMigration(0.1, 0, j, 0.5)
            assert e.time == 0.1
            assert e.source == 0
            assert e.dest == j
            assert e.proportion == 0.5
            e = msprime.MassMigration(time=0.1, source=0, destination=j, proportion=0.5)
            assert e.time == 0.1
            assert e.source == 0
            assert e.dest == j
            assert e.proportion == 0.5

    def test_population_parameters_population_id(self):
        with pytest.raises(ValueError):
            msprime.PopulationParametersChange(
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
            assert e.time == 0.1
            assert e.initial_size == 1.1
            assert e.growth_rate == 0.1
            assert e.population == j
            e = msprime.PopulationParametersChange(0.1, 1.1, 0.1, j)
            assert e.time == 0.1
            assert e.initial_size == 1.1
            assert e.growth_rate == 0.1
            assert e.population == j
            e = msprime.PopulationParametersChange(
                time=0.1, initial_size=1.1, growth_rate=0.1, population_id=j
            )
            assert e.time == 0.1
            assert e.initial_size == 1.1
            assert e.growth_rate == 0.1
            assert e.population == j

    def test_epoch_demographic_events(self):
        demography = msprime.Demography.isolated_model([10] * 2)
        demography.add_migration_rate_change(0.1, rate=1)
        dbg = demography.debug()
        assert dbg.num_epochs == 2
        for epoch in dbg.epochs:
            assert epoch.demographic_events == epoch.events


class TestEventsWithoutDemography:
    """
    Tests for "raw" DemographicEvent objects that are not associated
    with a Demography. These are supported for legacy reasons - downstream
    code like stdpopsim relied on them.
    """

    def test_convert_id(self):
        event = msprime.DemographicEvent(0.0)
        for j in range(5):
            assert event._convert_id(j) == j
            numpy_val = np.array([j])[0]
            assert event._convert_id(numpy_val) == j

        assert event._convert_id(-1) == -1
        assert event._convert_id(None) == -1

        with pytest.raises(ValueError, match="legacy-only operation"):
            event._convert_id("pop_0")


class TestLowLevelRepresentation:
    """
    Tests for the low-level demographic events interface.

    Some of these tests are needed to make sure that the _convert_id workaround
    works when we don't have a demography object specified, but we can
    get rid of these when we get stdpopsim fixed up.

    https://github.com/tskit-dev/msprime/issues/1037
    """

    def test_size_change(self):
        g = 100
        new_size = 512
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_population_parameters_change(
            time=g, initial_size=new_size, population="pop_0"
        )
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": 0,
            "initial_size": new_size,
        }
        assert event.get_ll_representation() == ll_event

    def test_growth_rate_change(self):
        g = 512
        growth_rate = 1
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_population_parameters_change(
            time=g, growth_rate=growth_rate, population="pop_1"
        )
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": 1,
            "growth_rate": growth_rate,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_population_parameters_change(
            time=g,
            growth_rate=growth_rate,
            population=1,
        )
        assert event.get_ll_representation() == ll_event

    def test_growth_rate_and_size_change(self):
        g = 1024
        growth_rate = 2
        initial_size = 8192
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_population_parameters_change(
            time=g, initial_size=initial_size, growth_rate=growth_rate, population=1
        )
        ll_event = {
            "type": "population_parameters_change",
            "time": g,
            "population": 1,
            "initial_size": initial_size,
            "growth_rate": growth_rate,
        }
        assert event.get_ll_representation() == ll_event

    def test_migration_rate_change(self):
        g = 1024
        migration_rate = 0.125
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_migration_rate_change(
            time=g, source="pop_0", dest="pop_1", rate=migration_rate
        )
        ll_event = {
            "type": "migration_rate_change",
            "time": g,
            "source": 0,
            "dest": 1,
            "migration_rate": migration_rate,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_migration_rate_change(
            time=g, source=0, dest=1, rate=migration_rate
        )
        assert event.get_ll_representation() == ll_event

    def test_migration_rate_change_all_pops(self):
        g = 1024
        migration_rate = 0.125
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_migration_rate_change(
            time=g, source=None, dest=None, rate=migration_rate
        )
        ll_event = {
            "type": "migration_rate_change",
            "time": g,
            "source": -1,
            "dest": -1,
            "migration_rate": migration_rate,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_migration_rate_change(time=g, rate=migration_rate)
        assert event.get_ll_representation() == ll_event

    def test_symmetric_migration_rate_change(self):
        g = 1024
        migration_rate = 0.125
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_symmetric_migration_rate_change(
            time=g, populations=["pop_0", "pop_1"], rate=migration_rate
        )
        ll_event = {
            "type": "symmetric_migration_rate_change",
            "time": g,
            "populations": [0, 1],
            "rate": migration_rate,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_symmetric_migration_rate_change(
            time=g, populations=[0, 1], rate=migration_rate
        )
        assert event.get_ll_representation() == ll_event

    def test_mass_migration(self):
        g = 1234
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_mass_migration(
            time=g, source="pop_0", dest="pop_1", proportion=0.5
        )
        ll_event = {
            "type": "mass_migration",
            "time": g,
            "source": 0,
            "dest": 1,
            "proportion": 0.5,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_mass_migration(time=g, source=0, dest=1, proportion=0.5)
        assert event.get_ll_representation() == ll_event

    def test_population_split(self):
        g = 1234
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_population_split(
            time=g, derived=["pop_0"], ancestral="pop_1"
        )
        ll_event = {
            "type": "population_split",
            "time": g,
            "derived": [0],
            "ancestral": 1,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_population_split(time=g, derived=[0], ancestral="pop_1")
        assert event.get_ll_representation() == ll_event

    def test_admixture(self):
        g = 1234
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_admixture(
            time=g,
            derived="pop_0",
            ancestral=["pop_1", "pop_2"],
            proportions=[0.25, 0.75],
        )
        ll_event = {
            "type": "admixture",
            "time": g,
            "derived": 0,
            "ancestral": [1, 2],
            "proportions": [0.25, 0.75],
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_admixture(
            time=g, derived=0, ancestral=[1, 2], proportions=[0.25, 0.75]
        )
        assert event.get_ll_representation() == ll_event

    def test_simple_bottleneck(self):
        g = 1234
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_simple_bottleneck(
            time=g, population="pop_0", proportion=0.5
        )
        ll_event = {
            "type": "simple_bottleneck",
            "time": g,
            "population": 0,
            "proportion": 0.5,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_simple_bottleneck(time=g, population=0, proportion=0.5)
        assert event.get_ll_representation() == ll_event

    def test_instantaneous_bottleneck(self):
        g = 1234
        demography = msprime.Demography.isolated_model([10] * 10)
        event = demography.add_instantaneous_bottleneck(
            time=g, population="pop_0", strength=0.5
        )
        ll_event = {
            "type": "instantaneous_bottleneck",
            "time": g,
            "population": 0,
            "strength": 0.5,
        }
        assert event.get_ll_representation() == ll_event
        event = demography.add_instantaneous_bottleneck(
            time=g, population=0, strength=0.5
        )
        assert event.get_ll_representation() == ll_event


class TestDemographyDebugger:
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
        assert dd.num_epochs == times.shape[0]
        assert dd.num_epochs == pop_size.shape[1]
        assert dd.num_populations == pop_size.shape[0]
        for j in range(dd.num_epochs):
            assert dd.epochs[j].start_time == times[j]
            for k in range(dd.num_populations):
                assert dd.epochs[j].populations[k].start_size == pop_size[k, j]

    def test_equal_after(self):
        population_configurations = [
            msprime.PopulationConfiguration(sample_size=10),
            msprime.PopulationConfiguration(sample_size=20),
        ]
        msprime.DemographyDebugger(population_configurations=population_configurations)
        assert population_configurations[0].sample_size == 10
        assert population_configurations[1].sample_size == 20

    def test_model_change_events(self):
        population_configurations = [msprime.PopulationConfiguration(sample_size=10)]
        demographic_events = [msprime.SimulationModelChange(1, "hudson")]
        with pytest.raises(TypeError):
            msprime.DemographyDebugger(
                population_configurations=population_configurations,
                demographic_events=demographic_events,
            )

    def test_one_pop_zero_events(self):
        dd = msprime.DemographyDebugger(
            population_configurations=[msprime.PopulationConfiguration()]
        )
        assert len(dd.epochs) == 1
        self.verify_arrays(dd)
        e = dd.epochs[0]
        assert e.start_time == 0
        assert dd.epoch_times[0] == 0
        assert dd.population_size_history.shape[0] == 1
        assert math.isinf(e.end_time)
        assert len(e.events) == 0
        assert len(e.populations) == 1
        assert e.migration_matrix == [[0]]
        pop = e.populations[0]
        assert pop.growth_rate == 0
        assert pop.start_size == 1
        assert pop.end_size == 1

    def test_two_pop_different_sizes(self):
        dd = msprime.DemographyDebugger(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=10),
                msprime.PopulationConfiguration(initial_size=20),
            ]
        )
        self.verify_arrays(dd)
        assert len(dd.epochs) == 1
        e = dd.epochs[0]
        assert e.start_time == 0
        assert dd.population_size_history.shape[0] == 2
        assert dd.population_size_history[0][0] == 10
        assert dd.population_size_history[1][0] == 20
        assert math.isinf(e.end_time)
        assert len(e.events) == 0
        assert len(e.populations) == 2
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
        for pop in e.populations:
            assert pop.growth_rate == 0
            assert pop.start_size == pop.end_size
        assert e.populations[0].start_size == 10
        assert e.populations[1].start_size == 20

    def test_two_pop_different_growth_rates(self):
        g1 = 0.1
        g2 = 0.5
        p0_end_size = 10 * math.exp(-g1 * 10)
        p1_end_size = 20 * math.exp(-g2 * 10)
        demog = msprime.Demography.from_old_style(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=10, growth_rate=g1),
                msprime.PopulationConfiguration(initial_size=20, growth_rate=g2),
            ],
            demographic_events=[
                msprime.PopulationParametersChange(time=10, growth_rate=0)
            ],
        )
        dd = demog.debug()
        self.verify_arrays(dd)
        # Make sure we're testing the __repr__ paths.
        s = repr(dd)
        assert len(s) > 0
        assert len(dd.epochs) == 2
        e = dd.epochs[0]
        assert e.start_time == 0
        assert e.end_time == 10
        assert len(e.events) == 0
        assert len(e.populations) == 2
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
        assert e.populations[0].start_size == 10
        assert e.populations[0].end_size == p0_end_size
        assert e.populations[1].start_size == 20
        assert e.populations[1].end_size == p1_end_size

        e = dd.epochs[1]
        assert e.start_time == 10
        assert math.isinf(e.end_time)
        assert len(e.events) == 1
        d = e.events[0]
        assert d.time == 10
        assert d.growth_rate == 0
        assert d.initial_size is None
        assert d.population == -1
        assert len(e.populations) == 2
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
        for pop in e.populations:
            assert pop.growth_rate == 0
            assert pop.start_size == pop.end_size
        assert e.populations[0].start_size == p0_end_size
        assert e.populations[1].start_size == p1_end_size

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
        assert len(dd.epochs) == 3
        e = dd.epochs[0]
        assert e.start_time == 0
        assert e.end_time == 20
        np.testing.assert_array_equal(e.migration_matrix, [[0, 0.25], [0, 0]])
        for pop in e.populations:
            assert pop.growth_rate == 0
            assert pop.start_size == pop.end_size
            assert pop.start_size == 100

        e = dd.epochs[1]
        assert e.start_time == 20
        assert e.end_time, 22
        assert len(e.events) == 1
        d = e.events[0]
        assert d.source == -1
        assert d.dest == -1
        assert d.time == 20
        assert d.rate == 1
        assert len(e.populations) == 2
        np.testing.assert_array_equal(e.migration_matrix, [[0, 1], [1, 0]])
        for pop in e.populations:
            assert pop.growth_rate == 0
            assert pop.start_size == pop.end_size
            assert pop.start_size == 100

        e = dd.epochs[2]
        assert e.start_time == 22
        assert math.isinf(e.end_time)
        assert len(e.events) == 1
        d = e.events[0]
        assert d.source == 0
        assert d.dest == 1
        assert d.time == 22
        assert d.rate == 1.7
        assert len(e.populations) == 2
        np.testing.assert_array_equal(e.migration_matrix, [[0, 1.7], [1, 0]])
        for pop in e.populations:
            assert pop.growth_rate == 0
            assert pop.start_size == pop.end_size
            assert pop.start_size == 100

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

            assert len(dd.epochs) == 4
            e = dd.epochs[0]
            assert e.start_time == 0
            assert e.end_time == t1
            assert dd.epoch_times[0] == 0
            assert dd.epoch_times[1] == t1
            assert len(e.events) == 0
            assert len(e.populations) == 2
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
            assert e.populations[0].start_size == N0
            n0 = N0 * math.exp(-alpha * t1)
            assert e.populations[0].end_size == n0
            assert e.populations[1].start_size == N1
            assert e.populations[1].end_size == N1

            e = dd.epochs[1]
            assert e.start_time == t1
            assert e.end_time == t2
            assert len(e.events) == 1
            assert len(e.populations) == 2
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
            assert e.populations[0].start_size == n0
            n0 = N0 * math.exp(-alpha * t2)
            assert e.populations[0].end_size == n0
            assert e.populations[1].start_size == N1
            n1 = N1 * math.exp(-alpha * (t2 - t1))
            assert e.populations[1].end_size == n1

            e = dd.epochs[2]
            assert e.start_time == t2
            assert e.end_time == t3
            assert len(e.events) == 1
            assert len(e.populations) == 2
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
            assert e.populations[0].start_size == n0
            n0 = n0 * math.exp(alpha * (t3 - t2))
            assert e.populations[0].end_size == n0
            assert e.populations[1].start_size == n1
            n1 = N1 * math.exp(-alpha * (t3 - t1))
            assert e.populations[1].end_size == n1

            e = dd.epochs[3]
            assert e.start_time == t3
            assert math.isinf(e.end_time)
            assert len(e.events) == 1
            assert len(e.populations) == 2
            np.testing.assert_array_equal(e.migration_matrix, [[0, 0], [0, 0]])
            assert e.populations[0].start_size == n0
            assert e.populations[0].end_size == n0
            assert e.populations[1].start_size == n1
            assert e.populations[1].end_size == n1

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
            assert len(w) == 1

    def test_no_misspecification_no_warning(self):
        with warnings.catch_warnings(record=True) as w:
            self.check_model_misspecification_warning(misspecify=False)
            assert len(w) == 0


class TestDemographicEventMessages:
    def test_population_parameters_change(self):
        event = msprime.PopulationParametersChange(1.0, population=1, initial_size=2.0)
        assert event._parameters() == "population=1, initial_size=2.0"
        assert event._effect() == "initial_size → 2 for population 1"

        event = msprime.PopulationParametersChange(
            1.0, population="XX", growth_rate=2.0
        )
        assert event._parameters() == "population=XX, growth_rate=2.0"
        assert event._effect() == "growth_rate → 2 for population XX"

        event = msprime.PopulationParametersChange(
            1.0, population=0, initial_size=3, growth_rate=2.0
        )
        assert event._parameters() == "population=0, initial_size=3, growth_rate=2.0"
        assert (
            event._effect() == "initial_size → 3 and growth_rate → 2 for population 0"
        )

        for pop in [None, -1]:
            event = msprime.PopulationParametersChange(
                1.0, population=pop, growth_rate=2.0
            )
            assert event._parameters() == "population=-1, growth_rate=2.0"
            assert event._effect() == "growth_rate → 2 for all populations"

    def test_migration_rate_change(self):
        event = msprime.MigrationRateChange(time=1, rate=2)
        assert event._parameters() == "source=-1, dest=-1, rate=2"
        assert (
            event._effect() == "Backwards-time migration rate for all populations → 2"
        )

        event = msprime.MigrationRateChange(source=0, dest=1, time=1, rate=6)
        assert event._parameters() == "source=0, dest=1, rate=6"
        assert event._effect() == "Backwards-time migration rate from 0 to 1 → 6"

    def test_symmetric_migration_rate_change(self):
        event = msprime.SymmetricMigrationRateChange(time=1, populations=[0, 1], rate=2)
        assert event._parameters() == "populations=[0, 1], rate=2"
        assert event._effect() == (
            "Sets the symmetric migration rate between 0 and 1 to 2 per generation"
        )
        event = msprime.SymmetricMigrationRateChange(
            time=1, populations=[0, 1, 2], rate=2
        )
        assert event._parameters() == "populations=[0, 1, 2], rate=2"
        assert event._effect() == (
            "Sets the symmetric migration rate between all pairs of populations "
            "in [0, 1, 2] to 2 per generation"
        )

    def test_population_split(self):
        event = msprime.PopulationSplit(time=1, derived=[0], ancestral=2)
        assert event._parameters() == "derived=[0], ancestral=2"
        assert event._effect() == (
            "Moves all lineages from the '0' derived population to the "
            "ancestral '2' population. Also set '0' to inactive, and "
            "all migration rates to and from the derived population to zero."
        )

        event = msprime.PopulationSplit(time=1, derived=[0, 1], ancestral=2)
        assert event._parameters() == "derived=[0, 1], ancestral=2"
        assert event._effect() == (
            "Moves all lineages from derived populations '0' and '1' to the "
            "ancestral '2' population. Also set the derived populations to inactive, "
            "and all migration rates to and from the derived populations to zero."
        )

        event = msprime.PopulationSplit(time=1, derived=[0, 1, 2], ancestral=3)
        assert event._parameters() == "derived=[0, 1, 2], ancestral=3"
        assert event._effect() == (
            "Moves all lineages from derived populations [0, 1, 2] to the "
            "ancestral '3' population. Also set the derived populations to inactive, "
            "and all migration rates to and from the derived populations to zero."
        )

    def test_admixture(self):
        event = msprime.Admixture(time=1, derived=0, ancestral=[1], proportions=[1])
        assert event._parameters() == "derived=0 ancestral=[1] proportions=[1.00]"
        assert event._effect() == (
            "Moves all lineages from admixed population '0' to ancestral population. "
            "Lineages move to '1' with proba 1. Set '0' to inactive, and all "
            "migration rates to and from '0' to zero."
        )
        event = msprime.Admixture(
            time=1, derived=0, ancestral=[1, 2], proportions=[1 / 4, 3 / 4]
        )
        assert (
            event._parameters() == "derived=0 ancestral=[1, 2] proportions=[0.25, 0.75]"
        )
        assert event._effect() == (
            "Moves all lineages from admixed population '0' to ancestral populations. "
            "Lineages move to '1' with proba 0.25; '2' with proba 0.75. "
            "Set '0' to inactive, and all migration rates to and from '0' to zero."
        )

    def test_mass_migration(self):
        event = msprime.MassMigration(time=1, proportion=0.5, source=0, dest=1)
        event._parameters() == "source=0, dest=1, proportion=0.5"
        effect = (
            "Lineages currently in population 0 move to 1 with probability 0.5 "
            "(equivalent to individuals migrating from 1 to 0 forwards in time)"
        )
        assert event._effect() == effect

        event = msprime.MassMigration(time=1, proportion=1.0, source=0, dest=1)
        event._parameters() == "source=0, dest=1, proportion=1.0"
        effect = (
            "All lineages currently in population 0 move to 1 "
            "(equivalent to individuals migrating from 1 to 0 forwards in time)"
        )
        assert event._effect() == effect

    def test_simple_bottleneck(self):
        event = msprime.SimpleBottleneck(time=1, population=1, proportion=0.5)
        assert event._parameters() == "population=1, proportion=0.5"
        assert event._effect() == (
            "Lineages in population 1 coalesce with probability 0.5"
        )

    def test_instantaneous_bottleneck(self):
        event = msprime.InstantaneousBottleneck(time=1, population=1, strength=1.5)
        assert event._parameters() == "population=1, strength=1.5"
        assert event._effect() == "Equivalent to 1.5 generations of the coalescent"

    def test_census(self):
        event = msprime.CensusEvent(time=1)
        assert event._parameters() == ""
        assert event._effect() == (
            "Insert census nodes to record the location of all lineages"
        )


class DebugOutputBase:
    def test_zero_samples_old_style(self):
        population_configurations = [msprime.PopulationConfiguration()]
        self.verify(msprime.Demography.from_old_style(population_configurations))

    def test_one_population(self):
        demography = msprime.Demography.isolated_model([10])
        demography.add_population_parameters_change(0.1, initial_size=2),
        demography.add_population_parameters_change(0.1, growth_rate=10),
        self.verify(demography)

    def test_no_events(self):
        demography = msprime.Demography.isolated_model([10, 11])
        self.verify(demography)

    def test_island_model(self):
        demography = msprime.Demography.island_model([100, 100], migration_rate=0.1)
        self.verify(demography)

    def test_primates_species_tree(self):
        demography = msprime.Demography.from_species_tree(
            "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)",
            initial_size=10_000,
            time_units="myr",
            generation_time=28,
        )
        self.verify(demography)

    def test_all_events(self):
        demography = all_events_example_demography(integer_ids=True)
        self.verify(demography)

    def test_all_events_string_names(self):
        demography = all_events_example_demography()
        self.verify(demography)

    def test_ooa_model(self):
        self.verify(msprime.Demography._ooa_model())

    def test_american_admixture_model(self):
        self.verify(msprime.Demography._american_admixture_model())


class TestDemographyHtml(DebugOutputBase):
    def verify(self, demography):
        html = demography._repr_html_()
        root = xml.etree.ElementTree.fromstring(html)
        assert root.tag == "div"
        children = list(root)
        assert len(children) == 3
        for child in children:
            assert list(child)[2].tag == "table"
        pop_table = list(children[0])[2]
        rows = list(pop_table.find("tbody"))
        assert len(rows) == demography.num_populations
        migration_matrix_table = list(children[1])[2]
        rows = list(migration_matrix_table.find("tbody"))
        assert len(rows) == demography.num_populations
        for row in rows:
            assert len(row) == demography.num_populations + 1
        events_table = list(children[2])[2]
        rows = list(events_table.find("tbody"))
        assert len(rows) == len(demography.events)
        # TODO add more tests when the output format is finalised.


class TestDemographyDebuggerHtml(DebugOutputBase):
    def verify(self, demography):
        debugger = demography.debug()
        html = debugger._repr_html_()
        root = xml.etree.ElementTree.fromstring(html)
        assert root.tag == "div"
        root_children = list(root)
        assert len(root_children) == len(debugger.epochs)
        for details, epoch in zip(root_children, debugger.epochs):
            assert details.tag == "details"
            children = list(details)
            assert children[0].tag == "summary"
            assert children[1].tag == "div"
            assert children[0].text == epoch._title_text()


class TestDemographyText(DebugOutputBase):
    def verify(self, demography):
        text = str(demography)
        assert text.startswith("Demography")
        sections = text.split("╟")
        assert len(sections) == 4
        assert sections[0].startswith("Demography")
        sec_demography = sections[1].splitlines()
        assert len(sec_demography) == demography.num_populations + 5
        sec_migration_matrix = sections[2].splitlines()
        assert len(sec_migration_matrix) == demography.num_populations + 5
        # We don't really know how many lines will be in the events section.


class TestDemographyDebuggerText(DebugOutputBase):
    def verify(self, demography):
        debugger = demography.debug()
        text = str(debugger)
        assert text.startswith("DemographyDebugger")
        sections = text.split("Epoch")
        assert len(sections) - 1 == len(debugger.epochs)
        # Check the deprecated interface.
        buff = io.StringIO()
        debugger.print_history(buff)
        assert text == buff.getvalue()


class TestDemographyTextExamples:
    def test_one_population(self):
        demography = msprime.Demography.isolated_model([10])
        out = textwrap.dedent(
            """\
        Demography
        ╟  Populations
        ║  ┌────────────────────────────────────────────────────────────────────────────────────────┐
        ║  │ id │name   │description  │initial_size  │ growth_rate │  sampling_time│extra_metadata  │
        ║  ├────────────────────────────────────────────────────────────────────────────────────────┤
        ║  │ 0  │pop_0  │             │10.0          │    0.00     │              0│{}              │
        ║  └────────────────────────────────────────────────────────────────────────────────────────┘
        ╟  Migration Matrix
        ║  ┌───────────────┐
        ║  │       │ pop_0 │
        ║  ├───────────────┤
        ║  │  pop_0│   0   │
        ║  └───────────────┘
        ╟  Events
        ║  ┌───────────────────────────────────┐
        ║  │  time│type  │parameters  │effect  │
        ║  ├───────────────────────────────────┤
        ║  └───────────────────────────────────┘
        """  # noqa: B950
        )
        assert out == str(demography)

    def test_two_populations(self):
        demography = msprime.Demography.isolated_model([10, 20], growth_rate=[1, 2])
        demography.migration_matrix[0, 1] = 0.1
        demography.migration_matrix[1, 0] = 0.2
        out = textwrap.dedent(
            """\
        Demography
        ╟  Populations
        ║  ┌────────────────────────────────────────────────────────────────────────────────────────┐
        ║  │ id │name   │description  │initial_size  │ growth_rate │  sampling_time│extra_metadata  │
        ║  ├────────────────────────────────────────────────────────────────────────────────────────┤
        ║  │ 0  │pop_0  │             │10.0          │    1.00     │              0│{}              │
        ║  │ 1  │pop_1  │             │20.0          │    2.00     │              0│{}              │
        ║  └────────────────────────────────────────────────────────────────────────────────────────┘
        ╟  Migration Matrix
        ║  ┌───────────────────────┐
        ║  │       │ pop_0 │ pop_1 │
        ║  ├───────────────────────┤
        ║  │  pop_0│   0   │  0.1  │
        ║  │  pop_1│  0.2  │   0   │
        ║  └───────────────────────┘
        ╟  Events
        ║  ┌───────────────────────────────────┐
        ║  │  time│type  │parameters  │effect  │
        ║  ├───────────────────────────────────┤
        ║  └───────────────────────────────────┘
        """  # noqa: B950
        )
        assert out == str(demography)

    @pytest.mark.skip("Skipping until events are stable.")
    def test_all_events(self):
        demography = all_events_example_demography(integer_ids=True)

        out = textwrap.dedent(
            """\
        Demography
        ╟  Populations
        ║  ┌────────────────────────────────────────────────────────────────────────┐
        ║  │ id │name   │description  │initial_size  │ growth_rate │extra_metadata  │
        ║  ├────────────────────────────────────────────────────────────────────────┤
        ║  │ 0  │pop_0  │             │1.0           │     0.0     │{}              │
        ║  │ 1  │pop_1  │             │1.0           │     0.0     │{}              │
        ║  └────────────────────────────────────────────────────────────────────────┘
        ╟  Migration Matrix
        ║  ┌───────────────────────┐
        ║  │       │ pop_0 │ pop_1 │
        ║  ├───────────────────────┤
        ║  │  pop_0│   0   │   0   │
        ║  │  pop_1│   0   │   0   │
        ║  └───────────────────────┘
        ╟  Events
        ║  ┌──────────────────────────────────────────────────────────────────────────────────────┐
        ║  │  time│type            │parameters           │effect                                  │
        ║  ├──────────────────────────────────────────────────────────────────────────────────────┤
        ║  │   0.1│Population      │population=-1,       │initial_size → 2 for all populations    │
        ║  │      │parameter       │initial_size=2       │                                        │
        ║  │      │change          │                     │                                        │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │   0.1│Population      │population=-1,       │growth_rate → 10 for all populations    │
        ║  │      │parameter       │growth_rate=10       │                                        │
        ║  │      │change          │                     │                                        │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │   0.1│Population      │population=-1,       │initial_size → 1 and growth_rate → 10   │
        ║  │      │parameter       │initial_size=1,      │for all populations                     │
        ║  │      │change          │growth_rate=10       │                                        │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │   0.2│Migration rate  │source=0, dest=1,    │Backwards-time migration rate from 0    │
        ║  │      │change          │rate=1               │to 1 → 1                                │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │   0.2│Migration rate  │source=1, dest=0,    │Backwards-time migration rate from 1    │
        ║  │      │change          │rate=1               │to 0 → 1                                │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │   0.4│Mass Migration  │source=1, dest=0,    │Lineages currently in population 1      │
        ║  │      │                │proportion=0.9       │move to 0 with probability 0.9          │
        ║  │      │                │                     │(equivalent to individuals migrating    │
        ║  │      │                │                     │from 0 to 1 forwards in time)           │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │   0.4│Migration rate  │source=-1, dest=-1,  │Backwards-time migration rate for all   │
        ║  │      │change          │rate=0               │populations → 0                         │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │   0.5│Instantaneous   │population=0,        │Equivalent to 100 generations of the    │
        ║  │      │Bottleneck      │strength=100         │coalescent                              │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │  0.55│Census          │                     │Insert census nodes to record the       │
        ║  │      │                │                     │location of all lineages                │
        ║  │┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈│
        ║  │  0.56│Simple          │population=1,        │Lineages in population 1 coalesce with  │
        ║  │      │Bottleneck      │proportion=0.1       │probability 0.1                         │
        ║  └──────────────────────────────────────────────────────────────────────────────────────┘
        """  # noqa: B950
        )
        assert out == str(demography)


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
            assert np.isnan(rA)
            assert pA == 0

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
        assert abs(m - truth) / truth < 5e-3
        # test passing steps
        m = ddb.mean_coalescence_time(num_samples=[2], steps=[0, 100, 1000])
        assert abs(m - truth) / truth < 5e-3

    def test_logging(self):
        ddb = self.one_pop_example()
        with mock.patch("msprime.demography.logger.debug") as mocked_debug:
            ddb.mean_coalescence_time([2])
        assert mocked_debug.call_count == 3

    def test_mean_errors(self):
        ddb = self.one_pop_example()
        with pytest.raises(ValueError):
            ddb.mean_coalescence_time(
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
        assert abs(mA - tA) / tA < 5e-3
        assert abs(mB - tB) / tB < 5e-3
        assert abs(mC - tC) / tC < 5e-3

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
        assert max(np.abs(rates - rates2[::2])) < 1e-6
        assert max(np.abs(P - P2[::2])) < 1e-6

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
        with pytest.warns(UserWarning):
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
        with pytest.raises(ValueError):
            ddb.coalescence_rate_trajectory(steps=steps, num_samples=[2, 0])
        # Test that when steps are not strictly increasing values, we throw error.
        with pytest.raises(ValueError):
            ddb.coalescence_rate_trajectory(
                steps=np.flip(steps, axis=0), num_samples=[2]
            )
        # Test that when steps are negative, we throw error
        with pytest.raises(ValueError):
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

    @pytest.mark.slow
    def test_random_example(self):
        ddb = self.get_random_example()
        num_samples = list(range(ddb.num_populations))
        rates, P = ddb.coalescence_rate_trajectory(
            steps=np.linspace(0, 200, 2001), num_samples=num_samples
        )
        assert np.all(rates >= 0)
        assert np.all(P >= 0)
        assert np.all(P <= 1)
        assert np.all(np.diff(P) <= 0)
        coaltime = ddb.mean_coalescence_time(num_samples=num_samples)
        assert coaltime > 0
        coaltime2 = ddb.mean_coalescence_time(
            num_samples=num_samples, steps=np.linspace(0, 200, 501)
        )
        assert abs(coaltime - coaltime2) < 2


class TestMatrixExponential:
    """
    Test cases for the matrix exponential function.
    """

    def verify(self, A):
        E1 = scipy.linalg.expm(A)
        E2 = msprime.demography._matrix_exponential(A)
        assert E1.shape == E2.shape
        assert np.allclose(E1, E2)

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
            assert A == B

    def test_identity_exp(self):
        # (-1) * np.eye(k), compared to exp(-1) * np.eye(k)
        for k in range(2, 5):
            A = msprime.demography._matrix_exponential((-1) * np.eye(k))
            B = np.exp(-1) * np.eye(k)
            assert np.allclose(A, B)


class TestEventTimes:
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
            assert len(migrations) == 2
            assert migrations[0].time == start_time
            assert migrations[1].time == start_time
            nodes = list(ts.nodes())
            assert nodes[0].population == 0
            assert nodes[1].population == 0
            assert nodes[2].population == 1

    def test_negative_times(self):
        with pytest.raises(ValueError):
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
                assert ts.first().num_roots == 1


class TestCoalescenceLocations:
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
        assert tree.root == 2
        assert tree.time(2) > t
        assert tree.population(0) == 0
        assert tree.population(1) == 1
        assert tree.population(2) == 2
        assert ts.node(0).population == 0
        assert ts.node(1).population == 1
        assert list(ts.samples()) == [0, 1]
        assert list(ts.samples(0)) == [0]
        assert list(ts.samples(1)) == [1]
        assert ts.num_populations == 3

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
        assert tree.root == 2 * n - 2
        assert tree.time(tree.root) > t
        for j in range(n // 2):
            assert tree.population(j) == 0
            assert tree.population(n // 2 + j) == 1
            assert ts.get_population(j) == 0
            assert ts.get_population(n // 2 + j) == 1
        assert tree.population(tree.root) == 2

        assert np.array_equal(ts.samples(0), np.arange(n // 2, dtype=np.int32))
        assert np.array_equal(ts.samples(1), np.arange(n // 2, n, dtype=np.int32))
        assert np.array_equal(ts.samples(2), np.array([], dtype=np.int32))
        assert ts.num_populations == 3

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
        assert tree.root == 2 * n - 2
        assert tree.time(tree.root) > t
        for j in range(n // 3):
            assert tree.population(j) == 0
            assert tree.population(n // 3 + j) == 1
            assert tree.population(2 * (n // 3) + j) == 2
            assert ts.get_population(j) == 0
            assert ts.get_population(n // 3 + j) == 1
            assert ts.get_population(2 * (n // 3) + j) == 2
        # The MRCAs of 0, 1 and 3 must have occured in deme 0
        assert tree.population(tree.get_mrca(0, n // 3)) == 0
        assert tree.population(tree.get_mrca(0, 2 * (n // 3))) == 0
        # The MRCAs of all the samples within each deme must have
        # occured within that deme
        for k in range(3):
            deme_samples = range(k * (n // 3), (k + 1) * (n // 3))
            for u, v in itertools.combinations(deme_samples, 2):
                mrca_pop = tree.population(tree.get_mrca(u, v))
                assert k == mrca_pop
        assert np.array_equal(ts.samples(0), np.arange(n // 3, dtype=np.int32))
        assert np.array_equal(
            ts.samples(1), np.arange(n // 3, 2 * (n // 3), dtype=np.int32)
        )
        assert np.array_equal(ts.samples(2), np.arange(2 * (n // 3), n, dtype=np.int32))

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
            assert tree.population(j) == j
            assert ts.get_population(j) == j
            assert ts.samples(j) == [j]
        # The MRCA of 0 and 1 should happen in 1 at time > t1, and < t2
        u = tree.get_mrca(0, 1)
        assert u == 4
        assert tree.population(u) == 1
        g = tree.time(u) * 4
        assert t1 < g < t2
        # The MRCA of 0, 1 and 2 should happen in 2 at time > t2 and < t3
        u = tree.get_mrca(0, 2)
        assert u == 5
        assert tree.population(u) == 2
        assert t2 < tree.time(u) < t3
        # The MRCA of 0, 1, 2 and 3 should happen in 3 at time > t3
        u = tree.get_mrca(0, 3)
        assert u == 6
        assert tree.population(u) == 3
        assert tree.time(u) > t3

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
        assert tree.population(0) == 0
        assert tree.population(1) == 3
        assert ts.node(0).population == 0
        assert ts.node(1).population == 3
        assert list(ts.samples(0)) == [0]
        assert list(ts.samples(1)) == []
        assert list(ts.samples(2)) == []
        assert list(ts.samples(3)) == [1]
        # The MRCA of 0, 1 in 3 at time > t3
        u = tree.get_mrca(0, 1)
        assert u == 2
        assert tree.population(u) == 3
        g = tree.time(u) * 4
        assert g > t3

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
        assert tree.population(0) == 0
        assert tree.population(1) == 3
        assert ts.node(0).population == 0
        assert ts.node(1).population == 3
        assert list(ts.samples(0)) == [0]
        assert list(ts.samples(1)) == []
        assert list(ts.samples(2)) == []
        assert list(ts.samples(3)) == [1]
        # The MRCA of 0, 1 in 3 at time > t3
        u = tree.get_mrca(0, 1)
        assert u == 2
        assert tree.population(u) == 3
        g = tree.time(u) * 4
        assert g > t3
        assert ts.num_populations == 4

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
        assert tree.root == 2
        assert tree.time(2) > t / 4
        assert tree.population(0) == 0
        assert tree.population(1) == 1
        assert tree.population(2) == 2
        assert ts.node(0).population == 0
        assert ts.node(1).population == 1
        assert list(ts.samples()) == [0, 1]
        assert list(ts.samples(0)) == [0]
        assert list(ts.samples(1)) == [1]

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
        assert tree.root == 2
        assert tree.time(2) > t / 4
        assert tree.population(0) == 0
        assert tree.population(1) == 1
        assert tree.population(2) == 2
        assert ts.node(0).population == 0
        assert ts.node(1).population == 1
        assert list(ts.samples()) == [0, 1]
        assert list(ts.samples(0)) == [0]
        assert list(ts.samples(1)) == [1]

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
        assert tree.root == 2
        assert tree.time(2) > t
        assert tree.population(0) == 0
        assert tree.population(1) == num_demes - 1
        assert tree.population(2) == num_demes - 1
        assert ts.node(0).population == 0
        assert ts.node(1).population == num_demes - 1

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
        assert tree.root == 2
        assert tree.time(2) > t
        assert tree.population(0) == 0
        assert tree.population(1) == num_demes - 1
        assert tree.population(2) == num_demes - 1
        assert ts.node(0).population == 0
        assert ts.node(1).population == num_demes - 1

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
        assert tree.time(tree.root) > t4
        assert tree.population(tree.root) == 0
        # The parent of all the samples from each deme should be in that deme.
        for pop in range(3):
            parents = [tree.get_parent(u) for u in ts.samples(population=pop)]
            for v in parents:
                assert tree.population(v) == pop
        assert ts.num_populations == 3


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
        assert ts.num_migrations == len(migrations)
        oldest_t = max(node.time for node in ts.nodes())
        for mig in migrations:
            assert mig.left >= 0
            assert mig.right <= ts.sequence_length
            assert mig.source != mig.dest
            assert 0 <= mig.node < ts.num_nodes
            assert 0 <= mig.time < oldest_t
        # Migrations must be listed in non-decreasing time order.
        for j in range(1, len(migrations)):
            assert migrations[j - 1].time <= migrations[j].time

    def verify_two_pops_single_sample(self, ts, t):
        self.verify_migrations(ts)
        migrations = list(ts.migrations())
        assert len(migrations) == 2
        m0 = migrations[0]
        m1 = migrations[1]
        assert m0.left == 0
        assert m0.right == 1
        assert m0.source == 0
        assert m0.dest == 2
        assert m0.time == t
        assert m1.left == 0
        assert m1.right == 1
        assert m1.source == 1
        assert m1.dest == 2
        assert m1.time == t

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
        assert len(migrations) > 0
        for mig in migrations:
            assert mig.time > 0
            assert ts.node(mig.node).population == 1
            assert mig.source == 1
            assert mig.dest == 0
            assert mig.left == 0
            assert mig.right == 1
        # All lineages should end up in population 1
        for tree in ts.trees():
            node = ts.node(tree.root)
            assert node.population == 0

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
        assert ts.num_trees > 1
        migrations = list(ts.migrations())
        assert len(migrations) > 0
        for mig in migrations:
            assert mig.time > 0
            assert ts.node(mig.node).population == 1
            assert mig.source == 1
            assert mig.dest == 0
            assert mig.left >= 0
            assert mig.right <= 1

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
        assert ts.num_trees > 10
        migrations = list(ts.migrations())
        assert len(migrations) > 0
        for mig in migrations:
            assert mig.time > 0
            assert mig.node == 1
            assert mig.source == 1
            assert mig.dest == 0
            assert mig.left >= 0
            assert mig.right <= 1


class TestMigrationRecordsHudson(MigrationRecordsMixin):
    model = "hudson"


class TestMigrationRecordsSmc(MigrationRecordsMixin):
    model = "smc"


class TestMigrationRecordsSmcPrime(MigrationRecordsMixin):
    model = "smc_prime"


class TestMigrationRecordsDtwf(MigrationRecordsMixin):
    model = "dtwf"


class TestFullArgMigration:
    """
    Tests for migration with the full ARG.
    """

    def verify_two_pops_full_arg(self, ts):
        migrations = ts.tables.migrations
        edges = ts.tables.edges
        nodes = ts.tables.nodes
        for mig in migrations:
            assert nodes[mig.node].flags == msprime.NODE_IS_MIG_EVENT
            assert nodes[mig.node].time == mig.time
            assert nodes[mig.node].population == mig.dest
            e1 = np.where(edges.parent == mig.node)
            e2 = np.where(edges.left == mig.left)
            e = np.intersect1d(e1[0], e2[0])
            assert len(e) == 1
            e = e.item()
            assert edges[e].right == mig.right
            assert nodes[edges[e].child].population == mig.source
        for edge in edges:
            if nodes[edge.parent].flags == msprime.NODE_IS_MIG_EVENT:
                m1 = np.where(migrations.node == edge.parent)
                m2 = np.where(migrations.left == edge.left)
                m = np.intersect1d(m1[0], m2[0])
                assert len(m) == 1
                m = m.item()
                assert migrations[m].right == edge.right
                assert migrations[m].time == nodes[edge.parent].time
                assert migrations[m].source == nodes[edge.child].population
                assert migrations[m].dest == nodes[edge.parent].population

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
            assert u == 2
            assert g <= tree.time(u)


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
            assert len(migrations) > 0
            for mr in migrations:
                self.assertAlmostEqual(g, mr.time)
            tree = next(ts.trees())
            u = tree.get_mrca(0, 1)
            assert u == 2
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


class TestTimeUnitsWrightFisher(TimeUnitsMixin):
    model = "dtwf"


class HistoricalSamplingMixin:
    """
    Tests to make sure historical sampling works correctly.
    """

    def test_two_diploid_samples(self):
        N = 100
        sampling_time = 1.01 * N
        ts = msprime.sim_ancestry(
            demography=msprime.Demography.island_model([N, N], migration_rate=1),
            samples=[
                msprime.SampleSet(1, population=0),
                msprime.SampleSet(1, population=1, time=sampling_time),
            ],
            ploidy=2,
            model=self.model,
            random_seed=3,
        )
        for t in ts.trees():
            assert t.get_time(0) == 0
            assert t.get_time(1) == 0
            assert t.get_time(2) == sampling_time
            assert t.get_time(3) == sampling_time

    def test_two_samples(self):
        N = 100
        sampling_time = 1.01 * N
        for recombination_rate in [0, 1]:
            ts = msprime.simulate(
                Ne=N,
                model=self.model,
                recombination_rate=recombination_rate,
                length=1,
                samples=[msprime.Sample(0, 0), msprime.Sample(0, sampling_time)],
                random_seed=3,
            )
            for t in ts.trees():
                assert t.get_time(0) == 0
                assert t.get_time(1) == sampling_time
                assert t.get_parent(0) == t.get_parent(1)
                assert t.get_parent(1) == t.get_parent(0)
                assert t.get_time(t.get_parent(0)) > sampling_time

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
            assert ts.num_nodes == 3
            assert nodes[0].time == 0
            assert nodes[1].time == sampling_time
            assert nodes[2].time > sampling_time
            assert nodes[2].time > start_time

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
        assert t.get_time(0) == 0
        assert t.get_time(1) == st1
        assert t.get_time(2) == st2
        assert t.get_time(3) == st3
        assert t.get_time(t.get_parent(1)) > st1
        assert t.get_time(t.get_parent(2)) > st2
        assert t.get_time(t.get_parent(3)) > st3
        assert t.get_time(t.get_root()) > st3

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
            assert time[j] == sampling_time
        assert time[n - 1] == 0
        # Allow it to be within 10 coalescent time units.
        assert time[-1] < sampling_time + 10 * N

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
                assert time[0] == 0
                assert time[1] == 0
                assert time[2] > start_time
                if offset is None:
                    offset = time[2] - start_time
                else:
                    self.assertAlmostEqual(offset, time[2] - start_time)

    def test_negative_start_time(self):
        ts = msprime.simulate(2, Ne=10, start_time=-1, model=self.model, random_seed=2)
        tables = ts.tables
        assert tables.nodes[0].time == 0
        assert tables.nodes[1].time == 0
        assert len(tables.edges) == 2

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
            assert tables.nodes[0].time == 10
            assert tables.nodes[1].time == 10
            assert len(tables.edges) == 2

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
        assert t.get_time(0) == 0
        assert t.get_time(1) == sampling_time
        assert t.get_time(2) >= migration_time
        assert t.get_population(0) == 0
        assert t.get_population(1) == 1
        assert t.get_population(2) == 0

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
        assert len(tables.edges) == 2
        assert len(tables.migrations) > 1

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
        assert t.get_time(0) == 0
        assert t.get_time(1) == t1
        assert t.get_time(2) == t2
        assert t.get_time(3) == t3
        assert t.get_population(0) == 0
        assert t.get_population(1) == 1
        assert t.get_population(2) == 2
        assert t.get_population(3) == 3
        assert t.get_population(4) == 1
        assert t.get_population(5) == 2
        assert t.get_population(6) == 3
        assert t1 < t.get_time(4) < t2
        assert t2 < t.get_time(5) < t3
        assert t.get_time(6) > t3


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
                assert time[0] == sampling_time
                assert time[1] == 0
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
        assert time[0] == 0
        assert time[1] == 1.1
        assert time[2] == 1.2


class EndTimeMixin:
    """
    Tests for the max_time parameter.
    """

    def verify_empty_tree_sequence(self, n, ts):
        assert ts.num_edges == 0
        assert ts.num_trees == 1
        assert ts.num_nodes == n
        assert ts.num_samples == n
        tree = ts.first()
        assert tree.num_roots == n

    def verify_incomplete_tree_sequence(self, n, max_time, ts):
        assert ts.num_samples == n
        time = ts.tables.nodes.time
        for tree in ts.trees():
            # Every sample with time <= max_time will end on a path
            # with time == max_time
            for u in tree.samples():
                if time[u] <= max_time:
                    while tree.parent(u) != tskit.NULL:
                        u = tree.parent(u)
                    assert ts.node(u).time == max_time
                else:
                    assert tree.parent(u) == tskit.NULL
        max_roots = max(tree.num_roots for tree in ts.trees())
        assert max_roots > 1

    def test_zero_time(self):
        n = 10
        for n in [2, 10, 100]:
            ts = msprime.simulate(n, end_time=0, model=self.model)
            self.verify_empty_tree_sequence(n, ts)

    def test_negative(self):
        with pytest.raises(ValueError):
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
        assert tables1 == tables2

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
        assert ts.num_migrations > 0
        assert np.all(ts.tables.migrations.time < max_time)

    def test_ancient_samples(self):
        n = 40
        samples = [msprime.Sample(time=j, population=0) for j in range(n)]
        max_time = 20
        ts = msprime.simulate(
            samples=samples, Ne=10, end_time=max_time, model=self.model, random_seed=100
        )
        self.verify_incomplete_tree_sequence(n, max_time, ts)
        nodes = ts.tables.nodes
        assert np.array_equal(nodes.time[:n], np.arange(n))
        assert len(nodes) > n
        tree = ts.first()
        assert tree.num_roots > 1

    def test_all_ancient_samples(self):
        n = 40
        samples = [msprime.Sample(time=j + 1, population=0) for j in range(n)]
        max_time = 20
        ts = msprime.simulate(
            samples=samples, Ne=10, end_time=max_time, model=self.model, random_seed=100
        )
        self.verify_incomplete_tree_sequence(n, max_time, ts)
        nodes = ts.tables.nodes
        assert np.array_equal(nodes.time[:n], np.arange(n) + 1)
        assert len(nodes) > n
        tree = ts.first()
        assert tree.num_roots > 1

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
        assert ts.num_nodes == 3
        assert ts.num_edges == 1

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


class TestEndTimeHudson(EndTimeMixin):
    model = "hudson"


class TestEndTimeWrightFisher(EndTimeMixin):
    model = "dtwf"


class TestEventsBetweenGenerationsWrightFisher:
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
            assert node.time == int(node.time)


class TestOldStylePopulationMetadata:
    """
    Tests for the metadata behaviour on PopulationConfiguration objects.
    """

    def test_simple_case(self):
        md = {"x": "y"}
        ts = msprime.simulate(
            population_configurations=[msprime.PopulationConfiguration(2, metadata=md)],
            random_seed=1,
        )
        assert ts.num_populations == 1
        pop = ts.population(0)
        assert md == json.loads(pop.metadata)

    def test_old_style_metadata_via_sim_ancestry(self):
        md = {"x": "y"}
        demography = msprime.Demography.from_old_style(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=1, metadata=md)
            ]
        )
        ts = msprime.sim_ancestry(2, demography=demography, random_seed=1)
        assert ts.num_populations == 1
        pop = ts.population(0)
        expected = {
            "name": "pop_0",
            "description": "",
            **md,
        }
        assert expected == pop.metadata

    def test_old_style_metadata_description_is_merged(self):
        md = {"description": "y"}
        demography = msprime.Demography.from_old_style(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=1, metadata=md)
            ]
        )
        ts = msprime.sim_ancestry(2, demography=demography, random_seed=1)
        assert ts.num_populations == 1
        pop = ts.population(0)
        expected = {
            "name": "pop_0",
            "description": "y",
        }
        assert expected == pop.metadata

    def test_old_style_metadata_name_conflicts(self):
        demography = msprime.Demography()
        md = {"name": "x"}
        pc = msprime.PopulationConfiguration(initial_size=1, metadata=md)
        with pytest.raises(ValueError, match="doesn't match"):
            demography._add_population_from_old_style(pc, name="y")

    def test_old_style_metadata_name(self):
        demography = msprime.Demography()
        md = {"name": "x"}
        pc = msprime.PopulationConfiguration(initial_size=1, metadata=md)
        pop = demography._add_population_from_old_style(pc)
        assert pop.name == "x"

        pc = msprime.PopulationConfiguration(initial_size=1)
        pop = demography._add_population_from_old_style(pc, name="y")
        assert pop.name == "y"

    def test_old_style_metadata_name_is_merged(self):
        md = {"name": "y"}
        demography = msprime.Demography.from_old_style(
            population_configurations=[
                msprime.PopulationConfiguration(initial_size=1, metadata=md)
            ]
        )
        ts = msprime.sim_ancestry(2, demography=demography, random_seed=1)
        assert ts.num_populations == 1
        pop = ts.population(0)
        expected = {
            "name": "y",
            "description": "",
        }
        assert expected == pop.metadata

    def test_default(self):
        ts = msprime.simulate(
            population_configurations=[msprime.PopulationConfiguration(2)],
            random_seed=1,
        )

        assert ts.num_populations == 1
        pop = ts.population(0)
        assert pop.metadata == b""

        # An explicit value of None also gives us an empty user metadata dict.
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(2, metadata=None)
            ],
            random_seed=1,
        )
        assert ts.num_populations == 1
        pop = ts.population(0)
        assert pop.metadata == b""

    def test_errors(self):
        for bad_metadata in [b"asdf", Exception]:
            pop_conf = msprime.PopulationConfiguration(2, metadata=bad_metadata)
            with pytest.raises(TypeError):
                msprime.simulate(population_configurations=[pop_conf])

    def test_multi_population(self):
        for num_pops in range(1, 10):
            pop_configs = [
                msprime.PopulationConfiguration(2, metadata={"x": "x" * j})
                for j in range(num_pops)
            ]
            ts = msprime.simulate(
                population_configurations=pop_configs, random_seed=1, end_time=1
            )
            assert ts.num_populations == num_pops
            for j in range(num_pops):
                pop = ts.population(j)
                assert pop_configs[j].metadata == json.loads(pop.metadata)


class TestPopulationMetadata:
    """
    Tests for the metadata behaviour on Population objects.
    """

    def test_defaults(self):
        demography = msprime.Demography.isolated_model([10])
        ts = msprime.sim_ancestry(2, demography=demography, random_seed=2)
        assert ts.num_populations == 1
        metadata = ts.population(0).metadata
        assert metadata["name"] == "pop_0"
        assert metadata["description"] == ""
        assert len(metadata) == 2

    def test_extra_metadata(self):
        demography = msprime.Demography.isolated_model([10])
        md = {"x": "y", "z": "z", "abc": 1234}
        demography.populations[0].extra_metadata = md
        ts = msprime.sim_ancestry(2, demography=demography, random_seed=2)
        assert ts.num_populations == 1
        metadata = ts.population(0).metadata
        assert metadata["name"] == "pop_0"
        assert metadata["description"] == ""
        assert len(metadata) == len(md) + 2
        for key, value in md.items():
            assert metadata[key] == value

    def test_extra_metadata_overwrite_standard(self):
        demography = msprime.Demography.isolated_model([10])
        md = {"x": 1234, "name": "abc"}
        demography.populations[0].extra_metadata = md
        with pytest.raises(ValueError) as record:
            msprime.sim_ancestry(2, demography=demography, random_seed=2)
        assert str(record.value).startswith(
            "Cannot set standard metadata key(s) ['name']"
        )

        md = {"x": 1234, "name": "abc", "description": "sdf"}
        demography.populations[0].extra_metadata = md
        with pytest.raises(ValueError) as record:
            msprime.sim_ancestry(2, demography=demography, random_seed=2)
        assert str(record.value).startswith(
            "Cannot set standard metadata key(s) ['description', 'name']"
        )

    def test_island_model(self):
        demography = msprime.Demography.island_model([1] * 10, migration_rate=0)
        for j, pop in enumerate(demography.populations):
            pop.extra_metadata = {"x": "y", "z": "z" * j, "abc": 1234}
            pop.description = f"Island model node {j}"
        ts = msprime.sim_ancestry({0: 2}, demography=demography, random_seed=2)
        assert ts.num_populations == demography.num_populations
        for j in range(demography.num_populations):
            metadata = ts.population(j).metadata
            assert metadata["name"] == f"pop_{j}"
            assert metadata["description"] == demography.populations[j].description
            assert len(metadata) == 2 + len(pop.extra_metadata)


class TestCensusEvent:
    """
    Tests of the census demographic event.
    """

    def verify(self, ts, census_time):
        """
        Verifies that a census event has been added correctly.
        """
        census_ids = np.where(ts.tables.nodes.flags == msprime.NODE_IS_CEN_EVENT)[0]
        for u in census_ids:
            assert ts.tables.nodes.time[u] == census_time
        assert len(census_ids) > 1
        # Check that all samples have a census ancestor on each tree.
        for tree in ts.trees():
            leaves = []
            census_nodes = [u for u in census_ids if u in list(tree.nodes())]
            for node in census_nodes:
                assert len(tree.children(node)) == 1
                le = list(tree.leaves(node))
                leaves += le
            leaves.sort()
            assert leaves == [i for i in range(0, ts.num_samples)]

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
                assert nodes.population[row.parent] == nodes.population[row.child]

    def test_census_at_existing_node_time(self):
        with pytest.raises(_msprime.LibraryError):
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
                assert nodes.population[row.parent] == nodes.population[row.child]
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
                assert nodes.population[row.parent] == nodes.population[row.child]

    def test_no_census_nodes_above_root_nodes(self):
        ts = msprime.simulate(sample_size=2, random_seed=525)
        assert all(ts.tables.nodes.flags) != msprime.NODE_IS_CEN_EVENT
        tsc = msprime.simulate(
            sample_size=2,
            random_seed=525,
            demographic_events=[msprime.CensusEvent(time=2000)],
        )
        assert ts.tables.nodes == tsc.tables.nodes


class TestPossibleLineages:
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
        assert len(lineages) == 4
        assert np.all(lineages[(0, 50)] == [1, 1, 1, 1])
        assert np.all(lineages[(50, 100)] == [1, 1, 1, 0])
        assert np.all(lineages[(100, 150)] == [1, 1, 0, 0])
        assert np.all(lineages[(150, np.inf)] == [1, 0, 0, 0])

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
        assert len(lineages) == 1
        assert np.all(lineages[(0, np.inf)] == [1, 1])

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
        assert len(lineages) == 4
        assert np.all(lineages[(0, 50)] == [1, 1])
        assert np.all(lineages[(50, 100)] == [1, 0])
        assert np.all(lineages[(100, 150)] == [1, 1])
        assert np.all(lineages[(150, np.inf)] == [1, 0])

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
        assert len(lineages) == 7
        assert np.all(lineages[(0, 100)] == [1, 0, 1, 1, 0])
        assert np.all(lineages[(100, 200)] == [1, 1, 1, 1, 0])
        assert np.all(lineages[(200, 300)] == [1, 1, 1, 0, 0])
        assert np.all(lineages[(300, 400)] == [1, 0, 1, 0, 0])
        assert np.all(lineages[(400, 600)] == [1, 0, 1, 0, 1])
        assert np.all(lineages[(600, 700)] == [1, 0, 0, 0, 1])
        assert np.all(lineages[(700, np.inf)] == [1, 0, 0, 0, 0])


class TestLineageProbabilities:
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
        samples = []
        for pop_id in range(dd.num_populations):
            for time in dd.epoch_times:
                samples.append(
                    msprime.SampleSet(1, population=pop_id, ploidy=1, time=time)
                )
        reps = msprime.sim_ancestry(
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
                    assert probs[j, pop, n.population] > 0.0
                    for epoch in locs:
                        if n.time >= epoch[0] and n.time < epoch[1]:
                            assert locs[epoch][n.population]

    def test_two_pop(self):
        for _b in [2, 0]:
            dd, f = self.two_pop_example(1, 2)
            times = np.linspace(0, 10, 21)
            for st in [0.0, 2.5]:
                P = dd.lineage_probabilities(times + st, sample_time=st)
                for j, t in enumerate(times):
                    assert np.allclose(P[j, :, :], f(t))
            self.verify_simulation(dd)

    @pytest.mark.slow
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
        assert np.all([np.sum(P) == len(pop_config) for P in P_out])
        assert np.all(np.diag(P_out[0]) == [1, 1, 1, 1])
        assert np.all(np.diag(P_out[1]) == [1, 1, 1, 1])
        assert np.all(probs == [1, 0, 0, 0] for probs in P_out[5])
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
        assert np.allclose(P_out[0], [[1, 0], [f_pulse, 1 - f_pulse]])
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
        assert np.all(P_out[0] == np.eye(len(pop_config)))
        assert np.all(P_out[1] > 0)
        assert np.all(P_out[1] > 0)
        # checking if close because of precision of _matrix_exponential function
        assert np.all(np.isclose(P_out[3], [[1, 0], [1, 0]]))
        self.verify_simulation(dd)

        mig_mat = [[0, 0.01], [0, 0]]
        dd = msprime.DemographyDebugger(
            demographic_events=dem_events,
            population_configurations=pop_config,
            migration_matrix=mig_mat,
        )
        P_out = dd.lineage_probabilities([0, 50, 100, 150])
        assert np.all(P_out[0] == np.eye(len(pop_config)))
        assert abs(P_out[1][1][0]) < np.finfo(float).eps
        assert abs(P_out[2][1][0]) < np.finfo(float).eps
        # machine precision instead of zero because of _matrix_exponential function
        assert np.all(np.isclose(P_out[3], [[1, 0], [1, 0]]))
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
        assert np.all(P_out[0] == 0)
        P_out = dd.lineage_probabilities([0, 1, 2], sample_time=1)
        assert np.all(P_out[0] == 0)
        assert np.all(P_out[1] == [[1, 0], [0, 1]])
        assert np.allclose(np.sum(P_out[2], axis=1), 1)
        P_out = dd.lineage_probabilities([99, 100, 101], sample_time=1)
        assert np.all(P_out[0] > 0)
        assert np.all(P_out[1] > 0)
        assert np.allclose(P_out[2], [[1, 0], [1, 0]])
        self.verify_simulation(dd)


class TestPreCannedModels:
    """
    Tests for the specialised models returned by static methods
    on the Demography.
    """

    def assertZeroDiagonal(self, A):
        assert np.all(np.diagonal(A)) == 0


class TestIslandModel(TestPreCannedModels):
    def test_errors(self):
        for bad_N in [-1, -1e6]:
            with pytest.raises(ValueError):
                msprime.Demography.island_model([bad_N], 0.1)
        for bad_m in [-1, -1e5]:
            with pytest.raises(ValueError):
                msprime.Demography.island_model([1], bad_m)

    def test_one_pop(self):
        model = msprime.Demography.island_model([1], migration_rate=1)
        assert len(model.populations) == 1
        assert len(model.migration_matrix) == 1
        ts = msprime.sim_ancestry(samples={0: 2}, demography=model, random_seed=1)
        assert ts.num_populations == 1

    def test_migration(self):
        for N in [1, 2, 5]:
            model = msprime.Demography.island_model([1] * N, 0.1)
            assert len(model.populations) == N
            assert model.migration_matrix.shape == (N, N)
            self.assertZeroDiagonal(model.migration_matrix)
            assert np.all(model.migration_matrix[~np.eye(N, dtype=bool)] == 0.1)
            ts = msprime.sim_ancestry(
                samples={j: 1 for j in range(N)}, demography=model, random_seed=1
            )
            assert ts.num_populations == N
            assert ts.num_samples == 2 * N

    def test_initial_size(self):
        for Ne in [0.1, 1, 10]:
            model = msprime.Demography.island_model([Ne, Ne], 0.1)
            assert model.populations[0].initial_size == Ne
            assert model.populations[1].initial_size == Ne

    def test_migration_zero(self):
        initial_size = [1, 2, 3, 4]
        growth_rate = [0.1, 0.2, 0.3, 0.4]
        model1 = msprime.Demography.island_model(
            initial_size, growth_rate=growth_rate, migration_rate=0
        )
        model2 = msprime.Demography.isolated_model(
            initial_size, growth_rate=growth_rate
        )
        assert model1 == model2


class TestSteppingStoneModel(TestPreCannedModels):
    def test_errors(self):
        for bad_N in [-1, -1000, np.nan]:
            with pytest.raises(ValueError):
                msprime.Demography.stepping_stone_model([bad_N], 0.1)
        for bad_m in [-1, -1e5]:
            with pytest.raises(ValueError):
                msprime.Demography.stepping_stone_model([1], bad_m)

    def test_migration_zero(self):
        initial_size = [1, 2, 3, 4]
        growth_rate = [0.1, 0.2, 0.3, 0.4]
        model1 = msprime.Demography.stepping_stone_model(
            initial_size, growth_rate=growth_rate, migration_rate=0
        )
        model2 = msprime.Demography.isolated_model(
            initial_size, growth_rate=growth_rate
        )
        assert model1 == model2

    def test_one_pop(self):
        for boundaries in [True, False]:
            model = msprime.Demography.stepping_stone_model(
                [1], 1, boundaries=boundaries
            )
            assert len(model.populations) == 1
            assert len(model.migration_matrix) == 1
            ts = msprime.sim_ancestry(samples={0: 2}, demography=model, random_seed=1)
            assert ts.num_populations == 1

    def test_migration_circular(self):
        m = 0.3
        for N in [2, 3, 5]:
            model = msprime.Demography.stepping_stone_model([1] * N, m)
            # Circular is the default
            assert model == msprime.Demography.stepping_stone_model(
                [1] * N, m, boundaries=False
            )
            assert len(model.populations) == N
            assert model.migration_matrix.shape == (N, N)
            self.assertZeroDiagonal(model.migration_matrix)
            for j in range(N):
                adjacent = [(j - 1) % N, (j + 1) % N]
                for k in range(N):
                    if k in adjacent:
                        assert model.migration_matrix[j, k] == m
                    else:
                        assert model.migration_matrix[j, k] == 0
            ts = msprime.sim_ancestry(
                samples={j: 2 for j in range(N)}, demography=model, random_seed=1
            )
            assert ts.num_populations == N
            assert ts.num_samples == 4 * N

    def test_migration_line_two_pops(self):
        m = 1
        model = msprime.Demography.stepping_stone_model([1, 1], m, boundaries=True)
        assert len(model.populations) == 2
        assert model.migration_matrix.shape == (2, 2)
        assert np.all(model.migration_matrix == 0)

    def test_migration_line(self):
        m = 0.3
        for N in [3, 4, 5]:
            model = msprime.Demography.stepping_stone_model([1] * N, m, boundaries=True)
            assert len(model.populations) == N
            assert model.migration_matrix.shape == (N, N)
            self.assertZeroDiagonal(model.migration_matrix)
            for j in range(N):
                adjacent = []
                if j > 0:
                    adjacent.append(j - 1)
                if j < N - 1:
                    adjacent.append(j + 1)
                for k in range(N):
                    if k in adjacent:
                        assert model.migration_matrix[j, k] == m
                    else:
                        assert model.migration_matrix[j, k] == 0
            ts = msprime.sim_ancestry(
                samples={j: 2 for j in range(N)}, demography=model, random_seed=1
            )
            assert ts.num_populations == N
            assert ts.num_samples == 4 * N

    def test_Ne(self):
        model = msprime.Demography.stepping_stone_model([1, 1], 0.1)
        assert model.populations[0].initial_size == 1
        assert model.populations[1].initial_size == 1
        for Ne in [0.1, 1, 10]:
            model = msprime.Demography.stepping_stone_model([Ne, Ne], 0.1)
            assert model.populations[0].initial_size == Ne
            assert model.populations[1].initial_size == Ne


class TestDemographicEventBase:
    def test_str_methods_not_implemented(self):
        de = msprime.DemographicEvent(0)
        with pytest.raises(NotImplementedError):
            de._parameters()
        with pytest.raises(NotImplementedError):
            de._effect()


class TestDemographyObject:
    """
    Basic tests for the demography object.
    """

    def test_equality(self):
        m1 = msprime.Demography.island_model([1, 1], 1 / 3)
        m2 = msprime.Demography.island_model([1, 1], 1 / 3)
        assert m1 == m2
        assert m2 == m1
        assert m1 == m1
        assert not (m1 != m2)
        assert not (m1 != m1)
        assert m1 is not None
        assert m1 != []
        # Validation fills out defaults, so won't be strictly equal
        assert m1.validate() != m1
        assert m1.validate() == m1.validate()
        assert m1.validate() == m2.validate()

        m3 = msprime.Demography.island_model([1, 1], 1 / 3 + 0.001)
        assert m1 != m3
        assert m3 != m1

        assert m1 != msprime.Demography.isolated_model([1])
        assert m1 != msprime.Demography.isolated_model([1, 1])
        assert m1 != msprime.Demography.island_model([2, 1], 1 / 3)

        m1.add_event(msprime.SymmetricMigrationRateChange(1, [0, 1], 0.1))
        assert m1 != m2
        m2.add_event(msprime.SymmetricMigrationRateChange(1, [0, 1], 0.1))
        assert m1 == m2
        m1.events[0].rate = 0.01
        assert m1 != m2

    def test_debug(self):
        model = msprime.Demography.island_model([1, 1], 1 / 3)
        dbg1 = model.debug()
        assert dbg1.demography == model.validate()
        dbg2 = msprime.DemographyDebugger(demography=model)
        assert dbg1.demography == dbg2.demography
        assert str(dbg1) == str(dbg2)

    def test_population_name(self):

        demography = msprime.Demography.isolated_model([1])
        assert demography.populations[0].name == "pop_0"

        demography.populations[0].name = None
        with pytest.raises(ValueError) as excinfo:
            demography.validate()
        assert "A population name must be set." in str(excinfo.value)

        for bad_identifier in ["", "x ", "x y"]:
            demography.populations[0].name = bad_identifier
            with pytest.raises(ValueError) as excinfo:
                demography.validate()
            msg = "A population name must be a valid Python identifier"
            assert msg in str(excinfo.value)

    def test_duplicate_population_name(self):
        demography = msprime.Demography.isolated_model([1, 1])
        demography.populations[1].name = "pop_0"
        with pytest.raises(ValueError, match="Duplicate population name"):
            demography.validate()

        # It's also an error to do this at initialisation time.
        with pytest.raises(ValueError, match="Duplicate population name"):
            msprime.Demography(
                [msprime.Population(10, name="A"), msprime.Population(11, name="A")]
            )

    def test_duplicate_populations(self):
        pop = msprime.Population(10)
        with pytest.raises(ValueError, match="must be distinct"):
            msprime.Demography([pop] * 2)

    def test_population_ids_set_on_init(self):
        pop = msprime.Population(10, id=0)
        with pytest.raises(ValueError, match="ID should not be set"):
            msprime.Demography([pop])

    def test_add_population_error(self):
        model = msprime.Demography.isolated_model([1])
        with pytest.raises(ValueError, match="Duplicate population name"):
            model.add_population(name="pop_0")

    def test_add_population_properties(self):
        model = msprime.Demography()
        pop = model.add_population(
            initial_size=1234,
            growth_rate=123,
            name="XYZ",
            description="asdf",
            sampling_time=0.1234,
            extra_metadata={"x": "y"},
        )
        assert pop.initial_size == 1234
        assert pop.growth_rate == 123
        assert pop.name == "XYZ"
        assert pop.description == "asdf"
        assert pop.sampling_time == 0.1234
        assert pop.extra_metadata == {"x": "y"}

    def test_sidestepping_add_population(self):
        # It's possible for users to sidestep the add_population method.
        # Better to detect mistakes here than to put in opaque infrastructure
        # to try to prevent it.
        model = msprime.Demography.isolated_model([1])
        model.populations.append(msprime.Population(1))
        with pytest.raises(ValueError, match="Incorrect population ID"):
            model.validate()
        model = msprime.Demography.isolated_model([1])
        model.populations.append(msprime.Population(1, id=1, name="pop_0"))
        with pytest.raises(ValueError, match="Duplicate population name"):
            model.validate()

    def test_add_population_migration_matrix(self):
        model = msprime.Demography.island_model([1, 1], 0.1)
        M = np.array([[0, 0.1], [0.1, 0]])
        assert np.array_equal(model.migration_matrix, M)
        model.add_population(initial_size=1)
        M = np.array([[0, 0.1, 0], [0.1, 0, 0], [0, 0, 0]])
        assert np.array_equal(model.migration_matrix, M)

        model = msprime.Demography.island_model([1, 1, 1], 0.1)
        M = np.array([[0, 0.1, 0.1], [0.1, 0, 0.1], [0.1, 0.1, 0]])
        assert np.array_equal(model.migration_matrix, M)
        model.add_population(initial_size=1)
        M = np.array(
            [[0, 0.1, 0.1, 0], [0.1, 0, 0.1, 0], [0.1, 0.1, 0, 0], [0, 0, 0, 0]]
        )
        assert np.array_equal(model.migration_matrix, M)

    def test_add_population_sequential(self):
        model = msprime.Demography()
        for j in range(5):
            model.add_population(initial_size=1)
            model.validate()
            assert model.num_populations == j + 1
            assert np.all(model.migration_matrix == 0)
        for j, pop in enumerate(model.populations):
            assert pop.id == j
            assert pop.name == f"pop_{j}"

    def test_population_lookup(self):
        demography = msprime.Demography.isolated_model([1, 1])
        assert demography.populations[0].name == "pop_0"
        assert demography.populations[1].name == "pop_1"
        for bad_pop in [-1, -2, 2, "x", "AAA"]:
            with pytest.raises(KeyError):
                demography[bad_pop]
            assert bad_pop not in demography
        for bad_type in [b"sdf", 1.0]:
            with pytest.raises(TypeError):
                demography[bad_type]
            with pytest.raises(TypeError):
                bad_type in demography
        # String name lookup works
        assert demography["pop_0"] == demography.populations[0]
        assert "pop_0" in demography
        assert demography["pop_1"] == demography.populations[1]
        assert "pop_1" in demography
        # As does integer lookup
        assert demography[0] == demography.populations[0]
        assert 0 in demography
        assert demography[1] == demography.populations[1]
        assert 1 in demography
        # Numpy integer types are OK too.
        int_array = np.array([0, 1], dtype=np.int8)
        assert demography[int_array[0]] == demography.populations[0]
        assert demography[int_array[1]] == demography.populations[1]

    def test_set_migration_rate(self):
        demography = msprime.Demography.isolated_model([10] * 5)
        demography.set_migration_rate(0, 1, 1)
        assert demography.migration_matrix[0, 1] == 1
        demography.migration_matrix[0, 1] = 0
        assert np.all(demography.migration_matrix == 0)
        demography.set_migration_rate("pop_1", "pop_2", 1)
        assert demography.migration_matrix[1, 2] == 1
        demography.migration_matrix[1, 2] = 0
        assert np.all(demography.migration_matrix == 0)

    def test_set_migration_rate_errors(self):
        demography = msprime.Demography.isolated_model([10] * 5)
        with pytest.raises(KeyError):
            demography.set_migration_rate(-1, 0, 0)
        with pytest.raises(KeyError):
            demography.set_migration_rate("pop_100", 0, 0)
        with pytest.raises(ValueError, match="must be different"):
            demography.set_migration_rate(0, 0, 0)
        with pytest.raises(ValueError, match="must be different"):
            demography.set_migration_rate("pop_0", 0, 0)

    def test_set_symmetric_migration_rate_errors(self):
        demography = msprime.Demography.isolated_model([10] * 5)
        with pytest.raises(ValueError, match="at least two"):
            demography.set_symmetric_migration_rate([1], 0.0)
        with pytest.raises(KeyError):
            demography.set_symmetric_migration_rate([-1, 0], 0.0)

    def test_set_symmetric_migration_rate(self):
        m = 0.25
        demography = msprime.Demography.isolated_model([10] * 5)
        demography.set_symmetric_migration_rate([0, 1], m)
        assert demography.migration_matrix[0, 1] == m
        assert demography.migration_matrix[1, 0] == m
        demography.migration_matrix[0, 1] = 0
        demography.migration_matrix[1, 0] = 0
        assert np.all(demography.migration_matrix == 0)

        # Can use names as well.
        demography = msprime.Demography.isolated_model([10] * 5)
        demography.set_symmetric_migration_rate(["pop_0", "pop_1"], m)
        assert demography.migration_matrix[0, 1] == m
        assert demography.migration_matrix[1, 0] == m
        demography.migration_matrix[0, 1] = 0
        demography.migration_matrix[1, 0] = 0
        assert np.all(demography.migration_matrix == 0)

    @pytest.mark.parametrize("N", [2, 3, 7])
    def test_set_symmetric_migration_rate_all(self, N):
        demography1 = msprime.Demography.isolated_model([10] * N)
        demography1.set_symmetric_migration_rate(range(N), 0.1)
        demography2 = msprime.Demography.island_model([10] * N, 0.1)
        assert np.array_equal(
            demography1.migration_matrix, demography2.migration_matrix
        )

    def test_events_out_of_order(self):
        demography = msprime.Demography.isolated_model([10] * 2)
        demography.add_population_parameters_change(2, initial_size=1)
        demography.add_population_parameters_change(1, initial_size=5)
        demography.add_population_parameters_change(1, growth_rate=0.5)
        with pytest.raises(ValueError, match="Events must be time-sorted"):
            demography.debug()
        demography.sort_events()
        dbg = demography.debug()
        assert dbg.num_epochs == 3

        other = msprime.Demography.isolated_model([10] * 2)
        other.add_population_parameters_change(1, initial_size=5)
        other.add_population_parameters_change(1, growth_rate=0.5)
        other.add_population_parameters_change(2, initial_size=1)
        assert other == demography

    def test_sidestepping_add_event(self):
        demography = msprime.Demography.isolated_model([10] * 2)
        demography.events = [None]
        with pytest.raises(TypeError, match="DemographicEvent instances"):
            demography.debug()

    def test_isolated_model(self):
        demography = msprime.Demography.isolated_model([2])
        assert demography.num_populations == 1
        assert demography.num_events == 0
        assert demography.populations[0].initial_size == 2
        assert demography.populations[0].growth_rate == 0

        demography = msprime.Demography.isolated_model([2], growth_rate=[3])
        assert demography.num_populations == 1
        assert demography.populations[0].initial_size == 2
        assert demography.populations[0].growth_rate == 3

        demography = msprime.Demography.isolated_model([5, 6])
        assert demography.num_populations == 2
        assert demography.populations[0].initial_size == 5
        assert demography.populations[0].growth_rate == 0
        assert demography.populations[1].initial_size == 6
        assert demography.populations[1].growth_rate == 0

        demography = msprime.Demography.isolated_model([5, 6], growth_rate=[7, 8])
        assert demography.num_populations == 2
        assert demography.populations[0].initial_size == 5
        assert demography.populations[0].growth_rate == 7
        assert demography.populations[1].initial_size == 6
        assert demography.populations[1].growth_rate == 8

    def test_isolated_model_errors(self):
        with pytest.raises(ValueError):
            msprime.Demography.isolated_model([2, None])
        with pytest.raises(ValueError):
            msprime.Demography.isolated_model([2], growth_rate=[np.nan])
        with pytest.raises(ValueError):
            msprime.Demography.isolated_model(2)
        with pytest.raises(ValueError):
            msprime.Demography.isolated_model([2], growth_rate=3)
        with pytest.raises(ValueError):
            msprime.Demography.isolated_model([[], []])
        with pytest.raises(ValueError):
            msprime.Demography.isolated_model([1], growth_rate=[[], []])
        with pytest.raises(ValueError):
            msprime.Demography.isolated_model([1], growth_rate=[])

    def test_from_species_tree(self):
        # basic checks here - indepth testing in the test_species_tree_parsing.py file.
        demography = msprime.Demography.from_species_tree(
            "(popA:10.0,popB:10.0)", initial_size=1000
        )
        assert isinstance(demography, msprime.Demography)
        assert len(demography.populations) == 3
        assert demography.populations[0].name == "popA"
        assert demography.populations[0].initial_size == 1000
        assert demography.populations[1].name == "popB"
        assert demography.populations[0].initial_size == 1000
        assert demography.populations[2].name == "pop_2"
        assert demography.populations[2].initial_size == 1000
        assert np.all(demography.migration_matrix == 0)
        assert demography.num_events == len(demography.events)
        assert len(demography.events) == 1
        assert demography.events[0].time == 10
        assert demography.events[0].ancestral == "pop_2"
        assert demography.events[0].derived == ["popA", "popB"]

    def test_from_starbeast(self):
        with open("tests/data/species_trees/91genes_species_rev.tre") as f:
            nexus = f.read()
        demography = msprime.Demography.from_starbeast(nexus, 1)
        assert isinstance(demography, msprime.Demography)
        assert demography.populations[0].name == "spc12"

    def test_validate_resolves_defaults(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="C", initial_size=100)
        demography.set_symmetric_migration_rate(["A", "B"], 0.1)
        demography.add_population_split(10, derived=["A", "B"], ancestral="C")
        assert demography["A"].sampling_time is None
        assert demography["B"].sampling_time is None
        assert demography["C"].sampling_time == 10
        assert demography["A"].initially_active is None
        assert demography["B"].initially_active is None
        assert not demography["C"].initially_active
        validated = demography.validate()
        assert validated["A"].sampling_time == 0
        assert validated["B"].sampling_time == 0
        assert validated["C"].sampling_time == 10
        assert validated["A"].initially_active
        assert validated["B"].initially_active
        assert not validated["C"].initially_active


class TestDemographyFromOldStyle:
    """
    Tests the method for creating a demography object from the old
    style population_configurations, migration_matrix and demographic_events
    parameters.
    """

    def test_defaults(self):
        demog = msprime.Demography.from_old_style()
        assert demog.num_populations == 1
        assert list(demog.migration_matrix) == [[0]]
        assert list(demog.events) == []

    def test_pop_configs_defaults(self):
        for n in range(1, 5):
            pop_configs = [msprime.PopulationConfiguration() for _ in range(n)]
            demog = msprime.Demography.from_old_style(
                population_configurations=pop_configs
            )
            assert demog.num_populations == n
            np.testing.assert_array_equal(demog.migration_matrix, np.zeros((n, n)))
            assert list(demog.events) == []

    def test_ignore_sample_size(self):
        pop_configs = [msprime.PopulationConfiguration(sample_size=1)]
        with pytest.raises(ValueError, match="You have specified a `sample_size`"):
            msprime.Demography.from_old_style(population_configurations=pop_configs)
        demog = msprime.Demography.from_old_style(
            population_configurations=pop_configs, ignore_sample_size=True
        )
        assert demog.num_populations == 1

    def test_migration_matrix(self):
        for n in range(1, 5):
            pop_configs = [msprime.PopulationConfiguration() for _ in range(n)]
            M = np.ones((n, n))
            demog = msprime.Demography.from_old_style(
                population_configurations=pop_configs, migration_matrix=M
            )
            assert demog.num_populations == n
            np.testing.assert_array_equal(demog.migration_matrix, M)
            assert list(demog.events) == []

    def test_demographic_events(self):
        events = [
            msprime.PopulationParametersChange(time=j, initial_size=2)
            for j in range(10)
        ]
        demog = msprime.Demography.from_old_style(demographic_events=events)
        assert demog.num_populations == 1
        assert list(demog.migration_matrix) == [[0]]
        assert events == demog.events

    def test_population_defaults(self):
        pop_config = msprime.PopulationConfiguration()
        demog = msprime.Demography.from_old_style([pop_config])
        pop = demog.populations[0]
        assert pop.initial_size == 1
        assert pop_config.growth_rate == pop.growth_rate

    def test_Ne(self):
        pop_config = msprime.PopulationConfiguration()
        demog = msprime.Demography.from_old_style([pop_config], Ne=1234)
        pop = demog.populations[0]
        assert pop.initial_size == 1234
        assert pop_config.growth_rate == pop.growth_rate

    def test_values(self):
        pop_config = msprime.PopulationConfiguration(
            initial_size=1234, growth_rate=5678
        )
        demog = msprime.Demography.from_old_style([pop_config], Ne=1234)
        pop = demog.populations[0]
        assert pop_config.initial_size == pop.initial_size
        assert pop_config.growth_rate == pop.growth_rate


class TestPopulationNamesInEvents:
    """
    Test that using population names in demographic events has the
    same effect as using the IDs.
    """

    def test_mass_migration(self):
        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_mass_migration(1, source=0, dest=1, proportion=1)
        ts1 = msprime.sim_ancestry(
            {0: 1, 1: 1}, demography=demography, random_seed=1234
        )

        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_mass_migration(1, source="pop_0", dest="pop_1", proportion=1)
        ts2 = msprime.sim_ancestry(
            {0: 1, 1: 1}, demography=demography, random_seed=1234
        )
        assert ts1.equals(ts2, ignore_provenance=True)

    def test_admixture(self):
        demography = msprime.Demography.isolated_model([1000, 1000, 1000, 1000])
        demography.add_admixture(1, derived=0, ancestral=[1, 2], proportions=[0.5, 0.5])
        demography.add_population_split(2, derived=[1, 2], ancestral=3)
        ts1 = msprime.sim_ancestry({0: 10}, demography=demography, random_seed=1234)

        demography = msprime.Demography.isolated_model([1000, 1000, 1000, 1000])
        demography.add_admixture(
            1, derived="pop_0", ancestral=["pop_1", "pop_2"], proportions=[0.5, 0.5]
        )
        demography.add_population_split(2, derived=[1, 2], ancestral=3)
        ts2 = msprime.sim_ancestry({0: 10}, demography=demography, random_seed=1234)
        assert ts1.equals(ts2, ignore_provenance=True)

    def test_population_split(self):
        demography = msprime.Demography.isolated_model([1000, 1000, 1000])
        demography.add_population_split(1, derived=[0, 1], ancestral=2)
        ts1 = msprime.sim_ancestry(
            {0: 1, 1: 1}, demography=demography, random_seed=1234
        )

        demography = msprime.Demography.isolated_model([1000, 1000, 1000])
        demography.add_population_split(
            1, derived=["pop_0", "pop_1"], ancestral="pop_2"
        )
        ts2 = msprime.sim_ancestry(
            {0: 1, 1: 1}, demography=demography, random_seed=1234
        )
        assert ts1.equals(ts2, ignore_provenance=True)

    def test_migration_rate_change(self):
        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_migration_rate_change(1, source=0, dest=1, rate=1)
        ts1 = msprime.sim_ancestry(
            {0: 1, 1: 1}, demography=demography, random_seed=1234
        )

        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_migration_rate_change(1, source="pop_0", dest="pop_1", rate=1)
        ts2 = msprime.sim_ancestry(
            {0: 1, 1: 1}, demography=demography, random_seed=1234
        )
        assert ts1.equals(ts2, ignore_provenance=True)

    def test_population_parameters_change(self):
        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_population_parameters_change(1, population=0, initial_size=100)
        ts1 = msprime.sim_ancestry({0: 1}, demography=demography, random_seed=1234)

        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_population_parameters_change(
            1, population="pop_0", initial_size=100
        )
        ts2 = msprime.sim_ancestry({0: 1}, demography=demography, random_seed=1234)
        assert ts1.equals(ts2, ignore_provenance=True)

    def test_simple_bottleneck(self):
        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_simple_bottleneck(1, population=0, proportion=1)
        ts1 = msprime.sim_ancestry({0: 1}, demography=demography, random_seed=1234)

        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_simple_bottleneck(1, population="pop_0", proportion=1)
        ts2 = msprime.sim_ancestry({0: 1}, demography=demography, random_seed=1234)
        assert ts1.equals(ts2, ignore_provenance=True)

    def test_instantaneous_bottleneck(self):
        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_instantaneous_bottleneck(1, population=0, strength=1)
        ts1 = msprime.sim_ancestry({0: 1}, demography=demography, random_seed=1234)

        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.add_instantaneous_bottleneck(1, population="pop_0", strength=1)
        ts2 = msprime.sim_ancestry({0: 1}, demography=demography, random_seed=1234)
        assert ts1.equals(ts2, ignore_provenance=True)


class TestPopulationSplit:
    """
    Tests that the PopulationSplit event has the correct effects.
    """

    def test_two_pop_tree(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="AB", initial_size=100)
        demography.set_symmetric_migration_rate(["A", "B"], 0.1)
        demography.add_population_split(10, derived=["A", "B"], ancestral="AB")
        assert demography["A"].sampling_time is None
        assert demography["B"].sampling_time is None
        assert demography["AB"].sampling_time == 10
        dbg = demography.debug()
        assert len(dbg.epochs) == 2
        assert dbg.epochs[1].start_time == 10
        assert np.all(dbg.epochs[1].migration_matrix == 0)
        for derived in dbg.epochs[1].populations[:2]:
            assert derived.start_size == 0
            assert derived.end_size == 0
            assert derived.growth_rate == 0
            assert not derived.active
        ancestral = dbg.epochs[1].populations[-1]
        assert ancestral.start_size == 100
        assert ancestral.end_size == 100
        assert ancestral.growth_rate == 0
        assert ancestral.active
        ts = msprime.sim_ancestry(
            {"A": 1, "B": 1}, demography=demography, random_seed=32
        )
        assert ts.tables.nodes.time[-1] > 10

    def test_two_pop_merge_initially_active(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100, initially_active=True)
        demography.add_population(name="B", initial_size=100)
        demography.add_population_split(1, derived=["B"], ancestral="A")
        debug = demography.debug()
        assert debug.epochs[0].populations[0].active
        assert debug.epochs[0].populations[1].active
        assert debug.epochs[1].populations[0].active
        assert not debug.epochs[1].populations[1].active

    def test_two_pop_tree_growth_rates(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=10000, growth_rate=0.01)
        demography.add_population(name="B", initial_size=10000, growth_rate=0.01)
        demography.add_population(name="AB", initial_size=100000, growth_rate=0.005)
        demography.set_symmetric_migration_rate(["A", "B"], 0.1)
        demography.add_population_split(2, derived=["A", "B"], ancestral="AB")
        assert demography["A"].sampling_time is None
        assert demography["B"].sampling_time is None
        assert demography["AB"].sampling_time == 2
        dbg = demography.debug()
        assert len(dbg.epochs) == 2
        assert dbg.epochs[1].start_time == 2
        assert np.all(dbg.epochs[1].migration_matrix == 0)
        for derived in dbg.epochs[1].populations[:2]:
            assert derived.start_size == 0
            assert derived.end_size == 0
            assert derived.growth_rate == 0
        ancestral = dbg.epochs[1].populations[-1]
        assert ancestral.growth_rate == 0.005
        ts = msprime.sim_ancestry(
            {"A": 1, "B": 1}, demography=demography, random_seed=32
        )
        assert ts.tables.nodes.time[-1] > 10

    def test_three_pop_tree(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="C", initial_size=100)
        demography.add_population(name="ABC", initial_size=100)
        demography.set_symmetric_migration_rate(["A", "B", "C"], 0.1)
        demography.add_population_split(10, derived=["A", "B", "C"], ancestral="ABC")
        dbg = demography.debug()
        assert len(dbg.epochs) == 2
        assert dbg.epochs[1].start_time == 10
        assert np.all(dbg.epochs[1].migration_matrix == 0)
        for derived in dbg.epochs[1].populations[:3]:
            assert derived.start_size == 0
            assert derived.end_size == 0
            assert derived.growth_rate == 0
        ancestral = dbg.epochs[1].populations[-1]
        assert ancestral.start_size == 100
        assert ancestral.end_size == 100
        assert ancestral.growth_rate == 0
        ts = msprime.sim_ancestry(
            {"A": 1, "B": 1, "C": 1}, demography=demography, random_seed=32
        )
        assert ts.tables.nodes.time[-1] > 10

    def test_three_pop_binary_tree(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="C", initial_size=100)
        demography.add_population(name="AB", initial_size=100)
        demography.add_population(name="ABC", initial_size=100)
        demography.set_symmetric_migration_rate(["A", "B", "C"], 0.1)
        demography.add_population_split(10, derived=["A", "B"], ancestral="AB")
        demography.add_symmetric_migration_rate_change(
            10, populations=["C", "AB"], rate=0.1
        )
        demography.add_population_split(20, derived=["AB", "C"], ancestral="ABC")
        dbg = demography.debug()
        assert len(dbg.epochs) == 3
        assert dbg.epochs[1].start_time == 10
        M = np.zeros((5, 5))
        M[2, 3] = 0.1
        M[3, 2] = 0.1
        assert np.array_equal(M, dbg.epochs[1].migration_matrix)
        assert dbg.epochs[2].start_time == 20
        assert np.all(dbg.epochs[2].migration_matrix == 0)
        ts = msprime.sim_ancestry(
            {"A": 1, "B": 1, "C": 1}, demography=demography, random_seed=32
        )
        assert ts.tables.nodes.time[-1] > 20

    def test_four_pop_binary_tree(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="C", initial_size=100)
        demography.add_population(name="D", initial_size=100)
        demography.add_population(name="AB", initial_size=100)
        demography.add_population(name="CD", initial_size=100)
        demography.add_population(name="ABCD", initial_size=100)
        demography.set_symmetric_migration_rate(["A", "B", "C", "D"], 0.1)

        demography.add_population_split(10, derived=["A", "B"], ancestral="AB")
        demography.add_population_split(20, derived=["C", "D"], ancestral="CD")
        demography.add_symmetric_migration_rate_change(
            20, populations=["CD", "AB"], rate=0.1
        )
        demography.add_population_split(30, derived=["AB", "CD"], ancestral="ABCD")
        dbg = demography.debug()
        assert len(dbg.epochs) == 4
        assert dbg.epochs[1].start_time == 10
        M = np.zeros((7, 7))
        M[2, 3] = 0.1
        M[3, 2] = 0.1
        assert np.array_equal(M, dbg.epochs[1].migration_matrix)
        assert dbg.epochs[2].start_time == 20
        M = np.zeros((7, 7))
        M[4, 5] = 0.1
        M[5, 4] = 0.1
        assert np.array_equal(M, dbg.epochs[2].migration_matrix)
        assert np.all(dbg.epochs[3].migration_matrix == 0)
        ts = msprime.sim_ancestry(
            {"A": 1, "B": 1, "C": 1, "D": 1}, demography=demography, random_seed=32
        )
        assert ts.tables.nodes.time[-1] > 30
        assert ts.tables.nodes.population[-1] == 6

    def test_lineages_move_into_derived(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="AB", initial_size=100)
        demography.set_symmetric_migration_rate(["A", "B"], 0.1)
        demography.add_population_split(10, derived=["A", "B"], ancestral="AB")
        demography.add_mass_migration(10, source="AB", dest="A", proportion=1)
        with pytest.raises(_msprime.LibraryError, match="inactive population"):
            msprime.sim_ancestry(
                {"A": 1, "B": 1}, demography=demography, random_seed=32
            )

    def test_sample_derived(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="AB", initial_size=100)
        demography.set_symmetric_migration_rate(["A", "B"], 0.1)
        demography.add_population_split(10, derived=["A", "B"], ancestral="AB")
        with pytest.raises(_msprime.LibraryError, match="inactive population"):
            msprime.sim_ancestry(
                samples=[
                    msprime.SampleSet(1, "A"),
                    msprime.SampleSet(1, "B"),
                    msprime.SampleSet(1, "A", time=11),
                ],
                demography=demography,
                random_seed=32,
            )

    def test_sample_ancestral(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="AB", initial_size=100)
        demography.add_population_split(10, derived=["A", "B"], ancestral="AB")
        with pytest.raises(_msprime.LibraryError, match="inactive population"):
            msprime.sim_ancestry(
                samples=[
                    msprime.SampleSet(1, "A"),
                    msprime.SampleSet(1, "B"),
                    msprime.SampleSet(1, "AB", time=5),  # Before split
                ],
                demography=demography,
                random_seed=32,
            )

    @pytest.mark.parametrize("n", range(2, 13))
    def test_random_population_tree(self, n):
        ts = msprime.sim_ancestry(n, ploidy=1, random_seed=52)
        demography = msprime.Demography.isolated_model([10] * ts.num_nodes)
        tree = ts.first()
        for node in tree.nodes(order="timeasc"):
            if tree.is_internal(node):
                demography.add_population_split(
                    time=tree.time(node), derived=tree.children(node), ancestral=node
                )

        dbg = demography.debug()
        epoch = dbg.epochs[0]
        assert epoch.start_time == 0
        assert len(epoch.events) == 0
        for u in tree.nodes():
            if tree.is_leaf(u):
                assert epoch.populations[u].active
            else:
                assert not epoch.populations[u].active

        sm = msprime.PopulationStateMachine
        for epoch in dbg.epochs[1:]:
            ancestral = epoch.events[0].ancestral
            assert tree.time(ancestral) == epoch.start_time
            assert epoch.populations[ancestral].active
            assert epoch.populations[ancestral].state == sm.ACTIVE
            # All the populations older than ancestal should be
            # inactive
            u = tree.parent(ancestral)
            while u != tskit.NULL:
                assert not epoch.populations[u].active
                assert epoch.populations[u].state == sm.INACTIVE
                u = tree.parent(u)
            # All populations below ancestral should be previously active.
            for u in tree.nodes(ancestral):
                if u != ancestral:
                    assert not epoch.populations[u].active
                    assert epoch.populations[u].state == sm.PREVIOUSLY_ACTIVE
        assert math.isinf(epoch.end_time)

    def test_record_migrations(self):
        demography = msprime.Demography.isolated_model([10] * 2)
        demography.add_population_split(0.01, derived=[0], ancestral=1)
        ts = msprime.sim_ancestry(
            {0: 10}, demography=demography, record_migrations=True, random_seed=1234
        )
        assert ts.tables.nodes.population[-1] == 1
        assert ts.num_migrations > 1
        for migration in ts.migrations():
            assert migration.source == 0
            assert migration.dest == 1
            assert migration.time == 0.01

    def test_multi_split_3_pop(self):
        #   | <-
        #   |  |
        #   2  |
        #   ^  |
        #   |  |
        #   0  1
        demography = msprime.Demography.isolated_model([10] * 3)
        demography.add_population_split(time=1, ancestral=2, derived=[0])
        demography.add_population_split(time=2, ancestral=2, derived=[1])
        dbg = demography.debug()
        assert dbg.num_epochs == 3
        ts = msprime.sim_ancestry({0: 1, 1: 1}, demography=demography, random_seed=1234)
        assert ts.tables.nodes.population[-1] == 2

    def test_multi_split_4_pop(self):
        #   | <---
        #   |  | |
        #   3  | |
        #   ^  | |
        #   |  | |
        #   0  1 2
        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_population_split(time=1, ancestral=3, derived=[0])
        demography.add_population_split(time=2, ancestral=3, derived=[1, 2])
        dbg = demography.debug()
        assert dbg.num_epochs == 3
        ts = msprime.sim_ancestry(
            {0: 1, 1: 1, 2: 1}, demography=demography, random_seed=1234
        )
        assert ts.tables.nodes.population[-1] == 3

    def test_multi_split_sample_pop(self):
        #   | <--
        #   |   |
        #   |   3
        #   |  | |
        #   |  | |
        #   0  1 2
        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_population_split(time=1, ancestral=3, derived=[1, 2])
        demography.add_population_split(time=2, ancestral=0, derived=[3])
        dbg = demography.debug()
        assert dbg.demography[0].sampling_time == 2
        assert dbg.num_epochs == 3
        assert not dbg.epochs[0].populations[0].active
        assert dbg.epochs[0].populations[1].active
        assert dbg.epochs[0].populations[2].active
        # We can't sample from 0 because it is assumed to be inactive
        with pytest.raises(_msprime.InputError, match="sample .* inactive population"):
            msprime.sim_ancestry(
                [msprime.SampleSet(1, population=0, time=0)],
                demography=demography,
                random_seed=1234,
            )
        ts = msprime.sim_ancestry(
            {0: 1, 1: 1, 2: 1}, demography=demography, random_seed=1234
        )
        assert ts.tables.nodes.population[-1] == 0

    def test_multi_split_4_pop_sequential(self):
        #   | <---
        #   | <- |
        #   |  | |
        #   3  | |
        #   ^  | |
        #   |  | |
        #   0  1 2
        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_population_split(time=1, ancestral=3, derived=[0])
        demography.add_population_split(time=2, ancestral=3, derived=[1])
        demography.add_population_split(time=3, ancestral=3, derived=[2])
        dbg = demography.debug()
        assert dbg.num_epochs == 4
        ts = msprime.sim_ancestry(
            {0: 1, 1: 1, 2: 1}, demography=demography, random_seed=1234
        )
        assert ts.tables.nodes.population[-1] == 3

    def test_sampling_time_not_overwritten_in_ancestral(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=100)
        demography.add_population(name="B", initial_size=100)
        demography.add_population(name="C", initial_size=100, sampling_time=12)
        demography.add_population_split(10, derived=["A", "B"], ancestral="C")
        assert demography["A"].sampling_time is None
        assert demography["B"].sampling_time is None
        assert demography["C"].sampling_time == 12
        validated = demography.validate()
        assert validated["A"].sampling_time == 0
        assert validated["B"].sampling_time == 0
        assert validated["C"].sampling_time == 12


class TestAdmixture:
    """
    Tests for the admixture functionality.
    """

    def test_proportions_sum_to_one(self):
        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_admixture(
            1, derived=0, ancestral=[1, 2], proportions=[1 / 3, 2 / 3]
        )
        assert math.isclose(sum(demography.events[0].proportions), 1)

        for N in range(1, 100):
            demography = msprime.Demography.isolated_model([10] * (N + 1))
            demography.add_admixture(
                1, derived=N, ancestral=range(N), proportions=[1 / N] * N
            )
            assert math.isclose(sum(demography.events[0].proportions), 1)

        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_admixture(
            1, derived=0, ancestral=[1, 2, 3], proportions=[0.1, 0.1, 0.8]
        )
        assert math.isclose(sum(demography.events[0].proportions), 1)

        demography = msprime.Demography.isolated_model([10] * 4)
        with pytest.raises(ValueError, match="Sum of the admixture proportions"):
            # Arguably this is too strict, but we can relax later if needs be.
            # However, probably more useful to allow the N - 1 form, and compute
            # the remaining value, than to have a high tolerance for error.
            demography.add_admixture(
                1, derived=0, ancestral=[1, 2, 3], proportions=[0.33, 0.33, 0.33]
            )

    def test_4_pop_diamond(self):
        #     3
        #    / \
        #   1   2
        #    \ /
        #     0
        # Equal probabilities of going either direction
        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_admixture(
            1, derived=0, ancestral=[1, 2], proportions=[1 / 2, 1 / 2]
        )
        demography.add_population_split(2, derived=[1, 2], ancestral=3)
        ts = msprime.sim_ancestry({0: 100}, demography=demography, random_seed=1234)
        # We should have coalescences in all 4 populations
        assert set(ts.tables.nodes.population) == set(range(4))
        assert ts.tables.nodes.population[-1] == 3

        # Only go via pop 1
        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_admixture(1, derived=0, ancestral=[1, 2], proportions=[1, 0])
        demography.add_population_split(2, derived=[1, 2], ancestral=3)
        ts = msprime.sim_ancestry({0: 100}, demography=demography, random_seed=1234)
        # We should have no coalescences in 2
        assert set(ts.tables.nodes.population) == {0, 1, 3}
        assert ts.tables.nodes.population[-1] == 3

        # Only go via pop 2
        demography = msprime.Demography.isolated_model([10] * 4)
        demography.add_admixture(1, derived=0, ancestral=[1, 2], proportions=[0, 1])
        demography.add_population_split(2, derived=[1, 2], ancestral=3)
        ts = msprime.sim_ancestry({0: 100}, demography=demography, random_seed=1234)
        # We should have no coalescences in 1
        assert set(ts.tables.nodes.population) == {0, 2, 3}
        assert ts.tables.nodes.population[-1] == 3

    def test_multi_admixture(self):
        #     4
        #    / \
        #   2   3
        #   | x |
        #   0   1
        # Two different admixed populations from the same sources.
        demography = msprime.Demography.isolated_model([10] * 5)
        demography.add_admixture(
            1, derived=0, ancestral=[2, 3], proportions=[1 / 2, 1 / 2]
        )
        demography.add_admixture(
            1, derived=1, ancestral=[2, 3], proportions=[2 / 3, 1 / 3]
        )
        demography.add_population_split(3, derived=[2, 3], ancestral=4)
        ts = msprime.sim_ancestry(
            {0: 100, 1: 100}, demography=demography, random_seed=1234
        )
        # We should have coalescences in all 5 populations
        assert set(ts.tables.nodes.population) == set(range(5))
        assert ts.tables.nodes.population[-1] == 4

    def test_admix_into_ancestral(self):
        #          7
        #         / \
        #   |------->
        #   |--> 5   6
        #   | x / \ / \
        #   0   1 2 3 4
        # Amixture where ancestral populations are from population splits.
        demography = msprime.Demography.isolated_model([10] * 8)
        demography.add_population_split(1, derived=[3, 4], ancestral=6)
        demography.add_population_split(1, derived=[1, 2], ancestral=5)
        demography.add_admixture(
            2, derived=0, ancestral=[5, 6], proportions=[1 / 2, 1 / 2]
        )
        demography.add_population_split(3, derived=[5, 6], ancestral=7)
        ts = msprime.sim_ancestry(
            {0: 100, 1: 100}, demography=demography, random_seed=1234
        )
        assert ts.tables.nodes.population[-1] == 7

    def test_1_pop_admix(self):
        demography = msprime.Demography.isolated_model([10] * 2)
        demography.add_admixture(1, derived=0, ancestral=[1], proportions=[1])
        ts = msprime.sim_ancestry({0: 100}, demography=demography, random_seed=1234)
        assert ts.tables.nodes.population[-1] == 1

    def test_record_migrations(self):
        demography = msprime.Demography.isolated_model([10] * 2)
        demography.add_admixture(0.01, derived=0, ancestral=[1], proportions=[1])
        ts = msprime.sim_ancestry(
            {0: 10}, demography=demography, record_migrations=True, random_seed=1234
        )
        assert ts.tables.nodes.population[-1] == 1
        assert ts.num_migrations > 1
        for migration in ts.migrations():
            assert migration.source == 0
            assert migration.dest == 1
            assert migration.time == 0.01


class TestPopulationLoops:
    def verify_loop_error(self, demography):
        msg = "derived population in a population split must be active"
        with pytest.raises(
            _msprime.InputError, match="sample a lineage from an inactive"
        ):
            msprime.sim_ancestry({0: 1}, demography=demography)
        with pytest.raises(_msprime.LibraryError, match=msg):
            demography.debug()

    def test_ABA_loop(self):
        demography = msprime.Demography()
        demography.add_population(name="A")
        demography.add_population(name="B")
        demography.add_population_split(1, derived=["A"], ancestral="B")
        demography.add_population_split(2, derived=["B"], ancestral="A")
        self.verify_loop_error(demography)

    def test_ABA_admixture_loop(self):
        demography = msprime.Demography()
        demography.add_population(name="A")
        demography.add_population(name="B")
        demography.add_admixture(0.1, derived="A", ancestral=["B"], proportions=[1])
        demography.add_admixture(0.2, derived="B", ancestral=["A"], proportions=[1])
        msg = "All ancestral populations in admixture must already be active"
        with pytest.raises(_msprime.LibraryError, match=msg):
            demography.debug()
        with pytest.raises(_msprime.LibraryError, match=msg):
            msprime.sim_ancestry(
                {"A": 10, "B": 10}, demography=demography, random_seed=123
            )

    def test_admixture_derived_inactive(self):
        demography = msprime.Demography()
        demography.add_population(name="A")
        demography.add_population(name="B")
        demography.add_population(name="C")
        demography.add_population_split(0.1, derived=["A"], ancestral="B")
        demography.add_admixture(0.2, derived="A", ancestral=["C"], proportions=[1])
        msg = "derived population in an admixture must be active"
        with pytest.raises(_msprime.LibraryError, match=msg):
            demography.debug()
        with pytest.raises(_msprime.LibraryError, match=msg):
            msprime.sim_ancestry(
                {"A": 10, "B": 10}, demography=demography, random_seed=123
            )

    def test_admixture_ancestral_inactive(self):
        demography = msprime.Demography()
        demography.add_population(name="A")
        demography.add_population(name="B")
        demography.add_population(name="C")
        demography.add_population_split(0.1, derived=["A"], ancestral="B")
        demography.add_admixture(0.2, derived="C", ancestral=["A"], proportions=[1])
        msg = "ancestral populations in admixture must already be active"
        with pytest.raises(_msprime.LibraryError, match=msg):
            demography.debug()
        with pytest.raises(_msprime.LibraryError, match=msg):
            msprime.sim_ancestry(
                {"A": 10, "C": 10}, demography=demography, random_seed=123
            )

    def test_ABCA_loop(self):
        demography = msprime.Demography()
        demography.add_population(name="A")
        demography.add_population(name="B")
        demography.add_population(name="C")
        demography.add_population_split(1, derived=["A"], ancestral="B")
        demography.add_population_split(2, derived=["B"], ancestral="C")
        demography.add_population_split(3, derived=["C"], ancestral="A")
        self.verify_loop_error(demography)

    @pytest.mark.parametrize("n", range(2, 8))
    def test_line_loop(self, n):
        demography = msprime.Demography.isolated_model([1] * n)
        for j in range(n):
            if j < n - 1:
                demography.add_population_split(j + 1, derived=[j], ancestral=j + 1)

        demography.add_population_split(n, derived=[n - 1], ancestral=0)
        self.verify_loop_error(demography)

    def test_ABCDB_loop(self):
        demography = msprime.Demography()
        demography.add_population(name="A")
        demography.add_population(name="B")
        demography.add_population(name="C")
        demography.add_population(name="D")
        demography.add_population_split(1, derived=["A"], ancestral="B")
        demography.add_population_split(2, derived=["B"], ancestral="C")
        demography.add_population_split(3, derived=["C"], ancestral="D")
        demography.add_population_split(3, derived=["D"], ancestral="B")
        msg = "Attempt to set a previously active population to active"
        with pytest.raises(_msprime.LibraryError, match=msg):
            msprime.sim_ancestry({0: 1, 2: 1}, demography=demography, random_seed=1)
        with pytest.raises(_msprime.LibraryError, match=msg):
            demography.debug()


class TestNormaliseLineageMovements:
    """
    Tests for the function to normalise events that result in lineage
    movements into equivalent MassMigrations.
    """

    def test_single_mass_migration(self):
        demography = msprime.Demography.isolated_model([10] * 10)
        mm = msprime.MassMigration(time=1, source=1, dest=0, proportion=0.2)
        assert demography._normalise_lineage_movements([mm]) == {
            1: [msprime.LineageMovement(1, 0, 0.2)]
        }

    def test_filter_other_events(self):
        demography = msprime.Demography.isolated_model([10] * 10)
        events = [
            msprime.CensusEvent(0.1),
            msprime.PopulationParametersChange(0.1, growth_rate=10, initial_size=1),
            msprime.MigrationRateChange(0.1, rate=0),
            msprime.SymmetricMigrationRateChange(0.1, [0, 1], rate=0),
            msprime.InstantaneousBottleneck(0.1, population=0, strength=100),
            msprime.SimpleBottleneck(0.1, population=1, proportion=0.1),
            # Just include the mass migration here.
            msprime.MassMigration(time=0.1, source=1, dest=0, proportion=0.2),
        ]
        # Only lineage movements are included.
        assert demography._normalise_lineage_movements(events) == {
            1: [msprime.LineageMovement(1, 0, 0.2)]
        }

    def test_convert_population_split(self):
        demography = msprime.Demography.isolated_model([10] * 10)
        events = [msprime.PopulationSplit(0.1, [0], 2)]
        equiv = {
            0: [msprime.LineageMovement(0, 2, 1.0)],
        }
        assert equiv == demography._normalise_lineage_movements(events)

        events = [msprime.PopulationSplit(0.1, [0, 1], 2)]
        equiv = {
            0: [msprime.LineageMovement(0, 2, 1.0)],
            1: [msprime.LineageMovement(1, 2, 1.0)],
        }
        assert equiv == demography._normalise_lineage_movements(events)

        events = [msprime.PopulationSplit(0.1, [1, 0], 2)]
        assert equiv == demography._normalise_lineage_movements(events)

        events = [msprime.PopulationSplit(0.1, [1, 0, 3], 2)]
        equiv = {
            0: [msprime.LineageMovement(0, 2, 1.0)],
            1: [msprime.LineageMovement(1, 2, 1.0)],
            3: [msprime.LineageMovement(3, 2, 1.0)],
        }
        assert equiv == demography._normalise_lineage_movements(events)

    def test_convert_admixture(self):
        demography = msprime.Demography.isolated_model([10] * 10)
        demography.add_admixture(time=0.1, derived=0, ancestral=[1], proportions=[1])
        equiv = {
            0: [msprime.LineageMovement(0, 1, 1.0)],
        }
        assert equiv == demography._normalise_lineage_movements(demography.events)

        demography = msprime.Demography.isolated_model([10] * 10)
        demography.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2], proportions=[0, 1]
        )
        equiv = {
            0: [msprime.LineageMovement(0, 1, 0.0), msprime.LineageMovement(0, 2, 1.0)]
        }
        assert equiv == demography._normalise_lineage_movements(demography.events)

        demography = msprime.Demography.isolated_model([10] * 10)
        demography.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[0.25, 0.25, 0.5]
        )
        equiv = {
            0: [
                msprime.LineageMovement(0, 1, 0.25),
                msprime.LineageMovement(0, 2, 1 / 3),
                msprime.LineageMovement(0, 3, 1),
            ]
        }
        assert equiv == demography._normalise_lineage_movements(demography.events)

        # We reorder the populations by ID
        demography = msprime.Demography.isolated_model([10] * 10)
        demography.add_admixture(
            time=0.1, derived=0, ancestral=[3, 2, 1], proportions=[0.5, 0.25, 0.25]
        )
        equiv = {
            0: [
                msprime.LineageMovement(0, 1, 0.25),
                msprime.LineageMovement(0, 2, 1 / 3),
                msprime.LineageMovement(0, 3, 1),
            ]
        }
        assert equiv == demography._normalise_lineage_movements(demography.events)


class TestProportionConversion:
    """
    Tests for the conversion of absolute and sequential proportions.
    """

    def verify_round_trip(self, P):
        C = msprime.demography._proportions_to_sequential(P)
        # assert math.isclose(C[-1], total)
        assert np.all(np.array(C) >= 0)
        Pp = msprime.demography._sequential_to_proportions(C)
        np.allclose(P, Pp)

    def test_simple(self):
        assert msprime.demography._proportions_to_sequential([1]) == [1]
        assert msprime.demography._sequential_to_proportions([1]) == [1]
        assert msprime.demography._proportions_to_sequential([0.25, 0.25, 0.5]) == [
            0.25,
            1 / 3,
            1,
        ]
        assert msprime.demography._sequential_to_proportions([0.25, 1 / 3, 1]) == [
            0.25,
            0.25,
            0.5,
        ]

    @pytest.mark.parametrize(
        "P",
        [
            [1],
            [1 / 2, 1 / 2],
            [1 / 3] * 3,
            [1 / 6] * 6,
            [2 / 6, 1 / 6, 2 / 6, 1 / 6],
            np.ones(100) / 100,
        ],
    )
    def test_simple_examples(self, P):
        self.verify_round_trip(P)

    @pytest.mark.parametrize("n", [1, 10, 100, 1000])
    def test_random(self, n):
        rng = np.random.default_rng(n)
        z = rng.random(n)
        self.verify_round_trip(z / np.sum(z))

    @pytest.mark.parametrize(
        "P",
        [
            [0.1],
            [1 / 2, 1 / 4],
            [1 / 3] * 2,
            [1 / 6] * 5,
            [1 / 6, 1 / 6, 2 / 6, 1 / 6],
            np.ones(100) / 100,
        ],
    )
    def test_examples_less_than_1(self, P):
        self.verify_round_trip(P)


class TestDemographyEquivalent:
    """
    Tests for the function to normalise events that result in lineage
    movements into equivalent MassMigrations.
    """

    def test_different_population_parameters(self):
        d1 = msprime.Demography.isolated_model([1, 1])
        d2 = msprime.Demography.isolated_model([1, 1])
        d1.assert_equivalent(d2)
        assert d1.is_equivalent(d2)
        d2 = msprime.Demography.isolated_model([1, 2])
        assert not d1.is_equivalent(d2)
        with pytest.raises(AssertionError, match="1.0 ≠ 2.0"):
            d1.assert_equivalent(d2)
        d2 = msprime.Demography.isolated_model([1, 1], growth_rate=[0.1, 0.1])
        assert not d1.is_equivalent(d2)
        with pytest.raises(AssertionError, match="growth_rate not equal"):
            d1.assert_equivalent(d2)

    def test_different_migration_rates(self):
        d1 = msprime.Demography.stepping_stone_model([1, 1], migration_rate=0.1)
        d2 = msprime.Demography.stepping_stone_model([1, 1], migration_rate=0.1)
        d1.assert_equivalent(d2)
        assert d1.is_equivalent(d2)
        d2 = msprime.Demography.stepping_stone_model([1, 1], migration_rate=0.2)
        assert not d1.is_equivalent(d2)
        with pytest.raises(AssertionError, match="Migration matrices"):
            d1.assert_equivalent(d2)

        d1 = msprime.Demography.stepping_stone_model([1, 1, 1], migration_rate=0.1)
        d2 = msprime.Demography.stepping_stone_model([1, 1, 1], migration_rate=0.1)
        d2.set_migration_rate(source=0, dest=1, rate=0.5)
        with pytest.raises(AssertionError, match="pop_0, pop_1"):
            d1.assert_equivalent(d2)

    def test_different_names(self):
        d1 = msprime.Demography()
        d1.add_population(name="x")
        d2 = msprime.Demography()
        d2.add_population(name="y")
        with pytest.raises(AssertionError, match="names differ"):
            d1.assert_equivalent(d2)

    def test_population_split(self):
        d1 = msprime.Demography.isolated_model([1, 1, 1])
        d1.add_population_split(1.0, derived=[0, 1], ancestral=2)
        assert d1.is_equivalent(d1)

        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_split(1.0, derived=[0, 1], ancestral=2)
        assert d1.is_equivalent(d2)

        # The order of derived populations in a split isn't significant.
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_split(1.0, derived=[1, 0], ancestral=2)
        d1.assert_equivalent(d2)

        # The equivilant mass migrations aren't the same as the population
        # states don't match.
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_mass_migration(1.0, source=0, dest=2, proportion=1.0)
        d2.add_mass_migration(1.0, source=1, dest=2, proportion=1.0)
        assert not d1.is_equivalent(d2)

    def test_different_num_epochs(self):
        d1 = msprime.Demography.isolated_model([1, 1, 1])
        d1.add_migration_rate_change(1.0, rate=1)
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        assert not d1.is_equivalent(d2)
        with pytest.raises(AssertionError, match="Number of epochs not equal: 2 ≠ 1"):
            d1.assert_equivalent(d2)

    def test_different_num_populations(self):
        d1 = msprime.Demography.isolated_model([1, 1, 1])
        d2 = msprime.Demography.isolated_model([1, 1])
        assert not d1.is_equivalent(d2)
        with pytest.raises(
            AssertionError, match="Number of populations not equal: 3 ≠ 2"
        ):
            d1.assert_equivalent(d2)

    def test_different_sampling_times(self):
        d1 = msprime.Demography.isolated_model([1, 1])
        d2 = msprime.Demography.isolated_model([1, 1])
        d1[0].sampling_time = 0.01
        assert not d1.is_equivalent(d2)
        with pytest.raises(
            AssertionError, match="Sampling times not equal for pop_0: 0.01 ≠ 0"
        ):
            d1.assert_equivalent(d2)

    def test_different_epoch_times(self):
        d1 = msprime.Demography.isolated_model([1, 1, 1])
        d1.add_population_parameters_change(1.0, initial_size=2)
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_parameters_change(2.0, initial_size=2)
        assert not d1.is_equivalent(d2)
        with pytest.raises(AssertionError, match="at different times: 1.0 ≠ 2.0"):
            d1.assert_equivalent(d2)

        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_parameters_change(1.00001, initial_size=2)
        assert not d1.is_equivalent(d2)

        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_parameters_change(1 + 1e-10, initial_size=2)
        d1.assert_equivalent(d2)
        with pytest.raises(AssertionError, match="at different times"):
            d1.assert_equivalent(d2, rel_tol=1e-12)

    def test_different_population_size_epoch(self):
        d1 = msprime.Demography.isolated_model([1, 1, 1])
        d1.add_population_parameters_change(time=1.0, initial_size=2)
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_parameters_change(time=1.0, initial_size=3)
        with pytest.raises(AssertionError, match="start_size not equal"):
            d1.assert_equivalent(d2)
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_parameters_change(time=1.0, initial_size=2 + 1e-10)
        d1.assert_equivalent(d2)

    def test_different_growth_rate_epoch(self):
        d1 = msprime.Demography.isolated_model([1, 1, 1])
        d1.add_population_parameters_change(time=1.0, growth_rate=2)
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_parameters_change(time=1.0, growth_rate=3)
        with pytest.raises(AssertionError, match="growth_rate not equal"):
            d1.assert_equivalent(d2)
        d2 = msprime.Demography.isolated_model([1, 1, 1])
        d2.add_population_parameters_change(time=1.0, growth_rate=2 + 1e-10)
        d1.assert_equivalent(d2)

    def test_state_changes_unsupported(self):
        d1 = msprime.Demography.isolated_model([1, 1])
        d1.add_population_parameters_change(0.1, initial_size=1)
        d2 = msprime.Demography.isolated_model([1, 1])
        d2.add_simple_bottleneck(0.1, population=0, proportion=1)
        with pytest.raises(ValueError, match="State change events"):
            d1.assert_equivalent(d2)

        d1 = msprime.Demography.isolated_model([1, 1])
        d1.add_instantaneous_bottleneck(0.1, population=0, strength=1)
        d2 = msprime.Demography.isolated_model([1, 1])
        d2.add_population_parameters_change(0.1, initial_size=1)
        with pytest.raises(ValueError, match="State change events"):
            d1.assert_equivalent(d2)

    def test_admixture(self):
        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[0.1, 0.2, 0.7]
        )
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[0.1, 0.2, 0.7]
        )
        d1.assert_equivalent(d2)

        # Robust to permutations
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_admixture(
            time=0.1, derived=0, ancestral=[3, 2, 1], proportions=[0.7, 0.2, 0.1]
        )
        d1.assert_equivalent(d2)

        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_admixture(
            time=0.1, derived=0, ancestral=[2, 3, 1], proportions=[0.2, 0.7, 0.1]
        )
        d1.assert_equivalent(d2)

        # We detect significantly different model
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[0.1, 0.21, 0.69]
        )
        assert not d1.is_equivalent(d2)

    def test_admixture_equivalent_mass_migrations(self):
        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[1 / 3, 1 / 6, 1 / 2]
        )
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_mass_migration(time=0.1, source=0, dest=1, proportion=1 / 3)
        d2.add_mass_migration(time=0.1, source=0, dest=2, proportion=1 / 4)
        d2.add_mass_migration(time=0.1, source=0, dest=3, proportion=1)
        with pytest.raises(AssertionError, match="State mismatch"):
            d1.assert_equivalent(d2)

    def test_different_lineage_movement_dest(self):
        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[1 / 3, 1 / 6, 1 / 2]
        )
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 4], proportions=[1 / 3, 1 / 6, 1 / 2]
        )
        with pytest.raises(AssertionError, match="movement destination"):
            d1.assert_equivalent(d2)

        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_mass_migration(time=0.1, source=0, dest=1, proportion=1)
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_mass_migration(time=0.1, source=0, dest=2, proportion=1)
        with pytest.raises(AssertionError, match="movement destination"):
            d1.assert_equivalent(d2)

    def test_different_lineage_movement_source_pops(self):
        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_mass_migration(time=0.1, source=1, dest=0, proportion=1)
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_mass_migration(time=0.1, source=2, dest=0, proportion=1)
        with pytest.raises(AssertionError, match="set of populations"):
            d1.assert_equivalent(d2)

    def test_different_lineage_movement_numbers_out(self):
        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_mass_migration(time=0.1, source=1, dest=0, proportion=0.5)
        d1.add_mass_migration(time=0.1, source=1, dest=0, proportion=0.5)
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_mass_migration(time=0.1, source=1, dest=0, proportion=0.5)
        with pytest.raises(
            AssertionError, match="number of normalised lineage movements"
        ):
            d1.assert_equivalent(d2)

    def test_different_lineage_movement_proportion(self):
        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_mass_migration(time=0.1, source=1, dest=0, proportion=0.5)
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_mass_migration(time=0.1, source=1, dest=0, proportion=0.51)
        with pytest.raises(AssertionError, match="movement proportions"):
            d1.assert_equivalent(d2)

    def test_different_admixture_proportions(self):
        d1 = msprime.Demography.isolated_model([10] * 10)
        d1.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[1 / 3, 1 / 6, 1 / 2]
        )
        d2 = msprime.Demography.isolated_model([10] * 10)
        d2.add_admixture(
            time=0.1, derived=0, ancestral=[1, 2, 3], proportions=[1 / 2, 1 / 6, 1 / 3]
        )
        with pytest.raises(AssertionError, match="movement proportions"):
            d1.assert_equivalent(d2)


class TestFromOldStyleMap:
    """
    Tests for the conversion of populations into a new-style demography
    using a map of population IDs, describing active populations in
    each epoch and how they map to original populations.
    """

    def test_single_epoch(self):
        d = msprime.Demography.from_old_style(
            [msprime.PopulationConfiguration(initial_size=1)], population_map=[{"A": 0}]
        )
        assert d.num_populations == 1
        assert d.debug().num_epochs == 1
        assert "A" in d

    def test_one_pop_size_changes(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[
                msprime.PopulationParametersChange(1, initial_size=1),
                msprime.PopulationParametersChange(2, initial_size=2),
                msprime.PopulationParametersChange(3, initial_size=3),
            ],
            population_map=[{"A": 0}, {"A": 0}, {"A": 0}, {"A": 0}],
        )
        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population_parameters_change(1, population="A", initial_size=1)
        d2.add_population_parameters_change(2, population="A", initial_size=2)
        d2.add_population_parameters_change(3, population="A", initial_size=3)
        d1.assert_equivalent(d2)

    def test_initial_migration_rates(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            migration_matrix=[
                [0, 0.1],
                [0.2, 0],
            ],
            population_map=[{"A": 0, "B": 1}],
        )
        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population(initial_size=1, name="B")
        d2.set_migration_rate(source="A", dest="B", rate=0.1)
        d2.set_migration_rate(source="B", dest="A", rate=0.2)
        d1.assert_equivalent(d2)

    def test_migration_rate_change(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[
                msprime.MigrationRateChange(1, source=0, dest=1, rate=0.1),
                msprime.MigrationRateChange(1, source=1, dest=0, rate=0.2),
            ],
            population_map=[{"A": 0, "B": 1}, {"A": 0, "B": 1}],
        )
        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population(initial_size=1, name="B")
        d2.add_migration_rate_change(time=1, source="A", dest="B", rate=0.1)
        d2.add_migration_rate_change(time=1, source="B", dest="A", rate=0.2)
        d1.assert_equivalent(d2)

    def test_two_pop_tree(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[msprime.MassMigration(1, source=1, dest=0)],
            population_map=[{"A": 0, "B": 1}, {"C": 0}],
        )
        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population(initial_size=1, name="B")
        d2.add_population(initial_size=1, name="C")
        d2.add_population_split(1, derived=["A", "B"], ancestral="C")
        d1.assert_equivalent(d2)

        # Also same if we map C to 1
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[msprime.MassMigration(1, source=0, dest=1)],
            population_map=[{"A": 0, "B": 1}, {"C": 1}],
        )
        d1.assert_equivalent(d2)

    def test_two_pop_merge_into_first(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[msprime.MassMigration(1, source=1, dest=0)],
            population_map=[{"A": 0, "B": 1}, {"A": 0}],
        )
        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A", initially_active=True)
        d2.add_population(initial_size=1, name="B")
        d2.add_population_split(1, derived=["B"], ancestral="A")
        d1.assert_equivalent(d2)

        # Try the other way around also
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[msprime.MassMigration(1, source=0, dest=1)],
            population_map=[{"A": 0, "B": 1}, {"B": 1}],
        )
        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population(initial_size=1, name="B", initially_active=True)
        d2.add_population_split(1, derived=["A"], ancestral="B")
        d1.assert_equivalent(d2)

    def test_two_pop_tree_no_mass_migration(self):
        with pytest.raises(ValueError, match="Insufficient MassMigrations"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                demographic_events=[
                    msprime.MigrationRateChange(1, source=1, dest=0, rate=0)
                ],
                population_map=[{"A": 0, "B": 1}, {"C": 0}],
            )

    def test_two_pop_tree_incorrect_mass_migration(self):
        with pytest.raises(ValueError, match="MassMigration associated with pop"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                demographic_events=[
                    msprime.MassMigration(1, source=0, dest=1, proportion=0.1)
                ],
                population_map=[{"A": 0, "B": 1}, {"C": 0}],
            )

    def test_three_pop_multiple_splits(self):
        # 1 and 2 sequentially merge into 0
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[
                msprime.MassMigration(1, source=2, dest=0),
                msprime.MassMigration(2, source=1, dest=0),
            ],
            population_map=[{"A": 0, "B": 1, "C": 2}, {"A": 0, "B": 1}, {"A": 0}],
        )
        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A", initially_active=True)
        d2.add_population(initial_size=1, name="B")
        d2.add_population(initial_size=1, name="C")
        d2.add_population_split(1, derived=["C"], ancestral="A")
        d2.add_population_split(2, derived=["B"], ancestral="A")
        d1.assert_equivalent(d2)

    def test_two_pop_admix(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[
                msprime.MassMigration(1, source=2, dest=0, proportion=0.5),
                msprime.MassMigration(1, source=2, dest=1, proportion=1),
            ],
            population_map=[{"A": 0, "B": 1, "C": 2}, {"A": 0, "B": 1}],
        )

        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population(initial_size=1, name="B")
        d2.add_population(initial_size=1, name="C")
        d2.add_admixture(1, derived="C", ancestral=["A", "B"], proportions=[0.5, 0.5])
        d1.assert_equivalent(d2)

    def test_two_pop_admix_not_summing_to_one(self):
        with pytest.raises(ValueError, match="don't sum to 1"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                demographic_events=[
                    msprime.MassMigration(1, source=2, dest=0, proportion=0.5),
                    msprime.MassMigration(1, source=2, dest=1, proportion=0.5),
                ],
                population_map=[{"A": 0, "B": 1, "C": 2}, {"A": 0, "B": 1}],
            )

    def test_two_pop_admix_plus_simultaneous_pulse(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[
                msprime.MassMigration(1, source=2, dest=0, proportion=0.5),
                msprime.MassMigration(1, source=2, dest=1, proportion=1),
                # Pulse of migration from 3 into 0
                msprime.MassMigration(1, source=3, dest=0, proportion=0.5),
            ],
            population_map=[{"A": 0, "B": 1, "C": 2, "D": 3}, {"A": 0, "B": 1, "D": 3}],
        )

        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population(initial_size=1, name="B")
        d2.add_population(initial_size=1, name="C")
        d2.add_population(initial_size=1, name="D")
        d2.add_admixture(1, derived="C", ancestral=["A", "B"], proportions=[0.5, 0.5])
        d2.add_mass_migration(1, source="D", dest="A", proportion=0.5)
        d1.assert_equivalent(d2)

    def test_two_pop_admix_plus_simultaneous_pulse_with_split(self):
        d1 = msprime.Demography.from_old_style(
            [
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
                msprime.PopulationConfiguration(initial_size=1),
            ],
            demographic_events=[
                msprime.MassMigration(1, source=2, dest=0, proportion=0.5),
                msprime.MassMigration(1, source=2, dest=1, proportion=1),
                # Pulse of migration from 3 into 0
                msprime.MassMigration(1, source=3, dest=0, proportion=0.5),
                # A, B, D merge into E
                msprime.MassMigration(2, source=1, dest=0, proportion=1),
                msprime.MassMigration(2, source=3, dest=0, proportion=1),
            ],
            population_map=[
                {"A": 0, "B": 1, "C": 2, "D": 3},
                {"A": 0, "B": 1, "D": 3},
                {"E": 0},
            ],
        )

        d2 = msprime.Demography()
        d2.add_population(initial_size=1, name="A")
        d2.add_population(initial_size=1, name="B")
        d2.add_population(initial_size=1, name="C")
        d2.add_population(initial_size=1, name="D")
        d2.add_population(initial_size=1, name="E")
        d2.add_admixture(1, derived="C", ancestral=["A", "B"], proportions=[0.5, 0.5])
        d2.add_mass_migration(1, source="D", dest="A", proportion=0.5)
        d2.add_population_split(2, derived=["A", "B", "D"], ancestral="E")
        d1.assert_equivalent(d2)

    def test_epoch_mismatch(self):
        with pytest.raises(ValueError, match="number of epochs"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                demographic_events=[],
                population_map=[{"A": 0, "B": 1}, {"C": 0}],
            )

        with pytest.raises(ValueError, match="number of epochs"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                demographic_events=[
                    msprime.MassMigration(1, source=1, dest=0, proportion=1),
                    msprime.PopulationParametersChange(1.1, initial_size=2),
                ],
                population_map=[{"A": 0, "B": 1}, {"C": 0}],
            )

    def test_bad_population_ref(self):
        with pytest.raises(ValueError, match="Bad population reference"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                population_map=[{"A": 0, "B": 3}],
            )

        with pytest.raises(ValueError, match="Bad population reference"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                population_map=[{"A": -1}],
            )

    def test_repeated_ids(self):
        with pytest.raises(ValueError, match="must be unique"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                population_map=[{"A": 0, "B": 0}],
            )

    def test_mismatch_populations(self):
        with pytest.raises(ValueError, match="entries for all"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                population_map=[{"A": 0}],
            )

    def test_unsupported_events(self):
        unsupported_events = [
            msprime.SimpleBottleneck(time=1, population=0, proportion=1),
            msprime.InstantaneousBottleneck(time=1, population=0, strength=1),
            msprime.CensusEvent(time=1),
        ]
        for event in unsupported_events:
            with pytest.raises(ValueError, match="Only MassMigration"):
                msprime.Demography.from_old_style(
                    [
                        msprime.PopulationConfiguration(initial_size=1),
                    ],
                    demographic_events=[event],
                    population_map=[{"A": 0}, {"A": 0}],
                )

    def test_migration_into_inactive_error(self):
        with pytest.raises(ValueError, match="Non zero migration from an active"):
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                migration_matrix=[[0, 1], [0, 0]],
                demographic_events=[
                    msprime.MassMigration(1, source=1, dest=0),
                ],
                population_map=[{"A": 0, "B": 1}, {"A": 0}],
            )

    def test_migration_out_of_inactive_warning(self):
        with warnings.catch_warnings(record=True) as w:
            msprime.Demography.from_old_style(
                [
                    msprime.PopulationConfiguration(initial_size=1),
                    msprime.PopulationConfiguration(initial_size=1),
                ],
                migration_matrix=[[0, 0], [1, 0]],
                demographic_events=[
                    msprime.MassMigration(1, source=1, dest=0),
                ],
                population_map=[{"A": 0, "B": 1}, {"A": 0}],
            )
        assert len(w) == 1
        assert str(w[0].message).startswith("Migration out of inactive")


class TestLineageMovementEvents:
    def test_as_lineage_movement_abstract(self):
        e = msprime.LineageMovementEvent(1234)
        with pytest.raises(NotImplementedError):
            e._as_lineage_movements()


class TestStdpopsimModels:
    def stdpopsim_browning_admixture_model(self):
        """
        The currently released version of the AmericanAdmixture_4B11 model in
        stdpopsim has arbitrary time deltas put in between different events,
        making it incomparable with the models we define here. This has been
        fixed upstream, but while we're waiting for it to be released, here
        is the model copied verbatim.
        """

        # Model code was ported from Supplementary File 1.
        N0 = 7310  # initial population size
        Thum = 5920  # time (gens) of advent of modern humans
        Naf = 14474  # size of african population
        Tooa = 2040  # number of generations back to Out of Africa
        Nb = 1861  # size of out of Africa population
        mafb = 1.5e-4  # migration rate Africa and Out-of-Africa
        Teu = 920  # number generations back to Asia-Europe split
        Neu = 1032  # bottleneck population sizes
        Nas = 554
        mafeu = 2.5e-5  # mig. rates
        mafas = 7.8e-6
        meuas = 3.11e-5
        reu = 0.0038  # growth rate per generation in Europe
        ras = 0.0048  # growth rate per generation in Asia
        Tadmix = 12  # time of admixture
        Nadmix = 30000  # initial size of admixed population
        radmix = 0.05  # growth rate of admixed population
        # pop0 is Africa, pop1 is Europe, pop2 is Asia,  pop3 is admixed

        population_configurations = [
            msprime.PopulationConfiguration(
                initial_size=Naf,
                growth_rate=0.0,
            ),
            msprime.PopulationConfiguration(
                initial_size=Neu * math.exp(reu * Teu),
                growth_rate=reu,
            ),
            msprime.PopulationConfiguration(
                initial_size=Nas * math.exp(ras * Teu),
                growth_rate=ras,
            ),
            msprime.PopulationConfiguration(
                initial_size=Nadmix * math.exp(radmix * Tadmix),
                growth_rate=radmix,
            ),
        ]

        migration_matrix = [
            [0, mafeu, mafas, 0],
            [mafeu, 0, meuas, 0],
            [mafas, meuas, 0, 0],
            [0, 0, 0, 0],
        ]
        # Admixture event, 1/6 Africa, 2/6 Europe, 3/6 Asia
        admixture_event = [
            msprime.MassMigration(
                time=Tadmix, source=3, destination=0, proportion=1.0 / 6.0
            ),
            msprime.MassMigration(
                time=Tadmix, source=3, destination=1, proportion=2.0 / 5.0
            ),
            msprime.MassMigration(time=Tadmix, source=3, destination=2, proportion=1.0),
        ]
        # Asia and Europe split
        eu_event = [
            msprime.MigrationRateChange(time=Teu, rate=0.0),
            msprime.MassMigration(time=Teu, source=2, destination=1, proportion=1.0),
            msprime.PopulationParametersChange(
                time=Teu, initial_size=Nb, growth_rate=0.0, population_id=1
            ),
            msprime.MigrationRateChange(time=Teu, rate=mafb, matrix_index=(0, 1)),
            msprime.MigrationRateChange(time=Teu, rate=mafb, matrix_index=(1, 0)),
        ]
        # Out of Africa event
        ooa_event = [
            msprime.MigrationRateChange(time=Tooa, rate=0.0),
            msprime.MassMigration(time=Tooa, source=1, destination=0, proportion=1.0),
        ]
        # initial population size
        init_event = [
            msprime.PopulationParametersChange(
                time=Thum, initial_size=N0, population_id=0
            )
        ]
        demographic_events = admixture_event + eu_event + ooa_event + init_event

        return (
            population_configurations,
            migration_matrix,
            demographic_events,
        )

    def test_browning_admixture(self):
        # Compare to stdpopsim model
        (
            pop_configs,
            migration_matrix,
            events,
        ) = self.stdpopsim_browning_admixture_model()
        population_map = [
            {"AFR": 0, "EUR": 1, "EAS": 2, "ADMIX": 3},
            {"AFR": 0, "EUR": 1, "EAS": 2},
            {"AFR": 0, "OOA": 1},
            {"AMH": 0},
            {"ANC": 0},
        ]
        demog_sps = msprime.Demography.from_old_style(
            pop_configs,
            migration_matrix=migration_matrix,
            demographic_events=events,
            population_map=population_map,
        )
        demog_local = msprime.Demography._american_admixture_model()
        demog_local.assert_equivalent(demog_sps, rel_tol=1e-5)

    def test_ooa_remap(self):
        # This test is temporary while we are updating stdpopsim to use the
        # msprime APIs. See the nodes in the _ooa_model() code.
        demog_local = msprime.Demography._ooa_model()
        model_sps = stdpopsim.get_species("HomSap").get_demographic_model(
            "OutOfAfrica_3G09"
        )

        # Map from local population names into the equivalent in the stdpopsim
        # model, per epoch.
        epoch_pop_map = [
            {"YRI": 0, "CEU": 1, "CHB": 2},
            {"YRI": 0, "OOA": 1},
            {"AMH": 0},
            {"ANC": 0},
        ]

        remapped_demog = msprime.Demography.from_old_style(
            model_sps.population_configurations,
            migration_matrix=model_sps.migration_matrix,
            demographic_events=model_sps.demographic_events,
            population_map=epoch_pop_map,
        )
        demog_local.assert_equivalent(remapped_demog)

    def test_ooa_archaic(self):
        demog_local = msprime.Demography._ooa_archaic_model()
        model_sps = stdpopsim.get_species("HomSap").get_demographic_model(
            "OutOfAfricaArchaicAdmixture_5R19"
        )
        demog_sps = msprime.Demography.from_old_style(
            model_sps.population_configurations,
            migration_matrix=model_sps.migration_matrix,
            demographic_events=model_sps.demographic_events,
            population_map=[
                # Initial populations
                {"AFR": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
                # Archaic migrations start
                {"AFR": 0, "CEU": 1, "CHB": 2, "Neanderthal": 3, "ArchaicAFR": 4},
                # CEU, CHB -> OOA
                {"AFR": 0, "OOA": 1, "Neanderthal": 3, "ArchaicAFR": 4},
                # OOA -> AFR
                {"AFR": 0, "Neanderthal": 3, "ArchaicAFR": 4},
                # Archaic migrations stop
                {"AFR": 0, "Neanderthal": 3, "ArchaicAFR": 4},
                # Ancestral popsize -> 3600
                {"AFR": 0, "Neanderthal": 3, "ArchaicAFR": 4},
                # ArchaicAFR -> AFR
                {"AFR": 0, "Neanderthal": 3},
                # Neanderthal -> AFR
                {"AFR": 0},
            ],
        )
        demog_local.assert_equivalent(demog_sps)

    def test_ooa_manual(self):
        demog_local = msprime.Demography._ooa_model()
        debug_local = demog_local.debug()
        model_sps = stdpopsim.get_species("HomSap").get_demographic_model(
            "OutOfAfrica_3G09"
        )
        demog_sps = msprime.Demography.from_old_style(
            model_sps.population_configurations,
            migration_matrix=model_sps.migration_matrix,
            demographic_events=model_sps.demographic_events,
        )
        debug_sps = demog_sps.debug()

        # We map the populations the other way around in this test
        epoch_pop_map = [
            {0: "YRI", 1: "CEU", 2: "CHB"},
            {0: "YRI", 1: "OOA"},
            {0: "AMH"},
            {0: "ANC"},
        ]

        assert len(epoch_pop_map) == debug_sps.num_epochs
        assert debug_local.num_epochs == debug_sps.num_epochs
        # In epoch 0 the top corner of the migration matrix should be
        # the same.
        assert np.array_equal(
            debug_local.epochs[0].migration_matrix[:3, :3],
            debug_sps.epochs[0].migration_matrix,
        )
        assert np.all(debug_local.epochs[0].migration_matrix[3:, 3:] == 0)
        # There's only migration between OOA and YRI in epoch 1
        M = debug_local.epochs[1].migration_matrix.copy()
        assert M[0, 3] == debug_sps.epochs[1].migration_matrix[0, 1]
        assert M[3, 0] == debug_sps.epochs[1].migration_matrix[1, 0]
        M[0, 3] = 0
        M[3, 0] = 0
        assert np.all(M == 0)

        for pop_map, epoch_local, epoch_sps in zip(
            epoch_pop_map, debug_local.epochs, debug_sps.epochs
        ):
            assert epoch_local.start_time == epoch_sps.start_time
            assert epoch_local.end_time == epoch_sps.end_time
            for pop_id_sps, local_pop_name in pop_map.items():
                pop_id_local = demog_local[local_pop_name].id
                pop_local = epoch_local.populations[pop_id_local]
                pop_sps = epoch_sps.populations[pop_id_sps]
                assert pop_local.active
                assert pop_local.start_size == pop_sps.start_size
                assert pop_local.end_size == pop_sps.end_size
                assert pop_local.growth_rate == pop_sps.growth_rate
            if len(pop_map) == 1:
                assert np.all(epoch_local.migration_matrix == 0)
                assert np.all(epoch_sps.migration_matrix == 0)
