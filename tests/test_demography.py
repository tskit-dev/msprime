#
# Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

import math
import tempfile
import unittest

import msprime


class TestGrowthRates(unittest.TestCase):
    """
    Tests to see the growth rates we calculate give us the
    right values when we go through to the low-level debugging
    interface.
    """
    def test_single_growth_rate(self):
        # Set out our values in units of generations and absolute sizes.
        Ne = 1000
        growth_rate = -0.01
        end_time = 20
        end_size = Ne * math.exp(-growth_rate * end_time)
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=2, initial_size=Ne, growth_rate=growth_rate)]
        demographic_events = [
            msprime.PopulationParametersChange(time=end_time, growth_rate=0)]
        simulator = msprime.simulator_factory(
            Ne=Ne,
            population_configurations=population_configurations,
            demographic_events=demographic_events)
        ll_sim = simulator.create_ll_instance()
        ll_end_time = ll_sim.debug_demography()
        self.assertEqual(end_time, ll_end_time * 4 * Ne)
        populations = [
            msprime.Population(Ne=Ne, **d)
            for d in ll_sim.get_population_configuration()]
        self.assertEqual(len(populations), 1)
        pop = populations[0]
        self.assertEqual(pop.growth_rate, growth_rate)
        self.assertEqual(pop.initial_size, Ne)
        self.assertEqual(pop.get_size(end_time), end_size)
        # Now fast forward to the next time slice.
        ll_end_time = ll_sim.debug_demography()
        self.assertTrue(math.isinf(ll_end_time))
        populations = [
            msprime.Population(Ne=Ne, **d)
            for d in ll_sim.get_population_configuration()]
        pop = populations[0]
        self.assertEqual(pop.growth_rate, 0)
        self.assertEqual(pop.initial_size, end_size)
        self.assertEqual(pop.get_size(10), end_size)

    def test_symmetric_growth_rates(self):
        # Test a symmetric model where we start with a negative growth
        # rate and then increase back to the same value.
        Ne = 10001
        growth_rate = 0.0125
        delta_t = 50
        end_size = Ne * math.exp(-growth_rate * delta_t)
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=2, initial_size=Ne, growth_rate=growth_rate)]
        demographic_events = [
            msprime.PopulationParametersChange(
                time=delta_t, growth_rate=-growth_rate),
            msprime.PopulationParametersChange(
                time=2 * delta_t, growth_rate=0)]
        simulator = msprime.simulator_factory(
            Ne=Ne,
            population_configurations=population_configurations,
            demographic_events=demographic_events)
        ll_sim = simulator.create_ll_instance()
        ll_end_time = ll_sim.debug_demography()
        t = delta_t
        self.assertEqual(t, ll_end_time * 4 * Ne)
        populations = [
            msprime.Population(Ne=Ne, **d)
            for d in ll_sim.get_population_configuration()]
        pop = populations[0]
        self.assertEqual(pop.growth_rate, growth_rate)
        self.assertEqual(pop.initial_size, Ne)
        self.assertEqual(pop.get_size(delta_t), end_size)
        # Now fast forward to the next time slice.
        t += delta_t
        ll_end_time = ll_sim.debug_demography()
        self.assertEqual(t, ll_end_time * 4 * Ne)
        pop = [
            msprime.Population(Ne=Ne, **d)
            for d in ll_sim.get_population_configuration()][0]
        self.assertEqual(pop.growth_rate, -growth_rate)
        self.assertEqual(pop.initial_size, end_size)
        self.assertEqual(pop.get_size(delta_t), Ne)
        # Now fast forward to the next time slice.
        ll_end_time = ll_sim.debug_demography()
        self.assertTrue(math.isinf(ll_end_time))
        populations = [
            msprime.Population(Ne=Ne, **d)
            for d in ll_sim.get_population_configuration()]
        pop = populations[0]
        self.assertEqual(pop.growth_rate, 0)
        self.assertEqual(pop.initial_size, Ne)

    def test_single_growth_rate_size_change(self):
        # Set out our values in units of generations and absolute sizes.
        Ne = 1000
        growth_rate = -0.01
        end_time = 20
        end_size = Ne * math.exp(-growth_rate * end_time)
        new_size = 4 * Ne
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=2, initial_size=Ne, growth_rate=growth_rate)]
        demographic_events = [
            msprime.PopulationParametersChange(
                time=end_time, initial_size=new_size, growth_rate=0)]
        simulator = msprime.simulator_factory(
            Ne=Ne,
            population_configurations=population_configurations,
            demographic_events=demographic_events)
        ll_sim = simulator.create_ll_instance()
        ll_end_time = ll_sim.debug_demography()
        self.assertEqual(end_time, ll_end_time * 4 * Ne)
        populations = [
            msprime.Population(Ne=Ne, **d)
            for d in ll_sim.get_population_configuration()]
        self.assertEqual(len(populations), 1)
        pop = populations[0]
        self.assertEqual(pop.growth_rate, growth_rate)
        self.assertEqual(pop.initial_size, Ne)
        self.assertEqual(pop.get_size(end_time), end_size)
        # Now fast forward to the next time slice.
        ll_end_time = ll_sim.debug_demography()
        self.assertTrue(math.isinf(ll_end_time))
        populations = [
            msprime.Population(Ne=Ne, **d)
            for d in ll_sim.get_population_configuration()]
        pop = populations[0]
        self.assertEqual(pop.growth_rate, 0)
        self.assertEqual(pop.initial_size, new_size)
        self.assertEqual(pop.get_size(10), new_size)


class TestRateConversions(unittest.TestCase):
    """
    Tests for the demographic events interface.
    """
    def test_size_change(self):
        g = 100
        Ne = 1024
        new_size = 512
        event = msprime.PopulationParametersChange(
            time=g, initial_size=new_size)
        ll_event = {
            "type": "population_parameters_change",
            "time": g / (4 * Ne),
            "population_id": -1,
            "initial_size": new_size / Ne
        }
        self.assertEqual(
            event.get_ll_representation(1, Ne), ll_event)

    def test_growth_rate_change(self):
        g = 512
        Ne = 4096
        growth_rate = 1
        event = msprime.PopulationParametersChange(
            time=g, growth_rate=growth_rate, population_id=1)
        ll_event = {
            "type": "population_parameters_change",
            "time": g / (4 * Ne),
            "population_id": 1,
            "growth_rate": growth_rate * (4 * Ne)
        }
        self.assertEqual(
            event.get_ll_representation(1, Ne), ll_event)

    def test_growth_rate_and_size_change(self):
        g = 1024
        Ne = 4096
        growth_rate = 2
        initial_size = 8192
        event = msprime.PopulationParametersChange(
            time=g, initial_size=initial_size,
            growth_rate=growth_rate, population_id=1)
        ll_event = {
            "type": "population_parameters_change",
            "time": g / (4 * Ne),
            "population_id": 1,
            "initial_size": initial_size / Ne,
            "growth_rate": growth_rate * (4 * Ne)
        }
        self.assertEqual(
            event.get_ll_representation(1, Ne), ll_event)


class TestDemographyPrinter(unittest.TestCase):
    """
    Tests for the demography printer interface.
    """

    def verify_debug(
            self, population_configurations, migration_matrix,
            demographic_events):
        with tempfile.TemporaryFile("w+") as f:
            dp = msprime.DemographyPrinter(
                population_configurations, migration_matrix,
                demographic_events=demographic_events, file=f)
            dp.debug_history()
            f.seek(0)
            debug_output = f.read()
        # TODO when there is better output, write some tests to
        # verify its format.
        self.assertGreater(len(debug_output), 0)

    def test_one_population(self):
        population_configurations = [
            msprime.PopulationConfiguration(10)]
        migration_matrix = [[0]]
        demographic_events = [
            msprime.PopulationParametersChange(0.1, initial_size=2),
            msprime.PopulationParametersChange(0.1, growth_rate=10)]
        self.verify_debug(
            population_configurations, migration_matrix,
            demographic_events)

    def test_no_events(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10)]
        migration_matrix = [[0, 0], [0, 0]]
        self.verify_debug(
            population_configurations, migration_matrix, [])

    def test_demographic_events(self):
        population_configurations = [
            msprime.PopulationConfiguration(10),
            msprime.PopulationConfiguration(10)]
        migration_matrix = [[0, 0], [0, 0]]
        demographic_events = [
            msprime.PopulationParametersChange(0.1, initial_size=2),
            msprime.PopulationParametersChange(0.1, growth_rate=10),
            msprime.MassMigration(0.2, source=1, destination=0),
            msprime.MigrationRateChange(0.2, rate=0),
            msprime.MigrationRateChange(0.4, matrix_index=(0, 1), rate=1),
            msprime.MigrationRateChange(0.4, matrix_index=(1, 0), rate=1)]
        self.verify_debug(
            population_configurations, migration_matrix, demographic_events)
