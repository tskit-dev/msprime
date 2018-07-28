#
# Copyright (C) 2018 University of Oxford
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
Test cases for the simulate-from functionality.
"""
from __future__ import print_function
from __future__ import division

import unittest

import msprime
import tests.tsutil as tsutil


class TestBasicFunctionality(unittest.TestCase):
    """
    Basic tests for the from_ts argument for msprime.simulate.
    """
    def verify_from_tables(self, from_ts, final_ts):
        from_tables = from_ts.dump_tables()
        final_tables = final_ts.dump_tables()
        # Populations and individuals should be equal.
        self.assertEqual(from_tables.populations, final_tables.populations)
        self.assertEqual(from_tables.individuals, final_tables.individuals)
        # Other tables should be equal up to the from_tables.
        final_tables.nodes.truncate(len(from_tables.nodes))
        self.assertEqual(final_tables.nodes, from_tables.nodes)
        final_tables.edges.truncate(len(from_tables.edges))
        self.assertEqual(final_tables.edges, from_tables.edges)
        final_tables.sites.truncate(len(from_tables.sites))
        self.assertEqual(final_tables.sites, from_tables.sites)
        final_tables.mutations.truncate(len(from_tables.mutations))
        self.assertEqual(final_tables.mutations, from_tables.mutations)
        final_tables.provenances.truncate(len(from_tables.provenances))
        self.assertEqual(final_tables.provenances, from_tables.provenances)

    def verify_simulation_completed(self, ts):
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 1)

    def test_from_single_locus_decapitated(self):
        ts = msprime.simulate(10, random_seed=5)
        from_ts = tsutil.decapitate(ts, ts.num_edges // 2)
        final_ts = msprime.simulate(from_ts=from_ts, random_seed=2)
        self.verify_from_tables(from_ts, final_ts)
        self.verify_simulation_completed(final_ts)

    def test_single_locus_max_time(self):
        from_ts = msprime.simulate(20, __tmp_max_time=1, random_seed=5)
        self.assertGreater(max(tree.num_roots for tree in from_ts.trees()), 1)
        final_ts = msprime.simulate(from_ts=from_ts, random_seed=2)
        self.verify_from_tables(from_ts, final_ts)
        self.verify_simulation_completed(final_ts)

    def test_from_multi_locus_decapitated(self):
        ts = msprime.simulate(10, recombination_rate=2, random_seed=5)
        self.assertGreater(ts.num_trees, 1)
        from_ts = tsutil.decapitate(ts, ts.num_edges // 2)
        final_ts = msprime.simulate(from_ts=from_ts, random_seed=2)
        self.verify_from_tables(from_ts, final_ts)
        self.verify_simulation_completed(final_ts)

    def test_from_multi_locus_max_time(self):
        from_ts = msprime.simulate(
            10, recombination_rate=2, random_seed=5, __tmp_max_time=1)
        self.assertGreater(from_ts.num_trees, 1)
        final_ts = msprime.simulate(from_ts=from_ts, random_seed=2)
        self.verify_from_tables(from_ts, final_ts)
        self.verify_simulation_completed(final_ts)


class TestErrors(unittest.TestCase):
    """
    Basic tests for the from_ts argument for msprime.simulate.
    """
    def get_example_base(self):
        ts = msprime.simulate(10, random_seed=5)
        return tsutil.decapitate(ts, ts.num_edges // 2)

    def test_samples(self):
        base_ts = self.get_example_base()
        self.assertRaises(ValueError, msprime.simulate, 2, from_ts=base_ts)
        self.assertRaises(ValueError, msprime.simulate, sample_size=2, from_ts=base_ts)
        self.assertRaises(
            ValueError, msprime.simulate,
            samples=[msprime.Sample(0, 0) for _ in range(10)], from_ts=base_ts)
        self.assertRaises(
            ValueError, msprime.simulate,
            population_configurations=[
                msprime.PopulationConfiguration(sample_size=2)],
            from_ts=base_ts)

    def test_bad_type(self):
        for bad_type in [{}, 1, "asd"]:
            self.assertRaises(TypeError, msprime.simulate, from_ts=bad_type)
