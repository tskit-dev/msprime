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
Test cases for format conversions in msprime.
"""
from __future__ import print_function
from __future__ import division

import json
import tempfile
import unittest

import msprime


class TestRoundTrip(unittest.TestCase):
    """
    Tests if we can round trip convert a tree sequence in memory
    through a V2 file format and a V3 format.
    """

    def verify_round_trip(self, ts):
        with tempfile.NamedTemporaryFile(prefix="msp_ff_") as f:
            msprime.dump_legacy(ts, f.name)
            tsp = msprime.load_legacy(f.name)
        self.assertEqual(ts.get_num_records(), tsp.get_num_records())
        self.assertEqual(ts.get_sample_size(), tsp.get_sample_size())
        self.assertEqual(ts.get_sequence_length(), tsp.get_sequence_length())
        self.assertEqual(ts.get_num_mutations(), tsp.get_num_mutations())
        j = 0
        for r1, r2 in zip(ts.records(), tsp.records()):
            self.assertEqual(r1.left, r2.left)
            self.assertEqual(r1.right, r2.right)
            self.assertEqual(r1.node, r2.node)
            self.assertEqual(r1.time, r2.time)
            self.assertEqual(r1.children, r2.children)
            self.assertEqual(r1.population, r2.population)
            j += 1
        self.assertEqual(j, ts.get_num_records())
        self.assertEqual(ts.get_num_trees(), tsp.get_num_trees())
        j = 0
        for m1, m2 in zip(ts.mutations(), tsp.mutations()):
            self.assertEqual(m1.position, m2.position)
            self.assertEqual(m1.node, m2.node)
            j += 1
        self.assertEqual(ts.get_num_nodes(), tsp.get_num_nodes())
        for u in range(ts.get_sample_size()):
            self.assertEqual(ts.get_population(u), tsp.get_population(u))
        provenance = tsp.get_provenance()
        if ts.get_num_mutations() > 0:
            self.assertEqual(len(provenance), 3)
        else:
            self.assertEqual(len(provenance), 2)
        for p in provenance:
            self.assertIsInstance(json.loads(p), dict)

    def test_single_locus_no_mutation(self):
        ts = msprime.simulate(10)
        self.verify_round_trip(ts)

    def test_single_locus_with_mutation(self):
        ts = msprime.simulate(10, mutation_rate=10)
        self.verify_round_trip(ts)

    def test_multi_locus_with_mutation(self):
        ts = msprime.simulate(
            10, recombination_rate=1, length=10, mutation_rate=10)
        self.verify_round_trip(ts)

    def test_migration_example(self):
        n = 10
        t = 100
        population_configurations = [
            msprime.PopulationConfiguration(n // 2),
            msprime.PopulationConfiguration(n // 2),
            msprime.PopulationConfiguration(0),
        ]
        demographic_events = [
            msprime.MassMigration(time=t, source=0, destination=2),
            msprime.MassMigration(time=t, source=1, destination=2),
        ]
        ts = msprime.simulate(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            random_seed=1)
        self.verify_round_trip(ts)


class TestErrors(unittest.TestCase):
    """
    Test various API errors.
    """
    def test_v2_non_binary_records(self):
        demographic_events = [
            msprime.Bottleneck(time=0.01, proportion=1)
        ]
        ts = msprime.simulate(
            sample_size=10,
            demographic_events=demographic_events,
            random_seed=1)
        with tempfile.NamedTemporaryFile() as f:
            self.assertRaises(
                ValueError, msprime.dump_legacy, ts, f.name)

    def test_unsupported_format(self):
        ts = msprime.simulate(10)
        with tempfile.NamedTemporaryFile() as f:
            self.assertRaises(
                ValueError, msprime.dump_legacy, ts, f.name, version=3)
            # We refuse to read version 3 also
            # TODO FIXME
            # ts.dump(f.name, skip_h5close=True)
            # self.assertRaises(ValueError, msprime.load_legacy, f.name)
