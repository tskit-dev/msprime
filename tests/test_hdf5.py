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
Test cases for the HDF5 format in msprime.
"""
from __future__ import print_function
from __future__ import division

import json
import os
import sys
import tempfile
import unittest

import h5py

import msprime
import _msprime


def single_locus_no_mutation_example():
    return msprime.simulate(10, random_seed=10)


def single_locus_with_mutation_example():
    return msprime.simulate(10, mutation_rate=10, random_seed=11)


def multi_locus_with_mutation_example():
    return msprime.simulate(
        10, recombination_rate=1, length=10, mutation_rate=10,
        random_seed=2)


def migration_example():
    n = 10
    t = 1
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
    return ts


def bottleneck_example():
    return msprime.simulate(
        random_seed=1,
        sample_size=100,
        recombination_rate=0.1,
        length=10,
        demographic_events=[
            msprime.Bottleneck(time=0.01, proportion=0.75)])


def historical_sample_example():
    return msprime.simulate(
        recombination_rate=0.1,
        length=10,
        random_seed=1,
        samples=[(0, j) for j in range(10)])


class TestRoundTrip(unittest.TestCase):
    """
    Tests if we can round trip convert a tree sequence in memory
    through a V2 file format and a V3 format.
    """

    def verify_round_trip(self, ts):
        with tempfile.NamedTemporaryFile(prefix="msp_ff_") as f:
            tmp = sys.stderr
            try:
                with open("/dev/null", "w") as devnull:
                    sys.stderr = devnull
                    # We silence stderr here because h5py dumps out some
                    # spurious # error messages. See
                    # https://github.com/h5py/h5py/issues/390
                    msprime.dump_legacy(ts, f.name)
                    tsp = msprime.load_legacy(f.name)
            finally:
                sys.stderr = tmp
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
        self.verify_round_trip(single_locus_no_mutation_example())

    def test_single_locus_with_mutation(self):
        self.verify_round_trip(single_locus_with_mutation_example())

    def test_multi_locus_with_mutation(self):
        self.verify_round_trip(multi_locus_with_mutation_example())

    def test_migration_example(self):
        self.verify_round_trip(migration_example())


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
            ts.dump(f.name)
            self.assertRaises(ValueError, msprime.load_legacy, f.name)


class TestHdf5Format(unittest.TestCase):
    """
    Tests on the HDF5 file format.
    """

    def verify_tree_dump_format(self, ts, outfile):
        uint8 = "uint8"
        uint32 = "uint32"
        float64 = "float64"
        ts.dump(outfile.name)
        self.assertTrue(os.path.exists(outfile.name))
        self.assertGreater(os.path.getsize(outfile.name), 0)
        root = h5py.File(outfile.name, "r")
        # Check the basic root attributes
        format_version = root.attrs['format_version']
        self.assertEqual(format_version[0], 3)
        self.assertEqual(format_version[1], 1)
        keys = set(root.keys())
        self.assertLessEqual(keys, set(["mutations", "trees", "provenance"]))
        self.assertIn("trees", keys)
        self.assertIn("provenance", keys)
        if ts.get_num_mutations() == 0:
            self.assertNotIn("mutations", keys)
        else:
            self.assertIn("mutations", keys)
            g = root["mutations"]
            fields = [("node", uint32), ("position", float64)]
            self.assertEqual(set(g.keys()), set([name for name, _ in fields]))
            for name, dtype in fields:
                self.assertEqual(len(g[name].shape), 1)
                self.assertEqual(g[name].shape[0], ts.get_num_mutations())
                self.assertEqual(g[name].dtype, dtype)
            node = list(g["node"])
            position = list(g["position"])
            self.assertEqual(len(node), ts.get_num_mutations())
            self.assertEqual(len(position), ts.get_num_mutations())
            for j, mutation in enumerate(ts.mutations()):
                self.assertEqual(mutation.node, node[j])
                self.assertEqual(mutation.position, position[j])

        trees_group = root["trees"]
        self.assertEqual(
            set(trees_group.keys()),
            {"records", "nodes", "breakpoints", "indexes"})
        # Breakpoints must be equal to the sorted list of distinct left and
        # right values.
        breakpoints = set()
        for record in ts.records():
            breakpoints.add(record.left)
            breakpoints.add(record.right)
        breakpoints = sorted(list(breakpoints))
        self.assertEqual(list(trees_group["breakpoints"]), breakpoints)
        self.assertEqual(trees_group["breakpoints"].dtype, float64)

        nodes_group = trees_group["nodes"]
        self.assertEqual(set(nodes_group.keys()), {"population", "time"})
        self.assertEqual(nodes_group["population"].dtype, uint8)
        self.assertEqual(nodes_group["time"].dtype, float64)
        population = [0 for j in range(ts.get_num_nodes())]
        time = [0 for j in range(ts.get_num_nodes())]
        for u in range(ts.get_sample_size()):
            time[u] = ts.get_time(u)
            population[u] = ts.get_population(u)
        for record in ts.records():
            time[record.node] = record.time
            population[record.node] = record.population
        self.assertEqual(time, list(nodes_group["time"]))
        self.assertEqual(population, list(nodes_group["population"]))

        records_group = trees_group["records"]
        self.assertEqual(
            set(records_group.keys()),
            {"children", "num_children", "left", "right", "node"})
        for field in records_group.keys():
            self.assertEqual(records_group[field].dtype, uint32)
        left = list(records_group["left"])
        right = list(records_group["right"])
        children = list(records_group["children"])
        num_children = list(records_group["num_children"])
        node = list(records_group["node"])
        self.assertEqual(len(left), ts.get_num_records())
        self.assertEqual(len(right), ts.get_num_records())
        self.assertEqual(len(num_children), ts.get_num_records())
        self.assertEqual(len(node), ts.get_num_records())
        breakpoint_map = {breakpoints[j]: j for j in range(len(breakpoints))}
        children_offset = 0
        for j, record in enumerate(ts.records()):
            self.assertEqual(time[record.node], record.time)
            self.assertEqual(population[record.node], record.population)
            self.assertEqual(breakpoint_map[record.left], left[j])
            self.assertEqual(breakpoint_map[record.right], right[j])
            self.assertEqual(
                list(record.children),
                children[children_offset: children_offset + num_children[j]])
            children_offset += num_children[j]

        indexes_group = trees_group["indexes"]
        self.assertEqual(
            set(indexes_group.keys()), {"insertion_order", "removal_order"})
        for field in indexes_group.keys():
            self.assertEqual(indexes_group[field].dtype, uint32)
        I = sorted(
            range(ts.get_num_records()),
            key=lambda j: (left[j], time[node[j]]))
        O = sorted(
            range(ts.get_num_records()),
            key=lambda j: (right[j], -time[node[j]]))
        self.assertEqual(I, list(indexes_group["insertion_order"]))
        self.assertEqual(O, list(indexes_group["removal_order"]))
        root.close()

    def test_single_locus_no_mutation(self):
        with tempfile.NamedTemporaryFile() as f:
            self.verify_tree_dump_format(
                single_locus_no_mutation_example(), f)

    def test_single_locus_with_mutation(self):
        with tempfile.NamedTemporaryFile() as f:
            self.verify_tree_dump_format(
                single_locus_with_mutation_example(), f)

    def test_multi_locus_with_mutation(self):
        with tempfile.NamedTemporaryFile() as f:
            self.verify_tree_dump_format(
                multi_locus_with_mutation_example(), f)

    def test_migration_example(self):
        with tempfile.NamedTemporaryFile() as f:
            self.verify_tree_dump_format(migration_example(), f)

    def test_bottleneck_example(self):
        with tempfile.NamedTemporaryFile() as f:
            self.verify_tree_dump_format(bottleneck_example(), f)

    def test_historical_sample_example(self):
        with tempfile.NamedTemporaryFile() as f:
            self.verify_tree_dump_format(historical_sample_example(), f)

    def test_optional_provenance(self):
        ts = single_locus_no_mutation_example()
        with tempfile.NamedTemporaryFile() as f:
            ts.dump(f.name)
            hfile = h5py.File(f.name, "r+")
            del hfile["provenance"]
            hfile.close()
            del hfile
            other_ts = msprime.load(f.name)
            self.assertEqual(other_ts.get_provenance(), [])


class TestHdf5FormatErrors(unittest.TestCase):
    """
    Tests for errors in the HDF5 format.
    """

    def verify_fields(self, ts):
        names = []

        def visit(name):
            if name not in ["provenance", "mutations"]:
                names.append(name)
        with tempfile.NamedTemporaryFile() as f:
            ts.dump(f.name)
            hfile = h5py.File(f.name, "r")
            hfile.visit(visit)
            hfile.close()
        for name in names:
            with tempfile.NamedTemporaryFile() as f:
                ts.dump(f.name)
                hfile = h5py.File(f.name, "r+")
                del hfile[name]
                hfile.close()
                del hfile
                self.assertRaises(_msprime.LibraryError, msprime.load, f.name)

    def test_mandatory_fields_no_mutation(self):
        self.verify_fields(single_locus_no_mutation_example())

    def test_mandatory_fields_with_mutation(self):
        self.verify_fields(single_locus_with_mutation_example())

    def test_load_malformed_hdf5(self):
        with tempfile.NamedTemporaryFile(prefix="malformed_") as f:
            hfile = h5py.File(f.name, "w")
            # First try the empty hdf5 file.
            hfile.close()
            self.assertRaises(
                _msprime.LibraryError, msprime.load, f.name)

    def test_version_load_error(self):
        ts = msprime.simulate(10)
        for bad_version in [(0, 1), (0, 8), (2, 0)]:
            with tempfile.NamedTemporaryFile() as f:
                ts.dump(f.name)
                hfile = h5py.File(f.name, "r+")
                hfile.attrs['format_version'] = bad_version
                hfile.close()
                other_ts = _msprime.TreeSequence()
                self.assertRaises(
                    _msprime.LibraryError, other_ts.load, f.name)

    def test_load_bad_formats(self):
        # try loading a bunch of files in various formats.
        with tempfile.NamedTemporaryFile("wb") as f:
            # First, check the emtpy file.
            self.assertRaises(_msprime.LibraryError, msprime.load, f.name)
            # Now some ascii text
            f.write(b"Some ASCII text")
            f.flush()
            self.assertRaises(_msprime.LibraryError, msprime.load, f.name)
            f.seek(0)
            # Now write 8k of random bytes
            f.write(os.urandom(8192))
            f.flush()
            self.assertRaises(_msprime.LibraryError, msprime.load, f.name)
