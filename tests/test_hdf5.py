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

try:
    import itertools.izip as zip
except ImportError:
    pass

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
            msprime.SimpleBottleneck(time=0.01, proportion=0.75)])


def historical_sample_example():
    return msprime.simulate(
        recombination_rate=0.1,
        length=10,
        random_seed=1,
        samples=[(0, j) for j in range(10)])


class TestHdf5(unittest.TestCase):
    """
    Superclass of HDF5 tests.
    """
    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msp_hdf5_test_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)


class TestRoundTrip(TestHdf5):
    """
    Tests if we can round trip convert a tree sequence in memory
    through a V2 file format and a V3 format.
    """
    def verify_tree_sequences_equal(self, ts, tsp):
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
            self.assertEqual(m1.nodes, m2.nodes)
            j += 1
        self.assertEqual(ts.get_num_nodes(), tsp.get_num_nodes())
        for u in range(ts.get_sample_size()):
            self.assertEqual(ts.get_population(u), tsp.get_population(u))
        self.assertEqual(ts.num_trees, tsp.num_trees)
        num_trees = 0
        for t1, t2 in zip(ts.trees(), tsp.trees()):
            self.assertEqual(t1.parent_dict, t2.parent_dict)
            num_trees += 1
        self.assertEqual(num_trees, ts.num_trees)

        provenance = tsp.get_provenance()
        self.assertGreater(len(provenance), 1)
        for p in provenance:
            self.assertIsInstance(json.loads(p.decode()), dict)

    def verify_round_trip(self, ts, version):
        tmp = sys.stderr
        try:
            with open(os.devnull, "w") as devnull:
                sys.stderr = devnull
                # We silence stderr here because h5py dumps out some
                # spurious # error messages. See
                # https://github.com/h5py/h5py/issues/390
                msprime.dump_legacy(ts, self.temp_file, version=version)
                tsp = msprime.load_legacy(self.temp_file)
        finally:
            sys.stderr = tmp
        self.verify_tree_sequences_equal(ts, tsp)

    def verify_malformed_json_v2(self, ts, group_name, attr, bad_json):
        tmp = sys.stderr
        try:
            with open(os.devnull, "w") as devnull:
                sys.stderr = devnull
                # We silence stderr here because h5py dumps out some
                # spurious error messages. See
                # https://github.com/h5py/h5py/issues/390
                msprime.dump_legacy(ts, self.temp_file, 2)
                # Write some bad JSON to the provenance string.
                root = h5py.File(self.temp_file, "r+")
                group = root[group_name]
                group.attrs[attr] = bad_json
                root.close()
                tsp = msprime.load_legacy(self.temp_file)
        finally:
            sys.stderr = tmp
        self.verify_tree_sequences_equal(ts, tsp)

    def test_malformed_json_v2(self):
        ts = multi_locus_with_mutation_example()
        for group_name in ["trees", "mutations"]:
            for attr in ["environment", "parameters"]:
                for bad_json in ["", "{", "{},"]:
                    self.verify_malformed_json_v2(ts, group_name, attr, bad_json)

    def test_single_locus_no_mutation(self):
        self.verify_round_trip(single_locus_no_mutation_example(), 2)
        self.verify_round_trip(single_locus_no_mutation_example(), 3)

    def test_single_locus_with_mutation(self):
        self.verify_round_trip(single_locus_with_mutation_example(), 2)
        self.verify_round_trip(single_locus_with_mutation_example(), 3)

    def test_multi_locus_with_mutation(self):
        self.verify_round_trip(multi_locus_with_mutation_example(), 2)
        self.verify_round_trip(multi_locus_with_mutation_example(), 3)

    def test_migration_example(self):
        self.verify_round_trip(migration_example(), 2)
        self.verify_round_trip(migration_example(), 3)

    def test_bottleneck_example(self):
        self.verify_round_trip(migration_example(), 3)


class TestErrors(TestHdf5):
    """
    Test various API errors.
    """
    def test_v2_non_binary_records(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=0.01, proportion=1)
        ]
        ts = msprime.simulate(
            sample_size=10,
            demographic_events=demographic_events,
            random_seed=1)
        self.assertRaises(ValueError, msprime.dump_legacy, ts, self.temp_file, 2)

    def test_unsupported_format(self):
        ts = msprime.simulate(10)
        self.assertRaises(ValueError, msprime.dump_legacy, ts, self.temp_file, version=4)
        # We refuse to read current version also
        ts.dump(self.temp_file)
        self.assertRaises(ValueError, msprime.load_legacy, self.temp_file)


class TestHdf5Format(TestHdf5):
    """
    Tests on the HDF5 file format.
    """

    def verify_tree_dump_format(self, ts):
        uint32 = "<u4"
        float64 = "<f8"
        ts.dump(self.temp_file)
        self.assertTrue(os.path.exists(self.temp_file))
        self.assertGreater(os.path.getsize(self.temp_file), 0)
        root = h5py.File(self.temp_file, "r")
        # Check the basic root attributes
        format_version = root.attrs['format_version']
        self.assertEqual(format_version[0], 4)
        self.assertEqual(format_version[1], 0)
        keys = set(root.keys())
        self.assertLessEqual(keys, set(["mutations", "trees", "provenance"]))
        self.assertIn("trees", keys)
        self.assertIn("provenance", keys)
        self.assertIn("mutations", keys)
        g = root["mutations"]
        fields = [("nodes", uint32), ("num_nodes", uint32), ("position", float64)]
        if ts.num_mutations > 0:
            self.assertEqual(set(g.keys()), set([name for name, _ in fields]))
            for name, dtype in fields:
                self.assertEqual(len(g[name].shape), 1)
                self.assertEqual(g[name].shape[0], ts.get_num_mutations())
                self.assertEqual(g[name].dtype, dtype)
            flat_nodes = list(g["nodes"])
            num_nodes = list(g["num_nodes"])
            position = list(g["position"])
            nodes = []
            offset = 0
            for k in num_nodes:
                nodes.append(tuple(flat_nodes[offset: offset + k]))
                offset += k
            self.assertEqual(len(num_nodes), ts.get_num_mutations())
            self.assertEqual(len(position), ts.get_num_mutations())
            for j, mutation in enumerate(ts.mutations()):
                self.assertEqual(mutation.nodes, nodes[j])
                self.assertEqual(mutation.position, position[j])
        else:
            self.assertEqual(0, len(list(g.keys())))

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
        self.assertEqual(nodes_group["population"].dtype, uint32)
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
        self.verify_tree_dump_format(single_locus_no_mutation_example())

    def test_single_locus_with_mutation(self):
        self.verify_tree_dump_format(single_locus_with_mutation_example())

    def test_multi_locus_with_mutation(self):
        self.verify_tree_dump_format(multi_locus_with_mutation_example())

    def test_migration_example(self):
        self.verify_tree_dump_format(migration_example())

    def test_bottleneck_example(self):
        self.verify_tree_dump_format(bottleneck_example())

    def test_historical_sample_example(self):
        self.verify_tree_dump_format(historical_sample_example())

    def test_optional_provenance(self):
        ts = single_locus_no_mutation_example()
        ts.dump(self.temp_file)
        hfile = h5py.File(self.temp_file, "r+")
        del hfile["provenance"]
        hfile.close()
        del hfile
        other_ts = msprime.load(self.temp_file)
        self.assertEqual(other_ts.get_provenance(), [])


class TestHdf5FormatErrors(TestHdf5):
    """
    Tests for errors in the HDF5 format.
    """

    def verify_fields(self, ts):
        names = []

        def visit(name):
            # The only dataset we can delete on its own is provenance
            if name != "provenance":
                names.append(name)
        ts.dump(self.temp_file)
        hfile = h5py.File(self.temp_file, "r")
        hfile.visit(visit)
        hfile.close()
        # Delete each field in turn; this should cause a LibraryError
        for name in names:
            ts.dump(self.temp_file)
            hfile = h5py.File(self.temp_file, "r+")
            del hfile[name]
            hfile.close()
            del hfile
            self.assertRaises(_msprime.LibraryError, msprime.load, self.temp_file)

    def test_mandatory_fields_no_mutation(self):
        self.verify_fields(single_locus_no_mutation_example())

    def test_mandatory_fields_with_mutation(self):
        self.verify_fields(single_locus_with_mutation_example())

    def test_load_malformed_hdf5(self):
        hfile = h5py.File(self.temp_file, "w")
        # First try the empty hdf5 file.
        hfile.close()
        self.assertRaises(_msprime.LibraryError, msprime.load, self.temp_file)

    def test_version_load_error(self):
        ts = msprime.simulate(10)
        for bad_version in [(0, 1), (0, 8), (2, 0)]:
            ts.dump(self.temp_file)
            hfile = h5py.File(self.temp_file, "r+")
            hfile.attrs['format_version'] = bad_version
            hfile.close()
            other_ts = _msprime.TreeSequence()
            self.assertRaises(_msprime.LibraryError, other_ts.load, self.temp_file)

    def test_load_bad_formats(self):
        # try loading a bunch of files in various formats.
        # First, check the emtpy file.
        self.assertRaises(_msprime.LibraryError, msprime.load, self.temp_file)
        # Now some ascii text
        with open(self.temp_file, "wb") as f:
            f.write(b"Some ASCII text")
        self.assertRaises(_msprime.LibraryError, msprime.load, self.temp_file)
        # Now write 8k of random bytes
        with open(self.temp_file, "wb") as f:
            f.write(os.urandom(8192))
        self.assertRaises(_msprime.LibraryError, msprime.load, self.temp_file)
