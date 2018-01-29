#
# Copyright (C) 2016-2017 University of Oxford
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

import contextlib
import os
import sys
import tempfile
import unittest

import h5py
import numpy as np

import msprime
import _msprime
import tests.tsutil as tsutil


@contextlib.contextmanager
def silence_stderr():
    """
    Context manager to silence stderr. We do this here because h5py dumps out some
    spurious error messages. See https://github.com/h5py/h5py/issues/390
    """
    tmp = sys.stderr
    try:
        with open(os.devnull, "w") as devnull:
            sys.stderr = devnull
            yield
    finally:
        sys.stderr = tmp


def single_locus_no_mutation_example():
    return msprime.simulate(10, random_seed=10)


def single_locus_with_mutation_example():
    return msprime.simulate(10, mutation_rate=10, random_seed=11)


def multi_locus_with_mutation_example():
    return msprime.simulate(
        10, recombination_rate=1, length=10, mutation_rate=10,
        random_seed=2)


def recurrent_mutation_example():
    ts = msprime.simulate(10, recombination_rate=1, length=10, random_seed=2)
    return tsutil.insert_branch_mutations(ts)


def general_mutation_example():
    ts = msprime.simulate(10, recombination_rate=1, length=10, random_seed=2)
    nodes = msprime.NodeTable()
    edges = msprime.EdgeTable()
    ts.dump_tables(nodes=nodes, edges=edges)
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    sites.add_row(position=0, ancestral_state="A", metadata=b"{}")
    sites.add_row(position=1, ancestral_state="C", metadata=b"{'id':1}")
    mutations.add_row(site=0, node=0, derived_state="T")
    mutations.add_row(site=1, node=0, derived_state="G")
    return msprime.load_tables(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations)


def multichar_mutation_example():
    ts = msprime.simulate(10, recombination_rate=1, length=10, random_seed=2)
    return tsutil.insert_multichar_mutations(ts)


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
            msprime.SimpleBottleneck(time=0.01, population=0, proportion=0.75)])


def historical_sample_example():
    return msprime.simulate(
        recombination_rate=0.1,
        length=10,
        random_seed=1,
        samples=[(0, j) for j in range(10)])


def node_metadata_example():
    ts = msprime.simulate(
        sample_size=100, recombination_rate=0.1, length=10, random_seed=1)
    nodes = msprime.NodeTable()
    edges = msprime.EdgeTable()
    ts.dump_tables(nodes=nodes, edges=edges)
    new_nodes = msprime.NodeTable()
    metadatas = ["n_{}".format(u) for u in range(ts.num_nodes)]
    packed, offset = msprime.pack_strings(metadatas)
    new_nodes.set_columns(
        metadata=packed, metadata_offset=offset, flags=nodes.flags, time=nodes.time)
    return msprime.load_tables(nodes=new_nodes, edges=edges)


class TestHdf5(unittest.TestCase):
    """
    Superclass of HDF5 tests.
    """
    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msp_hdf5_test_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)


class TestLoadLegacyExamples(TestHdf5):
    """
    Tests using the saved legacy file examples to ensure we can load them.
    """
    def verify_tree_sequence(self, ts):
        # Just some quick checks to make sure the tree sequence makes sense.
        self.assertGreater(ts.sample_size, 0)
        self.assertGreater(ts.num_edges, 0)
        self.assertGreater(ts.num_sites, 0)
        self.assertGreater(ts.num_mutations, 0)
        self.assertGreater(ts.sequence_length, 0)
        for t in ts.trees():
            left, right = t.interval
            self.assertGreater(right, left)
            for site in t.sites():
                self.assertTrue(left <= site.position < right)
                for mut in site.mutations:
                    self.assertEqual(mut.site, site.id)

    def test_msprime_v_0_4_0(self):
        ts = msprime.load_legacy("tests/data/hdf5-formats/msprime-0.4.0_v3.1.hdf5")
        self.verify_tree_sequence(ts)

    def test_msprime_v_0_3_0(self):
        ts = msprime.load_legacy("tests/data/hdf5-formats/msprime-0.3.0_v2.0.hdf5")
        self.verify_tree_sequence(ts)


class TestRoundTrip(TestHdf5):
    """
    Tests if we can round trip convert a tree sequence in memory
    through a V2 file format and a V3 format.
    """
    def verify_tree_sequences_equal(self, ts, tsp):
        t1 = ts.tables
        # We need to sort and squash the edges in the new format because it
        # has gone through an edgesets representation. Simplest way to do this
        # is to call simplify.
        t2 = tsp.simplify().tables
        self.assertEqual(t1.nodes, t2.nodes)
        self.assertEqual(t1.edges, t2.edges)
        self.assertEqual(t1.sites, t2.sites)
        self.assertEqual(t1.mutations, t2.mutations)

    def verify_round_trip(self, ts, version):
        msprime.dump_legacy(ts, self.temp_file, version=version)
        with silence_stderr():
            tsp = msprime.load_legacy(self.temp_file)
        self.verify_tree_sequences_equal(ts, tsp)
        tsp.dump(self.temp_file)
        tsp = msprime.load(self.temp_file)
        self.verify_tree_sequences_equal(ts, tsp)

    def verify_malformed_json_v2(self, ts, group_name, attr, bad_json):
        msprime.dump_legacy(ts, self.temp_file, 2)
        # Write some bad JSON to the provenance string.
        root = h5py.File(self.temp_file, "r+")
        group = root[group_name]
        group.attrs[attr] = bad_json
        root.close()
        with silence_stderr():
            tsp = msprime.load_legacy(self.temp_file)
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

    def test_recurrent_mutation_example(self):
        ts = recurrent_mutation_example()
        for version in [2, 3]:
            self.assertRaises(
                ValueError, msprime.dump_legacy, ts, self.temp_file, version)

    def test_general_mutation_example(self):
        ts = general_mutation_example()
        for version in [2, 3]:
            self.assertRaises(
                ValueError, msprime.dump_legacy, ts, self.temp_file, version)

    def test_v2_no_samples(self):
        ts = multi_locus_with_mutation_example()
        msprime.dump_legacy(ts, self.temp_file, version=2)
        root = h5py.File(self.temp_file, "r+")
        del root['samples']
        root.close()
        with silence_stderr():
            tsp = msprime.load_legacy(self.temp_file)
        self.verify_tree_sequences_equal(ts, tsp)

    def test_duplicate_mutation_positions_single_value(self):
        ts = multi_locus_with_mutation_example()
        for version in [2, 3]:
            msprime.dump_legacy(ts, self.temp_file, version=version)
            root = h5py.File(self.temp_file, "r+")
            root['mutations/position'][:] = 0
            root.close()
            with silence_stderr():
                self.assertRaises(
                    msprime.DuplicatePositionsError, msprime.load_legacy, self.temp_file)
                tsp = msprime.load_legacy(
                    self.temp_file, remove_duplicate_positions=True)
            self.assertEqual(tsp.num_sites, 1)
            sites = list(tsp.sites())
            self.assertEqual(sites[0].position, 0)

    def test_duplicate_mutation_positions(self):
        ts = multi_locus_with_mutation_example()
        for version in [2, 3]:
            msprime.dump_legacy(ts, self.temp_file, version=version)
            root = h5py.File(self.temp_file, "r+")
            position = np.array(root['mutations/position'])
            position[0] = position[1]
            root['mutations/position'][:] = position
            root.close()
            with silence_stderr():
                self.assertRaises(
                    msprime.DuplicatePositionsError, msprime.load_legacy, self.temp_file)
                tsp = msprime.load_legacy(
                    self.temp_file, remove_duplicate_positions=True)
            self.assertEqual(tsp.num_sites, position.shape[0] - 1)
            position_after = list(s.position for s in tsp.sites())
            self.assertEqual(list(position[1:]), position_after)


class TestErrors(TestHdf5):
    """
    Test various API errors.
    """
    def test_v2_non_binary_records(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=0.01, population=0, proportion=1)
        ]
        ts = msprime.simulate(
            sample_size=10,
            demographic_events=demographic_events,
            random_seed=1)
        self.assertRaises(ValueError, msprime.dump_legacy, ts, self.temp_file, 2)

    def test_unsupported_version(self):
        ts = msprime.simulate(10)
        self.assertRaises(ValueError, msprime.dump_legacy, ts, self.temp_file, version=4)
        # We refuse to read current version also
        ts.dump(self.temp_file)
        self.assertRaises(ValueError, msprime.load_legacy, self.temp_file)

    def test_no_version_number(self):
        root = h5py.File(self.temp_file, "w")
        root.attrs["x"] = 0
        root.close()
        self.assertRaises(ValueError, msprime.load_legacy, self.temp_file)


class TestHdf5Format(TestHdf5):
    """
    Tests on the HDF5 file format.
    """

    def verify_metadata(self, group, num_rows):

        self.assertEqual(group["metadata_offset"].dtype, np.uint32)
        metadata_offset = list(group["metadata_offset"])
        metadata_length = 0
        if metadata_offset[-1] > 0:
            self.assertEqual(group["metadata"].dtype, np.int8)
            metadata = list(group["metadata"])
            metadata_length = len(metadata)
            self.assertEqual(metadata_offset[-1], metadata_length)
        else:
            self.assertNotIn("metadata", group)
        self.assertEqual(len(metadata_offset), num_rows + 1)

    def verify_tree_dump_format(self, ts):
        int8 = "<i1"
        int32 = "<i4"
        uint32 = "<u4"
        float64 = "<f8"
        ts.dump(self.temp_file)
        self.assertTrue(os.path.exists(self.temp_file))
        self.assertGreater(os.path.getsize(self.temp_file), 0)
        root = h5py.File(self.temp_file, "r")
        # Check the basic root attributes
        format_version = root.attrs['format_version']
        self.assertEqual(format_version[0], 10)
        self.assertEqual(format_version[1], 0)
        sequence_length = root.attrs['sequence_length']
        self.assertGreater(sequence_length, 0)
        keys = set(root.keys())
        self.assertIn("nodes", keys)
        self.assertIn("edges", keys)
        self.assertIn("sites", keys)
        self.assertIn("mutations", keys)
        self.assertIn("provenances", keys)
        # Not filled in yet, but the group should be present for forward compatability.
        self.assertIn("migrations", keys)

        g = root["sites"]
        fields = [
            ("position", float64), ("ancestral_state", int8),
            ("ancestral_state_offset", uint32)]
        self.verify_metadata(g, ts.num_sites)
        ancestral_state_offset = g["ancestral_state_offset"]
        if ts.num_sites == 0:
            self.assertEqual(ancestral_state_offset.shape, (1,))
            self.assertNotIn("ancestral_state", list(g.keys()))
        else:
            for name, dtype in fields:
                self.assertEqual(len(g[name].shape), 1)
                self.assertEqual(g[name].dtype, dtype)
            position = list(g["position"])
            self.assertEqual(len(position), ts.num_sites)
            self.assertEqual(len(ancestral_state_offset), ts.num_sites + 1)
            ancestral_states = msprime.unpack_strings(
                g["ancestral_state"], ancestral_state_offset)
            for j, site in enumerate(ts.sites()):
                self.assertEqual(position[j], site.position)
                self.assertEqual(ancestral_states[j], site.ancestral_state)

        g = root["mutations"]
        fields = [
            ("site", int32), ("node", int32), ("parent", int32),
            ("derived_state", int8), ("derived_state_offset", uint32)]
        derived_state_offset = g["derived_state_offset"]
        self.verify_metadata(g, ts.num_sites)
        if ts.num_mutations == 0:
            self.assertEqual(derived_state_offset.shape, (1,))
            self.assertNotIn("derived_state", list(g.keys()))
        else:
            for name, dtype in fields:
                self.assertEqual(len(g[name].shape), 1)
                self.assertEqual(g[name].dtype, dtype)
            self.assertEqual(derived_state_offset.shape[0], ts.num_mutations + 1)
            site = g["site"]
            node = g["node"]
            parent = g["parent"]
            for col in [site, node, parent]:
                self.assertEqual(col.shape[0], ts.num_mutations)
            derived_state = msprime.unpack_strings(
                g["derived_state"], derived_state_offset)
            j = 0
            for s in ts.sites():
                for mutation in s.mutations:
                    self.assertEqual(site[j], s.id)
                    self.assertEqual(mutation.site, site[j])
                    self.assertEqual(mutation.node, node[j])
                    self.assertEqual(mutation.parent, parent[j])
                    self.assertEqual(mutation.derived_state, derived_state[j])
                    j += 1

        # TODO some of these fields should be optional.
        nodes_group = root["nodes"]
        self.assertEqual(nodes_group["flags"].dtype, uint32)
        self.assertEqual(nodes_group["population"].dtype, int32)
        self.assertEqual(nodes_group["time"].dtype, float64)
        self.verify_metadata(nodes_group, ts.num_nodes)
        population = [0 for _ in range(ts.num_nodes)]
        time = [0 for _ in range(ts.num_nodes)]
        flags = [0 for _ in range(ts.num_nodes)]
        for i, node in enumerate(ts.nodes()):
            flags[i] = int(node.is_sample())
            time[i] = node.time
            population[i] = node.population
        self.assertEqual(time, list(nodes_group["time"]))
        self.assertEqual(population, list(nodes_group["population"]))
        self.assertEqual(flags, list(nodes_group["flags"]))

        edges_group = root["edges"]
        self.assertEqual(
            set(edges_group.keys()),
            {"indexes", "child", "left", "right", "parent"})

        self.assertEqual(edges_group["left"].dtype, float64)
        self.assertEqual(edges_group["right"].dtype, float64)
        self.assertEqual(edges_group["parent"].dtype, int32)
        self.assertEqual(edges_group["child"].dtype, int32)
        left = list(edges_group["left"])
        right = list(edges_group["right"])
        parent = list(edges_group["parent"])
        child = list(edges_group["child"])

        self.assertEqual(len(left), ts.num_edges)
        self.assertEqual(len(right), ts.num_edges)
        self.assertEqual(len(parent), ts.num_edges)
        self.assertEqual(len(child), ts.num_edges)
        for j, record in enumerate(ts.edges()):
            self.assertEqual(record.left, left[j])
            self.assertEqual(record.right, right[j])
            self.assertEqual(record.parent, parent[j])
            self.assertEqual(record.child, child[j])

        indexes_group = edges_group["indexes"]
        self.assertEqual(
            set(indexes_group.keys()), {"insertion_order", "removal_order"})
        for field in indexes_group.keys():
            self.assertEqual(indexes_group[field].dtype, int32, child[j])
        in_order = sorted(
            range(ts.num_edges),
            key=lambda j: (left[j], time[parent[j]]))
        out_order = sorted(
            range(ts.num_edges),
            key=lambda j: (right[j], -time[parent[j]], -child[j]))
        self.assertEqual(in_order, list(indexes_group["insertion_order"]))
        self.assertEqual(out_order, list(indexes_group["removal_order"]))

        provenances_group = root["provenances"]
        timestamp_offset = list(provenances_group["timestamp_offset"])
        record_offset = list(provenances_group["record_offset"])
        self.assertEqual(provenances_group["timestamp_offset"].dtype, uint32)
        self.assertEqual(provenances_group["record_offset"].dtype, uint32)
        self.assertEqual(len(timestamp_offset), ts.num_provenances + 1)
        self.assertEqual(len(record_offset), ts.num_provenances + 1)
        if timestamp_offset[-1] > 0:
            timestamp = list(provenances_group["timestamp"])
            self.assertEqual(provenances_group["timestamp"].dtype, int8)
            self.assertEqual(len(timestamp), timestamp_offset[-1])
        if record_offset[-1] > 0:
            record = list(provenances_group["record"])
            self.assertEqual(provenances_group["record"].dtype, int8)
            self.assertEqual(len(record), record_offset[-1])

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

    def test_node_metadata_example(self):
        self.verify_tree_dump_format(node_metadata_example())

    def test_general_mutation_example(self):
        self.verify_tree_dump_format(general_mutation_example())

    def test_multichar_mutation_example(self):
        self.verify_tree_dump_format(multichar_mutation_example())


class TestHdf5FormatErrors(TestHdf5):
    """
    Tests for errors in the HDF5 format.
    """

    def verify_fields(self, ts):
        names = []

        def visit(name):
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
