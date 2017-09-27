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
Test cases for the HDF5 format in msprime.
"""
from __future__ import print_function
from __future__ import division

import contextlib
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
import numpy as np

import msprime
import _msprime


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
    ts = msprime.simulate(
        10, recombination_rate=1, length=10, random_seed=2)
    sites = [msprime.Site(
        position=j, index=j, ancestral_state="0", mutations=[
            msprime.Mutation(site=j, node=k, derived_state="1") for k in range(j + 1)])
        for j in range(ts.sample_size)]
    return ts.copy(sites)


def general_mutation_example():
    ts = msprime.simulate(10, recombination_rate=1, length=10, random_seed=2)
    nodes = msprime.NodeTable()
    edges = msprime.EdgeTable()
    ts.dump_tables(nodes=nodes, edges=edges)
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    sites.add_row(position=0, ancestral_state="A")
    sites.add_row(position=1, ancestral_state="C")
    mutations.add_row(site=0, node=0, derived_state="T")
    mutations.add_row(site=1, node=0, derived_state="G")
    return msprime.load_tables(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations)


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


def node_name_example():
    ts = msprime.simulate(
        sample_size=100, recombination_rate=0.1, length=10, random_seed=1)
    nodes = msprime.NodeTable()
    edges = msprime.EdgeTable()
    ts.dump_tables(nodes=nodes, edges=edges)
    new_nodes = msprime.NodeTable()
    names = ["n_{}".format(u) for u in range(ts.num_nodes)]
    packed, length = msprime.pack_strings(names)
    new_nodes.set_columns(
        name=packed, name_length=length, flags=nodes.flags, time=nodes.time)
    return msprime.load_tables(
        nodes=new_nodes, edges=edges, provenance_strings=[b"sdf"])


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
            l, r = t.interval
            self.assertGreater(r, l)
            for site in t.sites():
                self.assertTrue(l <= site.position < r)
                for mut in site.mutations:
                    self.assertEqual(mut.site, site.index)

    def test_msprime_v_0_4_0(self):
        ts = msprime.load_legacy("tests/data/hdf5-formats/msprime-0.4.0_v3.1.hdf5")
        self.verify_tree_sequence(ts)

    def test_msprime_v_0_3_0(self):
        ts = msprime.load_legacy("tests/data/hdf5-formats/msprime-0.3.0_v2.0.hdf5")
        self.verify_tree_sequence(ts)


@unittest.skip("Diffs broken")
class TestRoundTrip(TestHdf5):
    """
    Tests if we can round trip convert a tree sequence in memory
    through a V2 file format and a V3 format.
    """
    def verify_tree_sequences_equal(self, ts, tsp):
        self.assertEqual(ts.num_edges, tsp.num_edges)
        self.assertEqual(ts.get_sample_size(), tsp.get_sample_size())
        self.assertEqual(ts.get_sequence_length(), tsp.get_sequence_length())
        self.assertEqual(ts.get_num_mutations(), tsp.get_num_mutations())
        self.assertEqual(ts.get_num_sites(), tsp.get_num_sites())
        self.assertEqual(ts.get_num_trees(), tsp.get_num_trees())
        self.assertEqual(list(ts.sites()), list(tsp.sites()))
        self.assertEqual(list(ts.nodes()), list(tsp.nodes()))
        self.assertEqual(list(ts.edges()), list(tsp.edges()))
        self.assertEqual(ts.get_num_nodes(), tsp.get_num_nodes())
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
            msprime.SimpleBottleneck(time=0.01, proportion=1)
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
        self.assertEqual(format_version[0], 7)
        self.assertEqual(format_version[1], 0)
        keys = set(root.keys())
        self.assertIn("nodes", keys)
        self.assertIn("edges", keys)
        self.assertIn("sites", keys)
        self.assertIn("mutations", keys)
        # Not filled in yet, but the group should be present for forward compatability.
        self.assertIn("migrations", keys)

        if "provenance"in keys:
            # TODO verify provenance
            pass
        g = root["sites"]
        fields = [
            ("position", float64), ("ancestral_state", int8),
            ("ancestral_state_length", uint32)]
        if ts.num_sites > 0:
            for name, dtype in fields:
                self.assertEqual(len(g[name].shape), 1)
                self.assertEqual(g[name].dtype, dtype)
            ancestral_state_length = g["ancestral_state_length"]
            self.assertEqual(ancestral_state_length.shape[0], ts.num_sites)
            position = list(g["position"])
            ancestral_states = msprime.unpack_strings(
                g["ancestral_state"], ancestral_state_length)
            for j, site in enumerate(ts.sites()):
                self.assertEqual(position[j], site.position)
                self.assertEqual(ancestral_states[j], site.ancestral_state)

        g = root["mutations"]
        fields = [
            ("site", int32), ("node", int32),
            ("derived_state", int8), ("derived_state_length", uint32)]
        if ts.num_mutations > 0:
            for name, dtype in fields:
                self.assertEqual(len(g[name].shape), 1)
                self.assertEqual(g[name].shape[0], ts.get_num_mutations())
                self.assertEqual(g[name].dtype, dtype)
            derived_state_length = g["derived_state_length"]
            self.assertEqual(derived_state_length.shape[0], ts.num_mutations)
            site = g["site"]
            node = g["node"]
            derived_state = msprime.unpack_strings(
                g["derived_state"], derived_state_length)
            j = 0
            for s in ts.sites():
                for mutation in s.mutations:
                    self.assertEqual(site[j], s.index)
                    self.assertEqual(mutation.site, site[j])
                    self.assertEqual(mutation.node, node[j])
                    self.assertEqual(mutation.derived_state, derived_state[j])
                    j += 1
        else:
            self.assertEqual(0, len(list(g.keys())))

        # TODO some of these fields should be optional.
        nodes_group = root["nodes"]
        self.assertEqual(nodes_group["flags"].dtype, uint32)
        self.assertEqual(nodes_group["population"].dtype, int32)
        self.assertEqual(nodes_group["time"].dtype, float64)
        self.assertEqual(nodes_group["name_length"].dtype, uint32)
        stored_fields = set(nodes_group.keys())
        total_name_length = sum(nodes_group["name_length"])
        if "name" in stored_fields:
            self.assertEqual(nodes_group["name"].dtype, int8)
            self.assertGreater(total_name_length, 0)
        else:
            self.assertEqual(total_name_length, 0)
        population = [0 for _ in range(ts.get_num_nodes())]
        time = [0 for _ in range(ts.get_num_nodes())]
        # FIXME this will break when we have samples not in 0...n-1
        flags = [0 for _ in range(ts.get_num_nodes())]
        for u in range(ts.get_sample_size()):
            time[u] = ts.get_time(u)
            population[u] = ts.get_population(u)
            flags[u] = 1  # FIXME
        for i, node in enumerate(ts.nodes()):
            # TODO - so this means we assume nodes are densely numbered in
            # increasing order
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
        I = sorted(
            range(ts.num_edges),
            key=lambda j: (left[j], time[parent[j]]))
        O = sorted(
            range(ts.num_edges),
            key=lambda j: (right[j], -time[parent[j]], -child[j]))
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

    def test_node_names_example(self):
        self.verify_tree_dump_format(node_name_example())

    def test_optional_provenance(self):
        ts = single_locus_no_mutation_example()
        ts.dump(self.temp_file)
        hfile = h5py.File(self.temp_file, "r+")
        del hfile["provenance"]
        hfile.close()
        del hfile
        other_ts = msprime.load(self.temp_file)
        self.assertEqual(other_ts.get_provenance(), [])
        self.verify_tree_dump_format(other_ts)

    def test_general_mutation_example(self):
        self.verify_tree_dump_format(general_mutation_example())


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
