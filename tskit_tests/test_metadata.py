# -*- coding: utf-8 -*-
"""
Tests for metadata handling.
"""
from __future__ import print_function
from __future__ import division

import json
import os
import tempfile
import unittest
import pickle

import numpy as np
import python_jsonschema_objects as pjs
import six
import msprime

import tskit


class TestMetadataHdf5RoundTrip(unittest.TestCase):
    """
    Tests that we can encode metadata under various formats and this will
    successfully round-trip through the HDF5 format.
    """
    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msp_hdf5meta_test_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)

    def test_json(self):
        ts = msprime.simulate(10, random_seed=1)
        tables = ts.dump_tables()
        nodes = tables.nodes
        # For each node, we create some Python metadata that can be JSON encoded.
        metadata = [
            {"one": j, "two": 2 * j, "three": list(range(j))} for j in range(len(nodes))]
        encoded, offset = tskit.pack_strings(map(json.dumps, metadata))
        nodes.set_columns(
            flags=nodes.flags, time=nodes.time, population=nodes.population,
            metadata_offset=offset, metadata=encoded)
        self.assertTrue(np.array_equal(nodes.metadata_offset, offset))
        self.assertTrue(np.array_equal(nodes.metadata, encoded))
        ts1 = tables.tree_sequence()
        for j, node in enumerate(ts1.nodes()):
            decoded_metadata = json.loads(node.metadata.decode())
            self.assertEqual(decoded_metadata, metadata[j])
        ts1.dump(self.temp_file)
        ts2 = tskit.load(self.temp_file)
        self.assertEqual(ts1.tables.nodes, ts2.tables.nodes)

    def test_pickle(self):
        ts = msprime.simulate(10, random_seed=1)
        tables = ts.dump_tables()
        # For each node, we create some Python metadata that can be pickled
        metadata = [
            {"one": j, "two": 2 * j, "three": list(range(j))}
            for j in range(ts.num_nodes)]
        encoded, offset = tskit.pack_bytes(list(map(pickle.dumps, metadata)))
        tables.nodes.set_columns(
            flags=tables.nodes.flags, time=tables.nodes.time,
            population=tables.nodes.population,
            metadata_offset=offset, metadata=encoded)
        self.assertTrue(np.array_equal(tables.nodes.metadata_offset, offset))
        self.assertTrue(np.array_equal(tables.nodes.metadata, encoded))
        ts1 = tables.tree_sequence()
        for j, node in enumerate(ts1.nodes()):
            decoded_metadata = pickle.loads(node.metadata)
            self.assertEqual(decoded_metadata, metadata[j])
        ts1.dump(self.temp_file)
        ts2 = tskit.load(self.temp_file)
        self.assertEqual(ts1.tables.nodes, ts2.tables.nodes)


class ExampleMetadata(object):
    """
    Simple class that we can pickle/unpickle in metadata.
    """
    def __init__(self, one=None, two=None):
        self.one = one
        self.two = two


class TestMetadataPickleDecoding(unittest.TestCase):
    """
    Tests in which use pickle.pickle to decode metadata in nodes, sites and mutations.
    """

    def test_nodes(self):
        tables = tskit.TableCollection(sequence_length=1)
        metadata = ExampleMetadata(one="node1", two="node2")
        pickled = pickle.dumps(metadata)
        tables.nodes.add_row(time=0.125, metadata=pickled)
        ts = tables.tree_sequence()
        node = ts.node(0)
        self.assertEqual(node.time, 0.125)
        self.assertEqual(node.metadata, pickled)
        unpickled = pickle.loads(node.metadata)
        self.assertEqual(unpickled.one, metadata.one)
        self.assertEqual(unpickled.two, metadata.two)

    def test_sites(self):
        tables = tskit.TableCollection(sequence_length=1)
        metadata = ExampleMetadata(one="node1", two="node2")
        pickled = pickle.dumps(metadata)
        tables.sites.add_row(position=0.1, ancestral_state="A", metadata=pickled)
        ts = tables.tree_sequence()
        site = ts.site(0)
        self.assertEqual(site.position, 0.1)
        self.assertEqual(site.ancestral_state, "A")
        self.assertEqual(site.metadata, pickled)
        unpickled = pickle.loads(site.metadata)
        self.assertEqual(unpickled.one, metadata.one)
        self.assertEqual(unpickled.two, metadata.two)

    def test_mutations(self):
        tables = tskit.TableCollection(sequence_length=1)
        metadata = ExampleMetadata(one="node1", two="node2")
        pickled = pickle.dumps(metadata)
        tables.nodes.add_row(time=0)
        tables.sites.add_row(position=0.1, ancestral_state="A")
        tables.mutations.add_row(site=0, node=0, derived_state="T", metadata=pickled)
        ts = tables.tree_sequence()
        mutation = ts.site(0).mutations[0]
        self.assertEqual(mutation.site, 0)
        self.assertEqual(mutation.node, 0)
        self.assertEqual(mutation.derived_state, "T")
        self.assertEqual(mutation.metadata, pickled)
        unpickled = pickle.loads(mutation.metadata)
        self.assertEqual(unpickled.one, metadata.one)
        self.assertEqual(unpickled.two, metadata.two)


class TestJsonSchemaDecoding(unittest.TestCase):
    """
    Tests in which use json-schema to decode the metadata.
    """
    schema = """{
        "title": "Example Metadata",
        "type": "object",
        "properties": {
            "one": {"type": "string"},
            "two": {"type": "string"}
        },
        "required": ["one", "two"]
    }"""

    def test_nodes(self):
        tables = tskit.TableCollection(sequence_length=1)
        builder = pjs.ObjectBuilder(json.loads(self.schema))
        ns = builder.build_classes()
        metadata = ns.ExampleMetadata(one="node1", two="node2")
        encoded = json.dumps(metadata.as_dict()).encode()
        tables.nodes.add_row(time=0.125, metadata=encoded)
        ts = tables.tree_sequence()
        node = ts.node(0)
        self.assertEqual(node.time, 0.125)
        self.assertEqual(node.metadata, encoded)
        decoded = ns.ExampleMetadata.from_json(node.metadata.decode())
        self.assertEqual(decoded.one, metadata.one)
        self.assertEqual(decoded.two, metadata.two)


class TestLoadTextMetadata(unittest.TestCase):
    """
    Tests that use the load_text interface.
    """

    def test_individuals(self):
        individuals = six.StringIO("""\
        id  flags location     metadata
        0   1     0.0,1.0,0.0  abc
        1   1     1.0,2.0      XYZ+
        2   0     2.0,3.0,0.0  !@#$%^&*()
        """)
        i = tskit.parse_individuals(
            individuals, strict=False, encoding='utf8', base64_metadata=False)
        expected = [(1, [0.0, 1.0, 0.0], 'abc'),
                    (1, [1.0, 2.0], 'XYZ+'),
                    (0, [2.0, 3.0, 0.0], '!@#$%^&*()')]
        for a, b in zip(expected, i):
            self.assertEqual(a[0], b.flags)
            self.assertEqual(len(a[1]), len(b.location))
            for x, y in zip(a[1], b.location):
                self.assertEqual(x, y)
            self.assertEqual(a[2].encode('utf8'),
                             b.metadata)

    def test_nodes(self):
        nodes = six.StringIO("""\
        id  is_sample   time    metadata
        0   1           0   abc
        1   1           0   XYZ+
        2   0           1   !@#$%^&*()
        """)
        n = tskit.parse_nodes(
            nodes, strict=False, encoding='utf8', base64_metadata=False)
        expected = ['abc', 'XYZ+', '!@#$%^&*()']
        for a, b in zip(expected, n):
            self.assertEqual(a.encode('utf8'),
                             b.metadata)

    def test_sites(self):
        sites = six.StringIO("""\
        position    ancestral_state metadata
        0.1 A   abc
        0.5 C   XYZ+
        0.8 G   !@#$%^&*()
        """)
        s = tskit.parse_sites(
            sites, strict=False, encoding='utf8', base64_metadata=False)
        expected = ['abc', 'XYZ+', '!@#$%^&*()']
        for a, b in zip(expected, s):
            self.assertEqual(a.encode('utf8'),
                             b.metadata)

    def test_mutations(self):
        mutations = six.StringIO("""\
        site    node    derived_state   metadata
        0   2   C   mno
        0   3   G   )(*&^%$#@!
        """)
        m = tskit.parse_mutations(
            mutations, strict=False, encoding='utf8', base64_metadata=False)
        expected = ['mno', ')(*&^%$#@!']
        for a, b in zip(expected, m):
            self.assertEqual(a.encode('utf8'),
                             b.metadata)

    def test_populations(self):
        populations = six.StringIO("""\
        id    metadata
        0     mno
        1     )(*&^%$#@!
        """)
        p = tskit.parse_populations(
            populations, strict=False, encoding='utf8', base64_metadata=False)
        expected = ['mno', ')(*&^%$#@!']
        for a, b in zip(expected, p):
            self.assertEqual(a.encode('utf8'),
                             b.metadata)
