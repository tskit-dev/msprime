# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 University of Oxford
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

import msprime


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
        encoded, offset = msprime.pack_strings(map(json.dumps, metadata))
        nodes.set_columns(
            flags=nodes.flags, time=nodes.time, population=nodes.population,
            metadata_offset=offset, metadata=encoded)
        self.assertTrue(np.array_equal(nodes.metadata_offset, offset))
        self.assertTrue(np.array_equal(nodes.metadata, encoded))
        ts1 = msprime.load_tables(nodes=nodes, edges=tables.edges)
        for j, node in enumerate(ts1.nodes()):
            decoded_metadata = json.loads(node.metadata.decode())
            self.assertEqual(decoded_metadata, metadata[j])
        ts1.dump(self.temp_file)
        ts2 = msprime.load(self.temp_file)
        self.assertEqual(ts1.tables.nodes, ts2.tables.nodes)

    def test_pickle(self):
        ts = msprime.simulate(10, random_seed=1)
        tables = ts.dump_tables()
        nodes = tables.nodes
        # For each node, we create some Python metadata that can be JSON encoded.
        metadata = [
            {"one": j, "two": 2 * j, "three": list(range(j))} for j in range(len(nodes))]
        encoded, offset = msprime.pack_bytes(map(pickle.dumps, metadata))
        nodes.set_columns(
            flags=nodes.flags, time=nodes.time, population=nodes.population,
            metadata_offset=offset, metadata=encoded)
        self.assertTrue(np.array_equal(nodes.metadata_offset, offset))
        self.assertTrue(np.array_equal(nodes.metadata, encoded))
        ts1 = msprime.load_tables(nodes=nodes, edges=tables.edges)
        for j, node in enumerate(ts1.nodes()):
            decoded_metadata = pickle.loads(node.metadata)
            self.assertEqual(decoded_metadata, metadata[j])
        ts1.dump(self.temp_file)
        ts2 = msprime.load(self.temp_file)
        self.assertEqual(ts1.tables.nodes, ts2.tables.nodes)
