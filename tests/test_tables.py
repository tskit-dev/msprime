#
# Copyright (C) 2017 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
Test cases for the low-level tables used to transfer information
between simulations and the tree sequence.
"""
from __future__ import print_function
from __future__ import division

import unittest

import numpy as np

import msprime


class TestTable(unittest.TestCase):

    def test_defaults(self):
        table = msprime.NodeTable()
        self.assertEqual(table.max_rows_increment, 1024)
        self.assertEqual(table.num_rows, 0)
        self.assertEqual(table.time.shape, (0,))
        self.assertEqual(table.flags.shape, (0,))

    def test_constructor(self):
        for bad_type in [Exception, msprime]:
            self.assertRaises(TypeError, msprime.NodeTable, flags=[0], time=bad_type)
            self.assertRaises(TypeError, msprime.NodeTable, flags=bad_type, time=[0])
        for bad_value in ["qwer", [0, "sd"]]:
            self.assertRaises(ValueError, msprime.NodeTable, flags=[0], time=bad_value)
            self.assertRaises(ValueError, msprime.NodeTable, flags=bad_value, time=[0])
        # Must specify both time and flags.
        self.assertRaises(TypeError, msprime.NodeTable, time=[1, 2])
        self.assertRaises(TypeError, msprime.NodeTable, flags=[1, 2])
        # Dimensions must be equal
        self.assertRaises(ValueError, msprime.NodeTable, time=[1, 2], flags=[1])
        self.assertRaises(ValueError, msprime.NodeTable, time=[1], flags=[1, 2])

    def test_max_rows_increment(self):
        table = msprime.NodeTable()
        for bad_value in [-1, 0, -2**10]:
            with self.assertRaises(ValueError):
                table.max_rows_increment = bad_value
        for v in [1, 100, 256]:
            table.max_rows_increment = v
            self.assertEqual(table.max_rows_increment, v)

    def test_set_read_only_attributes(self):
        table = msprime.NodeTable()
        with self.assertRaises(AttributeError):
            table.num_rows = 10
        with self.assertRaises(AttributeError):
            table.time = np.zeros(5)
        with self.assertRaises(AttributeError):
            table.flags = np.zeros(5)
        self.assertEqual(table.num_rows, 0)

    def test_times(self):
        time = np.array([0.1, 0.2, 0.3])
        flags = [0, 0, 0]
        table = msprime.NodeTable(flags=flags, time=time)
        self.assertEqual(table.max_rows_increment, time.shape[0])
        self.assertEqual(table.num_rows, time.shape[0])
        stored_time = table.time
        self.assertTrue(np.all(time == table.time))
        # We should have different objects each time we get the array.
        self.assertNotEqual(id(stored_time), id(table.time))

    def test_flags(self):
        time = np.array([0.1, 0.2, 0.3])
        flags = np.array([0, 1, 2], dtype=np.uint32)
        table = msprime.NodeTable(flags=flags, time=time)
        self.assertEqual(table.max_rows_increment, flags.shape[0])
        self.assertEqual(table.num_rows, flags.shape[0])
        stored_flags = table.flags
        self.assertTrue(np.all(flags == table.flags))
        # We should have different objects each flags we get the array.
        self.assertNotEqual(id(stored_flags), id(table.flags))
