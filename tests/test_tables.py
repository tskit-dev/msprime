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

    columns = None
    table_class = None

    def verify_defaults(self):
        table = self.table_class()
        self.assertEqual(table.max_rows_increment, 1024)
        self.assertEqual(table.num_rows, 0)
        for colname in self.columns:
            array = getattr(table, colname)
            self.assertEqual(array.shape, (0,))

    def verify_max_rows_increment(self):
        for bad_value in [-1, 0, -2**10]:
            self.assertRaises(ValueError, self.table_class, max_rows_increment=bad_value)
        for v in [1, 100, 256]:
            table = self.table_class(max_rows_increment=v)
            self.assertEqual(table.max_rows_increment, v)

    def verify_set_columns_interface(self):
        valid_input = [0]
        kwargs = {c: valid_input for c in self.columns}
        # Make sure this works.
        table = self.table_class()
        table.set_columns(**kwargs)
        for focal_col in self.columns:
            table = self.table_class()
            for bad_type in [Exception, msprime]:
                error_kwargs = dict(kwargs)
                error_kwargs[focal_col] = bad_type
                self.assertRaises(TypeError, table.set_columns, **error_kwargs)
            for bad_value in ["qwer", [0, "sd"]]:
                error_kwargs = dict(kwargs)
                error_kwargs[focal_col] = bad_value
                self.assertRaises(ValueError, table.set_columns, **error_kwargs)

    def verify_set_columns_data(self):
        for num_rows in [0, 10, 100, 1000]:
            input_data = {
                col: np.arange(num_rows, dtype=np.uint32) for col in self.columns}
            table = self.table_class()
            table.set_columns(**input_data)
            for col, input_array in input_data.items():
                output_array = getattr(table, col)
                self.assertEqual(input_array.shape, output_array.shape)
                if not np.all(input_array == output_array):
                    print("not equal:", col)
                    print(input_array)
                    print(output_array)
                # self.assertTrue(np.all(input_array == output_array))

    def verify_set_columns_input_sizes(self, equal_len_cols):
        num_rows = 100
        input_data = {col: np.arange(num_rows, dtype=np.uint32) for col in self.columns}
        table = self.table_class()
        table.set_columns(**input_data)
        for col in equal_len_cols:
            kwargs = dict(input_data)
            kwargs[col] = np.zeros(1, dtype=np.uint32)
            self.assertRaises(ValueError, table.set_columns, **kwargs)

    def verify_constructor_interface(self):
        for bad_type in ["1", None, []]:
            self.assertRaises(TypeError, self.table_class, max_rows_increment=bad_type)

    def verify_set_read_only_attributes(self):
        table = self.table_class()
        with self.assertRaises(AttributeError):
            table.num_rows = 10
        with self.assertRaises(AttributeError):
            table.max_rows_increment = 2
        for col in self.columns:
            with self.assertRaises(AttributeError):
                setattr(table, col, np.zeros(5))
        self.assertEqual(table.num_rows, 0)


class TestNodeTable(TestTable):

    columns = ["flags", "time"]
    table_class = msprime.NodeTable

    def test_defaults(self):
        self.verify_defaults()

    def test_constructor(self):
        self.verify_constructor_interface()

    def test_set_read_only_attributes(self):
        self.verify_set_read_only_attributes()

    def test_set_columns_interface(self):
        self.verify_set_columns_interface()
        table = msprime.NodeTable()
        # Must specify both time and flags.
        self.assertRaises(TypeError, table.set_columns, time=[1, 2])
        self.assertRaises(TypeError, table.set_columns, flags=[1, 2])
        # Dimensions must be equal
        self.assertRaises(ValueError, table.set_columns, time=[1, 2], flags=[1])
        self.assertRaises(ValueError, table.set_columns, time=[1], flags=[1, 2])

    def test_set_columns_data(self):
        self.verify_set_columns_data()

    def test_max_rows_increment(self):
        self.verify_max_rows_increment()

    def test_set_columns_input_sizes(self):
        equal_len_cols = ["time", "flags"]
        self.verify_set_columns_input_sizes(equal_len_cols)


class TestEdgesetTable(TestTable):
    columns = ["left", "right", "parent", "num_children", "children"]
    table_class = msprime.EdgesetTable

    def test_defaults(self):
        self.verify_defaults()
        table = msprime.EdgesetTable()
        self.assertEqual(table.max_total_children_increment, 1024)

    def test_max_rows_increment(self):
        self.verify_max_rows_increment()

    def test_increment_values(self):
        for bad_value in [-1, 0, -2**10]:
            self.assertRaises(ValueError, self.table_class, max_rows_increment=bad_value)
            self.assertRaises(
                ValueError, self.table_class, max_total_children_increment=bad_value)
        for v in [1, 100, 256]:
            table = self.table_class(max_rows_increment=v)
            self.assertEqual(table.max_rows_increment, v)
            table = self.table_class(max_total_children_increment=v)
            self.assertEqual(table.max_total_children_increment, v)

    def test_set_read_only_attributes(self):
        self.verify_set_read_only_attributes()
        table = self.table_class()
        with self.assertRaises(AttributeError):
            table.max_total_children_increment = 1

    def test_constructor(self):
        for bad_type in ["1", None, []]:
            self.assertRaises(
                TypeError, msprime.NodeTable, max_rows_increment=bad_type)
            self.assertRaises(
                TypeError, msprime.NodeTable, max_total_children_increment=bad_type)

    def test_set_columns_interface(self):
        self.verify_set_columns_interface()

    def test_set_columns_data(self):
        self.verify_set_columns_data()

    def test_set_columns_input_sizes(self):
        equal_len_cols = ["left", "right", "parent", "num_children"]
        self.verify_set_columns_input_sizes(equal_len_cols)
