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


class CommonTestsMixin(object):
    """
    Abstract base class for common table tests. Because of the design of unittest,
    we have to make this a mixin.
    """

    def test_max_rows_increment(self):
        for bad_value in [-1, 0, -2**10]:
            self.assertRaises(ValueError, self.table_class, max_rows_increment=bad_value)
        for v in [1, 100, 256]:
            table = self.table_class(max_rows_increment=v)
            self.assertEqual(table.max_rows_increment, v)

    def test_input_parameters_errors(self):
        self.assertGreater(len(self.input_parameters), 0)
        for param in self.input_parameters:
            for bad_value in [-1, 0, -2**10]:
                self.assertRaises(ValueError, self.table_class, **{param: bad_value})
            for bad_type in [None, ValueError, "ser"]:
                self.assertRaises(TypeError, self.table_class, **{param: bad_type})

    def test_input_parameter_values(self):
        self.assertGreater(len(self.input_parameters), 0)
        for param in self.input_parameters:
            for v in [1, 100, 256]:
                table = self.table_class(**{param: v})
                self.assertEqual(getattr(table, param), v)

    def test_set_columns_interface(self):
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

    def test_set_columns_input_sizes(self):
        num_rows = 100
        input_data = {col: np.arange(num_rows, dtype=np.uint32) for col in self.columns}
        table = self.table_class()
        table.set_columns(**input_data)
        for equal_len_col_set in self.equal_len_columns:
            for col in equal_len_col_set:
                kwargs = dict(input_data)
                kwargs[col] = np.zeros(1, dtype=np.uint32)
                self.assertRaises(ValueError, table.set_columns, **kwargs)

    def test_set_read_only_attributes(self):
        table = self.table_class()
        with self.assertRaises(AttributeError):
            table.num_rows = 10
        for param in self.input_parameters:
            with self.assertRaises(AttributeError):
                setattr(table, param, 2)
        for col in self.columns:
            with self.assertRaises(AttributeError):
                setattr(table, col, np.zeros(5))
        self.assertEqual(table.num_rows, 0)

    def test_defaults(self):
        table = self.table_class()
        self.assertEqual(table.num_rows, 0)
        for param in self.input_parameters:
            self.assertEqual(getattr(table, param), 1024)
        for colname in self.columns:
            array = getattr(table, colname)
            self.assertEqual(array.shape, (0,))

    def test_set_columns_data(self):
        for num_rows in [0, 10, 100, 1000]:
            input_data = {
                col: np.arange(num_rows, dtype=np.uint32) for col in self.columns}
            table = self.table_class()
            table.set_columns(**input_data)
            for col, input_array in input_data.items():
                output_array = getattr(table, col)
                self.assertEqual(input_array.shape, output_array.shape)
                self.assertTrue(np.all(input_array == output_array))


class TestNodeTable(unittest.TestCase, CommonTestsMixin):

    columns = ["flags", "time"]
    input_parameters = ["max_rows_increment"]
    equal_len_columns = [["time", "flags"]]
    table_class = msprime.NodeTable


class TestEdgesetTable(unittest.TestCase, CommonTestsMixin):

    columns = ["left", "right", "parent", "num_children", "children"]
    equal_len_columns = [["left", "right", "parent", "num_children"]]
    input_parameters = ["max_rows_increment", "max_total_children_increment"]
    table_class = msprime.EdgesetTable


class TestMutationsTable(unittest.TestCase, CommonTestsMixin):
    columns = ["position", "num_nodes", "nodes"]
    equal_len_columns = [["position", "num_nodes"]]
    input_parameters = ["max_rows_increment", "max_total_nodes_increment"]
    table_class = msprime.MutationTable


# class TestMigrationsTable(CommonTestsMixin):
#     columns = []
#     table_class = msprime.MigrationTable

#     def test_defaults(self):
#         self.verify_defaults()
