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


class Column(object):
    def __init__(self, name):
        self.name = name


class Int32Column(Column):
    def get_input(self, n):
        return np.arange(n, dtype=np.int32)


class UInt32Column(Column):
    def get_input(self, n):
        return np.arange(n, dtype=np.uint32)


class CharColumn(Column):
    def get_input(self, n):
        return np.zeros(n, dtype=np.int8)


class DoubleColumn(Column):
    def get_input(self, n):
        return np.arange(n, dtype=np.float64)


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
        kwargs = {c.name: c.get_input(1) for c in self.columns}
        # Make sure this works.
        table = self.table_class()
        table.set_columns(**kwargs)
        for focal_col in self.columns:
            table = self.table_class()
            for bad_type in [Exception, msprime]:
                error_kwargs = dict(kwargs)
                error_kwargs[focal_col.name] = bad_type
                self.assertRaises(TypeError, table.set_columns, **error_kwargs)
            for bad_value in ["qwer", [0, "sd"]]:
                error_kwargs = dict(kwargs)
                error_kwargs[focal_col.name] = bad_value
                self.assertRaises(ValueError, table.set_columns, **error_kwargs)

    def test_set_columns_input_sizes(self):
        num_rows = 100
        input_data = {col.name: col.get_input(num_rows) for col in self.columns}
        col_map = {col.name: col for col in self.columns}
        table = self.table_class()
        table.set_columns(**input_data)
        for equal_len_col_set in self.equal_len_columns:
            if len(equal_len_col_set) > 1:
                for col in equal_len_col_set:
                    kwargs = dict(input_data)
                    kwargs[col] = col_map[col].get_input(1)
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
                setattr(table, col.name, np.zeros(5))
        self.assertEqual(table.num_rows, 0)

    def test_defaults(self):
        table = self.table_class()
        self.assertEqual(table.num_rows, 0)
        for param in self.input_parameters:
            self.assertEqual(getattr(table, param), 1024)
        for col in self.columns:
            array = getattr(table, col.name)
            self.assertEqual(array.shape, (0,))

    def test_set_columns_data(self):
        for num_rows in [0, 10, 100, 1000]:
            input_data = {
                col.name: col.get_input(num_rows) for col in self.columns}
            table = self.table_class()
            table.set_columns(**input_data)
            for colname, input_array in input_data.items():
                output_array = getattr(table, colname)
                self.assertEqual(input_array.shape, output_array.shape)
                self.assertTrue(np.all(input_array == output_array))


class TestNodeTable(unittest.TestCase, CommonTestsMixin):

    columns = [
        UInt32Column("flags"),
        DoubleColumn("time"),
        Int32Column("population"),
        CharColumn("name")]
    input_parameters = ["max_rows_increment"]
    equal_len_columns = [["time", "flags", "population"]]
    table_class = msprime.NodeTable

    def test_variable_stuff(self):
        flags = np.arange(3, dtype=np.uint32)
        time = np.arange(3)
        names = [b"one", b"two", b"three"]
        packed = np.frombuffer(b'\0'.join(names + [b""]), dtype=np.int8)
        table = msprime.NodeTable()
        table.set_columns(flags=flags, time=time, name=packed)
        self.assertTrue(np.all(table.flags == flags))
        self.assertTrue(np.all(table.time == time))
        self.assertTrue(np.all(table.name == packed))
        unpacked = table.name.tostring().split(b"\0")[:-1]
        self.assertEqual(unpacked, names)


class TestEdgesetTable(unittest.TestCase, CommonTestsMixin):

    columns = [
        DoubleColumn("left"),
        DoubleColumn("right"),
        Int32Column("parent"),
        Int32Column("children")]
    equal_len_columns = [["left", "right", "parent"]]
    input_parameters = ["max_rows_increment", "max_children_length_increment"]
    table_class = msprime.EdgesetTable


@unittest.skip("set_column")
class TestMutationTypesTable(unittest.TestCase, CommonTestsMixin):
    columns = [
        CharColumn("ancestral_state"),
        CharColumn("derived_state")]
    equal_len_columns = [["ancestral_state", "derived_state"]]
    input_parameters = ["max_rows_increment"]
    table_class = msprime.MutationTypeTable


class TestMutationsTable(unittest.TestCase, CommonTestsMixin):
    columns = [
        DoubleColumn("position"),
        Int32Column("nodes")]
    equal_len_columns = [["position"]]
    input_parameters = ["max_rows_increment", "max_nodes_length_increment"]
    table_class = msprime.MutationTable


class TestMigrationsTable(unittest.TestCase, CommonTestsMixin):
    columns = [
        DoubleColumn("left"),
        DoubleColumn("right"),
        Int32Column("node"),
        Int32Column("source"),
        Int32Column("dest"),
        DoubleColumn("time")]
    input_parameters = ["max_rows_increment"]
    equal_len_columns = [["left", "right", "node", "source", "dest", "time"]]
    table_class = msprime.MigrationTable
