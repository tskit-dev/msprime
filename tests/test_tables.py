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

import random
import string
import unittest

import numpy as np

import msprime


def random_string(max_length):
    """
    Returns a random character string of the specified maximum length.
    """
    length = random.randint(0, max_length)
    return "".join(random.choice(string.printable) for _ in range(length))


class Column(object):
    def __init__(self, name):
        self.name = name


class Int32Column(Column):
    def get_input(self, n):
        return np.arange(n, dtype=np.int32)


class UInt8Column(Column):
    def get_input(self, n):
        return np.arange(n, dtype=np.uint8)


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
        for param, _ in self.input_parameters:
            for bad_value in [-1, 0, -2**10]:
                self.assertRaises(ValueError, self.table_class, **{param: bad_value})
            for bad_type in [None, ValueError, "ser"]:
                self.assertRaises(TypeError, self.table_class, **{param: bad_type})

    def test_input_parameter_values(self):
        self.assertGreater(len(self.input_parameters), 0)
        for param, _ in self.input_parameters:
            for v in [1, 100, 256]:
                table = self.table_class(**{param: v})
                self.assertEqual(getattr(table, param), v)

    def test_set_columns_string_errors(self):
        inputs = {c.name: c.get_input(1) for c in self.columns}
        for list_col, length_col in self.ragged_list_columns:
            value = list_col.get_input(1)
            inputs[list_col.name] = value
            inputs[length_col.name] = [1]
        # Make sure this works.
        table = self.table_class()
        table.set_columns(**inputs)
        for list_col, length_col in self.ragged_list_columns:
            kwargs = dict(inputs)
            del kwargs[list_col.name]
            self.assertRaises(TypeError, table.set_columns, **kwargs)
            kwargs = dict(inputs)
            del kwargs[length_col.name]
            self.assertRaises(TypeError, table.set_columns, **kwargs)

    def test_set_columns_interface(self):
        kwargs = {c.name: c.get_input(1) for c in self.columns}
        for list_col, length_col in self.ragged_list_columns:
            value = list_col.get_input(1)
            kwargs[list_col.name] = value
            kwargs[length_col.name] = [1]
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
        for list_col, length_col in self.ragged_list_columns:
            value = list_col.get_input(num_rows)
            input_data[list_col.name] = value
            input_data[length_col.name] = np.ones(num_rows, dtype=np.uint32)
            col_map[list_col.name] = list_col
            col_map[length_col.name] = length_col
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
        for param, default in self.input_parameters:
            with self.assertRaises(AttributeError):
                setattr(table, param, 2)
        for col in self.columns:
            with self.assertRaises(AttributeError):
                setattr(table, col.name, np.zeros(5))
        self.assertEqual(table.num_rows, 0)

    def test_defaults(self):
        table = self.table_class()
        self.assertEqual(table.num_rows, 0)
        for param, default in self.input_parameters:
            self.assertEqual(getattr(table, param), default)
        for col in self.columns:
            array = getattr(table, col.name)
            self.assertEqual(array.shape, (0,))

    def test_set_columns_data(self):
        for num_rows in [0, 10, 100, 1000]:
            input_data = {
                col.name: col.get_input(num_rows) for col in self.columns}
            for list_col, length_col in self.ragged_list_columns:
                value = list_col.get_input(num_rows)
                input_data[list_col.name] = value
                input_data[length_col.name] = np.ones(num_rows, dtype=np.uint32)
            table = self.table_class()
            for _ in range(5):
                table.set_columns(**input_data)
                for colname, input_array in input_data.items():
                    output_array = getattr(table, colname)
                    self.assertEqual(input_array.shape, output_array.shape)
                    self.assertTrue(np.all(input_array == output_array))
                table.reset()
                self.assertEqual(table.num_rows, 0)
                for colname in input_data.keys():
                    self.assertEqual(list(getattr(table, colname)), [])


class TestNodeTable(unittest.TestCase, CommonTestsMixin):

    columns = [
        UInt32Column("flags"),
        DoubleColumn("time"),
        Int32Column("population")]
    ragged_list_columns = [(CharColumn("name"),  UInt32Column("name_length"))]
    input_parameters = [("max_rows_increment", 1024)]
    equal_len_columns = [["time", "flags", "population", "name_length"]]
    table_class = msprime.NodeTable

    def test_optional_population(self):
        for num_rows in [0, 10, 100]:
            names = [str(j) for j in range(num_rows)]
            name, name_length = msprime.pack_strings(names)
            flags = list(range(num_rows))
            time = list(range(num_rows))
            table = msprime.NodeTable()
            table.set_columns(
                name=name, name_length=name_length, flags=flags, time=time)
            self.assertEqual(list(table.population), [-1 for _ in range(num_rows)])
            self.assertEqual(list(table.flags), flags)
            self.assertEqual(list(table.time), time)
            self.assertEqual(list(table.name), list(name))
            self.assertEqual(list(table.name_length), list(name_length))

    def test_random_names(self):
        for num_rows in [0, 10, 100]:
            names = [random_string(10) for _ in range(num_rows)]
            name, name_length = msprime.pack_strings(names)
            flags = list(range(num_rows))
            time = list(range(num_rows))
            table = msprime.NodeTable()
            table.set_columns(
                name=name, name_length=name_length, flags=flags, time=time)
            self.assertEqual(list(table.flags), flags)
            self.assertEqual(list(table.time), time)
            self.assertEqual(list(table.name), list(name))
            self.assertEqual(list(table.name_length), list(name_length))
            unpacked_names = msprime.unpack_strings(table.name, table.name_length)
            self.assertEqual(names, unpacked_names)

    def test_optional_names(self):
        for num_rows in [0, 10, 100]:
            flags = list(range(num_rows))
            time = list(range(num_rows))
            table = msprime.NodeTable()
            table.set_columns(flags=flags, time=time)
            self.assertEqual(len(list(table.name)), 0)
            self.assertEqual(list(table.name_length), [0 for _ in range(num_rows)])


class TestEdgesetTable(unittest.TestCase, CommonTestsMixin):

    columns = [
        DoubleColumn("left"),
        DoubleColumn("right"),
        Int32Column("parent")]
    ragged_list_columns = [(Int32Column("children"), UInt32Column("children_length"))]
    equal_len_columns = [["left", "right", "parent", "children_length"]]
    input_parameters = [
        ("max_rows_increment", 1024),
        ("max_children_length_increment", 1024)]
    table_class = msprime.EdgesetTable


class TestSiteTable(unittest.TestCase, CommonTestsMixin):
    columns = [DoubleColumn("position")]
    ragged_list_columns = [
        (CharColumn("ancestral_state"), UInt32Column("ancestral_state_length"))]
    equal_len_columns = [["position", "ancestral_state_length"]]
    input_parameters = [
        ("max_rows_increment", 1024),
        ("max_total_ancestral_state_length_increment", 1024)]
    table_class = msprime.SiteTable


class TestMutationTable(unittest.TestCase, CommonTestsMixin):
    columns = [
        Int32Column("site"),
        Int32Column("node")]
    ragged_list_columns = [
        (CharColumn("derived_state"), UInt32Column("derived_state_length"))]
    equal_len_columns = [["site", "node", "derived_state_length"]]
    input_parameters = [
        ("max_rows_increment", 1024),
        ("max_total_derived_state_length_increment", 1024)]
    table_class = msprime.MutationTable


class TestMigrationTable(unittest.TestCase, CommonTestsMixin):
    columns = [
        DoubleColumn("left"),
        DoubleColumn("right"),
        Int32Column("node"),
        Int32Column("source"),
        Int32Column("dest"),
        DoubleColumn("time")]
    ragged_list_columns = []
    input_parameters = [("max_rows_increment", 1024)]
    equal_len_columns = [["left", "right", "node", "source", "dest", "time"]]
    table_class = msprime.MigrationTable


class TestStringPacking(unittest.TestCase):
    """
    Tests the code for packing and unpacking strings into numpy arrays.
    """

    def test_simple_case(self):
        strings = ["hello", "world"]
        packed, length = msprime.pack_strings(strings)
        self.assertEqual(list(length), [5, 5])
        self.assertEqual(packed.shape, (10,))
        returned = msprime.unpack_strings(packed, length)
        self.assertEqual(returned, strings)

    def verify_packing(self, strings):
        packed, length = msprime.pack_strings(strings)
        self.assertEqual(packed.dtype, np.int8)
        self.assertEqual(length.dtype, np.uint32)
        self.assertEqual(list(length), [len(s) for s in strings])
        self.assertEqual(packed.shape[0], np.sum(length))
        returned = msprime.unpack_strings(packed, length)
        self.assertEqual(strings, returned)

    def test_regular_cases(self):
        for n in range(10):
            strings = ["a" * j for j in range(n)]
            self.verify_packing(strings)

    def test_random_cases(self):
        for n in range(100):
            strings = [random_string(10) for _ in range(n)]
            self.verify_packing(strings)
