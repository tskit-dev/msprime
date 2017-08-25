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
import _msprime


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
        for bad_value in [-1, -2**10]:
            self.assertRaises(ValueError, self.table_class, max_rows_increment=bad_value)
        for v in [1, 100, 256]:
            table = self.table_class(max_rows_increment=v)
            self.assertEqual(table.max_rows_increment, v)
        # Setting zero or not argument both denote the default.
        table = self.table_class()
        self.assertEqual(table.max_rows_increment, 1024)
        table = self.table_class(max_rows_increment=0)
        self.assertEqual(table.max_rows_increment, 1024)

    def test_input_parameters_errors(self):
        self.assertGreater(len(self.input_parameters), 0)
        for param, _ in self.input_parameters:
            for bad_value in [-1, -2**10]:
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
        table.append_columns(**kwargs)
        for focal_col in self.columns:
            table = self.table_class()
            for bad_type in [Exception, msprime]:
                error_kwargs = dict(kwargs)
                error_kwargs[focal_col.name] = bad_type
                self.assertRaises(TypeError, table.set_columns, **error_kwargs)
                self.assertRaises(TypeError, table.append_columns, **error_kwargs)
            for bad_value in ["qwer", [0, "sd"]]:
                error_kwargs = dict(kwargs)
                error_kwargs[focal_col.name] = bad_value
                self.assertRaises(ValueError, table.set_columns, **error_kwargs)
                self.assertRaises(ValueError, table.append_columns, **error_kwargs)

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
        table.append_columns(**input_data)
        for equal_len_col_set in self.equal_len_columns:
            if len(equal_len_col_set) > 1:
                for col in equal_len_col_set:
                    kwargs = dict(input_data)
                    kwargs[col] = col_map[col].get_input(1)
                    self.assertRaises(ValueError, table.set_columns, **kwargs)
                    self.assertRaises(ValueError, table.append_columns, **kwargs)

    def test_set_read_only_attributes(self):
        table = self.table_class()
        with self.assertRaises(AttributeError):
            table.num_rows = 10
        with self.assertRaises(AttributeError):
            table.max_rows = 10
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

    def test_append_columns_data(self):
        for num_rows in [0, 10, 100, 1000]:
            input_data = {
                col.name: col.get_input(num_rows) for col in self.columns}
            for list_col, length_col in self.ragged_list_columns:
                value = list_col.get_input(num_rows)
                input_data[list_col.name] = value
                input_data[length_col.name] = np.ones(num_rows, dtype=np.uint32)
            table = self.table_class()
            for j in range(1, 10):
                table.append_columns(**input_data)
                for colname, values in input_data.items():
                    output_array = getattr(table, colname)
                    input_array = np.hstack([values for _ in range(j)])
                    self.assertEqual(input_array.shape, output_array.shape)
                    self.assertTrue(np.array_equal(input_array, output_array))
                self.assertEqual(table.num_rows, j * num_rows)

    def test_append_columns_max_rows(self):
        for num_rows in [0, 10, 100, 1000]:
            input_data = {
                col.name: col.get_input(num_rows) for col in self.columns}
            for list_col, length_col in self.ragged_list_columns:
                value = list_col.get_input(num_rows)
                input_data[list_col.name] = value
                input_data[length_col.name] = np.ones(num_rows, dtype=np.uint32)
            for max_rows in [0, 1, 8192]:
                table = self.table_class(max_rows_increment=max_rows)
                for j in range(1, 10):
                    table.append_columns(**input_data)
                    self.assertEqual(table.num_rows, j * num_rows)
                    self.assertGreater(table.max_rows, table.num_rows)
                    if table.num_rows < max_rows:
                        self.assertEqual(table.max_rows, max_rows)

    def test_str(self):
        for num_rows in [0, 10]:
            input_data = {
                col.name: col.get_input(num_rows) for col in self.columns}
            for list_col, length_col in self.ragged_list_columns:
                value = list_col.get_input(num_rows)
                input_data[list_col.name] = value
                input_data[length_col.name] = np.ones(num_rows, dtype=np.uint32)
            table = self.table_class()
            table.set_columns(**input_data)
            s = str(table)
            self.assertEqual(len(s.splitlines()), num_rows + 1)

    def test_copy(self):
        for num_rows in [0, 10]:
            input_data = {
                col.name: col.get_input(num_rows) for col in self.columns}
            for list_col, length_col in self.ragged_list_columns:
                value = list_col.get_input(num_rows)
                input_data[list_col.name] = value
                input_data[length_col.name] = np.ones(num_rows, dtype=np.uint32)
            table = self.table_class()
            table.set_columns(**input_data)
            for _ in range(10):
                copy = table.copy()
                self.assertNotEqual(id(copy), id(table))
                self.assertIsInstance(copy, self.table_class)
                self.assertEqual(copy, table)
                table = copy

    def test_equality(self):
        for num_rows in [1, 10, 100]:
            input_data = {
                col.name: col.get_input(num_rows) for col in self.columns}
            for list_col, length_col in self.ragged_list_columns:
                value = list_col.get_input(num_rows)
                input_data[list_col.name] = value
                input_data[length_col.name] = np.ones(num_rows, dtype=np.uint32)
            t1 = self.table_class()
            t2 = self.table_class()
            self.assertEqual(t1, t1)
            self.assertEqual(t1, t2)
            t1.set_columns(**input_data)
            self.assertEqual(t1, t1)
            self.assertNotEqual(t1, t2)
            self.assertNotEqual(t2, t1)
            t2.set_columns(**input_data)
            self.assertEqual(t1, t2)
            self.assertEqual(t2, t2)
            t2.reset()
            self.assertNotEqual(t1, t2)
            self.assertNotEqual(t2, t1)
            # Check each column in turn to see if we are correctly checking values.
            for col in self.columns:
                col_copy = np.copy(input_data[col.name])
                input_data_copy = dict(input_data)
                input_data_copy[col.name] = col_copy
                t2.set_columns(**input_data_copy)
                self.assertEqual(t1, t2)
                col_copy += 1
                t2.set_columns(**input_data_copy)
                self.assertNotEqual(t1, t2)
                self.assertNotEqual(t2, t1)
            for list_col, length_col in self.ragged_list_columns:
                value = list_col.get_input(num_rows)
                input_data_copy = dict(input_data)
                input_data_copy[list_col.name] = value + 1
                t2.set_columns(**input_data_copy)
                self.assertNotEqual(t1, t2)
                value = list_col.get_input(num_rows + 1)
                input_data_copy = dict(input_data)
                input_data_copy[list_col.name] = value
                input_data_copy[length_col.name] = np.ones(num_rows, dtype=np.uint32)
                input_data_copy[length_col.name][0] = 2
                t2.set_columns(**input_data_copy)
                self.assertNotEqual(t1, t2)
                self.assertNotEqual(t2, t1)
            # Different types should always be unequal.
            self.assertNotEqual(t1, None)
            self.assertNotEqual(t1, [])


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


class TestSortTables(unittest.TestCase):
    """
    Tests for the sort_tables method.
    """
    random_seed = 12345

    def verify_randomise_tables(self, ts):
        tables = ts.dump_tables()
        nodes = tables.nodes
        edgesets = tables.edgesets
        sites = tables.sites
        mutations = tables.mutations
        # TODO deal with migrations.

        # Randomise the tables.
        random.seed(self.random_seed)
        randomised_edgesets = list(ts.edgesets())
        random.shuffle(randomised_edgesets)
        new_edgesets = msprime.EdgesetTable()
        for e in randomised_edgesets:
            randomised_children = list(e.children)
            random.shuffle(randomised_children)
            new_edgesets.add_row(e.left, e.right, e.parent, tuple(randomised_children))
        # Verify that import fails for randomised edgesets
        self.assertRaises(
            _msprime.LibraryError, ts.load_tables, nodes=nodes, edgesets=new_edgesets)

        randomised_sites = list(ts.sites())
        random.shuffle(randomised_sites)
        new_sites = msprime.SiteTable()
        new_mutations = msprime.MutationTable()
        for s in randomised_sites:
            new_sites.add_row(s.position, ancestral_state=s.ancestral_state)
            randomised_mutations = list(s.mutations)
            random.shuffle(randomised_mutations)
            for m in randomised_mutations:
                new_mutations.add_row(
                    site=s.index, node=m.node, derived_state=m.derived_state)
        if ts.num_sites > 1:
            # Verify that import fails for randomised sites
            self.assertRaises(
                _msprime.LibraryError, ts.load_tables, nodes=nodes, edgesets=edgesets,
                sites=new_sites, mutations=new_mutations)

        msprime.sort_tables(
            nodes, new_edgesets, sites=new_sites, mutations=new_mutations)
        # Verify the new and old edgesets are equal.
        self.assertEqual(list(edgesets.left), list(new_edgesets.left))
        self.assertEqual(list(edgesets.right), list(new_edgesets.right))
        self.assertEqual(list(edgesets.parent), list(new_edgesets.parent))
        self.assertEqual(list(edgesets.children), list(new_edgesets.children))
        self.assertEqual(
            list(edgesets.children_length), list(new_edgesets.children_length))
        # sites
        self.assertEqual(list(sites.position), list(new_sites.position))
        self.assertEqual(list(sites.ancestral_state), list(new_sites.ancestral_state))
        self.assertEqual(
            list(sites.ancestral_state_length), list(new_sites.ancestral_state_length))
        # mutations
        self.assertEqual(list(mutations.site), list(new_mutations.site))
        self.assertEqual(
            list(mutations.derived_state), list(new_mutations.derived_state))
        self.assertEqual(
            list(mutations.derived_state_length),
            list(new_mutations.derived_state_length))

        # make sure we can import a tree sequence both with and without the sites.
        ts_new = msprime.load_tables(nodes=nodes, edgesets=new_edgesets)
        self.assertEqual(ts_new.num_edgesets, ts.num_edgesets)
        self.assertEqual(ts_new.num_trees, ts.num_trees)
        ts_new = msprime.load_tables(
            nodes=nodes, edgesets=new_edgesets, sites=new_sites, mutations=new_mutations)
        self.assertEqual(ts_new.num_edgesets, ts.num_edgesets)
        self.assertEqual(ts_new.num_trees, ts.num_trees)
        self.assertEqual(ts_new.num_sites, ts.num_sites)
        self.assertEqual(ts_new.num_mutations, ts.num_mutations)

    def test_single_tree_no_mutations(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        self.verify_randomise_tables(ts)

    def test_many_trees_no_mutations(self):
        ts = msprime.simulate(10, recombination_rate=2, random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        self.verify_randomise_tables(ts)

    def test_single_tree_mutations(self):
        ts = msprime.simulate(10, mutation_rate=2, random_seed=self.random_seed)
        self.assertGreater(ts.num_sites, 2)
        self.verify_randomise_tables(ts)

    def test_many_trees_mutations(self):
        ts = msprime.simulate(
            10, recombination_rate=2, mutation_rate=2, random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        self.assertGreater(ts.num_sites, 2)
        self.verify_randomise_tables(ts)

    def get_nonbinary_example(self, mutation_rate):
        ts = msprime.simulate(
            sample_size=20, recombination_rate=10, random_seed=self.random_seed,
            mutation_rate=mutation_rate, demographic_events=[
                msprime.SimpleBottleneck(time=0.5, proportion=1)])
        # Make sure this really has some non-binary nodes
        found = False
        for r in ts.edgesets():
            if len(r.children) > 2:
                found = True
        self.assertTrue(found)
        return ts

    def test_nonbinary_trees(self):
        ts = self.get_nonbinary_example(mutation_rate=0)
        self.assertGreater(ts.num_trees, 2)
        self.verify_randomise_tables(ts)

    def test_nonbinary_trees_mutations(self):
        ts = self.get_nonbinary_example(mutation_rate=2)
        self.assertGreater(ts.num_trees, 2)
        self.assertGreater(ts.num_sites, 2)
        self.verify_randomise_tables(ts)

    def test_nonbinary_mutations(self):
        # Test the sorting behaviour when we have ragged entries in the ancestral
        # and derived states columns.
        nodes = msprime.NodeTable()
        nodes.add_row(time=0)
        for num_mutations in [1, 10, 50, 100, 8192]:
            edgesets = msprime.EdgesetTable()
            sites = msprime.SiteTable()
            mutations = msprime.MutationTable()
            # Create some awkward length ancestral and derived states.
            random.seed(self.random_seed)
            ancestral_states = []
            derived_states = []
            for j in range(num_mutations):
                s = "".join(random.choice("ACTG") for _ in range(random.randint(0, 8)))
                ancestral_states.append(s.encode())
                sites.add_row(
                    position=num_mutations - j, ancestral_state=ancestral_states[-1])
                s = "".join(random.choice("ACTG") for _ in range(random.randint(0, 8)))
                derived_states.append(s.encode())
                mutations.add_row(site=j, node=0, derived_state=derived_states[-1])
            msprime.sort_tables(nodes, edgesets, sites=sites, mutations=mutations)

            self.assertEqual(len(ancestral_states), sites.num_rows)
            ancestral_states.reverse()
            sorted_ancestral_state = sites.ancestral_state
            length = sites.ancestral_state_length
            offset = 0
            for j in range(sites.num_rows):
                s = sorted_ancestral_state[offset: offset + length[j]].tostring()
                self.assertEqual(s, ancestral_states[j])
                offset += length[j]

            self.assertEqual(len(derived_states), mutations.num_rows)
            derived_states.reverse()
            sorted_derived_state = mutations.derived_state
            length = mutations.derived_state_length
            offset = 0
            for j in range(sites.num_rows):
                s = sorted_derived_state[offset: offset + length[j]].tostring()
                self.assertEqual(s, derived_states[j])
                offset += length[j]

    def test_incompatible_edgesets(self):
        ts1 = msprime.simulate(10, random_seed=self.random_seed)
        ts2 = msprime.simulate(20, random_seed=self.random_seed)
        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()
        # The edgesets in tables2 will refer to nodes that don't exist.
        self.assertRaises(
            IndexError, msprime.sort_tables, tables1.nodes, tables2.edgesets)

    def test_incompatible_sites(self):
        ts1 = msprime.simulate(10, random_seed=self.random_seed)
        ts2 = msprime.simulate(10, mutation_rate=2, random_seed=self.random_seed)
        self.assertGreater(ts2.num_sites, 1)
        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()
        # The mutations in tables2 will refer to sites that don't exist.
        self.assertRaises(
            IndexError, msprime.sort_tables, tables1.nodes, tables1.edgesets,
            sites=tables1.sites, mutations=tables2.mutations)

    def test_incompatible_mutation_nodes(self):
        ts1 = msprime.simulate(2, random_seed=self.random_seed)
        ts2 = msprime.simulate(10, mutation_rate=2, random_seed=self.random_seed)
        self.assertGreater(ts2.num_sites, 1)
        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()
        # The mutations in tables2 will refer to nodes that don't exist.
        self.assertRaises(
            IndexError, msprime.sort_tables, tables1.nodes, tables1.edgesets,
            sites=tables2.sites, mutations=tables2.mutations)

    def test_empty_tables(self):
        nodes = msprime.NodeTable()
        edgesets = msprime.EdgesetTable()
        msprime.sort_tables(nodes, edgesets)
        self.assertEqual(nodes.num_rows, 0)
        self.assertEqual(edgesets.num_rows, 0)
        sites = msprime.SiteTable()
        mutations = msprime.MutationTable()
        msprime.sort_tables(nodes, edgesets, sites=sites, mutations=mutations)
        self.assertEqual(sites.num_rows, 0)
        self.assertEqual(mutations.num_rows, 0)
        migrations = msprime.MigrationTable()
        msprime.sort_tables(
            nodes, edgesets, sites=sites, mutations=mutations, migrations=migrations)
        self.assertEqual(migrations.num_rows, 0)
        msprime.sort_tables(nodes, edgesets, migrations=migrations)
        self.assertEqual(migrations.num_rows, 0)

    def test_sort_interface(self):
        self.assertRaises(TypeError, msprime.sort_tables)
        self.assertRaises(TypeError, msprime.sort_tables, nodes=msprime.NodeTable())
        self.assertRaises(
            TypeError, msprime.sort_tables, edgesets=msprime.EdgesetTable())
        self.assertRaises(
            TypeError, msprime.sort_tables, nodes=msprime.NodeTable(), edgesets=None)
        self.assertRaises(
            TypeError, msprime.sort_tables, nodes=None, edgesets=msprime.EdgesetTable())
        nodes = msprime.NodeTable()
        edgesets = msprime.EdgesetTable()
        # Verify that nodes and edgesets are OK
        msprime.sort_tables(nodes=nodes, edgesets=edgesets)
        for bad_type in [None, "", 1]:
            self.assertRaises(
                TypeError, msprime.sort_tables, nodes=None,
                edgesets=msprime.EdgesetTable(), sites=bad_type)
            self.assertRaises(
                TypeError, msprime.sort_tables, nodes=None,
                edgesets=msprime.EdgesetTable(), mutations=bad_type)
            self.assertRaises(
                TypeError, msprime.sort_tables, nodes=None,
                edgesets=msprime.EdgesetTable(), migrations=bad_type)
        # Must specify sites and mutations together.
        self.assertRaises(
            TypeError, msprime.sort_tables, nodes=nodes, edgesets=edgesets,
            sites=msprime.SiteTable())
        self.assertRaises(
            TypeError, msprime.sort_tables, nodes=nodes, edgesets=edgesets,
            mutations=msprime.MutationTable())
        sites = msprime.SiteTable()
        mutations = msprime.MutationTable()
        # Verify that tables are OK.
        msprime.sort_tables(
            nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)
        for bad_type in [None, "", 1]:
            self.assertRaises(
                TypeError, msprime.sort_tables,
                nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations,
                migrations=bad_type)


class TestSimplifyTables(unittest.TestCase):
    """
    Tests for the simplify_tables function.
    """
    random_seed = 42

    def test_full_samples(self):
        for n in [2, 10, 100, 1000]:
            ts = msprime.simulate(
                n, recombination_rate=1, mutation_rate=1, random_seed=self.random_seed)
            tables = ts.dump_tables()
            nodes_before = tables.nodes.copy()
            edgesets_before = tables.edgesets.copy()
            sites_before = tables.sites.copy()
            mutations_before = tables.mutations.copy()
            msprime.simplify_tables(
                samples=list(ts.samples()),
                nodes=tables.nodes, edgesets=tables.edgesets, sites=tables.sites,
                mutations=tables.mutations)
            self.assertEqual(nodes_before, tables.nodes)
            self.assertEqual(edgesets_before, tables.edgesets)
            self.assertEqual(sites_before, tables.sites)
            self.assertEqual(mutations_before, tables.mutations)

    def test_bad_samples(self):
        n = 10
        ts = msprime.simulate(n, random_seed=self.random_seed)
        tables = ts.dump_tables()
        for bad_node in [-1, n, n + 1, ts.num_nodes - 1, ts.num_nodes, 2**31 - 1]:
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, bad_node], nodes=tables.nodes, edgesets=tables.edgesets)

    def test_bad_edgeset_ordering(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        tables = ts.dump_tables()
        edgesets = tables.edgesets
        # Reversing the edgesets violates the ordering constraints.
        edgesets.set_columns(
            left=edgesets.left[::-1], right=edgesets.right[::-1],
            parent=edgesets.parent[::-1], children=edgesets.children[::-1],
            children_length=edgesets.children_length[::-1])
        self.assertRaises(
            _msprime.LibraryError, msprime.simplify_tables,
            samples=[0, 1], nodes=tables.nodes, edgesets=edgesets)

    def test_bad_edgesets(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        for bad_node in [-1, ts.num_nodes, ts.num_nodes + 1, 2**31 - 1]:
            # Bad parent node
            tables = ts.dump_tables()
            edgesets = tables.edgesets
            parent = edgesets.parent
            parent[0] = bad_node
            edgesets.set_columns(
                left=edgesets.left, right=edgesets.right, parent=parent,
                children=edgesets.children, children_length=edgesets.children_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=edgesets)
            # Bad child node
            tables = ts.dump_tables()
            edgesets = tables.edgesets
            children = edgesets.children
            children[0] = bad_node
            edgesets.set_columns(
                left=edgesets.left, right=edgesets.right, parent=edgesets.parent,
                children=children, children_length=edgesets.children_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=edgesets)
            # child == parent
            tables = ts.dump_tables()
            edgesets = tables.edgesets
            children = edgesets.children
            children[0] = edgesets.parent[0]
            edgesets.set_columns(
                left=edgesets.left, right=edgesets.right, parent=edgesets.parent,
                children=children, children_length=edgesets.children_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=edgesets)
            # left == right
            tables = ts.dump_tables()
            edgesets = tables.edgesets
            left = edgesets.left
            left[0] = edgesets.right[0]
            edgesets.set_columns(
                left=left, right=edgesets.right, parent=edgesets.parent,
                children=edgesets.children, children_length=edgesets.children_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=edgesets)
            # left > right
            tables = ts.dump_tables()
            edgesets = tables.edgesets
            left = edgesets.left
            left[0] = edgesets.right[0] + 1
            edgesets.set_columns(
                left=left, right=edgesets.right, parent=edgesets.parent,
                children=edgesets.children, children_length=edgesets.children_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=edgesets)

    def test_bad_mutation_nodes(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, mutation_rate=1)
        self.assertGreater(ts.num_mutations, 0)
        for bad_node in [-1, ts.num_nodes, 2**31 - 1]:
            tables = ts.dump_tables()
            mutations = tables.mutations
            node = mutations.node
            node[0] = bad_node
            mutations.set_columns(
                site=mutations.site, node=node, derived_state=mutations.derived_state,
                derived_state_length=mutations.derived_state_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=tables.edgesets,
                sites=tables.sites, mutations=mutations)

    def test_bad_mutation_sites(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, mutation_rate=1)
        self.assertGreater(ts.num_mutations, 0)
        for bad_site in [-1, ts.num_sites, 2**31 - 1]:
            tables = ts.dump_tables()
            mutations = tables.mutations
            site = mutations.site
            site[0] = bad_site
            mutations.set_columns(
                site=site, node=mutations.node, derived_state=mutations.derived_state,
                derived_state_length=mutations.derived_state_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=tables.edgesets,
                sites=tables.sites, mutations=mutations)

    def test_bad_site_positions(self):
        ts = msprime.simulate(10, random_seed=self.random_seed, mutation_rate=1)
        self.assertGreater(ts.num_mutations, 0)
        for bad_position in [-1, ts.sequence_length, ts.sequence_length + 0.01]:
            tables = ts.dump_tables()
            sites = tables.sites
            position = sites.position
            position[0] = bad_position
            sites.set_columns(
                position=position, ancestral_state=sites.ancestral_state,
                ancestral_state_length=sites.ancestral_state_length)
            self.assertRaises(
                _msprime.LibraryError, msprime.simplify_tables,
                samples=[0, 1], nodes=tables.nodes, edgesets=tables.edgesets,
                sites=sites, mutations=tables.mutations)

    def test_samples_interface(self):
        tables = msprime.simulate(50, random_seed=1).dump_tables()
        nodes = tables.nodes
        edgesets = tables.edgesets
        for good_form in [[0, 1], (0, 1), np.array([0, 1], dtype=np.int32)]:
            msprime.simplify_tables(good_form, nodes, edgesets)
        for short_samples in [[], (1,), [2], np.array([1], dtype=np.int32)]:
            self.assertRaises(
                ValueError, msprime.simplify_tables, short_samples, nodes, edgesets)
        for bad_type in [None, {}]:
            self.assertRaises(
                TypeError, msprime.simplify_tables, bad_type, nodes, edgesets)
        for bad_numbers in [["s"], ["abc"], [b"1234"]]:
            self.assertRaises(
                ValueError, msprime.simplify_tables, bad_numbers, nodes, edgesets)
        # We only accept numpy arrays of the right type
        for bad_dtype in [np.uint32, np.int64, np.float64]:
            self.assertRaises(
                TypeError, msprime.simplify_tables,
                np.array([0, 1], dtype=bad_dtype), nodes, edgesets)
        bad_samples = np.array([[0, 1], [2, 3]], dtype=np.int32)
        self.assertRaises(
            ValueError, msprime.simplify_tables, bad_samples, nodes, edgesets)

    def test_tables_interface(self):
        self.assertRaises(TypeError, msprime.simplify_tables)
        self.assertRaises(TypeError, msprime.simplify_tables, samples=[0, 1])
        self.assertRaises(
            TypeError, msprime.simplify_tables, samples=[0, 1],
            nodes=msprime.NodeTable())
        self.assertRaises(
            TypeError, msprime.simplify_tables, samples=[0, 1],
            edgesets=msprime.EdgesetTable())
        self.assertRaises(
            TypeError, msprime.simplify_tables, samples=[0, 1],
            nodes=msprime.NodeTable(), edgesets=None)
        self.assertRaises(
            TypeError, msprime.simplify_tables, samples=[0, 1],
            nodes=None, edgesets=msprime.EdgesetTable())
        tables = msprime.simulate(2, random_seed=1).dump_tables()
        nodes = tables.nodes
        edgesets = tables.edgesets
        samples = [0, 1]
        # Verify that samples, nodes and edgesets are OK
        msprime.simplify_tables(samples=samples, nodes=nodes, edgesets=edgesets)
        for bad_type in [None, "", 1]:
            self.assertRaises(
                TypeError, msprime.simplify_tables, samples=samples, nodes=None,
                edgesets=msprime.EdgesetTable(), sites=bad_type)
            self.assertRaises(
                TypeError, msprime.simplify_tables, samples=samples, nodes=None,
                edgesets=msprime.EdgesetTable(), mutations=bad_type)
            self.assertRaises(
                TypeError, msprime.simplify_tables, samples=samples, nodes=None,
                edgesets=msprime.EdgesetTable(), migrations=bad_type)
        # Must specify sites and mutations together.
        self.assertRaises(
            TypeError, msprime.simplify_tables, samples=samples, nodes=nodes,
            edgesets=edgesets, sites=msprime.SiteTable())
        self.assertRaises(
            TypeError, msprime.simplify_tables, samples=samples, nodes=nodes,
            edgesets=edgesets, mutations=msprime.MutationTable())
        sites = msprime.SiteTable()
        mutations = msprime.MutationTable()
        # Verify that tables are OK.
        msprime.simplify_tables(
            samples=samples, nodes=nodes, edgesets=edgesets, sites=sites,
            mutations=mutations)
        for bad_type in [None, "", 1]:
            self.assertRaises(
                TypeError, msprime.simplify_tables,
                nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations,
                migrations=bad_type)
        # Trying to supply migrations fails for now.
        self.assertRaises(
            ValueError, msprime.simplify_tables, samples=samples,
            nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations,
            migrations=msprime.MigrationTable())

    def test_node_table_empty_name_bug(self):
        # Issue #236. Calling simplify on copied tables unexpectedly fails.
        ts = msprime.simulate(20, random_seed=1)
        tables = ts.dump_tables()
        nodes = tables.nodes.copy()
        edgesets = tables.edgesets.copy()
        msprime.simplify_tables(
                samples=ts.samples(), nodes=nodes, edgesets=edgesets)
        self.assertEqual(nodes, tables.nodes)
        self.assertEqual(edgesets, tables.edgesets)
