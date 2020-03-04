#
# Copyright (C) 2018 University of Oxford
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
Tests for the recombination map functionality, mapping continuous physical
coordinates to discrete genetic loci and vice versa.
"""
import unittest
import tempfile
import os
import gzip
import warnings

import numpy as np

import msprime


class TestConstructorAndGetters(unittest.TestCase):

    def verify_warning(self, f):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            f()
            self.assertEqual(len(w), 1)

    def test_warn_on_num_loci_equal_seq_len(self):
        self.verify_warning(
                lambda: msprime.RecombinationMap([0, 100], [0.1, 0], num_loci=100))
        self.verify_warning(
                lambda: msprime.RecombinationMap.uniform_map(100, 0.1, num_loci=100))

    def test_unsupported_methods(self):
        recomb_map = msprime.RecombinationMap([0, 10], [.2, 0])
        self.assertRaises(ValueError, recomb_map.get_num_loci)
        self.assertRaises(ValueError, recomb_map.physical_to_genetic, 8)
        self.assertRaises(ValueError, recomb_map.physical_to_discrete_genetic, 8)
        self.assertRaises(ValueError, recomb_map.genetic_to_physical, 8)
        self.assertRaises(ValueError, recomb_map.get_per_locus_recombination_rate)

    def test_total_recombination_rate(self):
        recomb_map = msprime.RecombinationMap([0, 10], [.1, 0])
        self.assertEqual(recomb_map.get_total_recombination_rate(), 1)


class TestReadHapmap(unittest.TestCase):
    """
    Tests file reading code.
    """
    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(suffix="msp_recomb_map")
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.temp_file)
        except Exception:
            pass

    def test_read_hapmap_simple(self):
        with open(self.temp_file, "w+") as f:
            print("HEADER", file=f)
            print("chr1 0 1", file=f)
            print("chr1 1 5 x", file=f)
            print("s    2 0 x x x", file=f)
        rm = msprime.RecombinationMap.read_hapmap(self.temp_file)
        self.assertEqual(rm.get_positions(), [0, 1, 2])
        self.assertEqual(rm.get_rates(), [1e-8, 5e-8, 0])

    def test_read_hapmap_nonzero_start(self):
        with open(self.temp_file, "w+") as f:
            print("HEADER", file=f)
            print("chr1 1 5 x", file=f)
            print("s    2 0 x x x", file=f)
        rm = msprime.RecombinationMap.read_hapmap(self.temp_file)
        self.assertEqual(rm.get_positions(), [0, 1, 2])
        self.assertEqual(rm.get_rates(), [0, 5e-8, 0])

    def test_read_hapmap_nonzero_end(self):
        with open(self.temp_file, "w+") as f:
            print("HEADER", file=f)
            print("chr1 0 5 x", file=f)
            print("s    2 1 x x x", file=f)
        self.assertRaises(
            ValueError, msprime.RecombinationMap.read_hapmap, self.temp_file)

    def test_read_hapmap_gzipped(self):
        try:
            filename = self.temp_file + ".gz"
            with gzip.open(filename, "w+") as f:
                f.write(b"HEADER\n")
                f.write(b"chr1 0 1\n")
                f.write(b"chr1 1 5.5\n")
                f.write(b"s    2 0\n")
            rm = msprime.RecombinationMap.read_hapmap(filename)
            self.assertEqual(rm.get_positions(), [0, 1, 2])
            self.assertEqual(rm.get_rates(), [1e-8, 5.5e-8, 0])
        finally:
            os.unlink(filename)


class TestSlice(unittest.TestCase):
    def test_slice(self):
        # test RecombinationMap.slice(..., trim=False)
        a = msprime.RecombinationMap([0, 100, 200, 300, 400], [0, 1, 2, 3, 0])
        b = a.slice()
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(a.get_positions(), b.get_positions()))
        self.assertTrue(np.array_equal(a.get_rates(), b.get_rates()))

        b = a.slice(start=50)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 3, 0], b.get_rates()))

        b = a.slice(start=100)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 3, 0], b.get_rates()))

        b = a.slice(start=150)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 150, 200, 300, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 3, 0], b.get_rates()))

        b = a.slice(end=300)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(end=250)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 250, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(start=50, end=300)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(start=150, end=250)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 150, 200, 250, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(start=150, end=300)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 150, 200, 300, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(start=150, end=160)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 150, 160, 400], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 0, 0], b.get_rates()))

    def test_slice_with_floats(self):
        #  test RecombinationMap.slice(..., trim=False) with floats
        a = msprime.RecombinationMap(
                [np.pi*x for x in [0, 100, 200, 300, 400]], [0, 1, 2, 3, 0])
        b = a.slice(start=50*np.pi)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(a.get_positions(), b.get_positions()))
        self.assertTrue(np.array_equal(a.get_rates(), b.get_rates()))

        b = a.slice(start=150*np.pi)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(
            [np.pi*x for x in [0, 150, 200, 300, 400]], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 3, 0], b.get_rates()))

        b = a.slice(end=300*np.pi)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(
            [np.pi*x for x in [0, 100, 200, 300, 400]], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(end=250*np.pi)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(
            [np.pi*x for x in [0, 100, 200, 250, 400]], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(start=50*np.pi, end=300*np.pi)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(
            [np.pi*x for x in [0, 100, 200, 300, 400]], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0, 0], b.get_rates()))

        b = a.slice(start=150*np.pi, end=160*np.pi)
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(
            [np.pi*x for x in [0, 150, 160, 400]], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 0, 0], b.get_rates()))

    def test_slice_error(self):
        recomb_map = msprime.RecombinationMap([0, 100], [1, 0])
        with self.assertRaises(IndexError):
            recomb_map.slice(start=-1)
        with self.assertRaises(IndexError):
            recomb_map.slice(end=-1)
        with self.assertRaises(IndexError):
            recomb_map.slice(start=200)
        with self.assertRaises(IndexError):
            recomb_map.slice(end=200)
        with self.assertRaises(IndexError):
            recomb_map.slice(start=20, end=10)

    def test_getitem_slice(self):
        # test RecombinationMap slice syntax
        a = msprime.RecombinationMap([0, 100, 200, 300, 400], [0, 1, 2, 3, 0])
        b = a[:]
        self.assertEqual(a.get_sequence_length(), b.get_sequence_length())
        self.assertTrue(np.array_equal(a.get_positions(), b.get_positions()))
        self.assertTrue(np.array_equal(a.get_rates(), b.get_rates()))

        b = a[50:]
        self.assertEqual(350, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 50, 150, 250, 350], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 3, 0], b.get_rates()))

        b = a[100:]
        self.assertEqual(300, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 300], b.get_positions()))
        self.assertTrue(np.array_equal([1, 2, 3, 0], b.get_rates()))

        b = a[150:]
        self.assertEqual(250, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 50, 150, 250], b.get_positions()))
        self.assertTrue(np.array_equal([1, 2, 3, 0], b.get_rates()))

        b = a[:300]
        self.assertEqual(300, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 300], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.get_rates()))

        b = a[:250]
        self.assertEqual(250, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200, 250], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.get_rates()))

        b = a[50:300]
        self.assertEqual(250, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 50, 150, 250], b.get_positions()))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.get_rates()))

        b = a[100:300]
        self.assertEqual(200, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 100, 200], b.get_positions()))
        self.assertTrue(np.array_equal([1, 2, 0], b.get_rates()))

        b = a[150:250]
        self.assertEqual(100, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 50, 100], b.get_positions()))
        self.assertTrue(np.array_equal([1, 2, 0], b.get_rates()))

        b = a[150:160]
        self.assertEqual(10, b.get_sequence_length())
        self.assertTrue(np.array_equal([0, 10], b.get_positions()))
        self.assertTrue(np.array_equal([1, 0], b.get_rates()))

    def test_getitem_slice_with_negative_indexes_and_floats(self):
        # test RecombinationMap slice syntax with negative indexes and floats
        a = msprime.RecombinationMap([0, 100, 200, 300, 400], [0, 1, 2, 3, 0])

        b = a[150:250]
        c = a[150:-150]
        self.assertTrue(np.array_equal(b.get_positions(), c.get_positions()))
        self.assertTrue(np.array_equal(b.get_rates(), c.get_rates()))

        b = a[150:250]
        c = a[-250:250]
        self.assertTrue(np.array_equal(b.get_positions(), c.get_positions()))
        self.assertTrue(np.array_equal(b.get_rates(), c.get_rates()))

        b = a[150:250]
        c = a[-250:-150]
        self.assertTrue(np.array_equal(b.get_positions(), c.get_positions()))
        self.assertTrue(np.array_equal(b.get_rates(), c.get_rates()))

        b = a[:-np.pi]
        c = a[:400 - np.pi]
        self.assertTrue(np.array_equal(b.get_positions(), c.get_positions()))
        self.assertTrue(np.array_equal(b.get_rates(), c.get_rates()))

        b = a[-50*np.pi:-np.pi]
        c = a[400 - 50*np.pi:400 - np.pi]
        self.assertTrue(np.array_equal(b.get_positions(), c.get_positions()))
        self.assertTrue(np.array_equal(b.get_rates(), c.get_rates()))

    def test_getitem_slice_errors(self):
        recomb_map = msprime.RecombinationMap([0, 100], [1, 0])
        with self.assertRaises(TypeError):
            recomb_map["foo"]
        with self.assertRaises(TypeError):
            recomb_map[50]
        with self.assertRaises(IndexError):
            recomb_map[200:]
        with self.assertRaises(IndexError):
            recomb_map[:200]
        with self.assertRaises(IndexError):
            recomb_map[20:10]
        with self.assertRaises(IndexError):
            recomb_map[-10:-20]
        with self.assertRaises(IndexError):
            recomb_map[-101:]
        with self.assertRaises(IndexError):
            recomb_map[:-101]
