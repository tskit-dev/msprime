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
import random
import math

import numpy as np

import msprime


class PythonRecombinationMap(object):
    """
    A Python implementation of the RecombinationMap interface.
    """
    def __init__(self, positions, rates, discrete=False):
        assert len(positions) == len(rates)
        assert len(positions) >= 2
        assert sorted(positions) == positions
        assert positions[0] == 0
        self._positions = positions
        self._sequence_length = positions[-1]
        self._rates = rates
        self._discrete = discrete
        self._cumulative = np.insert(np.cumsum(np.diff(positions) * rates[:-1]), 0, 0)

    def get_total_recombination_rate(self):
        """
        Returns the effective recombination rate for this genetic map.
        This is the weighted mean of the rates across all intervals.
        """
        return self._cumulative[-1]

    def _genetic_to_physical_zero_rate(self, v):
        """
        If we have a zero recombination rate throughout then everything except
        L maps to 0.
        """
        if v >= self.get_total_recombination_rate():
            return self._sequence_length

        return 0

    def physical_to_genetic(self, x):
        return np.interp(x, self._positions, self._cumulative)

    def genetic_to_physical(self, v):
        if self.get_total_recombination_rate() == 0:
            return self._genetic_to_physical_zero_rate(v)
        if v == 0:
            return self._positions[0]

        index = np.searchsorted(self._cumulative, v) - 1
        y = self._positions[index] + (v - self._cumulative[index]) / self._rates[index]
        return math.floor(y) if self._discrete else y


class TestCoordinateConversion(unittest.TestCase):
    """
    Tests that we convert coordinates correctly.
    """

    def verify_coordinate_conversion(self, positions, rates, discrete=False):
        """
        Verifies coordinate conversions by the specified RecombinationMap
        instance.
        """
        L = positions[-1]
        rm = msprime.RecombinationMap(positions, rates, discrete=discrete)
        other_rm = PythonRecombinationMap(positions, rates, discrete=discrete)

        self.assertEqual(
            rm.get_total_recombination_rate(),
            other_rm.get_total_recombination_rate())
        num_random_trials = 10
        num_systematic_trials = 10
        values = [L * random.random() for j in range(num_random_trials)]
        for j in range(num_systematic_trials):
            values.append(L * j / num_systematic_trials)
        values += positions
        for x in values:
            # x is a physical coordinate
            y = rm.physical_to_genetic(x)
            self.assertAlmostEqual(y, other_rm.physical_to_genetic(x), delta=1e-10)
            self.assertTrue(0 <= y <= rm.get_total_recombination_rate())
            # Check if we can round trip approximately in real coordinates.
            xp = rm.genetic_to_physical(y)
            self.assertAlmostEqual(x, xp)
            # The different implementations might differ by very small amounts.
            self.assertAlmostEqual(xp, other_rm.genetic_to_physical(y))

    def test_zero_rate_two_intervals(self):
        # When we have a zero rate in some interval we no longer have a
        # bijective function, since all the physical coordinates in this
        # interval map to a single genetic coordinate.
        positions = [0, 0.25, 0.5, 0.75, 1]
        rates = [200, 0, 200, 0, 0]
        maps = [
            msprime.RecombinationMap(positions, rates),
            PythonRecombinationMap(positions, rates)]
        for rm in maps:
            total_recomb = rm.get_total_recombination_rate()
            self.assertEqual(100, total_recomb)
            # Between 0 and 0.25 and 0.5 and 0.75 we should be able to map 1-1
            # in physical coordinates.
            for x in [0, 0.125, 0.25, 0.50001, 0.66, 0.74999]:
                y = rm.physical_to_genetic(x)
                self.assertTrue(0 <= y <= total_recomb)
                z = rm.genetic_to_physical(y)
                self.assertAlmostEqual(x, z)

    def test_zero_rate_start(self):
        positions = [0, 50, 100]
        rates = [0, 1, 0]
        maps = [
            msprime.RecombinationMap(positions, rates, discrete=True),
            PythonRecombinationMap(positions, rates, discrete=True)]
        for rm in maps:
            # Anything <= 50 maps to 0
            for x in [0, 10, 49, 50]:
                self.assertEqual(0, rm.physical_to_genetic(x))
            self.assertEqual(0, rm.genetic_to_physical(0))
            # values > 50 should map to x - 50
            for x in [51, 55, 99, 100]:
                genetic_x = x - 50
                self.assertEqual(genetic_x, rm.physical_to_genetic(x))
                self.assertEqual(rm.genetic_to_physical(genetic_x), x)

    def test_zero_rate_end(self):
        positions = [0, 50, 100]
        rates = [1, 0, 0]
        maps = [
            msprime.RecombinationMap(positions, rates),
            PythonRecombinationMap(positions, rates)]
        for rm in maps:
            # Anything < 50 maps to x
            for x in [0, 10, 49]:
                self.assertEqual(x, rm.physical_to_genetic(x))
                self.assertEqual(x, rm.genetic_to_physical(x))
            # values >= 50 should map to 50
            for x in [50, 51, 55, 99, 100]:
                self.assertEqual(50, rm.physical_to_genetic(x))
            self.assertEqual(50, rm.genetic_to_physical(50))

    def test_one_rate(self):
        for rate in [0.1, 1.0, 10]:
            for L in [0.1, 1, 10, 1024, 1e6]:
                positions = [0, L]
                rates = [rate, 0]
                rm = msprime.RecombinationMap(positions, rates)
                self.assertEqual(rate * L, rm.get_total_recombination_rate())
                self.verify_coordinate_conversion(positions, rates)

    def test_simple_map(self):
        positions = [0, 0.25, 0.5, 0.75, 1]
        rates = [0.125, 0.25, 0.5, 0.75, 0]
        self.verify_coordinate_conversion(positions, rates)

    def test_random_map(self):
        for size in [2, 3, 4, 100]:
            positions = [0] + sorted(
                random.random() for _ in range(size - 2)) + [1]
            rates = [random.random() for _ in range(size - 1)] + [0]
            self.verify_coordinate_conversion(positions, rates)

    def test_simple_examples(self):
        rm = msprime.RecombinationMap([0, 0.9, 1], [2, 1, 0])
        self.assertAlmostEqual(rm.get_total_recombination_rate(), 1.9)
        rm = msprime.RecombinationMap([0, 0.5, 0.6, 1], [2, 1, 2, 0])
        self.assertAlmostEqual(rm.get_total_recombination_rate(), 1.9)

    def test_integer_round_trip(self):
        for L in [1, 10, 100]:
            maps = [
                msprime.RecombinationMap.uniform_map(L, 1, discrete=True),
                PythonRecombinationMap([0, L], [1, 0], discrete=True)]
            for rm in maps:
                for x in range(L + 1):
                    self.assertAlmostEqual(x, rm.genetic_to_physical(x))


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
        self.assertRaises(ValueError, recomb_map.physical_to_discrete_genetic, 8)
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
