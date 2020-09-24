#
# Copyright (C) 2020 University of Oxford
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
Test cases for the intervals module.
"""
import unittest

import numpy as np

import msprime


class TestRateMap(unittest.TestCase):
    def test_bad_input(self):
        bad_inputs = [
            ([], []),
            ([0], []),
            ([0], [0]),
            ([1, 2], [0]),
            ([0, -1], [0]),
            ([0, 1], [-1]),
        ]

        for pos, rate in bad_inputs:
            with self.assertRaises(ValueError):
                msprime.RateMap(pos, rate)


class TestSlice(unittest.TestCase):
    def test_slice(self):
        # test RateMap.slice(..., trim=False)
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 3])
        b = a.slice()
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal(a.position, b.position))
        self.assertTrue(np.array_equal(a.rate, b.rate))

        b = a.slice(start=50)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 3], b.rate))

        b = a.slice(start=100)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 3], b.rate))

        b = a.slice(start=150)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 150, 200, 300, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 3], b.rate))

        b = a.slice(end=300)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(end=250)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 250, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(start=50, end=300)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 300, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(start=150, end=250)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 150, 200, 250, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(start=150, end=300)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 150, 200, 300, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(start=150, end=160)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal([0, 150, 160, 400], b.position))
        self.assertTrue(np.array_equal([0, 1, 0], b.rate))

        # If we take an end-slice into a trailing zero-rate region,
        # we should recover the same map.
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 0])
        b = a.slice(end=350)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal(a.position, b.position))
        self.assertTrue(np.array_equal(a.rate, b.rate))

        b = a.slice(end=300)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal(a.position, b.position))
        self.assertTrue(np.array_equal(a.rate, b.rate))

    def test_slice_with_floats(self):
        #  test RateMap.slice(..., trim=False) with floats
        a = msprime.RateMap([np.pi * x for x in [0, 100, 200, 300, 400]], [0, 1, 2, 3])
        b = a.slice(start=50 * np.pi)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal(a.position, b.position))
        self.assertTrue(np.array_equal(a.rate, b.rate))

        b = a.slice(start=150 * np.pi)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(
            np.array_equal([np.pi * x for x in [0, 150, 200, 300, 400]], b.position)
        )
        self.assertTrue(np.array_equal([0, 1, 2, 3], b.rate))

        b = a.slice(end=300 * np.pi)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(
            np.array_equal([np.pi * x for x in [0, 100, 200, 300, 400]], b.position)
        )
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(end=250 * np.pi)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(
            np.array_equal([np.pi * x for x in [0, 100, 200, 250, 400]], b.position)
        )
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(start=50 * np.pi, end=300 * np.pi)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(
            np.array_equal([np.pi * x for x in [0, 100, 200, 300, 400]], b.position)
        )
        self.assertTrue(np.array_equal([0, 1, 2, 0], b.rate))

        b = a.slice(start=150 * np.pi, end=160 * np.pi)
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(
            np.array_equal([np.pi * x for x in [0, 150, 160, 400]], b.position)
        )
        self.assertTrue(np.array_equal([0, 1, 0], b.rate))

    def test_slice_error(self):
        recomb_map = msprime.RateMap([0, 100], [1])
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
        # test RateMap slice syntax
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 3])
        b = a[:]
        self.assertEqual(a.sequence_length, b.sequence_length)
        self.assertTrue(np.array_equal(a.position, b.position))
        self.assertTrue(np.array_equal(a.rate, b.rate))

        b = a[50:]
        self.assertEqual(350, b.sequence_length)
        self.assertTrue(np.array_equal([0, 50, 150, 250, 350], b.position))
        self.assertTrue(np.array_equal([0, 1, 2, 3], b.rate))

        b = a[100:]
        self.assertEqual(300, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 300], b.position))
        self.assertTrue(np.array_equal([1, 2, 3], b.rate))

        b = a[150:]
        self.assertEqual(250, b.sequence_length)
        self.assertTrue(np.array_equal([0, 50, 150, 250], b.position))
        self.assertTrue(np.array_equal([1, 2, 3], b.rate))

        b = a[:300]
        self.assertEqual(300, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 300], b.position))
        self.assertTrue(np.array_equal([0, 1, 2], b.rate))

        b = a[:250]
        self.assertEqual(250, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200, 250], b.position))
        self.assertTrue(np.array_equal([0, 1, 2], b.rate))

        b = a[50:300]
        self.assertEqual(250, b.sequence_length)
        self.assertTrue(np.array_equal([0, 50, 150, 250], b.position))
        self.assertTrue(np.array_equal([0, 1, 2], b.rate))

        b = a[100:300]
        self.assertEqual(200, b.sequence_length)
        self.assertTrue(np.array_equal([0, 100, 200], b.position))
        self.assertTrue(np.array_equal([1, 2], b.rate))

        b = a[150:250]
        self.assertEqual(100, b.sequence_length)
        self.assertTrue(np.array_equal([0, 50, 100], b.position))
        self.assertTrue(np.array_equal([1, 2], b.rate))

        b = a[150:160]
        self.assertEqual(10, b.sequence_length)
        self.assertTrue(np.array_equal([0, 10], b.position))
        self.assertTrue(np.array_equal([1], b.rate))

    def test_getitem_slice_with_negative_indexes_and_floats(self):
        # test RateMap slice syntax with negative indexes and floats
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 3])

        b = a[150:250]
        c = a[150:-150]
        self.assertTrue(np.array_equal(b.position, c.position))
        self.assertTrue(np.array_equal(b.rate, c.rate))

        b = a[150:250]
        c = a[-250:250]
        self.assertTrue(np.array_equal(b.position, c.position))
        self.assertTrue(np.array_equal(b.rate, c.rate))

        b = a[150:250]
        c = a[-250:-150]
        self.assertTrue(np.array_equal(b.position, c.position))
        self.assertTrue(np.array_equal(b.rate, c.rate))

        b = a[: -np.pi]
        c = a[: 400 - np.pi]
        self.assertTrue(np.array_equal(b.position, c.position))
        self.assertTrue(np.array_equal(b.rate, c.rate))

        b = a[-50 * np.pi : -np.pi]
        c = a[400 - 50 * np.pi : 400 - np.pi]
        self.assertTrue(np.array_equal(b.position, c.position))
        self.assertTrue(np.array_equal(b.rate, c.rate))

    def test_getitem_slice_errors(self):
        recomb_map = msprime.RateMap([0, 100], [1])
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
