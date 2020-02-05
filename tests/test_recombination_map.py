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
