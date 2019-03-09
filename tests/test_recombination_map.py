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
import random
import tempfile
import os
import gzip

import numpy as np

import msprime


class PythonRecombinationMap(object):
    """
    A Python implementation of the RecombinationMap interface.
    """
    def __init__(self, positions, rates, num_loci):
        assert len(positions) == len(rates)
        assert len(positions) >= 2
        assert sorted(positions) == positions
        assert positions[0] == 0
        self._positions = positions
        self._sequence_length = positions[-1]
        self._rates = rates
        self._num_loci = num_loci

    def get_total_recombination_rate(self):
        """
        Returns the effective recombination rate for this genetic map.
        This is the weighted mean of the rates across all intervals.
        """
        x = self._positions
        effective_rate = 0
        for j in range(len(x) - 1):
            length = (x[j + 1] - x[j])
            effective_rate += self._rates[j] * length
        return effective_rate

    def _physical_to_genetic_zero_rate(self, x):
        """
        If we have a zero recombination rate throughout, then we only have
        two possible values. Any value < L maps to 0, as this is the start
        of the interval. If x = L, then we map to num_loci.
        """
        ret = 0
        if x >= self._sequence_length:
            ret = self._num_loci
        return ret

    def _genetic_to_physical_zero_rate(self, v):
        """
        If we have a zero recombination rate throughout then everything except
        L maps to 0.
        """
        return self._sequence_length if v == self._num_loci else 0

    def physical_to_genetic(self, x):
        if self.get_total_recombination_rate() == 0:
            return self._physical_to_genetic_zero_rate(x)
        s = 0
        last_phys_x = 0
        j = 1
        while j < len(self._positions) - 1 and x > self._positions[j]:
            phys_x = self._positions[j]
            rate = self._rates[j - 1]
            s += (phys_x - last_phys_x) * rate
            j += 1
            last_phys_x = phys_x
        rate = self._rates[j - 1]
        s += (x - last_phys_x) * rate
        ret = s / self.get_total_recombination_rate()
        return ret * self._num_loci

    def physical_to_discrete_genetic(self, x):
        return int(round(self.physical_to_genetic(x)))

    def genetic_to_physical(self, v):
        if self.get_total_recombination_rate() == 0:
            return self._genetic_to_physical_zero_rate(v)
        if v == self._num_loci:
            # the code before will not map right ends to right end
            # if the last segments have rate zero.
            y = self._sequence_length
        else:
            # v is expressed in [0, m]. Rescale it back into the range
            # (0, total_mass).
            u = (v / self._num_loci) * self.get_total_recombination_rate()
            s = 0
            last_phys_x = 0
            rate = self._rates[0]
            j = 1
            while j < len(self._positions) and s < u:
                phys_x = self._positions[j]
                rate = self._rates[j - 1]
                s += (phys_x - last_phys_x) * rate
                j += 1
                last_phys_x = phys_x
            y = last_phys_x
            if rate != 0:
                y = last_phys_x - (s - u) / rate
        return y


class TestCoordinateConversion(unittest.TestCase):
    """
    Tests that we convert coordinates correctly.
    """

    def verify_coordinate_conversion(self, positions, rates, num_loci=10):
        """
        Verifies coordinate conversions by the specified RecombinationMap
        instance.
        """
        L = positions[-1]
        rm = msprime.RecombinationMap(positions, rates, num_loci)
        other_rm = PythonRecombinationMap(positions, rates, num_loci)
        if rm.get_size() == 2:
            # When we have very large numbers of loci, this calculations for
            # max distance is off by very small amounts, probably because of
            # machine precision. But if the expected diff is less than say, 10^-10
            # anyway, there's no point in worrying about it.
            max_discretisation_distance = max(1e-10, L / (2 * num_loci))
        else:
            # The above calculation works for a uniform map, but I couldn't
            # figure out how to generalise it. Cop out:
            max_discretisation_distance = L

        self.assertEqual(
            rm.get_total_recombination_rate(),
            other_rm.get_total_recombination_rate())
        num_random_trials = 10
        num_systematic_trials = 10
        values = [L * random.random() for j in range(num_random_trials)]
        for j in range(num_systematic_trials):
            values.append(L * j * 1 / num_systematic_trials)
        values += positions
        for x in values:
            # x is a physical coordinate
            y = rm.physical_to_genetic(x)
            self.assertEqual(y, other_rm.physical_to_genetic(x))
            self.assertTrue(0 <= y <= num_loci)
            # Check if we can round trip approximately in real coordinates.
            xp = rm.genetic_to_physical(y)
            self.assertAlmostEqual(x, xp)
            # The different implementations might differ by very small amounts.
            self.assertAlmostEqual(xp, other_rm.genetic_to_physical(y))

            # Verify the discrete coordinate conversion.
            k = other_rm.physical_to_discrete_genetic(x)
            if y != 0.5:
                # Yuck. Glibc and Python seem to disagree on which way to round
                # when the argument is 1/2. Easiest just skip.
                self.assertEqual(rm.physical_to_discrete_genetic(x), k)
            self.assertTrue(0 <= k <= num_loci)
            x_hat = other_rm.genetic_to_physical(k)
            delta = abs(x - x_hat)
            self.assertGreaterEqual(max_discretisation_distance, delta)

    def test_zero_rate_two_intervals(self):
        # When we have a zero rate in some interval we no longer have a
        # bijective function, since all the physical coordinates in this
        # interval map to a single genetic coordinate.
        positions = [0, 0.25, 0.5, 0.75, 1]
        rates = [1, 0, 1, 0, 0]
        num_loci = 100
        maps = [
            msprime.RecombinationMap(positions, rates, num_loci),
            PythonRecombinationMap(positions, rates, num_loci)]
        for rm in maps:
            self.assertEqual(0.5, rm.get_total_recombination_rate())
            # Between 0 and 0.25 and 0.5 and 0.75 we should be able to map 1-1
            # in physical coordinates.
            for x in [0, 0.125, 0.25, 0.50001, 0.66, 0.74999]:
                y = rm.physical_to_genetic(x)
                self.assertTrue(0 <= y <= num_loci)
                z = rm.genetic_to_physical(y)
                self.assertAlmostEqual(x, z)
            self.assertEqual(0, rm.physical_to_discrete_genetic(0))
            self.assertEqual(25, rm.physical_to_discrete_genetic(0.125))
            self.assertEqual(50, rm.physical_to_discrete_genetic(0.25))
            # Everything withinin the 0.25 to 0.5 interval should map to 50
            self.assertEqual(50, rm.physical_to_discrete_genetic(0.2500001))
            self.assertEqual(50, rm.physical_to_discrete_genetic(0.4))
            self.assertEqual(50, rm.physical_to_discrete_genetic(0.4999999))
            # The discretisation means that we can push values on one side of the
            # interval back to the other side.
            self.assertEqual(50, rm.physical_to_discrete_genetic(0.5000))
            self.assertEqual(51, rm.physical_to_discrete_genetic(0.505))
            # Anything above or equalto 0.75 should map to 100
            self.assertEqual(100, rm.physical_to_discrete_genetic(0.75))
            self.assertEqual(100, rm.physical_to_discrete_genetic(0.751))
            self.assertEqual(100, rm.physical_to_discrete_genetic(0.999))
            self.assertEqual(100, rm.physical_to_discrete_genetic(1.0))

            # All physical coordinates within the first 0 region should map
            # down to the first point, but in the last interval should map
            # to the last point:
            for start, end, point in [(0.25, 0.5, 0.25), (0.75, 1, 1)]:
                for x in [start + delta for delta in [0, 0.01, 0.1]] + [end]:
                    y = rm.physical_to_genetic(x)
                    self.assertTrue(0 <= y <= num_loci)
                    z = rm.genetic_to_physical(y)
                    self.assertEqual(z, point)
                    # We should map exactly in discrete space.
                    k = rm.physical_to_discrete_genetic(x)
                    self.assertEqual(rm.genetic_to_physical(k), point)

    def test_zero_rate_start(self):
        positions = [0, 50, 100]
        rates = [0, 1, 0]
        num_loci = 50
        maps = [
            msprime.RecombinationMap(positions, rates, num_loci),
            PythonRecombinationMap(positions, rates, num_loci)]
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
        num_loci = 50
        maps = [
            msprime.RecombinationMap(positions, rates, num_loci),
            PythonRecombinationMap(positions, rates, num_loci)]
        for rm in maps:
            # Anything < 50 maps to x
            for x in [0, 10, 49]:
                self.assertEqual(x, rm.physical_to_genetic(x))
                self.assertEqual(x, rm.genetic_to_physical(x))
            # values >= 50 should map to 50
            for x in [50, 51, 55, 99, 100]:
                self.assertEqual(50, rm.physical_to_genetic(x))
            self.assertEqual(100, rm.genetic_to_physical(50))

    def test_one_rate(self):
        for num_loci in [1, 10, 1024, 2**31 - 1]:
            for rate in [0.1, 1.0, 10]:
                for L in [0.1, 1, 10, 1024, 1e6]:
                    positions = [0, L]
                    rates = [rate, 0]
                    rm = msprime.RecombinationMap(positions, rates, num_loci)
                    self.assertEqual(rate * L, rm.get_total_recombination_rate())
                    self.verify_coordinate_conversion(positions, rates, num_loci)

    def test_simple_map(self):
        for num_loci in [1, 10, 100, 1025, 2**32 - 1]:
            positions = [0, 0.25, 0.5, 0.75, 1]
            rates = [0.125, 0.25, 0.5, 0.75, 0]
            self.verify_coordinate_conversion(positions, rates, num_loci)

    def test_random_map(self):
        for size in [2, 3, 4, 100]:
            positions = [0] + sorted(
                random.random() for _ in range(size - 2)) + [1]
            rates = [random.random() for _ in range(size - 1)] + [0]
            self.verify_coordinate_conversion(positions, rates)

    def test_simple_examples(self):
        rm = msprime.RecombinationMap([0, 0.9, 1], [2, 1, 0], 10)
        self.assertAlmostEqual(rm.get_total_recombination_rate(), 1.9)
        rm = msprime.RecombinationMap([0, 0.5, 0.6, 1], [2, 1, 2, 0], 100)
        self.assertAlmostEqual(rm.get_total_recombination_rate(), 1.9)

    def test_integer_round_trip(self):
        # We should be able to round trip integer coordinates exactly using
        # a flat recombination map with the right number of loci.
        for L in [1, 10, 100]:
            for rate in [0.1, 1, 100]:
                maps = [
                    msprime.RecombinationMap.uniform_map(L, rate, num_loci=L),
                    PythonRecombinationMap([0, L], [rate, 0], L)]
                for rm in maps:
                    for x in range(L + 1):
                        self.assertEqual(x, rm.physical_to_discrete_genetic(x))
                        self.assertAlmostEqual(x, rm.genetic_to_physical(x))

    def test_single_locus(self):
        eps = 1e-14
        for L in [0.1, 0.99, 1.0, 2, 3.3333, 1e6]:
            for rate in [0.1, 1, 100]:
                maps = [
                    msprime.RecombinationMap.uniform_map(L, rate, num_loci=1),
                    PythonRecombinationMap([0, L], [rate, 0], 1)]
                for rm in maps:
                    self.assertEqual(0, rm.physical_to_discrete_genetic(0))
                    self.assertEqual(0, rm.physical_to_discrete_genetic(eps))
                    self.assertEqual(0, rm.physical_to_discrete_genetic(L / 4))
                    self.assertEqual(1, rm.physical_to_discrete_genetic(L))
                    self.assertEqual(1, rm.physical_to_discrete_genetic(L - eps))
                    self.assertEqual(1, rm.physical_to_discrete_genetic(L - L / 4))

    def test_zero_recombination_rate(self):
        eps = 1e-10
        for L in [0.1, 0.99, 1.0, 2, 3.3333, 1e6]:
            maps = [
                msprime.RecombinationMap.uniform_map(L, 0, num_loci=1),
                PythonRecombinationMap([0, L], [0, 0], 1)]
            for rm in maps:
                self.assertEqual(0, rm.physical_to_discrete_genetic(0))
                self.assertEqual(0, rm.physical_to_discrete_genetic(eps))
                self.assertEqual(0, rm.physical_to_discrete_genetic(L / 4))
                self.assertEqual(1, rm.physical_to_discrete_genetic(L))
                # Even things that are closer to L are mapped down to zero
                # because interval is empty. L only maps to 1 because the interval
                # is half-open and so L is outside of it.
                self.assertEqual(0, rm.physical_to_discrete_genetic(L - L / 4))
                self.assertEqual(0, rm.physical_to_discrete_genetic(L - eps))

    def test_zero_rate_many_loci(self):
        for L in [0.125, 1, 100]:
            positions = [0, L]
            rates = [0, 0]
            for m in [1, 10, 16]:
                maps = [
                    msprime.RecombinationMap(positions, rates, m),
                    PythonRecombinationMap(positions, rates, m)]
                for rm in maps:
                    self.assertEqual(0.0, rm.get_total_recombination_rate())
                    # Any physical value < L should map to 0
                    for x in np.array([0, 0.24, 0.33, 0.99]) * L:
                        self.assertEqual(rm.physical_to_genetic(x), 0)
                        self.assertEqual(rm.physical_to_discrete_genetic(x), 0)
                    self.assertEqual(rm.physical_to_discrete_genetic(L), m)
                    self.assertEqual(rm.physical_to_genetic(L), m)
                    # Any genetic value from 0 to L - 1 should map to 0 in physical
                    # coordinates
                    for y in range(0, m):
                        self.assertEqual(rm.genetic_to_physical(y), 0)
                    self.assertEqual(rm.genetic_to_physical(m), L)


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
