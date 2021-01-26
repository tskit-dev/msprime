#
# Copyright (C) 2018-2020 University of Oxford
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
Tests for the legacy recombination map functionality.
"""
import io
import os
import random
import tempfile
import unittest
import warnings

import pytest

import msprime


class PythonRecombinationMap:
    """
    A Python implementation of the RecombinationMap interface.

    This uses a simple algorithm used in previous versions of msprime.
    """

    def __init__(self, positions, rates):
        assert len(positions) == len(rates)
        assert len(positions) >= 2
        assert sorted(positions) == positions
        assert positions[0] == 0
        self._positions = positions
        self._sequence_length = positions[-1]
        self._rates = rates

    def get_total_recombination_rate(self):
        """
        Returns the effective recombination rate for this genetic map.
        This is the weighted mean of the rates across all intervals.
        """
        x = self._positions
        effective_rate = 0
        for j in range(len(x) - 1):
            length = x[j + 1] - x[j]
            effective_rate += self._rates[j] * length
        return effective_rate

    def _genetic_to_physical_zero_rate(self, v):
        """
        If we have a zero recombination rate throughout then everything except
        L maps to 0.
        """
        if v > 0:
            return self._sequence_length
        return 0

    def _physical_to_genetic_zero_rate(self, x):
        """
        If we have a zero recombination rate throughout, then we only have
        two possible values. Any value < L maps to 0, as this is the start
        of the interval. If x = L, then we map to L.
        """
        ret = 0
        if x >= self._sequence_length:
            ret = self.get_total_recombination_rate()
        return ret

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
        return s

    def genetic_to_physical(self, v):
        if self.get_total_recombination_rate() == 0:
            return self._genetic_to_physical_zero_rate(v)
        u = v
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

    def verify_coordinate_conversion(self, positions, rates):
        """
        Verifies coordinate conversions by the specified RecombinationMap
        instance.
        """
        L = positions[-1]
        rm = msprime.RecombinationMap(positions, rates)
        other_rm = PythonRecombinationMap(positions, rates)

        assert (
            rm.get_total_recombination_rate() == other_rm.get_total_recombination_rate()
        )
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
            assert 0 <= y <= rm.get_total_recombination_rate()
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
            PythonRecombinationMap(positions, rates),
        ]
        for rm in maps:
            total_recomb = rm.get_total_recombination_rate()
            assert 100 == total_recomb
            # Between 0 and 0.25 and 0.5 and 0.75 we should be able to map 1-1
            # in physical coordinates.
            for x in [0, 0.125, 0.25, 0.50001, 0.66, 0.74999]:
                y = rm.physical_to_genetic(x)
                assert 0 <= y <= total_recomb
                z = rm.genetic_to_physical(y)
                self.assertAlmostEqual(x, z)

    def test_all_zero_rate(self):
        positions = [0, 100]
        rates = [0, 0]
        maps = [
            msprime.RecombinationMap(positions, rates),
            PythonRecombinationMap(positions, rates),
        ]
        for rm in maps:
            for x in [0, 50, 51, 55, 99, 100]:
                assert 0 == rm.physical_to_genetic(x)
                assert 0 == rm.genetic_to_physical(0)
                assert 100 == rm.genetic_to_physical(1)

    def test_zero_rate_start(self):
        positions = [0, 50, 100]
        rates = [0, 1, 0]
        maps = [
            msprime.RecombinationMap(positions, rates),
            PythonRecombinationMap(positions, rates),
        ]
        for rm in maps:
            # Anything <= 50 maps to 0
            for x in [0, 10, 49, 50]:
                assert 0 == rm.physical_to_genetic(x)
            assert 0 == rm.genetic_to_physical(0)
            # values > 50 should map to x - 50
            for x in [51, 55, 99, 100]:
                genetic_x = x - 50
                assert genetic_x == rm.physical_to_genetic(x)
                assert rm.genetic_to_physical(genetic_x) == x

    def test_zero_rate_end(self):
        positions = [0, 50, 100]
        rates = [1, 0, 0]
        maps = [
            msprime.RecombinationMap(positions, rates),
            PythonRecombinationMap(positions, rates),
        ]
        for rm in maps:
            # Anything < 50 maps to x
            for x in [0, 10, 49]:
                assert x == rm.physical_to_genetic(x)
                assert x == rm.genetic_to_physical(x)
            # values >= 50 should map to 50
            for x in [50, 51, 55, 99, 100]:
                assert 50 == rm.physical_to_genetic(x)
            assert 50 == rm.genetic_to_physical(50)

    def test_one_rate(self):
        for rate in [0.1, 1.0, 10]:
            for L in [0.1, 1, 10, 1024, 1e6]:
                positions = [0, L]
                rates = [rate, 0]
                rm = msprime.RecombinationMap(positions, rates)
                assert rate * L == rm.get_total_recombination_rate()
                self.verify_coordinate_conversion(positions, rates)

    def test_simple_map(self):
        positions = [0, 0.25, 0.5, 0.75, 1]
        rates = [0.125, 0.25, 0.5, 0.75, 0]
        self.verify_coordinate_conversion(positions, rates)

    def test_random_map(self):
        for size in [2, 3, 4, 100]:
            positions = [0] + sorted(random.random() for _ in range(size - 2)) + [1]
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
                msprime.RecombinationMap.uniform_map(L, 1),
                PythonRecombinationMap([0, L], [1, 0]),
            ]
            for rm in maps:
                for x in range(L + 1):
                    self.assertAlmostEqual(x, rm.genetic_to_physical(x))


class TestReadHapmap:
    """
    Tests file reading code.
    """

    def test_read_hapmap_deprecation_warning(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 1
            chr1 1 5 x
            chr1 2 0 x x x"""
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rm = msprime.RecombinationMap.read_hapmap(hapfile)
            assert len(w) == 1
        assert rm.get_positions() == [0, 1, 2]
        assert rm.get_rates() == [1e-8, 5e-8, 0]


class TestRecombinationMapInterface(unittest.TestCase):
    """
    Tests for the RecombinationMap interface.
    """

    def test_error_on_num_loci_not_equal_seq_len(self):
        with pytest.raises(ValueError):
            msprime.RecombinationMap.uniform_map(100, 0.1, num_loci=10)

    def test_unsupported_methods(self):
        recomb_map = msprime.RecombinationMap([0, 10], [0.2, 0])
        with pytest.raises(ValueError):
            recomb_map.get_num_loci()
        with pytest.raises(ValueError):
            recomb_map.physical_to_discrete_genetic(8)
        with pytest.raises(ValueError):
            recomb_map.get_per_locus_recombination_rate()

    def test_total_recombination_rate(self):
        recomb_map = msprime.RecombinationMap([0, 10], [0.1, 0])
        assert recomb_map.get_total_recombination_rate() == 1

    def test_basic_properties(self):
        recomb_map = msprime.RecombinationMap.uniform_map(10, 1)
        assert recomb_map.get_sequence_length() == 10
        assert recomb_map.get_length() == 10
        assert recomb_map.get_rates() == [1, 0]
        assert recomb_map.get_positions() == [0, 10]
        assert recomb_map.get_size() == 2

    def test_zero_recombination_map(self):
        # test that beginning and trailing zero recombination regions in the
        # recomb map are included in the sequence
        for n in range(3, 10):
            positions = list(range(n))
            rates = [0.0, 0.2] + [0.0] * (n - 2)
            recomb_map = msprime.RecombinationMap(positions, rates)
            ts = msprime.simulate(10, recombination_map=recomb_map)
            assert ts.sequence_length == n - 1
            assert min(ts.tables.edges.left) == 0.0
            assert max(ts.tables.edges.right) == n - 1.0

    def test_mean_recombination_rate(self):
        # Some quick sanity checks.
        recomb_map = msprime.RecombinationMap([0, 1], [1, 0])
        mean_rr = recomb_map.mean_recombination_rate
        assert mean_rr == 1.0

        recomb_map = msprime.RecombinationMap([0, 1, 2], [1, 0, 0])
        mean_rr = recomb_map.mean_recombination_rate
        assert mean_rr == 0.5

        recomb_map = msprime.RecombinationMap([0, 1, 2], [0, 0, 0])
        mean_rr = recomb_map.mean_recombination_rate
        assert mean_rr == 0.0

        # Test mean_recombination_rate is correct after reading from
        # a hapmap file. read_hapmap() ignores the cM
        # field, so here we test against using the cM field directly.
        def hapmap_rr(hapmap_file):
            first_pos = 0
            with open(hapmap_file) as f:
                next(f)  # skip header
                for line in f:
                    pos, rate, cM = map(float, line.split()[1:4])
                    if cM == 0:
                        first_pos = pos
            return cM / 100 / (pos - first_pos)

        hapmap = """chr pos        rate                    cM
                    1   4283592    3.79115663174456        0
                    1   4361401    0.0664276817058413      0.294986106359414
                    1   7979763   10.9082897515584         0.535345505591925
                    1   8007051    0.0976780648822495      0.833010916332456
                    1   8762788    0.0899929572085616      0.906829844052373
                    1   9477943    0.0864382908650907      0.971188757364862
                    1   9696341    4.76495005895746        0.990066707213216
                    1   9752154    0.0864316558730679      1.25601286485381
                    1   9881751    0.0                     1.26721414815999"""
        with tempfile.TemporaryDirectory() as temp_dir:
            hapfile = os.path.join(temp_dir, "hapmap.txt")
            with open(hapfile, "w") as f:
                f.write(hapmap)
            recomb_map = msprime.RateMap.read_hapmap(f.name)
            mean_rr = recomb_map.mean_rate
            mean_rr2 = hapmap_rr(hapfile)
        self.assertAlmostEqual(mean_rr, mean_rr2, places=15)
