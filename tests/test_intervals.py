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
import gzip
import io
import os
import platform

import numpy as np
import pytest

import msprime

IS_WINDOWS = platform.system() == "Windows"


class TestRateMap:
    """
    Test the underlying RateMap class.

    At the moment, much of the RateMap class is tested in test_recombination_map.py,
    as msprime.RecombinationMap contains a RateMap within it. We should probably move
    some of that testing into this class.
    """

    def test_rate_vs_cumulative(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 3])
        cumulative_rates = np.array([0, 1, 4])
        rate_map1 = msprime.RateMap(positions, rate=rates)
        rate_map2 = msprime.RateMap(positions, cumulative_rate=cumulative_rates)
        assert np.allclose(rate_map1.position, rate_map2.position)
        assert np.allclose(rate_map1.rate, rate_map2.rate)
        assert np.allclose(rate_map1.cumulative, rate_map2.cumulative)
        assert np.isclose(rate_map1.mean_rate, rate_map2.mean_rate)
        assert np.isclose(rate_map1.mean_rate, 2)

    def test_no_rate(self):
        positions = np.array([0, 1, 2])
        with pytest.raises(TypeError, match="One of"):
            msprime.RateMap(positions)

    def test_rate_and_cumulative(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 1])
        cumulative_rates = np.array([0, 1, 2])
        with pytest.raises(TypeError, match="both"):
            msprime.RateMap(positions, rates, cumulative_rates)

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
            with pytest.raises(ValueError):
                msprime.RateMap(pos, rate)

    def test_bad_length(self):
        positions = np.array([0, 1, 2])
        rates = np.array([0, 1, 2])
        with pytest.raises(ValueError, match="one less entry"):
            msprime.RateMap(positions, rates)

    def test_bad_first_pos(self):
        positions = np.array([1, 2, 3])
        rates = np.array([1, 1])
        with pytest.raises(ValueError, match="First position"):
            msprime.RateMap(positions, rates)

    def test_bad_first_rate(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, -1])
        with pytest.raises(ValueError, match="non-negative"):
            msprime.RateMap(positions, rates)

    def test_start_position(self):
        positions = np.array([0, 1, 2])
        rates = np.array([0, 1])
        rate_map = msprime.RateMap(positions, rates, start_position=0.5)
        assert np.isclose(rate_map.mean_rate, 1 / (1 + 0.5))

    def test_bad_start_position(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 1])
        with pytest.raises(ValueError, match="start_position"):
            msprime.RateMap(positions, rates, start_position=0.5)

    def test_end_position(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 0])
        rate_map = msprime.RateMap(positions, rates, end_position=1.5)
        assert np.isclose(rate_map.mean_rate, 1 / (1 + 0.5))

    def test_bad_end_position(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 1])
        with pytest.raises(ValueError, match="end_position"):
            msprime.RateMap(positions, rates, end_position=0.5)

    def test_read_only(self):
        positions = np.array([0, 0.25, 0.5, 0.75, 1])
        rates = np.array([0.125, 0.25, 0.5, 0.75])  # 1 shorter than positions
        rate_map = msprime.RateMap(positions, rates)
        assert np.all(rates == rate_map.rate)
        assert np.all(positions == rate_map.position)
        with pytest.raises(AttributeError):
            rate_map.rate = 2 * rate_map.rate
        with pytest.raises(AttributeError):
            rate_map.position = 2 * rate_map.position
        with pytest.raises(AttributeError):
            rate_map.cumulative = 2 * rate_map.cumulative
        with pytest.raises(ValueError):
            rate_map.rate[0] = 1
        with pytest.raises(ValueError):
            rate_map.position[0] = 1
        with pytest.raises(ValueError):
            rate_map.cumulative[0] = 1


class TestSlice:
    def test_slice(self):
        # test RateMap.slice(..., trim=False)
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 3])
        b = a.slice()
        assert a.sequence_length == b.sequence_length
        assert np.array_equal(a.position, b.position)
        assert np.array_equal(a.rate, b.rate)

        b = a.slice(start=50)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 100, 200, 300, 400], b.position)
        assert np.array_equal([0, 1, 2, 3], b.rate)

        b = a.slice(start=100)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 100, 200, 300, 400], b.position)
        assert np.array_equal([0, 1, 2, 3], b.rate)

        b = a.slice(start=150)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 150, 200, 300, 400], b.position)
        assert np.array_equal([0, 1, 2, 3], b.rate)

        b = a.slice(end=300)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 100, 200, 300, 400], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(end=250)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 100, 200, 250, 400], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(start=50, end=300)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 100, 200, 300, 400], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(start=150, end=250)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 150, 200, 250, 400], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(start=150, end=300)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 150, 200, 300, 400], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(start=150, end=160)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([0, 150, 160, 400], b.position)
        assert np.array_equal([0, 1, 0], b.rate)

        # If we take an end-slice into a trailing zero-rate region,
        # we should recover the same map.
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 0])
        b = a.slice(end=350)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal(a.position, b.position)
        assert np.array_equal(a.rate, b.rate)

        b = a.slice(end=300)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal(a.position, b.position)
        assert np.array_equal(a.rate, b.rate)

    def test_slice_with_floats(self):
        #  test RateMap.slice(..., trim=False) with floats
        a = msprime.RateMap([np.pi * x for x in [0, 100, 200, 300, 400]], [0, 1, 2, 3])
        b = a.slice(start=50 * np.pi)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal(a.position, b.position)
        assert np.array_equal(a.rate, b.rate)

        b = a.slice(start=150 * np.pi)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([np.pi * x for x in [0, 150, 200, 300, 400]], b.position)
        assert np.array_equal([0, 1, 2, 3], b.rate)

        b = a.slice(end=300 * np.pi)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([np.pi * x for x in [0, 100, 200, 300, 400]], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(end=250 * np.pi)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([np.pi * x for x in [0, 100, 200, 250, 400]], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(start=50 * np.pi, end=300 * np.pi)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([np.pi * x for x in [0, 100, 200, 300, 400]], b.position)
        assert np.array_equal([0, 1, 2, 0], b.rate)

        b = a.slice(start=150 * np.pi, end=160 * np.pi)
        assert a.sequence_length == b.sequence_length
        assert np.array_equal([np.pi * x for x in [0, 150, 160, 400]], b.position)
        assert np.array_equal([0, 1, 0], b.rate)

    def test_slice_error(self):
        recomb_map = msprime.RateMap([0, 100], [1])
        with pytest.raises(IndexError):
            recomb_map.slice(start=-1)
        with pytest.raises(IndexError):
            recomb_map.slice(end=-1)
        with pytest.raises(IndexError):
            recomb_map.slice(start=200)
        with pytest.raises(IndexError):
            recomb_map.slice(end=200)
        with pytest.raises(IndexError):
            recomb_map.slice(start=20, end=10)

    def test_getitem_slice(self):
        # test RateMap slice syntax
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 3])
        b = a[:]
        assert a.sequence_length == b.sequence_length
        assert np.array_equal(a.position, b.position)
        assert np.array_equal(a.rate, b.rate)

        b = a[50:]
        assert 350 == b.sequence_length
        assert np.array_equal([0, 50, 150, 250, 350], b.position)
        assert np.array_equal([0, 1, 2, 3], b.rate)

        b = a[100:]
        assert 300 == b.sequence_length
        assert np.array_equal([0, 100, 200, 300], b.position)
        assert np.array_equal([1, 2, 3], b.rate)

        b = a[150:]
        assert 250 == b.sequence_length
        assert np.array_equal([0, 50, 150, 250], b.position)
        assert np.array_equal([1, 2, 3], b.rate)

        b = a[:300]
        assert 300 == b.sequence_length
        assert np.array_equal([0, 100, 200, 300], b.position)
        assert np.array_equal([0, 1, 2], b.rate)

        b = a[:250]
        assert 250 == b.sequence_length
        assert np.array_equal([0, 100, 200, 250], b.position)
        assert np.array_equal([0, 1, 2], b.rate)

        b = a[50:300]
        assert 250 == b.sequence_length
        assert np.array_equal([0, 50, 150, 250], b.position)
        assert np.array_equal([0, 1, 2], b.rate)

        b = a[100:300]
        assert 200 == b.sequence_length
        assert np.array_equal([0, 100, 200], b.position)
        assert np.array_equal([1, 2], b.rate)

        b = a[150:250]
        assert 100 == b.sequence_length
        assert np.array_equal([0, 50, 100], b.position)
        assert np.array_equal([1, 2], b.rate)

        b = a[150:160]
        assert 10 == b.sequence_length
        assert np.array_equal([0, 10], b.position)
        assert np.array_equal([1], b.rate)

    def test_getitem_slice_with_negative_indexes_and_floats(self):
        # test RateMap slice syntax with negative indexes and floats
        a = msprime.RateMap([0, 100, 200, 300, 400], [0, 1, 2, 3])

        b = a[150:250]
        c = a[150:-150]
        assert np.array_equal(b.position, c.position)
        assert np.array_equal(b.rate, c.rate)

        b = a[150:250]
        c = a[-250:250]
        assert np.array_equal(b.position, c.position)
        assert np.array_equal(b.rate, c.rate)

        b = a[150:250]
        c = a[-250:-150]
        assert np.array_equal(b.position, c.position)
        assert np.array_equal(b.rate, c.rate)

        b = a[: -np.pi]
        c = a[: 400 - np.pi]
        assert np.array_equal(b.position, c.position)
        assert np.array_equal(b.rate, c.rate)

        b = a[-50 * np.pi : -np.pi]
        c = a[400 - 50 * np.pi : 400 - np.pi]
        assert np.array_equal(b.position, c.position)
        assert np.array_equal(b.rate, c.rate)

    def test_getitem_slice_errors(self):
        recomb_map = msprime.RateMap([0, 100], [1])
        with pytest.raises(TypeError):
            recomb_map["foo"]
        with pytest.raises(TypeError):
            recomb_map[50]
        with pytest.raises(IndexError):
            recomb_map[200:]
        with pytest.raises(IndexError):
            recomb_map[:200]
        with pytest.raises(IndexError):
            recomb_map[20:10]
        with pytest.raises(IndexError):
            recomb_map[-10:-20]
        with pytest.raises(IndexError):
            recomb_map[-101:]
        with pytest.raises(IndexError):
            recomb_map[:-101]


class TestReadHapmap:
    """
    Tests file reading code.
    """

    # Read the (normally unused) 3rd ('recombination rate') column from a file
    @staticmethod
    def get_direct_rates(hapfile):
        positions = []
        rates = []
        next(hapfile)  # skip header
        for line in hapfile:
            _, position, rate, _ = line.split()
            positions.append(float(position))
            rates.append(float(rate))
        return positions, np.array(rates) / 1e8  # convert to morgans / bp

    def test_simple(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 x 0
            chr1 1 x 1 x
            s    2 x 6 x x x"""
        )
        rm = msprime.read_hapmap(hapfile)
        np.testing.assert_array_equal(rm.position, [0, 1, 2])
        assert np.allclose(rm.rate, [1e-2, 5e-2])

    def test_sequence_length(self):
        """
        If sequence_length is specified it adds a 0-rate recombination region at the end
        """

    def test_nonzero_start_pos(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 1 x 0 x
            s    2 x 5 x x x"""
        )
        rm = msprime.read_hapmap(hapfile)
        np.testing.assert_array_equal(rm.position, [0, 1, 2])
        assert np.allclose(rm.rate, [0, 5e-2])

    def test_nonzero_start_cm_err(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 x 5
            s    2 x 1 x x"""
        )
        with pytest.raises(ValueError):
            msprime.read_hapmap(hapfile)

    def test_gzipped(self, tmp_path):
        hapfile = os.path.join(tmp_path, "hapmap.txt.gz")
        with gzip.GzipFile(hapfile, "wb") as gzfile:
            gzfile.write(b"HEADER\n")
            gzfile.write(b"chr1 0 x 0\n")
            gzfile.write(b"chr1 1000 x 0.001\n")
            gzfile.write(b"s    2000 x 0.0065\n")
        rm = msprime.read_hapmap(hapfile)
        np.testing.assert_array_equal(rm.position, [0, 1000, 2000])
        assert np.allclose(rm.rate, [1e-8, 5.5e-8])

    def test_example_files(self, hapmap_strings_fixture):
        # Test recombination_rate against some data from real files.
        for hapstring in hapmap_strings_fixture.values():
            hapfile = io.StringIO(hapstring)
            recomb_map = msprime.read_hapmap(hapfile)
            # read_hapmap() ignores the rate field, so test against the actual
            # precalculated rate field (converted to M/bp) directly.
            hapfile.seek(0)
            _, file_rates = self.get_direct_rates(hapfile)
            assert np.allclose(file_rates[:-1], recomb_map.rate[1:])

    def test_example_mean_rate(self, hapmap_strings_fixture):
        # This will only work on a hapmap file starting at 0 cM
        for name, hapstring in hapmap_strings_fixture.items():
            hapfile = io.StringIO(hapstring)
            recomb_map = msprime.read_hapmap(hapfile)
            hapfile.seek(0)
            pos, rates = self.get_direct_rates(hapfile)
            if "nonzero" in name:
                rates = np.insert(rates, 0, 0)  # insert initial 0
                pos = np.insert(pos, 0, 0)
            assert np.allclose(
                recomb_map.mean_rate, np.average(rates[:-1], weights=np.diff(pos))
            )

    def test_read_old_hapmap(self, hapmap_strings_fixture, tmp_path):
        """
        Check new read_hapmap code will read the write_hapmap output from the old code
        """
        for name, hapstring in hapmap_strings_fixture.items():
            oldmap = msprime.RecombinationMap.read_hapmap(io.StringIO(hapstring))
            newmap = msprime.read_hapmap(io.StringIO(hapstring))
            hapfile = os.path.join(tmp_path, "hapmap.txt")
            oldmap.write_hapmap(hapfile)
            new_via_oldmap = msprime.read_hapmap(hapfile)
            assert len(oldmap) == len(newmap) == len(new_via_oldmap)
            assert np.allclose(newmap.position, oldmap.position)
            assert np.allclose(newmap.position, new_via_oldmap.position)
            assert np.allclose(newmap.rate, oldmap.rate)
            assert np.allclose(newmap.rate, new_via_oldmap.rate)
            if "nonzero" in name:
                expected_cumul = newmap.cumulative - newmap.cumulative[1]
                expected_cumul[0] = 0
                assert np.allclose(expected_cumul, new_via_oldmap.cumulative)
                # Start position 0 in new version, as the first region deemed filled
                assert newmap.start_position == 0
                # We expect to get these wrong, as cumulative appears to start at 0
                assert np.isclose(new_via_oldmap.cumulative[0], 0)
                assert new_via_oldmap.start_position != 0
            else:
                assert np.allclose(newmap.cumulative, new_via_oldmap.cumulative)
                assert np.isclose(newmap.start_position, new_via_oldmap.start_position)
            assert np.isclose(newmap.end_position, oldmap.map.end_position)
            assert np.isclose(newmap.end_position, new_via_oldmap.end_position)
