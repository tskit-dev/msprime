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

import numpy as np
import pytest

import msprime


class TestRateMap:
    """
    Test the underlying RateMap class.

    NB: some of the RateMap class is tested in test_recombination_map.py,
    as msprime.RecombinationMap contains a RateMap within it. We should probably move
    some of that testing into this class.
    """

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

    def test_bad_rate(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, -1])
        with pytest.raises(ValueError, match="negative.*1"):
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
            rate_map.cumulative_mass = 2 * rate_map.cumulative_mass
        with pytest.raises(ValueError):
            rate_map.rate[0] = 1
        with pytest.raises(ValueError):
            rate_map.position[0] = 1
        with pytest.raises(ValueError):
            rate_map.cumulative_mass[0] = 1
        with pytest.raises(NotImplementedError, match="slice"):
            rate_map.start_position = 1
        with pytest.raises(NotImplementedError, match="slice"):
            rate_map.end_position = 1


class TestGetIntermediates:
    def test_get_rate(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(positions, rates)
        assert np.all(rate_map.get_rate([0.5, 1.5]) == rates)

    def test_get_bad_rate(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(positions, rates)
        with pytest.raises(ValueError, match="No rate exists"):
            rate_map.get_rate([1, -0.1])

    def test_get_cumulative_mass(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(positions, rates)
        assert np.allclose(rate_map.get_cumulative_mass([0.5, 1.5]), np.array([0.5, 3]))
        assert rate_map.get_cumulative_mass([2]) == rate_map.total_mass

    def test_get_bad_cumulative_mass(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(positions, rates)
        with pytest.raises(ValueError, match="physical positions"):
            rate_map.get_cumulative_mass([1, -0.1])
        with pytest.raises(ValueError, match="physical positions"):
            rate_map.get_cumulative_mass([1, 2.1])


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
    def test_read_hapmap_simple(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 1 x 0
            chr1 2 x 0.000001 x
            chr1 3 x 0.000006 x x"""
        )
        rm = msprime.RateMap.read_hapmap(hapfile)
        np.testing.assert_array_equal(rm.position, [0, 1, 2, 3])
        assert np.allclose(rm.rate, [0, 1e-8, 5e-8])
        assert rm.start_position == 1

    def test_read_hapmap_from_filename(self, tmp_path):
        with open(tmp_path / "hapfile.txt", "wt") as hapfile:
            hapfile.write(
                """\
                HEADER
                chr1 1 x 0
                chr1 2 x 0.000001 x
                chr1 3 x 0.000006 x x"""
            )
        rm = msprime.RateMap.read_hapmap(tmp_path / "hapfile.txt")
        np.testing.assert_array_equal(rm.position, [0, 1, 2, 3])
        assert np.allclose(rm.rate, [0, 1e-8, 5e-8])
        assert rm.start_position == 1

    @pytest.mark.filterwarnings("ignore:loadtxt")
    def test_read_hapmap_empty(self):
        hapfile = io.StringIO(
            """\
            HEADER"""
        )
        with pytest.raises(ValueError, match="Empty"):
            msprime.RateMap.read_hapmap(hapfile)

    def test_read_hapmap_col_pos(self):
        hapfile = io.StringIO(
            """\
            HEADER
            0 0
            0.000001 1 x
            0.000006 2 x x"""
        )
        rm = msprime.RateMap.read_hapmap(hapfile, position_col=1, map_col=0)
        np.testing.assert_array_equal(rm.position, [0, 1, 2])
        assert np.allclose(rm.rate, [1e-8, 5e-8])

    def test_read_hapmap_map_and_rate(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 0 0
            chr1 1 1 0.000001 x
            chr1 2 2 0.000006 x x"""
        )
        with pytest.raises(ValueError, match="both rate_col and map_col"):
            msprime.RateMap.read_hapmap(hapfile, rate_col=2, map_col=3)

    def test_read_hapmap_duplicate_pos(self):
        hapfile = io.StringIO(
            """\
            HEADER
            0 0
            0.000001 1 x
            0.000006 2 x x"""
        )
        with pytest.raises(ValueError, match="same columns"):
            msprime.RateMap.read_hapmap(hapfile, map_col=1)

    def test_read_hapmap_nonzero_rate_start(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 1 5 x
            chr1 2 0 x x x"""
        )
        rm = msprime.RateMap.read_hapmap(hapfile, rate_col=2)
        np.testing.assert_array_equal(rm.position, [0, 1, 2])
        np.testing.assert_array_equal(rm.rate, [0, 5e-8])

    def test_read_hapmap_nonzero_rate_end(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 5 x
            chr1 2 1 x x x"""
        )
        with pytest.raises(ValueError, match="last entry.*must be zero"):
            msprime.RateMap.read_hapmap(hapfile, rate_col=2)

    def test_read_hapmap_gzipped(self, tmp_path):
        hapfile = os.path.join(tmp_path, "hapmap.txt.gz")
        with gzip.GzipFile(hapfile, "wb") as gzfile:
            gzfile.write(b"HEADER\n")
            gzfile.write(b"chr1 0 1\n")
            gzfile.write(b"chr1 1 5.5\n")
            gzfile.write(b"chr1 2 0\n")
        rm = msprime.RateMap.read_hapmap(hapfile, rate_col=2)
        np.testing.assert_array_equal(rm.position, [0, 1, 2])
        np.testing.assert_array_equal(rm.rate, [1e-8, 5.5e-8])

    def test_read_hapmap_nonzero_map_start(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 1 x 0.000001
            chr1 2 x 0.000001 x
            chr1 3 x 0.000006 x x x"""
        )
        rm = msprime.RateMap.read_hapmap(hapfile)
        np.testing.assert_array_equal(rm.position, [0, 1, 2, 3])
        assert np.allclose(rm.rate, [1e-8, 0, 5e-8])

    def test_read_hapmap_bad_nonzero_map_start(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 x 0.0000005
            chr1 1 x 0.000001 x
            chr1 2 x 0.000006 x x x"""
        )
        with pytest.raises(ValueError, match="start.*must be zero"):
            msprime.RateMap.read_hapmap(hapfile)

    def test_sequence_length(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 x 0
            chr1 1 x 0.000001 x
            chr1 2 x 0.000006 x x x"""
        )
        # test identical seq len
        rm = msprime.RateMap.read_hapmap(hapfile, sequence_length=2)
        np.testing.assert_array_equal(rm.position, [0, 1, 2])
        assert np.allclose(rm.rate, [1e-8, 5e-8])
        assert rm.end_position == 2

        hapfile.seek(0)
        rm = msprime.RateMap.read_hapmap(hapfile, sequence_length=10)
        np.testing.assert_array_equal(rm.position, [0, 1, 2, 10])
        assert np.allclose(rm.rate, [1e-8, 5e-8, 0])
        assert rm.end_position == 2

    def test_bad_sequence_length(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 0 x 0
            chr1 1 x 0.000001 x
            chr1 2 x 0.000006 x x x"""
        )
        with pytest.raises(ValueError, match="sequence_length"):
            msprime.RateMap.read_hapmap(hapfile, sequence_length=1.999)

    def test_no_header(self):
        data = """\
            chr1 0 x 0
            chr1 1 x 0.000001 x
            chr1 2 x 0.000006 x x x"""
        hapfile_noheader = io.StringIO(data)
        hapfile_header = io.StringIO("chr pos rate cM\n" + data)
        with pytest.raises(ValueError):
            msprime.RateMap.read_hapmap(hapfile_header, has_header=False)
        rm1 = msprime.RateMap.read_hapmap(hapfile_header)
        rm2 = msprime.RateMap.read_hapmap(hapfile_noheader, has_header=False)
        assert np.array_equal(rm1.rate, rm2.rate)
        assert np.array_equal(rm1.position, rm2.position)

    def test_hapmap_fragment(self):
        hapfile = io.StringIO(
            """\
            chr pos        rate                    cM
            1   4283592    3.79115663174456        0
            1   4361401    0.0664276817058413      0.294986106359414
            1   7979763   10.9082897515584         0.535345505591925
            1   8007051    0.0976780648822495      0.833010916332456
            1   8762788    0.0899929572085616      0.906829844052373
            1   9477943    0.0864382908650907      0.971188757364862
            1   9696341    4.76495005895746        0.990066707213216
            1   9752154    0.0864316558730679      1.25601286485381
            1   9881751    0.0                     1.26721414815999"""
        )
        rm1 = msprime.RateMap.read_hapmap(hapfile)
        hapfile.seek(0)
        rm2 = msprime.RateMap.read_hapmap(hapfile, rate_col=2)
        assert np.allclose(rm1.rate, rm2.rate)
