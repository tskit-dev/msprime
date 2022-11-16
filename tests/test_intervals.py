#
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
import decimal
import fractions
import gzip
import io
import os
import pickle
import textwrap
import xml

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import msprime


class TestRateMapErrors:
    @pytest.mark.parametrize(
        ("position", "rate"),
        [
            ([], []),
            ([0], []),
            ([0], [0]),
            ([1, 2], [0]),
            ([0, -1], [0]),
            ([0, 1], [-1]),
        ],
    )
    def test_bad_input(self, position, rate):
        with pytest.raises(ValueError):
            msprime.RateMap(position=position, rate=rate)

    def test_zero_length_interval(self):
        with pytest.raises(ValueError, match=r"at indexes \[2 4\]"):
            msprime.RateMap(position=[0, 1, 1, 2, 2, 3], rate=[0, 0, 0, 0, 0])

    def test_bad_length(self):
        positions = np.array([0, 1, 2])
        rates = np.array([0, 1, 2])
        with pytest.raises(ValueError, match="one less entry"):
            msprime.RateMap(position=positions, rate=rates)

    def test_bad_first_pos(self):
        positions = np.array([1, 2, 3])
        rates = np.array([1, 1])
        with pytest.raises(ValueError, match="First position"):
            msprime.RateMap(position=positions, rate=rates)

    def test_bad_rate(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, -1])
        with pytest.raises(ValueError, match="negative.*1"):
            msprime.RateMap(position=positions, rate=rates)

    def test_bad_rate_with_missing(self):
        positions = np.array([0, 1, 2])
        rates = np.array([np.nan, -1])
        with pytest.raises(ValueError, match="negative.*1"):
            msprime.RateMap(position=positions, rate=rates)

    def test_read_only(self):
        positions = np.array([0, 0.25, 0.5, 0.75, 1])
        rates = np.array([0.125, 0.25, 0.5, 0.75])  # 1 shorter than positions
        rate_map = msprime.RateMap(position=positions, rate=rates)
        assert np.all(rates == rate_map.rate)
        assert np.all(positions == rate_map.position)
        with pytest.raises(AttributeError):
            rate_map.rate = 2 * rate_map.rate
        with pytest.raises(AttributeError):
            rate_map.position = 2 * rate_map.position
        with pytest.raises(AttributeError):
            rate_map.left = 1234
        with pytest.raises(AttributeError):
            rate_map.right = 1234
        with pytest.raises(AttributeError):
            rate_map.mid = 1234
        with pytest.raises(ValueError):
            rate_map.rate[0] = 1
        with pytest.raises(ValueError):
            rate_map.position[0] = 1
        with pytest.raises(ValueError):
            rate_map.left[0] = 1
        with pytest.raises(ValueError):
            rate_map.mid[0] = 1
        with pytest.raises(ValueError):
            rate_map.right[0] = 1


class TestGetRateAllKnown:
    examples = [
        msprime.RateMap(position=[0, 1], rate=[0]),
        msprime.RateMap(position=[0, 1], rate=[0.1]),
        msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2]),
        msprime.RateMap(position=[0, 1, 2], rate=[0, 0.2]),
        msprime.RateMap(position=[0, 1, 2], rate=[0.1, 1e-6]),
        msprime.RateMap(position=range(100), rate=range(99)),
    ]

    @pytest.mark.parametrize("rate_map", examples)
    def test_get_rate_mid(self, rate_map):
        rate = rate_map.get_rate(rate_map.mid)
        assert len(rate) == len(rate_map)
        for j in range(len(rate_map)):
            assert rate[j] == rate_map[rate_map.mid[j]]

    @pytest.mark.parametrize("rate_map", examples)
    def test_get_rate_left(self, rate_map):
        rate = rate_map.get_rate(rate_map.left)
        assert len(rate) == len(rate_map)
        for j in range(len(rate_map)):
            assert rate[j] == rate_map[rate_map.left[j]]

    @pytest.mark.parametrize("rate_map", examples)
    def test_get_rate_right(self, rate_map):
        rate = rate_map.get_rate(rate_map.right[:-1])
        assert len(rate) == len(rate_map) - 1
        for j in range(len(rate_map) - 1):
            assert rate[j] == rate_map[rate_map.right[j]]


class TestOperations:
    examples = [
        msprime.RateMap(position=[0, 1], rate=[0]),
        msprime.RateMap(position=[0, 1], rate=[0.1]),
        msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2]),
        msprime.RateMap(position=[0, 1, 2], rate=[0, 0.2]),
        msprime.RateMap(position=[0, 1, 2], rate=[0.1, 1e-6]),
        msprime.RateMap(position=range(100), rate=range(99)),
        # Missing data
        msprime.RateMap(position=[0, 1, 2], rate=[np.nan, 0]),
        msprime.RateMap(position=[0, 1, 2], rate=[0, np.nan]),
        msprime.RateMap(position=[0, 1, 2, 3], rate=[0, np.nan, 1]),
    ]

    @pytest.mark.parametrize("rate_map", examples)
    def test_num_intervals(self, rate_map):
        assert rate_map.num_intervals == len(rate_map.rate)
        assert rate_map.num_missing_intervals == np.sum(np.isnan(rate_map.rate))
        assert rate_map.num_non_missing_intervals == np.sum(~np.isnan(rate_map.rate))

    @pytest.mark.parametrize("rate_map", examples)
    def test_mask_arrays(self, rate_map):
        assert_array_equal(rate_map.missing, np.isnan(rate_map.rate))
        assert_array_equal(rate_map.non_missing, ~np.isnan(rate_map.rate))

    @pytest.mark.parametrize("rate_map", examples)
    def test_missing_intervals(self, rate_map):
        missing = []
        for left, right, rate in zip(rate_map.left, rate_map.right, rate_map.rate):
            if np.isnan(rate):
                missing.append([left, right])
        if len(missing) == 0:
            assert len(rate_map.missing_intervals()) == 0
        else:
            assert_array_equal(missing, rate_map.missing_intervals())

    @pytest.mark.parametrize("rate_map", examples)
    def test_mean_rate(self, rate_map):
        total_span = 0
        total_mass = 0
        for span, mass in zip(rate_map.span, rate_map.mass):
            if not np.isnan(mass):
                total_span += span
                total_mass += mass
        assert total_mass / total_span == rate_map.mean_rate

    @pytest.mark.parametrize("rate_map", examples)
    def test_total_mass(self, rate_map):
        assert rate_map.total_mass == np.nansum(rate_map.mass)

    @pytest.mark.parametrize("rate_map", examples)
    def test_get_cumulative_mass(self, rate_map):
        assert list(rate_map.get_cumulative_mass([0])) == [0]
        assert list(rate_map.get_cumulative_mass([rate_map.sequence_length])) == [
            rate_map.total_mass
        ]
        assert_array_equal(
            rate_map.get_cumulative_mass(rate_map.right), np.nancumsum(rate_map.mass)
        )

    @pytest.mark.parametrize("rate_map", examples)
    def test_get_rate(self, rate_map):
        assert_array_equal(rate_map.get_rate([0]), rate_map.rate[0])
        assert_array_equal(
            rate_map.get_rate([rate_map.sequence_length - 1e-9]), rate_map.rate[-1]
        )
        assert_array_equal(rate_map.get_rate(rate_map.left), rate_map.rate)

    @pytest.mark.parametrize("rate_map", examples)
    def test_map_semantics(self, rate_map):
        assert len(rate_map) == rate_map.num_non_missing_intervals
        assert_array_equal(list(rate_map.keys()), rate_map.mid[rate_map.non_missing])
        for x in rate_map.left[rate_map.missing]:
            assert x not in rate_map
        for x in rate_map.mid[rate_map.missing]:
            assert x not in rate_map


class TestFindIndex:
    def test_one_interval(self):
        rate_map = msprime.RateMap(position=[0, 10], rate=[0.1])
        for j in range(10):
            assert rate_map.find_index(j) == 0
        assert rate_map.find_index(0.0001) == 0
        assert rate_map.find_index(9.999) == 0

    def test_two_intervals(self):
        rate_map = msprime.RateMap(position=[0, 5, 10], rate=[0.1, 0.1])
        assert rate_map.find_index(0) == 0
        assert rate_map.find_index(0.0001) == 0
        assert rate_map.find_index(4.9999) == 0
        assert rate_map.find_index(5) == 1
        assert rate_map.find_index(5.1) == 1
        assert rate_map.find_index(7) == 1
        assert rate_map.find_index(9.999) == 1

    def test_three_intervals(self):
        rate_map = msprime.RateMap(position=[0, 5, 10, 15], rate=[0.1, 0.1, 0.1])
        assert rate_map.find_index(0) == 0
        assert rate_map.find_index(0.0001) == 0
        assert rate_map.find_index(4.9999) == 0
        assert rate_map.find_index(5) == 1
        assert rate_map.find_index(5.1) == 1
        assert rate_map.find_index(7) == 1
        assert rate_map.find_index(9.999) == 1
        assert rate_map.find_index(10) == 2
        assert rate_map.find_index(10.1) == 2
        assert rate_map.find_index(12) == 2
        assert rate_map.find_index(14.9999) == 2

    def test_out_of_bounds(self):
        rate_map = msprime.RateMap(position=[0, 10], rate=[0.1])
        for bad_value in [-1, -0.0001, 10, 10.0001, 1e9]:
            with pytest.raises(KeyError, match="out of bounds"):
                rate_map.find_index(bad_value)

    def test_input_types(self):
        rate_map = msprime.RateMap(position=[0, 10], rate=[0.1])
        assert rate_map.find_index(0) == 0
        assert rate_map.find_index(0.0) == 0
        assert rate_map.find_index(np.zeros(1)[0]) == 0


class TestSimpleExamples:
    def test_all_missing_one_interval(self):
        with pytest.raises(ValueError, match="missing data"):
            msprime.RateMap(position=[0, 10], rate=[np.nan])

    def test_all_missing_two_intervals(self):
        with pytest.raises(ValueError, match="missing data"):
            msprime.RateMap(position=[0, 5, 10], rate=[np.nan, np.nan])

    def test_count(self):
        rate_map = msprime.RateMap(position=[0, 5, 10], rate=[np.nan, 1])
        assert rate_map.num_intervals == 2
        assert rate_map.num_missing_intervals == 1
        assert rate_map.num_non_missing_intervals == 1

    def test_missing_arrays(self):
        rate_map = msprime.RateMap(position=[0, 5, 10], rate=[np.nan, 1])
        assert list(rate_map.missing) == [True, False]
        assert list(rate_map.non_missing) == [False, True]

    def test_missing_at_start_mean_rate(self):
        positions = np.array([0, 0.5, 1, 2])
        rates = np.array([np.nan, 0, 1])
        rate_map = msprime.RateMap(position=positions, rate=rates)
        assert np.isclose(rate_map.mean_rate, 1 / (1 + 0.5))

    def test_missing_at_end_mean_rate(self):
        positions = np.array([0, 1, 1.5, 2])
        rates = np.array([1, 0, np.nan])
        rate_map = msprime.RateMap(position=positions, rate=rates)
        assert np.isclose(rate_map.mean_rate, 1 / (1 + 0.5))

    def test_interval_properties_all_known(self):
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        assert list(rate_map.left) == [0, 1, 2]
        assert list(rate_map.right) == [1, 2, 3]
        assert list(rate_map.mid) == [0.5, 1.5, 2.5]
        assert list(rate_map.span) == [1, 1, 1]
        assert list(rate_map.mass) == [0.1, 0.2, 0.3]

    def test_pickle_non_missing(self):
        r1 = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        r2 = pickle.loads(pickle.dumps(r1))
        assert r1 == r2

    def test_pickle_missing(self):
        r1 = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, np.nan, 0.3])
        r2 = pickle.loads(pickle.dumps(r1))
        assert r1 == r2

    def test_get_cumulative_mass_all_known(self):
        rate_map = msprime.RateMap(position=[0, 10, 20, 30], rate=[0.1, 0.2, 0.3])
        assert list(rate_map.mass) == [1, 2, 3]
        assert list(rate_map.get_cumulative_mass([10, 20, 30])) == [1, 3, 6]

    def test_cumulative_mass_missing(self):
        rate_map = msprime.RateMap(position=[0, 10, 20, 30], rate=[0.1, np.nan, 0.3])
        assert list(rate_map.get_cumulative_mass([10, 20, 30])) == [1, 1, 4]


class TestDisplay:
    def test_str(self):
        rate_map = msprime.RateMap(position=[0, 10], rate=[0.1])
        s = """
        ┌──────────────────────────────────┐
        │left  │right  │  mid│  span│  rate│
        ├──────────────────────────────────┤
        │0     │10     │    5│    10│   0.1│
        └──────────────────────────────────┘
        """
        assert textwrap.dedent(s) == str(rate_map)

    def test_str_scinot(self):
        rate_map = msprime.RateMap(position=[0, 10], rate=[0.000001])
        s = """
        ┌───────────────────────────────────┐
        │left  │right  │  mid│  span│   rate│
        ├───────────────────────────────────┤
        │0     │10     │    5│    10│  1e-06│
        └───────────────────────────────────┘
        """
        assert textwrap.dedent(s) == str(rate_map)

    def test_repr(self):
        rate_map = msprime.RateMap(position=[0, 10], rate=[0.1])
        s = "RateMap(position=array([ 0., 10.]), rate=array([0.1]))"
        assert repr(rate_map) == s

    def test_repr_html(self):
        rate_map = msprime.RateMap(position=[0, 10], rate=[0.1])
        html = rate_map._repr_html_()
        root = xml.etree.ElementTree.fromstring(html)
        assert root.tag == "div"
        table = root.find("table")
        rows = list(table.find("tbody"))
        assert len(rows) == 1

    def test_long_table(self):
        n = 100
        rate_map = msprime.RateMap(position=range(n + 1), rate=[0.1] * n)
        headers, data = rate_map._display_table()
        assert len(headers) == 5
        assert len(data) == 21
        # check some left values
        assert int(data[0][0]) == 0
        assert int(data[-1][0]) == n - 1

    def test_short_table(self):
        n = 10
        rate_map = msprime.RateMap(position=range(n + 1), rate=[0.1] * n)
        headers, data = rate_map._display_table()
        assert len(headers) == 5
        assert len(data) == n
        # check some left values.
        assert int(data[0][0]) == 0
        assert int(data[-1][0]) == n - 1


class TestRateMapIsMapping:
    def test_items(self):
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        items = list(rate_map.items())
        assert items[0] == (0.5, 0.1)
        assert items[1] == (1.5, 0.2)
        assert items[2] == (2.5, 0.3)

    def test_keys(self):
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        assert list(rate_map.keys()) == [0.5, 1.5, 2.5]

    def test_values(self):
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        assert list(rate_map.values()) == [0.1, 0.2, 0.3]

    def test_in_points(self):
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        # Any point within the map are True
        for x in [0, 0.5, 1, 2.9999]:
            assert x in rate_map
        # Points outside the map are False
        for x in [-1, -0.0001, 3, 3.1]:
            assert x not in rate_map

    def test_in_slices(self):
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        # slices that are within the map are "in"
        for x in [slice(0, 0.5), slice(0, 1), slice(0, 2), slice(2, 3), slice(0, 3)]:
            assert x in rate_map
        # Any slice that doesn't fully intersect with the map "not in"
        assert slice(-0.001, 1) not in rate_map
        assert slice(0, 3.0001) not in rate_map
        assert slice(2.9999, 3.0001) not in rate_map
        assert slice(3, 4) not in rate_map
        assert slice(-2, -1) not in rate_map

    def test_other_types_not_in(self):
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        for other_type in [None, "sdf", "123", {}, [], Exception]:
            assert other_type not in rate_map

    def test_len(self):
        rate_map = msprime.RateMap(position=[0, 1], rate=[0.1])
        assert len(rate_map) == 1
        rate_map = msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2])
        assert len(rate_map) == 2
        rate_map = msprime.RateMap(position=[0, 1, 2, 3], rate=[0.1, 0.2, 0.3])
        assert len(rate_map) == 3

    def test_immutable(self):
        rate_map = msprime.RateMap(position=[0, 1], rate=[0.1])
        with pytest.raises(TypeError, match="item assignment"):
            rate_map[0] = 1
        with pytest.raises(TypeError, match="item deletion"):
            del rate_map[0]

    def test_eq(self):
        r1 = msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2])
        r2 = msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2])
        assert r1 == r1
        assert r1 == r2
        r2 = msprime.RateMap(position=[0, 1, 3], rate=[0.1, 0.2])
        assert r1 != r2
        assert msprime.RateMap(position=[0, 1], rate=[0.1]) != msprime.RateMap(
            position=[0, 1], rate=[0.2]
        )
        assert msprime.RateMap(position=[0, 1], rate=[0.1]) != msprime.RateMap(
            position=[0, 10], rate=[0.1]
        )

    def test_getitem_value(self):
        rate_map = msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2])
        assert rate_map[0] == 0.1
        assert rate_map[0.5] == 0.1
        assert rate_map[1] == 0.2
        assert rate_map[1.5] == 0.2
        assert rate_map[1.999] == 0.2
        # Try other types
        assert rate_map[np.array([1], dtype=np.float32)[0]] == 0.2
        assert rate_map[np.array([1], dtype=np.int32)[0]] == 0.2
        assert rate_map[np.array([1], dtype=np.float64)[0]] == 0.2
        assert rate_map[1 / 2] == 0.1
        assert rate_map[fractions.Fraction(1, 3)] == 0.1
        assert rate_map[decimal.Decimal(1)] == 0.2

    def test_getitem_slice(self):
        r1 = msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2])
        # The semantics of the slice() function are tested elsewhere.
        assert r1[:] == r1.copy()
        assert r1[:] is not r1
        assert r1[1:] == r1.slice(left=1)
        assert r1[:1.5] == r1.slice(right=1.5)
        assert r1[0.5:1.5] == r1.slice(left=0.5, right=1.5)

    def test_getitem_slice_step(self):
        r1 = msprime.RateMap(position=[0, 1, 2], rate=[0.1, 0.2])
        # Trying to set a "step" is a error
        with pytest.raises(TypeError, match="interval slicing"):
            r1[0:3:1]


class TestMappingMissingData:
    def test_get_missing(self):
        rate_map = msprime.RateMap(position=[0, 1, 2], rate=[np.nan, 0.2])
        with pytest.raises(KeyError, match="within a missing interval"):
            rate_map[0]
        with pytest.raises(KeyError, match="within a missing interval"):
            rate_map[0.999]

    def test_in_missing(self):
        rate_map = msprime.RateMap(position=[0, 1, 2], rate=[np.nan, 0.2])
        assert 0 not in rate_map
        assert 0.999 not in rate_map
        assert 1 in rate_map

    def test_keys_missing(self):
        rate_map = msprime.RateMap(position=[0, 1, 2], rate=[np.nan, 0.2])
        assert list(rate_map.keys()) == [1.5]


class TestGetIntermediates:
    def test_get_rate(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(position=positions, rate=rates)
        assert np.all(rate_map.get_rate([0.5, 1.5]) == rates)

    def test_get_rate_out_of_bounds(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(position=positions, rate=rates)
        with pytest.raises(ValueError, match="out of bounds"):
            rate_map.get_rate([1, -0.1])
        with pytest.raises(ValueError, match="out of bounds"):
            rate_map.get_rate([2])

    def test_get_cumulative_mass(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(position=positions, rate=rates)
        assert np.allclose(rate_map.get_cumulative_mass([0.5, 1.5]), np.array([0.5, 3]))
        assert rate_map.get_cumulative_mass([2]) == rate_map.total_mass

    def test_get_bad_cumulative_mass(self):
        positions = np.array([0, 1, 2])
        rates = np.array([1, 4])
        rate_map = msprime.RateMap(position=positions, rate=rates)
        with pytest.raises(ValueError, match="positions"):
            rate_map.get_cumulative_mass([1, -0.1])
        with pytest.raises(ValueError, match="positions"):
            rate_map.get_cumulative_mass([1, 2.1])


class TestSlice:
    def test_slice_no_params(self):
        # test RateMap.slice(..., trim=False)
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[0, 1, 2, 3])
        b = a.slice()
        assert a.sequence_length == b.sequence_length
        assert_array_equal(a.position, b.position)
        assert_array_equal(a.rate, b.rate)
        assert a == b

    def test_slice_left_examples(self):
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[0, 1, 2, 3])
        b = a.slice(left=50)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 50, 100, 200, 300, 400], b.position)
        assert_array_equal([np.nan, 0, 1, 2, 3], b.rate)

        b = a.slice(left=100)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 100, 200, 300, 400], b.position)
        assert_array_equal([np.nan, 1, 2, 3], b.rate)

        b = a.slice(left=150)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 150, 200, 300, 400], b.position)
        assert_array_equal([np.nan, 1, 2, 3], b.rate)

    def test_slice_right_examples(self):
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[0, 1, 2, 3])
        b = a.slice(right=300)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 100, 200, 300, 400], b.position)
        assert_array_equal([0, 1, 2, np.nan], b.rate)

        b = a.slice(right=250)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 100, 200, 250, 400], b.position)
        assert_array_equal([0, 1, 2, np.nan], b.rate)

    def test_slice_left_right_examples(self):
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[0, 1, 2, 3])
        b = a.slice(left=50, right=300)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 50, 100, 200, 300, 400], b.position)
        assert_array_equal([np.nan, 0, 1, 2, np.nan], b.rate)

        b = a.slice(left=150, right=250)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 150, 200, 250, 400], b.position)
        assert_array_equal([np.nan, 1, 2, np.nan], b.rate)

        b = a.slice(left=150, right=300)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 150, 200, 300, 400], b.position)
        assert_array_equal([np.nan, 1, 2, np.nan], b.rate)

        b = a.slice(left=150, right=160)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 150, 160, 400], b.position)
        assert_array_equal([np.nan, 1, np.nan], b.rate)

    def test_slice_right_missing(self):
        # If we take a right-slice into a trailing missing region,
        # we should recover the same map.
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[0, 1, 2, np.nan])
        b = a.slice(right=350)
        assert a.sequence_length == b.sequence_length
        assert_array_equal(a.position, b.position)
        assert_array_equal(a.rate, b.rate)

        b = a.slice(right=300)
        assert a.sequence_length == b.sequence_length
        assert_array_equal(a.position, b.position)
        assert_array_equal(a.rate, b.rate)

    def test_slice_left_missing(self):
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[np.nan, 1, 2, 3])
        b = a.slice(left=50)
        assert a.sequence_length == b.sequence_length
        assert_array_equal(a.position, b.position)
        assert_array_equal(a.rate, b.rate)

        b = a.slice(left=100)
        assert a.sequence_length == b.sequence_length
        assert_array_equal(a.position, b.position)
        assert_array_equal(a.rate, b.rate)

    def test_slice_with_floats(self):
        #  test RateMap.slice(..., trim=False) with floats
        a = msprime.RateMap(
            position=[np.pi * x for x in [0, 100, 200, 300, 400]], rate=[0, 1, 2, 3]
        )
        b = a.slice(left=50 * np.pi)
        assert a.sequence_length == b.sequence_length
        assert_array_equal([0, 50 * np.pi] + list(a.position[1:]), b.position)
        assert_array_equal([np.nan] + list(a.rate), b.rate)

    def test_slice_trim_left(self):
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[1, 2, 3, 4])
        b = a.slice(left=100, trim=True)
        assert b == msprime.RateMap(position=[0, 100, 200, 300], rate=[2, 3, 4])
        b = a.slice(left=50, trim=True)
        assert b == msprime.RateMap(position=[0, 50, 150, 250, 350], rate=[1, 2, 3, 4])

    def test_slice_trim_right(self):
        a = msprime.RateMap(position=[0, 100, 200, 300, 400], rate=[1, 2, 3, 4])
        b = a.slice(right=300, trim=True)
        assert b == msprime.RateMap(position=[0, 100, 200, 300], rate=[1, 2, 3])
        b = a.slice(right=350, trim=True)
        assert b == msprime.RateMap(position=[0, 100, 200, 300, 350], rate=[1, 2, 3, 4])

    def test_slice_error(self):
        recomb_map = msprime.RateMap(position=[0, 100], rate=[1])
        with pytest.raises(KeyError):
            recomb_map.slice(left=-1)
        with pytest.raises(KeyError):
            recomb_map.slice(right=-1)
        with pytest.raises(KeyError):
            recomb_map.slice(left=200)
        with pytest.raises(KeyError):
            recomb_map.slice(right=200)
        with pytest.raises(KeyError):
            recomb_map.slice(left=20, right=10)


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
        assert_array_equal(rm.position, [0, 1, 2, 3])
        assert np.allclose(rm.rate, [np.nan, 1e-8, 5e-8], equal_nan=True)

    def test_read_hapmap_from_filename(self, tmp_path):
        with open(tmp_path / "hapfile.txt", "w") as hapfile:
            hapfile.write(
                """\
                HEADER
                chr1 1 x 0
                chr1 2 x 0.000001 x
                chr1 3 x 0.000006 x x"""
            )
        rm = msprime.RateMap.read_hapmap(tmp_path / "hapfile.txt")
        assert_array_equal(rm.position, [0, 1, 2, 3])
        assert np.allclose(rm.rate, [np.nan, 1e-8, 5e-8], equal_nan=True)

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
        assert_array_equal(rm.position, [0, 1, 2])
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
        assert_array_equal(rm.position, [0, 1, 2])
        assert_array_equal(rm.rate, [np.nan, 5e-8])

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
        assert_array_equal(rm.position, [0, 1, 2])
        assert_array_equal(rm.rate, [1e-8, 5.5e-8])

    def test_read_hapmap_nonzero_map_start(self):
        hapfile = io.StringIO(
            """\
            HEADER
            chr1 1 x 0.000001
            chr1 2 x 0.000001 x
            chr1 3 x 0.000006 x x x"""
        )
        rm = msprime.RateMap.read_hapmap(hapfile)
        assert_array_equal(rm.position, [0, 1, 2, 3])
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
        assert_array_equal(rm.position, [0, 1, 2])
        assert np.allclose(rm.rate, [1e-8, 5e-8])

        hapfile.seek(0)
        rm = msprime.RateMap.read_hapmap(hapfile, sequence_length=10)
        assert_array_equal(rm.position, [0, 1, 2, 10])
        assert np.allclose(rm.rate, [1e-8, 5e-8, np.nan], equal_nan=True)

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
        assert_array_equal(rm1.rate, rm2.rate)
        assert_array_equal(rm1.position, rm2.position)

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
        assert np.allclose(rm1.position, rm2.position)
        assert np.allclose(rm1.rate, rm2.rate, equal_nan=True)
