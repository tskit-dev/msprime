#
# Copyright (C) 2015-2020 University of Oxford
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# utils.is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with utils.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Test cases for the high level interface to utils.
"""
import multiprocessing

import numpy as np

import msprime.core as core


# Convenience method for getting seeds in a subprocess.
def get_seed(x):
    return core.get_random_seed()


class TestDefaultRandomSeeds:
    """
    Tests for the default random seed generator.
    """

    def test_seed_generator_init(self):
        core.clear_seed_rng()
        seed = core.get_random_seed()
        assert seed > 0
        assert core.get_seed_rng() is not None

    def test_unique(self):
        n = 100
        core.clear_seed_rng()
        seeds1 = [core.get_random_seed() for _ in range(n)]
        assert len(set(seeds1)) == n
        seeds2 = [core.get_random_seed() for _ in range(n)]
        assert len(set(seeds2)) == n
        assert len(set(seeds2)) + len(set(seeds2)) == 2 * n

    def test_unique_multiple_processes_no_init(self):
        n = 100
        core.clear_seed_rng()
        with multiprocessing.Pool(5) as pool:
            seeds = pool.map(get_seed, range(n))
            assert len(set(seeds)) == n

    def test_unique_multiple_processes_init(self):
        n = 100
        core.get_random_seed()
        assert core.get_seed_rng() is not None
        with multiprocessing.Pool(5) as pool:
            seeds = pool.map(get_seed, range(n))
            assert len(set(seeds)) == n

    def test_set_seed_rng_seed(self):
        n = 10
        core.set_seed_rng_seed(42)
        seeds1 = [core.get_random_seed() for _ in range(n)]
        core.set_seed_rng_seed(42)
        seeds2 = [core.get_random_seed() for _ in range(n)]
        assert seeds1 == seeds2
        core.set_seed_rng_seed(1234)
        seeds3 = [core.get_random_seed() for _ in range(n)]
        assert seeds1 != seeds3


class TestIsInteger:
    """
    Tests for the function used to determine if a value can be interpreted
    as an integer.
    """

    def test_good_values(self):
        numpy_int_array = np.array([100], dtype=int)
        numpy_float_array = np.array([100], dtype=np.float64)
        good_values = [
            -1,
            0,
            10 ** 6,
            100_000,
            0x123,
            1.0,
            numpy_int_array[0],
            numpy_float_array[0],
        ]
        for good_value in good_values:
            assert core.isinteger(good_value)

    def test_bad_values(self):
        numpy_float_array = np.array([100.1], dtype=np.float64)
        bad_values = [
            [],
            None,
            {},
            np.array([10], dtype=int),
            numpy_float_array,
            numpy_float_array[0],
            1.1,
            "1",
            "100_000",
            "1.1",
            1e-3,
        ]
        for bad_value in bad_values:
            assert not core.isinteger(bad_value)
