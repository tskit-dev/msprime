#
# Copyright (C) 2015-2020 University of Oxford
#
# This file is part of utils.
#
# utils.is free software: you can redistribute it and/or modify
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
import sys
import unittest

import msprime.core as core


# Convenience method for getting seeds in a subprocess.
def get_seed(x):
    return core.get_random_seed()


class TestDefaultRandomSeeds(unittest.TestCase):
    """
    Tests for the default random seed generator.
    """

    def test_seed_generator_init(self):
        core.clear_seed_rng()
        seed = core.get_random_seed()
        self.assertGreater(seed, 0)
        self.assertIsNotNone(core.get_seed_rng())

    def test_unique(self):
        n = 100
        core.clear_seed_rng()
        seeds1 = [core.get_random_seed() for _ in range(n)]
        self.assertEqual(len(set(seeds1)), n)
        seeds2 = [core.get_random_seed() for _ in range(n)]
        self.assertEqual(len(set(seeds2)), n)
        self.assertEqual(len(set(seeds2)) + len(set(seeds2)), 2 * n)

    def test_unique_multiple_processes_no_init(self):
        n = 100
        core.clear_seed_rng()
        # Would use with block here, but not supported in Py < 3.3.
        pool = multiprocessing.Pool(5)
        seeds = pool.map(get_seed, range(n))
        self.assertEqual(len(set(seeds)), n)
        pool.terminate()
        pool.join()

    def test_unique_multiple_processes_init(self):
        n = 100
        core.get_random_seed()
        self.assertIsNotNone(core.get_seed_rng())
        # Would use with block here, but not supported in Py < 3.3.
        pool = multiprocessing.Pool(5)
        seeds = pool.map(get_seed, range(n))
        self.assertEqual(len(set(seeds)), n)
        pool.terminate()
        pool.join()


class TestAlmostEqual(unittest.TestCase):
    """
    Simple tests to ensure that the almost_equal() method is sensible.
    """

    def test_defaults(self):
        eps = sys.float_info.epsilon
        equal = [(1, 1), (0, 0), (1 + eps, 1), (1, 1 - eps), (10.000000000001, 10.0)]
        for a, b in equal:
            self.assertAlmostEqual(a, b)
            self.assertTrue(core.almost_equal(a, b))

    def test_near_zero(self):
        eps = sys.float_info.epsilon
        equal = [(0, 0), (eps, 0), (0, -eps), (-eps, eps)]
        for a, b in equal:
            self.assertAlmostEqual(a, b)
            self.assertTrue(core.almost_equal(a, b, abs_tol=1e-9))
        not_equal = [(0, 0.0000001), (-0.0000001, 0)]
        for a, b in not_equal:
            self.assertNotAlmostEqual(a, b)
            self.assertFalse(core.almost_equal(a, b, abs_tol=1e-9))
