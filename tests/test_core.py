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
