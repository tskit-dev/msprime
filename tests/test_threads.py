#
# Copyright (C) 2016-2017 University of Oxford
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
Test cases for threading enabled aspects of the API.
"""
import threading
import unittest
import platform

import msprime

IS_WINDOWS = platform.system() == "Windows"


def run_threads(worker, num_threads):
    results = [None for _ in range(num_threads)]
    threads = [
        threading.Thread(target=worker, args=(j, results))
        for j in range(num_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    return results


class TestSimulatorThreads(unittest.TestCase):
    """
    Tests that we can run simulations in separate threads and
    get the same results.
    """
    num_threads = 10

    def test_single_replicate_equality(self):

        def worker(thread_index, results):
            results[thread_index] = msprime.simulate(100, random_seed=10)

        results = run_threads(worker, self.num_threads)
        # Make sure the results are all the same.
        records = list(results[0].records())
        self.assertGreater(len(records), 0)
        for ts in results[1:]:
            self.assertEqual(records, list(ts.records()))

    def test_multiple_replicate_equality(self):
        num_replicates = 10

        def worker(thread_index, results):
            results[thread_index] = []
            iterator = msprime.simulate(
                10, random_seed=10, num_replicates=num_replicates)
            for ts in iterator:
                results[thread_index].append(list(ts.records()))

        results = run_threads(worker, self.num_threads)
        # Make sure the results are all the same.
        self.assertEqual(len(results[0]), num_replicates)
        self.assertGreater(len(results[0][0]), 0)
        for result in results[1:]:
            self.assertEqual(results[0], result)
