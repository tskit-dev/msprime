#
# Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
from __future__ import print_function
from __future__ import division

import unittest
import threading
import random

import numpy as np

import msprime


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


class TestLdCalculatorReplicates(unittest.TestCase):
    """
    Tests the LdCalculator object to ensure we get correct results
    when using threads.
    """
    num_test_mutations = 25

    def get_tree_sequence(self):
        ts = msprime.simulate(
            20, mutation_rate=10, recombination_rate=10, random_seed=10)
        mutations = random.sample(
            list(ts.mutations()), self.num_test_mutations)
        ts.set_mutations(sorted(mutations))
        return ts

    def test_get_r2_multiple_instances(self):
        # This is the nominal case where we have a separate LdCalculator
        # instance in each thread.
        ts = self.get_tree_sequence()
        ld_calc = msprime.LdCalculator(ts)
        A = ld_calc.get_r2_matrix()
        del ld_calc
        m = A.shape[0]

        def worker(thread_index, results):
            ld_calc = msprime.LdCalculator(ts)
            row = np.zeros(m)
            results[thread_index] = row
            for j in range(m):
                row[j] = ld_calc.get_r2(thread_index, j)

        results = run_threads(worker, m)
        for j in range(m):
            self.assertTrue(np.allclose(results[j], A[j]))

    def test_get_r2_single_instance(self):
        # This is the degenerate case where we have a single LdCalculator
        # instance shared by the threads. We should have only one thread
        # actually executing get_r2() at one time.
        ts = self.get_tree_sequence()
        ld_calc = msprime.LdCalculator(ts)
        A = ld_calc.get_r2_matrix()
        m = A.shape[0]

        def worker(thread_index, results):
            row = np.zeros(m)
            results[thread_index] = row
            for j in range(m):
                row[j] = ld_calc.get_r2(thread_index, j)

        results = run_threads(worker, m)
        for j in range(m):
            self.assertTrue(np.allclose(results[j], A[j]))

    def test_get_r2_array_multiple_instances(self):
        # This is the nominal case where we have a separate LdCalculator
        # instance in each thread.
        ts = self.get_tree_sequence()
        ld_calc = msprime.LdCalculator(ts)
        A = ld_calc.get_r2_matrix()
        m = A.shape[0]
        del ld_calc

        def worker(thread_index, results):
            ld_calc = msprime.LdCalculator(ts)
            results[thread_index] = np.array(
                ld_calc.get_r2_array(thread_index))

        results = run_threads(worker, m)
        for j in range(m):
            self.assertTrue(np.allclose(results[j], A[j, j + 1:]))

    def test_get_r2_array_single_instance(self):
        # This is the degenerate case where we have a single LdCalculator
        # instance shared by the threads. We should have only one thread
        # actually executing get_r2_array() at one time. Because the buffer
        # is shared by many different instances, we can't make any assertions
        # about the returned values --- they are essentially gibberish.
        # However, we shouldn't crash and burn, which is what this test
        # is here to check for.
        ts = self.get_tree_sequence()
        ld_calc = msprime.LdCalculator(ts)
        m = ts.get_num_mutations()

        def worker(thread_index, results):
            results[thread_index] = ld_calc.get_r2_array(thread_index).shape

        results = run_threads(worker, m)
        for j in range(m):
            self.assertEqual(results[j][0], m - j - 1)
