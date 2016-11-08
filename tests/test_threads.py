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
import _msprime


def run_threads(worker, num_threads):
    _exception_occured = [False for _ in range(num_threads)]

    def local_worker(thread_index, results):
        try:
            worker(thread_index, results)
        except Exception as e:
            _exception_occured[thread_index] = True
            raise e

    # TODO remove the mandatory results array here and just use local
    # variables instead.
    results = [None for _ in range(num_threads)]
    threads = [
        threading.Thread(target=local_worker, args=(j, results))
        for j in range(num_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if any(_exception_occured):
        raise Exception("Error occured in tests!!")
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


class TestTreeSequenceThreads(unittest.TestCase):
    """
    Tests to ensure that we obtain consistent results when running tree
    sequence methods with many threads.
    """
    num_threads = 15

    def get_tree_sequence(self):
        ts = msprime.simulate(
            50, mutation_rate=10, recombination_rate=5, random_seed=10)
        self.assertGreater(ts.num_mutations, 0)
        self.assertGreater(ts.num_trees, 1)
        return ts

    def test_set_mutations(self):
        ts = msprime.simulate(10)
        new_mutations = [(j / ts.sample_size, j) for j in range(ts.sample_size)]

        def worker(thread_index, unused_arg=None):
            ts.set_mutations(new_mutations)
        run_threads(worker, self.num_threads)
        self.assertEqual(
            [(mut.position, mut.node) for mut in ts.mutations()], new_mutations)

    def test_haplotypes(self):
        ts = self.get_tree_sequence()
        H = list(ts.haplotypes())

        def worker(thread_index, results):
            results[thread_index] = list(ts.haplotypes())

        results = run_threads(worker, self.num_threads)
        self.assertEqual(len(H), ts.sample_size)
        self.assertEqual(len(H[0]), ts.num_mutations)
        for result in results:
            self.assertEqual(H, result)

    def test_trees(self):
        ts = self.get_tree_sequence()
        T = []
        for t in ts.trees():
            T.append(t.parent_dict)

        def worker(thread_index, results):
            trees = []
            for t in ts.trees():
                trees.append(t.parent_dict)
            results[thread_index] = trees

        results = run_threads(worker, self.num_threads)
        for result in results:
            self.assertEqual(T, result)

    def test_variants(self):
        ts = self.get_tree_sequence()
        V = np.zeros((ts.num_mutations, ts.sample_size), dtype=int)
        for variant in ts.variants():
            V[variant.index] = variant.genotypes

        def worker(thread_index, results):
            variants = np.zeros((ts.num_mutations, ts.sample_size), dtype=int)
            for variant in ts.variants():
                variants[variant.index] = variant.genotypes
            results[thread_index] = variants

        results = run_threads(worker, self.num_threads)
        for result in results:
            self.assertTrue(np.all(result == V))


class TestSparseTreeThreads(unittest.TestCase):
    """
    Test that the sparse tree can be used across multiple threads.
    """
    num_threads = 20

    def test_getters(self):
        ts = msprime.simulate(10, mutation_rate=5)
        t1 = next(ts.trees())

        def worker(thread_index, results):
            t2 = next(ts.trees())
            # t2 is a thread private copy of the same tree.
            self.assertIsNot(t1, t2)
            self.assertEqual(t1.root, t2.root)
            self.assertEqual(t1.time(t1.root), t2.time(t2.root))
            self.assertEqual(list(t1.nodes()), list(t2.nodes()))
            for j in range(ts.sample_size):
                self.assertEqual(t1.mrca(0, j), t2.mrca(0, j))
                self.assertEqual(t1.tmrca(0, j), t2.tmrca(0, j))

        run_threads(worker, self.num_threads)

    def test_set_mutations(self):
        ts = msprime.simulate(10, recombination_rate=5, mutation_rate=5)
        self.assertGreater(ts.num_trees, 2)
        tree = next(ts.trees())

        def worker(thread_index, results):
            mutations = [(j * 0.1, j) for j in range(5)]
            self.assertRaises(_msprime.LibraryError, ts.set_mutations, mutations)
            j = 0
            for t in ts.trees():
                j += 1
            self.assertEqual(j, ts.num_trees)

        run_threads(worker, self.num_threads)

        del tree
        mutations = [(j * 0.1, j) for j in range(5)]
        ts.set_mutations(mutations)
        self.assertEqual(
            mutations, [(mut.position, mut.node) for mut in ts.mutations()])

    def test_refcounts(self):
        ts = msprime.simulate(10, recombination_rate=5, mutation_rate=5)
        self.assertGreater(ts.num_trees, 2)
        ll_ts = ts.get_ll_tree_sequence()
        self.assertEqual(ll_ts.get_reference_count(), 0)
        num_loops = 100

        def worker(thread_index, results):
            n = 10
            for _ in range(num_loops):
                trees = [next(ts.trees()) for _ in range(n)]
                ld_calcs = [msprime.LdCalculator(ts) for _ in range(n)]
                ld_calcs = []
                j = 0
                for t in ts.trees():
                    j += 1
                self.assertEqual(j, ts.num_trees)
                for u in ts.trees():
                    if u.index == 2:
                        break
                del trees
                del ld_calcs
                del t
                del u
        run_threads(worker, self.num_threads)
        self.assertEqual(ll_ts.get_reference_count(), 0)
