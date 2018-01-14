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
from __future__ import print_function
from __future__ import division

import sys
import threading
import unittest
import platform

import numpy as np

import msprime
import tests.tsutil as tsutil

IS_PY2 = sys.version_info[0] < 3
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


class TestLdCalculatorReplicates(unittest.TestCase):
    """
    Tests the LdCalculator object to ensure we get correct results
    when using threads.
    """
    num_test_sites = 25

    def get_tree_sequence(self):
        ts = msprime.simulate(20, mutation_rate=10, recombination_rate=10, random_seed=8)
        return tsutil.subsample_sites(ts, self.num_test_sites)

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


# @unittest.skipIf(IS_PY2, "Cannot test thread support on Py2.")
# Temporarily skipping these on windows too. See
# https://github.com/jeromekelleher/msprime/issues/344
@unittest.skipIf(IS_PY2 or IS_WINDOWS, "Cannot test thread support on Py2.")
class TestTables(unittest.TestCase):
    """
    Tests to ensure that attempts to access tables in threads correctly
    raise an exception.
    """
    def get_tables(self):
        # TODO include migrations here.
        ts = msprime.simulate(
            100, mutation_rate=10, recombination_rate=10, random_seed=8)
        return ts.tables

    def run_multiple_writers(self, writer, num_writers=32):
        barrier = threading.Barrier(num_writers)

        def writer_proxy(thread_index, results):
            barrier.wait()
            # Attempts to operate on a table while locked should raise a RuntimeError
            try:
                writer(thread_index, results)
                results[thread_index] = 0
            except RuntimeError:
                results[thread_index] = 1

        results = run_threads(writer_proxy, num_writers)
        failures = sum(results)
        successes = num_writers - failures
        # Note: we would like to insist that #failures is > 0, but this is too
        # stochastic to guarantee for test purposes.
        self.assertGreaterEqual(failures, 0)
        self.assertGreater(successes, 0)

    def run_failing_reader(self, writer, reader, num_readers=32):
        """
        Runs a test in which a single writer acceses some tables
        and a bunch of other threads try to read the data.
        """
        barrier = threading.Barrier(num_readers + 1)

        def writer_proxy():
            barrier.wait()
            writer()

        def reader_proxy(thread_index, results):
            barrier.wait()
            # Attempts to operate on a table while locked should raise a RuntimeError
            try:
                reader(thread_index, results)
                results[thread_index] = 0
            except RuntimeError:
                results[thread_index] = 1

        writer_thread = threading.Thread(target=writer_proxy)
        writer_thread.start()
        results = run_threads(reader_proxy, num_readers)
        writer_thread.join()

        failures = sum(results)
        successes = num_readers - failures
        # Note: we would like to insist that #failures is > 0, but this is too
        # stochastic to guarantee for test purposes.
        self.assertGreaterEqual(failures, 0)
        self.assertGreater(successes, 0)

    def test_many_simplify_nodes_edges(self):
        tables = self.get_tables()

        def writer(thread_index, results):
            msprime.simplify_tables([0, 1], nodes=tables.nodes, edges=tables.edges)

        self.run_multiple_writers(writer)

    def test_many_simplify_all_tables(self):
        tables = self.get_tables()

        def writer(thread_index, results):
            msprime.simplify_tables(
                [0, 1], nodes=tables.nodes, edges=tables.edges, sites=tables.sites,
                mutations=tables.mutations)

        self.run_multiple_writers(writer)

    def test_many_sort(self):
        tables = self.get_tables()

        def writer(thread_index, results):
            msprime.sort_tables(**tables.asdict())

        self.run_multiple_writers(writer)

    def run_simplify_access_table(self, table_name, col_name):
        tables = self.get_tables()

        def writer():
            msprime.simplify_tables(
                [0, 1], nodes=tables.nodes, edges=tables.edges,
                sites=tables.sites, mutations=tables.mutations)

        table = getattr(tables, table_name)

        def reader(thread_index, results):
            for j in range(100):
                x = getattr(table, col_name)
                assert x.shape[0] == len(table)

        self.run_failing_reader(writer, reader)

    def run_sort_access_table(self, table_name, col_name):
        tables = self.get_tables()

        def writer():
            msprime.sort_tables(**tables.asdict())

        table = getattr(tables, table_name)

        def reader(thread_index, results):
            for j in range(100):
                x = getattr(table, col_name)
                assert x.shape[0] == len(table)

        self.run_failing_reader(writer, reader)

    def test_simplify_access_nodes(self):
        self.run_simplify_access_table("nodes", "time")

    def test_simplify_access_edges(self):
        self.run_simplify_access_table("edges", "left")

    def test_simplify_access_sites(self):
        self.run_simplify_access_table("sites", "position")

    def test_simplify_access_mutations(self):
        self.run_simplify_access_table("mutations", "site")

    def test_sort_access_nodes(self):
        self.run_sort_access_table("nodes", "time")

    def test_sort_access_edges(self):
        self.run_sort_access_table("edges", "left")

    def test_sort_access_sites(self):
        self.run_sort_access_table("sites", "position")

    def test_sort_access_mutations(self):
        self.run_sort_access_table("mutations", "site")
