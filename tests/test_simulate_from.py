#
# Copyright (C) 2018 University of Oxford
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
Test cases for the simulate-from functionality.
"""
from __future__ import print_function
from __future__ import division

import unittest
import itertools

import numpy as np

import msprime
import _msprime
import tests.tsutil as tsutil


class TestBasicFunctionality(unittest.TestCase):
    """
    Basic tests for the from_ts argument for msprime.simulate.
    """
    def verify_from_tables(self, from_ts, final_ts, start_time):
        from_tables = from_ts.dump_tables()
        final_tables = final_ts.dump_tables()
        # Populations and individuals should be equal.
        self.assertEqual(from_tables.populations, final_tables.populations)
        self.assertEqual(from_tables.individuals, final_tables.individuals)
        # Time for new nodes > start_time
        new_time = final_tables.nodes.time[from_ts.num_nodes:]
        self.assertTrue(np.all(new_time > start_time))
        # Other tables should be equal up to the from_tables.
        final_tables.nodes.truncate(len(from_tables.nodes))
        self.assertEqual(final_tables.nodes, from_tables.nodes)
        final_tables.edges.truncate(len(from_tables.edges))
        self.assertEqual(final_tables.edges, from_tables.edges)
        # The mutation_rate parameter in simulate is not permitted, so we
        # should always have the same set of mutations before and after.
        self.assertEqual(final_tables.sites, from_tables.sites)
        self.assertEqual(final_tables.mutations, from_tables.mutations)
        final_tables.provenances.truncate(len(from_tables.provenances))
        self.assertEqual(final_tables.provenances, from_tables.provenances)

    def verify_simulation_completed(self, ts):
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 1)

    def test_from_single_locus_decapitated(self):
        ts = msprime.simulate(10, random_seed=5)
        from_ts = tsutil.decapitate(ts, ts.num_edges // 2)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_single_locus_max_time(self):
        from_ts = msprime.simulate(20, __tmp_max_time=1, random_seed=5)
        self.assertGreater(max(tree.num_roots for tree in from_ts.trees()), 1)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_single_locus_mutations(self):
        from_ts = msprime.simulate(
            20, __tmp_max_time=1, random_seed=5, mutation_rate=5)
        self.assertGreater(max(tree.num_roots for tree in from_ts.trees()), 1)
        self.assertGreater(from_ts.num_sites, 0)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_decapitated_mutations(self):
        ts = msprime.simulate(10, random_seed=5, mutation_rate=10)
        from_ts = tsutil.decapitate(ts, ts.num_edges // 2)
        self.assertGreater(from_ts.num_mutations, 0)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_from_multi_locus_decapitated(self):
        ts = msprime.simulate(10, recombination_rate=2, random_seed=5)
        self.assertGreater(ts.num_trees, 1)
        from_ts = tsutil.decapitate(ts, ts.num_edges // 2)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_from_multi_locus_max_time(self):
        from_ts = msprime.simulate(
            10, recombination_rate=2, random_seed=5, __tmp_max_time=1)
        self.assertTrue(any(tree.num_roots > 1 for tree in from_ts.trees()))
        self.assertGreater(from_ts.num_trees, 1)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_random_seeds_equal_outcome(self):
        from_ts = msprime.simulate(
            8, recombination_rate=2, random_seed=5, __tmp_max_time=1)
        self.assertGreater(from_ts.num_trees, 1)
        self.assertTrue(any(tree.num_roots > 1 for tree in from_ts.trees()))
        start_time = from_ts.tables.nodes.time.max()
        seed = 234
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=seed)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)
        final_tables = final_ts.dump_tables()
        final_tables.provenances.clear()
        for _ in range(10):
            other_ts = msprime.simulate(
                from_ts=from_ts, start_time=start_time, random_seed=seed)
            other_tables = other_ts.dump_tables()
            other_tables.provenances.clear()
            self.assertEqual(final_tables, other_tables)

    def test_individuals(self):
        from_ts = msprime.simulate(5, random_seed=5, __tmp_max_time=0.5)
        self.assertTrue(any(tree.num_roots > 1 for tree in from_ts.trees()))
        from_ts = tsutil.insert_random_ploidy_individuals(from_ts, seed=2)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_replicates(self):
        from_ts = msprime.simulate(
            15, recombination_rate=2, random_seed=5, __tmp_max_time=2)
        self.assertTrue(any(tree.num_roots > 1 for tree in from_ts.trees()))
        self.assertGreater(from_ts.num_trees, 1)
        start_time = from_ts.tables.nodes.time.max()
        replicates = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2, num_replicates=10)
        tables = []
        for final_ts in replicates:
            self.verify_from_tables(from_ts, final_ts, start_time)
            self.verify_simulation_completed(final_ts)
            tables.append(final_ts.dump_tables())
            tables[-1].provenances.clear()
        for a, b in itertools.combinations(tables, 2):
            self.assertNotEqual(a, b)

    def test_mutations_not_allowed(self):
        from_ts = msprime.simulate(15, random_seed=5, __tmp_max_time=2)
        start_time = from_ts.tables.nodes.time.max()
        with self.assertRaises(ValueError):
            msprime.simulate(
                from_ts=from_ts, start_time=start_time, mutation_rate=10)


class TestErrors(unittest.TestCase):
    """
    Basic tests for the from_ts argument for msprime.simulate.
    """
    def get_example_base(self):
        ts = msprime.simulate(10, random_seed=5)
        return tsutil.decapitate(ts, ts.num_edges // 2)

    def test_samples(self):
        base_ts = self.get_example_base()
        self.assertRaises(ValueError, msprime.simulate, 2, from_ts=base_ts)
        self.assertRaises(ValueError, msprime.simulate, sample_size=2, from_ts=base_ts)
        self.assertRaises(
            ValueError, msprime.simulate,
            samples=[msprime.Sample(0, 0) for _ in range(10)], from_ts=base_ts)
        self.assertRaises(
            ValueError, msprime.simulate,
            population_configurations=[
                msprime.PopulationConfiguration(sample_size=2)],
            from_ts=base_ts)

    def test_bad_type(self):
        for bad_type in [{}, 1, "asd"]:
            self.assertRaises(
                TypeError, msprime.simulate, from_ts=bad_type, start_time=1)

    def test_no_start_time(self):
        base_ts = self.get_example_base()
        self.assertRaises(ValueError, msprime.simulate, from_ts=base_ts)

    def test_start_time_less_than_zero(self):
        base_ts = self.get_example_base()
        with self.assertRaises(_msprime.InputError):
            msprime.simulate(from_ts=base_ts, start_time=-1)

    def test_start_time_less_than_base_nodes(self):
        base_ts = self.get_example_base()
        max_time = max(node.time for node in base_ts.nodes())
        for x in [0, max_time - 1, max_time - 1e-6]:
            with self.assertRaises(_msprime.InputError):
                msprime.simulate(from_ts=base_ts, start_time=x)

    def test_all_population_ids_null(self):
        base_ts = self.get_example_base()
        tables = base_ts.dump_tables()
        nodes = tables.nodes
        nodes.set_columns(
            flags=nodes.flags,
            time=nodes.time)
        with self.assertRaises(_msprime.InputError):
            msprime.simulate(from_ts=tables.tree_sequence(), start_time=nodes.time.max())
        nodes.set_columns(
            flags=nodes.flags,
            population=np.zeros_like(nodes.population),
            time=nodes.time)
        final_ts = msprime.simulate(
            from_ts=tables.tree_sequence(), start_time=nodes.time.max())
        self.assertEqual(
            sum(tree.num_roots for tree in final_ts.trees()), final_ts.num_trees)

    def test_single_population_id_null(self):
        base_ts = self.get_example_base()
        tables = base_ts.dump_tables()
        nodes = tables.nodes

        for j in range(base_ts.num_nodes):
            population = np.zeros_like(nodes.population)
            population[j] = -1
            nodes.set_columns(
                flags=nodes.flags,
                population=population,
                time=nodes.time)
            with self.assertRaises(_msprime.InputError):
                msprime.simulate(
                    from_ts=tables.tree_sequence(), start_time=nodes.time.max())
