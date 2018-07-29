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

    def test_single_locus_old_recombination(self):
        from_ts = msprime.simulate(20, __tmp_max_time=1, random_seed=5)
        self.assertGreater(max(tree.num_roots for tree in from_ts.trees()), 1)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2,
            recombination_rate=2)
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

    def test_from_multi_locus_old_recombination(self):
        ts = msprime.simulate(10, recombination_rate=2, random_seed=5)
        self.assertGreater(ts.num_trees, 1)
        from_ts = tsutil.decapitate(ts, ts.num_edges // 2)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2,
            recombination_rate=2)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_from_two_populations_max_time(self):
        from_ts = msprime.simulate(
            recombination_rate=2, random_seed=5, __tmp_max_time=1,
            population_configurations=[
                msprime.PopulationConfiguration(5),
                msprime.PopulationConfiguration(5)],
            migration_matrix=[[0, 1], [1, 0]])
        self.assertTrue(any(tree.num_roots > 1 for tree in from_ts.trees()))
        self.assertGreater(from_ts.num_trees, 1)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time,
            random_seed=2, recombination_rate=2,
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration()],
            migration_matrix=[[0, 1], [1, 0]])
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_from_two_populations_decapitated(self):
        from_ts = msprime.simulate(
            recombination_rate=2, random_seed=5,
            population_configurations=[
                msprime.PopulationConfiguration(5),
                msprime.PopulationConfiguration(5)],
            migration_matrix=[[0, 1], [1, 0]])
        from_ts = tsutil.decapitate(from_ts, from_ts.num_edges // 3)
        self.assertTrue(any(tree.num_roots > 1 for tree in from_ts.trees()))
        self.assertGreater(from_ts.num_trees, 1)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time,
            random_seed=2, recombination_rate=2,
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration()],
            migration_matrix=[[0, 1], [1, 0]])
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_many_populations(self):
        for N in range(1, 6):
            population_configurations = [
                msprime.PopulationConfiguration() for _ in range(N)]
            migration_matrix = np.ones((N, N))
            np.fill_diagonal(migration_matrix, 0)
            from_ts = msprime.simulate(
                samples=[msprime.Sample(0, 0) for _ in range(10)],
                recombination_rate=2, random_seed=15,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix)
            from_ts = tsutil.decapitate(from_ts, from_ts.num_edges // 3)
            self.assertTrue(any(tree.num_roots > 1 for tree in from_ts.trees()))
            self.assertGreater(from_ts.num_trees, 1)
            start_time = from_ts.tables.nodes.time.max()
            final_ts = msprime.simulate(
                from_ts=from_ts, start_time=start_time,
                random_seed=8, recombination_rate=2,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix)
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
        from_ts = msprime.simulate(25, random_seed=5, __tmp_max_time=0.5)
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

    def test_from_subclass(self):
        from_ts = msprime.simulate(20, __tmp_max_time=1, random_seed=5)

        class MockTreeSequence(msprime.TreeSequence):
            pass

        subclass_instance = MockTreeSequence(from_ts.ll_tree_sequence)
        self.assertTrue(type(subclass_instance), MockTreeSequence)
        self.assertIsInstance(subclass_instance, msprime.TreeSequence)
        self.assertGreater(max(tree.num_roots for tree in subclass_instance.trees()), 1)
        start_time = subclass_instance.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=subclass_instance, start_time=start_time, random_seed=2)
        self.verify_from_tables(subclass_instance, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    @unittest.skip("Loss of precision with sequence length, single locus.")
    def test_sequence_length(self):
        from_ts = msprime.simulate(
            5, __tmp_max_time=0.1, random_seed=5, length=5)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2, length=5)
        print("FROM")
        print(from_ts.tables.edges)
        for tree in from_ts.trees():
            print(tree.draw(format="unicode"))
        print("FINAL")
        for tree in final_ts.trees():
            print(tree.interval)
            print(tree.draw(format="unicode"))
        print(final_ts.tables.edges)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_sequence_length_recombination(self):
        from_ts = msprime.simulate(
            5, __tmp_max_time=0.1, random_seed=5, length=5, recombination_rate=5)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2, length=5,
            recombination_rate=5)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    @unittest.skip("Assertion fail")
    def test_tricky_recombination_map(self):
        # Zero rates within the map cause tricky mapping issues
        positions = [0, 0.25, 0.5, 0.75, 1]
        rates = [1, 0, 1, 0, 0]
        num_loci = 100
        recomb_map = msprime.RecombinationMap(positions, rates, num_loci)
        from_ts = msprime.simulate(
            5, __tmp_max_time=0.1, random_seed=5, recombination_map=recomb_map)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=2,
            recombination_map=recomb_map)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)

    def test_random_recombination_map(self):
        np.random.seed(10)
        k = 10
        position = np.random.random(k) * 10
        position[0] = 0
        position.sort()
        rate = np.random.random(k)
        rate[-1] = 0
        recomb_map = msprime.RecombinationMap(list(position), list(rate))
        from_ts = msprime.simulate(
            10, __tmp_max_time=0.1, random_seed=50, recombination_map=recomb_map)
        start_time = from_ts.tables.nodes.time.max()
        final_ts = msprime.simulate(
            from_ts=from_ts, start_time=start_time, random_seed=20,
            recombination_map=recomb_map)
        self.verify_from_tables(from_ts, final_ts, start_time)
        self.verify_simulation_completed(final_ts)


class TestBaseEquivalance(unittest.TestCase):
    """
    Check that it's equivalent to send a from_ts with no topology to running
    a straight simulation.
    """
    def verify_simple_model(
            self, n, seed=1, recombination_rate=None, length=None,
            recombination_map=None):
        ts1 = msprime.simulate(
            n, random_seed=seed, recombination_rate=recombination_rate,
            length=length, recombination_map=recombination_map)
        tables = msprime.TableCollection(ts1.sequence_length)
        tables.populations.add_row()
        for _ in range(n):
            tables.nodes.add_row(
                flags=msprime.NODE_IS_SAMPLE, time=0, population=0)
        ts2 = msprime.simulate(
            from_ts=tables.tree_sequence(), start_time=0, random_seed=seed,
            recombination_rate=recombination_rate,
            recombination_map=recombination_map)
        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()
        # print(tables1.edges)
        # print(tables2.edges)
        tables1.provenances.clear()
        tables2.provenances.clear()
        self.assertEqual(tables1, tables2)

    def test_single_locus_two_samples(self):
        for seed in range(1, 10):
            self.verify_simple_model(2, seed)

    def test_single_locus_five_samples(self):
        for seed in range(1, 10):
            self.verify_simple_model(5, seed)

    @unittest.skip("Major loss of precision in single locus case")
    def test_single_locus_sequence_length(self):
        for length in [0.1, 0.99, 5, 10, 33.333, 1000, 1e9]:
            self.verify_simple_model(5, 43, length=length)

    def test_multi_locus_two_samples(self):
        for seed in range(1, 10):
            self.verify_simple_model(2, seed, recombination_rate=1)

    def test_multi_locus_five_samples(self):
        for seed in range(1, 10):
            self.verify_simple_model(5, seed, recombination_rate=1)

    def test_multi_locus_sequence_length(self):
        for length in [0.1, 2.5, 4, 8, 33.33333]:
            self.verify_simple_model(5, 45, length=length, recombination_rate=1)

    def test_random_recombination_map(self):
        np.random.seed(100)
        k = 15
        position = np.random.random(k) * 100
        position[0] = 0
        position.sort()
        rate = np.random.random(k)
        rate[-1] = 0
        recomb_map = msprime.RecombinationMap(list(position), list(rate))
        self.verify_simple_model(10, 23, recombination_map=recomb_map)

    def test_random_recombination_map_small_num_loci(self):
        np.random.seed(100)
        k = 15
        position = np.random.random(k) * 100
        position[0] = 0
        position.sort()
        rate = np.random.random(k)
        rate[-1] = 0
        recomb_map = msprime.RecombinationMap(list(position), list(rate), num_loci=5)
        self.verify_simple_model(10, 23, recombination_map=recomb_map)

    def test_two_populations_migration(self):
        n = 10
        seed = 1234
        ts1 = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(n),
                msprime.PopulationConfiguration(0)],
            migration_matrix=[[0, 1], [1, 0]],
            record_migrations=True,
            random_seed=seed)
        tables = msprime.TableCollection(1)
        tables.populations.add_row()
        tables.populations.add_row()
        for _ in range(n):
            tables.nodes.add_row(
                flags=msprime.NODE_IS_SAMPLE, time=0, population=0)
        ts2 = msprime.simulate(
            from_ts=tables.tree_sequence(), start_time=0,
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration()],
            migration_matrix=[[0, 1], [1, 0]],
            record_migrations=True,
            random_seed=seed)
        tables1 = ts1.dump_tables()
        tables2 = ts2.dump_tables()
        tables1.provenances.clear()
        tables2.provenances.clear()
        self.assertEqual(tables1, tables2)


class TestErrors(unittest.TestCase):
    """
    Basic tests for the from_ts argument for msprime.simulate.
    """
    def get_example_base(self, num_populations=1, length=1):
        N = num_populations
        population_configurations = [
            msprime.PopulationConfiguration() for _ in range(N)]
        migration_matrix = np.ones((N, N))
        np.fill_diagonal(migration_matrix, 0)
        ts = msprime.simulate(
            samples=[msprime.Sample(0, 0) for _ in range(10)],
            length=length,
            random_seed=155,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix)
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

    def test_sequence_length_mismatch(self):
        base_ts = self.get_example_base(length=5)
        for bad_length in [1, 4.99, 5.01, 100]:
            with self.assertRaises(ValueError):
                msprime.simulate(from_ts=base_ts, start_time=100, length=bad_length)
            recomb_map = msprime.RecombinationMap.uniform_map(bad_length, 1)
            with self.assertRaises(ValueError):
                msprime.simulate(
                    from_ts=base_ts, start_time=100, recombination_map=recomb_map)

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

    def test_population_mismatch(self):
        for N in range(1, 5):
            base_ts = self.get_example_base(num_populations=N)
            start_time = max(node.time for node in base_ts.nodes())
            for k in range(1, N):
                if k != N:
                    with self.assertRaises(ValueError):
                        msprime.simulate(
                            from_ts=base_ts, start_time=start_time,
                            population_configurations=[
                                msprime.PopulationConfiguration() for _ in range(k)])

    def test_population_mismatch_no_population_configs(self):
        for N in range(2, 5):
            base_ts = self.get_example_base(num_populations=N)
            start_time = max(node.time for node in base_ts.nodes())
            with self.assertRaises(ValueError):
                msprime.simulate(from_ts=base_ts, start_time=start_time)


class TestSlimOutput(unittest.TestCase):
    """
    Verify that we can successfully simulate from SLiM output.
    """
    def finish_simulation(self, from_ts, recombination_rate=0, seed=1):
        population_configurations = [
            msprime.PopulationConfiguration() for _ in range(from_ts.num_populations)]
        return msprime.simulate(
            from_ts=from_ts, start_time=1,
            population_configurations=population_configurations,
            recombination_rate=recombination_rate,
            random_seed=seed)

    def verify_completed(self, from_ts, final_ts):
        from_tables = from_ts.dump_tables()
        final_tables = final_ts.dump_tables()
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
        print(from_ts.tables)
        print(final_ts.tables)
        self.assertEqual(max(tree.num_roots for tree in final_ts.trees()), 1)

    @unittest.skip("Single locus recomb map issue")
    def test_minimal_example_no_recombination(self):
        from_ts = msprime.load("tests/data/SLiM/minimal-example.trees")
        ts = self.finish_simulation(from_ts, recombination_rate=0, seed=1)
        self.verify_completed(from_ts, ts)

    def test_minimal_example_recombination(self):
        from_ts = msprime.load("tests/data/SLiM/minimal-example.trees")
        ts = self.finish_simulation(from_ts, recombination_rate=0.1, seed=1)
        self.verify_completed(from_ts, ts)
