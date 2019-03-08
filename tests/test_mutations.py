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
Test cases for the high level interface to msprime.
"""
import unittest
import json

import numpy as np

import msprime
from tests import tsutil
import tests.wright_fisher as wf


class TestMutateProvenance(unittest.TestCase):
    """
    Test that we correctly record the provenance for calls to mutate.
    """
    def test_mutation_rate(self):
        ts = msprime.simulate(10, random_seed=1)
        for mutation_rate in [0, 1, 1e-5]:
            mutated = msprime.mutate(ts, mutation_rate)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            self.assertEqual(record["parameters"]["command"], "mutate")
            self.assertEqual(record["parameters"]["rate"], mutation_rate)
            self.assertTrue(record["parameters"]["random_seed"] >= 0)

    def test_start_time(self):
        ts = msprime.simulate(10, random_seed=1)
        for start_time in [0, 1, -1]:
            mutated = msprime.mutate(ts, 1, start_time=start_time)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            self.assertEqual(record["parameters"]["command"], "mutate")
            self.assertEqual(record["parameters"]["start_time"], start_time)
            self.assertTrue(record["parameters"]["random_seed"] >= 0)

    def test_end_time(self):
        ts = msprime.simulate(10, random_seed=1)
        for end_time in [0, 1, 100]:
            mutated = msprime.mutate(ts, 1, start_time=-1, end_time=end_time)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            self.assertEqual(record["parameters"]["command"], "mutate")
            self.assertEqual(record["parameters"]["start_time"], -1)
            self.assertEqual(record["parameters"]["end_time"], end_time)
            self.assertTrue(record["parameters"]["random_seed"] >= 0)

    def test_seed(self):
        ts = msprime.simulate(10, random_seed=1)
        for seed in range(1, 10):
            mutated = msprime.mutate(ts, rate=1, random_seed=seed)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            self.assertEqual(record["parameters"]["command"], "mutate")
            self.assertEqual(record["parameters"]["rate"], 1)
            self.assertEqual(record["parameters"]["random_seed"], seed)

    def test_keep(self):
        ts = msprime.simulate(10, random_seed=1)
        for keep in [True, False]:
            mutated = msprime.mutate(ts, rate=1, keep=keep)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            self.assertEqual(record["parameters"]["command"], "mutate")
            self.assertEqual(record["parameters"]["keep"], keep)


class TestMutate(unittest.TestCase):
    """
    Tests the msprime.mutate function.
    """
    def verify_topology(self, source, dest):
        self.assertEqual(source.nodes, dest.nodes)
        self.assertEqual(source.edges, dest.edges)
        self.assertEqual(source.migrations, dest.migrations)
        self.assertEqual(source.populations, dest.populations)

    def verify_provenance(self, source, dest):
        before = list(source.provenances)
        after = list(dest.provenances)
        self.assertEqual(len(before) + 1, len(after))
        self.assertEqual(before, after[:-1])

    def test_zero_mutation_rate(self):
        ts = msprime.simulate(10, random_seed=1)
        mutated = msprime.mutate(ts, 0)
        t1 = ts.dump_tables()
        t2 = mutated.dump_tables()
        self.verify_topology(t1, t2)
        self.verify_provenance(t1, t2)
        self.assertEqual(t1.sites, t2.sites)
        self.assertEqual(t1.mutations, t2.mutations)

    def test_populations(self):
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(10),
                msprime.PopulationConfiguration(10)],
            migration_matrix=[[0, 1], [1, 0]],
            record_migrations=True,
            random_seed=1)
        mutated = msprime.mutate(ts, 0)
        t1 = ts.dump_tables()
        self.assertEqual(len(t1.populations), 2)
        self.assertGreater(len(t1.migrations), 0)
        t2 = mutated.dump_tables()
        self.verify_topology(t1, t2)
        self.verify_provenance(t1, t2)
        self.assertEqual(t1.sites, t2.sites)
        self.assertEqual(t1.mutations, t2.mutations)

    def test_mutation_overwrite(self):
        ts = msprime.simulate(10, mutation_rate=5, random_seed=2)
        self.assertGreater(ts.num_sites, 0)
        self.assertGreater(ts.num_mutations, 0)
        mutated = msprime.mutate(ts, 0)
        t1 = ts.dump_tables()
        self.assertEqual(len(t1.sites), ts.num_sites)
        t2 = mutated.dump_tables()
        self.verify_topology(t1, t2)
        self.verify_provenance(t1, t2)
        self.assertEqual(len(t2.sites), 0)
        self.assertEqual(len(t2.mutations), 0)

    def test_bad_rates(self):
        ts = msprime.simulate(2, random_seed=2)
        for bad_type in [{}, [], ts]:
            self.assertRaises(TypeError, msprime.mutate, ts, rate=bad_type)
        for bad_rate in ["abc", -1, -1e-6, -1e7]:
            self.assertRaises(ValueError, msprime.mutate, ts, bad_rate)

    def test_bad_models(self):
        ts = msprime.simulate(2, random_seed=2)
        for bad_type in [{}, "234"]:
            self.assertRaises(TypeError, msprime.mutate, ts, rate=0, model=bad_type)

    def test_bad_alphabet(self):
        for bad_alphabet in [-1, 2, 10**6, "1"]:
            self.assertRaises(ValueError, msprime.InfiniteSites, bad_alphabet)

    def test_bad_tree_sequence(self):
        for bad_type in [None, {}, "sdrf"]:
            self.assertRaises(ValueError, msprime.mutate, bad_type)

    def test_default_seeds(self):
        ts = msprime.simulate(20, random_seed=2)
        seeds = []
        for _ in range(10):
            mutated = msprime.mutate(ts, 0)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            seeds.append(record["parameters"]["random_seed"])
        self.assertEqual(len(seeds), len(set(seeds)))

    def test_identical_seed(self):
        ts = msprime.simulate(10, random_seed=2)
        mutated = [
            msprime.mutate(ts, rate=1, random_seed=2) for _ in range(1, 10)]
        self.assertGreater(mutated[0].num_sites, 0)
        self.assertGreater(mutated[0].num_mutations, 0)
        tables = [other_ts.dump_tables() for other_ts in mutated]
        self.assertTrue(all(tables[0].sites == t.sites for t in tables[1:]))
        self.assertTrue(all(tables[0].mutations == t.mutations for t in tables[1:]))

    def verify_binary_alphabet(self, ts):
        self.assertGreater(ts.num_sites, 0)
        self.assertTrue(all(site.ancestral_state == '0' for site in ts.sites()))
        self.assertTrue(
            all(mutation.derived_state == '1' for mutation in ts.mutations()))

    def verify_nucleotides_alphabet(self, ts):
        nucleotides = "ACGT"
        self.assertGreater(ts.num_sites, 0)
        self.assertTrue(all(site.ancestral_state in nucleotides for site in ts.sites()))
        self.assertTrue(
            all(mutation.derived_state in nucleotides for mutation in ts.mutations()))
        for site in ts.sites():
            self.assertNotEqual(site.ancestral_state, site.mutations[0].derived_state)

    def test_default_alphabet(self):
        ts = msprime.simulate(10, random_seed=2)
        mutated = msprime.mutate(ts, rate=1, random_seed=2)
        self.verify_binary_alphabet(mutated)

    def test_alphabet_binary(self):
        ts = msprime.simulate(10, random_seed=2)
        mutated = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.BINARY))
        self.verify_binary_alphabet(mutated)

    def test_alphabet_nucleotide(self):
        ts = msprime.simulate(10, random_seed=2)
        mutated = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.NUCLEOTIDES))
        self.verify_nucleotides_alphabet(mutated)

    def test_identical_seed_alphabets(self):
        ts = msprime.simulate(10, random_seed=2)
        binary = msprime.mutate(ts, rate=1, random_seed=2)
        nucleotides = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.NUCLEOTIDES))
        self.assertGreater(binary.num_sites, 0)
        self.assertGreater(binary.num_mutations, 0)
        self.assertEqual(binary.num_sites, nucleotides.num_sites)
        self.assertEqual(binary.num_mutations, nucleotides.num_mutations)
        for s1, s2 in zip(binary.sites(), nucleotides.sites()):
            self.assertEqual(s1.position, s2.position)
            self.assertEqual(s1.mutations[0].node, s2.mutations[0].node)


class TestInterval(unittest.TestCase):
    """
    Tests on the start_time and end_time parameters.
    """
    def test_errors(self):
        ts = msprime.simulate(10, random_seed=2)
        for start, end in [(-2, -3), (1, 0), (1e6, 1e5)]:
            with self.assertRaises(ValueError):
                msprime.mutate(ts, start_time=-2.0, end_time=-3.0)

    def test_stick_tree(self):
        tables = msprime.TableCollection(1.0)
        tables.nodes.add_row(flags=msprime.NODE_IS_SAMPLE, time=0)
        tables.nodes.add_row(flags=0, time=1)
        tables.nodes.add_row(flags=0, time=2)
        tables.edges.add_row(0, 1, 1, 0)
        tables.edges.add_row(0, 1, 2, 1)
        ts = tables.tree_sequence()

        tsm = msprime.mutate(ts, rate=100, end_time=1, random_seed=1)
        self.assertGreater(tsm.num_sites, 0)
        self.assertTrue(all(mut.node == 0 for mut in ts.mutations()))

        tsm = msprime.mutate(ts, rate=100, start_time=0, end_time=1, random_seed=1)
        self.assertGreater(tsm.num_sites, 0)
        self.assertTrue(all(mut.node == 0 for mut in ts.mutations()))

        tsm = msprime.mutate(ts, rate=100, start_time=0.5, end_time=1, random_seed=1)
        self.assertGreater(tsm.num_sites, 0)
        self.assertTrue(all(mut.node == 0 for mut in ts.mutations()))

        tsm = msprime.mutate(ts, rate=100, start_time=1, random_seed=1)
        self.assertGreater(tsm.num_sites, 0)
        self.assertTrue(all(mut.node == 1 for mut in ts.mutations()))

        tsm = msprime.mutate(ts, rate=100, start_time=1, end_time=2, random_seed=1)
        self.assertGreater(tsm.num_sites, 0)
        self.assertTrue(all(mut.node == 1 for mut in ts.mutations()))

        tsm = msprime.mutate(ts, rate=100, start_time=1.5, end_time=2, random_seed=1)
        self.assertGreater(tsm.num_sites, 0)
        self.assertTrue(all(mut.node == 0 for mut in ts.mutations()))

    def verify_mutations(self, ts, start, end):
        self.assertGreater(ts.num_sites, 0)
        time = ts.tables.nodes.time
        if start is None:
            start = np.min(time) - 1
        if end is None:
            end = np.min(time) + 1
        for tree in ts.trees():
            for site in tree.sites():
                for mutation in site.mutations:
                    top_node = tree.parent(mutation.node)
                    self.assertGreater(time[top_node], start)

    def verify(self, ts, rate=100):
        root_time = max(node.time for node in ts.nodes())
        leaf_time = min(node.time for node in ts.nodes())
        length = root_time - leaf_time

        end = root_time - length / 2
        tsm = msprime.mutate(ts, rate=rate, end_time=end)
        self.verify_mutations(tsm, None, end)

        start = leaf_time + length / 4
        end = root_time - length / 2
        tsm = msprime.mutate(ts, rate=rate, start_time=start, end_time=end)
        self.verify_mutations(tsm, start, end)

        start = root_time - length / 2
        end = root_time
        tsm = msprime.mutate(ts, rate=rate, start_time=start, end_time=end)
        self.verify_mutations(tsm, start, end)

        tsm = msprime.mutate(ts, rate=rate, start_time=start)
        self.verify_mutations(tsm, start, None)

    def test_coalescent_tree(self):
        ts = msprime.simulate(20, random_seed=2)
        self.verify(ts)

    def test_coalescent_trees(self):
        ts = msprime.simulate(20, recombination_rate=1, random_seed=2)
        self.verify(ts)

    def test_wright_fisher_trees(self):
        tables = wf.wf_sim(20, 10, seed=1, deep_history=False)
        tables.sort()
        ts = tables.tree_sequence()
        self.verify(ts, rate=10)
        ts = ts.simplify()
        self.verify(ts, rate=10)

    def test_negative_time(self):
        ts = msprime.simulate(10, recombination_rate=1, random_seed=2)
        tables = ts.dump_tables()
        time = tables.nodes.time
        max_time = np.max(time)
        for offset in [max_time / 2, max_time, 2 * max_time]:
            tables.nodes.set_columns(
                flags=tables.nodes.flags,
                time=time - offset)
            ts = tables.tree_sequence()
            self.verify(ts)


class TestKeep(unittest.TestCase):
    """
    Tests for the "keep" functionality in which we append new mutations
    to an existing set.
    """
    def verify(self, ts, rate, random_seed):
        no_keep = msprime.mutate(ts, rate=rate, random_seed=random_seed)
        self.assertGreater(no_keep.num_sites, 0)
        keep = msprime.mutate(ts, rate=rate, random_seed=random_seed, keep=True)
        # Can assume there's no collisions here, very unlikely.
        self.assertEqual(ts.num_sites + no_keep.num_sites, keep.num_sites)
        # Mutations are all infinite sites, so must be equal
        self.assertEqual(ts.num_mutations + no_keep.num_mutations, keep.num_mutations)
        old = {site.position for site in ts.sites()}
        new = {site.position for site in no_keep.sites()}
        both = {site.position for site in keep.sites()}
        self.assertEqual(old | new, both)
        self.verify_sites(ts, keep)

    def verify_sites(self, original, updated):
        site_map = {site.position: site for site in original.sites()}
        found = 0
        inserted_mutations = 0
        for site in updated.sites():
            if site.position in site_map:
                old_site = site_map[site.position]
                self.assertEqual(site.position, old_site.position)
                self.assertEqual(len(site.mutations), len(old_site.mutations))
                for mutation, old_mutation in zip(site.mutations, old_site.mutations):
                    self.assertEqual(mutation.metadata, old_mutation.metadata)
                    self.assertEqual(mutation.node, old_mutation.node)
                    self.assertEqual(mutation.derived_state, old_mutation.derived_state)
                    if old_mutation.parent == msprime.NULL_MUTATION:
                        self.assertEqual(mutation.parent, old_mutation.parent)
                    else:
                        self.assertEqual(
                            mutation.parent, old_mutation.parent + inserted_mutations)
                self.assertEqual(site.ancestral_state, old_site.ancestral_state)
                self.assertEqual(site.metadata, old_site.metadata)
                found += 1
            else:
                inserted_mutations += len(site.mutations)
        self.assertEqual(found, original.num_sites)

    def test_simple_binary(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=2)
        self.assertGreater(ts.num_sites, 0)
        self.verify(ts, 1, random_seed=2)

    def test_simple_nucleotide(self):
        ts = msprime.mutate(
            msprime.simulate(10, random_seed=2),
            rate=1,
            random_seed=2,
            model=msprime.InfiniteSites(msprime.NUCLEOTIDES))
        self.assertGreater(ts.num_sites, 0)
        self.verify(ts, 2, random_seed=3)

    def test_branch_mutations(self):
        ts = tsutil.insert_branch_mutations(
            msprime.simulate(10, recombination_rate=1, random_seed=2))
        self.assertGreater(ts.num_sites, 1)
        self.verify(ts, 3, random_seed=7)

    def test_multichar_mutations(self):
        ts = tsutil.insert_multichar_mutations(
            msprime.simulate(12, recombination_rate=4, random_seed=3))
        self.assertGreater(ts.num_sites, 5)
        self.verify(ts, 3, random_seed=7)

    def test_random_metadata(self):
        ts = tsutil.add_random_metadata(
            msprime.simulate(12, random_seed=3, mutation_rate=1))
        self.assertGreater(ts.num_sites, 5)
        self.verify(ts, 3, random_seed=7)

    def test_no_sites(self):
        ts = msprime.simulate(12, random_seed=3)
        self.assertEqual(ts.num_sites, 0)
        self.verify(ts, 3, random_seed=7)

    def test_same_seeds(self):
        ts = msprime.simulate(12, random_seed=3)
        self.assertEqual(ts.num_sites, 0)
        ts = msprime.mutate(ts, rate=1, random_seed=1)
        updated = msprime.mutate(ts, rate=1, random_seed=1, keep=True)
        # We should rejection sample away all the sites that we have already
        # generated so we just get another random sample.
        sites = set()
        for site in updated.sites():
            self.assertNotIn(site.position, sites)
            sites.add(site.position)
        for site in ts.sites():
            self.assertIn(site.position, sites)

    def test_keep_multichar_muts(self):
        ts = msprime.simulate(12, random_seed=3)
        ts = msprime.mutate(ts, rate=1, random_seed=1)
        self.assertGreater(ts.num_sites, 2)
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        for site in ts.sites():
            tables.sites.add_row(position=site.position, ancestral_state="A" * site.id)
            for mutation in site.mutations:
                tables.mutations.add_row(
                    site=site.id, node=mutation.node, derived_state="T" * site.id)
        original = tables.tree_sequence()
        updated = msprime.mutate(original, rate=1, random_seed=1, keep=True)
        self.verify_sites(original, updated)

    def test_keep_metadata(self):
        ts = msprime.simulate(12, random_seed=3)
        ts = msprime.mutate(ts, rate=1, random_seed=1)
        self.assertGreater(ts.num_sites, 2)
        # Set metadata on this ts so that we can be sure we keep the original
        # mutations.
        ts = tsutil.add_random_metadata(ts)
        other = msprime.mutate(ts, rate=1, random_seed=1, keep=True)
        self.verify_sites(ts, other)

    def test_keep_mutation_parent(self):
        ts = msprime.simulate(12, recombination_rate=3, random_seed=3)
        ts = tsutil.insert_branch_mutations(ts)
        self.assertGreater(ts.num_sites, 2)
        other = msprime.mutate(ts, rate=1, random_seed=1, keep=True)
        self.assertGreater(other.num_sites, ts.num_sites)
        self.verify_sites(ts, other)

    def test_keep_mutation_parent_zero_rate(self):
        ts = msprime.simulate(12, recombination_rate=3, random_seed=3)
        ts = tsutil.insert_branch_mutations(ts)
        self.assertGreater(ts.num_sites, 2)
        other = msprime.mutate(ts, rate=0, random_seed=1, keep=True)
        t1 = ts.dump_tables()
        t2 = other.dump_tables()
        t1.provenances.clear()
        t2.provenances.clear()
        self.assertEqual(t1, t2)
