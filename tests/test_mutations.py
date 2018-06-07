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
from __future__ import print_function
from __future__ import division

import unittest
import json

import msprime


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
        for bad_type in [None, {}, "234"]:
            self.assertRaises(TypeError, msprime.mutate, ts, rate=bad_type)
        for bad_rate in [-1, -1e-6, -1e7]:
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

    def test_provenance(self):
        ts = msprime.simulate(10, random_seed=1)
        for mutation_rate in [0, 1, 1e-5]:
            mutated = msprime.mutate(ts, mutation_rate)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            self.assertEqual(record["command"], "mutate")
            self.assertEqual(record["parameters"]["rate"], mutation_rate)
            self.assertTrue(record["parameters"]["random_seed"] >= 0)

        for seed in range(1, 10):
            mutated = msprime.mutate(ts, rate=1, random_seed=seed)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            self.assertEqual(record["command"], "mutate")
            self.assertEqual(record["parameters"]["rate"], 1)
            self.assertEqual(record["parameters"]["random_seed"], seed)

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
        tables = [ts.dump_tables() for ts in mutated]
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
