#
# Copyright (C) 2018-2020 University of Oxford
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
Test cases for mutation generation.
"""
import functools
import itertools
import json
import unittest

import attr
import numpy as np
import scipy.stats as stats
import tskit

import _msprime
import msprime
import tests.wright_fisher as wf
from tests import tsutil


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


class TestMutationModel(unittest.TestCase):
    def validate_model(self, model):
        num_alleles = len(model.alleles)
        self.assertEqual(num_alleles, len(model.root_distribution))
        self.assertEqual((num_alleles, num_alleles), model.transition_matrix.shape)
        self.assertEqual(sum(model.root_distribution), 1.0)
        for j in range(num_alleles):
            self.assertEqual(sum(model.transition_matrix[j]), 1.0)
        s = model.__str__()
        self.assertIsInstance(s, str)
        self.assertTrue(s.startswith("Mutation model with alleles"))

    def validate_stationary(self, model):
        # for most models, the root distribution should be the
        # stationary distribution of the transition matrix
        x = np.matmul(model.root_distribution, model.transition_matrix)
        self.assertSequenceEqual(list(x), list(model.root_distribution))

    def validate_asdict(self, model):
        m = msprime.MutationModel(**model.asdict())
        self.assertEqual(model, m)

    def test_bad_alleles(self):
        for alleles, err in [
            (0, TypeError),
            (b"xyz", TypeError),
            (["0", "1"], TypeError),
            ([b"0", "1"], TypeError),
            ([1, 2, 3], TypeError),
            ([b"a"], _msprime.LibraryError),
        ]:
            with self.assertRaises(err):
                msprime.MutationModel(
                    alleles=alleles,
                    root_distribution=[1] + [0] * (len(alleles) - 1),
                    transition_matrix=[[1 / len(alleles)] * len(alleles)]
                    * len(alleles),
                )

    def test_bad_root_distribution(self):
        alleles = [b"a", b"b", b"c"]
        transition_matrix = [[1 / len(alleles)] * len(alleles)] * len(alleles)
        for root_distribution, err in [
            ([1, 0], ValueError),
            ([1, 0, 0, 0], ValueError),
            (1.0, ValueError),
            ("xyz", ValueError),
            (None, ValueError),
            ([1.5, -0.5, 0], _msprime.LibraryError),
            ([0.5, 0, 0], _msprime.LibraryError),
        ]:
            with self.assertRaises(err):
                msprime.MutationModel(
                    alleles=alleles,
                    root_distribution=root_distribution,
                    transition_matrix=transition_matrix,
                )

    def test_bad_transition_matrix(self):
        alleles = [b"a", b"b", b"c"]
        root_distribution = [0.5, 0.25, 0.25]
        for transition_matrix, err in [
            ([[1, 0]] * 3, ValueError),
            ([[1, 0, 0]] * 2, ValueError),
            (1.0, ValueError),
            ("xyz", ValueError),
            (None, ValueError),
            ([[1.5, -0.5, 0]] * 3, _msprime.LibraryError),
            ([[0.5, 0, 0]] * 3, _msprime.LibraryError),
        ]:
            with self.assertRaises(err):
                msprime.MutationModel(
                    alleles=alleles,
                    root_distribution=root_distribution,
                    transition_matrix=transition_matrix,
                )

    def test_binary(self):
        model = msprime.BinaryMutations()
        self.validate_model(model)

    def test_jukes_cantor(self):
        model = msprime.JukesCantor()
        self.validate_model(model)
        self.validate_stationary(model)

    def test_new_model(self):
        alleles = [b"Alligator", b"Camel", b"Goat"]
        root_distribution = [0.5, 0.25, 0.25]
        transition_matrix = [[0.0, 0.5, 0.5], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        model = msprime.MutationModel(alleles, root_distribution, transition_matrix)
        self.assertListEqual(model.alleles, alleles)
        self.assertListEqual(list(model.root_distribution), root_distribution)
        for a, b in zip(model.transition_matrix, transition_matrix):
            self.assertListEqual(list(a), b)
        self.validate_model(model)
        self.validate_stationary(model)


class MutateMixin:
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

    def verify_binary_alphabet(self, ts):
        self.assertGreater(ts.num_sites, 0)
        self.assertTrue(all(site.ancestral_state == "0" for site in ts.sites()))
        self.assertTrue(
            all(mutation.derived_state == "1" for mutation in ts.mutations())
        )

    def verify_nucleotides_alphabet(self, ts):
        nucleotides = "ACGT"
        self.assertGreater(ts.num_sites, 0)
        self.assertTrue(all(site.ancestral_state in nucleotides for site in ts.sites()))
        self.assertTrue(
            all(mutation.derived_state in nucleotides for mutation in ts.mutations())
        )
        for site in ts.sites():
            self.assertNotEqual(site.ancestral_state, site.mutations[0].derived_state)


class TestMutate(unittest.TestCase, MutateMixin):
    """
    Tests the msprime.mutate function.
    """

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
                msprime.PopulationConfiguration(10),
            ],
            migration_matrix=[[0, 1], [1, 0]],
            record_migrations=True,
            random_seed=1,
        )
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
        for bad_rate in ["abc", "xxx"]:
            self.assertRaises(ValueError, msprime.mutate, ts, bad_rate)
        for bad_rate in [-1, -1e-6, -1e7]:
            self.assertRaises(_msprime.LibraryError, msprime.mutate, ts, bad_rate)

    def test_bad_models(self):
        ts = msprime.simulate(2, random_seed=2)
        for bad_type in [{}, "234"]:
            self.assertRaises(TypeError, msprime.mutate, ts, rate=0, model=bad_type)

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
        mutated = [msprime.mutate(ts, rate=1, random_seed=2) for _ in range(1, 10)]
        self.assertGreater(mutated[0].num_sites, 0)
        self.assertGreater(mutated[0].num_mutations, 0)
        tables = [other_ts.dump_tables() for other_ts in mutated]
        self.assertTrue(all(tables[0].sites == t.sites for t in tables[1:]))
        self.assertTrue(all(tables[0].mutations == t.mutations for t in tables[1:]))

    def test_default_alphabet(self):
        ts = msprime.simulate(10, random_seed=2)
        mutated = msprime.mutate(ts, rate=1, random_seed=2)
        self.verify_binary_alphabet(mutated)

    def test_deprecated_alphabet_binary(self):
        ts = msprime.simulate(10, random_seed=2)
        mutated = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.BINARY)
        )
        self.verify_binary_alphabet(mutated)

    def test_deprecated_alphabet_nucleotides(self):
        ts = msprime.simulate(10, random_seed=2)
        mutated = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.NUCLEOTIDES)
        )
        self.verify_nucleotides_alphabet(mutated)

    def test_deprecated_bad_alphabet(self):
        ts = msprime.simulate(10, random_seed=2)
        with self.assertRaises(ValueError):
            msprime.mutate(ts, rate=1, random_seed=2, model=msprime.InfiniteSites(-1))
        with self.assertRaises(ValueError):
            msprime.mutate(ts, rate=1, random_seed=2, model=msprime.InfiniteSites(2))

    def test_deprecated_identical_seed_alphabets(self):
        ts = msprime.simulate(10, random_seed=2)
        binary = msprime.mutate(ts, rate=1, random_seed=2)
        nucs = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.NUCLEOTIDES)
        )
        self.assertGreater(binary.num_sites, 0)
        self.assertGreater(binary.num_mutations, 0)
        self.assertEqual(binary.num_sites, nucs.num_sites)
        self.assertEqual(binary.num_mutations, nucs.num_mutations)
        for s1, s2 in zip(binary.sites(), nucs.sites()):
            self.assertEqual(s1.position, s2.position)
            self.assertEqual(s1.mutations[0].node, s2.mutations[0].node)


class TestFiniteSites(TestMutate):
    def verify_binary_alphabet(self, ts):
        binary = "01"
        self.assertGreater(ts.num_sites, 0)
        self.assertTrue(all(site.ancestral_state == "0" for site in ts.sites()))
        self.assertTrue(
            all(mutation.derived_state in binary for mutation in ts.mutations())
        )

    def verify_mutations(self, ts, discrete, check_probs=False, model=None):
        # Verify that mutation.parents are correct
        # and that mutations actually involve a change of state
        # that, if (not keep), has positive probability under the model
        if check_probs:
            assert model is not None
            alleles = [x.decode() for x in model.alleles]
        for site in ts.sites():
            if not discrete:
                self.assertEqual(len(site.mutations), 1)
            t = ts.at(site.position)
            parents = {}
            for mut in site.mutations:
                n = mut.node
                # TODO: once mutations have times
                # self.assertGreaterEqual(mut.time, ts.node(n).time)
                # self.assertLess(mut.time, ts.node(t.parent(n)).time)
                while n != tskit.NULL and n not in parents:
                    n = t.parent(n)
                if n == tskit.NULL:
                    self.assertEqual(mut.parent, tskit.NULL)
                    pa = site.ancestral_state
                else:
                    self.assertEqual(mut.parent, parents[n].id)
                    # TODO: once mutations have times
                    # self.assertLess(mut.time, parents[n].time)
                    pa = parents[n].derived_state
                self.assertNotEqual(mut.derived_state, pa)
                if check_probs:
                    self.assertTrue(pa in alleles)
                    self.assertTrue(mut.derived_state in alleles)
                    i = alleles.index(pa)
                    j = alleles.index(mut.derived_state)
                    self.assertGreater(model.transition_matrix[i][j], 0)
                # do this last since parents should be listed before children
                parents[mut.node] = mut

    def mutate(self, ts, model, rate=1.0, keep=False, discrete=False, seed=42):
        lib_ts = msprime.mutate(
            ts, rate=rate, random_seed=seed, model=model, keep=keep, discrete=discrete
        )
        self.verify_mutations(lib_ts, discrete, check_probs=not keep, model=model)
        py_ts = py_mutate(
            ts, rate=rate, random_seed=seed, model=model, keep=keep, discrete=discrete
        )
        lib_tables = lib_ts.dump_tables()
        lib_tables.provenances.clear()
        py_tables = py_ts.dump_tables()
        py_tables.provenances.clear()
        self.assertEqual(lib_tables, py_tables)
        return lib_ts

    def mutate_binary(
        self,
        ts,
        rate=1.0,
        keep=False,
        discrete=False,
        transition_matrix=None,
        root_distribution=None,
    ):
        alleles = [b"0", b"1"]
        if transition_matrix is None:
            transition_matrix = [[0.0, 1.0], [1.0, 0.0]]
        if root_distribution is None:
            root_distribution = [1.0, 0.0]
        model = msprime.MutationModel(alleles, root_distribution, transition_matrix)
        return self.mutate(ts, model, rate=rate, keep=keep, discrete=discrete)

    def mutate_nucleotides(
        self,
        ts,
        rate=1.0,
        keep=False,
        discrete=False,
        transition_matrix=None,
        root_distribution=None,
    ):
        alleles = [b"A", b"C", b"G", b"T"]
        if transition_matrix is None:
            transition_matrix = [[0.25] * 4] * 4
        if root_distribution is None:
            root_distribution = [0.25] * 4
        model = msprime.MutationModel(alleles, root_distribution, transition_matrix)
        model = msprime.MutationModel(alleles, root_distribution, transition_matrix)
        return self.mutate(ts, model, rate=rate, keep=keep, discrete=discrete)

    def test_alleles_binary(self):
        ts = msprime.simulate(10, random_seed=1)
        mutated = self.mutate_binary(ts)
        self.verify_binary_alphabet(mutated)

    def test_alleles_nucleotides(self):
        ts = msprime.simulate(10, random_seed=1)
        mutated = self.mutate_nucleotides(ts)
        self.verify_nucleotides_alphabet(mutated)

    def test_zero_mutation_rate(self):
        for keep in (True, False):
            ts = msprime.simulate(10, random_seed=1, mutation_rate=2.0)
            mutated = self.mutate_binary(ts, 0, keep=keep)
            t1 = ts.dump_tables()
            t2 = mutated.dump_tables()
            self.verify_topology(t1, t2)
            self.verify_provenance(t1, t2)
            if keep:
                self.assertEqual(t1.sites, t2.sites)
                self.assertEqual(t1.mutations, t2.mutations)
            else:
                self.assertEqual(t2.sites.num_rows, 0)

    def test_bad_mutate_order(self):
        ts = msprime.simulate(10, random_seed=1, recombination_rate=1, length=10)
        mutated = msprime.mutate(
            ts, 3, random_seed=5, start_time=0.0, end_time=0.5, discrete=True
        )
        with self.assertRaises(_msprime.LibraryError):
            msprime.mutate(
                mutated,
                3,
                random_seed=6,
                start_time=0.5,
                end_time=1.0,
                discrete=True,
                keep=True,
            )

    def test_one_way_mutation(self):
        for discrete in (True, False):
            ts = msprime.simulate(10, random_seed=1, recombination_rate=1.0, length=10)
            mut_matrix = [[0.0, 1.0], [0.0, 1.0]]
            mutated = self.mutate_binary(
                ts,
                rate=1.0,
                transition_matrix=mut_matrix,
                root_distribution=[1.0, 0.0],
                discrete=discrete,
            )
            self.assertGreater(mutated.num_sites, 0)
            if discrete:
                self.assertGreater(max([len(s.mutations) for s in mutated.sites()]), 1)
            for site in mutated.sites():
                self.assertEqual(site.ancestral_state, "0")
                self.assertGreaterEqual(len(site.mutations), 1)
                for mut in site.mutations:
                    self.assertEqual(mut.derived_state, "1")

    def test_flip_flop_mutation(self):
        nucleotides = "ACGT"
        for discrete in (True, False):
            ts = msprime.simulate(10, random_seed=1, recombination_rate=1.0, length=10)
            mut_matrix = [
                [0.0, 0.0, 0.5, 0.5],
                [0.0, 0.0, 0.5, 0.5],
                [0.5, 0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0, 0.0],
            ]
            mutated = self.mutate_nucleotides(
                ts, rate=5.0, transition_matrix=mut_matrix, discrete=discrete
            )
            self.assertGreater(mutated.num_sites, 0)
            if discrete:
                self.assertGreater(max([len(s.mutations) for s in mutated.sites()]), 1)
            for mut in ts.mutations():
                self.assertTrue(mut.derived_state in nucleotides)
                if mut.parent == -1:
                    pa = ts.site(mut.site).ancestral_state
                else:
                    pa = ts.mutation(mut.parent).derived_state
                if pa in "AC":
                    self.assertTrue(mut.derived_state in "GT")
                else:
                    self.assertTrue(mut.derived_state in "AC")

    def test_do_nothing_mutations(self):
        for discrete in (True, False):
            ts = msprime.simulate(10, random_seed=1, length=10)
            mut_matrix = [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
            mutated = self.mutate_nucleotides(
                ts, rate=5.0, transition_matrix=mut_matrix, discrete=discrete
            )
            self.assertEqual(mutated.num_sites, 0)
            t1 = ts.dump_tables()
            t2 = mutated.dump_tables()
            self.verify_topology(t1, t2)
            self.verify_provenance(t1, t2)
            self.assertEqual(t1.sites, t2.sites)
            self.assertEqual(t1.mutations, t2.mutations)

    def test_uniform_mutations(self):
        nucleotides = "ACGT"
        for discrete in (True, False):
            ts = msprime.simulate(10, random_seed=1, recombination_rate=1.0, length=10)
            mut_matrix = [
                [0.1, 0.3, 0.3, 0.3],
                [0.0, 0.0, 0.5, 0.5],
                [0.0, 0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5, 0.0],
            ]
            mutated = self.mutate_nucleotides(
                ts,
                rate=20.0,
                transition_matrix=mut_matrix,
                root_distribution=[1.0, 0.0, 0.0, 0.0],
                discrete=discrete,
            )
            self.assertGreater(mutated.num_sites, 0)
            if discrete:
                self.assertGreater(max([len(s.mutations) for s in mutated.sites()]), 1)
            num_nucs = {a: 0 for a in nucleotides}
            for site in mutated.sites():
                self.assertTrue(site.ancestral_state == "A")
            for mut in mutated.mutations():
                num_nucs[mut.derived_state] += 1
            self.assertEqual(num_nucs["A"], 0)
            self.assertGreater(num_nucs["C"], 0)
            self.assertGreater(num_nucs["G"], 0)
            self.assertGreater(num_nucs["T"], 0)

    def test_circular_mutations(self):
        nucleotides = "ACGT"
        for discrete in (True, False):
            ts = msprime.simulate(10, random_seed=1, recombination_rate=1.0, length=10)
            mut_matrix = [
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0, 0.0],
            ]
            mutated = self.mutate_nucleotides(
                ts,
                rate=2.0,
                transition_matrix=mut_matrix,
                root_distribution=[1.0, 0.0, 0.0, 0.0],
                discrete=discrete,
            )
            self.assertGreater(mutated.num_sites, 0)
            if discrete:
                self.assertGreater(max([len(s.mutations) for s in mutated.sites()]), 1)
            for site in mutated.sites():
                self.assertTrue(site.ancestral_state == "A")
                self.assertGreater(len(site.mutations), 0)
            for mut in mutated.mutations():
                if mut.parent == tskit.NULL:
                    self.assertEqual(mut.derived_state, "C")
                else:
                    pmut = mutated.mutation(mut.parent)
                    self.assertEqual(
                        (nucleotides.index(pmut.derived_state) + 1) % 4,
                        nucleotides.index(mut.derived_state),
                    )

    def test_integer_sites(self):
        ts = msprime.simulate(10, random_seed=5, length=10, recombination_rate=10.0)
        mutated = self.mutate_binary(ts, rate=5.0, discrete=True)
        self.assertEqual(mutated.site(0).position, 0.0)
        for site in mutated.sites():
            self.assertEqual(site.position, np.floor(site.position))
            self.assertLess(site.position, mutated.sequence_length)


class TestInterval(unittest.TestCase):
    """
    Tests on the start_time and end_time parameters.
    """

    def test_errors(self):
        ts = msprime.simulate(10, random_seed=2)
        for start, end in [(-2, -3), (1, 0), (1e6, 1e5)]:
            with self.assertRaises(ValueError):
                msprime.mutate(ts, start_time=start, end_time=end)

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
            tables.nodes.set_columns(flags=tables.nodes.flags, time=time - offset)
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
                            mutation.parent, old_mutation.parent + inserted_mutations
                        )
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

    def test_deprecated_simple_nucleotide(self):
        ts = msprime.mutate(
            msprime.simulate(10, random_seed=2),
            rate=1,
            random_seed=2,
            model=msprime.InfiniteSites(msprime.NUCLEOTIDES),
        )
        self.assertGreater(ts.num_sites, 0)
        self.verify(ts, 2, random_seed=3)

    def test_branch_mutations(self):
        ts = tsutil.insert_branch_mutations(
            msprime.simulate(10, recombination_rate=1, random_seed=2)
        )
        self.assertGreater(ts.num_sites, 1)
        self.verify(ts, 3, random_seed=7)

    def test_multichar_mutations(self):
        ts = tsutil.insert_multichar_mutations(
            msprime.simulate(12, recombination_rate=4, random_seed=3)
        )
        self.assertGreater(ts.num_sites, 5)
        self.verify(ts, 3, random_seed=7)

    def test_random_metadata(self):
        ts = tsutil.add_random_metadata(
            msprime.simulate(12, random_seed=3, mutation_rate=1)
        )
        self.assertGreater(ts.num_sites, 5)
        self.verify(ts, 3, random_seed=7)

    def test_no_sites(self):
        ts = msprime.simulate(12, random_seed=3)
        self.assertEqual(ts.num_sites, 0)
        self.verify(ts, 3, random_seed=7)

    def test_site_with_no_mutation(self):
        ts = tsutil.insert_site(msprime.simulate(12, random_seed=3))
        self.assertEqual(ts.num_sites, 1)
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
                    site=site.id, node=mutation.node, derived_state="T" * site.id
                )
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


class StatisticalTestMixin:

    p_threshold = 0.01

    def chisquare(self, observed, expected):
        obs = observed.reshape((np.product(observed.shape)),)
        exp = expected.reshape((np.product(observed.shape)),)
        not_zeros = exp > 0
        self.assertSequenceEqual(
            list(obs[np.logical_not(not_zeros)]),
            list(np.zeros((len(obs) - sum(not_zeros),))),
        )
        if sum(not_zeros) > 1:
            chisq = stats.chisquare(obs[not_zeros], exp[not_zeros])
            self.assertGreater(chisq.pvalue, self.p_threshold)

    def bonferroni(self, pvals):
        self.assertGreater(min(pvals), self.p_threshold / len(pvals))

    def sign_tst(self, diffs):
        np = sum(diffs > 0)
        nm = sum(diffs < 0)
        pp = stats.binom.cdf(np, n=np + nm, p=0.5)
        pm = stats.binom.cdf(nm, n=np + nm, p=0.5)
        self.assertGreater(min(pp, pm), self.p_threshold)


class TestMutationStatistics(unittest.TestCase, StatisticalTestMixin):
    def verify_model(self, model, verify_roots=True):
        ots = msprime.simulate(10, random_seed=5, recombination_rate=0.05, length=20)
        for discrete in (True, False):
            for rate in (1, 10):
                ts = msprime.mutate(
                    ots, random_seed=6, rate=rate, model=model, discrete=discrete
                )
                self.verify_model_ts(ts, model, discrete, verify_roots)

    def verify_model_ts(self, ts, model, discrete, verify_roots):
        alleles = [a.decode() for a in model.alleles]
        num_alleles = len(alleles)
        roots = np.zeros((num_alleles,))
        transitions = np.zeros((num_alleles, num_alleles))
        for site in ts.sites():
            aa = site.ancestral_state
            roots[alleles.index(aa)] += 1
        for mut in ts.mutations():
            if mut.parent == tskit.NULL:
                pa = ts.site(mut.site).ancestral_state
            else:
                pa = ts.mutation(mut.parent).derived_state
            da = mut.derived_state
            transitions[alleles.index(pa), alleles.index(da)] += 1
        num_muts = np.sum(roots)
        if verify_roots:
            exp_roots = num_muts * model.root_distribution
            self.chisquare(roots, exp_roots)
        for j, (row, p) in enumerate(zip(transitions, model.transition_matrix)):
            p[j] = 0
            p /= sum(p)
            self.chisquare(row, sum(row) * p)

    def verify_mutation_rates(self, model):
        # this test only works if the probability of dropping a mutation
        # doesn't depend on the previous state
        assert len(set(np.diag(model.transition_matrix))) == 1
        np.random.seed(23)
        ots = msprime.simulate(10, random_seed=5, recombination_rate=0.05, length=20)
        for discrete in (False, True):
            for rate in (1, 10):
                ts = msprime.mutate(
                    ots, random_seed=6, rate=rate, model=model, discrete=discrete
                )
                self.verify_mutation_rates_ts(ts, model, rate, discrete)

    def verify_mutation_rates_ts(self, ts, model, rate, discrete):
        mut_rate = rate * (1 - model.transition_matrix[0, 0])
        lower_pvals = np.zeros(ts.num_trees)
        upper_pvals = np.zeros(ts.num_trees)
        diffs = np.zeros(ts.num_trees)
        for j, t in enumerate(ts.trees()):
            if discrete:
                span = np.ceil(t.interval[1]) - np.ceil(t.interval[0])
            else:
                span = t.span
            mean = mut_rate * span * t.total_branch_length
            N = t.num_mutations
            lower_pvals[j] = stats.poisson.cdf(mean, N)
            upper_pvals[j] = 1.0 - stats.poisson.cdf(mean, N + 1)
            # if we draw an indepenent Poisson with the same mean
            # it should be greater than observed half the time it is different
            # NOTE: more powerful might be a Wilcoxon test?
            rand = stats.poisson.rvs(mean, 1)
            diffs[j] = rand - N
        self.sign_tst(diffs)
        self.bonferroni(lower_pvals)
        self.bonferroni(upper_pvals)

    def test_binary_model(self):
        model = msprime.BinaryMutations()
        self.verify_model(model)
        self.verify_mutation_rates(model)

    def test_jukes_cantor(self):
        model = msprime.JukesCantor()
        self.verify_model(model)
        self.verify_mutation_rates(model)

    def test_arbitrary_model(self):
        model = msprime.MutationModel(
            alleles=[b"abc", b"", b"x"],
            root_distribution=[0.8, 0.0, 0.2],
            transition_matrix=[[0.2, 0.4, 0.4], [0.1, 0.2, 0.7], [0.5, 0.3, 0.2]],
        )
        self.verify_model(model)
        self.verify_mutation_rates(model)


class TestPythonMutationGenerator(unittest.TestCase):
    """
    Tests for the python implementation of the mutation generator
    (but, note it is compared to the C version extensively above).
    """

    def verify(self, ts, **kwargs):
        rates = [0, 0.01, 0.2]
        discretes = [True, False]
        keeps = [True, False]
        for rate, keep, discrete in itertools.product(rates, keeps, discretes):
            ts1 = msprime.mutate(ts, rate=rate, keep=keep, discrete=discrete, **kwargs)
            ts2 = py_mutate(ts, rate=rate, keep=keep, discrete=discrete, **kwargs)
            tables1 = ts1.dump_tables()
            tables2 = ts2.dump_tables()
            tables1.provenances.clear()
            tables2.provenances.clear()
            self.assertEqual(tables1, tables2)

    def test_single_tree_no_mutations(self):
        ts = msprime.simulate(10, length=100, random_seed=1234)
        self.verify(ts, random_seed=234)

    def test_single_tree_mutations(self):
        ts = msprime.simulate(10, length=100, mutation_rate=0.1, random_seed=1234)
        self.assertGreater(ts.num_sites, 0)
        self.verify(ts, random_seed=34)

    def test_many_trees_no_mutations(self):
        ts = msprime.simulate(10, length=100, recombination_rate=0.1, random_seed=123)
        self.assertGreater(ts.num_trees, 1)
        self.verify(ts, random_seed=789)

    def test_many_trees_mutations(self):
        ts = msprime.simulate(
            10, length=100, mutation_rate=0.1, recombination_rate=2, random_seed=123
        )
        self.assertGreater(ts.num_trees, 1)
        self.assertGreater(ts.num_sites, 1)
        self.verify(ts, random_seed=789)


####################################################
# Python implementation of lib/mutgen.c algorithms #
####################################################


def py_mutate(ts, rate=None, random_seed=None, model=None, keep=False, discrete=False):
    """
    Same interface as mutations.mutate() and should provide identical results.
    """
    if rate is None:
        rate = 0
    if model is None:
        model = msprime.BinaryMutations()
    tables = ts.dump_tables()
    mutmap = msprime.MutationMap([0, ts.sequence_length], [rate, 0])
    mutgen = PythonMutationGenerator(mutmap, model)
    return mutgen.generate(tables, random_seed, keep=keep, discrete=discrete)


@attr.s
class Site:
    position = attr.ib()
    ancestral_state = attr.ib()
    metadata = attr.ib()
    mutations = attr.ib()
    new = attr.ib()

    def __str__(self):
        s = f"Position: {self.position}\t{self.ancestral_state}"
        s += f"\t{self.metadata}\t{self.new}\n"
        for mut in self.mutations:
            s += mut.__str__()
        return s

    def add_mutation(
        self,
        node,
        time,
        new,
        derived_state=None,
        metadata=b"",
        id=tskit.NULL,  # noqa: A002
    ):
        mutation = Mutation(
            node=node,
            derived_state=derived_state,
            parent=None,
            metadata=metadata,
            time=time,
            new=new,
            keep=True,
            id=id,
        )
        self.mutations.append(mutation)


@attr.s
class Mutation:
    node = attr.ib()
    derived_state = attr.ib()
    parent = attr.ib()
    metadata = attr.ib()
    time = attr.ib()
    new = attr.ib()
    keep = attr.ib()
    id = attr.ib()  # noqa: A003

    def __str__(self):
        if self.parent is None:
            parent_id = None
        else:
            parent_id = self.parent.id
        s = f"\t{self.id}\t\tnode: {self.node}\tparent: {parent_id}"
        s += f"\ttime: {self.time}\t{self.derived_state}\t{self.metadata}"
        s += f"\t(new: {self.new})\tkeep: {self.keep}"
        return s


def cmp_mutation(a, b):
    # Sort mutations by decreasing time and increasing parent,
    # but preserving order of any kept mutations (assumed to be
    # in order already). Kept mutations are given an id that is
    # their order in the initial tables, and new mutations have id -1.
    out = a.id * (not a.new) - b.id * (not b.new)
    if out == 0:
        out = (b.time > a.time) - (a.time > b.time)
    return out


class PythonMutationGenerator:
    """
    Python implementation of the C code which should produce precisely the same
    results.
    """

    def __init__(self, rate_map, model):
        """
        Defaults to all 0->1 mutations.
        """
        self.rate_map = rate_map
        # for compatability with the C code we provide alleles as bytes,
        # but we want them as strings here for simplicity.
        self.alleles = [allele.decode() for allele in model.alleles]
        self.transition_matrix = model.transition_matrix
        self.root_distribution = model.root_distribution
        self.sites = {}

    def print_state(self):
        print("alleles", self.alleles)
        positions = sorted(self.sites.keys())
        for pos in positions:
            print(self.sites[pos])

    def add_site(self, position, new, ancestral_state=None, metadata=b""):
        assert position not in self.sites
        site = Site(
            position=position,
            ancestral_state=ancestral_state,
            metadata=metadata,
            mutations=[],
            new=new,
        )
        self.sites[position] = site
        return site

    def initialise_sites(self, tables):
        nodes = tables.nodes
        mutation_rows = iter(tables.mutations)
        mutation_row = next(mutation_rows, None)
        j = 0
        for site_id, site_row in enumerate(tables.sites):
            site = self.add_site(
                position=site_row.position,
                new=False,
                ancestral_state=site_row.ancestral_state,
                metadata=site_row.metadata,
            )
            while mutation_row is not None and mutation_row.site == site_id:
                site.add_mutation(
                    node=mutation_row.node,
                    time=nodes.time[mutation_row.node],
                    new=False,
                    derived_state=mutation_row.derived_state,
                    metadata=mutation_row.metadata,
                    id=j,
                )
                j += 1
                mutation_row = next(mutation_rows, None)

    def populate_tables(self, tables):
        positions = sorted(self.sites.keys())
        site_id = 0
        for pos in positions:
            site = self.sites[pos]
            num_mutations = 0
            for mutation in site.mutations:
                if mutation.keep:
                    if mutation.parent is None:
                        parent_id = tskit.NULL
                    else:
                        parent_id = mutation.parent.id
                        assert parent_id >= 0
                    mutation_id = tables.mutations.add_row(
                        site_id,
                        mutation.node,
                        mutation.derived_state,
                        parent=parent_id,
                        metadata=mutation.metadata,
                    )
                    assert mutation_id > parent_id
                    mutation.id = mutation_id
                    num_mutations += 1

            if (not site.new) or num_mutations > 0:
                sid = tables.sites.add_row(
                    site.position, site.ancestral_state, site.metadata
                )
                assert sid == site_id
                site_id += 1

    def place_mutations(self, tables, discrete=False):
        # Insert a sentinel into the map for convenience.
        map_position = np.hstack([self.rate_map.position, [tables.sequence_length]])
        node_times = tables.nodes.time
        for edge in tables.edges:
            branch_start = node_times[edge.child]
            branch_end = node_times[edge.parent]
            branch_length = branch_end - branch_start
            index = np.searchsorted(map_position, edge.left)
            if map_position[index] > edge.left:
                index -= 1
            left = edge.left
            right = 0
            while right != edge.right:
                right = min(edge.right, map_position[index + 1])
                site_left = np.ceil(left) if discrete else left
                site_right = np.ceil(right) if discrete else right
                assert site_left <= site_right
                assert map_position[index] <= left
                assert right <= map_position[index + 1]
                assert right <= edge.right
                # Generate the mutations.
                rate = self.rate_map.rate[index]
                mu = rate * (site_right - site_left) * branch_length
                for _ in range(self.rng.poisson(mu)):
                    position = self.rng.flat(site_left, site_right)
                    if discrete:
                        position = np.floor(position)
                    assert edge.left <= position
                    assert position < edge.right
                    if position not in self.sites:
                        self.add_site(position=position, new=True)
                    site = self.sites[position]
                    time = self.rng.flat(branch_start, branch_end)
                    site.add_mutation(node=edge.child, time=time, new=True)
                index += 1
                left = right

    def choose_index(self, distribution):
        u = self.rng.flat(0, 1)
        j = 0
        while u > distribution[j]:
            u -= distribution[j]
            j += 1
        return j

    def pick_allele(self):
        return self.choose_index(self.root_distribution)

    def transition_allele(self, current_allele):
        return self.choose_index(self.transition_matrix[current_allele])

    def choose_alleles(self, tree_parent, site, mutation_id_offset):
        if site.new:
            ai = self.pick_allele()
            site.ancestral_state = self.alleles[ai]
        else:
            ai = self.alleles.index(site.ancestral_state)
        # sort mutations by (increasing id if both are not null,
        #  decreasing time, increasing insertion order)
        site.mutations.sort(key=functools.cmp_to_key(cmp_mutation))
        bottom_mutation = {}
        for mut in site.mutations:
            # Traverse up the tree to find the parent mutation
            # bottom_mutation[u] is the index in mutations of the most recent
            #    mutation seen on the edge above u so far, if any
            u = mut.node
            while u != tskit.NULL and u not in bottom_mutation:
                u = tree_parent[u]
            if u == tskit.NULL:
                # parent is the ancestral state
                pa = ai
                assert mut.parent is None
            else:
                assert u in bottom_mutation
                parent_mut = bottom_mutation[u]
                mut.parent = parent_mut
                assert mut.time <= parent_mut.time, "Parent after child mutation."
                if mut.time > parent_mut.time or (parent_mut.new and not mut.new):
                    raise ValueError(
                        "Generated mutation appears above "
                        "an existing mutation: cannot apply"
                        "finite sites mutations to a earlier"
                        "time period than where they already exist."
                    )
                if mut.new:
                    pa = self.alleles.index(parent_mut.derived_state)

            if mut.new:
                da = self.transition_allele(pa)
                if da == pa:
                    mut.keep = False
                else:
                    mut.derived_state = self.alleles[da]
            if mut.keep:
                bottom_mutation[mut.node] = mut

    def apply_mutations(self, tables):
        ts = tables.tree_sequence()
        positions = sorted(self.sites.keys())
        j = 0
        mutation_id_offset = 0
        tree_parent = np.repeat(tskit.NULL, tables.nodes.num_rows)
        for (_tree_left, tree_right), edges_out, edges_in in ts.edge_diffs():
            for edge in edges_out:
                tree_parent[edge.child] = tskit.NULL
            for edge in edges_in:
                tree_parent[edge.child] = edge.parent
            while j < len(positions) and positions[j] < tree_right:
                site = self.sites[positions[j]]
                self.choose_alleles(tree_parent, site, mutation_id_offset)
                num_mutations = len(site.mutations)
                mutation_id_offset += num_mutations
                j += 1

    def generate(self, tables, seed, keep=False, discrete=False):
        self.rng = _msprime.RandomGenerator(seed)
        if keep:
            self.initialise_sites(tables)
        tables.sites.clear()
        tables.mutations.clear()
        self.place_mutations(tables, discrete=discrete)
        self.apply_mutations(tables)
        self.populate_tables(tables)
        self.record_provenance(tables, seed, keep, discrete)
        return tables.tree_sequence()

    def record_provenance(self, tables, seed, keep, discrete):
        parameters = {
            "command": "mutate",
            "random_seed": seed,
            "keep": keep,
            "discrete": discrete,
            "rate_map": None,  # not working??
            "alleles": self.alleles,
            "transition_matrix": self.transition_matrix,
            "root_distribution": self.root_distribution,
        }
        provenance_dict = msprime.provenance.get_provenance_dict(parameters)
        encoded_provenance = msprime.provenance.json_encode_provenance(provenance_dict)
        tables.provenances.add_row(encoded_provenance)
