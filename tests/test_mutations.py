#
# Copyright (C) 2018-2021 University of Oxford
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
from __future__ import annotations

import dataclasses
import functools
import itertools
import json
import struct
from typing import Any
from typing import List

import numpy as np
import pytest
import scipy.stats as stats
import tskit

import msprime
import tests.wright_fisher as wf
from msprime import _msprime
from tests import tsutil


class TestDefaults:
    def test_sim_mutations(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        mts = msprime.sim_mutations(ts, rate=1, random_seed=1)
        # Finite sites, JC69 by default.
        assert mts.sequence_length == 1
        assert mts.num_sites == 1
        # High rate, so should have > 1 mutation at this site.
        assert mts.num_mutations > 1
        assert mts.site(0).position == 0
        assert mts.site(0).ancestral_state in "ACGT"
        for mutation in mts.mutations():
            assert mutation.derived_state in "ACGT"

    def test_sim_mutations_existing(self):
        ts = msprime.sim_ancestry(10, sequence_length=10, random_seed=1)
        mts = msprime.sim_mutations(ts, rate=0.1, random_seed=1)
        assert mts.num_sites > 1
        assert mts.num_mutations > 1
        # By default trying to put mutations on top of existing ones
        # raises an error.
        with pytest.raises(_msprime.LibraryError):
            msprime.sim_mutations(mts, rate=1, random_seed=1)
        # But it's OK if we set add_ancestral to True.
        extra_mts = msprime.sim_mutations(
            mts, rate=0.1, random_seed=2, add_ancestral=True
        )
        assert extra_mts.num_sites > 1
        assert extra_mts.num_mutations > mts.num_mutations

        # If we turn up the mutation rate too high we get some
        # bad ancestral mutations with silent state transitions.
        with pytest.raises(_msprime.LibraryError, match="silent transition"):
            msprime.sim_mutations(mts, rate=10, random_seed=2, add_ancestral=True)

    def test_mutate(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        mts = msprime.mutate(ts, rate=1, random_seed=1)
        # Continuous genome, 0/1 alleles
        assert mts.num_sites > 1
        assert mts.sequence_length == 1
        assert mts.num_mutations == mts.num_sites
        for site in mts.sites():
            assert site.ancestral_state == "0"
        for mutation in mts.mutations():
            assert mutation.derived_state == "1"

    def test_mutate_existing(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        mts = msprime.mutate(ts, rate=1, random_seed=1)
        # Continuous genome, 0/1 alleles
        assert mts.num_sites > 1
        extra_mts = msprime.mutate(ts, rate=0, random_seed=1)
        # Keep is false by default
        assert extra_mts.num_sites == 0
        extra_mts = msprime.mutate(ts, rate=1, random_seed=1)
        assert extra_mts.equals(mts, ignore_provenance=True)


class TestMutateProvenance:
    """
    Test that we correctly record the provenance for calls to mutate.
    """

    def test_mutation_rate(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        for mutation_rate in [0, 1, 1e-5]:
            mutated = msprime.sim_mutations(ts, mutation_rate)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            assert record["parameters"]["command"] == "sim_mutations"
            assert record["parameters"]["rate"] == mutation_rate
            assert record["parameters"]["random_seed"] >= 0

    def test_start_time(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        for start_time in [0, 1, -1]:
            mutated = msprime.sim_mutations(ts, 1, start_time=start_time)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            assert record["parameters"]["command"] == "sim_mutations"
            assert record["parameters"]["start_time"] == start_time
            assert record["parameters"]["random_seed"] >= 0

    def test_end_time(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        for end_time in [0, 1, 100]:
            mutated = msprime.sim_mutations(ts, 1, start_time=-1, end_time=end_time)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            assert record["parameters"]["command"] == "sim_mutations"
            assert record["parameters"]["start_time"] == -1
            assert record["parameters"]["end_time"] == end_time
            assert record["parameters"]["random_seed"] >= 0

    def test_seed(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        for seed in range(1, 10):
            mutated = msprime.sim_mutations(ts, rate=1, random_seed=seed)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            assert record["parameters"]["command"] == "sim_mutations"
            assert record["parameters"]["rate"] == 1
            assert record["parameters"]["random_seed"] == seed

    def test_keep(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        for keep in [True, False]:
            mutated = msprime.sim_mutations(ts, rate=1, keep=keep)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            assert record["parameters"]["command"] == "sim_mutations"
            assert record["parameters"]["keep"] == keep


class TestMatrixMutationModel:
    def validate_model(self, model):
        num_alleles = len(model.alleles)
        assert num_alleles == len(model.root_distribution)
        assert (num_alleles, num_alleles) == model.transition_matrix.shape
        assert np.allclose(sum(model.root_distribution), 1.0)
        for j in range(num_alleles):
            assert np.allclose(sum(model.transition_matrix[j]), 1.0)
        s = str(model)
        assert isinstance(s, str)
        assert s.startswith("Mutation model with alleles")

    def validate_stationary(self, model):
        # for most models, the root distribution should be the
        # stationary distribution of the transition matrix
        x = np.matmul(model.root_distribution, model.transition_matrix)
        assert np.allclose(list(x), list(model.root_distribution))

    def validate_asdict(self, model):
        m = msprime.MatrixMutationModel(**model.asdict())
        assert model == m

    def test_bad_alleles(self):
        for alleles, err in [
            (0, TypeError),
            ("xyz", TypeError),
            ([b"0", "1"], TypeError),
            ([1, 2, 3], TypeError),
            (["a"], _msprime.LibraryError),
        ]:
            with pytest.raises(err):
                msprime.MatrixMutationModel(
                    alleles=alleles,
                    root_distribution=[1] + [0] * (len(alleles) - 1),
                    transition_matrix=[[1 / len(alleles)] * len(alleles)]
                    * len(alleles),
                )

    def test_bad_root_distribution(self):
        alleles = ["a", "b", "c"]
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
            with pytest.raises(err):
                msprime.MatrixMutationModel(
                    alleles=alleles,
                    root_distribution=root_distribution,
                    transition_matrix=transition_matrix,
                )

    def test_bad_transition_matrix(self):
        alleles = ["a", "b", "c"]
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
            with pytest.raises(err):
                msprime.MatrixMutationModel(
                    alleles=alleles,
                    root_distribution=root_distribution,
                    transition_matrix=transition_matrix,
                )

    def test_binary(self):
        model = msprime.BinaryMutationModel()
        self.validate_model(model)

    def test_jukes_cantor(self):
        model = msprime.JC69MutationModel()
        self.validate_model(model)
        self.validate_stationary(model)

    def test_HKY_default(self):
        model = msprime.HKYMutationModel(0.75)
        self.validate_model(model)
        self.validate_stationary(model)

    def test_F84_default(self):
        model = msprime.F84MutationModel(0.75)
        self.validate_model(model)
        self.validate_stationary(model)

    def test_GTR_default(self):
        model = msprime.GTRMutationModel([1 / 6] * 6)
        self.validate_model(model)
        self.validate_stationary(model)

    def test_PAM(self):
        model = msprime.PAMMutationModel()
        self.validate_model(model)
        self.validate_stationary(model)

    def test_BLOSUM62(self):
        model = msprime.BLOSUM62MutationModel()
        self.validate_model(model)
        self.validate_stationary(model)

    def test_new_model(self):
        alleles = ["Alligator", "Camel", "Goat"]
        root_distribution = [0.5, 0.25, 0.25]
        transition_matrix = [[0.0, 0.5, 0.5], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        model = msprime.MatrixMutationModel(
            alleles, root_distribution, transition_matrix
        )
        assert model.alleles == alleles
        assert list(model.root_distribution) == root_distribution
        for a, b in zip(model.transition_matrix, transition_matrix):
            assert list(a) == b
        self.validate_model(model)
        self.validate_stationary(model)


class TestMutateRateMap:
    def test_no_muts(self):
        ts = msprime.sim_ancestry(10, sequence_length=10, random_seed=2)
        pos = [0, 3, 7, 10]
        rate = [1.0, 0.0, 1.0]
        mut_map = msprime.RateMap(position=pos, rate=rate)
        ts = msprime.sim_mutations(
            ts, rate=mut_map, random_seed=1, discrete_genome=True
        )
        site_pos = ts.tables.sites.position
        assert np.all(site_pos == [0, 1, 2, 7, 8, 9])

    def test_many_edges(self):
        ts = msprime.sim_ancestry(
            10, sequence_length=10, random_seed=6, recombination_rate=10
        )
        pos = [0, 3, 7, 10]
        rate = [1.0, 0.0, 1.0]
        mut_map = msprime.RateMap(position=pos, rate=rate)
        ts = msprime.sim_mutations(
            ts, rate=mut_map, random_seed=1, discrete_genome=True
        )
        site_pos = ts.tables.sites.position
        assert np.all(site_pos == [0, 1, 2, 7, 8, 9])

    def test_noninteger_positions_discrete(self):
        ts = msprime.sim_ancestry(
            10, sequence_length=10, random_seed=2, recombination_rate=1.0
        )
        pos = [0, 2.9, 4.3, 4.9, 9.5, 10]
        rate = [1.0, 0.0, 10.0, 0.0, 10.0]
        mut_map = msprime.RateMap(position=pos, rate=rate)
        ts = msprime.sim_mutations(
            ts, rate=mut_map, random_seed=2, discrete_genome=True
        )
        site_pos = ts.tables.sites.position
        assert np.all(site_pos == [0, 1, 2])

    def test_noninteger_positions_continuous(self):
        ts = msprime.sim_ancestry(
            10,
            sequence_length=10,
            random_seed=4,
            discrete_genome=False,
            recombination_rate=1.0,
        )
        pos = [0, 3, 4.3, 4.9, 10]
        rate = [1.0, 0.0, 1.0, 0.0]
        mut_map = msprime.RateMap(position=pos, rate=rate)
        ts = msprime.sim_mutations(
            ts, rate=mut_map, random_seed=2, discrete_genome=False
        )
        site_pos = ts.tables.sites.position
        for left, right, rate in zip(pos[:-1], pos[1:], rate):
            n = np.sum(np.logical_and(site_pos >= left, site_pos < right))
            assert (n > 0) == (rate > 0)

    def test_rate_map_applies_only_to_integers(self):
        ts = msprime.sim_ancestry(
            10,
            sequence_length=10,
            random_seed=3,
            discrete_genome=True,
            recombination_rate=1.0,
        )
        pos = [0, 0.5, 1.0, 4.5, 5.0, 9.5, 10]
        rate = [0, 10, 0, 10, 0, 10]
        mut_map = msprime.RateMap(position=pos, rate=rate)
        ts = msprime.sim_mutations(
            ts, rate=mut_map, random_seed=2, discrete_genome=True
        )
        assert ts.num_mutations == 0


class MutateMixin:
    def verify_topology(self, source, dest):
        assert source.nodes == dest.nodes
        assert source.edges == dest.edges
        assert source.migrations == dest.migrations
        assert source.populations == dest.populations

    def verify_provenance(self, source, dest):
        before = list(source.provenances)
        after = list(dest.provenances)
        assert len(before) + 1 == len(after)
        assert before == after[:-1]

    def verify_binary_alphabet(self, ts):
        assert ts.num_sites > 0
        assert all(site.ancestral_state == "0" for site in ts.sites())
        assert all(mutation.derived_state == "1" for mutation in ts.mutations())

    def verify_nucleotides_alphabet(self, ts):
        nucleotides = "ACGT"
        assert ts.num_sites > 0
        assert all(site.ancestral_state in nucleotides for site in ts.sites())
        assert all(mutation.derived_state in nucleotides for mutation in ts.mutations())
        for site in ts.sites():
            assert site.ancestral_state != site.mutations[0].derived_state


class TestMutate(MutateMixin):
    def test_default_alphabet(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        mutated = msprime.mutate(ts, rate=1, random_seed=2)
        self.verify_binary_alphabet(mutated)

    def test_continuous_genome(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        mutated = msprime.mutate(ts, rate=1, random_seed=2)
        position = mutated.tables.sites.position
        assert len(position) > 0
        assert np.all(np.floor(position) != position)

    def test_binary_alphabet(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        mutated = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.BINARY)
        )
        self.verify_binary_alphabet(mutated)

    def test_nucleotide_alphabet(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        mutated = msprime.mutate(
            ts, rate=1, random_seed=2, model=msprime.InfiniteSites(msprime.NUCLEOTIDES)
        )
        self.verify_nucleotides_alphabet(mutated)

    def test_bad_alphabet(self):
        ts = msprime.sim_ancestry(2, random_seed=2)
        for bad_alphabet in [-1, 3, "binary"]:
            with pytest.raises(ValueError):
                msprime.mutate(
                    ts, rate=1, random_seed=2, model=msprime.InfiniteSites(bad_alphabet)
                )

    def test_deprecated_identical_seed_alphabets(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        binary = msprime.mutate(ts, rate=1, random_seed=2)
        nucs = msprime.mutate(
            ts,
            rate=1,
            random_seed=2,
            model=msprime.InfiniteSites(msprime.NUCLEOTIDES),
        )
        assert binary.num_sites > 0
        assert binary.num_mutations > 0
        assert binary.num_sites == nucs.num_sites
        assert binary.num_mutations == nucs.num_mutations
        for s1, s2 in zip(binary.sites(), nucs.sites()):
            assert s1.position == s2.position
            assert s1.mutations[0].node == s2.mutations[0].node

    def test_new_models_unsupported(self):
        ts = msprime.sim_ancestry(2, random_seed=2)
        for model in [
            msprime.PAMMutationModel(),
            "jc69",
            msprime.BinaryMutationModel(),
        ]:
            with pytest.raises(ValueError):
                msprime.mutate(ts, rate=1, random_seed=2, model=model)


class TestSimMutations(MutateMixin):
    def test_unicode_alleles(self):
        alleles = ["🎄🌳", "💩" * 5]
        binary = msprime.BinaryMutationModel()
        model = msprime.MatrixMutationModel(
            alleles, binary.root_distribution, binary.transition_matrix
        )
        ts = msprime.sim_ancestry(8, random_seed=2)
        mts = msprime.sim_mutations(
            ts, rate=2, random_seed=1, model=model, discrete_genome=False
        )
        assert mts.num_sites > 0
        for tree in mts.trees():
            for site in tree.sites():
                assert site.ancestral_state == alleles[0]
                for mutation in site.mutations:
                    assert mutation.derived_state == alleles[1]

    def test_zero_mutation_rate(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        mutated = msprime.sim_mutations(ts, 0)
        t1 = ts.dump_tables()
        t2 = mutated.dump_tables()
        self.verify_topology(t1, t2)
        self.verify_provenance(t1, t2)
        assert t1.sites == t2.sites
        assert t1.mutations == t2.mutations

    def test_populations(self):
        demography = msprime.Demography.island_model([1, 1], migration_rate=1)
        ts = msprime.sim_ancestry(
            {0: 10, 1: 10},
            demography=demography,
            record_migrations=True,
            random_seed=1,
        )
        mutated = msprime.sim_mutations(ts, 0)
        t1 = ts.dump_tables()
        assert len(t1.populations) == 2
        assert len(t1.migrations) > 0
        t2 = mutated.dump_tables()
        self.verify_topology(t1, t2)
        self.verify_provenance(t1, t2)
        assert t1.sites == t2.sites
        assert t1.mutations == t2.mutations

    def test_mutation_overwrite(self):
        ts = msprime.simulate(10, mutation_rate=5, random_seed=2)
        assert ts.num_sites > 0
        assert ts.num_mutations > 0
        mutated = msprime.sim_mutations(ts, 0, keep=False)
        t1 = ts.dump_tables()
        assert len(t1.sites) == ts.num_sites
        t2 = mutated.dump_tables()
        self.verify_topology(t1, t2)
        self.verify_provenance(t1, t2)
        assert len(t2.sites) == 0
        assert len(t2.mutations) == 0

    def test_bad_rates(self):
        ts = msprime.sim_ancestry(2, random_seed=2)
        for bad_type in [{}, [], ts]:
            with pytest.raises(TypeError):
                msprime.sim_mutations(ts, rate=bad_type)
        for bad_rate in ["abc", "xxx"]:
            with pytest.raises(ValueError):
                msprime.sim_mutations(ts, bad_rate)
        for bad_rate in [-1, -1e-6, -1e7]:
            with pytest.raises(ValueError):
                msprime.sim_mutations(ts, bad_rate)

    def test_bad_models(self):
        ts = msprime.sim_ancestry(2, random_seed=2)
        for bad_type in [{}, True, 123]:
            with pytest.raises(TypeError):
                msprime.sim_mutations(ts, rate=0, model=bad_type)
        for bad_name in ["", "coffee"]:
            with pytest.raises(ValueError):
                msprime.sim_mutations(ts, rate=0, model=bad_name)

    def test_bad_tree_sequence(self):
        for bad_type in [None, {}, "sdrf"]:
            with pytest.raises(ValueError):
                msprime.sim_mutations(bad_type)

    def test_default_seeds(self):
        ts = msprime.sim_ancestry(20, random_seed=2)
        seeds = []
        for _ in range(10):
            mutated = msprime.sim_mutations(ts, 0)
            record = json.loads(mutated.provenance(mutated.num_provenances - 1).record)
            seeds.append(record["parameters"]["random_seed"])
        assert len(seeds) == len(set(seeds))

    def test_default_discrete_genome(self):
        ts = msprime.sim_ancestry(20, sequence_length=10, random_seed=2)
        mutated = msprime.sim_mutations(ts, 2, random_seed=1234)
        position = mutated.tables.sites.position
        assert len(position) > 0
        assert np.all(np.floor(position) == position)

    def test_identical_seed(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        mutated = [
            msprime.sim_mutations(ts, rate=1, random_seed=2) for _ in range(1, 10)
        ]
        assert mutated[0].num_sites > 0
        assert mutated[0].num_mutations > 0
        tables = [other_ts.dump_tables() for other_ts in mutated]
        assert all(tables[0].sites == t.sites for t in tables[1:])
        assert all(tables[0].mutations == t.mutations for t in tables[1:])

    def test_numpy_inputs(self):
        ts = msprime.sim_ancestry(10, sequence_length=10, random_seed=2)
        values = np.array([0, 1, 2.5])
        mutated1 = msprime.sim_mutations(
            ts,
            rate=values[2],
            random_seed=values[1],
            start_time=values[0],
            end_time=values[2],
        )
        mutated2 = msprime.sim_mutations(
            ts,
            rate=float(values[2]),
            random_seed=int(values[1]),
            start_time=float(values[0]),
            end_time=float(values[2]),
        )
        assert mutated1.equals(mutated2, ignore_timestamps=True)

    def test_boolean_inputs(self):
        ts = msprime.sim_ancestry(10, sequence_length=10, random_seed=2)
        # Falsey values aren't allowed
        with pytest.raises(TypeError):
            msprime.sim_mutations(ts, rate=10, discrete_genome=[])
        with pytest.raises(TypeError):
            msprime.sim_mutations(ts, rate=10, keep=[])
        with pytest.raises(TypeError):
            msprime.sim_mutations(ts, rate=10, add_ancestral=[])


class TestFiniteSites(MutateMixin):
    def verify_binary_alphabet(self, ts):
        binary = "01"
        assert ts.num_sites > 0
        assert all(site.ancestral_state == "0" for site in ts.sites())
        assert all(mutation.derived_state in binary for mutation in ts.mutations())

    def verify_mutations(self, ts, discrete_genome, check_probs=False, model=None):
        # Verify that mutation.parents are correct
        # and that mutations actually involve a change of state
        # that, if (not keep), has positive probability under the model
        if check_probs:
            assert model is not None
            alleles = model.alleles
        for site in ts.sites():
            if not discrete_genome:
                assert len(site.mutations) == 1
            t = ts.at(site.position)
            parents = {}
            for mut in site.mutations:
                n = mut.node
                assert mut.time >= ts.node(n).time
                assert mut.time < ts.node(t.parent(n)).time
                while n != tskit.NULL and n not in parents:
                    n = t.parent(n)
                if n == tskit.NULL:
                    assert mut.parent == tskit.NULL
                    pa = site.ancestral_state
                else:
                    assert mut.parent == parents[n].id
                    assert mut.time < parents[n].time
                    pa = parents[n].derived_state
                assert mut.derived_state != pa
                if check_probs:
                    assert pa in alleles
                    assert mut.derived_state in alleles
                    i = alleles.index(pa)
                    j = alleles.index(mut.derived_state)
                    assert model.transition_matrix[i][j] > 0
                # do this last since parents should be listed before children
                parents[mut.node] = mut

    def mutate(self, ts, model, rate=1.0, keep=False, discrete_genome=False, seed=42):
        lib_ts = msprime.sim_mutations(
            ts,
            rate=rate,
            random_seed=seed,
            model=model,
            keep=keep,
            discrete_genome=discrete_genome,
        )
        self.verify_mutations(
            lib_ts, discrete_genome, check_probs=not keep, model=model
        )
        py_ts = py_sim_mutations(
            ts,
            rate=rate,
            random_seed=seed,
            model=model,
            keep=keep,
            discrete_genome=discrete_genome,
        )
        lib_tables = lib_ts.dump_tables()
        lib_tables.provenances.clear()
        py_tables = py_ts.dump_tables()
        py_tables.provenances.clear()
        assert lib_tables == py_tables
        return lib_ts

    def mutate_binary(
        self,
        ts,
        rate=1.0,
        keep=False,
        discrete_genome=False,
        transition_matrix=None,
        root_distribution=None,
    ):
        alleles = ["0", "1"]
        if transition_matrix is None:
            transition_matrix = [[0.0, 1.0], [1.0, 0.0]]
        if root_distribution is None:
            root_distribution = [1.0, 0.0]
        model = msprime.MatrixMutationModel(
            alleles, root_distribution, transition_matrix
        )
        return self.mutate(
            ts, model, rate=rate, keep=keep, discrete_genome=discrete_genome
        )

    def mutate_nucleotides(
        self,
        ts,
        rate=1.0,
        keep=False,
        discrete_genome=False,
        transition_matrix=None,
        root_distribution=None,
    ):
        alleles = ["A", "C", "G", "T"]
        if transition_matrix is None:
            transition_matrix = [[0.25] * 4] * 4
        if root_distribution is None:
            root_distribution = [0.25] * 4
        model = msprime.MatrixMutationModel(
            alleles, root_distribution, transition_matrix
        )
        return self.mutate(
            ts, model, rate=rate, keep=keep, discrete_genome=discrete_genome
        )

    def test_alleles_binary(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
        mutated = self.mutate_binary(ts)
        self.verify_binary_alphabet(mutated)

    def test_alleles_nucleotides(self):
        ts = msprime.sim_ancestry(10, random_seed=1)
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
                assert t1.sites == t2.sites
                assert t1.mutations == t2.mutations
            else:
                assert t2.sites.num_rows == 0

    def test_bad_mutate_order(self):
        ts = msprime.sim_ancestry(
            10, random_seed=1, recombination_rate=1, sequence_length=10
        )
        mutated = msprime.sim_mutations(
            ts, 3, random_seed=5, start_time=0.0, end_time=0.5, discrete_genome=True
        )
        with pytest.raises(_msprime.LibraryError):
            msprime.sim_mutations(
                mutated,
                3,
                random_seed=6,
                start_time=0.5,
                end_time=1.0,
                discrete_genome=True,
                keep=True,
            )

    def test_one_way_mutation(self):
        for discrete_genome in (True, False):
            ts = msprime.sim_ancestry(
                10, random_seed=1, recombination_rate=1.0, sequence_length=10
            )
            mut_matrix = [[0.0, 1.0], [0.0, 1.0]]
            mutated = self.mutate_binary(
                ts,
                rate=1.0,
                transition_matrix=mut_matrix,
                root_distribution=[1.0, 0.0],
                discrete_genome=discrete_genome,
            )
            assert mutated.num_sites > 0
            if discrete_genome:
                assert max([len(s.mutations) for s in mutated.sites()]) > 1
            for site in mutated.sites():
                assert site.ancestral_state == "0"
                assert len(site.mutations) >= 1
                for mut in site.mutations:
                    assert mut.derived_state == "1"

    def test_flip_flop_mutation(self):
        nucleotides = "ACGT"
        for discrete_genome in (True, False):
            ts = msprime.sim_ancestry(
                10, random_seed=1, recombination_rate=1.0, sequence_length=10
            )
            mut_matrix = [
                [0.0, 0.0, 0.5, 0.5],
                [0.0, 0.0, 0.5, 0.5],
                [0.5, 0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0, 0.0],
            ]
            mutated = self.mutate_nucleotides(
                ts,
                rate=5.0,
                transition_matrix=mut_matrix,
                discrete_genome=discrete_genome,
            )
            assert mutated.num_sites > 0
            if discrete_genome:
                assert max([len(s.mutations) for s in mutated.sites()]) > 1
            for mut in ts.mutations():
                assert mut.derived_state in nucleotides
                if mut.parent == -1:
                    pa = ts.site(mut.site).ancestral_state
                else:
                    pa = ts.mutation(mut.parent).derived_state
                if pa in "AC":
                    assert mut.derived_state in "GT"
                else:
                    assert mut.derived_state in "AC"

    def test_do_nothing_mutations(self):
        for discrete_genome in (True, False):
            ts = msprime.sim_ancestry(10, random_seed=1, sequence_length=10)
            mut_matrix = [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
            mutated = self.mutate_nucleotides(
                ts,
                rate=5.0,
                transition_matrix=mut_matrix,
                discrete_genome=discrete_genome,
            )
            assert mutated.num_sites == 0
            t1 = ts.dump_tables()
            t2 = mutated.dump_tables()
            self.verify_topology(t1, t2)
            self.verify_provenance(t1, t2)
            assert t1.sites == t2.sites
            assert t1.mutations == t2.mutations

    def test_uniform_mutations(self):
        nucleotides = "ACGT"
        for discrete_genome in (True, False):
            ts = msprime.sim_ancestry(
                10, random_seed=1, recombination_rate=1.0, sequence_length=10
            )
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
                discrete_genome=discrete_genome,
            )
            assert mutated.num_sites > 0
            if discrete_genome:
                assert max([len(s.mutations) for s in mutated.sites()]) > 1
            num_nucs = {a: 0 for a in nucleotides}
            for site in mutated.sites():
                assert site.ancestral_state == "A"
            for mut in mutated.mutations():
                num_nucs[mut.derived_state] += 1
            assert num_nucs["A"] == 0
            assert num_nucs["C"] > 0
            assert num_nucs["G"] > 0
            assert num_nucs["T"] > 0

    def test_circular_mutations(self):
        nucleotides = "ACGT"
        for discrete_genome in (True, False):
            ts = msprime.sim_ancestry(
                10, random_seed=1, recombination_rate=1.0, sequence_length=10
            )
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
                discrete_genome=discrete_genome,
            )
            assert mutated.num_sites > 0
            if discrete_genome:
                assert max([len(s.mutations) for s in mutated.sites()]) > 1
            for site in mutated.sites():
                assert site.ancestral_state == "A"
                assert len(site.mutations) > 0
            for mut in mutated.mutations():
                if mut.parent == tskit.NULL:
                    assert mut.derived_state == "C"
                else:
                    pmut = mutated.mutation(mut.parent)
                    assert (
                        nucleotides.index(pmut.derived_state) + 1
                    ) % 4 == nucleotides.index(mut.derived_state)

    def test_integer_sites(self):
        ts = msprime.sim_ancestry(
            10, random_seed=5, sequence_length=10, recombination_rate=10.0
        )
        mutated = self.mutate_binary(ts, rate=5.0, discrete_genome=True)
        assert mutated.site(0).position == 0.0
        for site in mutated.sites():
            assert site.position == np.floor(site.position)
            assert site.position < mutated.sequence_length


class TestInterval:
    """
    Tests on the start_time and end_time parameters.
    """

    def test_errors(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        for start, end in [(-2, -3), (1, 0), (1e6, 1e5)]:
            with pytest.raises(ValueError):
                msprime.mutate(ts, start_time=start, end_time=end)

    def test_stick_tree(self):
        tables = tskit.TableCollection(1.0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
        tables.nodes.add_row(flags=0, time=1)
        tables.nodes.add_row(flags=0, time=2)
        tables.edges.add_row(0, 1, 1, 0)
        tables.edges.add_row(0, 1, 2, 1)
        ts = tables.tree_sequence()

        tsm = msprime.mutate(ts, rate=100, end_time=1, random_seed=1)
        assert tsm.num_sites > 0
        assert all(mut.node == 0 for mut in ts.mutations())

        tsm = msprime.mutate(ts, rate=100, start_time=0, end_time=1, random_seed=1)
        assert tsm.num_sites > 0
        assert all(mut.node == 0 for mut in ts.mutations())

        tsm = msprime.mutate(ts, rate=100, start_time=0.5, end_time=1, random_seed=1)
        assert tsm.num_sites > 0
        assert all(mut.node == 0 for mut in ts.mutations())

        tsm = msprime.mutate(ts, rate=100, start_time=1, random_seed=1)
        assert tsm.num_sites > 0
        assert all(mut.node == 1 for mut in ts.mutations())

        tsm = msprime.mutate(ts, rate=100, start_time=1, end_time=2, random_seed=1)
        assert tsm.num_sites > 0
        assert all(mut.node == 1 for mut in ts.mutations())

        tsm = msprime.mutate(ts, rate=100, start_time=1.5, end_time=2, random_seed=1)
        assert tsm.num_sites > 0
        assert all(mut.node == 0 for mut in ts.mutations())

    def verify_mutations(self, ts, start, end):
        assert ts.num_sites > 0
        time = ts.tables.nodes.time
        if start is None:
            start = np.min(time) - 1
        if end is None:
            end = np.min(time) + 1
        for tree in ts.trees():
            for site in tree.sites():
                for mutation in site.mutations:
                    top_node = tree.parent(mutation.node)
                    assert time[top_node] > start

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
        ts = msprime.sim_ancestry(20, random_seed=2)
        self.verify(ts)

    def test_coalescent_trees(self):
        ts = msprime.sim_ancestry(
            20, recombination_rate=1, sequence_length=10, random_seed=2
        )
        self.verify(ts)

    def test_wright_fisher_trees(self):
        tables = wf.wf_sim(20, 10, seed=1, deep_history=False)
        tables.sort()
        ts = tables.tree_sequence()
        self.verify(ts, rate=10)
        ts = ts.simplify()
        self.verify(ts, rate=10)

    def test_negative_time(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        tables = ts.dump_tables()
        time = tables.nodes.time
        max_time = np.max(time)
        for offset in [max_time / 2, max_time, 2 * max_time]:
            tables.nodes.set_columns(flags=tables.nodes.flags, time=time - offset)
            ts = tables.tree_sequence()
            self.verify(ts)


class TestKeep:
    """
    Tests for the "keep" functionality in which we append new mutations
    to an existing set.
    """

    def verify(self, ts, rate, random_seed):
        no_keep = msprime.sim_mutations(
            ts,
            rate=rate,
            random_seed=random_seed,
            discrete_genome=False,
            keep=False,
        )
        assert no_keep.num_sites > 0
        keep = msprime.sim_mutations(
            ts, rate=rate, random_seed=random_seed, keep=True, discrete_genome=False
        )
        # Can assume there's no collisions here, very unlikely.
        assert ts.num_sites + no_keep.num_sites == keep.num_sites
        # Mutations are all infinite sites, so must be equal
        assert ts.num_mutations + no_keep.num_mutations == keep.num_mutations
        old = {site.position for site in ts.sites()}
        new = {site.position for site in no_keep.sites()}
        both = {site.position for site in keep.sites()}
        assert old | new == both
        self.verify_sites(ts, keep)

    def verify_sites(self, original, updated):
        site_map = {site.position: site for site in original.sites()}
        found = 0
        inserted_mutations = 0
        for site in updated.sites():
            if site.position in site_map:
                old_site = site_map[site.position]
                assert site.position == old_site.position
                assert len(site.mutations) == len(old_site.mutations)
                for mutation, old_mutation in zip(site.mutations, old_site.mutations):
                    assert mutation.metadata == old_mutation.metadata
                    assert mutation.node == old_mutation.node
                    assert mutation.derived_state == old_mutation.derived_state
                    if old_mutation.parent == tskit.NULL:
                        assert mutation.parent == old_mutation.parent
                    else:
                        assert (
                            mutation.parent == old_mutation.parent + inserted_mutations
                        )
                assert site.ancestral_state == old_site.ancestral_state
                assert site.metadata == old_site.metadata
                found += 1
            else:
                inserted_mutations += len(site.mutations)
        assert found == original.num_sites

    def test_simple_binary(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=2)
        assert ts.num_sites > 0
        self.verify(ts, 1, random_seed=2)

    def test_branch_mutations(self):
        ts = tsutil.insert_branch_mutations(
            msprime.sim_ancestry(
                5, sequence_length=10, recombination_rate=1, random_seed=2
            )
        )
        assert ts.num_sites > 1
        self.verify(ts, 3, random_seed=7)

    def test_multichar_mutations(self):
        ts = tsutil.insert_multichar_mutations(
            msprime.sim_ancestry(
                12, sequence_length=10, recombination_rate=2, random_seed=3
            )
        )
        assert ts.num_sites > 5
        self.verify(ts, 3, random_seed=7)

    def test_random_metadata(self):
        ts = tsutil.add_random_metadata(
            msprime.simulate(12, random_seed=4, mutation_rate=1)
        )
        assert ts.num_sites > 5
        self.verify(ts, 3, random_seed=7)

    def test_no_sites(self):
        ts = msprime.sim_ancestry(12, random_seed=3)
        assert ts.num_sites == 0
        self.verify(ts, 3, random_seed=7)

    def test_site_with_no_mutation(self):
        ts = tsutil.insert_site(msprime.sim_ancestry(12, random_seed=3))
        assert ts.num_sites == 1
        self.verify(ts, 3, random_seed=7)

    def test_same_seeds(self):
        ts = msprime.sim_ancestry(12, random_seed=3)
        assert ts.num_sites == 0
        ts = msprime.sim_mutations(ts, rate=1, random_seed=1, discrete_genome=False)
        updated = msprime.sim_mutations(
            ts, rate=1, random_seed=1, keep=True, discrete_genome=False
        )
        # We should rejection sample away all the sites that we have already
        # generated so we just get another random sample.
        sites = set()
        for site in updated.sites():
            assert site.position not in sites
            sites.add(site.position)
        for site in ts.sites():
            assert site.position in sites

    def test_keep_multichar_muts(self):
        ts = msprime.sim_ancestry(12, random_seed=3)
        ts = msprime.sim_mutations(ts, rate=1, random_seed=1, discrete_genome=False)
        assert ts.num_sites > 2
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        for site in ts.sites():
            tables.sites.add_row(position=site.position, ancestral_state="A" * site.id)
            for mutation in site.mutations:
                tables.mutations.add_row(
                    site=site.id, node=mutation.node, derived_state="T" * site.id
                )
        tables.compute_mutation_times()
        original = tables.tree_sequence()
        updated = msprime.sim_mutations(
            original, rate=1, random_seed=1, keep=True, discrete_genome=False
        )
        self.verify_sites(original, updated)

    def test_keep_metadata(self):
        ts = msprime.sim_ancestry(12, random_seed=3)
        ts = msprime.sim_mutations(ts, rate=1, random_seed=1, discrete_genome=False)
        assert ts.num_sites > 2
        # Set metadata on this ts so that we can be sure we keep the original
        # mutations.
        ts = tsutil.add_random_metadata(ts)
        other = msprime.sim_mutations(
            ts, rate=1, random_seed=1, keep=True, discrete_genome=False
        )
        self.verify_sites(ts, other)

    def test_keep_mutation_parent(self):
        ts = msprime.sim_ancestry(
            6, sequence_length=10, recombination_rate=3, random_seed=3
        )
        ts = tsutil.insert_branch_mutations(ts)
        assert ts.num_sites > 2
        other = msprime.sim_mutations(
            ts, rate=1, random_seed=1, keep=True, discrete_genome=False
        )
        assert other.num_sites > ts.num_sites
        self.verify_sites(ts, other)

    def test_keep_mutation_parent_zero_rate(self):
        ts = msprime.sim_ancestry(
            6, recombination_rate=3, sequence_length=10, random_seed=3
        )
        ts = tsutil.insert_branch_mutations(ts)
        assert ts.num_sites > 2
        other = msprime.sim_mutations(
            ts, rate=0, random_seed=1, keep=True, discrete_genome=False
        )
        t1 = ts.dump_tables()
        t2 = other.dump_tables()
        t1.provenances.clear()
        t2.provenances.clear()
        assert t1 == t2

    def test_keep_unknown_time_muts(self):
        ts = msprime.sim_ancestry(6, random_seed=3)
        ts = msprime.sim_mutations(ts, rate=1, random_seed=1, discrete_genome=False)
        assert ts.num_sites > 2
        tables = ts.dump_tables()
        tables.mutations.set_columns(
            site=tables.mutations.site,
            node=tables.mutations.node,
            derived_state=tables.mutations.derived_state,
            derived_state_offset=tables.mutations.derived_state_offset,
        )
        ts = tables.tree_sequence()
        with pytest.raises(_msprime.LibraryError):
            msprime.sim_mutations(
                ts, rate=1, random_seed=1, keep=True, discrete_genome=False
            )

    def test_add_ancestral(self):
        ts = msprime.sim_ancestry(
            12, recombination_rate=3, random_seed=3, sequence_length=10
        )
        ts_mut = msprime.sim_mutations(ts, rate=0.1, random_seed=1)
        assert ts_mut.num_sites > 0
        with pytest.raises(_msprime.LibraryError):
            msprime.sim_mutations(
                ts_mut, rate=1, random_seed=1, keep=True, discrete_genome=True
            )
        ts_2mut = msprime.sim_mutations(
            ts_mut,
            rate=0.1,
            random_seed=3,
            discrete_genome=True,
            add_ancestral=True,
        )
        assert ts_2mut.num_mutations > ts_mut.num_mutations

    def test_keep_only_ancestral(self):
        # if timespan where mutations will be generated is younger than all
        # kept mutations, it shouldn't error out
        ts = msprime.sim_ancestry(
            12, recombination_rate=3, random_seed=3, sequence_length=10
        )
        ts_mut = msprime.sim_mutations(
            ts,
            rate=1,
            random_seed=1,
            discrete_genome=True,
            start_time=1.0,
            end_time=2.0,
        )
        ts_2mut = msprime.sim_mutations(
            ts_mut,
            rate=1,
            random_seed=3,
            discrete_genome=True,
            start_time=0.0,
            end_time=1.0,
        )
        assert ts_2mut.num_mutations > ts_mut.num_mutations

    def test_keep_ancestral_many_mutations(self):
        ts = msprime.sim_ancestry(5, random_seed=3)
        ts_mut1 = msprime.sim_mutations(ts, rate=1, random_seed=1)
        assert ts_mut1.num_sites == 1
        assert ts_mut1.num_mutations > 1
        with pytest.raises(_msprime.LibraryError, match="silent transition"):
            msprime.sim_mutations(ts_mut1, rate=1, random_seed=1, add_ancestral=True)

    def test_keep_ancestral_mutate(self):
        ts = msprime.sim_ancestry(5, random_seed=3)
        ts_mut1 = msprime.mutate(ts, rate=1, random_seed=1)
        ts_mut2 = msprime.mutate(ts_mut1, rate=1, random_seed=1, keep=True)
        assert set(ts_mut1.tables.sites.position) < set(ts_mut2.tables.sites.position)
        # Check that decoding variants succeeds
        assert len(list(ts_mut2.variants())) == ts_mut2.num_sites


class StatisticalTestMixin:

    p_threshold = 0.001

    def chisquare(self, observed, expected, p_th=p_threshold):
        obs = observed.reshape(np.product(observed.shape))
        exp = expected.reshape(np.product(observed.shape))
        assert np.all(obs[exp == 0] == 0)
        not_zeros = exp > 0
        if sum(not_zeros) > 1:
            chisq = stats.chisquare(obs[not_zeros], exp[not_zeros])
            assert chisq.pvalue > p_th

    def bonferroni(self, pvals):
        assert min(pvals) > self.p_threshold / len(pvals)

    def sign_tst(self, diffs):
        np = sum(diffs > 0)
        nm = sum(diffs < 0)
        pp = stats.binom.cdf(np, n=np + nm, p=0.5)
        pm = stats.binom.cdf(nm, n=np + nm, p=0.5)
        assert min(pp, pm) > self.p_threshold


@dataclasses.dataclass
class SlimMetadata:
    mutation_type_id: int
    selection_coeff: float
    subpop_index: int
    origin_generation: int
    nucleotide: int


class TestSLiMMutationModel:
    """
    Tests for the SLiM mutation generator.
    """

    def parse_slim_metadata(self, metadata):
        fmt = "<ifiib"
        md_size = 17
        assert struct.calcsize(fmt) == md_size
        assert len(metadata) % md_size == 0
        ret = []
        for j in range(len(metadata) // md_size):
            unpacked = struct.unpack(fmt, metadata[j * md_size : (j + 1) * md_size])
            ret.append(SlimMetadata(*unpacked))
        return ret

    def validate_slim_mutations(self, ts, mutation_type=0, slim_generation=1):
        # slim alleles should be lists of integers
        # and ancestral states the empty string
        for site in ts.sites():
            assert site.ancestral_state == ""
            alleles = {}
            for mutation in site.mutations:
                a = list(map(int, mutation.derived_state.split(",")))
                alleles[mutation.id] = a
                metadata = self.parse_slim_metadata(mutation.metadata)
                assert len(metadata) == len(a)
                if mutation.parent == tskit.NULL:
                    assert len(a) == 1
                else:
                    parent_allele = alleles[mutation.parent]
                    assert a[:-1] == parent_allele
                    parent_metadata = self.parse_slim_metadata(
                        ts.mutation(mutation.parent).metadata
                    )
                    assert metadata[:-1] == parent_metadata
                for md in metadata:
                    assert md.mutation_type_id == mutation_type
                    assert md.selection_coeff == 0
                    assert md.subpop_index == tskit.NULL
                    assert md.nucleotide == -1
                # time should only match on the last one
                assert md.origin_generation == slim_generation - int(mutation.time)

    def run_mutate(
        self,
        ts,
        rate=1,
        random_seed=42,
        mutation_type=0,
        mutation_id=0,
        slim_generation=1,
    ):

        model = msprime.SLiMMutationModel(
            type=mutation_type, next_id=mutation_id, slim_generation=slim_generation
        )
        mts1 = msprime.sim_mutations(
            ts, rate=rate, random_seed=random_seed, model=model, discrete_genome=True
        )
        assert mts1.num_mutations == model.next_id

        model = PythonSLiMMutationModel(
            mutation_type=mutation_type, next_id=mutation_id
        )
        mts2 = py_sim_mutations(
            ts, rate=rate, random_seed=random_seed, model=model, discrete_genome=True
        )

        t1 = mts1.dump_tables()
        t2 = mts2.dump_tables()
        assert t1.sites == t2.sites
        # Drop the mutation metadata - we're validating that elsewhere and
        # it's not worth complicating the Python generator with it.
        t1.mutations.set_columns(
            site=t1.mutations.site,
            node=t1.mutations.node,
            parent=t1.mutations.parent,
            time=t1.mutations.time,
            derived_state=t1.mutations.derived_state,
            derived_state_offset=t1.mutations.derived_state_offset,
        )
        assert t1.mutations == t2.mutations
        return mts1

    def test_slim_mutation_type(self):
        ts = msprime.sim_ancestry(4, sequence_length=2, random_seed=5)
        for mutation_type in range(1, 10):
            mts = self.run_mutate(
                ts, rate=5.0, random_seed=23, mutation_type=mutation_type
            )
            assert mts.num_mutations > 10
            self.validate_slim_mutations(mts, mutation_type=mutation_type)

    def test_slim_generation(self):
        ts = msprime.simulate(4, length=2, random_seed=5)
        for slim_generation in [-100, 0, 1, 256]:
            mts = self.run_mutate(
                ts, rate=5.0, random_seed=23, slim_generation=slim_generation
            )
            assert mts.num_mutations > 10
            self.validate_slim_mutations(mts, slim_generation=slim_generation)

    def test_binary_n_4_low_rate(self):
        ts = msprime.sim_ancestry(4, sequence_length=10, random_seed=5)
        mts = self.run_mutate(ts, rate=0.1, random_seed=23)
        assert mts.num_mutations > 1
        self.validate_slim_mutations(mts)

    def test_binary_n_4_high_rate(self):
        ts = msprime.sim_ancestry(4, sequence_length=2, random_seed=5)
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate_slim_mutations(mts)

    def test_binary_n_8_low_rate(self):
        ts = msprime.sim_ancestry(8, sequence_length=10, random_seed=50)
        mts = self.run_mutate(ts, rate=0.1, random_seed=342)
        assert mts.num_mutations > 1
        self.validate_slim_mutations(mts)

    def test_binary_n_8_high_rate(self):
        ts = msprime.sim_ancestry(8, sequence_length=10, random_seed=5)
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate_slim_mutations(mts)

    def test_binary_incomplete_trees(self):
        ts = msprime.sim_ancestry(8, sequence_length=5, random_seed=50, end_time=0.1)
        assert ts.first().num_roots > 1
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate_slim_mutations(mts)

    def test_binary_many_trees(self):
        ts = msprime.sim_ancestry(
            8,
            sequence_length=5,
            recombination_rate=5,
            random_seed=50,
            discrete_genome=False,
        )
        assert ts.num_trees > 20
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate_slim_mutations(mts)


class TestInfiniteAllelesMutationModel:
    """
    Tests for the Infinite alleles mutation model
    """

    def validate(self, ts, start_allele=0):
        allele = start_allele
        for site in ts.sites():
            assert int(site.ancestral_state) == allele
            allele += 1
            for mutation in site.mutations:
                assert int(mutation.derived_state) == allele
                allele += 1

    def validate_unique_alleles(self, ts):
        for site in ts.sites():
            alleles = [m.derived_state for m in site.mutations] + [site.ancestral_state]
            assert len(alleles) == len(set(alleles))

    def run_mutate(self, ts, rate=1, random_seed=42, start_allele=0):

        model = msprime.InfiniteAllelesMutationModel(start_allele=start_allele)
        mts1 = msprime.sim_mutations(
            ts, rate=rate, random_seed=random_seed, model=model, discrete_genome=True
        )
        model = PythonInfiniteAllelesMutationModel(start_allele=start_allele)
        mts2 = py_sim_mutations(
            ts, rate=rate, random_seed=random_seed, model=model, discrete_genome=True
        )

        t1 = mts1.dump_tables()
        t2 = mts2.dump_tables()
        assert t1.sites == t2.sites
        assert t1.mutations == t2.mutations
        return mts1

    def test_binary_n_4_low_rate(self):
        ts = msprime.sim_ancestry(4, sequence_length=10, random_seed=5)
        mts = self.run_mutate(ts, rate=0.1, random_seed=23)
        assert mts.num_mutations > 1
        self.validate(mts)

    def test_binary_n_4_high_rate(self):
        ts = msprime.sim_ancestry(4, sequence_length=2, random_seed=5)
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate(mts)

    def test_binary_n_8_low_rate(self):
        ts = msprime.sim_ancestry(8, sequence_length=10, random_seed=50)
        mts = self.run_mutate(ts, rate=0.1, random_seed=342)
        assert mts.num_mutations > 1
        self.validate(mts)

    def test_binary_n_8_high_rate(self):
        ts = msprime.sim_ancestry(8, sequence_length=10, random_seed=5)
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate(mts)

    def test_binary_incomplete_trees(self):
        ts = msprime.sim_ancestry(8, sequence_length=5, random_seed=50, end_time=0.1)
        assert ts.first().num_roots > 1
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate(mts)

    def test_binary_many_trees(self):
        ts = msprime.sim_ancestry(
            8, sequence_length=50, recombination_rate=5, random_seed=50
        )
        assert ts.num_trees > 20
        mts = self.run_mutate(ts, rate=2.0, random_seed=23)
        assert mts.num_mutations > 10
        self.validate(mts)

    def test_allele_overflow(self):
        ts = msprime.sim_ancestry(4, sequence_length=2, random_seed=5)
        start_allele = 2 ** 64 - 1
        model = msprime.InfiniteAllelesMutationModel(start_allele=start_allele)
        mts = msprime.sim_mutations(
            ts, rate=1, random_seed=32, model=model, discrete_genome=True
        )
        assert mts.num_sites > 0
        assert mts.num_mutations > 0
        allele = start_allele
        first_allele = True
        for site in mts.sites():
            assert int(site.ancestral_state) == allele
            if first_allele:
                allele = 0
                first_allele = False
            else:
                allele += 1
            for mutation in site.mutations:
                assert int(mutation.derived_state) == allele
                allele += 1

    def test_non_discrete_sites(self):
        ts = msprime.sim_ancestry(4, sequence_length=2, random_seed=5)
        model = msprime.InfiniteAllelesMutationModel()
        mts = msprime.sim_mutations(
            ts, rate=1, random_seed=32, model=model, discrete_genome=False
        )
        assert mts.num_sites > 0
        assert mts.num_mutations == mts.num_sites
        self.validate(mts)
        for site in mts.sites():
            assert len(site.mutations) == 1

    def test_keep_mutations(self):
        t = tskit.TableCollection(sequence_length=1)
        t.nodes.add_row(time=0)
        t.nodes.add_row(time=10)
        t.sites.add_row(ancestral_state="0", position=0)
        t.mutations.add_row(derived_state="1", node=1, site=0, time=10)
        t.edges.add_row(parent=1, child=0, left=0, right=1)
        model = msprime.InfiniteAllelesMutationModel(start_allele=2)
        ts = msprime.sim_mutations(
            t.tree_sequence(),
            rate=1,
            model=model,
            start_time=2,
            random_seed=1,
            keep=True,
            discrete_genome=True,
            add_ancestral=True,
        )
        self.validate_unique_alleles(ts)


class TestPythonMutationGenerator:
    """
    Tests for the python implementation of the mutation generator
    (but, note it is compared to the C version extensively above).
    """

    def verify(self, ts, **kwargs):
        rates = [0, 0.01, 0.2]
        discretes = [True, False]
        keeps = [True, False]
        for rate, keep, discrete_genome in itertools.product(rates, keeps, discretes):
            ts1 = msprime.sim_mutations(
                ts,
                rate=rate,
                keep=keep,
                discrete_genome=discrete_genome,
                add_ancestral=True,
                **kwargs,
            )
            ts2 = py_sim_mutations(
                ts, rate=rate, keep=keep, discrete_genome=discrete_genome, **kwargs
            )
            tables1 = ts1.dump_tables()
            tables2 = ts2.dump_tables()
            tables1.provenances.clear()
            tables2.provenances.clear()
            assert tables1 == tables2

    def test_single_tree_no_mutations(self):
        ts = msprime.sim_ancestry(10, sequence_length=100, random_seed=1234)
        self.verify(ts, random_seed=234)

    def test_single_tree_mutations(self):
        ts = msprime.simulate(10, length=100, mutation_rate=0.1, random_seed=1234)
        assert ts.num_sites > 0
        self.verify(ts, random_seed=34)

    def test_many_trees_no_mutations(self):
        ts = msprime.sim_ancestry(
            10, sequence_length=100, recombination_rate=0.1, random_seed=123
        )
        assert ts.num_trees > 1
        self.verify(ts, random_seed=789)

    def test_many_trees_mutations(self):
        ts = msprime.simulate(
            10, length=100, mutation_rate=0.1, recombination_rate=2, random_seed=123
        )
        assert ts.num_trees > 1
        assert ts.num_sites > 1
        self.verify(ts, random_seed=789)


####################################################
# Python implementation of lib/mutgen.c algorithms #
####################################################


def py_sim_mutations(
    ts,
    rate=None,
    random_seed=None,
    model=None,
    keep=False,
    discrete_genome=False,
    sequential_only=True,
):
    """
    Same interface as mutations.sim_mutations() and should provide identical results.
    """
    if rate is None:
        rate = 0
    if model is None:
        model = msprime.JC69MutationModel()
    if isinstance(model, PythonMutationModel):
        py_model = model
    else:
        py_model = PythonMutationMatrixModel(
            alleles=model.alleles,
            root_distribution=model.root_distribution,
            transition_matrix=model.transition_matrix,
        )
    tables = ts.dump_tables()
    mutmap = msprime.RateMap([0, ts.sequence_length], [rate])
    mutgen = PythonMutationGenerator(mutmap, py_model)
    return mutgen.generate(
        tables,
        random_seed,
        keep=keep,
        discrete_genome=discrete_genome,
        sequential_only=sequential_only,
    )


@dataclasses.dataclass
class Site:
    position: float
    ancestral_state: str
    metadata: bytes
    mutations: List[Mutation]
    new: bool

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


@dataclasses.dataclass
class Mutation:
    node: int
    derived_state: str
    parent: int
    metadata: bytes
    time: float
    new: bool
    keep: bool
    id: int  # noqa: A003

    def __str__(self):
        if self.parent is None:
            parent_id = None
        else:
            parent_id = self.parent.id
        s = f"\t{self.id}\t\tnode: {self.node}\tparent: {parent_id}"
        s += f"\ttime: {self.time}\t{self.derived_state}\t{self.metadata}"
        s += f"\t(new: {self.new})\tkeep: {self.keep}]"
        return s


class PythonMutationModel:
    # Base class for mutation models, which must define these methods:

    def root_allele(self, rng):
        pass

    def transition_allele(self, rng, current_allele):
        pass


@dataclasses.dataclass
class PythonSLiMMutationModel(PythonMutationModel):
    mutation_type: int = 0
    next_id: int = 0

    def root_allele(self, rng):
        return ""

    def transition_allele(self, rng, current_allele):
        out = current_allele
        if len(current_allele) > 0:
            out += ","
        out += str(self.next_id)
        self.next_id += 1
        return out


@dataclasses.dataclass
class PythonInfiniteAllelesMutationModel(PythonMutationModel):
    start_allele: int = 0
    next_allele: int = 0

    def __init__(self, start_allele):
        self.next_allele = start_allele

    def make_allele(self):
        ret = str(self.next_allele)
        self.next_allele += 1
        return ret

    def root_allele(self, rng):
        return self.make_allele()

    def transition_allele(self, rng, current_allele):
        return self.make_allele()


@dataclasses.dataclass
class PythonMutationMatrixModel(PythonMutationModel):
    # for compatability with the C code we provide alleles as bytes,
    # but we want them as strings here for simplicity.
    alleles: List[bytes]
    # Taking a short-cut here with the annotations
    root_distribution: Any
    transition_matrix: Any

    def choose_allele(self, rng, distribution):
        u = rng.flat(0, 1)
        j = 0
        while u > distribution[j]:
            u -= distribution[j]
            j += 1
        return self.alleles[j]

    def root_allele(self, rng):
        return self.choose_allele(rng, self.root_distribution)

    def transition_allele(self, rng, current_allele):
        j = self.alleles.index(current_allele)
        return self.choose_allele(rng, self.transition_matrix[j])


def cmp_mutation(a, b):
    # Sort mutations by decreasing time and increasing parent,
    # but preserving order of any kept mutations (assumed to be
    # in order already). Kept mutations are given an id that is
    # their order in the initial tables, and new mutations have id -1.
    out = a.id * (not a.new) - b.id * (not b.new)
    if out == 0:
        out = b.time - a.time
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
        self.model = model
        self.sites = {}

    def print_state(self):
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
                    time=mutation_row.time,
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
                        time=mutation.time,
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

    def place_mutations(self, tables, discrete_genome=False):
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
                site_left = np.ceil(left) if discrete_genome else left
                site_right = np.ceil(right) if discrete_genome else right
                assert site_left <= site_right
                assert map_position[index] <= left
                assert right <= map_position[index + 1]
                assert right <= edge.right
                # Generate the mutations.
                rate = self.rate_map.rate[index]
                mu = rate * (site_right - site_left) * branch_length
                for _ in range(self.rng.poisson(mu)[0]):
                    position = self.rng.flat(site_left, site_right)[0]
                    if discrete_genome:
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

    def choose_alleles(self, tree_parent, site, mutation_id_offset, sequential_only):
        if site.new:
            site.ancestral_state = self.model.root_allele(self.rng)
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
                pa = site.ancestral_state
                assert mut.parent is None
            else:
                assert u in bottom_mutation
                parent_mut = bottom_mutation[u]
                mut.parent = parent_mut
                assert mut.time <= parent_mut.time, "Parent after child mutation."
                if sequential_only and (
                    mut.time > parent_mut.time or (parent_mut.new and not mut.new)
                ):
                    raise ValueError(
                        "Generated mutation appears above "
                        "an existing mutation: cannot apply "
                        "finite sites mutations to a earlier "
                        "time period than where they already exist."
                    )
                if mut.new:
                    pa = parent_mut.derived_state

            if mut.new:
                da = self.model.transition_allele(self.rng, pa)
                if da == pa:
                    mut.keep = False
                else:
                    mut.derived_state = da
            if mut.keep:
                bottom_mutation[mut.node] = mut

    def apply_mutations(self, tables, sequential_only):
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
                self.choose_alleles(
                    tree_parent, site, mutation_id_offset, sequential_only
                )
                num_mutations = len(site.mutations)
                mutation_id_offset += num_mutations
                j += 1

    def generate(
        self, tables, seed, keep=False, discrete_genome=False, sequential_only=True
    ):
        self.rng = _msprime.RandomGenerator(seed)
        if keep:
            self.initialise_sites(tables)
        tables.sites.clear()
        tables.mutations.clear()
        self.place_mutations(tables, discrete_genome=discrete_genome)
        self.apply_mutations(tables, sequential_only=sequential_only)
        self.populate_tables(tables)
        self.record_provenance(tables, seed, keep, discrete_genome)
        return tables.tree_sequence()

    def record_provenance(self, tables, seed, keep, discrete_genome):
        parameters = {
            "command": "mutate",
            "random_seed": seed,
            "keep": keep,
            "discrete_genome": discrete_genome,
            "rate_map": None,  # not working??
            "model": None,  # TODO
        }
        provenance_dict = msprime.provenance.get_provenance_dict(parameters)
        encoded_provenance = msprime.provenance.json_encode_provenance(provenance_dict)
        tables.provenances.add_row(encoded_provenance)


class TestMutationModelFactory:
    """
    Tests that the mutation_model_factory function.
    """

    def test_bad_model_names(self):
        for bad_model in ["NOT", "", "MODEL", "gtr", "slim"]:
            with pytest.raises(ValueError):
                msprime.mutation_model_factory(bad_model)

    def test_named_model_variants(self):
        mutation_models = {
            "infinite_alleles": msprime.InfiniteAllelesMutationModel,
            "binary": msprime.BinaryMutationModel,
            "jc69": msprime.JC69MutationModel,
            "blosum62": msprime.BLOSUM62MutationModel,
            "pam": msprime.PAMMutationModel,
        }
        for name, model_class in mutation_models.items():

            model = msprime.mutation_model_factory(model=name.upper())
            assert isinstance(model, model_class)
            model = msprime.mutation_model_factory(model=name.title())
            assert isinstance(model, model_class)
            model = msprime.mutation_model_factory(model=name)
            assert isinstance(model, model_class)

    def test_bad_models(self):
        for bad_type in [1234, {}]:
            with pytest.raises(TypeError):
                msprime.mutation_model_factory(model=bad_type)

    def test_model_instances(self):
        models = [
            msprime.SLiMMutationModel(0, 0),
            msprime.InfiniteAllelesMutationModel(),
            msprime.BinaryMutationModel(),
            msprime.JC69MutationModel(),
            msprime.HKYMutationModel(0.75),
            msprime.F84MutationModel(0.75),
            msprime.GTRMutationModel([1 / 6] * 6),
            msprime.BLOSUM62MutationModel(),
            msprime.PAMMutationModel(),
        ]
        for model in models:
            new_model = msprime.mutation_model_factory(model=model)
            assert new_model is model
            assert new_model.__dict__ == model.__dict__


class TestModelClasses:
    def test_slim(self):
        m = msprime.SLiMMutationModel(type=1, next_id=2, slim_generation=9)
        assert m.next_id == 2
        assert m.slim_generation == 9
        assert (
            str(m) == "Mutation model for SLiM mutations of type m1\n" "  next ID: 2\n"
        )

    def test_infinite_alleles(self):
        m = msprime.InfiniteAllelesMutationModel(start_allele=1)
        assert (
            str(m) == "Infinite alleles mutation model, beginning with"
            " allele 1\n    next allele: 1\n"
        )


class TestDeprecatedApis:
    def test_simulate_mutate_keep(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=2)
        assert ts.num_sites > 0
        mts = msprime.mutate(ts, rate=1, random_seed=3, keep=True)
        assert set(mts.tables.sites.position) > set(ts.tables.sites.position)

    def test_simulate_sim_mutations(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=2)
        assert ts.num_sites > 0
        mts = msprime.sim_mutations(
            ts, rate=1, random_seed=3, keep=True, discrete_genome=False
        )
        assert set(mts.tables.sites.position) > set(ts.tables.sites.position)

    def test_sim_ancestry_mutate(self):
        ts = msprime.sim_ancestry(10, random_seed=2)
        assert ts.num_sites == 0
        mts = msprime.mutate(ts, rate=1, random_seed=3)
        assert mts.num_sites > 0
