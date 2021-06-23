import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
import pytest
import tskit

import msprime
from msprime import pedigrees


class TestPedigree(unittest.TestCase):
    """
    Tests for the wf_ped model.
    """

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="msp_ped_testcase_")
        self.temp_pedigree_text_file = os.path.join(self.temp_dir, "pedigree.txt")
        self.temp_pedigree_array_file = os.path.join(self.temp_dir, "pedigree.npy")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_simple_case(self):
        # Simple test to check that the current code is still roughly working.
        # Should be replace by better tests.
        individual = np.array([1, 2, 3, 4])
        parents = np.array([-1, -1, -1, -1, 0, 1, 1, 0]).reshape(-1, 2)
        time = np.array([1, 1, 0, 0])
        is_sample = np.array([0, 0, 1, 1])

        model = msprime.WrightFisherPedigree()
        ped = pedigrees.Pedigree(
            individual, parents, time, is_sample, sex=None, ploidy=2
        )
        ts = msprime.simulate(2, pedigree=ped, model=model, recombination_rate=1)
        table = ts.tables.individuals
        assert np.all(table[0].parents == [-1, -1])
        assert np.all(table[1].parents == [-1, -1])
        # FIXME these are currently being reversed.
        # assert np.all(table[2].parents == [0, 1])
        # assert np.all(table[3].parents == [1, 0])

    @unittest.skip("Currently broken")
    def test_pedigree_replicates(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        time = np.array([0, 0, 1, 1])
        is_sample = np.array([1, 1, 0, 0])

        model = msprime.WrightFisherPedigree()
        ped = pedigrees.Pedigree(
            individual, parents, time, is_sample, sex=None, ploidy=2
        )
        replicates = msprime.simulate(
            2, pedigree=ped, model=model, recombination_rate=1, num_replicates=100
        )
        for ts in replicates:
            assert ts is not None

    @unittest.skip("Currently broken")
    def test_pedigree_samples(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        time = np.array([0, 0, 1, 1])

        ped = pedigrees.Pedigree(individual, parents, time)
        # This interface is pretty poor - what we want is for the
        # samples argument to simulate to be interpreted as individual
        # IDs. For now, let's just leave it as it is, though.
        ped.set_samples(2)
        ts = msprime.simulate(2, pedigree=ped, model="wf_ped")
        assert ts.num_edges > 0

        ped.set_samples(sample_IDs=[1, 2])
        ts = msprime.simulate(2, pedigree=ped, model="wf_ped")
        assert ts is not None
        with pytest.raises(ValueError):
            ped.set_samples(sample_IDs=[1, 1])
        with pytest.raises(ValueError):
            ped.set_samples(sample_IDs=[1, 3])
        with pytest.raises(NotImplementedError):
            ped.set_samples(sample_IDs=[1, 3], probands_only=False)

        ped.set_samples(sample_IDs=[1, 2])
        assert ped.get_proband_indices() == [0, 1]

    def test_pedigree_read_write(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        time = np.array([0, 0, 1, 1])
        is_sample = np.array([1, 1, 0, 0])

        ped = pedigrees.Pedigree(
            individual, parents, time, is_sample, sex=None, ploidy=2
        )

        ped.save_txt(self.temp_pedigree_text_file)
        with pytest.raises(NotImplementedError):
            pedigrees.Pedigree.read_txt(
                self.temp_pedigree_text_file,
                sex_col=4,
            )
        # FIXME
        # The compute_times should be done automatically in this case .
        # ped_from_txt = pedigrees.Pedigree.read_txt(
        #     self.temp_pedigree_text_file, time_col=None
        # )
        # ts = msprime.simulate(2, pedigree=ped_from_txt, model="wf_ped")
        # self.assertTrue(ts is not None)
        # ped_from_txt = pedigrees.Pedigree.read_txt(
        #     self.temp_pedigree_text_file, time_col=3
        # )
        # ts = msprime.simulate(2, pedigree=ped_from_txt, model="wf_ped")
        # self.assertTrue(ts is not None)

        ped.save_npy(self.temp_pedigree_array_file)
        ped_from_npy = pedigrees.Pedigree.read_npy(self.temp_pedigree_array_file)
        # TODO compre this to the file above.
        assert isinstance(ped_from_npy, pedigrees.Pedigree)

    def test_pedigree_times(self):
        individual = np.array([1, 2, 3, 4])
        time = np.array([0, 0, 1, 1])

        parent_IDs = np.array([3, 4, 3, 4, 0, 0, 0, 0]).reshape(-1, 2)
        estimated_times = pedigrees.Pedigree.get_times(
            individual, parent_IDs=parent_IDs, check=True
        )
        assert (time == estimated_times).all()

    def test_pedigree_utils(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        parent_IDs = np.array([3, 4, 3, 4, 0, 0, 0, 0]).reshape(-1, 2)

        bad_individual = np.array([0, 1, 2, 3])
        with pytest.raises(ValueError):
            pedigrees.Pedigree.parent_ID_to_index(bad_individual, parent_IDs)
        assert (
            pedigrees.Pedigree.parent_ID_to_index(individual, parent_IDs) == parents
        ).all()
        assert (
            pedigrees.Pedigree.parent_index_to_ID(individual, parents) == parent_IDs
        ).all()

    def test_pedigree_sanity_checks(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        time = np.array([0, 0, 1, 1])
        is_sample = np.array([1, 1, 0, 0])

        with pytest.raises(NotImplementedError):
            pedigrees.Pedigree(
                individual=individual,
                parents=parents,
                time=time,
                ploidy=1,
            )
        bad_parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 4)
        with pytest.raises(ValueError):
            pedigrees.Pedigree(
                individual=individual,
                parents=bad_parents,
                time=time,
            )
        bad_individual = np.array([-1, 2, 3, 4])
        with pytest.raises(ValueError):
            pedigrees.Pedigree(
                individual=bad_individual,
                parents=parents,
                time=time,
            )

        bad_times = np.array([1, 1, 1, 1])
        ped = pedigrees.Pedigree(
            individual, parents, bad_times, is_sample, sex=None, ploidy=2
        )
        with pytest.raises(ValueError):
            ped.check_times(
                individual=ped.individual,
                parents=ped.parents,
                time=ped.time,
            )

        ped = pedigrees.Pedigree(individual, parents, time, sex=None, ploidy=2)
        with pytest.raises(ValueError):
            ped.set_samples()
        with pytest.raises(ValueError):
            ped.set_samples(num_samples=10)
        with pytest.raises(ValueError):
            ped.set_samples(num_samples=2, sample_IDs=[1, 2])

        with pytest.raises(ValueError):
            ped.get_times(individual=ped.individual)


def simulate_pedigree(
    num_founders,
    num_children_prob=(0, 0, 1),
    num_generations=3,
    sequence_length=1,
    random_seed=42,
) -> tskit.TableCollection:
    """
    Simulates pedigree.

    num_founders: Number of founders for pedigree
    num_children_prob: Array-like object of probabilities for number of children
        produced by each pair of mates. Index 0 corresponds to probability of
        0 children, index 1 corresponds to probability of 1 child, etc. Sum of
        entries must be 1. Defaults to (0, 0, 1), i.e., there are always exactly
        two children.
    num_generations: Number of generations to attempt to simulate
    sequence_length: The sequence_length of the output tables.
    random_seed: Random seed.
    """
    rng = np.random.RandomState(random_seed)
    tables = tskit.TableCollection(sequence_length=sequence_length)
    # Add single population for simplicity.
    tables.populations.add_row()

    def add_individual(generation, parents=(-1, -1)):
        ind_id = tables.individuals.add_row(parents=parents)
        time = num_generations - generation - 1
        for _ in parents:
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE if time == 0 else 0,
                time=time,
                population=0,
                individual=ind_id,
            )
        return ind_id

    curr_gen = [add_individual(0) for _ in range(num_founders)]
    for generation in range(1, num_generations):
        num_pairs = len(curr_gen) // 2
        if num_pairs == 0 and num_children_prob[0] != 1:
            raise Exception(
                f"Not enough people to make children in generation {generation}"
            )
        all_parents = rng.choice(curr_gen, size=(num_pairs, 2), replace=False)
        curr_gen = []
        for parents in all_parents:
            num_children = rng.choice(len(num_children_prob), p=num_children_prob)
            for _ in range(num_children):
                ind_id = add_individual(generation, parents=parents.astype(np.int32))
                curr_gen.append(ind_id)
    tables.build_index()
    return tables


class TestPedigreeSimulation:
    """
    Tests for simple pedigree simulator
    """

    def simple_sim(
        self, num_founders=2, num_children_prob=None, num_generations=2, random_seed=42
    ):
        if num_children_prob is None:
            num_children_prob = [0, 0, 1]
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=num_children_prob,
            num_generations=num_generations,
            random_seed=random_seed,
        )
        self.verify_nodes(tables)
        return tables.individuals

    def verify_nodes(self, tables):
        """
        Check that the nodes and individuals in the specified tables have the
        required properties.
        """
        ts = tables.tree_sequence()
        for ind in ts.individuals():
            assert len(ind.nodes) == 2
            assert len({ts.node(node).flags for node in ind.nodes}) == 1
            assert len({ts.node(node).time for node in ind.nodes}) == 1
            assert len({ts.node(node).population for node in ind.nodes}) == 1

    def test_one_generation_no_children(self):
        num_founders = 8
        tb = self.simple_sim(num_founders=num_founders, num_generations=1)
        assert len(tb) == num_founders

        # check that all parents of founders are missing
        assert all([np.array_equal(row.parents, [-1, -1]) for row in tb])

    def test_one_trio(self):
        tb = self.simple_sim(
            num_founders=2, num_children_prob=[0, 1], num_generations=2
        )
        assert len(tb) == 3
        assert np.array_equal(tb[2].parents, [1, 0])
        assert all([np.array_equal(tb[idx].parents, [-1, -1]) for idx in range(2)])

    def test_grandparents(self):
        tb = self.simple_sim(
            num_founders=4, num_children_prob=[0, 1], num_generations=3
        )
        assert len(tb) == 7
        assert all([np.array_equal(tb[idx].parents, [-1, -1]) for idx in range(4)])
        assert set(tb[4].parents).union(tb[5].parents) == set(range(4))
        assert set(tb[6].parents) == {4, 5}

    def test_insufficient_founders(self):
        with pytest.raises(Exception):
            self.simple_sim(num_founders=1, num_children_prob=[0, 1])
        with pytest.raises(Exception):
            self.simple_sim(num_founders=3, num_children_prob=[0, 1], num_generations=3)

    @pytest.mark.parametrize("num_children", range(5))
    def test_nonrandom_child_prob(self, num_children):
        tb = self.simple_sim(
            num_founders=2,
            num_children_prob=[0] * num_children + [1],
            num_generations=2,
        )
        assert len(tb) == 2 + num_children

    @pytest.mark.parametrize("num_children_prob", [[0.5, 0.5], [0, 0.5, 0.5]])
    def test_expected_num_children(self, num_children_prob):
        n_reps = 1000
        num_founders = 2
        num_children = []
        for rep in range(n_reps):
            tb = self.simple_sim(
                num_founders=num_founders,
                num_children_prob=num_children_prob,
                num_generations=2,
                random_seed=rep,
            )
            num_children.append(len(tb) - num_founders)
        expected_num_children = sum(
            prob * n
            for prob, n in zip(num_children_prob, range(len(num_children_prob)))
        )
        assert np.mean(num_children) == pytest.approx(expected_num_children, rel=0.1)

    def test_bad_num_children_prob(self):
        with pytest.raises(ValueError):
            self.simple_sim(num_children_prob=[2])
        with pytest.raises(ValueError):
            self.simple_sim(num_children_prob=[1, 1])

    def test_valid_pedigree(self):
        tb = self.simple_sim(
            num_founders=128,
            num_generations=10,
        )
        tc = tskit.TableCollection(1)
        tc.individuals.metadata_schema = tb.metadata_schema
        for row in tb:
            tc.individuals.append(row)
        tc.tree_sequence()  # creating tree sequence should succeed


if __name__ == "__main__":
    # Easy way to output a simulated pedigree for command line testing.
    # Run with python3 tests/test_pedigree.py NUM_FOUNDERS NUM_GENERATIONS > out.trees
    num_founders = int(sys.argv[1])
    num_generations = int(sys.argv[2])
    tables = simulate_pedigree(
        num_founders=num_founders,
        num_generations=num_generations,
    )
    ts = tables.tree_sequence()
    ts.dump(sys.stdout)
