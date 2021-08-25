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
                parents = np.sort(parents).astype(np.int32)
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
        print(tables.individuals)
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
        assert np.array_equal(tb[2].parents, [0, 1])
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


class TestSimulateThroughPedigree:
    def verify(self, input_tables):
        initial_state = input_tables.tree_sequence()
        ts = msprime.sim_ancestry(
            model="wf_ped",
            initial_state=initial_state,
            recombination_rate=0,
            population_size=1,
            random_seed=1,
        )
        output_tables = ts.tables

        input_tables.individuals.assert_equals(output_tables.individuals)
        input_tables.nodes.assert_equals(output_tables.nodes)

        for node in ts.nodes():
            assert node.individual != tskit.NULL

        # Make the set of ancestors for each individual
        ancestors = [set() for _ in ts.individuals()]
        for individual in ts.individuals():
            for parent_id in individual.parents:
                stack = [parent_id]
                while len(stack) > 0:
                    parent_id = stack.pop()
                    if parent_id != tskit.NULL:
                        ancestors[individual.id].add(parent_id)
                        for u in ts.individual(parent_id).parents:
                            stack.append(u)
        # Not used because certain tests fail asserts
        # tests failed: test_shallow[100,1000], test_deep, test_many_children,
        # test_half_sibs, test_double_inbreeding_sibs
        #
        # # founder nodes have correct time
        # founder_ids = {
        #     ind.id for ind in ts.individuals() if len(ancestors[ind.id]) == 0
        # }
        # for founder in founder_ids:
        #    for node_id in ts.individual(founder).nodes:
        #        node = ts.node(node_id)
        # for tree in ts.trees():
        #     for node_id in tree.nodes():
        #         node = ts.node(node_id)
        #         individual = ts.individual(node.individual)
        #         if tree.parent(node_id) != tskit.NULL:
        #             parent_node = ts.node(tree.parent(node_id))
        #             # tests failed: test_shallow[100,1000],
        #             # test_deep, test_many_children, test_half_sibs
        #             assert parent_node.individual in ancestors[individual.id]
        #         else:
        #             # tests failed: test_double_inbreeding_sibs
        #             assert node.individual in founder_ids

    def add_individual(self, tables, time, parents=(-1, -1)):
        parents = np.sort(parents).astype("int32")
        ind_id = tables.individuals.add_row(
            parents=parents, flags=len(tables.individuals)
        )
        for _ in parents:
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE if time == 0 else 0,
                time=time,
                population=0,
                individual=ind_id,
            )
        return ind_id

    @pytest.mark.parametrize("num_founders", [2, 3, 100, 1000])
    def test_shallow(self, num_founders):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0, 0, 1],
            num_generations=2,
            sequence_length=100,
        )
        self.verify(tables)

    @pytest.mark.parametrize("num_founders", [2, 3, 10, 20])
    def test_deep(self, num_founders):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0, 0, 1],
            num_generations=16,  # Takes a long time if this is increased
            sequence_length=100,
        )
        self.verify(tables)

    @pytest.mark.parametrize("num_founders", [2, 3, 10, 20])
    def test_many_children(self, num_founders):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0] * 99
            + [1],  # Each pair of parents will always have 100 children
            num_generations=2,
        )
        self.verify(tables)

    @unittest.skip("Currently broken")
    @pytest.mark.parametrize("num_founders", [2, 3, 10, 20])
    def test_unrelated(self, num_founders):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0],
            num_generations=1,
        )
        self.verify(tables)

    @unittest.skip("Currently broken")
    def test_duo(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        parent = self.add_individual(tables, time=1)
        self.add_individual(tables, time=0, parents=[parent, -1])

        tables.build_index()
        self.verify(tables)

    def test_trio(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        gen0 = [self.add_individual(tables, time=1) for _ in range(2)]
        self.add_individual(tables, time=0, parents=gen0)

        tables.build_index()
        self.verify(tables)

    def test_mixed_generation_parents(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        gen0 = [self.add_individual(tables, time=2) for _ in range(2)]
        gen1 = [self.add_individual(tables, time=1, parents=gen0)]
        self.add_individual(tables, time=0, parents=(gen0[0], gen1[0]))

        tables.build_index()
        self.verify(tables)

    def test_inbreeding_sibs(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        parents = [self.add_individual(tables, time=2) for _ in range(2)]
        sibs = [self.add_individual(tables, time=1, parents=parents) for _ in range(2)]
        self.add_individual(tables, time=0, parents=sibs)

        tables.build_index()
        self.verify(tables)

    def test_double_inbreeding_sibs(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        parents = [self.add_individual(tables, time=3) for _ in range(2)]
        sibs1 = [self.add_individual(tables, time=2, parents=parents) for _ in range(2)]
        sibs2 = [self.add_individual(tables, time=1, parents=sibs1) for _ in range(2)]
        self.add_individual(tables, time=0, parents=sibs2)

        tables.build_index()
        self.verify(tables)

    def test_half_sibs(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        parents = [self.add_individual(tables, time=1) for _ in range(3)]
        self.add_individual(tables, time=0, parents=parents[:2])
        self.add_individual(tables, time=0, parents=parents[1:])

        tables.build_index()
        print(tables.individuals)
        self.verify(tables)

    def test_inbreeding_half_sibs(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        parents = [self.add_individual(tables, time=2) for _ in range(3)]
        half_sib1 = self.add_individual(tables, time=1, parents=parents[:2])
        half_sib2 = self.add_individual(tables, time=1, parents=parents[1:])
        self.add_individual(tables, time=0, parents=(half_sib1, half_sib2))

        tables.build_index()
        self.verify(tables)

    def test_oedipus(self):
        tables = tskit.TableCollection(sequence_length=100)
        tables.populations.add_row()

        parents = [self.add_individual(tables, time=2) for _ in range(2)]
        offspring = self.add_individual(tables, time=1, parents=parents)
        self.add_individual(tables, time=0, parents=(parents[0], offspring))

        tables.build_index()
        self.verify(tables)


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
