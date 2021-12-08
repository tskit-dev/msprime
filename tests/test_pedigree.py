import dataclasses
import sys
import tempfile
from dataclasses import asdict

import numpy as np
import pytest
import tskit

import msprime
from msprime import _msprime
from msprime import pedigrees


# Note: these functions are left over from older code when we didn't
# have the PedigreeBuilder. New code should use this interface for
# building test pedigrees.
def get_base_tables(sequence_length, num_populations=1):
    demography = msprime.Demography.isolated_model([1] * num_populations)
    tables = pedigrees.PedigreeBuilder(demography).tables
    tables.sequence_length = sequence_length
    return tables


def add_pedigree_individual(
    tables, time, parents=(-1, -1), is_sample=None, population=0
):
    parents = np.sort(parents).astype("int32")
    ind_id = tables.individuals.add_row(parents=parents, flags=len(tables.individuals))
    flags = 0
    if is_sample is None:
        flags = tskit.NODE_IS_SAMPLE if time == 0 else 0
    elif is_sample:
        flags = tskit.NODE_IS_SAMPLE
    for _ in parents:
        tables.nodes.add_row(
            flags=flags,
            time=time,
            population=population,
            individual=ind_id,
        )
    return ind_id


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
    builder = msprime.PedigreeBuilder()

    time = num_generations - 1
    curr_gen = [builder.add_individual(time=time) for _ in range(num_founders)]
    for generation in range(1, num_generations):
        num_pairs = len(curr_gen) // 2
        if num_pairs == 0 and num_children_prob[0] != 1:
            raise Exception(
                f"Not enough people to make children in generation {generation}"
            )
        all_parents = rng.choice(curr_gen, size=(num_pairs, 2), replace=False)
        time -= 1
        curr_gen = []
        for parents in all_parents:
            num_children = rng.choice(len(num_children_prob), p=num_children_prob)
            for _ in range(num_children):
                parents = np.sort(parents).astype(np.int32)
                ind_id = builder.add_individual(time=time, parents=parents)
                curr_gen.append(ind_id)
    return builder.finalise(sequence_length)


class TestSimPedigree:
    """
    Tests for the simple Wright-Fisher pedigree simulator.
    """

    def verify(self, tables):
        tables.sequence_length = 10
        # Check that we can simulate with it.
        msprime.sim_ancestry(
            initial_state=tables, model="fixed_pedigree", random_seed=2
        )

    def test_zero_generations(self):
        tables = pedigrees.sim_pedigree(population_size=10, end_time=0, random_seed=1)
        assert len(tables.individuals) == 10
        # print(tables)
        assert np.all(tables.nodes.time) == 0
        assert np.all(tables.nodes.flags) == tskit.NODE_IS_SAMPLE
        assert np.all(tables.nodes.population) == 0
        self.verify(tables)

    @pytest.mark.parametrize("N", range(1, 5))
    @pytest.mark.parametrize("num_generations", range(1, 5))
    def test_n_generations(self, N, num_generations):
        tables = pedigrees.sim_pedigree(
            population_size=N, end_time=num_generations, random_seed=1
        )
        assert len(tables.individuals) == (num_generations + 1) * N

        tables.sequence_length = 1
        ts = tables.tree_sequence()

        t = num_generations
        # first N individuals are founders
        for ind in range(N):
            ind_obj = ts.individual(ind)
            assert list(ind_obj.parents) == [-1, -1]
            for node_id in ind_obj.nodes:
                node_obj = ts.node(node_id)
                assert node_obj.time == t

        n = N
        for _ in range(1, num_generations):
            t -= 1
            for ind in range(n, n + N):
                ind_obj = ts.individual(ind)
                for node_id in ind_obj.nodes:
                    node_obj = ts.node(node_id)
                    assert node_obj.time == t

                # parents are from the previous generation
                for parent in ind_obj.parents:
                    parent_obj = ts.individual(parent)
                    parent_node_obj = ts.node(parent_obj.nodes[0])
                    assert parent_node_obj.time == t + 1

            n += N
        # Last N nodes are all samples
        assert np.all(tables.nodes.flags[-2 * N :] == tskit.NODE_IS_SAMPLE)

        self.verify(tables)

    @pytest.mark.parametrize("num_generations", [1, 3, 7])
    def test_end_time_like_dtwf(self, num_generations):
        N = 10
        ped_tables = pedigrees.sim_pedigree(
            population_size=N, end_time=num_generations, random_seed=1
        )
        dtwf_ts = msprime.sim_ancestry(
            100,
            population_size=N,
            end_time=num_generations,
            random_seed=1,
            model="dtwf",
        )
        assert ped_tables.nodes.time[0] == num_generations
        assert dtwf_ts.tables.nodes.time[-1] == num_generations


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


def join_pedigrees(tables_list):
    """
    Join together pedigrees in the specified lists of tables and return
    the resulting pedigree.
    """
    tables = tables_list[0].copy()
    for other_tables in tables_list[1:]:
        num_individuals = len(tables.individuals)
        for node in other_tables.nodes:
            tables.nodes.append(
                node.replace(individual=node.individual + num_individuals)
            )
        for individual in other_tables.individuals:
            parents = individual.parents
            for j in range(len(parents)):
                if parents[j] != -1:
                    parents[j] += num_individuals
            tables.individuals.append(individual.replace(parents=parents))
    return tables.tree_sequence()


class TestJoinPedigrees:
    def test_two_pedigrees(self):
        tables1 = simulate_pedigree(
            num_founders=2,
            num_generations=2,
            sequence_length=100,
        )
        tables2 = simulate_pedigree(
            num_founders=5,
            num_generations=5,
            sequence_length=100,
        )
        joined = join_pedigrees([tables1, tables2])
        assert joined.num_nodes == len(tables1.nodes) + len(tables2.nodes)
        assert joined.num_individuals == len(tables1.individuals) + len(
            tables2.individuals
        )

    def test_three_pedigrees(self):

        tables1 = simulate_pedigree(
            num_founders=2,
            num_generations=2,
            sequence_length=100,
        )
        tables2 = simulate_pedigree(
            num_founders=5,
            num_generations=5,
            sequence_length=100,
        )
        tables3 = simulate_pedigree(
            num_founders=7,
            num_generations=1,
            sequence_length=100,
        )
        all_tables = [tables1, tables2, tables3]
        joined = join_pedigrees(all_tables)
        assert joined.num_nodes == sum(len(tables.nodes) for tables in all_tables)
        assert joined.num_individuals == sum(
            len(tables.individuals) for tables in all_tables
        )


class TestSimulateThroughPedigree:
    def verify(self, input_tables, recombination_rate=0):
        initial_state = input_tables.tree_sequence()

        # print()
        # print(list(input_tables.individuals.parents))
        # time = [
        #     input_tables.nodes.time[ind.nodes[0]]
        #     for ind in initial_state.individuals()
        # ]
        # print(time)

        # Using this low-level interface to make debugging easier. It's the
        # same effect as calling ts = msprime.sim_ancestry(...)
        sim = msprime.ancestry._parse_sim_ancestry(
            model="fixed_pedigree",
            initial_state=initial_state,
            recombination_rate=recombination_rate,
            random_seed=1,
        )
        sim.run()
        # print(sim)
        output_tables = tskit.TableCollection.fromdict(sim.tables.asdict())
        ts = output_tables.tree_sequence()
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
        founder_ids = {
            ind.id for ind in ts.individuals() if len(ancestors[ind.id]) == 0
        }

        for tree in ts.trees():
            for root in tree.roots:
                node = ts.node(root)
                # If this is a unary root it must be from a founder.
                if tree.num_children(root) == 1:
                    assert node.individual in founder_ids
            for node_id in tree.nodes():
                node = ts.node(node_id)
                individual = ts.individual(node.individual)
                if tree.parent(node_id) != tskit.NULL:
                    parent_node = ts.node(tree.parent(node_id))
                    assert parent_node.individual in ancestors[individual.id] | {
                        tskit.NULL
                    }
        return ts

    @pytest.mark.parametrize("num_founders", [2, 3, 5, 100])
    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_shallow(self, num_founders, recombination_rate):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0, 0, 1],
            num_generations=2,
            sequence_length=100,
        )
        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("num_founders", [2, 3, 10, 20])
    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_deep(self, num_founders, recombination_rate):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0, 0, 1],
            num_generations=6,  # Takes a long time if this is increased
            sequence_length=100,
        )
        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("num_founders", [2, 3])
    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_very_deep(self, num_founders, recombination_rate):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0, 0, 1],
            num_generations=16,  # Takes a long time if this is increased
            sequence_length=10,
        )
        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("num_founders", [2, 3, 10, 20])
    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_many_children(self, num_founders, recombination_rate):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0] * 99
            + [1],  # Each pair of parents will always have 100 children
            num_generations=2,
        )
        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("num_founders", [2, 3, 10, 20])
    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_unrelated(self, num_founders, recombination_rate):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0],
            num_generations=1,
        )
        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_two_pedigrees(self, recombination_rate):
        tables1 = simulate_pedigree(
            num_founders=5,
            num_generations=5,
            sequence_length=100,
        )
        tables2 = simulate_pedigree(
            num_founders=3,
            num_generations=8,
            sequence_length=100,
        )
        joined = join_pedigrees([tables1, tables2])
        ts = self.verify(joined.tables, recombination_rate)
        for tree in ts.trees():
            assert tree.num_roots >= 2

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_two_pedigrees_different_times(self, recombination_rate):
        tables1 = simulate_pedigree(
            num_founders=5,
            num_generations=5,
            sequence_length=100,
        )
        tables2 = simulate_pedigree(
            num_founders=3,
            num_generations=8,
            sequence_length=100,
        )
        # The pedigree in tables2 starts at 1 generation ago
        tables2.nodes.time += 1
        joined = join_pedigrees([tables1, tables2])
        ts = self.verify(joined.tables, recombination_rate)
        for tree in ts.trees():
            assert tree.num_roots >= 2
            for root in tree.roots:
                assert not tree.is_sample(root)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_two_pedigrees_non_overlapping_times(self, recombination_rate):
        tables1 = simulate_pedigree(
            num_founders=15,
            num_generations=3,
            sequence_length=100,
        )
        tables2 = simulate_pedigree(
            num_founders=6,
            num_generations=5,
            sequence_length=100,
        )
        # The pedigree in tables2 starts at 1 generation ago
        tables2.nodes.time += 5
        joined = join_pedigrees([tables1, tables2])
        joined.dump("joined_ped.ts")
        ts = self.verify(joined.tables, recombination_rate)
        for tree in ts.trees():
            assert tree.num_roots >= 2
            for root in tree.roots:
                assert not tree.is_sample(root)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_duo(self, recombination_rate):
        tables = get_base_tables(100)

        parent = add_pedigree_individual(tables, time=1)
        add_pedigree_individual(tables, time=0, parents=[parent, -1])

        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_trio(self, recombination_rate):
        tables = get_base_tables(100)

        gen0 = [add_pedigree_individual(tables, time=1) for _ in range(2)]
        add_pedigree_individual(tables, time=0, parents=gen0)

        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_mixed_generation_parents(self, recombination_rate):
        tables = get_base_tables(100)

        gen0 = [add_pedigree_individual(tables, time=2) for _ in range(2)]
        gen1 = [add_pedigree_individual(tables, time=1, parents=gen0)]
        add_pedigree_individual(tables, time=0, parents=(gen0[0], gen1[0]))

        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_inbreeding_sibs(self, recombination_rate):
        tables = get_base_tables(100)

        parents = [add_pedigree_individual(tables, time=2) for _ in range(2)]
        sibs = [
            add_pedigree_individual(tables, time=1, parents=parents) for _ in range(2)
        ]
        add_pedigree_individual(tables, time=0, parents=sibs)

        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_double_inbreeding_sibs(self, recombination_rate):
        tables = get_base_tables(100)

        parents = [add_pedigree_individual(tables, time=3) for _ in range(2)]
        sibs1 = [
            add_pedigree_individual(tables, time=2, parents=parents) for _ in range(2)
        ]
        sibs2 = [
            add_pedigree_individual(tables, time=1, parents=sibs1) for _ in range(2)
        ]
        add_pedigree_individual(tables, time=0, parents=sibs2)

        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_half_sibs(self, recombination_rate):
        tables = get_base_tables(100)

        parents = [add_pedigree_individual(tables, time=1) for _ in range(3)]
        add_pedigree_individual(tables, time=0, parents=parents[:2])
        add_pedigree_individual(tables, time=0, parents=parents[1:])

        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_inbreeding_half_sibs(self, recombination_rate):
        tables = get_base_tables(100)

        parents = [add_pedigree_individual(tables, time=2) for _ in range(3)]
        half_sib1 = add_pedigree_individual(tables, time=1, parents=parents[:2])
        half_sib2 = add_pedigree_individual(tables, time=1, parents=parents[1:])
        add_pedigree_individual(tables, time=0, parents=(half_sib1, half_sib2))

        self.verify(tables, recombination_rate)

    @pytest.mark.parametrize("recombination_rate", [0, 0.01])
    def test_oedipus(self, recombination_rate):
        tables = get_base_tables(100)

        parents = [add_pedigree_individual(tables, time=2) for _ in range(2)]
        offspring = add_pedigree_individual(tables, time=1, parents=parents)
        add_pedigree_individual(tables, time=0, parents=(parents[0], offspring))

        self.verify(tables, recombination_rate)

    def test_large_family(self):
        tables = get_base_tables(100)

        parents = [add_pedigree_individual(tables, time=1) for _ in range(2)]
        for _ in range(50):
            add_pedigree_individual(tables, time=0, parents=parents)
        ts = self.verify(tables)
        assert ts.num_trees == 1
        parent_nodes = {0, 1, 2, 3}
        assert set(ts.first().roots) <= parent_nodes

    def test_large_family_different_times(self):
        tables = get_base_tables(100)
        n = 20
        p1 = add_pedigree_individual(tables, time=n + 1)
        p2 = add_pedigree_individual(tables, time=n + 2)
        for j in range(n):
            add_pedigree_individual(tables, time=j, parents=[p1, p2], is_sample=True)
        ts = self.verify(tables)
        assert ts.num_trees == 1
        parent_nodes = {0, 1, 2, 3}
        assert set(ts.first().roots) <= parent_nodes

    def test_ancient_sample(self):
        tables = get_base_tables(100)
        parents = [add_pedigree_individual(tables, time=2) for _ in range(2)]
        add_pedigree_individual(tables, time=0, parents=parents)
        add_pedigree_individual(tables, time=1, parents=parents, is_sample=True)
        ts = self.verify(tables)
        assert ts.num_trees == 1
        parent_nodes = {0, 1, 2, 3}
        assert set(ts.first().roots) <= parent_nodes

    def test_no_time_zero_samples(self):
        tables = get_base_tables(100)
        parents = [add_pedigree_individual(tables, time=2) for _ in range(2)]
        add_pedigree_individual(tables, time=1, parents=parents, is_sample=True)
        add_pedigree_individual(tables, time=1, parents=parents, is_sample=True)
        ts = self.verify(tables)
        assert ts.num_trees == 1
        parent_nodes = {0, 1, 2, 3}
        assert set(ts.first().roots) <= parent_nodes

    def test_just_samples(self):
        tables = get_base_tables(100)
        add_pedigree_individual(tables, time=0, is_sample=True)
        add_pedigree_individual(tables, time=1, is_sample=True)
        ts = self.verify(tables)
        assert ts.num_trees == 1
        parent_nodes = {0, 1, 2, 3}
        assert set(ts.first().roots) == parent_nodes


class TestSimulateThroughPedigreeEventByEvent(TestSimulateThroughPedigree):
    def verify(self, input_tables, recombination_rate=0):
        ts1 = msprime.sim_ancestry(
            model="fixed_pedigree",
            initial_state=input_tables,
            recombination_rate=recombination_rate,
            random_seed=1,
        )
        sim = msprime.ancestry._parse_sim_ancestry(
            model="fixed_pedigree",
            initial_state=input_tables,
            recombination_rate=recombination_rate,
            random_seed=1,
        )
        sim.run(event_chunk=1)
        output_tables = tskit.TableCollection.fromdict(sim.tables.asdict())
        output_tables.assert_equals(ts1.tables, ignore_provenance=True)
        return ts1


class TestContinueSimulateThroughPedigree(TestSimulateThroughPedigree):
    """
    Can we extend the simulations of the pedigree as we'd expect?
    """

    def verify(self, input_tables, recombination_rate=0):

        ts1 = msprime.sim_ancestry(
            model="fixed_pedigree",
            initial_state=input_tables,
            recombination_rate=recombination_rate,
            random_seed=42,
        )
        # print(ts1.draw_text())
        ts2 = msprime.sim_ancestry(
            initial_state=ts1, population_size=1, random_seed=1234
        )
        for tree in ts2.trees():
            assert tree.num_roots == 1
            for node_id in tree.nodes():
                if tree.num_children(node_id) == 1:
                    # Any unary nodes should be associated with a pedigree
                    # founder.
                    node = ts2.node(node_id)
                    assert node.individual != tskit.NULL
                    individual = ts2.individual(node.individual)
                    assert list(individual.parents) == [tskit.NULL, tskit.NULL]
        tables1 = ts1.tables
        tables2 = ts2.tables
        tables1.individuals.assert_equals(
            tables2.individuals[: len(tables1.individuals)]
        )
        tables1.nodes.assert_equals(tables2.nodes[: len(tables1.nodes)])

        return ts1

    def test_family_different_times(self):
        # Check that we're picking up ancient samples correctly from the
        # pedigree simulation by having two different parents at very
        # different times, and a small population size. This should guarantee
        # that the nodes for the first parent coalesces first.
        #
        # 28.00┊   10    ┊
        #      ┊  ┏━┻━┓  ┊
        # 23.00┊  ┃   9  ┊
        #      ┊  ┃  ┏┻┓ ┊
        # 20.00┊  ┃  3 2 ┊
        #      ┊  ┃  ┃ ┃ ┊
        # 4.00 ┊  8  ┃ ┃ ┊
        #      ┊ ┏┻┓ ┃ ┃ ┊
        # 2.00 ┊ 0 1 ┃ ┃ ┊
        #      ┊ ┃ ┃ ┃ ┃ ┊
        # 1.00 ┊ ┃ 6 ┃ 7 ┊
        #      ┊ ┃   ┃   ┊
        # 0.00 ┊ 4   5   ┊
        #      0.00     100.00
        tables = get_base_tables(100)
        n = 2
        p1 = add_pedigree_individual(tables, time=n)
        p2 = add_pedigree_individual(tables, time=10 * n)
        for j in range(n):
            add_pedigree_individual(tables, time=j, parents=[p1, p2], is_sample=True)
        ts1 = self.verify(tables)
        ts2 = msprime.sim_ancestry(
            initial_state=ts1, population_size=5, model="dtwf", random_seed=1234
        )
        assert ts2.num_trees == 1
        tree = ts2.first()
        # The youngest parent nodes should coalesce
        assert tree.parent(0) == tree.parent(1)


class TestSimulateThroughPedigreeReplicates(TestSimulateThroughPedigree):
    def verify(self, input_tables, recombination_rate=0):
        num_replicates = 5
        replicates = list(
            msprime.sim_ancestry(
                model="fixed_pedigree",
                initial_state=input_tables,
                recombination_rate=recombination_rate,
                random_seed=42,
                num_replicates=num_replicates,
            )
        )

        ts1 = msprime.sim_ancestry(
            model="fixed_pedigree",
            initial_state=input_tables,
            recombination_rate=recombination_rate,
            random_seed=42,
        )
        ts1.tables.assert_equals(replicates[0].tables, ignore_provenance=True)

        assert len(replicates) == num_replicates
        for ts in replicates:
            output_tables = ts.tables
            input_tables.individuals.assert_equals(output_tables.individuals)
            input_tables.nodes.assert_equals(output_tables.nodes)
        return replicates[0]


class TestSimulateThroughPartialPedigree:
    def test_one_parent_removed_each_individual(self):
        tables = simulate_pedigree(
            num_founders=5,
            num_generations=5,
            sequence_length=100,
        )
        parents = tables.individuals.parents
        parents[0::2] = -1
        tables.individuals.parents = parents

        ts = msprime.sim_ancestry(
            model="fixed_pedigree",
            initial_state=tables,
            random_seed=42,
        )
        # We quickly reach deadends in the pedigree, just check that
        # we've done something reasonably sensible.
        assert ts.num_edges > 0

    def test_10_percent_parents_removed(self):
        tables = simulate_pedigree(
            num_founders=5,
            num_generations=10,
            sequence_length=100,
        )
        parents = tables.individuals.parents
        parents[0::10] = -1
        tables.individuals.parents = parents

        ts = msprime.sim_ancestry(
            model="fixed_pedigree",
            initial_state=tables,
            random_seed=42,
        )
        # We quickly reach deadends in the pedigree, just check that
        # we've done something reasonably sensible.
        assert ts.num_edges > 0


class TestSimulateThroughPedigreeErrors:
    def test_parents_are_samples(self):
        tables = get_base_tables(100)
        parents = [
            add_pedigree_individual(tables, time=1, is_sample=True) for _ in range(2)
        ]
        add_pedigree_individual(tables, parents=parents, time=0)

        with pytest.raises(_msprime.LibraryError, match="1855"):
            msprime.sim_ancestry(
                initial_state=tables,
                model="fixed_pedigree",
            )

    @pytest.mark.parametrize("num_parents", [0, 1, 3])
    def test_not_two_parents(self, num_parents):
        tables = get_base_tables(100)
        parents = [add_pedigree_individual(tables, time=1) for _ in range(num_parents)]
        add_pedigree_individual(tables, parents=parents, time=0)
        with pytest.raises(_msprime.InputError, match="exactly two parents"):
            msprime.sim_ancestry(initial_state=tables, model="fixed_pedigree")

    @pytest.mark.parametrize("num_nodes", [0, 1, 3])
    def test_not_two_nodes(self, num_nodes):
        tables = get_base_tables(100)
        ind = tables.individuals.add_row(parents=[-1, -1])
        for _ in range(num_nodes):
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind, population=0
            )
        with pytest.raises(_msprime.InputError, match="exactly two nodes"):
            msprime.sim_ancestry(initial_state=tables, model="fixed_pedigree")

    def test_node_times_disagree(self):
        tables = get_base_tables(100)
        ind = tables.individuals.add_row(parents=[-1, -1])
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind, population=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, time=1, individual=ind, population=0
        )
        with pytest.raises(_msprime.InputError, match="times for the two nodes"):
            msprime.sim_ancestry(initial_state=tables, model="fixed_pedigree")

    def test_node_populations_disagree(self):
        tables = get_base_tables(100)
        tables.populations.add_row({"name": "pop_1", "description": "Y"})
        ind = tables.individuals.add_row(parents=[-1, -1])
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind, population=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind, population=1
        )
        with pytest.raises(_msprime.InputError, match="populations for the two nodes"):
            msprime.sim_ancestry(initial_state=tables, model="fixed_pedigree")

    def test_no_samples(self):
        tables = get_base_tables(100)
        ind = tables.individuals.add_row(parents=[-1, -1])
        tables.nodes.add_row(flags=0, time=0, individual=ind, population=0)
        tables.nodes.add_row(flags=0, time=0, individual=ind, population=0)
        with pytest.raises(_msprime.InputError, match="samples"):
            msprime.sim_ancestry(initial_state=tables, model="fixed_pedigree")

    def test_pedigree_time_travel(self):
        tables = get_base_tables(100)
        parents = [add_pedigree_individual(tables, time=0) for _ in range(2)]
        add_pedigree_individual(tables, time=1, parents=parents, is_sample=True)
        with pytest.raises(_msprime.InputError, match="time for a parent must be"):
            msprime.sim_ancestry(initial_state=tables, model="fixed_pedigree")

    @pytest.mark.parametrize("num_founders", [2, 3, 10, 20])
    def test_all_samples(self, num_founders):
        tables = simulate_pedigree(
            num_founders=num_founders,
            num_children_prob=[0, 0, 1],
            num_generations=6,
            sequence_length=100,
        )
        flags = tables.nodes.flags
        flags[:] = tskit.NODE_IS_SAMPLE
        tables.nodes.flags = flags
        with pytest.raises(_msprime.LibraryError, match="1855"):
            msprime.sim_ancestry(initial_state=tables, model="fixed_pedigree")


class TestSimulateThroughPedigreeMultiplePops:
    def test_demography_raises_error(self):
        demography = msprime.Demography.isolated_model([10, 10])
        pb = msprime.PedigreeBuilder(demography)
        pb.add_individual(time=0, population=0)
        tables = pb.finalise(1)
        with pytest.raises(ValueError, match="Cannot specify demography"):
            msprime.sim_ancestry(
                initial_state=tables,
                demography=demography,
                model="fixed_pedigree",
                random_seed=1,
            )

    def test_trio_parents_different_pops(self):
        demography = msprime.Demography.isolated_model([10, 10])
        pb = msprime.PedigreeBuilder(demography)
        parents = [pb.add_individual(time=2, population=j) for j in range(2)]
        pb.add_individual(time=0, parents=parents, population=0)
        pedigree = pb.finalise(1)
        ts = msprime.sim_ancestry(
            initial_state=pedigree, model="fixed_pedigree", random_seed=1
        )
        # Founders are in different populations so cannot coalesce
        with pytest.raises(_msprime.LibraryError, match="Infinite waiting"):
            msprime.sim_ancestry(initial_state=ts, demography=demography, random_seed=2)

    def test_trio_parents_same_pop(self):
        demography = msprime.Demography.isolated_model([10, 10])
        pb = msprime.PedigreeBuilder(demography)
        parents = [pb.add_individual(time=2, population=1) for j in range(2)]
        pb.add_individual(time=0, parents=parents, population=0)
        pedigree = pb.finalise(1)

        ts = msprime.sim_ancestry(
            initial_state=pedigree, model="fixed_pedigree", random_seed=1
        )
        ts = msprime.sim_ancestry(
            initial_state=ts, demography=demography, random_seed=2
        )
        for j in range(4):
            assert ts.node(j).population == 1
        assert ts.node(4).population == 0
        assert ts.node(5).population == 0
        mrca = ts.node(6)
        assert mrca.individual == -1
        assert mrca.population == 1

    def test_trio_parents_different_pops_with_split(self):
        demography = msprime.Demography.isolated_model([10, 10, 10])
        demography.add_population_split(time=10, derived=[0, 1], ancestral=2)
        pb = msprime.PedigreeBuilder(demography)
        parents = [pb.add_individual(time=2, population=j) for j in range(2)]
        pb.add_individual(time=0, parents=parents, population=0)
        pedigree = pb.finalise(1)

        ts = msprime.sim_ancestry(
            initial_state=pedigree, model="fixed_pedigree", random_seed=1
        )
        # Should be able to coalesce due to the added population split
        ts = msprime.sim_ancestry(
            initial_state=ts, demography=demography, random_seed=2
        )
        assert ts.node(6).time > 2
        assert ts.node(6).population == 2

    def test_trio_parents_different_pops_with_migrations(self):
        demography = msprime.Demography.isolated_model([10, 10])
        demography.add_symmetric_migration_rate_change(
            time=0, populations=[0, 1], rate=0.1
        )
        pb = msprime.PedigreeBuilder(demography)
        parents = [pb.add_individual(time=2, population=j) for j in range(2)]
        pb.add_individual(time=0, parents=parents, population=0)
        pedigree = pb.finalise(1)

        ts = msprime.sim_ancestry(
            initial_state=pedigree, model="fixed_pedigree", random_seed=1
        )
        # Should be able to coalesce due to symmetric migration
        ts = msprime.sim_ancestry(
            initial_state=ts, demography=demography, random_seed=2
        )
        assert ts.node(6).time > 2


@dataclasses.dataclass
class FamEntry:
    fid: str = "0"
    iid: str = "0"
    pat: str = "0"
    mat: str = "0"
    sex: str = "0"
    phen: str = None

    def get_row(self, delimiter="\t"):
        return delimiter.join([x for x in asdict(self).values() if x is not None])


class TestParseFam:
    """
    Tests for the parse_fam function.
    """

    def get_parsed_fam(self, entries, delimiter="\t"):
        content = "\n".join([entry.get_row(delimiter=delimiter) for entry in entries])
        with tempfile.TemporaryFile() as f:
            f.write(bytes(content, "utf-8"))
            f.seek(0)
            return pedigrees.parse_fam(f)

    def test_empty_file(self):
        entries = []
        with pytest.warns(UserWarning):
            tb = self.get_parsed_fam(entries=entries)
        assert len(tb) == 0

    @pytest.mark.parametrize("iid", ["1", "a", "100", "abc"])
    def test_single_line(self, iid):
        entries = [FamEntry(iid=iid)]
        tb = self.get_parsed_fam(entries=entries)
        assert len(tb) == 1
        assert np.array_equal(tb[0].parents, [-1, -1])
        assert tb[0].metadata["plink_fid"] == "0"
        assert tb[0].metadata["plink_iid"] == str(iid)
        assert tb[0].metadata["sex"] == 0

    @pytest.mark.parametrize("iids", [("1", "2"), ("a", "b")])
    def test_multiple_line_file(self, iids):
        # test both integer and string IIDs
        iid1, iid2 = iids
        entries = [FamEntry(iid=iid1), FamEntry(iid=iid2)]
        tb = self.get_parsed_fam(entries=entries)
        assert len(tb) == 2
        for idx in range(2):
            assert np.array_equal(tb[idx].parents, [-1, -1])
            assert tb[idx].metadata["plink_fid"] == "0"
            assert tb[idx].metadata["plink_iid"] == str(entries[idx].iid)
            assert tb[idx].metadata["sex"] == 0

    @pytest.mark.parametrize("n_cols", range(1, 5))
    def test_insufficient_cols(self, n_cols):
        fields = list(asdict(FamEntry()))
        entry = FamEntry(iid="1")
        for field in fields[n_cols:]:
            entry.__setattr__(field, None)
        with pytest.raises(IndexError):
            self.get_parsed_fam(entries=[entry])

    def test_unrelated_duplicate_iids(self):
        # Individuals have the same IID, but are in different families
        entries = [FamEntry(iid="1"), FamEntry(fid="1", iid="1")]
        tb = self.get_parsed_fam(entries=entries)
        assert len(tb) == 2
        assert tb[0].metadata["plink_fid"] == "0"
        assert tb[0].metadata["plink_iid"] == "1"
        assert tb[1].metadata["plink_fid"] == "1"
        assert tb[1].metadata["plink_iid"] == "1"

    def test_duplicate_rows(self):
        entries = [FamEntry(iid="1"), FamEntry(iid="1")]
        with pytest.raises(ValueError):
            self.get_parsed_fam(entries=entries)

    def test_space_delimited(self):
        entries = [FamEntry(iid="1")]
        tb = self.get_parsed_fam(entries=entries, delimiter=" ")
        assert np.array_equal(tb[0].parents, [-1, -1])
        assert tb[0].metadata["plink_fid"] == "0"
        assert tb[0].metadata["plink_iid"] == "1"
        assert tb[0].metadata["sex"] == 0

    def test_missing_phen_col(self):
        entries = [FamEntry(iid="1", phen="1")]
        tb = self.get_parsed_fam(entries=entries)

        entries = [FamEntry(iid="1")]  # remove last column (PHEN column)
        tb_missing = self.get_parsed_fam(entries=entries)

        assert tb == tb_missing

    @pytest.mark.parametrize("sex", [-2, 3, "F"])
    def test_bad_sex_value(self, sex):
        entries = [FamEntry(iid="1", sex=str(sex))]
        with pytest.raises(ValueError):
            self.get_parsed_fam(entries=entries)

    def test_empty_sex_value(self):
        entries = [FamEntry(iid="1", sex="")]
        with pytest.raises(IndexError):
            self.get_parsed_fam(entries=entries)

    def test_single_family_map_parent_ids(self):
        # PAT is mapped if the individual exists in the dataset
        entries = [FamEntry(iid="1"), FamEntry(iid="2", pat="1")]
        tb = self.get_parsed_fam(entries=entries)
        assert np.array_equal(tb[1].parents, [0, -1])

        # MAT is mapped if the individual exists in the dataset
        entries = [FamEntry(iid="1"), FamEntry(iid="2", mat="1")]
        tb = self.get_parsed_fam(entries=entries)
        assert np.array_equal(tb[1].parents, [-1, 0])

        # both parent IDs are remapped if the both parents exist in the dataset
        entries = [
            FamEntry(iid="1"),
            FamEntry(iid="2"),
            FamEntry(iid="3", pat="1", mat="2"),
        ]
        tb = self.get_parsed_fam(entries=entries)
        assert np.array_equal(tb[2].parents, [0, 1])

    def test_missing_parent_id(self):
        # KeyError raised if at least one parent (PAT) does not exist in dataset
        entries = [
            FamEntry(iid="2"),
            FamEntry(iid="3", pat="1", mat="2"),
        ]
        with pytest.raises(KeyError):
            self.get_parsed_fam(entries=entries)

        # KeyError raised if at least one parent (MAT) does not exist in dataset
        entries = [
            FamEntry(iid="1"),
            FamEntry(iid="3", pat="1", mat="2"),
        ]
        with pytest.raises(KeyError):
            self.get_parsed_fam(entries=entries)

        # KeyError raised if both parents do not exist in dataset
        entries = [FamEntry(iid="1", pat="2", mat="3")]
        with pytest.raises(KeyError):
            self.get_parsed_fam(entries=entries)

    def test_multiple_family_map_parent_ids(self):
        # parents mapped correctly when the same parent ID is used in different families
        entries = [
            FamEntry(iid="2"),
            FamEntry(iid="1"),
            FamEntry(fid="1", iid="2"),
            FamEntry(fid="1", iid="1"),
            FamEntry(iid="3", pat="1", mat="2"),
            FamEntry(fid="1", iid="3", pat="1", mat="2"),
        ]
        tb = self.get_parsed_fam(entries=entries)
        for idx in range(4):
            assert np.array_equal(tb[idx].parents, [-1, -1])
        assert np.array_equal(tb[4].parents, [1, 0])
        assert np.array_equal(tb[5].parents, [3, 2])

        # KeyError raised when FID does not match, even if parent ID matches
        entries = [
            FamEntry(iid="2"),
            FamEntry(iid="1"),
            FamEntry(iid="3", pat="1", mat="2"),
            FamEntry(
                fid="1", iid="1", pat="2", mat="3"
            ),  # there is no parent with FID=1, IID=3
            FamEntry(fid="1", iid="2"),
        ]
        with pytest.raises(KeyError):
            self.get_parsed_fam(entries)

    def test_grandparents(self):
        entries = [
            FamEntry(iid="4"),
            FamEntry(iid="3"),
            FamEntry(iid="2"),
            FamEntry(iid="1"),
            FamEntry(iid="6", pat="3", mat="4"),
            FamEntry(iid="5", pat="1", mat="2"),
            FamEntry(iid="7", pat="5", mat="6"),
        ]
        tb = self.get_parsed_fam(entries=entries)
        assert np.array_equal(tb[4].parents, [1, 0])
        assert np.array_equal(tb[5].parents, [3, 2])
        assert np.array_equal(tb[6].parents, [5, 4])

    def test_children_before_parents(self, tmp_path):
        entries = [
            FamEntry(iid="1", pat="2", mat="3"),
            FamEntry(iid="2"),
            FamEntry(iid="3"),
        ]
        content = "\n".join([entry.get_row() for entry in entries])
        fam_path = f"{tmp_path}/test.fam"
        with open(fam_path, "w+") as f:
            f.write(content)
            f.seek(0)
            tb = pedigrees.parse_fam(f)

        tc = tskit.TableCollection(1)
        # Issue 1489 will make this better
        tc.individuals.metadata_schema = tb.metadata_schema
        for row in tb:
            tc.individuals.append(row)
        tc.tree_sequence()  # creating tree sequence should succeed

        assert np.array_equal(tb[0].parents, [1, 2])
        assert np.array_equal(tb[1].parents, [-1, -1])
        assert np.array_equal(tb[2].parents, [-1, -1])


class TestPedigreeBuilder:
    def test_is_sample_default(self):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0)
        pb.add_individual(time=0.0001)
        pb.add_individual(time=1)
        tables = pb.finalise(1)
        zero_nodes = tables.nodes.time == 0
        assert np.all(tables.nodes.flags[zero_nodes] == tskit.NODE_IS_SAMPLE)
        non_zero_nodes = tables.nodes.time != 0
        assert np.all(tables.nodes.flags[non_zero_nodes] == 0)

    def test_is_sample_default_override(self):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0, is_sample=True)
        pb.add_individual(time=0.0001, is_sample=True)
        pb.add_individual(time=1, is_sample=True)
        tables = pb.finalise(1)
        assert np.all(tables.nodes.flags == tskit.NODE_IS_SAMPLE)

    @pytest.mark.parametrize("parents", [[], [0], [0, 1, 2]])
    def test_bad_parents_length(self, parents):
        pb = pedigrees.PedigreeBuilder()
        with pytest.raises(ValueError, match="exactly two parents"):
            pb.add_individual(time=0, parents=parents)

    def test_default_parents(self):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0)
        tables = pb.finalise(1)
        np.testing.assert_array_equal(tables.individuals[0].parents, [-1, -1])

    @pytest.mark.parametrize(
        "parents", [[0, 1], [0.0, 1.0], np.array([0, 1], dtype=np.int8)]
    )
    def test_parents(self, parents):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=1)
        pb.add_individual(time=1)
        pb.add_individual(time=0, parents=parents)
        tables = pb.finalise(1)
        np.testing.assert_array_equal(tables.individuals[2].parents, parents)

    def test_default_demography(self):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0)
        ts = pb.finalise(1).tree_sequence()
        assert ts.num_populations == 1
        assert ts.population(0).metadata["name"] == "pop_0"
        assert all(node.population == 0 for node in ts.nodes())

    def test_one_pop_demography(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=1)
        pb = pedigrees.PedigreeBuilder(demography=demography)
        pb.add_individual(time=0)
        ts = pb.finalise(1).tree_sequence()
        assert ts.num_populations == 1
        assert ts.population(0).metadata["name"] == "A"
        assert all(node.population == 0 for node in ts.nodes())

    def test_two_pop_demography(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=1)
        demography.add_population(name="B", initial_size=1)
        pb = pedigrees.PedigreeBuilder(demography=demography)
        ts = pb.finalise(1).tree_sequence()
        assert ts.num_populations == 2
        assert ts.population(0).metadata["name"] == "A"
        assert ts.population(1).metadata["name"] == "B"

    def test_two_pop_demography_population_required(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=1)
        demography.add_population(name="B", initial_size=1)
        pb = pedigrees.PedigreeBuilder(demography=demography)
        with pytest.raises(ValueError, match="multi-population demography"):
            pb.add_individual(time=0)

    def test_two_pop_demography_populations_by_index(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=1)
        demography.add_population(name="B", initial_size=1)
        pb = pedigrees.PedigreeBuilder(demography=demography)
        pb.add_individual(time=0, population=0)
        pb.add_individual(time=0, population=1)
        ts = pb.finalise(1).tree_sequence()
        assert all(ts.node(u).population == 0 for u in [0, 1])
        assert all(ts.node(u).population == 1 for u in [2, 3])

    def test_two_pop_demography_populations_by_name(self):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=1)
        demography.add_population(name="B", initial_size=1)
        pb = pedigrees.PedigreeBuilder(demography=demography)
        pb.add_individual(time=0, population="A")
        pb.add_individual(time=0, population="B")
        ts = pb.finalise(1).tree_sequence()
        assert all(ts.node(u).population == 0 for u in [0, 1])
        assert all(ts.node(u).population == 1 for u in [2, 3])

    def test_node_individual_linkage(self):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0)
        pb.add_individual(time=0)
        tables = pb.finalise(1)
        assert len(tables.individuals) == 2
        assert len(tables.nodes) == 4
        assert tables.nodes.individual[0] == 0
        assert tables.nodes.individual[1] == 0
        assert tables.nodes.individual[2] == 1
        assert tables.nodes.individual[3] == 1

    @pytest.mark.parametrize("sequence_length", [0.1, 1, 10.11])
    def test_finalise_sequence_length(self, sequence_length):
        pb = pedigrees.PedigreeBuilder()
        tables = pb.finalise(sequence_length)
        assert tables.sequence_length == sequence_length

    def test_finalise_copy(self):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0)
        t1 = pb.finalise(1)
        t2 = pb.finalise(1)
        assert t1 == t2
        assert t1 is not t2

    @pytest.mark.parametrize("time", [0, 1, 10.11, np.array([12345.6])[0]])
    def test_time(self, time):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=time)
        t1 = pb.finalise(1)
        assert np.all(t1.nodes.time == time)

    @pytest.mark.parametrize(
        "metadata", [{}, {"A": "B"}, {"A": [0, 1, 2]}, {"A": {"B": "C"}}]
    )
    def test_metadata(self, metadata):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0, metadata=metadata)
        ts = pb.finalise(1).tree_sequence()
        assert ts.individual(0).metadata == metadata

    def test_add_individual_return_value(self):
        pb = pedigrees.PedigreeBuilder()
        n = 5
        for j in range(n):
            ind_id = pb.add_individual(time=0, metadata={"id": j})
            assert ind_id == j
        ts = pb.finalise(1).tree_sequence()
        # double check that these are the IDs in the tables
        for j in range(n):
            assert ts.individual(j).metadata == {"id": j}

    def test_indexed(self):
        pb = pedigrees.PedigreeBuilder()
        pb.add_individual(time=0)
        t1 = pb.finalise(1)
        assert t1.has_index()


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
