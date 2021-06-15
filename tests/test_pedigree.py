import os
import shutil
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
    tc: tskit.TableCollection,
    n_founders=8,
    n_children_prob=None,
    n_generations=3,
    random_sex=True,
):
    def choose_sex(random=True, sex=None):
        # male = 1, female = 2
        if random:
            sex = np.random.binomial(n=1, p=0.5) + 1
        else:  # non-random: alternate male female to ensure equal balance
            sex = (sex % 2) + 1
        return sex

    def choose_n_children(n_children_prob, n_families=1):
        return np.random.choice(
            a=len(n_children_prob), size=n_families, p=n_children_prob
        )

    tb = tc.individuals
    tb.metadata_schema = tskit.MetadataSchema(
        {
            "codec": "json",
            "type": "object",
            "properties": {
                "sex": {"type": "integer"},
            },
            "required": ["sex"],
            "additionalProperties": True,
        }
    )

    if n_children_prob is None:
        n_children_prob = [0, 0, 1]

    if random_sex:
        sex = None
    else:
        sex = 1  # male=1, female=2

    # first generation
    n_males = int(n_founders / 2)
    curr_gen = [range(n_males), range(n_males, n_founders)]
    for ind in range(n_founders):
        tb.add_row(parents=[-1, -1], metadata={"sex": (ind >= n_males) + 1})

    # next generations
    for gen_idx in range(2, n_generations + 1):
        next_gen = [[], []]
        avail_pat = np.random.permutation(curr_gen[0])
        avail_mat = np.random.permutation(curr_gen[1])
        n_pairs = min(len(curr_gen[0]), len(curr_gen[1]))
        if n_pairs == 0 and n_children_prob[0] != 1:
            raise Exception(
                f"Not enough people to make children in generation {gen_idx}"
            )
        pairs = zip(avail_pat[:n_pairs], avail_mat[:n_pairs])
        n_children_per_pair = choose_n_children(
            n_children_prob=n_children_prob, n_families=n_pairs
        )
        for (pat, mat), n_children in zip(pairs, n_children_per_pair):
            for _ in range(n_children):
                sex = choose_sex(random=random_sex, sex=sex)
                next_gen[sex - 1] += [len(tb)]
                tb.add_row(parents=[pat, mat], metadata={"sex": sex})
        curr_gen = next_gen


class TestPedigreeSimulation:
    """
    Tests for simple pedigree simulator
    """

    def simple_sim(
        self, n_founders=2, n_children_prob=None, n_generations=2, random_sex=True
    ):
        tc = tskit.TableCollection(0)
        if n_children_prob is None:
            n_children_prob = [0, 0, 1]
        simulate_pedigree(
            tc=tc,
            n_founders=n_founders,
            n_children_prob=n_children_prob,
            n_generations=n_generations,
            random_sex=random_sex,
        )
        return tc.individuals

    def test_one_generation_no_children(self):
        n_founders = 8
        tb = self.simple_sim(n_founders=n_founders, n_generations=1)
        assert len(tb) == n_founders

        # check that all parents of founders are missing
        assert all([np.array_equal(row.parents, [-1, -1]) for row in tb])

        n_males = len([row for row in tb if row.metadata["sex"] == 1])
        n_females = len([row for row in tb if row.metadata["sex"] == 2])
        assert n_males == n_females == n_founders / 2

    def test_one_trio(self):
        tb = self.simple_sim(n_founders=2, n_children_prob=[0, 1], n_generations=2)
        assert len(tb) == 3
        assert np.array_equal(tb[2].parents, [0, 1])
        assert all([np.array_equal(row.parents, [-1, -1]) for row in tb[:2]])

    def test_grandparents(self):
        tb = self.simple_sim(
            n_founders=4, n_children_prob=[0, 1], n_generations=3, random_sex=False
        )
        assert len(tb) == 7
        assert {row.parents[0] for row in tb[4:6]} == {0, 1}
        assert {row.parents[1] for row in tb[4:6]} == {2, 3}
        assert all([np.array_equal(row.parents, [-1, -1]) for row in tb[:4]])

    def test_insufficient_founders(self):
        with pytest.raises(Exception):
            self.simple_sim(n_founders=1, n_children_prob=[0, 1])
        with pytest.raises(Exception):
            self.simple_sim(n_founders=3, n_children_prob=[0, 1], n_generations=3)

    @pytest.mark.parametrize("n_children", range(5))
    def test_nonrandom_child_prob(self, n_children):
        tb = self.simple_sim(
            n_founders=2, n_children_prob=[0] * n_children + [1], n_generations=2
        )
        assert len(tb) == 2 + n_children

    @pytest.mark.parametrize("n_children_prob", [[0.5, 0.5], [0, 0.5, 0.5]])
    def test_expected_n_children(self, n_children_prob):
        n_reps = 100
        n_founders = 2
        n_children = []
        for _ in range(n_reps):
            tb = self.simple_sim(
                n_founders=n_founders, n_children_prob=n_children_prob, n_generations=2
            )
            n_children.append(len(tb) - n_founders)
        expected_n_children = sum(
            prob * n for prob, n in zip(n_children_prob, range(len(n_children_prob)))
        )
        assert np.mean(n_children) == pytest.approx(expected_n_children, rel=0.1)

    def test_bad_n_children_prob(self):
        with pytest.raises(ValueError):
            self.simple_sim(n_children_prob=[2])
        with pytest.raises(ValueError):
            self.simple_sim(n_children_prob=[1, 1])
