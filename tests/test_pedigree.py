import os
import shutil
import tempfile
import unittest

import numpy as np
import pytest

import msprime


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

    def test_pedigree_replicates(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        time = np.array([0, 0, 1, 1])
        is_sample = np.array([1, 1, 0, 0])

        model = msprime.WrightFisherPedigree()
        ped = msprime.Pedigree(individual, parents, time, is_sample, sex=None, ploidy=2)
        replicates = msprime.simulate(
            2, pedigree=ped, model=model, recombination_rate=1, num_replicates=100
        )
        for ts in replicates:
            assert ts is not None

    def test_pedigree_samples(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        time = np.array([0, 0, 1, 1])

        ped = msprime.Pedigree(individual, parents, time)
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

        ped = msprime.Pedigree(individual, parents, time, is_sample, sex=None, ploidy=2)

        ped.save_txt(self.temp_pedigree_text_file)
        with pytest.raises(NotImplementedError):
            msprime.Pedigree.read_txt(
                self.temp_pedigree_text_file,
                sex_col=4,
            )
        # FIXME
        # The compute_times should be done automatically in this case .
        # ped_from_txt = msprime.Pedigree.read_txt(
        #     self.temp_pedigree_text_file, time_col=None
        # )
        # ts = msprime.simulate(2, pedigree=ped_from_txt, model="wf_ped")
        # self.assertTrue(ts is not None)
        # ped_from_txt = msprime.Pedigree.read_txt(
        #     self.temp_pedigree_text_file, time_col=3
        # )
        # ts = msprime.simulate(2, pedigree=ped_from_txt, model="wf_ped")
        # self.assertTrue(ts is not None)

        ped.save_npy(self.temp_pedigree_array_file)
        ped_from_npy = msprime.Pedigree.read_npy(self.temp_pedigree_array_file)
        # TODO compre this to the file above.
        assert isinstance(ped_from_npy, msprime.Pedigree)

    def test_pedigree_times(self):
        individual = np.array([1, 2, 3, 4])
        time = np.array([0, 0, 1, 1])

        parent_IDs = np.array([3, 4, 3, 4, 0, 0, 0, 0]).reshape(-1, 2)
        estimated_times = msprime.Pedigree.get_times(
            individual, parent_IDs=parent_IDs, check=True
        )
        assert (time == estimated_times).all()

    def test_pedigree_utils(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        parent_IDs = np.array([3, 4, 3, 4, 0, 0, 0, 0]).reshape(-1, 2)

        bad_individual = np.array([0, 1, 2, 3])
        with pytest.raises(ValueError):
            msprime.Pedigree.parent_ID_to_index(bad_individual, parent_IDs)
        assert (
            msprime.Pedigree.parent_ID_to_index(individual, parent_IDs) == parents
        ).all()
        assert (
            msprime.Pedigree.parent_index_to_ID(individual, parents) == parent_IDs
        ).all()

    def test_pedigree_sanity_checks(self):
        individual = np.array([1, 2, 3, 4])
        parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        time = np.array([0, 0, 1, 1])
        is_sample = np.array([1, 1, 0, 0])

        with pytest.raises(NotImplementedError):
            msprime.Pedigree(
                individual=individual,
                parents=parents,
                time=time,
                ploidy=1,
            )
        bad_parents = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 4)
        with pytest.raises(ValueError):
            msprime.Pedigree(
                individual=individual,
                parents=bad_parents,
                time=time,
            )
        bad_individual = np.array([-1, 2, 3, 4])
        with pytest.raises(ValueError):
            msprime.Pedigree(
                individual=bad_individual,
                parents=parents,
                time=time,
            )

        bad_times = np.array([1, 1, 1, 1])
        ped = msprime.Pedigree(
            individual, parents, bad_times, is_sample, sex=None, ploidy=2
        )
        with pytest.raises(ValueError):
            ped.check_times(
                individual=ped.individual,
                parents=ped.parents,
                time=ped.time,
            )

        ped = msprime.Pedigree(individual, parents, time, sex=None, ploidy=2)
        with pytest.raises(ValueError):
            ped.set_samples()
        with pytest.raises(ValueError):
            ped.set_samples(num_samples=10)
        with pytest.raises(ValueError):
            ped.set_samples(num_samples=2, sample_IDs=[1, 2])

        with pytest.raises(ValueError):
            ped.get_times(individual=ped.individual)
