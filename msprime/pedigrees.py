#
# Copyright (C) 2019-2020 University of Oxford
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
Module responsible managing pedigree definitions.
"""
import collections
import inspect
import os

import numpy as np


# NOTE This functionality is preliminary and undocumented. It
# *will change*.

# TODO the main thing we need to do here is get rid of the 1-based
# IDs and the "individual" column so that it maps directly into the
# individual table in tskit. We can integrate with PED files by
# putting in a null individual if necessary, but it's horribly
# confusing working with external IDs here.

# Another major issue here is that the sampling strategy is
# bound up with the representation of the pedigree itself. We should
# remove any idea of sampling from this class and instead make the
# samples argument to simulate specify *individual* IDs, ie.
# ts = msprime.sim_ancestry([0, 1], pedigree=pedigree)
# sets the individuals 0 and 1 to be the samples.


class Pedigree:
    """
    Class representing a pedigree for simulations.

    :param ndarray individual: A 1D integer array containing strictly positive unique IDs
        for each individual in the pedigree
    :param ndarray parents: A 2D integer array containing the indices, in the
        individual array, of an individual's parents
    :param ndarray time: A 1D float array containing the time of each
        individual, in generations in the past
    :param int ploidy: The ploidy of individuals in the pedigree. Currently
        only ploidy of 2 is supported
    """

    def __init__(self, individual, parents, time, is_sample=None, sex=None, ploidy=2):
        if ploidy != 2:
            raise NotImplementedError("Ploidy != 2 not currently supported")

        if sex is not None:
            raise NotImplementedError(
                "Assigning individual sexes not currently supported"
            )

        if parents.shape[1] != ploidy:
            raise ValueError(
                "Ploidy {} conflicts with number of parents {}".format(
                    ploidy, parents.shape[1]
                )
            )

        if np.min(individual) <= 0:
            raise ValueError("Individual IDs must be > 0")

        self.individual = individual.astype(np.int32)
        self.num_individuals = len(individual)
        self.parents = parents.astype(np.int32)
        self.time = time.astype(np.float64)
        self.sex = sex
        self.ploidy = int(ploidy)

        self.is_sample = None
        if is_sample is not None:
            self.is_sample = is_sample.astype(np.uint32)
            self.samples = np.array(is_sample)[np.where(is_sample == 1)]
            self.num_samples = len(self.samples)

    def set_samples(self, num_samples=None, sample_IDs=None, probands_only=True):
        if probands_only is not True:
            raise NotImplementedError("Only probands may currently be set as samples.")

        if num_samples is None and sample_IDs is None:
            raise ValueError("Must specify one of num_samples of sample_IDs")

        if num_samples is not None and sample_IDs is not None:
            raise ValueError("Cannot specify both samples and num_samples.")

        self.is_sample = np.zeros((self.num_individuals), dtype=np.uint32)

        all_indices = range(len(self.individual))
        proband_indices = set(all_indices).difference(self.parents.ravel())

        if num_samples is not None:
            self.num_samples = num_samples
            if self.num_samples > len(proband_indices):
                raise ValueError(
                    (
                        "Cannot specify more samples ({}) than there are "
                        "probands in the pedigree ({}) "
                    ).format(self.num_samples, len(proband_indices))
                )
            sample_indices = np.random.choice(
                list(proband_indices), size=self.num_samples, replace=False
            )

        elif sample_IDs is not None:
            self.num_samples = len(sample_IDs)

            indices = all_indices
            if probands_only:
                indices = proband_indices

            sample_set = set(sample_IDs)
            sample_indices = [i for i in indices if self.individual[i] in sample_set]

        if len(sample_indices) != self.num_samples:
            raise ValueError(
                "Sample size mismatch - duplicate sample IDs or sample ID not "
                "in pedigree"
            )

        self.is_sample[sample_indices] = 1

    def get_proband_indices(self):
        all_indices = range(len(self.individual))
        proband_indices = set(all_indices).difference(self.parents.ravel())

        return sorted(proband_indices)

    def get_ll_representation(self):
        """
        Returns the low-level representation of this Pedigree.
        """
        return {
            "individual": self.individual,
            "parents": self.parents,
            "time": self.time,
            "is_sample": self.is_sample,
        }

    # FIXME this shouldn't be a static method, but should be called automatically
    # if a time isn't specified.
    @staticmethod
    def get_times(individual, parent_IDs=None, parents=None, check=False):
        """
        For pedigrees without specified times, crudely assigns times to
        all individuals.
        """
        if parents is None and parent_IDs is None:
            raise ValueError("Must specify either parent IDs or parent indices")

        if parents is None:
            parents = Pedigree.parent_ID_to_index(individual, parent_IDs)

        time = np.zeros(len(individual))
        all_indices = range(len(individual))
        proband_indices = set(all_indices).difference(parents.ravel())
        climber_indices = proband_indices

        t = 0
        while len(climber_indices) > 0:
            next_climbers = []
            for c_idx in climber_indices:
                if time[c_idx] < t:
                    time[c_idx] = t

                next_parents = [p for p in parents[c_idx] if p >= 0]
                next_climbers.extend(next_parents)

            climber_indices = list(set(next_climbers))
            t += 1

        if check:
            Pedigree.check_times(individual, parents, time)

        return time

    @staticmethod
    def check_times(individual, parents, time):
        for i, ind in enumerate(individual):
            for parent_ix in parents[i]:
                if parent_ix >= 0:
                    t1 = time[i]
                    t2 = time[parent_ix]
                    if t1 >= t2:
                        raise ValueError(
                            "Ind {} has time >= than parent {}".format(
                                ind, individual[parent_ix]
                            )
                        )

    @staticmethod
    def parent_ID_to_index(individual, parent_IDs):
        n_inds, n_parents = parent_IDs.shape
        parents = np.zeros(parent_IDs.shape, dtype=int)
        ind_to_index_dict = dict(zip(individual, range(n_inds)))

        if 0 in ind_to_index_dict:
            raise ValueError(
                "Invalid ID: 0 reserved to denote individual" "not in the genealogy"
            )
        ind_to_index_dict[0] = -1

        for i in range(n_inds):
            for j in range(n_parents):
                parent_ID = parent_IDs[i, j]
                parents[i, j] = ind_to_index_dict[parent_ID]

        return parents

    @staticmethod
    def parent_index_to_ID(individual, parents):
        n_inds, n_parents = parents.shape
        parent_IDs = np.zeros(parents.shape, dtype=int)

        for i in range(n_inds):
            for j in range(n_parents):
                parent_ID = 0
                if parents[i, j] >= 0:
                    parent_ID = individual[parents[i, j]]
                    parent_IDs[i, j] = parent_ID

        return parent_IDs

    @staticmethod
    def default_format():
        cols = {
            "individual": 0,
            "parents": [1, 2],
            "time": 3,
            "is_sample": None,
            "sexes": None,
        }

        return cols

    @staticmethod
    def read_txt(pedfile, time_col=None, sex_col=None, **kwargs):
        """
        Creates a Pedigree instance from a text file.
        """
        cols = Pedigree.default_format()
        cols["time"] = time_col
        cols["sexes"] = sex_col

        if sex_col:
            raise NotImplementedError("Specifying sex of individuals not yet supported")

        usecols = []
        for c in cols.values():
            if isinstance(c, collections.abc.Iterable):
                usecols.extend(c)
            elif c is not None:
                usecols.append(c)
        usecols = sorted(usecols)

        data = np.genfromtxt(pedfile, skip_header=1, usecols=usecols, dtype=float)

        individual = data[:, cols["individual"]].astype(int)
        parent_IDs = data[:, cols["parents"]].astype(int)
        parents = Pedigree.parent_ID_to_index(individual, parent_IDs)

        if cols["time"] is not None:
            time = data[:, cols["time"]]
        else:
            time = Pedigree.get_times(individual, parents=parents)

        return Pedigree(individual, parents, time, **kwargs)

    @staticmethod
    def read_npy(pedarray_file, **kwargs):
        """
        Reads pedigree from numpy .npy file with columns:

            ind ID, father array index, mother array index, time

        where time is given in generations.
        """
        basename, ext = os.path.split(pedarray_file)
        pedarray = np.load(pedarray_file)

        cols = Pedigree.default_format()
        if "cols" in kwargs:
            cols = kwargs["cols"]

        individual = pedarray[:, cols["individual"]]
        parents = np.stack([pedarray[:, i] for i in cols["parents"]], axis=1)
        parents = parents.astype(int)
        time = pedarray[:, cols["time"]]

        P = Pedigree(individual, parents, time, **kwargs)

        return P

    def build_array(self):
        cols = Pedigree.default_format()

        col_nums = []
        for v in cols.values():
            if isinstance(v, collections.abc.Iterable):
                col_nums.extend(v)
            elif v is not None:
                col_nums.append(v)

        n_cols = max(col_nums) + 1

        if max(np.diff(col_nums)) > 1:
            raise ValueError(f"Non-sequential columns in pedigree format: {col_nums}")

        pedarray = np.zeros((self.num_individuals, n_cols))
        pedarray[:, cols["individual"]] = self.individual
        pedarray[:, cols["parents"]] = self.parents
        pedarray[:, cols["time"]] = self.time

        return pedarray

    def save_txt(self, fname):
        """
        Saves pedigree in text format with columns:

            ind ID, father array index, mother array index, time

        where time is given in generations.
        """
        pedarray = self.build_array()
        cols = self.default_format()
        cols_to_save = [
            cols["individual"],
            cols["parents"][0],
            cols["parents"][1],
            cols["time"],
        ]
        pedarray = pedarray[cols_to_save]
        parent_IDs = Pedigree.parent_index_to_ID(self.individual, self.parents)
        pedarray[:, cols["parents"]] = parent_IDs

        with open(fname, "w") as f:
            header = "ind\tfather\tmother\ttime\n"
            f.write(header)
            for row in pedarray:
                f.write("\t".join([str(x) for x in row]) + "\n")

    def save_npy(self, fname):
        """
        Saves pedigree in numpy .npy format with columns:

            ind ID, father array index, mother array index, time

        where time is given in generations.
        """
        pedarray = self.build_array()
        np.save(fname, pedarray)

    def asdict(self):
        """
        Returns a dict of arguments to recreate this pedigree
        """
        return {
            key: getattr(self, key)
            for key in inspect.signature(self.__init__).parameters.keys()
            if hasattr(self, key)
        }
