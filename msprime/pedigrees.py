#
# Copyright (C) 2019-2021 University of Oxford
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
Pedigree utilities.
"""
import numpy as np
import tskit

from msprime import demography as demog_mod


class PedigreeBuilder:
    def __init__(self, demography=None):
        if demography is None:
            demography = demog_mod.Demography.isolated_model([1])
        self.demography = demography
        self.tables = tskit.TableCollection(0)
        demography.insert_populations(self.tables)

        assert len(self.tables.individuals) == 0
        self.tables.individuals.metadata_schema = tskit.MetadataSchema(
            {
                "codec": "json",
                "type": "object",
                "properties": {},
                "additionalProperties": True,
            }
        )

    def add_individual(
        self, *, time, is_sample=None, parents=None, population=None, metadata=None
    ):
        if is_sample is None:
            # By default, individuals at time 0 are marked as samples
            is_sample = time == 0

        if parents is None:
            parents = [-1, -1]
        if len(parents) != 2:
            raise ValueError("Must have exactly two parents")

        ind_id = self.tables.individuals.add_row(parents=parents, metadata=metadata)
        if population is None:
            if self.demography.num_populations == 1:
                population = 0
            else:
                raise ValueError(
                    "Must specify the population for each individual in a "
                    "multi-population demography"
                )
        population = self.demography[population].id
        flags = tskit.NODE_IS_SAMPLE if is_sample else 0
        for _ in range(2):
            self.tables.nodes.add_row(
                flags=flags, time=time, individual=ind_id, population=population
            )
        return ind_id

    def finalise(self, sequence_length):
        copy = self.tables.copy()
        copy.sequence_length = sequence_length
        # Not strictly necessary, but it's useful when the tables are written
        # to file and are loaded as as a tree sequence.
        copy.build_index()
        return copy


def sim_pedigree(
    *,
    population_size=None,
    sequence_length=None,
    random_seed=None,
    num_replicates=None,
    end_time=None,
):
    # Internal utility for generating pedigree data. This function is not
    # part of the public API and subject to arbitrary changes/removal
    # in the future.

    # We use the builder to get a standard set of tables, but bypass
    # it so that we can use numpy methods add to it more efficiently.
    tables = PedigreeBuilder().tables
    tables.sequence_length = -1 if sequence_length is None else sequence_length

    population = np.array([tskit.NULL], dtype=np.int32)
    rng = np.random.RandomState(random_seed)

    # To make the semantics compatible with dtwf, the end_time means the
    # *end* of generation end_time
    for time in reversed(range(end_time + 1)):
        N = population_size  # This could be derived from the Demography
        # Make the current generation
        progeny = len(tables.individuals) + np.arange(N, dtype=np.int32)

        # NB this is *with* replacement, so 1 / N chance of selfing
        parents = rng.choice(population, (N, 2))
        tables.individuals.append_columns(
            flags=np.zeros(N, dtype=np.uint32),
            parents=parents.reshape(N * 2),
            parents_offset=np.arange(N + 1, dtype=np.uint32) * 2,
        )
        node_individual = np.zeros(2 * N, dtype=np.int32)
        node_individual[0::2] = progeny
        node_individual[1::2] = progeny
        flags_value = 0 if time > 0 else tskit.NODE_IS_SAMPLE
        tables.nodes.append_columns(
            flags=np.full(2 * N, flags_value, dtype=np.uint32),
            time=np.full(2 * N, time, dtype=np.float64),
            population=np.full(2 * N, 0, dtype=np.int32),
            individual=node_individual,
        )
        population = progeny

    return tables


def parse_fam(fam_file):
    """
    Parse PLINK .fam file and convert to tskit IndividualTable.
    Assumes fam file contains five columns: FID, IID, PAT, MAT, SEX
    :param fam_file: PLINK .fam file object
    :param tskit.TableCollection tc: TableCollection with IndividualTable to
        which the individuals will be added
    """
    individuals = np.loadtxt(
        fname=fam_file,
        dtype=str,
        ndmin=2,  # read file as 2-D table
        usecols=(0, 1, 2, 3, 4),  # only keep FID, IID, PAT, MAT, SEX columns
    )  # requires same number of columns in each row, i.e. not ragged

    id_map = {}  # dict for translating PLINK ID to tskit IndividualTable ID
    for tskit_id, (plink_fid, plink_iid, _pat, _mat, _sex) in enumerate(individuals):
        # include space between strings to ensure uniqueness
        plink_id = f"{plink_fid} {plink_iid}"
        if plink_id in id_map:
            raise ValueError("Duplicate PLINK ID: {plink_id}")
        id_map[plink_id] = tskit_id
    id_map["0"] = -1  # -1 is used in tskit to denote "missing"

    tc = tskit.TableCollection(1)
    tb = tc.individuals
    tb.metadata_schema = tskit.MetadataSchema(
        {
            "codec": "json",
            "type": "object",
            "properties": {
                "plink_fid": {"type": "string"},
                "plink_iid": {"type": "string"},
                "sex": {"type": "integer"},
            },
            "required": ["plink_fid", "plink_iid", "sex"],
            "additionalProperties": True,
        }
    )
    for plink_fid, plink_iid, pat, mat, sex in individuals:
        sex = int(sex)
        if not (sex in range(3)):
            raise ValueError(
                "Sex must be one of the following: 0 (unknown), 1 (male), 2 (female)"
            )
        metadata_dict = {"plink_fid": plink_fid, "plink_iid": plink_iid, "sex": sex}
        pat_id = f"{plink_fid} {pat}" if pat != "0" else pat
        mat_id = f"{plink_fid} {mat}" if mat != "0" else mat
        tb.add_row(
            parents=[
                id_map[pat_id],
                id_map[mat_id],
            ],
            metadata=metadata_dict,
        )
    tc.sort()

    return tb
