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
from __future__ import annotations

import dataclasses
from typing import Any
from typing import Dict
from typing import List
from typing import Union

import numpy as np
import tskit

from msprime import demography as demog_mod


@dataclasses.dataclass
class PedColumn:
    name: str
    arg_name: str = ""

    def parse_value(self, value: str) -> str:
        return value


class IdColumn(PedColumn):
    pass


class NumericColumn(PedColumn):
    def parse_value(self, value: str) -> float:
        return float(value)


class BooleanColumn(PedColumn):
    def parse_value(self, value: str) -> bool:
        return bool(int(value))


def _parse_pedigree_header(header: str) -> list[PedColumn]:
    tokens = header.split()
    if tokens[0] != "#":
        raise ValueError("First line must be the header and start with #")
    required_cols = {
        "IID": IdColumn("IID"),
        "FID": IdColumn("FID"),
        "MID": IdColumn("MID"),
        "TIME": NumericColumn("TIME", "time"),
    }
    optional_cols = {
        "IS_SAMPLE": BooleanColumn("IS_SAMPLE", "is_sample"),
    }
    if len(tokens) < len(required_cols) + 1:
        raise ValueError(f"The {list(required_cols)} columns are required")
    tokens = tokens[1:]
    # Dicts are sorted since Python 3.7
    if tokens[: len(required_cols)] != list(required_cols):
        raise ValueError(f"The {list(required_cols)} columns are required")

    columns = [required_cols[colname] for colname in required_cols]
    for colname in tokens[len(required_cols) :]:
        if colname not in optional_cols:
            raise ValueError(f"Column {colname} not supported")
        columns.append(optional_cols[colname])
    return columns


@dataclasses.dataclass
class Individual:
    external_id: str | None = None
    parents: list[str | int] = dataclasses.field(default_factory=list)
    metadata: dict = dataclasses.field(default_factory=dict)


class PedigreeBuilder:
    def __init__(self, demography=None, individuals_metadata_schema=None):
        if demography is None:
            demography = demog_mod.Demography.isolated_model([1])
        self.demography = demography
        self.tables = tskit.TableCollection(0)
        self.tables.time_units = "generations"
        demography.insert_populations(self.tables)
        # Work around https://github.com/tskit-dev/tskit/issues/2080
        self.individuals = self.tables.individuals
        self.nodes = self.tables.nodes
        assert len(self.individuals) == 0
        if individuals_metadata_schema is not None:
            self.individuals.metadata_schema = individuals_metadata_schema

    def set_parents(self, ind_id, parents):
        row = self.tables.individuals[ind_id]
        self.tables.individuals[ind_id] = row.replace(parents=parents)

    def add_individual(
        self,
        *,
        time,
        external_id=None,
        is_sample=None,
        parents=None,
        population=None,
        metadata=None,
    ):
        if is_sample is None:
            # By default, individuals at time 0 are marked as samples
            is_sample = time == 0

        if parents is None:
            parents = [-1, -1]
        if len(parents) != 2:
            raise ValueError("Must have exactly two parents")

        if external_id is not None:
            metadata = {"external_id": external_id}

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

    def add_individuals(self, *, parents, time):
        population = 0
        flags = 0 if time > 0 else tskit.NODE_IS_SAMPLE

        N = parents.shape[0]
        ind_ids = np.arange(N, dtype=np.int32) + len(self.individuals)
        self.individuals.append_columns(
            flags=np.zeros(N, dtype=np.uint32),
            parents=parents.reshape(2 * N),
            parents_offset=np.arange(N + 1, dtype=np.uint32) * 2,
        )
        node_individual = np.zeros(2 * N, dtype=np.int32)
        node_individual[0::2] = ind_ids
        node_individual[1::2] = ind_ids
        self.nodes.append_columns(
            flags=np.full(2 * N, flags, dtype=np.uint32),
            time=np.full(2 * N, time, dtype=np.float64),
            population=np.full(2 * N, population, dtype=np.int32),
            individual=node_individual,
        )
        return ind_ids

    def _parse_pedigree_row(self, row: str, line_num: int, columns: list[PedColumn]):
        tokens = row.split()
        if len(tokens) != len(columns):
            raise ValueError(f"Incorrect number of columns at line {line_num}: {row}")
        external_id = tokens[0]
        parents = tokens[1:3]
        kwargs = {}
        for value, column in zip(tokens[3:], columns[3:]):
            # print(value, column)
            kwargs[column.arg_name] = column.parse_value(value)
            # print(column, value)
            # column.validate(a

        print(external_id, parents, kwargs)

    def parse_text(self, ped_file):
        """
        Parse a text describing a pedigree used for input to the
        :class:`.FixedPedigree` ancestry model.

        # iid p0 p1 time time
        """
        # First line must be the header and contain the column headers
        header = next(ped_file)
        columns = _parse_pedigree_header(header)
        print(columns)
        for line_num, row in enumerate(ped_file, 1):
            self._parse_pedigree_row(row, line_num, columns)

    def finalise(self, sequence_length):
        copy = self.tables.copy()
        copy.sequence_length = sequence_length
        # Not strictly necessary, but it's useful when the tables are written
        # to file and are loaded as as a tree sequence.
        copy.build_index()
        return copy


def sim_pedigree_backward(
    builder,
    rng,
    *,
    population_size,
    num_samples,
    end_time,
):
    ancestors = np.arange(num_samples, dtype=np.int32)
    parent_offset = num_samples

    # Because we don't have overlapping generations we can add the ancestors
    # to the pedigree once the previous generation has been produced
    for time in range(0, end_time):
        population = np.arange(population_size, dtype=np.int32)
        parents = rng.choice(population, (len(ancestors), 2))
        unique_parents = np.unique(parents)
        parent_ids = np.searchsorted(unique_parents, parents).astype(np.int32)
        assert np.all(unique_parents[parent_ids] == parents)
        a = builder.add_individuals(parents=parent_ids + parent_offset, time=time)
        assert np.array_equal(a, ancestors)
        ancestors = np.arange(len(unique_parents), dtype=np.int32) + parent_offset
        parent_offset += len(ancestors)
    # Add the founders
    builder.add_individuals(
        parents=np.full((len(ancestors), 2), tskit.NULL, dtype=np.int32), time=end_time
    )
    return builder.finalise(1)


def sim_pedigree_forward(
    builder,
    rng,
    *,
    population_size,
    end_time,
):
    population = np.array([tskit.NULL], dtype=np.int32)

    # To make the semantics compatible with dtwf, the end_time means the
    # *end* of generation end_time
    for time in reversed(range(end_time + 1)):
        N = population_size  # This could be derived from the Demography
        # NB this is *with* replacement, so 1 / N chance of selfing
        parents = rng.choice(population, (N, 2))
        population = builder.add_individuals(parents=parents, time=time)
    return builder.finalise(1)


def sim_pedigree(
    *,
    population_size=None,
    num_samples=None,
    sequence_length=None,
    random_seed=None,
    end_time=None,
    direction="forward",
):
    # Internal utility for generating pedigree data. This function is not
    # part of the public API and subject to arbitrary changes/removal
    # in the future.
    num_samples = population_size if num_samples is None else num_samples
    builder = PedigreeBuilder()
    rng = np.random.RandomState(random_seed)

    if direction == "forward":
        if num_samples != population_size:
            raise ValueError(
                "num_samples must be equal to population_size for forward simulation"
            )
        tables = sim_pedigree_forward(
            builder,
            rng,
            population_size=population_size,
            end_time=end_time,
        )
    elif direction == "backward":
        tables = sim_pedigree_backward(
            builder,
            rng,
            population_size=population_size,
            num_samples=num_samples,
            end_time=end_time,
        )
    else:
        raise ValueError("unknown time direction; choose 'backward' or 'forward'")

    tables.sequence_length = -1 if sequence_length is None else sequence_length
    return tables
