#
# Copyright (C) 2019-2022 University of Oxford
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

import numpy as np
import tskit

from msprime import demography as demog_mod


@dataclasses.dataclass
class Column:
    name: str

    def parse_value(self, value: str) -> str:
        return value


class IdColumn(Column):
    pass


class NumericColumn(Column):
    def parse_value(self, value: str) -> float:
        return float(value)


class BooleanColumn(Column):
    def parse_value(self, value: str) -> bool:
        return bool(int(value))


class StringOrIntColumn(Column):
    def parse_value(self, value: str) -> str | int:
        try:
            return int(value)
        except ValueError:
            pass
        return value


@dataclasses.dataclass
class TextFileIndividual:
    id: str  # noqa: A003
    parent0: str | None
    parent1: str | None
    time: float
    is_sample: bool | None = None
    population: int = 0


def _parse_pedigree_header(
    header: str, demography: None | demog_mod.Demography
) -> list[Column]:
    tokens = header.split()
    if tokens[0] != "#":
        raise ValueError("First line must be the header and start with #")
    required_cols = [
        IdColumn("id"),
        IdColumn("parent0"),
        IdColumn("parent1"),
        NumericColumn("time"),
    ]
    demography = (
        demog_mod.Demography.isolated_model([1]) if demography is None else demography
    )
    optional_cols = [
        BooleanColumn("is_sample"),
        # There is the potential for ambiguity here in the population column
        # if the string name of the population was an integer not equal to
        # its ID. However, the Demography class requires that population names
        # are valid Python identifiers, so this can't happen in practice.
        StringOrIntColumn("population"),
    ]
    col_names = [col.name for col in required_cols]
    tokens = tokens[1:]
    if tokens[: len(required_cols)] != col_names:
        raise ValueError(f"The {col_names} columns are required")

    columns = list(required_cols)
    optional_col_map = {col.name: col for col in optional_cols}
    for colname in tokens[len(required_cols) :]:
        if colname not in optional_col_map:
            raise ValueError(f"Column '{colname}' not supported")
        columns.append(optional_col_map[colname])
    return columns


def parse_pedigree(text_file, demography=None, sequence_length=None):
    """
    Parse a text describing a pedigree used for input to the
    :class:`.FixedPedigree` ancestry model.
    """
    # First line must be the header and contain the column headers
    header = next(text_file, None)
    if header is None:
        raise ValueError("Pedigree file must contain at least a header")
    columns = _parse_pedigree_header(header, demography=demography)
    individuals = []
    tsk_id_map = {".": -1}
    for tsk_id, row in enumerate(text_file):
        line = tsk_id + 2
        tokens = row.split()
        if len(tokens) != len(columns):
            raise ValueError(
                "Incorrect number of columns "
                f"(found={len(tokens)} need={len(columns)}) at line {line}: {row}"
            )
        kwargs = {}
        for value, column in zip(tokens, columns):
            try:
                kwargs[column.name] = column.parse_value(value)
            except ValueError as ve:
                raise ValueError(
                    f"Error parsing value for column '{column.name}' "
                    f"on line {line}: {ve}"
                )
        ind = TextFileIndividual(**kwargs)
        if ind.id == ".":
            raise ValueError(f"'.' cannot be used as an individual ID (line {line})")
        if ind.id in tsk_id_map:
            raise ValueError(f"Duplicate ID at line {line}")
        tsk_id_map[ind.id] = len(individuals)
        individuals.append(ind)

    # Convert the list of individuals to a tskit pedigree.
    builder = PedigreeBuilder(
        individuals_metadata_schema=tskit.MetadataSchema.permissive_json(),
        demography=demography,
    )
    for ind in individuals:
        parents = []
        for parent in [ind.parent0, ind.parent1]:
            if parent not in tsk_id_map:
                raise ValueError(
                    f"Parent ID '{parent}' not defined (line {tsk_id_map[ind.id] + 2})"
                )
            parents.append(tsk_id_map[parent])
        builder.add_individual(
            parents=parents,
            time=ind.time,
            is_sample=ind.is_sample,
            metadata={"file_id": ind.id},
            population=ind.population,
        )
    return builder.finalise(sequence_length)


def write_pedigree(ts, out):

    print("# id\tparent0\tparent1\ttime\tis_sample\tpopulation", file=out)
    for ind in ts.individuals():
        if len(ind.nodes) != 2:
            raise ValueError(
                "Invalid pedigree format: each individual must be associated with "
                "two nodes"
            )
        nodes = [ts.node(node) for node in ind.nodes]
        time = {node.time for node in nodes}
        is_sample = {node.is_sample() for node in nodes}
        population = {node.population for node in nodes}
        if len(time) != 1:
            raise ValueError("Pedigree individuals must have the same node time")
        if len(is_sample) != 1:
            raise ValueError("Pedigree individuals must have the same sample status")
        if len(population) != 1:
            raise ValueError("Pedigree individuals must have the same node population")

        time = time.pop()
        is_sample = is_sample.pop()
        population = population.pop()
        parents = []
        for parent in ind.parents:
            if parent == tskit.NULL:
                parents.append(".")
            else:
                parents.append(parent)
        print(
            f"{ind.id}\t{parents[0]}\t{parents[1]}\t{time}\t{is_sample}\t"
            f"{population}",
            file=out,
        )


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

    def add_individual(
        self,
        *,
        time,
        is_sample=None,
        parents=None,
        population=None,
        metadata=None,
    ):
        """
        Adds the specified individual, returning its ID.
        """
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

    def finalise(self, sequence_length=None):
        """
        Returns the table collection representing the defined pedigree.
        """
        copy = self.tables.copy()
        copy.sequence_length = 1
        # Not strictly necessary, but it's useful when the tables are written
        # to file and are loaded as as a tree sequence.
        copy.build_index()
        copy.sequence_length = -1 if sequence_length is None else sequence_length
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
