#
# Copyright (C) 2015-2021 University of Oxford
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
Module responsible for defining and debugging demographic models.
"""
from __future__ import annotations

import collections
import copy
import dataclasses
import enum
import inspect
import itertools
import logging
import math
import sys
import textwrap
import warnings
from typing import Any
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import numpy as np
import tskit

from . import ancestry
from . import core
from . import species_trees


logger = logging.getLogger(__name__)


def check_num_populations(num_populations):
    """
    Check if an input number of populations is valid.
    """
    if num_populations < 1:
        raise ValueError("Must have at least one population")


def check_migration_rate(migration_rate):
    """
    Check if an input migration rate makes sense.
    """
    if migration_rate < 0:
        raise ValueError("Migration rates must be non-negative")


@dataclasses.dataclass
class Population:
    """
    Define a :ref:`population <sec_demography_populations>` in a
    :class:`.Demography`.

    :ivar initial_size: The absolute size of the population at time zero.
    :vartype initial_size: float
    :var growth_rate: The exponential growth rate of the
        population per generation (forwards in time).
        Growth rates can be negative. This is zero for a
        constant population size, and positive for a population that has been
        growing. Defaults to 0.
    :vartype growth_rate: float
    :ivar name: The name of the population. If specified this must be a uniquely
        identifying string and must be a valid Python identifier (i.e., could be
        used as a variable name in Python code).
    :vartype name: str
    :ivar description: A short description of the population. Defaults to the
        empty string if not specified.
    :vartype description: str
    :ivar extra_metadata: A JSON-encodable dictionary of metadata items to be
        stored in the associated tskit population object. This dictionary
        must not contain keys for any of the pre-defined metadata items.
    :vartype extra_metadata: dict
    :ivar sampling_time: The default time at which samples are drawn from
        this population. See the :ref:`sec_ancestry_samples_sampling_time`
        section for more details.
    :vartype sampling_time: float
    :ivar id: The integer ID of this population within the parent
        :class:`.Demography`. This attribute is assigned by the Demography
        class and should not be set or changed by user code.
    :vartype id: int
    """

    initial_size: float = 0.0
    growth_rate: float = 0.0
    name: Union[str, None] = None
    description: str = ""
    extra_metadata: dict = dataclasses.field(default_factory=dict)
    sampling_time: Union[float, None] = None
    initially_active: Union[bool, None] = None

    # Keeping this as something we can init because this stops us
    # doing things like round-tripping through repr. The warning
    # above should suffice.
    id: Union[int, None] = dataclasses.field(default=None)  # noqa: A003

    def asdict(self):
        return dataclasses.asdict(self)

    def validate(self):
        # TODO more checks
        if self.initial_size < 0:
            raise ValueError("Negative population size")
        if self.name is None:
            raise ValueError("A population name must be set.")
        if not self.name.isidentifier():
            raise ValueError("A population name must be a valid Python identifier")


@dataclasses.dataclass
class Demography:
    """
    The definition of a demographic model for an msprime simulation,
    consisting of a set of populations, a migration matrix, and a list
    of demographic events. See the :ref:`sec_demography_definitions`
    section for precise mathematical definitions of these concepts.
    """

    populations: List[Population] = dataclasses.field(default_factory=list)
    events: List = dataclasses.field(default_factory=list)
    # Until we can use numpy type hints properly, it's not worth adding them
    # here. We still have to add in ignores below for indexed assignment errors.
    migration_matrix: Union[Any, None] = None

    def __post_init__(self):
        if self.migration_matrix is None:
            N = self.num_populations
            self.migration_matrix = np.zeros((N, N))

        # People might get cryptic errors from passing in copies of the same
        # population, so check for it.
        if len({id(pop) for pop in self.populations}) != len(self.populations):
            raise ValueError("Population objects must be distinct")

        # Assign the IDs and default names, if needed.
        for j, population in enumerate(self.populations):
            if population.id is not None:
                raise ValueError(
                    "Population ID should not be set before using to create "
                    "a Demography"
                )
            population.id = j
            if population.name is None:
                population.name = f"pop_{j}"
        self._validate_populations()

    def add_population(
        self,
        *,
        initial_size=None,
        growth_rate=None,
        name=None,
        description=None,
        extra_metadata=None,
        sampling_time=None,
        initially_active=None,
    ) -> Population:
        """
        TODO document
        """
        N = self.num_populations
        population = Population()
        population.id = N
        population.growth_rate = 0 if growth_rate is None else growth_rate
        population.initial_size = 1 if initial_size is None else initial_size
        population.name = f"pop_{population.id}" if name is None else name
        population.description = "" if description is None else description
        population.extra_metadata = {} if extra_metadata is None else extra_metadata
        population.sampling_time = sampling_time
        population.initially_active = initially_active
        self.populations.append(population)
        M = self.migration_matrix
        self.migration_matrix = np.zeros((N + 1, N + 1))
        self.migration_matrix[:N, :N] = M
        self._validate_populations()
        return population

    def _add_population_from_old_style(
        self, pop_config: PopulationConfiguration, name: Union[str, None] = None
    ) -> Population:
        population = self.add_population(
            name=name,
            initial_size=pop_config.initial_size,
            growth_rate=pop_config.growth_rate,
        )
        metadata = pop_config.metadata
        if metadata is not None and isinstance(metadata, collections.abc.Mapping):
            metadata = metadata.copy()
            if "name" in metadata:
                population.name = metadata.pop("name")
                if name is not None and name != population.name:
                    # Maybe this should be a warning, or just ignored entirely?
                    raise ValueError(
                        "Population name already set in old-style metadata and doesn't "
                        "match supplied name"
                    )
            if "description" in metadata:
                population.description = metadata.pop("description")
        population.extra_metadata = metadata
        return population

    def add_event(self, event: DemographicEvent) -> DemographicEvent:
        if not isinstance(event, DemographicEvent):
            raise TypeError("Events must be instances of DemographicEvent")
        event.demography = self
        self.events.append(event)
        return event

    def set_migration_rate(
        self, source: Union[str, int], dest: Union[str, int], rate: float
    ):
        """
        Sets the backwards-time rate of migration from the specified ``source``
        population to ``dest`` to the specified value. This has the effect of
        setting ``demography.migration_matrix[source, dest] = rate``.

        .. warning:: Note this is the **backwards time** migration rate and that
            ``source`` and ``dest`` are from the perspective of lineages in the
            coalescent process. See :ref:`sec_demography_migration` for more
            details and clarification on this vital point.

        The ``source`` and ``dest`` populations can be referred to either by
        their integer ``id`` or string ``name`` values.

        :param str,int source: The source population from which lineages originate
            in the backwards-time process.
        :param str,int dest: The destination population where lineages are move
            to in the backwards-time process.
        :param float rate: The per-generation migration rate.
        """
        source = self[source].id
        dest = self[dest].id
        if source == dest:
            raise ValueError("The source and dest populations must be different")
        self.migration_matrix[source, dest] = rate  # type: ignore

    def set_symmetric_migration_rate(
        self,
        populations: List[Union[str, int]],
        rate: float,
    ):
        """
        Sets the symmetric migration rate between all pairs of populations in
        the specified list to the specified value. For a given pair of population
        IDs ``j`` and ``k``, this sets ``demography.migration_matrix[j, k] = rate``
        and ``demography.migration_matrix[k, j] = rate``.

        Populations may be specified either by their integer IDs or by
        their string names.

        :param list populations: An iterable of population identifiers (integer
            IDs or string names).
        :param float rate: The value to set the migration matrix entries to.
        """
        # There's an argument for not checking this so that corner cases on
        # single population models can be handled. However, it's nearly always
        # going to be a user error where someone forgets the second population
        # so it seems better to raise an error to prevent hard-to-detect mistakes.
        if len(populations) < 2:
            raise ValueError("Must specify at least two populations")
        pop_ids = [self[identifier].id for identifier in populations]
        for pop_j, pop_k in itertools.combinations(pop_ids, 2):
            self.migration_matrix[pop_j, pop_k] = rate  # type: ignore
            self.migration_matrix[pop_k, pop_j] = rate  # type: ignore

    # Demographic events.

    def add_population_split(
        self, time: float, *, derived: List[Union[str, int]], ancestral: Union[str, int]
    ) -> PopulationSplit:
        """
        Adds a population split event at the specified time. In a population
        split event all lineages from the (more recent) derived populations
        move to the (more ancient) ancestral population. Forwards in time,
        this corresponds to the ancestral population splitting into the
        derived populations.

        .. todo:: Add some links here, to the documentation and to the
            other types of events we might want to support.

        In addition to moving lineages from the derived population(s) into the
        ancestral population, a population split has the following additional
        effects:

        - All migration rates to and from the derived populations are set to 0.
        - Population sizes and growth rates for the derived populations are set
          to 0.
        - The ``sampling_time`` of the ``ancestral`` :class:`.Population` is set
          to the time of this event, **if** the ``sampling_time`` for the
          ancestral population has not already been set.

        :param float time: The time at which this event occurs in generations.
        :param list(str, int) derived: The derived populations.
        :param str, int ancestral: The ancestral population.
        """
        pop = self[ancestral]
        if pop.initially_active is None:
            pop.initially_active = False
            if pop.sampling_time is None:
                pop.sampling_time = time
        return self.add_event(
            PopulationSplit(time=time, derived=derived, ancestral=ancestral)
        )

    def add_admixture(
        self,
        time: float,
        *,
        derived: Union[str, int],
        ancestral: List[Union[str, int]],
        proportions: List[float],
    ) -> Admixture:
        """
        Adds an admixture event at the specified time. In an admixture
        event all lineages from a (more recent) ``derived`` population
        move to a list of (more ancient) ``ancestral`` populations according
        to a list of ``proportions``, such that a given lineage has a
        probability ``proportions[j]`` of being moved to the population
        ``ancestral[j]``. This movement of lineages backwards in time
        corresponds to the initial state of the admixed derived population
        the specified ``time`` being composed of individuals from the
        specified ``ancestral`` populations in the specified ``proportions``.

        .. todo:: Add some links here, to the documentation and to the
            other types of events we might want to support.

        In addition to moving lineages from the derived population into the
        ancestral population(s), an admixture has the following additional
        effects:

        - All migration rates to and from the derived population are set to 0.
        - Population sizes and growth rates for the derived population are set
          to 0, and the poulation is marked as inactive.

        :param float time: The time at which this event occurs in generations.
        :param str, int derived: The derived population.
        :param list(str, int) ancestral: The ancestral populations.
        :param list(float) proportions: The proportion of the derived population
            from each of the ancestral populations at the time of the event.
        """
        # Useful feature here might be to support taking n - 1 proportion values
        # and computing 1 - sum for the last value. Could be tedious for users to
        # do this manually.
        if not math.isclose(sum(proportions), 1.0):
            raise ValueError("Sum of the admixture proportions must be approximately 1")
        return self.add_event(
            Admixture(
                time=time, derived=derived, ancestral=ancestral, proportions=proportions
            )
        )

    def add_mass_migration(self, time, *, source, dest, proportion):
        # TODO document. Not clear whether we document this with a warning
        # or we just leave it as population split, etc.
        return self.add_event(MassMigration(time, source, dest, proportion))

    def add_migration_rate_change(
        self,
        time: float,
        *,
        rate: float,
        source: Union[int, str, None] = None,
        dest: Union[int, str, None] = None,
    ) -> MigrationRateChange:
        """
        Changes the rate of migration from one deme to another to a new value at a
        specific time. Migration rates are specified in terms of the rate at which
        lineages move from population ``source`` to ``dest`` during the progress of
        the simulation.

        .. warning::
            Note that ``source`` and ``dest`` are from the perspective of the
            coalescent process, i.e. **backwards in time**; please see the
            :ref:`sec_demography_migration` section for more details on the
            interpretation of this migration model.

        By default, ``source=None`` and ``dest=None``, which results in all
        non-diagonal elements of the migration matrix being changed to the new
        rate. If ``source`` and ``dest`` are specified, they must refer to valid
        populations (either integer IDs or string names).

        :param float time: The time at which this event occurs in generations.
        :param float rate: The new per-generation migration rate.
        :param str, int source: The ID of the source population.
        :param str, int dest: The ID of the destination population.
        :param int source: The source population ID.
        """
        source = -1 if source is None else source
        dest = -1 if dest is None else dest
        return self.add_event(
            MigrationRateChange(time=time, source=source, dest=dest, rate=rate)
        )

    def add_symmetric_migration_rate_change(
        self, time: float, populations: List[Union[str, int]], rate: float
    ) -> SymmetricMigrationRateChange:
        """
        Sets the symmetric migration rate between all pairs of populations in
        the specified list to the specified value. For a given pair of population
        IDs ``j`` and ``k``, this sets ``migration_matrix[j, k] = rate``
        and ``migration_matrix[k, j] = rate``.

        Populations may be specified either by their integer IDs or by
        their string names.

        :param float time: The time at which this event occurs in generations.
        :param list populations: An sequence of population identifiers (integer
            IDs or string names).
        :param float rate: The new migration rate.
        """
        return self.add_event(
            SymmetricMigrationRateChange(time=time, populations=populations, rate=rate)
        )

    def add_population_parameters_change(
        self,
        time: float,
        *,
        initial_size: Union[float, None] = None,
        growth_rate: Union[float, None] = None,
        population: Union[int, None] = None,
    ) -> PopulationParametersChange:
        """
        Changes the size parameters of a population (or all populations)
        at a given time.

        :param float time: The length of time ago at which this event
            occurred.
        :param float initial_size: The absolute size of the population
            at the beginning of the time slice starting at ``time``. If None,
            the initial_size of the population is computed according to
            the initial population size and growth rate over the preceding
            time slice.
        :param float growth_rate: The new per-generation growth rate. If None,
            the growth rate is not changed. Defaults to None.
        :param str, int population: The ID of the population affected. If
            ``population`` is None, the changes affect all populations
            simultaneously.
        """
        event = PopulationParametersChange(
            time,
            initial_size=initial_size,
            growth_rate=growth_rate,
            population=population,
        )
        return self.add_event(event)

    def add_simple_bottleneck(
        self,
        time: float,
        population: Union[int, str],
        proportion: Union[float, None] = None,
    ) -> SimpleBottleneck:
        proportion = 1.0 if proportion is None else proportion
        return self.add_event(
            SimpleBottleneck(time=time, population=population, proportion=proportion)
        )

    def add_instantaneous_bottleneck(
        self, time: float, *, population: Union[str, int], strength: float
    ) -> InstantaneousBottleneck:
        return self.add_event(
            InstantaneousBottleneck(time=time, population=population, strength=strength)
        )

    def add_census(self, time: float) -> CensusEvent:
        """
        Add a "census" event at the specified time. In a census we add a node
        to each branch of every tree, thus recording the population that each
        lineage is in at the specified time.

        This may be used to record all ancestral haplotypes present at that
        time, and to extract other information related to these haplotypes: for
        instance to trace the local ancestry of a sample back to a set of
        contemporaneous ancestors, or to assess whether a subset of samples has
        coalesced more recently than the census time.

        See :ref:`sec_ancestry_census_events` for more details.
        """
        return self.add_event(CensusEvent(time))

    def _populations_table(self):
        cols = [
            ("id", ""),
            ("name", ""),
            ("description", ""),
            ("initial_size", ".1f"),
            ("growth_rate", ".2f"),
            ("sampling_time", ".2g"),
            ("extra_metadata", ""),
        ]
        data = [
            [f"{getattr(pop, attr):{fmt}}" for attr, fmt in cols]
            for pop in self.populations
        ]
        return [title for title, _ in cols], data

    def _populations_text(self):
        col_titles, data = self._populations_table()
        alignments = ["^", "<", "<", "<", "^", ">", "<"]
        data = [
            [
                [item.as_text() if isinstance(item, core.TableEntry) else item]
                for item in row
            ]
            for row in data
        ]
        return core.text_table(
            "Populations", [[title] for title in col_titles], alignments, data
        )

    def _populations_html(self):
        col_titles, data = self._populations_table()
        return core.html_table("Populations", col_titles, data)

    def _migration_rate_info(self, source, dest, rate):
        extra = None
        if source != dest:
            source_name = self.populations[source].name
            dest_name = self.populations[dest].name
            extra = (
                "Backwards in time migration rate from population "
                f"{source_name} to {dest_name} = {rate} per generation. "
                "Forwards in time, this is the expected number of migrants "
                f"moving from {dest_name} to {source_name} "
                f"per generation, divided by the size of {source_name}."
            )
        return core.TableEntry(f"{rate:.4g}", extra)

    def _migration_matrix_table(self):
        col_titles = [""] + [pop.name for pop in self.populations]
        data = []
        for j in range(self.num_populations):
            row = [self.populations[j].name] + [
                self._migration_rate_info(j, k, self.migration_matrix[j, k])
                for k in range(self.num_populations)
            ]
            data.append(row)
        return col_titles, data

    def _migration_matrix_text(self):
        col_titles, data = self._migration_matrix_table()
        alignments = ">" + "^" * self.num_populations
        data = [
            [
                [item.as_text() if isinstance(item, core.TableEntry) else item]
                for item in row
            ]
            for row in data
        ]
        return core.text_table(
            "Migration Matrix", [[title] for title in col_titles], alignments, data
        )

    def _migration_matrix_html(self):
        col_titles, data = self._migration_matrix_table()
        return core.html_table("Migration matrix", col_titles, data)

    def _events_text(self, events, title="Events"):
        col_titles = [["time"], ["type"], ["parameters"], ["effect"]]
        alignments = "><<<"
        data = []
        for event in events:
            type_text = textwrap.wrap(event._type_str, 15)
            description = textwrap.wrap(event._parameters(), 22)
            effect = textwrap.wrap(event._effect(), 38)
            row = [[f"{event.time:.4g}"], type_text, description, effect]
            data.append(row)
        return core.text_table(
            title, col_titles, alignments, data, internal_hlines=True
        )

    def _events_html(self, events, title="Events"):
        col_titles = ["time", "type", "parameters", "effect"]
        data = []
        for event in events:
            class_name = event.__class__.__name__
            # TODO change this to stable when 1.0 is released.
            type_html = (
                "<a href='https://tskit.dev/msprime/docs/latest/api.html#msprime."
                f"{class_name}'>{event._type_str}</a>"
            )
            row = [f"{event.time:.4g}", type_html, event._parameters(), event._effect()]
            data.append(row)
        return core.html_table(title, col_titles, data)

    def _repr_html_(self):
        resolved = self.validate()
        return (
            '<div style="margin-left:20px">'
            + resolved._populations_html()
            + resolved._migration_matrix_html()
            + resolved._events_html(self.events)
            + "</div>"
        )

    def __str__(self):
        resolved = self.validate()
        populations = resolved._populations_text()
        migration_matrix = resolved._migration_matrix_text()
        events = resolved._events_text(self.events)

        def indent(table):
            lines = table.splitlines()
            s = "╟  " + lines[0] + "\n"
            for line in lines[1:]:
                s += "║  " + line + "\n"
            return s

        s = (
            "Demography\n"
            + indent(populations)
            + indent(migration_matrix)
            + indent(events)
        )
        return s

    def __getitem__(self, identifier):
        """
        Returns the population with the specified ID or name.
        """
        if isinstance(identifier, str):
            for population in self.populations:
                if population.name == identifier:
                    return population
            else:
                raise KeyError(f"Population with name '{identifier}' not found")
        elif core.isinteger(identifier):
            # We don't support negative indexing here because -1 is used as
            # way to refer to *all* populations in demographic events, and
            # it would be too easy to introduce bugs in old code if we changed
            # the meaning of this.
            if identifier < 0 or identifier >= self.num_populations:
                raise KeyError(f"Population id {identifier} out of bounds")
            return self.populations[identifier]
        raise TypeError("Keys must be either string population names or integer IDs")

    def __contains__(self, identifier):
        """
        Support "in" lookups, i.e. ``if "YRI" in demography``.
        """
        try:
            self.__getitem__(identifier)
        except KeyError:
            return False
        return True

    @property
    def num_populations(self):
        return len(self.populations)

    @property
    def num_events(self):
        return len(self.events)

    def _validate_populations(self):
        names = set()
        for j, population in enumerate(self.populations):
            if population.id != j:
                raise ValueError(
                    "Incorrect population ID. ID values should not be updated "
                    "by users. Please use Demography.add_population to add extra "
                    "populations after initialisation."
                )
            population.validate()
            if population.name in names:
                raise ValueError(f"Duplicate population name: '{population.name}'")
            names.add(population.name)

    def validate(self):
        """
        Checks the demography looks sensible and raises errors/warnings
        appropriately, and return a copy in which all default values have
        been appropriately resolved.
        """
        self._validate_populations()
        migration_matrix = np.array(self.migration_matrix)
        N = self.num_populations
        if migration_matrix.shape != (N, N):
            raise ValueError(
                "migration matrix must be a N x N square matrix encoded "
                "as a list-of-lists or numpy array, where N is the number "
                "of populations. The diagonal "
                "elements of this matrix must be zero. For example, a "
                "valid matrix for a 3 population system is "
                "[[0, 1, 1], [1, 0, 1], [1, 1, 0]]"
            )
        last_event = None
        for event in self.events:
            if not isinstance(event, DemographicEvent):
                raise TypeError(
                    "Demographic events must be a list of DemographicEvent "
                    "instances sorted in non-decreasing order of time."
                )
            if last_event is not None:
                if last_event.time > event.time:
                    raise ValueError(
                        "Events must be time-sorted. Please use demography.sort_events()"
                        "if you add events out of order."
                    )
            last_event = event
        resolved = copy.deepcopy(self)
        for population in resolved.populations:
            if population.sampling_time is None:
                population.sampling_time = 0
            if population.initially_active is None:
                population.initially_active = True
        return resolved

    def sort_events(self):
        # Sort demographic events by time. Sorting is stable so the relative
        # order of events at the same time will be preserved.
        self.events.sort(key=lambda de: de.time)

    def insert_populations(self, tables):
        """
        Insert population definitions for this demography into the specified
        set of tables.

        :meta private:
        """
        metadata_schema = tskit.MetadataSchema(
            {
                "codec": "json",
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "description": {"type": ["string", "null"]},
                },
                # The name and description fields are always filled out by
                # msprime, so we tell downstream tools this by making them
                # "required" by the schema.
                "required": ["name", "description"],
                "additionalProperties": True,
            }
        )
        assert len(tables.populations) == 0
        tables.populations.metadata_schema = metadata_schema
        for population in self.populations:
            metadata = {
                "name": population.name,
                "description": population.description,
            }
            if population.extra_metadata is not None:
                intersection = set(population.extra_metadata.keys()) & set(
                    metadata.keys()
                )
                if len(intersection) > 0:
                    printed_list = list(sorted(intersection))
                    raise ValueError(
                        f"Cannot set standard metadata key(s) {printed_list} "
                        "using extra_metadata. Please set using the corresponding "
                        "property of the Population class."
                    )
                metadata.update(population.extra_metadata)
            tables.populations.add_row(metadata=metadata)

    def asdict(self):
        return {
            "populations": [pop.asdict() for pop in self.populations],
            "events": [event.asdict() for event in self.events],
            "migration_matrix": self.migration_matrix.tolist(),
        }

    def debug(self):
        """
        Returns a :class:`.DemographyDebugger` instance for this demography.

        :return: A DemographyDebugger object for this demography.
        :rtype: .DemographyDebugger
        """
        return DemographyDebugger(demography=self)

    def __eq__(self, other):
        try:
            self.assert_equal(other)
            return True
        except AssertionError:
            return False

    def assert_equal(self, other: Demography):
        """
        Compares this Demography with specified ``other`` and raises an
        AssertionError if they are not exactly equal.

        :param Demography other: The other demography to compare against.
        """
        # TODO we could potentially do better here with error messages
        # by showing a diff of the str() values for objects that differ.
        assert isinstance(other, Demography)
        assert self.num_populations == other.num_populations
        for p1, p2 in zip(self.populations, other.populations):
            assert p1 == p2, f"{p1} ≠ {p2}"
        assert np.array_equal(
            self.migration_matrix, other.migration_matrix
        )  # type: ignore
        assert self.num_events == other.num_events
        for e1, e2 in zip(self.events, other.events):
            assert e1 == e2, f"{e1} ≠ {e2}"

    def is_equivalent(self, other: Demography, rel_tol=None, abs_tol=None):
        """
        Compares this demography with the other and return True if they are
        equivalent up to the specified numerical tolerances. Two demographies
        are equivalent if, they have the same set of epochs defined by demographic
        events, and for each epoch:

        - The population's ``initial_size``, ``growth_rate`` and ``active``
          values are equal in all populations.
        - The migration matrices are equal
        - The same sequence of lineage movements through population splits, etc.

        All numerical comparisons are performed using :func:`py:math.isclose`.

        :param Demography other: The other demography to compare against.
        :param float rel_tol: The relative tolerance used by math.isclose.
        :param float abs_tol: The relative tolerance used by math.isclose.
        :return: True if this demography and other are equivalent up to numerical
            tolerances.
        :rtype bool: bool
        """
        try:
            self.assert_equivalent(other, rel_tol=rel_tol, abs_tol=abs_tol)
            return True
        except AssertionError:
            return False

    def assert_equivalent(
        self,
        other: Demography,
        rel_tol: Union[None, float] = None,
        abs_tol: Union[None, float] = None,
    ):
        # Same defaults as math.isclose
        rel_tol = 1e-9 if rel_tol is None else rel_tol
        abs_tol = 0 if abs_tol is None else abs_tol
        assert isinstance(other, Demography)
        self_dbg = self.debug()
        other_dbg = other.debug()
        if self.num_populations != other.num_populations:
            raise AssertionError(
                "Number of populations not equal: "
                f"{self.num_populations} ≠ {other.num_populations}"
            )
        # Compare the population attributes.
        # NB use the *resolved* versions from the debug objects
        for self_pop, other_pop in zip(
            self_dbg.demography.populations, other_dbg.demography.populations
        ):
            if self_pop.name != other_pop.name:
                raise AssertionError(
                    f"Population names differ: {self_pop.name} ≠ {other_pop.name}"
                )
            self_st = 0 if self_pop.sampling_time is None else self_pop.sampling_time
            other_st = 0 if other_pop.sampling_time is None else other_pop.sampling_time
            if not math.isclose(self_st, other_st):
                raise AssertionError(
                    f"Sampling times not equal for {self_pop.name}: "
                    f"{self_pop.sampling_time} ≠ {other_pop.sampling_time}"
                )

        if self_dbg.num_epochs != other_dbg.num_epochs:
            raise AssertionError(
                "Number of epochs not equal: "
                f"{self_dbg.num_epochs} ≠ {other_dbg.num_epochs}"
            )
        for j, (self_epoch, other_epoch) in enumerate(
            zip(self_dbg.epochs, other_dbg.epochs)
        ):
            if not math.isclose(
                self_epoch.start_time,
                other_epoch.start_time,
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            ):
                raise AssertionError(
                    f"Epoch[{j}] at different times: "
                    f"{self_epoch.start_time} ≠ {other_epoch.start_time}"
                )
            for self_pop, other_pop in zip(
                self_epoch.populations, other_epoch.populations
            ):
                if self_pop.state != other_pop.state:
                    raise AssertionError(
                        f"State mismatch in populations in epoch[{j}], {self_pop.name}: "
                        f"{self_pop.state} ≠ {other_pop.state}"
                    )
                if self_pop.state == PopulationStateMachine.ACTIVE:
                    if not math.isclose(
                        self_pop.start_size,
                        other_pop.start_size,
                        rel_tol=rel_tol,
                        abs_tol=abs_tol,
                    ):
                        raise AssertionError(
                            "Population start_size not equal to required precision "
                            f"in epoch[{j}], {self_pop.name}: "
                            f"{self_pop.start_size} ≠ {other_pop.start_size}"
                        )

                    if not math.isclose(
                        self_pop.growth_rate,
                        other_pop.growth_rate,
                        rel_tol=rel_tol,
                        abs_tol=abs_tol,
                    ):
                        raise AssertionError(
                            "Population growth_rate not equal to required precision "
                            f"in epoch[{j}], {self_pop.name}: "
                            f"{self_pop.growth_rate} ≠ {other_pop.growth_rate}"
                        )
            m_equal = np.isclose(
                self_epoch.migration_matrix,
                other_epoch.migration_matrix,
                rtol=rel_tol,
                atol=abs_tol,
            )
            if not np.all(m_equal):
                differs = "Differences between: "
                for source_id, dest_id in zip(*np.where(~m_equal)):
                    source = self[source_id].name
                    dest = self[dest_id].name
                    self_rate = self_epoch.migration_matrix[source_id, dest_id]
                    other_rate = other_epoch.migration_matrix[source_id, dest_id]
                    differs += f"({source}, {dest}: {self_rate} ≠ {other_rate}), "
                raise AssertionError(
                    f"Migration matrices in epoch[{j}] not equal: \n"
                    "self = \n"
                    f"{self_epoch.migration_matrix}\n"
                    "other = \n"
                    f"{other_epoch.migration_matrix}\n"
                    f"::{differs[:-2]}"
                )
            self._assert_lineage_movements_equivalent(
                self_epoch.events, other_epoch.events, rel_tol=rel_tol, abs_tol=abs_tol
            )

            self_state_change = [
                event
                for event in self_epoch.events
                if isinstance(event, StateChangeEvent)
            ]
            other_state_change = [
                event
                for event in other_epoch.events
                if isinstance(event, StateChangeEvent)
            ]
            if len(self_state_change) > 0 or len(other_state_change) > 0:
                raise ValueError(
                    "State change events not currently supported in equivalent. "
                    "Please open an issue on GitHub"
                )

    def _assert_lineage_movements_equivalent(
        self, self_events, other_events, *, rel_tol, abs_tol
    ):
        self_lineage_movements = self._normalise_lineage_movements(self_events)
        other_lineage_movements = self._normalise_lineage_movements(other_events)
        source_pops = set(self_lineage_movements.keys())
        if source_pops != set(other_lineage_movements.keys()):
            raise AssertionError(
                f"Mismatch in the set of populations affected by lineage movements: "
                f"{source_pops} ≠ {set(other_lineage_movements.keys())}"
            )
        for source in source_pops:
            self_out_movements = self_lineage_movements[source]
            other_out_movements = other_lineage_movements[source]
            if len(self_out_movements) != len(other_out_movements):
                raise AssertionError(
                    "Mismatch in number of normalised lineage movements out of "
                    f"{source}: {len(self_out_movements)} ≠ {len(other_out_movements)}"
                )
            for self_lm, other_lm in zip(self_out_movements, other_out_movements):
                if self_lm.dest != other_lm.dest:
                    raise AssertionError(
                        "Mismatch in lineage movement destination:"
                        f"{self_lm.dest} ≠ {other_lm.dest}"
                    )
                if not math.isclose(
                    self_lm.proportion,
                    other_lm.proportion,
                    rel_tol=rel_tol,
                    abs_tol=abs_tol,
                ):
                    raise AssertionError(
                        "Mismatch in normalised lineage movement proportions "
                        f"from {self_lm.source} to {self_lm.dest}: "
                        f"{self_lm.proportion} ≠ {other_lm.proportion}"
                    )

    def _normalise_lineage_movements(self, events: List[DemographicEvent]):
        """
        Extract the LineageMovementEvent instances from the specified list
        and normalise their effects into LineageMovement instances, and
        return a dictionary mapping source populations to the list of
        sequential lineage movements. For each source population we have a
        list of sequentual lineage movements, sorted by destination population ID.
        """
        assert len({event.time for event in events}) <= 1
        ret = collections.defaultdict(list)
        for event in events:
            if isinstance(event, LineageMovementEvent):
                for move in event._as_lineage_movements():
                    ret[move.source].append(move)
        for pop in list(ret.keys()):
            # We have a list of *conditional* lineage movements out of a
            # population. We canonicalise this by sorting by the
            # destination population, so we have to first convert back to
            # absolute proportions.
            assert all(lm.source == pop for lm in ret[pop])
            P = _sequential_to_proportions([pm.proportion for pm in ret[pop]])
            id_value_pairs = sorted([(lm.dest, p) for lm, p in zip(ret[pop], P)])
            S = _proportions_to_sequential([p for _, p in id_value_pairs])
            ret[pop] = [
                LineageMovement(source=pop, dest=id_value_pairs[j][0], proportion=S[j])
                for j in range(len(S))
            ]
        return ret

    @staticmethod
    def from_species_tree(
        tree,
        initial_size,
        *,
        time_units="gen",
        generation_time=None,
        growth_rate=None,
    ) -> Demography:
        """
        Parse a species tree in `Newick
        <https://en.wikipedia.org/wiki/Newick_format>`_ format and return the
        corresponding :class:`Demography` object. The tree is assumed to be
        rooted and ultrametric and branch lengths must be included and
        correspond to time, either in units of millions of years, years, or
        generations.

        The returned :class:`.Demography` object contains a
        :class:`.Population` for each node in the species tree. The
        population's ``name`` attribute will be either the corresponding node
        label from the newick tree, if it exists, or otherwise the name takes
        the form "pop_{j}", where j is the position of the given population in
        the list. Leaf populations are first in the list, and added in
        left-to-right order. Populations corresponding to the internal nodes
        are then added in a postorder traversal of the species tree. For each
        internal node a :class:`.PopulationSplit` event is added so that
        lineages move from its child populations at the appropriate time
        and rates of continuous migration to and from the child populations is
        set to zero. See the :ref:`sec_demography_events_population_split`
        section for more details.

        The initial sizes and growth rates for the populations in the model are
        set via the ``initial_size`` and ``growth_rate`` arguments. These can be
        specified in two ways: if a single number is provided, this is used
        for all populations. The argument may also be a mapping from population
        names to their respective values. For example:

        .. code-block:: python

            tree = "(A:10.0,B:10.0)C"
            initial_size = {"A": 1000, "B": 2000, "C": 100}
            demography = msprime.Demography.from_species_tree(tree, initial_size)

        Note that it is possible to have default population sizes for unnamed
        ancestral populations using a :class:`python:collections.defaultdict`, e.g.,

        .. code-block:: python

            tree = "(A:10.0,B:10.0)"
            initial_size = collections.defaultdict(lambda: 100)
            initial_size.update({"A": 1000, "B": 2000})
            demography = msprime.Demography.from_species_tree(tree, initial_size)

        :param str tree: The tree string in Newick format, with named leaves and branch
            lengths.
        :param initial_size: Each population's initial_size. May be a single number
            or a mapping from population names to their sizes.
        :param growth_rate: Each population's growth_rate. May be a single number
            or a mapping from population names to their exponential growth rates.
            Defaults to zero.
        :param str time_units: The units of time in which the species tree's
            branch lengths are measured. Allowed branch length units are millions of
            years, years, and generations; these should be specified with the strings
            ``"myr"``, ``"yr"``, or ``"gen"``, respectively. This defaults to
            ``"gen"``.
        :param float generation_time: The number of years per generation. If and only
            if the branch lengths are not in units of generations, the generation time
            must be specified. This defaults to `None`.
        :return: A Demography object representing the specified species tree.
        :rtype: .Demography
        """
        return species_trees.parse_species_tree(
            tree,
            initial_size=initial_size,
            growth_rate=growth_rate,
            time_units=time_units,
            generation_time=generation_time,
        )

    @staticmethod
    def from_starbeast(tree, generation_time, time_units="myr") -> Demography:
        """
        Parse a species tree produced by the program `TreeAnnotator
        <https://www.beast2.org/treeannotator>`_
        based on a posterior tree distribution generated with `StarBEAST
        <https://academic.oup.com/mbe/article/34/8/2101/3738283>`_  and return
        the corresponding Demography object.

        Species trees produced by TreeAnnotator are written in `Nexus
        <https://en.wikipedia.org/wiki/Nexus_file>`_ format and are rooted,
        bifurcating, and ultrametric. Branch lengths usually are in units of
        millions of years, but the use of other units is permitted by StarBEAST
        (and thus TreeAnnotator). This function allows branch length units of
        millions of years or years. Leaves must be named and the tree must
        include information on population sizes of leaf and ancestral species
        in the form of annotation with the "dmv" tag, which is the case for
        trees written by TreeAnnotator based on StarBEAST posterior tree
        distributions.

        The returned :class:`.Demography` object contains a
        :class:`.Population` for each node in the species tree. The
        population's ``name`` attribute will be either the corresponding node
        label from the newick tree, if it exists, or otherwise the name takes
        the form "pop_{j}", where j is the position of the given population in
        the list. Leaf populations are first in the list, and added in
        left-to-right order. Populations corresponding to the internal nodes
        are then added in a postorder traversal of the species tree. For each
        internal node a :class:`.PopulationSplit` event is added so that
        lineages move from its child populations at the appropriate time and
        rates of continuous migration to and from the child populations is set
        to zero. See the :ref:`sec_demography_events_population_split` section
        for more details.

        :param str tree: The tree string in Nexus format, with named leaves, branch
            lengths, and branch annotation. Typically, this string is the entire content
            of a file written by TreeAnnotator.
        :param float generation_time: The number of years per generation.
        :param str time_units: The units of time in which the species tree's
            branch lengths are measured. Allowed branch length units are millions of
            years, and years; these should be specified with the strings ``"myr"`` or
            ``"yr"``, respectively. This defaults to ``"myr"``.
        :return: A :class:`.Demography` instance that describing the information in the
            specified species tree.
        :rtype: .Demography
        """
        return species_trees.parse_starbeast(
            tree=tree,
            generation_time=generation_time,
            time_units=time_units,
        )

    def _from_old_style_map_populations(
        population_configurations: List[PopulationConfiguration],
        migration_matrix: List[List[float]],
        demographic_events: List[DemographicEvent],
        population_map: [List[Dict[int, str]]],
    ) -> Demography:
        direct_model = Demography._from_old_style_simple(
            population_configurations, migration_matrix, demographic_events
        )

        for id_map in population_map:
            if len(set(id_map.values())) != len(id_map):
                raise ValueError("Population IDs in old model must be unique")
            for old_id in id_map.values():
                if old_id < 0 or old_id >= direct_model.num_populations:
                    raise ValueError(
                        f"Bad population reference {old_id} in old style model"
                    )
        if len(population_map[0]) != len(population_configurations):
            raise ValueError(
                "The ID map for the first epoch must have entries for all "
                "populations in the old-style model"
            )

        # Suppress misspecification warnings; we'll give more specific messages
        # later.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dbg = direct_model.debug()

        if dbg.num_epochs != len(population_map):
            raise ValueError(
                "Mismatch in the number of epochs in the old style model "
                f"({dbg.num_epochs}) and specified in the population ID map "
                f"({len(population_map)})"
            )

        demography = Demography()
        # Set up the initial state of the model
        id_map = population_map[0]
        for new_name, old_id in id_map.items():
            pc = population_configurations[old_id]
            demography._add_population_from_old_style(pc, new_name)
        for epoch_id_map in population_map[1:]:
            for name in epoch_id_map.keys():
                if name not in demography:
                    demography.add_population(name=name, initial_size=0)

        # Fill out migration matrix
        for pop_1, pop_2 in itertools.combinations(id_map.keys(), 2):
            for pop_j, pop_k in [(pop_1, pop_2), (pop_2, pop_1)]:
                j = id_map[pop_j]
                k = id_map[pop_k]
                demography.set_migration_rate(
                    source=pop_j, dest=pop_k, rate=direct_model.migration_matrix[j, k]
                )

        # Set the state of populations and migration rates epoch by epoch.
        for epoch in dbg.epochs[1:]:
            epoch_id_map = population_map[epoch.index]
            # Set the population sizes and growth_rates for the extant populations.
            for new_id, old_id in epoch_id_map.items():
                pop_state = epoch.populations[old_id]
                demography.add_population_parameters_change(
                    epoch.start_time,
                    population=new_id,
                    initial_size=pop_state.start_size,
                    growth_rate=pop_state.growth_rate,
                )
            # Set the migration rates.
            for pop_1, pop_2 in itertools.combinations(epoch_id_map.keys(), 2):
                for pop_j, pop_k in [(pop_1, pop_2), (pop_2, pop_1)]:
                    j = epoch_id_map[pop_j]
                    k = epoch_id_map[pop_k]
                    demography.add_migration_rate_change(
                        time=epoch.start_time,
                        source=pop_j,
                        dest=pop_k,
                        rate=epoch.migration_matrix[j][k],
                    )

        # Add splits/admixtures/pulses according to the mass migrations in the
        # original model and ID maps.
        last_epoch_id_map = population_map[0]
        for epoch in dbg.epochs:
            logger.debug(
                f"Converting epoch[{epoch.index}] "
                f"{epoch.start_time:.2f}-{epoch.end_time:.2f}"
            )
            # New name -> old ID
            epoch_id_map = population_map[epoch.index]
            # Old ID -> new name
            old_id_map = {value: key for key, value in epoch_id_map.items()}
            last_epoch_pops = set(last_epoch_id_map.keys())
            epoch_pops = set(epoch_id_map.keys())
            derived = set()
            events = copy.deepcopy(epoch.events)
            old_derived_ids = []
            if last_epoch_pops != epoch_pops:
                ancestral = epoch_pops - last_epoch_pops
                derived = last_epoch_pops - epoch_pops
                old_derived_ids = []
                if len(ancestral) == 1:
                    logger.debug(
                        f"Adding split ancestral={ancestral} derived={derived}"
                    )
                    ancestral = ancestral.pop()
                    demography.add_population_split(
                        time=epoch.start_time,
                        derived=list(derived),
                        ancestral=ancestral,
                    )
                    old_derived_ids = [last_epoch_id_map[pop] for pop in derived]
                    derived_mass_migrations = 0
                    for event in events:
                        if isinstance(event, MassMigration):
                            if event.source in old_derived_ids:
                                derived_mass_migrations += 1
                            if event.proportion != 1:
                                raise ValueError(
                                    "MassMigration associated with population split "
                                    "with proportion != 1"
                                )
                    # This check is weak, but it will catch some errors at least.
                    # We're assuming the old model is well-formed, so this should
                    # help catch errors where the id map is slightly off.
                    if derived_mass_migrations < len(old_derived_ids) - 1:
                        raise ValueError(
                            "Insufficient MassMigrations found for population split"
                        )
                elif len(derived) == 1:
                    derived = derived.pop()
                    ancestral = []
                    sequential_proportions = []
                    old_derived_id = last_epoch_id_map[derived]
                    for event in events:
                        if event.source == old_derived_id:
                            ancestral.append(old_id_map[event.dest])
                            sequential_proportions.append(event.proportion)
                    proportions = _sequential_to_proportions(sequential_proportions)
                    if math.isclose(sum(proportions), 1):
                        if len(ancestral) == 1:
                            # This is a single split from the ancestral population,
                            # which continues to exist. We need to override the
                            # defaults for the behaviour we want here.
                            logger.debug(
                                f"Adding single split ancestral={ancestral[0]} "
                                f"derived={derived}"
                            )
                            pop = demography[ancestral[0]]
                            pop.initially_active = True
                            pop.sampling_time = None
                            demography.add_population_split(
                                time=epoch.start_time,
                                derived=[derived],
                                ancestral=ancestral[0],
                            )
                        else:
                            logger.debug(
                                f"Adding admixture ancestral={ancestral} "
                                f"derived={derived} proportions={proportions}"
                            )
                            demography.add_admixture(
                                time=epoch.start_time,
                                derived=derived,
                                ancestral=ancestral,
                                proportions=proportions,
                            )
                        old_derived_ids = [old_derived_id]
                    else:
                        raise ValueError(
                            "Admixture or single population split implied by ID map "
                            "but absolute population proportions don't sum to 1"
                        )
            for event in events:
                if isinstance(event, MassMigration):
                    if event.source not in old_derived_ids:
                        demography.add_mass_migration(
                            time=epoch.start_time,
                            source=old_id_map[event.source],
                            dest=old_id_map[event.dest],
                            proportion=event.proportion,
                        )
                elif not isinstance(
                    event, (MigrationRateChange, PopulationParametersChange)
                ):
                    raise ValueError(
                        "Only MassMigration, MigrationRateChange and "
                        "PopulationParametersChange events are supported"
                    )

            # Go through the inactive populations and check the migration
            # rates make sense. It's an error to migrate out of an inactive
            # population and a warning to migrate out of inactive pops.
            pairs = itertools.combinations(range(len(epoch.populations)), 2)
            active = set(epoch_id_map.values())
            for pop_1, pop_2 in pairs:
                for source, dest in [(pop_1, pop_2), (pop_2, pop_1)]:
                    if epoch.migration_matrix[source, dest] != 0:
                        if source in active and dest not in active:
                            raise ValueError(
                                "Non zero migration from an active population "
                                f"({source}) to inactive ({dest})"
                            )
                        if source not in active:
                            warnings.warn(
                                "Migration out of inactive population "
                                f"({source}) to ({dest}). This may be an "
                                "error in your model."
                            )

            last_epoch_id_map = epoch_id_map

        demography.sort_events()
        return demography

    @staticmethod
    def _from_old_style_simple(
        population_configurations=None,
        migration_matrix=None,
        demographic_events=None,
    ) -> Demography:
        """
        Creates a Demography object from the pre 1.0 style input parameters,
        reproducing the old semantics with respect to default values.
        """
        demography = Demography()
        for pop_config in population_configurations:
            demography._add_population_from_old_style(pop_config)
        if migration_matrix is not None:
            demography.migration_matrix = np.array(migration_matrix)
        if demographic_events is not None:
            for event in demographic_events:
                demography.add_event(copy.deepcopy(event))
        return demography

    @staticmethod
    def from_old_style(
        population_configurations=None,
        *,
        migration_matrix=None,
        demographic_events=None,
        Ne=1,
        ignore_sample_size=False,
        population_map: Union[[List[Dict[int, Union[str, int]]]], None] = None,
    ) -> Demography:
        """
        Creates a Demography object from the pre 1.0 style input parameters,
        reproducing the old semantics with respect to default values.

        No sample information is stored in the new-style :class:`.Demography`
        objects, and therefore if the ``sample_size`` attribute of any
        of the input :class:`.PopulationConfiguration` objects is set a
        ValueError will be raised by default. However, if the
        ``ignore_sample_size`` parameter is set to True, this check will
        not be performed and the sample sizes specified in the old-style
        :class:`.PopulationConfiguration` objects will be ignored.

        Please see the :ref:`sec_ancestry_samples` section for details on
        how to specify sample locations in :func:`.sim_ancestry`.

        .. todo:: Document the remaining parameters.
        """
        if population_configurations is None:
            pop_configs = [PopulationConfiguration(initial_size=Ne)]
        else:
            pop_configs = copy.deepcopy(population_configurations)
            for pop_config in pop_configs:
                if pop_config.initial_size is None:
                    pop_config.initial_size = Ne

                if pop_config.sample_size is not None and not ignore_sample_size:
                    raise ValueError(
                        "You have specified a `sample_size` in a "
                        "PopulationConfiguration object that is to be converted "
                        "into a new-style Demography object, "
                        "which does not contain any information about samples. "
                        "Please use the ``samples`` argument to sim_ancestry "
                        "instead, which provides flexible options for sampling "
                        "from different populations"
                    )

        if population_map is None:
            return Demography._from_old_style_simple(
                pop_configs, migration_matrix, demographic_events
            )
        else:
            return Demography._from_old_style_map_populations(
                pop_configs,
                migration_matrix,
                demographic_events,
                population_map,
            )

    @staticmethod
    def isolated_model(initial_size, *, growth_rate=None) -> Demography:
        """
        Returns a :class:`.Demography` object representing a collection of
        isolated populations with specified initial population sizes and
        growth rates. Please see :ref:`sec_demography` for more details on
        population sizes and growth rates.

        :param array_like initial_size: the ``initial_size`` value for each
            of the :class:`.Population` in the returned model. The length
            of the array corresponds to the number of populations.
            model.
        :param array_like growth_rate: The exponential growth rate for each
            population. Must be either None (the default, resulting a zero
            growth rate) or an array with the same length as ``initial_size``.
        :return: A Demography object representing this model, suitable as
            input to :func:`.simulate`.
        :rtype: .Demography
        """
        initial_size = np.array(initial_size, dtype=np.float64)
        if len(initial_size.shape) != 1:
            raise ValueError(
                "The initial_size argument must a 1D array of population size values"
            )
        if growth_rate is None:
            growth_rate = np.zeros_like(initial_size)
        else:
            growth_rate = np.array(growth_rate, dtype=np.float64)
        if initial_size.shape != growth_rate.shape:
            raise ValueError(
                "If growth_rate is specified it must be a 1D array of the same "
                "length as the population_size array"
            )
        if np.any(initial_size < 0):
            raise ValueError("population size values must be nonnegative.")
        if not np.all(np.isfinite(initial_size)):
            raise ValueError("population size values must be finite.")
        if not np.all(np.isfinite(growth_rate)):
            raise ValueError("growth_rate values must be finite.")
        populations = [
            Population(
                initial_size=initial_size[j],
                growth_rate=growth_rate[j],
            )
            for j in range(len(initial_size))
        ]
        return Demography(populations=populations)

    @staticmethod
    def island_model(initial_size, migration_rate, *, growth_rate=None) -> Demography:
        """
        Returns a :class:`.Demography` object representing a collection of
        populations with specified initial population sizes and growth
        rates, with symmetric migration between each pair of populations at the
        specified rate. Please see :ref:`sec_demography` for more details on
        population sizes and growth rates.

        :param array_like initial_size: the ``initial_size`` value for each
            of the :class:`.Population` in the returned model. The length
            of the array corresponds to the number of populations.
            model.
        :param float migration_rate: The migration rate between each pair of
            populations.
        :param array_like growth_rate: The exponential growth rate for each
            population. Must be either None (the default, resulting a zero
            growth rate) or an array with the same length as ``initial_size``.
        :return: A Demography object representing this model, suitable as
            input to :func:`.simulate`.
        :rtype: .Demography
        """
        model = Demography.isolated_model(initial_size, growth_rate=growth_rate)
        check_migration_rate(migration_rate)
        model.migration_matrix[:] = migration_rate
        np.fill_diagonal(model.migration_matrix, 0)
        return model

    @staticmethod
    def stepping_stone_model(
        initial_size, migration_rate, *, growth_rate=None, boundaries=False
    ) -> Demography:
        """
        Returns a :class:`.Demography` object representing a collection of
        populations with specified initial population sizes and growth
        rates, in which adjacent demes exchange migrants at the
        specified rate. Please see :ref:`sec_demography` for more details on
        population sizes and growth rates.

        .. note:: The current implementation only supports a one-dimensional
            stepping stone model, but higher dimensions could also be supported.
            Please open an issue on GitHub if this feature would be useful to you.

        :param array_like initial_size: the ``initial_size`` value for each
            of the :class:`.Population` in the returned model. The length
            of the array corresponds to the number of populations.
        :param float migration_rate: The migration rate between adjacent pairs
            of populations.
        :param array_like growth_rate: The exponential growth rate for each
            population. Must be either None (the default, resulting a zero
            growth rate) or an array with the same length as ``initial_size``.
        :param bool boundaries: If True the stepping stone model has boundary
            conditions imposed so that demes at either end of the chain do
            not exchange migrats. If False (the default), the set of
            populations is "circular" and migration takes place between the
            terminal demes.
        :return: A Demography object representing this model, suitable as
            input to :func:`.simulate`.
        :rtype: .Demography
        """
        initial_size = np.array(initial_size, dtype=np.float64)
        if len(initial_size.shape) > 1:
            raise ValueError(
                "Only 1D stepping stone models currently supported. Please open "
                "an issue on GitHub if you would like 2D (or more) models"
            )
        model = Demography.isolated_model(initial_size, growth_rate=growth_rate)
        check_migration_rate(migration_rate)
        if model.num_populations > 1:
            index1 = np.arange(model.num_populations, dtype=int)
            index2 = np.mod(index1 + 1, model.num_populations)
            model.migration_matrix[index1, index2] = migration_rate
            model.migration_matrix[index2, index1] = migration_rate
            if boundaries:
                model.migration_matrix[0, -1] = 0
                model.migration_matrix[-1, 0] = 0
        return model

    @staticmethod
    def _ooa_model():
        """
        Returns the Gutenkunst et al three population out-of-Africa model.

        This version is included here temporarily as a way to get some
        test coverage on the model compared with stdpopsim. Because we
        use this model in the documentation, we want make sure that it's
        doing what we think. We compare the model defined here then with
        the one presented in the docs, to ensure that no errors creep in.

        Once the upstream code in stdpopsim is updated to use msprime 1.0
        APIs we can remove this model and instead compare directly
        to the stdpopsim model with .is_equivalent() or whatever.
        """
        # Times are provided in years, so we convert into generations.
        generation_time = 25
        T_OOA = 21.2e3 / generation_time
        T_AMH = 140e3 / generation_time
        T_ANC = 220e3 / generation_time
        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        r_CEU = 0.004
        r_CHB = 0.0055
        N_CEU = 1000 / math.exp(-r_CEU * T_OOA)
        N_CHB = 510 / math.exp(-r_CHB * T_OOA)

        demography = Demography()
        demography.add_population(
            name="YRI",
            description="Yoruba in Ibadan, Nigeria",
            initial_size=12300,
        )
        demography.add_population(
            name="CEU",
            description=(
                "Utah Residents (CEPH) with Northern and Western European Ancestry"
            ),
            initial_size=N_CEU,
            growth_rate=r_CEU,
        )
        demography.add_population(
            name="CHB",
            description="Han Chinese in Beijing, China",
            initial_size=N_CHB,
            growth_rate=r_CHB,
        )
        demography.add_population(
            name="OOA",
            description="Bottleneck out-of-Africa population",
            initial_size=2100,
        )
        demography.add_population(
            name="AMH", description="Anatomically modern humans", initial_size=12300
        )
        demography.add_population(
            name="ANC",
            description="Ancestral equilibrium population",
            initial_size=7300,
        )

        # Set the migration rates between extant populations
        demography.set_symmetric_migration_rate(["CEU", "CHB"], 9.6e-5)
        demography.set_symmetric_migration_rate(["YRI", "CHB"], 1.9e-5)
        demography.set_symmetric_migration_rate(["YRI", "CEU"], 3e-5)

        demography.add_population_split(
            time=T_OOA, derived=["CEU", "CHB"], ancestral="OOA"
        )
        demography.add_symmetric_migration_rate_change(
            time=T_OOA, populations=["YRI", "OOA"], rate=25e-5
        )
        demography.add_population_split(
            time=T_AMH, derived=["YRI", "OOA"], ancestral="AMH"
        )
        demography.add_population_split(time=T_ANC, derived=["AMH"], ancestral="ANC")
        return demography

    @staticmethod
    def _ooa_archaic_model():
        """
        See notes for the _ooa model above.
        """
        # Implement the OutOfAfricaArchaicAdmixture_5R19 model
        # NOTE: this example isn't very well factored and needs more work.

        # Times are provided in years, so we convert into generations.
        generation_time = 29
        T_OOA = 36_000 / generation_time
        T_AMH = 60_700 / generation_time
        T_ANC = 300_000 / generation_time
        T_ArchaicAFR = 499_000 / generation_time
        T_Neanderthal = 559_000 / generation_time
        T_archaic_migration_start = 18_700 / generation_time
        T_archaic_migration_end = 125_000 / generation_time

        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        r_CEU = 0.00125
        r_CHB = 0.00372
        N_CEU = 2300 / math.exp(-r_CEU * T_OOA)
        N_CHB = 650 / math.exp(-r_CHB * T_OOA)

        demography = Demography()
        # This is the "trunk" population that we merge other populations into
        demography.add_population(
            name="AFR",
            description="African population",
            initial_size=13900,
            initially_active=True,
        )
        demography.add_population(
            name="CEU",
            description=(
                "Utah Residents (CEPH) with Northern and Western European Ancestry"
            ),
            initial_size=N_CEU,
            growth_rate=r_CEU,
        )
        demography.add_population(
            name="CHB",
            description="Han Chinese in Beijing, China",
            initial_size=N_CHB,
            growth_rate=r_CHB,
        )
        demography.add_population(
            name="Neanderthal",
            description="Putative Neanderthals",
            initial_size=3600,
        )
        demography.add_population(
            name="ArchaicAFR",
            description="Putative Archaic Africans",
            initial_size=3600,
        )
        demography.add_population(
            name="OOA",
            description="Bottleneck out-of-Africa population",
            initial_size=880,
        )

        # Set the migration rates between extant populations
        demography.set_symmetric_migration_rate(["CEU", "CHB"], 11.3e-5)
        demography.set_symmetric_migration_rate(["AFR", "CEU"], 2.48e-5)

        demography.add_symmetric_migration_rate_change(
            T_archaic_migration_start, ["CEU", "Neanderthal"], 0.825e-5
        )
        demography.add_symmetric_migration_rate_change(
            T_archaic_migration_start, ["CHB", "Neanderthal"], 0.825e-5
        )
        demography.add_symmetric_migration_rate_change(
            T_archaic_migration_start, ["ArchaicAFR", "AFR"], 1.98e-5
        )
        demography.add_migration_rate_change(T_archaic_migration_end, rate=0)

        demography.add_population_split(
            time=T_OOA, derived=["CEU", "CHB"], ancestral="OOA"
        )
        demography.add_symmetric_migration_rate_change(
            time=T_OOA, populations=["AFR", "OOA"], rate=52.2e-5
        )
        demography.add_symmetric_migration_rate_change(
            time=T_OOA, populations=["OOA", "Neanderthal"], rate=0.825e-5
        )
        demography.add_population_split(time=T_AMH, derived=["OOA"], ancestral="AFR")
        demography.add_symmetric_migration_rate_change(
            T_AMH, ["ArchaicAFR", "AFR"], 1.98e-5
        )
        demography.add_population_parameters_change(
            time=T_AMH, population="AFR", initial_size=13900
        )
        demography.add_population_parameters_change(
            time=T_ANC, population="AFR", initial_size=3600
        )
        demography.add_population_split(
            time=T_ArchaicAFR, derived=["ArchaicAFR"], ancestral="AFR"
        )
        demography.add_population_split(
            time=T_Neanderthal, derived=["Neanderthal"], ancestral="AFR"
        )
        demography.sort_events()
        return demography

    @staticmethod
    def _american_admixture_model():
        # Implementation of AmericanAdmixture_4B11 model. See notes from the _ooa
        # model above as to why this is here.
        T_OOA = 920
        N_EUR = 34039
        r_EUR = 0.0038
        N_EAS = 45852
        r_EAS = 0.0048
        T_ADMIX = 12
        N_ADMIX = 54664
        r_ADMIX = 0.05

        demography = Demography()
        demography.add_population(
            name="AFR", description="African population", initial_size=14474
        )
        demography.add_population(
            name="EUR",
            description="European population",
            initial_size=N_EUR,
            growth_rate=r_EUR,
        )
        demography.add_population(
            name="EAS",
            description="East Asian population",
            initial_size=N_EAS,
            growth_rate=r_EAS,
        )
        demography.add_population(
            name="ADMIX",
            description="Admixed America",
            initial_size=N_ADMIX,
            growth_rate=r_ADMIX,
        )
        demography.add_admixture(
            T_ADMIX,
            derived="ADMIX",
            ancestral=["AFR", "EUR", "EAS"],
            proportions=[1 / 6, 2 / 6, 3 / 6],
        )
        demography.add_population(
            name="OOA",
            description="Bottleneck out-of-Africa population",
            initial_size=1861,
        )
        demography.add_population(
            name="AMH", description="Anatomically modern humans", initial_size=14474
        )
        demography.add_population(
            name="ANC",
            description="Ancestral equilibrium population",
            initial_size=7310,
        )
        demography.set_symmetric_migration_rate(["AFR", "EUR"], 2.5e-5)
        demography.set_symmetric_migration_rate(["AFR", "EAS"], 0.78e-5)
        demography.set_symmetric_migration_rate(["EUR", "EAS"], 3.11e-5)

        demography.add_population_split(T_OOA, derived=["EUR", "EAS"], ancestral="OOA")
        demography.add_symmetric_migration_rate_change(
            time=T_OOA, populations=["AFR", "OOA"], rate=15e-5
        )
        demography.add_population_split(2040, derived=["OOA", "AFR"], ancestral="AMH")
        demography.add_population_split(5920, derived=["AMH"], ancestral="ANC")
        return demography


# This was lifted out of older code as-is. No point in updating it
# to use dataclasses, since all we want to do is maintain compatability
# with older code.
class PopulationConfiguration:
    """
    The initial configuration of a population (or deme) in a simulation.

    .. todo:: This is deprecated. Document as such.

    :param int sample_size: The number of initial samples that are drawn
        from this population.
    :param float initial_size: The absolute size of the population at time
        zero. Defaults to the reference population size :math:`N_e`.
    :param float growth_rate: The forwards-time exponential growth rate of the
        population per generation. Growth rates can be negative. This is zero for a
        constant population size, and positive for a population that has been
        growing. Defaults to 0.
    :param dict metadata: A JSON-encodable dictionary of metadata to associate
        with the corresponding Population in the output tree sequence.
        If not specified or None, no metadata is stored (i.e., an empty bytes array).
        Note that this metadata is ignored when using the ``from_ts`` argument to
        :func:`simulate`, as the population definitions in the tree sequence that
        is used as the starting point take precedence.
    """

    def __init__(
        self, sample_size=None, initial_size=None, growth_rate=0.0, metadata=None
    ):
        if initial_size is not None and initial_size < 0:
            raise ValueError("Population size must be >= 0")
        if sample_size is not None and sample_size < 0:
            raise ValueError("Sample size must be >= 0")
        self.sample_size = sample_size
        self.initial_size = initial_size
        self.growth_rate = growth_rate
        self.metadata = metadata

    def asdict(self):
        return dict(
            sample_size=self.sample_size,
            initial_size=self.initial_size,
            growth_rate=self.growth_rate,
            metadata=self.metadata,
        )


def _list_str(a: List, fmt=None):
    """
    Returns the specified items rendered as a string without quotes.
    """
    if fmt is None:
        joined = ", ".join(str(item) for item in a)
    else:
        joined = ", ".join(f"{item:{fmt}}" for item in a)
    return f"[{joined}]"


@dataclasses.dataclass
class DemographicEvent:
    """
    Superclass of demographic events that occur during simulations.
    """

    time: float
    demography: Demography = dataclasses.field(
        init=False, compare=False, default=None, repr=False
    )

    def _parameters(self):
        raise NotImplementedError()

    def _effect(self):
        raise NotImplementedError()

    def asdict(self):
        return {
            key: getattr(self, key)
            for key in inspect.signature(self.__init__).parameters.keys()
            if hasattr(self, key)
        }

    def _convert_id(self, population_ref):
        """
        Converts the specified population reference into an integer,
        suitable for input into the low-level code. We treat -1 as a special
        case because it's used as meaning "all populations" by the events.
        """
        if population_ref in [-1, None]:
            # Both of these mean "all populations"
            return -1
        if self.demography is None:
            # We need to be able to handle Events that are not associated with
            # a Demography to support old code. However, these should only ever
            # happen with integer IDs.
            if not core.isinteger(population_ref):
                raise ValueError(
                    "Working with demographic events not associated with a "
                    "Demography object is a legacy-only operation. Population "
                    "references must be integer IDs"
                )
            return population_ref
        return self.demography[population_ref].id


class ParameterChangeEvent(DemographicEvent):
    """
    Superclass of events that change some parameters in the underlying
    simulation model but don't actually affect the state in any other
    way.
    """


@dataclasses.dataclass
class PopulationParametersChange(ParameterChangeEvent):
    """
    Changes the demographic parameters of a population at a given time.

    This event generalises the ``-eg``, ``-eG``, ``-en`` and ``-eN``
    options from ``ms``. Note that unlike ``ms`` we do not automatically
    set growth rates to zero when the population size is changed.

    :param float time: The length of time ago at which this event
        occurred.
    :param float initial_size: The absolute diploid size of the population
        at the beginning of the time slice starting at ``time``. If None,
        this is calculated according to the initial population size and
        growth rate over the preceding time slice.
    :param float growth_rate: The new per-generation growth rate. If None,
        the growth rate is not changed. Defaults to None.
    :param int population: The ID of the population affected. If
        ``population`` is None, the changes affect all populations
        simultaneously.
    """

    initial_size: Union[float, None] = None
    growth_rate: Union[float, None] = None
    # TODO change the default to -1 to match MigrationRateChange.
    population: Union[int, None] = None
    # Deprecated.
    # TODO add a formal deprecation notice
    population_id: Union[int, None] = dataclasses.field(default=None, repr=False)

    _type_str: ClassVar[str] = "Population parameter change"

    def __post_init__(self):

        if self.population_id is not None and self.population is not None:
            raise ValueError(
                "population_id and population are aliases; cannot supply both."
            )

        if self.population_id is not None:
            self.population = self.population_id
        if self.growth_rate is None and self.initial_size is None:
            raise ValueError("Must specify one or more of growth_rate and initial_size")
        if self.initial_size is not None and self.initial_size < 0:
            raise ValueError("Cannot have a population size < 0")
        self.population = -1 if self.population is None else self.population

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.2 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        ret = {
            "type": "population_parameters_change",
            "time": self.time,
            "population": self._convert_id(self.population),
        }
        if self.growth_rate is not None:
            ret["growth_rate"] = self.growth_rate
        if self.initial_size is not None:
            ret["initial_size"] = self.initial_size
        return ret

    def _parameters(self):
        s = f"population={self.population}, "
        if self.initial_size is not None:
            s += f"initial_size={self.initial_size}, "
        if self.growth_rate is not None:
            s += f"growth_rate={self.growth_rate}, "
        return s[:-2]

    def _effect(self):
        s = ""
        if self.initial_size is not None:
            s += f"initial_size → {self.initial_size:.2g} "
            if self.growth_rate is not None:
                s += "and "
        if self.growth_rate is not None:
            s += f"growth_rate → {self.growth_rate:.3g} "
        s += "for"
        if self.population == -1:
            s += " all populations"
        else:
            s += f" population {self.population}"
        return s


@dataclasses.dataclass
class MigrationRateChange(ParameterChangeEvent):
    """
    Changes the rate of migration from one deme to another to a new value at a
    specific time. Migration rates are specified in terms of the rate at which
    lineages move from population ``source`` to ``dest`` during the progress of
    the simulation. Note that ``source`` and ``dest`` are from the perspective
    of the coalescent process; please see the :ref:`sec_ancestry_models`
    section for more details on the interpretation of this migration model.

    By default, ``source=-1`` and ``dest=-1``, which results in all
    non-diagonal elements of the migration matrix being changed to the new
    rate. If ``source`` and ``dest`` are specified, they must refer to valid
    population IDs.

    :param float time: The time at which this event occurs in generations.
    :param float rate: The new per-generation migration rate.
    :param int source: The ID of the source population.
    :param int dest: The ID of the destination population.
    :param int source: The source population ID.
    """

    rate: float
    source: int = -1
    dest: int = -1
    # Deprecated.
    # TODO add a formal deprecation notice
    matrix_index: Union[tuple, None] = dataclasses.field(default=None, repr=False)

    _type_str: ClassVar[str] = "Migration rate change"

    def __post_init__(self):
        # If the deprecated form is used, it overwrites the values of source
        # and dest
        if self.matrix_index is not None:
            self.source = self.matrix_index[0]
            self.dest = self.matrix_index[1]

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "migration_rate_change",
            "time": self.time,
            # Note: We'd like to change the name here to "rate" but it's best
            # to leave this alone until stdpopsim has been moved away from
            # using this internal API.
            "migration_rate": self.rate,
            "source": self._convert_id(self.source),
            "dest": self._convert_id(self.dest),
        }

    def _parameters(self):
        return f"source={self.source}, dest={self.dest}, rate={self.rate}"

    def _effect(self):
        ret = "Backwards-time migration rate "
        if self.source == -1 and self.dest == -1:
            ret += "for all populations "
        else:
            ret += f"from {self.source} to {self.dest} "
        ret += f"→ {self.rate}"
        return ret


# TODO not clear we want to document this as part of the external API.
@dataclasses.dataclass
class SymmetricMigrationRateChange(ParameterChangeEvent):
    """
    Sets the symmetric migration rate between all pairs of populations in
    the specified list to the specified value. For a given pair of population
    IDs ``j`` and ``k``, this sets ``migration_matrix[j, k] = rate``
    and ``migration_matrix[k, j] = rate``.

    Populations may be specified either by their integer IDs or by
    their string names.

    :param float time: The time at which this event occurs in generations.
    :param list populations: An iterable of population identifiers (integer
        IDs or string names).
    :param float rate: The new migration rate.
    """

    populations: List[Union[int, str]]
    rate: float

    _type_str: ClassVar[str] = dataclasses.field(
        default="Symmetric migration rate change", repr=False
    )

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "symmetric_migration_rate_change",
            "time": self.time,
            "populations": [self._convert_id(pop) for pop in self.populations],
            "rate": self.rate,
        }

    def _parameters(self):
        return f"populations={_list_str(self.populations)}, rate={self.rate}"

    def _effect(self):
        s = "Sets the symmetric migration rate between "
        if len(self.populations) == 2:
            s += f"{self.populations[0]} and {self.populations[1]} "
        else:
            s += f"all pairs of populations in {_list_str(self.populations)} "
        s += f"to {self.rate} per generation"
        return s


@dataclasses.dataclass
class LineageMovement:
    """
    A single instantaneous movement of lineages from one population to
    another. Note that 'source' and 'dest' are in the backwards-in-time
    sense.
    """

    source: int
    dest: int
    proportion: float


class LineageMovementEvent(DemographicEvent):
    """
    Superclass of events that move lineages around between populations.
    """

    def _as_lineage_movements(self) -> List[LineageMovement]:
        """
        Returns the equivalent of this lineage movement event as a
        list of lineage movements.
        """
        raise NotImplementedError()


@dataclasses.dataclass
class MassMigration(LineageMovementEvent):
    """
    A mass migration event in which some fraction of the population in one deme
    (the ``source``) simultaneously move to another deme (``dest``) during the
    progress of the simulation. Each lineage currently present in the source
    population moves to the destination population with probability equal to
    ``proportion``. Note that ``source`` and ``dest`` are from the perspective
    of the coalescent process; please see the :ref:`sec_ancestry_models`
    section for more details on the interpretation of this migration model.

    This event class generalises the population split (``-ej``) and
    admixture (``-es``) events from ``ms``. Note that MassMigrations
    do *not* have any side effects on the migration matrix.

    :param float time: The time at which this event occurs in generations.
    :param int source: The ID of the source population.
    :param int dest: The ID of the destination population.
    :param float proportion: The probability that any given lineage within
        the source population migrates to the destination population.
    """

    source: int
    # dest only has a default because of the deprecated destination attr.
    dest: Union[None, int] = None
    proportion: float = 1.0
    # Deprecated.
    # TODO add a formal deprecation notice
    destination: Union[int, None] = dataclasses.field(default=None, repr=False)

    _type_str: ClassVar[str] = dataclasses.field(default="Mass Migration", repr=False)

    def __post_init__(self):
        if self.dest is not None and self.destination is not None:
            raise ValueError("dest and destination are aliases; cannot supply both")
        if self.destination is not None:
            self.dest = self.destination

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "mass_migration",
            "time": self.time,
            "source": self._convert_id(self.source),
            "dest": self._convert_id(self.dest),
            "proportion": self.proportion,
        }

    def _parameters(self):
        return (
            f"source={self.source}, dest={self.dest}, proportion={self.proportion:.3g}"
        )

    def _effect(self):
        if self.proportion == 1.0:
            ret = (
                f"All lineages currently in population {self.source} move "
                f"to {self.dest} "
            )
        else:
            ret = (
                f"Lineages currently in population {self.source} move to {self.dest} "
                f"with probability {self.proportion:.3g} "
            )
        ret += (
            "(equivalent to individuals "
            f"migrating from {self.dest} to {self.source} forwards in time)"
        )
        return ret

    def _as_lineage_movements(self):
        return [
            LineageMovement(
                source=self._convert_id(self.source),
                dest=self._convert_id(self.dest),
                proportion=self.proportion,
            )
        ]


@dataclasses.dataclass
class PopulationSplit(LineageMovementEvent):
    """

    :param float time: The time at which this event occurs in generations.
    :param list(int) derived: The ID(s) of the derived population(s).
    :param int ancestral: The ID of the ancestral population.
    """

    derived: List[Union[int, str]]
    ancestral: Union[int, str]

    _type_str: ClassVar[str] = dataclasses.field(default="Population Split", repr=False)

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "population_split",
            "time": self.time,
            "derived": [self._convert_id(pop) for pop in self.derived],
            "ancestral": self._convert_id(self.ancestral),
        }

    def _parameters(self):
        return f"derived={_list_str(self.derived)}, ancestral={self.ancestral}"

    def _effect(self):
        s = "Moves all lineages from "
        if len(self.derived) == 1:
            s += f"the '{self.derived[0]}' derived population "
        else:
            s += "derived populations "
            if len(self.derived) == 2:
                s += f"'{self.derived[0]}' and '{self.derived[1]}' "
            else:
                s += f"{_list_str(self.derived)} "
        s += f"to the ancestral '{self.ancestral}' population. "
        s += "Also set "
        if len(self.derived) == 1:
            s += f"'{self.derived[0]}' "
        else:
            s += "the derived populations "
        s += (
            "to inactive, and all migration rates to and from "
            f"the derived population{'s' if len(self.derived) > 1 else ''} "
            "to zero."
        )

        return s

    def _as_lineage_movements(self):
        ancestral = self._convert_id(self.ancestral)
        return [
            LineageMovement(source=self._convert_id(pop), dest=ancestral, proportion=1)
            for pop in self.derived
        ]


@dataclasses.dataclass
class Admixture(LineageMovementEvent):
    """

    :param float time: The time at which this event occurs in generations.
    """

    derived: Union[int, str]
    ancestral: List[Union[int, str]]
    proportions: List[float]

    _type_str: ClassVar[str] = dataclasses.field(default="Admixture", repr=False)

    @property
    def num_ancestral(self):
        return len(self.ancestral)

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "admixture",
            "time": self.time,
            "derived": self._convert_id(self.derived),
            "ancestral": [self._convert_id(pop) for pop in self.ancestral],
            "proportions": self.proportions,
        }

    def _parameters(self):
        return (
            f"derived={self.derived} ancestral={_list_str(self.ancestral)} "
            f"proportions={_list_str(self.proportions, '.2f')}"
        )

    def _effect(self):
        move_to = "; ".join(
            f"'{pop}' with proba {proba:.3g}"
            for pop, proba in zip(self.ancestral, self.proportions)
        )
        return (
            f"Moves all lineages from admixed population '{self.derived}' "
            f"to ancestral population{'s' if len(self.ancestral) > 1 else ''}. "
            f"Lineages move to {move_to}. Set '{self.derived}' to inactive, "
            f"and all migration rates to and from '{self.derived}' to zero."
        )

    def _as_lineage_movements(self):
        derived = self._convert_id(self.derived)
        ancestral = [self._convert_id(pop) for pop in self.ancestral]
        # Conditioned on having already distributed a fraction q of the
        # lineages, the we need a fraction p / (1 - q) of the remaining
        # lineages to get an overall proportion of p.
        S = _proportions_to_sequential(self.proportions)
        return [
            LineageMovement(source=derived, dest=ancestral[j], proportion=S[j])
            for j in range(self.num_ancestral)
        ]


class StateChangeEvent(DemographicEvent):
    """
    Superclass of events that change the state of the simulation in complex
    ways.
    """


# This is an unsupported/undocumented demographic event.
@dataclasses.dataclass
class SimpleBottleneck(StateChangeEvent):
    population: int
    proportion: float = 1.0

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "simple_bottleneck",
            "time": self.time,
            "population": self._convert_id(self.population),
            "proportion": self.proportion,
        }

    _type_str: ClassVar[str] = dataclasses.field(
        default="Simple Bottleneck", repr=False
    )

    def _parameters(self):
        return f"population={self.population}, proportion={self.proportion}"

    def _effect(self):
        return (
            f"Lineages in population {self.population} coalesce with "
            f"probability {self.proportion}"
        )


# TODO document
@dataclasses.dataclass
class InstantaneousBottleneck(StateChangeEvent):
    population: int
    strength: float = 1.0

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "instantaneous_bottleneck",
            "time": self.time,
            "population": self._convert_id(self.population),
            "strength": self.strength,
        }

    _type_str: ClassVar[str] = dataclasses.field(
        default="Instantaneous Bottleneck", repr=False
    )

    def _parameters(self):
        return f"population={self.population}, strength={self.strength}"

    def _effect(self):
        return f"Equivalent to {self.strength} generations of the coalescent"


@dataclasses.dataclass
class CensusEvent(DemographicEvent):
    """
    An event that adds a node to each branch of every tree at a given time
    during the simulation. This may be used to record all ancestral haplotypes
    present at that time, and to extract other information related to these
    haplotypes: for instance to trace the local ancestry of a sample back to a
    set of contemporaneous ancestors, or to assess whether a subset of samples
    has coalesced more recently than the census time.

    See :ref:`sec_ancestry_census_events` for more details.

    :param float time: The time at which this event occurs in generations.
    """

    _type_str: ClassVar[str] = dataclasses.field(default="Census", repr=False)

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "census_event",
            "time": self.time,
        }

    def _parameters(self):
        return ""

    def _effect(self):
        return "Insert census nodes to record the location of all lineages"


def _sequential_to_proportions(S):
    """
    Given a list of sequential lineage proportions out of a population,
    return the absolute proportions of the original population this
    corresponds to.
    """
    P = []
    for j in range(len(S)):
        P.append(S[j] * (1 - sum(P[:j])))
    return P


def _proportions_to_sequential(P):
    """
    Given a list of absolute proportions of lineages moving out of
    a population, return the sequential conditional movements required
    to give them same proportions.
    """
    # Conditioned on having already distributed a fraction q of the
    # lineages, the we need a fraction p / (1 - q) of the remaining
    # lineages to get an overall proportion of p
    C = [P[j] / (1 - sum(P[:j])) for j in range(len(P))]
    return C


def _matrix_exponential(A):
    """
    Returns the matrix exponential of A.
    https://en.wikipedia.org/wiki/Matrix_exponential
    Note: this is not a general purpose method and is only intended for use within
    msprime.
    """
    d, Y = np.linalg.eig(A)
    Yinv = np.linalg.pinv(Y)
    D = np.diag(np.exp(d))
    B = np.matmul(Y, np.matmul(D, Yinv))
    return np.real_if_close(B, tol=1000)


class PopulationStateMachine(enum.IntEnum):
    """
    During a simulation each population has three possible states described
    by this state machine. In general, a population follows:

    INACTIVE -> ACTIVE -> PREVIOUSLY_ACTIVE

    All populations are by default ACTIVE at the start of the simulation,
    except if they are they are "ancestral" in a population split event.
    In this case populations are initially INACTIVE. A population
    then transitions from INACTIVE -> ACTIVE when the corresponding
    population split event occurs.

    Populations transition from ACTIVE -> PREVIOUSLY_ACTIVE when they
    are "derived" in either population split or admixture events.

    No other transitions are possible.
    """

    INACTIVE = 0
    ACTIVE = 1
    PREVIOUSLY_ACTIVE = 2


@dataclasses.dataclass
class PopulationState:
    """
    Simple class to represent the state of a population in terms of its
    demographic parameters. Note: start and end here refer to time flowing
    *backwards*!
    """

    id: int  # noqa: A003
    name: str
    start_size: float
    end_size: float
    growth_rate: float
    state: int

    @property
    def active(self):
        return self.state == PopulationStateMachine.ACTIVE


@dataclasses.dataclass
class Epoch:
    """
    Represents a single epoch in the simulation within which the state
    of the demographic parameters are constant.
    """

    index: int
    start_time: float
    end_time: float
    populations: List[PopulationState]
    migration_matrix: list  # TODO numpy array
    events: List[DemographicEvent]

    def _title_text(self):
        return (
            f"Epoch[{self.index}]: "
            f"[{self.start_time:.3g}, {self.end_time:.3g}) generations"
        )

    def _population_state_text(self):
        return (
            f"Populations "
            f"(total={len(self.populations)} active={self.num_active_populations})"
        )

    @property
    def demographic_events(self):
        # For compatibility with msprime 0.x
        return self.events

    @property
    def active_populations(self):
        return [pop for pop in self.populations if pop.active]

    @property
    def num_active_populations(self):
        return len(self.active_populations)


class DemographyDebugger:
    """
    Utilities to compute and display information about the state of populations
    during the different simulation epochs defined by demographic events.

    .. warning:: This class is not intended to be instantiated directly using
        the contructor - please use :meth:`.Demography.debug()` to obtain
        a DemographyDebugger for a given :class:`.Demography` instead.
    """

    def __init__(
        self,
        # Deprecated pre-1.0 parameters.
        Ne=1,
        population_configurations=None,
        migration_matrix=None,
        demographic_events=None,
        model=None,
        *,
        demography=None,
    ):
        if demography is None:
            # Support the pre-1.0 syntax
            demography = Demography.from_old_style(
                population_configurations,
                migration_matrix=migration_matrix,
                demographic_events=demographic_events,
                Ne=Ne,
                ignore_sample_size=True,
            )
        self.demography = demography.validate()
        self.num_populations = demography.num_populations
        self._make_epochs()
        self._check_misspecification()

    def _make_epochs(self):
        self.epochs = []
        simulator = ancestry._parse_sim_ancestry(
            demography=self.demography, init_for_debugger=True
        )
        start_time = 0
        end_time = 0
        abs_tol = 1e-9
        event_index = 0
        all_events = self.demography.events
        while not math.isinf(end_time):
            events = []
            while event_index < len(all_events) and math.isclose(
                all_events[event_index].time, start_time, abs_tol=abs_tol
            ):
                events.append(all_events[event_index])
                event_index += 1
            end_time = simulator.debug_demography()
            migration_matrix = simulator.migration_matrix
            pop_conf = simulator.population_configuration
            populations = [
                PopulationState(
                    id=j,
                    name=self.demography.populations[j].name,
                    start_size=simulator.compute_population_size(j, start_time),
                    end_size=simulator.compute_population_size(j, end_time),
                    growth_rate=pop_conf[j]["growth_rate"],
                    state=PopulationStateMachine(pop_conf[j]["state"]),
                )
                for j in range(self.num_populations)
            ]
            epoch_index = len(self.epochs)
            self.epochs.append(
                Epoch(
                    epoch_index,
                    start_time,
                    end_time,
                    populations,
                    migration_matrix,
                    events,
                )
            )
            start_time = end_time

    def _check_misspecification(self):
        """
        Check for things that might indicate model misspecification.
        """
        merged_pops = set()
        for epoch in self.epochs:
            for de in epoch.events:
                if isinstance(de, MassMigration) and de.proportion == 1:
                    merged_pops.add(de.source)
            mm = epoch.migration_matrix
            for pop_k in merged_pops:
                k = self.demography[pop_k].id
                if any(mm[k, :] != 0) or any(mm[:, k] != 0):
                    warnings.warn(
                        "Non-zero migration rates exist after merging "
                        f"population {k}. This almost certainly indicates "
                        "demographic misspecification."
                    )

    def _populations_table(self, epoch, as_text=True):
        active_populations = epoch.active_populations
        column_titles = ["", "start", "end", "growth_rate"]
        if len(active_populations) > 1:
            column_titles += [pop.name for pop in active_populations]
        data = []
        for pop in active_populations:
            row = [
                pop.name,
                f"{pop.start_size: .1f}",
                f"{pop.end_size: .1f}",
                f"{pop.growth_rate: .3g}",
            ]
            if len(active_populations) > 1:
                for other_pop in active_populations:
                    item = self.demography._migration_rate_info(
                        pop.id,
                        other_pop.id,
                        epoch.migration_matrix[pop.id, other_pop.id],
                    )
                    if as_text:
                        row.append(item.as_text())
                    else:
                        row.append(item)

            data.append(row)
        return column_titles, data

    def _populations_html(self, epoch):
        column_titles, data = self._populations_table(epoch, as_text=False)
        return core.html_table(epoch._population_state_text(), column_titles, data)

    def _populations_text(self, epoch):
        column_titles, data = self._populations_table(epoch)
        alignments = ">>><"
        if epoch.num_active_populations > 1:
            alignments += "^" * epoch.num_active_populations
        # Repack the table items as lists
        column_titles = [[x] for x in column_titles]
        data = [[[x] for x in row] for row in data]
        return core.text_table(
            epoch._population_state_text(), column_titles, alignments, data
        )

    def _repr_html_(self):
        out = ""
        for epoch in self.epochs:
            if epoch.index > 0:
                assert len(epoch.events) > 0
                title = f"Events @ generation {epoch.start_time:.3g}"
                out += self.demography._events_html(epoch.events, title)
                out += "</div></details>"
            else:
                assert len(epoch.events) == 0
            title = epoch._title_text()
            out += f'<details open="true"><summary>{title}</summary>'
            # Indent the content div slightly
            out += '<div style="margin-left:20px">'
            out += self._populations_html(epoch)
        out += "</div>"
        out += "</details>"
        return f"<div>{out}</div>"

    def print_history(self, output=sys.stdout):
        """
        Prints a summary of the history of the populations.

        Deprecated since 1.0: use ``print(debugger)`` instead.
        """
        print(self, file=output, end="")

    def __str__(self):
        def indent(table, header_char="╟", depth=4):
            lines = table.splitlines()
            s = header_char + (" " * depth) + lines[0] + "\n"
            for line in lines[1:]:
                s += "║" + (" " * depth) + line + "\n"
            return s

        def box(title):
            N = len(title) + 2
            top = "╠" + ("═" * N) + "╗"
            bottom = "╠" + ("═" * N) + "╝"
            return f"{top}\n║ {title} ║\n{bottom}\n"

        out = "DemographyDebugger\n"
        for epoch in self.epochs:
            if epoch.index > 0:
                assert len(epoch.events) > 0
                title = f"Events @ generation {epoch.start_time:.3g}"
                out += indent(self.demography._events_text(epoch.events, title))
            out += box(epoch._title_text())
            out += indent(self._populations_text(epoch))
        return out

    def population_size_trajectory(self, steps):
        """
        This function returns an array of per-population effective population sizes,
        as defined by the demographic model. These are the `initial_size`
        parameters of the model, modified by any population growth rates.
        The sizes are computed at the time points given by `steps`.

        :param list steps: List of times ago at which the population
            size will be computed.
        :return: Returns a numpy array of population sizes, with one column per
            population, whose [i,j]th entry is the size of population
            j at time steps[i] ago.
        """
        num_pops = self.num_populations
        N_t = np.zeros([len(steps), num_pops])
        for j, t in enumerate(steps):
            N, _ = self._pop_size_and_migration_at_t(t)
            N_t[j] = N
        return N_t

    def lineage_probabilities(self, steps, sample_time=0):
        """
        Returns an array such that P[j, a, b] is the probability that a lineage that
        started in population a at time sample_time is in population b at time steps[j]
        ago.

        This function reports sampling probabilities _before_ mass migration events
        at a step time, if a mass migration event occurs at one of those times.
        Migrations will then effect the next time step.

        :param list steps: A list of times to compute probabilities.
        :param sample_time: The time of sampling of the lineage. For any times in steps
            that are more recent than sample_time, the probability of finding the
            lineage in any population is zero.
        :return: An array of dimension len(steps) by num pops by num_pops.
        """
        num_pops = self.num_populations
        # P[i, j] will be the probability that a lineage that started in i is now in j
        P = np.eye(num_pops)

        # epochs are defined by mass migration events or changes to population sizes
        # or migration rates, so we add the epoch interval times to the steps that we
        # need to account for
        epoch_breaks = [t for t in self.epoch_times if t not in steps]
        all_steps = np.concatenate([steps, epoch_breaks])

        sampling = []
        if sample_time not in all_steps:
            sampling.append(sample_time)
        all_steps = np.concatenate((all_steps, sampling))

        ix = np.argsort(all_steps)
        all_steps = all_steps[ix]
        # keep track of the steps to report in P_out
        keep_steps = np.concatenate(
            [
                np.repeat(True, len(steps)),
                np.repeat(False, len(epoch_breaks)),
                np.repeat(False, len(sampling)),
            ]
        )[ix]

        assert len(np.unique(all_steps)) == len(all_steps)
        assert np.all(steps == all_steps[keep_steps])
        P_out = np.zeros((len(all_steps), num_pops, num_pops))

        first_step = 0
        while all_steps[first_step] < sample_time:
            first_step += 1

        P_out[first_step] = P

        # get ordered mass migration events
        mass_migration_objects = []
        mass_migration_times = []
        for demo in self.demography.events:
            if isinstance(demo, MassMigration):
                mass_migration_objects.append(demo)
                mass_migration_times.append(demo.time)

        for jj in range(first_step, len(all_steps) - 1):
            t_j = all_steps[jj]

            # apply any mass migration events to P
            # so if we sample at this time, we do no account for the instantaneous
            # mass migration events that occur at the same time. that will show up
            # at the next step
            if t_j > sample_time:
                for mass_mig_t, mass_mig_e in zip(
                    mass_migration_times, mass_migration_objects
                ):
                    if mass_mig_t == t_j:
                        S = np.eye(num_pops, num_pops)
                        S[mass_mig_e.source, mass_mig_e.dest] = mass_mig_e.proportion
                        S[mass_mig_e.source, mass_mig_e.source] = (
                            1 - mass_mig_e.proportion
                        )
                        P = np.matmul(P, S)

            # get continuous migration matrix over next interval
            _, M = self._pop_size_and_migration_at_t(t_j)
            dt = all_steps[jj + 1] - all_steps[jj]
            dM = np.diag([sum(s) for s in M])
            # advance to next interval time (dt) taking into account continuous mig
            P = P.dot(_matrix_exponential(dt * (M - dM)))
            P_out[jj + 1] = P

        return P_out[keep_steps]

    def possible_lineage_locations(self, samples=None):
        """
        Given the sampling configuration, this function determines when lineages are
        possibly found within each population over epochs defined by demographic events
        and sampling times. If no sampling configuration is given, we assume we sample
        lineages from every population at time zero. The samples are specified by a list
        of msprime Sample objects, so that possible ancient samples may be accounted for.

        :param list samples: A list of msprime Sample objects, which specify their
            populations and times.
        :return: Returns a dictionary with epoch intervals as keys whose values are a
            list with length equal to the number of populations with True and False
            indicating which populations could possibly contain lineages over that
            epoch. The epoch intervals are given by tuples: (epoch start, epoch end).
            The first epoch necessarily starts at time 0, and the final epoch has end
            time of infinity.
        """
        # get configuration of sampling times from samples ({time:[pops_sampled_from]})
        if samples is None:
            sampling_times = {0: [i for i in range(self.num_populations)]}
        else:
            sampling_times = collections.defaultdict(list)
            for sample in samples:
                sampling_times[sample.time].append(sample.population)
            for t in sampling_times.keys():
                sampling_times[t] = list(set(sampling_times[t]))

        all_steps = sorted(
            list(set([t for t in self.epoch_times] + list(sampling_times.keys())))
        )

        epochs = [(x, y) for x, y in zip(all_steps[:-1], all_steps[1:])]
        epochs.append((all_steps[-1], np.inf))

        # need to go a bit beyond last step and into the final epoch that extends to inf
        all_steps.append(all_steps[-1] + 1)

        indicators = {e: np.zeros(self.num_populations, dtype=bool) for e in epochs}
        for sample_time, demes in sampling_times.items():
            P_out = self.lineage_probabilities(all_steps, sample_time=sample_time)
            for epoch, P in zip(epochs, P_out[1:]):
                if epoch[1] <= sample_time:
                    # samples shouldn't affect the epoch previous to the sampling time
                    continue
                for deme in demes:
                    indicators[epoch][P[deme] > 0] = True

        # join epochs if adjacent epochs have same set of possible live populations
        combined_indicators = {}
        skip = 0
        for ii, (epoch, inds) in enumerate(indicators.items()):
            if skip > 0:
                skip -= 1
                continue
            this_epoch = epoch
            while ii + skip + 1 < len(epochs) and np.all(
                indicators[epochs[ii + 1 + skip]] == inds
            ):
                this_epoch = (this_epoch[0], epochs[ii + 1 + skip][1])
                skip += 1
            combined_indicators[this_epoch] = inds

        return combined_indicators

    def mean_coalescence_time(
        self, num_samples, min_pop_size=1, steps=None, rtol=0.005, max_iter=12
    ):
        """
        Compute the mean time until coalescence between lineages of two samples drawn
        from the sample configuration specified in `num_samples`. This is done using
        :meth:`coalescence_rate_trajectory
        <.DemographyDebugger.coalescence_rate_trajectory>`
        to compute the probability that the lineages have not yet coalesced by time `t`,
        and using these to approximate :math:`E[T] = \\int_t^\\infty P(T > t) dt`,
        where :math:`T` is the coalescence time. See
        :meth:`coalescence_rate_trajectory
        <.DemographyDebugger.coalescence_rate_trajectory>`
        for more details.

        To compute this, an adequate time discretization must be arrived at
        by iteratively extending or refining the current discretization.
        Debugging information about numerical convergence of this procedure is
        logged using the Python :mod:`logging` infrastructure. To make it appear, using
        the :mod:`daiquiri` module, do for instance::

            import daiquiri

            daiquiri.setup(level="DEBUG")
            debugger.mean_coalescence_time([2])

        will print this debugging information to stderr. Briefly, this outputs
        iteration number, mean coalescence time, maximum difference in probabilty
        of not having coalesced yet, difference to last coalescence time,
        probability of not having coalesced by the final time point, and
        whether the last iteration was an extension or refinement.

        :param list num_samples: A list of the same length as the number
            of populations, so that `num_samples[j]` is the number of sampled
            chromosomes in subpopulation `j`.
        :param int min_pop_size: See :meth:`coalescence_rate_trajectory
            <.DemographyDebugger.coalescence_rate_trajectory>`.
        :param list steps: The time discretization to start out with (by default,
            picks something based on epoch times).
        :param float rtol: The relative tolerance to determine mean coalescence time
            to (used to decide when to stop subdividing the steps).
        :param int max_iter: The maximum number of times to subdivide the steps.
        :return: The mean coalescence time (a number).
        :rtype: float
        """

        def mean_time(steps, P):
            # Mean is int_0^infty P(T > t) dt, which we estimate by discrete integration
            # assuming that f(t) = P(T > t) is piecewise exponential:
            # if f(u) = a exp(bu) then b = log(f(t)/f(s)) / (t-s) for each s < t, so
            # \int_s^t f(u) du = (a/b) \int_s^t exp(bu) b du = (a/b)(exp(bt) - exp(bs))
            #    = (t - s) * (f(t) - f(s)) / log(f(t) / f(s))
            # unless b = 0, of course.
            assert steps[0] == 0
            dt = np.diff(steps)
            dP = np.diff(P)

            with np.errstate(divide="ignore", invalid="ignore"):
                dlogP = np.diff(np.log(P))
            nz = np.logical_and(dP < 0, P[1:] * P[:-1] > 0)
            const = dP == 0
            return np.sum(dt[const] * (P[:-1])[const]) + np.sum(
                dt[nz] * dP[nz] / dlogP[nz]
            )

        if steps is None:
            last_N = max(self.population_size_history[:, self.num_epochs - 1])
            last_epoch = max(self.epoch_times)
            steps = sorted(
                list(
                    set(np.linspace(0, last_epoch + 12 * last_N, 101)).union(
                        set(self.epoch_times)
                    )
                )
            )
        p_diff = m_diff = np.inf
        last_P = np.inf
        step_type = "none"
        n = 0
        logger.debug(
            "iter    mean    P_diff    mean_diff last_P    adjust_type"
            "num_steps  last_step"
        )
        # The factors of 20 here are probably not optimal: clearly, we need to
        # compute P accurately, but there's no good reason for this stopping rule.
        # If populations have picewise constant size then we shouldn't need this:
        # setting steps equal to the epoch boundaries should suffice; while if
        # there is very fast exponential change in some epochs caution is needed.
        while n < max_iter and (
            last_P > rtol or p_diff > rtol / 20 or m_diff > rtol / 20
        ):
            last_steps = steps
            _, P1 = self.coalescence_rate_trajectory(
                steps=last_steps,
                num_samples=num_samples,
                min_pop_size=min_pop_size,
                double_step_validation=False,
            )
            m1 = mean_time(last_steps, P1)
            if last_P > rtol:
                step_type = "extend"
                steps = np.concatenate(
                    [steps, np.linspace(steps[-1], steps[-1] * 1.2, 20)[1:]]
                )
            else:
                step_type = "refine"
                inter = steps[:-1] + np.diff(steps) / 2
                steps = np.concatenate([steps, inter])
                steps.sort()
            _, P2 = self.coalescence_rate_trajectory(
                steps=steps,
                num_samples=num_samples,
                min_pop_size=min_pop_size,
                double_step_validation=False,
            )
            m2 = mean_time(steps, P2)
            keep_steps = np.in1d(steps, last_steps)
            p_diff = max(np.abs(P1 - P2[keep_steps]))
            m_diff = np.abs(m1 - m2) / m2
            last_P = P2[-1]
            n += 1
            # Use the old-style string formatting as this is the logging default
            logger.debug(
                "%d %g %g %g %g %s %d %d",
                n,
                m2,
                p_diff,
                m_diff,
                last_P,
                step_type,
                len(steps),
                max(steps),
            )

        if n == max_iter:
            raise ValueError(
                "Did not converge on an adequate discretisation: "
                "Increase max_iter or rtol. Consult the log for "
                "debugging information"
            )
        return m2

    def coalescence_rate_trajectory(
        self, steps, num_samples, min_pop_size=1, double_step_validation=True
    ):
        """
        This function will calculate the mean coalescence rates and proportions
        of uncoalesced lineages between the lineages of the sample
        configuration provided in `num_samples`, at each of the times ago
        listed by steps, in this demographic model. The coalescence rate at
        time t in the past is the average rate of coalescence of
        as-yet-uncoalesed lineages, computed as follows: let :math:`p(t)` be
        the probability that the lineages of a randomly chosen pair of samples
        has not yet coalesced by time :math:`t`, let :math:`p(z,t)` be the
        probability that the lineages of a randomly chosen pair of samples has
        not yet coalesced by time :math:`t` *and* are both in population
        :math:`z`, and let :math:`N(z,t)` be the diploid effective population
        size of population :math:`z` at time :math:`t`. Then the mean
        coalescence rate at time :math:`t` is :math:`r(t) = (\\sum_z p(z,t) /
        (2 * N(z,t)) / p(t)`.

        The computation is done by approximating population size trajectories
        with piecewise constant trajectories between each of the steps. For
        this to be accurate, the distance between the steps must be small
        enough so that (a) short epochs (e.g., bottlenecks) are not missed, and
        (b) populations do not change in size too much over that time, if they
        are growing or shrinking. This function optionally provides a simple
        check of this approximation by recomputing the coalescence rates on a
        grid of steps twice as fine and throwing a warning if the resulting
        values do not match to a relative tolerance of 0.001.

        :param list steps: The times ago at which coalescence rates will be computed.
        :param list num_samples: A list of the same length as the number
            of populations, so that `num_samples[j]` is the number of sampled
            chromosomes in subpopulation `j`.
        :param int min_pop_size: The smallest allowed population size during
            computation of coalescent rates (i.e., coalescence rates are actually
            1 / (2 * max(min_pop_size, N(z,t))). Spurious very small population sizes
            can occur in models where populations grow exponentially but are unused
            before some time in the past, and lead to floating point error.
            This should be set to a value smaller than the smallest
            desired population size in the model.
        :param bool double_step_validation: Whether to perform the check that
            step sizes are sufficiently small, as described above. This is highly
            recommended, and will take at most four times the computation.
        :return: A tuple of arrays whose jth elements, respectively, are the
            coalescence rate at the jth time point (denoted r(t[j]) above),
            and the probablility that a randomly chosen pair of lineages has
            not yet coalesced (denoted p(t[j]) above).
        :rtype: (numpy.array, numpy.array)
        """
        num_pops = self.num_populations
        if not len(num_samples) == num_pops:
            raise ValueError(
                "`num_samples` must have the same length as the number of populations"
            )
        steps = np.array(steps)
        if not np.all(np.diff(steps) > 0):
            raise ValueError("`steps` must be a sequence of increasing times.")
        if np.any(steps < 0):
            raise ValueError("`steps` must be non-negative")
        r, p_t = self._calculate_coalescence_rate_trajectory(
            steps=steps, num_samples=num_samples, min_pop_size=min_pop_size
        )
        if double_step_validation:
            inter = steps[:-1] + np.diff(steps) / 2
            double_steps = np.concatenate([steps, inter])
            double_steps.sort()
            rd, p_td = self._calculate_coalescence_rate_trajectory(
                steps=double_steps, num_samples=num_samples, min_pop_size=min_pop_size
            )
            assert np.all(steps == double_steps[::2])
            r_prediction_close = np.allclose(r, rd[::2], rtol=1e-3, equal_nan=True)
            p_prediction_close = np.allclose(p_t, p_td[::2], rtol=1e-3, equal_nan=True)
            if not (r_prediction_close and p_prediction_close):
                warnings.warn(
                    "Doubling the number of steps has resulted in different "
                    " predictions, please re-run with smaller step sizes to ensure "
                    " numerical accuracy."
                )
        return r, p_t

    def _calculate_coalescence_rate_trajectory(self, steps, num_samples, min_pop_size):
        num_pops = self.num_populations
        P = np.zeros([num_pops ** 2, num_pops ** 2])
        IA = np.array(range(num_pops ** 2)).reshape([num_pops, num_pops])
        Identity = np.eye(num_pops)
        for x in range(num_pops):
            for y in range(num_pops):
                P[IA[x, y], IA[x, y]] = num_samples[x] * (num_samples[y] - (x == y))
        P = P / np.sum(P)
        # add epoch breaks if not there already but remember which steps they are
        epoch_breaks = list(
            set([0.0] + [t for t in self.epoch_times if t not in steps])
        )
        steps_b = np.concatenate([steps, epoch_breaks])
        ix = np.argsort(steps_b)
        steps_b = steps_b[ix]
        keep_steps = np.concatenate(
            [np.repeat(True, len(steps)), np.repeat(False, len(epoch_breaks))]
        )[ix]
        assert np.all(steps == steps_b[keep_steps])
        mass_migration_objects = []
        mass_migration_times = []
        for demo in self.demography.events:
            if type(demo) == MassMigration:
                mass_migration_objects.append(demo)
                mass_migration_times.append(demo.time)
        num_steps = len(steps_b)
        # recall that steps_b[0] = 0.0
        r = np.zeros(num_steps)
        p_t = np.zeros(num_steps)
        for j in range(num_steps - 1):
            time = steps_b[j]
            dt = steps_b[j + 1] - steps_b[j]
            N, M = self._pop_size_and_migration_at_t(time)
            C = np.zeros([num_pops ** 2, num_pops ** 2])
            for idx in range(num_pops):
                C[IA[idx, idx], IA[idx, idx]] = 1 / (2 * max(min_pop_size, N[idx]))
            dM = np.diag([sum(s) for s in M])
            if time in mass_migration_times:
                idx = mass_migration_times.index(time)
                a = mass_migration_objects[idx].source
                b = mass_migration_objects[idx].dest
                p = mass_migration_objects[idx].proportion
                S = np.eye(num_pops ** 2, num_pops ** 2)
                for x in range(num_pops):
                    if x == a:
                        S[IA[a, a], IA[a, b]] = S[IA[a, a], IA[b, a]] = p * (1 - p)
                        S[IA[a, a], IA[b, b]] = p ** 2
                        S[IA[a, a], IA[a, a]] = (1 - p) ** 2
                    else:
                        S[IA[x, a], IA[x, b]] = S[IA[a, x], IA[b, x]] = p
                        S[IA[x, a], IA[x, a]] = S[IA[a, x], IA[a, x]] = 1 - p
                P = np.matmul(P, S)
            p_notcoal = np.sum(P)
            p_t[j] = p_notcoal
            if p_notcoal > 0:
                r[j] = np.sum(np.matmul(P, C)) / np.sum(P)
            else:
                r[j] = np.nan
            G = (np.kron(M - dM, Identity) + np.kron(Identity, M - dM)) - C
            P = np.matmul(P, _matrix_exponential(dt * G))
        p_notcoal = np.sum(P)
        p_t[num_steps - 1] = p_notcoal
        if p_notcoal > 0:
            r[num_steps - 1] = np.sum(np.matmul(P, C)) / p_notcoal
        else:
            r[num_steps - 1] = np.nan
        return r[keep_steps], p_t[keep_steps]

    def _pop_size_and_migration_at_t(self, t):
        """
        Returns a tuple (N, M) of population sizes (N) and migration rates (M) at
        time t ago.

        Note: this isn't part of the external API as it is be better to provide
        separate methods to access the population size and migration rates, and
        needing both together is specialised for internal calculations.

        :param float t: The time ago.
        :return: A tuple of arrays, of the same form as the population sizes and
            migration rate arrays of the demographic model.
        """
        j = 0
        while self.epochs[j].end_time <= t:
            j += 1
        N = self.population_size_history[:, j]
        for i, pop in enumerate(self.epochs[j].populations):
            s = t - self.epochs[j].start_time
            g = pop.growth_rate
            N[i] *= np.exp(-1 * g * s)
        return N, self.epochs[j].migration_matrix

    @property
    def population_size_history(self):
        """
        Returns a (num_pops, num_epochs) numpy array giving the starting population size
        for each population in each epoch.
        """
        num_pops = len(self.epochs[0].populations)
        pop_size = np.zeros((num_pops, len(self.epochs)))
        for j, epoch in enumerate(self.epochs):
            for k, pop in enumerate(epoch.populations):
                pop_size[k, j] = pop.start_size
        return pop_size

    @property
    def epoch_times(self):
        """
        Returns array of epoch times defined by the demographic model
        """
        return np.array([x.start_time for x in self.epochs])

    @property
    def num_epochs(self):
        """
        Returns the number of epochs defined by the demographic model.
        """
        return len(self.epochs)
