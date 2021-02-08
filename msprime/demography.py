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
import dataclasses
import inspect
import logging
import math
import sys
import textwrap
import warnings
from typing import ClassVar
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


def check_population_size(Ne):
    """
    Check if an input population size makes sense.
    """
    if Ne is not None and Ne <= 0:
        raise ValueError("Population size must be positive")


@dataclasses.dataclass
class Demography:
    """
    A description of a demographic model for an msprime simulation.

    TODO document properly.

    Population structure is modelled by specifying a fixed number of
    subpopulations :math:`d`, and a :math:`d \\times d` matrix :math:`M` of
    per-generation migration rates. The :math:`(j,k)^{th}` entry of :math:`M`
    is the expected number of migrants moving from population :math:`k` to
    population :math:`j` per generation, divided by the size of population
    :math:`j`. In terms of the coalescent process, :math:`M_{j,k}` gives the
    rate at which an ancestral lineage moves from population :math:`j` to
    population :math:`k`, as one follows it back through time. In
    continuous-time models, when :math:`M_{j,k}` is close to zero, this rate is
    approximately equivalent to the fraction of population :math:`j` that is
    replaced each generation by migrants from population :math:`k`. In
    discrete-time models, the equivalence is exact and each row of :math:`M`
    has the constraint :math:`\\sum_{k \\neq j} M_{j,k} \\leq 1`. This differs
    from the migration matrix one usually uses in population demography: if
    :math:`m_{k,j}` is the proportion of individuals (in the usual sense; not
    lineages) in population :math:`k` that move to population :math:`j` per
    generation, then translating this proportion of population :math:`k` to a
    proportion of population :math:`j`, we have :math:`M_{j,k} = m_{k,j}
    \\times N_k / N_j`.

    Each subpopulation has an initial absolute population size :math:`s`
    and a per generation exponential growth rate :math:`\\alpha`. The size of a
    given population at time :math:`t` in the past (measured in generations) is
    therefore given by :math:`s e^{-\\alpha t}`. Demographic events that occur in
    the history of the simulated population alter some aspect of this population
    configuration at a particular time in the past.

    """

    populations: List[Population] = dataclasses.field(default_factory=list)
    events: List = dataclasses.field(default_factory=list)
    migration_matrix: Union[np.ndarray, None] = None

    def __post_init__(self):
        if self.migration_matrix is None:
            N = self.num_populations
            self.migration_matrix = np.zeros((N, N))

        # Sort demographic events by time.
        self.events.sort(key=lambda de: de.time)
        self.__name_id_map = None

    def _populations_table(self):
        col_titles = [
            "id",
            "name",
            "description",
            "initial_size",
            "growth_rate",
            "extra_metadata",
        ]
        data = []
        for j, pop in enumerate(self.populations):
            row = [str(j)] + [f"{getattr(pop, attr)}" for attr in col_titles[1:]]
            data.append(row)
        return col_titles, data

    def _populations_text(self):
        col_titles, data = self._populations_table()
        alignments = ["^", "<", "<", "<", "^", "<"]
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
            extra = (
                "Backwards in time migration rate from population "
                f"{self.populations[source].name} to {self.populations[dest].name} "
                f"= {rate} per generation. "
                "Equivalant to **IMPLEMENT ME** forwards in time"
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
            description = textwrap.wrap(event._parameters(), 20)
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
        return (
            "<p>"
            + self._populations_html()
            + self._migration_matrix_html()
            + self._events_html(self.events)
            + "</p>"
        )

    def __str__(self):
        populations = self._populations_text()
        migration_matrix = self._migration_matrix_text()
        events = self._events_text(self.events)

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

    def name_to_id(self, name):
        """
        Returns the integer ID (i.e., its position in the list of populations)
        of the population with the specified name. If the name does not exist,
        raise a KeyError.

        Note: this function will raise an error if called before the ``validate``
        method is called.

        :param str name: The name of the population we wish to look up.
        :return: The integer ID of the population.
        :rtype: int
        """
        if self.__name_id_map is None:
            raise ValueError("Cannot call name_to_id before calling validate()")
        if name not in self.__name_id_map:
            raise KeyError(f"Population with name '{name}' not found in demography")
        return self.__name_id_map[name]

    @property
    def num_populations(self):
        return len(self.populations)

    @property
    def num_events(self):
        return len(self.events)

    def validate(self):
        """
        Checks the demography looks sensible and raises errors/warnings
        appropriately.
        """
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

        for event in self.events:
            if not isinstance(event, DemographicEvent):
                raise TypeError(
                    "Demographic events must be a list of DemographicEvent "
                    "instances sorted in non-decreasing order of time."
                )

        self.__name_id_map = {}
        for j, population in enumerate(self.populations):
            population.validate()
            if population.name in self.__name_id_map:
                raise ValueError(f"Duplicate population name: '{population.name}'")
            self.__name_id_map[population.name] = j

    def insert_populations(self, tables):
        """
        Insert population definitions for this demography into the specified
        set of tables.
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
        return dataclasses.asdict(self)

    def debug(self):
        """
        Returns a :class:`.DemographyDebugger` instance for this demography.

        :return: A DemographyDebugger object for this demography.
        :rtype: .DemographyDebugger
        """
        return DemographyDebugger(demography=self)

    def __eq__(self, other):
        if isinstance(other, Demography):
            return (
                self.populations == other.populations
                and np.array_equal(self.migration_matrix, other.migration_matrix)
                and self.events == other.events
            )
        else:
            return super().__eq__(other)

    @staticmethod
    def from_species_tree(
        tree,
        initial_size,
        *,
        time_units="gen",
        generation_time=None,
        growth_rate=None,
    ):
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
        internal node, :class:`.MassMigration` events are added so that
        lineages move from its child populations at the appropriate time.

        :warning: If continuous migration is added to the returned Demography
            by updating the migration matrix, it is important to
            note that the MassMigration events used to move lineages do *not*
            alter migration rates, and it should be ensured that
            migration rates to source populations of mass migration events are zero
            after the mass migration (viewed backwards in time).

        :todo: Implement the PopulationSplit event and document its use here.
            We'll still need to put in a warning, so users know that the migration
            matrix will be updated at every node.

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
        ancestral populations using a `collections.defaultdict`, e.g.,

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
    def from_starbeast(tree, generation_time, time_units="myr"):
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

        The returned :class:`.Demography` object contains a :class:`.Population` for
        each node in the species tree. The population's ``name`` attribute will
        be either the corresponding node label from the newick tree, if it exists,
        or otherwise the name takes the form "pop_{j}", where j is the position
        of the given population in the list. Leaf populations are first in the
        list, and added in left-to-right order. Populations corresponding to the
        internal nodes are then added in a postorder traversal of the species
        tree. For each internal node, :class:`.MassMigration` events are added
        so that lineages move from its child populations at the appropriate time.

        :warning: If continuous migration is added to the returned Demography
            by updating the migration matrix, it is important to
            note that the MassMigration events used to move lineages do *not*
            alter migration rates, and it should be ensured that
            migration rates to source populations of mass migration events are zero
            after the mass migration (viewed backwards in time).

        :todo: Implement the PopulationSplit event and document its use here.
            We'll still need to put in a warning, so users know that the migration
            matrix will be updated at every node.

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

    @staticmethod
    def from_old_style(
        population_configurations=None,
        migration_matrix=None,
        demographic_events=None,
        Ne=1,
    ):
        """
        Creates a Demography object from the pre 1.0 style input parameters,
        reproducing the old semantics with respect to default values.
        """
        demography = Demography()
        if population_configurations is None:
            demography.populations = [Population(initial_size=Ne)]
        else:
            for pop_config in population_configurations:
                demography.populations.append(Population.from_old_style(pop_config, Ne))
        for j, population in enumerate(demography.populations):
            if population.name is None:
                population.name = f"pop_{j}"
        if migration_matrix is None:
            migration_matrix = np.zeros(
                (demography.num_populations, demography.num_populations)
            )
        demography.migration_matrix = migration_matrix
        if demographic_events is not None:
            demography.events = demographic_events
        return demography

    @staticmethod
    def isolated_model(initial_size, *, growth_rate=None):
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
                name=f"pop_{j}",
            )
            for j in range(len(initial_size))
        ]
        return Demography(populations=populations)

    @staticmethod
    def island_model(initial_size, migration_rate, *, growth_rate=None):
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
    ):
        """
        Returns a :class:`.Demography` object representing a collection of
        populations with specified initial population sizes and growth
        rates, in which adjacent demes exchange migrants at the
        specified rate. Please see :ref:`sec_demography` for more details on
        population sizes and growth rates.

        .. note:: The current implementation on supports a one-dimensional stepping
            stone model, but higher dimensions could also be supported. Please
            open an issue on GitHub if this feature would be useful to you.

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

    # TODO commenting these two out for now in the interest of getting the
    # main changes merged. We can update the code to include these models
    # ported in from stdpopsim as a follow up.

    # @staticmethod
    # def piecewise_constant_size(N0, *args):
    #     """
    #     Returns a piecewise constant size demographic model, which allows for
    #     instantaneous population size change over multiple epochs in a single
    #     population.

    #     :ivar N0: The initial effective population size
    #     :vartype N0: float
    #     :ivar args: Each subsequent argument is a tuple (t, N) which gives the
    #         time at which the size change takes place and the population size.

    #     The usage is best illustrated by an example:

    #     .. code-block:: python

    #         # One change
    #         model1 = msprime.Demography.piecewise_constant_size(N0, (t1, N1))
    #         # Two changes
    #         model2 = msprime.Demography.piecewise_constant_size(
    #             N0, (t1, N1), (t2, N2))
    #     """
    #     model = Demography()
    #     model.populations = [
    #         Population(
    #             initial_size=N0,
    #             name="pop_0",
    #             description="Population in piecewise constant model",
    #         )
    #     ]
    #     model.migration_matrix = [[0]]
    #     model.demographic_events = []
    #     for t, N in args:
    #         model.demographic_events.append(
    #             PopulationParametersChange(
    #                 time=t, initial_size=N, growth_rate=0, population=0
    #             )
    #         )
    #     return model

    # @staticmethod
    # def im_model(N0, N1, NA, T, M01, M10):
    #     """
    #     An isolation with migration model where a single ancestral
    #     population of size NA splits into two populations of constant size N0
    #     and N1 at time T generations ago, with migration rates M01 and M10 between
    #     the split populations. Sampling is disallowed in population index 2,
    #     as this is the ancestral population.

    #     .. fixme:: Not clear what direction the migration rates are here.

    #     :ivar N0: The effective population size of population 0
    #     :vartype N0: float
    #     :ivar N1: The effective population size of population 1
    #     :vartype N1: float
    #     :ivar NA: The ancestral effective population size (population 2).
    #     :vartype NA: float
    #     :ivar T: Time of split between populations 0 and 1 (in generations)
    #     :vartype T: float
    #     :ivar M01: Migration rate from population 0 to 1
    #     :vartype M01: float
    #     :ivar M10: Migration rate from population M10
    #     :vartype M10: float

    #     Example usage:

    #     .. code-block:: python

    #         model = msprime.Demography.im_model(N0, N1, NA, T, M12, M21)
    #     """
    #     model = Demography()
    #     model.populations = [
    #         Population(initial_size=N0),
    #         Population(initial_size=N1),
    #         Population(initial_size=NA),
    #     ]
    #     model.migration_matrix = [[0, M01, 0], [M10, 0, 0], [0, 0, 0]]
    #     model.demographic_events = [
    #         MassMigration(time=T, source=0, destination=2, proportion=1),
    #         MassMigration(time=T, source=1, destination=2, proportion=1),
    #     ]
    #     return model


# TODO fixup the documentation here. This definition is partly lifted
# from stdpopsim and also derived from PopulationConfiguration. The
# idea is that PopulationConfiguration is to be maintained
# indefinitely as the old-style interface, but we convert to Population
# internally. Getting rid of PopConfig makes sense for two reasons:
# 1) it's a clunky and confusing name and 2) The sample_size
# attribute makes things really awkward, as it conflates declaring
# structure and how you sample from it.


@dataclasses.dataclass
class Population:
    """
    Define a single population in a simulation.

    :ivar initial_size: The absolute size of the population at time zero.
    :vartype initial_size: float
    :var growth_rate: The exponential growth rate of the
        population per generation (forwards in time).
        Growth rates can be negative. This is zero for a
        constant population size, and positive for a population that has been
        growing. Defaults to 0.
    :vartype growth_rate: float
    :ivar name: The name of the population. If specified this must be a uniquely
        identifying string.
    :vartype name: str
    :ivar description: a short description of the population
    :vartype description: str
    """

    initial_size: float
    growth_rate: float = 0.0
    name: Union[str, None] = None
    description: Union[str, None] = None
    extra_metadata: dict = dataclasses.field(default_factory=dict)
    sampling_time: float = 0

    def asdict(self):
        return dataclasses.asdict(self)

    @staticmethod
    def from_old_style(pop_config, Ne=1):
        """
        Returns a Population object derived from the specified old-style
        PopulationConfiguration. The Ne value is used as the ``initial_size``
        if this is not provided in the PopulationConfiguration.
        """
        initial_size = (
            Ne if pop_config.initial_size is None else pop_config.initial_size
        )
        population = Population(
            initial_size=initial_size,
            growth_rate=pop_config.growth_rate,
        )
        metadata = pop_config.metadata
        if metadata is not None and isinstance(metadata, collections.abc.Mapping):
            metadata = metadata.copy()
            if "name" in metadata:
                population.name = metadata.pop("name")
            if "description" in metadata:
                population.description = metadata.pop("description")
        population.extra_metadata = metadata
        return population

    def validate(self):
        # TODO more checks
        if self.initial_size < 0:
            raise ValueError("Negative population size")
        if self.name is None:
            raise ValueError("A population name must be set.")
        if not self.name.isidentifier():
            raise ValueError("A population name must be a valid Python identifier")


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


@dataclasses.dataclass
class DemographicEvent:
    """
    Superclass of demographic events that occur during simulations.
    """

    time: float

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


@dataclasses.dataclass
class PopulationParametersChange(DemographicEvent):
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
            "population": self.population,
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
            s += f"initial_size → {self.initial_size} "
            if self.growth_rate is not None:
                s += "and "
        if self.growth_rate is not None:
            s += f"growth_rate → {self.growth_rate} "
        s += "for"
        if self.population == -1:
            s += " all populations"
        else:
            s += f" population {self.population}"
        return s


@dataclasses.dataclass
class MigrationRateChange(DemographicEvent):
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
            "migration_rate": self.rate,
            "source": self.source,
            "dest": self.dest,
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


@dataclasses.dataclass
class MassMigration(DemographicEvent):
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
            "source": self.source,
            "dest": self.dest,
            "proportion": self.proportion,
        }

    def _parameters(self):
        return f"source={self.source}, dest={self.dest}, proportion={self.proportion}"

    def _effect(self):
        return (
            f"Lineages currently in population {self.source} move to {self.dest} "
            f"with probability {self.proportion} (equivalent to individuals "
            f"migrating from {self.dest} to {self.source} forwards in time)"
        )


# This is an unsupported/undocumented demographic event.
@dataclasses.dataclass
class SimpleBottleneck(DemographicEvent):
    population: int
    proportion: float = 1.0

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "simple_bottleneck",
            "time": self.time,
            "population": self.population,
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
class InstantaneousBottleneck(DemographicEvent):
    population: int
    strength: float = 1.0

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "instantaneous_bottleneck",
            "time": self.time,
            "population": self.population,
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


@dataclasses.dataclass
class PopulationState:
    """
    Simple class to represent the state of a population in terms of its
    demographic parameters.
    """

    start_size: float
    end_size: float
    growth_rate: float


@dataclasses.dataclass
class Epoch:
    """
    Represents a single epoch in the simulation within which the state
    of the demographic parameters are constant.
    """

    start_time: float
    end_time: float
    populations: List[PopulationState]
    migration_matrix: list  # TODO numpy array
    demographic_events: List[DemographicEvent]


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


class DemographyDebugger:
    """
    A class to facilitate debugging of population parameters and migration
    rates in the past.
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
                migration_matrix,
                demographic_events,
                Ne=Ne,
            )
        self.demography = demography
        self.num_populations = demography.num_populations
        self._make_epochs()
        self._check_misspecification()

    def _make_epochs(self):
        self.epochs = []
        # We don't actually use the samples, this is just to get the simulator
        # correctly initialised.
        samples = {}
        for j, pop in enumerate(self.demography.populations):
            if pop.initial_size > 0 and pop.sampling_time == 0:
                samples[j] = 1
        if len(samples) == 0:
            raise ValueError("No population with non-zero initial size.")
        simulator = ancestry._parse_sim_ancestry(
            demography=self.demography, samples=samples
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
            growth_rates = [
                conf["growth_rate"] for conf in simulator.population_configuration
            ]
            populations = [
                PopulationState(
                    start_size=simulator.compute_population_size(j, start_time),
                    end_size=simulator.compute_population_size(j, end_time),
                    growth_rate=growth_rates[j],
                )
                for j in range(self.num_populations)
            ]
            self.epochs.append(
                Epoch(start_time, end_time, populations, migration_matrix, events)
            )
            start_time = end_time

    def _check_misspecification(self):
        """
        Check for things that might indicate model misspecification.
        """
        merged_pops = set()
        for epoch in self.epochs:
            for de in epoch.demographic_events:
                if isinstance(de, MassMigration) and de.proportion == 1:
                    merged_pops.add(de.source)
            mm = epoch.migration_matrix
            for k in merged_pops:
                if any(mm[k, :] != 0) or any(mm[:, k] != 0):
                    warnings.warn(
                        "Non-zero migration rates exist after merging "
                        f"population {k}. This almost certainly indicates "
                        "demographic misspecification."
                    )

    def _populations_html(self, epoch):
        column_titles = ["", "start", "end", "growth_rate"] + [
            pop.name for pop in self.demography.populations
        ]
        data = []
        for j, pop in enumerate(epoch.populations):
            row = [
                self.demography.populations[j].name,
                f"{pop.start_size: .3g}",
                f"{pop.end_size: .3g}",
                f"{pop.growth_rate: .3g}",
            ]
            for k in range(self.demography.num_populations):
                item = self.demography._migration_rate_info(
                    j, k, epoch.migration_matrix[j, k]
                )
                row.append(item)
            data.append(row)
        return core.html_table("", column_titles, data)

    def _populations_text(self, epoch):
        column_titles = [[""], ["start"], ["end"], ["growth_rate"]] + [
            [pop.name] for pop in self.demography.populations
        ]
        alignments = ">>><" + "^" * self.demography.num_populations
        data = []
        for j, pop in enumerate(epoch.populations):
            row = [
                [self.demography.populations[j].name],
                [f"{pop.start_size: .3g}"],
                [f"{pop.end_size: .3g}"],
                [f"{pop.growth_rate: .3g}"],
            ]
            for k in range(self.demography.num_populations):
                row.append([f"{epoch.migration_matrix[j, k]:.3g}"])
            data.append(row)
        return core.text_table("Population state", column_titles, alignments, data)

    def _repr_html_(self):
        out = ""
        for j, epoch in enumerate(self.epochs):
            if j > 0:
                if len(epoch.demographic_events) > 0:
                    title = f"Events @ generation {epoch.start_time}"
                    out += self.demography._events_html(epoch.demographic_events, title)
                out += "</p>"
            else:
                assert len(epoch.demographic_events) == 0
            epoch_title = f"Epoch: {epoch.start_time} -- {epoch.end_time} generations"
            out += "<p>"
            out += f"<h3>{epoch_title}</h3>"
            out += self._populations_html(epoch)
        out += "</p>"
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
        for j, epoch in enumerate(self.epochs):
            if j > 0:
                if len(epoch.demographic_events) > 0:
                    title = f"Events @ generation {epoch.start_time}"
                    out += indent(
                        self.demography._events_text(epoch.demographic_events, title)
                    )
            epoch_title = f"Epoch: {epoch.start_time} -- {epoch.end_time} generations"
            out += box(epoch_title)
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
