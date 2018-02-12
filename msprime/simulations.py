#
# Copyright (C) 2015-2017 University of Oxford
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
Module responsible for running simulations.
"""
from __future__ import division
from __future__ import print_function

import collections
import gzip
import json
import math
import random
import sys

import _msprime
import msprime.tables as tables
import msprime.trees as trees
import msprime.provenance as provenance

# Make the low-level generator appear like its from this module
from _msprime import RandomGenerator
from _msprime import MutationGenerator


Sample = collections.namedtuple(
    "Sample",
    ["population", "time"])


def _get_random_seed():
    return random.randint(1, 2**32 - 1)


def almost_equal(a, b, rel_tol=1e-9, abs_tol=0.0):
    """
    Returns true if the specified pair of numbers are equal to
    within the specified tolerances.

    The signature and implementation are taken from PEP 485,
    https://www.python.org/dev/peps/pep-0485/
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def harmonic_number(n):
    """
    Returns the nth Harmonic number.
    """
    EulerGamma = 0.5772156649
    return math.log(n) + EulerGamma


def _check_population_configurations(population_configurations):
    err = (
        "Population configurations must be a list of PopulationConfiguration instances")
    if not isinstance(population_configurations, collections.Iterable):
        raise TypeError(err)
    for config in population_configurations:
        if not isinstance(config, PopulationConfiguration):
            raise TypeError(err)


def _replicate_generator(sim, mutation_generator, num_replicates, provenance_dict):
    """
    Generator function for the many-replicates case of the simulate
    function.
    """
    # TODO We should encode the replicate index in here with the rest of the
    # parameters. This will provide sufficient information to reproduce the
    # simulation if necessary. Much simpler than encoding the details of
    # the random number generator.
    provenance_record = json.dumps(provenance_dict)

    # Should use range here, but Python 2 makes this awkward...
    j = 0
    while j < num_replicates:
        j += 1
        sim.run()
        tree_sequence = sim.get_tree_sequence(mutation_generator, provenance_record)
        yield tree_sequence
        sim.reset()


def simulator_factory(
        sample_size=None,
        Ne=1,
        random_generator=None,
        length=None,
        recombination_rate=None,
        recombination_map=None,
        population_configurations=None,
        migration_matrix=None,
        samples=None,
        demographic_events=[],
        model=None,
        record_migrations=False):
    """
    Convenience method to create a simulator instance using the same
    parameters as the `simulate` function. Primarily used for testing.
    """
    condition = (
        sample_size is None and
        population_configurations is None and
        samples is None)
    if condition:
        raise ValueError(
            "Either sample_size, population_configurations or samples "
            "be specified")
    if sample_size is not None:
        if samples is not None:
            raise ValueError(
                "Cannot specify sample size and samples simultaneously.")
        if population_configurations is not None:
            raise ValueError(
                "Cannot specify sample size and population_configurations "
                "simultaneously.")
        s = Sample(population=0, time=0.0)
        the_samples = [s for _ in range(sample_size)]
    # If we have population configurations we may have embedded sample_size
    # values telling us how many samples to take from each population.
    if population_configurations is not None:
        _check_population_configurations(population_configurations)
        if samples is None:
            the_samples = []
            for j, conf in enumerate(population_configurations):
                if conf.sample_size is not None:
                    the_samples += [(j, 0) for _ in range(conf.sample_size)]
        else:
            for conf in population_configurations:
                if conf.sample_size is not None:
                    raise ValueError(
                        "Cannot specify population configuration sample size"
                        "and samples simultaneously")
            the_samples = samples
    elif samples is not None:
        the_samples = samples

    if recombination_map is None:
        the_length = 1 if length is None else length
        the_rate = 0 if recombination_rate is None else recombination_rate
        if the_length <= 0:
            raise ValueError("Cannot provide non-positive sequence length")
        if the_rate < 0:
            raise ValueError("Cannot provide negative recombination rate")
        recomb_map = RecombinationMap.uniform_map(the_length, the_rate)
    else:
        if length is not None or recombination_rate is not None:
            raise ValueError(
                "Cannot specify length/recombination_rate along with "
                "a recombination map")
        recomb_map = recombination_map

    sim = Simulator(the_samples, recomb_map, model, Ne)
    sim.store_migrations = record_migrations
    rng = random_generator
    if rng is None:
        rng = RandomGenerator(_get_random_seed())
    sim.random_generator = rng
    if population_configurations is not None:
        sim.set_population_configurations(population_configurations)
    if migration_matrix is not None:
        sim.set_migration_matrix(migration_matrix)
    if demographic_events is not None:
        sim.set_demographic_events(demographic_events)
    return sim


def simulate(
        sample_size=None,
        Ne=1,
        length=None,
        recombination_rate=None,
        recombination_map=None,
        mutation_rate=None,
        population_configurations=None,
        migration_matrix=None,
        demographic_events=[],
        samples=None,
        model=None,
        record_migrations=False,
        random_seed=None,
        mutation_generator=None,
        num_replicates=None):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns the resulting :class:`.TreeSequence`.

    :param int sample_size: The number of individuals in our sample.
        If not specified or None, this defaults to the sum of the
        subpopulation sample sizes. Either ``sample_size``,
        ``population_configurations`` or ``samples`` must be specified.
    :param float Ne: The effective (diploid) population size for the reference
        population. This determines the factor by which the per-generation
        recombination and mutation rates are scaled in the simulation.
        This defaults to 1 if not specified.
    :param float length: The length of the simulated region in bases.
        This parameter cannot be used along with ``recombination_map``.
        Defaults to 1 if not specified.
    :param float recombination_rate: The rate of recombination per base
        per generation. This parameter cannot be used along with
        ``recombination_map``. Defaults to 0 if not specified.
    :param recombination_map: The map
        describing the changing rates of recombination along the simulated
        chromosome. This parameter cannot be used along with the
        ``recombination_rate`` or ``length`` parameters, as these
        values are encoded within the map. Defaults to a uniform rate as
        described in the ``recombination_rate`` parameter if not specified.
    :type recombination_map: :class:`.RecombinationMap`
    :param float mutation_rate: The rate of mutation per base per
        generation. If not specified, no mutations are generated.
    :param list population_configurations: The list of
        :class:`.PopulationConfiguration` instances describing the
        sampling configuration, relative sizes and growth rates of
        the populations to be simulated. If this is not specified,
        a single population with a sample of size ``sample_size``
        is assumed.
    :type population_configurations: list or None.
    :param list migration_matrix: The matrix describing the rates
        of migration between all pairs of populations. If :math:`N`
        populations are defined in the ``population_configurations``
        parameter, then the migration matrix must be an
        :math:`N\\times N` matrix consisting of :math:`N` lists of
        length :math:`N` or an :math`N\\times N` numpy array.
    :param list demographic_events: The list of demographic events to
        simulate. Demographic events describe changes to the populations
        in the past. Events should be supplied in non-decreasing
        order of time. Events with the same time value will be applied
        sequentially in the order that they were supplied before the
        simulation algorithm continues with the next time step.
    :param list samples: The list specifying the location and time of
        all samples. This parameter may be used to specify historical
        samples, and cannot be used in conjunction with the ``sample_size``
        parameter. Each sample is a (``population``, ``time``) pair
        such that the sample in position ``j`` in the list of samples
        is drawn in the specified population at the specfied time. Time
        is measured in generations, as elsewhere.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
    :param int num_replicates: The number of replicates of the specified
        parameters to simulate. If this is not specified or None,
        no replication is performed and a :class:`.TreeSequence` object
        returned. If :obj:`num_replicates` is provided, the specified
        number of replicates is performed, and an iterator over the
        resulting :class:`.TreeSequence` objects returned.
    :return: The :class:`.TreeSequence` object representing the results
        of the simulation if no replication is performed, or an
        iterator over the independent replicates simulated if the
        :obj:`num_replicates` parameter has been used.
    :rtype: :class:`.TreeSequence` or an iterator over
        :class:`.TreeSequence` replicates.
    :warning: If using replication, do not store the results of the
        iterator in a list! For performance reasons, the same
        underlying object may be used for every TreeSequence
        returned which will most likely lead to unexpected behaviour.
    """
    seed = random_seed
    if random_seed is None:
        seed = _get_random_seed()
    # To support numpy integer inputs here too we convert to integer.
    rng = RandomGenerator(int(seed))
    sim = simulator_factory(
        sample_size=sample_size, random_generator=rng,
        Ne=Ne, length=length,
        recombination_rate=recombination_rate,
        recombination_map=recombination_map,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        samples=samples, model=model, record_migrations=record_migrations)
    # The provenance API is very tentative, and only included now as a
    # pre-alpha feature.
    parameters = {"TODO": "encode simulation parameters"}
    provenance_dict = provenance.get_provenance_dict("simulate", parameters)
    if mutation_generator is None:
        mu = 0 if mutation_rate is None else mutation_rate
        mutation_generator = MutationGenerator(rng, mu)
    else:
        if mutation_rate is not None:
            raise ValueError(
                "Cannot specify both mutation_rate and mutation_generator")
    if num_replicates is None:
        return next(_replicate_generator(sim, mutation_generator, 1, provenance_dict))
    else:
        return _replicate_generator(
            sim, mutation_generator, num_replicates, provenance_dict)


class Simulator(object):
    """
    Class to simulate trees under a variety of population models.
    """
    def __init__(self, samples, recombination_map, model="hudson", Ne=0.25):
        if len(samples) < 2:
            raise ValueError("Sample size must be >= 2")
        if len(samples) >= 2**32:
            raise ValueError("sample_size must be < 2**32")
        self.sample_size = len(samples)
        self.samples = samples
        if not isinstance(recombination_map, RecombinationMap):
            raise TypeError("RecombinationMap instance required")
        self.set_model(model, Ne)
        self.recombination_map = recombination_map
        self.random_generator = None
        self.population_configurations = [
            PopulationConfiguration(initial_size=self.model.population_size)]
        self.migration_matrix = [[0]]
        self.demographic_events = []
        self.store_migrations = False
        # Set default block sizes to 64K objects.
        # TODO does this give good performance in a range of scenarios?
        block_size = 64 * 1024
        # We always need at least n segments, so no point in making
        # allocation any smaller than this.
        self.segment_block_size = max(block_size, self.sample_size)
        self.avl_node_block_size = block_size
        self.node_mapping_block_size = block_size
        self.node_block_size = block_size
        self.edge_block_size = block_size
        self.migration_block_size = block_size
        # TODO is it useful to bring back the API to set this? Mostly
        # the amount of memory required is tiny.
        self.max_memory = sys.maxsize
        self.ll_sim = None
        self.node_table = tables.NodeTable(block_size)
        self.edge_table = tables.EdgeTable(block_size)
        self.migration_table = tables.MigrationTable(block_size)
        self.site_table = tables.SiteTable()
        self.mutation_table = tables.MutationTable(block_size)
        self.provenance_table = tables.ProvenanceTable()

    @property
    def num_loci(self):
        return self.recombination_map.get_num_loci()

    @property
    def sample_configuration(self):
        return [conf.sample_size for conf in self.population_configurations]

    @property
    def num_breakpoints(self):
        return self.ll_sim.get_num_breakpoints()

    @property
    def breakpoints(self):
        """
        Returns the recombination breakpoints translated into physical
        coordinates.
        """
        return [
            self.recombination_map.genetic_to_physical(x)
            for x in self.ll_sim.get_breakpoints()]

    @property
    def used_memory(self):
        return self.ll_sim.get_used_memory()

    @property
    def time(self):
        return self.ll_sim.get_time()

    @property
    def num_avl_node_blocks(self):
        return self.ll_sim.get_num_avl_node_blocks()

    @property
    def num_node_blocks(self):
        return self.ll_sim.get_num_node_blocks()

    @property
    def num_edge_blocks(self):
        return self.ll_sim.get_num_edge_blocks()

    @property
    def num_node_mapping_blocks(self):
        return self.ll_sim.get_num_node_mapping_blocks()

    @property
    def num_segment_blocks(self):
        return self.ll_sim.get_num_segment_blocks()

    @property
    def num_common_ancestor_events(self):
        return self.ll_sim.get_num_common_ancestor_events()

    @property
    def num_rejected_common_ancestor_events(self):
        return self.ll_sim.get_num_rejected_common_ancestor_events()

    @property
    def num_recombination_events(self):
        return self.ll_sim.get_num_recombination_events()

    @property
    def num_populations(self):
        return len(self.population_configurations)

    @property
    def num_migration_events(self):
        N = self.num_populations
        matrix = [[0 for j in range(N)] for k in range(N)]
        flat = self.ll_sim.get_num_migration_events()
        for j in range(N):
            for k in range(N):
                matrix[j][k] = flat[j * N + k]
        return matrix

    @property
    def total_num_migration_events(self):
        return sum(self.ll_sim.get_num_migration_events())

    @property
    def num_multiple_recombination_events(self):
        return self.ll_sim.get_num_multiple_recombination_events()

    def get_configuration(self):
        return json.loads(self.ll_sim.get_configuration_json())

    def set_migration_matrix(self, migration_matrix):
        err = (
            "migration matrix must be a N x N square matrix encoded "
            "as a list-of-lists, where N is the number of populations "
            "defined in the population_configurations. The diagonal "
            "elements of this matrix must be zero. For example, a "
            "valid matrix for a 3 population system is "
            "[[0, 1, 1], [1, 0, 1], [1, 1, 0]]")
        if self.population_configurations is None:
            raise ValueError(
                "Cannot specify a migration matrix without also providing a "
                "population_configurations argument.")
        N = len(self.population_configurations)
        if not isinstance(migration_matrix, list):
            try:
                migration_matrix = [list(row) for row in migration_matrix]
            except TypeError:
                raise TypeError(err)
        if len(migration_matrix) != N:
            raise ValueError(err)
        for row in migration_matrix:
            if not isinstance(row, list):
                raise TypeError(err)
            if len(row) != N:
                raise ValueError(err)
        self.migration_matrix = migration_matrix

    def set_population_configurations(self, population_configurations):
        _check_population_configurations(population_configurations)
        self.population_configurations = population_configurations
        # For any populations configurations in which the initial size is None,
        # set it to the population size.
        for pop_conf in self.population_configurations:
            if pop_conf.initial_size is None:
                pop_conf.initial_size = self.model.population_size
        # Now set the default migration matrix.
        N = len(self.population_configurations)
        self.migration_matrix = [[0 for j in range(N)] for k in range(N)]

    def set_demographic_events(self, demographic_events):
        err = (
            "Demographic events must be a list of DemographicEvent instances "
            "sorted in non-decreasing order of time.")
        if not isinstance(demographic_events, collections.Iterable):
            raise TypeError(err)
        for event in demographic_events:
            if not isinstance(event, DemographicEvent):
                raise TypeError(err)
        self.demographic_events = demographic_events

    def set_model(self, model, population_size):
        """
        Sets the simulation model to the specified value. This may be either a string
        or a SimulationModel instance. If None, the default simulation model is used
        (i.e., Hudson's algorithm).
        """
        model_map = {
            "hudson": StandardCoalescent(population_size),
            "smc": SmcApproxCoalescent(population_size),
            "smc_prime": SmcPrimeApproxCoalescent(population_size),
            "dtwf": DiscreteTimeWrightFisher(population_size)
        }
        if model is None:
            model_instance = StandardCoalescent(population_size)
        elif isinstance(model, str):
            lower_model = model.lower()
            if lower_model not in model_map:
                raise ValueError("Model '{}' unknown. Choose from {}".format(
                    model, list(model_map.keys())))
            model_instance = model_map[lower_model]
        else:
            if not isinstance(model, SimulationModel):
                raise TypeError(
                    "Simulation model must be a string or an instance of "
                    "SimulationModel")
            model_instance = model
        self.model = model_instance

    def create_ll_instance(self):
        # Now, convert the high-level values into their low-level
        # counterparts.
        ll_simulation_model = self.model.get_ll_representation()
        d = len(self.population_configurations)
        # The migration matrix must be flattened.
        ll_migration_matrix = [0 for j in range(d**2)]
        for j in range(d):
            for k in range(d):
                ll_migration_matrix[j * d + k] = self.migration_matrix[j][k]
        ll_population_configuration = [
            conf.get_ll_representation() for conf in self.population_configurations]
        ll_demographic_events = [
            event.get_ll_representation(d) for event in self.demographic_events]
        ll_recomb_rate = self.recombination_map.get_per_locus_recombination_rate()
        ll_sim = _msprime.Simulator(
            samples=self.samples,
            random_generator=self.random_generator,
            model=ll_simulation_model,
            num_loci=self.recombination_map.get_num_loci(),
            recombination_rate=ll_recomb_rate,
            migration_matrix=ll_migration_matrix,
            population_configuration=ll_population_configuration,
            demographic_events=ll_demographic_events,
            store_migrations=self.store_migrations,
            max_memory=self.max_memory,
            segment_block_size=self.segment_block_size,
            avl_node_block_size=self.avl_node_block_size,
            node_mapping_block_size=self.node_mapping_block_size,
            node_block_size=self.node_block_size,
            edge_block_size=self.edge_block_size,
            migration_block_size=self.migration_block_size)
        return ll_sim

    def run(self):
        """
        Runs the simulation until complete coalescence has occurred.
        """
        if self.random_generator is None:
            raise ValueError("A random generator instance must be set")
        if self.ll_sim is None:
            self.ll_sim = self.create_ll_instance()
        self.ll_sim.run()

    def get_tree_sequence(self, mutation_generator=None, provenance_record=None):
        """
        Returns a TreeSequence representing the state of the simulation.
        """
        ll_recomb_map = self.recombination_map.get_ll_recombination_map()
        self.ll_sim.populate_tables(
            self.node_table, self.edge_table, self.migration_table,
            recombination_map=ll_recomb_map)
        if mutation_generator is not None:
            mutation_generator.generate(
                self.node_table, self.edge_table, self.site_table, self.mutation_table)
        self.provenance_table.reset()
        if provenance_record is not None:
            self.provenance_table.add_row(provenance_record)
        ll_tree_sequence = _msprime.TreeSequence()
        ll_tree_sequence.load_tables(
            self.node_table, self.edge_table, self.migration_table,
            self.site_table, self.mutation_table, self.provenance_table)
        return trees.TreeSequence(ll_tree_sequence)

    def reset(self):
        """
        Resets the simulation so that we can perform another replicate.
        """
        if self.ll_sim is not None:
            self.ll_sim.reset()


class RecombinationMap(object):
    """
    A RecombinationMap represents the changing rates of recombination
    along a chromosome. This is defined via two lists of numbers:
    ``positions`` and ``rates``, which must be of the same length.
    Given an index j in these lists, the rate of recombination
    per base per generation is ``rates[j]`` over the interval
    ``positions[j]`` to ``positions[j + 1]``. Consequently, the first
    position must be zero, and by convention the last rate value
    is also required to be zero (although it does not used).

    :param list positions: The positions (in bases) denoting the
        distinct intervals where recombination rates change. These can
        be floating point values.
    :param list rates: The list of rates corresponding to the supplied
        ``positions``. Recombination rates are specified per base,
        per generation.
    :param int num_loci: The maximum number of non-recombining loci
        in the underlying simulation. By default this is set to
        the largest possible value, allowing the maximum resolution
        in the recombination process. However, for a finite sites
        model this can be set to smaller values.
    """
    DEFAULT_NUM_LOCI = 2**32 - 1
    """
    The default number of non-recombining loci in a RecombinationMap.
    """
    def __init__(self, positions, rates, num_loci=None):
        m = self.DEFAULT_NUM_LOCI
        if num_loci is not None:
            m = num_loci
        self._ll_recombination_map = _msprime.RecombinationMap(
            m, positions, rates)

    @classmethod
    def uniform_map(cls, length, rate, num_loci=None):
        """
        Returns a :class:`.RecombinationMap` instance in which the recombination
        rate is constant over a chromosome of the specified length. The optional
        ``num_loci`` controls the number of discrete loci in the underlying
        simulation, and is by default large enough to be effectively be
        a continuous model.

        The following map can be used to simulate a true finite locus model
        with a fixed number of loci ``m``::

            >>> recomb_map = RecombinationMap.uniform_map(m, rate, num_loci=m)

        :param float length: The length of the chromosome.
        :param float rate: The rate of recombination per unit of sequence length
            along this chromosome.
        :param int num_loci: The number of discrete loci in the underlying
            simulation. By default this is set to a large number.
        """
        return cls([0, length], [rate, 0], num_loci)

    @classmethod
    def read_hapmap(cls, filename):
        """
        Parses the specified file in HapMap format. These files must contain
        a single header line (which is ignored), and then each subsequent
        line denotes a position/rate pair. Positions are in units of bases,
        and recombination rates in centimorgans/Megabase. The first column
        in this file is ignored, as are subsequence columns after the
        Position and Rate. A sample of this format is as follows::

            Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
            chr1	55550	2.981822	0.000000
            chr1	82571	2.082414	0.080572
            chr1	88169	2.081358	0.092229
            chr1	254996	3.354927	0.439456
            chr1	564598	2.887498	1.478148

        :param str filename: The name of the file to be parsed. This may be
            in plain text or gzipped plain text.
        """
        positions = []
        rates = []
        if filename.endswith(".gz"):
            f = gzip.open(filename)
        else:
            f = open(filename)
        try:
            # Skip the header line
            f.readline()
            for j, line in enumerate(f):
                pos, rate, = map(float, line.split()[1:3])
                if j == 0:
                    if pos != 0:
                        positions.append(0)
                        rates.append(0)
                positions.append(pos)
                # Rate is expressed in centimorgans per megabase, which
                # we convert to per-base rates
                rates.append(rate * 1e-8)
            if rate != 0:
                raise ValueError(
                    "The last rate provided in the recombination map must zero")
        finally:
            f.close()
        return cls(positions, rates)

    def get_ll_recombination_map(self):
        return self._ll_recombination_map

    def physical_to_genetic(self, physical_x):
        return self._ll_recombination_map.physical_to_genetic(physical_x)

    def genetic_to_physical(self, genetic_x):
        return self._ll_recombination_map.genetic_to_physical(genetic_x)

    def get_total_recombination_rate(self):
        return self._ll_recombination_map.get_total_recombination_rate()

    def get_per_locus_recombination_rate(self):
        return self._ll_recombination_map.get_per_locus_recombination_rate()

    def get_size(self):
        return self._ll_recombination_map.get_size()

    def get_num_loci(self):
        return self._ll_recombination_map.get_num_loci()

    def get_positions(self):
        return self._ll_recombination_map.get_positions()

    def get_length(self):
        # this is a wasteful creation of a list; push the method down
        # into the low-level API.
        return self._ll_recombination_map.get_positions()[-1]

    def get_rates(self):
        return self._ll_recombination_map.get_rates()


class PopulationConfiguration(object):
    """
    The initial configuration of a population (or deme) in a simulation.

    :param int sample_size: The number of initial samples that are drawn
        from this population.
    :param float initial_size: The absolute size of the population at time
        zero. Defaults to the reference population size :math:`N_e`.
    :param float growth_rate: The exponential growth rate of the population
        per generation. Growth rates can be negative. This is zero for a
        constant population size. Defaults to 0.
    """
    def __init__(self, sample_size=None, initial_size=None, growth_rate=0.0):
        if initial_size is not None and initial_size <= 0:
            raise ValueError("Population size must be > 0")
        if sample_size is not None and sample_size < 0:
            raise ValueError("Sample size must be >= 0")
        self.sample_size = sample_size
        self.initial_size = initial_size
        self.growth_rate = growth_rate

    def get_ll_representation(self):
        """
        Returns the low-level representation of this PopulationConfiguration.
        """
        return {
            "initial_size": self.initial_size,
            "growth_rate": self.growth_rate
        }


class DemographicEvent(object):
    """
    Superclass of demographic events that occur during simulations.
    """
    def __init__(self, type_, time):
        self.type = type_
        self.time = time

    def __repr__(self):
        return repr(self.__dict__)


class PopulationParametersChange(DemographicEvent):
    """
    Changes the demographic parameters of a population at a given time.

    This event generalises the ``-eg``, ``-eG``, ``-en`` and ``-eN``
    options from ``ms``. Note that unlike ``ms`` we do not automatically
    set growth rates to zero when the population size is changed.

    :param float time: The time at which this event occurs in generations.
    :param float initial_size: The absolute size of the population
        at the beginning of the time slice starting at ``time``. If None,
        this is calculated according to the initial population size and
        growth rate over the preceding time slice.
    :param float growth_rate: The new per-generation growth rate. If None,
        the growth rate is not changed. Defaults to None.
    :param int populationd: The ID of the population affected. If
        ``population`` is None, the changes affect all populations
        simultaneously.
    """
    def __init__(
            self, time, initial_size=None, growth_rate=None, population=None,
            population_id=None):
        super(PopulationParametersChange, self).__init__(
            "population_parameters_change", time)
        if population_id is not None and population is not None:
            raise ValueError(
                "population_id and population are aliases; cannot supply both.")
        if population_id is not None:
            population = population_id
        if growth_rate is None and initial_size is None:
            raise ValueError("Must specify one or more of growth_rate and initial_size")
        if initial_size is not None and initial_size <= 0:
            raise ValueError("Cannot have a population size <= 0")
        self.time = time
        self.growth_rate = growth_rate
        self.initial_size = initial_size
        self.population = -1 if population is None else population

    def get_ll_representation(self, num_populations):
        ret = {
            "type": self.type,
            "time": self.time,
            "population": self.population
        }
        if self.growth_rate is not None:
            ret["growth_rate"] = self.growth_rate
        if self.initial_size is not None:
            ret["initial_size"] = self.initial_size
        return ret

    def __str__(self):
        s = "Population parameter change for {}: ".format(self.population)
        if self.initial_size is not None:
            s += "initial_size -> {} ".format(self.initial_size)
        if self.growth_rate is not None:
            s += "growth_rate -> {} ".format(self.growth_rate)
        return s


class MigrationRateChange(DemographicEvent):
    """
    Changes the rate of migration to a new value at a specific time.

    :param float time: The time at which this event occurs in generations.
    :param float rate: The new per-generation migration rate.
    :param tuple matrix_index: A tuple of two population IDs descibing
        the matrix index of interest. If ``matrix_index`` is None, all
        non-diagonal entries of the migration matrix are changed
        simultaneously.
    """
    def __init__(self, time, rate, matrix_index=None):
        super(MigrationRateChange, self).__init__(
            "migration_rate_change", time)
        self.rate = rate
        self.matrix_index = matrix_index

    def get_ll_representation(self, num_populations):
        matrix_index = -1
        if self.matrix_index is not None:
            matrix_index = self.matrix_index[0] * num_populations + self.matrix_index[1]
        return {
            "type": self.type,
            "time": self.time,
            "migration_rate": self.rate,
            "matrix_index": matrix_index
        }

    def __str__(self):
        if self.matrix_index is None:
            ret = "Migration rate change to {} everywhere".format(self.rate)
        else:
            ret = "Migration rate change for {} to {}".format(
                self.matrix_index, self.rate)
        return ret


class MassMigration(DemographicEvent):
    """
    A mass migration event in which some fraction of the population in
    one deme simultaneously move to another deme, viewed backwards in
    time. For each lineage currently present in the source population,
    they move to the destination population with probability equal to
    ``proportion``.

    This event class generalises the population split (``-ej``) and
    admixture (``-es``) events from ``ms``. Note that MassMigrations
    do *not* have any side effects on the migration matrix.

    :param float time: The time at which this event occurs in generations.
    :param int source: The ID of the source population.
    :param int dest: The ID of the destination population.
    :param float proportion: The probability that any given lineage within
        the source population migrates to the destination population.
    """
    def __init__(self, time, source, dest=None, proportion=1.0, destination=None):
        super(MassMigration, self).__init__("mass_migration", time)
        if dest is not None and destination is not None:
            raise ValueError(
                "dest and destination are aliases; cannot supply both")
        if destination is not None:
            dest = destination
        self.source = source
        self.dest = dest
        self.proportion = proportion

    def get_ll_representation(self, num_populations):
        return {
            "type": self.type,
            "time": self.time,
            "source": self.source,
            "dest": self.dest,
            "proportion": self.proportion
        }

    def __str__(self):
        return (
            "Mass migration: lineages move from {} to {} with "
            "probability {}".format(
                self.source, self.dest, self.proportion))


class SimpleBottleneck(DemographicEvent):
    # This is an unsupported/undocumented demographic event.
    def __init__(self, time, population=None, proportion=1.0, population_id=None):
        super(SimpleBottleneck, self).__init__("simple_bottleneck", time)
        if population_id is not None and population is not None:
            raise ValueError(
                "population_id and population are aliases; cannot supply both.")
        if population_id is not None:
            population = population_id
        self.population = population
        self.proportion = proportion

    def get_ll_representation(self, num_populations):
        return {
            "type": self.type,
            "time": self.time,
            "population": self.population,
            "proportion": self.proportion
        }

    def __str__(self):
        return (
            "Simple bottleneck: lineages in population {} coalesce "
            "probability {}".format(self.population, self.proportion))


class InstantaneousBottleneck(DemographicEvent):
    # TODO document

    def __init__(self, time, population=None, strength=1.0, population_id=None):
        super(InstantaneousBottleneck, self).__init__("instantaneous_bottleneck", time)
        if population_id is not None and population is not None:
            raise ValueError(
                "population_id and population are aliases; cannot supply both.")
        if population_id is not None:
            population = population_id
        self.population = population
        self.strength = strength

    def get_ll_representation(self, num_populations):
        return {
            "type": self.type,
            "time": self.time,
            "population": self.population,
            "strength": self.strength,
        }

    def __str__(self):
        return (
            "Instantaneous bottleneck in population {}: equivalent to {} "
            "generations of the coalescent".format(
                self.population, self.strength))


class SimulationModel(object):
    """
    Superclass of all simulation models.
    """
    name = None

    def __init__(self, population_size=1):
        self.population_size = population_size

    def get_ll_representation(self):
        return {"name": self.name, "population_size": self.population_size}


class StandardCoalescent(SimulationModel):
    """
    The classical coalescent with recombination model (i.e., Hudson's algorithm).
    """
    name = "hudson"


class SmcApproxCoalescent(SimulationModel):
    # TODO document
    name = "smc"


class SmcPrimeApproxCoalescent(SimulationModel):
    # TODO document
    name = "smc_prime"


class DiscreteTimeWrightFisher(SimulationModel):
    """
    A discrete backwards-time Wright Fisher model, with back-and-forth
    recombination
    """
    name = 'dtwf'


class ParametricSimulationModel(SimulationModel):
    """
    The superclass of simulation models that require extra parameters.
    """
    def get_ll_representation(self):
        d = super(ParametricSimulationModel, self).get_ll_representation()
        d.update(self.__dict__)
        return d


class BetaCoalescent(ParametricSimulationModel):
    # TODO document.
    name = "beta"

    # TODO what is a meaningful value for this parameter? Ideally, the default
    # would be the equivalent of the Kingman coalescent or something similar.
    def __init__(self, population_size=1, alpha=1, truncation_point=None):
        self.population_size = population_size
        self.alpha = alpha
        if truncation_point is None:
            truncation_point = sys.float_info.max
        self.truncation_point = truncation_point


class DiracCoalescent(ParametricSimulationModel):
    # TODO document
    name = "dirac"

    # TODO What is a meaningful default for this value? See above.
    def __init__(self, population_size=1, psi=0.5, c=10.0):
        self.population_size = population_size
        self.psi = psi
        self.c = c


class Population(object):
    """
    Simple class to represent the state of a population in terms of its
    demographic parameters.
    """
    def __init__(self, start_size, end_size, growth_rate):
        self.start_size = start_size
        self.end_size = end_size
        self.growth_rate = growth_rate

    def __repr__(self):
        return repr(self.__dict__)


class Epoch(object):
    """
    Represents a single epoch in the simulation within which the state
    of the demographic parameters are constant.
    """
    def __init__(
            self, start_time=None, end_time=None, populations=None,
            migration_matrix=None, demographic_events=None):
        self.start_time = start_time
        self.end_time = end_time
        self.populations = populations
        self.migration_matrix = migration_matrix
        self.demographic_events = demographic_events

    def __repr__(self):
        return repr(self.__dict__)


class DemographyDebugger(object):
    """
    A class to facilitate debugging of population parameters and migration
    rates in the past.
    """
    def __init__(
            self, Ne=1, population_configurations=None, migration_matrix=None,
            demographic_events=[]):
        self._precision = 3
        # Make sure that we have a sample size of at least 2 so that we can
        # initialise the simulator.
        sample_size = None
        if population_configurations is None:
            sample_size = 2
        else:
            saved_sample_sizes = [
                pop_config.sample_size for pop_config in population_configurations]
            for pop_config in population_configurations:
                pop_config.sample_size = 2
        simulator = simulator_factory(
            sample_size=sample_size, Ne=Ne,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        self._make_epochs(
            simulator, sorted(demographic_events, key=lambda e: e.time))

        if population_configurations is not None:
            # Restore the saved sample sizes.
            for pop_config, sample_size in zip(
                    population_configurations, saved_sample_sizes):
                pop_config.sample_size = sample_size

    def _make_epochs(self, simulator, demographic_events):
        self.epochs = []
        ll_sim = simulator.create_ll_instance()
        N = simulator.num_populations
        start_time = 0
        end_time = 0
        abs_tol = 1e-9
        event_index = 0
        while not math.isinf(end_time):
            events = []
            while (
                    event_index < len(demographic_events) and
                    almost_equal(
                        demographic_events[event_index].time,
                        start_time, abs_tol=abs_tol)):
                events.append(demographic_events[event_index])
                event_index += 1
            end_time = ll_sim.debug_demography()
            m = ll_sim.get_migration_matrix()
            migration_matrix = [[m[j * N + k] for k in range(N)] for j in range(N)]
            growth_rates = [
                conf["growth_rate"] for conf in ll_sim.get_population_configuration()]
            populations = [
                Population(
                    start_size=ll_sim.compute_population_size(j, start_time),
                    end_size=ll_sim.compute_population_size(j, end_time),
                    growth_rate=growth_rates[j]) for j in range(N)]
            self.epochs.append(Epoch(
                start_time, end_time, populations, migration_matrix, events))
            start_time = end_time

    def _print_populations(self, epoch, output):
        field_width = self._precision + 6
        growth_rate_field_width = 14
        sep_str = " | "
        N = len(epoch.migration_matrix)
        fmt = (
            "{id:<2} "
            "{start_size:^{field_width}}"
            "{end_size:^{field_width}}"
            "{growth_rate:>{growth_rate_field_width}}")
        print(fmt.format(
            id="", start_size="start", end_size="end",
            growth_rate="growth_rate", field_width=field_width,
            growth_rate_field_width=growth_rate_field_width), end=sep_str,
            file=output)
        for k in range(N):
            print("{0:^{1}}".format(k, field_width), end="", file=output)
        print(file=output)
        h = "-" * (field_width - 1)
        print(
            fmt.format(
                id="", start_size=h, end_size=h, growth_rate=h,
                field_width=field_width,
                growth_rate_field_width=growth_rate_field_width),
            end=sep_str, file=output)
        for k in range(N):
            s = "-" * (field_width - 1)
            print("{0:<{1}}".format(s, field_width), end="", file=output)
        print(file=output)
        for j, pop in enumerate(epoch.populations):
            s = (
                "{id:<2}|"
                "{start_size:^{field_width}.{precision}g}"
                "{end_size:^{field_width}.{precision}g}"
                "{growth_rate:>{growth_rate_field_width}.{precision}g}"
                ).format(
                    id=j,
                    start_size=pop.start_size,
                    end_size=pop.end_size,
                    growth_rate=pop.growth_rate,
                    precision=self._precision, field_width=field_width,
                    growth_rate_field_width=growth_rate_field_width)
            print(s, end=sep_str, file=output)
            for k in range(N):
                x = epoch.migration_matrix[j][k]
                print("{0:^{1}.{2}g}".format(
                    x, field_width, self._precision), end="", file=output)
            print(file=output)

    def print_history(self, output=sys.stdout):
        """
        Prints a summary of the history of the populations.
        """
        for epoch in self.epochs:
            if len(epoch.demographic_events) > 0:
                print("Events @ generation {}".format(epoch.start_time), file=output)
            for event in epoch.demographic_events:
                print("   -", event, file=output)
            s = "Epoch: {} -- {} generations".format(epoch.start_time, epoch.end_time)
            print("=" * len(s), file=output)
            print(s, file=output)
            print("=" * len(s), file=output)
            self._print_populations(epoch, output)
            print(file=output)
