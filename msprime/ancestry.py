#
# Copyright (C) 2015-2020 University of Oxford
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
Module responsible for defining and running ancestry simulations.
"""
import bisect
import collections
import copy
import gzip
import inspect
import logging
import os
import warnings

import attr
import numpy as np
import tskit

import _msprime
from . import core
from . import demography
from . import mutations
from . import provenance

logger = logging.getLogger(__name__)


Sample = collections.namedtuple("Sample", ["population", "time"])


def model_factory(model, reference_size=1):
    """
    Returns a simulation model corresponding to the specified model and
    reference size.
    - If model is None, the default simulation model is returned with the
      specified reference size.
    - If model is a string, return the corresponding model instance
      with the specified reference size
    - If model is an instance of SimulationModel, return a copy of it.
      In this case, if the model's reference_size is None, set it to the
      parameter reference size.
    - Otherwise return a type error.
    """
    if reference_size <= 0:
        raise ValueError("Model reference size must be > 0")
    model_map = {
        "hudson": StandardCoalescent(reference_size),
        "smc": SmcApproxCoalescent(reference_size),
        "smc_prime": SmcPrimeApproxCoalescent(reference_size),
        "dtwf": DiscreteTimeWrightFisher(reference_size),
        "wf_ped": WrightFisherPedigree(reference_size),
    }
    if model is None:
        model_instance = StandardCoalescent(reference_size)
    elif isinstance(model, str):
        lower_model = model.lower()
        if lower_model not in model_map:
            raise ValueError(
                "Model '{}' unknown. Choose from {}".format(
                    model, list(model_map.keys())
                )
            )
        model_instance = model_map[lower_model]
    elif isinstance(model, SimulationModel):
        model_instance = copy.copy(model)
        if model_instance.reference_size is None:
            model_instance.reference_size = reference_size
    else:
        raise TypeError(
            "Simulation model must be a string or an instance of SimulationModel"
        )
    return model_instance


def _check_population_configurations(population_configurations):
    err = (
        "Population configurations must be a list of PopulationConfiguration instances"
    )
    for config in population_configurations:
        if not isinstance(config, demography.PopulationConfiguration):
            raise TypeError(err)
        if config.initial_size is not None and config.initial_size <= 0:
            raise ValueError("Population size must be > 0")
        if config.sample_size is not None and config.sample_size < 0:
            raise ValueError("Sample size must be >= 0")


def _replicate_generator(
    sim, mutation_generator, num_replicates, provenance_dict, end_time
):
    """
    Generator function for the many-replicates case of the simulate
    function.
    """

    encoded_provenance = None
    # The JSON is modified for each replicate to insert the replicate number.
    # To avoid repeatedly encoding the same JSON (which can take milliseconds)
    # we insert a replaceable string.
    placeholder = "@@_MSPRIME_REPLICATE_INDEX_@@"
    if provenance_dict is not None:
        provenance_dict["parameters"]["replicate_index"] = placeholder
        encoded_provenance = provenance.json_encode_provenance(
            provenance_dict, num_replicates
        )

    for j in range(num_replicates):
        sim.run(end_time)
        replicate_provenance = None
        if encoded_provenance is not None:
            replicate_provenance = encoded_provenance.replace(
                f'"{placeholder}"', str(j)
            )
        tree_sequence = sim.get_tree_sequence(mutation_generator, replicate_provenance)
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
    pedigree=None,
    migration_matrix=None,
    samples=None,
    demographic_events=(),
    model=None,
    record_migrations=False,
    from_ts=None,
    start_time=None,
    end_time=None,
    record_full_arg=False,
    num_labels=None,
    gene_conversion_rate=None,
    gene_conversion_track_length=None,
):
    """
    Convenience method to create a simulator instance using the same
    parameters as the `simulate` function. Primarily used for testing.
    """
    condition = (
        sample_size is None
        and population_configurations is None
        and samples is None
        and from_ts is None
    )
    if condition:
        raise ValueError(
            "Either sample_size, population_configurations, samples or from_ts must "
            "be specified"
        )
    the_samples = None
    if sample_size is not None:
        if samples is not None:
            raise ValueError("Cannot specify sample size and samples simultaneously.")
        if population_configurations is not None:
            raise ValueError(
                "Cannot specify sample size and population_configurations "
                "simultaneously."
            )
        s = Sample(population=0, time=0.0)
        # In pedigrees samples are diploid individuals
        if pedigree is not None:
            sample_size *= 2  # TODO: Update for different ploidy
        the_samples = [s for _ in range(sample_size)]
    # If we have population configurations we may have embedded sample_size
    # values telling us how many samples to take from each population.
    if population_configurations is not None:
        if pedigree is not None:
            raise NotImplementedError(
                "Cannot yet specify population configurations "
                "and pedigrees simultaneously"
            )
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
                        "and samples simultaneously"
                    )
            the_samples = samples
    elif samples is not None:
        the_samples = samples

    if start_time is not None and start_time < 0:
        raise ValueError("start_time cannot be negative")

    if num_labels is not None and num_labels < 1:
        raise ValueError("Must have at least one structured coalescent label")

    if from_ts is not None:
        if not isinstance(from_ts, tskit.TreeSequence):
            raise TypeError("from_ts must be a TreeSequence instance.")
        population_mismatch_message = (
            "Mismatch in the number of populations in from_ts and simulation "
            "parameters. The number of populations in the simulation must be "
            "equal to the number of populations in from_ts"
        )
        if population_configurations is None:
            if from_ts.num_populations != 1:
                raise ValueError(population_mismatch_message)
        else:
            if from_ts.num_populations != len(population_configurations):
                raise ValueError(population_mismatch_message)

    if recombination_map is None:
        # Default to 1 if no from_ts; otherwise default to the sequence length
        # of from_ts
        if from_ts is None:
            the_length = 1 if length is None else length
        else:
            the_length = from_ts.sequence_length if length is None else length
        the_rate = 0 if recombination_rate is None else recombination_rate
        if the_length <= 0:
            raise ValueError("Cannot provide non-positive sequence length")
        if the_rate < 0:
            raise ValueError("Cannot provide negative recombination rate")
        recomb_map = RecombinationMap.uniform_map(the_length, the_rate, discrete=False)
    else:
        if length is not None or recombination_rate is not None:
            raise ValueError(
                "Cannot specify length/recombination_rate along with "
                "a recombination map"
            )
        recomb_map = recombination_map

    # FIXME check the valid inputs for GC. Should we allow it when we
    # have a non-trivial genetic map?
    if gene_conversion_rate is None:
        gene_conversion_rate = 0
    else:
        if not recomb_map.discrete:
            raise ValueError(
                "Cannot specify gene_conversion_rate along with "
                "a nondiscrete recombination map"
            )
    if gene_conversion_track_length is None:
        gene_conversion_track_length = 1

    if from_ts is not None:
        if from_ts.sequence_length != recomb_map.get_length():
            raise ValueError(
                "The simulated sequence length must be the same as "
                "from_ts.sequence_length"
            )

    sim = Simulator(the_samples, recomb_map, model, Ne, from_ts)
    sim.store_migrations = record_migrations
    sim.store_full_arg = record_full_arg
    sim.start_time = start_time
    sim.num_labels = num_labels
    sim.gene_conversion_rate = gene_conversion_rate
    sim.gene_conversion_track_length = gene_conversion_track_length
    rng = random_generator
    if rng is None:
        rng = _msprime.RandomGenerator(core.get_random_seed())
    sim.random_generator = rng
    if population_configurations is not None:
        sim.set_population_configurations(population_configurations)
    if migration_matrix is not None:
        sim.set_migration_matrix(migration_matrix)
    if demographic_events is not None:
        sim.set_demographic_events(demographic_events)
    if pedigree is not None:
        if not isinstance(sim.model, WrightFisherPedigree):
            raise NotImplementedError("Pedigree can only be specified for wf_ped model")
        sim.set_pedigree(pedigree)
    return sim


def simulate(
    sample_size=None,
    Ne=1,
    length=None,
    recombination_rate=None,
    recombination_map=None,
    mutation_rate=None,
    population_configurations=None,
    pedigree=None,
    migration_matrix=None,
    demographic_events=(),
    samples=None,
    model=None,
    record_migrations=False,
    random_seed=None,
    mutation_generator=None,
    num_replicates=None,
    replicate_index=None,
    from_ts=None,
    start_time=None,
    end_time=None,
    record_full_arg=False,
    num_labels=None,
    record_provenance=True,
    # FIXME add documentation for these.
    gene_conversion_rate=None,
    gene_conversion_track_length=None,
):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns the resulting :class:`tskit.TreeSequence`. Note that
    Ne is the effective diploid population size (so the effective number
    of genomes in the population is 2*Ne), but ``sample_size`` is the
    number of (monoploid) genomes sampled.

    :param int sample_size: The number of sampled monoploid genomes.  If not
        specified or None, this defaults to the sum of the subpopulation sample
        sizes. Either ``sample_size``, ``population_configurations`` or
        ``samples`` must be specified.
    :param float Ne: The effective (diploid) population size for the reference
        population. This defaults to 1 if not specified.
        Please see the :ref:`sec_api_simulation_models` section for more details
        on specifying simulations models.
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
    :param float mutation_rate: The rate of infinite sites
        mutations per unit of sequence length per generation.
        If not specified, no mutations are generated. This option only
        allows for infinite sites mutations with a binary (i.e., 0/1)
        alphabet. For more control over the mutational process, please
        use the :func:`.mutate` function.
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
        :math:`N \\times N` matrix with 0 on the diagonal, consisting of
        :math:`N` lists of length :math:`N` or an :math:`N \\times N` numpy
        array, with the [j, k]th element giving the fraction of
        population j that consists of migrants from population k in each
        generation.
    :param list demographic_events: The list of demographic events to
        simulate. Demographic events describe changes to the populations
        in the past. Events should be supplied in non-decreasing
        order of time in the past. Events with the same time value will be
        applied sequentially in the order that they were supplied before the
        simulation algorithm continues with the next time step.
    :param list samples: The list specifying the location and time of
        all samples. This parameter may be used to specify historical
        samples, and cannot be used in conjunction with the ``sample_size``
        parameter. Each sample is a (``population``, ``time``) pair
        such that the sample in position ``j`` in the list of samples
        is drawn in the specified population at the specfied time. Time
        is measured in generations ago, as elsewhere.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
    :param int num_replicates: The number of replicates of the specified
        parameters to simulate. If this is not specified or None,
        no replication is performed and a :class:`tskit.TreeSequence` object
        returned. If :obj:`num_replicates` is provided, the specified
        number of replicates is performed, and an iterator over the
        resulting :class:`tskit.TreeSequence` objects returned.
    :param int replicate_index: Return only a specific tree
        sequence from the set of replicates. This is used to recreate a specific tree
        sequence from e.g. provenance. This argument only makes sense when used with
        `random seed`, and is not compatible with `num_replicates`. Note also that
        msprime will have to create and discard all the tree sequences up to this index.
    :param tskit.TreeSequence from_ts: If specified, initialise the simulation
        from the root segments of this tree sequence and return the
        completed tree sequence. Please see :ref:`here
        <sec_api_simulate_from>` for details on the required properties
        of this tree sequence and its interactions with other parameters.
        (Default: None).
    :param float start_time: If specified, set the initial time that the
        simulation starts to this value. If not specified, the start
        time is zero if performing a simulation of a set of samples,
        or is the time of the oldest node if simulating from an
        existing tree sequence (see the ``from_ts`` parameter).
    :param float end_time: If specified, terminate the simulation at the
        specified time. In the returned tree sequence, all rootward paths from
        samples with time < end_time will end in a node with one child with
        time equal to end_time. Sample nodes with time >= end_time will
        also be present in the output tree sequence. If not specified or ``None``,
        run the simulation until all samples have an MRCA at all positions in
        the genome.
    :param bool record_full_arg: If True, record all intermediate nodes
        arising from common ancestor and recombination events in the output
        tree sequence. This will result in unary nodes (i.e., nodes in marginal
        trees that have only one child). Defaults to False.
    :param model: The simulation model to use.
        This can either be a string (e.g., ``"smc_prime"``) or an instance of
        a simulation model class (e.g, ``msprime.DiscreteTimeWrightFisher(100)``.
        Please see the :ref:`sec_api_simulation_models` section for more details
        on specifying simulations models.
    :type model: str or simulation model instance
    :param bool record_provenance: If True, record all configuration and parameters
        required to recreate the tree sequence. These can be accessed
        via ``TreeSequence.provenances()``).
    :return: The :class:`tskit.TreeSequence` object representing the results
        of the simulation if no replication is performed, or an
        iterator over the independent replicates simulated if the
        :obj:`num_replicates` parameter has been used.
    :rtype: :class:`tskit.TreeSequence` or an iterator over
        :class:`tskit.TreeSequence` replicates.
    :warning: If using replication, do not store the results of the
        iterator in a list! For performance reasons, the same
        underlying object may be used for every TreeSequence
        returned which will most likely lead to unexpected behaviour.
    """

    seed = random_seed
    if random_seed is None:
        seed = core.get_random_seed()
    seed = int(seed)
    rng = _msprime.RandomGenerator(seed)

    provenance_dict = None
    if record_provenance:
        argspec = inspect.getargvalues(inspect.currentframe())
        # num_replicates is excluded as provenance is per replicate
        # replicate index is excluded as it is inserted for each replicate
        parameters = {
            "command": "simulate",
            **{
                arg: argspec.locals[arg]
                for arg in argspec.args
                if arg not in ["num_replicates", "replicate_index"]
            },
        }
        parameters["random_seed"] = seed
        provenance_dict = provenance.get_provenance_dict(parameters)

    sim = simulator_factory(
        sample_size=sample_size,
        random_generator=rng,
        Ne=Ne,
        length=length,
        recombination_rate=recombination_rate,
        recombination_map=recombination_map,
        population_configurations=population_configurations,
        pedigree=pedigree,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        samples=samples,
        model=model,
        record_migrations=record_migrations,
        from_ts=from_ts,
        start_time=start_time,
        record_full_arg=record_full_arg,
        num_labels=num_labels,
        gene_conversion_rate=gene_conversion_rate,
        gene_conversion_track_length=gene_conversion_track_length,
    )

    if mutation_generator is not None:
        # This error was added in version 0.6.1.
        raise ValueError(
            "mutation_generator is not longer supported. Please use "
            "msprime.mutate instead"
        )

    if mutation_rate is not None:
        # There is ambiguity in how we should throw mutations onto partially
        # built tree sequences: on the whole thing, or must the newly added
        # topology? Before or after start_time? We avoid this complexity by
        # asking the user to use mutate(), which should have the required
        # flexibility.
        if from_ts is not None:
            raise ValueError(
                "Cannot specify mutation rate combined with from_ts. Please use "
                "msprime.mutate on the final tree sequence instead"
            )
        # There is ambiguity in how the start_time argument should interact with
        # the mutation generator: should we throw mutations down on the whole
        # tree or just the (partial) edges after start_time? To avoid complicating
        # things here, make the user use mutate() which should have the flexibility
        # to do whatever is needed.
        if start_time is not None and start_time > 0:
            raise ValueError(
                "Cannot specify mutation rate combined with a non-zero "
                "start_time. Please use msprime.mutate on the returned "
                "tree sequence instead"
            )
    # TODO when the ``discrete`` parameter is added here, pass it through
    # to make it a property of the mutation generator.
    mutation_generator = mutations._simple_mutation_generator(
        mutation_rate, sim.sequence_length, sim.random_generator
    )
    if replicate_index is not None and random_seed is None:
        raise ValueError(
            "Cannot specify replicate_index without random_seed as this "
            "has the same effect as not specifying replicate_index i.e. a "
            "random tree sequence"
        )
    if replicate_index is not None and num_replicates is not None:
        raise ValueError(
            "Cannot specify replicate_index with num_replicates as only "
            "the replicate_index specified will be returned."
        )
    if num_replicates is None and replicate_index is None:
        replicate_index = 0
    if replicate_index is not None:
        iterator = _replicate_generator(
            sim, mutation_generator, replicate_index + 1, provenance_dict, end_time
        )
        # Return the last element of the iterator
        deque = collections.deque(iterator, maxlen=1)
        return deque.pop()
    else:
        return _replicate_generator(
            sim, mutation_generator, num_replicates, provenance_dict, end_time
        )


class Simulator:
    """
    Class to simulate trees under a variety of population models.
    """

    def __init__(
        self, samples, recombination_map, model="hudson", Ne=0.25, from_ts=None
    ):
        if from_ts is None:
            if len(samples) < 2:
                raise ValueError("Sample size must be >= 2")
            self.samples = samples
        else:
            if samples is not None and len(samples) > 0:
                raise ValueError("Cannot specify samples with from_ts")
            self.samples = []
        if not isinstance(recombination_map, RecombinationMap):
            raise TypeError("RecombinationMap instance required")
        self.ll_sim = None
        self.model = model_factory(model, Ne)
        self.recombination_map = recombination_map
        self.from_ts = from_ts
        self.start_time = None
        self.random_generator = None
        self.population_configurations = [
            demography.PopulationConfiguration(initial_size=self.model.reference_size)
        ]
        self.pedigree = None
        self.migration_matrix = [[0]]
        self.demographic_events = []
        self.model_change_events = []
        self.store_migrations = False
        self.store_full_arg = False
        self.num_labels = None
        # We always need at least n segments, so no point in making
        # allocation any smaller than this.
        num_samples = (
            len(self.samples) if self.samples is not None else from_ts.num_samples
        )
        block_size = 64 * 1024
        self.segment_block_size = max(block_size, num_samples)
        self.avl_node_block_size = block_size
        self.node_mapping_block_size = block_size
        self.end_time = None
        self.gene_conversion_rate = 0
        self.gene_conversion_track_length = 1

    @property
    def sequence_length(self):
        return self.recombination_map.get_sequence_length()

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
        return self.ll_sim.get_breakpoints()

    @property
    def time(self):
        return self.ll_sim.get_time()

    @property
    def num_avl_node_blocks(self):
        return self.ll_sim.get_num_avl_node_blocks()

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
    def num_gene_conversion_events(self):
        return self.ll_sim.get_num_gene_conversion_events()

    @property
    def num_populations(self):
        return len(self.population_configurations)

    @property
    def num_migration_events(self):
        return self.ll_sim.get_num_migration_events()

    @property
    def total_num_migration_events(self):
        return np.sum(self.num_migration_events)

    @property
    def num_multiple_recombination_events(self):
        return self.ll_sim.get_num_multiple_recombination_events()

    def set_migration_matrix(self, migration_matrix):
        err = (
            "migration matrix must be a N x N square matrix encoded "
            "as a list-of-lists or numpy array, where N is the number of populations "
            "defined in the population_configurations. The diagonal "
            "elements of this matrix must be zero. For example, a "
            "valid matrix for a 3 population system is "
            "[[0, 1, 1], [1, 0, 1], [1, 1, 0]]"
        )
        N = len(self.population_configurations)
        self.migration_matrix = np.array(migration_matrix)
        if self.migration_matrix.shape != (N, N):
            raise ValueError(err)

    def set_population_configurations(self, population_configurations):
        _check_population_configurations(population_configurations)
        self.population_configurations = population_configurations
        # For any populations configurations in which the initial size is None,
        # set it to the population size.
        for pop_conf in self.population_configurations:
            if pop_conf.initial_size is None:
                pop_conf.initial_size = self.model.reference_size
        # Now set the default migration matrix.
        N = len(self.population_configurations)
        self.migration_matrix = np.zeros((N, N))

    def set_pedigree(self, pedigree):
        if len(self.samples) % 2 != 0:
            raise ValueError(
                "In (diploid) pedigrees, must specify two " "lineages per individual."
            )

        if pedigree.is_sample is None:
            pedigree.set_samples(num_samples=len(self.samples) // 2)

        if sum(pedigree.is_sample) * 2 != len(self.samples):
            raise ValueError(
                "{} sample lineages to be simulated, but {} in pedigree".format(
                    len(self.samples), pedigree.num_samples * 2
                )
            )

        if not isinstance(self.model, WrightFisherPedigree):
            raise ValueError(
                "Can only specify pedigrees for the "
                "`msprime.WrightFisherPedigree` simulation model."
            )

        pedigree_max_time = np.max(pedigree.time)
        if len(self.demographic_events) > 0:
            de_min_time = min([x.time for x in self.demographic_events])
            if de_min_time <= pedigree_max_time:
                raise NotImplementedError(
                    "Demographic events must be older than oldest pedigree founder."
                )
        if len(self.model_change_events) > 0:
            mc_min_time = min([x.time for x in self.model_change_events])
            if mc_min_time < pedigree_max_time:
                raise NotImplementedError(
                    "Model change events earlier than founders of pedigree unsupported."
                )

        self.pedigree = pedigree.get_ll_representation()

    def set_demographic_events(self, demographic_events):
        err = (
            "Demographic events must be a list of DemographicEvent instances "
            "sorted in non-decreasing order of time."
        )
        self.demographic_events = []
        self.model_change_events = []
        for event in demographic_events:
            if not isinstance(event, demography.DemographicEvent):
                raise TypeError(err)
            if isinstance(event, demography.SimulationModelChange):
                # Take a copy so that we're not modifying our input params
                event = copy.copy(event)
                # Update the model so that we can parse strings and set
                # the reference size appropriately.
                event.model = model_factory(event.model, self.model.reference_size)
                self.model_change_events.append(event)
            else:
                self.demographic_events.append(event)

    def __choose_num_labels(self):
        """
        Choose the number of labels appropriately, given the simulation
        models that will be simulated.
        """
        self.num_labels = 1
        models = [self.model] + [event.model for event in self.model_change_events]
        for model in models:
            if isinstance(model, SweepGenicSelection):
                self.num_labels = 2

    def create_ll_instance(self):
        if self.num_labels is None:
            self.__choose_num_labels()
        # Now, convert the high-level values into their low-level
        # counterparts.
        ll_simulation_model = self.model.get_ll_representation()
        logger.debug("Setting initial model %s", ll_simulation_model)
        ll_population_configuration = [
            conf.asdict() for conf in self.population_configurations
        ]
        ll_demographic_events = [
            event.get_ll_representation() for event in self.demographic_events
        ]
        ll_recomb_map = self.recombination_map.get_ll_recombination_map()
        self.ll_tables = _msprime.LightweightTableCollection(
            self.recombination_map.get_sequence_length()
        )
        if self.from_ts is not None:
            from_ts_tables = self.from_ts.tables.asdict()
            self.ll_tables.fromdict(from_ts_tables)
        start_time = -1 if self.start_time is None else self.start_time
        ll_sim = _msprime.Simulator(
            samples=self.samples,
            recombination_map=ll_recomb_map,
            tables=self.ll_tables,
            start_time=start_time,
            random_generator=self.random_generator,
            model=ll_simulation_model,
            migration_matrix=self.migration_matrix,
            population_configuration=ll_population_configuration,
            pedigree=self.pedigree,
            demographic_events=ll_demographic_events,
            store_migrations=self.store_migrations,
            store_full_arg=self.store_full_arg,
            num_labels=self.num_labels,
            segment_block_size=self.segment_block_size,
            avl_node_block_size=self.avl_node_block_size,
            node_mapping_block_size=self.node_mapping_block_size,
            gene_conversion_rate=self.gene_conversion_rate,
            gene_conversion_track_length=self.gene_conversion_track_length,
        )
        return ll_sim

    def run(self, end_time=None):
        """
        Runs the simulation until complete coalescence has occurred.
        """
        if self.ll_sim is None:
            self.ll_sim = self.create_ll_instance()
        for event in self.model_change_events:
            # If the event time is a callable, we compute the end_time
            # as a function of the current simulation time.
            current_time = self.ll_sim.get_time()
            model_start_time = event.time
            if callable(event.time):
                model_start_time = event.time(current_time)
            # If model_start_time is None, we run until the current
            # model completes. Note that when event.time is a callable
            # it can also return None for this behaviour.
            if model_start_time is None:
                model_start_time = np.inf
            if model_start_time < current_time:
                raise ValueError(
                    "Model start times out of order or not computed correctly. "
                    f"current time = {current_time}; start_time = {model_start_time}"
                )
            logger.debug("Running simulation until maximum: %f", model_start_time)
            self.ll_sim.run(model_start_time)
            ll_new_model = event.model.get_ll_representation()
            logger.debug("Changing to model %s", ll_new_model)
            self.ll_sim.set_model(ll_new_model)
        end_time = np.inf if end_time is None else end_time
        logger.debug("Running simulation until maximum: %f", end_time)
        self.ll_sim.run(end_time)
        self.ll_sim.finalise_tables()

    def get_tree_sequence(self, mutation_generator=None, provenance_record=None):
        """
        Returns a TreeSequence representing the state of the simulation.
        """
        if mutation_generator is not None:
            mutation_generator.generate(self.ll_tables)
        tables = tskit.TableCollection.fromdict(self.ll_tables.asdict())
        if provenance_record is not None:
            tables.provenances.add_row(provenance_record)
        if self.from_ts is None:
            # Add the populations with metadata
            assert len(tables.populations) == len(self.population_configurations)
            tables.populations.clear()
            for pop_config in self.population_configurations:
                tables.populations.add_row(metadata=pop_config.encode_metadata())
        return tables.tree_sequence()

    def reset(self):
        """
        Resets the simulation so that we can perform another replicate.
        """
        if self.ll_sim is not None:
            self.ll_sim.reset()


class RecombinationMap:
    """
    A RecombinationMap represents the changing rates of recombination
    along a chromosome. This is defined via two lists of numbers:
    ``positions`` and ``rates``, which must be of the same length.
    Given an index j in these lists, the rate of recombination
    per base per generation is ``rates[j]`` over the interval
    ``positions[j]`` to ``positions[j + 1]``. Consequently, the first
    position must be zero, and by convention the last rate value
    is also required to be zero (although it is not used).

    .. warning::
        The ``num_loci`` parameter is deprecated. To set a discrete number of
        possible recombination sites along the sequence, scale ``positions``
        to the desired number of sites and set ``discrete=True`` to ensure
        recombination occurs only at integer values.

    :param list positions: The positions (in bases) denoting the
        distinct intervals where recombination rates change. These can
        be floating point values.
    :param list rates: The list of rates corresponding to the supplied
        ``positions``. Recombination rates are specified per base,
        per generation.
    :param int num_loci: **This parameter is deprecated**.
        The maximum number of non-recombining loci
        in the underlying simulation. By default this is set to
        the largest possible value, allowing the maximum resolution
        in the recombination process. However, for a finite sites
        model this can be set to smaller values.
    :param bool discrete: Whether recombination can occur only at integer
        positions. When ``False``, recombination sites can take continuous
        values. To simulate a fixed number of loci, set this parameter to
        ``True`` and scale ``positions`` to span the desired number of loci.
    """

    def __init__(self, positions, rates, num_loci=None, discrete=False, map_start=0):
        if num_loci is not None:
            if num_loci == positions[-1]:
                warnings.warn("num_loci is no longer supported and should not be used.")
            else:
                raise ValueError(
                    "num_loci does not match sequence length. "
                    "To set a discrete number of recombination sites, "
                    "scale positions to span the desired number of loci "
                    "and set discrete=True"
                )
        self._ll_recombination_map = _msprime.RecombinationMap(
            positions, rates, discrete
        )
        self.map_start = map_start

    @classmethod
    def uniform_map(cls, length, rate, num_loci=None, discrete=False):
        """
        Returns a :class:`.RecombinationMap` instance in which the recombination
        rate is constant over a chromosome of the specified length. The optional
        ``discrete`` controls whether recombination sites can occur only on integer
        positions or can take continuous values. The legacy ``num_loci`` option is
        no longer supported and should not be used.

        The following map can be used to simulate a true finite locus model
        with a fixed number of loci ``m``::

            >>> recomb_map = RecombinationMap.uniform_map(m, rate, discrete=True)

        :param float length: The length of the chromosome.
        :param float rate: The rate of recombination per unit of sequence length
            along this chromosome.
        :param int num_loci: This parameter is no longer supported.
        :param bool discrete: Whether recombination can occur only at integer
            positions. When ``False``, recombination sites can take continuous
            values. To simulate a fixed number of loci, set this parameter to
            ``True`` and set ``length`` to the desired number of loci.
        """
        return cls([0, length], [rate, 0], num_loci=num_loci, discrete=discrete)

    @classmethod
    def read_hapmap(cls, filename):
        """
        Parses the specified file in HapMap format. These files must be
        white-space-delimited, and contain a single header line (which is
        ignored), and then each subsequent line contains the starting position
        and recombination rate for the segment from that position (inclusive)
        to the starting position on the next line (exclusive). Starting
        positions of each segment are given in units of bases, and
        recombination rates in centimorgans/Megabase. The first column in this
        file is ignored, as are additional columns after the third (Position is
        assumed to be the second column, and Rate is assumed to be the third).
        If the first starting position is not equal to zero, then a
        zero-recombination region is inserted at the start of the chromosome.

        A sample of this format is as follows::

            Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
            chr1	55550	        2.981822	0.000000
            chr1	82571	        2.082414	0.080572
            chr1	88169	        2.081358	0.092229
            chr1	254996	        3.354927	0.439456
            chr1	564598	        2.887498	1.478148
            ...
            chr1	182973428	2.512769	122.832331
            chr1	183630013	0.000000	124.482178

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
                    map_start = pos
                    if pos != 0:
                        positions.append(0)
                        rates.append(0)
                positions.append(pos)
                # Rate is expressed in centimorgans per megabase, which
                # we convert to per-base rates
                rates.append(rate * 1e-8)
            if rate != 0:
                raise ValueError(
                    "The last rate provided in the recombination map must be zero"
                )
        finally:
            f.close()
        return cls(positions, rates, map_start=map_start)

    @property
    def mean_recombination_rate(self):
        """
        Return the weighted mean recombination rate
        across all windows of the entire recombination map.
        """
        chrom_length = self._ll_recombination_map.get_sequence_length()

        positions = self._ll_recombination_map.get_positions()
        positions_diff = self._ll_recombination_map.get_positions()[1:]
        positions_diff = np.append(positions_diff, chrom_length)
        window_sizes = positions_diff - positions

        weights = window_sizes / chrom_length
        if self.map_start != 0:
            weights[0] = 0
        rates = self._ll_recombination_map.get_rates()

        return np.average(rates, weights=weights)

    def slice(self, start=None, end=None, trim=False):  # noqa: A003
        """
        Returns a subset of this recombination map between the specified end
        points. If start is None, it defaults to 0. If end is None, it defaults
        to the end of the map. If trim is True, remove the flanking
        zero recombination rate regions such that the sequence length of the
        new recombination map is end - start.
        """
        positions = self.get_positions()
        rates = self.get_rates()

        if start is None:
            i = 0
            start = 0
        if end is None:
            end = positions[-1]
            j = len(positions)

        if (
            start < 0
            or end < 0
            or start > positions[-1]
            or end > positions[-1]
            or start > end
        ):
            raise IndexError(f"Invalid subset: start={start}, end={end}")

        if start != 0:
            i = bisect.bisect_left(positions, start)
            if start < positions[i]:
                i -= 1
        if end != positions[-1]:
            j = bisect.bisect_right(positions, end, lo=i)

        new_positions = list(positions[i:j])
        new_rates = list(rates[i:j])
        new_positions[0] = start
        if end > new_positions[-1]:
            new_positions.append(end)
            new_rates.append(0)
        else:
            new_rates[-1] = 0
        if trim:
            new_positions = [pos - start for pos in new_positions]
        else:
            if new_positions[0] != 0:
                if new_rates[0] == 0:
                    new_positions[0] = 0
                else:
                    new_positions.insert(0, 0)
                    new_rates.insert(0, 0.0)
            if new_positions[-1] != positions[-1]:
                new_positions.append(positions[-1])
                new_rates.append(0)
        return self.__class__(new_positions, new_rates, discrete=self.discrete)

    def __getitem__(self, key):
        """
        Use slice syntax for obtaining a recombination map subset. E.g.
            >>> recomb_map_4m_to_5m = recomb_map[4e6:5e6]
        """
        if not isinstance(key, slice) or key.step is not None:
            raise TypeError("Only interval slicing is supported")
        start, end = key.start, key.stop
        if start is not None and start < 0:
            start += self.get_sequence_length()
        if end is not None and end < 0:
            end += self.get_sequence_length()
        return self.slice(start=start, end=end, trim=True)

    def get_ll_recombination_map(self):
        return self._ll_recombination_map

    def physical_to_genetic(self, physical_x):
        return self._ll_recombination_map.position_to_mass(physical_x)

    def physical_to_discrete_genetic(self, physical_x):
        raise ValueError("Discrete genetic space is no longer supported")

    def genetic_to_physical(self, genetic_x):
        return self._ll_recombination_map.mass_to_position(genetic_x)

    def get_total_recombination_rate(self):
        return self._ll_recombination_map.get_total_recombination_rate()

    def get_per_locus_recombination_rate(self):
        raise ValueError("Genetic loci are no longer supported")

    def get_size(self):
        return self._ll_recombination_map.get_size()

    def get_num_loci(self):
        raise ValueError("num_loci is no longer supported")

    def get_positions(self):
        # For compatability with existing code we convert to a list
        return list(self._ll_recombination_map.get_positions())

    def get_rates(self):
        # For compatability with existing code we convert to a list
        return list(self._ll_recombination_map.get_rates())

    def get_sequence_length(self):
        return self._ll_recombination_map.get_sequence_length()

    def get_length(self):
        # Deprecated: use sequence_length instead
        return self.get_sequence_length()

    @property
    def discrete(self):
        return self._ll_recombination_map.get_discrete()

    def asdict(self):
        return {
            "positions": self.get_positions(),
            "rates": self.get_rates(),
            "discrete": self.discrete,
            "map_start": self.map_start,
        }


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
            if isinstance(c, collections.Iterable):
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
            if isinstance(v, collections.Iterable):
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
        np.save(os.path.expanduser(fname), pedarray)

    def asdict(self):
        """
        Returns a dict of arguments to recreate this pedigree
        """
        return {
            key: getattr(self, key)
            for key in inspect.signature(self.__init__).parameters.keys()
            if hasattr(self, key)
        }


@attr.s
class SimulationModel:
    """
    Abstract superclass of all simulation models.
    """

    name = None

    reference_size = attr.ib(default=None)

    def get_ll_representation(self):
        return {"name": self.name, "reference_size": self.reference_size}

    def asdict(self):
        return attr.asdict(self)


class StandardCoalescent(SimulationModel):
    """
    The classical coalescent with recombination model (i.e., Hudson's algorithm).
    The string ``"hudson"`` can be used to refer to this model.

    This is the default simulation model.
    """

    name = "hudson"


class SmcApproxCoalescent(SimulationModel):
    """
    The original SMC model defined by McVean and Cardin. This
    model is implemented using a naive rejection sampling approach
    and so it may not be any more efficient to simulate than the
    standard Hudson model.

    The string ``"smc"`` can be used to refer to this model.
    """

    name = "smc"


class SmcPrimeApproxCoalescent(SimulationModel):
    """
    The SMC' model defined by Marjoram and Wall as an improvement on the
    original SMC. model is implemented using a naive rejection sampling
    approach and so it may not be any more efficient to simulate than the
    standard Hudson model.

    The string ``"smc_prime"`` can be used to refer to this model.
    """

    name = "smc_prime"


class DiscreteTimeWrightFisher(SimulationModel):
    """
    A discrete backwards-time Wright-Fisher model, with diploid back-and-forth
    recombination. The string ``"dtwf"`` can be used to refer to this model.

    Wright-Fisher simulations are performed very similarly to coalescent
    simulations, with all parameters denoting the same quantities in both
    models. Because events occur at discrete times however, the order in which
    they occur matters. Each generation consists of the following ordered
    events:

    - Migration events. As in the Hudson coalescent, these move single extant
      lineages between populations. Because migration events occur before
      lineages choose parents, migrant lineages choose parents from their new
      population in the same generation.
    - Demographic events. All events with `previous_generation < event_time <=
      current_generation` are carried out here.
    - Lineages draw parents. Each (monoploid) extant lineage draws a parent
      from their current population.
    - Diploid recombination. Each parent is diploid, so all child lineages
      recombine back-and-forth into the same two parental genome copies. These
      become two independent lineages in the next generation.
    - Historical sampling events. All historical samples with
      `previous_generation < sample_time <= current_generation` are inserted.

    """

    name = "dtwf"


class WrightFisherPedigree(SimulationModel):
    # TODO Complete documentation.
    """
    Backwards-time simulations through a pre-specified pedigree, with diploid
    individuals and back-and-forth recombination. The string ``"wf_ped"`` can
    be used to refer to this model.
    """
    name = "wf_ped"


class ParametricSimulationModel(SimulationModel):
    """
    The superclass of simulation models that require extra parameters.
    """

    def get_ll_representation(self):
        d = super().get_ll_representation()
        d.update(self.__dict__)
        return d


@attr.s
class BetaCoalescent(ParametricSimulationModel):
    """
    A diploid Xi-coalescent with up to four simultaneous multiple mergers and
    crossover recombination.

    There are two main differences between the Beta-Xi-coalescent and the
    standard coalescent. Firstly, the number of lineages that take part in each
    common ancestor event is random, with distribution determined by moments of
    the :math:`Beta(2 - \\alpha, \\alpha)`-distribution. In particular, when there
    are :math:`n` lineages, each set of :math:`k \\leq n` of them participates in a
    common ancestor event at rate

    .. math::
        \\frac{8}{B(2 - \\alpha, \\alpha)}
        \\int_0^1 x^{k - \\alpha - 1} (1 - x)^{n - k + \\alpha - 1} dx,

    where :math:`B(2 - \\alpha, \\alpha)` is the Beta-function.

    In a common ancestor event, all participating lineages are randomly split
    into four groups, corresponding to the four parental chromosomes in a diploid,
    bi-parental reproduction event. All lineages within each group merge simultaneously.

    .. warning::
        The prefactor of 8 in the common ancestor event rate arises as the product
        of two terms. A factor of 4 compensates for the fact that only one quarter
        of binary common ancestor events result in a merger due to diploidy.
        A further factor of 2 is included for consistency with the implementation
        of the Hudson model, in which :math:`n` lineages undergo binary mergers at
        rate :math:`n (n - 1)`.

    Secondly, the time scale predicted by the Beta-Xi-coalescent is proportional
    to :math:`N_e^{\\alpha - 1}` generations, where :math:`N_e` is the effective
    population size. Specifically, one unit of coalescent time corresponds to
    a number of generations given by

    .. math::
        \\frac{m^{\\alpha} N_e^{\\alpha - 1}}{\\alpha B(2 - \\alpha, \\alpha)},

    where

    .. math::
        m = 2 + \\frac{2^{\\alpha}}{3^{\\alpha - 1} (\\alpha - 1)}.

    Note that the time scale depends both on the effective population size :math:`N_e`
    and :math:`\\alpha`, and can be dramatically shorter than the timescale of the
    standard coalescent. Thus, effective population sizes must often be many orders
    of magnitude larger than census population sizes. The per-generation recombination
    rate is rescaled similarly to obtain the population-rescaled recombination rate.

    See `Schweinsberg (2003)
    <https://www.sciencedirect.com/science/article/pii/S0304414903000280>`_
    for the derivation of the common ancestor event rate, as well as the time scaling.
    Note however that the model of Schweinsberg (2003) is haploid, so that
    all participating lineages merge in a common ancestor event without
    splitting into four groups.

    :param float alpha: Determines the degree of skewness in the family size
        distribution, and must satisfy :math:`1 < \\alpha < 2`. Smaller values of
        :math:`\\alpha` correspond to greater skewness, and :math:`\\alpha = 2`
        would coincide with the standard coalescent.
    :param float truncation_point: Determines the maximum fraction of the
        population replaced by offspring in one reproduction event, and must
        satisfy :math:`0 < \\tau \\leq 1`, where :math:`\\tau` is the truncation point.
        The default is :math:`\\tau = 1`, which corresponds to the standard
        Beta-Xi-coalescent. When :math:`\\tau < 1`, the number of lineages
        participating in a common ancestor event is determined by moments
        of the :math:`Beta(2 - \\alpha, \\alpha)` distribution conditioned on not
        exceeding :math:`\\tau`, and the Beta-function in the expression
        for the time scale is also replaced by the incomplete Beta function
        :math:`Beta(\\tau; 2 - \\alpha, \\alpha)`.
    """

    name = "beta"

    alpha = attr.ib(default=None)
    truncation_point = attr.ib(default=1)


@attr.s
class DiracCoalescent(ParametricSimulationModel):
    """
    A diploid Xi-coalescent with up to four simultaneous multiple mergers and
    crossover recombination.

    The Dirac-Xi-coalescent is an implementation of the model of
    `Blath et al. (2013) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3527250/>`_
    The simulation proceeds similarly to the standard coalescent.
    In addition to binary common ancestor events at rate :math:`n (n - 1)` when
    there are :math:`n` lineages, potential multiple merger events take place
    at rate :math:`2 c > 0`. Each lineage participates in each multiple merger
    event independently with probability :math:`0 < \\psi \\leq 1`. All participating
    lineages are randomly split into four groups, corresponding to the four
    parental chromosomes present in a diploid, bi-parental reproduction event,
    and the lineages within each group merge simultaneously.

    .. warning::
        The Dirac-Xi-coalescent is obtained as a scaling limit of Moran models,
        rather than Wright-Fisher models. As a consequence, one unit of coalescent
        time is proportional to :math:`N_e^2` generations,
        rather than :math:`N_e` generations as in the standard coalescent.
        However, the coalescent recombination rate is obtained from the
        per-generation recombination probability by rescaling with
        :math:`N_e`. See :ref:`sec_tutorial_multiple_mergers`
        for an illustration of how this affects simulation output in practice.

    :param float c: Determines the rate of potential multiple merger events.
        We require :math:`c > 0`.
    :param float psi: Determines the fraction of the population replaced by
        offspring in one large reproduction event, i.e. one reproduction event
        giving rise to potential multiple mergers when viewed backwards in time.
        We require :math:`0 < \\psi \\leq 1`.
    """

    name = "dirac"

    psi = attr.ib(default=None)
    c = attr.ib(default=None)


@attr.s
class SweepGenicSelection(ParametricSimulationModel):
    # TODO document and finalise the API
    name = "sweep_genic_selection"

    position = attr.ib(default=None)
    start_frequency = attr.ib(default=None)
    end_frequency = attr.ib(default=None)
    alpha = attr.ib(default=None)
    dt = attr.ib(default=None)
