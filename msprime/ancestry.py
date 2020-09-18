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
import collections.abc
import copy
import inspect
import logging
import math
import struct
import sys

import attr
import numpy as np
import tskit

from . import core
from . import demography as demog
from . import intervals
from . import mutations
from . import provenance
from msprime import _msprime

logger = logging.getLogger(__name__)


def _model_factory(model):
    """
    Returns a simulation model corresponding to the specified model.
    - If model is None, the default simulation model is returned.
    - If model is a string, return the corresponding model instance.
    - If model is an instance of SimulationModel, return a copy of it.
    - Otherwise raise a type error.
    """
    model_map = {
        "hudson": StandardCoalescent(),
        "smc": SmcApproxCoalescent(),
        "smc_prime": SmcPrimeApproxCoalescent(),
        "dtwf": DiscreteTimeWrightFisher(),
        "wf_ped": WrightFisherPedigree(),
    }
    if model is None:
        model_instance = StandardCoalescent()
    elif isinstance(model, str):
        lower_model = model.lower()
        if lower_model not in model_map:
            raise ValueError(
                "Model '{}' unknown. Choose from {}".format(
                    model, list(model_map.keys())
                )
            )
        model_instance = model_map[lower_model]
    elif not isinstance(model, SimulationModel):
        raise TypeError(
            "Simulation model must be a string or an instance of SimulationModel"
        )
    else:
        model_instance = model
    return model_instance


def _parse_model_change_events(events):
    """
    Parses the specified list of events provided in model_arg[1:] into
    SimulationModelChange events. There are two different forms supported,
    and model descriptions are anything supported by model_factory.
    """
    err = (
        "Simulation model change events must be either a two-tuple "
        "(time, model), describing the time of the model change and "
        "the new model or be an instance of SimulationModelChange."
    )
    model_change_events = []
    for event in events:
        if isinstance(event, (tuple, list)):
            if len(event) != 2:
                raise ValueError(err)
            t = event[0]
            if t is not None:
                try:
                    t = float(t)
                except (TypeError, ValueError):
                    raise ValueError(
                        "Model change times must be either a floating point "
                        "value or None"
                    )
            event = SimulationModelChange(t, _model_factory(event[1]))
        elif isinstance(event, SimulationModelChange):
            # We don't want to modify our inputs, so take a deep copy.
            event = copy.copy(event)
            event.model = _model_factory(event.model)
        else:
            raise TypeError(err)
        model_change_events.append(event)
    return model_change_events


def _parse_model_arg(model_arg):
    """
    Parses the specified model argument from the simulate function,
    returning the initial model and any model change events.
    """
    err = (
        "The model argument must be either (a) a value that can be "
        "interpreted as a simulation model or (b) a list in which "
        "the first element is a model description and the remaining "
        "elements are model change events. These can either be described "
        "by a (time, model) tuple or SimulationModelChange instances."
    )
    if isinstance(model_arg, (list, tuple)):
        if len(model_arg) < 1:
            raise ValueError(err)
        model = _model_factory(model_arg[0])
        model_change_events = _parse_model_change_events(model_arg[1:])
    else:
        model = _model_factory(model_arg)
        model_change_events = []
    return model, model_change_events


def _filter_events(demographic_events):
    """
    Returns a tuple (demographic_events, model_change_events) which separates
    out the SimulationModelChange events from the list. This is to support the
    pre-1.0 syntax for model changes, where they were included in the
    demographic_events parameter.
    """
    filtered_events = []
    model_change_events = []
    for event in demographic_events:
        if isinstance(event, SimulationModelChange):
            model_change_events.append(event)
        else:
            filtered_events.append(event)
    # Make sure any model references are resolved.
    model_change_events = _parse_model_change_events(model_change_events)
    return filtered_events, model_change_events


def _check_population_configurations(population_configurations):
    err = (
        "Population configurations must be a list of PopulationConfiguration instances"
    )
    for config in population_configurations:
        if not isinstance(config, demog.PopulationConfiguration):
            raise TypeError(err)


def _samples_factory(sample_size, samples, population_configurations):
    """
    Returns a list of Sample objects, given the specified inputs.
    """
    the_samples = []
    if sample_size is not None:
        if samples is not None:
            raise ValueError("Cannot specify sample size and samples simultaneously.")
        if population_configurations is not None:
            raise ValueError(
                "Cannot specify sample size and population_configurations "
                "simultaneously."
            )
        s = demog.Sample(population=0, time=0.0)
        the_samples = [s for _ in range(sample_size)]
    # If we have population configurations we may have embedded sample_size
    # values telling us how many samples to take from each population.
    if population_configurations is not None:
        _check_population_configurations(population_configurations)
        if samples is None:
            the_samples = []
            for j, conf in enumerate(population_configurations):
                if conf.sample_size is not None:
                    the_samples += [demog.Sample(j, 0) for _ in range(conf.sample_size)]
        else:
            for conf in population_configurations:
                if conf.sample_size is not None:
                    raise ValueError(
                        "Cannot specify population configuration sample size"
                        " and samples simultaneously"
                    )
            the_samples = samples
    elif samples is not None:
        the_samples = samples
    return the_samples


def _demography_factory(
    Ne, demography, population_configurations, migration_matrix, demographic_events
):
    if demography is not None:
        if population_configurations is not None:
            raise ValueError(
                "The demography and population_configurations options "
                "cannot be used together"
            )
        if migration_matrix is not None:
            raise ValueError(
                "The demography and migration_matrix options cannot be used together"
            )
        if demographic_events is not None:
            raise ValueError(
                "The demography and demographic_events options cannot be used together"
            )
        # Take a copy so that we don't modify the input parameters when
        # resolving defaults
        demography = copy.deepcopy(demography)
    else:
        demography = demog.Demography.from_old_style(
            population_configurations, migration_matrix, demographic_events
        )

    # For any populations in which the initial size is None set it to Ne
    for pop in demography.populations:
        if pop.initial_size is None:
            pop.initial_size = Ne

    demography.validate()
    return demography


def _build_initial_tables(*, sequence_length, samples, ploidy, demography, pedigree):
    tables = tskit.TableCollection(sequence_length)

    if pedigree is None:
        for index, (population, time) in enumerate(samples):
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE, time=time, population=population
            )
            if population < 0:
                raise ValueError(f"Negative population ID in sample at index {index}")
            if population >= demography.num_populations:
                raise ValueError(
                    f"Invalid population reference '{population}' in sample "
                    f"at index {index}"
                )
    else:
        for parents, time, is_sample in zip(
            pedigree.parents, pedigree.time, pedigree.is_sample
        ):
            # We encode the parents in the metadata for now, but see
            # https://github.com/tskit-dev/tskit/issues/852
            encoded_parents = struct.pack("=ii", *parents)
            ind_id = tables.individuals.add_row(0, metadata=encoded_parents)
            node_flags = tskit.NODE_IS_SAMPLE if is_sample else 0
            for _ in range(ploidy):
                tables.nodes.add_row(node_flags, time, population=0, individual=ind_id)

    for population in demography.populations:
        md = population.temporary_hack_for_encoding_old_style_metadata()
        tables.populations.add_row(metadata=md)
    return tables


def _parse_simulate(
    sample_size=None,
    *,
    Ne=1,
    random_generator=None,
    random_seed=None,
    length=None,
    discrete_genome=None,
    recombination_rate=None,
    recombination_map=None,
    population_configurations=None,
    pedigree=None,
    migration_matrix=None,
    samples=None,
    demographic_events=None,
    model=None,
    record_migrations=False,
    from_ts=None,
    start_time=None,
    end_time=None,
    record_full_arg=False,
    num_labels=None,
    gene_conversion_rate=None,
    gene_conversion_track_length=None,
    demography=None,
    ploidy=None,
):
    """
    Argument parser for the simulate frontend. Interprets all the parameters
    and returns an appropriate instance of Simulator.
    """
    if Ne <= 0:
        raise ValueError("Population size must be positive")

    samples_specified = (
        sample_size is None
        and population_configurations is None
        and samples is None
        and from_ts is None
    )
    if samples_specified:
        # TODO remove the population_configurations message here if we're
        # deprecating it?
        raise ValueError(
            "Either sample_size, samples, population_configurations or from_ts must "
            "be specified"
        )
    samples = _samples_factory(sample_size, samples, population_configurations)

    model, model_change_events = _parse_model_arg(model)
    if demographic_events is not None:
        demographic_events, old_style_model_change_events = _filter_events(
            demographic_events
        )
        if len(old_style_model_change_events) > 0:
            if len(model_change_events) > 0:
                raise ValueError(
                    "Cannot specify SimulationModelChange events using both new-style "
                    "and pre 1.0 syntax"
                )
            model_change_events = old_style_model_change_events

    demography = _demography_factory(
        Ne, demography, population_configurations, migration_matrix, demographic_events
    )

    # The logic for checking from_ts and recombination map is bound together
    # in a complicated way, so we can factor them out into separate functions.
    if from_ts is None:
        if len(samples) < 2:
            raise ValueError("Sample size must be >= 2")
    else:
        if len(samples) > 0:
            raise ValueError("Cannot specify samples with from_ts")
        if not isinstance(from_ts, tskit.TreeSequence):
            raise TypeError("from_ts must be a TreeSequence instance.")
        if demography.num_populations != from_ts.num_populations:
            raise ValueError(
                "Mismatch in the number of populations in from_ts and simulation "
                "parameters. The number of populations in the simulation must be "
                "equal to the number of populations in from_ts"
            )

    if discrete_genome is None:
        discrete_genome = False

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
        if discrete_genome and length != math.floor(length):
            raise ValueError(
                "Cannot specify non-integer sequence length when discrete_genome "
                "is True"
            )
        recombination_map = intervals.RateMap.uniform(the_length, the_rate)
    else:
        if isinstance(recombination_map, intervals.RecombinationMap):
            # Convert from the legacy RecombinationMap class
            recombination_map = recombination_map.map
        elif not isinstance(recombination_map, intervals.RateMap):
            raise TypeError("RateMap instance required.")
        if length is not None or recombination_rate is not None:
            raise ValueError(
                "Cannot specify length/recombination_rate along with "
                "a recombination map"
            )
    if from_ts is not None:
        if recombination_map.sequence_length != from_ts.sequence_length:
            raise ValueError(
                "Recombination map and from_ts must have identical " "sequence_length"
            )

    if num_labels is not None and num_labels < 1:
        raise ValueError("Must have at least one structured coalescent label")

    # TODO should we allow a rate map in the same way as sim_ancestry?
    if gene_conversion_rate is None:
        gene_conversion_rate = 0
    else:
        if not discrete_genome:
            raise ValueError(
                "Cannot specify gene_conversion_rate along with "
                "a non discrete genome."
            )
    if gene_conversion_track_length is None:
        gene_conversion_track_length = 1
    gene_conversion_map = intervals.RateMap.uniform(
        recombination_map.sequence_length, gene_conversion_rate
    )

    # For the simulate code-path the rng will already be set, but
    # for convenience we allow it to be null to help with writing
    # tests. We also provide the random_seed argument for convenience.
    if random_seed is not None:
        if random_generator is not None:
            raise ValueError("Cannot specify both random_seed and random_generator")
        random_generator = _msprime.RandomGenerator(random_seed)
    if random_generator is None:
        random_generator = _msprime.RandomGenerator(core.get_random_seed())

    if ploidy is None:
        ploidy = 2

    if from_ts is None:
        tables = _build_initial_tables(
            sequence_length=recombination_map.sequence_length,
            samples=samples,
            ploidy=ploidy,
            demography=demography,
            pedigree=pedigree,
        )
    else:
        tables = from_ts.tables

    sim = Simulator(
        tables=tables,
        recombination_map=recombination_map,
        gene_conversion_map=gene_conversion_map,
        gene_conversion_track_length=gene_conversion_track_length,
        discrete_genome=discrete_genome,
        ploidy=ploidy,
        random_generator=random_generator,
        model=model,
        store_migrations=record_migrations,
        store_full_arg=record_full_arg,
        start_time=start_time,
        end_time=end_time,
        num_labels=num_labels,
        demography=demography,
        model_change_events=model_change_events,
    )
    return sim


def _parse_random_seed(seed):
    """
    Parse the specified random seed value and return a RandomGenerator
    object. If no seed is provided, generate a high-quality random seed.
    """
    if seed is None:
        seed = core.get_random_seed()
    seed = int(seed)
    return _msprime.RandomGenerator(seed)


def _build_provenance(command, random_seed, frame):
    """
    Builds a provenance dictionary suitable for use as the basis
    of tree sequence provenance in replicate simulations. Uses the
    specified stack frame to determine the values of the arguments
    passed in, with a few exceptions.
    """
    argspec = inspect.getargvalues(frame)
    # num_replicates is excluded as provenance is per replicate
    # replicate index is excluded as it is inserted for each replicate
    parameters = {
        "command": command,
        **{
            arg: argspec.locals[arg]
            for arg in argspec.args
            if arg not in ["num_replicates", "replicate_index"]
        },
    }
    parameters["random_seed"] = random_seed
    return provenance.get_provenance_dict(parameters)


def simulate(
    sample_size=None,
    *,
    Ne=1,
    length=None,
    discrete_genome=None,
    recombination_rate=None,
    recombination_map=None,
    mutation_rate=None,
    population_configurations=None,
    pedigree=None,
    migration_matrix=None,
    demographic_events=None,
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
    demography=None,
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
    :param bool discrete_genome: If True, use discrete coordinates
        in simulations such that recombination breakpoints and mutational
        sites can only occur at integer positions along the genome.
        Multiple mutations can occur at the same site.
        If False (the default) the recombination breakpoints and
        coordinates of mutational sites are continuous floating point values.
        All sites in the returned tree sequence will have exactly one mutation.
        Please see the :func:`.mutate` function for a more powerful approach to
        simulating mutations on a tree sequence. It is an error to specify
        the ``discrete_genome`` parameter at the same time as the
        ``recombination_map`` argument.
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
    :param list migration_matrix: The matrix describing the rates of migration
        between all pairs of populations. If :math:`N` populations are defined
        in the ``population_configurations`` parameter, then the migration
        matrix must be an :math:`N \\times N` matrix with 0 on the diagonal,
        consisting of :math:`N` lists of length :math:`N` or an :math:`N
        \\times N` numpy array. The :math:`[j, k]^{th}` element of the
        migration matrix gives the expected number of migrants moving from
        population :math:`k` to population :math:`j` per generation, divided by
        the size of population :math:`j`.  When simulating from the
        discrete-time Wright-Fisher model (``model = "dtwf"``), the row sums of
        the migration matrix must not exceed 1. There are no sum constraints for
        migration rates in continuous-time models.
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
    rng = _parse_random_seed(random_seed)
    provenance_dict = None
    if record_provenance:
        frame = inspect.currentframe()
        provenance_dict = _build_provenance("simulate", rng.seed, frame)

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
        mutation_rate = float(mutation_rate)

    sim = _parse_simulate(
        sample_size=sample_size,
        random_generator=rng,
        Ne=Ne,
        length=length,
        discrete_genome=discrete_genome,
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
        end_time=end_time,
        record_full_arg=record_full_arg,
        num_labels=num_labels,
        gene_conversion_rate=gene_conversion_rate,
        gene_conversion_track_length=gene_conversion_track_length,
        demography=demography,
    )

    return _wrap_replicates(
        sim,
        random_seed=random_seed,
        replicate_index=replicate_index,
        num_replicates=num_replicates,
        provenance_dict=provenance_dict,
        mutation_rate=mutation_rate,
    )


def _wrap_replicates(
    simulator,
    *,
    random_seed,
    replicate_index,
    num_replicates,
    provenance_dict,
    mutation_rate=None,
):
    """
    Wrapper for the logic used to run replicate simulations for the two
    frontends.
    """

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
        # This is the default case where we just return one replicate.
        replicate_index = 0
        num_replicates = 1
    if replicate_index is not None:
        num_replicates = replicate_index + 1

    iterator = simulator.run_replicates(
        num_replicates, mutation_rate=mutation_rate, provenance_dict=provenance_dict,
    )
    if replicate_index is not None:
        # Return the last element of the iterator
        deque = collections.deque(iterator, maxlen=1)
        return deque.pop()
    else:
        return iterator


def _parse_rate_map(rate_param, sequence_length, name):
    """
    Parse the specified input rate parameter value into a rate map.
    """
    # Note: in the future we might have another clause here where we
    # allow for a different  map per population. This could be
    # accepted as either a list of N rate maps, or a dictionary mapping
    # population names to maps.
    # See https://github.com/tskit-dev/msprime/issues/1095

    msg_head = f"Error in parsing rate map for {name}: "
    if isinstance(rate_param, intervals.RateMap):
        rate_map = rate_param
        if rate_map.sequence_length != sequence_length:
            raise ValueError(msg_head + "sequence_length must match")
    else:
        rate_param = 0 if rate_param is None else float(rate_param)
        rate_map = intervals.RateMap.uniform(sequence_length, rate_param)
    return rate_map


def _insert_integer_samples(n, ploidy, tables):
    """
    Insert n individuals with the specified ploidy into the tables.
    """
    # Doing this with numpy is slighty obscure, but it's *much* faster
    # than simple loops. For 10^7 samples with ploidy=2 we go from
    # several minutes to a few seconds.
    ind_flags = np.zeros(n, dtype=np.uint32)
    tables.individuals.set_columns(flags=ind_flags)
    node_individual = np.repeat(np.arange(n, dtype=np.int32), ploidy)
    N = n * ploidy
    tables.nodes.set_columns(
        flags=np.full(N, tskit.NODE_IS_SAMPLE, dtype=np.uint32),
        time=np.zeros(N),
        population=np.zeros(N, dtype=np.int32),
        individual=node_individual,
    )


def _parse_samples(samples, ploidy, tables):
    """
    Parse the specified "samples" value and insert them into the specified
    tables. How the samples argument is interpreted depends on the state
    of the tables.
    """
    if isinstance(samples, collections.abc.Iterable):
        # We check if the samples can be interpreted as a list first because
        # we may want to extend the semantics here later to include specifying
        # a set of individuals from a pedigree defined in the initial state.
        # The most natural way to do this would be a numpy array of individual
        # IDs. By checking the iterable condition first, we make sure that we're
        # not interpreting input numpy arrays of length 1 as a numeric argument.
        # Interpret as a list of Sample objects.

        # TODO this should probably be recast to using numpy input types for
        # efficiency here. We could regard the input array as a numpy struct
        # array, which would be a lot more efficient.
        # See https://github.com/tskit-dev/msprime/issues/1211
        for sample in samples:
            if not isinstance(sample, demog.Sample):
                raise TypeError("msprime.Sample object required")
            if not core.isinteger(sample.population):
                raise TypeError("Sample population references must be integers")
            population = int(sample.population)
            if population < 0:
                raise ValueError("Negative population ID")
            if sample.population >= len(tables.populations):
                raise ValueError("Sample population ID out of bounds")
            ind_id = tables.individuals.add_row(flags=0)
            for _ in range(ploidy):
                tables.nodes.add_row(
                    flags=tskit.NODE_IS_SAMPLE,
                    time=sample.time,
                    population=population,
                    individual=ind_id,
                )

    elif core.isinteger(samples):
        n = int(samples)
        if n < 0:
            raise ValueError("Cannot have a negative number of samples")
        if len(tables.populations) != 1:
            raise ValueError(
                "Numeric samples can only be used in single population models. "
                "Please use Demography.sample() to generate a list of samples "
                "for your model, which can be used instead."
            )
        _insert_integer_samples(n, ploidy, tables)
    else:
        raise TypeError(
            "The samples argument must be either an integer (for single population "
            "models) or a list of msprime.Sample objects. The Demography.sample() "
            "function is useful for generating samples for complex demographic "
            "models."
        )


def _parse_flag(value, *, default):
    """
    Parses a boolean flag, which can be either True, False, or None.
    If the input value is None, return the default. Otherwise,
    check that the input value is a bool.

    Note that we do *not* cast to a bool as this would accept
    truthy values like the empty list, etc. In this case None
    would be converted to False, potentially conflicting with
    the default value.
    """
    assert isinstance(default, bool)
    if value is None:
        return default
    if not isinstance(value, bool):
        raise TypeError("Boolean flag must be True, False, or None (the default value)")
    return value


def _parse_sim_ancestry(
    samples=None,
    *,
    sequence_length=None,
    recombination_rate=None,
    gene_conversion_rate=None,
    gene_conversion_track_length=None,
    discrete_genome=None,
    population_size=None,
    demography=None,
    ploidy=None,
    model=None,
    initial_state=None,
    start_time=None,
    end_time=None,
    record_migrations=None,
    record_full_arg=None,
    num_labels=None,
    random_seed=None,
    random_generator=None,
    num_replicates=None,
    replicate_index=None,
):
    """
    Argument parser for the sim_ancestry frontend. Interprets all the parameters
    and returns an appropriate instance of Simulator.
    """
    # As a general rule we try to cast any input value to the required types
    # early and in a way that provides an interpretable traceback.

    # Simple defaults.
    start_time = 0 if start_time is None else float(start_time)
    end_time = sys.float_info.max if end_time is None else float(end_time)
    discrete_genome = _parse_flag(discrete_genome, default=True)
    record_full_arg = _parse_flag(record_full_arg, default=False)
    record_migrations = _parse_flag(record_migrations, default=False)

    if initial_state is not None:
        if isinstance(initial_state, tskit.TreeSequence):
            initial_state = initial_state.dump_tables()
        elif not isinstance(initial_state, tskit.TableCollection):
            raise TypeError(
                "initial_state must either be a TreeSequence or TableCollection instance"
            )

    if sequence_length is None:
        # These are all the cases in which we derive the sequence_length
        # from somewhere else.
        if initial_state is not None:
            sequence_length = initial_state.sequence_length
        elif recombination_rate is None and gene_conversion_rate is None:
            # In this case, we're doing single-locus simulations, so a sequence
            # length of 1 makes sense.
            sequence_length = 1
        elif isinstance(recombination_rate, intervals.RateMap):
            sequence_length = recombination_rate.sequence_length
        elif isinstance(gene_conversion_rate, intervals.RateMap):
            sequence_length = gene_conversion_rate.sequence_length
        else:
            raise ValueError(
                "A sequence_length value must be specified. This can be either "
                "via the the sequence_length parameter itself, of implicitly "
                "through using a RateMap instance for the recombination_rate "
                "or gene_conversion_rate parameters, or via the initial_state "
                "tables. "
            )
    else:
        sequence_length = float(sequence_length)
    assert sequence_length is not None

    if discrete_genome and math.floor(sequence_length) != sequence_length:
        raise ValueError("Must have integer sequence length with discrete_genome=True")

    recombination_map = _parse_rate_map(
        recombination_rate, sequence_length, "recombination"
    )
    gene_conversion_map = _parse_rate_map(
        gene_conversion_rate, sequence_length, "gene conversion"
    )
    if gene_conversion_track_length is None:
        if gene_conversion_rate is None:
            # It doesn't matter what the track_length is, just set a
            # value to keep the low-level code happy.
            gene_conversion_track_length = 1
        else:
            raise ValueError(
                "Must specify track length when simulating gene conversion"
            )
    else:
        if gene_conversion_rate is None:
            raise ValueError(
                "Must specify gene conversion rate along with track length"
            )
        gene_conversion_track_length = float(gene_conversion_track_length)

    # Default to diploid
    ploidy = 2 if ploidy is None else ploidy
    if not core.isinteger(ploidy):
        raise TypeError("ploidy must be an integer")
    ploidy = int(ploidy)
    if ploidy < 1:
        raise ValueError("ploidy must be >= 1")

    model, model_change_events = _parse_model_arg(model)
    is_dtwf = isinstance(model, DiscreteTimeWrightFisher)

    # Check the demography. If no demography is specified, we default to a
    # single-population model with a given population size.
    if demography is None:
        if is_dtwf:
            # A default size of 1 isn't so smart for DTWF and almost certainly
            # an error.
            if population_size is None:
                raise ValueError(
                    "When using the DTWF model, the population size must be set "
                    "explicitly, either using the population_size or demography "
                    "arguments."
                )
        population_size = 1 if population_size is None else float(population_size)
        demography = demog.Demography.simple_model(population_size)
    elif isinstance(demography, demog.Demography):
        if population_size is not None:
            raise ValueError("Cannot specify demography and population size")
    else:
        raise TypeError("demography argument must be an instance of msprime.Demography")
    demography.validate()

    if initial_state is None:
        if samples is None:
            raise ValueError(
                "Either the samples or initial_state arguments must be provided"
            )
        initial_state = tskit.TableCollection(sequence_length)
        demography.insert_populations(initial_state)
        _parse_samples(samples, ploidy, initial_state)
    else:
        if samples is not None:
            raise ValueError("Cannot specify both samples and initial_state")
        if sequence_length != initial_state.sequence_length:
            raise ValueError(
                "The initial_state sequence length must be consistent with the"
                "value derived from either the sequence_length, "
                "recombination_rate or gene_conversion_rate parameters."
            )
        if len(initial_state.populations) == 0:
            raise ValueError(
                "initial_state tables must define at least one population."
            )

    # For the standard code path, random_generator should always be set. But,
    # we can also use the parser code elsewhere, so it's useful to be able
    # to either provide a random seed, or get a good one.
    if random_generator is None:
        random_generator = _parse_random_seed(random_seed)
    elif random_seed is not None:
        raise ValueError("Cannot specify both random seed and random generator")

    return Simulator(
        tables=initial_state,
        random_generator=random_generator,
        recombination_map=recombination_map,
        gene_conversion_map=gene_conversion_map,
        gene_conversion_track_length=gene_conversion_track_length,
        discrete_genome=discrete_genome,
        ploidy=ploidy,
        demography=demography,
        model=model,
        model_change_events=model_change_events,
        store_migrations=record_migrations,
        store_full_arg=record_full_arg,
        start_time=start_time,
        end_time=end_time,
        num_labels=num_labels,
    )


def sim_ancestry(
    samples=None,
    *,
    sequence_length=None,
    recombination_rate=None,
    gene_conversion_rate=None,
    gene_conversion_track_length=None,
    discrete_genome=None,
    population_size=None,
    demography=None,
    ploidy=None,
    model=None,
    initial_state=None,
    start_time=None,
    end_time=None,
    record_migrations=None,
    record_full_arg=None,
    num_labels=None,
    random_seed=None,
    num_replicates=None,
    replicate_index=None,
    record_provenance=None,
):

    random_generator = _parse_random_seed(random_seed)
    record_provenance = True if record_provenance is None else record_provenance
    provenance_dict = None
    if record_provenance:
        frame = inspect.currentframe()
        provenance_dict = _build_provenance(
            "sim_ancestry", random_generator.seed, frame
        )

    sim = _parse_sim_ancestry(
        samples=samples,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        gene_conversion_rate=gene_conversion_rate,
        gene_conversion_track_length=gene_conversion_track_length,
        discrete_genome=discrete_genome,
        population_size=population_size,
        demography=demography,
        ploidy=ploidy,
        model=model,
        initial_state=initial_state,
        start_time=start_time,
        end_time=end_time,
        record_migrations=record_migrations,
        record_full_arg=record_full_arg,
        num_labels=num_labels,
        random_generator=random_generator,
    )
    return _wrap_replicates(
        sim,
        random_seed=random_seed,
        replicate_index=replicate_index,
        num_replicates=num_replicates,
        provenance_dict=provenance_dict,
    )


class Simulator(_msprime.Simulator):
    """
    Class to simulate trees under a variety of population models.

    Note: this class is not intended to be instantiated directly
    and is only for internal library use. The interface may change
    arbitrarily between versions.
    """

    def __init__(
        self,
        *,
        tables,
        recombination_map,
        gene_conversion_map,
        gene_conversion_track_length,
        discrete_genome,
        ploidy,
        random_generator,
        demography,
        model_change_events,
        model=None,
        store_migrations=False,
        store_full_arg=False,
        start_time=None,
        end_time=None,
        num_labels=None,
    ):
        # We always need at least n segments, so no point in making
        # allocation any smaller than this.
        num_samples = len(tables.nodes)
        block_size = 64 * 1024
        segment_block_size = max(block_size, num_samples)
        avl_node_block_size = block_size
        node_mapping_block_size = block_size

        if num_labels is None:
            num_labels = self._choose_num_labels(model, model_change_events)

        # Now, convert the high-level values into their low-level
        # counterparts.
        ll_simulation_model = model.get_ll_representation()
        ll_population_configuration = [pop.asdict() for pop in demography.populations]
        ll_demographic_events = [
            event.get_ll_representation() for event in demography.events
        ]
        ll_recomb_map = recombination_map.asdict()
        ll_tables = _msprime.LightweightTableCollection(tables.sequence_length)
        ll_tables.fromdict(tables.asdict())

        # FIXME support arbitrary gene conversion maps.
        # https://github.com/tskit-dev/msprime/issues/1212
        assert len(gene_conversion_map) == 1
        gene_conversion_rate = gene_conversion_map.rate[0]

        start_time = -1 if start_time is None else start_time
        super().__init__(
            tables=ll_tables,
            recombination_map=ll_recomb_map,
            start_time=start_time,
            random_generator=random_generator,
            model=ll_simulation_model,
            migration_matrix=demography.migration_matrix,
            population_configuration=ll_population_configuration,
            demographic_events=ll_demographic_events,
            store_migrations=store_migrations,
            store_full_arg=store_full_arg,
            num_labels=num_labels,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size,
            gene_conversion_rate=gene_conversion_rate,
            gene_conversion_track_length=gene_conversion_track_length,
            discrete_genome=discrete_genome,
            ploidy=ploidy,
        )
        # highlevel attributes used externally that have no lowlevel equivalent
        self.end_time = end_time
        self.model_change_events = model_change_events
        self.demography = demography
        # Temporary, until we add the low-level infrastructure for the gc map
        # when we'll take the same approach as the recombination map.
        self.gene_conversion_map = gene_conversion_map

    def copy_tables(self):
        """
        Returns a copy of the underlying table collection. This is useful
        for testing and avoids using the LightweightTableCollection object,
        which is returned by self.tables.
        """
        return tskit.TableCollection.fromdict(self.tables.asdict())

    @property
    def sample_configuration(self):
        """
        Returns a list of the number of samples in each of the populations.
        """
        tables = self.copy_tables()
        num_samples = [0 for _ in tables.populations]
        for node in tables.nodes:
            if (node.flags & tskit.NODE_IS_SAMPLE) != 0:
                num_samples[node.population] += 1
        return num_samples

    @property
    def recombination_map(self):
        return intervals.RateMap(**super().recombination_map)

    def _choose_num_labels(self, model, model_change_events):
        """
        Choose the number of labels appropriately, given the simulation
        models that will be simulated.
        """
        num_labels = 1
        models = [model] + [event.model for event in model_change_events]
        for model in models:
            if isinstance(model, SweepGenicSelection):
                num_labels = 2
        return num_labels

    def _run_until(self, end_time, event_chunk=None, debug_func=None):
        # This is a pretty big default event chunk so that we don't spend
        # too much time going back and forth into Python. We could imagine
        # doing something a bit more sophisticated where we try to tune the
        # number of events so that we end up with roughly 10 second slices
        # (say).
        if event_chunk is None:
            event_chunk = 10 ** 4
        if event_chunk <= 0:
            raise ValueError("Must have at least 1 event per chunk")
        logger.info("Running model %s until max time: %f", self.model, end_time)
        while super().run(end_time, event_chunk) == _msprime.EXIT_MAX_EVENTS:
            logger.debug("time=%g ancestors=%d", self.time, self.num_ancestors)
            if debug_func is not None:
                debug_func(self)

    def run(self, event_chunk=None, debug_func=None):
        """
        Runs the simulation until complete coalescence has occurred.
        """
        for event in self.model_change_events:
            # If the event time is a callable, we compute the end_time
            # as a function of the current simulation time.
            current_time = self.time
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
            self._run_until(model_start_time, event_chunk, debug_func)
            if self.time > model_start_time:
                raise NotImplementedError(
                    "The previously running model does not support ending early "
                    "and the requested model change cannot be performed. Please "
                    "open an issue on GitHub if this functionality is something "
                    "you require"
                )
            ll_new_model = event.model.get_ll_representation()
            self.model = ll_new_model
        end_time = np.inf if self.end_time is None else self.end_time
        self._run_until(end_time, event_chunk, debug_func)
        self.finalise_tables()
        logger.info(
            "Completed at time=%g nodes=%d edges=%d",
            self.time,
            self.num_nodes,
            self.num_edges,
        )

    def run_replicates(
        self, num_replicates, *, mutation_rate=None, provenance_dict=None,
    ):
        """
        Sequentially yield the specified number of simulation replicates.
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
            self.run()

            if mutation_rate is not None:
                mutations._simple_mutate(
                    self.tables,
                    self.random_generator,
                    sequence_length=self.sequence_length,
                    rate=mutation_rate,
                    discrete_sites=self.discrete_genome,
                )
            tables = tskit.TableCollection.fromdict(self.tables.asdict())
            replicate_provenance = None
            if encoded_provenance is not None:
                replicate_provenance = encoded_provenance.replace(
                    f'"{placeholder}"', str(j)
                )
                tables.provenances.add_row(replicate_provenance)
            yield tables.tree_sequence()
            self.reset()


# TODO update the documentation here to state that using this class is
# deprecated, and users should use the model=[...] notation instead.
@attr.s
class SimulationModelChange:
    """
    An event representing a change of underlying :ref:`simulation model
    <sec_api_simulation_models>`.

    :param float time: The time at which the simulation model changes
        to the new model, in generations. After this time, all internal
        tree nodes, edges and migrations are the result of the new model.
        If time is set to None (the default), the model change will occur
        immediately after the previous model has completed. If time is a
        callable, the time at which the simulation model changes is the result
        of calling this function with the time that the previous model
        started with as a parameter.
    :param model: The new simulation model to use.
        This can either be a string (e.g., ``"smc_prime"``) or an instance of
        a simulation model class (e.g, ``msprime.DiscreteTimeWrightFisher(100)``.
        Please see the :ref:`sec_api_simulation_models` section for more details
        on specifying simulations models. If the argument is a string, the
        reference population size is set from the top level ``Ne`` parameter
        to :func:`.simulate`. If this is None (the default) the model is
        changed to the standard coalescent with a reference_size of
        Ne (if model was not specified).
    :type model: str or simulation model instance
    """

    time = attr.ib(default=None)
    model = attr.ib(default=None)

    def asdict(self):
        return attr.asdict(self)


@attr.s
class SimulationModel:
    """
    Abstract superclass of all simulation models.
    """

    name = None

    def get_ll_representation(self):
        return {"name": self.name}

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
    # TODO Since the pedigree is a necessary parameter for this simulation
    # model and it cannot be used with any other model we should make it a
    # parametric model where the parameter is the pedigree. This would
    # streamline a bunch of logic.
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
    A Lambda-coalescent with multiple mergers in the haploid cases, or a
    Xi-coalescent with simultaneous multiple mergers in the polyploid case.

    There are two main differences between the Beta-coalescent and the
    standard coalescent. Firstly, the number of lineages that take part in each
    common ancestor event is random, with distribution determined by moments of
    the :math:`Beta(2 - \\alpha, \\alpha)`-distribution. In particular, when there
    are :math:`n` lineages, each set of :math:`k \\leq n` of them participates in a
    common ancestor event at rate

    .. math::
        \\frac{1}{B(2 - \\alpha, \\alpha)}
        \\int_0^1 x^{k - \\alpha - 1} (1 - x)^{n - k + \\alpha - 1} dx,

    where :math:`B(2 - \\alpha, \\alpha)` is the Beta-function.

    If ploidy = 1, then all participating lineages merge into one common ancestor,
    corresponding to haploid, single-parent reproduction.
    If ploidy = :math:`p > 1`, all participating lineages split randomly into
    :math:`2 p` groups, corresponding to two-parent reproduction with :math:`p` copies
    of each chromosome per parent. All lineages within each group merge simultaneously.

    Secondly, the number of generations between common ancestor events predicted by the
    Beta-coalescent is proportional to :math:`N_e^{\\alpha - 1}`, where :math:`N_e` is
    the effective population size. Specifically, the mean number of generations until
    two lineages undergo a common ancestor event is

    .. math::
        G = \\frac{m^{\\alpha} N_e^{\\alpha - 1}}{\\alpha B(2 - \\alpha, \\alpha)},

    if ploidy = 1, and

    .. math::
        G = \\frac{m^{\\alpha} (N_e / 2)^{\\alpha - 1}}
            {2 p \\alpha B(2 - \\alpha, \\alpha)},

    if ploidy = :math:`p > 1`, where :math:`m` is the mean number of juveniles per
    family, and is given by

    .. math::
        m = q + \\frac{2^{\\alpha}}{3^{\\alpha - 1} (\\alpha - 1)},

    where :math:`q = 1` ploidy = 1, and :math:`q = 2` if ploidy > 1.

    In the polyploid case we divide the effective population size :math:`N_e` by two
    because we assume the :math:`N_e` polyploid individuals form :math:`N_e / 2`
    two-parent families in which reproduction takes place.

    .. warning::
        The number of generations between common ancestor events :math:`G` depends
        both on the effective population size :math:`N_e` and :math:`\\alpha`,
        and can be dramatically shorter than in the case of the
        standard coalescent. For :math:`\\alpha \\approx 1` that is due to
        insensitivity of :math:`G` to :math:`N_e` --- see
        :ref:`sec_api_simulation_models_multiple_mergers` for an illustration.
        For :math:`\\alpha \\approx 2`, :math:`G` is almost linear in
        :math:`N_e`, but can nevertheless be small because
        :math:`B(2 - \\alpha, \\alpha) \\rightarrow \\infty` as
        :math:`\\alpha \\rightarrow 2`. As a result, effective population sizes
        must often be many orders of magnitude larger than census population sizes
        to obtain realistic amounts of diversity in simulated samples.

    See `Schweinsberg (2003)
    <https://www.sciencedirect.com/science/article/pii/S0304414903000280>`_
    for the derivation of the common ancestor event rate,
    as well as the number of generations between common ancestor events.
    Note however that Schweinsberg (2003) only covers the haploid case.
    For details of the diploid extension, see
    `Blath et al. (2013) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3527250/>`_.
    The general polyploid model is analogous to the diploid case, with
    :math:`2 p` available copies of parental chromsomes per common ancestor event,
    and hence up to :math:`2 p` simultaneous mergers.

    :param float alpha: Determines the degree of skewness in the family size
        distribution, and must satisfy :math:`1 < \\alpha < 2`. Smaller values of
        :math:`\\alpha` correspond to greater skewness, and :math:`\\alpha = 2`
        would coincide with the standard coalescent.
    :param float truncation_point: Determines the maximum fraction of the
        population replaced by offspring in one reproduction event, and must
        satisfy :math:`0 < \\tau \\leq 1`, where :math:`\\tau` is the truncation point.
        The default is :math:`\\tau = 1`, which corresponds to the standard
        Beta-coalescent. When :math:`\\tau < 1`, the number of lineages
        participating in a common ancestor event is determined by moments
        of the Beta:math:`(2 - \\alpha, \\alpha)` distribution conditioned on not
        exceeding :math:`\\tau`, and the Beta-function in the expression
        for :math:`G` is also replaced by the incomplete Beta-function
        :math:`B(\\tau; 2 - \\alpha, \\alpha)`.
    """

    name = "beta"

    alpha = attr.ib(default=None)
    truncation_point = attr.ib(default=1)


@attr.s
class DiracCoalescent(ParametricSimulationModel):
    """
    A Lambda-coalescent with multiple mergers in the haploid cases, or a
    Xi-coalescent with simultaneous multiple mergers in the polyploid case.

    The Dirac-coalescent is an implementation of the model of
    `Blath et al. (2013) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3527250/>`_
    The simulation proceeds similarly to the standard coalescent.
    In addition to binary common ancestor events at rate :math:`n (n - 1) / 2` when
    there are :math:`n` lineages, potential multiple merger events take place
    at rate :math:`c > 0`. Each lineage participates in each multiple merger
    event independently with probability :math:`0 < \\psi \\leq 1`.

    If ploidy = 1, then all participating lineages merge into one common ancestor,
    corresponding to haploid, single-parent reproduction.
    If ploidy = :math:`p > 1`, all participating lineages split randomly into
    :math:`2 p` groups, corresponding to two-parent reproduction with :math:`p` copies
    of each chromosome per parent. All lineages within each group merge simultaneously.

    .. warning::
        The Dirac-coalescent is obtained as a scaling limit of Moran models,
        rather than Wright-Fisher models. As a consequence, the number of generations
        between coalescence events is proportional to :math:`N_e^2`,
        rather than :math:`N_e` generations as in the standard coalescent.
        See :ref:`sec_tutorial_multiple_mergers` for an illustration of how this
        affects simulation output in practice.

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
