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
Module responsible for defining and running ancestry simulations.
"""
from __future__ import annotations

import collections.abc
import copy
import dataclasses
import enum
import json
import logging
import math
import sys
import tempfile
import warnings
from typing import Any
from typing import ClassVar

import numpy as np
import tskit

from . import core
from . import demography as demog
from . import intervals
from . import mutations
from . import provenance
from msprime import _msprime

logger: logging.Logger = logging.getLogger(__name__)

TIME_UNITS_GENERATIONS = "generations"


class TimeUnitsMismatchWarning(UserWarning):
    """
    Warning raise when the time units specified in different parts of a
    simulation do not match.
    """


def _model_factory(model: None | str | AncestryModel) -> AncestryModel:
    """
    Returns an AncestryModel corresponding to the specified model
    description.
    - If model is None, the default simulation model is returned.
    - If model is a string, return the corresponding model instance.
    - If model is an instance of AncestryModel, return it unchanged.
    - Otherwise raise a type error.
    """
    model_map = {
        "hudson": StandardCoalescent(),
        "smc": SmcApproxCoalescent(),
        "smc_prime": SmcPrimeApproxCoalescent(),
        "dtwf": DiscreteTimeWrightFisher(),
        "fixed_pedigree": FixedPedigree(),
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
    elif not isinstance(model, AncestryModel):
        raise TypeError(
            "Ancestry model must be a string or an instance of AncestryModel"
        )
    else:
        model_instance = model
    return model_instance


def _parse_model_arg(model_arg):
    """
    Parses the specified model argument from the sim_ancestry function,
    returning the list of models.
    """
    # TODO be more lenient about what we accept as input here. Ideally
    # we'd like to support generators and consume them during the
    # actual simulation, but that would raise complications for how
    # to deal with replication.
    if isinstance(model_arg, (list, tuple)):
        if len(model_arg) < 1:
            raise ValueError("Must specify at least one AncestryModel")
        models = [_model_factory(model_desc) for model_desc in model_arg]
    else:
        models = [_model_factory(model_arg)]
    return models


def _resolve_models(events):
    model_change_events = []
    for event in events:
        assert isinstance(event, SimulationModelChange)
        # We don't want to modify our inputs, so take a deep copy.
        event = copy.copy(event)
        event.model = _model_factory(event.model)
        model_change_events.append(event)
    return model_change_events


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
    model_change_events = _resolve_models(model_change_events)
    return filtered_events, model_change_events


def _check_population_configurations(population_configurations):
    err = (
        "Population configurations must be a list of PopulationConfiguration instances"
    )
    for config in population_configurations:
        if not isinstance(config, demog.PopulationConfiguration):
            raise TypeError(err)


# This class is only used in the 0.x interface.
Sample = collections.namedtuple("Sample", ["population", "time"])


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
                    the_samples += [Sample(j, 0) for _ in range(conf.sample_size)]
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
    Ne, population_configurations, migration_matrix, demographic_events
):
    demography = demog.Demography.from_old_style(
        population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        Ne=Ne,
        ignore_sample_size=True,
    )
    return demography.validate()


def _build_initial_tables(*, sequence_length, samples, ploidy, demography):
    # NOTE: this is only used in the simulate() codepath.
    tables = tskit.TableCollection(sequence_length)
    tables.time_units = TIME_UNITS_GENERATIONS

    for index, (population, time) in enumerate(samples):
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE,
            time=time,
            population=population,
        )
        if population < 0:
            raise ValueError(f"Negative population ID in sample at index {index}")
        if population >= demography.num_populations:
            raise ValueError(
                f"Invalid population reference '{population}' in sample "
                f"at index {index}"
            )

    # This is for the simulate() code path so we don't add metadata schemas
    # and insert the user metadata in directly as encoded JSON, as before.
    for population in demography.populations:
        encoded_metadata = b""
        if population.extra_metadata is not None:
            encoded_metadata = json.dumps(population.extra_metadata).encode()
        tables.populations.add_row(encoded_metadata)

    return tables


def _parse_simulate(
    sample_size=None,
    *,
    Ne=1,
    length=None,
    recombination_rate=None,
    recombination_map=None,
    population_configurations=None,
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
    random_seed=None,
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
        raise ValueError(
            "Either sample_size, samples, population_configurations or from_ts must "
            "be specified"
        )
    samples = _samples_factory(sample_size, samples, population_configurations)

    models = [_model_factory(model)]
    if demographic_events is not None:
        demographic_events, model_change_events = _filter_events(demographic_events)
        current_time = 0 if start_time is None else start_time
        for mce in model_change_events:
            models[-1].duration = mce.time - current_time
            models.append(mce.model)
            current_time = mce.time
    is_dtwf = isinstance(models[0], DiscreteTimeWrightFisher)
    is_pedigree = any(isinstance(model, FixedPedigree) for model in models)

    demography = _demography_factory(
        Ne, population_configurations, migration_matrix, demographic_events
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
        recombination_map = intervals.RateMap.uniform(the_length, the_rate)
    else:
        if isinstance(recombination_map, intervals.RecombinationMap):
            if recombination_map._is_discrete:
                logger.info("Emulating v0.x discrete sites simulation")
                discrete_genome = True
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

    if from_ts is None:
        tables = _build_initial_tables(
            sequence_length=recombination_map.sequence_length,
            samples=samples,
            ploidy=2,
            demography=demography,
        )
    else:
        tables = from_ts.dump_tables()

    if record_full_arg:
        if is_dtwf or is_pedigree:
            raise ValueError(
                "Full ARG recording not supported in FixedPedigree and DTWF simulation"
            )
        coalescing_segments_only = False
        additional_nodes = NodeType(
            NodeType.RECOMBINANT
            | NodeType.COMMON_ANCESTOR
            | NodeType.MIGRANT
            | NodeType.GENE_CONVERSION
        )
    else:
        coalescing_segments_only = True
        additional_nodes = None

    # It's useful to call _parse_simulate outside the context of the main
    # entry point - so we want to get good seeds in this case too.
    random_seed = _parse_random_seed(random_seed)
    random_generator = _msprime.RandomGenerator(random_seed)

    sim = Simulator(
        tables=tables,
        recombination_map=recombination_map,
        models=models,
        store_migrations=record_migrations,
        additional_nodes=additional_nodes,
        coalescing_segments_only=coalescing_segments_only,
        start_time=start_time,
        end_time=end_time,
        num_labels=num_labels,
        demography=demography,
        # Defaults for the values that are not supported through simulate()
        gene_conversion_map=intervals.RateMap.uniform(
            recombination_map.sequence_length, 0
        ),
        gene_conversion_tract_length=0,
        discrete_genome=discrete_genome,
        ploidy=2,
        random_generator=random_generator,
    )
    return sim


def _parse_random_seed(seed):
    """
    Parse the specified random seed value. If no seed is provided, generate a
    high-quality random seed.
    """
    if seed is None:
        seed = core.get_random_seed()
    if isinstance(seed, np.ndarray):
        seed = seed[0]
    seed = int(seed)
    return seed


def _parse_replicate_index(*, replicate_index, random_seed, num_replicates):
    """
    Parse the replicate_index value, and ensure that its value makes sense
    in the context of the other parameters.
    """
    if replicate_index is None:
        return None
    if random_seed is None:
        raise ValueError("Cannot specify the replicate_index without a random_seed")
    if num_replicates is not None:
        raise ValueError("Cannot specify the replicate_index as well as num_replicates")
    replicate_index = int(replicate_index)
    if replicate_index < 0:
        raise ValueError("Cannot specify negative replicate_index.")
    return replicate_index


def simulate(
    sample_size=None,
    *,
    Ne=1,
    length=None,
    recombination_rate=None,
    recombination_map=None,
    mutation_rate=None,
    population_configurations=None,
    migration_matrix=None,
    demographic_events=None,
    samples=None,
    model=None,
    record_migrations=False,
    random_seed=None,
    replicate_index=None,
    mutation_generator=None,
    num_replicates=None,
    from_ts=None,
    start_time=None,
    end_time=None,
    record_full_arg=False,
    num_labels=None,
    record_provenance=True,
):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns the resulting :class:`tskit.TreeSequence`. Note that
    Ne is the effective diploid population size (so the effective number
    of genomes in the population is 2*Ne), but ``sample_size`` is the
    number of (monoploid) genomes sampled.

    .. important::
        This function is deprecated (but supported indefinitely);
        please use :func:`.sim_ancestry` in new code.

    :param int sample_size: The number of sampled monoploid genomes.  If not
        specified or None, this defaults to the sum of the subpopulation sample
        sizes. Either ``sample_size``, ``population_configurations`` or
        ``samples`` must be specified.
    :param float Ne: The effective (diploid) population size. This defaults to
        1 if not specified.
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
    :type population_configurations: list or None
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
        is drawn in the specified population at the specified time. Time
        is measured in generations ago, as elsewhere.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
    :param int num_replicates: The number of replicates of the specified
        parameters to simulate. If this is not specified or None,
        no replication is performed and a :class:`tskit.TreeSequence` object
        returned. If `num_replicates` is provided, the specified
        number of replicates is performed, and an iterator over the
        resulting :class:`tskit.TreeSequence` objects returned.
    :param tskit.TreeSequence from_ts: If specified, initialise the simulation
        from the root segments of this tree sequence and return the
        updated tree sequence. Please see :ref:`here
        <sec_ancestry_initial_state>` for details on the required properties
        of this tree sequence and its interactions with other parameters.
        (Default: None).
    :param float start_time: If specified, set the initial time that the
        simulation starts to this value. If not specified, the start
        time is zero if performing a simulation of a set of samples,
        or is the time of the oldest node if simulating from an
        existing tree sequence (see the ``from_ts`` parameter).
    :param float end_time: If specified, terminate the simulation at the
        specified time. In the returned tree sequence, all rootward paths from
        samples with time <= end_time will end in a node with one child with
        time equal to end_time. Sample nodes with time > end_time will
        also be present in the output tree sequence. If not specified or ``None``,
        run the simulation until all samples have an MRCA at all positions in
        the genome.
    :param bool record_full_arg: If True, record all intermediate nodes
        arising from common ancestor and recombination events in the output
        tree sequence. This will result in unary nodes (i.e., nodes in marginal
        trees that have only one child). Defaults to False.
    :param model: The simulation model to use.
        This can either be a string (e.g., ``"smc_prime"``) or an instance of
        a ancestry model class (e.g, ``msprime.DiscreteTimeWrightFisher()``.
    :type model: str or AncestryModel
    :param bool record_provenance: If True, record all input parameters
        in the tree sequence :ref:`tskit:sec_provenance`.
    :return: The :class:`tskit.TreeSequence` object representing the results
        of the simulation if no replication is performed, or an
        iterator over the independent replicates simulated if the
        `num_replicates` parameter has been used.
    :rtype: :class:`tskit.TreeSequence` or an iterator over
        :class:`tskit.TreeSequence` replicates.
    """
    replicate_index = _parse_replicate_index(
        random_seed=random_seed,
        num_replicates=num_replicates,
        replicate_index=replicate_index,
    )
    random_seed = _parse_random_seed(random_seed)
    provenance_dict = None
    if record_provenance:
        parameters = dict(
            command="simulate",
            sample_size=sample_size,
            Ne=Ne,
            length=length,
            recombination_rate=recombination_rate,
            recombination_map=recombination_map,
            mutation_rate=mutation_rate,
            population_configurations=population_configurations,
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
            random_seed=random_seed,
            # num_replicates is excluded as provenance is per replicate
            # replicate index is excluded as it is inserted for each replicate
        )
        provenance_dict = provenance.get_provenance_dict(parameters)

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
        Ne=Ne,
        length=length,
        recombination_rate=recombination_rate,
        recombination_map=recombination_map,
        population_configurations=population_configurations,
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
        random_seed=random_seed,
    )
    return _wrap_replicates(
        sim,
        num_replicates=num_replicates,
        replicate_index=replicate_index,
        provenance_dict=provenance_dict,
        mutation_rate=mutation_rate,
    )


def _wrap_replicates(
    simulator,
    *,
    num_replicates,
    replicate_index,
    provenance_dict,
    mutation_rate=None,
):
    """
    Wrapper for the logic used to run replicate simulations for the two
    frontends.
    """
    if num_replicates is None and replicate_index is None:
        # Default single-replicate case.
        replicate_index = 0
    if replicate_index is not None:
        num_replicates = replicate_index + 1
    iterator = simulator.run_replicates(
        num_replicates,
        mutation_rate=mutation_rate,
        provenance_dict=provenance_dict,
    )
    if replicate_index is not None:
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


def _insert_sample_sets(sample_sets, demography, default_ploidy, tables):
    """
    Insert the samples described in the specified {population_id: num_samples}
    map into the specified set of tables.
    """
    for sample_set in sample_sets:
        n = sample_set.num_samples
        population = demography[sample_set.population]
        time = (
            population.default_sampling_time
            if sample_set.time is None
            else sample_set.time
        )
        ploidy = default_ploidy if sample_set.ploidy is None else sample_set.ploidy
        logger.info(
            f"Sampling {n} individuals with ploidy {ploidy} in population "
            f"{population.id} (name='{population.name}') at time {time}"
        )
        node_individual = len(tables.individuals) + np.repeat(
            np.arange(n, dtype=np.int32), ploidy
        )
        ind_flags = np.zeros(n, dtype=np.uint32)
        tables.individuals.append_columns(flags=ind_flags)
        N = n * ploidy
        tables.nodes.append_columns(
            flags=np.full(N, tskit.NODE_IS_SAMPLE, dtype=np.uint32),
            time=np.full(N, time),
            population=np.full(N, population.id, dtype=np.int32),
            individual=node_individual,
        )


def _parse_sample_sets(sample_sets, demography):
    # Don't modify the inputs.
    sample_sets = copy.deepcopy(sample_sets)
    for sample_set in sample_sets:
        if not isinstance(sample_set, SampleSet):
            raise TypeError("msprime.SampleSet object required")
        if not core.isinteger(sample_set.num_samples):
            raise TypeError(
                "The number of samples to draw from a population must be an integer"
            )
        sample_set.num_samples = int(sample_set.num_samples)
        if sample_set.num_samples < 0:
            raise ValueError("Number of samples cannot be negative")
        if sample_set.population is None:
            if demography.num_populations == 1:
                sample_set.population = 0
            else:
                raise ValueError(
                    "Must specify a SampleSet population in multipopulation models"
                )

    if sum(sample_set.num_samples for sample_set in sample_sets) == 0:
        raise ValueError("Zero samples specified")
    return sample_sets


def _parse_samples(samples, demography):
    """
    Parse the specified "samples" value for sim_ancestry and return the equivalent
    list of SampleSet instances.
    """
    if isinstance(samples, collections.abc.Sequence):
        sample_sets = samples
    elif isinstance(samples, collections.abc.Mapping):
        sample_sets = [
            SampleSet(num_samples, population)
            for population, num_samples in samples.items()
        ]
    elif core.isinteger(samples):
        if demography.num_populations != 1:
            raise ValueError(
                "Numeric samples can only be used in single population models. "
                "Please use either the mapping form {pop_id: num_samples} or "
                "a list of SampleSet instances."
            )
        sample_sets = [SampleSet(samples)]
    else:
        raise TypeError(
            f"The value '{samples}' cannot be interpreted as sample specification. "
            "Samples must either be a single integer, a dict that maps populations "
            "to the number of samples for that population, or a list of SampleSet "
            "objects. Please see the online documentation for more details on "
            "the different forms."
        )
    return _parse_sample_sets(sample_sets, demography)


def _parse_sim_ancestry(
    samples=None,
    *,
    sequence_length=None,
    recombination_rate=None,
    gene_conversion_rate=None,
    gene_conversion_tract_length=None,
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
    additional_nodes=None,
    coalescing_segments_only=None,
    num_labels=None,
    random_seed=None,
    init_for_debugger=False,
):
    """
    Argument parser for the sim_ancestry frontend. Interprets all the parameters
    and returns an appropriate instance of Simulator.
    """
    # As a general rule we try to cast any input value to the required types
    # early and in a way that provides an interpretable traceback.

    # Simple defaults.
    start_time = 0 if start_time is None else float(start_time)
    end_time = math.inf if end_time is None else float(end_time)
    discrete_genome = core._parse_flag(discrete_genome, default=True)
    record_full_arg = core._parse_flag(record_full_arg, default=False)

    models = _parse_model_arg(model)
    is_dtwf = isinstance(models[0], DiscreteTimeWrightFisher)
    is_pedigree = any(isinstance(model, FixedPedigree) for model in models)

    if record_full_arg:
        if coalescing_segments_only is not None:
            raise ValueError(
                "Cannot set both coalescing_segments_only and record_full_arg."
            )
        if additional_nodes is not None:
            raise ValueError("Cannot set both record_full_arg and additional_nodes.")
        coalescing_segments_only = False
        if is_dtwf or is_pedigree:
            raise ValueError(
                "Full ARG recording not supported in FixedPedigree and DTWF simulation"
            )
        else:
            additional_nodes = NodeType(
                NodeType.RECOMBINANT
                | NodeType.COMMON_ANCESTOR
                | NodeType.MIGRANT
                | NodeType.GENE_CONVERSION
            )
    coalescing_segments_only = core._parse_flag(coalescing_segments_only, default=True)
    if additional_nodes is None:
        additional_nodes = NodeType(0)
    else:
        if not isinstance(additional_nodes, NodeType):
            raise ValueError(
                "additional_nodes should be an instance of msprime.NodeType."
            )
        if additional_nodes.value > 0:
            if coalescing_segments_only:
                raise ValueError(
                    "Set coalescing_segments_only to False when recording"
                    " additional_nodes."
                )

        if additional_nodes.value & NodeType.MIGRANT.value:
            if is_dtwf:
                raise ValueError(
                    "Recording MIGRANT nodes is currently not supported in"
                    " DTWF simulation."
                )
            if is_pedigree:
                raise ValueError(
                    "Recording MIGRANT nodes is not supported in "
                    "FixedPedigree simulation."
                )
            if any(isinstance(model, SweepGenicSelection) for model in models):
                raise ValueError(
                    "Recording MIGRANT nodes is not supported in "
                    "SweepGenicSelection simulation."
                )

    record_migrations = core._parse_flag(record_migrations, default=False)

    if initial_state is not None:
        if isinstance(initial_state, tskit.TreeSequence):
            initial_state = initial_state.dump_tables()
        elif not isinstance(initial_state, tskit.TableCollection):
            raise TypeError(
                "initial_state must either be a TreeSequence or TableCollection instance"
            )
        if initial_state.time_units == tskit.TIME_UNITS_UNCALIBRATED:
            raise ValueError(
                "Cannot use a tree sequence with uncalibrated time_units as "
                "the initial state"
            )
        if initial_state.time_units != TIME_UNITS_GENERATIONS:
            message = (
                f"The initial_state has time_units={initial_state.time_units} but "
                "time is measured in generations in msprime. This may lead to "
                "significant discrepancies between the timescales. "
                "If you wish to suppress this warning, you can use, e.g., "
                "warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)"
            )
            warnings.warn(message, TimeUnitsMismatchWarning, stacklevel=2)

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
    if gene_conversion_tract_length is None:
        if gene_conversion_rate is None:
            # It doesn't matter what the tract_length is, just set a
            # value to keep the low-level code happy.
            gene_conversion_tract_length = 1
        else:
            raise ValueError(
                "Must specify tract length when simulating gene conversion"
            )
    else:
        if gene_conversion_rate is None:
            raise ValueError(
                "Must specify gene conversion rate along with tract length"
            )
        gene_conversion_tract_length = float(gene_conversion_tract_length)

    # Default to diploid
    ploidy = 2 if ploidy is None else ploidy
    if not core.isinteger(ploidy):
        raise TypeError("ploidy must be an integer")
    ploidy = int(ploidy)
    if ploidy < 1:
        raise ValueError("ploidy must be >= 1")

    if is_pedigree:
        if demography is not None:
            raise ValueError("Cannot specify demography for FixedPedigree simulation")
        if population_size is not None:
            raise ValueError(
                "Cannot specify population_size for FixedPedigree simulation"
            )
        if initial_state is None:
            raise ValueError(
                "Must specify an input pedigree using the ``initial_state`` argument "
                "for a FixedPedigree simulation"
            )
        if ploidy != 2:
            raise ValueError("Fixed pedigree simulations must have ploidy=2")
        if gene_conversion_map.total_mass != 0:
            raise ValueError(
                "Gene conversion not supported in FixedPedigree simulation"
            )
        if record_migrations:
            raise ValueError(
                "Migration recording not supported in FixedPedigree simulation"
            )
        if start_time > 0:
            raise ValueError("Cannot specify start_time in FixedPedigree simulation")
        if len(models) > 1:
            raise ValueError(
                "Cannot use FixedPedigree simulation in conjunction with other models"
            )
        demography = demog.Demography.from_tree_sequence(
            initial_state.tree_sequence(), initial_size=1
        )

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
        if initial_state is not None:
            if population_size is None:
                raise ValueError(
                    "Must specify either a demography object or a population_size "
                    "(for single population models) when providing an initial_state."
                )
            if len(initial_state.populations) > 1:
                raise ValueError(
                    "Must specify demography for initial_state with > 1 population"
                )
        population_size = 1 if population_size is None else float(population_size)
        demography = demog.Demography.isolated_model([population_size])
    elif isinstance(demography, demog.Demography):
        if population_size is not None:
            raise ValueError("Cannot specify demography and population size")
    else:
        raise TypeError("demography argument must be an instance of msprime.Demography")
    demography = demography.validate()

    if initial_state is None:
        if samples is None and not init_for_debugger:
            raise ValueError(
                "Either the samples or initial_state arguments must be provided"
            )
        initial_state = tskit.TableCollection(sequence_length)
        initial_state.time_units = TIME_UNITS_GENERATIONS
        demography.insert_populations(initial_state)
        if not init_for_debugger:
            sample_sets = _parse_samples(samples, demography)
            _insert_sample_sets(sample_sets, demography, ploidy, initial_state)
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
        # Make sure the names match-up in the input demography.
        demography_check = demog.Demography.from_tree_sequence(
            initial_state.tree_sequence()
        )
        if demography.num_populations < demography_check.num_populations:
            raise ValueError(
                "Input demography must have at least as many populations as the "
                "initial state population table: "
                f"{demography.num_populations} < {demography_check.num_populations}"
            )

        for pop1, pop2 in zip(demography.populations, demography_check.populations):
            if pop1.name != pop2.name:
                raise ValueError(
                    "Population names in the input demography and the initial "
                    f"state population table must be equal: {pop1.name} ≠ {pop2.name}"
                )
        demography.insert_extra_populations(initial_state)

    # It's useful to call _parse_sim_ancestry outside the context of the main
    # entry point - so we want to get good seeds in this case too.
    random_seed = _parse_random_seed(random_seed)
    random_generator = _msprime.RandomGenerator(random_seed)

    return Simulator(
        tables=initial_state,
        recombination_map=recombination_map,
        gene_conversion_map=gene_conversion_map,
        gene_conversion_tract_length=gene_conversion_tract_length,
        discrete_genome=discrete_genome,
        ploidy=ploidy,
        demography=demography,
        models=models,
        store_migrations=record_migrations,
        additional_nodes=additional_nodes,
        coalescing_segments_only=coalescing_segments_only,
        start_time=start_time,
        end_time=end_time,
        num_labels=num_labels,
        random_generator=random_generator,
    )


def sim_ancestry(
    samples=None,
    *,
    demography=None,
    sequence_length=None,
    discrete_genome=None,
    recombination_rate=None,
    gene_conversion_rate=None,
    gene_conversion_tract_length=None,
    population_size=None,
    ploidy=None,
    model=None,
    initial_state=None,
    start_time=None,
    end_time=None,
    record_migrations=None,
    record_full_arg=None,
    additional_nodes=None,
    coalescing_segments_only=None,
    num_labels=None,
    random_seed=None,
    num_replicates=None,
    replicate_index=None,
    record_provenance=None,
):
    """
    Simulates an ancestral process described by the specified model, demography and
    samples, and return a :class:`tskit.TreeSequence` (or a sequence of
    replicate tree sequences).

    :param samples: The sampled individuals as either an integer, specifying
        the number of individuals to sample in a single-population model;
        or a list of :class:`.SampleSet` objects defining the properties of
        groups of similar samples; or as a mapping in which the keys
        are population identifiers (either an integer ID or string name)
        and the values are the number of samples to take from the corresponding
        population at its default sampling time. It is important to note that
        samples correspond to *individuals* here, and each sampled individual
        is usually associated with :math:`k` sample *nodes* (or genomes) when
        ``ploidy`` = :math:`k`. See the :ref:`sec_ancestry_samples` section
        for further details.
        Either ``samples`` or ``initial_state`` must be specified.
    :param demography: The demographic model to simulate, describing the
        extant and ancestral populations, their population sizes and growth
        rates, their migration rates, and demographic events affecting the
        populations over time. See the :ref:`sec_demography` section for
        details on how to specify demographic models and
        :ref:`sec_ancestry_samples` for details on how to specify the
        populations that samples are drawn from. If not specified (or None) we
        default to a single population with constant size 1
        (see also the ``population_size`` parameter).
    :param int ploidy: The number of monoploid genomes per sample individual
        (Default=2). See the :ref:`sec_ancestry_ploidy` section for usage examples.
    :param float sequence_length: The length of the genome sequence to simulate.
        See the :ref:`sec_ancestry_sequence_length` section for usage examples
        for this parameter and how it interacts with other parameters.
    :param bool discrete_genome: If True (the default) simulation occurs
        in discrete genome coordinates such that recombination and
        gene conversion breakpoints always occur at integer positions.
        Thus, multiple (e.g.) recombinations can occur at the same
        genome position. If ``discrete_genome`` is False simulations
        are performed using continuous genome coordinates. In this
        case multiple events at precisely the same genome location are very
        unlikely (but technically possible).
        See the :ref:`sec_ancestry_discrete_genome` section for usage examples.
    :param recombination_rate: The rate of recombination along the sequence;
        can be either a single value (specifying a single rate over the entire
        sequence) or an instance of :class:`RateMap`.
        See the :ref:`sec_ancestry_recombination` section for usage examples
        for this parameter and how it interacts with other parameters.
    :param gene_conversion_rate: The rate of gene conversion along the sequence.
        If provided, a value for ``gene_conversion_tract_length`` must also be
        specified. See the :ref:`sec_ancestry_gene_conversion` section
        for usage examples for this parameter and how it interacts with
        other parameters.
    :param gene_conversion_tract_length: The mean length of the gene conversion
        tracts. For discrete genomes the tract lengths are geometrically
        distributed with mean ``gene_conversion_tract_length``, which must be
        greater than or equal to 1. For continuous genomes the tract lengths are
        exponentially distributed with mean ``gene_conversion_tract_length``,
        which must be larger than 0.
    :param population_size: The number of individuals of the default single population
        :class:`.Demography`. If not specified, defaults to 1. Cannot be specified
        along with the ``demography`` parameter. See the :ref:`sec_demography`
        section for more details on demographic models and population sizes
        and the :ref:`sec_ancestry_demography` section for usage examples.
    :param int random_seed: The random seed. If this is not specified or `None`,
        a high-quality random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
        See the :ref:`sec_randomness_seeds` section for usage examples.
    :param int num_replicates: The number of replicates of the specified
        parameters to simulate. If this is not specified or `None`,
        no replication is performed and a :class:`tskit.TreeSequence` object
        returned. If `num_replicates` is provided, the specified
        number of replicates is performed, and an iterator over the
        resulting :class:`tskit.TreeSequence` objects returned.
        See the :ref:`sec_randomness_replication` section for examples.
    :param bool record_full_arg: If True, record all intermediate nodes
        arising from common ancestor and recombination events in the output
        tree sequence. This will result in unary nodes (i.e., nodes in marginal
        trees that have only one child). Defaults to False.
        See the :ref:`sec_ancestry_full_arg` section for examples.
    :param NodeType additional_nodes: Retain all ancestry for any node
        of the specified type. This will result in unary nodes.
        Defaults to class `msprime.NodeType(0)`.
    :param bool coalescing_segments_only: If False, retain all ancestry for any
        nodes that are coalescent nodes anywhere in the sequence. This will result in
        unary nodes. Should be set to False when :ref:`recording additional nodes
        <sec_ancestry_additional_nodes>`. Defaults to True.
    :param bool record_migrations: If True, record all migration events
        that occur in the :ref:`tskit:sec_migration_table_definition` of
        the output tree sequence. Defaults to False.
        See the :ref:`sec_ancestry_record_migrations` section for examples.
    :param tskit.TreeSequence initial_state: If specified, initialise the
        simulation from the root segments of this tree sequence and return the
        completed tree sequence. Please see
        :ref:`sec_ancestry_initial_state` for details of the required
        properties of this tree sequence and its interactions with other parameters.
        All information in the ``initial_state`` tables is preserved
        (including metadata) and included in the returned tree sequence.
        (Default: None).
    :param float start_time: If specified, set the initial time that the
        simulation starts to this value. If not specified, the start
        time is zero if performing a simulation of a set of samples,
        or is the time of the oldest node if simulating from an
        existing tree sequence (see the ``initial_state`` parameter).
        See the :ref:`sec_ancestry_start_time` section for examples.
    :param float end_time: If specified, terminate the simulation at the
        specified time. In the returned tree sequence, all rootward paths from
        samples with time < ``end_time`` will end in a node with one child with
        time equal to end_time. Any sample nodes with time >= ``end_time`` will
        also be present in the output tree sequence. If not specified or ``None``,
        run the simulation until all samples have an MRCA at all positions in
        the genome. See the :ref:`sec_ancestry_end_time` section for examples.
    :param bool record_provenance: If True (the default), record all input
        parameters in the tree sequence :ref:`tskit:sec_provenance`.
    :param model: The ancestry model to use. This can be either a
        single instance of :class:`.AncestryModel` (or a string that can be
        interpreted as an ancestry model), or a list of :class:`.AncestryModel`
        instances. If the ``duration`` attribute of any of these models is
        set, the simulation will be run until at most :math:`t + t_m`, where
        :math:`t` is the simulation time when the model starts and :math:`t_m`
        is the model's ``duration``. If the ``duration`` is not set, the
        simulation will continue until the model completes, the overall
        ``end_time`` is reached, or overall coalescence. See
        the :ref:`sec_ancestry_models_specifying` section for more details,
        and the :ref:`sec_ancestry_models` section for the available models
        and examples.
    :type model: str or msprime.AncestryModel or list
    :return: The :class:`tskit.TreeSequence` object representing the results
        of the simulation if no replication is performed, or an
        iterator over the independent replicates simulated if the
        `num_replicates` parameter has been used.
    :rtype: :class:`tskit.TreeSequence` or an iterator over
        :class:`tskit.TreeSequence` replicates.
    """
    record_provenance = True if record_provenance is None else record_provenance
    replicate_index = _parse_replicate_index(
        random_seed=random_seed,
        num_replicates=num_replicates,
        replicate_index=replicate_index,
    )
    random_seed = _parse_random_seed(random_seed)
    provenance_dict = None
    if record_provenance:
        parameters = dict(
            command="sim_ancestry",
            samples=samples,
            demography=demography,
            sequence_length=sequence_length,
            discrete_genome=discrete_genome,
            recombination_rate=recombination_rate,
            gene_conversion_rate=gene_conversion_rate,
            gene_conversion_tract_length=gene_conversion_tract_length,
            population_size=population_size,
            ploidy=ploidy,
            model=model,
            initial_state=initial_state,
            start_time=start_time,
            end_time=end_time,
            record_migrations=record_migrations,
            record_full_arg=record_full_arg,
            additional_nodes=additional_nodes,
            coalescing_segments_only=coalescing_segments_only,
            num_labels=num_labels,
            random_seed=random_seed,
            # num_replicates is excluded as provenance is per replicate
            # replicate index is excluded as it is inserted for each replicate
        )
        provenance_dict = provenance.get_provenance_dict(parameters)
    sim = _parse_sim_ancestry(
        samples=samples,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        gene_conversion_rate=gene_conversion_rate,
        gene_conversion_tract_length=gene_conversion_tract_length,
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
        additional_nodes=additional_nodes,
        coalescing_segments_only=coalescing_segments_only,
        num_labels=num_labels,
        random_seed=random_seed,
    )
    return _wrap_replicates(
        sim,
        num_replicates=num_replicates,
        replicate_index=replicate_index,
        provenance_dict=provenance_dict,
    )


class ExitReason(enum.IntEnum):
    """
    The different reasons that the low-level simulation exits.
    """

    MAX_EVENTS = _msprime.EXIT_MAX_EVENTS
    """
    We ran for the specified maximum number of events. We usually
    run for a maximum number of events so that we return to Python
    regularly, to update logs and to check if CTRL-C has been hit, etc.
    """

    MAX_TIME = _msprime.EXIT_MAX_TIME
    """
    The simulation is equal to the specified maximum time.
    """

    COALESCENCE = _msprime.EXIT_COALESCENCE
    """
    The simulation is complete, and we have fully coalesced.
    """

    MODEL_COMPLETE = _msprime.EXIT_MODEL_COMPLETE
    """
    The model we have specified has run to completion **without**
    resulting in coalescence.
    """


class NodeType(enum.Flag):
    """
    Specify the type of node for which you want to track ancestry. Extends
    the :class:`python:enum.Flag` class. The NodeType can be specified by
    means of bitwise operators on the members of the enumeration.
    See :ref:`additional nodes <sec_ancestry_additional_nodes>` for a
    definition of each of the different node types.
    """

    RECOMBINANT = 1 << 17
    COMMON_ANCESTOR = 1 << 18
    MIGRANT = 1 << 19
    CENSUS = 1 << 20
    GENE_CONVERSION = 1 << 21
    PASS_THROUGH = 1 << 22

    def asdict(self):
        return {"__class__": "msprime.NodeType", "value": self.value}


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
        gene_conversion_tract_length,
        discrete_genome,
        ploidy,
        demography,
        random_generator,
        models=None,
        store_migrations=False,
        additional_nodes=None,
        coalescing_segments_only=True,
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
            num_labels = self._choose_num_labels(models)
        if additional_nodes is None:
            additional_nodes = NodeType(0)

        # Now, convert the high-level values into their low-level
        # counterparts.
        ll_population_configuration = [pop.asdict() for pop in demography.populations]
        ll_demographic_events = [
            event.get_ll_representation() for event in demography.events
        ]
        ll_recomb_map = self._resolve_missing_intervals(recombination_map)

        ll_tables = _msprime.LightweightTableCollection(tables.sequence_length)
        ll_tables.fromdict(tables.asdict())

        # FIXME support arbitrary gene conversion maps.
        # https://github.com/tskit-dev/msprime/issues/1212
        assert len(gene_conversion_map.rate) == 1
        gene_conversion_rate = gene_conversion_map.rate[0]

        start_time = -1 if start_time is None else start_time
        super().__init__(
            tables=ll_tables,
            recombination_map=ll_recomb_map,
            start_time=start_time,
            random_generator=random_generator,
            migration_matrix=demography.migration_matrix,
            population_configuration=ll_population_configuration,
            demographic_events=ll_demographic_events,
            store_migrations=store_migrations,
            additional_nodes=additional_nodes.value,
            coalescing_segments_only=coalescing_segments_only,
            num_labels=num_labels,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size,
            gene_conversion_rate=gene_conversion_rate,
            gene_conversion_tract_length=gene_conversion_tract_length,
            discrete_genome=discrete_genome,
            ploidy=ploidy,
        )
        # Highlevel attributes used externally that have no lowlevel equivalent
        self.end_time = np.inf if end_time is None else end_time
        self.models = models
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

    def _choose_num_labels(self, models):
        """
        Choose the number of labels appropriately, given the ancestry
        models that will be simulated.
        """
        num_labels = 1
        for model in models:
            if isinstance(model, SweepGenicSelection):
                num_labels = 2
        return num_labels

    def _resolve_missing_intervals(self, recombination_map):
        """
        Inspect the recombination map for unknown intervals, resolve
        the appropriate recombination rate to use in the actual simulation,
        and return the low-level recombination map representation.

        Also store the set of missing_intervals for later use in delete_intervals.

        For now we only support unknown values in the flanking regions,
        and insist that we can have at most two of them (so, no attempt to
        coalesce adjacent unknown intervals at the ends of the map).
        See https://github.com/tskit-dev/msprime/issues/1604 for plans
        and discussion for how to simulate chromosomes with unknown
        regions in the middle.
        """
        self.missing_intervals = recombination_map.missing_intervals()
        error_msg = (
            "Missing regions of the genome other than the flanks are currently "
            "not supported. Please see "
            "https://github.com/tskit-dev/msprime/issues/1604"
        )
        if self.missing_intervals.shape[0] > 2:
            raise ValueError(error_msg)
        for left, right in self.missing_intervals:
            if not (left == 0 or right == recombination_map.sequence_length):
                raise ValueError(error_msg)

        ll_recomb_map = recombination_map.asdict()
        missing = recombination_map.missing
        rate = recombination_map.rate.copy()
        rate[missing] = 0
        ll_recomb_map["rate"] = rate
        return ll_recomb_map

    def _run_until(self, end_time, event_chunk=None, debug_func=None):
        # This is a pretty big default event chunk so that we don't spend
        # too much time going back and forth into Python. We could imagine
        # doing something a bit more sophisticated where we try to tune the
        # number of events so that we end up with roughly 10 second slices
        # (say).
        if event_chunk is None:
            event_chunk = 10**4
        if event_chunk <= 0:
            raise ValueError("Must have at least 1 event per chunk")
        logger.info("Running model %s until max time: %f", self.model, end_time)
        ret = ExitReason.MAX_EVENTS
        while ret == ExitReason.MAX_EVENTS:
            ret = ExitReason(super().run(end_time, event_chunk))
            if self.time > end_time:
                # Currently the Pedigree and Sweeps models are "non-reentrant"
                # We can change this to an assertion once these have been fixed.
                raise RuntimeError(
                    f"Model {self.model['name']} does not support interruption. "
                    "Please open an issue on GitHub"
                )
            logger.debug(
                "time=%g ancestors=%d ret=%s", self.time, self.num_ancestors, ret
            )
            if debug_func is not None:
                debug_func(self)
        return ret

    def run(self, event_chunk=None, debug_func=None):
        """
        Runs the simulation until complete coalescence has occurred,
        end_time has been reached, or all model durations have
        elapsed.
        """
        for j, model in enumerate(self.models):
            self.model = model._as_lowlevel()
            logger.info(
                "model[%d] %s started at time=%g nodes=%d edges=%d",
                j,
                self.model,
                self.time,
                self.num_nodes,
                self.num_edges,
            )
            model_duration = np.inf if model.duration is None else model.duration
            if model_duration < 0:
                raise ValueError("Model durations must be >= 0")
            end_time = min(self.time + model_duration, self.end_time)
            exit_reason = self._run_until(end_time, event_chunk, debug_func)
            if exit_reason == ExitReason.COALESCENCE or self.time == self.end_time:
                logger.debug("Skipping remaining %d models", len(self.models) - j - 1)
                break
        self.finalise_tables()
        logger.info(
            "Completed at time=%g nodes=%d edges=%d",
            self.time,
            self.num_nodes,
            self.num_edges,
        )

    def run_replicates(
        self,
        num_replicates,
        *,
        mutation_rate=None,
        provenance_dict=None,
    ):
        """
        Sequentially yield the specified number of simulation replicates.
        """
        encoded_provenance = None
        # The JSON is modified for each replicate to insert the replicate number.
        # To avoid repeatedly encoding the same JSON (which can take milliseconds)
        # we insert a replaceable string.
        placeholder = "@@_REPLICATE_INDEX_@@"
        if provenance_dict is not None:
            provenance_dict["parameters"]["replicate_index"] = placeholder
            encoded_provenance = provenance.json_encode_provenance(
                provenance_dict, num_replicates
            )

        for replicate_index in range(num_replicates):
            logger.info("Starting replicate %d", replicate_index)
            self.run()
            if mutation_rate is not None:
                # This is only called from simulate() or the ms interface,
                # so does not need any further parameters.
                mutations._simple_mutate(
                    tables=self.tables,
                    random_generator=self.random_generator,
                    sequence_length=self.sequence_length,
                    rate=mutation_rate,
                )
            tables = tskit.TableCollection.fromdict(self.tables.asdict())
            if len(self.missing_intervals) > 0:
                tables.delete_intervals(
                    self.missing_intervals, simplify=False, record_provenance=False
                )
            replicate_provenance = None
            if encoded_provenance is not None:
                replicate_provenance = encoded_provenance.replace(
                    f'"{placeholder}"', str(replicate_index)
                )
                tables.provenances.add_row(replicate_provenance)
            ts = tables.tree_sequence()
            yield ts
            self.reset()

    def __str__(self):
        # Warning! This can be very big as it's a direct dump of the low-level
        # data structures. If you want to debug a large simulation use
        # simulator.print_state(sys.stdout).
        # It's also possible to trigger segfault/assertions by calling this on a
        # simulator that's in an errored state - so use this method with care!
        with tempfile.TemporaryFile() as f:
            self.print_state(f)
            f.seek(0)
            s = f.read().decode()
        return "================\nSimulator state:\n================\n" + s


@dataclasses.dataclass
class SampleSet:
    """
    Specify a set of exchangeable sample individuals with a given ploidy
    value from a population at a given time. See the
    :ref:`sec_ancestry_samples` section for details and examples.
    """

    num_samples: int
    """
    The number of k-ploid sample **individuals** to draw.
    """
    population: int | str | None = None
    """
    The population in which the samples are drawn. May be either a
    string name or integer ID (see
    :ref:`sec_demography_populations_identifiers` details).
    """
    time: float | None = None
    """
    The time at which these samples are drawn. If not specified or None,
    defaults to the :attr:`.Population.default_sampling_time`.
    """
    ploidy: int | None = None
    """
    The number of monoploid genomes to sample for each sample individual.
    See the :ref:`sec_ancestry_ploidy` section for more details and
    examples.
    """

    def asdict(self):
        return dataclasses.asdict(self)


@dataclasses.dataclass
class SimulationModelChange:
    """
    Demographic event denoting an change in ancestry model.

    .. important::
        This class is deprecated (but supported indefinitely);
        please use the ``model`` argument in :func:`sim_ancestry`
        to specify multiple models in new code.
    """

    time: float | None = None
    """
    The time at which the ancestry model changes to the new model, in
    generations. After this time, all internal tree nodes, edges and migrations
    are the result of the new model. If time is set to None (the default), the
    model change will occur immediately after the previous model has completed.
    """

    # Can't use typehints here because having a reference to AncestryModel
    # breaks autodoc, which wants to call it msprime.ancestry.AncestryModel
    # whereas we have it documented as msprime.AncestryModel. Annoying.
    model: Any = None
    """
    The new ancestry model to use. This can either be a string (e.g.,
    ``"smc_prime"``) or an instance of an ancestry model class (e.g,
    ``msprime.DiscreteTimeWrightFisher()``. Please see the
    :ref:`sec_ancestry_models` section for more details on specifying these
    models. If this is None (the default) the model is changed to the standard
    coalescent.
    """

    def asdict(self):
        return dataclasses.asdict(self)


@dataclasses.dataclass(init=False)
class AncestryModel:
    """
    Abstract superclass of all ancestry models.
    """

    duration: float | None
    """
    The time duration that this model should run for. If None, the model
    will run until completion (i.e., until the simulation coalesces
    or the model itself completes). Otherwise, this defines the maximum
    time duration which the model can run. See the
    :ref:`sec_ancestry_models_specifying` section for more details.
    """
    name: ClassVar[str]

    # We have to define an __init__ to enforce keyword-only behaviour
    def __init__(self, *, duration=None):
        self.duration = duration

    # We need to have a separate _as_lowlevel and asdict because the
    # asdict form can't have the name in the dictionary for the
    # provenance code.
    def _as_lowlevel(self):
        d = {"name": self.name}
        d.update(self.asdict())
        return d

    def asdict(self):
        return dataclasses.asdict(self)


class StandardCoalescent(AncestryModel):
    """
    The classical coalescent with recombination model (i.e., Hudson's algorithm).
    The string ``"hudson"`` can be used to refer to this model.

    This is a continuous time model in which the time to the next event
    is exponentially distributed with rates depending on the population size(s),
    migration rates, numbers of extant lineages and the amount of ancestral
    material currently present. See
    `Kelleher et al. (2016) <https://doi.org/10.1371/journal.pcbi.1004842>`_ for a
    detailed description of the model and further references.
    """

    name = "hudson"


class SmcApproxCoalescent(AncestryModel):
    """
    The Sequentially Markov Coalescent (SMC) model defined by
    `McVean and Cardin (2005) <https://dx.doi.org/10.1098%2Frstb.2005.1673>`_.
    In the SMC, only common ancestor events that result in marginal coalescences
    are possible. Under this approximation, the marginal trees along the
    genome depend only on the immediately previous tree (i.e. are Markovian).

    .. note::
        This model is implemented using a naive rejection sampling approach
        and so it may not be any more efficient to simulate than the
        standard Hudson model.

    The string ``"smc"`` can be used to refer to this model.
    """

    name = "smc"


class SmcPrimeApproxCoalescent(AncestryModel):
    """
    The SMC' model defined by
    `Marjoram and Wall (2006) <https://doi.org/10.1186/1471-2156-7-16>`_
    as a refinement of the :class:`SMC<SmcApproxCoalescent>`. The SMC'
    extends the SMC by additionally allowing common ancestor events that
    join contiguous tracts of ancestral material (as well as events that
    result in marginal coalescences).

    .. note::
        This model is implemented using a naive rejection sampling approach
        and so it may not be any more efficient to simulate than the
        standard Hudson model.

    The string ``"smc_prime"`` can be used to refer to this model.
    """

    name = "smc_prime"


class DiscreteTimeWrightFisher(AncestryModel):
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


class FixedPedigree(AncestryModel):
    """
    Backwards-time simulations through a pre-specified pedigree, with diploid
    individuals and back-and-forth recombination. The string ``"fixed_pedigree"`` can
    be used to refer to this model.

    The model has no parameters, but requires that the
    :ref:`initial_state<sec_ancestry_initial_state>` parameter is provided to
    :func:`.sim_ancestry` and that it contains a valid pedigree. See the
    :ref:`sec_pedigrees` section for more information on the required
    encoding and methods for translating existing pedigrees into this format.

    Unlike other :ref:`ancestry models<sec_ancestry_models>`, the
    :class:`FixedPedigree` cannot be combined with other simulation models, and
    special attention must be paid to the role of founders when
    :ref:`recapitating these
    simulations<sec_ancestry_models_fixed_pedigree_completing>`.

    See the :ref:`sec_ancestry_models_fixed_pedigree` section for more information
    for more information and examples of running simulations within a
    fixed pedigree.
    """

    name = "fixed_pedigree"


class ParametricAncestryModel(AncestryModel):
    """
    The superclass of ancestry models that require extra parameters.
    """


@dataclasses.dataclass
class BetaCoalescent(ParametricAncestryModel):
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
    Beta-coalescent is proportional to :math:`N^{\\alpha - 1}`, where :math:`N` is
    the population size. Specifically, the mean number of generations until
    two lineages undergo a common ancestor event is

    .. math::
        G = \\frac{m^{\\alpha} N^{\\alpha - 1}}{\\alpha B(2 - \\alpha, \\alpha)},

    if ploidy = 1, and

    .. math::
        G = \\frac{2 p m^{\\alpha} (N / 2)^{\\alpha - 1}}
            {\\alpha B(2 - \\alpha, \\alpha)},

    if ploidy = :math:`p > 1`, where :math:`m` is the mean number of juveniles per
    family given by

    .. math::
        m = 2 + \\frac{2^{\\alpha}}{3^{\\alpha - 1} (\\alpha - 1)},

    if ploidy > 1, and

    .. math::
        m = 1 + \\frac{1}{2^{\\alpha - 1} (\\alpha - 1)},

    if ploidy = 1.

    In the polyploid case we divide the population size :math:`N` by two
    because we assume the :math:`N` polyploid individuals form :math:`N / 2`
    two-parent families in which reproduction takes place.

    .. warning::
        The number of generations between common ancestor events :math:`G` depends
        both on the population size :math:`N` and :math:`\\alpha`,
        and can be dramatically shorter than in the case of the
        standard coalescent. For :math:`\\alpha \\approx 1` that is due to
        insensitivity of :math:`G` to :math:`N` --- see
        :ref:`sec_ancestry_models_multiple_mergers` for an illustration.
        For :math:`\\alpha \\approx 2`, :math:`G` is almost linear in
        :math:`N`, but can nevertheless be small because
        :math:`B(2 - \\alpha, \\alpha) \\rightarrow \\infty` as
        :math:`\\alpha \\rightarrow 2`. As a result, population sizes
        must often be many orders of magnitude larger than census population sizes
        to obtain realistic amounts of diversity in simulated samples.

    See `Schweinsberg (2003)
    <https://www.sciencedirect.com/science/article/pii/S0304414903000280>`_
    for the derivation of the common ancestor event rate,
    as well as the number of generations between common ancestor events.
    Note however that Schweinsberg (2003) only covers the haploid case.
    For details of the diploid extension, see
    `Blath et al. (2013) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3527250/>`_,
    and `Birkner et al. (2018) <https://projecteuclid.org/euclid.ejp/1527818427>`_
    for a diploid version of the Schweinsberg (2003) model specifically.
    The general polyploid model is analogous to the diploid case, with
    :math:`2 p` available copies of parental chromosomes per common ancestor event,
    and hence up to :math:`2 p` simultaneous mergers.

    :param float alpha: Determines the degree of skewness in the family size
        distribution, and must satisfy :math:`1 < \\alpha < 2`. Smaller values of
        :math:`\\alpha` correspond to greater skewness, and :math:`\\alpha = 2`
        would coincide with the standard coalescent.
    :param float truncation_point: The maximum number of juveniles :math:`K` born to
        one family as a fraction of the population size :math:`N`. Must satisfy
        :math:`0 < K \\leq \\inf`. Determines the maximum fraction of the population
        replaced by offspring in one reproduction event, :math:`\\tau`, via
        :math:`\\tau = K / (K + m)`, where :math:`m` is the mean juvenile number
        above. The default is :math:`K = \\inf`, which corresponds to the standard
        Beta-coalescent with :math:`\\tau = 1`. When :math:`K < \\inf`, the number of
        lineages participating in a common ancestor event is determined by moments
        of the Beta:math:`(2 - \\alpha, \\alpha)` distribution conditioned on not
        exceeding :math:`\\tau`, and the Beta-function in the expression
        for :math:`G` is replaced by the incomplete Beta-function
        :math:`B(\\tau; 2 - \\alpha, \\alpha)`.
    """

    name = "beta"

    alpha: float | None
    truncation_point: float

    # We have to define an __init__ to enforce keyword-only behaviour
    def __init__(
        self, *, duration=None, alpha=None, truncation_point=sys.float_info.max
    ):
        self.duration = duration
        self.alpha = alpha
        self.truncation_point = truncation_point


@dataclasses.dataclass
class DiracCoalescent(ParametricAncestryModel):
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
        between coalescence events is proportional to :math:`N^2`,
        rather than :math:`N` generations as in the standard coalescent.
        See :ref:`sec_ancestry_models_multiple_mergers` for an illustration
        of how this affects simulation output in practice.

    :param float c: Determines the rate of potential multiple merger events.
        We require :math:`c > 0`.
    :param float psi: Determines the fraction of the population replaced by
        offspring in one large reproduction event, i.e. one reproduction event
        giving rise to potential multiple mergers when viewed backwards in time.
        We require :math:`0 < \\psi \\leq 1`.
    """

    name = "dirac"

    psi: float | None
    c: float | None

    # We have to define an __init__ to enforce keyword-only behaviour
    def __init__(self, *, duration=None, psi=None, c=None):
        self.duration = duration
        self.psi = psi
        self.c = c


@dataclasses.dataclass
class SweepGenicSelection(ParametricAncestryModel):
    """
    A selective sweep that has occurred in the history of the sample.
    This will lead to a burst of rapid coalescence near the selected site.

    The strength of selection during the sweep is determined by the
    parameter :math:`s`. Here we define s such that the
    fitness of the three genotypes at our beneficial locus are
    :math:`W_{bb}=1`, :math:`W_{Bb}=1 + s/2`, :math:`W_{BB}=1 + s`.
    Thus fitness of the heterozygote is intermediate to the
    two homozygotes.

    .. warning:: Models where `start_frequency` is
       :math:`\\geq \\frac{1}{ploidy*N}`
       will not simulate the trajectory until its origination.
       In practical terms, this issue means
       that "sweeps from neutral standing genetic variation" are currently not
       possible. The tracking issue for this on GitHub is
       `here <https://github.com/tskit-dev/msprime/issues/1762>`_.

    The model is one of a structured coalescent where selective backgrounds are
    defined as in
    `Braverman et al. (1995) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1206652/>`_
    The implementation details here follow closely those in discoal
    `(Kern and Schrider, 2016)
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5167068/>`_, with
    the important difference that :math:`s` in msprime is half of
    of that in discoal (i.e. for equivalent results  use :math:`2s`).

    See :ref:`sec_ancestry_models_selective_sweeps` for examples and
    details on how to specify different types of sweeps.

    .. warning::
        Currently models with more than one population and a selective sweep
        are not implemented. Population size changes during the sweep
        are not yet possible in msprime.

    :param float position: the location of the beneficial allele along the
        chromosome.
    :param float start_frequency: population frequency of the beneficial
        allele at the start of the selective sweep. E.g., for a *de novo*
        allele in a diploid population of size N, start frequency would be
        :math:`1/2N`.
    :param float end_frequency: population frequency of the beneficial
        allele at the end of the selective sweep.
    :param float s: :math:`s` is the selection coefficient of the beneficial mutation.
    :param float dt: dt is the small increment of time for stepping through
        the sweep phase of the model. a good rule of thumb is for this to be
        approximately :math:`1/40N` or smaller.
    """

    name = "sweep_genic_selection"

    position: float | None
    start_frequency: float | None
    end_frequency: float | None
    s: float | None
    dt: float | None

    # We have to define an __init__ to enforce keyword-only behaviour
    def __init__(
        self,
        *,
        duration=None,
        position=None,
        start_frequency=None,
        end_frequency=None,
        s=None,
        dt=None,
    ):
        self.duration = duration
        self.position = position
        self.start_frequency = start_frequency
        self.end_frequency = end_frequency
        self.s = s
        self.dt = dt
