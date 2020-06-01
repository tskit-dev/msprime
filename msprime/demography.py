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
Module responsible for defining and debugging demographic models.
"""
import collections
import inspect
import json
import logging
import math
import sys
import warnings

import attr
import numpy as np

from . import ancestry


logger = logging.getLogger(__name__)


@attr.s
class PopulationConfiguration:
    """
    The initial configuration of a population (or deme) in a simulation.

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

    sample_size = attr.ib(default=None)
    initial_size = attr.ib(default=None)
    growth_rate = attr.ib(default=0.0)
    metadata = attr.ib(default=None)

    def encode_metadata(self):
        encoded_metadata = b""
        if self.metadata is not None:
            encoded_metadata = json.dumps(self.metadata).encode()
        return encoded_metadata

    def asdict(self):
        return attr.asdict(self)


@attr.s
class DemographicEvent:
    """
    Superclass of demographic events that occur during simulations.
    """

    time = attr.ib()

    def asdict(self):
        return {
            key: getattr(self, key)
            for key in inspect.signature(self.__init__).parameters.keys()
            if hasattr(self, key)
        }


@attr.s
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

    initial_size = attr.ib(default=None)
    growth_rate = attr.ib(default=None)
    # TODO change the default to -1 to match MigrationRateChange.
    population = attr.ib(default=None)
    # Deprecated.
    # TODO add a formal deprecation notice
    population_id = attr.ib(default=None, repr=False)

    def __attrs_post_init__(self):

        if self.population_id is not None and self.population is not None:
            raise ValueError(
                "population_id and population are aliases; cannot supply both."
            )

        if self.population_id is not None:
            self.population = self.population_id
        if self.growth_rate is None and self.initial_size is None:
            raise ValueError("Must specify one or more of growth_rate and initial_size")
        if self.initial_size is not None and self.initial_size <= 0:
            raise ValueError("Cannot have a population size <= 0")
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

    def __str__(self):
        s = f"Population parameter change for {self.population}: "
        if self.initial_size is not None:
            s += f"initial_size -> {self.initial_size} "
        if self.growth_rate is not None:
            s += f"growth_rate -> {self.growth_rate} "
        return s[:-1]


@attr.s
class MigrationRateChange(DemographicEvent):
    """
    Changes the rate of migration from one deme to another to a new value at a
    specific time. Migration rates are specified in terms of the rate at which
    lineages move from population ``source`` to ``dest`` during the progress of
    the simulation. Note that ``source`` and ``dest`` are from the perspective
    of the coalescent process; please see the :ref:`sec_api_simulation_model`
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

    rate = attr.ib()
    source = attr.ib(default=-1)
    dest = attr.ib(default=-1)
    # Deprecated.
    # TODO add a formal deprecation notice
    matrix_index = attr.ib(default=None, repr=False)

    def __attrs_post_init__(self):
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

    def __str__(self):
        if self.source == -1 and self.dest == -1:
            ret = f"Migration rate change to {self.rate} everywhere"
        else:
            matrix_index = (self.source, self.dest)
            ret = f"Migration rate change for {matrix_index} to {self.rate}"
        return ret


@attr.s
class MassMigration(DemographicEvent):
    """
    A mass migration event in which some fraction of the population in one deme
    (the ``source``) simultaneously move to another deme (``dest``) during the
    progress of the simulation. Each lineage currently present in the source
    population moves to the destination population with probability equal to
    ``proportion``. Note that ``source`` and ``dest`` are from the perspective
    of the coalescent process; please see the :ref:`sec_api_simulation_model`
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

    source = attr.ib()
    # dest only has a default because of the deprecated destination attr.
    dest = attr.ib(default=None)
    proportion = attr.ib(default=1.0)
    # Deprecated.
    # TODO add a formal deprecation notice
    destination = attr.ib(default=None, repr=False)

    def __attrs_post_init__(self):
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

    def __str__(self):
        return (
            "Mass migration: "
            "Lineages moved with probability {} backwards in time with "
            "source {} & dest {}"
            "\n                     "
            "(equivalent to migration from {} to {} forwards in time)".format(
                self.proportion, self.source, self.dest, self.dest, self.source
            )
        )


# This is an unsupported/undocumented demographic event.
@attr.s
class SimpleBottleneck(DemographicEvent):
    population = attr.ib()
    proportion = attr.ib(default=1.0)

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "simple_bottleneck",
            "time": self.time,
            "population": self.population,
            "proportion": self.proportion,
        }

    def __str__(self):
        return (
            "Simple bottleneck: lineages in population {} coalesce "
            "probability {}".format(self.population, self.proportion)
        )


# TODO document
@attr.s
class InstantaneousBottleneck(DemographicEvent):
    population = attr.ib(default=None)
    strength = attr.ib(default=1.0)

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "instantaneous_bottleneck",
            "time": self.time,
            "population": self.population,
            "strength": self.strength,
        }

    def __str__(self):
        return (
            "Instantaneous bottleneck in population {}: equivalent to {} "
            "generations of the coalescent".format(self.population, self.strength)
        )


@attr.s
class CensusEvent(DemographicEvent):
    """
    An event that adds a node to each branch of every tree at a given time
    during the simulation. This may be used to record all ancestral haplotypes
    present at that time, and to extract other information related to these
    haplotypes: for instance to trace the local ancestry of a sample back to a
    set of contemporaneous ancestors, or to assess whether a subset of samples
    has coalesced more recently than the census time.
    See the :ref:`tutorial<sec_tutorial_demography_census>` for an example.

    :param float time: The time at which this event occurs in generations.
    """

    def get_ll_representation(self, num_populations=None):
        # We need to keep the num_populations argument until stdpopsim 0.1 is out
        # https://github.com/tskit-dev/msprime/issues/1037
        return {
            "type": "census_event",
            "time": self.time,
        }

    def __str__(self):
        return "Census event"


# Implementation note: these are treated as demographic events for the
# sake of the high-level interface, but are treated differently at run
# time. There is no corresponding demographic event in the C layer, as
# this would add too much complexity to the main loops. Instead, we
# detect these events at the high level, and insert calls to set_model
# as appropriate.
@attr.s
class SimulationModelChange(DemographicEvent):
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

    # We redefine time because we want it to default to None here.
    time = attr.ib(default=None)
    model = attr.ib(default=None)

    def __str__(self):
        return f"Population model changes to {self.model}"


@attr.s
class PopulationParameters:
    """
    Simple class to represent the state of a population in terms of its
    demographic parameters.
    """

    start_size = attr.ib(default=None)
    end_size = attr.ib(default=None)
    growth_rate = attr.ib(default=None)


@attr.s
class Epoch:
    """
    Represents a single epoch in the simulation within which the state
    of the demographic parameters are constant.
    """

    start_time = attr.ib(default=None)
    end_time = attr.ib(default=None)
    populations = attr.ib(default=None)
    migration_matrix = attr.ib(default=None)
    demographic_events = attr.ib(default=None)


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
        Ne=1,
        population_configurations=None,
        migration_matrix=None,
        demographic_events=None,
        model="hudson",
    ):
        if demographic_events is None:
            demographic_events = []
        self.demographic_events = demographic_events
        self._precision = 3
        # Make sure that we have a sample size of at least 2 so that we can
        # initialise the simulator.
        sample_size = None
        if population_configurations is None:
            sample_size = 2
        else:
            saved_sample_sizes = [
                pop_config.sample_size for pop_config in population_configurations
            ]
            for pop_config in population_configurations:
                pop_config.sample_size = 2
        simulator = ancestry.simulator_factory(
            sample_size=sample_size,
            model=model,
            Ne=Ne,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events,
        )
        if len(simulator.model_change_events) > 0:
            raise ValueError(
                "Model changes not currently supported by the DemographyDebugger. "
                "Please open an issue on GitHub if this feature would be useful to you"
            )
        assert len(simulator.model_change_events) == 0
        self._make_epochs(simulator, sorted(demographic_events, key=lambda e: e.time))
        self.simulation_model = simulator.model

        if population_configurations is not None:
            # Restore the saved sample sizes.
            for pop_config, sample_size in zip(
                population_configurations, saved_sample_sizes
            ):
                pop_config.sample_size = sample_size
        self.migration_matrix = migration_matrix
        self.population_configurations = population_configurations
        self._check_misspecification()

    def _make_epochs(self, simulator, demographic_events):
        self.epochs = []
        self.num_populations = simulator.num_populations
        ll_sim = simulator.create_ll_instance()
        N = simulator.num_populations
        start_time = 0
        end_time = 0
        abs_tol = 1e-9
        event_index = 0
        while not math.isinf(end_time):
            events = []
            while event_index < len(demographic_events) and math.isclose(
                demographic_events[event_index].time, start_time, abs_tol=abs_tol
            ):
                events.append(demographic_events[event_index])
                event_index += 1
            end_time = ll_sim.debug_demography()
            migration_matrix = ll_sim.get_migration_matrix()
            growth_rates = [
                conf["growth_rate"] for conf in ll_sim.get_population_configuration()
            ]
            populations = [
                PopulationParameters(
                    start_size=ll_sim.compute_population_size(j, start_time),
                    end_size=ll_sim.compute_population_size(j, end_time),
                    growth_rate=growth_rates[j],
                )
                for j in range(N)
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

    def _print_populations(self, epoch, output):
        field_width = self._precision + 6
        growth_rate_field_width = 14
        sep_str = " | "
        N = len(epoch.migration_matrix)
        fmt = (
            "{id:<2} "
            "{start_size:^{field_width}}"
            "{end_size:^{field_width}}"
            "{growth_rate:>{growth_rate_field_width}}"
        )
        print(
            fmt.format(
                id="",
                start_size="start",
                end_size="end",
                growth_rate="growth_rate",
                field_width=field_width,
                growth_rate_field_width=growth_rate_field_width,
            ),
            end=sep_str,
            file=output,
        )
        for k in range(N):
            print("{0:^{1}}".format(k, field_width), end="", file=output)
        print(file=output)
        h = "-" * (field_width - 1)
        print(
            fmt.format(
                id="",
                start_size=h,
                end_size=h,
                growth_rate=h,
                field_width=field_width,
                growth_rate_field_width=growth_rate_field_width,
            ),
            end=sep_str,
            file=output,
        )
        for _ in range(N):
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
                precision=self._precision,
                field_width=field_width,
                growth_rate_field_width=growth_rate_field_width,
            )
            print(s, end=sep_str, file=output)
            for k in range(N):
                x = epoch.migration_matrix[j][k]
                print(
                    "{0:^{1}.{2}g}".format(x, field_width, self._precision),
                    end="",
                    file=output,
                )
            print(file=output)

    def print_history(self, output=sys.stdout):
        """
        Prints a summary of the history of the populations.
        """
        print("Model = ", self.simulation_model, file=output)
        for epoch in self.epochs:
            if len(epoch.demographic_events) > 0:
                print(f"Events @ generation {epoch.start_time}", file=output)
            for event in epoch.demographic_events:
                print("   -", event, file=output)
            s = f"Epoch: {epoch.start_time} -- {epoch.end_time} generations"
            print("=" * len(s), file=output)
            print(s, file=output)
            print("=" * len(s), file=output)
            self._print_populations(epoch, output)
            print(file=output)

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
        for demo in self.demographic_events:
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
            r_prediction_close = np.allclose(r, rd[::2], rtol=1e-3)
            p_prediction_close = np.allclose(p_t, p_td[::2], rtol=1e-3)
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
        for demo in self.demographic_events:
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
            p_t[j] = np.sum(P)
            r[j] = np.sum(np.matmul(P, C)) / np.sum(P)
            G = (np.kron(M - dM, Identity) + np.kron(Identity, M - dM)) - C
            P = np.matmul(P, _matrix_exponential(dt * G))
        p_t[num_steps - 1] = np.sum(P)
        r[num_steps - 1] = np.sum(np.matmul(P, C)) / np.sum(P)
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
