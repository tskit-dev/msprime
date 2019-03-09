#
# Copyright (C) 2015-2018 University of Oxford
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
Command line interfaces to the msprime library.
"""
import argparse
import hashlib
import json
import os
import random
import signal
import sys

import msprime
import tskit


def set_sigpipe_handler():
    if os.name == "posix":
        # Set signal handler for SIGPIPE to quietly kill the program.
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)


#######################################################
# mspms: the ms compatible interface
#######################################################

mscompat_description = (
    "mspms is an ms-compatible interface to the msprime library. "
    "It simulates the coalescent with recombination for a variety of "
    "demographic models and outputs the results in a text-based format. "
    "It supports a subset of the functionality available in ms and aims "
    "for full compatibility.")
mscompat_recombination_help = (
    "Recombination at rate rho=4*N0*r where r is the rate of recombination "
    "between the ends of the region being simulated; num_loci is the number "
    "of sites between which recombination can occur")
msprime_citation_text = """
If you use msprime in your work, please cite the following paper:
Jerome Kelleher, Alison M Etheridge and Gilean McVean (2016), "Efficient
Coalescent Simulation and Genealogical Analysis for Large Sample Sizes",
PLoS Comput Biol 12(5): e1004842. doi: 10.1371/journal.pcbi.1004842
"""


def positive_int(value):
    int_value = int(float(value))
    if int_value <= 0:
        msg = "{} in an invalid postive integer value".format(value)
        raise argparse.ArgumentTypeError(msg)
    return int_value


def add_sample_size_argument(parser):
    parser.add_argument(
        "sample_size", type=positive_int,
        help="The number of genomes in the sample")


def add_tree_sequence_argument(parser):
    parser.add_argument(
        "tree_sequence", help="The msprime tree sequence file")


def add_precision_argument(parser):
    parser.add_argument(
        "--precision", "-p", type=int, default=6,
        help="The number of decimal places to print in records")


def generate_seeds():
    """
    Generate seeds to seed the RNG and output on the command line.
    """
    # Pull three numbers from the SystemRandom generator
    rng = random.SystemRandom()
    return [rng.randint(0, 2**31) for _ in range(3)]


def get_single_seed(seeds):
    """
    Takes the specified command line seeds convert them to a single value
    that can be used to seed the python random number generator.
    """
    assert len(seeds) == 3
    m = hashlib.md5()
    for s in seeds:
        # Colon separate the values to ensure that we don't have
        # collisions in situations like 1:23:45, 12:3:45.
        m.update("{}:".format(s).encode())
    # Now take the integer value of this modulo 2^32, as this is
    # the largest seed value we'll accept.
    return int(m.hexdigest(), 16) % (2**32)


class SimulationRunner(object):
    """
    Class to run msprime simulation and output the results.
    """
    def __init__(
            self, sample_size=1, num_loci=1, scaled_recombination_rate=0,
            num_replicates=1, migration_matrix=None,
            population_configurations=None, demographic_events=None,
            scaled_mutation_rate=0, print_trees=False,
            precision=3, random_seeds=None):
        self._sample_size = sample_size
        self._num_loci = num_loci
        self._num_replicates = num_replicates
        self._recombination_rate = scaled_recombination_rate
        self._mutation_rate = scaled_mutation_rate
        # For strict ms-compability we want to have m non-recombining loci
        recomb_map = msprime.RecombinationMap.uniform_map(
            num_loci, self._recombination_rate, num_loci)
        # If we have specified any population_configurations we don't want
        # to give the overall sample size.
        sample_size = self._sample_size
        if population_configurations is not None:
            sample_size = None
        # msprime measure's time in units of generations, given a specific
        # Ne value whereas ms uses coalescent time. To be compatible with ms,
        # we therefore need to use an Ne value of 1/4.
        self._simulator = msprime.simulator_factory(
            Ne=0.25,
            sample_size=sample_size,
            recombination_map=recomb_map,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        self._precision = precision
        self._print_trees = print_trees
        # sort out the random seeds
        ms_seeds = random_seeds
        if random_seeds is None:
            ms_seeds = generate_seeds()
        seed = get_single_seed(ms_seeds)
        self._random_generator = msprime.RandomGenerator(seed)
        self._ms_random_seeds = ms_seeds
        self._simulator.random_generator = self._random_generator
        self._mutation_generator = msprime.MutationGenerator(
            self._random_generator, self._mutation_rate)

    def get_num_replicates(self):
        """
        Returns the number of replicates we are to run.
        """
        return self._num_replicates

    def get_simulator(self):
        """
        Returns the simulator instance for this simulation runner.
        """
        return self._simulator

    def get_mutation_rate(self):
        """
        Returns the per-base, per generation mutation rate used by
        msprime.
        """
        return self._mutation_rate

    def print_trees(self, tree_sequence, output):
        """
        Print out the trees in ms-format from the specified tree sequence.
        When 'invisible' recombinations occur ms prints out copies of the
        same tree. Therefore, we must keep track of all breakpoints from the
        simulation and write out a tree for each one.
        """
        breakpoints = self._simulator.breakpoints + [self._num_loci]
        if self._num_loci == 1:
            tree = next(tree_sequence.trees())
            newick = tree.newick(precision=self._precision)
            print(newick, file=output)
        else:
            j = 0
            for tree in tree_sequence.trees():
                newick = tree.newick(precision=self._precision)
                left, right = tree.interval
                while j < len(breakpoints) and breakpoints[j] <= right:
                    length = breakpoints[j] - left
                    j += 1
                    # Print these seperately to avoid the cost of creating
                    # another string.
                    print("[{}]".format(int(length)), end="", file=output)
                    print(newick, file=output)

    def run(self, output):
        """
        Runs the simulations and writes the output to the specified
        file handle.
        """
        # The first line of ms's output is the command line.
        print(" ".join(sys.argv), file=output)
        print(" ".join(str(s) for s in self._ms_random_seeds), file=output)
        for j in range(self._num_replicates):
            self._simulator.run()
            tree_sequence = self._simulator.get_tree_sequence(self._mutation_generator)
            print(file=output)
            print("//", file=output)
            if self._print_trees:
                self.print_trees(tree_sequence, output)
            if self._mutation_rate > 0:
                s = tree_sequence.get_num_mutations()
                print("segsites:", s, file=output)
                if s != 0:
                    print("positions: ", end="", file=output)
                    positions = [
                        mutation.position / self._num_loci for mutation in
                        tree_sequence.mutations()]
                    positions.sort()
                    for position in positions:
                        print(
                            "{0:.{1}f}".format(position, self._precision),
                            end=" ", file=output)
                    print(file=output)
                    for h in tree_sequence.haplotypes():
                        print(h, file=output)
                else:
                    print(file=output)
            self._simulator.reset()


def convert_int(value, parser):
    """
    Converts the specified value to an integer if possible. If
    conversion fails, exit by calling parser.error.
    """
    try:
        return int(value)
    except ValueError:
        parser.error("invalid int value '{}'".format(value))


def convert_float(value, parser):
    """
    Converts the specified value to a float if possible. If
    conversion fails, exit by calling parser.error.
    """
    try:
        return float(value)
    except ValueError:
        parser.error("invalid float value '{}'".format(value))


def convert_population_id(parser, population_id, num_populations):
    """
    Checks the specified population ID makes sense and returns
    it as an integer, and converted into msprime's internal population
    ID scheme (i.e., zero based).
    """
    pid = int(population_id)
    if pid != population_id:
        msg = "Bad population ID '{}': must be an integer"
        parser.error(msg.format(population_id))
    if pid < 1 or pid > num_populations:
        msg = "Bad population ID '{}': must be 1 to num_populations"
        parser.error(msg.format(pid))
    return pid - 1


def check_migration_rate(parser, rate):
    """
    Checks that the specified migration rate makes sense.
    """
    if rate < 0:
        parser.error("Migration rates must be non-negative.")


def check_event_time(parser, time):
    """
    Checks that the specified time for an event makes sense
    """
    if time < 0:
        parser.error("Event times must be non-negative.")


def convert_migration_matrix(parser, input_matrix, num_populations):
    """
    Converts the specified migration matrix into the internal format.
    """
    if len(input_matrix) != num_populations**2:
        parser.error(
            "Must be num_populations^2 migration matrix entries")
    migration_matrix = [[
        0 for j in range(num_populations)] for k in range(num_populations)]
    for j in range(num_populations):
        for k in range(num_populations):
            if j != k:
                rate = convert_float(
                    input_matrix[j * num_populations + k], parser)
                check_migration_rate(parser, rate)
                migration_matrix[j][k] = rate
    return migration_matrix


def raise_admixture_incompatability_error(parser, other_option):
    """
    Because of the state dependency within ms, it is messy for us to
    support options that affect all populations in conjunction with
    population splits. For now, raise an error.
    """
    parser.error(
        "Cannot currently use the -es and {} options together. "
        "Please open an issue on GitHub if this functionality is "
        "important to you.".format(other_option))


def create_simulation_runner(parser, arg_list):
    """
    Parses the arguments and returns a SimulationRunner instance.
    """
    args = parser.parse_args(arg_list)
    if args.mutation_rate == 0 and not args.trees:
        parser.error("Need to specify at least one of --theta or --trees")
    num_loci = int(args.recombination[1])
    if args.recombination[1] != num_loci:
        parser.error("Number of loci must be integer value")
    if args.recombination[0] != 0.0 and num_loci < 2:
        parser.error("Number of loci must > 1")
    r = 0.0
    # We don't scale recombination or mutation rates by the size
    # of the region.
    if num_loci > 1:
        r = args.recombination[0] / (num_loci - 1)
    mu = args.mutation_rate / num_loci

    # Check the structure format.
    symmetric_migration_rate = 0.0
    num_populations = 1
    population_configurations = [msprime.PopulationConfiguration(args.sample_size)]
    migration_matrix = [[0.0]]
    if args.structure is not None:
        num_populations = convert_int(args.structure[0], parser)
        # We must have at least num_population sample_configurations
        if len(args.structure) < num_populations + 1:
            parser.error("Must have num_populations sample sizes")
        population_configurations = [None for j in range(num_populations)]
        for j in range(num_populations):
            population_configurations[j] = msprime.PopulationConfiguration(
                convert_int(args.structure[j + 1], parser))
        total = sum(conf.sample_size for conf in population_configurations)
        if total != args.sample_size:
            parser.error("Population sample sizes must sum to sample_size")
        # We optionally have the overall migration_rate here
        if len(args.structure) == num_populations + 2:
            symmetric_migration_rate = convert_float(
                args.structure[num_populations + 1], parser)
            check_migration_rate(parser, symmetric_migration_rate)
        elif len(args.structure) > num_populations + 2:
            parser.error("Too many arguments to --structure/-I")
        if num_populations > 1:
            migration_matrix = [[
                symmetric_migration_rate / (num_populations - 1) * int(j != k)
                for j in range(num_populations)]
                for k in range(num_populations)]
    else:
        if len(args.migration_matrix_entry) > 0:
            parser.error(
                "Cannot specify migration matrix entries without "
                "first providing a -I option")
        if args.migration_matrix is not None:
            parser.error(
                "Cannot specify a migration matrix without "
                "first providing a -I option")
    if args.migration_matrix is not None:
        migration_matrix = convert_migration_matrix(
            parser, args.migration_matrix, num_populations)
    for matrix_entry in args.migration_matrix_entry:
        dest = convert_population_id(parser, matrix_entry[0], num_populations)
        source = convert_population_id(
            parser, matrix_entry[1], num_populations)
        rate = matrix_entry[2]
        if dest == source:
            parser.error("Cannot set diagonal elements in migration matrix")
        check_migration_rate(parser, rate)
        migration_matrix[dest][source] = rate

    # Set the initial demography
    demographic_events = []
    if args.growth_rate is not None:
        for config in population_configurations:
            config.growth_rate = args.growth_rate
    for population_id, growth_rate in args.population_growth_rate:
        pid = convert_population_id(parser, population_id, num_populations)
        population_configurations[pid].growth_rate = growth_rate
    for population_id, size in args.population_size:
        pid = convert_population_id(parser, population_id, num_populations)
        population_configurations[pid].initial_size = size

    # First we look at population split events. We do this differently
    # to ms, as msprime requires a fixed number of population. Therefore,
    # modify the number of populations to take into account populations
    # splits. This is a messy hack, and will probably need to be changed.
    for index, (t, population_id, proportion) in args.admixture:
        check_event_time(parser, t)
        pid = convert_population_id(parser, population_id, num_populations)
        if proportion < 0 or proportion > 1:
            parser.error("Proportion value must be 0 <= p <= 1.")
        # In ms, the probability of staying in source is p and the probabilty
        # of moving to the new population is 1 - p.
        event = (index, msprime.MassMigration(
            t, pid, num_populations, 1 - proportion))
        demographic_events.append(event)

        num_populations += 1
        # We add another element to each row in the migration matrix
        # along with an other row. All new entries are zero.
        for row in migration_matrix:
            row.append(0)
        migration_matrix.append([0 for j in range(num_populations)])
        # Add another PopulationConfiguration object with a sample size
        # of zero.
        population_configurations.append(msprime.PopulationConfiguration(0))

    # Add the demographic events
    for index, (t, alpha) in args.growth_rate_change:
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-eG")
        check_event_time(parser, t)
        demographic_events.append(
            (index, msprime.PopulationParametersChange(
                time=t, growth_rate=alpha)))
    for index, (t, population_id, alpha) in args.population_growth_rate_change:
        pid = convert_population_id(parser, population_id, num_populations)
        check_event_time(parser, t)
        demographic_events.append(
            (index, msprime.PopulationParametersChange(
                time=t, growth_rate=alpha, population_id=pid)))
    for index, (t, x) in args.size_change:
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-eN")
        check_event_time(parser, t)
        demographic_events.append(
            (index, msprime.PopulationParametersChange(
                time=t, initial_size=x, growth_rate=0)))
    for index, (t, population_id, x) in args.population_size_change:
        check_event_time(parser, t)
        pid = convert_population_id(parser, population_id, num_populations)
        demographic_events.append(
            (index, msprime.PopulationParametersChange(
                time=t, initial_size=x, growth_rate=0, population_id=pid)))
    for index, (t, source, dest) in args.population_split:
        check_event_time(parser, t)
        source_id = convert_population_id(parser, source, num_populations)
        dest_id = convert_population_id(parser, dest, num_populations)
        demographic_events.append(
            (index, msprime.MassMigration(t, source_id, dest_id, 1.0)))
        # Set the migration rates for source to 0
        for j in range(num_populations):
            if j != source_id:
                event = msprime.MigrationRateChange(t, 0.0, (j, source_id))
                demographic_events.append((index, event))

    # Demographic events that affect the migration matrix
    if num_populations == 1:
        condition = (
            len(args.migration_rate_change) > 0 or
            len(args.migration_matrix_entry_change) > 0 or
            len(args.migration_matrix_change) > 0)
        if condition:
            parser.error("Cannot change migration rates for 1 population")
    for index, (t, x) in args.migration_rate_change:
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-eM")
        check_migration_rate(parser, x)
        check_event_time(parser, t)
        event = msprime.MigrationRateChange(
            t, x / (num_populations - 1))
        demographic_events.append((index, event))
    for index, event in args.migration_matrix_entry_change:
        t = event[0]
        check_event_time(parser, t)
        dest = convert_population_id(parser, event[1], num_populations)
        source = convert_population_id(parser, event[2], num_populations)
        if dest == source:
            parser.error("Cannot set diagonal elements in migration matrix")
        rate = event[3]
        check_migration_rate(parser, rate)
        msp_event = msprime.MigrationRateChange(t, rate, (dest, source))
        demographic_events.append((index, msp_event))
    for index, event in args.migration_matrix_change:
        if len(event) < 3:
            parser.error("Need at least three arguments to -ma")
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-ema")
        t = convert_float(event[0], parser)
        check_event_time(parser, t)
        if convert_int(event[1], parser) != num_populations:
            parser.error(
                "num_populations must be equal for new migration matrix")
        matrix = convert_migration_matrix(parser, event[2:], num_populations)
        for j in range(num_populations):
            for k in range(num_populations):
                if j != k:
                    msp_event = msprime.MigrationRateChange(
                        t, matrix[j][k], (j, k))
                    demographic_events.append((index, msp_event))

    # We've created all the events and PopulationConfiguration objects. Because
    # msprime uses absolute population sizes we need to rescale these relative
    # to Ne, since this is what ms does.
    for _, msp_event in demographic_events:
        if isinstance(msp_event, msprime.PopulationParametersChange):
            if msp_event.initial_size is not None:
                msp_event.initial_size /= 4
    for config in population_configurations:
        if config.initial_size is not None:
            config.initial_size /= 4

    demographic_events.sort(key=lambda x: (x[0], x[1].time))
    time_sorted = sorted(demographic_events, key=lambda x: x[1].time)
    if demographic_events != time_sorted:
        parser.error(
            "Demographic events must be supplied in non-decreasing "
            "time order")
    runner = SimulationRunner(
        sample_size=args.sample_size,
        num_loci=num_loci,
        migration_matrix=migration_matrix,
        population_configurations=population_configurations,
        demographic_events=[event for _, event in demographic_events],
        num_replicates=args.num_replicates,
        scaled_recombination_rate=r,
        scaled_mutation_rate=mu,
        precision=args.precision,
        print_trees=args.trees,
        random_seeds=args.random_seeds)
    return runner


class IndexedAction(argparse._AppendAction):
    """
    Argparse action class that allows us to find the overall ordering
    across several different options. We use this for the demographic
    events, as the order in which the events are applied matters for
    ms compatability.
    """
    index = 0

    def __call__(self, parser, namespace, values, option_string=None):
        super().__call__(
            parser, namespace, (IndexedAction.index, values), option_string)
        IndexedAction.index += 1


def convert_arg_line_to_args(arg_line):
    # from the docs on argparse.ArgumentParser.convert_arg_line_to_args
    return arg_line.split()


def make_load_file_action(next_parser):
    """
    Argparse action class to allow passing a filename containing arguments
    on the command line (for super-long argument files).
    From
        http://stackoverflow.com/q/27433316
    and
        http://stackoverflow.com/q/40060571
    """
    class LoadFromFile(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            try:
                with open(values) as f:
                    # note parses with 'next_parser' *not* with parser that is
                    # passed in
                    next_parser.parse_args(f.read().split(), namespace)
            except IOError as ioe:
                parser.error(ioe)
    return LoadFromFile


def get_mspms_parser(error_handler=None):
    # Ensure that the IndexedAction counter is set to zero. This is useful
    # for testing where we'll be creating lots of these parsers.
    IndexedAction.index = 0
    # to allow `-f` options we'll need a parser that can do all the
    # arguments except the positional (nonoptional) arguments.  We'll
    # create this one first, then at the end make the parser that
    # includes the positional arguments.
    parser = argparse.ArgumentParser(
        description=mscompat_description,
        epilog=msprime_citation_text)
    parser.convert_arg_line_to_args = convert_arg_line_to_args

    group = parser.add_argument_group("Behaviour")
    group.add_argument(
        "--mutation-rate", "-t", type=float, metavar="theta",
        help="Mutation rate theta=4*N0*mu", default=0)
    group.add_argument(
        "--trees", "-T", action="store_true",
        help="Print out trees in Newick format")
    group.add_argument(
        "--recombination", "-r", type=float, nargs=2, default=(0, 1),
        metavar=("rho", "num_loci"), help=mscompat_recombination_help)

    group = parser.add_argument_group("Structure and migration")
    group.add_argument(
        "--structure", "-I", nargs='+', metavar="value",
        help=(
            "Sample from populations with the specified deme structure. "
            "The arguments are of the form 'num_populations "
            "n1 n2 ... [4N0m]', specifying the number of populations, "
            "the sample configuration, and optionally, the migration "
            "rate for a symmetric island model"))
    group.add_argument(
        "--migration-matrix-entry", "-m", action="append",
        metavar=("dest", "source", "rate"),
        nargs=3, type=float, default=[],
        help=(
            "Sets an entry M[dest, source] in the migration matrix to the "
            "specified rate. source and dest are (1-indexed) population "
            "IDs. Multiple options can be specified."))
    group.add_argument(
        "--migration-matrix", "-ma", nargs='+', default=None, metavar="entry",
        help=(
            "Sets the migration matrix to the specified value. The "
            "entries are in the order M[1,1], M[1, 2], ..., M[2, 1],"
            "M[2, 2], ..., M[N, N], where N is the number of populations."))
    group.add_argument(
        "--migration-rate-change", "-eM", nargs=2, action=IndexedAction,
        type=float, default=[], metavar=("t", "x"),
        help=(
            "Set the symmetric island model migration rate to "
            "x / (npop - 1) at time t"))
    group.add_argument(
        "--migration-matrix-entry-change", "-em", action=IndexedAction,
        metavar=("time", "dest", "source", "rate"),
        nargs=4, type=float, default=[],
        help=(
            "Sets an entry M[dest, source] in the migration matrix to the "
            "specified rate at the specified time. source and dest are "
            "(1-indexed) population IDs."))
    group.add_argument(
        "--migration-matrix-change", "-ema", nargs='+', default=[],
        action=IndexedAction, metavar="entry", help=(
            "Sets the migration matrix to the specified value at time t."
            "The entries are in the order M[1,1], M[1, 2], ..., M[2, 1],"
            "M[2, 2], ..., M[N, N], where N is the number of populations."))

    group = parser.add_argument_group("Demography")
    group.add_argument(
        "--growth-rate", "-G", metavar="alpha", type=float,
        help="Set the growth rate to alpha for all populations.")
    group.add_argument(
        "--population-growth-rate", "-g", action="append", default=[],
        nargs=2, metavar=("population_id", "alpha"), type=float,
        help="Set the growth rate to alpha for a specific population.")
    group.add_argument(
        "--population-size", "-n", action="append", default=[],
        nargs=2, metavar=("population_id", "size"), type=float,
        help="Set the size of a specific population to size*N0.")

    group.add_argument(
        "--growth-rate-change", "-eG", nargs=2, action=IndexedAction,
        type=float, default=[], metavar=("t", "alpha"),
        help="Set the growth rate for all populations to alpha at time t")
    group.add_argument(
        "--population-growth-rate-change", "-eg", nargs=3,
        action=IndexedAction, type=float, default=[],
        metavar=("t", "population_id", "alpha"),
        help=(
            "Set the growth rate for a specific population to "
            "alpha at time t"))
    group.add_argument(
        "--size-change", "-eN", nargs=2, action=IndexedAction,
        type=float, default=[], metavar=("t", "x"),
        help="Set the population size for all populations to x * N0 at time t")
    group.add_argument(
        "--population-size-change", "-en", nargs=3,
        action=IndexedAction, type=float, default=[],
        metavar=("t", "population_id", "x"),
        help=(
            "Set the population size for a specific population to "
            "x * N0 at time t"))
    group.add_argument(
        "--population-split", "-ej", nargs=3,
        action=IndexedAction, type=float, default=[],
        metavar=("t", "dest", "source"),
        help=(
            "Move all lineages in population dest to source at time t. "
            "Forwards in time, this corresponds to a population split "
            "in which lineages in source split into dest. All migration "
            "rates for population source are set to zero."))
    group.add_argument(
        "--admixture", "-es", nargs=3,
        action=IndexedAction, type=float, default=[],
        metavar=("t", "population_id", "proportion"),
        help=(
            "Split the specified population into a new population, such "
            "that the specified proportion of lineages remains in "
            "the population population_id. Forwards in time this "
            "corresponds to an admixture event. The new population has ID "
            "num_populations + 1. Migration rates to and from the new "
            "population are set to 0, and growth rate is 0 and the "
            "population size for the new population is N0."))

    group = parser.add_argument_group("Miscellaneous")
    group.add_argument(
        "--random-seeds", "-seeds", nargs=3, type=positive_int,
        metavar=("x1", "x2", "x3"),
        help="Random seeds (must be three integers)")
    group.add_argument(
        "--precision", "-p", type=positive_int, default=3,
        help="Number of values after decimal place to print")

    # now for the parser that gets called first
    init_parser = argparse.ArgumentParser(
        description=mscompat_description,
        epilog=msprime_citation_text,
        add_help=False,
        parents=[parser])
    init_parser.convert_arg_line_to_args = convert_arg_line_to_args

    add_sample_size_argument(init_parser)
    init_parser.add_argument(
        "num_replicates", type=positive_int,
        help="Number of independent replicates")
    init_parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(msprime.__version__))
    init_parser.add_argument(
        "-f", "--filename", action=make_load_file_action(parser),
        help="Insert commands from a file at this point in the command line.")

    # Set the optional error handler (used for testing)
    if error_handler is not None:
        parser.error = error_handler
        init_parser.error = error_handler

    return init_parser


def get_mspms_runner(arg_list):
    parser = get_mspms_parser()
    return create_simulation_runner(parser, arg_list)


def mspms_main(arg_list=None):
    set_sigpipe_handler()
    sr = get_mspms_runner(arg_list)
    sr.run(sys.stdout)


#######################################################
# msp: the command line interface for msprime
#######################################################


def exit(message):
    sys.exit(message)


def run_upgrade(args):
    try:
        tree_sequence = tskit.load_legacy(args.source, args.remove_duplicate_positions)
    except tskit.DuplicatePositionsError:
        exit(
            "Error: Duplicate mutation positions in the source file detected.\n\n"
            "This is not supported in the current file format. Running \"upgrade -d\" "
            "will remove these duplicate positions. However, this will result in loss "
            "of data from the original file!")
    tree_sequence.dump(args.destination)


def run_dump_newick(args):
    tree_sequence = tskit.load(args.tree_sequence)
    for tree in tree_sequence.trees():
        newick = tree.newick(precision=args.precision)
        print(newick)


def run_dump_haplotypes(args):
    tree_sequence = tskit.load(args.tree_sequence)
    for h in tree_sequence.haplotypes():
        print(h)


def run_dump_variants(args):
    tree_sequence = tskit.load(args.tree_sequence)
    for variant in tree_sequence.variants(as_bytes=True):
        print(variant.position, end="\t")
        print("{}".format(variant.genotypes.decode()))


def run_dump_nodes(args):
    tree_sequence = tskit.load(args.tree_sequence)
    tree_sequence.dump_text(nodes=sys.stdout, precision=args.precision)


def run_dump_edges(args):
    tree_sequence = tskit.load(args.tree_sequence)
    tree_sequence.dump_text(edges=sys.stdout, precision=args.precision)


def run_dump_sites(args):
    tree_sequence = tskit.load(args.tree_sequence)
    tree_sequence.dump_text(sites=sys.stdout, precision=args.precision)


def run_dump_mutations(args):
    tree_sequence = tskit.load(args.tree_sequence)
    tree_sequence.dump_text(mutations=sys.stdout, precision=args.precision)


def run_dump_provenances(args):
    tree_sequence = tskit.load(args.tree_sequence)
    if args.human:
        for provenance in tree_sequence.provenances():
            d = json.loads(provenance.record)
            print("id={}, timestamp={}, record={}".format(
                provenance.id, provenance.timestamp, json.dumps(d, indent=4)))
    else:
        tree_sequence.dump_text(provenances=sys.stdout)


def run_dump_vcf(args):
    tree_sequence = tskit.load(args.tree_sequence)
    tree_sequence.write_vcf(sys.stdout, args.ploidy)


def run_dump_macs(args):
    """
    Write a macs formatted file so we can import into pbwt.
    """
    tree_sequence = tskit.load(args.tree_sequence)
    n = tree_sequence.get_sample_size()
    m = tree_sequence.get_sequence_length()
    print("COMMAND:\tnot_macs {} {}".format(n, m))
    print("SEED:\tASEED")
    for variant in tree_sequence.variants(as_bytes=True):
        print(
            "SITE:", variant.index, variant.position / m, 0.0,
            "{}".format(variant.genotypes.decode()), sep="\t")


def run_simulate(args):
    tree_sequence = msprime.simulate(
        sample_size=int(args.sample_size),
        Ne=args.effective_population_size,
        length=args.length,
        recombination_rate=args.recombination_rate,
        mutation_rate=args.mutation_rate,
        random_seed=args.random_seed)
    tree_sequence.dump(args.tree_sequence, zlib_compression=args.compress)


def get_msp_parser():
    top_parser = argparse.ArgumentParser(
        description="Command line interface for msprime.",
        epilog=msprime_citation_text)
    top_parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(msprime.__version__))
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    parser = subparsers.add_parser(
        "simulate",
        help="Run the simulation")
    add_sample_size_argument(parser)
    add_tree_sequence_argument(parser)
    parser.add_argument(
        "--length", "-L", type=float, default=1,
        help="The length of the simulated region in base pairs.")
    parser.add_argument(
        "--recombination-rate", "-r", type=float, default=0,
        help="The recombination rate per base per generation")
    parser.add_argument(
        "--mutation-rate", "-u", type=float, default=0,
        help="The mutation rate per base per generation")
    parser.add_argument(
        "--effective-population-size", "-N", type=float, default=1,
        help="The diploid effective population size Ne")
    parser.add_argument(
        "--random-seed", "-s", type=int, default=None,
        help="The random seed. If not specified one is chosen randomly")
    parser.add_argument(
        "--compress", "-z", action="store_true",
        help="Enable zlib compression")
    parser.set_defaults(runner=run_simulate)

    parser = subparsers.add_parser(
        "vcf",
        help="Write the tree sequence out in VCF format.")
    add_tree_sequence_argument(parser)
    parser.add_argument(
        "--ploidy", "-P", type=int, default=1,
        help="The ploidy level of samples")
    parser.set_defaults(runner=run_dump_vcf)

    parser = subparsers.add_parser(
        "nodes",
        help="Dump nodes in tabular format.")
    add_tree_sequence_argument(parser)
    add_precision_argument(parser)
    parser.set_defaults(runner=run_dump_nodes)

    parser = subparsers.add_parser(
        "edges",
        help="Dump edges in tabular format.")
    add_tree_sequence_argument(parser)
    add_precision_argument(parser)
    parser.set_defaults(runner=run_dump_edges)

    parser = subparsers.add_parser(
        "sites",
        help="Dump sites in tabular format.")
    add_tree_sequence_argument(parser)
    add_precision_argument(parser)
    parser.set_defaults(runner=run_dump_sites)

    parser = subparsers.add_parser(
        "mutations",
        help="Dump mutations in tabular format.")
    add_tree_sequence_argument(parser)
    add_precision_argument(parser)
    parser.set_defaults(runner=run_dump_mutations)

    parser = subparsers.add_parser(
        "provenances",
        help="Dump provenance information in tabular format.")
    add_tree_sequence_argument(parser)
    parser.add_argument(
        "-H", "--human", action="store_true",
        help="Print out the provenances in a human readable format")
    parser.set_defaults(runner=run_dump_provenances)

    parser = subparsers.add_parser(
        "haplotypes",
        help="Dump haplotypes in text format.")
    add_tree_sequence_argument(parser)
    parser.set_defaults(runner=run_dump_haplotypes)

    parser = subparsers.add_parser(
        "variants",
        help="Dump variants in text format.")
    add_tree_sequence_argument(parser)
    parser.set_defaults(runner=run_dump_variants)

    parser = subparsers.add_parser(
        "macs",
        help="Dump results in MaCS format.")
    add_tree_sequence_argument(parser)
    parser.set_defaults(runner=run_dump_macs)

    parser = subparsers.add_parser(
        "newick",
        help="Dump results in newick format.")
    add_tree_sequence_argument(parser)
    parser.add_argument(
        "--precision", "-p", type=int, default=3,
        help="The number of decimal places in branch lengths")
    parser.set_defaults(runner=run_dump_newick)

    parser = subparsers.add_parser(
        "upgrade",
        help="Upgrade legacy tree sequence files to the latest version.")
    parser.add_argument(
        "source", help="The source msprime tree sequence file in legacy format")
    parser.add_argument(
        "destination", help="The filename of the upgraded copy.")
    parser.add_argument(
        "--remove-duplicate-positions", "-d", action="store_true", default=False,
        help="Remove any duplicated mutation positions in the source file. ")
    parser.set_defaults(runner=run_upgrade)

    return top_parser


def msp_main(arg_list=None):
    set_sigpipe_handler()
    parser = get_msp_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
