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
Command line interfaces to the msprime library.
"""
import argparse
import hashlib
import os
import random
import signal
import sys
import warnings

import tskit

import msprime
from . import ancestry


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
    "for full compatibility."
)
mscompat_recombination_help = (
    "Recombination at rate rho=4*N0*r where r is the rate of recombination "
    "between the ends of the region being simulated; num_loci is the number "
    "of sites between which recombination can occur"
)
mscompat_gene_conversion_help = (
    "Gene conversion at rate gamma where gamma depends on the defined "
    "recombination rate rho=4*N0*r. If rho > 0, gc_recomb_ratio defines the ratio "
    "g/r, where r is the probability per generation of crossing-over and g the "
    "corresponding gene conversion probability. Gene conversions are initiated at "
    "rate gamma=rho*gc_recomb_ratio = 4*N0*r*gc_recomb_ratio. If rho = 0 the gene "
    "conversion rate is given by gamma=gc_recomb_ratio=4*N0*c where c is the rate "
    "of gene conversion initiation between the ends of the simulated region of "
    "length num_loci. If the recombination rate is not specified, standard "
    "parameters are used, i.e. rho = 0 and num_loci = 1. The length of the gene "
    "conversion tracts is geometrically distributed with mean tract_length. "
    "The mean tract_length needs to be larger than or equal to 1 for discrete "
    "genomes and larger than 0 for continuous genomes."
)
mshotcompat_hotspot_help = (
    "Recombination hotspots defined according to the msHOT format. This is "
    "defined as a sequence: n (start stop scale)+ where n is the number of "
    "hotspots and each hotspot spans [start, stop) where the recombination "
    "rate is the background recombination rate times scale. Adjacent hotspots "
    "may stop and start at the same position but must otherwise be non-overlapping "
    "and specified in ascending order."
)

msprime_citation_text = """
If you use msprime in your work, please cite the following paper:
Jerome Kelleher, Alison M Etheridge and Gilean McVean (2016), "Efficient
Coalescent Simulation and Genealogical Analysis for Large Sample Sizes",
PLoS Comput Biol 12(5): e1004842. doi: 10.1371/journal.pcbi.1004842
"""


def positive_int(value):
    int_value = int(float(value))
    if int_value <= 0:
        msg = f"{value} in an invalid postive integer value"
        raise argparse.ArgumentTypeError(msg)
    return int_value


def add_sample_size_argument(parser):
    parser.add_argument(
        "sample_size", type=positive_int, help="The number of genomes in the sample"
    )


def add_tree_sequence_argument(parser):
    parser.add_argument("tree_sequence", help="The msprime tree sequence file")


def add_precision_argument(parser):
    parser.add_argument(
        "--precision",
        "-p",
        type=int,
        default=6,
        help="The number of decimal places to print in records",
    )


def add_mutation_rate_argument(parser):
    parser.add_argument(
        "--mutation-rate",
        "-u",
        type=float,
        default=0,
        help="The mutation rate per base per generation",
    )


def add_random_seed_argument(parser):
    parser.add_argument(
        "--random-seed",
        "-s",
        type=int,
        default=None,
        help="The random seed. If not specified one is chosen randomly",
    )


def generate_seeds():
    """
    Generate seeds to seed the RNG and output on the command line.
    """
    # Pull three numbers from the SystemRandom generator
    rng = random.SystemRandom()
    return [rng.randint(0, 2 ** 31) for _ in range(3)]


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
        m.update(f"{s}:".encode())
    # Now take the integer value of this modulo 2^32, as this is
    # the largest seed value we'll accept.
    return int(m.hexdigest(), 16) % (2 ** 32)


def hotspots_to_recomb_map(hotspots, background_rate, seq_length):
    """
    Translates the background recombination rate and recombination
    hotspots according to the msHOT cli spec to a msprime
    recombination map.

    hotspots is of the form: [n] + [start, stop, factor]+
    where n is the number of hotspots and each hotspot spans [start, stop)
    with a recombination rate of background_rate * factor. Intervals between
    hotspots have a recombination rate of background_rate.
    """
    assert len(hotspots) > 0
    n_hotspots = hotspots[0]
    assert len(hotspots[1:]) == 3 * n_hotspots
    positions = []
    rates = []
    if hotspots[1] != 0:
        positions.append(0)
        rates.append(background_rate)
    for i in range(1, len(hotspots) - 2, 3):
        [start, stop, factor] = hotspots[i : i + 3]
        # Beginning hotspot
        positions.append(start)
        rates.append(factor * background_rate)

        if i == len(hotspots) - 3 or stop != hotspots[i + 3]:
            # Ending hotspot, back to normal recombination rate
            positions.append(stop)
            if stop != seq_length:
                rates.append(background_rate)

    if positions[-1] != seq_length:
        positions.append(seq_length)

    return msprime.RateMap(positions, rates)


class SimulationRunner:
    """
    Class to run msprime simulation and output the results.
    """

    def __init__(
        self,
        num_samples,
        demography=None,
        num_loci=1,
        recombination_rate=0,
        num_replicates=1,
        mutation_rate=0,
        print_trees=False,
        precision=3,
        ms_random_seeds=None,
        gene_conversion_rate=0,
        gene_conversion_tract_length=1,
        hotspots=None,
    ):
        self.num_loci = num_loci
        self.num_replicates = num_replicates
        self.mutation_rate = mutation_rate
        self.precision = precision
        self.print_trees = print_trees
        if hotspots is None:
            recomb_map = msprime.RateMap.uniform(num_loci, recombination_rate)
        else:
            recomb_map = hotspots_to_recomb_map(hotspots, recombination_rate, num_loci)

        if demography is None:
            # This is just used for testing so values don't really matter.
            demography = msprime.Demography.isolated_model([1])

        self.ms_random_seeds = ms_random_seeds
        if ms_random_seeds is None:
            self.ms_random_seeds = generate_seeds()
        # The ms command line requires three integers. We combine these into
        # a single seed.
        random_seed = get_single_seed(self.ms_random_seeds)

        # We need to get direct access to the simulator here because of the
        # "invisible" recombination breakpoints, so we can't run simulations
        # the usual way via sim_ancestry.
        self.simulator = ancestry._parse_sim_ancestry(
            samples=dict(enumerate(num_samples)),
            demography=demography,
            recombination_rate=recomb_map,
            gene_conversion_rate=gene_conversion_rate,
            gene_conversion_tract_length=gene_conversion_tract_length,
            ploidy=1,
            random_seed=random_seed,
        )

    def _print_trees(self, tree_sequence, output):
        """
        Print out the trees in ms-format from the specified tree sequence.
        When 'invisible' recombinations occur ms prints out copies of the
        same tree. Therefore, we must keep track of all breakpoints from the
        simulation and write out a tree for each one.
        """
        if self.num_loci == 1:
            tree = next(tree_sequence.trees())
            newick = tree.newick(precision=self.precision)
            print(newick, file=output)
        else:
            breakpoints = list(self.simulator.breakpoints) + [self.num_loci]
            j = 0
            for tree in tree_sequence.trees():
                newick = tree.newick(precision=self.precision)
                left, right = tree.interval
                while j < len(breakpoints) and breakpoints[j] <= right:
                    length = int(breakpoints[j] - left)
                    left = breakpoints[j]
                    j += 1
                    print(f"[{length}]", end="", file=output)
                    print(newick, file=output)

    def run(self, output):
        """
        Runs the simulations and writes the output to the specified
        file handle.
        """
        # The first line of ms's output is the command line.
        print(" ".join(sys.argv), file=output)
        print(" ".join(str(s) for s in self.ms_random_seeds), file=output)
        replicates = self.simulator.run_replicates(
            self.num_replicates,
            mutation_rate=self.mutation_rate,
        )
        for ts in replicates:
            print(file=output)
            print("//", file=output)
            if self.print_trees:
                self._print_trees(ts, output)
            if self.mutation_rate > 0:
                assert ts.num_sites == ts.num_mutations
                s = ts.num_sites
                print("segsites:", s, file=output)
                if s != 0:
                    print("positions: ", end="", file=output)
                    for site in ts.sites():
                        x = site.position / self.num_loci
                        print(
                            "{0:.{1}f}".format(x, self.precision),
                            end=" ",
                            file=output,
                        )
                    print(file=output)
                    for h in ts.haplotypes():
                        print(h, file=output)
                else:
                    print(file=output)


def convert_int(value, parser):
    """
    Converts the specified value to an integer if possible. If
    conversion fails, exit by calling parser.error.
    """
    try:
        return int(value)
    except ValueError:
        parser.error(f"invalid int value '{value}'")


def convert_float(value, parser):
    """
    Converts the specified value to a float if possible. If
    conversion fails, exit by calling parser.error.
    """
    try:
        return float(value)
    except ValueError:
        parser.error(f"invalid float value '{value}'")


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
    if len(input_matrix) != num_populations ** 2:
        parser.error("Must be num_populations^2 migration matrix entries")
    migration_matrix = [
        [0 for j in range(num_populations)] for k in range(num_populations)
    ]
    for j in range(num_populations):
        for k in range(num_populations):
            if j != k:
                rate = convert_float(input_matrix[j * num_populations + k], parser)
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
        "important to you.".format(other_option)
    )


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

    # ms uses a ratio to define the GC rate, but if the recombination rate
    # is zero we define the gc rate directly.
    gc_param, gc_tract_length = args.gene_conversion
    gc_rate = 0
    if r == 0.0:
        if num_loci > 1:
            gc_rate = gc_param / (num_loci - 1)
    else:
        gc_rate = r * gc_param

    demography = msprime.Demography.isolated_model([1])
    # Check the structure format.
    symmetric_migration_rate = 0.0
    num_populations = 1
    migration_matrix = [[0.0]]
    num_samples = [args.sample_size]
    if args.structure is not None:
        num_populations = convert_int(args.structure[0], parser)
        # We must have at least num_population sample_configurations
        if len(args.structure) < num_populations + 1:
            parser.error("Must have num_populations sample sizes")
        demography = msprime.Demography.isolated_model([1] * num_populations)
        num_samples = [0] * num_populations
        for j in range(num_populations):
            num_samples[j] = convert_int(args.structure[j + 1], parser)
        if sum(num_samples) != args.sample_size:
            parser.error("Population sample sizes must sum to sample_size")
        # We optionally have the overall migration_rate here
        if len(args.structure) == num_populations + 2:
            symmetric_migration_rate = convert_float(
                args.structure[num_populations + 1], parser
            )
            check_migration_rate(parser, symmetric_migration_rate)
        elif len(args.structure) > num_populations + 2:
            parser.error("Too many arguments to --structure/-I")
        if num_populations > 1:
            migration_matrix = [
                [
                    symmetric_migration_rate / (num_populations - 1) * int(j != k)
                    for j in range(num_populations)
                ]
                for k in range(num_populations)
            ]
    else:
        if len(args.migration_matrix_entry) > 0:
            parser.error(
                "Cannot specify migration matrix entries without "
                "first providing a -I option"
            )
        if args.migration_matrix is not None:
            parser.error(
                "Cannot specify a migration matrix without "
                "first providing a -I option"
            )
    if args.migration_matrix is not None:
        migration_matrix = convert_migration_matrix(
            parser, args.migration_matrix, num_populations
        )
    for matrix_entry in args.migration_matrix_entry:
        pop_i = convert_population_id(parser, matrix_entry[0], num_populations)
        pop_j = convert_population_id(parser, matrix_entry[1], num_populations)
        rate = matrix_entry[2]
        if pop_i == pop_j:
            parser.error("Cannot set diagonal elements in migration matrix")
        check_migration_rate(parser, rate)
        migration_matrix[pop_i][pop_j] = rate

    # Set the initial demography
    if args.growth_rate is not None:
        for population in demography.populations:
            population.growth_rate = args.growth_rate
    for population_id, growth_rate in args.population_growth_rate:
        pid = convert_population_id(parser, population_id, num_populations)
        demography.populations[pid].growth_rate = growth_rate
    for population_id, size in args.population_size:
        pid = convert_population_id(parser, population_id, num_populations)
        demography.populations[pid].initial_size = size

    demographic_events = []
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
        event = (index, msprime.MassMigration(t, pid, num_populations, 1 - proportion))
        demographic_events.append(event)

        num_populations += 1
        # We add another element to each row in the migration matrix
        # along with an other row. All new entries are zero.
        for row in migration_matrix:
            row.append(0)
        migration_matrix.append([0 for j in range(num_populations)])
        demography.add_population(initial_size=1)
        num_samples.append(0)

    # Add the demographic events
    for index, (t, alpha) in args.growth_rate_change:
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-eG")
        check_event_time(parser, t)
        demographic_events.append(
            (index, msprime.PopulationParametersChange(time=t, growth_rate=alpha))
        )
    for index, (t, population_id, alpha) in args.population_growth_rate_change:
        pid = convert_population_id(parser, population_id, num_populations)
        check_event_time(parser, t)
        demographic_events.append(
            (
                index,
                msprime.PopulationParametersChange(
                    time=t, growth_rate=alpha, population_id=pid
                ),
            )
        )
    for index, (t, x) in args.size_change:
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-eN")
        check_event_time(parser, t)
        demographic_events.append(
            (
                index,
                msprime.PopulationParametersChange(
                    time=t, initial_size=x, growth_rate=0
                ),
            )
        )
    for index, (t, population_id, x) in args.population_size_change:
        check_event_time(parser, t)
        pid = convert_population_id(parser, population_id, num_populations)
        demographic_events.append(
            (
                index,
                msprime.PopulationParametersChange(
                    time=t, initial_size=x, growth_rate=0, population_id=pid
                ),
            )
        )
    for index, (t, pop_i, pop_j) in args.population_split:
        check_event_time(parser, t)
        pop_i = convert_population_id(parser, pop_i, num_populations)
        pop_j = convert_population_id(parser, pop_j, num_populations)
        demographic_events.append((index, msprime.MassMigration(t, pop_i, pop_j, 1.0)))
        # Migration rates from subpopulation i (M[k, i], k != i) are set to zero.
        for k in range(num_populations):
            if k != pop_i:
                event = msprime.MigrationRateChange(t, 0.0, matrix_index=(k, pop_i))
                demographic_events.append((index, event))

    # Demographic events that affect the migration matrix
    if num_populations == 1:
        condition = (
            len(args.migration_rate_change) > 0
            or len(args.migration_matrix_entry_change) > 0
            or len(args.migration_matrix_change) > 0
        )
        if condition:
            parser.error("Cannot change migration rates for 1 population")
    for index, (t, x) in args.migration_rate_change:
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-eM")
        check_migration_rate(parser, x)
        check_event_time(parser, t)
        event = msprime.MigrationRateChange(t, x / (num_populations - 1))
        demographic_events.append((index, event))
    for index, event in args.migration_matrix_entry_change:
        t = event[0]
        check_event_time(parser, t)
        pop_i = convert_population_id(parser, event[1], num_populations)
        pop_j = convert_population_id(parser, event[2], num_populations)
        if pop_i == pop_j:
            parser.error("Cannot set diagonal elements in migration matrix")
        rate = event[3]
        check_migration_rate(parser, rate)
        msp_event = msprime.MigrationRateChange(t, rate, matrix_index=(pop_i, pop_j))
        demographic_events.append((index, msp_event))
    for index, event in args.migration_matrix_change:
        if len(event) < 3:
            parser.error("Need at least three arguments to -ma")
        if len(args.admixture) != 0:
            raise_admixture_incompatability_error(parser, "-ema")
        t = convert_float(event[0], parser)
        check_event_time(parser, t)
        if convert_int(event[1], parser) != num_populations:
            parser.error("num_populations must be equal for new migration matrix")
        matrix = convert_migration_matrix(parser, event[2:], num_populations)
        for j in range(num_populations):
            for k in range(num_populations):
                if j != k:
                    msp_event = msprime.MigrationRateChange(
                        t, matrix[j][k], matrix_index=(j, k)
                    )
                    demographic_events.append((index, msp_event))

    demographic_events.sort(key=lambda x: (x[0], x[1].time))
    time_sorted = sorted(demographic_events, key=lambda x: x[1].time)
    if demographic_events != time_sorted:
        parser.error("Demographic events must be supplied in non-decreasing time order")
    for _, event in demographic_events:
        demography.add_event(event)
    demography.migration_matrix[:] = migration_matrix

    # Adjust the population sizes so that the timescales agree. In principle
    # we could correct this with a ploidy value=0.5, but what we have here
    # seems less awful.
    for msp_event in demography.events:
        if isinstance(msp_event, msprime.PopulationParametersChange):
            if msp_event.initial_size is not None:
                msp_event.initial_size /= 2
    for pop in demography.populations:
        pop.initial_size /= 2

    runner = SimulationRunner(
        num_samples,
        demography,
        num_loci=num_loci,
        num_replicates=args.num_replicates,
        recombination_rate=r,
        mutation_rate=mu,
        gene_conversion_rate=gc_rate,
        gene_conversion_tract_length=gc_tract_length,
        precision=args.precision,
        print_trees=args.trees,
        ms_random_seeds=args.random_seeds,
        hotspots=args.hotspots,
    )
    return runner


class IndexedAction(argparse._AppendAction):
    """
    Argparse action class that allows us to find the overall ordering
    across several different options. We use this for the demographic
    events, as the order in which the events are applied matters for
    ms compatibility.
    """

    index = 0

    def __call__(self, parser, namespace, values, option_string=None):
        super().__call__(
            parser, namespace, (IndexedAction.index, values), option_string
        )
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
            except OSError as ioe:
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
        description=mscompat_description, epilog=msprime_citation_text
    )
    parser.convert_arg_line_to_args = convert_arg_line_to_args

    group = parser.add_argument_group("Behaviour")
    group.add_argument(
        "--mutation-rate",
        "-t",
        type=float,
        metavar="theta",
        help="Mutation rate theta=4*N0*mu",
        default=0,
    )
    group.add_argument(
        "--trees", "-T", action="store_true", help="Print out trees in Newick format"
    )
    group.add_argument(
        "--recombination",
        "-r",
        type=float,
        nargs=2,
        default=(0, 1),
        metavar=("rho", "num_loci"),
        help=mscompat_recombination_help,
    )
    group.add_argument(
        "--gene-conversion",
        "-c",
        type=float,
        nargs=2,
        default=(0, 1),
        metavar=("gc_recomb_ratio", "tract_length"),
        help=mscompat_gene_conversion_help,
    )
    group.add_argument(
        "--hotspots",
        "-v",
        type=float,
        nargs="+",
        default=None,
        help=mshotcompat_hotspot_help,
    )

    group = parser.add_argument_group("Structure and migration")
    group.add_argument(
        "--structure",
        "-I",
        nargs="+",
        metavar="value",
        help=(
            "Sample from populations with the specified deme structure. "
            "The arguments are of the form 'num_populations "
            "n1 n2 ... [4N0m]', specifying the number of populations, "
            "the sample configuration, and optionally, the migration "
            "rate for a symmetric island model"
        ),
    )
    group.add_argument(
        "--migration-matrix-entry",
        "-m",
        action="append",
        metavar=("i", "j", "rate"),
        nargs=3,
        type=float,
        default=[],
        help=(
            "Sets an entry M[i, j] in the migration matrix to the "
            "specified rate. i and j are (1-indexed) population "
            "IDs. Multiple options can be specified."
        ),
    )
    group.add_argument(
        "--migration-matrix",
        "-ma",
        nargs="+",
        default=None,
        metavar="entry",
        help=(
            "Sets the migration matrix to the specified value. The "
            "entries are in the order M[1,1], M[1, 2], ..., M[2, 1],"
            "M[2, 2], ..., M[N, N], where N is the number of populations."
        ),
    )
    group.add_argument(
        "--migration-rate-change",
        "-eM",
        nargs=2,
        action=IndexedAction,
        type=float,
        default=[],
        metavar=("t", "x"),
        help=(
            "Set the symmetric island model migration rate to "
            "x / (npop - 1) at time t"
        ),
    )
    group.add_argument(
        "--migration-matrix-entry-change",
        "-em",
        action=IndexedAction,
        metavar=("time", "i", "j", "rate"),
        nargs=4,
        type=float,
        default=[],
        help=(
            "Sets an entry M[i, j] in the migration matrix to the "
            "specified rate at the specified time. i and j are "
            "(1-indexed) population IDs."
        ),
    )
    group.add_argument(
        "--migration-matrix-change",
        "-ema",
        nargs="+",
        default=[],
        action=IndexedAction,
        metavar="entry",
        help=(
            "Sets the migration matrix to the specified value at time t."
            "The entries are in the order M[1,1], M[1, 2], ..., M[2, 1],"
            "M[2, 2], ..., M[N, N], where N is the number of populations."
        ),
    )

    group = parser.add_argument_group("Demography")
    group.add_argument(
        "--growth-rate",
        "-G",
        metavar="alpha",
        type=float,
        help="Set the growth rate to alpha for all populations.",
    )
    group.add_argument(
        "--population-growth-rate",
        "-g",
        action="append",
        default=[],
        nargs=2,
        metavar=("population_id", "alpha"),
        type=float,
        help="Set the growth rate to alpha for a specific population.",
    )
    group.add_argument(
        "--population-size",
        "-n",
        action="append",
        default=[],
        nargs=2,
        metavar=("population_id", "size"),
        type=float,
        help="Set the size of a specific population to size*N0.",
    )

    group.add_argument(
        "--growth-rate-change",
        "-eG",
        nargs=2,
        action=IndexedAction,
        type=float,
        default=[],
        metavar=("t", "alpha"),
        help="Set the growth rate for all populations to alpha at time t",
    )
    group.add_argument(
        "--population-growth-rate-change",
        "-eg",
        nargs=3,
        action=IndexedAction,
        type=float,
        default=[],
        metavar=("t", "population_id", "alpha"),
        help=("Set the growth rate for a specific population to " "alpha at time t"),
    )
    group.add_argument(
        "--size-change",
        "-eN",
        nargs=2,
        action=IndexedAction,
        type=float,
        default=[],
        metavar=("t", "x"),
        help="Set the population size for all populations to x * N0 at time t",
    )
    group.add_argument(
        "--population-size-change",
        "-en",
        nargs=3,
        action=IndexedAction,
        type=float,
        default=[],
        metavar=("t", "population_id", "x"),
        help=(
            "Set the population size for a specific population to " "x * N0 at time t"
        ),
    )
    group.add_argument(
        "--population-split",
        "-ej",
        nargs=3,
        action=IndexedAction,
        type=float,
        default=[],
        metavar=("t", "i", "j"),
        help=(
            "Move all lineages in population i to j at time t. "
            "Forwards in time, this corresponds to a population split "
            "in which lineages in j split into i. All migration "
            "rates for population i are set to zero."
        ),
    )
    group.add_argument(
        "--admixture",
        "-es",
        nargs=3,
        action=IndexedAction,
        type=float,
        default=[],
        metavar=("t", "population_id", "proportion"),
        help=(
            "Split the specified population into a new population, such "
            "that the specified proportion of lineages remains in "
            "the population population_id. Forwards in time this "
            "corresponds to an admixture event. The new population has ID "
            "num_populations + 1. Migration rates to and from the new "
            "population are set to 0, and growth rate is 0 and the "
            "population size for the new population is N0."
        ),
    )

    group = parser.add_argument_group("Miscellaneous")
    group.add_argument(
        "--random-seeds",
        "-seeds",
        nargs=3,
        type=positive_int,
        metavar=("x1", "x2", "x3"),
        help="Random seeds (must be three integers)",
    )
    group.add_argument(
        "--precision",
        "-p",
        type=positive_int,
        default=3,
        help="Number of values after decimal place to print",
    )

    # now for the parser that gets called first
    init_parser = argparse.ArgumentParser(
        description=mscompat_description,
        epilog=msprime_citation_text,
        add_help=False,
        parents=[parser],
    )
    init_parser.convert_arg_line_to_args = convert_arg_line_to_args

    add_sample_size_argument(init_parser)
    init_parser.add_argument(
        "num_replicates", type=positive_int, help="Number of independent replicates"
    )
    init_parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {msprime.__version__}"
    )
    init_parser.add_argument(
        "-f",
        "--filename",
        action=make_load_file_action(parser),
        help="Insert commands from a file at this point in the command line.",
    )

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


def run_simulate(args):

    if args.compress:
        warnings.warn(
            "The --compress option is no longer supported and does nothing. "
            "Please use the tszip utility to compress the output instead."
        )
    tree_sequence = msprime.simulate(
        sample_size=int(args.sample_size),
        Ne=args.effective_population_size,
        length=args.length,
        recombination_rate=args.recombination_rate,
        mutation_rate=args.mutation_rate,
        random_seed=args.random_seed,
    )
    tree_sequence.dump(args.tree_sequence)


def run_mutate(args):
    tree_sequence = tskit.load(args.tree_sequence)
    tree_sequence = msprime.sim_mutations(
        tree_sequence=tree_sequence,
        rate=args.mutation_rate,
        random_seed=args.random_seed,
        keep=args.keep,
        start_time=args.start_time,
        end_time=args.end_time,
        discrete_genome=args.discrete_genome,
    )
    tree_sequence.dump(args.output_tree_sequence)


def get_msp_parser():
    top_parser = argparse.ArgumentParser(
        description="Command line interface for msprime.", epilog=msprime_citation_text
    )
    top_parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {msprime.__version__}"
    )
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    add_simulate_subcommand(subparsers)
    add_mutate_subcommand(subparsers)

    return top_parser


def add_mutate_subcommand(subparsers):
    parser = subparsers.add_parser("mutate", help="Add mutations to a tree sequence.")
    add_tree_sequence_argument(parser)
    add_mutation_rate_argument(parser)
    add_random_seed_argument(parser)
    parser.add_argument(
        "output_tree_sequence",
        help="The tree sequence output file containing the new mutations",
    )
    parser.add_argument(
        "--keep",
        "-k",
        action="store_true",
        default=False,
        help="Keep mutations in input tree sequence",
    )
    parser.add_argument(
        "--discrete-genome",
        action="store_true",
        default=False,
        help="Generate mutations at only integer positions along the genome. ",
    )
    parser.add_argument(
        "--start-time",
        type=float,
        default=None,
        help="The minimum time ago at which a mutation can occur.",
    )
    parser.add_argument(
        "--end-time",
        type=float,
        default=None,
        help="The maximum time ago at which a mutation can occur.",
    )
    parser.set_defaults(runner=run_mutate)


def add_simulate_subcommand(subparsers) -> None:
    parser = subparsers.add_parser("simulate", help="Run the simulation")
    add_sample_size_argument(parser)
    add_tree_sequence_argument(parser)
    parser.add_argument(
        "--length",
        "-L",
        type=float,
        default=1,
        help="The length of the simulated region in base pairs.",
    )
    parser.add_argument(
        "--recombination-rate",
        "-r",
        type=float,
        default=0,
        help="The recombination rate per base per generation",
    )
    add_mutation_rate_argument(parser)
    parser.add_argument(
        "--effective-population-size",
        "-N",
        type=float,
        default=1,
        help="The diploid effective population size Ne",
    )
    add_random_seed_argument(parser)
    parser.add_argument(
        "--compress",
        "-z",
        action="store_true",
        help="Deprecated option with no effect; please use the tszip utility instead.",
    )
    parser.set_defaults(runner=run_simulate)


def msp_main(arg_list=None):
    set_sigpipe_handler()
    parser = get_msp_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
