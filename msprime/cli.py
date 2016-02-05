#
# Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
from __future__ import division
from __future__ import print_function

import argparse
import os
import random
import signal
import struct
import sys

import msprime


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


def positive_int(value):
    int_value = int(float(value))
    if int_value <= 0:
        msg = "{0} in an invalid postive integer value".format(value)
        raise argparse.ArgumentTypeError(msg)
    return int_value


def add_sample_size_argument(parser):
    parser.add_argument(
        "sample_size", type=positive_int,
        help="The number of individuals in the sample")


def add_history_file_argument(parser):
    parser.add_argument(
        "history_file", help="The msprime history file in HDF5 format")


def add_header_argument(parser):
    parser.add_argument(
        "--header", "-H", action="store_true", default=False,
        help="Print a header line in the output.")


def add_max_memory_argument(parser):
    parser.add_argument(
        "--max-memory", "-M", default=None,
        help=(
            "Maximum memory to use. If the simulation exceeds this limit "
            "exit with error status. Supports K,M and G suffixes"))


def get_seeds(random_seeds):
    """
    Takes the specified command line seeds and truncates them to 16
    bits if necessary. Then convert them to a single value that can
    be used to seed the python random number generator, and return
    both values.
    """
    max_seed = 2**16 - 1
    if random_seeds is None:
        seeds = [random.randint(1, max_seed) for j in range(3)]
    else:
        # Follow ms behaviour and truncate back to shorts
        seeds = [s if s < max_seed else max_seed for s in random_seeds]
    # Combine the together to get a 64 bit number
    seed = struct.unpack(">Q", struct.pack(">HHHH", 0, *seeds))[0]
    return seed, seeds


class SimulationRunner(object):
    """
    Class to run msprime simulation and output the results.
    """
    def __init__(
            self, sample_size=1, num_loci=1, recombination_rate=0,
            num_replicates=1, migration_matrix=None,
            population_configurations=None, demographic_events=None,
            mutation_rate=0, print_trees=False, max_memory="16M",
            precision=3, random_seeds=None):
        self._sample_size = sample_size
        self._num_loci = num_loci
        self._num_replicates = num_replicates
        self._recombination = recombination_rate
        self._mutation_rate = mutation_rate
        self._simulator = msprime.TreeSimulator(sample_size)
        self._simulator.set_max_memory(max_memory)
        self._simulator.set_num_loci(num_loci)
        self._simulator.set_migration_matrix(migration_matrix)
        self._simulator.set_population_configurations(
            population_configurations)
        self._simulator.set_demographic_events(demographic_events)
        self._simulator.set_scaled_recombination_rate(recombination_rate)
        self._precision = precision
        self._print_trees = print_trees
        # sort out the random seeds
        python_seed, ms_seeds = get_seeds(random_seeds)
        self._ms_random_seeds = ms_seeds
        random.seed(python_seed)

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

    def run(self, output):
        """
        Runs the simulations and writes the output to the specified
        file handle.
        """
        # The first line of ms's output is the command line.
        print(" ".join(sys.argv), file=output)
        print(" ".join(str(s) for s in self._ms_random_seeds), file=output)
        for j in range(self._num_replicates):
            self._simulator.set_random_seed(random.randint(0, 2**30))
            self._simulator.run()
            tree_sequence = self._simulator.get_tree_sequence()
            breakpoints = self._simulator.get_breakpoints()
            print(file=output)
            print("//", file=output)
            if self._print_trees:
                iterator = tree_sequence.newick_trees(
                    self._precision, breakpoints)
                if self._num_loci == 1:
                    for l, ns in iterator:
                        print(ns, file=output)
                else:
                    for l, ns in iterator:
                        # Print these seperately to avoid the cost of creating
                        # another string.
                        print("[{0}]".format(l), end="", file=output)
                        print(ns, file=output)
            if self._mutation_rate > 0:
                # The mutation rate in ms is multiplied by the size of the
                # region
                seed = self._simulator.get_random_seed()
                tree_sequence.generate_mutations(self._mutation_rate, seed)
                hg = msprime.HaplotypeGenerator(tree_sequence)
                s = tree_sequence.get_num_mutations()
                print("segsites:", s, file=output)
                if s != 0:
                    print("positions: ", end="", file=output)
                    positions = [
                        p / self._num_loci for p, _ in
                        tree_sequence.get_mutations()]
                    positions.sort()
                    for position in positions:
                        print(
                            "{0:.{1}f}".format(position, self._precision),
                            end=" ", file=output)
                    print(file=output)
                    for h in hg.haplotypes():
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


def create_simulation_runner(parser, arg_list):
    """
    Parses the arguments and returns a SimulationRunner instance.
    """
    args = parser.parse_args(arg_list)
    if args.mutation_rate == 0 and not args.trees:
        parser.error("Need to specify at least one of --theta or --trees")
    num_loci = int(args.recombination[1])
    r = 0.0
    # We don't scale recombination or mutation rates by the size
    # of the region.
    if num_loci > 1:
        r = args.recombination[0] / (num_loci - 1)
    mu = args.mutation_rate / num_loci

    # Check the structure format.
    symmetric_migration_rate = 0.0
    num_populations = 1
    population_configurations = [
        msprime.PopulationConfiguration(args.sample_size)]
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
        elif len(args.structure) > num_populations + 2:
            parser.error("Too many arguments to --structure/-I")
        migration_matrix = [[
            symmetric_migration_rate / (num_populations - 1) * int(j != k)
            for j in range(num_populations)] for k in range(num_populations)]
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
        if len(args.migration_matrix) != num_populations**2:
            parser.error(
                "Must be num_populations^2 migration matrix entries")
        for j in range(num_populations):
            for k in range(num_populations):
                if j != k:
                    migration_matrix[j][k] = convert_float(
                        args.migration_matrix[j * num_populations + k],
                        parser)
    for matrix_entry in args.migration_matrix_entry:
        dest = int(matrix_entry[0])
        source = int(matrix_entry[1])
        rate = matrix_entry[2]
        msg = "Bad population ID '{}': must be an integer"
        if dest != matrix_entry[0]:
            parser.error(msg.format(matrix_entry[0]))
        if source != matrix_entry[1]:
            parser.error(msg.format(matrix_entry[1]))
        msg = "Bad population ID '{}': must be 1 to num_populations"
        if source < 1 or source > num_populations:
            parser.error(msg.format(source))
        if dest < 1 or dest > num_populations:
            parser.error(msg.format(dest))
        if dest == source:
            parser.error("Cannot set diagonal elements in migration matrix")
        migration_matrix[dest - 1][source - 1] = rate

    # Get the demography parameters
    # TODO for strict ms compatability, we need to resolve the command line
    # ordering of arguments when there are two or more events with the
    # same start time. We should resolve this.
    demographic_events = []
    if args.growth_rate is not None:
        for config in population_configurations:
            config.growth_rate = args.growth_rate
    for t, alpha in args.growth_event:
        demographic_events.append(
            msprime.GrowthRateChangeEvent(t, alpha))
    for t, x in args.size_event:
        demographic_events.append(
            msprime.SizeChangeEvent(t, x))
    demographic_events.sort(key=lambda config: config.time)

    runner = SimulationRunner(
        sample_size=args.sample_size,
        num_loci=num_loci,
        migration_matrix=migration_matrix,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        num_replicates=args.num_replicates,
        recombination_rate=r,
        mutation_rate=mu,
        precision=args.precision,
        max_memory=args.max_memory,
        print_trees=args.trees,
        random_seeds=args.random_seeds)
    return runner


def get_mspms_parser():
    parser = argparse.ArgumentParser(description=mscompat_description)
    add_sample_size_argument(parser)
    parser.add_argument(
        "num_replicates", type=positive_int,
        help="Number of independent replicates")
    parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(msprime.__version__))

    group = parser.add_argument_group("Behaviour")
    group.add_argument(
        "--mutation-rate", "-t", type=float, metavar="THETA",
        help="Mutation rate theta=4*N0*mu", default=0)
    group.add_argument(
        "--trees", "-T", action="store_true",
        help="Print out trees in Newick format")
    group.add_argument(
        "--recombination", "-r", type=float, nargs=2, default=(0, 1),
        metavar=("RHO", "NUM_LOCI"), help=mscompat_recombination_help)

    group = parser.add_argument_group("Structure and migration")
    group.add_argument(
        "--structure", "-I", nargs='+',
        help=(
            "Sample from populations with the specified deme structure. "
            "The STRUCTURE argument is of the form 'num_populations "
            "n1 n2 ... [4N0m]', specifying the number of populations, "
            "the sample configuration, and optionally, the migration "
            "rate for a symmetric island model"))
    group.add_argument(
        "--migration-matrix-entry", "-m", action="append",
        metavar=("DEST", "SOURCE", "RATE"),
        nargs=3, type=float, default=[],
        help=(
            "Sets an entry M[DEST, SOURCE] in the migration matrix to the "
            "specified ate. SOURCE and DEST are (1-indexed) population "
            "IDs. Multiple options can be specified."))
    group.add_argument(
        "--migration-matrix", "-ma", nargs='+', default=None,
        help=(
            "Sets the migration matrix to the specified value. The "
            "entries are in the order M[1,1], M[1, 2], ..., M[2, 1],"
            "M[2, 2], ..., M[N, N], where N is the number of populations."))

    group = parser.add_argument_group("Demography")
    group.add_argument(
        "--growth-rate", "-G", metavar="alpha", type=float,
        help="Population growth rate alpha.")
    group.add_argument(
        "--growth-event", "-eG", nargs=2, action="append",
        type=float, default=[], metavar=("t", "alpha"),
        help="Set the growth rate to alpha at time t")
    group.add_argument(
        "--size-event", "-eN", nargs=2, action="append",
        type=float, default=[], metavar=("t", "x"),
        help="Set the population size to x * N0 at time t")
    group = parser.add_argument_group("Miscellaneous")
    group.add_argument(
        "--random-seeds", "-seeds", nargs=3, type=positive_int,
        metavar=("x1", "x2", "x3"),
        help="Random seeds (must be three integers)")
    group.add_argument(
        "--precision", "-p", type=positive_int, default=3,
        help="Number of values after decimal place to print")
    add_max_memory_argument(group)
    return parser


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

def run_dump_newick(args):
    tree_sequence = msprime.load(args.history_file)
    for l, ns in tree_sequence.newick_trees(args.precision):
        print(ns)


def run_dump_haplotypes(args):
    tree_sequence = msprime.load(args.history_file)
    for h in tree_sequence.haplotypes():
        print(h)


def run_dump_records(args):
    tree_sequence = msprime.load(args.history_file)
    if args.header:
        print("l", "r", "u", "c1", "c2", "t", sep="\t")
    for l, r, u, c, t in tree_sequence.records():
        print(l, r, u, c[0], c[1], t, sep="\t")


def run_dump_mutations(args):
    tree_sequence = msprime.load(args.history_file)
    if args.header:
        print("x", "u", sep="\t")
    for position, node in tree_sequence.mutations():
        print(position, node, sep="\t")


def run_dump_macs(args):
    """
    Write a macs formatted file so we can import into pbwt.
    """
    tree_sequence = msprime.load(args.history_file)
    n = tree_sequence.get_sample_size()
    m = tree_sequence.get_num_loci()
    print("COMMAND:\tnot_macs {} {}".format(n, m))
    print("SEED:\tASEED")
    site = 0
    for tree in tree_sequence.trees():
        for position, node in tree.mutations():
            h = ['0' for _ in range(n)]
            for u in tree.leaves(node):
                h[u - 1] = '1'
            print(
                "SITE:", site, position / m, 0.0, "".join(h), sep="\t"
            )
            site += 1


def run_simulate(args):
    tree_sequence = msprime.simulate(
        sample_size=int(args.sample_size),
        num_loci=int(args.num_loci),
        scaled_recombination_rate=args.recombination_rate,
        scaled_mutation_rate=args.mutation_rate,
        max_memory=args.max_memory,
        random_seed=args.random_seed)
    tree_sequence.dump(args.history_file, zlib_compression=args.compress)


def get_msp_parser():
    parser = argparse.ArgumentParser(
        description="Command line interface for msprime.")
    parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(msprime.__version__))
    # This is required to get uniform behaviour in Python2 and Python3
    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    simulate_parser = subparsers.add_parser(
        "simulate",
        help="Run the simulation")
    add_sample_size_argument(simulate_parser)
    add_history_file_argument(simulate_parser)
    simulate_parser.add_argument(
        "--num-loci", "-m", type=float, default=1,
        help="The number of non-recombining loci")
    simulate_parser.add_argument(
        "--recombination-rate", "-r", type=float, default=0,
        help=(
            "The rate at which recombination occurs between adjacent loci "
            "in units of 4N generations"))
    simulate_parser.add_argument(
        "--mutation-rate", "-u", type=float, default=0,
        help=(
            "The rate at which mutations occur within a single locus "
            "in units of 4N generations"))
    simulate_parser.add_argument(
        "--random-seed", "-s", type=int, default=None,
        help="The random seed. If not specified one is chosen randomly")
    add_max_memory_argument(simulate_parser)
    simulate_parser.add_argument(
        "--compress", "-z", action="store_true",
        help="Enable HDF5's transparent zlib compression")
    simulate_parser.set_defaults(runner=run_simulate)

    records_parser = subparsers.add_parser(
        "records",
        help="Dump records in tabular format.")
    add_history_file_argument(records_parser)
    add_header_argument(records_parser)
    records_parser.set_defaults(runner=run_dump_records)

    mutations_parser = subparsers.add_parser(
        "mutations",
        help="Dump mutations in tabular format.")
    add_history_file_argument(mutations_parser)
    add_header_argument(mutations_parser)
    mutations_parser.set_defaults(runner=run_dump_mutations)

    haplotypes_parser = subparsers.add_parser(
        "haplotypes",
        help="Dump results in tabular format.")
    add_history_file_argument(haplotypes_parser)
    haplotypes_parser.set_defaults(runner=run_dump_haplotypes)

    macs_parser = subparsers.add_parser(
        "macs",
        help="Dump results in MaCS format.")
    add_history_file_argument(macs_parser)
    macs_parser.set_defaults(runner=run_dump_macs)

    newick_parser = subparsers.add_parser(
        "newick",
        help="Dump results in newick format.")
    add_history_file_argument(newick_parser)
    newick_parser.add_argument(
        "--precision", "-p", type=int, default=3,
        help="The number of decimal places in branch lengths")
    newick_parser.set_defaults(runner=run_dump_newick)
    return parser


def msp_main(arg_list=None):
    set_sigpipe_handler()
    parser = get_msp_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
