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
Test cases for the command line interfaces to msprime
"""
from __future__ import print_function
from __future__ import division

import io
import itertools
import os
import random
import sys
import tempfile
import unittest

import msprime
import msprime.cli as cli

# We're forced to do this because dendropy doesn't support Python 3.
_dendropy_available = True
try:
    import dendropy
except ImportError:
    _dendropy_available = False


def div4(matrix):
    """
    Returns the specified matrix divided by 4. We need this for
    the rate changes in migration matrices.
    """
    return [[x / 4 for x in row] for row in matrix]


def capture_output(func, *args, **kwargs):
    """
    Runs the specified function and arguments, and returns the
    tuple (stdout, stderr) as strings.
    """
    buffer_class = io.BytesIO
    if sys.version_info[0] == 3:
        buffer_class = io.StringIO
    stdout = sys.stdout
    sys.stdout = buffer_class()
    stderr = sys.stderr
    sys.stderr = buffer_class()

    try:
        func(*args, **kwargs)
        stdout_output = sys.stdout.getvalue()
        stderr_output = sys.stderr.getvalue()
    finally:
        sys.stdout.close()
        sys.stdout = stdout
        sys.stderr.close()
        sys.stderr = stderr
    return stdout_output, stderr_output


class TestRandomSeeds(unittest.TestCase):
    """
    Test the random seed generation for the ms compatability layer.
    """
    def test_seed_conversion(self):
        num_random_tests = 100
        for max_seed in [1024, 2**16 - 1, 2**32 - 1]:
            input_seeds = set()
            python_seeds = set()
            values = set()
            for j in range(num_random_tests):
                input_seed = tuple(
                    [random.randint(1, max_seed) for k in range(3)])
                python_seed = cli.get_single_seed(input_seed)
                self.assertNotIn(input_seed, input_seeds)
                self.assertNotIn(python_seed, python_seeds)
                input_seeds.add(input_seed)
                python_seeds.add(python_seed)
                self.assertGreater(python_seed, 0)
                # Make sure it's deterministic
                python_seed2 = cli.get_single_seed(input_seed)
                self.assertEqual(python_seed, python_seed2)
                # Make sure it results in a distinct first draw.
                rng = random.Random()
                rng.seed(python_seed)
                u = rng.random()
                self.assertNotIn(u, values)
                values.add(u)
            self.assertEqual(len(values), num_random_tests)
            self.assertEqual(len(python_seeds), num_random_tests)
            self.assertEqual(len(input_seeds), num_random_tests)

    def test_seed_generation(self):
        num_random_tests = 100
        seeds = set()
        for _ in range(num_random_tests):
            s = tuple(cli.generate_seeds())
            self.assertEqual(len(set(s)), 3)
            self.assertNotIn(s, seeds)
            seeds.add(s)
        self.assertEqual(len(seeds), num_random_tests)

    def test_seed_conversion_order(self):
        seeds = set()
        for p in itertools.permutations([1, 2, 3]):
            s = cli.get_single_seed(p)
            self.assertNotIn(s, seeds)
            seeds.add(s)


class TestMspmsRoundTrip(unittest.TestCase):
    """
    Tests the mspms argument parsing to ensure that we correctly round-trip
    to the ms arguments.
    """
    @unittest.skip
    # Skip these tests until we've got the ms cmd line generation code
    # reimplemented.
    def test_msdoc_examples(self):

        arg_list = ["4", "1", "-T"]
        simulator = cli.get_mspms_runner(arg_list).get_simulator()
        line = simulator.get_ms_command_line()
        parser = cli.get_mspms_parser()
        args = parser.parse_args(line[1:])
        self.assertEqual(args.sample_size, 4)
        self.assertEqual(args.num_replicates, 1)
        self.assertEqual(args.trees, True)

        arg_list = "15 1000 -t 2.0 -eN 1.0 .1 -eN 2.0 4.0".split()
        simulator = cli.get_mspms_runner(arg_list).get_simulator()
        line = simulator.get_ms_command_line(
            num_replicates=1000, scaled_mutation_rate=2.0)
        parser = cli.get_mspms_parser()
        args = parser.parse_args(line[1:])
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 2.0)
        self.assertEqual(
            args.size_change, [(0, [1.0, 0.1]), (1, [2.0, 4.0])])

        arg_list = "15 1000 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5".split()
        simulator = cli.get_mspms_runner(arg_list).get_simulator()
        line = simulator.get_ms_command_line(
            num_replicates=1000, scaled_mutation_rate=6.4)
        parser = cli.get_mspms_parser()
        args = parser.parse_args(line[1:])
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 6.4)
        self.assertEqual(args.growth_rate, 6.93)
        self.assertEqual(args.growth_rate_change, [(0, [0.2, 0.0])])
        self.assertEqual(args.size_change, [(1, [0.3, 0.5])])


class TestMspmsArgumentParser(unittest.TestCase):
    """
    Tests the parser to ensure it works correctly and is ms compatible.
    """
    def parse_args(self, args):
        parser = cli.get_mspms_parser()
        return parser.parse_args(args)

    def test_msdoc_examples(self):
        args = self.parse_args(["4", "2", "-t", "5.0"])
        self.assertEqual(args.sample_size, 4)
        self.assertEqual(args.num_replicates, 2)
        self.assertEqual(args.mutation_rate, 5.0)
        self.assertEqual(args.trees, False)

        args = self.parse_args(["4", "2", "-T"])
        self.assertEqual(args.sample_size, 4)
        self.assertEqual(args.num_replicates, 2)
        self.assertEqual(args.mutation_rate, 0.0)
        self.assertEqual(args.trees, True)

        args = self.parse_args("15 1000 -t 10.04 -r 100.0 2501".split())
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 10.04)
        self.assertEqual(args.trees, False)
        self.assertEqual(args.recombination, [100, 2501])

        args = self.parse_args(
            "15 1000 -t 2.0 -eN 1.0 .1 -eN 2.0 4.0".split())
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 2.0)
        self.assertEqual(
            args.size_change, [(0, [1.0, 0.1]), (1, [2.0, 4.0])])

        args = self.parse_args(
            "15 1000 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5".split())
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 6.4)
        self.assertEqual(args.growth_rate, 6.93)
        self.assertEqual(args.growth_rate_change, [(0, [0.2, 0.0])])
        self.assertEqual(args.size_change, [(1, [0.3, 0.5])])

    def test_positional_arguments(self):
        args = self.parse_args(["40", "20"])
        self.assertEqual(args.sample_size, 40)
        self.assertEqual(args.num_replicates, 20)

    def test_mutations(self):
        args = self.parse_args(["40", "20"])
        self.assertEqual(args.mutation_rate, 0.0)
        args = self.parse_args(["40", "20", "-t", "10"])
        self.assertEqual(args.mutation_rate, 10.0)
        args = self.parse_args(["40", "20", "--mutation-rate=10"])
        self.assertEqual(args.mutation_rate, 10.0)

    def test_trees(self):
        args = self.parse_args(["40", "20"])
        self.assertEqual(args.trees, False)
        args = self.parse_args(["40", "20", "-T"])
        self.assertEqual(args.trees, True)
        args = self.parse_args(["40", "20", "--trees"])
        self.assertEqual(args.trees, True)

    def test_size_changes(self):
        args = self.parse_args(["40", "20"])
        self.assertEqual(args.size_change, [])
        args = self.parse_args("10 1 -eN 2.0 0.5".split())
        self.assertEqual(args.size_change, [(0, [2.0, 0.5])])
        args = self.parse_args("10 1 -eN 1.0 0.5 -eN 2.0 5.0".split())
        self.assertEqual(args.size_change, [(0, [1.0, 0.5]), (1, [2.0, 5.0])])
        args = self.parse_args("10 1 -eN 1.0 0.5 -eN 2.0 5.0".split())
        self.assertEqual(args.size_change, [(0, [1.0, 0.5]), (1, [2.0, 5.0])])
        args = self.parse_args("10 1 -I 2 10 0 -en 1 2 3".split())
        self.assertEqual(
            args.population_size_change, [(0, [1, 2, 3])])

    def test_growth_rates(self):
        args = self.parse_args(["40", "20"])
        self.assertEqual(args.growth_rate, None)
        self.assertEqual(args.growth_rate_change, [])

        args = self.parse_args("15 1000 -G 5.25".split())
        self.assertEqual(args.growth_rate, 5.25)
        self.assertEqual(args.growth_rate_change, [])

        args = self.parse_args("15 1000 -eG 1.0 5.25".split())
        self.assertEqual(args.growth_rate, None)
        self.assertEqual(args.growth_rate_change, [(0, [1.0, 5.25])])
        args = self.parse_args("15 1000 -eG 1.0 5.25 -eG 2.0 10".split())
        self.assertEqual(args.growth_rate, None)
        self.assertEqual(
            args.growth_rate_change, [(0, [1.0, 5.25]), (1, [2.0, 10.0])])
        args = self.parse_args(
            "15 1000 -eG 1.0 5.25 -eG 2.0 10 -G 4".split())
        self.assertEqual(args.growth_rate, 4.0)
        self.assertEqual(
            args.growth_rate_change, [(0, [1.0, 5.25]), (1, [2.0, 10.0])])
        args = self.parse_args("10 1 -I 2 10 0 -eg 1 2 3".split())
        self.assertEqual(
            args.population_growth_rate_change, [(0, [1, 2, 3])])

    def test_migration_rates(self):
        args = self.parse_args("15 1 -I 2 15 0 -eM 2 3 ".split())
        self.assertEqual(
            args.migration_rate_change, [(0, [2, 3])])
        args = self.parse_args(
            "15 1 -I 2 15 0 -eM 2 3 -eG 3 4 -eM 4 5".split())
        self.assertEqual(
            args.migration_rate_change, [(0, [2, 3]), (2, [4, 5])])


class CustomExceptionForTesting(Exception):
    """
    This exception class is used to check that errors are correctly
    thrown.
    """


class TestMspmsCreateSimulationRunnerErrors(unittest.TestCase):
    """
    Tests for errors that can be thrown when creating the simulation runner.
    """

    def setUp(self):
        self.parser = cli.get_mspms_parser()

        def f(message):
            # print("error:", message)
            raise CustomExceptionForTesting()
        self.parser.error = f

    def assert_parser_error(self, command_line):
        self.assertRaises(
            CustomExceptionForTesting, cli.create_simulation_runner,
            self.parser, command_line.split())

    def test_trees_or_mutations(self):
        self.assert_parser_error("10 1")
        self.assert_parser_error("10 1 -G 1")

    def test_structure(self):
        self.assert_parser_error("2 1 -T -I")
        self.assert_parser_error("2 1 -T -I x")
        self.assert_parser_error("2 1 -T -I 2")
        self.assert_parser_error("2 1 -T -I 2 1")
        self.assert_parser_error("2 1 -T -I 2 1x")
        self.assert_parser_error("2 1 -T -I 2 1x 1")
        self.assert_parser_error("2 1 -T -I 2 1 100")
        # We can also optionally have a migration rate
        self.assert_parser_error("2 1 -T -I 2 1 100 sd")
        self.assert_parser_error("2 1 -T -I 2 1 1 0.1 1")
        # Check for some higher values
        self.assert_parser_error("10 1 -T -I 4 1 1")
        self.assert_parser_error("10 1 -T -I 5 1 1 1 1 6 1 1")
        # Negative migration rates not allowed
        self.assert_parser_error("2 1 -T -I 2 1 1 -1")

    def test_migration_matrix_entry(self):
        # -m without -I raises an error
        self.assert_parser_error("10 1 -T -m 1 1 1")
        # Non int values not allowed
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 1.1 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 1 1.1 1")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 0 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 1 0 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 3 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 1 3 1")
        # Diagonal elements cannot be set.
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 1 1 1")
        # Negative rates not allowed
        self.assert_parser_error("10 1 -T -I 2 10 0 -m 1 2 -1")

    def test_migration_matrix_entry_change(self):
        # -em without -I raises an error
        self.assert_parser_error("10 1 -T -em 1 1 1 1")
        # Non int values not allowed
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 1.1 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 1 1.1 1")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 0 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 1 0 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 3 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 1 3 1")
        # Diagonal elements cannot be set.
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 1 1 1")
        # Negative rates not allowed
        self.assert_parser_error("10 1 -T -I 2 10 0 -em 1 1 2 -1")

    def test_migration_matrix(self):
        # -ma without -I raises an error
        self.assert_parser_error("10 1 -T -ma 1 1 1")
        # Incorrect lengths
        self.assert_parser_error("10 1 -T -I 2 5 5 -ma ")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ma 0 0")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ma 0 0 0")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ma 0 0 0 0 0")
        # Non float values in non-diagonals not allowed
        self.assert_parser_error("10 1 -T -I 2 5 5 -ma 0 x 0 0")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ma 0 0 x 0")
        # Negative values
        self.assert_parser_error("10 1 -T -I 2 5 5 -ma 0 -1 0 0")

    def test_migration_matrix_change(self):
        # -ema without -I raises an error
        self.assert_parser_error("10 1 -T -ema 1 1 1 1")
        # Incorrect lengths
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema 1 ")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema 1 0 0")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema 1 0 0 0")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema 1 0 0 0 0 0")
        # Non float values in non-diagonals not allowed
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema 1 2 0 x 0 0")
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema 1 2 0 0 x 0")
        # Negative values
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema 1 2 0 -1 0 0")
        # Non float times.
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema x 2 0 0 0 0")
        # Change in migration matrix size.
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema x 1 0")
        self.assert_parser_error(
            "10 1 -T -I 2 5 5 -ema x 3 0 0 0 0 0 0 0 0 0")

    def test_migration_rate_change(self):
        # -eM without -I raises error
        self.assert_parser_error("10 1 -T -eM 1 1")
        # Non int values
        self.assert_parser_error("10 1 -T -I 2 10 0 -eM 1 x")
        self.assert_parser_error("10 1 -T -I 2 10 0 -eM x 1")
        # Wrong number of args
        self.assert_parser_error("10 1 -T -I 2 10 0 -eM 1 1 1")
        # Negative migration rates not allowed
        self.assert_parser_error("10 1 -T -I 2 10 0 -eM 1 -1")
        # Negative times are not allowed
        self.assert_parser_error("10 1 -T -I 2 10 0 -eM -1 1")

    def test_unordered_demographic_events(self):
        self.assert_parser_error("10 1 -T -eN 0.2 1 -eN 0.1 1")
        self.assert_parser_error("10 1 -T -eG 0.2 1 -eN 0.1 1")
        self.assert_parser_error("10 1 -T -eG 0.2 1 -eG 0.1 1")
        self.assert_parser_error("10 1 -T -eG 0.1 1 -eN 0.21 1 -eG 0.2 1")
        self.assert_parser_error("10 1 -T -eG 0.1 1 -eG 0.21 1 -eG 0.2 1")
        self.assert_parser_error(
            "10 1 -T -I 2 10 0 -eG 0.1 1 -eM 0.21 1 -eG 0.2 1")

    def test_recombination(self):
        self.assert_parser_error("10 1 -T -r x 20")
        # Cannot have non-integer numbers of loci
        self.assert_parser_error("10 1 -T -r 1 x")
        self.assert_parser_error("10 1 -T -r 1 x")
        self.assert_parser_error("10 1 -T -r 1 1.1")
        # Number of loci must be > 1
        self.assert_parser_error("10 1 -T -r 1 0")
        self.assert_parser_error("10 1 -T -r 1 1")
        self.assert_parser_error("10 1 -T -r 1 -1")

    def test_population_growth_rate(self):
        self.assert_parser_error("10 1 -T -I 2 10 0 -g 1 x")
        self.assert_parser_error("10 1 -T -I 2 10 0 -g x 1")
        # Non int values not allowed for pop_id
        self.assert_parser_error("10 1 -T -I 2 10 0 -g 1.1 1.1")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -g 0 1.1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -g -1 1.1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -g 3 1.1")
        self.assert_parser_error("10 1 -T -I 4 10 0 0 0 -g 5 1.1")

    def test_population_size(self):
        self.assert_parser_error("10 1 -T -I 2 10 0 -n 1 x")
        self.assert_parser_error("10 1 -T -I 2 10 0 -n x 1")
        # Non int values not allowed for pop_id
        self.assert_parser_error("10 1 -T -I 2 10 0 -n 1.1 1.1")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -n 0 1.1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -n -1 1.1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -n 3 1.1")
        self.assert_parser_error("10 1 -T -I 4 10 0 0 0 -n 5 1.1")

    def test_population_growth_rate_change(self):
        self.assert_parser_error("10 1 -T -I 2 10 0 -eg")
        self.assert_parser_error("10 1 -T -I 2 10 0 -eg 1 1 x")
        self.assert_parser_error("10 1 -T -I 2 10 0 -eg x 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -eg 1 x 1")
        # Non int values not allowed for pop_id
        self.assert_parser_error("10 1 -T -I 2 10 0 -eg 0.1 1.1 1.1")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -eg 0.1 -1 1.1")
        # Negative times raise an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -eg -1 1 1")

    def test_population_size_change(self):
        self.assert_parser_error("10 1 -T -I 2 10 0 -en")
        self.assert_parser_error("10 1 -T -I 2 10 0 -en 1 1 x")
        self.assert_parser_error("10 1 -T -I 2 10 0 -en x 1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -en 1 x 1")
        # Non int values not allowed for pop_id
        self.assert_parser_error("10 1 -T -I 2 10 0 -en 0.1 1.1 1.1")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -en 0.1 -1 1.1")
        # Negative times raise an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -en -1 1 1")

    def test_population_split(self):
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej")
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej 1 ")
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej 1 2")
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej 1 2 3 4")
        # Non int values not allowed for pop_id
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej 0.1 1.1 2")
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej 0.1 1 2.2")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej 0.1 -1 2")
        # Negative times raise an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -ej -1 1 2")

    def test_admixture(self):
        self.assert_parser_error("10 1 -T -es")
        self.assert_parser_error("10 1 -T -es 1")
        self.assert_parser_error("10 1 -T -es 1 1")
        self.assert_parser_error("10 1 -T -es 1 1 1 2")
        # Non int values not allowed for pop_id
        self.assert_parser_error("10 1 -T -I 2 10 0 -es 0.1 1.1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -es 0.1 x 1")
        # Out-of-bounds raises an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -es 0.1 -1 1")
        self.assert_parser_error("10 1 -T -I 2 10 0 -es 0.1 3 1")
        # After an -ej, num_populations is increased by one
        self.assert_parser_error("10 1 -T -I 2 10 0 -es 0.1 2 1 -en 0.2 4 1")
        # Negative times raise an error
        self.assert_parser_error("10 1 -T -I 2 10 0 -es -1 1 1")
        # We don't support -es and any options that affect all pops.
        self.assert_parser_error("10 1 -T -I 2 10 0 -es 1 1 1 -eM 1 2")
        self.assert_parser_error(
            "10 1 -T -I 2 10 0 -es 1 1 1 -ema 0.5 2 1 2 3 4")
        self.assert_parser_error("10 1 -t 2.0 -eG 0.001 5.0 -es 0.01 1 0.0")
        self.assert_parser_error("10 1 -t 2.0 -eN 0.001 5.0 -es 0.01 1 0.0")


class TestMspmsCreateSimulationRunner(unittest.TestCase):
    """
    Test that we correctly create a simulator instance based on the
    command line arguments.
    """

    def create_runner(self, command_line):
        parser = cli.get_mspms_parser()
        return cli.create_simulation_runner(parser, command_line.split())

    def create_simulator(self, command_line):
        return self.create_runner(command_line).get_simulator()

    def test_mutation_rates(self):
        # Mutation rates over a sequence length 1
        runner = self.create_runner("2 1 -t 1")
        self.assertEqual(runner.get_mutation_rate(), 1 / 4)
        runner = self.create_runner("2 1 -t 2")
        self.assertEqual(runner.get_mutation_rate(), 2 / 4)

        # Mutation rates over a sequence length > 1
        runner = self.create_runner("2 1 -t 2 -r 0 10")
        self.assertEqual(runner.get_mutation_rate(), (2 / 4) / 10)
        runner = self.create_runner("2 1 -t 0.2 -r 1 2")
        self.assertEqual(runner.get_mutation_rate(), (0.2 / 4) / 2)

    def test_structure_args(self):
        sim = self.create_simulator("2 1 -T")
        self.assertEqual(sim.get_sample_configuration(), [2])
        self.assertEqual(sim.get_migration_matrix(), [[0]])

        # Specifying 1 population is the same as the default.
        sim = self.create_simulator("2 1 -T -I 1 2")
        self.assertEqual(sim.get_sample_configuration(), [2])
        self.assertEqual(sim.get_migration_matrix(), [[0]])

        sim = self.create_simulator("2 1 -T -I 2 1 1")
        self.assertEqual(sim.get_sample_configuration(), [1, 1])
        self.assertEqual(sim.get_migration_matrix(), [[0, 0], [0, 0]])

        # Default migration matrix is zeros
        sim = self.create_simulator("2 1 -T -I 2 2 0")
        self.assertEqual(sim.get_migration_matrix(), [[0, 0], [0, 0]])
        self.assertEqual(sim.get_sample_configuration(), [2, 0])

        sim = self.create_simulator("2 1 -T -I 2 1 1 0.1")
        self.assertEqual(
            sim.get_migration_matrix(), div4([[0, 0.1], [0.1, 0]]))
        self.assertEqual(sim.get_sample_configuration(), [1, 1])

        # Initial migration matrix is M / (num_pops - 1)
        sim = self.create_simulator("3 1 -T -I 3 1 1 1 2")
        self.assertEqual(sim.get_sample_configuration(), [1, 1, 1])
        self.assertEqual(
            sim.get_migration_matrix(),
            div4([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
        sim = self.create_simulator("15 1 -T -I 6 5 4 3 2 1 0")
        self.assertEqual(sim.get_sample_configuration(), [5, 4, 3, 2, 1, 0])

    def test_migration_matrix_entry(self):
        sim = self.create_simulator("3 1 -T -I 2 3 0 -m 1 2 1.1 -m 2 1 9.0")
        self.assertEqual(
            sim.get_migration_matrix(), div4([[0, 1.1], [9.0, 0]]))
        sim = self.create_simulator("3 1 -T -I 3 3 0 0 -m 1 2 1.1 -m 2 1 9.0")
        self.assertEqual(
            sim.get_migration_matrix(),
            div4([[0, 1.1, 0], [9.0, 0, 0], [0, 0, 0]]))

    def test_migration_matrix(self):
        # All migration rates are divided by 4 to get per-generation rates.
        # Diagonal values are ignored
        sim = self.create_simulator("2 1 -T -I 2 2 0 -ma 0 1 2 3")
        self.assertEqual(sim.get_migration_matrix(), div4([[0, 1], [2, 0]]))
        sim = self.create_simulator("2 1 -T -I 2 2 0 -ma x 1 2 x")
        self.assertEqual(sim.get_migration_matrix(), div4([[0, 1], [2, 0]]))
        sim = self.create_simulator("3 1 -T -I 3 1 1 1 -ma 1 2 3 4 5 6 7 8 9")
        self.assertEqual(
            sim.get_migration_matrix(),
            div4([[0, 2, 3], [4, 0, 6], [7, 8, 0]]))

    def test_simultaneous_events(self):
        sim = self.create_simulator("2 1 -T -eN 1 2.0 -eG 1.0 3 -eN 1 4")
        events = sim.get_demographic_events()
        self.assertEqual(len(events), 3)
        for event in events:
            self.assertEqual(event.time, 1.0 * 4)
        self.assertEqual(events[0].type, "population_parameters_change")
        self.assertEqual(events[0].initial_size, 2.0)
        self.assertEqual(events[0].growth_rate, 0)
        self.assertEqual(events[1].type, "population_parameters_change")
        self.assertEqual(events[1].growth_rate, 3 / 4)
        self.assertEqual(events[1].initial_size, None)
        self.assertEqual(events[2].type, "population_parameters_change")
        self.assertEqual(events[2].initial_size, 4)
        self.assertEqual(events[2].growth_rate, 0)

    def test_population_growth_rate(self):
        def f(args):
            sim = self.create_simulator(args)
            return [
                (c.initial_size, c.growth_rate * 4)
                for c in sim.get_population_configurations()]
        self.assertEqual(
            f("2 1 -T -I 3 2 0 0 -g 1 -1"),
            [(None, -1), (None, 0), (None, 0)])
        self.assertEqual(
            f("2 1 -T -I 4 2 0 0 0 -g 1 1 -g 2 2 -g 3 3"),
            [(None, 1), (None, 2), (None, 3), (None, 0)])
        # A -g should override a -G
        self.assertEqual(
            f("2 1 -T -I 3 2 0 0 -g 1 2 -G -1"),
            [(None, 2), (None, -1), (None, -1)])
        # The last -g should be effective
        self.assertEqual(
            f("2 1 -T -I 3 2 0 0 -g 1 1 -g 1 -1"),
            [(None, -1), (None, 0), (None, 0)])

    def test_population_size(self):
        def f(args):
            sim = self.create_simulator(args)
            return [
                (c.initial_size, c.growth_rate * 4)
                for c in sim.get_population_configurations()]
        self.assertEqual(
            f("2 1 -T -I 3 2 0 0 -n 1 2"),
            [(2, 0), (None, 0), (None, 0)])
        self.assertEqual(
            f("2 1 -T -I 4 2 0 0 0 -n 1 1 -n 2 2 -n 3 3"),
            [(1, 0), (2, 0), (3, 0), (None, 0)])
        # The last -n should be effective
        self.assertEqual(
            f("2 1 -T -I 3 2 0 0 -n 1 1 -n 1 0.1"),
            [(0.1, 0), (None, 0), (None, 0)])
        self.assertEqual(
            f("2 1 -T -I 3 2 0 0 -g 1 2 -n 1 0.1"),
            [(0.1, 2), (None, 0), (None, 0)])

    def test_population_growth_rate_change(self):
        def f(args):
            sim = self.create_simulator(args)
            return sim.get_demographic_events()
        events = f("2 1 -T -eg 0.1 1 2")
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].type, "population_parameters_change")
        self.assertEqual(events[0].growth_rate, 2.0 / 4)
        self.assertEqual(events[0].time, 0.1 * 4)
        self.assertEqual(events[0].population_id, 0)
        events = f("2 1 -T -I 2 1 1 -eg 0.1 1 2 -eg 0.2 2 3")
        self.assertEqual(len(events), 2)
        self.assertEqual(events[0].type, "population_parameters_change")
        self.assertEqual(events[0].growth_rate, 2.0 / 4)
        self.assertEqual(events[0].time, 0.1 * 4)
        self.assertEqual(events[0].population_id, 0)
        self.assertEqual(events[1].type, "population_parameters_change")
        self.assertEqual(events[1].growth_rate, 3.0 / 4)
        self.assertEqual(events[1].time, 0.2 * 4)
        self.assertEqual(events[1].population_id, 1)

    def test_population_size_change(self):
        def f(args):
            sim = self.create_simulator(args)
            return sim.get_demographic_events()
        events = f("2 1 -T -en 0.1 1 2")
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].type, "population_parameters_change")
        self.assertEqual(events[0].initial_size, 2.0)
        self.assertEqual(events[0].growth_rate, 0)
        self.assertEqual(events[0].time, 0.1 * 4)
        self.assertEqual(events[0].population_id, 0)
        events = f("2 1 -T -I 2 1 1 -en 0.1 1 2 -en 0.2 2 3")
        self.assertEqual(len(events), 2)
        self.assertEqual(events[0].type, "population_parameters_change")
        self.assertEqual(events[0].initial_size, 2.0)
        self.assertEqual(events[0].growth_rate, 0)
        self.assertEqual(events[0].time, 0.1 * 4)
        self.assertEqual(events[0].population_id, 0)
        self.assertEqual(events[1].type, "population_parameters_change")
        self.assertEqual(events[1].initial_size, 3.0)
        self.assertEqual(events[1].growth_rate, 0)
        self.assertEqual(events[1].time, 0.2 * 4)
        self.assertEqual(events[1].population_id, 1)

    def test_migration_rate_change(self):
        def check(args, results):
            sim = self.create_simulator(args)
            events = sim.get_demographic_events()
            self.assertEqual(len(events), len(results))
            for event, result in zip(events, results):
                self.assertEqual(event.type, "migration_rate_change")
                self.assertEqual(event.time, result[0] * 4)
                # Need to divide by 4 to correct rates.
                self.assertEqual(event.rate, result[1] / 4)
                self.assertEqual(event.matrix_index, result[2])
        check("2 1 -T -I 3 2 0 0 -eM 2.2 2", [(2.2, 1, None)])
        check(
            "2 1 -T -I 3 2 0 0 -eM 2.2 2 -eM 3.3 4",
            [(2.2, 1, None), (3.3, 2, None)])

    def test_migration_matrix_entry_change(self):
        def check(args, results):
            sim = self.create_simulator(args)
            events = sim.get_demographic_events()
            self.assertEqual(len(events), len(results))
            for event, result in zip(events, results):
                self.assertEqual(event.type, "migration_rate_change")
                self.assertEqual(event.time, result[0] * 4)
                # Need to divide by 4 to correct rates.
                self.assertEqual(event.rate, result[1] / 4)
                self.assertEqual(event.matrix_index, result[2])
        check("2 1 -T -I 3 2 0 0 -em 2.2 1 2 2", [(2.2, 2, (0, 1))])
        check(
            "2 1 -T -I 3 2 0 0 -eM 2.2 2 -em 3.3 3 1 5.5",
            [(2.2, 1, None), (3.3, 5.5, (2, 0))])

    def test_migration_matrix_change(self):
        def check(args, results):
            sim = self.create_simulator(args)
            # Make sure we haven't changed the initial matrix.
            matrix = sim.get_migration_matrix()
            for row in matrix:
                for entry in row:
                    self.assertEqual(entry, 0.0)
            events = sim.get_demographic_events()
            self.assertEqual(len(events), len(results))
            for event, result in zip(events, results):
                self.assertEqual(event.type, "migration_rate_change")
                self.assertEqual(event.time, result[0] * 4)
                # Need to divide by 4 to correct rates.
                self.assertEqual(event.rate, result[1] / 4)
                self.assertEqual(event.matrix_index, result[2])
        check(
            "2 1 -T -I 2 2 0 -ema 2.2 2 x 1 2 x",
            [(2.2, 1, (0, 1)), (2.2, 2, (1, 0))])
        check(
            "2 1 -T -I 3 2 0 0 -ema 2.2 3 x 1 2 3 x 4 5 6 x",
            [
                (2.2, 1, (0, 1)), (2.2, 2, (0, 2)),
                (2.2, 3, (1, 0)), (2.2, 4, (1, 2)),
                (2.2, 5, (2, 0)), (2.2, 6, (2, 1)),
            ])

    def test_population_split(self):
        def check(N, args, results):
            sim = self.create_simulator(args)
            events = sim.get_demographic_events()
            self.assertEqual(len(events), len(results) * N)
            k = 0
            for result in results:
                event = events[k]
                source = event.source
                self.assertEqual(event.type, "mass_migration")
                self.assertEqual(event.time, result[0] * 4)
                self.assertEqual(event.source, result[1])
                self.assertEqual(event.destination, result[2])
                # We also have to set the migration rates to 0 for the
                # population that didn't exist before now.
                k += 1
                for j in range(N):
                    if j != source:
                        event = events[k]
                        self.assertEqual(event.type, "migration_rate_change")
                        self.assertEqual(event.time, result[0] * 4)
                        self.assertEqual(event.rate, 0.0)
                        self.assertEqual(event.matrix_index, (j, source))
                        k += 1
        check(3, "2 1 -T -I 3 2 0 0 -ej 2.2 1 2", [(2.2, 0, 1)])
        check(
            3, "2 1 -T -I 3 2 0 0 -ej 2.2 1 2 -ej 2.3 1 3",
            [(2.2, 0, 1), (2.3, 0, 2)])
        check(
            4, "2 1 -T -I 4 2 0 0 0 -ej 2.2 1 2 -ej 2.3 1 3",
            [(2.2, 0, 1), (2.3, 0, 2)])

    def test_admixture(self):
        def check(N, args, results):
            sim = self.create_simulator(args)
            events = sim.get_demographic_events()
            self.assertEqual(sim.get_num_populations(), N)
            self.assertEqual(len(events), len(results))
            matrix = [[0 for _ in range(N)] for _ in range(N)]
            self.assertEqual(sim.get_migration_matrix(), matrix)
            for result, event in zip(results, events):
                self.assertEqual(event.type, "mass_migration")
                self.assertEqual(event.time, result[0] * 4)
                self.assertEqual(event.source, result[1])
                self.assertEqual(event.destination, result[2])
                self.assertEqual(event.proportion, result[3])
        check(2, "2 1 -T -es 2.2 1 1", [(2.2, 0, 1, 0)])
        check(
            3, "2 1 -T -es 2.2 1 1 -es 3.3 2 0",
            [(2.2, 0, 1, 0), (3.3, 1, 2, 1.0)])
        check(
            4, "2 1 -T -I 2 2 0 -es 2.2 1 1 -es 3.3 2 0",
            [(2.2, 0, 2, 0), (3.3, 1, 3, 1.0)])


class TestMspmsOutput(unittest.TestCase):
    """
    Tests the output of the ms compatible CLI.
    """

    def verify_newick_tree(self, tree, sample_size, precision):
        """
        Verifies that the specified string is a valid newick tree.
        """
        self.assertEqual(tree[-1], ";")
        if _dendropy_available:
            parsed_tree = dendropy.Tree.get_from_string(tree, schema="newick")
            leaf_labels = set(
                int(ts.label) for ts in parsed_tree.taxon_namespace)
            self.assertEqual(leaf_labels, set(range(1, sample_size + 1)))
            if precision > 0:
                self.assertGreater(parsed_tree.length(), 0)
        # TODO test the branch length precision output.

    def verify_output(
            self, sample_size=2, num_loci=1, recombination_rate=0,
            num_replicates=1, mutation_rate=0.0, print_trees=True,
            precision=3, random_seeds=[1, 2, 3]):
        """
        Runs the UI for the specified parameters, and parses the output
        to ensure it's consistent.
        """
        # TODO there is a problem here when we have a zero recombination
        # rate, as we can't convert between physical and genetic coords
        # in this case.
        sr = cli.SimulationRunner(
            sample_size=sample_size, num_loci=num_loci,
            scaled_recombination_rate=recombination_rate,
            num_replicates=num_replicates, scaled_mutation_rate=mutation_rate,
            print_trees=print_trees, precision=precision,
            random_seeds=random_seeds)
        with tempfile.TemporaryFile("w+") as f:
            sr.run(f)
            f.seek(0)
            # The first line contains the command line.
            line = f.readline().rstrip()
            self.assertEqual(line, " ".join(sys.argv))
            # The second line is three integers, equal to the seeds
            s = list(map(int, f.readline().split()))
            self.assertEqual(len(s), 3)
            if random_seeds is not None:
                self.assertEqual(s, random_seeds)
            # Now we've got a bunch of replicates. Each one starts with //
            num_replicates_found = 0
            line = next(f, None)
            while line is not None:
                # The first line is blank
                self.assertEqual(line, "\n")
                line = next(f, None)
                self.assertEqual(line, "//\n")
                num_replicates_found += 1
                # if we're displaying trees, the next set of lines should
                # be trees
                line = next(f, None)
                num_trees = 0
                total_length = 0
                while line is not None and line[0] in "([":
                    num_trees += 1
                    if num_loci == 1:
                        total_length += 1
                        self.assertEqual(line[0], "(")
                        tree = line.rstrip()
                    else:
                        self.assertEqual(line[0], "[")
                        j = line.find("]")
                        length = int(line[1:j])
                        self.assertGreater(length, 0)
                        total_length += length
                        tree = line[j + 1:].rstrip()
                    self.verify_newick_tree(tree, sample_size, precision)
                    line = next(f, None)
                self.assertEqual(total_length, num_loci)
                # if we have a non-zero mutation rate, we should have more
                # output.
                if mutation_rate > 0:
                    self.assertTrue(line.startswith("segsites: "))
                    s = int(line.split(":")[1])
                    self.assertGreaterEqual(s, 0)
                    line = next(f, None)
                    if s == 0:
                        self.assertEqual(line, "\n")
                        line = next(f, None)
                    else:
                        self.assertTrue(line.startswith("positions: "))
                        positions = line.split(":")[1].split()
                        self.assertEqual(len(positions), s)
                        for p in positions:
                            j = p.find(".")
                            if precision == 0:
                                self.assertEqual(j, -1)
                            else:
                                self.assertEqual(precision, len(p) - j - 1)
                        values = list(map(float, positions))
                        self.assertEqual(values, sorted(values))
                        for position in values:
                            self.assertGreaterEqual(position, 0.0)
                            self.assertLessEqual(position, 1.0)
                        line = next(f, None)
                        sequences_found = 0
                        while line is not None and line[0] in "01":
                            sequences_found += 1
                            sequence = line.rstrip()
                            self.assertEqual(len(sequence), s)
                            line = next(f, None)
                        self.assertEqual(sequences_found, sample_size)
            self.assertEqual(num_replicates, num_replicates_found)

    def test_zero_recombination_rate(self):
        self.verify_output(
            sample_size=10, mutation_rate=1, num_loci=10,
            recombination_rate=0, num_replicates=2)

    def test_num_replicates(self):
        for j in range(1, 10):
            self.verify_output(
                sample_size=10, mutation_rate=0, num_replicates=j)
            self.verify_output(
                sample_size=10, mutation_rate=10, num_replicates=j)
            self.verify_output(
                sample_size=10, mutation_rate=0, num_loci=10,
                recombination_rate=100, num_replicates=j)
            self.verify_output(
                sample_size=10, mutation_rate=0, num_loci=1,
                recombination_rate=1, num_replicates=j)
            self.verify_output(
                sample_size=10, mutation_rate=10, num_loci=1,
                recombination_rate=1, num_replicates=j)
            self.verify_output(
                sample_size=10, mutation_rate=10, num_loci=10,
                recombination_rate=10, num_replicates=j)

    def test_mutation_output(self):
        for n in [2, 3, 10]:
            self.verify_output(sample_size=n, mutation_rate=0.0)
            self.verify_output(sample_size=n, mutation_rate=1e-6)
            self.verify_output(sample_size=n, mutation_rate=10)

    def test_precision(self):
        for p in range(10):
            self.verify_output(mutation_rate=10, precision=p)

    def test_tree_output(self):
        for n in [2, 3, 10]:
            self.verify_output(sample_size=n, print_trees=True)
            self.verify_output(
                sample_size=n, num_loci=10, recombination_rate=10,
                print_trees=True)
            self.verify_output(
                sample_size=n, num_loci=100, recombination_rate=10,
                print_trees=True)

    def test_seeds_output(self):
        self.verify_output(random_seeds=None)
        self.verify_output(random_seeds=[2, 3, 4])

    def test_correct_streams(self):
        args = "15 1 -r 0 1.0 -eG 1.0 5.25 -eG 2.0 10 -G 4 -eN 3.0 1.0 -T"
        stdout, stderr = capture_output(cli.mspms_main, args.split())
        self.assertEqual(len(stderr), 0)
        # We've already tested the output pretty thoroughly above so a
        # simple test is fine here.
        self.assertEqual(len(stdout.splitlines()), 5)

    def test_seed_equivalence(self):
        sample_size = 10
        mutation_rate = 10
        # Run without seeds to get automatically generated seeds
        sr = cli.SimulationRunner(
            sample_size=sample_size, scaled_mutation_rate=mutation_rate)
        with tempfile.TemporaryFile("w+") as f:
            sr.run(f)
            f.seek(0)
            output1 = f.read()
        # Get the seeds
        seeds = list(map(int, output1.splitlines()[1].split()))
        # Run with the same seeds to get the same output.
        sr = cli.SimulationRunner(
            sample_size=sample_size, scaled_mutation_rate=mutation_rate,
            random_seeds=seeds)
        with tempfile.TemporaryFile("w+") as f:
            sr.run(f)
            f.seek(0)
            output2 = f.read()
        self.assertEqual(output1, output2)


class TestMspArgumentParser(unittest.TestCase):
    """
    Tests for the argument parsers in msp.
    """

    def test_simulate_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args([cmd, "10", "out.hdf5"])
        self.assertEqual(args.sample_size, 10)
        self.assertEqual(args.history_file, "out.hdf5")
        self.assertEqual(args.recombination_rate, 0.0)
        self.assertEqual(args.mutation_rate, 0.0)
        self.assertEqual(args.length, 1)
        self.assertEqual(args.effective_population_size, 1)
        self.assertEqual(args.random_seed, None)
        self.assertEqual(args.compress, False)

    def test_simulate_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args([
            cmd, "100", "out2.hdf5", "-L", "1e3", "-r", "5", "-u", "2",
            "-s", "1234", "-z", "-N", "11"])
        self.assertEqual(args.sample_size, 100)
        self.assertEqual(args.history_file, "out2.hdf5")
        self.assertEqual(args.recombination_rate, 5)
        self.assertEqual(args.length, 1000)
        self.assertEqual(args.random_seed, 1234)
        self.assertEqual(args.compress, True)
        self.assertEqual(args.effective_population_size, 11)

    def test_simulate_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args([
            cmd, "1000", "out3.hdf5",
            "--length", "1e4",
            "--recombination-rate", "6",
            "--effective-population-size", "1e5",
            "--mutation-rate", "1",
            "--random-seed", "123",
            "--compress"])
        self.assertEqual(args.sample_size, 1000)
        self.assertEqual(args.history_file, "out3.hdf5")
        self.assertEqual(args.recombination_rate, 6)
        self.assertEqual(args.length, 10000)
        self.assertEqual(args.effective_population_size, 10**5)
        self.assertEqual(args.random_seed, 123)
        self.assertEqual(args.compress, True)

    def test_records_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "records"
        history_file = "test.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, False)
        self.assertEqual(args.precision, 6)

    def test_records_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "records"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "-H", "-p", "8"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)
        self.assertEqual(args.precision, 8)

    def test_records_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "records"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "--header", "--precision", "5"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)
        self.assertEqual(args.precision, 5)

    def test_vcf_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "vcf"
        history_file = "test.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.ploidy, 1)

    def test_vcf_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "vcf"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "-P", "2"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.ploidy, 2)

    def test_vcf_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "vcf"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "--ploidy", "5"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.ploidy, 5)

    def test_mutations_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "mutations"
        history_file = "test.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, False)
        self.assertEqual(args.precision, 6)

    def test_mutations_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "mutations"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "-H", "-p", "4"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)
        self.assertEqual(args.precision, 4)

    def test_mutations_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "mutations"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "--header", "--precision", "9"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)
        self.assertEqual(args.precision, 9)

    def test_haplotypes_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "haplotypes"
        history_file = "test1.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)

    def test_variants_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "variants"
        history_file = "test1.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)

    def test_macs_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "macs"
        history_file = "test2.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)

    def test_newick_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "newick"
        history_file = "test3.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.precision, 3)

    def test_newick_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "newick"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "-p", "10"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.precision, 10)

    def test_newick_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "newick"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "--precision=5"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.precision, 5)

    def test_upgrade_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "upgrade"
        source = "in.hdf5"
        destination = "out.hdf5"
        args = parser.parse_args([
            cmd, source, destination])
        self.assertEqual(args.source, source)
        self.assertEqual(args.destination, destination)


class TestMspSimulateOutput(unittest.TestCase):
    """
    Tests the output of msp to ensure it's correct.
    """
    def setUp(self):
        fd, self._history_file = tempfile.mkstemp(
            prefix="msp_cli", suffix=".hdf5")
        os.close(fd)

    def tearDown(self):
        os.unlink(self._history_file)

    def test_run_defaults(self):
        cmd = "simulate"
        sample_size = 10
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, str(sample_size), self._history_file])
        self.assertEqual(len(stderr), 0)
        self.assertEqual(len(stdout), 0)

        tree_sequence = msprime.load(self._history_file)
        self.assertEqual(tree_sequence.get_sample_size(), sample_size)
        self.assertEqual(tree_sequence.get_sequence_length(), 1)
        self.assertEqual(tree_sequence.get_num_mutations(), 0)

    def test_simulate_short_args(self):
        cmd = "simulate"
        stdout, stdearr = capture_output(cli.msp_main, [
            cmd, "100", self._history_file, "-L", "1e2", "-r", "5", "-u", "2"])
        tree_sequence = msprime.load(self._history_file)
        self.assertEqual(tree_sequence.get_sample_size(), 100)
        self.assertEqual(tree_sequence.get_sequence_length(), 100)
        self.assertGreater(tree_sequence.get_num_mutations(), 0)


class TestMspConversionOutput(unittest.TestCase):
    """
    Tests the output of msp to ensure it's correct.
    """
    @classmethod
    def setUpClass(cls):
        cls._tree_sequence = msprime.simulate(
            10, length=10, recombination_rate=10,
            mutation_rate=10, random_seed=1)
        fd, cls._history_file = tempfile.mkstemp(
            prefix="msp_cli", suffix=".hdf5")
        os.close(fd)
        cls._tree_sequence.dump(cls._history_file)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls._history_file)

    def verify_records(self, output_records, header, precision):
        with tempfile.TemporaryFile("w+") as f:
            self._tree_sequence.write_records(f, header, precision)
            f.seek(0)
            output = f.read().splitlines()
        self.assertEqual(output, output_records)

    def test_records(self):
        cmd = "records"
        precision = 6
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        output_records = stdout.splitlines()
        self.verify_records(output_records, False, precision)
        # check the header.
        precision = 8
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file, "-H", "-p", str(precision)])
        self.assertEqual(len(stderr), 0)
        output_records = stdout.splitlines()
        self.assertEqual(
            list(output_records[0].split()),
            ["left", "right", "node", "children", "time", "population"])
        self.verify_records(output_records, True, precision)

    def verify_vcf(self, output_vcf):
        with tempfile.TemporaryFile("w+") as f:
            self._tree_sequence.write_vcf(f)
            f.seek(0)
            vcf = f.read()
        self.assertEqual(output_vcf, vcf)

    def test_vcf(self):
        cmd = "vcf"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        self.verify_vcf(stdout)

    def verify_mutations(self, output_mutations, header, precision):
        with tempfile.TemporaryFile("w+") as f:
            self._tree_sequence.write_mutations(f, header, precision)
            f.seek(0)
            output = f.read().splitlines()
        self.assertEqual(output, output_mutations)

    def test_mutations(self):
        cmd = "mutations"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file, "-p", "8"])
        self.assertEqual(len(stderr), 0)
        output_mutations = stdout.splitlines()
        self.verify_mutations(output_mutations, False, 8)
        # check the header.
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file, "-H", "-p", "8"])
        self.assertEqual(len(stderr), 0)
        output_mutations = stdout.splitlines()
        self.assertEqual(
            list(output_mutations[0].split()),
            ["position", "node"])
        self.verify_mutations(output_mutations, True, 8)

    def verify_haplotypes(self, output_haplotypes):
        haplotypes = list(self._tree_sequence.haplotypes())
        self.assertEqual(len(haplotypes), len(output_haplotypes))
        for h, line in zip(haplotypes, output_haplotypes):
            self.assertEqual(h, line)

    def test_haplotypes(self):
        cmd = "haplotypes"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        output_haplotypes = stdout.splitlines()
        self.verify_haplotypes(output_haplotypes)

    def verify_variants(self, output_variants):
        variants = [
            (v.position, v.genotypes.decode())
            for v in self._tree_sequence.variants(as_bytes=True)]
        self.assertEqual(len(variants), len(output_variants))
        for (pos, v), line in zip(variants, output_variants):
            self.assertEqual("{}\t{}".format(pos, v), line)

    def test_variants(self):
        cmd = "variants"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        output_variants = stdout.splitlines()
        self.verify_variants(output_variants)

    def verify_newick(self, output_newick):
        newick = list(self._tree_sequence.newick_trees())
        self.assertEqual(len(newick), len(output_newick))
        for (l, tree), line in zip(newick, output_newick):
            self.assertEqual(tree, line)

    def test_newick(self):
        cmd = "newick"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        output_newick = stdout.splitlines()
        self.verify_newick(output_newick)

    def test_macs(self):
        cmd = "macs"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        output = stdout.splitlines()
        self.assertTrue(output[0].startswith("COMMAND:"))
        self.assertTrue(output[1].startswith("SEED:"))
        self.assertEqual(
            len(output), 2 + self._tree_sequence.get_num_mutations())
        n = self._tree_sequence.get_sample_size()
        m = self._tree_sequence.get_sequence_length()
        mutations = list(self._tree_sequence.mutations())
        haplotypes = list(self._tree_sequence.haplotypes())
        for site, line in enumerate(output[2:]):
            splits = line.split()
            self.assertEqual(splits[0], "SITE:")
            self.assertEqual(int(splits[1]), site)
            position = mutations[site][0] / m
            self.assertAlmostEqual(float(splits[2]), position)
            col = splits[4]
            self.assertEqual(len(col), n)
            for j in range(n):
                self.assertEqual(col[j], haplotypes[j][site])


class TestUpgrade(unittest.TestCase):
    """
    Tests the results of the upgrade operation to ensure they are
    correct.
    """
    def test_conversion(self):
        ts1 = msprime.simulate(10)
        with tempfile.NamedTemporaryFile("w+") as v2_file, \
                tempfile.NamedTemporaryFile("w+") as v3_file:
            msprime.dump_legacy(ts1, v2_file.name)
            stdout, stderr = capture_output(cli.msp_main, [
                "upgrade", v2_file.name, v3_file.name])
            ts2 = msprime.load(v3_file.name)
        self.assertEqual(stdout, "")
        # We get some cruft on stderr that comes from h5py. This only happens
        # because we're mixing h5py and msprime for this test, so we can ignore
        # it.
        # self.assertEqual(stderr, "")
        # Quick checks to ensure we have the right tree sequence.
        # More thorough checks are done elsewhere.
        self.assertEqual(ts1.get_sample_size(), ts2.get_sample_size())
        self.assertEqual(ts1.get_num_records(), ts2.get_num_records())
        self.assertEqual(ts1.get_num_trees(), ts2.get_num_trees())
