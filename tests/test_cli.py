#
# Copyright (C) 2015 University of Oxford
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
import io
import itertools
import os
import random
import sys
import tempfile
import unittest
from unittest import mock

import newick
import numpy as np
import pytest
import tskit

import msprime
import msprime.cli as cli


def capture_output(func, *args, **kwargs):
    """
    Runs the specified function and arguments, and returns the
    tuple (stdout, stderr) as strings.
    """
    buffer_class = io.StringIO
    stdout = sys.stdout
    sys.stdout = buffer_class()
    stderr = sys.stderr
    sys.stderr = buffer_class()

    try:
        # Recent versions of MacOS seem to have issues with us calling signal
        # during tests.
        with mock.patch("signal.signal"):
            func(*args, **kwargs)
        stdout_output = sys.stdout.getvalue()
        stderr_output = sys.stderr.getvalue()
    finally:
        sys.stdout.close()
        sys.stdout = stdout
        sys.stderr.close()
        sys.stderr = stderr
    return stdout_output, stderr_output


class TestRandomSeeds:
    """
    Test the random seed generation for the ms compatability layer.
    """

    def test_seed_conversion(self):
        num_random_tests = 100
        for max_seed in [1024, 2 ** 16 - 1, 2 ** 32 - 1]:
            input_seeds = set()
            python_seeds = set()
            values = set()
            for _ in range(num_random_tests):
                input_seed = tuple([random.randint(1, max_seed) for k in range(3)])
                python_seed = cli.get_single_seed(input_seed)
                assert input_seed not in input_seeds
                assert python_seed not in python_seeds
                input_seeds.add(input_seed)
                python_seeds.add(python_seed)
                assert python_seed > 0
                # Make sure it's deterministic
                python_seed2 = cli.get_single_seed(input_seed)
                assert python_seed == python_seed2
                # Make sure it results in a distinct first draw.
                rng = random.Random()
                rng.seed(python_seed)
                u = rng.random()
                assert u not in values
                values.add(u)
            assert len(values) == num_random_tests
            assert len(python_seeds) == num_random_tests
            assert len(input_seeds) == num_random_tests

    def test_seed_generation(self):
        num_random_tests = 100
        seeds = set()
        for _ in range(num_random_tests):
            s = tuple(cli.generate_seeds())
            assert len(set(s)) == 3
            assert s not in seeds
            seeds.add(s)
        assert len(seeds) == num_random_tests

    def test_seed_conversion_order(self):
        seeds = set()
        for p in itertools.permutations([1, 2, 3]):
            s = cli.get_single_seed(p)
            assert s not in seeds
            seeds.add(s)


class TestCli(unittest.TestCase):
    """
    Superclass of tests for the CLI needing temp files.
    """

    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msp_cli_testcase_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)


class TestMspmsArgumentParser:
    """
    Tests the parser to ensure it works correctly and is ms compatible.
    """

    def parse_args(self, args):
        parser = cli.get_mspms_parser()
        return parser.parse_args(args)

    def test_msdoc_examples(self):
        args = self.parse_args(["4", "2", "-t", "5.0"])
        assert args.sample_size == 4
        assert args.num_replicates == 2
        assert args.mutation_rate == 5.0
        assert not args.trees

        args = self.parse_args(["4", "2", "-T"])
        assert args.sample_size == 4
        assert args.num_replicates == 2
        assert args.mutation_rate == 0.0
        assert args.trees

        args = self.parse_args("15 1000 -t 10.04 -r 100.0 2501".split())
        assert args.sample_size == 15
        assert args.num_replicates == 1000
        assert args.mutation_rate == 10.04
        assert not args.trees
        assert args.recombination == [100, 2501]

        args = self.parse_args("15 1000 -t 2.0 -eN 1.0 .1 -eN 2.0 4.0".split())
        assert args.sample_size == 15
        assert args.num_replicates == 1000
        assert args.mutation_rate == 2.0
        assert args.size_change == [(0, [1.0, 0.1]), (1, [2.0, 4.0])]

        args = self.parse_args("15 1000 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5".split())
        assert args.sample_size == 15
        assert args.num_replicates == 1000
        assert args.mutation_rate == 6.4
        assert args.growth_rate == 6.93
        assert args.growth_rate_change == [(0, [0.2, 0.0])]
        assert args.size_change == [(1, [0.3, 0.5])]

    def test_mshot_examples(self):
        args = self.parse_args("15 1000 -t 10.04".split())
        assert args.sample_size == 15
        assert args.num_replicates == 1000
        assert args.mutation_rate == 10.04
        assert not args.trees
        assert args.hotspots is None

        arg_str = "15 1000 -t 10.04 -r 100.0 25001 -v 2 100 200 10 7000 8000 20"
        args = self.parse_args(arg_str.split())
        assert args.sample_size == 15
        assert args.num_replicates == 1000
        assert args.mutation_rate == 10.04
        assert not args.trees
        assert args.recombination == [100, 25001]
        assert args.hotspots == [2.0, 100.0, 200.0, 10.0, 7000.0, 8000.0, 20.0]

    def test_positional_arguments(self):
        args = self.parse_args(["40", "20"])
        assert args.sample_size == 40
        assert args.num_replicates == 20

    def test_mutations(self):
        args = self.parse_args(["40", "20"])
        assert args.mutation_rate == 0.0
        args = self.parse_args(["40", "20", "-t", "10"])
        assert args.mutation_rate == 10.0
        args = self.parse_args(["40", "20", "--mutation-rate=10"])
        assert args.mutation_rate == 10.0

    def test_trees(self):
        args = self.parse_args(["40", "20"])
        assert not args.trees
        args = self.parse_args(["40", "20", "-T"])
        assert args.trees
        args = self.parse_args(["40", "20", "--trees"])
        assert args.trees

    def test_size_changes(self):
        args = self.parse_args(["40", "20"])
        assert args.size_change == []
        args = self.parse_args("10 1 -eN 2.0 0.5".split())
        assert args.size_change == [(0, [2.0, 0.5])]
        args = self.parse_args("10 1 -eN 1.0 0.5 -eN 2.0 5.0".split())
        assert args.size_change == [(0, [1.0, 0.5]), (1, [2.0, 5.0])]
        args = self.parse_args("10 1 -eN 1.0 0.5 -eN 2.0 5.0".split())
        assert args.size_change == [(0, [1.0, 0.5]), (1, [2.0, 5.0])]
        args = self.parse_args("10 1 -I 2 10 0 -en 1 2 3".split())
        assert args.population_size_change == [(0, [1, 2, 3])]

    def test_growth_rates(self):
        args = self.parse_args(["40", "20"])
        assert args.growth_rate is None
        assert args.growth_rate_change == []

        args = self.parse_args("15 1000 -G 5.25".split())
        assert args.growth_rate == 5.25
        assert args.growth_rate_change == []

        args = self.parse_args("15 1000 -eG 1.0 5.25".split())
        assert args.growth_rate is None
        assert args.growth_rate_change == [(0, [1.0, 5.25])]
        args = self.parse_args("15 1000 -eG 1.0 5.25 -eG 2.0 10".split())
        assert args.growth_rate is None
        assert args.growth_rate_change == [(0, [1.0, 5.25]), (1, [2.0, 10.0])]
        args = self.parse_args("15 1000 -eG 1.0 5.25 -eG 2.0 10 -G 4".split())
        assert args.growth_rate == 4.0
        assert args.growth_rate_change == [(0, [1.0, 5.25]), (1, [2.0, 10.0])]
        args = self.parse_args("10 1 -I 2 10 0 -eg 1 2 3".split())
        assert args.population_growth_rate_change == [(0, [1, 2, 3])]

    def test_migration_rates(self):
        args = self.parse_args("15 1 -I 2 15 0 -eM 2 3 ".split())
        assert args.migration_rate_change == [(0, [2, 3])]
        args = self.parse_args("15 1 -I 2 15 0 -eM 2 3 -eG 3 4 -eM 4 5".split())
        assert args.migration_rate_change == [(0, [2, 3]), (2, [4, 5])]

    def test_gene_conversion(self):
        args = self.parse_args("10 1 -r 1 100 -c 5 12".split())
        assert args.recombination == [1, 100]
        assert args.gene_conversion == [5, 12]


class CustomExceptionForTesting(Exception):
    """
    This exception class is used to check that errors are correctly
    thrown.
    """


class TestHotspotsToRecombMap(TestCli):
    def verify_map(self, recomb_map, expected_positions, expected_rates):
        assert np.array_equal(recomb_map.position, expected_positions)
        assert np.array_equal(recomb_map.rate, expected_rates)

    def test_multiple_hotspots(self):
        seq_length = 1000
        rate = 0.1
        hotspots = [2, 100, 200, 10, 700, 900, 20]
        expected_positions = [0, 100, 200, 700, 900, 1000]
        expected_rates = [0.1, 1.0, 0.1, 2.0, 0.1]
        recomb_map = cli.hotspots_to_recomb_map(hotspots, rate, seq_length)
        self.verify_map(recomb_map, expected_positions, expected_rates)

    def test_adjacent_hotspots(self):
        seq_length = 1000
        rate = 0.1
        hotspots = [2, 100, 200, 10, 200, 900, 20]
        expected_positions = [0, 100, 200, 900, 1000]
        expected_rates = [0.1, 1.0, 2.0, 0.1]
        recomb_map = cli.hotspots_to_recomb_map(hotspots, rate, seq_length)
        self.verify_map(recomb_map, expected_positions, expected_rates)

    def test_hotspot_on_left_bound(self):
        seq_length = 1000
        rate = 0.1
        hotspots = [1, 0, 200, 10]
        expected_positions = [0, 200, 1000]
        expected_rates = [1.0, 0.1]
        recomb_map = cli.hotspots_to_recomb_map(hotspots, rate, seq_length)
        self.verify_map(recomb_map, expected_positions, expected_rates)

    def test_hotspot_on_right_bound(self):
        seq_length = 1000
        rate = 0.1
        hotspots = [1, 800, 1000, 10]
        expected_positions = [0, 800, 1000]
        expected_rates = [0.1, 1.0]
        recomb_map = cli.hotspots_to_recomb_map(hotspots, rate, seq_length)
        self.verify_map(recomb_map, expected_positions, expected_rates)

    def test_hotspot_covering_whole_sequence(self):
        seq_length = 1000
        rate = 0.1
        hotspots = [1, 0, 1000, 10]
        expected_positions = [0, 1000]
        expected_rates = [1.0]
        recomb_map = cli.hotspots_to_recomb_map(hotspots, rate, seq_length)
        self.verify_map(recomb_map, expected_positions, expected_rates)


class TestMspmsCreateSimulationRunnerErrors(TestCli):
    """
    Tests for errors that can be thrown when creating the simulation runner.
    """

    def setUp(self):
        super().setUp()

        def error_handler(message):
            raise CustomExceptionForTesting()

        self.parser = cli.get_mspms_parser(error_handler)

    def assert_parser_error(self, command_line):
        split_cmd = command_line.split()
        with pytest.raises(CustomExceptionForTesting):
            cli.create_simulation_runner(
                self.parser,
                split_cmd,
            )
        with open(self.temp_file, "w") as f:
            # We're assuming the first two args are always the sample size
            # and num_replicates here.
            f.write(" ".join(split_cmd[2:]))
        with pytest.raises(CustomExceptionForTesting):
            cli.create_simulation_runner(
                self.parser,
                split_cmd[:2] + ["-f", self.temp_file],
            )

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
        self.assert_parser_error("10 1 -T -I 2 5 5 -ema x 3 0 0 0 0 0 0 0 0 0")

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
        self.assert_parser_error("10 1 -T -I 2 10 0 -eG 0.1 1 -eM 0.21 1 -eG 0.2 1")

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
        self.assert_parser_error("10 1 -T -I 2 10 0 -es 1 1 1 -ema 0.5 2 1 2 3 4")
        self.assert_parser_error("10 1 -t 2.0 -eG 0.001 5.0 -es 0.01 1 0.0")
        self.assert_parser_error("10 1 -t 2.0 -eN 0.001 5.0 -es 0.01 1 0.0")


class TestMspmsCreateSimulationRunner:
    """
    Test that we correctly create a simulator instance based on the
    command line arguments.
    """

    def create_runner(self, command_line):
        parser = cli.get_mspms_parser()
        return cli.create_simulation_runner(parser, command_line.split())

    def create_simulator(self, command_line):
        return self.create_runner(command_line).simulator

    def test_mutation_rates(self):
        # Mutation rates over a sequence length 1
        runner = self.create_runner("2 1 -t 1")
        assert runner.mutation_rate == 1
        runner = self.create_runner("2 1 -t 2")
        assert runner.mutation_rate == 2

        # Mutation rates over a sequence length > 1
        runner = self.create_runner("2 1 -t 2 -r 0 10")
        assert runner.mutation_rate == 2 / 10
        runner = self.create_runner("2 1 -t 0.2 -r 1 2")
        assert runner.mutation_rate == 0.2 / 2

    def test_recomb_map(self):
        runner = self.create_runner("15 1000 -t 10.04 -r 100.0 2501")
        uniform = msprime.RateMap([0, 2501], [0.04])
        actual = runner.simulator.recombination_map
        assert np.array_equal(actual.position, uniform.position)
        assert np.array_equal(actual.rate, uniform.rate)

        args = "15 1000 -t 10.04 -r 100.0 25001 -v 2 100 200 10 7000 8000 20"
        runner = self.create_runner(args)
        positions = [0, 100, 200, 7000, 8000, 25001]
        rates = [0.004, 0.04, 0.004, 0.08, 0.004]
        actual = runner.simulator.recombination_map
        assert np.array_equal(actual.position, positions)
        assert np.array_equal(actual.rate, rates)

        args = "15 1000 -t 10.04 -r 100.0 25001 -v 2 100 200 10 200 300 20"
        runner = self.create_runner(args)
        positions = [0, 100, 200, 300, 25001]
        rates = [0.004, 0.04, 0.08, 0.004]
        actual = runner.simulator.recombination_map
        assert np.array_equal(actual.position, positions)
        assert np.array_equal(actual.rate, rates)

        args = "15 1000 -t 10.04 -r 100.0 25001 -v 1 0 25001 0"
        runner = self.create_runner(args)
        positions = [0, 25001]
        rates = [0]
        actual = runner.simulator.recombination_map
        assert np.array_equal(actual.position, positions)
        assert np.array_equal(actual.rate, rates)

    def test_structure_args(self):
        sim = self.create_simulator("2 1 -T")
        assert sim.sample_configuration == [2]
        assert sim.demography.migration_matrix == [[0]]

        # Specifying 1 population is the same as the default.
        sim = self.create_simulator("2 1 -T -I 1 2")
        assert sim.sample_configuration == [2]
        assert sim.demography.migration_matrix == [[0]]

        sim = self.create_simulator("2 1 -T -I 2 1 1")
        assert sim.sample_configuration == [1, 1]
        np.testing.assert_array_equal(sim.demography.migration_matrix, [[0, 0], [0, 0]])

        # Default migration matrix is zeros
        sim = self.create_simulator("2 1 -T -I 2 2 0")
        np.testing.assert_array_equal(sim.demography.migration_matrix, [[0, 0], [0, 0]])
        assert sim.sample_configuration == [2, 0]

        sim = self.create_simulator("2 1 -T -I 2 1 1 0.1")
        np.testing.assert_array_equal(
            sim.demography.migration_matrix, [[0, 0.1], [0.1, 0]]
        )
        assert sim.sample_configuration == [1, 1]

        # Initial migration matrix is M / (num_pops - 1)
        sim = self.create_simulator("3 1 -T -I 3 1 1 1 2")
        assert sim.sample_configuration == [1, 1, 1]
        np.testing.assert_array_equal(
            sim.demography.migration_matrix, [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        )
        sim = self.create_simulator("15 1 -T -I 6 5 4 3 2 1 0")
        assert sim.sample_configuration == [5, 4, 3, 2, 1, 0]

    def test_migration_matrix_entry(self):
        sim = self.create_simulator("3 1 -T -I 2 3 0 -m 1 2 1.1 -m 2 1 9.0")
        np.testing.assert_array_equal(
            sim.demography.migration_matrix, [[0, 1.1], [9.0, 0]]
        )
        sim = self.create_simulator("3 1 -T -I 3 3 0 0 -m 1 2 1.1 -m 2 1 9.0")
        np.testing.assert_array_equal(
            sim.demography.migration_matrix, [[0, 1.1, 0], [9.0, 0, 0], [0, 0, 0]]
        )

    def test_migration_matrix(self):
        sim = self.create_simulator("2 1 -T -I 2 2 0 -ma 0 1 2 3")
        np.testing.assert_array_equal(sim.demography.migration_matrix, [[0, 1], [2, 0]])
        sim = self.create_simulator("2 1 -T -I 2 2 0 -ma x 1 2 x")
        np.testing.assert_array_equal(sim.demography.migration_matrix, [[0, 1], [2, 0]])
        sim = self.create_simulator("3 1 -T -I 3 1 1 1 -ma 1 2 3 4 5 6 7 8 9")
        np.testing.assert_array_equal(
            sim.demography.migration_matrix, [[0, 2, 3], [4, 0, 6], [7, 8, 0]]
        )

    def test_simultaneous_events(self):
        sim = self.create_simulator("2 1 -T -eN 1 2.0 -eG 1.0 3 -eN 1 4")
        events = sim.demography.events
        assert len(events) == 3
        for event in events:
            assert event.time == 1.0
        assert isinstance(events[0], msprime.PopulationParametersChange)
        assert events[0].initial_size == 1
        assert events[0].growth_rate == 0
        assert isinstance(events[1], msprime.PopulationParametersChange)
        assert events[1].growth_rate == 3
        assert events[1].initial_size is None
        assert isinstance(events[2], msprime.PopulationParametersChange)
        assert events[2].initial_size == 2
        assert events[2].growth_rate == 0

    def test_population_growth_rate(self):
        def f(args):
            sim = self.create_simulator(args)
            return [
                (c.initial_size * 2, c.growth_rate) for c in sim.demography.populations
            ]

        assert f("2 1 -T -I 3 2 0 0 -g 1 -1") == [(1, -1), (1, 0), (1, 0)]
        assert f("2 1 -T -I 4 2 0 0 0 -g 1 1 -g 2 2 -g 3 3") == [
            (1, 1),
            (1, 2),
            (1, 3),
            (1, 0),
        ]
        # A -g should override a -G
        assert f("2 1 -T -I 3 2 0 0 -g 1 2 -G -1") == [(1, 2), (1, -1), (1, -1)]
        # The last -g should be effective
        assert f("2 1 -T -I 3 2 0 0 -g 1 1 -g 1 -1") == [(1, -1), (1, 0), (1, 0)]

    def test_population_size(self):
        def f(args):
            sim = self.create_simulator(args)
            return [
                (c.initial_size * 2, c.growth_rate) for c in sim.demography.populations
            ]

        assert f("2 1 -T -I 3 2 0 0 -n 1 2") == [(2, 0), (1, 0), (1, 0)]
        assert f("2 1 -T -I 4 2 0 0 0 -n 1 1 -n 2 2 -n 3 3") == [
            (1, 0),
            (2, 0),
            (3, 0),
            (1, 0),
        ]
        # The last -n should be effective
        assert f("2 1 -T -I 3 2 0 0 -n 1 1 -n 1 0.1") == [(0.1, 0), (1, 0), (1, 0)]
        assert f("2 1 -T -I 3 2 0 0 -g 1 2 -n 1 0.1") == [(0.1, 2), (1, 0), (1, 0)]

    def test_population_growth_rate_change(self):
        def f(args):
            sim = self.create_simulator(args)
            return sim.demography.events

        events = f("2 1 -T -eg 0.1 1 2")
        assert len(events) == 1
        assert isinstance(events[0], msprime.PopulationParametersChange)
        assert events[0].growth_rate == 2.0
        assert events[0].time == 0.1
        assert events[0].population == 0
        events = f("2 1 -T -I 2 1 1 -eg 0.1 1 2 -eg 0.2 2 3")
        assert len(events) == 2
        assert isinstance(events[0], msprime.PopulationParametersChange)
        assert events[0].growth_rate == 2.0
        assert events[0].time == 0.1
        assert events[0].population == 0
        assert isinstance(events[1], msprime.PopulationParametersChange)
        assert events[1].growth_rate == 3.0
        assert events[1].time == 0.2
        assert events[1].population == 1

    def test_population_size_change(self):
        def f(args):
            sim = self.create_simulator(args)
            return sim.demography.events

        events = f("2 1 -T -en 0.1 1 2")
        assert len(events) == 1
        assert isinstance(events[0], msprime.PopulationParametersChange)
        assert events[0].initial_size == 2.0 / 2
        assert events[0].growth_rate == 0
        assert events[0].time == 0.1
        assert events[0].population == 0
        events = f("2 1 -T -I 2 1 1 -en 0.1 1 2 -en 0.2 2 3")
        assert len(events) == 2
        assert isinstance(events[0], msprime.PopulationParametersChange)
        assert events[0].initial_size == 2.0 / 2
        assert events[0].growth_rate == 0
        assert events[0].time == 0.1
        assert events[0].population == 0
        assert isinstance(events[1], msprime.PopulationParametersChange)
        assert events[1].initial_size == 3.0 / 2
        assert events[1].growth_rate == 0
        assert events[1].time == 0.2
        assert events[1].population == 1

    def test_migration_rate_change(self):
        def check(args, results):
            sim = self.create_simulator(args)
            events = sim.demography.events
            assert len(events) == len(results)
            for event, result in zip(events, results):
                assert isinstance(event, msprime.MigrationRateChange)
                assert event.time == result[0]
                assert event.rate == result[1]
                assert event.source == -1
                assert event.dest == -1

        check("2 1 -T -I 3 2 0 0 -eM 2.2 2", [(2.2, 1)])
        check("2 1 -T -I 3 2 0 0 -eM 2.2 2 -eM 3.3 4", [(2.2, 1), (3.3, 2)])

    def test_migration_matrix_entry_change(self):
        def check(args, results):
            sim = self.create_simulator(args)
            events = sim.demography.events
            assert len(events) == len(results)
            for event, result in zip(events, results):
                assert isinstance(event, msprime.MigrationRateChange)
                assert event.time == result[0]
                assert event.rate == result[1]
                assert (event.source, event.dest) == result[2]

        check("2 1 -T -I 3 2 0 0 -em 2.2 1 2 2", [(2.2, 2, (0, 1))])
        check(
            "2 1 -T -I 3 2 0 0 -eM 2.2 2 -em 3.3 3 1 5.5",
            [(2.2, 1, (-1, -1)), (3.3, 5.5, (2, 0))],
        )

    def test_migration_matrix_change(self):
        def check(args, results):
            sim = self.create_simulator(args)
            # Make sure we haven't changed the initial matrix.
            matrix = sim.demography.migration_matrix
            for row in matrix:
                for entry in row:
                    assert entry == 0.0
            events = sim.demography.events
            assert len(events) == len(results)
            for event, result in zip(events, results):
                assert isinstance(event, msprime.MigrationRateChange)
                assert event.time == result[0]
                assert event.rate == result[1]
                assert (event.source, event.dest) == result[2]

        check(
            "2 1 -T -I 2 2 0 -ema 2.2 2 x 1 2 x", [(2.2, 1, (0, 1)), (2.2, 2, (1, 0))]
        )
        check(
            "2 1 -T -I 3 2 0 0 -ema 2.2 3 x 1 2 3 x 4 5 6 x",
            [
                (2.2, 1, (0, 1)),
                (2.2, 2, (0, 2)),
                (2.2, 3, (1, 0)),
                (2.2, 4, (1, 2)),
                (2.2, 5, (2, 0)),
                (2.2, 6, (2, 1)),
            ],
        )

    def test_population_split(self):
        def check(N, args, results):
            sim = self.create_simulator(args)
            events = sim.demography.events
            assert len(events) == len(results) * N
            k = 0
            for result in results:
                event = events[k]
                new_pop = event.source
                assert isinstance(event, msprime.MassMigration)
                assert event.time == result[0]
                assert event.source == result[1]
                assert event.dest == result[2]
                # We also have to set the migration rates to 0 for the
                # population that didn't exist before now.
                k += 1
                for j in range(N):
                    if j != new_pop:
                        event = events[k]
                        assert isinstance(event, msprime.MigrationRateChange)
                        assert event.time == result[0]
                        assert event.rate == 0.0
                        assert event.source == j
                        assert event.dest == new_pop
                        k += 1

        check(3, "2 1 -T -I 3 2 0 0 -ej 2.2 1 2", [(2.2, 0, 1)])
        check(
            3, "2 1 -T -I 3 2 0 0 -ej 2.2 1 2 -ej 2.3 1 3", [(2.2, 0, 1), (2.3, 0, 2)]
        )
        check(
            4, "2 1 -T -I 4 2 0 0 0 -ej 2.2 1 2 -ej 2.3 1 3", [(2.2, 0, 1), (2.3, 0, 2)]
        )

    def test_admixture(self):
        def check(N, args, results):
            sim = self.create_simulator(args)
            events = sim.demography.events
            assert sim.num_populations == N
            assert len(events) == len(results)
            matrix = [[0 for _ in range(N)] for _ in range(N)]
            np.testing.assert_array_equal(sim.demography.migration_matrix, matrix)
            for result, event in zip(results, events):
                assert isinstance(event, msprime.MassMigration)
                assert event.time == result[0]
                assert event.source == result[1]
                assert event.dest == result[2]
                assert event.proportion == result[3]

        check(2, "2 1 -T -es 2.2 1 1", [(2.2, 0, 1, 0)])
        check(3, "2 1 -T -es 2.2 1 1 -es 3.3 2 0", [(2.2, 0, 1, 0), (3.3, 1, 2, 1.0)])
        check(
            4,
            "2 1 -T -I 2 2 0 -es 2.2 1 1 -es 3.3 2 0",
            [(2.2, 0, 2, 0), (3.3, 1, 3, 1.0)],
        )


class TestMspmsArgsFromFile(TestCli):
    """
    Test that parsing command line arguments from a file results gives
    the same results.
    """

    # We need to keep the arguments grouped together because they must be
    # complete when split between the file and command line.
    cmd_lines = [
        ["10", "2", "-T"],
        ["10", "2", "-t 10"],
        ["2", "1", "-T", "-I 3 2 0 0", "-ema 2.2 3 x 1 2 3 x 4 5 6 x"],
        ["2", "1", "-T", "-eN 1 2.0", "-eG 1.0 3", "-eN 1 4"],
        ["2", "1", "-T", "-I 3 2 0 0", "-ej 2.2 1 2", "-ej 2.3 1 3"],
        ["3", "10", "-I 2 3 0", "-m 1 2 1.1", "-m 2 1 9.0", "-t 5"],
    ]

    def verify_parsing(self, cmd_line_args, file_args):
        parser = cli.get_mspms_parser()
        cmd_line_result = vars(
            parser.parse_args(cmd_line_args.split() + file_args.split())
        )
        parser = cli.get_mspms_parser()
        with open(self.temp_file, "w") as f:
            f.write(file_args)
            f.flush()
        file_result = vars(
            parser.parse_args(cmd_line_args.split() + ["-f", self.temp_file])
        )
        assert cmd_line_result == file_result

    def test_empty_file(self):
        for cmd_line in self.cmd_lines:
            self.verify_parsing(" ".join(cmd_line), "")

    def test_all_options_in_file(self):
        for cmd_line in self.cmd_lines:
            self.verify_parsing(" ".join(cmd_line[:2]), " ".join(cmd_line[2:]))

    def test_middle_split(self):
        for cmd_line in self.cmd_lines:
            k = max(2, len(cmd_line) // 2)
            self.verify_parsing(" ".join(cmd_line[:k]), " ".join(cmd_line[k:]))


class TestMspmsArgsFromFileErrors(TestCli):
    """
    Tests for errors that can be thrown when reading arguments from a file.
    """

    def assert_parser_error(self, command_line):
        def error_handler(message):
            raise CustomExceptionForTesting()

        parser = cli.get_mspms_parser(error_handler)
        with pytest.raises(CustomExceptionForTesting):
            cli.create_simulation_runner(
                parser,
                command_line.split(),
            )

    def test_file_arg_in_file(self):
        with open(self.temp_file, "w") as f:
            f.write("-f otherfile")
        self.assert_parser_error(f"10 1 -f {self.temp_file}")

    def test_missing_file(self):
        self.assert_parser_error("10 1 -f /does/not/exist")


class TestMspmsOutput(TestCli):
    """
    Tests the output of the ms compatible CLI.
    """

    def verify_newick_tree(self, tree, sample_size, precision):
        """
        Verifies that the specified string is a valid newick tree.
        """
        assert tree[-1] == ";"
        newick_tree = newick.loads(tree)[0]
        leaf_names = newick_tree.get_leaf_names()
        assert sorted(leaf_names) == sorted([str(u + 1) for u in range(sample_size)])

    def verify_output(
        self,
        sample_size=2,
        num_loci=1,
        recombination_rate=0,
        num_replicates=1,
        mutation_rate=0.0,
        print_trees=True,
        precision=3,
        random_seeds=(1, 2, 3),
    ):
        """
        Runs the UI for the specified parameters, and parses the output
        to ensure it's consistent.
        """
        sr = cli.SimulationRunner(
            num_samples=[sample_size],
            demography=msprime.Demography.isolated_model([1]),
            num_loci=num_loci,
            recombination_rate=recombination_rate,
            num_replicates=num_replicates,
            mutation_rate=mutation_rate,
            print_trees=print_trees,
            precision=precision,
            ms_random_seeds=random_seeds,
        )
        with open(self.temp_file, "w+") as f:
            sr.run(f)
            f.seek(0)
            # The first line contains the command line.
            line = f.readline().rstrip()
            assert line == " ".join(sys.argv)
            # The second line is three integers, equal to the seeds
            s = tuple(map(int, f.readline().split()))
            assert len(s) == 3
            if random_seeds is not None:
                assert s == random_seeds
            # Now we've got a bunch of replicates. Each one starts with //
            num_replicates_found = 0
            line = next(f, None)
            while line is not None:
                # The first line is blank
                assert line == "\n"
                line = next(f, None)
                assert line == "//\n"
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
                        assert line[0] == "("
                        tree = line.rstrip()
                    else:
                        assert line[0] == "["
                        j = line.find("]")
                        length = int(line[1:j])
                        assert length > 0
                        total_length += length
                        tree = line[j + 1 :].rstrip()
                    self.verify_newick_tree(tree, sample_size, precision)
                    line = next(f, None)
                assert total_length == num_loci
                # if we have a non-zero mutation rate, we should have more
                # output.
                if mutation_rate > 0:
                    assert line.startswith("segsites: ")
                    s = int(line.split(":")[1])
                    assert s >= 0
                    line = next(f, None)
                    if s == 0:
                        assert line == "\n"
                        line = next(f, None)
                    else:
                        assert line.startswith("positions: ")
                        positions = line.split(":")[1].split()
                        assert len(positions) == s
                        for p in positions:
                            j = p.find(".")
                            if precision == 0:
                                assert j == -1
                            else:
                                assert precision == len(p) - j - 1
                        values = list(map(float, positions))
                        assert values == sorted(values)
                        for position in values:
                            assert position >= 0.0
                            assert position <= 1.0
                        line = next(f, None)
                        sequences_found = 0
                        while line is not None and line[0] in "01":
                            sequences_found += 1
                            sequence = line.rstrip()
                            assert len(sequence) == s
                            line = next(f, None)
                        assert sequences_found == sample_size
            assert num_replicates == num_replicates_found

    def test_zero_recombination_rate(self):
        self.verify_output(
            sample_size=10,
            mutation_rate=1,
            num_loci=10,
            recombination_rate=0,
            num_replicates=2,
        )

    def test_invisible_recombinations(self):
        self.verify_output(
            sample_size=10,
            mutation_rate=0,
            num_loci=100,
            recombination_rate=1,
            num_replicates=1,
        )

    def test_num_replicates(self):
        for j in range(1, 10):
            self.verify_output(sample_size=10, mutation_rate=0, num_replicates=j)
            self.verify_output(sample_size=10, mutation_rate=10, num_replicates=j)
            self.verify_output(
                sample_size=10,
                mutation_rate=0,
                num_loci=10,
                recombination_rate=100,
                num_replicates=j,
            )
            self.verify_output(
                sample_size=10,
                mutation_rate=0,
                num_loci=1,
                recombination_rate=1,
                num_replicates=j,
            )
            self.verify_output(
                sample_size=10,
                mutation_rate=10,
                num_loci=1,
                recombination_rate=1,
                num_replicates=j,
            )
            self.verify_output(
                sample_size=10,
                mutation_rate=10,
                num_loci=10,
                recombination_rate=10,
                num_replicates=j,
            )

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
                sample_size=n, num_loci=10, recombination_rate=10, print_trees=True
            )
            self.verify_output(
                sample_size=n, num_loci=100, recombination_rate=10, print_trees=True
            )

    def test_seeds_output(self):
        self.verify_output(random_seeds=None)
        self.verify_output(random_seeds=(2, 3, 4))

    def test_correct_streams(self):
        args = "15 1 -r 0 1.0 -eG 1.0 5.25 -eG 2.0 10 -G 4 -eN 3.0 1.0 -T"
        stdout, stderr = capture_output(cli.mspms_main, args.split())
        assert len(stderr) == 0
        # We've already tested the output pretty thoroughly above so a
        # simple test is fine here.
        assert len(stdout.splitlines()) == 5

    def test_seed_equivalence(self):
        sample_size = 10
        mutation_rate = 10
        # Run without seeds to get automatically generated seeds
        sr = cli.SimulationRunner([sample_size], mutation_rate=mutation_rate)
        with tempfile.TemporaryFile("w+") as f:
            sr.run(f)
            f.seek(0)
            output1 = f.read()
        # Get the seeds
        seeds = list(map(int, output1.splitlines()[1].split()))
        # Run with the same seeds to get the same output.
        sr = cli.SimulationRunner(
            [sample_size],
            mutation_rate=mutation_rate,
            ms_random_seeds=seeds,
        )
        with tempfile.TemporaryFile("w+") as f:
            sr.run(f)
            f.seek(0)
            output2 = f.read()
        assert output1 == output2


class TestMspArgumentParser:
    """
    Tests for the argument parsers in msp.
    """

    def test_simulate_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args([cmd, "10", "out.trees"])
        assert args.sample_size == 10
        assert args.tree_sequence == "out.trees"
        assert args.recombination_rate == 0.0
        assert args.mutation_rate == 0.0
        assert args.length == 1
        assert args.effective_population_size == 1
        assert args.random_seed is None
        assert not args.compress

    def test_simulate_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args(
            [
                cmd,
                "100",
                "out2.trees",
                "-L",
                "1e3",
                "-r",
                "5",
                "-u",
                "2",
                "-s",
                "1234",
                "-z",
                "-N",
                "11",
            ]
        )
        assert args.sample_size == 100
        assert args.tree_sequence == "out2.trees"
        assert args.recombination_rate == 5
        assert args.length == 1000
        assert args.random_seed == 1234
        assert args.compress
        assert args.effective_population_size == 11

    def test_simulate_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args(
            [
                cmd,
                "1000",
                "out3.trees",
                "--length",
                "1e4",
                "--recombination-rate",
                "6",
                "--effective-population-size",
                "1e5",
                "--mutation-rate",
                "1",
                "--random-seed",
                "123",
                "--compress",
            ]
        )
        assert args.sample_size == 1000
        assert args.tree_sequence == "out3.trees"
        assert args.recombination_rate == 6
        assert args.length == 10000
        assert args.effective_population_size == 10 ** 5
        assert args.random_seed == 123
        assert args.compress


class TestMspSimulateOutput(unittest.TestCase):
    """
    Tests the output of msp to ensure it's correct.
    """

    def setUp(self):
        fd, self._tree_sequence = tempfile.mkstemp(prefix="msp_cli", suffix=".trees")
        os.close(fd)

    def tearDown(self):
        os.unlink(self._tree_sequence)

    def test_run_defaults(self):
        cmd = "simulate"
        sample_size = 10
        stdout, stderr = capture_output(
            cli.msp_main, [cmd, str(sample_size), self._tree_sequence]
        )
        assert len(stderr) == 0
        assert len(stdout) == 0

        tree_sequence = tskit.load(self._tree_sequence)
        assert tree_sequence.get_sample_size() == sample_size
        assert tree_sequence.get_sequence_length() == 1
        assert tree_sequence.get_num_mutations() == 0

    def test_simulate_short_args(self):
        cmd = "simulate"
        stdout, stderr = capture_output(
            cli.msp_main,
            [cmd, "100", self._tree_sequence, "-L", "1e2", "-r", "5", "-u", "2"],
        )
        tree_sequence = tskit.load(self._tree_sequence)
        assert len(stderr) == 0
        assert len(stdout) == 0
        assert tree_sequence.get_sample_size() == 100
        assert tree_sequence.get_sequence_length() == 100
        assert tree_sequence.get_num_mutations() > 0

    def test_compress_warns(self):
        cmd = "simulate"
        with pytest.warns(UserWarning):
            capture_output(cli.msp_main, [cmd, "10", self._tree_sequence, "--compress"])
        tree_sequence = tskit.load(self._tree_sequence)
        assert tree_sequence.get_sample_size() == 10


class TestMspConversionOutput(unittest.TestCase):
    """
    Tests the output of msp to ensure it's correct.
    """

    @classmethod
    def setUpClass(cls):
        cls._tree_sequence = msprime.simulate(
            10, length=10, recombination_rate=10, mutation_rate=10, random_seed=1
        )
        fd, cls._tree_sequence_file = tempfile.mkstemp(
            prefix="msp_cli", suffix=".trees"
        )
        os.close(fd)
        cls._tree_sequence.dump(cls._tree_sequence_file)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls._tree_sequence_file)

    def test_mutate_keep(self):
        fd, out_tree_sequence_file = tempfile.mkstemp(
            prefix="msp_cli_mutate", suffix=".trees"
        )
        os.close(fd)

        cmd = "mutate"
        mutation_rate = 10
        seed = 1
        stdout, stderr = capture_output(
            cli.msp_main,
            [
                cmd,
                self._tree_sequence_file,
                out_tree_sequence_file,
                "-u",
                str(mutation_rate),
                "-s",
                str(seed),
                "--keep",
            ],
        )
        assert len(stderr) == 0

        previous_ts = tskit.load(self._tree_sequence_file)
        tree_sequence = tskit.load(out_tree_sequence_file)

        assert tree_sequence.get_num_mutations() > previous_ts.get_num_mutations()

    def test_mutate_discrete_start_end_time(self):
        fd, out_tree_sequence_file = tempfile.mkstemp(
            prefix="msp_cli_mutate", suffix=".trees"
        )
        os.close(fd)

        cmd = "mutate"
        mutation_rate = 10
        seed = 1
        stdout, stderr = capture_output(
            cli.msp_main,
            [
                cmd,
                self._tree_sequence_file,
                out_tree_sequence_file,
                "-u",
                str(mutation_rate),
                "-s",
                str(seed),
                "--discrete",
                "--start-time",
                "0",
                "--end-time",
                "2",
            ],
        )
        assert len(stderr) == 0

        tree_sequence = tskit.load(out_tree_sequence_file)
        tables = tree_sequence.dump_tables()

        assert max(tables.nodes.time[tables.mutations.node]) <= 2
        assert min(tables.nodes.time[tables.mutations.node]) >= 0
        assert set(tables.sites.position) <= set(range(10))
        assert tree_sequence.get_num_mutations() > 0
