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

import random
import sys
import tempfile
import unittest

import msprime.cli as cli

# We're forced to do this because dendropy doesn't support Python 3.
_dendropy_available = True
try:
    import dendropy
except ImportError:
    _dendropy_available = False


class TestRandomSeeds(unittest.TestCase):
    """
    Test the random seed generation for the ms compatability layer.
    """
    def test_within_range(self):
        num_random_tests = 100
        max_seed = 2**16 - 1
        generated_seeds = {}
        for j in range(num_random_tests):
            seeds = [random.randint(1, max_seed) for k in range(3)]
            python_seed, ms_seeds = cli.get_seeds(seeds)
            self.assertEqual(ms_seeds, seeds)
            self.assertGreater(python_seed, 0)
            generated_seeds[tuple(seeds)] = python_seed
            # Make sure it's deterministic
            python_seed2, ms_seeds2 = cli.get_seeds(seeds)
            self.assertEqual(ms_seeds, ms_seeds2)
            self.assertEqual(python_seed, python_seed2)
        self.assertEqual(
            len(generated_seeds), len(set(generated_seeds.keys())))


class TestArgumentParser(unittest.TestCase):
    """
    Tests the parser to ensure it works correctly and is ms compatible.
    """

    def test_msdoc_examples(self):
        parser = cli.get_mspms_parser()

        args = parser.parse_args(["4", "2", "-t", "5.0"])
        self.assertEqual(args.sample_size, 4)
        self.assertEqual(args.num_replicates, 2)
        self.assertEqual(args.mutation_rate, 5.0)
        self.assertEqual(args.trees, False)

        args = parser.parse_args(["4", "2", "-T"])
        self.assertEqual(args.sample_size, 4)
        self.assertEqual(args.num_replicates, 2)
        self.assertEqual(args.mutation_rate, 0.0)
        self.assertEqual(args.trees, True)

        args = parser.parse_args("15 1000 -t 10.04 -r 100.0 2501".split())
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 10.04)
        self.assertEqual(args.trees, False)
        self.assertEqual(args.recombination, [100, 2501])

        args = parser.parse_args(
            "15 1000 -t 2.0 -eN 1.0 .1 -eN 2.0 4.0".split())
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 2.0)
        self.assertEqual(args.size_event, [[1.0, 0.1], [2.0, 4.0]])

        args = parser.parse_args(
            "15 1000 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5".split())
        self.assertEqual(args.sample_size, 15)
        self.assertEqual(args.num_replicates, 1000)
        self.assertEqual(args.mutation_rate, 6.4)
        self.assertEqual(args.growth_rate, 6.93)
        self.assertEqual(args.growth_event, [[0.2, 0.0]])
        self.assertEqual(args.size_event, [[0.3, 0.5]])

    def test_positional_arguments(self):
        parser = cli.get_mspms_parser()
        args = parser.parse_args(["40", "20"])
        self.assertEqual(args.sample_size, 40)
        self.assertEqual(args.num_replicates, 20)
        # TODO test errors here.

    def test_mutations(self):
        parser = cli.get_mspms_parser()
        args = parser.parse_args(["40", "20"])
        self.assertEqual(args.mutation_rate, 0.0)
        args = parser.parse_args(["40", "20", "-t", "10"])
        self.assertEqual(args.mutation_rate, 10.0)
        args = parser.parse_args(["40", "20", "--mutation-rate=10"])
        self.assertEqual(args.mutation_rate, 10.0)
        # TODO test errors

    def test_trees(self):
        parser = cli.get_mspms_parser()
        args = parser.parse_args(["40", "20"])
        self.assertEqual(args.trees, False)
        args = parser.parse_args(["40", "20", "-T"])
        self.assertEqual(args.trees, True)
        args = parser.parse_args(["40", "20", "--trees"])
        self.assertEqual(args.trees, True)
        # TODO test errors

    def test_size_events(self):
        parser = cli.get_mspms_parser()
        args = parser.parse_args(["40", "20"])
        self.assertEqual(args.size_event, [])
        args = parser.parse_args("10 1 -eN 2.0 0.5".split())
        self.assertEqual(args.size_event, [[2.0, 0.5]])
        args = parser.parse_args("10 1 -eN 2.0 0.5 -eN 1.0 5.0".split())
        self.assertEqual(args.size_event, [[2.0, 0.5], [1.0, 5.0]])
        # TODO test errors

    def test_growth_rates(self):
        parser = cli.get_mspms_parser()
        args = parser.parse_args(["40", "20"])
        self.assertEqual(args.growth_rate, None)
        self.assertEqual(args.growth_event, [])
        args = parser.parse_args("15 1000 -G 5.25".split())
        self.assertEqual(args.growth_rate, 5.25)
        self.assertEqual(args.growth_event, [])
        args = parser.parse_args("15 1000 -eG 1.0 5.25".split())
        self.assertEqual(args.growth_rate, None)
        self.assertEqual(args.growth_event, [[1.0, 5.25]])
        args = parser.parse_args("15 1000 -eG 1.0 5.25 -eG 2.0 10".split())
        self.assertEqual(args.growth_rate, None)
        self.assertEqual(args.growth_event, [[1.0, 5.25], [2.0, 10.0]])
        args = parser.parse_args(
            "15 1000 -eG 1.0 5.25 -eG 2.0 10 -G 4".split())
        self.assertEqual(args.growth_rate, 4.0)
        self.assertEqual(args.growth_event, [[1.0, 5.25], [2.0, 10.0]])
        # TODO test errors


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
            max_memory="16M", precision=3, population_models=[],
            random_seeds=[1, 2, 3]):
        """
        Runs the UI for the specified parameters, and parses the output
        to ensure it's consistent.
        """
        sr = cli.SimulationRunner(
            sample_size=sample_size, num_loci=num_loci,
            recombination_rate=recombination_rate,
            num_replicates=num_replicates, mutation_rate=mutation_rate,
            print_trees=print_trees, precision=precision,
            population_models=population_models,
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
                        line = next(f, None)
                        sequences_found = 0
                        while line is not None and line[0] in "01":
                            sequences_found += 1
                            sequence = line.rstrip()
                            self.assertEqual(len(sequence), s)
                            line = next(f, None)
                        self.assertEqual(sequences_found, sample_size)
            self.assertEqual(num_replicates, num_replicates_found)

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
                sample_size=10, mutation_rate=0, num_loci=10,
                recombination_rate=0, num_replicates=j)
            self.verify_output(
                sample_size=10, mutation_rate=10, num_loci=10,
                recombination_rate=0, num_replicates=j)
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
