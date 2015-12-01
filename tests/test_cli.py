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


class TestMspmsArgumentParser(unittest.TestCase):
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

    def test_correct_streams(self):
        args = "15 1 -r 0 1.0 -eG 1.0 5.25 -eG 2.0 10 -G 4 -eN 3.0 1.0 -T"
        stdout, stderr = capture_output(cli.mspms_main, args.split())
        self.assertEqual(len(stderr), 0)
        # We've already tested the output pretty thoroughly above so a
        # simple test is fine here.
        self.assertEqual(len(stdout.splitlines()), 5)


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
        self.assertEqual(args.num_loci, 1)
        self.assertEqual(args.random_seed, None)
        self.assertEqual(args.max_memory, "1G")
        self.assertEqual(args.compress, False)

    def test_simulate_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args([
            cmd, "100", "out2.hdf5", "-m", "1e3", "-r", "5", "-u", "2",
            "-s", "1234", "-M", "2G", "-z"])
        self.assertEqual(args.sample_size, 100)
        self.assertEqual(args.history_file, "out2.hdf5")
        self.assertEqual(args.recombination_rate, 5)
        self.assertEqual(args.num_loci, 1000)
        self.assertEqual(args.random_seed, 1234)
        self.assertEqual(args.max_memory, "2G")
        self.assertEqual(args.compress, True)

    def test_simulate_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "simulate"
        args = parser.parse_args([
            cmd, "1000", "out3.hdf5",
            "--num-loci", "1e4",
            "--recombination-rate", "6",
            "--mutation-rate", "1",
            "--random-seed", "123",
            "--max-memory", "20M",
            "--compress"])
        self.assertEqual(args.sample_size, 1000)
        self.assertEqual(args.history_file, "out3.hdf5")
        self.assertEqual(args.recombination_rate, 6)
        self.assertEqual(args.num_loci, 10000)
        self.assertEqual(args.random_seed, 123)
        self.assertEqual(args.max_memory, "20M")
        self.assertEqual(args.compress, True)

    def test_records_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "records"
        history_file = "test.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, False)

    def test_records_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "records"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "-H"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)

    def test_records_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "records"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "--header"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)

    def test_mutations_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "mutations"
        history_file = "test.hdf5"
        args = parser.parse_args([cmd, history_file])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, False)

    def test_mutations_short_args(self):
        parser = cli.get_msp_parser()
        cmd = "mutations"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "-H"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)

    def test_mutations_long_args(self):
        parser = cli.get_msp_parser()
        cmd = "mutations"
        history_file = "test.hdf5"
        args = parser.parse_args([
            cmd, history_file, "--header"])
        self.assertEqual(args.history_file, history_file)
        self.assertEqual(args.header, True)

    def test_haplotypes_default_values(self):
        parser = cli.get_msp_parser()
        cmd = "haplotypes"
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
        self.assertEqual(tree_sequence.get_num_loci(), 1)
        self.assertEqual(tree_sequence.get_num_mutations(), 0)

    def test_simulate_short_args(self):
        cmd = "simulate"
        stdout, stdearr = capture_output(cli.msp_main, [
            cmd, "100", self._history_file, "-m", "1e2", "-r", "5", "-u", "2"])
        tree_sequence = msprime.load(self._history_file)
        self.assertEqual(tree_sequence.get_sample_size(), 100)
        self.assertEqual(tree_sequence.get_num_loci(), 100)
        self.assertGreater(tree_sequence.get_num_mutations(), 0)


class TestMspConversionOutput(unittest.TestCase):
    """
    Tests the output of msp to ensure it's correct.
    """
    @classmethod
    def setUpClass(cls):
        cls._tree_sequence = msprime.simulate(
            10, 10, scaled_recombination_rate=10,
            scaled_mutation_rate=10, random_seed=1)
        fd, cls._history_file = tempfile.mkstemp(
            prefix="msp_cli", suffix=".hdf5")
        os.close(fd)
        cls._tree_sequence.dump(cls._history_file)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls._history_file)

    def verify_records(self, output_records):
        records = list(self._tree_sequence.records())
        self.assertEqual(len(records), len(output_records))
        for (l, r, u, c, t), line in zip(records, output_records):
            splits = line.split()
            self.assertEqual(l, int(splits[0]))
            self.assertEqual(r, int(splits[1]))
            self.assertEqual(u, int(splits[2]))
            self.assertEqual(c[0], int(splits[3]))
            self.assertEqual(c[1], int(splits[4]))
            self.assertAlmostEqual(t, float(splits[5]))

    def test_records(self):
        cmd = "records"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        output_records = stdout.splitlines()
        self.verify_records(output_records)
        # check the header.
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file, "-H"])
        self.assertEqual(len(stderr), 0)
        output_records = stdout.splitlines()
        self.assertEqual(
            list(output_records[0].split()),
            ["l", "r", "u", "c1", "c2", "t"])
        self.verify_records(output_records[1:])

    def verify_mutations(self, output_mutations):
        mutations = list(self._tree_sequence.mutations())
        self.assertEqual(len(mutations), len(output_mutations))
        for (x, u), line in zip(mutations, output_mutations):
            splits = line.split()
            self.assertAlmostEqual(x, float(splits[0]))
            self.assertEqual(u, int(splits[1]))

    def test_mutations(self):
        cmd = "mutations"
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file])
        self.assertEqual(len(stderr), 0)
        output_mutations = stdout.splitlines()
        self.verify_mutations(output_mutations)
        # check the header.
        stdout, stderr = capture_output(cli.msp_main, [
            cmd, self._history_file, "-H"])
        self.assertEqual(len(stderr), 0)
        output_mutations = stdout.splitlines()
        self.assertEqual(
            list(output_mutations[0].split()), ["x", "u"])
        self.verify_mutations(output_mutations[1:])

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
        m = self._tree_sequence.get_num_loci()
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
