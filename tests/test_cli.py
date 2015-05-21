"""
Test cases for the command line interfaces to msprime
"""
from __future__ import print_function
from __future__ import division

import os
import random
import sys
import tempfile
import unittest

import tests
import msprime.cli as cli

class TestRandomSeeds(unittest.TestCase):
    """
    Test the random seed generation for the ms compatability layer.
    """
    def test_within_range(self):
        num_random_tests = 100
        max_seed = 2**16 - 1
        generated_seeds = {}
        for j in range(100):
            seeds = [random.randint(1, max_seed) for k in range(3)]
            python_seed, ms_seeds = cli.get_seeds(seeds)
            self.assertEqual(ms_seeds, seeds)
            self.assertGreater(python_seed, 0)
            generated_seeds[tuple(seeds)] = python_seed
            # Make sure it's deterministic
            python_seed2, ms_seeds2 = cli.get_seeds(seeds)
            self.assertEqual(ms_seeds, ms_seeds2)
            self.assertEqual(python_seed, python_seed2)
        self.assertEqual(len(generated_seeds),
                len(set(generated_seeds.keys())))


class TestMspmsOutput(unittest.TestCase):
    """
    Tests the output of the ms compatible CLI.
    """
    def verify_output(self,
            sample_size, num_replicates, print_trees=False, num_loci=None,
            recombination_rate=None, mutation_rate=None):
        """
        Runs the UI for the specified parameters, and parses the output
        to ensure it's consistent.
        """
        cmdLine = [str(sample_size), str(num_replicates)]
        if print_trees:
            cmdLine.append("-T")

        stdout = sys.stdout
        try:
            with tempfile.TemporaryFile() as f:
                sys.stdout = f
                cli.mspms_main(cmdLine)
                f.seek(0)
                output = f.readlines()
        finally:
            sys.stdout = stdout
        # TODO write the tests...


    def test_tree_output(self):
        self.verify_output(10, 1, True)
