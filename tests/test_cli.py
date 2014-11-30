"""
Test cases for the command line interfaces to msprime
"""
from __future__ import print_function
from __future__ import division

import os
import random
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
        self.assertEqual(len(generated_seeds),
                len(set(generated_seeds.keys())))
