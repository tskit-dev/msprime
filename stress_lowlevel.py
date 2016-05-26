"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
from __future__ import print_function
from __future__ import division

import unittest
import random

import tests.test_lowlevel as test_lowlevel
import tests.test_highlevel as test_highlevel

import msprime


def main():
    while True:
        # We don't want any random variation in the amount of memory
        # used from test-to-test.
        random.seed(1)
        testloader = unittest.TestLoader()
        test_lowlevel.enable_h5py_tests = False
        suite = testloader.loadTestsFromModule(test_lowlevel)
        l = testloader.loadTestsFromModule(test_highlevel)
        suite.addTests(l)
        unittest.TextTestRunner(verbosity=0).run(suite)

        # Run a large number of replicate simulations to make sure we're
        # not leaking memory here.
        num_replicates = 10000
        j = 0
        replicates = msprime.simulate(
            sample_size=100,
            recombination_rate=10,
            num_replicates=num_replicates)
        for ts in replicates:
            j += 1
        assert j == num_replicates


if __name__ == "__main__":
    main()
