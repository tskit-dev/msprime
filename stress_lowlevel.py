"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
from __future__ import print_function
from __future__ import division

import unittest
import random

import tests.test_demography as test_demography
import tests.test_hdf5 as test_hdf5
import tests.test_highlevel as test_highlevel
import tests.test_lowlevel as test_lowlevel
import tests.test_vcf as test_vcf

# Currently skipping some CLI tests so won't work here.
# import tests.test_cli as test_cli

import msprime


def main():
    while True:
        # We don't want any random variation in the amount of memory
        # used from test-to-test.
        random.seed(1)
        testloader = unittest.TestLoader()
        suite = testloader.loadTestsFromModule(test_lowlevel)
        # suite = testloader.loadTestsFromModule(test_hdf5)
        modules = [
            test_highlevel, test_demography, test_vcf, test_hdf5]
        for mod in modules:
            l = testloader.loadTestsFromModule(mod)
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
