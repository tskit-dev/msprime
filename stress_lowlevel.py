"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
from __future__ import print_function
from __future__ import division

import unittest
import random

import tests.test_demography as test_demography
import tests.test_highlevel as test_highlevel
import tests.test_lowlevel as test_lowlevel
import tests.test_vcf as test_vcf

# Currently skipping some CLI tests so won't work here.
# import tests.test_cli as test_cli

import msprime
import _msprime


def main():
    while True:
        # We don't want any random variation in the amount of memory
        # used from test-to-test.
        random.seed(1)
        testloader = unittest.TestLoader()
        suite = testloader.loadTestsFromModule(test_lowlevel)
        modules = [
            test_highlevel, test_demography, test_vcf]
        for mod in modules:
            l = testloader.loadTestsFromModule(mod)
            suite.addTests(l)
        unittest.TextTestRunner(verbosity=0).run(suite)

        # When we have lots of error conditions we get a small memory
        # leak in HDF5. To counter this we call H5close after every set of
        # tests. This means that we cannot call the hdf5 tests in this loop
        # though, because h5py does not like h5close being called while it
        # is still in use.
        _msprime.h5close()

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
