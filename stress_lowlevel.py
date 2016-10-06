"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
from __future__ import print_function
from __future__ import division

import unittest
import random
import resource
import os
import sys
import time

import tests.test_demography as test_demography
import tests.test_highlevel as test_highlevel
import tests.test_lowlevel as test_lowlevel
import tests.test_vcf as test_vcf
import tests.test_threads as test_threads
import tests.test_stats as test_stats

# import tests.test_cli as test_cli

import msprime
import _msprime


def main():
    max_rss = 0
    min_rss = 1e100
    devnull = open(os.devnull, 'w')
    testloader = unittest.TestLoader()
    suite = testloader.loadTestsFromModule(test_lowlevel)
    modules = [
        test_highlevel, test_demography, test_vcf]
    for mod in modules:
        l = testloader.loadTestsFromModule(mod)
        suite.addTests(l)
    print("iter\ttests\terr\tfail\tskip\tRSS\tmin\tmax")
    iteration = 0
    last_print = time.time()
    while True:
        # We don't want any random variation in the amount of memory
        # used from test-to-test.
        random.seed(1)
        runner = unittest.TextTestRunner(verbosity=0, stream=devnull)
        result = runner.run(suite)

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

        # When we have lots of error conditions we get a small memory
        # leak in HDF5. To counter this we call H5close after every set of
        # tests. This means that we cannot call the hdf5 tests in this loop
        # though, because h5py does not like h5close being called while it
        # is still in use.
        _msprime.h5close()

        rusage = resource.getrusage(resource.RUSAGE_SELF)
        if max_rss < rusage.ru_maxrss:
            max_rss = rusage.ru_maxrss
        if min_rss > rusage.ru_maxrss:
            min_rss = rusage.ru_maxrss

        # We don't want to flood stdout, so we rate-limit to 1 per second.
        if time.time() - last_print > 1:
            print(
                iteration, result.testsRun, len(result.failures),
                len(result.errors), len(result.skipped),
                rusage.ru_maxrss, min_rss,  max_rss,
                sep="\t", end="\r")
            last_print = time.time()
            sys.stdout.flush()

        iteration += 1


if __name__ == "__main__":
    main()
