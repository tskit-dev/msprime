"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
from __future__ import print_function
from __future__ import division

import argparse
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
import tests.test_tables as test_tables
import tests.test_topology as test_topology

import _msprime


def main():
    modules = {
        "demography": test_demography,
        "highlevel": test_highlevel,
        "lowlevel": test_lowlevel,
        "vcf": test_vcf,
        "threads": test_threads,
        "stats": test_stats,
        "tables": test_tables,
        "topology": test_topology,
    }
    parser = argparse.ArgumentParser(
        description="Run tests in a loop to stress low-level interface")
    parser.add_argument(
        "-m", "--module", help="Run tests only on this module",
        choices=list(modules.keys()))
    args = parser.parse_args()
    test_modules = list(modules.values())
    if args.module is not None:
        test_modules = [modules[args.module]]

    print("iter\ttests\terr\tfail\tskip\tRSS\tmin\tmax\tmax@iter")
    max_rss = 0
    max_rss_iter = 0
    min_rss = 1e100
    iteration = 0
    last_print = time.time()
    devnull = open(os.devnull, 'w')
    while True:
        # We don't want any random variation in the amount of memory
        # used from test-to-test.
        random.seed(1)
        testloader = unittest.TestLoader()
        suite = testloader.loadTestsFromModule(test_modules[0])
        for mod in test_modules[1:]:
            suite.addTests(testloader.loadTestsFromModule(mod))
        runner = unittest.TextTestRunner(verbosity=0, stream=devnull)
        result = runner.run(suite)

        # When we have lots of error conditions we get a small memory
        # leak in HDF5. To counter this we call H5close after every set of
        # tests. This means that we cannot call the hdf5 tests in this loop
        # though, because h5py does not like h5close being called while it
        # is still in use.
        _msprime.h5close()

        rusage = resource.getrusage(resource.RUSAGE_SELF)
        if max_rss < rusage.ru_maxrss:
            max_rss = rusage.ru_maxrss
            max_rss_iter = iteration
        if min_rss > rusage.ru_maxrss:
            min_rss = rusage.ru_maxrss

        # We don't want to flood stdout, so we rate-limit to 1 per second.
        if time.time() - last_print > 1:
            print(
                iteration, result.testsRun, len(result.failures),
                len(result.errors), len(result.skipped),
                rusage.ru_maxrss, min_rss,  max_rss, max_rss_iter,
                sep="\t", end="\r")
            last_print = time.time()
            sys.stdout.flush()

        iteration += 1


if __name__ == "__main__":
    main()
