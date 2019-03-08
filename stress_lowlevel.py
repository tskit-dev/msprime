"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
import argparse
import unittest
import random
import resource
import os
import sys
import time
import logging

import tests.test_demography as test_demography
import tests.test_highlevel as test_highlevel
import tests.test_lowlevel as test_lowlevel
import tests.test_dict_encoding as test_dict_encoding


def main():
    modules = {
        "demography": test_demography,
        "highlevel": test_highlevel,
        "lowlevel": test_lowlevel,
        "dict_encoding": test_dict_encoding,
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

    # Need to do this to silence the errors from the file_format tests.
    logging.basicConfig(level=logging.ERROR)

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
