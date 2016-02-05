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


if __name__ == "__main__":
    main()
