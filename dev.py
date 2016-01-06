"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import os
import math
import itertools

import msprime

def mutations():
    recomb_rates = [(10, 0.05), (20, 0.1), (30, 0), (40, 0.05)]
    for x, rate in recomb_rates:
        print(x, rate)
    max_rate = max(rate for _, rate in recomb_rates)
    print("max_ rate = ", max_rate)

    tree_sequence = msprime.simulate(10, 100, max_rate, random_seed=1)
    for tree in tree_sequence.trees():
        print(tree.get_interval())

if __name__ == "__main__":
    mutations()

