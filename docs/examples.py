"""
The examples used in the tutorial section.
"""
from __future__ import print_function
from __future__ import division


import msprime

def single_locus_example():
    tree = msprime.simulate_tree(10)
    print(tree)


if __name__ == "__main__":
    single_locus_example()
