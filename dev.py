"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import msprime


def haplotype_example():
    tree_sequence = msprime.simulate(
        sample_size=10, num_loci=1000, scaled_recombination_rate=0.1,
        scaled_mutation_rate=0.01, random_seed=1)
    for h in tree_sequence.haplotypes():
        print(h)


def dump_example():
    tree_sequence = msprime.simulate(
        sample_size=10, num_loci=1000, scaled_recombination_rate=0.1,
        scaled_mutation_rate=0.01, random_seed=1)
    haplotypes = list(tree_sequence.haplotypes())
    tree_sequence.dump("example.hdf5")
    # Now, load another tree sequence instance from this file
    other_tree_sequence = msprime.load("example.hdf5")
    other_haplotypes = list(other_tree_sequence.haplotypes())
    assert haplotypes == other_haplotypes


if __name__ == "__main__":
    # haplotype_example()
    dump_example()
