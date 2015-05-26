"""
Example in which we reproduce the simulations in the GQT paper,
"Efficient compression and analysis of large genetic variation datasets"
by Layer et al.
"""
from __future__ import print_function
from __future__ import division

import time

import msprime


def main():
    before = time.clock()
    # Run the actual simulations
    tree_sequence = msprime.simulate(
        sample_size=10**5,
        num_loci=100 * 10**6,
        scaled_recombination_rate=0.001,
        scaled_mutation_rate=0.001,
        max_memory="5G",
        random_seed=1  # Arbitrary - make this reproducible.
    )
    duration = time.clock() - before
    print("Simulated 100k genomes in {0:.3f} seconds.".format(duration))

    # Write the results to file, which is small and can be quickly reloaded
    # to avoid the cost of re-running the simulation. We can reload the
    # file in a few seconds using msprime.load(filename).
    tree_sequence.dump("large-example.hdf5")

    # Now write the haplotypes to a file.
    # WARNING! This takes a lot of memory (>100G), so make sure you don't
    # crash your server. This memory requirement will be drastically reduced
    # in future versions.
    before = time.clock()
    with open("large-example-haplotypes.txt", "w") as f:
        for h in tree_sequence.haplotypes():
            print(h, file=f)
    duration = time.clock() - before
    print("Wrote 100k haplotypes to file in {0:.3f} seconds".format(duration))


if __name__ == "__main__":
    main()
