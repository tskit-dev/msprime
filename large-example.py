"""
Example in which we reproduce the simulations in the GQT paper,
"Efficient compression and analysis of large genetic variation datasets"
by Layer et al.
"""
import time

import msprime


def main():
    before = time.clock()
    # Run the actual simulations
    tree_sequence = msprime.simulate(
        sample_size=10**6,
        length=100 * 10**6,
        Ne=1e4,
        recombination_rate=1e-8,
        mutation_rate=1e-8,
        random_seed=1  # Arbitrary - make this reproducible.
    )
    duration = time.clock() - before
    print("Simulated 100k genomes in {:.3f} seconds.".format(duration))

    # Write the results to file, which is small and can be quickly reloaded
    # to avoid the cost of re-running the simulation. We can reload the
    # file in a few seconds using msprime.load(filename).
    tree_sequence.dump("tmp__NOBACKUP__/large-example_{}.hdf5".format(msprime.__version__))


if __name__ == "__main__":
    main()
