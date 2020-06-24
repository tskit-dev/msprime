import daiquiri

import msprime


def logging_debug_example():
    daiquiri.setup(level="DEBUG")
    msprime.simulate(
        10 ** 5, Ne=10000, recombination_rate=2e-8, length=1e6, random_seed=32
    )
