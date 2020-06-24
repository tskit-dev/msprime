import daiquiri

import msprime


def logging_info_example():
    daiquiri.setup(level="INFO")
    msprime.simulate(10, Ne=1000, model=["dtwf", (100, "hudson")])
