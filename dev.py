"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import _msprime

def main():
    treefile = "tmp.dat"
    for j in range(1000):
        sim = _msprime.Simulator(sample_size=4000, random_seed=j,
                coalescence_record_filename=treefile,
                num_loci=1000, recombination_rate=0.1,
                max_memory=1024**3)
        print(sim.run())
        print(sim.get_num_loci())
        print(sim.get_recombination_rate())
        print(sim.get_time())
        print(sim.get_random_seed())


    # tv = _msprime.TreeViewer(treefile)
    # for length, pi, tau in tv:
    #     print(length, pi, tau)



if __name__ == "__main__":
    main()
