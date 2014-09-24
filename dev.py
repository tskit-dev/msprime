"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import json
import _msprime

def main():
    treefile = "tmp.dat"
    while True:
        sim = _msprime.Simulator(sample_size=4000, random_seed=1,
                tree_file_name=treefile,
                num_loci=10000, recombination_rate=0.1,
                max_memory=1024**3)
        print(sim.run())
        print(sim.get_num_loci())
        print(sim.get_recombination_rate())
        print(sim.get_time())
        print(sim.get_random_seed())

        tr = _msprime.TreeReader(treefile)
        print(tr.get_num_loci())
        print(tr.get_sample_size())
        print(tr.get_num_trees())
        print(tr.get_metadata())
        s = json.loads(tr.get_metadata())
        print(s)
        for j in range(tr.get_num_trees()):
            breakpoint, pi, tau = tr.get_tree(j)
            #print(j, breakpoint)



    # tv = _msprime.TreeViewer(treefile)
    # for length, pi, tau in tv:
    #     print(length, pi, tau)



if __name__ == "__main__":
    main()
