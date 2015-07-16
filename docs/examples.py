"""
The examples used in the tutorial section.
"""
from __future__ import print_function
from __future__ import division


import msprime

def single_locus_example():
    tree = msprime.simulate_tree(5, random_seed=1)
    print(tree)
    tree.draw("_static/simple-tree.svg")
    u = 1
    while u != 0:
        print("node {}: time = {}".format(u, tree.get_time(u)))
        u = tree.get_parent(u)

    print(tree.get_branch_length(7))


if __name__ == "__main__":
    single_locus_example()
