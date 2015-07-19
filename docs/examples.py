"""
The examples used in the tutorial section.
"""
from __future__ import print_function
from __future__ import division


import msprime

def single_locus_example():
    tree = msprime.simulate_tree(5, random_seed=1)
    print(tree)
    tree.draw("_static/simple-tree.svg", show_times=True)
    u = 1
    while u != 0:
        print("node {}: time = {}".format(u, tree.get_time(u)))
        u = tree.get_parent(u)

    print(tree.get_branch_length(7))

def multi_locus_example():
    tree_sequence = msprime.simulate(
        5, num_loci=10, scaled_recombination_rate=0.1, random_seed=19)
    j = 0
    for tree in tree_sequence.sparse_trees():
        print(tree.get_interval(), str(tree), sep="\t")
        tree.draw("_static/simple-tree-sequence-{}.svg".format(j))
        j += 1

if __name__ == "__main__":
    single_locus_example()
    multi_locus_example()
