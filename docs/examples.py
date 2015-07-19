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

def mutations_example():
    tree_sequence = msprime.simulate(
        5, num_loci=10, scaled_recombination_rate=0.1,
        scaled_mutation_rate=0.2, random_seed=19)
    print("Total mutations = ", tree_sequence.get_num_mutations())
    j = 0
    for tree in tree_sequence.sparse_trees():
        print(tree.get_interval(), list(tree.mutations()), sep="\t")
        tree.draw("_static/mutations-tree-sequence-{}.svg".format(j))
        j += 1

    haplotypes = [None] + [
        ['0' for _ in range(tree_sequence.get_num_mutations())]
        for _ in range(tree_sequence.get_sample_size())]
    site = 0
    for tree in tree_sequence.sparse_trees():
        for _, node in tree.mutations():
            for u in tree.leaves(node):
                haplotypes[u][site] = '1'
            site += 1
    for j in range(1, tree_sequence.get_sample_size() + 1):
        print(j, "".join(haplotypes[j]), sep="\t")


if __name__ == "__main__":
    single_locus_example()
    multi_locus_example()
    mutations_example()
