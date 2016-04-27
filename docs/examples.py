"""
The examples used in the tutorial section.
"""
from __future__ import print_function
from __future__ import division


import msprime
import numpy as np

def segregating_sites_example(n, theta, num_replicates):
    S = np.zeros(num_replicates)
    replicates = msprime.simulate(
        sample_size=n,
        mutation_rate=theta / 4,
        num_replicates=num_replicates)
    for j, tree_sequence in enumerate(replicates):
        S[j] = tree_sequence.get_num_mutations()
    S_mean_a = np.sum(1 / np.arange(1, n)) * theta
    S_var_a = (
        theta * np.sum(1 / np.arange(1, n)) +
        theta**2 * np.sum(1 / np.arange(1, n)**2))
    print("              mean              variance")
    print("Observed      {}\t\t{}".format(np.mean(S), np.var(S)))
    print("Analytical    {:.5f}\t\t{:.5f}".format(S_mean_a, S_var_a))


def single_locus_example():
    tree_sequence = msprime.simulate(5, random_seed=1)
    tree = next(tree_sequence.trees())
    print(tree)
    tree.draw("_static/simple-tree.svg", show_times=True)
    u = 0
    while u != msprime.NULL_NODE:
        print("node {}: time = {}".format(u, tree.get_time(u)))
        u = tree.get_parent(u)
    print(tree.get_branch_length(6))

def multi_locus_example():
    tree_sequence = msprime.simulate(
        sample_size=5, length=10, recombination_rate=0.02, random_seed=19)
    j = 0
    for tree in tree_sequence.trees():
        print(tree.get_interval(), str(tree), sep="\t")
        tree.draw("_static/simple-tree-sequence-{}.svg".format(j))
        j += 1

def mutations_example():
    tree_sequence = msprime.simulate(
        sample_size=5, length=10, recombination_rate=0.02,
        mutation_rate=0.02, random_seed=19)
    print("Total mutations = ", tree_sequence.get_num_mutations())
    j = 0
    for tree in tree_sequence.trees():
        print(tree.get_interval(), list(tree.mutations()), sep="\t")
        tree.draw("_static/mutations-tree-sequence-{}.svg".format(j))
        j += 1

    for tree in tree_sequence.trees():
        for position, node in tree.mutations():
            print("Mutation @ position {} has frequency {}".format(
                position, tree.get_num_leaves(node) / tree.get_sample_size()))


if __name__ == "__main__":
    # segregating_sites_example(10, 5, 10000)
    single_locus_example()
    multi_locus_example()
    mutations_example()
