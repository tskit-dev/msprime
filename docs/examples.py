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


def migration_example():
    # M is the overall symmetric migration rate, and d is the number
    # of demes.
    M = 0.2
    d = 3
    # We rescale m into per-generation values for msprime.
    m = M / (4 * (d - 1))
    # Allocate the initial sample. Because we are interested in the
    # between deme coalescence times, we choose one sample each
    # from the first two demes.
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=1),
        msprime.PopulationConfiguration(sample_size=1),
        msprime.PopulationConfiguration(sample_size=0)]
    # Now we set up the migration matrix. Since this is a symmetric
    # island model, we have the same rate of migration between all
    # pairs of demes. Diagonal elements must be zero.
    migration_matrix = [
        [0, m, m],
        [m, 0, m],
        [m, m, 0]]
    # We pass these values to the simulate function, and ask it
    # to run the required number of replicates.
    num_replicates = 10000
    replicates = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        num_replicates=num_replicates)
    # And then iterate over these replicates
    T = np.zeros(num_replicates)
    for i, tree_sequence in enumerate(replicates):
        tree = next(tree_sequence.trees())
        T[i] = tree.get_time(tree.get_root())
    # Finally, calculate the analytical expectation and print
    # out the results
    analytical = d / 2 + (d - 1) / (2 * M)
    print("Observed  =", np.mean(T))
    print("Predicted =", analytical)

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
    structure_example()
