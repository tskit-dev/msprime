"""
The examples used in the tutorial section.
"""
from __future__ import print_function
from __future__ import division

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import math
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


def demography_example():
    generation_time = 25
    # The ancestral population size.
    N_A = 7300
    T_AF = 220e3 / (generation_time * 4 * N_A)
    T_B = 140e3 / (generation_time * 4 * N_A)
    T_EU_AS = 21.2e3 / (generation_time * 4 * N_A)
    N_AF = 12300 / N_A
    N_B = 2100 / N_A
    N_EU0 = 1000 / N_A
    N_AS0 = 510 / N_A
    r_EU = 0.004 * 4 * N_A
    r_AS = 0.0055 * 4 * N_A
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the popupulation
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [0,   m_AF_EU, 0],
        [m_AF_EU, 0,   0],
        [0, 0, 0, ],
    ]
    demographic_events = [
        msprime.MassMigrationEvent(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
    ]
    dp = msprime.DemographyPrinter(
        population_configurations, migration_matrix,
        demographic_events, Ne=N_A)
    dp.debug_history()



if __name__ == "__main__":
    # segregating_sites_example(10, 5, 10000)
    # single_locus_example()
    # multi_locus_example()
    # mutations_example()
    # structure_example()
    demography_example()
