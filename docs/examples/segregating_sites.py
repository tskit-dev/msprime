import numpy as np

import msprime


def segregating_sites(n, theta, num_replicates):
    S = np.zeros(num_replicates)
    replicates = msprime.simulate(
        sample_size=n, mutation_rate=theta / 4, num_replicates=num_replicates
    )
    for j, tree_sequence in enumerate(replicates):
        S[j] = tree_sequence.get_num_mutations()
    S_mean_a = np.sum(1 / np.arange(1, n)) * theta
    S_var_a = theta * np.sum(1 / np.arange(1, n)) + theta ** 2 * np.sum(
        1 / np.arange(1, n) ** 2
    )
    print("              mean              variance")
    print(f"Observed      {np.mean(S):.5f}\t\t{np.var(S):.5f}")
    print(f"Analytical    {S_mean_a:.5f}\t\t{S_var_a:.5f}")
