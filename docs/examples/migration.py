import numpy as np

import msprime


def migration_example(num_replicates=10 ** 5):
    # M is the overall symmetric migration rate, and d is the number
    # of subpopulations.
    M = 0.2
    d = 3
    # We rescale m into per-generation values for msprime.
    m = M / (4 * (d - 1))
    # Allocate the initial sample. Because we are interested in the
    # between subpopulation coalescence times, we choose one sample each
    # from the first two subpopulations.
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=1),
        msprime.PopulationConfiguration(sample_size=1),
        msprime.PopulationConfiguration(sample_size=0),
    ]
    # Now we set up the migration matrix. Since this is a symmetric
    # island model, we have the same rate of migration between all
    # pairs of subpopulations. Diagonal elements must be zero.
    migration_matrix = [[0, m, m], [m, 0, m], [m, m, 0]]
    # We pass these values to the simulate function, and ask it
    # to run the required number of replicates.
    replicates = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        num_replicates=num_replicates,
    )
    # And then iterate over these replicates
    T = np.zeros(num_replicates)
    for i, tree_sequence in enumerate(replicates):
        tree = tree_sequence.first()
        # Convert the TMRCA to coalecent units.
        T[i] = tree.get_time(tree.get_root()) / 4
    # Finally, calculate the analytical expectation and print
    # out the results
    analytical = d / 2 + (d - 1) / (2 * M)
    print("Observed  =", np.mean(T))
    print("Predicted =", analytical)
