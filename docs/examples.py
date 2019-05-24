"""
The examples used in the tutorial section.
"""
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import math
import msprime
import numpy as np
import scipy.stats
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.collections

import random
import tqdm
import threading


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
        msprime.PopulationConfiguration(sample_size=0)]
    # Now we set up the migration matrix. Since this is a symmetric
    # island model, we have the same rate of migration between all
    # pairs of subpopulations. Diagonal elements must be zero.
    migration_matrix = [
        [0, m, m],
        [m, 0, m],
        [m, m, 0]]
    # We pass these values to the simulate function, and ask it
    # to run the required number of replicates.
    num_replicates = 1e6
    replicates = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        num_replicates=num_replicates)
    # And then iterate over these replicates
    T = np.zeros(num_replicates)
    for i, tree_sequence in enumerate(replicates):
        tree = next(tree_sequence.trees())
        # Convert the TMRCA to coalecent units.
        T[i] = tree.get_time(tree.get_root()) / 4
    # Finally, calculate the analytical expectation and print
    # out the results
    analytical = d / 2 + (d - 1) / (2 * M)
    print("Observed  =", np.mean(T))
    print("Predicted =", analytical)

def single_locus_example():
    tree_sequence = msprime.simulate(sample_size=5, Ne=1000, random_seed=1)
    tree = next(tree_sequence.trees())
    print(tree)
    tree.draw("_static/simple-tree.svg", show_times=True)
    u = 0
    while u != msprime.NULL_NODE:
        print("node {}: time = {}".format(u, tree.get_time(u)))
        u = tree.get_parent(u)
    print(tree.get_branch_length(6))
    print(tree.get_total_branch_length())

def multi_locus_example():
    tree_sequence = msprime.simulate(
        sample_size=5, Ne=1000, length=1e4, recombination_rate=2e-8,
        random_seed=19)
    j = 0
    for tree in tree_sequence.trees():
        print(tree.get_interval(), str(tree), sep="\t")
        tree.draw("_static/simple-tree-sequence-{}.svg".format(j))
        j += 1

def mutations_example():

    tree_sequence = msprime.simulate(
        sample_size=5, Ne=1000, length=1e4, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=19)
    print("Total mutations = ", tree_sequence.get_num_mutations())
    j = 0
    for tree in tree_sequence.trees():
        print(tree.get_interval(), list(tree.mutations()), sep="\t")
        tree.draw("_static/mutations-tree-sequence-{}.svg".format(j))
        j += 1

    for tree in tree_sequence.trees():
        for mutation in tree.mutations():
            print("Mutation @ position {} has frequency {}".format(
                mutation.position,
                tree.get_num_leaves(mutation.node) / tree.get_sample_size()))


def out_of_africa():
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting population sizes based on the growth
    # rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the popupulation
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dp = msprime.DemographyDebugger(
        Ne=N_A,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dp.print_history()


def variable_recomb_example():
    infile = "../hapmap/genetic_map_GRCh37_chr22.txt"
    # Read in the recombination map using the read_hapmap method,
    recomb_map = msprime.RecombinationMap.read_hapmap(infile)

    # Now we get the positions and rates from the recombination
    # map and plot these using 500 bins.
    positions = np.array(recomb_map.get_positions()[1:])
    rates = np.array(recomb_map.get_rates()[1:])
    num_bins = 500
    v, bin_edges, _ = scipy.stats.binned_statistic(
        positions, rates, bins=num_bins)
    x = bin_edges[:-1][np.logical_not(np.isnan(v))]
    y = v[np.logical_not(np.isnan(v))]
    fig, ax1 = pyplot.subplots(figsize=(16, 6))
    ax1.plot(x, y, color="blue")
    ax1.set_ylabel("Recombination rate")
    ax1.set_xlabel("Chromosome position")

    # Now we run the simulation for this map. We assume Ne=10^4
    # and have a sample of 100 individuals
    tree_sequence = msprime.simulate(
        sample_size=100,
        Ne=10**4,
        recombination_map=recomb_map)
    # Now plot the density of breakpoints along the chromosome
    breakpoints = np.array(list(tree_sequence.breakpoints()))
    ax2 = ax1.twinx()
    v, bin_edges = np.histogram(breakpoints, num_bins, density=True)
    ax2.plot(bin_edges[:-1], v, color="green")
    ax2.set_ylabel("Breakpoint density")
    ax2.set_xlim(1.5e7, 5.3e7)
    fig.savefig("_static/hapmap_chr22.svg")


def ld_matrix_example():
    ts = msprime.simulate(100, recombination_rate=10, mutation_rate=20,
            random_seed=1)
    ld_calc = msprime.LdCalculator(ts)
    A = ld_calc.get_r2_matrix()
    # Now plot this matrix.
    x = A.shape[0] / pyplot.rcParams['savefig.dpi']
    x = max(x, pyplot.rcParams['figure.figsize'][0])
    fig, ax = pyplot.subplots(figsize=(x, x))
    fig.tight_layout(pad=0)
    im = ax.imshow(A, interpolation="none", vmin=0, vmax=1, cmap="Blues")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in 'top', 'bottom', 'left', 'right':
        ax.spines[s].set_visible(False)
    pyplot.gcf().colorbar(im, shrink=.5, pad=0)
    pyplot.savefig("_static/ld.svg")


def find_ld_sites(
        tree_sequence, focal_mutations, max_distance=1e6, r2_threshold=0.5,
        num_threads=8):
    """
    Finds all mutations within a given distance that are in approximate LD
    with a given set of mutations in a TreeSequence.
    """
    results = {}
    progress_bar = tqdm.tqdm(total=len(focal_mutations), ncols=90)
    num_threads = min(num_threads, len(focal_mutations))

    def thread_worker(thread_index):
        ld_calc = msprime.LdCalculator(tree_sequence)
        chunk_size = int(math.ceil(len(focal_mutations) / num_threads))
        start = thread_index * chunk_size
        for focal_mutation in focal_mutations[start: start + chunk_size]:
            a = ld_calc.get_r2_array(
                focal_mutation, max_distance=max_distance,
                direction=msprime.REVERSE)
            rev_indexes = focal_mutation - np.nonzero(a >= r2_threshold)[0] - 1
            a = ld_calc.get_r2_array(
                focal_mutation, max_distance=max_distance,
                direction=msprime.FORWARD)
            fwd_indexes = focal_mutation + np.nonzero(a >= r2_threshold)[0] + 1
            indexes = np.concatenate((rev_indexes[::-1], fwd_indexes))
            results[focal_mutation] = indexes
            progress_bar.update()

    threads = [
        threading.Thread(target=thread_worker, args=(j,))
        for j in range(num_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    progress_bar.close()
    return results

def threads_example():
    ts = msprime.simulate(
        sample_size=1000, Ne=1e4, length=1e7, recombination_rate=2e-8,
        mutation_rate=2e-8)
    counts = np.zeros(ts.get_num_mutations())
    for t in ts.trees():
        for mutation in t.mutations():
            counts[mutation.index] = t.get_num_leaves(mutation.node)
    doubletons = np.nonzero(counts == 2)[0]
    results = find_ld_sites(ts, doubletons, num_threads=8)
    print(
        "Found LD sites for", len(results), "doubleton mutations out of",
        ts.get_num_mutations())

def set_mutations_example():
    tree_sequence = msprime.simulate(
        sample_size=10000, Ne=1e4, length=1e7, recombination_rate=2e-8,
        mutation_rate=2e-8)
    print("Simulated ", tree_sequence.get_num_mutations(), "mutations")
    common_mutations = []
    for tree in tree_sequence.trees():
        for mutation in tree.mutations():
            p = tree.get_num_leaves(mutation.node) / tree.get_sample_size()
            if p >= 0.5:
                common_mutations.append(mutation)
    tree_sequence.set_mutations(common_mutations)
    print("Reduced to ", tree_sequence.get_num_mutations(), "common mutations")

def variants_example():
    tree_sequence = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=10)
    print("Simulated ", tree_sequence.get_num_mutations(), "mutations")
    for variant in tree_sequence.variants():
        print(variant.index, variant.position, variant.genotypes, sep="\t")

def variant_matrix_example():
    print("\nCreating full variant matrix")
    tree_sequence = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=10)
    shape = tree_sequence.get_num_mutations(), tree_sequence.get_sample_size()
    A = np.empty(shape, dtype="u1")
    for variant in tree_sequence.variants():
        A[variant.index] = variant.genotypes
    print(A)

def historical_samples_example():
    samples = [
        msprime.Sample(population=0, time=0),
        msprime.Sample(0, 0),  # Or, we can use positional arguments.
        msprime.Sample(0, 1.0)
    ]
    tree_seq = msprime.simulate(samples=samples, random_seed=5)
    tree = next(tree_seq.trees())
    for u in range(tree_seq.get_num_nodes()):
        print(u, tree.get_parent(u), tree.get_time(u), sep="\t")


def wright_fisher(N, T, L=100, random_seed=None):
    """
    Simulate a Wright-Fisher population of N haploid individuals with L
    discrete loci for T generations. Based on Algorithm W from
    https://www.biorxiv.org/content/biorxiv/early/2018/01/16/248500.full.pdf
    """
    random.seed(random_seed)
    tables = msprime.TableCollection(L)
    P = np.arange(N, dtype=int)
    # Mark the initial generation as samples so that we remember these nodes.
    for j in range(N):
        tables.nodes.add_row(time=T, flags=msprime.NODE_IS_SAMPLE)
    t = T
    while t > 0:
        t -= 1
        Pp = P.copy()
        for j in range(N):
            u = tables.nodes.add_row(time=t, flags=0)
            Pp[j] = u
            a = random.randint(0, N - 1)
            b = random.randint(0, N - 1)
            x = random.randint(1, L - 1)
            tables.edges.add_row(0, x, P[a], u)
            tables.edges.add_row(x, L, P[b], u)
        P = Pp

    # Now do some table manipulations to ensure that the tree sequence
    # that we output has the form that msprime needs to finish the
    # simulation. Much of the complexity here is caused by the tables API
    # not allowing direct access to memory, which will change soon.

    # Mark the extant population as samples also
    flags = tables.nodes.flags
    flags[P] = msprime.NODE_IS_SAMPLE
    tables.nodes.set_columns(flags=flags, time=tables.nodes.time)
    tables.sort()
    # Simplify with respect to the current generation, but ensuring we keep the
    # ancient nodes from the initial population.
    tables.simplify()
    # Unmark the initial generation as samples
    flags = tables.nodes.flags
    time = tables.nodes.time
    flags[:] = 0
    flags[time == 0] = msprime.NODE_IS_SAMPLE
    # The final tables must also have at least one population which
    # the samples are assigned to
    tables.populations.add_row()
    tables.nodes.set_columns(
        flags=flags, time=time,
        population=np.zeros_like(tables.nodes.population))
    return tables.tree_sequence()


def simulate_from_example():

    num_loci = 2
    wf_ts = wright_fisher(10, 5, L=num_loci, random_seed=3)
    for tree in wf_ts.trees():
        tree.draw(path="_static/simulate_from_wf_{}.svg".format(tree.index))

    recomb_map = msprime.RecombinationMap.uniform_map(num_loci, 1, num_loci)
    coalesced_ts = msprime.simulate(
        from_ts=wf_ts, recombination_map=recomb_map, random_seed=5)

    for tree in coalesced_ts.trees():
        tree.draw(path="_static/simulate_from_coalesced_{}.svg".format(tree.index))

    final_ts = coalesced_ts.simplify()

    for tree in final_ts.trees():
        print("interval = ", tree.interval)
        print(tree.draw(format="unicode"))


def full_arg_example():
    ts = msprime.simulate(
        sample_size=5, recombination_rate=0.1, record_full_arg=True,
        random_seed=42)
    print(ts.tables.nodes)
    print()
    for tree in ts.trees():
        print("interval:", tree.interval)
        print(tree.draw(format="unicode"))

def hybrid_sim_example():
    ts = msprime.simulate(
        sample_size=6, Ne=1000, model="dtwf", random_seed=2,
        demographic_events=[
            msprime.SimulationModelChange(time=500, model="hudson")])
    print(ts.tables.nodes)




if __name__ == "__main__":
    # single_locus_example()
    # multi_locus_example()
    # mutations_example()
    # set_mutations_example()
    # variants_example()
    # variant_matrix_example()
    # historical_samples_example()
    # segregating_sites_example(10, 5, 100000)
    # migration_example()
    # out_of_africa()
    # variable_recomb_example()
    # ld_matrix_example()
    # threads_example()
    # simulate_from_example()
    # full_arg_example()
    hybrid_sim_example()
