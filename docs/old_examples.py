"""
Examples that are used in the tutorial but haven't been ported to
using the new framework for test examples that get tested.
"""
import matplotlib.collections
import matplotlib.pyplot as pyplot
import numpy as np
import scipy.stats

import msprime

# Force matplotlib to not use any Xwindows backend.
matplotlib.use("Agg")  # noqa: E402


# TODO this example is a bit of a pain because we either need to put
# the hapmap example into the repo or we need to get it using stdpopsim.
# Getting it from stdpopsim is annoying because that'll mean we have to
# download every time for CI. This will slow CI down quite a bit.
# Do we get much from this example? If you want to simulate variable
# recomb from a map then you're a lot better off using stdpopsim.
def variable_recomb_example():
    infile = "../hapmap/genetic_map_GRCh37_chr22.txt"
    # Read in the recombination map using the read_hapmap method,
    recomb_map = msprime.RecombinationMap.read_hapmap(infile)

    # Now we get the positions and rates from the recombination
    # map and plot these using 500 bins.
    positions = np.array(recomb_map.get_positions()[1:])
    rates = np.array(recomb_map.get_rates()[1:])
    num_bins = 500
    v, bin_edges, _ = scipy.stats.binned_statistic(positions, rates, bins=num_bins)
    x = bin_edges[:-1][np.logical_not(np.isnan(v))]
    y = v[np.logical_not(np.isnan(v))]
    fig, ax1 = pyplot.subplots(figsize=(16, 6))
    ax1.plot(x, y, color="blue")
    ax1.set_ylabel("Recombination rate")
    ax1.set_xlabel("Chromosome position")

    # Now we run the simulation for this map. We assume Ne=10^4
    # and have a sample of 100 individuals
    tree_sequence = msprime.simulate(
        sample_size=100, Ne=10 ** 4, recombination_map=recomb_map
    )
    # Now plot the density of breakpoints along the chromosome
    breakpoints = np.array(list(tree_sequence.breakpoints()))
    ax2 = ax1.twinx()
    v, bin_edges = np.histogram(breakpoints, num_bins, density=True)
    ax2.plot(bin_edges[:-1], v, color="green")
    ax2.set_ylabel("Breakpoint density")
    ax2.set_xlim(1.5e7, 5.3e7)
    fig.savefig("_static/hapmap_chr22.svg")


def historical_samples_example():
    samples = [
        msprime.Sample(population=0, time=0),
        msprime.Sample(0, 0),  # Or, we can use positional arguments.
        msprime.Sample(0, 1.0),
    ]
    tree_seq = msprime.simulate(samples=samples, random_seed=5)
    tree = next(tree_seq.trees())
    for u in range(tree_seq.get_num_nodes()):
        print(u, tree.get_parent(u), tree.get_time(u), sep="\t")


# JK: commented this out for now to keep the linter happy.
# def simulate_from_example():

#     num_loci = 2
#     wf_ts = wright_fisher(10, 5, L=num_loci, random_seed=3)
#     for tree in wf_ts.trees():
#         tree.draw(path=f"_static/simulate_from_wf_{tree.index}.svg")

#     recomb_map = msprime.RecombinationMap.uniform_map(num_loci, 1, num_loci)
#     coalesced_ts = msprime.simulate(
#         from_ts=wf_ts, recombination_map=recomb_map, random_seed=5
#     )

#     for tree in coalesced_ts.trees():
#         tree.draw(path=f"_static/simulate_from_coalesced_{tree.index}.svg")

#     final_ts = coalesced_ts.simplify()

#     for tree in final_ts.trees():
#         print("interval = ", tree.interval)
#         print(tree.draw(format="unicode"))


def full_arg_example():
    ts = msprime.simulate(
        sample_size=5, recombination_rate=0.1, record_full_arg=True, random_seed=42
    )
    print(ts.tables.nodes)
    print()
    for tree in ts.trees():
        print("interval:", tree.interval)
        print(tree.draw(format="unicode"))


def hybrid_sim_example():
    ts = msprime.simulate(
        sample_size=6,
        Ne=1000,
        model="dtwf",
        random_seed=2,
        demographic_events=[msprime.SimulationModelChange(time=500, model="hudson")],
    )
    print(ts.tables.nodes)


if __name__ == "__main__":
    # historical_samples_example()
    # variable_recomb_example()
    # simulate_from_example()
    # full_arg_example()
    hybrid_sim_example()
