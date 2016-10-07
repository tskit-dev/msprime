"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import time
import collections
import math
import glob
import subprocess
import sys
import itertools
import threading
import random


import numpy as np
import numpy.ma as ma
import matplotlib
import scipy.stats
import pandas as pd
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import tqdm

import _msprime
import msprime

def mutations():
    n = 10
    # num_reps = 1000
    num_reps = 1
    num_loci = 10001
    # recomb_rates = [(1000, 0.005), (2000, 0.01), (3000, 0), (10001, 0.05)]
    recomb_rates = [(10001, 0.05)]
    last_pos = 0
    mean_rate = 0
    for pos, rate in recomb_rates:
        d = (pos - last_pos - 1) / (num_loci - 1)
        mean_rate += d * rate
        # print("mean_rate + ", d, rate)
        # print("rate = ", rate, rate / (4 * 10**4))
        last_pos = pos
    assert last_pos == num_loci
    print("mean_rate = ", mean_rate)
    num_trees = 0
    for j in range(num_reps):
        simulator = msprime.TreeSimulator(n)
        simulator.set_num_loci(num_loci)
        simulator.set_scaled_recombination_rate(mean_rate)
        # simulator.set_random_seed(j)
        simulator.run()
        num_trees += simulator.get_num_breakpoints()
        ts = simulator.get_tree_sequence()
        for t in ts.trees():
            print(t.get_interval()[0])

    # Construct the scrm command line. Use the first value as the background
    # rate
    simulator.set_scaled_recombination_rate(recomb_rates[0][-1])

    cmd = simulator.get_ms_command_line(
        "/home/jk/work/wt/papers/msprime/simulators/scrm",
        num_replicates=num_reps)
    for j in range(len(recomb_rates) - 1):
        pos = recomb_rates[j][0]
        # We still scale the recombination rate by the full locus length,
        # not the subset that we are working over.
        length = num_loci - 1
        rate = recomb_rates[j + 1][1]
        cmd += ["-sr", str(pos), str(rate * length)]
    # print(cmd)
    print(" ".join(cmd))
    result = subprocess.check_output(cmd)
    scrm_num_trees = 0
    for line in result.splitlines():
        # print(line)
        if line.startswith(b"["):
            scrm_num_trees += 1
    print(num_trees / num_reps, scrm_num_trees / num_reps)

    # tree_sequence = msprime.simulate(10, 100, mean_rate, random_seed=1)

    # for record in tree_sequence.records():
    #     print(record)
    # for tree in tree_sequence.trees():
    #     print(tree.get_interval())

def plot_distance_maps(recomb_rates):
    # Plot the piecewise map of physical distance to recombination rate
    x = np.zeros(2 * len(recomb_rates))
    y = np.copy(x)
    last_phys_x = 0
    j = 0
    for phys_x, recomb_rate in recomb_rates:
        x[j] = last_phys_x
        y[j] = recomb_rate
        j += 1
        x[j] = phys_x
        y[j] = recomb_rate
        last_phys_x = phys_x
        j += 1
    pyplot.plot(x, y)
    pyplot.ylim(-0.01, 1.01)
    pyplot.savefig("phys_recomb_rate.png")

    pyplot.clf()

    x = np.zeros(1 + len(recomb_rates))
    y = np.copy(x)
    j = 1
    s = 0
    last_phys_x = 0
    for phys_x, recomb_rate in recomb_rates:
        s += (phys_x - last_phys_x) * recomb_rate
        y[j] = s
        x[j] = phys_x
        j += 1
        last_phys_x = phys_x
    pyplot.plot(x, y)
    # physical_dist = 21.6
    # genetic_dist = physical_to_genetic(physical_dist, recomb_rates)
    genetic_dist = 4
    physical_dist = genetic_to_physical(genetic_dist, recomb_rates)
    pyplot.axvline(x=physical_dist, color="green")
    pyplot.axhline(y=genetic_dist, color="green")
    pyplot.savefig("phys_genetic_distance.png")


def plot_1kg_map():
    infile = "tmp__NOBACKUP__/genetic_map_b36/genetic_map_chr1_b36.txt.gz"

    import pandas as pd
    df = pd.read_csv(infile, delim_whitespace=True, compression="gzip",
            names=["pos", "rate", "distance"], header=0)
    # print(df.pos)
    physical_length = df.pos.iloc[-1]
    num_crossovers = df.distance.iloc[-1] / 100
    Ne = 10**4
    rate = 4 * Ne * num_crossovers / physical_length
    print("Overall rate = {:.2E}".format(rate))

    scaled_rate = np.array(4 * Ne * (df.rate / 100) / 10**6)[:-1]
    print(scaled_rate)

    lengths = np.diff(df.pos)
    print(lengths)

    print(lengths * scaled_rate)


    # print("overall rate = ",
    # print(df["pos"])
    pyplot.plot(df.pos, df.rate)
    pyplot.savefig("1kg.png")


def simulations():
    n = 10
    m = 1000
    recomb_map = msprime.RecombinationMap(
        m, [0, 0.5, 0.6, 0.7, 1], [0.1, 10, 0, 0.1, 0])
    sim = msprime.TreeSimulator(n)
    sim.set_random_seed(1)
    sim.set_num_loci(m)
    sim.set_recombination_map(recomb_map)
    # sim.set_scaled_recombination_rate(
    #     recomb_map.get_total_recombination_rate())
    sim.run()
    ts = sim.get_tree_sequence()
    size = 0
    for l, records_in, records_out in ts.diffs():
        # print(l, records_in, records_out)
        size += l
    print("size", size, ts.get_sequence_length())
    for t in ts.trees():
        l, r = t.get_interval()
        # print(l, r)
    for l, ns in ts.newick_trees():
        print(l, ns)
    # ts.generate_mutations(2, 1)
    # for t in ts.trees():
    #     l, r = t.get_interval()
    #     print("tree:", recomb_map.genetic_to_physical(l / m),
    #             recomb_map.genetic_to_physical(l / m))
    #     for pos, node in t.mutations():
    #         print("\t", node, pos, recomb_map.genetic_to_physical(pos / m),
    #                 sep="\t")

def convert_hdf5():
    in_filename = "tmp__NOBACKUP__/mutations.hdf5"
    out_filename = "tmp__NOBACKUP__/mutations_double_coords.hdf5"
    import h5py
    infile = h5py.File(in_filename, "r")
    outfile = h5py.File(out_filename, "w")
    # print(root)
    # g = root["trees"]
    # fields = [
    #     ("left", uint32, 1), ("right", uint32, 1),
    #     ("node", uint32, 1), ("children", uint32, 2),
    #     ("time", float64, 1)]
    #         self.assertEqual(g[name].shape[0], ts.get_num_records())

def read_1kg_map():
    infile = "tmp__NOBACKUP__/genetic_map_b36/genetic_map_chr1_b36.txt.gz"
    # infile = "genetic_map_chr22_b36.txt"
    infile = "tmp__NOBACKUP__/genetic_map_GRCh37_chr2.txt"
    pattern = "tmp__NOBACKUP__/genetic_map_GRCh37_chr*.txt"
    # pattern = "tmp__NOBACKUP__/genetic_map_GRCh37_chrX_par1.txt"
    for infile in glob.glob(pattern):
        name = infile.split("_")[-1].split(".")[0]
        print(infile, name)
        recomb_map = msprime.RecombinationMap.read_hapmap(infile)
        positions = np.array(recomb_map.get_positions())
        rates = np.array(recomb_map.get_rates())

        # tree_seq = msprime.simulate(10, recombination_map=recomb_map)
        n = 10
        before = time.clock()
        ts = msprime.simulate(
            n, Ne=10**4, recombination_map=recomb_map)
        print("Simulation ran in ", time.clock() - before)
        # for t in ts.trees():
        #     breakpoints.append(t.get_interval()[0])
        # b = np.array(breakpoints)

        # N = 500
        # fig, ax1 = pyplot.subplots(figsize=(16, 6))
        # v, bin_edges, bin_number = scipy.stats.binned_statistic(
        #     positions, rates, bins=N)
        # x = bin_edges[:-1][np.logical_not(np.isnan(v))]
        # y = v[np.logical_not(np.isnan(v))]
        # ax1.plot(x, y, "-")

        # ax2 = ax1.twinx()
        # v, bin_edges = np.histogram(b, N)
        # ax2.plot(bin_edges[:-1], v, color="green")

        # fig.savefig("tmp__NOBACKUP__/hapmap_{}.png".format(name))


    #     print(t.get_interval())
    # print(ts.get_num_records())


def genetic_to_phys(genetic_x, num_loci, positions, rates):
    total_recomb_rate = 0
    size = len(positions)
    for j in range(1, size):
        phys_length = positions[j] - positions[j - 1]
        total_recomb_rate += phys_length * rates[j - 1]
    if total_recomb_rate == 0:
        ret = (genetic_x / num_loci) * phys_length
    else:
        x = (genetic_x / num_loci) * total_recomb_rate
        ret = 0
        if x > 0:
            s = 0
            k = 0
            while s < x:
                s += (positions[k + 1] - positions[k]) * rates[k]
                k += 1
            excess = (s - x) / rates[k - 1]
            ret = positions[k] - excess
    return ret

def genetic_to_phys_bulk(values, num_loci, positions, rates):
    total_recomb_rate = 0
    size = len(positions)
    n = len(values)
    for j in range(1, size):
        phys_length = positions[j] - positions[j - 1]
        total_recomb_rate += phys_length * rates[j - 1]
    ret = list(values)
    if total_recomb_rate == 0:
        for j in range(n):
            ret[j] = genetic_to_phys(
                values[j], num_loci, positions, rates)
    else:
        # Get rid of zero values
        j = 0
        while values[j] == 0:
            j += 1
        s = 0
        k = 0
        while j < n:
            if j > 0 and values[j - 1] > values[j]:
                raise Exception("Input list not sorted")
            x = (values[j] / num_loci) * total_recomb_rate
            while s < x:
                s += (positions[k + 1] - positions[k]) * rates[k]
                k += 1
            excess = (s - x) / rates[k - 1]
            ret[j] = positions[k] - excess
            j += 1
    return ret


def map_stuff():
    num_loci = 1000
    positions = [0, 50, 80, 100]
    rates =     [0.2, 0.1, 0.0, 0]

    values = [0, 10, 50, 100, 900, 1000]
    bulk = genetic_to_phys_bulk(values, num_loci, positions, rates)

    for x, y in zip(values, bulk):
        phys = genetic_to_phys(x, num_loci, positions, rates)
        print(x, "\t", phys, "\t", y)


def new_api():
    # ts = msprime.simulate(10)

    infile = "hapmap/genetic_map_GRCh37_chr22.txt"
    recomb_map = msprime.RecombinationMap.read_hapmap(infile)
    ts = msprime.simulate(
        100, Ne=10**4,
        recombination_map=recomb_map,
        mutation_rate=1e-8)
    ts.dump("tmp__NOBACKUP__/chr22.hdf5")


def replicate_example():
    theta = 5
    R = 1000
    replicates = msprime.simulate(
        sample_size=100, recombination_rate=2, mutation_rate=theta/4,
        num_replicates=R, random_seed=None)
    S = np.zeros(R)
    T = np.zeros(R)
    for j, tree_sequence in enumerate(replicates):
        S[j] = tree_sequence.get_num_mutations()
        T[j] = tree_sequence.get_num_trees()
    print("theta =", theta, "mean(S) = ", np.mean(S))
    print(np.mean(T))

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


def segregating_sites_example(n, theta, num_replicates):
    S = np.zeros(num_replicates)
    replicates = msprime.simulate(
        sample_size=n,
        mutation_rate=theta / 4,
        num_replicates=num_replicates)
    for j, tree_sequence in enumerate(replicates):
        S[j] = tree_sequence.get_num_mutations()
    # Now, calculate the analytical predictions
    S_mean_a = np.sum(1 / np.arange(1, n)) * theta
    S_var_a = (
        theta * np.sum(1 / np.arange(1, n)) +
        theta**2 * np.sum(1 / np.arange(1, n)**2))
    print("              mean              variance")
    print("Observed      {}\t\t{}".format(np.mean(S), np.var(S)))
    print("Analytical    {:.5f}\t\t:.5f}".format(S_mean_a, S_var_a))
        # columns=["left", "right", "node", "children", "time", "population"])


def variable_recomb_example():
    infile = "hapmap/genetic_map_GRCh37_chr22.txt"
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
    fig.savefig("hapmap_chr22.svg")


def pop_example():
    if False:
        t = 100
        ts = msprime.simulate(
            Ne=10**4,
            population_configurations=[
                msprime.PopulationConfiguration(sample_size=1000),
                msprime.PopulationConfiguration(sample_size=1000),
                msprime.PopulationConfiguration(sample_size=1000),
                msprime.PopulationConfiguration(sample_size=1000),
                msprime.PopulationConfiguration(sample_size=1000)],
            demographic_events=[
                msprime.MassMigration(time=t, source=1, destination=0),
                msprime.MassMigration(time=t, source=2, destination=0),
                msprime.MassMigration(time=t, source=3, destination=0),
                msprime.MassMigration(time=t, source=4, destination=0)],
            length=100 * 1e6,
            recombination_rate=2e-8,
            mutation_rate=2e-8,
            random_seed=1)
        ts.dump("populations.hdf5")
        print(
            ts.get_sample_size(), ts.get_num_trees(),
            ts.get_num_mutations())
    else:
        ts = msprime.load("populations.hdf5")
        before = time.clock()
        R = 1
        for i in range(R):
            for j in range(5):
                samples = ts.get_samples(population_id=j)
                pi = ts.get_pairwise_diversity(samples)
                # pi2 = ts.get_pairwise_diversity2(samples)
                # print(j, pi, pi2, pi == pi2)
                # print(j, pi2)
        duration = time.clock() - before
        print("duration = ", duration, " per call = ", duration / (5 * R))


def vcf_example():

    # n = 6 # 3 diploid samples from each pop
    # t = 100
    # ts = msprime.simulate(
    #     Ne=10**4,
    #     population_configurations=[
    #         msprime.PopulationConfiguration(sample_size=n),
    #         msprime.PopulationConfiguration(sample_size=n),
    #         msprime.PopulationConfiguration(sample_size=n),
    #         msprime.PopulationConfiguration(sample_size=n),
    #         msprime.PopulationConfiguration(sample_size=n)],
    #     demographic_events=[
    #         msprime.MassMigration(time=t, source=1, destination=0),
    #         msprime.MassMigration(time=t, source=2, destination=0),
    #         msprime.MassMigration(time=t, source=3, destination=0),
    #         msprime.MassMigration(time=t, source=4, destination=0)],
    #     length=1 * 1e6,
    #     recombination_rate=2e-8,
    #     mutation_rate=2e-8,
    #     random_seed=1)
    # with open("test.vcf", "w") as f:

    #     ts.write_vcf(f, ploidy=2)

    ts = msprime.load("tmp__NOBACKUP__/populations.hdf5")
    before = time.clock()
    num_genotypes = 0
    for variant in ts.variants():
        num_genotypes += len(variant.genotypes)
    print(num_genotypes, ts.get_sample_size() * ts.get_num_mutations())
    duration = time.clock() - before
    print("Done in ", duration, " gives ",
            num_genotypes * 1e-6 / duration, " MGenotypes decoded per second")
    print(num_genotypes)


    before = time.clock()
    with open("tmp__NOBACKUP__/tmp_1.vcf", "w") as f:
        ts.write_vcf(f, ploidy=1)
        size = f.tell()
    duration = time.clock() - before
    print("wrote vcf in ", duration, "seconds", (size / 2**20) / duration, "MB/s")
    before = time.clock()
    with open("tmp__NOBACKUP__/tmp_2.vcf", "w") as f:
        ts.write_vcf(f, ploidy=2)
    duration = time.clock() - before
    print("wrote vcf in ", duration, "seconds", (size / 2**20) / duration, "MB/s")

def records_example():
    # filename = "records.txt"

    # ts = msprime.load("out.hdf5")
    # with open(filename, "w") as f:
    #     ts.write_records(f)

    # with open(filename, "r") as f:
    #     ts2 = msprime.TreeSequence.load_records(f)
    # for r1, r2 in zip(ts.records(), ts2.records()):
    #     print(r1.left, r2.left)

    ts = msprime.load_txt("example.txt")
    for t in ts.trees():
        print(t)

def stuff():
    before = time.clock()
    # Run the actual simulations
    tree_sequence = msprime.simulate(
        sample_size=10**5,
        length=100 * 10**6,
        Ne=1e4,
        demographic_events=[
            msprime.Bottleneck(time=100, proportion=0.1),
            msprime.Bottleneck(time=200, proportion=0.1),
            msprime.Bottleneck(time=300, proportion=0.1),
            msprime.Bottleneck(time=500, proportion=0.1)],
        recombination_rate=1e-8,
        mutation_rate=1e-8,
        random_seed=1  # Arbitrary - make this reproducible.
    )
    duration = time.clock() - before
    print("Simulated 100k genomes in {0:.3f} seconds.".format(duration))

    tree_sequence.dump("tmp__NOBACKUP__/bottleneck-example-new3.hdf5")

def examine():
    ts = msprime.load("tmp__NOBACKUP__/bottleneck-example.hdf5")
    print("num_records = ", ts.get_num_records())
    non_binary_records = 0
    max_record_length = 0
    for r in ts.records():
        if len(r.children) > 2:
            non_binary_records +=1
            max_record_length = max(max_record_length, len(r.children))
    print("non_binary_records = ", non_binary_records)
    print("max_record_length = ", max_record_length)
    num_nodes = collections.Counter()
    num_trees = 0
    for t in ts.trees():
        num_nodes[len(list(t.nodes(t.get_root())))] += 1
        num_trees += 1
    print("num_trees = ", num_trees)
    for k, v in num_nodes.items():
        print(k, "->", v)

def convert_dev():
    filename = "v2.hdf5"
    ts = msprime.read_legacy_hdf5(filename)
    msprime.write_legacy_hdf5(ts, "v2p.hdf5")
    # with msprime.Hdf5FileReader(filename) as reader:
    #     ts = reader.get_tree_sequence()
    # print("ts = ")
    # ll_ts = _msprime.TreeSequence()
    # ll_ts.load_records(records)
    # ts = msprime.TreeSequence(ll_ts)
    # ts.set_mutations(mutations)
    # print(ts.get_sample_size())
    ts.dump("v3.hdf5")


def ld_worker(ld_calc, start, stop, num_mutations, index, lock, progress):
    t = threading.current_thread()
    with lock:
        print("Thread ", t.name, "active")
    buff = bytearray(num_mutations * 8)
    for j in range(start, stop):
        v = ld_calc.get_r2(
            dest=buff, source_index=j, max_mutations=num_mutations - 1,
            max_distance=1e100)
        # print(list(buff[:v * 8]))
        a = np.frombuffer(buff, "d", v)
        x = np.mean(a)
        if j % 1000 == 0:
            with lock:
                progress[index] = (j - start) / (stop - start)
                s = "\t".join("{:.2f}".format(p) for p in progress)
                print(t.name, ":", s, end="\r")
                sys.stdout.flush()


def ld_dev():
    # ts = msprime.simulate(100, recombination_rate=10, mutation_rate=5,
    #         random_seed=1)
    num_threads = 10
    ts = msprime.load(sys.argv[1])
    print("num trees = ", ts.get_num_trees())
    print("num mutations  = ", ts.get_num_mutations())
    # num_mutations = min(ts.get_num_mutations(), 100000)
    # num_mutations = ts.get_num_mutations()
    num_mutations = 1000
    ld_calcs = [
        _msprime.LdCalculator(ts._ll_tree_sequence) for _ in range(num_threads)]
    k = ts.get_num_trees() // num_threads
    start = 0
    next_block = k
    intervals = []
    for t in ts.trees():
        if t.get_index() >= next_block:
            mutations = list(t.mutations())
            if len(mutations) > 0:
                stop = mutations[-1].index
                intervals.append((start, stop))
                start = stop
                next_block += k

    threads = []
    lock = threading.Lock()
    progress = [0 for j in range(num_threads)]
    for j in range(num_threads):
        start, stop = intervals[j]
        t = threading.Thread(
            name="ld_worker_{}".format(j), target=ld_worker,
            args=(ld_calcs[j], start, stop, num_mutations, j, lock, progress))
        t.start()
        threads.append(t)

    print("Main thread joining")
    for t in threads:
        t.join()
    print("Main thread done")


def ld_triangle_plot():
    ts = msprime.simulate(100, recombination_rate=10, mutation_rate=20,
            random_seed=1)

    print("num_mutations = ", ts.get_num_mutations())
    ld_calc = msprime.LdCalculator(ts)
    A = ld_calc.get_r2_matrix()

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
    pyplot.savefig("ld.png")


def find_ld_sites(
        tree_sequence, focal_mutations, max_distance=1e6, r2_threshold=0.5,
        num_threads=8):
    """
    Finds all mutations within a given distance that are in approximate LD
    with a given set of mutations in a TreeSequence.
    """
    results = {}
    progress_bar = tqdm.tqdm(total=len(focal_mutations))
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
    # ts = msprime.load(sys.argv[1])

    ts = msprime.simulate(
        sample_size=1000, Ne=1e4, length=1e7, recombination_rate=2e-8,
        mutation_rate=2e-8)
    np.random.seed(1)
    num_focal_mutations = 100
    print("num_mutations = ", ts.get_num_mutations())
    focal_mutations = np.sort(np.random.randint(
        ts.get_num_mutations(), size=num_focal_mutations))
    results = find_ld_sites(ts, focal_mutations, num_threads=8)
    print("found LD sites for", len(results), "random mutations")

    # for k, v in results.items():
    #     print(k, "has ", len(v), "mutation in LD")
    #     ld_calc = msprime.LdCalculator(ts)
    #     for j in v:
    #         print("\t", k, j, ld_calc.get_r2(k, j))


def api_stuff():
    ts = msprime.simulate(10, mutation_rate=10)

    for variant in ts.variants(as_bytes=True):
        print(variant.genotypes, type(variant.genotypes),
                variant.genotypes.decode())

if __name__ == "__main__":
    # mutations()

    # plot_distance_maps(
    #     [(10, 0.1), (11, 1), (20, 0.1), (21, 1), (30, 0.1)]
    # )
    # plot_1kg_map()

    # read_1kg_map()

    # simulations()
    # convert_hdf5()
    # map_stuff()
    # new_api()
    # replicate_example()
    # migration_example()
    # segregating_sites_example(2, 5, 10000)
    # variable_recomb_example()
    # pop_example()
    # vcf_example()
    # records_example()
    # stuff()
    # examine()
    # convert_dev()
    # ld_dev()
    # ld_triangle_plot()
    # threads_example()
    api_stuff()
