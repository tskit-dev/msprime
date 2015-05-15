"""
Script to automate verification of the msprime simulator against
Hudson's ms.
"""
from __future__ import print_function
from __future__ import division

import os
import math
import random
import tempfile
import subprocess

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot


import msprime

def harmonic_number(n):
    """
    Returns the nth Harmonic number.
    """
    return sum(1 / k for k in range(1, n + 1))

def get_scaled_recombination_rate(Ne, m, r):
    """
    Returns rho = 4 * Ne * (m - 1) * r, the scaled recombination rate.
    """
    return 4 * Ne * (m - 1) * r

class Simulator(object):
    """
    Superclass of coalescent simulator objects.
    """
    def __init__(self, n, m, Ne, r, models=[], mutation_rate=None):
        self.sample_size = n
        self.recombination_rate = r
        self.num_loci = m
        self.effective_population_size = Ne
        self.population_models = models
        self.mutation_rate = mutation_rate

class MsSimulator(Simulator):
    """
    Class representing Hudson's ms simulator. Takes care of running the
    simulations and collating results.
    """
    def get_executable(self):
        return ["./data/ms/ms"]

    def generate_trees(self):
        return True

    def get_command_line(self, replicates):
        executable = self.get_executable()
        rho = get_scaled_recombination_rate(self.effective_population_size,
                self.num_loci, self.recombination_rate)
        args = executable + [str(self.sample_size), str(replicates)]
        if self.generate_trees():
            args += ["-T"]
        if self.num_loci > 1:
            args += ["-r", str(rho), str(self.num_loci)]
        if self.mutation_rate is not None:
            args += ["-t", str(self.mutation_rate)]
        for model in self.population_models:
            if isinstance(model, msprime.ConstantPopulationModel):
                v = ["-eN", str(model.start_time), str(model.size)]
            elif isinstance(model, msprime.ExponentialPopulationModel):
                v = ["-eG", str(model.start_time), str(model.alpha)]
            else:
                raise ValueError("unknown population model")
            args.extend(v)
        return args


class MsCoalescentStatisticsSimulator(MsSimulator):
    """
    A modified version of ms in which we output statistics about the
    coalescent algorithm.
    """
    def get_executable(self):
        return ["./data/ms/ms_summary_stats"]

    def run(self, replicates):
        with tempfile.TemporaryFile() as f:
            args = self.get_command_line(replicates)
            print(" ".join(args))
            subprocess.call(args, stdout=f)
            f.seek(0)
            df = pd.read_table(f)
        return df

class MutationStatisticsSimulator(object):
    """
    A mixin to run the simulation and pass the results through Hudson's
    sample_stats program.
    """
    def generate_trees(self):
        return False

    def run(self, replicates):
        args = self.get_command_line(replicates)
        print(" ".join(args))
        p1 = subprocess.Popen(args, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["./data/ms/sample_stats"], stdin=p1.stdout,
                stdout=subprocess.PIPE)
        p1.stdout.close()
        output = p2.communicate()[0]
        with tempfile.TemporaryFile() as f:
            f.write(output)
            f.seek(0)
            df = pd.read_table(f)
        return df



class MsMutationStatisticsSimulator(MutationStatisticsSimulator, MsSimulator):
    """
    Runs ms with a given set of parameters, and returns the mutation
    statistics.
    """

class MsprimeMutationStatisticsSimulator(MutationStatisticsSimulator, MsSimulator):
    """
    Runs msprime with a given set of parameters, and returns the mutation
    statistics.
    """
    def get_executable(self):
        return ["python", "mspms_dev.py"]

class MsprimeCoalescentStatisticsSimulator(Simulator):
    """
    Class to simlify running the msprime simulator and getting summary
    stats over many replicates.
    """
    def run(self, replicates):
        num_trees = [0 for j in range(replicates)]
        time = [0 for j in range(replicates)]
        ca_events = [0 for j in range(replicates)]
        re_events = [0 for j in range(replicates)]
        for j in range(replicates):
            sim = msprime.TreeSimulator(self.sample_size)
            sim.set_scaled_recombination_rate(4
                * self.effective_population_size * self.recombination_rate)
            sim.set_num_loci(self.num_loci)
            sim.set_max_memory("10G")
            for m in self.population_models:
                sim.add_population_model(m)
            tree_sequence = sim.run()
            num_trees[j] = sim.get_num_breakpoints()
            time[j] = sim.get_time()
            ca_events[j] = sim.get_num_coancestry_events()
            re_events[j] = sim.get_num_recombination_events()
        df = pd.DataFrame({"t":time, "num_trees":num_trees,
                "ca_events":ca_events, "re_events":re_events})
        return df

def run_verify_mutations(n, m, Ne, r, models, num_replicates, mutation_rate,
        output_prefix):
    """
    Runs ms and msprime for the specified parameters, and filters the results
    through Hudson's sample_stats program to get distributions of the
    haplotype statistics.
    """
    ms = MsMutationStatisticsSimulator(n, m, r, Ne, models, mutation_rate)
    df_ms = ms.run(num_replicates)
    msp = MsprimeMutationStatisticsSimulator(n, m, r, Ne, models, mutation_rate)
    df_msp = msp.run(num_replicates)
    for stat in ["pi", "ss", "D", "thetaH", "H"]:
        v1 = df_ms[stat]
        v2 = df_msp[stat]
        # pyplot.hist(v1, 20, alpha=0.5, label="ms")
        # pyplot.hist(v2, 20, alpha=0.5, label="msp")
        # pyplot.legend(loc="upper left")
        sm.graphics.qqplot(v1)
        sm.qqplot_2samples(v1, v2, line="45")
        f = "{0}_{1}.png".format(output_prefix, stat)
        pyplot.savefig(f, dpi=72)
        pyplot.clf()

def run_verify_coalescent(n, m, Ne, r, models, num_replicates, output_prefix):
    """
    Runs ms and msprime on the specified parameters and outputs qqplots
    of the coalescent simulation summary statistics with the specified
    prefix.
    """
    ms = MsCoalescentStatisticsSimulator(n, m, r, Ne, models)
    df_ms = ms.run(num_replicates)
    msp = MsprimeCoalescentStatisticsSimulator(n, m, r, Ne, models)
    df_msp = msp.run(num_replicates)
    for stat in ["t", "num_trees", "re_events", "ca_events"]:
        v1 = df_ms[stat]
        v2 = df_msp[stat]
        # pyplot.hist(v1, 20, alpha=0.5, label="ms")
        # pyplot.hist(v2, 20, alpha=0.5, label="msp")
        # pyplot.legend(loc="upper left")
        sm.graphics.qqplot(v1)
        sm.qqplot_2samples(v1, v2, line="45")
        f = "{0}_{1}.png".format(output_prefix, stat)
        pyplot.savefig(f, dpi=72)
        pyplot.clf()
        # pyplot.hist(v1, 20, alpha=0.5, label="ms")
        # pyplot.legend(loc="upper left")
        # f = "{0}_{1}_1.png".format(output_prefix, stat)
        # pyplot.savefig(f, dpi=72)
        # pyplot.clf()
        # pyplot.hist(v2, 20, alpha=0.5, label="msp")
        # pyplot.legend(loc="upper left")
        # f = "{0}_{1}_2.png".format(output_prefix, stat)
        # pyplot.savefig(f, dpi=72)
        # pyplot.clf()


def verify_random(k):

    random.seed(k)
    for j in range(k):
        n = random.randint(2, 100)
        m = random.randint(1, 10000)
        Ne = random.uniform(100, 1e4)
        r = random.uniform(1e-9, 1e-6)
        theta = random.uniform(1, 100)
        num_replicates = 1000
        output_prefix = "tmp__NOBACKUP__/random_{0}".format(j)
        models = []
        t = 0
        for j in range(random.randint(0, 10)):
            t += random.uniform(0, 0.3)
            p = random.uniform(0.1, 2.0)
            if random.random() < 0.5:
                mod = msprime.ConstantPopulationModel(t, p)
            else:
                mod = msprime.ExponentialPopulationModel(t, p)
            models.append(mod)
            print(mod.get_ll_model())
        print("running for", n, m, Ne, r, 4 * Ne * r)
        run_verify_coalescent(n, m, Ne, r, models, num_replicates, output_prefix)
        run_verify_mutations(n, m, Ne, r, models, num_replicates, theta, output_prefix)

def verify_exponential_models():
    random.seed(4)
    n = 15
    m = 4550
    Ne = 7730.75967602
    r = 7.05807713707e-07
    num_replicates = 10000
    output_prefix = "tmp__NOBACKUP__/expo_models"
    models = []
    t = 0.0
    for j in range(3):
        t += 0.1
        p = 100 * t
        mod = msprime.ExponentialPopulationModel(t, p)
        models.append(mod)
        print(mod.get_ll_model())
    # params = [(0.05, 0.1), (0.1, 0.2), (0.11, 1000), (0.15, 0.0001)]
    # models = [msprime.ConstantPopulationModel(t, p) for t, p in params]
    print("running for", n, m, Ne, r, 4 * Ne * r)
    run_verify_coalescent(n, m, Ne, r, models, num_replicates, output_prefix)

def verify_scrm_example():
    # -eN 0.3 0.5 -eG .3 7.0
    num_replicates = 10000
    models = [
            msprime.ConstantPopulationModel(0.3, 0.5),
            msprime.ExponentialPopulationModel(0.3, 7.0)]
    output_prefix = "tmp__NOBACKUP__/scrm"
    run_verify_coalescent(5, 1, 1, 0, models, num_replicates, output_prefix)

def verify_zero_growth_example():
    num_replicates = 10000
    models = [
            msprime.ExponentialPopulationModel(0.0, 6.93),
            msprime.ExponentialPopulationModel(0.2, 0.0),
            msprime.ConstantPopulationModel(0.3, 0.5)]
    output_prefix = "tmp__NOBACKUP__/zero"
    run_verify_coalescent(5, 1, 1, 0, models, num_replicates, output_prefix)

def verify_simple():
    # default to humanish recombination rates and population sizes.
    n = 400
    m = 100000
    Ne = 1e4
    r = 1e-8
    num_replicates = 1000
    models = [
            msprime.ConstantPopulationModel(0.1, 2.0),
            msprime.ConstantPopulationModel(0.4, 0.5),
            msprime.ExponentialPopulationModel(0.5, 1.0)]
    output_prefix = "tmp__NOBACKUP__/simple_coalescent"
    run_verify_coalescent(n, m, Ne, r, models, num_replicates, output_prefix)
    output_prefix = "tmp__NOBACKUP__/simple_mutations"
    run_verify_mutations(n, m, Ne, r, models, num_replicates, 10, output_prefix)

def verify_recombination_events():
    """
    Simple check to see if the expected number of recombination events
    is correct for large simulations.
    """
    n = 10000
    Ne = 10**4
    r = 1e-8
    num_replicates = 10
    for k in range(1, 21):
        m = k * 10**7
        msp = MsprimeSimulator(n, m, r, Ne, [])
        df = msp.run(num_replicates)
        R = get_scaled_recombination_rate(Ne, m, r)
        expected = R * harmonic_number(n - 1)
        print(m, df["num_trees"].mean(), expected, sep="\t")

def verify_mutations():

    n = 9
    m = 7165
    Ne = 3717
    r = 5.05e-07
    theta = 100
    num_replicates = 1000
    output_prefix = "tmp__NOBACKUP__/mutations"
    models = []
    run_verify_mutations(n, m, Ne, r, models, num_replicates, theta, output_prefix)
    run_verify_coalescent(n, m, Ne, r, models, num_replicates, output_prefix)


def main():
    # verify_recombination_events()
    verify_random(10)
    # verify_exponential_models()
    # verify_simple()
    # verify_zero_growth_example()
    # verify_scrm_example()
    verify_mutations()


def verify_human_demographics():
    """
    Model: 1e6 now, increasing from 2e4 400 generations ago
    (12kya), then 2e3 2000 generations ago (60kya) then 2e4 again fixed size
    beyond that.
    """
    n = 100
    m = 500000
    r = 1e-8
    num_replicates = 10000
    # Calculate the models
    N0 = 1e6
    N1 = 2e4
    N2 = 2e3
    N3 = 2e4
    g1 = 400
    g2 = 2000
    t1 = g1 / (4 * N0)
    t2 = g2 / (4 * N0)
    # Calculate the growth rates.
    alpha1 = -math.log(N1 / N0) / t1
    alpha2 = -math.log(N2 / N1) / (t2 - t1)
    # print(t1, t2)
    # print(alpha1, N0 * math.exp(- alpha1 * t1))
    # print(alpha2, N1 * math.exp(- alpha2 * (t2 - t1)))
    # num_replicates = 1
    models = [
            msprime.ExponentialPopulationModel(0, alpha1),
            msprime.ExponentialPopulationModel(t1, alpha2),
            msprime.ConstantPopulationModel(t2, N3 / N0),
            ]
    output_prefix = "tmp__NOBACKUP__/simple"
    run_verify_coalescent(n, m, N0, r, models, num_replicates, output_prefix)

if __name__ == "__main__":
    main()
    # verify_human_demographics()
