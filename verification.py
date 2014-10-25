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

import pandas as pd
import statsmodels.api as sm
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot


import msprime

def get_scaled_recombination_rate(Ne, m, r):
    """
    Returns rho = 4 * Ne * (m - 1) * r, the scaled recombination rate.
    """
    return 4 * Ne * (m - 1) * r

class Simulator(object):
    """
    Superclass of coalescent simulator objects.
    """
    def __init__(self, n, m, Ne, r, models=[]):
        self.sample_size = n
        self.recombination_rate = r
        self.num_loci = m
        self.effective_population_size = Ne
        self.population_models = models

class MsSimulator(Simulator):
    """
    Class representing Hudson's ms simulator. Takes care of running the
    simulations and collating results.
    """
    def run(self, replicates):
        executable = "./data/ms/ms_summary_stats"
        if not os.path.exists(executable):
            raise ValueError("ms executable does not exist. "
                    + "go to data/ms directory and type make")
        rho = get_scaled_recombination_rate(self.effective_population_size,
                self.num_loci, self.recombination_rate)
        args = [executable, str(self.sample_size), str(replicates), "-T",
                "-r", str(rho), str(self.num_loci)]
        for model in self.population_models:
            if isinstance(model, msprime.ConstantPopulationModel):
                v = ["-eN", str(model.start_time), str(model.size)]
            elif isinstance(model, msprime.ExponentialPopulationModel):
                v = ["-eG", str(model.start_time), str(model.alpha)]
            else:
                raise ValueError("unknown population model")
            args.extend(v)
        with tempfile.TemporaryFile() as f:
            print(args)
            subprocess.call(args, stdout=f)
            f.seek(0)
            df = pd.read_table(f)
        return df

class MsprimeSimulator(Simulator):
    """
    Class to simlify running the msprime simulator and getting summary
    stats over many replicates.
    """
    def run(self, replicates):
        num_trees = [0 for j in range(replicates)]
        time = [0 for j in range(replicates)]
        ca_events = [0 for j in range(replicates)]
        re_events = [0 for j in range(replicates)]
        with tempfile.NamedTemporaryFile() as f:
            for j in range(replicates):
                sim = msprime.TreeSimulator(self.sample_size, f.name)
                sim.set_scaled_recombination_rate(4
                    * self.effective_population_size * self.recombination_rate)
                sim.set_num_loci(self.num_loci)
                sim.set_max_memory("1G")
                for m in self.population_models:
                    sim.add_population_model(m)
                sim.run()
                num_trees[j] = sim.get_num_trees()
                tf = msprime.TreeFile(f.name)
                time[j] = sim.get_time()
                ca_events[j] = sim.get_num_coancestry_events()
                re_events[j] = sim.get_num_recombination_events()
        df = pd.DataFrame({"t":time, "num_trees":num_trees,
                "ca_events":ca_events, "re_events":re_events})
        return df


def run_verify(n, m, Ne, r, models, num_replicates, output_prefix):
    """
    Runs ms and msprime on the specified parameters and outputs qqplots
    with the specified prefix.
    """
    ms = MsSimulator(n, m, r, Ne, models)
    df_ms = ms.run(num_replicates)
    msp = MsprimeSimulator(n, m, r, Ne, models)
    df_msp = msp.run(num_replicates)
    for stat in ["t", "num_trees", "re_events", "ca_events"]:
        v1 = df_ms[stat]
        v2 = df_msp[stat]
        sm.graphics.qqplot(v1)
        sm.qqplot_2samples(v1, v2, line="45")
        f = "{0}_{1}.png".format(output_prefix, stat)
        pyplot.savefig(f, dpi=72)
        pyplot.clf()

def verify_random(k):

    random.seed(k)
    for j in range(k):
        n = random.randint(1, 100)
        m = random.randint(1, 10000)
        Ne = random.uniform(100, 1e4)
        r = random.uniform(1e-9, 1e-6)
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
        run_verify(n, m, Ne, r, models, num_replicates, output_prefix)
        break

def main():
    # default to humanish recombination rates and population sizes.
    n = 400
    m = 100000
    Ne = 1e4
    r = 1e-8
    num_replicates = 1000
    # num_replicates = 1
    models = [
            msprime.ConstantPopulationModel(0.1, 2.0),
            msprime.ConstantPopulationModel(0.4, 0.5),
            msprime.ExponentialPopulationModel(0.5, 1.0)]
    output_prefix = "tmp__NOBACKUP__/simple"
    # run_verify(n, m, Ne, r, models, num_replicates, output_prefix)
    # TODO definite problem here using random parameters.
    # - don't change until this has been  fixed.
    verify_random(100)

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
    run_verify(n, m, N0, r, models, num_replicates, output_prefix)

if __name__ == "__main__":
    # main()
    verify_human_demographics()
