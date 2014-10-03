"""
Script to automate verification of the msprime simulator against
Hudson's ms.
"""
from __future__ import print_function
from __future__ import division

import os
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

class MsSimulator(object):
    """
    Class representing Hudson's ms simulator. Takes care of running the
    simulations and collating results.
    """
    def __init__(self, n, m, Ne, r):
        self.sample_size = n
        self.recombination_rate = r
        self.num_loci = m
        self.effective_population_size = Ne
        self.executable = "./data/ms/ms_summary_stats"
        self.population_models = []
        if not os.path.exists(self.executable):
            raise ValueError("ms executable does not exist. "
                    + "go to data/ms directory and type make")

    def run(self, replicates):
        rho = get_scaled_recombination_rate(self.effective_population_size,
                self.num_loci, self.recombination_rate)
        args = [self.executable, str(self.sample_size), str(replicates), "-T",
                "-r", str(rho), str(self.num_loci)]
        # md = {POP_MODEL_CONSTANT: "-eN", POP_MODEL_EXP: "-eG"}
        for time, model, param in self.population_models:
            mstr = md[model]
            args.extend([mstr, str(time), str(param)])
        with tempfile.TemporaryFile() as f:
            # print(args)
            subprocess.call(args, stdout=f)
            f.seek(0)
            df = pd.read_table(f)
        return df

class MsprimeSimulator(object):

    def __init__(self, n, m, Ne, r):
        self.sample_size = n
        self.recombination_rate = r
        self.num_loci = m
        self.effective_population_size = Ne

    def run(self, replicates):
        num_trees = [0 for j in range(replicates)]
        time = [0 for j in range(replicates)]
        ca_events = [0 for j in range(replicates)]
        re_events = [0 for j in range(replicates)]
        for j in range(replicates):
            with tempfile.NamedTemporaryFile() as f:
                sim = msprime.TreeSimulator(self.sample_size, f.name)
                sim.set_recombination_rate(4 * self.effective_population_size
                        * self.recombination_rate)
                sim.set_num_loci(self.num_loci)
                sim.run()
                tf = msprime.TreeFile(f.name)
                num_trees[j] = tf.get_num_trees()
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
    ms = MsSimulator(n, m, r, Ne)
    df_ms = ms.run(num_replicates)
    msp = MsprimeSimulator(n, m, r, Ne)
    df_msp = msp.run(num_replicates)
    for stat in ["t", "num_trees", "re_events", "ca_events"]:
        v1 = df_ms[stat]
        v2 = df_msp[stat]
        sm.graphics.qqplot(v1)
        sm.qqplot_2samples(v1, v2, line="45")
        f = "{0}_{1}.png".format(output_prefix, stat)
        pyplot.savefig(f, dpi=72)
        pyplot.clf()


def main():
    # default to humanish recombination rates and population sizes.
    n = 400
    m = 100000
    Ne = 1e4
    r = 1e-8
    num_replicates = 1000
    models = []
    output_prefix = "verification__NOBACKUP__/simple"
    run_verify(n, m, Ne, r, models, num_replicates, output_prefix)

if __name__ == "__main__":
    main()
