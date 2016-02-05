"""
Script to automate verification of the msprime simulator against
Hudson's ms.
"""
from __future__ import print_function
from __future__ import division

import math
import os
import random
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import statsmodels.api as sm
import allel
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot

import msprime
import msprime.cli as cli



# TODO remove these once they have been replaced.

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


#######################

class SimulationVerifier(object):
    """
    Class to compare msprime against ms to ensure that the same distributions
    of values are output under the same parameters.
    """
    def __init__(self, output_dir):
        self._output_dir = output_dir
        self._instances = {}

    def _run_sample_stats(self, args):
        print("\t", " ".join(args))
        p1 = subprocess.Popen(args, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["./data/ms/sample_stats"], stdin=p1.stdout,
            stdout=subprocess.PIPE)
        p1.stdout.close()
        output = p2.communicate()[0]
        with tempfile.TemporaryFile() as f:
            f.write(output)
            f.seek(0)
            df = pd.read_table(f)
        return df

    def _run_ms_mutation_stats(self, args):
        executable = "./data/ms/ms"
        return self._run_sample_stats([executable] + args.split())

    def _run_msprime_mutation_stats(self, args):
        executable = ["python", "mspms_dev.py"]
        return self._run_sample_stats(executable + args.split())

    def _run_ms_coalescent_stats(self, args):
        executable = ["./data/ms/ms_summary_stats"]
        with tempfile.TemporaryFile() as f:
            argList = executable + args.split()
            print("\t", " ".join(argList))
            subprocess.call(argList, stdout=f)
            f.seek(0)
            df = pd.read_table(f)
        return df

    def _run_msprime_coalescent_stats(self, args):
        print("\t msprime:", args)
        runner = cli.get_mspms_runner(args.split())
        sim = runner.get_simulator()
        num_populations = sim.get_num_populations()
        replicates = runner.get_num_replicates()
        num_trees = [0 for j in range(replicates)]
        time = [0 for j in range(replicates)]
        ca_events = [0 for j in range(replicates)]
        re_events = [0 for j in range(replicates)]
        mig_events = [None for j in range(replicates)]
        for j in range(replicates):
            sim.reset()
            sim.set_random_seed(j)
            sim.run()
            num_trees[j] = sim.get_num_breakpoints() + 1
            time[j] = sim.get_time()
            ca_events[j] = sim.get_num_common_ancestor_events()
            re_events[j] = sim.get_num_recombination_events()
            mig_events[j] = [r for row in sim.get_num_migration_events() for r in row]
        d = {
            "t":time, "num_trees":num_trees,
            "ca_events":ca_events, "re_events":re_events}
        for j in range(num_populations**2):
            events = [0 for j in range(replicates)]
            for k in range(replicates):
                events[k] = mig_events[k][j]
            d["mig_events_{}".format(j)] = events
        df = pd.DataFrame(d)
        return df

    def _build_filename(self, *args):
        return os.path.join(self._output_dir, "_".join(args))

    def _plot_stats(self, key, stats_type, df_msp, df_ms):
        assert set(df_ms.columns.values) == set(df_msp.columns.values)
        for stat in df_ms.columns.values:
            v1 = df_ms[stat]
            v2 = df_msp[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = self._build_filename(key, stats_type, stat)
            pyplot.savefig(f, dpi=72)
            pyplot.close('all')

    def _run_coalescent_stats(self, key, args):
        df_msp = self._run_msprime_coalescent_stats(args)
        df_ms = self._run_ms_coalescent_stats(args)
        self._plot_stats(key, "coalescent", df_ms, df_msp)

    def _run_mutation_stats(self, key, args):
        df_msp = self._run_msprime_mutation_stats(args)
        df_ms = self._run_ms_mutation_stats(args)
        self._plot_stats(key, "mutation", df_ms, df_msp)

    def run(self, keys=None):
        the_keys = sorted(self._instances.keys())
        if keys is not None:
            the_keys = keys

        for key in the_keys:
            args = self._instances[key]
            print(key, args)
            self._run_coalescent_stats(key, args)
            self._run_mutation_stats(key, args)

    def add_ms_instance(self, key, command_line):
        """
        Adds a test instance with the specified ms command line.
        """
        self._instances[key] = command_line


def main():
    verifier = SimulationVerifier("tmp__NOBACKUP__")

    # Try various options independently
    verifier.add_ms_instance(
        "size-change1", "10 10000 -t 2.0 -eN 0.1 2.0")
    verifier.add_ms_instance(
        "growth-rate-change1", "10 10000 -t 2.0 -eG 0.1 5.0")
    # Examples from ms documentation
    verifier.add_ms_instance(
        "msdoc-simple-ex", "4 20000 -t 5.0")
    verifier.add_ms_instance(
        "msdoc-recomb-ex", "15 1000 -t 10.04 -r 100.0 2501")
    verifier.add_ms_instance(
        "msdoc-structure-ex1", "15 1000 -t 2.0 -I 3 10 4 1 5.0")
    verifier.add_ms_instance(
        "msdoc-structure-ex2",
        "15 1000 -t 2.0 -I 3 10 4 1 5.0 -m 1 2 10.0 -m 2 1 9.0")
    verifier.add_ms_instance(
        "msdoc-structure-ex3",
        "15 1000 -t 10.0 -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x")
    # The order of simultaneous events matters in ms.
    verifier.add_ms_instance(
        "simultaneous-ex1", "10 10000 -t 2.0 -eN 0.3 0.5 -eG .3 7.0")
    # Add a bunch more instances...

    keys = None
    if len(sys.argv) > 1:
        keys = sys.argv[1:]

    verifier.run(keys)

if __name__ == "__main__":
    main()
