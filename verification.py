"""
Script to automate verification of the msprime simulator against
Hudson's ms.
"""
from __future__ import print_function
from __future__ import division

import os
import random
import subprocess
import sys
import tempfile

import pandas as pd
import statsmodels.api as sm
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot

import msprime.cli as cli


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
            mig_events[j] = [
                r for row in sim.get_num_migration_events() for r in row]
        d = {
            "t": time, "num_trees": num_trees,
            "ca_events": ca_events, "re_events": re_events}
        for j in range(num_populations**2):
            events = [mig_events[k][j] for k in range(replicates)]
            d["mig_events_{}".format(j)] = events
        df = pd.DataFrame(d)
        return df

    def _build_filename(self, *args):
        output_dir = os.path.join(self._output_dir, args[0])
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        return os.path.join(output_dir, "_".join(args[1:]))

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

    def add_random_instance(
            self, key, num_populations=1, num_replicates=1000,
            num_demographic_events=0):
        m = random.randint(1, 1000)
        r = random.uniform(0.01, 0.1) * m
        theta = random.uniform(1, 100)
        N = num_populations
        sample_sizes = [random.randint(2, 10) for _ in range(N)]
        migration_matrix = [
            random.random() * (j % (N + 1) != 0) for j in range(N**2)]
        structure = ""
        if num_populations > 1:
            structure = "-I {} {} -ma {}".format(
                num_populations, " ".join(str(s) for s in sample_sizes),
                " ".join(str(r) for r in migration_matrix))
        cmd = "{} {} -t {} -r {} {} {}".format(
            sum(sample_sizes), num_replicates, theta, r, m, structure)

        if N > 1:
            # Add some migration matrix changes
            t = 0
            for j in range(1, 6):
                t += 0.125
                u = random.random()
                if u < 0.33:
                    cmd += " -eM {} {}".format(t, random.random())
                elif u < 0.66:
                    j = random.randint(1, N)
                    k = j
                    while k == j:
                        k = random.randint(1, N)
                    r = random.random()
                    cmd += " -em {} {}".format(t, j, k, r)
                else:
                    migration_matrix = [
                        random.random() * (j % (N + 1) != 0)
                        for j in range(N**2)]
                    cmd += " -ema {} {} {}".format(
                        t, N, " ".join(str(r) for r in migration_matrix))

        # Set some initial growth rates, etc.
        if N == 1:
            if random.random() < 0.5:
                cmd += " -G {}".format(random.random())
            else:
                cmd += " -eN 0 {}".format(random.random())
        # Add some demographic events
        t = 0
        for j in range(num_demographic_events):
            t += 0.125
            if random.random() < 0.5:
                cmd += " -eG {} {}".format(t, random.random())
            else:
                cmd += " -eN {} {}".format(t, random.random())

        self._instances[key] = cmd


def main():
    random.seed(2)
    verifier = SimulationVerifier("tmp__NOBACKUP__")

    # Try various options independently
    verifier.add_ms_instance(
        "size-change1", "10 10000 -t 2.0 -eN 0.1 2.0")
    verifier.add_ms_instance(
        "growth-rate-change1", "10 10000 -t 2.0 -eG 0.1 5.0")
    verifier.add_ms_instance(
        "growth-rate-2-pops1", "10 10000 -t 2.0 -I 2 5 5 2.5 -G 5.0")
    verifier.add_ms_instance(
        "growth-rate-2-pops2", "10 10000 -t 2.0 -I 2 5 5 2.5 -G 5.0 -g 1 0.1")
    verifier.add_ms_instance(
        "growth-rate-2-pops3", "10 10000 -t 2.0 -I 2 5 5 2.5 -g 1 0.1")
    verifier.add_ms_instance(
        "growth-rate-2-pops4", "10 10000 -t 2.0 -I 2 5 5 2.5 -eg 1.0 1 5.0")
    verifier.add_ms_instance(
        "pop-size-2-pops1", "100 10000 -t 2.0 -I 2 50 50 2.5 -n 1 0.1")
    verifier.add_ms_instance(
        "pop-size-2-pops2", "100 10000 -t 2.0 -I 2 50 50 2.5 -g 1 2 -n 1 0.1")
    verifier.add_ms_instance(
        "pop-size-2-pops3", "100 10000 -t 2.0 -I 2 50 50 2.5 -eN 0.5 3.5")
    verifier.add_ms_instance(
        "pop-size-2-pops4", "100 10000 -t 2.0 -I 2 50 50 2.5 -en 0.5 1 3.5")
    verifier.add_ms_instance(
        "migration-rate-2-pops1", "100 10000 -t 2.0 -I 2 50 50 0 -eM 3 5")
    verifier.add_ms_instance(
        "migration-matrix-2-pops1",
        "100 10000 -t 2.0 -I 2 50 50 -ma x 10 0 x")
    verifier.add_ms_instance(
        "migration-matrix-2-pops2",
        "100 10000 -t 2.0 -I 2 50 50 -m 1 2 10 -m 2 1 50")
    verifier.add_ms_instance(
        "migration-rate-change-2-pops1",
        "100 10000 -t 2.0 -I 2 50 50 -eM 5 10")
    verifier.add_ms_instance(
        "migration-matrix-entry-change-2-pops1",
        "100 10000 -t 2.0 -I 2 50 50 -em 0.5 2 1 10")
    verifier.add_ms_instance(
        "migration-matrix-change-2-pops1",
        "100 10000 -t 2.0 -I 2 50 50 -ema 10.0 2 x 10 0 x")
    verifier.add_ms_instance(
        "migration-matrix-change-2-pops2",
        "100 10000 -t 2.0 -I 2 50 50 -ema 1.0 2 x 0.1 0 x "
        "-eN 1.1 0 -ema 10 2 x 0 10 x")
    verifier.add_ms_instance(
        "population-split-2-pops1",
        "100 10000 -t 2.0 -I 2 50 50 5.0 -ej 2.0 1 2")
    verifier.add_ms_instance(
        "population-split-4-pops1",
        "100 10000 -t 2.0 -I 4 50 50 0 0 2.0 -ej 0.5 2 1")
    verifier.add_ms_instance(
        "population-split-4-pops2",
        "100 10000 -t 2.0 -I 4 25 25 25 25 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1")
    verifier.add_ms_instance(
        "population-split-4-pops3", (
        "100 10000 -t 2.0 -I 4 25 25 25 25 -ej 1 2 1 -em 1.5 4 1 2 "
        "-ej 2 3 1 -ej 3 4 1"))
    verifier.add_ms_instance(
        "admixture-1-pop1", "1000 1000 -t 2.0 -es 0.1 1 0.5 -em 0.1 1 2 1")
    verifier.add_ms_instance(
        "admixture-1-pop2", "1000 1000 -t 2.0 -es 0.1 1 0.1 -em 0.1 1 2 1")
    verifier.add_ms_instance(
        "admixture-1-pop3", "1000 1000 -t 2.0 -es 0.01 1 0.1 -em 0.1 2 1 1")
    verifier.add_ms_instance(
        "admixture-1-pop4",
        "1000 1000 -t 2.0 -es 0.01 1 0.1 -es 0.1 2 0 -em 0.1 3 1 1")
    verifier.add_ms_instance(
        "admixture-1-pop5",
        "1000 1000 -t 2.0 -es 0.01 1 0.1 -ej 1 2 1")
    verifier.add_ms_instance(
        "admixture-1-pop6", "1000 1000 -t 2.0 -es 0.01 1 0.0 -eg 0.02 2 5.0 ")
    verifier.add_ms_instance(
        "admixture-1-pop7", "1000 1000 -t 2.0 -es 0.01 1 0.0 -en 0.02 2 5.0 ")
    verifier.add_ms_instance(
        "admixture-2-pop1",
        "1000 1000 -t 2.0 -I 2 500 500 1 -es 0.01 1 0.1 -ej 1 3 1")
    verifier.add_ms_instance(
        "admixture-2-pop2",
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -em 2.0 3 1 1")
    verifier.add_ms_instance(
        "admixture-2-pop3", (
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -G 5.0 "
        "-em 2.0 3 1 1"))
    verifier.add_ms_instance(
        "admixture-2-pop4", (
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -eg 0.02 1 5.0 "
        "-em 0.02 3 1 1"))

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
    verifier.add_ms_instance(
        "msdoc-outgroup-sequence", "11 1000 -t 2.0 -I 2 1 10 -ej 6.0 1 2")
    verifier.add_ms_instance(
        "msdoc-stepping-stone", (
        "15 10000 -t 3.0 -I 6 0 7 0 0 8 0 -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 "
        "-m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 -em 2.0 3 4 "
        "2.5 -em 2.0 4 3 2.5"))

    # The order of simultaneous events matters in ms.
    verifier.add_ms_instance(
        "simultaneous-ex1", "10 10000 -t 2.0 -eN 0.3 0.5 -eG .3 7.0")
    # Add a bunch more instances...
    verifier.add_ms_instance(
        "zero-growth-rate", "10 10000 -t 2.0 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5")
    # Some examples provided by Konrad Lohse
    verifier.add_ms_instance(
        "konrad-1", (
        "4 1000 -t 2508 -I 2 2 2 0 -n 2 2.59 -ma x 0 1.502 x -ej 0.9485 1 2 "
        "-r 23.76 3000"))
    verifier.add_ms_instance(
        "konrad-2", (
        "3 10000 -t 0.423 -I 3 1 1 1 -es 0.0786 1 0.946635 -ej 0.0786 4 3 "
        "-ej 0.189256 1 2 -ej 0.483492 2 3"))
    verifier.add_ms_instance(
        "konrad-3", (
        "100 100 -t 2 -I 10 10 10 10 10 10 10 10 10 10 10 0.001 "))

    # Add some random instances.
    verifier.add_random_instance("random1")
    verifier.add_random_instance("random2", num_demographic_events=10)
    # verifier.add_random_instance("random2", num_populations=3)

    keys = None
    if len(sys.argv) > 1:
        keys = sys.argv[1:]

    verifier.run(keys)

if __name__ == "__main__":
    main()
