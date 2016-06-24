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

import scipy.special
import pandas as pd
import numpy as np
import numpy.random
import statsmodels.api as sm
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot

import dendropy
import msprime.cli as cli

import msprime


class SimulationVerifier(object):
    """
    Class to compare msprime against ms to ensure that the same distributions
    of values are output under the same parameters.
    """
    def __init__(self, output_dir):
        self._output_dir = output_dir
        self._instances = {}
        self._ms_executable = ["./data/ms/ms"]
        self._mspms_executable = ["python", "mspms_dev.py"]

    def get_ms_seeds(self):
        max_seed = 2**16
        seeds = [random.randint(1, max_seed) for j in range(3)]
        return ["-seed"] + map(str, seeds)

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
        return self._run_sample_stats(
            self._ms_executable + args.split() + self.get_ms_seeds())

    def _run_msprime_mutation_stats(self, args):
        return self._run_sample_stats(
            self._mspms_executable + args.split() + self.get_ms_seeds())

    def _run_ms_coalescent_stats(self, args):
        executable = ["./data/ms/ms_summary_stats"]
        with tempfile.TemporaryFile() as f:
            argList = executable + args.split() + self.get_ms_seeds()
            print("\t", " ".join(argList))
            subprocess.call(argList, stdout=f)
            f.seek(0)
            df = pd.read_table(f)
        return df

    def _run_msprime_coalescent_stats(self, args):
        print("\t msprime:", args)
        runner = cli.get_mspms_runner(args.split())
        sim = runner.get_simulator()
        rng = msprime.RandomGenerator(random.randint(1, 2**32 - 1))
        sim.set_random_generator(rng)
        num_populations = sim.get_num_populations()
        replicates = runner.get_num_replicates()
        num_trees = [0 for j in range(replicates)]
        time = [0 for j in range(replicates)]
        ca_events = [0 for j in range(replicates)]
        re_events = [0 for j in range(replicates)]
        mig_events = [None for j in range(replicates)]
        for j in range(replicates):
            sim.reset()
            sim.run()
            num_trees[j] = sim.get_num_breakpoints() + 1
            time[j] = sim.get_time() / 4  # Convert to coalescent units
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
            runner = self._instances[key]
            runner()

    def add_ms_instance(self, key, command_line):
        """
        Adds a test instance with the specified ms command line.
        """
        def f():
            print(key, command_line)
            self._run_coalescent_stats(key, command_line)
            self._run_mutation_stats(key, command_line)
        self._instances[key] = f

    def get_pairwise_coalescence_time(self, cmd, R):
        # print("\t", " ".join(cmd))
        output = subprocess.check_output(cmd)
        T = np.zeros(R)
        j = 0
        for line in output.splitlines():
            if line.startswith("("):
                t = dendropy.Tree.get_from_string(line, schema="newick")
                a = t.calc_node_ages()
                T[j] = a[-1]
                j += 1
        return T

    def run_pairwise_island_model(self):
        """
        Runs the check for the pairwise coalscence times for within
        and between populations.
        """
        R = 10000
        M = 0.2
        basedir = "tmp__NOBACKUP__/analytical_pairwise_island"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        for d in range(2, 6):
            cmd = "2 {} -T -I {} 2 {} {}".format(R, d, "0 " * (d - 1), M)
            T_w_ms = self.get_pairwise_coalescence_time(
                self._ms_executable + cmd.split() + self.get_ms_seeds(), R)
            T_w_msp = self.get_pairwise_coalescence_time(
                self._mspms_executable + cmd.split() + self.get_ms_seeds(), R)

            cmd = "2 {} -T -I {} 1 1 {} {}".format(R, d, "0 " * (d - 2), M)
            T_b_ms = self.get_pairwise_coalescence_time(
                self._ms_executable + cmd.split() + self.get_ms_seeds(), R)
            T_b_msp = self.get_pairwise_coalescence_time(
                self._mspms_executable + cmd.split() + self.get_ms_seeds(), R)
            print(d, np.mean(T_w_ms), np.mean(T_w_msp), d / 2,
                    np.mean(T_b_ms), np.mean(T_b_msp), (d + (d - 1) / M) / 2,
                    sep="\t")

            sm.graphics.qqplot(T_w_ms)
            sm.qqplot_2samples(T_w_ms, T_w_msp, line="45")
            f = os.path.join(basedir, "within_{}.png".format(d))
            pyplot.savefig(f, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(T_b_ms)
            sm.qqplot_2samples(T_b_ms, T_b_msp, line="45")
            f = os.path.join(basedir, "between_{}.png".format(d))
            pyplot.savefig(f, dpi=72)
            pyplot.close('all')

    def get_segregating_sites_histogram(self, cmd):
        print("\t", " ".join(cmd))
        output = subprocess.check_output(cmd)
        max_s = 200
        hist = np.zeros(max_s)
        for line in output.splitlines():
            if line.startswith("segsites"):
                s = int(line.split()[1])
                if s <= max_s:
                    hist[s] += 1
        return hist / np.sum(hist)

    def get_S_distribution(self, k, n, theta):
        """
        Returns the probability of having k segregating sites in a sample of
        size n. Wakely pg 94.
        """
        s = 0.0
        for i in range(2, n + 1):
            t1 = (-1)**i
            t2 = scipy.special.binom(n - 1, i - 1)
            t3 = (i - 1) / (theta + i - 1)
            t4 = (theta / (theta + i - 1))**k
            s += t1 * t2 * t3 * t4
        return s

    def run_s_analytical_check(self):
        """
        Runs the check for the number of segregating sites against the
        analytical prediction.
        """
        R = 100000
        theta = 2
        max_s = 20
        basedir = "tmp__NOBACKUP__/analytical_s"
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        for n in range(2, 15):
            cmd = "{} {} -t {}".format(n, R, theta)
            S_ms = self.get_segregating_sites_histogram(
                self._ms_executable + cmd.split() + self.get_ms_seeds())
            S_msp = self.get_segregating_sites_histogram(
                self._mspms_executable + cmd.split() + self.get_ms_seeds())
            filename = os.path.join(basedir, "{}.png".format(n))

            fig, ax = pyplot.subplots()
            index = np.arange(10)
            S_analytical = [self.get_S_distribution(j, n, theta) for j in index]
            bar_width = 0.35
            rects1 = pyplot.bar(
                index, S_ms[index], bar_width, color='b', label="ms")
            rects2 = pyplot.bar(
                index + bar_width, S_msp[index], bar_width, color='r', label="msp")
            pyplot.plot(index + bar_width, S_analytical, "o", color='k')
            pyplot.legend()
            pyplot.xticks(index + bar_width, [str(j) for j in index])
            pyplot.tight_layout()
            pyplot.savefig(filename)

    def run_pi_analytical_check(self):
        """
        Runs the check for pi against analytical predictions.
        """
        R = 100000
        theta = 4.5
        basedir = "tmp__NOBACKUP__/analytical_pi"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        sample_size = np.arange(2, 15)
        mean = np.zeros_like(sample_size, dtype=float)
        var = np.zeros_like(sample_size, dtype=float)
        predicted_mean = np.zeros_like(sample_size, dtype=float)
        predicted_var = np.zeros_like(sample_size, dtype=float)

        for k, n in enumerate(sample_size):
            pi = np.zeros(R)
            replicates = msprime.simulate(
                sample_size=n,
                mutation_rate=theta/4,
                num_replicates=R)
            for j, ts in enumerate(replicates):
                pi[j] = ts.get_pairwise_diversity()
            # Predicted mean is is theta.
            predicted_mean[k] = theta
            # From Wakely, eqn (4.14), pg. 101
            predicted_var[k] = (
                (n + 1) * theta / (3 * (n - 1)) +
                2 * (n**2 + n + 3) * theta**2 / (9 * n * (n - 1)))
            mean[k] = np.mean(pi)
            var[k] = np.var(pi)
            print(
                n, theta, np.mean(pi), predicted_var[k], np.var(pi),
                sep="\t")

        filename = os.path.join(basedir, "mean.png")
        pyplot.plot(sample_size, predicted_mean, "-")
        pyplot.plot(sample_size, mean, "-")
        pyplot.savefig(filename)
        pyplot.close('all')

        filename = os.path.join(basedir, "var.png")
        pyplot.plot(sample_size, predicted_var, "-")
        pyplot.plot(sample_size, var, "-")
        pyplot.savefig(filename)
        pyplot.close('all')


    def get_tbl_distribution(self, n, R, executable):
        """
        Returns an array of the R total branch length values from
        the specified ms-like executable.
        """
        cmd = executable + "{} {} -T -p 10".format(n, R).split()
        cmd += self.get_ms_seeds()
        print("\t", " ".join(cmd))
        output = subprocess.check_output(cmd)
        tbl = np.zeros(R)
        j = 0
        for line in output.splitlines():
            if line.startswith("("):
                t = dendropy.Tree.get_from_string(line, schema="newick")
                tbl[j] = t.length()
                j += 1
        return tbl

    def get_analytical_tbl(self, n, t):
        """
        Returns the probabily density of the total branch length t with
        a sample of n lineages. Wakeley Page 78.
        """
        t1 = (n - 1) / 2
        t2 = math.exp(-t / 2)
        t3 = pow(1 - math.exp(-t / 2), n - 2)
        return t1 * t2 * t3

    def run_tbl_analytical_check(self):
        """
        Runs the check for the total branch length.
        """
        R = 10000
        basedir = "tmp__NOBACKUP__/analytical_tbl"
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        for n in range(2, 15):
            tbl_ms = self.get_tbl_distribution(n, R, self._ms_executable)
            tbl_msp = self.get_tbl_distribution(n, R, self._mspms_executable)

            sm.graphics.qqplot(tbl_ms)
            sm.qqplot_2samples(tbl_ms, tbl_msp, line="45")
            filename = os.path.join(basedir, "qqplot_{}.png".format(n))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            hist_ms, bin_edges = np.histogram(tbl_ms, 20, density=True)
            hist_msp, _ = np.histogram(tbl_msp, bin_edges, density=True)

            index = bin_edges[:-1]
            # We don't seem to have the analytical value quite right here,
            # but since the value is so very close to ms's, there doesn't
            # seem to be much point in trying to fix it.
            analytical = [self.get_analytical_tbl(n, x * 2) for x in index]
            fig, ax = pyplot.subplots()
            bar_width = 0.15
            rects1 = pyplot.bar(
                index, hist_ms, bar_width, color='b', label="ms")
            rects2 = pyplot.bar(
                index + bar_width, hist_msp, bar_width, color='r', label="msp")
            pyplot.plot(index + bar_width, analytical, "o", color='k')
            pyplot.legend()
            # pyplot.xticks(index + bar_width, [str(j) for j in index])
            pyplot.tight_layout()
            filename = os.path.join(basedir, "hist_{}.png".format(n))
            pyplot.savefig(filename)

    def add_s_analytical_check(self):
        """
        Adds a check for the analytical predictions about the distribution
        of S, the number of segregating sites.
        """
        self._instances["analytical_s"] = self.run_s_analytical_check


    def add_pi_analytical_check(self):
        """
        Adds a check for the analytical predictions about the pi,
        the pairwise site diversity.
        """
        self._instances["analytical_pi"] = self.run_pi_analytical_check

    def add_total_branch_length_analytical_check(self):
        """
        Adds a check for the analytical check for the total branch length.
        """
        self._instances["analytical_tbl"] = self.run_tbl_analytical_check

    def add_pairwise_island_model_analytical_check(self):
        """
        Adds a check for the analytical check for pairwise island model
        """
        self._instances[
            "analytical_pairwise_island"] = self.run_pairwise_island_model

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

        self.add_ms_instance(key, cmd)


def main():
    # random.seed(2)
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
        "msdoc-two-species", (
        "15 10000 -t 11.2 -I 2 3 12 -g 1 44.36 -n 2 0.125 -eg 0.03125 1 0.0 "
        "-en 0.0625 2 0.05 -ej 0.09375 2 1"))
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
    verifier.add_random_instance(
        "random2", num_replicates=10**4, num_demographic_events=10)
    # verifier.add_random_instance("random2", num_populations=3)

    # Add analytical checks
    verifier.add_s_analytical_check()
    verifier.add_pi_analytical_check()
    verifier.add_total_branch_length_analytical_check()
    verifier.add_pairwise_island_model_analytical_check()

    keys = None
    if len(sys.argv) > 1:
        keys = sys.argv[1:]

    verifier.run(keys)

if __name__ == "__main__":
    main()
