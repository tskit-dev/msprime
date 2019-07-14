"""
Script to automate verification of the msprime simulator against
known statistical results and benchmark programs such as Hudson's ms.
"""
import collections
import math
import os
import random
import subprocess
import sys
import tempfile
import time

import scipy.special
import pandas as pd
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
# Note this must be done before importing statsmodels.
matplotlib.use('Agg') # NOQA
from matplotlib import pyplot
from matplotlib import lines as mlines
import seaborn as sns
import statsmodels.api as sm
import dendropy
import tqdm
import argparse

import msprime.cli as cli
import msprime


def harmonic_number(n):
    return np.sum(1 / np.arange(1, n + 1))


def hk_f(n, z):
    """
    Returns Hudson and Kaplan's f_n(z) function. This is based on the exact
    value for n=2 and the approximations given in the 1985 Genetics paper.
    """
    ret = 0
    if n == 2:
        ret = (18 + z) / (z**2 + 13 * z + 18)
    else:
        ret = sum(1 / j**2 for j in range(1, n)) * hk_f(2, z)
    return ret


def get_predicted_variance(n, R):
    # We import this here as it's _very_ slow to import and we
    # only use it in this case.
    import scipy.integrate

    def g(z):
        return (R - z) * hk_f(n, z)
    res, err = scipy.integrate.quad(g, 0, R)
    return R * harmonic_number(n - 1) + 2 * res


def write_slim_script(outfile, format_dict):
    slim_str = """
    // set up a simple neutral simulation
    initialize()
    {{
        initializeTreeSeq(checkCoalescence=T);
        initializeMutationRate(0);
        initializeMutationType('m1', 0.5, 'f', 0.0);
        // g1 genomic element type: uses m1 for all mutations
        initializeGenomicElementType('g1', m1, 1.0);
        // uniform chromosome
        initializeGenomicElement(g1, 0, {NUM_LOCI});
        // uniform recombination along the chromosome
        initializeRecombinationRate({RHO});
    }}
    // create a population
    1
    {{
        {POP_STRS};
        sim.tag = 0;
    }}
    // run for set number of generations
    1: late()
    {{
        if (sim.tag == 0) {{
            if (sim.treeSeqCoalesced()) {{
                sim.tag = sim.generation;
                catn(sim.tag + ': COALESCED');
            }}
        }}
        if (sim.generation == sim.tag * 10) {{
            sim.simulationFinished();
            catn('Ran a further ' + sim.tag * 10 + ' generations');
            sim.treeSeqOutput('{OUTFILE}');
        }}
    }}
    100000 late() {{
        catn('No coalescence after 100000 generations!');
    }}
    """
    with open(outfile, 'w') as f:
        f.write(slim_str.format(**format_dict))


def subsample_simplify_slim_treesequence(ts, sample_sizes):
    tables = ts.dump_tables()
    samples = set(ts.samples())
    num_populations = len(set(tables.nodes.population))
    assert len(sample_sizes) == num_populations

    subsample = []
    for i, size in enumerate(sample_sizes):
        # Stride 2 to only sample one chrom per diploid SLiM individual
        ss = np.where(tables.nodes.population == i)[0][::2]
        ss = list(samples.intersection(ss))
        ss = np.random.choice(ss, replace=False, size=size)
        subsample.extend(ss)

    tables.nodes.individual = None
    tables.individuals.clear()
    tables.simplify(subsample)
    ts = tables.tree_sequence()

    return ts


class SimulationVerifier(object):
    """
    Class to compare msprime against ms to ensure that the same distributions
    of values are output under the same parameters.
    """
    def __init__(self, output_dir):
        self._output_dir = output_dir
        self._instances = {}
        self._ms_executable = ["./data/ms"]
        self._scrm_executable = ["./data/scrm"]
        self._slim_executable = ["./data/slim"]
        self._discoal_executable = ["./data/discoal"]
        self._mspms_executable = [sys.executable, "mspms_dev.py"]

    def check_slim_version(self):
        # This may not be robust but it's a start
        min_version = 3.1
        raw_str = subprocess.check_output(self._slim_executable + ["-version"])
        version_list = str.split(str(raw_str))
        for i in range(len(version_list)):
            if version_list[i].lower() == 'version':
                version_str = version_list[i+1]
                break
        version = float(version_str.strip(' ,')[0:3])
        assert version >= min_version, "Require SLiM >= 3.1!"

    def get_ms_seeds(self):
        max_seed = 2**16
        seeds = [random.randint(1, max_seed) for j in range(3)]
        return ["-seed"] + list(map(str, seeds))

    def get_discoal_seeds(self):
        max_seed = 2**16
        seeds = [random.randint(1, max_seed) for j in range(3)]
        return ["-d"] + list(map(str, seeds))

    def _run_sample_stats(self, args):
        print("\t", " ".join(args))
        p1 = subprocess.Popen(args, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["./data/sample_stats"], stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        output = p2.communicate()[0]
        p1.wait()
        if p1.returncode != 0:
            raise ValueError("Error occured in subprocess: ", p1.returncode)
        with tempfile.TemporaryFile() as f:
            f.write(output)
            f.seek(0)
            df = pd.read_table(f)
        return df

    def _run_discoal_mutation_stats(self, args):
        return self._run_sample_stats(
            self._discoal_executable + args.split() + self.get_discoal_seeds())

    def _run_ms_mutation_stats(self, args):
        return self._run_sample_stats(
            self._ms_executable + args.split() + self.get_ms_seeds())

    def _run_msprime_mutation_stats(self, args):
        return self._run_sample_stats(
            self._mspms_executable + args.split() + self.get_ms_seeds())

    def _run_ms_coalescent_stats(self, args):
        executable = ["./data/ms_summary_stats"]
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
        sim.random_generator = rng
        num_populations = sim.num_populations
        replicates = runner.get_num_replicates()
        num_trees = [0 for j in range(replicates)]
        time = [0 for j in range(replicates)]
        ca_events = [0 for j in range(replicates)]
        re_events = [0 for j in range(replicates)]
        mig_events = [None for j in range(replicates)]
        for j in range(replicates):
            sim.reset()
            sim.run()
            num_trees[j] = sim.num_breakpoints + 1
            time[j] = sim.time
            ca_events[j] = sim.num_common_ancestor_events
            re_events[j] = sim.num_recombination_events
            mig_events[j] = [r for row in sim.num_migration_events for r in row]
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
            print(key)
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

    def _discoal_str_to_ms(self, args):
        # convert discoal string to msprime string
        tokens = args.split(" ")
        # cut out sites param
        del(tokens[2])
        # adjust popIDs
        for i in range(len(tokens)):
            # pop size change case
            if tokens[i] == "-en":
                tokens[i+2] = str(int(tokens[i + 2]) + 1)
            # migration rate case
            if tokens[i] == "-m":
                tokens[i+1] = str(int(tokens[i + 1]) + 1)
                tokens[i+2] = str(int(tokens[i + 2]) + 1)
        msp_str = " ".join(tokens)
        return(msp_str)

    def _run_mutation_discoal_stats(self, key, args):
        msp_str = self._discoal_str_to_ms(args)
        df_msp = self._run_msprime_mutation_stats(msp_str)
        df_d = self._run_discoal_mutation_stats(args)
        self._plot_stats(key, "mutation", df_d, df_msp)

    def add_discoal_instance(self, key, command_line):
        """
        Adds a test instance with the specified discoal command line.
        """
        def f():
            print(key, command_line)
            self._run_mutation_discoal_stats(key, command_line)
        self._instances[key] = f

    def get_pairwise_coalescence_time(self, cmd, R):
        # print("\t", " ".join(cmd))
        output = subprocess.check_output(cmd)
        T = np.zeros(R)
        j = 0
        for line in output.splitlines():
            if line.startswith(b"("):
                t = dendropy.Tree.get_from_string(line.decode(), schema="newick")
                a = t.calc_node_ages()
                T[j] = a[-1]
                j += 1
        return T

    def run_arg_recording(self):
        basedir = "tmp__NOBACKUP__/arg_recording"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        ts_node_counts = np.array([])
        arg_node_counts = np.array([])
        ts_tree_counts = np.array([])
        arg_tree_counts = np.array([])
        ts_edge_counts = np.array([])
        arg_edge_counts = np.array([])

        reps = 10000
        leaves = 1000
        rho = 0.2

        for i in range(reps):
            ts = msprime.simulate(
                sample_size=leaves,
                recombination_rate=rho,
                random_seed=i+1)
            ts_node_counts = np.append(ts_node_counts, ts.num_nodes)
            ts_tree_counts = np.append(ts_tree_counts, ts.num_trees)
            ts_edge_counts = np.append(ts_edge_counts, ts.num_edges)
            arg = msprime.simulate(
                sample_size=leaves,
                recombination_rate=rho,
                random_seed=i + 1,
                record_full_arg=True)
            arg = arg.simplify()
            arg_node_counts = np.append(arg_node_counts, arg.num_nodes)
            arg_tree_counts = np.append(arg_tree_counts, arg.num_trees)
            arg_edge_counts = np.append(arg_edge_counts, arg.num_edges)

        pp_ts = sm.ProbPlot(ts_node_counts)
        pp_arg = sm.ProbPlot(arg_node_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "nodes.png")
        pyplot.savefig(f, dpi=72)

        pp_ts = sm.ProbPlot(ts_tree_counts)
        pp_arg = sm.ProbPlot(arg_tree_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "trees.png")
        pyplot.savefig(f, dpi=72)

        pp_ts = sm.ProbPlot(ts_edge_counts)
        pp_arg = sm.ProbPlot(arg_edge_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "edges.png")
        pyplot.savefig(f, dpi=72)
        pyplot.close('all')

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
            print(
                d, np.mean(T_w_ms), np.mean(T_w_msp), d / 2,
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
            if line.startswith(b"segsites"):
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
            pyplot.bar(
                index, S_ms[index], bar_width, color='b', label="ms")
            pyplot.bar(
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

    def run_mean_coaltime_check(self):
        """
        Checks the mean coalescence time calculation against pi.
        """
        random.seed(5)
        num_models = 8
        num_reps = 8
        T = np.zeros((num_models, num_reps))
        U = np.zeros(num_models)
        print("coaltime: theory  mean  sd   z")
        for k in range(num_models):
            Ne = 100
            N = 4
            pop_sizes = [random.uniform(0.01, 10) * Ne for _ in range(N)]
            growth_rates = [random.uniform(-0.01, 0.01) for _ in range(N)]
            migration_matrix = [
                [random.random() * (i != j) for j in range(N)]
                for i in range(N)]
            sample_size = [random.randint(2, 10) for _ in range(N)]
            population_configurations = [
                    msprime.PopulationConfiguration(
                        initial_size=k,
                        sample_size=n,
                        growth_rate=r)
                    for k, n, r in zip(pop_sizes, sample_size, growth_rates)]
            demographic_events = []
            for i in [0, 1]:
                n = random.uniform(0.01, 10)
                r = 0
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=100, initial_size=n, growth_rate=r, population_id=i))
            for ij in [(0, 1), (2, 3), (0, 3)]:
                demographic_events.append(
                        msprime.MigrationRateChange(
                            180, random.random(),
                            matrix_index=ij))
            demographic_events.append(
                msprime.MassMigration(time=200, source=3, dest=0, proportion=0.3))
            for i in [1, 3]:
                n = random.uniform(0.01, 10)
                r = random.uniform(-0.01, 0.01)
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=210, initial_size=n, growth_rate=r, population_id=i))

            ddb = msprime.DemographyDebugger(
                    population_configurations=population_configurations,
                    demographic_events=demographic_events,
                    migration_matrix=migration_matrix)

            U[k] = ddb.mean_coalescence_time(num_samples=sample_size)

            mut_rate = 1e-8
            replicates = msprime.simulate(
                    length=1e7,
                    recombination_rate=1e-8,
                    mutation_rate=mut_rate,
                    population_configurations=population_configurations,
                    demographic_events=demographic_events,
                    migration_matrix=migration_matrix,
                    random_seed=5, num_replicates=num_reps)
            for j, ts in enumerate(replicates):
                T[k, j] = ts.get_pairwise_diversity()
                T[k, j] /= ts.sequence_length
                T[k, j] /= 2 * mut_rate
            mT = np.mean(T[k])
            sT = np.std(T[k])
            print("        {:.2f} {:.2f} {:.2f} {:.2f}".format(
                    U[k], mT, sT, (U[k] - mT)/(sT * np.sqrt(num_reps))))

        basedir = "tmp__NOBACKUP__/coaltime"
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        fig, ax = pyplot.subplots()
        ax.scatter(np.column_stack([U]*T.shape[1]), T)
        # where oh where is abline(0,1)
        line = mlines.Line2D([0, 1], [0, 1])
        line.set_transform(ax.transAxes)
        ax.add_line(line)
        ax.set_xlabel("calculated mean coaltime")
        ax.set_ylabel("pairwise diversity, scaled")
        filename = os.path.join(basedir, "mean_coaltimes.png")
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
            if line.startswith(b"("):
                t = dendropy.Tree.get_from_string(line.decode(), schema="newick")
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
            pyplot.bar(
                index, hist_ms, bar_width, color='b', label="ms")
            pyplot.bar(
                index + bar_width, hist_msp, bar_width, color='r', label="msp")
            pyplot.plot(index + bar_width, analytical, "o", color='k')
            pyplot.legend()
            # pyplot.xticks(index + bar_width, [str(j) for j in index])
            pyplot.tight_layout()
            filename = os.path.join(basedir, "hist_{}.png".format(n))
            pyplot.savefig(filename)

    def get_num_trees(self, cmd, R):
        print("\t", " ".join(cmd))
        output = subprocess.check_output(cmd)
        T = np.zeros(R)
        j = -1
        for line in output.splitlines():
            if line.startswith(b"//"):
                j += 1
            if line.startswith(b"["):
                T[j] += 1
        return T

    def get_scrm_num_trees(self, cmd, R):
        print("\t", " ".join(cmd))
        output = subprocess.check_output(cmd)
        T = np.zeros(R)
        j = -1
        for line in output.splitlines():
            if line.startswith(b"//"):
                j += 1
            if line.startswith(b"time"):
                T[j] += 1
        return T

    def get_scrm_oldest_time(self, cmd, R):
        print("\t", " ".join(cmd))
        output = subprocess.check_output(cmd)
        T = np.zeros(R)
        j = -1
        for line in output.splitlines():
            if line.startswith(b"//"):
                j += 1
            if line.startswith(b"time:"):
                T[j] = max(T[j], float(line.split()[1]))
        return T

    def run_cli_num_trees(self):
        """
        Runs the check for number of trees using the CLI.
        """
        r = 1e-8  # Per generation recombination rate.
        num_loci = np.linspace(100, 10**5, 10).astype(int)
        Ne = 10**4
        n = 100
        rho = r * 4 * Ne * (num_loci - 1)
        num_replicates = 100
        ms_mean = np.zeros_like(rho)
        msp_mean = np.zeros_like(rho)
        for j in range(len(num_loci)):
            cmd = "{} {} -T -r {} {}".format(
                n, num_replicates, rho[j], num_loci[j])
            T = self.get_num_trees(
                self._ms_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates)
            ms_mean[j] = np.mean(T)

            T = self.get_num_trees(
                self._mspms_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates)
            msp_mean[j] = np.mean(T)
        basedir = "tmp__NOBACKUP__/cli_num_trees"
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        pyplot.plot(rho, ms_mean, "o")
        pyplot.plot(rho, msp_mean, "^")
        pyplot.plot(rho, rho * harmonic_number(n - 1), "-")
        filename = os.path.join(basedir, "mean.png")
        pyplot.savefig(filename)
        pyplot.close('all')

    def run_smc_oldest_time(self):
        """
        Runs the check for number of trees using the CLI.
        """
        r = 1e-8  # Per generation recombination rate.
        num_loci = np.linspace(100, 10**5, 10).astype(int)
        Ne = 10**4
        n = 100
        rho = r * 4 * Ne * (num_loci - 1)
        num_replicates = 1000
        scrm_mean = np.zeros_like(rho)
        scrm_smc_mean = np.zeros_like(rho)
        msp_mean = np.zeros_like(rho)
        msp_smc_mean = np.zeros_like(rho)
        for j in range(len(num_loci)):

            cmd = "{} {} -L -r {} {} -p 14".format(
                n, num_replicates, rho[j], num_loci[j])
            T = self.get_scrm_oldest_time(
                self._scrm_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates)
            scrm_mean[j] = np.mean(T)

            cmd += " -l 0"
            T = self.get_scrm_oldest_time(
                self._scrm_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates)
            scrm_smc_mean[j] = np.mean(T)

            for dest, model in [(msp_mean, "hudson"), (msp_smc_mean, "smc_prime")]:
                replicates = msprime.simulate(
                    sample_size=n, length=num_loci[j],
                    recombination_rate=r, Ne=Ne, num_replicates=num_replicates,
                    model=model)
                T = np.zeros(num_replicates)
                for k, ts in enumerate(replicates):
                    for record in ts.records():
                        T[k] = max(T[k], record.time)
                # Normalise back to coalescent time.
                T /= 4 * Ne
                dest[j] = np.mean(T)
        basedir = "tmp__NOBACKUP__/smc_oldest_time"
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        pyplot.plot(rho, scrm_mean, "-", color="blue", label="scrm")
        pyplot.plot(rho, scrm_smc_mean, "-", color="red", label="scrm_smc")
        pyplot.plot(rho, msp_smc_mean, "--", color="red", label="msprime_smc")
        pyplot.plot(rho, msp_mean, "--", color="blue", label="msprime")
        pyplot.xlabel("rho")
        pyplot.ylabel("Mean oldest coalescence time")
        pyplot.legend(loc="lower right")
        filename = os.path.join(basedir, "mean.png")
        pyplot.savefig(filename)
        pyplot.close('all')

    def run_smc_num_trees(self):
        """
        Runs the check for number of trees in the SMC and full coalescent
        using the API. We compare this with scrm using the SMC as a check.
        """
        r = 1e-8  # Per generation recombination rate.
        L = np.linspace(100, 10**5, 10).astype(int)
        Ne = 10**4
        n = 100
        rho = r * 4 * Ne * (L - 1)
        num_replicates = 10000
        num_trees = np.zeros(num_replicates)
        mean_exact = np.zeros_like(rho)
        var_exact = np.zeros_like(rho)
        mean_smc = np.zeros_like(rho)
        var_smc = np.zeros_like(rho)
        mean_smc_prime = np.zeros_like(rho)
        var_smc_prime = np.zeros_like(rho)
        mean_scrm = np.zeros_like(rho)
        var_scrm = np.zeros_like(rho)

        for j in range(len(L)):
            # Run SCRM under the SMC to see if we get the correct variance.
            cmd = "{} {} -L -r {} {} -l 0".format(n, num_replicates, rho[j], L[j])
            T = self.get_scrm_num_trees(
                self._scrm_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates)
            mean_scrm[j] = np.mean(T)
            var_scrm[j] = np.var(T)
            # IMPORTANT!! We have to use the get_num_breakpoints method
            # on the simulator as there is a significant drop in the number
            # of trees if we use the tree sequence. There is a significant
            # number of common ancestor events that result in a recombination
            # being undone.
            exact_sim = msprime.simulator_factory(
                sample_size=n, recombination_rate=r, Ne=Ne, length=L[j])
            for k in range(num_replicates):
                exact_sim.run()
                num_trees[k] = exact_sim.num_breakpoints
                exact_sim.reset()
            mean_exact[j] = np.mean(num_trees)
            var_exact[j] = np.var(num_trees)

            smc_sim = msprime.simulator_factory(
                sample_size=n, recombination_rate=r, Ne=Ne, length=L[j],
                model="smc")
            for k in range(num_replicates):
                smc_sim.run()
                num_trees[k] = smc_sim.num_breakpoints
                smc_sim.reset()
            mean_smc[j] = np.mean(num_trees)
            var_smc[j] = np.var(num_trees)

            smc_prime_sim = msprime.simulator_factory(
                sample_size=n, recombination_rate=r, Ne=Ne, length=L[j],
                model="smc_prime")
            for k in range(num_replicates):
                smc_prime_sim.run()
                num_trees[k] = smc_prime_sim.num_breakpoints
                smc_prime_sim.reset()
            mean_smc_prime[j] = np.mean(num_trees)
            var_smc_prime[j] = np.var(num_trees)

        basedir = "tmp__NOBACKUP__/smc_num_trees"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        pyplot.plot(rho, mean_exact, "o", label="msprime (hudson)")
        pyplot.plot(rho, mean_smc, "^", label="msprime (smc)")
        pyplot.plot(rho, mean_smc_prime, "*", label="msprime (smc_prime)")
        pyplot.plot(rho, mean_scrm, "x", label="scrm")
        pyplot.plot(rho, rho * harmonic_number(n - 1), "-")
        pyplot.legend(loc="upper left")
        pyplot.xlabel("scaled recombination rate rho")
        pyplot.ylabel("Mean number of breakpoints")
        filename = os.path.join(basedir, "mean.png")
        pyplot.savefig(filename)
        pyplot.close('all')

        v = np.zeros(len(rho))
        for j in range(len(rho)):
            v[j] = get_predicted_variance(n, rho[j])
        pyplot.plot(rho, var_exact, "o", label="msprime (hudson)")
        pyplot.plot(rho, var_smc, "^", label="msprime (smc)")
        pyplot.plot(rho, var_smc_prime, "*", label="msprime (smc_prime)")
        pyplot.plot(rho, var_scrm, "x", label="scrm")
        pyplot.plot(rho, v, "-")
        pyplot.xlabel("scaled recombination rate rho")
        pyplot.ylabel("variance in number of breakpoints")
        pyplot.legend(loc="upper left")
        filename = os.path.join(basedir, "var.png")
        pyplot.savefig(filename)
        pyplot.close('all')

    def run_simulate_from_single_locus(self):
        num_replicates = 1000

        basedir = "tmp__NOBACKUP__/simulate_from_single_locus"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        for n in [10, 50, 100, 200]:
            print("running for n =", n)
            T1 = np.zeros(num_replicates)
            reps = msprime.simulate(n, num_replicates=num_replicates)
            for j, ts in enumerate(reps):
                T1[j] = np.max(ts.tables.nodes.time)

            for t in [0.5, 1, 1.5, 5]:
                T2 = np.zeros(num_replicates)
                reps = msprime.simulate(
                    n, num_replicates=num_replicates, end_time=t)
                for j, ts in enumerate(reps):
                    final_ts = msprime.simulate(
                        from_ts=ts, start_time=np.max(ts.tables.nodes.time))
                    final_ts = final_ts.simplify()
                    T2[j] = np.max(final_ts.tables.nodes.time)

                sm.graphics.qqplot(T1)
                sm.qqplot_2samples(T1, T2, line="45")
                filename = os.path.join(basedir, "T_mrca_n={}_t={}.png".format(n, t))
                pyplot.savefig(filename, dpi=72)
                pyplot.close('all')

    def run_simulate_from_multi_locus(self):
        num_replicates = 1000
        n = 100

        basedir = "tmp__NOBACKUP__/simulate_from_multi_locus"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        for m in [10, 50, 100, 1000]:
            print("running for m =", m)
            T1 = np.zeros(num_replicates)
            num_trees1 = np.zeros(num_replicates)
            recomb_map = msprime.RecombinationMap.uniform_map(1, 1, num_loci=m)
            reps = msprime.simulate(
                n, recombination_map=recomb_map, num_replicates=num_replicates)
            for j, ts in enumerate(reps):
                T1[j] = np.max(ts.tables.nodes.time)
                num_trees1[j] = ts.num_trees

            for t in [0.5, 1, 1.5, 5]:
                T2 = np.zeros(num_replicates)
                num_trees2 = np.zeros(num_replicates)
                reps = msprime.simulate(
                    n, num_replicates=num_replicates,
                    recombination_map=recomb_map, end_time=t)
                for j, ts in enumerate(reps):
                    final_ts = msprime.simulate(
                        from_ts=ts,
                        recombination_map=recomb_map,
                        start_time=np.max(ts.tables.nodes.time))
                    final_ts = final_ts.simplify()
                    T2[j] = np.max(final_ts.tables.nodes.time)
                    num_trees2[j] = final_ts.num_trees

                sm.graphics.qqplot(T1)
                sm.qqplot_2samples(T1, T2, line="45")
                filename = os.path.join(basedir, "T_mrca_m={}_t={}.png".format(m, t))
                pyplot.savefig(filename, dpi=72)
                pyplot.close('all')

                sm.graphics.qqplot(num_trees1)
                sm.qqplot_2samples(num_trees1, num_trees2, line="45")
                filename = os.path.join(basedir, "num_trees_m={}_t={}.png".format(m, t))
                pyplot.savefig(filename, dpi=72)
                pyplot.close('all')

    def run_simulate_from_recombination(self):
        num_replicates = 1000
        n = 100
        recombination_rate = 10

        basedir = "tmp__NOBACKUP__/simulate_from_recombination"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        T1 = np.zeros(num_replicates)
        num_trees1 = np.zeros(num_replicates)
        num_edges1 = np.zeros(num_replicates)
        num_nodes1 = np.zeros(num_replicates)
        reps = msprime.simulate(
            n, recombination_rate=recombination_rate, num_replicates=num_replicates)
        for j, ts in enumerate(reps):
            T1[j] = np.max(ts.tables.nodes.time)
            num_trees1[j] = ts.num_trees
            num_nodes1[j] = ts.num_nodes
            num_edges1[j] = ts.num_edges

        print(
            "original\tmean trees = ", np.mean(num_trees1),
            "\tmean nodes = ", np.mean(num_nodes1),
            "\tmean edges = ", np.mean(num_edges1))

        for t in [0.5, 1.0, 1.5, 5.0]:
            T2 = np.zeros(num_replicates)
            num_trees2 = np.zeros(num_replicates)
            num_nodes2 = np.zeros(num_replicates)
            num_edges2 = np.zeros(num_replicates)
            reps = msprime.simulate(
                n, num_replicates=num_replicates,
                recombination_rate=recombination_rate, end_time=t)
            for j, ts in enumerate(reps):
                final_ts = msprime.simulate(
                    from_ts=ts,
                    recombination_rate=recombination_rate,
                    start_time=np.max(ts.tables.nodes.time))
                assert max(t.num_roots for t in final_ts.trees()) == 1
                final_ts = final_ts.simplify()
                T2[j] = np.max(final_ts.tables.nodes.time)
                num_trees2[j] = final_ts.num_trees
                num_nodes2[j] = final_ts.num_nodes
                num_edges2[j] = final_ts.num_edges
            print(
                "t = ", t, "\tmean trees = ", np.mean(num_trees2),
                "\tmean nodes = ", np.mean(num_nodes2),
                "\tmean edges = ", np.mean(num_edges2))

            sm.graphics.qqplot(T1)
            sm.qqplot_2samples(T1, T2, line="45")
            filename = os.path.join(basedir, "T_mrca_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_trees1)
            sm.qqplot_2samples(num_trees1, num_trees2, line="45")
            filename = os.path.join(basedir, "num_trees_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_edges1)
            sm.qqplot_2samples(num_edges1, num_edges2, line="45")
            filename = os.path.join(basedir, "num_edges_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_nodes1)
            sm.qqplot_2samples(num_nodes1, num_nodes2, line="45")
            filename = os.path.join(basedir, "num_nodes_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

    def run_simulate_from_demography(self):
        # TODO this test is considerably complicated by the fact that we
        # can't compare migrations without having support in simplify.
        # When simplify with migrations support is added, also add a test
        # here to check that the number of migrations is equivalent.
        # It's  still a good check to have the underlying numbers of
        # events reported though, so keep these now that it's implemented.
        num_replicates = 1000
        n = 50
        recombination_rate = 10
        samples = [msprime.Sample(time=0, population=j % 2) for j in range(n)]
        population_configurations = [
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration()]
        migration_matrix = [[0, 1], [1, 0]]
        demographic_events = [
            msprime.SimpleBottleneck(time=5.1, population=0, proportion=0.4),
            msprime.SimpleBottleneck(time=10.1, population=1, proportion=0.4),
            msprime.SimpleBottleneck(time=15.1, population=1, proportion=0.4),
            msprime.SimpleBottleneck(time=25.1, population=0, proportion=0.4)]

        basedir = "tmp__NOBACKUP__/simulate_from_demography"
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        T1 = np.zeros(num_replicates)
        num_ca_events1 = np.zeros(num_replicates)
        num_re_events1 = np.zeros(num_replicates)
        num_mig_events1 = np.zeros(num_replicates)
        num_trees1 = np.zeros(num_replicates)
        num_edges1 = np.zeros(num_replicates)
        num_nodes1 = np.zeros(num_replicates)

        sim = msprime.simulator_factory(
            samples=samples,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events,
            recombination_rate=recombination_rate)
        print("t\ttrees\tnodes\tedges\tca\tre\tmig")
        for j in range(num_replicates):
            sim.run()
            ts = sim.get_tree_sequence()
            num_ca_events1[j] = sim.num_common_ancestor_events
            num_re_events1[j] = sim.num_recombination_events
            num_mig_events1[j] = sum(
                [r for row in sim.num_migration_events for r in row])
            T1[j] = np.max(ts.tables.nodes.time)
            num_trees1[j] = ts.num_trees
            num_nodes1[j] = ts.num_nodes
            num_edges1[j] = ts.num_edges
            sim.reset()

        print(
            "{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                -1,
                np.mean(num_trees1),
                np.mean(num_nodes1),
                np.mean(num_edges1),
                np.mean(num_ca_events1),
                np.mean(num_re_events1),
                np.mean(num_mig_events1)))

        for t in [5.0, 10.0, 15.0, 25.0]:
            T2 = np.zeros(num_replicates)
            num_trees2 = np.zeros(num_replicates)
            num_nodes2 = np.zeros(num_replicates)
            num_edges2 = np.zeros(num_replicates)
            num_ca_events2 = np.zeros(num_replicates)
            num_re_events2 = np.zeros(num_replicates)
            num_mig_events2 = np.zeros(num_replicates)
            sim = msprime.simulator_factory(
                samples=samples,
                end_time=t,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix,
                demographic_events=demographic_events,
                recombination_rate=recombination_rate)
            for j in range(num_replicates):
                sim.run()
                ts = sim.get_tree_sequence()
                num_ca_events2[j] = sim.num_common_ancestor_events
                num_re_events2[j] = sim.num_recombination_events
                num_mig_events2[j] = sum(
                    [r for row in sim.num_migration_events for r in row])
                sim.reset()

                max_time = max(node.time for node in ts.nodes())
                sim2 = msprime.simulator_factory(
                    from_ts=ts,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=[
                        e for e in demographic_events if e.time > max_time],
                    recombination_rate=recombination_rate)
                sim2.run()

                num_ca_events2[j] += sim2.num_common_ancestor_events
                num_re_events2[j] += sim2.num_recombination_events
                num_mig_events2[j] += sum(
                    [r for row in sim2.num_migration_events for r in row])

                final_ts = sim2.get_tree_sequence().simplify()
                T2[j] = np.max(final_ts.tables.nodes.time)
                num_trees2[j] = final_ts.num_trees
                num_nodes2[j] = final_ts.num_nodes
                num_edges2[j] = final_ts.num_edges
                sim.reset()

            print(
                "{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                    t,
                    np.mean(num_trees2),
                    np.mean(num_nodes2),
                    np.mean(num_edges2),
                    np.mean(num_ca_events2),
                    np.mean(num_re_events2),
                    np.mean(num_mig_events2)))

            sm.graphics.qqplot(T1)
            sm.qqplot_2samples(T1, T2, line="45")
            filename = os.path.join(basedir, "T_mrca_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_trees1)
            sm.qqplot_2samples(num_trees1, num_trees2, line="45")
            filename = os.path.join(basedir, "num_trees_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_edges1)
            sm.qqplot_2samples(num_edges1, num_edges2, line="45")
            filename = os.path.join(basedir, "num_edges_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_nodes1)
            sm.qqplot_2samples(num_nodes1, num_nodes2, line="45")
            filename = os.path.join(basedir, "num_nodes_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_ca_events1)
            sm.qqplot_2samples(num_ca_events1, num_ca_events2, line="45")
            filename = os.path.join(basedir, "num_ca_events_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_re_events1)
            sm.qqplot_2samples(num_re_events1, num_re_events2, line="45")
            filename = os.path.join(basedir, "num_re_events_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

            sm.graphics.qqplot(num_mig_events1)
            sm.qqplot_2samples(num_mig_events1, num_mig_events2, line="45")
            filename = os.path.join(basedir, "num_mig_events_t={}.png".format(t))
            pyplot.savefig(filename, dpi=72)
            pyplot.close('all')

    def run_simulate_from_benchmark(self):
        # A quick benchmark to show this running on a large example
        L = 50 * 10**6
        seed = 3
        for n in [10**3, 10**4, 10**5]:
            print("====================")
            print("n = ", n)
            print("====================")
            before = time.perf_counter()
            ts = msprime.simulate(
                n, recombination_rate=1e-8, Ne=10**4, length=L, random_seed=seed)
            duration = time.perf_counter() - before

            print("Full sim required {:.2f} sec".format(duration))

            before = time.perf_counter()
            t = ts.tables.nodes.time[-1] / 100
            ts = msprime.simulate(
                n, recombination_rate=1e-8, Ne=10**4, length=L, random_seed=seed,
                end_time=t)
            duration = time.perf_counter() - before
            print("Initial sim required {:.2f} sec".format(duration))
            roots = np.array([tree.num_roots for tree in ts.trees()])
            print("\t", roots.shape[0], "trees, mean roots = ", np.mean(roots))
            before = time.perf_counter()
            msprime.simulate(
                from_ts=ts, recombination_rate=1e-8, Ne=10**4, length=L,
                random_seed=seed)
            duration = time.perf_counter() - before
            print("Final sim required {:.2f} sec".format(duration))

    def run_dtwf_coalescent_comparison(self, test_name, **kwargs):
        df = pd.DataFrame()
        for model in ["hudson", "dtwf"]:
            kwargs["model"] = model
            print("Running: ", kwargs)
            data = collections.defaultdict(list)
            replicates = msprime.simulate(**kwargs)
            for ts in replicates:
                t_mrca = np.zeros(ts.num_trees)
                for tree in ts.trees():
                    t_mrca[tree.index] = tree.time(tree.root)
                data["tmrca_mean"].append(np.mean(t_mrca))
                data["num_trees"].append(ts.num_trees)
                data["model"].append(model)
            df = df.append(pd.DataFrame(data))

        basedir = os.path.join("tmp__NOBACKUP__", test_name)
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        df_hudson = df[df.model == "hudson"]
        df_dtwf = df[df.model == "dtwf"]
        for stat in ["tmrca_mean", "num_trees"]:
            v1 = df_hudson[stat]
            v2 = df_dtwf[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = os.path.join(basedir, "{}.png".format(stat))
            pyplot.savefig(f, dpi=72)
            pyplot.close('all')

    def add_dtwf_vs_coalescent_single_locus(self):
        """
        Checks the DTWF against the standard coalescent at a single locus.
        """
        def f():
            self.run_dtwf_coalescent_comparison(
                "dtwf_vs_coalescent_single_locus", sample_size=10, Ne=1000,
                num_replicates=300)
        self._instances["dtwf_vs_coalescent_single_locus"] = f

    def add_dtwf_vs_coalescent_low_recombination(self):
        """
        Checks the DTWF against the standard coalescent at a single locus.
        """
        def f():
            self.run_dtwf_coalescent_comparison(
                "dtwf_vs_coalescent_low_recombination", sample_size=10, Ne=1000,
                num_replicates=400, recombination_rate=0.01)
        self._instances["dtwf_vs_coalescent_low_recombination"] = f

    def add_dtwf_vs_coalescent(
            self, key, initial_sizes, sample_sizes, num_loci, recombination_rate,
            migration_matrix=None, growth_rates=None, num_replicates=None):
        """
        Generic test of DTWF vs hudson coalescent. Populations are not
        allowed to shrink to fewer than 100 individuals, and if starting with
        fewer than 100 have growth rate set to zero.
        """
        assert len(sample_sizes) == len(initial_sizes)
        num_pops = len(sample_sizes)

        if num_replicates is None:
            num_replicates = 200

        if growth_rates is None:
            default_growth_rate = 0.01
            growth_rates = [default_growth_rate] * num_pops

        population_configurations = []
        demographic_events = []

        for i in range(num_pops):
            if initial_sizes[i] > 100:
                # Growth rate set to zero at pop size 100
                t_100 = (np.log(initial_sizes[i]) - np.log(100)) / growth_rates[i]
                de = msprime.PopulationParametersChange(
                        t_100, growth_rate=0, population=i)
                demographic_events.append(de)

                growth_rate = growth_rates[i]
            else:
                # Enforce zero growth rate for small populations
                print("Warning - setting growth rate to zero for small",
                      "population of size", initial_sizes[i])
                growth_rate = 0

            population_configurations.append(
                    msprime.PopulationConfiguration(
                        sample_size=sample_sizes[i],
                        initial_size=initial_sizes[i],
                        growth_rate=growth_rate
                        )
                    )

        recombination_map = msprime.RecombinationMap(
                [0, num_loci], [recombination_rate, 0], num_loci=num_loci)

        if migration_matrix is None:
            default_mig_rate = 0.05
            migration_matrix = []
            for i in range(num_pops):
                row = [default_mig_rate] * num_pops
                row[i] = 0
                migration_matrix.append(row)

        def f():
            self.run_dtwf_coalescent_comparison(
                    key,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    num_replicates=num_replicates,
                    demographic_events=demographic_events,
                    recombination_map=recombination_map
                    )
        self._instances[key] = f

    def add_dtwf_vs_coalescent_2_pops_massmigration(self):
        population_configurations = [
            msprime.PopulationConfiguration(sample_size=10, initial_size=1000),
            msprime.PopulationConfiguration(sample_size=10, initial_size=1000)]
        recombination_map = msprime.RecombinationMap(
                [0, int(1e6)], [1e-8, 0], num_loci=int(1e8))
        demographic_events = [
            msprime.MassMigration(
                time=300, source=1, destination=0, proportion=1.0)]

        def f():
            self.run_dtwf_coalescent_comparison(
                "dtwf_vs_coalescent_2_pops_massmigrations",
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                # Ne=0.5,
                num_replicates=300,
                recombination_map=recombination_map)
        self._instances["dtwf_vs_coalescent_2_pops_massmigration"] = f

    def add_dtwf_vs_coalescent_2_pop_growth(self):
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=10, initial_size=1000, growth_rate=0.01)]
        recombination_map = msprime.RecombinationMap(
                [0, int(5e7)], [1e-8, 0], num_loci=int(5e7))

        def f():
            self.run_dtwf_coalescent_comparison(
                "dtwf_vs_coalescent_2_pop_growth",
                population_configurations=population_configurations,
                recombination_map=recombination_map,
                num_replicates=300)
        self._instances["dtwf_vs_coalescent_2_pop_growth"] = f

    def add_dtwf_vs_coalescent_2_pop_shrink(self):
        initial_size = 1000

        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=10, initial_size=initial_size, growth_rate=-0.01)]
        recombination_map = msprime.RecombinationMap(
            [0, int(1e7)], [1e-8, 0], num_loci=int(1e7))
        demographic_events = [
            msprime.PopulationParametersChange(
                time=200, initial_size=initial_size, growth_rate=0.01, population_id=0)
        ]

        def f():
            self.run_dtwf_coalescent_comparison(
                "dtwf_vs_coalescent_2_pop_shrink",
                population_configurations=population_configurations,
                recombination_map=recombination_map,
                demographic_events=demographic_events,
                num_replicates=300)
        self._instances["dtwf_vs_coalescent_2_pop_shrink"] = f

    def add_dtwf_vs_coalescent_multiple_bottleneck(self):
        population_configurations = [
            msprime.PopulationConfiguration(sample_size=5, initial_size=1000),
            msprime.PopulationConfiguration(sample_size=5, initial_size=1000)]
        recombination_map = msprime.RecombinationMap(
                [0, int(1e6)], [1e-8, 0], num_loci=int(1e8))
        # migration_matrix = [[0, 0.1], [0.1, 0]]

        demographic_events = [
            msprime.PopulationParametersChange(
                    time=100, initial_size=100, growth_rate=-0.01, population_id=0),
            msprime.PopulationParametersChange(
                    time=200, initial_size=100, growth_rate=-0.01, population_id=1),
            msprime.PopulationParametersChange(
                    time=300, initial_size=1000, growth_rate=0.01, population_id=0),
            msprime.PopulationParametersChange(
                    time=400, initial_size=1000, growth_rate=0.01, population_id=1),
            msprime.PopulationParametersChange(
                    time=500, initial_size=100, growth_rate=0, population_id=0),
            msprime.PopulationParametersChange(
                    time=600, initial_size=100, growth_rate=0, population_id=1),
            msprime.MigrationRateChange(
                    time=700, rate=0.1, matrix_index=(0, 1))
            ]

        def f():
            self.run_dtwf_coalescent_comparison(
                "dtwf_vs_coalescent_multiple_bottleneck",
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                # migration_matrix=migration_matrix,
                num_replicates=400,
                recombination_map=recombination_map)
        self._instances["dtwf_vs_coalescent_multiple_bottleneck"] = f

    def add_dtwf_vs_coalescent_random_instance(
            self, key, num_populations=1, num_replicates=200, num_demographic_events=0):

        N = num_populations
        num_loci = np.random.randint(1e5, 1e7)
        rho = 1e-8
        recombination_map = msprime.RecombinationMap(
                [0, num_loci], [rho, 0], num_loci=num_loci)

        population_configurations = []
        for i in range(N):
            population_configurations.append(
                msprime.PopulationConfiguration(
                    sample_size=np.random.randint(1, 10),
                    initial_size=int(1000 / N)))

        migration_matrix = []
        for i in range(N):
            migration_matrix.append(
                [random.uniform(0.05, 0.25) * (j != i) for j in range(N)])

        # Add demographic events and some migration rate changes
        t_max = 1000
        demographic_events = []
        times = sorted(np.random.randint(300, t_max, size=num_demographic_events))
        for t in times:
            initial_size = np.random.randint(500, 1000)
            # Setting growth_rate to 0 because it's too tricky to get
            # growth_rates in the DTWF which don't result in N going to 0.
            growth_rate = 0
            pop_id = np.random.randint(N)
            demographic_events.append(
                msprime.PopulationParametersChange(
                    time=t, initial_size=initial_size,
                    growth_rate=growth_rate, population_id=pop_id))

            if random.random() < 0.5 and N >= 2:
                rate = random.uniform(0.05, 0.25)
                index = tuple(
                    np.random.choice(range(num_populations), size=2, replace=False))
                demographic_events.append(
                    msprime.MigrationRateChange(time=t, rate=rate, matrix_index=index))

        # Collect all pops together to control coalescence times for DTWF
        for i in range(1, N):
            demographic_events.append(
                msprime.MassMigration(
                    time=t_max, source=i, destination=0, proportion=1.0))

        demographic_events.append(
            msprime.PopulationParametersChange(
                time=t_max, initial_size=100, growth_rate=0, population_id=0))

        def f():
            self.run_dtwf_coalescent_comparison(
                key,
                migration_matrix=migration_matrix,
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                num_replicates=num_replicates,
                recombination_map=recombination_map)
        self._instances[key] = f

    def run_dtwf_slim_comparison(self, test_name, slim_args, **kwargs):

        df = pd.DataFrame()

        kwargs["model"] = "dtwf"
        print("Running: ", kwargs)
        replicates = msprime.simulate(**kwargs)
        data = collections.defaultdict(list)
        for ts in replicates:
            t_mrca = np.zeros(ts.num_trees)
            for tree in ts.trees():
                t_mrca[tree.index] = tree.time(tree.root)
            data["tmrca_mean"].append(np.mean(t_mrca))
            data["num_trees"].append(ts.num_trees)
            data["model"].append("dtwf")

        basedir = os.path.join("tmp__NOBACKUP__", test_name)
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        slim_script = os.path.join(basedir, "slim_script.txt")
        outfile = os.path.join(basedir, "slim.trees")
        slim_args['OUTFILE'] = outfile
        write_slim_script(slim_script, slim_args)

        cmd = self._slim_executable + [slim_script]
        for _ in tqdm.tqdm(range(kwargs['num_replicates'])):
            subprocess.check_output(cmd)
            ts = msprime.load(outfile)
            ts = subsample_simplify_slim_treesequence(ts, slim_args['sample_sizes'])

            t_mrca = np.zeros(ts.num_trees)
            for tree in ts.trees():
                t_mrca[tree.index] = tree.time(tree.root)

            data["tmrca_mean"].append(np.mean(t_mrca))
            data["num_trees"].append(ts.num_trees)
            data["model"].append("slim")
        df = df.append(pd.DataFrame(data))

        df_slim = df[df.model == "slim"]
        df_dtwf = df[df.model == "dtwf"]
        for stat in ["tmrca_mean", "num_trees"]:
            v1 = df_slim[stat]
            v2 = df_dtwf[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = os.path.join(basedir, "{}.png".format(stat))
            pyplot.xlabel("DTWF")
            pyplot.ylabel("SLiM")
            pyplot.savefig(f, dpi=72)
            pyplot.close('all')

    def add_dtwf_vs_slim(
            self, key, initial_sizes, sample_sizes, num_loci,
            recombination_rate, migration_matrix=None, num_replicates=None):
        """
        Generic test of DTWF vs SLiM WF simulator, without growth rates
        """
        self.check_slim_version()
        assert len(sample_sizes) == len(initial_sizes)

        num_pops = len(sample_sizes)
        slim_args = {}

        if num_replicates is None:
            num_replicates = 200

        slim_args['sample_sizes'] = sample_sizes

        population_configurations = []
        slim_args['POP_STRS'] = ''
        for i in range(len(sample_sizes)):
            population_configurations.append(
                    msprime.PopulationConfiguration(
                        sample_size=sample_sizes[i],
                        initial_size=initial_sizes[i],
                        growth_rate=0
                        )
                    )
            slim_args['POP_STRS'] += "sim.addSubpop('p{i}', {N});\n".format(
                    i=i, N=initial_sizes[i])

        if migration_matrix is None:
            default_mig_rate = 0.01
            migration_matrix = []
            for i in range(num_pops):
                row = [default_mig_rate] * num_pops
                row[i] = 0
                migration_matrix.append(row)

        # SLiM rates are 'immigration' forwards in time, which matches
        # DTWF backwards-time 'emmigration'
        assert(len(migration_matrix) == num_pops)
        if num_pops > 1:
            for i in range(num_pops):
                row = migration_matrix[i]
                indices = [j for j in range(num_pops) if j != i]
                pop_names = ['p' + str(j) for j in indices]
                rates = [str(row[j]) for j in indices]

                to_pop_str = ','.join(pop_names)
                rate_str = ','.join(rates)

                mig_str = "p{}.setMigrationRates(c({}), c({}));\n".format(
                        i, to_pop_str, rate_str)
                slim_args['POP_STRS'] += mig_str

        num_loci = int(num_loci)
        recombination_map = msprime.RecombinationMap(
                [0, num_loci], [recombination_rate, 0], num_loci=num_loci)
        slim_args['RHO'] = recombination_rate
        slim_args['NUM_LOCI'] = num_loci

        def f():
            self.run_dtwf_slim_comparison(
                    key, slim_args,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    num_replicates=num_replicates,
                    recombination_map=recombination_map,
                    )
        self._instances[key] = f

    def run_xi_hudson_comparison(self, test_name, xi_model, **kwargs):
        df = pd.DataFrame()
        for model in ["hudson", xi_model]:
            kwargs["model"] = model
            model_str = "hudson"
            if model != "hudson":
                model_str = "Xi"
            print("Running: ", kwargs)
            replicates = msprime.simulate(**kwargs)
            data = collections.defaultdict(list)
            for ts in replicates:
                t_mrca = np.zeros(ts.num_trees)
                for tree in ts.trees():
                    t_mrca[tree.index] = tree.time(tree.root)
                data["tmrca_mean"].append(np.mean(t_mrca))
                data["num_trees"].append(ts.num_trees)
                data["num_nodes"].append(ts.num_nodes)
                data["num_edges"].append(ts.num_edges)
                data["model"].append(model_str)
            df = df.append(pd.DataFrame(data))

        basedir = os.path.join("tmp__NOBACKUP__", test_name)
        if not os.path.exists(basedir):
            os.mkdir(basedir)

        df_hudson = df[df.model == "hudson"]
        df_xi = df[df.model == "Xi"]
        for stat in ["tmrca_mean", "num_trees", "num_nodes", "num_edges"]:
            v1 = df_hudson[stat]
            v2 = df_xi[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = os.path.join(basedir, "{}.png".format(stat))
            pyplot.savefig(f, dpi=72)
            pyplot.close('all')

    def add_xi_dirac_vs_hudson_single_locus(self):
        """
        Checks Xi-dirac against the standard coalescent at a single locus.
        """
        def f():
            N = 100
            self.run_xi_hudson_comparison(
                "xi_dirac_vs_hudson_single_locus",
                msprime.DiracCoalescent(N, psi=0.99, c=0),
                sample_size=10, Ne=N, num_replicates=5000)
        self._instances["xi_dirac_vs_hudson_single_locus"] = f

    def add_xi_dirac_vs_hudson_recombination(self):
        """
        Checks Xi-dirac against the standard coalescent with recombination.
        """
        def f():
            N = 100
            self.run_xi_hudson_comparison(
                "xi_dirac_vs_hudson_recombination",
                msprime.DiracCoalescent(N, psi=0.99, c=0),
                sample_size=50, Ne=N, num_replicates=1000,
                recombination_rate=0.1)
        self._instances["xi_dirac_vs_hudson_recombination"] = f

    def compare_xi_dirac_sfs(self, sample_size, psi, c, sfs, num_replicates=1000):
        """
        Runs simulations of the xi dirac model and compares to the expected SFS.
        """
        print("running SFS for", sample_size, psi, c)
        reps = msprime.simulate(
            sample_size, num_replicates=num_replicates,
            model=msprime.DiracCoalescent(psi=psi, c=c))

        data = collections.defaultdict(list)
        tbl_sum = [0] * (sample_size - 1)
        for j, ts in enumerate(reps):
            for tree in ts.trees():
                tot_bl = 0.0
                tbl = [0] * (sample_size - 1)
                for node in tree.nodes():
                    if tree.parent(node) != msprime.NULL_NODE:
                        tbl[tree.num_samples(node)-1] = tbl[
                            tree.num_samples(node)-1] + tree.branch_length(node)
                        tot_bl = tot_bl + tree.branch_length(node)

                for xi in range(sample_size - 1):
                    rescaled_x = tbl[xi]/tot_bl
                    data["total_branch_length"].append(rescaled_x)
                    tbl_sum[xi] = tbl_sum[xi] + rescaled_x
                data["num_leaves"].extend(range(1, sample_size))

        basedir = os.path.join("tmp__NOBACKUP__", "xi_dirac_expected_sfs")
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        f = os.path.join(basedir, "n={}_psi={}.png".format(sample_size, psi))

        ax = sns.violinplot(
            data=data, x="num_leaves", y="total_branch_length", color="grey")
        ax.set_xlabel("num leaves")
        l1 = ax.plot(np.arange(sample_size - 1), sfs[::], "--", linewidth=3)
        l2 = ax.plot(
            np.arange(sample_size - 1), [x/num_replicates for x in tbl_sum],
            "--", linewidth=3)
        ax.legend((l1[0], l2[0]), ("Expected", "Observed"))
        pyplot.savefig(f, dpi=72)
        pyplot.close('all')

    def run_xi_dirac_expected_sfs(self):
        self.compare_xi_dirac_sfs(
            num_replicates=5000, sample_size=3, psi=0.01, c=1, sfs=[0.666667, 0.333333])
        self.compare_xi_dirac_sfs(
            num_replicates=5000, sample_size=3, psi=0.99, c=1,
            sfs=[0.6722604, 0.3277396])
        self.compare_xi_dirac_sfs(
            num_replicates=5000, sample_size=4, psi=0.01, c=1,
            sfs=[0.5457826, 0.2728913, 0.1813261])
        self.compare_xi_dirac_sfs(
            num_replicates=5000, sample_size=4, psi=0.99, c=1,
            sfs=[0.5611642, 0.2747103, 0.1641255])

        # MORE, NEED TO CHECK THESE VALUES

        self.compare_xi_dirac_sfs(
            num_replicates=1000,
            sample_size=13, psi=0.5, c=1,
            sfs=[
                0.418425, 0.121938, 0.092209, 0.070954, 0.056666, 0.047179,
                0.040545, 0.035631, 0.031841, 0.028832, 0.026796, 0.028985])

    def add_xi_dirac_expected_sfs(self):
        """
        Adds a check for xi_dirac matching expected SFS calculations.
        """
        self._instances["xi_dirac_expected_sfs"] = self.run_xi_dirac_expected_sfs

    def compare_xi_beta_sfs(self, sample_size, alpha, sfs, num_replicates=1000):
        """
        Runs simulations of the xi beta model and compares to the expected SFS.
        """
        print("running SFS for", sample_size, alpha)
        reps = msprime.simulate(
            sample_size, num_replicates=num_replicates,
            model=msprime.BetaCoalescent(alpha=alpha, truncation_point=1000))

        data = collections.defaultdict(list)
        tbl_sum = [0] * (sample_size - 1)
        for j, ts in enumerate(reps):
            for tree in ts.trees():
                tot_bl = 0.0
                tbl = [0] * (sample_size - 1)
                for node in tree.nodes():
                    if tree.parent(node) != msprime.NULL_NODE:
                        tbl[tree.num_samples(node)-1] = tbl[
                            tree.num_samples(node)-1] + tree.branch_length(node)
                        tot_bl = tot_bl + tree.branch_length(node)

                for xi in range(sample_size - 1):
                    rescaled_x = tbl[xi]/tot_bl
                    data["total_branch_length"].append(rescaled_x)
                    tbl_sum[xi] = tbl_sum[xi] + rescaled_x
                data["num_leaves"].extend(range(1, sample_size))

        basedir = os.path.join("tmp__NOBACKUP__", "xi_beta_expected_sfs")
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        f = os.path.join(basedir, "n={}_alpha={}.png".format(sample_size, alpha))

        ax = sns.violinplot(
            data=data, x="num_leaves", y="total_branch_length", color="grey")
        ax.set_xlabel("num leaves")
        l1 = ax.plot(np.arange(sample_size - 1), sfs[::], "--", linewidth=3)
        l2 = ax.plot(
            np.arange(sample_size - 1), [x/num_replicates for x in tbl_sum],
            "--", linewidth=3)
        ax.legend((l1[0], l2[0]), ("Expected", "Observed"))
        pyplot.savefig(f, dpi=72)
        pyplot.close('all')

    def run_xi_beta_expected_sfs(self):
        self.compare_xi_beta_sfs(
            num_replicates=5000,
            sample_size=3, alpha=1.01, sfs=[0.681653, 0.318347])

        self.compare_xi_beta_sfs(
            num_replicates=5000,
            sample_size=3, alpha=1.8, sfs=[0.6694913, 0.3305087])

        self.compare_xi_beta_sfs(
            num_replicates=5000,
            sample_size=4, alpha=1.01, sfs=[0.5684275, 0.2576535, 0.1739190])

        self.compare_xi_beta_sfs(
            num_replicates=5000,
            sample_size=4, alpha=1.8, sfs=[0.5501309, 0.2691312, 0.1807379])

        # MORE

        self.compare_xi_beta_sfs(
            num_replicates=1000,
            sample_size=13, alpha=1.01,
            sfs=[0.400253, 0.134518, 0.093954, 0.072698, 0.058500, 0.048636,
                 0.041617, 0.036404, 0.032334, 0.028913, 0.026112, 0.026060])

    def add_xi_beta_expected_sfs(self):
        """
        Adds a check for xi_beta matching expected SFS calculations.
        """
        self._instances["xi_beta_expected_sfs"] = self.run_xi_beta_expected_sfs

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

    def add_mean_coaltime_check(self):
        """
        Adds a check for the demography debugger predictions about
        mean coalescence time.
        """
        self._instances["mean_coaltime"] = self.run_mean_coaltime_check

    def add_total_branch_length_analytical_check(self):
        """
        Adds a check for the analytical check for the total branch length.
        """
        self._instances["analytical_tbl"] = self.run_tbl_analytical_check

    def add_pairwise_island_model_analytical_check(self):
        """
        Adds a check for the analytical check for pairwise island model
        """
        self._instances["analytical_pairwise_island"] = self.run_pairwise_island_model

    def add_arg_recording_check(self):
        """
        Adds a check that we get the right number of objects when we simplify
        a full arg.
        """
        self._instances["arg_recording"] = self.run_arg_recording

    def add_smc_num_trees_analytical_check(self):
        """
        Adds a check for the analytical number of trees under the SMC
        and the full coalescent.
        """
        self._instances["smc_num_trees"] = self.run_smc_num_trees

    def add_cli_num_trees_analytical_check(self):
        """
        Adds a check for the analytical number of trees using the CLI
        and comparing with ms.
        """
        self._instances["cli_num_trees"] = self.run_cli_num_trees

    def add_smc_oldest_time_check(self):
        """
        Adds a check the distribution of the oldest time of a
        coalescence in the smc using scrm.
        """
        self._instances["smc_oldest_time"] = self.run_smc_oldest_time

    def add_simulate_from_single_locus_check(self):
        """
        Check that the distributions are identitical when we run simulate_from
        at various time points.
        """
        self._instances[
            "simulate_from_single_locus"] = self.run_simulate_from_single_locus

    def add_simulate_from_multi_locus_check(self):
        """
        Check that the distributions are identitical when we run simulate_from
        at various time points.
        """
        self._instances[
            "simulate_from_multi_locus"] = self.run_simulate_from_multi_locus

    def add_simulate_from_recombination_check(self):
        """
        Check that the distributions are identitical when we run simulate_from
        at various time points.
        """
        self._instances[
            "simulate_from_recombination"] = self.run_simulate_from_recombination

    def add_simulate_from_demography_check(self):
        """
        Check that the distributions are identitical when we run simulate_from
        at various time points.
        """
        self._instances[
            "simulate_from_demography"] = self.run_simulate_from_demography

    def add_simulate_from_benchmark(self):
        """
        Check that the distributions are identitical when we run simulate_from
        at various time points.
        """
        self._instances["simulate_from_benchmark"] = self.run_simulate_from_benchmark

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


def run_tests(args):
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
        "-eN 1.1 0.001 -ema 10 2 x 0 10 x")
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
        "population-split-4-pops3",
        "100 10000 -t 2.0 -I 4 25 25 25 25 -ej 1 2 1 -em 1.5 4 1 2 "
        "-ej 2 3 1 -ej 3 4 1")
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
        "admixture-2-pop3",
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -G 5.0 "
        "-em 2.0 3 1 1")
    verifier.add_ms_instance(
        "admixture-2-pop4",
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -eg 0.02 1 5.0 "
        "-em 0.02 3 1 1")

    # # Examples from ms documentation
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
        "msdoc-two-species",
        "15 10000 -t 11.2 -I 2 3 12 -g 1 44.36 -n 2 0.125 -eg 0.03125 1 0.0 "
        "-en 0.0625 2 0.05 -ej 0.09375 2 1")
    verifier.add_ms_instance(
        "msdoc-stepping-stone",
        "15 10000 -t 3.0 -I 6 0 7 0 0 8 0 -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 "
        "-m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 -em 2.0 3 4 "
        "2.5 -em 2.0 4 3 2.5")

    # The order of simultaneous events matters in ms.
    verifier.add_ms_instance(
        "simultaneous-ex1", "10 10000 -t 2.0 -eN 0.3 0.5 -eG .3 7.0")
    # Add a bunch more instances...
    verifier.add_ms_instance(
        "zero-growth-rate", "10 10000 -t 2.0 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5")
    # Some examples provided by Konrad Lohse
    verifier.add_ms_instance(
        "konrad-1",
        "4 1000 -t 2508 -I 2 2 2 0 -n 2 2.59 -ma x 0 1.502 x -ej 0.9485 1 2 "
        "-r 23.76 3000")
    verifier.add_ms_instance(
        "konrad-2",
        "3 10000 -t 0.423 -I 3 1 1 1 -es 0.0786 1 0.946635 -ej 0.0786 4 3 "
        "-ej 0.189256 1 2 -ej 0.483492 2 3")
    verifier.add_ms_instance(
        "konrad-3", "100 100 -t 2 -I 10 10 10 10 10 10 10 10 10 10 10 0.001 ")

    # Simple discoal tests
    verifier.add_discoal_instance(
        "discoal-simple-ex",
        "15 1000 100 -t 5.0")
    verifier.add_discoal_instance(
        "discoal-size-change1",
        "15 1000 100 -t 10.0 -en 0.1 0 2.0")
    verifier.add_discoal_instance(
        "discoal-size-change2",
        "15 1000 100 -t 10.0 -en 0.1 0 0.1")
    verifier.add_discoal_instance(
        "discoal-size-change3",
        "50 1000 100 -t 10.0 -en 0.01 0 0.01")
    verifier.add_discoal_instance(
        "discoal-size-change4",
        "50 1000 100 -t 10.0 -en 0.01 0 0.5 -en 0.05 0 1.0")

    # Add some random instances.
    verifier.add_random_instance("random1")
    verifier.add_random_instance(
        "random2", num_replicates=10**4, num_demographic_events=10)
    # verifier.add_random_instance("random2", num_populations=3)

    # Add analytical checks
    verifier.add_s_analytical_check()
    verifier.add_pi_analytical_check()
    verifier.add_mean_coaltime_check()
    verifier.add_total_branch_length_analytical_check()
    verifier.add_pairwise_island_model_analytical_check()
    verifier.add_cli_num_trees_analytical_check()

    # ARG recording
    verifier.add_arg_recording_check()

    # Simulate-from checks.
    verifier.add_simulate_from_single_locus_check()
    verifier.add_simulate_from_multi_locus_check()
    verifier.add_simulate_from_recombination_check()
    verifier.add_simulate_from_demography_check()
    verifier.add_simulate_from_benchmark()

    # Add SMC checks against scrm.
    verifier.add_smc_num_trees_analytical_check()
    verifier.add_smc_oldest_time_check()

    # Add XiDirac checks against standard coalescent.
    verifier.add_xi_dirac_vs_hudson_single_locus()
    verifier.add_xi_dirac_vs_hudson_recombination()
    verifier.add_xi_dirac_expected_sfs()
    verifier.add_xi_beta_expected_sfs()

    # DTWF checks against coalescent.
    verifier.add_dtwf_vs_coalescent_single_locus()
    verifier.add_dtwf_vs_coalescent_low_recombination()
    verifier.add_dtwf_vs_coalescent_2_pops_massmigration()
    verifier.add_dtwf_vs_coalescent_2_pop_growth()
    verifier.add_dtwf_vs_coalescent_2_pop_shrink()
    verifier.add_dtwf_vs_coalescent_multiple_bottleneck()

    verifier.add_dtwf_vs_coalescent(
        'dtwf_vs_coalescent_long_region', [1000], [10], int(1e8), 1e-8)
    verifier.add_dtwf_vs_coalescent(
        'dtwf_vs_coalescent_short_region', [1000], [10], int(1e6), 1e-8)
    verifier.add_dtwf_vs_coalescent(
        'dtwf_vs_coalescent_2_pops', [500, 500], [5, 5], int(1e6), 1e-8,
        num_replicates=500)
    verifier.add_dtwf_vs_coalescent(
        'dtwf_vs_coalescent_3_pops', [500, 500, 500], [5, 2, 0], int(1e7), 1e-8)
    verifier.add_dtwf_vs_coalescent(
        'dtwf_vs_coalescent_4_pops', [1000, 1000, 1000, 1000], [0, 20, 0, 0],
        int(1e6), 1e-8, num_replicates=500)

    migration_matrix = [[0, 0.2, 0.1], [0.1, 0, 0.2], [0.2, 0.1, 0]]
    verifier.add_dtwf_vs_coalescent(
        'dtwf_vs_coalescent_3_pops_asymm_mig', [500, 500, 500], [20, 0, 0],
        int(1e6), 1e-8, migration_matrix=migration_matrix, num_replicates=500)

    migration_matrix = [[0, 0.5], [0.7, 0]]
    verifier.add_dtwf_vs_coalescent(
        'dtwf_vs_coalescent_2_pops_high_asymm_mig', [1000, 1000], [10, 10],
        int(1e6), 1e-8, migration_matrix=migration_matrix, num_replicates=200,
        growth_rates=[0.005, 0.005])

    # Random checks vs Hudson coalescent - extended only
    if args.extended is True:
        verifier.add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_1",
            num_populations=2, num_replicates=200, num_demographic_events=3)
        verifier.add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_2",
            num_populations=3, num_replicates=200, num_demographic_events=3)
        verifier.add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_3",
            num_populations=2, num_replicates=200, num_demographic_events=6)
        verifier.add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_4",
            num_populations=1, num_replicates=200, num_demographic_events=8)

    # DTWF checks against SLiM
    # TODO: Add back multi-pop tests of DTWF vs. SLiM when full diploid
    # simulations are implemented.
    verifier.add_dtwf_vs_slim(
        'dtwf_vs_slim_single_locus', [10], [10], 1, 0)
    verifier.add_dtwf_vs_slim(
        'dtwf_vs_slim_short_region', [100], [10], 1e7, 1e-8, num_replicates=200)
    verifier.add_dtwf_vs_slim(
        'dtwf_vs_slim_long_region', [50], [10], 1e8, 1e-8, num_replicates=200)

    keys = None
    if len(args.keys) > 0:
        keys = args.keys

    verifier.run(keys)


def add_simulator_arguments(parser):
    parser.add_argument("--extended", action="store_true")
    parser.add_argument("keys", nargs="*")


def main():
    parser = argparse.ArgumentParser()
    add_simulator_arguments(parser)
    args = parser.parse_args()
    run_tests(args)


if __name__ == "__main__":
    main()
