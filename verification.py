"""
Script to automate verification of the msprime simulator against
known statistical results and benchmark programs such as Hudson's ms.
"""
import argparse
import ast
import collections
import concurrent.futures
import itertools
import logging
import math
import os
import random
import subprocess
import sys
import tempfile
import time

import allel
import attr
import daiquiri
import dendropy
import matplotlib
import numpy as np
import pandas as pd
import pyvolve
import scipy.special
import scipy.stats
import seaborn as sns
import tqdm
import tskit
from matplotlib import lines as mlines
from matplotlib import pyplot

import msprime
import msprime.cli as cli
from msprime.demography import _matrix_exponential

# Force matplotlib to not use any Xwindows backend.
# Note this must be done before importing statsmodels.
matplotlib.use("Agg")
import statsmodels.api as sm  # noqa: E402


def flatten(l):
    return [x for sublist in l for x in sublist]


def scale_breakpoints(df, factor):
    def scale(points):
        return [factor * x for x in points]

    df["breakpoints"] = df["breakpoints"].map(scale)


def harmonic_number(n):
    return np.sum(1 / np.arange(1, n + 1))


def hk_f(n, z):
    """
    Returns Hudson and Kaplan's f_n(z) function. This is based on the exact
    value for n=2 and the approximations given in the 1985 Genetics paper.
    """
    ret = 0
    if n == 2:
        ret = (18 + z) / (z ** 2 + 13 * z + 18)
    else:
        ret = sum(1 / j ** 2 for j in range(1, n)) * hk_f(2, z)
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
    with open(outfile, "w") as f:
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


def run_dtwf_coalescent_stats(**kwargs):
    df = pd.DataFrame()

    for model in ["hudson", "dtwf"]:
        kwargs["model"] = model

        logging.debug(f"Running: {kwargs}")
        data = collections.defaultdict(list)
        replicates = msprime.simulate(**kwargs)
        for ts in replicates:
            t_mrca = np.zeros(ts.num_trees)
            t_intervals = []
            for tree in ts.trees():
                t_mrca[tree.index] = tree.time(tree.root)
                t_intervals.append(tree.interval)
            data["tmrca_mean"].append(np.mean(t_mrca))
            data["num_trees"].append(ts.num_trees)
            data["intervals"].append(t_intervals)
            data["model"].append(model)
        df = df.append(pd.DataFrame(data))
    return df


def plot_qq(v1, v2):
    sm.graphics.qqplot(v1)
    sm.qqplot_2samples(v1, v2, line="45")


def plot_breakpoints_hist(v1, v2, v1_name, v2_name):
    sns.kdeplot(v1, color="b", label=v1_name, shade=True, legend=False)
    sns.kdeplot(v2, color="r", label=v2_name, shade=True, legend=False)
    pyplot.legend(loc="upper right")


def all_breakpoints_in_replicates(replicates):
    return [right for intervals in replicates for left, right in intervals]


def plot_dtwf_coalescent_stats(basedir, df):
    df_hudson = df[df.model == "hudson"]
    df_dtwf = df[df.model == "dtwf"]
    for stat in ["tmrca_mean", "num_trees"]:
        plot_qq(df_hudson[stat], df_dtwf[stat])
        f = os.path.join(basedir, "{stat}.png")
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")

    hudson_breakpoints = all_breakpoints_in_replicates(df_hudson["intervals"])
    dtwf_breakpoints = all_breakpoints_in_replicates(df_dtwf["intervals"])
    if len(hudson_breakpoints) > 0 or len(dtwf_breakpoints) > 0:
        plot_breakpoints_hist(hudson_breakpoints, dtwf_breakpoints, "hudson", "dtwf")
        f = os.path.join(basedir, "breakpoints.png")
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")


def plot_tree_intervals(basedir, df):
    fig, ax_arr = pyplot.subplots(2, 1)
    for subplot_idx, model in enumerate(["hudson", "dtwf"]):
        intervals = df[df.model == model]["intervals"][0]
        for i, interval in enumerate(intervals):
            left, right = interval
            ax_arr[subplot_idx].set_title(model)
            ax_arr[subplot_idx].set_ylabel("tree index")
            ax_arr[subplot_idx].plot([left, right], [i, i], c="grey")

    ax_arr[1].set_xlabel("tree interval")
    pyplot.tight_layout()
    f = os.path.join(basedir, "breakpoints.png")
    pyplot.savefig(f, dpi=72)
    pyplot.close("all")


def make_test_dir(test_name):
    root = "tmp__NOBACKUP__"
    if not os.path.exists(root):
        os.mkdir(root)
    basedir = os.path.join(root, test_name)
    if not os.path.exists(basedir):
        os.mkdir(basedir)
    return basedir


@attr.s
class Test:
    """
    The superclass of all tests. Each test must define a name, a group
    and a run method.
    """

    name = attr.ib(type=str)
    group = attr.ib(type=str)
    _mspms_executable = attr.ib(init=False, default=[sys.executable, "mspms_dev.py"])
    _slim_executable = attr.ib(init=False, default=["./data/slim"])
    _ms_executable = attr.ib(init=False, default=["./data/ms"])

    def run(self, output_dir):
        raise NotImplementedError()

    def _run_sample_stats(self, args):
        logging.debug(f"\t {' '.join(args)}")
        p1 = subprocess.Popen(args, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["./data/sample_stats"], stdin=p1.stdout, stdout=subprocess.PIPE
        )
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

    def _build_filename(self, *args):
        output_dir = os.path.join(self._output_dir, args[0])
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        return os.path.join(output_dir, "_".join(args[1:]))

    def _plot_stats(self, key, stats_type, df1, df2, df1_name, df2_name):
        assert set(df1.columns.values) == set(df2.columns.values)
        for stat in df1.columns.values:
            v1 = df1[stat]
            v2 = df2[stat]
            if stat == "breakpoints":
                plot_breakpoints_hist(flatten(v1), flatten(v2), df1_name, df2_name)
            else:
                plot_qq(v1, v2)
            f = self._build_filename(key, stats_type, stat)
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def get_ms_seeds(self):
        max_seed = 2 ** 16
        seeds = [random.randint(1, max_seed) for j in range(3)]
        return ["-seed"] + list(map(str, seeds))

    def _run_msprime_mutation_stats(self, args):
        return self._run_sample_stats(
            self._mspms_executable + args.split() + self.get_ms_seeds()
        )


@attr.s
class MsTest(Test):

    # functions common in ms and random
    def _deserialize_breakpoints(self, df):
        breakpoints_strs = df["breakpoints"]
        breakpoints = [ast.literal_eval(l) for l in breakpoints_strs]
        df["breakpoints"] = breakpoints
        return df

    def _exec_coalescent_stats(self, executable, args, seeds=None):
        with tempfile.TemporaryFile() as f:
            argList = [executable] + args.split() + self.get_ms_seeds()
            logging.debug(f"\t {' '.join(argList)}")
            subprocess.call(argList, stdout=f)
            f.seek(0)
            df = pd.read_table(f)
        self._deserialize_breakpoints(df)
        return df

    def _run_ms_coalescent_stats(self, args):
        return self._exec_coalescent_stats("./data/ms_summary_stats", args)

    def _run_ms_mutation_stats(self, args):
        return self._run_sample_stats(
            self._ms_executable + args.split() + self.get_ms_seeds()
        )

    def _run_mutation_stats(self, key, args):
        df_ms = self._run_ms_mutation_stats(args)
        df_msp = self._run_msprime_mutation_stats(args)
        self._plot_stats(key, "mutation", df_ms, df_msp, "ms", "msp")

    def _run_mspms_coalescent_stats(self, args):
        logging.debug(f"\t mspms: {args}")
        runner = cli.get_mspms_runner(args.split())
        sim = runner.get_simulator()
        num_populations = sim.num_populations
        replicates = runner.get_num_replicates()
        num_trees = [0 for j in range(replicates)]
        time = [0 for j in range(replicates)]
        ca_events = [0 for j in range(replicates)]
        re_events = [0 for j in range(replicates)]
        gc_events = [0 for j in range(replicates)]
        mig_events = [None for j in range(replicates)]
        breakpoints = [[] for j in range(replicates)]
        for j in range(replicates):
            sim.reset()
            sim.run()
            num_trees[j] = sim.num_breakpoints + 1
            breakpoints[j] = sim.breakpoints
            time[j] = sim.time
            ca_events[j] = sim.num_common_ancestor_events
            re_events[j] = sim.num_recombination_events
            gc_events[j] = sim.num_gene_conversion_events
            mig_events[j] = [r for row in sim.num_migration_events for r in row]
        d = {
            "t": time,
            "num_trees": num_trees,
            "ca_events": ca_events,
            "re_events": re_events,
            "gc_events": gc_events,
        }
        for j in range(num_populations ** 2):
            events = [mig_events[k][j] for k in range(replicates)]
            d[f"mig_events_{j}"] = events
        d["breakpoints"] = breakpoints
        df = pd.DataFrame(d)
        return df

    def _run_coalescent_stats(self, key, args):
        df_msp = self._run_mspms_coalescent_stats(args)
        df_ms = self._run_ms_coalescent_stats(args)
        self._plot_stats(key, "coalescent", df_msp, df_ms, "msp", "ms")

    # end of tests common to MS and random
    def _run_variable_recombination_coalescent_stats(self, key, args):
        df_msp = self._run_mspms_coalescent_stats(args)
        df_mshot = self._run_mshot_coalescent_stats(args)
        self._plot_stats(key, "recomb map coalescent", df_msp, df_mshot, "msp", "msHOT")

    def _run_mshot_coalescent_stats(self, args):
        return self._exec_coalescent_stats("./data/msHOT_summary_stats", args)

    def _run_ms_mshot_stats(self, key, args):
        df_ms = self._run_ms_coalescent_stats(args)
        df_mshot = self._run_mshot_coalescent_stats(args)
        self._plot_stats(key, "ms mshot consistency", df_mshot, df_ms, "msHOT", "ms")


@attr.s
class MsTest1(MsTest):
    command = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        logging.debug(f"{self.name} {self.command}")
        self._run_coalescent_stats(self.name, self.command)
        self._run_mutation_stats(self.name, self.command)


def register_ms_tests_1(runner):
    def register(name, command):
        runner.register(MsTest1(name=name, group="ms", command=command))

    register("size-change1", "10 10000 -t 2.0 -eN 0.1 2.0")
    register("growth-rate-change1", "10 10000 -t 2.0 -eG 0.1 5.0")
    register("growth-rate-2-pops1", "10 10000 -t 2.0 -I 2 5 5 2.5 -G 5.0")
    register("growth-rate-2-pops2", "10 10000 -t 2.0 -I 2 5 5 2.5 -G 5.0 -g 1 0.1")
    register("growth-rate-2-pops3", "10 10000 -t 2.0 -I 2 5 5 2.5 -g 1 0.1")
    register("growth-rate-2-pops4", "10 10000 -t 2.0 -I 2 5 5 2.5 -eg 1.0 1 5.0")
    register("pop-size-2-pops1", "100 10000 -t 2.0 -I 2 50 50 2.5 -n 1 0.1")
    register("pop-size-2-pops2", "100 10000 -t 2.0 -I 2 50 50 2.5 -g 1 2 -n 1 0.1")
    register("pop-size-2-pops3", "100 10000 -t 2.0 -I 2 50 50 2.5 -eN 0.5 3.5")
    register("pop-size-2-pops4", "100 10000 -t 2.0 -I 2 50 50 2.5 -en 0.5 1 3.5")
    register("migration-rate-2-pops1", "100 10000 -t 2.0 -I 2 50 50 0 -eM 3 5")
    register("migration-matrix-2-pops1", "100 10000 -t 2.0 -I 2 50 50 -ma x 10 0 x")
    register(
        "migration-matrix-2-pops2", "100 10000 -t 2.0 -I 2 50 50 -m 1 2 10 -m 2 1 50"
    )
    register("migration-rate-change-2-pops1", "100 10000 -t 2.0 -I 2 50 50 -eM 5 10")
    register(
        "migration-matrix-entry-change-2-pops1",
        "100 10000 -t 2.0 -I 2 50 50 -em 0.5 2 1 10",
    )
    register(
        "migration-matrix-change-2-pops1",
        "100 10000 -t 2.0 -I 2 50 50 -ema 10.0 2 x 10 0 x",
    )
    cmd = """100 10000 -t 2.0 -I 2 50 50 -ema 1.0
      2 x 0.1 0 x -eN 1.1 0.001 -ema 10 2 x 0 10 x"""
    register(
        "migration-matrix-change-2-pops2", cmd,
    )
    register("population-split-2-pops1", "100 10000 -t 2.0 -I 2 50 50 5.0 -ej 2.0 1 2")
    register(
        "population-split-4-pops1", "100 10000 -t 2.0 -I 4 50 50 0 0 2.0 -ej 0.5 2 1"
    )
    register(
        "population-split-4-pops2",
        "100 10000 -t 2.0 -I 4 25 25 25 25 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1",
    )
    register(
        "population-split-4-pops3",
        "100 10000 -t 2.0 -I 4 25 25 25 25 -ej 1 2 1 -em 1.5 4 1 2 -ej 2 3 1 -ej 3 4 1",
    )
    register("admixture-1-pop1", "1000 1000 -t 2.0 -es 0.1 1 0.5 -em 0.1 1 2 1")
    register("admixture-1-pop2", "1000 1000 -t 2.0 -es 0.1 1 0.1 -em 0.1 1 2 1")
    register("admixture-1-pop3", "1000 1000 -t 2.0 -es 0.01 1 0.1 -em 0.1 2 1 1")
    register(
        "admixture-1-pop4", "1000 1000 -t 2.0 -es 0.01 1 0.1 -es 0.1 2 0 -em 0.1 3 1 1"
    )
    register("admixture-1-pop5", "1000 1000 -t 2.0 -es 0.01 1 0.1 -ej 1 2 1")
    register("admixture-1-pop6", "1000 1000 -t 2.0 -es 0.01 1 0.0 -eg 0.02 2 5.0 ")
    register("admixture-1-pop7", "1000 1000 -t 2.0 -es 0.01 1 0.0 -en 0.02 2 5.0 ")
    register(
        "admixture-2-pop1", "1000 1000 -t 2.0 -I 2 500 500 1 -es 0.01 1 0.1 -ej 1 3 1"
    )
    register(
        "admixture-2-pop2",
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -em 2.0 3 1 1",
    )
    register(
        "admixture-2-pop3",
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -G 5.0 " "-em 2.0 3 1 1",
    )
    register(
        "admixture-2-pop4",
        "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -eg 0.02 1 5.0 -em 0.02 3 1 1",
    )
    # FIXME disabling this test until GC with recombination rate=0
    # register("gene-conversion-1-r0", "100 10000 -t 5.0 -r 0 2501 -c 10 1")
    register("gene-conversion-1", "100 10000 -t 5.0 -r 0.01 2501 -c 1000 1")
    register("gene-conversion-2", "100 10000 -t 5.0 -r 10 2501 -c 2 1")
    register("gene-conversion-2-tl-10", "100 10000 -t 5.0 -r 10 2501 -c 2 10")
    register("gene-conversion-2-tl-100", "100 10000 -t 5.0 -r 10 2501 -c 2 100")
    register("msdoc-simple-ex", "4 20000 -t 5.0")
    register("msdoc-recomb-ex", "15 1000 -t 10.04 -r 100.0 2501")
    register("msdoc-structure-ex1", "15 1000 -t 2.0 -I 3 10 4 1 5.0")
    register(
        "msdoc-structure-ex2", "15 1000 -t 2.0 -I 3 10 4 1 5.0 -m 1 2 10.0 -m 2 1 9.0"
    )
    register(
        "msdoc-structure-ex3",
        "15 1000 -t 10.0 -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x",
    )
    register("msdoc-outgroup-sequence", "11 1000 -t 2.0 -I 2 1 10 -ej 6.0 1 2")
    cmd = "15 10000 -t 11.2 -I 2 3 12 -g 1 44.36 -n 2 \
     0.125 -eg 0.03125 1 0.0 -en 0.0625 2 0.05 -ej 0.09375 2 1"
    register(
        "msdoc-two-species", cmd,
    )
    cmd = "15 10000 -t 3.0 -I 6 0 7 0 0 8 0 -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 -m 3 \
     2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 -em 2.0 3 4 2.5 -em 2.0 4 3 2.5"
    register(
        "msdoc-stepping-stone", cmd,
    )
    register("simultaneous-ex1", "10 10000 -t 2.0 -eN 0.3 0.5 -eG .3 7.0")
    register("zero-growth-rate", "10 10000 -t 2.0 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5")
    cmd = "4 1000 -t 2508 -I 2 2 2 0 -n 2 2.59 \
     -ma x 0 1.502 x -ej 0.9485 1 2 -r 23.76 3000"
    register(
        "konrad-1", cmd,
    )
    cmd = "3 10000 -t 0.423 -I 3 1 1 1 -es 0.0786 1 0.946635 \
      -ej 0.0786 4 3 -ej 0.189256 1 2 -ej 0.483492 2 3"
    register(
        "konrad-2", cmd,
    )
    register("konrad-3", "100 100 -t 2 -I 10 10 10 10 10 10 10 10 10 10 10 0.001 ")


@attr.s
class MsTest2(MsTest):
    command = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        logging.debug(f"{self.name} {self.command}")
        self._run_variable_recombination_coalescent_stats(self.name, self.command)


def register_ms_tests_2(runner):
    def register(name, command):
        runner.register(MsTest2(name=name, group="ms", command=command))

    register(
        "mshotdoc-hotspot-ex",
        "10 1000 -t 10.4 -r 10.0 25000 -v 2 100 200 10 7000 8000 20",
    )
    register(
        "mshot-zero-recomb-interval", "10 1000 -t 10.4 -r 10.0 25000 -v 1 5000 13000 0"
    )
    register("mshot-zero-recomb", "10 1000 -t 10.4 -r 10.0 25000 -v 1 100 25000 0")
    hotspots = "4 1000 2000 0 7000 8000 20 12000 15000 10 20000 22000 0"
    register(
        "mshot-high-recomb-variance", f"10 1000 -t 10.4 -r 10.0 25000 -v {hotspots}"
    )


@attr.s
class MsTest3(MsTest):
    command = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        logging.debug(f"{self.name} {self.command}")
        self._run_ms_mshot_stats(self.name, self.command)


def register_ms_tests_3(runner):
    def register(name, command):
        runner.register(MsTest3(name=name, group="ms", command=command))

    register("ms-mshot-consistency-check", "10 1000 -t 10.4 -r 10.0 25000")


def register_ms_tests(runner):
    register_ms_tests_1(runner)
    register_ms_tests_2(runner)
    register_ms_tests_3(runner)


@attr.s
class DiscoalTest(Test):
    _discoal_executable = attr.ib(init=False, default=["./data/discoal"])

    def get_discoal_seeds(self):
        max_seed = 2 ** 16
        seeds = [random.randint(1, max_seed) for j in range(3)]
        return ["-d"] + list(map(str, seeds))

    def _discoal_str_to_ms(self, args):
        # convert discoal string to msprime string
        tokens = args.split(" ")
        # cut out sites param
        del tokens[2]
        # adjust popIDs
        for i in range(len(tokens)):
            # pop size change case
            if tokens[i] == "-en":
                tokens[i + 2] = str(int(tokens[i + 2]) + 1)
            # migration rate case
            if tokens[i] == "-m":
                tokens[i + 1] = str(int(tokens[i + 1]) + 1)
                tokens[i + 2] = str(int(tokens[i + 2]) + 1)
        msp_str = " ".join(tokens)
        return msp_str

    def _run_discoal_mutation_stats(self, args):
        return self._run_sample_stats(
            self._discoal_executable + args.split() + self.get_discoal_seeds()
        )

    def _run_mutation_discoal_stats(self, key, args):
        msp_str = self._discoal_str_to_ms(args)
        df_msp = self._run_msprime_mutation_stats(msp_str)
        df_d = self._run_discoal_mutation_stats(args)
        self._plot_stats(key, "mutation", df_d, df_msp, "discoal", "msp")

    def run_sweep_comparison(self, key, args):
        """
        does a comparison of sweep model vs discoal
        """
        # This is broken currently: https://github.com/tskit-dev/msprime/issues/942
        # Skip noisily
        logging.warning("Skipping sweep comparison due to known bug.")
        return

        # TODO We should be parsing the args string here to derive these values
        # from the input. There's no point in having it as a parameter otherwise.
        df = pd.DataFrame()
        refsize = 1e4
        seqlen = 1e4
        nreps = 1000
        mod = msprime.SweepGenicSelection(
            position=np.floor(seqlen / 2),
            start_frequency=1.0 / (2 * refsize),
            end_frequency=1.0 - (1.0 / (2 * refsize)),
            alpha=1000.0,
            dt=1.0 / (40 * refsize),
        )
        mu = 2.5e-4
        rho = 2.5e-4
        data = collections.defaultdict(list)
        replicates = msprime.simulate(
            10,
            model=mod,
            length=seqlen,
            recombination_rate=rho,
            mutation_rate=mu,
            num_replicates=nreps,
            # Change to Hudson after sweep finishes
            demographic_events=[msprime.SimulationModelChange()],
        )
        for ts in replicates:
            data["pi"].append(ts.diversity(span_normalise=False))
            data["D"].append(ts.Tajimas_D())
            data["ss"].append(ts.segregating_sites(span_normalise=False))
        data["pi"] = np.array(data["pi"]).flatten()
        data["D"] = np.array(data["D"]).flatten()
        data["ss"] = np.array(data["ss"]).flatten()
        df = pd.DataFrame.from_dict(data)
        df = df.fillna(0)
        df_d = self._run_discoal_mutation_stats(args)
        df_df = df_d[["pi", "D", "ss"]]
        logging.debug(f"msp pi mean: {df['pi'].mean()}")
        logging.debug(f"discoal pi mean: {df_df['pi'].mean()}")
        logging.debug(f"msp ss mean: {df['ss'].mean()}")
        logging.debug(f"discoal ss mean: {df_df['ss'].mean()}")
        logging.debug(f"msp D mean: {df['D'].mean()}")
        logging.debug(f"discoal D mean: {df_df['D'].mean()}")
        logging.debug(f"sample sizes msp: {len(df['pi'])} discoal: {len(df_df['pi'])}")
        self._plot_stats(key, f"mutation{df}, {df_df}")


@attr.s
class DiscoalTest1(DiscoalTest):
    command = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        logging.debug(f"{self.name} {self.command}")
        self._run_mutation_discoal_stats(self.name, self.command)


def register_discoal_tests_1(runner):
    def register(name, command):
        runner.register(DiscoalTest1(name=name, group="discoal", command=command))

    register("discoal-simple-ex", "15 1000 100 -t 5.0")
    register("discoal-size-change1", "10 10000 100 -t 10.0 -en 0.1 0 2.0")
    register("discoal-size-change2", "10 10000 100 -t 10.0 -en 0.1 0 0.1")
    register("discoal-size-change3", "10 10000 100 -t 10.0 -en 0.01 0 0.01")
    register(
        "discoal-size-change4", "10 10000 100 -t 10.0 -en 0.01 0 0.5 -en 0.05 0 1.0"
    )


@attr.s
class DiscoalTest2(DiscoalTest):
    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_sweep_comparison(
            "sweep",
            "10 1000 10000 -t 10.0 -r 10.0 -ws 0 -a \
                                  500 -x 0.5 -N 10000",
        )


def register_discoal_tests_2(runner):
    runner.register(DiscoalTest2(name="sweep_vs_discoal", group="discoal"))


def register_discoal_tests(runner):
    register_discoal_tests_1(runner)
    register_discoal_tests_2(runner)


@attr.s
class DtwfVsCoalescentTest(Test):
    _discoal_executable = attr.ib(init=False, default=["./data/discoal"])

    def run_dtwf_pedigree_comparison(self, test_name, **kwargs):
        df = pd.DataFrame()
        pedigree = kwargs["pedigree"]
        assert kwargs["sample_size"] % 2 == 0
        sample_size = kwargs["sample_size"]
        sample_size_diploid = sample_size // 2
        for model in ["wf_ped", "dtwf"]:
            kwargs["model"] = model
            kwargs["pedigree"] = None
            kwargs["sample_size"] = sample_size
            if model == "wf_ped":
                kwargs["sample_size"] = sample_size_diploid
                kwargs["pedigree"] = pedigree

                des = []
                if "demographic_events" in kwargs:
                    des = kwargs["demographic_events"]
                max_ped_time = max(pedigree.times)
                des.append(msprime.SimulationModelChange(max_ped_time, "dtwf"))
                des = sorted(des, key=lambda x: x.time)
                kwargs["demographic_events"] = des

            logging.debug(f"Running: {kwargs}")
            data = collections.defaultdict(list)
            try:
                replicates = msprime.simulate(**kwargs)
                for ts in replicates:
                    t_mrca = np.zeros(ts.num_trees)
                    for tree in ts.trees():
                        t_mrca[tree.index] = tree.time(tree.root)
                    data["tmrca_mean"].append(np.mean(t_mrca))
                    data["num_trees"].append(ts.num_trees)
                    data["model"].append(model)
            except Exception as e:
                logging.warning(f"TEST FAILED!!!: {e}")
                return
            df = df.append(pd.DataFrame(data))

        basedir = make_test_dir(test_name)

        df_wf_ped = df[df.model == "wf_ped"]
        df_dtwf = df[df.model == "dtwf"]
        for stat in ["tmrca_mean", "num_trees"]:
            v1 = df_wf_ped[stat]
            v2 = df_dtwf[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = os.path.join(basedir, f"{stat}.png")
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def add_dtwf_vs_pedigree_single_locus(self):
        """
        Checks the DTWF against the standard coalescent at a single locus.
        """
        pedigree_file = "tests/data/pedigrees/wf_100Ne_10000gens.txt"
        pedigree = msprime.Pedigree.read_txt(pedigree_file, time_col=3)

        def f():
            self.run_dtwf_pedigree_comparison(
                "dtwf_vs_pedigree_single_locus",
                sample_size=10,
                Ne=100,
                num_replicates=400,
                length=1,
                pedigree=pedigree,
                recombination_rate=0,
                mutation_rate=1e-8,
            )

        self._instances["dtwf_vs_pedigree_single_locus"] = f

    def add_dtwf_vs_pedigree_short_region(self):
        """
        Checks the DTWF against the standard coalescent at a single locus.
        """
        pedigree_file = "tests/data/pedigrees/wf_100Ne_10000gens.txt"
        pedigree = msprime.Pedigree.read_txt(pedigree_file, time_col=3)

        def f():
            self.run_dtwf_pedigree_comparison(
                "dtwf_vs_pedigree_short_region",
                sample_size=10,
                Ne=100,
                num_replicates=400,
                length=1e6,
                pedigree=pedigree,
                recombination_rate=1e-8,
                mutation_rate=1e-8,
            )

        self._instances["dtwf_vs_pedigree_short_region"] = f

    def add_dtwf_vs_pedigree_long_region(self):
        """
        Checks the DTWF against the standard coalescent at a single locus.
        """
        pedigree_file = "tests/data/pedigrees/wf_100Ne_10000gens.txt"
        pedigree = msprime.Pedigree.read_txt(pedigree_file, time_col=3)

        def f():
            self.run_dtwf_pedigree_comparison(
                "dtwf_vs_pedigree_long_region",
                sample_size=10,
                Ne=100,
                num_replicates=200,
                length=1e8,
                pedigree=pedigree,
                recombination_rate=1e-8,
                mutation_rate=1e-8,
            )

        self._instances["dtwf_vs_pedigree_long_region"] = f

    def run_dtwf_coalescent_comparison(self, test_name, **kwargs):
        basedir = make_test_dir(test_name)
        df = run_dtwf_coalescent_stats(**kwargs)
        plot_dtwf_coalescent_stats(basedir, df)

    def run_dtwf_coalescent_tree_interval_comparison(self, test_name, **kwargs):
        basedir = make_test_dir(test_name)
        df = run_dtwf_coalescent_stats(**kwargs)
        plot_tree_intervals(basedir, df)

    def run_dtwf_slim_comparison(self, test_name, slim_args, **kwargs):

        df = pd.DataFrame()

        kwargs["model"] = "dtwf"
        logging.debug(f"Running: {kwargs}")
        replicates = msprime.simulate(**kwargs)
        data = collections.defaultdict(list)
        for ts in replicates:
            t_mrca = np.zeros(ts.num_trees)
            for tree in ts.trees():
                t_mrca[tree.index] = tree.time(tree.root)
            data["tmrca_mean"].append(np.mean(t_mrca))
            data["num_trees"].append(ts.num_trees)
            data["model"].append("dtwf")

        basedir = make_test_dir(test_name)

        slim_script = os.path.join(basedir, "slim_script.txt")
        outfile = os.path.join(basedir, "slim.trees")
        slim_args["OUTFILE"] = outfile
        write_slim_script(slim_script, slim_args)

        cmd = self._slim_executable + [slim_script]
        for _ in tqdm.tqdm(range(kwargs["num_replicates"])):
            subprocess.check_output(cmd)
            ts = msprime.load(outfile)
            ts = subsample_simplify_slim_treesequence(ts, slim_args["sample_sizes"])

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
            f = os.path.join(basedir, f"{stat}.png")
            pyplot.xlabel("DTWF")
            pyplot.ylabel("SLiM")
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def check_slim_version(self):
        # This may not be robust but it's a start
        min_version = 3.1
        raw_str = subprocess.check_output(self._slim_executable + ["-version"])
        version_list = str.split(str(raw_str))
        for i in range(len(version_list)):
            if version_list[i].lower() == "version":
                version_str = version_list[i + 1]
                break
        version = float(version_str.strip(" ,")[0:3])
        assert version >= min_version, "Require SLiM >= 3.1!"


@attr.s
class DtwfvsCoalescentTest1(DtwfVsCoalescentTest):
    simulate_args = attr.ib(factory=dict)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        logging.debug(f"{self.name} {self.simulate_args}")
        if "recombination_map" in self.simulate_args.keys():
            if "uniform" in self.simulate_args.keys():
                self.simulate_args[
                    "recombination_map"
                ] = msprime.RecombinationMap.uniform_map(
                    **self.simulate_args["recombination_map"]
                )
                self.simulate_args.pop("uniform", None)
            else:
                self.simulate_args["recombination_map"] = msprime.RecombinationMap(
                    **self.simulate_args["recombination_map"]
                )
        self.run_dtwf_coalescent_comparison(self.name, **self.simulate_args)


def dtwf_vs_coalescent_single_locus(group):
    """
    Checks the DTWF against the standard coalescent at a single locus.
    """

    return DtwfvsCoalescentTest1(
        "dtwf_vs_coalescent_single_locus",
        group,
        simulate_args=dict(sample_size=10, Ne=1000, num_replicates=300,),
    )


def dtwf_vs_coalescent_recomb_discrete_hotspots(group):
    """
    Checks the DTWF against the standard coalescent with a
    discrete recombination map with variable rates.
    """
    test_name = "dtwf_vs_coalescent_discrete_hotspots"

    recombination_map = {
        "positions": [0, 100, 500, 900, 1200, 1500, 2000],
        "rates": [0.00001, 0, 0.0002, 0.00005, 0, 0.001, 0],
        "discrete": True,
    }
    return DtwfvsCoalescentTest1(
        test_name,
        group,
        simulate_args=dict(
            sample_size=10,
            Ne=1000,
            recombination_map=recombination_map,
            num_replicates=300,
        ),
    )


def dtwf_vs_coalescent_recomb_continuous_hotspots(group):
    """
    Checks the DTWF against the standard coalescent with a
    continuous recombination map with variable rates.
    """
    test_name = "dtwf_vs_coalescent_continuous_hotspots"

    recombination_map = {
        "positions": [0, 0.1, 0.5, 0.9, 1.2, 1.5, 2.0],
        "rates": [0.00001, 0, 0.0002, 0.00005, 0, 0.001, 0],
    }

    return DtwfvsCoalescentTest1(
        test_name,
        group,
        simulate_args=dict(
            sample_size=10,
            Ne=1000,
            recombination_map=recombination_map,
            num_replicates=300,
        ),
    )


def dtwf_vs_coalescent_single_forced_recombination(group):
    test_name = "dtwf_vs_coalescent_single_forced_recombination"

    recombination_map = {
        "positions": [0, 100, 101, 201],
        "rates": [0, 1, 0, 0],
        "discrete": True,
    }

    return DtwfvsCoalescentTest1(
        test_name,
        group,
        simulate_args=dict(
            sample_size=10,
            Ne=10,
            num_replicates=1,
            recombination_map=recombination_map,
        ),
    )


def dtwf_vs_coalescent_low_recombination(group):
    """
    Checks the DTWF against the standard coalescent at a single locus.
    """

    return DtwfvsCoalescentTest1(
        "dtwf_vs_coalescent_low_recombination",
        group,
        simulate_args=dict(
            sample_size=10, Ne=1000, num_replicates=400, recombination_rate=0.01,
        ),
    )


def dtwf_vs_coalescent_2_pops_massmigration(group):
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=10, initial_size=1000),
        msprime.PopulationConfiguration(sample_size=10, initial_size=1000),
    ]
    recombination_map = {"positions": [0, int(1e6)], "rates": [1e-8, 0]}
    demographic_events = [
        msprime.MassMigration(time=300, source=1, destination=0, proportion=1.0)
    ]

    test_name = "dtwf_vs_coalescent_2_pops_massmigration"

    return DtwfvsCoalescentTest1(
        test_name,
        group,
        simulate_args=dict(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            # Ne=0.5,
            num_replicates=300,
            recombination_map=recombination_map,
        ),
    )


def dtwf_vs_coalescent_2_pop_growth(group):
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=10, initial_size=1000, growth_rate=0.01
        )
    ]
    recombination_map = {
        "positions": [0, int(5e7)],
        "rates": [1e-8, 0],
        "discrete": True,
    }

    return DtwfvsCoalescentTest1(
        "dtwf_vs_coalescent_2_pop_growth",
        group,
        simulate_args=dict(
            population_configurations=population_configurations,
            recombination_map=recombination_map,
            num_replicates=300,
        ),
    )


def dtwf_vs_coalescent_2_pop_shrink(group):
    initial_size = 1000

    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=10, initial_size=initial_size, growth_rate=-0.01
        )
    ]
    recombination_map = {
        "positions": [0, int(1e7)],
        "rates": [1e-8, 0],
        "discrete": True,
    }
    demographic_events = [
        msprime.PopulationParametersChange(
            time=200, initial_size=initial_size, growth_rate=0.01, population_id=0
        )
    ]

    return DtwfvsCoalescentTest1(
        "dtwf_vs_coalescent_2_pop_shrink",
        group,
        simulate_args=dict(
            population_configurations=population_configurations,
            recombination_map=recombination_map,
            demographic_events=demographic_events,
            num_replicates=300,
        ),
    )


def dtwf_vs_coalescent_multiple_bottleneck(group):
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=5, initial_size=1000),
        msprime.PopulationConfiguration(sample_size=5, initial_size=1000),
    ]
    recombination_map = {"positions": [0, int(1e6)], "rates": [1e-8, 0]}
    # migration_matrix = [[0, 0.1], [0.1, 0]]

    demographic_events = [
        msprime.PopulationParametersChange(
            time=100, initial_size=100, growth_rate=-0.01, population_id=0
        ),
        msprime.PopulationParametersChange(
            time=200, initial_size=100, growth_rate=-0.01, population_id=1
        ),
        msprime.PopulationParametersChange(
            time=300, initial_size=1000, growth_rate=0.01, population_id=0
        ),
        msprime.PopulationParametersChange(
            time=400, initial_size=1000, growth_rate=0.01, population_id=1
        ),
        msprime.PopulationParametersChange(
            time=500, initial_size=100, growth_rate=0, population_id=0
        ),
        msprime.PopulationParametersChange(
            time=600, initial_size=100, growth_rate=0, population_id=1
        ),
        msprime.MigrationRateChange(time=700, rate=0.1, matrix_index=(0, 1)),
    ]

    return DtwfvsCoalescentTest1(
        "dtwf_vs_coalescent_multiple_bottleneck",
        group,
        simulate_args=dict(
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            # migration_matrix=migration_matrix,
            num_replicates=400,
            recombination_map=recombination_map,
        ),
    )


def register_dtwfvscoalescent_tests_1(runner):
    group = "dtwfvscoalescent"
    tests = [
        dtwf_vs_coalescent_single_locus(group),
        dtwf_vs_coalescent_recomb_discrete_hotspots(group),
        dtwf_vs_coalescent_recomb_continuous_hotspots(group),
        dtwf_vs_coalescent_single_forced_recombination(group),
        dtwf_vs_coalescent_low_recombination(group),
        dtwf_vs_coalescent_2_pops_massmigration(group),
        dtwf_vs_coalescent_2_pop_growth(group),
        dtwf_vs_coalescent_2_pop_shrink(group),
        dtwf_vs_coalescent_multiple_bottleneck(group),
    ]

    for test in tests:
        runner.register(test)


def add_dtwf_vs_coalescent(
    key,
    group,
    initial_sizes,
    sample_sizes,
    num_loci,
    recombination_rate,
    migration_matrix=None,
    growth_rates=None,
    num_replicates=None,
):
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
            de = msprime.PopulationParametersChange(t_100, growth_rate=0, population=i)
            demographic_events.append(de)

            growth_rate = growth_rates[i]
        else:
            # Enforce zero growth rate for small populations
            logging.warning(
                f"Warning - setting growth rate to zero for small \
                population of size {initial_sizes[i]}",
            )
            growth_rate = 0

        population_configurations.append(
            msprime.PopulationConfiguration(
                sample_size=sample_sizes[i],
                initial_size=initial_sizes[i],
                growth_rate=growth_rate,
            )
        )

    recombination_map = {
        "length": num_loci,
        "rate": recombination_rate,
        "discrete": True,
    }

    if migration_matrix is None:
        default_mig_rate = 0.05
        migration_matrix = []
        for i in range(num_pops):
            row = [default_mig_rate] * num_pops
            row[i] = 0
            migration_matrix.append(row)

    return DtwfvsCoalescentTest1(
        key,
        group,
        simulate_args=dict(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            num_replicates=num_replicates,
            demographic_events=demographic_events,
            recombination_map=recombination_map,
            uniform=True,
        ),
    )


def register_dtwfvscoalescent_tests_2(runner):
    group = "dtwfvscoalescent"
    tests = [
        add_dtwf_vs_coalescent(
            "dtwf_vs_coalescent_long_region", group, [1000], [10], int(1e8), 1e-8
        ),
        add_dtwf_vs_coalescent(
            "dtwf_vs_coalescent_short_region", group, [1000], [10], int(1e6), 1e-8
        ),
        add_dtwf_vs_coalescent(
            "dtwf_vs_coalescent_2_pops",
            group,
            [500, 500],
            [5, 5],
            int(1e6),
            1e-8,
            num_replicates=500,
        ),
        add_dtwf_vs_coalescent(
            "dtwf_vs_coalescent_3_pops",
            group,
            [500, 500, 500],
            [5, 2, 0],
            int(1e7),
            1e-8,
        ),
        add_dtwf_vs_coalescent(
            "dtwf_vs_coalescent_4_pops",
            group,
            [1000, 1000, 1000, 1000],
            [0, 20, 0, 0],
            int(1e6),
            1e-8,
            num_replicates=500,
        ),
    ]
    migration_matrix = [[0, 0.2, 0.1], [0.1, 0, 0.2], [0.2, 0.1, 0]]
    tests.append(
        add_dtwf_vs_coalescent(
            "dtwf_vs_coalescent_3_pops_asymm_mig",
            group,
            [500, 500, 500],
            [20, 0, 0],
            int(1e6),
            1e-8,
            migration_matrix=migration_matrix,
            num_replicates=500,
        )
    )

    migration_matrix = [[0, 0.5], [0.7, 0]]
    tests.append(
        add_dtwf_vs_coalescent(
            "dtwf_vs_coalescent_2_pops_high_asymm_mig",
            group,
            [1000, 1000],
            [10, 10],
            int(1e6),
            1e-8,
            migration_matrix=migration_matrix,
            num_replicates=200,
            growth_rates=[0.005, 0.005],
        )
    )

    for test in tests:
        runner.register(test)


@attr.s
class DtwfvsCoalescentTest2(DtwfVsCoalescentTest):
    slim_args = attr.ib(factory=dict)
    simulate_args = attr.ib(factory=dict)

    def run(self, output_dir):
        self.check_slim_version()
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        logging.debug(f"{self.name} {self.slim_args} {self.simulate_args}")
        if "recombination_map" in self.simulate_args.keys():
            self.simulate_args["recombination_map"] = msprime.RecombinationMap(
                **self.simulate_args["recombination_map"]
            )

        self.run_dtwf_slim_comparison(self.name, self.slim_args, **self.simulate_args)


def add_dtwf_vs_slim(
    key,
    group,
    initial_sizes,
    sample_sizes,
    num_loci,
    recombination_rate,
    migration_matrix=None,
    num_replicates=None,
):
    """
    Generic test of DTWF vs SLiM WF simulator, without growth rates
    """
    assert len(sample_sizes) == len(initial_sizes)

    num_pops = len(sample_sizes)
    slim_args = {}

    if num_replicates is None:
        num_replicates = 200

    slim_args["sample_sizes"] = sample_sizes

    population_configurations = []
    slim_args["POP_STRS"] = ""
    for i in range(len(sample_sizes)):
        population_configurations.append(
            msprime.PopulationConfiguration(
                sample_size=sample_sizes[i],
                initial_size=initial_sizes[i],
                growth_rate=0,
            )
        )
        slim_args["POP_STRS"] += "sim.addSubpop('p{i}', {N});\n".format(
            i=i, N=initial_sizes[i]
        )

    if migration_matrix is None:
        default_mig_rate = 0.01
        migration_matrix = []
        for i in range(num_pops):
            row = [default_mig_rate] * num_pops
            row[i] = 0
            migration_matrix.append(row)

    # SLiM rates are 'immigration' forwards in time, which matches
    # DTWF backwards-time 'emmigration'
    assert len(migration_matrix) == num_pops
    if num_pops > 1:
        for i in range(num_pops):
            row = migration_matrix[i]
            indices = [j for j in range(num_pops) if j != i]
            pop_names = ["p" + str(j) for j in indices]
            rates = [str(row[j]) for j in indices]

            to_pop_str = ",".join(pop_names)
            rate_str = ",".join(rates)

            mig_str = "p{}.setMigrationRates(c({}), c({}));\n".format(
                i, to_pop_str, rate_str
            )
            slim_args["POP_STRS"] += mig_str

    num_loci = int(num_loci)
    recombination_map = {
        "positions": [0, num_loci],
        "rates": [recombination_rate, 0],
        "discrete": True,
    }
    slim_args["RHO"] = recombination_rate
    slim_args["NUM_LOCI"] = num_loci

    return DtwfvsCoalescentTest2(
        key,
        group,
        slim_args,
        simulate_args=dict(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            num_replicates=num_replicates,
            recombination_map=recombination_map,
        ),
    )


def register_dtwfvscoalescent_tests_3(runner):
    group = "dtwfvscoalescent"
    tests = [
        add_dtwf_vs_slim("dtwf_vs_slim_single_locus", group, [10], [10], 1, 0),
        add_dtwf_vs_slim(
            "dtwf_vs_slim_short_region",
            group,
            [100],
            [10],
            1e7,
            1e-8,
            num_replicates=200,
        ),
        add_dtwf_vs_slim(
            "dtwf_vs_slim_long_region", group, [50], [10], 1e8, 1e-8, num_replicates=200
        ),
    ]

    for test in tests:
        runner.register(test)


# run only if args.extended is true
def add_dtwf_vs_coalescent_random_instance(
    key, group, num_populations=1, num_replicates=200, num_demographic_events=0
):

    N = num_populations
    num_loci = np.random.randint(1e5, 1e7)
    rho = 1e-8
    recombination_map = {
        "positions": [0, num_loci],
        "rates": [rho, 0],
        "discrete": True,
    }

    population_configurations = []
    for _ in range(N):
        population_configurations.append(
            msprime.PopulationConfiguration(
                sample_size=np.random.randint(1, 10), initial_size=int(1000 / N)
            )
        )

    migration_matrix = []
    for i in range(N):
        migration_matrix.append(
            [random.uniform(0.05, 0.25) * (j != i) for j in range(N)]
        )

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
                time=t,
                initial_size=initial_size,
                growth_rate=growth_rate,
                population_id=pop_id,
            )
        )

        if random.random() < 0.5 and N >= 2:
            rate = random.uniform(0.05, 0.25)
            index = tuple(
                np.random.choice(range(num_populations), size=2, replace=False)
            )
            demographic_events.append(
                msprime.MigrationRateChange(time=t, rate=rate, matrix_index=index)
            )

    # Collect all pops together to control coalescence times for DTWF
    for i in range(1, N):
        demographic_events.append(
            msprime.MassMigration(time=t_max, source=i, destination=0, proportion=1.0)
        )

    demographic_events.append(
        msprime.PopulationParametersChange(
            time=t_max, initial_size=100, growth_rate=0, population_id=0
        )
    )

    return DtwfvsCoalescentTest1(
        key,
        group,
        simulate_args=dict(
            migration_matrix=migration_matrix,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            num_replicates=num_replicates,
            recombination_map=recombination_map,
        ),
    )


def register_dtwfvscoalescent_tests_4(runner):
    group = "dtwfvscoalescent"
    tests = [
        add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_1",
            group,
            num_populations=2,
            num_replicates=200,
            num_demographic_events=3,
        ),
        add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_2",
            group,
            num_populations=3,
            num_replicates=200,
            num_demographic_events=3,
        ),
        add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_3",
            group,
            num_populations=2,
            num_replicates=200,
            num_demographic_events=6,
        ),
        add_dtwf_vs_coalescent_random_instance(
            "dtwf_vs_coalescent_random_4",
            group,
            num_populations=1,
            num_replicates=200,
            num_demographic_events=8,
        ),
    ]

    for test in tests:
        runner.register(test)


def register_dtwfvscoalescent_tests(runner, add_random_dtwf):
    if add_random_dtwf:
        register_dtwfvscoalescent_tests_4(runner)
    else:
        register_dtwfvscoalescent_tests_1(runner)
        register_dtwfvscoalescent_tests_2(runner)
        register_dtwfvscoalescent_tests_3(runner)


@attr.s
class XiHudsonTest(Test):
    def verify_breakpoint_distribution(
        self, basedir_name, name, sample_size, Ne, r, L, model, growth_rate=0
    ):
        """
        Verifies that the number of recombination breakpoints is proportional to
        the total branch length across all trees.
        """
        basedir = make_test_dir(basedir_name)
        ts = msprime.simulate(
            Ne=Ne,
            recombination_rate=r,
            length=L,
            population_configurations=[
                msprime.PopulationConfiguration(
                    sample_size=sample_size, initial_size=Ne, growth_rate=growth_rate
                )
            ],
            model=model,
        )
        empirical = []
        for tree in ts.trees():
            area = tree.total_branch_length * tree.span
            empirical.append(area)

        scipy.stats.probplot(empirical, dist=scipy.stats.expon(Ne * r), plot=pyplot)
        path = os.path.join(basedir, f"{name}_growth={growth_rate}.png")
        logging.debug(f"Writing {path}")
        pyplot.savefig(path)
        pyplot.close("all")

    # used in Xi and Hudson tests
    def verify_recombination(
        self, basedir_name, name, sample_size, Ne, r, m, L, model, growth_rate=0
    ):
        """
        Verifies that the number of recombination equals the number of mutation.
        """
        basedir = make_test_dir(basedir_name)
        empirical_theta = []
        empirical_rho = []
        for _ in range(1, 500):
            ts = msprime.simulate(
                Ne=Ne,
                recombination_rate=r,
                mutation_rate=m,
                length=L,
                population_configurations=[
                    msprime.PopulationConfiguration(
                        sample_size=sample_size,
                        initial_size=Ne,
                        growth_rate=growth_rate,
                    )
                ],
                model=model,
            )
            empirical_theta.append(ts.get_num_sites())
            ts = msprime.simulator_factory(
                Ne=Ne,
                recombination_rate=r,
                length=L,
                population_configurations=[
                    msprime.PopulationConfiguration(
                        sample_size=sample_size,
                        initial_size=Ne,
                        growth_rate=growth_rate,
                    )
                ],
                model=model,
            )
            ts.run()
            empirical_rho.append(ts.num_breakpoints)
        empirical_rho.sort()
        empirical_theta.sort()
        empirical_rho = np.array(empirical_rho)
        empirical_theta = np.array(empirical_theta)
        plot_qq(empirical_theta, empirical_rho)
        path = os.path.join(basedir, f"{name}_growth={growth_rate}_rec_check.png")
        logging.debug(f"Writing {path}")
        pyplot.savefig(path)
        pyplot.close("all")


@attr.s
class XiTest(XiHudsonTest):
    def run_xi_hudson_comparison(self, test_name, xi_model, **kwargs):
        df = pd.DataFrame()
        for model in ["hudson", xi_model]:
            kwargs["model"] = model
            model_str = "hudson"
            if model != "hudson":
                model_str = "Xi"
            logging.debug(f"Running: {kwargs}")
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

        basedir = make_test_dir(test_name)

        df_hudson = df[df.model == "hudson"]
        df_xi = df[df.model == "Xi"]
        for stat in ["tmrca_mean", "num_trees", "num_nodes", "num_edges"]:
            v1 = df_hudson[stat]
            v2 = df_xi[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = os.path.join(basedir, f"{stat}.png")
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def compare_xi_dirac_sfs(self, sample_size, psi, c, sfs, num_replicates=1000):
        """
        Runs simulations of the xi dirac model and calculates
        E[Bi/B] (Bi branch length having i leaves and B total branch length)
        and compares to the expected SFS.
        """
        logging.debug(f"running SFS for {sample_size} {psi} {c}")
        reps = msprime.simulate(
            sample_size,
            num_replicates=num_replicates,
            model=msprime.DiracCoalescent(psi=psi, c=c),
        )

        data = collections.defaultdict(list)
        tbl_sum = [0] * (sample_size - 1)
        for ts in reps:
            for tree in ts.trees():
                tot_bl = 0.0
                tbl = [0] * (sample_size - 1)
                for node in tree.nodes():
                    if tree.parent(node) != msprime.NULL_NODE:
                        tbl[tree.num_samples(node) - 1] = tbl[
                            tree.num_samples(node) - 1
                        ] + tree.branch_length(node)
                        tot_bl = tot_bl + tree.branch_length(node)

                for xi in range(sample_size - 1):
                    rescaled_x = tbl[xi] / tot_bl
                    data["total_branch_length"].append(rescaled_x)
                    tbl_sum[xi] = tbl_sum[xi] + rescaled_x
                data["num_leaves"].extend(range(1, sample_size))

        basedir = make_test_dir("xi_dirac_expected_sfs")
        f = os.path.join(basedir, f"n={sample_size}_psi={psi}_c={c}.png")
        ax = sns.violinplot(
            data=data, x="num_leaves", y="total_branch_length", color="grey"
        )
        ax.set_xlabel("num leaves")
        l1 = ax.plot(np.arange(sample_size - 1), sfs[::], "--", linewidth=3)
        l2 = ax.plot(
            np.arange(sample_size - 1),
            [x / num_replicates for x in tbl_sum],
            "--",
            linewidth=3,
        )
        ax.legend((l1[0], l2[0]), ("Expected", "Observed"))
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")

    def compare_normalised_xi_dirac_sfs(
        self, sample_size, psi, c, sfs, num_replicates=1000
    ):
        """
        Runs simulations of the xi dirac model and calculates
        E[Bi]/E[B] (Bi branch length having i leaves and B total branch length)
        and compares to the expected SFS.
        """
        logging.debug(f"running SFS for {sample_size} {psi} {c}")
        reps = msprime.simulate(
            sample_size,
            num_replicates=num_replicates,
            model=msprime.DiracCoalescent(psi=psi, c=c),
        )

        data = collections.defaultdict(list)
        tbl_sum = [0] * (sample_size - 1)
        tot_bl_sum = [0]
        for ts in reps:
            for tree in ts.trees():
                tot_bl = 0.0
                tbl = [0] * (sample_size - 1)
                for node in tree.nodes():
                    if tree.parent(node) != msprime.NULL_NODE:
                        tbl[tree.num_samples(node) - 1] = tbl[
                            tree.num_samples(node) - 1
                        ] + tree.branch_length(node)
                        tot_bl = tot_bl + tree.branch_length(node)

                for xi in range(sample_size - 1):
                    rescaled_x = tbl[xi]
                    data["total_branch_length"].append(rescaled_x / tot_bl)
                    tbl_sum[xi] = tbl_sum[xi] + rescaled_x
                tot_bl_sum[0] = tot_bl_sum[0] + tot_bl
                data["num_leaves"].extend(range(1, sample_size))

        basedir = make_test_dir("xi_dirac_expected_sfs")
        f = os.path.join(basedir, f"n={sample_size}_psi={psi}_c={c}.png")
        ax = sns.violinplot(
            data=data, x="num_leaves", y="total_branch_length", color="grey"
        )
        ax.set_xlabel("num leaves")
        l1 = ax.plot(np.arange(sample_size - 1), sfs[::], "--", linewidth=3)
        l2 = ax.plot(
            np.arange(sample_size - 1),
            [(x / num_replicates) / (tot_bl_sum[0] / num_replicates) for x in tbl_sum],
            "--",
            linewidth=3,
        )
        ax.legend((l1[0], l2[0]), ("Expected", "Observed"))
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")

    def compare_normalized_xi_beta_sfs(
        self, sample_size, alpha, sfs, num_replicates=1000
    ):
        """
        Runs simulations of the xi beta model and compares to the expected SFS.
        """
        logging.debug(f"running SFS for {sample_size} {alpha}")
        reps = msprime.simulate(
            sample_size,
            num_replicates=num_replicates,
            model=msprime.BetaCoalescent(alpha=alpha, truncation_point=1),
        )

        data = collections.defaultdict(list)
        tbl_sum = [0] * (sample_size - 1)
        tot_bl_sum = [0]
        for ts in reps:
            for tree in ts.trees():
                tot_bl = 0.0
                tbl = [0] * (sample_size - 1)
                for node in tree.nodes():
                    if tree.parent(node) != msprime.NULL_NODE:
                        tbl[tree.num_samples(node) - 1] = tbl[
                            tree.num_samples(node) - 1
                        ] + tree.branch_length(node)
                        tot_bl = tot_bl + tree.branch_length(node)

                for xi in range(sample_size - 1):
                    rescaled_x = tbl[xi]
                    data["total_branch_length"].append(rescaled_x / tot_bl)
                    tbl_sum[xi] = tbl_sum[xi] + rescaled_x
                tot_bl_sum[0] = tot_bl_sum[0] + tot_bl
                data["num_leaves"].extend(range(1, sample_size))
        basedir = make_test_dir("xi_beta_expected_sfs")
        f = os.path.join(basedir, f"n={sample_size}_alpha={alpha}.png")
        ax = sns.violinplot(
            data=data, x="num_leaves", y="total_branch_length", color="grey"
        )
        ax.set_xlabel("num leaves")
        l1 = ax.plot(np.arange(sample_size - 1), sfs[::], "--", linewidth=3)
        l2 = ax.plot(
            np.arange(sample_size - 1),
            [(x / num_replicates) / (tot_bl_sum[0] / num_replicates) for x in tbl_sum],
            "--",
            linewidth=3,
        )
        ax.legend((l1[0], l2[0]), ("Expected", "Observed"))
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")

    def run_xi_dirac_expected_sfs(self):

        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.1,
            c=1,
            sfs=[
                0.35352303,
                0.17672997,
                0.11781921,
                0.08836481,
                0.07069227,
                0.05891075,
                0.05049574,
                0.04418514,
                0.03927908,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.3,
            c=1,
            sfs=[
                0.35430737,
                0.17650201,
                0.11762438,
                0.08822363,
                0.07058696,
                0.05883259,
                0.05044232,
                0.04416277,
                0.03931799,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.5,
            c=1,
            sfs=[
                0.35655911,
                0.17596878,
                0.11711820,
                0.08785514,
                0.07030139,
                0.05860142,
                0.05025410,
                0.04402755,
                0.03931431,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.9,
            c=1,
            sfs=[
                0.36443828,
                0.17490683,
                0.11614708,
                0.08717119,
                0.06965759,
                0.05790491,
                0.04939935,
                0.04279132,
                0.03758346,
            ],
        )

        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=3,
            psi=0.1,
            c=10,
            sfs=[0.6667343, 0.3332657],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=3,
            psi=0.3,
            c=10,
            sfs=[0.6682113, 0.3317887],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=3,
            psi=0.5,
            c=10,
            sfs=[0.6721853, 0.3278147],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=3,
            psi=0.9,
            c=10,
            sfs=[0.6852703, 0.3147297],
        )

        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.1,
            c=10,
            sfs=[
                0.35385062,
                0.17661522,
                0.11773706,
                0.08830646,
                0.07064941,
                0.05887993,
                0.05047626,
                0.04418035,
                0.03930470,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.3,
            c=10,
            sfs=[
                0.36053858,
                0.17456975,
                0.11610005,
                0.08713599,
                0.06977685,
                0.05822906,
                0.05002797,
                0.04398723,
                0.03963453,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.5,
            c=10,
            sfs=[
                0.37556917,
                0.17015781,
                0.11285655,
                0.08495119,
                0.06808802,
                0.05683977,
                0.04886055,
                0.04309158,
                0.03958537,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.9,
            c=10,
            sfs=[
                0.41154361,
                0.15908770,
                0.10852899,
                0.08341563,
                0.06647774,
                0.05471783,
                0.04592602,
                0.03818041,
                0.03212207,
            ],
        )

        ##########################################################################
        # Compare SFS when c=10000 to the expected SFS whetre c tend to infinity #
        ##########################################################################

        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.1,
            c=10000,
            sfs=[
                0.36939374,
                0.17057448,
                0.11408360,
                0.08571572,
                0.06874076,
                0.05749423,
                0.04958115,
                0.04390987,
                0.04050644,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.3,
            c=10000,
            sfs=[
                0.39876239,
                0.15840021,
                0.10834860,
                0.08165271,
                0.06562863,
                0.05508280,
                0.04777344,
                0.04280604,
                0.04154517,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.5,
            c=10000,
            sfs=[
                0.42603419,
                0.14512841,
                0.10505636,
                0.07956441,
                0.06368639,
                0.05328134,
                0.04595869,
                0.04078814,
                0.04050205,
            ],
        )
        self.compare_normalised_xi_dirac_sfs(
            num_replicates=10000,
            sample_size=10,
            psi=0.9,
            c=10000,
            sfs=[
                0.47543921,
                0.11338801,
                0.10691661,
                0.08342993,
                0.06358921,
                0.05162311,
                0.04334855,
                0.03416865,
                0.02809671,
            ],
        )

    def compare_xi_beta_sfs(self, sample_size, alpha, sfs, num_replicates=1000):
        """
        Runs simulations of the xi beta model and compares to the expected SFS.
        """
        logging.debug(f"running SFS for {sample_size} {alpha}")
        reps = msprime.simulate(
            sample_size,
            num_replicates=num_replicates,
            model=msprime.BetaCoalescent(alpha=alpha, truncation_point=1),
        )

        data = collections.defaultdict(list)
        tbl_sum = [0] * (sample_size - 1)
        for ts in reps:
            for tree in ts.trees():
                tot_bl = 0.0
                tbl = [0] * (sample_size - 1)
                for node in tree.nodes():
                    if tree.parent(node) != msprime.NULL_NODE:
                        tbl[tree.num_samples(node) - 1] = tbl[
                            tree.num_samples(node) - 1
                        ] + tree.branch_length(node)
                        tot_bl = tot_bl + tree.branch_length(node)

                for xi in range(sample_size - 1):
                    rescaled_x = tbl[xi] / tot_bl
                    data["total_branch_length"].append(rescaled_x)
                    tbl_sum[xi] = tbl_sum[xi] + rescaled_x
                data["num_leaves"].extend(range(1, sample_size))

        basedir = make_test_dir("xi_beta_expected_sfs")
        f = os.path.join(basedir, f"n={sample_size}_alpha={alpha}.png")
        ax = sns.violinplot(
            data=data, x="num_leaves", y="total_branch_length", color="grey"
        )
        ax.set_xlabel("num leaves")
        l1 = ax.plot(np.arange(sample_size - 1), sfs[::], "--", linewidth=3)
        l2 = ax.plot(
            np.arange(sample_size - 1),
            [x / num_replicates for x in tbl_sum],
            "--",
            linewidth=3,
        )
        ax.legend((l1[0], l2[0]), ("Expected", "Observed"))
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")

    def run_xi_beta_breakpoints(self):
        basedir_name = "xi_beta_breakpoints"
        Ne = 10 ** 4
        for alpha in [1.1, 1.3, 1.6, 1.9]:
            self.verify_breakpoint_distribution(
                basedir_name,
                f"n=100_alpha={alpha}",
                sample_size=100,
                Ne=Ne,
                r=1e-7,
                L=10 ** 6,
                model=msprime.BetaCoalescent(alpha=alpha),
            )
            # Add a growth rate with a higher recombination rate so
            # we still get decent numbers of trees
            self.verify_breakpoint_distribution(
                basedir_name,
                f"n=100_alpha={alpha}",
                sample_size=100,
                Ne=Ne,
                r=1e-7,
                L=10 ** 6,
                model=msprime.BetaCoalescent(alpha=alpha),
                growth_rate=0.05,
            )

    def run_xi_dirac_breakpoints(self):
        basedir_name = "xi_dirac_breakpoints"
        Ne = 10 ** 2
        for psi in [0.1, 0.3, 0.6, 0.9]:
            for c in [1, 10]:
                self.verify_breakpoint_distribution(
                    basedir_name,
                    f"n=100_psi={psi}_c={c}",
                    sample_size=100,
                    Ne=Ne,
                    r=1e-8,
                    L=10 ** 6,
                    model=msprime.DiracCoalescent(psi=psi, c=c),
                )
                # Add a growth rate with a higher recombination rate so
                # we still get decent numbers of trees
                self.verify_breakpoint_distribution(
                    basedir_name,
                    f"n=100_psi={psi}_c={c}",
                    sample_size=100,
                    Ne=Ne,
                    r=1e-7,
                    L=10 ** 6,
                    model=msprime.DiracCoalescent(psi=psi, c=c),
                    growth_rate=0.05,
                )

    def run_xi_beta_recombinations(self):
        basedir_name = "xi_beta_recombinations"
        Ne = 10000
        for alpha in [1.1, 1.3, 1.5, 1.9]:
            self.verify_recombination(
                basedir_name,
                f"n=100_alpha={alpha}",
                sample_size=100,
                Ne=Ne,
                r=1e-8,
                m=1e-8,
                L=10 ** 6,
                model=msprime.BetaCoalescent(alpha=alpha),
            )

    def run_xi_dirac_recombinations(self):
        basedir_name = "xi_dirac_recombinations"
        Ne = 100
        for psi in [0.1, 0.3, 0.5, 0.9]:
            for c in [1, 10, 100]:
                # The Dirac coalescent has branch lengths proportional to Ne^2
                # and recombination rate proportional to Ne so need to divide
                # the mutation rate by Ne to make numbers of mutations and
                # recombinations comparable.
                self.verify_recombination(
                    basedir_name,
                    f"n=100_psi={psi}_c={c}",
                    sample_size=100,
                    Ne=Ne,
                    r=1e-8,
                    m=1e-8,
                    L=10 ** 6,
                    model=msprime.DiracCoalescent(psi=psi, c=c),
                )

    def run_xi_beta_expected_sfs(self):

        self.compare_normalized_xi_beta_sfs(
            num_replicates=100000,
            sample_size=10,
            alpha=1.1,
            sfs=[
                0.40838865,
                0.15645421,
                0.10765060,
                0.08178884,
                0.06548874,
                0.05455910,
                0.04672861,
                0.04082172,
                0.03811953,
            ],
        )

        self.compare_normalized_xi_beta_sfs(
            num_replicates=100000,
            sample_size=10,
            alpha=1.3,
            sfs=[
                0.39612917,
                0.16173072,
                0.10932728,
                0.08270507,
                0.06630221,
                0.05534012,
                0.04754038,
                0.04182775,
                0.03909731,
            ],
        )

        self.compare_normalized_xi_beta_sfs(
            num_replicates=100000,
            sample_size=10,
            alpha=1.5,
            sfs=[
                0.38395732,
                0.16650213,
                0.11136301,
                0.08395003,
                0.06731437,
                0.05622960,
                0.04837457,
                0.04268961,
                0.03961935,
            ],
        )

        self.compare_normalized_xi_beta_sfs(
            num_replicates=100000,
            sample_size=10,
            alpha=1.9,
            sfs=[
                0.35961114,
                0.17486018,
                0.11638771,
                0.08734266,
                0.06992360,
                0.05832611,
                0.05007349,
                0.04396363,
                0.03951149,
            ],
        )


@attr.s
class XiTest1(XiTest):
    model = attr.ib(type=type(msprime.DiracCoalescent()))
    simulate_args = attr.ib(factory=dict)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        logging.debug(f"{self.name} {self.model} {self.simulate_args}")
        self.run_xi_hudson_comparison(self.name, self.model, **self.simulate_args)


@attr.s
class XiTest2(XiTest):
    # Adds a check for xi_beta recombination breakpoints

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_xi_beta_breakpoints()


@attr.s
class XiTest3(XiTest):
    # Adds a check for xi_dirac recombination breakpoints

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_xi_dirac_breakpoints()


@attr.s
class XiTest4(XiTest):
    # Adds a check for xi_dirac matching expected SFS calculations.

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_xi_dirac_expected_sfs()


@attr.s
class XiTest5(XiTest):
    # Adds a check for xi_beta matching expected SFS calculations.

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_xi_beta_expected_sfs()


@attr.s
class XiTest6(XiTest):
    # Adds a check for xi_beta recombination breakpoints

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_xi_beta_recombinations()


@attr.s
class XiTest7(XiTest):
    # Adds a check for xi_dirac recombination breakpoints

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_xi_dirac_recombinations()


def register_xi_tests(runner):
    group = "xi"
    tests = [
        # add_xi_dirac_vs_hudson_recombination
        # Checks Xi-dirac against the standard coalescent with recombination.
        XiTest1(
            "xi_dirac_vs_hudson_recombination",
            group,
            msprime.DiracCoalescent(psi=0.99, c=0),
            simulate_args=dict(
                sample_size=50, Ne=50, num_replicates=1000, recombination_rate=0.1,
            ),
        ),
        # add_xi_dirac_vs_hudson_single_locus
        # Checks Xi-dirac against the standard coalescent at a single locus.
        XiTest1(
            "xi_dirac_vs_hudson_single_locus",
            group,
            msprime.DiracCoalescent(psi=0.99, c=0),
            simulate_args=dict(sample_size=10, Ne=100, num_replicates=5000,),
        ),
        XiTest2(name="xi_beta_bp", group=group),
        XiTest3(name="xi_dirac_bp", group=group),
        XiTest4(name="xi_dirac_sfs", group=group),
        XiTest5(name="xi_beta_sfs", group=group),
        XiTest6(name="xi_beta_recomb", group=group),
        XiTest7(name="xi_beta_recomb", group=group),
    ]

    for test in tests:
        runner.register(test)


@attr.s
class HudsonTest(XiHudsonTest):
    def run_hudson_breakpoints(self):
        basedir_name = "hudson_breakpoints"
        self.verify_breakpoint_distribution(
            basedir_name,
            "single_pop_n_50",
            sample_size=50,
            Ne=10 ** 4,
            r=1e-8,
            L=10 ** 6,
            model="hudson",
        )
        self.verify_breakpoint_distribution(
            basedir_name,
            "single_pop_n_100",
            sample_size=100,
            Ne=10 ** 4,
            r=1e-8,
            L=10 ** 6,
            model="hudson",
        )
        # Add a growth rate with a higher recombination rate so
        # we still get decent numbers of trees
        self.verify_breakpoint_distribution(
            basedir_name,
            "single_pop_n_100_growth",
            sample_size=100,
            Ne=10 ** 4,
            r=1e-7,
            L=10 ** 6,
            model="hudson",
            growth_rate=0.05,
        )

    def run_Hudson_recombinations(self):
        basedir_name = "hudson_recombinations"
        self.verify_recombination(
            basedir_name,
            f"n=100_hudson",
            sample_size=100,
            Ne=10000,
            r=1e-8,
            m=1e-8,
            L=10 ** 6,
            model="hudson",
        )


@attr.s
class HudsonTest1(HudsonTest):
    """
    Adds a check for hudson recombination breakpoints
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_Hudson_recombinations()


@attr.s
class HudsonTest2(HudsonTest):
    """
    Adds a check for hudson recombination breakpoints
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_hudson_breakpoints()


def register_hudson_tests(runner):
    group = "hudson"
    tests = [
        HudsonTest1(name="hudson recombination", group=group),
        HudsonTest2(name="hudson breakpoints", group=group),
    ]

    for test in tests:
        runner.register(test)


@attr.s
class ContDiscreteTest(Test):
    def _run_msprime_coalescent_stats(self, **kwargs):
        logging.debug(f"\t msprime: {kwargs}")
        if "num_replicates" in kwargs:
            replicates = kwargs["num_replicates"]
            num_trees = [0 for i in range(replicates)]
            breakpoints = [0 for i in range(replicates)]
            for i, ts in enumerate(msprime.simulate(**kwargs)):
                num_trees[i] = ts.num_trees
                breakpoints[i] = list(ts.breakpoints())
        else:
            ts = msprime.simulate(**kwargs)
            num_trees = [ts.num_trees]
            breakpoints = [list(ts.breakpoints)]

        d = {"num_trees": num_trees, "breakpoints": breakpoints}
        df = pd.DataFrame(d)
        return df

    def run_cont_discrete_comparison(
        self, key, model, discrete_recomb_map, cont_recomb_map
    ):
        sample_size = 10
        num_replicates = 400
        df_discrete = self._run_msprime_coalescent_stats(
            num_replicates=num_replicates,
            sample_size=sample_size,
            model=model,
            recombination_map=discrete_recomb_map,
        )
        df_cont = self._run_msprime_coalescent_stats(
            num_replicates=num_replicates,
            sample_size=sample_size,
            model=model,
            recombination_map=cont_recomb_map,
        )

        discrete_length = discrete_recomb_map.get_sequence_length()
        cont_length = cont_recomb_map.get_sequence_length()
        scale_breakpoints(df_cont, discrete_length / cont_length)
        self._plot_stats(
            key,
            "compare continuous and discrete coordinates",
            df_discrete,
            df_cont,
            "discrete",
            "continuous",
        )

    def run_uniform_recomb_cont_discrete_comparison(self, key, model):
        discrete_recomb_map = msprime.RecombinationMap.uniform_map(
            2000000, 1e-5, discrete=True
        )
        cont_recomb_map = msprime.RecombinationMap.uniform_map(
            1, 2000000 * 1e-5, discrete=False
        )

        self.run_cont_discrete_comparison(
            key, model, discrete_recomb_map, cont_recomb_map
        )

    def run_variable_recomb_cont_discrete_comparison(self, key, model):
        r = 1e-5
        discrete_positions = [0, 10000, 50000, 150000, 200000]
        discrete_rates = [0.0, r, 5 * r, r / 2, 0.0]
        cont_positions = [x / 200000 for x in discrete_positions]
        cont_rates = [x * 200000 for x in discrete_rates]

        discrete_recomb_map = msprime.RecombinationMap(
            discrete_positions, discrete_rates, discrete=True
        )

        cont_recomb_map = msprime.RecombinationMap(
            cont_positions, cont_rates, discrete=False
        )

        self.run_cont_discrete_comparison(
            key, model, discrete_recomb_map, cont_recomb_map
        )

    def run_continuous_discrete_same_scale(self, key, model):
        discrete_recomb_map = msprime.RecombinationMap.uniform_map(
            2000000, 1e-5, discrete=True
        )
        cont_recomb_map = msprime.RecombinationMap.uniform_map(
            2000000, 1e-5, discrete=False
        )
        self.run_cont_discrete_comparison(
            key, model, discrete_recomb_map, cont_recomb_map
        )


@attr.s
class ContDiscreteTest1(ContDiscreteTest):
    model = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_continuous_discrete_same_scale(self.name, self.model)


@attr.s
class ContDiscreteTest2(ContDiscreteTest):
    """
    Adds checks comparing equivalent simulations in discrete space
    and scaled up continuous space.
    """

    model = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_uniform_recomb_cont_discrete_comparison(
            self.name + "_uniform", self.model
        )
        self.run_variable_recomb_cont_discrete_comparison(
            self.name + "_variable", self.model
        )


def register_contdiscrete_tests(runner):
    group = "contdiscrete"
    tests = [
        ContDiscreteTest1("hudson_cont_discrete_same_scale", group, "hudson"),
        ContDiscreteTest1("dtwf_cont_discrete_same_scale", group, "dtwf"),
        ContDiscreteTest2("hudson_recomb_cont_discrete", group, "hudson"),
        ContDiscreteTest2("dtwf_recomb_cont_discrete", group, "dtwf"),
    ]
    for test in tests:
        runner.register(test)


@attr.s
class ArgRecordTest(Test):
    def run_arg_recording(self):
        basedir = make_test_dir("arg_recording")

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
                sample_size=leaves, recombination_rate=rho, random_seed=i + 1
            )
            ts_node_counts = np.append(ts_node_counts, ts.num_nodes)
            ts_tree_counts = np.append(ts_tree_counts, ts.num_trees)
            ts_edge_counts = np.append(ts_edge_counts, ts.num_edges)
            arg = msprime.simulate(
                sample_size=leaves,
                recombination_rate=rho,
                random_seed=i + 1,
                record_full_arg=True,
            )
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
        pyplot.close("all")

    def run_multiple_merger_arg_recording(self):
        basedir = make_test_dir("mmc_arg_recording")

        ts_node_counts = np.array([])
        arg_node_counts = np.array([])
        ts_tree_counts = np.array([])
        arg_tree_counts = np.array([])
        ts_edge_counts = np.array([])
        arg_edge_counts = np.array([])

        reps = 10000
        leaves = 100
        rho = 2

        for i in range(reps):
            ts = msprime.simulate(
                sample_size=leaves,
                recombination_rate=rho,
                random_seed=i + 1,
                model=msprime.BetaCoalescent(alpha=1.1, truncation_point=1),
            )
            ts_node_counts = np.append(ts_node_counts, ts.num_nodes)
            ts_tree_counts = np.append(ts_tree_counts, ts.num_trees)
            ts_edge_counts = np.append(ts_edge_counts, ts.num_edges)
            arg = msprime.simulate(
                sample_size=leaves,
                recombination_rate=rho,
                random_seed=i + 1,
                model=msprime.BetaCoalescent(alpha=1.1, truncation_point=1),
                record_full_arg=True,
            )
            arg = arg.simplify()
            arg_node_counts = np.append(arg_node_counts, arg.num_nodes)
            arg_tree_counts = np.append(arg_tree_counts, arg.num_trees)
            arg_edge_counts = np.append(arg_edge_counts, arg.num_edges)

        pp_ts = sm.ProbPlot(ts_node_counts)
        pp_arg = sm.ProbPlot(arg_node_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "beta_nodes.png")
        pyplot.savefig(f, dpi=72)

        pp_ts = sm.ProbPlot(ts_tree_counts)
        pp_arg = sm.ProbPlot(arg_tree_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "beta_trees.png")
        pyplot.savefig(f, dpi=72)

        pp_ts = sm.ProbPlot(ts_edge_counts)
        pp_arg = sm.ProbPlot(arg_edge_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "beta_edges.png")
        pyplot.savefig(f, dpi=72)

        ts_node_counts = np.array([])
        arg_node_counts = np.array([])
        ts_tree_counts = np.array([])
        arg_tree_counts = np.array([])
        ts_edge_counts = np.array([])
        arg_edge_counts = np.array([])

        for i in range(reps):
            ts = msprime.simulate(
                sample_size=leaves,
                recombination_rate=rho,
                random_seed=i + 1,
                model=msprime.DiracCoalescent(psi=0.9, c=1),
            )
            ts_node_counts = np.append(ts_node_counts, ts.num_nodes)
            ts_tree_counts = np.append(ts_tree_counts, ts.num_trees)
            ts_edge_counts = np.append(ts_edge_counts, ts.num_edges)
            arg = msprime.simulate(
                sample_size=leaves,
                recombination_rate=rho,
                random_seed=i + 1,
                model=msprime.DiracCoalescent(psi=0.9, c=1),
                record_full_arg=True,
            )
            arg = arg.simplify()
            arg_node_counts = np.append(arg_node_counts, arg.num_nodes)
            arg_tree_counts = np.append(arg_tree_counts, arg.num_trees)
            arg_edge_counts = np.append(arg_edge_counts, arg.num_edges)

        pp_ts = sm.ProbPlot(ts_node_counts)
        pp_arg = sm.ProbPlot(arg_node_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "dirac_nodes.png")
        pyplot.savefig(f, dpi=72)

        pp_ts = sm.ProbPlot(ts_tree_counts)
        pp_arg = sm.ProbPlot(arg_tree_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "dirac_trees.png")
        pyplot.savefig(f, dpi=72)

        pp_ts = sm.ProbPlot(ts_edge_counts)
        pp_arg = sm.ProbPlot(arg_edge_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        f = os.path.join(basedir, "dirac_edges.png")
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")


@attr.s
class ArgRecordTest1(ArgRecordTest):
    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_arg_recording()


@attr.s
class ArgRecordTest2(ArgRecordTest):
    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_multiple_merger_arg_recording()


def register_argrecord_tests(runner):
    group = "argrecord"
    tests = [
        # Check that we get the right number of objects when we simplify
        # a full arg.
        ArgRecordTest1(name="simple check", group=group),
        ArgRecordTest2(name="merge check", group=group),
    ]

    for test in tests:
        runner.register(test)


@attr.s
class AnalyticalTest(Test):
    def get_segregating_sites_histogram(self, cmd):
        logging.debug(f"\t {' '.join(cmd)}")
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
            t1 = (-1) ** i
            t2 = scipy.special.binom(n - 1, i - 1)
            t3 = (i - 1) / (theta + i - 1)
            t4 = (theta / (theta + i - 1)) ** k
            s += t1 * t2 * t3 * t4
        return s

    def run_s_analytical_check(self):
        """
        Runs the check for the number of segregating sites against the
        analytical prediction.
        """
        R = 100000
        theta = 2
        basedir = make_test_dir("analytical_s")
        for n in range(2, 15):
            cmd = f"{n} {R} -t {theta}"
            S_ms = self.get_segregating_sites_histogram(
                self._ms_executable + cmd.split() + self.get_ms_seeds()
            )
            S_msp = self.get_segregating_sites_histogram(
                self._mspms_executable + cmd.split() + self.get_ms_seeds()
            )
            filename = os.path.join(basedir, f"{n}.png")

            fig, ax = pyplot.subplots()
            index = np.arange(10)
            S_analytical = [self.get_S_distribution(j, n, theta) for j in index]
            bar_width = 0.35
            pyplot.bar(index, S_ms[index], bar_width, color="b", label="ms")
            pyplot.bar(
                index + bar_width, S_msp[index], bar_width, color="r", label="msp"
            )
            pyplot.plot(index + bar_width, S_analytical, "o", color="k")
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
        basedir = make_test_dir("analytical_pi")

        sample_size = np.arange(2, 15)
        mean = np.zeros_like(sample_size, dtype=float)
        var = np.zeros_like(sample_size, dtype=float)
        predicted_mean = np.zeros_like(sample_size, dtype=float)
        predicted_var = np.zeros_like(sample_size, dtype=float)

        for k, n in enumerate(sample_size):
            pi = np.zeros(R)
            replicates = msprime.simulate(
                sample_size=n, mutation_rate=theta / 4, num_replicates=R
            )
            for j, ts in enumerate(replicates):
                pi[j] = ts.get_pairwise_diversity()
            # Predicted mean is is theta.
            predicted_mean[k] = theta
            # From Wakely, eqn (4.14), pg. 101
            predicted_var[k] = (n + 1) * theta / (3 * (n - 1)) + 2 * (
                n ** 2 + n + 3
            ) * theta ** 2 / (9 * n * (n - 1))
            mean[k] = np.mean(pi)
            var[k] = np.var(pi)

            logging.debug(
                f"{n}\t{theta}\t{np.mean(pi)}\t{predicted_var[k]}\t{np.var(pi)}"
            )

        filename = os.path.join(basedir, "mean.png")
        pyplot.plot(sample_size, predicted_mean, "-")
        pyplot.plot(sample_size, mean, "-")
        pyplot.savefig(filename)
        pyplot.close("all")

        filename = os.path.join(basedir, "var.png")
        pyplot.plot(sample_size, predicted_var, "-")
        pyplot.plot(sample_size, var, "-")
        pyplot.savefig(filename)
        pyplot.close("all")

    def run_correlation_between_trees_analytical_check(self):
        """
        Runs the check for the probability of same tree at two sites against
        analytical predictions.
        """
        R = 1000
        basedir = make_test_dir("analytical_corr_same_tree")

        sample_size = 2
        gc_length_rate_ratio = np.array([0.05, 0.5, 5.0])
        gc_length = np.array([100, 50, 20])
        gc_rate = 1.0 / (gc_length_rate_ratio * gc_length)
        seq_length = 500
        predicted_prob = np.zeros([gc_length_rate_ratio.size, seq_length], dtype=float)
        empirical_prob_first = np.zeros(
            [gc_length_rate_ratio.size, seq_length], dtype=float
        )
        empirical_prob_mid = np.zeros(
            [gc_length_rate_ratio.size, seq_length], dtype=float
        )
        empirical_prob_last = np.zeros(
            [gc_length_rate_ratio.size, seq_length], dtype=float
        )

        for k, l in enumerate(gc_length):
            same_root_count_first = np.zeros(seq_length)
            same_root_count_mid = np.zeros(seq_length)
            same_root_count_last = np.zeros(seq_length)
            low_recombination_rate = 0.000001
            recomb_map = msprime.RecombinationMap.uniform_map(
                seq_length, low_recombination_rate, discrete=True
            )
            replicates = msprime.simulate(
                sample_size=sample_size,
                recombination_map=recomb_map,
                gene_conversion_rate=gc_rate[k],
                gene_conversion_track_length=gc_length[k],
                num_replicates=R,
            )
            for ts in replicates:
                firstroot = ts.first().root
                lastroot = ts.last().root
                for tree in ts.trees():
                    left, right = tree.interval
                    if left <= seq_length / 2 < right:
                        midroot = tree.root
                for tree in ts.trees():
                    left, right = map(int, tree.interval)
                    if firstroot == tree.root:
                        same_root_count_first[left:right] += 1
                    if lastroot == tree.root:
                        same_root_count_last[left:right] += 1
                    if midroot == tree.root:
                        same_root_count_mid[left:right] += 1
            empirical_prob_first[k, :] = same_root_count_first / R
            empirical_prob_last[k, :] = same_root_count_last / R
            empirical_prob_mid[k, :] = same_root_count_mid / R
            # Predicted prob
            # From Wiuf, Hein, 2000, eqn (15), pg. 457
            rG = (
                2 / gc_length_rate_ratio[k] * (1.0 - np.exp(-np.arange(seq_length) / l))
            )
            predicted_prob[k, :] = (18.0 + rG) / (18.0 + 13.0 * rG + rG * rG)

        x = np.arange(500) + 1
        filename = os.path.join(basedir, "prob_first.png")
        pyplot.plot(x, predicted_prob[0], "--")
        pyplot.plot(x, empirical_prob_first[0], "-")
        pyplot.plot(x, predicted_prob[1], "--")
        pyplot.plot(x, empirical_prob_first[1], "-")
        pyplot.plot(x, predicted_prob[2], "--")
        pyplot.plot(x, empirical_prob_first[2], "-")
        pyplot.savefig(filename)
        pyplot.close("all")

        filename = os.path.join(basedir, "prob_last.png")
        pyplot.plot(x, predicted_prob[0, ::-1], "--")
        pyplot.plot(x, empirical_prob_last[0], "-")
        pyplot.plot(x, predicted_prob[1, ::-1], "--")
        pyplot.plot(x, empirical_prob_last[1], "-")
        pyplot.plot(x, predicted_prob[2, ::-1], "--")
        pyplot.plot(x, empirical_prob_last[2], "-")
        pyplot.savefig(filename)
        pyplot.close("all")

        filename = os.path.join(basedir, "prob_mid.png")
        pyplot.plot(
            x,
            np.concatenate((predicted_prob[0, 249::-1], predicted_prob[0, :250])),
            "--",
        )
        pyplot.plot(x, empirical_prob_mid[0], "-")
        pyplot.plot(
            x,
            np.concatenate((predicted_prob[1, 249::-1], predicted_prob[1, :250])),
            "--",
        )
        pyplot.plot(x, empirical_prob_mid[1], "-")
        pyplot.plot(
            x,
            np.concatenate((predicted_prob[2, 249::-1], predicted_prob[2, :250])),
            "--",
        )
        pyplot.plot(x, empirical_prob_mid[2], "-")
        pyplot.savefig(filename)
        pyplot.close("all")

        filename = os.path.join(basedir, "prob_first_zoom.png")
        x = np.arange(10) + 1
        pyplot.plot(x, predicted_prob[0, range(10)], "--")
        pyplot.plot(x, empirical_prob_first[0, range(10)], "-")
        pyplot.plot(x, predicted_prob[1, range(10)], "--")
        pyplot.plot(x, empirical_prob_first[1, range(10)], "-")
        pyplot.plot(x, predicted_prob[2, range(10)], "--")
        pyplot.plot(x, empirical_prob_first[2, range(10)], "-")
        pyplot.savefig(filename)
        pyplot.close("all")

    def run_mean_coaltime_check(self):
        """
        Checks the mean coalescence time calculation against pi.
        """
        random.seed(5)
        num_models = 8
        num_reps = 8
        T = np.zeros((num_models, num_reps))
        U = np.zeros(num_models)
        logging.debug("coaltime: theory  mean  sd   z")
        for k in range(num_models):
            Ne = 100
            N = 4
            pop_sizes = [random.uniform(0.01, 10) * Ne for _ in range(N)]
            growth_rates = [random.uniform(-0.01, 0.01) for _ in range(N)]
            migration_matrix = [
                [random.random() * (i != j) for j in range(N)] for i in range(N)
            ]
            sample_size = [random.randint(2, 10) for _ in range(N)]
            population_configurations = [
                msprime.PopulationConfiguration(
                    initial_size=k, sample_size=n, growth_rate=r
                )
                for k, n, r in zip(pop_sizes, sample_size, growth_rates)
            ]
            demographic_events = []
            for i in [0, 1]:
                n = random.uniform(0.01, 10)
                r = 0
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=100, initial_size=n, growth_rate=r, population_id=i
                    )
                )
            for ij in [(0, 1), (2, 3), (0, 3)]:
                demographic_events.append(
                    msprime.MigrationRateChange(180, random.random(), matrix_index=ij)
                )
            demographic_events.append(
                msprime.MassMigration(time=200, source=3, dest=0, proportion=0.3)
            )
            for i in [1, 3]:
                n = random.uniform(0.01, 10)
                r = random.uniform(-0.01, 0.01)
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=210, initial_size=n, growth_rate=r, population_id=i
                    )
                )

            ddb = msprime.DemographyDebugger(
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                migration_matrix=migration_matrix,
            )

            U[k] = ddb.mean_coalescence_time(num_samples=sample_size)

            mut_rate = 1e-8
            replicates = msprime.simulate(
                length=1e7,
                recombination_rate=1e-8,
                mutation_rate=mut_rate,
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                migration_matrix=migration_matrix,
                random_seed=5,
                num_replicates=num_reps,
            )
            for j, ts in enumerate(replicates):
                T[k, j] = ts.get_pairwise_diversity()
                T[k, j] /= ts.sequence_length
                T[k, j] /= 2 * mut_rate
            mT = np.mean(T[k])
            sT = np.std(T[k])
            logging.debug(
                "        {:.2f} {:.2f} {:.2f} {:.2f}".format(
                    U[k], mT, sT, (U[k] - mT) / (sT * np.sqrt(num_reps))
                )
            )

        basedir = make_test_dir("coaltime")
        fig, ax = pyplot.subplots()
        ax.scatter(np.column_stack([U] * T.shape[1]), T)
        # where oh where is abline(0,1)
        line = mlines.Line2D([0, 1], [0, 1])
        line.set_transform(ax.transAxes)
        ax.add_line(line)
        ax.set_xlabel("calculated mean coaltime")
        ax.set_ylabel("pairwise diversity, scaled")
        filename = os.path.join(basedir, "mean_coaltimes.png")
        pyplot.savefig(filename)
        pyplot.close("all")

    def get_tbl_distribution(self, n, R, executable):
        """
        Returns an array of the R total branch length values from
        the specified ms-like executable.
        """
        cmd = executable + f"{n} {R} -T -p 10".split()
        cmd += self.get_ms_seeds()
        logging.debug(f"\t {' '.join(cmd)}")
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
        basedir = make_test_dir("analytical_tbl")
        for n in range(2, 15):
            tbl_ms = self.get_tbl_distribution(n, R, self._ms_executable)
            tbl_msp = self.get_tbl_distribution(n, R, self._mspms_executable)

            sm.graphics.qqplot(tbl_ms)
            sm.qqplot_2samples(tbl_ms, tbl_msp, line="45")
            filename = os.path.join(basedir, f"qqplot_{n}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            hist_ms, bin_edges = np.histogram(tbl_ms, 20, density=True)
            hist_msp, _ = np.histogram(tbl_msp, bin_edges, density=True)

            index = bin_edges[:-1]
            # We don't seem to have the analytical value quite right here,
            # but since the value is so very close to ms's, there doesn't
            # seem to be much point in trying to fix it.
            analytical = [self.get_analytical_tbl(n, x * 2) for x in index]
            fig, ax = pyplot.subplots()
            bar_width = 0.15
            pyplot.bar(index, hist_ms, bar_width, color="b", label="ms")
            pyplot.bar(index + bar_width, hist_msp, bar_width, color="r", label="msp")
            pyplot.plot(index + bar_width, analytical, "o", color="k")
            pyplot.legend()
            # pyplot.xticks(index + bar_width, [str(j) for j in index])
            pyplot.tight_layout()
            filename = os.path.join(basedir, f"hist_{n}.png")
            pyplot.savefig(filename)

    def get_num_trees(self, cmd, R):
        logging.debug(f"\t {' '.join(cmd)}")
        output = subprocess.check_output(cmd)
        T = np.zeros(R)
        j = -1
        for line in output.splitlines():
            if line.startswith(b"//"):
                j += 1
            if line.startswith(b"["):
                T[j] += 1
        return T

    def run_cli_num_trees(self):
        """
        Runs the check for number of trees using the CLI.
        """
        r = 1e-8  # Per generation recombination rate.
        num_loci = np.linspace(100, 10 ** 5, 10).astype(int)
        Ne = 10 ** 4
        n = 100
        rho = r * 4 * Ne * (num_loci - 1)
        num_replicates = 100
        ms_mean = np.zeros_like(rho)
        msp_mean = np.zeros_like(rho)
        for j in range(len(num_loci)):
            cmd = "{} {} -T -r {} {}".format(n, num_replicates, rho[j], num_loci[j])
            T = self.get_num_trees(
                self._ms_executable + cmd.split() + self.get_ms_seeds(), num_replicates
            )
            ms_mean[j] = np.mean(T)

            T = self.get_num_trees(
                self._mspms_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates,
            )
            msp_mean[j] = np.mean(T)
        basedir = make_test_dir("cli_num_trees")
        pyplot.plot(rho, ms_mean, "o")
        pyplot.plot(rho, msp_mean, "^")
        pyplot.plot(rho, rho * harmonic_number(n - 1), "-")
        filename = os.path.join(basedir, "mean.png")
        pyplot.savefig(filename)
        pyplot.close("all")

    def get_pairwise_coalescence_time(self, cmd, R):
        # logging.debug(f"\t {' '.join(cmd)}")
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

    def run_pairwise_island_model(self):
        """
        Runs the check for the pairwise coalscence times for within
        and between populations.
        """
        R = 10000
        M = 0.2
        basedir = make_test_dir("analytical_pairwise_island")

        for d in range(2, 6):
            cmd = "2 {} -T -I {} 2 {} {}".format(R, d, "0 " * (d - 1), M)
            T_w_ms = self.get_pairwise_coalescence_time(
                self._ms_executable + cmd.split() + self.get_ms_seeds(), R
            )
            T_w_msp = self.get_pairwise_coalescence_time(
                self._mspms_executable + cmd.split() + self.get_ms_seeds(), R
            )

            cmd = "2 {} -T -I {} 1 1 {} {}".format(R, d, "0 " * (d - 2), M)
            T_b_ms = self.get_pairwise_coalescence_time(
                self._ms_executable + cmd.split() + self.get_ms_seeds(), R
            )
            T_b_msp = self.get_pairwise_coalescence_time(
                self._mspms_executable + cmd.split() + self.get_ms_seeds(), R
            )
            logging.debug(
                f"\
                {d}\t\
                {np.mean(T_w_ms)}\t\
                {np.mean(T_w_msp)}\t\
                {d / 2}\t\
                {np.mean(T_b_ms)}\t\
                {np.mean(T_b_msp)}\t\
                {(d + (d - 1) / M) / 2}"
            )

            sm.graphics.qqplot(T_w_ms)
            sm.qqplot_2samples(T_w_ms, T_w_msp, line="45")
            f = os.path.join(basedir, f"within_{d}.png")
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(T_b_ms)
            sm.qqplot_2samples(T_b_ms, T_b_msp, line="45")
            f = os.path.join(basedir, f"between_{d}.png")
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")


@attr.s
class AnalyticalTest1(AnalyticalTest):
    """
    Adds a check for the analytical predictions about the distribution
    of S, the number of segregating sites.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_s_analytical_check()


@attr.s
class AnalyticalTest2(AnalyticalTest):
    """
    Adds a check for the analytical predictions about the pi,
    the pairwise site diversity.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_pi_analytical_check()


@attr.s
class AnalyticalTest3(AnalyticalTest):
    """
    Adds a check for the analytical predictions about the correlation between
    trees in the case of gene conversion.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_correlation_between_trees_analytical_check()


@attr.s
class AnalyticalTest4(AnalyticalTest):
    """
    Adds a check for the demography debugger predictions about
    mean coalescence time.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_mean_coaltime_check()


@attr.s
class AnalyticalTest5(AnalyticalTest):
    """
    Adds a check for the analytical check for the total branch length.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_tbl_analytical_check()


@attr.s
class AnalyticalTest6(AnalyticalTest):
    """
    Adds a check for the analytical check for pairwise island model
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_pairwise_island_model()


@attr.s
class AnalyticalTest7(AnalyticalTest):
    """
    Adds a check for the analytical number of trees using the CLI
    and comparing with ms.add_s_analytical_check
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_cli_num_trees()


def register_analytical_tests(runner):
    group = "analytical"
    tests = [
        AnalyticalTest1(name="s", group=group),
        AnalyticalTest2(name="pi", group=group),
        AnalyticalTest3(name="corr_trees", group=group),
        AnalyticalTest4(name="mean coaltime", group=group),
        AnalyticalTest5(name="branch length", group=group),
        AnalyticalTest6(name="island", group=group),
        AnalyticalTest7(name="cli num trees", group=group),
    ]

    for test in tests:
        runner.register(test)


@attr.s
class SmcTest(Test):
    _scrm_executable = attr.ib(init=False, default=["./data/scrm"])

    def get_scrm_num_trees(self, cmd, R):
        logging.debug(f"\t {' '.join(cmd)}")
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
        logging.debug(f"\t {' '.join(cmd)}")
        output = subprocess.check_output(cmd)
        T = np.zeros(R)
        j = -1
        for line in output.splitlines():
            if line.startswith(b"//"):
                j += 1
            if line.startswith(b"time:"):
                T[j] = max(T[j], float(line.split()[1]))
        return T

    def run_smc_oldest_time(self):
        """
        Runs the check for number of trees using the CLI.
        """
        r = 1e-8  # Per generation recombination rate.
        num_loci = np.linspace(100, 10 ** 5, 10).astype(int)
        Ne = 10 ** 4
        n = 100
        rho = r * 4 * Ne * (num_loci - 1)
        num_replicates = 1000
        scrm_mean = np.zeros_like(rho)
        scrm_smc_mean = np.zeros_like(rho)
        msp_mean = np.zeros_like(rho)
        msp_smc_mean = np.zeros_like(rho)
        for j in range(len(num_loci)):

            cmd = "{} {} -L -r {} {} -p 14".format(
                n, num_replicates, rho[j], num_loci[j]
            )
            T = self.get_scrm_oldest_time(
                self._scrm_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates,
            )
            scrm_mean[j] = np.mean(T)

            cmd += " -l 0"
            T = self.get_scrm_oldest_time(
                self._scrm_executable + cmd.split() + self.get_ms_seeds(),
                num_replicates,
            )
            scrm_smc_mean[j] = np.mean(T)

            for dest, model in [(msp_mean, "hudson"), (msp_smc_mean, "smc_prime")]:
                replicates = msprime.simulate(
                    sample_size=n,
                    length=num_loci[j],
                    recombination_rate=r,
                    Ne=Ne,
                    num_replicates=num_replicates,
                    model=model,
                )
                T = np.zeros(num_replicates)
                for k, ts in enumerate(replicates):
                    for record in ts.records():
                        T[k] = max(T[k], record.time)
                # Normalise back to coalescent time.
                T /= 4 * Ne
                dest[j] = np.mean(T)
        basedir = make_test_dir("smc_oldest_time")
        pyplot.plot(rho, scrm_mean, "-", color="blue", label="scrm")
        pyplot.plot(rho, scrm_smc_mean, "-", color="red", label="scrm_smc")
        pyplot.plot(rho, msp_smc_mean, "--", color="red", label="msprime_smc")
        pyplot.plot(rho, msp_mean, "--", color="blue", label="msprime")
        pyplot.xlabel("rho")
        pyplot.ylabel("Mean oldest coalescence time")
        pyplot.legend(loc="lower right")
        filename = os.path.join(basedir, "mean.png")
        pyplot.savefig(filename)
        pyplot.close("all")

    def run_smc_num_trees(self):
        """
        Runs the check for number of trees in the SMC and full coalescent
        using the API. We compare this with scrm using the SMC as a check.
        """
        r = 1e-8  # Per generation recombination rate.
        L = np.linspace(100, 10 ** 5, 10).astype(int)
        Ne = 10 ** 4
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
                num_replicates,
            )
            mean_scrm[j] = np.mean(T)
            var_scrm[j] = np.var(T)
            # IMPORTANT!! We have to use the get_num_breakpoints method
            # on the simulator as there is a significant drop in the number
            # of trees if we use the tree sequence. There is a significant
            # number of common ancestor events that result in a recombination
            # being undone.
            exact_sim = msprime.simulator_factory(
                sample_size=n, recombination_rate=r, Ne=Ne, length=L[j]
            )
            for k in range(num_replicates):
                exact_sim.run()
                num_trees[k] = exact_sim.num_breakpoints
                exact_sim.reset()
            mean_exact[j] = np.mean(num_trees)
            var_exact[j] = np.var(num_trees)

            smc_sim = msprime.simulator_factory(
                sample_size=n, recombination_rate=r, Ne=Ne, length=L[j], model="smc"
            )
            for k in range(num_replicates):
                smc_sim.run()
                num_trees[k] = smc_sim.num_breakpoints
                smc_sim.reset()
            mean_smc[j] = np.mean(num_trees)
            var_smc[j] = np.var(num_trees)

            smc_prime_sim = msprime.simulator_factory(
                sample_size=n,
                recombination_rate=r,
                Ne=Ne,
                length=L[j],
                model="smc_prime",
            )
            for k in range(num_replicates):
                smc_prime_sim.run()
                num_trees[k] = smc_prime_sim.num_breakpoints
                smc_prime_sim.reset()
            mean_smc_prime[j] = np.mean(num_trees)
            var_smc_prime[j] = np.var(num_trees)

        basedir = make_test_dir("smc_num_trees")

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
        pyplot.close("all")

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
        pyplot.close("all")


@attr.s
class SmcTest1(SmcTest):
    """
    Adds a check for the analytical number of trees under the SMC
    and the full coalescent.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_smc_num_trees()


@attr.s
class SmcTest2(SmcTest):
    """
    Adds a check the distribution of the oldest time of a
    coalescence in the smc using scrm.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_smc_oldest_time()


def register_smc_tests(runner):
    group = "smc"
    tests = [
        SmcTest1(name="num trees", group=group),
        SmcTest2(name="time", group=group),
    ]
    for test in tests:
        runner.register(test)


@attr.s
class SimTest(Test):
    def run_simulate_from_single_locus(self):
        num_replicates = 1000

        basedir = make_test_dir("simulate_from_single_locus")

        for n in [10, 50, 100, 200]:
            logging.debug(f"running for n = {n}")
            T1 = np.zeros(num_replicates)
            reps = msprime.simulate(n, num_replicates=num_replicates)
            for j, ts in enumerate(reps):
                T1[j] = np.max(ts.tables.nodes.time)

            for t in [0.5, 1, 1.5, 5]:
                T2 = np.zeros(num_replicates)
                reps = msprime.simulate(n, num_replicates=num_replicates, end_time=t)
                for j, ts in enumerate(reps):
                    final_ts = msprime.simulate(
                        from_ts=ts, start_time=np.max(ts.tables.nodes.time)
                    )
                    final_ts = final_ts.simplify()
                    T2[j] = np.max(final_ts.tables.nodes.time)

                sm.graphics.qqplot(T1)
                sm.qqplot_2samples(T1, T2, line="45")
                filename = os.path.join(basedir, f"T_mrca_n={n}_t={t}.png")
                pyplot.savefig(filename, dpi=72)
                pyplot.close("all")

    def run_simulate_from_multi_locus(self):
        num_replicates = 1000
        n = 100

        basedir = make_test_dir("simulate_from_multi_locus")

        for m in [10, 50, 100, 1000]:
            logging.debug(f"running for m = {m}")
            T1 = np.zeros(num_replicates)
            num_trees1 = np.zeros(num_replicates)
            recomb_map = msprime.RecombinationMap.uniform_map(m, 1 / m, discrete=True)
            reps = msprime.simulate(
                n, recombination_map=recomb_map, num_replicates=num_replicates
            )
            for j, ts in enumerate(reps):
                T1[j] = np.max(ts.tables.nodes.time)
                num_trees1[j] = ts.num_trees

            for t in [0.5, 1, 1.5, 5]:
                T2 = np.zeros(num_replicates)
                num_trees2 = np.zeros(num_replicates)
                reps = msprime.simulate(
                    n,
                    num_replicates=num_replicates,
                    recombination_map=recomb_map,
                    end_time=t,
                )
                for j, ts in enumerate(reps):
                    final_ts = msprime.simulate(
                        from_ts=ts,
                        recombination_map=recomb_map,
                        start_time=np.max(ts.tables.nodes.time),
                    )
                    final_ts = final_ts.simplify()
                    T2[j] = np.max(final_ts.tables.nodes.time)
                    num_trees2[j] = final_ts.num_trees

                sm.graphics.qqplot(T1)
                sm.qqplot_2samples(T1, T2, line="45")
                filename = os.path.join(basedir, f"T_mrca_m={m}_t={t}.png")
                pyplot.savefig(filename, dpi=72)
                pyplot.close("all")

                sm.graphics.qqplot(num_trees1)
                sm.qqplot_2samples(num_trees1, num_trees2, line="45")
                filename = os.path.join(basedir, f"num_trees_m={m}_t={t}.png")
                pyplot.savefig(filename, dpi=72)
                pyplot.close("all")

    def run_simulate_from_recombination(self):
        num_replicates = 1000
        n = 100
        recombination_rate = 10

        basedir = make_test_dir("simulate_from_recombination")

        T1 = np.zeros(num_replicates)
        num_trees1 = np.zeros(num_replicates)
        num_edges1 = np.zeros(num_replicates)
        num_nodes1 = np.zeros(num_replicates)
        reps = msprime.simulate(
            n, recombination_rate=recombination_rate, num_replicates=num_replicates
        )
        for j, ts in enumerate(reps):
            T1[j] = np.max(ts.tables.nodes.time)
            num_trees1[j] = ts.num_trees
            num_nodes1[j] = ts.num_nodes
            num_edges1[j] = ts.num_edges

        logging.debug(
            "original\tmean trees ={}\
             \tmean nodes = {}\
             \tmean edges =  {}".format(
                np.mean(num_trees1), np.mean(num_nodes1), np.mean(num_edges1)
            )
        )

        for t in [0.5, 1.0, 1.5, 5.0]:
            T2 = np.zeros(num_replicates)
            num_trees2 = np.zeros(num_replicates)
            num_nodes2 = np.zeros(num_replicates)
            num_edges2 = np.zeros(num_replicates)
            reps = msprime.simulate(
                n,
                num_replicates=num_replicates,
                recombination_rate=recombination_rate,
                end_time=t,
            )
            for j, ts in enumerate(reps):
                final_ts = msprime.simulate(
                    from_ts=ts,
                    recombination_rate=recombination_rate,
                    start_time=np.max(ts.tables.nodes.time),
                )
                assert max(t.num_roots for t in final_ts.trees()) == 1
                final_ts = final_ts.simplify()
                T2[j] = np.max(final_ts.tables.nodes.time)
                num_trees2[j] = final_ts.num_trees
                num_nodes2[j] = final_ts.num_nodes
                num_edges2[j] = final_ts.num_edges
            logging.debug(
                "t = {}\
                \tmean trees = {}\
                \tmean nodes = {}\
                \tmean edges = {}".format(
                    t, np.mean(num_trees2), np.mean(num_nodes2), np.mean(num_edges2)
                )
            )

            sm.graphics.qqplot(T1)
            sm.qqplot_2samples(T1, T2, line="45")
            filename = os.path.join(basedir, f"T_mrca_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_trees1)
            sm.qqplot_2samples(num_trees1, num_trees2, line="45")
            filename = os.path.join(basedir, f"num_trees_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_edges1)
            sm.qqplot_2samples(num_edges1, num_edges2, line="45")
            filename = os.path.join(basedir, f"num_edges_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_nodes1)
            sm.qqplot_2samples(num_nodes1, num_nodes2, line="45")
            filename = os.path.join(basedir, f"num_nodes_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

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
            msprime.PopulationConfiguration(),
        ]
        migration_matrix = [[0, 1], [1, 0]]
        demographic_events = [
            msprime.SimpleBottleneck(time=5.1, population=0, proportion=0.4),
            msprime.SimpleBottleneck(time=10.1, population=1, proportion=0.4),
            msprime.SimpleBottleneck(time=15.1, population=1, proportion=0.4),
            msprime.SimpleBottleneck(time=25.1, population=0, proportion=0.4),
        ]

        basedir = make_test_dir("simulate_from_demography")

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
            recombination_rate=recombination_rate,
        )
        logging.debug("t\ttrees\tnodes\tedges\tca\tre\tmig")
        for j in range(num_replicates):
            sim.run()
            ts = sim.get_tree_sequence()
            num_ca_events1[j] = sim.num_common_ancestor_events
            num_re_events1[j] = sim.num_recombination_events
            num_mig_events1[j] = sum(
                [r for row in sim.num_migration_events for r in row]
            )
            T1[j] = np.max(ts.tables.nodes.time)
            num_trees1[j] = ts.num_trees
            num_nodes1[j] = ts.num_nodes
            num_edges1[j] = ts.num_edges
            sim.reset()

        logging.debug(
            "{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                -1,
                np.mean(num_trees1),
                np.mean(num_nodes1),
                np.mean(num_edges1),
                np.mean(num_ca_events1),
                np.mean(num_re_events1),
                np.mean(num_mig_events1),
            )
        )

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
                recombination_rate=recombination_rate,
            )
            for j in range(num_replicates):
                sim.run()
                ts = sim.get_tree_sequence()
                num_ca_events2[j] = sim.num_common_ancestor_events
                num_re_events2[j] = sim.num_recombination_events
                num_mig_events2[j] = sum(
                    [r for row in sim.num_migration_events for r in row]
                )
                sim.reset()

                max_time = max(node.time for node in ts.nodes())
                sim2 = msprime.simulator_factory(
                    from_ts=ts,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=[
                        e for e in demographic_events if e.time > max_time
                    ],
                    recombination_rate=recombination_rate,
                )
                sim2.run()

                num_ca_events2[j] += sim2.num_common_ancestor_events
                num_re_events2[j] += sim2.num_recombination_events
                num_mig_events2[j] += sum(
                    [r for row in sim2.num_migration_events for r in row]
                )

                final_ts = sim2.get_tree_sequence().simplify()
                T2[j] = np.max(final_ts.tables.nodes.time)
                num_trees2[j] = final_ts.num_trees
                num_nodes2[j] = final_ts.num_nodes
                num_edges2[j] = final_ts.num_edges
                sim.reset()

            logging.debug(
                "{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                    t,
                    np.mean(num_trees2),
                    np.mean(num_nodes2),
                    np.mean(num_edges2),
                    np.mean(num_ca_events2),
                    np.mean(num_re_events2),
                    np.mean(num_mig_events2),
                )
            )

            sm.graphics.qqplot(T1)
            sm.qqplot_2samples(T1, T2, line="45")
            filename = os.path.join(basedir, f"T_mrca_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_trees1)
            sm.qqplot_2samples(num_trees1, num_trees2, line="45")
            filename = os.path.join(basedir, f"num_trees_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_edges1)
            sm.qqplot_2samples(num_edges1, num_edges2, line="45")
            filename = os.path.join(basedir, f"num_edges_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_nodes1)
            sm.qqplot_2samples(num_nodes1, num_nodes2, line="45")
            filename = os.path.join(basedir, f"num_nodes_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_ca_events1)
            sm.qqplot_2samples(num_ca_events1, num_ca_events2, line="45")
            filename = os.path.join(basedir, f"num_ca_events_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_re_events1)
            sm.qqplot_2samples(num_re_events1, num_re_events2, line="45")
            filename = os.path.join(basedir, f"num_re_events_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_mig_events1)
            sm.qqplot_2samples(num_mig_events1, num_mig_events2, line="45")
            filename = os.path.join(basedir, f"num_mig_events_t={t}.png")
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

    def run_simulate_from_benchmark(self):
        # A quick benchmark to show this running on a large example
        L = 50 * 10 ** 6
        seed = 3
        for n in [10 ** 3, 10 ** 4, 10 ** 5]:
            logging.debug("====================")
            logging.debug(f"n = {n}")
            logging.debug("====================")
            before = time.perf_counter()
            ts = msprime.simulate(
                n, recombination_rate=1e-8, Ne=10 ** 4, length=L, random_seed=seed
            )
            duration = time.perf_counter() - before

            logging.debug(f"Full sim required {duration:.2f} sec")

            before = time.perf_counter()
            t = ts.tables.nodes.time[-1] / 100
            ts = msprime.simulate(
                n,
                recombination_rate=1e-8,
                Ne=10 ** 4,
                length=L,
                random_seed=seed,
                end_time=t,
            )
            duration = time.perf_counter() - before
            logging.debug(f"Initial sim required {duration:.2f} sec")
            roots = np.array([tree.num_roots for tree in ts.trees()])
            logging.debug("\t", roots.shape[0], "trees, mean roots = ", np.mean(roots))
            before = time.perf_counter()
            msprime.simulate(
                from_ts=ts,
                recombination_rate=1e-8,
                Ne=10 ** 4,
                length=L,
                random_seed=seed,
            )
            duration = time.perf_counter() - before
            logging.debug(f"Final sim required {duration:.2f} sec")


@attr.s
class SimTest1(SimTest):
    """
    Check that the distributions are identitical when we run simulate_from
    at various time points.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_simulate_from_single_locus()


@attr.s
class SimTest2(SimTest):
    """
    Check that the distributions are identitical when we run simulate_from
    at various time points.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_simulate_from_multi_locus()


@attr.s
class SimTest3(SimTest):
    """
    Check that the distributions are identitical when we run simulate_from
    at various time points.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_simulate_from_recombination()


@attr.s
class SimTest4(SimTest):
    """
    Check that the distributions are identitical when we run simulate_from
    at various time points.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self.run_simulate_from_demography()


@attr.s
class SimTest5(SimTest):
    """
    Check that the distributions are identitical when we run simulate_from
    at various time points.
    """

    def run(self, output_dir):
        self._output_dir = output_dir
        self.run_simulate_from_benchmark()
        logging.info(f"running {self.group} {self.name}")


def register_sim_tests(runner):
    group = "sim"
    tests = [
        SimTest1(name="single locus", group=group),
        SimTest2(name="multi locus", group=group),
        SimTest3(name="recombination", group=group),
        SimTest4(name="demography", group=group),
        SimTest5(name="benchmark", group=group),
    ]
    for test in tests:
        runner.register(test)


def random_instance(
    key, group, num_populations=1, num_replicates=1000, num_demographic_events=0
):
    m = random.randint(1, 1000)
    r = random.uniform(0.01, 0.1) * m
    theta = random.uniform(1, 100)
    N = num_populations
    sample_sizes = [random.randint(2, 10) for _ in range(N)]
    migration_matrix = [random.random() * (j % (N + 1) != 0) for j in range(N ** 2)]
    structure = ""
    if num_populations > 1:
        structure = "-I {} {} -ma {}".format(
            num_populations,
            " ".join(str(s) for s in sample_sizes),
            " ".join(str(r) for r in migration_matrix),
        )
    cmd = "{} {} -t {} -r {} {} {}".format(
        sum(sample_sizes), num_replicates, theta, r, m, structure
    )

    if N > 1:
        # Add some migration matrix changes
        t = 0
        for j in range(1, 6):
            t += 0.125
            u = random.random()
            if u < 0.33:
                cmd += f" -eM {t} {random.random()}"
            elif u < 0.66:
                j = random.randint(1, N)
                k = j
                while k == j:
                    k = random.randint(1, N)
                r = random.random()
                cmd += f" -em {t} {j}"
            else:
                migration_matrix = [
                    random.random() * (j % (N + 1) != 0) for j in range(N ** 2)
                ]
                cmd += " -ema {} {} {}".format(
                    t, N, " ".join(str(r) for r in migration_matrix)
                )

    # Set some initial growth rates, etc.
    if N == 1:
        if random.random() < 0.5:
            cmd += f" -G {random.random()}"
        else:
            cmd += f" -eN 0 {random.random()}"
    # Add some demographic events
    t = 0
    for _ in range(num_demographic_events):
        t += 0.125
        if random.random() < 0.5:
            cmd += f" -eG {t} {random.random()}"
        else:
            cmd += f" -eN {t} {random.random()}"

    return MsTest1(key, group, cmd)


def register_msrandom_tests(runner):
    def register(name, **kwargs):
        runner.register(random_instance(name, "random", **kwargs))

    register("random1")
    register("random2", num_replicates=10 ** 4, num_demographic_events=10)


@attr.s
class MutationTest(Test):
    def _transition_matrix_chi_sq(self, transitions, transition_matrix):
        tm_chisq = []
        for row, p in zip(transitions, transition_matrix):
            not_zeros = p > 0
            if sum(not_zeros) > 1:
                chisq = scipy.stats.chisquare(row[not_zeros], p[not_zeros])
                tm_chisq.append(chisq.statistic)
            else:
                tm_chisq.append(None)

        return tm_chisq

    def _transitions(self, sequences, ts, alleles, mutation_rate, Q):
        num_alleles = len(alleles)
        transitions = np.zeros((num_alleles, num_alleles), dtype=int)
        expected = np.zeros((num_alleles, num_alleles))

        for edge in ts.edges():
            for idx in range(int(ts.sequence_length)):
                p = sequences[edge.parent][idx]
                c = sequences[edge.child][idx]
                transitions[alleles.index(p), alleles.index(c)] += 1
                j = alleles.index(p)
                expected[j, :] += _matrix_exponential(
                    ts.first().branch_length(edge.child) * mutation_rate * Q
                )[j, :]

        return (transitions, expected)

    def get_allele_counts(self, ts):
        if ts.num_sites == 0:
            df_ts = allel.HaplotypeArray(np.zeros((2, ts.num_samples), dtype=int))
        else:
            df_ts = allel.HaplotypeArray(ts.genotype_matrix())
        return df_ts.count_alleles()

    def get_transition_stats(self, ts, alleles, mutation_rate, Q):
        num_alleles = len(alleles)
        observed_transitions_ts = np.zeros((num_alleles, num_alleles))
        expected_ts = np.zeros((num_alleles, num_alleles))

        corr = ts.sequence_length / ts.num_sites
        # -> for this method to perform optimally, corr==1
        # at least one  mutation on each site
        assert ts.num_trees == 1
        tree = ts.first()
        for v in ts.variants(samples=range(ts.num_nodes), impute_missing_data=True):
            for n in tree.nodes():
                pn = tree.parent(n)
                if pn != tskit.NULL:
                    pa = v.alleles[v.genotypes[pn]]
                else:
                    pa = v.site.ancestral_state
                da = v.alleles[v.genotypes[n]]
                observed_transitions_ts[alleles.index(pa), alleles.index(da)] += 1
                j = alleles.index(pa)
                expected_ts[j, :] += _matrix_exponential(
                    tree.branch_length(n) * mutation_rate * corr * Q
                )[j, :]
        return observed_transitions_ts, expected_ts

    def plot_stats(self, df_test, df_msprime, alleles, test_prog, model):
        test_key = f"{test_prog}-{model}"
        # plot results
        for name in ["pi", "root_distribution"]:
            sg_results = sm.ProbPlot(df_test[name].dropna())
            ts_results = sm.ProbPlot(df_msprime[name].dropna())
            sm.qqplot_2samples(
                sg_results,
                ts_results,
                ylabel=f"quantiles {test_prog}",
                xlabel="quantiles msprime",
                line="45",
            )
            outfile = self._build_filename(test_key, name)
            pyplot.savefig(outfile)
        pyplot.clf()
        if len(alleles) == 4:
            rows, columns = 2, 2
        else:
            rows, columns = 5, 4
        fig, axs = pyplot.subplots(rows, columns, figsize=(12, 12))
        for i, co in enumerate(itertools.product(range(rows), range(columns))):
            a = alleles[i]
            size = min(df_test[a].dropna().size, df_msprime[a].dropna().size)
            temp_test = sm.ProbPlot(df_test[a].dropna()[:size])
            temp_msprime = sm.ProbPlot(df_msprime[a].dropna()[:size])
            sm.qqplot_2samples(
                temp_test,
                temp_msprime,
                ylabel=f"quantiles {test_prog}",
                xlabel="quantiles msprime",
                line="45",
                ax=axs[co],
            )
            axs[co].set_title(a)
        outfile = self._build_filename(test_key, "alleles")
        pyplot.savefig(outfile)


class SeqGenTest(MutationTest):
    _seq_gen_executable = ["./data/seq-gen"]

    def _run_seq_gen(self, tree, args, model, alleles, num_sites, mutation_rate, Q):

        ts = tree.tree_sequence
        newick = tree.newick()
        cmd = self._seq_gen_executable + args
        num_sequences = 2 * ts.num_samples - 1
        with tempfile.TemporaryFile("w+") as in_file, tempfile.TemporaryFile(
            "w+"
        ) as out_file:
            in_file.write(newick)
            in_file.seek(0)
            subprocess.call(cmd, stdin=in_file, stdout=out_file)
            out_file.seek(0)
            sequences = {}
            # Skip the first line
            out_file.readline()
            for line, node in zip(out_file, ts.first().nodes()):
                sample_id, sequence = line.split()
                sequences[node] = sequence
                assert len(sequence) == ts.sequence_length
        assert len(sequences) == num_sequences

        num_alleles = len(alleles)

        ancestral_sequence = sequences[len(sequences) - 1]
        observed_ancestral_sg = np.zeros((num_alleles,))
        for idx in np.random.choice(int(ts.sequence_length), num_sites, replace=False):
            b = ancestral_sequence[idx]
            observed_ancestral_sg[alleles.index(b)] += 1

        def replace_variants(variants):
            u = np.unique(variants)
            repl = [i for i in range(len(u))]
            return np.array([dict(zip(u, repl))[i] for i in variants])

        ord_sequences = {
            key: [ord(element) % 32 for element in value]
            for key, value in sequences.items()
        }
        transitions_sg, expected = self._transitions(
            sequences, ts, alleles, mutation_rate, Q
        )

        sg_sequences = np.transpose(
            np.array([ord_sequences[key] for key in range(ts.num_samples)])
        )
        sg_reduced = np.apply_along_axis(replace_variants, 1, sg_sequences)
        sg_genotypes = allel.HaplotypeArray(sg_reduced)
        sg_counts = sg_genotypes.count_alleles()

        return (sg_counts, transitions_sg, observed_ancestral_sg, expected)

    def _run_seq_gen_msprime_stats(self, model, length=20, num_samples=10):
        """
        Runs a comparison between mutations generated by msprime and seq_gen
        for the specified model and returns a tuple of data frames ready
        for plotting.
        """
        model_dict = {
            "JC69": {"model_id": msprime.JukesCantor(), "par": ["-m", "HKY"]},
            "HKY": {
                "model_id": msprime.HKY(
                    kappa=1.5, equilibrium_frequencies=[0.2, 0.3, 0.1, 0.4]
                ),
                "par": ["-m", "HKY", "-t", "0.75", "-f", "0.2,0.3,0.1,0.4"],
            },
            "F84": {
                "model_id": msprime.F84(
                    kappa=1.0, equilibrium_frequencies=[0.3, 0.25, 0.2, 0.25]
                ),
                "par": ["-m", "F84", "-t", "0.5", "-f", "0.3,0.25,0.2,0.25"],
            },
            "GTR": {
                "model_id": msprime.GTR(
                    relative_rates=[0.4, 0.1, 0.4, 0.2, 0.4, 0.4],
                    equilibrium_frequencies=[0.3, 0.2, 0.3, 0.2],
                ),
                "par": [
                    "-m",
                    "GTR",
                    "-r",
                    "0.4,0.1,0.4,0.2,0.4,0.4",
                    "-f",
                    "0.3,0.2,0.3,0.2",
                ],
            },
            "PAM": {"model_id": msprime.PAM(), "par": ["-m", "PAM"]},
            "BLOSUM62": {"model_id": msprime.BLOSUM62(), "par": ["-m", "BLOSUM"]},
        }

        num_replicates = 250
        sg_results = collections.defaultdict(list)
        ts_results = collections.defaultdict(list)
        pos = [i for i in range(1, length + 1)]
        transition_matrix = model_dict[model]["model_id"].transition_matrix
        root_distribution = model_dict[model]["model_id"].root_distribution
        alleles = model_dict[model]["model_id"].alleles
        num_alleles = len(alleles)
        mutation_rate = 1e-4 if num_alleles == 4 else 1.5e-3
        Q = transition_matrix.copy()
        Q -= np.eye(num_alleles)
        mut_rate_seq_gen = np.sum(-Q.diagonal() * root_distribution) * mutation_rate
        args = ["-q", "-s", str(mut_rate_seq_gen), "-l", str(length), "-wa"]
        args += model_dict[model]["par"]

        Ne = 1e4

        for _ in tqdm.tqdm(range(num_replicates)):
            ts = msprime.simulate(num_samples, Ne=Ne, length=length)
            ts_mutated = msprime.mutate(
                ts,
                rate=mutation_rate,
                model=model_dict[model]["model_id"],
                discrete=True,
            )
            num_sites = ts_mutated.num_sites
            t = ts_mutated.first()
            t_span = np.ceil(t.interval[1] - np.ceil(t.interval[0]))
            # expected number of ancestral alleles for sites
            expected_ancestral_states_ts = np.zeros(num_alleles)
            change_probs = transition_matrix.sum(axis=1) - np.diag(transition_matrix)
            expected_ancestral_states_ts += (
                root_distribution
                * t_span
                * (1 - np.exp(-mutation_rate * t.total_branch_length * change_probs))
            )

            # observed number of ancestral alleles
            obs_ancestral_states_ts = np.zeros((num_alleles,))
            for site in ts_mutated.sites():
                aa = site.ancestral_state
                obs_ancestral_states_ts[alleles.index(aa)] += 1

            # expected and observed number of transitions ts
            # root distribution == equilibrium freqs for these tests,
            # as is the case in seq-gen
            observed_transitions_ts, expected_ts = self.get_transition_stats(
                ts_mutated, alleles, mutation_rate, Q
            )

            # run Seq-gen and calculate statistics
            (
                c_sg,
                observed_transitions_sg,
                observed_ancestral_sg,
                expected_sg,
            ) = self._run_seq_gen(
                t,
                args,
                model_dict[model]["model_id"],
                alleles,
                num_sites,
                mutation_rate,
                Q,
            )

            c_ts = self.get_allele_counts(ts_mutated)
            # Compute pi
            pi_sg = allel.sequence_diversity(pos, c_sg)
            sg_results["pi"].append(pi_sg)
            pi_ts = allel.sequence_diversity(pos, c_ts)
            ts_results["pi"].append(pi_ts)

            # Compute chisquare stats.
            tm_chisq_sg = self._transition_matrix_chi_sq(
                observed_transitions_sg, expected_sg
            )
            # in Seq-Gen the ancestral sequence is determined first
            expected_num_ancestral_states_sg = root_distribution * num_sites
            root_chisq_sg = scipy.stats.chisquare(
                observed_ancestral_sg, expected_num_ancestral_states_sg
            ).statistic

            tm_chisq_ts = self._transition_matrix_chi_sq(
                observed_transitions_ts, expected_ts
            )
            root_chisq_ts = scipy.stats.chisquare(
                obs_ancestral_states_ts, expected_ancestral_states_ts
            ).statistic
            ts_results["root_distribution"].append(root_chisq_ts)
            sg_results["root_distribution"].append(root_chisq_sg)

            for idx, a in enumerate(alleles):
                sg_results[a].append(tm_chisq_sg[idx])
                ts_results[a].append(tm_chisq_ts[idx])

        df_sg = pd.DataFrame.from_dict(sg_results)
        df_ts = pd.DataFrame.from_dict(ts_results)
        return df_sg, df_ts, alleles

    def _run_seq_gen_msprime_comparison(self, model, length=20, num_samples=10):
        df_sg, df_ts, alleles = self._run_seq_gen_msprime_stats(
            model, length, num_samples
        )
        self.plot_stats(df_sg, df_ts, alleles, "seqgen", model)


@attr.s
class SeqGenTest1(SeqGenTest):
    model = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self._run_seq_gen_msprime_comparison(self.model)


def register_seqgen_tests(runner):
    def register(model):
        group = "seqgen"
        name = f"{group}_{model}"
        runner.register(SeqGenTest1(name=name, group=group, model=model))

    register("JC69")
    register("HKY")
    register("F84")
    register("GTR")
    register("PAM")
    register("BLOSUM62")


@attr.s
class PyvolveTest(MutationTest):
    def _run_pyvolve(
        self, tree, py_model, model, alleles, num_sites, mutation_rate, ts_mutrate, Q
    ):
        ts = tree.tree_sequence
        seq_length = int(ts.sequence_length)
        node_labels = {u: str(u) for u in ts.samples()}
        newick = tree.newick(node_labels=node_labels)
        pyvolve_tree = pyvolve.read_tree(tree=newick, scale_tree=mutation_rate)
        pyvolve_model = pyvolve.Partition(models=py_model, size=seq_length)
        sim = pyvolve.Evolver(tree=pyvolve_tree, partitions=pyvolve_model)
        sim(ratefile=None, infofile=None, seqfile=None)
        seqs = sim.get_sequences(anc=True)  # seq-dict is sorted in pre-order
        sequences = {}
        for key, node in zip(seqs.keys(), ts.first().nodes()):
            sequences[node] = seqs[key]
            assert len(seqs[key]) == ts.sequence_length
        assert len(sequences) == 2 * ts.num_samples - 1

        num_alleles = len(alleles)
        ancestral_sequence = sequences[len(sequences) - 1]
        roots_d_py = np.zeros((num_alleles,))
        for idx in np.random.choice(int(ts.sequence_length), num_sites, replace=False):
            b = ancestral_sequence[idx]
            roots_d_py[alleles.index(b)] += 1

        def replace_variants(variants):
            u = np.unique(variants)
            repl = [i for i in range(len(u))]
            return np.array([dict(zip(u, repl))[i] for i in variants])

        ord_sequences = {
            key: [ord(element) % 32 for element in value]
            for key, value in sequences.items()
        }
        transitions_py, expected = self._transitions(
            sequences, ts, alleles, ts_mutrate, Q
        )
        py_sequences = np.transpose(
            np.array([ord_sequences[key] for key in range(ts.num_samples)])
        )
        py_reduced = np.apply_along_axis(replace_variants, 1, py_sequences)
        py_genotypes = allel.HaplotypeArray(py_reduced)
        py_counts = py_genotypes.count_alleles()

        return (py_counts, transitions_py, roots_d_py, expected)

    def _run_pyvolve_stats(self, model, length=20, num_samples=10):

        model_dict = {
            "JC69": {
                "model_id": msprime.JukesCantor(),
                "pyvolve_model": pyvolve.Model("nucleotide"),
            },
            "HKY": {
                "model_id": msprime.HKY(
                    kappa=1.5, equilibrium_frequencies=[0.2, 0.3, 0.1, 0.4]
                ),
                "pyvolve_model": pyvolve.Model(
                    "nucleotide", {"kappa": 1.5, "state_freqs": [0.2, 0.3, 0.1, 0.4]}
                ),
            },
            "PAM": {
                "model_id": msprime.PAM(),
                "pyvolve_model": pyvolve.Model("DAYHOFFDCMUT"),
            },
            "BLOSUM62": {
                "model_id": msprime.BLOSUM62(),
                "pyvolve_model": pyvolve.Model("BLOSUM62"),
            },
        }

        num_replicates = 250
        py_results = collections.defaultdict(list)
        ts_results = collections.defaultdict(list)
        pos = [i for i in range(1, length + 1)]
        alleles = model_dict[model]["model_id"].alleles
        num_alleles = len(alleles)
        mutation_rate = 1e-4 if num_alleles == 4 else 1.5e-3
        transition_matrix = model_dict[model]["model_id"].transition_matrix
        root_distribution = model_dict[model]["model_id"].root_distribution
        Q = transition_matrix.copy()
        Q -= np.eye(num_alleles)
        mut_rate_pyvolve = np.sum(-Q.diagonal() * root_distribution) * mutation_rate

        for _ in tqdm.tqdm(range(num_replicates)):
            ts = msprime.simulate(num_samples, Ne=1e4, length=length)
            ts_mutated = msprime.mutate(
                ts,
                rate=mutation_rate,
                model=model_dict[model]["model_id"],
                discrete=True,
            )

            num_sites = ts_mutated.num_sites

            t = ts_mutated.first()
            t_span = np.ceil(t.interval[1] - np.ceil(t.interval[0]))
            # expected number of ancestral alleles for sites
            expected_ancestral_states_ts = np.zeros(num_alleles)
            change_probs = transition_matrix.sum(axis=1) - np.diag(transition_matrix)
            expected_ancestral_states_ts += (
                root_distribution
                * t_span
                * (1 - np.exp(-mutation_rate * t.total_branch_length * change_probs))
            )

            # observed number of ancestral alleles
            obs_ancestral_states_ts = np.zeros((num_alleles,))
            for site in ts_mutated.sites():
                aa = site.ancestral_state
                obs_ancestral_states_ts[alleles.index(aa)] += 1
            observed_transitions_ts, expected = self.get_transition_stats(
                ts_mutated, alleles, mutation_rate, Q
            )

            # run pyvolve and calculate statistics
            (
                c_py,
                observed_transitions_py,
                observed_ancestral_py,
                expected_py,
            ) = self._run_pyvolve(
                t,
                model_dict[model]["pyvolve_model"],
                model_dict[model]["model_id"],
                alleles,
                num_sites,
                mut_rate_pyvolve,
                mutation_rate,
                Q,
            )
            pi_py = allel.sequence_diversity(pos, c_py)

            tm_chisq_py = self._transition_matrix_chi_sq(
                observed_transitions_py, expected_py
            )

            expected_num_ancestral_states_py = root_distribution * num_sites
            root_chisq_py = scipy.stats.chisquare(
                observed_ancestral_py, expected_num_ancestral_states_py
            ).statistic

            tm_chisq_ts = self._transition_matrix_chi_sq(
                observed_transitions_ts, expected
            )

            root_chisq_ts = scipy.stats.chisquare(
                obs_ancestral_states_ts, expected_ancestral_states_ts
            ).statistic

            c_ts = self.get_allele_counts(ts_mutated)
            pi_ts = allel.sequence_diversity(pos, c_ts)

            ts_results["pi"].append(pi_ts)
            ts_results["root_distribution"].append(root_chisq_ts)
            py_results["pi"].append(pi_py)
            py_results["root_distribution"].append(root_chisq_py)
            for idx, a in enumerate(alleles):
                ts_results[a].append(tm_chisq_ts[idx])
                py_results[a].append(tm_chisq_py[idx])

        df_py = pd.DataFrame.from_dict(py_results)
        df_ts = pd.DataFrame.from_dict(ts_results)
        return df_py, df_ts, alleles

    def _run_pyvolve_comparison(self, model, length=20, num_samples=10):
        df_py, df_ts, alleles = self._run_pyvolve_stats(model, length, num_samples)
        self.plot_stats(df_py, df_ts, alleles, "pyvolve", model)


@attr.s
class PyvolveTest1(PyvolveTest):
    model = attr.ib(type=str)

    def run(self, output_dir):
        self._output_dir = output_dir
        logging.info(f"running {self.group} {self.name}")
        self._run_pyvolve_comparison(self.model)


def register_pyvolve_tests(runner):
    def register(model):
        group = "pyvolve"
        name = f"{group}_{model}"
        runner.register(PyvolveTest1(name=name, group=group, model=model))

    register("JC69")
    register("HKY")
    register("PAM")
    register("BLOSUM62")


@attr.s
class TestRunner:
    """
    Class responsible for registering all known tests and running
    them.
    """

    tests = attr.ib(init=False, default=[])
    groups = attr.ib(init=False, default=set())

    def register(self, test):
        self.tests.append(test)
        self.groups.add(test.group)

    def run(self, output_dir, threads=1, names=None, group=None):
        if names is not None and group is not None:
            raise ValueError("Cannot specify test names and group at the same time")

        if names is not None:
            tests = [test for test in self.tests if test.name in names]
        elif group is not None:
            tests = [test for test in self.tests if test.group == group]
        else:
            tests = self.tests

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(test.run, output_dir) for test in tests]

            assert len(futures) == len(tests)

            pbar = tqdm.tqdm(total=len(futures), desc="Running tests", leave=True)
            for future in concurrent.futures.as_completed(futures):
                pbar.update(1)
                logging.debug(f"{future} done")


def setup_logging(args):
    log_level = "WARN"
    if args.verbose == 1:
        log_level = "INFO"
    elif args.verbose >= 2:
        log_level = "DEBUG"

    daiquiri.setup(level=log_level)
    msprime_logger = daiquiri.getLogger("msprime")
    msprime_logger.setLevel("WARN")


def add_simulator_arguments(parser, groups):
    parser.add_argument(
        "--extended",
        action="store_true",
        help="Run extended tests in dtwf_vs_coalescent",
    )
    parser.add_argument(
        "--group",
        "-g",
        default=None,
        choices=groups,
        help="Run all tests for specified group",
    )
    parser.add_argument("tests", nargs="*", help="Run specific tests")
    parser.add_argument(
        "--num-threads", "-t", type=int, default=1, help="Specify number of threads"
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Set the log level. 0:WARN, 1:INFO, 2:DEBUG",
    )


def run_tests(args, runner):
    output_dir = "tmp__NOBACKUP__"
    if args.extended:
        register_dtwfvscoalescent_tests(runner, add_random_dtwf=True)

    if len(args.tests) > 0:
        runner.run(output_dir, threads=args.num_threads, names=args.tests)
    elif args.group is not None:
        runner.run(output_dir, threads=args.num_threads, group=args.group)
    else:
        runner.run(output_dir, threads=args.num_threads)


def main():
    runner = TestRunner()

    register_ms_tests(runner)
    register_discoal_tests(runner)
    register_contdiscrete_tests(runner)
    register_dtwfvscoalescent_tests(runner, add_random_dtwf=False)
    register_xi_tests(runner)
    register_hudson_tests(runner)
    register_argrecord_tests(runner)
    register_msrandom_tests(runner)
    register_sim_tests(runner)
    register_analytical_tests(runner)
    register_smc_tests(runner)
    register_seqgen_tests(runner)
    register_pyvolve_tests(runner)

    parser = argparse.ArgumentParser()
    add_simulator_arguments(parser, runner.groups)
    args = parser.parse_args()
    setup_logging(args)
    run_tests(args, runner)


if __name__ == "__main__":
    main()
