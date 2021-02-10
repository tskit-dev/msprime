"""
Script to automate verification of msprime against known statistical
results and benchmark programs such as ms and Seq-Gen.

Tests are structured in a similar way to Python unittests. Tests
are organised into classes of similar tests. Ideally, each test
in the class is a simple call to a general method with
different parameters (this is called ``_run``, by convention).
Tests must be *independent* and not depend on any shared
state within the test class, other than the ``self.output_dir``
variable which is guaranteed to be set when the method is called.

The output directory is <output-dir>/<class name>/<test name>.
Each test should output one or more diagnostic plots, which have
a clear interpretation as "correct" or "incorrect". QQ-plots
are preferred, where possible. Numerical results can also be
output by using ``logging.debug()``, where appropriate; to
view these, append ``--debug`` to the comand line running
your tests.

Test classes must be a subclass of the ``Test`` class defined
in this module.

To run the tests, first get some help from the CLI:

    python3 verification.py --help

This will output some basic help on the tests. Use

    python3 verification.py --list

to show all the available tests.

If you run without any arguments, this will run all the tests
sequentially. The progress bar and output behaviour can be
controlled using command line parameters, and running over
multiple processes is possible.

If you wish to run a specific tests, you can provide the
test names as positional arguments, i.e.,

    python3 verification.py test_msdoc_outgroup_sequence test_msdoc_recomb_ex

will just run these two specific tests.

Using the ``-c`` option allows you to run all tests in a
given class.

Gotchas:
- Any test superclasses must be abstract. That is, you cannot
  inherit from a test class that contains any tests.
- Test method names must be unique across *all* classes.

"""
import argparse
import ast
import collections
import concurrent.futures
import inspect
import itertools
import json
import logging
import math
import pathlib
import pickle
import random
import subprocess
import sys
import tempfile
import warnings

import allel
import attr
import daiquiri
import dendropy
import matplotlib
import numpy as np
import pandas as pd
import pyslim
import pyvolve
import scipy.special
import scipy.stats
import seaborn as sns
import tqdm
import tskit
from matplotlib import pyplot

import msprime
import msprime.cli as cli
from msprime.demography import _matrix_exponential

# Force matplotlib to not use any Xwindows backend.
# Note this must be done before importing statsmodels.
matplotlib.use("Agg")
import statsmodels.api as sm  # noqa: E402


_mspms_executable = [sys.executable, "mspms_dev.py"]
_slim_executable = ["./data/slim"]
_ms_executable = ["./data/ms"]
_discoal_executable = ["./data/discoal"]
_scrm_executable = ["./data/scrm"]
_msms_executable = ["java", "-Xmx1G", "-jar", "data/msms.jar"]


def flatten(li):
    return [x for sublist in li for x in sublist]


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


def write_sweep_slim_script(outfile, format_dict):
    slim_str = """
        initialize() {{
        initializeTreeSeq();
        initializeMutationRate(0);
        initializeMutationType('m1', 0.5, 'f', 0.0);
        initializeMutationType('m2', 0.5, 'f', {s});
        initializeGenomicElementType('g1', m1, 1.0);
        initializeGenomicElement(g1, 0, {NUMLOCI});
        initializeRecombinationRate({r});
        }}
        s1 200000 late() {{
                sim.treeSeqOutput('{OUTFILE}');
                    sim.simulationFinished();
        }}

        1 {{
            // save this run's identifier, used to save and restore
            defineConstant("simID", getSeed());
            sim.addSubpop("p1", {POPSIZE});
            sim.setValue("flag",0);
        }}

        2 late() {{
            // save the state of the simulation
            sim.treeSeqOutput("/tmp/slim_" + simID + ".trees");
            target = sample(p1.genomes, 1);
            target.addNewDrawnMutation(m2, {SWEEPPOS});
        }}
        2:2000 late() {{
            if (sim.countOfMutationsOfType(m2) == 0)
            {{
                fixed = (sum(sim.substitutions.mutationType == m2) == 1);
                if (fixed){{
                    sim.setValue("flag", sim.getValue("flag") + 1);
                    }}
                if (fixed)
                {{
                    if (sim.getValue("flag") == 1){{
                        sim.rescheduleScriptBlock(s1,
                        start=sim.generation+{TAU}, end=sim.generation+{TAU});
                    }}
                }}
                else
                {{
                    sim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");
                    setSeed(rdunif(1, 0, asInteger(2^62) - 1));
                    target = sample(p1.genomes, 1);
                    target.addNewDrawnMutation(m2, {SWEEPPOS});
                }}
            }}
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


def plot_qq(v1, v2):
    sm.graphics.qqplot(v1)
    sm.qqplot_2samples(v1, v2, line="45")


def plot_stat_hist(v1, v2, v1_name, v2_name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sns.kdeplot(v1, color="b", shade=True, label=v1_name, legend=False)
        sns.kdeplot(v2, color="r", shade=True, label=v2_name, legend=False)
        pyplot.legend(loc="upper right")


def plot_breakpoints_hist(v1, v2, v1_name, v2_name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sns.kdeplot(v1, color="b", label=v1_name, shade=True, legend=False)
        sns.kdeplot(v2, color="r", label=v2_name, shade=True, legend=False)
        pyplot.legend(loc="upper right")


def all_breakpoints_in_replicates(replicates):
    return [right for intervals in replicates for left, right in intervals]


@attr.s
class Test:
    """
    The superclass of all tests. The only attribute defined is the output
    directory for the test, which is guaranteed to exist when the
    test method is called.
    """

    output_dir = attr.ib(type=str, default=None)

    def _run_sample_stats(self, args):
        logging.debug(f"{' '.join(args)}")
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
            df = pd.read_csv(f, sep="\t")
        return df

    def _build_filename(self, *args):
        return self.output_dir / "_".join(args[1:])

    def _plot_stats(self, stats_type, df1, df2, df1_name, df2_name):
        assert set(df1.columns.values) == set(df2.columns.values)
        for stat in df1.columns.values:
            v1 = df1[stat]
            v2 = df2[stat]
            if stat == "breakpoints":
                plot_breakpoints_hist(flatten(v1), flatten(v2), df1_name, df2_name)
                pyplot.xlabel("genome")
                f = self._build_filename(stats_type, stat)
                pyplot.savefig(f, dpi=72)
            else:
                plot_qq(v1, v2)
                pyplot.xlabel(df1_name)
                pyplot.ylabel(df2_name)
                f = self._build_filename(stats_type, stat)
                pyplot.savefig(f, dpi=72)
                pyplot.close("all")
                plot_stat_hist(v1, v2, df1_name, df2_name)
                f = self._build_filename(stats_type, stat)
                f = str(f) + ".hist.png"
                pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def get_ms_seeds(self):
        max_seed = 2 ** 16
        seeds = [random.randint(1, max_seed) for j in range(3)]
        return ["-seed"] + list(map(str, seeds))

    def _run_msprime_mutation_stats(self, args):
        return self._run_sample_stats(
            _mspms_executable + args.split() + self.get_ms_seeds()
        )


class MsTest(Test):
    """
    Superclass of tests that perform comparisons with ms. Provides some
    infrastructure for common operations.
    """

    def _deserialize_breakpoints(self, df):
        breakpoints_strs = df["breakpoints"]
        breakpoints = [ast.literal_eval(literal) for literal in breakpoints_strs]
        df["breakpoints"] = breakpoints
        return df

    def _exec_coalescent_stats(self, executable, args, seeds=None):
        with tempfile.TemporaryFile() as f:
            argList = [executable] + args.split() + self.get_ms_seeds()
            logging.debug(f"{' '.join(argList)}")
            subprocess.call(argList, stdout=f)
            f.seek(0)
            df = pd.read_table(f)
        self._deserialize_breakpoints(df)
        return df

    def _run_ms_coalescent_stats(self, args):
        return self._exec_coalescent_stats("./data/ms_summary_stats", args)

    def _run_ms_mutation_stats(self, args):
        return self._run_sample_stats(
            _ms_executable + args.split() + self.get_ms_seeds()
        )

    def _run_mutation_stats(self, args):
        df_ms = self._run_ms_mutation_stats(args)
        df_msp = self._run_msprime_mutation_stats(args)
        self._plot_stats("mutation", df_ms, df_msp, "ms", "msp")

    def _run_mspms_coalescent_stats(self, args):
        logging.debug(f"mspms: {args}")
        runner = cli.get_mspms_runner(args.split())
        sim = runner.simulator
        num_populations = sim.num_populations
        replicates = runner.num_replicates
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

    def _run_coalescent_stats(self, args):
        df_msp = self._run_mspms_coalescent_stats(args)
        df_ms = self._run_ms_coalescent_stats(args)
        self._plot_stats("coalescent", df_msp, df_ms, "msp", "ms")

    # end of tests common to MS and random
    def _run_variable_recombination_coalescent_stats(self, args):
        df_msp = self._run_mspms_coalescent_stats(args)
        df_mshot = self._run_mshot_coalescent_stats(args)
        self._plot_stats("recomb map coalescent", df_msp, df_mshot, "msp", "msHOT")

    def _run_mshot_coalescent_stats(self, args):
        return self._exec_coalescent_stats("./data/msHOT_summary_stats", args)

    def _run(self, cmd):
        self._run_coalescent_stats(cmd)
        self._run_mutation_stats(cmd)


class MsDemography(MsTest):
    def test_size_change_1(self):
        self._run("10 10000 -t 2.0 -eN 0.1 2.0")

    def test_growth_rate_change_1(self):
        self._run("10 10000 -t 2.0 -eG 0.1 5.0")

    def test_growth_rate_change1(self):
        self._run("10 10000 -t 2.0 -eG 0.1 5.0")

    def test_growth_rate_2_pops1(self):
        self._run("10 10000 -t 2.0 -I 2 5 5 2.5 -G 5.0")

    def test_growth_rate_2_pops2(self):
        self._run("10 10000 -t 2.0 -I 2 5 5 2.5 -G 5.0 -g 1 0.1")

    def test_growth_rate_2_pops3(self):
        self._run("10 10000 -t 2.0 -I 2 5 5 2.5 -g 1 0.1")

    def test_growth_rate_2_pops4(self):
        self._run("10 10000 -t 2.0 -I 2 5 5 2.5 -eg 1.0 1 5.0")

    def test_pop_size_2_pops1(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 2.5 -n 1 0.1")

    def test_pop_size_2_pops2(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 2.5 -g 1 2 -n 1 0.1")

    def test_pop_size_2_pops3(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 2.5 -eN 0.5 3.5")

    def test_pop_size_2_pops4(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 2.5 -en 0.5 1 3.5")

    def test_migration_rate_2_pops1(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 0 -eM 3 5")

    def test_migration_matrix_2_pops1(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 -ma x 10 0 x")

    def test_migration_matrix_2_pops2(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 -m 1 2 10 -m 2 1 50")

    def test_migration_rate_change_2_pops1(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 -eM 5 10")

    def test_migration_matrix_entry_change_2_pops1(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 -em 0.5 2 1 10")

    def test_migration_matrix_change_2_pops1(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 -ema 10.0 2 x 10 0 x")

    def migration_matrix_change_2_pops2(self):
        cmd = """100 10000 -t 2.0 -I 2 50 50 -ema 1.0
          2 x 0.1 0 x -eN 1.1 0.001 -ema 10 2 x 0 10 x"""
        self._run(cmd)

    def test_population_split_2_pops1(self):
        self._run("100 10000 -t 2.0 -I 2 50 50 5.0 -ej 2.0 1 2")

    def test_population_split_4_pops1(self):
        self._run("100 10000 -t 2.0 -I 4 50 50 0 0 2.0 -ej 0.5 2 1")

    def test_population_split_4_pops2(self):
        self._run("100 10000 -t 2.0 -I 4 25 25 25 25 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1")

    def test_population_split_4_pops3(self):
        cmd = (
            "100 10000 -t 2.0 -I 4 25 25 25 25 -ej 1 2 1 "
            "-em 1.5 4 1 2 -ej 2 3 1 -ej 3 4 1"
        )
        self._run(cmd)

    def test_admixture_1_pop1(self):
        self._run("1000 1000 -t 2.0 -es 0.1 1 0.5 -em 0.1 1 2 1")

    def test_admixture_1_pop2(self):
        self._run("1000 1000 -t 2.0 -es 0.1 1 0.1 -em 0.1 1 2 1")

    def test_admixture_1_pop3(self):
        self._run("1000 1000 -t 2.0 -es 0.01 1 0.1 -em 0.1 2 1 1")

    def test_admixture_1_pop4(self):
        self._run("1000 1000 -t 2.0 -es 0.01 1 0.1 -es 0.1 2 0 -em 0.1 3 1 1")

    def test_admixture_1_pop5(self):
        self._run("1000 1000 -t 2.0 -es 0.01 1 0.1 -ej 1 2 1")

    def test_admixture_1_pop6(self):
        self._run("1000 1000 -t 2.0 -es 0.01 1 0.0 -eg 0.02 2 5.0 ")

    def test_admixture_1_pop7(self):
        self._run("1000 1000 -t 2.0 -es 0.01 1 0.0 -en 0.02 2 5.0 ")

    def test_admixture_2_pop1(self):
        self._run("1000 1000 -t 2.0 -I 2 500 500 1 -es 0.01 1 0.1 -ej 1 3 1")

    def test_admixture_2_pop2(self):
        self._run("1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -em 2.0 3 1 1")

    def test_admixture_2_pop3(self):
        self._run(
            "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 -G 5.0 " "-em 2.0 3 1 1"
        )

    def test_admixture_2_pop4(self):
        cmd = (
            "1000 1000 -t 2.0 -I 2 500 500 2 -es 0.01 1 0.75 "
            "-eg 0.02 1 5.0 -em 0.02 3 1 1"
        )
        self._run(cmd)


class MsGeneConversion(MsTest):
    def _run(self, cmd):
        # The mutation stats are a waste of time for GC, they tell us basically
        # nothing.
        self._run_coalescent_stats(cmd)

    def test_gene_conversion_c10_r0(self):
        self._run("100 10000 -t 5.0 -r 0 2501 -c 10 1")

    def test_gene_conversion_c100_tl1000_r0(self):
        self._run("100 10000 -t 5.0 -r 0 2501 -c 100 1000")

    def test_gene_conversion_c1000_tl_1(self):
        self._run("100 10000 -t 5.0 -r 0.01 2501 -c 1000 1")

    def test_gene_conversion_c1000_tl_1000(self):
        self._run("100 10000 -t 5.0 -r 0.01 2501 -c 1000 1000")

    def test_gene_conversion_c2_r10(self):
        self._run("100 10000 -t 5.0 -r 10 2501 -c 2 1")

    def test_gene_conversion_c2_tl_10_r10(self):
        self._run("100 10000 -t 5.0 -r 10 2501 -c 2 10")

    def test_gene_conversion_c2_tl_100(self):
        self._run("100 10000 -t 5.0 -r 10 2501 -c 2 100")

    def test_gene_conversion_c2_tl_100_r0(self):
        self._run("100 10000 -t 5.0 -r 0 2501 -c 2 100")

    def test_gene_conversion_c20_tl_1000_r0(self):
        self._run("100 10000 -t 5.0 -r 0 2501 -c 20 1000")


class MsDocExamples(MsTest):
    def test_msdoc_simple_ex(self):
        self._run("4 20000 -t 5.0")

    def test_msdoc_recomb_ex(self):
        self._run("15 1000 -t 10.04 -r 100.0 2501")

    def test_msdoc_structure_ex1(self):
        self._run("15 1000 -t 2.0 -I 3 10 4 1 5.0")

    def test_msdoc_structure_ex2(self):
        self._run("15 1000 -t 2.0 -I 3 10 4 1 5.0 -m 1 2 10.0 -m 2 1 9.0")

    def test_msdoc_structure_ex3(self):
        self._run("15 1000 -t 10.0 -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x")

    def test_msdoc_outgroup_sequence(self):
        self._run("11 1000 -t 2.0 -I 2 1 10 -ej 6.0 1 2")

    def test_msdoc_two_species(self):
        cmd = (
            "15 10000 -t 11.2 -I 2 3 12 -g 1 44.36 -n 2 "
            "0.125 -eg 0.03125 1 0.0 -en 0.0625 2 0.05 -ej 0.09375 2 1"
        )
        self._run(cmd)

    def test_msdoc_stepping_stone(self):
        cmd = (
            "15 10000 -t 3.0 -I 6 0 7 0 0 8 0 -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 "
            "-m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 "
            "-em 2.0 3 4 2.5 -em 2.0 4 3 2.5"
        )
        self._run(cmd)


class MsMiscExamples(MsTest):
    """
    Miscellaneous examples that have been good for finding bugs.
    """

    def test_simultaneous_ex1(self):
        self._run("10 10000 -t 2.0 -eN 0.3 0.5 -eG .3 7.0")

    def test_zero_growth_rate(self):
        self._run("10 10000 -t 2.0 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5")

    def test_konrad_1(self):
        cmd = (
            "4 1000 -t 2508 -I 2 2 2 0 -n 2 2.59 "
            "-ma x 0 1.502 x -ej 0.9485 1 2 -r 23.76 3000"
        )
        self._run(cmd)

    def test_konrad_2(self):
        cmd = (
            "3 10000 -t 0.423 -I 3 1 1 1 -es 0.0786 1 0.946635 "
            "-ej 0.0786 4 3 -ej 0.189256 1 2 -ej 0.483492 2 3"
        )
        self._run(cmd)

    def test_konrad_3(self):
        self._run("100 100 -t 2 -I 10 10 10 10 10 10 10 10 10 10 10 0.001 ")


class MsRandom(MsTest):
    """
    Some tests made by generating random parameters.
    """

    def _run(self, num_populations=1, num_replicates=1000, num_demographic_events=0):
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

        super()._run(cmd)

    def test_ms_random_1(self):
        self._run()

    def test_ms_random_2(self):
        self._run(num_replicates=10 ** 4, num_demographic_events=10)

    def test_ms_random_2_pops1(self):
        self._run(num_populations=3)


class MsHotTest(MsTest):
    def _run(self, cmd):
        self._run_variable_recombination_coalescent_stats(cmd)

    def test_mshotdoc_hotspot_ex(self):
        self._run("10 1000 -t 10.4 -r 10.0 25000 -v 2 100 200 10 7000 8000 20")

    def test_mshot_zero_recomb_interval(self):
        self._run("10 1000 -t 10.4 -r 10.0 25000 -v 1 5000 13000 0")

    def test_mshot_zero_recomb(self):
        self._run("10 1000 -t 10.4 -r 10.0 25000 -v 1 100 25000 0")

    def test_mshot_high_recomb_variants(self):
        hotspots = "4 1000 2000 0 7000 8000 20 12000 15000 10 20000 22000 0"
        cmd = f"10 1000 -t 10.4 -r 10.0 25000 -v {hotspots}"
        self._run(cmd)


class DiscoalTest(Test):
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
            _discoal_executable + args.split() + self.get_discoal_seeds()
        )

    def _run_mutation_discoal_stats(self, args):
        msp_str = self._discoal_str_to_ms(args)
        df_msp = self._run_msprime_mutation_stats(msp_str)
        df_d = self._run_sample_stats(
            _discoal_executable + args.split() + self.get_discoal_seeds()
        )
        self._plot_stats("mutation", df_d, df_msp, "discoal", "msp")

    def _discoal_str_to_simulation(self, args):
        # takes discoal command line as input
        # and returns msprime run treeseqs

        tokens = args.split(" ")
        # positional args
        sample_size = int(tokens[0])
        nreps = int(tokens[1])
        seq_length = int(tokens[2])
        # parse discoal command line for params
        # init ones we definitely need for comparison
        theta = rho = alpha = sweep_site = sweep_mod_time = None
        refsize = 1e6
        for i in range(3, len(tokens)):
            # pop size change case
            if tokens[i] == "-en":
                raise ValueError(
                    "sweeps with population size changes remain unimplemented"
                )
            # migration rate case
            if (tokens[i] == "-m") or (tokens[i] == "-p"):
                raise ValueError(
                    "sweeps with multiple populations remain unimplemented"
                )
            # split or admixture case
            if (tokens[i] == "-ea") or (tokens[i] == "-ed"):
                raise ValueError("sweeps with splits or admixture not supported")
            # sweep params
            if tokens[i] == "-x":
                sweep_site = float(tokens[i + 1])
            if (tokens[i] == "-ws") or (tokens[i] == "-wd") or (tokens[i] == "-wn"):
                sweep_mod_time = float(tokens[i + 1])
            if tokens[i] == "-a":
                alpha = float(tokens[i + 1])
            if tokens[i] == "-N":
                refsize = float(tokens[i + 1])
            # coalescent params
            if tokens[i] == "-t":
                theta = float(tokens[i + 1])
            if tokens[i] == "-r":
                rho = float(tokens[i + 1])
        mod_list = [("hudson")]
        if alpha is not None:
            # sweep model
            s = alpha / (2 * refsize)
            mod = msprime.SweepGenicSelection(
                position=np.floor(sweep_site * seq_length),
                start_frequency=1.0 / (2 * refsize),
                end_frequency=1.0 - (1.0 / (2 * refsize)),
                s=s * 2,  # discoal fitness model is 1, 1+s, 1+2s
                dt=1e-6,
            )
            mod_list.append((sweep_mod_time, mod))
            # if an event is defined from discoal line
            # best thing to do is rescale to Ne=0.25
            # so that time scale are consistent
            # see note at msprime/cli.py line 626
            # and following for alternate solution
            if sweep_mod_time > 0:
                refsize = 0.25
                mod.s = alpha / refsize
        # append final model
        mod_list.append((None, "hudson"))
        # scale theta and rho and create recomb_map
        recomb_map = msprime.RecombinationMap.uniform_map(
            seq_length, rho / 4 / refsize / (seq_length - 1)
        )
        mu = theta / 4 / refsize / seq_length
        replicates = msprime.simulate(
            sample_size,
            Ne=refsize,
            model=mod_list,
            recombination_map=recomb_map,
            mutation_rate=mu,
            num_replicates=nreps,
        )
        return replicates


class DiscoalCompatibility(DiscoalTest):
    """
    Basic tests to make sure that we have correctly set up the
    discoal interface.
    """

    def _run(self, cmd):
        self._run_mutation_discoal_stats(cmd)

    def test_discoal_simple_ex(self):
        self._run("15 1000 100 -t 5.0")

    def test_discoal_size_change1(self):
        self._run("10 10000 100 -t 10.0 -en 0.1 0 2.0")

    def test_discoal_size_change2(self):
        self._run("10 10000 100 -t 10.0 -en 0.1 0 0.1")

    def test_discoal_size_change3(self):
        self._run("10 10000 100 -t 10.0 -en 0.01 0 0.01")

    def test_discoal_size_change4(self):
        self._run("10 10000 100 -t 10.0 -en 0.01 0 0.5 -en 0.05 0 1.0")


# TODO we need to fix this test and to add a good number of examples.


class DiscoalSweeps(DiscoalTest):
    """
    Compare the result of sweeps in msprime and discoal.
    """

    def _run(self, args):
        df = pd.DataFrame()
        data = collections.defaultdict(list)
        replicates = self._discoal_str_to_simulation(args)
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
        self._plot_stats("mutation", df, df_df, "msp", "discoal")

    def test_sweep_ex0(self):
        cmd = "10 1000 10000 -t 10.0 -r 10.0"
        self._run(cmd)

    def test_sweep_no_rec_ex1(self):
        cmd = "10 1000 10000 -t 10.0 -r 0.0 -ws 0 -a 100 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_no_rec_ex2(self):
        cmd = "10 1000 10000 -t 10.0 -r 0.0 -ws 0 -a 200 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_rec_ex1(self):
        cmd = "10 1000 10000 -t 10.0 -r 10.0 -ws 0 -a 1000 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_rec_ex2(self):
        cmd = "10 1000 10000 -t 10.0 -r 20.0 -ws 0 -a 1000 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_rec_ex3(self):
        cmd = "10 1000 10000 -t 10.0 -r 100.0 -ws 0 -a 1000 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_rec_ex4(self):
        cmd = "10 1000 10000 -t 10.0 -r 400.0 -ws 0 -a 2000 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_rec_ex5(self):
        cmd = "10 1000 10000 -t 100.0 -r 100.0 -ws 0 -a 250 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_tau_ex1(self):
        cmd = "10 1000 10000 -t 10.0 -r 20.0 -ws 0.001 -a 250 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_tau_ex2(self):
        cmd = "10 1000 10000 -t 10.0 -r 20.0 -ws 0.01 -a 250 -x 0.5 -N 10000"
        self._run(cmd)

    def test_sweep_tau_ex3(self):
        cmd = "10 1000 10000 -t 10.0 -r 20.0 -ws 1.0 -a 250 -x 0.5 -N 10000"
        self._run(cmd)


def sample_recap_simplify(slim_ts, sample_size, Ne, r, mu):
    """
    takes a ts from slim and samples, recaps, simplifies
    """
    recap = msprime.sim_ancestry(
        initial_state=slim_ts,
        population_size=Ne,
        recombination_rate=r,
        start_time=slim_ts.metadata["SLiM"]["generation"],
    )
    rts = pyslim.SlimTreeSequence(recap)
    logging.debug(f"pyslim: slim generation:{slim_ts.metadata['SLiM']['generation']}")
    alive_inds = rts.individuals_alive_at(0)
    keep_indivs = np.random.choice(alive_inds, sample_size, replace=False)
    keep_nodes = []
    for i in keep_indivs:
        keep_nodes.extend(rts.individual(i).nodes)
    logging.debug(f"before simplify {rts.num_nodes} nodes")
    sts = rts.simplify(keep_nodes)
    logging.debug(f"after simplify {sts.num_nodes} nodes")
    logging.debug(f"after simplify {sts.num_trees} trees")
    return pyslim.SlimTreeSequence(msprime.mutate(sts, rate=mu))


class SweepVsSlim(Test):
    """
    Tests where we compare the msprime sweeps with SLiM simulations.
    """

    def run_sweep_slim_comparison(self, slim_args, **kwargs):
        df = pd.DataFrame()

        kwargs["model"] = "msp"
        logging.debug(f"Running: {kwargs}")
        seq_length = kwargs.get("seq_length")
        pop_size = kwargs.get("pop_size")
        s = kwargs.get("s")
        tau = kwargs.get("tau")
        sample_size = kwargs.get("sample_size")
        recombination_rate = kwargs.get("recombination_rate")
        num_replicates = kwargs.get("num_replicates")
        sweep = msprime.SweepGenicSelection(
            position=seq_length / 2,
            start_frequency=1.0 / (2 * pop_size),
            end_frequency=1.0 - (1.0 / (2 * pop_size)),
            s=s,
            dt=1e-6,
        )
        replicates = msprime.sim_ancestry(
            sample_size,
            population_size=pop_size,
            model=[None, (tau, sweep), (None, "hudson")],
            recombination_rate=recombination_rate,
            sequence_length=seq_length,
            num_replicates=num_replicates,
        )
        wins = range(0, int(seq_length + 1), int(seq_length / 20))
        mids = np.zeros(len(wins) - 1)
        for i in range(len(wins) - 1):
            mids[i] = (wins[i + 1] + wins[i]) / 2
        msp_win_pis = []
        slim_win_pis = []
        data = collections.defaultdict(list)
        for ts in replicates:
            t_mrca = np.zeros(ts.num_trees)
            for tree in ts.trees():
                t_mrca[tree.index] = tree.time(tree.root)
            data["tmrca_mean"].append(np.mean(t_mrca))
            data["num_trees"].append(ts.num_trees)
            mutated_ts = msprime.sim_mutations(ts, rate=1e-8)
            data["pi"].append(mutated_ts.diversity().reshape((1,))[0])
            data["model"].append("msp")
            msp_num_samples = ts.num_samples
            msp_win_pis.append(mutated_ts.diversity(windows=wins))
        slim_script = self.output_dir / "slim_script.txt"
        outfile = self.output_dir / "slim.trees"
        slim_args["OUTFILE"] = str(outfile)
        write_sweep_slim_script(slim_script, slim_args)

        cmd = _slim_executable + [slim_script]
        for _ in range(kwargs["num_replicates"]):
            subprocess.check_output(cmd)
            ts = pyslim.load(outfile)
            rts = sample_recap_simplify(
                ts, sample_size, pop_size, recombination_rate, 1e-8
            )
            assert rts.num_samples == msp_num_samples

            t_mrca = np.zeros(rts.num_trees)
            for tree in rts.trees():
                t_mrca[tree.index] = tree.time(tree.root)

            data["tmrca_mean"].append(np.mean(t_mrca))
            data["num_trees"].append(rts.num_trees)
            slim_win_pis.append(rts.diversity(windows=wins))
            data["pi"].append(rts.diversity().reshape((1,))[0])
            data["model"].append("slim")
        df = df.append(pd.DataFrame(data))

        df_slim = df[df.model == "slim"]
        df_msp = df[df.model == "msp"]
        for stat in ["tmrca_mean", "num_trees", "pi"]:
            v1 = df_slim[stat]
            v2 = df_msp[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            pyplot.xlabel("msp")
            pyplot.ylabel("SLiM")
            f = self.output_dir / f"{stat}.png"
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")
            plot_stat_hist(v1, v2, "slim", "msp")
            f = self.output_dir / f"{stat}.hist.png"
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")
        pyplot.plot(mids, np.array(msp_win_pis).mean(axis=0), label="msp")
        pyplot.plot(mids, np.array(slim_win_pis).mean(axis=0), label="slim")
        pyplot.title(f"tau: {tau}")
        pyplot.xlabel("location (bp)")
        pyplot.ylabel("pairwise diversity")
        pyplot.legend()
        f = self.output_dir / "pi_wins.png"
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")

    def _run(
        self,
        sample_size,
        seq_length,
        pop_size,
        recombination_rate,
        s,
        tau,
        num_replicates=None,
    ):
        """
        basic tests for sweeps vs slim
        """
        slim_args = {}

        if num_replicates is None:
            num_replicates = 20

        # These are *diploid* samples in msprime
        slim_args["sample_size"] = 2 * sample_size
        slim_args["r"] = recombination_rate
        slim_args["NUMLOCI"] = int(seq_length - 1)
        slim_args["POPSIZE"] = int(pop_size)
        slim_args["TAU"] = tau
        slim_args["s"] = s
        slim_args["SWEEPPOS"] = int(seq_length / 2)
        self.run_sweep_slim_comparison(
            slim_args,
            pop_size=pop_size,
            sample_size=sample_size,
            num_replicates=num_replicates,
            seq_length=seq_length,
            tau=tau,
            s=s,
            recombination_rate=recombination_rate,
        )

    def test_sweep_vs_slim_ex1(self):
        self._run(10, 1e6, 1e3, 1e-7, 0.25, 1, num_replicates=10)

    def test_sweep_vs_slim_ex2(self):
        self._run(10, 1e6, 1e3, 1e-7, 0.25, 200, num_replicates=10)

    def test_sweep_vs_slim_ex3(self):
        self._run(10, 1e6, 1e3, 1e-7, 0.25, 1000, num_replicates=10)

    def test_sweep_vs_slim_ex4(self):
        self._run(10, 1e6, 1e3, 1e-7, 0.25, 2000, num_replicates=10)

    def test_sweep_vs_slim_ex5(self):
        self._run(10, 1e6, 1e3, 1e-7, 0.25, 5000, num_replicates=10)


# FIXME disabling these for now because they are unreliable and result in
# errors (root: EXCEPTION:No columns to parse from file). This is probably
# harmless as it means there's no data simulated, but it stops the rest
# of the tests from running so we need to deal with it somehow.


class MsmsSweeps:
    """
    Compare msms with msprime/discoal for selective sweeps.

    NOTE:
    1. Msms allows user to specify selection starting time/frequency (-SI), or,
    alternatively, specify selection ending time/frequency (-SF); msprime is
    able to simulate selection similar to the '-SF' option in msms
    2. Msms allows user to specify different selection coefficients for AA and
    Aa, but in msprime/disocal only the selection coefficient for aA can be
    specified, and use h=0.5 to calculate that for AA.
    """

    def _msms_str_to_parameters(self, msms_cmd):
        """
        Parse msms cmdline arguments into a dictionary. This method is called
        by `_run_msp_sample_stats`

        msms cmdline pattern:
            nsam nrep -t theta -r rho num_sites -SF end_time end_frequency \
            -SAA sAA -SaA saA -Sp sel_pos -N refsize -seed rand_seed

        eg. "5 1 -t 200 -r 200 500000 -SF 0.002 0.9 -Sp 0.5"\
            " -SaA 5000 -SAA 10000 -N 10000 -seed 1"
        """
        # initialize local variables
        end_time_lst = []  # use list for multiple sweeps
        end_frequency_lst = []
        num_sweeps = 0
        sAA = saA = sel_pos = -1.0
        saA = -0.5
        refsize = 1
        rand_seed = (random.randint(1, 2 ** 16),)

        # parse arguments
        tokens = msms_cmd.split(" ")
        for ind in range(len(tokens)):
            if ind == 0:
                nsam = int(tokens[ind])
                nrep = int(tokens[ind + 1])
            elif tokens[ind] == "-t":
                theta = float(tokens[ind + 1])
            elif tokens[ind] == "-r":
                rho = float(tokens[ind + 1])
                num_sites = int(tokens[ind + 2])
            elif tokens[ind] == "-SF":
                num_sweeps += 1
                end_time_lst.append(float(tokens[ind + 1]))
                end_frequency_lst.append(float(tokens[ind + 2]))
            elif tokens[ind] == "-Sp":
                sel_pos = float(tokens[ind + 1])
            elif tokens[ind] == "-SAA":
                sAA = float(tokens[ind + 1])
            elif tokens[ind] == "-SaA":
                saA = float(tokens[ind + 1])
            elif tokens[ind] == "-N":
                refsize = int(tokens[ind + 1])
            elif tokens[ind] == "-seed":
                rand_seed = float(tokens[ind + 1])
            else:
                pass

        # check if h = 0.5
        if abs(saA * 2 - sAA) > 1e-5:
            logging.warning(
                "If 2 * saA is not equal to sAA, saA is set to sAA / 2,"
                "that is, h can only be 0.5 in msprime"
            )
            saA = sAA / 2.0

        return {
            "nsam": nsam,
            "nrep": nrep,
            "num_sweeps": num_sweeps,
            "end_time_lst": end_time_lst,
            "end_frequency_lst": end_frequency_lst,
            "refsize": refsize,
            "alpha": saA,
            "theta": theta,
            "rho": rho,
            "num_sites": num_sites,
            "sel_pos": sel_pos,
            "rand_seed": rand_seed,
        }

    def _update_msms_cmd_to_match_discoal(self, msms_cmd):
        """
        NOTE: discoal does not have options to specify allele frequencies
        and instead it calculates the frequency internally according to
        refsize. When msp or msms is compared with discoal, the
        "end_frequency_lst" from msms command arguments will be replaced by the
        following calculations.
        """
        msms_params = self._msms_str_to_parameters(msms_cmd)
        if msms_params["num_sweeps"] == 0:
            return msms_cmd

        logging.warning(
            "When compared with discoal, selected allele frequency options are"
            " recalculated following discoal's way"
        )

        # recalculate frequencies
        refsize = msms_params["refsize"]
        end_frequency_lst = msms_params["end_frequency_lst"]
        end_frequency_lst = [1 - 0.5 / refsize for _ in end_frequency_lst]

        # construct new msms cmd
        new_cmd = [
            str(msms_params["nsam"]),
            str(msms_params["nrep"]),
            "-t",
            str(msms_params["theta"]),
            "-r",
            str(msms_params["rho"]),
            str(msms_params["num_sites"]),
        ]
        for i in range(len(end_frequency_lst)):
            new_cmd += [
                "-SF",
                str(msms_params["end_time_lst"][i]),
                str(end_frequency_lst[i]),
                "-Sp",
                str(msms_params["sel_pos"]),
            ]
        new_cmd += [
            "-SaA",
            str(msms_params["alpha"]),
            "-SAA",
            str(msms_params["alpha"] * 2),
            "-N",
            str(msms_params["refsize"]),
        ]
        new_msms_cmd = " ".join(new_cmd)
        return new_msms_cmd

    def _msms_params_to_run_msp(self, params):
        """
        Run simulation for a single sample and return a tree sequence. This
        method is called by `_run_msp_sample_stats` in a loop to generate nrep
        samples.
        """
        if params["num_sweeps"] > 0:
            model = [None]
            for i in range(params["num_sweeps"]):
                temp_model = msprime.SweepGenicSelection(
                    position=params["sel_pos"],
                    end_frequency=params["end_frequency_lst"][i],
                    start_frequency=0.5 / params["refsize"],
                    alpha=params["alpha"],
                    dt=1.0 / (40 * params["refsize"]),
                )
                t_start = params["end_time_lst"][i]
                model.append((t_start, temp_model))
            model.append((None, None))

        else:
            model = "hudson"

        scale_factor = params["num_sites"]
        recombination_rate = params["rho"] / (scale_factor - 1)
        mutation_rate = params["theta"] / scale_factor

        ts = msprime.simulate(
            sample_size=params["nsam"],
            Ne=0.25,
            length=params["num_sites"],
            mutation_rate=mutation_rate,
            recombination_rate=recombination_rate,
            model=model,
        )

        return ts

    def _run_msp_sample_stats(self, msms_cmd):
        """
        Call methods to parse cmdline options and run simulation,
        and then output in ms format, pipe thru sample_stats and finally return
        stats dataframe.
        """
        temp_file = tempfile.gettempdir() + "/tmp_msp_out"
        output = open(temp_file, "w")

        # run simulation and print ms format data into a file
        msms_params = self._msms_str_to_parameters(msms_cmd)
        num_replicates = msms_params["nrep"]
        print("ms " + msms_cmd, file=output)  # needed by sample_stat tools
        self._ms_random_seeds = msms_params["rand_seed"] = self.get_ms_seeds()

        for _ in range(num_replicates):
            tree_sequence = self._msms_params_to_run_msp(msms_params)
            print(file=output)
            print("//", file=output)
            if msms_params["theta"] > 0:
                s = tree_sequence.get_num_mutations()
                print("segsites:", s, file=output)

                if s != 0:
                    print("positions: ", end="", file=output)
                    positions = [
                        mutation.position / msms_params["num_sites"]
                        for mutation in tree_sequence.mutations()
                    ]
                    positions.sort()
                    for position in positions:
                        print("{0:.{1}f}".format(position, 8), end=" ", file=output)
                    print(file=output)
                    for h in tree_sequence.haplotypes():
                        print(h, file=output)

                else:
                    print(file=output)
        output.close()

        # pipe ms format output to sample_stats
        p1 = subprocess.Popen(["cat", temp_file], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["./data/sample_stats"], stdin=p1.stdout, stdout=subprocess.PIPE
        )
        p1.stdout.close()
        output = p2.communicate()[0]
        p1.wait()

        # read into pandas frame and return it
        with tempfile.TemporaryFile() as f:
            f.write(output)
            f.seek(0)
            df = pd.read_csv(f, sep="\t")
        return df

    def _run_msms_sample_stats(self, cmd):
        return self._run_sample_stats(_msms_executable + cmd.split(" "))

    def _convert_to_discoal_cmd(self, msms_cmd):
        """
        called by _run_discoal_sample_stats to convert msms cmdline args
        to discoal cmdline args

        NOTE: if -N option is not specified, discoal internally use N=1,000,000
        """
        params = self._msms_str_to_parameters(msms_cmd)
        return "%d %d %d -t %f -r %f -ws %f -a %f -x %f -N %f" % (
            params["nsam"],
            params["nrep"],
            params["num_sites"],
            params["theta"],
            params["rho"],
            params["end_time_lst"][0],
            params["alpha"],
            params["sel_pos"],
            params["refsize"],
        )

    def _run_discoal_sample_stats(self, msms_cmd):
        discoal_cmd = self._convert_to_discoal_cmd(msms_cmd)
        return self._run_sample_stats(_discoal_executable + discoal_cmd.split(" "))

    def _cmp_msms_vs_msp(self, cmd):
        df_msp = self._run_msp_sample_stats(cmd)
        df_msms = self._run_msms_sample_stats(cmd)
        self._plot_stats("msp_msms", df_msp, df_msms, "msp", "msms")

    def _cmp_discoal_vs_msp_via_msms_cmd(self, cmd):
        cmd = self._update_msms_cmd_to_match_discoal(cmd)
        df_discoal = self._run_discoal_sample_stats(cmd)
        df_msp = self._run_msp_sample_stats(cmd)
        self._plot_stats("msp_discoal", df_msp, df_discoal, "msp", "discoal")

    def _cmp_msms_vs_discoal(self, cmd):
        cmd = self._update_msms_cmd_to_match_discoal(cmd)
        df_discoal = self._run_discoal_sample_stats(cmd)
        df_msms = self._run_msms_sample_stats(cmd)
        self._plot_stats("discoal_msms", df_discoal, df_msms, "discoal", "msms")

    def test_neutral_msms_vs_msp(self):
        self._cmp_msms_vs_msp("100 300 -t 200 -r 200 500000 -N 10000")

    def test_selective_discoal_vs_msp(self):
        self._cmp_discoal_vs_msp_via_msms_cmd(
            "100 300 -t 20 -r 20 50000"
            " -SF 0 0.99995 -Sp 0.5 -SaA 5000 -SAA 10000 -N 10000"
        )

    def test_selective_msms_vs_msp(self):
        self._cmp_msms_vs_msp(
            "100 300 -t 200 -r 200 500000"
            " -SF 0 0.9 -Sp 0.5 -SaA 5000 -SAA 10000 -N 10000"
        )

    def test_selective_msms_vs_msp_small_s(self):
        self._cmp_msms_vs_msp(
            "100 300 -t 200 -r 200 500000 -SF 0 0.9 -Sp 0.5 -SaA 1 -SAA 2 -N 10000"
        )

    def test_selective_msms_vs_msp_multiple_sweeps(self):
        self._cmp_msms_vs_msp(
            "100 300 -t 200 -r 200 500000"
            " -SF 0 0.9 -Sp 0.5"
            " -SF 0.1 0.9 -Sp 0.5 -SaA 5000 -SAA 10000 -N 10000"
        )

    def _test_selective_msp_50Mb(self):
        """
        Test runtime of msprime for long chromosomes
        """
        self._cmp_msp_sample_stats(
            "1000 1 -t 20000 -r 20000 50000000"
            " -SF 0 0.9 -Sp 0.5 -SaA 5000 -SAA 10000 -N 10000"
        )

    def test_selective_msms_vs_discoal(self):
        self._cmp_msms_vs_discoal(
            # "100 300 -t 20 -r 20 50000"
            "100 300 -t 20 -r 20 5000"
            " -SF 0 0.9 -Sp 0.5 -SaA 5000 -SAA 10000 -N 10000"
        )

    def test_selective_msms_vs_msp_use_discoal_paper_param(self):
        self._cmp_msms_vs_msp(
            "100 300 -t 100 -r 100 250000"
            " -SF 0 0.99995 -Sp 0.5 -SaA 2000 -SAA 4000 -N 10000"
        )

    def test_selective_msms_vs_discoal_use_discoal_paper_param(self):
        """
        NOTE: tests calling discoal will take a much longer time to finish
        especially when large num_sites are used. Use the lines commented out
        instead if we want to reproduce the results posted in issue # 1173
        """
        self._cmp_msms_vs_discoal(
            # "100 300 -t 100 -r 100 250000"
            "100 300 -t 100 -r 100 2500"
            " -SF 0 0.99995 -Sp 0.5 -SaA 2000 -SAA 4000 -N 10000"
        )

    def test_selective_msms_vs_discoal_random_param(self):
        self._cmp_msms_vs_discoal(
            # "100 300 -t 40 -r 40 50000"
            "100 300 -t 40 -r 40 5000"
            " -SF 0 0.99995 -Sp 0.5 -SaA 1000 -SAA 2000 -N 10000"
        )

    def test_selective_discoal_vs_msp_use_discoal_paper_param(self):
        self._cmp_discoal_vs_msp_via_msms_cmd(
            # "100 300 -t 100 -r 100 250000"
            "100 300 -t 100 -r 100 2500"
            " -SF 0 0.99995 -Sp 0.5 -SaA 2000 -SAA 4000 -N 10000"
        )

    def test_selective_discoal_vs_msp_random_param(self):
        self._cmp_discoal_vs_msp_via_msms_cmd(
            # "100 300 -t 40 -r 40 50000"
            "100 300 -t 40 -r 40 5000"
            " -SF 0 0.99995 -Sp 0.5 -SaA 1000 -SAA 2000 -N 10000"
        )


class SweepAnalytical(Test):
    """
    Analytical comparisons wrt to sweeps
    """

    def hermissonPennings_exp_sojourn(self, alpha):
        """
        analytic expectation of sojourn time
        equation A.17 from Hermisson and Pennings
        """
        inner = np.log(alpha) + np.euler_gamma - (1.0 / alpha)
        return 4.0 / alpha * inner

    def charlesworth_exp_sojourn(self, alpha, s):
        """
        same as above but scaled in number of gens
        """
        inner = np.log(alpha) + np.euler_gamma - (1.0 / alpha)
        return 4.0 / s * inner

    def test_sojourn_time(self):
        alphas = np.arange(5e-3, 5e-2, 5e-3)
        refsize = 1e4
        nreps = 500
        seqlen = 1e4
        mu = 2.5e-8
        rho = 0
        p0 = 1.0 / (2 * refsize)
        p1 = 1 - p0
        dt = 1.0 / (400 * refsize)
        pos = np.floor(seqlen / 2)
        df = pd.DataFrame()
        data = collections.defaultdict(list)
        for a in alphas:
            mod = msprime.SweepGenicSelection(
                start_frequency=p0, end_frequency=p1, s=a, dt=dt, position=pos
            )
            s = a / 2 / refsize
            replicates = msprime.simulate(
                10,
                Ne=refsize,
                model=[mod],
                length=seqlen,
                num_labels=2,
                recombination_rate=rho,
                mutation_rate=mu,
                num_replicates=nreps,
            )

            reptimes = np.zeros(nreps)
            i = 0
            for x in replicates:
                tree_times = np.zeros(x.num_trees)
                j = 0
                for tree in x.trees():
                    tree_times[j] = np.max([tree.time(root) for root in tree.roots])
                    j += 1
                reptimes[i] = np.max(tree_times)
                i += 1
            data["alpha_means"].append(np.mean(reptimes))
            data["exp_means"].append(self.charlesworth_exp_sojourn(a, s))
        df = pd.DataFrame.from_dict(data)
        df = df.fillna(0)
        sm.qqplot_2samples(df["exp_means"], df["alpha_means"], line="45")
        pyplot.xlabel("expected sojourn time")
        pyplot.ylabel("simulated sojourn time")
        f = self.output_dir / "sojourn.png"
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")

    def test_sojourn_time2(self):
        alpha = 1000
        refsizes = [0.25, 0.5, 1.0]
        selrefsize = 1000
        nreps = 500
        seqlen = 1e4
        mu = 2.5e-8
        rho = 0
        p0 = 1.0 / (2 * selrefsize)
        p1 = 1 - p0
        dt = 1.0 / (400 * selrefsize)
        pos = np.floor(seqlen / 2)
        df = pd.DataFrame()
        data = collections.defaultdict(list)
        for n in refsizes:
            s = alpha / (2 * n)
            mod = msprime.SweepGenicSelection(
                start_frequency=p0, end_frequency=p1, s=s, dt=dt, position=pos
            )
            replicates = msprime.simulate(
                10,
                Ne=n,
                model=[mod],
                length=seqlen,
                num_labels=2,
                recombination_rate=rho,
                mutation_rate=mu,
                num_replicates=nreps,
            )

            reptimes = np.zeros(nreps)
            i = 0
            for x in replicates:
                tree_times = np.zeros(x.num_trees)
                j = 0
                for tree in x.trees():
                    tree_times[j] = np.max([tree.time(root) for root in tree.roots])
                    j += 1
                reptimes[i] = np.max(tree_times)
                i += 1
            data["alpha_means"].append(np.mean(reptimes))
            data["exp_means"].append(self.hermissonPennings_exp_sojourn(alpha) * 2 * n)
        df = pd.DataFrame.from_dict(data)
        df = df.fillna(0)
        sm.qqplot_2samples(df["exp_means"], df["alpha_means"], line="45")
        pyplot.xlabel("expected sojourn time")
        pyplot.ylabel("simulated sojourn time")
        f = self.output_dir / "sojourn.png"
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")


# FIXME disabling these for now because the pedigree file that
# they depend on doesn't exist. (Tests won't be picked up unless
# they subclass Test.)


class DtwfPedigreeVsCoalescent:
    def run_dtwf_pedigree_comparison(self, **kwargs):
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
            replicates = msprime.simulate(**kwargs)
            for ts in replicates:
                t_mrca = np.zeros(ts.num_trees)
                for tree in ts.trees():
                    t_mrca[tree.index] = tree.time(tree.root)
                data["tmrca_mean"].append(np.mean(t_mrca))
                data["num_trees"].append(ts.num_trees)
                data["model"].append(model)
            df = df.append(pd.DataFrame(data))

        df_wf_ped = df[df.model == "wf_ped"]
        df_dtwf = df[df.model == "dtwf"]
        for stat in ["tmrca_mean", "num_trees"]:
            v1 = df_wf_ped[stat]
            v2 = df_dtwf[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = self.output_dir / f"{stat}.png"
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def test_dtwf_vs_pedigree_single_locus(self):
        pedigree_file = "tests/data/pedigrees/wf_100Ne_10000gens.txt"
        pedigree = msprime.Pedigree.read_txt(pedigree_file, time_col=3)

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

    def test_dtwf_vs_pedigree_short_region(self):
        pedigree_file = "tests/data/pedigrees/wf_100Ne_10000gens.txt"
        pedigree = msprime.Pedigree.read_txt(pedigree_file, time_col=3)

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

    def test_dtwf_vs_pedigree_long_region(self):
        pedigree_file = "tests/data/pedigrees/wf_100Ne_10000gens.txt"
        pedigree = msprime.Pedigree.read_txt(pedigree_file, time_col=3)

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


class DtwfVsCoalescent(Test):
    """
    Tests where we compare the DTWF with coalescent simulations.
    """

    def run_dtwf_coalescent_stats(self, **kwargs):
        df = pd.DataFrame()

        for model in ["hudson", "dtwf"]:
            kwargs["model"] = model

            logging.debug(f"Running: {kwargs}")
            data = collections.defaultdict(list)
            replicates = msprime.sim_ancestry(**kwargs)
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

    def plot_dtwf_coalescent_stats(self, df):
        df_hudson = df[df.model == "hudson"]
        df_dtwf = df[df.model == "dtwf"]
        for stat in ["tmrca_mean", "num_trees"]:
            plot_qq(df_hudson[stat], df_dtwf[stat])
            f = self.output_dir / f"{stat}.png"
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

        hudson_breakpoints = all_breakpoints_in_replicates(df_hudson["intervals"])
        dtwf_breakpoints = all_breakpoints_in_replicates(df_dtwf["intervals"])
        if len(hudson_breakpoints) > 0 or len(dtwf_breakpoints) > 0:
            plot_breakpoints_hist(
                hudson_breakpoints, dtwf_breakpoints, "hudson", "dtwf"
            )
            pyplot.savefig(self.output_dir / "breakpoints.png", dpi=72)
            pyplot.close("all")

    def plot_tree_intervals(self, df):
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
        pyplot.savefig(self.output_dir / "intervals.png", dpi=72)
        pyplot.close("all")

    def _run(self, **kwargs):
        df = self.run_dtwf_coalescent_stats(**kwargs)
        self.plot_dtwf_coalescent_stats(df)
        self.plot_tree_intervals(df)


class DtwfVsCoalescentSimple(DtwfVsCoalescent):
    """
    Straightforward tests where we pass through simulate args directly.
    """

    def test_dtwf_vs_coalescent_single_locus(self):
        self._run(samples=10, population_size=1000, num_replicates=300)

    def test_dtwf_vs_coalescent_recomb_discrete_hotspots(self):
        """
        Checks the DTWF against the standard coalescent with a
        discrete recombination map with variable rates.
        """
        recombination_map = msprime.RateMap(
            position=[0, 100, 500, 900, 1200, 1500, 2000],
            rate=[0.00001, 0, 0.0002, 0.00005, 0, 0.001],
        )
        self._run(
            samples=10,
            population_size=1000,
            recombination_rate=recombination_map,
            num_replicates=300,
            discrete_genome=True,
        )

    def test_dtwf_vs_coalescent_recomb_continuous_hotspots(self):
        """
        Checks the DTWF against the standard coalescent with a
        continuous recombination map with variable rates.
        """
        recombination_map = msprime.RateMap(
            position=[0, 0.1, 0.5, 0.9, 1.2, 1.5, 2.0],
            rate=[0.00001, 0, 0.0002, 0.00005, 0, 0.001],
        )
        self._run(
            samples=10,
            population_size=1000,
            recombination_rate=recombination_map,
            num_replicates=300,
            discrete_genome=False,
        )

    def test_dtwf_vs_coalescent_single_forced_recombination(self):
        recombination_map = msprime.RateMap(position=[0, 100, 101, 201], rate=[0, 1, 0])
        self._run(
            samples=10,
            population_size=10,
            num_replicates=1,
            discrete_genome=True,
            recombination_rate=recombination_map,
        )

    def test_dtwf_vs_coalescent_low_recombination(self):
        self._run(
            samples=10,
            population_size=1000,
            num_replicates=400,
            recombination_rate=0.01,
            sequence_length=5,
        )

    def test_dtwf_vs_coalescent_2_pops_massmigration(self):
        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.events.append(
            msprime.MassMigration(time=300, source=1, destination=0, proportion=1.0)
        )
        self._run(
            samples={0: 10, 1: 10},
            demography=demography,
            sequence_length=10 ** 6,
            num_replicates=300,
            recombination_rate=1e-8,
        )

    def test_dtwf_vs_coalescent_1_pop_growth(self):
        self._run(
            samples=10,
            demography=msprime.Demography.isolated_model([1000], growth_rate=[0.01]),
            recombination_rate=1e-8,
            sequence_length=5e7,
            num_replicates=300,
            discrete_genome=True,
        )

    def test_dtwf_vs_coalescent_1_pop_shrink(self):
        initial_size = 1000
        demography = msprime.Demography.isolated_model(
            [initial_size], growth_rate=[-0.01]
        )
        demography.events.append(
            msprime.PopulationParametersChange(
                time=200, initial_size=initial_size, growth_rate=0.01, population=0
            )
        )
        self._run(
            samples=10,
            demography=demography,
            recombination_rate=1e-8,
            sequence_length=5e7,
            num_replicates=300,
            discrete_genome=True,
        )

    def test_dtwf_vs_coalescent_multiple_bottleneck(self):
        demography = msprime.Demography.isolated_model([1000, 1000])
        demography.events = [
            msprime.PopulationParametersChange(
                time=100, initial_size=100, growth_rate=-0.01, population=0
            ),
            msprime.PopulationParametersChange(
                time=200, initial_size=100, growth_rate=-0.01, population=1
            ),
            msprime.PopulationParametersChange(
                time=300, initial_size=1000, growth_rate=0.01, population=0
            ),
            msprime.PopulationParametersChange(
                time=400, initial_size=1000, growth_rate=0.01, population=1
            ),
            msprime.PopulationParametersChange(
                time=500, initial_size=100, growth_rate=0, population=0
            ),
            msprime.PopulationParametersChange(
                time=600, initial_size=100, growth_rate=0, population=1
            ),
            msprime.MigrationRateChange(time=700, rate=0.1, matrix_index=(0, 1)),
        ]
        self._run(
            samples={0: 5, 1: 5},
            demography=demography,
            num_replicates=400,
            recombination_rate=1e-8,
            sequence_length=5e7,
        )


class DtwfVsCoalescentHighLevel(DtwfVsCoalescent):
    """
    Tests for the DTWF and coalescent when we use a slightly more
    high-level intervace.
    """

    def _run(
        self,
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

        demography = msprime.Demography.isolated_model(
            initial_sizes, growth_rate=growth_rates
        )

        for i in range(num_pops):
            if initial_sizes[i] > 100:
                # Growth rate set to zero at pop size 100
                t_100 = (np.log(initial_sizes[i]) - np.log(100)) / growth_rates[i]
                de = msprime.PopulationParametersChange(
                    t_100, growth_rate=0, population=i
                )
                demography.events.append(de)

            else:
                # Enforce zero growth rate for small populations
                logging.warning(
                    f"Warning - setting growth rate to zero for small \
                    population of size {initial_sizes[i]}",
                )
                demography.populations[i].growth_rate = 0

        if migration_matrix is None:
            default_mig_rate = 0.05
            migration_matrix = []
            for i in range(num_pops):
                row = [default_mig_rate] * num_pops
                row[i] = 0
                migration_matrix.append(row)

        demography.migration_matrix = migration_matrix

        super()._run(
            samples={j: sample_size for j, sample_size in enumerate(sample_sizes)},
            demography=demography,
            num_replicates=num_replicates,
            sequence_length=num_loci,
            recombination_rate=recombination_rate,
            discrete_genome=True,
        )

    def test_dtwf_vs_coalescent_long_region(self):
        self._run([1000], [10], int(1e8), 1e-8)

    def test_dtwf_vs_coalescent_short_region(self):
        self._run([1000], [10], int(1e6), 1e-8)

    def test_dtwf_vs_coalescent_2_pops(self):
        self._run(
            [500, 500],
            [5, 5],
            int(1e6),
            1e-8,
            num_replicates=500,
        )

    def test_dtwf_vs_coalescent_3_pops(self):
        self._run(
            [500, 500, 500],
            [5, 2, 0],
            int(1e7),
            1e-8,
        )

    def test_dtwf_vs_coalescent_4_pops(self):
        self._run(
            [1000, 1000, 1000, 1000],
            [0, 20, 0, 0],
            int(1e6),
            1e-8,
            num_replicates=500,
        )

    def test_dtwf_vs_coalescent_3_pops_asymm_mig(self):
        migration_matrix = [[0, 0.2, 0.1], [0.1, 0, 0.2], [0.2, 0.1, 0]]
        self._run(
            [500, 500, 500],
            [20, 0, 0],
            int(1e6),
            1e-8,
            migration_matrix=migration_matrix,
            num_replicates=500,
        )

    def test_dtwf_vs_coalescent_2_pops_high_asymm_mig(self):

        migration_matrix = [[0, 0.5], [0.7, 0]]
        self._run(
            [1000, 1000],
            [10, 10],
            int(1e6),
            1e-8,
            migration_matrix=migration_matrix,
            num_replicates=200,
            growth_rates=[0.005, 0.005],
        )


class DtwfVsSlim(Test):
    """
    Tests where we compare the DTWF with SLiM simulations.
    """

    def run_dtwf_slim_comparison(self, slim_args, **kwargs):
        df = pd.DataFrame()

        kwargs["model"] = "dtwf"
        logging.debug(f"Running: {kwargs}")
        replicates = msprime.sim_ancestry(**kwargs)
        data = collections.defaultdict(list)
        for ts in replicates:
            t_mrca = np.zeros(ts.num_trees)
            for tree in ts.trees():
                t_mrca[tree.index] = tree.time(tree.root)
            data["tmrca_mean"].append(np.mean(t_mrca))
            data["num_trees"].append(ts.num_trees)
            data["model"].append("dtwf")
            msp_num_samples = ts.num_samples

        slim_script = self.output_dir / "slim_script.txt"
        outfile = self.output_dir / "slim.trees"
        slim_args["OUTFILE"] = str(outfile)
        write_slim_script(slim_script, slim_args)

        cmd = _slim_executable + [slim_script]
        for _ in range(kwargs["num_replicates"]):
            subprocess.check_output(cmd)
            ts = tskit.load(outfile)
            ts = subsample_simplify_slim_treesequence(ts, slim_args["sample_sizes"])
            assert ts.num_samples == msp_num_samples

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
            pyplot.xlabel("DTWF")
            pyplot.ylabel("SLiM")
            f = self.output_dir / f"{stat}.png"
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def check_slim_version(self):
        # This may not be robust but it's a start
        min_version = 3.1
        raw_str = subprocess.check_output(_slim_executable + ["-version"])
        version_list = str.split(str(raw_str))
        for i in range(len(version_list)):
            if version_list[i].lower() == "version":
                version_str = version_list[i + 1]
                break
        version = float(version_str.strip(" ,")[0:3])
        assert version >= min_version, "Require SLiM >= 3.1!"

    def _run(
        self,
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

        sample_sizes = np.array(sample_sizes)
        num_pops = len(sample_sizes)
        slim_args = {}

        if num_replicates is None:
            num_replicates = 200

        # These are *diploid* samples in msprime
        slim_args["sample_sizes"] = 2 * sample_sizes
        demography = msprime.Demography.isolated_model(initial_sizes)

        slim_args["POP_STRS"] = ""
        for i in range(num_pops):
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
        demography.migration_matrix = migration_matrix

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

        slim_args["RHO"] = recombination_rate
        slim_args["NUM_LOCI"] = int(num_loci)

        self.run_dtwf_slim_comparison(
            slim_args,
            samples={j: sample_size for j, sample_size in enumerate(sample_sizes)},
            demography=demography,
            num_replicates=num_replicates,
            sequence_length=num_loci,
            recombination_rate=recombination_rate,
            discrete_genome=True,
        )

    def test_dtwf_vs_slim_single_locus(self):
        self._run([100], [10], 1, 0)

    def test_dtwf_vs_slim_single_locus_2_pops(self):
        self._run([20, 20], [5, 5], 1, 0)

    def test_dtwf_vs_slim_short_region(self):
        self._run([100], [10], 1e7, 1e-8, num_replicates=200)

    def test_dtwf_vs_slim_long_region(self):
        self._run([50], [10], 1e8, 1e-8, num_replicates=200)


class DtwfVsCoalescentRandom(DtwfVsCoalescent):
    """
    Runs randomly generated test parameters.
    """

    def _run(self, num_populations=1, num_replicates=200, num_demographic_events=0):

        # Make this deterministic
        np.random.seed(42)
        random.seed(42)

        N = num_populations
        num_loci = np.random.randint(1e5, 1e7)
        num_samples = np.random.randint(2, 10, size=num_populations)
        demography = msprime.Demography.isolated_model([1000 / N] * num_populations)

        migration_matrix = []
        for i in range(N):
            migration_matrix.append(
                [random.uniform(0.05, 0.25) * (j != i) for j in range(N)]
            )
        demography.migration_matrix = migration_matrix

        # Add demographic events and some migration rate changes
        t_max = 1000
        times = sorted(np.random.randint(300, t_max, size=num_demographic_events))
        for t in times:
            initial_size = np.random.randint(500, 1000)
            # Setting growth_rate to 0 because it's too tricky to get
            # growth_rates in the DTWF which don't result in N going to 0.
            growth_rate = 0
            pop_id = np.random.randint(N)
            demography.events.append(
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
                demography.events.append(
                    msprime.MigrationRateChange(time=t, rate=rate, matrix_index=index)
                )

        # Collect all pops together to control coalescence times for DTWF
        for i in range(1, N):
            demography.events.append(
                msprime.MassMigration(
                    time=t_max, source=i, destination=0, proportion=1.0
                )
            )

        demography.events.append(
            msprime.PopulationParametersChange(
                time=t_max, initial_size=100, growth_rate=0, population_id=0
            )
        )
        super()._run(
            samples={j: sample_size for j, sample_size in enumerate(num_samples)},
            demography=demography,
            num_replicates=num_replicates,
            sequence_length=num_loci,
            recombination_rate=1e-8,
            discrete_genome=True,
        )

    def test_dtwf_vs_coalescent_random_1(self):
        self._run(num_populations=2, num_replicates=200, num_demographic_events=3)

    def test_dtwf_vs_coalescent_random_2(self):
        self._run(num_populations=3, num_replicates=200, num_demographic_events=3)

    def test_dtwf_vs_coalescent_random_3(self):
        self._run(num_populations=2, num_replicates=200, num_demographic_events=6)

    def test_dtwf_vs_coalescent_random_4(self):
        self._run(num_populations=1, num_replicates=200, num_demographic_events=8)


class RecombinationBreakpointTest(Test):
    """
    Verifies that the number of recombination breakpoints is proportional to
    the total branch length across all trees.
    """

    def verify_breakpoint_distribution(
        self, name, sample_size, Ne, r, L, ploidy, model, growth_rate=0
    ):
        ts = msprime.sim_ancestry(
            samples=sample_size,
            demography=msprime.Demography.isolated_model(
                [Ne], growth_rate=[growth_rate]
            ),
            ploidy=ploidy,
            sequence_length=L,
            recombination_rate=r,
            model=model,
        )
        area = [tree.total_branch_length * tree.span for tree in ts.trees()]
        scipy.stats.probplot(area, dist=scipy.stats.expon(Ne * r), plot=pyplot)
        path = self.output_dir / f"{name}_growth={growth_rate}_ploidy={ploidy}.png"
        logging.debug(f"Writing {path}")
        pyplot.savefig(path)
        pyplot.close("all")

    def test_xi_beta_breakpoints(self):
        Ne = 10 ** 4
        for alpha in [1.1, 1.3, 1.6, 1.9]:
            for p in [1, 2]:
                self.verify_breakpoint_distribution(
                    f"n=100_alpha={alpha}",
                    sample_size=100,
                    Ne=Ne,
                    r=1e-7,
                    L=10 ** 6,
                    ploidy=p,
                    model=msprime.BetaCoalescent(alpha=alpha),
                )
                # Add a growth rate with a higher recombination rate so
                # we still get decent numbers of trees
                self.verify_breakpoint_distribution(
                    f"growth_n=100_alpha={alpha}",
                    sample_size=100,
                    Ne=Ne,
                    r=1e-7,
                    L=10 ** 6,
                    ploidy=p,
                    model=msprime.BetaCoalescent(alpha=alpha),
                    growth_rate=0.05,
                )

    def test_xi_dirac_breakpoints(self):
        Ne = 10 ** 2
        for psi in [0.1, 0.3, 0.6, 0.9]:
            for c in [1, 10]:
                for p in [1, 2]:
                    self.verify_breakpoint_distribution(
                        f"n=100_psi={psi}_c={c}",
                        sample_size=100,
                        Ne=Ne,
                        r=1e-8,
                        L=10 ** 6,
                        ploidy=p,
                        model=msprime.DiracCoalescent(psi=psi, c=c),
                    )
                    # Add a growth rate with a higher recombination rate so
                    # we still get decent numbers of trees
                    self.verify_breakpoint_distribution(
                        f"growth_n=100_psi={psi}_c={c}",
                        sample_size=100,
                        Ne=Ne,
                        r=1e-7,
                        L=10 ** 6,
                        ploidy=p,
                        model=msprime.DiracCoalescent(psi=psi, c=c),
                        growth_rate=0.05,
                    )

    def test_hudson_breakpoints(self):
        for p in [1, 2]:
            self.verify_breakpoint_distribution(
                "single_pop_n_50",
                sample_size=50,
                Ne=10 ** 4,
                r=1e-8,
                L=10 ** 6,
                ploidy=p,
                model="hudson",
            )
            self.verify_breakpoint_distribution(
                "single_pop_n_100",
                sample_size=100,
                Ne=10 ** 4,
                r=1e-8,
                L=10 ** 6,
                ploidy=p,
                model="hudson",
            )
            self.verify_breakpoint_distribution(
                "single_pop_n_100_growth",
                sample_size=100,
                Ne=10 ** 4,
                r=1e-7,
                L=10 ** 6,
                ploidy=p,
                model="hudson",
                growth_rate=0.05,
            )


class RecombinationMutationTest(Test):
    """
    Verifies that the number of recombinations equals the number of mutations
    since both should be proportional to the total branch lenght of the
    trees.
    """

    def verify_recombination(
        self, name, sample_size, Ne, r, m, L, ploidy, model, growth_rate=0
    ):
        num_replicates = 500
        empirical_theta = []
        empirical_rho = []
        for _ in range(num_replicates):
            sim = msprime.ancestry._parse_sim_ancestry(
                samples=[msprime.SampleSet(sample_size, ploidy=1)],
                recombination_rate=r,
                sequence_length=L,
                ploidy=ploidy,
                demography=msprime.Demography.isolated_model(
                    [Ne], growth_rate=[growth_rate]
                ),
                model=model,
            )
            ts = next(sim.run_replicates(1))
            empirical_rho.append(sim.num_breakpoints)
            ts = msprime.sim_mutations(ts, rate=m)
            empirical_theta.append(ts.get_num_sites())
        empirical_rho.sort()
        empirical_theta.sort()
        empirical_rho = np.array(empirical_rho)
        empirical_theta = np.array(empirical_theta)
        plot_qq(empirical_theta, empirical_rho)
        path = (
            self.output_dir
            / f"{name}_growth={growth_rate}_ploidy={ploidy}_rec_check.png"
        )
        logging.debug(f"Writing {path}")
        pyplot.savefig(path)
        pyplot.close("all")

    def test_xi_beta_recombinations(self):
        Ne = 10000
        for alpha in [1.1, 1.3, 1.5, 1.9]:
            for p in [1, 2]:
                self.verify_recombination(
                    f"n=100_alpha={alpha}",
                    sample_size=100,
                    Ne=Ne,
                    r=1e-8,
                    m=1e-8,
                    L=10 ** 6,
                    ploidy=p,
                    model=msprime.BetaCoalescent(alpha=alpha),
                )

    def test_xi_dirac_recombinations(self):
        Ne = 100
        for psi in [0.1, 0.5, 0.9]:
            for c in [1, 10]:
                for p in [1, 2]:
                    self.verify_recombination(
                        f"n=100_psi={psi}_c={c}",
                        sample_size=100,
                        Ne=Ne,
                        r=1e-8,
                        m=1e-8,
                        L=10 ** 6,
                        ploidy=p,
                        model=msprime.DiracCoalescent(psi=psi, c=c),
                    )

    def test_hudson_recombinations(self):
        for p in [1, 2]:
            self.verify_recombination(
                "n=100_hudson",
                sample_size=100,
                Ne=10000,
                r=1e-8,
                m=1e-8,
                L=10 ** 6,
                ploidy=p,
                model="hudson",
            )


class XiVsHudsonTest(Test):
    """
    Test that Xi dirac coalescent is equivalent to the Hudson model in the
    appropriate regime.
    """

    def _run(self, xi_model, num_replicates, num_samples, **kwargs):
        df = pd.DataFrame()
        for model in ["hudson", xi_model]:
            simulate_args = dict(kwargs)
            simulate_args["model"] = model
            model_str = "hudson"
            if model != "hudson":
                model_str = "Xi"
                # The Xi Dirac coalescent scales differently than the Hudson model.
                # (Ne² for Dirac and 2Ne for Hudson).
                # We need NeDirac= square_root(2NeHudson).
                simulate_args["population_size"] = math.sqrt(
                    int(simulate_args["ploidy"]) * int(simulate_args["population_size"])
                )
            logging.debug(f"Running: {simulate_args}")
            sim = msprime.ancestry._parse_sim_ancestry(
                samples=[msprime.SampleSet(num_samples, ploidy=1)],
                sequence_length=1,
                discrete_genome=False,
                **simulate_args,
            )
            replicates = sim.run_replicates(num_replicates)
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

        df_hudson = df[df.model == "hudson"]
        df_xi = df[df.model == "Xi"]
        p = int(simulate_args["ploidy"])
        for stat in ["tmrca_mean", "num_trees", "num_nodes", "num_edges"]:
            v1 = df_hudson[stat]
            v2 = df_xi[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            f = self.output_dir / f"{stat}_ploidy={p}.png"
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def test_xi_dirac_vs_hudson_recombination(self):
        self._run(
            msprime.DiracCoalescent(psi=0.99, c=0),
            num_replicates=1000,
            num_samples=50,
            population_size=10000,
            recombination_rate=0.001,
            ploidy=1,
        )
        self._run(
            msprime.DiracCoalescent(psi=0.99, c=0),
            num_replicates=1000,
            num_samples=50,
            population_size=10000,
            recombination_rate=0.001,
            ploidy=2,
        )

    def test_xi_dirac_vs_hudson_single_locus(self):
        self._run(
            msprime.DiracCoalescent(psi=0.99, c=0),
            num_replicates=5000,
            num_samples=10,
            population_size=10000,
            ploidy=1,
        )
        self._run(
            msprime.DiracCoalescent(psi=0.99, c=0),
            num_replicates=5000,
            num_samples=10,
            population_size=10000,
            ploidy=2,
        )


class KnownSFS(Test):
    """
    Compare the simulated SFS to precomputed known values.
    """

    def compare_sfs(self, sample_size, ploidy, model, num_replicates, sfs, name):
        data = collections.defaultdict(list)
        tbl_sum = [0] * (sample_size - 1)
        tot_bl_sum = [0]
        # Because we have input cases which sample_size % ploidy != 0 we can't
        # use the standard approach for specifying the samples. Sidestep this
        # by setting the initial state directly.
        tables = tskit.TableCollection(1)
        tables.populations.add_row()
        for _ in range(sample_size):
            tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0, population=0)
        replicates = msprime.sim_ancestry(
            initial_state=tables,
            ploidy=ploidy,
            model=model,
            num_replicates=num_replicates,
        )
        for ts in replicates:
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

        f = self.output_dir / f"{name}.png"
        ax = sns.violinplot(
            data=data, x="num_leaves", y="total_branch_length", color="grey"
        )
        ax.set_xlabel("num leaves")
        l1 = ax.plot(np.arange(sample_size - 1), sfs[::], ":", linewidth=3, marker="^")
        l2 = ax.plot(
            np.arange(sample_size - 1),
            [(x / num_replicates) / (tot_bl_sum[0] / num_replicates) for x in tbl_sum],
            "--",
            marker="o",
            linewidth=2,
        )
        ax.legend((l1[0], l2[0]), ("Expected", "Observed"))
        pyplot.savefig(f, dpi=72)
        pyplot.close("all")


class DiracSFS(KnownSFS):
    def _run(
        self, sample_size=10, ploidy=2, psi=None, c=None, sfs=None, num_replicates=10000
    ):
        """
        Runs simulations of the xi dirac model and calculates
        E[Bi]/E[B] (Bi branch length having i leaves and B total branch length)
        and compares to the expected SFS.
        """
        logging.debug(f"running SFS for {sample_size} {psi} {c}")
        model = msprime.DiracCoalescent(psi=psi, c=c)
        name = f"n={sample_size}_psi={psi}_c={c}_ploidy={ploidy}"
        self.compare_sfs(sample_size, ploidy, model, num_replicates, sfs, name)

    def test_xi_dirac_expected_sfs_psi_0_1_c_1(self):

        self._run(
            psi=0.1,
            c=1,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_psi_0_3_c_1(self):
        self._run(
            psi=0.3,
            c=1,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_psi_0_5_c_1(self):
        self._run(
            psi=0.5,
            c=1,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_psi_0_9_c_1(self):
        self._run(
            psi=0.9,
            c=1,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_n3(self):
        self._run(sample_size=3, ploidy=2, psi=0.1, c=10, sfs=[0.6667343, 0.3332657])
        self._run(sample_size=3, ploidy=2, psi=0.3, c=10, sfs=[0.6682113, 0.3317887])
        self._run(sample_size=3, ploidy=2, psi=0.5, c=10, sfs=[0.6721853, 0.3278147])
        self._run(sample_size=3, ploidy=2, psi=0.9, c=10, sfs=[0.6852703, 0.3147297])
        self._run(sample_size=3, ploidy=1, psi=0.1, c=10000, sfs=[0.678571, 0.321429])
        self._run(sample_size=3, ploidy=1, psi=0.3, c=10000, sfs=[0.708333, 0.291667])
        self._run(sample_size=3, ploidy=1, psi=0.5, c=10000, sfs=[0.750000, 0.250000])
        self._run(sample_size=3, ploidy=1, psi=0.9, c=10000, sfs=[0.916667, 0.083333])

    def test_xi_dirac_expected_sfs_psi_0_1_c_10(self):
        self._run(
            psi=0.1,
            c=10,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_psi_0_3_c_10(self):
        self._run(
            psi=0.3,
            c=10,
            ploidy=2,
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
        self._run(
            num_replicates=10000,
            sample_size=10,
            psi=0.5,
            c=10,
            ploidy=2,
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
        self._run(
            num_replicates=10000,
            sample_size=10,
            psi=0.9,
            c=10,
            ploidy=2,
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

    # Compare SFS when c=10000 to the expected SFS where c tends to infinity

    def test_xi_dirac_expected_sfs_psi_0_1_c_10000(self):
        self._run(
            psi=0.1,
            c=10000,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_psi_0_3_c_10000(self):
        self._run(
            psi=0.3,
            c=10000,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_psi_0_5_c_10000(self):
        self._run(
            psi=0.5,
            c=10000,
            ploidy=2,
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

    def test_xi_dirac_expected_sfs_psi_0_9_c_10000(self):
        self._run(
            psi=0.9,
            c=10000,
            ploidy=2,
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

    def test_dirac_expected_sfs_psi_0_1_c_10000(self):
        self._run(
            psi=0.1,
            c=10000,
            ploidy=1,
            sfs=[
                0.422312,
                0.148277,
                0.101947,
                0.077241,
                0.062498,
                0.052964,
                0.046659,
                0.043069,
                0.045033,
            ],
        )

    def test_dirac_expected_sfs_psi_0_3_c_10000(self):
        self._run(
            psi=0.3,
            c=10000,
            ploidy=1,
            sfs=[
                0.570300,
                0.083920,
                0.067942,
                0.056251,
                0.047302,
                0.041406,
                0.038521,
                0.039844,
                0.054512,
            ],
        )

    def test_dirac_expected_sfs_psi_0_5_c_10000(self):
        self._run(
            psi=0.5,
            c=10000,
            ploidy=1,
            sfs=[
                0.710037,
                0.036594,
                0.031667,
                0.031557,
                0.032135,
                0.031557,
                0.031667,
                0.036594,
                0.058192,
            ],
        )

    def test_dirac_expected_sfs_psi_0_9_c_10000(self):
        self._run(
            psi=0.9,
            c=10000,
            ploidy=1,
            sfs=[
                0.927920,
                0.001810,
                0.000476,
                0.000096,
                0.000148,
                0.001040,
                0.005356,
                0.018413,
                0.044742,
            ],
        )


class BetaSFS(KnownSFS):
    def _run(self, sample_size, ploidy, alpha, sfs, num_replicates=1000):
        """
        Runs simulations of the xi beta model and compares to the expected SFS.
        """
        logging.debug(f"running Beta SFS for {sample_size} {alpha}")
        model = msprime.BetaCoalescent(alpha=alpha)
        name = f"n={sample_size}_alpha={alpha}_ploidy={ploidy}"
        self.compare_sfs(sample_size, ploidy, model, num_replicates, sfs, name)

    def test_xi_beta_expected_sfs_alpha1_1(self):

        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.1,
            ploidy=2,
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

    def test_xi_beta_expected_sfs_alpha1_3(self):
        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.3,
            ploidy=2,
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

    def test_xi_beta_expected_sfs_alpha1_5(self):
        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.5,
            ploidy=2,
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

    def test_xi_beta_expected_sfs_alpha1_9(self):
        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.9,
            ploidy=2,
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

    def test_beta_expected_sfs_alpha1_1(self):
        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.1,
            ploidy=1,
            sfs=[
                0.580175,
                0.119103,
                0.066440,
                0.047197,
                0.038166,
                0.033879,
                0.032796,
                0.035382,
                0.046863,
            ],
        )

    def test_beta_expected_sfs_alpha1_3(self):
        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.3,
            ploidy=1,
            sfs=[
                0.521296,
                0.137166,
                0.078487,
                0.056070,
                0.045115,
                0.039481,
                0.037258,
                0.038479,
                0.046649,
            ],
        )

    def test_beta_expected_sfs_alpha1_5(self):
        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.5,
            ploidy=1,
            sfs=[
                0.467491,
                0.152216,
                0.090245,
                0.065103,
                0.052216,
                0.045067,
                0.041436,
                0.040898,
                0.045330,
            ],
        )

    def test_beta_expected_sfs_alpha1_9(self):
        self._run(
            num_replicates=100000,
            sample_size=10,
            alpha=1.9,
            ploidy=1,
            sfs=[
                0.374086,
                0.173264,
                0.112565,
                0.083644,
                0.066914,
                0.056165,
                0.048856,
                0.043826,
                0.040681,
            ],
        )


class XiGrowth(Test):
    def compare_tmrca(
        self, pop_size, growth_rate, model, num_replicates, a, b, ploidy, name
    ):
        demography = msprime.Demography.isolated_model(
            initial_size=[pop_size], growth_rate=[growth_rate]
        )

        replicates = msprime.ancestry.sim_ancestry(
            2,
            demography=demography,
            model=model,
            ploidy=ploidy,
            num_replicates=num_replicates,
        )
        T1 = np.array([ts.first().tmrca(0, 1) for ts in replicates])
        sm.graphics.qqplot(
            T1, dist=scipy.stats.gompertz, distargs=(a / b,), scale=1 / b, line="45"
        )
        filename = self.output_dir / f"{name}.png"
        pyplot.savefig(filename, dpi=72)
        pyplot.close("all")


class BetaGrowth(XiGrowth):
    def _run(self, pop_size, alpha, growth_rate, num_replicates=10000):
        logging.debug(f"running Beta growth for {pop_size} {alpha} {growth_rate}")
        b = growth_rate * (alpha - 1)
        model = (msprime.BetaCoalescent(alpha=alpha),)
        ploidy = 2
        a = 1 / (2 * ploidy * self.compute_beta_timescale(pop_size, alpha, ploidy))
        name = f"N={pop_size}_alpha={alpha}_growth_rate={growth_rate}_ploidy={ploidy}"
        self.compare_tmrca(
            pop_size, growth_rate, model, num_replicates, a, b, ploidy, name
        )
        ploidy = 1
        a = 1 / self.compute_beta_timescale(pop_size, alpha, ploidy)
        name = f"N={pop_size}_alpha={alpha}_growth_rate={growth_rate}_ploidy={ploidy}"
        self.compare_tmrca(
            pop_size, growth_rate, model, num_replicates, a, b, ploidy, name
        )

    def compute_beta_timescale(self, pop_size, alpha, ploidy):
        if ploidy > 1:
            N = pop_size / 2
            m = 2 + np.exp(
                alpha * np.log(2) + (1 - alpha) * np.log(3) - np.log(alpha - 1)
            )
        else:
            N = pop_size
            m = 1 + np.exp((1 - alpha) * np.log(2) - np.log(alpha - 1))
        ret = np.exp(
            alpha * np.log(m)
            + (alpha - 1) * np.log(N)
            - np.log(alpha)
            - scipy.special.betaln(2 - alpha, alpha)
        )
        return ret

    def test_10_15_01(self):
        self._run(pop_size=10, alpha=1.5, growth_rate=0.1)

    def test_1000_19_0001(self):
        self._run(pop_size=1000, alpha=1.9, growth_rate=0.001)

    def test_100000_11_001(self):
        self._run(pop_size=100000, alpha=1.1, growth_rate=0.01)


class DiracGrowth(XiGrowth):
    def _run(self, pop_size, c, psi, growth_rate, num_replicates=10000):
        logging.debug(f"running Dirac growth for {pop_size} {c} {psi} {growth_rate}")
        b = growth_rate
        model = (msprime.DiracCoalescent(psi=psi, c=c),)
        p = 2
        a = (1 + c * psi * psi / (2 * p)) / (pop_size * pop_size)
        name = f"N={pop_size}_c={c}_psi={psi}_growth_rate={growth_rate}_ploidy={p}"
        self.compare_tmrca(pop_size, growth_rate, model, num_replicates, a, b, p, name)
        p = 1
        a = (1 + c * psi * psi) / (pop_size * pop_size)
        name = f"N={pop_size}_c={c}_psi={psi}_growth_rate={growth_rate}_ploidy={p}"
        self.compare_tmrca(pop_size, growth_rate, model, num_replicates, a, b, p, name)

    def test_1_01_05_01(self):
        self._run(pop_size=1, c=0.1, psi=0.5, growth_rate=0.1)

    def test_10_05_07_0001(self):
        self._run(pop_size=10, c=0.5, psi=0.7, growth_rate=0.001)

    def test_100_1_09_001(self):
        self._run(pop_size=100, c=1, psi=0.9, growth_rate=0.01)

    def test_10_5_03_01(self):
        self._run(pop_size=10, c=5, psi=0.3, growth_rate=0.1)


class ContinuousVsDiscreteRecombination(Test):
    def _run_msprime_coalescent_stats(self, **kwargs):
        logging.debug(f"\t msprime: {kwargs}")
        if "num_replicates" in kwargs:
            replicates = kwargs["num_replicates"]
            num_trees = [0 for i in range(replicates)]
            breakpoints = [0 for i in range(replicates)]
            for i, ts in enumerate(msprime.sim_ancestry(**kwargs)):
                num_trees[i] = ts.num_trees
                breakpoints[i] = list(ts.breakpoints())
        else:
            ts = msprime.sim_ancestry(**kwargs)
            num_trees = [ts.num_trees]
            breakpoints = [list(ts.breakpoints)]

        d = {"num_trees": num_trees, "breakpoints": breakpoints}
        df = pd.DataFrame(d)
        return df

    def run_cont_discrete_comparison(self, model, recomb_map):
        sample_size = 10
        num_replicates = 400
        N = 100
        df_discrete = self._run_msprime_coalescent_stats(
            num_replicates=num_replicates,
            samples=sample_size,
            population_size=N,
            model=model,
            recombination_rate=recomb_map,
            discrete_genome=True,
        )
        df_cont = self._run_msprime_coalescent_stats(
            num_replicates=num_replicates,
            samples=sample_size,
            model=model,
            population_size=N,
            recombination_rate=recomb_map,
            discrete_genome=False,
        )
        self._plot_stats(
            "compare continuous and discrete coordinates",
            df_discrete,
            df_cont,
            "discrete",
            "continuous",
        )


class UniformRecombination(ContinuousVsDiscreteRecombination):
    def _run(self, model):
        recomb_map = msprime.RateMap.uniform(2000000, 1e-6)
        self.run_cont_discrete_comparison(model, recomb_map)

    def test_hudson_cont_discrete_uniform(self):
        self._run("hudson")

    def test_dtwf_cont_discrete_uniform(self):
        self._run("dtwf")


class VariableRecombination(ContinuousVsDiscreteRecombination):
    def _run(self, model):
        r = 1e-6
        positions = [0, 10000, 50000, 150000, 200000]
        rates = [0.0, r, 5 * r, r / 2]

        recomb_map = msprime.RateMap(positions, rates)

        self.run_cont_discrete_comparison(model, recomb_map)

    def test_hudson_cont_discrete_variable(self):
        self._run("hudson")

    def test_dtwf_cont_discrete_variable(self):
        self._run("dtwf")


class ArgRecordTest(Test):
    """
    Check that we get the same distributions of nodes and edges when
    we simplify an ARG as we get in a direct simulation.
    """

    def _run(self, num_replicates=1000, **kwargs):

        ts_node_counts = np.array([])
        arg_node_counts = np.array([])
        ts_tree_counts = np.array([])
        arg_tree_counts = np.array([])
        ts_edge_counts = np.array([])
        arg_edge_counts = np.array([])

        for ts in msprime.simulate(num_replicates=num_replicates, **kwargs):
            ts_node_counts = np.append(ts_node_counts, ts.num_nodes)
            ts_tree_counts = np.append(ts_tree_counts, ts.num_trees)
            ts_edge_counts = np.append(ts_edge_counts, ts.num_edges)

        reps = msprime.simulate(
            num_replicates=num_replicates, record_full_arg=True, **kwargs
        )
        for arg in reps:
            arg = arg.simplify()
            arg_node_counts = np.append(arg_node_counts, arg.num_nodes)
            arg_tree_counts = np.append(arg_tree_counts, arg.num_trees)
            arg_edge_counts = np.append(arg_edge_counts, arg.num_edges)

        pp_ts = sm.ProbPlot(ts_node_counts)
        pp_arg = sm.ProbPlot(arg_node_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        pyplot.savefig(self.output_dir / "nodes.png", dpi=72)

        pp_ts = sm.ProbPlot(ts_tree_counts)
        pp_arg = sm.ProbPlot(arg_tree_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        pyplot.savefig(self.output_dir / "num_trees.png", dpi=72)

        pp_ts = sm.ProbPlot(ts_edge_counts)
        pp_arg = sm.ProbPlot(arg_edge_counts)
        sm.qqplot_2samples(pp_ts, pp_arg, line="45")
        pyplot.savefig(self.output_dir / "edges.png", dpi=72)
        pyplot.close("all")

    def test_arg_hudson_n10_rho_20(self):
        self._run(sample_size=10, recombination_rate=20)

    def test_arg_hudson_n1000_rho_0_2(self):
        self._run(sample_size=1000, recombination_rate=0.2)

    def test_arg_beta_n100_rho_2(self):
        model = msprime.BetaCoalescent(alpha=1.1)
        self._run(sample_size=100, recombination_rate=2, model=model)

    def test_arg_dirac_n100_rho_2(self):
        model = msprime.DiracCoalescent(psi=0.9, c=1)
        self._run(sample_size=100, recombination_rate=2, model=model)


class HudsonAnalytical(Test):
    """
    Miscellaneous tests for the hudson model where we verify against
    analytical results.
    """

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

    def test_analytical_segsites(self):
        """
        Runs the check for the number of segregating sites against the
        analytical prediction. We also compare against ms.
        """
        R = 100000
        theta = 2
        for n in range(2, 15):
            logging.debug(f"Running n = {n}")
            cmd = f"{n} {R} -t {theta}"
            S_ms = self.get_segregating_sites_histogram(
                _ms_executable + cmd.split() + self.get_ms_seeds()
            )
            S_msp = self.get_segregating_sites_histogram(
                _mspms_executable + cmd.split() + self.get_ms_seeds()
            )

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
            pyplot.savefig(self.output_dir / f"{n}.png")

    def test_analytical_pi(self):
        """
        Runs the check for pi against analytical predictions.
        """
        R = 100000
        theta = 4.5

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
            # Predicted mean is theta.
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

        filename = self.output_dir / "mean.png"
        pyplot.plot(sample_size, predicted_mean, "-")
        pyplot.plot(sample_size, mean, "-")
        pyplot.savefig(filename)
        pyplot.close("all")

        filename = self.output_dir / "var.png"
        pyplot.plot(sample_size, predicted_var, "-")
        pyplot.plot(sample_size, var, "-")
        pyplot.savefig(filename)
        pyplot.close("all")

    def test_gc_correlation_between_trees(self):
        """
        Runs the check for the probability of same tree at two sites against
        analytical predictions.
        """
        R = 1000
        sample_size = 1  # 2 diploids
        gc_length_rate_ratio = np.array([0.05, 0.5, 5.0])
        gc_length = np.array([100, 50, 20])
        gc_rate = 0.25 / (gc_length_rate_ratio * gc_length)
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
            replicates = msprime.sim_ancestry(
                samples=sample_size,
                sequence_length=seq_length,
                gene_conversion_rate=gc_rate[k],
                gene_conversion_tract_length=gc_length[k],
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
        pyplot.plot(x, predicted_prob[0], "--", label="prediction")
        pyplot.plot(x, empirical_prob_first[0], "-", label="simulation")
        pyplot.plot(x, predicted_prob[1], "--")
        pyplot.plot(x, empirical_prob_first[1], "-")
        pyplot.plot(x, predicted_prob[2], "--")
        pyplot.plot(x, empirical_prob_first[2], "-")
        pyplot.xlabel("chromosome positon")
        pyplot.ylabel("fraction of trees identical to first position tree")
        pyplot.legend(loc="upper right")
        pyplot.savefig(self.output_dir / "prob_first.png")
        pyplot.close("all")

        pyplot.plot(x, predicted_prob[0, ::-1], "--", label="prediction")
        pyplot.plot(x, empirical_prob_last[0], "-", label="simulation")
        pyplot.plot(x, predicted_prob[1, ::-1], "--")
        pyplot.plot(x, empirical_prob_last[1], "-")
        pyplot.plot(x, predicted_prob[2, ::-1], "--")
        pyplot.plot(x, empirical_prob_last[2], "-")
        pyplot.xlabel("chromosome positon")
        pyplot.ylabel("fraction of trees identical to last position tree")
        pyplot.legend(loc="upper left")
        pyplot.savefig(self.output_dir / "prob_last.png")
        pyplot.close("all")

        pyplot.plot(
            x,
            np.concatenate((predicted_prob[0, 249::-1], predicted_prob[0, :250])),
            "--",
            label="prediction",
        )
        pyplot.plot(x, empirical_prob_mid[0], "-", label="simulation")
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
        pyplot.xlabel("chromosome positon")
        pyplot.ylabel("fraction of trees identical to middle position tree")
        pyplot.legend(loc="upper right")
        pyplot.savefig(self.output_dir / "prob_mid.png")
        pyplot.close("all")

        x = np.arange(10) + 1
        pyplot.plot(x, predicted_prob[0, range(10)], "--", label="prediction")
        pyplot.plot(x, empirical_prob_first[0, range(10)], "-", label="simulation")
        pyplot.plot(x, predicted_prob[1, range(10)], "--")
        pyplot.plot(x, empirical_prob_first[1, range(10)], "-")
        pyplot.plot(x, predicted_prob[2, range(10)], "--")
        pyplot.plot(x, empirical_prob_first[2, range(10)], "-")
        pyplot.xlabel("chromosome positon")
        pyplot.ylabel("fraction of trees identical to first position tree")
        pyplot.legend(loc="upper right")
        pyplot.savefig(self.output_dir / "prob_first_zoom.png")
        pyplot.close("all")

    def test_gc_tract_length_expectation(self):
        """
        Runs the check for the mean length of gene conversion tracts.
        """
        num_replicates = 100
        n = 10
        gene_conversion_rate = 5
        gc_tract_lengths = np.append(np.arange(1, 5.25, 0.25), [10, 50])

        for discrete_genome in [True, False]:
            data_to_plot = []

            for k, l in enumerate(gc_tract_lengths):
                num_gc_events = np.zeros(num_replicates)
                num_internal_gc_events = np.zeros(num_replicates)
                sum_internal_gc_tract_lengths = np.zeros(num_replicates)

                sim = msprime.ancestry._parse_sim_ancestry(
                    samples=n,
                    sequence_length=100,
                    gene_conversion_rate=gene_conversion_rate,
                    gene_conversion_tract_length=gc_tract_lengths[k],
                    discrete_genome=discrete_genome,
                    ploidy=1,
                )
                for j, _ts in enumerate(sim.run_replicates(num_replicates)):
                    num_gc_events[j] = sim.num_gene_conversion_events
                    num_internal_gc_events[j] = sim.num_internal_gene_conversion_events
                    sum_internal_gc_tract_lengths[j] = sim.sum_internal_gc_tract_lengths
                    sim.reset()
                data_to_plot.append(
                    sum_internal_gc_tract_lengths / num_internal_gc_events / l
                )
            pyplot.boxplot(data_to_plot, labels=gc_tract_lengths)
            pyplot.xlabel("tl: mean tract length specified")
            pyplot.ylabel("average internal tract length / tl")
            filename = f"mean_gc_tract_lengths_discrete={int(discrete_genome)}.png"
            pyplot.savefig(self.output_dir / filename)
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

    def test_analytical_tbl(self):
        """
        Runs the check for the total branch length.
        """
        R = 10000
        for n in range(2, 15):
            logging.debug(f"Running for n = {n}")
            tbl_ms = self.get_tbl_distribution(n, R, _ms_executable)
            tbl_msp = self.get_tbl_distribution(n, R, _mspms_executable)

            sm.graphics.qqplot(tbl_ms)
            sm.qqplot_2samples(tbl_ms, tbl_msp, line="45")
            pyplot.savefig(self.output_dir / f"qqplot_{n}.png", dpi=72)
            pyplot.close("all")

            hist_ms, bin_edges = np.histogram(tbl_ms, 20, density=True)
            hist_msp, _ = np.histogram(tbl_msp, bin_edges, density=True)

            index = bin_edges[:-1]
            # NOTE We don't to have the analytical value quite right here,
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
            pyplot.savefig(self.output_dir / f"hist_{n}.png")

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

    def test_analytical_num_trees(self):
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
                _ms_executable + cmd.split() + self.get_ms_seeds(), num_replicates
            )
            ms_mean[j] = np.mean(T)

            T = self.get_num_trees(
                _mspms_executable + cmd.split() + self.get_ms_seeds(), num_replicates
            )
            msp_mean[j] = np.mean(T)
        pyplot.plot(rho, ms_mean, "o")
        pyplot.plot(rho, msp_mean, "^")
        pyplot.plot(rho, rho * harmonic_number(n - 1), "-")
        pyplot.savefig(self.output_dir / "mean.png")
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

    def test_analytical_pairwise_island_model(self):
        """
        Runs the check for the pairwise coalscence times for within
        and between populations.
        """
        R = 10000
        M = 0.2

        for d in range(2, 6):
            cmd = "2 {} -T -I {} 2 {} {}".format(R, d, "0 " * (d - 1), M)
            T_w_ms = self.get_pairwise_coalescence_time(
                _ms_executable + cmd.split() + self.get_ms_seeds(), R
            )
            T_w_msp = self.get_pairwise_coalescence_time(
                _mspms_executable + cmd.split() + self.get_ms_seeds(), R
            )

            cmd = "2 {} -T -I {} 1 1 {} {}".format(R, d, "0 " * (d - 2), M)
            T_b_ms = self.get_pairwise_coalescence_time(
                _ms_executable + cmd.split() + self.get_ms_seeds(), R
            )
            T_b_msp = self.get_pairwise_coalescence_time(
                _mspms_executable + cmd.split() + self.get_ms_seeds(), R
            )
            t_within = d / 2
            t_between = (d + (d - 1) / M) / 2
            logging.debug(
                f"d={d} within=({np.mean(T_w_msp):.2f},{t_within}) "
                f"between=({np.mean(T_b_msp):.2f}, {t_between})"
            )

            sm.graphics.qqplot(T_w_ms)
            sm.qqplot_2samples(T_w_ms, T_w_msp, line="45")
            pyplot.savefig(self.output_dir / f"within_{d}.png", dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(T_b_ms)
            sm.qqplot_2samples(T_b_ms, T_b_msp, line="45")
            pyplot.savefig(self.output_dir / f"between_{d}.png", dpi=72)
            pyplot.close("all")


class DemographyDebugger(Test):
    """
    Tests for the demography debugger methods.
    """

    def verify_ddb_mean_coaltime(self, model_factory, name):
        """
        Checks the mean coalescence time calculation against pi.
        """
        num_reps = 20
        T = []
        U = []
        logging.debug("coaltime: theory  mean  sd   z")
        for k, model in enumerate(model_factory()):
            ddb = msprime.DemographyDebugger(
                population_configurations=model["population_configurations"],
                demographic_events=model["demographic_events"],
                migration_matrix=model["migration_matrix"],
            )
            u = ddb.mean_coalescence_time(num_samples=model["sample_size"], max_iter=18)
            U.append(u)

            mut_rate = 1e-7
            replicates = msprime.simulate(
                length=1e7,
                recombination_rate=1e-8,
                mutation_rate=mut_rate,
                population_configurations=model["population_configurations"],
                demographic_events=model["demographic_events"],
                migration_matrix=model["migration_matrix"],
                random_seed=5 + k,
                num_replicates=num_reps,
            )
            TT = np.zeros(num_reps)
            for j, ts in enumerate(replicates):
                TT[j] = ts.diversity(ts.samples())
                TT[j] /= 2 * mut_rate
            T.append(TT)
            mT = np.mean(TT)
            sT = np.std(TT)
            logging.debug(
                "        {:.2f} {:.2f} {:.2f} {:.2f}".format(
                    u, mT, sT, (u - mT) / (sT / np.sqrt(num_reps))
                )
            )

        U = np.array(U)
        T = np.array(T)
        fig, ax = pyplot.subplots()
        ax.scatter(np.column_stack([U] * T.shape[1]), T)
        ax.scatter(U, np.mean(T, axis=1))
        # where oh where is abline(0,1)
        x_vals = np.array(ax.get_xlim())
        ax.plot(x_vals, x_vals, "--")
        ax.set_xlabel("calculated mean coaltime")
        ax.set_ylabel("pairwise diversity, scaled")
        pyplot.savefig(self.output_dir / f"{name}_mean_coaltimes.png")
        pyplot.close("all")

    def random_model_factory(self):
        """
        Checks the mean coalescence time calculation against pi.
        """
        random.seed(5)
        num_models = 20
        for _ in range(num_models):
            Ne = 100
            npops = 4
            pop_sizes = [random.uniform(0.1, 1) * Ne for _ in range(npops)]
            growth_rates = [random.uniform(-0.001, 0.001) for _ in range(npops)]
            migration_matrix = [
                [random.random() * (i != j) for j in range(npops)] for i in range(npops)
            ]
            sample_size = [random.randint(2, 10) for _ in range(npops)]
            population_configurations = [
                msprime.PopulationConfiguration(
                    initial_size=j, sample_size=n, growth_rate=r
                )
                for j, n, r in zip(pop_sizes, sample_size, growth_rates)
            ]
            demographic_events = []
            for i in [0, 1]:
                n = random.uniform(0.1, 10) * Ne
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
                n = random.uniform(0.1, 10) * Ne
                r = random.uniform(-0.01, 0.01)
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=210, initial_size=n, growth_rate=r, population_id=i
                    )
                )
            for i in [1, 2, 3]:
                n = random.uniform(0.1, 10) * Ne
                r = random.uniform(0.0, 0.01)
                demographic_events.append(
                    msprime.PopulationParametersChange(
                        time=250, initial_size=n, growth_rate=r, population_id=i
                    )
                )
            yield {
                "population_configurations": population_configurations,
                "demographic_events": demographic_events,
                "migration_matrix": migration_matrix,
                "sample_size": sample_size,
            }

    def migration_model_factory(self):
        random.seed(5)
        Ne = 100
        npops = 3
        num_models = 10
        for k in range(num_models):
            pop_sizes = [Ne] * (npops - 1) + [Ne * (2 ** k)]
            migration_matrix = [
                [2 ** (k - 4) * ((i - j) % npops == 1) / Ne for j in range(npops)]
                for i in range(npops)
            ]
            sample_size = [1 + j for j in range(npops)]
            population_configurations = [
                msprime.PopulationConfiguration(initial_size=j, sample_size=n)
                for j, n in zip(pop_sizes, sample_size)
            ]
            demographic_events = []
            yield {
                "population_configurations": population_configurations,
                "demographic_events": demographic_events,
                "migration_matrix": migration_matrix,
                "sample_size": sample_size,
            }

    def popsize_change_model_factory(self):
        random.seed(5)
        Ne = 100
        npops = 3
        num_models = 16
        change_times = [j * Ne / 4 for j in range(8)]
        for k in range(num_models):
            pop_sizes = [Ne] * (npops - 1) + [Ne * (2 ** k)]
            migration_matrix = [
                [10 * ((i - j) % npops == 1) / Ne for j in range(npops)]
                for i in range(npops)
            ]
            sample_size = [1 + j for j in range(npops)]
            population_configurations = [
                msprime.PopulationConfiguration(initial_size=j, sample_size=n)
                for j, n in zip(pop_sizes, sample_size)
            ]
            demographic_events = []
            for t in change_times:
                pop_sizes = pop_sizes[1:] + pop_sizes[:1]
                r = 0
                for i, n in enumerate(pop_sizes):
                    demographic_events.append(
                        msprime.PopulationParametersChange(
                            time=t, initial_size=n, growth_rate=r, population_id=i
                        )
                    )
            yield {
                "population_configurations": population_configurations,
                "demographic_events": demographic_events,
                "migration_matrix": migration_matrix,
                "sample_size": sample_size,
            }

    def test_random_mean_coaltime(self):
        """
        Checks the mean coalescence time calculation against pi.
        """
        self.verify_ddb_mean_coaltime(self.random_model_factory, "random")

    def test_popsize_change_mean_coaltime(self):
        """
        Checks the mean coalescence time calculation against pi for some models
        with population size changes.
        """
        self.verify_ddb_mean_coaltime(
            self.popsize_change_model_factory, "popsize_change"
        )

    def test_migration_mean_coaltime(self):
        """
        Checks the mean coalescence time calculation against pi
        for some models with migration.
        """
        self.verify_ddb_mean_coaltime(self.migration_model_factory, "migration")


class SmcTest(Test):
    """
    Tests for the SMC model against scrm.
    """

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

    def test_smc_oldest_time(self):
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
                _scrm_executable + cmd.split() + self.get_ms_seeds(), num_replicates
            )
            scrm_mean[j] = np.mean(T)

            cmd += " -l 0"
            T = self.get_scrm_oldest_time(
                _scrm_executable + cmd.split() + self.get_ms_seeds(), num_replicates
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
        pyplot.plot(rho, scrm_mean, "-", color="blue", label="scrm")
        pyplot.plot(rho, scrm_smc_mean, "-", color="red", label="scrm_smc")
        pyplot.plot(rho, msp_smc_mean, "--", color="red", label="msprime_smc")
        pyplot.plot(rho, msp_mean, "--", color="blue", label="msprime")
        pyplot.xlabel("rho")
        pyplot.ylabel("Mean oldest coalescence time")
        pyplot.legend(loc="lower right")
        pyplot.savefig(self.output_dir / "mean.png")
        pyplot.close("all")

    def test_smc_num_trees(self):
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
                _scrm_executable + cmd.split() + self.get_ms_seeds(), num_replicates
            )
            mean_scrm[j] = np.mean(T)
            var_scrm[j] = np.var(T)
            # IMPORTANT!! We have to use the get_num_breakpoints method
            # on the simulator as there is a significant drop in the number
            # of trees if we use the tree sequence. There is a significant
            # number of common ancestor events that result in a recombination
            # being undone.
            exact_sim = msprime.ancestry._parse_simulate(
                sample_size=n, recombination_rate=r, Ne=Ne, length=L[j]
            )
            for k in range(num_replicates):
                exact_sim.run()
                num_trees[k] = exact_sim.num_breakpoints
                exact_sim.reset()
            mean_exact[j] = np.mean(num_trees)
            var_exact[j] = np.var(num_trees)

            smc_sim = msprime.ancestry._parse_simulate(
                sample_size=n, recombination_rate=r, Ne=Ne, length=L[j], model="smc"
            )
            for k in range(num_replicates):
                smc_sim.run()
                num_trees[k] = smc_sim.num_breakpoints
                smc_sim.reset()
            mean_smc[j] = np.mean(num_trees)
            var_smc[j] = np.var(num_trees)

            smc_prime_sim = msprime.ancestry._parse_simulate(
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

        pyplot.plot(rho, mean_exact, "o", label="msprime (hudson)")
        pyplot.plot(rho, mean_smc, "^", label="msprime (smc)")
        pyplot.plot(rho, mean_smc_prime, "*", label="msprime (smc_prime)")
        pyplot.plot(rho, mean_scrm, "x", label="scrm")
        pyplot.plot(rho, rho * harmonic_number(n - 1), "-")
        pyplot.legend(loc="upper left")
        pyplot.xlabel("scaled recombination rate rho")
        pyplot.ylabel("Mean number of breakpoints")
        pyplot.savefig(self.output_dir / "mean.png")
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
        pyplot.savefig(self.output_dir / "var.png")
        pyplot.close("all")


class SimulateFrom(Test):
    def test_simulate_from_single_locus(self):
        num_replicates = 1000

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
                filename = self.output_dir / f"T_mrca_n={n}_t={t}.png"
                pyplot.savefig(filename, dpi=72)
                pyplot.close("all")

    def test_simulate_from_multi_locus(self):
        num_replicates = 1000
        n = 50

        for m in [10, 50, 100, 1000]:
            logging.debug(f"running for m = {m}")
            T1 = np.zeros(num_replicates)
            num_trees1 = np.zeros(num_replicates)
            recomb_rate = 1 / m
            reps = msprime.sim_ancestry(
                n,
                recombination_rate=recomb_rate,
                sequence_length=m,
                num_replicates=num_replicates,
            )
            for j, ts in enumerate(reps):
                T1[j] = np.max(ts.tables.nodes.time)
                num_trees1[j] = ts.num_trees

            for t in [0.5, 1, 1.5, 5]:
                T2 = np.zeros(num_replicates)
                num_trees2 = np.zeros(num_replicates)
                reps = msprime.sim_ancestry(
                    n,
                    num_replicates=num_replicates,
                    recombination_rate=recomb_rate,
                    sequence_length=m,
                    end_time=t,
                )
                for j, ts in enumerate(reps):
                    final_ts = msprime.sim_ancestry(
                        initial_state=ts,
                        recombination_rate=recomb_rate,
                        start_time=np.max(ts.tables.nodes.time),
                    )
                    final_ts = final_ts.simplify()
                    T2[j] = np.max(final_ts.tables.nodes.time)
                    num_trees2[j] = final_ts.num_trees

                sm.graphics.qqplot(T1)
                sm.qqplot_2samples(T1, T2, line="45")
                filename = self.output_dir / f"T_mrca_m={m}_t={t}.png"
                pyplot.savefig(filename, dpi=72)
                pyplot.close("all")

                sm.graphics.qqplot(num_trees1)
                sm.qqplot_2samples(num_trees1, num_trees2, line="45")
                filename = self.output_dir / f"num_trees_m={m}_t={t}.png"
                pyplot.savefig(filename, dpi=72)
                pyplot.close("all")

    def test_simulate_from_recombination(self):
        num_replicates = 1000
        n = 100
        recombination_rate = 10

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
            "original mean: trees={:.2f} nodes={:.2f} edges={:.2f}".format(
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
                "t = {} mean: trees={:.2f} nodes={:.2f} edges={:.2f}".format(
                    t, np.mean(num_trees2), np.mean(num_nodes2), np.mean(num_edges2)
                )
            )

            sm.graphics.qqplot(T1)
            sm.qqplot_2samples(T1, T2, line="45")
            filename = self.output_dir / f"T_mrca_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_trees1)
            sm.qqplot_2samples(num_trees1, num_trees2, line="45")
            filename = self.output_dir / f"num_trees_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_edges1)
            sm.qqplot_2samples(num_edges1, num_edges2, line="45")
            filename = self.output_dir / f"num_edges_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_nodes1)
            sm.qqplot_2samples(num_nodes1, num_nodes2, line="45")
            filename = self.output_dir / f"num_nodes_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

    def test_simulate_from_demography(self):
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

        T1 = np.zeros(num_replicates)
        num_ca_events1 = np.zeros(num_replicates)
        num_re_events1 = np.zeros(num_replicates)
        num_mig_events1 = np.zeros(num_replicates)
        num_trees1 = np.zeros(num_replicates)
        num_edges1 = np.zeros(num_replicates)
        num_nodes1 = np.zeros(num_replicates)

        sim = msprime.ancestry._parse_simulate(
            samples=samples,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events,
            recombination_rate=recombination_rate,
        )
        logging.debug("t\ttrees\tnodes\tedges\tca\tre\tmig")
        for j, ts in enumerate(sim.run_replicates(num_replicates)):
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
            "{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
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
            sim = msprime.ancestry._parse_simulate(
                samples=samples,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix,
                demographic_events=demographic_events,
                recombination_rate=recombination_rate,
                end_time=t,
            )
            for j, ts in enumerate(sim.run_replicates(num_replicates)):
                num_ca_events2[j] = sim.num_common_ancestor_events
                num_re_events2[j] = sim.num_recombination_events
                num_mig_events2[j] = sum(
                    [r for row in sim.num_migration_events for r in row]
                )

                max_time = max(node.time for node in ts.nodes())
                sim2 = msprime.ancestry._parse_simulate(
                    from_ts=ts,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=[
                        e for e in demographic_events if e.time > max_time
                    ],
                    recombination_rate=recombination_rate,
                )
                final_ts = next(sim2.run_replicates(1)).simplify()

                num_ca_events2[j] += sim2.num_common_ancestor_events
                num_re_events2[j] += sim2.num_recombination_events
                num_mig_events2[j] += sum(
                    [r for row in sim2.num_migration_events for r in row]
                )

                T2[j] = np.max(final_ts.tables.nodes.time)
                num_trees2[j] = final_ts.num_trees
                num_nodes2[j] = final_ts.num_nodes
                num_edges2[j] = final_ts.num_edges
                sim.reset()

            logging.debug(
                "{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
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
            filename = self.output_dir / f"T_mrca_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_trees1)
            sm.qqplot_2samples(num_trees1, num_trees2, line="45")
            filename = self.output_dir / f"num_trees_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_edges1)
            sm.qqplot_2samples(num_edges1, num_edges2, line="45")
            filename = self.output_dir / f"num_edges_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_nodes1)
            sm.qqplot_2samples(num_nodes1, num_nodes2, line="45")
            filename = self.output_dir / f"num_nodes_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_ca_events1)
            sm.qqplot_2samples(num_ca_events1, num_ca_events2, line="45")
            filename = self.output_dir / f"num_ca_events_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_re_events1)
            sm.qqplot_2samples(num_re_events1, num_re_events2, line="45")
            filename = self.output_dir / f"num_re_events_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")

            sm.graphics.qqplot(num_mig_events1)
            sm.qqplot_2samples(num_mig_events1, num_mig_events2, line="45")
            filename = self.output_dir / f"num_mig_events_t={t}.png"
            pyplot.savefig(filename, dpi=72)
            pyplot.close("all")


class MutationStatsTest(Test):
    def plot_relative_error(self, x_values, observed, expected, name):
        x = np.array([np.full(o.shape, xv) for xv, o in zip(x_values, observed)])
        observed = np.array(observed)
        expected = np.array(expected)
        outfile = self._build_filename(None, name)
        if not np.all(observed[expected == 0] == 0):
            raise ValueError("Impossible mutations occurred!")
        nonzero = expected > 0
        rel_err = (observed[nonzero] - expected[nonzero]) / expected[nonzero]
        unique_x = np.unique(x)
        x_index = np.searchsorted(unique_x, x[nonzero])
        mean_rel_err = np.bincount(x_index, weights=np.abs(rel_err))
        mean_rel_err /= np.bincount(x_index)
        n_expected = np.repeat(-1, len(unique_x))
        for j, exp in enumerate(expected):
            n_expected[j] = 1 / np.mean(1 / exp[exp > 0])

        fig, (ax1, ax2) = pyplot.subplots(1, 2, figsize=(12, 6))
        ax1.scatter(x[nonzero], rel_err)
        ax1.plot([0, max(unique_x)], [0, 0], linestyle=":")
        ax1.set_xlabel("sample size")
        ax1.set_ylabel("relative error")
        ax2.plot(unique_x, mean_rel_err, label="mean relative error")
        ax2.plot(
            unique_x,
            1 / np.sqrt(1 + n_expected),
            linestyle=":",
            label="rough expected behaviour",
        )
        ax2.plot([0, max(unique_x)], [0, 0], linestyle=":")
        ax2.set_ylim(0, max(0.001, max(mean_rel_err)))
        ax2.set_xlabel("sample size")
        ax2.set_ylabel("mean relative error")
        ax2.legend()
        pyplot.savefig(outfile, dpi=72)
        pyplot.close(fig)

    def plot_uniform(self, x, name):
        outfile = self._build_filename(None, name)
        x = np.array(sorted(x))
        fig, ax = pyplot.subplots(1, 1, figsize=(8, 8))
        ax.scatter(np.linspace(0, 1, len(x)), x)
        ax.plot([-0.1, 1.1], [-0.1, 1.1], "r-", linewidth=2)
        ax.set_xlabel("expected relative mutation spacings")
        ax.set_ylabel("observed relative mutation spacings")
        pyplot.savefig(outfile, dpi=72)
        pyplot.close(fig)

    def plot_y_equals_x(self, x, y, name):
        x = np.array(x).flatten()
        y = np.array(y).flatten()
        xx = np.linspace(1, 1.1 * max(x), 51)
        outfile = self._build_filename(None, name)
        fig, ax = pyplot.subplots(1, 1, figsize=(8, 8))
        ax.scatter(x, y)
        ax.plot(
            [0, 1.1 * np.max(x)], [0, 1.1 * np.max(x)], "r-", linewidth=2, label="y = x"
        )
        ax.plot(
            xx, xx + 4 * np.sqrt(xx), "r:", linewidth=2, label="rough expected bounds"
        )
        ax.plot(xx, xx - 4 * np.sqrt(xx), "r:", linewidth=2)
        ax.legend()
        ax.set_xlabel("expected")
        ax.set_ylabel("observed")
        pyplot.savefig(outfile, dpi=72)
        pyplot.close(fig)

    def verify_model(self, model, name, verify_rates=False, state_independent=False):
        L = 100000
        ots = msprime.sim_ancestry(
            8,
            random_seed=7,
            recombination_rate=3 / L,
            sequence_length=L,
            discrete_genome=True,
        )
        for discrete_genome in (True, False):
            verify_times = (not discrete_genome) or state_independent
            x = []
            observed = []
            expected = []
            observed_roots = []
            expected_roots = []
            observed_rates = []
            expected_rates = []
            observed_times = []
            for nmuts in (100, 500, 1000, 1500, 2500, 3500, 5000, 7500, 10000, 20000):
                rate = nmuts / L
                ts = msprime.sim_mutations(
                    ots,
                    random_seed=6,
                    rate=rate,
                    model=model,
                    discrete_genome=discrete_genome,
                )
                x.append(nmuts)
                # transitions
                obs, exp = self.verify_transitions(ts, model, discrete_genome, rate)
                observed.append(obs)
                expected.append(exp)
                # root distributions
                obs_roots, exp_roots = self.verify_roots(
                    ts, model, discrete_genome, rate
                )
                observed_roots.append(obs_roots)
                expected_roots.append(exp_roots)
                # mutation rates
                obs_rates, exp_rates = self.verify_mutation_rates(
                    ts, model, rate, discrete_genome
                )
                observed_rates.append(obs_rates)
                expected_rates.append(exp_rates)
                if verify_times:
                    obs_times = self.verify_mutation_times(ts)
                    observed_times.extend(obs_times)
            if discrete_genome:
                pname = f"{name}_discrete"
            else:
                pname = f"{name}_continuous"
            if name != "binary":
                self.plot_relative_error(
                    x, observed=observed, expected=expected, name=pname + "_transitions"
                )
            self.plot_relative_error(
                x,
                observed=observed_roots,
                expected=expected_roots,
                name=pname + "_roots",
            )
            # check mutation times
            if verify_times:
                self.plot_uniform(observed_times, name=pname + "_times")
            if verify_rates:
                # this test only works if the probability of dropping a mutation
                # doesn't depend on the previous state
                assert len(set(np.diag(model.transition_matrix))) == 1
                self.plot_y_equals_x(
                    observed_rates, expected_rates, name=pname + "_rates"
                )

    def verify_transitions(self, ts, model, discrete_genome, mutation_rate):
        alleles = model.alleles
        num_alleles = len(alleles)

        observed = np.zeros((num_alleles, num_alleles))
        expected = np.zeros((num_alleles, num_alleles))
        for mut in ts.mutations():
            if mut.parent == tskit.NULL:
                pa = ts.site(mut.site).ancestral_state
            else:
                pa = ts.mutation(mut.parent).derived_state
            da = mut.derived_state
            observed[alleles.index(pa), alleles.index(da)] += 1
        for j, (row, p) in enumerate(zip(observed, model.transition_matrix)):
            p[j] = 0
            p /= sum(p)
            expected[j, :] = sum(row) * p
        return observed, expected

    def verify_roots(self, ts, model, discrete_genome, mutation_rate):
        alleles = model.alleles
        num_alleles = len(alleles)
        observed = np.zeros((num_alleles,))

        for site in ts.sites():
            aa = site.ancestral_state
            observed[alleles.index(aa)] += 1

        expected = np.zeros(num_alleles)
        change_probs = model.transition_matrix.sum(axis=1) - np.diag(
            model.transition_matrix
        )
        for t in ts.trees():
            if discrete_genome:
                t_span = np.ceil(t.interval[1] - np.ceil(t.interval[0]))
                expected += (
                    model.root_distribution
                    * t_span
                    * (
                        1
                        - np.exp(-mutation_rate * t.total_branch_length * change_probs)
                    )
                )
            else:
                t_span = t.span
                expected += (
                    model.root_distribution
                    * mutation_rate
                    * t.total_branch_length
                    * t_span
                    * change_probs
                )

        return observed, expected

    def verify_mutation_rates(self, ts, model, rate, discrete_genome):
        mut_rate = rate * (1 - model.transition_matrix[0, 0])
        observed = np.zeros(ts.num_trees)
        expected = np.zeros(ts.num_trees)
        for j, t in enumerate(ts.trees()):
            if discrete_genome:
                span = np.ceil(t.interval[1]) - np.ceil(t.interval[0])
            else:
                span = t.span
            mean = mut_rate * span * t.total_branch_length
            observed[j] = t.num_mutations
            # if we draw an indepenent Poisson with the same mean
            # it should be greater than observed half the time it is different
            expected[j] = scipy.stats.poisson.rvs(mean, 1)
        return observed, expected

    def verify_mutation_times(self, ts):
        start_time = np.full(ts.num_mutations, -1, dtype=np.float32)
        end_time = np.full(ts.num_mutations, -1, dtype=np.float32)
        mut_time = np.full(ts.num_mutations, -1, dtype=np.float32)
        for t in ts.trees():
            for mut in t.mutations():
                mut_time[mut.id] = mut.time
                end_time[mut.id] = t.time(t.parent(mut.node))
                start_time[mut.id] = t.time(mut.node)
                if mut.parent != tskit.NULL:
                    end_time[mut.id] = min(
                        end_time[mut.id], ts.mutation(mut.parent).time
                    )
                    start_time[mut.parent] = max(start_time[mut.parent], mut.time)
        return (mut_time - start_time) / (end_time - start_time)

    def test_binary_model_stats(self):
        model = msprime.BinaryMutationModel()
        self.verify_model(
            model, name="binary", state_independent=True, verify_rates=True
        )

    def test_jukes_cantor_stats(self):
        model = msprime.JC69MutationModel()
        self.verify_model(
            model, name="jukes_cantor", state_independent=True, verify_rates=True
        )

    def test_HKY_stats(self):
        equilibrium_frequencies = [0.3, 0.2, 0.3, 0.2]
        kappa = 0.75
        model = msprime.HKYMutationModel(
            kappa=kappa, equilibrium_frequencies=equilibrium_frequencies
        )
        self.verify_model(model, name="HKY")

    def test_F84_stats(self):
        equilibrium_frequencies = [0.4, 0.1, 0.1, 0.4]
        kappa = 0.75
        model = msprime.F84MutationModel(
            kappa=kappa, equilibrium_frequencies=equilibrium_frequencies
        )
        self.verify_model(model, name="F84")

    def test_GTR_stats(self):
        relative_rates = [0.2, 0.1, 0.7, 0.5, 0.3, 0.4]
        equilibrium_frequencies = [0.3, 0.4, 0.2, 0.1]
        model = msprime.GTRMutationModel(
            relative_rates=relative_rates,
            equilibrium_frequencies=equilibrium_frequencies,
        )
        self.verify_model(model, name="GTR")

    def test_PAM_stats(self):
        model = msprime.PAMMutationModel()
        self.verify_model(model, name="PAM")

    def test_BLOSUM62_stats(self):
        model = msprime.BLOSUM62MutationModel()
        self.verify_model(model, name="BLOSUM62")

    def test_arbitrary_model_stats(self):
        model = msprime.MatrixMutationModel(
            alleles=["abc", "", "x"],
            root_distribution=[0.8, 0.0, 0.2],
            transition_matrix=[[0.2, 0.4, 0.4], [0.1, 0.2, 0.7], [0.5, 0.3, 0.2]],
        )
        self.verify_model(
            model, name="arbitrary", state_independent=True, verify_rates=True
        )


class MutationRateMapTest(Test):
    def verify_subdivided(self, ts, rate_map, discrete=False):
        outfile = self._build_filename(None, "mutation_counts")
        sub_pos = np.unique(
            np.sort(
                np.concatenate(
                    [
                        rate_map.position,
                        np.linspace(
                            rate_map.start_position, rate_map.end_position, 101
                        ),
                    ]
                )
            )
        )
        sub_rate = [rate_map.get_rate(p) for p in sub_pos[:-1]]
        sub_map = msprime.RateMap(position=sub_pos, rate=sub_rate)
        t0 = msprime.sim_mutations(ts, rate=rate_map, discrete_genome=discrete).tables
        t1 = msprime.sim_mutations(ts, rate=sub_map, discrete_genome=discrete).tables
        if discrete:
            # make an equivalent map with breaks on integers
            # to use in calculating expected values
            int_pos = np.unique(np.ceil(rate_map.position))
            int_rate = [rate_map.get_rate(p) for p in int_pos[:-1]]
            rate_map = msprime.RateMap(int_pos, int_rate)

        bins = np.linspace(0, ts.sequence_length, min(101, int(ts.sequence_length + 1)))
        breaks = np.unique(np.sort(np.concatenate([bins, rate_map.position])))
        segsites = ts.segregating_sites(windows=breaks, mode="branch")
        expected = np.bincount(
            np.searchsorted(bins, breaks[:-1], "right") - 1,
            weights=segsites,
            minlength=len(bins) - 1,
        )
        for j in range(len(expected)):
            left = bins[j]
            right = bins[j + 1]
            mass = rate_map.get_cumulative_mass(right) - rate_map.get_cumulative_mass(
                left
            )
            expected[j] *= mass
        lower = scipy.stats.poisson.ppf(0.025, expected)
        upper = scipy.stats.poisson.ppf(1 - 0.025, expected)
        counts0 = np.bincount(
            np.digitize(t0.sites.position[t0.mutations.site], bins) - 1,
            minlength=len(bins) - 1,
        )
        counts1 = np.bincount(
            np.digitize(t1.sites.position[t1.mutations.site], bins) - 1,
            minlength=len(bins) - 1,
        )

        fig, ax = pyplot.subplots(1, 1, figsize=(8, 6))
        ax.scatter(bins[:-1], counts0, label="coarse map")
        ax.scatter(bins[:-1], counts1, label="fine map")
        ax.plot(bins[:-1], expected, label="expected number")
        ax.plot(bins[:-1], lower, "r:", linewidth=2, label="rough expected bounds")
        ax.plot(bins[:-1], upper, "r:", linewidth=2)
        ax.set_xlabel("Position")
        ax.set_ylabel("Num mutations in bin")
        ax.legend()
        pyplot.savefig(outfile, dpi=72)
        pyplot.close(fig)

    def test_subdivide(self):
        ts = msprime.sim_ancestry(
            1000,
            sequence_length=1e6,
            recombination_rate=1e-8,
            population_size=10000,
            random_seed=1,
        )
        rate_map = msprime.RateMap([0, 1e6], [1e-8])
        self.verify_subdivided(ts, rate_map)

    def test_varying_rate(self):
        ts = msprime.sim_ancestry(
            1000,
            sequence_length=1e6,
            recombination_rate=1e-8,
            population_size=10000,
            random_seed=1,
        )
        rate_map = msprime.RateMap([0, 3e5, 6e5, 1e6], [2e-8, 1e-9, 1e-8])
        self.verify_subdivided(ts, rate_map)

    def test_shorter_chromosome(self):
        ts = msprime.sim_ancestry(
            1000,
            sequence_length=20,
            recombination_rate=0.05,
            population_size=100,
            random_seed=12,
        )
        rate_map = msprime.RateMap(
            position=[0, 1.1, 10, 11.5, 13.8, 15.2, 15.9, 20],
            rate=[0.1, 0.2, 0.0, 0.2, 0.0, 100, 0.0],
        )
        self.verify_subdivided(ts, rate_map, discrete=True)


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
        for v in ts.variants(samples=range(ts.num_nodes), isolated_as_missing=False):
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
            "JC69": {"model_id": msprime.JC69MutationModel(), "par": ["-m", "HKY"]},
            "HKY": {
                "model_id": msprime.HKYMutationModel(
                    kappa=1.5, equilibrium_frequencies=[0.2, 0.3, 0.1, 0.4]
                ),
                "par": ["-m", "HKY", "-t", "0.75", "-f", "0.2,0.3,0.1,0.4"],
            },
            "F84": {
                "model_id": msprime.F84MutationModel(
                    kappa=1.0, equilibrium_frequencies=[0.3, 0.25, 0.2, 0.25]
                ),
                "par": ["-m", "F84", "-t", "0.5", "-f", "0.3,0.25,0.2,0.25"],
            },
            "GTR": {
                "model_id": msprime.GTRMutationModel(
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
            "PAM": {"model_id": msprime.PAMMutationModel(), "par": ["-m", "PAM"]},
            "BLOSUM62": {
                "model_id": msprime.BLOSUM62MutationModel(),
                "par": ["-m", "BLOSUM"],
            },
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

        for _ in range(num_replicates):
            ts = msprime.simulate(num_samples, Ne=Ne, length=length)
            ts_mutated = msprime.sim_mutations(
                ts,
                rate=mutation_rate,
                model=model_dict[model]["model_id"],
                discrete_genome=True,
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

    # Test methods

    def test_JC69(self):
        self._run_seq_gen_msprime_comparison("JC69")

    def test_HKY(self):
        self._run_seq_gen_msprime_comparison("HKY")

    def test_F84(self):
        self._run_seq_gen_msprime_comparison("F84")

    def test_GTR(self):
        self._run_seq_gen_msprime_comparison("GTR")

    def test_PAM(self):
        self._run_seq_gen_msprime_comparison("PAM")

    def test_BLOSUM62(self):
        self._run_seq_gen_msprime_comparison("BLOSUM62")


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
                "model_id": msprime.JC69MutationModel(),
                "pyvolve_model": pyvolve.Model("nucleotide"),
            },
            "HKY": {
                "model_id": msprime.HKYMutationModel(
                    kappa=1.5, equilibrium_frequencies=[0.2, 0.3, 0.1, 0.4]
                ),
                "pyvolve_model": pyvolve.Model(
                    "nucleotide", {"kappa": 1.5, "state_freqs": [0.2, 0.3, 0.1, 0.4]}
                ),
            },
            "PAM": {
                "model_id": msprime.PAMMutationModel(),
                "pyvolve_model": pyvolve.Model("DAYHOFFDCMUT"),
            },
            "BLOSUM62": {
                "model_id": msprime.BLOSUM62MutationModel(),
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

        for _ in range(num_replicates):
            ts = msprime.simulate(num_samples, Ne=1e4, length=length)
            ts_mutated = msprime.sim_mutations(
                ts,
                rate=mutation_rate,
                model=model_dict[model]["model_id"],
                discrete_genome=True,
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

    def test_pyv_JC69(self):
        self._run_pyvolve_comparison("JC69")

    def test_pyv_HKY(self):
        self._run_pyvolve_comparison("HKY")

    def test_pyv_PAM(self):
        self._run_pyvolve_comparison("PAM")

    def test_pyv_BLOSUM62(self):
        self._run_pyvolve_comparison("BLOSUM62")


class OlderMsprimeTest(Test):
    """
    Run tests against older versions of msprime.
    """

    def _run_in_venv(self, num_replicates, **kwargs):
        """
        Runs the specified simulation in the older version of msprime
        using a venv.
        """
        with tempfile.TemporaryDirectory(dir=self.output_dir.resolve()) as tempdir:
            tempdir = pathlib.Path(tempdir)
            params_file = tempdir / "params.pkl"
            output_prefix = tempdir / "output"
            with open(params_file, "wb") as f:
                pickle.dump(kwargs, f)
            cmd = (
                "cd data && ./msprime-0.7.4/bin/python run_old_msprime.py "
                f"{num_replicates} {params_file} {output_prefix}"
            )
            subprocess.run(cmd, shell=True, check=True)
            count = 0
            for trees_file in tempdir.glob("*.trees"):
                ts = tskit.load(trees_file)
                prov = json.loads(ts.provenance(0).record)
                assert prov["software"] == {"name": "msprime", "version": "0.7.4"}
                yield ts
                count += 1
            assert count == num_replicates

    def _run(self, num_replicates, **kwargs):
        logging.debug(f"Running: {num_replicates} replicates of {kwargs}")
        data = collections.defaultdict(list)
        old_version = "0.7.4"
        new_version = msprime.__version__
        iter1 = self._run_in_venv(num_replicates, **kwargs)
        iter2 = msprime.simulate(num_replicates=num_replicates, **kwargs)
        for ts1, ts2 in zip(iter1, iter2):
            assert ts1.sequence_length == ts2.sequence_length
            assert ts1.num_samples == ts2.num_samples
            for ts, version in [(ts1, old_version), (ts2, new_version)]:
                t_mrca = np.zeros(ts.num_trees)
                for tree in ts.trees():
                    t_mrca[tree.index] = tree.time(tree.root)
                data["tmrca_mean"].append(np.mean(t_mrca))
                data["num_trees"].append(ts.num_trees)
                data["num_nodes"].append(ts.num_nodes)
                data["num_edges"].append(ts.num_edges)
                data["version"].append(version)
        df = pd.DataFrame(data)

        df_new = df[df.version == new_version]
        df_old = df[df.version == old_version]
        for stat in ["tmrca_mean", "num_trees", "num_nodes", "num_edges"]:
            v1 = df_new[stat]
            v2 = df_old[stat]
            sm.graphics.qqplot(v1)
            sm.qqplot_2samples(v1, v2, line="45")
            pyplot.xlabel(new_version)
            pyplot.ylabel(old_version)
            f = self.output_dir / f"{stat}.png"
            pyplot.savefig(f, dpi=72)
            pyplot.close("all")

    def test_msprime_n1e2_no_recomb(self):
        self._run(10000, sample_size=100)

    def test_msprime_n1e4_no_recomb(self):
        self._run(1000, sample_size=10 ** 4)

    def test_msprime_n1e3_long_genome(self):
        self._run(
            1000, sample_size=10 ** 2, Ne=10 ** 4, recombination_rate=1e-8, length=1e6
        )

    def test_msprime_n1e2_long_genome(self):
        self._run(
            2000, sample_size=10 ** 2, Ne=10 ** 4, recombination_rate=1e-8, length=1e6
        )

    def test_msprime_n10_long_genome(self):
        self._run(1000, sample_size=10, Ne=10 ** 4, recombination_rate=1e-8, length=1e6)

    def test_msprime_n2_long_genome(self):
        self._run(1000, sample_size=2, Ne=10 ** 4, recombination_rate=1e-8, length=1e7)


###############################################
# Infrastructure for running the tests and CLI
###############################################


@attr.s
class TestInstance:
    """
    A single test instance, that consists of the test class and the test method
    name.
    """

    test_class = attr.ib()
    method_name = attr.ib()

    def run(self, basedir):
        logging.info(f"Running {self}")
        output_dir = pathlib.Path(basedir) / self.test_class / self.method_name
        output_dir.mkdir(parents=True, exist_ok=True)

        instance = getattr(sys.modules[__name__], self.test_class)(output_dir)
        method = getattr(instance, self.method_name)
        method()


@attr.s
class TestSuite:
    """
    Class responsible for registering all known tests.
    """

    tests = attr.ib(init=False, factory=dict)
    classes = attr.ib(init=False, factory=set)

    def register(self, test_class, method_name):
        test_instance = TestInstance(test_class, method_name)
        if method_name in self.tests:
            raise ValueError(f"Test name {method_name} already used.")
        self.tests[method_name] = test_instance
        self.classes.add(test_class)

    def get_tests(self, names=None, test_class=None):
        if names is not None:
            tests = [self.tests[name] for name in names]
        elif test_class is not None:
            tests = [
                test for test in self.tests.values() if test.test_class == test_class
            ]
        else:
            tests = list(self.tests.values())
        return tests


@attr.s
class TestRunner:
    """
    Class responsible for running test instances.
    """

    def __run_sequential(self, tests, basedir, progress):
        for test in tests:
            test.run(basedir)
            progress.update()

    def __run_parallel(self, tests, basedir, num_threads, progress):
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=num_threads
        ) as executor:
            futures = [executor.submit(test.run, basedir) for test in tests]
            exception = None
            for future in concurrent.futures.as_completed(futures):
                exception = future.exception()
                if exception is not None:
                    # At least tell the user that we've had an exception.
                    # Other stuff will still keep going, though.
                    logging.error("EXCEPTION:%s", exception)
                    break
                progress.update()
            if exception is not None:
                # Try to clear out as much work as we can, but it'll still run a
                # lot of stuff before we finish
                for future in futures:
                    future.cancel()
                raise exception

    def run(self, tests, basedir, num_threads, show_progress):
        progress = tqdm.tqdm(total=len(tests), disable=not show_progress)
        logging.info(f"running {len(tests)} tests using {num_threads} processes")
        if num_threads <= 1:
            self.__run_sequential(tests, basedir, progress)
        else:
            self.__run_parallel(tests, basedir, num_threads, progress)
        progress.close()


def setup_logging(args):
    log_level = "INFO"
    if args.quiet:
        log_level = "WARN"
    if args.debug:
        log_level = "DEBUG"

    daiquiri.setup(level=log_level)
    msprime_logger = daiquiri.getLogger("msprime")
    msprime_logger.setLevel("WARN")
    mpl_logger = daiquiri.getLogger("matplotlib")
    mpl_logger.setLevel("WARN")


def run_tests(suite, args):

    setup_logging(args)
    runner = TestRunner()

    if len(args.tests) > 0:
        tests = suite.get_tests(names=args.tests)
    elif args.test_class is not None:
        tests = suite.get_tests(test_class=args.test_class)
    else:
        tests = suite.get_tests()

    runner.run(tests, args.output_dir, args.num_threads, not args.no_progress)


def make_suite():
    suite = TestSuite()

    for cls_name, cls in inspect.getmembers(sys.modules[__name__]):
        if inspect.isclass(cls) and issubclass(cls, Test):
            test_class_instance = cls()
            for name, thing in inspect.getmembers(test_class_instance):
                if inspect.ismethod(thing):
                    if name.startswith("test_"):
                        suite.register(cls_name, name)
    return suite


def main():
    suite = make_suite()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test-class",
        "-c",
        default=None,
        choices=sorted(suite.classes),
        help="Run all tests for specified test class",
    )
    parser.add_argument(
        "tests",
        nargs="*",
        help="Run specific tests. Use the --list option to see those available",
    )
    parser.add_argument(
        "--output-dir",
        "-d",
        type=str,
        default="tmp__NOBACKUP__",
        help="specify the base output directory",
    )
    parser.add_argument(
        "--num-threads", "-t", type=int, default=1, help="Specify number of threads"
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--quiet", "-q", action="store_true", help="Do not write any output"
    )
    group.add_argument(
        "--debug", "-D", action="store_true", help="Write out debug output"
    )
    parser.add_argument(
        "--no-progress", "-n", action="store_true", help="Do not show progress bar"
    )
    parser.add_argument(
        "--list", "-l", action="store_true", help="List available checks and exit"
    )
    args = parser.parse_args()
    if args.list:
        print("All available tests")
        for test in suite.tests.values():
            print(test.test_class, test.method_name, sep="\t")
    else:
        run_tests(suite, args)


if __name__ == "__main__":
    main()
