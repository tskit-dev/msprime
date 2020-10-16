# MIT License
#
# Copyright (c) 2020 Tskit Developers
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Benchmarks for msprime using airspeed velocity. Please see the developer
documentation for details on how to run these and how to develop your
own benchmarks.
"""
try:
    import stdpopsim

    stdpopsim_available = True
except TypeError:
    stdpopsim_available = False

import msprime


class LargeSimulationBenchmark:

    # ASV is designed to accurately time functions that execute in a fraction
    # of a second. But, we're interested in profiling large simulations that
    # run in 10s of seconds (at least). We want to run the target function
    # *exactly* once, and we need to set all these variables to do that.
    warmup_time = 0
    processes = 1
    repeat = 1
    number = 1
    rounds = 1
    min_run_count = 1
    timeout = 120

    def setup(self):
        # Stuff that depends on the recomb_map_chr22 will fail
        if stdpopsim_available:
            species = stdpopsim.get_species("HomSap")
            genetic_map = species.get_genetic_map("HapMapII_GRCh37")
            self.recomb_map_chr22 = genetic_map.get_chromosome_map("chr22")


class Hudson(LargeSimulationBenchmark):
    def _run_large_sample_size(self):
        msprime.simulate(
            sample_size=10 ** 6,
            length=1e7,
            Ne=10 ** 4,
            recombination_rate=1e-8,
            random_seed=42,
        )

    def time_large_sample_size(self):
        self._run_large_sample_size()

    def peakmem_large_sample_size(self):
        self._run_large_sample_size()

    def _run_long_sequence_length(self):
        msprime.simulate(
            sample_size=100,
            length=1e8,
            Ne=10 ** 4,
            recombination_rate=1e-8,
            random_seed=42,
        )

    def time_long_sequence_length(self):
        self._run_long_sequence_length()

    def peakmem_long_sequence_length(self):
        self._run_long_sequence_length()

    def _run_long_sequence_length_gene_conversion(self):
        msprime.sim_ancestry(
            sample_size=100,
            length=1e8,
            Ne=10 ** 4,
            gene_conversion_rate=1e-8,
            # 100Kb tract length.
            gene_conversion_tract_length=100 * 1e3,
            random_seed=43,
        )

    def time_long_sequence_length_gene_conversion(self):
        self._run_long_sequence_length()

    def peakmem_long_sequence_length_gene_conversion(self):
        self._run_long_sequence_length()

    def _run_human_chr22(self):
        msprime.simulate(
            sample_size=100,
            Ne=10 ** 4,
            recombination_map=self.recomb_map_chr22,
            random_seed=234,
        )

    def time_human_chr22(self):
        self._run_human_chr22()

    def peakmem_human_chr22(self):
        self._run_human_chr22()

    def _run_many_replicates(self):
        for _ in msprime.simulate(10, num_replicates=10 ** 5, random_seed=1234):
            pass

    def time_many_replicates(self):
        self._run_many_replicates()

    def peakmem_many_replicates(self):
        self._run_many_replicates()

    # 2 populations, high migration.
    # Lots of populations, 1D stepping stone.


class DTWF(LargeSimulationBenchmark):
    def _run_large_population_size(self):
        msprime.simulate(
            sample_size=1000,
            Ne=10 ** 6,
            length=1e5,
            recombination_rate=1e-8,
            random_seed=42,
            model="dtwf",
            end_time=1000,
        )

    def time_large_population_size(self):
        self._run_large_population_size()

    def peakmem_large_population_size(self):
        self._run_large_population_size()

    def _run_long_sequence_length(self):
        msprime.simulate(
            sample_size=100,
            Ne=10 ** 4,
            length=1e7,
            recombination_rate=1e-8,
            random_seed=42,
            model="dtwf",
            # Tuning this to give ~30s runtime.
            end_time=5e4,
        )

    def time_long_sequence_length(self):
        self._run_long_sequence_length()

    def peakmem_long_sequence_length(self):
        self._run_long_sequence_length()

    def _run_human_chr22(self):
        msprime.simulate(
            sample_size=100,
            Ne=10 ** 4,
            recombination_map=self.recomb_map_chr22,
            random_seed=234,
            end_time=10000,
            model="dtwf",
        )

    def time_human_chr22(self):
        self._run_human_chr22()

    def peakmem_human_chr22(self):
        self._run_human_chr22()

    def _run_many_replicates(self):
        reps = msprime.simulate(
            10,
            Ne=100,
            num_replicates=10 ** 5,
            random_seed=1234,
            model="dtwf",
            end_time=100,
        )
        for _ in reps:
            pass

    def time_many_replicates(self):
        self._run_many_replicates()

    def peakmem_many_replicates(self):
        self._run_many_replicates()
