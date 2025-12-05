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


class HudsonlargeSampleSize(LargeSimulationBenchmark):

    def _get_params(self):
        return {
            "samples": 0.5 * (10**6),
            "sequence_length": 1e7,
            "population_size": 10**4,
            "recombination_rate": 1e-8,
            "random_seed": 42,
        }

    def run(self):
        return msprime.sim_ancestry(**self._get_params())

    def time_test(self):
        self.run()

    def peakmem_test(self):
        self.run()


class HudsonlargeSampleSizeOverRoot(HudsonlargeSampleSize):
    def _get_params(self):
        return {
            **super()._get_params(),
            "stop_at_local_mrca": False,
        }


class HudsonLongSequenceLength(HudsonlargeSampleSize):
    def _get_params(self):
        return {
            **super()._get_params(),
            "sequence_length": 1e8,
            "samples": 50,
        }


class HudsonLongSequenceLengthGeneConversion(HudsonlargeSampleSize):
    def _get_params(self):
        return {
            "sequence_length": 1e8,
            "samples": 50,
            "gene_conversion_rate": 1e-8,
            "gene_conversion_tract_length": 100 * 1e3,
            "random_seed": 43,
        }


class HudsonHumanChr22(HudsonlargeSampleSize):

    def _get_params(self):
        return {
            **super()._get_params(),
            "sequence_length": None,
            "samples": 50,
            "recombination_rate": self.recomb_map_chr22,
        }


class HudsonManyReplicates(HudsonlargeSampleSize):

    def run(self):
        params = {"samples": 10, "num_replicates": 10**5, "random_seed": 1234}
        for _ in msprime.sim_ancestry(**params):
            pass


class HudsonHumanChr22OverRoot(HudsonlargeSampleSize):
    def _get_params(self):
        return {
            **super()._get_params(),
            "sequence_length": None,
            "samples": 50,
            "recombination_rate": self.recomb_map_chr22,
            "stop_at_local_mrca": False,
        }


class DTWFLargePopulationSize(LargeSimulationBenchmark):
    def _get_params(self):
        return {
            "samples": 500,
            "sequence_length": 1e5,
            "population_size": 10**6,
            "recombination_rate": 1e-8,
            "random_seed": 42,
            "model": "dtwf",
            "end_time": 1000,
        }

    def run(self):
        return msprime.sim_ancestry(**self._get_params())

    def time_test(self):
        self.run()

    def peakmem_test(self):
        self.run()


class DTWFLongSequenceLength(DTWFLargePopulationSize):
    def _get_params(self):
        return {
            **super()._get_params(),
            "sequence_length": 1e7,
            "samples": 50,
            "end_time": 5e4,
            "population_size": 10**4,
        }


class DTWFHumanChr22(DTWFLargePopulationSize):

    def _get_params(self):
        return {
            **super()._get_params(),
            "sequence_length": None,
            "samples": 50,
            "recombination_rate": self.recomb_map_chr22,
            "end_time": 10000,
            "population_size": 10**4,
        }


class DTWFManyReplicates(DTWFLargePopulationSize):
    def run(self):
        params = {
            "samples": 5,
            "population_size": 100,
            "num_replicates": 10**5,
            "random_seed": 1234,
            "model": "dtwf",
            "end_time": 100,
        }
        for _ in msprime.sim_ancestry(**params):
            pass
