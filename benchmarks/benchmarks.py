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
Benchmarks for msprime - use asv to run
"""
import stdpopsim

import msprime

DEFAULTS = {
    "sample_size": 1000,
    "length": 1e6,
    "Ne": 10 ** 4,
    "recombination_rate": 1e-8,
    "random_seed": 42,
}


class HudsonLength:
    version = "1"
    param_names = ["length"]
    params = [10 ** j for j in range(2, 8)]

    def time_hudson_length(self, length):
        msprime.simulate(**{**DEFAULTS, "length": length})


class HudsonSampleSize:
    version = "1"
    param_names = ["sample_size"]
    params = [10 ** j for j in range(1, 7)]

    def time_hudson_sample_size(self, sample_size):
        msprime.simulate(**{**DEFAULTS, "sample_size": sample_size})


class HumanRecombinationMapSuite:
    version = "1"

    def setup(self):
        species = stdpopsim.get_species("HomSap")
        genetic_map = species.get_genetic_map("HapMapII_GRCh37")
        self.recomb_map_chr22 = genetic_map.get_chromosome_map("chr22")

    def time_hudson_chr22(self):
        msprime.simulate(
            100,
            Ne=10 ** 4,
            recombination_map=self.recomb_map_chr22,
            random_seed=234,
            end_time=1000,
        )

    def time_dtwf_chr22(self):
        msprime.simulate(
            100,
            Ne=10 ** 4,
            recombination_map=self.recomb_map_chr22,
            random_seed=234,
            model="dtwf",
            end_time=1000,
        )


# class SimulateNe:
#     version = "1"
#     param_names = ["sample_size", "Ne"]
#     params = [SAMPLE_SIZES, [1, 1e3, 2.5e3, 5e3, 1e4]]

#     def time_simulate_Ne(self, sample_size, Ne):
#         msprime.simulate(sample_size, **{**DEFAULTS, "Ne": Ne})


# class SimulateRecombination:
#     version = "1"
#     param_names = ["sample_size", "recombination_rate"]
#     params = [SAMPLE_SIZES, [0, 1e-9, 1e-8, 5e-8, 1e-7]]

#     def time_simulate_recombination(self, sample_size, recombination_rate):
#         msprime.simulate(
#             sample_size, **{**DEFAULTS, "recombination_rate": recombination_rate}
#         )


# class SimulateMutation:
#     version = "1"
#     param_names = ["sample_size", "mutation_rate"]
#     params = [SAMPLE_SIZES, [0, 1e-8, 1e-7, 1e-6, 1e-5]]

#     def time_simulate_mutation(self, sample_size, mutation_rate):
#         msprime.simulate(sample_size, **{**DEFAULTS, "mutation_rate": mutation_rate})


# class MemSuite:
#     param_names = ['sample_size', 'length']
#     params = [
#                 [2, 1e3, 1e6],
#                 [10, 1e4, 1e6]
#     ]
#
#     def peakmem_simulate(self, sample_size, length):
#         msprime.simulate(sample_size, length=length, Ne=1000, recombination_rate=2e-8,
#                          random_seed=42)
