# MIT License
#
# Copyright (c) 2018-2020 Tskit Developers
# Copyright (c) 2015-2018 University of Oxford
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
import msprime

SAMPLE_SIZES = [2]  # , 1000, 10_000, 50_000, 100_000, 200_000]
DEFAULTS = {
    "length": 1e6,
    "Ne": 1000,
    "recombination_rate": 1e-8,
    "mutation_rate": 1e-8,
    "random_seed": 42,
}


class SimulateLength:
    version = "1"
    param_names = ["sample_size", "length"]
    params = [SAMPLE_SIZES, [10, 1e4, 1e5, 5e5, 1e6]]

    def time_simulate_length(self, sample_size, length):
        msprime.simulate(sample_size, **{**DEFAULTS, "length": length})


class SimulateNe:
    version = "1"
    param_names = ["sample_size", "Ne"]
    params = [SAMPLE_SIZES, [1, 1e3, 2.5e3, 5e3, 1e4]]

    def time_simulate_Ne(self, sample_size, Ne):
        msprime.simulate(sample_size, **{**DEFAULTS, "Ne": Ne})


class SimulateRecombination:
    version = "1"
    param_names = ["sample_size", "recombination_rate"]
    params = [SAMPLE_SIZES, [0, 1e-9, 1e-8, 5e-8, 1e-7]]

    def time_simulate_recombination(self, sample_size, recombination_rate):
        msprime.simulate(
            sample_size, **{**DEFAULTS, "recombination_rate": recombination_rate}
        )


class SimulateMutation:
    version = "1"
    param_names = ["sample_size", "mutation_rate"]
    params = [SAMPLE_SIZES, [0, 1e-8, 1e-7, 1e-6, 1e-5]]

    def time_simulate_mutation(self, sample_size, mutation_rate):
        msprime.simulate(sample_size, **{**DEFAULTS, "mutation_rate": mutation_rate})


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
