#
# Copyright (C) 2014 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msprime is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msprime.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Module responsible to generating and reading tree files.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import json
import platform
import random
import sys

import _msprime
from _msprime import InputError
from _msprime import LibraryError
from . import __version__


def harmonic_number(n):
    """
    Returns the nth Harmonic number.
    """
    return sum(1 / k for k in range(1, n + 1))

def simulate_trees(sample_size, num_loci, scaled_recombination_rate,
        population_models=[], random_seed=None, max_memory="10M"):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns an iterator over the resulting trees.
    """
    sim = TreeSimulator(sample_size)
    sim.set_num_loci(num_loci)
    sim.set_scaled_recombination_rate(scaled_recombination_rate)
    if random_seed is not None:
        sim.set_random_seed(random_seed)
    sim.set_max_memory(max_memory)
    for m in population_models:
        sim.add_population_model(m)
    tree_sequence = sim.run()
    return tree_sequence.sparse_trees()

def simulate_tree(sample_size, population_models=[], random_seed=None,
        max_memory="10M"):
    """
    Simulates the coalescent at a single locus for the specified sample size
    under the specified list of population models.
    """
    iterator = simulate_trees(sample_size, 1, 0, population_models,
            random_seed, max_memory)
    l, pi, tau = next(iterator)
    return pi, tau

class TreeSimulator(object):
    """
    Class to simulate trees under the standard neutral coalescent with
    recombination.
    """
    def __init__(self, sample_size):
        self._sample_size = sample_size
        self._scaled_recombination_rate = 1.0
        self._num_loci = 1
        self._population_models = []
        self._random_seed = None
        self._segment_block_size = None
        self._avl_node_block_size = None
        self._node_mapping_block_size = None
        self._coalescence_record_block_size = None
        self._max_memory = None
        self._ll_sim = None

    def get_sample_size(self):
        return self._sample_size

    def get_scaled_recombination_rate(self):
        return self._scaled_recombination_rate

    def get_num_loci(self):
        return self._num_loci

    def get_random_seed(self):
        return self._random_seed

    def get_population_models(self):
        return self._population_models

    def get_num_breakpoints(self):
        return self._ll_sim.get_num_breakpoints()

    def get_breakpoints(self):
        return self._ll_sim.get_breakpoints()

    def get_used_memory(self):
        return self._ll_sim.get_used_memory()

    def get_time(self):
        return self._ll_sim.get_time()

    def get_num_avl_node_blocks(self):
        return self._ll_sim.get_num_avl_node_blocks()

    def get_num_coalescence_record_blocks(self):
        return self._ll_sim.get_num_coalescence_record_blocks()

    def get_num_node_mapping_blocks(self):
        return self._ll_sim.get_num_node_mapping_blocks()

    def get_num_segment_blocks(self):
        return self._ll_sim.get_num_segment_blocks()

    def get_num_coancestry_events(self):
        return self._ll_sim.get_num_coancestry_events()

    def get_num_recombination_events(self):
        return self._ll_sim.get_num_recombination_events()

    def get_max_memory(self):
        return self._ll_sim.get_max_memory()

    def add_population_model(self, pop_model):
        self._population_models.append(pop_model)

    def set_num_loci(self, num_loci):
        self._num_loci = num_loci

    def set_scaled_recombination_rate(self, scaled_recombination_rate):
        self._scaled_recombination_rate = scaled_recombination_rate

    def set_effective_population_size(self, effective_population_size):
        self._effective_population_size = effective_population_size

    def set_random_seed(self, random_seed):
        self._random_seed = random_seed

    def set_segment_block_size(self, segment_block_size):
        self._segment_block_size = segment_block_size

    def set_avl_node_block_size(self, avl_node_block_size):
        self._avl_node_block_size = avl_node_block_size

    def set_node_mapping_block_size(self, node_mapping_block_size):
        self._node_mapping_block_size = node_mapping_block_size

    def set_coalescence_record_block_size(self, coalescence_record_block_size):
        self._coalescence_record_block_size = coalescence_record_block_size

    def set_max_memory(self, max_memory):
        """
        Sets the approximate maximum memory used by the simulation
        to the specified value.  This can be suffixed with
        K, M or G to specify units of Kibibytes, Mibibytes or Gibibytes.
        """
        s = max_memory
        d = {"K":2**10, "M":2**20, "G":2**30}
        multiplier = 1
        value = s
        if s.endswith(tuple(d.keys())):
            value = s[:-1]
            multiplier = d[s[-1]]
        n = int(value)
        self._max_memory = n * multiplier

    def _set_environment_defaults(self):
        """
        Sets sensible default values for the memory usage parameters.
        """
        # Set the block sizes using our estimates.
        n = self._sample_size
        m = self._num_loci
        # First check to make sure they are sane.
        if not isinstance(n, int):
            raise TypeError("Sample size must be an integer")
        if not isinstance(m, int):
            raise TypeError("Number of loci must be an integer")
        if n < 2:
            raise ValueError("Sample size must be >= 2")
        if m < 1:
            raise ValueError("Postive number of loci required")
        rho = 4 * self._scaled_recombination_rate * (m - 1)
        num_trees = min(m // 2, rho * harmonic_number(n - 1))
        b = 10 # Baseline maximum
        num_trees = max(b, int(num_trees))
        num_avl_nodes = max(b, 4 * n + num_trees)
        # TODO This is probably much too large now.
        num_segments = max(b, int(0.0125 * n  * rho))
        if self._avl_node_block_size is None:
            self._avl_node_block_size = num_avl_nodes
        if self._segment_block_size is None:
            self._segment_block_size = num_segments
        if self._node_mapping_block_size is None:
            self._node_mapping_block_size = num_trees
        if self._coalescence_record_block_size is None:
            memory = 16 * 2**10  # 16M
            # Each coalescence record is 32bytes
            self._coalescence_record_block_size = memory // 32
        if self._random_seed is None:
            self._random_seed = random.randint(0, 2**31 - 1)
        if self._max_memory is None:
            self._max_memory = 10 * 1024 * 1024 # 10MiB by default

    def run(self):
        """
        Runs the simulation until complete coalescence has occured.
        """
        # Sort the models by start time
        models = sorted(self._population_models, key=lambda m: m.start_time)
        models = [m.get_ll_model() for m in models]
        assert self._ll_sim is None
        self._set_environment_defaults()
        self._ll_sim = _msprime.Simulator(
            sample_size=self._sample_size,
            num_loci=self._num_loci,
            population_models=models,
            scaled_recombination_rate=self._scaled_recombination_rate,
            random_seed=self._random_seed,
            max_memory=self._max_memory,
            segment_block_size=self._segment_block_size,
            avl_node_block_size=self._avl_node_block_size,
            node_mapping_block_size=self._node_mapping_block_size,
            coalescence_record_block_size=self._coalescence_record_block_size)
        self._ll_sim.run()
        ll_tree_sequence = _msprime.TreeSequence()
        ll_tree_sequence.create(self._ll_sim)
        ts = TreeSequence(ll_tree_sequence)
        return ts

    def reset(self):
        """
        Resets the simulation so that we can perform another replicate.
        """
        self._ll_sim = None

class TreeSequence(object):

    def __init__(self, ll_tree_sequence):
        self._ll_tree_sequence = ll_tree_sequence

    def get_ll_tree_sequence(self):
        return self._ll_tree_sequence

    def print_state(self):
        print("TODO")
        # print("parameters = ")
        # print(json.dumps(self._parameters, sort_keys=True, indent=4))
        # print("environment = ")
        # print(json.dumps(self._environment, sort_keys=True, indent=4))
        # for j in range(self._num_records):
        #     print(self._left[j], self._right[j], self._children[j],
        #             self._parent[j], self._time[j], sep="\t")

    def dump(self, path, zlib_compression=False):
        """
        Writes the tree sequence to the specified file path.
        """
        self._ll_tree_sequence.dump(path, zlib_compression)

    @classmethod
    def load(cls, path):
        ts = _msprime.TreeSequence()
        ts.load(path)
        return TreeSequence(ts)

    def get_sample_size(self):
        return self._ll_tree_sequence.get_sample_size()

    def get_num_loci(self):
        return self._ll_tree_sequence.get_num_loci()

    def get_num_records(self):
        return self._ll_tree_sequence.get_num_records()

    def get_num_mutations(self):
        return self._ll_tree_sequence.get_num_mutations()

    def get_breakpoints(self):
        return self._ll_tree_sequence.get_breakpoints()

    def get_mutations(self):
        return self._ll_tree_sequence.get_mutations()

    def records(self):
        for j in range(self.get_num_records()):
            yield self._ll_tree_sequence.get_record(j)

    def diffs(self, all_breaks=False):
        iterator = _msprime.TreeDiffIterator(self._ll_tree_sequence, all_breaks)
        return iterator

    def sparse_trees(self):
        n = self.get_sample_size()
        pi = {}
        tau = {j:0 for j in range(1, n + 1)}
        iterator = self.diffs()
        length, records_out, records_in = next(iterator)
        assert len(records_out) == 0
        for node, children, time in records_in:
            tau[node] = time
            pi[node] = 0
            for c in children:
                pi[c] = node
        yield length, pi, tau
        root = node
        del pi[root]

        for length, records_out, records_in in iterator:
            # print("ROOT = ", root)
            for node, children, time in records_out:
                # print("OUT:", children, node,time, sep="\t")
                del tau[node]
                for c in children:
                    del pi[c]
            for node, children, time in records_in:
                # print("IN :", children, node,time, sep="\t")
                tau[node] = time
                for c in children:
                    pi[c] = node
            # TODO this is a O(h) operation per tree, which seems
            # unneccessary. However, I can't seem to see a clean way to
            # keep track of the root using these records...
            v = 1
            while v in pi:
                v = pi[v]
            pi[v] = 0
            yield length, pi, tau
            del pi[v]

    def newick_trees(self, precision=3, breakpoints=None):
        iterator = _msprime.NewickConverter(self._ll_tree_sequence, precision)
        if breakpoints is None:
            for length, tree in iterator:
                yield length, tree
        else:
            trees_covered = 0
            j = 0
            for length, tree in iterator:
                trees_covered += length
                while breakpoints[j] < trees_covered:
                    j += 1
                    yield breakpoints[j] - breakpoints[j - 1], tree


    def generate_mutations(self, scaled_mutation_rate, random_seed=None):
        """
        Generates mutation according to the infinite sites model. This
        method over-writes any existing mutations stored in the tree
        sequence.
        """
        seed = random_seed
        if random_seed is None:
            seed = random.randint(0, 2**31)
        self._ll_tree_sequence.generate_mutations(scaled_mutation_rate, seed)


class HaplotypeGenerator(object):

    def __init__(self, tree_sequence):
        self._tree_sequence = tree_sequence
        ts = self._tree_sequence.get_ll_tree_sequence()
        self._ll_haplotype_generator = _msprime.HaplotypeGenerator(ts)

    def get_haplotype(self, sample_id):
        return self._ll_haplotype_generator.get_haplotype(sample_id)

    def haplotypes(self):
        for j in range(1, self._tree_sequence.get_sample_size() + 1):
            yield self.get_haplotype(j)


class PopulationModel(object):
    """
    Superclass of simulation population models.
    """
    def __init__(self, start_time):
        self.start_time = start_time

    def get_ll_model(self):
        """
        Returns the low-level model corresponding to this population
        model.
        """
        return self.__dict__

class ConstantPopulationModel(PopulationModel):
    """
    Class representing a constant-size population model. The size of this
    is expressed relative to the size of the population at sampling time.
    """
    def __init__(self, start_time, size):
        super(ConstantPopulationModel, self).__init__(start_time)
        self.size = size
        self.type = _msprime.POP_MODEL_CONSTANT


class ExponentialPopulationModel(PopulationModel):
    """
    Class representing an exponentially growing or shrinking population.
    TODO document model.
    """
    def __init__(self, start_time, alpha):
        super(ExponentialPopulationModel, self).__init__(start_time)
        self.alpha = alpha
        self.type = _msprime.POP_MODEL_EXPONENTIAL


