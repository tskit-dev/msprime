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
Msprime is a reimplementation of Hudson's classical ms simulator for
modern datasets.
"""
from __future__ import print_function
from __future__ import division

import os
import random
import tempfile

import _msprime

__version__ = '1.0.0a1'

def simulate_trees(sample_size, num_loci, recombination_rate,
        population_models=[], max_memory="10M"):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns an iterator over the resulting trees.
    """
    fd, tf = tempfile.mkstemp(prefix="msp_", suffix=".dat")
    os.close(fd)
    try:
        sim = TreeSimulator(sample_size, tf)
        sim.set_num_loci(num_loci)
        sim.set_recombination_rate(recombination_rate)
        sim.set_max_memory(max_memory)
        for m in population_models:
            sim.add_population_model(m)
        sim.run()
        # We are done with the simulator so free up the memory
        del sim
        tree_file = TreeFile(tf)
        tree_file.sort()
        for l, pi, tau in tree_file:
            yield l, pi, tau
    finally:
        os.unlink(tf)

def simulate_tree(sample_size, population_models=[]):
    """
    Simulates the coalescent at a single locus for the specified sample size
    under the specified list of population models.
    """
    iterator = simulate_trees(sample_size, 1, 0, population_models)
    l, pi, tau = next(iterator)
    return pi, tau

class TreeSimulator(object):
    """
    Class to simulate trees under the standard neutral coalescent with
    recombination.
    """
    def __init__(self, sample_size, tree_file_name):
        self._sample_size = sample_size
        self._tree_file_name = tree_file_name
        self._recombination_rate = 0.5
        self._num_loci = 1
        self._population_models = []
        self._random_seed = None
        self._segment_block_size = None
        self._avl_node_block_size = None
        self._node_mapping_block_size = None
        self._max_memory = None
        self._ll_sim = None

    def get_num_trees(self):
        return self._ll_sim.get_num_trees()

    def get_time(self):
        return self._ll_sim.get_time()

    def get_num_coancestry_events(self):
        return self._ll_sim.get_num_coancestry_events()

    def get_num_recombination_events(self):
        return self._ll_sim.get_num_recombination_events()

    def add_population_model(self, pop_model):
        self._population_models.append(pop_model)

    def set_num_loci(self, num_loci):
        self._num_loci = num_loci

    def set_recombination_rate(self, recombination_rate):
        self._recombination_rate = recombination_rate

    def set_random_seed(self, random_seed):
        self._random_seed = random_seed

    def set_segment_block_size(self, segment_block_size):
        self._segment_block_size = segment_block_size

    def set_avl_node_block_size(self, avl_node_block_size):
        self._avl_node_block_size = avl_node_block_size

    def set_node_mapping_block_size(self, node_mapping_block_size):
        self._node_mapping_block_size = node_mapping_block_size

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
        if self._max_memory is None:
            self._max_memory = 10 * 1024 * 1024 # Very small by default.
        # TODO use the estimates to set these depending on the other parameter
        # values.
        if self._avl_node_block_size is None:
            self._avl_node_block_size = 1024
        if self._segment_block_size is None:
            self._segment_block_size = 1024
        if self._node_mapping_block_size is None:
            self._node_mapping_block_size = 1024
        if self._random_seed is None:
            self._random_seed = random.randint(0, 2**31 - 1)

    def run(self, max_time=None):
        """
        Runs the simulation until complete coalescence has occured or until
        the (optional) max_time.
        """
        assert self._ll_sim == None
        self._set_environment_defaults()
        self._ll_sim = _msprime.Simulator(sample_size=self._sample_size,
                num_loci=self._num_loci,
                recombination_rate=self._recombination_rate,
                random_seed=self._random_seed,
                tree_file_name=self._tree_file_name,
                max_memory=self._max_memory,
                segment_block_size=self._segment_block_size,
                avl_node_block_size=self._avl_node_block_size,
                node_mapping_block_size=self._node_mapping_block_size)
        if max_time == None:
            self._ll_sim.run()
        else:
            self._ll_sim.run(max_time)

class TreeFile(object):
    """
    Class to read trees from a tree file generated by the msprime simulation.
    """
    def __init__(self, tree_file_name):
        self._ll_tree_file = _msprime.TreeFile(tree_file_name)

    def sort(self):
        self._ll_tree_file.sort()

    def get_sample_size(self):
        return self._ll_tree_file.get_sample_size()

    def get_num_loci(self):
        return self._ll_tree_file.get_num_loci()

    def get_num_trees(self):
        return self._ll_tree_file.get_num_trees()

    def get_metadata(self):
        return self._ll_tree_file.get_metadata()

    def __iter__(self):
        n = 2 * self.get_sample_size()
        pi = [0 for j in range(n)]
        tau = [0 for j in range(n)]
        b = 1
        for l, r, c1, c2, p, t in self._ll_tree_file:
            if l != b:
                yield l - b, pi, tau
                b = l
            pi[c1] = p
            pi[c2] = p
            tau[p] = t
        yield self.get_num_loci() - l, pi, tau








