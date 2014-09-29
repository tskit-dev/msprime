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


class TreeSimulator(object):
    """
    Class to simulate trees under the standard neutral coalescent with
    recombination.
    """
    def __init__(self, sample_size):
        self._sample_size = sample_size
        self._recombination_rate = 0.5
        self._num_loci = 1
        self._population_models = []
        self._random_seed = None
        self._tree_file_name = None
        self._delete_tree_file = False
        self._segment_block_size = None
        self._avl_node_block_size = None
        self._node_mapping_block_size = None
        self._max_memory = None
        self._ll_sim = None

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
        if self._tree_file_name is None:
            fd, self._tree_file_name = tempfile.mkstemp(prefix="msp_",
                    suffix=".dat")
            os.close(fd)
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

    def get_trees(self):
        """
        Returns an iterator over the trees in tree file.
        """

