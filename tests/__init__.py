#
# Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
""""
Common code for the msprime test cases.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import collections
import random
import unittest

import _msprime


def setUp():
    # Make random tests reproducible.
    random.seed(210)


class PythonSparseTree(object):
    """
    Presents the same interface as the SparseTree object for testing. This
    is tightly coupled with the PythonTreeSequence object below which updates
    the internal structures during iteration.
    """
    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        self.parent = {}
        self.children = {}
        self.time = {}
        self.left = 0
        self.right = 0
        self.root = 0

    def get_num_nodes(self):
        return self.num_nodes

    def get_sample_size(self):
        return self.sample_size

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def get_parent(self, node):
        return self.parent[node]

    def get_children(self, node):
        return self.children[node]

    def get_time(self, node):
        return self.time[node]

    def get_root(self):
        return self.root

    def get_parent_dict(self):
        return self.parent

    def get_time_dict(self):
        return self.time


class PythonTreeSequence(object):
    """
    A python implementation of the TreeDiffIterator algorithm.
    """
    def __init__(self, tree_sequence, breakpoints=None):
        self._tree_sequence = tree_sequence
        self._sample_size = tree_sequence.get_sample_size()
        self._breakpoints = breakpoints

    def records(self):
        for j in range(self._tree_sequence.get_num_records()):
            yield self._tree_sequence.get_record(j, _msprime.MSP_ORDER_LEFT)

    def _diffs(self):
        left = 0
        used_records = collections.defaultdict(list)
        records_in = []
        for l, r, node, children, t in self.records():
            if l != left:
                yield l - left, used_records[left], records_in
                del used_records[left]
                records_in = []
                left = l
            used_records[r].append((node, children, t))
            records_in.append((node, children, t))
        yield r - left, used_records[left], records_in

    def _diffs_with_breaks(self):
        k = 1
        x = 0
        b = self._breakpoints
        for length, records_out, records_in in self._diffs():
            x += length
            yield b[k] - b[k - 1], records_out, records_in
            while self._breakpoints[k] != x:
                k += 1
                yield b[k] - b[k - 1], [], []
            k += 1

    def diffs(self, all_breaks=False):
        if all_breaks:
            return self._diffs_with_breaks()
        else:
            return self._diffs()

    def sparse_trees(self):
        st = PythonSparseTree(self._tree_sequence.get_num_nodes())
        st.sample_size = self._tree_sequence.get_sample_size()
        st.left = 0
        st.time = {j: 0 for j in range(1, st.sample_size + 1)}
        for length, records_out, records_in in self.diffs():
            for node, children, t in records_out:
                del st.time[node]
                del st.children[node]
                for c in children:
                    del st.parent[c]
            for node, children, t in records_in:
                st.time[node] = t
                st.children[node] = children
                for c in children:
                    st.parent[c] = node
            # Insert the root
            root = 1
            while root in st.parent:
                root = st.parent[root]
            st.parent[root] = 0
            st.root = root
            st.right += length
            assert len(st.parent) == 2 * st.sample_size - 1
            yield st
            del st.parent[root]
            st.left = st.right


class MsprimeTestCase(unittest.TestCase):
    """
    Superclass of all tests msprime simulator test cases.
    """
