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

NULL_NODE = -1


def setUp():
    # Make random tests reproducible.
    random.seed(210)


class PythonSparseTree(object):
    """
    Presents the same interface as the SparseTree object for testing. This
    is tightly coupled with the PythonTreeSequence object below which updates
    the internal structures during iteration.
    """
    def __init__(self):
        self.parent = {}
        self.children = {}
        self.time = {}
        self.left = 0
        self.right = 0
        self.root = 0
        self.index = -1
        self.sample_size = 0
        # We need a mutations function, so this name is taken.
        self.mutation_list = []

    @classmethod
    def from_sparse_tree(self, sparse_tree):
        ret = PythonSparseTree()
        ret.root = sparse_tree.get_root()
        ret.sample_size = sparse_tree.get_sample_size()
        ret.left, ret.right = sparse_tree.get_interval()
        ret.mutation_list = list(sparse_tree.mutations())
        ret.index = sparse_tree.get_index()
        # Traverse the tree and update the details as we go
        # We don't use the traversal method here because this
        # is used to test them.
        stack = [sparse_tree.get_root()]
        while len(stack) > 0:
            u = stack.pop()
            ret.time[u] = sparse_tree.get_time(u)
            if sparse_tree.is_internal(u):
                c = sparse_tree.get_children(u)
                stack.extend(c)
                for child in c:
                    ret.parent[child] = u
                ret.children[u] = c
        ret.parent[sparse_tree.get_root()] = NULL_NODE
        assert ret == sparse_tree
        return ret

    def _preorder_nodes(self, u, l):
        l.append(u)
        if u in self.children:
            for c in self.children[u]:
                self._preorder_nodes(c, l)

    def nodes(self, root=None, order="preorder"):
        u = root
        if root is None:
            u = self.root
        if order == "preorder":
            l = []
            self._preorder_nodes(u, l)
            return iter(l)
        else:
            raise ValueError("order not supported")

    def get_sample_size(self):
        return self.sample_size

    def get_interval(self):
        return self.left, self.right

    def get_parent(self, node):
        return self.parent[node]

    def get_children(self, node):
        return self.children[node]

    def get_time(self, node):
        return self.time[node]

    def get_root(self):
        return self.root

    def get_index(self):
        return self.index

    def get_parent_dict(self):
        return self.parent

    def get_time_dict(self):
        return self.time

    def mutations(self):
        return iter(self.mutation_list)


class PythonTreeSequence(object):
    """
    A python implementation of the TreeDiffIterator algorithm.
    """
    def __init__(self, tree_sequence, breakpoints=None):
        self._tree_sequence = tree_sequence
        self._sample_size = tree_sequence.get_sample_size()
        self._breakpoints = breakpoints
        self._mutations = tree_sequence.get_mutations()

    def records(self):
        for j in range(self._tree_sequence.get_num_records()):
            yield self._tree_sequence.get_record(j, _msprime.MSP_ORDER_LEFT)

    def _diffs(self):
        left = 0
        used_records = collections.defaultdict(list)
        records_in = []
        for l, r, node, children, t, _ in self.records():
            if l != left:
                # Out records must be sorted in reverse time order.
                records_out = sorted(used_records[left], key=lambda x: -x[2])
                yield l - left, records_out, records_in
                del used_records[left]
                records_in = []
                left = l
            used_records[r].append((node, children, t))
            records_in.append((node, children, t))
        records_out = sorted(used_records[left], key=lambda x: -x[2])
        yield r - left, records_out, records_in

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

    def trees(self):
        st = PythonSparseTree()
        st.sample_size = self._tree_sequence.get_sample_size()
        st.left = 0
        st.time = {j: 0 for j in range(st.sample_size)}
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
            root = 0
            while root in st.parent:
                root = st.parent[root]
            st.parent[root] = NULL_NODE
            st.root = root
            st.index += 1
            st.right += length
            # Add in all the mutations
            st.mutation_list = [
                (p, u, j) for (p, u, j) in self._mutations
                if st.left <= p < st.right
            ]
            yield st
            del st.parent[root]
            st.left = st.right


class PythonRecombinationMap(object):
    """
    A Python implementation of the RecombinationMap interface.
    """
    def __init__(self, positions, rates, num_loci):
        assert len(positions) == len(rates)
        assert len(positions) >= 2
        assert sorted(positions) == positions
        assert positions[0] == 0
        assert positions[-1] == 1
        self._positions = positions
        self._rates = rates
        self._num_loci = num_loci

    def get_total_recombination_rate(self):
        """
        Returns the effective recombination rate for this genetic map.
        This is the weighted mean of the rates across all intervals.
        """
        x = self._positions
        effective_rate = 0
        for j in range(len(x) - 1):
            length = (x[j + 1] - x[j])
            effective_rate += self._rates[j] * length
        return effective_rate

    def physical_to_genetic(self, x):
        if self.get_total_recombination_rate() == 0:
            ret = x
        else:
            s = 0
            last_phys_x = 0
            j = 1
            while j < len(self._positions) - 1 and x > self._positions[j]:
                phys_x = self._positions[j]
                rate = self._rates[j - 1]
                s += (phys_x - last_phys_x) * rate
                j += 1
                last_phys_x = phys_x
            rate = self._rates[j - 1]
            s += (x - last_phys_x) * rate
            ret = 0
            if self.get_total_recombination_rate() > 0:
                ret = s / self.get_total_recombination_rate()
        return ret * self._num_loci

    def genetic_to_physical(self, v):
        if self.get_total_recombination_rate() == 0:
            return v / self._num_loci
        # v is expressed in [0, m]. Rescale it back into the range
        # (0, total_mass).
        u = (v / self._num_loci) * self.get_total_recombination_rate()
        s = 0
        last_phys_x = 0
        rate = self._rates[0]
        j = 1
        while j < len(self._positions) and s < u:
            phys_x = self._positions[j]
            rate = self._rates[j - 1]
            s += (phys_x - last_phys_x) * rate
            j += 1
            last_phys_x = phys_x
        y = last_phys_x - (s - u) / rate
        return y


class MRCACalculator(object):
    """
    Class to that allows us to compute the nearest common ancestor of arbitrary
    nodes in an oriented forest.

    This is an implementation of Schieber and Vishkin's nearest common ancestor
    algorithm from TAOCP volume 4A, pg.164-167 [K11]_. Preprocesses the
    input tree into a sideways heap in O(n) time and processes queries for the
    nearest common ancestor between an arbitary pair of nodes in O(1) time.

    :param oriented_forest: the input oriented forest
    :type oriented_forest: list of integers
    """
    LAMBDA = 0

    def __init__(self, oriented_forest):
        # We turn this oriened forest into a 1 based array by adding 1
        # to everything
        converted = [0] + [x + 1 for x in oriented_forest]
        self.__preprocess(converted)

    def __preprocess(self, oriented_forest):
        """
        Preprocess the oriented forest, so that we can answer mrca queries
        in constant time.
        """
        n = len(oriented_forest)
        child = [self.LAMBDA for i in range(n)]
        parent = [self.LAMBDA for i in range(n)]
        sib = [self.LAMBDA for i in range(n)]
        self.__lambda = [0 for i in range(n)]
        self.__pi = [0 for i in range(n)]
        self.__tau = [0 for i in range(n)]
        self.__beta = [0 for i in range(n)]
        self.__alpha = [0 for i in range(n)]
        for u in range(n):
            v = oriented_forest[u]
            sib[u] = child[v]
            child[v] = u
            parent[u] = v
        p = child[self.LAMBDA]
        n = 0
        self.__lambda[0] = -1
        while p != self.LAMBDA:
            notDone = True
            while notDone:
                n += 1
                self.__pi[p] = n
                self.__tau[n] = self.LAMBDA
                self.__lambda[n] = 1 + self.__lambda[n >> 1]
                if child[p] != self.LAMBDA:
                    p = child[p]
                else:
                    notDone = False
            self.__beta[p] = n
            notDone = True
            while notDone:
                self.__tau[self.__beta[p]] = parent[p]
                if sib[p] != self.LAMBDA:
                    p = sib[p]
                    notDone = False
                else:
                    p = parent[p]
                    if p != self.LAMBDA:
                        h = self.__lambda[n & -self.__pi[p]]
                        self.__beta[p] = ((n >> h) | 1) << h
                    else:
                        notDone = False
        # Begin the second traversal
        self.__lambda[0] = self.__lambda[n]
        self.__pi[self.LAMBDA] = 0
        self.__beta[self.LAMBDA] = 0
        self.__alpha[self.LAMBDA] = 0
        p = child[self.LAMBDA]
        while p != self.LAMBDA:
            notDone = True
            while notDone:
                a = (
                    self.__alpha[parent[p]] |
                    (self.__beta[p] & -self.__beta[p])
                )
                self.__alpha[p] = a
                if child[p] != self.LAMBDA:
                    p = child[p]
                else:
                    notDone = False
            notDone = True
            while notDone:
                if sib[p] != self.LAMBDA:
                    p = sib[p]
                    notDone = False
                else:
                    p = parent[p]
                    notDone = p != self.LAMBDA

    def get_mrca(self, x, y):
        """
        Returns the most recent common ancestor of the nodes x and y,
        or -1 if the nodes belong to different trees.

        :param x: the first node
        :param y: the second node
        :return: the MRCA of nodes x and y
        """
        # WE need to rescale here because SV expects 1-based arrays.
        return self._sv_mrca(x + 1, y + 1) - 1

    def _sv_mrca(self, x, y):
        if self.__beta[x] <= self.__beta[y]:
            h = self.__lambda[self.__beta[y] & -self.__beta[x]]
        else:
            h = self.__lambda[self.__beta[x] & -self.__beta[y]]
        k = self.__alpha[x] & self.__alpha[y] & -(1 << h)
        h = self.__lambda[k & -k]
        j = ((self.__beta[x] >> h) | 1) << h
        if j == self.__beta[x]:
            xhat = x
        else:
            l = self.__lambda[self.__alpha[x] & ((1 << h) - 1)]
            xhat = self.__tau[((self.__beta[x] >> l) | 1) << l]
        if j == self.__beta[y]:
            yhat = y
        else:
            l = self.__lambda[self.__alpha[y] & ((1 << h) - 1)]
            yhat = self.__tau[((self.__beta[y] >> l) | 1) << l]
        if self.__pi[xhat] <= self.__pi[yhat]:
            z = xhat
        else:
            z = yhat
        return z


class MsprimeTestCase(unittest.TestCase):
    """
    Superclass of all tests msprime simulator test cases.
    """
