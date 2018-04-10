#
# Copyright (C) 2015-2018 University of Oxford
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

import random
import sys
import unittest
import base64

try:
    # We run some tests on the CLI to make sure that we can work in a minimal
    # sense without numpy. We should extract the PythonSimplifier (and other
    # algorithms) out into their own module so we don't pull in numpy for
    # all tests.
    import numpy as np
except ImportError:
    pass

import msprime

NULL_NODE = -1


def setUp():
    # Make random tests reproducible.
    random.seed(210)


class MsprimeTestCase(unittest.TestCase):
    """
    Superclass of all tests msprime simulator test cases.
    """


class PythonSparseTree(object):
    """
    Presents the same interface as the SparseTree object for testing. This
    is tightly coupled with the PythonTreeSequence object below which updates
    the internal structures during iteration.
    """
    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        self.parent = [msprime.NULL_NODE for _ in range(num_nodes)]
        self.left_child = [msprime.NULL_NODE for _ in range(num_nodes)]
        self.right_child = [msprime.NULL_NODE for _ in range(num_nodes)]
        self.left_sib = [msprime.NULL_NODE for _ in range(num_nodes)]
        self.right_sib = [msprime.NULL_NODE for _ in range(num_nodes)]
        self.above_sample = [False for _ in range(num_nodes)]
        self.is_sample = [False for _ in range(num_nodes)]
        self.left = 0
        self.right = 0
        self.root = 0
        self.index = -1
        self.left_root = -1
        # We need a sites function, so this name is taken.
        self.site_list = []

    @classmethod
    def from_sparse_tree(cls, sparse_tree):
        ret = PythonSparseTree(sparse_tree.num_nodes)
        ret.left, ret.right = sparse_tree.get_interval()
        ret.site_list = list(sparse_tree.sites())
        ret.index = sparse_tree.get_index()
        ret.left_root = sparse_tree.left_root
        for u in range(ret.num_nodes):
            ret.parent[u] = sparse_tree.parent(u)
            ret.left_child[u] = sparse_tree.left_child(u)
            ret.right_child[u] = sparse_tree.right_child(u)
            ret.left_sib[u] = sparse_tree.left_sib(u)
            ret.right_sib[u] = sparse_tree.right_sib(u)
        assert ret == sparse_tree
        return ret

    @property
    def roots(self):
        u = self.left_root
        roots = []
        while u != msprime.NULL_NODE:
            roots.append(u)
            u = self.right_sib[u]
        return roots

    def children(self, u):
        v = self.left_child[u]
        ret = []
        while v != msprime.NULL_NODE:
            ret.append(v)
            v = self.right_sib[v]
        return ret

    def _preorder_nodes(self, u, l):
        l.append(u)
        for c in self.children(u):
            self._preorder_nodes(c, l)

    def _postorder_nodes(self, u, l):
        for c in self.children(u):
            self._postorder_nodes(c, l)
        l.append(u)

    def _inorder_nodes(self, u, l):
        children = self.children(u)
        if len(children) > 0:
            mid = len(children) // 2
            for v in children[:mid]:
                self._inorder_nodes(v, l)
            l.append(u)
            for v in children[mid:]:
                self._inorder_nodes(v, l)
        else:
            l.append(u)

    def _levelorder_nodes(self, u, l, level):
        l[level].append(u) if level < len(l) else l.append([u])
        for c in self.children(u):
            self._levelorder_nodes(c, l, level + 1)

    def nodes(self, root=None, order="preorder"):
        roots = [root]
        if root is None:
            roots = self.roots
        for u in roots:
            node_list = []
            if order == "preorder":
                self._preorder_nodes(u, node_list)
            elif order == "inorder":
                self._inorder_nodes(u, node_list)
            elif order == "postorder":
                self._postorder_nodes(u, node_list)
            elif order == "levelorder" or order == "breadthfirst":
                # Returns nodes in their respective levels
                # Nested list comprehension flattens node_list in order
                self._levelorder_nodes(u, node_list, 0)
                node_list = iter([i for level in node_list for i in level])
            else:
                raise ValueError("order not supported")
            for v in node_list:
                yield v

    def get_interval(self):
        return self.left, self.right

    def get_parent(self, node):
        return self.parent[node]

    def get_children(self, node):
        return self.children[node]

    def get_index(self):
        return self.index

    def get_parent_dict(self):
        d = {
            u: self.parent[u] for u in range(self.num_nodes)
            if self.parent[u] != msprime.NULL_NODE}
        return d

    def sites(self):
        return iter(self.site_list)

    def __eq__(self, other):
        return (
            self.get_parent_dict() == other.get_parent_dict() and
            self.get_interval() == other.get_interval() and
            self.roots == other.roots and
            self.get_index() == other.get_index() and
            list(self.sites()) == list(other.sites()))

    def __ne__(self, other):
        return not self.__eq__(other)

    def newick(self, precision=0, time_scale=0):
        # We only support 0 branch lengths here because this information isn't
        # immediately available.
        assert time_scale == 0 and precision == 0
        assert len(self.roots) == 1
        return self._build_newick(self.left_root) + ";"

    def _build_newick(self, node):
        if self.left_child[node] == msprime.NULL_NODE:
            s = "{0}".format(node + 1)
        else:
            s = "("
            for child in self.children(node):
                s += self._build_newick(child) + ":0,"
            s = s[:-1] + ")"
        return s


class PythonTreeSequence(object):
    """
    A python implementation of the TreeSequence object.
    """
    def __init__(self, tree_sequence, breakpoints=None):
        self._tree_sequence = tree_sequence
        self._num_samples = tree_sequence.get_num_samples()
        self._breakpoints = breakpoints
        self._sites = []

        def make_mutation(id_):
            site, node, derived_state, parent, metadata = tree_sequence.get_mutation(id_)
            return msprime.Mutation(
                id_=id_, site=site, node=node, derived_state=derived_state,
                parent=parent, metadata=metadata)
        for j in range(tree_sequence.get_num_sites()):
            pos, ancestral_state, ll_mutations, id_, metadata = tree_sequence.get_site(j)
            self._sites.append(msprime.Site(
                id_=id_, position=pos, ancestral_state=ancestral_state,
                mutations=[make_mutation(ll_mut) for ll_mut in ll_mutations],
                metadata=metadata))

    def edge_diffs(self):
        M = self._tree_sequence.get_num_edges()
        sequence_length = self._tree_sequence.get_sequence_length()
        edges = [msprime.Edge(*self._tree_sequence.get_edge(j)) for j in range(M)]
        time = [self._tree_sequence.get_node(edge.parent)[1] for edge in edges]
        in_order = sorted(range(M), key=lambda j: (
            edges[j].left, time[j], edges[j].parent, edges[j].child))
        out_order = sorted(range(M), key=lambda j: (
            edges[j].right, -time[j], -edges[j].parent, -edges[j].child))
        j = 0
        k = 0
        left = 0
        while j < M or left < sequence_length:
            e_out = []
            e_in = []
            while k < M and edges[out_order[k]].right == left:
                h = out_order[k]
                e_out.append(edges[h])
                k += 1
            while j < M and edges[in_order[j]].left == left:
                h = in_order[j]
                e_in.append(edges[h])
                j += 1
            right = sequence_length
            if j < M:
                right = min(right, edges[in_order[j]].left)
            if k < M:
                right = min(right, edges[out_order[k]].right)
            yield (left, right), e_out, e_in
            left = right

    def trees(self):
        M = self._tree_sequence.get_num_edges()
        sequence_length = self._tree_sequence.get_sequence_length()
        edges = [
            msprime.Edge(*self._tree_sequence.get_edge(j)) for j in range(M)]
        t = [
            self._tree_sequence.get_node(j)[1]
            for j in range(self._tree_sequence.get_num_nodes())]
        in_order = sorted(
            range(M), key=lambda j: (
                edges[j].left, t[edges[j].parent], edges[j].parent, edges[j].child))
        out_order = sorted(
            range(M), key=lambda j: (
                edges[j].right, -t[edges[j].parent], -edges[j].parent, -edges[j].child))
        j = 0
        k = 0
        N = self._tree_sequence.get_num_nodes()
        st = PythonSparseTree(N)

        samples = list(self._tree_sequence.get_samples())
        for l in range(len(samples)):
            if l < len(samples) - 1:
                st.right_sib[samples[l]] = samples[l + 1]
            if l > 0:
                st.left_sib[samples[l]] = samples[l - 1]
            st.above_sample[samples[l]] = True
            st.is_sample[samples[l]] = True

        st.left_root = msprime.NULL_NODE
        if len(samples) > 0:
            st.left_root = samples[0]

        u = st.left_root
        roots = []
        while u != -1:
            roots.append(u)
            v = st.right_sib[u]
            if v != -1:
                assert st.left_sib[v] == u
            u = v

        st.left = 0
        while j < M or st.left < sequence_length:
            while k < M and edges[out_order[k]].right == st.left:
                p = edges[out_order[k]].parent
                c = edges[out_order[k]].child
                k += 1

                lsib = st.left_sib[c]
                rsib = st.right_sib[c]
                if lsib == msprime.NULL_NODE:
                    st.left_child[p] = rsib
                else:
                    st.right_sib[lsib] = rsib
                if rsib == msprime.NULL_NODE:
                    st.right_child[p] = lsib
                else:
                    st.left_sib[rsib] = lsib
                st.parent[c] = msprime.NULL_NODE
                st.left_sib[c] = msprime.NULL_NODE
                st.right_sib[c] = msprime.NULL_NODE

                # If c is not above a sample then we have nothing to do as we
                # cannot affect the status of any roots.
                if st.above_sample[c]:
                    # Compute the new above sample status for the nodes from
                    # p up to root.
                    v = p
                    above_sample = False
                    while v != msprime.NULL_NODE and not above_sample:
                        above_sample = st.is_sample[v]
                        u = st.left_child[v]
                        while u != msprime.NULL_NODE:
                            above_sample = above_sample or st.above_sample[u]
                            u = st.right_sib[u]
                        st.above_sample[v] = above_sample
                        root = v
                        v = st.parent[v]

                    if not above_sample:
                        # root is no longer above samples. Remove it from the root list.
                        lroot = st.left_sib[root]
                        rroot = st.right_sib[root]
                        st.left_root = msprime.NULL_NODE
                        if lroot != msprime.NULL_NODE:
                            st.right_sib[lroot] = rroot
                            st.left_root = lroot
                        if rroot != msprime.NULL_NODE:
                            st.left_sib[rroot] = lroot
                            st.left_root = rroot
                        st.left_sib[root] = msprime.NULL_NODE
                        st.right_sib[root] = msprime.NULL_NODE

                    # Add c to the root list.
                    # print("Insert ", c, "into root list")
                    if st.left_root != msprime.NULL_NODE:
                        lroot = st.left_sib[st.left_root]
                        if lroot != msprime.NULL_NODE:
                            st.right_sib[lroot] = c
                        st.left_sib[c] = lroot
                        st.left_sib[st.left_root] = c
                    st.right_sib[c] = st.left_root
                    st.left_root = c

            while j < M and edges[in_order[j]].left == st.left:
                p = edges[in_order[j]].parent
                c = edges[in_order[j]].child
                j += 1

                # print("insert ", c, "->", p)
                st.parent[c] = p
                u = st.right_child[p]
                lsib = st.left_sib[c]
                rsib = st.right_sib[c]
                if u == msprime.NULL_NODE:
                    st.left_child[p] = c
                    st.left_sib[c] = msprime.NULL_NODE
                    st.right_sib[c] = msprime.NULL_NODE
                else:
                    st.right_sib[u] = c
                    st.left_sib[c] = u
                    st.right_sib[c] = msprime.NULL_NODE
                st.right_child[p] = c

                if st.above_sample[c]:
                    v = p
                    above_sample = False
                    while v != msprime.NULL_NODE and not above_sample:
                        above_sample = st.above_sample[v]
                        st.above_sample[v] = st.above_sample[v] or st.above_sample[c]
                        root = v
                        v = st.parent[v]
                    # print("root = ", root, st.above_sample[root])

                    if not above_sample:
                        # Replace c with root in root list.
                        # print("replacing", root, "with ", c ," in root list")
                        if lsib != msprime.NULL_NODE:
                            st.right_sib[lsib] = root
                        if rsib != msprime.NULL_NODE:
                            st.left_sib[rsib] = root
                        st.left_sib[root] = lsib
                        st.right_sib[root] = rsib
                        st.left_root = root
                    else:
                        # Remove c from root list.
                        # print("remove ", c ," from root list")
                        st.left_root = msprime.NULL_NODE
                        if lsib != msprime.NULL_NODE:
                            st.right_sib[lsib] = rsib
                            st.left_root = lsib
                        if rsib != msprime.NULL_NODE:
                            st.left_sib[rsib] = lsib
                            st.left_root = rsib

            st.right = sequence_length
            if j < M:
                st.right = min(st.right, edges[in_order[j]].left)
            if k < M:
                st.right = min(st.right, edges[out_order[k]].right)
            assert st.left_root != msprime.NULL_NODE
            while st.left_sib[st.left_root] != msprime.NULL_NODE:
                st.left_root = st.left_sib[st.left_root]
            st.index += 1
            # Add in all the sites
            st.site_list = [
                site for site in self._sites if st.left <= site.position < st.right]
            yield st
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
            ell = self.__lambda[self.__alpha[x] & ((1 << h) - 1)]
            xhat = self.__tau[((self.__beta[x] >> ell) | 1) << ell]
        if j == self.__beta[y]:
            yhat = y
        else:
            ell = self.__lambda[self.__alpha[y] & ((1 << h) - 1)]
            yhat = self.__tau[((self.__beta[y] >> ell) | 1) << ell]
        if self.__pi[xhat] <= self.__pi[yhat]:
            z = xhat
        else:
            z = yhat
        return z


def overlapping_segments(segments):
    """
    Returns an iterator over the (left, right, X) tuples describing the
    distinct overlapping segments in the specified set.
    """
    S = sorted(segments, key=lambda x: x.left)
    n = len(S)
    # Insert a sentinel at the end for convenience.
    S.append(Segment(sys.float_info.max, 0))
    right = S[0].left
    X = []
    j = 0
    while j < n:
        # Remove any elements of X with right <= left
        left = right
        X = [x for x in X if x.right > left]
        if len(X) == 0:
            left = S[j].left
        while j < n and S[j].left == left:
            X.append(S[j])
            j += 1
        j -= 1
        right = min(x.right for x in X)
        right = min(right, S[j + 1].left)
        yield left, right, X
        j += 1

    while len(X) > 0:
        left = right
        X = [x for x in X if x.right > left]
        if len(X) > 0:
            right = min(x.right for x in X)
            yield left, right, X


class Segment(object):
    """
    A class representing a single segment. Each segment has a left and right,
    denoting the loci over which it spans, a node and a next, giving the next
    in the chain.

    The node it records is the *output* node ID.
    """
    def __init__(self, left=None, right=None, node=None, next=None):
        self.left = left
        self.right = right
        self.node = node
        self.next = next

    def __str__(self):
        s = "({}-{}->{}:next={})".format(
            self.left, self.right, self.node, repr(self.next))
        return s

    def __repr__(self):
        return repr((self.left, self.right, self.node))

    def __lt__(self, other):
        return (self.left, self.right, self.node) < (other.left, other.right, self.node)


class Simplifier(object):
    """
    Simplifies a tree sequence to its minimal representation given a subset
    of the leaves.
    """
    def __init__(self, ts, sample, filter_zero_mutation_sites=True):
        self.ts = ts
        self.n = len(sample)
        self.sequence_length = ts.sequence_length
        self.filter_zero_mutation_sites = filter_zero_mutation_sites
        self.num_mutations = ts.num_mutations
        self.input_sites = list(ts.sites())
        self.A_head = [None for _ in range(ts.num_nodes)]
        self.A_tail = [None for _ in range(ts.num_nodes)]
        # Use a head sentinel on the linked lists for convenience.
        for j in range(ts.num_nodes):
            self.A_head[j] = Segment(0, 0, -1)
            self.A_tail[j] = self.A_head[j]
        self.mutation_table = msprime.MutationTable(ts.num_mutations)
        self.node_table = msprime.NodeTable(ts.num_nodes)
        self.edge_table = msprime.EdgeTable(ts.num_edges)
        self.site_table = msprime.SiteTable(ts.num_sites)
        self.mutation_table = msprime.MutationTable(ts.num_mutations)
        self.edge_buffer = []
        self.node_id_map = np.zeros(ts.num_nodes, dtype=np.int32) - 1
        self.mutation_node_map = [-1 for _ in range(self.num_mutations)]
        self.samples = set(sample)
        for sample_id in sample:
            output_id = self.record_node(sample_id, is_sample=True)
            self.add_ancestry(sample_id, 0, self.sequence_length, output_id)
        # We keep a map of input nodes to mutations.
        self.mutation_map = [[] for _ in range(ts.num_nodes)]
        position = ts.tables.sites.position
        site = ts.tables.mutations.site
        node = ts.tables.mutations.node
        for mutation_id in range(ts.num_mutations):
            site_position = position[site[mutation_id]]
            self.mutation_map[node[mutation_id]].append((site_position, mutation_id))

    def record_node(self, input_id, is_sample=False):
        """
        Adds a new node to the output table corresponding to the specified input
        node ID.
        """
        node = self.ts.node(input_id)
        flags = node.flags
        # Need to zero out the sample flag
        flags &= ~msprime.NODE_IS_SAMPLE
        if is_sample:
            flags |= msprime.NODE_IS_SAMPLE
        output_id = self.node_table.add_row(
            flags=flags, time=node.time, population=node.population,
            metadata=node.metadata)
        self.node_id_map[input_id] = output_id
        return output_id

    def flush_edges(self):
        """
        Flush the edges to the output table after sorting and squashing
        any redundant records.
        """
        if len(self.edge_buffer) > 0:
            self.edge_buffer.sort(key=lambda e: (e.child, e.left))
            parent = self.edge_buffer[0].parent
            left = self.edge_buffer[0].left
            right = self.edge_buffer[0].right
            child = self.edge_buffer[0].child
            for e in self.edge_buffer[1:]:
                assert e.parent == parent
                if e.left != right or e.child != child:
                    self.edge_table.add_row(left, right, parent, child)
                    left = e.left
                    child = e.child
                right = e.right
            self.edge_table.add_row(left, right, parent, child)
            self.edge_buffer = []

    def record_edge(self, left, right, parent, child):
        """
        Adds an edge to the output list.
        """
        self.edge_buffer.append(
            msprime.Edge(left=left, right=right, parent=parent, child=child))

    def print_state(self):
        print(".................")
        print("Ancestors: ")
        num_nodes = len(self.A_tail)
        for j in range(num_nodes):
            print("\t", j, "->", end="")
            x = self.A_head[j].next
            while x is not None:
                print("({}-{}->{})".format(x.left, x.right, x.node), end="")
                x = x.next
            print()
        print("Mutation map:")
        for u in range(len(self.mutation_map)):
            v = self.mutation_map[u]
            if len(v) > 0:
                print("\t", u, "->", v)
        print("Node ID map: (input->output)")
        for input_id, output_id in enumerate(self.node_id_map):
            print("\t", input_id, "->", output_id)
        print("Mutation node map")
        for j in range(self.num_mutations):
            print("\t", j, "->", self.mutation_node_map[j])
        print("Output nodes:")
        print(self.node_table)
        print("Output Edges: ")
        print(self.edge_table)
        print("Output sites:")
        print(self.site_table)
        print("Output mutations: ")
        print(self.mutation_table)
        self.check_state()

    def add_ancestry(self, input_id, left, right, node):
        tail = self.A_tail[input_id]
        if tail.right == left and tail.node == node:
            tail.right = right
        else:
            x = Segment(left, right, node)
            tail.next = x
            self.A_tail[input_id] = x

    def merge_labeled_ancestors(self, S, input_id):
        """
        All ancestry segments in S come together into a new parent.
        The new parent must be assigned and any overlapping segments coalesced.
        """
        output_id = self.node_id_map[input_id]
        is_sample = output_id != -1
        if is_sample:
            # Free up the existing ancestry mapping.
            x = self.A_tail[input_id]
            assert x.left == 0 and x.right == self.sequence_length
            self.A_tail[input_id] = self.A_head[input_id]
            self.A_head[input_id].next = None

        prev_right = 0
        for left, right, X in overlapping_segments(S):
            if len(X) == 1:
                ancestry_node = X[0].node
                if is_sample:
                    self.record_edge(left, right, output_id, X[0].node)
                    ancestry_node = output_id
            else:
                if output_id == -1:
                    output_id = self.record_node(input_id)
                ancestry_node = output_id
                for x in X:
                    self.record_edge(left, right, output_id, x.node)
            if is_sample and left != prev_right:
                # Fill in any gaps in the ancestry for the sample
                self.add_ancestry(input_id, prev_right, left, output_id)
            self.add_ancestry(input_id, left, right, ancestry_node)
            prev_right = right

        if is_sample and prev_right != self.sequence_length:
            # If a trailing gap exists in the sample ancestry, fill it in.
            self.add_ancestry(input_id, prev_right, self.sequence_length, output_id)
        self.flush_edges()

    def process_parent_edges(self, edges):
        """
        Process all of the edges for a given parent.
        """
        assert len(set(e.parent for e in edges)) == 1
        parent = edges[0].parent
        S = []
        for edge in edges:
            x = self.A_head[edge.child].next
            while x is not None:
                if x.right > edge.left and edge.right > x.left:
                    y = Segment(max(x.left, edge.left), min(x.right, edge.right), x.node)
                    S.append(y)
                x = x.next
        self.merge_labeled_ancestors(S, parent)
        self.check_state()
        # self.print_state()

    def finalise_sites(self):
        # Build a map from the old mutation IDs to new IDs. Any mutation that
        # has not been mapped to a node in the new tree sequence will be removed.
        mutation_id_map = [-1 for _ in range(self.num_mutations)]
        num_output_mutations = 0

        for site in self.ts.sites():
            num_output_site_mutations = 0
            for mut in site.mutations:
                mapped_node = self.mutation_node_map[mut.id]
                mapped_parent = -1
                if mut.parent != -1:
                    mapped_parent = mutation_id_map[mut.parent]
                if mapped_node != -1:
                    keep = True
                    if mapped_parent == -1 and site.ancestral_state == mut.derived_state:
                        keep = False
                    if keep:
                        mutation_id_map[mut.id] = num_output_mutations
                        num_output_mutations += 1
                        num_output_site_mutations += 1
            output_site = True
            if self.filter_zero_mutation_sites and num_output_site_mutations == 0:
                output_site = False

            if output_site:
                for mut in site.mutations:
                    if mutation_id_map[mut.id] != -1:
                        mapped_parent = -1
                        if mut.parent != -1:
                            mapped_parent = mutation_id_map[mut.parent]
                        self.mutation_table.add_row(
                            site=len(self.site_table),
                            node=self.mutation_node_map[mut.id],
                            parent=mapped_parent,
                            derived_state=mut.derived_state,
                            metadata=mut.metadata)
                self.site_table.add_row(
                    position=site.position, ancestral_state=site.ancestral_state,
                    metadata=site.metadata)

    def map_mutation_nodes(self):
        for input_node in range(len(self.mutation_map)):
            mutations = self.mutation_map[input_node]
            seg = self.A_head[input_node].next
            m_index = 0
            while seg is not None and m_index < len(mutations):
                x, mutation_id = mutations[m_index]
                if seg.left <= x < seg.right:
                    self.mutation_node_map[mutation_id] = seg.node
                    m_index += 1
                elif x >= seg.right:
                    seg = seg.next
                else:
                    assert x < seg.left
                    m_index += 1

    def simplify(self):
        # print("START")
        # self.print_state()
        if self.ts.num_edges > 0:
            all_edges = list(self.ts.edges())
            edges = all_edges[:1]
            for e in all_edges[1:]:
                if e.parent != edges[0].parent:
                    self.process_parent_edges(edges)
                    edges = []
                edges.append(e)
            self.process_parent_edges(edges)
        # self.print_state()
        self.map_mutation_nodes()
        self.finalise_sites()
        ts = msprime.load_tables(
            nodes=self.node_table, edges=self.edge_table,
            sites=self.site_table, mutations=self.mutation_table,
            sequence_length=self.sequence_length)
        return ts, self.node_id_map

    def check_state(self):
        num_nodes = len(self.A_head)
        for j in range(num_nodes):
            head = self.A_head[j]
            tail = self.A_tail[j]
            assert head.left == 0
            assert head.right == 0
            assert head.node == -1
            x = head
            while x.next is not None:
                x = x.next
            assert x == tail
            x = head.next
            while x is not None:
                assert x.left < x.right
                if x.next is not None:
                    assert x.right <= x.next.left
                    # We should also not have any squashable segments.
                    if x.right == x.next.left:
                        assert x.node != x.next.node
                x = x.next


def base64_encode(metadata):
    """
    Returns the specified metadata bytes object encoded as an ASCII-safe
    string.
    """
    return base64.b64encode(metadata).decode('utf8')


if __name__ == "__main__":
    # Simple CLI for running simplifier above.
    ts = msprime.load(sys.argv[1])
    samples = list(map(int, sys.argv[2:]))
    s = Simplifier(ts, samples)
    # s.print_state()
    tss, _ = s.simplify()
    tables = tss.dump_tables()
    print("Output:")
    print(tables.nodes)
    print(tables.edges)
    print(tables.sites)
    print(tables.mutations)
