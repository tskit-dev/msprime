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
"""
Common code for the msprime test cases.
"""
import random
import unittest

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
        ret.sparse_tree = sparse_tree
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

    def newick(self, root=None, precision=16, node_labels=None):
        if node_labels is None:
            node_labels = {u: str(u + 1) for u in self.sparse_tree.leaves()}
        if root is None:
            root = self.left_root
        return self._build_newick(root, precision, node_labels) + ";"

    def _build_newick(self, node, precision, node_labels):
        label = node_labels.get(node, "")
        if self.left_child[node] == msprime.NULL_NODE:
            s = label
        else:
            s = "("
            for child in self.children(node):
                branch_length = self.sparse_tree.branch_length(child)
                subtree = self._build_newick(child, precision, node_labels)
                s += subtree + ":{0:.{1}f},".format(branch_length, precision)
            s = s[:-1] + label + ")"
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
