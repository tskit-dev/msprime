#
# Copyright (C) 2015-2017 University of Oxford
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
import heapq
import random
import sys
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
    def __init__(self):
        self.parent = {}
        self.children = {}
        self.time = {}
        self.left = 0
        self.right = 0
        self.root = 0
        self.index = -1
        self.sample_size = 0
        # We need a sites function, so this name is taken.
        self.site_list = []

    @classmethod
    def from_sparse_tree(self, sparse_tree):
        ret = PythonSparseTree()
        ret.root = sparse_tree.get_root()
        ret.sample_size = sparse_tree.get_sample_size()
        ret.left, ret.right = sparse_tree.get_interval()
        ret.site_list = list(sparse_tree.sites())
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
        assert ret == sparse_tree
        return ret

    def _preorder_nodes(self, u, l):
        l.append(u)
        if u in self.children:
            for c in self.children[u]:
                self._preorder_nodes(c, l)

    def _postorder_nodes(self, u, l):
        if u in self.children:
            for c in self.children[u]:
                self._postorder_nodes(c, l)
        l.append(u)

    def _inorder_nodes(self, u, l):
        if u in self.children:
            self._inorder_nodes(self.children[u][0], l)
            l.append(u)
            self._inorder_nodes(self.children[u][1], l)
        else:
            l.append(u)

    def _levelorder_nodes(self, u, l, level):
        l[level].append(u) if level < len(l) else l.append([u])
        if u in self.children:
            for c in self.children[u]:
                self._levelorder_nodes(c, l, level + 1)

    def nodes(self, root=None, order="preorder"):
        u = root
        l = []
        if root is None:
            u = self.root
        if order == "preorder":
            self._preorder_nodes(u, l)
            return iter(l)
        elif order == "inorder":
            self._inorder_nodes(u, l)
            return iter(l)
        elif order == "postorder":
            self._postorder_nodes(u, l)
            return iter(l)
        elif order == "levelorder" or order == "breadthfirst":
            # Returns nodes in their respective levels
            # Nested list comprehension flattens l in order
            self._levelorder_nodes(u, l, 0)
            return iter([i for level in l for i in level])
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

    def sites(self):
        return iter(self.site_list)


class PythonTreeSequence(object):
    """
    A python implementation of the TreeSequence object.
    """
    def __init__(self, tree_sequence, breakpoints=None):
        self._tree_sequence = tree_sequence
        self._sample_size = tree_sequence.get_sample_size()
        self._breakpoints = breakpoints
        self._sites = []
        _Site = collections.namedtuple(
            "Site",
            ["position", "ancestral_state", "index", "mutations"])
        _Mutation = collections.namedtuple(
            "Mutation",
            ["site", "node", "derived_state"])
        for j in range(tree_sequence.get_num_sites()):
            pos, ancestral_state, mutations, index = tree_sequence.get_site(j)
            self._sites.append(_Site(
                position=pos, ancestral_state=ancestral_state, index=index,
                mutations=[_Mutation(*mut) for mut in mutations]))

    def _diffs(self):
        M = self._tree_sequence.get_num_edgesets()
        records = [self._tree_sequence.get_record(j) for j in range(M)]
        l = [record[0] for record in records]
        r = [record[1] for record in records]
        u = [record[2] for record in records]
        c = [record[3] for record in records]
        t = [record[4] for record in records]
        I = sorted(range(M), key=lambda j: (l[j], t[j]))
        O = sorted(range(M), key=lambda j: (r[j], -t[j]))
        j = 0
        k = 0
        while j < M:
            r_out = []
            r_in = []
            x = l[I[j]]
            while r[O[k]] == x:
                h = O[k]
                r_out.append((u[h], c[h], t[h]))
                k += 1
            while j < M and l[I[j]] == x:
                h = I[j]
                r_in.append((u[h], c[h], t[h]))
                j += 1
            yield r[O[k]] - x, r_out, r_in

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
        M = self._tree_sequence.get_num_edgesets()
        records = [self._tree_sequence.get_record(j) for j in range(M)]
        l = [record[0] for record in records]
        r = [record[1] for record in records]
        u = [record[2] for record in records]
        c = [record[3] for record in records]
        t = [record[4] for record in records]
        I = sorted(range(M), key=lambda j: (l[j], t[j]))
        O = sorted(range(M), key=lambda j: (r[j], -t[j]))
        j = 0
        k = 0
        st = PythonSparseTree()
        st.sample_size = self._tree_sequence.get_sample_size()
        st.left = 0
        st.time = {j: 0 for j in range(st.sample_size)}
        while j < M:
            x = l[I[j]]
            while r[O[k]] == x:
                h = O[k]
                del st.children[u[h]]
                del st.time[u[h]]
                for q in c[h]:
                    del st.parent[q]
                k += 1
            while j < M and l[I[j]] == x:
                h = I[j]
                st.children[u[h]] = c[h]
                st.time[u[h]] = t[h]
                for q in c[h]:
                    st.parent[q] = u[h]
                j += 1
            st.left = x
            st.right = r[O[k]]
            # Insert the root
            root = 0
            while root in st.parent:
                root = st.parent[root]
            st.root = root
            st.index += 1
            # Add in all the sites
            st.site_list = [
                site for site in self._sites if st.left <= site.position < st.right]
            yield st


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
        s = "({}-{}->{}:{}:: next={})".format(
            self.left, self.right, self.node, self.mutations, repr(self.next))
        return s

    def __lt__(self, other):
        return (self.left, self.right, self.node) < (other.left, other.right, self.node)


class SortedMap(dict):
    """
    Simple implementation of a sorted mapping. Based on the API for bintrees.AVLTree.
    We don't use bintrees here because it is not available on Windows.
    """
    def floor_key(self, k):
        ret = None
        for key in sorted(self.keys()):
            if key <= k:
                ret = key
            if key > k:
                break
            ret = key
        return ret

    def succ_key(self, k):
        ret = None
        for key in sorted(self.keys()):
            if key >= k:
                ret = key
            if key > k:
                break
        return ret


class Simplifier(object):
    """
    Simplifies a tree sequence to its minimal representation given a subset
    of the leaves.
    """
    def __init__(self, ts, sample):
        self.ts = ts
        self.n = len(sample)
        self.m = ts.sequence_length
        self.input_sites = list(ts.sites())
        # A maps input node IDs to the extant ancestor chain. Once the algorithm
        # has processed the ancestors, they are are removed from the map.
        self.A = {}
        # Output tables
        self.node_table = msprime.NodeTable(ts.num_nodes)
        self.edgeset_table = msprime.EdgesetTable(ts.num_edgesets)
        self.site_table = msprime.SiteTable(max(1, ts.num_sites))
        self.mutation_table = msprime.MutationTable(max(1, ts.num_mutations))
        self.last_edgeset = None
        self.num_output_nodes = 0
        self.output_sites = {}
        # Keep track of then number of segments we alloc and free to ensure we
        # don't leak.
        self.num_used_segments = 0
        for j, sample_id in enumerate(sample):
            # segment label (j) is the output node ID
            x = self.alloc_segment(0, self.m, j)
            # and the label in A is the input node ID
            self.A[sample_id] = x
            self.record_node(sample_id)
        self.S = SortedMap()
        self.S[0] = self.n
        self.S[self.m] = -1
        # We keep a sorted map of mutations for each input node.
        self.mutation_map = [SortedMap() for _ in range(ts.num_nodes)]
        for site in self.ts.sites():
            for mut in site.mutations:
                self.mutation_map[mut.node][site.position] = mut

    def get_mutations(self, input_id, left, right):
        """
        Returns all mutations for the specified input ID over the specified
        interval.
        """
        mutations = self.mutation_map[input_id]
        ret = SortedMap()
        pos = mutations.succ_key(left)
        while pos is not None and pos < right:
            mut = mutations.pop(pos)
            ret[pos] = mut
            assert left <= pos < right
            pos = mutations.succ_key(pos)
        # print("GET_MUTATIONS", input_id, left, right, "::", ret)
        return ret

    def alloc_segment(self, left, right, node, next=None):
        """
        Allocates a new segment with the specified values.
        """
        s = Segment(left, right, node, next)
        s.mutations = SortedMap()
        self.num_used_segments += 1
        return s

    def free_segment(self, u):
        """
        Frees the specified segment.

        Note: this method is only here to ensure that we are not leaking segments
        in the C implementation.
        """
        self.num_used_segments -= 1

    def record_node(self, input_id):
        """
        Adds a new node to the output table corresponding to the specified input
        node ID.
        """
        node = self.ts.node(input_id)
        self.node_table.add_row(
            flags=node.flags, time=node.time, population=node.population)
        self.num_output_nodes += 1

    def record_edgeset(self, left, right, parent, children):
        """
        Adds an edgeset to the output list. This method used the ``last_edgeset``
        variable to check for adjacent records that may be squashed. Thus, the
        last edgeset will not be entered in the table, which must be done manually.
        """
        sorted_children = tuple(sorted(children))
        if self.last_edgeset is None:
            self.last_edgeset = left, right, parent, sorted_children
        else:
            last_left, last_right, last_parent, last_children = self.last_edgeset
            squash_condition = (
                last_parent == parent and
                last_children == sorted_children and
                last_right == left)
            if squash_condition:
                self.last_edgeset = last_left, right, parent, sorted_children
            else:
                # Flush the last edgeset
                self.edgeset_table.add_row(
                    left=last_left, right=last_right, parent=last_parent,
                    children=last_children)
                self.last_edgeset = left, right, parent, sorted_children

    def segment_chain_str(self, segment):
        u = segment
        s = ""
        while u is not None:
            s += "({0}-{1}->{2}:{3})".format(
                    u.left, u.right, u.node, u.mutations)
            u = u.next
        return s

    def print_heaps(self, L):
        copy = list(L)
        ordered = [heapq.heappop(copy) for _ in L]
        print("H = ")
        for l, x in ordered:
            print("\t", l, ":", self.segment_chain_str(x))

    def print_state(self):
        print(".................")
        print("Ancestors: ", len(self.A))
        for x in self.A.keys():
            s = str(x) + ": " + self.segment_chain_str(self.A[x])
            print("\t\t" + s)
        print("Overlap counts", len(self.S))
        for k in sorted(self.S.keys()):
            x = self.S[k]
            print("\t", k, "\t:\t", x)
        print("Mutation map:")
        for u in range(len(self.mutation_map)):
            v = self.mutation_map[u]
            if len(v) > 0:
                print("\t", u, "->", v)

        # print("Output nodes:")
        # print(self.node_table)
        # print("Output Edgesets: ")
        # print(self.edgeset_table)
        # print("Output sites and mutations: ")
        # for site in self.sites:
        #     print(site)

    def simplify(self):
        # print("START")
        # self.print_state()
        the_parents = [
            (node.time, input_id) for input_id, node in enumerate(self.ts.nodes())]
        # need to deal with parents in order by birth time-ago
        the_parents.sort()
        for time, input_id in the_parents:
            # print()
            # print("---> doing parent: ", input_id, "at time", time)
            # self.print_state()
            if len(self.A) == 0:
                break
            # inefficent way to pull all edges corresponding to a given parent
            edgesets = [x for x in self.ts.edgesets() if x.parent == input_id]
            # print("edgesets = ", edgesets)
            if len(edgesets) > 0:
                # pull out the ancestry segments that will be merged
                # print("before = ")
                H = []
                for edgeset in edgesets:
                    for child in edgeset.children:
                        if child in self.A:
                            self.remove_ancestry(edgeset.left, edgeset.right, child, H)
                self.check_state()
                # print("---- will merge these segments (H):")
                # self.print_heaps(H)
                # print("---- State before merging:")
                # self.print_state()
                self.merge_labeled_ancestors(H, input_id)
                # print("MERGE DONE")
                # print("---- merged: ", input_id, "->", parent.index)
                # self.print_state()
                self.check_state()
        # print("------ done!")
        # self.print_state()
        extant_segments = 0
        for x in self.A.values():
            while x is not None:
                extant_segments += 1
                x = x.next
        # assert self.num_used_segments == extant_segments

        # Flush the last edgeset to the table and create the new tree sequence.
        left, right, parent, children = self.last_edgeset
        self.edgeset_table.add_row(
            left=left, right=right, parent=parent, children=children)

        # Add in the sites and mutations.
        for j, position in enumerate(sorted(self.output_sites.keys())):
            site = self.output_sites[position]
            self.site_table.add_row(
                position=site.position, ancestral_state=site.ancestral_state)
            # Order the mutations within a site by node for compatability with
            # existing implementation. Since nodes are ordered by time this is
            # equivalent to sorting by time.
            for mutation in site.mutations:
                self.mutation_table.add_row(
                    site=j, node=mutation.node, derived_state=mutation.derived_state)

#         print(self.node_table)
#         print(self.edgeset_table)
#         print(self.site_table)
#         print(self.mutation_table)
        return msprime.load_tables(
            nodes=self.node_table, edgesets=self.edgeset_table,
            sites=self.site_table, mutations=self.mutation_table)

    def record_mutation(self, node, mutation):
        position = self.input_sites[mutation.site].position
        if position not in self.output_sites:
            site = msprime.Site(
                position=position, ancestral_state="0", mutations=[], index=None)
            self.output_sites[position] = site
        else:
            site = self.output_sites[position]
        site.mutations.append(
            msprime.Mutation(site=None, node=node, derived_state=mutation.derived_state))

    def remove_ancestry(self, left, right, input_id, H):
        """
        Remove the ancestry for the specified input node over the specified interval
        by snipping out elements of the segment chain and modifying any segments
        that overlap the edges. Update the specified heapq H of (x.left, x)
        tuples, where x is the head of a linked list of ancestral segments that we
        remove from the chain for input_id.
        """
        head = self.A[input_id]
        # Add mutations to the segments
        x = head
        while x is not None:
            mutations = self.get_mutations(input_id, x.left, x.right)
            for pos, mut in mutations.items():
                x.mutations[self.input_sites[mut.site].position] = mut
            x = x.next
        # print("REMOVE ANCESTRY:", input_id, left, right)
        x = head
        last = None
        # Skip the leading segments before left.
        while x is not None and x.right <= left:
            last = x
            x = x.next
        if x is not None and x.left < left:
            # The left edge of x overhangs. Insert a new segment for the excess.
            y = self.alloc_segment(x.left, left, x.node, None)
            x.left = left
            if last is not None:
                last.next = y
            last = y
            if x == head:
                head = last
            # Remove the mutations in x with pos < left and add to y.
            pos = x.mutations.succ_key(-1)
            while pos is not None and pos < left:
                mut = x.mutations.pop(pos)
                y.mutations[pos] = mut
                pos = x.mutations.succ_key(pos)

        if x is not None and x.left < right:
            # x is the first segment within the target interval, so add it to the
            # output heapq.
            heapq.heappush(H, (x.left, x))
            # Skip over segments strictly within the interval
            while x is not None and x.right <= right:
                x_prev = x
                x = x.next
            if x is not None and x.left < right:
                # We have an overhang on the right hand side. Create a new
                # segment for the overhang and terminate the output chain.
                y = self.alloc_segment(right, x.right, x.node, x.next)
                x.right = right
                x.next = None
                # Remove the mutations in x with pos >= right and add to y.
                pos = x.mutations.succ_key(right)
                while pos is not None:
                    mut = x.mutations.pop(pos)
                    y.mutations[pos] = mut
                    pos = x.mutations.succ_key(pos)
                x = y
            elif x_prev is not None:
                x_prev.next = None
        # x is the first segment in the new chain starting after right.
        if last is None:
            head = x
        else:
            last.next = x
        if head is None:
            del self.A[input_id]
        else:
            self.A[input_id] = head

    def merge_labeled_ancestors(self, H, input_id):
        '''
        All ancestry segments in H come together into a new parent.
        The new parent must be assigned;
        any overlapping segments coalesced;
        and node IDs in the mutation table remapped.
        '''
        # H is a heapq of (x.left, x) tuples,
        # with x an ancestor, i.e., a list of segments.
        coalescence = False
        alpha = None
        z = None
        while len(H) > 0:
            # self.print_heaps(H)
            alpha = None
            l = H[0][0]
            X = []
            r_max = self.m + 1
            while len(H) > 0 and H[0][0] == l:
                x = heapq.heappop(H)[1]
                X.append(x)
                r_max = min(r_max, x.right)
            if len(H) > 0:
                r_max = min(r_max, H[0][0])
            if len(X) == 1:
                x = X[0]
                if len(H) > 0 and H[0][0] < x.right:
                    alpha = self.alloc_segment(x.left, H[0][0], x.node)
                    x.left = H[0][0]
                    heapq.heappush(H, (x.left, x))
                else:
                    if x.next is not None:
                        y = x.next
                        heapq.heappush(H, (y.left, y))
                    alpha = x
                    alpha.next = None
            else:
                if not coalescence:
                    coalescence = True
                    self.record_node(input_id)
                # output node ID
                u = self.num_output_nodes - 1
                # We must also break if the next left value is less than
                # any of the right values in the current overlap set.
                if l not in self.S:
                    j = self.S.floor_key(l)
                    self.S[l] = self.S[j]
                if r_max not in self.S:
                    j = self.S.floor_key(r_max)
                    self.S[r_max] = self.S[j]
                # Update the number of extant segments.
                if self.S[l] == len(X):
                    self.S[l] = 0
                    r = self.S.succ_key(l)
                else:
                    r = l
                    while r < r_max and self.S[r] != len(X):
                        self.S[r] -= len(X) - 1
                        r = self.S.succ_key(r)
                    alpha = self.alloc_segment(l, r, u)
                # Update the heaps and make the record.
                children = []
                for x in X:
                    # Record and remove mutations for the coalescing segment
                    pos = x.mutations.succ_key(-1)
                    while pos is not None and pos < r:
                        mut = x.mutations.pop(pos)
                        self.record_mutation(x.node, mut)
                        pos = x.mutations.succ_key(pos)
                    children.append(x.node)
                    if x.right == r:
                        self.free_segment(x)
                        if x.next is not None:
                            y = x.next
                            heapq.heappush(H, (y.left, y))
                    elif x.right > r:
                        x.left = r
                        heapq.heappush(H, (x.left, x))
                self.record_edgeset(l, r, u, children)

            # loop tail; update alpha and integrate it into the state.
            if alpha is not None:
                if z is None:
                    # Add a new mapping for the input_id to the segment chain starting
                    # with alpha.
                    self.A[input_id] = alpha
                else:
                    z.next = alpha
                z = alpha

    def check_state(self):
        # print("CHECK_STATE")
        # self.print_state()
        for input_id, x in self.A.items():
            while x is not None:
                assert x.left < x.right
                for pos, mut in x.mutations.items():
                    assert pos == self.input_sites[mut.site].position
                    assert x.left <= pos < x.right
                if x.next is not None:
                    assert x.right <= x.next.left
                x = x.next


if __name__ == "__main__":
    # Simple CLI for running simplifier above.
    ts = msprime.load(sys.argv[1])
    samples = list(map(int, sys.argv[2:]))
    s = Simplifier(ts, samples)
    tss = s.simplify()
    tables = tss.dump_tables()
    print("Output:")
    print(tables.nodes)
    print(tables.edgesets)
    print(tables.sites)
    print(tables.mutations)
