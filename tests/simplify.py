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
Python implementation of the simplify algorithm.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import sys

try:
    # We run some tests on the CLI to make sure that we can work in a minimal
    # sense without numpy. We should extract the PythonSimplifier (and other
    # algorithms) out into their own module so we don't pull in numpy for
    # all tests.
    import numpy as np
except ImportError:
    pass

import msprime


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
        self.tables = msprime.TableCollection(sequence_length=ts.sequence_length)
        # We don't touch populations, so add them straigtht in.
        t = ts.tables
        self.tables.populations.set_columns(
            metadata=t.populations.metadata,
            metadata_offset=t.populations.metadata_offset)
        # For now we don't remove individuals that have no nodes referring to
        # them, so we can just copy the table.
        self.tables.individuals.set_columns(
            flags=t.individuals.flags,
            location=t.individuals.location,
            location_offset=t.individuals.location_offset,
            metadata=t.individuals.metadata,
            metadata_offset=t.individuals.metadata_offset)
        self.edge_buffer = {}
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
        output_id = self.tables.nodes.add_row(
            flags=flags, time=node.time, population=node.population,
            metadata=node.metadata, individual=node.individual)
        self.node_id_map[input_id] = output_id
        return output_id

    def flush_edges(self):
        """
        Flush the edges to the output table after sorting and squashing
        any redundant records.
        """
        for child in sorted(self.edge_buffer.keys()):
            for edge in self.edge_buffer[child]:
                self.tables.edges.add_row(edge.left, edge.right, edge.parent, edge.child)
        self.edge_buffer.clear()

    def record_edge(self, left, right, parent, child):
        """
        Adds an edge to the output list.
        """
        if child not in self.edge_buffer:
            self.edge_buffer[child] = [msprime.Edge(left, right, parent, child)]
        else:
            last = self.edge_buffer[child][-1]
            if last.right == left:
                last.right = right
            else:
                self.edge_buffer[child].append(msprime.Edge(left, right, parent, child))

    def print_state(self):
        print(".................")
        print("Ancestors: ")
        num_nodes = len(self.A_tail)
        for j in range(num_nodes):
            print("\t", j, "->", end="")
            x = self.A_head[j]
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
        print("Output:")
        print(self.tables)
        self.check_state()

    def add_ancestry(self, input_id, left, right, node):
        tail = self.A_tail[input_id]
        if tail is None:
            x = Segment(left, right, node)
            self.A_head[input_id] = x
            self.A_tail[input_id] = x
        else:
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
            self.A_tail[input_id] = None
            self.A_head[input_id] = None

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
            x = self.A_head[edge.child]
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
                        self.tables.mutations.add_row(
                            site=len(self.tables.sites),
                            node=self.mutation_node_map[mut.id],
                            parent=mapped_parent,
                            derived_state=mut.derived_state,
                            metadata=mut.metadata)
                self.tables.sites.add_row(
                    position=site.position, ancestral_state=site.ancestral_state,
                    metadata=site.metadata)

    def map_mutation_nodes(self):
        for input_node in range(len(self.mutation_map)):
            mutations = self.mutation_map[input_node]
            seg = self.A_head[input_node]
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
            sequence_length=self.sequence_length,
            **self.tables.asdict())
        return ts, self.node_id_map

    def check_state(self):
        num_nodes = len(self.A_head)
        for j in range(num_nodes):
            head = self.A_head[j]
            tail = self.A_tail[j]
            if head is None:
                assert tail is None
            else:
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
