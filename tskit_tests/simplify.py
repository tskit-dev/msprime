"""
Python implementation of the simplify algorithm.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import sys

import numpy as np

import tskit


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
    def __init__(
            self, ts, sample, reduce_to_site_topology=False, filter_sites=True,
            filter_populations=True, filter_individuals=True):
        self.ts = ts
        self.n = len(sample)
        self.reduce_to_site_topology = reduce_to_site_topology
        self.sequence_length = ts.sequence_length
        self.filter_sites = filter_sites
        self.filter_populations = filter_populations
        self.filter_individuals = filter_individuals
        self.num_mutations = ts.num_mutations
        self.input_sites = list(ts.sites())
        self.A_head = [None for _ in range(ts.num_nodes)]
        self.A_tail = [None for _ in range(ts.num_nodes)]
        self.tables = tskit.TableCollection(sequence_length=ts.sequence_length)
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
        self.position_lookup = None
        if self.reduce_to_site_topology:
            self.position_lookup = np.hstack([[0], position, [self.sequence_length]])

    def record_node(self, input_id, is_sample=False):
        """
        Adds a new node to the output table corresponding to the specified input
        node ID.
        """
        node = self.ts.node(input_id)
        flags = node.flags
        # Need to zero out the sample flag
        flags &= ~tskit.NODE_IS_SAMPLE
        if is_sample:
            flags |= tskit.NODE_IS_SAMPLE
        output_id = self.tables.nodes.add_row(
            flags=flags, time=node.time, population=node.population,
            metadata=node.metadata, individual=node.individual)
        self.node_id_map[input_id] = output_id
        return output_id

    def rewind_node(self, input_id, output_id):
        """
        Remove the mapping for the specified input and output node pair. This is
        done because there are no edges referring to the node.
        """
        assert output_id == len(self.tables.nodes) - 1
        assert output_id == self.node_id_map[input_id]
        self.tables.nodes.truncate(output_id)
        self.node_id_map[input_id] = -1

    def flush_edges(self):
        """
        Flush the edges to the output table after sorting and squashing
        any redundant records.
        """
        num_edges = 0
        for child in sorted(self.edge_buffer.keys()):
            for edge in self.edge_buffer[child]:
                self.tables.edges.add_row(edge.left, edge.right, edge.parent, edge.child)
                num_edges += 1
        self.edge_buffer.clear()
        return num_edges

    def record_edge(self, left, right, parent, child):
        """
        Adds an edge to the output list.
        """
        if self.reduce_to_site_topology:
            X = self.position_lookup
            left_index = np.searchsorted(X, left)
            right_index = np.searchsorted(X, right)
            # Find the smallest site position index greater than or equal to left
            # and right, i.e., slide each endpoint of an interval to the right
            # until they hit a site position. If both left and right map to the
            # the same position then we discard this edge. We also discard an
            # edge if left = 0 and right is less than the first site position.
            if left_index == right_index or (left_index == 0 and right_index == 1):
                return
            # Remap back to zero if the left end maps to the first site.
            if left_index == 1:
                left_index = 0
            left = X[left_index]
            right = X[right_index]
        if child not in self.edge_buffer:
            self.edge_buffer[child] = [tskit.Edge(left, right, parent, child)]
        else:
            last = self.edge_buffer[child][-1]
            if last.right == left:
                last.right = right
            else:
                self.edge_buffer[child].append(tskit.Edge(left, right, parent, child))

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
                    self.record_edge(left, right, output_id, ancestry_node)
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
        if output_id != -1:
            num_edges = self.flush_edges()
            if num_edges == 0 and not is_sample:
                self.rewind_node(input_id, output_id)

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
            if self.filter_sites and num_output_site_mutations == 0:
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

    def finalise_references(self):
        input_populations = self.ts.tables.populations
        population_id_map = np.arange(len(input_populations) + 1, dtype=np.int32)
        # Trick to ensure the null population gets mapped to itself.
        population_id_map[-1] = -1
        input_individuals = self.ts.tables.individuals
        individual_id_map = np.arange(len(input_individuals) + 1, dtype=np.int32)
        # Trick to ensure the null individual gets mapped to itself.
        individual_id_map[-1] = -1

        population_ref_count = np.ones(len(input_populations), dtype=int)
        if self.filter_populations:
            population_ref_count[:] = 0
            population_id_map[:] = -1
        individual_ref_count = np.ones(len(input_individuals), dtype=int)
        if self.filter_individuals:
            individual_ref_count[:] = 0
            individual_id_map[:] = -1

        for node in self.tables.nodes:
            if self.filter_populations and node.population != tskit.NULL:
                population_ref_count[node.population] += 1
            if self.filter_individuals and node.individual != tskit.NULL:
                individual_ref_count[node.individual] += 1

        for input_id, count in enumerate(population_ref_count):
            if count > 0:
                row = input_populations[input_id]
                output_id = self.tables.populations.add_row(metadata=row.metadata)
                population_id_map[input_id] = output_id
        for input_id, count in enumerate(individual_ref_count):
            if count > 0:
                row = input_individuals[input_id]
                output_id = self.tables.individuals.add_row(
                    flags=row.flags, location=row.location, metadata=row.metadata)
                individual_id_map[input_id] = output_id

        # Remap the population ID references for nodes.
        nodes = self.tables.nodes
        nodes.set_columns(
            flags=nodes.flags,
            time=nodes.time,
            metadata=nodes.metadata,
            metadata_offset=nodes.metadata_offset,
            individual=individual_id_map[nodes.individual],
            population=population_id_map[nodes.population])

        # We don't support migrations for now. We'll need to remap these as well.
        assert self.ts.num_migrations == 0

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
        self.finalise_references()
        ts = self.tables.tree_sequence()
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
    ts = tskit.load(sys.argv[1])
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
