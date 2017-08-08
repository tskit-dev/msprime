"""
Python versions of the algorithms from the paper.
"""
from __future__ import print_function
from __future__ import division

import sys
import random
import tempfile
import argparse
import heapq
import math

import msprime
import numpy as np
import bintrees

from six.moves import StringIO


class Segment(object):
    """
    A class representing a single segment. Each segment has a left and right,
    denoting the loci over which it spans, a node and a next, giving the next
    in the chain.

    The node it records is the *output* node ID.
    """
    def __init__(self):
        self.left = None
        self.right = None
        self.node = None
        self.prev = None
        self.next = None

    def __str__(self):
        s = "({}-{}->{}: prev={} next={})".format(
            self.left, self.right, self.node, repr(self.prev),
            repr(self.next))
        return s

    def __lt__(self, other):
        return (self.left, self.right, self.node) < (other.left, other.right, self.node)


class Simplifier(object):
    """
    Modified from Simulator().
    """
    def __init__(self, ts, sample):
        sample_size = len(ts.samples())
        self.ts = ts
        self.sample = sample
        self.n = len(sample)
        self.m = ts.sequence_length
        # A maps input node IDs to the extant ancestor chain. Once the algorithm
        # has processed the ancestors, they are are removed from the map.
        self.A = {}
        # Output tables
        self.node_table = msprime.NodeTable(ts.num_nodes)
        self.edgeset_table = msprime.EdgesetTable(ts.num_edgesets)
        # will keep here output sites, and associated mutations,
        # with keys equal to position
        self.sites = {}
        self.last_edgeset = None
        self.num_output_nodes = 0
        self.num_output_mutations = 0
        self.num_output_sites = 0
        # Keep track of then number of segments we alloc and free to ensure we
        # don't leak.
        self.num_used_segments = 0

        j = 0
        for j, sample_id in enumerate(sample):
            # segment label (j) is the output node ID
            x = self.alloc_segment(0, self.m, j)
            # and the label in A is the input node ID
            self.A[sample_id] = x
            self.record_sample_node(sample_id)
            self.record_mutations(input_id = sample_id,
                                  output_id = self.num_output_nodes - 1,
                                  left = 0.0, right = self.m)

        self.S = bintrees.AVLTree()
        self.S[0] = self.n
        self.S[self.m] = -1

    def alloc_segment(self, left, right, node, prev=None, next=None):
        """
        Allocates a new segment with the specified values.
        """
        assert left < right
        s = Segment()
        s.left = left
        s.right = right
        s.node = node
        s.next = next
        s.prev = prev
        self.num_used_segments += 1
        return s

    def free_segment(self, u):
        """
        Frees the specified segment.

        Note: this method is only here to ensure that we are not leaking segments
        in the C implementation.
        """
        self.num_used_segments -= 1

    def record_sample_node(self, input_id):
        """
        Adds a new node to the output table corresponding to the specified input
        node ID (which will be `self.num_output_nodes - 1`).
        """
        # If we were to keep track of the full output to input mapping in M:
        # self.M[self.num_output_nodes] = input_id
        node = self.ts.node(input_id)
        flags = node.flags & msprime.NODE_IS_SAMPLE
        self.node_table.add_row(
            flags=flags, time=node.time, population=node.population)
        self.num_output_nodes += 1

    def check_or_record_node(self, input_id):
        """
        Adds a new node to the output table corresponding to the specified input
        node ID, *unless* the `input_id` corresponds to a sample, which has
        already been added.  In either case, returns the output node ID.

        Since this is not called to add sample nodes, set flags to be not samples here.
        """
        # If we were to keep track of the full output to input mapping in M:
        # self.M[self.num_output_nodes] = input_id
        if input_id not in self.sample:
            node = self.ts.node(input_id)
            flags = node.flags & ~msprime.NODE_IS_SAMPLE
            self.node_table.add_row(
                flags=node.flags, time=node.time, population=node.population)
            self.num_output_nodes += 1
            output_id = self.num_output_nodes - 1
        else:
            output_id = self.sample.index(input_id)
        # print("id", input_id, "maps to", output_id)
        return output_id

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

    def record_mutations(self, input_id, left, right, output_id):
        """
        For each input mutation associated with node `input_id` between `left`
        and `right`, add a new mutation to the output table associated with
        node `output_id`.  Leaves `index` of sites and `site` of mutations
        as `None` because site ordering will be determined on output.
        """
        # inefficiently find all matching mutations
        for site in self.ts.sites():
            if (site.position >= left) and (site.position < right):
                for mutation in site.mutations:
                    if mutation.node == input_id:
                        # insert site in output if it doesn't exist
                        if site.position not in self.sites:
                            new_site = msprime.Site(
                                        position = site.position,
                                        ancestral_state = site.ancestral_state,
                                        index = None,
                                        mutations = [])
                            self.num_output_sites += 1
                        else:
                            new_site = self.sites[site.position]
                        # and insert mutation in output with output_id
                        new_mutation = msprime.Mutation(
                                        site=None,
                                        node = output_id,
                                        derived_state = mutation.derived_state)
                        new_site.mutations.append(new_mutation)
                        self.num_output_mutations += 1

    def update_ancestral_state(self, input_id, left, right):
        """
        This function is called when it is discovered that the unversal MRCA of
        the samples is the input node `input_id` on the segment `[left,
        right)`, which may be different than in the original tree sequence.
        For this reason, the ancestral states of any sites in that region must
        be updated.
        """
        # inefficiently...
        for site in self.ts.sites():
            if (site.position in self.sites) \
                    and (site.position >= left) and (site.position < right):
                # find the most recent mutation on the path from input_id back
                # to the root, if any
                #   DO THIS SOMEHOW:
                new_ancestral_state = self.ts.allele_of_this_individual(input_id)
                self.sites[site.position].ancestral_state = new_ancestral_state

    def segment_chain_str(self, segment):
        u = segment
        s = ""
        while u is not None:
            s += "({0}-{1}->{2})".format(u.left, u.right, u.node)
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
        for k, x in self.S.items():
            print("\t", k, "\t:\t", x)
        print("Output nodes:")
        print(self.node_table)
        print("Output Edgesets: ")
        print(self.edgeset_table)
        print("Output sites and mutations: ")
        for site in self.sites:
            print(site)

    def simplify(self):
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
                    # print("remove edgeset:", edgeset)
                    # self.print_state()
                    for child in edgeset.children:
                        if child in self.A:
                            self.remove_ancestry(edgeset.left, edgeset.right, child, H)
                # print("---- will merge these segments (H):")
                # self.print_heaps(H)
                # print("---- State before merging:")
                # self.print_state()
                self.merge_labeled_ancestors(H, input_id)
                # print("---- merged: ", input_id)
                # self.print_state()
        # print("------ done!")
        # self.print_state()
        # assert self.num_used_segments == 0

        # Flush the last edgeset to the table and create the new tree sequence.
        left, right, parent, children = self.last_edgeset
        self.edgeset_table.add_row(
            left=left, right=right, parent=parent, children=children)
        # construct Site and Mutation tables
        site_table = msprime.SiteTable()
        mutation_table = msprime.MutationTable()
        for k, pos in enumerate(sorted(self.sites.keys())):
            site = self.sites[pos]
            site_table.add_row(position = site.position,
                               ancestral_state = self.ancestral_state)
            for mut in site.mutations:
                mutation_table.add_row(site=k, node=mut.node,
                                       derived_state=mut.derived_state)

        return msprime.load_tables(nodes=self.node_table, edgesets=self.edgeset_table,
                                   sites=site_table, mutations=mutation_table)


    def remove_ancestry(self, left, right, input_id, H):
        """
        Remove the ancestry for the specified input node over the specified interval
        by snipping out elements of the segment chain and modifying any segments
        that overlap the edges. Update the specified heapq H of (x.left, x)
        tuples, where x is the head of a linked list of ancestral segments that we
        remove from the chain for input_id.
        """
        x = self.A[input_id]
        head = x
        last = None
        while x is not None and x.right <= left:
            last = x
            x = x.next
        if x is not None and x.left < left:
            # The left edge of x overhangs. Insert a new segment for the excess.
            y = self.alloc_segment(x.left, left, x.node, x.prev, None)
            x.left = left
            if last is not None:
                last.next = y
            last = y
            if x == head:
                head = last
        if x is not None and x.left < right:
            # x is the first segment within the target interval, so add it to the
            # output heapq.
            x.prev = None
            heapq.heappush(H, (x.left, x))
            # Skip over segments strictly within the interval
            x_prev = x
            while x is not None and x.right <= right:
                x_prev = x
                x = x.next
            if x is not None and x.left < right:
                # We have an overhang on the right hand side. Create a new
                # segment for the overhang and terminate the output chain.
                y = self.alloc_segment(right, x.right, x.node, None, x.next)
                x.right = right
                x.next = None
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
        u = None
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
                    # output node ID
                    u = self.check_or_record_node(input_id)
                assert u is not None
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
                    # all done with this segment, so check ancestral states
                    self.update_ancestral_state(input_id, l, r)
                else:
                    r = l
                    while r < r_max and self.S[r] != len(X):
                        self.S[r] -= len(X) - 1
                        r = self.S.succ_key(r)
                    alpha = self.alloc_segment(l, r, u)
                # Update the heaps and make the record.
                children = []
                for x in X:
                    if x.node is not u:
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
                alpha.prev = z
                z = alpha
                # and copy over any mutations on this segment to the output
                self.record_mutations(input_id, alpha.left, alpha.right, alpha.node)


def run_simplify(args):
    """
    Runs simplify on the tree sequence.
    """
    ts = msprime.load(args.tree_sequence)
    random.seed(args.random_seed)
    sample = random.sample(ts.samples(), args.sample_size)
    s = Simplifier(ts, sample)
    new_ts = s.simplify()
    print("Input:")
    for t in ts.trees():
        print(t)
    print("Output:")
    for t in new_ts.trees():
        print(t)
    # process_trees(new_ts)


def add_simplifier_arguments(parser):
    parser.add_argument("tree_sequence", type=str)
    parser.add_argument("sample_size", type=int)
    parser.add_argument(
        "--random_seed", "-s", type=int, default=1)


def main():
    parser = argparse.ArgumentParser()
    # This is required to get uniform behaviour in Python2 and Python3
    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    simplify_parser = subparsers.add_parser(
        "simplify",
        help="Simplify the tree sequence to fewer samples..")
    add_simplifier_arguments(simplify_parser)
    simplify_parser.set_defaults(runner=run_simplify)

    args = parser.parse_args()
    args.runner(args)


if __name__ == "__main__":
    main()
