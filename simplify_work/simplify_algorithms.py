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

from six.moves import StringIO

from algorithms import *


class SearchablePopulation(object):
    """
    Class representing a population in the simulation.
    """
    def __init__(self, id_):
        self._id = id_
        # we will use this to identify ancestors with the input node IDs:
        #   keys are input node IDs
        #   and values are their head ancestral segments
        self._ancestors = {}

    def __contains__(self, key):
        return key in self._ancestors

    def __getitem__(self, key):
        return self._ancestors[key]

    def __setitem__(self, key, value):
        self._ancestors[key] = value

    def print_state(self):
        print("Population ", self._id)
        for x in self._ancestors:
            s = str(x) + ": "
            u = self._ancestors[x]
            while u is not None:
                s += "({0}-{1}->{2}({3}))".format(
                    u.left, u.right, u.node, u.index)
                u = u.next
            print("\t\t" + s)

    def get_num_ancestors(self):
        return len(self._ancestors)

    def remove(self, index):
        """
        Removes and returns the individual at the specified index.
        """
        return self._ancestors.pop(index)

    def add_with_id(self, ind, individual):
        """
        Inserts the specified individual into this population.
        """
        self._ancestors[ind] = individual

    def __iter__(self):
        return iter(self._ancestors)


class Simplifier(Simulator):
    """
    Modified from Simulator().
    """
    def __init__(self, ts, sample, max_segments=100):
        # since we will have no migration events,
        # need to use one population but copy over population information
        N = 1
        sample_size = len(ts.samples())
        num_loci = ts.sequence_length

        self.ts = ts
        self.n = len(sample)
        self.m = num_loci
        self.max_segments = max_segments
        self.segment_stack = []
        self.segments = [None for j in range(self.max_segments + 1)]
        for j in range(self.max_segments):
            s = Segment(j + 1)
            self.segments[j + 1] = s
            self.segment_stack.append(s)
        # records from `ts` refer to IDs that we must associate with ancestors
        #   this mapping is recorded here
        self.P = [SearchablePopulation(id_) for id_ in range(N)]
        self.C = []
        self.L = FenwickTree(self.max_segments)
        self.S = bintrees.AVLTree()

        # unused stuff included to use other Simulator methods
        self.r = 0.0
        self.migration_matrix = [0.0]
        self.modifier_events = []

        # need to record this to allow for samples not at time 0.0
        self.sample_times = [ts.time(u) for u in ts.samples()]

        # set this as a constant to make code clear below
        self.pop_index = 0
        self.A = self.P[self.pop_index]
        j = 0
        for k in ts.samples():
            if k in sample:
                # segment label (j) is the output node ID
                x = self.alloc_segment(0, self.m, j, self.pop_index)
                self.L.set_value(x.index, self.m - 1)
                # and the label in A is the input node ID
                self.A.add_with_id(k, x)
                j += 1
        self.S[0] = self.n
        self.S[self.m] = -1
        self.t = 0
        # this (w) gives the next output node ID to be assigned
        #   when coalescent events occur:
        self.w = self.n
        self.num_ca_events = 0
        self.num_re_events = 0

    def print_state(self):
        print(".................")
        print("State @ time ", self.t)
        print("Links = ", self.L.get_total())
        print("Node mappings: ")
        self.A.print_state()
        print("Overlap counts", len(self.S))
        for k, x in self.S.items():
            print("\t", k, "\t:\t", x)
        print("Coalescence records: ")
        for rec in self.C:
            print("\t", rec)

    def simplify(self):
        # need to deal with parents in order by birth time-ago
        the_parents = [(parent.time, input_id) for input_id, parent in enumerate(self.ts.nodes())]
        the_parents.sort()
        for parent_time, input_id in the_parents:
            # print("---> doing parent: ", input_id, "at time", parent_time)
            # self.print_state()
            # inefficent way to pull all edges corresponding to a given parent
            edges = [x for x in self.ts.edgesets() if x.parent == input_id]
            if len(edges) > 0:
                self.t = parent_time
                # pull out the ancestry segments that will be merged
                H = self.remove_ancestry(edges)
                # print("---- will merge these segments (H):")
                # self.print_heaps(H)
                # print("---- State before merging:")
                # self.print_state()
                if len(H) > 0:
                    # and merge them: just like merge_ancestors but needs to return the index
                    # of the first segment of the parent to update P with
                    # or returns None if that parent is empty
                    parent = self.merge_labeled_ancestors(H, self.pop_index)
                    if parent is not None:
                        # this replaces pop.add() in merge_ancestors
                        self.A.add_with_id(input_id, parent)
                        # print("---- merged: ", input_id, "->", parent.index)
                    # self.print_state()
        # print("------ done!")
        # self.print_state()

    def get_ancestor(self, u):
        if u in self.A:
            out = self.A[u]
        else:
            out = None
        return out

    def remove_ancestry(self, edges):
        """
        Remove (modifying in place) and return the subset of the ancestors 
        lying within all intervals (left, right) for each of the children
        for each edge in edges. Modified from paint_simplify::remove_paint().
        The output, H, is a heapq of (x.left, x) tuples, where x is the head of
        an linked list of ancestral segments.
        """
        H = []
        for edge in edges:
            # print("remove edge:", edge)
            # self.print_state()
            for child in edge.children:
                if child in self.A:
                    x = self.get_ancestor(child)
                    # y will be the last segment to the left of edge, if any,
                    #   which we may need to make sure links to the next one after
                    y = None
                    # and z will be the first segment after edge, if any
                    z = None
                    # and w will be the previous segment sent to output
                    w = None
                    while x is not None and edge.right > x.left:
                        # print("begin     x: " + x.__str__())
                        # print("begin     y: " + y.__str__())
                        # print("begin     z: " + z.__str__())
                        # print("begin     w: " + w.__str__())
                        # intervals are half-open: [left, right)
                        #  so that the left coordinate is inclusive and the right
                        if edge.left < x.right and edge.right > x.left:
                            # we have overlap
                            seg_right = x.right
                            out_left = max(edge.left, x.left)
                            out_right = min(edge.right, x.right)
                            overhang_left = (x.left < out_left)
                            overhang_right = (x.right > out_right)
                            if overhang_left:
                                # this means x will be the first before removed segment
                                y = x
                                # revise x to be the left part
                                x.right = out_left
                                # the remaining segment will be sent to output
                                # with the previously output segment w as the previous one
                                next_w = self.alloc_segment(
                                    out_left, out_right, x.node, x.population, w, None)
                            else:
                                # remove x, and use it as next_w
                                x.prev = w
                                x.right = out_right
                                next_w = x
                            if w is None:
                                # then we're at the head of an ancestor that we are outputting to H
                                heapq.heappush(H, (next_w.left, next_w))
                            else:
                                w.next = next_w
                            w = next_w
                            if overhang_right:
                                # add new segment for right overhang, which will be the last one
                                # remaining in this ancestor after the removed segment
                                z = self.alloc_segment(
                                    out_right, seg_right, x.node, x.population, y, x.next)
                                # y.next is updated below
                                if x.next is not None:
                                    x.next.prev = z
                                break
                        else:
                            # maybe THIS segment was the first one before edge
                            y = x
                        # move on to the next segment
                        x = x.next
                    # don't do wrap-up if we haven't actually done anything
                    if w is not None:
                        w.next = None
                        if not overhang_right:
                            z = x
                        if y is not None:
                            y.next = z
                        if z is not None:
                            z.prev = y
                        if y is None:
                            # must update P[child]
                            if z is None:
                                self.A.remove(child)
                            else:
                                self.A[child] = z
                    # print("end     x:" + x.__str__())
                    # print("end     y:" + y.__str__())
                    # print("end     z:" + z.__str__())
                    # print("end     w:" + w.__str__())
            # print(" ... state of H while in removing loop ...")
            # self.print_heaps(H)
        return H

    def merge_labeled_ancestors(self, H, pop_id):
        '''
        As merge_ancestors but returning the index of the resulting head ancestry segment.
        '''
        # H is a heapq of (x.left, x) tuples,
        # with x an ancestor, i.e., a list of segments.
        # This will merge everyone in H and add them to population pop_id
        pop = self.P[pop_id]
        defrag_required = False
        coalescence = False
        alpha = None
        z = None
        out = None
        while len(H) > 0:
            # print("LOOP HEAD")
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
                    alpha = self.alloc_segment(
                        x.left, H[0][0], x.node, x.population)
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
                    self.w += 1
                # output node ID
                u = self.w - 1
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
                    alpha = self.alloc_segment(l, r, u, pop_id)
                # Update the heaps and make the record.
                children = []
                for x in X:
                    children.append(x.node)
                    if x.right == r:
                        self.free_segment(x)
                        if x.next is not None:
                            y = x.next
                            heapq.heappush(H, (y.left, y))
                    elif x.right > r:
                        x.left = r
                        heapq.heappush(H, (x.left, x))
                self.C.append((l, r, u, children, self.t))

            # loop tail; update alpha and integrate it into the state.
            if alpha is not None:
                if z is None:
                    # the only place where this differs from merge_ancestors():
                    out = alpha
                    # pop.add(alpha)
                    self.L.set_value(alpha.index, alpha.right - alpha.left - 1)
                else:
                    defrag_required |= (
                        z.right == alpha.left and z.node == alpha.node)
                    z.next = alpha
                    self.L.set_value(alpha.index, alpha.right - z.right)
                alpha.prev = z
                z = alpha
        if defrag_required:
            self.defrag_segment_chain(z)
        if coalescence:
            self.defrag_breakpoints()
        return out

    def write_text(self, nodes_file, edgesets_file):
        """
        Writes the records out as text.  Modified to allow samples from nonzero times.
        """
        num_nodes = max(r[2] for r in self.C) + 1
        time = self.sample_times + [0 for _ in range(num_nodes - self.n)]
        print("is_sample\ttime", file=nodes_file)
        print("left\tright\tparent\tchildren", file=edgesets_file)
        for left, right, u, children, t in self.C:
            time[u] = t
            print(
                left, right, u, ",".join(str(c) for c in sorted(children)),
                sep="\t", file=edgesets_file)
        for u in range(num_nodes):
            print(
                int(u < self.n), time[u], sep="\t", file=nodes_file)



def run_simplify(args):
    """
    Runs simplify on the tree sequence.
    """
    ts = msprime.load(args.tree_sequence)
    sample = random.sample(ts.samples(), args.sample_size)
    random.seed(args.random_seed)
    s = Simplifier(ts, sample)
    s.simplify()
    nodes_file = StringIO()
    edgesets_file = StringIO()
    s.write_text(nodes_file, edgesets_file)
    nodes_file.seek(0)
    edgesets_file.seek(0)
    print(nodes_file.getvalue())
    print(edgesets_file.getvalue())
    new_ts = msprime.load_text(nodes_file, edgesets_file)
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
