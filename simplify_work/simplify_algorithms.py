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


class Simplifier(object):
    """
    Modified from Simulator().
    """
    def __init__(self, ts, sample, max_segments=100):
        sample_size = len(ts.samples())
        self.ts = ts
        self.n = len(sample)
        self.m = ts.sequence_length
        self.max_segments = max_segments
        self.segment_stack = []
        self.segments = [None for j in range(self.max_segments + 1)]
        for j in range(self.max_segments):
            s = Segment(j + 1)
            self.segments[j + 1] = s
            self.segment_stack.append(s)
        # A maps input node IDs to the extant ancestor chain. Once the algorithm
        # has processed the ancestors, they are are removed from the map.
        self.A = {}
        # M maps output node IDs to their corresponding input nodes
        self.M = [-1 for _ in range(ts.num_nodes)]
        # Output edgesets
        self.E = []
        self.S = bintrees.AVLTree()

        j = 0
        for j, sample_id in enumerate(sample):
            # segment label (j) is the output node ID
            x = self.alloc_segment(0, self.m, j)
            # and the label in A is the input node ID
            self.A[sample_id] = x
            self.M[j] = sample_id

        self.S[0] = self.n
        self.S[self.m] = -1
        # this (w) gives the next output node ID to be assigned
        #   when coalescent events occur:
        self.w = self.n

    def alloc_segment(self, left, right, node, prev=None, next=None):
        """
        Pops a new segment off the stack and sets its properties.
        """
        s = self.segment_stack.pop()
        s.left = left
        s.right = right
        s.node = node
        s.next = next
        s.prev = prev
        return s

    def free_segment(self, u):
        """
        Frees the specified segment making it ready for reuse and
        setting its weight to zero.
        """
        self.segment_stack.append(u)

    def print_heaps(self, L):
        copy = list(L)
        ordered = [heapq.heappop(copy) for _ in L]
        print("L = ")
        for l, x in ordered:
            print("\t", l, ":", end="")
            u = x
            s = ""
            while u is not None:
                s += "({0}-{1}->{2}({3}))".format(
                    u.left, u.right, u.node, u.index)
                u = u.next
            print(s)

    def print_state(self):
        print(".................")
        print("Ancestors: ", len(self.A))
        for x in self.A.keys():
            s = str(x) + ": "
            u = self.A[x]
            while u is not None:
                s += "({0}-{1}->{2}({3}))".format(
                    u.left, u.right, u.node, u.index)
                u = u.next
            print("\t\t" + s)
        print("Node mappings: (output->input)")
        for j in range(self.w):
            print("\t", j, self.M[j])
        print("Overlap counts", len(self.S))
        for k, x in self.S.items():
            print("\t", k, "\t:\t", x)
        print("Output Edgesets: ")
        for e in self.E:
            print("\t", e)

    def simplify(self):
        the_parents = [
            (node.time, input_id) for input_id, node in enumerate(self.ts.nodes())]
        # need to deal with parents in order by birth time-ago
        the_parents.sort()
        for time, input_id in the_parents:
            # print("---> doing parent: ", input_id, "at time", time)
            # self.print_state()
            if len(self.A) == 0:
                break
            # inefficent way to pull all edges corresponding to a given parent
            edgesets = [x for x in self.ts.edgesets() if x.parent == input_id]
            # print("edgesets = ", edgesets)
            if len(edgesets) > 0:
                # pull out the ancestry segments that will be merged
                H = self.remove_ancestry(edgesets)
                # print("---- will merge these segments (H):")
                # self.print_heaps(H)
                # print("---- State before merging:")
                # self.print_state()
                self.merge_labeled_ancestors(H, input_id)
                # print("---- merged: ", input_id, "->", parent.index)
                # self.print_state()
        # print("------ done!")
        # self.print_state()

    def remove_ancestry(self, edgesets):
        """
        Remove (modifying in place) and return the subset of the ancestors
        lying within all intervals (left, right) for each of the children
        for each edgeset in edgesets. Modified from paint_simplify::remove_paint().
        The output, H, is a heapq of (x.left, x) tuples, where x is the head of
        an linked list of ancestral segments.
        """
        H = []
        for edgeset in edgesets:
            # print("remove edge:", edge)
            # self.print_state()
            for child in edgeset.children:
                if child in self.A:
                    x = self.A[child]
                    # y will be the last segment to the left of edgeset, if any,
                    #   which we may need to make sure links to the next one after
                    y = None
                    # and z will be the first segment after edgeset, if any
                    z = None
                    # and w will be the previous segment sent to output
                    w = None
                    while x is not None and edgeset.right > x.left:
                        # print("begin     x: " + x.__str__())
                        # print("begin     y: " + y.__str__())
                        # print("begin     z: " + z.__str__())
                        # print("begin     w: " + w.__str__())
                        # intervals are half-open: [left, right)
                        #  so that the left coordinate is inclusive and the right
                        if edgeset.left < x.right and edgeset.right > x.left:
                            # we have overlap
                            seg_right = x.right
                            out_left = max(edgeset.left, x.left)
                            out_right = min(edgeset.right, x.right)
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
                                    out_left, out_right, x.node, w, None)
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
                                    out_right, seg_right, x.node, y, x.next)
                                # y.next is updated below
                                if x.next is not None:
                                    x.next.prev = z
                                break
                        else:
                            # maybe THIS segment was the first one before edgeset
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
                                del self.A[child]
                            else:
                                self.A[child] = z
                    # print("end     x:" + x.__str__())
                    # print("end     y:" + y.__str__())
                    # print("end     z:" + z.__str__())
                    # print("end     w:" + w.__str__())
            # print(" ... state of H while in removing loop ...")
            # self.print_heaps(H)
        return H

    def merge_labeled_ancestors(self, H, input_id):
        '''
        As merge_ancestors but returning the index of the resulting head ancestry segment.
        '''
        # H is a heapq of (x.left, x) tuples,
        # with x an ancestor, i.e., a list of segments.
        # This will merge everyone in H and add them to population pop_id
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
                    # Allocate a new output node ID and map it back to the input ID.
                    self.M[self.w] = input_id
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
                    alpha = self.alloc_segment(l, r, u)
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
                self.E.append((l, r, u, children))

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

    def write_text(self, nodes_file, edgesets_file):
        """
        Writes the records out as text.  Modified to allow samples from nonzero times.
        """
        num_nodes = self.w
        input_nodes = [self.ts.node(self.M[u]) for u in range(num_nodes)]
        print("is_sample\tpopulation\ttime", file=nodes_file)
        for node in input_nodes:
            print(node.flags, node.population, node.time, sep="\t", file=nodes_file)

        print("left\tright\tparent\tchildren", file=edgesets_file)
        # collapse adjacent identical ones
        left, right, parent, _, = self.E[0]
        children = sorted(self.E[0][3])
        k = 1
        while k < len(self.E):
            nleft, nright, nparent, _ = self.E[k]
            nchildren = sorted(self.E[k][3])
            if (right == nleft and len(children) == len(nchildren) and parent == nparent and
                all([a == b for a, b in zip(children, nchildren)])):
                    # squash this record into the last
                    right = nright
            else:
                print(
                    left, right, parent, ",".join(str(c) for c in children),
                    sep="\t", file=edgesets_file)
                left, right, parent, children = nleft, nright, nparent, nchildren
                children = nchildren
            k += 1
        print(
            left, right, parent, ",".join(str(c) for c in children),
            sep="\t", file=edgesets_file)


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
