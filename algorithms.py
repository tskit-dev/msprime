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

import bintrees
import msprime
import numpy as np
import statsmodels.api as sm
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot


class FenwickTree(object):
    """
    A Fenwick Tree to represent cumulative frequency tables over
    integers. Each index from 1 to max_index initially has a
    zero frequency.

    This is an implementation of the Fenwick tree (also known as a Binary
    Indexed Tree) based on "A new data structure for cumulative frequency
    tables", Software Practice and Experience, Vol 24, No 3, pp 327 336 Mar
    1994. This implementation supports any non-negative frequencies, and the
    search procedure always returns the smallest index such that its cumulative
    frequency <= f. This search procedure is a slightly modified version of
    that presented in Tech Report 110, "A new data structure for cumulative
    frequency tables: an improved frequency-to-symbol algorithm." available at
    https://www.cs.auckland.ac.nz/~peter-f/FTPfiles/TechRep110.ps
    """
    def __init__(self, max_index):
        assert max_index > 0
        self.__max_index = max_index
        self.__tree = [0 for j in range(max_index + 1)]
        # Compute the binary logarithm of max_index
        u = self.__max_index
        while u != 0:
            self.__log_max_index = u
            u -= (u & -u)

    def get_total(self):
        """
        Returns the total cumulative frequency over all indexes.
        """
        return self.get_cumulative_frequency(self.__max_index)

    def increment(self, index, v):
        """
        Increments the frequency of the specified index by the specified
        value.
        """
        assert 0 < index <= self.__max_index
        j = index
        while j <= self.__max_index:
            self.__tree[j] += v
            j += (j & -j)

    def set_value(self, index, v):
        """
        Sets the frequency at the specified index to the specified value.
        """
        f = self.get_frequency(index)
        self.increment(index, v - f)

    def get_cumulative_frequency(self, index):
        """
        Returns the cumulative frequency of the specified index.
        """
        assert 0 < index <= self.__max_index
        j = index
        s = 0
        while j > 0:
            s += self.__tree[j]
            j -= (j & -j)
        return s

    def get_frequency(self, index):
        """
        Returns the frequency of the specified index.
        """
        assert 0 < index <= self.__max_index
        j = index
        v = self.__tree[j]
        p = j & (j - 1)
        j -= 1
        while p != j:
            v -= self.__tree[j]
            j = j & (j - 1)
        return v

    def find(self, v):
        """
        Returns the smallest index with cumulative sum >= v.
        """
        j = 0
        s = v
        half = self.__log_max_index
        while half > 0:
            # Skip non-existant entries
            while j + half > self.__max_index:
                half >>= 1
            k = j + half
            if s > self.__tree[k]:
                j = k
                s -= self.__tree[j]
            half >>= 1
        return j + 1


class Segment(object):
    """
    A class representing a single segment. Each segment has a left
    and right, denoting the loci over which it spans, a node and a
    next, giving the next in the chain.
    """
    def __init__(self, index):
        self.left = None
        self.right = None
        self.node = None
        self.prev = None
        self.next = None
        self.population = None
        self.index = index

    def __str__(self):
        s = "({0}:{1}-{2}->{3}: prev={4} next={5})".format(
            self.index, self.left, self.right, self.node, repr(self.prev),
            repr(self.next))
        return s


class Population(object):
    """
    Class representing a population in the simulation.
    """
    def __init__(self, id_):
        self._id = id_
        self._start_time = 0
        self._start_size = 1.0
        self._growth_rate = 0
        # We'd like to use an AVLTree here for P but the API doesn't quite
        # do what we need. Lists are inefficient here and should not be
        # used in a real implementation.
        self._ancestors = []

    def print_state(self):
        print("Population ", self._id)
        print("\tstart_size = ", self._start_size)
        print("\tgrowth_rate = ", self._growth_rate)
        print("\tAncestors: ", len(self._ancestors))
        for u in self._ancestors:
            s = ""
            while u is not None:
                s += "({0}-{1}->{2}({3}))".format(
                    u.left, u.right, u.node, u.index)
                u = u.next
            print("\t\t" + s)

    def set_growth_rate(self, growth_rate, time):
        # TODO This doesn't work because we need to know what the time
        # is so we can set the start size accordingly. Need to look at
        # ms's model carefully to see what it actually does here.
        new_size = self.get_size(time)
        self._start_size = new_size
        self._start_time = time
        self._growth_rate = growth_rate

    def set_start_size(self, start_size):
        self._start_size = start_size
        self._growth_rate = 0

    def get_num_ancestors(self):
        return len(self._ancestors)

    def get_size(self, t):
        """
        Returns the size of this population at time t.
        """
        dt = t - self._start_time
        return self._start_size * math.exp(-self._growth_rate * dt)

    def get_common_ancestor_waiting_time(self, t):
        """
        Returns the random waiting time until a common ancestor event
        occurs within this population.
        """
        ret = sys.float_info.max
        k = len(self._ancestors)
        if k > 1:
            u = random.expovariate(k * (k - 1))
            if self._growth_rate == 0:
                ret = self._start_size * u
            else:
                dt = t - self._start_time
                z = (1 + self._growth_rate * self._start_size *
                    math.exp(-self._growth_rate * dt) * u)
                if z > 0:
                    ret = math.log(z) / self._growth_rate
        return ret

    def remove(self, index):
        """
        Removes and returns the individual at the specified index.
        """
        return self._ancestors.pop(index)

    def add(self, individual):
        """
        Inserts the specified individual into this population.
        """
        self._ancestors.append(individual)

    def __iter__(self):
        return iter(self._ancestors)

class Simulator(object):
    """
    A reference implementation of the multi locus simulation algorithm.
    """
    def __init__(
            self, sample_size, num_loci, recombination_rate, migration_matrix,
            sample_configuration, population_growth_rates, population_sizes,
            population_growth_rate_changes, population_size_changes,
            migration_matrix_element_changes, bottlenecks, max_segments=100):
        # Must be a square matrix.
        N = len(migration_matrix)
        assert len(sample_configuration) == N
        assert len(population_growth_rates) == N
        assert len(population_sizes) == N
        for j in range(N):
            assert N == len(migration_matrix[j])
            assert migration_matrix[j][j] == 0
        assert sum(sample_configuration) == sample_size

        self.n = sample_size
        self.m = num_loci
        self.r = recombination_rate
        self.migration_matrix = migration_matrix
        self.max_segments = max_segments
        self.segment_stack = []
        self.segments = [None for j in range(self.max_segments + 1)]
        for j in range(self.max_segments):
            s = Segment(j + 1)
            self.segments[j + 1] = s
            self.segment_stack.append(s)
        self.P = [Population(id_) for id_ in range(N)]
        self.C = []
        self.L = FenwickTree(self.max_segments)
        self.S = bintrees.AVLTree()
        j = 0
        for pop_index in range(N):
            sample_size = sample_configuration[pop_index]
            self.P[pop_index].set_start_size(population_sizes[pop_index])
            self.P[pop_index].set_growth_rate(
                population_growth_rates[pop_index], 0)
            for k in range(sample_size):
                x = self.alloc_segment(0, self.m, j, pop_index)
                self.L.set_value(x.index, self.m - 1)
                self.P[pop_index].add(x)
                j += 1
        self.S[0] = self.n
        self.S[self.m] = -1
        self.t = 0
        self.w = self.n
        self.num_ca_events = 0
        self.num_re_events = 0
        self.modifier_events = [(sys.float_info.max, None, None)]
        for time, pop_id, new_size in population_size_changes:
            self.modifier_events.append(
                (time, self.change_population_size, (int(pop_id), new_size)))
        for time, pop_id, new_rate in population_growth_rate_changes:
            self.modifier_events.append(
                (time, self.change_population_growth_rate,
                    (int(pop_id), new_rate, time)))
        for time, pop_i, pop_j, new_rate in migration_matrix_element_changes:
            self.modifier_events.append(
                (time, self.change_migration_matrix_element,
                    (int(pop_i), int(pop_j), new_rate)))
        for time, pop_id, intensity in bottlenecks:
            self.modifier_events.append(
                (time, self.bottleneck_event, (int(pop_id), intensity)))
        self.modifier_events.sort()

    def change_population_size(self, pop_id, size):
        print("Changing pop size to ", size)
        self.P[pop_id].set_start_size(size)

    def change_population_growth_rate(self, pop_id, rate, time):
        print("Changing growth rate to ", rate)
        self.P[pop_id].set_growth_rate(rate, time)

    def change_migration_matrix_element(self, pop_i, pop_j, rate):
        print("Changing migration rate", pop_i, pop_j, rate)
        self.migration_matrix[pop_i][pop_j] = rate

    def alloc_segment(
            self, left, right, node, pop_index, prev=None, next=None):
        """
        Pops a new segment off the stack and sets its properties.
        """
        s = self.segment_stack.pop()
        s.left = left
        s.right = right
        s.node = node
        s.population = pop_index
        s.next = next
        s.prev = prev
        return s

    def free_segment(self, u):
        """
        Frees the specified segment making it ready for reuse and
        setting its weight to zero.
        """
        self.L.set_value(u.index, 0)
        self.segment_stack.append(u)

    def simulate(self):
        """
        Simulates the algorithm until all loci have coalesced.
        """
        infinity = sys.float_info.max
        while sum(pop.get_num_ancestors() for pop in self.P) != 0:
            self.verify()
            rate = self.r * self.L.get_total()
            t_re = infinity
            if rate != 0:
                t_re = random.expovariate(rate)
            # Common ancestor events occur within demes.
            t_ca = infinity
            for index, pop in enumerate(self.P):
                t = pop.get_common_ancestor_waiting_time(self.t)
                if t < t_ca:
                    t_ca = t
                    ca_population = index
            t_mig = infinity
            # Migration events happen at the rates in the matrix.
            for j in range(len(self.P)):
                source_size = self.P[j].get_num_ancestors()
                for k in range(len(self.P)):
                    rate = source_size * self.migration_matrix[j][k]
                    if rate > 0:
                        t = random.expovariate(rate)
                        if t < t_mig:
                            t_mig = t
                            mig_source = j
                            mig_dest = k
            min_time = min(t_re, t_ca, t_mig)
            assert min_time != infinity
            if self.t + min_time > self.modifier_events[0][0]:
                t, func, args = self.modifier_events.pop(0)
                self.t = t
                func(*args)
            else:
                self.t += min_time
                if min_time == t_re:
                    # print("RE EVENT")
                    self.recombination_event()
                elif min_time == t_ca:
                    # print("CA EVENT")
                    self.common_ancestor_event(ca_population)
                else:
                    # print("MIG EVENT")
                    self.migration_event(mig_source, mig_dest)

    def migration_event(self, j, k):
        """
        Migrates an individual from population j to population k.
        """
        # print("Migrating ind from ", j, " to ", k)
        # print("Population sizes:", [len(pop) for pop in self.P])
        index = random.randint(0, self.P[j].get_num_ancestors() - 1)
        x = self.P[j].remove(index)
        self.P[k].add(x)
        # Set the population id for each segment also.
        u = x
        while u is not None:
            u.population = k
            u = u.next
        # print("AFTER Population sizes:", [len(pop) for pop in self.P])

    def recombination_event(self):
        """
        Implements a recombination event.
        """
        self.num_re_events += 1
        h = random.randint(1, self.L.get_total())
        # Get the segment containing the h'th link
        y = self.segments[self.L.find(h)]
        k = y.right - self.L.get_cumulative_frequency(y.index) + h - 1
        x = y.prev
        if y.left < k:
            # Make new segment
            z = self.alloc_segment(
                k, y.right, y.node, y.population, None, y.next)
            if y.next is not None:
                y.next.prev = z
            y.next = None
            y.right = k
            self.L.increment(y.index, k - z.right)
        else:
            # split the link between x and y.
            x.next = None
            y.prev = None
            z = y
        self.L.set_value(z.index, z.right - z.left - 1)
        self.P[z.population].add(z)

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

    def bottleneck_event(self, pop_id, intensity):
        # self.print_state()
        # Merge some of the ancestors.
        pop = self.P[pop_id]
        H = []
        for _ in range(pop.get_num_ancestors()):
            if random.random() < intensity:
                x = pop.remove(0)
                heapq.heappush(H, (x.left, x))
        self.merge_ancestors(H, pop_id)

    def merge_ancestors(self, H, pop_id):
        pop = self.P[pop_id]
        defrag_required = False
        coalescence = False
        alpha = None
        z = None
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
                    pop.add(alpha)
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

    def defrag_segment_chain(self, z):
        y = z
        while y.prev is not None:
            x = y.prev
            if x.right == y.left and x.node == y.node:
                x.right = y.right
                x.next = y.next
                if y.next is not None:
                    y.next.prev = x
                self.L.increment(x.index, y.right - y.left)
                self.free_segment(y)
            y = x

    def defrag_breakpoints(self):
        # Defrag the breakpoints set
        j = 0
        k = 0
        while k < self.m:
            k = self.S.succ_key(j)
            if self.S[j] == self.S[k]:
                del self.S[k]
            else:
                j = k

    def common_ancestor_event(self, population_index):
        """
        Implements a coancestry event.
        """
        pop = self.P[population_index]
        self.num_ca_events += 1
        # Choose two ancestors uniformly.
        j = random.randint(0, pop.get_num_ancestors() - 1)
        x = pop.remove(j)
        j = random.randint(0, pop.get_num_ancestors() - 1)
        y = pop.remove(j)
        pop = self.P[population_index]
        z = None
        coalescence = False
        defrag_required = False
        while x is not None or y is not None:
            alpha = None
            if x is None or y is None:
                if x is not None:
                    alpha = x
                    x = None
                if y is not None:
                    alpha = y
                    y = None
            else:
                if y.left < x.left:
                    beta = x
                    x = y
                    y = beta
                if x.right <= y.left:
                    alpha = x
                    x = x.next
                    alpha.next = None
                elif x.left != y.left:
                    alpha = self.alloc_segment(
                        x.left, y.left, x.node, x.population)
                    x.left = y.left
                else:
                    if not coalescence:
                        coalescence = True
                        self.w += 1
                    u = self.w - 1
                    # Put in breakpoints for the outer edges of the coalesced
                    # segment
                    l = x.left
                    r_max = min(x.right, y.right)
                    if l not in self.S:
                        j = self.S.floor_key(l)
                        self.S[l] = self.S[j]
                    if r_max not in self.S:
                        j = self.S.floor_key(r_max)
                        self.S[r_max] = self.S[j]
                    # Update the number of extant segments.
                    if self.S[l] == 2:
                        self.S[l] = 0
                        r = self.S.succ_key(l)
                    else:
                        r = l
                        while r < r_max and self.S[r] != 2:
                            self.S[r] -= 1
                            r = self.S.succ_key(r)
                        alpha = self.alloc_segment(l, r, u, population_index)
                    self.C.append((l, r, u, (x.node, y.node), self.t))
                    # Now trim the ends of x and y to the right sizes.
                    if x.right == r:
                        self.free_segment(x)
                        x = x.next
                    else:
                        x.left = r
                    if y.right == r:
                        self.free_segment(y)
                        y = y.next
                    else:
                        y.left = r

            # loop tail; update alpha and integrate it into the state.
            if alpha is not None:
                if z is None:
                    pop.add(alpha)
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

    def print_state(self):
        print("State @ time ", self.t)
        print("Links = ", self.L.get_total())
        print("Modifier events = ")
        for t, f, args in self.modifier_events:
            print("\t", t, f, args)
        print("Population sizes:", [pop.get_num_ancestors() for pop in self.P])
        print("Migration Matrix:")
        for row in self.migration_matrix:
            print("\t", row)
        for population in self.P:
            population.print_state()
        print("Overlap counts", len(self.S))
        for k, x in self.S.items():
            print("\t", k, "\t:\t", x)
        print("Fenwick tree:", self.L.get_total())
        for j in range(1, self.max_segments + 1):
            s = self.L.get_frequency(j)
            if s != 0:
                print(
                    "\t", j, "->", s, self.L.get_cumulative_frequency(j))
        print("Coalescence records: ")
        for rec in self.C:
            print("\t", rec)
        self.verify()

    def verify(self):
        """
        Checks that the state of the simulator is consistent.
        """
        q = 0
        for pop_index, pop in enumerate(self.P):
            for u in pop:
                assert u.prev is None
                left = u.left
                right = u.left
                while u is not None:
                    assert u.population == pop_index
                    assert u.left <= u.right
                    if u.prev is not None:
                        s = u.right - u.prev.right
                    else:
                        s = u.right - u.left - 1
                    assert s == self.L.get_frequency(u.index)
                    right = u.right
                    v = u.next
                    if v is not None:
                        assert v.prev == u
                        if u.right > v.left:
                            print("ERROR", u, v)
                        assert u.right <= v.left
                    u = v
                q += right - left - 1
        assert q == self.L.get_total()

        assert self.S[self.m] == -1
        # Check the ancestry tracking.
        A = bintrees.AVLTree()
        A[0] = 0
        A[self.m] = -1
        for pop in self.P:
            for u in pop:
                while u is not None:
                    if u.left not in A:
                        k = A.floor_key(u.left)
                        A[u.left] = A[k]
                    if u.right not in A:
                        k = A.floor_key(u.right)
                        A[u.right] = A[k]
                    k = u.left
                    while k < u.right:
                        A[k] += 1
                        k = A.succ_key(k)
                    u = u.next
        # Now, defrag A
        j = 0
        k = 0
        while k < self.m:
            k = A.succ_key(j)
            if A[j] == A[k]:
                del A[k]
            else:
                j = k
        assert list(A.items()) == list(self.S.items())

    def verify_end(self):
        """
        Verify the state of the simulation at the end.
        """
        # Check the coalescence records to make sure they correctly cover
        # space
        left_coords = sorted(set(r[0] for r in self.C))
        for k in left_coords:
            c = 0
            for l, r, _, _, _, _ in self.C:
                if l <= k < r:
                    c += 1
            assert c == self.n - 1

    def write_records(self, out):
        """
        Writes the records out as text.
        """
        for left, right, u, children, time in self.C:
            print(
                left, right, u, ",".join(str(c) for c in sorted(children)),
                time, 0, sep="\t", file=out)


class Tree(object):
    """
    A tree in a tree sequence. This class represents a single point
    in the seqeunce of trees and has next() and prev() methods to
    change the state of the tree into other trees in the sequence.
    """
    def __init__(self, left, right, node, children, time):
        self.__left = left
        self.__right = right
        self.__node = node
        self.__children = children
        self.__time = time
        self.__num_records = len(left)
        self.__left_order = sorted(
            range(self.__num_records), key=lambda j: (left[j], time[j]))
        self.__right_order = sorted(
            range(self.__num_records), key=lambda j: (right[j], -time[j]))
        self.__direction = -1
        self.parent = [-1 for j in range(max(node) + 1)]
        self.num_trees = len(set(left))

    def first(self):
        self.__left_index = 0
        self.__right_index = 0
        self.__direction = 1
        self.index = -1
        self.next()

    def last(self):
        self.__left_index = self.__num_records - 1
        self.__right_index = self.__num_records - 1
        self.__direction = -1
        self.index = self.num_trees
        self.prev()

    def next(self):
        if self.index == self.num_trees - 1:
            raise ValueError("Cannot call next on last tree")
        j = int(self.__direction != 1)
        self.__direction = 1
        il = self.__left_index + j
        ir = self.__right_index + j
        x = self.__left[self.__left_order[il]]
        M = self.__num_records
        while self.__right[self.__right_order[ir]] == x:
            h = self.__right_order[ir]
            for q in self.__children[h]:
                self.parent[q] = -1
            ir += 1
        while il < M and self.__left[self.__left_order[il]] == x:
            h = self.__left_order[il]
            for q in self.__children[h]:
                self.parent[q] = self.__node[h]
            il += 1
        self.__left_index = il
        self.__right_index = ir
        self.index += 1

    def prev(self):
        if self.index == 0:
            raise ValueError("Cannot call prev on the first tree")
        j = int(self.__direction != -1)
        self.__direction = -1
        il = self.__left_index - j
        ir = self.__right_index - j
        x = self.__right[self.__right_order[ir]]
        while self.__left[self.__left_order[il]] == x:
            h = self.__left_order[il]
            for q in self.__children[h]:
                self.parent[q] = -1
            il -= 1
        while ir >= 0 and self.__right[self.__right_order[ir]] == x:
            h = self.__right_order[ir]
            for q in self.__children[h]:
                self.parent[q] = self.__node[h]
            ir -= 1
        self.__left_index = il
        self.__right_index = ir
        self.index -= 1


def generate_trees(l, r, u, c, t):
    """
    Algorithm T. Sequentially visits all trees in the specified
    tree sequence.
    """
    # Calculate the index vectors
    M = len(l)
    I = sorted(range(M), key=lambda j: (l[j], t[j]))
    O = sorted(range(M), key=lambda j: (r[j], -t[j]))
    pi = [-1 for j in range(max(u) + 1)]
    j = 0
    k = 0
    while j < M:
        x = l[I[j]]
        while r[O[k]] == x:
            h = O[k]
            for q in c[h]:
                pi[q] = -1
            k = k + 1
        while j < M and l[I[j]] == x:
            h = I[j]
            for q in c[h]:
                pi[q] = u[h]
            j += 1
        yield pi

def reverse_generate_trees(l, r, u, c, t):
    """
    Reversed version of Algorithm T. Sequentially visits all trees
    in the specified tree sequence in right-to-left order.
    """
    # Calculate the index vectors
    M = len(l)
    O = sorted(range(M), key=lambda j: (l[j], t[j]))
    I = sorted(range(M), key=lambda j: (r[j], -t[j]))
    pi = [-1 for j in range(max(u) + 1)]
    j = M - 1
    k = M - 1
    while j >= 0:
        x = r[I[j]]
        while l[O[k]] == x:
            h = O[k]
            for q in c[h]:
                pi[q] = -1
            k -= 1
        while j >= 0 and r[I[j]] == x:
            h = I[j]
            for q in c[h]:
                pi[q] = u[h]
            j -= 1
        yield pi

def count_leaves(l, r, u, c, t, S):
    """
    Algorithm L. Sequentially visits all trees in the specified
    tree sequence and maintain a count of the leaf nodes in the
    specified set for each node.
    """
    # Calculate the index vectors
    M = len(l)
    I = sorted(range(M), key=lambda j: (l[j], t[j]))
    O = sorted(range(M), key=lambda j: (r[j], -t[j]))
    pi = [-1 for j in range(max(u) + 1)]
    beta = [0 for j in range(max(u) + 1)]
    for j in S:
        beta[j] = 1
    j = 0
    k = 0
    while j < M:
        x = l[I[j]]
        while r[O[k]] == x:
            h = O[k]
            b = 0
            for q in c[h]:
                pi[q] = -1
                b += beta[q]
            k += 1
            v = u[h]
            while v != -1:
                beta[v] -= b
                v = pi[v]
        while j < M and l[I[j]] == x:
            h = I[j]
            b = 0
            for q in c[h]:
                pi[q] = u[h]
                b += beta[q]
            j = j + 1
            v = u[h]
            while v != -1:
                beta[v] += b
                v = pi[v]
        yield pi, beta

class LeafListNode(object):
    def __init__(self, value, next=None):
        self.value = value
        self.next = next

    def __str__(self):
        next = -1 if self.next is None else self.next.value
        return "{}->{}".format(self.value, next)

def propagate_leaf_loss(u, pi, xi, head, tail):
    # Invalidate the head and tail pointers above u that depend
    # on this node.
    head[u] = None
    tail[u] = None


def propagate_leaf_gain(u, pi, xi, head, tail):
    num_children = len(xi[u])
    for j in range(1, num_children):
        tail[xi[u][j - 1]].next = head[xi[u][j]]
    head[u] = head[xi[u][0]]
    tail[u] = tail[xi[u][-1]]

    v = u
    w = pi[v]
    while w != -1:
        j = xi[w].index(v)
        if j != 0:
            break
        head[w] = head[u]
        v = w
        w = pi[w]
    v = u
    w = pi[v]
    while w != -1:
        j = xi[w].index(v)
        if j !=len(xi[w]) - 1:
            break
        tail[w] = tail[u]
        v = w
        w = pi[w]

def post_propagate_leaf_gain(u, pi, xi, head, tail):
    v = u
    w = pi[v]
    while w != -1:
        j = xi[w].index(v)
        if j < len(xi[w]) - 1:
            tail[v].next = head[xi[w][j + 1]]
        if j > 0:
            tail[xi[w][j - 1]].next = head[v]
        v = w
        w = pi[w]


def leaf_sets(l, r, u, c, t, S):
    """
    Sequentially visits all trees in the specified
    tree sequence and maintain the leaf sets for all leaves in
    specified set for each node.
    """
    # Calculate the index vectors
    M = len(l)
    I = sorted(range(M), key=lambda j: (l[j], t[j]))
    O = sorted(range(M), key=lambda j: (r[j], -t[j]))
    pi = [-1 for j in range(max(u) + 1)]
    xi = [[] for j in range(max(u) + 1)]
    head = [None for j in range(max(u) + 1)]
    tail = [None for j in range(max(u) + 1)]
    for j in S:
        node = LeafListNode(j)
        head[j] = node
        tail[j] = node
    j = 0
    k = 0
    while j < M:
        x = l[I[j]]
        while r[O[k]] == x:
            h = O[k]
            propagate_leaf_loss(u[h], pi, xi, head, tail)
            for q in c[h]:
                pi[q] = -1
            xi[u[h]] = []
            k += 1
        before = j
        while j < M and l[I[j]] == x:
            h = I[j]
            for q in c[h]:
                pi[q] = u[h]
            xi[u[h]] = c[h]
            propagate_leaf_gain(u[h], pi, xi, head, tail)
            j += 1
        j = before
        while j < M and l[I[j]] == x:
            h = I[j]
            post_propagate_leaf_gain(u[h], pi, xi, head, tail)
            j += 1
        yield pi, xi, head, tail


def check_consistency(n, pi, xi, head, tail):
    """
    Checks the consistency of the specified parent list, child list
    and head and tail leaf list pointers.
    """
    root = 0
    while pi[root] != -1:
        root = pi[root]
    all_leaves = list(leaves(root, xi))
    assert set(all_leaves) == set(range(n))
    for u in nodes(root, xi):
        node_leaves = list(leaves(u, xi))
        if node_leaves[0] != head[u].value:
            print("HERROR: head incorrect:", head[u].value)
        if node_leaves[-1] != tail[u].value:
            print("TERROR: tail incorrect:", tail[u].value)
        list_leaves = []
        x = head[u]
        while True:
            if x.value in list_leaves:
                print("ERROR!!!", x.value, "already in leaf list at index",
                        list_leaves.index(x.value), "len = ", len(list_leaves))
                break
            list_leaves.append(x.value)
            if x == tail[u]:
                break
            x = x.next
        if list_leaves != node_leaves:
            print("ERROR")
            print(list_leaves)
            print(node_leaves)
        # assert list_leaves == node_leaves
        # print(list_leaves)
        # print(list_leaves == node_leaves)

    print("TREE")
    print_tree(root, xi, head, tail, 0)

def print_tree(u, xi, head, tail, depth):
    indent = "  " * depth
    leaf_list = []
    # x = head[u]
    # while True:
    #     leaf_list.append(x.value)
    #     if x == tail[u]:
    #         break
    #     x = x.next
    print("{}[{}:{}]\thead = {}; tail = {}\t{}".format(
        indent, u, len(xi[u]), head[u], tail[u], leaf_list))
    for v in xi[u]:
        print_tree(v, xi, head, tail, depth + 1)

def nodes(root, xi):
    """
    Returns an iterator over the nodes in the specified tree
    """
    stack = [root]
    while len(stack) > 0:
        u = stack.pop()
        stack.extend(reversed(xi[u]))
        yield u

def leaves(u, xi):
    """
    Returns an iterator over the leaves below the specified node.
    """
    for v in nodes(u, xi):
        if len(xi[v]) == 0:
            yield v


def run_trees(args):
    process_trees(args.history_file)


def process_trees(records_file):
    # Read in the records
    l = []
    r = []
    u = []
    c = []
    t = []
    with open(records_file) as f:
        for line in f:
            toks = line.split()
            l.append(float(toks[0]))
            r.append(float(toks[1]))
            u.append(int(toks[2]))
            children = list(map(int, toks[3].split(",")))
            c.append(children)
            t.append(float(toks[4]))
    # N = len(l)
    print("Trees:")
    forward = []
    for pi in generate_trees(l, r, u, c, t):
        forward.append(list(pi))
    T = len(forward)
    print("START")
    tree = Tree(l, r, u, c, t)
    tree.first()
    for j in range(T):
        assert tree.index == j
        assert tree.parent == forward[j]
        if j < T - 1:
            tree.next()

    tree = Tree(l, r, u, c, t)
    tree.last()
    for j in range(T):
        print(j, tree.index, T)
        assert tree.index == T - j - 1
        assert tree.parent == forward[T - j - 1]
        if j < T - 1:
            tree.prev()

    tree = Tree(l, r, u, c, t)
    tree.first()
    print("first done")
    assert tree.parent == forward[0]
    assert tree.index == 0
    for _ in range(1):
        for _ in range(T // 2):
            tree.next()
            print("\t", tree.index)
            assert tree.parent == forward[tree.index]
        for j in range(10):
            print("Reverse")
            for _ in range(T // 3):
                tree.prev()
                print("\t", tree.index)
                assert tree.parent == forward[tree.index]
            print("Forward")
            for _ in range(T // 3):
                tree.next()
                print("\t", tree.index)
                assert tree.parent == forward[tree.index]


    # n = min(u)
    # S = set(range(n))
    # # print("Counts:")
    # # for pi, beta in count_leaves(l, r, u, c, t, S):
    # #     print("\t", beta)
    # print("Counts:")
    # for pi, xi, head, tail in leaf_sets(l, r, u, c, t, S):
    #     check_consistency(n, pi, xi, head, tail)
    #     print(pi)
    #     print(xi)

def run_simulate(args):
    """
    Runs the simulation and outputs the results in text.
    """
    n = args.sample_size
    m = args.num_loci
    rho = args.recombination_rate
    num_populations = args.num_populations
    migration_matrix = [
        [args.migration_rate * int(j != k) for j in range(num_populations)]
        for k in range(num_populations)]
    sample_configuration = [0 for j in range(num_populations)]
    population_growth_rates = [0 for j in range(num_populations)]
    population_sizes = [1 for j in range(num_populations)]
    sample_configuration[0] = n
    if args.sample_configuration is not None:
        sample_configuration = args.sample_configuration
    if args.population_growth_rates is not None:
        population_growth_rates = args.population_growth_rates
    if args.population_sizes is not None:
        population_sizes = args.population_sizes
    random.seed(args.random_seed)
    s = Simulator(
        n, m, rho, migration_matrix,
        sample_configuration, population_growth_rates,
        population_sizes, args.population_growth_rate_change,
        args.population_size_change,
        args.migration_matrix_element_change,
        args.bottleneck, 10000)
    s.simulate()
    # TEMP
    with tempfile.NamedTemporaryFile(prefix="msp_alg") as f:
        s.write_records(f)
        f.flush()
        process_trees(f.name)

def add_simulator_arguments(parser):
    parser.add_argument("sample_size", type=int)
    parser.add_argument(
        "--random-seed", "-s", type=int, default=1)
    parser.add_argument(
        "--num-loci", "-m", type=int, default=100)
    parser.add_argument(
        "--num-replicates", "-R", type=int, default=1000)
    parser.add_argument(
        "--recombination-rate", "-r", type=float, default=0.1)
    parser.add_argument(
        "--num-populations", "-p", type=int, default=1)
    parser.add_argument(
        "--migration-rate", "-g", type=float, default=1)
    parser.add_argument(
        "--sample-configuration", type=int, nargs="+", default=None)
    parser.add_argument(
        "--population-growth-rates", type=float, nargs="+", default=None)
    parser.add_argument(
        "--population-sizes", type=float, nargs="+", default=None)
    parser.add_argument(
        "--population-size-change", type=float, nargs=3, action="append",
        default=[])
    parser.add_argument(
        "--population-growth-rate-change", type=float, nargs=3,
        action="append", default=[])
    parser.add_argument(
        "--migration-matrix-element-change", type=float, nargs=4,
        action="append", default=[])
    parser.add_argument(
        "--bottleneck", type=float, nargs=3, action="append", default=[])


def main():
    parser = argparse.ArgumentParser()
    # This is required to get uniform behaviour in Python2 and Python3
    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    simulate_parser = subparsers.add_parser(
        "simulate",
        help="Simulate the process and output the results in text")
    add_simulator_arguments(simulate_parser)
    simulate_parser.set_defaults(runner=run_simulate)

    trees_parser = subparsers.add_parser(
        "trees",
        help="Shows the trees from an text records file")
    trees_parser.add_argument("history_file")

    trees_parser.set_defaults(runner=run_trees)

    args = parser.parse_args()
    args.runner(args)


if __name__ == "__main__":
    main()
