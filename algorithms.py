"""
Python version of the simulation algorithm.
"""
from __future__ import print_function
from __future__ import division

import sys
import random
import argparse
import heapq
import math
import numpy as np

import bintrees
import msprime


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

    def __lt__(self, other):
        return ((self.left, self.right, self.population, self.node)
                < (other.left, other.right, other.population, self.node))


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
                z = (
                    1 + self._growth_rate * self._start_size
                    * math.exp(-self._growth_rate * dt) * u)
                if z > 0:
                    ret = math.log(z) / self._growth_rate
        return ret

    def get_ind_range(self, t):
        """ Returns ind labels at time t """
        first_ind = np.sum([self.get_size(t_prev) for t_prev in range(0, t)])
        last_ind = first_ind + self.get_size(t)

        return range(int(first_ind), int(last_ind)+1)

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
            migration_matrix_element_changes, bottlenecks, model='hudson',
	    max_segments=100):
        # Must be a square matrix.
        N = len(migration_matrix)
        assert len(sample_configuration) == N
        assert len(population_growth_rates) == N
        assert len(population_sizes) == N
        for j in range(N):
            assert N == len(migration_matrix[j])
            assert migration_matrix[j][j] == 0
        assert sum(sample_configuration) == sample_size

        self.model = model
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
        self.L = FenwickTree(self.max_segments)
        self.S = bintrees.AVLTree()
        # The output tree sequence.
        self.nodes = msprime.NodeTable()
        self.edges = msprime.EdgeTable()
        self.edge_buffer = []
        for pop_index in range(N):
            sample_size = sample_configuration[pop_index]
            self.P[pop_index].set_start_size(population_sizes[pop_index])
            self.P[pop_index].set_growth_rate(
                population_growth_rates[pop_index], 0)
            for k in range(sample_size):
                j = len(self.nodes)
                x = self.alloc_segment(0, self.m, j, pop_index)
                self.L.set_value(x.index, self.m - 1)
                self.P[pop_index].add(x)
                self.nodes.add_row(
                    flags=msprime.NODE_IS_SAMPLE, time=0, population=pop_index)
                j += 1
        self.S[0] = self.n
        self.S[self.m] = -1
        self.t = 0
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

    def store_node(self, population):
        self.flush_edges()
        self.nodes.add_row(time=self.t, population=population)

    def flush_edges(self):
        """
        Flushes the edges in the edge buffer to the table, squashing any adjacent edges.
        """
        if len(self.edge_buffer) > 0:
            parent = len(self.nodes) - 1
            self.edge_buffer.sort(key=lambda e: (e.child, e.left))
            left = self.edge_buffer[0].left
            right = self.edge_buffer[0].right
            child = self.edge_buffer[0].child
            assert self.edge_buffer[0].parent == parent
            for e in self.edge_buffer[1:]:
                assert e.parent == parent
                if e.left != right or e.child != child:
                    self.edges.add_row(left, right, parent, child)
                    left = e.left
                    child = e.child
                right = e.right
            self.edges.add_row(left, right, parent, child)
            self.edge_buffer = []

    def store_edge(self, left, right, parent, child):
        """
        Stores the specified edge to the output tree sequence.
        """
        self.edge_buffer.append(
            msprime.Edge(left=left, right=right, parent=parent, child=child))

    def finalise(self):
        """
        Finalises the simulation returns an msprime tree sequence object.
        """
        self.flush_edges()
        return msprime.load_tables(nodes=self.nodes, edges=self.edges)

    def simulate(self, model='hudson'):
        if self.model == 'hudson':
            self.hudson_simulate()
        elif self.model == 'dtwf':
            self.dtwf_simulate()
        else:
            print("Error: bad model specification -", self.model)
            raise ValueError

        return self.finalise()

    def hudson_simulate(self):
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
                    self.hudson_recombination_event()
                elif min_time == t_ca:
                    # print("CA EVENT")
                    self.common_ancestor_event(ca_population)
                else:
                    # print("MIG EVENT")
                    self.migration_event(mig_source, mig_dest)
        return self.finalise()

    def dtwf_simulate(self):
        """
        Simulates the algorithm until all loci have coalesced.
        """
        while sum(pop.get_num_ancestors() for pop in self.P) != 0:
            self.t += 1
            self.verify()

            self.dtwf_generation()

    def dtwf_generation(self):
        """
        Evolves one generation of a Wright Fisher population
        """
        for pop_idx, pop in enumerate(self.P):
            ## Cluster haploid inds by parent
            cur_inds = pop.get_ind_range(self.t)
            offspring = bintrees.AVLTree()
            for i in range(pop.get_num_ancestors()-1, -1, -1):
                ## Popping every ancestor every generation is inefficient.
                ## In the C implementation we store a pointer to the 
                ## ancestor so we can pop only if we need to merge
                anc = pop.remove(i)
                parent = np.random.choice(cur_inds)
                if parent not in offspring:
                    offspring[parent] = []
                offspring[parent].append(anc)

            ## Draw recombinations in children and sort segments by
            ## inheritance direction
            for children in offspring.values():
                H = [[], []]
                for child in children:
                    segs_pair = self.dtwf_recombine(child)

                    ## Collect segments inherited from the same individual
                    for i, seg in enumerate(segs_pair):
                        if seg is None:
                            continue
                        assert seg.prev is None
                        heapq.heappush(H[i], (seg.left, seg))

                ## Merge segments
                for h in H:
                    self.merge_ancestors(h, pop_idx)

        ## Migration events happen at the rates in the matrix.
        for j in range(len(self.P)):
            source_size = self.P[j].get_num_ancestors()
            for k in range(len(self.P)):
                if j == k:
                    continue
                mig_rate = source_size * self.migration_matrix[j][k]
                num_migs = min(source_size, np.random.poisson(mig_rate))
                for _ in range(num_migs):
                    mig_source = j
                    mig_dest = k
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

    def hudson_recombination_event(self):
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

    def dtwf_recombine(self, x):
        """
        Chooses breakpoints and returns segments sorted by inheritance
        direction, by iterating through segment chain starting with x
        """
        u = self.alloc_segment(-1, -1, -1, -1, None, None)
        v = self.alloc_segment(-1, -1, -1, -1, None, None)
        seg_tails = [u, v]

        if self.r > 0:
            mu = 1. / self.r
            k = 1. + x.left + np.random.exponential(mu)
        else:
            mu = np.inf
            k = np.inf

        ix = np.random.randint(2)
        seg_tails[ix].next = x
        seg_tails[ix] = x

        while x is not None:
            seg_tails[ix] = x
            y = x.next

            if x.right > k:
                assert x.left <= k
                self.num_re_events += 1
                ix = (ix + 1) % 2
                # Make new segment
                z = self.alloc_segment(
                    k, x.right, x.node, x.population, seg_tails[ix], x.next)
                if x.next is not None:
                    x.next.prev = z
                seg_tails[ix].next = z
                seg_tails[ix] = z
                x.next = None
                x.right = k
                x = z
                k = 1 + k + np.random.exponential(mu)
            elif x.right <= k and y is not None and y.left >= k:
                ## Recombine between segment and the next
                assert seg_tails[ix] == x
                x.next = None
                y.prev = None
                while y.left > k:
                    self.num_re_events += 1
                    ix = (ix + 1) % 2
                    k = 1 + k + np.random.exponential(1. / self.r)
                seg_tails[ix].next = y
                y.prev = seg_tails[ix]
                seg_tails[ix] = y
                x = y
            else:
                ## No recombination between x.right and y.left
                x = y

        ## Remove sentinal segments - this can be handled more simply
        ## with pointers in C implemetation
        if u.next is not None:
            u.next.prev = None
        s = u
        u = s.next
        self.free_segment(s)

        if v.next is not None:
            v.next.prev = None
        s = v
        v = s.next
        self.free_segment(s)

        return u, v

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
                    self.store_node(pop_id)
                u = len(self.nodes) - 1
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
                for x in X:
                    self.store_edge(l, r, u, x.node)
                    if x.right == r:
                        self.free_segment(x)
                        if x.next is not None:
                            y = x.next
                            heapq.heappush(H, (y.left, y))
                    elif x.right > r:
                        x.left = r
                        heapq.heappush(H, (x.left, x))

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
                        self.store_node(population_index)
                    u = len(self.nodes) - 1
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
                    self.store_edge(l, r, u, x.node)
                    self.store_edge(l, r, u, y.node)
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
        print("nodes")
        print(self.nodes)
        print("edges")
        print(self.edges)
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
                    if self.model != 'dtwf':
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
        if self.model != 'dtwf':
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
        args.bottleneck, args.model, 10000)
    ts = s.simulate()
    ts.dump(args.output_file)
    if args.verbose:
        s.print_state()


def add_simulator_arguments(parser):
    parser.add_argument("sample_size", type=int)
    parser.add_argument("output_file")
    parser.add_argument(
            "-v", "--verbose", help="increase output verbosity", action="store_true")
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
    parser.add_argument(
        "--model", default='hudson')


def main():
    parser = argparse.ArgumentParser()
    add_simulator_arguments(parser)
    args = parser.parse_args()
    run_simulate(args)


if __name__ == "__main__":
    main()
