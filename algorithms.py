"""
Python version of the simulation algorithm.
"""
import argparse
import heapq
import logging
import math
import random
import sys

import bintrees
import daiquiri
import numpy as np
import tskit

import msprime


logger = daiquiri.getLogger()


class FenwickTree:
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
        self.__value = [0 for j in range(max_index + 1)]
        # Compute the binary logarithm of max_index
        u = self.__max_index
        while u != 0:
            self.__log_max_index = u
            u -= u & -u

    def get_total(self):
        """
        Returns the total cumulative frequency over all indexes.
        """
        return self.get_cumulative_sum(self.__max_index)

    def increment(self, index, v):
        """
        Increments the frequency of the specified index by the specified
        value.
        """
        assert 0 < index <= self.__max_index
        self.__value[index] += v
        j = index
        while j <= self.__max_index:
            self.__tree[j] += v
            j += j & -j

    def set_value(self, index, v):
        """
        Sets the frequency at the specified index to the specified value.
        """
        f = self.get_value(index)
        self.increment(index, v - f)

    def get_cumulative_sum(self, index):
        """
        Returns the cumulative frequency of the specified index.
        """
        assert 0 < index <= self.__max_index
        j = index
        s = 0
        while j > 0:
            s += self.__tree[j]
            j -= j & -j
        return s

    def get_value(self, index):
        """
        Returns the frequency of the specified index.
        """
        return self.__value[index]

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


class Segment:
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
        self.label = 0
        self.index = index

    @staticmethod
    def show_chain(seg):
        s = ""
        while seg is not None:
            s += f"[{seg.left}, {seg.right}: {seg.node}], "
            seg = seg.next
        return s[:-2]

    def __lt__(self, other):
        return (self.left, self.right, self.population, self.node) < (
            other.left,
            other.right,
            other.population,
            self.node,
        )


class Population:
    """
    Class representing a population in the simulation.
    """

    def __init__(self, id_, num_labels=1):
        self.id = id_
        self.start_time = 0
        self.start_size = 1.0
        self.growth_rate = 0
        # Keep a list of each label.
        # We'd like to use AVLTrees here for P but the API doesn't quite
        # do what we need. Lists are inefficient here and should not be
        # used in a real implementation.
        self._ancestors = [[] for _ in range(num_labels)]

    def print_state(self):
        print("Population ", self.id)
        print("\tstart_size = ", self.start_size)
        print("\tgrowth_rate = ", self.growth_rate)
        print("\tAncestors: ", len(self._ancestors))
        for label, ancestors in enumerate(self._ancestors):
            print("\tLabel = ", label)
            for u in ancestors:
                print("\t\t" + Segment.show_chain(u))

    def set_growth_rate(self, growth_rate, time):
        # TODO This doesn't work because we need to know what the time
        # is so we can set the start size accordingly. Need to look at
        # ms's model carefully to see what it actually does here.
        new_size = self.get_size(time)
        self.start_size = new_size
        self.start_time = time
        self.growth_rate = growth_rate

    def set_start_size(self, start_size):
        self.start_size = start_size
        self.growth_rate = 0

    def get_num_ancestors(self, label=None):
        if label is None:
            return sum(len(label_ancestors) for label_ancestors in self._ancestors)
        else:
            return len(self._ancestors[label])

    def get_size(self, t):
        """
        Returns the size of this population at time t.
        """
        dt = t - self.start_time
        return self.start_size * math.exp(-self.growth_rate * dt)

    def get_common_ancestor_waiting_time(self, t):
        """
        Returns the random waiting time until a common ancestor event
        occurs within this population.
        """
        ret = sys.float_info.max
        k = self.get_num_ancestors()
        if k > 1:
            u = random.expovariate(k * (k - 1))
            if self.growth_rate == 0:
                ret = self.start_size * u
            else:
                dt = t - self.start_time
                z = (
                    1
                    + self.growth_rate
                    * self.start_size
                    * math.exp(-self.growth_rate * dt)
                    * u
                )
                if z > 0:
                    ret = math.log(z) / self.growth_rate
        return ret

    def get_ind_range(self, t):
        """ Returns ind labels at time t """
        first_ind = np.sum([self.get_size(t_prev) for t_prev in range(0, int(t))])
        last_ind = first_ind + self.get_size(t)

        return range(int(first_ind), int(last_ind) + 1)

    def remove(self, index, label=0):
        """
        Removes and returns the individual at the specified index.
        """
        return self._ancestors[label].pop(index)

    def remove_individual(self, individual, label=0):
        """
        Removes the given individual from its population.
        """
        return self._ancestors[label].remove(individual)

    def add(self, individual, label=0):
        """
        Inserts the specified individual into this population.
        """
        assert individual.label == label
        self._ancestors[label].append(individual)

    def __iter__(self):
        # will default to label 0
        # inter_label() extends behavior
        return iter(self._ancestors[0])

    def iter_label(self, label):
        """
        Iterates ancestors in popn from a label
        """
        return iter(self._ancestors[label])

    def iter_ancestors(self):
        """
        Iterates over all ancestors in a population over all labels.
        """
        for ancestors in self._ancestors:
            yield from ancestors

    def find_indv(self, indv):
        """
        find the index of an ancestor in population
        """
        return self._ancestors[indv.label].index(indv)


# TODO this needs to be refactored a bit and cleaned up. We can integrate
# better with tskit for this.
class Pedigree:
    """
    Class representing a pedigree for use with the DTWF model, as implemented
    in C library
    """

    def __init__(self, num_individuals, ploidy):
        self.ploidy = ploidy
        self.num_individuals = num_individuals
        self.inds = []
        self.samples = []
        self.num_samples = 0
        self.ind_heap = []
        self.is_climbing = False

        # Stores most recently merged segment
        self.merged_segment = None

    def set_pedigree(self, inds, parents, times, is_sample):
        self.ploidy = parents.shape[1]
        self.is_sample = is_sample
        self.inds = [
            Individual(ploidy=self.ploidy) for i in range(self.num_individuals)
        ]

        num_samples = 0
        for i in range(self.num_individuals):
            ind = self.inds[i]
            ind.id = inds[i]
            assert ind.id > 0
            ind.time = times[i]
            for j, parent in enumerate(parents[i]):
                if parent != tskit.NULL:
                    assert parent >= 0
                    ind.parents[j] = self.inds[parent]

            if is_sample[i] != 0:
                assert is_sample[i] == 1
                self.samples.append(ind)
                num_samples += 1

        self.num_samples = num_samples

    def load_pop(self, pop):
        """
        Loads segments from a given pop into the pedigree samples
        """
        if self.num_sample_lineages() != pop.get_num_ancestors():
            err_str = (
                "Ped samples: "
                + str(self.num_sample_lineages())
                + " Samples: "
                + str(pop.get_num_ancestors())
                + " - must be equal!"
            )
            raise ValueError(err_str)

        # for i, anc in enumerate(pop):
        for i in range(pop.get_num_ancestors() - 1, -1, -1):
            anc = pop.remove(i)
            # Each individual gets 'ploidy' lineages
            ind = self.samples[i // self.ploidy]
            parent_ix = i % self.ploidy
            ind.add_segment(anc, parent_ix=parent_ix)

        # Add samples to queue to prepare for climbing - might be better if
        # included in previous loop
        for ind in self.samples:
            self.push_ind(ind)

    def assign_times(self, check=False):
        """
        For pedigrees without specified times, crudely assigns times to
        all individuals.
        """
        if len(self.samples) == 0:
            self.set_samples()
        assert len(self.samples) > 0

        climbers = [s for s in self.samples]
        t = 0
        while len(climbers) > 0:
            next_climbers = []
            for climber in climbers:
                if climber.time < t:
                    climber.time = t
                for parent in climber.parents:
                    if parent is not None:
                        next_climbers.append(parent)
            climbers = next_climbers
            t += 1

        if check:
            for ind in self.inds:
                for parent in ind.parents:
                    if parent is not None:
                        assert ind.time < parent.time

    def build_ind_queue(self):
        """
        Set up heap queue of samples, so most recent can be popped for merge.
        Heapify in case samples are not all at t=0.
        """
        self.ind_heap = [(ind.time, ind) for ind in self.samples]
        heapq.heapify(self.ind_heap)

    def push_ind(self, ind):
        """
        Adds an individual to the heap queue
        """
        assert ind.queued is False
        ind.queued = True
        heapq.heappush(self.ind_heap, ind)

    def pop_ind(self):
        """
        Pops the most recent individual off the heap queue
        """
        ind = heapq.heappop(self.ind_heap)
        assert ind.queued
        ind.queued = False

        return ind

    def num_sample_lineages(self):
        return len(self.samples) * self.ploidy

    def print_samples(self):
        for s in self.samples:
            print(s)


class Individual:
    """
    Class representing a diploid individual in the DTWF model. Trying to make
    arbitrary ploidy possible at some point in the future.
    """

    def __init__(self, ploidy=2):
        self.id = None  # This is the index of the individual in pedigree.inds
        self.ploidy = ploidy
        self.parents = [None for i in range(ploidy)]
        self.segments = [[] for i in range(ploidy)]
        self.sex = None
        self.time = -1
        self.queued = False

        # For debugging - to ensure we only merge once
        self.merged = False

    def __str__(self):
        parents = []
        for p in self.parents:
            if p is not None:
                parents.append(str(p.id))
            else:
                parents.append("None")

        parents_str = ",".join(parents)

        return f"(ID: {self.id}, time: {self.time}, parents: {parents_str})"

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return self.time < other.time

    def add_segment(self, seg, parent_ix):
        """
        Adds a segment to a parental segment heap, which allows easy merging
        later.
        """
        heapq.heappush(self.segments[parent_ix], (seg.left, seg))

    def num_lineages(self):
        return sum([len(s) for s in self.segments])


class TrajectorySimulator:
    """
    Class to simulate an allele frequency trajectory on which to condition
    the coalescent simulation.
    """

    def __init__(self, initial_freq, end_freq, alpha, time_slice):
        self._initial_freq = initial_freq
        self._end_freq = end_freq
        self._alpha = alpha
        self._time_slice = time_slice
        self._reset()

    def _reset(self):
        self._allele_freqs = []
        self._times = []

    def _genic_selection_stochastic_forwards(self, dt, freq, alpha):
        ux = (alpha * freq * (1 - freq)) / np.tanh(alpha * freq)
        sign = 1 if random.random() < 0.5 else -1
        freq += (ux * dt) + sign * np.sqrt(freq * (1.0 - freq) * dt)
        return freq

    def _simulate(self):
        """
        Proposes a sweep trajectory and returns the acceptance probability.
        """
        x = self._end_freq  # backward time
        current_size = 1
        t_inc = self._time_slice
        t = 0
        while x > self._initial_freq:
            # print("x: ",x)
            self._allele_freqs.append(max(x, self._initial_freq))
            self._times.append(t)
            # just a note below
            # current_size = self._size_calculator(t)
            #
            x = 1.0 - self._genic_selection_stochastic_forwards(
                t_inc, 1.0 - x, self._alpha * current_size
            )
            t += self._time_slice
        # will want to return current_size / N_max
        # for prototype this always equals 1
        return 1

    def run(self):
        while random.random() > self._simulate():
            self.reset()
        return self._allele_freqs, self._times


class RateMap:
    def __init__(self, positions, rates):
        self.positions = positions
        self.rates = rates
        self.cumulative = RateMap.recomb_mass(positions, rates)

    @staticmethod
    def recomb_mass(positions, rates):
        recomb_mass = 0
        cumulative = [recomb_mass]
        for i in range(1, len(positions)):
            recomb_mass += (positions[i] - positions[i - 1]) * rates[i - 1]
            cumulative.append(recomb_mass)
        return cumulative

    @property
    def sequence_length(self):
        return self.positions[-1]

    @property
    def total_mass(self):
        return self.cumulative[-1]

    @property
    def mean_rate(self):
        return self.total_mass / self.sequence_length

    def mass_between(self, left, right):
        left_mass = self.position_to_mass(left)
        right_mass = self.position_to_mass(right)
        return right_mass - left_mass

    def position_to_mass(self, pos):
        if pos == self.positions[0]:
            return 0
        if pos >= self.positions[-1]:
            return self.cumulative[-1]

        index = self._search(self.positions, pos)
        assert index > 0
        index -= 1
        offset = pos - self.positions[index]
        return self.cumulative[index] + offset * self.rates[index]

    def mass_to_position(self, recomb_mass):
        if recomb_mass == 0:
            return 0
        index = self._search(self.cumulative, recomb_mass)
        assert index > 0
        index -= 1
        mass_in_interval = recomb_mass - self.cumulative[index]
        pos = self.positions[index] + (mass_in_interval / self.rates[index])
        return pos

    def shift_by_mass(self, pos, mass):
        result_mass = self.position_to_mass(pos) + mass
        return self.mass_to_position(result_mass)

    def _search(self, values, query):
        left = 0
        right = len(values) - 1
        while left < right:
            m = (left + right) // 2
            if values[m] < query:
                left = m + 1
            else:
                right = m
        return left


class OverlapCounter:
    def __init__(self, seq_length):
        self.seq_length = seq_length
        self.overlaps = self._make_segment(0, seq_length, 0)

    def overlaps_at(self, pos):
        assert 0 <= pos < self.seq_length
        curr_interval = self.overlaps
        while curr_interval is not None:
            if curr_interval.left <= pos < curr_interval.right:
                return curr_interval.node
            curr_interval = curr_interval.next
        raise ValueError("Bad overlap count chain")

    def increment_interval(self, left, right):
        """
        Increment the count that spans the interval
        [left, right), creating additional intervals in overlaps
        if necessary.
        """
        curr_interval = self.overlaps
        while left < right:
            if curr_interval.left == left:
                if curr_interval.right <= right:
                    curr_interval.node += 1
                    left = curr_interval.right
                    curr_interval = curr_interval.next
                else:
                    self._split(curr_interval, right)
                    curr_interval.node += 1
                    break
            else:
                if curr_interval.right < left:
                    curr_interval = curr_interval.next
                else:
                    self._split(curr_interval, left)
                    curr_interval = curr_interval.next

    def _split(self, seg, bp):  # noqa: A002
        """
        Split the segment at breakpoint and add in another segment
        from breakpoint to seg.right. Set the original segment's
        right endpoint to breakpoint
        """
        right = self._make_segment(bp, seg.right, seg.node)
        if seg.next is not None:
            seg.next.prev = right
            right.next = seg.next
        right.prev = seg
        seg.next = right
        seg.right = bp

    def _make_segment(self, left, right, count):
        seg = Segment(0)
        seg.left = left
        seg.right = right
        seg.node = count
        return seg


class Simulator:
    """
    A reference implementation of the multi locus simulation algorithm.
    """

    def __init__(
        self,
        *,
        tables,
        recombination_map,
        migration_matrix,
        population_growth_rates,
        population_sizes,
        population_growth_rate_changes,
        population_size_changes,
        migration_matrix_element_changes,
        bottlenecks,
        census_times,
        model="hudson",
        max_segments=100,
        num_labels=1,
        sweep_trajectory=None,
        full_arg=False,
        time_slice=None,
        gene_conversion_rate=0.0,
        gene_conversion_length=1,
        pedigree=None,
        discrete_genome=True,
    ):
        # Must be a square matrix.
        N = len(migration_matrix)
        assert len(tables.populations) == N
        assert len(population_growth_rates) == N
        assert len(population_sizes) == N
        for j in range(N):
            assert N == len(migration_matrix[j])
            assert migration_matrix[j][j] == 0
        assert gene_conversion_length >= 1

        self.tables = tables
        self.model = model
        self.L = tables.sequence_length
        self.recomb_map = recombination_map
        self.gc_map = RateMap([0, self.L], [gene_conversion_rate, 0])
        self.tract_length = gene_conversion_length
        self.discrete_genome = discrete_genome
        self.migration_matrix = migration_matrix
        self.num_labels = num_labels
        self.num_populations = N
        self.max_segments = max_segments
        self.full_arg = full_arg
        self.segment_stack = []
        self.segments = [None for j in range(self.max_segments + 1)]
        for j in range(self.max_segments):
            s = Segment(j + 1)
            self.segments[j + 1] = s
            self.segment_stack.append(s)
        self.P = [Population(id_, num_labels) for id_ in range(N)]
        if self.recomb_map.total_mass == 0:
            self.recomb_mass_index = None
        else:
            self.recomb_mass_index = [
                FenwickTree(self.max_segments) for j in range(num_labels)
            ]
        if self.gc_map.total_mass == 0:
            self.gc_mass_index = None
        else:
            self.gc_mass_index = [
                FenwickTree(self.max_segments) for j in range(num_labels)
            ]
        self.S = bintrees.AVLTree()
        for pop in self.P:
            pop.set_start_size(population_sizes[pop.id])
            pop.set_growth_rate(population_growth_rates[pop.id], 0)
        self.edge_buffer = []
        self.pedigree = pedigree

        self.initialise(tables.tree_sequence())

        if pedigree is not None:
            assert N == 1  # <- only support single pop/pedigree for now
            pop = self.P[0]
            pedigree.load_pop(pop)

        self.num_ca_events = 0
        self.num_re_events = 0
        self.num_gc_events = 0

        # Sweep variables
        self.sweep_site = (self.L // 2) - 1  # need to add options here
        self.sweep_trajectory = sweep_trajectory
        self.time_slice = time_slice

        self.modifier_events = [(sys.float_info.max, None, None)]
        for time, pop_id, new_size in population_size_changes:
            self.modifier_events.append(
                (time, self.change_population_size, (int(pop_id), new_size))
            )
        for time, pop_id, new_rate in population_growth_rate_changes:
            self.modifier_events.append(
                (
                    time,
                    self.change_population_growth_rate,
                    (int(pop_id), new_rate, time),
                )
            )
        for time, pop_i, pop_j, new_rate in migration_matrix_element_changes:
            self.modifier_events.append(
                (
                    time,
                    self.change_migration_matrix_element,
                    (int(pop_i), int(pop_j), new_rate),
                )
            )
        for time, pop_id, intensity in bottlenecks:
            self.modifier_events.append(
                (time, self.bottleneck_event, (int(pop_id), 0, intensity))
            )
        for time in census_times:
            self.modifier_events.append((time[0], self.census_event, time))
        self.modifier_events.sort()

    def initialise(self, ts):
        root_time = np.max(self.tables.nodes.time)
        self.t = root_time

        root_segments_head = [None for _ in range(ts.num_nodes)]
        root_segments_tail = [None for _ in range(ts.num_nodes)]
        last_S = -1
        for tree in ts.trees():
            left, right = tree.interval
            S = 0 if tree.num_roots == 1 else tree.num_roots
            if S != last_S:
                self.S[left] = S
                last_S = S
            # If we have 1 root this is a special case and we don't add in
            # any ancestral segments to the state.
            if tree.num_roots > 1:
                for root in tree.roots:
                    population = ts.node(root).population
                    if root_segments_head[root] is None:
                        seg = self.alloc_segment(left, right, root, population)
                        root_segments_head[root] = seg
                        root_segments_tail[root] = seg
                    else:
                        tail = root_segments_tail[root]
                        if tail.right == left:
                            tail.right = right
                        else:
                            seg = self.alloc_segment(
                                left, right, root, population, tail
                            )
                            tail.next = seg
                            root_segments_tail[root] = seg
        self.S[self.L] = -1

        # Insert the segment chains into the algorithm state.
        for node in range(ts.num_nodes):
            seg = root_segments_head[node]
            if seg is not None:
                self.P[seg.population].add(seg)
                while seg is not None:
                    self.set_segment_mass(seg)
                    seg = seg.next

    def ancestors_remain(self):
        """
        Returns True if the simulation is not finished, i.e., there is some ancestral
        material that has not fully coalesced.
        """
        return sum(pop.get_num_ancestors() for pop in self.P) != 0

    def change_population_size(self, pop_id, size):
        self.P[pop_id].set_start_size(size)

    def change_population_growth_rate(self, pop_id, rate, time):
        self.P[pop_id].set_growth_rate(rate, time)

    def change_migration_matrix_element(self, pop_i, pop_j, rate):
        self.migration_matrix[pop_i][pop_j] = rate

    def alloc_segment(
        self,
        left,
        right,
        node,
        population,
        prev=None,
        next=None,  # noqa: A002
        label=0,
    ):
        """
        Pops a new segment off the stack and sets its properties.
        """
        s = self.segment_stack.pop()
        s.left = left
        s.right = right
        s.node = node
        s.population = population
        s.next = next
        s.prev = prev
        s.label = label
        return s

    def copy_segment(self, segment):
        return self.alloc_segment(
            left=segment.left,
            right=segment.right,
            node=segment.node,
            population=segment.population,
            next=segment.next,
            prev=segment.prev,
            label=segment.label,
        )

    def free_segment(self, u):
        """
        Frees the specified segment making it ready for reuse and
        setting its weight to zero.
        """
        if self.recomb_mass_index is not None:
            self.recomb_mass_index[u.label].set_value(u.index, 0)
        if self.gc_mass_index is not None:
            self.gc_mass_index[u.label].set_value(u.index, 0)
        self.segment_stack.append(u)

    def store_node(self, population, flags=0):
        self.flush_edges()
        self.tables.nodes.add_row(time=self.t, flags=flags, population=population)

    def flush_edges(self):
        """
        Flushes the edges in the edge buffer to the table, squashing any adjacent edges.
        """
        if len(self.edge_buffer) > 0:
            parent = len(self.tables.nodes) - 1
            self.edge_buffer.sort(key=lambda e: (e.child, e.left))
            left = self.edge_buffer[0].left
            right = self.edge_buffer[0].right
            child = self.edge_buffer[0].child
            assert self.edge_buffer[0].parent == parent
            for e in self.edge_buffer[1:]:
                assert e.parent == parent
                if e.left != right or e.child != child:
                    self.tables.edges.add_row(left, right, parent, child)
                    left = e.left
                    child = e.child
                right = e.right
            self.tables.edges.add_row(left, right, parent, child)
            self.edge_buffer = []

    def store_edge(self, left, right, parent, child):
        """
        Stores the specified edge to the output tree sequence.
        """
        self.edge_buffer.append(
            tskit.Edge(left=left, right=right, parent=parent, child=child)
        )

    def finalise(self):
        """
        Finalises the simulation returns an msprime tree sequence object.
        """
        self.flush_edges()
        ts = self.tables.tree_sequence()
        return ts

    def simulate(self, end_time):
        self.verify()
        if self.model == "hudson":
            self.hudson_simulate(end_time)
        elif self.model == "dtwf":
            self.dtwf_simulate()
        elif self.model == "wf_ped":
            self.pedigree_simulate()
        elif self.model == "single_sweep":
            self.single_sweep_simulate()
        else:
            print("Error: bad model specification -", self.model)
            raise ValueError
        return self.finalise()

    def get_potential_destinations(self):
        """
        For each population return the set of populations for which it has a
        non-zero migration into.
        """
        N = len(self.P)
        potential_destinations = [set() for _ in range(N)]
        for j in range(N):
            for k in range(N):
                if self.migration_matrix[j][k] > 0:
                    potential_destinations[j].add(k)
        return potential_destinations

    def get_total_recombination_rate(self, label):
        total_rate = 0
        if self.recomb_mass_index is not None:
            total_rate = self.recomb_mass_index[label].get_total()
        return total_rate

    def get_total_gc_rate(self, label):
        total_rate = 0
        if self.gc_mass_index is not None:
            total_rate = self.gc_mass_index[label].get_total()
        return total_rate

    def get_total_gc_left_rate(self, label):
        gc_left_total = self.get_total_gc_left(label)
        return gc_left_total

    def get_total_gc_left(self, label):
        gc_left_total = 0
        num_ancestors = sum(pop.get_num_ancestors() for pop in self.P)
        mean_gc_rate = self.gc_map.mean_rate
        gc_left_total = num_ancestors * mean_gc_rate * self.tract_length
        return gc_left_total

    def find_cleft_individual(self, label, cleft_value):
        mean_gc_rate = self.gc_map.mean_rate
        individual_index = math.floor(cleft_value / (mean_gc_rate * self.tract_length))
        for pop in self.P:
            num_ancestors = pop.get_num_ancestors()
            if individual_index < num_ancestors:
                return pop._ancestors[label][individual_index]
            individual_index -= num_ancestors
        raise AssertionError()

    def hudson_simulate(self, end_time):
        """
        Simulates the algorithm until all loci have coalesced.
        """
        infinity = sys.float_info.max
        non_empty_pops = {pop.id for pop in self.P if pop.get_num_ancestors() > 0}
        potential_destinations = self.get_potential_destinations()

        # only worried about label 0 below
        while len(non_empty_pops) > 0:
            self.verify()
            if self.t >= end_time:
                break
            # self.print_state()
            re_rate = self.get_total_recombination_rate(label=0)
            t_re = infinity
            if re_rate > 0:
                t_re = random.expovariate(re_rate)

            # Gene conversion can occur within segments ..
            gc_rate = self.get_total_gc_rate(label=0)
            t_gcin = infinity
            if gc_rate > 0:
                t_gcin = random.expovariate(gc_rate)
            # ... or to the left of the first segment.
            gc_left_rate = self.get_total_gc_left_rate(label=0)
            t_gc_left = infinity
            if gc_left_rate > 0:
                t_gc_left = random.expovariate(gc_left_rate)

            # Common ancestor events occur within demes.
            t_ca = infinity
            for index in non_empty_pops:
                pop = self.P[index]
                assert pop.get_num_ancestors() > 0
                t = pop.get_common_ancestor_waiting_time(self.t)
                if t < t_ca:
                    t_ca = t
                    ca_population = index
            t_mig = infinity
            # Migration events happen at the rates in the matrix.
            for j in non_empty_pops:
                source_size = self.P[j].get_num_ancestors()
                assert source_size > 0
                # for k in range(len(self.P)):
                for k in potential_destinations[j]:
                    rate = source_size * self.migration_matrix[j][k]
                    assert rate > 0
                    t = random.expovariate(rate)
                    if t < t_mig:
                        t_mig = t
                        mig_source = j
                        mig_dest = k
            min_time = min(t_re, t_ca, t_gcin, t_gc_left, t_mig)
            assert min_time != infinity
            if self.t + min_time > self.modifier_events[0][0]:
                t, func, args = self.modifier_events.pop(0)
                self.t = t
                func(*args)
                # Don't bother trying to maintain the non-zero lists
                # through demographic events, just recompute them.
                non_empty_pops = {
                    pop.id for pop in self.P if pop.get_num_ancestors() > 0
                }
                potential_destinations = self.get_potential_destinations()
                event = "MOD"
            else:
                self.t += min_time
                if min_time == t_re:
                    event = "RE"
                    self.hudson_recombination_event(0)
                elif min_time == t_gcin:
                    event = "GCI"
                    self.wiuf_gene_conversion_within_event(0)
                elif min_time == t_gc_left:
                    event = "GCL"
                    self.wiuf_gene_conversion_left_event(0)
                elif min_time == t_ca:
                    event = "CA"
                    self.common_ancestor_event(ca_population, 0)
                    if self.P[ca_population].get_num_ancestors() == 0:
                        non_empty_pops.remove(ca_population)
                else:
                    event = "MIG"
                    self.migration_event(mig_source, mig_dest)
                    if self.P[mig_source].get_num_ancestors() == 0:
                        non_empty_pops.remove(mig_source)
                    assert self.P[mig_dest].get_num_ancestors() > 0
                    non_empty_pops.add(mig_dest)
            logger.info(
                "%s time=%f n=%d",
                event,
                self.t,
                sum(pop.get_num_ancestors() for pop in self.P),
            )

            X = {pop.id for pop in self.P if pop.get_num_ancestors() > 0}
            assert non_empty_pops == X
        return self.finalise()

    def single_sweep_simulate(self):
        """
        Does a structed coalescent until end_freq is reached, using
        information in self.weep_trajectory.

        """
        allele_freqs, times = self.sweep_trajectory
        sweep_traj_step = 0
        x = allele_freqs[sweep_traj_step]

        assert self.num_populations == 1

        # go through segments and assign labels
        # a bit ugly with the two loops because
        # of dealing with the pops
        indices = []
        for idx, u in enumerate(self.P[0].iter_label(0)):
            if random.random() < x:
                self.set_labels(u, 1)
                indices.append(idx)
            else:
                assert u.label == 0
        popped = 0
        for i in indices:
            tmp = self.P[0].remove(i - popped, 0)
            popped += 1
            self.P[0].add(tmp, 1)

        # main loop time
        t_inc_orig = self.time_slice
        e_time = 0.0
        while self.ancestors_remain() and sweep_traj_step < len(times) - 1:
            self.verify()
            event_prob = 1.0
            while event_prob > random.random() and sweep_traj_step < len(times) - 1:
                sweep_traj_step += 1
                x = allele_freqs[sweep_traj_step]
                e_time += times[sweep_traj_step]
                # self.t = self.t + times[sweep_traj_step]
                sweep_pop_sizes = [
                    self.P[0].get_num_ancestors(label=0),
                    self.P[0].get_num_ancestors(label=1),
                ]
                # print(sweep_pop_sizes)
                p_rec_b = self.get_total_recombination_rate(0) * t_inc_orig
                p_rec_B = self.get_total_recombination_rate(1) * t_inc_orig

                # JK NOTE: We should probably factor these pop size calculations
                # into a method in Population like get_common_ancestor_waiting_time().
                # That way we can handle exponentially growing populations as well?
                p_coal_b = (
                    (sweep_pop_sizes[0] * (sweep_pop_sizes[0] - 1))
                    / (1.0 - x)
                    * t_inc_orig
                    / self.P[0].start_size
                )
                p_coal_B = (
                    (sweep_pop_sizes[1] * (sweep_pop_sizes[1] - 1))
                    / x
                    * t_inc_orig
                    / self.P[0].start_size
                )
                sweep_pop_tot_rate = p_rec_b + p_rec_B + p_coal_b + p_coal_B

                total_rate = sweep_pop_tot_rate
                if total_rate == 0:
                    break
                event_prob *= 1.0 - total_rate

            if total_rate == 0:
                break
            if self.t + e_time > self.modifier_events[0][0]:
                t, func, args = self.modifier_events.pop(0)
                self.t = t
                func(*args)
            else:
                self.t += e_time
                # choose which event happened
                # print("event time: "+str(self.t))
                if random.random() < sweep_pop_tot_rate / total_rate:
                    # even in sweeping pop, choose which kind
                    r = random.random()
                    e_sum = p_coal_B
                    if r < e_sum / sweep_pop_tot_rate:
                        # coalescent in B
                        self.common_ancestor_event(0, 1)
                    else:
                        e_sum += p_coal_b
                        if r < e_sum / sweep_pop_tot_rate:
                            # coalescent in b
                            self.common_ancestor_event(0, 0)
                        else:
                            e_sum += p_rec_B
                            if r < e_sum / sweep_pop_tot_rate:
                                # recomb in B
                                self.hudson_recombination_event_sweep_phase(
                                    1, self.sweep_site, x
                                )
                            else:
                                # recomb in b
                                self.hudson_recombination_event_sweep_phase(
                                    0, self.sweep_site, 1.0 - x
                                )
        # clean up the labels at end
        for idx, u in enumerate(self.P[0].iter_label(1)):
            tmp = self.P[0].remove(idx, u.label)
            self.set_labels(u, 0)
            self.P[0].add(tmp)

    def pedigree_simulate(self):
        """
        Simulates through the provided pedigree, stopping at the top.
        """
        assert self.pedigree is not None
        self.dtwf_climb_pedigree()

        # Complete simulations under dtwf model
        # self.dtwf_simulate

    def dtwf_simulate(self):
        """
        Simulates the algorithm until all loci have coalesced.
        """
        while self.ancestors_remain():
            self.t += 1
            self.verify()

            self.dtwf_generation()

    def dtwf_generation(self):
        """
        Evolves one generation of a Wright Fisher population
        """
        for pop_idx, pop in enumerate(self.P):
            # Cluster haploid inds by parent
            parent_inds = pop.get_ind_range(self.t)
            offspring = bintrees.AVLTree()
            for anc in pop.iter_label(0):
                parent_index = np.random.choice(parent_inds.stop - parent_inds.start)
                parent = parent_inds.start + parent_index
                if parent not in offspring:
                    offspring[parent] = []
                offspring[parent].append(anc)

            # Draw recombinations in children and sort segments by
            # inheritance direction
            for children in offspring.values():
                H = [[], []]
                for child in children:
                    segs_pair = self.dtwf_recombine(child)
                    for seg in segs_pair:
                        if seg is not None and seg.index != child.index:
                            pop.add(seg)

                    self.verify()
                    # Collect segments inherited from the same individual
                    for i, seg in enumerate(segs_pair):
                        if seg is None:
                            continue
                        assert seg.prev is None
                        heapq.heappush(H[i], (seg.left, seg))

                # Merge segments
                for h in H:
                    segments_to_merge = len(h)
                    if segments_to_merge == 1:
                        h = []
                    elif segments_to_merge >= 2:
                        for _, individual in h:
                            pop.remove_individual(individual)
                        if segments_to_merge == 2:
                            self.merge_two_ancestors(pop_idx, 0, h[0][1], h[1][1])
                        else:
                            self.merge_ancestors(h, pop_idx, 0)  # label 0 only
            self.verify()

        # Migration events happen at the rates in the matrix.
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

    def dtwf_climb_pedigree(self):
        """
        Simulates transmission of ancestral material through a pre-specified
        pedigree
        """
        assert len(self.pedigree.ind_heap) > 0
        assert self.num_populations == 1  # Single pop/pedigree for now
        self.pedigree.is_climbing = True

        # Store founders for adding back to population when climbing is done
        # TODO: Store as heap to add to population for simultaneous DTWF sims?
        founder_lineages = []

        while len(self.pedigree.ind_heap) > 0:
            next_ind = self.pedigree.pop_ind()
            self.t = next_ind.time
            assert next_ind.num_lineages() > 0
            assert next_ind.merged is False

            for segments, parent in zip(next_ind.segments, next_ind.parents):
                # This parent may not have contributed any ancestral material
                # to the samples.
                if len(segments) == 0:
                    continue

                # Merge segments inherited from this ind and recombine
                self.merge_ancestors(segments, 0, 0)
                merged_segment = self.pedigree.merged_segment
                self.pedigree.merged_segment = None
                assert merged_segment.prev is None

                # If parent is None, we are at a pedigree founder and we add
                # to founder lineages.
                if parent is None:
                    founder_lineages.append(merged_segment)
                    continue

                # Recombine and climb segments to parents.
                segs_pair = self.dtwf_recombine(merged_segment)
                for i, seg in enumerate(segs_pair):
                    if seg is None:
                        continue
                    assert seg.prev is None
                    parent.add_segment(seg, parent_ix=i)

                if parent.queued is False:
                    self.pedigree.push_ind(parent)

            next_ind.merged = True

        self.pedigree.is_climbing = False

        # Add lineages back to population.
        for lineage in founder_lineages:
            self.P[0].add(lineage)

        self.verify()

    def store_arg_edges(self, segment):
        u = len(self.tables.nodes) - 1
        # Store edges pointing to current node to the left
        x = segment
        while x is not None:
            if x.node != u:
                self.store_edge(x.left, x.right, u, x.node)
            x.node = u
            x = x.prev
        # Store edges pointing to current node to the right
        x = segment
        while x is not None:
            if x.node != u:
                self.store_edge(x.left, x.right, u, x.node)
            x.node = u
            x = x.next

    def migration_event(self, j, k):
        """
        Migrates an individual from population j to population k.
        Only does label 0
        """
        label = 0
        source = self.P[j]
        dest = self.P[k]
        index = random.randint(0, source.get_num_ancestors(label) - 1)
        x = source.remove(index, label)
        dest.add(x, label)
        if self.full_arg:
            self.store_node(k, flags=msprime.NODE_IS_MIG_EVENT)
            self.store_arg_edges(x)
        # Set the population id for each segment also.
        u = x
        while u is not None:
            u.population = k
            u = u.next

    def get_recomb_left_bound(self, seg):
        """
        Returns the left bound for genomic region over which the specified
        segment represents recombination events.
        """
        if seg.prev is None:
            left_bound = seg.left + 1 if self.discrete_genome else seg.left
        else:
            left_bound = seg.prev.right
        return left_bound

    def get_gc_left_bound(self, seg):
        # TODO remove me
        return self.get_recomb_left_bound(seg)

    def set_segment_mass(self, seg):
        """
        Sets the mass for the specified segment. All links *must* be
        appropriately set before calling this function.
        """
        if self.recomb_mass_index is not None:
            mass_index = self.recomb_mass_index[seg.label]
            recomb_left_bound = self.get_recomb_left_bound(seg)
            recomb_mass = self.recomb_map.mass_between(recomb_left_bound, seg.right)
            mass_index.set_value(seg.index, recomb_mass)
        if self.gc_mass_index is not None:
            mass_index = self.gc_mass_index[seg.label]
            gc_left_bound = self.get_gc_left_bound(seg)
            gc_mass = self.gc_map.mass_between(gc_left_bound, seg.right)
            mass_index.set_value(seg.index, gc_mass)

    def set_labels(self, segment, new_label):
        """
        Move the specified segment to the specified label.
        """
        mass_indexes = [self.recomb_mass_index, self.gc_mass_index]
        while segment is not None:
            masses = []
            for mass_index in mass_indexes:
                if mass_index is not None:
                    masses.append(mass_index[segment.label].get_value(segment.index))
                    mass_index[segment.label].set_value(segment.index, 0)
            segment.label = new_label
            for mass, mass_index in zip(masses, mass_indexes):
                if mass_index is not None:
                    mass_index[segment.label].set_value(segment.index, mass)
            segment = segment.next

    def choose_breakpoint(self, mass_index, rate_map):
        assert mass_index.get_total() > 0
        random_mass = random.uniform(0, mass_index.get_total())
        y = self.segments[mass_index.find(random_mass)]
        y_cumulative_mass = mass_index.get_cumulative_sum(y.index)
        y_right_mass = rate_map.position_to_mass(y.right)
        bp_mass = y_right_mass - (y_cumulative_mass - random_mass)
        bp = rate_map.mass_to_position(bp_mass)
        if self.discrete_genome:
            bp = math.floor(bp)
        return y, bp

    def hudson_recombination_event(self, label, return_heads=False):
        """
        Implements a recombination event.
        """
        self.num_re_events += 1
        y, bp = self.choose_breakpoint(self.recomb_mass_index[label], self.recomb_map)
        x = y.prev
        if y.left < bp:
            #   x         y
            # =====  ===|====  ...
            #          bp
            # becomes
            #   x     y
            # =====  ===  α        (LHS)
            #           =====  ... (RHS)
            alpha = self.copy_segment(y)
            alpha.left = bp
            alpha.prev = None
            if y.next is not None:
                y.next.prev = alpha
            y.next = None
            y.right = bp
            self.set_segment_mass(y)
            lhs_tail = y
        else:
            #   x            y
            # =====  |   =========  ...
            #
            # becomes
            #   x
            # =====          α          (LHS)
            #            =========  ... (RHS)
            x.next = None
            y.prev = None
            alpha = y
            lhs_tail = x
        self.set_segment_mass(alpha)
        self.P[alpha.population].add(alpha, label)
        if self.full_arg:
            self.store_node(lhs_tail.population, flags=msprime.NODE_IS_RE_EVENT)
            self.store_arg_edges(lhs_tail)
            self.store_node(alpha.population, flags=msprime.NODE_IS_RE_EVENT)
            self.store_arg_edges(alpha)
        ret = None
        if return_heads:
            x = lhs_tail
            # Seek back to the head of the x chain
            while x.prev is not None:
                x = x.prev
            ret = x, alpha
        return ret

    def generate_gc_tract_length(self):
        # generate tract length
        if self.discrete_genome:
            tl = np.random.geometric(1 / self.tract_length)
        else:
            tl = np.random.exponential(self.tract_length)
        return tl

    def wiuf_gene_conversion_within_event(self, label):
        """
        Implements a gene conversion event that starts within a segment
        """
        # TODO This is more complicated than it needs to be now because
        # we're not trying to simulate the full GC process with this
        # one event anymore. Look into what bits can be dropped now
        # that we're simulating gc_left separately again.
        y, left_breakpoint = self.choose_breakpoint(
            self.gc_mass_index[label], self.gc_map
        )
        x = y.prev
        # generate tract_length
        tl = self.generate_gc_tract_length()
        assert tl > 0
        right_breakpoint = left_breakpoint + tl
        if y.left >= right_breakpoint:
            #                  y
            # ...  |   |   ========== ...
            #     lbp rbp
            return None
        self.num_gc_events += 1

        # Process left break
        insert_alpha = True
        if left_breakpoint <= y.left:
            #  x             y
            # =====  |  ==========
            #       lbp
            #
            # becomes
            #  x
            # =====         α
            #           ==========
            if x is None:
                # In this case we *don't* insert alpha because it is already
                # the head of a segment chain
                insert_alpha = False
            else:
                x.next = None
            y.prev = None
            alpha = y
            tail = x
        else:
            #  x             y
            # =====     ====|=====
            #              lbp
            #
            # becomes
            #  x         y
            # =====     ====   α
            #               ======
            alpha = self.copy_segment(y)
            alpha.left = left_breakpoint
            alpha.prev = None
            if y.next is not None:
                y.next.prev = alpha
            y.next = None
            y.right = left_breakpoint
            self.set_segment_mass(y)
            tail = y
        self.set_segment_mass(alpha)

        # Find the segment z that the right breakpoint falls in
        z = alpha
        while z is not None and right_breakpoint >= z.right:
            z = z.next

        head = None
        # Process the right break
        if z is not None:
            if z.left < right_breakpoint:
                #   tail             z
                # ======
                #       ...  ===|==========
                #              rbp
                #
                # becomes
                #  tail              head
                # =====         ===========
                #      ...   ===
                #             z
                head = self.copy_segment(z)
                head.left = right_breakpoint
                if z.next is not None:
                    z.next.prev = head
                z.right = right_breakpoint
                z.next = None
                self.set_segment_mass(z)
            else:
                #   tail             z
                # ======
                #   ...   |   =============
                #        rbp
                #
                # becomes
                #  tail             z
                # ======      =============
                #  ...
                if z.prev is not None:
                    z.prev.next = None
                head = z
            if tail is not None:
                tail.next = head
            head.prev = tail
            self.set_segment_mass(head)

        #        y            z
        #  |  ========== ... ===== |
        # lbp                     rbp
        # When y and z are the head and tail of the segment chains, then
        # this GC event does nothing. This logic takes are of this situation.
        new_individual_head = None
        if insert_alpha:
            new_individual_head = alpha
        elif head is not None:
            new_individual_head = head
        if new_individual_head is not None:
            self.P[new_individual_head.population].add(
                new_individual_head, new_individual_head.label
            )

    def wiuf_gene_conversion_left_event(self, label):
        """
        Implements a gene conversion event that started left of a first segment.
        """
        random_gc_left = random.uniform(0, self.get_total_gc_left(label))
        # Get segment where gene conversion starts from left
        y = self.find_cleft_individual(label, random_gc_left)
        assert y is not None

        # generate tract_length
        tl = self.generate_gc_tract_length()

        bp = y.left + tl
        while y is not None and y.right <= bp:
            y = y.next

        # if the gene conversion spans the whole individual nothing happens
        if y is None:
            #    last segment
            # ... ==========   |
            #                  bp
            # stays in current state
            return None

        self.num_gc_events += 1
        x = y.prev
        if y.left < bp:
            #  x          y
            # =====   =====|====
            #              bp
            # becomes
            #  x         y
            # =====   =====
            #              =====
            #                α
            alpha = self.copy_segment(y)
            alpha.left = bp
            alpha.prev = None
            if alpha.next is not None:
                alpha.next.prev = alpha
            y.next = None
            y.right = bp
            self.set_segment_mass(y)
        else:
            #  x          y
            # ===== |  =========
            #       bp
            # becomes
            #  x
            # =====
            #          =========
            #              α
            # split the link between x and y.
            x.next = None
            y.prev = None
            alpha = y
        self.set_segment_mass(alpha)
        assert alpha.prev is None
        self.P[alpha.population].add(alpha, label)

    def hudson_recombination_event_sweep_phase(self, label, sweep_site, pop_freq):
        """
        Implements a recombination event in during a selective sweep.
        """
        lhs, rhs = self.hudson_recombination_event(label, return_heads=True)

        r = random.random()
        if sweep_site < rhs.left:
            if r < 1.0 - pop_freq:
                # move rhs to other population
                t_idx = self.P[rhs.population].find_indv(rhs)
                self.P[rhs.population].remove(t_idx, rhs.label)
                self.set_labels(rhs, 1 - label)
                self.P[rhs.population].add(rhs, rhs.label)
        else:
            if r < 1.0 - pop_freq:
                # move lhs to other population
                t_idx = self.P[lhs.population].find_indv(lhs)
                self.P[lhs.population].remove(t_idx, lhs.label)
                self.set_labels(lhs, 1 - label)
                self.P[lhs.population].add(lhs, lhs.label)

    def dtwf_generate_breakpoint(self, start):
        left_bound = start + 1 if self.discrete_genome else start
        mass_to_next_recomb = np.random.exponential(1.0)
        bp = self.recomb_map.shift_by_mass(left_bound, mass_to_next_recomb)
        if self.discrete_genome:
            bp = math.floor(bp)
        return bp

    def dtwf_recombine(self, x):
        """
        Chooses breakpoints and returns segments sorted by inheritance
        direction, by iterating through segment chain starting with x
        """
        u = self.alloc_segment(-1, -1, -1, -1, None, None)
        v = self.alloc_segment(-1, -1, -1, -1, None, None)
        seg_tails = [u, v]

        # TODO Should this be the recombination rate going foward from x.left?
        if self.recomb_map.total_mass > 0:
            k = self.dtwf_generate_breakpoint(x.left)
        else:
            k = np.inf

        ix = np.random.randint(2)
        seg_tails[ix].next = x

        while x is not None:
            seg_tails[ix] = x
            y = x.next

            if x.right > k:
                assert x.left < k
                self.num_re_events += 1
                ix = (ix + 1) % 2
                # Make new segment
                if seg_tails[ix] is u or seg_tails[ix] is v:
                    tail = None
                else:
                    tail = seg_tails[ix]
                z = self.copy_segment(x)
                z.left = k
                z.prev = tail
                self.set_segment_mass(z)
                if x.next is not None:
                    x.next.prev = z
                seg_tails[ix].next = z
                seg_tails[ix] = z
                x.next = None
                x.right = k
                self.set_segment_mass(x)
                x = z
                k = self.dtwf_generate_breakpoint(k)
            elif x.right <= k and y is not None and y.left >= k:
                # Recombine between segment and the next
                assert seg_tails[ix] == x
                x.next = None
                y.prev = None
                while y.left >= k:
                    self.num_re_events += 1
                    ix = (ix + 1) % 2
                    k = self.dtwf_generate_breakpoint(k)
                seg_tails[ix].next = y
                if seg_tails[ix] is u or seg_tails[ix] is v:
                    tail = None
                else:
                    tail = seg_tails[ix]
                y.prev = tail
                self.set_segment_mass(y)
                seg_tails[ix] = y
                x = y
            else:
                # No recombination between x.right and y.left
                x = y

        # Remove sentinel segments - this can be handled more simply
        # with pointers in C implemetation
        s = u
        u = s.next
        self.free_segment(s)

        s = v
        v = s.next
        self.free_segment(s)

        return u, v

    def census_event(self, time):
        for pop in self.P:
            for ancestor in pop.iter_ancestors():
                seg = ancestor
                self.flush_edges()
                u = self.tables.nodes.add_row(
                    time=time, flags=msprime.NODE_IS_CEN_EVENT, population=pop.id
                )
                while seg is not None:
                    # Add an edge joining the segment to the new node.
                    self.store_edge(seg.left, seg.right, u, seg.node)
                    seg.node = u
                    seg = seg.next

    def bottleneck_event(self, pop_id, label, intensity):
        # self.print_state()
        # Merge some of the ancestors.
        pop = self.P[pop_id]
        H = []
        for _ in range(pop.get_num_ancestors()):
            if random.random() < intensity:
                x = pop.remove(0)
                heapq.heappush(H, (x.left, x))
        self.merge_ancestors(H, pop_id, label)

    def merge_ancestors(self, H, pop_id, label):
        pop = self.P[pop_id]
        defrag_required = False
        coalescence = False
        alpha = None
        z = None
        while len(H) > 0:
            alpha = None
            left = H[0][0]
            X = []
            r_max = self.L
            while len(H) > 0 and H[0][0] == left:
                x = heapq.heappop(H)[1]
                X.append(x)
                r_max = min(r_max, x.right)
            if len(H) > 0:
                r_max = min(r_max, H[0][0])
            if len(X) == 1:
                x = X[0]
                if len(H) > 0 and H[0][0] < x.right:
                    alpha = self.alloc_segment(x.left, H[0][0], x.node, x.population)
                    alpha.label = label
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
                u = len(self.tables.nodes) - 1
                # We must also break if the next left value is less than
                # any of the right values in the current overlap set.
                if left not in self.S:
                    j = self.S.floor_key(left)
                    self.S[left] = self.S[j]
                if r_max not in self.S:
                    j = self.S.floor_key(r_max)
                    self.S[r_max] = self.S[j]
                # Update the number of extant segments.
                if self.S[left] == len(X):
                    self.S[left] = 0
                    right = self.S.succ_key(left)
                else:
                    right = left
                    while right < r_max and self.S[right] != len(X):
                        self.S[right] -= len(X) - 1
                        right = self.S.succ_key(right)
                    alpha = self.alloc_segment(left, right, u, pop_id)
                # Update the heaps and make the record.
                for x in X:
                    self.store_edge(left, right, u, x.node)
                    if x.right == right:
                        self.free_segment(x)
                        if x.next is not None:
                            y = x.next
                            heapq.heappush(H, (y.left, y))
                    elif x.right > right:
                        x.left = right
                        heapq.heappush(H, (x.left, x))

            # loop tail; update alpha and integrate it into the state.
            if alpha is not None:
                if z is None:
                    # Pedigrees don't currently track lineages in Populations,
                    # so keep reference to merged segments instead.
                    if self.pedigree is not None and self.pedigree.is_climbing:
                        assert self.pedigree.merged_segment is None
                        self.pedigree.merged_segment = alpha
                    else:
                        pop.add(alpha, label)
                else:
                    if self.full_arg:
                        defrag_required |= z.right == alpha.left
                    else:
                        defrag_required |= (
                            z.right == alpha.left and z.node == alpha.node
                        )
                    z.next = alpha
                alpha.prev = z
                self.set_segment_mass(alpha)
                z = alpha
        if self.full_arg:
            if not coalescence:
                self.store_node(pop_id, flags=msprime.NODE_IS_CA_EVENT)
            self.store_arg_edges(z)
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
                self.set_segment_mass(x)
                self.free_segment(y)
            y = x

    def defrag_breakpoints(self):
        # Defrag the breakpoints set
        j = 0
        k = 0
        while k < self.L:
            k = self.S.succ_key(j)
            if self.S[j] == self.S[k]:
                del self.S[k]
            else:
                j = k

    def common_ancestor_event(self, population_index, label):
        """
        Implements a coancestry event.
        """
        pop = self.P[population_index]
        # Choose two ancestors uniformly.
        j = random.randint(0, pop.get_num_ancestors(label) - 1)
        x = pop.remove(j, label)
        j = random.randint(0, pop.get_num_ancestors(label) - 1)
        y = pop.remove(j, label)
        self.merge_two_ancestors(population_index, label, x, y)

    def merge_two_ancestors(self, population_index, label, x, y):
        pop = self.P[population_index]
        self.num_ca_events += 1
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
                    alpha = self.copy_segment(x)
                    alpha.prev = None
                    alpha.next = None
                    alpha.right = y.left
                    x.left = y.left
                else:
                    if not coalescence:
                        coalescence = True
                        self.store_node(population_index)
                    u = len(self.tables.nodes) - 1
                    # Put in breakpoints for the outer edges of the coalesced
                    # segment
                    left = x.left
                    r_max = min(x.right, y.right)
                    if left not in self.S:
                        j = self.S.floor_key(left)
                        self.S[left] = self.S[j]
                    if r_max not in self.S:
                        j = self.S.floor_key(r_max)
                        self.S[r_max] = self.S[j]
                    # Update the number of extant segments.
                    if self.S[left] == 2:
                        self.S[left] = 0
                        right = self.S.succ_key(left)
                    else:
                        right = left
                        while right < r_max and self.S[right] != 2:
                            self.S[right] -= 1
                            right = self.S.succ_key(right)
                        alpha = self.alloc_segment(
                            left=left,
                            right=right,
                            node=u,
                            population=population_index,
                            label=label,
                        )
                    self.store_edge(left, right, u, x.node)
                    self.store_edge(left, right, u, y.node)
                    # Now trim the ends of x and y to the right sizes.
                    if x.right == right:
                        self.free_segment(x)
                        x = x.next
                    else:
                        x.left = right
                    if y.right == right:
                        self.free_segment(y)
                        y = y.next
                    else:
                        y.left = right

            # loop tail; update alpha and integrate it into the state.
            if alpha is not None:
                if z is None:
                    pop.add(alpha, label)
                else:
                    if self.full_arg:
                        defrag_required |= z.right == alpha.left
                    else:
                        defrag_required |= (
                            z.right == alpha.left and z.node == alpha.node
                        )
                    z.next = alpha
                alpha.prev = z
                self.set_segment_mass(alpha)
                z = alpha

        if self.full_arg:
            if not coalescence:
                self.store_node(population_index, flags=msprime.NODE_IS_CA_EVENT)
            self.store_arg_edges(z)
        if defrag_required:
            self.defrag_segment_chain(z)
        if coalescence:
            self.defrag_breakpoints()

    def print_state(self, verify=False):
        print("State @ time ", self.t)
        for label in range(self.num_labels):
            print(
                "Recomb mass = ",
                0
                if self.recomb_mass_index is None
                else self.recomb_mass_index[label].get_total(),
            )
            print(
                "GC mass = ",
                0
                if self.gc_mass_index is None
                else self.gc_mass_index[label].get_total(),
            )
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
        for label in range(self.num_labels):
            if self.recomb_mass_index is not None:
                print(
                    "recomb_mass_index [%d]: %d"
                    % (label, self.recomb_mass_index[label].get_total())
                )
                for j in range(1, self.max_segments + 1):
                    s = self.recomb_mass_index[label].get_value(j)
                    if s != 0:
                        seg = self.segments[j]
                        left_bound = self.get_recomb_left_bound(seg)
                        sp = self.recomb_map.mass_between(left_bound, seg.right)
                        print("\t", j, "->", s, sp)
            if self.gc_mass_index is not None:
                print(
                    "gc_mass_index [%d]: %d"
                    % (label, self.gc_mass_index[label].get_total())
                )
                for j in range(1, self.max_segments + 1):
                    s = self.gc_mass_index[label].get_value(j)
                    if s != 0:
                        seg = self.segments[j]
                        left_bound = self.get_gc_left_bound(seg)
                        sp = self.gc_map.mass_between(left_bound, seg.right)
                        print("\t", j, "->", s, sp)
        print("nodes")
        print(self.tables.nodes)
        print("edges")
        print(self.tables.edges)
        if verify:
            self.verify()

    def verify_segments(self):
        for pop in self.P:
            for label in range(self.num_labels):
                for head in pop.iter_label(label):
                    assert head.prev is None
                    prev = head
                    u = head.next
                    while u is not None:
                        assert prev.next is u
                        assert u.prev is prev
                        assert u.left >= prev.right
                        assert u.label == head.label
                        assert u.population == head.population
                        prev = u
                        u = u.next

    def verify_overlaps(self):
        overlap_counter = OverlapCounter(self.L)
        for pop in self.P:
            for label in range(self.num_labels):
                for u in pop.iter_label(label):
                    while u is not None:
                        overlap_counter.increment_interval(u.left, u.right)
                        u = u.next

        for pos, count in self.S.items():
            if pos != self.L:
                assert count == overlap_counter.overlaps_at(pos)

        assert self.S[self.L] == -1
        # Check the ancestry tracking.
        A = bintrees.AVLTree()
        A[0] = 0
        A[self.L] = -1
        for pop in self.P:
            for label in range(self.num_labels):
                for u in pop.iter_label(label):
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
        while k < self.L:
            k = A.succ_key(j)
            if A[j] == A[k]:
                del A[k]
            else:
                j = k
        assert list(A.items()) == list(self.S.items())

    def verify_mass_index(self, label, mass_index, rate_map, compute_left_bound):
        assert mass_index is not None
        total_mass = 0
        alt_total_mass = 0
        for pop_index, pop in enumerate(self.P):
            for u in pop.iter_label(label):
                assert u.prev is None
                left = compute_left_bound(u)
                while u is not None:
                    assert u.population == pop_index
                    assert u.left < u.right
                    left_bound = compute_left_bound(u)
                    s = rate_map.mass_between(left_bound, u.right)
                    right = u.right
                    index_value = mass_index.get_value(u.index)
                    total_mass += index_value
                    assert math.isclose(s, index_value, abs_tol=1e-6)
                    v = u.next
                    if v is not None:
                        assert v.prev == u
                        assert u.right <= v.left
                    u = v

                s = rate_map.mass_between(left, right)
                alt_total_mass += s
        assert math.isclose(total_mass, mass_index.get_total(), abs_tol=1e-6)
        assert math.isclose(total_mass, alt_total_mass, abs_tol=1e-6)

    def verify(self):
        """
        Checks that the state of the simulator is consistent.
        """
        self.verify_segments()
        if self.model != "wf_ped":
            # The wf_ped model doesn't maintain a bunch of stuff. It would probably
            # be simpler if it did.
            self.verify_overlaps()
            for label in range(self.num_labels):
                if self.recomb_mass_index is None:
                    assert self.recomb_map.total_mass == 0
                else:
                    self.verify_mass_index(
                        label,
                        self.recomb_mass_index[label],
                        self.recomb_map,
                        self.get_recomb_left_bound,
                    )

                if self.gc_mass_index is None:
                    assert self.gc_map.total_mass == 0
                else:
                    self.verify_mass_index(
                        label,
                        self.gc_mass_index[label],
                        self.gc_map,
                        self.get_gc_left_bound,
                    )


def run_simulate(args):
    """
    Runs the simulation and outputs the results in text.
    """
    n = args.sample_size
    m = args.sequence_length
    rho = args.recombination_rate
    gc_rate = args.gene_conversion_rate[0]
    mean_tract_length = args.gene_conversion_rate[1]
    num_populations = args.num_populations
    migration_matrix = [
        [args.migration_rate * int(j != k) for j in range(num_populations)]
        for k in range(num_populations)
    ]
    sample_configuration = [0 for j in range(num_populations)]
    population_growth_rates = [0 for j in range(num_populations)]
    population_sizes = [1 for j in range(num_populations)]
    if args.pedigree_file is not None:
        n *= 2
    sample_configuration[0] = n
    if args.sample_configuration is not None:
        sample_configuration = args.sample_configuration
    if args.population_growth_rates is not None:
        population_growth_rates = args.population_growth_rates
    if args.population_sizes is not None:
        population_sizes = args.population_sizes
    if args.recomb_positions is None or args.recomb_rates is None:
        recombination_map = RateMap([0, m], [rho, 0])
    else:
        positions = args.recomb_positions
        rates = args.recomb_rates
        recombination_map = RateMap(positions, rates)
    num_labels = 1
    sweep_trajectory = None
    if args.model == "single_sweep":
        if num_populations > 1:
            raise ValueError("Multiple populations not currently supported")
        # Compute the trajectory
        if args.trajectory is None:
            raise ValueError("Must provide trajectory (init_freq, end_freq, alpha)")
        init_freq, end_freq, alpha = args.trajectory
        traj_sim = TrajectorySimulator(init_freq, end_freq, alpha, args.time_slice)
        sweep_trajectory = traj_sim.run()
        num_labels = 2
    pedigree = None
    if args.pedigree_file is not None:
        # TODO we should be ignoring the sample number and instead reading
        # the is_sample status from the file.
        if n % 2 != 0:
            raise ValueError(
                "Must specify an even number of sample lineages for "
                "diploid simulations"
            )

        num_diploid_individuals = n // 2
        py_pedigree = msprime.Pedigree.read_txt(args.pedigree_file)
        py_pedigree.set_samples(num_diploid_individuals)
        ll_pedigree = py_pedigree.get_ll_representation()
        pedigree = Pedigree(py_pedigree.num_individuals, py_pedigree.ploidy)
        pedigree.set_pedigree(
            ll_pedigree["individual"],
            ll_pedigree["parents"],
            ll_pedigree["time"],
            ll_pedigree["is_sample"],
        )
    random.seed(args.random_seed)
    np.random.seed(args.random_seed + 1)

    if args.from_ts is None:
        tables = tskit.TableCollection(m)
        if pedigree is None:
            for pop_id, sample_count in enumerate(sample_configuration):
                tables.populations.add_row()
                for _ in range(sample_count):
                    tables.nodes.add_row(
                        flags=tskit.NODE_IS_SAMPLE, time=0, population=pop_id
                    )
        else:
            # Assume one population for now.
            tables.populations.add_row()
            for is_sample, ind in zip(pedigree.is_sample, pedigree.inds):
                # TODO add information about the individual.
                ind_id = tables.individuals.add_row()
                if is_sample:
                    for _ in range(ind.ploidy):
                        tables.nodes.add_row(
                            flags=tskit.NODE_IS_SAMPLE,
                            time=ind.time,
                            population=0,
                            individual=ind_id,
                        )
    else:
        from_ts = tskit.load(args.from_ts)
        tables = from_ts.dump_tables()

    s = Simulator(
        tables=tables,
        recombination_map=recombination_map,
        migration_matrix=migration_matrix,
        model=args.model,
        population_growth_rates=population_growth_rates,
        population_sizes=population_sizes,
        population_growth_rate_changes=args.population_growth_rate_change,
        population_size_changes=args.population_size_change,
        migration_matrix_element_changes=args.migration_matrix_element_change,
        bottlenecks=args.bottleneck,
        census_times=args.census_time,
        max_segments=100000,
        num_labels=num_labels,
        full_arg=args.full_arg,
        sweep_trajectory=sweep_trajectory,
        time_slice=args.time_slice,
        gene_conversion_rate=gc_rate,
        gene_conversion_length=mean_tract_length,
        pedigree=pedigree,
        discrete_genome=args.discrete,
    )
    ts = s.simulate(args.end_time)
    ts.dump(args.output_file)
    if args.verbose:
        s.print_state()


def add_simulator_arguments(parser):
    parser.add_argument("sample_size", type=int)
    parser.add_argument("output_file")
    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )
    parser.add_argument(
        "-l",
        "--log-level",
        type=int,
        help="Set log-level to the specified value",
        default=logging.WARNING,
    )
    parser.add_argument("--random-seed", "-s", type=int, default=1)
    parser.add_argument("--sequence-length", "-L", type=int, default=100)
    parser.add_argument("--discrete", "-d", action="store_true")
    parser.add_argument("--num-replicates", "-R", type=int, default=1000)
    parser.add_argument("--recombination-rate", "-r", type=float, default=0.01)
    parser.add_argument("--recomb-positions", type=float, nargs="+", default=None)
    parser.add_argument("--recomb-rates", type=float, nargs="+", default=None)
    parser.add_argument(
        "--gene-conversion-rate", "-c", type=float, nargs=2, default=[0, 3]
    )
    parser.add_argument("--num-populations", "-p", type=int, default=1)
    parser.add_argument("--migration-rate", "-g", type=float, default=1)
    parser.add_argument("--sample-configuration", type=int, nargs="+", default=None)
    parser.add_argument(
        "--population-growth-rates", type=float, nargs="+", default=None
    )
    parser.add_argument("--population-sizes", type=float, nargs="+", default=None)
    parser.add_argument(
        "--population-size-change", type=float, nargs=3, action="append", default=[]
    )
    parser.add_argument(
        "--population-growth-rate-change",
        type=float,
        nargs=3,
        action="append",
        default=[],
    )
    parser.add_argument(
        "--migration-matrix-element-change",
        type=float,
        nargs=4,
        action="append",
        default=[],
    )
    parser.add_argument(
        "--bottleneck", type=float, nargs=3, action="append", default=[]
    )
    parser.add_argument(
        "--census-time", type=float, nargs=1, action="append", default=[]
    )
    parser.add_argument(
        "--trajectory",
        type=float,
        nargs=3,
        default=None,
        help="Parameters for the allele frequency trajectory simulation",
    )
    parser.add_argument(
        "--full-arg",
        action="store_true",
        default=False,
        help="Store the full ARG with all recombination and common ancestor nodes",
    )
    parser.add_argument(
        "--time-slice",
        type=float,
        default=1e-6,
        help="The delta_t value for selective sweeps",
    )
    parser.add_argument("--model", default="hudson")
    parser.add_argument("--pedigree-file", default=None)
    parser.add_argument(
        "--from-ts",
        "-F",
        default=None,
        help=(
            "Specify the tree sequence to complete. The sample_size argument "
            "is ignored if this is provided"
        ),
    )
    parser.add_argument(
        "--end-time", type=float, default=np.inf, help="The end for simulations."
    )


def main(args=None):
    parser = argparse.ArgumentParser()
    add_simulator_arguments(parser)
    args = parser.parse_args(args)
    log_output = daiquiri.output.Stream(
        sys.stdout,
        formatter=daiquiri.formatter.ColorFormatter(fmt="[%(levelname)s] %(message)s"),
    )
    daiquiri.setup(level=args.log_level, outputs=[log_output])
    run_simulate(args)


if __name__ == "__main__":
    main()
