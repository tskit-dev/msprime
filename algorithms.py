"""
Python version of the simulation algorithm.
"""
import argparse
import heapq
import math
import random
import sys

import bintrees
import numpy as np
import tskit

import msprime


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
        # Compute the binary logarithm of max_index
        u = self.__max_index
        while u != 0:
            self.__log_max_index = u
            u -= u & -u

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
            j += j & -j

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
            j -= j & -j
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


class Segment:
    """
    A class representing a single segment. Each segment has a left
    and right, denoting the loci over which it spans, a node and a
    next, giving the next in the chain.
    """

    def __init__(self, index):
        self.left = None
        self.right = None
        self.left_mass = None
        self.right_mass = None
        self.node = None
        self.prev = None
        self.next = None
        self.population = None
        self.label = 0
        self.index = index

    def __str__(self):
        s = "({}:{}-{}->{}: prev={} next={})".format(
            self.index,
            self.left,
            self.right,
            self.node,
            repr(self.prev),
            repr(self.next),
        )
        return s

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
                s = ""
                while u is not None:
                    s += "({}-{}->{}({});lab:{})".format(
                        u.left, u.right, u.node, u.index, u.label
                    )
                    u = u.next
                print("\t\t" + s)

    def get_cleft(self, tracklength):
        cleft = 0
        for ancestors in self._ancestors:
            for u in ancestors:
                left = u.left
                while u.next is not None:
                    u = u.next
                right = u.right
                dist = right - left
                cleft += 1 - ((tracklength - 1) / tracklength) ** (dist - 1)
        return cleft

    def find_cleft(self, rvalue, tracklength):
        for ancestors in self._ancestors:
            for u in ancestors:
                left = u.left
                index = u.index
                while u.next is not None:
                    u = u.next
                right = u.right
                dist = right - left
                rvalue -= 1 - ((tracklength - 1) / tracklength) ** (dist - 1)
                if rvalue <= 0:
                    break
            return rvalue, index, dist

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
        first_ind = np.sum([self.get_size(t_prev) for t_prev in range(0, t)])
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
        Iterates over all ancestors in a population.
        """
        for ancestors in self._ancestors:
            yield from ancestors

    def find_indv(self, indv):
        """
        find the index of an ancestor in population
        """
        return self._ancestors[indv.label].index(indv)


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


class RecombinationMap:
    def __init__(self, positions, rates, discrete):
        self.positions = positions
        self.rates = rates
        self.discrete = discrete
        self.cumulative = RecombinationMap.recomb_mass(positions, rates)

    @staticmethod
    def recomb_mass(positions, rates):
        recomb_mass = 0
        cumulative = [recomb_mass]
        for i in range(1, len(positions)):
            recomb_mass += (positions[i] - positions[i - 1]) * rates[i - 1]
            cumulative.append(recomb_mass)
        return cumulative

    @property
    def total_recombination_rate(self):
        return self.cumulative[-1]

    def mass_between(self, left, right):
        left_mass = self.position_to_mass(left)
        right_mass = self.position_to_mass(right)
        return right_mass - left_mass

    def mass_between_left_exclusive(self, left, right):
        left_bound = left + 1 if self.discrete else left
        return self.mass_between(left_bound, right)

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
        return math.floor(pos) if self.discrete else pos

    def shift_by_mass(self, pos, mass):
        result_mass = self.position_to_mass(pos) + mass
        return self.mass_to_position(result_mass)

    def sample_poisson(self, start):
        left_bound = start + 1 if self.discrete else start
        mass_to_next_recomb = np.random.exponential(1.0)
        return self.shift_by_mass(left_bound, mass_to_next_recomb)

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

    def _split(self, seg, breakpoint):  # noqa: A002
        """
        Split the segment at breakpoint and add in another segment
        from breakpoint to seg.right. Set the original segment's
        right endpoint to breakpoint
        """
        right = self._make_segment(breakpoint, seg.right, seg.node)
        if seg.next is not None:
            seg.next.prev = right
            right.next = seg.next
        right.prev = seg
        seg.next = right
        seg.right = breakpoint

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
        sample_size,
        num_loci,
        recombination_rate,
        recombination_map,
        migration_matrix,
        sample_configuration,
        population_growth_rates,
        population_sizes,
        population_growth_rate_changes,
        population_size_changes,
        migration_matrix_element_changes,
        bottlenecks,
        census_times,
        model="hudson",
        from_ts=None,
        max_segments=100,
        num_labels=1,
        sweep_trajectory=None,
        full_arg=False,
        time_slice=None,
        gene_conversion_rate=0.0,
        gene_conversion_length=1,
        pedigree=None,
    ):
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
        self.recomb_map = recombination_map
        self.g = gene_conversion_rate
        self.tracklength = gene_conversion_length
        self.pc = (self.tracklength - 1) / self.tracklength
        if self.tracklength == 1:
            self.lnpc = -math.inf
        else:
            self.lnpc = math.log(1.0 - 1.0 / self.tracklength)
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
        self.L = [FenwickTree(self.max_segments) for j in range(num_labels)]
        self.S = bintrees.AVLTree()
        for pop in self.P:
            pop.set_start_size(population_sizes[pop.id])
            pop.set_growth_rate(population_growth_rates[pop.id], 0)
        self.edge_buffer = []
        self.from_ts = from_ts
        self.pedigree = pedigree

        if from_ts is None:
            self.tables = msprime.TableCollection(sequence_length=num_loci)
            for pop_index in range(N):
                self.tables.populations.add_row()
                sample_size = sample_configuration[pop_index]
                for _ in range(sample_size):
                    j = len(self.tables.nodes)
                    x = self.alloc_segment(
                        0,
                        self.m,
                        0,
                        self.recomb_map.position_to_mass(self.m),
                        j,
                        pop_index,
                    )
                    self.set_single_segment_mass(x)
                    self.P[pop_index].add(x)
                    self.tables.nodes.add_row(
                        flags=msprime.NODE_IS_SAMPLE, time=0, population=pop_index
                    )
                    j += 1
            self.S[0] = self.n
            self.S[self.m] = -1
            self.t = 0
        else:
            ts = msprime.load(from_ts)
            if ts.sequence_length != self.m:
                raise ValueError("Sequence length in from_ts must match")
            if ts.num_populations != N:
                raise ValueError("Number of populations in from_ts must match")
            self.initialise_from_ts(ts)

        if pedigree is not None:
            assert N == 1  # <- only support single pop/pedigree for now
            pop = self.P[0]
            pedigree.load_pop(pop)

        self.num_ca_events = 0
        self.num_re_events = 0
        self.num_gc_events = 0

        # Sweep variables
        self.sweep_site = (self.m // 2) - 1  # need to add options here
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

    def initialise_from_ts(self, ts):
        self.tables = ts.dump_tables()
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
                        seg = self.alloc_segment(
                            left,
                            right,
                            self.recomb_map.position_to_mass(left),
                            self.recomb_map.position_to_mass(right),
                            root,
                            population,
                        )
                        root_segments_head[root] = seg
                        root_segments_tail[root] = seg
                    else:
                        tail = root_segments_tail[root]
                        if tail.right == left:
                            tail.right = right
                        else:
                            seg = self.alloc_segment(
                                left,
                                right,
                                self.recomb_map.position_to_mass(left),
                                self.recomb_map.position_to_mass(right),
                                root,
                                population,
                                tail,
                            )
                            tail.next = seg
                            root_segments_tail[root] = seg
        self.S[self.m] = -1

        # Insert the segment chains into the algorithm state.
        for node in range(ts.num_nodes):
            seg = root_segments_head[node]
            if seg is not None:
                self.L.set_value(seg.index, seg.right - seg.left - 1)
                self.P[seg.population].add(seg)
                prev = seg
                seg = seg.next
                while seg is not None:
                    self.L.set_value(seg.index, seg.right - prev.right)
                    prev = seg
                    seg = seg.next

    def ancestors_remain(self):
        """
        Returns True if the simulation is not finished, i.e., there is some ancestral
        material that has not fully coalesced.
        """
        return sum(pop.get_num_ancestors() for pop in self.P) != 0

    def change_population_size(self, pop_id, size):
        print("Changing pop size to ", size)
        for i in range(self.num_labels):
            self.P[i][pop_id].set_start_size(size)

    def change_population_growth_rate(self, pop_id, rate, time):
        print("Changing growth rate to ", rate)
        for i in range(self.num_labels):
            self.P[i][pop_id].set_growth_rate(rate, time)

    def change_migration_matrix_element(self, pop_i, pop_j, rate):
        print("Changing migration rate", pop_i, pop_j, rate)
        self.migration_matrix[pop_i][pop_j] = rate

    def get_cleft_total(self, tracklength):
        cleft = 0
        for pop in self.P:
            cleft += pop.get_cleft(tracklength)
        return cleft

    def find_cleft_individual(self, rvalue, tracklength):
        for pop in self.P:
            if rvalue > 0:
                rvalue, index, distance = pop.find_cleft(rvalue, tracklength)
        return index, distance

    def alloc_segment(
        self,
        left,
        right,
        left_mass,
        right_mass,
        node,
        pop_index,
        prev=None,
        next=None,  # noqa: A002
    ):
        """
        Pops a new segment off the stack and sets its properties.
        """
        s = self.segment_stack.pop()
        s.left = left
        s.right = right
        s.left_mass = left_mass
        s.right_mass = right_mass
        s.node = node
        s.population = pop_index
        s.next = next
        s.prev = prev
        s.label = 0
        return s

    def free_segment(self, u):
        """
        Frees the specified segment making it ready for reuse and
        setting its weight to zero.
        """
        self.L[u.label].set_value(u.index, 0)
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
            msprime.Edge(left=left, right=right, parent=parent, child=child)
        )

    def finalise(self):
        """
        Finalises the simulation returns an msprime tree sequence object.
        """
        self.flush_edges()
        ts = self.tables.tree_sequence()
        return ts

    def simulate(self, model="hudson"):
        if self.model == "hudson":
            self.hudson_simulate()
        elif self.model == "dtwf":
            self.dtwf_simulate()
        elif self.model == "wf_ped":
            self.pedigree_simulate()
        elif self.model == "single_sweep":
            # self.print_state()
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

    def hudson_simulate(self):
        """
        Simulates the algorithm until all loci have coalesced.
        """
        infinity = sys.float_info.max
        non_empty_pops = {pop.id for pop in self.P if pop.get_num_ancestors() > 0}
        potential_destinations = self.get_potential_destinations()

        # only worried about label 0 below
        while len(non_empty_pops) > 0:
            self.verify()
            recomb_mass = self.L[0].get_total()
            rate = recomb_mass
            t_re = infinity
            if rate != 0:
                t_re = random.expovariate(rate)
            # Gene conversion can occur within segments ..
            rate = self.g * self.recomb_map.mass_to_position(recomb_mass)
            t_gcin = infinity
            if rate != 0:
                t_gcin = random.expovariate(rate)
            # .. or left of the first segment
            cleft = self.get_cleft_total(self.tracklength)
            assert cleft <= sum(pop.get_num_ancestors() for pop in self.P)
            rate = self.g * self.tracklength * cleft
            t_gcleft = infinity
            if rate != 0:
                t_gcleft = random.expovariate(rate)
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
            min_time = min(t_re, t_ca, t_gcin, t_gcleft, t_mig)
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
            else:
                self.t += min_time
                if min_time == t_re:
                    # print("RE EVENT")
                    self.hudson_recombination_event(0)
                elif min_time == t_gcin:
                    # print("GCI EVENT")
                    self.wiuf_geneconversion_within_event(0)
                elif min_time == t_gcleft:
                    # print("GCL EVENT")
                    self.wiuf_geneconversion_left_event(0)
                elif min_time == t_ca:
                    # print("CA EVENT")
                    self.common_ancestor_event(ca_population, 0)
                    if self.P[ca_population].get_num_ancestors() == 0:
                        non_empty_pops.remove(ca_population)
                else:
                    # print("MIG EVENT")
                    self.migration_event(mig_source, mig_dest)
                    if self.P[mig_source].get_num_ancestors() == 0:
                        non_empty_pops.remove(mig_source)
                    assert self.P[mig_dest].get_num_ancestors() > 0
                    non_empty_pops.add(mig_dest)

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
                p_rec_b = self.L[0].get_total() * t_inc_orig
                p_rec_B = self.L[1].get_total() * t_inc_orig

                # JK NOTE: We should probably factor these pop size calculations
                # into a method in Population like get_common_ancestor_waiting_time().
                # That way we can handle exponentially growing populations as well?
                p_coal_b = (
                    (sweep_pop_sizes[0] * (sweep_pop_sizes[0] - 1))
                    / (1.0 - x)
                    * t_inc_orig
                    / self.P[0]._start_size
                )
                p_coal_B = (
                    (sweep_pop_sizes[1] * (sweep_pop_sizes[1] - 1))
                    / x
                    * t_inc_orig
                    / self.P[0]._start_size
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

    def set_segment_left_endpoint(self, seg, pos):
        seg.left = pos
        seg.left_mass = self.recomb_map.position_to_mass(pos)

    def add_segment_mass_between(self, seg1, l_mass, r_mass):
        self.L[seg1.label].increment(seg1.index, r_mass - l_mass)

    # Add the mass subtended by the endpoints of seg2 to that of seg1
    def add_segment_mass(self, seg1, seg2):
        mass = seg2.right_mass - seg2.left_mass
        self.L[seg1.label].increment(seg1.index, mass)

    # Subtract the mass subtended by the endpoints of seg2 to that of seg1
    def subtract_segment_mass(self, seg1, seg2):
        mass = seg2.left_mass - seg2.right_mass
        self.L[seg1.label].increment(seg1.index, mass)

    def set_segment_mass(self, seg, tail_seg):
        mass = seg.right_mass - tail_seg.right_mass
        self.L[seg.label].set_value(seg.index, mass)

    def set_single_segment_mass(self, seg):
        mass = self.recomb_map.mass_between_left_exclusive(seg.left, seg.right)
        self.L[seg.label].set_value(seg.index, mass)

    def pick_segments_and_breakpoint(self, label):
        h = random.uniform(0, self.L[label].get_total())
        y = self.segments[self.L[label].find(h)]

        t = self.L[label].get_cumulative_frequency(y.index)
        k = self.recomb_map.mass_to_position(y.right_mass - (t - h))
        if k == y.left and y.prev is None:
            return self.pick_segments_and_breakpoint(label)

        return y.prev, y, k

    def hudson_recombination_event(self, label, return_heads=False):
        """
        Implements a recombination event.
        """
        self.num_re_events += 1
        x, y, k = self.pick_segments_and_breakpoint(label)
        k_mass = self.recomb_map.position_to_mass(k)
        if y.left < k:
            # Make new segment
            z = self.alloc_segment(
                k, y.right, k_mass, y.right_mass, y.node, y.population, None, y.next
            )
            if y.next is not None:
                y.next.prev = z
            y.next = None
            y.right = k
            y.right_mass = k_mass
            self.subtract_segment_mass(y, z)
            lhs_tail = y
        else:
            # split the link between x and y.
            x.next = None
            y.prev = None
            z = y
            lhs_tail = x
        z.label = label
        self.set_single_segment_mass(z)
        self.P[z.population].add(z, label)
        if self.full_arg:
            self.store_node(lhs_tail.population, flags=msprime.NODE_IS_RE_EVENT)
            self.store_arg_edges(lhs_tail)
            self.store_node(z.population, flags=msprime.NODE_IS_RE_EVENT)
            self.store_arg_edges(z)
        ret = None
        if return_heads:
            x = lhs_tail
            # Seek back to the head of the x chain
            while x.prev is not None:
                x = x.prev
            ret = x, z
        return ret

    def cut_right_break(self, lhs_tail, y, new_segment, track_end, label):
        assert lhs_tail is not None
        lhs_tail.next = new_segment
        self.set_segment_mass(new_segment, lhs_tail)
        if y.next is not None:
            y.next.prev = new_segment
        y.next = None
        y.right = track_end
        y.right_mass = self.recomb_map.position_to_mass(track_end)
        self.add_segment_mass_between(y, new_segment.right_mass, y.right_mass)

    def wiuf_geneconversion_within_event(self, label, return_heads=False):
        """
        Implements a gene conversion event that starts within a segment
        """
        h = random.uniform(0, self.L[label].get_total())
        # generate tracklength
        tl = np.random.geometric(1 / self.tracklength)
        # Get the segment containing the h'th link
        y = self.segments[self.L[label].find(h)]
        t = self.L[label].get_cumulative_frequency(y.index)
        k = self.recomb_map.shift_by_mass(y.right, h - t)
        k_plus_tl_mass = self.recomb_map.position_to_mass(k + tl)
        k_mass = self.recomb_map.position_to_mass(k)
        # check if the gene conversion falls between segments --> no effect
        if y.left >= k + tl:
            # print("noneffective GCI EVENT")
            return None
        self.num_gc_events += 1
        x = y.prev
        # both breaks are within the same segment
        if k + tl < y.right:
            if k <= y.left:
                y.prev = None
                z2 = self.alloc_segment(
                    k + tl,
                    y.right,
                    k_plus_tl_mass,
                    y.right_mass,
                    y.node,
                    y.population,
                    x,
                    y.next,
                )
                lhs_tail = x
                self.cut_right_break(lhs_tail, y, z2, k + tl, label)
                z = y
            elif k > y.left:
                z = self.alloc_segment(
                    k, k + tl, k_mass, k_plus_tl_mass, y.node, y.population, None, None
                )
                z2 = self.alloc_segment(
                    k + tl,
                    y.right,
                    k_plus_tl_mass,
                    y.right_mass,
                    y.node,
                    y.population,
                    y,
                    y.next,
                )
                if y.next is not None:
                    y.next.prev = z2
                y.next = z2
                y.right = k
                y.right_mass = k_mass
                self.set_segment_mass(z2, y)
                self.subtract_segment_mass(y, z)
                self.subtract_segment_mass(y, z2)
                lhs_tail = y
        # breaks are in separate segments
        else:
            # Get the segment y2 containing the end of the conversion tract
            y2 = y
            while y2 is not None and k + tl >= y2.right:
                y2 = y2.next
            # process left break
            if k <= y.left:
                if x is not None:
                    x.next = None
                y.prev = None
                z = y
                lhs_tail = x
            elif k > y.left:
                z = self.alloc_segment(
                    k, y.right, k_mass, y.right_mass, y.node, y.population, None, y.next
                )
                self.set_single_segment_mass(z)
                if y.next is not None:
                    y.next.prev = z
                y.next = None
                y.right = k
                y.right_mass = k_mass
                self.subtract_segment_mass(y, z)
                lhs_tail = y
            # process right break
            if y2 is not None:
                if y2.left < k + tl:
                    z2 = self.alloc_segment(
                        k + tl,
                        y2.right,
                        k_plus_tl_mass,
                        y2.right_mass,
                        y2.node,
                        y2.population,
                        lhs_tail,
                        y2.next,
                    )
                    self.cut_right_break(lhs_tail, y2, z2, k + tl, label)
                    if z2.prev is None:
                        z = z2
                elif y2.left >= k + tl:
                    lhs_tail.next = y2
                    y2.prev.next = None
                    y2.prev = lhs_tail
                    self.set_segment_mass(y2, lhs_tail)
        # update population
        z.label = label
        self.set_single_segment_mass(z)
        self.P[z.population].add(z, label)
        # TODO check what needs to be added for full arg
        ret = None
        if return_heads:
            x = lhs_tail
            # Seek back to the head of the x chain
            while x.prev is not None:
                x = x.prev
            ret = x, z
        return ret

    def wiuf_geneconversion_left_event(self, label, return_heads=False):
        """
        Implements a gene conversion event that started left of a first segment.
        """
        self.num_gc_events += 1
        h = random.uniform(0, self.get_cleft_total(self.tracklength))
        # Get segment where gene conversion starts from left and length of the individual
        index, distance = self.find_cleft_individual(h, self.tracklength)
        y = self.segments[index]
        # generate tracklength
        k = y.left + math.floor(
            1.0
            + math.log(1.0 - random.random() * (1.0 - (self.pc) ** (distance - 1)))
            / self.lnpc
        )
        k_mass = self.recomb_map.position_to_mass(k)
        while y.right <= k:
            y = y.next
        x = y.prev
        if y.left < k:
            # Make new segment
            z = self.alloc_segment(
                k, y.right, k_mass, y.right_mass, y.node, y.population, None, y.next
            )
            if y.next is not None:
                y.next.prev = z
            y.next = None
            y.right = k
            y.right_mass = k_mass
            self.subtract_segment_mass(y, z)
            lhs_tail = y
        else:
            # split the link between x and y.
            x.next = None
            y.prev = None
            z = y
            lhs_tail = x
        z.label = label
        self.set_single_segment_mass(z)
        self.P[z.population].add(z, label)
        # TODO check what needs to be added for full arg
        ret = None
        if return_heads:
            x = lhs_tail
            # Seek back to the head of the x chain
            while x.prev is not None:
                x = x.prev
            ret = x, z
        return ret

    def set_labels(self, segment, new_label):
        while segment is not None:
            recomb_mass = self.L[segment.label].get_frequency(segment.index)
            self.L[new_label].set_value(segment.index, recomb_mass)
            segment.label = new_label
            segment = segment.next

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
        k = self.recomb_map.sample_poisson(start)
        assert k > start
        return k

    def dtwf_recombine(self, x):
        """
        Chooses breakpoints and returns segments sorted by inheritance
        direction, by iterating through segment chain starting with x
        """
        u = self.alloc_segment(-1, -1, -1, -1, -1, -1, None, None)
        v = self.alloc_segment(-1, -1, -1, -1, -1, -1, None, None)
        seg_tails = [u, v]

        # TODO Should this be the recombination rate going foward from x.left?
        if self.recomb_map.total_recombination_rate > 0:
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
                k_mass = self.recomb_map.position_to_mass(k)
                self.num_re_events += 1
                ix = (ix + 1) % 2
                # Make new segment
                if seg_tails[ix] is u or seg_tails[ix] is v:
                    tail = None
                else:
                    tail = seg_tails[ix]
                z = self.alloc_segment(
                    k, x.right, k_mass, x.right_mass, x.node, x.population, tail, x.next
                )
                if z.prev is None:
                    self.set_single_segment_mass(z)
                else:
                    self.set_segment_mass(z, z.prev)
                if x.next is not None:
                    x.next.prev = z
                seg_tails[ix].next = z
                seg_tails[ix] = z
                x.next = None
                x.right = k
                x.right_mass = k_mass
                self.subtract_segment_mass(x, z)
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
                if y.prev is None:
                    self.set_single_segment_mass(y)
                else:
                    self.set_segment_mass(y, y.prev)
                seg_tails[ix] = y
                x = y
            else:
                # No recombination between x.right and y.left
                x = y

        # Remove sentinal segments - this can be handled more simply
        # with pointers in C implemetation
        s = u
        u = s.next
        self.free_segment(s)

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
                s += f"({u.left}-{u.right}->{u.node}({u.index}))"
                u = u.next
            print(s)

    def census_event(self, time):
        for pop in self.P:
            for ancestor in pop.iter_ancestors():
                seg = ancestor
                self.flush_edges()
                u = self.tables.nodes.add_row(
                    time=time, flags=msprime.NODE_IS_CEN_EVENT, population=pop._id
                )
                while seg is not None:
                    # Add an edge joining the segment to the new node.
                    self.store_edge(seg.left, seg.right, u, seg.node)
                    seg.node = u
                    seg = seg.next

    def bottleneck_event(self, pop_id, label, intensity):
        # self.print_state()
        # Merge some of the ancestors.
        pop = self.P[label][pop_id]
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
            r_max = self.m
            while len(H) > 0 and H[0][0] == left:
                x = heapq.heappop(H)[1]
                X.append(x)
                r_max = min(r_max, x.right)
            if len(H) > 0:
                r_max = min(r_max, H[0][0])
            if len(X) == 1:
                x = X[0]
                if len(H) > 0 and H[0][0] < x.right:
                    next_l_mass = H[0][1].left_mass
                    alpha = self.alloc_segment(
                        x.left, H[0][0], x.left_mass, next_l_mass, x.node, x.population
                    )
                    alpha.label = label
                    x.left = H[0][0]
                    x.left_mass = next_l_mass
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
                    alpha = self.alloc_segment(
                        left,
                        right,
                        self.recomb_map.position_to_mass(left),
                        self.recomb_map.position_to_mass(right),
                        u,
                        pop_id,
                    )
                # Update the heaps and make the record.
                for x in X:
                    self.store_edge(left, right, u, x.node)
                    if x.right == right:
                        self.free_segment(x)
                        if x.next is not None:
                            y = x.next
                            heapq.heappush(H, (y.left, y))
                    elif x.right > right:
                        self.set_segment_left_endpoint(x, right)
                        heapq.heappush(H, (x.left, x))

            # loop tail; update alpha and integrate it into the state.
            if alpha is not None:
                if z is None:
                    self.set_single_segment_mass(alpha)
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
                    self.set_segment_mass(alpha, z)
                alpha.prev = z
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
                x.right_mass = y.right_mass
                x.next = y.next
                if y.next is not None:
                    y.next.prev = x
                self.add_segment_mass(x, y)
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
                    alpha = self.alloc_segment(
                        x.left, y.left, x.left_mass, y.left_mass, x.node, x.population
                    )
                    x.left = y.left
                    x.left_mass = y.left_mass
                    alpha.label = x.label
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
                            left,
                            right,
                            self.recomb_map.position_to_mass(left),
                            self.recomb_map.position_to_mass(right),
                            u,
                            population_index,
                        )
                        alpha.label = label
                    self.store_edge(left, right, u, x.node)
                    self.store_edge(left, right, u, y.node)
                    # Now trim the ends of x and y to the right sizes.
                    if x.right == right:
                        self.free_segment(x)
                        x = x.next
                    else:
                        self.set_segment_left_endpoint(x, right)
                    if y.right == right:
                        self.free_segment(y)
                        y = y.next
                    else:
                        self.set_segment_left_endpoint(y, right)

            # loop tail; update alpha and integrate it into the state.
            if alpha is not None:
                if z is None:
                    pop.add(alpha, label)
                    self.set_single_segment_mass(alpha)
                else:
                    if self.full_arg:
                        defrag_required |= z.right == alpha.left
                    else:
                        defrag_required |= (
                            z.right == alpha.left and z.node == alpha.node
                        )
                    z.next = alpha
                    self.set_segment_mass(alpha, z)
                alpha.prev = z
                z = alpha

        if self.full_arg:
            if not coalescence:
                self.store_node(population_index, flags=msprime.NODE_IS_CA_EVENT)
            self.store_arg_edges(z)

        if defrag_required:
            self.defrag_segment_chain(z)
        if coalescence:
            self.defrag_breakpoints()

    def print_state(self):
        print("State @ time ", self.t)
        for l in range(self.num_labels):
            print("Links = ", self.L[l].get_total())
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
        for l in range(self.num_labels):
            print("Fenwick tree[%d]: %d" % (l, self.L[l].get_total()))
            for j in range(1, self.max_segments + 1):
                s = self.L[l].get_frequency(j)
                if s != 0:
                    print("\t", j, "->", s, self.L[l].get_cumulative_frequency(j))
        print("nodes")
        print(self.tables.nodes)
        print("edges")
        print(self.tables.edges)
        self.verify()

    def verify_overlaps(self):
        overlap_counter = OverlapCounter(self.m)
        for pop in self.P:
            for l in range(self.num_labels):
                for u in pop.iter_label(l):
                    while u is not None:
                        overlap_counter.increment_interval(u.left, u.right)
                        u = u.next

        for pos, count in self.S.items():
            if pos != self.m:
                assert count == overlap_counter.overlaps_at(pos)

    def verify(self):
        """
        Checks that the state of the simulator is consistent.
        """
        self.verify_overlaps()
        q = 0
        for l in range(self.num_labels):
            total_mass = 0
            alt_total_mass = 0
            for pop_index, pop in enumerate(self.P):
                for u in pop.iter_label(l):
                    assert u.prev is None
                    left = u.left
                    while u is not None:
                        assert u.population == pop_index
                        assert u.left < u.right
                        l_mass = self.recomb_map.position_to_mass(u.left)
                        r_mass = self.recomb_map.position_to_mass(u.right)
                        assert math.isclose(u.left_mass, l_mass, abs_tol=1e-6)
                        assert math.isclose(u.right_mass, r_mass, abs_tol=1e-6)
                        if u.prev is not None:
                            s = self.recomb_map.mass_between(u.prev.right, u.right)
                            assert u.prev.label == u.label
                        else:
                            s = self.recomb_map.mass_between_left_exclusive(
                                u.left, u.right
                            )
                        right = u.right
                        if self.model != "wf_ped":
                            freq = self.L[l].get_frequency(u.index)
                            total_mass += freq
                            assert math.isclose(s, freq, abs_tol=1e-6)
                        v = u.next
                        if v is not None:
                            assert v.prev == u
                            if u.right > v.left:
                                print("ERROR", u, v)
                            assert u.right <= v.left
                        u = v
                    s = self.recomb_map.mass_between_left_exclusive(left, right)
                    q += s
                    alt_total_mass += s
            assert math.isclose(total_mass, self.L[l].get_total(), abs_tol=1e-6)
            assert math.isclose(total_mass, alt_total_mass, abs_tol=1e-6)
        # add check for dealing with labels
        lab_tot = 0
        for l in range(self.num_labels):
            lab_tot += self.L[l].get_total()
        if self.model != "wf_ped":
            assert math.isclose(q, lab_tot, abs_tol=1e-6)

        assert self.S[self.m] == -1
        # Check the ancestry tracking.
        A = bintrees.AVLTree()
        A[0] = 0
        A[self.m] = -1
        for pop in self.P:
            for l in range(self.num_labels):
                for u in pop.iter_label(l):
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
    if rho == 0:
        gamma = args.gene_conversion_rate[0]
    else:
        gamma = args.gene_conversion_rate[0] * rho
    mean_tracklength = args.gene_conversion_rate[1]
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
        recombination_map = RecombinationMap([0, m], [rho, 0], True)
    else:
        positions = args.recomb_positions
        rates = args.recomb_rates
        recombination_map = RecombinationMap(positions, rates, True)
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
    s = Simulator(
        n,
        m,
        rho,
        recombination_map,
        migration_matrix,
        sample_configuration,
        population_growth_rates,
        population_sizes,
        args.population_growth_rate_change,
        args.population_size_change,
        args.migration_matrix_element_change,
        args.bottleneck,
        args.census_time,
        args.model,
        from_ts=args.from_ts,
        max_segments=100000,
        num_labels=num_labels,
        full_arg=args.full_arg,
        sweep_trajectory=sweep_trajectory,
        time_slice=args.time_slice,
        gene_conversion_rate=gamma,
        gene_conversion_length=mean_tracklength,
        pedigree=pedigree,
    )
    ts = s.simulate()
    ts.dump(args.output_file)
    if args.verbose:
        s.print_state()


def add_simulator_arguments(parser):
    parser.add_argument("sample_size", type=int)
    parser.add_argument("output_file")
    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )
    parser.add_argument("--random-seed", "-s", type=int, default=1)
    parser.add_argument("--num-loci", "-m", type=int, default=100)
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


def main():
    parser = argparse.ArgumentParser()
    add_simulator_arguments(parser)
    args = parser.parse_args()
    run_simulate(args)


if __name__ == "__main__":
    main()
