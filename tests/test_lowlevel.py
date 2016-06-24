# Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
Test cases for the low level C interface to msprime.
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

import collections
import heapq
import itertools
import json
import math
import os.path
import random
import tempfile
import unittest

import tests
import _msprime
# We don't want to import all of the high-level library
# here. All we need is the version.
from msprime import __version__ as _library_version

# We use this to disable hdf5 tests when we're looping through the
# tests to stress test the low-level memory management, etc. For
# obscure reasons we cannot import h5py in the normal way (see
# below for more discussion on this)
enable_h5py_tests = True

# Root node marker
NULL_NODE = -1
NULL_POPULATION = -1


def uniform_recombination_map(sim):
    """
    Returns a uniform recombination map for the specified simulator.
    range 0 to scale.
    """
    return _msprime.RecombinationMap(
        sim.get_num_loci(),
        [0, sim.get_num_loci()],
        [sim.get_scaled_recombination_rate(), 0])


def get_population_configuration(
        sample_size=0, growth_rate=0.0, initial_size=1.0):
    """
    Returns a population configuration dictionary suitable for passing
    to the low-level API.
    """
    return {
        "sample_size": sample_size,
        "growth_rate": growth_rate,
        "initial_size": initial_size
    }


def get_population_parameters_change_event(
        time=0.0, population_id=-1, initial_size=None, growth_rate=None):
    """
    Returns a population change event for the specified values.
    """
    ret = {
        "type": "population_parameters_change",
        "time": time,
        "population_id": population_id
    }
    if initial_size is not None:
        ret["initial_size"] = initial_size
    if growth_rate is not None:
        ret["growth_rate"] = growth_rate
    return ret


def get_size_change_event(time=0.0, size=1.0, population_id=-1):
    """
    Returns a size change demographic event.
    """
    return get_population_parameters_change_event(
        time, population_id, initial_size=size, growth_rate=0)


def get_growth_rate_change_event(time=0.0, growth_rate=1.0, population_id=-1):
    """
    Returns a growth_rate change demographic event.
    """
    return get_population_parameters_change_event(
        time, population_id, growth_rate=growth_rate)


def get_migration_rate_change_event(
        time=0.0, migration_rate=1.0, matrix_index=-1):
    """
    Returns a migration_rate change demographic event.
    """
    return {
        "type": "migration_rate_change",
        "migration_rate": migration_rate,
        "time": time,
        "matrix_index": matrix_index
    }


def get_mass_migration_event(
        time=0.0, source=0, destination=1, proportion=1):
    """
    Returns a mass_migration demographic event.
    """
    return {
        "type": "mass_migration",
        "time": time,
        "source": source,
        "destination": destination,
        "proportion": proportion
    }


def get_migration_matrix(num_populations, value=1.0):
    """
    Returns a simple migration matrix.
    """
    return [
        value * (j != k) for j in range(num_populations)
        for k in range(num_populations)]


def get_random_demographic_events(num_populations, num_events):
    """
    Return some random demographic events for the specified number
    of populations. Note: we return num_events of *each type*.
    """
    events = []
    for j in range(num_events):
        events.append(get_size_change_event(
            time=random.random(), size=random.random(),
            population_id=random.randint(-1, num_populations - 1)))
        events.append(get_growth_rate_change_event(
            time=random.random(), growth_rate=random.random(),
            population_id=random.randint(-1, num_populations - 1)))
        if num_populations > 1:
            matrix_index = 0
            if random.random() < 0.5:
                matrix_index = -1
            else:
                # Don't pick diagonal elements
                while matrix_index % (num_populations + 1) == 0:
                    matrix_index = random.randint(0, num_populations**2 - 1)
            events.append(get_migration_rate_change_event(
                time=random.random(), migration_rate=random.random(),
                matrix_index=matrix_index))
            # Add a mass migration event.
            source = random.randint(0, num_populations - 1)
            destination = source
            while destination == source:
                destination = random.randint(0, num_populations - 1)
            # We use proportion of 0 or 1 so that we can test deterministically
            events.append(get_mass_migration_event(
                time=random.random(), proportion=random.choice([0, 1]),
                source=source, destination=destination))
    sorted_events = sorted(events, key=lambda x: x["time"])
    return sorted_events


def get_leaf_counts(st):
    """
    Returns a list of the leaf node counts for the specfied sparse tree.
    """
    nu = [0 for j in range(st.get_num_nodes())]
    for j in range(st.get_sample_size()):
        u = j
        while u != NULL_NODE:
            nu[u] += 1
            u = st.get_parent(u)
    return nu


def get_tracked_leaf_counts(st, tracked_leaves):
    """
    Returns a list giving the number of leaves in the specified list
    that are in the subtree rooted at each node.
    """
    nu = [0 for j in range(st.get_num_nodes())]
    for j in tracked_leaves:
        # Duplicates not permitted.
        assert nu[j] == 0
        u = j
        while u != NULL_NODE:
            nu[u] += 1
            u = st.get_parent(u)
    return nu


def oriented_forests(n):
    """
    Implementation of Algorithm O from TAOCP section 7.2.1.6.
    Generates all canonical n-node oriented forests.
    """
    p = [k - 1 for k in range(0, n + 1)]
    k = 1
    while k != 0:
        yield p
        if p[n] > 0:
            p[n] = p[p[n]]
            yield p
        k = n
        while k > 0 and p[k] == 0:
            k -= 1
        if k != 0:
            j = p[k]
            d = k - j
            not_done = True
            while not_done:
                if p[k - d] == p[j]:
                    p[k] = p[j]
                else:
                    p[k] = p[k - d] + d
                if k == n:
                    not_done = False
                else:
                    k += 1


def get_mrca(pi, x, y):
    """
    Returns the most recent common ancestor of nodes x and y in the
    oriented forest pi.
    """
    x_parents = [x]
    j = x
    while j != 0:
        j = pi[j]
        x_parents.append(j)
    y_parents = {y: None}
    j = y
    while j != 0:
        j = pi[j]
        y_parents[j] = None
    # We have the complete list of parents for x and y back to root.
    mrca = 0
    j = 0
    while x_parents[j] not in y_parents:
        j += 1
    mrca = x_parents[j]
    return mrca


class TestMRCACalculator(unittest.TestCase):
    """
    Class to test the Schieber-Vishkin algorithm.

    These tests are included here as we use the MRCA calculator below in
    our tests.
    """
    def test_all_oriented_forests(self):
        # Runs through all possible oriented forests and checks all possible
        # node pairs using an inferior algorithm.
        for n in range(2, 9):
            for pi in oriented_forests(n):
                sv = tests.MRCACalculator(pi)
                for j in range(1, n + 1):
                    for k in range(1, j + 1):
                        mrca = get_mrca(pi, j, k)
                        self.assertEqual(mrca, sv.get_mrca(j, k))


class TestModule(tests.MsprimeTestCase):
    """
    Tests for module level stuff.
    """
    def test_library_versions(self):
        major, minor = _msprime.get_gsl_version()
        self.assertIsInstance(major, int)
        self.assertGreater(major, 0)
        self.assertIsInstance(minor, int)
        major, minor, revision = _msprime.get_hdf5_version()
        self.assertIsInstance(major, int)
        self.assertGreater(major, 0)
        self.assertIsInstance(minor, int)
        self.assertIsInstance(revision, int)
        version_str = _msprime.get_library_version_str()
        self.assertEqual(version_str, _library_version)


def get_random_population_models(n):
    """
    Returns n random population models.
    """
    t = 0.0
    models = []
    for j in range(n):
        t += random.uniform(0, 0.3)
        p = random.uniform(0.1, 2.0)
        mod = {"start_time": t}
        if random.random() < 0.5:
            mod["type"] = _msprime.POP_MODEL_EXPONENTIAL
            mod["alpha"] = p
        else:
            mod["type"] = _msprime.POP_MODEL_CONSTANT
            mod["size"] = p
        models.append(mod)
    return models


class LowLevelTestCase(tests.MsprimeTestCase):
    """
    Superclass of tests for the low-level interface.
    """
    def verify_sparse_tree_dict(self, n, pi, tau):
        """
        Verifies that the specified sparse tree in dict format is a
        consistent coalescent history for a sample of size n.
        """
        self.assertEqual(set(pi.keys()), set(tau.keys()))
        self.assertEqual(len(pi), 2 * n - 1)
        self.assertEqual(len(tau), 2 * n - 1)
        # NULL_NODE should not be a node
        self.assertNotIn(NULL_NODE, pi)
        self.assertNotIn(NULL_NODE, tau)
        # verify the root is equal for all leaves
        root = 0
        while pi[root] != NULL_NODE:
            root = pi[root]
        for j in range(n):
            k = j
            while pi[k] != NULL_NODE:
                k = pi[k]
            self.assertEqual(k, root)
        self.assertIn(root, tau)
        # 0 to n - 1 inclusive should always be nodes
        for j in range(n):
            self.assertIn(j, pi)
            self.assertIn(j, tau)
        num_children = collections.defaultdict(int)
        for j in pi.keys():
            num_children[pi[j]] += 1
        # nodes 0 to n are leaves.
        for j in range(n):
            self.assertNotEqual(pi[j], 0)
            self.assertEqual(tau[j], 0)
            self.assertEqual(num_children[j], 0)
        # All non-leaf nodes should be binary with non-zero times.
        taup = {}
        for j in pi.keys():
            if j > n:
                self.assertEqual(num_children[j], 2)
                self.assertGreater(tau[j], 0.0)
                taup[j] = tau[j]
        # times of non leaves should be distinct
        self.assertEqual(len(set(taup)), len(taup))
        # Times of leaves should be zero, and increasing up the tree
        for j in range(n):
            self.assertEqual(tau[j], 0.0)
            last_time = -1
            k = j
            while k in pi:
                self.assertNotEqual(k, pi[k])
                self.assertGreater(tau[k], last_time)
                last_time = tau[k]
                k = pi[k]

    def get_tree_sequence(
            self, sample_size=10, num_loci=100, mutation_rate=10,
            random_seed=1):
        rho = 1.0
        rng = _msprime.RandomGenerator(random_seed)
        sim = _msprime.Simulator(
            sample_size, rng, num_loci=num_loci,
            scaled_recombination_rate=rho)
        sim.run()
        ts = _msprime.TreeSequence()
        recomb_map = uniform_recombination_map(sim)
        ts.create(sim, recomb_map)
        ts.generate_mutations(mutation_rate, rng)
        return ts

    def get_example_tree_sequences(self):
        for n in [2, 3, 100]:
            for m in [1, 2, 100]:
                for mu in [0, 10]:
                    yield self.get_tree_sequence(n, m, mu)

    def verify_iterator(self, iterator):
        """
        Checks that the specified non-empty iterator implements the
        iterator protocol correctly.
        """
        l = list(iterator)
        self.assertGreater(len(l), 0)
        for j in range(10):
            self.assertRaises(StopIteration, next, iterator)


class TestSimulationState(LowLevelTestCase):
    """
    Longer running tests that examine the simulation state in detail for a
    variety of parameters.
    """
    def verify_running_simulation(self, sim):
        """
        Verifies the state of the specified simulation that has run
        for at least one event.
        """
        self.assertRaises(_msprime.LibraryError, sim.debug_demography)
        self.assertGreaterEqual(sim.get_num_breakpoints(), 0)
        self.assertGreater(sim.get_time(), 0.0)
        self.assertGreater(sim.get_num_ancestors(), 1)
        self.assertGreater(sim.get_used_memory(), 0)
        events = sim.get_num_common_ancestor_events()
        events += sim.get_num_recombination_events()
        events += sum(sim.get_num_migration_events())
        self.assertGreater(events, 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
        n = sim.get_sample_size()
        m = sim.get_num_loci()
        N = sim.get_num_populations()
        ancestors = sim.get_ancestors()
        self.assertEqual(len(ancestors), sim.get_num_ancestors())
        for ind in ancestors:
            the_pop_id = ind[0][-1]
            for l, r, node, pop_id in ind:
                self.assertTrue(0 <= l < m)
                self.assertTrue(1 <= r <= m)
                self.assertGreaterEqual(node, 0)
                self.assertTrue(0 <= pop_id < N)
                self.assertEqual(the_pop_id, pop_id)
        breakpoints = [0] + sim.get_breakpoints() + [m]
        self.assertEqual(len(breakpoints), sim.get_num_breakpoints() + 2)
        self.assertEqual(breakpoints, sorted(breakpoints))
        records = sim.get_coalescence_records()
        self.assertEqual(len(records), sim.get_num_coalescence_records())
        for l, r, p, children, t, pop in records:
            self.assertEqual(children, tuple(sorted(children)))
            self.assertTrue(0 <= l < m)
            self.assertTrue(1 <= r <= m)
            self.assertGreater(t, 0.0)
            self.assertTrue(0 <= pop < sim.get_num_populations())
            self.assertIn(l, breakpoints)
            self.assertIn(r, breakpoints)
        # The amount of ancestral material in the coalescence records and
        # the extant segments over all intervals should be either n (if
        # the given interval has not fully coalesced yet) or n - 1 (if
        # full coalescence has occured).
        segments_am = [0 for b in breakpoints[:-1]]
        for ind in ancestors:
            for l, r, _, _ in ind:
                j = breakpoints.index(l)
                while breakpoints[j] < r:
                    segments_am[j] += 1
                    j += 1
        records_am = [0 for b in breakpoints[:-1]]
        for record in records:
            l, r = record[0:2]
            j = breakpoints.index(l)
            while breakpoints[j] < r:
                records_am[j] += 1
                j += 1
        for segment_am, record_am in zip(segments_am, records_am):
            if segment_am == 0:
                self.assertEqual(record_am, n - 1)
            else:
                self.assertEqual(segment_am + record_am, n)

    def verify_trees_equal(self, n, pi, tau, pop, sparse_tree):
        """
        Verifies that the specified maps are equivalent to the specified
        sparse tree object.
        """
        self.assertEqual(n, sparse_tree.get_sample_size())
        pi_p = {}
        tau_p = {}
        pop_p = {}
        for j in range(n):
            u = j
            while u != NULL_NODE and u not in pi_p:
                pop_p[u] = sparse_tree.get_population(u)
                pi_p[u] = sparse_tree.get_parent(u)
                tau_p[u] = sparse_tree.get_time(u)
                u = pi_p[u]
        self.assertEqual(pi_p, pi)
        self.assertEqual(pop_p, pop)
        self.assertEqual(tau_p, tau)
        self.assertEqual(pi_p[sparse_tree.get_root()], NULL_NODE)

    def verify_mrcas(self, sparse_tree):
        """
        Verifies that the MRCAs for the specified tree are computed correctly.
        """
        oriented_forest = [sparse_tree.get_parent(j) for j in range(
            sparse_tree.get_num_nodes())]
        self.assertEqual(
            sum(j > 0 for j in oriented_forest),
            2 * sparse_tree.get_sample_size() - 2)
        mrca_calc = tests.MRCACalculator(oriented_forest)
        nodes = range(sparse_tree.get_num_nodes())
        for u, v in itertools.combinations(nodes, 2):
            self.assertEqual(
                mrca_calc.get_mrca(u, v), sparse_tree.get_mrca(u, v))

    def verify_trees(self, sim, sorted_records):
        """
        Verifies that the specified set of sorted coalescence records
        corresponds to correct trees for the specified simulation.
        """
        ts = _msprime.TreeSequence()
        ts.create(sim, uniform_recombination_map(sim))
        st = _msprime.SparseTree(ts)
        st_iter = _msprime.SparseTreeIterator(ts, st)
        n = sim.get_sample_size()
        pi = {}
        tau = {j: 0 for j in range(n)}
        pops = {j: 0 for j in range(n)}
        # We only support single population here; tests for population
        # assignment are done in the demography tests.
        self.assertEqual(sim.get_num_populations(), 1)
        last_l = 0
        last_t = 0
        num_trees = 0
        live_segments = []
        for l, r, node, children, t, pop in sorted_records:
            if last_l != l:
                last_l = l
                last_t = 0
                for j in range(n):
                    assert j in pi
                # insert the root
                v = 0
                while v in pi:
                    v = pi[v]
                pi[v] = -1
                self.verify_sparse_tree_dict(n, pi, tau)
                # Make sure this is equal to the sparse tree we get from
                # the iterator.
                st = next(st_iter)
                self.verify_trees_equal(n, pi, tau, pops, st)
                self.verify_mrcas(st)
                del pi[v]
                num_trees += 1
            else:
                last_t = t
            heapq.heappush(live_segments, (r, (children, node)))
            while live_segments[0][0] <= l:
                x, (other_children, p) = heapq.heappop(live_segments)
                for c in other_children:
                    del pi[c]
                del tau[p]
                del pops[p]
            pi[children[0]] = node
            pi[children[1]] = node
            tau[node] = t
            pops[node] = pop
            # Ensure that records are sorted by time within a block
            self.assertLessEqual(last_t, t)
        for j in range(n):
            assert j in pi
        # Insert the root branch.
        v = 0
        while v in pi:
            v = pi[v]
        pi[v] = -1
        self.verify_sparse_tree_dict(n, pi, tau)
        st = next(st_iter)
        self.verify_trees_equal(n, pi, tau, pops, st)
        num_trees += 1
        self.assertLessEqual(num_trees, sim.get_num_breakpoints() + 2)
        self.assertRaises(StopIteration, next, st_iter)

    def verify_squashed_records(self, sorted_records):
        """
        Checks to see if there were any unsquashed records in the specified
        set of time sorted records.
        """
        u = sorted_records[0]
        for v in sorted_records[1:]:
            # An unsquashed record would be two adjacent records with the
            # same c1, c2, p and t values.
            self.assertFalse(all(u[j] == v[j] for j in range(1, 5)))
            u = v

    def verify_completed_simulation(self, sim):
        """
        Verifies the state of the specified completed simulation.
        """
        self.assertEqual(sim.get_ancestors(), [])
        self.assertEqual(sim.get_num_ancestors(), 0)
        self.assertGreaterEqual(sim.get_num_breakpoints(), 0)
        self.assertGreater(sim.get_num_coalescence_records(), 0)
        self.assertGreater(sim.get_time(), 0.0)
        events = sim.get_num_common_ancestor_events()
        events += sim.get_num_recombination_events()
        events += sum(sim.get_num_migration_events())
        self.assertGreater(events, 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
        self.assertGreater(sim.get_used_memory(), 0)

        records = sim.get_coalescence_records()
        self.assertGreater(len(records), 0)
        self.assertEqual(len(records), sim.get_num_coalescence_records())
        # Records should be in nondecreasing time order
        times = [t for l, r, node, children, t, pop in records]
        # Children should always be sorted in order.
        for _, _, _, children, _, _ in records:
            self.assertEqual(children, tuple(sorted(children)))
        self.assertEqual(times, sorted(times))
        self.assertEqual(times[-1], sim.get_time())
        self.verify_squashed_records(records)
        left_sorted_records = sorted(records, key=lambda r: (r[0], r[-2]))
        right_sorted_records = sorted(records, key=lambda r: (r[1], -r[-2]))
        self.verify_trees(sim, left_sorted_records)
        # Check the TreeSequence. Ensure we get the records back in the
        # correct orders.
        ts = _msprime.TreeSequence()
        recomb_map = uniform_recombination_map(sim)
        ts.create(sim, recomb_map)
        ts_records = [
            ts.get_record(j, _msprime.MSP_ORDER_TIME)
            for j in range(len(records))]
        self.assertEqual(ts_records, records)
        ts_records = [
            ts.get_record(j, _msprime.MSP_ORDER_LEFT)
            for j in range(len(records))]
        self.assertEqual(ts_records, left_sorted_records)
        ts_records = [
            ts.get_record(j, _msprime.MSP_ORDER_RIGHT)
            for j in range(len(records))]
        self.assertEqual(ts_records, right_sorted_records)
        # Check the simulation parameters
        self.assertEqual(
            json.loads(ts.get_simulation_parameters()),
            json.loads(sim.get_configuration_json())),

    def verify_random_parameters(self):
        mb = 1024 * 1024
        n = random.randint(2, 100)
        m = random.randint(1, 10**6)
        rho = random.uniform(0, 1000)
        N = random.randint(1, 4)
        migration_matrix = [
            random.random() * (j != k) for j in range(N) for k in range(N)]
        population_configuration = [
            get_population_configuration(
                n * (j == 0), random.random(), random.random())
            for j in range(N)]
        demographic_events = get_random_demographic_events(
            N, random.randint(1, 5))
        random_seed = random.randint(0, 2**31)
        max_memory = random.randint(10 * mb, 100 * mb)
        segment_block_size = random.randint(1, 100)
        node_mapping_block_size = random.randint(1, 100)
        avl_node_block_size = random.randint(1, 100)
        coalescence_record_block_size = random.randint(1, 100)
        sim = _msprime.Simulator(
            sample_size=n,
            random_generator=_msprime.RandomGenerator(random_seed),
            num_loci=m,
            scaled_recombination_rate=rho,
            population_configuration=population_configuration,
            demographic_events=demographic_events,
            migration_matrix=migration_matrix,
            max_memory=max_memory,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size,
            coalescence_record_block_size=coalescence_record_block_size)
        for _ in range(3):
            # Check initial state
            self.assertEqual(0, sim.get_num_breakpoints())
            self.assertEqual(0.0, sim.get_time())
            self.assertEqual(n, sim.get_num_ancestors())
            self.assertEqual(0, sim.get_num_common_ancestor_events())
            self.assertEqual(0, sim.get_num_recombination_events())
            self.assertEqual(0, sum(sim.get_num_migration_events()))
            self.assertGreater(sim.get_num_avl_node_blocks(), 0)
            self.assertGreater(sim.get_num_segment_blocks(), 0)
            self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
            self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
            self.assertEqual(sim.get_sample_size(), n)
            self.assertEqual(sim.get_num_loci(), m)
            self.assertEqual(n, len(sim.get_ancestors()))
            self.assertEqual(sim.get_migration_matrix(), migration_matrix)
            self.assertEqual(
                sim.get_population_configuration(), population_configuration)
            # Check the configuration json
            config = {
                "sample_size": n, "num_loci": m,
                "scaled_recombination_rate": rho,
                "migration_matrix": migration_matrix,
                "demographic_events": demographic_events,
                "population_configuration": population_configuration
            }
            self.assertEqual(json.loads(sim.get_configuration_json()), config)
            a = 0
            nodes = set()
            pop_sizes = [0 for _ in population_configuration]
            for ind in sim.get_ancestors():
                self.assertEqual(len(ind), 1)
                l, r, node, pop_id = ind[0]
                pop_sizes[pop_id] += 1
                self.assertEqual(l, 0)
                self.assertEqual(r, m)
                self.assertFalse(node in nodes)
                nodes.add(node)
                a += 1
            for pop_config, pop_size in zip(
                    population_configuration, pop_sizes):
                self.assertEqual(pop_config["sample_size"], pop_size)
            self.assertEqual(a, n)
            for j in range(3):
                # Check the getters to ensure we've got the right values.
                self.assertEqual(n, sim.get_sample_size())
                self.assertEqual(m, sim.get_num_loci())
                self.assertEqual(rho, sim.get_scaled_recombination_rate())
                self.assertEqual(max_memory, sim.get_max_memory())
                self.assertEqual(
                    segment_block_size, sim.get_segment_block_size())
                self.assertEqual(
                    avl_node_block_size, sim.get_avl_node_block_size())
                self.assertEqual(
                    node_mapping_block_size, sim.get_node_mapping_block_size())
                self.assertEqual(
                    coalescence_record_block_size,
                    sim.get_coalescence_record_block_size())
                self.assertEqual(
                    json.loads(sim.get_configuration_json()), config)
                # Run this for a tiny amount of time and check the state
                self.assertFalse(sim.run(1e-8))
                self.verify_running_simulation(sim)
            sim.reset()

    def verify_tree_diffs(self, tree_sequence):
        n = tree_sequence.get_sample_size()
        L = tree_sequence.get_sequence_length()
        # Check some basic properties of the diffs.
        diffs = list(_msprime.TreeDiffIterator(tree_sequence))
        length, records_out, records_in = diffs[0]
        self.assertGreaterEqual(length, 0)
        self.assertEqual(len(records_out), 0)
        self.assertEqual(len(records_in), n - 1)
        # Make sure in records are in increasing time order.
        time_sorted = sorted(records_in, key=lambda x: x[2])
        self.assertEqual(time_sorted, records_in)
        self.assertEqual(sum([l for l, _, _ in diffs]), L)
        for l, records_out, records_in in diffs[1:]:
            self.assertGreaterEqual(l, 0)
            self.assertGreaterEqual(l, 0)
            self.assertEqual(len(records_out), len(records_in))
            for node, children, time in records_out + records_in:
                for c in children:
                    self.assertGreaterEqual(c, 0)
                    self.assertGreater(node, c)
                    self.assertGreater(time, 0.0)
            # Make sure in records are in increasing time order.
            time_sorted = sorted(records_in, key=lambda x: x[2])
            self.assertEqual(time_sorted, records_in)
            # Make sure out records are in decreasing time order.
            time_sorted = sorted(records_out, key=lambda x: -x[2])
            self.assertEqual(time_sorted, records_out)
        # Compare with the Python implementation.
        pts = tests.PythonTreeSequence(tree_sequence)
        python_diffs = list(pts.diffs())
        self.assertGreaterEqual(len(python_diffs), 0)
        self.assertEqual(diffs, python_diffs)

    def verify_leaf_counts(self, tree_sequence):
        st = _msprime.SparseTree(tree_sequence, [])
        for _ in _msprime.SparseTreeIterator(tree_sequence, st):
            self.assertEqual(st.get_count_leaves(), True)
            nu = get_leaf_counts(st)
            nu_prime = [
                st.get_num_leaves(j) for j in
                range(st.get_num_nodes())]
            self.assertEqual(nu, nu_prime)

    def verify_simulation(self, n, m, r):
        """
        Runs the specified simulation and verifies its state.
        """
        # These tests don't work for n == 2
        assert n > 2
        mb = 1024 * 1024
        random_seed = random.randint(0, 2**31)
        sim = _msprime.Simulator(
            sample_size=n, num_loci=m,
            scaled_recombination_rate=r,
            random_generator=_msprime.RandomGenerator(random_seed),
            max_memory=10 * mb, segment_block_size=1000,
            avl_node_block_size=1000, node_mapping_block_size=1000,
            coalescence_record_block_size=1000)
        for _ in range(3):
            # Run the sim for a tiny amount of time and check.
            self.assertFalse(sim.run(1e-8))
            self.verify_running_simulation(sim)
            increment = 0.01
            t = sim.get_time() + increment
            while not sim.run(t):
                self.assertGreaterEqual(sim.get_time(), t)
                self.verify_running_simulation(sim)
                t += increment
            self.verify_completed_simulation(sim)
            # Check the tree sequence.
            tree_sequence = _msprime.TreeSequence()
            recomb_map = uniform_recombination_map(sim)
            tree_sequence.create(sim, recomb_map)
            self.verify_tree_diffs(tree_sequence)
            self.verify_leaf_counts(tree_sequence)
            sim.reset()

    def test_random_sims(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            self.verify_random_parameters()

    def test_small_sims(self):
        self.verify_simulation(3, 1, 0.0)
        self.verify_simulation(3, 100, 0.0)
        self.verify_simulation(3, 10, 1000.0)
        self.verify_simulation(5, 10, 10.0)

    def test_event_by_event(self):
        n = 10
        m = 100
        sim = _msprime.Simulator(
            sample_size=n, num_loci=m, scaled_recombination_rate=1,
            random_generator=_msprime.RandomGenerator(1))
        # We run until time -1 to for initialisation
        sim.run(-1)
        self.assertEqual(sim.get_time(), 0)
        ancestors = sim.get_ancestors()
        self.assertEqual(len(ancestors), n)
        nodes = []
        for ancestor in ancestors:
            self.assertEqual(len(ancestor), 1)
            for l, r, node, pop_id in ancestor:
                self.assertEqual(l, 0)
                self.assertEqual(r, m)
                self.assertEqual(pop_id, 0)
                nodes.append(node)
        self.assertEqual(sorted(nodes), list(range(n)))
        events = 0
        while not sim.run_event():
            events += 1
            total_events = (
                sim.get_num_common_ancestor_events() +
                sim.get_num_recombination_events() +
                sum(sim.get_num_migration_events()))
            self.assertEqual(events, total_events)

    def test_demographic_events(self):
        n = 10
        N = 3
        migration_matrix = [
            random.random() * (j != k) for j in range(N) for k in range(N)]
        population_configuration = [
            get_population_configuration(
                n * (j == 0), random.random(), random.random())
            for j in range(N)]
        demographic_events = get_random_demographic_events(N, 10)
        start_times = [0 for j in range(N)]
        # Rescale time back to very small values so we know that they
        # will definitely happen
        for event in demographic_events:
            event["time"] *= 1e-6
        sim = _msprime.Simulator(
            n, _msprime.RandomGenerator(1), migration_matrix=migration_matrix,
            population_configuration=population_configuration,
            demographic_events=demographic_events)
        # Use a second instance to track the demographic events debugger.
        sim2 = _msprime.Simulator(
            n, _msprime.RandomGenerator(1), migration_matrix=migration_matrix,
            population_configuration=population_configuration,
            demographic_events=demographic_events)
        self.assertEqual(sim.get_migration_matrix(), migration_matrix)
        self.assertEqual(
            sim.get_population_configuration(), population_configuration)

        # Now run the demography debug forward.
        next_event_time = sim2.debug_demography()
        self.assertEqual(sim2.get_migration_matrix(), migration_matrix)
        self.assertEqual(
            sim2.get_population_configuration(), population_configuration)
        # For each event we now run the simulator forward until this time
        # and make sure that the internal state is what it should be.
        for event in demographic_events:
            t = event["time"]
            event_type = event["type"]
            self.assertEqual(next_event_time, t)
            self.assertEqual(sim2.get_migration_matrix(), migration_matrix)
            sim.run(t)
            self.assertEqual(sim.get_time(), t)
            if event_type == "migration_rate_change":
                matrix_index = event["matrix_index"]
                rate = event["migration_rate"]
                if matrix_index == -1:
                    for j in range(N):
                        for k in range(N):
                            if j != k:
                                migration_matrix[j * N + k] = rate
                else:
                    migration_matrix[matrix_index] = rate
            elif event_type == "mass_migration":
                source = event["source"]
                proportion = event["proportion"]
                pop_sizes = [0 for j in range(N)]
                for ind in sim.get_ancestors():
                    _, _, _, pop_id = ind[0]
                    pop_sizes[pop_id] += 1
                if proportion == 1:
                    self.assertEqual(pop_sizes[source], 0)
            else:
                population_id = event["population_id"]
                indexes = [population_id]
                if population_id == -1:
                    indexes = range(N)
                initial_size = event.get("initial_size", None)
                growth_rate = event.get("growth_rate", None)
                for j in indexes:
                    s = population_configuration[j]["initial_size"]
                    a = population_configuration[j]["growth_rate"]
                    size = initial_size
                    if initial_size is None:
                        # Calculate the new size
                        size = s * math.exp(-a * (t - start_times[j]))
                    alpha = a
                    if growth_rate is not None:
                        alpha = growth_rate
                    population_configuration[j]["initial_size"] = size
                    population_configuration[j]["growth_rate"] = alpha
                    start_times[j] = t
            self.assertEqual(sim.get_migration_matrix(), migration_matrix)
            self.assertEqual(
                sim.get_population_configuration(), population_configuration)
            next_event_time = sim2.debug_demography()
        self.assertTrue(math.isinf(next_event_time))


class TestSimulator(LowLevelTestCase):
    """
    Tests for the low-level interface to the simulator.
    """
    def test_bad_parameters(self):
        rng = _msprime.RandomGenerator(1)

        def f(sample_size=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                sample_size, _msprime.RandomGenerator(random_seed), **kwargs)
        # sample_size and random_seed are mandatory
        self.assertRaises(TypeError, _msprime.Simulator)
        self.assertRaises(TypeError, _msprime.Simulator, sample_size=10)
        self.assertRaises(TypeError, _msprime.Simulator, random_generator=rng)
        # check types
        for bad_type in ["1", None, [], {}, int]:
            self.assertRaises(TypeError, f, sample_size=bad_type)
            self.assertRaises(TypeError, f, random_generator=bad_type)
            self.assertRaises(TypeError, f, scaled_recombination_rate=bad_type)
            self.assertRaises(TypeError, f, max_memory=bad_type)
            self.assertRaises(TypeError, f, avl_node_block_size=bad_type)
            self.assertRaises(TypeError, f, segment_block_size=bad_type)
            self.assertRaises(TypeError, f, node_mapping_block_size=bad_type)
            self.assertRaises(
                TypeError, f, coalescence_record_block_size=bad_type)
        # Check for bad values.
        self.assertRaises(_msprime.InputError, f, sample_size=0)
        self.assertRaises(_msprime.InputError, f, sample_size=1)
        self.assertRaises(_msprime.InputError, f, num_loci=0)
        self.assertRaises(_msprime.InputError, f, scaled_recombination_rate=-1)
        self.assertRaises(_msprime.InputError, f, max_memory=0)
        self.assertRaises(_msprime.InputError, f, avl_node_block_size=0)
        self.assertRaises(_msprime.InputError, f, segment_block_size=0)
        self.assertRaises(_msprime.InputError, f, node_mapping_block_size=0)
        self.assertRaises(
            _msprime.InputError, f, coalescence_record_block_size=0)
        # Check for other type specific errors.
        self.assertRaises(OverflowError, f, max_memory=2**65)

    def test_deleting_rng(self):
        rng = _msprime.RandomGenerator(1)
        sim = _msprime.Simulator(10, rng)
        del rng
        sim.run()

    def test_defaults(self):
        n = 10
        sim = _msprime.Simulator(n, _msprime.RandomGenerator(1))
        self.assertEqual(sim.get_migration_matrix(), [0.0])
        self.assertEqual(
            sim.get_population_configuration(),
            [get_population_configuration(n)])
        self.assertEqual(
            json.loads(sim.get_configuration_json())["demographic_events"],
            [])

    def test_bad_population_configurations(self):
        def f(population_configuration):
            return _msprime.Simulator(
                2, _msprime.RandomGenerator(1),
                population_configuration=population_configuration)
        self.assertRaises(TypeError, f, "")
        self.assertRaises(TypeError, f, [""])
        self.assertRaises(ValueError, f, [{}])
        # We must supply all three parameters.
        parameters = ["sample_size", "initial_size", "growth_rate"]
        for k in range(1, 3):
            for t in itertools.combinations(parameters, k):
                d = {k: 2 for k in t}
                self.assertRaises(ValueError, f, [d])
                self.assertRaises(ValueError, f, [d, d])
        self.assertRaises(ValueError, f, [{"start_time": 1, "type": -1000}])
        for bad_number in ["", None, [], {}]:
                self.assertRaises(
                    TypeError, f,
                    [get_population_configuration(sample_size=bad_number)])
                self.assertRaises(
                    TypeError, f,
                    [get_population_configuration(initial_size=bad_number)])
                self.assertRaises(
                    TypeError, f,
                    [get_population_configuration(growth_rate=bad_number)])
        # Cannot have negative sample size or initial_size
        self.assertRaises(
            ValueError, f, [get_population_configuration(sample_size=-1)])
        self.assertRaises(
            ValueError, f, [get_population_configuration(initial_size=-1)])

    def test_bad_sample_configurations(self):
        rng = _msprime.RandomGenerator(1)

        def f(sample_size, pop_sample_sizes):
            population_configuration = [
                get_population_configuration(n) for n in pop_sample_sizes]
            N = len(pop_sample_sizes)
            migration_matrix = [0 for j in range(N) for k in range(N)]
            _msprime.Simulator(
                sample_size, rng,
                population_configuration=population_configuration,
                migration_matrix=migration_matrix)
        for bad_type in [{}, None, 2, [""], [[]], [None]]:
            self.assertRaises(
                TypeError, _msprime.Simulator, 2, rng,
                population_configuration=bad_type)
        # Negative numbers are ValueErrors.
        self.assertRaises(ValueError, f, 2, [-1, 3])
        # Cannot have empty list
        self.assertRaises(ValueError, f, 2, [])
        # Sample sizes that don't add up are input errors.
        self.assertRaises(_msprime.InputError, f, 2, [3])
        self.assertRaises(_msprime.InputError, f, 2, [2, 1])
        self.assertRaises(_msprime.InputError, f, 5, [2, 2])
        # Must provide population_configuration if a migration_matrix
        # is supplied.
        self.assertRaises(
            ValueError, _msprime.Simulator, 2, rng,
            migration_matrix=[0, 0, 0, 0])

    def test_get_population_configurations(self):
        def f(sample_size, conf_tuples):
            population_configuration = [
                get_population_configuration(
                    sample_size=n, initial_size=p, growth_rate=a)
                for n, p, a in conf_tuples]
            N = len(population_configuration)
            migration_matrix = [0 for j in range(N) for k in range(N)]
            s = _msprime.Simulator(
                sample_size, _msprime.RandomGenerator(1),
                population_configuration=population_configuration,
                migration_matrix=migration_matrix)
            conf_dicts = s.get_population_configuration()
            self.assertEqual(len(conf_dicts), len(conf_tuples))
            for conf_dict, conf_tuple in zip(conf_dicts, conf_tuples):
                self.assertEqual(len(conf_dict), 3)
                self.assertEqual(conf_dict["sample_size"], conf_tuple[0])
                self.assertEqual(conf_dict["initial_size"], conf_tuple[1])
                self.assertEqual(conf_dict["growth_rate"], conf_tuple[2])
        f(2, [(2, 1, 1)])
        f(2, [(1, 2, 0), (1, 0.5, 0.1)])
        f(5, [(0, 1, 0.25), (1, 2, 0.5), (2, 3, 0.75), (2, 4, 1)])

    def test_bad_migration_matrix(self):
        def f(num_populations, migration_matrix):
            population_configuration = [
                get_population_configuration()
                for j in range(num_populations)]
            population_configuration[0]["sample_size"] = 2
            return _msprime.Simulator(
                2, _msprime.RandomGenerator(1),
                population_configuration=population_configuration,
                migration_matrix=migration_matrix)
        for bad_type in ["", {}, None, 2, [""], [[]], [None]]:
            self.assertRaises(TypeError, f, 1, bad_type)
        for bad_value in [[1, 2], [-1], [1, 2, 3]]:
            self.assertRaises(ValueError, f, 1, bad_value)

        # Providing the wrong number of populations provokes a ValueError
        self.assertRaises(ValueError, f, 1, [1, 1])
        self.assertRaises(ValueError, f, 2, [1, 1])
        self.assertRaises(ValueError, f, 2, [1, 1, 1])
        self.assertRaises(ValueError, f, 2, [1, 1, 1, 1, 1])
        # Negative values also provoke ValueError
        self.assertRaises(ValueError, f, 2, [0, 1, -1, 0])

        bad_matrices = [
            # Non-zero diagonal gives a InputError
            [[1, 1],
             [1, 1]],
            [[0, 1],
             [1, 1]],
            [[0, 1, 1],
             [1, 0, 1],
             [1, 0, 1]],
        ]
        for matrix in bad_matrices:
            num_populations = len(matrix[0])
            flattened = [v for row in matrix for v in row]
            self.assertRaises(
                _msprime.InputError, f, num_populations, flattened)
        # Must provide migration_matrix when a population config is
        # provided.
        for N in range(1, 5):
            pop_conf = [get_population_configuration(2)] + [
                get_population_configuration(0) for j in range(N - 1)]
            self.assertRaises(
                ValueError, _msprime.Simulator, 2, _msprime.RandomGenerator(1),
                population_configuration=pop_conf)

    def test_get_migration_matrix(self):
        for N in range(1, 10):
            population_configuration = [get_population_configuration(2)] + [
                get_population_configuration(0) for _ in range(N - 1)]
            random_matrix = [
                random.random() * (j != k) for j in range(N)
                for k in range(N)]
            # Deliberately stress the JSON encoding code.
            nasty_matrix = [
                random.random() * 1e9 * (j != k) for j in range(N)
                for k in range(N)]
            matrices = [random_matrix, nasty_matrix]
            for migration_matrix in matrices:
                sim = _msprime.Simulator(
                    2, _msprime.RandomGenerator(1),
                    migration_matrix=migration_matrix,
                    population_configuration=population_configuration)
                self.assertEqual(migration_matrix, sim.get_migration_matrix())
                json_config = json.loads(sim.get_configuration_json())
                json_matrix = json_config["migration_matrix"]
                self.assertEqual(migration_matrix, json_matrix)

    def test_bad_demographic_event_types(self):
        def f(events):
            _msprime.Simulator(
                2, _msprime.RandomGenerator(1), demographic_events=events)
        for bad_type in [None, {}, "", 1]:
            self.assertRaises(TypeError, f, bad_type)
        for bad_type in [None, "", 1, []]:
            self.assertRaises(TypeError, f, [bad_type])
        for bad_type in [None, [], 0]:
            size_change_event = get_size_change_event()
            size_change_event["type"] = bad_type
            self.assertRaises(ValueError, f, [size_change_event])
        for bad_event in [b'', b'1', b'x' * 1000, b'Size_change', 2, 'none']:
            size_change_event = get_size_change_event()
            size_change_event["type"] = bad_event
            self.assertRaises(ValueError, f, [size_change_event])
        for bad_type in [[], "", {}]:
            size_change_event = get_size_change_event(time=bad_type)
            self.assertRaises(TypeError, f, [size_change_event])
            size_change_event = get_size_change_event(population_id=bad_type)
            self.assertRaises(TypeError, f, [size_change_event])
            size_change_event = get_size_change_event(size=bad_type)
            self.assertRaises(TypeError, f, [size_change_event])

    def test_bad_demographic_event_values(self):
        def f(events, num_populations=1):
            population_configuration = [get_population_configuration(2)] + [
                get_population_configuration(0)
                for _ in range(num_populations - 1)]
            _msprime.Simulator(
                2, _msprime.RandomGenerator(1), demographic_events=events,
                population_configuration=population_configuration,
                migration_matrix=get_migration_matrix(num_populations))
        event_generators = [
            get_size_change_event, get_growth_rate_change_event,
            get_migration_rate_change_event,
            get_mass_migration_event]
        for event_generator in event_generators:
            # Negative times not allowed.
            event = event_generator(time=-1)
            self.assertRaises(ValueError, f, [event])
        event_generators = [
            get_size_change_event, get_growth_rate_change_event]
        # Check for bad population ids.
        for event_generator in event_generators:
            for bad_pop_id in [-2, 1, 10**6]:
                event = event_generator(
                    population_id=bad_pop_id)
                self.assertRaises(_msprime.InputError, f, [event])
            for k in range(1, 4):
                event = event_generator(population_id=k)
                self.assertRaises(_msprime.InputError, f, [event], k)
                events = [event_generator(), event]
                self.assertRaises(_msprime.InputError, f, events, k)
        for bad_pop_id in [-2, 1, 10**6]:
            event = get_mass_migration_event(source=bad_pop_id)
            self.assertRaises(_msprime.InputError, f, [event])
            event = get_mass_migration_event(destination=bad_pop_id)
            self.assertRaises(_msprime.InputError, f, [event])
        # Negative size values not allowed
        size_change_event = get_size_change_event(size=-5)
        self.assertRaises(_msprime.InputError, f, [size_change_event])
        # population changes where we specify neither initial_size
        # or growth rate are illegal.
        event = get_population_parameters_change_event()
        self.assertRaises(_msprime.InputError, f, [event])
        # Check for bad matrix indexes
        event_generator = get_migration_rate_change_event
        for bad_index in [-2, 1, 10**6]:
            event = event_generator(matrix_index=bad_index)
            self.assertRaises(_msprime.InputError, f, [event])
        for k in range(1, 4):
            event = event_generator(matrix_index=k * k)
            self.assertRaises(_msprime.InputError, f, [event], k)
            events = [event_generator(), event]
            self.assertRaises(_msprime.InputError, f, events, k)
            # Diagonal matrix values are not accepted.
            for j in range(k):
                event = event_generator(matrix_index=j * k + j)
                self.assertRaises(_msprime.InputError, f, [event], k)
        # Tests specific for mass_migration
        event = get_mass_migration_event(source=0, destination=0)
        self.assertRaises(_msprime.InputError, f, [event])
        for bad_proportion in [-1, 1.1, 1e7]:
            event = get_mass_migration_event(proportion=bad_proportion)
            self.assertRaises(_msprime.InputError, f, [event])

    def test_unsorted_demographic_events(self):
        event_generators = [
            get_size_change_event, get_growth_rate_change_event,
            get_migration_rate_change_event, get_mass_migration_event]
        events = []
        for event_generator in event_generators:
            for _ in range(3):
                events.append(event_generator(time=random.random()))
        sorted_events = sorted(events, key=lambda x: x["time"])
        self.assertNotEqual(events, sorted_events)
        self.assertRaises(
            _msprime.InputError, _msprime.Simulator, 10,
            _msprime.RandomGenerator(1), demographic_events=events)

    def test_get_demographic_events(self):
        # If no events are set, we should get the empty list.
        sim = _msprime.Simulator(2, _msprime.RandomGenerator(1))

        def get_events(s):
            d = json.loads(s.get_configuration_json())
            return d["demographic_events"]

        self.assertEqual([], get_events(sim))
        num_populations = 3
        events = get_random_demographic_events(3, 3)
        population_configuration = [
            get_population_configuration(2),
            get_population_configuration(0),
            get_population_configuration(0)]
        sim = _msprime.Simulator(
            2, _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=[0 for j in range(num_populations**2)],
            demographic_events=events)
        self.assertEqual(events, get_events(sim))

    def test_seed_equality(self):
        simulations = [
            {
                "sample_size": 10,
            }, {
                "sample_size": 10, "num_loci": 100,
                "scaled_recombination_rate": 0.1,
            }, {
                "sample_size": 10,
                "population_configuration": [
                    get_population_configuration(sample_size=3),
                    get_population_configuration(sample_size=3),
                    get_population_configuration(sample_size=4)],
                "migration_matrix": get_migration_matrix(3, 3),
                "demographic_events": [
                    get_size_change_event(0.1, 0.5),
                    get_migration_rate_change_event(0.2, 2),
                    get_growth_rate_change_event(0.3, 2),
                    get_mass_migration_event(0.4, 0, 1, 0.5)]
            }
        ]
        seed = 10
        for params in simulations:
            params["random_generator"] = _msprime.RandomGenerator(seed)
            sim1 = _msprime.Simulator(**params)
            params["random_generator"] = _msprime.RandomGenerator(seed)
            sim2 = _msprime.Simulator(**params)
            sim1.run()
            sim2.run()
            self.assertEqual(
                sim1.get_num_coalescence_records(),
                sim2.get_num_coalescence_records())
            self.assertEqual(
                sim1.get_coalescence_records(),
                sim2.get_coalescence_records())
            self.assertEqual(sim1.get_breakpoints(), sim2.get_breakpoints())

    def test_infinite_waiting_time(self):
        # With no migration we should have an infinite waiting time.
        population_configuration = [
            get_population_configuration(5),
            get_population_configuration(5)]
        sim = _msprime.Simulator(
            10, _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=[0, 0, 0, 0])
        self.assertRaises(_msprime.LibraryError, sim.run)

    def test_simple_event_counters(self):
        for n in [2, 10, 20]:
            sim = _msprime.Simulator(
                n, _msprime.RandomGenerator(1), scaled_recombination_rate=0)
            sim.run()
            self.assertEqual(n - 1, sim.get_num_common_ancestor_events())
            self.assertEqual(0, sim.get_num_recombination_events())
            self.assertEqual([0], sim.get_num_migration_events())

    def test_simple_migration_event_counters(self):
        n = 10
        # No migration at all
        population_configuration = [
            get_population_configuration(n),
            get_population_configuration(0)]
        sim = _msprime.Simulator(
            n, _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=[0.0, 0.0, 0.0, 0.0])
        sim.run()
        self.assertEqual(n - 1, sim.get_num_common_ancestor_events())
        self.assertEqual(0, sim.get_num_recombination_events())
        self.assertEqual([0, 0, 0, 0], sim.get_num_migration_events())
        # Migration between only pops 0 and 1
        matrix = [
            [0, 5, 0],
            [5, 0, 0],
            [0, 0, 0]]
        flattened = [x for row in matrix for x in row]
        population_configuration = [
            get_population_configuration(5),
            get_population_configuration(5),
            get_population_configuration(0)]
        sim = _msprime.Simulator(
            n, _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=flattened)
        sim.run()
        self.assertEqual(n - 1, sim.get_num_common_ancestor_events())
        self.assertEqual(0, sim.get_num_recombination_events())
        migration_events = sim.get_num_migration_events()
        for rate, num_events in zip(flattened, migration_events):
            if rate == 0:
                self.assertEqual(num_events, 0)
            else:
                self.assertGreater(num_events, 0)

    def test_large_migration_matrix_counters(self):
        # Put in linear migration into a larger matrix of populations.
        n = 10
        num_populations = 10
        migration_matrix = [
            [0 for j in range(num_populations)]
            for k in range(num_populations)]
        active_pops = [3, 5, 7]
        for index in range(len(active_pops) - 1):
            j = active_pops[index]
            k = active_pops[index + 1]
            migration_matrix[j][k] = 2
            migration_matrix[k][j] = 2
        population_configuration = [
            get_population_configuration() for _ in range(num_populations)]
        population_configuration[
            active_pops[0]] = get_population_configuration(2)
        population_configuration[
            active_pops[1]] = get_population_configuration(2)
        population_configuration[
            active_pops[2]] = get_population_configuration(6)
        flattened = [x for row in migration_matrix for x in row]
        sim = _msprime.Simulator(
            n, _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=flattened)
        sim.run()
        self.assertEqual(n - 1, sim.get_num_common_ancestor_events())
        self.assertEqual(0, sim.get_num_recombination_events())
        migration_events = sim.get_num_migration_events()
        self.assertEqual(len(migration_events), len(flattened))
        for rate, num_events in zip(flattened, migration_events):
            if rate == 0:
                self.assertEqual(num_events, 0)
            else:
                self.assertGreater(num_events, 0)

    def test_mass_migration(self):
        n = 10
        t = 0.01
        dt = 0.0000001
        sim = _msprime.Simulator(
            n, _msprime.RandomGenerator(1),
            population_configuration=[
                get_population_configuration(n),
                get_population_configuration(0)],
            demographic_events=[
                get_migration_rate_change_event(t),
                get_mass_migration_event(
                    t + dt, source=0, destination=1, proportion=1)],
            migration_matrix=[0, 0, 0, 0])
        sim.run(t)
        pop_sizes_before = [0, 0]
        for ind in sim.get_ancestors():
            for _, _, _, pop_id in ind:
                pop_sizes_before[pop_id] += 1
        sim.run(t + dt)
        pop_sizes_after = [0, 0]
        for ind in sim.get_ancestors():
            for _, _, _, pop_id in ind:
                pop_sizes_after[pop_id] += 1
        self.assertEqual(pop_sizes_before[0], pop_sizes_after[1])

    def test_recombination_event_counters(self):
        n = 10
        # No migration
        population_configuration = [
            get_population_configuration(n),
            get_population_configuration(0)]
        sim = _msprime.Simulator(
            n, _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=[0.0, 0.0, 0.0, 0.0], num_loci=10,
            scaled_recombination_rate=10)
        sim.run()
        self.assertLessEqual(n - 1, sim.get_num_common_ancestor_events())
        self.assertLess(0, sim.get_num_recombination_events())
        self.assertEqual([0, 0, 0, 0], sim.get_num_migration_events())

    def test_single_sink_population_counters(self):
        n = 10
        # Migration only into population 2.
        matrix = [
            [0, 0, 1],
            [0, 0, 1],
            [0, 0, 0]]
        flattened = [x for row in matrix for x in row]
        population_configuration = [
            get_population_configuration(5),
            get_population_configuration(5),
            get_population_configuration(0)]
        sim = _msprime.Simulator(
            n, _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=flattened)
        sim.run()
        self.assertEqual(n - 1, sim.get_num_common_ancestor_events())
        self.assertEqual(0, sim.get_num_recombination_events())
        migration_events = sim.get_num_migration_events()
        for rate, num_events in zip(flattened, migration_events):
            if rate == 0:
                self.assertEqual(num_events, 0)
            else:
                self.assertGreater(num_events, 0)

    def test_reset(self):
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        times = set()
        for _ in range(10):
            sim.run()
            t = sim.get_time()
            self.assertNotIn(t, times)
            times.add(t)
            sim.reset()
            self.assertEqual(sim.get_time(), 0)


class TestTreeSequence(LowLevelTestCase):
    """
    Tests for the low-level interface for the TreeSequence.
    """

    def test_seed_equality(self):
        rng = _msprime.RandomGenerator(10)
        simulations = [
            {
                "sample_size": 10,
                "random_generator": rng,
            }, {
                "sample_size": 10, "num_loci": 100,
                "scaled_recombination_rate": 0.1,
                "random_generator": rng,
            },
        ]
        for params in simulations:
            sim1 = _msprime.Simulator(**params)
            sim2 = _msprime.Simulator(**params)
            sim1.run()
            sim2.run()
            t1 = _msprime.TreeSequence()
            t2 = _msprime.TreeSequence()
            t1.create(sim1, uniform_recombination_map(sim1))
            t2.create(sim1, uniform_recombination_map(sim1))
            self.assertEqual(t1.get_num_records(), t2.get_num_records())
            r1 = [t1.get_record(j) for j in range(t1.get_num_records())]
            r2 = [t2.get_record(j) for j in range(t2.get_num_records())]
            self.assertEqual(r1, r2)
            self.assertEqual(t1.get_mutations(), t2.get_mutations())
            rng1 = _msprime.RandomGenerator(1)
            rng2 = _msprime.RandomGenerator(1)
            t1.generate_mutations(1, rng1)
            t2.generate_mutations(1, rng2)
            self.assertEqual(t1.get_mutations(), t2.get_mutations())

    def test_create_empty_tree_sequence(self):
        sim = _msprime.Simulator(
            sample_size=10, random_generator=_msprime.RandomGenerator(10))
        ts = _msprime.TreeSequence()
        recomb_map = uniform_recombination_map(sim)
        self.assertRaises(ValueError, ts.create, sim, recomb_map)
        sim.run(0.001)
        self.assertRaises(ValueError, ts.create, sim, recomb_map)
        sim.run()
        ts.create(sim, recomb_map)
        self.assertEqual(
            sim.get_num_coalescence_records(), ts.get_num_records())

    def test_file_errors(self):
        ts1 = self.get_tree_sequence()

        def loader(*args):
            ts2 = _msprime.TreeSequence()
            ts2.load(*args)

        for func in [ts1.dump, loader]:
            self.assertRaises(TypeError, func)
            for bad_type in [1, None, [], {}]:
                self.assertRaises(TypeError, func, bad_type)
            # Try to dump/load files we don't have access to or don't exist.
            for f in ["/", "/test.hdf5", "/dir_does_not_exist/x.hdf5"]:
                self.assertRaises(_msprime.LibraryError, func, f)
                try:
                    func(f)
                except _msprime.LibraryError as e:
                    message = str(e)
                    self.assertGreater(len(message), 0)
            # use a long filename and make sure we don't overflow error
            # buffers
            f = "/" + 4000 * "x"
            self.assertRaises(_msprime.LibraryError, func, f)
            try:
                func(f)
            except _msprime.LibraryError as e:
                message = str(e)
                self.assertLess(len(message), 1024)

    def test_initial_state(self):
        # Check the initial state for the dump and load operations.
        ts = _msprime.TreeSequence()
        self.assertRaises(ValueError, ts.dump, "/tmp/file")
        ts = self.get_tree_sequence()
        self.assertRaises(ValueError, ts.load, "/tmp/file")

    def test_not_initialised_after_init_error(self):
        ts = _msprime.TreeSequence()
        self.assertRaises(TypeError, ts.create)
        self.assertRaises(ValueError, ts.dump, "filename")
        self.assertRaises(_msprime.LibraryError, ts.load, "/dev/null")
        self.assertRaises(ValueError, ts.dump, "filename")

    def test_num_nodes(self):
        for ts in self.get_example_tree_sequences():
            max_node = 0
            for j in range(ts.get_num_records()):
                _, _, node, _, _, _ = ts.get_record(j)
                if node > max_node:
                    max_node = node
            self.assertEqual(max_node + 1, ts.get_num_nodes())

    def test_get_population(self):
        for ts in self.get_example_tree_sequences():
            for bad_type in ["1", None, []]:
                self.assertRaises(TypeError, ts.get_population, bad_type)
            self.assertRaises(_msprime.LibraryError, ts.get_population, -1)
            self.assertRaises(
                _msprime.LibraryError, ts.get_population, ts.get_sample_size())
            self.assertRaises(
                _msprime.LibraryError, ts.get_population,
                ts.get_sample_size() + 1)
            for j in range(ts.get_sample_size()):
                # We only check for a single population here. Multi population
                # tests are done in test_demography.
                self.assertEqual(ts.get_population(j), 0)

    def verify_dump_equality(self, ts, outfile):
        """
        Verifies that we can dump a copy of the specified tree sequence
        to the specified file, and load an identical copy.
        """
        ts.dump(outfile.name, skip_h5close=True)
        ts2 = _msprime.TreeSequence()
        ts2.load(outfile.name)
        self.assertEqual(ts.get_sample_size(), ts2.get_sample_size())
        self.assertEqual(ts.get_sequence_length(), ts2.get_sequence_length())
        self.assertEqual(ts.get_num_mutations(), ts2.get_num_mutations())
        self.assertEqual(ts.get_num_nodes(), ts2.get_num_nodes())
        records1 = [
            ts.get_record(j) for j in range(ts.get_num_records())]
        records2 = [
            ts2.get_record(j) for j in range(ts2.get_num_records())]
        self.assertEqual(records1, records2)
        self.assertEqual(ts.get_mutations(), ts2.get_mutations())
        self.assertEqual(
            ts.get_simulation_parameters(), ts2.get_simulation_parameters())
        self.assertEqual(
            ts.get_mutation_parameters(), ts2.get_mutation_parameters())

    def test_dump_equality(self):
        for ts in self.get_example_tree_sequences():
            with tempfile.NamedTemporaryFile() as f:
                self.verify_dump_equality(ts, f)

    def test_generate_mutations_interface(self):
        ts = _msprime.TreeSequence()
        # This hasn't been initialised, so should fail.
        self.assertRaises(ValueError, ts.generate_mutations, ts)
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        sim.run()
        ts.create(sim, uniform_recombination_map(sim))
        self.assertRaises(TypeError, ts.generate_mutations)
        self.assertRaises(TypeError, ts.generate_mutations, mutation_rate=1.0)
        self.assertRaises(TypeError, ts.generate_mutations, random_seed=1.0)
        self.assertRaises(
            TypeError, ts.generate_mutations, mutation_rate=1.0,
            random_seed=10, recombination_map=None)
        self.assertRaises(
            TypeError, ts.generate_mutations, 10, 1, None)
        self.assertRaises(
            TypeError, ts.generate_mutations, mutation_rate=10,
            random_seed=1.0, invalid_param=7)
        self.assertIsNone(ts.get_mutation_parameters())

    def verify_mutations(self, ts):
        mutations = ts.get_mutations()
        self.assertGreater(ts.get_num_mutations(), 0)
        self.assertEqual(len(mutations), ts.get_num_mutations())
        # Check the form of the mutations
        for position, node in mutations:
            self.assertIsInstance(node, int)
            self.assertGreaterEqual(node, 0)
            self.assertLessEqual(node, ts.get_num_nodes())
            self.assertIsInstance(position, float)
            self.assertGreater(position, 0)
            self.assertLess(position, ts.get_sequence_length())
        for j in range(3):
            self.assertEqual(mutations, ts.get_mutations())
        # mutations must be sorted by position order.
        self.assertEqual(mutations, sorted(mutations))

    def test_mutations(self):
        # A mutation rate of 0 should give 0 mutations
        for ts in self.get_example_tree_sequences():
            ts.generate_mutations(0.0, _msprime.RandomGenerator(1))
            self.assertEqual(ts.get_num_mutations(), 0)
            self.assertEqual(len(ts.get_mutations()), 0)
            self.assertIsNone(ts.get_mutation_parameters())
            # A non-zero mutation rate will give more than 0 mutations.
            ts.generate_mutations(10.0, _msprime.RandomGenerator(2))
            self.verify_mutations(ts)
            json_str = ts.get_mutation_parameters()
            self.assertIsNotNone(json_str)
            params = json.loads(json_str)
            self.assertEqual(params["mutation_rate"], 10.0)

    def test_mutation_persistence(self):
        ts = self.get_tree_sequence(mutation_rate=0.0)
        self.assertEqual(ts.get_num_mutations(), 0)
        last_mutations = ts.get_mutations()
        # We should be able to over-write mutations as many times as we like.
        for j in range(10):
            rng = _msprime.RandomGenerator(j + 1)
            ts.generate_mutations(10, rng)
            mutations = ts.get_mutations()
            self.assertNotEqual(mutations, last_mutations)
            last_mutations = mutations

    def test_set_mutations(self):
        ts = self.get_tree_sequence(mutation_rate=0.0)
        for x in [None, "", {}, tuple(), 1]:
            self.assertRaises(TypeError, ts.set_mutations, x)
        invalid_mutations = ["", [1, 1], {1, 1}, None]
        for mutation in invalid_mutations:
            self.assertRaises(TypeError, ts.set_mutations, [mutation])
        invalid_mutations = [tuple(), (1,), (1, 2, 3)]
        for mutation in invalid_mutations:
            self.assertRaises(ValueError, ts.set_mutations, [mutation])
        invalid_mutations = [("1", 0), (0, "1"), (None, 0), ([], 0)]
        for mutation in invalid_mutations:
            self.assertRaises(TypeError, ts.set_mutations, [mutation])
        invalid_mutations = [
            (-1, 1), (ts.get_sequence_length() + 1, 1), (2**32, 1),
            (0, -1), (0, ts.get_num_nodes())]
        for mutation in invalid_mutations:
            self.assertRaises(
                _msprime.LibraryError, ts.set_mutations, [mutation])
        valid_mutations = [[], [(0.1, 1)], [(0.1, 1), (0.2, 2)]]
        for mutations in valid_mutations:
            ts.set_mutations(mutations)
            self.assertEqual(ts.get_mutations(), mutations)
            if len(mutations) == 0:
                self.assertIsNone(ts.get_mutation_parameters())
            else:
                self.assertEqual("{}", ts.get_mutation_parameters())
            # Test dumping the mutations
            with tempfile.NamedTemporaryFile() as f:
                self.verify_dump_equality(ts, f)

    def test_constructor_interface(self):
        tree_sequence = _msprime.TreeSequence()
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        recomb_map = uniform_recombination_map(sim)
        for x in [None, "", {}, [], 1]:
            self.assertRaises(TypeError, tree_sequence.create, x, recomb_map)
            self.assertRaises(TypeError, tree_sequence.create, sim, x)
        # Creating iterators or running method should fail as we
        # haven't initialised it.
        self.assertRaises(ValueError, tree_sequence.get_record, 0)
        self.assertRaises(ValueError, tree_sequence.get_num_records)
        self.assertRaises(ValueError, _msprime.TreeDiffIterator, tree_sequence)
        self.assertRaises(ValueError, _msprime.SparseTree, tree_sequence)
        sparse_tree = _msprime.SparseTree(self.get_tree_sequence(2))
        self.assertRaises(
            ValueError, _msprime.SparseTreeIterator, tree_sequence,
            sparse_tree)
        sim.run()
        tree_sequence.create(sim, recomb_map)
        self.assertRaises(ValueError, tree_sequence.create, sim)
        num_records = sim.get_num_coalescence_records()
        self.assertEqual(num_records, tree_sequence.get_num_records())
        sim_records = sim.get_coalescence_records()
        ts_records = [tree_sequence.get_record(j) for j in range(num_records)]
        self.assertEqual(sim_records, ts_records)

    def verify_get_record_interface(self, ts):
        num_records = ts.get_num_records()
        # We don't accept Python negative indexes here.
        self.assertRaises(IndexError, ts.get_record, -1)
        for j in [0, 10, 10**6]:
            self.assertRaises(IndexError, ts.get_record, num_records + j)
            self.assertRaises(
                IndexError, ts.get_record, num_records + j,
                _msprime.MSP_ORDER_TIME)
            self.assertRaises(
                IndexError, ts.get_record, num_records + j,
                _msprime.MSP_ORDER_LEFT)
            self.assertRaises(
                IndexError, ts.get_record, num_records + j,
                _msprime.MSP_ORDER_RIGHT)
        # Make sure we don't accept bad values for order or index.
        for x in [None, "", {}, []]:
            self.assertRaises(TypeError, ts.get_record, x)
            self.assertRaises(TypeError, ts.get_record, 0, x)
        for x in [-1, 3, 5, 10**6]:
            self.assertRaises(_msprime.LibraryError, ts.get_record, 0, x)

    def test_get_record_interface(self):
        for ts in self.get_example_tree_sequences():
            self.verify_get_record_interface(ts)

    def test_create_interface(self):
        tree_sequence = _msprime.TreeSequence()
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        recomb_map = uniform_recombination_map(sim)
        # First two params are mandatory.
        self.assertRaises(TypeError, tree_sequence.create)
        self.assertRaises(TypeError, tree_sequence.create, sim)
        self.assertRaises(TypeError, tree_sequence.create, recomb_map, sim)
        for bad_type in ["", {}, [], None]:
            self.assertRaises(TypeError, tree_sequence.create, sim, bad_type)
            self.assertRaises(
                TypeError, tree_sequence.create, bad_type, recomb_map)
            # Ne argument is optional, must be a number
            self.assertRaises(
                TypeError, tree_sequence.create, sim, recomb_map, bad_type)

    def test_record_scaling(self):
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        sim.run()
        recomb_map = uniform_recombination_map(sim)
        for Ne in [0.25, 1, 10, 1e6]:
            tree_sequence = _msprime.TreeSequence()
            tree_sequence.create(sim, recomb_map, Ne)
            self.assertEqual(
                sim.get_num_coalescence_records(),
                tree_sequence.get_num_records())
            sim_times = [
                r[-2] for r in sim.get_coalescence_records()]
            for j in range(tree_sequence.get_num_records()):
                generation = tree_sequence.get_record(j)[-2]
                self.assertEqual(generation, sim_times[j] * 4 * Ne)

    def test_pairwise_diversity(self):
        for ts in self.get_example_tree_sequences():
            for bad_type in ["", None, {}]:
                self.assertRaises(
                    TypeError, ts.get_pairwise_diversity, bad_type)
            self.assertRaises(
                ValueError, ts.get_pairwise_diversity, [])
            self.assertRaises(
                ValueError, ts.get_pairwise_diversity, [0])
            self.assertRaises(
                ValueError, ts.get_pairwise_diversity,
                [0, ts.get_sample_size()])
            samples = list(range(ts.get_sample_size()))
            pi1 = ts.get_pairwise_diversity(samples)
            self.assertGreaterEqual(pi1, 0)


class TestHdf5Format(LowLevelTestCase):
    """
    Tests on the HDF5 file format.
    """

    def verify_mutation_parameters_json(self, json_str):
        parameters = json.loads(json_str.decode())
        self.assertIn("mutation_rate", parameters)
        self.assertIsInstance(parameters["mutation_rate"], int, float)

    def verify_tree_parameters_json(self, json_str):
        parameters = json.loads(json_str.decode())
        self.assertIn("scaled_recombination_rate", parameters)
        self.assertIn("sample_size", parameters)
        self.assertIn("num_loci", parameters)
        self.assertIn("population_configuration", parameters)
        self.assertIn("migration_matrix", parameters)
        self.assertIn("demographic_events", parameters)
        self.assertIsInstance(
            parameters["scaled_recombination_rate"], int, float)
        self.assertIsInstance(parameters["sample_size"], int)
        self.assertIsInstance(parameters["num_loci"], int)
        self.assertIsInstance(parameters["population_configuration"], list)
        self.assertIsInstance(parameters["migration_matrix"], list)
        self.assertIsInstance(parameters["demographic_events"], list)

    def verify_environment_json(self, json_str):
        environment = json.loads(json_str.decode())
        self.assertIn("msprime_version", environment)
        version = environment["msprime_version"]
        self.assertEqual(version, _library_version)
        self.assertIn("hdf5_version", environment)
        version = list(map(int, environment["hdf5_version"].split(".")))
        self.assertEqual(len(version), 3)
        self.assertIn("gsl_version", environment)
        version = list(map(int, environment["gsl_version"].split(".")))
        self.assertEqual(len(version), 2)
        uname_keys = [
            "kernel_name", "kernel_release", "kernel_version",
            "hardware_identifier"]
        for key in uname_keys:
            self.assertIn(key, environment)

    def verify_tree_dump_format(self, ts, outfile):
        uint8 = "uint8"
        uint32 = "uint32"
        float64 = "float64"
        # This is an ugly hack here, but we have to do it to
        # be able to use h5py, as it keeps some global state
        # open, and we nuke this when we call h5close() to clean up.
        import h5py
        ts.dump(outfile.name, skip_h5close=True)
        self.assertTrue(os.path.exists(outfile.name))
        self.assertGreater(os.path.getsize(outfile.name), 0)
        root = h5py.File(outfile.name, "r")
        keys = set(root.keys())
        self.assertLessEqual(keys, set(["mutations", "trees", "samples"]))
        self.assertIn("trees", keys)
        self.assertIn("samples", keys)
        if ts.get_num_mutations() == 0:
            self.assertNotIn("mutations", keys)
        else:
            self.assertIn("mutations", keys)
        # Check the basic root attributes
        format_version = root.attrs['format_version']
        self.assertEqual(format_version[0], 2)
        self.assertEqual(format_version[1], 0)
        self.assertEqual(root.attrs["sample_size"], ts.get_sample_size())
        self.assertEqual(
            root.attrs["sequence_length"], ts.get_sequence_length())
        if ts.get_num_mutations() > 0:
            g = root["mutations"]
            # Check the parameters and environment attributes
            self.verify_mutation_parameters_json(g.attrs["parameters"])
            self.verify_environment_json(g.attrs["environment"])
            fields = [("node", uint32), ("position", float64)]
            self.assertEqual(set(g.keys()), set([name for name, _ in fields]))
            for name, dtype in fields:
                self.assertEqual(len(g[name].shape), 1)
                self.assertEqual(g[name].shape[0], ts.get_num_mutations())
                self.assertEqual(g[name].dtype, dtype)
        g = root["trees"]
        self.verify_tree_parameters_json(g.attrs["parameters"])
        self.verify_environment_json(g.attrs["environment"])
        fields = [
            ("left", float64, 1), ("right", float64, 1),
            ("node", uint32, 1), ("children", uint32, 2),
            ("population", uint8, 1), ("time", float64, 1)]
        self.assertEqual(set(g.keys()), set([name for name, _, _ in fields]))
        for name, dtype, dims in fields:
            self.assertEqual(len(g[name].shape), dims)
            self.assertEqual(g[name].shape[0], ts.get_num_records())
            if dims == 2:
                self.assertEqual(g[name].shape[1], 2)
            self.assertEqual(g[name].dtype, dtype)
        g = root["samples"]
        fields = [("population", uint8)]
        self.assertEqual(set(g.keys()), set([name for name, _, in fields]))
        for name, dtype in fields:
            self.assertEqual(len(g[name].shape), 1)
            self.assertEqual(g[name].shape[0], ts.get_sample_size())
            self.assertEqual(g[name].dtype, dtype)
        root.close()

    def test_dump_format(self):
        if enable_h5py_tests:
            for ts in self.get_example_tree_sequences():
                with tempfile.NamedTemporaryFile() as f:
                    self.verify_tree_dump_format(ts, f)

    # @unittest.skip("Skipping due to weird h5py behaviour")
    def test_load_malformed_hdf5(self):
        if enable_h5py_tests:
            # See above for why we import h5py here.
            import h5py
            ts = _msprime.TreeSequence()
            with tempfile.NamedTemporaryFile() as f:
                hfile = h5py.File(f.name, "w")
                # First try the empty hdf5 file.
                hfile.close()
                self.assertRaises(
                    _msprime.LibraryError, ts.load, f.name, skip_h5close=True)

    @unittest.skip("Skipping due to weird h5py behaviour")
    # This test works when it's run on its own, but fails when other tests
    # are run.
    def test_version_load_error(self):
        if enable_h5py_tests:
            # See above for why we import h5py here.
            import h5py
            ts = list(self.get_example_tree_sequences())[0]
            for bad_version in [(0, 1), (0, 8), (2, 0)]:
                with tempfile.NamedTemporaryFile() as f:
                    ts.dump(f.name, skip_h5close=True)
                    hfile = h5py.File(f.name, "r+")
                    hfile.attrs['format_version'] = bad_version
                    hfile.close()
                    other_ts = _msprime.TreeSequence()
                    self.assertRaises(
                        _msprime.LibraryError, other_ts.load, f.name)

    @unittest.skip("Skipping due to weird h5py behaviour")
    def test_optional_population(self):
        if enable_h5py_tests:
            ts = list(self.get_example_tree_sequences())[-1]
            # See above for why we import h5py here.
            import h5py
            with tempfile.NamedTemporaryFile() as f:
                ts.dump(f.name, skip_h5close=True)
                hfile = h5py.File(f.name, "r+")
                del hfile["trees"]["population"]
                hfile.close()
                del hfile
                other_ts = _msprime.TreeSequence()
                other_ts.load(f.name)
                self.assertEqual(
                    ts.get_num_records(), other_ts.get_num_records())
                for j in range(other_ts.get_num_records()):
                    record = other_ts.get_record(j)
                    pop = record[-1]
                    self.assertEqual(pop, NULL_POPULATION)
                    old_record = ts.get_record(j)
                    self.assertEqual(record[:-1], old_record[:-1])

    def test_optional_samples(self):
        if enable_h5py_tests:
            ts = list(self.get_example_tree_sequences())[-1]
            # See above for why we import h5py here.
            import h5py
            with tempfile.NamedTemporaryFile() as f:
                ts.dump(f.name, skip_h5close=True)
                hfile = h5py.File(f.name, "r+")
                del hfile["samples"]["population"]
                hfile.close()
                del hfile
                other_ts = _msprime.TreeSequence()
                other_ts.load(f.name)
                self.assertEqual(
                    ts.get_num_records(), other_ts.get_num_records())
                for j in range(ts.get_num_records()):
                    record = other_ts.get_record(j)
                    old_record = ts.get_record(j)
                    self.assertEqual(record, old_record)
                self.assertEqual(
                    ts.get_mutations(), other_ts.get_mutations())
                # Check that we get the NULL for leaf nodes.
                st = _msprime.SparseTree(other_ts)
                for _ in _msprime.SparseTreeIterator(other_ts, st):
                    for j in range(ts.get_sample_size()):
                        self.assertEqual(
                            st.get_population(j), NULL_POPULATION)

    def test_load_bad_formats(self):
        # try loading a bunch of files in various formats.
        ts = _msprime.TreeSequence()
        with tempfile.NamedTemporaryFile("wb") as f:
            # First, check the emtpy file.
            self.assertRaises(_msprime.LibraryError, ts.load, f.name)
            # Now some ascii text
            f.write(b"Some ASCII text")
            f.flush()
            self.assertRaises(_msprime.LibraryError, ts.load, f.name)
            f.seek(0)
            # Now write 8k of random bytes
            f.write(os.urandom(8192))
            f.flush()
            self.assertRaises(_msprime.LibraryError, ts.load, f.name)


class TestNewickConverter(LowLevelTestCase):
    """
    Tests for the low-level newick converter.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.NewickConverter)
        self.assertRaises(TypeError, _msprime.NewickConverter, None)
        ts = _msprime.TreeSequence()
        # This hasn't been initialised, so should fail.
        self.assertRaises(ValueError, _msprime.NewickConverter, ts)
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        sim.run()
        ts.create(sim, uniform_recombination_map(sim))
        for bad_type in [None, "", [], {}]:
            self.assertRaises(
                TypeError, _msprime.NewickConverter, ts, precision=bad_type)
            self.assertRaises(
                TypeError, _msprime.NewickConverter, ts, Ne=bad_type)
        before = list(_msprime.NewickConverter(ts))
        self.assertGreater(len(before), 0)
        iterator = _msprime.NewickConverter(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = list(iterator)
        self.assertEqual(before, after)
        # make sure the basic form of the output is correct.
        for length, tree in before:
            self.assertIsInstance(length, float)
            self.assertIsInstance(tree, str)

    def test_iterator(self):
        ts = self.get_tree_sequence()
        ncs = [
            _msprime.NewickConverter(ts), _msprime.NewickConverter(ts, 1),
            _msprime.NewickConverter(ts, 0)]
        for nc in ncs:
            self.verify_iterator(nc)

    def get_times(self, tree):
        """
        Returns the time strings from the specified newick tree.
        """
        ret = []
        current_time = None
        for c in tree:
            if c == ":":
                current_time = ""
            elif c in [",", ")"]:
                ret.append(current_time)
                current_time = None
            elif current_time is not None:
                current_time += c
        return ret

    def test_precision(self):
        ts = self.get_tree_sequence()
        self.assertRaises(ValueError, _msprime.NewickConverter, ts, -1)
        self.assertRaises(ValueError, _msprime.NewickConverter, ts, 17)
        self.assertRaises(ValueError, _msprime.NewickConverter, ts, 100)
        for precision in range(17):
            for l, tree in _msprime.NewickConverter(ts, precision=precision):
                times = self.get_times(tree)
                for t in times:
                    if precision == 0:
                        self.assertNotIn(".", t)
                    else:
                        point = t.find(".")
                        self.assertEqual(precision, len(t) - point - 1)


class TestTreeDiffIterator(LowLevelTestCase):
    """
    Tests for the low-level tree diff iterator.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.TreeDiffIterator)
        self.assertRaises(TypeError, _msprime.TreeDiffIterator, None)
        ts = _msprime.TreeSequence()
        # This hasn't been initialised, so should fail.
        self.assertRaises(ValueError, _msprime.TreeDiffIterator, ts)
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        sim.run()
        ts.create(sim, uniform_recombination_map(sim))
        before = list(_msprime.TreeDiffIterator(ts))
        iterator = _msprime.TreeDiffIterator(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = list(iterator)
        self.assertEqual(before, after)

    def test_iterator(self):
        ts = self.get_tree_sequence()
        self.verify_iterator(_msprime.TreeDiffIterator(ts))


class TestSparseTreeIterator(LowLevelTestCase):
    """
    Tests for the low-level sparse tree iterator.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.SparseTreeIterator)
        self.assertRaises(TypeError, _msprime.SparseTreeIterator, None)
        ts = _msprime.TreeSequence()
        # This hasn't been initialised, so should fail.
        other_ts = self.get_tree_sequence(10)
        tree = _msprime.SparseTree(other_ts)
        self.assertRaises(ValueError, _msprime.SparseTreeIterator, ts, tree)
        sim = _msprime.Simulator(10, _msprime.RandomGenerator(1))
        sim.run()
        ts.create(sim, uniform_recombination_map(sim))
        tree = _msprime.SparseTree(ts)
        n_before = 0
        parents_before = []
        for t in _msprime.SparseTreeIterator(ts, tree):
            n_before += 1
            self.assertIs(t, tree)
            pi = {}
            for j in range(t.get_num_nodes()):
                pi[j] = t.get_parent(j)
            parents_before.append(pi)
        self.assertEqual(n_before, len(list(_msprime.TreeDiffIterator(ts))))
        # If we remove the objects, we should get the same results.
        iterator = _msprime.SparseTreeIterator(ts, tree)
        del tree
        del ts
        n_after = 0
        parents_after = []
        for index, t in enumerate(iterator):
            n_after += 1
            self.assertIsInstance(t, _msprime.SparseTree)
            pi = {}
            for j in range(t.get_num_nodes()):
                pi[j] = t.get_parent(j)
            parents_after.append(pi)
            self.assertEqual(index, t.get_index())
        self.assertEqual(parents_before, parents_after)

    def test_iterator(self):
        ts = self.get_tree_sequence()
        tree = _msprime.SparseTree(ts)
        self.verify_iterator(_msprime.SparseTreeIterator(ts, tree))

    def test_root_bug(self):
        # Reproduce a simulation that provoked a root calculation bug.
        params = {
            "random_generator": _msprime.RandomGenerator(878638576),
            "sample_size": 64,
            "num_loci": 57,
            "scaled_recombination_rate": 0.192184324155680
        }
        sim = _msprime.Simulator(**params)
        sim.run()
        ts = _msprime.TreeSequence()
        ts.create(sim, uniform_recombination_map(sim))
        st = _msprime.SparseTree(ts)
        for st in _msprime.SparseTreeIterator(ts, st):
            root = 0
            while st.get_parent(root) != NULL_NODE:
                root = st.get_parent(root)
            self.assertEqual(st.get_root(), root)


class TestHaplotypeGenerator(LowLevelTestCase):
    """
    Tests for the low-level haplotype generator.
    """
    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.HaplotypeGenerator)
        ts = _msprime.TreeSequence()
        # This hasn't been initialised, so should fail.
        self.assertRaises(ValueError, _msprime.HaplotypeGenerator, ts)
        ts = self.get_tree_sequence(num_loci=10)

        for bad_type in ["", {}, [], None]:
            self.assertRaises(
                TypeError, _msprime.HaplotypeGenerator, ts, bad_type)
        n = ts.get_sample_size()
        hg = _msprime.HaplotypeGenerator(ts)
        before = list(hg.get_haplotype(j) for j in range(n))
        hg = _msprime.HaplotypeGenerator(ts)
        num_mutations = ts.get_num_mutations()
        del ts
        # We should keep a reference to the tree sequence.
        after = list(hg.get_haplotype(j) for j in range(n))
        self.assertEqual(before, after)
        # make sure the basic form of the output is correct.
        for h in before:
            self.assertGreater(len(h), 0)
            self.assertIsInstance(h, str)
            self.assertEqual(len(h), num_mutations)


class TestVariantGenerator(LowLevelTestCase):
    """
    Tests for the low-level variant generator.
    """
    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.VariantGenerator)
        ts = _msprime.TreeSequence()
        # This hasn't been initialised, so should fail.
        self.assertRaises(ValueError, _msprime.VariantGenerator, ts)
        ts = self.get_tree_sequence(num_loci=10)

        for bad_type in ["", {}, [], None]:
            self.assertRaises(
                TypeError, _msprime.VariantGenerator, ts, bad_type)
        vg = _msprime.VariantGenerator(ts)
        before = list(vg)
        vg = _msprime.VariantGenerator(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = list(vg)
        self.assertEqual(before, after)

    def test_form(self):
        ts = self.get_tree_sequence(num_loci=10)
        variants = list(_msprime.VariantGenerator(ts))
        self.assertGreater(len(variants), 0)
        self.assertEqual(len(variants), ts.get_num_mutations())
        positions = []
        for pos, variant in variants:
            positions.append(pos)
            self.assertEqual(len(variant), ts.get_sample_size())
        self.assertEqual(
            positions,
            [pos for pos, _ in ts.get_mutations()])


class TestSparseTree(LowLevelTestCase):
    """
    Tests on the low-level sparse tree interface.
    """
    def test_mutations(self):
        for ts in self.get_example_tree_sequences():
            st = _msprime.SparseTree(ts)
            all_mutations = ts.get_mutations()
            all_tree_mutations = []
            for st in _msprime.SparseTreeIterator(ts, st):
                tree_mutations = st.get_mutations()
                self.assertEqual(st.get_num_mutations(), len(tree_mutations))
                all_tree_mutations.extend(tree_mutations)
                for position, node in tree_mutations:
                    self.assertTrue(st.get_left() <= position < st.get_right())
                    self.assertNotEqual(st.get_parent(node), 0)
            self.assertEqual(all_tree_mutations, all_mutations)

    def test_wrong_size(self):
        ts1 = self.get_tree_sequence(sample_size=10)
        ts2 = self.get_tree_sequence(sample_size=2)
        st1 = _msprime.SparseTree(ts1)
        st2 = _msprime.SparseTree(ts2)
        self.assertRaises(
            _msprime.LibraryError, _msprime.SparseTreeIterator, ts1, st2)
        self.assertRaises(
            _msprime.LibraryError, _msprime.SparseTreeIterator, ts2, st1)
        ts1 = self.get_tree_sequence(sample_size=10, num_loci=100)
        ts2 = self.get_tree_sequence(sample_size=10, num_loci=1)
        st1 = _msprime.SparseTree(ts1)
        st2 = _msprime.SparseTree(ts2)
        self.assertRaises(
            _msprime.LibraryError, _msprime.SparseTreeIterator, ts1, st2)
        self.assertRaises(
            _msprime.LibraryError, _msprime.SparseTreeIterator, ts2, st1)
        ts1 = self.get_tree_sequence(
            sample_size=10, num_loci=10, mutation_rate=0)
        ts2 = self.get_tree_sequence(
            sample_size=10, num_loci=10, mutation_rate=10)
        st1 = _msprime.SparseTree(ts1)
        st2 = _msprime.SparseTree(ts2)
        self.assertRaises(
            _msprime.LibraryError, _msprime.SparseTreeIterator, ts1, st2)
        self.assertRaises(
            _msprime.LibraryError, _msprime.SparseTreeIterator, ts2, st1)

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.SparseTree)
        for bad_type in ["", {}, [], None, 0]:
            self.assertRaises(
                TypeError, _msprime.SparseTree, bad_type)
        ts = self.get_tree_sequence()
        for bad_type in ["", {}, True, 1, None]:
            self.assertRaises(TypeError, _msprime.SparseTree, ts, bad_type)
        for n in range(1, 10):
            ts = self.get_tree_sequence(num_loci=n)
            st = _msprime.SparseTree(ts)
            self.assertEqual(st.get_num_nodes(), ts.get_num_nodes())
            self.assertEqual(st.get_sample_size(), ts.get_sample_size())
            # An uninitialised sparse tree should always be zero.
            self.assertEqual(st.get_root(), 0)
            self.assertEqual(st.get_left(), 0)
            self.assertEqual(st.get_right(), 0)
            for j in range(n):
                self.assertEqual(st.get_parent(j), NULL_NODE)
                self.assertEqual(st.get_population(j), NULL_POPULATION)
                self.assertEqual(st.get_children(j), (NULL_NODE, NULL_NODE))
                self.assertEqual(st.get_time(j), 0)

    def test_bad_tracked_leaves(self):
        ts = self.get_tree_sequence()
        for bad_type in ["", {}, [], None]:
            self.assertRaises(TypeError, _msprime.SparseTree, ts, [bad_type])
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, [1, bad_type])
        for bad_leaf in [ts.get_sample_size(), 10**6, -1e6]:
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts, [bad_leaf])
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts, [1, bad_leaf])
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts,
                [1, bad_leaf, 1])

    def test_count_all_leaves(self):
        ts = self.get_tree_sequence(num_loci=10)
        st = _msprime.SparseTree(ts)
        # Without initialisation we should be 0 leaves for every node
        # that is not a leaf.
        for j in range(st.get_num_nodes()):
            l = 1 if j < st.get_sample_size() else 0
            self.assertEqual(st.get_num_leaves(j), l)
            self.assertEqual(st.get_num_tracked_leaves(j), 0)
        # Now, try this for a tree sequence.
        for st in _msprime.SparseTreeIterator(ts, st):
            nu = get_leaf_counts(st)
            nu_prime = [
                st.get_num_leaves(j) for j in
                range(st.get_num_nodes())]
            self.assertEqual(nu, nu_prime)
            # For tracked leaves, this should be all zeros.
            nu = [
                st.get_num_tracked_leaves(j) for j in
                range(st.get_num_nodes())]
            self.assertEqual(nu, list([0 for _ in nu]))

    def test_count_tracked_leaves(self):
        ts = self.get_tree_sequence(sample_size=5, num_loci=10)
        leaves = [j for j in range(ts.get_sample_size())]
        powerset = itertools.chain.from_iterable(
            itertools.combinations(leaves, r) for r in range(len(leaves) + 1))
        for subset in map(list, powerset):
            # Ordering shouldn't make any different.
            random.shuffle(subset)
            st = _msprime.SparseTree(ts, subset)
            for st in _msprime.SparseTreeIterator(ts, st):
                nu = get_tracked_leaf_counts(st, subset)
                nu_prime = [
                    st.get_num_tracked_leaves(j) for j in
                    range(st.get_num_nodes())]
                self.assertEqual(nu, nu_prime)
        # Passing duplicated values should have no effect.
        leaf = 1
        for j in range(1, 20):
            tracked_leaves = [leaf for _ in range(j)]
            st = _msprime.SparseTree(ts, tracked_leaves)
            for st in _msprime.SparseTreeIterator(ts, st):
                nu = get_tracked_leaf_counts(st, [leaf])
                nu_prime = [
                    st.get_num_tracked_leaves(j) for j in
                    range(st.get_num_nodes())]
                self.assertEqual(nu, nu_prime)

    def test_bounds_checking(self):
        for m in range(1, 10):
            ts = self.get_tree_sequence(num_loci=m)
            n = ts.get_num_nodes()
            st = _msprime.SparseTree(ts)
            for v in [-100, -1, n + 1, n + 100, n * 100]:
                self.assertRaises(ValueError, st.get_parent, v)
                self.assertRaises(ValueError, st.get_children, v)
                self.assertRaises(ValueError, st.get_time, v)
                self.assertRaises(
                    ValueError, _msprime.LeafListIterator, st, v)

    def test_mrca_interface(self):
        for num_loci in range(1, 10):
            for sample_size in range(2, 5):
                ts = self.get_tree_sequence(
                    num_loci=num_loci, sample_size=sample_size)
                num_nodes = ts.get_num_nodes()
                st = _msprime.SparseTree(ts)
                for v in [num_nodes, 10**6, NULL_NODE]:
                    self.assertRaises(ValueError, st.get_mrca, v, v)
                    self.assertRaises(ValueError, st.get_mrca, v, 1)
                    self.assertRaises(ValueError, st.get_mrca, 1, v)
                # All the mrcas for an uninitialised tree should be NULL_NODE
                for u, v in itertools.combinations(range(num_nodes), 2):
                    self.assertEqual(st.get_mrca(u, v), NULL_NODE)

    def test_index(self):
        for num_loci in [1, 100]:
            ts = self.get_tree_sequence(num_loci=num_loci, sample_size=10)
            st = _msprime.SparseTree(ts)
            for index, st in enumerate(_msprime.SparseTreeIterator(ts, st)):
                self.assertEqual(index, st.get_index())

    def test_bad_mutations(self):
        ts = self.get_tree_sequence(2, num_loci=200, random_seed=1)
        # Anything bigger than num_nodes should be caught immediately.
        u = ts.get_num_nodes()
        for bad_node in [u, u + 1, 2 * u, -1]:
            self.assertRaises(
                _msprime.LibraryError, ts.set_mutations, [(0.1, bad_node)])
        # We shouldn't be able to assign mutations to the root node
        st = _msprime.SparseTree(ts)
        for st in _msprime.SparseTreeIterator(ts, st):
            x = st.get_left()
            # For more subtle issues where we put mutations on nodes not in
            # the tree, we have to wait until later to detect it.
            other_ts = self.get_tree_sequence(2, num_loci=200, random_seed=1)
            for u in range(ts.get_num_nodes()):
                if st.get_parent(u) == NULL_NODE:
                    other_ts.set_mutations([(x, u)])
                    self.assertRaises(
                        _msprime.LibraryError, _msprime.HaplotypeGenerator,
                        other_ts)
                    self.assertRaises(
                        _msprime.LibraryError, list,
                        _msprime.VariantGenerator(other_ts))


class TestLeafListIterator(LowLevelTestCase):
    """
    Tests for the low-level leaf list iterator.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.LeafListIterator)
        self.assertRaises(TypeError, _msprime.LeafListIterator, None)
        ts = self.get_tree_sequence()
        tree = _msprime.SparseTree(ts)
        for bad_type in [None, "1", []]:
            self.assertRaises(
                TypeError, _msprime.LeafListIterator, tree, bad_type)
        for bad_node in [-1, tree.get_num_nodes() + 1]:
            self.assertRaises(
                ValueError, _msprime.LeafListIterator, tree, bad_node)
        # Do nasty things...
        iterator = _msprime.LeafListIterator(tree, 1)
        del tree
        del ts
        self.assertEqual(list(iterator), [1])

    def test_iterator(self):
        ts = self.get_tree_sequence()
        tree = _msprime.SparseTree(ts)
        for tree in _msprime.SparseTreeIterator(ts, tree):
            self.verify_iterator(_msprime.LeafListIterator(tree, 1))
            self.verify_iterator(
                _msprime.LeafListIterator(tree, tree.get_root()))

    def test_leaf_list(self):
        ts = self.get_tree_sequence()
        st = _msprime.SparseTree(ts)
        for t in _msprime.SparseTreeIterator(ts, st):
            # All leaf nodes should have themselves.
            for j in range(t.get_sample_size()):
                leaves = list(_msprime.LeafListIterator(t, j))
                self.assertEqual(leaves, [j])
            # All non-tree nodes should have 0
            for j in range(t.get_num_nodes()):
                if t.get_parent(j) == 0 and j != t.get_root():
                    leaves = list(_msprime.LeafListIterator(t, j))
                    self.assertEqual(len(leaves), 0)
            # The root should have all leaves.
            leaves = list(_msprime.LeafListIterator(t, t.get_root()))
            self.assertEqual(
                sorted(leaves), list(range(t.get_sample_size())))


class TestRecombinationMap(LowLevelTestCase):
    """
    Tests for the low-level Recombination Map.
    """
    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.RecombinationMap)
        self.assertRaises(TypeError, _msprime.RecombinationMap, 1, None)
        self.assertRaises(TypeError, _msprime.RecombinationMap, 1, None, None)
        self.assertRaises(TypeError, _msprime.RecombinationMap, 1, {}, {})
        self.assertRaises(TypeError, _msprime.RecombinationMap, "1", [], [])
        self.assertRaises(
            ValueError, _msprime.RecombinationMap, 1, [0, 0.1], [1, 2, 3])
        self.assertRaises(
            ValueError, _msprime.RecombinationMap, 1, [], [0, 0])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap, 1, [], [])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap,
            1, [0, -2], [0, 0])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap,
            1, [1, 0], [0, 0])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap,
            1, [0, 1, 0.5], [0, 0, 0])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap,
            0, [0, 1], [0.1, 0])

    def verify_random_recomb_map(self, size, num_random_checks):
        positions = [0] + sorted(
            random.random() for j in range(size - 2)) + [1]
        rates = [random.uniform(0, 5) for j in range(size - 1)] + [0]
        num_loci = 1000
        rm = _msprime.RecombinationMap(
            num_loci=num_loci, positions=positions, rates=rates)
        for bad_type in [{}, None, "10", []]:
            self.assertRaises(TypeError, rm.genetic_to_physical, bad_type)
        for bad_value in [-1, num_loci + 0.001, 1e7, 2**32]:
            self.assertRaises(ValueError, rm.genetic_to_physical, bad_value)
        self.assertEqual(rm.get_size(), size)
        self.assertEqual(rm.get_positions(), positions)
        self.assertEqual(rm.get_rates(), rates)
        self.assertEqual(rm.genetic_to_physical(0), 0)
        self.assertEqual(rm.genetic_to_physical(num_loci), 1)
        total_rate = rm.get_total_recombination_rate()
        self.assertGreater(total_rate, 0)
        self.assertEqual(
            total_rate / (num_loci - 1),
            rm.get_per_locus_recombination_rate())
        for j in range(num_random_checks):
            x = random.uniform(0, num_loci)
            y = rm.genetic_to_physical(x)
            self.assertTrue(0 <= y <= 1)
            z = rm.physical_to_genetic(y)
            self.assertAlmostEqual(x, z)

    def test_coordinate_conversions(self):
        num_random_checks = 100
        for size in [2, 3, 4, 5, 100]:
            self.verify_random_recomb_map(size, num_random_checks)

    def test_total_rate(self):
        m = 1000
        rm = _msprime.RecombinationMap(m, [0, 10], [0.25, 0])
        self.assertEqual(rm.get_total_recombination_rate(), 2.5)
        rm = _msprime.RecombinationMap(m, [0, 10, 20], [0.25, 0.5, 0])
        self.assertEqual(rm.get_total_recombination_rate(), 7.5)

    def test_zero_rate(self):
        for m in [1, 10, 1000]:
            rm = _msprime.RecombinationMap(m, [0, m], [0, 0])
            self.assertEqual(rm.get_total_recombination_rate(), 0)
            self.assertEqual(rm.genetic_to_physical(0), 0)
            self.assertEqual(rm.get_per_locus_recombination_rate(), 0)
            self.assertEqual(rm.get_size(), 2)
            for j in range(m + 1):
                self.assertEqual(rm.genetic_to_physical(j), j)

    def test_uniform_rate(self):
        for m in [1, 10, 100]:
            rm = _msprime.RecombinationMap(m, [0, m], [0.001, 0])
            # For a uniform map with an equal number of loci and physical
            # length we should map 1-1 exactly.
            for v in range(0, m + 1):
                self.assertEqual(v, rm.genetic_to_physical(v))
                # We don't bother special casing for physical to genetic as
                # this isn't really used.
                self.assertAlmostEqual(v, rm.physical_to_genetic(v))


class TestRandomGenerator(unittest.TestCase):
    """
    Tests for the random generator class.
    """
    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.RandomGenerator)
        for bad_type in ["x", 1.0, {}]:
            self.assertRaises(TypeError, _msprime.RandomGenerator, bad_type)
        for bad_value in [-1, 0, 2**32]:
            self.assertRaises(ValueError, _msprime.RandomGenerator, bad_value)

    def test_seed(self):
        for s in [1, 10, 2**32 - 1]:
            rng = _msprime.RandomGenerator(s)
            self.assertEqual(rng.get_seed(), s)


class TestDemographyDebugger(unittest.TestCase):
    """
    Tests for the demography debugging interface.
    """
    def get_simulator(self, events):
        return _msprime.Simulator(
            2, _msprime.RandomGenerator(1), demographic_events=events)

    def test_zero_events(self):
        sim = self.get_simulator([])
        self.assertTrue(math.isinf(sim.debug_demography()))

    def test_state_machine_errors(self):
        sim = self.get_simulator([])
        # It's an error to call debug_demography after run()
        sim.run(1e-9)
        self.assertRaises(_msprime.LibraryError, sim.debug_demography)
        self.assertRaises(_msprime.LibraryError, sim.debug_demography)
        sim.run()
        # It's an error run after debug_demography
        sim = self.get_simulator([])
        self.assertTrue(math.isinf(sim.debug_demography()))
        self.assertRaises(_msprime.LibraryError, sim.run)
        self.assertRaises(_msprime.LibraryError, sim.run)
        self.assertTrue(math.isinf(sim.debug_demography()))
