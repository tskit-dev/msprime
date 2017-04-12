# Copyright (C) 2015-2016 University of Oxford
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

import array
import collections
import heapq
import itertools
import math
import os
import platform
import random
import sys
import tempfile
import unittest

import tests
import _msprime

# We don't want to import all of the high-level library
# here. All we need is the version.
from msprime import __version__ as _library_version

# Root node marker
NULL_NODE = -1
NULL_POPULATION = -1

IS_PY2 = sys.version_info[0] < 3
IS_WINDOWS = platform.system() == "Windows"


def uniform_recombination_map(sim):
    """
    Returns a uniform recombination map for the specified simulator.
    range 0 to scale.
    """
    return _msprime.RecombinationMap(
        sim.get_num_loci(),
        [0, sim.get_num_loci()],
        [sim.get_scaled_recombination_rate(), 0])


def get_simulation_model(name="hudson", **kwargs):
    """
    Returns simulation model dictionary suitable for passing to the low-level API.
    """
    d = {"name": name}
    d.update(kwargs)
    return d


def get_population_configuration(growth_rate=0.0, initial_size=1.0):
    """
    Returns a population configuration dictionary suitable for passing
    to the low-level API.
    """
    return {
        "growth_rate": growth_rate,
        "initial_size": initial_size
    }


def get_samples(sample_size):
    """
    Returns a sample list for the specified size.
    """
    t = (0, 0)
    return [t for _ in range(sample_size)]


def get_population_samples(*args):
    """
    Returns a sample list for the specified list of subpopulation
    sample sizes.
    """
    ret = []
    for j, n in enumerate(args):
        t = (j, 0)
        for _ in range(n):
            ret.append(t)
    return ret


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


def get_simple_bottleneck_event(
        time=0.0, population_id=0, proportion=1):
    """
    Returns a simple bottleneck demographic event.
    """
    return {
        "type": "simple_bottleneck",
        "time": time,
        "population_id": population_id,
        "proportion": proportion
    }


def get_instantaneous_bottleneck_event(
        time=0.0, population_id=0, strength=10):
    """
    Returns a instantaneous bottleneck demographic event.
    """
    return {
        "type": "instantaneous_bottleneck",
        "time": time,
        "population_id": population_id,
        "strength": strength
    }


def get_migration_matrix(num_populations, value=1.0):
    """
    Returns a simple migration matrix.
    """
    return [
        value * (j != k) for j in range(num_populations)
        for k in range(num_populations)]


def get_example_simulator(
        sample_size=10, random_seed=1, num_populations=1,
        store_migrations=False):
    samples = [(j % num_populations, 0) for j in range(sample_size)]
    migration_matrix = [1 for _ in range(num_populations**2)]
    for j in range(num_populations):
        migration_matrix[j * num_populations + j] = 0
    population_configuration = [
        get_population_configuration() for j in range(num_populations)]
    sim = _msprime.Simulator(
        samples, _msprime.RandomGenerator(random_seed),
        population_configuration=population_configuration,
        migration_matrix=migration_matrix,
        store_migrations=store_migrations)
    return sim


def populate_tree_sequence(
        sim, Ne=0.25, mutation_generator=None, provenance_strings=[]):
    nodes = _msprime.NodeTable()
    edgesets = _msprime.EdgesetTable()
    migrations = _msprime.MigrationTable()
    sites = _msprime.SiteTable()
    mutations = _msprime.MutationTable()
    ts = _msprime.TreeSequence()
    sim.populate_tables(nodes, edgesets, migrations, Ne=Ne)
    if mutation_generator is not None:
        mutation_generator.generate(nodes, edgesets, sites, mutations)
    ts.load_tables(
        nodes, edgesets, migrations, sites, mutations, provenance_strings)
    return ts


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
            # Add some bottlenecks
            events.append(get_simple_bottleneck_event(
                time=random.random(), proportion=random.uniform(0, 0.25),
                population_id=random.randint(0, num_populations - 1)))
            events.append(get_instantaneous_bottleneck_event(
                time=random.random(), strength=random.uniform(0, 0.01),
                population_id=random.randint(0, num_populations - 1)))

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
        self.assertLessEqual(len(pi), 2 * n - 1)
        self.assertEqual(len(tau), len(pi))
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
                self.assertGreaterEqual(num_children[j], 2)
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
            random_seed=1, demographic_events=[], num_provenance_strings=5):
        rho = 1.0
        rng = _msprime.RandomGenerator(random_seed)
        sim = _msprime.Simulator(
            get_samples(sample_size), rng, num_loci=num_loci,
            scaled_recombination_rate=rho,
            demographic_events=demographic_events)
        sim.run()
        nodes = _msprime.NodeTable()
        edgesets = _msprime.EdgesetTable()
        migrations = _msprime.MigrationTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        ts = _msprime.TreeSequence()
        mutgen = _msprime.MutationGenerator(rng, mutation_rate)
        provenance_strings = [
            b"xxxxxxx" * (j + 1) for j in range(num_provenance_strings)]
        sim.populate_tables(nodes, edgesets, migrations)
        mutgen.generate(nodes, edgesets, sites, mutations)
        ts.load_tables(
            nodes, edgesets, migrations, sites, mutations,
            provenance_strings=provenance_strings)
        self.assertEqual(ts.get_num_nodes(), nodes.num_rows)
        self.assertEqual(ts.get_num_mutations(), mutations.num_rows)
        self.assertEqual(ts.get_num_edgesets(), edgesets.num_rows)
        return ts

    def get_nonbinary_tree_sequence(self):
        bottlenecks = [
            get_simple_bottleneck_event(0.1, 0, 0.1),
            get_instantaneous_bottleneck_event(0.11, 0, 5)]
        return self.get_tree_sequence(demographic_events=bottlenecks)

    def get_example_tree_sequences(self):
        bottlenecks = [
            get_simple_bottleneck_event(0.1, 0, 0.1),
            get_instantaneous_bottleneck_event(0.11, 0, 5)]
        for n in [2, 3, 100]:
            for m in [1, 2, 100]:
                for mu in [0, 10]:
                    yield self.get_tree_sequence(n, m, mu)
                    yield self.get_tree_sequence(
                        n, m, mu, demographic_events=bottlenecks)

    def get_example_migration_tree_sequence(self):
        ts = _msprime.TreeSequence()
        n = 20
        sim = get_example_simulator(n, num_populations=3, store_migrations=True)
        sim.run()
        nodes = _msprime.NodeTable()
        edgesets = _msprime.EdgesetTable()
        migrations = _msprime.MigrationTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        ts = _msprime.TreeSequence()
        sim.populate_tables(nodes, edgesets, migrations)
        mutation_rate = 10
        mutgen = _msprime.MutationGenerator(_msprime.RandomGenerator(1), mutation_rate)
        mutgen.generate(nodes, edgesets, sites, mutations)
        ts.load_tables(nodes, edgesets, migrations, sites, mutations)
        self.assertGreater(sim.get_num_migrations(), 0)
        self.assertGreater(ts.get_num_mutations(), 0)
        return ts

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
        self.assertGreater(sim.get_num_migration_blocks(), 0)
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
        migrations = sim.get_migrations()
        self.assertEqual(len(migrations), sim.get_num_migrations())
        if sim.get_store_migrations():
            self.assertGreaterEqual(
                sim.get_num_migrations(), sum(sim.get_num_migration_events()))

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
        ts = populate_tree_sequence(sim)
        st = _msprime.SparseTree(ts)
        st_iter = _msprime.SparseTreeIterator(st)
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
            for c in children:
                pi[c] = node
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
        self.assertEqual(ts.get_num_trees(), num_trees)
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
        events += sim.get_num_rejected_common_ancestor_events()
        events += sum(sim.get_num_migration_events())
        self.assertGreater(events, 0)
        if sim.get_model() == "hudson":
            self.assertEqual(sim.get_num_rejected_common_ancestor_events(), 0)
        elif sim.get_model() in ["smc", "smc_prime"]:
            self.assertGreaterEqual(sim.get_num_rejected_common_ancestor_events(), 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
        self.assertGreater(sim.get_num_migration_blocks(), 0)
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
        self.verify_trees(sim, left_sorted_records)
        # Check the TreeSequence. Ensure we get the records back in the
        # correct orders.
        ts = populate_tree_sequence(sim)
        ts_records = [ts.get_record(j) for j in range(len(records))]
        self.assertEqual(ts_records, records)

    def verify_random_parameters(self):
        mb = 1024 * 1024
        n = random.randint(2, 100)
        m = random.randint(1, 10**6)
        rho = random.uniform(0, 1000)
        N = random.randint(1, 4)
        store_migrations = random.choice([True, False])
        migration_matrix = [
            random.random() * (j != k) for j in range(N) for k in range(N)]
        population_configuration = [
            get_population_configuration(random.random(), random.random())
            for j in range(N)]
        demographic_events = get_random_demographic_events(
            N, random.randint(1, 5))
        sample_sizes = [0 for j in range(N)]
        sample_sizes[0] = n
        random_seed = random.randint(0, 2**31)
        max_memory = random.randint(10 * mb, 100 * mb)
        segment_block_size = random.randint(1, 100)
        node_mapping_block_size = random.randint(1, 100)
        avl_node_block_size = random.randint(1, 100)
        coalescence_record_block_size = random.randint(1, 100)
        migration_block_size = random.randint(1, 100)
        sim = _msprime.Simulator(
            samples=get_population_samples(*sample_sizes),
            random_generator=_msprime.RandomGenerator(random_seed),
            num_loci=m,
            store_migrations=store_migrations,
            scaled_recombination_rate=rho,
            population_configuration=population_configuration,
            demographic_events=demographic_events,
            migration_matrix=migration_matrix,
            max_memory=max_memory,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size,
            coalescence_record_block_size=coalescence_record_block_size,
            migration_block_size=migration_block_size)
        for _ in range(3):
            # Check initial state
            self.assertEqual(0, sim.get_num_breakpoints())
            self.assertEqual(0.0, sim.get_time())
            self.assertEqual(n, sim.get_num_ancestors())
            self.assertEqual(0, sim.get_num_common_ancestor_events())
            self.assertEqual(0, sim.get_num_rejected_common_ancestor_events())
            self.assertEqual(0, sim.get_num_recombination_events())
            self.assertEqual(0, sum(sim.get_num_migration_events()))
            self.assertGreater(sim.get_num_avl_node_blocks(), 0)
            self.assertGreater(sim.get_num_segment_blocks(), 0)
            self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
            self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
            self.assertGreater(sim.get_num_migration_blocks(), 0)
            self.assertEqual(sim.get_sample_size(), n)
            self.assertEqual(sim.get_num_loci(), m)
            self.assertEqual(n, len(sim.get_ancestors()))
            self.assertEqual(sim.get_migration_matrix(), migration_matrix)
            self.assertEqual(
                sim.get_population_configuration(), population_configuration)
            a = 0
            nodes = set()
            samples = sim.get_samples()
            sample_pop_sizes = [0 for _ in population_configuration]
            for sample in samples:
                sample_pop_sizes[sample[0]] += 1
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
            for n1, n2 in zip(sample_pop_sizes, pop_sizes):
                self.assertEqual(n1, n2)
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
                    migration_block_size, sim.get_migration_block_size())
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
        nodes_in = sum(len(children) - 1 for _, children, _ in records_in)
        self.assertEqual(nodes_in, n - 1)
        # Make sure in records are in increasing time order.
        time_sorted = sorted(records_in, key=lambda x: x[2])
        self.assertEqual(time_sorted, records_in)
        self.assertEqual(sum([l for l, _, _ in diffs]), L)
        for l, records_out, records_in in diffs[1:]:
            self.assertGreaterEqual(l, 0)
            self.assertGreaterEqual(l, 0)
            nodes_out = 0
            for _, children, _ in records_out:
                nodes_out += len(children) - 1
            nodes_in = 0
            for _, children, _ in records_in:
                nodes_in += len(children) - 1
            self.assertEqual(nodes_out, nodes_in)
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
        st = _msprime.SparseTree(
            tree_sequence, flags=_msprime.LEAF_COUNTS, tracked_leaves=[])
        for _ in _msprime.SparseTreeIterator(st):
            self.assertEqual(st.get_flags(), _msprime.LEAF_COUNTS)
            nu = get_leaf_counts(st)
            nu_prime = [
                st.get_num_leaves(j) for j in
                range(st.get_num_nodes())]
            self.assertEqual(nu, nu_prime)

    def verify_simulation(
            self, n, m, r, demographic_events=[], model=get_simulation_model()):
        """
        Runs the specified simulation and verifies its state.
        """
        # These tests don't work for n == 2
        assert n > 2
        mb = 1024 * 1024
        random_seed = random.randint(0, 2**31)
        sim = _msprime.Simulator(
            samples=get_samples(n), num_loci=m,
            scaled_recombination_rate=r,
            random_generator=_msprime.RandomGenerator(random_seed),
            demographic_events=demographic_events,
            max_memory=10 * mb, segment_block_size=1000,
            avl_node_block_size=1000, node_mapping_block_size=1000,
            coalescence_record_block_size=1000, model=model)
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
            tree_sequence = populate_tree_sequence(sim)
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
        self.verify_simulation(
            5, 10, 10.0,
            demographic_events=[
                get_simple_bottleneck_event(time=0.2, proportion=1)])
        self.verify_simulation(3, 10, 1.0, model=get_simulation_model("smc"))
        self.verify_simulation(4, 10, 2.0, model=get_simulation_model("smc_prime"))

    def test_event_by_event(self):
        n = 10
        m = 100
        sim = _msprime.Simulator(
            samples=get_samples(n), num_loci=m, scaled_recombination_rate=1,
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
        random.seed(11)
        n = 10
        N = 3
        migration_matrix = [
            random.random() * (j != k) for j in range(N) for k in range(N)]
        population_configuration = [
            get_population_configuration(random.random(), random.random())
            for j in range(N)]
        demographic_events = get_random_demographic_events(N, 10)
        start_times = [0 for j in range(N)]
        # Rescale time back to very small values so we know that they
        # will definitely happen
        for event in demographic_events:
            event["time"] *= 1e-6
        sim = _msprime.Simulator(
            get_population_samples(n, 0, 0), _msprime.RandomGenerator(1),
            migration_matrix=migration_matrix,
            population_configuration=population_configuration,
            demographic_events=demographic_events)
        # Use a second instance to track the demographic events debugger.
        sim2 = _msprime.Simulator(
            get_samples(n), _msprime.RandomGenerator(1),
            migration_matrix=migration_matrix,
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
            completed = sim.run(t)
            self.assertFalse(completed)
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
            elif event_type in ["simple_bottleneck", "instantaneous_bottleneck"]:
                # Not much we can test for here...
                pass
            else:
                self.assertEqual(event_type, "population_parameters_change")
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
                get_samples(sample_size),
                _msprime.RandomGenerator(random_seed), **kwargs)
        # samples and random_seed are mandatory
        self.assertRaises(TypeError, _msprime.Simulator)
        self.assertRaises(
            TypeError, _msprime.Simulator, get_samples(10))
        self.assertRaises(TypeError, _msprime.Simulator, random_generator=rng)
        # check types
        for bad_type in ["1", None, {}, int]:
            self.assertRaises(TypeError, f, samples=bad_type)
            self.assertRaises(TypeError, f, random_generator=bad_type)
            self.assertRaises(TypeError, f, scaled_recombination_rate=bad_type)
            self.assertRaises(TypeError, f, max_memory=bad_type)
            self.assertRaises(TypeError, f, avl_node_block_size=bad_type)
            self.assertRaises(TypeError, f, segment_block_size=bad_type)
            self.assertRaises(TypeError, f, node_mapping_block_size=bad_type)
            self.assertRaises(
                TypeError, f, coalescence_record_block_size=bad_type)
        # Check for bad values.
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

    def test_non_parametric_simulation_models(self):

        def f(sample_size=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(sample_size),
                _msprime.RandomGenerator(random_seed), **kwargs)
        for bad_type in [0, None, str]:
            self.assertRaises(TypeError, f, model=bad_type)
        for bad_dict in [{}, {"noname": 1}]:
            self.assertRaises(ValueError, f, model=bad_dict)
        for bad_model in ["", "SMC", "ABC", "hud", None, 1234, {}, []]:
            self.assertRaises(ValueError, f, model=get_simulation_model(bad_model))
        for name in ["hudson", "smc", "smc_prime"]:
            model = get_simulation_model(name)
            sim = f(model=model)
            self.assertEqual(sim.get_model(), model)

    def test_dirac_simulation_model(self):

        def f(sample_size=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(sample_size),
                _msprime.RandomGenerator(random_seed), **kwargs)
        for bad_type in [None, str, "sdf"]:
            model = get_simulation_model("dirac", psi=bad_type, c=1.0)
            self.assertRaises(TypeError, f, model=model)
            model = get_simulation_model("dirac", psi=0.5, c=bad_type)
            self.assertRaises(TypeError, f, model=model)
        self.assertRaises(ValueError, f, model=get_simulation_model("dirac"))
        # TODO check for bad values when range checking is done.
        for psi in [0.99, 0.2, 1e-4]:
            for c in [5.0, 1e2, 1e-4]:
                model = get_simulation_model("dirac", psi=psi, c=c)
                sim = f(model=model)
                self.assertEqual(sim.get_model(), model)

    def test_beta_simulation_model(self):

        def f(sample_size=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(sample_size),
                _msprime.RandomGenerator(random_seed), **kwargs)
        for bad_type in [None, str, "sdf"]:
            model = get_simulation_model("beta", alpha=bad_type, truncation_point=1)
            self.assertRaises(TypeError, f, model=model)
            model = get_simulation_model("beta", alpha=1, truncation_point=bad_type)
            self.assertRaises(TypeError, f, model=model)
        model = get_simulation_model("beta", alpha=1)
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model("beta", truncation_point=1)
        self.assertRaises(ValueError, f, model=model)
        # TODO check for bad values when range checking is done.
        for alpha in [1.0, 2.2, 1e-4]:
            for truncation_point in [0.1, 1.5, 1e9]:
                model = get_simulation_model(
                    "beta", alpha=alpha, truncation_point=truncation_point)
                sim = f(model=model)
                self.assertEqual(sim.get_model(), model)

    def test_store_migrations(self):
        def f(sample_size=10, random_seed=1, **kwargs):
            samples = [(j % 2, 0) for j in range(sample_size)]
            migration_matrix = [0, 1, 1, 0]
            population_configuration = [
                get_population_configuration() for j in range(2)]
            return _msprime.Simulator(
                samples, _msprime.RandomGenerator(random_seed),
                population_configuration=population_configuration,
                migration_matrix=migration_matrix, **kwargs)
        for bad_type in [[], "False", None, {}, str]:
            self.assertRaises(TypeError, f, store_migrations=bad_type)
        for store_records in [True, False]:
            sim = f(store_migrations=store_records)
            self.assertEqual(sim.get_store_migrations(), store_records)
            sim.run()
            if store_records:
                self.assertGreater(sim.get_num_migrations(), 0)
            else:
                self.assertEqual(sim.get_num_migrations(), 0)
            records = sim.get_migrations()
            self.assertEqual(len(records), sim.get_num_migrations())
            for l, r, node, source, dest, time in records:
                self.assertTrue(0 <= l < sim.get_num_loci())
                self.assertTrue(0 < r <= sim.get_num_loci())
                self.assertTrue(0 <= source < 2)
                self.assertTrue(0 <= dest < 2)
                self.assertGreaterEqual(node, 0)
                self.assertGreater(time, 0)
        # By default, this should be off.
        sim = f()
        self.assertFalse(sim.get_store_migrations())

    def test_bad_samples(self):
        rng = _msprime.RandomGenerator(1)

        def f(samples):
            return _msprime.Simulator(samples, rng)

        for bad_type in [None, {}, _msprime.Simulator]:
            self.assertRaises(TypeError, f, bad_type)
            self.assertRaises(TypeError, f, [(0, 0), bad_type])
            self.assertRaises(TypeError, f, [(0, 0), (bad_type, 0)])
            self.assertRaises(TypeError, f, [(0, 0), (0, bad_type)])
        self.assertRaises(_msprime.InputError, f, [])
        self.assertRaises(_msprime.InputError, f, [(0, 0)])
        self.assertRaises(ValueError, f, [(0, 0), (0, 0, 0)])
        self.assertRaises(ValueError, f, [(0, 0), (-1, 0)])
        self.assertRaises(ValueError, f, [(0, 0), (0, -1)])
        # Only tuples are supported.
        self.assertRaises(TypeError, f, [(0, 0), [0, 0]])

    def test_get_samples(self):
        N = 4
        samples = [
            (random.randint(0, N - 1), random.random()) for _ in range(10)]
        # There must be at least one sample at the present time.
        samples[-1] = (0, 0)
        rng = _msprime.RandomGenerator(1)
        sim = _msprime.Simulator(
            samples, rng,
            population_configuration=[
                get_population_configuration() for _ in range(N)],
            migration_matrix=[0 for j in range(N * N)])
        self.assertEqual(samples, sim.get_samples())

    def test_deleting_rng(self):
        rng = _msprime.RandomGenerator(1)
        sim = _msprime.Simulator(get_samples(10), rng)
        del rng
        sim.run()

    def test_defaults(self):
        n = 10
        sim = _msprime.Simulator(get_samples(n), _msprime.RandomGenerator(1))
        self.assertEqual(sim.get_migration_matrix(), [0.0])
        self.assertEqual(
            sim.get_population_configuration(),
            [get_population_configuration()])

    def test_bad_population_configurations(self):
        def f(population_configuration):
            return _msprime.Simulator(
                get_samples(2), _msprime.RandomGenerator(1),
                population_configuration=population_configuration)
        self.assertRaises(TypeError, f, "")
        self.assertRaises(TypeError, f, [""])
        self.assertRaises(ValueError, f, [{}])
        # We must supply bot parameters.
        parameters = ["initial_size", "growth_rate"]
        for k in range(1, 3):
            for t in itertools.combinations(parameters, k):
                d = {k: 2 for k in t}
                self.assertRaises(ValueError, f, [d])
                self.assertRaises(ValueError, f, [d, d])
        self.assertRaises(ValueError, f, [{"start_time": 1, "type": -1000}])
        for bad_number in ["", None, [], {}]:
                self.assertRaises(
                    TypeError, f,
                    [get_population_configuration(initial_size=bad_number)])
                self.assertRaises(
                    TypeError, f,
                    [get_population_configuration(growth_rate=bad_number)])
        # Cannot have negative for initial_size
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
                get_samples(sample_size), rng,
                population_configuration=population_configuration,
                migration_matrix=migration_matrix)
        for bad_type in [{}, None, 2, [""], [[]], [None]]:
            self.assertRaises(
                TypeError, _msprime.Simulator, get_samples(2), rng,
                population_configuration=bad_type)
        # Cannot have empty list
        self.assertRaises(ValueError, f, 2, [])
        # Must provide population_configuration if a migration_matrix
        # is supplied.
        self.assertRaises(
            ValueError, _msprime.Simulator, get_samples(2), rng,
            migration_matrix=[0, 0, 0, 0])

    def test_get_population_configurations(self):
        def f(sample_size, conf_tuples):
            population_configuration = [
                get_population_configuration(initial_size=p, growth_rate=a)
                for p, a in conf_tuples]
            N = len(population_configuration)
            migration_matrix = [0 for j in range(N) for k in range(N)]
            s = _msprime.Simulator(
                get_samples(sample_size), _msprime.RandomGenerator(1),
                population_configuration=population_configuration,
                migration_matrix=migration_matrix)
            conf_dicts = s.get_population_configuration()
            self.assertEqual(len(conf_dicts), len(conf_tuples))
            for conf_dict, conf_tuple in zip(conf_dicts, conf_tuples):
                self.assertEqual(len(conf_dict), 2)
                self.assertEqual(conf_dict["initial_size"], conf_tuple[0])
                self.assertEqual(conf_dict["growth_rate"], conf_tuple[1])
        f(2, [(1, 1)])
        f(2, [(2, 0), (0.5, 0.1)])
        f(5, [(1, 0.25), (2, 0.5), (3, 0.75), (4, 1)])

    def test_bad_migration_matrix(self):
        def f(num_populations, migration_matrix):
            population_configuration = [
                get_population_configuration()
                for j in range(num_populations)]
            population_configuration[0]["sample_size"] = 2
            return _msprime.Simulator(
                get_samples(2), _msprime.RandomGenerator(1),
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
                ValueError, _msprime.Simulator, get_samples(2),
                _msprime.RandomGenerator(1),
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
                    get_samples(2), _msprime.RandomGenerator(1),
                    migration_matrix=migration_matrix,
                    population_configuration=population_configuration)
                self.assertEqual(migration_matrix, sim.get_migration_matrix())

    def test_bad_demographic_event_types(self):
        def f(events):
            _msprime.Simulator(
                get_samples(2), _msprime.RandomGenerator(1),
                demographic_events=events)
        event_generators = [
            get_size_change_event, get_growth_rate_change_event,
            get_migration_rate_change_event,
            get_mass_migration_event, get_simple_bottleneck_event,
            get_instantaneous_bottleneck_event]
        for bad_type in [None, {}, "", 1]:
            self.assertRaises(TypeError, f, bad_type)
        for bad_type in [None, "", 1, []]:
            self.assertRaises(TypeError, f, [bad_type])
        for bad_type in [None, [], 0]:
            for generator in event_generators:
                event = generator()
                event["type"] = bad_type
                self.assertRaises(ValueError, f, [event])
        for bad_event in [b'', b'1', b'x' * 1000, b'Size_change', 2, 'none']:
            for generator in event_generators:
                event["type"] = bad_event
                self.assertRaises(ValueError, f, [event])
        for bad_type in [[], "", {}]:
            for generator in event_generators:
                event = generator(time=bad_type)
                self.assertRaises(TypeError, f, [event])
            event = get_size_change_event(population_id=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_size_change_event(size=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_instantaneous_bottleneck_event(population_id=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_instantaneous_bottleneck_event(strength=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_simple_bottleneck_event(population_id=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_simple_bottleneck_event(proportion=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_mass_migration_event(source=bad_type, destination=0)
            self.assertRaises(TypeError, f, [event])
            event = get_mass_migration_event(source=0, destination=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_mass_migration_event(
                source=0, destination=1, proportion=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_migration_rate_change_event(matrix_index=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_migration_rate_change_event(
                matrix_index=0, migration_rate=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_growth_rate_change_event(population_id=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_growth_rate_change_event(population_id=0, growth_rate=bad_type)
            self.assertRaises(TypeError, f, [event])

    def test_bad_demographic_event_values(self):
        def f(events, num_populations=1):
            population_configuration = [get_population_configuration(2)] + [
                get_population_configuration(0)
                for _ in range(num_populations - 1)]
            _msprime.Simulator(
                get_samples(2), _msprime.RandomGenerator(1),
                demographic_events=events,
                population_configuration=population_configuration,
                migration_matrix=get_migration_matrix(num_populations))
        event_generators = [
            get_size_change_event, get_growth_rate_change_event,
            get_migration_rate_change_event,
            get_mass_migration_event, get_simple_bottleneck_event,
            get_instantaneous_bottleneck_event]
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
            event = get_simple_bottleneck_event(population_id=bad_pop_id)
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
        # Tests specific for mass_migration and bottleneck
        event = get_mass_migration_event(source=0, destination=0)
        self.assertRaises(_msprime.InputError, f, [event])
        for bad_proportion in [-1, 1.1, 1e7]:
            event = get_mass_migration_event(proportion=bad_proportion)
            self.assertRaises(_msprime.InputError, f, [event])
            event = get_simple_bottleneck_event(proportion=bad_proportion)
            self.assertRaises(_msprime.InputError, f, [event])

    def test_unsorted_demographic_events(self):
        event_generators = [
            get_size_change_event, get_growth_rate_change_event,
            get_migration_rate_change_event, get_mass_migration_event,
            get_simple_bottleneck_event]
        events = []
        for event_generator in event_generators:
            for _ in range(3):
                events.append(event_generator(time=random.random()))
        sorted_events = sorted(events, key=lambda x: x["time"])
        self.assertNotEqual(events, sorted_events)
        self.assertRaises(
            _msprime.InputError, _msprime.Simulator, get_samples(10),
            _msprime.RandomGenerator(1), demographic_events=events)

    def test_seed_equality(self):
        simulations = [
            {
                "samples": get_samples(10),
            }, {
                "samples": get_samples(100),
                "demographic_events": [
                    get_simple_bottleneck_event(0.01, 0, 1.0)],
            }, {
                "samples": get_samples(10), "num_loci": 100,
                "scaled_recombination_rate": 0.1,
            }, {
                "samples": get_population_samples(3, 3, 4),
                "population_configuration": [
                    get_population_configuration(),
                    get_population_configuration(),
                    get_population_configuration()],
                "migration_matrix": get_migration_matrix(3, 3),
                "demographic_events": [
                    get_size_change_event(0.1, 0.5),
                    get_migration_rate_change_event(0.2, 2),
                    get_growth_rate_change_event(0.3, 2),
                    get_mass_migration_event(0.4, 0, 1, 0.5),
                    get_simple_bottleneck_event(0.5, 0.5)]
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
            get_population_configuration(),
            get_population_configuration()]
        sim = _msprime.Simulator(
            get_population_samples(5, 5), _msprime.RandomGenerator(1),
            population_configuration=population_configuration,
            migration_matrix=[0, 0, 0, 0])
        self.assertRaises(_msprime.LibraryError, sim.run)

    def test_simple_event_counters(self):
        for n in [2, 10, 20]:
            sim = _msprime.Simulator(
                get_samples(n), _msprime.RandomGenerator(1),
                scaled_recombination_rate=0)
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
            get_samples(n), _msprime.RandomGenerator(1),
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
            get_samples(n), _msprime.RandomGenerator(1),
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
        sample_sizes = [0 for _ in range(num_populations)]
        sample_sizes[active_pops[0]] = 2
        sample_sizes[active_pops[1]] = 2
        sample_sizes[active_pops[2]] = 6
        flattened = [x for row in migration_matrix for x in row]
        sim = _msprime.Simulator(
            get_population_samples(*sample_sizes),
            _msprime.RandomGenerator(1),
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
            get_samples(n), _msprime.RandomGenerator(1),
            population_configuration=[
                get_population_configuration(),
                get_population_configuration()],
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

    def test_bottleneck(self):
        n = 10
        t1 = 0.01
        t2 = 0.02
        t3 = 0.03
        sim = _msprime.Simulator(
            get_population_samples(n, n),
            _msprime.RandomGenerator(1),
            population_configuration=[
                get_population_configuration(),
                get_population_configuration()],
            demographic_events=[
                get_simple_bottleneck_event(t1, population_id=0, proportion=1),
                get_simple_bottleneck_event(t2, population_id=1, proportion=1),
                get_mass_migration_event(
                    t3, source=0, destination=1, proportion=1)],
            migration_matrix=[0, 0, 0, 0])
        sim.run(t1)
        pop_sizes = [0, 0]
        for ind in sim.get_ancestors():
            for _, _, _, pop_id in ind:
                pop_sizes[pop_id] += 1
        self.assertEqual(pop_sizes[0], 1)
        sim.run(t2)
        pop_sizes = [0, 0]
        for ind in sim.get_ancestors():
            for _, _, _, pop_id in ind:
                pop_sizes[pop_id] += 1
        self.assertEqual(pop_sizes[1], 1)
        sim.run()
        self.assertGreater(sim.get_time(), t3)

    def test_recombination_event_counters(self):
        n = 10
        # No migration
        population_configuration = [
            get_population_configuration(n),
            get_population_configuration(0)]
        sim = _msprime.Simulator(
            get_samples(n), _msprime.RandomGenerator(1),
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
            get_population_configuration(),
            get_population_configuration(),
            get_population_configuration()]
        sim = _msprime.Simulator(
            get_population_samples(5, 5, 0), _msprime.RandomGenerator(1),
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
        sim = _msprime.Simulator(get_samples(10), _msprime.RandomGenerator(1))
        times = set()
        for _ in range(10):
            sim.run()
            t = sim.get_time()
            self.assertNotIn(t, times)
            times.add(t)
            sim.reset()
            self.assertEqual(sim.get_time(), 0)

    def test_populate_tables_interface(self):
        node_table = _msprime.NodeTable()
        edgeset_table = _msprime.EdgesetTable()
        migration_table = _msprime.MigrationTable()
        sim = _msprime.Simulator(get_samples(10), _msprime.RandomGenerator(1))
        sim.run()

        recomb_map = uniform_recombination_map(sim)
        # nodes, edgesets and migrations are mandatory.
        self.assertRaises(TypeError, sim.populate_tables)
        self.assertRaises(TypeError, sim.populate_tables, nodes=node_table)
        self.assertRaises(
            TypeError, sim.populate_tables, nodes=node_table, migrations=migration_table)
        self.assertRaises(
            TypeError, sim.populate_tables, nodes=node_table, edgesets=edgeset_table)

        kwargs = {
            "nodes": node_table,
            "edgesets": edgeset_table,
            "migrations": migration_table
        }
        for bad_type in ["", {}, [], None]:
            for k in kwargs.keys():
                d = dict(kwargs)
                d[k] = bad_type
                self.assertRaises(TypeError, sim.populate_tables, **d)
            self.assertRaises(
                TypeError, sim.populate_tables, recombination_map=bad_type, **kwargs)
            self.assertRaises(
                TypeError, sim.populate_tables, Ne=bad_type, **kwargs)

        sim.populate_tables(**kwargs)
        self.assertEqual(edgeset_table.num_rows, sim.get_num_coalescence_records())
        kwargs["recombination_map"] = recomb_map
        sim.populate_tables(**kwargs)
        self.assertEqual(edgeset_table.num_rows, sim.get_num_coalescence_records())
        kwargs["Ne"] = 0.25
        sim.populate_tables(**kwargs)
        self.assertEqual(edgeset_table.num_rows, sim.get_num_coalescence_records())


class TestTreeSequence(LowLevelTestCase):
    """
    Tests for the low-level interface for the TreeSequence.
    """
    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msp_ll_ts_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)

    def test_seed_equality(self):
        simulations = [
            {
                "samples": get_samples(10),
            }, {
                "samples": get_samples(10), "num_loci": 100,
                "scaled_recombination_rate": 0.1,
            },
        ]
        for params in simulations:
            rng1 = _msprime.RandomGenerator(1)
            rng2 = _msprime.RandomGenerator(1)
            sim1 = _msprime.Simulator(random_generator=rng1, **params)
            sim2 = _msprime.Simulator(random_generator=rng2, **params)
            sim1.run()
            sim2.run()
            mutgen1 = _msprime.MutationGenerator(rng1, 1)
            mutgen2 = _msprime.MutationGenerator(rng2, 1)
            t1 = populate_tree_sequence(sim1, mutation_generator=mutgen1)
            t2 = populate_tree_sequence(sim2, mutation_generator=mutgen2)
            self.assertEqual(t1.get_num_edgesets(), t2.get_num_edgesets())
            r1 = [t1.get_record(j) for j in range(t1.get_num_edgesets())]
            r2 = [t2.get_record(j) for j in range(t2.get_num_edgesets())]
            self.assertEqual(r1, r2)
            m1 = [t1.get_mutation(j) for j in range(t1.get_num_mutations())]
            m2 = [t2.get_mutation(j) for j in range(t2.get_num_mutations())]
            self.assertEqual(m1, m2)

    @unittest.skipIf(IS_WINDOWS, "File permissions on Windows")
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
        # Check the initial state to make sure that it is empty.
        ts = _msprime.TreeSequence()
        self.assertEqual(ts.get_sample_size(), 0)
        self.assertEqual(ts.get_sequence_length(), 0)
        self.assertEqual(ts.get_num_trees(), 0)
        self.assertEqual(ts.get_num_edgesets(), 0)
        self.assertEqual(ts.get_num_mutations(), 0)
        self.assertEqual(ts.get_num_migrations(), 0)
        self.verify_dump_equality(ts)

    def test_num_nodes(self):
        for ts in self.get_example_tree_sequences():
            max_node = 0
            for j in range(ts.get_num_edgesets()):
                _, _, node, _, _, _ = ts.get_record(j)
                if node > max_node:
                    max_node = node
            self.assertEqual(max_node + 1, ts.get_num_nodes())

    def verify_dump_equality(self, ts):
        """
        Verifies that we can dump a copy of the specified tree sequence
        to the specified file, and load an identical copy.
        """
        ts.dump(self.temp_file)
        ts2 = _msprime.TreeSequence()
        ts2.load(self.temp_file)
        self.assertEqual(ts.get_sample_size(), ts2.get_sample_size())
        self.assertEqual(ts.get_sequence_length(), ts2.get_sequence_length())
        self.assertEqual(ts.get_num_mutations(), ts2.get_num_mutations())
        self.assertEqual(ts.get_num_nodes(), ts2.get_num_nodes())
        records1 = [
            ts.get_record(j) for j in range(ts.get_num_edgesets())]
        records2 = [
            ts2.get_record(j) for j in range(ts2.get_num_edgesets())]
        self.assertEqual(records1, records2)
        mutations1 = [ts.get_mutation(j) for j in range(ts.get_num_mutations())]
        mutations2 = [ts2.get_mutation(j) for j in range(ts2.get_num_mutations())]
        self.assertEqual(mutations1, mutations2)
        self.assertEqual(
            ts.get_provenance_strings(), ts2.get_provenance_strings())

    def test_dump_equality(self):
        for ts in self.get_example_tree_sequences():
            self.verify_dump_equality(ts)

    def verify_mutations(self, ts):
        mutations = [ts.get_mutation(j) for j in range(ts.get_num_mutations())]
        self.assertGreater(ts.get_num_mutations(), 0)
        self.assertEqual(len(mutations), ts.get_num_mutations())
        # Check the form of the mutations
        for j, (position, nodes, index) in enumerate(mutations):
            self.assertEqual(j, index)
            for node in nodes:
                self.assertIsInstance(node, int)
                self.assertGreaterEqual(node, 0)
                self.assertLessEqual(node, ts.get_num_nodes())
            self.assertIsInstance(position, float)
            self.assertGreater(position, 0)
            self.assertLess(position, ts.get_sequence_length())
        # mutations must be sorted by position order.
        self.assertEqual(mutations, sorted(mutations))

    def test_get_record_interface(self):
        for ts in self.get_example_tree_sequences():
            num_records = ts.get_num_edgesets()
            # We don't accept Python negative indexes here.
            self.assertRaises(IndexError, ts.get_record, -1)
            for j in [0, 10, 10**6]:
                self.assertRaises(IndexError, ts.get_record, num_records + j)
            for x in [None, "", {}, []]:
                self.assertRaises(TypeError, ts.get_record, x)

    def test_get_node_interface(self):
        for ts in self.get_example_tree_sequences():
            num_nodes = ts.get_num_nodes()
            # We don't accept Python negative indexes here.
            self.assertRaises(IndexError, ts.get_node, -1)
            for j in [0, 10, 10**6]:
                self.assertRaises(IndexError, ts.get_node, num_nodes + j)
            for x in [None, "", {}, []]:
                self.assertRaises(TypeError, ts.get_node, x)

    def test_get_edgeset_interface(self):
        for ts in self.get_example_tree_sequences():
            num_edgesets = ts.get_num_edgesets()
            # We don't accept Python negative indexes here.
            self.assertRaises(IndexError, ts.get_edgeset, -1)
            for j in [0, 10, 10**6]:
                self.assertRaises(IndexError, ts.get_edgeset, num_edgesets + j)
            for x in [None, "", {}, []]:
                self.assertRaises(TypeError, ts.get_edgeset, x)

    def test_get_migration_interface(self):
        ts = self.get_example_migration_tree_sequence()
        for bad_type in ["", None, {}]:
            self.assertRaises(TypeError, ts.get_migration, bad_type)
        num_records = ts.get_num_migrations()
        # We don't accept Python negative indexes here.
        self.assertRaises(IndexError, ts.get_migration, -1)
        for j in [0, 10, 10**6]:
            self.assertRaises(IndexError, ts.get_migration, num_records + j)

    def test_migrations(self):
        sim = get_example_simulator(10, num_populations=3, store_migrations=True)
        sim.run()
        ts = populate_tree_sequence(sim)
        self.assertGreater(sim.get_num_migrations(), 0)
        self.assertEqual(
            sim.get_num_migrations(), ts.get_num_migrations())
        sim_records = sim.get_migrations()
        ts_records = [
            ts.get_migration(j) for j in range(ts.get_num_migrations())]
        self.assertEqual(len(sim_records), len(ts_records))
        for sim_r, ts_r in zip(sim_records, ts_records):
            self.assertEqual(len(sim_r), len(ts_r))
            self.assertEqual(len(sim_r), 6)
            # Left and right have been remapped, so might be slightly different.
            self.assertAlmostEqual(sim_r[0], ts_r[0])
            self.assertAlmostEqual(sim_r[1], ts_r[1])
            for j in range(2, len(ts_r)):
                self.assertEqual(sim_r[j], ts_r[j])

    def test_record_scaling(self):
        sim = get_example_simulator(10, num_populations=2, store_migrations=True)
        sim.run()
        for Ne in [0.25, 1, 10, 1e6]:
            tree_sequence = populate_tree_sequence(sim, Ne=Ne)
            self.assertEqual(
                sim.get_num_coalescence_records(),
                tree_sequence.get_num_edgesets())
            sim_times = [
                r[-2] for r in sim.get_coalescence_records()]
            for j in range(tree_sequence.get_num_edgesets()):
                generation = tree_sequence.get_record(j)[-2]
                self.assertAlmostEqual(generation, sim_times[j] * 4 * Ne)
            self.assertGreater(sim.get_num_migrations(), 0)
            self.assertEqual(
                sim.get_num_migrations(),
                tree_sequence.get_num_migrations())
            sim_times = [
                r[-1] for r in sim.get_migrations()]
            for j in range(tree_sequence.get_num_migrations()):
                generation = tree_sequence.get_migration(j)[-1]
                self.assertAlmostEqual(generation, sim_times[j] * 4 * Ne)

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
            self.assertRaises(
                _msprime.LibraryError, ts.get_pairwise_diversity, [0, 0])
            samples = list(range(ts.get_sample_size()))
            pi1 = ts.get_pairwise_diversity(samples)
            self.assertGreaterEqual(pi1, 0)

    @unittest.skip("SKIP simplify")
    def test_simplify(self):
        out = _msprime.TreeSequence()
        for ts in self.get_example_tree_sequences():
            for bad_type in ["", None, {}]:
                self.assertRaises(TypeError, ts.simplify)
                self.assertRaises(TypeError, ts.simplify, out, bad_type)
            self.assertRaises(ValueError, ts.simplify, out, [])
            self.assertRaises(ValueError, ts.simplify, out, [0])
            self.assertRaises(ValueError, ts.simplify, out, [0, ts.get_sample_size()])
            self.assertRaises(_msprime.LibraryError, ts.simplify, out, [0, 0])
            s1 = _msprime.TreeSequence()
            ts.simplify(s1, [0, 1])
            s2 = _msprime.TreeSequence()
            ts.simplify(s2, [1, 0])
            self.assertEqual(ts.get_sequence_length(), s1.get_sequence_length())
            self.assertEqual(s1.get_sample_size(), s2.get_sample_size())
            self.assertEqual(s1.get_num_edgesets(), s2.get_num_edgesets())
            self.assertEqual(s1.get_num_mutations(), s2.get_num_mutations())
            self.assertEqual(s1.get_num_trees(), s2.get_num_trees())

    def test_provenance_strings_populate(self):
        rng = _msprime.RandomGenerator(1)
        sim = _msprime.Simulator(get_samples(10), rng)
        sim.run()
        for bad_type in [{}, rng, None, 5]:
            provenance_strings = [bad_type]
            self.assertRaises(
                TypeError, populate_tree_sequence, sim,
                provenance_strings=provenance_strings)

        for j in range(10):
            strings = []
            for k in range(j):
                # Use some fairly big strings to stress things out
                strings.append(b"x" * (k + 1) * 8192)
            ts = populate_tree_sequence(sim, provenance_strings=strings)
            self.assertEqual(ts.get_provenance_strings(), strings)
            ts.dump(self.temp_file)
            ts2 = _msprime.TreeSequence()
            ts2.load(self.temp_file)
            self.assertEqual(ts2.get_provenance_strings(), strings)


class TestNewickConverter(LowLevelTestCase):
    """
    Tests for the low-level newick converter.
    """
    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        conv = _msprime.NewickConverter(ts)
        self.assertEqual(list(conv), [])

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.NewickConverter)
        self.assertRaises(TypeError, _msprime.NewickConverter, None)
        sim = _msprime.Simulator(get_samples(10), _msprime.RandomGenerator(1))
        sim.run()
        ts = populate_tree_sequence(sim)
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

    def test_nonbinary_trees(self):
        ts = self.get_nonbinary_tree_sequence()

        def f():
            return list(_msprime.NewickConverter(ts))
        self.assertRaises(_msprime.LibraryError, f)


class TestVcfConverter(LowLevelTestCase):
    """
    Tests for the low-level vcf converter.
    """
    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        converter = _msprime.VcfConverter(ts)
        header = converter.get_header()
        self.assertGreater(len(header), 0)
        self.assertTrue(header.startswith("##fileformat"))
        self.assertEqual(list(converter), [])

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.VcfConverter)
        self.assertRaises(TypeError, _msprime.VcfConverter, None)
        rng = _msprime.RandomGenerator(1)
        sim = _msprime.Simulator(get_samples(10), rng)
        mutgen = _msprime.MutationGenerator(rng, 10)
        ts = _msprime.TreeSequence()
        sim.run()
        ts = populate_tree_sequence(sim, mutation_generator=mutgen)
        self.assertGreater(ts.get_num_mutations(), 0)
        for bad_type in [None, "", [], {}]:
            self.assertRaises(
                TypeError, _msprime.VcfConverter, ts, ploidy=bad_type)
        converter = _msprime.VcfConverter(ts)
        before = converter.get_header() + "".join(converter)
        self.assertGreater(len(before), 0)
        converter = _msprime.VcfConverter(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = converter.get_header() + "".join(converter)
        self.assertEqual(before, after)

    def test_ploidy(self):
        ts = self.get_tree_sequence(sample_size=10)
        for bad_ploidy in [-1, 3, 4, 11, 20, 10**6]:
            self.assertRaises(
                _msprime.LibraryError, _msprime.VcfConverter, ts,
                bad_ploidy)

    def test_positions(self):
        # Make sure we successfully discretise the positions when
        # we have lots of mutations
        ts = self.get_tree_sequence(sample_size=10, mutation_rate=10)
        converter = _msprime.VcfConverter(ts)
        positions = []
        for row in converter:
            pos = float(row.split()[1])
            self.assertEqual(pos, int(pos))
            positions.append(pos)
        self.assertEqual(sorted(positions), positions)
        self.assertEqual(len(set(positions)), len(positions))

    def test_rows(self):
        ts = self.get_tree_sequence()
        num_mutations = ts.get_num_mutations()
        converter = _msprime.VcfConverter(ts)
        self.assertGreater(num_mutations, 0)
        num_rows = 0
        for row in converter:
            num_rows += 1
        self.assertEqual(num_rows, num_mutations)

    def test_header(self):
        ts = self.get_tree_sequence()
        converter = _msprime.VcfConverter(ts)
        header = converter.get_header()
        # We do more indepth testing elsewhere, just make sure
        # it roughly makes sense.
        self.assertGreater(len(header), 0)
        self.assertTrue(header.startswith("##fileformat"))

    def test_iterator(self):
        ts = self.get_tree_sequence()
        iters = [
            _msprime.VcfConverter(ts), _msprime.VcfConverter(ts, 1),
            _msprime.VcfConverter(ts, 2)]
        for iterator in iters:
            self.verify_iterator(iterator)


class TestTreeDiffIterator(LowLevelTestCase):
    """
    Tests for the low-level tree diff iterator.
    """
    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        iterator = _msprime.TreeDiffIterator(ts)
        self.assertEqual(list(iterator), [])

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.TreeDiffIterator)
        self.assertRaises(TypeError, _msprime.TreeDiffIterator, None)
        sim = _msprime.Simulator(get_samples(10), _msprime.RandomGenerator(1))
        sim.run()
        ts = populate_tree_sequence(sim)
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
    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        tree = _msprime.SparseTree(ts)
        iterator = _msprime.SparseTreeIterator(tree)
        self.assertEqual(list(iterator), [])

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.SparseTreeIterator)
        self.assertRaises(TypeError, _msprime.SparseTreeIterator, None)
        ts = _msprime.TreeSequence()
        self.assertRaises(TypeError, _msprime.SparseTreeIterator, ts)
        sim = _msprime.Simulator(get_samples(10), _msprime.RandomGenerator(1))
        sim.run()
        ts = populate_tree_sequence(sim)
        tree = _msprime.SparseTree(ts)
        n_before = 0
        parents_before = []
        for t in _msprime.SparseTreeIterator(tree):
            n_before += 1
            self.assertIs(t, tree)
            pi = {}
            for j in range(t.get_num_nodes()):
                pi[j] = t.get_parent(j)
            parents_before.append(pi)
        self.assertEqual(n_before, len(list(_msprime.TreeDiffIterator(ts))))
        # If we remove the objects, we should get the same results.
        iterator = _msprime.SparseTreeIterator(tree)
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
        self.verify_iterator(_msprime.SparseTreeIterator(tree))

    def test_root_bug(self):
        # Reproduce a simulation that provoked a root calculation bug.
        params = {
            "random_generator": _msprime.RandomGenerator(878638576),
            "samples": get_samples(64),
            "num_loci": 57,
            "scaled_recombination_rate": 0.192184324155680
        }
        sim = _msprime.Simulator(**params)
        sim.run()
        ts = populate_tree_sequence(sim)
        st = _msprime.SparseTree(ts)
        for st in _msprime.SparseTreeIterator(st):
            root = 0
            while st.get_parent(root) != NULL_NODE:
                root = st.get_parent(root)
            self.assertEqual(st.get_root(), root)


class TestHaplotypeGenerator(LowLevelTestCase):
    """
    Tests for the low-level haplotype generator.
    """
    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        hg = _msprime.HaplotypeGenerator(ts)
        self.assertRaises(IndexError, hg.get_haplotype, 0)

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.HaplotypeGenerator)
        ts = self.get_tree_sequence(num_loci=10)
        for bad_type in ["", {}, [], None]:
            self.assertRaises(TypeError, _msprime.HaplotypeGenerator, ts, bad_type)
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
    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        vg = _msprime.VariantGenerator(ts, bytearray())
        self.assertEqual(list(vg), [])

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.VariantGenerator)
        ts = self.get_tree_sequence(num_loci=10)
        buff = bytearray(ts.get_sample_size())
        for bad_type in ["", {}, [], None]:
            self.assertRaises(TypeError, _msprime.VariantGenerator, bad_type, buff)
            self.assertRaises(TypeError, _msprime.VariantGenerator, ts, bad_type)
            self.assertRaises(TypeError, _msprime.VariantGenerator, ts, buff, bad_type)
        for size in [0, 1, ts.get_sample_size() - 1]:
            buff = bytearray(size)
            self.assertRaises(BufferError, _msprime.VariantGenerator, ts, buff)

        buff = bytearray(ts.get_sample_size())
        vg = _msprime.VariantGenerator(ts, buff)
        before = list(vg)
        vg = _msprime.VariantGenerator(ts, buff)
        del ts
        # We should keep a reference to the tree sequence.
        after = list(vg)
        self.assertEqual(before, after)

    def test_buffer_nastiness(self):
        ts = self.get_tree_sequence(num_loci=10)
        buff = bytearray(ts.get_sample_size())
        variants = list(_msprime.VariantGenerator(ts, buff))
        del buff
        j = 0
        for _ in enumerate(variants):
            j += 1
        self.assertEqual(j, ts.get_num_mutations())

        buff = bytearray(ts.get_sample_size())
        variants = list(_msprime.VariantGenerator(ts, buff))
        buff.extend(0 for j in range(1000))
        j = 0
        for _ in enumerate(variants):
            j += 1
        self.assertEqual(j, ts.get_num_mutations())

        buff = bytearray(ts.get_sample_size())
        variants = list(_msprime.VariantGenerator(ts, buff))
        buff.pop()
        j = 0
        for _ in enumerate(variants):
            j += 1
        self.assertEqual(j, ts.get_num_mutations())

    def test_form(self):
        ts = self.get_tree_sequence(num_loci=10)
        buff = bytearray(ts.get_sample_size())
        variants = list(_msprime.VariantGenerator(ts, buff))
        self.assertGreater(len(variants), 0)
        self.assertEqual(len(variants), ts.get_num_sites())
        sites = [ts.get_site(j) for j in range(ts.get_num_sites())]
        self.assertEqual(variants, sites)
        for _ in _msprime.VariantGenerator(ts, buff):
            self.assertEqual(len(buff), ts.get_sample_size())
            for b in buff:
                self.assertIn(b, [0, 1])
        for _ in _msprime.VariantGenerator(ts, buff, False):
            self.assertEqual(len(buff), ts.get_sample_size())
            for b in buff:
                self.assertIn(b, [0, 1])
        for _ in _msprime.VariantGenerator(ts, buff, True):
            self.assertEqual(len(buff), ts.get_sample_size())
            for b in buff:
                self.assertIn(b, [ord('0'), ord('1')])

    def test_iterator(self):
        ts = self.get_tree_sequence()
        buff = bytearray(ts.get_sample_size())
        variants = _msprime.VariantGenerator(ts, buff)
        self.verify_iterator(variants)


class TestSparseTree(LowLevelTestCase):
    """
    Tests on the low-level sparse tree interface.
    """
    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        st = _msprime.SparseTree(ts)
        self.assertEqual(list(_msprime.SparseTreeIterator(st)), [])

    def test_flags(self):
        ts = self.get_tree_sequence()
        st = _msprime.SparseTree(ts)
        self.assertEqual(st.get_flags(), 0)
        # We should still be able to count the leaves, just inefficiently.
        self.assertEqual(st.get_num_leaves(0), 1)
        self.assertRaises(_msprime.LibraryError, st.get_num_tracked_leaves, 0)
        self.assertRaises(_msprime.LibraryError, _msprime.LeafListIterator, st, 0)
        all_flags = [
            0, _msprime.LEAF_COUNTS, _msprime.LEAF_LISTS,
            _msprime.LEAF_COUNTS | _msprime.LEAF_LISTS]
        for flags in all_flags:
            st = _msprime.SparseTree(ts, flags=flags)
            self.assertEqual(st.get_flags(), flags)
            self.assertEqual(st.get_num_leaves(0), 1)
            if flags & _msprime.LEAF_COUNTS:
                self.assertEqual(st.get_num_tracked_leaves(0), 0)
            else:
                self.assertRaises(_msprime.LibraryError, st.get_num_tracked_leaves, 0)
            if flags & _msprime.LEAF_LISTS:
                leaves = list(_msprime.LeafListIterator(st, 0))
                self.assertEqual(leaves, [0])
            else:
                self.assertRaises(
                    _msprime.LibraryError, _msprime.LeafListIterator, st, 0)

    def test_sites(self):
        for ts in self.get_example_tree_sequences():
            st = _msprime.SparseTree(ts)
            all_sites = [ts.get_site(j) for j in range(ts.get_num_sites())]
            all_tree_sites = []
            j = 0
            for st in _msprime.SparseTreeIterator(st):
                tree_sites = st.get_sites()
                self.assertEqual(st.get_num_sites(), len(tree_sites))
                all_tree_sites.extend(tree_sites)
                for position, ancestral_state, mutations, index in tree_sites:
                    self.assertTrue(st.get_left() <= position < st.get_right())
                    self.assertEqual(index, j)
                    for site, node, derived_state in mutations:
                        self.assertEqual(site, index)
                        self.assertNotEqual(st.get_parent(node), NULL_NODE)
                    j += 1
            self.assertEqual(all_tree_sites, all_sites)

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.SparseTree)
        for bad_type in ["", {}, [], None, 0]:
            self.assertRaises(
                TypeError, _msprime.SparseTree, bad_type)
        ts = self.get_tree_sequence()
        for bad_type in ["", {}, True, 1, None]:
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, tracked_leaves=bad_type)
        for bad_type in ["", {}, None, []]:
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, flags=bad_type)
        for n in range(1, 10):
            ts = self.get_tree_sequence(sample_size=10, num_loci=n)
            st = _msprime.SparseTree(ts)
            self.assertEqual(st.get_num_nodes(), ts.get_num_nodes())
            self.assertEqual(st.get_sample_size(), ts.get_sample_size())
            # An uninitialised sparse tree should always be zero.
            self.assertEqual(st.get_root(), 0)
            self.assertEqual(st.get_left(), 0)
            self.assertEqual(st.get_right(), 0)
            for j in range(n):
                self.assertEqual(st.get_parent(j), NULL_NODE)
                self.assertEqual(st.get_population(j), 0)
                self.assertEqual(st.get_children(j), tuple())
                self.assertEqual(st.get_time(j), 0)

    def test_memory_error(self):
        # This provokes a bug where we weren't reference counting
        # the tree sequence properly, and the underlying memory for a
        # sparse tree was getting corrupted.
        for n in range(1, 10):
            ts = self.get_tree_sequence(sample_size=100, num_loci=n)
            num_nodes = ts.get_num_nodes()
            st = _msprime.SparseTree(ts)
            # deleting the tree sequence should still give a well formed
            # sparse tree.
            st_iter = _msprime.SparseTreeIterator(st)
            next(st_iter)
            del ts
            del st_iter
            # Do a quick traversal just to exercise the tree
            stack = [st.get_root()]
            while len(stack) > 0:
                u = stack.pop()
                self.assertLess(u, num_nodes)
                stack.extend(st.get_children(u))

    def test_bad_tracked_leaves(self):
        ts = self.get_tree_sequence()
        flags = _msprime.LEAF_COUNTS
        for bad_type in ["", {}, [], None]:
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, flags=flags,
                tracked_leaves=[bad_type])
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, flags=flags,
                tracked_leaves=[1, bad_type])
        for bad_leaf in [ts.get_sample_size(), 10**6, -1e6]:
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts, flags=flags,
                tracked_leaves=[bad_leaf])
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts, flags=flags,
                tracked_leaves=[1, bad_leaf])
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts,
                tracked_leaves=[1, bad_leaf, 1])

    def test_count_all_leaves(self):
        for ts in self.get_example_tree_sequences():
            self.verify_iterator(_msprime.TreeDiffIterator(ts))
            st = _msprime.SparseTree(ts, flags=_msprime.LEAF_COUNTS)
            # Without initialisation we should be 0 leaves for every node
            # that is not a leaf.
            for j in range(st.get_num_nodes()):
                l = 1 if j < st.get_sample_size() else 0
                self.assertEqual(st.get_num_leaves(j), l)
                self.assertEqual(st.get_num_tracked_leaves(j), 0)
            # Now, try this for a tree sequence.
            for st in _msprime.SparseTreeIterator(st):
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
        examples = [
            self.get_tree_sequence(sample_size=5, num_loci=10),
            self.get_tree_sequence(
                sample_size=12,
                demographic_events=[
                    get_simple_bottleneck_event(0.2, proportion=1.0)])]
        # Ensure that the second example does have some non-binary records
        ts = examples[1]
        found = False
        for j in range(ts.get_num_edgesets()):
            r = ts.get_record(j)
            if len(r[3]) > 2:
                found = True
        self.assertTrue(found)
        for ts in examples:
            leaves = [j for j in range(ts.get_sample_size())]
            powerset = itertools.chain.from_iterable(
                itertools.combinations(leaves, r)
                for r in range(len(leaves) + 1))
            for subset in map(list, powerset):
                # Ordering shouldn't make any different.
                random.shuffle(subset)
                st = _msprime.SparseTree(
                    ts, flags=_msprime.LEAF_COUNTS, tracked_leaves=subset)
                for st in _msprime.SparseTreeIterator(st):
                    nu = get_tracked_leaf_counts(st, subset)
                    nu_prime = [
                        st.get_num_tracked_leaves(j) for j in
                        range(st.get_num_nodes())]
                    self.assertEqual(nu, nu_prime)
            # Passing duplicated values should raise an error
            leaf = 1
            for j in range(2, 20):
                tracked_leaves = [leaf for _ in range(j)]
                self.assertRaises(
                    _msprime.LibraryError, _msprime.SparseTree,
                    ts, flags=_msprime.LEAF_COUNTS,
                    tracked_leaves=tracked_leaves)

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
            for index, st in enumerate(_msprime.SparseTreeIterator(st)):
                self.assertEqual(index, st.get_index())

    def test_bad_mutations(self):
        ts = self.get_tree_sequence(2, num_loci=200, random_seed=1)
        node_table = _msprime.NodeTable()
        edgeset_table = _msprime.EdgesetTable()
        site_table = _msprime.SiteTable()
        mutation_table = _msprime.MutationTable()
        migration_table = _msprime.MigrationTable()
        ts.dump_tables(nodes=node_table, edgesets=edgeset_table)

        def f(mutations):
            position = []
            node = []
            site = []
            ancestral_state = []
            ancestral_state_length = []
            derived_state = []
            derived_state_length = []
            for j, (p, n) in enumerate(mutations):
                site.append(j)
                position.append(p)
                ancestral_state.append("0")
                ancestral_state_length.append(1)
                derived_state.append("1")
                derived_state_length.append(1)
                node.append(n)
            site_table.set_columns(
                position=position, ancestral_state=ancestral_state,
                ancestral_state_length=ancestral_state_length)
            mutation_table.set_columns(
                site=site, node=node, derived_state=derived_state,
                derived_state_length=derived_state_length)
            ts2 = _msprime.TreeSequence()
            ts2.load_tables(
                nodes=node_table, edgesets=edgeset_table, migrations=migration_table,
                sites=site_table, mutations=mutation_table)
        self.assertRaises(_msprime.LibraryError, f, [(0.1, -1)])
        l = ts.get_sequence_length()
        u = ts.get_num_nodes()
        for bad_node in [u, u + 1, 2 * u]:
            self.assertRaises(_msprime.LibraryError, f, [(0.1, bad_node)])
        for bad_pos in [-1, l, l + 1]:
            self.assertRaises(_msprime.LibraryError, f, [(l, 0)])

    def test_free(self):
        ts = self.get_tree_sequence()
        t = _msprime.SparseTree(
            ts, flags=_msprime.LEAF_COUNTS | _msprime.LEAF_LISTS)
        no_arg_methods = [
            t.get_root, t.get_sample_size, t.get_index, t.get_left, t.get_right,
            t.get_num_sites, t.get_flags, t.get_sites, t.get_num_nodes]
        node_arg_methods = [
            t.get_parent, t.get_population, t.get_children, t.get_num_leaves,
            t.get_num_tracked_leaves]
        two_node_arg_methods = [t.get_mrca]
        for method in no_arg_methods:
            method()
        for method in node_arg_methods:
            method(0)
        for method in two_node_arg_methods:
            method(0, 0)
        t.free()
        self.assertRaises(RuntimeError, t.free)
        for method in no_arg_methods:
            self.assertRaises(RuntimeError, method)
        for method in node_arg_methods:
            self.assertRaises(RuntimeError, method, 0)
        for method in two_node_arg_methods:
            self.assertRaises(RuntimeError, method, 0, 0)


class TestTablesInterface(LowLevelTestCase):
    """
    Tests for the dump/load tables interface.
    """

    def verify_node_table(self, nodes, ts):
        """
        Verifies that the specified tree sequence and nodes table contain
        the same data.
        """
        self.assertEqual(nodes.num_rows, ts.get_num_nodes())
        self.assertGreater(nodes.num_rows, 0)
        # flags = nodes.flags
        population = nodes.population
        time = nodes.time
        for j in range(ts.get_num_nodes()):
            t = ts.get_node(j)
            # FIXME check flags
            # self.assertEqual(flags[j], t[0])
            self.assertEqual(time[j], t[1])
            self.assertEqual(population[j], t[2])

    def verify_edgeset_table(self, edgesets, ts):
        """
        Verifies that the specified tree sequence and edgesets table contain
        the same data.
        """
        self.assertEqual(edgesets.num_rows, ts.get_num_edgesets())
        self.assertGreater(edgesets.num_rows, 0)
        left = edgesets.left
        right = edgesets.right
        parent = edgesets.parent
        children = edgesets.children
        children_length = edgesets.children_length
        offset = 0
        for j in range(ts.get_num_edgesets()):
            t = ts.get_edgeset(j)
            self.assertEqual(t[0], left[j])
            self.assertEqual(t[1], right[j])
            self.assertEqual(t[2], parent[j])
            for child in t[3]:
                self.assertEqual(child, children[offset])
                offset += 1
            self.assertEqual(len(t[3]), children_length[j])

    def verify_migration_table(self, migrations, ts):
        """
        Verifies that the specified tree sequence and migrations table contain
        the same data.
        """
        self.assertEqual(migrations.num_rows, ts.get_num_migrations())
        self.assertGreater(migrations.num_rows, 0)
        left = migrations.left
        right = migrations.right
        node = migrations.node
        source = migrations.source
        dest = migrations.dest
        time = migrations.time
        for j in range(ts.get_num_migrations()):
            t = ts.get_migration(j)
            self.assertEqual(t[0], left[j])
            self.assertEqual(t[1], right[j])
            self.assertEqual(t[2], node[j])
            self.assertEqual(t[3], source[j])
            self.assertEqual(t[4], dest[j])
            self.assertEqual(t[5], time[j])

    def verify_site_table(self, sites, ts):
        """
        Verifies that the specified tree sequence and sites table contain
        the same data.
        """
        self.assertEqual(sites.num_rows, ts.get_num_sites())
        self.assertGreater(sites.num_rows, 0)
        position = sites.position
        ancestral_state = sites.ancestral_state
        for j in range(ts.get_num_sites()):
            t = ts.get_site(j)
            self.assertEqual(t[0], position[j])
            self.assertEqual(t[1], chr(ancestral_state[j]))

    def verify_mutation_table(self, mutations, ts):
        """
        Verifies that the specified tree sequence and mutations table contain
        the same data.
        """
        self.assertEqual(mutations.num_rows, ts.get_num_mutations())
        self.assertGreater(mutations.num_rows, 0)
        site = mutations.site
        node = mutations.node
        derived_state = mutations.derived_state
        derived_state_length = mutations.derived_state_length
        offset = 0
        for j in range(ts.get_num_mutations()):
            t = ts.get_mutation(j)
            self.assertEqual(t[0], site[j])
            self.assertEqual(t[1], node[j])
            self.assertEqual(len(t[2]), derived_state_length[j])
            for c in t[2]:
                self.assertEqual(ord(c), derived_state[offset])
                offset += 1

    def test_dump_tables(self):
        ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edgesets = _msprime.EdgesetTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()
        kwargs = {
            "nodes": nodes,
            "edgesets": edgesets,
            "sites": sites,
            "mutations": mutations,
            "migrations": migrations,
        }
        ts.dump_tables(**kwargs)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)
        self.verify_migration_table(migrations, ts)
        self.verify_mutation_table(mutations, ts)

    def test_dump_tables_errors(self):
        ts = self.get_example_migration_tree_sequence()
        kwargs = {
            "nodes": _msprime.NodeTable(),
            "edgesets": _msprime.EdgesetTable(),
            "sites": _msprime.SiteTable(),
            "mutations": _msprime.MutationTable(),
            "migrations": _msprime.MigrationTable()
        }

        for bad_type in [None, "", []]:
            for parameter in kwargs.keys():
                copy = dict(kwargs)
                copy[parameter] = bad_type
                self.assertRaises(TypeError, ts.dump_tables, **copy)
        self.assertRaises(TypeError, ts.dump_tables)
        self.assertRaises(TypeError, ts.dump_tables, nodes=_msprime.NodeTable())
        self.assertRaises(TypeError, ts.dump_tables, edgesets=_msprime.EdgesetTable())
        # Both mutations and mutation types must be specified together
        self.assertRaises(
            TypeError, ts.dump_tables, nodes=_msprime.NodeTable(),
            edgesets=_msprime.EdgesetTable(), mutations=_msprime.MutationTable())
        self.assertRaises(
            TypeError, ts.dump_tables, nodes=_msprime.NodeTable(),
            edgesets=_msprime.EdgesetTable(),
            sites=_msprime.SiteTable())

    def test_dump_tables_optional_parameters(self):
        ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edgesets = _msprime.EdgesetTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()

        ts.dump_tables(nodes=nodes, edgesets=edgesets)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)

        ts.dump_tables(nodes=nodes, edgesets=edgesets, migrations=migrations)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)
        self.verify_migration_table(migrations, ts)

        ts.dump_tables(
            nodes=nodes, edgesets=edgesets, sites=sites,
            mutations=mutations)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)
        self.verify_site_table(sites, ts)
        self.verify_mutation_table(mutations, ts)

    def test_load_tables(self):
        ex_ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edgesets = _msprime.EdgesetTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()
        provenance = [b"sdf", b"wer"]
        kwargs = {
            "nodes": nodes,
            "edgesets": edgesets,
            "sites": sites,
            "mutations": mutations,
            "migrations": migrations,
        }
        ex_ts.dump_tables(**kwargs)
        ts = _msprime.TreeSequence()
        kwargs["provenance_strings"] = provenance
        ts.load_tables(**kwargs)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)
        self.verify_migration_table(migrations, ts)
        self.verify_site_table(sites, ts)
        self.verify_mutation_table(mutations, ts)
        self.assertEqual(ts.get_provenance_strings(), provenance)

    def test_load_tables_errors(self):
        ex_ts = self.get_example_migration_tree_sequence()
        kwargs = {
            "nodes": _msprime.NodeTable(),
            "edgesets": _msprime.EdgesetTable(),
            "sites": _msprime.SiteTable(),
            "mutations": _msprime.MutationTable(),
            "migrations": _msprime.MigrationTable(),
        }
        ex_ts.dump_tables(**kwargs)
        ts = _msprime.TreeSequence()
        kwargs["provenance_strings"] = ["sdf", "wer"]

        for bad_type in [None, "", tuple()]:
            for parameter in kwargs.keys():
                copy = dict(kwargs)
                copy[parameter] = bad_type
                self.assertRaises(TypeError, ts.load_tables, **copy)
        self.assertRaises(TypeError, ts.load_tables)
        self.assertRaises(TypeError, ts.load_tables, nodes=_msprime.NodeTable())
        self.assertRaises(TypeError, ts.load_tables, edgesets=_msprime.EdgesetTable())
        # Both mutations and mutation types must be specified together
        self.assertRaises(
            TypeError, ts.load_tables, nodes=_msprime.NodeTable(),
            edgesets=_msprime.EdgesetTable(), mutations=_msprime.MutationTable())
        self.assertRaises(
            TypeError, ts.load_tables, nodes=_msprime.NodeTable(),
            edgesets=_msprime.EdgesetTable(),
            sites=_msprime.SiteTable())

    def test_load_tables_optional_parameters(self):
        ex_ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edgesets = _msprime.EdgesetTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()
        ex_ts.dump_tables(
            nodes=nodes, edgesets=edgesets, sites=sites,
            mutations=mutations, migrations=migrations)

        ts = _msprime.TreeSequence()
        ts.load_tables(nodes=nodes, edgesets=edgesets)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)
        self.assertEqual(ts.get_num_migrations(), 0)
        self.assertEqual(ts.get_num_mutations(), 0)

        ts = _msprime.TreeSequence()
        ts.load_tables(nodes=nodes, edgesets=edgesets, migrations=migrations)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)
        self.verify_migration_table(migrations, ts)
        self.assertEqual(ts.get_num_mutations(), 0)

        ts = _msprime.TreeSequence()
        ts.load_tables(
            nodes=nodes, edgesets=edgesets, sites=sites,
            mutations=mutations)
        self.verify_node_table(nodes, ts)
        self.verify_edgeset_table(edgesets, ts)
        self.verify_site_table(sites, ts)
        self.verify_mutation_table(mutations, ts)
        self.assertEqual(ts.get_num_migrations(), 0)

    def test_node_table_add_row(self):
        table = _msprime.NodeTable()
        table.add_row()
        self.assertEqual(table.num_rows, 1)
        self.assertEqual(table.population, [-1])
        self.assertEqual(table.flags, [0])
        self.assertEqual(table.time, [0])
        self.assertEqual(list(table.name), [])
        self.assertEqual(list(table.name_length), [0])

        name = "abcde"
        table = _msprime.NodeTable()
        table.add_row(flags=5, population=10, time=1.23, name=name)
        self.assertEqual(table.num_rows, 1)
        self.assertEqual(table.population, [10])
        self.assertEqual(table.flags, [5])
        self.assertEqual(table.time, [1.23])
        self.assertEqual(list(table.name), [ord(c) for c in name])
        self.assertEqual(list(table.name_length), [len(name)])

    def test_node_table_add_row_errors(self):
        table = _msprime.NodeTable()
        for bad_type in ["3423", None, []]:
            self.assertRaises(TypeError, table.add_row, flags=bad_type)
            self.assertRaises(TypeError, table.add_row, population=bad_type)
            self.assertRaises(TypeError, table.add_row, time=bad_type)
        for bad_type in [234, None, []]:
            self.assertRaises(TypeError, table.add_row, name=bad_type)

    def test_edgeset_table_add_row(self):
        table = _msprime.EdgesetTable()
        table = _msprime.EdgesetTable()
        table.add_row(left=0, right=1, parent=2, children=(0, 1))
        self.assertEqual(table.num_rows, 1)
        self.assertEqual(table.left, [0])
        self.assertEqual(table.right, [1])
        self.assertEqual(table.parent, [2])
        self.assertEqual(list(table.children), [0, 1])
        self.assertEqual(list(table.children_length), [2])

    def test_edgeset_table_add_row_errors(self):
        table = _msprime.EdgesetTable()
        for bad_type in ["3423", None, []]:
            self.assertRaises(
                TypeError, table.add_row, left=bad_type, right=1, parent=1,
                children=(0, 1))
            self.assertRaises(
                TypeError, table.add_row, left=0, right=bad_type, parent=1,
                children=(0, 1))
            self.assertRaises(
                TypeError, table.add_row, left=0, right=1, parent=bad_type,
                children=(0, 1))
            self.assertRaises(
                TypeError, table.add_row, left=0, right=1, parent=1, children=bad_type)
            self.assertRaises(
                TypeError, table.add_row, left=0, right=1, parent=1,
                children=(0, bad_type))
        for bad_type in [234, None, []]:
            self.assertRaises(TypeError, table.add_row, name=bad_type)

    def test_add_row_data(self):
        nodes = _msprime.NodeTable()
        edgesets = _msprime.EdgesetTable()
        mutations = _msprime.MutationTable()
        sites = _msprime.SiteTable()
        ts = self.get_example_migration_tree_sequence()
        ts.dump_tables(
            nodes=nodes, edgesets=edgesets, sites=sites,
            mutations=mutations)
        self.assertGreater(sites.num_rows, 0)
        self.assertGreater(mutations.num_rows, 0)

        new_nodes = _msprime.NodeTable()
        for j in range(ts.get_num_nodes()):
            flags, time, population, name = ts.get_node(j)
            new_nodes.add_row(
                flags=flags, time=time, population=population, name=name)
        self.assertEqual(list(nodes.time), list(new_nodes.time))
        self.assertEqual(list(nodes.flags), list(new_nodes.flags))
        self.assertEqual(list(nodes.population), list(new_nodes.population))
        self.assertEqual(list(nodes.name), list(new_nodes.name))

        new_edgesets = _msprime.EdgesetTable()
        for j in range(ts.get_num_edgesets()):
            left, right, parent, children = ts.get_edgeset(j)
            new_edgesets.add_row(left, right, parent, children)
        self.assertEqual(list(edgesets.left), list(new_edgesets.left))
        self.assertEqual(list(edgesets.right), list(new_edgesets.right))
        self.assertEqual(list(edgesets.parent), list(new_edgesets.parent))
        self.assertEqual(list(edgesets.children), list(new_edgesets.children))

        new_sites = _msprime.SiteTable()
        for j in range(ts.get_num_sites()):
            t = ts.get_site(j)
            new_sites.add_row(t[0], t[1])
        self.assertEqual(list(new_sites.position), list(sites.position))
        self.assertEqual(list(new_sites.ancestral_state), list(sites.ancestral_state))

        new_mutations = _msprime.MutationTable()
        for j in range(ts.get_num_mutations()):
            site, node, derived_state = ts.get_mutation(j)
            new_mutations.add_row(site, node, derived_state)
        self.assertEqual(list(new_mutations.site), list(mutations.site))
        self.assertEqual(list(new_mutations.node), list(mutations.node))
        self.assertEqual(
            list(new_mutations.derived_state), list(mutations.derived_state))


class TestLeafListIterator(LowLevelTestCase):
    """
    Tests for the low-level leaf list iterator.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.LeafListIterator)
        self.assertRaises(TypeError, _msprime.LeafListIterator, None)
        ts = self.get_tree_sequence()
        flags = _msprime.LEAF_COUNTS | _msprime.LEAF_LISTS
        tree = _msprime.SparseTree(ts, flags=flags)
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
        flags = _msprime.LEAF_COUNTS | _msprime.LEAF_LISTS
        tree = _msprime.SparseTree(ts, flags=flags)
        for tree in _msprime.SparseTreeIterator(tree):
            self.verify_iterator(_msprime.LeafListIterator(tree, 1))
            self.verify_iterator(
                _msprime.LeafListIterator(tree, tree.get_root()))

    def test_leaf_list(self):
        examples = [
            self.get_tree_sequence(sample_size=5, num_loci=10),
            self.get_tree_sequence(
                demographic_events=[
                    get_simple_bottleneck_event(0.2, proportion=1.0)])]
        flags = _msprime.LEAF_COUNTS | _msprime.LEAF_LISTS
        for ts in examples:
            st = _msprime.SparseTree(ts, flags=flags)
            for t in _msprime.SparseTreeIterator(st):
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
        self.assertRaises(ValueError, _msprime.RecombinationMap, 1, [0, 0.1], [1, 2, 3])
        self.assertRaises(ValueError, _msprime.RecombinationMap, 1, [], [0, 0])
        self.assertRaises(_msprime.LibraryError, _msprime.RecombinationMap, 1, [], [])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap, 1, [0, -2], [0, 0])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap, 1, [1, 0], [0, 0])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap, 1, [0, 1, 0.5], [0, 0, 0])
        self.assertRaises(
            _msprime.LibraryError, _msprime.RecombinationMap, 0, [0, 1], [0.1, 0])

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
        self.assertAlmostEqual(rm.genetic_to_physical(num_loci), 1)
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


class TestMutationGenerator(unittest.TestCase):
    """
    Tests for the mutation generator class.
    """
    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.MutationGenerator)
        rng = _msprime.RandomGenerator(1)
        for bad_type in ["x", {}, None]:
            self.assertRaises(TypeError, _msprime.MutationGenerator, rng, bad_type)
        for bad_value in [-1, -1e-3]:
            self.assertRaises(ValueError, _msprime.MutationGenerator, rng, bad_value)

    def test_mutation_rate(self):
        rng = _msprime.RandomGenerator(1)
        for rate in [0, 1e-12, 1e12, 1000, 0.01]:
            mutgen = _msprime.MutationGenerator(rng, rate)
            self.assertEqual(mutgen.get_mutation_rate(), rate)


class TestDemographyDebugger(unittest.TestCase):
    """
    Tests for the demography debugging interface.
    """
    def get_simulator(self, events):
        return _msprime.Simulator(
            get_samples(2), _msprime.RandomGenerator(1),
            demographic_events=events)

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


class TestLdCalculator(LowLevelTestCase):
    """
    Tests for the low-level leaf list iterator.
    """

    def get_array(self, buff, length=None):
        """
        Returns an array of double from a buffer.
        """
        code = b'd' if IS_PY2 else 'd'
        if length is None:
            return array.array(code, bytes(buff))
        else:
            return array.array(code, bytes(buff[:length * 8]))

    def get_buffer(self, num_values):
        """
        Returns a buffer suitable for storing the specified number
        of doubles.
        """
        return bytearray(8 * num_values)

    def test_empty_tree_sequence(self):
        ts = _msprime.TreeSequence()
        ldc = _msprime.LdCalculator(ts)
        buff = self.get_buffer(0)
        self.assertRaises(IndexError, ldc.get_r2_array, buff, 0)
        self.assertRaises(IndexError, ldc.get_r2, 0, 1)

    def test_constructor(self):
        for bad_type in [None, "1", []]:
            self.assertRaises(TypeError, _msprime.LdCalculator, bad_type)

    def test_refcounts(self):
        ts = self.get_tree_sequence()
        ldc = _msprime.LdCalculator(ts)
        buff = self.get_buffer(1)
        v = ldc.get_r2_array(dest=buff, source_index=0)
        self.assertEqual(v, 1)
        a = self.get_array(buff)
        # Delete the underlying tree sequence to ensure that nothing
        # nasty happens
        del ts
        buff = self.get_buffer(1)
        v = ldc.get_r2_array(dest=buff, source_index=0)
        self.assertEqual(v, 1)
        b = self.get_array(buff)
        self.assertEqual(a, b)

    def test_bad_mutation_indexes(self):
        ts = self.get_tree_sequence()
        m = ts.get_num_mutations()
        self.assertGreater(m, 0)
        bad_indexes = [-1, m, m + 1, m + 100]
        ldc = _msprime.LdCalculator(ts)
        for bad_index in bad_indexes:
            self.assertRaises(
                IndexError, ldc.get_r2_array, self.get_buffer(1), bad_index)
            self.assertRaises(IndexError, ldc.get_r2, bad_index, 0)
            self.assertRaises(IndexError, ldc.get_r2, 0, bad_index)

    def test_get_r2_interface(self):
        ts = self.get_tree_sequence()
        ldc = _msprime.LdCalculator(ts)
        self.assertRaises(TypeError, ldc.get_r2)
        self.assertRaises(TypeError, ldc.get_r2, 0)
        for bad_type in [None, "1", []]:
            self.assertRaises(TypeError, ldc.get_r2, 0, bad_type)
            self.assertRaises(TypeError, ldc.get_r2, bad_type, 0)

    def test_get_r2_array_interface(self):
        ts = self.get_tree_sequence()
        ldc = _msprime.LdCalculator(ts)
        for bad_type in [None, "1", []]:
            self.assertRaises(TypeError, ldc.get_r2_array, bad_type, 0)
            self.assertRaises(
                TypeError, ldc.get_r2_array, self.get_buffer(1), bad_type)
            self.assertRaises(
                TypeError, ldc.get_r2_array, self.get_buffer(1), 0,
                direction=bad_type)
            self.assertRaises(
                TypeError, ldc.get_r2_array, self.get_buffer(1), 0,
                max_mutations=bad_type)
            self.assertRaises(
                TypeError, ldc.get_r2_array, self.get_buffer(1), 0,
                max_distance=bad_type)
        buffers = [b'bytes', bytes()]
        for bad_buff in buffers:
            self.assertRaises(BufferError, ldc.get_r2_array, bad_buff, 0)
        for j in range(min(10, ts.get_num_mutations())):
            buff = self.get_buffer(j)
            # If we pass a buffer of a given size we should get back this
            # number of values.
            v = ldc.get_r2_array(buff, 0)
            self.assertEqual(v, j)
            v = ldc.get_r2_array(buff, 0, max_mutations=j)
            self.assertEqual(v, j)
            # If we set max_mutations to > size of the buffer this should
            # be an error.
            for k in range(1, 5):
                self.assertRaises(
                    BufferError, ldc.get_r2_array, buff, 0,
                    max_mutations=j + k)
        for bad_direction in [0, -2, 2, 10**6]:
            self.assertRaises(
                ValueError, ldc.get_r2_array, self.get_buffer(0), 0,
                direction=bad_direction)
        for bad_distance in [-1e-6, -1, -1e6]:
            self.assertRaises(
                ValueError, ldc.get_r2_array, self.get_buffer(0), 0,
                max_distance=bad_distance)

    def test_get_r2_array_from_new(self):
        ts = self.get_tree_sequence()
        self.assertGreater(ts.get_num_trees(), 1)
        self.assertGreater(ts.get_num_mutations(), 3)
        m = ts.get_num_mutations()
        buff = self.get_buffer(m)
        # We create a new instance of ldc each time to make sure we get
        # the correct behaviour on a new instance.
        ldc = _msprime.LdCalculator(ts)
        v = ldc.get_r2_array(buff, 0, direction=_msprime.FORWARD, max_mutations=1)
        self.assertEqual(v, 1)
        ldc = _msprime.LdCalculator(ts)
        v = ldc.get_r2_array(buff, 0, direction=_msprime.REVERSE, max_mutations=1)
        self.assertEqual(v, 0)
        ldc = _msprime.LdCalculator(ts)
        v = ldc.get_r2_array(buff, m - 1, direction=_msprime.FORWARD, max_mutations=1)
        self.assertEqual(v, 0)
        ldc = _msprime.LdCalculator(ts)
        v = ldc.get_r2_array(buff, m - 1, direction=_msprime.REVERSE, max_mutations=1)
        self.assertEqual(v, 1)
        ldc = _msprime.LdCalculator(ts)
        v = ldc.get_r2_array(buff, m // 2, direction=_msprime.FORWARD, max_mutations=1)
        self.assertEqual(v, 1)
        ldc = _msprime.LdCalculator(ts)
        v = ldc.get_r2_array(buff, m // 2, direction=_msprime.REVERSE, max_mutations=1)
        self.assertEqual(v, 1)

    def test_get_r2_array_random_seeks(self):
        num_start_positions = 100
        num_retries = 3
        ts = self.get_tree_sequence()
        self.assertGreater(ts.get_num_trees(), 1)
        self.assertGreater(ts.get_num_mutations(), 3)
        m = ts.get_num_mutations()
        buff = self.get_buffer(m)
        start_positions = [
            random.randint(0, m) for _ in range(num_start_positions)]
        directions = [
            random.choice([_msprime.FORWARD, _msprime.REVERSE])
            for _ in range(num_start_positions)]
        results = [[] for j in range(num_start_positions)]
        ldc = _msprime.LdCalculator(ts)
        for _ in range(num_retries):
            for j in range(num_start_positions):
                v = ldc.get_r2_array(
                    buff, start_positions[j], direction=directions[j],
                    max_mutations=10)
                self.assertLessEqual(v, 10)
                results[j].append(list(self.get_array(buff, v)))
        for result in results:
            self.assertEqual(len(result), num_retries)
            for rp in result[1:]:
                self.assertEqual(result[0], rp)
