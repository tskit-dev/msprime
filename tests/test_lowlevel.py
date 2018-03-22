# Copyright (C) 2015-2017 University of Oxford
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
        [sim.get_recombination_rate(), 0])


def get_simulation_model(name="hudson", population_size=0.25, **kwargs):
    """
    Returns simulation model dictionary suitable for passing to the low-level API.
    """
    d = {"name": name, "population_size": population_size}
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


def get_samples(num_samples):
    """
    Returns a sample list for the specified size.
    """
    t = (0, 0)
    return [t for _ in range(num_samples)]


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
        time=0.0, population=-1, initial_size=None, growth_rate=None):
    """
    Returns a population change event for the specified values.
    """
    ret = {
        "type": "population_parameters_change",
        "time": time,
        "population": population
    }
    if initial_size is not None:
        ret["initial_size"] = initial_size
    if growth_rate is not None:
        ret["growth_rate"] = growth_rate
    return ret


def get_size_change_event(time=0.0, size=1.0, population=-1):
    """
    Returns a size change demographic event.
    """
    return get_population_parameters_change_event(
        time, population, initial_size=size, growth_rate=0)


def get_growth_rate_change_event(time=0.0, growth_rate=1.0, population=-1):
    """
    Returns a growth_rate change demographic event.
    """
    return get_population_parameters_change_event(
        time, population, growth_rate=growth_rate)


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


def get_mass_migration_event(time=0.0, source=0, dest=1, proportion=1):
    """
    Returns a mass_migration demographic event.
    """
    return {
        "type": "mass_migration",
        "time": time,
        "source": source,
        "dest": dest,
        "proportion": proportion
    }


def get_simple_bottleneck_event(
        time=0.0, population=0, proportion=1):
    """
    Returns a simple bottleneck demographic event.
    """
    return {
        "type": "simple_bottleneck",
        "time": time,
        "population": population,
        "proportion": proportion
    }


def get_instantaneous_bottleneck_event(
        time=0.0, population=0, strength=10):
    """
    Returns a instantaneous bottleneck demographic event.
    """
    return {
        "type": "instantaneous_bottleneck",
        "time": time,
        "population": population,
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
        num_samples=10, Ne=0.25, random_seed=1, num_populations=1,
        store_migrations=False):
    samples = [(j % num_populations, 0) for j in range(num_samples)]
    migration_matrix = [1 for _ in range(num_populations**2)]
    for j in range(num_populations):
        migration_matrix[j * num_populations + j] = 0
    population_configuration = [
        get_population_configuration() for j in range(num_populations)]
    sim = _msprime.Simulator(
        samples, _msprime.RandomGenerator(random_seed),
        population_configuration=population_configuration,
        migration_matrix=migration_matrix,
        store_migrations=store_migrations,
        model=get_simulation_model(population_size=Ne))
    return sim


def populate_tree_sequence(sim, mutation_generator=None, provenances=[]):
    nodes = _msprime.NodeTable()
    edges = _msprime.EdgeTable()
    migrations = _msprime.MigrationTable()
    sites = _msprime.SiteTable()
    mutations = _msprime.MutationTable()
    provenance_table = _msprime.ProvenanceTable()
    ts = _msprime.TreeSequence()
    sim.populate_tables(nodes, edges, migrations)
    if mutation_generator is not None:
        mutation_generator.generate(nodes, edges, sites, mutations)
    for timestamp, record in provenances:
        provenance_table.add_row(timestamp=timestamp, record=record)
    ts.load_tables(nodes, edges, migrations, sites, mutations, provenance_table)
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
            population=random.randint(-1, num_populations - 1)))
        events.append(get_growth_rate_change_event(
            time=random.random(), growth_rate=random.random(),
            population=random.randint(-1, num_populations - 1)))
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
            dest = source
            while dest == source:
                dest = random.randint(0, num_populations - 1)
            # We use proportion of 0 or 1 so that we can test deterministically
            events.append(get_mass_migration_event(
                time=random.random(), proportion=random.choice([0, 1]),
                source=source, dest=dest))
            # Add some bottlenecks
            events.append(get_simple_bottleneck_event(
                time=random.random(), proportion=random.uniform(0, 0.25),
                population=random.randint(0, num_populations - 1)))
            events.append(get_instantaneous_bottleneck_event(
                time=random.random(), strength=random.uniform(0, 0.01),
                population=random.randint(0, num_populations - 1)))

    sorted_events = sorted(events, key=lambda x: x["time"])
    return sorted_events


def get_sample_counts(tree_sequence, st):
    """
    Returns a list of the sample node counts for the specfied sparse tree.
    """
    nu = [0 for j in range(st.get_num_nodes())]
    for j in range(tree_sequence.get_num_samples()):
        u = j
        while u != NULL_NODE:
            nu[u] += 1
            u = st.get_parent(u)
    return nu


def get_tracked_sample_counts(st, tracked_samples):
    """
    Returns a list giving the number of samples in the specified list
    that are in the subtree rooted at each node.
    """
    nu = [0 for j in range(st.get_num_nodes())]
    for j in tracked_samples:
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
    def verify_sparse_tree_dict(self, n, pi):
        """
        Verifies that the specified sparse tree in dict format is a
        consistent coalescent history for a sample of size n.
        """
        self.assertLessEqual(len(pi), 2 * n - 1)
        # NULL_NODE should not be a node
        self.assertNotIn(NULL_NODE, pi)
        # verify the root is equal for all samples
        root = 0
        while pi[root] != NULL_NODE:
            root = pi[root]
        for j in range(n):
            k = j
            while pi[k] != NULL_NODE:
                k = pi[k]
            self.assertEqual(k, root)
        # 0 to n - 1 inclusive should always be nodes
        for j in range(n):
            self.assertIn(j, pi)
        num_children = collections.defaultdict(int)
        for j in pi.keys():
            num_children[pi[j]] += 1
        # nodes 0 to n are samples.
        for j in range(n):
            self.assertNotEqual(pi[j], 0)
            self.assertEqual(num_children[j], 0)
        # All non-sample nodes should be binary
        for j in pi.keys():
            if j > n:
                self.assertGreaterEqual(num_children[j], 2)

    def get_tree_sequence(
            self, num_samples=10, num_loci=100, mutation_rate=10,
            random_seed=1, demographic_events=[], num_provenance_records=5):
        rho = 1.0
        rng = _msprime.RandomGenerator(random_seed)
        sim = _msprime.Simulator(
            get_samples(num_samples), rng, num_loci=num_loci,
            recombination_rate=rho,
            demographic_events=demographic_events)
        sim.run()
        nodes = _msprime.NodeTable()
        edges = _msprime.EdgeTable()
        migrations = _msprime.MigrationTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        provenances = _msprime.ProvenanceTable()
        ts = _msprime.TreeSequence()
        mutgen = _msprime.MutationGenerator(rng, mutation_rate)
        for j in range(num_provenance_records):
            provenances.add_row(timestamp="y" * j, record="x" * j)
        sim.populate_tables(nodes, edges, migrations)
        mutgen.generate(nodes, edges, sites, mutations)
        ts.load_tables(nodes, edges, migrations, sites, mutations, provenances)
        self.assertEqual(ts.get_num_nodes(), nodes.num_rows)
        self.assertEqual(ts.get_num_mutations(), mutations.num_rows)
        self.assertEqual(ts.get_num_edges(), edges.num_rows)
        self.assertEqual(ts.get_num_provenances(), provenances.num_rows)
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
        edges = _msprime.EdgeTable()
        migrations = _msprime.MigrationTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        provenances = _msprime.ProvenanceTable()
        for j in range(4):
            provenances.add_row(timestamp="y" * (j + 1), record="x" * j)
        ts = _msprime.TreeSequence()
        sim.populate_tables(nodes, edges, migrations)
        mutation_rate = 10
        mutgen = _msprime.MutationGenerator(_msprime.RandomGenerator(1), mutation_rate)
        mutgen.generate(nodes, edges, sites, mutations)
        ts.load_tables(nodes, edges, migrations, sites, mutations, provenances)
        self.assertGreater(sim.get_num_migrations(), 0)
        self.assertGreater(ts.get_num_mutations(), 0)
        return ts

    def verify_iterator(self, iterator):
        """
        Checks that the specified non-empty iterator implements the
        iterator protocol correctly.
        """
        list_ = list(iterator)
        self.assertGreater(len(list_), 0)
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
        self.assertGreater(sim.get_num_node_blocks(), 0)
        self.assertGreater(sim.get_num_edge_blocks(), 0)
        self.assertGreater(sim.get_num_migration_blocks(), 0)
        n = sim.get_num_samples()
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
        nodes = sim.get_nodes()
        self.assertEqual(len(nodes), sim.get_num_nodes())
        for j, (flags, t, pop, metadata) in enumerate(nodes):
            if j < sim.get_num_samples():
                self.assertEqual(t, 0.0)
                self.assertEqual(flags, 1)
            else:
                self.assertGreater(t, 0.0)
                self.assertEqual(flags, 0)
            self.assertTrue(0 <= pop < sim.get_num_populations())
            self.assertEqual(len(metadata), 0)

        edges = sim.get_edges()
        self.assertEqual(len(edges), sim.get_num_edges())
        for l, r, p, c in edges:
            self.assertTrue(0 <= l < m)
            self.assertTrue(1 <= r <= m)
            self.assertIn(l, breakpoints)
            self.assertIn(r, breakpoints)
        # The amount of ancestral material in the edges and
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
        for l, r, _, _ in edges:
            j = breakpoints.index(l)
            while breakpoints[j] < r:
                # We are assuming a binary coalescent here. We would need to
                # do something more complicated for more general models.
                records_am[j] += 0.5
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

    def verify_trees_equal(self, n, pi, sparse_tree):
        """
        Verifies that the specified parent map is equivalent to the specified
        sparse tree object.
        """
        pi_p = {}
        for j in range(n):
            u = j
            while u != NULL_NODE and u not in pi_p:
                pi_p[u] = sparse_tree.get_parent(u)
                u = pi_p[u]
        self.assertEqual(pi_p, pi)
        u = sparse_tree.get_left_root()
        while u != NULL_NODE:
            self.assertEqual(pi_p[u], NULL_NODE)
            u = sparse_tree.get_right_sib(u)

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

    def verify_trees(self, sim, sorted_edges):
        """
        Verifies that the specified set of (left, parent) sorted edges
        corresponds to correct trees for the specified simulation.
        """
        ts = populate_tree_sequence(sim)
        st = _msprime.SparseTree(ts)
        st_iter = _msprime.SparseTreeIterator(st)
        n = sim.get_num_samples()
        pi = {}
        last_l = 0
        num_trees = 0
        live_segments = []
        for l, r, parent, child in sorted_edges:
            if last_l != l:
                last_l = l
                for j in range(n):
                    assert j in pi
                # insert the root
                v = 0
                while v in pi:
                    v = pi[v]
                pi[v] = -1
                self.verify_sparse_tree_dict(n, pi)
                # Make sure this is equal to the sparse tree we get from
                # the iterator.
                st = next(st_iter)
                self.verify_trees_equal(n, pi, st)
                self.verify_mrcas(st)
                del pi[v]
                num_trees += 1
            heapq.heappush(live_segments, (r, (child, parent)))
            while live_segments[0][0] <= l:
                x, (other_child, p) = heapq.heappop(live_segments)
                del pi[other_child]
            pi[child] = parent
        for j in range(n):
            assert j in pi
        # Insert the root branch.
        v = 0
        while v in pi:
            v = pi[v]
        pi[v] = -1
        self.verify_sparse_tree_dict(n, pi)
        st = next(st_iter)
        self.verify_trees_equal(n, pi, st)
        num_trees += 1
        self.assertEqual(ts.get_num_trees(), num_trees)
        self.assertLessEqual(num_trees, sim.get_num_breakpoints() + 2)
        self.assertRaises(StopIteration, next, st_iter)

    def verify_squashed_edges(self, sorted_edges):
        """
        Checks to see if there were any unsquashed edges in the specified
        set of time sorted edges.
        """
        _, last_r, last_p, last_c = sorted_edges[0]
        for l, r, p, c in sorted_edges[1:]:
            if p == last_p and c == last_c:
                self.assertNotEqual(last_r, l)
            last_r, last_p, last_c = r, p, c

    def verify_completed_simulation(self, sim):
        """
        Verifies the state of the specified completed simulation.
        """
        self.assertEqual(sim.get_ancestors(), [])
        self.assertEqual(sim.get_num_ancestors(), 0)
        self.assertGreaterEqual(sim.get_num_breakpoints(), 0)
        self.assertGreater(sim.get_num_edges(), 0)
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
        self.assertGreater(sim.get_num_node_blocks(), 0)
        self.assertGreater(sim.get_num_edge_blocks(), 0)
        self.assertGreater(sim.get_num_migration_blocks(), 0)
        self.assertGreater(sim.get_used_memory(), 0)

        edges = sim.get_edges()
        self.assertGreater(len(edges), 0)
        self.assertEqual(len(edges), sim.get_num_edges())
        # Edges should be returned in canonical order
        self.assertEqual(edges, sorted(edges, key=lambda e: (e[2], e[3], e[0])))
        # Nodes should be in nondecreasing time order
        times = [node[1] for node in sim.get_nodes()]
        self.assertEqual(times, sorted(times))
        self.assertEqual(times[-1], sim.get_time())
        self.verify_squashed_edges(edges)

        left_sorted_edges = sorted(edges, key=lambda r: (r[0], r[2]))
        self.verify_trees(sim, left_sorted_edges)
        # Check the TreeSequence. Ensure we get the edges back in the
        # correct orders.
        ts = populate_tree_sequence(sim)
        ts_edges = [ts.get_edge(j) for j in range(ts.get_num_edges())]
        j = 0
        for edge in edges:
            left, right, parent, child = edge
            self.assertEqual(left, ts_edges[j][0])
            self.assertEqual(right, ts_edges[j][1])
            self.assertEqual(parent, ts_edges[j][2])
            self.assertEqual(child, ts_edges[j][3])
            j += 1
        assert j == len(ts_edges)

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
        num_sampless = [0 for j in range(N)]
        num_sampless[0] = n
        random_seed = random.randint(0, 2**31)
        max_memory = random.randint(10 * mb, 100 * mb)
        segment_block_size = random.randint(1, 100)
        node_mapping_block_size = random.randint(1, 100)
        avl_node_block_size = random.randint(1, 100)
        node_block_size = random.randint(1, 100)
        edge_block_size = random.randint(1, 100)
        migration_block_size = random.randint(1, 100)
        sim = _msprime.Simulator(
            samples=get_population_samples(*num_sampless),
            random_generator=_msprime.RandomGenerator(random_seed),
            num_loci=m,
            store_migrations=store_migrations,
            recombination_rate=rho,
            population_configuration=population_configuration,
            demographic_events=demographic_events,
            migration_matrix=migration_matrix,
            max_memory=max_memory,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size,
            node_block_size=node_block_size,
            edge_block_size=edge_block_size,
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
            self.assertGreater(sim.get_num_node_blocks(), 0)
            self.assertGreater(sim.get_num_edge_blocks(), 0)
            self.assertGreater(sim.get_num_migration_blocks(), 0)
            self.assertEqual(sim.get_num_samples(), n)
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
                self.assertEqual(n, sim.get_num_samples())
                self.assertEqual(m, sim.get_num_loci())
                self.assertEqual(rho, sim.get_recombination_rate())
                self.assertEqual(max_memory, sim.get_max_memory())
                self.assertEqual(segment_block_size, sim.get_segment_block_size())
                self.assertEqual(avl_node_block_size, sim.get_avl_node_block_size())
                self.assertEqual(
                    node_mapping_block_size, sim.get_node_mapping_block_size())
                self.assertEqual(node_block_size, sim.get_node_block_size())
                self.assertEqual(edge_block_size, sim.get_edge_block_size())
                self.assertEqual(migration_block_size, sim.get_migration_block_size())
                # Run this for a tiny amount of time and check the state
                self.assertFalse(sim.run(1e-8))
                self.verify_running_simulation(sim)
            sim.reset()

    def verify_tree_diffs(self, tree_sequence):
        n = tree_sequence.get_num_samples()
        L = tree_sequence.get_sequence_length()
        t = [tree_sequence.get_node(u)[1] for u in range(tree_sequence.get_num_nodes())]
        # Check some basic properties of the diffs.
        diffs = list(_msprime.TreeDiffIterator(tree_sequence))
        (left, right), edges_out, edges_in = diffs[0]
        self.assertEqual(left, 0)
        self.assertGreater(right, 0)
        self.assertEqual(len(edges_out), 0)
        self.assertLessEqual(len(edges_in), 2 * n - 2)
        last_right = right
        for (left, right), edges_out, edges_in in diffs[1:]:
            self.assertEqual(last_right, left)
            last_right = right
            self.assertGreater(right, left)
            self.assertLessEqual(right, L)
            for l, r, p, c in edges_in:
                self.assertEqual(l, left)
            for l, r, p, c in edges_out:
                self.assertEqual(r, left)
            # Make sure in edges are in increasing time order.
            time_sorted = sorted(edges_in, key=lambda x: t[x[2]])
            self.assertEqual(time_sorted, edges_in)
            # Make sure out edges are in decreasing time order.
            time_sorted = sorted(edges_out, key=lambda x: -t[x[2]])
            self.assertEqual(time_sorted, edges_out)
        # Compare with the Python implementation.
        pts = tests.PythonTreeSequence(tree_sequence)
        python_diffs = list(pts.edge_diffs())
        self.assertGreaterEqual(len(python_diffs), 0)
        self.assertEqual(len(diffs), len(python_diffs))
        for diff, py_diff in zip(diffs, python_diffs):
            self.assertEqual(diff[0], py_diff[0])
            self.assertEqual(len(diff[1]), len(py_diff[1]))
            self.assertEqual(len(diff[2]), len(py_diff[2]))
            for edge, py_edge in zip(diff[1], py_diff[1]):
                self.assertEqual(
                    edge, (py_edge.left, py_edge.right, py_edge.parent, py_edge.child))
            for edge, py_edge in zip(diff[2], py_diff[2]):
                self.assertEqual(
                    edge, (py_edge.left, py_edge.right, py_edge.parent, py_edge.child))

    def verify_sample_counts(self, tree_sequence):
        st = _msprime.SparseTree(
            tree_sequence, flags=_msprime.SAMPLE_COUNTS, tracked_samples=[])
        for _ in _msprime.SparseTreeIterator(st):
            self.assertEqual(st.get_flags(), _msprime.SAMPLE_COUNTS)
            nu = get_sample_counts(tree_sequence, st)
            nu_prime = [
                st.get_num_samples(j) for j in
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
            recombination_rate=r,
            random_generator=_msprime.RandomGenerator(random_seed),
            demographic_events=demographic_events,
            max_memory=10 * mb, segment_block_size=1000,
            avl_node_block_size=1000, node_mapping_block_size=1000,
            node_block_size=1000, edge_block_size=1000, model=model)
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
            self.verify_sample_counts(tree_sequence)
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
            samples=get_samples(n), num_loci=m, recombination_rate=1,
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
        self.assertRaises(_msprime.LibraryError, sim.compute_population_size, N + 1, 0)
        # For each event we now run the simulator forward until this time
        # and make sure that the internal state is what it should be.
        for event in demographic_events:
            t = event["time"]
            event_type = event["type"]
            self.assertEqual(next_event_time, t)
            self.assertEqual(sim2.get_migration_matrix(), migration_matrix)
            completed = sim.run(t)
            for j in range(N):
                s = sim2.compute_population_size(j, t)
                self.assertGreater(s, 0)
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
                population = event["population"]
                indexes = [population]
                if population == -1:
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

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
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
            self.assertRaises(TypeError, f, recombination_rate=bad_type)
            self.assertRaises(TypeError, f, max_memory=bad_type)
            self.assertRaises(TypeError, f, avl_node_block_size=bad_type)
            self.assertRaises(TypeError, f, segment_block_size=bad_type)
            self.assertRaises(TypeError, f, node_mapping_block_size=bad_type)
            self.assertRaises(TypeError, f, node_block_size=bad_type)
            self.assertRaises(TypeError, f, edge_block_size=bad_type)
        # Check for bad values.
        self.assertRaises(_msprime.InputError, f, num_loci=0)
        self.assertRaises(_msprime.InputError, f, recombination_rate=-1)
        self.assertRaises(_msprime.InputError, f, max_memory=0)
        self.assertRaises(_msprime.InputError, f, avl_node_block_size=0)
        self.assertRaises(_msprime.InputError, f, segment_block_size=0)
        self.assertRaises(_msprime.InputError, f, node_mapping_block_size=0)
        self.assertRaises(_msprime.InputError, f, node_block_size=0)
        self.assertRaises(_msprime.InputError, f, edge_block_size=0)
        # Check for other type specific errors.
        self.assertRaises(OverflowError, f, max_memory=2**65)

    def test_non_parametric_simulation_models(self):

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
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

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
                _msprime.RandomGenerator(random_seed), **kwargs)
        for bad_type in [None, str, "sdf"]:
            model = get_simulation_model("dirac", psi=bad_type, c=1.0)
            self.assertRaises(TypeError, f, model=model)
            model = get_simulation_model("dirac", psi=0.5, c=bad_type)
            self.assertRaises(TypeError, f, model=model)
        self.assertRaises(ValueError, f, model=get_simulation_model("dirac"))
        for bad_psi in [-1, 0, -1e-6, 1, 1e6]:
            self.assertRaises(
                ValueError, f, model=get_simulation_model("dirac", c=1, psi=bad_psi))
        for bad_c in [-1, -1e-6]:
            self.assertRaises(
                ValueError, f, model=get_simulation_model("dirac", psi=0.5, c=bad_c))
        for psi in [0.99, 0.2, 1e-4]:
            for c in [5.0, 1e2, 1e-4]:
                model = get_simulation_model("dirac", psi=psi, c=c)
                sim = f(model=model)
                self.assertEqual(sim.get_model(), model)

    def test_beta_simulation_model(self):

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
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
                    "beta", 2, alpha=alpha, truncation_point=truncation_point)
                sim = f(model=model)
                self.assertEqual(sim.get_model(), model)

    def test_store_migrations(self):
        def f(num_samples=10, random_seed=1, **kwargs):
            samples = [(j % 2, 0) for j in range(num_samples)]
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
            [get_population_configuration(initial_size=0.25)])

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
        for bad_size in [-1, 0]:
            self.assertRaises(
                _msprime.InputError, f,
                [get_population_configuration(initial_size=bad_size)])

    def test_bad_sample_configurations(self):
        rng = _msprime.RandomGenerator(1)

        def f(num_samples, pop_num_sampless):
            population_configuration = [
                get_population_configuration(n) for n in pop_num_sampless]
            N = len(pop_num_sampless)
            migration_matrix = [0 for j in range(N) for k in range(N)]
            _msprime.Simulator(
                get_samples(num_samples), rng,
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
        def f(num_samples, conf_tuples):
            population_configuration = [
                get_population_configuration(initial_size=p, growth_rate=a)
                for p, a in conf_tuples]
            N = len(population_configuration)
            migration_matrix = [0 for j in range(N) for k in range(N)]
            s = _msprime.Simulator(
                get_samples(num_samples), _msprime.RandomGenerator(1),
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
            population_configuration[0]["num_samples"] = 2
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
            event = get_size_change_event(population=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_size_change_event(size=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_instantaneous_bottleneck_event(population=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_instantaneous_bottleneck_event(strength=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_simple_bottleneck_event(population=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_simple_bottleneck_event(proportion=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_mass_migration_event(source=bad_type, dest=0)
            self.assertRaises(TypeError, f, [event])
            event = get_mass_migration_event(source=0, dest=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_mass_migration_event(
                source=0, dest=1, proportion=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_migration_rate_change_event(matrix_index=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_migration_rate_change_event(
                matrix_index=0, migration_rate=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_growth_rate_change_event(population=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_growth_rate_change_event(population=0, growth_rate=bad_type)
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
                    population=bad_pop_id)
                self.assertRaises(_msprime.InputError, f, [event])
            for k in range(1, 4):
                event = event_generator(population=k)
                self.assertRaises(_msprime.InputError, f, [event], k)
                events = [event_generator(), event]
                self.assertRaises(_msprime.InputError, f, events, k)
        for bad_pop_id in [-2, 1, 10**6]:
            event = get_mass_migration_event(source=bad_pop_id)
            self.assertRaises(_msprime.InputError, f, [event])
            event = get_mass_migration_event(dest=bad_pop_id)
            self.assertRaises(_msprime.InputError, f, [event])
            event = get_simple_bottleneck_event(population=bad_pop_id)
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
        event = get_mass_migration_event(source=0, dest=0)
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
                "recombination_rate": 0.1,
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
            self.assertEqual(sim1.get_num_edges(), sim2.get_num_edges())
            self.assertEqual(sim1.get_edges(), sim2.get_edges())
            self.assertEqual(sim1.get_num_nodes(), sim2.get_num_nodes())
            self.assertEqual(sim1.get_nodes(), sim2.get_nodes())
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
                recombination_rate=0)
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
        num_sampless = [0 for _ in range(num_populations)]
        num_sampless[active_pops[0]] = 2
        num_sampless[active_pops[1]] = 2
        num_sampless[active_pops[2]] = 6
        flattened = [x for row in migration_matrix for x in row]
        sim = _msprime.Simulator(
            get_population_samples(*num_sampless),
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
                    t + dt, source=0, dest=1, proportion=1)],
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
                get_simple_bottleneck_event(t1, population=0, proportion=1),
                get_simple_bottleneck_event(t2, population=1, proportion=1),
                get_mass_migration_event(t3, source=0, dest=1, proportion=1)],
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
            recombination_rate=10)
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
        edge_table = _msprime.EdgeTable()
        migration_table = _msprime.MigrationTable()
        sim = _msprime.Simulator(get_samples(10), _msprime.RandomGenerator(1))
        sim.run()

        recomb_map = uniform_recombination_map(sim)
        # nodes, edges and migrations are mandatory.
        self.assertRaises(TypeError, sim.populate_tables)
        self.assertRaises(TypeError, sim.populate_tables, nodes=node_table)
        self.assertRaises(
            TypeError, sim.populate_tables, nodes=node_table, migrations=migration_table)
        self.assertRaises(
            TypeError, sim.populate_tables, nodes=node_table, edges=edge_table)

        kwargs = {
            "nodes": node_table,
            "edges": edge_table,
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
        self.assertEqual(edge_table.num_rows, sim.get_num_edges())
        kwargs["recombination_map"] = recomb_map
        sim.populate_tables(**kwargs)
        self.assertEqual(edge_table.num_rows, sim.get_num_edges())


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
                "recombination_rate": 0.1,
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
            self.assertEqual(t1.get_num_edges(), t2.get_num_edges())
            r1 = [t1.get_edge(j) for j in range(t1.get_num_edges())]
            r2 = [t2.get_edge(j) for j in range(t2.get_num_edges())]
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
        self.assertRaises(ValueError, ts.get_num_samples)
        self.assertRaises(ValueError, ts.get_sequence_length)
        self.assertRaises(ValueError, ts.get_num_trees)
        self.assertRaises(ValueError, ts.get_num_edges)
        self.assertRaises(ValueError, ts.get_num_mutations)
        self.assertRaises(ValueError, ts.get_num_migrations)
        self.assertRaises(ValueError, ts.get_num_migrations)
        self.assertRaises(ValueError, ts.get_genotype_matrix)
        self.assertRaises(ValueError, ts.dump)

    def test_num_nodes(self):
        for ts in self.get_example_tree_sequences():
            max_node = 0
            for j in range(ts.get_num_edges()):
                _, _, parent, child = ts.get_edge(j)
                for node in [parent, child]:
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
        self.assertEqual(ts.get_num_samples(), ts2.get_num_samples())
        self.assertEqual(ts.get_sequence_length(), ts2.get_sequence_length())
        self.assertEqual(ts.get_num_mutations(), ts2.get_num_mutations())
        self.assertEqual(ts.get_num_nodes(), ts2.get_num_nodes())
        records1 = [ts.get_edge(j) for j in range(ts.get_num_edges())]
        records2 = [ts2.get_edge(j) for j in range(ts2.get_num_edges())]
        self.assertEqual(records1, records2)
        mutations1 = [ts.get_mutation(j) for j in range(ts.get_num_mutations())]
        mutations2 = [ts2.get_mutation(j) for j in range(ts2.get_num_mutations())]
        self.assertEqual(mutations1, mutations2)
        provenances1 = [ts.get_provenance(j) for j in range(ts.get_num_provenances())]
        provenances2 = [ts2.get_provenance(j) for j in range(ts2.get_num_provenances())]
        self.assertEqual(provenances1, provenances2)

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

    def test_get_edge_interface(self):
        for ts in self.get_example_tree_sequences():
            num_edges = ts.get_num_edges()
            # We don't accept Python negative indexes here.
            self.assertRaises(IndexError, ts.get_edge, -1)
            for j in [0, 10, 10**6]:
                self.assertRaises(IndexError, ts.get_edge, num_edges + j)
            for x in [None, "", {}, []]:
                self.assertRaises(TypeError, ts.get_edge, x)

    def test_get_node_interface(self):
        for ts in self.get_example_tree_sequences():
            num_nodes = ts.get_num_nodes()
            # We don't accept Python negative indexes here.
            self.assertRaises(IndexError, ts.get_node, -1)
            for j in [0, 10, 10**6]:
                self.assertRaises(IndexError, ts.get_node, num_nodes + j)
            for x in [None, "", {}, []]:
                self.assertRaises(TypeError, ts.get_node, x)

    def test_get_genotype_matrix_interface(self):
        for ts in self.get_example_tree_sequences():
            num_samples = ts.get_num_samples()
            num_sites = ts.get_num_sites()
            G = ts.get_genotype_matrix()
            self.assertEqual(G.shape, (num_sites, num_samples))

    def test_get_migration_interface(self):
        ts = self.get_example_migration_tree_sequence()
        for bad_type in ["", None, {}]:
            self.assertRaises(TypeError, ts.get_migration, bad_type)
        num_records = ts.get_num_migrations()
        # We don't accept Python negative indexes here.
        self.assertRaises(IndexError, ts.get_migration, -1)
        for j in [0, 10, 10**6]:
            self.assertRaises(IndexError, ts.get_migration, num_records + j)

    def test_get_samples(self):
        ts = self.get_example_migration_tree_sequence()
        # get_samples takes no arguments.
        self.assertRaises(TypeError, ts.get_samples, 0)
        self.assertEqual(list(range(ts.get_num_samples())), ts.get_samples())

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
        for Ne in [0.25, 1, 10, 1e6]:
            sim = get_example_simulator(
                10, Ne=Ne, num_populations=2, store_migrations=True)
            sim.run()
            tree_sequence = populate_tree_sequence(sim)
            sim_times = [node[1] for node in sim.get_nodes()]
            for j in range(tree_sequence.get_num_nodes()):
                _, generation, _, _ = tree_sequence.get_node(j)
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
                [0, ts.get_num_samples()])
            self.assertRaises(
                _msprime.LibraryError, ts.get_pairwise_diversity, [0, 0])
            samples = list(range(ts.get_num_samples()))
            pi1 = ts.get_pairwise_diversity(samples)
            self.assertGreaterEqual(pi1, 0)

    def test_provenance_populate(self):
        rng = _msprime.RandomGenerator(1)
        sim = _msprime.Simulator(get_samples(10), rng)
        sim.run()
        for j in range(10):
            pr = []
            for k in range(j):
                # Use some fairly big strings to stress things out
                timestamp = "x" * (k + 1) * 8192
                record = "y" * (k + 1) * 8192
                pr.append((timestamp, record))
            ts = populate_tree_sequence(sim, provenances=pr)
            records = [ts.get_provenance(j) for j in range(ts.get_num_provenances())]
            self.assertEqual(records, pr)

#             ts.dump(self.temp_file)
#             ts2 = _msprime.TreeSequence()
#             ts2.load(self.temp_file)
#             self.assertEqual(ts2.get_provenance_strings(), strings)


class TestVcfConverter(LowLevelTestCase):
    """
    Tests for the low-level vcf converter.
    """
    def test_uninitialised_tree_sequence(self):
        ts = _msprime.TreeSequence()
        self.assertRaises(ValueError, _msprime.VcfConverter, ts)

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
        for bad_type in [None, 0, [], {}]:
            self.assertRaises(
                TypeError, _msprime.VcfConverter, ts, contig_id=bad_type)
        converter = _msprime.VcfConverter(ts)
        before = converter.get_header() + "".join(converter)
        self.assertGreater(len(before), 0)
        converter = _msprime.VcfConverter(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = converter.get_header() + "".join(converter)
        self.assertEqual(before, after)

    def test_ploidy(self):
        ts = self.get_tree_sequence(num_samples=10)
        for bad_ploidy in [-1, 3, 4, 11, 20, 10**6]:
            self.assertRaises(
                _msprime.LibraryError, _msprime.VcfConverter, ts, ploidy=bad_ploidy)
        # 0 is caught as a ValueError before getting to the library code.
        self.assertRaises(ValueError, _msprime.VcfConverter, ts, ploidy=0)

    def test_contig_id(self):
        ts = self.get_tree_sequence(num_samples=10)
        self.assertRaises(ValueError, _msprime.VcfConverter, ts, contig_id="")
        for contig_id in ["1", "chr1", 8192 * "x"]:
            converter = _msprime.VcfConverter(ts, contig_id=contig_id)
            for row in converter:
                self.assertTrue(row.startswith(contig_id))
        # Try to be nasty and break things.
        s = 1024 * "1234"
        converter = _msprime.VcfConverter(ts, contig_id=s)
        del s
        del ts
        for row in converter:
            self.assertTrue(row.startswith("1234"))

    def test_positions(self):
        # Make sure we successfully discretise the positions when
        # we have lots of mutations
        ts = self.get_tree_sequence(num_samples=10, mutation_rate=10)
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
    def test_uninitialised_tree_sequence(self):
        ts = _msprime.TreeSequence()
        self.assertRaises(ValueError, _msprime.TreeDiffIterator, ts)

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
    def test_uninitialised_tree_sequence(self):
        ts = _msprime.TreeSequence()
        self.assertRaises(ValueError, _msprime.SparseTree, ts)

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
            "recombination_rate": 0.192184324155680
        }
        sim = _msprime.Simulator(**params)
        sim.run()
        ts = populate_tree_sequence(sim)
        st = _msprime.SparseTree(ts)
        for st in _msprime.SparseTreeIterator(st):
            root = 0
            while st.get_parent(root) != NULL_NODE:
                root = st.get_parent(root)
            self.assertEqual(st.get_left_root(), root)


class TestHaplotypeGenerator(LowLevelTestCase):
    """
    Tests for the low-level haplotype generator.
    """
    def test_uninitialised_tree_sequence(self):
        ts = _msprime.TreeSequence()
        self.assertRaises(ValueError, _msprime.HaplotypeGenerator, ts)

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.HaplotypeGenerator)
        ts = self.get_tree_sequence(num_loci=10)
        for bad_type in ["", {}, [], None]:
            self.assertRaises(TypeError, _msprime.HaplotypeGenerator, ts, bad_type)
        n = ts.get_num_samples()
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
    def test_uninitialised_tree_sequence(self):
        ts = _msprime.TreeSequence()
        self.assertRaises(ValueError, _msprime.VariantGenerator, ts)

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.VariantGenerator)
        ts = self.get_tree_sequence(num_loci=10)
        for bad_type in ["", {}, [], None]:
            self.assertRaises(TypeError, _msprime.VariantGenerator, bad_type)
            self.assertRaises(TypeError, _msprime.VariantGenerator, ts, bad_type)

        vg = _msprime.VariantGenerator(ts)
        before = [
            (site, genotypes.tobytes(), alleles) for site, genotypes, alleles in vg]
        vg = _msprime.VariantGenerator(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = [
            (site, genotypes.tobytes(), alleles) for site, genotypes, alleles in vg]
        self.assertEqual(before, after)

    def test_form(self):
        ts = self.get_tree_sequence(num_loci=10)
        variants = list(_msprime.VariantGenerator(ts))
        self.assertGreater(len(variants), 0)
        self.assertEqual(len(variants), ts.get_num_sites())
        sites = list(range(ts.get_num_sites()))
        self.assertEqual([site for site, _, _ in variants], sites)
        for _, genotypes, alleles in _msprime.VariantGenerator(ts):
            self.assertEqual(genotypes.shape, (ts.get_num_samples(), ))
            self.assertEqual(alleles, ('0', '1'))

    def test_iterator(self):
        ts = self.get_tree_sequence()
        variants = _msprime.VariantGenerator(ts)
        self.verify_iterator(variants)


class TestSparseTree(LowLevelTestCase):
    """
    Tests on the low-level sparse tree interface.
    """

    def test_flags(self):
        ts = self.get_tree_sequence()
        st = _msprime.SparseTree(ts)
        self.assertEqual(st.get_flags(), 0)
        # We should still be able to count the samples, just inefficiently.
        self.assertEqual(st.get_num_samples(0), 1)
        self.assertRaises(_msprime.LibraryError, st.get_num_tracked_samples, 0)
        self.assertRaises(_msprime.LibraryError, _msprime.SampleListIterator, st, 0)
        all_flags = [
            0, _msprime.SAMPLE_COUNTS, _msprime.SAMPLE_LISTS,
            _msprime.SAMPLE_COUNTS | _msprime.SAMPLE_LISTS]
        for flags in all_flags:
            st = _msprime.SparseTree(ts, flags=flags)
            self.assertEqual(st.get_flags(), flags)
            self.assertEqual(st.get_num_samples(0), 1)
            if flags & _msprime.SAMPLE_COUNTS:
                self.assertEqual(st.get_num_tracked_samples(0), 0)
            else:
                self.assertRaises(_msprime.LibraryError, st.get_num_tracked_samples, 0)
            if flags & _msprime.SAMPLE_LISTS:
                samples = list(_msprime.SampleListIterator(st, 0))
                self.assertEqual(samples, [0])
            else:
                self.assertRaises(
                    _msprime.LibraryError, _msprime.SampleListIterator, st, 0)

    def test_sites(self):
        for ts in self.get_example_tree_sequences():
            st = _msprime.SparseTree(ts)
            all_sites = [ts.get_site(j) for j in range(ts.get_num_sites())]
            all_tree_sites = []
            j = 0
            mutation_id = 0
            for st in _msprime.SparseTreeIterator(st):
                tree_sites = st.get_sites()
                self.assertEqual(st.get_num_sites(), len(tree_sites))
                all_tree_sites.extend(tree_sites)
                for position, ancestral_state, mutations, index, metadata in tree_sites:
                    self.assertTrue(st.get_left() <= position < st.get_right())
                    self.assertEqual(index, j)
                    self.assertEqual(metadata, b"")
                    for mut_id in mutations:
                        site, node, derived_state, parent, metadata = \
                            ts.get_mutation(mut_id)
                        self.assertEqual(site, index)
                        self.assertEqual(mutation_id, mut_id)
                        self.assertNotEqual(st.get_parent(node), NULL_NODE)
                        self.assertEqual(metadata, b"")
                        mutation_id += 1
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
                TypeError, _msprime.SparseTree, ts, tracked_samples=bad_type)
        for bad_type in ["", {}, None, []]:
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, flags=bad_type)
        for n in range(1, 10):
            ts = self.get_tree_sequence(num_samples=10, num_loci=n)
            st = _msprime.SparseTree(ts)
            self.assertEqual(st.get_num_nodes(), ts.get_num_nodes())
            # An uninitialised sparse tree should always be zero.
            self.assertEqual(st.get_left_root(), 0)
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
            ts = self.get_tree_sequence(num_samples=100, num_loci=n)
            num_nodes = ts.get_num_nodes()
            st = _msprime.SparseTree(ts)
            # deleting the tree sequence should still give a well formed
            # sparse tree.
            st_iter = _msprime.SparseTreeIterator(st)
            next(st_iter)
            del ts
            del st_iter
            # Do a quick traversal just to exercise the tree
            stack = [st.get_left_root()]
            while len(stack) > 0:
                u = stack.pop()
                self.assertLess(u, num_nodes)
                stack.extend(st.get_children(u))

    def test_bad_tracked_samples(self):
        ts = self.get_tree_sequence()
        flags = _msprime.SAMPLE_COUNTS
        for bad_type in ["", {}, [], None]:
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, flags=flags,
                tracked_samples=[bad_type])
            self.assertRaises(
                TypeError, _msprime.SparseTree, ts, flags=flags,
                tracked_samples=[1, bad_type])
        for bad_sample in [10**6, -1e6]:
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts, flags=flags,
                tracked_samples=[bad_sample])
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts, flags=flags,
                tracked_samples=[1, bad_sample])
            self.assertRaises(
                ValueError, _msprime.SparseTree, ts,
                tracked_samples=[1, bad_sample, 1])

    def test_count_all_samples(self):
        for ts in self.get_example_tree_sequences():
            self.verify_iterator(_msprime.TreeDiffIterator(ts))
            st = _msprime.SparseTree(ts, flags=_msprime.SAMPLE_COUNTS)
            # Without initialisation we should be 0 samples for every node
            # that is not a sample.
            for j in range(st.get_num_nodes()):
                count = 1 if j < ts.get_num_samples() else 0
                self.assertEqual(st.get_num_samples(j), count)
                self.assertEqual(st.get_num_tracked_samples(j), 0)
            # Now, try this for a tree sequence.
            for st in _msprime.SparseTreeIterator(st):
                nu = get_sample_counts(ts, st)
                nu_prime = [
                    st.get_num_samples(j) for j in
                    range(st.get_num_nodes())]
                self.assertEqual(nu, nu_prime)
                # For tracked samples, this should be all zeros.
                nu = [
                    st.get_num_tracked_samples(j) for j in
                    range(st.get_num_nodes())]
                self.assertEqual(nu, list([0 for _ in nu]))

    def test_count_tracked_samples(self):
        examples = [
            self.get_tree_sequence(num_samples=5, num_loci=10),
            self.get_tree_sequence(
                num_samples=12,
                demographic_events=[
                    get_simple_bottleneck_event(0.2, proportion=1.0)])]
        # Ensure that there are some non-binary nodes.
        non_binary = False
        for ts in examples:
            st = _msprime.SparseTree(ts)
            for st in _msprime.SparseTreeIterator(st):
                for u in range(ts.get_num_nodes()):
                    if len(st.get_children(u)) > 1:
                        non_binary = True
            samples = [j for j in range(ts.get_num_samples())]
            powerset = itertools.chain.from_iterable(
                itertools.combinations(samples, r)
                for r in range(len(samples) + 1))
            for subset in map(list, powerset):
                # Ordering shouldn't make any different.
                random.shuffle(subset)
                st = _msprime.SparseTree(
                    ts, flags=_msprime.SAMPLE_COUNTS, tracked_samples=subset)
                for st in _msprime.SparseTreeIterator(st):
                    nu = get_tracked_sample_counts(st, subset)
                    nu_prime = [
                        st.get_num_tracked_samples(j) for j in
                        range(st.get_num_nodes())]
                    self.assertEqual(nu, nu_prime)
            # Passing duplicated values should raise an error
            sample = 1
            for j in range(2, 20):
                tracked_samples = [sample for _ in range(j)]
                self.assertRaises(
                    _msprime.LibraryError, _msprime.SparseTree,
                    ts, flags=_msprime.SAMPLE_COUNTS,
                    tracked_samples=tracked_samples)
        self.assertTrue(non_binary)

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
                    ValueError, _msprime.SampleListIterator, st, v)

    def test_mrca_interface(self):
        for num_loci in range(1, 10):
            for num_samples in range(2, 5):
                ts = self.get_tree_sequence(
                    num_loci=num_loci, num_samples=num_samples)
                num_nodes = ts.get_num_nodes()
                st = _msprime.SparseTree(ts)
                for v in [num_nodes, 10**6, NULL_NODE]:
                    self.assertRaises(ValueError, st.get_mrca, v, v)
                    self.assertRaises(ValueError, st.get_mrca, v, 1)
                    self.assertRaises(ValueError, st.get_mrca, 1, v)
                # All the mrcas for an uninitialised tree should be NULL_NODE
                for u, v in itertools.combinations(range(num_nodes), 2):
                    self.assertEqual(st.get_mrca(u, v), NULL_NODE)

    def test_newick_precision(self):

        def get_times(tree):
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

        ts = self.get_tree_sequence(num_samples=10, num_loci=5)
        st = _msprime.SparseTree(ts)
        for st in _msprime.SparseTreeIterator(st):
            self.assertRaises(ValueError, st.get_newick, precision=-1)
            self.assertRaises(ValueError, st.get_newick, precision=17)
            self.assertRaises(ValueError, st.get_newick, precision=100)
            for precision in range(17):
                tree = st.get_newick(precision=precision).decode()
                times = get_times(tree)
                self.assertGreater(len(times), ts.get_num_samples())
                for t in times:
                    if precision == 0:
                        self.assertNotIn(".", t)
                    else:
                        point = t.find(".")
                        self.assertEqual(precision, len(t) - point - 1)

    @unittest.skip("Correct initialisation for sparse tree.")
    def test_newick_interface(self):
        ts = self.get_tree_sequence(num_loci=10, num_samples=10)
        st = _msprime.SparseTree(ts)
        # TODO this will break when we correctly handle multiple roots.
        self.assertEqual(st.get_newick(), b"1;")
        for bad_type in [None, "", [], {}]:
            self.assertRaises(TypeError, st.get_newick, precision=bad_type)
            self.assertRaises(TypeError, st.get_newick, ts, time_scale=bad_type)
        for st in _msprime.SparseTreeIterator(st):
            newick = st.get_newick()
            self.assertTrue(newick.endswith(b";"))

    def test_index(self):
        for num_loci in [1, 100]:
            ts = self.get_tree_sequence(num_loci=num_loci, num_samples=10)
            st = _msprime.SparseTree(ts)
            for index, st in enumerate(_msprime.SparseTreeIterator(st)):
                self.assertEqual(index, st.get_index())

    def test_bad_mutations(self):
        ts = self.get_tree_sequence(2, num_loci=200, random_seed=1)
        node_table = _msprime.NodeTable()
        edge_table = _msprime.EdgeTable()
        site_table = _msprime.SiteTable()
        mutation_table = _msprime.MutationTable()
        migration_table = _msprime.MigrationTable()
        ts.dump_tables(nodes=node_table, edges=edge_table)

        def f(mutations):
            position = []
            node = []
            site = []
            ancestral_state = []
            ancestral_state_offset = [0]
            derived_state = []
            derived_state_offset = [0]
            for j, (p, n) in enumerate(mutations):
                site.append(j)
                position.append(p)
                ancestral_state.append("0")
                ancestral_state_offset.append(ancestral_state_offset[-1] + 1)
                derived_state.append("1")
                derived_state_offset.append(derived_state_offset[-1] + 1)
                node.append(n)
            site_table.set_columns(
                position=position, ancestral_state=ancestral_state,
                ancestral_state_offset=ancestral_state_offset)
            mutation_table.set_columns(
                site=site, node=node, derived_state=derived_state,
                derived_state_offset=derived_state_offset)
            ts2 = _msprime.TreeSequence()
            ts2.load_tables(
                nodes=node_table, edges=edge_table, migrations=migration_table,
                sites=site_table, mutations=mutation_table)
        self.assertRaises(_msprime.LibraryError, f, [(0.1, -1)])
        length = ts.get_sequence_length()
        u = ts.get_num_nodes()
        for bad_node in [u, u + 1, 2 * u]:
            self.assertRaises(_msprime.LibraryError, f, [(0.1, bad_node)])
        for bad_pos in [-1, length, length + 1]:
            self.assertRaises(_msprime.LibraryError, f, [(length, 0)])

    def test_free(self):
        ts = self.get_tree_sequence()
        t = _msprime.SparseTree(
            ts, flags=_msprime.SAMPLE_COUNTS | _msprime.SAMPLE_LISTS)
        no_arg_methods = [
            t.get_left_root, t.get_index, t.get_left, t.get_right,
            t.get_num_sites, t.get_flags, t.get_sites, t.get_num_nodes]
        node_arg_methods = [
            t.get_parent, t.get_population, t.get_children, t.get_num_samples,
            t.get_num_tracked_samples]
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

    def verify_edge_table(self, edges, ts):
        """
        Verifies that the specified tree sequence and edges table contain
        the same data.
        """
        self.assertEqual(edges.num_rows, ts.get_num_edges())
        self.assertGreater(edges.num_rows, 0)
        left = edges.left
        right = edges.right
        parent = edges.parent
        child = edges.child
        for j in range(ts.get_num_edges()):
            t = ts.get_edge(j)
            self.assertEqual(t[0], left[j])
            self.assertEqual(t[1], right[j])
            self.assertEqual(t[2], parent[j])
            self.assertEqual(t[3], child[j])

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
        derived_state_offset = mutations.derived_state_offset
        for j in range(ts.get_num_mutations()):
            t = ts.get_mutation(j)
            self.assertEqual(t[0], site[j])
            self.assertEqual(t[1], node[j])
            length = derived_state_offset[j + 1] - derived_state_offset[j]
            self.assertEqual(len(t[2]), length)
            for c in t[2]:
                self.assertEqual(ord(c), derived_state[derived_state_offset[j]])

    def verify_provenance_table(self, provenances, ts):
        """
        Verifies that the specified tree sequence and provenances table contain
        the same data.
        """
        self.assertEqual(provenances.num_rows, ts.get_num_provenances())
        self.assertGreater(provenances.num_rows, 0)
        timestamp = provenances.timestamp
        timestamp_offset = provenances.timestamp_offset
        record = provenances.record
        record_offset = provenances.record_offset
        for j in range(ts.get_num_provenances()):
            t, r = ts.get_provenance(j)
            decoded = timestamp[timestamp_offset[j]: timestamp_offset[j + 1]]
            self.assertEqual("".join(chr(d) for d in decoded), t)
            decoded = record[record_offset[j]: record_offset[j + 1]]
            self.assertEqual("".join(chr(d) for d in decoded), r)

    def test_dump_tables(self):
        ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edges = _msprime.EdgeTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()
        provenances = _msprime.ProvenanceTable()
        kwargs = {
            "nodes": nodes,
            "edges": edges,
            "sites": sites,
            "mutations": mutations,
            "migrations": migrations,
            "provenances": provenances,
        }
        ts.dump_tables(**kwargs)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_migration_table(migrations, ts)
        self.verify_mutation_table(mutations, ts)
        self.verify_provenance_table(provenances, ts)

    def test_dump_tables_errors(self):
        ts = self.get_example_migration_tree_sequence()
        kwargs = {
            "nodes": _msprime.NodeTable(),
            "edges": _msprime.EdgeTable(),
            "sites": _msprime.SiteTable(),
            "mutations": _msprime.MutationTable(),
            "migrations": _msprime.MigrationTable(),
            "provenances": _msprime.ProvenanceTable()
        }

        for bad_type in [None, "", []]:
            for parameter in kwargs.keys():
                copy = dict(kwargs)
                copy[parameter] = bad_type
                self.assertRaises(TypeError, ts.dump_tables, **copy)
        self.assertRaises(TypeError, ts.dump_tables)
        self.assertRaises(TypeError, ts.dump_tables, nodes=_msprime.NodeTable())
        self.assertRaises(TypeError, ts.dump_tables, edges=_msprime.EdgeTable())
        # Both mutations and mutation types must be specified together
        self.assertRaises(
            TypeError, ts.dump_tables, nodes=_msprime.NodeTable(),
            edges=_msprime.EdgeTable(), mutations=_msprime.MutationTable())
        self.assertRaises(
            TypeError, ts.dump_tables, nodes=_msprime.NodeTable(),
            edges=_msprime.EdgeTable(),
            sites=_msprime.SiteTable())

    def test_dump_tables_optional_parameters(self):
        ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edges = _msprime.EdgeTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()
        provenances = _msprime.ProvenanceTable()

        ts.dump_tables(nodes=nodes, edges=edges)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)

        ts.dump_tables(nodes=nodes, edges=edges, migrations=migrations)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_migration_table(migrations, ts)

        ts.dump_tables(
            nodes=nodes, edges=edges, sites=sites,
            mutations=mutations)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_site_table(sites, ts)
        self.verify_mutation_table(mutations, ts)

        ts.dump_tables(nodes=nodes, edges=edges, provenances=provenances)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_provenance_table(provenances, ts)

    def test_load_tables(self):
        ex_ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edges = _msprime.EdgeTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()
        provenances = _msprime.ProvenanceTable()
        kwargs = {
            "nodes": nodes,
            "edges": edges,
            "sites": sites,
            "mutations": mutations,
            "migrations": migrations,
            "provenances": provenances,
        }
        ex_ts.dump_tables(**kwargs)
        ts = _msprime.TreeSequence()
        ts.load_tables(**kwargs)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_migration_table(migrations, ts)
        self.verify_site_table(sites, ts)
        self.verify_mutation_table(mutations, ts)
        self.verify_provenance_table(provenances, ts)

    def test_load_tables_errors(self):
        ex_ts = self.get_example_migration_tree_sequence()
        kwargs = {
            "nodes": _msprime.NodeTable(),
            "edges": _msprime.EdgeTable(),
            "sites": _msprime.SiteTable(),
            "mutations": _msprime.MutationTable(),
            "migrations": _msprime.MigrationTable(),
            "provenances": _msprime.ProvenanceTable(),
        }
        ex_ts.dump_tables(**kwargs)
        ts = _msprime.TreeSequence()

        for bad_type in [None, "", tuple()]:
            for parameter in kwargs.keys():
                copy = dict(kwargs)
                copy[parameter] = bad_type
                self.assertRaises(TypeError, ts.load_tables, **copy)
        self.assertRaises(TypeError, ts.load_tables)
        self.assertRaises(TypeError, ts.load_tables, nodes=_msprime.NodeTable())
        self.assertRaises(TypeError, ts.load_tables, edges=_msprime.EdgeTable())
        # Both mutations and mutation types must be specified together
        self.assertRaises(
            TypeError, ts.load_tables, nodes=_msprime.NodeTable(),
            edges=_msprime.EdgeTable(), mutations=_msprime.MutationTable())
        self.assertRaises(
            TypeError, ts.load_tables, nodes=_msprime.NodeTable(),
            edges=_msprime.EdgeTable(),
            sites=_msprime.SiteTable())

    def test_load_tables_optional_parameters(self):
        ex_ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edges = _msprime.EdgeTable()
        sites = _msprime.SiteTable()
        mutations = _msprime.MutationTable()
        migrations = _msprime.MigrationTable()
        provenances = _msprime.ProvenanceTable()
        ex_ts.dump_tables(
            nodes=nodes, edges=edges, sites=sites,
            mutations=mutations, migrations=migrations,
            provenances=provenances)

        ts = _msprime.TreeSequence()
        ts.load_tables(nodes=nodes, edges=edges)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.assertEqual(ts.get_num_migrations(), 0)
        self.assertEqual(ts.get_num_mutations(), 0)

        ts = _msprime.TreeSequence()
        ts.load_tables(nodes=nodes, edges=edges, migrations=migrations)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_migration_table(migrations, ts)
        self.assertEqual(ts.get_num_mutations(), 0)

        ts = _msprime.TreeSequence()
        ts.load_tables(
            nodes=nodes, edges=edges, sites=sites,
            mutations=mutations)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_site_table(sites, ts)
        self.verify_mutation_table(mutations, ts)
        self.assertEqual(ts.get_num_migrations(), 0)

        ts = _msprime.TreeSequence()
        ts.load_tables(nodes=nodes, edges=edges, provenances=provenances)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)
        self.verify_provenance_table(provenances, ts)
        self.assertEqual(ts.get_num_migrations(), 0)
        self.assertEqual(ts.get_num_mutations(), 0)

    def test_load_tables_sequence_length(self):
        ex_ts = self.get_example_migration_tree_sequence()
        nodes = _msprime.NodeTable()
        edges = _msprime.EdgeTable()
        ex_ts.dump_tables(nodes=nodes, edges=edges)
        L = ex_ts.get_sequence_length()

        ts = _msprime.TreeSequence()
        ts.load_tables(nodes=nodes, edges=edges)
        self.assertEqual(ts.get_sequence_length(), L)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)

        ts = _msprime.TreeSequence()
        ts.load_tables(nodes=nodes, edges=edges, sequence_length=L + 1)
        self.assertEqual(ts.get_sequence_length(), L + 1)
        self.verify_node_table(nodes, ts)
        self.verify_edge_table(edges, ts)

        # Zero is used as a special value to infer the value from the edges.
        for sequence_length in [0, L]:
            ts = _msprime.TreeSequence()
            ts.load_tables(
                nodes=nodes, edges=edges, sequence_length=sequence_length)
            self.assertEqual(ts.get_sequence_length(), L)
            self.verify_node_table(nodes, ts)
            self.verify_edge_table(edges, ts)

        ts = _msprime.TreeSequence()
        self.assertRaises(
            _msprime.LibraryError, ts.load_tables,
            nodes=nodes, edges=edges, sequence_length=L - 0.5)

        for bad_type in ["", None, []]:
            self.assertRaises(
                TypeError, ts.load_tables,
                nodes=nodes, edges=edges, sequence_length=bad_type)
        for bad_value in [-1, -100]:
            self.assertRaises(
                _msprime.LibraryError, ts.load_tables,
                nodes=nodes, edges=edges, sequence_length=bad_value)

    def test_node_table_add_row(self):
        table = _msprime.NodeTable()
        table.add_row()
        self.assertEqual(table.num_rows, 1)
        self.assertEqual(table.population, [-1])
        self.assertEqual(table.flags, [0])
        self.assertEqual(table.time, [0])
        self.assertEqual(list(table.metadata), [])
        self.assertEqual(list(table.metadata_offset), [0, 0])

        metadata = b"abcde"
        table = _msprime.NodeTable()
        table.add_row(flags=5, population=10, time=1.23, metadata=metadata)
        self.assertEqual(table.num_rows, 1)
        self.assertEqual(table.population, [10])
        self.assertEqual(table.flags, [5])
        self.assertEqual(table.time, [1.23])
        self.assertEqual(table.metadata.tobytes(), metadata)
        self.assertEqual(list(table.metadata_offset), [0, len(metadata)])

    def test_node_table_add_row_errors(self):
        table = _msprime.NodeTable()
        for bad_type in ["3423", None, []]:
            self.assertRaises(TypeError, table.add_row, flags=bad_type)
            self.assertRaises(TypeError, table.add_row, population=bad_type)
            self.assertRaises(TypeError, table.add_row, time=bad_type)
        for bad_type in [234, []]:
            self.assertRaises(TypeError, table.add_row, metadata=bad_type)

    def test_edge_table_add_row(self):
        table = _msprime.EdgeTable()
        table = _msprime.EdgeTable()
        table.add_row(left=0, right=1, parent=2, child=0)
        self.assertEqual(table.num_rows, 1)
        self.assertEqual(table.left, [0])
        self.assertEqual(table.right, [1])
        self.assertEqual(table.parent, [2])
        self.assertEqual(table.child, [0])

    def test_edge_table_add_row_errors(self):
        table = _msprime.EdgeTable()
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
            self.assertRaises(TypeError, table.add_row, metadata=bad_type)

    def test_add_row_data(self):
        nodes = _msprime.NodeTable()
        edges = _msprime.EdgeTable()
        mutations = _msprime.MutationTable()
        sites = _msprime.SiteTable()
        ts = self.get_example_migration_tree_sequence()
        ts.dump_tables(nodes=nodes, edges=edges, sites=sites, mutations=mutations)
        self.assertGreater(sites.num_rows, 0)
        self.assertGreater(mutations.num_rows, 0)

        new_nodes = _msprime.NodeTable()
        for j in range(ts.get_num_nodes()):
            flags, time, population, metadata = ts.get_node(j)
            new_nodes.add_row(
                flags=flags, time=time, population=population, metadata=metadata)
        self.assertEqual(list(nodes.time), list(new_nodes.time))
        self.assertEqual(list(nodes.flags), list(new_nodes.flags))
        self.assertEqual(list(nodes.population), list(new_nodes.population))
        self.assertEqual(list(nodes.metadata), list(new_nodes.metadata))

        new_edges = _msprime.EdgeTable()
        for j in range(ts.get_num_edges()):
            left, right, parent, children = ts.get_edge(j)
            new_edges.add_row(left, right, parent, children)
        self.assertEqual(list(edges.left), list(new_edges.left))
        self.assertEqual(list(edges.right), list(new_edges.right))
        self.assertEqual(list(edges.parent), list(new_edges.parent))
        self.assertEqual(list(edges.child), list(new_edges.child))

        new_sites = _msprime.SiteTable()
        for j in range(ts.get_num_sites()):
            t = ts.get_site(j)
            new_sites.add_row(t[0], t[1])
        self.assertEqual(list(new_sites.position), list(sites.position))
        self.assertEqual(list(new_sites.ancestral_state), list(sites.ancestral_state))

        new_mutations = _msprime.MutationTable()
        for j in range(ts.get_num_mutations()):
            site, node, derived_state, parent, metadata = ts.get_mutation(j)
            self.assertEqual(j, new_mutations.num_rows)
            self.assertEqual(metadata, b'')
            new_mutations.add_row(site, node, derived_state, parent)
        self.assertEqual(list(new_mutations.site), list(mutations.site))
        self.assertEqual(list(new_mutations.node), list(mutations.node))
        self.assertEqual(list(new_mutations.parent), list(mutations.parent))
        self.assertEqual(
            list(new_mutations.derived_state), list(mutations.derived_state))

    def test_site_table_add_row(self):
        table = _msprime.SiteTable()
        table.add_row(position=1, ancestral_state="0")
        self.assertEqual(table.num_rows, 1)
        self.assertEqual(table.position, [1])
        self.assertEqual(list(table.ancestral_state), [ord("0")])
        self.assertEqual(list(table.ancestral_state_offset), [0, 1])
        self.assertEqual(list(table.metadata), [])
        self.assertEqual(list(table.metadata_offset), [0, 0])

        metadata = b"abcde"
        table.add_row(position=2, ancestral_state="1", metadata=metadata)
        self.assertEqual(table.num_rows, 2)
        self.assertEqual(list(table.position), [1, 2])
        self.assertEqual(list(table.ancestral_state), [ord("0"), ord("1")])
        self.assertEqual(list(table.ancestral_state_offset), [0, 1, 2])
        self.assertEqual(table.metadata.tobytes(), metadata)
        self.assertEqual(list(table.metadata_offset), [0, 0, len(metadata)])

    def test_site_table_add_row_errors(self):
        table = _msprime.SiteTable()
        self.assertRaises(TypeError, table.add_row)
        for bad_type in [[], {}]:
            self.assertRaises(
                TypeError, table.add_row, position=bad_type, ancestral_state="1")
            self.assertRaises(
                TypeError, table.add_row, ancestral_state=bad_type, position=1)
            self.assertRaises(
                TypeError, table.add_row, ancestral_state="0", position=1,
                metadata=bad_type)


class TestSampleListIterator(LowLevelTestCase):
    """
    Tests for the low-level sample list iterator.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.SampleListIterator)
        self.assertRaises(TypeError, _msprime.SampleListIterator, None)
        ts = self.get_tree_sequence()
        flags = _msprime.SAMPLE_COUNTS | _msprime.SAMPLE_LISTS
        tree = _msprime.SparseTree(ts, flags=flags)
        for bad_type in [None, "1", []]:
            self.assertRaises(
                TypeError, _msprime.SampleListIterator, tree, bad_type)
        for bad_node in [-1, tree.get_num_nodes() + 1]:
            self.assertRaises(
                ValueError, _msprime.SampleListIterator, tree, bad_node)
        # Do nasty things...
        iterator = _msprime.SampleListIterator(tree, 1)
        del tree
        del ts
        self.assertEqual(list(iterator), [1])

    def test_iterator(self):
        ts = self.get_tree_sequence()
        flags = _msprime.SAMPLE_COUNTS | _msprime.SAMPLE_LISTS
        tree = _msprime.SparseTree(ts, flags=flags)
        for tree in _msprime.SparseTreeIterator(tree):
            self.verify_iterator(_msprime.SampleListIterator(tree, 1))
            self.verify_iterator(
                _msprime.SampleListIterator(tree, tree.get_left_root()))

    def test_sample_list(self):
        examples = [
            self.get_tree_sequence(num_samples=5, num_loci=10),
            self.get_tree_sequence(
                demographic_events=[
                    get_simple_bottleneck_event(0.2, proportion=1.0)])]
        flags = _msprime.SAMPLE_COUNTS | _msprime.SAMPLE_LISTS
        for ts in examples:
            st = _msprime.SparseTree(ts, flags=flags)
            for t in _msprime.SparseTreeIterator(st):
                # All sample nodes should have themselves.
                for j in range(ts.get_num_samples()):
                    samples = list(_msprime.SampleListIterator(t, j))
                    self.assertEqual(samples, [j])
                # All non-tree nodes should have 0
                for j in range(t.get_num_nodes()):
                    if t.get_parent(j) == NULL_NODE \
                            and t.get_left_child(j) == NULL_NODE:
                        samples = list(_msprime.SampleListIterator(t, j))
                        self.assertEqual(len(samples), 0)
                # The roots should have all samples.
                u = t.get_left_root()
                samples = []
                while u != NULL_NODE:
                    samples.extend(_msprime.SampleListIterator(t, u))
                    u = t.get_right_sib(u)
                self.assertEqual(sorted(samples), list(range(ts.get_num_samples())))


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
    Tests for the low-level sample list iterator.
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

    def test_uninitialised_tree_sequence(self):
        ts = _msprime.TreeSequence()
        self.assertRaises(ValueError, _msprime.LdCalculator, ts)

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
        rng = random.Random(6)
        start_positions = [rng.randint(0, m - 1) for _ in range(num_start_positions)]
        directions = [
            rng.choice([_msprime.FORWARD, _msprime.REVERSE])
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
