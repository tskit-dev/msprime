# Copyright (C) 2015-2018 University of Oxford
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
import collections
import heapq
import itertools
import math
import random
import unittest

import tskit

import tests
import _msprime
import _tskit

# Root node marker
NULL_NODE = -1


def uniform_recombination_map(num_loci=1, rate=0, L=None):
    """
    Returns a uniform recombination map for the specified number of loci
    and rate.
    """
    if L is None:
        L = num_loci
    return _msprime.RecombinationMap(
        num_loci=num_loci, positions=[0, L], rates=[rate, 0])


def get_simulation_model(name="hudson", reference_size=0.25, **kwargs):
    """
    Returns simulation model dictionary suitable for passing to the low-level API.
    """
    d = {"name": name, "reference_size": reference_size}
    d.update(kwargs)
    return d


def get_sweep_genic_selection_model(
        reference_size=0.25, position=0.5, start_frequency=0.1, end_frequency=0.9,
        alpha=0.1, dt=0.1):
    """
    Returns a sweep model for the specified parameters.
    """
    return get_simulation_model(
            name="sweep_genic_selection", reference_size=reference_size,
            position=position, start_frequency=start_frequency,
            end_frequency=end_frequency, alpha=alpha, dt=dt)


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
    tables = _msprime.LightweightTableCollection()
    samples = [(j % num_populations, 0) for j in range(num_samples)]
    migration_matrix = [1 for _ in range(num_populations**2)]
    for j in range(num_populations):
        migration_matrix[j * num_populations + j] = 0
    population_configuration = [
        get_population_configuration() for j in range(num_populations)]
    sim = _msprime.Simulator(
        samples,
        uniform_recombination_map(10, 0.1),
        _msprime.RandomGenerator(random_seed),
        tables=tables,
        population_configuration=population_configuration,
        migration_matrix=migration_matrix,
        store_migrations=store_migrations,
        model=get_simulation_model(reference_size=Ne))
    return sim, tables


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


class TestModule(tests.MsprimeTestCase):
    """
    Tests for module level stuff.
    """
    def test_library_versions(self):
        major, minor = _msprime.get_gsl_version()
        self.assertIsInstance(major, int)
        self.assertGreater(major, 0)
        self.assertIsInstance(minor, int)

    def test_gsl_error_handler(self):
        # Not a lot we can test here. Just make sure the error handler is unset when
        # we finish.
        _msprime.restore_gsl_error_handler()
        _msprime.unset_gsl_error_handler()
        _msprime.restore_gsl_error_handler()
        _msprime.unset_gsl_error_handler()


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
        self.assertGreaterEqual(sim.get_time(), 0.0)
        self.assertGreater(sim.get_num_ancestors(), 1)
        events = sim.get_num_common_ancestor_events()
        events += sim.get_num_recombination_events()
        events += sum(sim.get_num_migration_events())
        self.assertGreaterEqual(events, 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
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
        for j, (flags, t, pop, ind, metadata) in enumerate(nodes):
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
                pi_p[u] = sparse_tree.parent(u)
                u = pi_p[u]
        self.assertEqual(pi_p, pi)
        u = sparse_tree.left_root
        while u != NULL_NODE:
            self.assertEqual(pi_p[u], NULL_NODE)
            u = sparse_tree.right_sib(u)

    def verify_trees(self, sim, sorted_edges, ts):
        """
        Verifies that the specified set of (left, parent) sorted edges
        corresponds to correct trees for the specified simulation.
        """
        st_iter = ts.trees()
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

    def verify_completed_simulation(self, sim, tables):
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

        ts = tskit.TableCollection.fromdict(tables.asdict()).tree_sequence()
        left_sorted_edges = sorted(edges, key=lambda r: (r[0], r[2]))
        self.verify_trees(sim, left_sorted_edges, ts)
        ts_edges = list(ts.edges())
        j = 0
        for edge in edges:
            left, right, parent, child = edge
            self.assertEqual(left, ts_edges[j].left)
            self.assertEqual(right, ts_edges[j].right)
            self.assertEqual(parent, ts_edges[j].parent)
            self.assertEqual(child, ts_edges[j].child)
            j += 1
        assert j == len(ts_edges)

    def verify_random_parameters(self, seed):
        rng = random.Random(seed)
        n = rng.randint(2, 100)
        m = rng.randint(1, 10**6)
        rho = rng.uniform(0, 1000)
        N = rng.randint(1, 4)
        store_migrations = rng.choice([True, False])
        migration_matrix = [
            rng.random() * (j != k) for j in range(N) for k in range(N)]
        population_configuration = [
            get_population_configuration(rng.random(), rng.random())
            for j in range(N)]
        demographic_events = get_random_demographic_events(N, rng.randint(1, 5))
        start_time = rng.uniform(0, demographic_events[0]["time"])
        num_sampless = [0 for j in range(N)]
        num_sampless[0] = n
        random_seed = rng.randint(0, 2**31)
        segment_block_size = rng.randint(1, 100)
        node_mapping_block_size = rng.randint(1, 100)
        avl_node_block_size = rng.randint(1, 100)
        sim = _msprime.Simulator(
            samples=get_population_samples(*num_sampless),
            recombination_map=uniform_recombination_map(num_loci=m, rate=rho, L=m - 1),
            tables=_msprime.LightweightTableCollection(),
            random_generator=_msprime.RandomGenerator(random_seed),
            store_migrations=store_migrations,
            start_time=start_time,
            population_configuration=population_configuration,
            demographic_events=demographic_events,
            migration_matrix=migration_matrix,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size)
        for _ in range(3):
            # Check initial state
            self.assertEqual(0, sim.get_num_breakpoints())
            self.assertEqual(start_time, sim.get_time())
            self.assertEqual(n, sim.get_num_ancestors())
            self.assertEqual(0, sim.get_num_common_ancestor_events())
            self.assertEqual(0, sim.get_num_rejected_common_ancestor_events())
            self.assertEqual(0, sim.get_num_recombination_events())
            self.assertEqual(0, sum(sim.get_num_migration_events()))
            self.assertGreater(sim.get_num_avl_node_blocks(), 0)
            self.assertGreater(sim.get_num_segment_blocks(), 0)
            self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
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
                self.assertAlmostEqual(rho, sim.get_recombination_rate())
                self.assertEqual(segment_block_size, sim.get_segment_block_size())
                self.assertEqual(avl_node_block_size, sim.get_avl_node_block_size())
                self.assertEqual(
                    node_mapping_block_size, sim.get_node_mapping_block_size())
                # Run this for a tiny amount of time and check the state
                self.assertFalse(sim.run(start_time + 2e-8))
                self.verify_running_simulation(sim)
            sim.reset()

    def verify_tree_diffs(self, tree_sequence):
        n = tree_sequence.get_num_samples()
        L = tree_sequence.get_sequence_length()
        t = [tree_sequence.get_node(u)[1] for u in range(tree_sequence.get_num_nodes())]
        # Check some basic properties of the diffs.
        diffs = list(_tskit.TreeDiffIterator(tree_sequence))
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

    def verify_simulation(
            self, n, m, r, demographic_events=[], model=get_simulation_model()):
        """
        Runs the specified simulation and verifies its state.
        """
        # These tests don't work for n == 2
        assert n > 2
        random_seed = random.randint(0, 2**31)
        tables = _msprime.LightweightTableCollection()
        sim = _msprime.Simulator(
            samples=get_samples(n),
            recombination_map=uniform_recombination_map(num_loci=m, rate=r),
            random_generator=_msprime.RandomGenerator(random_seed),
            tables=tables, demographic_events=demographic_events,
            segment_block_size=1000, avl_node_block_size=1000,
            node_mapping_block_size=1000, model=model)
        for _ in range(3):
            # Run the sim for a tiny amount of time and check.
            if not sim.run():
                self.verify_running_simulation(sim)
                increment = 0.01
                t = sim.get_time() + increment
                while not sim.run(t):
                    self.assertLess(sim.get_time(), t)
                    self.verify_running_simulation(sim)
                    t += increment
            self.verify_completed_simulation(sim, tables)
            sim.reset()

    def test_random_sims(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            self.verify_random_parameters(j)

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
            samples=get_samples(n),
            recombination_map=uniform_recombination_map(num_loci=m, rate=1),
            random_generator=_msprime.RandomGenerator(1),
            tables=_msprime.LightweightTableCollection())
        # We run until time 0 to for initialisation
        sim.run(0)
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
            get_population_samples(n, 0, 0),
            uniform_recombination_map(),
            _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
            migration_matrix=migration_matrix,
            population_configuration=population_configuration,
            demographic_events=demographic_events)
        # Use a second instance to track the demographic events debugger.
        sim2 = _msprime.Simulator(
            get_samples(n), uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(), migration_matrix=migration_matrix,
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
            dt = 1e-10
            completed = sim.run(t + dt)
            for j in range(N):
                s = sim2.compute_population_size(j, t)
                self.assertGreater(s, 0)
            self.assertFalse(completed)
            self.assertEqual(sim.get_time(), t + dt)
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
        recomb_map = uniform_recombination_map()

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
                recomb_map, _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(), **kwargs)
        # samples recomb_map and random_seed are mandatory
        self.assertRaises(TypeError, _msprime.Simulator)
        self.assertRaises(TypeError, _msprime.Simulator, get_samples(10))
        self.assertRaises(
                TypeError, _msprime.Simulator, get_samples(10), random_generator=rng)
        self.assertRaises(
                TypeError, _msprime.Simulator, get_samples(10),
                recombination_map=recomb_map)
        self.assertRaises(
                TypeError, _msprime.Simulator, get_samples(10),
                recombination_map=recomb_map, random_generator=rng)
        # check types
        for bad_type in ["1", None, {}, int]:
            self.assertRaises(TypeError, f, tables=bad_type)
            self.assertRaises(TypeError, f, rng=bad_type)
            self.assertRaises(TypeError, f, samples=bad_type)
            self.assertRaises(TypeError, f, random_generator=bad_type)
            self.assertRaises(TypeError, f, recombination_rate=bad_type)
            self.assertRaises(TypeError, f, avl_node_block_size=bad_type)
            self.assertRaises(TypeError, f, segment_block_size=bad_type)
            self.assertRaises(TypeError, f, node_mapping_block_size=bad_type)
            self.assertRaises(TypeError, f, start_time=bad_type)
            self.assertRaises(TypeError, f, num_labels=bad_type)
        # Check for bad values.
        self.assertRaises(_msprime.InputError, f, avl_node_block_size=0)
        self.assertRaises(_msprime.InputError, f, segment_block_size=0)
        self.assertRaises(_msprime.InputError, f, node_mapping_block_size=0)
        self.assertRaises(_msprime.InputError, f, num_labels=0)
        self.assertRaises(_msprime.InputError, f, num_labels=-1)
        # Check for other type specific errors.
        self.assertRaises(OverflowError, f, avl_node_block_size=2**65)

    def test_num_labels(self):
        recomb_map = uniform_recombination_map()

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
                recomb_map, _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(), **kwargs)

        for num_labels in range(1, 10):
            sim = f(num_labels=num_labels)
            self.assertEqual(sim.get_num_labels(), num_labels)
            sim.run()
            self.assertEqual(sim.get_num_labels(), num_labels)

    def test_record_scaling(self):
        for Ne in [0.25, 1, 10, 1e6]:
            sim, ll_tables = get_example_simulator(
                10, Ne=Ne, num_populations=2, store_migrations=True)
            sim.run()
            tables = tskit.TableCollection.fromdict(ll_tables.asdict())
            sim_times = [node[1] for node in sim.get_nodes()]
            for j in range(len(tables.nodes)):
                node = tables.nodes[j]
                self.assertAlmostEqual(node.time, sim_times[j])
            self.assertGreater(sim.get_num_migrations(), 0)
            self.assertEqual(sim.get_num_migrations(), len(tables.migrations))
            sim_times = [r[-1] for r in sim.get_migrations()]
            for j in range(len(tables.migrations)):
                generation = tables.migrations[j].time
                self.assertAlmostEqual(generation, sim_times[j])

    def test_bad_run_args(self):
        sim, _ = get_example_simulator(10)
        for bad_type in ["x", []]:
            self.assertRaises(TypeError, sim.run, bad_type)
        for bad_value in [-1, -1e-6]:
            self.assertRaises(ValueError, sim.run, bad_value)

    def test_migrations(self):
        sim, ll_tables = get_example_simulator(
            10, num_populations=3, store_migrations=True)
        sim.run()
        tables = tskit.TableCollection.fromdict(ll_tables.asdict())
        self.assertGreater(sim.get_num_migrations(), 0)
        self.assertEqual(sim.get_num_migrations(), len(tables.migrations))
        sim_records = sim.get_migrations()
        self.assertEqual(len(sim_records), len(tables.migrations))
        for sim_r, ts_r in zip(sim_records, tables.migrations):
            self.assertEqual(len(sim_r), 6)
            # Left and right have been remapped, so might be slightly different.
            self.assertAlmostEqual(sim_r[0], ts_r.left)
            self.assertAlmostEqual(sim_r[1], ts_r.right)
            self.assertAlmostEqual(sim_r[2], ts_r.node)
            self.assertAlmostEqual(sim_r[3], ts_r.source)
            self.assertAlmostEqual(sim_r[4], ts_r.dest)
            self.assertAlmostEqual(sim_r[5], ts_r.time)

    def test_non_parametric_simulation_models(self):

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
                uniform_recombination_map(),
                _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(), **kwargs)
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
                uniform_recombination_map(),
                _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(),
                **kwargs)
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
                uniform_recombination_map(),
                _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(), **kwargs)
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

    def test_sweep_genic_selection_simulation_model_errors(self):
        L = 10

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
                uniform_recombination_map(L=L),
                _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(), **kwargs)

        # Partialy specified models.
        model = get_simulation_model("sweep_genic_selection")
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model("sweep_genic_selection", position=0.5)
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection", position=0.5, start_frequency=0.1)
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection", position=0.5, start_frequency=0.1,
            end_frequency=0.5)
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection", position=0.5, start_frequency=0.1,
            end_frequency=0.5, alpha=0.1)
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection", position=0.5, start_frequency=0.1,
            end_frequency=0.5, dt=0.1)
        self.assertRaises(ValueError, f, model=model)

        for bad_type in [None, str, "sdf"]:
            model = get_sweep_genic_selection_model(position=bad_type)
            self.assertRaises(TypeError, f, model=model)
            model = get_sweep_genic_selection_model(start_frequency=bad_type)
            self.assertRaises(TypeError, f, model=model)
            model = get_sweep_genic_selection_model(end_frequency=bad_type)
            self.assertRaises(TypeError, f, model=model)
            model = get_sweep_genic_selection_model(alpha=bad_type)
            self.assertRaises(TypeError, f, model=model)
            model = get_sweep_genic_selection_model(dt=bad_type)
            self.assertRaises(TypeError, f, model=model)

        for bad_position in [-1, L, L + 1]:
            model = get_sweep_genic_selection_model(position=bad_position)
            self.assertRaises(_msprime.InputError, f, model=model)

        # If we don't specify 2 labels we get an error.
        model = get_sweep_genic_selection_model()
        for bad_labels in [1, 3, 10]:
            sim = f(model=model, num_labels=bad_labels)
            self.assertRaises(_msprime.LibraryError, sim.run)

    def test_sweep_genic_selection_simulation_locus_round_trip(self):
        L = 10

        def f(num_samples=10, random_seed=1, **kwargs):
            return _msprime.Simulator(
                get_samples(num_samples),
                uniform_recombination_map(num_loci=L, rate=1, L=L),
                _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(), **kwargs)

        for j in range(L):
            model = get_sweep_genic_selection_model(position=j)
            sim = f(model=model)
            model = sim.get_model()
            # Because we have a flat map we should convert exactly.
            self.assertEqual(j, model["locus"])

    def test_store_migrations(self):
        def f(num_samples=10, random_seed=1, **kwargs):
            samples = [(j % 2, 0) for j in range(num_samples)]
            migration_matrix = [0, 1, 1, 0]
            population_configuration = [
                get_population_configuration() for j in range(2)]
            return _msprime.Simulator(
                samples, uniform_recombination_map(),
                _msprime.RandomGenerator(random_seed),
                _msprime.LightweightTableCollection(),
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
            return _msprime.Simulator(
                samples, uniform_recombination_map(), rng,
                _msprime.LightweightTableCollection())

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
            samples, uniform_recombination_map(), rng,
            _msprime.LightweightTableCollection(),
            population_configuration=[
                get_population_configuration() for _ in range(N)],
            migration_matrix=[0 for j in range(N * N)])
        self.assertEqual(samples, sim.get_samples())

    def test_deleting_tables(self):
        rng = _msprime.RandomGenerator(1)
        tables = _msprime.LightweightTableCollection()
        sim = _msprime.Simulator(
            get_samples(10), uniform_recombination_map(), rng, tables)
        del tables
        sim.run()

    def test_deleting_rng(self):
        rng = _msprime.RandomGenerator(1)
        sim = _msprime.Simulator(
            get_samples(10), uniform_recombination_map(), rng,
            _msprime.LightweightTableCollection())
        del rng
        sim.run()

    def test_deleting_recomb_map(self):
        rng = _msprime.RandomGenerator(1)
        recomb_map = uniform_recombination_map()
        sim = _msprime.Simulator(
            get_samples(10), recomb_map, rng, _msprime.LightweightTableCollection())
        del recomb_map
        sim.run()

    def test_defaults(self):
        n = 10
        sim = _msprime.Simulator(
            get_samples(n), uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection())
        self.assertEqual(sim.get_migration_matrix(), [0.0])
        self.assertEqual(
            sim.get_population_configuration(),
            [get_population_configuration(initial_size=0.25)])

    def test_bad_population_configurations(self):
        def f(population_configuration):
            return _msprime.Simulator(
                get_samples(2), uniform_recombination_map(), _msprime.RandomGenerator(1),
                _msprime.LightweightTableCollection(),
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
                get_samples(num_samples), uniform_recombination_map(), rng,
                _msprime.LightweightTableCollection(),
                population_configuration=population_configuration,
                migration_matrix=migration_matrix)
        for bad_type in [{}, None, 2, [""], [[]], [None]]:
            self.assertRaises(
                TypeError, _msprime.Simulator, get_samples(2),
                uniform_recombination_map(), rng, population_configuration=bad_type)
        # Cannot have empty list
        self.assertRaises(ValueError, f, 2, [])
        # Must provide population_configuration if a migration_matrix
        # is supplied.
        self.assertRaises(
            ValueError, _msprime.Simulator, get_samples(2),
            uniform_recombination_map(), rng, _msprime.LightweightTableCollection(),
            migration_matrix=[0, 0, 0, 0])

    def test_get_population_configurations(self):
        def f(num_samples, conf_tuples):
            population_configuration = [
                get_population_configuration(initial_size=p, growth_rate=a)
                for p, a in conf_tuples]
            N = len(population_configuration)
            migration_matrix = [0 for j in range(N) for k in range(N)]
            s = _msprime.Simulator(
                get_samples(num_samples),
                uniform_recombination_map(), _msprime.RandomGenerator(1),
                _msprime.LightweightTableCollection(),
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
                get_samples(2), uniform_recombination_map(), _msprime.RandomGenerator(1),
                _msprime.LightweightTableCollection(),
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
                ValueError, _msprime.Simulator,
                get_samples(2), uniform_recombination_map(), _msprime.RandomGenerator(1),
                _msprime.LightweightTableCollection(), population_configuration=pop_conf)

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
                    get_samples(2), uniform_recombination_map(),
                    _msprime.RandomGenerator(1),
                    _msprime.LightweightTableCollection(),
                    migration_matrix=migration_matrix,
                    population_configuration=population_configuration)
                self.assertEqual(migration_matrix, sim.get_migration_matrix())

    def test_bad_demographic_event_types(self):
        def f(events):
            _msprime.Simulator(
                get_samples(2), uniform_recombination_map(), _msprime.RandomGenerator(1),
                _msprime.LightweightTableCollection(), demographic_events=events)
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
                get_samples(2), uniform_recombination_map(), _msprime.RandomGenerator(1),
                _msprime.LightweightTableCollection(), demographic_events=events,
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
            uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
            demographic_events=events)

    def test_seed_equality(self):
        simulations = [
            {
                "samples": get_samples(10),
            }, {
                "samples": get_samples(100),
                "demographic_events": [
                    get_simple_bottleneck_event(0.01, 0, 1.0)],
            }, {
                "samples": get_samples(10),
                "recombination_map": uniform_recombination_map(num_loci=10, rate=1)
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
            params["recombination_map"] = uniform_recombination_map()
            params["tables"] = _msprime.LightweightTableCollection()
            sim1 = _msprime.Simulator(**params)
            params["random_generator"] = _msprime.RandomGenerator(seed)
            params["recombination_map"] = uniform_recombination_map()
            params["tables"] = _msprime.LightweightTableCollection()
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
            get_population_samples(5, 5), uniform_recombination_map(),
            _msprime.RandomGenerator(1), _msprime.LightweightTableCollection(),
            population_configuration=population_configuration,
            migration_matrix=[0, 0, 0, 0])
        self.assertRaises(_msprime.LibraryError, sim.run)

    def test_simple_event_counters(self):
        for n in [2, 10, 20]:
            sim = _msprime.Simulator(
                get_samples(n), uniform_recombination_map(), _msprime.RandomGenerator(1),
                _msprime.LightweightTableCollection())
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
            get_samples(n), uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
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
            get_samples(n), uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
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
            uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
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
            get_samples(n), uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
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
        sim.run(t + 2 * dt)
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
        dt = 1e-9
        sim = _msprime.Simulator(
            get_population_samples(n, n),
            uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
            population_configuration=[
                get_population_configuration(),
                get_population_configuration()],
            demographic_events=[
                get_simple_bottleneck_event(t1, population=0, proportion=1),
                get_simple_bottleneck_event(t2, population=1, proportion=1),
                get_mass_migration_event(t3, source=0, dest=1, proportion=1)],
            migration_matrix=[0, 0, 0, 0])
        sim.run(t1 + dt)
        pop_sizes = [0, 0]
        for ind in sim.get_ancestors():
            for _, _, _, pop_id in ind:
                pop_sizes[pop_id] += 1
        self.assertEqual(pop_sizes[0], 1)
        sim.run(t2 + dt)
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
            get_samples(n), uniform_recombination_map(num_loci=10, rate=10),
            _msprime.RandomGenerator(1), _msprime.LightweightTableCollection(),
            population_configuration=population_configuration,
            migration_matrix=[0.0, 0.0, 0.0, 0.0])
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
            get_population_samples(5, 5, 0),
            uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(),
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
        sim = _msprime.Simulator(
            get_samples(10), uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection())
        times = set()
        for _ in range(10):
            sim.run()
            t = sim.get_time()
            self.assertNotIn(t, times)
            times.add(t)
            sim.reset()
            self.assertEqual(sim.get_time(), 0)


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
            total_rate / (num_loci - 1), rm.get_per_locus_recombination_rate())
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

    def test_range_errors(self):
        rm = _msprime.RecombinationMap(100, [0, 10], [0.1, 0])
        self.assertRaises(ValueError, rm.physical_to_discrete_genetic, -1)
        self.assertRaises(ValueError, rm.physical_to_discrete_genetic, 10.1)
        self.assertRaises(ValueError, rm.physical_to_genetic, -1)
        self.assertRaises(ValueError, rm.physical_to_genetic, 10.1)
        self.assertRaises(ValueError, rm.genetic_to_physical, -1)
        self.assertRaises(ValueError, rm.genetic_to_physical, 100.1)


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
    def test_basic_constructor(self):
        self.assertRaises(TypeError, _msprime.MutationGenerator)
        mg = _msprime.MutationGenerator(_msprime.RandomGenerator(1), 0)
        self.assertEqual(mg.get_mutation_rate(), 0)

    def test_rng(self):
        for bad_type in ["x", {}, None]:
            self.assertRaises(
                TypeError, _msprime.MutationGenerator, random_generator=bad_type)

    def test_mutation_rate(self):
        rng = _msprime.RandomGenerator(1)
        for bad_type in ["x", {}, None]:
            self.assertRaises(TypeError, _msprime.MutationGenerator, rng, bad_type)
        for bad_value in [-1, -1e-3]:
            self.assertRaises(ValueError, _msprime.MutationGenerator, rng, bad_value)
        for rate in [0, 1e-12, 1e12, 1000, 0.01]:
            mutgen = _msprime.MutationGenerator(rng, rate)
            self.assertEqual(mutgen.get_mutation_rate(), rate)

    def test_time_interval(self):
        rng = _msprime.RandomGenerator(1)
        for bad_type in ["x", {}, None]:
            with self.assertRaises(TypeError):
                _msprime.MutationGenerator(rng, 0, start_time=bad_type)
            with self.assertRaises(TypeError):
                _msprime.MutationGenerator(rng, 0, end_time=bad_type)
        for start_time, end_time in [(1, 0), (-1, -2), (200, 100)]:
            with self.assertRaises(_msprime.LibraryError):
                _msprime.MutationGenerator(
                    rng, 0, start_time=start_time, end_time=end_time)

    def test_alphabet(self):
        rng = _msprime.RandomGenerator(1)
        for bad_type in ["x", {}, None]:
            self.assertRaises(
                TypeError, _msprime.MutationGenerator, random_generator=rng,
                mutation_rate=0, alphabet=bad_type)
        for bad_value in [-1, 2, 10**6]:
            self.assertRaises(
                ValueError, _msprime.MutationGenerator, random_generator=rng,
                mutation_rate=0, alphabet=bad_value)
        for alphabet in [0, 1]:
            mg = _msprime.MutationGenerator(
                random_generator=rng, mutation_rate=0, alphabet=alphabet)
            self.assertEqual(alphabet, mg.get_alphabet())


class TestDemographyDebugger(unittest.TestCase):
    """
    Tests for the demography debugging interface.
    """
    def get_simulator(self, events):
        return _msprime.Simulator(
            get_samples(2), uniform_recombination_map(), _msprime.RandomGenerator(1),
            _msprime.LightweightTableCollection(), demographic_events=events)

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
