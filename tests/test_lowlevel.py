#
# Copyright (C) 2015-2020 University of Oxford
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
import functools
import inspect
import io
import itertools
import math
import pathlib
import platform
import random
import tempfile
import unittest

import _tskit
import numpy as np
import tskit

import tests
from msprime import _msprime

# Root node marker
NULL_NODE = -1

IS_WINDOWS = platform.system() == "Windows"


def uniform_rate_map(L=1, rate=0):
    """
    Returns a uniform recombination map for the specified sequence length
    and rate.
    """
    return {"position": [0, L], "rate": [rate]}


def get_mutation_model(model=0):
    """
    Returns a simple mutation model instance.
    """
    if model == 0:
        alleles = ["0", "1"]
        root_distribution = [1, 0]
        transition_matrix = [[0, 1], [1, 0]]
    elif model == 1:
        alleles = ["A", "C", "T", "G"]
        root_distribution = [0.25, 0.25, 0.25, 0.25]
        transition_matrix = np.zeros((4, 4))
        transition_matrix[:, 0] = 1
    return _msprime.MatrixMutationModel(
        alleles=alleles,
        root_distribution=root_distribution,
        transition_matrix=transition_matrix,
    )


def get_simulation_model(name="hudson", **kwargs):
    """
    Returns simulation model dictionary suitable for passing to the low-level API.
    """
    d = {"name": name}
    d.update(kwargs)
    return d


def get_sweep_genic_selection_model(
    reference_size=0.25,
    position=0.5,
    start_frequency=0.1,
    end_frequency=0.9,
    alpha=0.1,
    dt=0.1,
):
    """
    Returns a sweep model for the specified parameters.
    """
    return get_simulation_model(
        name="sweep_genic_selection",
        position=position,
        start_frequency=start_frequency,
        end_frequency=end_frequency,
        alpha=alpha,
        dt=dt,
    )


def get_population_configuration(growth_rate=0.0, initial_size=1.0):
    """
    Returns a population configuration dictionary suitable for passing
    to the low-level API.
    """
    return {"growth_rate": growth_rate, "initial_size": initial_size}


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
    time=0.0, population=-1, initial_size=None, growth_rate=None
):
    """
    Returns a population change event for the specified values.
    """
    ret = {
        "type": "population_parameters_change",
        "time": time,
        "population": population,
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
        time, population, initial_size=size, growth_rate=0
    )


def get_growth_rate_change_event(time=0.0, growth_rate=1.0, population=-1):
    """
    Returns a growth_rate change demographic event.
    """
    return get_population_parameters_change_event(
        time, population, growth_rate=growth_rate
    )


def get_migration_rate_change_event(time=0.0, migration_rate=1.0, source=-1, dest=-1):
    """
    Returns a migration_rate change demographic event.
    """
    return {
        "type": "migration_rate_change",
        "migration_rate": migration_rate,
        "time": time,
        "source": source,
        "dest": dest,
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
        "proportion": proportion,
    }


def get_simple_bottleneck_event(time=0.0, population=0, proportion=1):
    """
    Returns a simple bottleneck demographic event.
    """
    return {
        "type": "simple_bottleneck",
        "time": time,
        "population": population,
        "proportion": proportion,
    }


def get_instantaneous_bottleneck_event(time=0.0, population=0, strength=10):
    """
    Returns a instantaneous bottleneck demographic event.
    """
    return {
        "type": "instantaneous_bottleneck",
        "time": time,
        "population": population,
        "strength": strength,
    }


def get_migration_matrix(num_populations, value=1.0):
    """
    Returns a simple migration matrix.
    """
    m = np.full((num_populations, num_populations), value)
    np.fill_diagonal(m, 0)
    return m


def make_minimal_tables():
    """
    Returns the simplest set of tables we can run a simulation from
    """
    tables = tskit.TableCollection(1)
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
    tables.populations.add_row()
    ll_tables = _msprime.LightweightTableCollection(1)
    ll_tables.fromdict(tables.asdict())
    return ll_tables


def make_sim(
    samples=2, *, sequence_length=1, num_populations=1, random_seed=42, **kwargs
):
    """
    Helper function to make a functioning simulator object.
    """
    tables = tskit.TableCollection(sequence_length)
    if isinstance(samples, int):
        for _ in range(samples):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=0)
    else:
        for population, time in samples:
            tables.nodes.add_row(
                flags=tskit.NODE_IS_SAMPLE, time=time, population=population
            )
    for j in range(num_populations):
        tables.populations.add_row(f"pop_{j}".encode())
    ll_tables = _msprime.LightweightTableCollection(tables.sequence_length)
    ll_tables.fromdict(tables.asdict())
    sim = _msprime.Simulator(
        ll_tables, random_generator=_msprime.RandomGenerator(random_seed), **kwargs
    )
    return sim


def get_example_simulator(
    num_samples=10, Ne=0.25, random_seed=1, num_populations=1, store_migrations=False
):
    samples = [(j % num_populations, 0) for j in range(num_samples)]
    migration_matrix = get_migration_matrix(num_populations, 1)
    population_configuration = [
        get_population_configuration() for j in range(num_populations)
    ]
    sim = make_sim(
        samples=samples,
        sequence_length=10,
        num_populations=num_populations,
        recombination_map=uniform_rate_map(10, 0.1),
        population_configuration=population_configuration,
        migration_matrix=migration_matrix,
        store_migrations=store_migrations,
        model=get_simulation_model(reference_size=Ne),
    )
    return sim


def get_random_demographic_events(num_populations, num_events, rng=None):
    """
    Return some random demographic events for the specified number
    of populations. Note: we return num_events of *each type*.
    """
    if rng is None:
        rng = random.Random(1234)
    events = []
    for _ in range(num_events):
        events.append(
            get_size_change_event(
                time=rng.random(),
                size=rng.random(),
                population=rng.randint(-1, num_populations - 1),
            )
        )
        events.append(
            get_growth_rate_change_event(
                time=rng.random(),
                growth_rate=rng.random(),
                population=rng.randint(-1, num_populations - 1),
            )
        )
        if num_populations > 1:
            source = -1
            dest = -1
            if rng.random() < 0.5:
                # Don't pick diagonal elements
                while source == dest:
                    source = rng.randint(0, num_populations - 1)
                    dest = rng.randint(0, num_populations - 1)
            events.append(
                get_migration_rate_change_event(
                    time=rng.random(),
                    migration_rate=rng.random(),
                    source=source,
                    dest=dest,
                )
            )
            # Add a mass migration event.
            source = rng.randint(0, num_populations - 1)
            dest = source
            while dest == source:
                dest = rng.randint(0, num_populations - 1)
            # We use proportion of 0 or 1 so that we can test deterministically
            events.append(
                get_mass_migration_event(
                    time=rng.random(),
                    proportion=rng.choice([0, 1]),
                    source=source,
                    dest=dest,
                )
            )
            # Add some bottlenecks
            events.append(
                get_simple_bottleneck_event(
                    time=rng.random(),
                    proportion=rng.uniform(0, 0.25),
                    population=rng.randint(0, num_populations - 1),
                )
            )
            events.append(
                get_instantaneous_bottleneck_event(
                    time=rng.random(),
                    strength=rng.uniform(0, 0.01),
                    population=rng.randint(0, num_populations - 1),
                )
            )

    sorted_events = sorted(events, key=lambda x: x["time"])
    return sorted_events


class TestModule(unittest.TestCase):
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
    for _ in range(n):
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


class LowLevelTestCase(tests.SequenceEqualityMixin, unittest.TestCase):
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
        self.assertGreaterEqual(sim.num_breakpoints, 0)
        self.assertGreaterEqual(sim.time, 0.0)
        self.assertGreater(sim.num_ancestors, 1)
        events = sim.num_common_ancestor_events
        events += sim.num_recombination_events
        events += np.sum(sim.num_migration_events)
        self.assertGreaterEqual(events, 0)
        self.assertGreater(sim.num_avl_node_blocks, 0)
        self.assertGreater(sim.num_segment_blocks, 0)
        self.assertGreaterEqual(sim.num_fenwick_rebuilds, 0)
        self.assertGreater(sim.num_node_mapping_blocks, 0)
        n = sum(sim.tables.asdict()["nodes"]["flags"] == tskit.NODE_IS_SAMPLE)
        L = sim.sequence_length
        N = sim.num_populations
        ancestors = sim.ancestors
        self.assertEqual(len(ancestors), sim.num_ancestors)
        for ind in ancestors:
            the_pop_id = ind[0][-1]
            for left, right, node, pop_id in ind:
                self.assertTrue(0 <= left < L)
                self.assertTrue(1 <= right <= L)
                self.assertGreaterEqual(node, 0)
                self.assertTrue(0 <= pop_id < N)
                self.assertEqual(the_pop_id, pop_id)
        breakpoints = [0] + list(sim.breakpoints) + [L]
        self.assertEqual(len(breakpoints), sim.num_breakpoints + 2)
        self.assertEqual(breakpoints, sorted(breakpoints))
        tables = tskit.TableCollection.fromdict(sim.tables.asdict())
        nodes = tables.nodes
        self.assertEqual(len(nodes), sim.num_nodes)
        for j, node in enumerate(nodes):
            if j < n:
                self.assertEqual(node.time, 0.0)
                self.assertEqual(node.flags, 1)
            else:
                self.assertGreater(node.time, 0.0)
                self.assertEqual(node.flags, 0)
            self.assertTrue(0 <= node.population < sim.num_populations)
            self.assertEqual(len(node.metadata), 0)

        edges = tables.edges
        self.assertEqual(len(edges), sim.num_edges)
        # The amount of ancestral material in the edges and
        # the extant segments over all intervals should be either n (if
        # the given interval has not fully coalesced yet) or n - 1 (if
        # full coalescence has occured).
        segments_am = [0 for b in breakpoints[:-1]]
        for ind in ancestors:
            for left, right, _, _ in ind:
                j = breakpoints.index(left)
                while breakpoints[j] < right:
                    segments_am[j] += 1
                    j += 1
        records_am = [0 for b in breakpoints[:-1]]
        for edge in edges:
            j = breakpoints.index(edge.left)
            while breakpoints[j] < edge.right:
                # We are assuming a binary coalescent here. We would need to
                # do something more complicated for more general models.
                records_am[j] += 0.5
                j += 1
        for segment_am, record_am in zip(segments_am, records_am):
            if segment_am == 0:
                self.assertEqual(record_am, n - 1)
            else:
                self.assertEqual(segment_am + record_am, n)
        self.assertEqual(len(tables.migrations), sim.num_migrations)
        if sim.store_migrations:
            self.assertGreaterEqual(
                sim.num_migrations, np.sum(sim.num_migration_events)
            )
        sim.verify(True)

    def verify_squashed_edges(self, sorted_edges):
        """
        Checks to see if there were any unsquashed edges in the specified
        set of time sorted edges.
        """
        last_edge = sorted_edges[0]
        for edge in sorted_edges[1:]:
            if edge.parent == last_edge.parent and edge.child == last_edge.child:
                self.assertNotEqual(last_edge.right, edge.left)
            last_edge = edge

    def verify_completed_simulation(self, sim):
        """
        Verifies the state of the specified completed simulation.
        """
        self.assertEqual(sim.ancestors, [])
        self.assertEqual(sim.num_ancestors, 0)
        self.assertGreaterEqual(sim.num_breakpoints, 0)
        self.assertGreater(sim.num_edges, 0)
        self.assertGreater(sim.time, 0.0)
        events = sim.num_common_ancestor_events
        events += sim.num_recombination_events
        events += sim.num_rejected_common_ancestor_events
        events += sum(sim.num_migration_events)
        self.assertGreater(events, 0)
        model_name = sim.model["name"]
        if model_name == "hudson":
            self.assertEqual(sim.num_rejected_common_ancestor_events, 0)
        elif model_name in ["smc", "smc_prime"]:
            self.assertGreaterEqual(sim.num_rejected_common_ancestor_events, 0)
        self.assertGreater(sim.num_avl_node_blocks, 0)
        self.assertGreater(sim.num_segment_blocks, 0)
        self.assertGreater(sim.num_node_mapping_blocks, 0)

        tables = tskit.TableCollection.fromdict(sim.tables.asdict())
        edges = list(tables.edges)
        self.assertGreater(len(edges), 0)
        self.assertEqual(len(edges), sim.num_edges)
        # Edges should be returned in canonical order
        self.assertEqual(
            edges, sorted(edges, key=lambda e: (e.parent, e.child, e.left))
        )
        # Nodes should be in nondecreasing time order
        times = [node.time for node in tables.nodes]
        self.assertEqual(times, sorted(times))
        self.assertEqual(times[-1], sim.time)
        self.verify_squashed_edges(edges)

        ts = tables.tree_sequence()
        ts_edges = list(ts.edges())
        j = 0
        for edge in edges:
            self.assertEqual(edge.left, ts_edges[j].left)
            self.assertEqual(edge.right, ts_edges[j].right)
            self.assertEqual(edge.parent, ts_edges[j].parent)
            self.assertEqual(edge.child, ts_edges[j].child)
            j += 1
        assert j == len(ts_edges)

    def verify_random_parameters(self, seed):
        rng = random.Random(seed)
        n = rng.randint(2, 100)
        m = rng.randint(1, 10 ** 6)
        rho = rng.uniform(0, 1000)
        N = rng.randint(1, 4)
        store_migrations = rng.choice([True, False])
        migration_matrix = [
            [rng.random() * (j != k) for j in range(N)] for k in range(N)
        ]
        population_configuration = [
            get_population_configuration(rng.random(), rng.random()) for j in range(N)
        ]
        demographic_events = get_random_demographic_events(N, rng.randint(1, 5), rng)
        start_time = rng.uniform(0, demographic_events[0]["time"])
        num_sampless = [0 for j in range(N)]
        num_sampless[0] = n
        random_seed = rng.randint(0, 2 ** 31)
        segment_block_size = rng.randint(1, 100)
        node_mapping_block_size = rng.randint(1, 100)
        avl_node_block_size = rng.randint(1, 100)
        sim = make_sim(
            samples=get_population_samples(*num_sampless),
            sequence_length=m,
            num_populations=N,
            random_seed=random_seed,
            recombination_map=uniform_rate_map(L=m, rate=rho),
            store_migrations=store_migrations,
            start_time=start_time,
            population_configuration=population_configuration,
            demographic_events=demographic_events,
            migration_matrix=migration_matrix,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size,
        )
        for _ in range(3):
            # Check initial state
            self.assertEqual(0, sim.num_breakpoints)
            self.assertEqual(start_time, sim.time)
            self.assertEqual(n, sim.num_ancestors)
            self.assertEqual(0, sim.num_common_ancestor_events)
            self.assertEqual(0, sim.num_rejected_common_ancestor_events)
            self.assertEqual(0, sim.num_recombination_events)
            self.assertEqual(0, np.sum(sim.num_migration_events))
            self.assertGreater(sim.num_avl_node_blocks, 0)
            self.assertGreater(sim.num_segment_blocks, 0)
            self.assertGreater(sim.num_node_mapping_blocks, 0)
            self.assertEqual(sim.sequence_length, m)
            self.assertEqual(n, len(sim.ancestors))
            self.assertTrue(np.array_equal(sim.migration_matrix, migration_matrix))
            self.assertEqual(sim.population_configuration, population_configuration)
            a = 0
            nodes = set()
            pop_sizes = [0 for _ in population_configuration]
            for ind in sim.ancestors:
                self.assertEqual(len(ind), 1)
                left, right, node, pop_id = ind[0]
                pop_sizes[pop_id] += 1
                self.assertEqual(left, 0)
                self.assertEqual(right, m)
                self.assertFalse(node in nodes)
                nodes.add(node)
                a += 1
            for n1, n2 in zip(num_sampless, pop_sizes):
                self.assertEqual(n1, n2)
            self.assertEqual(a, n)
            for _ in range(3):
                # Check the getters to ensure we've got the right values.
                self.assertEqual(m, sim.sequence_length)
                self.assertEqual(segment_block_size, sim.segment_block_size)
                self.assertEqual(avl_node_block_size, sim.avl_node_block_size)
                self.assertEqual(node_mapping_block_size, sim.node_mapping_block_size)
                # Run this for a tiny amount of time and check the state
                ret = sim.run(start_time + 2e-8)
                self.assertEqual(ret, _msprime.EXIT_MAX_TIME)
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
            for l_in, _r_in, _p, _c in edges_in:
                self.assertEqual(l_in, left)
            for _l_out, r_out, _p, _c in edges_out:
                self.assertEqual(r_out, left)
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
                    edge, (py_edge.left, py_edge.right, py_edge.parent, py_edge.child)
                )
            for edge, py_edge in zip(diff[2], py_diff[2]):
                self.assertEqual(
                    edge, (py_edge.left, py_edge.right, py_edge.parent, py_edge.child)
                )

    def verify_simulation(self, n, m, r, demographic_events=(), model=None):
        """
        Runs the specified simulation and verifies its state.
        """
        if model is None:
            model = get_simulation_model()
        # These tests don't work for n == 2
        assert n > 2
        random_seed = random.randint(0, 2 ** 31)
        sim = make_sim(
            samples=n,
            sequence_length=m,
            random_seed=random_seed,
            recombination_map=uniform_rate_map(L=m, rate=r),
            demographic_events=list(demographic_events),
            segment_block_size=1000,
            avl_node_block_size=1000,
            node_mapping_block_size=1000,
            model=model,
        )
        for _ in range(3):
            # Run the sim for a tiny amount of time and check.
            if sim.run(1e-9) != _msprime.EXIT_COALESCENCE:
                self.verify_running_simulation(sim)
                increment = 0.01
                t = sim.time + increment
                while sim.run(t) != _msprime.EXIT_COALESCENCE:
                    self.assertLessEqual(sim.time, t)
                    self.verify_running_simulation(sim)
                    t += increment
            self.verify_completed_simulation(sim)
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
            5,
            10,
            10.0,
            demographic_events=[get_simple_bottleneck_event(time=0.2, proportion=1)],
        )
        self.verify_simulation(3, 10, 1.0, model=get_simulation_model("smc"))
        self.verify_simulation(4, 10, 2.0, model=get_simulation_model("smc_prime"))

    def test_event_by_event(self):
        n = 10
        m = 100
        sim = make_sim(
            n, sequence_length=m, recombination_map=uniform_rate_map(L=m, rate=1),
        )
        # We run until time 0 to for initialisation
        sim.run(0)
        self.assertEqual(sim.time, 0)
        ancestors = sim.ancestors
        self.assertEqual(len(ancestors), n)
        nodes = []
        for ancestor in ancestors:
            self.assertEqual(len(ancestor), 1)
            for left, right, node, pop_id in ancestor:
                self.assertEqual(left, 0)
                self.assertEqual(right, m)
                self.assertEqual(pop_id, 0)
                nodes.append(node)
        self.assertEqual(sorted(nodes), list(range(n)))
        events = 0
        while sim.run(max_events=1) == _msprime.EXIT_MAX_EVENTS:
            events += 1
            total_events = (
                sim.num_common_ancestor_events
                + sim.num_recombination_events
                + np.sum(sim.num_migration_events)
            )
            self.assertEqual(events, total_events)
        self.verify_completed_simulation(sim)

    def test_demographic_events(self):
        rng = random.Random(11)
        n = 10
        N = 3
        migration_matrix = np.array(
            [[rng.random() * (j != k) for j in range(N)] for k in range(N)]
        )
        population_configuration = [
            get_population_configuration(rng.random(), rng.random()) for j in range(N)
        ]
        demographic_events = get_random_demographic_events(N, 10, rng)
        start_times = [0 for j in range(N)]
        # Rescale time back to very small values so we know that they
        # will definitely happen
        for event in demographic_events:
            event["time"] *= 1e-6
        sim = make_sim(
            get_population_samples(n, 0, 0),
            num_populations=3,
            migration_matrix=migration_matrix,
            population_configuration=population_configuration,
            demographic_events=demographic_events,
        )
        # Use a second instance to track the demographic events debugger.
        sim2 = make_sim(
            get_population_samples(n, 0, 0),
            num_populations=3,
            migration_matrix=migration_matrix,
            population_configuration=population_configuration,
            demographic_events=demographic_events,
        )
        self.assertTrue(np.array_equal(sim.migration_matrix, migration_matrix))
        self.assertEqual(sim.population_configuration, population_configuration)

        # Now run the demography debug forward.
        next_event_time = sim2.debug_demography()
        self.assertTrue(np.array_equal(sim2.migration_matrix, migration_matrix))
        self.assertEqual(sim2.population_configuration, population_configuration)
        self.assertRaises(_msprime.LibraryError, sim.compute_population_size, N + 1, 0)
        # For each event we now run the simulator forward until this time
        # and make sure that the internal state is what it should be.
        for event in demographic_events:
            t = event["time"]
            event_type = event["type"]
            self.assertEqual(next_event_time, t)
            self.assertTrue(np.array_equal(sim2.migration_matrix, migration_matrix))
            dt = 1e-10
            status = sim.run(t + dt)
            for j in range(N):
                s = sim2.compute_population_size(j, t)
                self.assertGreater(s, 0)
            self.assertEqual(status, _msprime.EXIT_MAX_TIME)
            self.assertEqual(sim.time, t + dt)
            if event_type == "migration_rate_change":
                source = event["source"]
                dest = event["dest"]
                rate = event["migration_rate"]
                if source == -1 and dest == -1:
                    migration_matrix[:] = rate
                    np.fill_diagonal(migration_matrix, 0)
                else:
                    assert source >= 0 and dest >= 0
                    migration_matrix[source, dest] = rate
            elif event_type == "mass_migration":
                source = event["source"]
                proportion = event["proportion"]
                pop_sizes = [0 for j in range(N)]
                for ind in sim.ancestors:
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
            self.assertTrue(np.array_equal(sim.migration_matrix, migration_matrix))
            self.assertEqual(sim.population_configuration, population_configuration)
            next_event_time = sim2.debug_demography()
        self.assertTrue(math.isinf(next_event_time))


class TestSimulator(LowLevelTestCase):
    """
    Tests for the low-level interface to the simulator.
    """

    def test_run(self):
        sim = make_sim(10)
        status = sim.run()
        self.assertEqual(status, _msprime.EXIT_COALESCENCE)

        sim.reset()
        for bad_type in ["sdf", [], {}]:
            with self.assertRaises(TypeError):
                sim.run(end_time=bad_type)
            with self.assertRaises(TypeError):
                sim.run(max_events=bad_type)

        with self.assertRaises(ValueError):
            sim.run(end_time=-1)
        with self.assertRaises(ValueError):
            sim.run(max_events=0)
        self.assertEqual(sim.time, 0)
        sim.run(max_events=1)
        self.assertGreater(sim.time, 0)

    def test_print_state(self):
        sim = make_sim(10)
        with tempfile.TemporaryFile("w+") as f:
            sim.print_state(f)
            f.seek(0)
            output = f.read()
        self.assertGreater(len(output), 0)
        self.assertTrue(output.startswith("simulation model"))

    def test_uninitialised(self):
        # Create a properly initialised instance so we can run inspect on it.
        sim = make_sim(10)
        attributes = []
        methods = []
        for name, value in inspect.getmembers(sim):
            if not name.startswith("__"):
                if inspect.isbuiltin(value):
                    methods.append(name)
                else:
                    attributes.append(name)
        uninitialised_sim = _msprime.Simulator.__new__(_msprime.Simulator)
        for attr in attributes:
            with self.assertRaises(SystemError):
                getattr(uninitialised_sim, attr)
        for method_name in methods:
            method = getattr(uninitialised_sim, method_name)
            with self.assertRaises(SystemError):
                method()

    def test_verify(self):
        sim = make_sim(10)
        sim.verify()
        with self.assertRaises(TypeError):
            sim.verify("asdg")

    def test_fenwick_drift(self):
        sim = make_sim(10)
        self.assertEqual(sim.fenwick_drift(0), 0)
        self.assertRaises(TypeError, sim.fenwick_drift, "sdf")
        for bad_label in [-1, 1, 100]:
            self.assertRaises(ValueError, sim.fenwick_drift, bad_label)

    def test_discrete_genome(self):
        def f(discrete_genome):
            return make_sim(10, discrete_genome=discrete_genome)

        for discrete_genome in [True, False]:
            sim = f(discrete_genome)
            self.assertEqual(sim.discrete_genome, discrete_genome)
        self.assertEqual(sim.fenwick_drift(0), 0)
        self.assertRaises(TypeError, f, "sdf")

    def test_ploidy(self):
        def f(ploidy):
            return make_sim(10, ploidy=ploidy)

        for ploidy in [1, 2, 3]:
            sim = f(ploidy)
            self.assertEqual(sim.ploidy, ploidy)
        for bad_type in ["sdf", [], 0.0]:
            with self.assertRaises(TypeError):
                f(bad_type)
        for bad_ploidy in [-1, -100, 0]:
            with self.assertRaises(_msprime.InputError):
                f(bad_ploidy)

    @unittest.skipIf(IS_WINDOWS, "windows IO is weird")
    def test_print_state_errors(self):
        sim = make_sim(10)
        with self.assertRaises(TypeError):
            sim.print_state()
        with self.assertRaises(io.UnsupportedOperation):
            sim.print_state(io.StringIO())
        with tempfile.TemporaryDirectory() as tmpdir:
            path = pathlib.Path(tmpdir) / "testfile"
            # Write something into the file
            with open(path, "w") as f:
                print("something", file=f)
            with open(path) as f:
                with self.assertRaises(OSError):
                    sim.print_state(f)

    def test_recombination_map(self):
        def f(recomb_map):
            return make_sim(
                2,
                sequence_length=recomb_map["position"][-1],
                recombination_map=recomb_map,
            )

        maps = [
            uniform_rate_map(),
            uniform_rate_map(100, 0.001),
            {"position": np.arange(100), "rate": np.arange(99)},
            {"position": np.arange(1000), "rate": np.ones(999)},
        ]
        for recomb_map in maps:
            sim = f(recomb_map)
            for _ in range(2):
                other_map = sim.recombination_map
                self.assertEqual(len(other_map), 2)
                self.assertEqual(other_map.keys(), recomb_map.keys())
                for key, array in recomb_map.items():
                    self.assertTrue(np.array_equal(array, other_map[key]))

    def test_bad_recombination_map(self):
        def f(recomb_map):
            return make_sim(2, recombination_map=recomb_map)

        for bad_type in ["", None, []]:
            with self.assertRaises(TypeError):
                f(bad_type)
        for missing_key in [{}, {"position": []}, {"rate": []}]:
            with self.assertRaises(ValueError):
                f(missing_key)
        for bad_array in ["sdrf", b"sxdf", None, [[], []]]:
            with self.assertRaises(ValueError):
                f({"position": bad_array, "rate": []})
            with self.assertRaises(ValueError):
                f({"position": [], "rate": bad_array})

        with self.assertRaises(ValueError):
            f({"position": [0, 1], "rate": []})
        with self.assertRaises(_msprime.InputError):
            f({"position": [0], "rate": []})
        with self.assertRaises(_msprime.InputError):
            f({"position": [1, 0], "rate": [0]})
        with self.assertRaises(_msprime.InputError):
            f({"position": [0, -1], "rate": [0]})
        with self.assertRaises(_msprime.InputError):
            f({"position": [0, 1], "rate": [-1]})

    def test_bad_parameters(self):
        tables = _msprime.LightweightTableCollection(1)
        rng = _msprime.RandomGenerator(1)

        # tables and random_generator are mandatory.
        self.assertRaises(TypeError, _msprime.Simulator)
        self.assertRaises(TypeError, _msprime.Simulator, tables)
        # Empty tables raises and input error
        self.assertRaises(
            _msprime.InputError, _msprime.Simulator, tables, random_generator=rng
        )

        # check types
        for bad_type in ["1", None, {}, int]:
            self.assertRaises(
                TypeError, _msprime.Simulator, tables=bad_type, random_generator=rng
            )
            self.assertRaises(
                TypeError, _msprime.Simulator, tables=tables, random_generator=bad_type
            )
            self.assertRaises(TypeError, make_sim, recombination_rate=bad_type)
            self.assertRaises(TypeError, make_sim, avl_node_block_size=bad_type)
            self.assertRaises(TypeError, make_sim, segment_block_size=bad_type)
            self.assertRaises(TypeError, make_sim, node_mapping_block_size=bad_type)
            self.assertRaises(TypeError, make_sim, start_time=bad_type)
            self.assertRaises(TypeError, make_sim, num_labels=bad_type)
            self.assertRaises(TypeError, make_sim, gene_conversion_rate=bad_type)
            self.assertRaises(
                TypeError, make_sim, gene_conversion_track_length=bad_type
            )
        # Check for bad values.
        self.assertRaises(_msprime.InputError, make_sim, avl_node_block_size=0)
        self.assertRaises(_msprime.InputError, make_sim, segment_block_size=0)
        self.assertRaises(_msprime.InputError, make_sim, node_mapping_block_size=0)
        self.assertRaises(_msprime.InputError, make_sim, num_labels=0)
        self.assertRaises(_msprime.InputError, make_sim, num_labels=-1)
        self.assertRaises(_msprime.InputError, make_sim, gene_conversion_rate=-1)
        # Track length is ignored if gene_conversion_rate is 0
        self.assertRaises(
            _msprime.InputError,
            make_sim,
            gene_conversion_rate=1,
            gene_conversion_track_length=-100,
        )

        # Check for other type specific errors.
        self.assertRaises(OverflowError, make_sim, avl_node_block_size=2 ** 65)

    def test_num_labels(self):
        for num_labels in range(1, 10):
            sim = make_sim(5, num_labels=num_labels)
            self.assertEqual(sim.num_labels, num_labels)
            sim.run()
            self.assertEqual(sim.num_labels, num_labels)

    def test_bad_run_args(self):
        sim = get_example_simulator(10)
        for bad_type in ["x", []]:
            self.assertRaises(TypeError, sim.run, bad_type)
        for bad_value in [-1, -1e-6]:
            self.assertRaises(ValueError, sim.run, bad_value)

    def test_migrations(self):
        sim = get_example_simulator(10, num_populations=3, store_migrations=True)
        sim.run()
        tables = tskit.TableCollection.fromdict(sim.tables.asdict())
        self.assertGreater(sim.num_migrations, 0)
        self.assertEqual(sim.num_migrations, len(tables.migrations))

    def test_non_parametric_simulation_models(self):
        for bad_type in [0, None, str]:
            self.assertRaises(TypeError, make_sim, model=bad_type)
        for bad_dict in [{}, {"noname": 1}]:
            self.assertRaises(ValueError, make_sim, model=bad_dict)
        for bad_model in ["", "SMC", "ABC", "hud", None, 1234, {}, []]:
            self.assertRaises(
                ValueError, make_sim, model=get_simulation_model(bad_model)
            )
        for name in ["hudson", "smc", "smc_prime"]:
            model = get_simulation_model(name)
            sim = make_sim(model=model)
            self.assertEqual(sim.model, model)

    def test_dirac_simulation_model(self):

        for bad_type in [None, str, "sdf"]:
            model = get_simulation_model("dirac", psi=bad_type, c=1.0)
            self.assertRaises(TypeError, make_sim, model=model)
            model = get_simulation_model("dirac", psi=0.5, c=bad_type)
            self.assertRaises(TypeError, make_sim, model=model)
        self.assertRaises(ValueError, make_sim, model=get_simulation_model("dirac"))
        for bad_psi in [-1, 0, -1e-6, 1, 1e6]:
            self.assertRaises(
                ValueError,
                make_sim,
                model=get_simulation_model("dirac", c=1, psi=bad_psi),
            )
        for bad_c in [-1, -1e-6]:
            self.assertRaises(
                ValueError,
                make_sim,
                model=get_simulation_model("dirac", psi=0.5, c=bad_c),
            )
        for psi in [0.99, 0.2, 1e-4]:
            for c in [5.0, 1e2, 1e-4]:
                model = get_simulation_model("dirac", psi=psi, c=c)
                sim = make_sim(model=model)
                self.assertEqual(sim.model, model)

    def test_beta_simulation_model(self):
        for bad_type in [None, str, "sdf"]:
            model = get_simulation_model("beta", alpha=bad_type, truncation_point=1)
            self.assertRaises(TypeError, make_sim, model=model)
            model = get_simulation_model("beta", alpha=1, truncation_point=bad_type)
            self.assertRaises(TypeError, make_sim, model=model)
        model = get_simulation_model("beta", alpha=1)
        self.assertRaises(ValueError, make_sim, model=model)
        model = get_simulation_model("beta", truncation_point=1)
        self.assertRaises(ValueError, make_sim, model=model)
        # should have 1 < alpha < 2 and 0 < truncation_point <= 1
        for alpha in np.arange(1.01, 2, 0.01):
            for truncation_point in np.arange(0, 1, 0.01) + 0.01:
                model = get_simulation_model(
                    "beta", alpha=alpha, truncation_point=truncation_point
                )
                sim = make_sim(model=model)
                self.assertEqual(sim.model, model)
        # bad values
        for alpha in (-1e9, -1, 0, 1, 2, 5, 1e9):
            model = get_simulation_model("beta", alpha=alpha, truncation_point=1)
            with self.assertRaises(_msprime.InputError):
                sim = make_sim(model=model)
        for truncation_point in [-1e9, -1, 0, 1.1, 2, 1e9]:
            model = get_simulation_model(
                "beta", alpha=1.5, truncation_point=truncation_point
            )
            with self.assertRaises(_msprime.InputError):
                sim = make_sim(model=model)

    def test_sweep_genic_selection_simulation_model_errors(self):
        L = 10

        def f(num_samples=10, **kwargs):
            return make_sim(
                num_samples,
                sequence_length=L,
                recombination_map=uniform_rate_map(L=L),
                **kwargs,
            )

        # Partialy specified models.
        model = get_simulation_model("sweep_genic_selection")
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model("sweep_genic_selection", position=0.5)
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection", position=0.5, start_frequency=0.1
        )
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection",
            position=0.5,
            start_frequency=0.1,
            end_frequency=0.5,
        )
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection",
            position=0.5,
            start_frequency=0.1,
            end_frequency=0.5,
            alpha=0.1,
        )
        self.assertRaises(ValueError, f, model=model)
        model = get_simulation_model(
            "sweep_genic_selection",
            position=0.5,
            start_frequency=0.1,
            end_frequency=0.5,
            dt=0.1,
        )
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

        def f(num_samples=10, **kwargs):
            return make_sim(
                num_samples,
                sequence_length=L,
                recombination_map=uniform_rate_map(L=L),
                **kwargs,
            )

        for j in range(L):
            model = get_sweep_genic_selection_model(position=j)
            sim = f(model=model)
            model = sim.model
            # Because we have a flat map we should convert exactly.
            self.assertEqual(j, model["locus"])

    def test_sweep_after_coalescence(self):
        sim = make_sim(
            10,
            sequence_length=10,
            recombination_map=uniform_rate_map(L=10, rate=1),
            num_labels=2,
        )
        status = sim.run()
        self.assertEqual(status, _msprime.EXIT_COALESCENCE)
        t_before = sim.time
        for _ in range(100):
            model = get_sweep_genic_selection_model(position=5)
            sim.model = model
            self.assertEqual(sim.run(), _msprime.EXIT_COALESCENCE)
            self.assertEqual(t_before, sim.time)

    def test_store_migrations(self):
        def f(num_samples=10, **kwargs):
            samples = [(j % 2, 0) for j in range(num_samples)]
            migration_matrix = [[0, 1], [1, 0]]
            population_configuration = [
                get_population_configuration() for j in range(2)
            ]
            return make_sim(
                samples,
                num_populations=2,
                population_configuration=population_configuration,
                migration_matrix=migration_matrix,
                **kwargs,
            )

        for bad_type in [[], "False", None, {}, str]:
            self.assertRaises(TypeError, f, store_migrations=bad_type)
        for store_records in [True, False]:
            sim = f(store_migrations=store_records)
            self.assertEqual(sim.store_migrations, store_records)
            sim.run()
            if store_records:
                self.assertGreater(sim.num_migrations, 0)
            else:
                self.assertEqual(sim.num_migrations, 0)
            tables = tskit.TableCollection.fromdict(sim.tables.asdict())
            self.assertEqual(len(tables.migrations), sim.num_migrations)
            for migration in tables.migrations:
                self.assertTrue(0 <= migration.left < sim.sequence_length)
                self.assertTrue(0 < migration.right <= sim.sequence_length)
                self.assertTrue(0 <= migration.source < 2)
                self.assertTrue(0 <= migration.dest < 2)
                self.assertGreaterEqual(migration.node, 0)
                self.assertGreater(migration.time, 0)
        # By default, this should be off.
        sim = f()
        self.assertFalse(sim.store_migrations)

    def test_deleting_tables(self):
        rng = _msprime.RandomGenerator(1)
        tables = make_minimal_tables()
        sim = _msprime.Simulator(tables, rng)
        del tables
        sim.run()

    def test_deleting_rng(self):
        rng = _msprime.RandomGenerator(1)
        tables = make_minimal_tables()
        sim = _msprime.Simulator(tables, rng)
        del rng
        sim.run()

    def test_bad_population_configurations(self):
        def f(population_configuration):
            return make_sim(
                num_populations=len(population_configuration),
                population_configuration=population_configuration,
            )

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
                TypeError, f, [get_population_configuration(initial_size=bad_number)]
            )
            self.assertRaises(
                TypeError, f, [get_population_configuration(growth_rate=bad_number)]
            )
        # Cannot have negative for initial_size
        for bad_size in [-1, -1e300]:
            self.assertRaises(
                _msprime.InputError,
                f,
                [get_population_configuration(initial_size=bad_size)],
            )

    def test_get_population_configurations(self):
        def f(num_samples, conf_tuples):
            population_configuration = [
                get_population_configuration(initial_size=p, growth_rate=a)
                for p, a in conf_tuples
            ]
            N = len(population_configuration)
            migration_matrix = [[0 for j in range(N)] for k in range(N)]
            s = make_sim(
                num_samples,
                num_populations=N,
                population_configuration=population_configuration,
                migration_matrix=migration_matrix,
            )
            conf_dicts = s.population_configuration
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
                get_population_configuration() for j in range(num_populations)
            ]
            return make_sim(
                2,
                num_populations=num_populations,
                population_configuration=population_configuration,
                migration_matrix=migration_matrix,
            )

        for bad_type in ["", {}, None, 2, [""], [[]], [None]]:
            self.assertRaises(ValueError, f, 1, bad_type)
        for bad_value in [[1, 2], [-1], [1, 2, 3]]:
            self.assertRaises(ValueError, f, 1, bad_value)

        # Providing the wrong number of populations provokes a ValueError
        self.assertRaises(ValueError, f, 1, np.zeros((2, 2)))
        self.assertRaises(ValueError, f, 2, np.zeros((1, 1)))
        self.assertRaises(ValueError, f, 2, np.zeros((3, 3)))
        # Negative values also provoke ValueError
        self.assertRaises(_msprime.InputError, f, 2, [[0, 1], [-1, 0]])

        bad_matrices = [
            # Non-zero diagonal gives a InputError
            [[1, 1], [1, 1]],
            [[0, 1], [1, 1]],
            [[0, 1, 1], [1, 0, 1], [1, 0, 1]],
        ]
        for matrix in bad_matrices:
            num_populations = len(matrix[0])
            self.assertRaises(_msprime.InputError, f, num_populations, matrix)
        # Must provide migration_matrix when a population config is
        # provided.
        for N in range(1, 5):
            pop_conf = [get_population_configuration(2)] + [
                get_population_configuration(0) for j in range(N - 1)
            ]
            with self.assertRaises(ValueError):
                make_sim(
                    2, num_populations=N, population_configuration=pop_conf,
                )

    def test_get_migration_matrix(self):
        for N in range(1, 10):
            population_configuration = [get_population_configuration(2)] + [
                get_population_configuration(0) for _ in range(N - 1)
            ]
            random_matrix = [
                [random.random() * (j != k) for j in range(N)] for k in range(N)
            ]
            # Deliberately stress the JSON encoding code.
            nasty_matrix = [
                [random.random() * 1e9 * (j != k) for j in range(N)] for k in range(N)
            ]
            matrices = [random_matrix, nasty_matrix]
            for migration_matrix in matrices:
                sim = make_sim(
                    num_populations=N,
                    migration_matrix=migration_matrix,
                    population_configuration=population_configuration,
                )
                self.assertTrue(np.array_equal(migration_matrix, sim.migration_matrix))

    def test_bad_demographic_event_types(self):
        def f(events):
            return make_sim(demographic_events=events)

        event_generators = [
            get_size_change_event,
            get_growth_rate_change_event,
            get_migration_rate_change_event,
            get_mass_migration_event,
            get_simple_bottleneck_event,
            get_instantaneous_bottleneck_event,
        ]
        for bad_type in [None, {}, "", 1]:
            self.assertRaises(TypeError, f, bad_type)
        for bad_type in [None, "", 1, []]:
            self.assertRaises(TypeError, f, [bad_type])
        for bad_type in [None, [], 0]:
            for generator in event_generators:
                event = generator()
                event["type"] = bad_type
                self.assertRaises(ValueError, f, [event])
        for bad_event in [b"", b"1", b"x" * 1000, b"Size_change", 2, "none"]:
            for generator in event_generators:
                event = generator()
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
            event = get_mass_migration_event(source=0, dest=1, proportion=bad_type)
            self.assertRaises(TypeError, f, [event])

            event = get_migration_rate_change_event(source=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_migration_rate_change_event(dest=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_migration_rate_change_event(
                source=0, dest=1, migration_rate=bad_type
            )
            self.assertRaises(TypeError, f, [event])

            event = get_growth_rate_change_event(population=bad_type)
            self.assertRaises(TypeError, f, [event])
            event = get_growth_rate_change_event(population=0, growth_rate=bad_type)
            self.assertRaises(TypeError, f, [event])

    def test_bad_demographic_event_values(self):
        def f(events, num_populations=1):
            population_configuration = [get_population_configuration(2)] + [
                get_population_configuration(0) for _ in range(num_populations - 1)
            ]
            return make_sim(
                num_populations=num_populations,
                demographic_events=events,
                population_configuration=population_configuration,
                migration_matrix=get_migration_matrix(num_populations),
            )

        event_generators = [
            get_size_change_event,
            get_growth_rate_change_event,
            get_migration_rate_change_event,
            get_mass_migration_event,
            get_simple_bottleneck_event,
            get_instantaneous_bottleneck_event,
        ]
        for event_generator in event_generators:
            # Negative times not allowed.
            event = event_generator(time=-1)
            self.assertRaises(ValueError, f, [event])
        event_generators = [get_size_change_event, get_growth_rate_change_event]
        # Check for bad population ids.
        for event_generator in event_generators:
            for bad_pop_id in [-2, 1, 10 ** 6]:
                event = event_generator(population=bad_pop_id)
                self.assertRaises(_msprime.InputError, f, [event])
            for k in range(1, 4):
                event = event_generator(population=k)
                self.assertRaises(_msprime.InputError, f, [event], k)
                events = [event_generator(), event]
                self.assertRaises(_msprime.InputError, f, events, k)
        for bad_pop_id in [-2, 1, 10 ** 6]:
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
        event = get_migration_rate_change_event(source=-1, dest=1)
        self.assertRaises(_msprime.InputError, f, [event])
        event = get_migration_rate_change_event(source=-1, dest=1)
        self.assertRaises(_msprime.InputError, f, [event])
        for k in range(1, 4):
            # Diagonal matrix values are not accepted.
            event = get_migration_rate_change_event(source=0, dest=0)
            self.assertRaises(_msprime.InputError, f, [event], k)

        # Tests specific for mass_migration and bottleneck
        event = get_mass_migration_event(source=0, dest=0)
        self.assertRaises(_msprime.InputError, f, [event])
        for bad_proportion in [-1, 1.1, 1e7]:
            event = get_mass_migration_event(proportion=bad_proportion)
            self.assertRaises(_msprime.InputError, f, [event])
            event = get_simple_bottleneck_event(proportion=bad_proportion)
            self.assertRaises(_msprime.InputError, f, [event])

    def test_bad_demographic_event_messages(self):
        def f(events, num_populations=1):
            population_configuration = [get_population_configuration(2)] + [
                get_population_configuration(0) for _ in range(num_populations - 1)
            ]
            return make_sim(
                num_populations=num_populations,
                demographic_events=events,
                population_configuration=population_configuration,
                migration_matrix=get_migration_matrix(num_populations),
            )

        # Make sure we're getting the right index in error message.
        events = [get_size_change_event(), get_size_change_event(size=-1)]
        try:
            f(events)
        except _msprime.InputError as e:
            message = str(e)
        self.assertTrue(message.startswith("Input error in demographic_events[1]"))

    def test_unsorted_demographic_events(self):
        event_generators = [
            get_size_change_event,
            get_growth_rate_change_event,
            get_migration_rate_change_event,
            get_mass_migration_event,
            get_simple_bottleneck_event,
        ]
        events = []
        for event_generator in event_generators:
            for _ in range(3):
                events.append(event_generator(time=random.random()))
        sorted_events = sorted(events, key=lambda x: x["time"])
        self.assertNotEqual(events, sorted_events)
        with self.assertRaises(_msprime.InputError):
            make_sim(demographic_events=events)

    def test_seed_equality(self):
        simulations = [
            {"samples": get_samples(10)},
            {
                "samples": get_samples(100),
                "demographic_events": [get_simple_bottleneck_event(0.01, 0, 1.0)],
            },
            {
                "samples": get_samples(10),
                "sequence_length": 10,
                "recombination_map": uniform_rate_map(L=10, rate=1),
            },
            {
                "samples": get_population_samples(3, 3, 4),
                "num_populations": 3,
                "population_configuration": [
                    get_population_configuration(),
                    get_population_configuration(),
                    get_population_configuration(),
                ],
                "migration_matrix": get_migration_matrix(3, 3),
                "demographic_events": [
                    get_size_change_event(0.1, 0.5),
                    get_migration_rate_change_event(0.2, 2),
                    get_growth_rate_change_event(0.3, 2),
                    get_mass_migration_event(0.4, 0, 1, 0.5),
                    get_simple_bottleneck_event(0.5, 0.5),
                ],
            },
        ]
        seed = 10
        for params in simulations:
            sim1 = make_sim(random_seed=seed, **params)
            sim2 = make_sim(random_seed=seed, **params)
            sim1.run()
            sim2.run()
            tables1 = tskit.TableCollection.fromdict(sim1.tables.asdict())
            tables2 = tskit.TableCollection.fromdict(sim2.tables.asdict())
            self.assertEqual(sim1.num_edges, sim2.num_edges)
            self.assertEqual(tables1.edges, tables2.edges)
            self.assertEqual(sim1.num_nodes, sim2.num_nodes)
            self.assertEqual(tables1.nodes, tables2.nodes)
            self.assertEqual(sim1.breakpoints, sim2.breakpoints)

    def test_infinite_waiting_time(self):
        # With no migration we should have an infinite waiting time.
        population_configuration = [
            get_population_configuration(),
            get_population_configuration(),
        ]
        sim = make_sim(
            get_population_samples(5, 5),
            num_populations=2,
            population_configuration=population_configuration,
            migration_matrix=[[0, 0], [0, 0]],
        )
        self.assertRaises(_msprime.LibraryError, sim.run)

    def test_simple_event_counters(self):
        for n in [2, 10, 20]:
            sim = make_sim(n)
            sim.run()
            self.assertEqual(n - 1, sim.num_common_ancestor_events)
            self.assertEqual(0, sim.num_recombination_events)
            self.assertEqual([0], sim.num_migration_events)

    def test_simple_migration_event_counters(self):
        n = 10
        # No migration at all
        population_configuration = [
            get_population_configuration(0),
            get_population_configuration(0),
        ]
        sim = make_sim(
            get_samples(n),
            num_populations=2,
            population_configuration=population_configuration,
            migration_matrix=[[0.0, 0.0], [0.0, 0.0]],
        )
        sim.run()
        self.assertEqual(n - 1, sim.num_common_ancestor_events)
        self.assertEqual(0, sim.num_recombination_events)
        self.assertTrue(np.all(sim.num_migration_events == 0))
        # Migration between only pops 0 and 1
        matrix = np.array([[0, 5, 0], [5, 0, 0], [0, 0, 0]])
        population_configuration = [
            get_population_configuration(0),
            get_population_configuration(0),
            get_population_configuration(0),
        ]
        sim = make_sim(
            get_population_samples(5, 5),
            num_populations=3,
            population_configuration=population_configuration,
            migration_matrix=matrix,
        )
        sim.run()
        self.assertEqual(n - 1, sim.num_common_ancestor_events)
        self.assertEqual(0, sim.num_recombination_events)
        migration_events = sim.num_migration_events
        self.assertTrue(np.all(migration_events[matrix == 0] == 0))
        self.assertTrue(np.all(migration_events[matrix > 0] > 0))

    def test_large_migration_matrix_counters(self):
        # Put in linear migration into a larger matrix of populations.
        n = 10
        num_populations = 10
        migration_matrix = np.zeros((num_populations, num_populations))
        active_pops = [3, 5, 7]
        for index in range(len(active_pops) - 1):
            j = active_pops[index]
            k = active_pops[index + 1]
            migration_matrix[j][k] = 2
            migration_matrix[k][j] = 2
        population_configuration = [
            get_population_configuration() for _ in range(num_populations)
        ]
        num_sampless = [0 for _ in range(num_populations)]
        num_sampless[active_pops[0]] = 2
        num_sampless[active_pops[1]] = 2
        num_sampless[active_pops[2]] = 6
        sim = make_sim(
            get_population_samples(*num_sampless),
            num_populations=num_populations,
            population_configuration=population_configuration,
            migration_matrix=migration_matrix,
        )
        sim.run()
        self.assertEqual(n - 1, sim.num_common_ancestor_events)
        self.assertEqual(0, sim.num_recombination_events)
        migration_events = sim.num_migration_events
        self.assertEqual(migration_matrix.shape, migration_events.shape)
        self.assertTrue(np.all(migration_events[migration_matrix == 0] == 0))
        self.assertTrue(np.all(migration_events[migration_matrix > 0] > 0))

    def test_mass_migration(self):
        n = 10
        t = 0.01
        dt = 0.0000001
        sim = make_sim(
            samples=n,
            num_populations=2,
            population_configuration=[
                get_population_configuration(),
                get_population_configuration(),
            ],
            demographic_events=[
                get_migration_rate_change_event(t),
                get_mass_migration_event(t + dt, source=0, dest=1, proportion=1),
            ],
            migration_matrix=[[0, 0], [0, 0]],
        )
        sim.run(t)
        pop_sizes_before = [0, 0]
        for ind in sim.ancestors:
            for _, _, _, pop_id in ind:
                pop_sizes_before[pop_id] += 1
        sim.run(t + 2 * dt)
        pop_sizes_after = [0, 0]
        for ind in sim.ancestors:
            for _, _, _, pop_id in ind:
                pop_sizes_after[pop_id] += 1
        self.assertEqual(pop_sizes_before[0], pop_sizes_after[1])

    def test_bottleneck(self):
        n = 10
        t1 = 0.01
        t2 = 0.02
        t3 = 0.03
        dt = 1e-9
        sim = make_sim(
            get_population_samples(n, n),
            num_populations=2,
            population_configuration=[
                get_population_configuration(),
                get_population_configuration(),
            ],
            demographic_events=[
                get_simple_bottleneck_event(t1, population=0, proportion=1),
                get_simple_bottleneck_event(t2, population=1, proportion=1),
                get_mass_migration_event(t3, source=0, dest=1, proportion=1),
            ],
            migration_matrix=[[0, 0], [0, 0]],
        )
        sim.run(t1 + dt)
        pop_sizes = [0, 0]
        for ind in sim.ancestors:
            for _, _, _, pop_id in ind:
                pop_sizes[pop_id] += 1
        self.assertEqual(pop_sizes[0], 1)
        sim.run(t2 + dt)
        pop_sizes = [0, 0]
        for ind in sim.ancestors:
            for _, _, _, pop_id in ind:
                pop_sizes[pop_id] += 1
        self.assertEqual(pop_sizes[1], 1)
        sim.run()
        self.assertGreater(sim.time, t3)

    def test_recombination_event_counters(self):
        n = 10
        # No migration
        population_configuration = [
            get_population_configuration(0),
            get_population_configuration(0),
        ]
        sim = make_sim(
            get_samples(n),
            num_populations=2,
            sequence_length=10,
            recombination_map=uniform_rate_map(L=10, rate=10),
            population_configuration=population_configuration,
            migration_matrix=[[0.0, 0.0], [0.0, 0.0]],
        )
        sim.run()
        self.assertLessEqual(n - 1, sim.num_common_ancestor_events)
        self.assertLess(0, sim.num_recombination_events)
        self.assertTrue(np.array_equal(np.zeros((2, 2)), sim.num_migration_events))

    def test_single_sink_population_counters(self):
        n = 10
        # Migration only into population 2.
        matrix = np.array([[0, 0, 1], [0, 0, 1], [0, 0, 0]])
        population_configuration = [
            get_population_configuration(),
            get_population_configuration(),
            get_population_configuration(),
        ]
        sim = make_sim(
            get_population_samples(5, 5, 0),
            num_populations=3,
            population_configuration=population_configuration,
            migration_matrix=matrix,
        )
        sim.run()
        self.assertEqual(n - 1, sim.num_common_ancestor_events)
        self.assertEqual(0, sim.num_recombination_events)
        migration_events = sim.num_migration_events
        self.assertTrue(np.all(migration_events[matrix == 0] == 0))
        self.assertTrue(np.all(migration_events[matrix > 0] > 0))

    def test_reset(self):
        sim = make_sim(10)
        times = set()
        for _ in range(10):
            sim.run()
            t = sim.time
            self.assertNotIn(t, times)
            times.add(t)
            sim.reset()
            self.assertEqual(sim.time, 0)


class TestRandomGenerator(unittest.TestCase):
    """
    Tests for the random generator class.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.RandomGenerator)
        for bad_type in ["x", 1.0, {}]:
            self.assertRaises(TypeError, _msprime.RandomGenerator, bad_type)
        for bad_value in [-1, 0, 2 ** 32]:
            self.assertRaises(ValueError, _msprime.RandomGenerator, bad_value)

    def test_seed(self):
        for s in [1, 10, 2 ** 32 - 1]:
            rng = _msprime.RandomGenerator(s)
            self.assertEqual(rng.get_seed(), s)

    def test_uninitialised(self):
        uninitialised_rng = _msprime.RandomGenerator.__new__(_msprime.RandomGenerator)
        with self.assertRaises(SystemError):
            uninitialised_rng.get_seed()
        with self.assertRaises(SystemError):
            uninitialised_rng.flat()
        with self.assertRaises(SystemError):
            uninitialised_rng.poisson()
        with self.assertRaises(SystemError):
            uninitialised_rng.uniform_int()

    def test_flat_errors(self):
        rng = _msprime.RandomGenerator(1)
        self.assertRaises(TypeError, rng.flat)
        self.assertRaises(TypeError, rng.flat, 0)
        for bad_type in ["as", [], None]:
            self.assertRaises(TypeError, rng.flat, bad_type, 1)
            self.assertRaises(TypeError, rng.flat, 0, bad_type)

    def test_flat(self):
        for seed in [1, 2, 2 ** 32 - 1]:
            rng1 = _msprime.RandomGenerator(seed)
            rng2 = _msprime.RandomGenerator(seed)
            values = [0, 1, 10, -10, 1e200, -1e200]
            for a, b in itertools.product(values, repeat=2):
                self.assertEqual(rng1.flat(a, b), rng2.flat(a, b))

    def test_poisson_errors(self):
        rng = _msprime.RandomGenerator(1)
        self.assertRaises(TypeError, rng.poisson)
        for bad_type in ["as", [], None]:
            self.assertRaises(TypeError, rng.poisson, bad_type)

    def test_poisson(self):
        for seed in [1, 2, 2 ** 32 - 1]:
            rng1 = _msprime.RandomGenerator(seed)
            rng2 = _msprime.RandomGenerator(seed)
            values = [0.001, 1e-6, 0, 1, 10, -10, 100]
            for mu in values:
                self.assertEqual(rng1.poisson(mu), rng2.poisson(mu))

    def test_uniform_int_errors(self):
        rng = _msprime.RandomGenerator(1)
        self.assertRaises(TypeError, rng.uniform_int)
        for bad_type in ["as", [], None]:
            self.assertRaises(TypeError, rng.uniform_int, bad_type)

    def test_uniform_int(self):
        for seed in [1, 2, 2 ** 32 - 1]:
            rng1 = _msprime.RandomGenerator(seed)
            rng2 = _msprime.RandomGenerator(seed)
            values = [-1, 0, 1, 2, 10, 100, 2 ** 31]
            for n in values:
                self.assertEqual(rng1.uniform_int(n), rng2.uniform_int(n))


class TestMatrixMutationModel(unittest.TestCase):
    """
    Tests for the mutation model class.
    """

    def test_constructor(self):
        self.assertRaises(TypeError, _msprime.MatrixMutationModel)
        self.assertRaises(TypeError, _msprime.MatrixMutationModel, [])
        self.assertRaises(TypeError, _msprime.MatrixMutationModel, [], [])

        for bad_type in [None, "x", 123]:
            self.assertRaises(TypeError, _msprime.MatrixMutationModel, bad_type, [], [])

    def test_bad_matrix(self):
        alleles = ["0", "1"]
        dist = [1, 0]
        bad_matrixes = [
            [],
            [[], []],
            [[0, 1], [0]],
            [[0, 1], [0, 0], [0, 0]],
            ["wwer"],
            None,
            {},
            [{}, {}],
        ]
        for bad_matrix in bad_matrixes:
            with self.assertRaises(ValueError):
                _msprime.MatrixMutationModel(alleles, dist, bad_matrix)

    def test_bad_lengths(self):
        for num_alleles in [2, 4, 6]:
            alleles = [str(j) for j in range(num_alleles)]
            distribution = np.zeros(num_alleles)
            matrix = np.zeros((num_alleles, num_alleles))
            self.assertRaises(
                ValueError, _msprime.MatrixMutationModel, alleles, [], matrix
            )
            self.assertRaises(
                ValueError,
                _msprime.MatrixMutationModel,
                alleles,
                distribution[:-1],
                matrix,
            )
            self.assertRaises(
                ValueError, _msprime.MatrixMutationModel, alleles, distribution, []
            )
            self.assertRaises(
                ValueError,
                _msprime.MatrixMutationModel,
                alleles,
                distribution,
                matrix[:, :-1],
            )
            self.assertRaises(
                ValueError,
                _msprime.MatrixMutationModel,
                alleles,
                distribution,
                matrix[:-1, :-1],
            )

    def test_bad_alleles(self):
        for alleles in [[b"A", b"B"], ["0", "1", b"2"], [b"x", None]]:
            n = len(alleles)
            distribution = np.zeros(n)
            matrix = np.zeros((n, n))
            with self.assertRaises(TypeError):
                _msprime.MatrixMutationModel(alleles, distribution, matrix)

        for alleles in [[], ["a"], ["asdfsadg"]]:
            n = len(alleles)
            distribution = np.zeros(n)
            matrix = np.zeros((n, n))
            with self.assertRaises(_msprime.LibraryError):
                _msprime.MatrixMutationModel(alleles, distribution, matrix)

    def test_uninitialised(self):
        model = _msprime.MatrixMutationModel.__new__(_msprime.MatrixMutationModel)
        with self.assertRaises(SystemError):
            model.root_distribution
        with self.assertRaises(SystemError):
            model.alleles
        with self.assertRaises(SystemError):
            model.transition_matrix

    def test_good_alleles(self):
        for n in range(2, 10):
            alleles = [f"{j}" for j in range(n)]
            dist = np.zeros(n)
            dist[0] = 1
            matrix = np.zeros((n, n))
            matrix[:, 0] = 1
            mm = _msprime.MatrixMutationModel(alleles, dist, matrix)
            self.assertEqual(mm.alleles, alleles)
            self.assertTrue(np.array_equal(mm.root_distribution, dist))
            self.assertTrue(np.array_equal(mm.transition_matrix, matrix))

    def test_unicode_alleles(self):
        unicode_strings = ["á", "þ÷ý", "🌳", "🎄🌳"]
        for allele in unicode_strings:
            alleles = ["", allele]
            dist = [1, 0]
            matrix = np.zeros((2, 2))
            matrix[:, 0] = 1
            mm = _msprime.MatrixMutationModel(alleles, dist, matrix)
            self.assertEqual(mm.alleles, alleles)
            self.assertTrue(np.array_equal(mm.root_distribution, dist))
            self.assertTrue(np.array_equal(mm.transition_matrix, matrix))

    def test_multichar_alleles(self):
        for n in range(2, 10):
            alleles = ["x" * j for j in range(n)]
            dist = np.zeros(n)
            dist[:] = 1 / n
            matrix = np.zeros((n, n))
            matrix[:] = 1 / n
            mm = _msprime.MatrixMutationModel(alleles, dist, matrix)
            self.assertEqual(mm.alleles, alleles)
            self.assertTrue(np.array_equal(mm.root_distribution, dist))
            self.assertTrue(np.array_equal(mm.transition_matrix, matrix))

    def test_bad_probabilities(self):
        alleles = ["0", "1"]
        matrix = np.zeros((2, 2))
        matrix[:] = 0.5
        for bad_root in [[1, 1], [-1, 2], [0, 100]]:
            with self.assertRaises(_msprime.LibraryError):
                _msprime.MatrixMutationModel(alleles, bad_root, matrix)

        matrix[:] = 1
        with self.assertRaises(_msprime.LibraryError):
            _msprime.MatrixMutationModel(alleles, [0, 1], matrix)
        matrix[:] = -1
        with self.assertRaises(_msprime.LibraryError):
            _msprime.MatrixMutationModel(alleles, [0, 1], matrix)
        matrix = [[0, 0], [1, 1]]
        with self.assertRaises(_msprime.LibraryError):
            _msprime.MatrixMutationModel(alleles, [0, 1], matrix)


class TestSLiMMutationModel(unittest.TestCase):
    """
    Tests for the slim mutation model class.
    """

    def test_constructor_errors(self):
        self.assertRaises(TypeError, _msprime.SLiMMutationModel)

        for bad_type in ["sdr", 0.222, None]:
            with self.assertRaises(TypeError):
                _msprime.SLiMMutationModel(type=bad_type, next_id=234)
                _msprime.SLiMMutationModel(type=0, next_id=bad_type)
                _msprime.SLiMMutationModel(type=0, next_id=0, block_size=bad_type)

        with self.assertRaises(_msprime.LibraryError):
            _msprime.SLiMMutationModel(type=-1, next_id=0)
            _msprime.SLiMMutationModel(type=1, next_id=-1)

    def test_uninitialised(self):
        model = _msprime.SLiMMutationModel.__new__(_msprime.SLiMMutationModel)
        with self.assertRaises(SystemError):
            model.type
        with self.assertRaises(SystemError):
            model.next_id

    def test_type(self):
        for mutation_type in [0, 10, 2 ** 31 - 1]:
            model = _msprime.SLiMMutationModel(type=mutation_type)
            self.assertEqual(model.type, mutation_type)
            self.assertEqual(model.next_id, 0)

    def test_next_id(self):
        for next_id in [0, 10, 2 ** 63 - 1]:
            model = _msprime.SLiMMutationModel(0, next_id=next_id)
            self.assertEqual(model.type, 0)
            self.assertEqual(model.next_id, next_id)


class TestInfiniteAllelesMutationModel(unittest.TestCase):
    """
    Tests for the infinite alleles mutation model class.
    """

    def test_constructor_errors(self):
        for bad_type in ["sdr", 0.222, None]:
            with self.assertRaises(TypeError):
                _msprime.InfiniteAllelesMutationModel(start_allele=bad_type)

    def test_defaults(self):
        iamm = _msprime.InfiniteAllelesMutationModel()
        self.assertEqual(iamm.start_allele, 0)
        self.assertEqual(iamm.next_allele, 0)

    def test_uninitialised(self):
        model = _msprime.InfiniteAllelesMutationModel.__new__(
            _msprime.InfiniteAllelesMutationModel
        )
        with self.assertRaises(SystemError):
            model.start_allele
        with self.assertRaises(SystemError):
            model.next_allele

    def test_start_allele(self):
        for start_allele in [0, 10, 2 ** 64 - 1]:
            iamm = _msprime.InfiniteAllelesMutationModel(start_allele)
            self.assertEqual(iamm.start_allele, start_allele)
            self.assertEqual(iamm.next_allele, start_allele)

    def test_start_allele_overflow(self):
        for start_allele in [2 ** 64, 2 ** 65 + 1]:
            iamm = _msprime.InfiniteAllelesMutationModel(start_allele)
            self.assertEqual(iamm.start_allele, start_allele % 2 ** 64)


class TestSimMutations(unittest.TestCase):
    """
    Tests for the sim_mutations function.
    """

    def test_mandatory_args(self):
        tables = _msprime.LightweightTableCollection(1.0)
        rate_map = uniform_rate_map(1, 1)
        model = get_mutation_model()
        self.assertRaises(TypeError, _msprime.sim_mutations)
        self.assertRaises(TypeError, _msprime.sim_mutations, tables)
        rng = _msprime.RandomGenerator(1)
        _msprime.sim_mutations(tables, rng, rate_map, model)
        for bad_type in [None, "x"]:
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables=bad_type,
                    random_generator=rng,
                    rate_map=rate_map,
                    model=model,
                )
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables=tables,
                    random_generator=bad_type,
                    rate_map=rate_map,
                    model=model,
                )
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables=tables, random_generator=rng, rate_map=bad_type, model=model
                )
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables=tables,
                    random_generator=rng,
                    rate_map=rate_map,
                    model=bad_type,
                )
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables=tables,
                    random_generator=rng,
                    rate_map=rate_map,
                    model=model,
                    discrete_sites=bad_type,
                )
        self.assertIsNone(
            _msprime.sim_mutations(
                tables=tables, random_generator=rng, rate_map=rate_map, model=model,
            )
        )

    def test_optional_args(self):
        tables = _msprime.LightweightTableCollection(1.0)
        rate_map = uniform_rate_map(1, 1)
        model = get_mutation_model()
        rng = _msprime.RandomGenerator(1)
        generate = functools.partial(
            _msprime.sim_mutations,
            tables=tables,
            random_generator=rng,
            rate_map=rate_map,
            model=model,
        )
        generate()
        for bad_type in [[], None, "asdf"]:
            with self.assertRaises(TypeError):
                generate(discrete_sites=bad_type)
            with self.assertRaises(TypeError):
                generate(keep=bad_type)
            with self.assertRaises(TypeError):
                generate(start_time=bad_type)
            with self.assertRaises(TypeError):
                generate(end_time=bad_type)

    def test_tables(self):
        imap = uniform_rate_map(1)
        rng = _msprime.RandomGenerator(1)
        model = get_mutation_model()
        for bad_type in ["x", {}, None]:
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables=bad_type, random_generator=rng, rate_map=imap, model=model,
                )
        tables = _msprime.LightweightTableCollection.__new__(
            _msprime.LightweightTableCollection
        )
        with self.assertRaises(SystemError):
            _msprime.sim_mutations(
                tables=tables, random_generator=rng, rate_map=imap, model=model,
            )

    def test_rng(self):
        imap = uniform_rate_map(1)
        tables = _msprime.LightweightTableCollection(1.0)
        for bad_type in ["x", {}, None]:
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables=tables,
                    random_generator=bad_type,
                    rate_map=imap,
                    model=get_mutation_model(),
                )
        rng = _msprime.RandomGenerator.__new__(_msprime.RandomGenerator)
        # This doesn't have __init__ called, so will have NULL as the rng.
        with self.assertRaises(SystemError):
            _msprime.sim_mutations(
                tables=tables,
                random_generator=rng,
                rate_map=imap,
                model=get_mutation_model(),
            )

    def test_mutation_map(self):
        rng = _msprime.RandomGenerator(1)
        tables = _msprime.LightweightTableCollection(1.0)
        for bad_type in ["x", None, [[], []]]:
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables, rng, rate_map=bad_type, model=get_mutation_model(),
                )
        for bad_mutation_map in [{}, {"position": [], "rate": []}]:
            with self.assertRaises(ValueError):
                _msprime.sim_mutations(
                    tables, rng, rate_map=bad_mutation_map, model=get_mutation_model(),
                )

    def test_model(self):
        rng = _msprime.RandomGenerator(1)
        imap = uniform_rate_map(1)
        tables = _msprime.LightweightTableCollection(1.0)
        for bad_type in ["x", {}, None, [[], []]]:
            with self.assertRaises(TypeError):
                _msprime.sim_mutations(
                    tables, rng, rate_map=imap, model=bad_type,
                )
        model_classes = [
            _msprime.SLiMMutationModel,
            _msprime.InfiniteAllelesMutationModel,
            _msprime.MatrixMutationModel,
        ]
        for cls in model_classes:
            uninitialised_model = cls.__new__(cls)
            with self.assertRaises(SystemError):
                _msprime.sim_mutations(
                    tables, rng, rate_map=imap, model=uninitialised_model,
                )

    def test_slim_model(self):
        rng = _msprime.RandomGenerator(1)
        imap = uniform_rate_map(1)
        tables = _msprime.LightweightTableCollection(1.0)
        model = _msprime.SLiMMutationModel(1234, 5678)
        _msprime.sim_mutations(tables, rng, imap, model)
        self.assertEqual(model.next_id, 5678)

    def test_infinite_alleles_model(self):
        rng = _msprime.RandomGenerator(1)
        imap = uniform_rate_map(1, 1)
        tables = _msprime.LightweightTableCollection(1.0)
        model = _msprime.InfiniteAllelesMutationModel(1234)
        _msprime.sim_mutations(tables, rng, imap, model)
        self.assertEqual(model.start_allele, 1234)
        self.assertEqual(model.next_allele, 1234)

    def test_time_interval(self):
        rng = _msprime.RandomGenerator(1)
        rate_map = uniform_rate_map(1, 1)
        tables = _msprime.LightweightTableCollection(1)
        mutgen = functools.partial(
            _msprime.sim_mutations, tables, rng, rate_map, get_mutation_model()
        )
        for bad_type in ["x", {}, None]:
            with self.assertRaises(TypeError):
                mutgen(start_time=bad_type)
            with self.assertRaises(TypeError):
                mutgen(end_time=bad_type)
        for start_time, end_time in [(1, 0), (-1, -2), (200, 100)]:
            with self.assertRaises(_msprime.LibraryError):
                mutgen(start_time=start_time, end_time=end_time)

    def verify_block_size(self, tables):
        rng = _msprime.RandomGenerator(1)
        ll_tables = _msprime.LightweightTableCollection()
        ll_tables.fromdict(tables.asdict())
        rate_map = uniform_rate_map(tables.sequence_length, 2)
        _msprime.sim_mutations(
            ll_tables, rng, rate_map, get_mutation_model(1), keep=True
        )
        after_tables = tskit.TableCollection.fromdict(ll_tables.asdict())
        self.assertEqual(tables, after_tables)

    def test_keep_mutations_block_size_mutations(self):
        tables = tskit.TableCollection(1)
        tables.nodes.add_row(0)
        tables.sites.add_row(0, "a")
        for _ in range(8192):
            tables.mutations.add_row(0, node=0, derived_state="b")
        tables.build_index()
        tables.compute_mutation_parents()
        self.verify_block_size(tables)

    def test_keep_mutations_block_size_ancestral_state(self):
        tables = tskit.TableCollection(1)
        big_alloc = 64 * 1024
        tables.nodes.add_row(0)
        tables.sites.add_row(0, "a" * big_alloc)
        self.verify_block_size(tables)

    def test_keep_mutations_block_size_derived_state(self):
        tables = tskit.TableCollection(1)
        tables.nodes.add_row(0)
        tables.sites.add_row(0, "a")
        big_alloc = 64 * 1024
        tables.mutations.add_row(0, node=0, derived_state="b" * big_alloc)
        self.verify_block_size(tables)

    def test_keep_mutations_block_size_metadata(self):
        tables = tskit.TableCollection(1)
        big_alloc = 64 * 1024
        tables.nodes.add_row(0)
        tables.sites.add_row(0, "a", metadata=b"x" * big_alloc)
        self.verify_block_size(tables)
        tables.mutations.add_row(0, node=0, derived_state="b" * 2 * big_alloc)
        self.verify_block_size(tables)


class TestDemographyDebugger(unittest.TestCase):
    """
    Tests for the demography debugging interface.
    """

    def get_simulator(self, events):
        return make_sim(samples=2, demographic_events=events,)

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


class TestLikelihood(unittest.TestCase):
    """
    Tests for the low-level likelihood calculation interface.
    """

    def get_arg(self):
        L = 20
        rate_map = uniform_rate_map(L=L, rate=2)
        sim = make_sim(
            samples=5,
            sequence_length=L,
            recombination_map=rate_map,
            store_full_arg=True,
        )
        sim.run()
        _msprime.sim_mutations(
            sim.tables, _msprime.RandomGenerator(1), rate_map, get_mutation_model()
        )
        t = tskit.TableCollection.fromdict(sim.tables.asdict())
        self.assertGreater(len(t.edges), 10)
        self.assertGreater(len(t.mutations), 10)
        return sim.tables

    def test_simple_example(self):
        tables = self.get_arg()
        log_lik = _msprime.log_likelihood_arg(tables, 1, 1)
        self.assertLess(log_lik, 0)

    def test_interface(self):
        tables = self.get_arg()
        self.assertRaises(TypeError, _msprime.log_likelihood_arg)
        self.assertRaises(TypeError, _msprime.log_likelihood_arg, tables)
        self.assertRaises(TypeError, _msprime.log_likelihood_arg, tables, 1)
        self.assertRaises(
            TypeError, _msprime.log_likelihood_arg, tables, recombination_rate=1
        )

        for bad_tables in [None, {}, "SDf"]:
            with self.assertRaises(TypeError):
                _msprime.log_likelihood_arg(bad_tables, 1, 1)

        for bad_float in [[], None, "SDF"]:
            with self.assertRaises(TypeError):
                _msprime.log_likelihood_arg(tables, bad_float, 1)
            with self.assertRaises(TypeError):
                _msprime.log_likelihood_arg(tables, 1, bad_float)

        for bad_rec_rate in [-0.1, -1]:
            with self.assertRaises(ValueError):
                _msprime.log_likelihood_arg(tables, 1, recombination_rate=bad_rec_rate)

        for bad_Ne in [0, -1]:
            with self.assertRaises(_msprime.LibraryError):
                _msprime.log_likelihood_arg(tables, bad_Ne, recombination_rate=1)

    def test_bad_tables(self):
        # Pass in a table collection that can't be made into a tree
        # sequence.
        lw_tables = _msprime.LightweightTableCollection(0)
        with self.assertRaises(_msprime.LibraryError):
            _msprime.log_likelihood_arg(lw_tables, 1, 1)
