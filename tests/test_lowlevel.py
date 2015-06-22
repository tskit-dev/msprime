#
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

import heapq
import json
import random
import os.path
import tempfile

import tests
import _msprime


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
    def get_tree_sequence(
            self, sample_size=10, num_loci=100, mutation_rate=10):
        models = get_random_population_models(3)
        sim = _msprime.Simulator(
            sample_size, random_seed=1, num_loci=num_loci,
            scaled_recombination_rate=1, population_models=models)
        sim.run()
        ts = _msprime.TreeSequence()
        ts.create(sim)
        ts.generate_mutations(mutation_rate, 1)
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
        self.assertGreater(sim.get_num_breakpoints(), 0)
        self.assertGreater(sim.get_time(), 0.0)
        self.assertGreater(sim.get_num_ancestors(), 1)
        self.assertGreater(sim.get_used_memory(), 0)
        events = sim.get_num_coancestry_events()
        events += sim.get_num_recombination_events()
        self.assertGreater(events, 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
        n = sim.get_sample_size()
        m = sim.get_num_loci()
        ancestors = sim.get_ancestors()
        self.assertEqual(len(ancestors), sim.get_num_ancestors())
        for ind in ancestors:
            for l, r, node in ind:
                self.assertTrue(0 <= l < m)
                self.assertTrue(1 <= r <= m)
                self.assertGreaterEqual(node, 1)
        breakpoints = sim.get_breakpoints()
        self.assertEqual(len(breakpoints), sim.get_num_breakpoints())
        self.assertEqual(breakpoints, sorted(breakpoints))
        self.assertEqual(breakpoints[0], 0)
        self.assertEqual(breakpoints[-1], m)
        records = sim.get_coalescence_records()
        self.assertEqual(len(records), sim.get_num_coalescence_records())
        for l, r, p, children, t in records:
            self.assertEqual(children, tuple(sorted(children)))
            self.assertTrue(0 <= l < m)
            self.assertTrue(1 <= r <= m)
            self.assertGreater(t, 0.0)
            self.assertIn(l, breakpoints)
            self.assertIn(r, breakpoints)
        # The amount of ancestral material in the coalescence records and
        # the extant segments over all intervals should be either n (if
        # the given interval has not fully coalesced yet) or n - 1 (if
        # full coalescence has occured).
        segments_am = [0 for b in breakpoints[:-1]]
        for ind in ancestors:
            for l, r, _ in ind:
                j = breakpoints.index(l)
                while breakpoints[j] < r:
                    segments_am[j] += 1
                    j += 1
        records_am = [0 for b in breakpoints[:-1]]
        for l, r, _, _, _ in records:
            j = breakpoints.index(l)
            while breakpoints[j] < r:
                records_am[j] += 1
                j += 1
        for segment_am, record_am in zip(segments_am, records_am):
            if segment_am == 0:
                self.assertEqual(record_am, n - 1)
            else:
                self.assertEqual(segment_am + record_am, n)

    def verify_trees(self, sim, sorted_records):
        """
        Verifies that the specified set of sorted coalescence records
        corresponds to correct trees for the specified simulation.
        """
        n = sim.get_sample_size()
        pi = {}
        tau = {j: 0 for j in range(1, n + 1)}
        last_l = 0
        last_t = 0
        num_trees = 0
        live_segments = []
        for l, r, node, children, t in sorted_records:
            if last_l != l:
                last_l = l
                last_t = 0
                for j in range(1, n + 1):
                    assert j in pi
                # insert the root
                v = 1
                while v in pi:
                    v = pi[v]
                pi[v] = 0
                self.verify_sparse_tree(n, pi, tau)
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
            pi[children[0]] = node
            pi[children[1]] = node
            tau[node] = t
            # Ensure that records are sorted by time within a block
            self.assertLessEqual(last_t, t)
        for j in range(1, n + 1):
            assert j in pi
        # Insert the root branch.
        v = 1
        while v in pi:
            v = pi[v]
        pi[v] = 0
        self.verify_sparse_tree(n, pi, tau)
        num_trees += 1
        self.assertLessEqual(num_trees, sim.get_num_breakpoints())

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
        self.assertGreater(sim.get_num_breakpoints(), 0)
        self.assertGreater(sim.get_num_coalescence_records(), 0)
        self.assertGreater(sim.get_time(), 0.0)
        events = sim.get_num_coancestry_events()
        events += sim.get_num_recombination_events()
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
        times = [t for l, r, node, children, t in records]
        # Children should always be sorted in order.
        for _, _, _, children, _ in records:
            self.assertEqual(children, tuple(sorted(children)))
        self.assertEqual(times, sorted(times))
        self.assertEqual(times[-1], sim.get_time())
        self.verify_squashed_records(records)
        left_sorted_records = sorted(records, key=lambda r: r[0])
        right_sorted_records = sorted(records, key=lambda r: r[1])
        self.verify_trees(sim, left_sorted_records)
        # Check the TreeSequence. Ensure we get the records back in the
        # correct orders.
        ts = _msprime.TreeSequence()
        ts.create(sim)
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
        json_str = ts.get_simulation_parameters()
        parameters = json.loads(json_str)
        self.assertEqual(parameters["random_seed"], sim.get_random_seed())
        self.assertEqual(parameters["sample_size"], sim.get_sample_size())
        self.assertEqual(parameters["num_loci"], sim.get_num_loci())
        self.assertEqual(
            parameters["scaled_recombination_rate"],
            sim.get_scaled_recombination_rate())
        models = parameters["population_models"]
        sim_models = sim.get_population_models()
        self.assertEqual(models, sim_models)

    def verify_random_parameters(self):
        mb = 1024 * 1024
        n = random.randint(2, 1000)
        m = random.randint(1, 10**6)
        rho = random.uniform(0, 1000)
        num_pop_models = random.randint(0, 10)
        models = get_random_population_models(num_pop_models)
        random_seed = random.randint(0, 2**31)
        max_memory = random.randint(10 * mb, 100 * mb)
        segment_block_size = random.randint(1, 100)
        node_mapping_block_size = random.randint(1, 100)
        avl_node_block_size = random.randint(1, 100)
        coalescence_record_block_size = random.randint(1, 100)
        sim = _msprime.Simulator(
            sample_size=n, num_loci=m, population_models=models,
            scaled_recombination_rate=rho,
            random_seed=random_seed, max_memory=max_memory,
            segment_block_size=segment_block_size,
            avl_node_block_size=avl_node_block_size,
            node_mapping_block_size=node_mapping_block_size,
            coalescence_record_block_size=coalescence_record_block_size)
        # Check initial state
        self.assertEqual(0, sim.get_num_breakpoints())
        self.assertEqual(0.0, sim.get_time())
        self.assertEqual(0, sim.get_num_ancestors())
        self.assertEqual(0, sim.get_num_coancestry_events())
        self.assertEqual(0, sim.get_num_recombination_events())
        self.assertEqual(sim.get_num_avl_node_blocks(), 0)
        self.assertEqual(sim.get_num_segment_blocks(), 0)
        self.assertEqual(sim.get_num_node_mapping_blocks(), 0)
        self.assertEqual(sim.get_num_coalescence_record_blocks(), 0)
        self.assertEqual(sim.get_sample_size(), n)
        self.assertEqual(sim.get_num_loci(), m)
        self.assertEqual(0, len(sim.get_ancestors()))
        # Force initialisaton of the simulation by running for 0 events.
        sim.run(-1.0)
        a = 0
        nodes = set()
        for ind in sim.get_ancestors():
            self.assertEqual(len(ind), 1)
            l, r, node = ind[0]
            self.assertEqual(l, 0)
            self.assertEqual(r, m)
            self.assertFalse(node in nodes)
            nodes.add(node)
            a += 1
        self.assertEqual(a, n)
        for j in range(3):
            # Check the getters to ensure we've got the right values.
            self.assertEqual(n, sim.get_sample_size())
            self.assertEqual(m, sim.get_num_loci())
            self.assertEqual(rho, sim.get_scaled_recombination_rate())
            self.assertEqual(random_seed, sim.get_random_seed())
            self.assertEqual(max_memory, sim.get_max_memory())
            self.assertEqual(segment_block_size, sim.get_segment_block_size())
            self.assertEqual(
                avl_node_block_size, sim.get_avl_node_block_size())
            self.assertEqual(
                node_mapping_block_size, sim.get_node_mapping_block_size())
            self.assertEqual(
                coalescence_record_block_size,
                sim.get_coalescence_record_block_size())
            self.verify_population_models(sim, models)
            # Run this for a tiny amount of time and check the state
            self.assertFalse(sim.run(1e-8))
            self.verify_running_simulation(sim)
            self.verify_population_models(sim, models)

    def verify_population_models(self, sim, models):
        """
        Verifies that the population models returned by the simulator
        match the specified list.
        """
        pop_models = sim.get_population_models()
        self.assertEqual(pop_models, models)

    def verify_tree_diffs(self, tree_sequence):
        n = tree_sequence.get_sample_size()
        m = tree_sequence.get_num_loci()
        # Check some basic properties of the diffs.
        diffs = list(_msprime.TreeDiffIterator(tree_sequence))
        length, records_out, records_in = diffs[0]
        self.assertGreaterEqual(length, 0)
        self.assertEqual(len(records_out), 0)
        self.assertEqual(len(records_in), n - 1)
        self.assertEqual(sum([l for l, _, _ in diffs]), m)
        for l, records_out, records_in in diffs[1:]:
            self.assertGreaterEqual(l, 0)
            self.assertGreaterEqual(l, 0)
            self.assertEqual(len(records_out), len(records_in))
            for node, children, time in records_out + records_in:
                for c in children:
                    self.assertGreater(c, 0)
                    self.assertGreater(node, c)
                    self.assertGreater(time, 0.0)
        # Compare with the Python implementation.
        pts = tests.PythonTreeSequence(tree_sequence)
        python_diffs = list(pts.diffs())
        self.assertGreaterEqual(len(python_diffs), 0)
        self.assertEqual(diffs, python_diffs)

    def verify_simulation(self, n, m, r, models):
        """
        Runs the specified simulation and verifies its state.
        """
        # These tests don't work for n == 2
        assert n > 2
        mb = 1024 * 1024
        random_seed = random.randint(0, 2**31)
        sim = _msprime.Simulator(
            sample_size=n, num_loci=m, population_models=models,
            scaled_recombination_rate=r, random_seed=random_seed,
            max_memory=10 * mb, segment_block_size=1000,
            avl_node_block_size=1000, node_mapping_block_size=1000,
            coalescence_record_block_size=1000)
        self.verify_population_models(sim, models)
        # Run the sim for a tiny amount of time and check.
        # self.assertFalse(sim.run(1e-8))
        self.assertFalse(sim.run(1e-8))
        self.verify_running_simulation(sim)
        increment = 0.01
        t = sim.get_time() + increment
        while not sim.run(t):
            self.assertGreaterEqual(sim.get_time(), t)
            self.verify_running_simulation(sim)
            t += increment
        self.verify_completed_simulation(sim)
        self.verify_population_models(sim, models)
        # Check the tree sequence.
        tree_sequence = _msprime.TreeSequence()
        tree_sequence.create(sim)
        self.verify_tree_diffs(tree_sequence)

    def test_random_sims(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            self.verify_random_parameters()

    def test_small_sims(self):
        self.verify_simulation(3, 1, 0.0, [])
        self.verify_simulation(3, 100, 0.0, [])
        self.verify_simulation(3, 10, 1000.0, [])
        self.verify_simulation(5, 10, 10.0, [])
        self.verify_simulation(10, 100, 1.0, [])
        self.verify_simulation(100, 100, 0.1, [])

    def test_population_models(self):
        exp_model = _msprime.POP_MODEL_EXPONENTIAL
        m1 = {"alpha": 0.0, "start_time": 0.0, "type": exp_model}
        m2 = {"alpha": 1.0, "start_time": 0.5, "type": exp_model}
        # TODO add constant pop models here too.
        for n in [5, 10, 100]:
            self.verify_simulation(n, 1, 0.0, [m1])
            self.verify_simulation(n, 1, 0.0, [m2])
            self.verify_simulation(n, 1, 0.0, [m1, m2])


class TestSimulator(LowLevelTestCase):
    """
    Tests for the low-level interface to the simulator.
    """
    def test_bad_parameters(self):
        def f(sample_size=10, random_seed=1, **kwargs):
            return _msprime.Simulator(sample_size, random_seed, **kwargs)
        # sample_size and random_seed are mandatory
        self.assertRaises(TypeError, _msprime.Simulator)
        self.assertRaises(TypeError, _msprime.Simulator, sample_size=10)
        self.assertRaises(TypeError, _msprime.Simulator, random_seed=10)
        # check types
        for bad_type in ["1", None, [], {}, int]:
            self.assertRaises(TypeError, f, sample_size=bad_type)
            self.assertRaises(TypeError, f, random_seed=bad_type)
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

    def test_bad_population_models(self):
        def f(population_models):
            return _msprime.Simulator(
                2, 1, population_models=population_models)
        self.assertRaises(TypeError, f, "")
        self.assertRaises(TypeError, f, [""])
        self.assertRaises(_msprime.InputError, f, [{}])
        self.assertRaises(_msprime.InputError, f, [{"start_time": 1}])
        self.assertRaises(
            _msprime.InputError, f, [{"start_time": 1, "type": -1000}])
        self.assertRaises(
            _msprime.InputError, f,
            [{"start_time": 1, "type": _msprime.POP_MODEL_EXPONENTIAL}])
        self.assertRaises(
            _msprime.InputError, f,
            [{"start_time": 1, "type": _msprime.POP_MODEL_CONSTANT}])
        m = [{
            "start_time": -1.0, "type": _msprime.POP_MODEL_EXPONENTIAL,
            "alpha": 1.0}]
        self.assertRaises(_msprime.InputError, f, m)
        m = get_random_population_models(10)
        sim = f(m)
        self.assertIsNot(sim, None)
        m.reverse()
        self.assertRaises(_msprime.InputError, f, m)


class TestTreeSequence(LowLevelTestCase):
    """
    Tests for the low-level interface for the TreeSequence.
    """

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

    def verify_mutation_parameters_json(self, json_str):
        parameters = json.loads(json_str.decode())
        self.assertIn("scaled_mutation_rate", parameters)
        self.assertIn("random_seed", parameters)
        self.assertIsInstance(parameters["scaled_mutation_rate"], float)
        self.assertIsInstance(parameters["random_seed"], int)

    def verify_tree_parameters_json(self, json_str):
        parameters = json.loads(json_str.decode())
        self.assertIn("scaled_recombination_rate", parameters)
        self.assertIn("random_seed", parameters)
        self.assertIn("sample_size", parameters)
        self.assertIn("num_loci", parameters)
        self.assertIn("population_models", parameters)
        self.assertIsInstance(parameters["scaled_recombination_rate"], float)
        self.assertIsInstance(parameters["random_seed"], int)
        self.assertIsInstance(parameters["sample_size"], int)
        self.assertIsInstance(parameters["num_loci"], int)
        self.assertIsInstance(parameters["population_models"], list)
        for m in parameters["population_models"]:
            self.assertIsInstance(m, dict)
            self.assertIn("start_time", m)
            self.assertIn("type", m)
            self.assertIsInstance(m["start_time"], float)
            t = m["type"]
            self.assertIsInstance(t, int)
            self.assertIn(t, [
                _msprime.POP_MODEL_EXPONENTIAL, _msprime.POP_MODEL_CONSTANT])
            if t == _msprime.POP_MODEL_EXPONENTIAL:
                self.assertIn("alpha", m)
                self.assertIsInstance(m["alpha"], float)
            else:
                self.assertIn("size", m)
                self.assertIsInstance(m["size"], float)

    def verify_environment_json(self, json_str):
        environment = json.loads(json_str.decode())
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
        self.assertLessEqual(keys, set(["mutations", "trees"]))
        self.assertIn("trees", keys)
        if ts.get_num_mutations() == 0:
            self.assertNotIn("mutations", keys)
        else:
            self.assertIn("mutations", keys)
        # Check the basic root attributes
        format_version = root.attrs['format_version']
        self.assertEqual(format_version[0], 0)
        self.assertEqual(format_version[1], 1)
        self.assertEqual(root.attrs["sample_size"], ts.get_sample_size())
        self.assertEqual(root.attrs["num_loci"], ts.get_num_loci())

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
            ("left", uint32, 1), ("right", uint32, 1),
            ("node", uint32, 1), ("children", uint32, 2),
            ("time", float64, 1)]
        self.assertEqual(set(g.keys()), set([name for name, _, _ in fields]))
        for name, dtype, dims in fields:
            self.assertEqual(len(g[name].shape), dims)
            self.assertEqual(g[name].shape[0], ts.get_num_records())
            if dims == 2:
                self.assertEqual(g[name].shape[1], 2)
            self.assertEqual(g[name].dtype, dtype)

    def test_dump_format(self):
        for ts in self.get_example_tree_sequences():
            with tempfile.NamedTemporaryFile() as f:
                self.verify_tree_dump_format(ts, f)

    def test_num_nodes(self):
        for ts in self.get_example_tree_sequences():
            num_nodes = 0
            for j in range(ts.get_num_records()):
                _, _, node, _, _ = ts.get_record(j)
                if node > num_nodes:
                    num_nodes = node
            self.assertEqual(num_nodes, ts.get_num_nodes())

    def verify_dump_equality(self, ts, outfile):
        """
        Verifies that we can dump a copy of the specified tree sequence
        to the specified file, and load an identical copy.
        """
        ts.dump(outfile.name)
        ts2 = _msprime.TreeSequence()
        ts2.load(outfile.name)
        self.assertEqual(ts.get_sample_size(), ts2.get_sample_size())
        self.assertEqual(ts.get_num_loci(), ts2.get_num_loci())
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

    def test_mutations(self):
        ts = _msprime.TreeSequence()
        # This hasn't been initialised, so should fail.
        self.assertRaises(ValueError, ts.generate_mutations, ts)
        sim = _msprime.Simulator(10, 1)
        sim.run()
        ts.create(sim)
        self.assertRaises(TypeError, ts.generate_mutations)
        self.assertRaises(TypeError, ts.generate_mutations, mutation_rate=1.0)
        self.assertRaises(TypeError, ts.generate_mutations, random_seed=1.0)
        self.assertRaises(
            TypeError, ts.generate_mutations, mutation_rate=10,
            random_seed=1.0, invalid_param=7)
        self.assertIsNone(ts.get_mutation_parameters())
        # A mutation rate of 0 should give 0 mutations
        ts.generate_mutations(0.0, random_seed=1)
        for j in range(3):
            self.assertEqual(ts.get_num_mutations(), 0)
            self.assertEqual(len(ts.get_mutations()), 0)
            self.assertIsNone(ts.get_mutation_parameters())
        # A non-zero mutation rate will give more than 0 mutations.
        ts.generate_mutations(10.0, random_seed=2)
        json_str = ts.get_mutation_parameters()
        self.assertIsNotNone(json_str)
        params = json.loads(json_str)
        self.assertEqual(params["scaled_mutation_rate"], 10.0)
        self.assertEqual(params["random_seed"], 2)
        mutations = ts.get_mutations()
        self.assertGreater(ts.get_num_mutations(), 0)
        self.assertEqual(len(mutations), ts.get_num_mutations())
        # Check the form of the mutations
        for node, position in mutations:
            self.assertIsInstance(node, int)
            self.assertGreater(node, 0)
            self.assertLessEqual(node, ts.get_num_nodes())
            self.assertIsInstance(position, float)
            self.assertGreater(position, 0)
            self.assertLess(position, ts.get_num_loci())
        for j in range(3):
            self.assertEqual(mutations, ts.get_mutations())

    def test_mutation_persistence(self):
        ts = self.get_tree_sequence(mutation_rate=0.0)
        self.assertEqual(ts.get_num_mutations(), 0)
        last_mutations = ts.get_mutations()
        # We should be able to over-write mutations as many times as we like.
        for j in range(10):
            ts.generate_mutations(10, j)
            mutations = ts.get_mutations()
            self.assertNotEqual(mutations, last_mutations)
            last_mutations = mutations

    def test_constructor_interface(self):
        tree_sequence = _msprime.TreeSequence()
        for x in [None, "", {}, [], 1]:
            self.assertRaises(TypeError, tree_sequence.create, x)
        # Creating iterators or running method should fail as we
        # haven't initialised it.
        self.assertRaises(ValueError, tree_sequence.get_record, 0)
        self.assertRaises(ValueError, tree_sequence.get_num_records)
        self.assertRaises(ValueError, _msprime.TreeDiffIterator, tree_sequence)
        sim = _msprime.Simulator(10, 1)
        sim.run()
        tree_sequence.create(sim)
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
        sim = _msprime.Simulator(10, 1)
        sim.run()
        ts.create(sim)
        for bad_type in [None, "", 2.3, [], {}]:
            self.assertRaises(
                TypeError, _msprime.NewickConverter, ts, precision=bad_type)
        before = list(_msprime.NewickConverter(ts))
        self.assertGreater(len(before), 0)
        iterator = _msprime.NewickConverter(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = list(iterator)
        self.assertEqual(before, after)
        # make sure the basic form of the output is correct.
        for length, tree in before:
            self.assertIsInstance(length, int)
            self.assertIsInstance(tree, str)

    def test_iterator(self):
        ts = self.get_tree_sequence()
        ncs = [
            _msprime.NewickConverter(ts), _msprime.NewickConverter(ts, 1),
            _msprime.NewickConverter(ts, 0)]
        for nc in ncs:
            self.verify_iterator(nc)

    def test_precision(self):
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
        ts = self.get_tree_sequence()
        self.assertRaises(ValueError, _msprime.NewickConverter, ts, -1)
        self.assertRaises(ValueError, _msprime.NewickConverter, ts, 17)
        self.assertRaises(ValueError, _msprime.NewickConverter, ts, 100)
        for precision in range(17):
            for l, tree in _msprime.NewickConverter(ts, precision=precision):
                times = get_times(tree)
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
        sim = _msprime.Simulator(10, 1)
        sim.run()
        ts.create(sim)
        before = list(_msprime.TreeDiffIterator(ts))
        iterator = _msprime.TreeDiffIterator(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = list(iterator)
        self.assertEqual(before, after)

    def test_iterator(self):
        ts = self.get_tree_sequence()
        self.verify_iterator(_msprime.TreeDiffIterator(ts))


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
        before = list(hg.get_haplotype(j) for j in range(1, n + 1))
        hg = _msprime.HaplotypeGenerator(ts)
        num_mutations = ts.get_num_mutations()
        del ts
        # We should keep a reference to the tree sequence.
        after = list(hg.get_haplotype(j) for j in range(1, n + 1))
        self.assertEqual(before, after)
        # make sure the basic form of the output is correct.
        for h in before:
            self.assertGreater(len(h), 0)
            self.assertIsInstance(h, str)
            self.assertEqual(len(h), num_mutations)
