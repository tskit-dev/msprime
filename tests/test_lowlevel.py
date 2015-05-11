"""
Test cases for the low level C interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import collections
import json
import heapq
import random
import os.path
import tempfile

import tests
import _msprime

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

class PythonTreeSequence(object):
    """
    A python implementation of the TreeDiffIterator algorithm.
    """
    def __init__(self, tree_sequence):
        self._tree_sequence = tree_sequence
        self._sample_size = tree_sequence.get_sample_size()
        self._breakpoints = tree_sequence.get_breakpoints()

    def records(self):
        for j in range(self._tree_sequence.get_num_records()):
            yield self._tree_sequence.get_record(j)

    def _diffs(self):
        n = self._sample_size
        left = 0
        used_records = collections.defaultdict(list)
        records_in = []
        for l, r, children, parent, t in self.records():
            if l != left:
                yield l - left, used_records[left], records_in
                del used_records[left]
                records_in = []
                left = l
            used_records[r].append((children, parent, t))
            records_in.append((children, parent, t))
        yield r - left, used_records[left], records_in

    def _diffs_with_breaks(self):
        k = 1
        x = 0
        b = self._breakpoints
        for length, records_out, records_in in self._diffs():
            x += length
            yield b[k] - b[k - 1], records_out, records_in
            while self._breakpoints[k] != x:
                k += 1
                yield b[k] - b[k - 1], [], []
            k += 1

    def diffs(self, all_breaks=False):
        if all_breaks:
            return self._diffs_with_breaks()
        else:
            return self._diffs()



class TestInterface(tests.MsprimeTestCase):
    """
    Test the low-level interface to make sure it is robust.
    """
    def test_gsl_version(self):
        s = _msprime.get_gsl_version()
        self.assertGreater(len(s), 0)
        self.assertIsInstance(s, str)

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
        for l, r, children, p, t in records:
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
        m = sim.get_num_loci()
        pi = {}
        tau = {j:0 for j in range(1, n + 1)}
        last_l = 0
        last_t = 0
        num_trees = 0
        live_segments = []
        for l, r, children, parent, t in sorted_records:
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
            heapq.heappush(live_segments, (r, (children, parent)))
            while live_segments[0][0] <= l:
                x, (other_children, p) = heapq.heappop(live_segments)
                for c in other_children:
                    del pi[c]
                del tau[p]
            pi[children[0]] = parent
            pi[children[1]] = parent
            tau[parent] = t
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
        times = [t for l, r, children, parent, t in records]
        # Children should always be sorted in order.
        for _, _, children, _, _ in records:
            self.assertEqual(children, tuple(sorted(children)))
        self.assertEqual(times, sorted(times))
        self.assertEqual(times[-1], sim.get_time())
        self.verify_squashed_records(records)
        sorted_records = sorted(records, key=lambda r: r[0])
        self.verify_trees(sim, sorted_records)
        # Check the TreeSequence
        ts = _msprime.TreeSequence()
        ts.create(sim)
        records = [ts.get_record(j) for j in range(len(records))]
        self.assertEqual(records, sorted_records)


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
        sim = _msprime.Simulator(sample_size=n, num_loci=m,
                population_models=models,
                scaled_recombination_rate=rho,
                random_seed=random_seed, max_memory=max_memory,
                segment_block_size=segment_block_size,
                avl_node_block_size=avl_node_block_size,
                node_mapping_block_size=node_mapping_block_size,
                coalescence_record_block_size=coalescence_record_block_size)
        # Check initial state
        self.assertEqual(2, sim.get_num_breakpoints())
        self.assertEqual(0.0, sim.get_time())
        self.assertEqual(n, sim.get_num_ancestors())
        self.assertEqual(0, sim.get_num_coancestry_events())
        self.assertEqual(0, sim.get_num_recombination_events())
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        self.assertGreater(sim.get_num_coalescence_record_blocks(), 0)
        self.assertEqual(sim.get_sample_size(), n)
        self.assertEqual(sim.get_num_loci(), m)
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
            self.assertEqual(avl_node_block_size, sim.get_avl_node_block_size())
            self.assertEqual(node_mapping_block_size, sim.get_node_mapping_block_size())
            self.assertEqual(coalescence_record_block_size,
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
        # TODO fix this API!
        pop_models = sim.get_population_models()
        # Exclude first and last elements.
        self.assertGreater(len(pop_models), 1)
        pop_models = pop_models[1:-1]
        self.assertEqual(pop_models, models)

    def verify_iterators(self, tree_sequence):
        """
        Verifies the iterator protocol is properly implemented for the
        low-level iterator classes.
        """
        iters = [_msprime.TreeDiffIterator(tree_sequence)]
        for iterator in iters:
            l = list(iterator)
            self.assertGreater(len(l), 0)
            for j in range(10):
                self.assertRaises(StopIteration, next, iterator)

    def verify_tree_diffs(self, tree_sequence):
        n = tree_sequence.get_sample_size()
        m = tree_sequence.get_num_loci()
        for all_breakpoints in [False, True]:
            # Check some basic properties of the diffs.
            diffs = list(_msprime.TreeDiffIterator(
                tree_sequence, all_breakpoints))
            length, records_out, records_in = diffs[0]
            self.assertGreaterEqual(length, 0)
            self.assertEqual(len(records_out), 0)
            self.assertEqual(len(records_in), n - 1)
            self.assertEqual(sum([l for l, _, _ in diffs]), m)
            for l, records_out, records_in in diffs[1:]:
                self.assertGreaterEqual(l, 0)
                self.assertGreaterEqual(l, 0)
                self.assertEqual(len(records_out), len(records_in))
                for children, parent, time in records_out + records_in:
                    for c in children:
                        self.assertGreater(c, 0)
                        self.assertGreater(parent, c)
                        self.assertGreater(time, 0.0)
            # Compare with the Python implementation.
            pts = PythonTreeSequence(tree_sequence)
            python_diffs = list(pts.diffs(all_breakpoints))
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
        sim = _msprime.Simulator(sample_size=n, num_loci=m,
                population_models=models, scaled_recombination_rate=r,
                random_seed=random_seed,
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
        self.verify_iterators(tree_sequence)
        self.verify_tree_diffs(tree_sequence)
        self.assertEqual(tree_sequence.get_breakpoints(), sim.get_breakpoints())

    def test_tree_sequence_interface(self):
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
        # We don't accept Python negative indexes here.
        self.assertRaises(IndexError, tree_sequence.get_record, -1)
        for j in [0, 10, 10**6]:
            self.assertRaises(IndexError, tree_sequence.get_record,
                    num_records + j)
        # Check out the diff iterator.
        iterator = _msprime.TreeDiffIterator(tree_sequence)
        records = list(iterator)
        self.assertGreater(len(records), 0)
        # Now we do something horrible. Allocate an iterator, but then delete
        # the underlying TreeSequence. We should still be able to iterate over
        # the records, as we should keep a reference to the TreeSequence in
        # the iterator, preventing it from being freed.
        iterator = _msprime.TreeDiffIterator(tree_sequence)
        del tree_sequence
        new_records = list(iterator)
        self.assertEqual(new_records, records)

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
        const_model = _msprime.POP_MODEL_CONSTANT

    def test_population_models(self):
        exp_model = _msprime.POP_MODEL_EXPONENTIAL
        m1 = {"alpha":0.0, "start_time":0.0, "type":exp_model}
        m2 = {"alpha":1.0, "start_time":0.5, "type":exp_model}
        # TODO add constant pop models here too.
        for n in [5, 10, 100]:
            self.verify_simulation(n, 1, 0.0, [m1])
            self.verify_simulation(n, 1, 0.0, [m2])
            self.verify_simulation(n, 1, 0.0, [m1, m2])

    def test_bad_parameters(self):
        def f(n=2, m=1, **kwargs):
            return _msprime.Simulator(n, m, **kwargs)
        self.assertRaises(TypeError, f, n=None)
        self.assertRaises(TypeError, f, m=None)
        self.assertRaises(OverflowError, f, max_memory=2**65)
        # TODO add more tests!

    def test_bad_population_models(self):
        def f(population_models):
            return _msprime.Simulator(2, 1,
                    population_models=population_models)
        self.assertRaises(TypeError, f, "")
        self.assertRaises(TypeError, f, [""])
        self.assertRaises(_msprime.InputError, f, [{}])
        self.assertRaises(_msprime.InputError, f, [{"start_time":1}])
        self.assertRaises(_msprime.InputError, f,
                [{"start_time":1, "type":-1000}])
        self.assertRaises(_msprime.InputError, f,
                [{"start_time":1, "type":_msprime.POP_MODEL_EXPONENTIAL}])
        self.assertRaises(_msprime.InputError, f,
                [{"start_time":1, "type":_msprime.POP_MODEL_CONSTANT}])
        m = get_random_population_models(10)
        m.reverse()
        # TODO This should really be an input error.
        self.assertRaises(_msprime.LibraryError, f, m)

