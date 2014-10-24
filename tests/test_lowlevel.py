"""
Test cases for the low level C interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import json
import random
import os.path

import tests
import _msprime

class TestInterface(tests.MsprimeTestCase):
    """
    Test the low-level interface to make sure it is robust.
    """
    def verify_running_simulation(self, sim):
        """
        Verifies the state of the specified simulation that has run
        for at least one event.
        """
        self.assertTrue(os.path.exists(sim.get_tree_file_name()))
        self.assertGreater(sim.get_num_trees(), 0)
        self.assertGreater(sim.get_time(), 0.0)
        self.assertGreater(sim.get_num_ancestors(), 1)
        events = sim.get_num_coancestry_events()
        events += sim.get_num_recombination_events()
        self.assertGreater(events, 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        n = sim.get_sample_size()
        m = sim.get_num_loci()
        for ind in sim.get_ancestors():
            for l, r, node in ind:
                self.assertTrue(1 <= l <= m)
                self.assertTrue(1 <= r <= m)
                self.assertTrue(1 <= node <= 2 * n)

    def verify_tree_file_information(self, sim, tf):
        """
        Verifies that the tree file corresponds to the specified
        simulation.
        """
        self.assertEqual(tf.get_sample_size(), sim.get_sample_size())
        self.assertEqual(tf.get_num_loci(), sim.get_num_loci())
        md = json.loads(tf.get_metadata())
        self.assertEqual(sim.get_sample_size(), md["sample_size"])
        self.assertEqual(sim.get_num_loci(), md["num_loci"])
        self.assertEqual(sim.get_random_seed(), md["random_seed"])
        self.assertEqual(sim.get_tree_file_name(), md["tree_file_name"])
        # TODO check rest of metadata.

    def verify_trees(self, sim, sorted_records):
        """
        Verifies that the specified set of sorted coalescence records
        corresponds to correct trees for the specified simulation.
        """
        n = sim.get_sample_size()
        m = sim.get_num_loci()
        pi = [0 for j in range(2 * n)]
        tau = [0.0 for j in range(2 * n)]
        last_l = 1
        for l, c1, c2, parent, t in sorted_records:
            if last_l != l:
                last_l = l
                self.verify_tree(n, pi, tau)
            pi[c1] = parent
            pi[c2] = parent
            tau[parent] = t
        self.verify_tree(n, pi, tau)


    def verify_completed_simulation(self, sim):
        """
        Verifies the state of the specified completed simulation.
        """
        self.assertEqual(sim.get_ancestors(), [])
        self.assertEqual(sim.get_num_ancestors(), 0)
        self.assertTrue(os.path.exists(sim.get_tree_file_name()))
        self.assertGreater(sim.get_num_trees(), 0)
        self.assertGreater(sim.get_time(), 0.0)
        events = sim.get_num_coancestry_events()
        events += sim.get_num_recombination_events()
        self.assertGreater(events, 0)
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        # Open the tree file and get the records
        tf = _msprime.TreeFile(sim.get_tree_file_name())
        self.verify_tree_file_information(sim, tf)
        self.assertTrue(tf.iscomplete())
        self.assertFalse(tf.issorted())
        records = [r for r in tf]
        self.assertGreater(len(records), 0)
        # Records should be in nondecreasing time order
        times = [t for l, c1, c2, parent, t in records]
        self.assertEqual(times, sorted(times))
        self.assertEqual(times[-1], sim.get_time())
        # The iterator has been used, so we should not be able
        # to go through it again.
        r2 = [r for r in tf]
        self.assertEqual(0, len(r2))
        # Sort the records
        _msprime.sort_tree_file(sim.get_tree_file_name())
        # Verify we can't resort the file
        self.assertRaises(ValueError, _msprime.sort_tree_file,
                sim.get_tree_file_name())
        # Read back the records and verify
        tf = _msprime.TreeFile(sim.get_tree_file_name())
        self.verify_tree_file_information(sim, tf)
        self.assertTrue(tf.iscomplete())
        self.assertTrue(tf.issorted())
        sorted_records = [r for r in tf]
        self.assertEqual(sorted_records, sorted(records, key=lambda r: r[0]))
        self.verify_trees(sim, sorted_records)

    def verify_random_paramters(self):
        mb = 1024 * 1024
        n = random.randint(2, 1000)
        m = random.randint(1, 10**6)
        rho = random.uniform(0, 1000)
        models = []
        random_seed = random.randint(0, 2**31)
        max_memory = random.randint(10 * mb, 100 * mb)
        segment_block_size = random.randint(1, 100)
        node_mapping_block_size = random.randint(1, 100)
        avl_node_block_size = random.randint(1, 100)
        sim = _msprime.Simulator(sample_size=n, num_loci=m,
                population_models=models, scaled_recombination_rate=rho,
                random_seed=random_seed, tree_file_name=self._treefile,
                max_memory=max_memory, segment_block_size=segment_block_size,
                avl_node_block_size=avl_node_block_size,
                node_mapping_block_size=node_mapping_block_size)
        # Check initial state
        self.assertEqual(1, sim.get_num_trees())
        self.assertEqual(0.0, sim.get_time())
        self.assertEqual(n, sim.get_num_ancestors())
        self.assertEqual(0, sim.get_num_coancestry_events())
        self.assertEqual(0, sim.get_num_recombination_events())
        self.assertGreater(sim.get_num_avl_node_blocks(), 0)
        self.assertGreater(sim.get_num_segment_blocks(), 0)
        self.assertGreater(sim.get_num_node_mapping_blocks(), 0)
        a = 0
        nodes = set()
        for ind in sim.get_ancestors():
            self.assertEqual(len(ind), 1)
            l, r, node  = ind[0]
            self.assertEqual(l, 1)
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
            self.assertEqual(self._treefile, sim.get_tree_file_name())
            # TODO fix population models!
            # self.assertEqual(models, sim.get_population_models())
            self.assertEqual(max_memory, sim.get_max_memory())
            self.assertEqual(segment_block_size, sim.get_segment_block_size())
            self.assertEqual(avl_node_block_size, sim.get_avl_node_block_size())
            self.assertEqual(node_mapping_block_size, sim.get_node_mapping_block_size())
            # Run this for a tiny amount of time and check the state
            self.assertFalse(sim.run(1e-8))
            self.verify_running_simulation(sim)

    def verify_simulation(self, n, m, r, models):
        """
        Runs the specified simulation and verifies its state.
        """
        mb = 1024 * 1024
        random_seed = random.randint(0, 2**31)
        sim = _msprime.Simulator(sample_size=n, num_loci=m,
                population_models=models, scaled_recombination_rate=r,
                random_seed=random_seed, tree_file_name=self._treefile,
                max_memory=10 * mb, segment_block_size=1000,
                avl_node_block_size=1000, node_mapping_block_size=1000)
        # Run the sim for a tiny amount of time and check.
        self.assertFalse(sim.run(1e-8))
        self.verify_running_simulation(sim)
        # Now run until coalescence
        self.assertTrue(sim.run())
        self.verify_completed_simulation(sim)

    def test_random_sims(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            self.verify_random_paramters()

    def test_small_sims(self):
        self.verify_simulation(3, 1, 0.0, [])
        self.verify_simulation(3, 100, 0.0, [])
        self.verify_simulation(5, 10, 10.0, [])
        self.verify_simulation(10, 100, 1.0, [])

class TestTreeFile(tests.MsprimeTestCase):
    """
    Test cases on the specifics of reading and writing tree files.
    """
    def get_simulator(self, filename):
        sim = _msprime.Simulator(sample_size=10, tree_file_name=filename,
                random_seed=1)
        return sim

    def test_file_errors(self):
        files = ["/no_permissions", "dir/does/not/exist"]
        for f in files:
            self.assertRaises(_msprime.LibraryError, self.get_simulator, f)
            self.assertRaises(_msprime.LibraryError, _msprime.TreeFile, f)
            self.assertRaises(_msprime.LibraryError,
                    _msprime.sort_tree_file, f)

    def test_bad_formats(self):
        # TODO add tests where we take a valid file and mangle it a bit.
        bad_contents = ["",
            "This is not a valid tree file",
        ]
        for contents in bad_contents:
            with open(self._treefile, "w") as f:
                f.write(contents)
            self.assertRaises(_msprime.LibraryError, _msprime.TreeFile,
                    self._treefile)
            self.assertRaises(_msprime.LibraryError, _msprime.sort_tree_file,
                    self._treefile)


