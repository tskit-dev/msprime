"""
Test cases for the low level C interface to msprime.
"""
from __future__ import print_function
from __future__ import division

import random

import tests
import _msprime

class TestInterface(tests.MsprimeTestCase):
    """
    Test the low-level interface to make sure it is robust.
    """
    def verify_random_paramters(self):
        mb = 1024 * 1024
        n = random.randint(2, 1000)
        m = random.randint(1, 10**6)
        r = random.uniform(0, 1000)
        models = []
        random_seed = random.randint(0, 2**31)
        max_memory = random.randint(10 * mb, 100 * mb)
        segment_block_size = random.randint(1, 100)
        node_mapping_block_size = random.randint(1, 100)
        avl_node_block_size = random.randint(1, 100)
        sim = _msprime.Simulator(sample_size=n, num_loci=m,
                population_models=models, scaled_recombination_rate=r,
                random_seed=random_seed, tree_file_name=self._treefile,
                max_memory=max_memory, segment_block_size=segment_block_size,
                avl_node_block_size=avl_node_block_size,
                node_mapping_block_size=node_mapping_block_size)
        # Check the getters to ensure we've got the right values.
        self.assertEqual(n, sim.get_sample_size())
        self.assertEqual(m, sim.get_num_loci())
        self.assertEqual(r, sim.get_scaled_recombination_rate())
        self.assertEqual(random_seed, sim.get_random_seed())
        self.assertEqual(self._treefile, sim.get_tree_file_name())
        # TODO fix population models!
        # self.assertEqual(models, sim.get_population_models())
        self.assertEqual(max_memory, sim.get_max_memory())
        self.assertEqual(segment_block_size, sim.get_segment_block_size())
        self.assertEqual(avl_node_block_size, sim.get_avl_node_block_size())
        self.assertEqual(node_mapping_block_size, sim.get_node_mapping_block_size())
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
        # BUG here!
        # self.assertFalse(sim.run(1e-8))

        # Now run until coalescence
        self.assertTrue(sim.run())

    def test_random_sims(self):
        num_random_sims = 10
        for j in range(num_random_sims):
            self.verify_random_paramters()

    def test_small_sims(self):
        self.verify_simulation(10, 100, 1.0, [])
