#
# Copyright (C) 2017 University of Oxford
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
Tests for the newick output feature.
"""
from __future__ import print_function
from __future__ import division

import unittest

import msprime

import newick


class NewickTest(unittest.TestCase):
    """
    Tests that the newick output has the properties that we need using
    external Newick parser.
    """
    random_seed = 155

    def verify_newick_topology(self, tree, time_scale=1):
        ns = tree.newick(precision=16, time_scale=time_scale)
        newick_tree = newick.loads(ns)[0]
        leaf_names = newick_tree.get_leaf_names()
        self.assertEqual(
            sorted(leaf_names),
            sorted([str(u + 1) for u in tree.leaves()]))
        for u in tree.leaves():
            name = str(u + 1)
            node = newick_tree.get_node(name)
            while tree.parent(u) != msprime.NULL_NODE:
                self.assertAlmostEqual(node.length, time_scale * tree.branch_length(u))
                node = node.ancestor
                u = tree.parent(u)
            self.assertIsNone(node.ancestor)

    def get_nonbinary_example(self):
        ts = msprime.simulate(
            sample_size=20, recombination_rate=10, random_seed=self.random_seed,
            demographic_events=[
                msprime.SimpleBottleneck(time=0.5, population=0, proportion=1)])
        # Make sure this really has some non-binary nodes
        found = False
        for e in ts.edgesets():
            if len(e.children) > 2:
                found = True
                break
        self.assertTrue(found)
        return ts

    def get_binary_example(self):
        ts = msprime.simulate(
            sample_size=25, recombination_rate=5, random_seed=self.random_seed)
        return ts

    def get_multiroot_example(self):
        ts = msprime.simulate(
            sample_size=50, recombination_rate=5, random_seed=self.random_seed)
        tables = ts.dump_tables()
        edges = tables.edges
        n = len(edges) // 2
        edges.set_columns(
            left=edges.left[:n], right=edges.right[:n],
            parent=edges.parent[:n], child=edges.child[:n])
        return msprime.load_tables(nodes=tables.nodes, edges=edges)

    def test_nonbinary_tree(self):
        ts = self.get_nonbinary_example()
        for t in ts.trees():
            self.verify_newick_topology(t)

    def test_binary_tree(self):
        ts = self.get_binary_example()
        for t in ts.trees():
            self.verify_newick_topology(t)

    def test_time_scale(self):
        tree = next(self.get_binary_example().trees())
        for scale in [0, -1, 0.5, 1e6]:
            self.verify_newick_topology(tree, scale)
