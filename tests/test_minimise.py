#
# Copyright (C) 2018 University of Oxford
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
Python implementation of the minimise algorithm.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import unittest

import numpy as np

import msprime


def squash_edges(ts):
    """
    Returns the edges in the tree sequence squashed.
    """
    t = ts.tables.nodes.time
    edges = list(ts.edges())
    edges.sort(key=lambda e: (t[e.parent], e.parent, e.child, e.left))
    if len(edges) == 0:
        return []

    squashed = []
    last_e = edges[0]
    for e in edges[1:]:
        condition = (
            e.parent != last_e.parent or
            e.child != last_e.child or
            e.left != last_e.right)
        if condition:
            squashed.append(last_e)
            last_e = e
        last_e.right = e.right
    squashed.append(last_e)
    return squashed


def minimise(ts):
    """
    Returns a tree sequence with the minimal information required to represent
    the tree topologies at its sites.
    """
    tables = ts.dump_tables()
    edge_map = {}

    def add_edge(left, right, parent, child):
        new_edge = msprime.Edge(left, right, parent, child)
        if child not in edge_map:
            edge_map[child] = new_edge
        else:
            edge = edge_map[child]
            if edge.right == left and edge.parent == parent:
                # Squash
                edge.right = right
            else:
                tables.edges.add_row(edge.left, edge.right, edge.parent, edge.child)
                edge_map[child] = new_edge

    tables.edges.clear()

    edge_buffer = []
    first_site = True
    for tree in ts.trees():
        if tree.num_sites > 0:
            sites = list(tree.sites())
            if first_site:
                x = 0
                first_site = False
            else:
                x = sites[0].position
            # Flush the edge buffer.
            for left, parent, child in edge_buffer:
                add_edge(left, x, parent, child)
            # Add edges for each node in the tree.
            edge_buffer.clear()
            for root in tree.roots:
                for u in tree.nodes(root):
                    if u != root:
                        edge_buffer.append((x, tree.parent(u), u))
    # Add the final edges.
    for left, parent, child in edge_buffer:
        add_edge(left, tables.sequence_length, parent, child)
    # Flush the remaining edges to the table
    for edge in edge_map.values():
        tables.edges.add_row(edge.left, edge.right, edge.parent, edge.child)

    tables.sort()
    # record = msprime.get_provenance_dict(command="minimise")
    # tables.provenances.add_row(record=json.dumps(record))
    return tables.tree_sequence()


class TestMinimise(unittest.TestCase):
    """
    Tests for the minimise function.
    """

    def verify(self, ts):
        source_tables = ts.tables
        positions = set(source_tables.sites.position)
        positions.add(0)
        positions.add(ts.sequence_length)
        mts = minimise(ts)
        for edge in mts.edges():
            self.assertIn(edge.left, positions)
            self.assertIn(edge.right, positions)
        minimised_trees = mts.trees()
        minimised_tree = next(minimised_trees)
        minimised_tree_sites = minimised_tree.sites()
        for tree in ts.trees():
            for site in tree.sites():
                minimised_site = next(minimised_tree_sites, None)
                if minimised_site is None:
                    minimised_tree = next(minimised_trees)
                    minimised_tree_sites = minimised_tree.sites()
                    minimised_site = next(minimised_tree_sites)
                self.assertEqual(site, minimised_site)
            if tree.num_sites > 0:
                self.assertEqual(tree.parent_dict, minimised_tree.parent_dict)
        self.assertTrue(np.array_equal(ts.genotype_matrix(), mts.genotype_matrix()))

        edges = list(mts.edges())
        squashed = squash_edges(mts)
        self.assertEqual(len(edges), len(squashed))
        self.assertEqual(edges, squashed)

    def test_simple_recombination(self):
        ts = msprime.simulate(5, random_seed=1, recombination_rate=2, mutation_rate=2)
        self.verify(ts)

    def test_large_recombination(self):
        ts = msprime.simulate(25, random_seed=12, recombination_rate=5, mutation_rate=15)
        self.verify(ts)

    def test_no_recombination(self):
        ts = msprime.simulate(5, random_seed=1, mutation_rate=2)
        self.verify(ts)

    def test_no_mutation(self):
        ts = msprime.simulate(5, random_seed=1)
        self.verify(ts)
