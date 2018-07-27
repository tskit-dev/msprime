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
Tests for the minimise algorithm which reduces a tree sequence down to
the topology visible at its sites.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import unittest

import numpy as np

import msprime
from tests.simplify import Simplifier
import tests.tsutil as tsutil


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
    the tree topologies at its sites. Uses a left-to-right algorithm.
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
        # print(tree.interval)
        # print(tree.draw(format="unicode"))
        if tree.num_sites > 0:
            sites = list(tree.sites())
            if first_site:
                x = 0
                # print("First site", sites)
                first_site = False
            else:
                x = sites[0].position
            # Flush the edge buffer.
            for left, parent, child in edge_buffer:
                add_edge(left, x, parent, child)
            # Add edges for each node in the tree.
            edge_buffer = []
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
    ts = tables.tree_sequence()
    # Now simplify to remove redundant nodes.
    return ts.simplify(map_nodes=True, filter_zero_mutation_sites=False)


class TestMinimise(unittest.TestCase):
    """
    Tests for the minimise function.
    """

    def verify(self, ts):
        source_tables = ts.tables
        X = source_tables.sites.position
        position_count = {x: 0 for x in X}
        position_count[0] = 0
        position_count[ts.sequence_length] = 0
        mts, node_map = minimise(ts)
        for edge in mts.edges():
            self.assertIn(edge.left, position_count)
            self.assertIn(edge.right, position_count)
            position_count[edge.left] += 1
            position_count[edge.right] += 1
        if ts.num_sites == 0:
            # We should have zero edges output.
            self.assertEqual(mts.num_edges, 0)
        elif X[0] != 0:
            # The first site (if it's not zero) should be mapped to zero so
            # this never occurs in edges.
            self.assertEqual(position_count[X[0]], 0)

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
                self.assertEqual(site.position, minimised_site.position)
                self.assertEqual(site.ancestral_state, minimised_site.ancestral_state)
                self.assertEqual(site.metadata, minimised_site.metadata)
                self.assertEqual(len(site.mutations), len(minimised_site.mutations))

                for mutation, minimised_mutation in zip(
                        site.mutations, minimised_site.mutations):
                    self.assertEqual(
                        mutation.derived_state, minimised_mutation.derived_state)
                    self.assertEqual(mutation.metadata, minimised_mutation.metadata)
                    self.assertEqual(mutation.parent, minimised_mutation.parent)
                    self.assertEqual(node_map[mutation.node], minimised_mutation.node)
            if tree.num_sites > 0:
                mapped_dict = {
                    node_map[u]: node_map[v] for u, v in tree.parent_dict.items()}
                self.assertEqual(mapped_dict, minimised_tree.parent_dict)
        self.assertTrue(np.array_equal(ts.genotype_matrix(), mts.genotype_matrix()))

        edges = list(mts.edges())
        squashed = squash_edges(mts)
        self.assertEqual(len(edges), len(squashed))
        self.assertEqual(edges, squashed)

        # Verify against simplify implementations.
        s = Simplifier(
            ts, ts.samples(), reduce_to_site_topology=True,
            filter_zero_mutation_sites=False)
        sts, _ = s.simplify()

        t1 = mts.tables
        t2 = sts.tables
        self.assertEqual(t1.nodes,  t2.nodes)
        self.assertEqual(t1.edges, t2.edges)
        self.assertEqual(t1.sites, t2.sites)
        self.assertEqual(t1.mutations, t2.mutations)
        self.assertEqual(t1.populations, t2.populations)
        self.assertEqual(t1.individuals, t2.individuals)
        return mts

    def test_no_recombination_one_site(self):
        ts = msprime.simulate(15, random_seed=1)
        tables = ts.dump_tables()
        tables.sites.add_row(position=0.25, ancestral_state="0")
        mts = self.verify(tables.tree_sequence())
        self.assertEqual(mts.num_trees, 1)

    def test_simple_recombination_one_site(self):
        ts = msprime.simulate(15, random_seed=1, recombination_rate=2)
        tables = ts.dump_tables()
        tables.sites.add_row(position=0.25, ancestral_state="0")
        mts = self.verify(tables.tree_sequence())
        self.assertEqual(mts.num_trees, 1)

    def test_simple_recombination_fixed_sites(self):
        ts = msprime.simulate(5, random_seed=1, recombination_rate=2)
        tables = ts.dump_tables()
        for x in [0.25, 0.5, 0.75]:
            tables.sites.add_row(position=x, ancestral_state="0")
        self.verify(tables.tree_sequence())

    def get_integer_edge_ts(self, n, m):
        recombination_map = msprime.RecombinationMap.uniform_map(m, 1, num_loci=m)
        ts = msprime.simulate(n, random_seed=1, recombination_map=recombination_map)
        self.assertGreater(ts.num_trees, 1)
        for edge in ts.edges():
            self.assertEqual(int(edge.left), edge.left)
            self.assertEqual(int(edge.right), edge.right)
        return ts

    def test_integer_edges_one_site(self):
        ts = self.get_integer_edge_ts(5, 10)
        tables = ts.dump_tables()
        tables.sites.add_row(position=1, ancestral_state="0")
        mts = self.verify(tables.tree_sequence())
        self.assertEqual(mts.num_trees, 1)

    def test_integer_edges_all_sites(self):
        ts = self.get_integer_edge_ts(5, 10)
        tables = ts.dump_tables()
        for x in range(10):
            tables.sites.add_row(position=x, ancestral_state="0")
        mts = self.verify(tables.tree_sequence())
        self.assertEqual(mts.num_trees, ts.num_trees)

    def test_simple_recombination_site_at_zero(self):
        ts = msprime.simulate(5, random_seed=1, recombination_rate=2)
        tables = ts.dump_tables()
        tables.sites.add_row(position=0, ancestral_state="0")
        mts = self.verify(tables.tree_sequence())
        self.assertEqual(mts.num_trees, 1)

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

    def test_many_roots(self):
        ts = msprime.simulate(25, random_seed=12, recombination_rate=2, length=10)
        tables = tsutil.decapitate(ts, ts.num_edges // 2).dump_tables()
        for x in range(10):
            tables.sites.add_row(x, "0")
        self.verify(tables.tree_sequence())

    def test_branch_sites(self):
        ts = msprime.simulate(15, random_seed=12, recombination_rate=2, length=10)
        ts = tsutil.insert_branch_sites(ts)
        self.verify(ts)

    def test_jiggled_samples(self):
        ts = msprime.simulate(8, random_seed=13, recombination_rate=2, length=10)
        ts = tsutil.jiggle_samples(ts)
        self.verify(ts)
