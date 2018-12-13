"""
Tests for the newick output feature.
"""
from __future__ import print_function
from __future__ import division

import unittest

import msprime

import newick


class TestNewick(unittest.TestCase):
    """
    Tests that the newick output has the properties that we need using
    external Newick parser.
    """
    random_seed = 155

    def verify_newick_topology(self, tree, root=None, node_labels=None):
        if root is None:
            root = tree.root
        ns = tree.newick(precision=16, root=root, node_labels=node_labels)
        if node_labels is None:
            leaf_labels = {u: str(u + 1) for u in tree.leaves(root)}
        else:
            leaf_labels = {u: node_labels[u] for u in tree.leaves(root)}
        newick_tree = newick.loads(ns)[0]
        leaf_names = newick_tree.get_leaf_names()
        self.assertEqual(sorted(leaf_names), sorted(leaf_labels.values()))
        for u in tree.leaves(root):
            name = leaf_labels[u]
            node = newick_tree.get_node(name)
            while u != root:
                self.assertAlmostEqual(node.length, tree.branch_length(u))
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
        ts = msprime.simulate(sample_size=50, random_seed=self.random_seed)
        tables = ts.dump_tables()
        edges = tables.edges
        n = len(edges) // 2
        edges.set_columns(
            left=edges.left[:n], right=edges.right[:n],
            parent=edges.parent[:n], child=edges.child[:n])
        return tables.tree_sequence()

    def test_nonbinary_tree(self):
        ts = self.get_nonbinary_example()
        for t in ts.trees():
            self.verify_newick_topology(t)

    def test_binary_tree(self):
        ts = self.get_binary_example()
        for t in ts.trees():
            self.verify_newick_topology(t)

    def test_multiroot(self):
        ts = self.get_multiroot_example()
        t = ts.first()
        self.assertRaises(ValueError, t.newick)
        for root in t.roots:
            self.verify_newick_topology(t, root=root)

    def test_all_nodes(self):
        ts = msprime.simulate(10, random_seed=5)
        tree = ts.first()
        for u in tree.nodes():
            self.verify_newick_topology(tree, root=u)

    def test_binary_leaf_labels(self):
        tree = self.get_binary_example().first()
        labels = {u: "x_{}".format(u) for u in tree.leaves()}
        self.verify_newick_topology(tree, node_labels=labels)

    def test_nonbinary_leaf_labels(self):
        ts = self.get_nonbinary_example()
        for t in ts.trees():
            labels = {u: str(u) for u in t.leaves()}
            self.verify_newick_topology(t, node_labels=labels)

    def test_all_node_labels(self):
        tree = msprime.simulate(5, random_seed=2).first()
        labels = {u: "x_{}".format(u) for u in tree.nodes()}
        ns = tree.newick(node_labels=labels)
        root = newick.loads(ns)[0]
        self.assertEqual(root.name, labels[tree.root])
        self.assertEqual(
            sorted([n.name for n in root.walk()]), sorted(labels.values()))

    def test_single_node_label(self):
        tree = msprime.simulate(5, random_seed=2).first()
        labels = {tree.root: "XXX"}
        ns = tree.newick(node_labels=labels)
        root = newick.loads(ns)[0]
        self.assertEqual(root.name, labels[tree.root])
        self.assertEqual(
            [n.name for n in root.walk()],
            [labels[tree.root]] + [None for _ in range(len(list(tree.nodes())) - 1)])
