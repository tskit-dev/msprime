"""
Test cases for branch length statistic computation.
"""
from __future__ import print_function
from __future__ import division


import unittest
import math
import random
import itertools

import six

import msprime

from simplify_algorithms import *


class SimplifyTestCase(unittest.TestCase):
    """
    Tests of the simplify implementation.
    """
    random_seed = 23


def single_childify(ts):
    """
    Builds a new equivalent tree sequence whose edgesets all have singleton children.
    """
    tables = ts.dump_tables()
    print(tables.edgesets)
    tables.edgesets.reset()
    for u in range(ts.num_nodes):
        parent_edges = [r for r in ts.edgesets() if u == r.parent]
        print(u, "__", parent_edges)
        children = []
        for r in parent_edges:
            children.extend(r.children)
        for child in set(children):
            edges = [r for r in parent_edges if child in r.children]
            lefts = [r.left for r in edges]
            rights = [r.right for r in edges]
            do_lefts = [lefts[0]]
            do_rights = []
            for k in range(len(lefts)-1):
                if lefts[k+1] != rights[k]:
                    do_lefts.append(lefts[k+1])
                    do_rights.append(rights[k])
            do_rights.append(rights[-1])
            print(child, " :: ", do_lefts, " ---- ", do_rights)
            assert len(do_lefts) == len(do_rights)
            for k in range(len(do_lefts)):
                tables.edgesets.add_row(
                    left=do_lefts[k], right=do_rights[k], children=(child,), parent=u)
    print(tables.edgesets)
    new_ts = msprime.load_tables(**tables._asdict())
    return new_ts


def do_simplify(ts, sample=None):
    if sample is None:
        sample = ts.samples()
    s = Simplifier(ts, sample)
    new_ts = s.simplify()
    return new_ts


def print_tables(ts):
    nodes = msprime.NodeTable()
    edges = msprime.EdgesetTable()
    ts.dump_tables(nodes=nodes, edgesets=edges)
    print(nodes)
    print(edges)


def reset_population(ts):
    """
    Set population to zero.
    """
    tables = ts.dump_tables()
    tables.nodes.reset()
    for x in ts.nodes():
        tables.nodes.add_row(flags=x.flags, population=0, time=x.time)
    new_ts = msprime.load_tables(**tables._asdict())
    return new_ts


class TestOverlappingRecords(SimplifyTestCase):
    """
    Tests that simplify() can deal with overlapping records.
    """

    def test_single_tree(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        ts_single = single_childify(ts)
        tss_nopops = do_simplify(ts_single)
        # simplify is not currently recording population: set to 0
        tss = reset_population(tss_nopops)
        self.assertEqual(list(tss.records()), list(ts.records()))

    def test_many_trees(self):
        ts = msprime.simulate(
                5, recombination_rate=5,
                random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        ts_single = single_childify(ts)
        tss_nopops = do_simplify(ts_single)
        # simplify is not currently recording population: set to 0
        tss = reset_population(tss_nopops)
        print("ts:")
        print_tables(ts)
        print("ts_single:")
        print_tables(ts_single)
        print("tss:")
        print_tables(tss)
        self.assertEqual(list(tss.records()), list(ts.records()))


class TestWithVisuals(SimplifyTestCase):
    """
    Some pedantic tests with ascii depictions of what's supposed to happen.
    """

    def verify_simplify_topology(self, ts, sample):
        # copies from test_topology.py
        print("## sample:")
        print(sample)
        new_ts = do_simplify(ts, sample)
        print("## ts:")
        print_tables(ts)
        print("## simplified:")
        print_tables(new_ts)
        sample_map = {k: j for j, k in enumerate(sample)}
        old_trees = ts.trees()
        old_tree = next(old_trees)
        self.assertGreaterEqual(ts.get_num_trees(), new_ts.get_num_trees())
        for new_tree in new_ts.trees():
            new_left, new_right = new_tree.get_interval()
            old_left, old_right = old_tree.get_interval()
            # Skip ahead on the old tree until new_left is within its interval
            while old_right <= new_left:
                old_tree = next(old_trees)
                old_left, old_right = old_tree.get_interval()
            # If the TMRCA of all pairs of samples is the same, then we have the
            # same information. We limit this to at most 500 pairs
            pairs = itertools.islice(itertools.combinations(sample, 2), 500)
            for pair in pairs:
                mapped_pair = [sample_map[u] for u in pair]
                mrca1 = old_tree.get_mrca(*pair)
                mrca2 = new_tree.get_mrca(*mapped_pair)
                self.assertEqual(old_tree.get_time(mrca1), new_tree.get_time(mrca2))
                self.assertEqual(
                    old_tree.get_population(mrca1), new_tree.get_population(mrca2))

    def test_partial_non_sample_external_nodes(self):
        # A somewhat more complicated test case with a partially specified,
        # non-sampled tip.
        #
        # Here is the situation:
        #
        # 1.0             7
        # 0.7            / \                                            6
        #               /   \                                          / \
        # 0.5          /     5                      5                 /   5
        #             /     / \                    / \               /   / \
        # 0.4        /     /   4                  /   4             /   /   4
        #           /     /   / \                /   / \           /   /   / \
        #          /     /   3   \              /   /   \         /   /   3   \
        #         /     /         \            /   /     \       /   /         \
        # 0.0    0     1           2          1   0       2     0   1           2
        #
        #          (0.0, 0.2),                 (0.2, 0.8),         (0.8, 1.0)

        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           0.2  # Non sample leaf
        4       0           0.4
        5       0           0.5
        6       0           0.7
        7       0           1.0
        """)
        edgesets = six.StringIO("""\
        left    right   parent  children
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.0     0.2     7       0,5
        """)
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6: -1, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: -1, 4: 5, 5: -1, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        ts = msprime.load_text(nodes=nodes, edgesets=edgesets)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(ts.num_nodes, 8)
        # check topologies agree:
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        # check .simplify() works here
        self.verify_simplify_topology(ts, [0, 1, 2])

    def test_partial_non_sample_external_nodes_2(self):
        # The same situation as above, but partial tip is labeled '7' not '3':
        #
        # 1.0          6
        # 0.7         / \                                       5
        #            /   \                                     / \
        # 0.5       /     4                 4                 /   4
        #          /     / \               / \               /   / \
        # 0.4     /     /   3             /   3             /   /   3
        #        /     /   / \           /   / \           /   /   / \
        #       /     /   7   \         /   /   \         /   /   7   \
        #      /     /         \       /   /     \       /   /         \
        # 0.0 0     1           2     1   0       2     0   1           2
        #
        #          (0.0, 0.2),         (0.2, 0.8),         (0.8, 1.0)
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           0.4
        4       0           0.5
        5       0           0.7
        6       0           1.0
        7       0           0    # Non sample leaf
        """)
        edgesets = six.StringIO("""\
        left    right   parent  children
        0.0     0.2     3       2,7
        0.2     0.8     3       0,2
        0.8     1.0     3       2,7
        0.0     0.2     4       1,3
        0.2     0.8     4       1,3
        0.8     1.0     4       1,3
        0.8     1.0     5       0,4
        0.0     0.2     6       0,4
        """)
        true_trees = [
            {0: 6, 1: 4, 2: 3, 3: 4, 4: 6, 5: -1, 6: -1, 7: 3},
            {0: 3, 1: 4, 2: 3, 3: 4, 4: -1, 5: -1, 6: -1, 7: -1},
            {0: 5, 1: 4, 2: 3, 3: 4, 4: 5, 5: -1, 6: -1, 7: 3}]
        ts = msprime.load_text(nodes=nodes, edgesets=edgesets)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        # sample size check works here since 7 > 3
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(ts.num_nodes, 8)
        # check topologies agree:
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        self.verify_simplify_topology(ts, [0, 1, 2])

    def test_single_offspring_records(self):
        # Here we have inserted a single-offspring record
        # (for 6 on the left segment):
        #
        # 1.0             7
        # 0.7            / 6                                                  6
        #               /   \                                                / \
        # 0.5          /     5                       5                      /   5
        #             /     / \                     / \                    /   / \
        # 0.4        /     /   4                   /   4                  /   /   4
        # 0.3       /     /   / \                 /   / \                /   /   / \
        #          /     /   3   \               /   /   \              /   /   3   \
        #         /     /         \             /   /     \            /   /         \
        # 0.0    0     1           2           1   0       2          0   1           2
        #
        #          (0.0, 0.2),               (0.2, 0.8),              (0.8, 1.0)
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   1           0
        3   0           0       # Non sample leaf
        4   0           0.4
        5   0           0.5
        6   0           0.7
        7   0           1.0
        """)
        edgesets = six.StringIO("""\
        left    right   parent  children
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.0     0.2     6       5
        0.0     0.2     7       0,6
        """)
        ts = msprime.load_text(nodes, edgesets)
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: 7, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: -1, 4: 5, 5: -1, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(ts.num_nodes, 8)
        # check topologies agree:
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        self.verify_simplify_topology(ts, [0, 1, 2])

    def test_many_single_offspring(self):
        # a more complex test with single offspring
        # With `(i,j,x)->k` denoting that individual `k` inherits from `i` on `[0,x)`
        #    and from `j` on `[x,1)`:
        # 1. Begin with an individual `3` (and another anonymous one) at `t=0`.
        # 2. `(3,?,1.0)->4` and `(3,?,1.0)->5` at `t=1`
        # 3. `(4,3,0.9)->6` and `(3,5,0.1)->7` and then `3` dies at `t=2`
        # 4. `(6,7,0.7)->8` at `t=3`
        # 5. `(8,6,0.8)->9` and `(7,8,0.2)->10` at `t=4`.
        # 6. `(3,9,0.6)->0` and `(9,10,0.5)->1` and `(10,4,0.4)->2` at `t=5`.
        # 7. We sample `0`, `1`, and `2`.
        # Here are the trees:
        # t                  |              |              |             |
        #
        # 0       --3--      |     --3--    |     --3--    |    --3--    |    --3--
        #        /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \
        # 1     4   |   5    |   4   |   5  |   4       5  |  4       5  |  4       5
        #       |\ / \ /|    |   |\   \     |   |\     /   |  |\     /   |  |\     /|
        # 2     | 6   7 |    |   | 6   7    |   | 6   7    |  | 6   7    |  | 6   7 |
        #       | |\ /| |    |   |  \  |    |   |  \  |    |  |  \       |  |  \    | ...
        # 3     | | 8 | |    |   |   8 |    |   |   8 |    |  |   8      |  |   8   |
        #       | |/ \| |    |   |  /  |    |   |  /  |    |  |  / \     |  |  / \  |
        # 4     | 9  10 |    |   | 9  10    |   | 9  10    |  | 9  10    |  | 9  10 |
        #       |/ \ / \|    |   |  \   \   |   |  \   \   |  |  \   \   |  |  \    |
        # 5     0   1   2    |   0   1   2  |   0   1   2  |  0   1   2  |  0   1   2
        #
        #                    |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5
        # ... continued:
        # t                  |             |             |             |
        #
        # 0         --3--    |    --3--    |    --3--    |    --3--    |    --3--
        #          /     \   |   /     \   |   /     \   |   /     \   |   /  |  \
        # 1       4       5  |  4       5  |  4       5  |  4       5  |  4   |   5
        #         |\     /|  |   \     /|  |   \     /|  |   \     /|  |     /   /|
        # 2       | 6   7 |  |    6   7 |  |    6   7 |  |    6   7 |  |    6   7 |
        #         |  \    |  |     \    |  |       /  |  |    |  /  |  |    |  /  |
        # 3  ...  |   8   |  |      8   |  |      8   |  |    | 8   |  |    | 8   |
        #         |  / \  |  |     / \  |  |     / \  |  |    |  \  |  |    |  \  |
        # 4       | 9  10 |  |    9  10 |  |    9  10 |  |    9  10 |  |    9  10 |
        #         |    /  |  |   /   /  |  |   /   /  |  |   /   /  |  |   /   /  |
        # 5       0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2
        #
        #         0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0

        true_trees = [
            {0: 4, 1: 9, 2: 10, 3: -1, 4: 3, 5: 3, 6: 4, 7: 3, 8: 6, 9: 8, 10: 7},
            {0: 4, 1: 9, 2: 10, 3: -1, 4: 3, 5: 3, 6: 4, 7: 5, 8: 6, 9: 8, 10: 7},
            {0: 4, 1: 9, 2: 10, 3: -1, 4: 3, 5: 3, 6: 4, 7: 5, 8: 6, 9: 8, 10: 8},
            {0: 4, 1: 9,  2: 5, 3: -1, 4: 3, 5: 3, 6: 4, 7: 5, 8: 6, 9: 8, 10: 8},
            {0: 4, 1: 10, 2: 5, 3: -1, 4: 3, 5: 3, 6: 4, 7: 5, 8: 6, 9: 8, 10: 8},
            {0: 9, 1: 10, 2: 5, 3: -1, 4: 3, 5: 3, 6: 4, 7: 5, 8: 6, 9: 8, 10: 8},
            {0: 9, 1: 10, 2: 5, 3: -1, 4: 3, 5: 3, 6: 4, 7: 5, 8: 7, 9: 8, 10: 8},
            {0: 9, 1: 10, 2: 5, 3: -1, 4: 3, 5: 3, 6: 4, 7: 5, 8: 7, 9: 6, 10: 8},
            {0: 9, 1: 10, 2: 5, 3: -1, 4: 3, 5: 3, 6: 3, 7: 5, 8: 7, 9: 6, 10: 8}
        ]
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           5
        4       0           4
        5       0           4
        6       0           3
        7       0           3
        8       0           2
        9       0           1
        10      0           1
        """)
        edgesets = six.StringIO("""\
        left    right   parent  children
        0.5     1.0     10      1
        0.0     0.4     10      2
        0.6     1.0     9       0
        0.0     0.5     9       1
        0.8     1.0     8       10
        0.2     0.8     8       9,10
        0.0     0.2     8       9
        0.7     1.0     7       8
        0.0     0.2     7       10
        0.8     1.0     6       9
        0.0     0.7     6       8
        0.4     1.0     5       2,7
        0.1     0.4     5       7
        0.6     0.9     4       6
        0.0     0.6     4       0,6
        0.9     1.0     3       4,5,6
        0.1     0.9     3       4,5
        0.0     0.1     3       4,5,7
        """)

        ts = msprime.load_text(nodes, edgesets)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, len(true_trees))
        self.assertEqual(ts.num_nodes, 11)
        self.assertEqual(len(list(ts.diffs())), ts.num_trees)
        # check topologies agree:
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        self.verify_simplify_topology(ts, [0, 1])
        self.verify_simplify_topology(ts, [1, 2])
        self.verify_simplify_topology(ts, [2, 0])

    def test_tricky_switches(self):
        # suppose the topology has:
        # left right parent children
        #  0.0   0.5      6      0,1
        #  0.5   1.0      6      4,5
        #  0.0   0.4      7      2,3
        #
        # --------------------------
        #
        #        12         .        12         .        12         .
        #       /  \        .       /  \        .       /  \        .
        #     11    \       .      /    \       .      /    \       .
        #     / \    \      .     /     10      .     /     10      .
        #    /   \    \     .    /     /  \     .    /     /  \     .
        #   6     7    8    .   6     9    8    .   6     9    8    .
        #  / \   / \   /\   .  / \   / \   /\   .  / \   / \   /\   .
        # 0   1 2   3 4  5  . 0   1 2   3 4  5  . 4   5 2   3 0  1  .
        #                   .                   .                   .
        # 0.0              0.4                 0.5                 1.0
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       1           0
        4       1           0
        5       1           0
        6       0           1
        7       0           1
        8       0           1
        9       0           1
        10      0           2
        11      0           3
        12      0           4
        """)
        edgesets = six.StringIO("""\
        left right parent children
        0.0  0.5   6      0,1
        0.5  1.0   6      4,5
        0.0  0.4   7      2,3
        0.0  0.5   8      4,5
        0.5  1.0   8      0,1
        0.4  1.0   9      2,3
        0.4  1.0   10     8,9
        0.0  0.4   11     6,7
        0.0  0.4   12     8,11
        0.4  1.0   12     6,10
        """)
        true_trees = [
                {0: 6, 1: 6, 2: 7, 3: 7, 4: 8, 5: 8, 6: 11,
                    7: 11, 8: 12, 9: -1, 10: -1, 11: 12, 12: -1},
                {0: 6, 1: 6, 2: 9, 3: 9, 4: 8, 5: 8, 6: 12,
                    7: -1, 8: 10, 9: 10, 10: 12, 11: -1, 12: -1},
                {0: 8, 1: 8, 2: 9, 3: 9, 4: 6, 5: 6, 6: 12,
                    7: -1, 8: 10, 9: 10, 10: 12, 11: -1, 12: -1}
        ]
        ts = msprime.load_text(nodes, edgesets)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 6)
        self.assertEqual(ts.num_trees, len(true_trees))
        self.assertEqual(ts.num_nodes, 13)
        self.assertEqual(len(list(ts.diffs())), ts.num_trees)
        # check topologies agree:
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        for samples in [[0, 2], [0, 4], [2, 4]]:
            print("Verifying", samples)
            self.verify_simplify_topology(ts, samples)

    def test_tricky_simplify(self):
        # Continue as above but invoke simplfy:
        #
        #         12         .          12         .
        #        /  \        .         /  \        .
        #      11    \       .       11    \       .
        #      / \    \      .       / \    \      .
        #    13   \    \     .      /  15    \     .
        #    / \   \    \    .     /   / \    \    .
        #   6  14   7    8   .    6  14   7    8   .
        #  / \     / \   /\  .   / \     / \   /\  .
        # 0   1   2   3 4  5 .  0   1   2   3 4  5 .
        #                    .                     .
        # 0.0               0.1                   0.4
        #
        #  .        12         .        12         .
        #  .       /  \        .       /  \        .
        #  .      /    \       .      /    \       .
        #  .     /     10      .     /     10      .
        #  .    /     /  \     .    /     /  \     .
        #  .   6     9    8    .   6     9    8    .
        #  .  / \   / \   /\   .  / \   / \   /\   .
        #  . 0   1 2   3 4  5  . 4   5 2   3 0  1  .
        #  .                   .                   .
        # 0.4                 0.5                 1.0
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       1           0
        4       1           0
        5       1           0
        6       0           1
        7       0           1
        8       0           1
        9       0           1
        10      0           2
        11      0           3
        12      0           4
        13      0           2
        14      0           1
        15      0           2
        """)
        edgesets = six.StringIO("""\
        left right parent children
        0.0  0.5   6      0,1
        0.5  1.0   6      4,5
        0.0  0.4   7      2,3
        0.0  0.5   8      4,5
        0.5  1.0   8      0,1
        0.4  1.0   9      2,3
        0.4  1.0   10     8,9
        0.0  0.1   13     6,14
        0.1  0.4   15     7,14
        0.0  0.1   11     7,13
        0.1  0.4   11     6,15
        0.0  0.4   12     8,11
        0.4  1.0   12     6,10
        """)
        true_trees = [
                {0: 6, 1: 6, 2: 7, 3: 7, 4: 8, 5: 8, 6: 11,
                    7: 11, 8: 12, 9: -1, 10: -1, 11: 12, 12: -1},
                {0: 6, 1: 6, 2: 9, 3: 9, 4: 8, 5: 8, 6: 12,
                    7: -1, 8: 10, 9: 10, 10: 12, 11: -1, 12: -1},
                {0: 8, 1: 8, 2: 9, 3: 9, 4: 6, 5: 6, 6: 12,
                    7: -1, 8: 10, 9: 10, 10: 12, 11: -1, 12: -1}
        ]
        big_ts = msprime.load_text(nodes, edgesets)
        self.assertEqual(big_ts.num_trees, 1+len(true_trees))
        self.assertEqual(big_ts.num_nodes, 16)
        ts = big_ts.simplify()
        self.assertEqual(ts.sample_size, 6)
        self.assertEqual(ts.num_nodes, 13)

    def test_ancestral_samples(self):
        # Check that specifying samples to be not at time 0.0 works.
        #
        # 1.0             7
        # 0.7            / \                      8                     6
        #               /   \                    / \                   / \
        # 0.5          /     5                  /   5                 /   5
        #             /     / \                /   / \               /   / \
        # 0.4        /     /   4              /   /   4             /   /   4
        #           /     /   / \            /   /   / \           /   /   / \
        # 0.2      /     /   3   \          3   /   /   \         /   /   3   \
        #         /     /    *    \         *  /   /     \       /   /    *    \
        # 0.0    0     1           2          1   0       2     0   1           2
        #              *           *          *           *         *           *
        #          (0.0, 0.2),                 (0.2, 0.8),         (0.8, 1.0)
        #
        # Simplified, keeping [1,2,3]
        #
        # 1.0
        # 0.7                                     5
        #                                        / \
        # 0.5                4                  /   4                     4
        #                   / \                /   / \                   / \
        # 0.4              /   3              /   /   3                 /   3
        #                 /   / \            /   /     \               /   / \
        # 0.2            /   2   \          2   /       \             /   2   \
        #               /    *    \         *  /         \           /    *    \
        # 0.0          0           1          0           1         0           1
        #              *           *          *           *         *           *
        #          (0.0, 0.2),                 (0.2, 0.8),         (0.8, 1.0)

        nodes = six.StringIO("""\
        id      is_sample   time
        0       0           0
        1       1           0
        2       1           0
        3       1           0.2
        4       0           0.4
        5       0           0.5
        6       0           0.7
        7       0           1.0
        8       0           0.8
        """)
        edgesets = six.StringIO("""\
        left    right   parent  children
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.2     0.8     8       3,5
        0.0     0.2     7       0,5
        """)
        first_ts = msprime.load_text(nodes=nodes, edgesets=edgesets)
        ts = do_simplify(first_ts)
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6: -1, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: 8, 4: 5, 5: 8, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        # maps [1,2,3] -> [0,1,2]
        true_simplified_trees = [
            {0: 4, 1: 3, 2: 3, 3: 4},
            {0: 4, 1: 4, 2: 5, 4: 5},
            {0: 4, 1: 3, 2: 3, 3: 4}]
        self.assertEqual(first_ts.sample_size, 3)
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(first_ts.num_trees, 3)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(first_ts.num_nodes, 9)
        self.assertEqual(ts.num_nodes, 6)
        self.assertEqual(first_ts.time(3), 0.2)
        self.assertEqual(ts.time(2), 0.2)
        # check topologies agree:
        tree_dicts = [t.parent_dict for t in first_ts.trees()]
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        tree_simplified_dicts = [t.parent_dict for t in ts.trees()]
        for a, t in zip(true_simplified_trees, tree_simplified_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        # check .simplify() works here
        self.verify_simplify_topology(first_ts, [1, 2, 3])

    def test_all_ancestral_samples(self):
        # Check that specifying samples all to be not at time 0.0 works.
        #
        # 1.0             7
        # 0.7            / \                      8                     6
        #               /   \                    / \                   / \
        # 0.5          /     5                  /   5                 /   5
        #             /     / \                /   / \               /   / \
        # 0.4        /     /   4              /   /   4             /   /   4
        #           /     /   / \            /   /   / \           /   /   / \
        # 0.2      /     /   3   \          3   /   /   \         /   /   3   \
        #         /     1    *    2         *  1   /     2       /   1    *    2
        # 0.0    0      *         *            *  0      *      0    *         *
        #
        #          (0.0, 0.2),                 (0.2, 0.8),         (0.8, 1.0)

        nodes = six.StringIO("""\
        id      is_sample   time
        0       0           0
        1       1           0.1
        2       1           0.1
        3       1           0.2
        4       0           0.4
        5       0           0.5
        6       0           0.7
        7       0           1.0
        8       0           0.8
        """)
        edgesets = six.StringIO("""\
        left    right   parent  children
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.2     0.8     8       3,5
        0.0     0.2     7       0,5
        """)
        ts = msprime.load_text(nodes=nodes, edgesets=edgesets)
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6: -1, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: 8, 4: 5, 5: 8, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(ts.num_nodes, 9)
        self.assertEqual(ts.time(0), 0.0)
        self.assertEqual(ts.time(1), 0.1)
        self.assertEqual(ts.time(2), 0.1)
        self.assertEqual(ts.time(3), 0.2)
        # check topologies agree:
        tree_dicts = [t.parent_dict for t in ts.trees()]
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        # check .simplify() works here
        self.verify_simplify_topology(ts, [1, 2, 3])

    def test_disconnected(self):
        # A complex situation with multiple roots.
        #
        # 6     5      7___    .   5       7___    .   5       7___    .   5       7     6 .
        #       |     / \  \   .   |      / \  \   .   |      / \  \   .   |      / \    | .
        # 5    12   13   10 11 .  12    13   10 11 .  12    13   10 11 .  12    13   10 11 .
        #          / |\    \   .   |   / |\    \   .   |\   /|\    \   .   |\   /|\    \   .
        # 4      14  | \   15  .   | 14  | \   15  .   | |14 | \   15  .   | |14 | \   15  .
        #        /\  |  \   |  .   |  \  |  \   |  .   | |   |  \   |  .   | |   |  \   |  .
        # 3    17 16 18  |  |  .  17  16 18  |  |  .  17 16  18  |  |  .  17 16  18  |  |  .
        #       |        |  |  .   |         |  |  .   |         |  |  .   |         |  |  .
        # 2    21       20 19  .   21       20 19  .   21       20 19  .   21       20 19  .
        #      /|        |  |  .   /|        |  |  .   /|        |  |  .   /|        |  |  .
        # 1  22 |        |  |  . 22 |        |  |  . 22 |        |  |  . 22 |        |  |  .
        #    /| |        |  |  . /| |        |  |  . /| |        |  |  . /| |        |  |  .
        # 0 4 2 3        1  0  .4 2 3        1  0  .4 2 3        1  0  .4 2 3        1  0  .
        #
        #   0.0 ------------- 0.18 -------------- 0.21 -------------- 0.34 ------------ 0.43
        #
        # ... continued:
        #
        # 6   .   5      6   7    .       5      6   7    .       5      6   7    .
        #     .   |     / \  |    .       |     / \  |    .       |     / \  |    .
        # 5   .  12    13 11 10   .      12    13 11 10   .      12    13 11 10   .
        #     .   |\   /|\    \   .     / |\   /|\    \   .     / |\   /|     \   .
        # 4   .   | |14 | \   15  .    /  | |14 | \   15  .    /  | |14 |     15  .
        #     .   | |   |  \   |  .   /   | |   |  \      .   /   | |   |         .
        # 3   .  17 16  18  |  |  .  |   17 16  18  |     .  /   17 16  18        .
        #     .   |         |  |  .  |    |         |     . /     |\              .
        # 2   .   21       20 19  . 19    21       20     .19    21 20            .
        #     .   /|        |  |  .  |    /|        |     . |    /|  |            .
        # 1   . 22 |        |  |  .  |  22 |        |     . |  22 |  |            .
        #     . /| |        |  |  .  |  /| |        |     . |  /| |  |            .
        # 0   .4 2 3        1  0  .  0 4 2 3        1     . 0 4 2 3  1            .
        #
        #   0.43 --- 0.47 ----- 0.75 ------------------ 0.98 ------------------ 1.0

        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   1           0
        3   1           0
        4   1           0
        5   0           6
        6   0           6
        7   0           6
        8   0           6
        9   0           6
        10  0           5
        11  0           5
        12  0           5
        13  0           5
        14  0           4
        15  0           4
        16  0           3
        17  0           3
        18  0           3
        19  0           2
        20  0           2
        21  0           2
        22  0           1
        """)
        edgesets = six.StringIO("""\
        left    right   parent  children
        0.0     1.0     22      2,4
        0.0     1.0     21      3,22
        0.0     1.0     20      1
        0.0     1.0     19      0
        0.0     0.34    17      21
        0.34    0.98    17      21
        0.98    1.0     17      20,21
        0.0     0.75    15      19
        0.0     0.18    14      16,17
        0.18    0.21    14      16
        0.0     0.98    13      14,18,20
        0.98    1.0     13      14,18
        0.18    0.21    12      17
        0.21    0.75    12      16,17
        0.75    1.0     12      16,17,19
        0.0     1.0     10      15
        0.0     0.34    7       10,11,13
        0.34    0.43    7       10,13
        0.43    1.0     7       10
        0.34    0.43    6       11
        0.43    1.0     6       11,13
        0.0     1.0     5       12
        """)
        ts = msprime.load_text(nodes=nodes, edgesets=edgesets)
        self.verify_simplify_topology(ts, [0, 1, 2, 3, 4])
        self.verify_simplify_topology(ts, [0, 1])
        self.verify_simplify_topology(ts, [1, 2])
        self.verify_simplify_topology(ts, [2, 0])
