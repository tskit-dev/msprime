#
# Copyright (C) 2016-2017 University of Oxford
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
Test cases for the supported topological variations and operations.
"""
from __future__ import print_function
from __future__ import division

try:
    # We use the zip as iterator functionality here.
    from future_builtins import zip
except ImportError:
    # This fails for Python 3.x, but that's fine.
    pass

import unittest
import itertools
import random

import six
import numpy as np

import msprime
import _msprime
import tests
import tests.tsutil as tsutil
import tests.test_wright_fisher as wf


def generate_segments(n, sequence_length=100, seed=None):
    rng = random.Random(seed)
    segs = []
    for j in range(n):
        left = rng.randint(0, sequence_length - 1)
        right = rng.randint(left + 1, sequence_length)
        assert left < right
        segs.append(tests.Segment(left, right, j))
    return segs


class TestOverlappingSegments(unittest.TestCase):
    """
    Tests for the overlapping segments algorithm required for simplify.
    This test probably belongs somewhere else.
    """

    def test_random(self):
        segs = generate_segments(10, 20, 1)
        for left, right, X in tests.overlapping_segments(segs):
            self.assertGreater(right, left)
            self.assertGreater(len(X), 0)

    def test_empty(self):
        ret = list(tests.overlapping_segments([]))
        self.assertEqual(len(ret), 0)

    def test_single_interval(self):
        for j in range(1, 10):
            segs = [tests.Segment(0, 1, j) for _ in range(j)]
            ret = list(tests.overlapping_segments(segs))
            self.assertEqual(len(ret), 1)
            left, right, X = ret[0]
            self.assertEqual(left, 0)
            self.assertEqual(right, 1)
            self.assertEqual(sorted(segs), sorted(X))

    def test_stairs_down(self):
        segs = [
            tests.Segment(0, 1, 0),
            tests.Segment(0, 2, 1),
            tests.Segment(0, 3, 2)]
        ret = list(tests.overlapping_segments(segs))
        self.assertEqual(len(ret), 3)

        left, right, X = ret[0]
        self.assertEqual(left, 0)
        self.assertEqual(right, 1)
        self.assertEqual(sorted(X), sorted(segs))

        left, right, X = ret[1]
        self.assertEqual(left, 1)
        self.assertEqual(right, 2)
        self.assertEqual(sorted(X), sorted(segs[1:]))

        left, right, X = ret[2]
        self.assertEqual(left, 2)
        self.assertEqual(right, 3)
        self.assertEqual(sorted(X), sorted(segs[2:]))

    def test_stairs_up(self):
        segs = [
            tests.Segment(0, 3, 0),
            tests.Segment(1, 3, 1),
            tests.Segment(2, 3, 2)]
        ret = list(tests.overlapping_segments(segs))
        self.assertEqual(len(ret), 3)

        left, right, X = ret[0]
        self.assertEqual(left, 0)
        self.assertEqual(right, 1)
        self.assertEqual(X, segs[:1])

        left, right, X = ret[1]
        self.assertEqual(left, 1)
        self.assertEqual(right, 2)
        self.assertEqual(sorted(X), sorted(segs[:2]))

        left, right, X = ret[2]
        self.assertEqual(left, 2)
        self.assertEqual(right, 3)
        self.assertEqual(sorted(X), sorted(segs))

    def test_pyramid(self):
        segs = [
            tests.Segment(0, 5, 0),
            tests.Segment(1, 4, 1),
            tests.Segment(2, 3, 2)]
        ret = list(tests.overlapping_segments(segs))
        self.assertEqual(len(ret), 5)

        left, right, X = ret[0]
        self.assertEqual(left, 0)
        self.assertEqual(right, 1)
        self.assertEqual(X, segs[:1])

        left, right, X = ret[1]
        self.assertEqual(left, 1)
        self.assertEqual(right, 2)
        self.assertEqual(sorted(X), sorted(segs[:2]))

        left, right, X = ret[2]
        self.assertEqual(left, 2)
        self.assertEqual(right, 3)
        self.assertEqual(sorted(X), sorted(segs))

        left, right, X = ret[3]
        self.assertEqual(left, 3)
        self.assertEqual(right, 4)
        self.assertEqual(sorted(X), sorted(segs[:2]))

        left, right, X = ret[4]
        self.assertEqual(left, 4)
        self.assertEqual(right, 5)
        self.assertEqual(sorted(X), sorted(segs[:1]))

    def test_gap(self):
        segs = [
            tests.Segment(0, 2, 0),
            tests.Segment(3, 4, 1)]
        ret = list(tests.overlapping_segments(segs))
        self.assertEqual(len(ret), 2)

        left, right, X = ret[0]
        self.assertEqual(left, 0)
        self.assertEqual(right, 2)
        self.assertEqual(X, segs[:1])

        left, right, X = ret[1]
        self.assertEqual(left, 3)
        self.assertEqual(right, 4)
        self.assertEqual(X, segs[1:])


class TopologyTestCase(unittest.TestCase):
    """
    Superclass of test cases containing common utilities.
    """
    random_seed = 123456

    def assert_haplotypes_equal(self, ts1, ts2):
        h1 = list(ts1.haplotypes())
        h2 = list(ts2.haplotypes())
        self.assertEqual(h1, h2)

    def assert_variants_equal(self, ts1, ts2):
        v1 = list(ts1.variants(as_bytes=True))
        v2 = list(ts2.variants(as_bytes=True))
        self.assertEqual(v1, v2)

    def check_num_samples(self, ts, x):
        """
        Compare against x, a list of tuples of the form
        `(tree number, parent, number of samples)`.
        """
        k = 0
        tss = ts.trees(sample_counts=True)
        t = next(tss)
        for j, node, nl in x:
            while k < j:
                t = next(tss)
                k += 1
            self.assertEqual(nl, t.num_samples(node))

    def check_num_tracked_samples(self, ts, tracked_samples, x):
        k = 0
        tss = ts.trees(sample_counts=True, tracked_samples=tracked_samples)
        t = next(tss)
        for j, node, nl in x:
            while k < j:
                t = next(tss)
                k += 1
            self.assertEqual(nl, t.num_tracked_samples(node))

    def check_sample_iterator(self, ts, x):
        """
        Compare against x, a list of tuples of the form
        `(tree number, node, sample ID list)`.
        """
        k = 0
        tss = ts.trees(sample_lists=True)
        t = next(tss)
        for j, node, samples in x:
            while k < j:
                t = next(tss)
                k += 1
            for u, v in zip(samples, t.samples(node)):
                self.assertEqual(u, v)


class TestZeroRoots(unittest.TestCase):
    """
    Tests that for the case in which we have zero samples and therefore
    zero roots in our trees.
    """
    def remove_samples(self, ts):
        tables = ts.dump_tables()
        tables.nodes.set_columns(
            flags=np.zeros_like(tables.nodes.flags),
            time=tables.nodes.time)
        return tables.tree_sequence()

    def verify(self, ts, no_root_ts):
        self.assertEqual(ts.num_trees, no_root_ts.num_trees)
        for tree, no_root in zip(ts.trees(), no_root_ts.trees()):
            self.assertEqual(no_root.num_roots, 0)
            self.assertEqual(no_root.left_root, msprime.NULL_NODE)
            self.assertEqual(no_root.roots, [])
            self.assertEqual(tree.parent_dict, no_root.parent_dict)

    def test_single_tree(self):
        ts = msprime.simulate(10, random_seed=1)
        no_root_ts = self.remove_samples(ts)
        self.assertEqual(ts.num_trees, 1)
        self.verify(ts, no_root_ts)

    def test_multiple_trees(self):
        ts = msprime.simulate(10, recombination_rate=2, random_seed=1)
        no_root_ts = self.remove_samples(ts)
        self.assertGreater(ts.num_trees, 1)
        self.verify(ts, no_root_ts)


class TestEmptyTreeSequences(TopologyTestCase):
    """
    Tests covering tree sequences that have zero edges.
    """
    def test_zero_nodes(self):
        nodes = msprime.NodeTable()
        edges = msprime.EdgeTable()
        # Without a sequence length this should fail.
        self.assertRaises(
            _msprime.LibraryError, msprime.load_tables, nodes=nodes, edges=edges)
        ts = msprime.load_tables(nodes=nodes, edges=edges, sequence_length=1)
        self.assertEqual(ts.sequence_length, 1)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 0)
        self.assertEqual(ts.num_edges, 0)
        t = next(ts.trees())
        self.assertEqual(t.index, 0)
        self.assertEqual(t.left_root, msprime.NULL_NODE)
        self.assertEqual(t.interval, (0, 1))
        self.assertEqual(t.roots, [])
        self.assertEqual(t.root, msprime.NULL_NODE)
        self.assertEqual(t.parent_dict, {})
        self.assertEqual(list(t.nodes()), [])
        self.assertEqual(list(ts.haplotypes()), [])
        self.assertEqual(list(ts.variants()), [])
        methods = [t.parent, t.left_child, t.right_child, t.left_sib, t.right_sib]
        for method in methods:
            for u in [-1, 0, 1, 100]:
                self.assertRaises(ValueError, method, u)
        tsp = ts.simplify()
        self.assertEqual(tsp.num_nodes, 0)
        self.assertEqual(tsp.num_edges, 0)

    def test_one_node_zero_samples(self):
        nodes = msprime.NodeTable()
        nodes.add_row(time=0, flags=0)
        edges = msprime.EdgeTable()
        # Without a sequence length this should fail.
        self.assertRaises(
            _msprime.LibraryError, msprime.load_tables, nodes=nodes, edges=edges)
        ts = msprime.load_tables(nodes=nodes, edges=edges, sequence_length=1)
        self.assertEqual(ts.sequence_length, 1)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 1)
        self.assertEqual(ts.sample_size, 0)
        self.assertEqual(ts.num_edges, 0)
        self.assertEqual(ts.num_sites, 0)
        self.assertEqual(ts.num_mutations, 0)
        t = next(ts.trees())
        self.assertEqual(t.index, 0)
        self.assertEqual(t.left_root, msprime.NULL_NODE)
        self.assertEqual(t.interval, (0, 1))
        self.assertEqual(t.roots, [])
        self.assertEqual(t.root, msprime.NULL_NODE)
        self.assertEqual(t.parent_dict, {})
        self.assertEqual(list(t.nodes()), [])
        self.assertEqual(list(ts.haplotypes()), [])
        self.assertEqual(list(ts.variants()), [])
        methods = [t.parent, t.left_child, t.right_child, t.left_sib, t.right_sib]
        for method in methods:
            self.assertEqual(method(0), msprime.NULL_NODE)
            for u in [-1, 1, 100]:
                self.assertRaises(ValueError, method, u)

    def test_one_node_zero_samples_sites(self):
        nodes = msprime.NodeTable()
        nodes.add_row(time=0, flags=0)
        edges = msprime.EdgeTable()
        sites = msprime.SiteTable()
        mutations = msprime.MutationTable()
        sites.add_row(position=0.5, ancestral_state='0')
        mutations.add_row(site=0, derived_state='1', node=0)
        ts = msprime.load_tables(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations,
            sequence_length=1)
        self.assertEqual(ts.sequence_length, 1)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 1)
        self.assertEqual(ts.sample_size, 0)
        self.assertEqual(ts.num_edges, 0)
        self.assertEqual(ts.num_sites, 1)
        self.assertEqual(ts.num_mutations, 1)
        t = next(ts.trees())
        self.assertEqual(t.index, 0)
        self.assertEqual(t.left_root, msprime.NULL_NODE)
        self.assertEqual(t.interval, (0, 1))
        self.assertEqual(t.roots, [])
        self.assertEqual(t.root, msprime.NULL_NODE)
        self.assertEqual(t.parent_dict, {})
        self.assertEqual(len(list(t.sites())), 1)
        self.assertEqual(list(t.nodes()), [])
        self.assertEqual(list(ts.haplotypes()), [])
        self.assertEqual(len(list(ts.variants())), 1)
        tsp = ts.simplify()
        self.assertEqual(tsp.num_nodes, 0)
        self.assertEqual(tsp.num_edges, 0)

    def test_one_node_one_sample(self):
        nodes = msprime.NodeTable()
        nodes.add_row(time=0, flags=msprime.NODE_IS_SAMPLE)
        edges = msprime.EdgeTable()
        # Without a sequence length this should fail.
        self.assertRaises(
            _msprime.LibraryError, msprime.load_tables, nodes=nodes, edges=edges)
        ts = msprime.load_tables(nodes=nodes, edges=edges, sequence_length=1)
        self.assertEqual(ts.sequence_length, 1)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 1)
        self.assertEqual(ts.sample_size, 1)
        self.assertEqual(ts.num_edges, 0)
        t = next(ts.trees())
        self.assertEqual(t.index, 0)
        self.assertEqual(t.left_root, 0)
        self.assertEqual(t.interval, (0, 1))
        self.assertEqual(t.roots, [0])
        self.assertEqual(t.root, 0)
        self.assertEqual(t.parent_dict, {})
        self.assertEqual(list(t.nodes()), [0])
        self.assertEqual(list(ts.haplotypes()), [""])
        self.assertEqual(list(ts.variants()), [])
        methods = [t.parent, t.left_child, t.right_child, t.left_sib, t.right_sib]
        for method in methods:
            self.assertEqual(method(0), msprime.NULL_NODE)
            for u in [-1, 1, 100]:
                self.assertRaises(ValueError, method, u)
        tsp = ts.simplify()
        self.assertEqual(tsp.num_nodes, 1)
        self.assertEqual(tsp.num_edges, 0)

    def test_one_node_one_sample_sites(self):
        nodes = msprime.NodeTable()
        nodes.add_row(time=0, flags=msprime.NODE_IS_SAMPLE)
        edges = msprime.EdgeTable()
        sites = msprime.SiteTable()
        mutations = msprime.MutationTable()
        sites.add_row(position=0.5, ancestral_state='0')
        mutations.add_row(site=0, derived_state='1', node=0)
        ts = msprime.load_tables(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations,
            sequence_length=1)
        self.assertEqual(ts.sequence_length, 1)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 1)
        self.assertEqual(ts.sample_size, 1)
        self.assertEqual(ts.num_edges, 0)
        self.assertEqual(ts.num_sites, 1)
        self.assertEqual(ts.num_mutations, 1)
        t = next(ts.trees())
        self.assertEqual(t.index, 0)
        self.assertEqual(t.left_root, 0)
        self.assertEqual(t.interval, (0, 1))
        self.assertEqual(t.roots, [0])
        self.assertEqual(t.root, 0)
        self.assertEqual(t.parent_dict, {})
        self.assertEqual(list(t.nodes()), [0])
        self.assertEqual(list(ts.haplotypes()), ["1"])
        self.assertEqual(len(list(ts.variants())), 1)
        methods = [t.parent, t.left_child, t.right_child, t.left_sib, t.right_sib]
        for method in methods:
            self.assertEqual(method(0), msprime.NULL_NODE)
            for u in [-1, 1, 100]:
                self.assertRaises(ValueError, method, u)
        tsp = ts.simplify(filter_sites=False)
        self.assertEqual(tsp.num_nodes, 1)
        self.assertEqual(tsp.num_edges, 0)
        self.assertEqual(tsp.num_sites, 1)


class TestHoleyTreeSequences(TopologyTestCase):
    """
    Tests for tree sequences in which we have partial (or no) trees defined
    over some of the sequence.
    """
    def verify_trees(self, ts, expected):
        observed = []
        for t in ts.trees():
            observed.append((t.interval, t.parent_dict))
        self.assertEqual(expected, observed)
        # Test simple algorithm also.
        observed = []
        for interval, parent in tsutil.algorithm_T(ts):
            parent_dict = {j: parent[j] for j in range(ts.num_nodes) if parent[j] >= 0}
            observed.append((interval, parent_dict))
        self.assertEqual(expected, observed)

    def verify_zero_roots(self, ts):
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 0)
            self.assertEqual(tree.left_root, msprime.NULL_NODE)
            self.assertEqual(tree.roots, [])

    def test_simple_hole(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       2       0
        2       3       2       0
        0       1       2       1
        2       3       2       1
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        expected = [
            ((0, 1), {0: 2, 1: 2}),
            ((1, 2), {}),
            ((2, 3), {0: 2, 1: 2})]
        self.verify_trees(ts, expected)

    def test_simple_hole_zero_roots(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   0           0
        1   0           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       2       0
        2       3       2       0
        0       1       2       1
        2       3       2       1
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        expected = [
            ((0, 1), {0: 2, 1: 2}),
            ((1, 2), {}),
            ((2, 3), {0: 2, 1: 2})]
        self.verify_trees(ts, expected)
        self.verify_zero_roots(ts)

    def test_initial_gap(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        1       2       2       0,1
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        expected = [
            ((0, 1), {}),
            ((1, 2), {0: 2, 1: 2})]
        self.verify_trees(ts, expected)

    def test_initial_gap_zero_roots(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   0           0
        1   0           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        1       2       2       0,1
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        expected = [
            ((0, 1), {}),
            ((1, 2), {0: 2, 1: 2})]
        self.verify_trees(ts, expected)
        self.verify_zero_roots(ts)

    def test_final_gap(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       2       2       0,1
        """)
        ts = msprime.load_text(nodes, edges, sequence_length=3, strict=False)
        expected = [
            ((0, 2), {0: 2, 1: 2}),
            ((2, 3), {})]
        self.verify_trees(ts, expected)

    def test_final_gap_zero_roots(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   0           0
        1   0           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       2       2       0,1
        """)
        ts = msprime.load_text(nodes, edges, sequence_length=3, strict=False)
        expected = [
            ((0, 2), {0: 2, 1: 2}),
            ((2, 3), {})]
        self.verify_trees(ts, expected)
        self.verify_zero_roots(ts)

    def test_initial_and_final_gap(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        1       2       2       0,1
        """)
        ts = msprime.load_text(nodes, edges, sequence_length=3, strict=False)
        expected = [
            ((0, 1), {}),
            ((1, 2), {0: 2, 1: 2}),
            ((2, 3), {})]
        self.verify_trees(ts, expected)

    def test_initial_and_final_gap_zero_roots(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   0           0
        1   0           0
        2   0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        1       2       2       0,1
        """)
        ts = msprime.load_text(nodes, edges, sequence_length=3, strict=False)
        expected = [
            ((0, 1), {}),
            ((1, 2), {0: 2, 1: 2}),
            ((2, 3), {})]
        self.verify_trees(ts, expected)
        self.verify_zero_roots(ts)


class TestTsinferExamples(TopologyTestCase):
    """
    Test cases on troublesome topology examples that arose from tsinfer.
    """
    def test_no_last_tree(self):
        # The last tree was not being generated here because of a bug in
        # the low-level tree generation code.
        nodes = six.StringIO("""\
        id      is_sample   population      time
        0       1       -1              3.00000000000000
        1       1       -1              2.00000000000000
        2       1       -1              2.00000000000000
        3       1       -1              2.00000000000000
        4       1       -1              2.00000000000000
        5       1       -1              1.00000000000000
        6       1       -1              1.00000000000000
        7       1       -1              1.00000000000000
        8       1       -1              1.00000000000000
        9       1       -1              1.00000000000000
        10      1       -1              1.00000000000000
        """)
        edges = six.StringIO("""\
        id      left            right           parent  child
        0       62291.41659631  79679.17408763  1       5
        1       62291.41659631  62374.60889677  1       6
        2       122179.36037089 138345.43104411 1       7
        3       67608.32330402  79679.17408763  1       8
        4       122179.36037089 138345.43104411 1       8
        5       62291.41659631  79679.17408763  1       9
        6       126684.47550333 138345.43104411 1       10
        7       23972.05905068  62291.41659631  2       5
        8       79679.17408763  82278.53390076  2       5
        9       23972.05905068  62291.41659631  2       6
        10      79679.17408763  110914.43816806 2       7
        11      145458.28890561 189765.31932273 2       7
        12      79679.17408763  110914.43816806 2       8
        13      145458.28890561 200000.00000000 2       8
        14      23972.05905068  62291.41659631  2       9
        15      79679.17408763  110914.43816806 2       9
        16      145458.28890561 145581.18329797 2       10
        17      4331.62138785   23972.05905068  3       6
        18      4331.62138785   23972.05905068  3       9
        19      110914.43816806 122179.36037089 4       7
        20      138345.43104411 145458.28890561 4       7
        21      110914.43816806 122179.36037089 4       8
        22      138345.43104411 145458.28890561 4       8
        23      110914.43816806 112039.30503475 4       9
        24      138345.43104411 145458.28890561 4       10
        25      0.00000000      200000.00000000 0       1
        26      0.00000000      200000.00000000 0       2
        27      0.00000000      200000.00000000 0       3
        28      0.00000000      200000.00000000 0       4
        """)
        ts = msprime.load_text(nodes, edges, sequence_length=200000, strict=False)
        pts = tests.PythonTreeSequence(ts.get_ll_tree_sequence())
        num_trees = 0
        for t in pts.trees():
            num_trees += 1
        self.assertEqual(num_trees, ts.num_trees)
        n = 0
        for pt, t in zip(pts.trees(), ts.trees()):
            self.assertEqual((pt.left, pt.right), t.interval)
            for j in range(ts.num_nodes):
                self.assertEqual(pt.parent[j], t.parent(j))
                self.assertEqual(pt.left_child[j], t.left_child(j))
                self.assertEqual(pt.right_child[j], t.right_child(j))
                self.assertEqual(pt.left_sib[j], t.left_sib(j))
                self.assertEqual(pt.right_sib[j], t.right_sib(j))
            n += 1
        self.assertEqual(n, num_trees)
        intervals = [t.interval for t in ts.trees()]
        self.assertEqual(intervals[0][0], 0)
        self.assertEqual(intervals[-1][-1], ts.sequence_length)


class TestRecordSquashing(TopologyTestCase):
    """
    Tests that we correctly squash adjacent equal records together.
    """
    def test_single_record(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       1       0
        1       2       1       0
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        tss, node_map = ts.simplify(map_nodes=True)
        self.assertEqual(list(node_map), [0, 1])
        self.assertEqual(tss.dump_tables().nodes, ts.dump_tables().nodes)
        simplified_edges = list(tss.edges())
        self.assertEqual(len(simplified_edges), 1)
        e = simplified_edges[0]
        self.assertEqual(e.left, 0)
        self.assertEqual(e.right, 2)

    def test_single_tree(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        ts_redundant = tsutil.insert_redundant_breakpoints(ts)
        tss = ts_redundant.simplify()
        self.assertEqual(tss.dump_tables().nodes, ts.dump_tables().nodes)
        self.assertEqual(tss.dump_tables().edges, ts.dump_tables().edges)

    def test_many_trees(self):
        ts = msprime.simulate(
            20, recombination_rate=5, random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        ts_redundant = tsutil.insert_redundant_breakpoints(ts)
        tss = ts_redundant.simplify()
        self.assertEqual(tss.dump_tables().nodes, ts.dump_tables().nodes)
        self.assertEqual(tss.dump_tables().edges, ts.dump_tables().edges)


class TestRedundantBreakpoints(TopologyTestCase):
    """
    Tests for dealing with redundant breakpoints within the tree sequence.
    These are records that may be squashed together into a single record.
    """
    def test_single_tree(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        ts_redundant = tsutil.insert_redundant_breakpoints(ts)
        self.assertEqual(ts.sample_size, ts_redundant.sample_size)
        self.assertEqual(ts.sequence_length, ts_redundant.sequence_length)
        self.assertEqual(ts_redundant.num_trees, 2)
        trees = [t.parent_dict for t in ts_redundant.trees()]
        self.assertEqual(len(trees), 2)
        self.assertEqual(trees[0], trees[1])
        self.assertEqual([t.parent_dict for t in ts.trees()][0], trees[0])

    def test_many_trees(self):
        ts = msprime.simulate(
            20, recombination_rate=5, random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        ts_redundant = tsutil.insert_redundant_breakpoints(ts)
        self.assertEqual(ts.sample_size, ts_redundant.sample_size)
        self.assertEqual(ts.sequence_length, ts_redundant.sequence_length)
        self.assertGreater(ts_redundant.num_trees, ts.num_trees)
        self.assertGreater(ts_redundant.num_edges, ts.num_edges)
        redundant_trees = ts_redundant.trees()
        redundant_t = next(redundant_trees)
        comparisons = 0
        for t in ts.trees():
            while redundant_t is not None and redundant_t.interval[1] <= t.interval[1]:
                self.assertEqual(t.parent_dict, redundant_t.parent_dict)
                comparisons += 1
                redundant_t = next(redundant_trees, None)
        self.assertEqual(comparisons, ts_redundant.num_trees)


class TestUnaryNodes(TopologyTestCase):
    """
    Tests for situations in which we have unary nodes in the tree sequence.
    """
    def test_simple_case(self):
        # Simple case where we have n = 2 and some unary nodes.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        3       0           1
        4       0           2
        5       0           3
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       2       0
        0       1       3       1
        0       1       4       2,3
        0       1       5       4
        """)
        sites = "position    ancestral_state\n"
        mutations = "site    node    derived_state\n"
        for j in range(5):
            position = j * 1 / 5
            sites += "{} 0\n".format(position)
            mutations += "{} {} 1\n".format(j, j)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=six.StringIO(sites),
            mutations=six.StringIO(mutations), strict=False)

        self.assertEqual(ts.sample_size, 2)
        self.assertEqual(ts.num_nodes, 6)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_sites, 5)
        self.assertEqual(ts.num_mutations, 5)
        self.assertEqual(len(list(ts.edge_diffs())), ts.num_trees)
        t = next(ts.trees())
        self.assertEqual(
            t.parent_dict, {0: 2, 1: 3, 2: 4, 3: 4, 4: 5})
        self.assertEqual(t.mrca(0, 1), 4)
        self.assertEqual(t.mrca(0, 2), 2)
        self.assertEqual(t.mrca(0, 4), 4)
        self.assertEqual(t.mrca(0, 5), 5)
        self.assertEqual(t.mrca(0, 3), 4)
        H = list(ts.haplotypes())
        self.assertEqual(H[0], "10101")
        self.assertEqual(H[1], "01011")

    def test_ladder_tree(self):
        # We have a single tree with a long ladder of unary nodes along a path
        num_unary_nodes = 30
        n = 2
        nodes = """\
            is_sample   time
            1           0
            1           0
        """
        edges = """\
            left right parent child
            0    1     2      0
        """
        for j in range(num_unary_nodes + 2):
            nodes += "0 {}\n".format(j + 2)
        for j in range(num_unary_nodes):
            edges += "0 1 {} {}\n".format(n + j + 1, n + j)
        root = num_unary_nodes + 3
        root_time = num_unary_nodes + 3
        edges += "0    1     {}      1,{}\n".format(root, num_unary_nodes + 2)
        ts = msprime.load_text(six.StringIO(nodes), six.StringIO(edges), strict=False)
        t = next(ts.trees())
        self.assertEqual(t.mrca(0, 1), root)
        self.assertEqual(t.tmrca(0, 1), root_time)
        ts_simplified, node_map = ts.simplify(map_nodes=True)
        test_map = [msprime.NULL_NODE for _ in range(ts.num_nodes)]
        test_map[0] = 0
        test_map[1] = 1
        test_map[root] = 2
        self.assertEqual(list(node_map), test_map)
        self.assertEqual(ts_simplified.num_edges, 2)
        t = next(ts_simplified.trees())
        self.assertEqual(t.mrca(0, 1), 2)
        self.assertEqual(t.tmrca(0, 1), root_time)

    def verify_unary_tree_sequence(self, ts):
        """
        Take the specified tree sequence and produce an equivalent in which
        unary records have been interspersed.
        """
        self.assertGreater(ts.num_trees, 2)
        self.assertGreater(ts.num_mutations, 2)
        tables = ts.dump_tables()
        next_node = ts.num_nodes
        node_times = {j: node.time for j, node in enumerate(ts.nodes())}
        edges = []
        for e in ts.edges():
            node = ts.node(e.parent)
            t = node.time - 1e-14  # Arbitrary small value.
            next_node = len(tables.nodes)
            tables.nodes.add_row(time=t, population=node.population)
            edges.append(msprime.Edge(
                left=e.left, right=e.right, parent=next_node, child=e.child))
            node_times[next_node] = t
            edges.append(msprime.Edge(
                left=e.left, right=e.right, parent=e.parent, child=next_node))
        edges.sort(key=lambda e: node_times[e.parent])
        tables.edges.reset()
        for e in edges:
            tables.edges.add_row(
                left=e.left, right=e.right, child=e.child, parent=e.parent)
        ts_new = msprime.load_tables(**tables.asdict())
        self.assertGreater(ts_new.num_edges, ts.num_edges)
        self.assert_haplotypes_equal(ts, ts_new)
        self.assert_variants_equal(ts, ts_new)
        ts_simplified = ts_new.simplify()
        self.assertEqual(list(ts_simplified.records()), list(ts.records()))
        self.assert_haplotypes_equal(ts, ts_simplified)
        self.assert_variants_equal(ts, ts_simplified)
        self.assertEqual(len(list(ts.edge_diffs())), ts.num_trees)

    def test_binary_tree_sequence_unary_nodes(self):
        ts = msprime.simulate(
            20, recombination_rate=5, mutation_rate=5, random_seed=self.random_seed)
        self.verify_unary_tree_sequence(ts)

    def test_nonbinary_tree_sequence_unary_nodes(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=1.0, population=0, proportion=0.95)]
        ts = msprime.simulate(
            20, recombination_rate=10, mutation_rate=5,
            demographic_events=demographic_events, random_seed=self.random_seed)
        found = False
        for r in ts.edgesets():
            if len(r.children) > 2:
                found = True
        self.assertTrue(found)
        self.verify_unary_tree_sequence(ts)


class TestGeneralSamples(TopologyTestCase):
    """
    Test cases in which we have samples at arbitrary nodes (i.e., not at
    {0,...,n - 1}).
    """
    def test_simple_case(self):
        # Simple case where we have n = 3 and samples starting at n.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       0           2
        1       0           1
        2       1           0
        3       1           0
        4       1           0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       1       2,3
        0       1       0       1,4
        """)
        sites = six.StringIO("""\
        position    ancestral_state
        0.1     0
        0.2     0
        0.3     0
        0.4     0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state
        0       2       1
        1       3       1
        2       4       1
        3       1       1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)

        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(list(ts.samples()), [2, 3, 4])
        self.assertEqual(ts.num_nodes, 5)
        self.assertEqual(ts.num_nodes, 5)
        self.assertEqual(ts.num_sites, 4)
        self.assertEqual(ts.num_mutations, 4)
        self.assertEqual(len(list(ts.edge_diffs())), ts.num_trees)
        t = next(ts.trees())
        self.assertEqual(t.root, 0)
        self.assertEqual(t.parent_dict, {1: 0, 2: 1, 3: 1, 4: 0})
        H = list(ts.haplotypes())
        self.assertEqual(H[0], "1001")
        self.assertEqual(H[1], "0101")
        self.assertEqual(H[2], "0010")

        tss, node_map = ts.simplify(map_nodes=True)
        self.assertEqual(list(node_map), [4, 3, 0, 1, 2])
        # We should have the same tree sequence just with canonicalised nodes.
        self.assertEqual(tss.sample_size, 3)
        self.assertEqual(list(tss.samples()), [0, 1, 2])
        self.assertEqual(tss.num_nodes, 5)
        self.assertEqual(tss.num_trees, 1)
        self.assertEqual(tss.num_sites, 4)
        self.assertEqual(tss.num_mutations, 4)
        self.assertEqual(len(list(ts.edge_diffs())), ts.num_trees)
        t = next(tss.trees())
        self.assertEqual(t.root, 4)
        self.assertEqual(t.parent_dict, {0: 3, 1: 3, 2: 4, 3: 4})
        H = list(tss.haplotypes())
        self.assertEqual(H[0], "1001")
        self.assertEqual(H[1], "0101")
        self.assertEqual(H[2], "0010")

    def verify_permuted_nodes(self, ts):
        """
        Take the specified tree sequence and permute the nodes, verifying that we
        get back a tree sequence with the correct properties.
        """
        # Mapping from the original nodes into nodes in the new tree sequence.
        node_map = list(range(ts.num_nodes))
        random.shuffle(node_map)
        # Change the permutation so that the relative order of samples is maintained.
        # Then, we should get back exactly the same tree sequence after simplify
        # and haplotypes and variants are also equal.
        samples = sorted(node_map[:ts.sample_size])
        node_map = samples + node_map[ts.sample_size:]
        permuted = tsutil.permute_nodes(ts, node_map)
        self.assertEqual(ts.sequence_length, permuted.sequence_length)
        self.assertEqual(list(permuted.samples()), samples)
        self.assertEqual(list(permuted.haplotypes()), list(ts.haplotypes()))
        self.assertEqual(
            [v.genotypes for v in permuted.variants(as_bytes=True)],
            [v.genotypes for v in ts.variants(as_bytes=True)])
        self.assertEqual(ts.num_trees, permuted.num_trees)
        j = 0
        for t1, t2 in zip(ts.trees(), permuted.trees()):
            t1_dict = {node_map[k]: node_map[v] for k, v in t1.parent_dict.items()}
            self.assertEqual(node_map[t1.root], t2.root)
            self.assertEqual(t1_dict, t2.parent_dict)
            for u1 in t1.nodes():
                u2 = node_map[u1]
                self.assertEqual(
                    sorted([node_map[v] for v in t1.samples(u1)]),
                    sorted(list(t2.samples(u2))))
            j += 1
        self.assertEqual(j, ts.num_trees)

        # The simplified version of the permuted tree sequence should be in canonical
        # form, and identical to the original.
        simplified, s_node_map = permuted.simplify(map_nodes=True)

        original_tables = ts.dump_tables()
        simplified_tables = simplified.dump_tables()
        original_tables.provenances.clear()
        simplified_tables.provenances.clear()

        self.assertEqual(
            original_tables.sequence_length, simplified_tables.sequence_length)
        self.assertEqual(original_tables.nodes, simplified_tables.nodes)
        self.assertEqual(original_tables.edges, simplified_tables.edges)
        self.assertEqual(original_tables.sites, simplified_tables.sites)
        self.assertEqual(original_tables.mutations, simplified_tables.mutations)
        self.assertEqual(original_tables.individuals, simplified_tables.individuals)
        self.assertEqual(original_tables.populations, simplified_tables.populations)

        self.assertEqual(original_tables, simplified_tables)
        self.assertEqual(ts.sequence_length, simplified.sequence_length)
        for tree in simplified.trees():
            pass

        for u, v in enumerate(node_map):
            self.assertEqual(s_node_map[v], u)
        self.assertTrue(np.array_equal(simplified.samples(), ts.samples()))
        self.assertEqual(list(simplified.nodes()), list(ts.nodes()))
        self.assertEqual(list(simplified.edges()), list(ts.edges()))
        self.assertEqual(list(simplified.sites()), list(ts.sites()))
        self.assertEqual(list(simplified.haplotypes()), list(ts.haplotypes()))
        self.assertEqual(
            list(simplified.variants(as_bytes=True)), list(ts.variants(as_bytes=True)))

    def test_single_tree_permuted_nodes(self):
        ts = msprime.simulate(10,  mutation_rate=5, random_seed=self.random_seed)
        self.verify_permuted_nodes(ts)

    def test_binary_tree_sequence_permuted_nodes(self):
        ts = msprime.simulate(
            20, recombination_rate=5, mutation_rate=5, random_seed=self.random_seed)
        self.verify_permuted_nodes(ts)

    def test_nonbinary_tree_sequence_permuted_nodes(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=1.0, population=0, proportion=0.95)]
        ts = msprime.simulate(
            20, recombination_rate=10, mutation_rate=5,
            demographic_events=demographic_events, random_seed=self.random_seed)
        found = False
        for e in ts.edgesets():
            if len(e.children) > 2:
                found = True
        self.assertTrue(found)
        self.verify_permuted_nodes(ts)


class TestSimplifyExamples(TopologyTestCase):
    """
    Tests for simplify where we write out the input and expected output
    or we detect expected errors.
    """
    def verify_simplify(
            self, samples, filter_sites=True,
            nodes_before=None, edges_before=None, sites_before=None,
            mutations_before=None, nodes_after=None, edges_after=None,
            sites_after=None, mutations_after=None, debug=False):
        """
        Verifies that if we run simplify on the specified input we get the
        required output.
        """
        ts = msprime.load_text(
            nodes=six.StringIO(nodes_before),
            edges=six.StringIO(edges_before),
            sites=six.StringIO(sites_before) if sites_before is not None else None,
            mutations=(
                six.StringIO(mutations_before)
                if mutations_before is not None else None),
            strict=False)
        before = ts.dump_tables()

        ts = msprime.load_text(
            nodes=six.StringIO(nodes_after),
            edges=six.StringIO(edges_after),
            sites=six.StringIO(sites_after) if sites_after is not None else None,
            mutations=(
                six.StringIO(mutations_after)
                if mutations_after is not None else None),
            strict=False,
            sequence_length=before.sequence_length)
        after = ts.dump_tables()

        ts = before.tree_sequence()
        # Make sure it's a valid topology. We want to be sure we evaluate the
        # whole iterator
        for t in ts.trees():
            self.assertTrue(t is not None)
        before.simplify(
            samples=samples, filter_sites=filter_sites)
        if debug:
            print("before")
            print(before)
            print("after")
            print(after)
        self.assertEqual(before, after)

    def test_unsorted_edges(self):
        # We have two nodes at the same time and interleave edges for
        # these nodes together. This is an error because all edges for
        # a given parent must be contigous.
        nodes_before = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        3       0           1
        """
        edges_before = """\
        left    right   parent  child
        0       1       2       0,1
        0       1       3       0,1
        1       2       2       0,1
        1       2       3       0,1
        """
        nodes = msprime.parse_nodes(six.StringIO(nodes_before), strict=False)
        edges = msprime.parse_edges(six.StringIO(edges_before), strict=False)
        self.assertRaises(
            _msprime.LibraryError, msprime.simplify_tables,
            samples=[0, 1], nodes=nodes, edges=edges)

    def test_single_binary_tree(self):
        #
        # 2        4
        #         / \
        # 1      3   \
        #       / \   \
        # 0   (0)(1)  (2)
        nodes_before = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           1
        4       0           2
        """
        edges_before = """\
        left    right   parent  child
        0       1       3       0,1
        0       1       4       2,3
        """
        # We sample 0 and 2, so we get
        nodes_after = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           2
        """
        edges_after = """\
        left    right   parent  child
        0       1       2       0,1
        """
        self.verify_simplify(
            samples=[0, 2],
            nodes_before=nodes_before, edges_before=edges_before,
            nodes_after=nodes_after, edges_after=edges_after)

    def test_single_binary_tree_internal_sample(self):
        #
        # 2        4
        #         / \
        # 1     (3)  \
        #       / \   \
        # 0   (0)  1  (2)
        nodes_before = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           0
        3       1           1
        4       0           2
        """
        edges_before = """\
        left    right   parent  child
        0       1       3       0,1
        0       1       4       2,3
        """
        # We sample 0 and 3, so we get
        nodes_after = """\
        id      is_sample   time
        0       1           0
        1       1           1
        """
        edges_after = """\
        left    right   parent  child
        0       1       1       0
        """
        self.verify_simplify(
            samples=[0, 3],
            nodes_before=nodes_before, edges_before=edges_before,
            nodes_after=nodes_after, edges_after=edges_after)

    def test_single_binary_tree_internal_sample_meet_at_root(self):
        # 3          5
        #           / \
        # 2        4  (6)
        #         / \
        # 1     (3)  \
        #       / \   \
        # 0   (0)  1   2
        nodes_before = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           0
        3       1           1
        4       0           2
        5       0           3
        6       1           2
        """
        edges_before = """\
        left    right   parent  child
        0       1       3       0,1
        0       1       4       2,3
        0       1       5       4,6
        """
        # We sample 0 and 3 and 6, so we get
        nodes_after = """\
        id      is_sample   time
        0       1           0
        1       1           1
        2       1           2
        3       0           3
        """
        edges_after = """\
        left    right   parent  child
        0       1       1       0
        0       1       3       1,2
        """
        self.verify_simplify(
            samples=[0, 3, 6],
            nodes_before=nodes_before, edges_before=edges_before,
            nodes_after=nodes_after, edges_after=edges_after)

    def test_single_binary_tree_simple_mutations(self):
        # 3          5
        #           / \
        # 2        4   \
        #         / \   s0
        # 1      3   s1  \
        #       / \   \   \
        # 0   (0) (1)  2  (6)
        nodes_before = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           0
        3       0           1
        4       0           2
        5       0           3
        6       1           0
        """
        edges_before = """\
        left    right   parent  child
        0       1       3       0,1
        0       1       4       2,3
        0       1       5       4,6
        """
        sites_before = """\
        id  position    ancestral_state
        0   0.1         0
        1   0.2         0
        """
        mutations_before = """\
        site    node    derived_state
        0       6       1
        1       2       1
        """

        # We sample 0 and 2 and 6, so we get
        nodes_after = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           1
        3       0           3
        """
        edges_after = """\
        left    right   parent  child
        0       1       3       0,1
        0       1       4       2,3
        """
        sites_after = """\
        id  position    ancestral_state
        0   0.1         0
        """
        mutations_after = """\
        site    node    derived_state
        0       2       1
        """
        self.verify_simplify(
            samples=[0, 1, 6],
            nodes_before=nodes_before, edges_before=edges_before,
            sites_before=sites_before, mutations_before=mutations_before,
            nodes_after=nodes_after, edges_after=edges_after,
            sites_after=sites_after, mutations_after=mutations_after)
        # If we don't filter the fixed sites, we should get the same
        # mutations and the original sites table back.
        self.verify_simplify(
            samples=[0, 1, 6], filter_sites=False,
            nodes_before=nodes_before, edges_before=edges_before,
            sites_before=sites_before, mutations_before=mutations_before,
            nodes_after=nodes_after, edges_after=edges_after,
            sites_after=sites_before, mutations_after=mutations_after)

    def test_overlapping_edges(self):
        nodes = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        """
        edges_before = """\
        left    right   parent  child
        0       2       2       0
        1       3       2       1
        """
        # We resolve the overlapping edges here. Since the flanking regions
        # have no interesting edges, these are left out of the output.
        edges_after = """\
        left    right   parent  child
        1       2       2       0,1
        """
        self.verify_simplify(
            samples=[0, 1],
            nodes_before=nodes, edges_before=edges_before,
            nodes_after=nodes, edges_after=edges_after)

    def test_overlapping_edges_internal_samples(self):
        nodes = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           1
        """
        edges = """\
        left    right   parent  child
        0       2       2       0
        1       3       2       1
        """
        self.verify_simplify(
            samples=[0, 1, 2],
            nodes_before=nodes, edges_before=edges, nodes_after=nodes, edges_after=edges)

    def test_unary_edges_no_overlap(self):
        nodes_before = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        """
        edges_before = """\
        left    right   parent  child
        0       2       2       0
        2       3       2       1
        """
        # Because there is no overlap between the samples, we just get an
        # empty set of output edges.
        nodes_after = """\
        id      is_sample   time
        0       1           0
        1       1           0
        """
        edges_after = """\
        left    right   parent  child
        """
        self.verify_simplify(
            samples=[0, 1],
            nodes_before=nodes_before, edges_before=edges_before,
            nodes_after=nodes_after, edges_after=edges_after)

    def test_unary_edges_no_overlap_internal_sample(self):
        nodes_before = """\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           1
        """
        edges_before = """\
        left    right   parent  child
        0       1       2       0
        1       2       2       1
        """
        self.verify_simplify(
            samples=[0, 1, 2],
            nodes_before=nodes_before, edges_before=edges_before,
            nodes_after=nodes_before, edges_after=edges_before)


class TestNonSampleExternalNodes(TopologyTestCase):
    """
    Tests for situations in which we have tips that are not samples.
    """
    def test_simple_case(self):
        # Simplest case where we have n = 2 and external non-sample nodes.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        3       0           0
        4       0           0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       2       0,1,3,4
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.1         0
        1   0.2         0
        2   0.3         0
        3   0.4         0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state
        0       0       1
        1       1       1
        2       3       1
        3       4       1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        self.assertEqual(ts.sample_size, 2)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 5)
        self.assertEqual(ts.num_sites, 4)
        self.assertEqual(ts.num_mutations, 4)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 2, 1: 2, 3: 2, 4: 2})
        self.assertEqual(t.root, 2)
        ts_simplified, node_map = ts.simplify(map_nodes=True)
        self.assertEqual(list(node_map), [0, 1, 2, -1, -1])
        self.assertEqual(ts_simplified.num_nodes, 3)
        self.assertEqual(ts_simplified.num_trees, 1)
        t = next(ts_simplified.trees())
        self.assertEqual(t.parent_dict, {0: 2, 1: 2})
        self.assertEqual(t.root, 2)
        # We should have removed the two non-sample mutations.
        self.assertEqual([s.position for s in t.sites()], [0.1, 0.2])

    def test_unary_non_sample_external_nodes(self):
        # Take an ordinary tree sequence and put a bunch of external non
        # sample nodes on it.
        ts = msprime.simulate(
            15, recombination_rate=5, random_seed=self.random_seed, mutation_rate=5)
        self.assertGreater(ts.num_trees, 2)
        self.assertGreater(ts.num_mutations, 2)
        tables = ts.dump_tables()
        next_node = ts.num_nodes
        tables.edges.reset()
        for e in ts.edges():
            tables.edges.add_row(e.left, e.right, e.parent, e.child)
            tables.edges.add_row(e.left, e.right, e.parent, next_node)
            tables.nodes.add_row(time=0)
            next_node += 1
        tables.sort()
        ts_new = tables.tree_sequence()
        self.assertEqual(ts_new.num_nodes, next_node)
        self.assertEqual(ts_new.sample_size, ts.sample_size)
        self.assert_haplotypes_equal(ts, ts_new)
        self.assert_variants_equal(ts, ts_new)
        ts_simplified = ts_new.simplify()
        self.assertEqual(ts_simplified.num_nodes, ts.num_nodes)
        self.assertEqual(ts_simplified.sample_size, ts.sample_size)
        self.assertEqual(list(ts_simplified.records()), list(ts.records()))
        self.assert_haplotypes_equal(ts, ts_simplified)
        self.assert_variants_equal(ts, ts_simplified)


class TestMultipleRoots(TopologyTestCase):
    """
    Tests for situations where we have multiple roots for the samples.
    """

    def test_simplest_degenerate_case(self):
        # Simplest case where we have n = 2 and no edges.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.1         0
        1   0.2         0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state
        0       0         1
        1       1         1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations,
            sequence_length=1, strict=False)
        self.assertEqual(ts.num_nodes, 2)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_sites, 2)
        self.assertEqual(ts.num_mutations, 2)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {})
        self.assertEqual(sorted(t.roots), [0, 1])
        self.assertEqual(list(ts.haplotypes()), ["10", "01"])
        self.assertEqual(
            [v.genotypes for v in ts.variants(as_bytes=True)], [b"10", b"01"])
        simplified = ts.simplify()
        t1 = ts.dump_tables()
        t2 = simplified.dump_tables()
        self.assertEqual(t1.nodes, t2.nodes)
        self.assertEqual(t1.edges, t2.edges)

    def test_simplest_non_degenerate_case(self):
        # Simplest case where we have n = 4 and two trees.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       1           0
        4       0           1
        5       0           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       4       0,1
        0       1       5       2,3
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.1         0
        1   0.2         0
        2   0.3         0
        3   0.4         0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state
        0       0       1
        1       1       1
        2       2       1
        3       3       1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        self.assertEqual(ts.num_nodes, 6)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_sites, 4)
        self.assertEqual(ts.num_mutations, 4)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 4, 2: 5, 3: 5})
        self.assertEqual(list(ts.haplotypes()), ["1000", "0100", "0010", "0001"])
        self.assertEqual(
            [v.genotypes for v in ts.variants(as_bytes=True)],
            [b"1000", b"0100", b"0010", b"0001"])
        self.assertEqual(t.mrca(0, 1), 4)
        self.assertEqual(t.mrca(0, 4), 4)
        self.assertEqual(t.mrca(2, 3), 5)
        self.assertEqual(t.mrca(0, 2), msprime.NULL_NODE)
        self.assertEqual(t.mrca(0, 3), msprime.NULL_NODE)
        self.assertEqual(t.mrca(2, 4), msprime.NULL_NODE)
        ts_simplified, node_map = ts.simplify(map_nodes=True)
        for j in range(4):
            self.assertEqual(node_map[j], j)
        self.assertEqual(ts_simplified.num_nodes, 6)
        self.assertEqual(ts_simplified.num_trees, 1)
        self.assertEqual(ts_simplified.num_sites, 4)
        self.assertEqual(ts_simplified.num_mutations, 4)
        t = next(ts_simplified.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 4, 2: 5, 3: 5})

    def test_two_reducable_trees(self):
        # We have n = 4 and two trees, with some unary nodes and non-sample leaves
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       1           0
        4       0           1
        5       0           1
        6       0           2
        7       0           3
        8       0           0   # Non sample leaf
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1      4         0
        0       1      5         1
        0       1      6         4,5
        0       1      7         2,3,8
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.1         0
        1   0.2         0
        2   0.3         0
        3   0.4         0
        4   0.5         0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state
        0       0       1
        1       1       1
        2       2       1
        3       3       1
        4       8       1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        self.assertEqual(ts.num_nodes, 9)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_sites, 5)
        self.assertEqual(ts.num_mutations, 5)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 5, 2: 7, 3: 7, 4: 6, 5: 6, 8: 7})
        self.assertEqual(list(ts.haplotypes()), ["10000", "01000", "00100", "00010"])
        self.assertEqual(
            [v.genotypes for v in ts.variants(as_bytes=True)],
            [b"1000", b"0100", b"0010", b"0001", b"0000"])
        self.assertEqual(t.mrca(0, 1), 6)
        self.assertEqual(t.mrca(2, 3), 7)
        self.assertEqual(t.mrca(2, 8), 7)
        self.assertEqual(t.mrca(0, 2), msprime.NULL_NODE)
        self.assertEqual(t.mrca(0, 3), msprime.NULL_NODE)
        self.assertEqual(t.mrca(0, 8), msprime.NULL_NODE)
        ts_simplified, node_map = ts.simplify(map_nodes=True)
        for j in range(4):
            self.assertEqual(node_map[j], j)
        self.assertEqual(ts_simplified.num_nodes, 6)
        self.assertEqual(ts_simplified.num_trees, 1)
        t = next(ts_simplified.trees())
        # print(ts_simplified.tables)
        self.assertEqual(
            list(ts_simplified.haplotypes()), ["1000", "0100", "0010", "0001"])
        self.assertEqual(
            [v.genotypes for v in ts_simplified.variants(as_bytes=True)],
            [b"1000", b"0100", b"0010", b"0001"])
        # The site over the non-sample external node should have been discarded.
        sites = list(t.sites())
        self.assertEqual(sites[-1].position, 0.4)
        self.assertEqual(t.parent_dict, {0: 4, 1: 4, 2: 5, 3: 5})

    def test_one_reducable_tree(self):
        # We have n = 4 and two trees. One tree is reducable and the other isn't.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       1           0
        4       0           1
        5       0           1
        6       0           2
        7       0           3
        8       0           0   # Non sample leaf
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1      4         0
        0       1      5         1
        0       1      6         4,5
        0       1      7         2,3,8
        """)
        ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
        self.assertEqual(ts.num_nodes, 9)
        self.assertEqual(ts.num_trees, 1)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 5, 2: 7, 3: 7, 4: 6, 5: 6, 8: 7})
        self.assertEqual(t.mrca(0, 1), 6)
        self.assertEqual(t.mrca(2, 3), 7)
        self.assertEqual(t.mrca(2, 8), 7)
        self.assertEqual(t.mrca(0, 2), msprime.NULL_NODE)
        self.assertEqual(t.mrca(0, 3), msprime.NULL_NODE)
        self.assertEqual(t.mrca(0, 8), msprime.NULL_NODE)
        ts_simplified = ts.simplify()
        self.assertEqual(ts_simplified.num_nodes, 6)
        self.assertEqual(ts_simplified.num_trees, 1)
        t = next(ts_simplified.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 4, 2: 5, 3: 5})

    # NOTE: This test has not been checked since updating to the text representation
    # so there might be other problems with it.
    def test_mutations_over_roots(self):
        # Mutations over root nodes should be ok when we have multiple roots.
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           0
        3       0           1
        4       0           2
        5       0           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       3       0,1
        0       1       4       3
        0       1       5       2
        """)
        sites = six.StringIO("""\
        id  position    ancestral_state
        0   0.1         0
        1   0.2         0
        2   0.3         0
        3   0.4         0
        4   0.5         0
        5   0.6         0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state
        0       0       1
        1       1       1
        2       3       1
        3       4       1
        4       2       1
        5       5       1
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        self.assertEqual(ts.num_nodes, 6)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_sites, 6)
        self.assertEqual(ts.num_mutations, 6)
        t = next(ts.trees())
        self.assertEqual(len(list(t.sites())), 6)
        haplotypes = ["101100", "011100", "000011"]
        variants = [b"100", b"010", b"110", b"110", b"001", b"001"]
        self.assertEqual(list(ts.haplotypes()), haplotypes)
        self.assertEqual([v.genotypes for v in ts.variants(as_bytes=True)], variants)
        ts_simplified = ts.simplify(filter_sites=False)
        self.assertEqual(list(ts_simplified.haplotypes()), haplotypes)
        self.assertEqual(
            [v.genotypes for v in ts_simplified.variants(as_bytes=True)], variants)

    def test_break_single_tree(self):
        # Take a single largish tree from msprime, and remove the oldest record.
        # This breaks it into two subtrees.
        ts = msprime.simulate(20, random_seed=self.random_seed, mutation_rate=4)
        self.assertGreater(ts.num_mutations, 5)
        tables = ts.dump_tables()
        tables.edges.set_columns(
            left=tables.edges.left[:-1],
            right=tables.edges.right[:-1],
            parent=tables.edges.parent[:-1],
            child=tables.edges.child[:-1])
        ts_new = msprime.load_tables(**tables.asdict())
        self.assertEqual(ts.sample_size, ts_new.sample_size)
        self.assertEqual(ts.num_edges, ts_new.num_edges + 1)
        self.assertEqual(ts.num_trees, ts_new.num_trees)
        self.assert_haplotypes_equal(ts, ts_new)
        self.assert_variants_equal(ts, ts_new)
        roots = set()
        t_new = next(ts_new.trees())
        for u in ts_new.samples():
            while t_new.parent(u) != msprime.NULL_NODE:
                u = t_new.parent(u)
            roots.add(u)
        self.assertEqual(len(roots), 2)
        self.assertEqual(sorted(roots), sorted(t_new.roots))


class TestWithVisuals(TopologyTestCase):
    """
    Some pedantic tests with ascii depictions of what's supposed to happen.
    """

    def verify_simplify_topology(self, ts, sample, haplotypes=False):
        # copies from test_highlevel.py
        new_ts, node_map = ts.simplify(sample, map_nodes=True)
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
                mapped_pair = [node_map[u] for u in pair]
                mrca1 = old_tree.get_mrca(*pair)
                mrca2 = new_tree.get_mrca(*mapped_pair)
                self.assertEqual(mrca2, node_map[mrca1])
        if haplotypes:
            orig_haps = list(ts.haplotypes())
            simp_haps = list(new_ts.haplotypes())
            for i, j in enumerate(sample):
                self.assertEqual(orig_haps[j], simp_haps[i])

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
        edges = six.StringIO("""\
        left    right   parent  child
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
        ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
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
        edges = six.StringIO("""\
        left    right   parent  child
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
        ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
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
        edges = six.StringIO("""\
        left    right   parent  child
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.0     0.2     6       5
        0.0     0.2     7       0,6
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
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
        # 1     4   |   5    |   4   *   5  |   4       5  |  4       5  |  4       5
        #       |\ / \ /|    |   |\   \     |   |\     /   |  |\     /   |  |\     /|
        # 2     | 6   7 |    |   | 6   7    |   | 6   7    |  | 6   7    |  | 6   7 |
        #       | |\ /| |    |   |  \  *    |   |  \  |    |  |  *       |  |  *    | ...
        # 3     | | 8 | |    |   |   8 |    |   *   8 *    |  |   8      |  |   8   |
        #       | |/ \| |    |   |  /  |    |   |  /  |    |  |  * *     |  |  / \  |
        # 4     | 9  10 |    |   | 9  10    |   | 9  10    |  | 9  10    |  | 9  10 |
        #       |/ \ / \|    |   |  \   *   |   |  \   \   |  |  \   *   |  |  \    |
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
        true_haplotypes = ['0100', '0001', '1110']
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
        edges = six.StringIO("""\
        left    right   parent  child
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
        sites = six.StringIO("""\
        position    ancestral_state
        0.05        0
        0.15        0
        0.25        0
        0.4         0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state   parent
        0       7       1               -1
        0      10       0               0
        0       2       1               1
        1       0       1               -1
        1      10       1               -1
        2       8       1               -1
        2       9       0               5
        2      10       0               5
        2       2       1               7
        3       8       1               -1
        """)
        ts = msprime.load_text(nodes, edges, sites, mutations, strict=False)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, len(true_trees))
        self.assertEqual(ts.num_nodes, 11)
        self.assertEqual(len(list(ts.edge_diffs())), ts.num_trees)
        # check topologies agree:
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        for j, x in enumerate(ts.haplotypes()):
            self.assertEqual(x, true_haplotypes[j])
        self.verify_simplify_topology(ts, [0, 1, 2], haplotypes=True)
        self.verify_simplify_topology(ts, [1, 0, 2], haplotypes=True)
        self.verify_simplify_topology(ts, [0, 1], haplotypes=False)
        self.verify_simplify_topology(ts, [1, 2], haplotypes=False)
        self.verify_simplify_topology(ts, [2, 0], haplotypes=False)

    def test_tricky_switches(self):
        # suppose the topology has:
        # left right parent child
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
        edges = six.StringIO("""\
        left right parent child
        0.0  0.5   6      0
        0.0  0.5   6      1
        0.5  1.0   6      4
        0.5  1.0   6      5
        0.0  0.4   7      2,3
        0.5  1.0   8      0
        0.5  1.0   8      1
        0.0  0.5   8      4
        0.0  0.5   8      5
        0.4  1.0   9      2,3
        0.4  1.0   10     8,9
        0.0  0.4   11     6,7
        0.4  1.0   12     6
        0.0  0.4   12     8
        0.4  1.0   12     10
        0.0  0.4   12     11
        """)
        true_trees = [
            {0: 6, 1: 6, 2: 7, 3: 7, 4: 8, 5: 8, 6: 11,
                7: 11, 8: 12, 9: -1, 10: -1, 11: 12, 12: -1},
            {0: 6, 1: 6, 2: 9, 3: 9, 4: 8, 5: 8, 6: 12,
                7: -1, 8: 10, 9: 10, 10: 12, 11: -1, 12: -1},
            {0: 8, 1: 8, 2: 9, 3: 9, 4: 6, 5: 6, 6: 12,
                7: -1, 8: 10, 9: 10, 10: 12, 11: -1, 12: -1}
        ]
        ts = msprime.load_text(nodes, edges, strict=False)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 6)
        self.assertEqual(ts.num_trees, len(true_trees))
        self.assertEqual(ts.num_nodes, 13)
        self.assertEqual(len(list(ts.edge_diffs())), ts.num_trees)
        # check topologies agree:
        for a, t in zip(true_trees, tree_dicts):
            for k in a.keys():
                if k in t.keys():
                    self.assertEqual(t[k], a[k])
                else:
                    self.assertEqual(a[k], msprime.NULL_NODE)
        self.verify_simplify_topology(ts, [0, 2])
        self.verify_simplify_topology(ts, [0, 4])
        self.verify_simplify_topology(ts, [2, 4])

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
        edges = six.StringIO("""\
        left right parent child
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
        big_ts = msprime.load_text(nodes, edges, strict=False)
        self.assertEqual(big_ts.num_trees, 1 + len(true_trees))
        self.assertEqual(big_ts.num_nodes, 16)
        ts, node_map = big_ts.simplify(map_nodes=True)
        self.assertEqual(list(node_map[:6]), list(range(6)))
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
        edges = six.StringIO("""\
        left    right   parent  child
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.2     0.8     8       3,5
        0.0     0.2     7       0,5
        """)
        first_ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
        ts, node_map = first_ts.simplify(map_nodes=True)
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6: -1, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: 8, 4: 5, 5: 8, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        # maps [1,2,3] -> [0,1,2]
        self.assertEqual(node_map[1], 0)
        self.assertEqual(node_map[2], 1)
        self.assertEqual(node_map[3], 2)
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
        self.assertEqual(first_ts.node(3).time, 0.2)
        self.assertEqual(ts.node(2).time, 0.2)
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
        edges = six.StringIO("""\
        left    right   parent  child
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.2     0.8     8       3,5
        0.0     0.2     7       0,5
        """)
        ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6: -1, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: 8, 4: 5, 5: 8, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(ts.num_nodes, 9)
        self.assertEqual(ts.node(0).time, 0.0)
        self.assertEqual(ts.node(1).time, 0.1)
        self.assertEqual(ts.node(2).time, 0.1)
        self.assertEqual(ts.node(3).time, 0.2)
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

    def test_internal_sampled_node(self):
        # 1.0             7
        # 0.7            / \                      8                     6
        #               /   \                    / \                   / \
        # 0.5          /     5                  /   5                 /   5
        #             /     /*\                /   /*\               /   /*\
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
        5       1           0.5
        6       0           0.7
        7       0           1.0
        8       0           0.8
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.0     0.2     4       2,3
        0.2     0.8     4       0,2
        0.8     1.0     4       2,3
        0.0     1.0     5       1,4
        0.8     1.0     6       0,5
        0.2     0.8     8       3,5
        0.0     0.2     7       0,5
        """)
        ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6: -1, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: 8, 4: 5, 5: 8, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        self.assertEqual(ts.sample_size, 4)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(ts.num_nodes, 9)
        self.assertEqual(ts.node(0).time, 0.0)
        self.assertEqual(ts.node(1).time, 0.1)
        self.assertEqual(ts.node(2).time, 0.1)
        self.assertEqual(ts.node(3).time, 0.2)
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
        self.check_num_samples(
            ts,
            [(0, 5, 4), (0, 2, 1), (0, 7, 4), (0, 4, 2),
             (1, 4, 1), (1, 5, 3), (1, 8, 4), (1, 0, 0),
             (2, 5, 4), (2, 1, 1)])
        self.check_num_tracked_samples(
            ts, [1, 2, 5],
            [(0, 5, 3), (0, 2, 1), (0, 7, 3), (0, 4, 1),
             (1, 4, 1), (1, 5, 3), (1, 8, 3), (1, 0, 0),
             (2, 5, 3), (2, 1, 1)])
        self.check_sample_iterator(
            ts,
            [(0, 0, []), (0, 5, [5, 1, 2, 3]), (0, 4, [2, 3]),
             (1, 5, [5, 1, 2]), (2, 4, [2, 3])])
        # pedantically check the SparseTree methods on the second tree
        tst = ts.trees()
        t = next(tst)
        t = next(tst)
        self.assertEqual(t.branch_length(1), 0.4)
        self.assertEqual(t.is_internal(0), False)
        self.assertEqual(t.is_leaf(0), True)
        self.assertEqual(t.is_sample(0), False)
        self.assertEqual(t.is_internal(1), False)
        self.assertEqual(t.is_leaf(1), True)
        self.assertEqual(t.is_sample(1), True)
        self.assertEqual(t.is_internal(5), True)
        self.assertEqual(t.is_leaf(5), False)
        self.assertEqual(t.is_sample(5), True)
        self.assertEqual(t.is_internal(4), True)
        self.assertEqual(t.is_leaf(4), False)
        self.assertEqual(t.is_sample(4), False)
        self.assertEqual(t.root, 8)
        self.assertEqual(t.mrca(0, 1), 5)
        self.assertEqual(t.sample_size, 4)


class TestBadTrees(unittest.TestCase):
    """
    Tests for bad tree sequence topologies that can only be detected when we
    try to create trees.
    """
    def test_simplest_contradictory_children(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        3       0           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.0     1.0     2       0
        0.0     1.0     3       0
        """)
        ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
        self.assertRaises(_msprime.LibraryError, list, ts.trees())

    def test_partial_overlap_contradictory_children(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        3       0           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.0     1.0     2       0,1
        0.5     1.0     3       0
        """)
        ts = msprime.load_text(nodes=nodes, edges=edges, strict=False)
        self.assertRaises(_msprime.LibraryError, list, ts.trees())


class TestSimplify(unittest.TestCase):
    """
    Tests that the implementations of simplify() do what they are supposed to.
    """
    random_seed = 23
    #
    #          8
    #         / \
    #        /   \
    #       /     \
    #      7       \
    #     / \       6
    #    /   5     / \
    #   /   / \   /   \
    #  4   0   1 2     3
    small_tree_ex_nodes = """\
    id      is_sample   population      time
    0       1       0               0.00000000000000
    1       1       0               0.00000000000000
    2       1       0               0.00000000000000
    3       1       0               0.00000000000000
    4       1       0               0.00000000000000
    5       0       0               0.14567111023387
    6       0       0               0.21385545626353
    7       0       0               0.43508024345063
    8       0       0               1.60156352971203
    """
    small_tree_ex_edges = """\
    id      left            right           parent  child
    0       0.00000000      1.00000000      5       0,1
    1       0.00000000      1.00000000      6       2,3
    2       0.00000000      1.00000000      7       4,5
    3       0.00000000      1.00000000      8       6,7
    """

    def do_simplify(
            self, ts, samples=None, compare_lib=True, filter_sites=True,
            filter_populations=True, filter_individuals=True):
        """
        Runs the Python test implementation of simplify.
        """
        if samples is None:
            samples = ts.samples()
        s = tests.Simplifier(
            ts, samples, filter_sites=filter_sites,
            filter_populations=filter_populations, filter_individuals=filter_individuals)
        new_ts, node_map = s.simplify()
        if compare_lib:
            sts, lib_node_map1 = ts.simplify(
                samples,
                filter_sites=filter_sites,
                filter_individuals=filter_individuals,
                filter_populations=filter_populations,
                map_nodes=True)
            lib_tables1 = sts.dump_tables()

            lib_tables2 = ts.dump_tables()
            lib_node_map2 = lib_tables2.simplify(
                samples,
                filter_sites=filter_sites,
                filter_individuals=filter_individuals,
                filter_populations=filter_populations)

            py_tables = new_ts.dump_tables()
            for lib_tables, lib_node_map in [
                    (lib_tables1, lib_node_map1), (lib_tables2, lib_node_map2)]:
                # print("lib = ")
                # print(lib_tables.nodes)
                # print(lib_tables.edges)
                # print("py = ")
                # print(py_tables.nodes)
                # print(py_tables.edges)

                self.assertEqual(lib_tables.nodes, py_tables.nodes)
                self.assertEqual(lib_tables.edges, py_tables.edges)
                self.assertEqual(lib_tables.migrations, py_tables.migrations)
                self.assertEqual(lib_tables.sites, py_tables.sites)
                self.assertEqual(lib_tables.mutations, py_tables.mutations)
                self.assertEqual(lib_tables.individuals, py_tables.individuals)
                self.assertEqual(lib_tables.populations, py_tables.populations)
                self.assertTrue(all(node_map == lib_node_map))
        return new_ts, node_map

    def verify_single_childified(self, ts):
        """
        Modify the specified tree sequence so that it has lots of unary
        nodes. Run simplify and verify we get the same tree sequence back.
        """
        ts_single = tsutil.single_childify(ts)
        tss, node_map = self.do_simplify(ts_single)
        # All original nodes should still be present.
        for u in range(ts.num_samples):
            self.assertEqual(u, node_map[u])
        # All introduced nodes should be mapped to null.
        for u in range(ts.num_samples, ts_single.num_samples):
            self.assertEqual(node_map[u], msprime.NULL_NODE)
        t1 = ts.dump_tables()
        t2 = tss.dump_tables()
        self.assertEqual(t1.nodes, t2.nodes)
        self.assertEqual(t1.edges, t2.edges)
        self.assertEqual(t1.sites, t2.sites)
        self.assertEqual(t1.mutations, t2.mutations)

    def verify_multiroot_internal_samples(self, ts):
        ts_multiroot = tsutil.decapitate(ts, ts.num_edges // 2)
        ts1 = tsutil.jiggle_samples(ts_multiroot)
        ts2, node_map = self.do_simplify(ts1)
        self.assertGreaterEqual(ts1.num_trees, ts2.num_trees)
        trees2 = ts2.trees()
        t2 = next(trees2)
        for t1 in ts1.trees():
            self.assertTrue(t2.interval[0] <= t1.interval[0])
            self.assertTrue(t2.interval[1] >= t1.interval[1])
            pairs = itertools.combinations(ts1.samples(), 2)
            for pair in pairs:
                mapped_pair = [node_map[u] for u in pair]
                mrca1 = t1.get_mrca(*pair)
                mrca2 = t2.get_mrca(*mapped_pair)
                if mrca1 == msprime.NULL_NODE:
                    assert mrca2 == msprime.NULL_NODE
                else:
                    self.assertEqual(node_map[mrca1], mrca2)
            if t2.interval[1] == t1.interval[1]:
                t2 = next(trees2, None)

    def test_single_tree(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        self.verify_single_childified(ts)
        self.verify_multiroot_internal_samples(ts)

    def test_single_tree_mutations(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=self.random_seed)
        self.assertGreater(ts.num_sites, 1)
        self.do_simplify(ts)
        self.verify_single_childified(ts)

    def test_many_trees_mutations(self):
        ts = msprime.simulate(
            10, recombination_rate=1, mutation_rate=10, random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        self.assertGreater(ts.num_sites, 2)
        self.do_simplify(ts)
        self.verify_single_childified(ts)

    def test_many_trees(self):
        ts = msprime.simulate(5, recombination_rate=4, random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        self.verify_single_childified(ts)
        self.verify_multiroot_internal_samples(ts)

    def test_small_tree_internal_samples(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        nodes = tables.nodes
        flags = nodes.flags
        # The parent of samples 0 and 1 is 5. Change this to an internal sample
        # and set 0 and 1 to be unsampled.
        flags[0] = 0
        flags[1] = 0
        flags[5] = msprime.NODE_IS_SAMPLE
        nodes.set_columns(flags=flags, time=nodes.time)
        ts = msprime.load_tables(nodes=nodes, edges=tables.edges)
        self.assertEqual(ts.sample_size, 4)
        tss, node_map = self.do_simplify(ts, [3, 5])
        self.assertEqual(node_map[3], 0)
        self.assertEqual(node_map[5], 1)
        self.assertEqual(tss.num_nodes, 3)
        self.assertEqual(tss.num_edges, 2)

    def test_small_tree_linear_samples(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        nodes = tables.nodes
        flags = nodes.flags
        # 7 is above 0. These are the only two samples
        flags[:] = 0
        flags[0] = msprime.NODE_IS_SAMPLE
        flags[7] = msprime.NODE_IS_SAMPLE
        nodes.set_columns(flags=flags, time=nodes.time)
        ts = msprime.load_tables(nodes=nodes, edges=tables.edges)
        self.assertEqual(ts.sample_size, 2)
        tss, node_map = self.do_simplify(ts, [0, 7])
        self.assertEqual(node_map[0], 0)
        self.assertEqual(node_map[7], 1)
        self.assertEqual(tss.num_nodes, 2)
        self.assertEqual(tss.num_edges, 1)
        t = next(tss.trees())
        self.assertEqual(t.parent_dict, {0: 1})

    def test_small_tree_internal_and_external_samples(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        nodes = tables.nodes
        flags = nodes.flags
        # 7 is above 0 and 1.
        flags[:] = 0
        flags[0] = msprime.NODE_IS_SAMPLE
        flags[1] = msprime.NODE_IS_SAMPLE
        flags[7] = msprime.NODE_IS_SAMPLE
        nodes.set_columns(flags=flags, time=nodes.time)
        ts = msprime.load_tables(nodes=nodes, edges=tables.edges)
        self.assertEqual(ts.sample_size, 3)
        tss, node_map = self.do_simplify(ts, [0, 1, 7])
        self.assertEqual(node_map[0], 0)
        self.assertEqual(node_map[1], 1)
        self.assertEqual(node_map[7], 2)
        self.assertEqual(tss.num_nodes, 4)
        self.assertEqual(tss.num_edges, 3)
        t = next(tss.trees())
        self.assertEqual(t.parent_dict, {0: 3, 1: 3, 3: 2})

    def test_small_tree_mutations(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        # Add some simple mutations here above the nodes we're keeping.
        tables.sites.add_row(position=0.25, ancestral_state="0")
        tables.sites.add_row(position=0.5, ancestral_state="0")
        tables.sites.add_row(position=0.75, ancestral_state="0")
        tables.sites.add_row(position=0.8, ancestral_state="0")
        tables.mutations.add_row(site=0, node=0, derived_state="1")
        tables.mutations.add_row(site=1, node=2, derived_state="1")
        tables.mutations.add_row(site=2, node=7, derived_state="1")
        tables.mutations.add_row(site=3, node=0, derived_state="1")
        ts = tables.tree_sequence()
        self.assertEqual(ts.num_sites, 4)
        self.assertEqual(ts.num_mutations, 4)
        tss = self.do_simplify(ts, [0, 2])[0]
        self.assertEqual(tss.sample_size, 2)
        self.assertEqual(tss.num_mutations, 4)
        self.assertEqual(list(tss.haplotypes()), ["1011", "0100"])

    def test_small_tree_filter_zero_mutations(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        ts = tsutil.insert_branch_sites(ts)
        self.assertEqual(ts.num_sites, 8)
        self.assertEqual(ts.num_mutations, 8)
        tss, _ = self.do_simplify(ts, [4, 0, 1], filter_sites=True)
        self.assertEqual(tss.num_sites, 5)
        self.assertEqual(tss.num_mutations, 5)
        tss, _ = self.do_simplify(ts, [4, 0, 1], filter_sites=False)
        self.assertEqual(tss.num_sites, 8)
        self.assertEqual(tss.num_mutations, 5)

    def test_small_tree_fixed_sites(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        # Add some simple mutations that will be fixed after simplify
        tables.sites.add_row(position=0.25, ancestral_state="0")
        tables.sites.add_row(position=0.5, ancestral_state="0")
        tables.sites.add_row(position=0.75, ancestral_state="0")
        tables.mutations.add_row(site=0, node=2, derived_state="1")
        tables.mutations.add_row(site=1, node=3, derived_state="1")
        tables.mutations.add_row(site=2, node=6, derived_state="1")
        ts = tables.tree_sequence()
        self.assertEqual(ts.num_sites, 3)
        self.assertEqual(ts.num_mutations, 3)
        tss, _ = self.do_simplify(ts, [4, 1])
        self.assertEqual(tss.sample_size, 2)
        self.assertEqual(tss.num_mutations, 0)
        self.assertEqual(list(tss.haplotypes()), ["", ""])

    def test_small_tree_mutations_over_root(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        tables.sites.add_row(position=0.25, ancestral_state="0")
        tables.mutations.add_row(site=0, node=8, derived_state="1")
        ts = tables.tree_sequence()
        self.assertEqual(ts.num_sites, 1)
        self.assertEqual(ts.num_mutations, 1)
        for filt in [True, False]:
            tss, _ = self.do_simplify(ts, [0, 1], filter_sites=filt)
            self.assertEqual(tss.num_sites, 1)
            self.assertEqual(tss.num_mutations, 1)

    def test_small_tree_recurrent_mutations(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        # Add recurrent mutation on the root branches
        tables.sites.add_row(position=0.25, ancestral_state="0")
        tables.mutations.add_row(site=0, node=6, derived_state="1")
        tables.mutations.add_row(site=0, node=7, derived_state="1")
        ts = tables.tree_sequence()
        self.assertEqual(ts.num_sites, 1)
        self.assertEqual(ts.num_mutations, 2)
        tss = self.do_simplify(ts, [4, 3])[0]
        self.assertEqual(tss.sample_size, 2)
        self.assertEqual(tss.num_sites, 1)
        self.assertEqual(tss.num_mutations, 2)
        self.assertEqual(list(tss.haplotypes()), ["1", "1"])

    def best_small_tree_back_mutations(self):
        ts = msprime.load_text(
            nodes=six.StringIO(self.small_tree_ex_nodes),
            edges=six.StringIO(self.small_tree_ex_edges), strict=False)
        tables = ts.dump_tables()
        # Add a chain of mutations
        tables.sites.add_row(position=0.25, ancestral_state="0")
        tables.mutations.add_row(site=0, node=7, derived_state="1")
        tables.mutations.add_row(site=0, node=5, derived_state="0")
        tables.mutations.add_row(site=0, node=1, derived_state="1")
        ts = tables.tree_sequence()
        self.assertEqual(ts.num_sites, 1)
        self.assertEqual(ts.num_mutations, 3)
        self.assertEqual(list(ts.haplotypes()), ["0", "1", "0", "0", "1"])
        # First check if we simplify for all samples and keep original state.
        tss = self.do_simplify(ts, [0, 1, 2, 3, 4])
        self.assertEqual(tss.sample_size, 5)
        self.assertEqual(tss.num_sites, 1)
        self.assertEqual(tss.num_mutations, 3)
        self.assertEqual(list(tss.haplotypes()), ["0", "1", "0", "0", "1"])

        # The ancestral state above 5 should be 0.
        tss = self.do_simplify(ts, [0, 1])
        self.assertEqual(tss.sample_size, 2)
        self.assertEqual(tss.num_sites, 1)
        self.assertEqual(tss.num_mutations, 1)
        self.assertEqual(list(tss.haplotypes()), ["0", "1"])

        # The ancestral state above 7 should be 1.
        tss = self.do_simplify(ts, [4, 0, 1])
        self.assertEqual(tss.sample_size, 3)
        self.assertEqual(tss.num_sites, 1)
        self.assertEqual(tss.num_mutations, 2)
        self.assertEqual(list(tss.haplotypes()), ["1", "0", "1"])

    def test_overlapping_unary_edges(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       2       2       0
        1       3       2       1
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        self.assertEqual(ts.sample_size, 2)
        self.assertEqual(ts.num_trees, 3)
        self.assertEqual(ts.sequence_length, 3)
        tss, node_map = self.do_simplify(ts)
        self.assertEqual(list(node_map), [0, 1, 2])
        trees = [{}, {0: 2, 1: 2}, {}]
        for t in tss.trees():
            self.assertEqual(t.parent_dict, trees[t.index])

    def test_overlapping_unary_edges_internal_samples(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           0
        2       1           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       2       2       0
        1       3       2       1
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        self.assertEqual(ts.sample_size, 3)
        self.assertEqual(ts.num_trees, 3)
        trees = [{0: 2}, {0: 2, 1: 2}, {1: 2}]
        for t in ts.trees():
            self.assertEqual(t.parent_dict, trees[t.index])
        tss, node_map = self.do_simplify(ts)
        self.assertEqual(list(node_map), [0, 1, 2])

    def test_isolated_samples(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       1           1
        2       1           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        """)
        ts = msprime.load_text(nodes, edges, sequence_length=1, strict=False)
        self.assertEqual(ts.num_samples, 3)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 3)
        tss, node_map = self.do_simplify(ts, compare_lib=True)
        self.assertEqual(ts.tables.nodes, tss.tables.nodes)
        self.assertEqual(ts.tables.edges, tss.tables.edges)
        self.assertEqual(list(node_map), [0, 1, 2])

    def test_internal_samples(self):
        nodes = six.StringIO("""\
        id      is_sample   population      time
        0       1       -1              1.00000000000000
        1       0       -1              1.00000000000000
        2       1       -1              1.00000000000000
        3       0       -1              1.31203521181726
        4       0       -1              2.26776380586006
        5       1       -1              0.00000000000000

        """)
        edges = six.StringIO("""\
        id      left            right           parent  child
        0       0.62185118      1.00000000      1       5
        1       0.00000000      0.62185118      2       5
        2       0.00000000      1.00000000      3       0,2
        3       0.00000000      1.00000000      4       1,3
        """)

        ts = msprime.load_text(nodes, edges, strict=False)
        tss, node_map = self.do_simplify(ts, [5, 2, 0], compare_lib=True)
        self.assertEqual(node_map[5], 0)
        self.assertEqual(node_map[2], 1)
        self.assertEqual(node_map[0], 2)
        self.assertEqual(node_map[1], -1)
        self.assertEqual(node_map[3], 3)
        self.assertEqual(node_map[4], 4)
        self.assertEqual(tss.sample_size, 3)
        self.assertEqual(tss.num_trees, 2)
        trees = [{0: 1, 1: 3, 2: 3}, {0: 4, 1: 3, 2: 3, 3: 4}]
        for t in tss.trees():
            self.assertEqual(t.parent_dict, trees[t.index])

    def test_many_mutations_over_single_sample_ancestral_state(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       1       0
        """)
        sites = six.StringIO("""\
        position    ancestral_state
        0           0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state   parent
        0       0       1               -1
        0       0       0               0
        """)
        ts = msprime.load_text(
            nodes, edges, sites=sites, mutations=mutations, strict=False)
        self.assertEqual(ts.sample_size, 1)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_sites, 1)
        self.assertEqual(ts.num_mutations, 2)
        tss, node_map = self.do_simplify(ts)
        self.assertEqual(tss.num_sites, 1)
        self.assertEqual(tss.num_mutations, 2)
        self.assertEqual(list(tss.haplotypes()), ["0"])

    def test_many_mutations_over_single_sample_derived_state(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       1           0
        1       0           1
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       1       0
        """)
        sites = six.StringIO("""\
        position    ancestral_state
        0           0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state   parent
        0       0       1               -1
        0       0       0               0
        0       0       1               1
        """)
        ts = msprime.load_text(
            nodes, edges, sites=sites, mutations=mutations, strict=False)
        self.assertEqual(ts.sample_size, 1)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_sites, 1)
        self.assertEqual(ts.num_mutations, 3)
        tss, node_map = self.do_simplify(ts)
        self.assertEqual(tss.num_sites, 1)
        self.assertEqual(tss.num_mutations, 3)
        self.assertEqual(list(tss.haplotypes()), ["1"])

    def test_many_trees_filter_zero_mutations(self):
        ts = msprime.simulate(5, recombination_rate=1, random_seed=10)
        self.assertGreater(ts.num_trees, 3)
        ts = tsutil.insert_branch_sites(ts)
        self.assertEqual(ts.num_sites, ts.num_mutations)
        self.assertGreater(ts.num_sites, ts.num_trees)
        for filter_sites in [True, False]:
            tss, _ = self.do_simplify(
                ts, samples=None, filter_sites=filter_sites)
            self.assertEqual(ts.num_sites, tss.num_sites)
            self.assertEqual(ts.num_mutations, tss.num_mutations)

    def test_many_trees_filter_zero_multichar_mutations(self):
        ts = msprime.simulate(5, recombination_rate=1, random_seed=10)
        self.assertGreater(ts.num_trees, 3)
        ts = tsutil.insert_multichar_mutations(ts)
        self.assertEqual(ts.num_sites, ts.num_trees)
        self.assertEqual(ts.num_mutations, ts.num_trees)
        for filter_sites in [True, False]:
            tss, _ = self.do_simplify(ts, samples=None, filter_sites=filter_sites)
            self.assertEqual(ts.num_sites, tss.num_sites)
            self.assertEqual(ts.num_mutations, tss.num_mutations)

    def test_simple_population_filter(self):
        ts = msprime.simulate(10, random_seed=2)
        tables = ts.dump_tables()
        tables.populations.add_row(metadata=b"unreferenced")
        self.assertEqual(len(tables.populations), 2)
        tss, _ = self.do_simplify(tables.tree_sequence(), filter_populations=True)
        self.assertEqual(tss.num_populations, 1)
        tss, _ = self.do_simplify(tables.tree_sequence(), filter_populations=False)
        self.assertEqual(tss.num_populations, 2)

    def test_interleaved_populations_filter(self):
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration(10),
                msprime.PopulationConfiguration(),
                msprime.PopulationConfiguration()],
            random_seed=2)
        self.assertEqual(ts.num_populations, 4)
        tables = ts.dump_tables()
        # Edit the populations so we can identify the rows.
        tables.populations.clear()
        for j in range(4):
            tables.populations.add_row(metadata=bytes([j]))
        ts = tables.tree_sequence()
        id_map = np.array([-1, 0, -1, -1], dtype=np.int32)
        tss, _ = self.do_simplify(ts, filter_populations=True)
        self.assertEqual(tss.num_populations, 1)
        population = tss.population(0)
        self.assertEqual(population.metadata, bytes([1]))
        self.assertTrue(np.array_equal(
            id_map[ts.tables.nodes.population], tss.tables.nodes.population))
        tss, _ = self.do_simplify(ts, filter_populations=False)
        self.assertEqual(tss.num_populations, 4)

    def test_removed_node_population_filter(self):
        tables = msprime.TableCollection(1)
        tables.populations.add_row(metadata=bytes(0))
        tables.populations.add_row(metadata=bytes(1))
        tables.populations.add_row(metadata=bytes(2))
        tables.nodes.add_row(flags=1, population=0)
        # Because flags=0 here, this node will be simplified out and the node
        # will disappear.
        tables.nodes.add_row(flags=0, population=1)
        tables.nodes.add_row(flags=1, population=2)
        tss, _ = self.do_simplify(tables.tree_sequence(), filter_populations=True)
        self.assertEqual(tss.num_nodes, 2)
        self.assertEqual(tss.num_populations, 2)
        self.assertEqual(tss.population(0).metadata, bytes(0))
        self.assertEqual(tss.population(1).metadata, bytes(2))
        self.assertEqual(tss.node(0).population, 0)
        self.assertEqual(tss.node(1).population, 1)

        tss, _ = self.do_simplify(
            tables.tree_sequence(), filter_populations=False)
        self.assertEqual(tss.tables.populations, tables.populations)

    def test_simple_individual_filter(self):
        tables = msprime.TableCollection(1)
        tables.individuals.add_row(flags=0)
        tables.individuals.add_row(flags=1)
        tables.nodes.add_row(flags=1, individual=0)
        tables.nodes.add_row(flags=1, individual=0)
        tss, _ = self.do_simplify(tables.tree_sequence(), filter_individuals=True)
        self.assertEqual(tss.num_nodes, 2)
        self.assertEqual(tss.num_individuals, 1)
        self.assertEqual(tss.individual(0).flags, 0)

        tss, _ = self.do_simplify(tables.tree_sequence(), filter_individuals=False)
        self.assertEqual(tss.tables.individuals, tables.individuals)

    def test_interleaved_individual_filter(self):
        tables = msprime.TableCollection(1)
        tables.individuals.add_row(flags=0)
        tables.individuals.add_row(flags=1)
        tables.individuals.add_row(flags=2)
        tables.nodes.add_row(flags=1, individual=1)
        tables.nodes.add_row(flags=1, individual=-1)
        tables.nodes.add_row(flags=1, individual=1)
        tss, _ = self.do_simplify(tables.tree_sequence(), filter_individuals=True)
        self.assertEqual(tss.num_nodes, 3)
        self.assertEqual(tss.num_individuals, 1)
        self.assertEqual(tss.individual(0).flags, 1)

        tss, _ = self.do_simplify(tables.tree_sequence(), filter_individuals=False)
        self.assertEqual(tss.tables.individuals, tables.individuals)

    def test_removed_node_individual_filter(self):
        tables = msprime.TableCollection(1)
        tables.individuals.add_row(flags=0)
        tables.individuals.add_row(flags=1)
        tables.individuals.add_row(flags=2)
        tables.nodes.add_row(flags=1, individual=0)
        # Because flags=0 here, this node will be simplified out and the node
        # will disappear.
        tables.nodes.add_row(flags=0, individual=1)
        tables.nodes.add_row(flags=1, individual=2)
        tss, _ = self.do_simplify(tables.tree_sequence(), filter_individuals=True)
        self.assertEqual(tss.num_nodes, 2)
        self.assertEqual(tss.num_individuals, 2)
        self.assertEqual(tss.individual(0).flags, 0)
        self.assertEqual(tss.individual(1).flags, 2)
        self.assertEqual(tss.node(0).individual, 0)
        self.assertEqual(tss.node(1).individual, 1)

        tss, _ = self.do_simplify(tables.tree_sequence(), filter_individuals=False)
        self.assertEqual(tss.tables.individuals, tables.individuals)

    def verify_simplify_haplotypes(self, ts, samples):
        sub_ts, node_map = self.do_simplify(ts, samples, filter_sites=False)
        self.assertEqual(ts.num_sites, sub_ts.num_sites)
        sub_haplotypes = list(sub_ts.haplotypes())
        all_samples = list(ts.samples())
        k = 0
        for j, h in enumerate(ts.haplotypes()):
            if k == len(samples):
                break
            if samples[k] == all_samples[j]:
                self.assertEqual(h, sub_haplotypes[k])
                k += 1

    def test_single_tree_recurrent_mutations(self):
        ts = msprime.simulate(6, random_seed=10)
        for mutations_per_branch in [1, 2, 3]:
            ts = tsutil.insert_branch_mutations(ts, mutations_per_branch)
            for num_samples in range(1, ts.num_samples):
                for samples in itertools.combinations(ts.samples(), num_samples):
                    self.verify_simplify_haplotypes(ts, samples)

    def test_many_trees_recurrent_mutations(self):
        ts = msprime.simulate(5, recombination_rate=1, random_seed=10)
        self.assertGreater(ts.num_trees, 3)
        for mutations_per_branch in [1, 2, 3]:
            ts = tsutil.insert_branch_mutations(ts, mutations_per_branch)
            for num_samples in range(1, ts.num_samples):
                for samples in itertools.combinations(ts.samples(), num_samples):
                    self.verify_simplify_haplotypes(ts, samples)

    def test_single_multiroot_tree_recurrent_mutations(self):
        ts = msprime.simulate(6, random_seed=10)
        ts = tsutil.decapitate(ts, ts.num_edges // 2)
        for mutations_per_branch in [1, 2, 3]:
            ts = tsutil.insert_branch_mutations(ts, mutations_per_branch)
            for num_samples in range(1, ts.num_samples):
                for samples in itertools.combinations(ts.samples(), num_samples):
                    self.verify_simplify_haplotypes(ts, samples)

    def test_many_multiroot_trees_recurrent_mutations(self):
        ts = msprime.simulate(7, recombination_rate=1, random_seed=10)
        self.assertGreater(ts.num_trees, 3)
        ts = tsutil.decapitate(ts, ts.num_edges // 2)
        for mutations_per_branch in [1, 2, 3]:
            ts = tsutil.insert_branch_mutations(ts, mutations_per_branch)
            for num_samples in range(1, ts.num_samples):
                for samples in itertools.combinations(ts.samples(), num_samples):
                    self.verify_simplify_haplotypes(ts, samples)

    def test_single_tree_recurrent_mutations_internal_samples(self):
        ts = msprime.simulate(6, random_seed=10)
        ts = tsutil.jiggle_samples(ts)
        for mutations_per_branch in [1, 2, 3]:
            ts = tsutil.insert_branch_mutations(ts, mutations_per_branch)
            for num_samples in range(1, ts.num_samples):
                for samples in itertools.combinations(ts.samples(), num_samples):
                    self.verify_simplify_haplotypes(ts, samples)

    def test_many_trees_recurrent_mutations_internal_samples(self):
        ts = msprime.simulate(5, recombination_rate=1, random_seed=10)
        ts = tsutil.jiggle_samples(ts)
        self.assertGreater(ts.num_trees, 3)
        for mutations_per_branch in [1, 2, 3]:
            ts = tsutil.insert_branch_mutations(ts, mutations_per_branch)
            for num_samples in range(1, ts.num_samples):
                for samples in itertools.combinations(ts.samples(), num_samples):
                    self.verify_simplify_haplotypes(ts, samples)


class TestMutationParent(unittest.TestCase):
    """
    Tests that mutation parent is correctly specified, and that we correctly
    recompute it with compute_mutation_parent.
    """
    seed = 42

    def verify_parents(self, ts):
        parent = tsutil.compute_mutation_parent(ts)
        tables = ts.tables
        self.assertTrue(np.array_equal(parent, tables.mutations.parent))
        mutations = tables.mutations
        mutations.set_columns(
            site=mutations.site, node=mutations.node,
            derived_state=mutations.derived_state,
            derived_state_offset=mutations.derived_state_offset)
        self.assertTrue(np.all(mutations.parent == msprime.NULL_MUTATION))
        tables.compute_mutation_parents()
        self.assertTrue(np.array_equal(parent, tables.mutations.parent))

    def test_example(self):
        nodes = six.StringIO("""\
        id      is_sample   time
        0       0           2.0
        1       0           1.0
        2       0           1.0
        3       1           0
        4       1           0
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0.0    0.5   2  3
        0.0    0.8   2  4
        0.5    1.0   1  3
        0.0    1.0   0  1
        0.0    1.0   0  2
        0.8    1.0   0  4
        """)
        sites = six.StringIO("""\
        position    ancestral_state
        0.1     0
        0.5     0
        0.9     0
        """)
        mutations = six.StringIO("""\
        site    node    derived_state   parent
        0       1       1               -1
        0       2       1               -1
        0       3       2               1
        1       0       1               -1
        1       1       1               3
        1       3       2               4
        1       2       1               3
        1       4       2               6
        2       0       1               -1
        2       1       1               8
        2       2       1               8
        2       4       1               8
        """)
        ts = msprime.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False)
        self.verify_parents(ts)

    def test_single_muts(self):
        ts = msprime.simulate(10, random_seed=self.seed, mutation_rate=3.0,
                              recombination_rate=1.0)
        self.verify_parents(ts)

    def test_with_jukes_cantor(self):
        ts = msprime.simulate(10, random_seed=self.seed, mutation_rate=0.0,
                              recombination_rate=1.0)
        # make *lots* of recurrent mutations
        mut_ts = tsutil.jukes_cantor(ts, num_sites=10, mu=1,
                                     multiple_per_node=False, seed=self.seed)
        self.verify_parents(mut_ts)

    def test_with_jukes_cantor_multiple_per_node(self):
        ts = msprime.simulate(10, random_seed=self.seed, mutation_rate=0.0,
                              recombination_rate=1.0)
        # make *lots* of recurrent mutations
        mut_ts = tsutil.jukes_cantor(ts, num_sites=10, mu=1,
                                     multiple_per_node=True, seed=self.seed)
        self.verify_parents(mut_ts)

    def verify_branch_mutations(self, ts, mutations_per_branch):
        ts = tsutil.insert_branch_mutations(ts, mutations_per_branch)
        self.assertGreater(ts.num_mutations, 1)
        self.verify_parents(ts)

    def test_single_tree_one_mutation_per_branch(self):
        ts = msprime.simulate(6, random_seed=10)
        self.verify_branch_mutations(ts, 1)

    def test_single_tree_two_mutations_per_branch(self):
        ts = msprime.simulate(10, random_seed=9)
        self.verify_branch_mutations(ts, 2)

    def test_single_tree_three_mutations_per_branch(self):
        ts = msprime.simulate(8, random_seed=9)
        self.verify_branch_mutations(ts, 3)

    def test_single_multiroot_tree_recurrent_mutations(self):
        ts = msprime.simulate(6, random_seed=10)
        ts = tsutil.decapitate(ts, ts.num_edges // 2)
        for mutations_per_branch in [1, 2, 3]:
            self.verify_branch_mutations(ts, mutations_per_branch)

    def test_many_multiroot_trees_recurrent_mutations(self):
        ts = msprime.simulate(7, recombination_rate=1, random_seed=10)
        self.assertGreater(ts.num_trees, 3)
        ts = tsutil.decapitate(ts, ts.num_edges // 2)
        for mutations_per_branch in [1, 2, 3]:
            self.verify_branch_mutations(ts, mutations_per_branch)


class TestSimpleTreeAlgorithm(unittest.TestCase):
    """
    Tests for the direct implementation of Algorithm T in tsutil.py.

    See TestHoleyTreeSequences above for further tests on wacky topologies.
    """
    def test_zero_nodes(self):
        nodes = msprime.NodeTable()
        edges = msprime.EdgeTable()
        ts = msprime.load_tables(nodes=nodes, edges=edges, sequence_length=1)
        self.assertEqual(ts.sequence_length, 1)
        self.assertEqual(ts.num_trees, 1)
        # Test the simple tree iterator.
        trees = list(tsutil.algorithm_T(ts))
        self.assertEqual(len(trees), 1)
        (left, right), parent = trees[0]
        self.assertEqual(left, 0)
        self.assertEqual(right, 1)
        self.assertEqual(parent, [])

    def test_one_node(self):
        nodes = msprime.NodeTable()
        edges = msprime.EdgeTable()
        nodes.add_row()
        ts = msprime.load_tables(nodes=nodes, edges=edges, sequence_length=1)
        self.assertEqual(ts.sequence_length, 1)
        self.assertEqual(ts.num_trees, 1)
        # Test the simple tree iterator.
        trees = list(tsutil.algorithm_T(ts))
        self.assertEqual(len(trees), 1)
        (left, right), parent = trees[0]
        self.assertEqual(left, 0)
        self.assertEqual(right, 1)
        self.assertEqual(parent, [-1])

    def test_single_coalescent_tree(self):
        ts = msprime.simulate(10, random_seed=1, length=10)
        tree = ts.first()
        p1 = [tree.parent(j) for j in range(ts.num_nodes)]
        interval, p2 = next(tsutil.algorithm_T(ts))
        self.assertEqual(interval, tree.interval)
        self.assertEqual(p1, p2)

    def test_coalescent_trees(self):
        ts = msprime.simulate(8, recombination_rate=5, random_seed=1, length=2)
        self.assertGreater(ts.num_trees, 2)
        new_trees = tsutil.algorithm_T(ts)
        for tree in ts.trees():
            interval, p2 = next(new_trees)
            p1 = [tree.parent(j) for j in range(ts.num_nodes)]
            self.assertEqual(interval, tree.interval)
            self.assertEqual(p1, p2)
        self.assertRaises(StopIteration, next, new_trees)


class TestSampleLists(unittest.TestCase):
    """
    Tests for the sample lists algorithm.
    """
    def verify(self, ts):
        tree1 = tsutil.LinkedTree(ts)
        s = str(tree1)
        self.assertIsNotNone(s)
        trees = ts.trees(sample_lists=True)
        for left, right in tree1.sample_lists():
            tree2 = next(trees)
            assert (left, right) == tree2.interval
            for u in tree2.nodes():
                self.assertEqual(tree1.left_sample[u], tree2.left_sample(u))
                self.assertEqual(tree1.right_sample[u], tree2.right_sample(u))
            for j in range(ts.num_samples):
                self.assertEqual(tree1.next_sample[j], tree2.next_sample(j))
        assert right == ts.sequence_length

        tree1 = tsutil.LinkedTree(ts)
        trees = ts.trees(sample_lists=False)
        sample_index_map = ts.samples()
        for left, right in tree1.sample_lists():
            tree2 = next(trees)
            for u in range(ts.num_nodes):
                samples2 = list(tree2.samples(u))
                samples1 = []
                index = tree1.left_sample[u]
                if index != msprime.NULL:
                    self.assertEqual(
                        sample_index_map[tree1.left_sample[u]], samples2[0])
                    self.assertEqual(
                        sample_index_map[tree1.right_sample[u]], samples2[-1])
                    stop = tree1.right_sample[u]
                    while True:
                        assert index != -1
                        samples1.append(sample_index_map[index])
                        if index == stop:
                            break
                        index = tree1.next_sample[index]
                self.assertEqual(samples1, samples2)
        assert right == ts.sequence_length

    def test_single_coalescent_tree(self):
        ts = msprime.simulate(10, random_seed=1, length=10)
        self.verify(ts)

    def test_coalescent_trees(self):
        ts = msprime.simulate(8, recombination_rate=5, random_seed=1, length=2)
        self.assertGreater(ts.num_trees, 2)
        self.verify(ts)

    def test_coalescent_trees_internal_samples(self):
        ts = msprime.simulate(8, recombination_rate=5, random_seed=10, length=2)
        self.assertGreater(ts.num_trees, 2)
        self.verify(tsutil.jiggle_samples(ts))

    def test_coalescent_trees_all_samples(self):
        ts = msprime.simulate(8, recombination_rate=5, random_seed=10, length=2)
        self.assertGreater(ts.num_trees, 2)
        tables = ts.dump_tables()
        tables.nodes.set_columns(
            flags=np.zeros_like(tables.nodes.flags) + msprime.NODE_IS_SAMPLE,
            time=tables.nodes.time)
        self.verify(tables.tree_sequence())

    def test_wright_fisher_trees_unsimplified(self):
        tables = wf.wf_sim(10, 5, deep_history=False, seed=2)
        tables.sort()
        ts = tables.tree_sequence()
        self.verify(ts)

    def test_wright_fisher_trees_simplified(self):
        tables = wf.wf_sim(10, 5, deep_history=False, seed=1)
        tables.sort()
        ts = tables.tree_sequence()
        ts = ts.simplify()
        self.verify(ts)

    def test_wright_fisher_trees_simplified_one_gen(self):
        tables = wf.wf_sim(10, 1, deep_history=False, seed=1)
        tables.sort()
        ts = tables.tree_sequence()
        ts = ts.simplify()
        self.verify(ts)

    def test_nonbinary_trees(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=1.0, population=0, proportion=0.95)]
        ts = msprime.simulate(
            20, recombination_rate=10, mutation_rate=5,
            demographic_events=demographic_events, random_seed=7)
        found = False
        for e in ts.edgesets():
            if len(e.children) > 2:
                found = True
        self.assertTrue(found)
        self.verify(ts)

    def test_many_multiroot_trees(self):
        ts = msprime.simulate(7, recombination_rate=1, random_seed=10)
        self.assertGreater(ts.num_trees, 3)
        ts = tsutil.decapitate(ts, ts.num_edges // 2)
        self.verify(ts)


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


def reduce_topology(ts):
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
    return ts.simplify(map_nodes=True, filter_sites=False)


class TestReduceTopology(unittest.TestCase):
    """
    Tests to ensure that reduce topology in simplify is equivalent to the
    reduce_topology function above.
    """

    def verify(self, ts):
        source_tables = ts.tables
        X = source_tables.sites.position
        position_count = {x: 0 for x in X}
        position_count[0] = 0
        position_count[ts.sequence_length] = 0
        mts, node_map = reduce_topology(ts)
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
        s = tests.Simplifier(
            ts, ts.samples(), reduce_to_site_topology=True, filter_sites=False)
        sts1, _ = s.simplify()
        sts2 = ts.simplify(reduce_to_site_topology=True, filter_sites=False)
        t1 = mts.tables
        for sts in [sts2, sts2]:
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

    def test_zero_sites(self):
        ts = msprime.simulate(5, random_seed=2)
        self.assertEqual(ts.num_sites, 0)
        mts = ts.simplify(reduce_to_site_topology=True)
        self.assertEqual(mts.num_trees, 1)
        self.assertEqual(mts.num_edges, 0)

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


def search_sorted(a, v):
    """
    Implementation of searchsorted based on binary search with the same
    semantics as numpy's searchsorted. Used as the basis of the C
    implementation which we use in the simplify algorithm.
    """
    upper = len(a)
    if upper == 0:
        return 0
    lower = 0
    while upper - lower > 1:
        mid = (upper + lower) // 2
        if (v >= a[mid]):
            lower = mid
        else:
            upper = mid
    offset = 0
    if a[lower] < v:
        offset = 1
    return lower + offset


class TestSearchSorted(unittest.TestCase):
    """
    Tests for the basic implementation of search_sorted.
    """
    def verify(self, a):
        a = np.array(a)
        start, end = a[0], a[-1]
        # Check random values.
        np.random.seed(43)
        for v in np.random.uniform(start, end, 10):
            self.assertEqual(search_sorted(a, v), np.searchsorted(a, v))
        # Check equal values.
        for v in a:
            self.assertEqual(search_sorted(a, v), np.searchsorted(a, v))
        # Check values outside bounds.
        for v in [start - 2, start - 1, end, end + 1, end + 2]:
            self.assertEqual(search_sorted(a, v), np.searchsorted(a, v))

    def test_range(self):
        for j in range(1, 20):
            self.verify(range(j))

    def test_negative_range(self):
        for j in range(1, 20):
            self.verify(-1 * np.arange(j)[::-1])

    def test_random_unit_interval(self):
        np.random.seed(143)
        for size in range(1, 100):
            a = np.random.random(size=size)
            a.sort()
            self.verify(a)

    def test_random_interval(self):
        np.random.seed(143)
        for _ in range(10):
            interval = np.random.random(2) * 10
            interval.sort()
            a = np.random.uniform(*interval, size=100)
            a.sort()
            self.verify(a)

    def test_random_negative(self):
        np.random.seed(143)
        for _ in range(10):
            interval = np.random.random(2) * 5
            interval.sort()
            a = -1 * np.random.uniform(*interval, size=100)
            a.sort()
            self.verify(a)

    def test_edge_cases(self):
        for v in [0, 1]:
            self.assertEqual(search_sorted([], v), np.searchsorted([], v))
            self.assertEqual(search_sorted([1], v), np.searchsorted([1], v))
