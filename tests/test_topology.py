#
# Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

import unittest
import itertools

import msprime
import _msprime


def cr(left=0, right=1, node=None, children=None, time=None, population=0):
    """
    Convenience short-hand function for defining records.
    """
    return msprime.CoalescenceRecord(
        left=left, right=right, node=node, children=children, time=time,
        population=population)


def build_tree_sequence(records, mutations=[]):
    ts = _msprime.TreeSequence()
    # FIXME this is a hack to get the sample size. This API will be removed in favour
    # of specifying nodes, so no point in cleaning it up.
    n = max(2, min(r[2] for r in records))
    samples = [(0, 0) for _ in range(n)]
    ts.load_records(samples, records, mutations)
    return msprime.TreeSequence(ts)


def insert_redundant_breakpoints(ts):
    """
    Builds a new tree sequence containing redundant breakpoints.
    """
    new_records = []
    for r in ts.records():
        x = r.left + (r.right - r.left) / 2
        new_records.append(cr(
            left=r.left, right=x, children=r.children, node=r.node,
            population=r.population, time=r.time))
        new_records.append(cr(
            left=x, right=r.right, children=r.children, node=r.node,
            population=r.population, time=r.time))
    new_ts = build_tree_sequence(new_records)
    assert new_ts.num_records == 2 * ts.num_records
    return new_ts


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


class TestRecordSquashing(TopologyTestCase):
    """
    Tests that we correctly squash adjacent equal records together.
    """
    def test_single_record(self):
        records = [
            cr(left=0, right=1, node=2, children=(0, 1), time=1),
            cr(left=1, right=2, node=2, children=(0, 1), time=1),
        ]
        ts = build_tree_sequence(records)
        self.assertEqual(list(ts.records()), records)
        tss = ts.simplify()
        simplified_records = list(tss.records())
        self.assertEqual(len(simplified_records), 1)
        r = simplified_records[0]
        self.assertEqual(r.left, 0)
        self.assertEqual(r.right, 2)

    def test_single_tree(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        ts_redundant = insert_redundant_breakpoints(ts)
        tss = ts_redundant.simplify()
        self.assertEqual(list(tss.records()), list(ts.records()))

    def test_many_trees(self):
        ts = msprime.simulate(
                20, recombination_rate=5,
                random_seed=self.random_seed)
        self.assertGreater(ts.num_trees, 2)
        ts_redundant = insert_redundant_breakpoints(ts)
        tss = ts_redundant.simplify()
        self.assertEqual(list(tss.records()), list(ts.records()))


class TestRedundantBreakpoints(TopologyTestCase):
    """
    Tests for dealing with redundant breakpoints within the tree sequence.
    These are records that may be squashed together into a single record.
    """
    def test_single_tree(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        ts_redundant = insert_redundant_breakpoints(ts)
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
        ts_redundant = insert_redundant_breakpoints(ts)
        self.assertEqual(ts.sample_size, ts_redundant.sample_size)
        self.assertEqual(ts.sequence_length, ts_redundant.sequence_length)
        self.assertGreater(ts_redundant.num_trees, ts.num_trees)
        self.assertGreater(ts_redundant.num_records, ts.num_records)
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
        records = [
            cr(left=0, right=1, node=2, children=(0,), time=1),
            cr(left=0, right=1, node=3, children=(1,), time=1),
            cr(left=0, right=1, node=4, children=(2, 3), time=2),
            cr(left=0, right=1, node=5, children=(4,), time=3),
        ]
        mutations = [(j * 1 / 5, (j,)) for j in range(5)]
        ts = build_tree_sequence(records, mutations)
        self.assertEqual(ts.sample_size, 2)
        self.assertEqual(ts.num_nodes, 6)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(len(list(ts.diffs())), ts.num_trees)
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
        num_unary_nodes = 50
        n = 2
        records = [
            cr(left=0, right=1, node=n, children=(0,), time=1)]
        for j in range(num_unary_nodes):
            records.append(cr(node=n + j + 1, children=(n + j,), time=j + 2))
        root = num_unary_nodes + 3
        root_time = num_unary_nodes + 3
        records.append(cr(node=root, children=(1, num_unary_nodes + 2), time=root_time))
        ts = build_tree_sequence(records)
        t = next(ts.trees())
        self.assertEqual(t.mrca(0, 1), root)
        self.assertEqual(t.tmrca(0, 1), root_time)
        ts_simplified = ts.simplify()
        self.assertEqual(ts_simplified.num_records, 1)
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
        new_records = []
        next_node = ts.num_nodes
        for r in ts.records():
            u = r.node
            t = r.time - 1e-14  # Arbitrary small value.
            children = []
            for v in r.children:
                new_records.append(cr(
                    left=r.left, right=r.right, population=r.population,
                    node=next_node, children=(v,), time=t))
                children.append(next_node)
                next_node += 1
            new_records.append(cr(
                left=r.left, right=r.right, population=r.population,
                node=u, children=tuple(children), time=r.time))
        new_records.sort(key=lambda r: r.time)
        ts_new = build_tree_sequence(new_records, list(ts.mutations()))
        self.assertGreater(ts_new.num_records, ts.num_records)
        self.assert_haplotypes_equal(ts, ts_new)
        self.assert_variants_equal(ts, ts_new)
        ts_simplified = ts_new.simplify()
        self.assertEqual(list(ts_simplified.records()), list(ts.records()))
        self.assert_haplotypes_equal(ts, ts_simplified)
        self.assert_variants_equal(ts, ts_simplified)
        self.assertEqual(len(list(ts.diffs())), ts.num_trees)

    def test_binary_tree_sequence_unary_nodes(self):
        ts = msprime.simulate(
            20, recombination_rate=5, mutation_rate=5, random_seed=self.random_seed)
        self.verify_unary_tree_sequence(ts)

    def test_nonbinary_tree_sequence_unary_nodes(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=1.0, proportion=0.95)]
        ts = msprime.simulate(
            20, recombination_rate=10, mutation_rate=5,
            demographic_events=demographic_events, random_seed=self.random_seed)
        found = False
        for r in ts.records():
            if len(r.children) > 2:
                found = True
        self.assertTrue(found)
        self.verify_unary_tree_sequence(ts)


class TestNonSampleExternalNodes(TopologyTestCase):
    """
    Tests for situations in which we have tips that are not samples.
    """
    def test_simple_case(self):
        # Simplest case where we have n = 2 and external non-sample nodes.
        records = [cr(node=2, children=(0, 1, 3, 4), time=1)]
        mutations = [
            msprime.Mutation(index=0, position=0.1, nodes=(0,)),
            msprime.Mutation(index=1, position=0.2, nodes=(1,)),
            msprime.Mutation(index=2, position=0.3, nodes=(3,)),
            msprime.Mutation(index=3, position=0.4, nodes=(4,)),
        ]
        ts = build_tree_sequence(records, mutations)
        self.assertEqual(ts.sample_size, 2)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(ts.num_nodes, 5)
        t = next(ts.trees())
        self.assertEqual(list(t.mutations()), mutations)
        self.assertEqual(t.parent_dict, {0: 2, 1: 2, 3: 2, 4: 2})
        self.assertEqual(t.time_dict, {0: 0, 1: 0, 3: 0, 4: 0, 2: 1})
        self.assertEqual(t.root, 2)
        ts_simplified = ts.simplify()
        self.assertEqual(ts_simplified.num_nodes, 3)
        self.assertEqual(ts_simplified.num_trees, 1)
        t = next(ts_simplified.trees())
        self.assertEqual(t.parent_dict, {0: 2, 1: 2})
        self.assertEqual(t.time_dict, {0: 0, 1: 0, 2: 1})
        self.assertEqual(t.root, 2)
        # We should have removed the two non-sample mutations.
        self.assertEqual(list(t.mutations()), mutations[:-2])

    def test_unary_non_sample_external_nodes(self):
        # Take an ordinary tree sequence and put a bunch of external non
        # sample nodes on it.
        ts = msprime.simulate(
            15, recombination_rate=5, random_seed=self.random_seed, mutation_rate=5)
        self.assertGreater(ts.num_trees, 2)
        self.assertGreater(ts.num_mutations, 2)
        new_records = []
        next_node = ts.num_nodes
        for r in ts.records():
            children = tuple(list(r.children) + [next_node])
            new_records.append(cr(
                left=r.left, right=r.right, node=r.node, time=r.time,
                population=r.population, children=children))
            next_node += 1
        ts_new = build_tree_sequence(new_records, list(ts.mutations()))
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
        # Simplest case where we have n = 2 and two unary records.
        # This cannot be simplified, since there are no trees to recover.
        records = [
            cr(left=0, right=1, node=2, children=(0,), time=1),
            cr(left=0, right=1, node=3, children=(1,), time=1),
        ]
        mutations = [
            msprime.Mutation(index=0, position=0.1, nodes=(0,)),
            msprime.Mutation(index=1, position=0.2, nodes=(1,)),
        ]
        ts = build_tree_sequence(records, mutations)
        self.assertEqual(ts.num_nodes, 4)
        self.assertEqual(ts.num_trees, 1)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 2, 1: 3})
        self.assertEqual(t.time_dict, {0: 0, 1: 0, 2: 1, 3: 1})
        self.assertEqual(list(t.mutations()), mutations)
        self.assertEqual(list(ts.haplotypes()), ["10", "01"])
        self.assertEqual(
            [v.genotypes for v in ts.variants(as_bytes=True)], [b"10", b"01"])
        self.assertRaises(_msprime.LibraryError, ts.simplify)

    def test_simplest_non_degenerate_case(self):
        # Simplest case where we have n = 4 and two trees.
        records = [
            cr(left=0, right=1, node=4, children=(0, 1), time=1),
            cr(left=0, right=1, node=5, children=(2, 3), time=2),
        ]
        mutations = [
            msprime.Mutation(index=0, position=0.1, nodes=(0,)),
            msprime.Mutation(index=1, position=0.2, nodes=(1,)),
            msprime.Mutation(index=2, position=0.3, nodes=(2,)),
            msprime.Mutation(index=3, position=0.4, nodes=(3,)),
        ]
        ts = build_tree_sequence(records, mutations)
        self.assertEqual(ts.num_nodes, 6)
        self.assertEqual(ts.num_trees, 1)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 4, 2: 5, 3: 5})
        self.assertEqual(t.time_dict, {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 2})
        self.assertEqual(list(t.mutations()), mutations)
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
        ts_simplified = ts.simplify()
        self.assertEqual(ts_simplified.num_nodes, 6)
        self.assertEqual(ts_simplified.num_trees, 1)
        t = next(ts_simplified.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 4, 2: 5, 3: 5})
        self.assertEqual(t.time_dict, {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 2})
        self.assertEqual(list(t.mutations()), mutations)

    def test_two_reducable_trees(self):
        # We have n = 4 and two trees, with some unary nodes and non-sample leaves
        records = [
            cr(left=0, right=1, node=4, children=(0,), time=1),
            cr(left=0, right=1, node=5, children=(1,), time=1),
            cr(left=0, right=1, node=6, children=(4, 5), time=2),
            cr(left=0, right=1, node=7, children=(2, 3, 8), time=3),
        ]
        mutations = [
            msprime.Mutation(index=0, position=0.1, nodes=(0,)),
            msprime.Mutation(index=1, position=0.2, nodes=(1,)),
            msprime.Mutation(index=2, position=0.3, nodes=(2,)),
            msprime.Mutation(index=3, position=0.4, nodes=(3,)),
            msprime.Mutation(index=4, position=0.5, nodes=(8,)),
        ]
        ts = build_tree_sequence(records, mutations)
        self.assertEqual(ts.num_nodes, 9)
        self.assertEqual(ts.num_trees, 1)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 5, 2: 7, 3: 7, 4: 6, 5: 6, 8: 7})
        self.assertEqual(
            t.time_dict, {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 1, 6: 2, 7: 3, 8: 0})
        self.assertEqual(list(t.mutations()), mutations)
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
        ts_simplified = ts.simplify()
        self.assertEqual(ts_simplified.num_nodes, 6)
        self.assertEqual(ts_simplified.num_trees, 1)
        t = next(ts_simplified.trees())
        self.assertEqual(
            list(ts_simplified.haplotypes()), ["1000", "0100", "0010", "0001"])
        self.assertEqual(
            [v.genotypes for v in ts_simplified.variants(as_bytes=True)],
            [b"1000", b"0100", b"0010", b"0001"])
        # The mutation over the non-sample external node should have been discarded.
        self.assertEqual(list(t.mutations()), mutations[:-1])
        self.assertEqual(t.parent_dict, {0: 4, 1: 4, 2: 5, 3: 5})
        self.assertEqual(t.time_dict, {0: 0, 1: 0, 2: 0, 3: 0, 4: 2, 5: 3})

    def test_one_reducable_tree(self):
        # We have n = 3 and two trees. One tree is reducable and the other isn't.
        records = [
            cr(left=0, right=1, node=4, children=(0,), time=1),
            cr(left=0, right=1, node=5, children=(1,), time=1),
            cr(left=0, right=1, node=6, children=(4, 5), time=2),
            cr(left=0, right=1, node=7, children=(2, 3, 8), time=3),
        ]
        ts = build_tree_sequence(records)
        self.assertEqual(ts.num_nodes, 9)
        self.assertEqual(ts.num_trees, 1)
        t = next(ts.trees())
        self.assertEqual(t.parent_dict, {0: 4, 1: 5, 2: 7, 3: 7, 4: 6, 5: 6, 8: 7})
        self.assertEqual(
            t.time_dict, {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 1, 6: 2, 7: 3, 8: 0})
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
        self.assertEqual(t.time_dict, {0: 0, 1: 0, 2: 0, 3: 0, 4: 2, 5: 3})

    def test_mutations_over_roots(self):
        # Mutations over root nodes should be ok when we have multiple roots.
        records = [
            cr(left=0, right=1, node=3, children=(0, 1), time=1),
            cr(left=0, right=1, node=4, children=(3,), time=2),
            cr(left=0, right=1, node=5, children=(2,), time=2),
        ]
        mutations = [
            msprime.Mutation(index=0, position=0.1, nodes=(0,)),
            msprime.Mutation(index=1, position=0.2, nodes=(1,)),
            msprime.Mutation(index=2, position=0.3, nodes=(3,)),
            msprime.Mutation(index=3, position=0.4, nodes=(4,)),
            msprime.Mutation(index=4, position=0.5, nodes=(2,)),
            msprime.Mutation(index=5, position=0.6, nodes=(5,)),
        ]
        ts = build_tree_sequence(records, mutations)
        self.assertEqual(ts.num_nodes, 6)
        self.assertEqual(ts.num_trees, 1)
        self.assertEqual(list(ts.mutations()), mutations)
        t = next(ts.trees())
        self.assertEqual(list(t.mutations()), mutations)
        haplotypes = ["101100", "011100", "000011"]
        variants = [b"100", b"010", b"110", b"110", b"001", b"001"]
        self.assertEqual(list(ts.haplotypes()), haplotypes)
        self.assertEqual([v.genotypes for v in ts.variants(as_bytes=True)], variants)
        ts_simplified = ts.simplify(filter_root_mutations=False)
        self.assertEqual(list(ts_simplified.haplotypes()), haplotypes)
        self.assertEqual(
            [v.genotypes for v in ts_simplified.variants(as_bytes=True)], variants)
        ts_simplified = ts.simplify(filter_root_mutations=True)
        self.assertEqual(list(ts_simplified.haplotypes()), ["10", "01", "00"])
        self.assertEqual(
            [v.genotypes for v in ts_simplified.variants(as_bytes=True)],
            [b"100", b"010"])

    def test_break_single_tree(self):
        # Take a single largish tree from msprime, and remove the oldest record.
        # This breaks it into two subtrees.
        ts = msprime.simulate(20, random_seed=self.random_seed, mutation_rate=4)
        self.assertGreater(ts.num_mutations, 5)
        records = list(ts.records())
        ts_new = build_tree_sequence(records[:-1], list(ts.mutations()))
        self.assertEqual(ts.sample_size, ts_new.sample_size)
        self.assertEqual(ts.num_records, ts_new.num_records + 1)
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
        self.assertIn(t_new.root, roots)


class TestWithVisuals(TopologyTestCase):
    """
    Some pedantic tests with ascii depictions of what's supposed to happen.
    """

    def verify_simplify_topology(self, ts, sample):
        # copies from test_highlevel.py
        new_ts = ts.simplify(sample)
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

        records = [
            cr(left=0.0, right=0.2, node=4, children=(2, 3), time=0.4),
            cr(left=0.2, right=0.8, node=4, children=(0, 2), time=0.4),
            cr(left=0.8, right=1.0, node=4, children=(2, 3), time=0.4),
            cr(left=0.0, right=1.0, node=5, children=(1, 4), time=0.5),
            cr(left=0.8, right=1.0, node=6, children=(0, 5), time=0.7),
            cr(left=0.0, right=0.2, node=7, children=(0, 5), time=1.0),
        ]
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6: -1, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: -1, 4: 5, 5: -1, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        ts = build_tree_sequence(records)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 4)
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

        records = [
            cr(left=0.0, right=0.2, node=3, children=(2, 7), time=0.4),
            cr(left=0.2, right=0.8, node=3, children=(0, 2), time=0.4),
            cr(left=0.8, right=1.0, node=3, children=(2, 7), time=0.4),
            cr(left=0.0, right=0.2, node=4, children=(1, 3), time=0.5),
            cr(left=0.2, right=0.8, node=4, children=(1, 3), time=0.5),
            cr(left=0.8, right=1.0, node=4, children=(1, 3), time=0.5),
            cr(left=0.8, right=1.0, node=5, children=(0, 4), time=0.7),
            cr(left=0.0, right=0.2, node=6, children=(0, 4), time=1.0),
        ]
        true_trees = [
            {0: 6, 1: 4, 2: 3, 3: 4, 4: 6, 5: -1, 6: -1, 7: 3},
            {0: 3, 1: 4, 2: 3, 3: 4, 4: -1, 5: -1, 6: -1, 7: -1},
            {0: 5, 1: 4, 2: 3, 3: 4, 4: 5, 5: -1, 6: -1, 7: 3}]
        ts = build_tree_sequence(records)
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
        records = [
            cr(left=0.0, right=0.2, node=4, children=(2, 3), time=0.4),
            cr(left=0.2, right=0.8, node=4, children=(0, 2), time=0.4),
            cr(left=0.8, right=1.0, node=4, children=(2, 3), time=0.4),
            cr(left=0.0, right=1.0, node=5, children=(1, 4), time=0.5),
            cr(left=0.8, right=1.0, node=6, children=(0, 5), time=0.7),
            cr(left=0.0, right=0.2, node=6, children=(5,), time=0.7),
            cr(left=0.0, right=0.2, node=7, children=(0, 6), time=1.0),
        ]
        true_trees = [
            {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: 7, 7: -1},
            {0: 4, 1: 5, 2: 4, 3: -1, 4: 5, 5: -1, 6: -1, 7: -1},
            {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: -1, 7: -1}]
        ts = build_tree_sequence(records)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        self.assertEqual(ts.sample_size, 4)
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

        records = [
            cr(left=0.5, right=1.0, node=10, children=(1,), time=5.0-4.0),
            cr(left=0.0, right=0.4, node=10, children=(2,), time=5.0-4.0),
            cr(left=0.6, right=1.0, node=9, children=(0,), time=5.0-4.0),
            cr(left=0.0, right=0.5, node=9, children=(1,), time=5.0-4.0),
            cr(left=0.8, right=1.0, node=8, children=(10,), time=5.0-3.0),
            cr(left=0.2, right=0.8, node=8, children=(9, 10), time=5.0-3.0),
            cr(left=0.0, right=0.2, node=8, children=(9,), time=5.0-3.0),
            cr(left=0.7, right=1.0, node=7, children=(8,), time=5.0-2.0),
            cr(left=0.0, right=0.2, node=7, children=(10,), time=5.0-2.0),
            cr(left=0.8, right=1.0, node=6, children=(9,), time=5.0-2.0),
            cr(left=0.0, right=0.7, node=6, children=(8,), time=5.0-2.0),
            cr(left=0.4, right=1.0, node=5, children=(2, 7), time=5.0-1.0),
            cr(left=0.1, right=0.4, node=5, children=(7,), time=5.0-1.0),
            cr(left=0.6, right=0.9, node=4, children=(6,), time=5.0-1.0),
            cr(left=0.0, right=0.6, node=4, children=(0, 6), time=5.0-1.0),
            cr(left=0.9, right=1.0, node=3, children=(4, 5, 6), time=5.0-0.0),
            cr(left=0.1, right=0.9, node=3, children=(4, 5), time=5.0-0.0),
            cr(left=0.0, right=0.1, node=3, children=(4, 5, 7), time=5.0-0.0),
        ]

        ts = build_tree_sequence(records)
        tree_dicts = [t.parent_dict for t in ts.trees()]
        # self.assertEqual(ts.sample_size, 3)
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
