#
# Copyright (C) 2016 University of Oxford
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
Test cases for stats calculations in msprime.
"""
from __future__ import print_function
from __future__ import division

import unittest
import sys

import numpy as np

import msprime
import _msprime
import tests.tsutil as tsutil
import tests.test_wright_fisher as wf


IS_PY2 = sys.version_info[0] < 3


def get_r2_matrix(ts):
    """
    Returns the matrix for the specified tree sequence. This is computed
    via a straightforward Python algorithm.
    """
    n = ts.get_sample_size()
    m = ts.get_num_mutations()
    A = np.zeros((m, m), dtype=float)
    for t1 in ts.trees():
        for sA in t1.sites():
            assert len(sA.mutations) == 1
            mA = sA.mutations[0]
            A[sA.id, sA.id] = 1
            fA = t1.get_num_samples(mA.node) / n
            samples = list(t1.samples(mA.node))
            for t2 in ts.trees(tracked_samples=samples):
                for sB in t2.sites():
                    assert len(sB.mutations) == 1
                    mB = sB.mutations[0]
                    if sB.position > sA.position:
                        fB = t2.get_num_samples(mB.node) / n
                        fAB = t2.get_num_tracked_samples(mB.node) / n
                        D = fAB - fA * fB
                        r2 = D * D / (fA * fB * (1 - fA) * (1 - fB))
                        A[sA.id, sB.id] = r2
                        A[sB.id, sA.id] = r2
    return A


class TestLdCalculator(unittest.TestCase):
    """
    Tests for the LdCalculator class.
    """

    num_test_sites = 50

    def verify_matrix(self, ts):
        m = ts.get_num_sites()
        ldc = msprime.LdCalculator(ts)
        A = ldc.get_r2_matrix()
        self.assertEqual(A.shape, (m, m))
        B = get_r2_matrix(ts)
        self.assertTrue(np.allclose(A, B))

        # Now look at each row in turn, and verify it's the same
        # when we use get_r2 directly.
        for j in range(m):
            a = ldc.get_r2_array(j, direction=msprime.FORWARD)
            b = A[j, j + 1:]
            self.assertEqual(a.shape[0], m - j - 1)
            self.assertEqual(b.shape[0], m - j - 1)
            self.assertTrue(np.allclose(a, b))
            a = ldc.get_r2_array(j, direction=msprime.REVERSE)
            b = A[j, :j]
            self.assertEqual(a.shape[0], j)
            self.assertEqual(b.shape[0], j)
            self.assertTrue(np.allclose(a[::-1], b))

        # Now check every cell in the matrix in turn.
        for j in range(m):
            for k in range(m):
                self.assertAlmostEqual(ldc.get_r2(j, k), A[j, k])

    def verify_max_distance(self, ts):
        """
        Verifies that the max_distance parameter works as expected.
        """
        mutations = list(ts.mutations())
        ldc = msprime.LdCalculator(ts)
        A = ldc.get_r2_matrix()
        j = len(mutations) // 2
        for k in range(j):
            x = mutations[j + k].position - mutations[j].position
            a = ldc.get_r2_array(j, max_distance=x)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j + 1: j + 1 + k], a))
            x = mutations[j].position - mutations[j - k].position
            a = ldc.get_r2_array(j, max_distance=x, direction=msprime.REVERSE)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j - k: j], a[::-1]))
        L = ts.get_sequence_length()
        m = len(mutations)
        a = ldc.get_r2_array(0, max_distance=L)
        self.assertEqual(a.shape[0], m - 1)
        self.assertTrue(np.allclose(A[0, 1:], a))
        a = ldc.get_r2_array(m - 1, max_distance=L, direction=msprime.REVERSE)
        self.assertEqual(a.shape[0], m - 1)
        self.assertTrue(np.allclose(A[m - 1, :-1], a[::-1]))

    def verify_max_mutations(self, ts):
        """
        Verifies that the max mutations parameter works as expected.
        """
        mutations = list(ts.mutations())
        ldc = msprime.LdCalculator(ts)
        A = ldc.get_r2_matrix()
        j = len(mutations) // 2
        for k in range(j):
            a = ldc.get_r2_array(j, max_mutations=k)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j + 1: j + 1 + k], a))
            a = ldc.get_r2_array(j, max_mutations=k, direction=msprime.REVERSE)
            self.assertEqual(a.shape[0], k)
            self.assertTrue(np.allclose(A[j, j - k: j], a[::-1]))

    def test_single_tree_simulated_mutations(self):
        ts = msprime.simulate(20, mutation_rate=10, random_seed=15)
        ts = tsutil.subsample_sites(ts, self.num_test_sites)
        self.verify_matrix(ts)
        self.verify_max_distance(ts)

    def test_deprecated_aliases(self):
        ts = msprime.simulate(20, mutation_rate=10, random_seed=15)
        ts = tsutil.subsample_sites(ts, self.num_test_sites)
        ldc = msprime.LdCalculator(ts)
        A = ldc.get_r2_matrix()
        B = ldc.r2_matrix()
        self.assertTrue(np.array_equal(A, B))
        a = ldc.get_r2_array(0)
        b = ldc.r2_array(0)
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(ldc.get_r2(0, 1), ldc.r2(0, 1))

    def test_single_tree_regular_mutations(self):
        ts = msprime.simulate(self.num_test_sites, length=self.num_test_sites)
        ts = tsutil.insert_branch_mutations(ts)
        # We don't support back mutations, so this should fail.
        self.assertRaises(_msprime.LibraryError, self.verify_matrix, ts)
        self.assertRaises(_msprime.LibraryError, self.verify_max_distance, ts)

    def test_tree_sequence_regular_mutations(self):
        ts = msprime.simulate(
            self.num_test_sites, recombination_rate=1,
            length=self.num_test_sites)
        self.assertGreater(ts.get_num_trees(), 10)
        t = ts.dump_tables()
        t.sites.reset()
        t.mutations.reset()
        for j in range(self.num_test_sites):
            site_id = len(t.sites)
            t.sites.add_row(position=j, ancestral_state="0")
            t.mutations.add_row(site=site_id, derived_state="1", node=j)
        ts = msprime.load_tables(**t.asdict())
        self.verify_matrix(ts)
        self.verify_max_distance(ts)

    def test_tree_sequence_simulated_mutations(self):
        ts = msprime.simulate(20, mutation_rate=10, recombination_rate=10)
        self.assertGreater(ts.get_num_trees(), 10)
        ts = tsutil.subsample_sites(ts, self.num_test_sites)
        self.verify_matrix(ts)
        self.verify_max_distance(ts)
        self.verify_max_mutations(ts)


def set_partitions(collection):
    """
    Returns an ierator over all partitions of the specified set.

    From https://stackoverflow.com/questions/19368375/set-partitions-in-python
    """
    if len(collection) == 1:
        yield [collection]
    else:
        first = collection[0]
        for smaller in set_partitions(collection[1:]):
            for n, subset in enumerate(smaller):
                yield smaller[:n] + [[first] + subset] + smaller[n + 1:]
            yield [[first]] + smaller


def naive_mean_descendants(ts, reference_sets):
    """
    Straightforward implementation of mean sample ancestry by iterating
    over the trees and nodes in each tree.
    """
    # TODO generalise this to allow arbitrary nodes, not just samples.
    C = np.zeros((ts.num_nodes, len(reference_sets)))
    T = np.zeros(ts.num_nodes)
    tree_iters = [ts.trees(tracked_samples=sample_set) for sample_set in reference_sets]
    for _ in range(ts.num_trees):
        trees = [next(tree_iter) for tree_iter in tree_iters]
        left, right = trees[0].interval
        length = right - left
        for node in trees[0].nodes():
            num_samples = trees[0].num_samples(node)
            if num_samples > 0:
                for j, tree in enumerate(trees):
                    C[node, j] += length * tree.num_tracked_samples(node)
                T[node] += length
    for node in range(ts.num_nodes):
        if T[node] > 0:
            C[node] /= T[node]
    return C


class TestMeanDescendants(unittest.TestCase):
    """
    Tests the TreeSequence.mean_descendants method.
    """
    def verify(self, ts, reference_sets):
        C1 = naive_mean_descendants(ts, reference_sets)
        C2 = tsutil.mean_descendants(ts, reference_sets)
        C3 = ts.mean_descendants(reference_sets)
        self.assertEqual(C1.shape, C2.shape)
        self.assertTrue(np.allclose(C1, C2))
        self.assertTrue(np.allclose(C1, C3))
        return C1

    def test_two_populations_high_migration(self):
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(8),
                msprime.PopulationConfiguration(8)],
            migration_matrix=[[0, 1], [1, 0]],
            recombination_rate=3,
            random_seed=5)
        self.assertGreater(ts.num_trees, 1)
        self.verify(ts, [ts.samples(0), ts.samples(1)])

    def test_single_tree(self):
        ts = msprime.simulate(6, random_seed=1)
        S = [range(3), range(3, 6)]
        C = self.verify(ts, S)
        for j, samples in enumerate(S):
            tree = next(ts.trees(tracked_samples=samples))
            for u in tree.nodes():
                self.assertEqual(tree.num_tracked_samples(u), C[u, j])

    def test_single_tree_partial_samples(self):
        ts = msprime.simulate(6, random_seed=1)
        S = [range(3), range(3, 4)]
        C = self.verify(ts, S)
        for j, samples in enumerate(S):
            tree = next(ts.trees(tracked_samples=samples))
            for u in tree.nodes():
                self.assertEqual(tree.num_tracked_samples(u), C[u, j])

    def test_single_tree_all_sample_sets(self):
        ts = msprime.simulate(6, random_seed=1)
        for S in set_partitions(list(range(ts.num_samples))):
            C = self.verify(ts, S)
            for j, samples in enumerate(S):
                tree = next(ts.trees(tracked_samples=samples))
                for u in tree.nodes():
                    self.assertEqual(tree.num_tracked_samples(u), C[u, j])

    def test_many_trees_all_sample_sets(self):
        ts = msprime.simulate(6, recombination_rate=2, random_seed=1)
        self.assertGreater(ts.num_trees, 2)
        for S in set_partitions(list(range(ts.num_samples))):
            self.verify(ts, S)

    def test_wright_fisher_unsimplified_all_sample_sets(self):
        tables = wf.wf_sim(
            4, 5, seed=1, deep_history=False, initial_generation_samples=False,
            num_loci=10)
        tables.sort()
        ts = tables.tree_sequence()
        for S in set_partitions(list(ts.samples())):
            self.verify(ts, S)

    def test_wright_fisher_unsimplified(self):
        tables = wf.wf_sim(
            20, 15, seed=1, deep_history=False, initial_generation_samples=False,
            num_loci=20)
        tables.sort()
        ts = tables.tree_sequence()
        samples = ts.samples()
        self.verify(ts, [samples[:10], samples[10:]])

    def test_wright_fisher_simplified(self):
        tables = wf.wf_sim(
            30, 10, seed=1, deep_history=False, initial_generation_samples=False,
            num_loci=5)
        tables.sort()
        ts = tables.tree_sequence()
        samples = ts.samples()
        self.verify(ts, [samples[:10], samples[10:]])


def naive_genealogical_nearest_neighbours(ts, focal, reference_sets):
    # Make sure everyhing is a sample so we can use the tracked_samples option.
    # This is a limitation of the current API.
    tables = ts.dump_tables()
    tables.nodes.set_columns(
        flags=np.ones_like(tables.nodes.flags),
        time=tables.nodes.time)
    ts = tables.tree_sequence()

    A = np.zeros((len(focal), len(reference_sets)))
    L = np.zeros(len(focal))
    reference_set_map = np.zeros(ts.num_nodes, dtype=int) - 1
    for k, ref_set in enumerate(reference_sets):
        for u in ref_set:
            reference_set_map[u] = k
    tree_iters = [
        ts.trees(tracked_samples=reference_nodes) for reference_nodes in reference_sets]
    for _ in range(ts.num_trees):
        trees = list(map(next, tree_iters))
        length = trees[0].interval[1] - trees[0].interval[0]
        for j, u in enumerate(focal):
            v = trees[0].parent(u)
            while v != msprime.NULL_NODE:
                total = sum(tree.num_tracked_samples(v) for tree in trees)
                if total > 1:
                    break
                v = trees[0].parent(v)
            if v != msprime.NULL_NODE:
                focal_node_set = reference_set_map[u]
                for k, tree in enumerate(trees):
                    # If the focal node is in the current set, we subtract its
                    # contribution from the numerator
                    n = tree.num_tracked_samples(v) - (k == focal_node_set)
                    # If the focal node is in *any* reference set, we subtract its
                    # contribution from the demoninator.
                    A[j, k] += length * n / (total - int(focal_node_set != -1))
                L[j] += length
    # Normalise by the accumulated value for each focal node.
    index = L > 0
    L = L[index]
    L = L.reshape((L.shape[0], 1))
    A[index, :] /= L
    return A


class TestGenealogicalNearestNeighbours(unittest.TestCase):
    """
    Tests the TreeSequence.genealogical_nearest_neighbours method.
    """
    def verify(self, ts, reference_sets, focal=None):
        if focal is None:
            focal = [u for refset in reference_sets for u in refset]
        A1 = naive_genealogical_nearest_neighbours(ts, focal, reference_sets)
        A2 = tsutil.genealogical_nearest_neighbours(ts, focal, reference_sets)
        A3 = ts.genealogical_nearest_neighbours(focal, reference_sets)
        if IS_PY2:
            # Threads not supported on PY2
            self.assertRaises(
                ValueError, ts.genealogical_nearest_neighbours, focal,
                reference_sets, num_threads=3)
        else:
            A4 = ts.genealogical_nearest_neighbours(focal, reference_sets, num_threads=3)
            self.assertTrue(np.array_equal(A3, A4))
        self.assertEqual(A1.shape, A2.shape)
        self.assertEqual(A1.shape, A3.shape)
        self.assertTrue(np.allclose(A1, A2))
        self.assertTrue(np.allclose(A1, A3))
        if all(ts.node(u).is_sample() for u in focal):
            # When the focal nodes are samples, we can assert some stronger properties.
            fully_rooted = True
            for tree in ts.trees():
                if tree.num_roots > 1:
                    fully_rooted = False
                    break
            if fully_rooted:
                self.assertTrue(np.allclose(np.sum(A1, axis=1), 1))
            else:
                all_refs = [u for refset in reference_sets for u in refset]
                # Any node that hits a root before meeting a descendent of the reference
                # nodes must have total zero.
                coalescence_found = np.array([False for _ in all_refs])
                for tree in ts.trees(tracked_samples=all_refs):
                    for j, u in enumerate(focal):
                        while u != msprime.NULL_NODE:
                            if tree.num_tracked_samples(u) > 1:
                                coalescence_found[j] = True
                                break
                            u = tree.parent(u)
                self.assertTrue(np.allclose(np.sum(A1[coalescence_found], axis=1), 1))
                # Anything where there's no coalescence, ever is zero by convention.
                self.assertTrue(
                    np.allclose(
                        np.sum(A1[np.logical_not(coalescence_found)], axis=1), 0))
        return A1

    def test_two_populations_high_migration(self):
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(18),
                msprime.PopulationConfiguration(18)],
            migration_matrix=[[0, 1], [1, 0]],
            recombination_rate=8,
            random_seed=5)
        self.assertGreater(ts.num_trees, 1)
        self.verify(ts, [ts.samples(0), ts.samples(1)])

    def test_single_tree(self):
        ts = msprime.simulate(6, random_seed=1)
        S = [range(3), range(3, 6)]
        self.verify(ts, S)

    def test_single_tree_internal_reference_sets(self):
        ts = msprime.simulate(10, random_seed=1)
        tree = ts.first()
        S = [[u] for u in tree.children(tree.root)]
        self.verify(ts, S, ts.samples())

    def test_single_tree_all_nodes(self):
        ts = msprime.simulate(10, random_seed=1)
        S = [np.arange(ts.num_nodes, dtype=np.int32)]
        self.verify(ts, S, np.arange(ts.num_nodes, dtype=np.int32))

    def test_single_tree_partial_samples(self):
        ts = msprime.simulate(6, random_seed=1)
        S = [range(3), range(3, 4)]
        self.verify(ts, S)

    def test_single_tree_all_sample_sets(self):
        ts = msprime.simulate(6, random_seed=1)
        for S in set_partitions(list(range(ts.num_samples))):
            self.verify(ts, S)

    def test_many_trees_all_sample_sets(self):
        ts = msprime.simulate(6, recombination_rate=2, random_seed=1)
        self.assertGreater(ts.num_trees, 2)
        for S in set_partitions(list(range(ts.num_samples))):
            self.verify(ts, S)

    def test_many_trees_sequence_length(self):
        for L in [0.5, 1.5, 3.3333]:
            ts = msprime.simulate(6, length=L, recombination_rate=2, random_seed=1)
            self.verify(ts, [range(3), range(3, 6)])

    def test_many_trees_all_nodes(self):
        ts = msprime.simulate(6, length=4, recombination_rate=2, random_seed=1)
        S = [np.arange(ts.num_nodes, dtype=np.int32)]
        self.verify(ts, S, np.arange(ts.num_nodes, dtype=np.int32))

    def test_wright_fisher_unsimplified_all_sample_sets(self):
        tables = wf.wf_sim(
            4, 5, seed=1, deep_history=True, initial_generation_samples=False,
            num_loci=10)
        tables.sort()
        ts = tables.tree_sequence()
        for S in set_partitions(list(ts.samples())):
            self.verify(ts, S)

    def test_wright_fisher_unsimplified(self):
        tables = wf.wf_sim(
            20, 15, seed=1, deep_history=True, initial_generation_samples=False,
            num_loci=20)
        tables.sort()
        ts = tables.tree_sequence()
        samples = ts.samples()
        self.verify(ts, [samples[:10], samples[10:]])

    def test_wright_fisher_initial_generation(self):
        tables = wf.wf_sim(
            20, 15, seed=1, deep_history=True, initial_generation_samples=True,
            num_loci=20)
        tables.sort()
        tables.simplify()
        ts = tables.tree_sequence()
        samples = ts.samples()
        founders = [u for u in samples if ts.node(u).time > 0]
        samples = [u for u in samples if ts.node(u).time == 0]
        self.verify(ts, [founders[:10], founders[10:]], samples)

    def test_wright_fisher_initial_generation_no_deep_history(self):
        tables = wf.wf_sim(
            20, 15, seed=2, deep_history=False, initial_generation_samples=True,
            num_loci=20)
        tables.sort()
        tables.simplify()
        ts = tables.tree_sequence()
        samples = ts.samples()
        founders = [u for u in samples if ts.node(u).time > 0]
        samples = [u for u in samples if ts.node(u).time == 0]
        A = self.verify(ts, [founders[:10], founders[10:]], samples)
        # Because the founders are all isolated, the stat must be zero.
        self.assertTrue(np.all(A == 0))

    def test_wright_fisher_unsimplified_multiple_roots(self):
        tables = wf.wf_sim(
            20, 15, seed=1, deep_history=False, initial_generation_samples=False,
            num_loci=20)
        tables.sort()
        ts = tables.tree_sequence()
        samples = ts.samples()
        self.verify(ts, [samples[:10], samples[10:]])

    def test_wright_fisher_simplified(self):
        tables = wf.wf_sim(
            31, 10, seed=1, deep_history=True, initial_generation_samples=False,
            num_loci=5)
        tables.sort()
        ts = tables.tree_sequence().simplify()
        samples = ts.samples()
        self.verify(ts, [samples[:10], samples[10:]])

    def test_wright_fisher_simplified_multiple_roots(self):
        tables = wf.wf_sim(
            31, 10, seed=1, deep_history=False, initial_generation_samples=False,
            num_loci=5)
        tables.sort()
        ts = tables.tree_sequence()
        samples = ts.samples()
        self.verify(ts, [samples[:10], samples[10:]])

    def test_empty_ts(self):
        tables = msprime.TableCollection(1.0)
        tables.nodes.add_row(1, 0)
        tables.nodes.add_row(1, 0)
        ts = tables.tree_sequence()
        self.verify(ts, [[0], [1]])


def exact_genealogical_nearest_neighbours(ts, focal, reference_sets):
    # Same as above, except we return the per-tree value for a single node.

    # Make sure everyhing is a sample so we can use the tracked_samples option.
    # This is a limitation of the current API.
    tables = ts.dump_tables()
    tables.nodes.set_columns(
        flags=np.ones_like(tables.nodes.flags),
        time=tables.nodes.time)
    ts = tables.tree_sequence()

    A = np.zeros((len(reference_sets), ts.num_trees))
    L = np.zeros(ts.num_trees)
    reference_set_map = np.zeros(ts.num_nodes, dtype=int) - 1
    for k, ref_set in enumerate(reference_sets):
        for u in ref_set:
            reference_set_map[u] = k
    tree_iters = [
        ts.trees(tracked_samples=reference_nodes) for reference_nodes in reference_sets]
    u = focal
    for _ in range(ts.num_trees):
        trees = list(map(next, tree_iters))
        v = trees[0].parent(u)
        while v != msprime.NULL_NODE:
            total = sum(tree.num_tracked_samples(v) for tree in trees)
            if total > 1:
                break
            v = trees[0].parent(v)
        if v != msprime.NULL_NODE:
            # The length is only reported where the statistic is defined.
            L[trees[0].index] = trees[0].interval[1] - trees[0].interval[0]
            focal_node_set = reference_set_map[u]
            for k, tree in enumerate(trees):
                # If the focal node is in the current set, we subtract its
                # contribution from the numerator
                n = tree.num_tracked_samples(v) - (k == focal_node_set)
                # If the focal node is in *any* reference set, we subtract its
                # contribution from the demoninator.
                A[k, tree.index] = n / (total - int(focal_node_set != -1))
    return A, L


class TestExactGenealogicalNearestNeighbours(TestGenealogicalNearestNeighbours):

    def verify(self, ts, reference_sets, focal=None):
        if focal is None:
            focal = [u for refset in reference_sets for u in refset]
        A = ts.genealogical_nearest_neighbours(focal, reference_sets)

        for j, u in enumerate(focal):
            T, L = exact_genealogical_nearest_neighbours(ts, u, reference_sets)
            # Ignore the cases where the node has no GNNs
            if np.sum(L) > 0:
                mean = np.sum(T * L, axis=1) / np.sum(L)
                self.assertTrue(np.allclose(mean, A[j]))
        return A
