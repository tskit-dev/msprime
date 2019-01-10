"""
Test cases for the low level C interface to tskit.
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

import collections
import itertools
import os
import platform
import random
import sys
import tempfile
import unittest

import msprime

import _tskit

IS_PY2 = sys.version_info[0] < 3
IS_WINDOWS = platform.system() == "Windows"


def get_tracked_sample_counts(st, tracked_samples):
    """
    Returns a list giving the number of samples in the specified list
    that are in the subtree rooted at each node.
    """
    nu = [0 for j in range(st.get_num_nodes())]
    for j in tracked_samples:
        # Duplicates not permitted.
        assert nu[j] == 0
        u = j
        while u != _tskit.NULL:
            nu[u] += 1
            u = st.get_parent(u)
    return nu


def get_sample_counts(tree_sequence, st):
    """
    Returns a list of the sample node counts for the specfied sparse tree.
    """
    nu = [0 for j in range(st.get_num_nodes())]
    for j in range(tree_sequence.get_num_samples()):
        u = j
        while u != _tskit.NULL:
            nu[u] += 1
            u = st.get_parent(u)
    return nu


class LowLevelTestCase(unittest.TestCase):
    """
    Superclass of tests for the low-level interface.
    """
    def verify_tree_dict(self, n, pi):
        """
        Verifies that the specified sparse tree in dict format is a
        consistent coalescent history for a sample of size n.
        """
        self.assertLessEqual(len(pi), 2 * n - 1)
        # _tskit.NULL should not be a node
        self.assertNotIn(_tskit.NULL, pi)
        # verify the root is equal for all samples
        root = 0
        while pi[root] != _tskit.NULL:
            root = pi[root]
        for j in range(n):
            k = j
            while pi[k] != _tskit.NULL:
                k = pi[k]
            self.assertEqual(k, root)
        # 0 to n - 1 inclusive should always be nodes
        for j in range(n):
            self.assertIn(j, pi)
        num_children = collections.defaultdict(int)
        for j in pi.keys():
            num_children[pi[j]] += 1
        # nodes 0 to n are samples.
        for j in range(n):
            self.assertNotEqual(pi[j], 0)
            self.assertEqual(num_children[j], 0)
        # All non-sample nodes should be binary
        for j in pi.keys():
            if j > n:
                self.assertGreaterEqual(num_children[j], 2)

    def get_example_tree_sequence(self):
        ts = msprime.simulate(10, recombination_rate=0.1, random_seed=1)
        return ts.ll_tree_sequence

    def get_example_tree_sequences(self):
        yield self.get_example_tree_sequence()
        yield self.get_example_migration_tree_sequence()

    def get_example_migration_tree_sequence(self):
        pop_configs = [msprime.PopulationConfiguration(5) for _ in range(2)]
        migration_matrix = [[0, 1], [1, 0]]
        ts = msprime.simulate(
            population_configurations=pop_configs,
            migration_matrix=migration_matrix,
            mutation_rate=1,
            record_migrations=True,
            random_seed=1)
        return ts.ll_tree_sequence

    def verify_iterator(self, iterator):
        """
        Checks that the specified non-empty iterator implements the
        iterator protocol correctly.
        """
        list_ = list(iterator)
        self.assertGreater(len(list_), 0)
        for j in range(10):
            self.assertRaises(StopIteration, next, iterator)


class TestTableCollection(LowLevelTestCase):
    """
    Tests for the low-level TableCollection class
    """

    def test_reference_deletion(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        tc = ts.tables.ll_tables
        # Get references to all the tables
        tables = [
            tc.individuals, tc.nodes, tc.edges, tc.migrations, tc.sites, tc.mutations,
            tc.populations, tc.provenances]
        del tc
        for _ in range(10):
            for table in tables:
                self.assertGreater(len(str(table)), 0)


class TestTreeSequence(LowLevelTestCase):
    """
    Tests for the low-level interface for the TreeSequence.
    """
    def setUp(self):
        fd, self.temp_file = tempfile.mkstemp(prefix="msp_ll_ts_")
        os.close(fd)

    def tearDown(self):
        os.unlink(self.temp_file)

    @unittest.skipIf(IS_WINDOWS, "File permissions on Windows")
    def test_file_errors(self):
        ts1 = self.get_example_tree_sequence()

        def loader(*args):
            ts2 = _tskit.TreeSequence()
            ts2.load(*args)

        for func in [ts1.dump, loader]:
            self.assertRaises(TypeError, func)
            for bad_type in [1, None, [], {}]:
                self.assertRaises(TypeError, func, bad_type)
            # Try to dump/load files we don't have access to or don't exist.
            for f in ["/", "/test.trees", "/dir_does_not_exist/x.trees"]:
                self.assertRaises(_tskit.FileFormatError, func, f)
                try:
                    func(f)
                except _tskit.FileFormatError as e:
                    message = str(e)
                    self.assertGreater(len(message), 0)
            # use a long filename and make sure we don't overflow error
            # buffers
            f = "/" + 4000 * "x"
            self.assertRaises(_tskit.FileFormatError, func, f)
            try:
                func(f)
            except _tskit.FileFormatError as e:
                message = str(e)
                self.assertLess(len(message), 1024)

    def test_initial_state(self):
        # Check the initial state to make sure that it is empty.
        ts = _tskit.TreeSequence()
        self.assertRaises(ValueError, ts.get_num_samples)
        self.assertRaises(ValueError, ts.get_sequence_length)
        self.assertRaises(ValueError, ts.get_num_trees)
        self.assertRaises(ValueError, ts.get_num_edges)
        self.assertRaises(ValueError, ts.get_num_mutations)
        self.assertRaises(ValueError, ts.get_num_migrations)
        self.assertRaises(ValueError, ts.get_num_migrations)
        self.assertRaises(ValueError, ts.get_genotype_matrix)
        self.assertRaises(ValueError, ts.dump)

    def test_num_nodes(self):
        for ts in self.get_example_tree_sequences():
            max_node = 0
            for j in range(ts.get_num_edges()):
                _, _, parent, child = ts.get_edge(j)
                for node in [parent, child]:
                    if node > max_node:
                        max_node = node
            self.assertEqual(max_node + 1, ts.get_num_nodes())

    def verify_dump_equality(self, ts):
        """
        Verifies that we can dump a copy of the specified tree sequence
        to the specified file, and load an identical copy.
        """
        ts.dump(self.temp_file)
        ts2 = _tskit.TreeSequence()
        ts2.load(self.temp_file)
        self.assertEqual(ts.get_num_samples(), ts2.get_num_samples())
        self.assertEqual(ts.get_sequence_length(), ts2.get_sequence_length())
        self.assertEqual(ts.get_num_mutations(), ts2.get_num_mutations())
        self.assertEqual(ts.get_num_nodes(), ts2.get_num_nodes())
        records1 = [ts.get_edge(j) for j in range(ts.get_num_edges())]
        records2 = [ts2.get_edge(j) for j in range(ts2.get_num_edges())]
        self.assertEqual(records1, records2)
        mutations1 = [ts.get_mutation(j) for j in range(ts.get_num_mutations())]
        mutations2 = [ts2.get_mutation(j) for j in range(ts2.get_num_mutations())]
        self.assertEqual(mutations1, mutations2)
        provenances1 = [ts.get_provenance(j) for j in range(ts.get_num_provenances())]
        provenances2 = [ts2.get_provenance(j) for j in range(ts2.get_num_provenances())]
        self.assertEqual(provenances1, provenances2)

    def test_dump_equality(self):
        for ts in self.get_example_tree_sequences():
            self.verify_dump_equality(ts)

    def verify_mutations(self, ts):
        mutations = [ts.get_mutation(j) for j in range(ts.get_num_mutations())]
        self.assertGreater(ts.get_num_mutations(), 0)
        self.assertEqual(len(mutations), ts.get_num_mutations())
        # Check the form of the mutations
        for j, (position, nodes, index) in enumerate(mutations):
            self.assertEqual(j, index)
            for node in nodes:
                self.assertIsInstance(node, int)
                self.assertGreaterEqual(node, 0)
                self.assertLessEqual(node, ts.get_num_nodes())
            self.assertIsInstance(position, float)
            self.assertGreater(position, 0)
            self.assertLess(position, ts.get_sequence_length())
        # mutations must be sorted by position order.
        self.assertEqual(mutations, sorted(mutations))

    def test_get_edge_interface(self):
        for ts in self.get_example_tree_sequences():
            num_edges = ts.get_num_edges()
            # We don't accept Python negative indexes here.
            self.assertRaises(IndexError, ts.get_edge, -1)
            for j in [0, 10, 10**6]:
                self.assertRaises(IndexError, ts.get_edge, num_edges + j)
            for x in [None, "", {}, []]:
                self.assertRaises(TypeError, ts.get_edge, x)

    def test_get_node_interface(self):
        for ts in self.get_example_tree_sequences():
            num_nodes = ts.get_num_nodes()
            # We don't accept Python negative indexes here.
            self.assertRaises(IndexError, ts.get_node, -1)
            for j in [0, 10, 10**6]:
                self.assertRaises(IndexError, ts.get_node, num_nodes + j)
            for x in [None, "", {}, []]:
                self.assertRaises(TypeError, ts.get_node, x)

    def test_get_genotype_matrix_interface(self):
        for ts in self.get_example_tree_sequences():
            num_samples = ts.get_num_samples()
            num_sites = ts.get_num_sites()
            G = ts.get_genotype_matrix()
            self.assertEqual(G.shape, (num_sites, num_samples))

    def test_get_migration_interface(self):
        ts = self.get_example_migration_tree_sequence()
        for bad_type in ["", None, {}]:
            self.assertRaises(TypeError, ts.get_migration, bad_type)
        num_records = ts.get_num_migrations()
        # We don't accept Python negative indexes here.
        self.assertRaises(IndexError, ts.get_migration, -1)
        for j in [0, 10, 10**6]:
            self.assertRaises(IndexError, ts.get_migration, num_records + j)

    def test_get_samples(self):
        ts = self.get_example_migration_tree_sequence()
        # get_samples takes no arguments.
        self.assertRaises(TypeError, ts.get_samples, 0)
        self.assertEqual(list(range(ts.get_num_samples())), ts.get_samples())

    def test_pairwise_diversity(self):
        for ts in self.get_example_tree_sequences():
            for bad_type in ["", None, {}]:
                self.assertRaises(
                    TypeError, ts.get_pairwise_diversity, bad_type)
            self.assertRaises(
                ValueError, ts.get_pairwise_diversity, [])
            self.assertRaises(
                ValueError, ts.get_pairwise_diversity, [0])
            self.assertRaises(
                ValueError, ts.get_pairwise_diversity,
                [0, ts.get_num_samples()])
            self.assertRaises(
                _tskit.LibraryError, ts.get_pairwise_diversity, [0, 0])
            samples = list(range(ts.get_num_samples()))
            pi1 = ts.get_pairwise_diversity(samples)
            self.assertGreaterEqual(pi1, 0)

    def test_genealogical_nearest_neighbours(self):
        for ts in self.get_example_tree_sequences():
            self.assertRaises(TypeError, ts.genealogical_nearest_neighbours)
            self.assertRaises(
                TypeError, ts.genealogical_nearest_neighbours, focal=None)
            self.assertRaises(
                TypeError, ts.genealogical_nearest_neighbours, focal=ts.get_samples(),
                reference_sets={})
            self.assertRaises(
                ValueError, ts.genealogical_nearest_neighbours, focal=ts.get_samples(),
                reference_sets=[])

            bad_array_values = ["", {}, "x", [[[0], [1, 2]]]]
            for bad_array_value in bad_array_values:
                self.assertRaises(
                    ValueError, ts.genealogical_nearest_neighbours,
                    focal=bad_array_value, reference_sets=[[0], [1]])
                self.assertRaises(
                    ValueError, ts.genealogical_nearest_neighbours,
                    focal=ts.get_samples(), reference_sets=[[0], bad_array_value])
                self.assertRaises(
                    ValueError, ts.genealogical_nearest_neighbours,
                    focal=ts.get_samples(), reference_sets=[bad_array_value])
            focal = ts.get_samples()
            A = ts.genealogical_nearest_neighbours(focal, [focal[2:], focal[:2]])
            self.assertEqual(A.shape, (len(focal), 2))

    def test_mean_descendants(self):
        for ts in self.get_example_tree_sequences():
            self.assertRaises(TypeError, ts.mean_descendants)
            self.assertRaises(TypeError, ts.mean_descendants, reference_sets={})
            self.assertRaises(ValueError, ts.mean_descendants, reference_sets=[])

            bad_array_values = ["", {}, "x", [[[0], [1, 2]]]]
            for bad_array_value in bad_array_values:
                self.assertRaises(
                    ValueError, ts.mean_descendants,
                    reference_sets=[[0], bad_array_value])
                self.assertRaises(
                    ValueError, ts.mean_descendants, reference_sets=[bad_array_value])
            focal = ts.get_samples()
            A = ts.mean_descendants([focal[2:], focal[:2]])
            self.assertEqual(A.shape, (ts.get_num_nodes(), 2))


class TestTreeDiffIterator(LowLevelTestCase):
    """
    Tests for the low-level tree diff iterator.
    """
    def test_uninitialised_tree_sequence(self):
        ts = _tskit.TreeSequence()
        self.assertRaises(ValueError, _tskit.TreeDiffIterator, ts)

    def test_constructor(self):
        self.assertRaises(TypeError, _tskit.TreeDiffIterator)
        self.assertRaises(TypeError, _tskit.TreeDiffIterator, None)
        ts = self.get_example_tree_sequence()
        before = list(_tskit.TreeDiffIterator(ts))
        iterator = _tskit.TreeDiffIterator(ts)
        del ts
        # We should keep a reference to the tree sequence.
        after = list(iterator)
        self.assertEqual(before, after)

    def test_iterator(self):
        ts = self.get_example_tree_sequence()
        self.verify_iterator(_tskit.TreeDiffIterator(ts))


class TestTreeIterator(LowLevelTestCase):
    """
    Tests for the low-level sparse tree iterator.
    """
    def test_uninitialised_tree_sequence(self):
        ts = _tskit.TreeSequence()
        self.assertRaises(ValueError, _tskit.Tree, ts)

    def test_constructor(self):
        self.assertRaises(TypeError, _tskit.TreeIterator)
        self.assertRaises(TypeError, _tskit.TreeIterator, None)
        ts = _tskit.TreeSequence()
        self.assertRaises(TypeError, _tskit.TreeIterator, ts)
        ts = self.get_example_tree_sequence()
        tree = _tskit.Tree(ts)
        n_before = 0
        parents_before = []
        for t in _tskit.TreeIterator(tree):
            n_before += 1
            self.assertIs(t, tree)
            pi = {}
            for j in range(t.get_num_nodes()):
                pi[j] = t.get_parent(j)
            parents_before.append(pi)
        self.assertEqual(n_before, len(list(_tskit.TreeDiffIterator(ts))))
        # If we remove the objects, we should get the same results.
        iterator = _tskit.TreeIterator(tree)
        del tree
        del ts
        n_after = 0
        parents_after = []
        for index, t in enumerate(iterator):
            n_after += 1
            self.assertIsInstance(t, _tskit.Tree)
            pi = {}
            for j in range(t.get_num_nodes()):
                pi[j] = t.get_parent(j)
            parents_after.append(pi)
            self.assertEqual(index, t.get_index())
        self.assertEqual(parents_before, parents_after)

    def test_iterator(self):
        ts = self.get_example_tree_sequence()
        tree = _tskit.Tree(ts)
        self.verify_iterator(_tskit.TreeIterator(tree))


class TestTree(LowLevelTestCase):
    """
    Tests on the low-level sparse tree interface.
    """

    def test_flags(self):
        ts = self.get_example_tree_sequence()
        st = _tskit.Tree(ts)
        self.assertEqual(st.get_flags(), 0)
        # We should still be able to count the samples, just inefficiently.
        self.assertEqual(st.get_num_samples(0), 1)
        self.assertRaises(_tskit.LibraryError, st.get_num_tracked_samples, 0)
        all_flags = [
            0, _tskit.SAMPLE_COUNTS, _tskit.SAMPLE_LISTS,
            _tskit.SAMPLE_COUNTS | _tskit.SAMPLE_LISTS]
        for flags in all_flags:
            st = _tskit.Tree(ts, flags=flags)
            self.assertEqual(st.get_flags(), flags)
            self.assertEqual(st.get_num_samples(0), 1)
            if flags & _tskit.SAMPLE_COUNTS:
                self.assertEqual(st.get_num_tracked_samples(0), 0)
            else:
                self.assertRaises(_tskit.LibraryError, st.get_num_tracked_samples, 0)
            if flags & _tskit.SAMPLE_LISTS:
                self.assertEqual(0, st.get_left_sample(0))
                self.assertEqual(0, st.get_right_sample(0))
            else:
                self.assertRaises(ValueError, st.get_left_sample, 0)
                self.assertRaises(ValueError, st.get_right_sample, 0)
                self.assertRaises(ValueError, st.get_next_sample, 0)

    def test_sites(self):
        for ts in self.get_example_tree_sequences():
            st = _tskit.Tree(ts)
            all_sites = [ts.get_site(j) for j in range(ts.get_num_sites())]
            all_tree_sites = []
            j = 0
            mutation_id = 0
            for st in _tskit.TreeIterator(st):
                tree_sites = st.get_sites()
                self.assertEqual(st.get_num_sites(), len(tree_sites))
                all_tree_sites.extend(tree_sites)
                for position, ancestral_state, mutations, index, metadata in tree_sites:
                    self.assertTrue(st.get_left() <= position < st.get_right())
                    self.assertEqual(index, j)
                    self.assertEqual(metadata, b"")
                    for mut_id in mutations:
                        site, node, derived_state, parent, metadata = \
                            ts.get_mutation(mut_id)
                        self.assertEqual(site, index)
                        self.assertEqual(mutation_id, mut_id)
                        self.assertNotEqual(st.get_parent(node), _tskit.NULL)
                        self.assertEqual(metadata, b"")
                        mutation_id += 1
                    j += 1
            self.assertEqual(all_tree_sites, all_sites)

    def test_constructor(self):
        self.assertRaises(TypeError, _tskit.Tree)
        for bad_type in ["", {}, [], None, 0]:
            self.assertRaises(
                TypeError, _tskit.Tree, bad_type)
        ts = self.get_example_tree_sequence()
        for bad_type in ["", {}, True, 1, None]:
            self.assertRaises(
                TypeError, _tskit.Tree, ts, tracked_samples=bad_type)
        for bad_type in ["", {}, None, []]:
            self.assertRaises(
                TypeError, _tskit.Tree, ts, flags=bad_type)
        for ts in self.get_example_tree_sequences():
            st = _tskit.Tree(ts)
            self.assertEqual(st.get_num_nodes(), ts.get_num_nodes())
            # An uninitialised sparse tree should always be zero.
            self.assertEqual(st.get_left_root(), 0)
            self.assertEqual(st.get_left(), 0)
            self.assertEqual(st.get_right(), 0)
            for j in range(ts.get_num_samples()):
                self.assertEqual(st.get_parent(j), _tskit.NULL)
                self.assertEqual(st.get_children(j), tuple())
                self.assertEqual(st.get_time(j), 0)

    def test_memory_error(self):
        # This provokes a bug where we weren't reference counting
        # the tree sequence properly, and the underlying memory for a
        # sparse tree was getting corrupted.
        for ts in self.get_example_tree_sequences():
            num_nodes = ts.get_num_nodes()
            st = _tskit.Tree(ts)
            # deleting the tree sequence should still give a well formed
            # sparse tree.
            st_iter = _tskit.TreeIterator(st)
            next(st_iter)
            del ts
            del st_iter
            # Do a quick traversal just to exercise the tree
            stack = [st.get_left_root()]
            while len(stack) > 0:
                u = stack.pop()
                self.assertLess(u, num_nodes)
                stack.extend(st.get_children(u))

    def test_bad_tracked_samples(self):
        ts = self.get_example_tree_sequence()
        flags = _tskit.SAMPLE_COUNTS
        for bad_type in ["", {}, [], None]:
            self.assertRaises(
                TypeError, _tskit.Tree, ts, flags=flags,
                tracked_samples=[bad_type])
            self.assertRaises(
                TypeError, _tskit.Tree, ts, flags=flags,
                tracked_samples=[1, bad_type])
        for bad_sample in [10**6, -1e6]:
            self.assertRaises(
                ValueError, _tskit.Tree, ts, flags=flags,
                tracked_samples=[bad_sample])
            self.assertRaises(
                ValueError, _tskit.Tree, ts, flags=flags,
                tracked_samples=[1, bad_sample])
            self.assertRaises(
                ValueError, _tskit.Tree, ts,
                tracked_samples=[1, bad_sample, 1])

    def test_count_all_samples(self):
        for ts in self.get_example_tree_sequences():
            self.verify_iterator(_tskit.TreeDiffIterator(ts))
            st = _tskit.Tree(ts, flags=_tskit.SAMPLE_COUNTS)
            # Without initialisation we should be 0 samples for every node
            # that is not a sample.
            for j in range(st.get_num_nodes()):
                count = 1 if j < ts.get_num_samples() else 0
                self.assertEqual(st.get_num_samples(j), count)
                self.assertEqual(st.get_num_tracked_samples(j), 0)
            # Now, try this for a tree sequence.
            for st in _tskit.TreeIterator(st):
                nu = get_sample_counts(ts, st)
                nu_prime = [
                    st.get_num_samples(j) for j in
                    range(st.get_num_nodes())]
                self.assertEqual(nu, nu_prime)
                # For tracked samples, this should be all zeros.
                nu = [
                    st.get_num_tracked_samples(j) for j in
                    range(st.get_num_nodes())]
                self.assertEqual(nu, list([0 for _ in nu]))

    def test_count_tracked_samples(self):
        # Ensure that there are some non-binary nodes.
        non_binary = False
        for ts in self.get_example_tree_sequences():
            st = _tskit.Tree(ts)
            for st in _tskit.TreeIterator(st):
                for u in range(ts.get_num_nodes()):
                    if len(st.get_children(u)) > 1:
                        non_binary = True
            samples = [j for j in range(ts.get_num_samples())]
            powerset = itertools.chain.from_iterable(
                itertools.combinations(samples, r)
                for r in range(len(samples) + 1))
            for subset in map(list, powerset):
                # Ordering shouldn't make any different.
                random.shuffle(subset)
                st = _tskit.Tree(
                    ts, flags=_tskit.SAMPLE_COUNTS, tracked_samples=subset)
                for st in _tskit.TreeIterator(st):
                    nu = get_tracked_sample_counts(st, subset)
                    nu_prime = [
                        st.get_num_tracked_samples(j) for j in
                        range(st.get_num_nodes())]
                    self.assertEqual(nu, nu_prime)
            # Passing duplicated values should raise an error
            sample = 1
            for j in range(2, 20):
                tracked_samples = [sample for _ in range(j)]
                self.assertRaises(
                    _tskit.LibraryError, _tskit.Tree,
                    ts, flags=_tskit.SAMPLE_COUNTS,
                    tracked_samples=tracked_samples)
        self.assertTrue(non_binary)

    def test_bounds_checking(self):
        for ts in self.get_example_tree_sequences():
            n = ts.get_num_nodes()
            st = _tskit.Tree(
                ts, flags=_tskit.SAMPLE_COUNTS | _tskit.SAMPLE_LISTS)
            for v in [-100, -1, n + 1, n + 100, n * 100]:
                self.assertRaises(ValueError, st.get_parent, v)
                self.assertRaises(ValueError, st.get_children, v)
                self.assertRaises(ValueError, st.get_time, v)
                self.assertRaises(ValueError, st.get_left_sample, v)
                self.assertRaises(ValueError, st.get_right_sample, v)
            n = ts.get_num_samples()
            for v in [-100, -1, n + 1, n + 100, n * 100]:
                self.assertRaises(ValueError, st.get_next_sample, v)

    def test_mrca_interface(self):
        for ts in self.get_example_tree_sequences():
            num_nodes = ts.get_num_nodes()
            st = _tskit.Tree(ts)
            for v in [num_nodes, 10**6, _tskit.NULL]:
                self.assertRaises(ValueError, st.get_mrca, v, v)
                self.assertRaises(ValueError, st.get_mrca, v, 1)
                self.assertRaises(ValueError, st.get_mrca, 1, v)
            # All the mrcas for an uninitialised tree should be _tskit.NULL
            for u, v in itertools.combinations(range(num_nodes), 2):
                self.assertEqual(st.get_mrca(u, v), _tskit.NULL)

    def test_newick_precision(self):

        def get_times(tree):
            """
            Returns the time strings from the specified newick tree.
            """
            ret = []
            current_time = None
            for c in tree:
                if c == ":":
                    current_time = ""
                elif c in [",", ")"]:
                    ret.append(current_time)
                    current_time = None
                elif current_time is not None:
                    current_time += c
            return ret

        ts = self.get_example_tree_sequence()
        st = _tskit.Tree(ts)
        for st in _tskit.TreeIterator(st):
            self.assertRaises(ValueError, st.get_newick, root=0, precision=-1)
            self.assertRaises(ValueError, st.get_newick, root=0, precision=17)
            self.assertRaises(ValueError, st.get_newick, root=0, precision=100)
            for precision in range(17):
                tree = st.get_newick(
                    root=st.get_left_root(), precision=precision).decode()
                times = get_times(tree)
                self.assertGreater(len(times), ts.get_num_samples())
                for t in times:
                    if precision == 0:
                        self.assertNotIn(".", t)
                    else:
                        point = t.find(".")
                        self.assertEqual(precision, len(t) - point - 1)

    @unittest.skip("Correct initialisation for sparse tree.")
    def test_newick_interface(self):
        ts = self.get_tree_sequence(num_loci=10, num_samples=10)
        st = _tskit.Tree(ts)
        # TODO this will break when we correctly handle multiple roots.
        self.assertEqual(st.get_newick(), b"1;")
        for bad_type in [None, "", [], {}]:
            self.assertRaises(TypeError, st.get_newick, precision=bad_type)
            self.assertRaises(TypeError, st.get_newick, ts, time_scale=bad_type)
        for st in _tskit.TreeIterator(st):
            newick = st.get_newick()
            self.assertTrue(newick.endswith(b";"))

    def test_index(self):
        for ts in self.get_example_tree_sequences():
            st = _tskit.Tree(ts)
            for index, st in enumerate(_tskit.TreeIterator(st)):
                self.assertEqual(index, st.get_index())

    def test_bad_mutations(self):
        ts = self.get_example_tree_sequence()
        tables = _tskit.TableCollection()
        ts.dump_tables(tables)

        def f(mutations):
            position = []
            node = []
            site = []
            ancestral_state = []
            ancestral_state_offset = [0]
            derived_state = []
            derived_state_offset = [0]
            for j, (p, n) in enumerate(mutations):
                site.append(j)
                position.append(p)
                ancestral_state.append("0")
                ancestral_state_offset.append(ancestral_state_offset[-1] + 1)
                derived_state.append("1")
                derived_state_offset.append(derived_state_offset[-1] + 1)
                node.append(n)
            tables.sites.set_columns(dict(
                position=position, ancestral_state=ancestral_state,
                ancestral_state_offset=ancestral_state_offset,
                metadata=None, metadata_offset=None))
            tables.mutations.set_columns(dict(
                site=site, node=node, derived_state=derived_state,
                derived_state_offset=derived_state_offset,
                parent=None, metadata=None, metadata_offset=None))
            ts2 = _tskit.TreeSequence()
            ts2.load_tables(tables)
        self.assertRaises(_tskit.LibraryError, f, [(0.1, -1)])
        length = ts.get_sequence_length()
        u = ts.get_num_nodes()
        for bad_node in [u, u + 1, 2 * u]:
            self.assertRaises(_tskit.LibraryError, f, [(0.1, bad_node)])
        for bad_pos in [-1, length, length + 1]:
            self.assertRaises(_tskit.LibraryError, f, [(length, 0)])

    def test_free(self):
        ts = self.get_example_tree_sequence()
        t = _tskit.Tree(
            ts, flags=_tskit.SAMPLE_COUNTS | _tskit.SAMPLE_LISTS)
        no_arg_methods = [
            t.get_left_root, t.get_index, t.get_left, t.get_right,
            t.get_num_sites, t.get_flags, t.get_sites, t.get_num_nodes]
        node_arg_methods = [
            t.get_parent, t.get_population, t.get_children, t.get_num_samples,
            t.get_num_tracked_samples]
        two_node_arg_methods = [t.get_mrca]
        for method in no_arg_methods:
            method()
        for method in node_arg_methods:
            method(0)
        for method in two_node_arg_methods:
            method(0, 0)
        t.free()
        self.assertRaises(RuntimeError, t.free)
        for method in no_arg_methods:
            self.assertRaises(RuntimeError, method)
        for method in node_arg_methods:
            self.assertRaises(RuntimeError, method, 0)
        for method in two_node_arg_methods:
            self.assertRaises(RuntimeError, method, 0, 0)

    def test_sample_list(self):
        flags = _tskit.SAMPLE_COUNTS | _tskit.SAMPLE_LISTS
        # Note: we're assuming that samples are 0-n here.
        for ts in self.get_example_tree_sequences():
            st = _tskit.Tree(ts, flags=flags)
            for t in _tskit.TreeIterator(st):
                # All sample nodes should have themselves.
                for j in range(ts.get_num_samples()):
                    self.assertEqual(t.get_left_sample(j), j)
                    self.assertEqual(t.get_right_sample(j), j)

                # All non-tree nodes should have 0
                for j in range(t.get_num_nodes()):
                    if t.get_parent(j) == _tskit.NULL \
                            and t.get_left_child(j) == _tskit.NULL:
                        self.assertEqual(t.get_left_sample(j), _tskit.NULL)
                        self.assertEqual(t.get_right_sample(j), _tskit.NULL)
                # The roots should have all samples.
                u = t.get_left_root()
                samples = []
                while u != _tskit.NULL:
                    sample = t.get_left_sample(u)
                    end = t.get_right_sample(u)
                    while True:
                        samples.append(sample)
                        if sample == end:
                            break
                        sample = t.get_next_sample(sample)
                    u = t.get_right_sib(u)
                self.assertEqual(sorted(samples), list(range(ts.get_num_samples())))


class TestModuleFunctions(unittest.TestCase):
    """
    Tests for the module level functions.
    """
    def test_kastore_version(self):
        version = _tskit.get_kastore_version()
        self.assertEqual(version, (0, 1, 0))
