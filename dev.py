"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import os
import sys
import math
import json
import time
import random
import tempfile
import unittest
import collections
import subprocess

import numpy as np
import numpy.random

import statsmodels.api as sm

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot

import _msprime
import msprime

def print_sim(sim):
    print("sample_size = ", sim.get_sample_size())
    print("num_loci = ", sim.get_num_loci())
    print("recombination_rate = ", sim.get_scaled_recombination_rate())
    print("random_seed = ", sim.get_random_seed())
    print("time = ", sim.get_time())
    print("num_ancestors = ", sim.get_num_ancestors())
    print("num_breakpoints = ", sim.get_num_breakpoints())
    print("num_coalescence_records = ", sim.get_num_coalescence_records())
    print("num_coalescence_record_blocks = ", sim.get_num_coalescence_record_blocks())
    for segs in sim.get_ancestors():
        print(segs)
    print("breakpoints = ", sim.get_breakpoints())
    for cr in sim.get_coalescence_records():
        print(cr)
    print("population models = ")
    for model in sim.get_population_models():
        print(model)

def check_sim(sim):
    n = sim.get_sample_size()
    ancestors = sim.get_ancestors()
    records = sim.get_coalescence_records()
    breakpoints = sim.get_breakpoints()
    # The amount of ancestral material in the coalescence records and
    # the extant segments over all intervals should be n.
    segments_am = [0 for b in breakpoints[:-1]]
    for ind in ancestors:
        for l, r, _ in ind:
            j = breakpoints.index(l)
            while breakpoints[j] < r:
                segments_am[j] += 1
                j += 1
    records_am = [0 for b in breakpoints[:-1]]
    for l, r, _, _, _ in zip(*records):
        j = breakpoints.index(l)
        while breakpoints[j] < r:
            records_am[j] += 1
            j += 1
    for segment_am, record_am in zip(segments_am, records_am):
        if segment_am == 0:
            assert record_am == n - 1
        else:
            assert segment_am + record_am == n


def ll_main():

    j = 0
    while True:
    # if True:
        j += 1
        models = [{"type":_msprime.POP_MODEL_CONSTANT, "start_time":0.3, "size":0.2},
                {"type":_msprime.POP_MODEL_EXPONENTIAL, "start_time":0.5, "alpha":5}]
        sim = _msprime.Simulator(sample_size=4, random_seed=j,
                num_loci=1000, scaled_recombination_rate=0.1,
                max_memory=1024**3, segment_block_size=10**6,
                coalescence_record_block_size=1000,
                population_models=models)
        for t in [0.001, 0.1, 0.15, 0.5, 1000]:
            before = time.time()
            # print(sim.run())
            sim.run(t)
            duration = time.time() - before
            # print("Ran in", duration)
            # print_sim(sim)
            segs = sim.get_ancestors()
            crs = sim.get_coalescence_records()
            c = 0
            for r in crs:
                c += 1
            bps = sim.get_breakpoints()
            # print(len(segs), len(crs), len(bps))
            # check_sim(sim)
        # crs  sim.get_coalescence_records()
        ts = _msprime.TreeSequence()
        ts.create(sim)
        iterator = _msprime.TreeDiffIterator(ts, True)
        c = 0
        for r in iterator:
            c += 1
        iterator = _msprime.TreeDiffIterator(ts, False)
        del ts
        c = 0
        for r in iterator:
            c += 1
            # print(r)
        # assert c == len(crs[0])
        # # print(iterator)
        ts = _msprime.TreeSequence()
        ts.create(sim)
        nc = _msprime.NewickConverter(ts, 4, False)
        c = 0
        for l, t in nc:
            c += l + len(t)
        nc = _msprime.NewickConverter(ts, 4, False)
        x = next(nc)
        files = ["/", "/nofile", "/" + "x" * 4192]
        for f in files:
            try:
                ts.dump(f)
            except _msprime.LibraryError as e:
                pass
        # f = "tmp__NOBACKUP__/dev.hdf5"
        # ts.dump(f)
        # ts2 = _msprime.TreeSequence()
        # ts2.load(f)
        # records1 = [ts.get_record(j) for j in range(ts.get_num_records())]
        # records2 = [ts2.get_record(j) for j in range(ts2.get_num_records())]
        # assert records1 == records2
        # # ts.generate_mutations(10, 1)
        # f = "tmp__NOBACKUP__/dev-mutaions.hdf5"
        # ts.dump(f)
        # n = ts.get_sample_size()
        # hg = _msprime.HaplotypeGenerator(ts)
        # h1 = [hg.get_haplotype(j) for j in range(1, n + 1)]
        # # ts2 = _msprime.TreeSequence()
        # # ts2.load(f)
        # ts2 = ts
        # hg = _msprime.HaplotypeGenerator(ts2)
        # h2 = [hg.get_haplotype(j) for j in range(1, n + 1)]
        # assert h1 == h2

        # for r in ts.records():
        #     print(r)
        # for length, nodes_out, nodes_in in ts.

        # for j in range(tr.get_num_trees()):
        #     breakpoint, pi, tau = tr.get_tree(j)
        #     #print(j, breakpoint)
        # duration = time.time() - before
        # print("tree loop Ran in", duration)
        # print(sim.run())
        # print_sim(sim)

    # tv = _msprime.TreeViewer(treefile)
    # for length, pi, tau in tv:
    #     print(length, pi, tau)

def new_sparse_trees(tree_sequence):
    """
    prototype implementation of list-based sparse trees.
    """
    # TODO add get num node to tree sequence
    n = tree_sequence.get_sample_size()
    R = max(node for _, _, node, _, _ in tree_sequence.records())
    tau  = [0.0 for j in range(R + 1)]
    pi = [0 for j in range(R + 1)]
    position = 0
    records = list(tree_sequence.records())
    right_sorted = sorted(records, key=lambda x: x[1])
    j = 0
    # for t1, t2 in zip(records, right_sorted):
    #     print(t1[:-1], "|", t2[:-1], sep="\t")
    for j in range(n - 1):
        l, r, node, children, time = records[j]
        assert l == 0
        tau[node] = time
        for c in children:
            pi[c] = node
    position = right_sorted[0][1]
    yield position, pi, tau
    j += 1
    while j < len(records):
        # print("j = ", j, j - n + 1)
        k = j - n + 1
        while right_sorted[k][1] == position:
            # print("Removing", right_sorted[k])
            for c in right_sorted[k][3]:
                pi[c] = 0
            k += 1
        while j < len(records) and records[j][0] == position:

            _, _, node, children, time = records[j]
            # print("inserting ", records[j])
            tau[node] = time
            for c in children:
                pi[c] = node
            j += 1
        yield right_sorted[k][1] - position, pi, tau
        position = right_sorted[k][1]
        # print("Tree @ position", position)
        # for x in range(1, R + 1):
        #     print(x, "\t", pi[x])
        # print("yielding", pi)
        # print("unsetting", j - n, right_sorted[j - n])
        # for c in right_sorted[j - n][3]:
        #     pi[c] = 0


def hl_main():

    random.seed(1)
    # pi, tau = msprime.simulate_tree(4)
    # print(pi, tau)
    # for l, pi, tau in msprime.simulate_trees(3, 100, 0.1):
    #     print(l, pi, tau)
    treefile = "tmp__NOBACKUP__/tmp.hdf5"
    n = 3000
    sim = msprime.TreeSimulator(n)
    sim.set_random_seed(1)
    sim.set_num_loci(1000)
    sim.set_scaled_recombination_rate(0.1)
    sim.set_max_memory("10G")

    models = [
            msprime.ExponentialPopulationModel(0, 1),
            msprime.ExponentialPopulationModel(0.1, 2),
            msprime.ConstantPopulationModel(0.5, 2.0),
    ]
    for m in models:
        sim.add_population_model(m)

    tree_sequence = sim.run()
    st = tree_sequence.sparse_trees()
    for l, pi, tau in new_sparse_trees(tree_sequence):
        l2, pi2, tau2 = next(st)
        for j in range(1, n + 1):
            path1, path2 = [], []
            u = j
            while u != 0:
                path1.append(u)
                u = pi[u]
            u = j
            while u != 0:
                path2.append(u)
                u = pi2[u]
            # print("new", path1)
            # print("old", path2)
            assert path1 == path2
            if l != l2:
                print(l, l2)
            assert l == l2


      # for _, pi, _ in tree_sequence.sparse_trees():
    #     print(pi)
    #     path = []
    #     j = 1
    #     while j != 0:
    #         path.append(j)
    #         j = pi[j]
    #     print("\t", path)
    # print()
    # for j in range(1, n + 1):
    #     print(j, ":", tree_sequence.sample_nodes(j))

    # hg = msprime.HaplotypeGenerator(tree_sequence, 116.1,
    #         random_seed=1)
    # for h in hg.haplotypes():
    #     print(h)
    # print()
    # print("locations:", hg.get_locations())
    # # hg = msprime.CHaplotypeGenerator(tree_sequence, 116.1,
    # #         random_seed=1)
    # # for h in hg.haplotypes():
    # #     print(h)

    """
    N = 10000
    s = np.zeros(N)
    cs = np.zeros(N)
    for j in range(N):
        hg = msprime.CHaplotypeGenerator(tree_sequence, 116.1,
                random_seed=j)
        cs[j] = hg.get_num_segregating_sites()
        hg = msprime.HaplotypeGenerator(tree_sequence, 116.1,
                random_seed=j)
        s[j] = hg.get_num_segregating_sites()
    print("mean:", np.mean(cs), np.mean(s))
    print("var :", np.var(cs), np.var(s))
    sm.graphics.qqplot(s)
    sm.qqplot_2samples(s, cs, line="45")
    f = "tmp__NOBACKUP__/s.png"
    pyplot.savefig(f, dpi=72)
    pyplot.clf()
    """

def analyse_records(treefile):
    tree_sequence = msprime.TreeSequence.load(treefile)
    n = tree_sequence.get_sample_size()
    num_records = []
    iterator = tree_sequence.diffs()
    length, r_out, r_in = next(iterator)
    parents = {}
    for children, parent, time in r_in:
        for j in children:
            parents[j] = parent
    assert len(r_out) == 0

    num_trees = 0
    num_leaf_diffs = 0
    for length, r_out, r_in in iterator:
        num_trees  += 1
        k = len(r_in)
        root = 1
        while root in parents:
            root = parents[root]
        print("TREE", len(r_in))
        print("ROOT", root)
        print("OUT:")
        analyse_subtrees(r_out)
        for children, parent, time in r_out:
            v = parent
            d = 0
            while v in parents:
                v = parents[v]
                d += 1
            leaves = 0
            for j in children:
                if j <= n:
                    leaves += 1
            print("\t", d, "\t", leaves, "\t", children, parent, time)
            if k == 1:

                assert parent == root
        if leaves > 0:
            num_leaf_diffs += 1
        for children, parent, time in r_out:
            for j in children:
                del parents[j]
        print("IN:")
        analyse_subtrees(r_in)
        for children, parent, time in r_in:
            for j in children:
                parents[j] = parent
            print("\t", children, parent, time)
        num_records.append(len(r_in))

    print(num_trees, num_leaf_diffs / num_trees)
    pyplot.hist(num_records, range(1, 10), normed=True)
    f = "tmp__NOBACKUP__/num_records.png"
    pyplot.savefig(f, dpi=72)
    pyplot.clf()

def analyse_subtrees(records):

    parents = {}
    internal_nodes = set()
    for children, parent, time in records:
        internal_nodes.add(parent)
        for c in children:
            parents[c] = parent
    leaves = set()
    for children, parent, time in records:
        for c in children:
            if c not in internal_nodes:
                leaves.add(c)
    subtrees = collections.defaultdict(list)
    for leaf in leaves:
        v = leaf
        while v in parents:
            v = parents[v]
        subtrees[v].append(leaf)
    print("subtrees:", subtrees)
    print("internal_nodes :", internal_nodes)
    print("leaves = ", leaves)
    print("roots = ", list(subtrees.keys()))


def print_newick(filename):
    ts = msprime.TreeSequence.load(filename)
    j = 0
    s = 0
    for l, newick in ts.newick_trees(10):
        j += 1
        s += l
        # # print()
        print(l, newick)
        # if j % 10000 == 0:
        #     print(j, s)
        #     break

def large_sim():
    n = 10**3
    # m = 20 * 10**6
    m = 20 * 10**6
    r = 1e-8
    # Calculate the models
    N0 = 1e6
    N1 = 2e4
    N2 = 2e3
    N3 = 2e4
    g1 = 400
    g2 = 2000
    t1 = g1 / (4 * N0)
    t2 = g2 / (4 * N0)
    # Calculate the growth rates.
    alpha1 = -math.log(N1 / N0) / t1
    alpha2 = -math.log(N2 / N1) / (t2 - t1)
    models = [
            msprime.ExponentialPopulationModel(0, alpha1),
            msprime.ExponentialPopulationModel(t1, alpha2),
            msprime.ConstantPopulationModel(t2, N3 / N0),
    ]
    treefile = "tmp__NOBACKUP__/large_tree_models.dat"
    ts = msprime.TreeSimulator(n, treefile)
    ts.set_random_seed(1)
    ts.set_num_loci(m)
    ts.set_scaled_recombination_rate(4 * N0 * r)
    ts.set_squash_records(True)
    ts.set_max_memory("245G")
    # ts.set_segment_block_size(int(1e8))
    for m in models:
        ts.add_population_model(m)
    try:
        ts.run()
    except KeyboardInterrupt:
        pass
    except Exception as e:
        print("Died!", e)
    print("num_avl_node_blocks = ", ts.get_num_avl_node_blocks())
    print("num_node_mapping_blocks = ", ts.get_num_node_mapping_blocks())
    print("num_segment_blocks = ", ts.get_num_segment_blocks())

    # tf = msprime.TreeFile(treefile)
    # for l, r, c1, c2, p, t in tf.records():
    #     print(l, r, c1, c2, p, t, sep="\t")
    # msprime.sort_tree_file()

def memory_test():
    while True:
        n = random.randint(2, 100)
        m = random.randint(10, 1000)
        r = random.random()
        treefile = "tmp__NOBACKUP__/treefile.dat"
        # TODO add some different n and m values here.
        ts = msprime.TreeSimulator(n, treefile)
        # todo verify all the setters.
        # self.assertEqual(ts.get_sample_size(), n)
        ts.set_scaled_recombination_rate(r)
        ts.set_num_loci(m)
        ts.run()
        msprime.sort_tree_file(treefile)
        tf = msprime.TreeFile(treefile)
        l = [t for t in tf]


def example1():
    treefile = "tmp__NOBACKUP__/treefile.dat"
    # TODO add some different n and m values here.
    ts = msprime.TreeSimulator(5, treefile)
    ts.set_scaled_recombination_rate(0.1)
    ts.set_num_loci(60)
    ts.set_random_seed(10)
    ts.set_squash_records(True)
    ts.run()
    msprime.sort_tree_file(treefile)
    tf = msprime.TreeFile(treefile)
    for l, pi, tau in tf.sparse_trees():
        print(l, pi, tau)
    # for l, records_in, records_out in tf.get_tree_diffs():
    #     print(l, records_in, records_out, sep="\t")

    # pi = {}
    # tau = {}
    # for l, records_in, records_out in tf.get_tree_diffs():
    #     for c1, c2, p in records_out:
    #         del pi[c1]
    #         del pi[c2]
    #         del tau[p]
    #     for c1, c2, p, t in records_in:
    #         pi[c1] = p
    #         pi[c2] = p
    #         tau[p] = t
    #     print(l, pi, tau, sep="\t")
    # print()
    # tf = msprime.TreeFile(treefile)
    # for l, pi, tau in tf.trees():
    #     print(l, pi, tau, sep="\t")
    # # for l, newick in tf.newick_trees():
    #     print(l, ":", newick)
    # models = [
    #     msprime.ConstantPopulationModel(0.5, 2.0),
    #     msprime.ConstantPopulationModel(1.0, 0.5)
    # ]
    # pi, tau = msprime.simulate_tree(10, models)

def isdescendent(u, v, pi):
    """
    Returns True if the node u is a descendent of the node v in the
    specified oriented forest pi.
    """
    j = u
    while j != v and j != 0:
        j = pi[j]
    return j == v

def mutation_dev():
    random.seed(1)
    np.random.seed(1)
    n = 4
    m = 100
    r = 0.1
    mu = 2.0
    treefile = "tmp__NOBACKUP__/treefile.dat"
    ts = msprime.TreeSimulator(n, treefile)
    ts.set_scaled_recombination_rate(r)
    ts.set_num_loci(m)
    ts.run()
    msprime.sort_tree_file(treefile)
    tf = msprime.TreeFile(treefile)
    pi = [0 for j in range(2 * n)]
    tau = [0 for j in range(2 * n)]
    b = 1
    sequences = ["" for j in range(n + 1)]
    # This algorithm basically works, we just need to actually prove
    # it. It'll all hang on the fact that records are sorted in time
    # order, so that we never get the wrong time in between trees.
    for l, c1, c2, p, t in tf.records():
        print(l, c1, c2, p, t)
        for c in [c1, c2]:
            k = np.random.poisson(mu * (t - tau[c]))
            if k > 0:
                print(k, "mutations happened on ", c, "->", p)
                for j in range(1, n + 1):
                    v = str(int(isdescendent(j, c, pi)))
                    sequences[j] += v

        if l != b:
            print("new tree:",  l - b, pi, tau)
            b = l
        pi[c1] = p
        pi[c2] = p
        tau[p] = t

    for s in sequences:
        print(s)

def edit_visualisation():
    """
    Create some visualisation of the edit structure of the tree file.
    """
    random.seed(1)
    np.random.seed(1)
    n = 5
    m = 10
    r = 0.5
    mu = 2.0
    treefile = "tmp__NOBACKUP__/treefile.dat"
    ts = msprime.TreeSimulator(n, treefile)
    ts.set_scaled_recombination_rate(r)
    ts.set_num_loci(m)
    ts.run()
    msprime.sort_tree_file(treefile)
    tf = msprime.TreeFile(treefile)
    max_time = 0
    for l, pi, tau in tf:
        max_time = max(max_time, *tau)
    tf = msprime.TreeFile(treefile)
    for j, (l, pi, tau) in enumerate(tf):
        output_file = "tmp__NOBACKUP__/viz_{0}.pdf".format(j)
        fd, filename = tempfile.mkstemp(".asy")
        out = os.fdopen(fd, "w")
        make_tree_visualisation(pi, tau, out, max_time)
        out.close()
        args = ["asy", filename, "-f", "pdf", "-o", output_file]
        subprocess.check_call(args)
        os.unlink(filename)
    tf = msprime.TreeFile(treefile)
    for r in tf.records():
        print(r)

def make_tree_visualisation(pi, tau, out, max_time):
    print(pi, tau, max_time)
    print("size(8cm, 0);", file=out)
    print("defaultpen(fontsize(6pt));", file=out)
    n = len(pi) // 2
    nodes = {}
    for j in range(1, n + 1):
        x = (j, 0)
        nodes[j] = x
    for j in range(n + 1, 2 * n):
        k = j - n
        if j % 2 == 0:
            x0 = j % n + 0.5 * k
        else:
            x0 = n - (j % n - 1) + 0.5 * k
        x = (x0, k)
        nodes[j] = x
    for j, x in nodes.items():
        s = 'dot({0});'.format(x)
        l = "S" if j <= n else "NW"
        print(s, file=out)
        s = 'label("{0}, t={3:.3f}", {1}, {2});'.format(j, x, l, tau[j])
        print(s, file=out)
    m = [0 for j in range(2 * n)]
    for j in range(1, n + 1):
        u = j
        while pi[u] != 0 and m[u] == 0:
            m[u] = 1
            x = nodes[u]
            y = nodes[pi[u]]
            s = 'draw({0}--{1});'.format(x, y)
            print(s, file=out)
            u = pi[u]

def print_tree_file(filename):
    ts = msprime.TreeSequence.load(filename)
    max_node = 0
    for l, r, node, (c1, c2), t in ts.records():
        if node > max_node:
            max_node = node
        # print(l, r, node, c1, c2, t, sep="\t")
    print(max_node)

def dump_simulation(filename, n=10, m=100):
    sim = msprime.TreeSimulator(n)
    sim.set_num_loci(m)
    sim.set_scaled_recombination_rate(1)
    ts = sim.run()
    for l, r, node, children, time in ts.records():
        print(l, r, node, children, time, sep="\t")
    # ts.dump(filename)


class VerifyTrees(unittest.TestCase):

    def verify_sparse_tree(self, n, pi, tau):
        """
        Verifies that the specified sparse_tree is a consistent coalescent history
        for a sample of size n.
        """
        self.assertEqual(set(pi.keys()), set(tau.keys()))
        self.assertEqual(len(pi), 2 * n - 1)
        self.assertEqual(len(tau), 2 * n - 1)
        # Zero should not be a node
        self.assertNotIn(0, pi)
        self.assertNotIn(0, tau)
        # 1 to n inclusive should always be nodes
        for j in range(1, n + 1):
            self.assertIn(j, pi)
            self.assertIn(j, tau)
        num_children = collections.defaultdict(int)
        roots = 0
        for j in pi.keys():
            num_children[pi[j]] += 1
            roots += pi[j] == 0
        if roots != 1:
            print(pi)
        self.assertEqual(roots, 1)
        # nodes 1 to n are leaves.
        for j in range(1, n + 1):
            self.assertNotEqual(pi[j], 0)
            self.assertEqual(tau[j], 0)
            self.assertEqual(num_children[j], 0)
        # All non-leaf nodes should be binary with non-zero times.
        taup = {}
        for j in pi.keys():
            if j > n:
                self.assertEqual(num_children[j], 2)
                self.assertGreater(tau[j], 0.0)
                taup[j] = tau[j]
        # times of non leaves should be distinct
        self.assertEqual(len(set(taup)), len(taup))
        self.verify_general_tree(n, pi, tau)

    def verify_general_tree(self, n, pi, tau):
        """
        Verifies properties that should hold for both sparse and
        dense trees.
        """
        # Times of leaves should be zero, and increasing up the tree
        for j in range(1, n + 1):
            self.assertEqual(tau[j], 0.0)
            last_time = -1
            k = j
            while k != 0:
                self.assertNotEqual(k, pi[k])
                self.assertGreater(tau[k], last_time)
                last_time = tau[k]
                k = pi[k]

    def test_large_tree(self):
        tree_file_name = "tmp__NOBACKUP__/human-treefile-squashed-5-5.dat"
        tf = msprime.TreeFile(tree_file_name)
        if not tf.issorted():
            msprime.sort_tree_file(tree_file_name)
            tf = msprime.TreeFile(tree_file_name)
        metadata = json.loads(tf.get_metadata())
        n = metadata["sample_size"]
        count = 0
        for l, pi, tau in tf.sparse_trees():
            count += 1
            self.verify_sparse_tree(n, pi, tau)
            if count % 10 == 0:
                print(".", end="")
                sys.stdout.flush()
        print()
        print("Checked", count, "trees")


if __name__ == "__main__":
    # unittest.main()
    # dump_simulation(sys.argv[1])
    # print_tree_file(sys.argv[1])
    # analyse_records(sys.argv[1])
    # edit_visualisation()
    # mutation_dev()
    # example1()
    # hl_main()
    ll_main()
    # print_newick(sys.argv[1])
    # memory_test()
    # large_sim()

