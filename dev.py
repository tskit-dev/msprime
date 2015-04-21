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
            # print(sim.run(t))
            duration = time.time() - before
            # print("Ran in", duration)
            # print_sim(sim)
            segs = sim.get_ancestors()
            crs = sim.get_coalescence_records()
            bps = sim.get_breakpoints()
            # print(len(segs), len(crs), len(bps))
            check_sim(sim)

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

def hl_main():

    random.seed(1)
    # pi, tau = msprime.simulate_tree(4)
    # print(pi, tau)
    # for l, pi, tau in msprime.simulate_trees(3, 100, 0.1):
    #     print(l, pi, tau)
    treefile = "tmp__NOBACKUP__/tmp.hdf5"
    sim = msprime.TreeSimulator(10)
    sim.set_random_seed(1)
    sim.set_num_loci(1000)
    sim.set_scaled_recombination_rate(0.1)
    models = [
            msprime.ExponentialPopulationModel(0, 1),
            msprime.ExponentialPopulationModel(0.1, 2),
            msprime.ConstantPopulationModel(0.5, 2.0),
    ]
    for m in models:
        sim.add_population_model(m)

    tree_sequence = sim.run()
    # tree_sequence.dump(treefile)
    # tree_sequence.print_state()
    for newick in tree_sequence.newick_trees():
        print(newick)
    # for l, records_out, records_in in tree_sequence.diffs():
    #     print("\tin  = ", records_in)
    #     print("\tout = ", records_out)

    # ts = msprime.TreeSequence.load(treefile)
    # print("after")
    # ts.print_state()

def print_newick(filename):
    ts = msprime.TreeSequence.load(filename)
    j = 0
    s = 0
    for l, newick in ts.newick_trees():
        j += 1
        s += l
        # newick.tofile(sys.stdout)
        # print()
        # print(l, newick)
        if j % 100000 == 0:
            print(j, s)
            break

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

def print_tree_records(treefile):
    tf = msprime.TreeFile(treefile)
    for l, r, c1, c2, p, t in tf.records():
        print(l, r, c1, c2, p, t, sep="\t")

def sort_tree(f):
    tf = msprime.TreeFile(f)
    print(tf.issorted())
    msprime.sort_tree_file(f)

def print_tree(f):
    tf = msprime.TreeFile(f)
    print(tf.issorted())
    print(tf.get_metadata())
    j = 0
    t = 0
    for l, pi, tau in tf:
        j += 1
        t += l
        # print(l, pi, tau)
    print(j)
    print(t)


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

def print_tree_file(tree_file_name):
    tf = msprime.TreeFile(tree_file_name)
    if not tf.issorted():
        msprime.sort_tree_file(tree_file_name)
        tf = msprime.TreeFile(tree_file_name)
    metadata = json.loads(tf.get_metadata())
    # print(metadata)
    n = metadata["sample_size"]
    # for l, r, c1, c2, p, t in tf.records():
    #     print(l, r, c1, c2, p, t, sep="\t")
    count = 0
    # for l, pi, tau in tf.sparse_trees():
    for l, r, c1, c2, p, t in tf.records():
        # print(l, pi, tau)
        # print(l, len(pi))
        count += 1
    print(count)

    # for r in msprime.TreeFile(tree_file_name).records():
    #     print(r)

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
    # print_tree_file(sys.argv[1])
    # edit_visualisation()
    # mutation_dev()
    # example1()
    # hl_main()
    # ll_main()
    print_newick(sys.argv[1])
    # memory_test()
    # large_sim()
    # print_tree_records(sys.argv[1])
    # sort_tree("tmp__NOBACKUP__/large_tree.dat")
    # print_tree("tmp__NOBACKUP__/large_tree.dat")

