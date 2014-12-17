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
import subprocess

import numpy as np
import numpy.random

import _msprime
import msprime

def print_sim(sim):
    print("sample_size = ", sim.get_sample_size())
    print("num_loci = ", sim.get_num_loci())
    print("recombination_rate = ", sim.get_recombination_rate())
    print("random_seed = ", sim.get_random_seed())
    print("time = ", sim.get_time())
    print("tree_file = ", sim.get_tree_file())
    print("num_ancestors = ", sim.get_num_ancestors())
    for segs in sim.get_ancestors():
        print(segs)
    print("population models = ")
    for model in sim.get_population_models():
        print(model)
    # print("X = ", sim.get_X())
    # print("X = ", sim.get_X())
    # print("X = ", sim.get_X())

def ll_main():
    treefile = "tmp__NOBACKUP__/tmp2.dat"
    j = 0
    while True:
        j += 1
        models = [{"type":_msprime.POP_MODEL_CONSTANT, "time":0.3, "size":0.2},
                {"type":_msprime.POP_MODEL_EXPONENTIAL, "time":0.5, "alpha":5}]
        sim = _msprime.Simulator(sample_size=400, random_seed=j,
                tree_file_name=treefile,
                num_loci=1000, recombination_rate=0.1,
                max_memory=1024**3, segment_block_size=10**6,
                population_models=models)
        before = time.time()
        print(sim.run())
        #print(sim.run(0.5))
        duration = time.time() - before
        print("Ran in", duration)
        print_sim(sim)

        before = time.time()
        tr = _msprime.TreeFile(treefile)
        tr.sort()
        duration = time.time() - before
        print("create tree_reader Ran in", duration)
        print(tr.get_num_loci())
        print(tr.get_sample_size())
        print(tr.get_num_trees())
        print(tr.get_metadata())
        s = json.loads(tr.get_metadata())
        # print(s)
        before = time.time()
        j = 0
        for r in tr:
            j += 1
        print(j)

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
    pi, tau = msprime.simulate_tree(4)
    print(pi, tau)
    for l, pi, tau in msprime.simulate_trees(3, 100, 0.1):
        print(l, pi, tau)

def large_sim():
    n = 10**3
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
    ts.set_num_loci(m)
    ts.set_scaled_recombination_rate(4 * N0 * r)
    ts.set_max_memory("245G")
    ts.set_segment_block_size(int(1e8))
    try:
        ts.run()
    except KeyboardInterrupt:
        pass
    except Exception as e:
        print("Died!", e)
    print("num_avl_node_blocks = ", ts.get_num_avl_node_blocks())
    print("num_node_mapping_blocks = ", ts.get_num_node_mapping_blocks())
    print("num_segment_blocks = ", ts.get_num_segment_blocks())

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
    ts.set_num_loci(20)
    ts.set_random_seed(10)
    ts.run()
    msprime.sort_tree_file(treefile)
    tf = msprime.TreeFile(treefile)
    for l, newick in tf.newick_trees():
        print(l, ":", newick)
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
    for l, pi, tau in msprime.TreeFile(tree_file_name):
        print(l, pi, tau)
    for r in msprime.TreeFile(tree_file_name).records():
        print(r)

if __name__ == "__main__":
    # print_tree_file(sys.argv[1])
    # edit_visualisation()
    # mutation_dev()
    example1()
    #hl_main()
    # ll_main()
    # memory_test()
    # large_sim()
    # sort_tree("tmp__NOBACKUP__/large_tree.dat")
    # print_tree("tmp__NOBACKUP__/large_tree.dat")
