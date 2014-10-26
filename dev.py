"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import math
import json
import time
import random

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


if __name__ == "__main__":
    #hl_main()
    # ll_main()
    # memory_test()
    large_sim()
    # sort_tree("tmp__NOBACKUP__/large_tree.dat")
    # print_tree("tmp__NOBACKUP__/large_tree.dat")
