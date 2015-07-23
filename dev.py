"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import msprime


def haplotype_example():
    tree_sequence = msprime.simulate(
        sample_size=10, num_loci=1000, scaled_recombination_rate=0.1,
        scaled_mutation_rate=0.01, random_seed=1)
    for h in tree_sequence.haplotypes():
        print(h)


def dump_example():
    tree_sequence = msprime.simulate(
        sample_size=10, num_loci=1000, scaled_recombination_rate=0.1,
        scaled_mutation_rate=0.01, random_seed=1)
    haplotypes = list(tree_sequence.haplotypes())
    tree_sequence.dump("example.hdf5")
    # Now, load another tree sequence instance from this file
    other_tree_sequence = msprime.load("example.hdf5")
    other_haplotypes = list(other_tree_sequence.haplotypes())
    assert haplotypes == other_haplotypes

def newick_example():
    tree_sequence = msprime.load("example.hdf5")
    with open("example.newick", "w") as f:
        iterator = tree_sequence.newick_trees(8)
        for l, ns in iterator:
            print("[{0}]".format(l), end="", file=f)
            print(ns, file=f)

def tree_example():
    tree_sequence = msprime.simulate(
        sample_size=10, num_loci=1000, scaled_recombination_rate=0.1,
        random_seed=1)
    for tree in tree_sequence.sparse_trees():
        print(tree)


def dump_file(filename):
    tree_sequence = msprime.load(filename)
    for r in tree_sequence.records():
        print(r)

class AnonymousObject():
    pass

def generate_trees(tree_sequence):
    N = tree_sequence.get_num_nodes() + 1
    M = tree_sequence.get_num_records()
    s = AnonymousObject()
    s.left = [0 for _ in range(M)]
    s.right = [0 for _ in range(M)]
    s.node = [0 for _ in range(M)]
    s.children = [None for _ in range(M)]
    s.time = [0 for _ in range(M)]
    for j, record in enumerate(tree_sequence.records()):
        s.left[j] = record[0]
        s.right[j] = record[1]
        s.node[j] = record[2]
        s.children[j] = record[3]
        s.time[j] = record[4]
    L = sorted(range(M), key=lambda j: (s.left[j], s.time[j]))
    R = sorted(range(M), key=lambda j: (s.right[j], -s.time[j]))
    t = AnonymousObject()
    t.parent = [0 for _ in range(N)]
    t.time = [0 for _ in range(N)]
    t.children = [None for _ in range(N)]
    t.left = 0
    t.right = s.right[R[0]]
    t.root = 0
    j = 0
    k = 0
    while True:
        while j < M and s.left[L[j]] == t.left:
            v = L[j]
            t.parent[s.children[v][0]] = s.node[v]
            t.parent[s.children[v][1]] = s.node[v]
            t.time[s.node[v]] = s.time[v]
            t.children[s.node[v]] = s.children[v]
            if s.node[v] > t.root:
                t.root = s.node[v]
            j += 1
        yield t
        if j == M:
            break
        while s.right[R[k]] == t.right:
            v = R[k]
            t.parent[s.children[v][0]] = 0
            t.parent[s.children[v][1]] = 0
            t.children[s.node[v]] = None
            t.time[s.node[v]] = 0
            if s.node[v] == t.root:
                t.root = max(s.children[v])
            k += 1
        t.left = t.right
        t.right = s.right[R[k]]

def tree_algorithm():
    tree_sequence = msprime.simulate(
        sample_size=40, num_loci=1000, scaled_recombination_rate=0.1,
        random_seed=1)
    for t1, t2 in zip(generate_trees(tree_sequence), tree_sequence.sparse_trees()):
        assert t1.parent == list(t2.parent)
        assert t1.time == list(t2.time)
        assert t1.left == t2.left
        assert t1.right == t2.right
        assert t1.root == t2.root

def large_example():
    msprime.simulate(10**6, 100 * 10**6, scaled_recombination_rate=0.001,
            random_seed=1, max_memory="10G")

def small_example():
    ts = msprime.simulate(5, 10, scaled_recombination_rate=0.1, random_seed=1)
    for sp in ts.sparse_trees():
        print(sp)
    ts.dump("example.hdf5")


def draw_tree():
    tree = msprime.simulate_tree(5, random_seed=1, scaled_mutation_rate=0.1)
    print(list(tree.mutations()))
    print(list(tree.nodes()))
    tree.draw("example.svg")

def draw_trees():
    # tree = msprime.simulate_trees(5, random_seed=1, scaled_mutation_rate=0.1)
    tree_sequence = msprime.simulate(
        10, 100, scaled_recombination_rate=1, random_seed=1)
    for j, tree in enumerate(tree_sequence.trees()):
        tree.draw("tmp__NOBACKUP__/example-{}.svg".format(j))


if __name__ == "__main__":
    # haplotype_example()
    # dump_example()
    # newick_example()
    # tree_example()
    # dump_file("example.hdf5")
    # tree_algorithm()
    # large_example()
    # draw_tree()
    # small_example()
    draw_trees()
