"""
Simple client code for development purposes.
"""

from __future__ import print_function
from __future__ import division

import os
import math
import itertools

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
    for tree in tree_sequence.trees():
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
    for sp in ts.trees():
        print(sp)
    ts.dump("example.hdf5")

def diffs_example():
    ts = msprime.simulate(5, 10, scaled_recombination_rate=0.1, random_seed=1)
    for length, records_out, records_in in ts.diffs():
        print(length, records_out, records_in)

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

def large_leaf_count_example():
    n = 1000
    ts = msprime.simulate(n, 100000, scaled_recombination_rate=0.1, random_seed=1)
    num_trees = 0
    potential_nodes = []
    for tree, diff in itertools.izip(ts.trees(True), ts.diffs()):
        num_trees += 1
        _, _, records_in = diff
        # Go through all the nodes that are introduced in this tree and
        # see how many leaves they subtend.
        for node, _, _ in records_in:
            if tree.get_num_leaves(node) == 500:
                potential_nodes.append((node, tree.get_interval()))
    print(potential_nodes)
    print(num_trees)


def count_leaves(pi, n):
    nu = [0 for j in pi]
    for u in range(1, n + 1):
        v = u
        while v != 0:
            nu[v] += 1
            v = pi[v]
    return nu

def leaf_count_example():
    # This is an implementation of the algorithm to count the number
    # of leaves under each node.
    n = 20
    ts = msprime.simulate(n, 10000, scaled_recombination_rate=0.1, random_seed=1)
    pi = [0 for j in range(ts.get_num_nodes() + 1)]
    nu = [0 for j in range(ts.get_num_nodes() + 1)]
    for j in range(1, n + 1):
        nu[j] = 1
    for l, records_out, records_in in ts.diffs():
        for node, children, time in records_out:
            # print("out:", node, children)
            for c in children:
                pi[c] = 0
            leaves_lost = nu[node]
            u = node
            while u != 0:
                nu[u] -= leaves_lost
                # print("Setting nu[", u, "] = ", nu[u])
                u = pi[u]
        for node, children, time in records_in:
            # print("in:", node, children)
            num_leaves = 0
            for c in children:
                pi[c] = node
                num_leaves += nu[c]
            u = node
            while u != 0:
                nu[u] += num_leaves
                # print("setting nu[", u, "] = ", nu[u])
                u = pi[u]
        nup = count_leaves(pi, n)
        assert nup == nu
        # print("pi = ", pi)
        # print("nu = ", nu)
        # print("nup= ", nup)

class LeafSetNode(object):
    """
    Class representing a single leaf in the leaf set.
    """
    def __init__(self, value):
        self.value = value
        self.next = None
        self.prev = None

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return self.__str__()

def leaf_set_example():
    """
    Development for the leaf set algorithm.
    """
    n = 1000

    ts = msprime.simulate(n, 10000, scaled_recombination_rate=0.1, random_seed=1)
    head = [None for j in range(ts.get_num_nodes() + 1)]
    tail = [None for j in range(ts.get_num_nodes() + 1)]
    pi = [0 for j in range(ts.get_num_nodes() + 1)]
    chi = [None for j in range(ts.get_num_nodes() + 1)]
    for j in range(1, n + 1):
        set_node = LeafSetNode(j)
        head[j] = set_node
        tail[j] = set_node
    # print(ts)
    # for st in ts.trees():
    #     print(st)

    for st, (l, records_out, records_in) in itertools.izip(ts.trees(), ts.diffs()):
        for node, (c1, c2), time in records_out:
            print("out:", node, c1, c2)
            pi[c1] = 0
            pi[c2] = 0
            chi[node] = None
            # Break links in the leaf chain.
            tail[c1].next = None
            tail[c2].next = None
            # Clear out head and tail pointers.
            head[node] = None
            tail[node] = None

        for node, (c1, c2), time in records_in:
            print("in:", node, c1, c2)
            pi[c1] = node
            pi[c2] = node
            chi[node] = c1, c2
            u = node
            while u != 0:
                d1, d2 = chi[u]
                print("\tsettting head/tail", u, d1, d2)
                head[u] = head[d1]
                tail[u] = tail[d2]
                u = pi[u]

        for node, _, _ in records_in:
            u = node
            print("propagating", u)
            while u != 0:
                d1, d2 = chi[u]
                print("\tpropogating", u, d1, d2)
                print("\tbefore tail[d1] = ", tail[d1], "next = ", tail[d1].next)
                if tail[d1] is None:
                    break
                tail[d1].next = head[d2]
                print("\tafter  tail[d1] = ", tail[d1], "next = ", tail[d1].next)
                u = pi[u]


        print("node", "pi", "chi", "head", "tail", sep="\t")
        for j in range(ts.get_num_nodes() + 1):
            if st.get_parent(j) != 0 or j == st.get_root():
                print(
                    j, st.get_parent(j), st.get_children(j),
                    head[j], tail[j], sep="\t")
        # for u in st.nodes():
        for u in [st.get_root()]:
            leaves = list(st.leaves(u))
            leaf_list = []
            v = head[u]
            while v is not tail[u]:
                # print("\t", v.value, v.next)
                leaf_list.append(v.value)
                assert len(leaf_list) < n
                v = v.next
            leaf_list.append(v.value)
            # print(u, leaf_list, leaves, sep="\t")
            if leaf_list != leaves:
                print(leaf_list)
                print(leaves)
            assert leaf_list == leaves
        # print()

def allele_frequency_example():
    # n = 10000
    # ts = msprime.simulate(
    #     n, 100000, scaled_recombination_rate=0.1, scaled_mutation_rate=0.1,
    #     random_seed=1)
    ts = msprime.load("tmp__NOBACKUP__/gqt.hdf5")
    n = ts.get_sample_size()
    num_mutations = 0
    min_frequency = 0.0001
    num_trees = 0
    for tree in ts.trees():
        num_trees += 1
        for pos, node in tree.mutations():
            if tree.get_num_leaves(node) / n < min_frequency:
                num_mutations += 1
    print("num_mutatinos = ", num_mutations, "\t", num_mutations / ts.get_num_mutations())
    print("total_mutations = ", ts.get_num_mutations())
    print("num_trees = ", num_trees)

def write_ped(ts, num_cases, prefix):
    import numpy as np
    # Write the ped file.
    s = ts.get_num_mutations()
    indices = np.arange(0, s)
    with open(prefix + ".ped", "w") as f:
        for j, h in enumerate(ts.haplotypes(), 1):
            family_id = j
            ind_id = j
            paternal_id = 0
            maternal_id = 0
            sex = 1
            # 2 == case, 1 == control
            phenotype = 2 if j <= num_cases else 1
            print(
                family_id, ind_id, paternal_id, maternal_id, sex, phenotype,
                end=" ", file=f)
            # Neat trick for converting numerical strings to numpy array:
            # interpret data as raw chars, and subtrace ascii code.
            nph = np.fromstring(h, np.int8) - ord('0')
            # Plink genotypes are 1/2, not 0/1
            nph += 1
            genotypes = np.zeros(4 * s, dtype="S1")
            genotypes[indices * 4] = nph
            genotypes[indices * 4 + 1] = " "
            genotypes[indices * 4 + 2] = nph
            genotypes[indices * 4 + 3] = " "
            genotypes.tofile(f)
            # Remove trailing space
            f.seek(-1, os.SEEK_CUR)
            print(file=f)
    # Make a map file. We're using haploid data, so we put it on the
    # male X chromosome.
    with open(prefix + ".map", "w") as f:
        for j in range(ts.get_num_mutations()):
            print("X", "rs{}".format(j), 0, j + 1, file=f)



def gwas_example():
    # n = 100
    # ts = msprime.simulate(
    #     n, 1000, scaled_recombination_rate=0.1, scaled_mutation_rate=10,
    #     random_seed=1)
    ts = msprime.load("tmp__NOBACKUP__/gqt.hdf5")
    n = ts.get_sample_size()
    num_cases = n // 2
    # write_ped(ts, num_cases, "tmp__NOBACKUP__/test")
    # write_ped(ts, num_cases, "tmp__NOBACKUP__/plink/gqt")
    write_plink_assoc(ts, num_cases)

def write_plink_assoc(ts, num_cases):
    site = 0
    n = ts.get_sample_size()
    cases = list(range(1, num_cases + 1))
    for tree in ts.trees(cases):
        for pos, node in tree.mutations():
            num_leaves = tree.get_num_leaves(node)
            cases_with_mut = tree.get_num_tracked_leaves(node)
            controls_with_mut = tree.get_num_leaves(node) - cases_with_mut
            f_cases = cases_with_mut / num_cases
            f_controls = controls_with_mut / (n - num_cases)
            if num_leaves >= n / 2:
                # The mutation is the major allele
                a1 = 1
                a2 = 2
                fa = 1 - f_cases
                fu = 1 - f_controls
            else:
                # The mutation is the minor allele
                a1 = 2
                a2 = 1
                fa = f_cases
                fu = f_controls
            case_odds = fa / (1 - fa)
            control_odds = fu / (1 - fu)
            odds_ratio = "NA" if control_odds == 0 else case_odds / control_odds
            print(site + 1, a1, fa, fu, a2, odds_ratio, sep="\t")

def convert_haplotypes():
    import sys
    ts = msprime.load(sys.argv[1])
    ts.generate_mutations(0.0001, 1)
    print("Generated mutations", ts.get_num_mutations())
    c = 0
    for h in ts.haplotypes():
        c += 1
        # print(h)
    print("generated ", c, "haplotypes")

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
    # draw_trees()
    # leaf_count_example()
    # large_leaf_count_example()
    # diffs_example()
    # leaf_set_example()
    # allele_frequency_example()
    # gwas_example()
    convert_haplotypes()
