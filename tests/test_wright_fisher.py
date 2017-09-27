#
# Copyright (C) 2017 University of Oxford
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
Test various functions using messy tables output by a forwards-time simulator.
"""
from __future__ import print_function
from __future__ import division

import itertools
import random
import unittest

import numpy as np
import math

import msprime


def random_breakpoint():
    return min(1.0, max(0.0, 2 * random.random() - 0.5))


def random_mutations(rate, infinite_sites=True):
    nmuts = np.random.poisson(lam=rate)
    if not infinite_sites:
        muts = [math.floor(100*random.random())/100 for _ in range(nmuts)]
    else:
        muts = [random.random() for _ in range(nmuts)]
    return muts


def random_allele():
    return random.choice(['A', 'C', 'G', 'T'])


def wf_sim(
        N, ngens, survival=0.0, mutation_rate=0.0, deep_history=True, debug=False,
        seed=None):
    """
    SIMPLE simulation of a bisexual, haploid Wright-Fisher population of size N
    for ngens generations, in which each individual survives with probability
    survival and only those who die are replaced.  The chromosome is 1.0
    Morgans long, and the mutation rate is in units of mutations/Morgan/generation.
    """
    if seed is not None:
        random.seed(seed)
    nodes = msprime.NodeTable()
    edgesets = msprime.EdgesetTable()
    migrations = msprime.MigrationTable()
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    mut_positions = {}
    if deep_history:
        # initial population
        init_ts = msprime.simulate(N, recombination_rate=1.0)
        init_ts.dump_tables(nodes=nodes, edgesets=edgesets)
        nodes.set_columns(time=nodes.time + ngens, flags=nodes.flags)
    else:
        for _ in range(N):
            nodes.add_row(time=ngens)

    pop = list(range(N))
    for t in range(ngens - 1, -1, -1):
        if debug:
            print("t:", t)
            print("pop:", pop)

        dead = [random.random() > survival for k in pop]
        # sample these first so that all parents are from the previous gen
        new_parents = [
            (random.choice(pop), random.choice(pop)) for k in range(sum(dead))]
        k = 0
        if debug:
            print("Replacing", sum(dead), "individuals.")
        for j in range(N):
            if dead[j]:
                # this is: offspring ID, lparent, rparent, breakpoint
                offspring = nodes.num_rows
                nodes.add_row(time=t)
                lparent, rparent = new_parents[k]
                k += 1
                bp = random_breakpoint()
                muts = random_mutations(mutation_rate)
                if debug:
                    print("--->", offspring, lparent, rparent, bp)
                pop[j] = offspring
                if bp > 0.0:
                    edgesets.add_row(
                        left=0.0, right=bp, parent=lparent, children=(offspring,))
                if bp < 1.0:
                    edgesets.add_row(
                        left=bp, right=1.0, parent=rparent, children=(offspring,))
                for mut in muts:
                    if mut not in mut_positions:
                        mut_positions[mut] = sites.num_rows
                        sites.add_row(position=mut, ancestral_state=random_allele())
                    mutations.add_row(
                            site=mut_positions[mut], node=offspring,
                            derived_state=random_allele())

    if debug:
        print("Done! Final pop:")
        print(pop)
    flags = [
        (msprime.NODE_IS_SAMPLE if u in pop else 0) for u in range(nodes.num_rows)]
    nodes.set_columns(time=nodes.time, flags=flags)
    if debug:
        print("Done.")
        print("Nodes:")
        print(nodes)
        print("Edgesets:")
        print(edgesets)
        print("Sites:")
        print(sites)
        print("Mutations:")
        print(mutations)
        # print("Migrations:")
        # print(tables.migrations)

    return msprime.TableTuple(nodes, edgesets, migrations, sites, mutations)


def get_tree_roots(ts, tree):
    """
    Returns the set of roots in the specified tree in the specified tree
    sequence. This is defined as the unique end points of paths from
    the samples.

    NOTE: This should be removed when the appropriate method has been
    implemented in the Tree class.
    """
    roots = set()
    for u in ts.samples():
        while tree.parent(u) != msprime.NULL_NODE:
            u = tree.parent(u)
        roots.add(u)
    return roots


class TestSimulation(unittest.TestCase):
    """
    Tests that the simulations produce the output we expect.
    """
    random_seed = 5678

    def test_non_overlapping_generations(self):
        tables = wf_sim(N=10, ngens=10, survival=0.0, seed=self.random_seed)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edgesets.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edgesets = tables.edgesets
        msprime.sort_tables(nodes=nodes, edgesets=edgesets)
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.simplify_tables(samples=samples, nodes=nodes, edgesets=edgesets)
        # We should have a complete tree sequence with no unary nodes.
        ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
        for e in ts.edgesets():
            self.assertGreater(len(e.children), 1)
        # All trees should have exactly one root and the leaves should be the samples.
        for tree in ts.trees():
            leaves = set(tree.leaves(tree.root))
            self.assertEqual(leaves, set(ts.samples()))
            roots = get_tree_roots(ts, tree)
            self.assertEqual(roots, set([tree.root]))

    def test_overlapping_generations(self):
        tables = wf_sim(N=30, ngens=10, survival=0.85, seed=self.random_seed)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edgesets.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edgesets = tables.edgesets
        msprime.sort_tables(nodes=nodes, edgesets=edgesets)
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.simplify_tables(samples=samples, nodes=nodes, edgesets=edgesets)
        ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
        for tree in ts.trees():
            roots = get_tree_roots(ts, tree)
            self.assertEqual(roots, {tree.root})

    @unittest.skip("Skipping simplify bug.")
    def test_one_generation_no_deep_history(self):
        N = 20
        tables = wf_sim(N=N, ngens=1, deep_history=False, seed=self.random_seed)
        self.assertEqual(tables.nodes.num_rows, 2 * N)
        self.assertGreater(tables.edgesets.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edgesets = tables.edgesets
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.sort_tables(nodes=nodes, edgesets=edgesets)
        msprime.simplify_tables(samples=samples, nodes=nodes, edgesets=edgesets)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edgesets.num_rows, 0)
        ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
        for tree in ts.trees():
            roots = get_tree_roots(ts, tree)
            all_samples = set()
            for root in roots:
                root_samples = set(tree.samples(root))
                self.assertEqual(len(root_samples & all_samples), 0)
                all_samples |= root_samples
            self.assertEqual(all_samples, set(ts.samples()))

    def test_many_generations_no_deep_history(self):
        N = 10
        ngens = 100
        tables = wf_sim(N=N, ngens=ngens, deep_history=False, seed=self.random_seed)
        self.assertEqual(tables.nodes.num_rows, N * (ngens + 1))
        self.assertGreater(tables.edgesets.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edgesets = tables.edgesets
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.sort_tables(nodes=nodes, edgesets=edgesets)
        msprime.simplify_tables(samples=samples, nodes=nodes, edgesets=edgesets)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edgesets.num_rows, 0)
        # We are assuming that everything has coalesced and we have single-root trees
        ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
        for tree in ts.trees():
            num_roots = len(get_tree_roots(ts, tree))
            self.assertEqual(num_roots, 1)

    def test_with_mutations(self):
        N = 10
        ngens = 100
        tables = wf_sim(N=N, ngens=ngens, mutation_rate=1.0,
                        deep_history=False, seed=self.random_seed)
        self.assertGreater(tables.sites.num_rows, 0)
        self.assertGreater(tables.mutations.num_rows, 0)
        nodes = tables.nodes
        edgesets = tables.edgesets
        sites = tables.sites
        mutations = tables.mutations
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.sort_tables(nodes=nodes, edgesets=edgesets,
                            sites=sites, mutations=mutations)
        msprime.simplify_tables(samples=samples, nodes=nodes, edgesets=edgesets,
                                sites=sites, mutations=mutations)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edgesets.num_rows, 0)
        self.assertGreater(tables.sites.num_rows, 0)
        self.assertGreater(tables.mutations.num_rows, 0)
        # We are assuming that everything has coalesced and we have single-root trees
        ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
        self.assertEqual(ts.sample_size, N)
        print(sites)
        print(mutations)
        for hap in ts.haplotypes():
            self.assertGreaterEqual(len(hap), 0)
            self.assertLessEqual(len(hap), tables.sites.num_rows)


class TestSimplify(unittest.TestCase):

    def get_wf_sims(self, seed):
        """
        Returns an iterator of example tree sequences produced by
        the WF simulator.
        """
        for N in [5, 10, 20]:
            for surv in [0.0, 0.5, 0.9]:
                for mut in [0.0, 0.5]:
                    tables = wf_sim(N=N, ngens=N, survival=surv,
                                    mutation_rate=mut, seed=seed)
                    msprime.sort_tables(
                        nodes=tables.nodes, edgesets=tables.edgesets,
                        sites=tables.sites, mutations=tables.mutations)
                    self.verify_simulation(tables, ngens=N)
                    yield tables

    def verify_simulation(self, tables, ngens):
        """
        Verify that in the full set of returned tables there is parentage
        information for every individual, except those initially present.
        """
        ts = msprime.load_tables(
            nodes=tables.nodes, edgesets=tables.edgesets,
            sites=tables.sites, mutations=tables.mutations)
        for u in range(tables.nodes.num_rows):
            if tables.nodes.time[u] <= ngens:
                lefts = []
                rights = []
                k = 0
                for edgeset in ts.edgesets():
                    if u in edgeset.children:
                        lefts.append(edgeset.left)
                        rights.append(edgeset.right)
                    k += 1
                lefts.sort()
                rights.sort()
                self.assertEqual(lefts[0], 0.0)
                self.assertEqual(rights[-1], 1.0)
                for k in range(len(lefts) - 1):
                    self.assertEqual(lefts[k + 1], rights[k])

    def verify_simplify(self, ts, new_ts, samples=None):
        """
        Check that trees in `ts` match `new_ts` after mapping `samples` in `ts`
        to `range(len(samples))` in `new_ts`.  Modified from
        `verify_simplify_topology`.
        """
        if samples is None:
            samples = ts.samples()
        sample_map = {k: j for j, k in enumerate(samples)}
        # check trees agree at these points
        locs = [random.random() for _ in range(20)]
        locs += random.sample(list(ts.breakpoints())[:-1], min(20, ts.num_trees))
        locs.sort()
        old_trees = ts.trees()
        new_trees = new_ts.trees()
        old_right = -1
        new_right = -1
        for loc in locs:
            while old_right <= loc:
                old_tree = next(old_trees)
                old_left, old_right = old_tree.get_interval()
            assert old_left <= loc < old_right
            while new_right <= loc:
                new_tree = next(new_trees)
                new_left, new_right = new_tree.get_interval()
            assert new_left <= loc < new_right
            pairs = itertools.islice(itertools.combinations(samples, 2), 500)
            for pair in pairs:
                mapped_pair = [sample_map[u] for u in pair]
                mrca1 = old_tree.get_mrca(*pair)
                self.assertTrue(mrca1 != -1)
                tmrca1 = old_tree.get_time(mrca1)
                mrca2 = new_tree.get_mrca(*mapped_pair)
                self.assertTrue(mrca2 != -1)
                tmrca2 = new_tree.get_time(mrca2)
                self.assertEqual(tmrca1, tmrca2)

    def verify_haplotypes(self, ts, samples=None):
        """
        Check that haplotypes are unchanged by simplify, except for removing
        invariant sites.
        """
        if samples is None:
            samples = ts.samples()
        new_ts = ts.simplify(samples=samples)
        sample_map = {k: j for j, k in enumerate(samples)}
        haps = list(ts.haplotypes())
        new_haps = list(new_ts.haplotypes())
        pos = [s.position for s in ts.sites()]
        new_pos = [s.position for s in new_ts.sites()]
        invariant = [x not in new_pos for x in pos]
        reference_hap = haps[samples[0]]
        print(pos)
        print(new_pos)
        print(invariant)
        print(sample_map)
        print(reference_hap)
        print([(haps[k], new_haps[sample_map[k]]) for k in sample_map])
        for k in sample_map:
            print("k:", k)
            assert k in ts.samples()
            i = 0  # position in new haplotype
            for j, ch in enumerate(haps[k]):
                print("i:", i, "j:", j, invariant[j])
                if invariant[j]:
                    self.assertEqual(reference_hap[j], ch)
                else:
                    self.assertEqual(new_haps[sample_map[k]][i], ch)
                    i += 1

    def test_simplify(self):
        #  check that simplify(big set) -> simplify(subset) equals simplify(subset)
        seed = 23
        random.seed(seed)
        for tables in self.get_wf_sims(seed=seed):
            ts = msprime.load_tables(
                nodes=tables.nodes, edgesets=tables.edgesets,
                sites=tables.sites, mutations=tables.mutations)
            for nsamples in [2, 5, 10]:
                big_ts = ts.simplify(samples=ts.samples())
                sub_samples = random.sample(
                        big_ts.samples(), min(nsamples, len(big_ts.samples())))
                small_ts = ts.simplify(samples=[ts.samples()[k] for k in sub_samples])
                self.verify_simplify(big_ts, small_ts, samples=sub_samples)
                self.verify_haplotypes(big_ts, samples=sub_samples)

    def test_simplify_tables(self):
        seed = 23
        for tables in self.get_wf_sims(seed=seed):
            ts = msprime.load_tables(
                nodes=tables.nodes, edgesets=tables.edgesets,
                sites=tables.sites, mutations=tables.mutations)
            for nsamples in [2, 5, 10]:
                nodes = tables.nodes.copy()
                edgesets = tables.edgesets.copy()
                sites = tables.sites.copy()
                mutations = tables.mutations.copy()
                msprime.simplify_tables(
                    samples=ts.samples(), nodes=nodes, edgesets=edgesets,
                    sites=sites, mutations=mutations)
                big_ts = msprime.load_tables(
                    nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)
                sub_samples = random.sample(
                    big_ts.samples(), min(nsamples, len(big_ts.samples())))
                msprime.simplify_tables(
                    samples=sub_samples, nodes=nodes, edgesets=edgesets,
                    sites=sites, mutations=mutations)
                small_ts = msprime.load_tables(
                    nodes=nodes, edgesets=edgesets, sites=sites, mutations=mutations)
                self.verify_simplify(big_ts, small_ts, samples=sub_samples)
