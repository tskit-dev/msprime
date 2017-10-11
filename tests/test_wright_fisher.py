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

import msprime
import tests


def random_breakpoint():
    return min(1.0, max(0.0, 2 * random.random() - 0.5))


def random_mutations(rate):
    nmuts = np.random.poisson(lam=rate)
    return [random.random() for _ in range(nmuts)]


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
    edges = msprime.EdgeTable()
    migrations = msprime.MigrationTable()
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    if deep_history:
        # initial population
        init_ts = msprime.simulate(N, recombination_rate=1.0)
        init_ts.dump_tables(nodes=nodes, edges=edges)
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
                # muts = random_mutations(mutation_rate)
                if debug:
                    print("--->", offspring, lparent, rparent, bp)
                pop[j] = offspring
                if bp > 0.0:
                    edges.add_row(left=0.0, right=bp, parent=lparent, child=offspring)
                if bp < 1.0:
                    edges.add_row(left=bp, right=1.0, parent=rparent, child=offspring)
                # for mut in muts:
                #     # Pretty sure this won't work.
                #     if mut not in mut_positions:
                #         mut_positions[mut] = sites.num_rows
                #         sites.add_row(site=mut, ancestral_state=random_allele())
                #     mutations.add_row(
                #         site=mut_positions[mut], node=offspring,
                # derived_state=random_allele())

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
        print("Edges:")
        print(edges)
        # print("Sites:")
        # print(tables.sites)
        # print("Mutations:")
        # print(tables.mutations)
        # print("Migrations:")
        # print(tables.migrations)

    return msprime.TableTuple(nodes, edges, migrations, sites, mutations)


class TestSimulation(unittest.TestCase):
    """
    Tests that the simulations produce the output we expect.
    """
    random_seed = 5678

    def test_non_overlapping_generations(self):
        tables = wf_sim(N=10, ngens=10, survival=0.0, seed=self.random_seed)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edges.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edges = tables.edges
        msprime.sort_tables(nodes=nodes, edges=edges)
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.simplify_tables(samples=samples, nodes=nodes, edges=edges)
        ts = msprime.load_tables(nodes=nodes, edges=edges)
        # All trees should have exactly one root and the leaves should be the samples,
        # and all internal nodes should have arity > 1
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 1)
            leaves = set(tree.leaves(tree.root))
            self.assertEqual(leaves, set(ts.samples()))
            for u in tree.nodes():
                if tree.is_internal(u):
                    self.assertGreater(len(tree.children(u)), 1)

    def test_overlapping_generations(self):
        tables = wf_sim(N=30, ngens=10, survival=0.85, seed=self.random_seed)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edges.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edges = tables.edges
        msprime.sort_tables(nodes=nodes, edges=edges)
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.simplify_tables(samples=samples, nodes=nodes, edges=edges)
        ts = msprime.load_tables(nodes=nodes, edges=edges)
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 1)

    def test_one_generation_no_deep_history(self):
        N = 20
        tables = wf_sim(N=N, ngens=1, deep_history=False, seed=self.random_seed)
        self.assertEqual(tables.nodes.num_rows, 2 * N)
        self.assertGreater(tables.edges.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edges = tables.edges
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.sort_tables(nodes=nodes, edges=edges)
        msprime.simplify_tables(samples=samples, nodes=nodes, edges=edges)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edges.num_rows, 0)
        ts = msprime.load_tables(nodes=nodes, edges=edges)
        for tree in ts.trees():
            all_samples = set()
            for root in tree.roots:
                root_samples = set(tree.samples(root))
                self.assertEqual(len(root_samples & all_samples), 0)
                all_samples |= root_samples
            self.assertEqual(all_samples, set(ts.samples()))

    def test_many_generations_no_deep_history(self):
        N = 10
        ngens = 100
        tables = wf_sim(N=N, ngens=ngens, deep_history=False, seed=self.random_seed)
        self.assertEqual(tables.nodes.num_rows, N * (ngens + 1))
        self.assertGreater(tables.edges.num_rows, 0)
        self.assertEqual(tables.sites.num_rows, 0)
        self.assertEqual(tables.mutations.num_rows, 0)
        self.assertEqual(tables.migrations.num_rows, 0)
        nodes = tables.nodes
        edges = tables.edges
        samples = np.where(nodes.flags == msprime.NODE_IS_SAMPLE)[0].astype(np.int32)
        msprime.sort_tables(nodes=nodes, edges=edges)
        msprime.simplify_tables(samples=samples, nodes=nodes, edges=edges)
        self.assertGreater(tables.nodes.num_rows, 0)
        self.assertGreater(tables.edges.num_rows, 0)
        # We are assuming that everything has coalesced and we have single-root trees
        ts = msprime.load_tables(nodes=nodes, edges=edges)
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 1)


class TestSimplify(unittest.TestCase):
    """
    Tests for simplify on cases generated by the Wright-Fisher simulator.
    """
    def assertTreeSequencesEqual(self, ts1, ts2):
        self.assertEqual(list(ts1.samples()), list(ts2.samples()))
        self.assertEqual(ts1.sequence_length, ts2.sequence_length)
        ts1_tables = ts1.dump_tables()
        ts2_tables = ts2.dump_tables()
        self.assertEqual(ts1_tables.nodes, ts2_tables.nodes)
        self.assertEqual(ts1_tables.edges, ts2_tables.edges)
        self.assertEqual(ts1_tables.sites, ts2_tables.sites)
        self.assertEqual(ts1_tables.mutations, ts2_tables.mutations)

    def get_wf_sims(self, seed):
        """
        Returns an iterator of example tree sequences produced by the WF simulator.
        """
        for N in [5, 10, 20]:
            for surv in [0.0, 0.5, 0.9]:
                tables = wf_sim(N=N, ngens=N, survival=surv, seed=seed)
                msprime.sort_tables(
                    nodes=tables.nodes, edges=tables.edges,
                    sites=tables.sites, mutations=tables.mutations)
                ts = msprime.load_tables(
                    nodes=tables.nodes, edges=tables.edges,
                    sites=tables.sites, mutations=tables.mutations)
                self.verify_simulation(ts, ngens=N)
                yield ts

    def verify_simulation(self, ts, ngens):
        """
        Verify that in the full set of returned tables there is parentage
        information for every individual, except those initially present.
        """
        tables = ts.dump_tables()
        for u in range(tables.nodes.num_rows):
            if tables.nodes.time[u] <= ngens:
                lefts = []
                rights = []
                k = 0
                for edge in ts.edges():
                    if u == edge.child:
                        lefts.append(edge.left)
                        rights.append(edge.right)
                    k += 1
                lefts.sort()
                rights.sort()
                self.assertEqual(lefts[0], 0.0)
                self.assertEqual(rights[-1], 1.0)
                for k in range(len(lefts) - 1):
                    self.assertEqual(lefts[k + 1], rights[k])

    def verify_simplify(self, ts, new_ts, sample_map):
        """
        Check that trees in `ts` match `new_ts` using the specified sample_map.
        Modified from `verify_simplify_topology`.
        """
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
            # print("comparing trees")
            # print("interval:", old_tree.interval)
            # print(old_tree.draw(format="unicode"))
            # print("interval:", new_tree.interval)
            # print(new_tree.draw(format="unicode"))
            pairs = itertools.islice(itertools.combinations(sample_map.keys(), 2), 500)
            for pair in pairs:
                mapped_pair = [sample_map[u] for u in pair]
                mrca1 = old_tree.get_mrca(*pair)
                self.assertTrue(mrca1 != -1)
                tmrca1 = old_tree.get_time(mrca1)
                mrca2 = new_tree.get_mrca(*mapped_pair)
                self.assertTrue(mrca2 != -1)
                tmrca2 = new_tree.get_time(mrca2)
                self.assertEqual(tmrca1, tmrca2)

    def test_simplify(self):
        #  check that simplify(big set) -> simplify(subset) equals simplify(subset)
        seed = 23
        random.seed(seed)
        for ts in self.get_wf_sims(seed=seed):
            s = tests.Simplifier(ts, ts.samples())
            py_full_ts, py_full_map = s.simplify()
            full_ts, full_map = ts.simplify(ts.samples())
            self.assertEqual(py_full_map, full_map)
            self.assertTreeSequencesEqual(full_ts, py_full_ts)

            for nsamples in [2, 5, 10]:
                sub_samples = random.sample(ts.samples(), min(nsamples, ts.sample_size))
                s = tests.Simplifier(ts, sub_samples)
                py_small_ts, py_small_map = s.simplify()
                small_ts, small_map = ts.simplify(samples=sub_samples)
                self.assertTreeSequencesEqual(small_ts, py_small_ts)
                self.verify_simplify(ts, small_ts, small_map)

    def test_simplify_tables(self):
        seed = 71
        for ts in self.get_wf_sims(seed=seed):
            tables = ts.dump_tables()
            for nsamples in [2, 5, 10]:
                nodes = tables.nodes.copy()
                edges = tables.edges.copy()
                sites = tables.sites.copy()
                mutations = tables.mutations.copy()
                sub_samples = random.sample(ts.samples(), min(nsamples, ts.num_samples))

                new_sample_ids = np.zeros(len(sub_samples), dtype=np.int32)
                msprime.simplify_tables(
                    samples=sub_samples, nodes=nodes, edges=edges,
                    sites=sites, mutations=mutations, sample_map=new_sample_ids)
                sample_map = dict(zip(sub_samples, new_sample_ids))
                small_ts = msprime.load_tables(
                    nodes=nodes, edges=edges, sites=sites, mutations=mutations)
                self.verify_simplify(ts, small_ts, sample_map)
