#
# Copyright (C) 2015-2017 University of Oxford
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

try:
    # We use the zip as iterator functionality here.
    from future_builtins import zip
except ImportError:
    # This fails for Python 3.x, but that's fine.
    pass

import itertools
import math
import os
import random
import sys
import six
import unittest

import numpy as np

import msprime
import tests


def get_wf_sims(seed):
    """
    Returns an iterator of example tree sequences produced by 
    the WF simulator.
    """
    for N in [5, 10]:
        for surv in [0.0, 0.5]:
            tables = tests.wf_sim(N=N, ngens=N, survival=surv, seed=seed)
            msprime.sort_tables(nodes=tables.nodes, edgesets=tables.edgesets,
                                sites=tables.sites, mutations=tables.mutations)
            verify_simulation(tables, ngens=N)
            yield tables


def verify_simulation(tables, ngens):
    """
    Verify that in the full set of returned tables there is parentage
    information for every individual, except those initially present.
    """
    ts = msprime.load_tables(nodes=tables.nodes, edgesets=tables.edgesets,
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
            assert lefts[0] == 0.0
            assert rights[-1] == 1.0
            for k in range(len(lefts)-1):
                assert lefts[k+1] == rights[k]

class TestWFsim(unittest.TestCase):
    """
    Tests that tables produced by the WF sim are properly dealt with.
    """

    def verify_simplify(self, ts, new_ts, samples=None):
        '''
        Check that trees in `ts` match `new_ts` after mapping `samples` in `ts
        to `range(len(samples))` in `new_ts`.  Modified from
        `verify_simplify_topology`.
        '''
        if samples is None:
            samples = ts.samples()
        #
        for x in ts.dump_tables():
            print(x)
        for x in new_ts.dump_tables():
            print(x)
        #
        sample_map = {k: j for j, k in enumerate(samples)}
        # check trees agree at these points
        locs = [random.random() for _ in range(20)]
        locs += random.sample(list(ts.breakpoints())[:-1], 
                              min(20, ts.num_trees))
        locs.sort()
        old_trees = ts.trees()
        new_trees = new_ts.trees()
        old_right = -1
        new_right = -1
        for loc in locs:
            while old_right <= loc:
                old_tree = next(old_trees)
                old_left, old_right = old_tree.get_interval()
            while new_right <= loc:
                new_tree = next(new_trees)
                new_left, new_right = new_tree.get_interval()
            pairs = itertools.islice(itertools.combinations(samples, 2), 500)
            print("Loc:", loc)
            print(old_tree.parent_dict)
            print(new_tree.parent_dict)
            for pair in pairs:
                mapped_pair = [sample_map[u] for u in pair]
                print(pair, mapped_pair)
                mrca1 = old_tree.get_mrca(*pair)
                self.assertTrue(mrca1 != -1)
                tmrca1 = old_tree.get_time(mrca1)
                mrca2 = new_tree.get_mrca(*mapped_pair)
                self.assertTrue(mrca2 != -1)
                tmrca2 = new_tree.get_time(mrca2)
                self.assertEqual(tmrca1, tmrca2)

    def test_simplify(self):
        """
        check that simplify(big set) -> simplify(subset)
        equals simplify(subset)
        """
        seed = 23
        random.seed(seed)
        for tables in get_wf_sims(seed=seed):
            for x in tables:
                print(x)
            ts = msprime.load_tables(nodes=tables.nodes, edgesets=tables.edgesets,
                                     sites=tables.sites, mutations=tables.mutations)
            for nsamples in [2, 5, 10]:
                big_ts = ts.simplify(samples=ts.samples())
                sub_samples = random.sample(big_ts.samples(), min(nsamples, len(big_ts.samples())))
                print("samples:", nsamples, sub_samples)
                small_ts = ts.simplify(samples=[ts.samples()[k] for k in sub_samples])
                self.verify_simplify(big_ts, small_ts, samples=sub_samples)

    @unittest.skip("error at second MSP_BAD_PARAM_VALUE of MSP_ERR_BAD_PARAM_VALUE")
    def dont_test_simplify_tables(self):
        seed = 23
        for tables in get_wf_sims(seed=seed):
            ts = msprime.load_tables(nodes=tables.nodes, edgesets=tables.edgesets,
                                     sites=tables.sites, mutations=tables.mutations)
            for nsamples in [2, 5, 10]:
                nodes = tables.nodes.copy()
                edgesets = tables.edgesets.copy()
                sites = tables.sites.copy()
                mutations = tables.mutations.copy()
                msprime.simplify_tables(samples=ts.samples(),
                                        nodes=nodes, edgesets=edgesets,
                                        sites=sites, mutations=mutations)
                big_ts = msprime.load_tables(nodes=nodes, edgesets=edgesets,
                                             sites=sites, mutations=mutations)
                sub_samples = random.sample(big_ts.samples(), min(nsamples, len(big_ts.samples())))
                msprime.simplify_tables(samples=sub_samples,
                                        nodes=nodes, edgesets=edgesets,
                                        sites=sites, mutations=mutations)
                small_ts = msprime.load_tables(nodes=nodes, edgesets=edgesets,
                                             sites=sites, mutations=mutations)
                self.verify_simplify(big_ts, small_ts, samples=sub_samples)
