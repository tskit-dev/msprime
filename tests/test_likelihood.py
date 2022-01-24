#
# Copyright (C) 2019-2020 University of Oxford
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
Test cases for likelihood functions in msprime.
"""
import collections
import math
import unittest

import numpy as np
import pytest
import tskit

import msprime

# Python implementations of the likelihood calculations


def log_arg_likelihood(arg, recombination_rate, Ne=1):
    if recombination_rate < 0:
        raise ValueError("Recombination rate must be >= 0")
    tables = arg.tables
    number_of_lineages = arg.num_samples
    number_of_links = number_of_lineages * arg.sequence_length
    number_of_edges = len(tables.edges)
    edges_above = collections.defaultdict(list)
    for edge in tables.edges:
        edges_above[edge.child].append(edge)
    edge = 0
    time = 0
    ret = 0
    while edge < number_of_edges:
        parent = tables.edges[edge].parent
        rate = (
            number_of_lineages * (number_of_lineages - 1) / (4 * Ne)
            + number_of_links * recombination_rate
        )
        ret -= rate * (tables.nodes[parent].time - time)
        time = tables.nodes[parent].time
        child = tables.edges[edge].child
        if tables.nodes[parent].flags == msprime.NODE_IS_RE_EVENT:
            if recombination_rate == 0:
                ret = -float("inf")
                break
            while tables.edges[edge].child == child:
                if tables.edges[edge].parent != tables.edges[edge - 1].parent:
                    gap = tables.edges[edge].left - tables.edges[edge - 1].right
                edge += 1
            number_of_links -= gap
            number_of_lineages += 1
            if gap == 0:
                # we evaluate the density rather than the probability if there is no gap
                gap = 1
            ret += math.log(recombination_rate * gap)
        else:
            ret -= math.log(2 * Ne)
            segment_length_in_children = -tables.edges.left[edge]
            while edge < number_of_edges and tables.edges[edge].child == child:
                edge += 1
            segment_length_in_children += (
                tables.edges.right[edge - 1] - tables.edges.left[edge]
            )
            child = tables.edges[edge].child
            while edge < number_of_edges and tables.edges[edge].child == child:
                edge += 1
            segment_length_in_children += tables.edges.right[edge - 1]
            if parent in edges_above:
                segment_length_in_parent = (
                    edges_above[parent][-1].right - edges_above[parent][0].left
                )
                number_of_lineages -= 1
                number_of_links -= segment_length_in_children - segment_length_in_parent
            else:
                number_of_lineages -= 2
                number_of_links -= segment_length_in_children
    return ret


class TestKnownExamples:
    """
    Tests for likelihood evaluation with the full ARG in cases where we've
    calculated the exact values beforehand.
    """

    def test_log_likelihoods(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )

        tables.edges.add_row(left=0, right=0.5, parent=3, child=2)
        tables.edges.add_row(left=0.5, right=1, parent=4, child=2)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )

        tables.edges.add_row(left=0, right=1, parent=5, child=1)
        tables.edges.add_row(left=0, right=0.5, parent=5, child=3)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.25)

        tables.edges.add_row(left=0, right=1, parent=6, child=0)
        tables.edges.add_row(left=0.5, right=1, parent=6, child=4)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.5)

        tables.edges.add_row(left=0, right=1, parent=7, child=5)
        tables.edges.add_row(left=0, right=1, parent=7, child=6)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=1)

        tables.mutations.add_row(site=0, node=0, derived_state="1")
        tables.mutations.add_row(site=1, node=1, derived_state="1")
        tables.mutations.add_row(site=2, node=3, derived_state="1")
        tables.mutations.add_row(site=3, node=0, derived_state="1")
        tables.mutations.add_row(site=4, node=5, derived_state="1")

        tables.sites.add_row(0.1, "0")
        tables.sites.add_row(0.2, "0")
        tables.sites.add_row(0.3, "0")
        tables.sites.add_row(0.4, "0")
        tables.sites.add_row(0.45, "0")

        tables.populations.add_row()

        arg = tables.tree_sequence()

        rho = np.arange(0.1, 10, 0.1)
        for r in rho:
            log_arg_likelihood_exact = math.log(r) - (3 + 3 * r) * 0.1
            log_arg_likelihood_exact -= (6 + 3 * r) * 0.15
            log_arg_likelihood_exact -= (3 + 2.5 * r) * 0.25
            log_arg_likelihood_exact -= (1 + 2 * r) * 0.5
            assert math.isclose(
                log_arg_likelihood_exact, msprime.log_arg_likelihood(arg, r, 0.5)
            )

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 19 / 8
        for t in theta:
            mutation_ll_exact = 5 * math.log(tree_length * t) - tree_length * t
            mutation_ll_exact -= 2 * math.log(4 * tree_length)
            mutation_ll_exact -= 2 * math.log(tree_length)
            mutation_ll_exact += math.log(3 / (4 * tree_length))
            assert math.isclose(
                mutation_ll_exact,
                msprime.log_mutation_likelihood(arg, t),
            )

    def test_multiple_mrcas(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )

        tables.edges.add_row(left=0, right=0.5, parent=2, child=1)
        tables.edges.add_row(left=0.5, right=1, parent=3, child=1)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )

        tables.edges.add_row(left=0, right=0.5, parent=4, child=0)
        tables.edges.add_row(left=0.5, right=1, parent=5, child=0)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.15
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.15
        )

        tables.edges.add_row(left=0, right=0.5, parent=6, child=2)
        tables.edges.add_row(left=0, right=0.5, parent=6, child=4)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.5)

        tables.edges.add_row(left=0.5, right=1, parent=7, child=3)
        tables.edges.add_row(left=0.5, right=1, parent=7, child=5)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=1)

        tables.mutations.add_row(site=0, node=1, derived_state="1")
        tables.mutations.add_row(site=1, node=4, derived_state="1")
        tables.mutations.add_row(site=2, node=3, derived_state="1")

        tables.sites.add_row(0.1, "0")
        tables.sites.add_row(0.2, "0")
        tables.sites.add_row(0.7, "0")

        tables.populations.add_row()

        arg = tables.tree_sequence()

        rho = np.arange(0.1, 10, 0.1)
        for r in rho:
            log_arg_likelihood_exact = math.log(r) - (1 + 2 * r) * 0.1
            log_arg_likelihood_exact += math.log(r) - (3 + 2 * r) * 0.05
            log_arg_likelihood_exact -= (6 + 2 * r) * 0.35
            log_arg_likelihood_exact -= (1 + r) * 0.5
            assert math.isclose(
                log_arg_likelihood_exact, msprime.log_arg_likelihood(arg, r, 0.5)
            )

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 1.5
        for t in theta:
            mutation_ll_exact = 3 * math.log(tree_length * t) - tree_length * t
            mutation_ll_exact -= math.log(tree_length)
            mutation_ll_exact -= 2 * math.log(2 * tree_length)
            assert math.isclose(
                mutation_ll_exact,
                msprime.log_mutation_likelihood(arg, t),
            )

    def test_merger_with_overhang(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )

        tables.edges.add_row(left=0, right=0.5, parent=2, child=1)
        tables.edges.add_row(left=0.5, right=1, parent=3, child=1)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )

        tables.edges.add_row(left=0, right=0.7, parent=4, child=0)
        tables.edges.add_row(left=0.7, right=1, parent=5, child=0)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.15
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.15
        )

        tables.edges.add_row(left=0, right=0.5, parent=6, child=2)
        tables.edges.add_row(left=0, right=0.7, parent=6, child=4)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.5)

        tables.edges.add_row(left=0.5, right=1, parent=7, child=3)
        tables.edges.add_row(left=0.7, right=1, parent=7, child=5)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=1)

        tables.edges.add_row(left=0.5, right=0.7, parent=8, child=6)
        tables.edges.add_row(left=0.5, right=0.7, parent=8, child=7)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=1.3)

        tables.mutations.add_row(site=0, node=1, derived_state="1")
        tables.mutations.add_row(site=1, node=4, derived_state="1")
        tables.mutations.add_row(site=2, node=3, derived_state="1")

        tables.sites.add_row(0.1, "0")
        tables.sites.add_row(0.6, "0")
        tables.sites.add_row(0.75, "0")

        tables.populations.add_row()

        arg = tables.tree_sequence()

        rho = np.arange(0.1, 10, 0.1)
        for r in rho:
            log_arg_likelihood_exact = math.log(r) - (1 + 2 * r) * 0.1
            log_arg_likelihood_exact += math.log(r) - (3 + 2 * r) * 0.05
            log_arg_likelihood_exact -= (6 + 2 * r) * 0.35
            log_arg_likelihood_exact -= (3 + r) * 0.5
            log_arg_likelihood_exact -= (1 + 0.4 * r) * 0.3
            assert math.isclose(
                log_arg_likelihood_exact, msprime.log_arg_likelihood(arg, r, 0.5)
            )

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 81 / 50
        for t in theta:
            mutation_ll_exact = 3 * math.log(tree_length * t) - tree_length * t
            mutation_ll_exact -= math.log(2 * tree_length)
            mutation_ll_exact += math.log(1.3 / tree_length)
            mutation_ll_exact -= math.log(tree_length)
            assert math.isclose(
                mutation_ll_exact,
                msprime.log_mutation_likelihood(arg, t),
            )

    def test_gap_in_ancestral_material(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )

        tables.edges.add_row(left=0, right=0.3, parent=2, child=0)
        tables.edges.add_row(left=0.3, right=1, parent=3, child=0)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )

        tables.edges.add_row(left=0.3, right=0.5, parent=4, child=3)
        tables.edges.add_row(left=0.5, right=1, parent=5, child=3)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.2
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.2
        )

        tables.edges.add_row(left=0, right=0.3, parent=6, child=2)
        tables.edges.add_row(left=0.5, right=1, parent=6, child=5)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.3)

        tables.edges.add_row(left=0, right=1, parent=7, child=1)
        tables.edges.add_row(left=0, right=0.3, parent=7, child=6)
        tables.edges.add_row(left=0.5, right=1, parent=7, child=6)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.4)

        tables.edges.add_row(left=0.3, right=0.5, parent=8, child=4)
        tables.edges.add_row(left=0.3, right=0.5, parent=8, child=7)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.5)

        tables.mutations.add_row(site=0, node=2, derived_state="1")
        tables.mutations.add_row(site=1, node=3, derived_state="1")
        tables.mutations.add_row(site=2, node=7, derived_state="1")
        tables.mutations.add_row(site=3, node=6, derived_state="1")
        tables.mutations.add_row(site=4, node=5, derived_state="1")

        tables.sites.add_row(0.1, "0")
        tables.sites.add_row(0.35, "0")
        tables.sites.add_row(0.4, "0")
        tables.sites.add_row(0.6, "0")
        tables.sites.add_row(0.8, "0")

        tables.populations.add_row()

        arg = tables.tree_sequence()

        rho = np.arange(0.1, 10, 0.1)
        for r in rho:
            log_arg_likelihood_exact = math.log(r) - (1 + 2 * r) * 0.1
            log_arg_likelihood_exact += math.log(r) - (3 + 2 * r) * 0.1
            log_arg_likelihood_exact -= (6 + 2 * r) * 0.1
            log_arg_likelihood_exact -= (3 + 2.2 * r) * 0.1
            log_arg_likelihood_exact -= (1 + 0.4 * r) * 0.1
            assert math.isclose(
                log_arg_likelihood_exact, msprime.log_arg_likelihood(arg, r, 0.5)
            )

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 0.84
        for t in theta:
            mutation_ll_exact = 5 * math.log(tree_length * t) - tree_length * t
            mutation_ll_exact += 3 * math.log(0.4 / tree_length)
            mutation_ll_exact += 2 * math.log(0.5 / tree_length)
            assert math.isclose(
                mutation_ll_exact,
                msprime.log_mutation_likelihood(arg, t),
            )

    def test_recombination_in_material_gap(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE, population=0, individual=-1, time=0
        )

        tables.edges.add_row(left=0, right=0.3, parent=2, child=0)
        tables.edges.add_row(left=0.3, right=1, parent=3, child=0)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.1
        )

        tables.edges.add_row(left=0.3, right=0.5, parent=4, child=3)
        tables.edges.add_row(left=0.5, right=1, parent=5, child=3)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.2
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.2
        )

        tables.edges.add_row(left=0, right=0.3, parent=6, child=2)
        tables.edges.add_row(left=0.5, right=1, parent=6, child=5)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.3)

        tables.edges.add_row(left=0, right=0.3, parent=7, child=6)
        tables.edges.add_row(left=0.5, right=1, parent=8, child=6)
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.4
        )
        tables.nodes.add_row(
            flags=msprime.NODE_IS_RE_EVENT, population=0, individual=-1, time=0.4
        )

        tables.edges.add_row(left=0.3, right=0.5, parent=9, child=4)
        tables.edges.add_row(left=0.5, right=1, parent=9, child=8)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.5)

        tables.edges.add_row(left=0, right=1, parent=10, child=1)
        tables.edges.add_row(left=0, right=0.3, parent=10, child=7)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.6)

        tables.edges.add_row(left=0.3, right=1, parent=11, child=9)
        tables.edges.add_row(left=0.3, right=1, parent=11, child=10)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.7)

        tables.mutations.add_row(site=0, node=2, derived_state="1")
        tables.mutations.add_row(site=1, node=3, derived_state="1")
        tables.mutations.add_row(site=2, node=6, derived_state="1")
        tables.mutations.add_row(site=3, node=5, derived_state="1")

        tables.sites.add_row(0.1, "0")
        tables.sites.add_row(0.35, "0")
        tables.sites.add_row(0.6, "0")
        tables.sites.add_row(0.8, "0")

        tables.populations.add_row()

        arg = tables.tree_sequence()

        rho = np.arange(0.1, 10, 0.1)
        for r in rho:
            log_arg_likelihood_exact = math.log(r) - (1 + 2 * r) * 0.1
            log_arg_likelihood_exact += math.log(r) - (3 + 2 * r) * 0.1
            log_arg_likelihood_exact -= (6 + 2 * r) * 0.1
            log_arg_likelihood_exact += math.log(0.2 * r) - (3 + 2.2 * r) * 0.1
            log_arg_likelihood_exact -= (6 + 2 * r) * 0.1
            log_arg_likelihood_exact -= (3 + 2 * r) * 0.1
            log_arg_likelihood_exact -= (1 + 1.4 * r) * 0.1
            assert math.isclose(
                log_arg_likelihood_exact, msprime.log_arg_likelihood(arg, r, 0.5)
            )

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 1.34
        for t in theta:
            mutation_ll_exact = 4 * math.log(tree_length * t) - tree_length * t
            mutation_ll_exact += math.log(0.6 / tree_length)
            mutation_ll_exact += 3 * math.log(0.7 / tree_length)
            assert math.isclose(
                mutation_ll_exact,
                msprime.log_mutation_likelihood(arg, t),
            )


class TestSimulatedExamples(unittest.TestCase):
    """
    Given some simulated ARGs test that the Python likelihood implementation
    computes the same likelihoods as the C code.
    """

    # TODO Add mutation rate as parameter here.
    def verify(self, ts, recombination_rate, Ne):
        l1 = msprime.log_arg_likelihood(
            ts, recombination_rate=recombination_rate, Ne=Ne
        )
        l2 = log_arg_likelihood(ts, recombination_rate, Ne)
        self.assertAlmostEqual(l1, l2)

    def test_no_recombination(self):
        ts = msprime.simulate(10, random_seed=2)
        self.verify(ts, recombination_rate=1, Ne=0.5)
        self.verify(ts, recombination_rate=0.5, Ne=0.5)
        self.verify(ts, recombination_rate=2, Ne=0.5)
        self.verify(ts, recombination_rate=1, Ne=1)
        self.verify(ts, recombination_rate=1, Ne=2)

    def test_small_arg_no_mutation(self):
        ts = msprime.simulate(
            5, recombination_rate=1, random_seed=12, record_full_arg=True
        )
        assert ts.num_edges > 10
        self.verify(ts, recombination_rate=1, Ne=0.5)
        self.verify(ts, recombination_rate=0.5, Ne=0.5)
        self.verify(ts, recombination_rate=2, Ne=0.5)
        self.verify(ts, recombination_rate=1, Ne=1)
        self.verify(ts, recombination_rate=1, Ne=2)

    def test_negative_rec_rate(self):
        ts = msprime.simulate(
            5, recombination_rate=1, random_seed=12, record_full_arg=True
        )
        with pytest.raises(ValueError):
            msprime.log_arg_likelihood(ts, recombination_rate=-1)
        with pytest.raises(ValueError):
            log_arg_likelihood(ts, recombination_rate=-1)

    def test_zero_mut_rate(self):
        # No mutations
        ts = msprime.simulate(
            5, recombination_rate=1, random_seed=12, record_full_arg=True
        )
        lik = msprime.log_mutation_likelihood(ts, 0)
        assert lik == 0

        # With mutations
        ts = msprime.simulate(
            5,
            recombination_rate=1,
            mutation_rate=1,
            random_seed=12,
            record_full_arg=True,
        )
        lik = msprime.log_mutation_likelihood(ts, 0)
        assert lik == float("-inf")


class TestOddTopologies(unittest.TestCase):
    """
    Tests that we give sensible results when we run on weird topologies.
    """

    def verify(self, ts):
        for r in [0.001, 1]:
            l1 = msprime.log_arg_likelihood(ts, r)
            l2 = log_arg_likelihood(ts, r)
            self.assertAlmostEqual(l1, l2)

    def test_zero_edges(self):
        tables = tskit.TableCollection(1)
        for _ in range(2):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE)
        self.verify(tables.tree_sequence())

    def test_no_edges_mutations(self):
        tables = tskit.TableCollection(1)
        for _ in range(2):
            tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE)
        tables.sites.add_row(0, "A")
        tables.mutations.add_row(0, 0, "T")
        self.verify(tables.tree_sequence())
