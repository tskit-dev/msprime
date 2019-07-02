#
# Copyright (C) 2019 University of Oxford
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
import math
import msprime
import numpy as np
import tskit
import unittest


class TestLikelihood(unittest.TestCase):
    """
    Tests for likelihood evaluation with the full ARG.
    """
    def test_log_likelihoods(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)

        tables.edges.add_row(left=0, right=0.5, parent=3, child=2)
        tables.edges.add_row(left=0.5, right=1, parent=4, child=2)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)

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
            self.assertTrue(math.isclose(log_arg_likelihood_exact,
                                         msprime.log_arg_likelihood(arg, r)))

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 19 / 8
        for t in theta:
            unnormalised_mutation_ll_exact = (5 * math.log(tree_length * t) -
                                              tree_length * t)
            unnormalised_mutation_ll_exact -= 2 * math.log(4 * tree_length)
            unnormalised_mutation_ll_exact -= 2 * math.log(tree_length)
            unnormalised_mutation_ll_exact += math.log(3 / (4 * tree_length))
            self.assertTrue(math.isclose(
                            unnormalised_mutation_ll_exact,
                            msprime.unnormalised_log_mutation_likelihood(arg, t)))

    def test_multiple_mrcas(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)

        tables.edges.add_row(left=0, right=0.5, parent=2, child=1)
        tables.edges.add_row(left=0.5, right=1, parent=3, child=1)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)

        tables.edges.add_row(left=0, right=0.5, parent=4, child=0)
        tables.edges.add_row(left=0.5, right=1, parent=5, child=0)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.15)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.15)

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
            self.assertTrue(math.isclose(log_arg_likelihood_exact,
                                         msprime.log_arg_likelihood(arg, r)))

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 1.5
        for t in theta:
            unnormalised_mutation_ll_exact = (3 * math.log(tree_length * t) -
                                              tree_length * t)
            unnormalised_mutation_ll_exact -= math.log(tree_length)
            unnormalised_mutation_ll_exact -= 2 * math.log(2 * tree_length)
            self.assertTrue(math.isclose(
                            unnormalised_mutation_ll_exact,
                            msprime.unnormalised_log_mutation_likelihood(arg, t)))

    def test_merger_with_overhang(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)

        tables.edges.add_row(left=0, right=0.5, parent=2, child=1)
        tables.edges.add_row(left=0.5, right=1, parent=3, child=1)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)

        tables.edges.add_row(left=0, right=0.7, parent=4, child=0)
        tables.edges.add_row(left=0.7, right=1, parent=5, child=0)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.15)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.15)

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
            self.assertTrue(math.isclose(log_arg_likelihood_exact,
                                         msprime.log_arg_likelihood(arg, r)))

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 81 / 50
        for t in theta:
            unnormalised_mutation_ll_exact = (3 * math.log(tree_length * t) -
                                              tree_length * t)
            unnormalised_mutation_ll_exact -= math.log(2 * tree_length)
            unnormalised_mutation_ll_exact += math.log(1.3 / tree_length)
            unnormalised_mutation_ll_exact -= math.log(tree_length)
            self.assertTrue(math.isclose(
                            unnormalised_mutation_ll_exact,
                            msprime.unnormalised_log_mutation_likelihood(arg, t)))

    def test_gap_in_ancestral_material(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)

        tables.edges.add_row(left=0, right=0.3, parent=2, child=0)
        tables.edges.add_row(left=0.3, right=1, parent=3, child=0)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)

        tables.edges.add_row(left=0.3, right=0.5, parent=4, child=3)
        tables.edges.add_row(left=0.5, right=1, parent=5, child=3)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.2)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.2)

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
            self.assertTrue(math.isclose(log_arg_likelihood_exact,
                                         msprime.log_arg_likelihood(arg, r)))

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 0.84
        for t in theta:
            unnormalised_mutation_ll_exact = (5 * math.log(tree_length * t) -
                                              tree_length * t)
            unnormalised_mutation_ll_exact += 3 * math.log(0.4 / tree_length)
            unnormalised_mutation_ll_exact += 2 * math.log(0.5 / tree_length)
            self.assertTrue(math.isclose(
                            unnormalised_mutation_ll_exact,
                            msprime.unnormalised_log_mutation_likelihood(arg, t)))

    def test_recombination_in_material_gap(self):
        tables = tskit.TableCollection(sequence_length=1)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, population=0,
                             individual=-1, time=0)

        tables.edges.add_row(left=0, right=0.3, parent=2, child=0)
        tables.edges.add_row(left=0.3, right=1, parent=3, child=0)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.1)

        tables.edges.add_row(left=0.3, right=0.5, parent=4, child=3)
        tables.edges.add_row(left=0.5, right=1, parent=5, child=3)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.2)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.2)

        tables.edges.add_row(left=0, right=0.3, parent=6, child=2)
        tables.edges.add_row(left=0.5, right=1, parent=6, child=5)
        tables.nodes.add_row(flags=0, population=0, individual=-1, time=0.3)

        tables.edges.add_row(left=0, right=0.3, parent=7, child=6)
        tables.edges.add_row(left=0.5, right=1, parent=8, child=6)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.4)
        tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, population=0,
                             individual=-1, time=0.4)

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
            self.assertTrue(math.isclose(log_arg_likelihood_exact,
                                         msprime.log_arg_likelihood(arg, r)))

        theta = np.arange(0.1, 10, 0.1)
        tree_length = 1.34
        for t in theta:
            unnormalised_mutation_ll_exact = (4 * math.log(tree_length * t) -
                                              tree_length * t)
            unnormalised_mutation_ll_exact += math.log(0.6 / tree_length)
            unnormalised_mutation_ll_exact += 3 * math.log(0.7 / tree_length)
            self.assertTrue(math.isclose(
                            unnormalised_mutation_ll_exact,
                            msprime.unnormalised_log_mutation_likelihood(arg, t)))
