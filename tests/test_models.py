#
# Copyright (C) 2016 University of Oxford
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
Test cases for simulation models to see if they have the correct
basic properties.
"""
from __future__ import print_function
from __future__ import division

import sys
import unittest

import _msprime
import msprime


class TestRejectedCommonAncestorEventCounts(unittest.TestCase):
    """
    Tests to see if we get the correct number of rejected commone ancestor
    events from the various models.
    """
    def test_hudson(self):
        threshold = 20
        sim = msprime.simulator_factory(sample_size=10, recombination_rate=5)
        sim.run()
        self.assertGreater(sim.get_num_common_ancestor_events(), threshold)
        self.assertGreater(sim.get_num_recombination_events(), threshold)
        self.assertEqual(sim.get_num_rejected_common_ancestor_events(), 0)

        sim = msprime.simulator_factory(
            sample_size=10, recombination_rate=5, model="hudson")
        sim.run()
        self.assertGreater(sim.get_num_common_ancestor_events(), threshold)
        self.assertGreater(sim.get_num_recombination_events(), threshold)
        self.assertEqual(sim.get_num_rejected_common_ancestor_events(), 0)

    def test_smc_variants(self):
        for model in ["smc", "smc_prime"]:
            threshold = 20
            sim = msprime.simulator_factory(
                sample_size=10, recombination_rate=5, model=model)
            sim.run()
            self.assertGreater(sim.get_num_common_ancestor_events(), threshold)
            self.assertGreater(sim.get_num_recombination_events(), threshold)
            self.assertGreater(sim.get_num_rejected_common_ancestor_events(), 0)


class TestCoalescenceRecords(unittest.TestCase):
    """
    Tests that the coalescence records have the correct properties.
    """
    def test_gaps(self):
        # SMC simulations should never have adjacent coalescence records with
        # a non-zero distance between them and the same time/node value.
        # First we do a simulation with the standard model to make sure
        # we have plausible parameter values.
        sample_size = 10
        recombination_rate = 20
        random_seed = 1

        ts = msprime.simulate(
            sample_size=sample_size, recombination_rate=recombination_rate,
            random_seed=random_seed)
        records = list(ts.records())
        num_found = 0
        for j in range(1, len(records)):
            r = records[j - 1]
            s = records[j]
            if r.right != s.left and r.node == s.node:
                num_found += 1
        self.assertGreater(num_found, 10)  # Make a reasonable threshold

        # Now do the same for SMC and SMC'.
        for model in ["smc", "smc_prime"]:
            ts = msprime.simulate(
                sample_size=sample_size, recombination_rate=recombination_rate,
                random_seed=random_seed, model=model)
            records = list(ts.records())
            num_found = 0
            for j in range(1, len(records)):
                r = records[j - 1]
                s = records[j]
                if r.right != s.left and r.node == s.node:
                    num_found += 1
            self.assertEqual(num_found, 0)


class TestModelParsing(unittest.TestCase):
    """
    Tests the parsing code for model strings.
    """
    def test_bad_models(self):
        for bad_model in ["NOT", "",  "MODEL"]:
            self.assertRaises(ValueError, msprime.simulate, 10, model=bad_model)

    def test_named_model_variants(self):
        simulation_models = [
            ("hudson", msprime.StandardCoalescent),
            ("smc", msprime.SmcApproxCoalescent),
            ("smc_prime", msprime.SmcPrimeApproxCoalescent)
        ]
        for name, model in simulation_models:
            sim = msprime.simulator_factory(sample_size=10, model=name.upper())
            self.assertIsInstance(sim.get_model(), model)
            sim = msprime.simulator_factory(sample_size=10, model=name.title())
            self.assertIsInstance(sim.get_model(), model)

    def test_model_instances(self):
        for bad_type in [1234, {}]:
            self.assertRaises(
                TypeError, msprime.simulator_factory, sample_size=2, model=bad_type)
        models = [
            msprime.StandardCoalescent(),
            msprime.SmcApproxCoalescent(),
            msprime.SmcPrimeApproxCoalescent(),
            msprime.BetaCoalescent(),
            msprime.DiracCoalescent(),
        ]
        for model in models:
            sim = msprime.simulator_factory(sample_size=10, model=model)
            self.assertEqual(sim.get_model(), model)


class TestParametricModels(unittest.TestCase):
    """
    Tests for the parametric simulation models.
    """
    def test_beta_coalescent_parameters(self):
        dbl_max = sys.float_info.max
        for alpha in [-1, 0, 1.1]:
            model = msprime.BetaCoalescent(alpha)
            self.assertEqual(model.alpha, alpha)
            self.assertEqual(model.truncation_point, dbl_max)
            d = model.get_ll_representation()
            self.assertEqual(d, {
                "name": "beta",
                "alpha": alpha,
                "truncation_point": dbl_max})
        alpha = 2
        for truncation_point in [0, 3, 1e6]:
            model = msprime.BetaCoalescent(alpha, truncation_point)
            self.assertEqual(model.alpha, alpha)
            self.assertEqual(model.truncation_point, truncation_point)
            d = model.get_ll_representation()
            self.assertEqual(d, {
                "name": "beta",
                "alpha": alpha,
                "truncation_point": truncation_point})

    def test_dirac_coalescent_parameters(self):
        for psi in [0.01, 0.5, 0.99]:
            for c in [1e-6, 1.0, 1e2]:
                model = msprime.DiracCoalescent(psi, c)
                self.assertEqual(model.psi, psi)
                self.assertEqual(model.c, c)
                d = model.get_ll_representation()
                self.assertEqual(d, {"name": "dirac", "psi": psi, "c": c})


class TestMultipleMergerModels(unittest.TestCase):
    """
    Runs tests on the multiple merger coalescent models.
    """
    def test_dirac_coalescent(self):
        model = msprime.DiracCoalescent(0.3, 10)
        ts = msprime.simulate(sample_size=10, model=model)
        # TODO real tests
        self.assertTrue(ts is not None)

    def test_beta_coalescent(self):
        model = msprime.BetaCoalescent(5, 10)
        ts = msprime.simulate(sample_size=10, model=model)
        # TODO real tests
        self.assertTrue(ts is not None)


class TestUnsupportedDemographicEvents(unittest.TestCase):
    """
    Some demographic events are not supported until specific models.
    """
    def test_smc_bottlenecks(self):
        # TODO we should have a better exception here.
        for model in ["smc", "smc_prime"]:
            self.assertRaises(
                _msprime.InputError, msprime.simulate, 10, model=model,
                demographic_events=[msprime.SimpleBottleneck(1)])
            self.assertRaises(
                _msprime.InputError, msprime.simulate, 10, model=model,
                demographic_events=[msprime.InstantaneousBottleneck(1)])
