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
import sys
import unittest

import numpy as np

import _msprime
import msprime


class TestRejectedCommonAncestorEventCounts(unittest.TestCase):
    """
    Tests to see if we get the correct number of rejected commone ancestor
    events from the various models.
    """
    def test_hudson(self):
        threshold = 20
        sim = msprime.simulator_factory(sample_size=10, recombination_rate=10)
        sim.random_generator = msprime.RandomGenerator(2)
        sim.run()
        self.assertGreater(sim.num_common_ancestor_events, threshold)
        self.assertGreater(sim.num_recombination_events, threshold)
        self.assertEqual(sim.num_rejected_common_ancestor_events, 0)

        sim2 = msprime.simulator_factory(
            sample_size=10, recombination_rate=10, model="hudson")
        sim2.random_generator = msprime.RandomGenerator(2)
        sim2.run()
        self.assertEqual(
            sim2.num_common_ancestor_events, sim.num_common_ancestor_events)
        self.assertEqual(
            sim2.num_recombination_events, sim.num_recombination_events)
        self.assertEqual(sim2.num_rejected_common_ancestor_events, 0)

    def test_smc_variants(self):
        for model in ["smc", "smc_prime"]:
            threshold = 20
            sim = msprime.simulator_factory(
                sample_size=10, recombination_rate=5, model=model)
            sim.random_generator = msprime.RandomGenerator(3)
            sim.run()
            self.assertGreater(sim.num_common_ancestor_events, threshold)
            self.assertGreater(sim.num_recombination_events, threshold)
            self.assertGreater(sim.num_rejected_common_ancestor_events, 0)


class TestEdges(unittest.TestCase):
    """
    Tests that the edges have the correct properties.
    """
    def test_gaps(self):
        # SMC simulations should never have adjacent edgesets with
        # a non-zero distance between them and the same parent.
        # First we do a simulation with the standard model to make sure
        # we have plausible parameter values.
        sample_size = 10
        recombination_rate = 20
        random_seed = 1

        ts = msprime.simulate(
            sample_size=sample_size, recombination_rate=recombination_rate,
            random_seed=random_seed)
        edgesets = sorted(ts.edgesets(), key=lambda e: (e.parent, e.left))
        num_found = 0
        for j in range(1, len(edgesets)):
            r = edgesets[j - 1]
            s = edgesets[j]
            if r.right != s.left and r.parent == s.parent:
                num_found += 1
        self.assertGreater(num_found, 10)  # Make a reasonable threshold

        # Now do the same for SMC and SMC'.
        for model in ["smc", "smc_prime"]:
            ts = msprime.simulate(
                sample_size=sample_size, recombination_rate=recombination_rate,
                random_seed=random_seed, model=model)
            edgesets = sorted(ts.edgesets(), key=lambda e: (e.parent, e.left))
            num_found = 0
            for j in range(1, len(edgesets)):
                r = edgesets[j - 1]
                s = edgesets[j]
                if r.right != s.left and r.parent == s.parent:
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
            ("smc_prime", msprime.SmcPrimeApproxCoalescent),
            ("dtwf", msprime.DiscreteTimeWrightFisher)
        ]
        for name, model in simulation_models:
            sim = msprime.simulator_factory(sample_size=10, model=name.upper())
            self.assertIsInstance(sim.model, model)
            sim = msprime.simulator_factory(sample_size=10, model=name.title())
            self.assertIsInstance(sim.model, model)

    def test_model_instances(self):
        for bad_type in [1234, {}]:
            self.assertRaises(
                TypeError, msprime.simulator_factory, sample_size=2, model=bad_type)
        models = [
            msprime.StandardCoalescent(),
            msprime.SmcApproxCoalescent(),
            msprime.SmcPrimeApproxCoalescent(),
            msprime.DiscreteTimeWrightFisher(),
            msprime.BetaCoalescent(),
            msprime.DiracCoalescent(),
        ]
        for model in models:
            sim = msprime.simulator_factory(sample_size=10, model=model)
            self.assertEqual(sim.model, model)


class TestParametricModels(unittest.TestCase):
    """
    Tests for the parametric simulation models.
    """
    def test_beta_coalescent_parameters(self):
        N = 1000
        dbl_max = sys.float_info.max
        for alpha in [1.01, 1.5, 1.99]:
            model = msprime.BetaCoalescent(N, alpha)
            self.assertEqual(model.population_size, N)
            self.assertEqual(model.alpha, alpha)
            self.assertEqual(model.truncation_point, dbl_max)
            d = model.get_ll_representation()
            self.assertEqual(d, {
                "name": "beta",
                "population_size": N,
                "alpha": alpha,
                "truncation_point": dbl_max})
        alpha = 1.5
        for truncation_point in [0, 3, 1e6]:
            model = msprime.BetaCoalescent(N, alpha, truncation_point)
            self.assertEqual(model.population_size, N)
            self.assertEqual(model.alpha, alpha)
            self.assertEqual(model.truncation_point, truncation_point)
            d = model.get_ll_representation()
            self.assertEqual(d, {
                "name": "beta",
                "population_size": N,
                "alpha": alpha,
                "truncation_point": truncation_point})

    def test_dirac_coalescent_parameters(self):
        N = 10
        for psi in [0.01, 0.5, 0.99]:
            for c in [1e-6, 1.0, 1e2]:
                model = msprime.DiracCoalescent(N, psi, c)
                self.assertEqual(model.population_size, N)
                self.assertEqual(model.psi, psi)
                self.assertEqual(model.c, c)
                d = model.get_ll_representation()
                self.assertEqual(d, {
                    "name": "dirac", "population_size": N, "psi": psi, "c": c})


class TestMultipleMergerModels(unittest.TestCase):
    """
    Runs tests on the multiple merger coalescent models.
    """
    def test_dirac_coalescent_kingman_regime(self):
        # When c=0, we should a kingman coalescent and no multiple mergers.
        model = msprime.DiracCoalescent(psi=0.5, c=0)
        ts = msprime.simulate(sample_size=10, model=model, random_seed=2)
        for t in ts.trees():
            for u in t.nodes():
                if t.is_internal(u):
                    self.assertEqual(len(t.children(u)), 2)

    def verify_non_binary(self, ts):
        non_binary = False
        for e in ts.edgesets():
            if len(e.children) > 2:
                non_binary = True
                break
        self.assertTrue(non_binary)

    def test_dirac_coalescent_lambda_regime(self):
        # With large c and psi ~ 1, we should be guaranteed some multiple mergers.
        model = msprime.DiracCoalescent(psi=0.999, c=1000)
        ts = msprime.simulate(sample_size=100, model=model, random_seed=4)
        self.verify_non_binary(ts)

    def test_dirac_coalescent_lambda_regime_recombination(self):
        model = msprime.DiracCoalescent(psi=0.999, c=1)
        ts = msprime.simulate(
            sample_size=100, recombination_rate=100, model=model, random_seed=3)
        self.verify_non_binary(ts)

    def test_dirac_coalescent(self):
        model = msprime.DiracCoalescent(100, 0.3, 10)
        ts = msprime.simulate(sample_size=10, model=model)
        # TODO real tests
        self.assertTrue(ts is not None)

    def test_beta_coalescent(self):
        model = msprime.BetaCoalescent(
            population_size=5, alpha=1.5, truncation_point=10)
        ts = msprime.simulate(sample_size=10, model=model)
        # TODO real tests
        self.assertTrue(ts is not None)

    def test_beta_coalescent_integration_fails(self):
        model = msprime.BetaCoalescent(
            population_size=5, alpha=1e-10, truncation_point=10)
        with self.assertRaises(_msprime.LibraryError):
            msprime.simulate(sample_size=10, model=model)

    def test_dtwf(self):
        model = msprime.DiscreteTimeWrightFisher()
        ts = msprime.simulate(sample_size=10, model=model)
        self.assertTrue(ts is not None)
        self.verify_non_binary(ts)


class TestUnsupportedFullArg(unittest.TestCase):
    """
    Full ARG recording isn't supported on anything except standard coalescent.
    """
    def test_smc(self):
        for model in ["smc", "smc_prime"]:
            self.assertRaises(
                _msprime.LibraryError, msprime.simulate, 10, model=model,
                record_full_arg=True)

    def test_dtwf(self):
        for model in [msprime.DiscreteTimeWrightFisher(10)]:
            self.assertRaises(
                _msprime.LibraryError, msprime.simulate, 10, model=model,
                record_full_arg=True)

    def test_multiple_mergers(self):
        for model in [msprime.BetaCoalescent(10), msprime.DiracCoalescent(10)]:
            self.assertRaises(
                _msprime.LibraryError, msprime.simulate, 10, model=model,
                record_full_arg=True)


class TestMixedModels(unittest.TestCase):
    """
    Tests that we can run mixed simulation models.
    """
    def test_wf_hudson_single_locus(self):
        Ne = 100
        t = 10
        ts = msprime.simulate(
            sample_size=10,
            model=msprime.DiscreteTimeWrightFisher(Ne),
            demographic_events=[
                msprime.SimulationModelChange(t, msprime.StandardCoalescent(Ne))],
            random_seed=2)
        tree = ts.first()
        self.assertEqual(tree.num_roots, 1)
        times = ts.tables.nodes.time
        dtwf_times = times[np.logical_and(times > 0, times < t)]
        self.assertGreater(dtwf_times.shape[0], 0)
        self.assertTrue(np.all(dtwf_times == np.floor(dtwf_times)))
        coalescent_times = times[times > t]
        self.assertGreater(coalescent_times.shape[0], 0)
        self.assertTrue(np.all(coalescent_times != np.floor(coalescent_times)))

    def test_wf_hudson_recombinatation(self):
        Ne = 100
        t = 100
        ts = msprime.simulate(
            sample_size=10,
            model=msprime.DiscreteTimeWrightFisher(Ne),
            recombination_rate=0.1,
            demographic_events=[
                msprime.SimulationModelChange(t, msprime.StandardCoalescent(Ne))],
            random_seed=2)
        tree = ts.first()
        self.assertEqual(tree.num_roots, 1)
        times = ts.tables.nodes.time
        dtwf_times = times[np.logical_and(times > 0, times < t)]
        self.assertGreater(dtwf_times.shape[0], 0)
        self.assertTrue(np.all(dtwf_times == np.floor(dtwf_times)))
        coalescent_times = times[times > t]
        self.assertGreater(coalescent_times.shape[0], 0)
        self.assertTrue(np.all(coalescent_times != np.floor(coalescent_times)))

    def test_wf_hudson_back_and_forth(self):
        Ne = 100
        t1 = 100
        t2 = 200
        ts = msprime.simulate(
            sample_size=10,
            model=msprime.DiscreteTimeWrightFisher(Ne),
            recombination_rate=0.1,
            demographic_events=[
                msprime.SimulationModelChange(t1, msprime.StandardCoalescent(Ne)),
                msprime.SimulationModelChange(t2, msprime.DiscreteTimeWrightFisher(Ne))],
            random_seed=2)
        tree = ts.first()
        self.assertEqual(tree.num_roots, 1)
        times = ts.tables.nodes.time
        dtwf_times = times[np.logical_and(times > 0, times < t1, times > t2)]
        self.assertGreater(dtwf_times.shape[0], 0)
        self.assertTrue(np.all(dtwf_times == np.floor(dtwf_times)))
        coalescent_times = times[np.logical_and(times > t1, times < t2)]
        self.assertGreater(coalescent_times.shape[0], 0)
        self.assertTrue(np.all(coalescent_times != np.floor(coalescent_times)))

    def test_many_models(self):
        Ne = 10000
        ts = msprime.simulate(
            Ne=Ne,
            sample_size=10,
            recombination_rate=0.1,
            demographic_events=[
                msprime.SimulationModelChange(10, msprime.StandardCoalescent(Ne)),
                msprime.SimulationModelChange(20, msprime.SmcApproxCoalescent(Ne)),
                msprime.SimulationModelChange(30, msprime.SmcPrimeApproxCoalescent(Ne)),
                msprime.SimulationModelChange(
                    40, msprime.DiscreteTimeWrightFisher(100)),
                msprime.SimulationModelChange(
                    50, msprime.BetaCoalescent(population_size=10)),
                msprime.SimulationModelChange(60, msprime.StandardCoalescent(0.1))],
            random_seed=10)
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 1)


class TestSingleSweep(unittest.TestCase):
    """
    Tests for the single sweep model.
    """
    def test_incorrect_num_labels(self):
        trajectory = msprime.SweepTrajectory([0, 0.1], [0.1, 0.11])
        model = msprime.SingleSweep(position=0.5, trajectory=trajectory)
        for num_labels in [1, 3, 10]:
            with self.assertRaises(_msprime.LibraryError):
                msprime.simulate(
                    10, recombination_rate=1, model=model, num_labels=num_labels)

    @unittest.skip("Assertion `overlaps != NULL' fails.")
    def test_simple_no_recomb(self):
        trajectory = msprime.SweepTrajectory([0, 0.1], [0.1, 0.11])
        model = msprime.SingleSweep(position=0, trajectory=trajectory)
        msprime.simulate(10, model=model, num_labels=2)

    @unittest.skip("Assertion `overlaps != NULL' fails.")
    def test_simple_recomb(self):
        trajectory = msprime.SweepTrajectory([0, 0.1], [0.1, 0.11])
        model = msprime.SingleSweep(position=0.5, trajectory=trajectory)
        msprime.simulate(10, recombination_rate=1, model=model, num_labels=2)

    @unittest.skip("Assertion `overlaps != NULL' fails.")
    def test_trajectory_no_recomb(self):
        trajectory = msprime.simulate_trajectory(
            start_frequency=0.1, end_frequency=0.8, alpha=0.1, time_slice=0.01,
            random_seed=32)
        model = msprime.SingleSweep(position=0.5, trajectory=trajectory)
        msprime.simulate(10, model=model, num_labels=2)

    def test_simple_example(self):
        trajectory = msprime.SweepTrajectory([0, 0.1, 0.2, 0.3], [0.95, 0.85, 0.75, 0.5])
        recomb_map = msprime.RecombinationMap.uniform_map(
            length=100, num_loci=100, rate=0.1)
        model = msprime.SingleSweep(
            population_size=100, position=5, trajectory=trajectory)
        ts = msprime.simulate(
            sample_size=10,
            recombination_map=recomb_map, model=model, num_labels=2,
            demographic_events=[
                msprime.SimulationModelChange(0.3, msprime.StandardCoalescent(100))])
        for tree in ts.trees():
            self.assertEqual(tree.num_roots, 1)


class TestSimulateTrajectory(unittest.TestCase):
    """
    Tests the trajectory simulator.
    """
    def test_same_random_seed_equal(self):
        trajectories = [msprime.simulate_trajectory(
            start_frequency=0.125, end_frequency=0.75, alpha=0.125, time_slice=1 / 16,
            random_seed=42) for _ in range(10)]
        for trajectory in trajectories:
            self.assertTrue(np.array_equal(trajectory.time, trajectories[0].time))
            self.assertTrue(
                np.array_equal(
                    trajectory.allele_frequency, trajectories[0].allele_frequency))

    # TODO more tests as we define functionality here.
