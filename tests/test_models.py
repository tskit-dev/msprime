#
# Copyright (C) 2016-2020 University of Oxford
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
import numpy as np
import pytest

import msprime
from msprime import _msprime
from msprime import ancestry


class TestIntrospectionInterface:
    """
    Tests that we have meaningful repr and str functions for all the
    classes used in the model hierarchy.
    """

    def test_standard_coalescent(self):
        model = msprime.StandardCoalescent()
        repr_s = "StandardCoalescent()"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_smc_models(self):
        model = msprime.SmcApproxCoalescent()
        repr_s = "SmcApproxCoalescent()"
        assert repr(model) == repr_s
        assert str(model) == repr_s

        model = msprime.SmcPrimeApproxCoalescent()
        repr_s = "SmcPrimeApproxCoalescent()"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_dtwf(self):
        model = msprime.DiscreteTimeWrightFisher()
        repr_s = "DiscreteTimeWrightFisher()"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_wf_pedigree(self):
        model = msprime.WrightFisherPedigree()
        repr_s = "WrightFisherPedigree()"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_beta_coalescent(self):
        model = msprime.BetaCoalescent(alpha=1, truncation_point=2)
        repr_s = "BetaCoalescent(alpha=1, truncation_point=2)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_dirac_coalescent(self):
        model = msprime.DiracCoalescent(psi=123, c=2)
        repr_s = "DiracCoalescent(psi=123, c=2)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_sweep_genic_selection(self):
        model = msprime.SweepGenicSelection(
            position=1, start_frequency=0.5, end_frequency=0.9, s=0.1, dt=0.01
        )
        repr_s = (
            "SweepGenicSelection(position=1, start_frequency=0.5, "
            "end_frequency=0.9, s=0.1, dt=0.01)"
        )
        assert repr(model) == repr_s
        assert str(model) == repr_s


class TestModelFactory:
    """
    Tests that the model_factory function, which we use for instantiating model
    objects in various different ways, works correctly.
    """

    def test_bad_model_names(self):
        for bad_model in ["NOT", "", "MODEL"]:
            with pytest.raises(ValueError):
                ancestry._model_factory(model=bad_model)

    def test_named_model_variants(self):
        simulation_models = [
            ("hudson", msprime.StandardCoalescent),
            ("smc", msprime.SmcApproxCoalescent),
            ("smc_prime", msprime.SmcPrimeApproxCoalescent),
            ("dtwf", msprime.DiscreteTimeWrightFisher),
            ("wf_ped", msprime.WrightFisherPedigree),
        ]
        for name, model_class in simulation_models:
            model = ancestry._model_factory(model=name.upper())
            assert isinstance(model, model_class)
            model = ancestry._model_factory(model=name.title())
            assert isinstance(model, model_class)
            model = ancestry._model_factory(model=name)
            assert isinstance(model, model_class)

    def test_named_parametric_models_fail(self):
        parametric_models = ["beta", "dirac"]
        for name in parametric_models:
            with pytest.raises(ValueError):
                ancestry._model_factory(model=name)

    def test_bad_models(self):
        for bad_type in [1234, {}]:
            with pytest.raises(TypeError):
                ancestry._model_factory(model=bad_type)

    def test_model_instances(self):
        models = [
            msprime.StandardCoalescent(),
            msprime.SmcApproxCoalescent(),
            msprime.SmcPrimeApproxCoalescent(),
            msprime.DiscreteTimeWrightFisher(),
            msprime.WrightFisherPedigree(),
            msprime.SweepGenicSelection(
                position=0.5,
                start_frequency=0.1,
                end_frequency=0.9,
                s=0.01,
                dt=0.01,
            ),
            msprime.BetaCoalescent(alpha=2),
            msprime.DiracCoalescent(psi=1, c=1),
        ]
        for model in models:
            new_model = ancestry._model_factory(model=model)
            assert new_model is model
            assert new_model.__dict__ == model.__dict__


class TestParseModel:
    """
    Tests for the input model parsing code.
    """

    def test_none(self):
        model, events = ancestry._parse_model_arg(None)
        assert model == msprime.StandardCoalescent()
        assert events == []

    def test_single_model(self):
        model, events = ancestry._parse_model_arg("hudson")
        assert model == msprime.StandardCoalescent()
        assert events == []

        model, events = ancestry._parse_model_arg(msprime.StandardCoalescent())
        assert model == msprime.StandardCoalescent()
        assert events == []

        model, events = ancestry._parse_model_arg("dtwf")
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == []

        model, events = ancestry._parse_model_arg(msprime.DiscreteTimeWrightFisher())
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == []

    def test_single_model_list(self):
        model, events = ancestry._parse_model_arg([None])
        assert model == msprime.StandardCoalescent()
        assert events == []

        # Tuples are also accepted as input.
        model, events = ancestry._parse_model_arg((None,))
        assert model == msprime.StandardCoalescent()
        assert events == []

        model, events = ancestry._parse_model_arg(["hudson"])
        assert model == msprime.StandardCoalescent()
        assert events == []

        model, events = ancestry._parse_model_arg([msprime.StandardCoalescent()])
        assert model == msprime.StandardCoalescent()
        assert events == []

        model, events = ancestry._parse_model_arg(["dtwf"])
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == []

        model, events = ancestry._parse_model_arg([msprime.DiscreteTimeWrightFisher()])
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == []

    def test_one_event(self):
        expected_event = msprime.AncestryModelChange(
            time=1.33, model=msprime.StandardCoalescent()
        )
        model, events = ancestry._parse_model_arg(["dtwf", (1.33, "hudson")])
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == [expected_event]

        model, events = ancestry._parse_model_arg(["dtwf", (1.33, None)])
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == [expected_event]

        model, events = ancestry._parse_model_arg(
            ["dtwf", (1.33, msprime.StandardCoalescent())]
        )
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == [expected_event]

        model, events = ancestry._parse_model_arg(["dtwf", expected_event])
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == [expected_event]
        # We should take a copy of the event.
        assert events[0] is not expected_event

        model, events = ancestry._parse_model_arg(["dtwf", (None, None)])
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == [
            msprime.AncestryModelChange(time=None, model=msprime.StandardCoalescent())
        ]

    def test_two_events(self):
        expected_events = [
            msprime.AncestryModelChange(time=1, model=msprime.StandardCoalescent()),
            msprime.AncestryModelChange(time=2, model=msprime.SmcApproxCoalescent()),
        ]
        model, events = ancestry._parse_model_arg(["dtwf", (1, "hudson"), (2, "smc")])
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == expected_events

        model, events = ancestry._parse_model_arg(
            ["dtwf", (1, None), (2, msprime.SmcApproxCoalescent())]
        )
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == expected_events

        model, events = ancestry._parse_model_arg(
            ["dtwf", expected_events[0], (2, msprime.SmcApproxCoalescent())]
        )
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == expected_events

        model, events = ancestry._parse_model_arg(
            ["dtwf", expected_events[0], (2, msprime.SmcApproxCoalescent())]
        )
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == expected_events
        assert events[0] is not expected_events[0]

        model, events = ancestry._parse_model_arg(["dtwf"] + expected_events)
        assert model == msprime.DiscreteTimeWrightFisher()
        assert events == expected_events
        assert events[0] is not expected_events[0]
        assert events[1] is not expected_events[1]

    def test_errors(self):
        with pytest.raises(ValueError):
            ancestry._parse_model_arg([])
        with pytest.raises(ValueError):
            ancestry._parse_model_arg("X")
        # Anything that's not a list or tuple is interpreted as a model
        with pytest.raises(TypeError):
            ancestry._parse_model_arg({})

        for bad_model_change_type in [None, "str", {}]:
            with pytest.raises(TypeError):
                ancestry._parse_model_arg([None, bad_model_change_type])

        for bad_model_change_tuple in [[], [1, None, None]]:
            with pytest.raises(ValueError):
                ancestry._parse_model_arg(["hudson", bad_model_change_tuple])

        for bad_time in ["sdf", [], {}]:
            with pytest.raises(ValueError):
                ancestry._parse_model_arg(["hudson", (bad_time, "hudson")])

        for bad_model_type in [[], {}]:
            with pytest.raises(TypeError):
                ancestry._parse_model_arg(["hudson", (1, bad_model_type)])


class TestRejectedCommonAncestorEventCounts:
    """
    Tests to see if we get the correct number of rejected commone ancestor
    events from the various models.
    """

    def test_hudson(self):
        threshold = 20
        sim = ancestry._parse_simulate(
            sample_size=10,
            recombination_rate=10,
            random_seed=32,
        )
        sim.run()
        assert sim.num_common_ancestor_events > threshold
        assert sim.num_recombination_events > threshold
        assert sim.num_rejected_common_ancestor_events == 0

        sim2 = ancestry._parse_simulate(
            sample_size=10,
            recombination_rate=10,
            model="hudson",
            random_seed=32,
        )
        sim2.run()
        assert sim2.num_common_ancestor_events == sim.num_common_ancestor_events
        assert sim2.num_recombination_events == sim.num_recombination_events
        assert sim2.num_rejected_common_ancestor_events == 0

    def test_smc_variants(self):
        for model in ["smc", "smc_prime"]:
            threshold = 20
            sim = ancestry._parse_simulate(
                sample_size=10, recombination_rate=5, model=model, random_seed=432
            )
            sim.run()
            assert sim.num_rejected_common_ancestor_events > 0
            assert sim.num_common_ancestor_events > threshold
            assert sim.num_recombination_events > threshold


class TestEdges:
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
            sample_size=sample_size,
            recombination_rate=recombination_rate,
            random_seed=random_seed,
        )
        edgesets = sorted(ts.edgesets(), key=lambda e: (e.parent, e.left))
        num_found = 0
        for j in range(1, len(edgesets)):
            r = edgesets[j - 1]
            s = edgesets[j]
            if r.right != s.left and r.parent == s.parent:
                num_found += 1
        assert num_found > 10  # Make a reasonable threshold

        # Now do the same for SMC and SMC'.
        for model in ["smc", "smc_prime"]:
            ts = msprime.simulate(
                sample_size=sample_size,
                recombination_rate=recombination_rate,
                random_seed=random_seed,
                model=model,
            )
            edgesets = sorted(ts.edgesets(), key=lambda e: (e.parent, e.left))
            num_found = 0
            for j in range(1, len(edgesets)):
                r = edgesets[j - 1]
                s = edgesets[j]
                if r.right != s.left and r.parent == s.parent:
                    num_found += 1
            assert num_found == 0


class TestParametricModels:
    """
    Tests for the parametric simulation models.
    """

    def test_beta_coalescent_parameters(self):
        for alpha in [1.01, 1.5, 1.99]:
            model = msprime.BetaCoalescent(alpha, truncation_point=1)
            assert model.alpha == alpha
            assert model.truncation_point == 1
            d = model.get_ll_representation()
            assert d == {"name": "beta", "alpha": alpha, "truncation_point": 1}
        alpha = 1.5
        for truncation_point in [0.01, 0.5, 1]:
            model = msprime.BetaCoalescent(alpha, truncation_point)
            assert model.alpha == alpha
            assert model.truncation_point == truncation_point
            d = model.get_ll_representation()
            assert d == {
                "name": "beta",
                "alpha": alpha,
                "truncation_point": truncation_point,
            }

    def test_dirac_coalescent_parameters(self):
        for psi in [0.01, 0.5, 0.99]:
            for c in [1e-6, 1.0, 1e2]:
                model = msprime.DiracCoalescent(psi, c)
                assert model.psi == psi
                assert model.c == c
                d = model.get_ll_representation()
                assert d == {"name": "dirac", "psi": psi, "c": c}


class TestMultipleMergerModels:
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
                    assert len(t.children(u)) == 2

    def verify_non_binary(self, ts):
        non_binary = False
        for e in ts.edgesets():
            if len(e.children) > 2:
                non_binary = True
                break
        assert non_binary

    def test_dirac_coalescent_lambda_regime(self):
        # With large c and psi ~ 1, we should be guaranteed some multiple mergers.
        model = msprime.DiracCoalescent(psi=0.999, c=1000)
        ts = msprime.simulate(sample_size=100, model=model, random_seed=4)
        self.verify_non_binary(ts)

    def test_dirac_coalescent_lambda_regime_recombination(self):
        model = msprime.DiracCoalescent(psi=0.9, c=100)
        ts = msprime.simulate(
            sample_size=100, recombination_rate=100, model=model, random_seed=3
        )
        self.verify_non_binary(ts)

    def test_dirac_coalescent(self):
        model = msprime.DiracCoalescent(0.3, 10)
        ts = msprime.simulate(Ne=100, sample_size=10, model=model)
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_beta_coalescent(self):
        model = msprime.BetaCoalescent(alpha=1.5)
        ts = msprime.simulate(Ne=5, sample_size=10, model=model)
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_dtwf(self):
        model = msprime.DiscreteTimeWrightFisher()
        ts = msprime.simulate(sample_size=10, model=model)
        assert ts is not None
        self.verify_non_binary(ts)

    def test_wf_ped(self):
        inds = np.array([1, 2, 3, 4])
        parent_indices = np.array([2, 3, 2, 3, -1, -1, -1, -1]).reshape(-1, 2)
        times = np.array([0, 0, 1, 1])
        is_sample = np.array([1, 1, 0, 0])

        model = msprime.WrightFisherPedigree()
        ped = msprime.Pedigree(
            inds, parent_indices, times, is_sample, sex=None, ploidy=2
        )
        ts = msprime.simulate(2, pedigree=ped, model=model)
        assert ts is not None


class TestDtwf:
    """
    Tests for the DTWF model.
    """

    def test_low_recombination(self):
        # https://github.com/tskit-dev/msprime/issues/831
        ts = msprime.simulate(
            10, Ne=1e2, model="dtwf", recombination_rate=1e-9, random_seed=2
        )
        assert ts.num_trees == 1

    def test_very_low_recombination(self):
        ts = msprime.simulate(
            10, Ne=1e2, model="dtwf", recombination_rate=1e-300, random_seed=2
        )
        assert ts.num_trees == 1

    def test_single_recombination(self):
        recombination_map = msprime.RateMap([0, 100, 101, 200], [0, 1, 0])
        ts = msprime.sim_ancestry(
            10,
            population_size=10,
            model="dtwf",
            random_seed=2,
            recombination_rate=recombination_map,
            discrete_genome=True,
        )
        assert ts.num_trees == 2
        trees = ts.trees()
        assert next(trees).interval == (0, 100)
        assert next(trees).interval == (100, 200)

    def test_no_recombination_interval(self):
        positions = [0, 50, 100, 150, 200]
        rates = [0.01, 0.0, 0.1, 0.005, 0.0]
        recombination_map = msprime.RecombinationMap(positions, rates)
        ts = msprime.simulate(
            10, Ne=10, model="dtwf", random_seed=2, recombination_map=recombination_map
        )
        for tree in ts.trees():
            left, right = tree.interval
            assert left < 50 or left > 100
            assert right < 50 or right > 100


class TestUnsupportedFullArg:
    """
    Full ARG recording isn't supported on the discrete time Wright-Fisher model
    """

    def test_dtwf(self):
        for model in [msprime.DiscreteTimeWrightFisher()]:
            with pytest.raises(_msprime.LibraryError):
                msprime.simulate(
                    10,
                    model=model,
                    record_full_arg=True,
                )


class TestMixedModels:
    """
    Tests that we can run mixed simulation models.
    """

    def test_ped_wf_single_locus(self):
        inds = np.array([1, 2, 3, 4, 5, 6])
        parent_indices = np.array([4, 5, 4, 5, 4, 5, 4, 5, -1, -1, -1, -1]).reshape(
            -1, 2
        )
        times = np.array([0, 0, 0, 0, 1, 1])
        is_sample = np.array([1, 1, 1, 1, 0, 0])
        t = max(times)

        model = msprime.WrightFisherPedigree()
        ped = msprime.Pedigree(
            inds, parent_indices, times, is_sample, sex=None, ploidy=2
        )
        ts = msprime.simulate(
            sample_size=4, Ne=2, pedigree=ped, model=(model, (t, "dtwf"))
        )
        tree = ts.first()
        assert tree.num_roots == 1
        all_times = ts.tables.nodes.time
        ped_times = all_times[np.logical_and(all_times > 0, all_times <= t)]
        assert ped_times.shape[0] > 0
        assert np.all(ped_times == np.floor(ped_times))
        wf_times = all_times[all_times > t]
        assert wf_times.shape[0] > 0

    def test_pedigree_unsupported_events(self):
        inds = np.array([1, 2, 3, 4, 5, 6])
        parent_indices = np.array([4, 5, 4, 5, 4, 5, 4, 5, -1, -1, -1, -1]).reshape(
            -1, 2
        )
        times = np.array([0, 0, 0, 0, 1, 1])
        is_sample = np.array([1, 1, 1, 1, 0, 0])
        t = max(times)

        ped = msprime.Pedigree(
            inds, parent_indices, times, is_sample, sex=None, ploidy=2
        )

        bad_model_change = msprime.SimulationModelChange(
            0.5, msprime.DiscreteTimeWrightFisher()
        )
        with pytest.raises(NotImplementedError):
            msprime.simulate(
                4,
                pedigree=ped,
                demographic_events=[bad_model_change],
                model="wf_ped",
            )
        bad_demographic_event = msprime.PopulationParametersChange(t, initial_size=2)
        with pytest.raises(_msprime.LibraryError):
            msprime.simulate(
                4,
                pedigree=ped,
                demographic_events=[bad_demographic_event],
                model="wf_ped",
            )

    def test_ped_wf_recombination(self):
        inds = np.array([1, 2, 3, 4, 5, 6])
        parent_indices = np.array([4, 5, 4, 5, 4, 5, 4, 5, -1, -1, -1, -1]).reshape(
            -1, 2
        )
        times = np.array([0, 0, 0, 0, 1, 1])
        is_sample = np.array([1, 1, 1, 1, 0, 0])
        t = max(times)

        model = msprime.WrightFisherPedigree()
        ped = msprime.Pedigree(
            inds, parent_indices, times, is_sample, sex=None, ploidy=2
        )
        ts = msprime.simulate(
            sample_size=4,
            pedigree=ped,
            recombination_rate=0.1,
            model=[model, (1, "dtwf")],
        )
        tree = ts.first()
        assert tree.num_roots == 1
        all_times = ts.tables.nodes.time
        ped_times = all_times[np.logical_and(all_times > 0, all_times <= t)]
        assert ped_times.shape[0] > 0
        assert np.all(ped_times == np.floor(ped_times))
        wf_times = all_times[all_times > t]
        assert wf_times.shape[0] > 0

    def test_wf_hudson_single_locus(self):
        Ne = 100
        t = 10
        ts = msprime.simulate(
            sample_size=10, Ne=Ne, model=["dtwf", (t, "hudson")], random_seed=2
        )
        tree = ts.first()
        assert tree.num_roots == 1
        times = ts.tables.nodes.time
        dtwf_times = times[np.logical_and(times > 0, times < t)]
        assert dtwf_times.shape[0] > 0
        assert np.all(dtwf_times == np.floor(dtwf_times))
        coalescent_times = times[times > t]
        assert coalescent_times.shape[0] > 0

    def test_wf_hudson_ancient_samples(self):
        Ne = 10
        t = 10
        n = 20
        ts = msprime.simulate(
            samples=[msprime.Sample(time=j, population=0) for j in range(n)],
            Ne=Ne,
            model=["dtwf", (t, "hudson")],
            random_seed=2,
        )
        tree = ts.first()
        assert tree.num_roots == 1
        times = ts.tables.nodes.time[ts.tables.nodes.flags == 0]
        dtwf_times = times[np.logical_and(times > 0, times < t)]
        assert dtwf_times.shape[0] > 0
        assert np.all(dtwf_times == np.floor(dtwf_times))
        coalescent_times = times[times > t]
        assert coalescent_times.shape[0] > 0
        assert np.all(coalescent_times != np.floor(coalescent_times))

    def test_wf_hudson_recombinatation(self):
        Ne = 100
        t = 100
        ts = msprime.simulate(
            sample_size=10,
            Ne=Ne,
            model=["dtwf", (t, "hudson")],
            recombination_rate=0.1,
            random_seed=2,
        )
        tree = ts.first()
        assert tree.num_roots == 1
        times = ts.tables.nodes.time
        dtwf_times = times[np.logical_and(times > 0, times < t)]
        assert dtwf_times.shape[0] > 0
        assert np.all(dtwf_times == np.floor(dtwf_times))
        coalescent_times = times[times > t]
        assert coalescent_times.shape[0] > 0
        assert np.all(coalescent_times != np.floor(coalescent_times))

    def test_wf_hudson_different_specifications(self):
        Ne = 100
        t = 100
        ts1 = msprime.simulate(
            sample_size=10,
            Ne=Ne,
            model=["dtwf", (t, "hudson")],
            recombination_rate=0.1,
            random_seed=2,
        )
        ts2 = msprime.simulate(
            sample_size=10,
            recombination_rate=0.1,
            Ne=Ne,
            model="dtwf",
            demographic_events=[msprime.SimulationModelChange(t, "hudson")],
            random_seed=2,
        )
        ts3 = msprime.simulate(
            sample_size=10,
            recombination_rate=0.1,
            Ne=Ne,
            model="dtwf",
            demographic_events=[msprime.SimulationModelChange(t)],
            random_seed=2,
        )
        t1 = ts1.dump_tables()
        t2 = ts2.dump_tables()
        t3 = ts3.dump_tables()
        t1.provenances.clear()
        t2.provenances.clear()
        t3.provenances.clear()
        assert t1 == t2
        assert t1 == t3

    def test_wf_hudson_back_and_forth(self):
        Ne = 100
        t1 = 100
        t2 = 200
        ts = msprime.simulate(
            sample_size=10,
            Ne=Ne,
            model=["dtwf", (t1, "hudson"), (t2, "dtwf")],
            recombination_rate=0.1,
            random_seed=2,
        )
        tree = ts.first()
        assert tree.num_roots == 1
        times = ts.tables.nodes.time
        dtwf_times = times[np.logical_and(times > 0, times < t1, times > t2)]
        assert dtwf_times.shape[0] > 0
        assert np.all(dtwf_times == np.floor(dtwf_times))
        coalescent_times = times[np.logical_and(times > t1, times < t2)]
        assert coalescent_times.shape[0] > 0
        assert np.all(coalescent_times != np.floor(coalescent_times))

    def test_many_models_simulate(self):
        Ne = 10000
        ts = msprime.simulate(
            Ne=Ne,
            sample_size=10,
            # Use the old-style SimulationModelChange
            model=[
                "hudson",
                msprime.SimulationModelChange(10, msprime.StandardCoalescent()),
                msprime.SimulationModelChange(20, msprime.SmcApproxCoalescent()),
                msprime.SimulationModelChange(30, msprime.SmcPrimeApproxCoalescent()),
                msprime.SimulationModelChange(40, msprime.DiscreteTimeWrightFisher()),
                msprime.SimulationModelChange(50, msprime.BetaCoalescent(alpha=1.1)),
                msprime.SimulationModelChange(60, msprime.StandardCoalescent()),
            ],
            random_seed=10,
        )
        for tree in ts.trees():
            assert tree.num_roots == 1

    def test_many_models_sim_ancestry(self):
        ts = msprime.sim_ancestry(
            samples=10,
            population_size=10_000,
            model=[
                "hudson",
                msprime.AncestryModelChange(10, msprime.StandardCoalescent()),
                msprime.AncestryModelChange(20, msprime.SmcApproxCoalescent()),
                msprime.AncestryModelChange(30, msprime.SmcPrimeApproxCoalescent()),
                msprime.AncestryModelChange(40, msprime.DiscreteTimeWrightFisher()),
                msprime.AncestryModelChange(50, msprime.BetaCoalescent(alpha=1.1)),
                msprime.AncestryModelChange(60, msprime.StandardCoalescent()),
            ],
            random_seed=10,
        )
        for tree in ts.trees():
            assert tree.num_roots == 1

    def test_too_many_models(self):
        # What happens when we have loads of models
        model_names = ["hudson", "smc"]
        models = ["hudson"]
        for j in range(1000):
            models.append(
                msprime.SimulationModelChange(time=j, model=model_names[j % 2])
            )
        ts = msprime.simulate(10, model=models, random_seed=2)
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_models_demographic_events(self):
        Ne = 10000
        ts = msprime.simulate(
            Ne=Ne,
            sample_size=10,
            recombination_rate=0.1,
            model=[
                None,
                msprime.SimulationModelChange(10, msprime.StandardCoalescent()),
            ],
            demographic_events=[
                msprime.SimpleBottleneck(11, population=0, proportion=1.0),
            ],
            random_seed=10,
        )
        for tree in ts.trees():
            assert tree.num_roots == 1
            assert ts.node(tree.root).time == 11

    def test_models_out_of_order(self):
        with pytest.raises(ValueError):
            msprime.simulate(
                Ne=10 ** 6,
                sample_size=10,
                model=[
                    "hudson",
                    msprime.SimulationModelChange(10, "hudson"),
                    msprime.SimulationModelChange(8, "hudson"),
                ],
                random_seed=2,
            )

    def test_model_change_negative_time(self):
        with pytest.raises(ValueError):
            msprime.simulate(
                Ne=10,
                sample_size=10,
                model=[None, msprime.SimulationModelChange(-10, "hudson")],
            )

    def test_model_change_time_bad_func(self):
        def bad_func(t):
            return t - 1

        with pytest.raises(ValueError):
            msprime.simulate(
                Ne=10,
                sample_size=10,
                model=[
                    None,
                    msprime.SimulationModelChange(1, "hudson"),
                    msprime.SimulationModelChange(bad_func, "hudson"),
                ],
            )


class TestSweepGenicSelection:
    """
    Tests for the single sweep model.
    """

    def assertTreeSequencesEqual(self, ts1, ts2):
        t1 = ts1.dump_tables()
        t2 = ts2.dump_tables()
        t1.provenances.clear()
        t2.provenances.clear()
        assert t1 == t2

    def test_incorrect_num_labels(self):
        model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.1, end_frequency=0.9, s=0.01, dt=0.01
        )
        for num_labels in [1, 3, 10]:
            with pytest.raises(_msprime.LibraryError):
                msprime.simulate(
                    10, recombination_rate=1, model=model, num_labels=num_labels
                )

    def test_sweep_coalescence_no_recomb(self):
        N = 1e6
        model = msprime.SweepGenicSelection(
            position=0.5,
            start_frequency=1.0 / (2 * N),
            end_frequency=1.0 - (1.0 / (2 * N)),
            s=0.1,
            dt=1e-6,
        )
        ts = msprime.simulate(
            10, model=model, Ne=1000, recombination_rate=0, random_seed=2
        )
        assert ts.num_trees == 1
        for tree in ts.trees():
            assert tree.num_roots == 1

    def test_sweep_coalescence_recomb(self):
        N = 1e6
        model = msprime.SweepGenicSelection(
            position=0.5,
            start_frequency=1.0 / (2 * N),
            end_frequency=1.0 - (1.0 / (2 * N)),
            s=0.1,
            dt=1e-6,
        )
        ts = msprime.simulate(
            10, model=model, Ne=1000, recombination_rate=1, random_seed=2
        )
        assert ts.num_trees > 1

    def test_sweep_coalescence_same_seed(self):
        model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.6, end_frequency=0.7, s=0.1, dt=1e-6
        )
        ts1 = msprime.simulate(5, model=model, random_seed=2)
        ts2 = msprime.simulate(5, model=model, random_seed=2)
        self.assertTreeSequencesEqual(ts1, ts2)

    def test_sweep_start_time_complete(self):
        sweep_model = msprime.SweepGenicSelection(
            position=0.5,
            start_frequency=0.01,
            end_frequency=0.99,
            s=0.25,
            dt=1e-6,
        )
        t_start = 0.1
        ts = msprime.simulate(
            10,
            Ne=1000,
            recombination_rate=2,
            model=[None, (t_start, sweep_model), (None, None)],
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_sweep_start_time_incomplete(self):
        # Short sweep that doesn't make complete coalescence.
        sweep_model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.69, end_frequency=0.7, s=0.01, dt=1e-6
        )
        t_start = 0.1
        ts = msprime.simulate(
            10,
            Ne=1000,
            recombination_rate=2,
            model=["hudson", (t_start, sweep_model)],
            random_seed=2,
        )
        assert any(tree.num_roots > 1 for tree in ts.trees())

    def test_sweep_model_change_time_complete(self):
        # Short sweep that doesn't coalesce followed
        # by Hudson phase to finish up coalescent
        sweep_model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.69, end_frequency=0.72, s=0.1, dt=1e-6
        )
        ts = msprime.simulate(
            10,
            Ne=1000,
            recombination_rate=2,
            model=[sweep_model, (None, None)],
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())

        # Returning None from a function should be identical
        ts2 = msprime.simulate(
            10,
            Ne=1000,
            recombination_rate=2,
            model=[
                sweep_model,
                msprime.SimulationModelChange(lambda t: None, "hudson"),
            ],
            random_seed=2,
        )
        self.assertTreeSequencesEqual(ts, ts2)

        # Make sure that the Hudson phase did something.
        ts = msprime.simulate(
            10, Ne=1000, recombination_rate=2, model=sweep_model, random_seed=2
        )
        assert any(tree.num_roots > 1 for tree in ts.trees())

    def test_many_sweeps(self):
        sweep_models = [
            msprime.SweepGenicSelection(
                position=j, start_frequency=0.69, end_frequency=0.7, s=0.001, dt=1e-6
            )
            for j in range(10)
        ]
        ts = msprime.simulate(
            10,
            Ne=1000,
            length=10,
            recombination_rate=0.2,
            model=[None, msprime.SimulationModelChange(0.01, sweep_models[0])]
            + [msprime.SimulationModelChange(None, model) for model in sweep_models]
            + [msprime.SimulationModelChange()],
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_many_sweeps_regular_times_model_change(self):
        events = []
        for j in range(10):
            sweep_model = msprime.SweepGenicSelection(
                position=j,
                start_frequency=0.69,
                end_frequency=0.7,
                s=0.1,
                dt=1e-6,
            )
            # Start the sweep after 0.01 generations of Hudson
            events.append(
                msprime.SimulationModelChange(
                    time=lambda t: t + 0.01, model=sweep_model
                )
            )
            # Revert back to Hudson until the next sweep
            events.append(msprime.SimulationModelChange())
        # The old-style approach
        ts = msprime.simulate(
            10,
            Ne=1000,
            length=10,
            recombination_rate=0.2,
            demographic_events=events,
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())
        # The recommended way.
        ts2 = msprime.simulate(
            10,
            Ne=1000,
            model=["hudson"] + events,
            length=10,
            recombination_rate=0.2,
            random_seed=2,
        )
        self.assertTreeSequencesEqual(ts2, ts)

    def test_too_many_sweeps(self):
        # What happens when we have loads of sweeps
        events = []
        for _ in range(1000):
            sweep_model = msprime.SweepGenicSelection(
                position=5,
                start_frequency=0.69,
                end_frequency=0.7,
                s=0.1,
                dt=1e-6,
            )
            # Start the sweep after 0.1 generations of Hudson
            events.append(
                msprime.SimulationModelChange(time=lambda t: t + 0.1, model=sweep_model)
            )
            # Revert back to Hudson until the next sweep
            events.append(msprime.SimulationModelChange())
        # Old style
        ts = msprime.simulate(
            10,
            Ne=100,
            length=10,
            recombination_rate=0.2,
            demographic_events=events,
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())
        # New style
        ts2 = msprime.simulate(
            10,
            model=[None] + events,
            Ne=100,
            length=10,
            recombination_rate=0.2,
            random_seed=2,
        )
        self.assertTreeSequencesEqual(ts2, ts)
