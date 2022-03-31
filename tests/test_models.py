#
# Copyright (C) 2016-2021 University of Oxford
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


nonparametric_model_classes = [
    msprime.StandardCoalescent,
    msprime.SmcApproxCoalescent,
    msprime.SmcPrimeApproxCoalescent,
    msprime.DiscreteTimeWrightFisher,
    msprime.FixedPedigree,
]


class TestIntrospectionInterface:
    """
    Tests that we have meaningful repr and str functions for all the
    classes used in the model hierarchy.
    """

    def test_standard_coalescent(self):
        model = msprime.StandardCoalescent()
        repr_s = "StandardCoalescent(duration=None)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_smc_models(self):
        model = msprime.SmcApproxCoalescent()
        repr_s = "SmcApproxCoalescent(duration=None)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

        model = msprime.SmcPrimeApproxCoalescent()
        repr_s = "SmcPrimeApproxCoalescent(duration=None)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_dtwf(self):
        model = msprime.DiscreteTimeWrightFisher()
        repr_s = "DiscreteTimeWrightFisher(duration=None)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_fixed_pedigreeigree(self):
        model = msprime.FixedPedigree()
        repr_s = "FixedPedigree(duration=None)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_beta_coalescent(self):
        model = msprime.BetaCoalescent(alpha=1, truncation_point=2)
        repr_s = "BetaCoalescent(duration=None, alpha=1, truncation_point=2)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_dirac_coalescent(self):
        model = msprime.DiracCoalescent(psi=123, c=2)
        repr_s = "DiracCoalescent(duration=None, psi=123, c=2)"
        assert repr(model) == repr_s
        assert str(model) == repr_s

    def test_sweep_genic_selection(self):
        model = msprime.SweepGenicSelection(
            position=1, start_frequency=0.5, end_frequency=0.9, s=0.1, dt=0.01
        )
        repr_s = (
            "SweepGenicSelection(duration=None, position=1, start_frequency=0.5, "
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
            ("fixed_pedigree", msprime.FixedPedigree),
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
            msprime.FixedPedigree(),
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
        models = ancestry._parse_model_arg(None)
        assert models == [msprime.StandardCoalescent()]

    def test_single_model(self):
        models = ancestry._parse_model_arg("hudson")
        assert models == [msprime.StandardCoalescent()]

        models = ancestry._parse_model_arg(msprime.StandardCoalescent())
        assert models == [msprime.StandardCoalescent()]

        models = ancestry._parse_model_arg("dtwf")
        assert models == [msprime.DiscreteTimeWrightFisher()]

        models = ancestry._parse_model_arg(msprime.DiscreteTimeWrightFisher())
        assert models == [msprime.DiscreteTimeWrightFisher()]

    def test_single_model_list(self):
        models = ancestry._parse_model_arg([None])
        assert models == [msprime.StandardCoalescent()]

        # Tuples are also accepted as input.
        models = ancestry._parse_model_arg((None,))
        assert models == [msprime.StandardCoalescent()]

        models = ancestry._parse_model_arg(["hudson"])
        assert models == [msprime.StandardCoalescent()]

        models = ancestry._parse_model_arg([msprime.StandardCoalescent()])
        assert models == [msprime.StandardCoalescent()]

        models = ancestry._parse_model_arg(["dtwf"])
        assert models == [msprime.DiscreteTimeWrightFisher()]

        models = ancestry._parse_model_arg([msprime.DiscreteTimeWrightFisher()])
        assert models == [msprime.DiscreteTimeWrightFisher()]

    def test_one_change(self):
        models = ancestry._parse_model_arg(
            [msprime.DiscreteTimeWrightFisher(duration=1.33), "hudson"]
        )
        assert models == [
            msprime.DiscreteTimeWrightFisher(duration=1.33),
            msprime.StandardCoalescent(),
        ]

        models = ancestry._parse_model_arg(
            [msprime.DiscreteTimeWrightFisher(duration=1.33), None]
        )
        assert models == [
            msprime.DiscreteTimeWrightFisher(duration=1.33),
            msprime.StandardCoalescent(),
        ]

    def test_two_changes(self):
        models = ancestry._parse_model_arg(
            [
                msprime.DiscreteTimeWrightFisher(duration=1),
                msprime.StandardCoalescent(duration=2),
                msprime.SmcApproxCoalescent(duration=3),
            ]
        )
        assert models == [
            msprime.DiscreteTimeWrightFisher(duration=1),
            msprime.StandardCoalescent(duration=2),
            msprime.SmcApproxCoalescent(duration=3),
        ]

    def test_errors(self):
        with pytest.raises(ValueError, match="at least one AncestryModel"):
            ancestry._parse_model_arg([])
        for bad_model in ["X", ["X"]]:
            with pytest.raises(ValueError, match="Model 'X' unknown"):
                ancestry._parse_model_arg(bad_model)

        # Anything that's not a list or tuple is interpreted as a model
        with pytest.raises(TypeError, match="Ancestry model must be a string or"):
            ancestry._parse_model_arg({})


class TestClassesKeywordArgs:
    @pytest.mark.parametrize("cls", nonparametric_model_classes)
    def test_non_parametric(self, cls):
        model = cls()
        assert model.duration is None
        model = cls(duration=1)
        assert model.duration == 1
        with pytest.raises(TypeError, match="takes 1 positional"):
            cls(1)

    def test_beta_coalescent(self):
        model = msprime.BetaCoalescent(alpha=1, truncation_point=2)
        assert model.duration is None
        assert model.alpha == 1
        assert model.truncation_point == 2

        model = msprime.BetaCoalescent(alpha=2, truncation_point=1, duration=3)
        assert model.duration == 3
        assert model.alpha == 2
        assert model.truncation_point == 1

        model = msprime.BetaCoalescent()
        assert model.duration is None
        assert model.alpha is None

        with pytest.raises(TypeError, match="takes 1 positional"):
            msprime.BetaCoalescent(1)

    def test_dirac_coalescent(self):
        model = msprime.DiracCoalescent(psi=123, c=2)
        assert model.duration is None
        assert model.psi == 123
        assert model.c == 2

        model = msprime.DiracCoalescent(psi=123, c=2, duration=3)
        assert model.duration == 3
        assert model.psi == 123
        assert model.c == 2

        model = msprime.DiracCoalescent()
        assert model.duration is None
        assert model.psi is None

        with pytest.raises(TypeError, match="takes 1 positional"):
            msprime.DiracCoalescent(1)

    def test_sweep_genic_selection(self):
        model = msprime.SweepGenicSelection(
            position=1, start_frequency=0.5, end_frequency=0.9, s=0.1, dt=0.01
        )
        assert model.duration is None
        assert model.position == 1
        assert model.start_frequency == 0.5
        assert model.end_frequency == 0.9
        assert model.s == 0.1
        assert model.dt == 0.01

        model = msprime.SweepGenicSelection(
            position=2,
            start_frequency=0.5,
            end_frequency=0.9,
            s=0.1,
            dt=0.01,
            duration=1234,
        )
        assert model.duration == 1234
        assert model.position == 2
        assert model.start_frequency == 0.5
        assert model.end_frequency == 0.9
        assert model.s == 0.1
        assert model.dt == 0.01

        model = msprime.SweepGenicSelection()
        assert model.duration is None
        assert model.position is None
        assert model.start_frequency is None
        assert model.end_frequency is None
        assert model.s is None
        assert model.dt is None

        with pytest.raises(TypeError, match="takes 1 positional"):
            msprime.SweepGenicSelection(1)


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
            model = msprime.BetaCoalescent(alpha=alpha, truncation_point=1)
            assert model.alpha == alpha
            assert model.truncation_point == 1
            d = model._as_lowlevel()
            assert d == {
                "name": "beta",
                "duration": None,
                "alpha": alpha,
                "truncation_point": 1,
            }
        alpha = 1.5
        for truncation_point in [0.01, 0.5, 1]:
            model = msprime.BetaCoalescent(
                alpha=alpha, truncation_point=truncation_point
            )
            assert model.alpha == alpha
            assert model.truncation_point == truncation_point
            d = model._as_lowlevel()
            assert d == {
                "name": "beta",
                "alpha": alpha,
                "truncation_point": truncation_point,
                "duration": None,
            }

    def test_dirac_coalescent_parameters(self):
        for psi in [0.01, 0.5, 0.99]:
            for c in [1e-6, 1.0, 1e2]:
                model = msprime.DiracCoalescent(psi=psi, c=c)
                assert model.psi == psi
                assert model.c == c
                d = model._as_lowlevel()
                assert d == {"name": "dirac", "psi": psi, "c": c, "duration": None}


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
        model = msprime.DiracCoalescent(psi=0.3, c=10)
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
        recombination_map = msprime.RateMap(position=[0, 100, 101, 200], rate=[0, 1, 0])
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

    def test_wf_hudson_single_locus(self):
        t = 10
        ts = msprime.sim_ancestry(
            10,
            population_size=100,
            model=[msprime.DiscreteTimeWrightFisher(duration=t), "hudson"],
            random_seed=3,
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
        t = 10
        n = 20
        ts = msprime.sim_ancestry(
            [msprime.SampleSet(1, time=j, population=0) for j in range(n)],
            population_size=10,
            model=[msprime.DiscreteTimeWrightFisher(duration=t), "hudson"],
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
        t = 100
        ts = msprime.sim_ancestry(
            10,
            population_size=100,
            model=[msprime.DiscreteTimeWrightFisher(duration=t), "hudson"],
            recombination_rate=0.1,
            sequence_length=10,
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
        ts1 = msprime.sim_ancestry(
            samples=5,
            population_size=Ne,
            model=[msprime.DiscreteTimeWrightFisher(duration=t), "hudson"],
            recombination_rate=0.1,
            sequence_length=1,
            discrete_genome=False,
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
        # Not worth trying to puzzle out the slight differences in tables
        # between the old and new form. The edges are the same, good enough.
        assert ts1.tables.edges == ts2.tables.edges
        assert ts2.equals(ts3, ignore_provenance=True)

    def test_wf_hudson_back_and_forth(self):
        t1 = 100
        t2 = 200
        ts = msprime.sim_ancestry(
            5,
            population_size=100,
            model=[
                msprime.DiscreteTimeWrightFisher(duration=t1),
                msprime.StandardCoalescent(duration=t2 - t1),
                "dtwf",
            ],
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
            model="hudson",
            demographic_events=[
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
                msprime.StandardCoalescent(duration=10),
                msprime.SmcApproxCoalescent(duration=10),
                msprime.SmcPrimeApproxCoalescent(duration=10),
                msprime.DiscreteTimeWrightFisher(duration=10),
                msprime.BetaCoalescent(alpha=1.1, duration=10),
                msprime.StandardCoalescent(),
            ],
            random_seed=10,
        )
        for tree in ts.trees():
            assert tree.num_roots == 1

    def test_many_models(self):
        # What happens when we have loads of models
        models = [
            msprime.StandardCoalescent(duration=0.1),
            msprime.SmcApproxCoalescent(duration=0.1),
        ]
        ts = msprime.sim_ancestry(10, model=models * 1000, random_seed=2)
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_models_demographic_events(self):
        Ne = 10000
        ts = msprime.simulate(
            Ne=Ne,
            sample_size=10,
            recombination_rate=0.1,
            model="smc",
            demographic_events=[
                msprime.SimulationModelChange(10, msprime.StandardCoalescent()),
                msprime.SimpleBottleneck(11, population=0, proportion=1.0),
            ],
            random_seed=10,
        )
        for tree in ts.trees():
            assert tree.num_roots == 1
            assert ts.node(tree.root).time == 11

    def test_models_out_of_order(self):
        with pytest.raises(ValueError, match="durations must be >= 0"):
            msprime.simulate(
                Ne=10**6,
                sample_size=10,
                model="hudson",
                demographic_events=[
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
                model=None,
                demographic_events=[msprime.SimulationModelChange(-10, "hudson")],
            )


class TestSweepGenicSelection:
    """
    Tests for the single sweep model.
    """

    def test_model_end_broken(self):
        # Checking that we're correctly detecting the fact that
        # sweeps are non renentrant.
        model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.1, end_frequency=0.9, s=0.01, dt=0.01
        )
        with pytest.raises(RuntimeError, match="does not support interruption"):
            msprime.sim_ancestry(10, model=model, end_time=0.0001)

    def test_incorrect_num_labels(self):
        model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.1, end_frequency=0.9, s=0.01, dt=0.01
        )
        for num_labels in [1, 3, 10]:
            # Not the best error, but this shouldn't be exposed to the user anyway.
            with pytest.raises(
                _msprime.LibraryError, match="configuration is not supported"
            ):
                msprime.sim_ancestry(
                    10,
                    model=model,
                    num_labels=num_labels,
                )

    def test_float_position_discrete_genome(self):
        # We can still specify a floating point position if we want, doesn't
        # make any difference.
        model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.1, end_frequency=0.9, s=0.01, dt=0.01
        )
        ts = msprime.sim_ancestry(
            10, sequence_length=1, recombination_rate=20, model=model
        )
        assert ts.num_trees == 1

    def test_sweep_coalescence_no_recomb(self):
        N = 1e6
        model = msprime.SweepGenicSelection(
            position=0.5,
            start_frequency=1.0 / (2 * N),
            end_frequency=1.0 - (1.0 / (2 * N)),
            s=0.1,
            dt=1e-6,
        )
        ts = msprime.sim_ancestry(10, model=model, population_size=1000, random_seed=2)
        assert ts.num_trees == 1
        for tree in ts.trees():
            assert tree.num_roots == 1

    def test_sweep_coalescence_recomb(self):
        N = 1e6
        model = msprime.SweepGenicSelection(
            position=5,
            start_frequency=1.0 / (2 * N),
            end_frequency=1.0 - (1.0 / (2 * N)),
            s=0.1,
            dt=1e-6,
        )
        ts = msprime.sim_ancestry(
            5,
            model=model,
            population_size=1000,
            recombination_rate=1,
            sequence_length=10,
            random_seed=2,
        )
        assert ts.num_trees > 1

    def test_sweep_coalescence_same_seed(self):
        model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.6, end_frequency=0.7, s=0.1, dt=1e-6
        )
        ts1 = msprime.sim_ancestry(5, model=model, random_seed=2)
        ts2 = msprime.sim_ancestry(5, model=model, random_seed=2)
        assert ts1.equals(ts2, ignore_provenance=True)

    def test_sweep_start_time_complete(self):
        sweep_model = msprime.SweepGenicSelection(
            position=0.5,
            start_frequency=0.01,
            end_frequency=0.99,
            s=0.25,
            dt=1e-6,
        )
        ts = msprime.sim_ancestry(
            10,
            population_size=1000,
            recombination_rate=2,
            sequence_length=10,
            model=[msprime.StandardCoalescent(duration=0.1), sweep_model, "hudson"],
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_sweep_start_time_incomplete(self):
        # Short sweep that doesn't make complete coalescence.
        sweep_model = msprime.SweepGenicSelection(
            position=0.5, start_frequency=0.69, end_frequency=0.7, s=0.01, dt=1e-6
        )
        ts = msprime.sim_ancestry(
            10,
            population_size=1000,
            recombination_rate=2,
            sequence_length=10,
            model=[msprime.StandardCoalescent(duration=0.1), sweep_model],
            random_seed=3,
        )
        assert any(tree.num_roots > 1 for tree in ts.trees())

    def test_sweep_model_change_time_complete(self):
        # Short sweep that doesn't coalesce followed
        # by Hudson phase to finish up coalescent
        sweep_model = msprime.SweepGenicSelection(
            position=5, start_frequency=0.69, end_frequency=0.72, s=0.1, dt=1e-6
        )
        ts = msprime.sim_ancestry(
            10,
            population_size=1000,
            sequence_length=10,
            recombination_rate=2,
            model=[sweep_model, "hudson"],
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())

        ts = msprime.sim_ancestry(
            10,
            population_size=1000,
            recombination_rate=2,
            model=sweep_model,
            random_seed=2,
            sequence_length=10,
            discrete_genome=False,
        )
        assert any(tree.num_roots > 1 for tree in ts.trees())

    def test_many_sweeps(self):
        sweep_models = [
            msprime.SweepGenicSelection(
                position=j, start_frequency=0.69, end_frequency=0.7, s=0.001, dt=1e-6
            )
            for j in range(10)
        ]
        ts = msprime.sim_ancestry(
            10,
            population_size=1000,
            sequence_length=10,
            recombination_rate=0.2,
            model=sweep_models + ["hudson"],
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_many_sweeps_regular_times_model_change(self):
        models = []
        for j in range(0, 10):
            models.extend(
                [
                    # Start each sweep after 0.01 generations of Hudson
                    msprime.StandardCoalescent(duration=0.01),
                    msprime.SweepGenicSelection(
                        position=j,
                        start_frequency=0.69,
                        end_frequency=0.7,
                        s=0.1,
                        dt=1e-6,
                    ),
                ]
            )
        # Complete the simulation with Hudson
        models.append("hudson")
        ts = msprime.sim_ancestry(
            3,
            population_size=1000,
            sequence_length=10,
            recombination_rate=0.2,
            model=models,
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())

    def test_lots_of_sweeps(self):
        # What happens when we have loads of sweeps
        sweep_model = msprime.SweepGenicSelection(
            position=5,
            start_frequency=0.69,
            end_frequency=0.7,
            s=0.1,
            dt=1e-6,
        )
        # Start each sweep after 0.1 generations of Hudson
        models = [msprime.StandardCoalescent(duration=0.1), sweep_model] * 1000
        ts = msprime.sim_ancestry(
            10,
            population_size=100,
            sequence_length=10,
            recombination_rate=0.2,
            model=models,
            random_seed=2,
        )
        assert all(tree.num_roots == 1 for tree in ts.trees())
