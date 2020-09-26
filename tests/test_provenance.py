#
# Copyright (C) 2018-2020 University of Oxford
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
Tests for the provenance information attached to tree sequences.
"""
import base64
import inspect
import json
import marshal
import unittest

import numpy as np
import python_jsonschema_objects as pjs
import tskit

import msprime
from msprime import _msprime
from msprime import ancestry


class TestProvenance(unittest.TestCase):
    """
    Basic tests for the provenance dict function.
    """

    def test_libraries(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = json.loads(ts.provenance(0).record)
        libs = prov["environment"]["libraries"]
        self.assertEqual(
            libs["gsl"], {"version": ".".join(map(str, _msprime.get_gsl_version()))}
        )
        self.assertEqual(libs["tskit"], {"version": tskit.__version__})


class TestBuildProvenance(unittest.TestCase):
    """
    Tests for the provenance dictionary building. This dictionary is used
    to encode the parameters for the msprime simulations.
    """

    def test_basic(self):
        def somefunc(a, b):
            frame = inspect.currentframe()
            return ancestry._build_provenance("cmd", 1234, frame)

        d = somefunc(42, 43)
        tskit.validate_provenance(d)
        params = d["parameters"]
        self.assertEqual(params["command"], "cmd")
        self.assertEqual(params["random_seed"], 1234)
        self.assertEqual(params["a"], 42)
        self.assertEqual(params["b"], 43)

    def test_replicates(self):
        def somefunc(*, a, b, num_replicates, replicate_index):
            frame = inspect.currentframe()
            return ancestry._build_provenance("the_cmd", 42, frame)

        d = somefunc(b="b", a="a", num_replicates=100, replicate_index=1234)
        tskit.validate_provenance(d)
        params = d["parameters"]
        self.assertEqual(params["command"], "the_cmd")
        self.assertEqual(params["random_seed"], 42)
        self.assertEqual(params["a"], "a")
        self.assertEqual(params["b"], "b")
        self.assertFalse("num_replicates" in d)
        self.assertFalse("replicate_index" in d)


class ValidateSchemas(unittest.TestCase):
    """
    Check that the schemas we produce in msprime are valid.
    """

    def test_simulation(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = json.loads(ts.provenance(0).record)
        tskit.validate_provenance(prov)
        self.assertEqual(prov["parameters"]["command"], "simulate")

    def test_mutate(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(ts, rate=1, random_seed=1)
        prov = json.loads(ts.provenance(1).record)
        tskit.validate_provenance(prov)
        self.assertEqual(prov["parameters"]["command"], "mutate")

    def test_simplify(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = ts.simplify()
        prov = json.loads(ts.provenance(1).record)
        tskit.validate_provenance(prov)
        self.assertEqual(prov["parameters"]["command"], "simplify")


class TestBuildObjects(unittest.TestCase):
    """
    Check that we can build objects from the json schema as we'd expect.
    """

    def decode(self, prov):
        builder = pjs.ObjectBuilder(tskit.provenance.get_schema())
        ns = builder.build_classes()
        return ns.TskitProvenance.from_json(prov)

    def test_simulation(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = ts.provenance(0).record
        decoded = self.decode(prov)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "simulate")
        self.assertEqual(decoded.parameters.random_seed, 1)

    def test_no_provenance(self):
        ts = msprime.simulate(5, record_provenance=False)
        self.assertEqual(len(list(ts.provenances())), 0)

    def test_encode_simulation_models(self):
        simple_model = ["hudson", [10, "dtwf"], [20, "smc"], [None, None]]
        ts = msprime.simulate(10, model=simple_model)
        decoded = self.decode(ts.provenance(0).record)
        parameters = decoded.parameters
        self.assertEqual(parameters.sample_size, 10)
        self.assertEqual(list(parameters.model), simple_model)

        model_instances = [
            msprime.StandardCoalescent(),
            msprime.SimulationModelChange(10, msprime.DiscreteTimeWrightFisher()),
            msprime.SimulationModelChange(20, msprime.SmcApproxCoalescent()),
            msprime.SimulationModelChange(30, msprime.BetaCoalescent(alpha=1.1)),
        ]
        ts = msprime.simulate(10, model=model_instances)
        decoded = self.decode(ts.provenance(0).record)
        parameters = decoded.parameters
        self.assertEqual(parameters.sample_size, 10)
        self.assertEqual(
            parameters.model[0], {"__class__": "msprime.ancestry.StandardCoalescent"}
        )
        self.assertDictEqual(
            parameters.model[1],
            {
                "__class__": "msprime.ancestry.SimulationModelChange",
                "model": {"__class__": "msprime.ancestry.DiscreteTimeWrightFisher"},
                "time": 10,
            },
        )
        self.assertDictEqual(
            parameters.model[2],
            {
                "__class__": "msprime.ancestry.SimulationModelChange",
                "model": {"__class__": "msprime.ancestry.SmcApproxCoalescent"},
                "time": 20,
            },
        )
        self.assertDictEqual(
            parameters.model[3],
            {
                "__class__": "msprime.ancestry.SimulationModelChange",
                "model": {
                    "__class__": "msprime.ancestry.BetaCoalescent",
                    "alpha": 1.1,
                    "truncation_point": 1.0,
                },
                "time": 30,
            },
        )

    def test_encode_numpy_functions_and_classes(self):
        seeds = np.ones(1, dtype=int)
        pop_configs = [msprime.PopulationConfiguration(5) for _ in range(2)]

        def time(t):
            return t + 0.01

        def bad_func(t):
            if msprime:
                pass
            return t

        with self.assertLogs("msprime.provenance", level="WARNING") as cm:
            ts = msprime.simulate(
                random_seed=seeds[0],
                population_configurations=pop_configs,
                migration_matrix=[[0, 1], [1, 0]],
                demographic_events=[
                    msprime.SimulationModelChange(time=time),
                    msprime.SimulationModelChange(time=bad_func),
                    msprime.SimulationModelChange(time=lambda t: t + seeds[0]),
                ],
            )
        self.assertIn(
            "WARNING:msprime.provenance:Configuration function bad_func refers"
            " to global variables ('msprime',) which are not currently"
            " serialized",
            cm.output,
        )
        self.assertIn(
            "WARNING:msprime.provenance:Configuration function <lambda> refers"
            " to outer scope variables which are not currently serialized",
            cm.output,
        )
        prov = ts.provenance(0).record
        decoded = self.decode(prov)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "simulate")
        parameters = decoded.parameters
        self.assertEqual(parameters.random_seed, 1)
        self.assertEqual(
            parameters.population_configurations,
            [
                {
                    "sample_size": 5,
                    "initial_size": None,
                    "growth_rate": 0.0,
                    "metadata": None,
                    "__class__": "msprime.demography.PopulationConfiguration",
                },
                {
                    "sample_size": 5,
                    "initial_size": None,
                    "growth_rate": 0.0,
                    "metadata": None,
                    "__class__": "msprime.demography.PopulationConfiguration",
                },
            ],
        )
        self.assertEqual(parameters.migration_matrix, [[0, 1], [1, 0]])
        self.assertDictEqual(
            parameters.demographic_events[0],
            {
                "time": {
                    "__function__": base64.b64encode(
                        marshal.dumps(time.__code__)
                    ).decode("utf-8"),
                    "defaults": None,
                },
                "model": None,
                "__class__": "msprime.ancestry.SimulationModelChange",
            },
        )
        self.assertDictEqual(
            parameters.demographic_events[1],
            {
                "time": {
                    "__error__": "Configuration function bad_func "
                    "refers to global variables"
                    " ('msprime',) which are not currently serialized"
                },
                "model": None,
                "__class__": "msprime.ancestry.SimulationModelChange",
            },
        ),
        self.assertDictEqual(
            parameters.demographic_events[2],
            {
                "time": {
                    "__error__": "Configuration function <lambda> refers to outer scope"
                    " variables which are not currently serialized"
                },
                "model": None,
                "__class__": "msprime.ancestry.SimulationModelChange",
            },
        )

    def test_unencodable_function(self):
        class BadClass:
            def __call__(self, t):
                return t

        with self.assertRaises(TypeError):
            msprime.simulate(
                5,
                random_seed=1,
                demographic_events=[msprime.SimulationModelChange(time=BadClass)],
            )

    def test_replicate_index(self):
        for i, ts in enumerate(msprime.simulate(5, num_replicates=3)):
            decoded = self.decode(ts.provenance(0).record)
            self.assertEqual(decoded.parameters["replicate_index"], i)

    def test_large_provenance_warning(self):
        with self.assertLogs("msprime.provenance", level="WARNING") as cm:
            msprime.provenance.json_encode_provenance(
                {"test": np.zeros(shape=(500_000,))}, 5
            )
        self.assertIn(
            "WARNING:msprime.provenance:The provenance information for the resulting"
            " tree sequence is 2.38MB. This is nothing to worry about as"
            " provenance is a good thing to have, but if you want to save this"
            " memory/storage space you can disable provenance recording by setting"
            " record_provenance=False",
            cm.output,
        )

    def test_mutate(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(
            ts, rate=2, random_seed=1, start_time=0, end_time=100, keep=False
        )
        decoded = self.decode(ts.provenance(1).record)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "mutate")
        self.assertEqual(decoded.parameters.random_seed, 1)
        self.assertEqual(decoded.parameters.rate, 2)
        self.assertEqual(decoded.parameters.start_time, 0)
        self.assertEqual(decoded.parameters.end_time, 100)
        self.assertEqual(decoded.parameters.keep, False)
        self.assertEqual(
            decoded.parameters.model["__class__"],
            "msprime.mutations.BinaryMutationModel",
        )

    def test_mutate_model(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(ts, model="jc69")
        decoded = self.decode(ts.provenance(1).record)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "mutate")
        self.assertEqual(
            decoded.parameters.model["__class__"], "msprime.mutations.JC69MutationModel"
        )

    def test_mutate_map(self):
        ts = msprime.simulate(5, random_seed=1)
        rate_map = msprime.RateMap(position=[0, 0.5, 1], rate=[0, 1])
        ts = msprime.mutate(ts, rate=rate_map)
        decoded = self.decode(ts.provenance(1).record)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "mutate")
        self.assertEqual(
            decoded.parameters.rate["__class__"], "msprime.intervals.RateMap"
        )
        self.assertEqual(
            decoded.parameters.rate["position"]["__ndarray__"], list(rate_map.position)
        )
        self.assertEqual(
            decoded.parameters.rate["rate"]["__ndarray__"], list(rate_map.rate)
        )

    def test_mutate_numpy(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(
            ts,
            rate=np.array([2])[0],
            random_seed=np.array([1])[0],
            start_time=np.array([0])[0],
            end_time=np.array([100][0]),
            keep=np.array([False][0]),
        )
        decoded = self.decode(ts.provenance(1).record)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "mutate")
        self.assertEqual(decoded.parameters.random_seed, 1)
        self.assertEqual(decoded.parameters.rate, 2)
        self.assertEqual(decoded.parameters.start_time, 0)
        self.assertEqual(decoded.parameters.end_time, 100)
        self.assertEqual(decoded.parameters.keep, False)


class TestParseProvenance(unittest.TestCase):
    def test_parse_non_msprime(self):
        # Most testing of parsing is in the round tripping below
        with self.assertRaises(ValueError):
            ts = msprime.simulate(5)
            prov = json.loads(ts.provenance(0).record)
            prov["software"]["name"] = "skynet"
            tables = ts.tables
            tables.provenances.clear()
            tables.provenances.add_row(json.dumps(prov))
            command, prov = msprime.provenance.parse_provenance(
                tables.tree_sequence().provenance(0), None
            )

    def test_current_ts(self):
        ts1 = msprime.simulate(5, random_seed=1)
        ts2 = msprime.mutate(ts1)
        command, prov = msprime.provenance.parse_provenance(ts2.provenance(1), ts1)
        self.assertEquals(command, "mutate")
        self.assertEquals(prov["tree_sequence"], ts1)


class TestRoundTrip(unittest.TestCase):
    """
    Check that we can recreate tree sequences from their provenances
    """

    def verify_equal(self, ts1, ts2):
        ts1_tables = ts1.dump_tables()
        ts2_tables = ts2.dump_tables()
        # We clear the provenances as they are not equal due to timestamps
        ts1_tables.provenances.clear()
        ts2_tables.provenances.clear()
        self.assertEqual(ts1_tables, ts2_tables)

    def verify(self, ts, recurse=True):
        recreated = None
        for provenance in ts.provenances():
            command, parsed_prov = msprime.provenance.parse_provenance(
                provenance, recreated
            )
            recreated = getattr(msprime, command)(**parsed_prov)
        self.verify_equal(ts, recreated)
        # We verify again to check that the recreated provenace is also good.
        if recurse:
            self.verify(recreated, recurse=False)


class TestSimulateRoundTrip(TestRoundTrip):
    def test_sample_size(self):
        ts = msprime.simulate(10)
        self.verify(ts)

    def test_random_seed(self):
        ts = msprime.simulate(10, random_seed=1234)
        self.verify(ts)

    def test_population_configuration(self):
        pop_configs = [msprime.PopulationConfiguration(5) for _ in range(2)]
        ts = msprime.simulate(
            population_configurations=pop_configs,
            migration_matrix=[[0, 1], [1, 0]],
            demographic_events=[msprime.SimulationModelChange(time=lambda t: t + 0.1)],
        )
        self.verify(ts)

    def test_simulation_models(self):
        simple_model = ["hudson", [10, "dtwf"], [20, "smc"]]
        ts = msprime.simulate(10, model=simple_model)
        self.verify(ts)

        model_instances = [
            msprime.StandardCoalescent(),
            msprime.SimulationModelChange(10, msprime.DiscreteTimeWrightFisher()),
            msprime.SimulationModelChange(20, msprime.SmcApproxCoalescent()),
            msprime.SimulationModelChange(30, msprime.BetaCoalescent(alpha=1.1)),
        ]
        ts = msprime.simulate(10, model=model_instances)
        self.verify(ts)

    def test_pedigree(self):
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
            demographic_events=[
                msprime.SimulationModelChange(t, msprime.DiscreteTimeWrightFisher())
            ],
            model=model,
        )
        self.verify(ts)

    def test_from_ts(self):
        from_ts = msprime.simulate(10, random_seed=5, end_time=0.5)
        ts = msprime.simulate(from_ts=from_ts, random_seed=2)
        self.verify(ts)


class TestMutateRoundTrip(TestRoundTrip):
    def test_mutate_round_trip(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(
            ts, rate=2, random_seed=1, start_time=0, end_time=100, keep=False
        )
        self.verify(ts)

    def test_mutate_rate_map(self):
        ts = msprime.simulate(5, random_seed=1)
        rate_map = msprime.RateMap(position=[0, 0.5, 1], rate=[0, 1])
        ts = msprime.mutate(ts, rate=rate_map)
        self.verify(ts)


class TestRetainsProvenance(unittest.TestCase):
    def test_simulate_retains_provenance(self):
        from_ts = msprime.simulate(10, random_seed=5, end_time=0.5)
        ts = msprime.simulate(from_ts=from_ts, random_seed=2)
        self.assertEqual(ts.num_provenances, from_ts.num_provenances + 1)
        self.assertEqual(ts.provenance(0), from_ts.provenance(0))
