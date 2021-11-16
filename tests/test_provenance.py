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
import json
import logging

import numpy as np
import pytest
import python_jsonschema_objects as pjs
import tskit

import msprime
from msprime import _msprime


class TestProvenance:
    """
    Basic tests for the provenance dict function.
    """

    def test_libraries(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = json.loads(ts.provenance(0).record)
        libs = prov["environment"]["libraries"]
        assert libs["gsl"] == {
            "version": ".".join(map(str, _msprime.get_gsl_version()))
        }
        assert libs["tskit"] == {"version": tskit.__version__}


class ValidateSchemas:
    """
    Check that the schemas we produce in msprime are valid.
    """

    def test_simulation(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = json.loads(ts.provenance(0).record)
        tskit.validate_provenance(prov)
        assert prov["parameters"]["command"] == "simulate"

    def test_mutate(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(ts, rate=1, random_seed=1)
        prov = json.loads(ts.provenance(1).record)
        tskit.validate_provenance(prov)
        assert prov["parameters"]["command"] == "mutate"

    def test_simplify(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = ts.simplify()
        prov = json.loads(ts.provenance(1).record)
        tskit.validate_provenance(prov)
        assert prov["parameters"]["command"] == "simplify"


class TestBuildObjects:
    """
    Check that we can build objects from the json schema as we'd expect.
    """

    def decode(self, prov):
        builder = pjs.ObjectBuilder(tskit.provenance.get_schema())
        ns = builder.build_classes()
        return ns.TskitProvenance.from_json(prov)

    def test_simulate(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = ts.provenance(0).record
        decoded = self.decode(prov)
        assert decoded.schema_version == "1.0.0"
        assert decoded.parameters.command == "simulate"
        assert decoded.parameters.random_seed == 1

    def test_sim_ancestry(self):
        ts = msprime.sim_ancestry(5, random_seed=1)
        prov = ts.provenance(0).record
        decoded = self.decode(prov)
        assert decoded.schema_version == "1.0.0"
        assert decoded.parameters.command == "sim_ancestry"
        assert decoded.parameters.random_seed == 1

    def test_no_provenance(self):
        ts = msprime.simulate(5, record_provenance=False)
        assert len(list(ts.provenances())) == 0
        ts = msprime.sim_ancestry(5, record_provenance=False)
        assert len(list(ts.provenances())) == 0

    def test_encode_simulation_models(self):
        models = [
            msprime.StandardCoalescent(duration=10),
            msprime.DiscreteTimeWrightFisher(duration=10),
            msprime.SmcApproxCoalescent(duration=10),
            msprime.StandardCoalescent(),
        ]
        ts = msprime.sim_ancestry(10, model=models, random_seed=1234)
        decoded = self.decode(ts.provenance(0).record)
        parameters = decoded.parameters
        assert parameters.model[0] == {
            "__class__": "msprime.ancestry.StandardCoalescent",
            "duration": 10,
        }
        assert parameters.model[1] == {
            "__class__": "msprime.ancestry.DiscreteTimeWrightFisher",
            "duration": 10,
        }
        assert parameters.model[2] == {
            "__class__": "msprime.ancestry.SmcApproxCoalescent",
            "duration": 10,
        }
        assert parameters.model[3] == {
            "__class__": "msprime.ancestry.StandardCoalescent",
            "duration": None,
        }

    def test_unencodable_function(self):
        class BadClass:
            def __call__(self, t):
                return t

        with pytest.raises(TypeError):
            msprime.simulate(
                5,
                random_seed=1,
                demographic_events=[msprime.SimulationModelChange(time=BadClass)],
            )

    def test_replicate_index(self):
        for i, ts in enumerate(msprime.simulate(5, num_replicates=3)):
            decoded = self.decode(ts.provenance(0).record)
            assert decoded.parameters["replicate_index"] == i

        for i, ts in enumerate(msprime.sim_ancestry(5, num_replicates=3)):
            decoded = self.decode(ts.provenance(0).record)
            assert decoded.parameters["replicate_index"] == i

    def test_large_provenance_warning(self, caplog):
        with caplog.at_level(logging.WARNING):
            msprime.provenance.json_encode_provenance(
                {"test": np.zeros(shape=(500_000,))}, 5
            )
        assert (
            "The provenance information for the resulting"
            " tree sequence is 2.38MB. This is nothing to worry about as"
            " provenance is a good thing to have, but if you want to save this"
            " memory/storage space you can disable provenance recording by setting"
            " record_provenance=False" in caplog.text
        )

    def test_sim_mutations(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.sim_mutations(
            ts, rate=2, random_seed=1, start_time=0, end_time=100, keep=False
        )
        decoded = self.decode(ts.provenance(1).record)
        assert decoded.schema_version == "1.0.0"
        assert decoded.parameters.command == "sim_mutations"
        assert decoded.parameters.random_seed == 1
        assert decoded.parameters.rate == 2
        assert decoded.parameters.start_time == 0
        assert decoded.parameters.end_time == 100
        assert not decoded.parameters.keep

    def test_mutate_model(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.sim_mutations(ts, model="pam")
        decoded = self.decode(ts.provenance(1).record)
        assert decoded.schema_version == "1.0.0"
        assert decoded.parameters.command == "sim_mutations"
        assert decoded.parameters.model == "pam"

    def test_mutate_map(self):
        ts = msprime.simulate(5, random_seed=1)
        rate_map = msprime.RateMap(position=[0, 0.5, 1], rate=[0, 1])
        ts = msprime.sim_mutations(ts, rate=rate_map)
        decoded = self.decode(ts.provenance(1).record)
        assert decoded.schema_version == "1.0.0"
        assert decoded.parameters.command == "sim_mutations"
        assert decoded.parameters.rate["__class__"] == "msprime.intervals.RateMap"
        assert decoded.parameters.rate["position"]["__ndarray__"] == list(
            rate_map.position
        )
        assert decoded.parameters.rate["rate"]["__ndarray__"] == list(rate_map.rate)

    def test_mutate_numpy(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.sim_mutations(
            ts,
            rate=np.array([2])[0],
            random_seed=np.array([1])[0],
            start_time=np.array([0])[0],
            end_time=np.array([100][0]),
        )
        decoded = self.decode(ts.provenance(1).record)
        assert decoded.schema_version == "1.0.0"
        assert decoded.parameters.command == "sim_mutations"
        assert decoded.parameters.random_seed == 1
        # The dtype values change depending on platform, so not much
        # point in trying to test exactly.
        assert decoded.parameters.rate["__npgeneric__"] == "2"
        assert decoded.parameters.start_time["__npgeneric__"] == "0"
        assert decoded.parameters.end_time["__ndarray__"] == 100


class TestParseProvenance:
    def test_parse_non_msprime(self):
        # Most testing of parsing is in the round tripping below
        with pytest.raises(ValueError):
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
        ts1 = msprime.sim_ancestry(5, random_seed=1)
        ts2 = msprime.sim_mutations(ts1)
        command, prov = msprime.provenance.parse_provenance(ts2.provenance(1), ts1)
        assert command == "sim_mutations"
        assert prov["tree_sequence"] == ts1


class TestRoundTrip:
    """
    Check that we can recreate tree sequences from their provenances
    """

    def verify_equal(self, ts1, ts2):
        ts1_tables = ts1.dump_tables()
        ts2_tables = ts2.dump_tables()
        # We clear the provenances as they are not equal due to timestamps
        ts1_tables.provenances.clear()
        ts2_tables.provenances.clear()
        assert ts1_tables == ts2_tables

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
            demographic_events=[msprime.SimulationModelChange(time=10)],
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

    def test_mutate_numpy(self):
        ts = msprime.sim_ancestry(5, random_seed=1)
        ts = msprime.sim_mutations(
            ts,
            rate=np.array([2])[0],
            random_seed=np.array([1])[0],
            start_time=np.array([0])[0],
            end_time=np.array([100][0]),
        )
        self.verify(ts)


class TestRetainsProvenance:
    def test_simulate_retains_provenance(self):
        from_ts = msprime.simulate(10, random_seed=5, end_time=0.5)
        ts = msprime.simulate(from_ts=from_ts, random_seed=2)
        assert ts.num_provenances == from_ts.num_provenances + 1
        assert ts.provenance(0) == from_ts.provenance(0)
