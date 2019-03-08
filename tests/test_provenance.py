#
# Copyright (C) 2018-2019 University of Oxford
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
import unittest
import json

import tskit
import numpy as np
import python_jsonschema_objects as pjs

import _msprime
import msprime


class TestProvenance(unittest.TestCase):
    """
    Basic tests for the provenance dict function.
    """
    def test_libraries(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = json.loads(ts.provenance(0).record)
        libs = prov["environment"]["libraries"]
        self.assertEqual(libs["gsl"], {
            "version": ".".join(map(str, _msprime.get_gsl_version()))})
        self.assertEqual(libs["tskit"], {"version": tskit.__version__})


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

    def test_simulation_numpy(self):
        seeds = np.ones(1, dtype=int)
        ts = msprime.simulate(5, random_seed=seeds[0])
        prov = ts.provenance(0).record
        decoded = self.decode(prov)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "simulate")
        self.assertEqual(decoded.parameters.random_seed, 1)

    def test_mutate(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(
            ts, rate=2, random_seed=1, start_time=0, end_time=100, keep=False)
        decoded = self.decode(ts.provenance(1).record)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "mutate")
        self.assertEqual(decoded.parameters.random_seed, 1)
        self.assertEqual(decoded.parameters.rate, 2)
        self.assertEqual(decoded.parameters.start_time, 0)
        self.assertEqual(decoded.parameters.end_time, 100)
        self.assertEqual(decoded.parameters.keep, False)

    def test_mutate_numpy(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(
            ts,
            rate=np.array([2])[0],
            random_seed=np.array([1])[0],
            start_time=np.array([0])[0],
            end_time=np.array([100][0]),
            keep=np.array([False][0]))
        decoded = self.decode(ts.provenance(1).record)
        self.assertEqual(decoded.schema_version, "1.0.0")
        self.assertEqual(decoded.parameters.command, "mutate")
        self.assertEqual(decoded.parameters.random_seed, 1)
        self.assertEqual(decoded.parameters.rate, 2)
        self.assertEqual(decoded.parameters.start_time, 0)
        self.assertEqual(decoded.parameters.end_time, 100)
        self.assertEqual(decoded.parameters.keep, False)
