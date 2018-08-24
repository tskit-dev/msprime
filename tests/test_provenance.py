#
# Copyright (C) 2018 University of Oxford
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
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import unittest
import json

import python_jsonschema_objects as pjs
import numpy as np

import _msprime
import msprime
import msprime.provenance as provenance


class TestProvenance(unittest.TestCase):
    """
    Basic tests for the provenance dict function.
    """

    def test_libraries(self):
        d = provenance.get_provenance_dict()
        libs = d["environment"]["libraries"]
        self.assertEqual(libs["gsl"], {
            "version": ".".join(map(str, _msprime.get_gsl_version()))})
        self.assertEqual(libs["kastore"], {"version": "0.1.0"})

    # TODO more tests when we finalise the format of these dictionaries.


class TestEnvironment(unittest.TestCase):
    """
    Basic tests for the provenance dict function.
    """
    def test_cache(self):
        d = provenance.get_environment()
        self.assertIn("os", d)
        self.assertIs(d, provenance.get_environment())


def get_provenance(
        software_name="x", software_version="y", schema_version="1", environment=None,
        parameters=None):
    document = {
        "schema_version": schema_version,
        "software": {
            "name": software_name,
            "version": software_version,
        },
        "environment": {} if environment is None else environment,
        "parameters": {} if parameters is None else parameters,
    }
    return document


class TestSchema(unittest.TestCase):
    """
    Tests for schema validation.
    """
    def test_empty(self):
        with self.assertRaises(msprime.ProvenanceValidationError):
            msprime.validate_provenance({})

    def test_missing_keys(self):
        minimal = get_provenance()
        msprime.validate_provenance(minimal)
        for key in minimal.keys():
            copy = dict(minimal)
            del copy[key]
            with self.assertRaises(msprime.ProvenanceValidationError):
                msprime.validate_provenance(copy)
        copy = dict(minimal)
        del copy["software"]["name"]
        with self.assertRaises(msprime.ProvenanceValidationError):
            msprime.validate_provenance(copy)
        copy = dict(minimal)
        del copy["software"]["version"]
        with self.assertRaises(msprime.ProvenanceValidationError):
            msprime.validate_provenance(copy)

    def test_software_types(self):
        for bad_type in [0, [1, 2, 3], {}]:
            doc = get_provenance(software_name=bad_type)
            with self.assertRaises(msprime.ProvenanceValidationError):
                msprime.validate_provenance(doc)
            doc = get_provenance(software_version=bad_type)
            with self.assertRaises(msprime.ProvenanceValidationError):
                msprime.validate_provenance(doc)

    def test_schema_version_empth(self):
        doc = get_provenance(schema_version="")
        with self.assertRaises(msprime.ProvenanceValidationError):
            msprime.validate_provenance(doc)

    def test_software_empty_strings(self):
        doc = get_provenance(software_name="")
        with self.assertRaises(msprime.ProvenanceValidationError):
            msprime.validate_provenance(doc)
        doc = get_provenance(software_version="")
        with self.assertRaises(msprime.ProvenanceValidationError):
            msprime.validate_provenance(doc)

    def test_minimal(self):
        minimal = {
            "schema_version": "1",
            "software": {
                "name": "x",
                "version": "y",
            },
            "environment": {},
            "parameters": {}
        }
        msprime.validate_provenance(minimal)

    def test_extra_stuff(self):
        extra = {
            "you": "can",
            "schema_version": "1",
            "software": {
                "put": "anything",
                "name": "x",
                "version": "y",
            },
            "environment": {"extra": ["you", "want"]},
            "parameters": {"so": ["long", "its", "JSON", 0]}
        }
        msprime.validate_provenance(extra)


class ValidateSchemas(unittest.TestCase):
    """
    Check that the schemas we produce in msprime are valid.
    """
    def test_simulation(self):
        ts = msprime.simulate(5, random_seed=1)
        prov = json.loads(ts.provenance(0).record)
        msprime.validate_provenance(prov)
        self.assertEqual(prov["parameters"]["command"], "simulate")

    def test_mutate(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = msprime.mutate(ts, rate=1, random_seed=1)
        prov = json.loads(ts.provenance(1).record)
        msprime.validate_provenance(prov)
        self.assertEqual(prov["parameters"]["command"], "mutate")

    def test_simplify(self):
        ts = msprime.simulate(5, random_seed=1)
        ts = ts.simplify()
        prov = json.loads(ts.provenance(1).record)
        msprime.validate_provenance(prov)
        self.assertEqual(prov["parameters"]["command"], "simplify")


class TestBuildObjects(unittest.TestCase):
    """
    Check that we can build objects from the json schema as we'd expect.
    """
    def decode(self, prov):
        schema_file = "msprime/provenance.schema.json"
        with open(schema_file) as f:
            schema = f.read()
        builder = pjs.ObjectBuilder(json.loads(schema))
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
