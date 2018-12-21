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

import tskit

import _msprime
import msprime
import msprime.provenance as provenance


class TestProvenance(unittest.TestCase):
    """
    Basic tests for the provenance dict function.
    """

    @unittest.skip("FIX gsl version")
    def test_libraries(self):
        d = provenance.get_provenance_dict()
        libs = d["environment"]["libraries"]
        self.assertEqual(libs["gsl"], {
            "version": ".".join(map(str, _msprime.get_gsl_version()))})
        self.assertEqual(libs["kastore"], {"version": "0.1.0"})

    # TODO more tests when we finalise the format of these dictionaries.


@unittest.skip("TMP")
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
