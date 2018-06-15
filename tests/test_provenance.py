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

import _msprime

import msprime.provenance as provenance


class TestProvenance(unittest.TestCase):
    """
    Basic tests for the provenance dict function.
    """

    def test_libraries(self):
        d = provenance.get_provenance_dict("test")
        libs = d["environment"]["libraries"]
        self.assertEqual(libs["gsl"], {"version": _msprime.get_gsl_version()})
        self.assertEqual(libs["kastore"], {"version": (0, 1, 0)})

    # TODO more tests when we finalise the format of these dictionaries.
