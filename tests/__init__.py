#
# Copyright (C) 2015-2018 University of Oxford
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
Common code for the msprime test cases.
"""
import unittest

import numpy as np


class SequenceEqualityMixin:
    """
    Overwrites unittest.TestCase.assertEqual to work with numpy arrays.

    Note: unittest.TestCase.assertSequenceEqual also fails to work with
    numpy arrays, and assertEqual works with ordinary lists/tuples anyway.
    """

    def assertEqual(self, it1, it2, msg=None):
        if isinstance(it1, np.ndarray):
            it1 = list(it1)
        if isinstance(it2, np.ndarray):
            it2 = list(it2)
        unittest.TestCase.assertEqual(self, it1, it2, msg=msg)
