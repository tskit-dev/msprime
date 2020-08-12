#
# Copyright (C) 2020 University of Oxford
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
Test cases for the high level interface to msprime.
"""
import os
import tempfile
import unittest

import msprime


class TestRecombinationMap(unittest.TestCase):
    """
    Tests for the RecombinationMap class.
    """

    def test_zero_recombination_map(self):
        # test that beginning and trailing zero recombination regions in the
        # recomb map are included in the sequence
        for n in range(3, 10):
            positions = list(range(n))
            rates = [0.0, 0.2] + [0.0] * (n - 2)
            recomb_map = msprime.RecombinationMap(positions, rates)
            ts = msprime.simulate(10, recombination_map=recomb_map)
            self.assertEqual(ts.sequence_length, n - 1)
            self.assertEqual(min(ts.tables.edges.left), 0.0)
            self.assertEqual(max(ts.tables.edges.right), n - 1.0)

    def test_mean_recombination_rate(self):
        # Some quick sanity checks.
        recomb_map = msprime.RecombinationMap([0, 1], [1, 0])
        mean_rr = recomb_map.mean_recombination_rate
        self.assertEqual(mean_rr, 1.0)

        recomb_map = msprime.RecombinationMap([0, 1, 2], [1, 0, 0])
        mean_rr = recomb_map.mean_recombination_rate
        self.assertEqual(mean_rr, 0.5)

        recomb_map = msprime.RecombinationMap([0, 1, 2], [0, 0, 0])
        mean_rr = recomb_map.mean_recombination_rate
        self.assertEqual(mean_rr, 0.0)

        # Test mean_recombination_rate is correct after reading from
        # a hapmap file. RecombinationMap.read_hapmap() ignores the cM
        # field, so here we test against using the cM field directly.
        def hapmap_rr(hapmap_file):
            first_pos = 0
            with open(hapmap_file) as f:
                next(f)  # skip header
                for line in f:
                    pos, rate, cM = map(float, line.split()[1:4])
                    if cM == 0:
                        first_pos = pos
            return cM / 100 / (pos - first_pos)

        hapmap = """chr pos        rate                    cM
                    1   4283592    3.79115663174456        0
                    1   4361401    0.0664276817058413      0.294986106359414
                    1   7979763   10.9082897515584         0.535345505591925
                    1   8007051    0.0976780648822495      0.833010916332456
                    1   8762788    0.0899929572085616      0.906829844052373
                    1   9477943    0.0864382908650907      0.971188757364862
                    1   9696341    4.76495005895746        0.990066707213216
                    1   9752154    0.0864316558730679      1.25601286485381
                    1   9881751    0.0                     1.26721414815999"""
        with tempfile.TemporaryDirectory() as temp_dir:
            hapfile = os.path.join(temp_dir, "hapmap.txt")
            with open(hapfile, "w") as f:
                f.write(hapmap)
            recomb_map = msprime.RecombinationMap.read_hapmap(f.name)
            mean_rr = recomb_map.mean_recombination_rate
            mean_rr2 = hapmap_rr(hapfile)
        self.assertAlmostEqual(mean_rr, mean_rr2, places=15)
