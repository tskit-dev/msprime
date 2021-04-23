#
# Copyright (C) 2021 University of Oxford
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
Test cases for demes support.
"""
import demes

import msprime


class TestDemes:
    def test_ooa_example(self):
        b = demes.Builder(
            description="Gutenkunst et al. (2009) three-population model.",
            doi=["10.1371/journal.pgen.1000695"],
            time_units="years",
            generation_time=25,
        )
        b.add_deme("ANC", epochs=[dict(end_time=220e3, start_size=7300)])
        b.add_deme(
            "AMH",
            ancestors=["ANC"],
            epochs=[dict(end_time=140e3, start_size=12300)],
        )
        b.add_deme(
            "OOA", ancestors=["AMH"], epochs=[dict(end_time=21.2e3, start_size=2100)]
        )
        b.add_deme("YRI", ancestors=["AMH"], epochs=[dict(start_size=12300)])
        # Use floating point values here to make the comparisons simpler.
        b.add_deme(
            "CEU",
            ancestors=["OOA"],
            epochs=[dict(start_size=1000, end_size=29725.34354)],
        )
        b.add_deme(
            "CHB",
            ancestors=["OOA"],
            epochs=[dict(start_size=510, end_size=54090.33108)],
        )
        b.add_migration(demes=["YRI", "OOA"], rate=25e-5)
        b.add_migration(demes=["YRI", "CEU"], rate=3e-5)
        b.add_migration(demes=["YRI", "CHB"], rate=1.9e-5)
        b.add_migration(demes=["CEU", "CHB"], rate=9.6e-5)

        g = b.resolve()

        ooa1 = msprime.Demography.from_demes(g)
        ooa2 = msprime.Demography._ooa_model().copy([d.name for d in g.demes])
        ooa2.assert_equivalent(ooa1)
