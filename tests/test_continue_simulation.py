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
Test cases for continue_simulation.
"""
import pytest

import msprime


class TestContinueSimulation:
    def get_oldest_time(self, ts):
        return max(
            [
                ts.node(ts.first().roots[i]).time
                for i in range(0, len(ts.first().roots), 1)
            ]
        )

    def test_basic_function(self):
        ts = msprime.simulate(1000, random_seed=37)
        time = 1
        cts = msprime.continue_simulation(ts, time, random_seed=37)
        assert ts.num_mutations == cts.num_mutations
        assert ts.num_sites == cts.num_sites
        assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
        assert round(self.get_oldest_time(cts) - time, 10) == round(
            self.get_oldest_time(ts), 10
        )

    def test_with_recombination(self):
        ts = msprime.simulate(1000, recombination_rate=1e-8, random_seed=72)
        time = 1
        cts = msprime.continue_simulation(
            ts, time, recombination_rate=1e-8, random_seed=72
        )
        assert ts.num_mutations == cts.num_mutations
        assert ts.num_sites == cts.num_sites
        assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
        assert round(self.get_oldest_time(cts) - time, 10) == round(
            self.get_oldest_time(ts), 10
        )

    def test_large_simulation(self):
        ts = msprime.simulate(10000, Ne=5000, recombination_rate=1e-8, random_seed=72)
        time = 1
        cts = msprime.continue_simulation(
            ts,
            time,
            sample_size=10000,
            Ne=5000,
            recombination_rate=1e-8,
            random_seed=72,
        )
        assert ts.num_mutations == cts.num_mutations
        assert ts.num_sites == cts.num_sites
        assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
        assert round(self.get_oldest_time(cts) - time, 10) == round(
            self.get_oldest_time(ts), 10
        )

    def test_large_simulation_large_time(self):
        ts = msprime.simulate(10000, Ne=5000, recombination_rate=1e-8, random_seed=72)
        time = 1000
        cts = msprime.continue_simulation(
            ts, time, Ne=5000, recombination_rate=1e-8, random_seed=72
        )
        assert ts.num_mutations == cts.num_mutations
        assert ts.num_sites == cts.num_sites
        assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
        assert round(self.get_oldest_time(cts) - time, 10) == round(
            self.get_oldest_time(ts), 10
        )

    def test_large_simulation_large_time_different_sample_sizes(self):
        ts = msprime.simulate(10000, Ne=5000, recombination_rate=1e-8, random_seed=72)
        time = 1000
        cts = msprime.continue_simulation(
            ts, time, Ne=5000, sample_size=2000, recombination_rate=1e-8, random_seed=72
        )
        assert ts.num_mutations == cts.num_mutations
        assert ts.num_sites == cts.num_sites
        assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
        assert round(self.get_oldest_time(cts) - time, 10) == round(
            self.get_oldest_time(ts), 10
        )

    def test_large_simulation_large_time_with_mutation(self):
        ts = msprime.simulate(
            10000, Ne=5000, mutation_rate=1e-8, recombination_rate=1e-8, random_seed=72
        )
        time = 1000
        cts = msprime.continue_simulation(
            ts, time, Ne=5000, recombination_rate=1e-8, random_seed=72
        )
        assert ts.num_mutations == cts.num_mutations
        assert ts.num_sites == cts.num_sites
        assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
        assert round(self.get_oldest_time(cts) - time, 10) == round(
            self.get_oldest_time(ts), 10
        )

    def test_continue_nodes(self):
        ts = msprime.simulate(
            10000, Ne=5000, mutation_rate=1e-8, recombination_rate=1e-8, random_seed=72
        )
        continue_nodes = ts.samples()[1::2]
        time = 1000
        cts = msprime.continue_simulation(
            ts,
            time,
            continue_nodes=continue_nodes,
            Ne=5000,
            recombination_rate=1e-8,
            random_seed=72,
        )
        assert ts.num_mutations == cts.num_mutations
        assert ts.num_sites == cts.num_sites
        assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
        assert round(self.get_oldest_time(cts) - time, 10) == round(
            self.get_oldest_time(ts), 10
        )

    def test_no_samples(self):
        ts = msprime.simulate(10, random_seed=15)
        with pytest.raises(RuntimeError):
            msprime.continue_simulation(ts, 10000, sample_size=10000, random_seed=15)

    def test_too_many_samples(self):
        ts = msprime.simulate(sample_size=2, random_seed=23)
        with pytest.raises(RuntimeError):
            msprime.continue_simulation(ts, 1, sample_size=100000, random_seed=23)
