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
from nose.tools import raises

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

    @raises(RuntimeError)
    def test_no_samples(self):
        ts = msprime.simulate(10, random_seed=15)
        msprime.continue_simulation(ts, 10000, sample_size=10000, random_seed=15)

    @raises(RuntimeError)
    def test_too_many_samples(self):
        ts = msprime.simulate(sample_size=2, random_seed=23)
        msprime.continue_simulation(ts, 1, sample_size=100000, random_seed=23)

    """def test_bug_instance1(self):
        ts = msprime.simulate(10, recombination_rate=10, end_time=0.5, random_seed=103)
        print("test_bug_instance1")
        self.verify(ts)

    def test_bug_instance2(self):
        ts = msprime.simulate(10, recombination_rate=10, end_time=1.0, random_seed=61)
        print("test_bug_instance2")
        self.verify(ts)

    def test_large_recombination(self):
        ts = msprime.simulate(15, recombination_rate=1.0, random_seed=2, end_time=0.25)
        print("test_large_recombination")
        self.verify(ts)

    def test_simple_recombination(self):
        ts = msprime.simulate(10, recombination_rate=0.1, random_seed=1, end_time=0.5)
        print("test_simple_recombination")
        self.verify(ts)

    def test_discrete_loci(self):
        ts = msprime.sim_ancestry(
            10,
            sequence_length=10,
            recombination_rate=1,
            random_seed=1,
            end_time=0.5,
        )
        print("test_discrete_loci")
        self.verify(ts)

    def test_simple_recombination_time_zero(self):
        ts = msprime.simulate(10, recombination_rate=0.1, random_seed=4, end_time=0.0)
        print("test_simple_recombination_time_zero")
        self.verify(ts)

    def test_no_recombination(self):
        ts = msprime.simulate(10, random_seed=2, end_time=0.5)
        print("test_no_recombination")
        self.verify(ts)

    def test_no_recombination_time_zero(self):
        ts = msprime.simulate(10, random_seed=3, end_time=0.0)
        print("test_no_recombination_time_zero")
        self.verify(ts)

    def test_dtwf_recombination(self):
        ts = msprime.simulate(
            10, Ne=100, model="dtwf", random_seed=2, end_time=100, recombination_rate=10
        )
        assert ts.num_trees > 1
        print("test_dtwf_recombination")
        self.verify(ts)

    def test_dtwf_no_recombination(self):
        ts = msprime.simulate(10, Ne=100, model="dtwf", random_seed=2, end_time=100)
        print("test_dtwf_no_recombination")
        self.verify(ts)

    def test_dtwf_no_recombination_time_zero(self):
        ts = msprime.simulate(10, Ne=100, model="dtwf", random_seed=2, end_time=0)
        print("test_dtwf_no_recombination_time_zero")
        self.verify(ts)

    def test_mixed_models_no_recombination(self):
        ts = msprime.simulate(
            10,
            Ne=10,
            model="dtwf",
            random_seed=2,
            end_time=100,
            demographic_events=[
                msprime.SimulationModelChange(10, msprime.StandardCoalescent())
            ],
        )
        print("test_mixed_models_no_recombination")
        self.verify(ts)"""
