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
    def verify(self, ts, new_ts, time, continue_nodes=None):
        assert round(new_ts.max_root_time - time, 10) == round(ts.max_root_time, 10)
        assert ts.num_mutations == new_ts.num_mutations
        assert ts.num_sites == new_ts.num_sites
        assert list(ts.tables.sites.position) == list(new_ts.tables.sites.position)
        if continue_nodes is None:
            continue_nodes = ts.samples()
        section = list(set(continue_nodes).intersection(new_ts.samples()))
        assert len(section) == 0
        for tree in new_ts.trees():
            assert new_ts.node(tree.root).time > time
        # check that subsetting to all nodes present before time gets back orig. tables
        # sub_ts = new_ts.subset(np.where(new_ts.tables.nodes.time > time)[0])
        # assert sub_ts.tables.equals(ts.tables)

    def test_basic_function(self):
        ts = msprime.sim_ancestry(samples=1000, sequence_length=1.0, random_seed=37)
        time = 1
        cts = msprime.continue_simulation(ts, time, sequence_length=1.0, random_seed=37)
        self.verify(ts, cts, time)

    def test_with_recombination(self):
        ts = msprime.sim_ancestry(
            samples=1000, recombination_rate=1e-8, sequence_length=1.0, random_seed=72
        )
        time = 1
        cts = msprime.continue_simulation(
            ts, time, recombination_rate=1e-8, sequence_length=1.0, random_seed=72
        )
        self.verify(ts, cts, time)

    def test_large_simulation(self):
        ts = msprime.sim_ancestry(
            samples=10000,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        time = 1
        cts = msprime.continue_simulation(
            ts,
            time,
            sample_size=10000,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        self.verify(ts, cts, time)

    def test_large_simulation_large_time(self):
        ts = msprime.sim_ancestry(
            samples=10000,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        time = 1000
        cts = msprime.continue_simulation(
            ts,
            time,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        self.verify(ts, cts, time)

    def test_large_simulation_large_time_different_sample_sizes(self):
        ts = msprime.sim_ancestry(
            samples=10000,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        time = 1000
        cts = msprime.continue_simulation(
            ts,
            time,
            population_size=5000,
            sample_size=2000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        self.verify(ts, cts, time)

    def test_large_simulation_large_time_with_mutation(self):
        ts = msprime.sim_ancestry(
            samples=10000,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        ts = msprime.sim_mutations(ts, rate=1e-8)
        time = 1000
        cts = msprime.continue_simulation(
            ts,
            time,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        self.verify(ts, cts, time)

    def test_continue_nodes(self):
        ts = msprime.sim_ancestry(
            samples=10000,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        ts = msprime.sim_mutations(ts, rate=1e-8)
        continue_nodes = ts.samples()[0::2]
        time = 1000
        cts = msprime.continue_simulation(
            ts,
            time,
            continue_nodes=continue_nodes,
            population_size=5000,
            recombination_rate=1e-8,
            sequence_length=1e4,
            random_seed=72,
        )
        self.verify(ts, cts, time, continue_nodes=continue_nodes)

    def test_no_uncoalesced_samples(self):
        ts = msprime.sim_ancestry(samples=10, sequence_length=1.0, random_seed=15)
        time = 10000
        cts = msprime.continue_simulation(
            ts, time, sample_size=10000, sequence_length=1.0, random_seed=15
        )
        self.verify(ts, cts, time)

    def test_too_many_samples(self):
        ts = msprime.sim_ancestry(samples=2, sequence_length=1.0, random_seed=23)
        with pytest.raises(RuntimeError):
            msprime.continue_simulation(
                ts, 1, sample_size=100000, sequence_length=1.0, random_seed=23
            )
