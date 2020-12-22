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
Test cases for the continue_simulation and _adjust_tables_time methods.
"""
# import itertools
# import numpy as np
# import pytest
# import tskit
import msprime

# import tests.tsutil as tsutil
# import tests.wright_fisher as wf
# from msprime import _msprime


class TestContinueSimulation:
    def verify(self, ts):
        if ts.num_populations == 1:
            # recomb_rate = 1.0 / ts.sequence_length
            times = [1, 1.0, 100, 100.0, 100.5, 1000, 1000.0, 1000.5]
            for time in times:
                print(len(ts.samples()))
                cts = msprime.continue_simulation(
                    ts,
                    time,
                    Ne=20000,
                    mutation_rate=1e-8,
                    recombination_rate=1e-8,
                    sample_size=10000,
                )
                assert ts.num_mutations == cts.num_mutations
                assert ts.num_sites == cts.num_sites
                assert list(ts.tables.sites.position) == list(cts.tables.sites.position)
                # need to get oldest root for each, indexes change
                max_old = max(
                    [
                        ts.node(ts.first().roots[i]).time
                        for i in range(0, len(ts.first().roots), 1)
                    ]
                )
                max_new = max(
                    [
                        cts.node(cts.first().roots[i]).time
                        for i in range(0, len(cts.first().roots), 1)
                    ]
                )
                # fr = cts.node(cts.first().roots[0])
                # r = ts.node(ts.first().roots[0])
                print(max_old, max_new)
                assert round(max_new - time, 10) == round(max_old, 10)

    def test_bug_instance1(self):
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
        self.verify(ts)
