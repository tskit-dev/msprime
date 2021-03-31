# flake8: NOQA
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
Check that importing * from msprime works as expected.
"""
from msprime import *


def test_demography():
    demog = Demography.isolated_model([1, 2])
    assert demog.num_populations == 2


def test_sim_ancestry():
    ts = sim_ancestry(2, random_seed=1)
    assert ts.num_samples == 4


def test_sim_ancestry_model():
    ts = sim_ancestry(2, random_seed=1, model=SmcApproxCoalescent())
    assert ts.num_samples == 4


def test_sim_mutations():
    ts = sim_ancestry(2, random_seed=1)
    ts = sim_mutations(ts, rate=1, random_seed=1)
    assert ts.num_mutations > 0


def test_sim_mutation_model():
    ts = sim_ancestry(2, random_seed=1)
    ts = sim_mutations(ts, rate=1, model=JC69(), random_seed=1)
    assert ts.num_mutations > 0


def test_simulate():
    ts = simulate(2, random_seed=1)
    assert ts.num_samples == 2


def test_mutate():
    ts = simulate(2, random_seed=1)
    ts = mutate(ts, rate=1, random_seed=1)
    assert ts.num_mutations > 0


def test_dir():
    import msprime

    assert len(dir(msprime)) > 0
