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
import textwrap

import demes
import hypothesis as hyp
import numpy as np
import pytest

import msprime
from .demes_delete_me import graphs


def validate_demes_demography(graph, demography):
    """
    Checks that the specified Demes graph and msprime are consistent.
    """
    assert len(graph.demes) == len(demography.populations)
    assert {deme.name for deme in graph.demes} == set(demography.keys())
    for deme in graph.demes:
        population = demography[deme.name]
        # The default sampling time should always be the demes end_time
        assert population.default_sampling_time == deme.end_time
        # TODO once we have manual management of population states in #1679
        # and changed to use this in the demes converter we can turn back
        # on this check.
        # start_time = deme.end_time
        # end_time = deme.start_time
        # for epoch in dbg.epochs:
        #     pop = epoch.populations[pop_id]
        #     if end_time <= epoch.start_time:
        #         assert pop.state == PopulationStateMachine.PREVIOUSLY_ACTIVE
        #     elif start_time >= epoch.end_time:
        #         assert pop.state == PopulationStateMachine.INACTIVE
        #     else:
        #         assert pop.state == PopulationStateMachine.ACTIVE


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

    @hyp.settings(
        deadline=None, print_blob=True, suppress_health_check=[hyp.HealthCheck.too_slow]
    )
    @hyp.given(graphs())
    def test_random_graph(self, graph):
        demography = msprime.Demography.from_demes(graph)
        validate_demes_demography(graph, demography)


class TestYamlExamples:
    def from_yaml(self, yaml):
        graph = demes.loads(textwrap.dedent(yaml))
        demography = msprime.Demography.from_demes(graph)
        validate_demes_demography(graph, demography)
        return demography

    def get_active_populations(self, dbg):
        return {
            (e.start_time, e.end_time): [p.name for p in e.active_populations]
            for e in dbg.epochs
        }

    def test_one_pop_two_epoch(self):
        yaml = """\
        description: Single-population two-epoch demography.
        time_units: generations
        demes:

        - name: deme0
          description: A deme that doubles in size 100 generations ago.
          epochs:
          - start_size: 1000
            end_time: 100
          - start_size: 2000
            end_time: 0
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 1
        assert d.populations[0].name == "deme0"
        dbg = d.debug()
        assert dbg.num_epochs == 2
        assert dbg.epochs[0].populations[0].start_size == 2000
        assert dbg.epochs[1].populations[0].start_size == 1000
        assert dbg.epochs[1].start_time == 100
        active_pops = self.get_active_populations(dbg)
        assert "deme0" in active_pops[(0, 100)]
        assert "deme0" in active_pops[(100, np.inf)]

    def test_zigzag(self):
        yaml = """\
        time_units: generations
        demes:
          - name: generic
            epochs:
            - {end_time: 34133.31, start_size: 7156}
            - {end_time: 8533.33, end_size: 71560}
            - {end_time: 2133.33, end_size: 7156}
            - {end_time: 533.33, end_size: 71560}
            - {end_time: 133.33, end_size: 7156}
            - {end_time: 33.333, end_size: 71560}
            - {end_time: 0, end_size: 71560}
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 1
        dbg = d.debug()
        assert dbg.epochs[0].populations[0].start_size == pytest.approx(71560)
        assert dbg.epochs[1].populations[0].start_size == pytest.approx(71560)
        assert dbg.epochs[2].populations[0].start_size == pytest.approx(7156)
        assert dbg.epochs[3].populations[0].start_size == pytest.approx(71560)
        assert dbg.epochs[4].populations[0].start_size == pytest.approx(7156)
        assert dbg.epochs[5].populations[0].start_size == pytest.approx(71560)
        assert dbg.epochs[6].populations[0].start_size == pytest.approx(7156)
        assert dbg.epochs[0].start_time == 0
        assert dbg.epochs[1].start_time == pytest.approx(33.333)
        assert dbg.epochs[2].start_time == pytest.approx(133.33)
        assert dbg.epochs[3].start_time == pytest.approx(533.33)
        assert dbg.epochs[4].start_time == pytest.approx(2133.33)
        assert dbg.epochs[5].start_time == pytest.approx(8533.33)
        assert dbg.epochs[6].start_time == pytest.approx(34133.31)

    def test_split(self):
        yaml = """\
        time_units: generations
        demes:
          - name: X
            epochs:
              - end_time: 1000
                start_size: 2000
          - name: A
            ancestors:
              - X
            epochs:
              - start_size: 2000
          - name: B
            ancestors:
              - X
            epochs:
              - start_size: 2000
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 3
        assert d.populations[0].name == "X"
        assert d.populations[1].name == "A"
        assert d.populations[2].name == "B"
        assert d.populations[0].initial_size == 0
        assert d.populations[1].initial_size == 2000
        assert d.populations[2].initial_size == 2000
        assert d.num_events == 2
        assert isinstance(d.events[1], msprime.demography.PopulationSplit)

        dbg = d.debug()
        assert dbg.num_epochs == 2
        epoch = dbg.epochs[0]
        assert not epoch.populations[0].active
        assert epoch.populations[1].active
        assert epoch.populations[1].start_size == 2000
        assert epoch.populations[2].active
        assert epoch.populations[2].start_size == 2000

        epoch = dbg.epochs[1]
        assert epoch.start_time == 1000
        assert epoch.populations[0].active
        assert epoch.populations[0].start_size == 2000
        assert not epoch.populations[1].active
        assert not epoch.populations[2].active

        x = dbg.possible_lineage_locations(["A", "B"])
        assert len(x) == 2
        assert list(x[(0, 1000)]) == [False, True, True]
        assert list(x[(1000, np.inf)]) == [True, False, False]

        active_pops = self.get_active_populations(dbg)
        assert "A" in active_pops[(0, 1000)]
        assert "B" in active_pops[(0, 1000)]
        assert "X" in active_pops[(1000, np.inf)]
        assert "X" not in active_pops[(0, 1000)]
        assert "A" not in active_pops[(1000, np.inf)]
        assert "B" not in active_pops[(1000, np.inf)]

    def test_branch(self):
        yaml = """\
        time_units: generations
        demes:
          - name: X
            epochs:
              - start_size: 2000
          - name: A
            ancestors: [X]
            start_time: 1000
            epochs:
              - start_size: 2000
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 2

        dbg = d.debug()
        assert dbg.num_epochs == 2
        for epoch in dbg.epochs:
            assert epoch.populations[0].start_size == 2000
        assert dbg.epochs[0].populations[1].start_size == 2000
        assert dbg.epochs[1].populations[1].start_size == 0

        x = dbg.possible_lineage_locations(["X", "A"])
        assert len(x) == 2
        assert list(x[(0, 1000)]) == [True, True]
        assert list(x[(1000, np.inf)]) == [True, False]

        active_pops = self.get_active_populations(dbg)
        assert "A" in active_pops[(0, 1000)]
        assert "X" in active_pops[(0, 1000)]
        assert "X" in active_pops[(1000, np.inf)]
        assert "A" not in active_pops[(1000, np.inf)]

    def test_pulses(self):
        yaml = """\
        time_units: generations
        demes:
          - name: X
            epochs:
              - start_size: 2000
          - name: A
            ancestors: [X]
            start_time: 1000
            epochs:
              - start_size: 2000
        pulses:
          - source: X
            dest: A
            time: 500
            proportion: 0.1
          - source: A
            dest: X
            time: 100
            proportion: 0.2
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 2
        assert len(d.events) == 3
        dbg = d.debug()
        assert dbg.num_epochs == 4
        assert len(dbg.epochs[1].events) == 1
        assert len(dbg.epochs[2].events) == 1
        assert isinstance(dbg.epochs[1].events[0], msprime.demography.MassMigration)
        assert isinstance(dbg.epochs[2].events[0], msprime.demography.MassMigration)
        assert dbg.epochs[1].events[0].time == 100
        assert dbg.epochs[1].events[0].dest == "A"
        assert dbg.epochs[1].events[0].source == "X"
        assert dbg.epochs[1].events[0].proportion == 0.2
        assert dbg.epochs[2].events[0].time == 500
        assert dbg.epochs[2].events[0].dest == "X"
        assert dbg.epochs[2].events[0].source == "A"
        assert dbg.epochs[2].events[0].proportion == 0.1

        active_pops = self.get_active_populations(dbg)
        assert "A" in active_pops[(0, 100)]
        assert "X" in active_pops[(0, 100)]
        assert "A" in active_pops[(100, 500)]
        assert "X" in active_pops[(100, 500)]
        assert "A" in active_pops[(500, 1000)]
        assert "X" in active_pops[(500, 1000)]
        assert "X" in active_pops[(1000, np.inf)]
        assert "A" not in active_pops[(1000, np.inf)]

    def test_merger(self):
        yaml = """\
        time_units: generations
        demes:
          - name: X
            epochs:
              - start_size: 2000
                end_time: 1000
          - name: A
            ancestors: [X]
            epochs:
              - start_size: 2000
                end_time: 100
          - name: B
            ancestors: [X]
            epochs:
              - start_size: 2000
                end_time: 100
          - name: C
            ancestors: [A, B]
            proportions: [0.2, 0.8]
            start_time: 100
            epochs:
              - start_size: 2000
                end_time: 0
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 4
        assert len(d.events) == 5
        dbg = d.debug()
        assert (
            msprime.demography.Admixture(
                time=100, derived="C", ancestral=["A", "B"], proportions=[0.2, 0.8]
            )
            in d.events
        )
        x = dbg.possible_lineage_locations(["X", "A", "B", "C"])
        assert len(x) == 3
        assert list(x[(0, 100)]) == [False, False, False, True]
        assert list(x[(100, 1000)]) == [False, True, True, False]
        assert list(x[(1000, np.inf)]) == [True, False, False, False]

        active_pops = self.get_active_populations(dbg)
        assert "C" in active_pops[(0, 100)]
        assert "A" in active_pops[(100, 1000)]
        assert "B" in active_pops[(100, 1000)]
        assert "X" in active_pops[(1000, np.inf)]
        assert "C" not in active_pops[(100, 1000)]
        assert "C" not in active_pops[(1000, np.inf)]
        assert "A" not in active_pops[(0, 100)]
        assert "B" not in active_pops[(0, 100)]

    def test_admixture(self):
        yaml = """\
        time_units: generations
        demes:
          - name: X
            epochs:
              - start_size: 2000
                end_time: 1000
          - name: A
            ancestors: [X]
            epochs:
              - start_size: 2000
                end_time: 0
          - name: B
            ancestors: [X]
            epochs:
              - start_size: 2000
                end_time: 0
          - name: C
            ancestors: [A, B]
            proportions: [0.2, 0.8]
            start_time: 100
            epochs:
              - start_size: 2000
                end_time: 0
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 4
        assert len(d.events) == 3
        dbg = d.debug()
        assert (
            msprime.demography.Admixture(
                time=100, derived="C", ancestral=["A", "B"], proportions=[0.2, 0.8]
            )
            in d.events
        )
        x = dbg.possible_lineage_locations(["X", "A", "B", "C"])
        assert len(x) == 3
        assert list(x[(0, 100)]) == [False, True, True, True]
        assert list(x[(100, 1000)]) == [False, True, True, False]
        assert list(x[(1000, np.inf)]) == [True, False, False, False]

    def test_non_zero_end_time(self):
        yaml = """\
        time_units: generations
        demes:
          - name: A
            epochs:
              - start_size: 2000
                end_time: 100
        """
        self.from_yaml(yaml)
