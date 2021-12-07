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
# import copy
import math
import textwrap

import demes
import numpy as np
import pytest

import msprime


def validate_demes_demography(graph, demography):
    """
    Checks that the specified Demes graph and msprime are consistent.
    """
    graph = graph.in_generations()
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


class TestFromDemes:
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

    def test_selfing_rate_unsupported(self):
        b = demes.Builder()
        b.add_deme("A", epochs=[dict(start_size=1, selfing_rate=0.1)])
        graph = b.resolve()
        with pytest.raises(ValueError, match="selfing_rate"):
            msprime.Demography.from_demes(graph)

    def test_cloning_rate_unsupported(self):
        b = demes.Builder()
        b.add_deme("A", epochs=[dict(start_size=1, cloning_rate=0.1)])
        graph = b.resolve()
        with pytest.raises(ValueError, match="cloning_rate"):
            msprime.Demography.from_demes(graph)


class TestFromYamlExamples:
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
          - sources: [X]
            dest: A
            time: 500
            proportions: [0.1]
          - sources: [A]
            dest: X
            time: 100
            proportions: [0.2]
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

    @pytest.mark.filterwarnings(
        # demes warns about multiple pulses with the same time
        "ignore:Multiple pulses:UserWarning"
    )
    def test_pulses_ordering(self):
        yaml = """\
        time_units: generations
        defaults:
          epoch:
            start_size: 2000
        demes:
          - name: A
          - name: B
          - name: C
          - name: D
          - name: E
        pulses:
          - {sources: [A], dest: B, time: 100, proportions: [0.1]}
          - {sources: [B], dest: C, time: 100, proportions: [0.2]}
          - {sources: [C], dest: D, time: 100, proportions: [0.3]}
          - {sources: [A], dest: E, time: 200, proportions: [0.4]}
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 5
        assert len(d.events) == 4
        for event in d.events:
            assert isinstance(event, msprime.demography.MassMigration)
        assert d.events[0].time == 100
        assert d.events[0].source == "D"
        assert d.events[0].dest == "C"
        assert d.events[0].proportion == 0.3
        assert d.events[1].time == 100
        assert d.events[1].source == "C"
        assert d.events[1].dest == "B"
        assert d.events[1].proportion == 0.2
        assert d.events[2].time == 100
        assert d.events[2].source == "B"
        assert d.events[2].dest == "A"
        assert d.events[2].proportion == 0.1
        assert d.events[3].time == 200
        assert d.events[3].source == "E"
        assert d.events[3].dest == "A"
        assert d.events[3].proportion == 0.4

    def test_multipulse(self):
        b = demes.Builder(defaults=dict(epoch=dict(start_size=1)))
        b.add_deme("a")
        b.add_deme("b")
        b.add_deme("c")
        b.add_pulse(sources=["a", "b"], proportions=[0.1, 0.2], dest="c", time=1)
        g = b.resolve()
        d = msprime.Demography.from_demes(g)
        assert len(d.events) == 2
        assert isinstance(d.events[0], msprime.MassMigration)
        assert isinstance(d.events[1], msprime.MassMigration)
        assert d.events[0].source == "c"
        assert d.events[0].dest == "a"
        assert d.events[0].proportion == 0.1
        assert d.events[1].source == "c"
        assert d.events[1].dest == "b"
        assert d.events[1].proportion == 0.2 / (1 - 0.1)

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

    def test_symmetric_migration(self):
        yaml = """\
        time_units: generations
        defaults:
          epoch:
            start_size: 2000
        demes:
          - name: A
          - name: B
          - name: C
        migrations:
          - {demes: [A, B, C], rate: 0.1}
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 3
        assert len(d.events) == 0
        for j in range(2):
            for k in range(2):
                rate = d.migration_matrix[j, k]
                if j != k:
                    assert rate == 0.1
                else:
                    assert rate == 0

    def check_asymmetric_matrix(self, migration_matrix, expected_rate, j, k):
        N = len(migration_matrix)
        for jj in range(N):
            for kk in range(N):
                rate = migration_matrix[jj][kk]
                if j == jj and k == kk:
                    assert rate == expected_rate
                else:
                    assert rate == 0

    def test_asymmetric_migration(self):
        yaml = """\
        time_units: generations
        defaults:
          epoch:
            start_size: 2000
        demes:
          - name: A
          - name: B
        migrations:
          - {source: A, dest: B, rate: 0.1}
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 2
        assert len(d.events) == 0
        self.check_asymmetric_matrix(d.migration_matrix, 0.1, d["B"].id, d["A"].id)

    def test_asymmetric_migration_multiple(self):
        yaml = """\
        time_units: generations
        defaults:
          epoch:
            start_size: 2000
        demes:
          - name: A
          - name: B
        migrations:
          - {source: A, dest: B, rate: 0.1, start_time: 200, end_time: 100}
          - {source: A, dest: B, rate: 0.2, start_time: 100, end_time: 50}
          - {source: A, dest: B, rate: 0.3, start_time: 50, end_time: 0}
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 2
        dbg = d.debug()
        assert dbg.num_epochs == 4

        j = d["B"].id
        k = d["A"].id
        self.check_asymmetric_matrix(dbg.epochs[0].migration_matrix, 0.3, j, k)
        self.check_asymmetric_matrix(dbg.epochs[1].migration_matrix, 0.2, j, k)
        self.check_asymmetric_matrix(dbg.epochs[2].migration_matrix, 0.1, j, k)
        self.check_asymmetric_matrix(dbg.epochs[3].migration_matrix, 0, j, k)

    def test_activate_ghost_split(self):
        yaml = """\
        time_units: generations
        demes:
          - name: A
            epochs: [{start_size: 100, end_time: 100}]
          - name: B
            ancestors: [A]
            epochs: [{start_size: 100}]
          - name: C
            ancestors: [A]
            epochs: [{start_size: 100, end_time: 50}]
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 3
        dbg = d.debug()
        assert dbg.num_epochs == 3

        assert np.all(dbg.epoch_start_time == np.array([0, 50, 100]))

        assert len(dbg.epochs[0].active_populations) == 1
        assert "B" == dbg.epochs[0].active_populations[0].name
        assert len(dbg.epochs[1].active_populations) == 2
        active_pops = [_.name for _ in dbg.epochs[1].active_populations]
        assert "B" in active_pops and "C" in active_pops
        assert len(dbg.epochs[2].active_populations) == 1
        assert "A" in dbg.epochs[2].active_populations[0].name

    def test_activate_ghost_admixture(self):
        yaml = """\
        description: two populations admix, with one persisting but extinct before zero
        time_units: generations
        demes:
          - name: A
            epochs: [{start_size: 100, end_time: 100}]
          - name: B
            ancestors: [A]
            epochs: [{start_size: 100, end_time: 25}]
          - name: C
            ancestors: [A]
            epochs: [{start_size: 100, end_time: 50}]
          - name: D
            ancestors: [B, C]
            proportions: [0.3, 0.7]
            start_time: 50
            epochs: [{start_size: 100, end_time: 0}]
        """
        d = self.from_yaml(yaml)
        assert d.num_populations == 4
        dbg = d.debug()
        assert dbg.num_epochs == 4

        assert np.all(dbg.epoch_start_time == np.array([0, 25, 50, 100]))

        assert len(dbg.epochs[0].active_populations) == 1
        assert "D" == dbg.epochs[0].active_populations[0].name
        assert len(dbg.epochs[1].active_populations) == 2
        active_pops = [_.name for _ in dbg.epochs[1].active_populations]
        assert "B" in active_pops and "D" in active_pops
        assert len(dbg.epochs[2].active_populations) == 2
        active_pops = [_.name for _ in dbg.epochs[2].active_populations]
        assert "B" in active_pops and "C" in active_pops
        assert len(dbg.epochs[3].active_populations) == 1
        assert "A" in dbg.epochs[3].active_populations[0].name


class TestToDemes:
    def test_ooa_example(self):
        ooa1 = msprime.Demography._ooa_model()
        graph = ooa1.to_demes()
        assert len(graph.demes) == len(ooa1)
        assert sorted(d.name for d in graph.demes) == sorted(ooa1)
        assert graph["CEU"].ancestors == ["OOA"]
        assert graph["CHB"].ancestors == ["OOA"]
        assert graph["OOA"].ancestors == ["AMH"]
        assert graph["YRI"].ancestors == ["AMH"]
        assert graph["AMH"].ancestors == ["ANC"]
        assert graph["ANC"].ancestors == []
        assert len(graph.migrations) == 8
        assert len(graph.pulses) == 0

    def test_american_admixture_example(self):
        aa1 = msprime.Demography._american_admixture_model()
        graph = aa1.to_demes()
        assert len(graph.demes) == len(aa1)
        assert sorted(d.name for d in graph.demes) == sorted(aa1)
        assert sorted(graph["ADMIX"].ancestors) == sorted(["EUR", "EAS", "AFR"])
        assert graph["EUR"].ancestors == ["OOA"]
        assert graph["EAS"].ancestors == ["OOA"]
        assert graph["AFR"].ancestors == ["AMH"]
        assert graph["OOA"].ancestors == ["AMH"]
        assert graph["AMH"].ancestors == ["ANC"]
        assert graph["ANC"].ancestors == []
        assert len(graph.migrations) == 8
        assert len(graph.pulses) == 0

    @pytest.mark.parametrize("num_pops", [1, 2, 3])
    def test_isolated(self, num_pops):
        demog = msprime.Demography.isolated_model([1000] * num_pops)
        graph = demog.to_demes()
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 0
        assert len(graph.demes) == num_pops
        for j in range(num_pops):
            name = f"pop_{j}"
            assert name in graph
            assert len(graph[name].epochs) == 1
            epoch = graph[name].epochs[0]
            assert math.isinf(epoch.start_time)
            assert epoch.end_time == 0
            assert epoch.start_size == 1000
            assert epoch.end_size == 1000

    @pytest.mark.parametrize("num_derived", [1, 2, 3])
    def test_split(self, num_derived):
        demog = msprime.Demography.isolated_model([1000] * (num_derived + 1))
        demog.add_population_split(
            100,
            derived=[f"pop_{j}" for j in range(1, num_derived + 1)],
            ancestral="pop_0",
        )
        graph = demog.to_demes()
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 0
        assert len(graph.demes) == num_derived + 1
        assert graph["pop_0"].ancestors == []
        assert graph["pop_0"].end_time == 100
        for j in range(1, num_derived + 1):
            name = f"pop_{j}"
            assert graph[name].ancestors == ["pop_0"]
            assert graph[name].start_time == 100
            assert len(graph[name].epochs) == 1
            epoch = graph[name].epochs[0]
            assert epoch.end_time == 0
            assert epoch.start_size == 1000
            assert epoch.end_size == 1000

    def test_split_via_mass_migration(self):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.add_mass_migration(100, source="pop_1", dest="pop_0", proportion=1)
        graph = demog.to_demes()
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 0
        assert len(graph.demes) == 2
        assert graph["pop_0"].ancestors == []
        assert graph["pop_0"].end_time == 0
        assert graph["pop_1"].ancestors == ["pop_0"]
        assert graph["pop_1"].start_time == 100

    @pytest.mark.parametrize("num_ancestral", [1, 2, 3])
    def test_admixture(self, num_ancestral):
        ancestral = [f"pop_{j}" for j in range(1, num_ancestral + 1)]
        proportions = [1 / num_ancestral for _ in range(num_ancestral)]
        demog = msprime.Demography.isolated_model([1000] * (num_ancestral + 1))
        demog.add_admixture(
            100,
            derived="pop_0",
            ancestral=ancestral,
            proportions=proportions,
        )
        graph = demog.to_demes()
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 0
        assert len(graph.demes) == num_ancestral + 1
        assert tuple(
            sorted(zip(graph["pop_0"].ancestors, graph["pop_0"].proportions))
        ) == tuple(zip(ancestral, proportions))
        assert graph["pop_0"].start_time == 100
        assert graph["pop_0"].end_time == 0
        for j in range(1, num_ancestral + 1):
            name = f"pop_{j}"
            assert graph[name].ancestors == []
            assert graph[name].end_time == 0
            assert len(graph[name].epochs) == 1
            epoch = graph[name].epochs[0]
            assert epoch.end_time == 0
            assert epoch.start_size == 1000
            assert epoch.end_size == 1000

    def test_admixture_via_mass_migration(self):
        demog = msprime.Demography.isolated_model([1000] * 4)
        demog.add_mass_migration(100, source="pop_0", dest="pop_1", proportion=0.5)
        demog.add_mass_migration(100, source="pop_0", dest="pop_2", proportion=0.5)
        demog.add_mass_migration(100, source="pop_0", dest="pop_3", proportion=1)
        graph = demog.to_demes()
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 0
        assert len(graph.demes) == 4
        assert sorted(graph["pop_0"].ancestors) == ["pop_1", "pop_2", "pop_3"]
        assert tuple(
            sorted(zip(graph["pop_0"].ancestors, graph["pop_0"].proportions))
        ) == tuple(zip(["pop_1", "pop_2", "pop_3"], [0.5, 0.25, 0.25]))
        assert graph["pop_0"].start_time == 100
        assert graph["pop_0"].end_time == 0
        for j in range(1, len(demog)):
            name = f"pop_{j}"
            assert graph[name].ancestors == []
            assert graph[name].end_time == 0
            assert len(graph[name].epochs) == 1
            epoch = graph[name].epochs[0]
            assert epoch.end_time == 0
            assert epoch.start_size == 1000
            assert epoch.end_size == 1000

    @pytest.mark.parametrize(
        "time_lo,time_hi", [(0, math.inf), (0, 20), (10, math.inf), (10, 20)]
    )
    def test_asymmetric_migration_via_rate_change(self, time_lo, time_hi):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.add_migration_rate_change(time_lo, rate=0.1, source="pop_0", dest="pop_1")
        if not math.isinf(time_hi):
            demog.add_migration_rate_change(time_hi, rate=0)
        graph = demog.to_demes()
        assert len(graph.demes) == 2
        assert len(graph.migrations) == 1
        assert len(graph.pulses) == 0
        migration = graph.migrations[0]
        assert migration.source == "pop_1"
        assert migration.dest == "pop_0"
        assert migration.rate == 0.1
        if math.isinf(time_hi):
            assert math.isinf(migration.start_time)
        else:
            assert migration.start_time == time_hi
        assert migration.end_time == time_lo

    def test_asymmetric_migration_via_initial_matrix(self):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.set_migration_rate(source="pop_0", dest="pop_1", rate=0.1)
        graph = demog.to_demes()
        assert len(graph.demes) == 2
        assert len(graph.migrations) == 1
        assert len(graph.pulses) == 0
        migration = graph.migrations[0]
        assert migration.source == "pop_1"
        assert migration.dest == "pop_0"
        assert migration.rate == 0.1
        assert math.isinf(migration.start_time)
        assert migration.end_time == 0

    def test_asymmetric_migration_multiple(self):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.add_migration_rate_change(10, rate=0.1, source="pop_0", dest="pop_1")
        demog.add_migration_rate_change(15, rate=0.1, source="pop_0", dest="pop_1")
        demog.add_migration_rate_change(20, rate=0.2, source="pop_0", dest="pop_1")
        demog.add_migration_rate_change(25, rate=0.2, source="pop_0", dest="pop_1")
        demog.add_migration_rate_change(30, rate=0, source="pop_0", dest="pop_1")
        graph = demog.to_demes()
        assert len(graph.demes) == 2
        assert len(graph.migrations) == 2
        assert len(graph.pulses) == 0
        assert graph.migrations[0].source == "pop_1"
        assert graph.migrations[0].dest == "pop_0"
        assert graph.migrations[0].start_time == 30
        assert graph.migrations[0].end_time == 20
        assert graph.migrations[0].rate == 0.2
        assert graph.migrations[1].source == "pop_1"
        assert graph.migrations[1].dest == "pop_0"
        assert graph.migrations[1].start_time == 20
        assert graph.migrations[1].end_time == 10
        assert graph.migrations[1].rate == 0.1

    @pytest.mark.parametrize("num_pops", [2, 5, 10])
    def test_symmetric_migration(self, num_pops):
        demog = msprime.Demography.island_model([1000] * num_pops, 0.1)
        graph = demog.to_demes()
        assert len(graph.demes) == num_pops
        assert len(graph.migrations) == num_pops * (num_pops - 1)
        assert len(graph.pulses) == 0
        for migration in graph.migrations:
            assert math.isinf(migration.start_time)
            assert migration.end_time == 0
            assert migration.rate == 0.1

    @pytest.mark.parametrize("num_pops", [2, 5, 10])
    def test_symmetric_migration_time_limited(self, num_pops):
        demog = msprime.Demography.isolated_model([1000] * num_pops)
        demog.add_symmetric_migration_rate_change(
            100, populations=list(demog), rate=0.1
        )
        demog.add_symmetric_migration_rate_change(200, populations=list(demog), rate=0)
        graph = demog.to_demes()
        assert len(graph.demes) == num_pops
        assert len(graph.migrations) == num_pops * (num_pops - 1)
        assert len(graph.pulses) == 0
        for migration in graph.migrations:
            assert migration.start_time == 200
            assert migration.end_time == 100
            assert migration.rate == 0.1

    @pytest.mark.filterwarnings("ignore:Non-zero migration.*after merging:UserWarning")
    def test_bad_migration1(self):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.add_mass_migration(100, source="pop_0", dest="pop_1", proportion=1)
        # Migration after pop_0 is merged.
        demog.add_migration_rate_change(200, source="pop_1", dest="pop_0", rate=0.1)
        with pytest.raises(ValueError, match="invalid migration"):
            demog.to_demes()

    @pytest.mark.filterwarnings("ignore:Non-zero migration.*after merging:UserWarning")
    def test_bad_migration2(self):
        demog = msprime.Demography.island_model([1000] * 2, 0.1)
        demog.add_mass_migration(100, source="pop_0", dest="pop_1", proportion=1)
        # Forgot to turn off migration after merging pop_0 into pop_1.
        with pytest.raises(ValueError, match="invalid migration"):
            demog.to_demes()

    def test_pulse(self):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.add_mass_migration(100, source="pop_0", dest="pop_1", proportion=0.1)
        graph = demog.to_demes()
        assert len(graph.demes) == 2
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 1
        pulse = graph.pulses[0]
        assert pulse.time == 100
        assert pulse.sources[0] == "pop_1"
        assert pulse.dest == "pop_0"
        assert pulse.proportions[0] == 0.1

    @pytest.mark.filterwarnings(
        # demes warns about multiple pulses with the same time
        "ignore:Multiple pulses:UserWarning"
    )
    def test_pulse_ordering(self):
        demog = msprime.Demography.isolated_model([1000] * 5)
        demog.add_mass_migration(100, source="pop_0", dest="pop_1", proportion=0.1)
        demog.add_mass_migration(100, source="pop_1", dest="pop_2", proportion=0.2)
        demog.add_mass_migration(100, source="pop_2", dest="pop_3", proportion=0.3)
        demog.add_mass_migration(200, source="pop_3", dest="pop_4", proportion=0.4)
        graph = demog.to_demes()
        assert len(graph.demes) == 5
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 4
        pulses = graph.pulses
        assert pulses[0].sources[0] == "pop_4"
        assert pulses[0].dest == "pop_3"
        assert pulses[0].proportions[0] == 0.4
        assert pulses[0].time == 200
        assert pulses[1].sources[0] == "pop_3"
        assert pulses[1].dest == "pop_2"
        assert pulses[1].proportions[0] == 0.3
        assert pulses[1].time == 100
        assert pulses[2].sources[0] == "pop_2"
        assert pulses[2].dest == "pop_1"
        assert pulses[2].proportions[0] == 0.2
        assert pulses[2].time == 100
        assert pulses[3].sources[0] == "pop_1"
        assert pulses[3].dest == "pop_0"
        assert pulses[3].proportions[0] == 0.1
        assert pulses[3].time == 100

    def test_bad_pulse(self):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.add_mass_migration(100, source="pop_0", dest="pop_1", proportion=1)
        demog.add_mass_migration(200, source="pop_0", dest="pop_1", proportion=0.1)
        with pytest.raises(ValueError, match="invalid pulse"):
            demog.to_demes()

    @pytest.mark.parametrize(
        "sizes,times",
        (([1000], [0]), ([1000, 500], [0, 20]), ([1000, 50, 1000], [0, 20, 30])),
    )
    def test_piecewise_constant_sizes(self, sizes, times):
        assert len(sizes) == len(times)
        demog = msprime.Demography.isolated_model([1234])
        for size, time in zip(sizes, times):
            demog.add_population_parameters_change(time, initial_size=size)
        graph = demog.to_demes()
        assert len(graph.demes) == 1
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 0
        epochs = graph.demes[0].epochs
        assert len(epochs) == len(sizes)
        for j in range(len(epochs)):
            assert epochs[j].size_function == "constant"
            assert epochs[j].start_size == sizes[len(sizes) - j - 1]
            assert epochs[j].end_size == sizes[len(sizes) - j - 1]
            if j == 0:
                assert math.isinf(epochs[j].start_time)
            else:
                assert epochs[j].start_time == times[len(times) - j]
            assert epochs[j].end_time == times[len(times) - j - 1]

    @pytest.mark.parametrize(
        "sizes,growth_rates,times",
        (
            ([1000, 500], [0.1, 0], [0, 20]),
            ([1000, None, 1000], [0.1, -0.1, 0], [0, 20, 30]),
            ([1000, 5000, 1000], [0.1, None, 0], [0, 20, 30]),
            ([1000, 5000, None], [0.1, None, 0], [0, 20, 30]),
        ),
    )
    def test_piecewise_exponential_sizes(self, sizes, growth_rates, times):
        assert growth_rates[-1] == 0  # can't convert to demes graph otherwise
        assert len(sizes) == len(growth_rates) == len(times)
        demog = msprime.Demography.isolated_model([1234])
        for size, growth_rate, time in zip(sizes, growth_rates, times):
            demog.add_population_parameters_change(
                time, initial_size=size, growth_rate=growth_rate
            )
        graph = demog.to_demes()
        assert len(graph.demes) == 1
        assert len(graph.migrations) == 0
        assert len(graph.pulses) == 0
        epochs = graph.demes[0].epochs
        assert len(epochs) == len(sizes)
        size = 4321
        growth_rate = 0
        start_times = times[1:] + [math.inf]
        for j, (start_time, end_time) in enumerate(zip(start_times, times)):
            epoch = epochs[len(epochs) - j - 1]
            if sizes[j] is not None:
                end_size = sizes[j]
            if growth_rates[j] is not None:
                growth_rate = growth_rates[j]

            if j == len(sizes) - 1:
                assert epoch.size_function == "constant"
                start_size = end_size
            else:
                assert epoch.size_function == "exponential"
                dt = start_time - end_time
                start_size = end_size * math.exp(dt * -growth_rate)

            assert math.isclose(epoch.start_size, start_size)
            assert math.isclose(epoch.end_size, end_size)

            if math.isinf(start_time):
                assert math.isinf(epoch.start_time)
            else:
                assert epoch.start_time == start_time
            assert epoch.end_time == end_time

            end_size = start_size

    def test_bad_growth_rate(self):
        demog = msprime.Demography.isolated_model([1000], growth_rate=[-1e-5])
        with pytest.raises(
            ValueError, match="growth rate for infinite-length epoch is invalid"
        ):
            demog.to_demes()

    def test_bad_population_size_change(self):
        demog = msprime.Demography.isolated_model([1000] * 2)
        demog.add_mass_migration(100, source="pop_0", dest="pop_1", proportion=1)
        demog.add_population_parameters_change(
            200, initial_size=500, population="pop_0"
        )
        with pytest.raises(ValueError, match="outside pop_0's existence interval"):
            demog.to_demes()

    def test_unsupported_event_simple_bottleneck(self):
        demog = msprime.Demography.isolated_model([1000])
        demog.add_simple_bottleneck(100, population=0, proportion=0.5)
        with pytest.raises(ValueError, match="Cannot convert"):
            demog.to_demes()

    def test_unsupported_event_instantaneous_bottleneck(self):
        demog = msprime.Demography.isolated_model([1000])
        demog.add_instantaneous_bottleneck(100, population=0, strength=100)
        with pytest.raises(ValueError, match="Cannot convert"):
            demog.to_demes()

    def test_census_event_is_ignored(self):
        demog = msprime.Demography.isolated_model([1000])
        demog.add_census(100)
        demog.to_demes()

    @pytest.mark.filterwarnings(
        # demes warns about multiple pulses with the same time
        "ignore:Multiple pulses:UserWarning"
    )
    def test_multipulse_roundtrip_two_sources(self):
        b = demes.Builder(defaults=dict(epoch=dict(start_size=1)))
        b.add_deme("a")
        b.add_deme("b")
        b.add_deme("c")
        b.add_pulse(sources=["a", "b"], proportions=[0.1, 0.2], dest="c", time=1)
        g = b.resolve()
        d = msprime.Demography.from_demes(g)
        g2 = d.to_demes()
        assert len(g2.pulses) == 2
        assert g2.pulses[0].sources[0] == "b"
        assert g2.pulses[0].proportions[0] == 0.2 / (1 - 0.1)
        assert g2.pulses[1].sources[0] == "a"
        assert g2.pulses[1].proportions[0] == 0.1

    @pytest.mark.filterwarnings(
        # demes warns about multiple pulses with the same time
        "ignore:Multiple pulses:UserWarning"
    )
    def test_multipulse_roundtrip_three_sources(self):
        b = demes.Builder(defaults=dict(epoch=dict(start_size=1)))
        b.add_deme("a")
        b.add_deme("b")
        b.add_deme("c")
        b.add_deme("d")
        b.add_pulse(
            sources=["a", "b", "c"], proportions=[0.1, 0.2, 0.3], dest="d", time=1
        )
        g = b.resolve()
        d = msprime.Demography.from_demes(g)
        g2 = d.to_demes()
        assert len(g2.pulses) == 3
        assert g2.pulses[0].sources[0] == "c"
        assert g2.pulses[0].proportions[0] == 0.3 / (1 - 0.1 - 0.2)
        assert g2.pulses[1].sources[0] == "b"
        assert g2.pulses[1].proportions[0] == 0.2 / (1 - 0.1)
        assert g2.pulses[2].sources[0] == "a"
        assert g2.pulses[2].proportions[0] == 0.1
