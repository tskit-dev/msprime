#
# Copyright (C) 2015-2021 University of Oxford
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
Test cases for basic ancestry simulation operations.
"""
import datetime
import json
import logging
import math
import random
import warnings

import numpy as np
import pytest
import tskit

import msprime
from msprime import _msprime
from msprime import ancestry


def tree_sequences_equal(ts1, ts2):
    """
    Returns True if the specified tree sequences are equal, ignoring
    their provenances.
    """
    t1 = ts1.dump_tables()
    t2 = ts2.dump_tables()
    t1.provenances.clear()
    t2.provenances.clear()
    return t1 == t2


def has_discrete_genome(ts):
    """
    Returns True if the specified tree sequence has discrete genome coordinates.
    """
    tables = ts.tables
    edges_left = np.all(tables.edges.left == np.floor(tables.edges.left))
    edges_right = np.all(tables.edges.right == np.floor(tables.edges.right))
    migrations_left = np.all(tables.migrations.left == np.floor(tables.migrations.left))
    migrations_right = np.all(
        tables.migrations.right == np.floor(tables.migrations.right)
    )
    sites = np.all(tables.sites.position == np.floor(tables.sites.position))
    return edges_left and edges_right and migrations_left and migrations_right and sites


def get_bottleneck_examples():
    """
    Returns an iterator of example tree sequences with nonbinary
    trees.
    """
    bottlenecks = [
        msprime.SimpleBottleneck(0.01, 0, proportion=0.05),
        msprime.SimpleBottleneck(0.02, 0, proportion=0.25),
        msprime.SimpleBottleneck(0.03, 0, proportion=1),
    ]
    for n in [3, 10, 100]:
        ts = msprime.simulate(
            n,
            length=100,
            recombination_rate=1,
            demographic_events=bottlenecks,
            random_seed=n,
        )
        yield ts


class TestFullArg:
    """
    Tests for recording the full ARG.
    """

    def verify(self, sim, multiple_mergers=False):
        sim.random_generator.seed = 1234
        sim.run()
        tree_sequence = next(sim.run_replicates(1))
        # Check if we have multiple merger somewhere.
        found = False
        for edgeset in tree_sequence.edgesets():
            if len(edgeset.children) > 2:
                found = True
                break
        assert multiple_mergers == found

        flags = tree_sequence.tables.nodes.flags
        time = tree_sequence.tables.nodes.time
        # TODO add checks for migrations.
        re_nodes = np.where(flags == msprime.NODE_IS_RE_EVENT)[0]
        ca_nodes = np.where(flags == msprime.NODE_IS_CA_EVENT)[0]
        coal_nodes = np.where(flags == 0)[0]
        # There should be two recombination nodes for every event
        assert np.array_equal(time[re_nodes[::2]], time[re_nodes[1::2]])  # Odd indexes
        assert re_nodes.shape[0] / 2 == sim.num_recombination_events
        if not multiple_mergers:
            assert (
                ca_nodes.shape[0] + coal_nodes.shape[0]
                == sim.num_common_ancestor_events
            )
        # After simplification, all the RE and CA nodes should be gone.
        ts_simplified = tree_sequence.simplify()
        new_flags = ts_simplified.tables.nodes.flags
        new_time = ts_simplified.tables.nodes.time
        assert np.sum(new_flags == msprime.NODE_IS_RE_EVENT) == 0
        assert np.sum(new_flags == msprime.NODE_IS_CA_EVENT) == 0
        # All coal nodes from the original should be identical to the originals
        assert np.array_equal(time[coal_nodes], new_time[new_flags == 0])
        assert ts_simplified.num_nodes <= tree_sequence.num_nodes
        assert ts_simplified.num_edges <= tree_sequence.num_edges
        return tree_sequence

    def test_no_recombination(self):
        sim = ancestry._parse_simulate(10, record_full_arg=True)
        ts = self.verify(sim)
        ts_simplified = ts.simplify()
        t1 = ts.tables
        t2 = ts_simplified.tables
        assert t1.nodes == t2.nodes
        assert t1.edges == t2.edges

    def test_recombination_n25(self):
        sim = ancestry._parse_simulate(
            25,
            recombination_rate=1,
            record_full_arg=True,
        )
        self.verify(sim)

    def test_recombination_n5(self):
        sim = ancestry._parse_simulate(
            5,
            recombination_rate=10,
            record_full_arg=True,
        )
        self.verify(sim)

    def test_recombination_n50(self):
        sim = ancestry._parse_simulate(
            50,
            recombination_rate=2,
            record_full_arg=True,
        )
        self.verify(sim)

    def test_recombination_n100(self):
        sim = ancestry._parse_simulate(
            100,
            recombination_rate=0.2,
            record_full_arg=True,
        )
        self.verify(sim)

    def test_multimerger(self):
        sim = ancestry._parse_simulate(
            100,
            recombination_rate=0.1,
            record_full_arg=True,
            demographic_events=[
                msprime.InstantaneousBottleneck(time=0.1, population=0, strength=5)
            ],
        )
        self.verify(sim, multiple_mergers=True)


class TestSimulator:
    """
    Runs tests on the underlying Simulator object.
    """

    def verify_simulation(self, n, m, r):
        """
        Verifies a simulation for the specified parameters.
        """
        sim = ancestry._parse_sim_ancestry(
            n,
            sequence_length=m,
            recombination_rate=r,
            random_seed=32,
        )
        sim.run()
        assert sim.num_breakpoints == len(sim.breakpoints)
        assert sim.time > 0
        assert sim.num_avl_node_blocks > 0
        assert sim.num_segment_blocks > 0
        assert sim.num_node_mapping_blocks > 0
        tree_sequence = next(sim.run_replicates(1))
        t = 0.0
        for record in tree_sequence.nodes():
            if record.time > t:
                t = record.time
        assert sim.time == t
        assert sim.num_common_ancestor_events > 0
        assert sim.num_recombination_events >= 0
        assert np.sum(sim.num_migration_events) >= 0
        assert sim.num_multiple_recombination_events >= 0

    def test_random_parameters(self):
        num_random_sims = 10
        for _ in range(num_random_sims):
            n = random.randint(2, 100)
            m = random.randint(10, 100)
            r = random.random()
            self.verify_simulation(n, m, r)

    def test_perf_parameters(self):
        sim = ancestry._parse_simulate(10)
        sim.run()
        assert sim.avl_node_block_size > 0
        assert sim.segment_block_size > 0
        assert sim.node_mapping_block_size > 0

    def test_event_chunk(self):
        sim = ancestry._parse_simulate(10)
        for bad_chunk in [-(2 ** 32), -1, 0]:
            with pytest.raises(ValueError):
                sim.run(event_chunk=bad_chunk)
        sim.reset()
        sim.run(event_chunk=2 ** 32 + 1)
        sim.reset()
        sim.run(event_chunk=2 ** 64 + 1)

    def test_debug_func(self):
        sim = ancestry._parse_simulate(10)
        count = 0

        def f(sim):
            nonlocal count
            count += 1

        sim.run(event_chunk=1, debug_func=f)
        assert count > 0

    def test_info_logging(self, caplog):
        sim = ancestry._parse_simulate(10)
        sim.random_generator.seed = 42
        with caplog.at_level(logging.INFO):
            sim.run()
        assert len(caplog.records) == 2
        assert (
            caplog.messages[0] == "Running model {'name': 'hudson'}"
            " until max time: inf"
        )

        assert caplog.messages[1].startswith("Completed at time")

    def test_debug_logging(self, caplog):
        sim = ancestry._parse_simulate(3)
        sim.random_generator.seed = 42
        with caplog.at_level(logging.DEBUG):
            sim.run(event_chunk=1)
        assert len(caplog.records) == 3
        assert (
            caplog.messages[0] == "Running model {'name': 'hudson'}"
            " until max time: inf"
        )
        assert caplog.messages[-1].startswith("Completed at time")
        assert caplog.messages[1] == "time=0.312845 ancestors=2"

    def test_debug_logging_dtwf(self, caplog):
        sim = ancestry._parse_simulate(3, Ne=10, model="dtwf")
        sim.random_generator.seed = 42
        with caplog.at_level(logging.DEBUG):
            sim.run(event_chunk=1)
            assert len(caplog.records) == 5
            assert (
                caplog.messages[0] == "Running model {'name': 'dtwf'}"
                " until max time: inf"
            )
            assert caplog.messages[-1].startswith("Completed at time")
            assert caplog.messages[1] == "time=1 ancestors=3"


class TestParseRandomSeed:
    """
    Tests for parsing the random seed values.
    """

    def test_default(self):
        # Make sure we get different random seeds when calling sequentially.
        seeds = [ancestry._parse_random_seed(None) for _ in range(100)]
        assert len(set(seeds)) == len(seeds)
        assert all(isinstance(seed, int) for seed in seeds)

    def test_numpy(self):
        seed = 12345
        ret_seed = ancestry._parse_random_seed(np.array([seed], dtype=int)[0])
        assert ret_seed == seed
        ret_seed = ancestry._parse_random_seed(np.array([seed], dtype=int))
        assert ret_seed == seed
        assert isinstance(ret_seed, int)

    def test_ints(self):
        # Anything that can be cast to an int is fine.
        for seed in [1234, 12.0, "12"]:
            ret_seed = ancestry._parse_random_seed(seed)
            assert ret_seed == int(seed)


class TestParseReplicateIndex:
    def test_none(self):
        assert (
            ancestry._parse_replicate_index(
                replicate_index=None, random_seed=None, num_replicates=None
            )
            is None
        )
        assert (
            ancestry._parse_replicate_index(
                replicate_index=None, random_seed=1, num_replicates=None
            )
            is None
        )
        assert (
            ancestry._parse_replicate_index(
                replicate_index=None, random_seed=1, num_replicates=1
            )
            is None
        )

    def test_random_seed(self):
        with pytest.raises(ValueError):
            ancestry._parse_replicate_index(
                replicate_index=0, random_seed=None, num_replicates=None
            )
        assert (
            ancestry._parse_replicate_index(
                replicate_index=1234, random_seed=1, num_replicates=None
            )
            == 1234
        )

    def test_num_replicates(self):
        with pytest.raises(ValueError):
            ancestry._parse_replicate_index(
                replicate_index=0, random_seed=1, num_replicates=1
            )

    def test_values(self):
        assert (
            ancestry._parse_replicate_index(
                replicate_index=0, random_seed=1, num_replicates=None
            )
            == 0
        )
        assert (
            ancestry._parse_replicate_index(
                replicate_index="1234", random_seed=1, num_replicates=None
            )
            == 1234
        )
        with pytest.raises(ValueError):
            ancestry._parse_replicate_index(
                replicate_index=-1, random_seed=1, num_replicates=None
            )


class TestParseSimAncestry:
    """
    Tests that the front-end for the sim_ancestry function correctly
    creates simulators with the required parameter values.
    """

    def test_sequence_length(self):
        # a single locus simulation will have sequence_length = 1
        sim = ancestry._parse_sim_ancestry(10)
        assert sim.sequence_length == 1
        assert sim.copy_tables().sequence_length == 1
        assert sim.recombination_map.total_mass == 0
        assert sim.gene_conversion_map.total_mass == 0

        # if we specify a rate_map for either GC or recomb this defines length.
        rate_map = msprime.RateMap.uniform(101, 0)
        sim = ancestry._parse_sim_ancestry(10, recombination_rate=rate_map)
        assert sim.sequence_length == rate_map.sequence_length
        assert sim.copy_tables().sequence_length == sim.sequence_length

        sim = ancestry._parse_sim_ancestry(
            10, gene_conversion_rate=rate_map, gene_conversion_tract_length=1
        )
        assert sim.sequence_length == rate_map.sequence_length
        assert sim.copy_tables().sequence_length == sim.sequence_length

        # If we have an initial_state this defines sequence_length
        initial_state = tskit.TableCollection(1234)
        initial_state.populations.add_row()
        sim = ancestry._parse_sim_ancestry(initial_state=initial_state)
        assert sim.sequence_length == 1234
        assert sim.copy_tables().sequence_length == sim.sequence_length

    def test_sequence_length_errors(self):
        # scaler rate values with no squence length is an error
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, recombination_rate=1)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10, gene_conversion_rate=1, gene_conversion_tract_length=1
            )

        # A rate map with a value that disagrees with sequence length
        rate_map = msprime.RateMap.uniform(101, 0)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10, recombination_rate=rate_map, sequence_length=1
            )
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10,
                gene_conversion_rate=rate_map,
                gene_conversion_tract_length=1,
                sequence_length=1,
            )

        # A different rate map with a sequence_length that disagrees
        other_rate_map = msprime.RateMap.uniform(1, 0)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10,
                recombination_rate=other_rate_map,
                gene_conversion_rate=rate_map,
                gene_conversion_tract_length=1,
            )
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10,
                recombination_rate=other_rate_map,
                gene_conversion_rate=rate_map,
                gene_conversion_tract_length=1,
                sequence_length=other_rate_map.sequence_length,
            )
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10,
                recombination_rate=other_rate_map,
                gene_conversion_rate=rate_map,
                gene_conversion_tract_length=1,
                sequence_length=rate_map.sequence_length,
            )
        # Both maps disagree with sequence_length
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10,
                recombination_rate=other_rate_map,
                gene_conversion_rate=rate_map,
                gene_conversion_tract_length=1,
                sequence_length=56789,
            )

        # An initial state with a sequence_length that disagrees.
        initial_state = tskit.TableCollection(1234).tree_sequence()
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(initial_state=initial_state, sequence_length=1)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                initial_state=initial_state, recombination_rate=rate_map
            )
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                initial_state=initial_state,
                gene_conversion_rate=rate_map,
                gene_conversion_tract_length=1,
            )
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                initial_state=initial_state,
                recombination_rate=other_rate_map,
                gene_conversion_rate=rate_map,
                gene_conversion_tract_length=1,
            )

    def test_sequence_length_discrete_genome(self):
        # Can't have floating point sequence_length with discrete_genome
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, sequence_length=1.1, discrete_genome=True)
        # Anything goes if we have a continuous genome, though.
        sim = ancestry._parse_sim_ancestry(
            10, sequence_length=0.1, discrete_genome=False
        )
        assert sim.sequence_length == 0.1

    def test_sequence_length_bad_arguments(self):
        for bad_value in ["x", b"sdf"]:
            with pytest.raises(ValueError):
                ancestry._parse_sim_ancestry(10, sequence_length=bad_value)

        for bad_type in [[], {}]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(10, sequence_length=bad_type)

    def test_gene_conversion_simple(self):
        for rate in ["1234", 1234, 1234.0]:
            sim = ancestry._parse_sim_ancestry(
                10,
                sequence_length=10,
                gene_conversion_rate=rate,
                gene_conversion_tract_length=5,
            )
            assert sim.sequence_length == 10
            gc_map = sim.gene_conversion_map
            assert gc_map.sequence_length == 10
            assert len(gc_map.rate) == 1
            assert gc_map.rate[0] == 1234
            assert sim.gene_conversion_tract_length == 5

    def test_gene_conversion_errors(self):
        # No tract length is an error
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10, sequence_length=10, gene_conversion_rate=1234
            )
        # Specifying a tract_length and no map is also an error.
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(
                10, sequence_length=10, gene_conversion_tract_length=5.5
            )

        for bad_type in [[], {}]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(
                    10,
                    sequence_length=10,
                    gene_conversion_rate=bad_type,
                    gene_conversion_tract_length=1,
                )

            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(
                    10,
                    sequence_length=10,
                    gene_conversion_rate=1,
                    gene_conversion_tract_length=bad_type,
                )

    def test_discrete_genome(self):
        # default is True
        sim = ancestry._parse_sim_ancestry(10, sequence_length=10)
        assert sim.discrete_genome
        for discrete_genome in [True, False]:
            sim = ancestry._parse_sim_ancestry(
                10, sequence_length=10, discrete_genome=discrete_genome
            )
            assert sim.discrete_genome == discrete_genome
        # Falsey values are not OK
        for discrete_genome in ["", []]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(
                    10, sequence_length=10, discrete_genome=discrete_genome
                )

    def test_start_time(self):
        # default is 0
        sim = ancestry._parse_sim_ancestry(10)
        assert sim.start_time == 0
        for start_time in [1234, 1234.34, "1", "1.234"]:
            sim = ancestry._parse_sim_ancestry(10, start_time=start_time)
            assert sim.start_time == float(start_time)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, start_time="bad value")
        with pytest.raises(TypeError):
            ancestry._parse_sim_ancestry(10, start_time=[])

    def test_end_time(self):
        # default is inf
        sim = ancestry._parse_sim_ancestry(10)
        assert sim.end_time == math.inf
        for end_time in [1234, 1234.34, "1", "1.234"]:
            sim = ancestry._parse_sim_ancestry(10, end_time=end_time)
            assert sim.end_time == float(end_time)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, end_time="bad value")
        with pytest.raises(TypeError):
            ancestry._parse_sim_ancestry(10, end_time=[])

    def test_record_migrations(self):
        # default is False
        sim = ancestry._parse_sim_ancestry(10)
        assert not sim.record_migrations
        for record_migrations in [True, False]:
            sim = ancestry._parse_sim_ancestry(10, record_migrations=record_migrations)
            assert sim.record_migrations == bool(record_migrations)
        for truthy in [0, []]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(10, record_migrations=truthy)

    def test_record_full_arg(self):
        # default is False
        sim = ancestry._parse_sim_ancestry(10)
        assert not sim.record_full_arg
        for record_full_arg in [True, False]:
            sim = ancestry._parse_sim_ancestry(10, record_full_arg=record_full_arg)
            assert sim.record_full_arg == bool(record_full_arg)

        for truthy in [0, []]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(10, record_full_arg=truthy)

    def test_ploidy(self):
        # default is 2
        sim = ancestry._parse_sim_ancestry(10)
        assert sim.ploidy == 2
        for ploidy in [1, 2, np.array([33])[0]]:
            sim = ancestry._parse_sim_ancestry(10, ploidy=ploidy)
            assert sim.ploidy == int(ploidy)

        for bad_ploidy in [0, -1]:
            with pytest.raises(ValueError):
                ancestry._parse_sim_ancestry(10, ploidy=bad_ploidy)

        for bad_ploidy in ["1", "0.1", 0.1, np.array([0.1])[0]]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(10, ploidy=bad_ploidy)

    def test_population_size(self):
        # default is 1
        sim = ancestry._parse_sim_ancestry(10)
        assert sim.demography.num_populations == 1
        assert sim.demography.populations[0].initial_size == 1
        assert sim.demography.populations[0].growth_rate == 0
        for pop_size in [2, 0.1, 100, 1e6, "100"]:
            sim = ancestry._parse_sim_ancestry(10, population_size=pop_size)
            assert sim.demography.num_populations == 1
            assert sim.demography.populations[0].initial_size == float(pop_size)
            assert sim.demography.populations[0].growth_rate == 0
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, population_size=-1)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, population_size="bad value")
        with pytest.raises(TypeError):
            ancestry._parse_sim_ancestry(10, population_size=[])

        # Cannot specify a population_size and demography args.
        demography = msprime.Demography.stepping_stone_model([1], 0)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, demography=demography, population_size=1)

    def test_demography(self):
        demography = msprime.Demography.stepping_stone_model([1] * 5, 0.1)
        samples = {0: 2}
        sim = ancestry._parse_sim_ancestry(samples, demography=demography)
        assert sim.demography is not demography
        assert sim.num_populations == demography.num_populations
        assert np.array_equal(sim.migration_matrix, demography.migration_matrix)
        # Numeric samples fail here as we have more than 1 population
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(5, demography=demography)
        with pytest.raises(TypeError):
            ancestry._parse_sim_ancestry(samples, demography="not a demography")

        demography.populations[0].initial_size = -1
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(samples, demography=demography)

    def test_model(self):
        # Extensive testing of the model parsing is done elsewhere.
        sim = ancestry._parse_sim_ancestry(10, model="smc")
        assert sim.model["name"] == "smc"
        sim = ancestry._parse_sim_ancestry(10, model=("smc", (10, "hudson")))
        assert sim.model["name"] == "smc"
        assert len(sim.model_change_events) == 1
        assert sim.model_change_events[0].time == 10

    def test_dtwf_population_size(self):
        # It's an error to not specify a pop size for dtwf.
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, model="dtwf", ploidy=2)

    def test_initial_state_errors(self):
        tables = tskit.TableCollection(1)
        tables.populations.add_row()
        # sequence_length doesn't match.
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(initial_state=tables, sequence_length=100)
        # Must have at least one population
        tables = tskit.TableCollection(1)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(initial_state=tables)
        for bad_type in [[], "sdf", {}]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(initial_state=bad_type)

    def test_initial_state(self):
        ts = msprime.sim_ancestry(10, end_time=0.01, random_seed=2)
        # Same if we use either the tables or tree sequence object.
        sim = ancestry._parse_sim_ancestry(initial_state=ts)
        assert sim.copy_tables() == ts.tables
        sim = ancestry._parse_sim_ancestry(initial_state=ts.tables)
        assert sim.copy_tables() == ts.tables

    def test_num_labels(self):
        for num_labels in [1, 2, 10]:
            sim = ancestry._parse_sim_ancestry(10, num_labels=num_labels)
            assert sim.num_labels == num_labels

    def test_samples_and_initial_state(self):
        # If we specify neither samples or initial_state we get an error
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(samples=None)

        # Specifying both is also an error.
        tables = tskit.TableCollection(1)
        tables.populations.add_row()
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(2, initial_state=tables)


class TestSimAncestrySamples:
    def test_negative_samples(self):
        for ploidy in [1, 2, 5]:
            with pytest.raises(ValueError):
                ancestry._parse_sim_ancestry(-1, ploidy=ploidy)

    def test_numeric_samples(self):
        for n in [1, 3, 10]:
            for ploidy in [1, 2, 3]:
                sim = ancestry._parse_sim_ancestry(n, ploidy=ploidy)
                assert sim.ploidy == ploidy
                tables = sim.copy_tables()
                assert len(tables.individuals) == n
                for individual in tables.individuals:
                    assert individual.flags == 0
                    assert len(individual.location) == 0
                assert len(tables.nodes) == n * ploidy
                assert len(tables.populations) == 1
                for node_id, node in enumerate(tables.nodes):
                    assert node.individual == node_id // ploidy
                    assert node.time == 0
                    assert node.flags == tskit.NODE_IS_SAMPLE
                    assert node.population == 0

    def test_sample_sets_override_ploidy(self):
        samples = [
            msprime.SampleSet(1, ploidy=2),
            msprime.SampleSet(1),  # No ploidy, should use default
            msprime.SampleSet(1, ploidy=1),
        ]
        sim = ancestry._parse_sim_ancestry(samples, ploidy=3)
        tables = sim.copy_tables()
        assert len(tables.individuals) == 3
        assert len(tables.nodes) == 6
        assert np.all(tables.nodes.individual[0:2] == 0)
        assert np.all(tables.nodes.individual[2:5] == 1)
        assert tables.nodes.individual[5] == 2

    def test_sample_sets_override_time(self):
        demography = msprime.Demography.isolated_model([1])
        demography.populations[0].sampling_time = 10
        samples = [
            msprime.SampleSet(1, ploidy=2, time=2),
            msprime.SampleSet(1, ploidy=3),  # No time, should use default
            msprime.SampleSet(1, ploidy=1, time=0),
        ]
        sim = ancestry._parse_sim_ancestry(samples, demography=demography)
        tables = sim.copy_tables()
        assert len(tables.individuals) == 3
        assert len(tables.nodes) == 6
        assert np.all(tables.nodes.time[0:2] == 2)
        assert np.all(tables.nodes.time[2:5] == 10)
        assert tables.nodes.time[5] == 0

    def test_numeric_samples_types(self):
        # Make sure the various different ways we can specify a numeric
        # set of samples all give the same answer
        values = [10, 10.0, np.array([10], dtype=int)[0]]
        for value in values:
            sim = ancestry._parse_sim_ancestry(value)
            assert len(sim.copy_tables().individuals) == 10

    def test_numeric_samples_only_simple_demography(self):
        # A simple 1-population model is fine.
        demography = msprime.Demography.isolated_model([1])
        sim = ancestry._parse_sim_ancestry(10, demography=demography)
        assert len(sim.copy_tables().individuals) == 10

        demography = msprime.Demography.stepping_stone_model([1, 1], 0)
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry(10, demography=demography)

    def test_population_name_sampling_errors(self):
        model = msprime.Demography.island_model([1] * 2, migration_rate=1)
        model.populations[0].name = "A"
        model.populations[1].name = "B"
        with pytest.raises(ValueError):
            ancestry._parse_sim_ancestry({}, demography=model)

        for bad_sample in [{"A": 1, "B": -1}, {"A": -1}, {"A": 0, "B": -10}]:
            with pytest.raises(ValueError):
                ancestry._parse_sim_ancestry(bad_sample, demography=model)
        for bad_type in [6.6, "a"]:
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry({"A": bad_type}, demography=model)
        for bad_name in ["C", "AC", "", "X" * 100]:
            with pytest.raises(KeyError):
                ancestry._parse_sim_ancestry({bad_name: 1}, demography=model)

    def test_population_name_samples_two_populations(self):
        model = msprime.Demography.island_model([1] * 2, migration_rate=1)
        model.populations[0].name = "A"
        model.populations[1].name = "B"
        sim = ancestry._parse_sim_ancestry({"A": 1}, demography=model)
        tables = sim.copy_tables()
        assert len(tables.nodes) == 2
        assert np.all(tables.nodes.population == 0)

        sim = ancestry._parse_sim_ancestry({"B": 1}, demography=model)
        tables = sim.copy_tables()
        assert len(tables.nodes) == 2
        assert np.all(tables.nodes.population == 1)

        sim = ancestry._parse_sim_ancestry({"A": 1, "B": 1}, demography=model)
        tables = sim.copy_tables()
        assert len(tables.nodes) == 4
        assert np.all(tables.nodes.population[:2] == 0)
        assert np.all(tables.nodes.population[2:] == 1)

        sim = ancestry._parse_sim_ancestry({"B": 1, "A": 1}, demography=model)
        tables = sim.copy_tables()
        assert len(tables.nodes) == 4
        assert np.all(tables.nodes.population[:2] == 1)
        assert np.all(tables.nodes.population[2:] == 0)

        sim = ancestry._parse_sim_ancestry({0: 1, "B": 1}, demography=model)
        tables = sim.copy_tables()
        assert len(tables.nodes) == 4
        assert np.all(tables.nodes.population[:2] == 0)
        assert np.all(tables.nodes.population[2:] == 1)

        sim = ancestry._parse_sim_ancestry({1: 1, 0: 1}, demography=model)
        tables = sim.copy_tables()
        assert len(tables.nodes) == 4
        assert np.all(tables.nodes.population[:2] == 1)
        assert np.all(tables.nodes.population[2:] == 0)

        # Referring to the same population by name and ID means we do the
        # sampling twice.
        sim = ancestry._parse_sim_ancestry({0: 1, "A": 1}, demography=model)
        tables = sim.copy_tables()
        assert len(tables.nodes) == 4
        assert np.all(tables.nodes.population[:4] == 0)

    def verify_samples_map(self, samples, demography, ploidy):
        demography = demography.validate()
        sim = ancestry._parse_sim_ancestry(
            samples=samples, demography=demography, ploidy=ploidy
        )
        n = sum(samples.values())
        assert sim.ploidy == ploidy
        tables = sim.copy_tables()
        assert len(tables.individuals) == n
        assert len(tables.nodes) == n * ploidy
        assert len(tables.populations) == demography.num_populations
        node_id = 0
        ind_id = 0
        for pop_id, num_samples in samples.items():
            for _ in range(num_samples):
                individual = tables.individuals[ind_id]
                for _ in range(ploidy):
                    node = tables.nodes[node_id]
                    assert node.individual == ind_id
                    assert node.time == demography.populations[pop_id].sampling_time
                    assert node.population == pop_id
                    assert node.flags == tskit.NODE_IS_SAMPLE
                    node_id += 1
                assert individual.flags == 0
                assert len(individual.location) == 0
                ind_id += 1

    def test_sample_demography(self):
        demography = msprime.Demography.isolated_model([1])
        self.verify_samples_map({0: 10}, demography, ploidy=1)
        self.verify_samples_map({0: 10}, demography, ploidy=2)

        demography = msprime.Demography.stepping_stone_model([1] * 5, 0)
        samples = {j: j for j in range(5)}
        self.verify_samples_map(samples, demography, ploidy=1)
        self.verify_samples_map(samples, demography, ploidy=2)

        samples = {4: 15}
        self.verify_samples_map(samples, demography, ploidy=1)
        self.verify_samples_map(samples, demography, ploidy=2)

    def test_sample_time_with_map(self):
        demography = msprime.Demography.stepping_stone_model([1, 1], 0)
        demography.populations[0].sampling_time = 1
        demography.populations[1].sampling_time = 2
        samples = {0: 5, 1: 5}
        self.verify_samples_map(samples, demography, ploidy=1)
        self.verify_samples_map(samples, demography, ploidy=2)

    def test_bad_list_samples_population(self):
        demography = msprime.Demography.stepping_stone_model([1, 1], 0)
        for bad_pop in [-1, 2]:
            samples = [msprime.SampleSet(2, time=0, population=bad_pop)]
            with pytest.raises(KeyError):
                ancestry._parse_sim_ancestry(samples=samples, demography=demography)
        for bad_pop in [1.1, ValueError]:
            samples = [msprime.SampleSet(2, time=0, population=bad_pop)]
            with pytest.raises(TypeError):
                ancestry._parse_sim_ancestry(samples=samples, demography=demography)

    def test_list_samples_no_population(self):
        # No population value is fine when we have one population
        samples = [msprime.SampleSet(2)]
        ts = ancestry.sim_ancestry(samples)
        assert ts.num_populations == 1

        # But not if we have 2 or more pops.
        demography = msprime.Demography.stepping_stone_model([1, 1], 0)
        with pytest.raises(ValueError, match="SampleSet population"):
            ancestry._parse_sim_ancestry(samples=samples, demography=demography)

    def test_sample_list_types(self):
        bad_sample_types = [ValueError, msprime.Sample(0, 0), "asdrf"]
        for bad_sample_type in bad_sample_types:
            with pytest.raises(TypeError, match="SampleSet object"):
                ancestry._parse_sim_ancestry([msprime.SampleSet(1), bad_sample_type])

    def test_bad_top_level_sample_types(self):
        bad_sample_types = [ValueError]
        for bad_sample_type in bad_sample_types:
            with pytest.raises(TypeError, match="different forms"):
                ancestry._parse_sim_ancestry(bad_sample_type)


class TestParseSimulate:
    """
    Tests that the front-end for the simulate function correctly
    creates simulators with the required parameter values.
    """

    def test_length(self):
        for bad_length in [-1, 0, -1e-6]:
            with pytest.raises(ValueError):
                ancestry._parse_simulate(10, length=bad_length)

    def test_num_labels(self):
        for bad_value in [-1, 0, 0.1]:
            with pytest.raises(ValueError):
                ancestry._parse_simulate(10, num_labels=bad_value)

    def test_sample_size(self):
        with pytest.raises(ValueError):
            ancestry._parse_simulate()
        with pytest.raises(ValueError):
            ancestry._parse_simulate(1)
        with pytest.raises(ValueError):
            ancestry._parse_simulate(sample_size=1)
        for n in [2, 100, 1000]:
            sim = ancestry._parse_simulate(n)
            tables = sim.copy_tables()
            assert len(tables.populations) == 1
            assert len(tables.individuals) == 0
            assert len(tables.edges) == 0
            assert len(tables.nodes) == n
            assert np.all(tables.nodes.flags == tskit.NODE_IS_SAMPLE)
            assert np.all(tables.nodes.time == 0)
            assert np.all(tables.nodes.individual == tskit.NULL)

    def test_effective_population_size(self):
        def f(Ne):
            return ancestry._parse_simulate(10, Ne=Ne)

        for bad_value in [-1, -1e16, 0]:
            with pytest.raises(ValueError):
                f(bad_value)
        for Ne in [1, 10, 1e5]:
            sim = f(Ne)
            assert sim.demography.populations[0].initial_size == Ne
        # Test the default.
        sim = ancestry._parse_simulate(10)
        assert sim.demography.populations[0].initial_size == 1

    def test_population_configurations(self):
        def f(configs):
            return ancestry._parse_simulate(population_configurations=configs)

        for bad_type in [10, ["sdf"], "sdfsd"]:
            with pytest.raises(TypeError):
                f(bad_type)
        # Just test the basic equalities here. The actual
        # configuration options are tested elewhere.
        for N in range(1, 10):
            pop_configs = [
                msprime.PopulationConfiguration(5, initial_size=5) for _ in range(N)
            ]
            sample_size = 5 * N
            sim = ancestry._parse_simulate(population_configurations=pop_configs)
            assert len(sim.demography.populations) == len(pop_configs)
            for pop, pop_config in zip(sim.demography.populations, pop_configs):
                assert pop.initial_size == pop_config.initial_size
                assert pop.growth_rate == pop_config.growth_rate
            tables = sim.copy_tables()
            assert len(tables.nodes) == sample_size
            assert len(sim.population_configuration) == N
        # The default is a single population
        sim = ancestry._parse_simulate(10)
        assert len(sim.population_configuration) == 1

    def test_sample_size_population_configuration(self):
        for d in range(1, 5):
            # Zero sample size is always an error
            configs = [msprime.PopulationConfiguration(0) for _ in range(d)]
            with pytest.raises(ValueError):
                ancestry._parse_simulate(population_configurations=configs)
            configs = [msprime.PopulationConfiguration(2) for _ in range(d)]
            sim = ancestry._parse_simulate(population_configurations=configs)
            tables = sim.copy_tables()
            assert len(tables.nodes) == 2 * d
            i = 0
            for j in range(d):
                for _ in range(2):
                    node = tables.nodes[i]
                    assert node.population == j
                    assert node.time == 0
                    assert node.flags == tskit.NODE_IS_SAMPLE
                    i += 1

    def test_migration_matrix(self):
        # Cannot specify a migration matrix without population
        # configurations
        with pytest.raises(ValueError):
            ancestry._parse_simulate(10, migration_matrix=[])
        for N in range(1, 10):
            pop_configs = [msprime.PopulationConfiguration(5) for _ in range(N)]
            sim = ancestry._parse_simulate(population_configurations=pop_configs)
            # If we don't specify a matrix, it's 0 everywhere.
            matrix = np.zeros((N, N))
            np.testing.assert_array_equal(sim.migration_matrix, matrix)

            def f(matrix):
                return ancestry._parse_simulate(
                    population_configurations=pop_configs, migration_matrix=matrix
                )

            matrix = [[(j + k) * int(j != k) for j in range(N)] for k in range(N)]
            sim = f(matrix)
            np.testing.assert_array_equal(sim.demography.migration_matrix, matrix)
            # Try with equivalent numpy array.
            sim = f(np.array(matrix))
            np.testing.assert_array_equal(sim.demography.migration_matrix, matrix)
            np.testing.assert_array_equal(sim.migration_matrix, matrix)
            for bad_type in [{}, "", 234, 1.2]:
                with pytest.raises(ValueError):
                    f(bad_type)
            # Now check for the structure of the matrix.
            matrix[0][0] = "bad value"
            with pytest.raises(TypeError):
                f(matrix)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                matrix[0] = None
                with pytest.raises(ValueError):
                    f(matrix)
                matrix[0] = []
                with pytest.raises(ValueError):
                    f(matrix)
            # Simple numpy array.
            matrix = np.ones((N, N))
            np.fill_diagonal(matrix, 0)
            sim = f(matrix)
            np.testing.assert_array_equal(
                np.array(sim.demography.migration_matrix), matrix
            )
            sim.run()
            events = np.array(sim.num_migration_events)
            assert events.shape == (N, N)
            assert np.all(events >= 0)

    def test_default_migration_matrix(self):
        sim = ancestry._parse_simulate(10)
        assert sim.migration_matrix == [0.0]

    def test_demographic_events(self):
        for bad_type in ["sdf", 234, [12], [None]]:
            with pytest.raises(TypeError):
                ancestry._parse_simulate(2, demographic_events=bad_type)
        # TODO test for bad values.

    def test_recombination_rate(self):
        def f(recomb_rate):
            return ancestry._parse_simulate(10, recombination_rate=recomb_rate)

        for bad_type in ["", {}, []]:
            with pytest.raises(TypeError):
                f(bad_type)
        for bad_value in [-1, -1e15]:
            with pytest.raises(ValueError):
                f(bad_value)
        for rate in [0, 1e-3, 10]:
            sim = f(rate)
            recomb_map = sim.recombination_map
            assert list(recomb_map.position) == [0, 1]
            assert list(recomb_map.rate) == [rate]
            assert sim.sequence_length == recomb_map.sequence_length

    def test_recombination_map(self):
        def f(recomb_map):
            return ancestry._parse_simulate(10, recombination_map=recomb_map)

        with pytest.raises(TypeError):
            f("wrong type")
        for n in range(2, 10):
            positions = list(range(n))
            rates = [0.1 * j for j in range(n - 1)]
            # Use the old-form RecombinationMap
            recomb_map = msprime.RecombinationMap(positions, rates + [0.0])
            sim = ancestry._parse_simulate(10, recombination_map=recomb_map)
            other_map = sim.recombination_map
            assert list(other_map.position) == positions
            assert list(other_map.rate) == rates
            assert sim.sequence_length == other_map.sequence_length
            # Use the new-form RateMap
            rate_map = msprime.RateMap(positions, rates)
            sim = ancestry._parse_simulate(10, recombination_map=rate_map)
            other_map = sim.recombination_map
            assert list(other_map.position) == positions
            assert list(other_map.rate) == rates
            assert sim.sequence_length == other_map.sequence_length

    def test_combining_recomb_map_and_rate_length(self):
        recomb_map = msprime.RecombinationMap([0, 1], [1, 0])
        with pytest.raises(ValueError):
            ancestry._parse_simulate(
                10,
                recombination_map=recomb_map,
                length=1,
            )
        with pytest.raises(ValueError):
            ancestry._parse_simulate(
                10,
                recombination_map=recomb_map,
                recombination_rate=100,
            )
        with pytest.raises(ValueError):
            ancestry._parse_simulate(
                10,
                recombination_map=recomb_map,
                length=1,
                recombination_rate=1,
            )

    def test_sample_combination_errors(self):
        # Make sure that the various ways we can specify the samples
        # operate correctly.
        s = msprime.Sample(time=0.0, population=0)
        with pytest.raises(ValueError):
            ancestry._parse_simulate()
        # Cannot provide sample_size with either population configurations
        # or samples
        with pytest.raises(ValueError):
            ancestry._parse_simulate(sample_size=2, samples=[s, s])
        pop_configs = [msprime.PopulationConfiguration(sample_size=2)]
        with pytest.raises(ValueError):
            ancestry._parse_simulate(
                sample_size=2,
                population_configurations=pop_configs,
            )
        # If we provide samples and population_configurations we cannot
        # have a sample size for the config.
        pop_configs = [msprime.PopulationConfiguration(sample_size=2)]
        with pytest.raises(ValueError):
            ancestry._parse_simulate(
                samples=[s, s],
                population_configurations=pop_configs,
            )
        pop_configs = [
            msprime.PopulationConfiguration(sample_size=None),
            msprime.PopulationConfiguration(sample_size=2),
        ]
        with pytest.raises(ValueError):
            ancestry._parse_simulate(
                samples=[s, s],
                population_configurations=pop_configs,
            )

    def test_samples(self):
        pop_configs = [
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
            msprime.PopulationConfiguration(),
        ]
        samples = [
            msprime.Sample(population=0, time=0),
            msprime.Sample(population=1, time=1),
            msprime.Sample(population=2, time=2),
        ]
        sim = ancestry._parse_simulate(
            samples=samples, population_configurations=pop_configs
        )
        tables = sim.copy_tables()
        assert len(tables.nodes) == len(samples)
        for node, sample in zip(tables.nodes, samples):
            assert node.population == sample.population
            assert node.time == sample.time
            assert node.flags == tskit.NODE_IS_SAMPLE

    def test_sample_is_tuple(self):
        s = msprime.Sample(population=1, time=2)
        assert s.population == 1
        assert s.time == 2
        assert s == (1, 2)

    def test_new_old_style_model_changes_equal(self):
        models = [
            msprime.SweepGenicSelection(
                position=j, start_frequency=j, end_frequency=j, s=j, dt=j
            )
            for j in range(1, 10)
        ]
        # Old style
        sim = ancestry._parse_simulate(
            sample_size=2,
            Ne=10,
            demographic_events=[
                msprime.SimulationModelChange(None, model) for model in models
            ],
        )
        assert len(sim.model_change_events) == len(models)
        for event, model in zip(sim.model_change_events, models):
            assert event.model == model

        sim2 = ancestry._parse_simulate(
            sample_size=2,
            Ne=10,
            model=[None]
            + [msprime.SimulationModelChange(None, model) for model in models],
        )
        assert sim.model_change_events == sim2.model_change_events

    def test_model_change_old_style(self):
        main_model = msprime.SmcApproxCoalescent()
        sim = ancestry._parse_simulate(
            Ne=100,
            sample_size=2,
            model=main_model,
            demographic_events=[
                msprime.SimulationModelChange(1, msprime.DiscreteTimeWrightFisher()),
                msprime.SimulationModelChange(2, None),
            ],
        )
        assert len(sim.model_change_events) == 2
        assert sim.model_change_events[0].time == 1
        # When model=None we change to the standard coalescent
        assert sim.model_change_events[1].time == 2
        assert sim.model_change_events[1].model.name == "hudson"

        # This should be the same in new notation
        sim = ancestry._parse_simulate(
            Ne=100, sample_size=2, model=[main_model, (1, "dtwf"), (2, None)]
        )
        assert len(sim.model_change_events) == 2
        assert sim.model_change_events[0].time == 1
        # When model=None we change to the standard coalescent
        assert sim.model_change_events[1].time == 2
        assert sim.model_change_events[1].model.name == "hudson"

    def test_bad_sample_population_reference(self):
        # What happens when we reference a population that doesn't exist?
        with pytest.raises(
            ValueError, match="Invalid population reference '1' in sample at index 1"
        ):
            msprime.simulate(
                samples=[
                    msprime.Sample(population=0, time=0),
                    msprime.Sample(population=1, time=0),
                ]
            )

        with pytest.raises(
            ValueError, match="Negative population ID in sample at index 2"
        ):
            msprime.simulate(
                samples=[
                    msprime.Sample(population=0, time=0),
                    msprime.Sample(population=0, time=0),
                    msprime.Sample(population=-1, time=0),
                ]
            )


class TestSimAncestryInterface:
    """
    Some simple tests cases for the sim_ancestry interface.
    """

    def test_defaults(self):
        n = 10
        # Diploid sim by default.
        ts = msprime.sim_ancestry(n)
        assert ts.num_samples == 2 * n
        assert ts.num_individuals == n
        assert ts.num_trees == 1
        assert ts.num_sites == 0
        assert ts.sequence_length == 1

    def test_ploidy(self):
        n = 10
        for k in [1, 2, 3, 4]:
            ts = msprime.sim_ancestry(n, ploidy=k)
            assert ts.num_samples == k * n
            assert ts.num_trees == 1
            assert ts.num_sites == 0
            assert ts.sequence_length == 1
            # TODO check for individuals

    def test_hudson_time_scale(self):
        n = 10
        seed = 1234
        for ploidy in [1, 2, 3, 7]:
            ts1 = msprime.sim_ancestry(n * ploidy, ploidy=1, random_seed=seed)
            ts2 = msprime.sim_ancestry(n, ploidy=ploidy, random_seed=seed)
            t1 = ts1.tables
            t2 = ts2.tables
            assert np.allclose(t1.nodes.time * ploidy, t2.nodes.time)
            assert t1.edges == t2.edges

    def test_ploidy_demography(self):
        n = 2
        demography = msprime.Demography.island_model([100] * 2, 0.1)
        for k in [1, 2, 3, 4]:
            ts = msprime.sim_ancestry(
                samples={0: n, 1: n}, ploidy=k, demography=demography
            )
            assert ts.num_samples == 2 * n * k
            assert ts.num_trees == 1
            assert ts.num_sites == 0
            assert ts.sequence_length == 1
            assert ts.num_populations == 2
            assert ts.num_individuals == 2 * n

    def test_random_seed(self):
        ts1 = msprime.sim_ancestry(10, random_seed=1)
        ts2 = msprime.sim_ancestry(10, random_seed=1)
        assert tree_sequences_equal(ts1, ts2)

        ts2 = msprime.sim_ancestry(10, random_seed=2)
        assert not tree_sequences_equal(ts1, ts2)

    def test_population_size(self):
        ts1 = msprime.sim_ancestry(10, population_size=1, random_seed=2)
        # Defaults to 1
        ts2 = msprime.sim_ancestry(10, random_seed=2)
        assert tree_sequences_equal(ts1, ts2)
        ts2 = msprime.sim_ancestry(10, population_size=100, random_seed=2)
        # Acts as a simple scaling factor on times.
        assert ts1.tables.edges == ts2.tables.edges
        assert np.allclose(100 * ts1.tables.nodes.time, ts2.tables.nodes.time)

    def test_replicates(self):
        ts = msprime.sim_ancestry(10)
        assert isinstance(ts, tskit.TreeSequence)
        for n in [0, 1, 2, 5]:
            ts_list = list(msprime.sim_ancestry(10, num_replicates=n))
            assert len(ts_list) == n
            for ts in ts_list:
                assert isinstance(ts, tskit.TreeSequence)
                assert ts.num_individuals == 10
                assert ts.num_samples == 20
                assert ts.num_trees == 1

    def test_recombination_rate(self):
        ts = msprime.sim_ancestry(10, recombination_rate=1, sequence_length=10)
        assert ts.num_samples == 20
        assert ts.sequence_length == 10
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)
        # A non-zero recombination_rate and no sequence length is an error
        with pytest.raises(ValueError):
            msprime.sim_ancestry(10, recombination_rate=1)
        # But if we specify a rate map, that's OK.
        rate_map = msprime.RateMap.uniform(sequence_length=10, rate=1)
        ts = msprime.sim_ancestry(10, recombination_rate=rate_map)
        assert ts.num_samples == 20
        assert ts.sequence_length == 10
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)

        # We should get precisely the same ts if we have the same seed
        ts1 = msprime.sim_ancestry(
            10, recombination_rate=1, sequence_length=10, random_seed=1
        )
        ts2 = msprime.sim_ancestry(10, recombination_rate=rate_map, random_seed=1)
        assert tree_sequences_equal(ts1, ts2)

    def test_gene_conversion_rate(self):
        ts = msprime.sim_ancestry(
            10,
            gene_conversion_rate=1,
            gene_conversion_tract_length=2,
            sequence_length=10,
            random_seed=14,
        )
        assert ts.num_trees > 1
        for tree in ts.trees():
            assert tree.num_roots == 1

    def test_gene_conversion_simple_map(self):
        n = 10
        ts = msprime.sim_ancestry(
            n,
            gene_conversion_rate=1,
            gene_conversion_tract_length=1,
            sequence_length=10,
            recombination_rate=1,
            discrete_genome=True,
        )
        assert isinstance(ts, tskit.TreeSequence)
        assert ts.num_samples == 2 * n
        assert ts.num_trees > 1

    def test_gene_conversion_continuous(self):
        with pytest.raises(ValueError):
            msprime.sim_ancestry(
                10,
                gene_conversion_rate=1,
                gene_conversion_tract_length=1,
                discrete_genome=False,
            )

    def test_gene_conversion_default_map(self):
        n = 10
        ts = msprime.sim_ancestry(
            n,
            sequence_length=10,
            gene_conversion_rate=1,
            gene_conversion_tract_length=1,
            discrete_genome=True,
        )
        assert isinstance(ts, tskit.TreeSequence)
        assert ts.num_samples == 2 * n
        assert ts.num_trees > 1

    def test_model(self):
        ts1 = msprime.sim_ancestry(10, population_size=100, random_seed=2)
        ts2 = msprime.sim_ancestry(
            10, population_size=100, model="hudson", random_seed=2
        )
        assert tree_sequences_equal(ts1, ts2)
        ts2 = msprime.sim_ancestry(
            10, population_size=100, model="dtwf", random_seed=2, ploidy=2
        )
        assert not tree_sequences_equal(ts1, ts2)

    def test_continuous_genome(self):
        ts = msprime.sim_ancestry(
            10, recombination_rate=10, sequence_length=1, discrete_genome=False
        )
        assert ts.num_samples == 20
        assert ts.sequence_length == 1
        assert ts.num_trees > 1
        assert not has_discrete_genome(ts)

    def test_discrete_genome(self):
        # Default to discrete_genome=True
        ts = msprime.sim_ancestry(
            10, recombination_rate=10, sequence_length=10, random_seed=2
        )
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)

        ts = msprime.sim_ancestry(
            10,
            recombination_rate=10,
            sequence_length=10,
            random_seed=2,
            discrete_genome=True,
        )
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)

        ts = msprime.sim_ancestry(
            10,
            recombination_rate=1,
            sequence_length=10,
            random_seed=2,
            discrete_genome=False,
        )
        assert ts.num_trees > 1
        assert not has_discrete_genome(ts)

    def test_discrete_genome_migrations(self):
        def run_sim(discrete_genome=None):
            demography = msprime.Demography.stepping_stone_model([1, 1], 0.1)
            return msprime.sim_ancestry(
                samples={0: 5, 1: 5},
                demography=demography,
                sequence_length=5,
                recombination_rate=1,
                discrete_genome=discrete_genome,
                record_migrations=True,
                random_seed=2134,
            )

        ts_discrete = run_sim(True)
        assert ts_discrete.num_trees > 1
        assert ts_discrete.num_migrations > 1
        assert has_discrete_genome(ts_discrete)

        ts_continuous = run_sim(False)
        assert ts_continuous.num_trees > 1
        assert ts_continuous.num_migrations > 1
        assert not has_discrete_genome(ts_continuous)

    def test_numpy_random_seed(self):
        seed = np.array([12345], dtype=np.int64)[0]
        assert seed.dtype == np.int64
        ts1 = msprime.sim_ancestry(10, random_seed=seed)
        ts2 = msprime.sim_ancestry(10, random_seed=seed)
        assert ts1.tables.nodes == ts2.tables.nodes

    def test_record_provenance(self):
        # The content of the provenances is tested elsewhere.
        ts = msprime.sim_ancestry(10, random_seed=2)
        assert ts.num_provenances == 1
        ts = msprime.sim_ancestry(10, random_seed=2, record_provenance=False)
        assert ts.num_provenances == 0

    def test_dtwf(self):
        ts = msprime.sim_ancestry(
            10, population_size=100, model="dtwf", ploidy=2, random_seed=1234
        )
        assert ts.num_trees == 1
        assert ts.first().num_roots == 1
        # All node times should be integers
        time = ts.tables.nodes.time
        assert np.all(time == np.floor(time))

    def test_dtwf_non_diploid(self):
        for ploidy in [1, 3, 7]:
            with pytest.raises(_msprime.LibraryError):
                msprime.sim_ancestry(
                    10, population_size=100, model="dtwf", ploidy=ploidy
                )

    def test_sweep_coalescence(self):
        N = 1e6
        model = msprime.SweepGenicSelection(
            position=0.5,
            start_frequency=1.0 / (2 * N),
            end_frequency=1.0 - (1.0 / (2 * N)),
            s=1000,
            dt=1e-6,
        )
        ts = msprime.sim_ancestry(10, model=model)
        assert ts.num_trees == 1
        assert ts.first().num_roots == 1

    def test_start_time(self):
        ts = msprime.sim_ancestry(10, ploidy=1, random_seed=42, start_time=100)
        assert np.all(ts.tables.nodes.time[10:] > 100)

    def test_end_time(self):
        ts = msprime.sim_ancestry(10, random_seed=42, end_time=0.01)
        assert np.all(ts.tables.nodes.time <= 0.01)
        assert ts.first().num_roots > 1

    def test_record_migrations(self):
        demography = msprime.Demography.stepping_stone_model([1, 1], 0.1)
        samples = {0: 2, 1: 2}
        ts = msprime.sim_ancestry(samples, demography=demography, random_seed=42)
        assert ts.first().num_roots == 1
        # Migrations are off by default
        assert ts.num_migrations == 0

        ts = msprime.sim_ancestry(
            samples, demography=demography, random_seed=42, record_migrations=True
        )
        assert ts.first().num_roots == 1
        # Migrations are off by default
        assert ts.num_migrations > 0

        ts = msprime.sim_ancestry(
            samples, demography=demography, random_seed=42, record_migrations=False
        )
        assert ts.first().num_roots == 1
        assert ts.num_migrations == 0

    def test_record_full_arg(self):
        ts = msprime.sim_ancestry(
            4,
            recombination_rate=1,
            random_seed=2,
            sequence_length=10,
            record_full_arg=True,
        )
        flags = ts.tables.nodes.flags
        assert np.sum(flags == msprime.NODE_IS_RE_EVENT) > 0
        for record_full_arg in [None, False]:
            ts = msprime.sim_ancestry(
                4,
                recombination_rate=1,
                random_seed=2,
                sequence_length=10,
                record_full_arg=record_full_arg,
            )
            flags = ts.tables.nodes.flags
            assert np.sum(flags == msprime.NODE_IS_RE_EVENT) == 0

    def test_initial_tables_recapitate(self):
        # Simple recapitate scenario
        ts = msprime.sim_ancestry(5, end_time=0.5, random_seed=53)
        assert ts.first().num_roots > 1
        recapitated1 = msprime.sim_ancestry(initial_state=ts, random_seed=234)
        assert recapitated1.num_trees == 1
        assert recapitated1.first().num_roots == 1

        # We should get the same answer from the providing the tables argument.
        recapitated2 = msprime.sim_ancestry(initial_state=ts.tables, random_seed=234)
        assert tree_sequences_equal(recapitated1, recapitated2)

    def test_replicate_index(self):
        reps = msprime.sim_ancestry(12, random_seed=52, num_replicates=10)
        for index, ts1 in enumerate(reps):
            ts2 = msprime.sim_ancestry(12, random_seed=52, replicate_index=index)
            t1 = ts1.dump_tables()
            t2 = ts2.dump_tables()
            t1.provenances.clear()
            t2.provenances.clear()
            assert t1 == t2


class TestSimulateInterface:
    """
    Some simple test cases for the simulate() interface.
    """

    def test_defaults(self):
        n = 10
        ts = msprime.simulate(n)
        assert isinstance(ts, tskit.TreeSequence)
        assert ts.get_sample_size() == n
        assert ts.get_num_trees() == 1
        assert ts.get_num_mutations() == 0
        assert ts.get_sequence_length() == 1
        assert len(list(ts.provenances())) == 1

    def test_positional_args_not_allowed(self):
        with pytest.raises(TypeError):
            msprime.simulate(2, 100)

    def verify_provenance(self, provenance):
        """
        Checks that the specified provenance object has the right sort of
        properties.
        """
        # Generate the ISO 8601 time for now, without the high precision suffix,
        # and compare the prefixes.
        today = datetime.date.today().isoformat()
        k = len(today)
        assert provenance.timestamp[:k] == today
        assert provenance.timestamp[k] == "T"
        d = json.loads(provenance.record)
        assert len(d) > 0

    def test_provenance(self):
        ts = msprime.simulate(10)
        assert ts.num_provenances == 1
        self.verify_provenance(ts.provenance(0))
        for ts in msprime.simulate(10, num_replicates=10):
            assert ts.num_provenances == 1
            self.verify_provenance(ts.provenance(0))

    def test_end_time(self):
        ts = msprime.simulate(15, recombination_rate=2, random_seed=8, end_time=0.1)
        for tree in ts.trees():
            for root in tree.roots:
                assert tree.time(root) == 0.1

    def test_replicates(self):
        n = 20
        num_replicates = 10
        count = 0
        for ts in msprime.simulate(n, num_replicates=num_replicates):
            count += 1
            assert isinstance(ts, tskit.TreeSequence)
            assert ts.get_sample_size() == n
            assert ts.get_num_trees() == 1
        assert num_replicates == count

    def test_mutations(self):
        n = 10
        ts = msprime.simulate(n, mutation_rate=10)
        assert isinstance(ts, tskit.TreeSequence)
        assert ts.get_sample_size() == n
        assert ts.get_num_trees() == 1
        assert ts.get_num_mutations() > 0

    def test_no_mutations_with_start_time(self):
        with pytest.raises(ValueError):
            msprime.simulate(10, mutation_rate=10, start_time=3)
        # But fine if we set start_time = None
        ts = msprime.simulate(10, mutation_rate=10, start_time=None, random_seed=1)
        assert ts.num_sites > 0

    def test_mutation_generator_unsupported(self):
        with pytest.raises(ValueError):
            msprime.simulate(10, mutation_generator="some non-None value")

    def test_mutation_interface(self):
        for bad_type in [{}, [], self]:
            with pytest.raises(TypeError):
                msprime.simulate(10, mutation_rate=bad_type)
        for bad_value in ["x", "234x"]:
            with pytest.raises(ValueError):
                msprime.simulate(10, mutation_rate=bad_value)

    def test_recombination(self):
        n = 10
        ts = msprime.simulate(n, recombination_rate=10)
        assert isinstance(ts, tskit.TreeSequence)
        assert ts.sample_size == n
        assert ts.num_trees > 1
        assert ts.num_mutations == 0

    def test_num_labels(self):
        # Running simulations with different numbers of labels in the default
        # setting should have no effect.
        tables = [
            msprime.simulate(10, num_labels=num_labels, random_seed=1).tables
            for num_labels in range(1, 5)
        ]
        for t in tables:
            t.provenances.clear()
        for t in tables:
            assert t == tables[0]

    def test_replicate_index(self):
        reps = msprime.simulate(10, random_seed=52, num_replicates=10)
        for index, ts1 in enumerate(reps):
            ts2 = msprime.simulate(10, random_seed=52, replicate_index=index)
            t1 = ts1.dump_tables()
            t2 = ts2.dump_tables()
            t1.provenances.clear()
            t2.provenances.clear()
            assert t1 == t2


class TestDiscreteGenomeSimulate:
    """
    Test that the num_loci==sequence_length trick for doing a discrete genome
    simulation in simulate() still works.
    """

    def test_simple_example(self):
        rmap = msprime.RecombinationMap.uniform_map(100, 0.1, num_loci=100)
        ts = msprime.simulate(10, recombination_map=rmap, random_seed=1)
        assert has_discrete_genome(ts)
        assert ts.sequence_length == 100
        assert ts.num_trees > 1

    def test_single_locus(self):
        rmap = msprime.RecombinationMap.uniform_map(1, 100, num_loci=1)
        ts = msprime.simulate(10, recombination_map=rmap, random_seed=1)
        assert has_discrete_genome(ts)
        assert ts.sequence_length == 1
        assert ts.num_trees == 1

    @pytest.mark.parametrize("m", [2, 5, 7])
    def test_small_examples(self, m):
        rmap = msprime.RecombinationMap.uniform_map(m, 100, num_loci=m)
        ts = msprime.simulate(10, recombination_map=rmap, random_seed=1)
        assert has_discrete_genome(ts)
        assert ts.sequence_length == m
        assert ts.num_trees == m

    def test_simulate_from_example(self):
        # Example from the docs on a WF simulation
        N = 10
        num_loci = 2
        recomb_map = msprime.RecombinationMap.uniform_map(
            num_loci, 1 / num_loci, num_loci
        )
        ts = msprime.simulate(5, Ne=N / 2, recombination_map=recomb_map, random_seed=5)
        assert has_discrete_genome(ts)
        assert ts.sequence_length == num_loci
        assert ts.num_trees > 1


class TestReprRoundTrip:
    """
    Tests that we can eval the repr of objects to round trip them.
    """

    def assert_repr_round_trip(self, obj_list):
        for obj in obj_list:
            obj_copy = eval(repr(obj), globals(), msprime.__dict__)
            assert obj_copy == obj
            assert not (obj_copy is obj)

    def test_population(self):
        examples = [
            msprime.Population(initial_size=2),
            msprime.Population(2, growth_rate=5),
            msprime.Population(initial_size=234, growth_rate=10),
        ]
        self.assert_repr_round_trip(examples)

    def test_population_parameters_change(self):
        examples = [
            msprime.PopulationParametersChange(time=1, initial_size=1),
            msprime.PopulationParametersChange(time=1, growth_rate=2),
            msprime.PopulationParametersChange(time=1, growth_rate=1, population=2),
            msprime.PopulationParametersChange(
                time=3, initial_size=3, growth_rate=1, population=2
            ),
        ]
        self.assert_repr_round_trip(examples)

    def test_migration_rate_change(self):
        examples = [
            msprime.MigrationRateChange(time=1, rate=1),
            msprime.MigrationRateChange(time=1, rate=1, source=1, dest=2),
        ]
        self.assert_repr_round_trip(examples)

    def test_mass_migration(self):
        examples = [
            msprime.MassMigration(time=1, source=1, dest=2),
            msprime.MassMigration(time=1, source=1, dest=2, proportion=0.2),
        ]
        self.assert_repr_round_trip(examples)

    def test_population_split(self):
        examples = [
            msprime.PopulationSplit(time=1, derived=[1, 2], ancestral=3),
        ]
        self.assert_repr_round_trip(examples)

    def test_simulation_model_change(self):
        examples = [
            msprime.SimulationModelChange(),
            msprime.SimulationModelChange(model="hudson"),
            msprime.SimulationModelChange(model=msprime.DiscreteTimeWrightFisher()),
            msprime.SimulationModelChange(
                model=msprime.BetaCoalescent(alpha=1, truncation_point=2)
            ),
        ]
        self.assert_repr_round_trip(examples)

    def test_simple_bottleneck(self):
        examples = [
            msprime.SimpleBottleneck(time=10, population=2),
            msprime.SimpleBottleneck(time=10, population=2, proportion=0.5),
        ]
        self.assert_repr_round_trip(examples)

    def test_instantaneous_bottleneck(self):
        examples = [
            msprime.InstantaneousBottleneck(time=10, population=1),
            msprime.InstantaneousBottleneck(time=10, population=1, strength=10),
        ]
        self.assert_repr_round_trip(examples)

    def test_census_event(self):
        examples = [
            msprime.CensusEvent(time=10),
        ]
        self.assert_repr_round_trip(examples)

    def test_simulation_models(self):
        examples = [
            msprime.StandardCoalescent(),
            msprime.SmcApproxCoalescent(),
            msprime.SmcPrimeApproxCoalescent(),
            msprime.DiscreteTimeWrightFisher(),
            msprime.WrightFisherPedigree(),
            msprime.BetaCoalescent(),
            msprime.BetaCoalescent(alpha=1, truncation_point=10),
            msprime.DiracCoalescent(),
            msprime.DiracCoalescent(psi=1234, c=56),
            msprime.SweepGenicSelection(
                position=1, start_frequency=0.5, end_frequency=0.9, s=0.1, dt=1e-4
            ),
        ]
        self.assert_repr_round_trip(examples)
