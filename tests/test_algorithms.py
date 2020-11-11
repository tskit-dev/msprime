"""
Tests for the algorithms.py script.
"""
import pathlib
import platform
import tempfile
import unittest.mock

import numpy as np
import pytest
import tskit

import algorithms
import msprime
import tests.test_cli as test_cli


IS_WINDOWS = platform.system() == "Windows"


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


@pytest.mark.skipif(IS_WINDOWS, reason="Bintrees isn't availble on windows")
class TestAlgorithms:
    def run_script(self, cmd):
        # print("RUN", cmd)
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile = pathlib.Path(tmpdir) / "out.trees"
            # Avoid dairquiri mucking up the logging setup for unittest.
            with unittest.mock.patch("daiquiri.setup"):
                algorithms.main(cmd.split() + [str(outfile)])
            return tskit.load(outfile)

    def test_defaults(self):
        ts = self.run_script("10")
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert not has_discrete_genome(ts)
        assert ts.sequence_length == 100

    def test_discrete(self):
        ts = self.run_script("10 -d")
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)

    def test_dtwf(self):
        ts = self.run_script("10 --model=dtwf")
        assert ts.num_trees > 1
        assert not has_discrete_genome(ts)
        assert ts.sequence_length == 100

    def test_dtwf_migration(self):
        ts = self.run_script("10 -r 0 --model=dtwf -p 2 -g 0.1")
        assert ts.num_trees == 1
        assert ts.sequence_length == 100
        assert ts.num_populations == 2

    def test_dtwf_discrete(self):
        ts = self.run_script("10 -d --model=dtwf")
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)

    def test_full_arg(self):
        ts = self.run_script("30 -L 200 --full-arg")
        assert ts.num_trees > 1
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_RE_EVENT) > 0
        assert np.sum(node_flags == msprime.NODE_IS_CA_EVENT) > 0

    def test_migration_full_arg(self):
        ts = self.run_script("10 -p 3 -g 0.1 --full-arg")
        assert ts.num_trees > 1
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_MIG_EVENT) > 0

    def test_gc(self):
        ts = self.run_script("10 -c 0.4 2 -d")
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)
        assert ts.sequence_length == 100

    def test_recomb_map(self):
        # Need the final -- to break the end of the recomb_rates option
        ts = self.run_script(
            "10 -d --recomb-positions 0 10 100 --recomb-rates 0 0.1 0 --"
        )
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)
        assert ts.sequence_length == 100

    def test_census_event(self):
        ts = self.run_script("10 --census-time 0.01")
        assert ts.num_trees > 1
        node_time = ts.tables.nodes.time
        assert np.sum(node_time == 0.01) > 0

    def test_single_sweep(self):
        ts = self.run_script(
            "10 --model=single_sweep --trajectory 0.1 0.9 0.01 --time-slice=0.1"
        )
        assert ts.num_trees > 1

    def test_single_sweep_shorter(self):
        # Try to catch some situtations not covered in other test.
        ts = self.run_script(
            "10 --model=single_sweep --trajectory 0.1 0.5 0.1 --time-slice=0.01"
        )
        assert ts.num_trees > 1

    def test_bottleneck(self):
        ts = self.run_script("10 -r 0 --bottleneck 0.1 0 2")
        assert ts.num_trees == 1
        tree = ts.first()
        found = False
        for node in tree.nodes():
            if len(tree.children(node)) > 2:
                found = True
        assert found

    def test_pop_size_change_event(self):
        ts = self.run_script("10 -r 0 --population-size-change 1 0 1")
        assert ts.num_trees == 1

    def test_pop_growth_rate_change_event(self):
        ts = self.run_script("10 -r 0 --population-growth-rate-change 1 0 0.1")
        assert ts.num_trees == 1

    def test_migration_rate_change_event(self):
        ts = self.run_script("10 -r 0 -p 2 --migration-matrix-element-change 1 0 1 0.1")
        assert ts.num_trees == 1

    def test_verbose(self):
        stdout, stderr = test_cli.capture_output(
            self.run_script, "10 -c 0.4 2 -v --end-time=0.1"
        )
        assert len(stdout) > 0
        assert len(stderr) == 0

    def test_from_ts(self):
        from_ts = msprime.sim_ancestry(
            10,
            sequence_length=100,
            discrete_genome=True,
            recombination_rate=2,
            random_seed=2,
            end_time=0.1,
        )
        tree = from_ts.first()
        assert tree.num_roots > 1
        with tempfile.TemporaryDirectory() as tmpdir:
            from_ts_path = pathlib.Path(tmpdir) / "from.trees"
            from_ts.dump(from_ts_path)
            ts = self.run_script(f"0 --from-ts {from_ts_path}")
            assert ts.num_trees > 1
            assert ts.sequence_length == 100

    def test_pedigree(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ped_path = pathlib.Path(tmpdir) / "tmp.ped"
            individual = np.array([1, 2, 3, 4], dtype=int)
            parents = np.array([2, 3, 2, 3, -1, -1, -1, -1], dtype=int).reshape(-1, 2)
            time = np.array([0, 0, 1, 1])
            is_sample = np.array([1, 1, 0, 0], dtype=int)
            ped = msprime.Pedigree(
                individual, parents, time, is_sample, sex=None, ploidy=2
            )
            ped.save_txt(ped_path)
            ts = self.run_script(f"2 --pedigree-file {ped_path} --model=wf_ped")
            for tree in ts.trees():
                assert tree.num_roots > 1
            assert ts.num_individuals == len(individual)
            assert ts.num_samples == 4
