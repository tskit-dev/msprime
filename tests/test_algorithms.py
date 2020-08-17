"""
Tests for the algorithms.py script.
"""
import pathlib
import platform
import tempfile
import unittest

import numpy as np
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


@unittest.skipIf(IS_WINDOWS, "Bintrees isn't availble on windows")
class TestAlgorithms(unittest.TestCase):
    def run_script(self, cmd):
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile = pathlib.Path(tmpdir) / "out.trees"
            algorithms.main(cmd.split() + [str(outfile)])
            return tskit.load(outfile)

    def test_defaults(self):
        ts = self.run_script("10")
        self.assertEqual(ts.num_samples, 10)
        self.assertGreater(ts.num_trees, 1)
        self.assertFalse(has_discrete_genome(ts))
        self.assertEqual(ts.sequence_length, 100)

    def test_discrete(self):
        ts = self.run_script("10 -d")
        self.assertGreater(ts.num_trees, 1)
        self.assertTrue(has_discrete_genome(ts))

    def test_dtwf(self):
        ts = self.run_script("10 --model=dtwf")
        self.assertGreater(ts.num_trees, 1)
        self.assertFalse(has_discrete_genome(ts))
        self.assertEqual(ts.sequence_length, 100)

    def test_dtwf_discrete(self):
        ts = self.run_script("10 -d --model=dtwf")
        self.assertGreater(ts.num_trees, 1)
        self.assertTrue(has_discrete_genome(ts))

    def test_full_arg(self):
        ts = self.run_script("30 -L 200 --full-arg")
        self.assertGreater(ts.num_trees, 1)
        node_flags = ts.tables.nodes.flags
        self.assertGreater(np.sum(node_flags == msprime.NODE_IS_RE_EVENT), 0)
        self.assertGreater(np.sum(node_flags == msprime.NODE_IS_CA_EVENT), 0)

    def test_gc(self):
        ts = self.run_script("10 -c 0.1 2 -d")
        self.assertGreater(ts.num_trees, 1)
        self.assertTrue(has_discrete_genome(ts))
        self.assertEqual(ts.sequence_length, 100)

    def test_verbose(self):
        stdout, stderr = test_cli.capture_output(self.run_script, "10 -v")
        self.assertGreater(len(stdout), 0)
        self.assertEqual(len(stderr), 0)
