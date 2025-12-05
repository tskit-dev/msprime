# Tests for the algorithms.py script.
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
from tests.test_pedigree import simulate_pedigree


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


def verify_unary(ts):
    min_children = np.zeros(ts.num_nodes, dtype=int)
    max_children = np.zeros_like(min_children)

    for tree in ts.trees():
        for i in range(ts.num_nodes):
            n = tree.num_children_array[i]
            if n > 0:
                if min_children[i] > 0:
                    min_children[i] = min(min_children[i], n)
                else:
                    min_children[i] = n
                max_children[i] = max(max_children[i], n)

    assert np.any(min_children == 1)
    for minc, maxc in zip(min_children, max_children):
        if minc == 1:
            assert maxc >= 2


def verify_pedigree_unary(ts):
    direct_ancestors = [set() for _ in range(ts.num_nodes)]

    # collect all direct ancestors for each node
    for node in ts.nodes():
        individual = ts.individual(node.individual)
        for parent_id in individual.parents:
            if parent_id != tskit.NULL:
                parent = ts.individual(parent_id)
                for parent_node in parent.nodes:
                    direct_ancestors[node.id].add(parent_node)

    # assert that in each marginal tree all nodes a parent from
    # the constructed set
    for tree in ts.trees():
        for node_id in tree.nodes():
            if tree.parent(node_id) != tskit.NULL:
                direct_ancestor = direct_ancestors[node_id]
                assert tree.parent(node_id) in direct_ancestor


def verify_dtwf_unary(ts):
    # assert all nodes in each tree and subtree only differ one generation
    for tree in ts.trees():
        queue = [tree.root]
        while queue:
            parent = queue.pop(0)
            time = tree.time(parent) - 1
            num_children = 0
            for child in tree.children(parent):
                assert tree.time(child) == time
                queue.append(child)
                num_children += 1
            if num_children == 0:
                assert tree.time(parent) == 0


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

    @pytest.mark.parametrize("options", ["", "--continue-after-local-mrca"])
    def test_dtwf(self, options):
        ts = self.run_script(f"10 --model=dtwf {options}")
        assert ts.num_trees > 1
        assert not has_discrete_genome(ts)
        assert ts.sequence_length == 100

    @pytest.mark.parametrize("options", ["", "--continue-after-local-mrca"])
    def test_dtwf_migration(self, options):
        ts = self.run_script(f"10 -r 0 --model=dtwf -p 2 -g 0.1 {options}")
        assert ts.num_trees == 1
        assert ts.sequence_length == 100
        assert ts.num_populations == 2

    def test_dtwf_discrete(self):
        ts = self.run_script("10 -d --model=dtwf")
        assert ts.num_trees > 1
        assert has_discrete_genome(ts)

    @pytest.mark.parametrize("options", ["", "--continue-after-local-mrca"])
    def test_full_arg(self, options):
        node_value = sum(2**i for i in (17, 18, 19, 21))
        ts = self.run_script(
            f"30 -L 200 --additional-nodes {node_value} --all-segments {options}"
        )
        assert ts.num_trees > 1
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_RE_EVENT) > 0
        assert np.sum(node_flags == msprime.NODE_IS_CA_EVENT) > 0

    def test_migration_full_arg(self):
        node_value = sum(2**i for i in (17, 18, 19, 21))
        ts = self.run_script(
            f"10 -p 3 -g 0.1 --additional-nodes {node_value} --all-segments"
        )
        assert ts.num_trees > 1
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_MIG_EVENT) > 0

    def test_store_unary(self):
        ts = self.run_script("10 --all-segments")
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert not has_discrete_genome(ts)
        assert ts.sequence_length == 100
        verify_unary(ts)

    @pytest.mark.parametrize("options", ["", "--continue-after-local-mrca"])
    def test_store_unary_dtwf(self, options):
        node_value = 1 << 18 | 1 << 22
        cmd = (
            f"10 --all-segments -d --additional-nodes {node_value} "
            f"--model=dtwf {options}"
        )
        ts = self.run_script(cmd)
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert ts.sequence_length == 100
        verify_dtwf_unary(ts)

    @pytest.mark.parametrize("options", ["", "--continue-after-local-mrca"])
    def test_store_unary_dtwf_re(self, options):
        node_value = 1 << 17 | 1 << 18 | 1 << 22
        cmd = (
            f"10 --all-segments -d --additional-nodes {node_value} "
            f"--model=dtwf {options}"
        )
        ts = self.run_script(cmd)
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert ts.sequence_length == 100
        assert (
            np.sum(np.bitwise_and(ts.tables.nodes.flags, msprime.NODE_IS_RE_EVENT) > 0)
            > 0
        )
        verify_dtwf_unary(ts)

    def test_store_unary_pedigree(self):
        tables = simulate_pedigree(num_founders=4, num_generations=10)
        node_value = 1 << 18 | 1 << 22
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            tables.dump(ts_path)
            ts = self.run_script(
                f"0 --from-ts {ts_path} -r 1 --all-segments  "
                f"--additional-nodes {node_value} --model=fixed_pedigree"
            )
        verify_pedigree_unary(ts)

    def test_store_unary_pedigree_re(self):
        tables = simulate_pedigree(num_founders=4, num_generations=10)
        node_value = 1 << 17 | 1 << 18 | 1 << 22
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            tables.dump(ts_path)
            ts = self.run_script(
                f"0 --from-ts {ts_path} -r 1 --all-segments  "
                f"--additional-nodes {node_value} --model=fixed_pedigree"
            )
        assert (
            np.sum(np.bitwise_and(ts.tables.nodes.flags, msprime.NODE_IS_RE_EVENT) > 0)
            > 0
        )
        verify_pedigree_unary(ts)

    def test_store_unary_pedigree_small(self):
        pb = msprime.PedigreeBuilder()
        mom_id = pb.add_individual(time=1)
        dad_id = pb.add_individual(time=1)
        pb.add_individual(time=0, parents=[mom_id, dad_id], is_sample=True)
        pedigree = pb.finalise(sequence_length=100)
        node_value = 1 << 18 | 1 << 22
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            pedigree.dump(ts_path)
            ts = self.run_script(
                f"0 --from-ts {ts_path} -r 0.1 --all-segments  "
                f"--additional-nodes {node_value} --model=fixed_pedigree "
            )
        verify_pedigree_unary(ts)

    def test_store_common_ancestor(self):
        node_value = 1 << 18
        ts = self.run_script(
            f"10 -r 0.5 --additional-nodes {node_value} --all-segments"
        )
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert ts.sequence_length == 100
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_CA_EVENT) > 0

    def test_store_recombinant(self):
        node_value = 1 << 17
        ts = self.run_script(f"10 --additional-nodes {node_value} --all-segments")
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert ts.sequence_length == 100
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_RE_EVENT) > 0

    def test_store_migrant(self):
        node_value = 1 << 19
        ts = self.run_script(
            f"10 -p 3 -g 0.1 --additional-nodes {node_value} --all-segments"
        )
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert ts.sequence_length == 100
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_MIG_EVENT) > 0

    def test_store_ca_re(self):
        node_value = 1 << 17 | 1 << 18
        ts = self.run_script(
            f"10 -p 3 -g 0.1 -r 0.1 --additional-nodes {node_value} --all-segments"
        )
        assert ts.num_samples == 10
        assert ts.num_trees > 1
        assert ts.sequence_length == 100
        node_flags = ts.tables.nodes.flags
        assert np.sum(node_flags == msprime.NODE_IS_CA_EVENT) > 0
        assert np.sum(node_flags == msprime.NODE_IS_RE_EVENT) > 0

    @pytest.mark.parametrize("model", ["hudson", "smc_k", "dtwf"])
    def test_continue_after_local_mrca_multi_tree(self, model):
        r = 0.1
        ts = self.run_script(f"10 --continue-after-local-mrca -r {r} --model={model}")
        assert ts.num_trees >= 2
        # All roots should have the same time
        root_times = []
        for tree in ts.trees():
            root_times.append(tree.time(tree.root))
        assert len(set(root_times)) == 1

    @pytest.mark.parametrize("model", ["hudson", "smc_k", "dtwf"])
    def test_continue_after_local_mrca_single_tree(self, model):
        ts = self.run_script(f"10 --continue-after-local-mrca -r 0 --model={model}")
        root_times = [tree.time(tree.root) for tree in ts.trees()]
        assert len(set(root_times)) == 1

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

    @pytest.mark.parametrize(
        ["num_founders", "num_generations", "r"],
        [
            (4, 5, 0),
            (4, 5, 0.1),
            (4, 5, 1),
            (4, 10, 0.1),
            (10, 5, 0.1),
            (100, 2, 0.1),
        ],
    )
    def test_pedigree(self, num_founders, num_generations, r):
        tables = simulate_pedigree(
            num_founders=num_founders, num_generations=num_generations
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            tables.dump(ts_path)
            ts = self.run_script(f"0 --from-ts {ts_path} -r {r} --model=fixed_pedigree")
        for node in ts.nodes():
            assert node.individual != tskit.NULL
        assert ts.num_samples == 2 * num_founders

        # Make the set of ancestors for each individual
        ancestors = [set() for _ in ts.individuals()]
        for individual in ts.individuals():
            for parent_id in individual.parents:
                stack = [parent_id]
                while len(stack) > 0:
                    parent_id = stack.pop()
                    if parent_id != tskit.NULL:
                        ancestors[individual.id].add(parent_id)
                        for u in ts.individual(parent_id).parents:
                            stack.append(u)
        # founder nodes have correct time
        founder_ids = {
            ind.id for ind in ts.individuals() if len(ancestors[ind.id]) == 0
        }
        for founder in founder_ids:
            for node_id in ts.individual(founder).nodes:
                node = ts.node(node_id)
                assert node.time == num_generations - 1
        for tree in ts.trees():
            for root in tree.roots:
                node = ts.node(root)
                # If this is a unary root it must be from a founder.
                if tree.num_children(root) == 1:
                    assert node.individual in founder_ids
            for node_id in tree.nodes():
                node = ts.node(node_id)
                individual = ts.individual(node.individual)
                if tree.parent(node_id) != tskit.NULL:
                    parent_node = ts.node(tree.parent(node_id))
                    assert parent_node.individual in ancestors[individual.id] | {
                        tskit.NULL
                    }
                else:
                    assert node.individual in founder_ids

    @pytest.mark.parametrize("r", [0, 0.1, 1])
    def test_pedigree_trio(self, r):
        input_tables = simulate_pedigree(
            num_founders=2, num_children_prob=(0, 1), num_generations=2
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            input_tables.dump(ts_path)
            ts = self.run_script(f"0 --from-ts {ts_path} -r {r} --model=fixed_pedigree")
        output_tables = ts.dump_tables()
        input_tables.individuals.assert_equals(output_tables.individuals)
        input_tables.nodes.assert_equals(output_tables.nodes[: len(input_tables.nodes)])
        assert len(output_tables.edges) >= 2

    @pytest.mark.parametrize("r", [0, 0.1, 1])
    def test_pedigree_internal_sample(self, r):
        input_tables = simulate_pedigree(num_founders=4, num_generations=4)
        flags = input_tables.nodes.flags
        flags[8:12] = tskit.NODE_IS_SAMPLE  # Make inds 4 and 5 samples
        input_tables.nodes.flags = flags
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            input_tables.dump(ts_path)
            ts = self.run_script(f"0 --from-ts {ts_path} -r {r} --model=fixed_pedigree")
        output_tables = ts.dump_tables()
        input_tables.individuals.assert_equals(output_tables.individuals)
        input_tables.nodes.assert_equals(output_tables.nodes[: len(input_tables.nodes)])
        assert len(output_tables.edges) >= 2

    @pytest.mark.parametrize("r", [0, 0.1, 1])
    def test_pedigree_only_internal_samples(self, r):
        input_tables = simulate_pedigree(
            num_founders=4, num_generations=4, sample_gen=[1]
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            input_tables.dump(ts_path)
            ts = self.run_script(f"0 --from-ts {ts_path} -r {r} --model=fixed_pedigree")
        output_tables = ts.dump_tables()
        input_tables.individuals.assert_equals(output_tables.individuals)
        input_tables.nodes.assert_equals(output_tables.nodes[: len(input_tables.nodes)])
        assert np.all(ts.nodes_time[ts.samples()] == 1)
        assert len(output_tables.edges) >= 2

    @pytest.mark.parametrize("num_founders", [1, 2, 20])
    def test_one_gen_pedigree(self, num_founders):
        tables = simulate_pedigree(num_founders=num_founders, num_generations=1)
        with tempfile.TemporaryDirectory() as tmpdir:
            ts_path = pathlib.Path(tmpdir) / "pedigree.trees"
            tables.dump(ts_path)
            ts = self.run_script(f"0 --from-ts {ts_path} -r 1 --model=fixed_pedigree")
        assert len(ts.dump_tables().edges) == 0

    @pytest.mark.parametrize(
        "cmd",
        [
            "10 -L 1000 -d -r 0.01 --model smc_k",
            "10 -L 1000 -r 0.01 --model smc_k",
            "10 -L 1000 -r 0.01 --model smc_k --offset 0.50",
            "10 -L 1000 -d -r 0.01 --model smc_k -p 2 -g 0.1",
            "10 -L 1000 -d -c 0.04 2  --model smc_k",
            "10 -L 1000 -c 0.04 2  --model smc_k --offset 0.75",
        ],
    )
    def test_smck(self, cmd):
        ts = self.run_script(cmd)
        assert ts.num_trees > 1
        for tree in ts.trees():
            assert tree.num_roots == 1
