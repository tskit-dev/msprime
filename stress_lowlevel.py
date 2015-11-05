"""
Code to stress the low-level API as much as possible to expose
any memory leaks or error handling issues.
"""
from __future__ import print_function
from __future__ import division

import tempfile

import _msprime


class StressTester(object):
    """
    Class to stress the low-level interface by repeatedly running
    simulations and being as nasty as possible to the API.
    """
    def __init__(
            self, sample_size=2, num_loci=1, scaled_recombination_rate=1.0,
            random_seed=1):
        self._sample_size = sample_size
        self._num_loci = num_loci
        self._scaled_recombination_rate = scaled_recombination_rate
        self._random_seed = random_seed

    def run_module_functions(self):
        """
        Runs the module level functions.
        """
        major, minor = _msprime.get_gsl_version()
        major, minor, revision = _msprime.get_hdf5_version()

    def check_simulator_errors(self):
        """
        Provoke some errors in the simulator initialisation code
        """
        try:
            _msprime.Simulator()
            assert False
        except TypeError:
            pass
        try:
            _msprime.Simulator(sample_size=0, random_seed=1)
            assert False
        except _msprime.InputError:
            pass
        try:
            _msprime.Simulator(
                sample_size=10, random_seed=1, segment_block_size=0)
            assert False
        except _msprime.InputError:
            pass
        try:
            _msprime.Simulator(
                sample_size=10, random_seed=1, population_models=[{}])
            assert False
        except _msprime.InputError:
            pass

    def get_simulator(self):
        models = [
            {
                "type": _msprime.POP_MODEL_CONSTANT, "start_time": 0.3,
                "size": 0.2},
            {
                "type": _msprime.POP_MODEL_EXPONENTIAL, "start_time": 0.5,
                "alpha": 5}]
        sim = _msprime.Simulator(
            sample_size=self._sample_size, random_seed=self._random_seed,
            num_loci=self._num_loci,
            scaled_recombination_rate=self._scaled_recombination_rate,
            max_memory=1024**3, segment_block_size=10**6,
            coalescence_record_block_size=1000,
            population_models=models)
        return sim

    def check_uninitialised_tree_sequence(self):
        ts = _msprime.TreeSequence()
        sim = _msprime.Simulator(sample_size=5, random_seed=2)
        sim.run()
        other_ts = _msprime.TreeSequence()
        other_ts.create(sim)
        tree = _msprime.SparseTree(other_ts)
        functions = [
            ts.dump, ts.generate_mutations,
            lambda: _msprime.NewickConverter(ts),
            lambda: _msprime.TreeDiffIterator(ts),
            lambda: _msprime.SparseTreeIterator(ts, tree),
            lambda: _msprime.HaplotypeGenerator(ts)]
        for f in functions:
            try:
                f()
                assert False
            except ValueError:
                pass

    def check_tree_sequence_iterators(self, sim):
        # Check the tree diff iterator
        ts = _msprime.TreeSequence()
        ts.create(sim)
        diffs = list(_msprime.TreeDiffIterator(ts))
        iterator = _msprime.TreeDiffIterator(ts)
        del ts
        assert list(iterator) == diffs
        # Check the newick generator
        ts = _msprime.TreeSequence()
        ts.create(sim)
        trees = list(_msprime.NewickConverter(ts, 4))
        iterator = _msprime.NewickConverter(ts, 4)
        del ts
        assert list(iterator) == trees

        # Check the SparseTreeIterator
        ts = _msprime.TreeSequence()
        ts.create(sim)
        tree = _msprime.SparseTree(ts)
        parents_before = []
        for t in _msprime.SparseTreeIterator(ts, tree):
            pi = {}
            for j in range(tree.get_num_nodes()):
                pi[j] = tree.get_parent(j)
            parents_before.append(pi)
        iterator = _msprime.SparseTreeIterator(ts, tree)
        del tree
        del ts
        parents_after = []
        for t in iterator:
            pi = {}
            for j in range(t.get_num_nodes()):
                pi[j] = t.get_parent(j)
            parents_after.append(pi)
        assert parents_before == parents_after

        # Check the LeafListIterator
        ts = _msprime.TreeSequence()
        ts.create(sim)
        tree = _msprime.SparseTree(ts)
        tree_iterator = _msprime.SparseTreeIterator(ts, tree)
        del tree
        del ts
        count = 0
        for t in tree_iterator:
            count += 1
            leaf_iterator = _msprime.LeafListIterator(t, 1)
            del t
            assert list(leaf_iterator) == [1]
        assert count > 0

    def check_tree_sequence_file_errors(self, sim):
        ts = _msprime.TreeSequence()
        ts.create(sim)
        files = ["/", "/nofile", "/" + "x" * 4192]
        for f in files:
            try:
                ts.dump(f)
                assert False
            except _msprime.LibraryError:
                pass
        files = ["/", "/nofile", "/etc/passwd", "/dev/urandom"]
        for f in files:
            try:
                ts = _msprime.TreeSequence()
                ts.load(f)
                assert False
            except _msprime.LibraryError:
                pass

    def check_tree_sequence_dump_equality(self, sim):
        ts1 = _msprime.TreeSequence()
        ts1.create(sim)
        ts2 = _msprime.TreeSequence()
        ts2.create(sim)
        ts2.generate_mutations(10, 1)
        for ts in [ts1, ts2]:
            with tempfile.NamedTemporaryFile() as f:
                ts.dump(f.name)
                ts2 = _msprime.TreeSequence()
                ts2.load(f.name)
                records1 = [
                    ts.get_record(j) for j in range(ts.get_num_records())]
                records2 = [
                    ts2.get_record(j) for j in range(ts2.get_num_records())]
                assert records1 == records2
                assert ts.get_mutations() == ts2.get_mutations()

    def check_haplotype_generator(self, sim):
        ts = _msprime.TreeSequence()
        ts.create(sim)
        hg = _msprime.HaplotypeGenerator(ts)
        try:
            hg.get_haplotype(0)
            assert False
        except _msprime.LibraryError:
            pass
        n = sim.get_sample_size()
        h1 = [hg.get_haplotype(j) for j in range(1, n + 1)]
        ts.generate_mutations(10, 1)
        hg2 = _msprime.HaplotypeGenerator(ts)
        h2 = [hg2.get_haplotype(j) for j in range(1, n + 1)]
        assert h1 != h2

    def check_provenance_strings(self, sim):
        ts = _msprime.TreeSequence()
        ts.create(sim)
        s = ts.get_simulation_parameters()
        assert s is not None
        s = ts.get_mutation_parameters()
        assert s is None
        ts.generate_mutations(10, 1)
        s = ts.get_mutation_parameters()
        assert s is not None

    def check_sparse_tree(self, sim):
        try:
            _msprime.SparseTree()
            assert False
        except TypeError:
            pass
        ts = _msprime.TreeSequence()
        try:
            _msprime.SparseTree(ts)
            assert False
        except ValueError:
            pass
        ts.create(sim)
        _msprime.SparseTree(ts)

    def check_multiple_tree_sequences(self):
        params = {"sample_size": 10, "random_seed": 1}
        sim1 = _msprime.Simulator(**params)
        sim2 = _msprime.Simulator(**params)
        sim1.run()
        sim2.run()
        t1 = _msprime.TreeSequence()
        t2 = _msprime.TreeSequence()
        t1.create(sim1)
        t2.create(sim1)
        assert t1.get_num_records() == t2.get_num_records()
        r1 = [t1.get_record(j) for j in range(t1.get_num_records())]
        r2 = [t2.get_record(j) for j in range(t2.get_num_records())]
        assert r1 == r2

    def check_mutations(self, sim):
        ts = _msprime.TreeSequence()
        ts.create(sim)
        ts.generate_mutations(10, 1)
        st = _msprime.SparseTree(ts)
        all_mutations = ts.get_mutations()
        all_tree_mutations = []
        for st in _msprime.SparseTreeIterator(ts, st):
            tree_mutations = st.get_mutations()
            all_tree_mutations.extend(tree_mutations)
        assert all_tree_mutations == all_mutations
        for mutations in [[], [(0, 1), (0, 1)]]:
            ts.set_mutations(mutations)
            assert ts.get_mutations() == mutations
            params = None if len(mutations) == 0 else "{}"
            assert ts.get_mutation_parameters() == params
            with tempfile.NamedTemporaryFile() as f:
                ts.dump(f.name)

    def check_count_leaves(self, sim):
        ts = _msprime.TreeSequence()
        ts.create(sim)
        st = _msprime.SparseTree(ts)
        for st in _msprime.SparseTreeIterator(ts, st):
            for j in range(0, ts.get_num_nodes() + 1):
                assert st.get_num_tracked_leaves(j) == 0
        # provoke some error conditions to make sure we're not leaking
        # memory.
        leaves = list(range(1, 10**6))
        try:
            st = _msprime.SparseTree(ts, leaves)
        except _msprime.LibraryError:
            pass
        leaves[-1] = {}
        try:
            st = _msprime.SparseTree(ts, leaves)
        except TypeError:
            pass
        leaves = list(range(1, ts.get_sample_size() + 1))
        st = _msprime.SparseTree(ts, leaves)
        for st in _msprime.SparseTreeIterator(ts, st):
            for j in range(0, ts.get_num_nodes() + 1):
                assert (
                    st.get_num_tracked_leaves(j) == st.get_num_leaves(j))

    def run(self):
        self.run_module_functions()
        self.check_simulator_errors()
        sim = self.get_simulator()
        for t in [0.001, 0.1, 0.15, 0.5, 1000]:
            sim.run(t)
            segs = sim.get_ancestors()
            assert len(segs) == sim.get_num_ancestors()
            crs = list(sim.get_coalescence_records())
            assert len(crs) == sim.get_num_coalescence_records()
            bps = list(sim.get_breakpoints())
            assert len(bps) == sim.get_num_breakpoints()
        self.check_uninitialised_tree_sequence()
        self.check_sparse_tree(sim)
        self.check_tree_sequence_iterators(sim)
        self.check_tree_sequence_file_errors(sim)
        self.check_tree_sequence_dump_equality(sim)
        self.check_haplotype_generator(sim)
        self.check_provenance_strings(sim)
        self.check_multiple_tree_sequences()
        self.check_mutations(sim)
        self.check_count_leaves(sim)


def main():
    tester = StressTester(
        sample_size=10, num_loci=1000, scaled_recombination_rate=0.1)
    while True:
        tester.run()


if __name__ == "__main__":
    main()
