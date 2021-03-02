#
# Copyright (C) 2020-2021 University of Oxford
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
Tests for the parsing of species trees in newick and StarBEAST format.
"""
import collections

import numpy as np
import pytest

import msprime
import msprime.species_trees as species_trees


def get_non_binary_tree(n):
    demographic_events = [
        msprime.SimpleBottleneck(time=0.1, population=0, proportion=0.5)
    ]
    ts = msprime.simulate(n, demographic_events=demographic_events, random_seed=3)
    tree = ts.first()
    found = False
    for u in tree.nodes():
        if len(tree.children(u)) > 2:
            found = True
    assert found
    return tree


class TestIsNumber:
    """
    Test the is_number function.
    """

    def test_good_examples(self):
        for x in ["2", "2.0", "1000", "-1e3", "1e-6"]:
            assert species_trees.is_number(x)

    def test_bad_examples(self):
        for x in ["", "x2.0", "1000x", ";-1e3", ";;"]:
            assert not species_trees.is_number(x)


class TestParseNumberOrMapping:
    @pytest.mark.parametrize("N", [1, 1000, 0.01, np.array([1])[0]])
    def test_single_value(self, N):
        res = species_trees.parse_number_or_mapping(N, "")
        # This is a defaultdict, so we should get the right value regardless
        # of the input.
        for key in ["x", "", None, 123]:
            assert res[key] == N

    def test_simple_mapping(self):
        arg = {"x": 1, "y": 2}
        res = species_trees.parse_number_or_mapping(arg, "")
        assert arg is res

    def test_defaultdict(self):
        arg = collections.defaultdict(lambda x: 1234)
        res = species_trees.parse_number_or_mapping(arg, "")
        assert arg is res

    @pytest.mark.parametrize("value", [[], None])
    def test_errors(self, value):
        message = "this is a message"
        with pytest.raises(TypeError):
            species_trees.parse_number_or_mapping(value, message)


class TestSpeciesNamedInternalNodes:
    def test_initial_size_full_map(self):
        tree = "(A:10.0,B:10.0)C"
        initial_size = {"A": 234, "B": 567, "C": 8}
        demography = msprime.Demography.from_species_tree(tree, initial_size)
        assert demography.populations[0].name == "A"
        assert demography.populations[0].initial_size == 234
        assert demography.populations[1].name == "B"
        assert demography.populations[1].initial_size == 567
        assert demography.populations[2].name == "C"
        assert demography.populations[2].initial_size == 8

    def test_initial_size_partial_map(self):
        tree = "(A:10.0,B:10.0)C"
        initial_size = collections.defaultdict(lambda: 100)
        initial_size.update({"A": 234, "C": 8})
        demography = msprime.Demography.from_species_tree(tree, initial_size)
        assert demography.populations[0].name == "A"
        assert demography.populations[0].initial_size == 234
        assert demography.populations[1].name == "B"
        assert demography.populations[1].initial_size == 100
        assert demography.populations[2].name == "C"
        assert demography.populations[2].initial_size == 8

    def test_growth_rate_full_map(self):
        tree = "(A:10.0,B:10.0)C"
        growth_rate = {"A": 234, "B": 567, "C": 8}
        demography = msprime.Demography.from_species_tree(
            tree, 10, growth_rate=growth_rate
        )
        assert demography.populations[0].name == "A"
        assert demography.populations[0].initial_size == 10
        assert demography.populations[0].growth_rate == 234
        assert demography.populations[1].name == "B"
        assert demography.populations[1].growth_rate == 567
        assert demography.populations[1].initial_size == 10
        assert demography.populations[2].name == "C"
        assert demography.populations[2].growth_rate == 8
        assert demography.populations[2].initial_size == 10

    def test_growth_rate_partial_map(self):
        tree = "(A:10.0,B:10.0)C"
        growth_rate = collections.defaultdict(lambda: 100)
        growth_rate.update({"A": 234, "C": 8})
        demography = msprime.Demography.from_species_tree(
            tree, 10, growth_rate=growth_rate
        )
        assert demography.populations[0].name == "A"
        assert demography.populations[0].initial_size == 10
        assert demography.populations[0].growth_rate == 234
        assert demography.populations[1].name == "B"
        assert demography.populations[1].initial_size == 10
        assert demography.populations[1].growth_rate == 100
        assert demography.populations[2].name == "C"
        assert demography.populations[2].initial_size == 10
        assert demography.populations[2].growth_rate == 8


class TestSpeciesTreeRoundTrip:
    """
    Tests that we get what we expect when we parse trees produced from
    msprime/tskit.
    """

    def make_newick(self, tree):
        """
        Return a standard newick encoding we can use to get node IDs.
        """
        return tree.newick(node_labels={u: f"node_{u}" for u in tree.nodes()})

    def verify_non_ultrametric(self, tree):
        newick = tree.newick()
        with pytest.raises(ValueError):
            species_trees.parse_species_tree(newick, initial_size=1)

    def verify(
        self,
        tree,
        newick=None,
        initial_size=1,
        time_units="gen",
        generation_time=None,
    ):
        if newick is None:
            newick = self.make_newick(tree)
        demography = species_trees.parse_species_tree(
            newick,
            initial_size=initial_size,
            time_units=time_units,
            generation_time=generation_time,
        )
        assert demography.num_populations == tree.num_nodes
        for pop in demography.populations:
            assert pop.initial_size == initial_size
            assert pop.growth_rate == 0
            assert pop.name is not None

        # Population IDs are mapped to leaves first, and then to the internal nodes
        # in postorder
        pop_id_map = {}
        k = 0
        for u in tree.leaves():
            pop_id_map[u] = k
            k += 1

        for u in tree.nodes(order="postorder"):
            if tree.is_internal(u):
                pop_id_map[u] = k
                k += 1

        for u in tree.nodes():
            pop = demography.populations[pop_id_map[u]]
            assert pop.growth_rate == 0
            if tree.is_leaf(u):
                # Assuming we're using the make_newick function above
                assert pop.name == f"node_{u}"

        # We should have demographic events for every internal node, and
        # events should be output in increasing time order.
        j = 0
        for node in tree.nodes(order="timeasc"):
            if tree.is_internal(node):
                event = demography.events[j]
                j += 1
                assert isinstance(event, msprime.PopulationSplit)
                assert event.time == pytest.approx(tree.time(node))
                assert event.ancestral == demography[pop_id_map[node]].name
                assert event.derived == [
                    demography[pop_id_map[child]].name for child in tree.children(node)
                ]
        assert j == len(demography.events)

    def test_n2_binary(self):
        tree = msprime.simulate(2, random_seed=2).first()
        self.verify(tree)

    def test_n2_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, 0), (0, 1)], random_seed=2)
        self.verify_non_ultrametric(ts.first())

    def test_n5_binary(self):
        ts = msprime.simulate(5, random_seed=2)
        tree = ts.first()
        self.verify(tree, initial_size=1)

    def test_n5_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, j) for j in range(5)], random_seed=2)
        self.verify_non_ultrametric(ts.first())

    def test_n7_binary(self):
        ts = msprime.simulate(7, random_seed=2)
        tree = ts.first()
        self.verify(tree, initial_size=11)

    def test_n7_binary_embedded_whitespace(self):
        # Check for embedded whitespace in the newick string
        tree = msprime.simulate(7, random_seed=2).first()
        newick = self.make_newick(tree)
        self.verify(tree, newick="    " + newick)
        self.verify(tree, newick=newick + "        ")
        self.verify(tree, newick=newick + "\n")
        self.verify(tree, newick=newick.replace("(", "( "))
        self.verify(tree, newick=newick.replace("(", "(\n"))
        self.verify(tree, newick=newick.replace(")", ") "))
        self.verify(tree, newick=newick.replace(")", ")\n"))
        self.verify(tree, newick=newick.replace(":", " : "))
        self.verify(tree, newick=newick.replace(":", "\n:\n"))
        self.verify(tree, newick=newick.replace(":", "\t:\t"))
        self.verify(tree, newick=newick.replace(",", "\t,\t"))
        self.verify(tree, newick=newick.replace(",", "\n,\n"))
        self.verify(tree, newick=newick.replace(",", "     ,  "))

    def test_n100_binary(self):
        ts = msprime.simulate(100, random_seed=2)
        tree = ts.first()
        self.verify(tree, initial_size=11)

    def test_n10_non_binary(self):
        tree = get_non_binary_tree(10)
        self.verify(tree, initial_size=3.1234)

    def test_n10_binary_years(self):
        ts = msprime.simulate(10, random_seed=2)
        generation_time = 5
        tree = ts.first()
        self.verify(tree, initial_size=1, time_units="yr", generation_time=1)
        tables = ts.dump_tables()
        times = tables.nodes.time
        flags = tables.nodes.flags
        scaled_times = [generation_time * time for time in times]
        tables.nodes.set_columns(flags=flags, time=scaled_times)
        ts = tables.tree_sequence()
        scaled_tree = ts.first()
        self.verify(
            tree,
            newick=self.make_newick(scaled_tree),
            time_units="yr",
            generation_time=generation_time,
        )

    def test_n10_binary_million_years(self):
        ts = msprime.simulate(10, random_seed=2)
        generation_time = 5
        tree = ts.first()
        tables = ts.dump_tables()
        times = tables.nodes.time
        flags = tables.nodes.flags
        scaled_times = [time / (1e6 / generation_time) for time in times]
        tables.nodes.set_columns(flags=flags, time=scaled_times)
        ts = tables.tree_sequence()
        scaled_tree = ts.first()
        self.verify(
            tree,
            newick=self.make_newick(scaled_tree),
            time_units="myr",
            generation_time=generation_time,
        )


def make_nexus(tree, pop_size_map):
    """
    Returns the specified tree formatted as StarBEAST compatible nexus.
    """
    node_labels = {}
    leaf_names = []
    count = 0
    for u in tree.nodes():
        name = ""
        if tree.is_leaf(u):
            count += 1
            name = str(u)
            leaf_names.append(name)
            node_labels[u] = f"{count}[&dmv={{{pop_size_map[u]}}},"
            node_labels[u] += "dmv1=0.260,dmv1_95%_HPD={0.003,0.625},"
            node_labels[u] += "dmv1_median=0.216,dmv1_range={0.001,1.336},"
            node_labels[u] += "height=1.310E-15,height_95%_HPD={0.0,3.552E-15},"
            node_labels[u] += "height_median=0.0,height_range={0.0,7.105E-15},"
            node_labels[u] += "length=2.188,length_95%_HPD={1.725,2.634},"
            node_labels[u] += "length_median=2.182,length_range={1.307,3.236}]"
        else:
            node_labels[u] = f"[&dmv={{{pop_size_map[u]}}},"
            node_labels[u] += "dmv1=0.260,dmv1_95%_HPD={0.003,0.625},"
            node_labels[u] += "dmv1_median=0.216,dmv1_range={0.001,1.336},"
            node_labels[u] += "height=1.310E-15,height_95%_HPD={0.0,3.552E-15},"
            node_labels[u] += "height_median=0.0,height_range={0.0,7.105E-15},"
            node_labels[u] += "length=2.188,length_95%_HPD={1.725,2.634},"
            node_labels[u] += "length_median=2.182,length_range={1.307,3.236}]"
    newick = tree.newick(node_labels=node_labels)
    out = "#NEXUS\n\n"
    out += "Begin taxa;\n"
    out += "    Dimensions ntax=" + str(len(leaf_names)) + ";\n"
    out += "    Taxlabels\n"
    for name in leaf_names:
        out += "        spc" + str(name) + "\n"
    out += "        ;\n"
    out += "End;\n"
    out += "Begin trees;\n"
    out += "    Translate\n"
    count = 0
    for name in leaf_names:
        count += 1
        out += "             " + str(count) + " spc" + name + ",\n"
    out = out[:-2]
    out += "\n;\n"
    out += "tree TREE1 = " + newick + "\n"
    out += "End;\n"
    return out


class TestStarbeastRoundTrip:
    """
    Tests that we get what we expect when we parse trees produced from
    msprime/tskit.
    """

    def verify_non_ultrametric(self, tree, pop_size_map):
        nexus = make_nexus(tree, pop_size_map)
        with pytest.raises(ValueError):
            species_trees.parse_starbeast(nexus, 10)

    def verify(
        self,
        tree,
        pop_size_map,
        nexus=None,
        time_units="yr",
        generation_time=1,
    ):
        if nexus is None:
            nexus = make_nexus(tree, pop_size_map)
        demography = species_trees.parse_starbeast(nexus, generation_time, time_units)
        assert demography.num_populations == tree.num_nodes
        for pop in demography.populations:
            assert pop.growth_rate == 0

        # Population IDs are mapped to leaves first, and then to the internal nodes
        # in postorder
        pop_id_map = {}
        k = 0
        for u in tree.leaves():
            pop_id_map[u] = k
            k += 1

        for u in tree.nodes(order="postorder"):
            if tree.is_internal(u):
                pop_id_map[u] = k
                k += 1

        for u in tree.nodes():
            pop = demography.populations[pop_id_map[u]]
            assert pop.growth_rate == 0
            assert pop.initial_size == pop_size_map[u]
            if tree.is_leaf(u):
                # Note: we're assuming the default newick here in tskit that labels
                # nodes as their id + 1.
                assert pop.name == f"spc{u}"
            else:
                assert pop.name == f"pop_{pop_id_map[u]}"

        # We should have demographic events for every internal node, and
        # events should be output in increasing time order.
        j = 0
        for node in tree.nodes(order="timeasc"):
            if tree.is_internal(node):
                event = demography.events[j]
                j += 1
                assert isinstance(event, msprime.PopulationSplit)
                assert event.time == pytest.approx(tree.time(node))
                assert event.ancestral == demography[pop_id_map[node]].name
                assert event.derived == [
                    demography[pop_id_map[child]].name for child in tree.children(node)
                ]
        assert j == len(demography.events)

    def test_n2_binary(self):
        tree = msprime.simulate(2, random_seed=2).first()
        self.verify(tree, {u: 1 for u in tree.nodes()})

    def test_n2_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, 0), (0, 1)], random_seed=2)
        tree = ts.first()
        self.verify_non_ultrametric(tree, {u: 2.123 for u in tree.nodes()})

    def test_n5_binary(self):
        ts = msprime.simulate(5, random_seed=2)
        tree = ts.first()
        self.verify(tree, {u: 1 + u for u in tree.nodes()})

    def test_n5_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, j) for j in range(5)], random_seed=2)
        tree = ts.first()
        self.verify_non_ultrametric(tree, {u: 1 / (1 + u) for u in tree.nodes()})

    def test_n7_binary(self):
        ts = msprime.simulate(7, random_seed=2)
        tree = ts.first()
        self.verify(tree, {u: 7 for u in tree.nodes()})

    def test_n100_binary(self):
        ts = msprime.simulate(100, random_seed=2)
        tree = ts.first()
        self.verify(tree, {u: 1e-4 for u in tree.nodes()})

    def test_n10_non_binary(self):
        tree = get_non_binary_tree(10)
        self.verify(tree, {u: 0.1 for u in tree.nodes()})

    def test_n10_binary_million_years(self):
        ts = msprime.simulate(10, random_seed=2)
        generation_time = 5
        tree = ts.first()
        pop_size_map = {u: 0.1 for u in tree.nodes()}
        nexus = make_nexus(tree, pop_size_map)
        tables = ts.dump_tables()
        times = tables.nodes.time
        flags = tables.nodes.flags
        scaled_times = [time * (1e6 / generation_time) for time in times]
        tables.nodes.set_columns(flags=flags, time=scaled_times)
        ts = tables.tree_sequence()
        scaled_tree = ts.first()
        scaled_pop_size_map = {u: 0.1 * (1e6 / generation_time) for u in pop_size_map}
        self.verify(
            scaled_tree,
            nexus=nexus,
            pop_size_map=scaled_pop_size_map,
            time_units="myr",
            generation_time=generation_time,
        )


class TestSpeciesTreeParsingErrors:
    """
    Tests for parsing of species trees in newick format.
    """

    def test_bad_params(self):
        with pytest.raises(TypeError):
            species_trees.parse_species_tree()
        with pytest.raises(TypeError):
            species_trees.parse_species_tree(tree="()")
        with pytest.raises(TypeError):
            species_trees.parse_species_tree(initial_size=1)

    def test_unequal_branch_lengths(self):
        with pytest.raises(ValueError):
            species_trees.parse_species_tree(
                tree="(popA:100.0,popB:10.0)", initial_size=1000
            )

    def test_duplicate_name(self):
        with pytest.raises(ValueError, match="Duplicate population name"):
            species_trees.parse_species_tree(
                tree="(popA:100.0,popA:100.0)", initial_size=1
            )

    def test_bad_tree(self):
        bad_trees = [
            "",
            ";",
            "abcd",
            ";;;",
            "___",
            "∞",
            "(",
            ")",
            "()",
            "( )",
            "(()())",
            "((3:0.39,5:0.39]:1.39,(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
            "((3:0.39,5:0.39(:1.39,(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
            "((3:0.39,5:0.39,:1.39,(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
            "(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
        ]
        for bad_tree in bad_trees:
            with pytest.raises(ValueError):
                species_trees.parse_species_tree(tree=bad_tree, initial_size=1)

    def test_bad_parameter(self):
        good_tree = "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)"
        good_time_units = "myr"
        good_ne = 10000
        good_generation_time = 5
        for bad_time_units in [-3, "asdf", ["myr"]]:
            with pytest.raises(ValueError):
                species_trees.parse_species_tree(
                    good_tree,
                    time_units=bad_time_units,
                    initial_size=good_ne,
                    generation_time=good_generation_time,
                )

        with pytest.raises(TypeError):
            species_trees.parse_species_tree(good_tree, None)

        for bad_ne in [-3, "x"]:
            with pytest.raises(ValueError):
                species_trees.parse_species_tree(
                    good_tree,
                    time_units=good_time_units,
                    initial_size=bad_ne,
                    generation_time=good_generation_time,
                )
        for bad_generation_time in [None, -3, "x"]:
            with pytest.raises(ValueError):
                species_trees.parse_species_tree(
                    good_tree,
                    time_units=good_time_units,
                    initial_size=good_ne,
                    generation_time=bad_generation_time,
                )
        for bad_time_units in ["gen"]:
            with pytest.raises(ValueError):
                species_trees.parse_species_tree(
                    good_tree,
                    time_units=bad_time_units,
                    initial_size=good_ne,
                    generation_time=good_generation_time,
                )


class TestSpeciesTreeExamples:
    """
    Tests that we get the expected value in simple examples.
    """

    def test_4_species_parse(self):
        good_tree = "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)"
        good_time_units = "myr"
        good_ne = 10000
        good_generation_time = 20
        demography = species_trees.parse_species_tree(
            good_tree,
            time_units=good_time_units,
            initial_size=good_ne,
            generation_time=good_generation_time,
        )
        assert isinstance(demography.populations, list)
        assert len(demography.populations) == 7
        for pop in demography.populations:
            assert isinstance(pop, msprime.demography.Population)
        assert isinstance(demography.events, list)
        assert len(demography.events) == 3
        for mm in demography.events:
            assert isinstance(mm, msprime.demography.PopulationSplit)

    def test_4_species_run(self):
        species_tree = (
            "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)"
        )
        spec = species_trees.parse_species_tree(
            species_tree,
            time_units="myr",
            initial_size=10000,
            generation_time=20,
        )

        # Take one sample from each population
        ts = msprime.sim_ancestry(
            samples={j: 1 for j in range(4)}, demography=spec, ploidy=1
        )

        assert ts.num_trees == 1
        assert ts.num_samples == 4
        assert ts.num_populations == 7
        for j, u in enumerate(ts.samples()):
            assert ts.node(u).population == j

        pops = list(ts.populations())
        assert pops[0].metadata["name"] == "human"
        assert pops[1].metadata["name"] == "chimpanzee"
        assert pops[2].metadata["name"] == "gorilla"
        assert pops[3].metadata["name"] == "orangutan"
        assert pops[4].metadata["name"] == "pop_4"
        assert pops[5].metadata["name"] == "pop_5"
        assert pops[6].metadata["name"] == "pop_6"

        # Use the population names to get the samples
        samples = dict(human=4, gorilla=2)
        ts = msprime.sim_ancestry(samples=samples, demography=spec)
        assert ts.num_trees == 1
        assert ts.num_samples == 12
        for j, u in enumerate(ts.samples()):
            pop = 0 if j < 8 else 2
            assert ts.node(u).population == pop

        # Order of keywords is respected
        ts = msprime.sim_ancestry(samples={"gorilla": 2, "human": 4}, demography=spec)
        assert ts.num_trees == 1
        assert ts.num_samples == 12
        for j, u in enumerate(ts.samples()):
            pop = 2 if j < 4 else 0
            assert ts.node(u).population == pop


class TestStarbeastParsingErrors:
    """
    Tests for parsing of species trees in nexus format, written by
    StarBEAST.
    """

    def test_bad_tree(self):
        bad_trees = []
        tree_file = "tests/data/species_trees/101g_nucl_conc_unconst.combined.nwk.tre"
        with open(tree_file) as f:
            bad_trees.append(f.read())
        good_nexus = "#NEXUS\n\n"
        good_nexus += "Begin taxa;\n"
        good_nexus += "    Dimensions ntax=3;\n"
        good_nexus += "    Taxlabels\n"
        good_nexus += "           spc01\n"
        good_nexus += "           spc02\n"
        good_nexus += "           spc03\n"
        good_nexus += "           ;\n"
        good_nexus += "End;\n"
        good_nexus += "Begin trees;\n"
        good_nexus += "    Translate\n"
        good_nexus += "     1 spc01,\n"
        good_nexus += "     2 spc02,\n"
        good_nexus += "     3 spc03\n"
        good_nexus += "     ;\n"
        good_nwk = "tree TREE1 = ((1[&dmv={0.1}]:1,2[&dmv={0.2}]:1)[&dmv={0.3}]"
        good_nexus += "End;\n"
        bad_trees.append(good_nexus.replace("#NEXUS", "#NEXU"))
        bad_trees.append(good_nexus.replace("#NEXUS", "NEXUS"))
        bad_trees.append(good_nexus.replace("tree TREE1", "tre TREE1"))
        bad_trees.append(good_nexus.replace("End;", ""))
        bad_trees.append(good_nexus.replace("Translate", "T"))
        bad_trees.append(good_nexus.replace("2 spc02,", "2 spc02"))
        bad_trees.append(good_nexus.replace("2 spc02,", "2 spc02 asdf,"))
        bad_trees.append(good_nexus.replace("2 spc02,", "2 spc03,"))
        bad_trees.append(good_nexus.replace("2 spc02,", "spc02 2,"))
        bad_trees.append(good_nexus.replace("2 spc02,", "asdf2 spc02,"))
        bad_trees.append(good_nexus.replace("2 spc02,", "spc02; 2,"))
        bad_trees.append(good_nexus.replace(";\n", ""))
        bad_trees.append(good_nexus.replace("Taxlabels", "Begin trees;"))
        bad_trees.append(good_nexus.replace("dmv", "emv"))
        bad_trees.append(good_nexus.replace("[", ""))
        bad_trees.append(good_nexus.replace("[", "").replace("]", ""))
        bad_trees.append(good_nexus.replace("=", ""))
        bad_trees.append(good_nexus.replace("Begin taxa", "Begin trees"))
        bad_trees.append(good_nexus.replace("Begin trees", "Begin taxa"))
        bad_trees.append(good_nexus.replace("[&dmv={0.5}]", ""))
        bad_trees.append(good_nexus.replace("[&dmv={0.1}]", ""))
        bad_trees.append(good_nexus.replace("[&dmv={0.1}]", "[&dmv={asdf}]"))
        bad_trees.append(good_nexus.replace(":1,2[&dmv", ":1, 2[&dmv"))
        bad_trees.append(good_nexus.replace(good_nwk, good_nwk + good_nwk))
        good_generation_time = 5
        for bad_tree in bad_trees:
            with pytest.raises(ValueError):
                species_trees.parse_starbeast(
                    tree=bad_tree, generation_time=good_generation_time
                )

    def test_bad_annotations(self):
        good = "((1[&dmv={0.1}]:1,2[&dmv={0.2}]:1)[&dmv={0.3}])"
        assert species_trees.strip_extra_annotations(good) == good
        bad_examples = [
            # No annotations
            "((1:1,2:1)",
            # Mismatched annotations
            "((1[]:1,2[]:1)[]]",
            "((1[]:1,2[]:1)[",
            "((1[]:1,2[]:1)]",
            # Missing all dmvs
            "((1[]:1,2[]:1)[]",
            # Missing closing }
            "((1[&dmv={]:1,2[]:1)[]",
        ]
        for example in bad_examples:
            with pytest.raises(ValueError):
                species_trees.strip_extra_annotations(example)

    def test_bad_annotations_in_tree(self):
        name_map = {f"n{j}": f"n{j}" for j in range(3)}
        good = "(n1[&dmv={1}]:1.14,n2[&dmv={1}]:1.14)[&dmv={1}]"
        spec = species_trees.process_starbeast_tree(good, 1, name_map)
        assert len(spec.populations) == 3
        assert len(spec.events) == 1
        bad_examples = [
            # Missing one dmv
            "(n1[&dmv={1}]:1.14,n2[&dmv={1}]:1.14)[&={1}]",
            # No annotation
            "(n1[&dmv={1}]:1.14,n2[&dmv={1}]:1.14)",
        ]
        for example in bad_examples:
            with pytest.raises(ValueError):
                species_trees.process_starbeast_tree(example, 1, name_map)

    def test_bad_translation(self):
        good = "translate 1 spc1, 2 spc2, 3 spc3"
        assert species_trees.parse_translate_command(good) == {
            "1": "spc1",
            "2": "spc2",
            "3": "spc3",
        }
        bad_examples = [
            "translate 1,",
            "translate 1 spc1 more, 2 spc2",
            "translate 1 spc1, 1 spc2",
            "translate 1 spc1, 2 spc1",
        ]
        for example in bad_examples:
            with pytest.raises(ValueError):
                species_trees.parse_translate_command(example)

    def test_bad_parameter(self):
        with open("tests/data/species_trees/91genes_species_rev.tre") as f:
            good_tree = f.read()
            good_time_units = "myr"
            for bad_time_units in [-3, "asdf", ["myr"], "gen"]:
                with pytest.raises(ValueError):
                    species_trees.parse_starbeast(
                        tree=f.read(),
                        time_units=bad_time_units,
                        generation_time=5,
                    )
            for bad_generation_time in [-3, "sdf"]:
                with pytest.raises(ValueError):
                    species_trees.parse_starbeast(
                        tree=good_tree,
                        time_units=good_time_units,
                        generation_time=bad_generation_time,
                    )
            for bad_generation_time in [None, {}]:
                with pytest.raises(TypeError):
                    species_trees.parse_starbeast(
                        tree=good_tree,
                        time_units=good_time_units,
                        generation_time=bad_generation_time,
                    )


class TestStarbeastExamples:
    """
    Tests for known examples in starbeast format.
    """

    def test_12_species(self):
        with open("tests/data/species_trees/91genes_species_rev.tre") as f:
            good_tree = f.read()
            good_time_units = "myr"
            good_generation_time = 5
            spec = species_trees.parse_starbeast(
                tree=good_tree,
                time_units=good_time_units,
                generation_time=good_generation_time,
            )
            assert len(spec.populations) == 23
            for pop in spec.populations[:12]:
                species_name = pop.name
                assert species_name.startswith("spc")
                assert species_name[3:].isnumeric()
            assert len(spec.events) == 11
            for mm in spec.events:
                assert isinstance(mm, msprime.PopulationSplit)
