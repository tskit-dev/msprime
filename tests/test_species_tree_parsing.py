#
# Copyright (C) 2020 University of Oxford
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
Tests for the parsing of species trees in newick and starbeast format.
"""
import unittest

import msprime


def get_non_binary_tree(n):
    demographic_events = [
        msprime.SimpleBottleneck(time=0.1, population=0, proportion=0.5)]
    ts = msprime.simulate(
        n, demographic_events=demographic_events, random_seed=1)
    tree = ts.first()
    found = False
    for u in tree.nodes():
        if len(tree.children(u)) > 2:
            found = True
    assert found
    return tree


class TestSpeciesTreeRoundTrip(unittest.TestCase):
    """
    Tests that we get what we expect when we parse trees produced from
    msprime/tskit.
    """
    def verify(self, tree, newick=None, Ne=1, branch_length_units="gen",
               generation_time=None):
        if newick is None:
            newick = tree.newick()
        population_configurations, demographic_events = msprime.parse_species_tree(
            newick, Ne=Ne, branch_length_units=branch_length_units,
            generation_time=generation_time)
        self.assertEqual(len(population_configurations), tree.num_samples())
        for pop_config in population_configurations:
            self.assertEqual(pop_config.initial_size, Ne)
            self.assertEqual(pop_config.growth_rate, 0)
            self.assertIn("species_name", pop_config.metadata)
            self.assertTrue(pop_config.metadata["species_name"].isdigit())

        # Population IDs are mapped to leaves as they are encountered in a postorder
        # traversal.
        pop_id_map = {}
        k = 0
        for u in tree.nodes(order="postorder"):
            if tree.is_leaf(u):
                pop_id_map[u] = k
                k += 1
            else:
                pop_id_map[u] = pop_id_map[tree.left_child(u)]

        # We should have demographic events for every non-unary internal node, and
        # events should be output in increasing time order.
        j = 0
        for node in [u for u in tree.nodes(order="timeasc")]:
            children = tree.children(node)
            if len(children) > 1:
                self.assertEqual(node, tree.mrca(children[0], children[1]))
                dest = pop_id_map[node]
                for child in children[1:]:
                    event = demographic_events[j]
                    j += 1
                    self.assertIsInstance(event, msprime.MassMigration)
                    self.assertAlmostEqual(event.time, tree.time(node))
                    source = pop_id_map[child]
                    self.assertEqual(event.source, source)
                    self.assertEqual(event.dest, dest)

        self.assertEqual(j, len(demographic_events))

    def test_n2_binary(self):
        tree = msprime.simulate(2, random_seed=2).first()
        self.verify(tree)

    def test_n2_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, 0), (0, 1)], random_seed=2)
        self.verify(ts.first(), Ne=5)

    def test_n5_binary(self):
        ts = msprime.simulate(5, random_seed=2)
        tree = ts.first()
        self.verify(tree, Ne=1)

    def test_n5_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, j) for j in range(5)], random_seed=2)
        self.verify(ts.first(), Ne=10)

    def test_n7_binary(self):
        ts = msprime.simulate(7, random_seed=2)
        tree = ts.first()
        self.verify(tree, Ne=11)

    def test_n7_binary_embedded_whitespace(self):
        # Check for embedded whitespace in the newick string
        tree = msprime.simulate(7, random_seed=2).first()
        newick = tree.newick()
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
        self.verify(tree, Ne=11)

    def test_n10_non_binary(self):
        tree = get_non_binary_tree(10)
        self.verify(tree, Ne=3.1234)

    def test_n10_binary_years(self):
        ts = msprime.simulate(10, random_seed=2)
        generation_time = 5
        tree = ts.first()
        self.verify(tree, Ne=1, branch_length_units="yr", generation_time=1)
        tables = ts.dump_tables()
        times = tables.nodes.time
        flags = tables.nodes.flags
        scaled_times = [generation_time * time for time in times]
        tables.nodes.set_columns(flags=flags, time=scaled_times)
        ts = tables.tree_sequence()
        scaled_tree = ts.first()
        self.verify(tree, newick=scaled_tree.newick(), Ne=1,
                    branch_length_units="yr", generation_time=generation_time)

    def test_n10_binary_million_years(self):
        ts = msprime.simulate(10, random_seed=2)
        generation_time = 5
        tree = ts.first()
        tables = ts.dump_tables()
        times = tables.nodes.time
        flags = tables.nodes.flags
        scaled_times = [generation_time * time * 1E-6 for time in times]
        tables.nodes.set_columns(flags=flags, time=scaled_times)
        ts = tables.tree_sequence()
        scaled_tree = ts.first()
        self.verify(tree, newick=scaled_tree.newick(), Ne=1,
                    branch_length_units="myr", generation_time=generation_time)


def make_nexus(tree, pop_size_map):
    """
    Returns the specified tree formatted as starbeast compatable nexus.
    """
    node_labels = {}
    for u in tree.nodes():
        name = ""
        if tree.is_leaf(u):
            name = str(u)
        node_labels[u] = f"{name}[&dmv={{{pop_size_map[u]}}}]"
    newick = tree.newick(node_labels=node_labels)
    out = "#NEXUS\n"
    out += "Begin trees;\n"
    out += "tree TREE1 = " + newick
    out += "End;"
    return out


class TestStarbeastRoundTrip(unittest.TestCase):
    """
    Tests that we get what we expect when we parse trees produced from
    msprime/tskit.
    """
    def verify(self, tree, pop_size_map, nexus=None):
        nexus = make_nexus(tree, pop_size_map)

        population_configurations, demographic_events = msprime.parse_starbeast(
            nexus, "yr", 1)
        self.assertEqual(len(population_configurations), tree.num_samples())

        # Population IDs are mapped to leaves as they are encountered in a postorder
        # traversal.
        pop_id_map = {}
        k = 0
        for u in tree.nodes(order="postorder"):
            if tree.is_leaf(u):
                pop_id_map[u] = k
                k += 1
            else:
                pop_id_map[u] = pop_id_map[tree.left_child(u)]

        for u in tree.leaves():
            pop_config = population_configurations[pop_id_map[u]]
            self.assertEqual(pop_config.initial_size, pop_size_map[u])
            self.assertEqual(pop_config.growth_rate, 0)
            # TODO check the metadata: we should be encoding the species name

        # We should have demographic events for every non-unary internal node, and
        # events should be output in increasing time order.
        j = 0
        for node in [u for u in tree.nodes(order="timeasc")]:
            children = tree.children(node)
            if len(children) > 1:
                dest = pop_id_map[node]
                for child in children[1:]:
                    event = demographic_events[j]
                    j += 1
                    self.assertIsInstance(event, msprime.MassMigration)
                    self.assertAlmostEqual(event.time, tree.time(node))
                    source = pop_id_map[child]
                    self.assertEqual(event.source, source)
                    self.assertEqual(event.dest, dest)
                event = demographic_events[j]
                j += 1
                self.assertIsInstance(event, msprime.PopulationParametersChange)
                self.assertAlmostEqual(event.time, tree.time(node))
                self.assertEqual(event.population, dest)
                self.assertEqual(event.growth_rate, None)
                self.assertEqual(event.initial_size, pop_size_map[node])

        self.assertEqual(j, len(demographic_events))

    def test_n2_binary(self):
        tree = msprime.simulate(2, random_seed=2).first()
        self.verify(tree, {u: 1 for u in tree.nodes()})

    def test_n2_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, 0), (0, 1)], random_seed=2)
        tree = ts.first()
        self.verify(tree, {u: 2.123 for u in tree.nodes()})

    def test_n5_binary(self):
        ts = msprime.simulate(5, random_seed=2)
        tree = ts.first()
        self.verify(tree, {u: 1 + u for u in tree.nodes()})

    def test_n5_binary_non_ultrametric(self):
        ts = msprime.simulate(samples=[(0, j) for j in range(5)], random_seed=2)
        tree = ts.first()
        self.verify(tree, {u: 1 / (1 + u) for u in tree.nodes()})

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


class TestSpeciesTreeParsingErrors(unittest.TestCase):
    """
    Tests for parsing of species trees in newick format.
    """
    def test_bad_tree(self):
        for bad_type in [None, {}, 123]:
            with self.assertRaises(TypeError):
                msprime.parse_species_tree(tree=bad_type, Ne=1)

        bad_trees = [
            "", ";",  "abcd", ";;;", "___", "∞",
            "(", ")", "()", "( )", "(()())",
            "((3:0.39,5:0.39]:1.39,(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
            "((3:0.39,5:0.39(:1.39,(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
            "((3:0.39,5:0.39,:1.39,(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
            "(4:0.47,(1:0.18,2:0.18):0.29):1.31);",
        ]
        for bad_tree in bad_trees:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(tree=bad_tree, Ne=1)

    def test_bad_parameter(self):
        good_tree = "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)"
        good_branch_length_units = "myr"
        good_ne = 10000
        good_generation_time = 5
        for bad_branch_length_units in [-3, "asdf", ["myr"]]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    good_tree,
                    branch_length_units=bad_branch_length_units,
                    Ne=good_ne,
                    generation_time=good_generation_time
                    )
        for bad_ne in [None, -3]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    good_tree,
                    branch_length_units=good_branch_length_units,
                    Ne=bad_ne,
                    generation_time=good_generation_time
                    )
        for bad_generation_time in [None, -3]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    good_tree,
                    branch_length_units=good_branch_length_units,
                    Ne=good_ne,
                    generation_time=bad_generation_time
                    )
        for bad_branch_length_units in ["gen"]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    good_tree,
                    branch_length_units=bad_branch_length_units,
                    Ne=good_ne,
                    generation_time=good_generation_time
                    )

    def test_output(self):
        good_tree = "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)"
        good_branch_length_units = "myr"
        good_ne = 10000
        good_generation_time = 20
        parsed_tuple = msprime.parse_species_tree(
                good_tree,
                branch_length_units=good_branch_length_units,
                Ne=good_ne,
                generation_time=good_generation_time
                )
        self.assertEqual(len(parsed_tuple), 2)
        self.assertIs(type(parsed_tuple[0]), list)
        self.assertEqual(len(parsed_tuple[0]), 4)
        for pc in parsed_tuple[0]:
            self.assertIs(type(pc), msprime.simulations.PopulationConfiguration)
        self.assertIs(type(parsed_tuple[1]), list)
        self.assertEqual(len(parsed_tuple[1]), 3)
        for mm in parsed_tuple[1]:
            self.assertIs(type(mm), msprime.simulations.MassMigration)


class TestStarBeastParsingErrors(unittest.TestCase):
    """
    Tests for parsing of species trees in nexus format, written by
    StarBEAST.
    """
    def test_bad_tree(self):
        tree_file = "tests/data/species_trees/101g_nucl_conc_unconst.combined.nwk.tre"
        with open(tree_file) as f:
            nwk = f.read()
        good_generation_time = 5
        for bad_tree in [None, {}, 123, nwk]:
            with self.assertRaises(ValueError):
                msprime.parse_starbeast(
                        tree=bad_tree,
                        generation_time=good_generation_time
                        )

    def test_bad_parameter(self):
        with open("tests/data/species_trees/91genes_species_rev.tre") as f:
            good_tree = f.read()
            good_branch_length_units = "myr"
            for bad_branch_length_units in [-3, "asdf", ["myr"], "gen"]:
                with self.assertRaises(ValueError):
                    msprime.parse_starbeast(
                            tree=f.read(),
                            branch_length_units=bad_branch_length_units,
                            generation_time=5)
            for bad_generation_time in [-3]:
                with self.assertRaises(ValueError):
                    msprime.parse_starbeast(
                        tree=good_tree,
                        branch_length_units=good_branch_length_units,
                        generation_time=bad_generation_time
                        )
            for bad_generation_time in [None, {}]:
                with self.assertRaises(TypeError):
                    msprime.parse_starbeast(
                        tree=good_tree,
                        branch_length_units=good_branch_length_units,
                        generation_time=bad_generation_time
                        )

    def test_output(self):
        with open("tests/data/species_trees/91genes_species_rev.tre") as f:
            good_tree = f.read()
            good_branch_length_units = "myr"
            good_generation_time = 5
            parsed_tuple = msprime.parse_starbeast(
                    tree=good_tree,
                    branch_length_units=good_branch_length_units,
                    generation_time=good_generation_time
                    )
            self.assertEqual(len(parsed_tuple), 2)
            self.assertIs(type(parsed_tuple[0]), list)
            self.assertEqual(len(parsed_tuple[0]), 12)
            for pc in parsed_tuple[0]:
                self.assertIs(type(pc), msprime.simulations.PopulationConfiguration)
            self.assertIs(type(parsed_tuple[1]), list)
            self.assertEqual(len(parsed_tuple[1]), 22)
            event_types = [msprime.simulations.MassMigration]
            event_types.append(msprime.simulations.PopulationParametersChange)
            for mm in parsed_tuple[1]:
                self.assertIn(type(mm), event_types)
