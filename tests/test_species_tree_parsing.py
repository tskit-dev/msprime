#
# Copyright (C) 2018 University of Oxford
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


class TestSpeciesTreeParsingErrors(unittest.TestCase):
    """
    Tests for parsing of species trees in newick format.
    """
    def test_bad_tree(self):
        good_ne = 10000
        for bad_tree in [None, {}, 123, "(((h:5,c:5):3,g:8)\n:10,o:18)", "asdf"]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    species_tree=bad_tree,
                    Ne=good_ne
                    )

    def test_bad_parameter(self):
        good_tree = "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)"
        good_branch_length_units = "myr"
        good_ne = 10000
        good_generation_time = 5
        for bad_branch_length_units in [-3, "asdf", ["myr"]]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    species_tree=good_tree,
                    branch_length_units=bad_branch_length_units,
                    Ne=good_ne,
                    generation_time=good_generation_time
                    )
        for bad_ne in [None, -3]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    species_tree=good_tree,
                    branch_length_units=good_branch_length_units,
                    Ne=bad_ne,
                    generation_time=good_generation_time
                    )
        for bad_generation_time in [None, -3]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    species_tree=good_tree,
                    branch_length_units=good_branch_length_units,
                    Ne=good_ne,
                    generation_time=bad_generation_time
                    )
        for bad_branch_length_units in ["gen"]:
            with self.assertRaises(ValueError):
                msprime.parse_species_tree(
                    species_tree=good_tree,
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
                species_tree=good_tree,
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
                        species_tree=bad_tree,
                        generation_time=good_generation_time
                        )

    def test_bad_parameter(self):
        with open("tests/data/species_trees/91genes_species_rev.tre") as f:
            good_tree = f.read()
            good_branch_length_units = "myr"
            for bad_branch_length_units in [-3, "asdf", ["myr"], "gen"]:
                with self.assertRaises(ValueError):
                    msprime.parse_starbeast(
                            species_tree=f.read(),
                            branch_length_units=bad_branch_length_units,
                            generation_time=5)
            for bad_generation_time in [-3]:
                with self.assertRaises(ValueError):
                    msprime.parse_starbeast(
                        species_tree=good_tree,
                        branch_length_units=good_branch_length_units,
                        generation_time=bad_generation_time
                        )
            for bad_generation_time in [None, {}]:
                with self.assertRaises(TypeError):
                    msprime.parse_starbeast(
                        species_tree=good_tree,
                        branch_length_units=good_branch_length_units,
                        generation_time=bad_generation_time
                        )

    def test_output(self):
        with open("tests/data/species_trees/91genes_species_rev.tre") as f:
            good_tree = f.read()
            good_branch_length_units = "myr"
            good_generation_time = 5
            parsed_tuple = msprime.parse_starbeast(
                    species_tree=good_tree,
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
