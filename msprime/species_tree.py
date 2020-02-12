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
Module responsible for parsing species trees.
"""
import sys
import re

import msprime

class Tree(object):
    """
    The tree object, consisting of multiple edges. Trees are checked for whether they are
    ultrametric; only ultrametric trees are allowed. Species IDs are stored as node IDs of
    extant edges.
    """
    def __init__(self, newick_string):
        self.newick_string = newick_string
        self.edges = []

    def get_newick_string(self):
        return self.newick_string

    def parse_newick_string(self):
        working_newick_string = self.newick_string
        number_of_internal_nodes = 0
        number_of_edges = 0

        # Remove comments from the tree string.
        pattern = re.compile("\[.*?\]")
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(working_newick_string)
            if hit != None:
                working_newick_string = working_newick_string.replace(hit.group(0),"")

        # Check whether a branch above the root is present, and if so, remove it.
        if working_newick_string[0:2] == "((" and working_newick_string[-1] == ")" and working_newick_string[-2] != ")":
            level = 0
            newick_string_tail_start_pos = 0
            newick_string_tail = ""
            for pos in range(len(working_newick_string)):
                if working_newick_string[pos] == "(":
                    level += 1
                if working_newick_string[pos] == ")":
                    level -= 1
                if level == 1 and pos > 1:
                    newick_string_tail_start_pos = pos
                    newick_string_tail = working_newick_string[pos+1:]
                    break
            if newick_string_tail.count(",") == 0:
                working_newick_string = working_newick_string[1:newick_string_tail_start_pos+1]

        # Parse the bifurcating part of the tree.
        assert ":" in working_newick_string, 'It appears that the tree string does not include branch lengths!'
        pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?):([\d\.Ee-]+?),([a-zA-Z0-9_\.\-]+?):([\d\.Ee-]+?)\)")
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(working_newick_string)
            if hit != None:
                number_of_internal_nodes += 1
                number_of_edges += 2
                internal_node_id = "internalNode" + str(number_of_internal_nodes) + "X"
                edge1 = Edge("edge" + str(number_of_edges-1) + "X")
                edge1.set_node_ids([internal_node_id, hit.group(1)])
                edge1.set_length(float(hit.group(2)))
                edge2 = Edge("edge" + str(number_of_edges) + "X")
                edge2.set_node_ids([internal_node_id, hit.group(3)])
                edge2.set_length(float(hit.group(4)))
                self.edges.append(edge1)
                self.edges.append(edge2)
                working_newick_string = working_newick_string.replace(hit.group(0), internal_node_id)

        # Make sure the remaining string includes a single node and use this node id to determine root edges.
        pattern_rooted = re.compile("^internalNode\d+X$")
        hit_rooted = pattern_rooted.search(working_newick_string)
        err = 'The species tree string could not be parsed! '
        err += 'This can occur when the species tree is not strictly bifurcating and contains polytomies. '
        err += 'The remaining unparsed part of the newick string: \"{}\"'.format(working_newick_string)
        assert hit_rooted, err
        root_node_id = hit_rooted.group(0)
        for edge in self.get_edges():
            if edge.get_node_ids()[0] == root_node_id:
                edge.set_is_root_edge(True)
            else:
                edge.set_is_root_edge(False)

    def parse_extended_newick_string(self):
        working_newick_string = self.newick_string
        number_of_internal_nodes = 0
        number_of_edges = 0

        # Check whether a branch above the root is present, and if so, remove it.
        if working_newick_string[0:2] == "((" and working_newick_string[-1] == ")" and working_newick_string[-2] != ")":
            level = 0
            newick_string_tail_start_pos = 0
            newick_string_tail = ""
            for pos in range(len(working_newick_string)):
                if working_newick_string[pos] == "(":
                    level += 1
                if working_newick_string[pos] == ")":
                    level -= 1
                if level == 1 and pos > 1:
                    newick_string_tail_start_pos = pos
                    newick_string_tail = working_newick_string[pos+1:]
                    break
            if newick_string_tail.count(",") == 0:
                working_newick_string = working_newick_string[1:newick_string_tail_start_pos+1]

        # Parse the bifurcating part of the tree.
        assert ":" in working_newick_string, 'It appears that the tree string does not include branch lengths!'
        assert "[" in working_newick_string, 'It appears that the tree string does not include annotation!'
        pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?)\[\&dmv=\{([\d\.]+?)\}\]:([\d\.Ee-]+?),([a-zA-Z0-9_\.\-]+?)\[\&dmv=\{([\d\.]+?)\}\]:([\d\.Ee-]+?)\)")
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(working_newick_string)
            if hit != None:
                number_of_internal_nodes += 1
                number_of_edges += 2
                internal_node_id = "internalNode" + str(number_of_internal_nodes) + "X"
                edge1 = Edge("edge" + str(number_of_edges-1) + "X")
                edge1.set_node_ids([internal_node_id, hit.group(1)])
                edge1.set_dmv(float(hit.group(2)))
                edge1.set_length(float(hit.group(3)))
                edge2 = Edge("edge" + str(number_of_edges) + "X")
                edge2.set_node_ids([internal_node_id, hit.group(4)])
                edge2.set_dmv(float(hit.group(5)))
                edge2.set_length(float(hit.group(6)))
                self.edges.append(edge1)
                self.edges.append(edge2)
                working_newick_string = working_newick_string.replace(hit.group(0), internal_node_id)

        # Make sure the remaining string includes a single node and use this node id to determine root edges.
        pattern_rooted = re.compile("^internalNode\d+X$")
        hit_rooted = pattern_rooted.search(working_newick_string)
        assert hit_rooted, 'The newick tree string could not be parsed! The remaining unparsed part of the newick string: \"{}\"'.format(working_newick_string)
        root_node_id = hit_rooted.group(0)
        for edge in self.get_edges():
            if edge.get_node_ids()[0] == root_node_id:
                edge.set_is_root_edge(True)
            else:
                edge.set_is_root_edge(False)

    def get_edges(self):
        return self.edges

    def get_number_of_edges(self):
        return len(self.edges)

    def get_number_of_extant_edges(self):
        number_of_extant_edges = 0
        for edge in self.edges:
            if edge.get_is_extant():
                number_of_extant_edges += 1
        return number_of_extant_edges

    def set_extant_progeny_ids(self):
        for edge in self.get_edges():
            if edge.get_is_extant():
                species_id = edge.get_node_ids()[1]
                this_edge = edge
                species_id_added_to_root_edge = False
                while species_id_added_to_root_edge == False:
                    this_edge.add_extant_progeny_id(species_id)
                    if this_edge.get_is_root_edge():
                        species_id_added_to_root_edge = True
                    else:
                        for other_edge in self.get_edges():
                            if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                                parent_edge = other_edge
                                break
                    this_edge = parent_edge

    def set_times(self):
        # Get the durations between root and extant edges.
        total_edge_lengths = []
        for edge in self.get_edges():
            if edge.get_is_extant():
                total_edge_length = edge.get_length()
                if edge.get_is_root_edge():
                    total_edge_lengths.append(total_edge_length)
                else:
                    root_edge_found = False
                    this_edge = edge
                    while root_edge_found == False:
                        for other_edge in self.get_edges():
                            if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                                parent_edge = other_edge
                                break
                        assert parent_edge, 'The parent edge for edge {} could not be found'.format(this_edge.get_id())
                        total_edge_length += parent_edge.get_length()
                        if parent_edge.get_is_root_edge():
                            root_edge_found = True
                            total_edge_lengths.append(total_edge_length)
                        else:
                            this_edge = parent_edge
        # Make sure that the tree is at least roughly ultrametric. A maximum difference of 1% between the shortest and longest root-to-tip distances is allowed.
        max_total_edge_length = max(total_edge_lengths)
        if max_total_edge_length - min(total_edge_lengths) > 0.01 * max_total_edge_length:
            raise NotImplementedError(
                    'The tree appears not to be ultrametric. '
                    'The use of non-ultrametric trees is not yet supported.')
        # Extend terminal edges if necessary.
        for edge in self.get_edges():
            if edge.get_is_extant():
                total_edge_length = edge.get_length()
                if edge.get_is_root_edge():
                    edge.set_length(edge.get_length() + max_total_edge_length - total_edge_length)
                else:
                    root_edge_found = False
                    this_edge = edge
                    while root_edge_found == False:
                        for other_edge in self.get_edges():
                            if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                                parent_edge = other_edge
                                break
                        assert parent_edge, 'The parent edge for edge {} could not be found'.format(this_edge.get_id())
                        total_edge_length += parent_edge.get_length()
                        if parent_edge.get_is_root_edge():
                            root_edge_found = True
                            edge.set_length(round(edge.get_length() + max_total_edge_length - total_edge_length,8))
                        else:
                            this_edge = parent_edge
        # First specify the edges for which the parents still need to be identified.
        for edge in self.get_edges():
            if edge.get_is_root_edge():
                edge.set_parent_needs_times(False)
            else:
                edge.set_parent_needs_times(True)
        # Set the times of all edges.
        for edge in self.get_edges():
            if edge.get_is_extant() == True:
                edge.set_termination(0.0)
                edge.set_origin(edge.get_length())
                this_edge = edge
                while this_edge.get_parent_needs_times():
                    for other_edge in self.get_edges():
                        if other_edge.get_node_ids()[1] == this_edge.get_node_ids()[0]:
                            parent_edge = other_edge
                            break
                    assert parent_edge, 'The parent edge for edge {} could not be found'.format(this_edge.get_id())
                    parent_edge.set_termination(this_edge.get_origin())
                    parent_edge.set_origin(this_edge.get_origin() + parent_edge.get_length())
                    this_edge.set_parent_needs_times = False
                    this_edge = parent_edge

    def get_origin(self):
        for edge in self.get_edges():
            if edge.get_is_root_edge():
                return edge.get_origin()

    def info(self):
        info_string = ''
        info_string += 'Tree'.ljust(20)
        info_string += '\n'
        info_string += 'Number of edges:'.ljust(20)
        info_string += str(self.get_number_of_edges())
        info_string += '\n'
        return info_string


class Edge(object):
    """
    The edge object, with two IDs corresponding to the nodes at both ends of the edge,
    and ages of origin and termination.
    """
    def __init__(self, id):
        self.id = id
        self.node_ids = []
        self.extant_progeny_ids = []
        self.origin = None
        self.termination = None
        self.dmv = None
        self.parent_needs_times = True

    def get_id(self):
        return self.id

    def set_node_ids(self, node_ids):
        self.node_ids = node_ids

    def get_node_ids(self):
        return self.node_ids

    def get_is_extant(self):
        if self.node_ids[1][0:12] == 'internalNode':
            return False
        else:
            return True

    def set_length(self, length):
        self.length = length

    def get_length(self):
        return self.length

    def set_is_root_edge(self, is_root_edge):
        self.is_root_edge = is_root_edge

    def get_is_root_edge(self):
        return self.is_root_edge

    def add_extant_progeny_id(self, extant_progeny_id):
        self.extant_progeny_ids.append(extant_progeny_id)

    def get_extant_progeny_ids(self):
        return self.extant_progeny_ids

    def set_termination(self, termination):
        self.termination = termination

    def get_termination(self):
        return self.termination

    def set_origin(self, origin):
        self.origin = origin

    def get_origin(self):
        return self.origin

    def set_parent_needs_times(self, parent_needs_times):
        self.parent_needs_times = parent_needs_times

    def get_parent_needs_times(self):
        return self.parent_needs_times

    def set_dmv(self, dmv):
        self.dmv = dmv

    def get_pop_size(self, generations_per_branch_length_unit):
        return self.dmv * generations_per_branch_length_unit

    def info(self):
        info_string = ''
        info_string += 'Edge id:'.ljust(28)
        info_string += self.id
        info_string += '\n'
        info_string += 'Edge node 1 id:'.ljust(28)
        info_string += self.node_ids[0]
        info_string += '\n'
        info_string += 'Edge node 2 id:'.ljust(28)
        info_string += self.node_ids[1]
        info_string += '\n'
        info_string += 'Edge length:'.ljust(28)
        info_string += str(self.length)
        info_string += '\n'
        if self.dmv != None:
            info_string += 'Edge dmv:'.ljust(28)
            info_string += str(self.dmv)
            info_string += '\n'            
        info_string += 'Edge origin:'.ljust(28)
        info_string += str(self.origin)
        info_string += '\n'
        info_string += 'Edge termination:'.ljust(28)
        info_string += str(self.termination)
        info_string += '\n'
        info_string += 'Edge is extant:'.ljust(28)
        info_string += str(self.get_is_extant())
        info_string += '\n'
        info_string += 'Edge is root edge:'.ljust(28)
        info_string += str(self.is_root_edge)
        info_string += '\n'
        info_string += 'Edge extant progeny ids:'.ljust(28)
        for extant_progeny_id in self.extant_progeny_ids:
            info_string += '{}, '.format(extant_progeny_id)
        info_string = info_string[:-2]
        info_string += '\n'
        return info_string


def parse_species_tree(
        filename=None,
        branch_length_units="gen",
        sample_size=None,
        Ne=None,
        generation_time=None):
    """
    Method to parse species trees from files in Newick https://en.wikipedia.org/wiki/Newick_format format or 
    Nexus (https://en.wikipedia.org/wiki/Nexus_file) format with embedded Newick strings.

    Trees must be rooted, binary, and ultrametric and branch lengths must be included and correspond to time,
    either in units of millions of years ("myr"), years ("yr"), or generations ("gen"; default). Leafs must be
    named. An example for an accepted tree string in Newick format is this:
    (((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)
    The tree string can end with a semi-colon, but this is not required.

    We allow the following input:
    1) filename must be specified.
    2) If and only if the tree is not in StarBEAST format, Ne should be specified.
    3) If and only if the branch lengths are not in units of generations, the generation time should be specified.

    As msprime does not name populations and instead only assigns numbers to populations (starting with 0), the
    species names stored in the species tree are not used when defining instances of PopulationConfiguration.
    However, these instances will be sorted according to the alphabetical order of the corresponding species names,
    allowing the user to link species from the species tree with populations in msprime.

    The sample size is not used except that it is stored in the resulting instances of PopulationConfiguration.

    Introgression events, as stored in the extended Newick format by PhyloNet and possibly other programs, are not
    parsed and therefore not used.
    Rate matrices of continuous migration between co-existing lineages can not be stored in (extended) Newick or
    Nexus format and are therefore also not used here. However, these can be added to the model after population
    configurations and demographic events are defined with this function.
    """
    
    # Make sure a filename is specified.
    assert filename, 'A filename must be specified when calling function parse_species_tree().'

    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if not branch_length_units in allowed_branch_lenth_units:
        err = 'The specified units for branch lengths ("{}") are not accepted. '.format(branch_length_units)
        err += 'Accepted units are "myr" (millions of years), "yr" (years), and "gen" (generations).'
        raise ValueError(err)

    # Make sure that the sample size is either None or non-negative.
    if sample_size is not None and sample_size < 0:
        raise ValueError("Sample size must be >= 0")

    # Make sure that the population size is either None or positive.
    if Ne is not None and Ne <= 0:
        raise ValueError("Population size Ne must be > 0")

    # Make sure that the generation time is either None or positive.
    if generation_time is not None and generation_time <= 0:
        raise ValueError("Generation time must be > 0")

    # Make sure that the generation time is specified if and only if branch lengths are not in units of generations.
    if branch_length_units == "gen":
        assert generation_time == None, 'With branch lengths in units of generations ("gen"), generation_time should not be specified.'
    else:
        assert generation_time, 'With branch lengths in units of "{}", generation_time must be specified.'.format(branch_length_units)

    # Get the number of generations per branch length unit.
    if branch_length_units == "gen":
        generations_per_branch_length_unit = 1
    elif branch_length_units == "myr":
        generations_per_branch_length_unit = 1000000/generation_time
    else:
        generations_per_branch_length_unit = 1/generation_time

    # Read the input file.
    tree_string = None
    starbeast_format = False
    starbeast_root_pop_size = None
    infile = open(filename)
    inlines = infile.readlines()
    assert len(inlines) > 0, 'File {} is empty.'.format(filename)
    if inlines[0][0:6].lower() == '#nexus':
        # Assume the input is in nexus format. Maximally one tree string is read.
        in_tree = False
        translate_string = ""
        for line in inlines:
            if line.strip().lower() == 'begin trees;':
                in_tree = True
            elif line.strip().lower() == 'end;':
                in_tree = False
            elif in_tree and line.strip() is not '':
                clean_line = line.strip()
                if clean_line.lower().split()[0] == "translate":
                    in_translation = True
                    translate_string += clean_line[9:]
                elif in_translation:
                    if ";" in clean_line:
                        translate_string += clean_line.split(";")[0]
                        in_translation = False
                    else:
                        translate_string += clean_line
                if clean_line[0:4].lower() == "tree":
                    if '[' in clean_line and ']' in clean_line:
                        if "dmv=" in clean_line:
                            starbeast_format = True
                            # Remove all annotation except for dmv.
                            tree_string_raw = ''
                            in_comment = False
                            in_dmv = False
                            for x in range(len(clean_line)):
                                if clean_line[x] == '[':
                                    in_comment = True
                                    tree_string_raw += clean_line[x]
                                elif clean_line[x] == ']':
                                    in_comment = False
                                    tree_string_raw += clean_line[x]
                                elif in_comment:
                                    if clean_line[x-1] == "[" and clean_line[x] == "&":
                                        tree_string_raw += "&"
                                    if clean_line[x-5:x] == "dmv={":
                                        in_dmv = True
                                        tree_string_raw += "dmv={"
                                    if in_dmv:
                                        tree_string_raw += clean_line[x]
                                    if clean_line[x] == "}":
                                        in_dmv = False
                                else:
                                    tree_string_raw += clean_line[x]
                        else:
                            # Remove all annotation.
                            tree_string_raw = ''
                            in_comment = False
                            for letter in clean_line:
                                if letter == '[':
                                    in_comment = True
                                elif letter == ']':
                                    in_comment = False
                                elif in_comment == False:
                                    tree_string_raw += letter
                    if starbeast_format:
                        tree_patterns = re.search('(\(.+\))\[\&dmv=\{([\d\.]+?)\}\]',tree_string_raw)
                        tree_string = tree_patterns.group(1)
                        starbeast_root_pop_size = float(tree_patterns.group(2)) * generations_per_branch_length_unit
                    else:
                        tree_patterns = re.search('\(.+\)',tree_string_raw)
                        tree_string = tree_patterns.group(0)
                    # Use the translation block (if present) to back-translate species IDs in the tree string.
                    if translate_string != "":
                        translation_originals = []
                        translation_replaceds = []
                        translation_list = translate_string.split(",")
                        for item in translation_list:
                            item_list = item.split()
                            assert len(item_list) > 1, 'The translation block of file {} is malformed.'.format(filename)
                            assert len(item_list) == 2, 'Species IDs in the translation block of file {} appear to include whitespace. This is not supported.'.format(filename)
                            translation_originals.append(item_list[1])
                            translation_replaceds.append(item_list[0])
                        assert len(translation_originals) == len(set(translation_originals)), 'One or more species names in the translation block of file {} are duplicated.'.format(filename)
                        assert len(translation_replaceds) == len(set(translation_replaceds)), 'One or more translations in the translation block of file {} are duplicated.'.format(filename)
                        translation_originals_and_replaceds = translation_originals + translation_replaceds
                        err = 'One or more of the species names are identical to one or more of the translations in the translation block of file {}.'.format(filename)
                        assert len(translation_originals_and_replaceds) == len(set(translation_originals_and_replaceds)), err
                        for x in range(len(translation_originals)):
                            find_string = ',{}['.format(translation_replaceds[x])
                            replace_string = ',{}['.format(translation_originals[x])
                            tree_string = tree_string.replace(find_string,replace_string)
                            find_string = '({}['.format(translation_replaceds[x])
                            replace_string = '({}['.format(translation_originals[x])
                            tree_string = tree_string.replace(find_string,replace_string)
                    break

    else:
        # Test if the input is in newick format.
        assert len(inlines) == 1, 'File {} does not appear to be in Nexus format but contains more than one line.'.format(filename)
        assert inlines[0][0] == '(', 'Unexpected tree format in file {}. The first character should be "(".'.format(filename)
        assert inlines[0].count('(') == inlines[0].count(')'), 'The number of opening and closing parentheses differ in the tree string!'
        tree_string_raw = inlines[0]
        tree_patterns = re.search('\(.+\)',tree_string_raw)
        tree_string = tree_patterns.group(0)

    # Make sure a tree string is found.
    assert tree_string, 'No tree could be found in file {}.'.format(filename)

    # Make sure that between the beginning and the end of the tree string, there are always more parentheses opened than closed.
    parentheses_open = 0
    for char in tree_string[:-2]:
        if char == "(":
            parentheses_open += 1
        elif char == ")":
            parentheses_open -= 1
        assert parentheses_open > 0, 'The tree string in file {} appears to be malformed (possibly the first and last parentheses are missing).'.format(filename)

    # Now that we know if the tree is in starbeast format or not, make sure that Ne is specified if and only if the tree is not in starbeast format.
    if starbeast_format:
        assert Ne == None, 'With species trees in StarBEAST format, Ne should not be specified.'
    else:
        assert Ne, 'Ne should be specified.'
    
    # Use the tree string to generate a tree object.
    species_tree = Tree(tree_string)

    # Parse the newick tree string.
    if starbeast_format:
        species_tree.parse_extended_newick_string()
    else:
        species_tree.parse_newick_string()

    # Set the extant progeny ids for each branch.
    species_tree.set_extant_progeny_ids()

    # Set the origin and termination times of all branches.
    species_tree.set_times()

    # Get the age of the species tree.
    tree_origin = species_tree.get_origin()

    # Get a list of species IDs along with population sizes.
    species_ids = []
    terminal_pop_sizes = []
    for edge in species_tree.get_edges():
        if edge.get_is_extant():
            species_ids.append(edge.get_node_ids()[1])
            if starbeast_format:
                terminal_pop_sizes.append(edge.get_pop_size(generations_per_branch_length_unit))
            else:
                terminal_pop_sizes.append(Ne)
    lists_are_sorted = False
    while lists_are_sorted == False:
        lists_are_sorted = True
        for x in range(len(species_ids)-1):
            if species_ids[x] > species_ids[x+1]:
                species_ids[x], species_ids[x+1] = species_ids[x+1], species_ids[x]
                terminal_pop_sizes[x], terminal_pop_sizes[x+1] = terminal_pop_sizes[x+1], terminal_pop_sizes[x]
                lists_are_sorted = False
                break

    # Determine at which time which populations should merge.
    indices1 = []
    indices2 = []
    internal_pop_sizes = []
    divergence_times = []
    for edge in species_tree.get_edges():
        if edge.get_is_extant() == False:
            d1_species = []
            d2_species = []
            for other_edge in species_tree.get_edges():
                if edge.get_node_ids()[1] == other_edge.get_node_ids()[0]:
                    if d1_species == []:
                        d1_species = other_edge.get_extant_progeny_ids()
                    elif d2_species == []:
                        d2_species = other_edge.get_extant_progeny_ids()
                        break
            d1_species.sort()
            d2_species.sort()
            index1 = species_ids.index(d1_species[0])
            index2 = species_ids.index(d2_species[0])
            if index1 > index2:
                index1, index2 = index2, index1
            indices1.append(index1)
            indices2.append(index2)
            divergence_times.append(edge.get_termination())
            if starbeast_format:
                internal_pop_sizes.append(edge.get_pop_size(generations_per_branch_length_unit))
            else:
                internal_pop_sizes.append(Ne)
    # Do the same for the root.
    d1_species = []
    d2_species = []
    for edge in species_tree.get_edges():
        if edge.get_is_root_edge():
            if d1_species == []:
                d1_species = edge.get_extant_progeny_ids()
            elif d2_species == []:
                d2_species = edge.get_extant_progeny_ids()
                divergence_times.append(edge.get_origin())
                if starbeast_format:
                    internal_pop_sizes.append(starbeast_root_pop_size)
                else:
                    internal_pop_sizes.append(Ne)
                break
    d1_species.sort()
    d2_species.sort()
    index1 = species_ids.index(d1_species[0])
    index2 = species_ids.index(d2_species[0])
    if index1 > index2:
        index1, index2 = index2, index1
    indices1.append(index1)
    indices2.append(index2)

    # Sort the arrays indices1, indices2, divergence_times, and population sizes according to divergence_time.
    lists_are_sorted = False
    while lists_are_sorted == False:
        lists_are_sorted = True
        for x in range(len(indices1)-1):
            if divergence_times[x] > divergence_times[x+1]:
                indices1[x], indices1[x+1] = indices1[x+1], indices1[x]
                indices2[x], indices2[x+1] = indices2[x+1], indices2[x]
                divergence_times[x], divergence_times[x+1] = divergence_times[x+1], divergence_times[x]
                internal_pop_sizes[x], internal_pop_sizes[x+1] = internal_pop_sizes[x+1], internal_pop_sizes[x]
                lists_are_sorted = False
                break

    # Define the species/population tree for msprime.
    population_configurations = []
    if starbeast_format:
        for x in range(species_tree.get_number_of_extant_edges()):
            population_configurations.append(msprime.PopulationConfiguration(sample_size=sample_size, initial_size=terminal_pop_sizes[x]))
        demographic_events = []
        for x in range(len(indices1)-1):
            demographic_events.append(msprime.MassMigration(time=divergence_times[x]*generations_per_branch_length_unit, source=indices2[x], destination=indices1[x], proportion=1.0))
            demographic_events.append(msprime.PopulationParametersChange(time=divergence_times[x]*generations_per_branch_length_unit, initial_size=internal_pop_sizes[x], population_id=indices1[x]))
        demographic_events.append(msprime.MassMigration(time=divergence_times[-1]*generations_per_branch_length_unit, source=indices2[-1], destination=indices1[-1], proportion=1.0))
        demographic_events.append(msprime.PopulationParametersChange(time=divergence_times[-1]*generations_per_branch_length_unit, initial_size=internal_pop_sizes[-1], population_id=indices1[-1]))
    else:
        for _ in range(species_tree.get_number_of_extant_edges()):
            population_configurations.append(msprime.PopulationConfiguration(sample_size=sample_size, initial_size=Ne))
        demographic_events = []
        for x in range(len(indices1)-1):
            demographic_events.append(msprime.MassMigration(time=divergence_times[x]*generations_per_branch_length_unit, source=indices2[x], destination=indices1[x], proportion=1.0))
        demographic_events.append(msprime.MassMigration(time=divergence_times[-1]*generations_per_branch_length_unit, source=indices2[-1], destination=indices1[-1], proportion=1.0))

    # Return a tuple of population_configurations and demographic_events.
    return (population_configurations, demographic_events)
