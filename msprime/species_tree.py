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
        if ":" in working_newick_string:
            pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?):([\d\.Ee-]+?),([a-zA-Z0-9_\.\-]+?):([\d\.Ee-]+?)\)")
        else:
            print("ERROR: It appears that the tree string does not include branch lengths!")
            sys.exit(1)
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
        if hit_rooted == None:
            print('ERROR: The newick tree string could not be parsed!')
            print('The remaining unparsed part of the newick string: \"{}\"'.format(working_newick_string))
            sys.exit(1)
        else:
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
        if ":" in working_newick_string:
            if "[" in working_newick_string:
                pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?)\[\&dmv=\{([\d\.]+?)\}\]:([\d\.Ee-]+?),([a-zA-Z0-9_\.\-]+?)\[\&dmv=\{([\d\.]+?)\}\]:([\d\.Ee-]+?)\)")
            else:
                print("ERROR: It appears that the tree string does not include annotation!")
                sys.exit(1)
        else:
            print("ERROR: It appears that the tree string does not include branch lengths!")
            sys.exit(1)
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
        if hit_rooted == None:
            print('ERROR: The newick tree string could not be parsed!')
            print('The remaining unparsed part of the newick string: \"{}\"'.format(working_newick_string))
            sys.exit(1)
        else:
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
        # These should be approximately similar, if not produce a warning.
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
                        if parent_edge == None:
                            print('ERROR: The parent edge for edge {} could not be found'.format(this_edge.get_id()))
                            sys.exit(1)
                        total_edge_length += parent_edge.get_length()
                        if parent_edge.get_is_root_edge():
                            root_edge_found = True
                            total_edge_lengths.append(total_edge_length)
                        else:
                            this_edge = parent_edge
        max_total_edge_length = max(total_edge_lengths)
        if max_total_edge_length - min(total_edge_lengths) > 0.1:
            print('WARNING: The tree appears not to be ultrametric. Some terminal branches will be extended so that they all end at the same time.')
            print('')
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
                        if parent_edge == None:
                            print('ERROR: The parent edge for edge {} could not be found'.format(this_edge.get_id()))
                            sys.exit(1)
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
                    if parent_edge == None:
                        print('ERROR: The parent edge for edge {} could not be found'.format(this_edge.get_id()))
                        sys.exit(1)
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

    def get_pop_size(self, generation_time):
        return self.dmv * (1000000/float(generation_time))

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
        branch_lengths="myr",
        pop_size=None,
        generation_time=None,
        diploid=True):
    """
    Method to parse species trees from files in Newick format or Nexus format with embedded Newick strings.
    Trees must be ultrametric and branch lengths must correspond to time, either in units of millions of years
    ("myr"; default), years ("yr"), or generations ("gen"). XXXTODO: Implement options for branch lengths.
    """
    # Read the input file.
    species_tree = None
    starbeast_format = False
    infile = open(filename)
    inlines = infile.readlines()
    if inlines[0][0:6].lower() == '#nexus':
        # Assume the input is in nexus format. Maximally one tree string is read.
        in_tree = False
        for line in inlines:
            if line.strip().lower() == 'begin trees;':
                in_tree = True
            elif line.strip().lower() == 'end;':
                in_tree = False
            elif in_tree and line.strip() is not '':
                clean_line = line.strip()
                if clean_line[0:4].lower() == "tree":
                    if '[' in clean_line and ']' in clean_line:
                        if "dmv=" in clean_line:
                            # Turn on star-c-genie mode
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
                        root_pop_size = float(tree_patterns.group(2)) * (1000000/generation_time)
                    else:
                        tree_patterns = re.search('\(.+\)',tree_string_raw)
                        tree_string = tree_patterns.group(0)
                    species_tree = Tree(tree_string)
                    break
    else:
        # Test if the input is in newick format.
        if inlines[0][0] == '(':
            if inlines[0].count('(') != inlines[0].count(')'):
                print('ERROR: The number of opening and closing parentheses differ in the tree string!')
                sys.exit(1)
            else:
                tree_string_raw = inlines[0]
                tree_patterns = re.search('\(.+\)',tree_string_raw)
                tree_string = tree_patterns.group(0)
                species_tree = Tree(tree_string)
        else:
            print('ERROR: Unexpected tree format in file {}'.format(filename))

    # Make sure a tree is found.
    if species_tree is None:
        print('\nERROR: No tree could be found in file {}'.format(filename))
        sys.exit(1)

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

    # Get a list of species Ids.
    species_ids = []
    terminal_pop_sizes = []
    for edge in species_tree.get_edges():
        if edge.get_is_extant():
            species_ids.append(edge.get_node_ids()[1])
            if starbeast_format:
                terminal_pop_sizes.append(edge.get_pop_size(generation_time))
            else:
                terminal_pop_sizes.append(pop_size)
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
                internal_pop_sizes.append(edge.get_pop_size(generation_time))
            else:
                internal_pop_sizes.append(pop_size)
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
                    internal_pop_sizes.append(root_pop_size)
                else:
                    internal_pop_sizes.append(pop_size)
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

    # Feedback.
    if starbeast_format:
        print('      -> newick string annotation will be used to inform per-branch population sizes.')

    # Define the species/population tree for msprime.
    population_configurations = []
    if diploid:
        population_configurations_diploid = []
    if starbeast_format:
        for x in range(species_tree.get_number_of_extant_edges()):
            population_configurations.append(msprime.PopulationConfiguration(sample_size=1, initial_size=terminal_pop_sizes[x]))
            if diploid:
                population_configurations_diploid.append(msprime.PopulationConfiguration(sample_size=2, initial_size=terminal_pop_sizes[x]))
        demographic_events = []
        for x in range(len(indices1)-1):
            demographic_events.append(msprime.MassMigration(time=(divergence_times[x]*1000000)/generation_time, source=indices2[x], destination=indices1[x], proportion=1.0))
            demographic_events.append(msprime.PopulationParametersChange(time=(divergence_times[x]*1000000)/generation_time, initial_size=internal_pop_sizes[x], population_id=indices1[x]))
        demographic_events.append(msprime.MassMigration(time=(divergence_times[-1]*1000000)/generation_time, source=indices2[-1], destination=indices1[-1], proportion=1.0))
        demographic_events.append(msprime.PopulationParametersChange(time=(divergence_times[-1]*1000000)/generation_time, initial_size=internal_pop_sizes[-1], population_id=indices1[-1]))
    else:
        for _ in range(species_tree.get_number_of_extant_edges()):
            population_configurations.append(msprime.PopulationConfiguration(sample_size=1, initial_size=pop_size))
            if diploid:
                population_configurations_diploid.append(msprime.PopulationConfiguration(sample_size=2, initial_size=pop_size))
        demographic_events = []
        for x in range(len(indices1)-1):
            demographic_events.append(msprime.MassMigration(time=(divergence_times[x]*1000000)/generation_time, source=indices2[x], destination=indices1[x], proportion=1.0))
        demographic_events.append(msprime.MassMigration(time=(divergence_times[-1]*1000000)/generation_time, source=indices2[-1], destination=indices1[-1], proportion=1.0))

    # Return a tuple of population_configurations and demographic_events.
    return (population_configurations, demographic_events)
