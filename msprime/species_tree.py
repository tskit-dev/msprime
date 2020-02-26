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
import re
import newick
import sys # XXX Remove

import msprime


def parse_species_tree(
        filename=None,
        branch_length_units="gen",
        Ne=None,
        generation_time=None):
    """
    Method to parse species trees from files in Newick
    https://en.wikipedia.org/wiki/Newick_format format or Nexus
    (https://en.wikipedia.org/wiki/Nexus_file) format with embedded
    Newick strings.

    Trees must be rooted, binary, and ultrametric and branch lengths
    must be included and correspond to time, either in units of millions
    of years ("myr"), years ("yr"), or generations ("gen"; default).
    Leafs must be named. An example for an accepted tree string in
    Newick format is this:
    (((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)
    The tree string can end with a semi-colon, but this is not required.

    We allow the following input:
    1) filename must be specified.
    2) If and only if the tree is not in StarBEAST format,
        Ne should be specified.
    3) If and only if the branch lengths are not in units of
        generations, the generation time should be specified.

    As msprime does not name populations and instead only assigns
    numbers to populations (starting with 0), the species names stored
    in the species tree are not used when defining instances of
    PopulationConfiguration. However, these instances will be sorted
    according to the alphabetical order of the corresponding species
    names, allowing the user to link species from the species tree with
    populations in msprime.

    Introgression events, as stored in the extended Newick format by
    PhyloNet and possibly other programs, are not parsed and therefore
    not used.
    Rate matrices of continuous migration between co-existing lineages
    can not be stored in (extended) Newick or Nexus format and are
    therefore also not used here. However, these can be added to the
    model after population configurations and demographic events are
    defined with this function.
    """

    # Make sure a filename is specified.
    err = 'A filename must be specified when calling function parse_species_tree().'
    assert filename, err

    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = 'The specified units for branch lengths ('
        err += '"{}") are not accepted. '.format(branch_length_units)
        err += 'Accepted units are "myr" (millions of years), "yr" (years), '
        err += 'and "gen" (generations).'
        raise ValueError(err)

    # Make sure that the population size is either None or positive.
    if Ne is not None and Ne <= 0:
        raise ValueError("Population size Ne must be > 0")

    # Make sure that the generation time is either None or positive.
    if generation_time is not None and generation_time <= 0:
        raise ValueError("Generation time must be > 0")

    # Make sure that the generation time is specified if and only if
    # branch lengths are not in units of generations.
    if branch_length_units == "gen":
        err = 'With branch lengths in units of generations ("gen"), generation_time '
        err += 'should not be specified.'
        assert generation_time is None, err
    else:
        err = 'With branch lengths in units of '
        err += '"{}", generation_time must be specified.'.format(branch_length_units)
        assert generation_time, err

    # Get the number of generations per branch length unit.
    if branch_length_units == "gen":
        generations_per_branch_length_unit = 1
    elif branch_length_units == "myr":
        generations_per_branch_length_unit = 10**6/generation_time
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
            elif in_tree and line.strip() != '':
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
                                elif in_comment is False:
                                    tree_string_raw += letter
                    if starbeast_format:
                        find_pattern = '(\\(.+\\))\\[\\&dmv=\\{([\\d\\.]+?)\\}\\]'
                        tree_patterns = re.search(find_pattern, tree_string_raw)
                        tree_string = tree_patterns.group(1)
                        starbeast_root_pop_size = float(tree_patterns.group(2))
                        starbeast_root_pop_size *= generations_per_branch_length_unit
                    else:
                        tree_patterns = re.search('\\(.+\\)', tree_string_raw)
                        tree_string = tree_patterns.group(0)
                    # Use the translation block (if present) to
                    # back-translate species IDs in the tree string.
                    if translate_string != "":
                        transl_originals = []
                        transl_replaceds = []
                        translation_list = translate_string.split(",")
                        for item in translation_list:
                            item_list = item.split()
                            err = 'The translation block of file '
                            err += '{} is malformed.'.format(filename)
                            assert len(item_list) > 1, err
                            err = 'Species IDs in the translation '
                            err += 'block of file {} '.format(filename)
                            err += 'appear to include whitespace. This '
                            err += 'is not supported.'
                            assert len(item_list) == 2, err
                            transl_originals.append(item_list[1])
                            transl_replaceds.append(item_list[0])
                        err = 'One or more species names in the translation '
                        err += 'block of file {} are duplicated.'.format(filename)
                        assert len(transl_originals) == len(set(transl_originals)), err
                        err = 'One or more translations in the translation '
                        err += 'block of file {} are duplicated.'.format(filename)
                        assert len(transl_replaceds) == len(set(transl_replaceds)), err
                        transl_orireps = transl_originals + transl_replaceds
                        err = 'One or more of the species names are identical '
                        err += 'to one or more of the translations in the '
                        err += 'translation block of file {}.'.format(filename)
                        assert len(transl_orireps) == len(set(transl_orireps)), err
                        for x in range(len(transl_originals)):
                            find_str = ',{}['.format(transl_replaceds[x])
                            replace_str = ',{}['.format(transl_originals[x])
                            tree_string = tree_string.replace(find_str, replace_str)
                            find_str = '({}['.format(transl_replaceds[x])
                            replace_str = '({}['.format(transl_originals[x])
                            tree_string = tree_string.replace(find_str, replace_str)
                    break

    else:
        # Test if the input is in newick format.
        err = 'File {} does not appear to be in Nexus '.format(filename)
        err += 'format but contains more than one line.'
        assert len(inlines) == 1, err
        err = 'Unexpected tree format in file {}. '.format(filename)
        err += 'The first character should be "(".'
        assert inlines[0][0] == '(', err
        err = 'The number of opening and closing parentheses differ '
        err += 'in the tree string!'
        assert inlines[0].count('(') == inlines[0].count(')'), err
        tree_string_raw = inlines[0]
        tree_patterns = re.search('\\(.+\\)', tree_string_raw)
        tree_string = tree_patterns.group(0)

    # Make sure a tree string is found.
    assert tree_string, 'No tree could be found in file {}.'.format(filename)

    # Make sure that between the beginning and the end of the tree string,
    # there are always more parentheses opened than closed.
    parentheses_open = 0
    for char in tree_string[:-2]:
        if char == "(":
            parentheses_open += 1
        elif char == ")":
            parentheses_open -= 1
        err = 'The tree string in file {} '.format(filename)
        err += 'appears to be malformed (possibly the first and '
        err += 'last parentheses are missing).'
        assert parentheses_open > 0, err

    # Now that we know if the tree is in starbeast format or not, make
    # sure that Ne is specified if and only if the tree is not in
    # starbeast format.
    if starbeast_format:
        err = 'With species trees in StarBEAST format, Ne should not be specified.'
        assert Ne is None, err
    else:
        assert Ne, 'Ne should be specified.'

    # Parse the newick tree string.
    root = newick.loads(tree_string)[0]

    # Remove nodes that have only a single child if there should be any.
    root.remove_redundant_nodes()

    # Resolve polytomies by injecting zero-length branches.
    root.resolve_polytomies()

    # Set node depths (distances from root).
    max_depth = 0
    for node in root.walk():
        depth = 0
        moving_node = node
        while moving_node.ancestor:
            depth += moving_node.length
            moving_node = moving_node.ancestor
        node.depth = depth
        if depth > max_depth:
            max_depth = depth

    # Set node heights (distances from present).
    for node in root.walk():
        node.height = max_depth - node.depth

    # Get a list of species IDs along with population sizes.
    species_ids = []
    terminal_nes = []
    for leaf in root.get_leaves():
        if starbeast_format:
            species_id = leaf.name.split("[")[0]
            species_ids.append(species_id)
            edge_dmv = float(leaf.name.split("{")[1].split("}")[0])
            assert edge_dmv > 0, 'Parsed Ne is zero.'
            edge_ne = edge_dmv * generations_per_branch_length_unit
            terminal_nes.append(edge_ne)
        else:
            species_ids.append(leaf.name)
            terminal_nes.append(Ne)
    lists_are_sorted = False
    while lists_are_sorted is False:
        lists_are_sorted = True
        for x in range(len(species_ids)-1):
            if species_ids[x] > species_ids[x+1]:
                species_ids[x], species_ids[x+1] = species_ids[x+1], species_ids[x]
                terminal_nes[x], terminal_nes[x+1] = terminal_nes[x+1], terminal_nes[x]
                lists_are_sorted = False
                break

    # Determine at which time which populations should merge.
    indices1 = []
    indices2 = []
    internal_nes = []
    divergence_times = []
    for node in root.walk():
        if node.name and node.is_leaf is False:
            if starbeast_format:
                d1_species = [l.split("[")[0] for l in node.descendants[0].get_leaf_names()]
                d2_species = [l.split("[")[0] for l in node.descendants[1].get_leaf_names()]
            else:
                d1_species = node.descendants[0].get_leaf_names()
                d2_species = node.descendants[1].get_leaf_names()
            index1 = species_ids.index(sorted(d1_species)[0])
            index2 = species_ids.index(sorted(d2_species)[0])
            if index1 > index2:
                index1, index2 = index2, index1
            indices1.append(index1)
            indices2.append(index2)
            divergence_times.append(node.height)
            if starbeast_format:
                edge_dmv = float(node.name.split("{")[1].split("}")[0])
                assert edge_dmv > 0, 'Parsed Ne is zero.'
                edge_ne = edge_dmv * generations_per_branch_length_unit
                internal_nes.append(edge_ne)
            else:
                internal_nes.append(Ne)

    # Add information for the root.
    if starbeast_format:
        d1_species = [l.split("[")[0] for l in root.descendants[0].get_leaf_names()]
        d2_species = [l.split("[")[0] for l in root.descendants[1].get_leaf_names()]
    else:
        d1_species = root.descendants[0].get_leaf_names()
        d2_species = root.descendants[1].get_leaf_names()
    index1 = species_ids.index(sorted(d1_species)[0])
    index2 = species_ids.index(sorted(d2_species)[0])
    if index1 > index2:
        index1, index2 = index2, index1
    indices1.append(index1)
    indices2.append(index2)
    divergence_times.append(root.height)
    if starbeast_format:
        internal_nes.append(starbeast_root_pop_size)
    else:
        internal_nes.append(Ne)

    # Sort the arrays indices1, indices2, divergence_times, and
    # population sizes according to divergence_time.
    s = sorted(zip(divergence_times,indices1,indices2,internal_nes))
    divergence_times,indices1,indices2,internal_nes = map(list, zip(*s))

    # Define the species/population tree for msprime.
    population_configurations = []
    if starbeast_format:
        for x in range(len(root.get_leaves())):
            population_configurations.append(
                msprime.PopulationConfiguration(
                    initial_size=terminal_nes[x]))
        demographic_events = []
        for x in range(len(indices1)-1):
            demographic_events.append(
                msprime.MassMigration(
                    time=divergence_times[x]*generations_per_branch_length_unit,
                    source=indices2[x],
                    destination=indices1[x],
                    proportion=1.0))
            demographic_events.append(
                msprime.PopulationParametersChange(
                    time=divergence_times[x]*generations_per_branch_length_unit,
                    initial_size=internal_nes[x],
                    population_id=indices1[x]))
        demographic_events.append(
            msprime.MassMigration(
                time=divergence_times[-1]*generations_per_branch_length_unit,
                source=indices2[-1],
                destination=indices1[-1],
                proportion=1.0))
        demographic_events.append(
            msprime.PopulationParametersChange(
                time=divergence_times[-1]*generations_per_branch_length_unit,
                initial_size=internal_nes[-1],
                population_id=indices1[-1]))
    else:
        for _ in range(len(root.get_leaves())):
            population_configurations.append(
                msprime.PopulationConfiguration(
                    initial_size=Ne))
        demographic_events = []
        for x in range(len(indices1)-1):
            demographic_events.append(
                msprime.MassMigration(
                    time=divergence_times[x]*generations_per_branch_length_unit,
                    source=indices2[x],
                    destination=indices1[x],
                    proportion=1.0))
        demographic_events.append(
            msprime.MassMigration(
                time=divergence_times[-1]*generations_per_branch_length_unit,
                source=indices2[-1],
                destination=indices1[-1],
                proportion=1.0))

    # Return a tuple of population_configurations and demographic_events.
    return population_configurations, demographic_events
