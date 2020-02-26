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

import msprime


def get_generations_per_branch_length_unit(
        branch_length_units=None,
        generation_time=None):
    """
    Method to calculate the number of generations per branch length
    unit, given the branch length unit and a generation time.
    """
    if branch_length_units == "gen":
        generations_per_branch_length_unit = 1
    elif branch_length_units == "myr":
        generations_per_branch_length_unit = 10**6/generation_time
    else:
        generations_per_branch_length_unit = 1/generation_time
    return generations_per_branch_length_unit


def get_species_tree_string(species_tree_line=None):
    """
    Method to check the Newick tree string.
    """
    tree_patterns = re.search('\\(.+\\)', species_tree_line)
    # Make sure the pattern is found.
    if tree_patterns is None:
        raise ValueError("No species tree string found.")
    species_tree_string = tree_patterns.group(0)
    # Make sure that the numbers of opening and closing parentheses are
    # identical.
    if species_tree_string.count('(') != species_tree_string.count(')'):
        err = "The numbers of opening and closing parentheses in the species "
        err += "tree string must be identical."
        raise ValueError(err)
    # Make sure that between the beginning and the end of the tree string,
    # there are always more parentheses opened than closed.
    parentheses_open = 0
    for char in species_tree_string[:-2]:
        if char == "(":
            parentheses_open += 1
        elif char == ")":
            parentheses_open -= 1
        if parentheses_open <= 0:
            raise ValueError("The string seems to include multiple trees.")
    return species_tree_string


def parse_species_tree(
        species_tree=None,
        branch_length_units="gen",
        Ne=None,
        generation_time=None):
    """
    Method to parse species trees in Newick
    (https://en.wikipedia.org/wiki/Newick_format) format.

    Trees are assumed to be rooted and ultrametric and branch lengths
    must be included and correspond to time, either in units of millions
    of years ("myr"), years ("yr"), or generations ("gen"; default).
    Leafs must be named. An example for an accepted tree string in
    Newick format is:
    (((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)
    The tree string can end with a semi-colon, but this is not required.

    - An estimate of the effective population size Ne should be
        specified.
    - If and only if the branch lengths are not in units of
        generations, the generation time should be specified.
    """

    # Make sure a species tree is specified.
    if species_tree is None:
        raise ValueError("A species tree must be specified.")

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
        if generation_time is not None:
            err = 'With branch lengths in units of generations ("gen"), '
            err += 'a generation time should not be specified additionally.'
            raise ValueError(err)
    else:
        if generation_time is None:
            err = 'With branch lengths in units of '
            err += '"{}", a generation time must be '.format(branch_length_units)
            err += 'specified additionally.'
            raise ValueError(err)

    # Make sure that a population size is specified.
    if Ne is None:
        raise ValueError("Ne should be specified.")

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time
        )

    # Read the input file.
    species_tree_lines = species_tree.splitlines(False)
    if len(species_tree_lines) > 1:
        raise ValueError("The species tree has multiple lines.")
    species_tree_string = get_species_tree_string(species_tree_lines[0])

    # Parse the newick tree string.
    root = newick.loads(species_tree_string)[0]

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

    # Get a list of species IDs along with terminal population sizes.
    species_ids = []
    terminal_nes = []
    for leaf in root.get_leaves():
        species_ids.append(leaf.name)
        terminal_nes.append(Ne)
    s = sorted(zip(species_ids, terminal_nes))
    species_ids, terminal_nes = map(list, zip(*s))

    # Determine at which time which populations should merge.
    indices1 = []
    indices2 = []
    internal_nes = []
    divergence_times = []
    for node in root.walk():
        if node is not root and node.is_leaf is False:
            d1_species = node.descendants[0].get_leaf_names()
            d2_species = node.descendants[1].get_leaf_names()
            index1 = species_ids.index(sorted(d1_species)[0])
            index2 = species_ids.index(sorted(d2_species)[0])
            if index1 > index2:
                index1, index2 = index2, index1
            indices1.append(index1)
            indices2.append(index2)
            divergence_times.append(node.height)
            internal_nes.append(Ne)

    # Add information for the root.
    d1_species = root.descendants[0].get_leaf_names()
    d2_species = root.descendants[1].get_leaf_names()
    index1 = species_ids.index(sorted(d1_species)[0])
    index2 = species_ids.index(sorted(d2_species)[0])
    if index1 > index2:
        index1, index2 = index2, index1
    indices1.append(index1)
    indices2.append(index2)
    divergence_times.append(root.height)
    internal_nes.append(Ne)

    # Sort the arrays indices1, indices2, divergence_times, and
    # population sizes according to divergence_time.
    s = sorted(zip(divergence_times, indices1, indices2, internal_nes))
    divergence_times, indices1, indices2, internal_nes = map(list, zip(*s))

    # Define the species/population tree for msprime.
    population_configurations = []
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


def parse_starbeast(
        species_tree=None,
        branch_length_units="myr",
        generation_time=None):
    """
    Method to parse species trees produced by the program StarBEAST,
    (https://academic.oup.com/mbe/article/34/8/2101/3738283)
    written in Nexus (https://en.wikipedia.org/wiki/Nexus_file) format.

    Species trees written by StarBEAST are rooted, bi-furcating, and
    ultrametric. Branch lengths usually are in units of millions
    of years ("myr"), but the use of years ("yr") as units is also
    possible. Leafs must be named. The population size is not required
    as input since StarBEAST stores estimates for population sizes
    in branch annotation. However, to convert these estimates to
    absolute population sizes, a generation time is required.
    """

    # Make sure a species tree is specified.
    if species_tree is None:
        raise ValueError("A species tree must be specified.")

    # Make sure that branch length units are either "myr" or "yr".
    allowed_branch_lenth_units = ["myr", "yr"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = 'The specified units for branch lengths ('
        err += '"{}") are not accepted. '.format(branch_length_units)
        err += 'Accepted units are "myr" (millions of years) or "yr" (years).'
        raise ValueError(err)

    # Make sure that the generation time is positive.
    if generation_time <= 0:
        raise ValueError("Generation time must be > 0")

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time
        )

    # Read the input file.
    species_tree_lines = species_tree.splitlines(False)
    if species_tree_lines[0][0:6].lower() != '#nexus':
        raise ValueError("The species tree does not appear to be in Nexus format")
    in_tree = False
    translate_string = ""
    for line in species_tree_lines:
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
                species_tree_string = get_species_tree_string(clean_line)
                break

    # Make sure the annotation can be read.
    if '[' not in species_tree_string or ']' not in species_tree_string:
        raise ValueError("Could not read annotation from species tree")
    if "dmv=" not in species_tree_string:
        raise ValueError("Could not find dmv tag in annotation.")

    # Remove all annotation except for dmv.
    clean_species_tree_string = ''
    in_comment = False
    in_dmv = False
    for x in range(len(species_tree_string)):
        if species_tree_string[x] == '[':
            in_comment = True
            clean_species_tree_string += species_tree_string[x]
        elif species_tree_string[x] == ']':
            in_comment = False
            clean_species_tree_string += species_tree_string[x]
        elif in_comment:
            if species_tree_string[x-1] == "[" and species_tree_string[x] == "&":
                clean_species_tree_string += "&"
            if species_tree_string[x-5:x] == "dmv={":
                in_dmv = True
                clean_species_tree_string += "dmv={"
            if in_dmv:
                clean_species_tree_string += species_tree_string[x]
            if species_tree_string[x] == "}":
                in_dmv = False
        else:
            clean_species_tree_string += species_tree_string[x]

    # Get the population size at the root.
    find_pattern = '(\\(.+\\))\\[\\&dmv=\\{([\\d\\.]+?)\\}\\]'
    tree_patterns = re.search(find_pattern, clean_species_tree_string)
    starbeast_root_pop_size = float(tree_patterns.group(2))
    starbeast_root_pop_size *= generations_per_branch_length_unit

    # Use the translation block (if present) to
    # back-translate species IDs in the tree string.
    if translate_string != "":
        transl_originals = []
        transl_replaceds = []
        translation_list = translate_string.split(",")
        for item in translation_list:
            item_list = item.split()
            assert len(item_list) > 1, "The translation block is malformed."
            if len(item_list) != 2:
                err = "Species IDs in the translation block appear to include "
                err += "whitespace. This is not supported."
                raise ValueError(err)
            transl_originals.append(item_list[1])
            transl_replaceds.append(item_list[0])

        # Make sure both the species names and their translations
        # are unique.
        if len(transl_originals) != len(set(transl_originals)):
            err = 'One or more species names in the translation '
            err += 'block are duplicated.'
            raise ValueError(err)
        if len(transl_replaceds) != len(set(transl_replaceds)):
            err = 'One or more translations in the translation '
            err += 'block are duplicated.'
            raise ValueError(err)

        # Make sure no species names are identical with any
        # translations.
        transl_orireps = transl_originals + transl_replaceds
        if len(transl_orireps) != len(set(transl_orireps)):
            err = 'One or more of the species names are identical '
            err += 'to one or more of the translations in the '
            err += 'translation block of file {}.'.format(filename)
            raise ValueError(err)

        # Backtranslate to species names.
        work_string = clean_species_tree_string
        for x in range(len(transl_originals)):
            find_str = ',{}['.format(transl_replaceds[x])
            replace_str = ',{}['.format(transl_originals[x])
            work_string = work_string.replace(find_str, replace_str)
            find_str = '({}['.format(transl_replaceds[x])
            replace_str = '({}['.format(transl_originals[x])
            work_string = work_string.replace(find_str, replace_str)
        species_tree_string = work_string

    # Parse the newick tree string.
    root = newick.loads(species_tree_string)[0]

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

    # Get a list of species IDs along with terminal population sizes.
    species_ids = []
    terminal_nes = []
    for leaf in root.get_leaves():
        species_id = leaf.name.split("[")[0]
        species_ids.append(species_id)
        edge_dmv = float(leaf.name.split("{")[1].split("}")[0])
        assert edge_dmv > 0, 'Parsed Ne is zero.'
        edge_ne = edge_dmv * generations_per_branch_length_unit
        terminal_nes.append(edge_ne)
    s = sorted(zip(species_ids, terminal_nes))
    species_ids, terminal_nes = map(list, zip(*s))

    # Determine at which time which populations should merge.
    indices1 = []
    indices2 = []
    internal_nes = []
    divergence_times = []
    for node in root.walk():
        if node is not root and node.is_leaf is False:
            d1_species = [l.split("[")[0] for l in node.descendants[0].get_leaf_names()]
            d2_species = [l.split("[")[0] for l in node.descendants[1].get_leaf_names()]
            index1 = species_ids.index(sorted(d1_species)[0])
            index2 = species_ids.index(sorted(d2_species)[0])
            if index1 > index2:
                index1, index2 = index2, index1
            indices1.append(index1)
            indices2.append(index2)
            divergence_times.append(node.height)
            edge_dmv = float(node.name.split("{")[1].split("}")[0])
            assert edge_dmv > 0, 'Parsed Ne is zero.'
            edge_ne = edge_dmv * generations_per_branch_length_unit
            internal_nes.append(edge_ne)

    # Add information for the root.
    d1_species = [l.split("[")[0] for l in root.descendants[0].get_leaf_names()]
    d2_species = [l.split("[")[0] for l in root.descendants[1].get_leaf_names()]
    index1 = species_ids.index(sorted(d1_species)[0])
    index2 = species_ids.index(sorted(d2_species)[0])
    if index1 > index2:
        index1, index2 = index2, index1
    indices1.append(index1)
    indices2.append(index2)
    divergence_times.append(root.height)
    internal_nes.append(starbeast_root_pop_size)

    # Sort the arrays indices1, indices2, divergence_times, and
    # population sizes according to divergence_time.
    s = sorted(zip(divergence_times, indices1, indices2, internal_nes))
    divergence_times, indices1, indices2, internal_nes = map(list, zip(*s))

    # Define the species/population tree for msprime.
    population_configurations = []
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

    # Return a tuple of population_configurations and demographic_events.
    return population_configurations, demographic_events
