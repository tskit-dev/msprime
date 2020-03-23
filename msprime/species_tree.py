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


def parse_newick(tree, branch_length_multiplier):
    """
    Parses the newick tree and annotates the resulting nodes with their
    time values, appropriately scaled.
    """
    # Parse the newick tree string.
    parsed = newick.loads(tree)
    if len(parsed) == 0:
        raise ValueError(f"Not a valid newick tree: '{tree}'")
    root = parsed[0]

    # Set node depths (distances from root).
    stack = [(root, 0)]
    num_nodes = 0
    max_depth = 0
    while len(stack) > 0:
        node, depth = stack.pop()
        if depth > max_depth:
            max_depth = depth
        num_nodes += 1
        node.depth = depth
        for child in node.descendants:
            stack.append((child, depth + child.length))
    if num_nodes < 3:
        raise ValueError("Newick tree must have at least three nodes")

    # Set node times (distances from present).
    for node in root.walk():
        node.time = (max_depth - node.depth) * branch_length_multiplier

    return root


def parse_species_tree(
        tree=None,
        Ne=None,
        branch_length_units="gen",
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
    if not isinstance(tree, str):
        raise TypeError("The species tree must be a string")

    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = 'The specified units for branch lengths ('
        err += '"{}") are not accepted. '.format(branch_length_units)
        err += 'Accepted units are "myr" (millions of years), "yr" (years), '
        err += 'and "gen" (generations).'
        raise ValueError(err)

    # Make sure that Ne is specified and positive.
    if Ne is None:
        raise ValueError("Ne should be specified.")
    else:
        try:
            Ne = float(Ne)
        except ValueError:
            raise ValueError("Population size Ne must be numeric.")
        if Ne <= 0:
            raise ValueError("Population size Ne must be > 0.")

    # Make sure that the generation time is either None or positive.
    if generation_time is not None:
        try:
            generation_time = float(generation_time)
        except ValueError:
            raise ValueError("Generation time must be numeric.")
        if generation_time <= 0:
            raise ValueError("Generation time must be > 0.")

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

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time)

    root = parse_newick(tree, generations_per_branch_length_unit)

    population_configurations = []
    demographic_events = []
    # The destination is the left-most leaf for each node. Because we are
    # using the n leaf populations, we map each node back to the leaf
    # population that it corresponds to.
    leaf_map = {}
    for node in root.walk("postorder"):
        if len(node.descendants) == 0:
            population_configurations.append(
                msprime.PopulationConfiguration(initial_size=Ne, metadata={"species_name": node.name.strip()}))
            leaf_map[node] = len(population_configurations) - 1
        else:
            # The parent species maps implicitly to the left-most child
            # species, so we don't generate any MassMigrations for that.
            leaf_map[node] = leaf_map[node.descendants[0]]
            # For each child species after the left-most one, we create
            # a MassMigration into the left-most species.
            for child in node.descendants[1:]:
                demographic_events.append(
                    msprime.MassMigration(
                        source=leaf_map[child], dest=leaf_map[node], time=node.time))

    demographic_events.sort(key=lambda de: de.time)

    return population_configurations, demographic_events


def parse_starbeast(
        tree=None,
        branch_length_units="myr",
        generation_time=None):
    """
    Method to parse species trees produced by the program StarBEAST,
    (https://academic.oup.com/mbe/article/34/8/2101/3738283)
    written in Nexus (https://en.wikipedia.org/wiki/Nexus_file) format.

    Species trees written by StarBEAST are rooted, bi-furcating, and
    ultrametric. Branch lengths usually are in units of millions
    of years ("myr"), but the use of years ("yr") as units is also
    possible. Leafs must be named. However, to convert these estimates to
    absolute population sizes, a generation time is required.
    """

    # Make sure a species tree is specified.
    if type(tree) is not str:
        raise ValueError("A species tree must be specified.")

    # Make sure that branch length units are either "myr" or "yr".
    allowed_branch_lenth_units = ["myr", "yr"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = 'The specified units for branch lengths ('
        err += '"{}") are not accepted. '.format(branch_length_units)
        err += 'Accepted units are "myr" (millions of years) or "yr" (years).'
        raise ValueError(err)

    # Make sure that the generation time is positive.
    if generation_time is None:
        raise TypeError("Generation time must not be None.")
    try:
        generation_time = float(generation_time)
    except ValueError:
        raise ValueError("Generation time must be numeric.")
    if generation_time <= 0:
        raise ValueError("Generation time must be > 0.")

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time)

    # Read the input file.
    tree_lines = tree.splitlines(False)
    if tree_lines[0][0:6].lower() != '#nexus':
        raise ValueError("The species tree does not appear to be in Nexus format")
    in_translation = False
    in_tree = False
    translate_string = ""
    for line in tree_lines:
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
                break

    # Skip forward until the newick string starts
    tree_string = clean_line[clean_line.find("("):]

    # TODO these need to be made into ValueErrors and tested.

    # Make sure the annotation can be read.
    assert '[' in tree_string, "Could not read annotation."
    assert ']' in tree_string, "Could not read annotation."
    assert "dmv=" in tree_string, "Could not find dmv tag in annotation."

    # Because the newick module doesn't support parsing extended newick attributes
    # in general, we have to clean thing up manually before parsing the tree. Here,
    # we get rid of all annotations except for dmv.
    clean_tree_string = ''
    in_comment = False
    in_dmv = False
    for x in range(len(tree_string)):
        if tree_string[x] == '[':
            in_comment = True
            clean_tree_string += tree_string[x]
        elif tree_string[x] == ']':
            in_comment = False
            clean_tree_string += tree_string[x]
        elif in_comment:
            if tree_string[x-1] == "[" and tree_string[x] == "&":
                clean_tree_string += "&"
            if tree_string[x-5:x] == "dmv={":
                in_dmv = True
                clean_tree_string += "dmv={"
            if in_dmv:
                clean_tree_string += tree_string[x]
            if tree_string[x] == "}":
                in_dmv = False
        else:
            clean_tree_string += tree_string[x]

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
            err += 'translation block.'
            raise ValueError(err)

        # Backtranslate to species names.
        work_string = clean_tree_string
        for x in range(len(transl_originals)):
            find_str = ',{}['.format(transl_replaceds[x])
            replace_str = ',{}['.format(transl_originals[x])
            work_string = work_string.replace(find_str, replace_str)
            find_str = '({}['.format(transl_replaceds[x])
            replace_str = '({}['.format(transl_originals[x])
            work_string = work_string.replace(find_str, replace_str)
        tree_string = work_string

    root = parse_newick(tree_string, generations_per_branch_length_unit)

    population_configurations = []
    demographic_events = []
    # The destination is the left-most leaf for each node. Because we are
    # using the n leaf populations, we map each node back to the leaf
    # population that it corresponds to.
    leaf_map = {}
    for node in root.walk("postorder"):
        find_pattern = '\\&dmv=\\{([\\d\\.]+?)\\}'
        dmv_patterns = re.search(find_pattern, node.name)
        if dmv_patterns is None:
            raise ValueError("No dmv annotation for node")
        pop_size = float(dmv_patterns.group(1)) * generations_per_branch_length_unit

        if len(node.descendants) == 0:
            population_configurations.append(
                msprime.PopulationConfiguration(initial_size=pop_size))
            leaf_map[node] = len(population_configurations) - 1
        else:
            # The parent species maps implicitly to the left-most child
            # species, so we don't generate any MassMigrations for that.
            leaf_map[node] = leaf_map[node.descendants[0]]
            # For each child species after the left-most one, we create
            # a MassMigration into the left-most species.
            for child in node.descendants[1:]:
                demographic_events.append(
                    msprime.MassMigration(
                        source=leaf_map[child], dest=leaf_map[node], time=node.time))
            demographic_events.append(
                msprime.PopulationParametersChange(
                    node.time,
                    initial_size=pop_size,
                    population_id=leaf_map[node]))
    demographic_events.sort(key=lambda de: de.time)
    return population_configurations, demographic_events
