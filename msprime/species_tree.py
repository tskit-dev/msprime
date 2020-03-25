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
import newick
import re

import msprime


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


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
    The newick library (https://pypi.org/project/newick/) is used for
    tree parsing. The tree string can end with a semi-colon, but this is
    not required.

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

    # Parse the tree with the newick library.
    root = parse_newick(tree, generations_per_branch_length_unit)

    # Define population configuratons and demographic events according to the
    # specified population size and the divergence times in the species tree.
    # Per divergence event (node in the tree), a mass migration with a proportion
    # of 1 of the population is used. The destination is the left-most leaf for
    # each node. Because we are using the n leaf populations, we map each node back
    # to the leaf population that it corresponds to.
    population_configurations = []
    demographic_events = []
    leaf_map = {}
    for node in root.walk("postorder"):
        if len(node.descendants) == 0:
            # Per extant species (= leaf node) in the tree, add a population with
            # size Ne. Species names are stored as metadata with the "species_name"
            # tag.
            population_configurations.append(msprime.PopulationConfiguration(
                    initial_size=Ne, metadata={"species_name": node.name.strip()}))
            leaf_map[node] = len(population_configurations) - 1
        else:
            # Per internal node, add one (if the node is bifurcating) or multiple
            # (if the node is multi-furcating) MassMigrations. The parent species
            # maps implicitly to the left-most child species, so we don't generate
            # any MassMigrations for that.
            leaf_map[node] = leaf_map[node.descendants[0]]
            # For each child species after the left-most one, we create
            # a MassMigration into the left-most species.
            for child in node.descendants[1:]:
                demographic_events.append(
                    msprime.MassMigration(
                        source=leaf_map[child], dest=leaf_map[node], time=node.time))

    # Sort demographic events by time.
    demographic_events.sort(key=lambda de: de.time)

    return population_configurations, demographic_events


def parse_starbeast(
        tree=None,
        branch_length_units="myr",
        generation_time=None):
    """
    Method to parse species trees produced by the program TreeAnnotator
    based on posterior tree distributions generated with StarBEAST
    (https://academic.oup.com/mbe/article/34/8/2101/3738283).
    Species trees produced by TreeAnnotator are written in Nexus
    (https://en.wikipedia.org/wiki/Nexus_file) format.

    Species trees written by TreeAnnotator are rooted, bi-furcating, and
    ultrametric. Branch lengths usually are in units of millions
    of years ("myr"), but the use of other units is permitted by
    StarBEAST (and thus TreeAnnotator). Here, we allow only units of
    millions of years ("myr") or years ("yr"). Leafs must be named and
    the tree must include information on population sizes of leaf and
    ancestral species in the form of annotation with the "dmv" tag,
    which is the case for trees written by StarBEAST and TreeAnnotator.
    However, to convert these "dmv" values into absolute population
    sizes, a generation time is required and must be provided separately
    to the tree.
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

    # Make sure that the generation time is specified and positive.
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

    # From the Nexus format definition (Maddison et al. 1997):
    # "For the most part, whitespace, inluding newline characters, is ignored,
    # with two exceptions: (1) whitespace indicates boundaries between words;
    # (2) in interleaved matrices, newline characters indicate the boundary
    # between data of different taxa."
    # As we do not parse matrices (we're only interested in taxa and trees
    # blocks), we ignore (2), replace newline characters with spaces and
    # replace multiple whitespaces with a single one.
    nexus_string = tree.replace("\n", " ")
    nexus_string = ' '.join(nexus_string.split())

    # From the Nexus format definition (Maddison et al. 1997):
    # "Commands or subcommands that differ only in case are homonymous."
    # We turn the whole nexus string into lowercase.
    nexus_string = nexus_string.lower()

    # Make sure that the string is in Nexus format.
    if nexus_string[0:6] != '#nexus':
        raise ValueError("The species tree does not appear to be in Nexus format.")

    # From the Nexus format definition (Maddison et al. 1997):
    # "Blocks are series of commands, beginning with a Begin command and ending
    # with an End command."
    # As semicolons are used only to terminate commands, potentially present
    # whitespace before semicolons has no meaning; we remove it for easier
    # parsing.
    # Then we identify the trees block and raise a ValueError if none is found.
    # This could be done with a regexp instead.
    nexus_string = nexus_string.replace(" ;", ";")
    tree_block_string = ""
    in_tree_block = False
    for x in range(len(nexus_string)):
        if nexus_string[x:x+12] == "begin trees;":
            in_tree_block = True
            tree_block_string += nexus_string[x]
        elif in_tree_block and nexus_string[x:x+4] == "end;":
            tree_block_string += nexus_string[x:x+4]
            break
        elif in_tree_block:
            tree_block_string += nexus_string[x]
    if tree_block_string == "" or tree_block_string[-4:] != "end;":
        raise ValueError("The Nexus string does not include a complete trees block.")

    # From the Nexus format definition (Maddison et al. 1997):
    # "Commands follow a simple format: the first token in the command
    # is the command name, which is followed by a series of tokens and
    # whitespace; the command is terminated by a semicolon."
    # Get the commands from the tree block, ignoring the begin and end
    # statements of the block.
    tree_block_commands = tree_block_string.split(";")
    tree_block_commands = [c.strip() for c in tree_block_commands]
    assert tree_block_commands[0] == "begin trees", "Tree block malformed"
    assert tree_block_commands[-1] == "", "Tree block malformed"
    assert tree_block_commands[-2] == "end", "Tree block malformed"
    assert len(tree_block_commands) > 3, "Tree block malformed"
    tree_block_commands = tree_block_commands[1:-2]

    # Ensure that exactly one of the commands is a translate command and
    # exactly one is a tree command, which is the case when the Nexus file
    # is written with TreeAnnotator based on a posterior tree distribution
    # generated with StarBEAST.
    translate_commands = []
    tree_commands = []
    for command in tree_block_commands:
        command_list = command.split()
        if command_list[0] == "translate":
            translate_commands.append(command)
        elif command_list[0] == "tree":
            tree_commands.append(command)
    if len(translate_commands) != 1:
        err = "The Nexus string does not contain exactly one translate command."
        raise ValueError(err)
    if len(tree_commands) != 1:
        err = "The Nexus string does not contain exactly one tree command."
        raise ValueError(err)
    translate_command = translate_commands[0]
    tree_command = tree_commands[0]
    assert "(" in tree_command, "No parentheses in tree string"

    # Skip forward until the newick string starts.
    tree_string = tree_command[tree_command.find("("):]

    # Make sure that annotation is included and has the dmv tag, which specifies
    # StarBEAST's mean estimate of the population size scaled by the number of
    # generation per time unit.
    assert "[" in tree_string, "Could not read annotation"
    assert "]" in tree_string, "Could not read annotation"
    assert tree_string.count("[") == tree_string.count("]"), "Could not read annotation"
    assert "&dmv={" in tree_string, "Could not find dmv tag in annotation"
    assert "}" in tree_string, "Could not read dmv annotation"
    assert tree_string.count("{") == tree_string.count("}"), "Could not read annotation"

    # Make sure that each substring that begins with an opening square bracket and ends
    # with a closing square bracket does not contain any further square or round brackets
    # in it and that it does include the dmv tag.
    in_annotation = False
    annotation_string = ""
    for x in range(len(tree_string)):
        if tree_string[x] == "[":
            in_annotation = True
            annotation_string += tree_string[x]
        elif tree_string[x] == "]":
            in_annotation = False
            annotation_string += tree_string[x]
            assert "[" not in annotation_string[1:-1], "Square bracket in annotation"
            assert "]" not in annotation_string[1:-1], "Square bracket in annotation"
            assert annotation_string.count("&dmv=") == 1, "Multiple or no dmv tags"
            # Make sure that the dmv tag is followed by a closing curly bracket.
            # Also ensure that the dmv tag is the first in the annotation.
            dmv_string = ""
            in_dmv = False
            for y in range(len(annotation_string)):
                dmv_string += annotation_string[y]
                if annotation_string[y] == "}":
                    break
            err = "Uneven curly parentheses in dmv annotation"
            assert dmv_string.count("{") == dmv_string.count("}"), err
            assert dmv_string.count("{") == 1, "Multiple or no values in dmv annotation"
            annotation_string = ""
            # Make sure that a positive number is found between curly brackets.
            clean_dmv_string = dmv_string.split("{")[1][0:-1]
            assert is_number(clean_dmv_string), "dmv annotation is not a number"
            assert float(clean_dmv_string) >= 0, "dmv annotation is negative"
        elif in_annotation:
            annotation_string += tree_string[x]

    # Because the newick module doesn't support parsing extended newick attributes
    # in general, we have to clean thing up manually before parsing the tree. Here,
    # we get rid of all annotations except for dmv.
    clean_tree_string = ""
    in_annotation = False
    in_dmv = False
    for x in range(len(tree_string)):
        if tree_string[x] == "[":
            in_annotation = True
            clean_tree_string += tree_string[x]
        elif tree_string[x] == "]":
            in_annotation = False
            clean_tree_string += tree_string[x]
        elif in_annotation:
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

    # Make sure that only dmv annotation remains in the tree string.
    in_annotation = False
    annotation_string = ""
    for x in range(len(clean_tree_string)):
        if clean_tree_string[x] == "[":
            in_annotation = True
            annotation_string += clean_tree_string[x]
        elif clean_tree_string[x] == "]":
            in_annotation = False
            annotation_string += clean_tree_string[x]
            assert annotation_string[0:7] == "[&dmv={", "Annotation could not be read"
            assert annotation_string[-2:] == "}]", "Annotation could not be read"
            assert is_number(annotation_string[7:-2]), "dmv annotation is not a number"
            assert float(annotation_string[7:-2]) >= 0, "dmv annotation is negative"
            annotation_string = ""
        elif in_annotation:
            annotation_string += clean_tree_string[x]

    # First, we trim the "translate" tag from the beginning.
    assert translate_command[0:10] == "translate "
    translate_command = translate_command[10:]

    # Use the translation block to back-translate species IDs in
    # the tree string.
    transl_originals = []
    transl_replaceds = []
    translation_list = translate_command.split(",")
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
        find_str1 = ',{}['.format(transl_replaceds[x])
        replace_str1 = ',{}['.format(transl_originals[x])
        find_str2 = '({}['.format(transl_replaceds[x])
        replace_str2 = '({}['.format(transl_originals[x])
        # Make sure that exactly one of the two find patterns is found.
        count1 = work_string.count(find_str1)
        count2 = work_string.count(find_str2)
        assert count1 + count2 == 1, "Ambiguous back-translation"
        if count1 == 1:
            work_string = work_string.replace(find_str1, replace_str1)
        else:
            work_string = work_string.replace(find_str2, replace_str2)
    clean_tree_string = work_string

    root = parse_newick(clean_tree_string, generations_per_branch_length_unit)

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
            # Per extant species (= leaf node) in the tree, add a population with
            # size pop_size. Species names are stored as metadata with the "species_name"
            # tag.
            metadata = {"species_name": node.name.strip()}
            population_configurations.append(
                msprime.PopulationConfiguration(initial_size=pop_size,
                                                metadata=metadata))
            leaf_map[node] = len(population_configurations) - 1
        else:
            # Per internal node, add one (if the node is bifurcating) or multiple
            # (if the node is multi-furcating) MassMigrations. The parent species
            # maps implicitly to the left-most child species, so we don't generate
            # any MassMigrations for that.
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
