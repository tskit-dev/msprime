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

from . import demography as demog


def parse_starbeast(tree, generation_time, branch_length_units="myr"):
    """
    This function parses a species tree produced by the program `TreeAnnotator
    <https://www.beast2.org/treeannotator>`_
    based on a posterior tree distribution generated with `StarBEAST
    <https://academic.oup.com/mbe/article/34/8/2101/3738283>`_  and defines a
    simulation model according to the species tree. Species trees produced by
    TreeAnnotator are written in `Nexus
    <https://en.wikipedia.org/wiki/Nexus_file>`_ format and are rooted,
    bifurcating, and ultrametric. Branch lengths usually are in units of millions
    of years, but the use of other units is permitted by StarBEAST (and thus
    TreeAnnotator). This function allows branch length units of millions of years
    or years. Leaves must be named and the tree must include information on
    population sizes of leaf and ancestral species in the form of annotation with
    the "dmv" tag, which is the case for trees written by TreeAnnotator based on
    StarBEAST posterior tree distributions.

    After reading the input tree, this function defines a
    :class:`.PopulationConfiguration` instance for each terminal node in the tree,
    corresponding to extant species. These population configurations store the
    species' name and population size, both according to information from the input
    tree. Additionally, a :class:`.MassMigration` instance is defined for each
    internal node, with the time of the mass migration set according to the age of
    the node in the species tree. For each internal node, the left one of the two
    descendant populations is arbitrarily selected as the destination in the mass
    migration defined for that node. A :class:`.PopulationParametersChange`
    instance is also added for each internal node to adjust the population
    size of the destination population according to the information given in the
    tree for the population size of the species that is ancestral to the node. Like
    the mass migration event defined for the same node, the time of the population
    parameter change is also set according to the age of the node.

    :param str tree: The tree string in Nexus format, with named leaves, branch
        lengths, and branch annotation. Typically, this string is the entire content
        of a file written by TreeAnnotator.
    :param float generation_time: The number of years per generation.
    :param str branch_length_units: The units of time in which the species tree's
        branch lengths are measured. Allowed branch length units are millions of
        years, and years; these should be specified with the strings ``"myr"`` or
        ``"yr"``, respectively. This defaults to ``"myr"``.
    :return: A tuple of two lists of which the first contains
        :class:`.PopulationConfiguration` instances and the second contains
        :class:`.MassMigration` and :class:`.PopulationParametersChange` instances.
        The population configurations specify the size of each population according
        to the information from the input species tree and the species name
        corresponding to each population. Species names are stored as metadata in
        each :class:`.PopulationConfiguration` instance, with the metadata tag
        "species_name". Sampling configurations and growth rates are not specified
        in the population configurations. The list of population configurations is
        ordered according to the order of the corresponding extant species in a
        `post-order tree traversal
        <https://en.wikipedia.org/wiki/Tree_traversal#Post-order_(LRN)>`_.
        The list of mass migration events and population parameter changes is
        ordered by the time of the events, from young to old events.
    :rtype: (list, list)
    :warning: This function does not modify migration matrices. When the population
        configurations and mass migration events returned by this function are used
        to simulate with the :func:`.simulate` function, it should be ensured that
        migration rates to source populations of mass migration events are zero
        after the mass migration (viewed backwards in time).
    """

    # Make sure that branch length units are either "myr" or "yr".
    allowed_branch_lenth_units = ["myr", "yr"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = "The specified units for branch lengths ("
        err += f'"{branch_length_units}") are not accepted. '
        err += 'Accepted units are "myr" (millions of years) or "yr" (years).'
        raise ValueError(err)

    generation_time = check_generation_time(generation_time)
    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time
    )

    translate_string, tree_string = parse_nexus(tree)
    species_name_map = parse_translate_command(translate_string)
    clean_tree_string = strip_extra_annotations(tree_string)
    return process_starbeast_tree(
        clean_tree_string, generations_per_branch_length_unit, species_name_map
    )


def parse_species_tree(tree, Ne, branch_length_units="gen", generation_time=None):
    """
    This function parses a species tree in
    `Newick <https://en.wikipedia.org/wiki/Newick_format>`_ format and defines a
    simulation model according to the species tree. The tree is assumed to be
    rooted and ultrametric and branch lengths must be included and correspond to
    time, either in units of millions of years, years, or generations. Leaves must
    be named.

    After reading the input tree, this function defines a
    :class:`.PopulationConfiguration` instance for each terminal node in the tree,
    corresponding to extant species. These population configurations store the
    species' name and population size. The specified Ne is used as the size of all
    populations. Additionally, one or more :class:`.MassMigration` instances are
    defined for each internal node, with the time of the mass migration set
    according to the age of the node in the species tree. :math:`n - 1` mass
    migration events are defined for internal nodes with `n` descendants, meaning
    that a single event is defined for bifurcating nodes. For each internal node,
    the left-most of the descendant populations is arbitrarily selected as the
    destination in all mass migrations defined for that node.

    :param str tree: The tree string in Newick format, with named leaves and branch
        lengths.
    :param float Ne: The effective population size.
    :param str branch_length_units: The units of time in which the species tree's
        branch lengths are measured. Allowed branch length units are millions of
        years, years, and generations; these should be specified with the strings
        ``"myr"``, ``"yr"``, or ``"gen"``, respectively. This defaults to
        ``"gen"``.
    :param float generation_time: The number of years per generation. If and only
        if the branch lengths are not in units of generations, the generation time
        must be specified. This defaults to `None`.
    :type generation_time: float or None
    :return: A tuple of two lists of which the first contains
        :class:`.PopulationConfiguration` instances and the second contains
        :class:`.MassMigration` instances. The population configurations specify
        the size of each population and the species name corresponding to each
        population. Species names are stored as metadata in each
        :class:`.PopulationConfiguration` instance, with the metadata tag
        "species_name". Sampling configurations and growth rates are not specified
        in the population configurations. The list of population configurations is
        ordered according to the order of the corresponding extant species in a
        `post-order tree traversal
        <https://en.wikipedia.org/wiki/Tree_traversal#Post-order_(LRN)>`_. The list
        of mass migration events is ordered by the time of the events, from young
        to old events.
    :rtype: (list, list)
    :warning: This function does not modify migration matrices. When the population
        configurations and mass migration events returned by this function are used
        to simulate with the :func:`.simulate` function, it should be ensured that
        migration rates to source populations of mass migration events are zero
        after the mass migration (viewed backwards in time).
    """

    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = "The specified units for branch lengths ("
        err += f'"{branch_length_units}") are not accepted. '
        err += 'Accepted units are "myr" (millions of years), "yr" (years), '
        err += 'and "gen" (generations).'
        raise ValueError(err)

    try:
        Ne = float(Ne)
    except ValueError:
        raise ValueError("Population size Ne must be numeric.")
    if Ne <= 0:
        raise ValueError("Population size Ne must be > 0.")

    # Make sure that the generation time is either None or positive.
    if generation_time is not None:
        generation_time = check_generation_time(generation_time)

    # Make sure that the generation time is specified if and only if
    # branch lengths are not in units of generations.
    if branch_length_units == "gen":
        if generation_time is not None:
            err = 'With branch lengths in units of generations ("gen"), '
            err += "a generation time should not be specified additionally."
            raise ValueError(err)
    else:
        if generation_time is None:
            err = "With branch lengths in units of "
            err += f'"{branch_length_units}", a generation time must be '
            err += "specified additionally."
            raise ValueError(err)

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time
    )

    # Parse the tree with the newick library.
    root = parse_newick(tree, generations_per_branch_length_unit)

    # Define populations and demographic events according to the
    # specified population size and the divergence times in the species tree.
    # Per divergence event (node in the tree), a mass migration with a proportion
    # of 1 of the population is used. The destination is the left-most leaf for
    # each node. Because we are using the n leaf populations, we map each node back
    # to the leaf population that it corresponds to.
    populations = []
    events = []
    leaf_map = {}
    for node in root.walk("postorder"):
        if len(node.descendants) == 0:
            # Per extant species (= leaf node) in the tree, add a population with
            # size Ne. Species names are stored as metadata with the "species_name"
            # tag.
            populations.append(
                demog.Population(initial_size=Ne, name=node.name.strip())
            )
            leaf_map[node] = len(populations) - 1
        else:
            # Per internal node, add one (if the node is bifurcating) or multiple
            # (if the node is multi-furcating) MassMigrations. The parent species
            # maps implicitly to the left-most child species, so we don't generate
            # any MassMigrations for that.
            leaf_map[node] = leaf_map[node.descendants[0]]
            # For each child species after the left-most one, we create
            # a MassMigration into the left-most species.
            for child in node.descendants[1:]:
                events.append(
                    demog.MassMigration(
                        source=leaf_map[child], dest=leaf_map[node], time=node.time
                    )
                )
    return demog.Demography(populations=populations, events=events)


def process_starbeast_tree(
    tree_string, generations_per_branch_length_unit, species_name_map
):
    """
    Process the specified starbeast newick string with embedded dmv annotations
    (but no others) and return the resulting population_configurations and
    demographic_events.
    """
    root = parse_newick(tree_string, generations_per_branch_length_unit)
    populations = []
    events = []
    # The destination is the left-most leaf for each node. Because we are
    # using the n leaf populations, we map each node back to the leaf
    # population that it corresponds to.
    leaf_map = {}
    for node in root.walk("postorder"):
        if node.name is None:
            raise ValueError("Annotation missing for one or more nodes.")
        find_pattern = "\\&dmv=\\{([\\d\\.]+?)\\}"
        dmv_patterns = re.search(find_pattern, node.name)
        if dmv_patterns is None:
            raise ValueError("No dmv annotation for node")
        pop_size = float(dmv_patterns.group(1)) * generations_per_branch_length_unit

        if len(node.descendants) == 0:
            # Per extant species (= leaf node) in the tree, add a population with
            # size pop_size. Species names are stored as metadata with the "species_name"
            # tag.
            newick_id = node.name.strip().split("[")[0]
            species_name = species_name_map[newick_id]
            populations.append(
                demog.Population(initial_size=pop_size, name=species_name)
            )
            leaf_map[node] = len(populations) - 1
        else:
            # Per internal node, add one (if the node is bifurcating) or multiple
            # (if the node is multi-furcating) MassMigrations. The parent species
            # maps implicitly to the left-most child species, so we don't generate
            # any MassMigrations for that.
            leaf_map[node] = leaf_map[node.descendants[0]]
            # For each child species after the left-most one, we create
            # a MassMigration into the left-most species.
            for child in node.descendants[1:]:
                events.append(
                    demog.MassMigration(
                        source=leaf_map[child], dest=leaf_map[node], time=node.time
                    )
                )
            events.append(
                demog.PopulationParametersChange(
                    node.time, initial_size=pop_size, population_id=leaf_map[node]
                )
            )
    return demog.Demography(populations=populations, events=events)


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def check_generation_time(generation_time):
    try:
        generation_time = float(generation_time)
    except ValueError:
        raise ValueError("Generation time must be numeric.")
    if generation_time <= 0:
        raise ValueError("Generation time must be > 0.")
    return generation_time


def get_generations_per_branch_length_unit(branch_length_units, generation_time):
    """
    Method to calculate the number of generations per branch length
    unit, given the branch length unit and a generation time.
    """
    if branch_length_units == "gen":
        generations_per_branch_length_unit = 1
    elif branch_length_units == "myr":
        generations_per_branch_length_unit = 10 ** 6 / generation_time
    else:
        generations_per_branch_length_unit = 1 / generation_time
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


def strip_extra_annotations(tree_string):
    """
    Takes the input newick string and strips all extended newick annotations
    other than the dmv attribute, returning the simplified newick string.
    """
    if "[" not in tree_string:
        raise ValueError("No annotation in tree string.")
    if tree_string.count("[") != tree_string.count("]"):
        raise ValueError("Unbalanced square brackets in annotation.")
    if "&dmv={" not in tree_string:
        raise ValueError("No dmv tag in annotation.")
    if "}" not in tree_string:
        raise ValueError("No closing curly brackets in annotation.")

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
            if tree_string[x - 1] == "[" and tree_string[x] == "&":
                clean_tree_string += "&"
            if tree_string[x - 5 : x] == "dmv={":
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

    return clean_tree_string


def parse_translate_command(translate_command):
    """
    Parses the species IDs used in a nexus newick string to their
    more verbose species names. Returns a dictionary mapping the newick
    values to the species names.
    """
    # Use the translation block to back-translate species IDs in
    # the tree string. The Nexus format definition does not
    # define the format of the translate block. Here, we only
    # allow comma-separated translation pairs with the species
    # name used in the tree string to the left and the translation
    # to the right.
    # An example of an allowed translation is:
    # "translate 1 spc1, 2 spc2, 3 spc3;"

    # First, we trim the "translate" tag from the beginning.
    assert translate_command[0:10] == "translate "
    translate_command = translate_command[10:]
    mapping = {}
    for item in translate_command.split(","):
        item_list = item.split()
        if len(item_list) <= 1:
            raise ValueError("Missing translation in translation block.")
        if len(item_list) != 2:
            err = "Species IDs in the translation block appear to include "
            err += "whitespace. This is not supported."
            raise ValueError(err)
        newick_id, species_name = item_list
        if newick_id in mapping:
            raise ValueError(
                f"Newick ID {newick_id} defined multiple times in translation"
            )
        mapping[newick_id] = species_name
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Duplicate species names in translation")
    return mapping


def parse_nexus(nexus):
    """
    Parse the specified nexus string, returning translate and tree comand
    strings.

    NOTE because we're assuming that the data is generated by starbeast we
    aren't exhaustive in checking for malformed input. We try to catch
    a lot of errors and to give good error messages in these cases. We also
    put a lot of assertions to make sure that we're not silently
    accepting malformed input data. Nonetheless, this is definitely not
    a general purpose Nexus parser and should not be expected to work on
    anything other than input that closely resembles starbeast output.
    """
    # From the Nexus format definition (Maddison et al. 1997):
    # "For the most part, whitespace, inluding newline characters, is ignored,
    # with two exceptions: (1) whitespace indicates boundaries between words;
    # (2) in interleaved matrices, newline characters indicate the boundary
    # between data of different taxa."
    # As we do not parse matrices (we're only interested in taxa and trees
    # blocks), we ignore (2), replace newline characters with spaces and
    # replace multiple whitespaces with a single one.
    nexus_string = nexus.replace("\n", " ")
    nexus_string = " ".join(nexus_string.split())

    # From the Nexus format definition (Maddison et al. 1997):
    # "Commands or subcommands that differ only in case are homonymous."
    # We turn the whole nexus string into lowercase.
    nexus_string = nexus_string.lower()

    # Make sure that the string is in Nexus format.
    if nexus_string[0:6] != "#nexus":
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
        if nexus_string[x : x + 12] == "begin trees;":
            in_tree_block = True
            tree_block_string += nexus_string[x]
        elif in_tree_block and nexus_string[x : x + 4] == "end;":
            tree_block_string += nexus_string[x : x + 4]
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
    tree_string = tree_command[tree_command.find("(") :]
    return translate_command, tree_string
