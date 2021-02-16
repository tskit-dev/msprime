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
import collections
import re

import newick

from . import demography as demog


def parse_starbeast(tree, generation_time, time_units="myr"):
    """
    Parse a nexusencoded species tree into a Demography object. See the
    documentation of :class:`.Demography.from_starbeast` (the public interface)
    for details.
    """

    # Make sure that branch length units are either "myr" or "yr".
    allowed_branch_lenth_units = ["myr", "yr"]
    if time_units not in allowed_branch_lenth_units:
        err = "The specified units for branch lengths ("
        err += f'"{time_units}") are not accepted. '
        err += 'Accepted units are "myr" (millions of years) or "yr" (years).'
        raise ValueError(err)

    generation_time = check_generation_time(generation_time)
    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        time_units, generation_time
    )

    translate_string, tree_string = parse_nexus(tree)
    species_name_map = parse_translate_command(translate_string)
    clean_tree_string = strip_extra_annotations(tree_string)
    return process_starbeast_tree(
        clean_tree_string, generations_per_branch_length_unit, species_name_map
    )


def parse_number_or_mapping(value, message):
    """
    Interpret the specified value as either a single floating point value,
    or a mapping and returns a mapping.
    """
    try:
        x = float(value)
        value = collections.defaultdict(lambda: x)
    except TypeError:
        if not isinstance(value, collections.abc.Mapping):
            raise TypeError(message)
    return value


def parse_initial_size(initial_size):
    error_message = (
        "initial_size argument must be a single number or a mapping from "
        "species names to their population sizes."
    )
    return parse_number_or_mapping(initial_size, error_message)


def parse_growth_rate(growth_rate):
    error_message = (
        "growth_rate argument must be a single number or a mapping from "
        "species names to their exponential growth rates."
    )
    return parse_number_or_mapping(growth_rate, error_message)


def parse_species_tree(
    tree,
    initial_size,
    *,
    time_units="gen",
    generation_time=None,
    growth_rate=None,
):
    """
    Parse a newick encoded species tree into a Demography object. See the
    documentation of :class:`.Demography.from_species_tree` (the public interface)
    for details.
    """
    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if time_units not in allowed_branch_lenth_units:
        err = "The specified units for branch lengths ("
        err += f'"{time_units}") are not accepted. '
        err += 'Accepted units are "myr" (millions of years), "yr" (years), '
        err += 'and "gen" (generations).'
        raise ValueError(err)

    initial_size = parse_initial_size(initial_size)
    if growth_rate is None:
        growth_rate = 0
    growth_rate = parse_growth_rate(growth_rate)

    # Make sure that the generation time is either None or positive.
    if generation_time is not None:
        generation_time = check_generation_time(generation_time)

    # Make sure that the generation time is specified if and only if
    # branch lengths are not in units of generations.
    if time_units == "gen":
        if generation_time is not None:
            err = 'With branch lengths in units of generations ("gen"), '
            err += "a generation time should not be specified additionally."
            raise ValueError(err)
    else:
        if generation_time is None:
            err = "With branch lengths in units of "
            err += f'"{time_units}", a generation time must be '
            err += "specified additionally."
            raise ValueError(err)

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        time_units, generation_time
    )

    # Parse the tree with the newick library.
    root = parse_newick(tree, generations_per_branch_length_unit)

    # Define populations and demographic events according to the
    # specified population size and the divergence times in the species tree.
    # Each divergence event (node in the tree) corresponds to an ancestral
    # population, and mass migration events with proportion 1 move all lineages
    # from the child populations into this new populationl

    population_id_map = {}
    demography = demog.Demography()

    def add_population(node):
        name = None
        if node.name is not None:
            stripped = node.name.strip()
            if len(stripped) > 0:
                name = stripped
        population = demography.add_population(
            initial_size=initial_size[name],
            growth_rate=growth_rate[name],
            name=name,
        )
        population_id_map[node] = population.name
        return population.name

    # Add in the leaf populations first so that they get IDs 0..n - 1
    for node in root.walk():
        if len(node.descendants) == 0:
            add_population(node)

    # Now add in the internal node populations and the mass migration events
    # joining them.
    for node in root.walk("postorder"):
        if len(node.descendants) > 0:
            population_id = add_population(node)
            child_pops = [population_id_map[child] for child in node.descendants]
            demography.add_population_split(
                time=node.time, derived=child_pops, ancestral=population_id
            )
    demography.sort_events()
    demography.validate()
    return demography


def process_starbeast_tree(
    tree_string, generations_per_branch_length_unit, species_name_map
):
    """
    Process the specified starbeast newick string with embedded dmv annotations
    (but no others) and return the resulting population_configurations and
    demographic_events.
    """
    root = parse_newick(tree_string, generations_per_branch_length_unit)
    demography = demog.Demography()
    population_size_map = {}
    population_id_map = {}

    # The process here follows the same basic logic as parse_species_tree above
    # but with some extra elaborations to account for changing population sizes
    # and details of the extended newick annotations.

    def add_population(node):
        name = None
        newick_id = node.name.strip().split("[")[0]
        if len(newick_id) > 0:
            name = species_name_map[newick_id]

        population = demography.add_population(
            initial_size=population_size_map[node], name=name
        )

        population_id_map[node] = population.name
        return population.name

    for node in root.walk():
        if node.name is None:
            raise ValueError("Annotation missing for one or more nodes.")
        find_pattern = "\\&dmv=\\{([\\d\\.]+?)\\}"
        dmv_patterns = re.search(find_pattern, node.name)
        if dmv_patterns is None:
            raise ValueError("No dmv annotation for node")
        pop_size = float(dmv_patterns.group(1)) * generations_per_branch_length_unit
        population_size_map[node] = pop_size
        if len(node.descendants) == 0:
            add_population(node)

    for node in root.walk("postorder"):
        if len(node.descendants) > 0:
            population_id = add_population(node)
            demography.add_population_split(
                time=node.time,
                ancestral=population_id,
                derived=[population_id_map[child] for child in node.descendants],
            )
    demography.sort_events()
    demography.validate()
    return demography


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


def get_generations_per_branch_length_unit(time_units, generation_time):
    """
    Method to calculate the number of generations per branch length
    unit, given the branch length unit and a generation time.
    """
    if time_units == "gen":
        generations_per_branch_length_unit = 1
    elif time_units == "myr":
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
        # We don't allow non ultrametric trees for now because it's unclear
        # how we should deal with taking samples in this case. The code
        # all works perfectly well other than this, though.
        if node.is_leaf:
            if abs(node.time) > 1e-8:  # Arbitrary cutoff
                raise ValueError(
                    f"All leaf populations must be at time 0: time={node.time}"
                )
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
