"""
Examples illustrating the parsing of species trees from files in Newick
or Nexus format, and the definition of tuples combining population
configurations and demographic events.
"""
import msprime
import sys

# Parse the tree from file primates.tre, written in plain Newick format,
# generate a tree sequence based on this species tree, and inspect the
# demography.
print("Simple newick tree with branch lengths in units of millions of years:")
parsed_tuple = msprime.parse_species_tree(
        species_tree="(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)",
        branch_length_units="myr",
        Ne=10000,
        generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()

# Do the same with the species tree from file primates_years.tre, a
# modified version of primates.tre in which branch lengths are in units
# of years instead of millions of years.
print("\n\nSimple newick tree with branch lengths in units of years:")
parsed_tuple = msprime.parse_species_tree(
        species_tree="(((human:5600000,chimpanzee:5600000):3000000,gorilla:8600000):9400000,orangutan:18000000)",
        branch_length_units="yr",
        Ne=10000,
        generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()

# Do the same with the species tree from file primates_simultaneous.tre,
# a modified version of primates.tre in which two divergence events are
# simultaneous.
print("\n\nSimple newick tree with two simultaneous divergence events:")
parsed_tuple = msprime.parse_species_tree(
        species_tree="((human:5,chimpanzee:5):3,(gorilla:5,orangutan:5):3)",
        branch_length_units="myr",
        Ne=10000,
        generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()


# Parse the species tree primates_polytomy.tre with a polytomy.
print("\n\nSimple newick tree with polytomy:")
parsed_tuple = msprime.parse_species_tree(
        species_tree="((human:8.6,chimpanzee:8.6,gorilla:8.6):9.4,orangutan:18.0)",
        branch_length_units="myr",
        Ne=10000,
        generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()

# Parse the non-ultrametric species tree in file.
print("\n\nNon-ultrametric newick tree:")
parsed_tuple = msprime.parse_species_tree(
        species_tree="(((human:5.6,chimpanzee:5.6):3.0,gorilla:7.6):9.4,orangutan:18.0)",
        branch_length_units="myr",
        Ne=10000,
        generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()

# Parse the tree from file 91genes_species_rev.tre, written in StarBEAST
# format and including a translation block and population sizes for each
# branch.
# Generate a tree sequence based on the species tree and draw the first
# tree from the tree sequence.
# As the population sizes are encoded in the tree, Ne does not need to be specified.
print("\n\nStarBEAST species tree with annotation for population sizes:")
with open("tests/data/species_trees/91genes_species_rev.tre", "r") as f:
    parsed_tuple = msprime.parse_starbeast(
            species_tree=f.read(),
            branch_length_units="myr",
            generation_time=5)
population_configurations = parsed_tuple[0]
demographic_events=parsed_tuple[1]
for n in population_configurations:
    n.sample_size = 2
tree_sequence = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        recombination_rate=1e-7,
        length=100)
tree = tree_sequence.first()
print(tree.draw(format="unicode"))
print("number of trees: ", tree_sequence.num_trees)


# Parse the large species tree with over 100 species from
# 101g_nucl_conc_unconst.combined.tre, written in Nexus format but
# without population sizes per branch.
# This tree is from Musilova et al. (2019), available from
# http://evoinformatics.eu/opsin_evolution.htm.
print("\n\nLarge newick tree with over 100 species:")
with open("tests/data/species_trees/101g_nucl_conc_unconst.combined.nwk.tre", "r") as f:
    parsed_tuple = msprime.parse_species_tree(
            species_tree=f.read(),
            branch_length_units="myr",
            Ne=1000,
            generation_time=5)
population_configurations = parsed_tuple[0]
demographic_events=parsed_tuple[1]
for n in population_configurations:
    n.sample_size = 2
tree_sequence = msprime.simulate(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1],
        recombination_rate=1e-8,
        length=400)
print("number of trees: ", tree_sequence.num_trees)
for tree in tree_sequence.trees():
        print("interval: ", tree.interval)
