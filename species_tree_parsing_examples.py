"""
Examples illustrating the parsing of species trees from files in Newick
or Nexus format, and the definition of tuples combining population
configurations and demographic events.
"""
import msprime

# Parse the tree from file primates.tre, written in plain Newick format, generate a tree sequence based on this species tree, and inspect the demography.
print("\n\nprimates.tre")
parsed_tuple = msprime.parse_species_tree(filename="tests/data/species_trees/primates.tre", branch_length_units="myr", sample_size=3, Ne=10000, generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()


# Do the same with the species tree from file primates_years.tre, a modified version of primates.tre in which branch lengths are in units of years instead of millions of years.
print("\n\nprimates_years.tre")
parsed_tuple = msprime.parse_species_tree(filename="tests/data/species_trees/primates_years.tre", branch_length_units="yr", sample_size=3, Ne=10000, generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()


# Do the same with the species tree from file primates_simultaneous.tre, a modified version of primates.tre in which two divergence events are simultaneous.
print("\n\nprimates_simultaneous.tre")
parsed_tuple = msprime.parse_species_tree(filename="tests/data/species_trees/primates_simultaneous.tre", branch_length_units="myr", sample_size=3, Ne=10000, generation_time=28)
dd = msprime.DemographyDebugger(
        population_configurations=parsed_tuple[0],
        demographic_events=parsed_tuple[1])
dd.print_history()


# Try parsing the species tree primates_polytomy.tre with a polytomy (this should result in an error).
print("\n\nprimates_polytomy.tre")
try:
	parsed_tuple = msprime.parse_species_tree(filename="tests/data/species_trees/primates_polytomy.tre", branch_length_units="myr", sample_size=3, Ne=10000, generation_time=28)
except:
	print("The species tree primates_polytomy.tre could not be parsed.")


# Try parsing the non-ultrametric species tree in file  (this should result in an error).
print("\n\nprimates_nonultrametric.tre")
try:
	parsed_tuple = msprime.parse_species_tree(filename="tests/data/species_trees/primates_nonultrametric.tre", branch_length_units="myr", sample_size=3, Ne=10000, generation_time=28)
except:
	print("The species tree primates_nonultrametric.tre could not be parsed.")


# Parse the tree from file 91genes_species_rev.tre, written in StarBEAST format and including a translation block and population sizes for each branch.
# Generate a tree sequence based on the species tree and draw the first tree from the tree sequence.
# As the population sizes are encoded in the tree, Ne does not need to be specified.
print("\n\n91genes_species_rev.tre")
parsed_tuple = msprime.parse_species_tree(filename="tests/data/species_trees/91genes_species_rev.tre", branch_length_units="myr", sample_size=2, generation_time=5)
tree_sequence = msprime.simulate(population_configurations=parsed_tuple[0], demographic_events=parsed_tuple[1], recombination_rate=1e-7, length=100)
tree = tree_sequence.first()
print(tree.draw(format="unicode"))
print("number of trees: ", tree_sequence.num_trees)


# Parse the large species tree with over 100 species from 101g_nucl_conc_unconst.combined.tre, written in Nexus format but without population sizes per branch.
# This tree is from Musilova et al. (2019), available from http://evoinformatics.eu/opsin_evolution.htm.
print("\n\n101g_nucl_conc_unconst.combined.tre")
parsed_tuple = msprime.parse_species_tree(filename="tests/data/species_trees/101g_nucl_conc_unconst.combined.tre", branch_length_units="myr", sample_size=2, Ne=1000, generation_time=5)
tree_sequence = msprime.simulate(population_configurations=parsed_tuple[0], demographic_events=parsed_tuple[1], recombination_rate=1e-8, length=400)
print("number of trees: ", tree_sequence.num_trees)
for tree in tree_sequence.trees():
    print("interval: ", tree.interval)
