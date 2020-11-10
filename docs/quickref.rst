.. _sec_quickref:

===================
API quick reference
===================

+---------------------------------------+-------------------------------------------+
| **Main simulation functions**                                                     |
+---------------------------------------+-------------------------------------------+
| :func:`.sim_ancestry`                 |  Simulate ancestral topology              |
+---------------------------------------+-------------------------------------------+
| :func:`.sim_mutations`                |  Simulate mutations on a given topology   |
+---------------------------------------+-------------------------------------------+
| **Demography**                                                                    |
+---------------------------------------+-------------------------------------------+
| :class:`.Demography`                  |  Description of demographic model         |
+---------------------------------------+-------------------------------------------+
| :class:`.DemographyDebugger`          |  Debugger for demographic models          |
+---------------------------------------+-------------------------------------------+
| :meth:`.Demography.island_model`      |  Island model demography                  |
+---------------------------------------+-------------------------------------------+
| :meth:`.Demography.stepping_stone_1d` |  1D stepping stone model demography       |
+---------------------------------------+-------------------------------------------+
| :meth:`.Demography.parse_species_tree`|  Demography from newick species tree      |
+---------------------------------------+-------------------------------------------+
| :meth:`.Demography.parse_starbeast`   |  Demography from StarBeast output         |
+---------------------------------------+-------------------------------------------+
| **Ancestry models**                                                               |
+---------------------------------------+-------------------------------------------+
| :class:`.StandardCoalescent`          | Coalescent with recombination ("hudson")  |
+---------------------------------------+-------------------------------------------+
| :class:`.SmcApproxCoalescent`         | Sequentially Markov Coalescent ("smc")    |
+---------------------------------------+-------------------------------------------+
| :class:`.SmcPrimeApproxCoalescent`    | SMC'("smc_prime")                         |
+---------------------------------------+-------------------------------------------+
| :class:`.DiscreteTimeWrightFisher`    | Generation-by-generation Wright-Fisher    |
+---------------------------------------+-------------------------------------------+
| :class:`.BetaCoalescent`              | Beta coalescent multiple-merger           |
+---------------------------------------+-------------------------------------------+
| :class:`.DiracCoalescent`             | Dirac coalescent multiple-merger          |
+---------------------------------------+-------------------------------------------+
| :class:`.SweepGenicSelection`         | Selective sweep at a linked locus         |
+---------------------------------------+-------------------------------------------+
| **Mutation models**                                                               |
+---------------------------------------+-------------------------------------------+
| :class:`.BinaryMutationModel`         | 0/1 flip-flopping alleles                 |
+---------------------------------------+-------------------------------------------+
| :class:`.JC69MutationModel`           | Jukes & Cantor '69, nucleotides           |
+---------------------------------------+-------------------------------------------+
| :class:`.HKYMutationModel`            | Hasegawa, Kishino & Yano '85, nucleotides |
+---------------------------------------+-------------------------------------------+
| :class:`.F84MutationModel`            | Felsenstein '84, nucleotides              |
+---------------------------------------+-------------------------------------------+
| :class:`.GTRMutationModel`            | general time-reversible, nucleotides      |
+---------------------------------------+-------------------------------------------+
| :class:`.BLOSUM62MutationModel`       | amino acids                               |
+---------------------------------------+-------------------------------------------+
| :class:`.PAMMutationModel`            | amino acids                               |
+---------------------------------------+-------------------------------------------+
| :class:`.MatrixMutationModel`         | general finite-state mutations            |
+---------------------------------------+-------------------------------------------+
| :class:`.InfiniteAllelesMutationModel`| a generic infinite-alleles model          |
+---------------------------------------+-------------------------------------------+
| :class:`.SLiMMutationModel`           | SLiM compatible mutations                 |
+---------------------------------------+-------------------------------------------+
| **Computing likelihoods**                                                         |
+---------------------------------------+-------------------------------------------+
| :func:`.log_arg_likelihood`           |  Likelihood of an ARG topology            |
+---------------------------------------+-------------------------------------------+
| :func:`.log_mutation_likelihood`      |  Likelihood of a set of mutations         |
+---------------------------------------+-------------------------------------------+
| **Deprecated (pre 1.0) simulation functions**                                     |
+---------------------------------------+-------------------------------------------+
| :func:`.simulate`                     |  Simulate ancestral topology & mutations  |
+---------------------------------------+-------------------------------------------+
| :func:`.mutate`                       |  Simulate mutations on a given topology   |
+---------------------------------------+-------------------------------------------+


