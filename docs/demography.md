(sec_demography)=

# Demography

:::{warning}

This documentation is under heavy construction. Please note
any outstanding TODOs before opening issues.

:::

```{eval-rst}
.. todo:: This section needs some intro material
```

{class}`.Demography`                       
: Description of a demographic model       

**Constructing simple models**                                                       

{meth}`.Demography.isolated_model`        
: A model of isolated populations           

{meth}`.Demography.island_model`          
: An interconnected web of populations      

{meth}`.Demography.stepping_stone_model`  
: Simple spatial model                      

**Constructing models from existing definitions**                                    

{meth}`.Demography.from_species_tree`     
: Parse a newick species tree               

{meth}`.Demography.from_starbeast`        
: Parse StarBeast output                    

{meth}`.Demography.from_old_style`        
: Demography from msprime legacy 0.x inputs        

**Debugging**                                                                        

{meth}`.Demography.debug`                  
: Get a Debugger for a demography          
 
{class}`.DemographyDebugger`               
: Debugger for demographic models          

---

(sec_demography_migration)=

## Migration rates

(sec_demography_events)=

## Demographic Events

## Species trees

```{eval-rst}
.. todo:: Review this documentation which is a combination of the old API and
    tutorial docs.
```

Species trees hold information about the sequence and the times at which species
diverged from each other. Viewed backwards in time, divergence events are equivalent
to mass migration events in which all lineages from one population move to another
population. The history of a set of populations can thus be modelled according to
a given species tree. To faciliate the specification of the model,
{meth}`.Demography.from_species_tree` parses a species tree and returns
the corresponding Demography object.

When species trees are estimated with a program like [StarBEAST](<https://academic.oup.com/mbe/article/34/8/2101/3738283>) they can further
contain estimates on the population sizes of extant and ancestral species.
{meth}`.Demography.from_starbeast` parses species trees estimated with StarBEAST and uses
these estimates to define the population configurations.

Note that when the species tree has branch lengths not in units of generations but
in units of years or millions of years (which is common), a generation time in years
is required for parsing.

```{eval-rst}
.. todo:: refactor this content to use the new conventions and jupyter.
```

Models used in `msprime` for simulation can be designed to approximate the
diversification history of a group of diverging species, by defining, for each species
divergence, a mass migration event at which all lineages from one population move
into another population. To faciliate the specification of these mass migration events
it is possible to parse a species tree and generate the set of mass migration events
automatically.

To be parseable, a species tree must be encoded in
[Newick](<https://en.wikipedia.org/wiki/Newick_format>) format, with named leaves and
branch lengths. One example of a parseable species tree in Newick format is
`(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)`. When visualized
in a software like [FigTree](<http://tree.bio.ed.ac.uk/software/figtree/>), this tree
appears as shown below:

```{image} _static/primates.svg
:width: 400px
:alt: A species tree with branch lengths.

```

In the above figure, numbers written on branches indicate the lengths of these branches.
In this case, the units of the branch lengths are millions of years, which is common
for species trees; however, trees with branch lengths in units of years or generations
can also be parsed. When the branch lengths are in units of millions of years or in units
of years (i.e., not in units of generations), a generation time in number of years must
be provided so that the simulation model can be set up. In addition, a population size
is required. With the species tree, a generation time, and the population size, the model
can be generated:

```{code-block} python

import msprime
pop_configs, demographic_events = msprime.parse_species_tree(
        tree="(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)",
        Ne=10000,
        branch_length_units="myr",
        generation_time=28)

```

The `parse_species_tree()` method returns a tuple of two lists:

```{code-block} python

print(type(pop_configs))
# <class 'list'>
print(type(demographic_events))
# <class 'list'>

```

The first of these two lists contains instances of {class}`.PopulationConfiguration` and
the second contains instances of {class}`.MassMigration`:

```{code-block} python

print(type(pop_configs[0]))
# <class 'msprime.simulations.PopulationConfiguration'>
print(type(demographic_events[0]))
# <class 'msprime.simulations.MassMigration'>

```

The mass migration events are ordered by the time of the event and they thus specify
the sequence in which lineages are moved from source to destination populations:

```{code-block} python

print(demographic_events[0])
# Mass migration: Lineages moved with probability 1.0 backwards in time with source 1 & dest 0
                     (equivalent to migration from 0 to 1 forwards in time)
print(demographic_events[1])
# Mass migration: Lineages moved with probability 1.0 backwards in time with source 2 & dest 0
                     (equivalent to migration from 0 to 2 forwards in time)
print(demographic_events[2])
# Mass migration: Lineages moved with probability 1.0 backwards in time with source 3 & dest 0
                     (equivalent to migration from 0 to 3 forwards in time)

```

The above output indicates that --- viewed backwards in time --- lineages from
populations 1, 2, and 3 are consecutively moved into population 0. Viewed forwards in
time instead, this means that population 3 is the first to diverge, followed by
population 2 and finally the divergence between populations 0 and 1. This
sequence of divergences corresponds to the species tree if population 3 is orangutan
and populations 2, 1, and 0 are gorilla, chimpanzee, and human, respectively. While
the parsed species names are not used as population labels, they are included in the
population configurations in the form of metadata, with the "species_name" tag:

```{code-block} python

print([pop_config.metadata for pop_config in pop_configs])
# [{'species_name': 'human'}, {'species_name': 'chimpanzee'}, {'species_name': 'gorilla'}, {'species_name': 'orangutan'}]

```

As the above output shows, the information on the topology of the species tree is
fully included in the set of population configurations and mass migration events.
It also illustrates that it is always the left one (when viewed from the root
towards the tips of the species tree in the tree figure above) of the two populations
descending from a divergence event that is used as the destination population in
mass migration events.

The population configurations also define the population size:

```{code-block} python

print([pop_config.initial_size for pop_config in pop_configs])
# [10000.0, 10000.0, 10000.0, 10000.0]

```

The population size is 10,000 because this value was specified for Ne when
calling {func}`.parse_species_tree`. The growth rates are zero for all populations,
meaning that they all have constant population sizes:

```{code-block} python

print([pop_config.growth_rate for pop_config in pop_configs])
# [0.0, 0.0, 0.0, 0.0]

```

To simulate under the model corresponding to the species tree, the population
configurations and mass migration events are used as input for
{func}`.simulate`. We can specify the genomes to sample either by using the
`samples` parameter or by setting a sample size for each population:

```{code-block} python

for pop_config in pop_configs:
    pop_config.sample_size = 2
# pop_configs now has 2 genomes sampled from each population at time 0. Equivalent to
# msprime.simulate(samples=[(0,0), (0,0), (1,0), (1,0), (2,0), (2,0), (3,0), (3,0)], ...)
tree_sequence = msprime.simulate(
        population_configurations=pop_configs,
        demographic_events=demographic_events)
tree = tree_sequence.first()
print(tree.draw(format="unicode"))
#    14
#  ┏━━┻━━━┓
#  ┃     13
#  ┃   ┏━━┻━━┓
#  ┃   ┃    12
#  ┃   ┃   ┏━┻━┓
#  ┃   ┃   ┃  11
#  ┃   ┃   ┃  ┏┻┓
# 10   ┃   ┃  ┃ ┃
# ┏┻┓  ┃   ┃  ┃ ┃
# ┃ ┃  ┃   9  ┃ ┃
# ┃ ┃  ┃  ┏┻┓ ┃ ┃
# ┃ ┃  8  ┃ ┃ ┃ ┃
# ┃ ┃ ┏┻┓ ┃ ┃ ┃ ┃
# 6 7 4 5 2 3 0 1

```

The correspondence between the model and the species tree can also be verified
by using the demography debugger:

```{code-block} python

dd = msprime.DemographyDebugger(
    population_configurations=pop_configs,
    demographic_events=demographic_events)
dd.print_history()
# Model =  hudson(reference_size=1)
# ================================
# Epoch: 0 -- 200000.0 generations
# ================================
#      start     end      growth_rate |     0        1        2        3
#    -------- --------       -------- | -------- -------- -------- --------
# 0 |  1e+04    1e+04               0 |     0        0        0        0
# 1 |  1e+04    1e+04               0 |     0        0        0        0
# 2 |  1e+04    1e+04               0 |     0        0        0        0
# 3 |  1e+04    1e+04               0 |     0        0        0        0
#
# Events @ generation 200000.0
#    - Mass migration: Lineages moved with probability 1.0 backwards in time with source 1 & dest 0
#                      (equivalent to migration from 0 to 1 forwards in time)
# =================================================
# Epoch: 200000.0 -- 307142.85714285716 generations
# =================================================
#      start     end      growth_rate |     0        1        2        3
#    -------- --------       -------- | -------- -------- -------- --------
# 0 |  1e+04    1e+04               0 |     0        0        0        0
# 1 |  1e+04    1e+04               0 |     0        0        0        0
# 2 |  1e+04    1e+04               0 |     0        0        0        0
# 3 |  1e+04    1e+04               0 |     0        0        0        0
#
# Events @ generation 307142.85714285716
#    - Mass migration: Lineages moved with probability 1.0 backwards in time with source 2 & dest 0
#                      (equivalent to migration from 0 to 2 forwards in time)
# =========================================================
# Epoch: 307142.85714285716 -- 642857.142857143 generations
# =========================================================
#      start     end      growth_rate |     0        1        2        3
#    -------- --------       -------- | -------- -------- -------- --------
# 0 |  1e+04    1e+04               0 |     0        0        0        0
# 1 |  1e+04    1e+04               0 |     0        0        0        0
# 2 |  1e+04    1e+04               0 |     0        0        0        0
# 3 |  1e+04    1e+04               0 |     0        0        0        0
#
# Events @ generation 642857.142857143
#    - Mass migration: Lineages moved with probability 1.0 backwards in time with source 3 & dest 0
#                      (equivalent to migration from 0 to 3 forwards in time)
# ==========================================
# Epoch: 642857.142857143 -- inf generations
# ==========================================
#      start     end      growth_rate |     0        1        2        3
#    -------- --------       -------- | -------- -------- -------- --------
# 0 |  1e+04    1e+04               0 |     0        0        0        0
# 1 |  1e+04    1e+04               0 |     0        0        0        0
# 2 |  1e+04    1e+04               0 |     0        0        0        0
# 3 |  1e+04    1e+04               0 |     0        0        0        0

```

The epoch boundaries 200000, 307142.9, and 642857.1 correspond to the species
divergence times 5.6, 8.6, and 18.0 after converting the branch length units
of the species tree from millions of years to generations with the specified
generation time of 28 years.


