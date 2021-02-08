---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

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

(sec_demography_populations)=

## Populations

(sec_demography_populations_initial_size)=

### Initial size



(sec_demography_migration)=

## Migration rates


```{code-cell}

import msprime

demography = msprime.Demography.stepping_stone_model([10] * 3, migration_rate=0.1)
demography
```

(sec_demography_events)=

## Demographic Events

## Importing model definitions

### Species trees

```{eval-rst}
.. todo:: This section needs another pass once the rest of the sections
    have a draft. We also need to improve the visualisation tools here
    to help explain stuff.
```

Species trees hold information about the sequence and the times at which species
diverged from each other. Viewed backwards in time, divergence events are equivalent
to mass migration events in which all lineages from one population move to another
population. The {meth}`.Demography.from_species_tree` method parses a
species tree in [Newick format](https://en.wikipedia.org/wiki/Newick_format)
and returns a {class}`.Demography` object. These species trees do not contain
any information on the sizes of the relevant populations, however, and so these
must be specified separately using the ``initial_size`` argument
(see `ref`{sec_demography_populations_initial_size}).
When species trees are estimated with a program like
[StarBEAST](<https://academic.oup.com/mbe/article/34/8/2101/3738283>) they can
contain estimates on the population sizes of extant and ancestral species.
The {meth}`.Demography.from_starbeast` method parses species trees estimated
with StarBEAST and sets the population sizes accordingly.

:::{note}
When a species tree has branch lengths in units of years or millions of years
rather than generations (which is common), a generation time in years
must be specified.
:::

Species trees must be encoded in
[Newick](<https://en.wikipedia.org/wiki/Newick_format>) format, with named leaves and
branch lengths. Consider the following species tree of four primates, for example:

```
(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)
```

When visualized in software like
[FigTree](http://tree.bio.ed.ac.uk/software/figtree/), the tree
looks like this:

```{figure} _static/primates.svg
:width: 400px
:alt: A species tree with branch lengths.
```

The numbers written on branches indicate the lengths of these branches,
which in this case is expressed in millions of years. We can then
convert this species tree into a {class}`.Demography` using the
{meth}`.Demography.from_species_tree`:

```{code-cell}

import msprime

demography = msprime.Demography.from_species_tree(
    "(((human:5.6,chimpanzee:5.6):3.0,gorilla:8.6):9.4,orangutan:18.0)",
    initial_size=10_000,
    time_units="myr",
    generation_time=28)
demography.debug()
```

Because the species tree does not contain complete information about the
demographic model, we must provide some extra information. The
``initial_size`` parameter lets us specify the size of each of the
populations; here we give all the populations a fixed size of 10,000
(there is much more flexibility here, however). Because the branch
lengths in the species tree are given in millions of years, we also
need provide a ``time_units`` parameter and a ``generation_time``.

The epoch boundaries 200000, 307142.9, and 642857.1 correspond to the species
divergence times 5.6, 8.6, and 18.0 after converting the branch length units
of the species tree from millions of years to generations with the specified
generation time of 28 years.


Running the simulation is then straightforward:

```{code-cell}

ts = msprime.sim_ancestry(
    {"human": 2, "orangutan": 2}, demography=demography, random_seed=2)
ts
```

Note that the names of the populations are recorded in the population
metadata:

```{code-cell}
for population in ts.populations():
    print(population.metadata)
```
