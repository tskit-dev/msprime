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


```{code-cell}
:tags: [remove-cell]

    import msprime
    from IPython.display import SVG

```

(sec_demography)=

# Demography

By default, msprime assumes a single randomly mating population of a fixed
size for {ref}`ancestry simulations <sec_ancestry>`, which is unrealistic
for most purposes. To enable more realistic and complex simulations,
msprime models population structure
by defining a set of discrete {ref}`sec_demography_populations`,
with {ref}`sec_demography_migration` between these populations occuring
at some rate. The parameters of this model can change over time, via
{ref}`sec_demography_events`.  Please see the {ref}`sec_demography_definitions`
section for mathematical details.

The information required to define a demographic model is encapsulated
by the {class}`.Demography` object. To run
{ref}`ancestry simulations <sec_ancestry>` for a given demography
we use the `demography` parameter to {func}`.sim_ancestry` and define
what populations our {ref}`samples <sec_ancestry_samples>` are drawn from.
For example, here we create a three-population island model (each
with population size 100), and take one diploid sample from each
of the three populations:

```{code-cell}
demography = msprime.Demography.island_model([100, 100, 100], 0.1)
ts = msprime.sim_ancestry(
    samples={pop_id: 1 for pop_id in range(demography.num_populations)},
    demography=demography,
    random_seed=1234
)
ts
```

---

## Quick reference

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

**Numerical methods**

{meth}`.DemographyDebugger.coalescence_rate_trajectory`
: Compute mean coalescence rate and fraction of uncoalesced lineages

{meth}`.DemographyDebugger.mean_coalescence_time`
: Compute the mean time until coalescence

{meth}`.DemographyDebugger.lineage_probabilities`
: Probability of the location of lineages over time

---

(sec_demography_definitions)=

## Definitions

Formally, population structure in msprime is modelled by specifying a fixed
number of subpopulations {math}`d`.
Each population has an initial absolute population size {math}`s`
and a per generation exponential growth rate {math}`\alpha`. The size of a
given population at time {math}`t` in the past (measured in generations) is
therefore given by {math}`s e^{-\alpha t}`. Demographic events that occur in
the history of the simulated population alter some aspect of this population
configuration at a particular time in the past.

Continuous migration between populations is modelled by a
{math}`d \times d` matrix {math}`M` of per-generation migration rates.
The {math}`(j,k)^{th}` entry
of {math}`M` is the expected number of migrants moving from population
{math}`k` to population {math}`j` per generation, divided by the size of
population {math}`j`. In terms of the coalescent process, {math}`M_{j,k}`
gives the rate at which an ancestral lineage moves from population
{math}`j` to population {math}`k`, as one follows it back through time. In
continuous-time models, when {math}`M_{j,k}` is close to zero, this rate is
approximately equivalent to the fraction of population {math}`j` that is
replaced each generation by migrants from population {math}`k`. In
discrete-time models, the equivalence is exact and each row of {math}`M`
has the constraint {math}`\sum_{k \neq j} M_{j,k} \leq 1`. This differs
from the migration matrix one usually uses in population demography: if
{math}`m_{k,j}` is the proportion of individuals (in the usual sense; not
lineages) in population {math}`k` that move to population {math}`j` per
generation, then translating this proportion of population {math}`k` to a
proportion of population {math}`j`, we have
{math}`M_{j,k} = m_{k,j} \times N_k / N_j`.

The details of population structure in msprime closely follow the model
used in the classical
[ms](https://academic.oup.com/bioinformatics/article/18/2/337/225783)
program.

---

(sec_demography_populations)=

## Populations

A {class}`.Demography` contains a list of {class}`.Population` objects. Populations
essentially have two purposes:

1. To define the state of the population at the start of the simulation (i.e.,
    the {ref}`sec_demography_populations_initial_size` and
    the {ref}`sec_demography_populations_growth_rate`)
2. To define {ref}`sec_demography_populations_metadata` which is
    associated with the corresponding {class}`tskit.Population` objects
    in simulated tree sequences.

When used as part of the {class}`.Demography`, each population has an integer
ID (its zero-based index in the `populations` list) and a ``name`` attribute.
By default, msprime assigns the name `pop_j` to the population at index
`j` in a demographic model. These default names can be overridden, and
users should give populations more descriptive names when building complex models

:::{note}
Population names must be unique within a single {class}`.Demography` and be valid
[Python identifiers](https://docs.python.org/3/reference/lexical_analysis.html#identifiers).
For example, this means that names like "my_pop_1" is valid, but "my-pop-1" and "my pop 1"
are not.
:::

(sec_demography_populations_initial_size)=

### Initial size

The `initial_size` of a population is its size at the start of a simulation
(looking backwards in time). If the population's
{ref}`sec_demography_populations_growth_rate` is zero, then
the population will have the same size for all time (unless there are
some {ref}`sec_demography_events` that change these parameters).

When using methods like {meth}`.Demography.island_model` to describe
multipopulation models, the ``initial_size`` parameter takes
a **list** of size values which defines both the number of
populations and their sizes. For example, here we define two isolated
populations (i.e., which have no migration) with sizes 100 and 200:

```{code-cell}
demography = msprime.Demography.isolated_model([100, 200])
demography
```

When we have multiple populations that are the same size we can use
some Python tricks to avoid code duplication. For example, here
we create three populations with size 100:

```{code-cell}
demography = msprime.Demography.isolated_model([100] * 3)
demography
```

```{warning}
We use the {meth}`.Demography.isolated_model` function here as a
convenient means of creating example demographies. However, it's
important to note that simulations in which we sample from
multiple isolated populations will fail with an error because
the lineages in question can never coalesce.
```

(sec_demography_populations_growth_rate)=

### Growth rate

Each population has an expoential growth rate so that the
size of a population with initial size {math}`s`
and growth rate {math}`\alpha` is {math}`s e^{-\alpha t}`
at time {math}`t` generations in the past (see the
{ref}`sec_demography_definitions` section for more details).

Growth rates for functions like {meth}`.Demography.island_model`
that construct a {class}`.Demography` are specified in a similar
way to {ref}`sec_demography_populations_initial_size`: we provide
a list of size equal to the number of populations. For
example, here we define 2 populations with different
population sizes and growth rates:

```{code-cell}
demography = msprime.Demography.isolated_model([100, 200], growth_rate=[0.01, 0.02])
demography
```

```{note}
The {class}`DemographyDebugger` is a useful tool to help understand
how population sizes change over time. See

```

(sec_demography_populations_metadata)=

### Metadata

In [tskit](https://tskit.dev/tskit) the {class}`Population<tskit.Population>` class
largely exists as a container for {ref}`tskit:sec_metadata`. Metadata is
important and useful: it let's us associated information about our simulated
populations with the output tree sequence; we can then use this information
later when we are analysing the data.

Msprime associates two pieces of mandatory metadata with every population:
their `name` and `description`. For example,

```{code-cell}
demography = msprime.Demography.stepping_stone_model([100, 100], migration_rate=0.01)
ts = msprime.sim_ancestry({0: 1, 1: 1}, demography=demography)
print([population.metadata for population in ts.populations()])
```

Here we have two populations in the output tree sequence, and the metadata
for each population is a dictionary containing the keys ``name`` and
``description``. These correspond to the same attributes on the msprime
{class}`.Population` class. We don't have to just use the defaults
for these values: we can set them to (more or less) anything we like.
For example,

```{code-cell}
demography = msprime.Demography.stepping_stone_model([100, 100], migration_rate=0.01)
demography.populations[0].name = "awesome_pop"
demography.populations[0].description = "This population is totally awesome"
demography.populations[1].name = "more_awesome_pop"
demography.populations[1].description = "This population is even more awesome"
ts = msprime.sim_ancestry({0: 1, 1: 1}, demography=demography)
for population in ts.populations():
    print(f"id={population.id}: {population.metadata}")
```

As well as the default ``name`` and ``description`` metadata keys
we can also associate additional metadata with population objects using
the ``extra_metadata`` attribute of the msprime {class}`.Population`
object. For example,

```{code-cell}
demography = msprime.Demography.stepping_stone_model([100, 100], migration_rate=0.01)
demography.populations[0].name = "awesome_pop"
demography.populations[0].extra_metadata = {"emoji": "👍"}
demography.populations[1].name = "more_awesome_pop"
demography.populations[1].extra_metadata = {"emoji": "🤘"}
ts = msprime.sim_ancestry({0: 1, 1: 1}, demography=demography)
for pop in ts.populations():
    print(pop.id, "\t", pop.metadata["emoji"], pop.metadata["name"])
```

(sec_demography_migration)=

## Migration

Migration is the process of lineages moving from one population to another
during the course of the simulation.
This either happens through continuous migration, where a rate of
migration between each pair of populations is defined, or through
{ref}`sec_demography_events_mass_migration` events. In this section
we are concerned with continuous migration.

As described in the {ref}`sec_demography_definitions` section,
continuous migration between populations is modelled by a matrix
of rates, so that `M[j, k]` is the rate at which lineages move from
population `j` to population `k` in the **coalescent process**,
that is, **backwards in time**. Lineages that move from population
`j` to `k` backwards in time actually correspond to individuals
migrating from population `k` to `j` **forwards in time**.

:::{note}
If you're confused at this point, don't worry. Everyone is confused by
this.
:::

Let's take an example to clarify. Suppose we have a two population model
in which we have migration at rate 0.1 from population `0` to `1` and no
migration in the other direction. We'll then take one haploid sample
(for simplicity) from each population:

```{code-cell}
demography = msprime.Demography.isolated_model([100, 100])
# Rate of backwards migration from population 0 to 1
demography.migration_matrix[0, 1] = 0.1
ts = msprime.sim_ancestry(
    samples={0: 1, 1: 1},
    demography=demography,
    ploidy=1,
    random_seed=12345)
ts.tables.nodes
```

Because we have only two samples and no recombination we have
only one coalescent event in the simulation, and therefore
only one ancestral node. We are interested in what
populations nodes are associated with, which is shown in
the ``population`` column of the node table. We can see
that our samples, nodes 0 and 1, are associated with populations 0 and 1,
respectively. Node 2, which is the ancestor of nodes
0 and 1, is from population 1 which means that the **lineage**
for node 0 moved from population 0 to population 1 as we traced
its history back through time, so that it could ultimately
coalesce with lineage 1 in population 1. However, **forwards** in
time, this must mean that one of the individuals along node 0's lineage
(each lineage is a succession of individuals, passing on genetic
information through the generations) must have migrated
from population ``1`` to population ``0``.

:::{note}
If you're still confused, don't worry, it's still OK. Just remember
that migration rates are confusing and come back to check the
documentation whenever you need to work with them.
:::

```{eval-rst}
.. todo:: Put in plug for demes here.
```

(sec_demography_events)=

## Demographic Events

Setting the population parameters and migration matrix in the {class}`.Demography`
object define the state of demographic model at the start of the simulation
(i.e., the present). We are often interested in population models in which
these parameters change over time; these are implemented through a set
of demographic events, which make some instantaneous changes to the state
of the simulation.

The {class}`.DemographyDebugger` is very useful for inspecting the
state of the simulation at each of the epochs defined by these events.

### Population parameters change

The {class}`.PopulationParametersChange` event is used to change the
``initial_size`` or ``growth_rate`` of a
{ref}`population <sec_demography_populations>` (or populations)
at a given time. For example,
here we create a two-population island model and add two events that
change the parameters of these populations over time. The first
event changes the size of population 0 to 200, 10 generations
in the past. The second event changes the size of *both* populations
to 10 after 20 generations using the ``population=None`` shorthand.

```{code-cell}
demography = msprime.Demography.island_model([100, 100], migration_rate=0.01)
demography.events = [
    msprime.PopulationParametersChange(time=10, population=0, initial_size=200),
    msprime.PopulationParametersChange(time=20, population=None, initial_size=10),
]
demography.debug()
```

### Migration rate change

{class}`.MigrationRateChange` events are used to change the state of the
{ref}`migration matrix <sec_demography_migration>` at a given time.
Here we create a two population model with no migration. We then
add migration (**backwards in time**; see the
{ref}`sec_demography_migration` section) from population 0 to population 1
at rate 0.1 at generation 10. Then, at generation 20 we set the
migration rate between all pairs of populations to 0.2.

```{code-cell}
demography = msprime.Demography.isolated_model([100, 100])
demography.events = [
    msprime.MigrationRateChange(time=10, source=0, dest=1, rate=0.1),
    msprime.MigrationRateChange(time=20, rate=0.2)
]
demography.debug()
```


(sec_demography_events_mass_migration)=

### Mass Migration

{class}`.MassMigration` events move a proportion of the lineages currently
present in one population (the ``source``) to another (the ``dest``).

:::{warning}
There are two things it's vitally important to realise about mass migrations:

1. Like {ref}`continuous migration <sec_demography_migration>` the source
   and destination populations are from the perspective of the coalescent
   process, that is **backwards in time**. So, if a lineage is moved from
   `source` population ``A`` to `dest` population ``B`` by msprime, this
   corresponds to an individual moving **from** ``B`` **to** ``A``,
   forwards in time. (Yes, this is
   {ref}`confusing <sec_demography_migration>`!)

2. Mass migration events **do not** alter the migration rates between
   populations. So, even if all lineages are moved out of a
   particular population (i.e., if ``proportion=1``) that population
   may still have lineages migrating into it after the event. This
   can easily lead to errors when modelling population splits using
   {class}`.MassMigration` events. This error was present in the
   documentation for an older version of msprime, which lead to an
   incorrect version of an example demographic model being used
   in a number of publications. See
   [this paper](http://dx.doi.org/10.1016/j.ajhg.2020.08.017) for
   more details.
:::

The effect of mass migration events are summarised in the
{class}`.DemographyDebugger` output:

```{code-cell}
demography = msprime.Demography.island_model([100, 100], migration_rate=0.1)
demography.events = [
    msprime.MassMigration(time=10, source=0, dest=1, proportion=0.5),
]
demography.debug()
```


```{eval-rst}
.. todo:: Finish this section when we have the PopulationSplit event
    implemented.
```

## Defining empirical models


```{code-cell}
import math

# First we set out the maximum likelihood values of the various parameters
# given in Table 1.
N_A = 7300
N_B = 2100
N_AF = 12300
N_EU0 = 1000
N_AS0 = 510
# Times are provided in years, so we convert into generations.
generation_time = 25
T_AF = 220e3 / generation_time
T_B = 140e3 / generation_time
T_EU_AS = 21.2e3 / generation_time
# We need to work out the starting population sizes based on the growth
# rates provided for these two populations
r_EU = 0.004
r_AS = 0.0055
N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
# Migration rates during the various epochs.
m_AF_B = 25e-5
m_AF_EU = 3e-5
m_AF_AS = 1.9e-5
m_EU_AS = 9.6e-5
# Population IDs correspond to their indexes in the popupulation
# configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
# initially.

populations = [
    msprime.Population(name="YRI", initial_size=N_AF),
    msprime.Population(name="CEU", initial_size=N_EU, growth_rate=r_EU),
    msprime.Population(name="CHB", initial_size=N_AS, growth_rate=r_AS),
]
migration_matrix = [
    [0, m_AF_EU, m_AF_AS],
    [m_AF_EU, 0, m_EU_AS],
    [m_AF_AS, m_EU_AS, 0],
]
events = [
    # CEU and CHB merge into B with rate changes at T_EU_AS
    msprime.MassMigration(time=T_EU_AS, source=2, destination=1, proportion=1.0),
    msprime.MigrationRateChange(time=T_EU_AS, rate=0),
    msprime.MigrationRateChange(time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
    msprime.MigrationRateChange(time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
    msprime.PopulationParametersChange(
        time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1
    ),
    # Population B merges into YRI at T_B
    msprime.MassMigration(time=T_B, source=1, destination=0, proportion=1.0),
    msprime.MigrationRateChange(time=T_B, rate=0),
    # Size changes to N_A at T_AF
    msprime.PopulationParametersChange(
        time=T_AF, initial_size=N_A, population_id=0
    ),
]
demography = msprime.Demography(
    populations=populations,
    migration_matrix=migration_matrix,
    events=events)
demography.debug()
```


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

## Numerical predictions

```{eval-rst}
.. todo:: Add a section explaining how to use
    :meth:`.DemographyDebugger.lineage_probabilities` and friends.
