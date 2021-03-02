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

:::{warning}
This page is under heavy development and some parts are inconsistent.
See the {ref}`sec_demography_defining_empirical_models` section for
idiomatic examples of using the API to define models from the literature.
:::

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
demography.add_population_parameters_change(
    time=10, population=0, initial_size=200)
demography.add_population_parameters_change(
    time=20, population=None, initial_size=10)
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
demography.add_migration_rate_change(time=10, source=0, dest=1, rate=0.1)
demography.add_migration_rate_change(time=20, rate=0.2)
demography.debug()
```

(sec_demography_events_population_split)=

### Population split

```{eval-rst}
.. todo:: Write this section
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
demography.add_mass_migration(time=10, source=0, dest=1, proportion=0.5)
demography.debug()
```


```{eval-rst}
.. todo:: Finish this section when we have the PopulationSplit event
    implemented.
```

(sec_demography_defining_empirical_models)=

## Defining empirical models

### Splitting populations

```{eval-rst}
.. todo:: Update this text a bit once we've finalised the API.
```

To illustrate ``msprime``'s demography API on a real example from the
literature, we implement the
[Gutenkunst et al.](http://dx.doi.org/10.1371/journal.pgen.1000695)
out-of-Africa model.
The parameter values used are taken from
[Table 1](http://dx.doi.org/10.1371/journal.pgen.1000695.t001).
Here is an illustration of the model using the [demography
package](https://github.com/apragsdale/demography>)
(see also [Figure 2B](http://dx.doi.org/10.1371/journal.pgen.1000695.g002)
of the Gutenkunst et. al paper):

```{image} _static/Gutenkunst_OOA_diagram.svg
:width: 500px
:align: center
:alt: Schematic of Gutenkunst et al. (2009) out-of-Africa model.
```

The code below is provided as an example to help you develop your
own models. If you want to use this precise model in your analyses
we strongly recommend using {ref}`stdpopsim <stdpopsim:sec_introduction>`,
which provides a community maintained {ref}`catalog <stdpopsim:sec_catalog>`
of species information and demographic models for simulation. The
model given here is identical to the
{ref}`HomSam/OutOfAfrica_3G09 <stdpopsim:sec_catalog_homsap_models_outofafrica_3g09>`
model.

```{eval-rst}
.. todo:: The model isn't actually identical any more because we're adding
    the new populations for OOA and ancestral. We should explain this
    distinction, but hopefully the stdpopsim model will be updated to
    use the demes-like approach soon too.
```

:::{warning}
The version of this model in this documentatino from 31 May 2016 to 29 May 2020
(on the stable branch) was **incorrect**. Specifically, it mistakenly
allowed for migration to continue beyond the merger of the African and
Eurasian bottleneck populations. This has now been fixed, but if you had
copied this model from the tutorial for your own analyses, you should
update your model code or use the implementation that has been verified in
{ref}`stdpopsim project <stdpopsim:sec_introduction>`. See
[here](https://github.com/jeromekelleher/msprime-model-errors) for more
details on the faulty model and its likely effects on downstream analyses.
:::

```{code-cell}
import math

# Times are provided in years, so we convert into generations.
generation_time = 25
T_OOA = 21.2e3 / generation_time
T_AMH = 140e3 / generation_time
T_ANC = 220e3 / generation_time
# We need to work out the starting (diploid) population sizes based on
# the growth rates provided for these two populations
r_CEU = 0.004
r_CHB = 0.0055
N_CEU = 1000 / math.exp(-r_CEU * T_OOA)
N_CHB = 510 / math.exp(-r_CHB * T_OOA)

demography = msprime.Demography()
demography.add_population(
    name="YRI",
    description="Yoruba in Ibadan, Nigeria",
    initial_size=12300,
)
demography.add_population(
    name="CEU",
    description=(
        "Utah Residents (CEPH) with Northern and Western European Ancestry"
    ),
    initial_size=N_CEU,
    growth_rate=r_CEU,
)
demography.add_population(
    name="CHB",
    description="Han Chinese in Beijing, China",
    initial_size=N_CHB,
    growth_rate=r_CHB,
)
demography.add_population(
    name="OOA",
    description="Bottleneck out-of-Africa population",
    initial_size=2100,
)
demography.add_population(
    name="AMH", description="Anatomically modern humans", initial_size=12300
)
demography.add_population(
    name="ANC",
    description="Ancestral equilibrium population",
    initial_size=7300,
)

# Set the migration rates between extant populations
demography.set_symmetric_migration_rate(["CEU", "CHB"], 9.6e-5)
demography.set_symmetric_migration_rate(["YRI", "CHB"], 1.9e-5)
demography.set_symmetric_migration_rate(["YRI", "CEU"], 3e-5)

demography.add_population_split(
    time=T_OOA, derived=["CEU", "CHB"], ancestral="OOA"
)
demography.add_symmetric_migration_rate_change(
    time=T_OOA, populations=["YRI", "OOA"], rate=25e-5
)
demography.add_population_split(
    time=T_AMH, derived=["YRI", "OOA"], ancestral="AMH"
)
demography.add_population_split(time=T_ANC, derived=["AMH"], ancestral="ANC")

demography.debug()

```

```{code-cell}
:tags: [remove-cell]

# Make sure we don't insert any errors in this version of the model.
# Once we have an equivalent version in stdpopsim we can update here
# to compare against that instead and remove the Demography._ooa_model()

demography.assert_equal(msprime.Demography._ooa_model())

```
### Admixture


```{code-cell}
import math

# Implementation of the stdpopsim AmericanAdmixture_4B11 model.
T_OOA = 920
N_EUR = 34039
r_EUR = 0.0038
N_EAS = 45852
r_EAS = 0.0048
T_ADMIX = 12
N_ADMIX = 54664
r_ADMIX = 0.05

demography = msprime.Demography()
demography.add_population(
    name="AFR", description="African population", initial_size=14474
)
demography.add_population(
    name="EUR",
    description="European population",
    initial_size=N_EUR,
    growth_rate=r_EUR,
)
demography.add_population(
    name="EAS",
    description="East Asian population",
    initial_size=N_EAS,
    growth_rate=r_EAS,
)
demography.add_population(
    name="ADMIX",
    description="Admixed America",
    initial_size=N_ADMIX,
    growth_rate=r_ADMIX,
)
demography.add_admixture(
    T_ADMIX,
    derived="ADMIX",
    ancestral=["AFR", "EUR", "EAS"],
    proportions=[1 / 6, 2 / 6, 3 / 6],
)
demography.add_population(
    name="OOA",
    description="Bottleneck out-of-Africa population",
    initial_size=1861,
)
demography.add_population(
    name="AMH", description="Anatomically modern humans", initial_size=14474
)
demography.add_population(
    name="ANC",
    description="Ancestral equilibrium population",
    initial_size=7310,
)
demography.set_symmetric_migration_rate(["AFR", "EUR"], 2.5e-5)
demography.set_symmetric_migration_rate(["AFR", "EAS"], 0.78e-5)
demography.set_symmetric_migration_rate(["EUR", "EAS"], 3.11e-5)

demography.add_population_split(T_OOA, derived=["EUR", "EAS"], ancestral="OOA")
demography.add_symmetric_migration_rate_change(
    time=T_OOA, populations=["AFR", "OOA"], rate=15e-5
)
demography.add_population_split(2040, derived=["OOA", "AFR"], ancestral="AMH")
demography.add_population_split(5920, derived=["AMH"], ancestral="ANC")

demography.debug()
```

```{code-cell}
:tags: [remove-cell]

# Make sure we don't insert any errors in this version of the model.
# Once we have an equivalent version in stdpopsim we can update here
# to compare against that instead and remove the local model

demography.assert_equal(msprime.Demography._american_admixture_model())

```

### Trunk population models

For many empirical models we want to sequentially merge derived populations
into a "trunk" population.

```{eval-rst}
.. todo::
    It would be better to illustrate the merging-to-a-trunk idea with the OOA
    model above, since that's what we already have an illustration of above
    and it's a good bit simpler.
```

```{code-cell}
import math

# Times are provided in years, so we convert into generations.
generation_time = 29
T_OOA = 36_000 / generation_time
T_AMH = 60_700 / generation_time
T_ANC = 300_000 / generation_time
T_ArchaicAFR = 499_000 / generation_time
T_Neanderthal = 559_000 / generation_time
T_archaic_migration_start = 18_700 / generation_time
T_archaic_migration_end = 125_000 / generation_time

# We need to work out the starting (diploid) population sizes based on
# the growth rates provided for these two populations
r_CEU = 0.00125
r_CHB = 0.00372
N_CEU = 2300 / math.exp(-r_CEU * T_OOA)
N_CHB = 650 / math.exp(-r_CHB * T_OOA)

demography = msprime.Demography()
# This is the "trunk" population that we merge other populations into
demography.add_population(
    name="AFR",
    description="African population",
    initial_size=13900,
    initially_active=True,
)
demography.add_population(
    name="CEU",
    description=(
        "Utah Residents (CEPH) with Northern and Western European Ancestry"
    ),
    initial_size=N_CEU,
    growth_rate=r_CEU,
)
demography.add_population(
    name="CHB",
    description="Han Chinese in Beijing, China",
    initial_size=N_CHB,
    growth_rate=r_CHB,
)
demography.add_population(
    name="Neanderthal",
    description="Putative Neanderthals",
    initial_size=3600,
)
demography.add_population(
    name="ArchaicAFR",
    description="Putative Archaic Africans",
    initial_size=3600,
)
demography.add_population(
    name="OOA",
    description="Bottleneck out-of-Africa population",
    initial_size=880,
)

# Set the migration rates between extant populations
demography.set_symmetric_migration_rate(["CEU", "CHB"], 11.3e-5)
demography.set_symmetric_migration_rate(["AFR", "CEU"], 2.48e-5)

demography.add_symmetric_migration_rate_change(
    T_archaic_migration_start, ["CEU", "Neanderthal"], 0.825e-5
)
demography.add_symmetric_migration_rate_change(
    T_archaic_migration_start, ["CHB", "Neanderthal"], 0.825e-5
)
demography.add_symmetric_migration_rate_change(
    T_archaic_migration_start, ["ArchaicAFR", "AFR"], 1.98e-5
)
demography.add_migration_rate_change(T_archaic_migration_end, rate=0)

demography.add_population_split(
    time=T_OOA, derived=["CEU", "CHB"], ancestral="OOA"
)
demography.add_symmetric_migration_rate_change(
    time=T_OOA, populations=["AFR", "OOA"], rate=52.2e-5
)
demography.add_symmetric_migration_rate_change(
    time=T_OOA, populations=["OOA", "Neanderthal"], rate=0.825e-5
)
demography.add_population_split(time=T_AMH, derived=["OOA"], ancestral="AFR")
demography.add_symmetric_migration_rate_change(
    T_AMH, ["ArchaicAFR", "AFR"], 1.98e-5
)
demography.add_population_parameters_change(
    time=T_AMH, population="AFR", initial_size=13900
)
demography.add_population_parameters_change(
    time=T_ANC, population="AFR", initial_size=3600
)
demography.add_population_split(
    time=T_ArchaicAFR, derived=["ArchaicAFR"], ancestral="AFR"
)
demography.add_population_split(
    time=T_Neanderthal, derived=["Neanderthal"], ancestral="AFR"
)
demography.sort_events()

demography.debug()
```

```{code-cell}
:tags: [remove-cell]

# Make sure we don't insert any errors in this version of the model.
# Once we have an equivalent version in stdpopsim we can update here
# to compare against that instead and remove the local model

demography.assert_equal(msprime.Demography._ooa_archaic_model())

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
