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


(sec_legacy_0x)=
# Legacy (version 0.x) APIs


```{eval-rst}
.. todo:: This page is under construction, and needs work to ensure
    that the target audiences are able to find the content they need.
```


```{eval-rst}
.. todo:: An overview of why we added the new APIs, what they are
  and some reassurance that old code will continue to work.
```


With a few exceptions, you **should not need to change your code** and
it should remain working indefinitely.


## Upgrading code

This section is to help legacy 0.x users of msprime get up to speed quickly, summarising
the new APIs and their main differences to what you are used to in the 0.x versions.

The main change is that there are two new functions, {func}`.sim_ancestry` and
{func}`.sim_mutations` which correspond to the 0.x functions {func}`.simulate`
and {func}`.mutate`. The 0.x functions are **deprecated** but **will continue
to be supported indefinitely**.

### Backwards compatibility

All existing simulations should work as before, *except* for simulations relying on
the detailed properties of RecombinationMaps. If your code uses the `num_loci`
property, then it may need to be updated. The reason for this is that `msprime`
has changed to simulate directly in physical coordinates internally (which has
greatly simplified the code and solved many thorny issues) and this is fundamentally
incompatible with the approach taken in 0.x. In the vast majority of cases, this
will have no effect.

If you are using `num_loci` to simulate a discrete genome, it may be simplest to
convert your code to use the new {func}`.sim_ancestry` method. If you were following
a recipe to simulate multiple chromosomes under the DTWF model, please see
the {ref}`updated recipe <sec_ancestry_multiple_chromosomes>`.

### Ancestry

The new {func}`.sim_ancestry` function replaces the 0.x {func}`.simulate`
function and is very similar. There are some important differences though:

* The `samples` parameter now refers to the **number of individuals**
  rather than **the number of nodes** (i.e. monoploid genomes).
  Because the default {ref}`ploidy <sec_ancestry_ploidy>`
  is 2 (see the next point) the upshot is that `sim_ancestry(2)` will
  result in a tree sequence with *four* sample nodes, not two. (It is
  possible to override this behaviour using the list of {class}`.SampleSet`
  objects parameter to `samples`.)
* The `Ne` parameter in  0.x {func}`.simulate` function has been replaced
  with the `population_size` parameter.
* There is now a {ref}`sec_ancestry_ploidy` parameter, which has
  two effects:

  1. Sets the default number of sample nodes per *individual*
  2. Changes the timescale of the coalescent process (TODO link to a section
     that explains this effect.) By default `ploidy` is 2 and
     time is scaled scaled in units of 4N generations, which is the same as
     msprime 0.x.
* Rather than two parameters `num_samples` and `samples`, the
  {func}`.sim_ancestry` function has a single parameter `samples` which
  has different behaviour depending on the type of parameters provided.
  See {ref}`sec_ancestry_samples` for details.

  Note in particular that a list of `Sample` objects is **not** supported.
* Similarly, there is now one parameter `recombination_rate` which can
  be either a single value or a {class}`.RateMap` object. Note that the
  0.x {class}`.RecombinationMap` is deprecated and not supported as input
  to {func}`.sim_ancestry`. See {ref}`sec_ancestry_recombination` for more
  details.
* Simulations are peformed on a **discrete** genome by default. To get the
  0.x behaviour of a continuous genome, set `discrete_genome=False`.
  See {ref}`sec_ancestry_discrete_genome` for more details.
* The `from_ts` parameter used has been renamed to `initial_state` and
  accepts either a {class}`tskit.TableCollection` or {class}`tskit.TreeSequence`
  parameter. See {ref}`sec_ancestry_initial_state` for details.
* There is **no** `mutation_rate` parameter to {func}`.sim_ancestry`: use
  {func}`.sim_mutations` instead.
* The `population_configurations`, `migration_matrix` and `demographic_events`
  parameters have been replace with a single parameter `demography`, which must take
  a {class}`.Demography` instance. (See the next section for more details.)

### Demography

* A new {class}`.Demography` object has been added for version 1.0 which
  encapsulates the functionality needed to define and debug demographic models
  in msprime. Demographic models can only be specified to `sim_ancestry`
  using an instance of this class.
* It is easy to create a {class}`.Demography` from the 0.x
  `population_configurations`, `migration_matrix` and `demographic_events`
  values using the {meth}`.Demography.from_old_style` method.
* The {class}`.DemographyDebugger` class should no longer be instantiated
  directly; instead use the {meth}`.Demography.debug` method.

### Mutations

* For symmetry with the {func}`.sim_ancestry` function, there is now a {func}`.sim_mutations`
  function. The 0.x {func}`.mutate` function is **deprecated**.
* The {func}`.sim_mutations` function works on a **discrete** genome by default.
* There are now also many new mutation models supported by {func}`.sim_mutations`;
  see {ref}`sec_mutations` for details. These are *not* supported in the deprecated
  {func}`.mutate` function.
* The simulated mutations now have a simulated ``time`` value, which specifies the
  precise time that the mutation occurred. Note that this value is also provided in the
  returned tables for the deprecated ``simulate()`` and ``mutate()`` functions,
  which may lead to some compatibility issues. (We considered removing the simulated
  mutation times for these 0.x functions for strict compatibility, but this would
  have broken any code using the ``keep`` option in mutate.)

(sec_legacy_0x_genome_discretisation)=
### Genome discretisation

In msprime 0.x, recombination was implemented internally using a discrete
number of genetic loci. That is, the simulation was performed in
*genetic* coordinates, which were then mapped back to *physical* coordinates
at the end of simulation. This had the significant advantage that
recombination could be implemented during the simulation as a uniform process
over these discrete loci. However, it led to
a number of different numerical issues encountered when mapping back and
forth between genetic and physical coordinates as well as limiting
what could be done in terms of gene converion and other processes.
We therefore changed to using physical coordinates throughout the simulation
for msprime 1.0.

The simulations in 0.x and 1.x are almost entirely compatible and everything
should work as expected when running 0.x code on msprime 1.0 or later. However,
there is one (hopefully obscure) case in which code written for msprime 0.x
will no longer work.

The ``num_loci`` argument to the 0.x class {class}`.RecombinationMap`
was used to control the number of discrete genetic loci in the simulation. By
default, this was set to a large number ({math}`\sim 2^{32}`), effectively
giving a continuous coordinate space when mapped back into physical units.
By setting the ``num_loci`` equal
to the sequence length of the RecombinationMap, we could also specify
discrete physical loci. Specifying whether we simulate in discrete or continuous
genome coordinates is now done using the ``discrete_genome`` argument
to {func}`.sim_ancestry` (see the {ref}`sec_ancestry_discrete_genome`
section for more details). The {class}`.RateMap` class is now used to
specify varying rates of recombination along the genome and no longer
has any concept of genetic "loci" --- the choice of coordinate space
is now decoupled from our specification of the recombination process.

Both the cases of discrete and fully continuous genomes are well
supported in msprime 1.x  and so nearly all existing code
should continue to work as expected.
What is no longer supported is specifying the "granularity" of the
continuous space via the ``num_loci`` parameter, and if we try
to set ``num_loci`` to anything other than the sequence length
we get an error:

```{code-cell}
:tags: [raises-exception]
import msprime

# Here we try to make a sequence length of 10 with 5 discrete loci
recomb_map = msprime.RecombinationMap(positions=[0, 10], rates=[0.1, 0], num_loci=5)
```

If you get this error, please check whether specifying a
number of loci like this was actually what you intended. Almost
certainly you actually wanted to simulate a continuous genome
(omit the ``num_loci`` parameter) or a discrete genome
with the breaks occurring integer boundaries (set ``num_loci``
equal to the sequence length).

If not, please let us know your use case and we may be able
to accommodate it in the new code. Until then, you will need
to downgrade msprime to 0.7.x for your simulations to run.

## API Reference


### Ancestry

```{eval-rst}
.. autofunction:: msprime.simulate()
```

```{eval-rst}
.. autoclass:: msprime.PopulationConfiguration
```

```{eval-rst}
.. autoclass:: msprime.PopulationParametersChange

```

```{eval-rst}
.. autoclass:: msprime.MigrationRateChange

```

```{eval-rst}
.. autoclass:: msprime.MassMigration
```

```{eval-rst}
.. autoclass:: msprime.CensusEvent
```

```{eval-rst}
.. autoclass:: msprime.Sample
```

```{eval-rst}
.. autoclass:: msprime.SimulationModelChange
```

### Recombination maps

```{eval-rst}
.. autoclass:: msprime.RecombinationMap
    :members:
```

### Mutations

```{eval-rst}
.. autofunction:: msprime.mutate
```

```{eval-rst}
.. autoclass:: msprime.InfiniteSites
```


