
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

This section is to help 0.x users of the get up to speed quickly, summarising the new
APIs and their main differences to what you are used to.

The main change is that there are two new functions, {func}`.sim_ancestry` and
{func}`.sim_mutations` which correspond to the 0.x functions {func}`.simulate`
and {func}`.mutate`. The 0.x functions are **deprecated** but **will continue
to be supported indefinitely**.

## Backwards compatibility

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

## Ancestry

The new {func}`.sim_ancestry` function replaces the 0.x {func}`.simulate`
function and is very similar. There are some important differences though:

* The {}`samples` argument now refers to the **number of individuals**
  rather than **the number of nodes** (i.e. monoploid genomes).
  Because the default {ref}`ploidy <sec_ancestry_ploidy>`
  is 2 (see the next point) the upshot is that `sim_ancestry(2)` will
  result in a tree sequence with *four* sample nodes, not two. (It is
  possible to override this behaviour using the list of {class}`.SampleSet`
  objects argument to `samples`.)
* There is now a {ref}`sec_ancestry_ploidy` argument, which has
  two effects:

  1. Sets the default number of sample nodes per *individual*
  2. Changes the timescale of the coalescent process (TODO link to a section
     that explains this effect.) By default `ploidy` is 2 and
     time is scaled scaled in units of 4N generations, which is the same as
     msprime 0.x.
* Rather than two arguments `num_samples` and `samples`, the
  {func}`.sim_ancestry` function has a single argument `samples` which
  has different behaviour depending on the type of arguments provided.
  See {ref}`sec_ancestry_samples` for details.

  Note in particular that a list of `Sample` objects is **not** supported.
* Similarly, there is now one argument `recombination_rate` which can
  be either a single value or a {class}`.RateMap` object. Note that the
  0.x {class}`.RecombinationMap` is deprecated and not supported as input
  to {func}`.sim_ancestry`. See {ref}`sec_ancestry_recombination` for more
  details.
* Simulations are peformed on a **discrete** genome by default. To get the
  0.x behaviour of a continuous genome, set `discrete_genome=False`.
  See {ref}`sec_ancestry_discrete_genome` for more details.
* The `from_ts` argument used has been renamed to `initial_state` and
  accepts either a {class}`tskit.TableCollection` or {class}`tskit.TreeSequence`
  argument. See {ref}`sec_ancestry_initial_state` for details.
* There is **no** `mutation_rate` argument to {func}`.sim_ancestry`: use
  {func}`.sim_mutations` instead.
* The `population_configurations`, `migration_matrix` and `demographic_events`
  arguments have been replace with a single argument `demography`, which must take
  a {class}`.Demography` instance. (See the next section for more details.)

## Demography

* A new {class}`.Demography` object has been added for version 1.0 which
  encapsulates the functionality needed to define and debug demographic models
  in msprime. Demographic models can only be specified to `sim_ancestry`
  using an instance of this class.
* It is easy to create a {class}`.Demography` from the 0.x
  `population_configurations`, `migration_matrix` and `demographic_events`
  values using the {meth}`.Demography.from_old_style` method.
* The {class}`.DemographyDebugger` class should no longer be instantiated
  directly; instead use the {meth}`.Demography.debug` method.

## Mutations

* For symmetry with the {func}`.sim_ancestry` function, there is now a {func}`.sim_mutations`
  function. The 0.x {func}`.mutate` function is **deprecated**.
* The {func}`.sim_mutations` function works on a **discrete** genome by default.
* There are now also many new mutation models supported by :func:.sim_mutations;
  see {ref}`sec_mutations` for details. These are *not* supported in the deprecated
  {func}`.mutate` function.
* The simulated mutations now have a simulated ``time`` value, which specifies the 
  precise time that the mutation occurred. Note that this value is also provided in the
  returned tables for the deprecated ``simulate()`` and ``mutate()`` functions,
  which may lead to some compatibility issues. (We considered removing the simulated
  mutation times for these 0.x functions for strict compatibility, but this would
  have broken any code using the ``keep`` option in mutate.)

## Utilities

* The 0.x class {class}`.RecombinationMap` has been **deprecated** in favour of the new
  {class}`.RateMap`. This was to (a) generalise the interface to accomodate varying
  rates of mutation and gene conversion along the genome; and (b) convert to a
  more modern numpy-based API.


## API Reference


```{eval-rst}
.. autoclass:: msprime.PopulationConfiguration
```

```{eval-rst}
.. autoclass:: msprime.RecombinationMap
    :members:
```

```{eval-rst}
.. autofunction:: msprime.simulate()
```

```{eval-rst}
.. autoclass:: msprime.SimulationModelChange
```

```{eval-rst}
.. autofunction:: msprime.mutate
```
