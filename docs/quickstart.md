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

(sec_quickstart)=

# Quickstart

This page gives some simple examples of how to use the major features
of {program}`msprime`, with links to more detailed documentation
and tutorial content.

See the {ref}`sec_installation` page for instructions on installing 
{program}`msprime` (short version: ``pip install msprime`` or 
``conda install -c conda-forge msprime`` will work for most users).

(sec_quickstart_ancestry)=

## Ancestry

Msprime simulates ancestral histories for a set of sample genomes 
using backwards-in-time population genetic models.
Here we run a simple simulation of a short recombining sequence under
human-like parameters:

```{code-cell}

    import msprime
    from IPython.display import SVG

    # Simulate an ancestral history for 3 diploid samples under the coalescent
    # with recombination on a 5kb region with human-like parameters.
    ts = msprime.sim_ancestry(
        samples=3, 
        recombination_rate=1e-8, 
        sequence_length=5_000,
        population_size=10_000, 
        random_seed=123456)
    # Visualise the simulated ancestral history.
    SVG(ts.draw_svg())
```

In this example we simulate the ancestral history of three diploid 
individuals (see {ref}`sec_ancestry_samples` and {ref}`sec_ancestry_ploidy`)
for a 5kb sequence with a 
[recombination](<https://en.wikipedia.org/wiki/Genetic_recombination>)
rate of {math}`10^{-8}` 
(see {ref}`sec_ancestry_genome_properties`)
from a population with a constant size of 10,000 (see 
the {ref}`sec_quickstart_demography` section below)
under the default 
[coalescent](<https://en.wikipedia.org/wiki/Coalescent_theory>)
ancestry model (see the {ref}`sec_ancestry_models` for details on 
other available models).
To ensure that 
the output of this example is predictable, we set a random seed 
(see {ref}`sec_ancestry_random_seed`).

When recombination is present, the ancestry of a sample of DNA sequences
cannot be represented by a single genealogical tree relating the 
samples to their genetic ancestors; there is instead
a *sequence* of highly correlated trees along the genome.
The result of our simulation is therefore a [tree sequence](https://tskit.dev)
object from the {ref}`tskit <tskit:sec_introduction>` library,
which provides a rich suite of operations for 
analysing these genealogical histories: see the 
{ref}`tutorials:sec_tskit_getting_started` tutorial for help. 
In this example we show a visualisation
of the four different trees along the 5kb region 
(see the {ref}`tutorials:sec_tskit_viz` tutorial for more 
examples).  Because we have specified three diploid sample 
*individuals*, each of these trees has 6 "sample" nodes 
(the "leaves" or "tips"), because each diploid individual
has two monoploid genomes (see {ref}`sec_ancestry_samples`).

See the {ref}`sec_ancestry` section for more details on 
ancestry simulations.

## Mutations

```{eval-rst}
.. todo:: This is a WIP
```

The {func}`.sim_ancestry` function generates a simulated ancestral
history for some samples. 
If we want [genome sequence](<https://en.wikipedia.org/wiki/Genome>)
we must also simulate some
[mutations](<https://en.wikipedia.org/wiki/Mutation>) on these trees.
However, it's important to note that it's not always necessary to 
simulate mutations in order to use the simulations; often, it's 
better *not to*; see the 
% TODO enable this once the tutorials build is fixed
% {ref}`tutorials:sec_tskit_no_mutations` 
tutorial for more information.

```{code-cell}
mutated_ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=54321)
mutated_ts.tables.sites
```

```{code-cell}
mutated_ts.tables.mutations
```

```{code-cell}
SVG(mutated_ts.draw_svg())
```

```{code-cell}
for variant in mutated_ts.variants():
    print(variant)
```

```{eval-rst}
.. todo:: Not sure how much detail we want to get into here. The salient points
    we want to get across are that mutations are not automatically part of the 
    simulation, that we only output sites with mutations and we have efficient
    ways to work with the results. A lot of this should be pointing to the 
    "getting started with tskit" tutorial.

```

```{eval-rst}
.. todo:: List of pointers to the relevant sections of the documentation.

```

(sec_quickstart_demography)=

## Demography

By default ancestry simulations assume an extremely simple
population structure in which a single randomly mating population
of a fixed size exists for all time. For most simulations this
is an unrealistic assumption, and so msprime provides a way
to describe more complex demographic histories.

```{code-cell}

    # Create a 1D stepping stone model of demograpy
    demography = msprime.Demography.stepping_stone_model([100] * 10, migration_rate=0.1)
    # Take one diploid sample each from the first and last demes
    samples = {0: 1, 9: 1}
    # Simulate an ancestral history for this demography and sample.
    ts = msprime.sim_ancestry(samples=samples, demography=demography)
    ts.tables.nodes
```

```{eval-rst}
.. todo:: Links into more detailed documentation
```

