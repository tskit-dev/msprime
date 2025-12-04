% This is the recommended link for external projects to link in to
% msprime, so don't change it!

(sec_intro)=
# Introduction

This is the manual for {program}`msprime`, a population genetics
simulator of ancestry and DNA sequence evolution based on
[tskit](https://tskit.dev). {program}`msprime` can simulate
{ref}`ancestral histories<sec_ancestry>` for a sample of individuals,
consistent with a given {ref}`demography<sec_demography>` under a
range of different models and evolutionary processes. It
can also simulate {ref}`mutations <sec_mutations>` on a given ancestral
history (which can be produced by {program}`msprime` ancestry
simulations or other programs supporting [tskit](https://tskit.dev))
under a variety of different {ref}`models <sec_mutations_models>`
of genome sequence evolution.

Besides this manual, there are a number of other resources
available for learning about {program}`tskit` and {program}`msprime`:

- The [tskit tutorials](https://tskit.dev/tutorials) site provides an
  introduction to the tree sequence format used by {program}`msprime`,
  as well as in-depth tutorials
  on more {ref}`advanced msprime topics<tutorials:sec_msprime>` and
  examples of how to {ref}`analyse<tutorials:sec_analysing_tree_sequences>`
  simulated {program}`tskit` tree sequences.

  - For a formal introduction to tree sequences and how they represent
    ancestral recombination graphs,
    see the [2024 Genetics paper](https://doi.org/10.1093/genetics/iyae100).

- Our [Discussions board](https://github.com/tskit-dev/msprime/discussions)
  is a great place to ask questions like "how do I do X" or "what's the best
  way to do Y". Please make questions as clear as possible, and be respectful,
  helpful, and kind.

- The book chapter
  [Coalescent simulation with msprime](https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_9)
  is a comprehensive introduction to running coalescent simulations with
  {program}`msprime`, and provides many examples
  of how to run and use coalescent simulations. **Note however** that
  the chapter uses the deprecated {ref}`legacy 0.x API<sec_legacy_0x>`,
  and so does not follow current best practices.

- If you would like to understand more about the underlying algorithms
  for {program}`msprime`, please see the
  [2016 PLoS Computational Biology paper](https://doi.org/10.1371/journal.pcbi.1004842).
  For more information on the {ref}`sec_ancestry_models_dtwf` model,
  please see the [2020 PLoS Genetics paper](https://doi.org/10.1371/journal.pgen.1008619).
  For a formal introduction to recent capabilities, see the
  [2022 Genetics paper](https://doi.org/10.1093/genetics/iyab229)
  describing the {program}`msprime` 1.0 release.

```{important}
If you use {program}`msprime` in your work, please remember to
cite it appropriately: see the {ref}`citations<sec_citation>` page
for details.
```

## Contents

```{tableofcontents}
```
