% This is the recommended link for external projects to link in to
% msprime, so don't change it!

(sec_intro)=
# Introduction

This is the manual for msprime, a population genetics simulator
based on [tskit](https://tskit.dev). Msprime can simulate
{ref}`ancestral histories<sec_ancestry>` for a sample of individuals,
consistent with a given {ref}`demography<sec_demography>` under a
range of different models and evolutionary processes. Msprime can
also simulate {ref}`mutations <sec_mutations>` on a given ancestral
history (which can be produced by msprime ancestry simulations
or other programs supporting [tskit](https://tskit.dev)) under
a variety of diffent {ref}`models <sec_mutations_models>`
of genome sequence evolution.

Besides this manual, there are a number of other resources
available for learning about tskit and msprime:

- The [tskit tutorials](https://tskit.dev/tutorials) site contains
  in-depth tutorials on different aspects of msprime simulations
  as well as how to analyse simulated tree sequences.

- The book chapter
  [Coalescent simulation with msprime](https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_9)
  is a comprehensive introduction to
  running coalescent simulations with msprime, and provides many examples
  of how to run and use coalescent simulations. **Note however** that
  the chapter uses the deprecated {ref}`legacy 0.x API<sec_legacy_0x>`,
  and so does not follow current best practices.


```{eval-rst}
.. todo::
    Some top level content
```

```{tableofcontents}
```

