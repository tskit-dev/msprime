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


(sec_likelihood)=

# Computing likelihoods

In general, the likelihood function is defined as the probability (or
probability density) assigned to an observed outcome by a particular model.
It is often alternately referred to as the sampling probability, and is a
cornerstone of many statistical inference procedures.

`msprime` provides functions for evaluating two sampling probabilities:
that of a stored tree sequence for a given diploid population size
{math}`N_e` and per-link, per-generation recombination probability {math}`r`
under the standard ancestral recombination graph (ARG); and that of a pattern
of mutations given a tree sequence and per-site, per-generation mutation
probability {math}`\mu` under the infinite sites model.

:::{note}
Generic analysis tools for arbitrary tree sequences, allowing efficient calculation
of genetic diversity, allele frequency spectra, F statistics, and so on, are
provided within [tskit](https://tskit.dev/tskit/docs/stable/) and described in the
[statistics documentation](https://tskit.dev/tskit/docs/stable/stats.html).
The likelihood functions below are provided by `msprime` because they are specific
to the ARG and infinite sites models. 
:::

:::{warning}
Evaluating these log likelihoods requires knowledge of a whole ARG,
rather than just the more compact tree sequence representation. Additionally
the log likelihood implementations are only available for continuous genomes.
Hence, the underlying tree sequence has to conform to the
`record_full_arg = True` and `discrete_genome = False` options of the
{func}`.sim_ancestry` function.
:::

---

## Quick reference

{func}`.log_arg_likelihood`                                               
: Log likelihood of an ARG topology and branch lengths                    

{func}`.log_mutation_likelihood`                                          
: Log likelihood of a pattern of mutations arising from a given ARG       

---

(sec_likelihood_topology)=

## ARG sampling probability

The following example simulates an ARG with 5 diploid samples and evaluates
the likelihood of the realisation for four combinations of parameters.
Note that one combination coincides with the parameters used to simulate
the ARG, while the other combinations are quite different.

```{code-cell}
    import msprime

    ts = msprime.sim_ancestry(
        5, recombination_rate=1, record_full_arg=True,
        sequence_length=1, discrete_genome=False, random_seed=42)
    print(msprime.log_arg_likelihood(ts, recombination_rate=0, Ne=1))
    print(msprime.log_arg_likelihood(ts, recombination_rate=0.1, Ne=1))
    print(msprime.log_arg_likelihood(ts, recombination_rate=1, Ne=1))
    print(msprime.log_arg_likelihood(ts, recombination_rate=10, Ne=10))
```

In this example, the simulated ARG contains at least one recombination,
which is an event of probability 0 when the `recombination_rate = 0`.
Hence, the log likelihood for recombination rate zero returns a numerical
representation of negative infinity, i.e. the logarithm of zero. The other
three combinations of parameters all result in positive likelihoods,
or finite log likelihoods.

The data was generated from a recombination rate of one and a default
population size of one, and these parameters give rise to a relatively high
log likelihood. The same ARG would have been an unlikely realisation under
very different parameters, which thus result in more negative values of the
log likelihood.

## Mutation sampling probability

The next example adds random mutations to the tree sequence
generated above in {ref}`sec_likelihood_topology` and evaluates the
unnormalised log probability of the mutation realisation, given the
tree sequence and a prescribed mutation rate.

```{code-cell}
    ts = msprime.mutate(ts, rate=1, random_seed=42)
    print(msprime.log_mutation_likelihood(ts, mutation_rate=0))
    print(msprime.log_mutation_likelihood(ts, mutation_rate=0.1))
    print(msprime.log_mutation_likelihood(ts, mutation_rate=1))
    print(msprime.log_mutation_likelihood(ts, mutation_rate=10))
```

Since there is at least one mutation in the realisation, its probability
given `mutation_rate = 0` is 0, resulting in a log likelihood of negative
infinity. The mutation realisation is a typical outcome when
`mutation_rate = 1`, which is the value used to simulate it, and which thus
results in a relatively high log likelihood. Mutation rates which are
significantly higher or lower result in more negative log likelihoods
because generating the same realisation using those rates is unlikely.
