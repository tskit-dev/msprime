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

{func}`.log_arg_likelihood`                
: Likelihood of an ARG topology            

{func}`.log_mutation_likelihood`           
: Likelihood of a set of mutations         


`msprime` provides the capability to evaluate two sampling probabilities:
that of a stored tree sequence for a given diploid effective population size
{math}`N_e` and per-link, per-generation recombination probability {math}`r`
under the standard ancestral recombination graph; and that of a pattern of
mutations given a tree sequence and per-site, per-generation mutation
probability {math}`\mu` under the infinite sites model. In both cases,
the tree sequence must conform to the `record_full_arg` option of
the {func}`.simulate` function.

The following sections illustrate the evaluation of these log likelihoods.

(sec_likelihood_topology)=

## Topology

The following example simulates a tree sequence with 5 diploid samples
and evaluates the likelihood of the realisation for various parameter
combinations.

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

## Mutations

The following example adds random mutations to the tree sequence
generated in {ref}`sec_likelihood_topology` and evaluates the
unnormalised log likelihood of the mutations given the tree sequence.

```{code-cell}
    ts = msprime.mutate(ts, rate=1, random_seed=42)
    print(msprime.log_mutation_likelihood(ts, mutation_rate=0))
    print(msprime.log_mutation_likelihood(ts, mutation_rate=0.1))
    print(msprime.log_mutation_likelihood(ts, mutation_rate=1))
    print(msprime.log_mutation_likelihood(ts, mutation_rate=10))
```


