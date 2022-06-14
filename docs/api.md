# API Reference

This lists the detailed reference documentation for the msprime Python API.
The reference documentation aims to be concise, precise and exhaustive;
as such, it is not the best place to start if you are new to
a particular piece of functionality. Please see the top-level documentation
for {ref}`sec_ancestry`, {ref}`sec_mutations` or {ref}`sec_demography`
for discussion and examples of individual features.

## Summary

```{eval-rst}
.. currentmodule:: msprime
```

### Ancestry

```{eval-rst}
.. autosummary::

  sim_ancestry
  SampleSet
  StandardCoalescent
  SmcApproxCoalescent
  SmcPrimeApproxCoalescent
  DiscreteTimeWrightFisher
  FixedPedigree
  BetaCoalescent
  DiracCoalescent
  SweepGenicSelection
```

### Mutations

```{eval-rst}
.. autosummary::

  sim_mutations
  JC69
  HKY
  F84
  GTR
  BLOSUM62
  PAM
  BinaryMutationModel
  MatrixMutationModel
  MicrosatMutationModel
  SMM
  TPM
  EL2
  InfiniteAlleles
  SLiMMutationModel
```

### Demography

```{eval-rst}
.. autosummary::

  Population
  Demography
  DemographyDebugger
```

### Utilities

```{eval-rst}
.. autosummary::

  RateMap
  PedigreeBuilder
  parse_pedigree


```

## Reference documentation

### Ancestry

```{eval-rst}
.. autofunction:: msprime.sim_ancestry
```

```{eval-rst}
.. autoclass:: msprime.SampleSet
    :members:
```

#### Models

```{eval-rst}
.. autoclass:: msprime.AncestryModel
    :members:
```

```{eval-rst}
.. autoclass:: msprime.StandardCoalescent
```

```{eval-rst}
.. autoclass:: msprime.SmcApproxCoalescent
```

```{eval-rst}
.. autoclass:: msprime.SmcPrimeApproxCoalescent
```

```{eval-rst}
.. autoclass:: msprime.DiscreteTimeWrightFisher
```

```{eval-rst}
.. autoclass:: msprime.FixedPedigree
```

```{eval-rst}
.. autoclass:: msprime.BetaCoalescent
```

```{eval-rst}
.. autoclass:: msprime.DiracCoalescent
```

```{eval-rst}
.. autoclass:: msprime.SweepGenicSelection
```

### Mutations

```{eval-rst}
.. autofunction:: msprime.sim_mutations
```

#### Models

```{eval-rst}
.. autoclass:: msprime.MutationModel
```

```{eval-rst}
.. autoclass:: msprime.MatrixMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.MicrosatMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.SMM()
```

```{eval-rst}
.. autoclass:: msprime.TPM()
```

```{eval-rst}
.. autoclass:: msprime.EL2()
```

```{eval-rst}
.. autoclass:: msprime.BinaryMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.JC69()
```

```{eval-rst}
.. autoclass:: msprime.HKY()
```

```{eval-rst}
.. autoclass:: msprime.F84()
```

```{eval-rst}
.. autoclass:: msprime.GTR()
```

```{eval-rst}
.. autoclass:: msprime.BLOSUM62()
```

```{eval-rst}
.. autoclass:: msprime.PAM()

```

```{eval-rst}
.. autoclass:: msprime.InfiniteAlleles()
```

```{eval-rst}
.. autoclass:: msprime.SLiMMutationModel()
```

(sec_api_node_flags)=

### Node flags

In the tskit {ref}`tskit:sec_node_table_definition` node flags specify
particular properties about nodes. Msprime follows the standard approach
of setting the {data}`tskit.NODE_IS_SAMPLE` flag for all sample nodes,
with all other nodes having a flags value of 0.

Msprime defines some extra flags that help us to identify particular
nodes in some situations:

```{data} msprime.NODE_IS_RE_EVENT

The node is an ARG recombination event. Each recombination event is marked
with two nodes, one identifying the individual providing the genetic
material to the left of the breakpoint and the other providing the genetic
material the right. Only present if the ``record_full_arg`` option is
specified.

```

```{data} msprime.NODE_IS_CA_EVENT

The node is an ARG common ancestor event that did not result in
marginal coalescence. Only present if the ``record_full_arg`` option is
specified.

```

```{data} msprime.NODE_IS_MIG_EVENT

The node is an ARG migration event identifying the individual that migrated.
Can be used in combination with the ``record_migrations`` option.
Only present if the ``record_full_arg`` option is
specified.

```

```{data} msprime.NODE_IS_CEN_EVENT

The node was created by a census event. Please see the
{ref}`sec_ancestry_census_events` section for more details.

```

### Rate maps

```{eval-rst}
.. autoclass:: msprime.RateMap
    :members:
```

### Demography

```{eval-rst}
.. autoclass:: msprime.Demography()
    :members:
```

```{eval-rst}
.. autoclass:: msprime.Population
    :members:
```

```{eval-rst}
.. autoclass:: msprime.DemographyDebugger
    :members:

```

### Likelihoods

```{eval-rst}
.. autofunction:: msprime.log_arg_likelihood
```

```{eval-rst}
.. autofunction:: msprime.log_mutation_likelihood
```

### Pedigrees

```{eval-rst}
.. autoclass:: msprime.PedigreeBuilder
    :members:
```

```{eval-rst}
.. autofunction:: msprime.parse_pedigree
```
