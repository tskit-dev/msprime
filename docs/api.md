# API Reference


## Summary 

### Simulation functions
 
```{eval-rst}
.. currentmodule:: msprime

.. autosummary::

  sim_ancestry
  sim_mutations
  
```

### Ancestry Models


```{eval-rst}
.. autosummary::

  StandardCoalescent
  SmcApproxCoalescent
  SmcPrimeApproxCoalescent
  DiscreteTimeWrightFisher
  BetaCoalescent
  DiracCoalescent
  SweepGenicSelection
```

### Mutation Models

```{eval-rst}
.. autosummary::

  BinaryMutationModel
  JC69MutationModel
  HKYMutationModel
  F84MutationModel
  GTRMutationModel
  BLOSUM62MutationModel
  PAMMutationModel
  MatrixMutationModel
  InfiniteAllelesMutationModel
  SLiMMutationModel
```

### Demography

```{eval-rst}
.. autosummary::

  Demography
  DemographyDebugger
  Population
  MassMigration 
  PopulationParametersChange
  MigrationRateChange
```


```{eval-rst}
.. todo:: Section for utililies like RateMap
```

## Ancestry

```{eval-rst}
.. autofunction:: msprime.sim_ancestry
```

```{eval-rst}
.. autoclass:: msprime.SampleSet
```

```{eval-rst}
.. autoclass:: msprime.AncestryModelChange
```

### Models

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
.. autoclass:: msprime.BetaCoalescent
```

```{eval-rst}
.. autoclass:: msprime.DiracCoalescent
```

```{eval-rst}
.. autoclass:: msprime.SweepGenicSelection
```

## Mutations


### Models
```{eval-rst}
.. autoclass:: msprime.BinaryMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.JC69MutationModel()
```

```{eval-rst}
.. autoclass:: msprime.HKYMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.F84MutationModel()
```

```{eval-rst}
.. autoclass:: msprime.GTRMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.BLOSUM62MutationModel()
```

```{eval-rst}
.. autoclass:: msprime.PAMMutationModel()

```

```{eval-rst}
.. autoclass:: msprime.MatrixMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.InfiniteAllelesMutationModel()
```

```{eval-rst}
.. autoclass:: msprime.SLiMMutationModel()
```


(sec_api_node_flags)=

## Node flags

For standard coalescent simulations, all samples are marked with the
{data}`tskit.NODE_IS_SAMPLE` flag; internal nodes all have a flags value of 0.

<!---
todo link these up with the examples sections below where they are used.
-->

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
Can be used in combination with the ``record_migrations``.
Only present if the ``record_full_arg`` option is
specified.

```

```{data} msprime.NODE_IS_CEN_EVENT

The node was created by a :class:`msprime.CensusEvent`.

```

## Utilities


### Rate maps

```{eval-rst}
.. todo:: Add some high-level content here.
```

```{eval-rst}
.. autoclass:: msprime.RateMap
    :members:
```

## Demography


```{eval-rst}
.. autoclass:: msprime.Demography
    :members:
```

```{eval-rst}
.. autoclass:: msprime.Population
    :members:
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
.. autoclass:: msprime.DemographyDebugger
    :members:

```

## Likelihoods

```{eval-rst}
.. autofunction:: msprime.log_arg_likelihood
```

```{eval-rst}
.. autofunction:: msprime.log_mutation_likelihood
```

## Mutations

```{eval-rst}
.. autofunction:: msprime.sim_mutations
```
