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

import io
import msprime
import tskit
from IPython.display import SVG, set_matplotlib_formats
from matplotlib import pyplot as plt
import numpy as np
import networkx as nx

set_matplotlib_formats("svg")
```

(sec_pedigrees)=

# Pedigrees

Pedigrees describe parent-offspring relationships between individual organisms,
and can be provided as input to contrain simulations of genetic ancestry
(see the {ref}`sec_ancestry_models_fixed_pedigree` section for details).
In this section we describe the data structures used to encode pedigree
data in msprime, and the utilities used to create input pedigrees.

:::{todo}
This page is incomplete. Needs another pass through.
https://github.com/tskit-dev/msprime/issues/2008
:::

(sec_pedigrees_encoding)=

## Pedigree encoding

Msprime uses the tskit
{ref}`node<tskit:sec_node_table_definition>` and
{ref}`individual<tskit:sec_individual_table_definition>` tables to encode pedigree
information.
Each {ref}`individual<sec_data_model_definitions_individual>`
is defined by a row in the individual table, and
this row index is the individual's **id**. Pedigree relationships
are defined using the ``parents`` column, which contains the IDs
of an individual's parents. Further information is stored in the
{ref}`nodes<sec_data_model_definitions_node>` associated with the individual.

:::{note}
It is important to note that sample status, time and population information
is associated with an individual's nodes, which is an artefact of
tskit's node-centric design.
:::

It is possible to specify a pedigree model directly using tskit APIs,
but it is simpler to use the {class}`.PedigreeBuilder` utility
class (or the simplified text {ref}`sec_pedigrees_file_format`).
In the following example we build a simple trio:

```{code-cell}
pb = msprime.PedigreeBuilder()
mom_id  = pb.add_individual(time=1)
dad_id  = pb.add_individual(time=1)
pb.add_individual(time=0, parents=[mom_id, dad_id], is_sample=True)
pedigree = pb.finalise()
# TODO replace with display(pedigree) when its implemented in tskit
# https://github.com/tskit-dev/tskit/issues/2093
print(pedigree)
```

The pedigree returned by the {meth}`~.PedigreeBuilder.finalise` method
contains the pedigree information defined by the calls to
{meth}`~.PedigreeBuilder.add_individual`. For this trio, we began by adding
the parents, and because they are founders, we don't provide any information
about their parents. Each call to returns the integer ID of
newly added individual.

### Requirements



(sec_pedigrees_file_format)=

## Simplified file format

The methods described in {ref}`sec_pedigrees_encoding` are general and
allow for arbitrary metadata to be associated with individuals. It
is often convenient to work with a text based representation of
the pedigree, which is supported by the {func}`.parse_pedigree`.

### Basic structure

```{code-cell}
txt = """\
# id parent0 parent1 time
child mom dad 0
mom NA NA 1
dad NA NA 1
"""
pedigree = msprime.parse_pedigree(io.StringIO(txt))
display(pedigree.individuals)
display(pedigree.nodes)
```


(sec_pedigrees_visualisation)=

## Visualising pedigrees

```{code-cell}
import io
txt = """\
# id parent0 parent1 time
child mom dad 0
mom NA NA 1
dad NA NA 1
"""
pedigree = msprime.parse_pedigree(io.StringIO(txt))
display(pedigree.individuals)

display(pedigree.nodes)
```

```{code-cell}
pedigree.sequence_length = 100
ts = msprime.sim_ancestry(initial_state=pedigree, model="fixed_pedigree",
    random_seed=1)

SVG(ts.draw_svg())
```



## Visualising pedigrees

```{code-cell}


def draw_pedigree(ped_ts):

    G = nx.DiGraph()
    for ind in ped_ts.individuals():
        time = ped_ts.node(ind.nodes[0]).time
        G.add_node(ind.id, time=time)
        for p in ind.parents:
            if p != tskit.NULL:
                G.add_edge(ind.id, p)
    pos = nx.multipartite_layout(G, subset_key="time", align="horizontal")
    nx.draw_networkx(G, pos, with_labels=True)

tables = msprime.pedigrees.sim_pedigree(
    population_size=20, end_time=3, num_samples=5, direction="backward", random_seed=1)
tables.individuals
tables.sequence_length = 1


draw_pedigree(tables.tree_sequence())

```


