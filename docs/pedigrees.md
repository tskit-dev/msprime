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

Pedigrees describe parent-offspring relationships between individuals,
and can be provided as input to constrain simulations of genetic ancestry
(see the {ref}`sec_ancestry_models_fixed_pedigree` section for details).
In this section we describe the data structures used to encode pedigrees
in msprime, and the utilities used to create input pedigrees.

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
It is important to note that sample status, time, and population information
is associated with an individual's nodes (see the {ref}`sec_ancestry_ploidy`
section for details), which is an artefact of tskit's node-centric design.
:::

It is possible to specify a pedigree model directly using tskit APIs,
but it is simpler to use the {class}`.PedigreeBuilder` utility
class (or the {ref}`sec_pedigrees_file_format`).
In the following example we build a simple trio:

```{code-cell}
pb = msprime.PedigreeBuilder()
mom_id = pb.add_individual(time=1)
dad_id = pb.add_individual(time=1)
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
about their parents. Each call to {meth}`~.PedigreeBuilder.add_individual`
returns the integer ID of newly added individual.
### Requirements

This section lists the detailed requirements of the low-level encoding
for the tskit node and individual tables.

:::{warning}
Most of this information is only
of interest if you are constructing these tables by hand. For most
purposes it's best to use either the {class}`.PedigreeBuilder`
API or the {ref}`sec_pedigrees_file_format`.
:::

- Each individual must have exactly two parents. Unknown parents are
  indicated using ``-1`` (thus, a founder individual will have
  ``parents=[-1, -1]``.
- Each individual must be associated with exactly two nodes.
- All nodes associated with an individual must have the same ``time``
  and ``population`` values. (These are referred to as the individual's
  time and population, as a shorthand.)
- An individual's time must be less than all its parent's times.


(sec_pedigrees_metadata)=

## Metadata

It is often useful to associate metadata with individuals when building
a pedigree. This allows us to either directly store information about the
individuals for later analysis, or to store an identifier to facilitate
matching with existing records.

In this example we update our simple trio by associating a string
identifier with each individual. We use the
{meth}`tskit.MetadataSchema.permissive_json` function here to create
a schema, as this is an easy and quick way to get a functional
metadata schema. Much more sophisticated approaches are possible, however: please
see the {ref}`tskit metadata documentation<tskit:sec_metadata>` for
more information.

```{code-cell}
pb = msprime.PedigreeBuilder(
    individuals_metadata_schema=tskit.MetadataSchema.permissive_json())
mom_id = pb.add_individual(time=1, metadata={"name": "mom"})
dad_id = pb.add_individual(time=1, metadata={"name": "dad"})
pb.add_individual(
    time=0, parents=[mom_id, dad_id], is_sample=True, metadata={"name": "child"})
pedigree = pb.finalise()
display(pedigree.individuals)
```

(sec_pedigrees_file_format)=

## Simplified file format

The methods described in {ref}`sec_pedigrees_encoding` are general and
allow for arbitrary metadata to be associated with individuals. It
is often convenient to work with a text based representation of
the pedigree, which is supported by the {func}`.parse_pedigree` function.
The {ref}`sec_pedigrees_file_format_definition` section defines the columns and
requirements for this file format.

:::{seealso}
See the
 {ref}`sec_pedigrees_file_format_basic_example` and subsequent sections
for examples.
:::


(sec_pedigrees_file_format_definition)=

### Format definition

This section describes the detailed rules for the pedigree text file format.

The first line of the file must be a header starting with a ``#`` character,
listing the included columns separated by white space.

Columns can be listed in any order. Columns from the required columns list
below must be included, and only recognised columns can be included in the
file.

All other lines in the file are rows defining a particular individual.
Each row must contain data for all columns defined in the header,
separated by whitespace.

#### Required columns

id
: A unique string identifier for this individual.

parent0
: The parent that will be assigned to ``parents[0]`` in the
  {ref}`tskit encoding<sec_pedigrees_encoding>`. This must
  be either a identifier defined in the ``id`` column, or
  the value ``.`` (representing missing data).

parent1
: The parent that will be assigned to ``parents[1]`` in the
  {ref}`tskit encoding<sec_pedigrees_encoding>`. The value
  requirements are identical to the ``parent0`` column.

time
: The time value to associate with the individual.

:::{note}
   The ``time`` column may become optional in later releases.
:::

#### Optional columns

is_sample
: The sample status of the individual. If equal to ``1``, the
  nodes associated with this individual will be marked as samples.
  If ``0`` they will not. No other values are supported.

population
: The population to associate with the individual. Values must
  correspond to
  {ref}`population identifiers<sec_demography_populations_identifiers>`
  in a supplied {class}`.Demography` object


(sec_pedigrees_file_format_basic_example)=

### Basic example

In this example we encode a trio, ``child``, ``mom`` and ``dad``:

```{code-cell}
txt = """\
# id parent0 parent1 time
child mom dad 0
mom    .   .  1
dad    .   .  1
"""
pedigree = msprime.parse_pedigree(io.StringIO(txt))
display(pedigree.individuals)
display(pedigree.nodes)
```

We have three individuals in our input pedigree, and there is therefore
three rows in the individual table. In the tskit encoding, individuals
are referred to by their integer ID (the corresponding row in the
individual table). Individuals are added to the table in the same order
they appear in the file, and therefore ``child`` corresponds to ID 0,
``mom`` ID 1 and ``dad`` ID 2. These ID values are then used in the
``parents`` column to describe the pedigree relationships.

Time information about individuals is stored in the node table,
where the time of an individual is associated with its two nodes.
For example, nodes 0 and 1 correspond to the two genomes of
individual 0 (``child``), and these are both at time 0.

The individual table also includes the original string ``file_id`` as
metadata to faciliate joining with existing data sources. We can
access these IDs as follows:

```{code-cell}
for id, ind in enumerate(pedigree.individuals):
    print(id, "->", ind.metadata["file_id"])
```

```{seealso}
See the tskit {ref}`metadata documentation<tskit:sec_metadata>` for more information
on how to use metadata.
```

### Specifying samples

In the previous example we did not specify which of our individuals were
the samples ("probands"). In this case, {func}`.parse_pedigree` assumes
that any individual at time 0 is a sample. This is show in the example
above in the node table, where we can see that the flags value for
both of ``child``'s nodes is 1 (corresponding to {data}`tskit.NODE_IS_SAMPLE`)

We can override this behaviour by providing an ``is_sample`` column:

```{code-cell}
import io
txt = """\
# id parent0 parent1 time is_sample
child1 mom1 dad1 0 1
mom1     .   .   1 0
dad1     .   .   1 0
child2 mom2 dad2 1 1
mom2     .   .   2 0
dad2     .   .   2 0
"""
pedigree = msprime.parse_pedigree(io.StringIO(txt))
display(pedigree.individuals)
display(pedigree.nodes)
```

Here we have two trios, where the child in the second trio is from the
same generation as the parents in the second. We use the ``is_sample``
column to specify that ``child2`` is a sample as well as ``child1``.


### Demography information

If the founder individuals of the pedigree belong to different populations
and we wish to simulate ancestry
{ref}`beyond the pedigree<sec_ancestry_models_fixed_pedigree_completing>`,
we must define a {ref}`demographic model<sec_demography>` and define
the populations that the individuals belong to.

:::{seealso}
See the {ref}`sec_ancestry_models_fixed_pedigree_demography` section for
a full discussion of how demography interacts with fixed pedigree simulations
and important caveats.
:::

To do this we must pass the {class}`.Demography` instance defining the
demographic model to {func}`.parse_pedigree`, and include
a ``population`` column in the file.

```{code-cell}
import io
txt = """\
# id parent0 parent1 time is_sample population
child1 mom1 dad1 0 1 A
mom1     .   .   1 0 A
dad1     .   .   1 0 A
child2 mom2 dad2 1 1 B
mom2     .   .   2 0 B
dad2     .   .   2 0 B
"""

demography = msprime.Demography()
demography.add_population(name="A", initial_size=10)
demography.add_population(name="B", initial_size=20)
demography.add_population(name="C", initial_size=100)
demography.add_population_split(time=10, derived=["A", "B"], ancestral="C");

pedigree = msprime.parse_pedigree(io.StringIO(txt), demography=demography)
display(pedigree.populations)
display(pedigree.individuals)
display(pedigree.nodes)
```

(sec_pedigrees_visualisation)=

## Visualising pedigrees

It is often useful to visualise a (small) pedigree. The following function
is used in this documentation, which may be a useful recipe for others:

```{code-cell}
def draw_pedigree(ped_ts):

    G = nx.DiGraph()
    for ind in ped_ts.individuals():
        time = ped_ts.node(ind.nodes[0]).time
        pop = ped_ts.node(ind.nodes[0]).population
        G.add_node(ind.id, time=time, population=pop)
        for p in ind.parents:
            if p != tskit.NULL:
                G.add_edge(ind.id, p)
    pos = nx.multipartite_layout(G, subset_key="time", align="horizontal")
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    node_colours = [colours[node_attr["population"]] for node_attr in G.nodes.values()]
    nx.draw_networkx(G, pos, with_labels=True, node_color=node_colours)
    plt.show()
```

See the {ref}`sec_ancestry_models_fixed_pedigree_demography` section
for an example of this function drawing pedigrees in a multi-population
model.
