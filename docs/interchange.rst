.. _sec_interchange:

#########################
Tree sequence interchange
#########################

The correlated genealogical trees that describe the shared ancestry of a set of
samples are stored concisely in ``msprime`` as a collection of
easy-to-understand tables. These are output by coalescent simulation in
``msprime`` or can be read in from another source. This page documents
the structure of the tables, and the different methods of interchanging
genealogical data to and from the msprime API. We begin by defining
the basic concepts that we need and the structure of the tables in the
`Data model`_ section. We then describe the tabular text formats that can
be used as simple interchange mechanism for small amounts of data in the
`Text file formats`_ section. The `Binary interchange`_ section then describes
the efficient Python API for table interchange using numpy arrays. Finally,
we describe the binary format used by msprime to efficiently
store tree sequences on disk in the `Tree sequence file format`_ section.


.. _sec_data_model:

**********
Data model
**********

To begin, here are definitions of some key ideas encountered later.

tree
    A "gene tree", i.e., the genealogical tree describing how a collection of
    genomes (usually at the tips of the tree) are related to each other at some
    chromosomal location. See :ref:`sec_nodes_or_individuals` for discussion
    of what a "genome" is.

tree sequence
    A "succinct tree sequence" (or tree sequence, for brevity) is an efficient
    encoding of a sequence of correlated trees, such as one encounters looking
    at the gene trees along a genome. A tree sequence efficiently captures the
    structure shared by adjacent trees, (essentially) storing only what differs
    between them.

node
    Each branching point in each tree is associated with a particular genome
    in a particular ancestor, called "nodes".  Since each node represents a
    specific genome it has a unique ``time``, thought of as its birth time,
    which determines the height of any branching points it is associated with.
    See :ref:`sec_nodes_or_individuals` for discussion of what a "node" is.

individual

    In certain situations we are interested in how nodes (representing
    individual homologous genomes) are grouped together into individuals
    (e.g., two nodes per diploid individual). For example, when we are working
    with polyploid samples it is useful to associate metadata with a specific
    individual rather than duplicate this information on the constituent nodes.
    See :ref:`sec_nodes_or_individuals` for more discussion on this point.

sample
    The focal nodes of a tree sequence, usually thought of as those that we
    have obtained data from. The specification of these affects various
    methods: (1) :meth:`TreeSequence.variants` and
    :meth:`TreeSequence.haplotypes` will output the genotypes of the samples,
    and :meth:`SparseTree.roots` only return roots ancestral to at least one
    sample. (See the :ref:`node table definitions <sec_node_table_definition>`
    for information on how the sample
    status a node is encoded in the ``flags`` column.)

edge
    The topology of a tree sequence is defined by a set of **edges**. Each
    edge is a tuple ``(left, right, parent, child)``, which records a
    parent-child relationship among a pair of nodes on the
    on the half-open interval of chromosome ``[left, right)``.

site
    Tree sequences can define the mutational state of nodes as well as their
    topological relationships. A **site** is thought of as some position along
    the genome at which variation occurs. Each site is associated with
    a unique position and ancestral state.

mutation
    A mutation records the change of state at a particular site 'above'
    a particular node (more precisely, along the branch between the node
    in question and its parent). Each mutation is associated with a specific
    site (which defines the position along the genome), a node (which defines
    where it occurs within the tree at this position), and a derived state
    (which defines the mutational state inherited by all nodes in the subtree
    rooted at the focal node). In more complex situations in which we have
    back or recurrent mutations, a mutation must also specify its 'parent'
    mutation.

migration
    An event at which a parent and child node were born in different populations.

population
    A grouping of nodes, e.g., by sampling location.

provenance
    An entry recording the origin and history of the data encoded in a tree sequence.

ID
    In the set of interconnected tables that we define here, we refer
    throughout to the IDs of particular entities. The ID of an
    entity (e.g., a node) is defined by the position of the corresponding
    row in the table. These positions are zero indexed. For example, if we
    refer to node with ID zero, this corresponds to the node defined by the
    first row in the node table.

Sequence length
    This value defines the coordinate space in which the edges and site positions
    are defined. This is most often assumed to be equal to the largest
    ``right`` coordinate in the edge table, but there are situations in which
    we might wish to specify the sequence length explicitly.


A tree sequence can be stored in a collection of eight tables:
:ref:`Node <sec_node_table_definition>`,
:ref:`Edge <sec_edge_table_definition>`,
:ref:`Individual <sec_individual_table_definition>`,
:ref:`Site <sec_site_table_definition>`,
:ref:`Mutation <sec_mutation_table_definition>`,
:ref:`Migration <sec_migration_table_definition>`,
:ref:`Population <sec_population_table_definition>`, and
:ref:`Provenance <sec_provenance_table_definition>`.
The Node and Edge tables store the genealogical
relationships that define the trees, and the Individual table
describes how multiple genomes are grouped within individuals;
the Site and Mutation tables describe where mutations fall
on the trees; the Migration table describes how lineages move across space;
and the Provenance table contains information on where the data came from.
Only Node and Edge tables are necessary to encode the genealogical trees;
Sites and Mutations are optional but necessary to encode polymorphism
(sequence) data; the remainder are optional.
In the following sections we define these components of a tree sequence in
more detail.


.. _sec_nodes_or_individuals:

Nodes, Genomes, or Individuals?
===============================

The natural unit of biological analysis is (usually) the *individual*. However,
many organisms we study are diploid, and so each individual contains *two*
homologous copies of the entire genome, separately inherited from the two
parental individuals. Since each monoploid copy of the genome is inherited separately,
each diploid individual lies at the end of two distinct lineages, and so will
be represented by *two* places in any given genealogical tree. This makes it
difficult to precisely discuss tree sequences for diploids, as we have no
simple way to refer to the bundle of chromosomes that make up the "copy of the
genome inherited from one particular parent". For this reason, in this
documentation we use the non-descriptive term "node" to refer to this concept
-- and so, a diploid individual is composed of two nodes -- although we use the
term "genome" at times, for concreteness.

Several properties naturally associated with individuals are in fact assigned
to nodes in what follows: birth time and population. This is for two reasons:
First, since coalescent simulations naturally lack a notion of polyploidy, earlier
versions of ``msprime`` lacked the notion of an individual. Second, ancestral
nodes are not naturally grouped together into individuals -- we know they must have
existed, but have no way of inferring this grouping, so in fact many nodes in
an empirically-derived tree sequence will not be associated with individuals,
even though their birth times might be inferred.


.. _sec_table_definitions:

*****************
Table definitions
*****************

.. _sec_table_types_definitions:

Table types
===========

.. _sec_node_table_definition:

Node Table
----------

A **node** defines a monoploid set of chromosomes (a "genome") of a specific
individual that was born at some time in the past: the set of
chromosomes inherited from a particular one of the individual's parents.
(See :ref:`sec_nodes_or_individuals` for more discussion.)
Every vertex in the marginal trees of a tree sequence corresponds
to exactly one node, and a node may be present in many trees. The
node table contains five columns, of which ``flags`` and ``time`` are
mandatory:

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
flags               uint32              Bitwise flags.
time                double              Birth time of node
population          int32               Birth population of node.
individual          int32               The individual the node belongs to.
metadata            binary              Node :ref:`sec_metadata_definition`
================    ==============      ===========

The ``time`` column records the birth time of the individual in question,
and is a floating point value. Similarly,
the ``population`` column records the ID of the population where this
individual was born. If not provided, ``population`` defaults to the
null ID (-1). Otherwise, the population ID must refer to a row in the
:ref:`sec_population_table_definition`.
The ``individual`` column records the ID of the
:ref:`Individual <sec_individual_table_definition>`
individual that this node belongs to. If specified, the ID must refer
to a valid individual. If not provided, ``individual``
defaults to the null ID (-1).

The ``flags`` column stores information about a particular node, and
is composed of 32 bitwise boolean values. Currently, the only flag defined
is ``IS_SAMPLE = 1``, which defines the *sample* status of nodes. Marking
a particular node as a "sample" means, for example, that the mutational state
of the node will be included in the genotypes produced by
:meth:`.TreeSequence.variants`.

Bits 0-15 (inclusive) of the ``flags`` column are reserved for internal use by
``tskit`` and should not be used by applications for anything other
than the purposes documented here. Bits 16-31 (inclusive) are free for applications
to use for any purpose and will not be altered or interpreteted by
``tskit``.

See the :ref:`sec_node_requirements` section for details on the properties
required for a valid set of nodes.

For convenience, the :ref:`text format <sec_text_file_format>` for nodes
decomposes the ``flags`` value into its separate values. Thus, in the
text format we have a column for ``is_sample``, which corresponds to the
``flags`` column in the underlying table. As more flags values are
defined, these will be added to the text file format.

The ``metadata`` column provides a location for client code to store
information about each node. See the :ref:`sec_metadata_definition` section for
more details on how metadata columns should be used.

.. note::
    The distinction between ``flags`` and ``metadata`` is that flags
    holds information about a node that the library understands, whereas
    metadata holds information about a node that the library *does not*
    understand. Metadata is for storing auxiliarly information that is
    not necessary for the core tree sequence algorithms.

.. _sec_individual_table_definition:


Individual Table
----------------

An **individual** defines how nodes (which can be seen
as representing single chromosomes) group together in a polyploid individual.
The individual table contains three columns, of which only ``flags`` is mandatory.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
flags               uint32              Bitwise flags.
location            double              Location in arbitrary dimensions
metadata            binary              Individual :ref:`sec_metadata_definition`
================    ==============      ===========

See the :ref:`sec_individual_requirements` section for details on the properties
required for a valid set of individuals.

The ``flags`` column stores information about a particular individual, and
is composed of 32 bitwise boolean values. Currently, no flags are
defined.

Bits 0-15 (inclusive) of the ``flags`` column are reserved for internal use by
``tskit`` and should not be used by applications for anything other
than the purposes documented here. Bits 16-31 (inclusive) are free for applications
to use for any purpose and will not be altered or interpreteted by
``tskit``.

The ``location`` column stores the location of an individual in arbitrary
dimensions. This column is :ref:`ragged <sec_encoding_ragged_columns>`, and
so different individuals can have locations with different dimensions (i.e.,
one individual may have location ``[]`` and another ``[0, 1, 0]``. This could
therefore be used to store other quantities (e.g., phenotype).

The ``metadata`` column provides a location for client code to store
information about each individual. See the :ref:`sec_metadata_definition` section for
more details on how metadata columns should be used.

.. note::
    The distinction between ``flags`` and ``metadata`` is that flags
    holds information about a individual that the library understands, whereas
    metadata holds information about a individual that the library *does not*
    understand. Metadata is for storing auxiliarly information that is
    not necessary for the core tree sequence algorithms.


.. _sec_edge_table_definition:

Edge Table
----------

An **edge** defines a parent-child relationship between a pair of nodes
over a specific sequence interval. The edge table contains four columns,
all of which are mandatory:

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
left                double              Left coordinate of the edge (inclusive).
right               double              Right coordinate of the edge (exclusive).
parent              int32               Parent node ID.
child               int32               Child node ID.
================    ==============      ===========

Each row in an edge table describes a half-open genomic interval ``[left, right)``
over which the ``child`` inherited from the given ``parent``.
The ``left`` and ``right`` columns are defined using double precision
floating point values. The ``parent`` and ``child``
columns specify integer IDs in the associated :ref:`sec_node_table_definition`.

See the :ref:`sec_edge_requirements` section for details on the properties
required for a valid set of edges.

.. _sec_site_table_definition:

Site Table
----------

A **site** defines a particular location along the genome in which
we are interested in observing the allelic state. The site table
contains three columns, of which ``position`` and ``ancestral_state``
are mandatory.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
position            double              Position of site in genome coordinates.
ancestral_state     text                The state at the root of the tree.
metadata            binary              Site :ref:`sec_metadata_definition`.
================    ==============      ===========

The ``position`` column is a floating point value defining the location
of the site in question along the genome.

The ``ancestral_state`` column specifies the allelic state at the root
of the tree, thus defining the state that nodes inherit if no mutations
intervene. The column stores text character data of arbitrary length.

The ``metadata`` column provides a location for client code to store
information about each site. See the :ref:`sec_metadata_definition` section for
more details on how metadata columns should be used.

See the :ref:`sec_site_requirements` section for details on the properties
required for a valid set of sites.

.. _sec_mutation_table_definition:

Mutation Table
--------------

A **mutation** defines a change of allelic state on a tree at a particular site.
The mutation table contains five columns, of which ``site``, ``node`` and
``derived_state`` are mandatory.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
site                int32               The ID of the site the mutation occurs at.
node                int32               The node this mutation occurs at.
parent              int32               The ID of the parent mutation.
derived_state       char                The allelic state resulting from the mutation.
metadata            binary              Mutation :ref:`sec_metadata_definition`.
================    ==============      ===========

The ``site`` column is an integer value defining the ID of the
:ref:`site <sec_site_table_definition>` at which this mutation occurred.

The ``node`` column is an integer value defining the ID of the
first :ref:`node <sec_node_table_definition>` in the tree below this mutation.

The ``derived_state`` column specifies the allelic state resulting from the mutation,
thus defining the state that the ``node`` and any descendant nodes in the
subtree inherit unless further mutations occur. The column stores text
character data of arbitrary length.

The ``parent`` column is an integer value defining the ID of the mutation whose
allelic state this mutation replaced. If there is no mutation at the
site in question on the path back to root, then this field is set to the
null ID (-1). (The ``parent`` column is only required in situations
where there are multiple mutations at a given site. For
"infinite sites" mutations, it can be ignored.)

The ``metadata`` column provides a location for client code to store
information about each site. See the :ref:`sec_metadata_definition` section for
more details on how metadata columns should be used.

See the :ref:`sec_mutation_requirements` section for details on the properties
required for a valid set of mutations.

.. _sec_migration_table_definition:

Migration Table
---------------

In simulations, trees can be thought of as spread across space, and it is
helpful for inferring demographic history to record this history.
Migrations are performed by individual ancestors, but most likely not by an
individual whose genome is tracked as a ``node`` (as in a discrete-deme model they are
unlikely to be both a migrant and a most recent common ancestor).  So,
``msprime`` records when a segment of ancestry has moved between
populations. This table is not required, even if different nodes come from
different populations.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
left                double              Left coordinate of the migrating segment (inclusive).
right               double              Right coordinate of the migrating segment (exclusive).
node                int32               Node ID.
source              int32               Source population ID.
dest                int32               Destination population ID.
time                double              Time of migration event.
================    ==============      ===========


The ``left`` and ``right`` columns are floating point values defining the
half-open segment of genome affected. The ``source`` and ``dest`` columns
record the IDs of the respective populations. The ``node`` column records the
ID of the node that was associated with the ancestry segment in question
at the time of the migration event. The ``time`` column is holds floating
point values recording the time of the event.

See the :ref:`sec_migration_requirements` section for details on the properties
required for a valid set of mutations.

.. _sec_population_table_definition:

Population Table
----------------

A **population** defines a grouping of individuals that a node can
be said to belong to.

The population table contains one column, ``metadata``.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
metadata            binary              Population :ref:`sec_metadata_definition`.
================    ==============      ===========


The ``metadata`` column provides a location for client code to store
information about each population. See the :ref:`sec_metadata_definition` section for
more details on how metadata columns should be used.

See the :ref:`sec_population_requirements` section for details on the properties
required for a valid set of populations.


.. _sec_provenance_table_definition:

Provenance Table
----------------

.. todo::
    Document the provenance table.


================    ==============      ===========
Column              Type                Description
================    ==============      ===========
timestamp           char                Timestamp in `ISO-8601 <https://en.wikipedia.org/wiki/ISO_8601>`_ format.
record              char                Provenance record.
================    ==============      ===========



.. _sec_metadata_definition:

Metadata
========

Users of the tables API sometimes need to store auxiliary information for
the various entities defined here. For example, in a forwards-time simulation,
the simulation engine may wish to store the time at which a particular mutation
arose or some other pertinent information. If we are representing real data,
we may wish to store information derived from a VCF INFO field, or associate
information relating to samples or populations. The columns defined in tables
here are deliberately minimal: we define columns only for information which
the library itself can use. All other information is considered to be
**metadata**, and is stored in the ``metadata`` columns of the various
tables.

Arbitrary binary data can be stored in ``metadata`` columns, and the
``msprime`` library makes no attempt to interpret this information. How the
information held in this field is encoded is entirely the choice of client code.

To ensure that metadata can be safely interchanged using the :ref:`sec_text_file_format`,
each row is `base 64 encoded <https://en.wikipedia.org/wiki/Base64>`_. Thus,
binary information can be safely printed and exchanged, but may not be
human readable.

.. todo::
    We plan on providing more sophisticated tools for working with metadata
    in future, including the auto decoding metadata via pluggable
    functions and the ability to store metadata schemas so that metadata
    is self-describing.


.. _sec_valid_tree_sequence_requirements:

Valid tree sequence requirements
================================

Arbitrary data can be stored in tables using the classes in the
:ref:`sec_tables_api`. However, only a :class:`.TableCollection`
that fulfils a set of requirements represents
a valid :class:`.TreeSequence` object which can be obtained
using the :meth:`.TableCollection.tree_sequence` method. In this
section we list these requirements, and explain their rationale.
Violations of most of these requirements are detected when the
user attempts to load a tree sequence via :func:`.load` or
:meth:`.TableCollection.tree_sequence`, raising an informative
error message. Some more complex requirements may not be detectable at load-time,
and errors may not occur until certain operations are attempted.
These are documented below.
We also provide tools that can transform a collection of tables into a valid
collection of tables, so long as they are logically consistent,
as described in :ref:`sec_table_transformations`.

.. _sec_individual_requirements:

Individual requirements
-----------------------

Individuals are a basic type in a tree sequence and are not defined with
respect to any other tables. Therefore, there are no requirements on
individuals.

There are no requirements regarding the ordering of individuals.
Sorting a set of tables using :meth:`.TableCollection.sort` has
no effect on the individuals.

.. _sec_node_requirements:

Node requirements
-----------------

Given a valid set of individuals and populations, the requirements for
each node are:

- ``population`` must either be null (-1) or refer to a valid population ID;
- ``individual`` must either be null (-1) or refer to a valid individual ID.

An ID refers to a zero-indexed row number in the relevant table,
and so is "valid" if is between 0 and one less than the number of rows in the relevant table.

There are no requirements regarding the ordering of nodes with respect to time.

For simplicity and algorithmic efficiency, all nodes referring to the same
(non-null) individual must be contiguous.

Sorting a set of tables using :meth:`.TableCollection.sort`
has no effect on nodes.

.. _sec_edge_requirements:

Edge requirements
-----------------

Given a valid set of nodes and a sequence length :math:`L`, the simple
requirements for each edge are:

- We must have :math:`0 \leq` ``left`` :math:`<` ``right`` :math:`\leq L`;
- ``parent`` and ``child`` must be valid node IDs;
- ``time[parent]`` > ``time[child]``;
- edges must be unique (i.e., no duplicate edges are allowed).

The first requirement simply ensures that the interval makes sense. The
third requirement ensures that we cannot have loops, since time is
always increasing as we ascend the tree.

To ensure a valid tree sequence there is one further requirement:

- The set of intervals on which each node is a child must be disjoint.

This guarantees that we cannot have contradictory edges (i.e.,
where a node ``a`` is a child of both ``b`` and ``c``), and ensures that
at each point on the sequence we have a well-formed forest of trees.
Because this is a more complex semantic requirement, it is **not** detected
at load time. This error is detected during tree traversal, via, e.g.,
the :meth:`.TreeSequence.trees` iterator.

In the interest of algorithmic efficiency, edges must have the following
sortedness properties:

- All edges for a given parent must be contiguous;
- Edges must be listed in nondecreasing order of ``parent`` time;
- Within the edges for a given ``parent``, edges must be sorted
  first by ``child`` ID and then by ``left`` coordinate.

Violations of these requirements are detected at load time.
The :meth:`.TableCollection.sort` method will ensure that these sortedness
properties are fulfilled.

.. _sec_site_requirements:

Site requirements
-----------------

Given a valid set of nodes and a sequence length :math:`L`, the simple
requirements for a valid set of sites are:

- We must have :math:`0 \leq` ``position`` :math:`< L`;
- ``position`` values must be unique.

For simplicity and algorithmic efficiency, sites must also:

- Be sorted in increasing order of ``position``.

Violations of these requirements are detected at load time.
The :meth:`.TableCollection.sort` method ensures that sites are sorted
according to these criteria.

.. _sec_mutation_requirements:

Mutation requirements
---------------------

Given a valid set of nodes, edges and sites, the
requirements for a valid set of mutations are:

- ``site`` must refer to a valid site ID;
- ``node`` must refer to a valid node ID;
- ``parent`` must either be the null ID (-1) or a valid mutation ID within the
  current table

Furthermore,

- If another mutation occurs on the tree above the mutation in
  question, its ID must be listed as the ``parent``.

For simplicity and algorithmic efficiency, mutations must also:

- be sorted by site ID;
- when there are multiple mutations per site, parent mutations must occur
  **before** their children (i.e. if a mutation with ID :math:`x` has
  ``parent`` with ID :math:`y`, then we must have :math:`y < x`).

Violations of these sorting requirements are detected at load time.
The :meth:`.TableCollection.sort` method ensures that mutations are sorted
according site ID, but does not at present enforce that mutations occur
after their parent mutations.

Mutations also have the requirement that they must result in a
change of state. For example, if we have a site with ancestral state
of "A" and a single mutation with derived state "A", then this
mutation does not result in any change of state. This error is
raised at run-time when we reconstruct sample genotypes, for example
in the :meth:`.TreeSequence.variants` iterator.

.. _sec_migration_requirements:

Migration requirements
----------------------

Given a valid set of nodes and edges, the requirements for a value set of
migrations are:

- ``left`` and ``right`` must lie within the tree sequence coordinate space (i.e.,
  from 0 to ``sequence_length``).
- ``time`` must be strictly between the time of its ``node`` and the time of any
  ancestral node from which that node inherits on the segment ``[left, right)``.
- The ``population`` of any such ancestor matching ``source``, if another
  ``migration`` does not intervene.

To enable efficient processing, migrations must also be:

- Sorted by nondecreasing ``time`` value.

Note in particular that there is no requirement that adjacent migration records
should be "squashed". That is, we can have two records ``m1`` and ``m2``
such that ``m1.right`` = ``m2.left`` and with the ``node``, ``source``,
``dest`` and ``time`` fields equal. This is because such records will usually
represent two independent ancestral segments migrating at the same time, and
as such squashing them into a single record would result in a loss of information.


.. _sec_population_requirements:

Population requirements
-----------------------

There are no requirements on a population table.

.. _sec_provenance_requirements:

Provenance requirements
-----------------------

The `timestamp` column of a provenance table should be in
`ISO-8601 <https://en.wikipedia.org/wiki/ISO_8601>`_ format.

The `record` should be valid JSON with structure defined in the Provenance
Schema section (TODO).


.. _sec_table_transformations:

Table transformation methods
============================

The following methods operate *in place* on a :class:`TableCollection`,
transforming them while preserving information.
In some applications, tables may most naturally be produced in a way that is
logically consistent, but not meeting all the requirements for validity that
are established for algorithmic and efficiency reasons.
These methods (while having other uses), can be used to make such a set of
tables valid, and thus ready to be loaded into a tree sequence.

This section is best skipped unless you are writing a program that records
tables directly.

Simplification
--------------

Simplification of a tree sequence is in fact a transformation method applied
to the underlying tables: the method :meth:`TreeSequence.simplify` calls
:meth:`TableCollection.simplify` on the tables, and loads a new tree sequence.
The main purpose of this method is to remove redundant information,
only retaining the minimal tree sequence necessary to describe the genealogical
history of the ``samples`` provided.

Furthermore, ``simplify`` is guaranteed to:

- preserve relative ordering of any rows in the Site and Mutation tables
  that are not discarded.

The :meth:`TableCollection.simplify` method can be applied to a collection of
tables that does not have the ``mutations.parent`` entries filled in, as long
as all other validity requirements are satisfied.

Sorting
-------

The validity requirements for a set of tables to be loaded into a tree sequence
listed in :ref:`sec_table_definitions` are of two sorts: logical consistency,
and sortedness. The :meth:`TableCollection.sort` method can be used to make
completely valid a set of tables that satisfies all requirements other than
sortedness.

This method can also be used on slightly more general collections of tables:
it is not required that ``site`` positions be unique in the table collection to
be sorted. The method has two additional properties:

- it preserves relative ordering between sites at the same position, and
- it preserves relative ordering between mutations at the same site.

:meth:`TableCollection.sort` does not check the validity of the `parent`
property of the mutation table. However, because the method preserves mutation
order among mutations at the same site, if mutations are already sorted so that
each mutation comes after its parent (e.g., if they are ordered by time of
appearance), then this property is preserved, even if the `parent` properties
are not specified.


Removing duplicate sites
------------------------

The :meth:`TableCollection.deduplicate_sites` method can be used to save a tree
sequence recording method the bother of checking to see if a given site already
exists in the site table. If there is more than one site with the same
position, all but the first is removed, and all mutations referring to the
removed sites are edited to refer to the first (and remaining) site. Order is
preserved.


Computing mutation parents
--------------------------

If each edge had at most only a single mutation, then the ``parent`` property
of the mutation table would be easily inferred from the tree at that mutation's
site. If mutations are entered into the mutation table ordered by time of
appearance, then this sortedness allows us to infer the parent of each mutation
even for mutations occurring on the same branch. The
:meth:`TableCollection.compute_mutation_parents` method will take advantage
of this fact to compute the ``parent`` column of a mutation table, if all
other information is valid.


Recording tables in forwards time
---------------------------------

The above methods enable the following scheme for recording site and mutation
tables during a forwards-time simulation. Whenever a new mutation is
encountered:

1. Add a new ``site`` to the site table at this position.
2. Add a new ``mutation`` to the mutation table at the newly created site.

This is lazy and wrong, because:

a. There might have already been sites in the site table with the same position,
b. and/or a mutation (at the same position) that this mutation should record as
   its ``parent``.

But, it's all OK because here's what we do:

1. Add rows to the mutation and site tables as described above.
2. Periodically, ``sort``, ``deduplicate_sites``,  and ``simplify``, then
   return to (1.), except that
3. Sometimes, to output the tables, ``sort``, ``compute_mutation_parents``,
    (optionally ``simplify``), and dump these out to a file.

*Note:* as things are going along we do *not* have to
``compute_mutation_parents``, which is nice, because this is a nontrivial step
that requires construction all the trees along the sequence. Computing mutation
parents only has to happen before the final (output) step.

This is OK as long as the forwards-time simulation outputs things in order by when
they occur, because these operations have the following properties:

1. Mutations appear in the mutation table ordered by time of appearance, so
   that a mutation will always appear after the one that it replaced (i.e.,
   its parent).
2. Furthermore, if mutation B appears after mutation A, but at the same site,
   then mutation B's site will appear after mutation A's site in the site
   table.
3. ``sort`` sorts sites by position, and then by ID, so that the relative
   ordering of sites at the same position is maintained, thus preserving
   property (2).
4. ``sort`` sorts mutations by site, and then by ID, thus preserving property
   (1); if the mutations are at separate sites (but the same position), this
   fact is thanks to property (2).
5. ``simplify`` also preserves ordering of any rows in the site and mutation
   tables that do not get discarded.
6. ``deduplicate_sites`` goes through and collapses all sites at the same
   position to only one site, maintaining order otherwise.
7. ``compute_mutation_parents`` fills in the ``parent`` information by using
    property (1).


.. _sec_text_file_format:

*****************
Text file formats
*****************

The tree sequence text file format is based on a simple whitespace
delimited approach. Each table corresponds to a single file, and is
composed of a number of whitespace delimited columns. The first
line of each file must be a **header** giving the names of each column.
Subsequent rows must contain data for each of these columns, following
the usual conventions. Each table has a set of mandatory and optional columns which are
described below. The columns can be provided in any order, and extra columns
can be included in the file. Note, in particular, that this means that
an ``id`` column may be present in any of these files, but it will be
ignored (IDs are always determined by the position of the row in a table).

We present the text format below using the following very simple tree
sequence, with four nodes, two trees, and three mutations at two sites,
both on the first tree::

    time ago
    --------
      3            3
                ┏━━┻━━┓
                ╋     ╋         2
                ┃     ╋      ┏━━┻━━┓
      0         0     1      0     1

    position  0           7          10

A deletion from AT to A has occurred at position 2 on the branch leading to
node 0, and two mutations have occurred at position 4 on the branch leading to
node 1, first from A to T, then a back mutation to A. The genotypes of our two
samples, nodes 0 and 1, are therefore AA and ATA.


.. _sec_individual_text_format:

Individual text format
======================

The individual text format must contain a ``flags`` column.
Optionally, there may also be a ``location`` and
``metadata`` columns. See the :ref:`individual table definitions
<sec_individual_table_definition>` for details on these columns.

Note that there are currently no globally defined ``flags``, but the column
is still required; a value of ``0`` means that there are no flags set.

The ``location`` column should be a sequence of comma-separated numeric
values. They do not all have to be the same length.

An example individual table::

    flags   location
    0           0.5,1.2
    0           1.0,3.4
    0
    0           1.2
    0           3.5,6.3
    0           0.5,0.5
    0           0.5
    0           0.7,0.6,0.0
    0           0.5,0.0


.. _sec_node_text_format:

Node text format
================

The node text format must contain the columns ``is_sample`` and
``time``. Optionally, there may also be ``population``, ``individual``, and
``metadata`` columns. See the :ref:`node table definitions
<sec_node_table_definition>` for details on these columns.

Note that we do not have a ``flags`` column in the text file format, but
instead use ``is_sample`` (which may be 0 or 1). Currently, ``IS_SAMPLE`` is
the only flag value defined for nodes, and as more flags are defined we will
allow for extra columns in the text format.

An example node table::

    is_sample   individual   time
    1           0            0.0
    1           0            0.0
    0           -1           1.0
    0           -1           3.0

.. _sec_edge_text_format:

Edge text format
================

The edge text format must contain the columns ``left``,
``right``, ``parent`` and ``child``.
See the :ref:`edge table definitions <sec_edge_table_definition>`
for details on these columns.

An example edge table::

    left   right   parent  child
    0.0    7.0     2       0
    0.0    7.0     2       1
    7.0    10.0    3       0
    7.0    10.0    3       1


.. _sec_site_text_format:

Site text format
================

The site text format must contain the columns ``position`` and
``ancestral_state``. The ``metadata`` column may also be optionally
present. See the
:ref:`site table definitions <sec_site_table_definition>`
for details on these columns.

sites::

    position      ancestral_state
    2.0           AT
    4.0           A

.. _sec_mutation_text_format:

Mutation text format
====================

The mutation text format must contain the columns ``site``,
``node`` and ``derived_state``. The ``parent`` and ``metadata`` columns
may also be optionally present (but ``parent`` must be specified if
more than one mutation occurs at the same site). See the
:ref:`mutation table definitions <sec_site_table_definition>`
for details on these columns.

mutations::

    site   node    derived_state    parent
    0      0       A                -1
    1      0       T                -1
    1      1       A                1



Population text format
======================

Population tables only have a ``metadata`` column, so the text format for
a population table requires there to be a ``metadata`` column. See the
:ref:`population table definitions <sec_population_table_definition>` for
details.

An example population table::

    id   metadata
    0    cG9wMQ==
    1    cG9wMg==

The ``metadata`` contains base64-encoded data (in this case, the strings
``pop1`` and ``pop1``).


.. _sec_binary_interchange:

******************
Binary interchange
******************

In this section we describe the high-level details of the API for interchanging
table data via numpy arrays. Please see the :ref:`sec_tables_api` for detailed
description of the functions and methods.

The tables API is based on **columnar** storage of the data. In memory, each
table is organised as a number of blocks of contiguous storage, one for
each column. There are many advantages to this approach, but the key
property for us is that allows for very efficient transfer of data
in and out of tables. Rather than inserting data into tables row-by-row
(which can be done using the ``add_row`` methods), it is much more
efficient to add many rows at the same time by providing pointers to blocks of
contiguous memory. By taking
this approach, we can work with tables containing gigabytes of data very
efficiently.

We use the `numpy Array API <https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.html>`_
to allow us to define and work with numeric arrays of the required types.
Node IDs, for example, are defined using 32 bit integers. Thus, the
``parent`` column of an :ref:`sec_edge_table_definition`'s with ``n`` rows
is a block ``4n`` bytes.

This approach is very straightforward for columns in which each row contains
a fixed number of values. However, dealing with columns containing a
**variable** number of values is more problematic.

.. _sec_encoding_ragged_columns:

Encoding ragged columns
=======================

A **ragged** column is a column in which the rows are not of a fixed length.
For example, :ref:`sec_metadata_definition` columns contain binary of data of arbitrary
length. To encode such columns in the tables API, we store **two** columns:
one contains the flattened array of data and another stores the **offsets**
of each row into this flattened array. Consider an example::

    >>> s = msprime.SiteTable()
    >>> s.add_row(0, "A")
    >>> s.add_row(0, "")
    >>> s.add_row(0, "TTT")
    >>> s.add_row(0, "G")
    >>> print(s)
    id      position        ancestral_state metadata
    0       0.00000000      A
    1       0.00000000
    2       0.00000000      TTT
    3       0.00000000      G
    >>> s.ancestral_state
    array([65, 84, 84, 84, 71], dtype=int8)
    >>> s.ancestral_state.tobytes()
    b'ATTTG'
    >>> s.ancestral_state_offset
    array([0, 1, 1, 4, 5], dtype=uint32)
    >>> s.ancestral_state[s.ancestral_state_offset[2]: s.ancestral_state_offset[3]].tobytes()
    b'TTT'

In this example we create a :ref:`sec_site_table_definition` with four rows,
and then print out this table. We can see that the second row has the
empty string as its ``ancestral_state``, and the third row's
``ancestral_state`` is ``TTT``. When we print out the tables ``ancestral_state``
column, we see that its a numpy array of length 5: this is the
flattened array of `ASCII encoded <https://en.wikipedia.org/wiki/ASCII>`_
values for these rows. When we decode these bytes using the
numpy ``tobytes`` method, we get the string 'ATTTG'. This flattened array
can now be transferred efficiently in memory like any other column. We
then use the ``ancestral_state_offset`` column to allow us find the
individual rows. For a row ``j``::

    ancestral_state[ancestral_state_offset[j]: ancestral_state_offset[j + 1]]

gives us the array of bytes for the ancestral state in that row.

For a table with ``n`` rows, any offset column must have ``n + 1``
values, the first of which is always ``0``. The values in this column must be
nondecreasing, and cannot exceed the length of the ragged column in question.

.. _sec_tree_sequence_file_format:

**************************
Tree sequence file format
**************************

To make tree sequence data as efficient and easy as possible to use, we store the
data on file in a columnar, binary format. The format is based on the
`kastore <https://pypi.org/project/kastore/>`_ package, which is a simple
key-value store for numerical data. There is a one-to-one correspondence
between the tables described above and the arrays stored in these files.

By convention, these files are given the ``.trees`` suffix (although this
is not enforced in any way), and we will sometimes refer to them as ".trees"
files. We also refer to them as "tree sequence files".

.. todo::
    Link to the documentation for kastore, and describe the arrays that are
    stored as well as the top-level metadata.

.. The root group contains two attributes, ``format_version`` and ``sequence_length``.
.. The ``format_version`` is a pair ``(major, minor)`` describing the file format version.
.. This document describes version 10.0. The ``sequence_length`` attribute defines the
.. coordinate space over which edges and sites are defined. This must be present
.. and be greater than or equal to the largest coordinate present.

.. ================    ==============      ======      ===========
.. Path                Type                Dim         Description
.. ================    ==============      ======      ===========
.. /format_version     H5T_STD_U32LE       2           The (major, minor) file format version.
.. /sequence_length    H5T_IEEE_F64LE      1           The maximum value of a sequence coordinate.
.. ================    ==============      ======      ===========

.. Nodes group
.. ===========

.. The ``/nodes`` group stores the :ref:`sec_node_table_definition`.

.. =======================     ==============
.. Path                        Type
.. =======================     ==============
.. /nodes/flags                H5T_STD_U32LE
.. /nodes/population           H5T_STD_I32LE
.. /nodes/time                 H5T_IEEE_F64LE
.. /nodes/metadata             H5T_STD_I8LE
.. /nodes/metadata_offset      H5T_STD_U32LE
.. =======================     ==============

.. Edges group
.. ===========

.. The ``/edges`` group stores the :ref:`sec_edge_table_definition`.

.. ===================       ==============
.. Path                      Type
.. ===================       ==============
.. /edges/left               H5T_IEEE_F64LE
.. /edges/right              H5T_IEEE_F64LE
.. /edges/parent             H5T_STD_I32LE
.. /edges/child              H5T_STD_I32LE
.. ===================       ==============

.. Indexes group
.. -------------

.. The ``/edges/indexes`` group records information required to efficiently
.. reconstruct the individual trees from the tree sequence. The
.. ``insertion_order`` dataset contains the order in which records must be applied
.. and the ``removal_order`` dataset the order in which records must be
.. removed for a left-to-right traversal of the trees.

.. ==============================     ==============
.. Path                               Type
.. ==============================     ==============
.. /edges/indexes/insertion_order     H5T_STD_I32LE
.. /edges/indexes/removal_order       H5T_STD_I32LE
.. ==============================     ==============

.. Sites group
.. ===========

.. The sites group stores the :ref:`sec_site_table_definition`.

.. =============================   ==============
.. Path                            Type
.. =============================   ==============
.. /sites/position                 H5T_IEEE_F64LE
.. /sites/ancestral_state          H5T_STD_I8LE
.. /sites/ancestral_state_offset   H5T_STD_U32LE
.. /sites/metadata                 H5T_STD_I8LE
.. /sites/metadata_offset          H5T_STD_U32LE
.. =============================   ==============

.. Mutations group
.. ===============

.. The mutations group stores the :ref:`sec_mutation_table_definition`.

.. ===============================  ==============
.. Path                             Type
.. ===============================  ==============
.. /mutations/site                  H5T_STD_I32LE
.. /mutations/node                  H5T_STD_I32LE
.. /mutations/parent                H5T_STD_I32LE
.. /mutations/derived_state         H5T_STD_I8LE
.. /mutations/derived_state_offset  H5T_STD_U32LE
.. /mutations/metadata              H5T_STD_I8LE
.. /mutations/metadata_offset       H5T_STD_U32LE
.. ===============================  ==============

.. Migrations group
.. ================

.. The ``/migrations`` group stores the :ref:`sec_migration_table_definition`.

.. ===================       ==============
.. Path                      Type
.. ===================       ==============
.. /migrations/left          H5T_IEEE_F64LE
.. /migrations/right         H5T_IEEE_F64LE
.. /migrations/node          H5T_STD_I32LE
.. /migrations/source        H5T_STD_I32LE
.. /migrations/dest          H5T_STD_I32LE
.. /migrations/time          H5T_IEEE_F64LE
.. ===================       ==============

.. Provenances group
.. =================

.. The provenances group stores the :ref:`sec_provenance_table_definition`.

.. ===============================  ==============
.. Path                             Type
.. ===============================  ==============
.. /provenances/timestamp           H5T_STD_I8LE
.. /provenances/timestamp_offset    H5T_STD_U32LE
.. /provenances/record              H5T_STD_I8LE
.. /provenances/record_offset       H5T_STD_U32LE
.. ===============================  ==============


Legacy Versions
===============

Tree sequence files written by older versions of msprime are not readable by
newer versions of msprime. For major releases of msprime, :ref:`sec_msp_upgrade`
will convert older tree sequence files to the latest version.

File formats from version 11 onwards are based on
`kastore <https://pypi.org/project/kastore/>`_;
previous to this, the file format was based on HDF5.

However many changes to the tree sequence format are not part of major
releases. The table below gives these versions.

.. to obtain hashes where versions were changed:
        git log --oneline -L40,41:lib/msprime.h
   then on each hash, to obtain the parent where a merge occured:
        git log --merges --pretty=format:"%h" fc17dbd | head -n 1
   in some cases this didn't work so required hand manipulation. checks were
   done (after checkign out and rebuilding) with:
        python msp_dev.py simulate 10 tmp.trees && h5dump tmp.trees | head
   For versions 11 and onwards, use kastore to get the version:
        kastore dump format/version tmp.trees

=======    =================
Version    Commit Short Hash
=======    =================
11.0       5646cd3
10.0       e4396a7
9.0        e504abd
8.0        299ddc9
7.0        ca9c0c5
6.0        6310725
5.0        62659fb
4.0        a586646
3.2        8f44bed
3.1        d69c059
3.0        7befdcf
2.1        a26a227
2.0        7c507f3
1.1        c143dd9
1.0        04722d8
0.3        f42215e
0.1        34ac742
=======    =================
