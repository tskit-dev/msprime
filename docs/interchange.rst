.. _sec_interchange:

#########################
Tree sequence interchange
#########################

The correlated genealogical trees that describe the shared ancestry of set of
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
we describe the HDF5-based file format using by msprime to efficiently
store tree sequences in the `HDF5 file format`_ section.


.. _sec_data_model:

**********
Data model
**********

To begin, here are definitions of some key ideas encountered later.

tree
    A "gene tree", i.e., the genealogical tree describing how each of the
    individuals at the tips of the tree are related to each other.

tree sequence
    A "succinct tree sequence" (or tree sequence, for brevity) is an efficient
    encoding of a sequence of correlated trees. A tree sequence efficiently
    captures the structure shared by adjacent trees, (essentially) storing only
    what differs between them.

node
    Each branching point in each tree is associated with a particular ancestor,
    called "nodes".  Since each node represents a certain ancestor, it has a
    unique ``time``, thought of as her birth time, which determines the height
    of any branching points she is associated with.

sample
    Those nodes in the tree that we have obtained data from.  These are
    distinguished from other nodes by the fact that a tree sequence *must*
    describe the genealogical history of all samples at every point on the
    genome. (See :ref:`sec_node_table_definition` for information on how the sample
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
    back or recurrent mutations, a mutation must also specify it's 'parent'
    mutation.

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

.. todo:: Define migration and provenance types.

A tree sequence can be stored in a collection of six tables:
:ref:`Node <sec_node_table_definition>`,
:ref:`Edge <sec_edge_table_definition>`,
:ref:`Site <sec_site_table_definition>`,
:ref:`Mutation <sec_mutation_table_definition>`,
:ref:`Migration <sec_migration_table_definition>`, and
:ref:`Provenance <sec_provenance_table_definition>`.
The first two store the genealogical
relationships that define the trees; the next two describe where mutations fall
on those trees; the Migration table describes how lineages move across space;
and the Provenance table contains information on where the data came from.
In the following sections we define these components of a tree sequence in
more detail.

.. _sec_table_definitions:

Table definitions
=================


.. _sec_node_table_definition:

Node Table
----------

A **node** defines a specific ancestor that was born at some time in
the past. Every vertex in the marginal trees of a tree sequence corresponds
to exactly one node, and a node may be present in many trees. The
node table contains four columns, of which ``flags`` and ``time`` are
mandatory:

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
flags               uint32              Bitwise flags.
time                double              Birth time of node
population          int32               Birth population of node.
metadata            binary              Node :ref:`sec_metadata_definition`
================    ==============      ===========

The ``time`` column records the birth time of the individual in question,
and is a floating point value. Similarly,
the ``population`` column records the ID of the population where this
individual was born. If not provided, ``population`` defaults to the
null ID (-1).

The ``flags`` column stores information about a particular node, and
is composed of 32 bitwise boolean values. Currently, the only flag defined
is ``IS_SAMPLE = 1``, which defines the sample status of nodes. Marking
a particular node as a sample means, for example, that the mutational state
of the node will be included in the genotypes produced by
:meth:`.TreeSequence.variants`.

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

Each row in an edge table describes the half-open genomic interval
affected ``[left, right)``, and the ``parent`` and ``child`` on that interval.
The ``left`` and ``right`` columns are defined using double precision
floating point values for flexibility. The ``parent`` and ``child``
columns specify integer IDs in the associated :ref:`sec_node_table_definition`.

See the :ref:`sec_edge_requirements` section for details on the properties
required for a valid set of edges.

.. _sec_site_table_definition:

Site Table
----------

A **site** defines a particular location along the genome in which
we are interested in observing the mutational state. The site table
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

The ``ancestral_state`` column specifies the mutational state at the root
of the tree, thus defining the state that nodes inherit (unless mutations
occur). The column stores text character data of arbitrary length.

The ``metadata`` column provides a location for client code to store
information about each site. See the :ref:`sec_metadata_definition` section for
more details on how metadata columns should be used.

See the :ref:`sec_site_requirements` section for details on the properties
required for a valid set of sites.

.. _sec_mutation_table_definition:

Mutation Table
--------------

A **mutation** defines a change of mutational state on a tree at a particular site.
The mutation table contains five columns, of which ``site``, ``node`` and
``derived_state`` are mandatory.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
site                int32               The ID of the site the mutation occurs at.
node                int32               The node this mutation occurs at.
parent              int32               The ID of the parent mutation.
derived_state       char                The mutational state at the defined node.
metadata            char                Mutation :ref:`sec_metadata_definition`.
================    ==============      ===========

The ``site`` column is an integer value defining the ID of the
:ref:`site <sec_site_table_definition>` at which this mutation occured.

The ``node`` column is an integer value defining the ID of the
first :ref:`node <sec_node_table_definition>` in the tree to inherit this mutation.

The ``derived_state`` column specifies the mutational state at the specified node,
thus defining the state that nodes in the subtree inherit (unless further mutations
occur). The column stores text character data of arbitrary length.

The ``parent`` column is an integer value defining the ID of the
mutation from which this mutation inherits. If there is no mutation at the
site in question on the path back to root, then this field is set to the
null ID (-1). (The ``parent`` column is only required in situations
where there are multiple mutations at a given site. For
simple infinite sites mutations, it can be ignored.)

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
individual tracked as a ``node`` (as in a discrete-deme model they are
unlikely to be both a migrant and a most recent common ancestor).  So,
``msprime`` records when a segment of ancestry has moved between
populations.

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
half-open segment of genome affected.


.. todo::
    Document the remaining fields

See the :ref:`sec_migration_requirements` section for details on the properties
required for a valid set of mutations.

..     This ``migration`` records that the ancestor who was alive 2.1 time units
..     in the past from which ``node`` 3 inherited the segment of genome between
..     0.0 and 0.3 migrated from population 0 to population 1.


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
:ref:`sec_tables_api`. However, only tables that fulfil a set of
requirements represent a valid :class:`.TreeSequence` object, and
can be loaded using :func:`.load_tables`. In this
section we list these requirements, and explain their rationale.
Violations of most of these requirements are detected when the
user attempts to load a tree sequence via :func:`.load` or
:func:`.load_tables`, raising an informative error message. Some
more complex requirements may not be detectable at load-time,
and errors may not occur until certain operations are attempted.
These are documented below.


.. _sec_node_requirements:

Node requirements
-----------------

Nodes are the most basic type in a tree sequence, and are not defined with
respect to any other tables. Therefore, the requirements for nodes are
trivial.

- Node times must be non-negative.

There are no requirements regarding the ordering of nodes with respect to time
or any other field. Sorting a set of tables using :func:`.sort_tables` has
no effect on the nodes.

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

Semantically, to ensure a valid tree sequence there is one further requirement:

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
The :func:`.sort_tables` function will ensure that these sortedness
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
The :func:`.sort_tables` function ensures that sites are sorted
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

For simplicity and algorithmic efficiency, mutations must also:

- be sorted by site ID;
- when there are multiple mutations per site, parent mutations must occur
  **before** their children (i.e. if a mutation with ID :math:`x` has
  ``parent`` with ID :math:`y`, then we must have :math:`y < x`).

Violations of these sorting requirements are detected at load time.
The :func:`.sort_tables` function ensures that mutationsare sorted
according to these criteria.

Mutations also have the requirement that they must result in a
change of state. For example, if we have a site with ancestral state
of "A" and a single mutation with derived state "A", then this
mutation does not result in any change of state. This error is
raised at run-time when we reconstruct sample genotypes, for example
in the :meth:`.TreeSequence.variants` iterator.

.. _sec_migration_requirements:

Migration requirements
----------------------

.. todo::
    Add requirements for valid migrations.

.. A valid ``migration``:

.. 1. Has ``time`` strictly between the time of its ``node`` and the time of any
..    ancestral node from which that node inherits on the segment ``[left,
..    right)``.
.. 2. Has the ``population`` of any such ancestor matching ``source``, if another
..    ``migration`` does not intervene.


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

.. todo::
    Update the examples in this section to be a very simple tree sequence
    with (say) 4 nodes and two trees, and include a picture. This
    example can also be used in the binary interchange section also.

.. _sec_node_text_format:

Node text format
================

The node text format must contain the columns ``is_sample`` and
``time``. Optionally, there may also be a ``population`` and
``metadata`` columns. See the :ref:`node table definitions
<sec_node_table_definition>` for details on these columns.

Note that we do not have a ``flags`` column in the text file format, but
instead use ``is_sample`` (which may be 0 or 1). Currently, ``IS_SAMPLE`` is
the only flag value defined for nodes, and as more flags are defined we will
allow for extra columns in the text format.

An example node table::

    is_sample   time    population
    1           0.0     0
    1           0.0     0
    1           0.0     0
    1           0.0     0
    0           0.071   0
    0           0.090   0
    0           0.170   0
    0           0.202   0
    0           0.253   0


.. _sec_edge_text_format:

Edge text format
================

The edge text format must contain the columns ``left``,
``right``, ``parent`` and ``child``.
See the :ref:`edge table definitions <sec_edge_table_definition>`
for details on these columns.

An example edge table::

    left    right   parent  child
    2       10      4       2
    0       10      5       1
    0       7       6       0
    7       10      7       0
    0       2       8       2


.. _sec_site_text_format:

Site text format
================

The site text format must contain the columns ``position`` and
``ancestral_state``. The ``metadata`` column may also be optionally
present. See the
:ref:`site table definitions <sec_site_table_definition>`
for details on these columns.

sites::

    position    ancestral_state
    0.1         A
    8.5         AT

.. _sec_mutation_text_format:

Mutation text format
====================

The mutation text format must contain the columns ``site``,
``node`` and ``derived_state``. The ``parent`` and ``metadata`` columns
may also be optionally present. See the
:ref:`mutation table definitions <sec_site_table_definition>`
for details on these columns.

mutations::

    site    node    derived_state
    0       3       G
    1       6       T
    1       0       A


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
contigous memory. By taking
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

Note that for a table with ``n`` rows, any offset column must have ``n + 1``
values. The values in this column must be nondecreasing, and cannot exceed
the length of the ragged column in question.

.. _sec_hdf5_file_format:

****************
HDF5 file format
****************

To make tree sequence data as efficient and easy as possible to use, we store the
data on disk in a `HDF5 <https://www.hdfgroup.org/HDF5/>`_ based file format.
Using the specification defined here, it should be straightforward to access tree
sequence information produced by ``msprime`` in any language with `HDF5 support
<https://en.wikipedia.org/wiki/Hierarchical_Data_Format#Interfaces>`_.

The file format is broken into a number of groups, and each group
corresponds to one of the tables above (possibly including some extra
information for efficiency). In general, each group will contain a dataset
corresponding to a column in the table in question. All groups must be
present.

To work around limitations in some versions of the HDF5 library, empty
columns are **not** stored. For example, if there is no metadata associated
with nodes, the ``metadata`` column in the node table will be empty, and
the corresponding ``metadata`` dataset will not be present in the HDF5 file.

Variable length data is handled in the same manner as the
:ref:`Tables API <sec_encoding_ragged_columns>`
above: we store two arrays, one containing the flattened data, and another
storing offsets into this array.

The root group contains two attributes, ``format_version`` and ``sequence_length``.
The ``format_version`` is a pair ``(major, minor)`` describing the file format version.
This document describes version 10.0. The ``sequence_length`` attribute defines the
coordinate space over which edges and sites are defined. This must be present
and be greater than or equal to the largest coordinate present.

================    ==============      ======      ===========
Path                Type                Dim         Description
================    ==============      ======      ===========
/format_version     H5T_STD_U32LE       2           The (major, minor) file format version.
/sequence_length    H5T_IEEE_F64LE      1           The maximum value of a sequence coordinate.
================    ==============      ======      ===========

Nodes group
===========

The ``/nodes`` group stores the :ref:`sec_node_table_definition`.

=======================     ==============
Path                        Type
=======================     ==============
/nodes/flags                H5T_STD_U32LE
/nodes/population           H5T_STD_I32LE
/nodes/time                 H5T_IEEE_F64LE
/nodes/metadata             H5T_STD_I8LE
/nodes/metadata_offset      H5T_STD_U32LE
=======================     ==============

Edges group
===========

The ``/edges`` group stores the :ref:`sec_edge_table_definition`.

===================       ==============
Path                      Type
===================       ==============
/edges/left               H5T_IEEE_F64LE
/edges/right              H5T_IEEE_F64LE
/edges/parent             H5T_STD_I32LE
/edges/child              H5T_STD_I32LE
===================       ==============

Indexes group
-------------

The ``/edges/indexes`` group records information required to efficiently
reconstruct the individual trees from the tree sequence. The
``insertion_order`` dataset contains the order in which records must be applied
and the ``removal_order`` dataset the order in which records must be
removed for a left-to-right traversal of the trees.

==============================     ==============
Path                               Type
==============================     ==============
/edges/indexes/insertion_order     H5T_STD_I32LE
/edges/indexes/removal_order       H5T_STD_I32LE
==============================     ==============

Sites group
===========

The sites group stores the :ref:`sec_site_table_definition`.

=============================   ==============
Path                            Type
=============================   ==============
/sites/position                 H5T_IEEE_F64LE
/sites/ancestral_state          H5T_STD_I8LE
/sites/ancestral_state_offset   H5T_STD_U32LE
/sites/metadata                 H5T_STD_I8LE
/sites/metadata_offset          H5T_STD_U32LE
=============================   ==============

Mutations group
===============

The mutations group stores the :ref:`sec_mutation_table_definition`.

===============================  ==============
Path                             Type
===============================  ==============
/mutations/site                  H5T_STD_I32LE
/mutations/node                  H5T_STD_I32LE
/mutations/parent                H5T_STD_I32LE
/mutations/derived_state         H5T_STD_I8LE
/mutations/derived_state_offset  H5T_STD_U32LE
/mutations/metadata              H5T_STD_I8LE
/mutations/metadata_offset       H5T_STD_U32LE
===============================  ==============

Migrations group
================

The ``/migrations`` group stores the :ref:`sec_migration_table_definition`.

===================       ==============
Path                      Type
===================       ==============
/migrations/left          H5T_IEEE_F64LE
/migrations/right         H5T_IEEE_F64LE
/migrations/node          H5T_STD_I32LE
/migrations/source        H5T_STD_I32LE
/migrations/dest          H5T_STD_I32LE
/migrations/time          H5T_IEEE_F64LE
===================       ==============

Provenances group
=================

The provenances group stores the :ref:`sec_provenance_table_definition`.

===============================  ==============
Path                             Type
===============================  ==============
/provenances/timestamp           H5T_STD_I8LE
/provenances/timestamp_offset    H5T_STD_U32LE
/provenances/record              H5T_STD_I8LE
/provenances/record_offset       H5T_STD_U32LE
===============================  ==============


Legacy Versions
===============

Tree sequence files written by older versions of msprime are not readable by
newer versions of msprime. For major releases of msprime, :ref:`sec_msp_upgrade`
will convert older tree sequence files to the latest version.

However many changes to the tree sequence format are not part of major
releases. The table below gives these versions (contained in the root group
attribute, ``format_version`` as a pair ``(major, minor)``).

.. to obtain hashes where versions were changed:
        git log --oneline -L40,41:lib/msprime.h
   then on each hash, to obtain the parent where a merge occured:
        git log --merges --pretty=format:"%h" fc17dbd | head -n 1
   in some cases this didn't work so required hand manipulation. checks were
   done (after checkign out and rebuilding) with:
        python msp_dev.py simulate 10 tmp.hdf5 && h5dump tmp.hdf5 | head

=======    =================
Version    Commit Short Hash
=======    =================
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
