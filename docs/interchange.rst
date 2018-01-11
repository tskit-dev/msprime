.. _sec-interchange:

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


.. _sec-data-model:

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
    genome. (See :ref:`sec-node-table-definition` for information on how the sample
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

A tree sequence can be stored in a collection of six tables: Node, Edge, Site,
Mutation, Migration, and Provenance. The first two store the genealogical
relationships that define the trees; the next two describe where mutations fall
on those trees; the Migration table describes how lineages move across space;
and the Provenance table contains information on where the data came from.
In the following sections we define these components of a tree sequence in
more detail.

Table definitions
=================


.. _sec-node-table-definition:

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
metadata            binary              Node :ref:`sec-metadata-definition`
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
:meth:``TreeSequence.variants``.

For convenience, the :ref:`text format <sec-text-file-format>` for nodes
decomposes the ``flags`` value into it's separate values. Thus, in the
text format we have a column for ``is_sample``, which corresponds to the
the ``flags`` column in the underlying table. As more flags values are
defined, these will be added to the text file format.

The ``metadata`` column provides a location for client code to store
information about each node. See the :ref:`sec-metadata-definition` section for
more details on how metadata columns should be used.

.. note::
    The distinction between ``flags`` and ``metadata`` is that flags
    holds information about a node that the library understands, whereas
    metadata holds information about a node that the library *does not*
    understand. Metadata is for storing auxiliarly information that is
    not necessary for the core tree sequence algorithms.


.. _sec-edge-table-definition:

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
affected ``[left, right)``, the ``parent`` and the ``child`` on that interval.
The ``left`` and ``right`` columns are defined using double precision
floating point values for flexibility. The ``parent`` and ``child``
columns specify integer IDs in the associated :ref:`sec-node-table-definition`.


.. _sec-site-table-definition:

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
metadata            binary              Site :ref:`sec-metadata-definition`.
================    ==============      ===========

The ``position`` column is a floating point value defining the location
of the site in question along the genome.

The ``ancestral_state`` column specifies the mutational state at the root
of the tree, thus defining the state that nodes inherit (unless mutations
occur). The column stores text character data of arbitrary length.

The ``metadata`` column provides a location for client code to store
information about each site. See the :ref:`sec-metadata-definition` section for
more details on how metadata columns should be used.


.. _sec-mutation-table-definition:

Mutation Table
--------------

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
site                int32               The ID of the site the mutation occurs at.
node                int32               The node this mutation occurs at.
parent              int32               The ID of the parent mutation.
derived_state       char                The mutational state at the defined node.
metadata            char                Mutation :ref:`sec-metadata-definition`.
================    ==============      ===========


.. _sec-migration-table-definition:

Migration Table
---------------

In simulations, trees can be thought of as spread across space, and it is
helpful for inferring demographic history to record this history.  This is
stored using the following type.

migration
    Migrations are performed by individual ancestors, but most likely not by an
    individual tracked as a ``node`` (as in a discrete-deme model they are
    unlikely to be both a migrant and a most recent common ancestor).  So,
    ``msprime`` records when a segment of ancestry has moved between
    populations::

        left    right   node    source  dest    time
        0.0     0.3     3       0       1       2.1

    This ``migration`` records that the ancestor who was alive 2.1 time units
    in the past from which ``node`` 3 inherited the segment of genome between
    0.0 and 0.3 migrated from population 0 to population 1.

A valid ``migration``:

1. Has ``time`` strictly between the time of its ``node`` and the time of any
   ancestral node from which that node inherits on the segment ``[left,
   right)``.
2. Has the ``population`` of any such ancestor matching ``source``, if another
   ``migration`` does not intervene.

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



.. _sec-provenance-table-definition:

Provenance Table
----------------

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
timestamp           char                Timestamp in `ISO-8601 <https://en.wikipedia.org/wiki/ISO_8601>`_ format.
record              char                Provenance record.
================    ==============      ===========



.. _sec-metadata-definition:

Metadata
========

.. _sec-valid-tree-sequence-requirements:

Valid tree sequence requirements
================================

**Explain and list the requirements for a set of tables to form a valid tree
sequence**.

.. _sec-structural-requirements:

Structural requirements
-----------------------


1. All birth times must be greater than or equal to zero.

To disallow time travel and multiple inheritance:

1. Offspring must be born after their parents (and hence, no loops).
2. The set of intervals on which each individual is a child must be disjoint.

and for algorithmic reasons:

3. The leftmost endpoint of each chromosome is 0.0.
4. Node times must be strictly greater than zero.


.. _sec-ordering-requirements:

Ordering requirements
---------------------

Edges are ordered by

- time of parent, then
- parent node ID, then
- child node ID, then
- left endpoint.

Sites are ordered by position, and Mutations are ordered by site.

5. Edges must be sorted in nondecreasing time order.
6. The set of intervals on which each individual is a parent must be disjoint.

A set of tables satisfying requirements 1-4 can be transformed into a completely
valid set of tables by applying first ``sort_tables()`` (which ensures 5)
and then ``simplify_tables()`` (which ensures 6).

Note that since each node time is equal to the (birth) time of the
corresponding parent, time is measured in clock time (not meioses).


To allow for efficent algorithms, it is required that

8. Sites are sorted by increasing position,
9. and mutations are sorted by site.

.. _sec-text-file-format:

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

.. _sec-node-text-format:

Node text format
================

The node text format must contain the columns ``is_sample`` and
``time``. Optionally, there may also be a ``population`` and
``metadata`` columns. See the :ref:`node table definitions
<sec-node-table-definition>` for details on these columns.

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


.. _sec-edge-text-format:

Edge text format
================

The edge text format must contain the columns ``left``,
``right``, ``parent`` and ``child``.
See the :ref:`edge table definitions
<sec-edge-table-definition>` for details on these columns.

An example edge table::

    left    right   parent  child
    2       10      4       2
    0       10      5       1
    0       7       6       0
    7       10      7       0
    0       2       8       2


.. _sec-site-text-format:

Site text format
================

The site text format must contain the columns ``position`` and
``ancestral_state``. The ``metadata`` column may also be optionally
present. See the :ref:`site table definitions
<sec-site-table-definition>` for details on these columns.

sites::

    position    ancestral_state
    0.1         A
    8.5         AT

.. _sec-mutation-text-format:

Mutation text format
====================

The mutation text format must contain the columns ``site``,
``node`` and ``derived_state``. The ``parent`` and ``metadata`` columns
may also be optionally present. See the :ref:`mutation table definitions
<sec-site-table-definition>` for details on these columns.

mutations::

    site    node    derived_state
    0       3       G
    1       6       T
    1       0       A


.. _sec-binary-interchange:

******************
Binary interchange
******************

In this section we describe the high-level details of the API for interchanging
table data via numpy arrays. Please see the :ref:`sec-tables-api` for detailed
description of the functions and methods.


.. _sec-encoding-ragged-columns:

Encoding ragged columns
=======================

    **todo: This section will define how to work with ragged columns. It's not clear
    This has been referred to from elsewhere, but we should probably rename it**

.. _sec-variable-length-columns:

Variable length columns
=======================

.. Keeping this for now as we're referring to it below. Merge these two into one
   section and get rid of duplicate refs.



.. Sorting and simplifying tables
.. ==============================

.. Tables that are noncontradictory but do not satisfy all algorithmic requirements
.. listed above may be converted to a TreeSequence by first sorting, then simplifying
.. them (both operate on the tables **in place**):

.. .. autofunction:: msprime.sort_tables(nodes, edges[, migrations, sites, mutations, edge_start])

.. **Note:** the following function is more general than
.. ``TreeSequence.simplify()``, since it can be applied to tables not satisfying
.. all criteria above (and that hence could not be loaded into a TreeSequence).



.. NodeTable
.. =========

.. .. autoclass:: msprime.NodeTable


.. EdgeTable
.. ============

.. .. autoclass:: msprime.EdgeTable


.. SiteTable
.. =========

.. .. autoclass:: msprime.SiteTable


.. MutationTable
.. =============

.. .. autoclass:: msprime.MutationTable


.. Import and export
.. =================

.. This section describes how to extract tables from a ``TreeSequence``, and how
.. to construct a ``TreeSequence`` from tables.  Since tree sequences are
.. immutible, often the best way to modify a ``TreeSequence`` is something along
.. the lines of (for ``ts`` a ``TreeSequence``)::

..     nodes = msprime.NodeTable()
..     edges = msprime.EdgeTable()
..     ts.dump_tables(nodes=nodes, edges=edges)
..     # (modify nodes and edges)
..     ts.load_tables(nodes=nodes, edges=edges)


.. .. automethod:: msprime.TreeSequence.load_tables

.. .. automethod:: msprime.TreeSequence.dump_tables
..    :noindex:


.. _sec-hdf5-file-format:

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
:ref:`Tables API <sec-variable-length-columns>`
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

The ``/nodes`` group stores the :ref:`sec-node-table-definition`.

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

The ``/edges`` group stores the :ref:`sec-edge-table-definition`.

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

The sites group stores the :ref:`sec-site-table-definition`.

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

The mutations group stores the :ref:`sec-mutation-table-definition`.

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

The ``/migrations`` group stores the :ref:`sec-migration-table-definition`.

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

The provenances group stores the :ref:`sec-provenance-table-definition`.

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
newer versions of msprime. For major releases of msprime, :ref:`sec-msp-upgrade`
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
