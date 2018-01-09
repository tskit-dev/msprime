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
`Text file formats`_ section. The `Tables API`_ section then describes
the efficient Python API for table interchange using numpy arrays. Finally,
we describe the HDF5-based file format using by msprime to efficiently
store tree sequences in the `HDF5 file format`_ section.


.. _sec-data-model:

**********
Data model
**********

To begin, here are definitions of some key ideas encountered later.  This will
define the terminology, as well as giving properties of the tables that these
are stored in.


.. These are properties that can be assumed when writing methods
.. that operate on an ``msprime`` tree sequence; the function ``sort_tables`` is
.. provided to put unsorted tables in the proper order.

Defintions
==========

First are those that describe genealogical relationships:

tree
    A "gene tree", i.e., the genealogical tree describing how each of the
    individuals at the tips of the tree are related to each other.  A "tree
    sequence" contains information sufficient to reconstruct the genealogical
    tree relating all samples to each other at any point along the genome.

node
    Each branching point in each tree is associated with a particular ancestor,
    called "nodes".  Since each node represents a certain ancestor, it has a
    unique ``time``, thought of as her birth time, which determines the height
    of any branching points she is associated with.  A given node will be
    associated with branching points of all trees across a region if that node
    is the most recent common ancestor to the subtending tips across that
    region.  For each node, we record::

        (flags, population, time)

    where ``flags`` records information about the ancestor; ``population`` is
    the integer ID of the ancestor's (birth) population, and ``time`` is how
    long ago the ancestor was born.  Each node also has a unique (integer) ID,
    but this is *not* recorded explicitly - rather, the individual's ID is
    given by their position in the tree sequence's node table.

samples
    Those nodes in the tree that we have obtained data from.  These are
    distinguished from other nodes by the fact that a tree sequence *must*
    describe the genealogical history of all samples at every point on the
    genome.  These are a special kind of node, having ``flags`` set to 1 (as a
    binary mask).

edge
    Tree sequences are constructed by specifying over which segments of genome
    which nodes inherit from which other nodes.  This information is stored by
    recording::

        (left, right, parent, child)

    where each node in ``child`` inherits from the node ``parent``
    on the half-open interval of chromosome ``[left, right)``.


Here are the formal requirements for a set of nodes and edges to make sense,
and to allow ``msprime``'s algorithms to work properly.

.. _sec-encoding-ragged-columns:

Encoding ragged columns
=======================

    **todo: This section will define how to work with ragged columns. It's not clear
    yet where it should be placed.**


.. _sec-node-table-definition:

Node Table
==========

.. todo Clear up distinction between flags and is_sample.

Example table:

===    =========  ==========   ====
id     is_sample  population   time
===    =========  ==========   ====
0      1          0            0.0
1      1          1            0.0
2      0          0            0.0
3      1          0            0.5
4      0          2            2.1
===    =========  ==========   ====

Node IDs are *not* recorded; rather the `id` column shows the row index, so
that the `k`-th row describes the node whose ID is `k`.  `is_sample`
records whether the node is a sample (=1) or not (=0).  `population` is an
integer population ID, and `time` is the time since that individual was
born, as a float.

Requirements:

1. All birth times must be greater than or equal to zero.

It is not required that the `time` column be ordered or that all samples
must be at the top.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
flags               uint32              Bitwise flags.
time                double              Birth time of node
population          int32               Birth population of node.
metadata            char                Node :ref:`sec-metadata-definition`
================    ==============      ===========

.. _sec-edge-table-definition:

Edge Table
==========

=====   =====   ======  =====
left	right	parent	child
=====   =====   ======  =====
0.0     0.4     3       0
0.0     0.4     3       2
0.4     1.0     3       0
0.4     1.0     3       1
0.4     1.0     3       2
0.0     0.4     4       1
0.0     0.4     4       3
=====   =====   ======  =====

Each row in an edge table describes the half-open genomic interval
affected `[left, right)`, the `parent` and the `child` on that interval.


================    ==============      ===========
Column              Type                Description
================    ==============      ===========
left                double              Left coordinate of the edge (inclusive).
right               double              Right coordinate of the edge (exclusive).
parent              int32               Parent node ID.
child               int32               Child node ID.
================    ==============      ===========


.. _sec-migration-table-definition:

Migration Table
===============

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


.. _sec-site-table-definition:

Site Table
==========

Rather than storing a position on the genome directly, a ``mutation``
stores the index of a ``site``, that describes that position.  This is to
allow efficient processing of multiple mutations at the same genomic
position.  A ``site`` records a position on the genome where a mutation has
occurred along with the ancestral state (i.e., the state at the root of the
tree at that position)::

    id	position	ancestral_state
    0	0.1	        0

As with nodes, the ``id`` is not stored directly, but is implied by its
index in the site table.


To allow for efficent algorithms, it is required that

8. Sites are sorted by increasing position,
9. and mutations are sorted by site.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
position            double              Position of site in genome coordinates.
ancestral_state     char                The state at the root of the tree.
metadata            char                Site :ref:`sec-metadata-definition`.
================    ==============      ===========


.. _sec-mutation-table-definition:

Mutation Table
==============

This type records a mutation that has occurred at some point in the
genealogical history.  Each mutation is associated with a particular
``node`` (i.e., a particular ancestor), so that any sample which inherits
from that node will also inherit that mutation, unless another mutation
intervenes.  The type records::

    site	node	derived_state
    0	    14	    1

Here ``site`` is the index of the ``site`` at which the mutation occurred,
``node`` records the ID of the ancestral node associated with the mutation,
and ``derived_state`` is the allele that any sample inheriting from that
node at this site will have if another mutation does not intervene.  The
``node`` is not necessarily the ancestor in whom the mutation occurred, but
rather the ancestor at the bottom of the branch in the tree at that site on
which the mutation occurred.

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
site                int32               The ID of the site the mutation occurs at.
node                int32               The node this mutation occurs at.
parent              int32               The ID of the parent mutation.
derived_state       char                The mutational state at the defined node.
metadata            char                Site :ref:`sec-metadata-definition`.
================    ==============      ===========


.. _sec-provenance-table-definition:

Provenance Table
================

================    ==============      ===========
Column              Type                Description
================    ==============      ===========
timestamp           char                Timestamp in `ISO-8601 <https://en.wikipedia.org/wiki/ISO_8601>`_ format.
record              char                Provenance record.
================    ==============      ===========


.. todo: move this to somewhere else.
.. In addition to genealogical relationships, ``msprime`` generates and stores
.. mutations.  Associating these with nodes means that a variant shared by many
.. individuals need only be stored once, allowing retrieval and processing of
.. variant information much more efficiently than if every individual's genotype
.. was stored directly.

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



.. _sec-text-file-format:

*****************
Text file formats
*****************


An example of a simple tree sequence for four samples with
three distinct trees is as follows.

nodes::

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

edges::

    left    right   node    children
    2       10      4       2,3
    0       2       5       1,3
    2       10      5       1,4
    0       7       6       0,5
    7       10      7       0,5
    0       2       8       2,6


This example is equivalent to the tree sequence illustrated in Figure 4 of
the `PLoS Computational Biology paper
<http://dx.doi.org/10.1371/journal.pcbi.1004842>`_. Nodes are given here in
time order (since this is a backwards-in-time tree sequence), but they may
be allocated in any order. In particular, left-to-right tree sequences are
fully supported.

An example of a ``sites`` and ``mutations`` file for the tree sequence
defined in the previous example is as follows.

sites::

    position    ancestral_state
    0.1         0
    8.5         0

mutations::

    site    node    derived_state
    0       3       1
    1       6       1
    1       0       0


.. _sec-tables-api:

**********
Tables API
**********


.. _sec-variable-length-columns:

Variable length columns
=======================

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
