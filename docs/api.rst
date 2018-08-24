.. _sec_api:

=================
API Documentation
=================

This is the API documentation for ``msprime``, and provides detailed information
on the Python programming interface. See the :ref:`sec_tutorial` for an
introduction to using this API to run simulations and analyse the results.

****************
Simulation model
****************

The simulation model in ``msprime`` closely follows the classical ``ms``
program. Unlike ``ms``, however, time is measured in generations rather than
in units of :math:`4 N_e` generations, i.e., "coalescent units".
This means that when simulating a population with diploid effective size :math:`N_e`,
the mean time to coalescence between two samples
in an ``msprime`` simulation will be around :math:`2 N_e`,
while in an ``ms`` simulation, the mean time will be around :math:`0.5`.
Internally, ``msprime`` uses the same algorithm as ``ms``,
and so the ``Ne`` parameter to the :func:`.simulate` function
still acts as a time scaling, and can be set to ``0.5`` to match many theoretical results,
or to ``0.25`` to match ``ms``. Population sizes for each
subpopulation and for past demographic events are also defined as absolute values, **not**
scaled by ``Ne``. All migration rates and growth rates are also per generation.

When running simulations we define the length :math:`L` of the sequence in
question using the ``length`` parameter. This defines the coordinate space
within which trees and mutations are defined. :math:`L` is a continuous value,
so units are arbitrary, and coordinates can take any value from :math:`0` to
:math:`L`. (So, although we recommend setting the units of length to be
analogous to "bases", events can occur at fractional positions.)
Mutations occur in an infinite sites process along this sequence,
and mutation rates are specified per generation, per unit of sequence length.
Thus, given the per-generation mutation rate :math:`\mu`, the rate of mutation
over the entire sequence in coalescent time units is :math:`\theta = 4 N_e \mu
L`. It is important to remember these scaling factors when comparing with
analytical results!

Similarly, recombination rates are per unit of sequence length and per
generation in ``msprime``. Thus, given the per generation crossover rate
:math:`r`, the overall rate of recombination between the ends of the sequence
in coalescent time units is :math:`\rho = 4 N_e r L`. Although breakpoints do
not necessarily occur at integer locations, the underlying recombination model
is finite, and the behaviour of a small number of loci can be modelled using
the :class:`RecombinationMap` class. However, this is considered an advanced
feature and the majority of cases should be well served with the default
recombination model and number of loci.

Population structure is modelled by specifying a fixed number of subpopulations
:math:`d`, and a :math:`d \times d` matrix :math:`M` of per generation
migration rates. Each element of the matrix :math:`M_{j,k}` defines
the fraction of population :math:`j` that consists of migrants from
population :math:`k` in each generation.
Each subpopulation has an initial absolute population size :math:`s`
and a per generation exponential growth rate :math:`\alpha`. The size of a
given population at time :math:`t` in the past (measured in generations) is
therefore given by :math:`s e^{-\alpha t}`. Demographic events that occur in
the history of the simulated population alter some aspect of this population
configuration at a particular time in the past.

.. warning:: This parameterisation of recombination, mutation and
    migration rates is different to :program:`ms`, which states these
    rates over the entire region and in coalescent time units. The
    motivation for this is to allow the user change the size of the simulated
    region without having to rescale the recombination and mutation rates,
    and to also allow users directly state times and rates in units of
    generations. However, the ``mspms`` command line application is
    fully :program:`ms` compatible.

*******************
Running simulations
*******************

The :func:`.simulate` function provides the primary interface to running
coalescent simulations in msprime.

.. autofunction:: msprime.simulate

++++++++++++++++++++
Population structure
++++++++++++++++++++

Population structure is modelled in ``msprime`` by specifying a fixed number of
subpopulations, with the migration rates between those subpopulations defined by a migration
matrix. Each subpopulation has an ``initial_size`` that defines its absolute diploid size at
time zero and a per-generation ``growth_rate`` which specifies the exponential
growth rate of the sub-population. We must also define the number of genomes to
sample from each subpopulation. The number of populations and their initial
configuration is defined using the ``population_configurations`` parameter to
:func:`.simulate`, which takes a list of :class:`.PopulationConfiguration`
instances. Population IDs are zero indexed, and correspond to their position in
the list.

Samples are drawn sequentially from populations in increasing order of
population ID. For example, if we specified an overall sample size of 6, and
specify that 2 samples are drawn from population 0 and 4 from population 1,
then samples 0 and 1 will be initially located in population 0, and
samples 2, 3, 4, and 5 will be drawn from population 2.

Given :math:`N` populations, migration matrices are specified using an :math:`N
\times N` matrix of between-subpopulation migration rates. See the
documentation for :func:`.simulate` and the `Simulation model`_ section for
more details on the migration rates.

.. autoclass:: msprime.PopulationConfiguration

++++++++++++++++++
Demographic Events
++++++++++++++++++

Demographic events change some aspect of the population configuration
at some time in the past, and are specified using the ``demographic_events``
parameter to :func:`.simulate`. Each element of this list must be an
instance of one of the following demographic events
that are currently supported. Note that all times are measured in
generations, all sizes are absolute (i.e., *not* relative to :math:`N_e`),
and all rates are per-generation.

.. autoclass:: msprime.PopulationParametersChange
.. autoclass:: msprime.MigrationRateChange
.. autoclass:: msprime.MassMigration

++++++++++++++++++++++++++++
Debugging demographic models
++++++++++++++++++++++++++++

.. warning:: The ``DemographyDebugger`` class is preliminary, and the API
    is likely to change in the future.

.. autoclass:: msprime.DemographyDebugger
    :members:

++++++++++++++++++++++++++++
Variable recombination rates
++++++++++++++++++++++++++++

.. autoclass:: msprime.RecombinationMap
    :members:

.. _sec_api_simulate_from:

+++++++++++++++++++++++++++++++++++++++++++++
Initialising simulations from a tree sequence
+++++++++++++++++++++++++++++++++++++++++++++

By default ``msprime`` simulations are initialised by specifying a set of samples,
using the ``sample_size`` or  ``samples`` parameters to :func:`.simulate`. This
initialises the simulation with segments of ancestral material covering the
whole sequence. Simulation then proceeds backwards in time until a most recent
common ancestor has been found at all points along this sequence. We can
also start simulations from different initial conditions by using the
``from_ts`` argument to :func:`.simulate`. Informally, we take an 'unfinished'
tree sequence as a parameter to simulate, initialise the simulation
from the state of this tree sequence and then run the simulation until
coalescence. The returned tree sequence is then the result of taking the
input tree sequence and completing the trees using the coalescent.

This is useful for forwards-time simulators such as
`SLiM <https://messerlab.org/slim/>`_ that can output tree sequences. By running
forward-time simulation for a certain number of generations we obtain a
tree sequence, but these trees may not have had sufficient time to
reach a most recent common ancestor. By using the ``from_ts`` argument
to :func:`.simulate` we can combine the best of both forwards- and
backwards-time simulators. The recent past can be simulated forwards
in time and the ancient past by the coalescent. The coalescent
simulation is initialised by the root segments of the
input tree sequence, ensuring that the minimal amount of ancestral
material possible is simulated.

Please see the :ref:`tutorial <sec_tutorial_simulate_from>` for an example of how to use this
feature with a simple forwards-time Wright-Fisher simulator

------------------
Input requirements
------------------

Any tree sequence can be provided as input to this process, but there is a
specific topological requirement that must be met for the simulations to be
statistically correct. To ensure that ancestral segments are correctly associated within chromosomes
when constructing the initial conditions for the coalescent simulation,
forward-time simulators **must** retain the nodes corresponding to the
initial generation. Furthermore, for every sample in the final generation
(i.e. the extant population at the present time) there must be a path
to one of the founder population nodes. (Please see the :ref:`tutorial <sec_tutorial_simulate_from>`
for further explanation of this point and an example.)

-----------------------------
Recombination map limitations
-----------------------------

Because of the way that ``msprime`` handles recombination internally, care must
be taken when specifying recombination when using the ``from_ts`` argument.
If recombination positions are generated in the same way in both the initial
tree sequence and the coalescent simulation, then everything should work.
However, the fine scale details of the underlying recombination model matter,
so matching nonuniform recombination maps between simulators may not be
possible at present. (To make it work, we must ensure that every recombination
breakpoint in ``from_ts`` matches exactly to a possible recombination
breakpoint in msprime's recombination map, which is not guaranteed because of
msprime's discrete recombination model.)

One case in which it is guaranteed to work is if ``from_ts`` has integer
coordinates, and we want to simulate a coalescent with a uniform recombination
rate. In this case, to have a uniform recombination rate ``r`` use::

    L = int(from_ts.sequence_length)
    recomb_map = msprime.RecombinationMap.uniform_map(L, r, L)
    final_ts = mpsrime.simulate(from_ts=from_ts, recomb_map=recomb_map)


++++++++++++++++++++
Simulating mutations
++++++++++++++++++++

When running coalescent simulations it's usually most convenient to use the
``mutation_rate`` argument to the :func:`.simulate` function to throw neutral
mutations down on the trees. However, sometimes we wish to throw mutations
down on an existing tree sequence: for example, if we want to see the outcome
of different random mutational processes on top of a single simulated topology,
or if we have obtained the tree sequence from another program and wish to
overlay neutral mutations on this tree sequence.

.. autoclass:: msprime.InfiniteSites

.. data:: msprime.BINARY == 0

    The binary mutation alphabet where ancestral states are always "0" and
    derived states "1".

.. data:: msprime.NUCLEOTIDES == 1

    The nucleotides mutation alphabet in which ancestral and derived states are
    chosen from the characters "A", "C", "G" and "T".

.. autofunction:: msprime.mutate


************************
Using simulation results
************************

The :class:`.TreeSequence` class represents a sequence of correlated trees
output by a simulation. The :class:`.SparseTree` class represents a single
tree in this sequence.
These classes are the interfaces used to interact with the trees
and mutational information stored in a tree sequence returned from a simulation.
There are also methods for loading data into these objects, either from the native
format using :func:`msprime.load`, or from another sources
using :func:`msprime.load_text` or :meth:`.TableCollection.tree_sequence`.

+++++++++++++++++
Top level-classes
+++++++++++++++++


.. autoclass:: msprime.TreeSequence()
    :members:

.. autoclass:: msprime.SparseTree()
    :members:


+++++++++
Constants
+++++++++

These constants are used in :class:`.TreeSequence` and :class:`.SparseTree`
instances to define special values of various variables.

.. data:: msprime.NULL_NODE == -1

    Special reserved value, representing the null node. If the parent of a
    given node is null, then this node is either a root (the ancestor of
    one or more samples), or the node is not present in the current tree.

.. data:: msprime.NULL_POPULATION == -1

    Special reserved value, representing the null population ID. This value is
    returned if the population associated with a particular node is not
    defined, or population information was not available in the underlying tree
    sequence.

.. data:: msprime.NULL_INDIVIDUAL == -1

    Special reserved value, representing the null individual ID. This value is
    returned if the individual associated with a particular node is not
    defined, or individual information was not available in the underlying tree
    sequence.

.. data:: msprime.NULL_MUTATION == -1

    Special reserved value, representing the null mutation ID. If there is a
    single mutation at a site, or if a mutation at a given site does not have
    an ancestal mutations it inherited from, this is the value returned by
    ``mutation.parent``. See also :class:`.Mutation`.


.. data:: msprime.FORWARD == 1

    Constant representing the forward direction of travel (i.e.,
    increasing genomic coordinate values).

.. data:: msprime.REVERSE == -1

    Constant representing the reverse direction of travel (i.e.,
    decreasing genomic coordinate values).


++++++++++++++++++++++++
Simple container classes
++++++++++++++++++++++++

These classes are simple shallow containers representing the entities defined
in the :ref:`sec_data_model`. These classes are not intended to be instantiated
directly, but are the return types for the various iterators provided by the
:class:`.TreeSequence` and :class:`.SparseTree` classes.

.. autoclass:: msprime.Individual()
    :members:

.. autoclass:: msprime.Node()
    :members:

.. autoclass:: msprime.Edge()
    :members:

.. autoclass:: msprime.Site()
    :members:

.. autoclass:: msprime.Mutation()
    :members:

.. autoclass:: msprime.Variant()
    :members:

.. autoclass:: msprime.Migration()
    :members:

.. autoclass:: msprime.Population()
    :members:

++++++++++++
Loading data
++++++++++++

There are several methods for loading data into a :class:`.TreeSequence`
instance. The simplest and most convenient is the use the :func:`msprime.load`
function to load a :ref:`tree sequence file <sec_tree_sequence_file_format>`. For small
scale data and debugging, it is often convenient to use the
:func:`msprime.load_text` to read data in the :ref:`text file format
<sec_text_file_format>`. The :meth:`.TableCollection.tree_sequence` function
efficiently creates a :class:`.TreeSequence` object from a set of tables
using the :ref:`Tables API <sec_tables_api>`.


.. autofunction:: msprime.load

.. autofunction:: msprime.load_text


**********************
Calculating statistics
**********************

The ``msprime`` API provides methods for efficiently calculating
population genetics statistics from a given :class:`.TreeSequence`.

.. autoclass:: msprime.LdCalculator(tree_sequence)
    :members:


.. _sec_tables_api:

***********
Tables API
***********

The :ref:`tables API <sec_binary_interchange>` provides an efficient way of working
with and interchanging :ref:`tree sequence data <sec_data_model>`. Each table
class (e.g, :class:`.NodeTable`, :class:`.EdgeTable`) has a specific set of
columns with fixed types, and a set of methods for setting and getting the data
in these columns. The number of rows in the table ``t`` is given by ``len(t)``.
Each table supports accessing the data either by row or column. To access the
row ``j`` in table ``t`` simply use ``t[j]``. The value returned by such an
access is an instance of :func:`collections.namedtuple`, and therefore supports
either positional or named attribute access. To access the data in
a column, we can use standard attribute access which will return a numpy array
of the data. For example::

    >>> import msprime
    >>> t = msprime.EdgeTable()
    >>> t.add_row(left=0, right=1, parent=10, child=11)
    0
    >>> t.add_row(left=1, right=2, parent=9, child=11)
    1
    >>> print(t)
    id      left            right           parent  child
    0       0.00000000      1.00000000      10      11
    1       1.00000000      2.00000000      9       11
    >>> t[0]
    EdgeTableRow(left=0.0, right=1.0, parent=10, child=11)
    >>> t[-1]
    EdgeTableRow(left=1.0, right=2.0, parent=9, child=11)
    >>> t.left
    array([ 0.,  1.])
    >>> t.parent
    array([10,  9], dtype=int32)
    >>> len(t)
    2
    >>>

Tables also support the :mod:`pickle` protocol, and so can be easily
serialised and deserialised (for example, when performing parallel
computations using the :mod:`multiprocessing` module). ::

    >>> serialised = pickle.dumps(t)
    >>> t2 = pickle.loads(serialised)
    >>> print(t2)
    id      left            right           parent  child
    0       0.00000000      1.00000000      10      11
    1       1.00000000      2.00000000      9       11

However, pickling will not be as efficient as storing tables
in the native :ref:`format <sec_tree_sequence_file_format>`.

Tables support the equality operator ``==`` based on the data
held in the columns::

    >>> t == t2
    True
    >>> t is t2
    False
    >>> t2.add_row(0, 1, 2, 3)
    2
    >>> print(t2)
    id      left            right           parent  child
    0       0.00000000      1.00000000      10      11
    1       1.00000000      2.00000000      9       11
    2       0.00000000      1.00000000      2       3
    >>> t == t2
    False



.. _sec_tables_api_text_columns:

++++++++++++
Text columns
++++++++++++

As described in the :ref:`sec_encoding_ragged_columns`, working with
variable length columns is somewhat more involved. Columns
encoding text data store the **encoded bytes** of the flattened
strings, and the offsets into this column in two separate
arrays.

Consider the following example::

    >>> t = msprime.SiteTable()
    >>> t.add_row(0, "A")
    >>> t.add_row(1, "BB")
    >>> t.add_row(2, "")
    >>> t.add_row(3, "CCC")
    >>> print(t)
    id      position        ancestral_state metadata
    0       0.00000000      A
    1       1.00000000      BB
    2       2.00000000
    3       3.00000000      CCC
    >>> t[0]
    SiteTableRow(position=0.0, ancestral_state='A', metadata=b'')
    >>> t[1]
    SiteTableRow(position=1.0, ancestral_state='BB', metadata=b'')
    >>> t[2]
    SiteTableRow(position=2.0, ancestral_state='', metadata=b'')
    >>> t[3]
    SiteTableRow(position=3.0, ancestral_state='CCC', metadata=b'')

Here we create a :class:`.SiteTable` and add four rows, each with a different
``ancestral_state``. We can then access this information from each
row in a straightforward manner. Working with the data in the columns
is a little trickier, however::

    >>> t.ancestral_state
    array([65, 66, 66, 67, 67, 67], dtype=int8)
    >>> t.ancestral_state_offset
    array([0, 1, 3, 3, 6], dtype=uint32)
    >>> msprime.unpack_strings(t.ancestral_state, t.ancestral_state_offset)
    ['A', 'BB', '', 'CCC']

Here, the ``ancestral_state`` array is the UTF8 encoded bytes of the flattened
strings, and the ``ancestral_state_offset`` is the offset into this array
for each row. The :func:`.unpack_strings` function, however, is a convient
way to recover the original strings from this encoding. We can also use the
:func:`.pack_strings` to insert data using this approach::

    >>> a, off = msprime.pack_strings(["0", "12", ""])
    >>> t.set_columns(position=[0, 1, 2], ancestral_state=a, ancestral_state_offset=off)
    >>> print(t)
    id      position        ancestral_state metadata
    0       0.00000000      0
    1       1.00000000      12
    2       2.00000000

When inserting many rows with standard infinite sites mutations (i.e.,
ancestral state is "0"), it is more efficient to construct the
numpy arrays directly than to create a list of strings and use
:func:`.pack_strings`. When doing this, it is important to note that
it is the **encoded** byte values that are stored; by default, we
use UTF8 (which corresponds to ASCII for simple printable characters).::

    >>> t_s = msprime.SiteTable()
    >>> m = 10
    >>> a = ord("0") + np.zeros(m, dtype=np.int8)
    >>> off = np.arange(m + 1, dtype=np.uint32)
    >>> t_s.set_columns(position=np.arange(m), ancestral_state=a, ancestral_state_offset=off)
    >>> print(t_s)
    id      position        ancestral_state metadata
    0       0.00000000      0
    1       1.00000000      0
    2       2.00000000      0
    3       3.00000000      0
    4       4.00000000      0
    5       5.00000000      0
    6       6.00000000      0
    7       7.00000000      0
    8       8.00000000      0
    9       9.00000000      0
    >>> t_s.ancestral_state
    array([48, 48, 48, 48, 48, 48, 48, 48, 48, 48], dtype=int8)
    >>> t_s.ancestral_state_offset
    array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10], dtype=uint32)

Here we create 10 sites at regular positions, each with ancestral state equal to
"0". Note that we use ``ord("0")`` to get the ASCII code for "0" (48), and create
10 copies of this by adding it to an array of zeros. We have done this for
illustration purposes: it is equivalent (though slower for large examples) to do
``a, off = msprime.pack_strings(["0"] * m)``.

Mutations can be handled similarly::

    >>> t_m = msprime.MutationTable()
    >>> site = np.arange(m, dtype=np.int32)
    >>> d, off = msprime.pack_strings(["1"] * m)
    >>> node = np.zeros(m, dtype=np.int32)
    >>> t_m.set_columns(site=site, node=node, derived_state=d, derived_state_offset=off)
    >>> print(t_m)
    id      site    node    derived_state   parent  metadata
    0       0       0       1       -1
    1       1       0       1       -1
    2       2       0       1       -1
    3       3       0       1       -1
    4       4       0       1       -1
    5       5       0       1       -1
    6       6       0       1       -1
    7       7       0       1       -1
    8       8       0       1       -1
    9       9       0       1       -1
    >>>


.. _sec_tables_api_binary_columns:

++++++++++++++
Binary columns
++++++++++++++

Columns storing binary data take the same approach as
:ref:`sec_tables_api_text_columns` to encoding
:ref:`variable length data <sec_encoding_ragged_columns>`.
The difference between the two is
only raw :class:`bytes` values are accepted: no character encoding or
decoding is done on the data. Consider the following example::


    >>> t = msprime.NodeTable()
    >>> t.add_row(metadata=b"raw bytes")
    >>> t.add_row(metadata=pickle.dumps({"x": 1.1}))
    >>> t[0].metadata
    b'raw bytes'
    >>> t[1].metadata
    b'\x80\x03}q\x00X\x01\x00\x00\x00xq\x01G?\xf1\x99\x99\x99\x99\x99\x9as.'
    >>> pickle.loads(t[1].metadata)
    {'x': 1.1}
    >>> print(t)
    id      flags   population      time    metadata
    0       0       -1      0.00000000000000        cmF3IGJ5dGVz
    1       0       -1      0.00000000000000        gAN9cQBYAQAAAHhxAUc/8ZmZmZmZmnMu
    >>> t.metadata
    array([ 114,   97,  119,   32,   98,  121,  116,  101,  115, -128,    3,
            125,  113,    0,   88,    1,    0,    0,    0,  120,  113,    1,
             71,   63,  -15, -103, -103, -103, -103, -103, -102,  115,   46], dtype=int8)
    >>> t.metadata_offset
    array([ 0,  9, 33], dtype=uint32)


Here we add two rows to a :class:`.NodeTable`, with different
:ref:`metadata <sec_metadata_definition>`. The first row contains a simple
byte string, and the second contains a Python dictionary serialised using
:mod:`pickle`. We then show several different (and seemingly incompatible!)
different views on the same data.

When we access the data in a row (e.g., ``t[0].metadata``) we are returned
a Python bytes object containing precisely the bytes that were inserted.
The pickled dictionary is encoded in 24 bytes containing unprintable
characters, and when we unpickle it using :func:`pickle.loads`, we obtain
the original dictionary.

When we print the table, however, we see some data which is seemingly
unrelated to the original contents. This is because the binary data is
`base64 encoded <https://en.wikipedia.org/wiki/Base64>`_ to ensure
that it is print-safe (and doesn't break your terminal). (See the
:ref:`sec_metadata_definition` section for more information on the
use of base64 encoding.).

Finally, when we print the ``metadata`` column, we see the raw byte values
encoded as signed integers. As for :ref:`sec_tables_api_text_columns`,
the ``metadata_offset`` column encodes the offsets into this array. So, we
see that the first metadata value is 9 bytes long and the second is 24.

The :func:`pack_bytes` and :func:`unpack_bytes` functions are also useful
for encoding data in these columns.

+++++++++++++
Table classes
+++++++++++++

This section describes the methods and variables available for each
table class. For description and definition of each table's meaning
and use, see :ref:`the table definitions <sec_table_definitions>`.

.. Overriding the default signatures for the tables here as they will be
.. confusing to most users.

.. autoclass:: msprime.IndividualTable()
    :members:
    :inherited-members:

.. autoclass:: msprime.NodeTable()
    :members:
    :inherited-members:

.. autoclass:: msprime.EdgeTable()
    :members:
    :inherited-members:

.. autoclass:: msprime.MigrationTable()
    :members:
    :inherited-members:

.. autoclass:: msprime.SiteTable()
    :members:
    :inherited-members:

.. autoclass:: msprime.MutationTable()
    :members:
    :inherited-members:

.. autoclass:: msprime.PopulationTable()
    :members:
    :inherited-members:

.. autoclass:: msprime.ProvenanceTable()
    :members:
    :inherited-members:

+++++++++++++++++
Table Collections
+++++++++++++++++

Each of the table classes defines a different aspect of the structure of
a tree sequence. It is convenient to be able to refer to a set of these
tables which together define a tree sequence. We
refer to this grouping of related tables as a ``TableCollection``.
The :class:`.TableCollection` and :class:`.TreeSequence` classes are
deeply related. A ``TreeSequence`` instance is based on the information
encoded in a ``TableCollection``. Tree sequences are **immutable**, and
provide methods for obtaining trees from the sequence. A ``TableCollection``
is **mutable**, and does not have any methods for obtaining trees.
The ``TableCollection`` class essentially exists to allow the
dynamic creation of tree sequences.

.. autoclass:: msprime.TableCollection(sequence_length=0)
    :members:


+++++++++++++++
Table functions
+++++++++++++++

.. autofunction:: msprime.parse_nodes

.. autofunction:: msprime.parse_edges

.. autofunction:: msprime.parse_sites

.. autofunction:: msprime.parse_mutations

.. autofunction:: msprime.pack_strings

.. autofunction:: msprime.unpack_strings

.. autofunction:: msprime.pack_bytes

.. autofunction:: msprime.unpack_bytes


.. _sec_provenance_api:

**************
Provenance API
**************

We provide some preliminary support for validating JSON documents against the
:ref:`provenance schema <sec_provenance>`. Programmatic access to provenance
information is planned for future versions.

.. autofunction:: msprime.validate_provenance

.. autoexception:: msprime.ProvenanceValidationError



