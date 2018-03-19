.. _sec-api:

=================
API Documentation
=================

This is the API documentation for ``msprime``, and provides detailed information
on the Python programming interface. See the :ref:`sec-tutorial` for an
introduction to using this API to run simulations and analyse the results.

****************
Simulation model
****************

The simulation model in ``msprime`` closely follows the classical ``ms``
program. Unlike ``ms``, however, time is measured in generations rather than
"coalescent units". Internally the same simulation algorithm is used, but
``msprime`` provides a translation layer to allow the user input times and
rates in generations. Similarly, the times associated with the trees produced
by ``msprime`` are in measured generations. To enable this translation from
generations into coalescent units and vice-versa, a reference effective
population size must be provided, which is given by the ``Ne`` parameter in the
:func:`.simulate` function. (Note that we assume diploid population sizes
thoughout, since we scale by :math:`4 N_e`.) Population sizes for individual
demes and for past demographic events are defined as absolute values, **not**
scaled by ``Ne``. All migration rates and growth rates are also per generation.

When running simulations we define the length in bases :math:`L` of the
sequence in question using the ``length`` parameter. This defines the
coordinate space within which trees and mutations are defined. :math:`L` is a
continuous value, and coordinates can take any value from :math:`0` to
:math:`L`. Mutations occur in an infinite sites process along this sequence,
and mutation rates are specified per generation, per unit of sequence length.
Thus, given the per-generation mutation rate :math:`\mu`, the rate of mutation
over the entire sequence in coalescent time units is :math:`\theta = 4 N_e \mu
L`. It is important to remember these scaling factors when comparing with
analytical results!

Similarly, recombination rates are per base, per generation in ``msprime``.
Thus, given the per generation crossover rate :math:`r`, the overall rate
of recombination between the ends of the sequence in coalescent time units
is :math:`\rho = 4 N_e r L`. Recombination events occur in a continuous
coordinate space, so that breakpoints do not necessarily occur at integer
locations. However, the underlying recombination model is finite, and the
behaviour of a small number of loci can be modelled using the
:class:`RecombinationMap` class. However, this is considered an advanced
feature and the majority of cases should be well served with the default
recombination model and number of loci.

Population structure is modelled by specifying a fixed number of demes
:math:`d`, and a :math:`d \times d` matrix :math:`M` of per generation
migration rates. Each element of the matrix :math:`M_{j,k}` defines
the fraction of population :math:`j` that consists of migrants from
population :math:`k` in each generation.
Each deme has an initial absolute population size :math:`s`
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
demes, with the migration rates between those demes defined by a migration
matrix. Each deme has an ``initial_size`` that defines its absolute size at
time zero and a per-generation ``growth_rate`` which specifies the exponential
growth rate of the sub-population. We must also define the size of the sample
to draw from each deme. The number of populations and their initial
configuration is defined using the ``population_configurations`` parameter to
:func:`.simulate`, which takes a list of :class:`.PopulationConfiguration`
instances. Population IDs are zero indexed, and correspond to their position in
the list.

Samples are drawn sequentially from populations in increasing order of
population ID. For example, if we specified an overall sample size of 5, and
specify that 2 samples are drawn from population 0 and 3 from population 1,
then individuals 0 and 1 will be initially located in population 0, and
individuals 2, 3 and 4 will be drawn from population 2.

Given :math:`N` populations, migration matrices are specified using an :math:`N
\times N` matrix of deme-to-deme migration rates. See the documentation for
:func:`.simulate` and the `Simulation model`_ section for more details on the
migration rates.

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


*******************
Processing results
*******************

The :class:`.TreeSequence` class represents a sequence of correlated trees
output by a simulation. The :class:`.SparseTree` class represents a single
tree in this sequence.

+++++++++
Constants
+++++++++

.. data:: msprime.NULL_NODE = -1

    Special reserved value, representing the null node. If the parent of a
    given node is null, then this node is a root. Similarly, if the children of
    a node are null, this node is a leaf.

.. data:: msprime.NULL_POPULATION = -1

    Special reserved value, representing the null population ID. If the
    population associated with a particular tree node is not defined,
    or population information was not available in the underlying
    tree sequence, then this value will be returned by
    :meth:`.SparseTree.get_population`.

.. data:: msprime.FORWARD = 1

    Constant representing the forward direction of travel (i.e.,
    increasing coordinate values).

.. data:: msprime.REVERSE = -1

    Constant representing the reverse direction of travel (i.e.,
    decreasing coordinate values).

++++++++++++
Loading data
++++++++++++

There are several methods for loading data into the msprime API. The simplest
and most convenient is the use the :func:`msprime.load` function to load
a :ref:`HDF ancestry file <sec-hdf5-file-format>`. For small scale data
and debugging, it is often convenient to use the :func:`msprime.load_text`
to read data in the :ref:`text file format <sec-text-file-format>`.
The :func:`msprime.load_tables` function efficiently loads large volumes
of data using the :ref:`Tables API <sec-tables-api>`.


.. autofunction:: msprime.load

.. autofunction:: msprime.load_text

.. autofunction:: msprime.load_tables

++++++++++++++++++
TreeSequence class
++++++++++++++++++

.. autoclass:: msprime.TreeSequence()
    :members:

++++++++++++++++
SparseTree class
++++++++++++++++

.. autoclass:: msprime.SparseTree()
    :members:

**********************
Calculating statistics
**********************

The ``msprime`` API provides methods for efficiently calculating
population genetics statistics from a given :class:`.TreeSequence`.

.. autoclass:: msprime.LdCalculator(tree_sequence)
    :members:


.. _sec-tables-api:

***********
Tables API
***********

The :ref:`tables API <sec-binary-interchange>` provides an efficient way of working
with and interchanging :ref:`tree sequence data <sec-data-model>`. Each table
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
in the native :ref:`HDF5 format <sec-hdf5-file-format>`.

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

+++++++++++++
Table classes
+++++++++++++

.. autoclass:: msprime.NodeTable
    :members:

.. autoclass:: msprime.EdgeTable
    :members:

.. autoclass:: msprime.MigrationTable
    :members:

.. autoclass:: msprime.SiteTable
    :members:

.. autoclass:: msprime.MutationTable
    :members:

.. autoclass:: msprime.ProvenanceTable
    :members:

+++++++++++++++
Table functions
+++++++++++++++

.. autofunction:: msprime.sort_tables

.. autofunction:: msprime.simplify_tables

.. autofunction:: msprime.parse_nodes

.. autofunction:: msprime.parse_edges

.. autofunction:: msprime.parse_sites

.. autofunction:: msprime.parse_mutations
