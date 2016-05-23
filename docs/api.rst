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

We assume a standard diploid coalescent model in which time is measured in
units of :math:`2N_e` generations, where :math:`N_e` is the Wright-Fisher
population size of our reference population at time zero. As a convenience
when converting from per-generation rates, we can specify :math:`N_e` using
the simulation parameter ``Ne``. All other population
sizes are defined with respect to the this population size. When running
simulations we defined the length in bases :math:`L` of the sequence in question
using the ``length`` parameter. This defines the coordinate space within
which trees and mutations are defined. :math:`L` is a continuous value, and
coordinates can take any value from :math:`0` to :math:`L`.

Mutations occur in an infinite sites process along this sequence, and
mutation rates are specified per generation, per unit of sequence length.
Thus, given the per-generation mutation rate :math:`\mu`, the rate of mutation
over the entire sequence in coalescent time units is :math:`\theta = 4 N_e \mu L`.
It is important to remember these scaling factors when comparing with
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
migration rates. Therefore, the rate of migration in coalescent time units
between a given pair of demes is :math:`4 N_e M_{j, k}`. Each deme has an
initial size :math:`s` that is measured relative to :math:`N_e`, and an
exponential growth rate :math:`\alpha`. The size of a given population at time
:math:`t` in the past (assuming no demographic events have occurred) is
therefore given by :math:`s N_e e^{-\alpha t}`. Demographic events that occur
in the history of the simulated population alter some aspect of this population
configuration at a particular time in the past.


.. warning:: This parameterisation of recombination, mutation and
    migration rates is different to :program:`ms`, which states these
    rates over the entire region and in coalescent time units. The
    motivation for this is to allow the user change the size of the simulated
    region without having to rescale the recombination and mutation rates,
    and to allow users to easily state per-generation recombination rates.
    However, the ``mspms`` command line application is fully :program:`ms`
    compatible.


*******************
Running simulations
*******************

The :func:`.simulate` function provides the primary interface to running
coalescent simulations in msprime.

.. autofunction:: msprime.simulate

++++++++++++++++++++
Population structure
++++++++++++++++++++

Population structure is modelled in ``msprime`` by specifying a fixed
number of demes, with the migration rates between those demes defined by
a migration matrix. Each deme has an ``initial_size`` that defines its
size relative to :math:`N_e`, and a ``growth_rate`` which specifies exponential
growth rate of the sub-population. We must also define the size of the
sample to draw from each deme. The number of populations and their
initial configuration is defined using the ``population_configuration``
parameter to :func:`.simulate`, which takes a list of
:class:`.PopulationConfiguration` instances. Population
IDs are zero indexed, and correspond to their position in the list.

Samples are drawn sequentially from populations in increasing order of
population ID. For example, if we specified an overall sample size of
5, and specify that 2 samples are drawn from population 0 and 3 from
population 1, then individuals 0 and 1 will be initially located in
population 0, and individuals 2, 3 and 4 will be drawn from
population 2.

Given :math:`N` populations, migration matrices are specified using an :math:`N
x N` matrix of deme-to-deme migration rates. See the documentation for
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
that are currently supported.

.. autoclass:: msprime.GrowthRateChangeEvent
.. autoclass:: msprime.SizeChangeEvent
.. autoclass:: msprime.MigrationRateChangeEvent
.. autoclass:: msprime.MassMigrationEvent


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

.. data:: msprime.NULL_NODE = - 1

    Special reserved value, representing the null node. If the parent of a
    given node is null, then this node is a root. Similarly, if the children of
    a node are null, this node is a leaf.

.. autofunction:: msprime.load

.. autoclass:: msprime.TreeSequence()
    :members:

.. autoclass:: msprime.SparseTree()
    :members:

