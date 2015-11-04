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
units of :math:`2N_0` generations, where :math:`N_0` is the Wright-Fisher
population size at time zero. Recombination is modelled via :math:`m` discrete
loci, between which recombination occurs. These loci as can be viewed as
classical genetic loci consisting of several kilobases of non-recombining
sequence or as single nucleotides, depending on what we wish to model. The
:attr:`scaled_recombination_rate` specifies the rate at which recombination
happens between *adjacent pairs* of loci per :math:`4N_0` generations. The
:attr:`num_loci` parameter determines the coordinate space for a simulated tree
sequence; trees always begin and end at discrete locations in this space.

Mutations follow the infinite sites model, and fall in a continuous coordinate
space from :math:`0` to :math:`m`. The mutations rate is specified by the
:attr:`scaled_mutation_rate` parameter, which gives the rate at which mutations
fall per unit length per :math:`4N_0` generations.

.. warning:: This parameterisation of recombination and mutation rate is
    different to :program:`ms`, which states these rates over the entire region. The
    motivation for this is to allow the user change the size of the simulated
    region without having to rescale the recombination and mutation rates.
    However, the ``mspms`` command line application is fully :program:`ms`
    compatible.

++++++++++++++++++
Demographic events
++++++++++++++++++

Demographic events in the history of the population are specified via
:class:`PopulationModel` instances. Each population model has a starting time
(which may be zero), and some associated parameters. Arbitrarily complex
demographic histories can specified by providing a list of these population
model to :func:`simulate`. This method of specifying demographic events
closely follows the ``-eG`` and ``-eN`` options to :program:`ms`, and the models
are equivalent to these options.

.. autoclass:: msprime.ConstantPopulationModel

.. autoclass:: msprime.ExponentialPopulationModel


*******************
Running simulations
*******************

The :func:`.simulate` function provides the basic interface to running the
coalescent with recombination simulation. The :func:`.simulate_tree` function
provides a simplified interface for the special case in which we have no
recombination.

.. autofunction:: msprime.simulate

.. autofunction:: msprime.simulate_tree


*******************
Processing results
*******************

The :class:`.TreeSequence` class represents a sequence of correlated trees
output by a simulation. The :class:`.SparseTree` class represents a single
tree in this sequence.

.. autofunction:: msprime.load

.. autoclass:: msprime.TreeSequence()
    :members:

.. autoclass:: msprime.SparseTree()
    :members:

