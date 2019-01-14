.. _sec_api:

=================
API Documentation
=================

This is the API documentation for ``msprime``, and provides detailed information
on the Python programming interface. See the :ref:`sec_tutorial` for an
introduction to using this API to run simulations.
See the `tskit documentation <https://tskit.readthedocs.io/en/stable>`_ for
information on how to use the
`tskit Python API <https://tskit.readthedocs.io/en/stable/python-api.html>`_
to analyse simulation results.

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
the :class:`.RecombinationMap` class. However, this is considered an advanced
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

********************
Population structure
********************

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

******************
Demographic Events
******************

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

****************************
Variable recombination rates
****************************

.. autoclass:: msprime.RecombinationMap
    :members:

.. _sec_api_simulate_from:

*********************************************
Initialising simulations from a tree sequence
*********************************************

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

++++++++++++++++++
Input requirements
++++++++++++++++++

Any tree sequence can be provided as input to this process, but there is a
specific topological requirement that must be met for the simulations to be
statistically correct. To ensure that ancestral segments are correctly associated within chromosomes
when constructing the initial conditions for the coalescent simulation,
forward-time simulators **must** retain the nodes corresponding to the
initial generation. Furthermore, for every sample in the final generation
(i.e. the extant population at the present time) there must be a path
to one of the founder population nodes. (Please see the :ref:`tutorial <sec_tutorial_simulate_from>`
for further explanation of this point and an example.)

+++++++++++++++++++++++++++++
Recombination map limitations
+++++++++++++++++++++++++++++

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


********************
Simulating mutations
********************

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
