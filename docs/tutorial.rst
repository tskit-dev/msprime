.. _sec_tutorial:

========
Tutorial
========

This is the tutorial for the Python interface to the ``msprime``
library. Detailed :ref:`sec_api` is also available for this
library. An :program:`ms` compatible :ref:`command line interface <sec_cli>`
is also available if you wish to use ``msprime`` directly within
an existing work flow.


****************
Simulating trees
****************

Running simulations is very straightforward in ``msprime``::


    >>> import msprime
    >>> tree_sequence = msprime.simulate(sample_size=6, Ne=1000)
    >>> tree = tree_sequence.first()
    >>> print(tree.draw(format="unicode"))

        10
     ┏━━┻━┓
     ┃    9
     ┃  ┏━┻━┓
     8  ┃   ┃
    ┏┻┓ ┃   ┃
    ┃ ┃ ┃   7
    ┃ ┃ ┃ ┏━┻┓
    ┃ ┃ ┃ ┃  6
    ┃ ┃ ┃ ┃ ┏┻┓
    3 5 0 4 1 2


Here, we simulate the coalescent for a sample of size six
with an effective population size of 1000 diploids,
and then print out a depiction of the resulting tree.
The :func:`.simulate` function returns a
:class:`.TreeSequence` object, which provides a very
efficient way to access the correlated trees in simulations
involving recombination. In this example we know that
there can only be one tree because we have not provided
a value for ``recombination_rate``, and it
defaults to zero.
Therefore, we access the only tree in the
sequence using the call ``tree_sequence.first()``.
Finally, we draw a simple depiction of the tree to the terminal
using the :meth:`~.SparseTree.draw` method.

Genealogical trees record the lines of descent along which genomes
have been inherited. Since diploids have two copies of each autosomal
chromosome, diploid individuals contain two such lines of descent:
the simulation above provides the genealogical history of only three diploids.

Trees are represented within ``msprime`` in a slightly unusual way. In
the majority of libraries dealing with trees, each node is
represented as an object in memory and the relationship
between nodes as pointers between these objects. In ``msprime``,
however, nodes are *integers*.
In the tree above, we can see that the leaves of the tree
are labelled with 0 to 5, and all the internal nodes of the tree
are also integers with the root of the tree being 10.

We can easily trace our path
back to the root for a particular sample using the
:meth:`~.SparseTree.parent` method::

    >>> u = 2
    >>> while u != msprime.NULL_NODE:
    >>>     print("node {}: time = {}".format(u, tree.time(u)))
    >>>     u = tree.parent(u)
    node 2: time = 0.0
    node 6: time = 11.59282234272971
    node 7: time = 129.57841077196494
    node 9: time = 1959.4591339636365
    node 10: time = 5379.737460469677


In this code chunk we iterate up the tree starting at node 0 and
stop when we get to the root. We know that a node is a root
if its parent is :const:`msprime.NULL_NODE`, which is a special
reserved node. (The value of the null node is -1, but we recommend
using the symbolic constant to make code more readable.) We also use
the :meth:`~.SparseTree.time` method to get the time
for each node, which corresponds to the time in generations
at which the coalescence event happened during the simulation.
We can also obtain the length of a branch joining a node to
its parent using the :meth:`~.SparseTree.branch_length`
method::

    >>> print(tree.branch_length(6))
    117.98558842923524

The branch length for node 6 is about 118 generations, since
the birth times of node 6 was 11 generations ago, and the birth time of its
parent, node 7, was around 129 generations ago.
It is also
often useful to obtain the total branch length of the tree, i.e.,
the sum of the lengths of all branches::

    >>> print(tree.total_branch_length)
    13238.125493096279

*************
Recombination
*************

Simulating the history of a single locus is a very useful, but we are most
often interesting in simulating the history of our sample across large genomic
regions under the influence of recombination. The ``msprime`` API is
specifically designed to make this common requirement both easy and efficient.
To model genomic sequences under the influence of recombination we have
two parameters to the :func:`.simulate()` function.
The ``length`` parameter specifies the length of the simulated sequence,
and is a floating point number, so recombination (and mutation) can
occur at any location along the sequence (the units are arbitrary).
If ``length`` is not supplied, it is assumed to be 1.0. The ``recombination_rate``
parameter specifies the rate of crossing over per unit of length per generation,
and is zero by default. See the :ref:`sec_api` for a discussion of the precise
recombination model used.

Here, we simulate the trees across over a 10kb region with a recombination
rate of :math:`2 \times 10^{-8}` per base per generation, with a diploid
effective population size of 1000::

    >>> tree_sequence = msprime.simulate(
    ...     sample_size=6, Ne=1000, length=1e4, recombination_rate=2e-8)
    >>> for tree in tree_sequence.trees():
    ...     print("-" * 20)
    ...     print("tree {}: interval = {}".format(tree.index, tree.interval))
    ...     print(tree.draw(format="unicode"))
    --------------------
    tree 0: interval = (0.0,  6016.224463474058)
       11
    ┏━━┻━━┓
    ┃     10
    ┃  ┏━━┻━┓
    ┃  ┃    9
    ┃  ┃  ┏━┻┓
    ┃  7  ┃  ┃
    ┃ ┏┻┓ ┃  ┃
    ┃ ┃ ┃ ┃  6
    ┃ ┃ ┃ ┃ ┏┻┓
    3 0 1 2 4 5

    --------------------
    tree 1: interval = (6016.224463474058, 10000.0)
         10
      ┏━━┻━━┓
      9     ┃
    ┏━┻┓    ┃
    ┃  ┃    8
    ┃  ┃  ┏━┻┓
    ┃  ┃  ┃  7
    ┃  ┃  ┃ ┏┻┓
    ┃  6  ┃ ┃ ┃
    ┃ ┏┻┓ ┃ ┃ ┃
    2 4 5 3 0 1

In this example, we use the :meth:`~.TreeSequence.trees`
method to iterate over the trees in the sequence. For each tree
we print out its index (i.e., its position in the sequence) and
the interval the tree covers (i.e., the genomic
coordinates which all share precisely this tree) using the
:attr:`~.SparseTree.index` and :attr:`~.SparseTree.interval` attributes.
Thus, the first tree covers the
first 6kb of sequence and the second tree covers the remaining 4kb.
We can see
that these trees share a great deal of their structure, but that there are
also important differences between the trees.

.. warning:: Do not store the values returned from the
    :meth:`~.TreeSequence.trees` iterator in a list and operate
    on them afterwards! For efficiency reasons ``msprime`` uses the same
    instance of :class:`.SparseTree` for each tree in the sequence
    and updates the internal state for each new tree. Therefore, if you store
    the trees returned from the iterator in a list, they will all refer
    to the same tree.


*********
Mutations
*********

Mutations are generated in ``msprime`` by throwing mutations down
on the branches of trees at a particular rate. The mutations are
generated under the infinite sites model, and so each mutation
occurs at a unique (floating point) point position along the
genomic interval occupied by a tree. The mutation rate for simulations
is specified using the ``mutation_rate`` parameter of
:func:`.simulate`. For example, the following chunk simulates 50kb
of nonrecombining sequence with a mutation rate of :math:`1 \times 10^{-8}`
per base per generation::

    >>> tree_sequence = msprime.simulate(
    ...    sample_size=6, Ne=1000, length=50e3, mutation_rate=1e-8, random_seed=30)
    >>> tree = tree_sequence.first()
    >>> for site in tree.sites():
    ...     for mutation in site.mutations:
    ...         print("Mutation @ position {:.2f} over node {}".format(
    ...             site.position, mutation.node))
    Mutation @ position 1556.54 over node 9
    Mutation @ position 4485.17 over node 6
    Mutation @ position 9788.56 over node 6
    Mutation @ position 11759.03 over node 6
    Mutation @ position 11949.32 over node 6
    Mutation @ position 14321.77 over node 9
    Mutation @ position 31454.99 over node 6
    Mutation @ position 45125.69 over node 9
    Mutation @ position 49709.68 over node 6

    >>> print(tree.draw(format="unicode"))
        10
     ┏━━┻━━┓
     ┃     9
     ┃   ┏━┻━┓
     ┃   ┃   8
     ┃   ┃  ┏┻┓
     ┃   7  ┃ ┃
     ┃  ┏┻┓ ┃ ┃
     6  ┃ ┃ ┃ ┃
    ┏┻┓ ┃ ┃ ┃ ┃
    0 4 2 5 1 3


********
Variants
********

We are often interesting in accessing the sequence data that results from
simulations directly. The most efficient way to do this is by using
the :meth:`.TreeSequence.variants` method, which returns an iterator
over all the :class:`.Variant` objects arising from the trees and mutations.
Each variant contains a reference to the site object, as well as the
alleles and the observed sequences for each sample in the ``genotypes``
field::

    >>> tree_sequence = msprime.simulate(
    ...     sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
    ...     mutation_rate=2e-8, random_seed=10)
    >>> for variant in tree_sequence.variants():
    ...     print(
    ...         variant.site.id, variant.site.position,
    ...         variant.alleles, variant.genotypes, sep="\t")
    0       2432.768327416852       ('0', '1')      [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
    1       2577.6939414924095      ('0', '1')      [1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1]
    2       2844.682702049562       ('0', '1')      [0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
    3       4784.266628557816       ('0', '1')      [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]

In this example we simulate some data and then print out the observed
sequences. We loop through each variant and print out the observed state of
each sample as an array of zeros and ones, along with the index and position
of the corresponding mutation.  In this example, the
alleles are always ``'0'`` (the ancestral state) and ``'1'``
(the derived state), because we are simulating with the infinite sites mutation
model, in which each mutation occurs at a unique position in the genome.
More complex models are possible, however.

This way of working with the sequence data is quite efficient because we
do not need to keep the entire genotype matrix in memory at once. However, if
we do want the full genotype matrix it is simple to obtain::

    >>> A = tree_sequence.genotype_matrix()
    >>> A
    array([[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
           [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=uint8)

In this example, we run the same simulation but this time
store the entire variant matrix in a two-dimensional numpy array.
This is useful for integrating with tools such as
`scikit allel <https://scikit-allel.readthedocs.io/en/latest/>`_.

******************
Historical samples
******************

Simulating coalescent histories in which some of the samples are not
from the present time is straightforward in ``msprime``.
By using the ``samples`` argument to :meth:`msprime.simulate`
we can specify the location and time at which all samples are made.

.. code-block:: python

    def historical_samples_example():
        samples = [
            msprime.Sample(population=0, time=0),
            msprime.Sample(0, 0),  # Or, we can use positional arguments.
            msprime.Sample(0, 1.0),
            msprime.Sample(0, 1.0)
        ]
        tree_seq = msprime.simulate(samples=samples)
        tree = tree_seq.first()
        for u in tree.nodes():
            print(u, tree.parent(u), tree.time(u), sep="\t")
        print(tree.draw(format="unicode"))

In this example we create four samples, two taken at the present time
and two taken 1.0 generations in the past, as might represent one modern
and one ancient diploid individual. There are a number of
different ways in which we can describe the samples using the
``msprime.Sample`` object (samples can be provided as plain tuples also
if more convenient). Running this example, we get::


    >>> historical_samples_example()
    6    -1    2.8240255501413247
    4    6    0.0864109319103291
    0    4    0.0
    1    4    0.0
    5    6    1.9249243960710336
    2    5    1.0
    3    5    1.0
       6
     ┏━┻━┓
     ┃   5
     ┃  ┏┻┓
     ┃  2 3
     ┃
     4
    ┏┻┓
    0 1


Because nodes ``0`` and ``1`` were sampled at time 0, their times in the tree
are both 0. Nodes ``2`` and ``3`` were sampled at time 1.0, and so their times are recorded
as 1.0 in the tree.

***********
Replication
***********

A common task for coalescent simulations is to check the accuracy of analytical
approximations to statistics of interest. To do this, we require many independent
replicates of a given simulation. ``msprime`` provides a simple and efficient
API for replication: by providing the ``num_replicates`` argument to the
:func:`.simulate` function, we can iterate over the replicates
in a straightforward manner. Here is an example where we compare the
analytical results for the number of segregating sites with simulations:

.. code-block:: python

    import msprime
    import numpy as np

    def segregating_sites_example(n, theta, num_replicates):
        S = np.zeros(num_replicates)
        replicates = msprime.simulate(
            Ne=0.5,
            sample_size=n,
            mutation_rate=theta / 2,
            num_replicates=num_replicates)
        for j, tree_sequence in enumerate(replicates):
            S[j] = tree_sequence.num_sites
        # Now, calculate the analytical predictions
        S_mean_a = np.sum(1 / np.arange(1, n)) * theta
        S_var_a = (
            theta * np.sum(1 / np.arange(1, n)) +
            theta**2 * np.sum(1 / np.arange(1, n)**2))
        print("              mean              variance")
        print("Observed      {}\t\t{}".format(np.mean(S), np.var(S)))
        print("Analytical    {:.5f}\t\t{:.5f}".format(S_mean_a, S_var_a))

Running this code, we get::

    >>> segregating_sites_example(10, 5, 100000)
              mean              variance
    Observed      14.17893          53.0746740551
    Analytical    14.14484          52.63903


Note that in this example we set :math:`N_e = 0.5` and
the mutation rate to :math:`\theta / 2` when calling :func:`.simulate`.
This works because ``msprime`` simulates Kingman's coalescent,
for which :math:`N_e` is only a time scaling;
since :math:`N_e` is the diploid effective population size,
setting :math:`N_e = 0.5` means that the mean time for two samples to coalesce
is equal to one time unit in the resulting trees.
This is helpful for converting the diploid per-generation time units
of msprime into the haploid coalescent units used in many
theoretical results. However, it is important to note that conventions
vary widely, and great care is needed with such factor-of-two
rescalings.

********************
Population structure
********************


Population structure in ``msprime`` closely follows the model used in the
``ms`` simulator: we have :math:`N` subpopulations with an :math:`N\times N`
matrix describing the migration rates between these subpopulations. The
sample sizes, population sizes and growth rates of all subpopulations
can be specified independently. Migration rates are specified using
a migration matrix. Unlike ``ms`` however, all times and rates are specified
in generations and all populations sizes are absolute (that is, not
multiples of :math:`N_e`).

In the following example, we calculate the mean coalescence time for
a pair of lineages sampled in different subpopulations in a symmetric island
model, and compare this with the analytical expectation.

.. code-block:: python

    import msprime
    import numpy as np

    def migration_example(num_replicates=10**4):
        # M is the overall symmetric migration rate, and d is the number
        # of subpopulations.
        M = 0.2
        d = 3
        m = M / (2 * (d - 1))
        # Allocate the initial sample. Because we are interested in the
        # between-subpopulation coalescence times, we choose one sample each
        # from the first two subpopulations.
        population_configurations = [
            msprime.PopulationConfiguration(sample_size=1),
            msprime.PopulationConfiguration(sample_size=1),
            msprime.PopulationConfiguration(sample_size=0)]
        # Now we set up the migration matrix. Since this is a symmetric
        # island model, we have the same rate of migration between all
        # pairs of subpopulations. Diagonal elements must be zero.
        migration_matrix = [
            [0, m, m],
            [m, 0, m],
            [m, m, 0]]
        # We pass these values to the simulate function, and ask it
        # to run the required number of replicates.
        replicates = msprime.simulate(Ne=0.5,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            num_replicates=num_replicates)
        # And then iterate over these replicates
        T = np.zeros(num_replicates)
        for i, tree_sequence in enumerate(replicates):
            tree = tree_sequence.first()
            T[i] = tree.time(tree.root) / 4
        # Finally, calculate the analytical expectation and print
        # out the results
        analytical = d / 4 + (d - 1) / (4 * M)
        print("Observed  =", np.mean(T))
        print("Predicted =", analytical)

Again, we set :math:`N_e = 0.5` to agree with convention in theoretical results,
where usually one coalescent time unit is, in generations, the effective number of *haploid* individuals.
Running this example we get::


    >>> migration_example()
    Observed  = 3.254904176088153
    Predicted = 3.25


**********
Demography
**********

Msprime provides a flexible and simple way to model past demographic events
in arbitrary combinations. Here is an example describing the
`Gutenkunst et al. <http://dx.doi.org/10.1371/journal.pgen.1000695>`_
out-of-Africa model. See
`Figure 2B <http://dx.doi.org/10.1371/journal.pgen.1000695.g002>`_
for a schematic of this model, and
`Table 1 <http://dx.doi.org/10.1371/journal.pgen.1000695.t001>`_ for
the values used.
Coalescent simulation moves from the present back into the past,
so times are in units of generations *ago*, and we build the model
with most recent events first.

.. todo:: Add a diagram of the model for convenience.

.. code-block:: python

    import math
    def out_of_africa():
        # First we set out the maximum likelihood values of the various parameters
        # given in Table 1.
        N_A = 7300
        N_B = 2100
        N_AF = 12300
        N_EU0 = 1000
        N_AS0 = 510
        # Times are provided in years, so we convert into generations.
        generation_time = 25
        T_AF = 220e3 / generation_time
        T_B = 140e3 / generation_time
        T_EU_AS = 21.2e3 / generation_time
        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        r_EU = 0.004
        r_AS = 0.0055
        N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
        N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
        # Migration rates during the various epochs.
        m_AF_B = 25e-5
        m_AF_EU = 3e-5
        m_AF_AS = 1.9e-5
        m_EU_AS = 9.6e-5
        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
        # initially.
        population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=0, initial_size=N_AF),
            msprime.PopulationConfiguration(
                sample_size=1, initial_size=N_EU, growth_rate=r_EU),
            msprime.PopulationConfiguration(
                sample_size=1, initial_size=N_AS, growth_rate=r_AS)
        ]
        migration_matrix = [
            [      0, m_AF_EU, m_AF_AS],
            [m_AF_EU,       0, m_EU_AS],
            [m_AF_AS, m_EU_AS,       0],
        ]
        demographic_events = [
            # CEU and CHB merge into B with rate changes at T_EU_AS
            msprime.MassMigration(
                time=T_EU_AS, source=2, destination=1, proportion=1.0),
            msprime.MigrationRateChange(time=T_EU_AS, rate=0),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
            # Population B merges into YRI at T_B
            msprime.MassMigration(
                time=T_B, source=1, destination=0, proportion=1.0),
            # Size changes to N_A at T_AF
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        dd.print_history()


The :class:`.DemographyDebugger` provides a method to debug the history that
you have described so that you can be sure that the migration rates, population
sizes and growth rates are all as you intend during each epoch::

    =============================
    Epoch: 0 -- 848.0 generations
    =============================
         start     end      growth_rate |     0        1        2
       -------- --------       -------- | -------- -------- --------
    0 |1.23e+04 1.23e+04              0 |     0      3e-05   1.9e-05
    1 |2.97e+04   1e+03           0.004 |   3e-05      0     9.6e-05
    2 |5.41e+04    510           0.0055 |  1.9e-05  9.6e-05     0

    Events @ generation 848.0
       - Mass migration: lineages move from 2 to 1 with probability 1.0
       - Migration rate change to 0 everywhere
       - Migration rate change for (0, 1) to 0.00025
       - Migration rate change for (1, 0) to 0.00025
       - Population parameter change for 1: initial_size -> 2100 growth_rate -> 0
    ==================================
    Epoch: 848.0 -- 5600.0 generations
    ==================================
         start     end      growth_rate |     0        1        2
       -------- --------       -------- | -------- -------- --------
    0 |1.23e+04 1.23e+04              0 |     0     0.00025     0
    1 | 2.1e+03  2.1e+03              0 |  0.00025     0        0
    2 |   510   2.27e-09         0.0055 |     0        0        0

    Events @ generation 5600.0
       - Mass migration: lineages move from 1 to 0 with probability 1.0
    ===================================
    Epoch: 5600.0 -- 8800.0 generations
    ===================================
         start     end      growth_rate |     0        1        2
       -------- --------       -------- | -------- -------- --------
    0 |1.23e+04 1.23e+04              0 |     0     0.00025     0
    1 | 2.1e+03  2.1e+03              0 |  0.00025     0        0
    2 |2.27e-09 5.17e-17         0.0055 |     0        0        0

    Events @ generation 8800.0
       - Population parameter change for 0: initial_size -> 7300
    ================================
    Epoch: 8800.0 -- inf generations
    ================================
         start     end      growth_rate |     0        1        2
       -------- --------       -------- | -------- -------- --------
    0 | 7.3e+03  7.3e+03              0 |     0     0.00025     0
    1 | 2.1e+03  2.1e+03              0 |  0.00025     0        0
    2 |5.17e-17     0            0.0055 |     0        0        0

.. warning:: The output of the :meth:`.DemographyDebugger.print_history` method
    is intended only for debugging purposes, and is not meant to be machine
    readable. The format is also preliminary; if there is other information
    that you think would be useful, please `open an issue on GitHub
    <https://github.com/tskit-dev/msprime/issues>`_

Once you are satisfied that the demographic history that you have built
is correct, it can then be simulated by calling the :func:`.simulate`
function.

******************
Recombination maps
******************

The ``msprime`` API allows us to quickly and easily simulate data from an
arbitrary recombination map. In this example we read a recombination
map for human chromosome 22, and simulate a single replicate. After
the simulation is completed, we plot histograms of the recombination
rates and the simulated breakpoints. These show that density of
breakpoints follows the recombination rate closely.

.. code-block:: python

    import numpy as np
    import scipy.stats
    import matplotlib.pyplot as pyplot

    def variable_recomb_example():
        infile = "hapmap/genetic_map_GRCh37_chr22.txt"
        # Read in the recombination map using the read_hapmap method,
        recomb_map = msprime.RecombinationMap.read_hapmap(infile)

        # Now we get the positions and rates from the recombination
        # map and plot these using 500 bins.
        positions = np.array(recomb_map.get_positions()[1:])
        rates = np.array(recomb_map.get_rates()[1:])
        num_bins = 500
        v, bin_edges, _ = scipy.stats.binned_statistic(
            positions, rates, bins=num_bins)
        x = bin_edges[:-1][np.logical_not(np.isnan(v))]
        y = v[np.logical_not(np.isnan(v))]
        fig, ax1 = pyplot.subplots(figsize=(16, 6))
        ax1.plot(x, y, color="blue")
        ax1.set_ylabel("Recombination rate")
        ax1.set_xlabel("Chromosome position")

        # Now we run the simulation for this map. We simulate
        # 50 diploids (100 sampled genomes) in a population with Ne=10^4.
        tree_sequence = msprime.simulate(
            sample_size=100,
            Ne=10**4,
            recombination_map=recomb_map)
        # Now plot the density of breakpoints along the chromosome
        breakpoints = np.array(list(tree_sequence.breakpoints()))
        ax2 = ax1.twinx()
        v, bin_edges = np.histogram(breakpoints, num_bins, density=True)
        ax2.plot(bin_edges[:-1], v, color="green")
        ax2.set_ylabel("Breakpoint density")
        ax2.set_xlim(1.5e7, 5.3e7)
        fig.savefig("hapmap_chr22.svg")


.. image:: _static/hapmap_chr22.svg
   :width: 800px
   :alt: Density of breakpoints along the chromosome.

**************
Calculating LD
**************

The ``msprime`` API provides methods to efficiently calculate
population genetics statistics. For example, the :class:`.LdCalculator`
class allows us to compute pairwise `linkage disequilibrium
<https://en.wikipedia.org/wiki/Linkage_disequilibrium>`_ coefficients.
Here we use the :meth:`.get_r2_matrix` method to easily make an
LD plot using `matplotlib <http://matplotlib.org/>`_. (Thanks to
the excellent `scikit-allel
<http://scikit-allel.readthedocs.io/en/latest/index.html>`_
for the basic `plotting code
<http://scikit-allel.readthedocs.io/en/latest/_modules/allel/stats/ld.html#plot_pairwise_ld>`_
used here.)

.. code-block:: python

    import msprime
    import matplotlib.pyplot as pyplot

    def ld_matrix_example():
        ts = msprime.simulate(100, recombination_rate=10, mutation_rate=20,
                random_seed=1)
        ld_calc = msprime.LdCalculator(ts)
        A = ld_calc.r2_matrix()
        # Now plot this matrix.
        x = A.shape[0] / pyplot.rcParams['figure.dpi']
        x = max(x, pyplot.rcParams['figure.figsize'][0])
        fig, ax = pyplot.subplots(figsize=(x, x))
        fig.tight_layout(pad=0)
        im = ax.imshow(A, interpolation="none", vmin=0, vmax=1, cmap="Blues")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in 'top', 'bottom', 'left', 'right':
            ax.spines[s].set_visible(False)
        pyplot.gcf().colorbar(im, shrink=.5, pad=0)
        pyplot.savefig("ld.svg")


.. image:: _static/ld.svg
   :width: 800px
   :alt: An example LD matrix plot.

.. _sec_tutorial_threads:

********************
Working with threads
********************

When performing large calculations it's often useful to split the
work over multiple processes or threads. The msprime API can
be used without issues across multiple processes, and the Python
:mod:`multiprocessing` module often provides a very effective way to
work with many replicate simulations in parallel.

When we wish to work with a single very large dataset, however, threads can
offer better resource usage because of the shared memory space. The Python
:mod:`threading` library gives a very simple interface to lightweight CPU
threads and allows us to perform several CPU intensive tasks in parallel. The
``msprime`` API is designed to allow multiple threads to work in parallel when
CPU intensive tasks are being undertaken.

.. note:: In the CPython implementation the `Global Interpreter Lock
   <https://wiki.python.org/moin/GlobalInterpreterLock>`_ ensures that
   only one thread executes Python bytecode at one time. This means that
   Python code does not parallelise well across threads, but avoids a large
   number of nasty pitfalls associated with multiple threads updating
   data structures in parallel. Native C extensions like ``numpy`` and ``msprime``
   release the GIL while expensive tasks are being performed, therefore
   allowing these calculations to proceed in parallel.

In the following example we wish to find all mutations that are in approximate
LD (:math:`r^2 \geq 0.5`) with a given set of mutations. We parallelise this
by splitting the input array between a number of threads, and use the
:meth:`.LdCalculator.r2_array` method to compute the :math:`r^2` value
both up and downstream of each focal mutation, filter out those that
exceed our threshold, and store the results in a dictionary. We also
use the very cool `tqdm <https://pypi.python.org/pypi/tqdm>`_ module to give us a
progress bar on this computation.

.. code-block:: python

    import threading
    import numpy as np
    import tqdm
    import msprime

    def find_ld_sites(
            tree_sequence, focal_mutations, max_distance=1e6, r2_threshold=0.5,
            num_threads=8):
        results = {}
        progress_bar = tqdm.tqdm(total=len(focal_mutations))
        num_threads = min(num_threads, len(focal_mutations))

        def thread_worker(thread_index):
            ld_calc = msprime.LdCalculator(tree_sequence)
            chunk_size = int(math.ceil(len(focal_mutations) / num_threads))
            start = thread_index * chunk_size
            for focal_mutation in focal_mutations[start: start + chunk_size]:
                a = ld_calc.r2_array(
                    focal_mutation, max_distance=max_distance,
                    direction=msprime.REVERSE)
                rev_indexes = focal_mutation - np.nonzero(a >= r2_threshold)[0] - 1
                a = ld_calc.r2_array(
                    focal_mutation, max_distance=max_distance,
                    direction=msprime.FORWARD)
                fwd_indexes = focal_mutation + np.nonzero(a >= r2_threshold)[0] + 1
                indexes = np.concatenate((rev_indexes[::-1], fwd_indexes))
                results[focal_mutation] = indexes
                progress_bar.update()

        threads = [
            threading.Thread(target=thread_worker, args=(j,))
            for j in range(num_threads)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        progress_bar.close()
        return results

    def threads_example():
        ts = msprime.simulate(
            sample_size=1000, Ne=1e4, length=1e7, recombination_rate=2e-8,
            mutation_rate=2e-8)
        counts = np.zeros(ts.num_sites)
        for tree in ts.trees():
            for site in tree.sites():
                assert len(site.mutations) == 1
                mutation = site.mutations[0]
                counts[site.id] = tree.num_samples(mutation.node)
        doubletons = np.nonzero(counts == 2)[0]
        results = find_ld_sites(ts, doubletons, num_threads=8)
        print(
            "Found LD sites for", len(results), "doubleton sites out of",
            ts.num_sites)

In this example, we first simulate 1000 samples of 10 megabases and find all
doubleton mutations in the resulting tree sequence. We then call the
``find_ld_sites()`` function to find all mutations that are within 1 megabase
of these doubletons and have an :math:`r^2` statistic of greater than 0.5.

The ``find_ld_sites()`` function performs these calculations in parallel using
8 threads. The real work is done in the nested ``thread_worker()`` function,
which is called once by each thread. In the thread worker, we first allocate an
instance of the :class:`.LdCalculator` class. (It is **critically important**
that each thread has its own instance of :class:`.LdCalculator`, as the threads
will not work efficiently otherwise.) After this, each thread works out the
slice of the input array that it is responsible for, and then iterates over
each focal mutation in turn. After the :math:`r^2` values have been calculated,
we then find the indexes of the mutations corresponding to values greater than
0.5 using :func:`numpy.nonzero`. Finally, the thread stores the resulting array
of mutation indexes in the ``results`` dictionary, and moves on to the next
focal mutation.


Running this example we get::

    >>> threads_example()
    100%|████████████████████████████████████████████████| 4045/4045 [00:09<00:00, 440.29it/s]
    Found LD sites for 4045 doubleton mutations out of 60100

**********************
Editing tree sequences
**********************

Sometimes we wish to make some minor modifications to a tree sequence that has
been generated by a simulation. However, tree sequence objects are **immutable**
and so we cannot edit a them in place. To modify a tree sequence, we need to
extract the underlying :ref:`tables <sec_table_definitions>` of information, edit these tables,
and then create a new tree sequence from them.
These tables succinctly store everything we need to know
about a tree sequence, and can be manipulated using the :ref:`sec_tables_api`.
In the following example, we use this approach
to remove all singleton sites from a given tree sequence.

.. code-block:: python

    def strip_singletons(ts):
        tables = ts.dump_tables()
        tables.sites.clear()
        tables.mutations.clear()
        for tree in ts.trees():
            for site in tree.sites():
                assert len(site.mutations) == 1  # Only supports infinite sites muts.
                mut = site.mutations[0]
                if tree.num_samples(mut.node) > 1:
                    site_id = tables.sites.add_row(
                        position=site.position,
                        ancestral_state=site.ancestral_state)
                    tables.mutations.add_row(
                        site=site_id, node=mut.node, derived_state=mut.derived_state)
        return tables.tree_sequence()


This function takes a tree sequence containing some infinite sites mutations as
input, and returns a copy in which all singleton sites have been removed.
The approach is very simple: we get a copy of the underlying
table data in a :class:`.TableCollection` object, and first clear the
site and mutation tables. We then consider each site in turn,
and if the allele frequency of
the mutation is greater than one, we add the site and mutation to our
output tables using :meth:`.SiteTable.add_row` and :meth:`.MutationTable.add_row`.
(In this case we consider only simple infinite sites mutations,
where we cannot have back or recurrent mutations. These would require a slightly
more involved approach where we keep a map of mutation IDs so that
mutation ``parent`` values could be computed. We have also omitted the
site and mutation metadata in the interest of simplicity.)

After considering each site, we then create a new tree sequence using
the :meth:`.TableCollection.tree_sequence` method on our updated tables.
Using this function then, we get::

    >>> ts = msprime.simulate(10, mutation_rate=10)
    >>> ts.num_sites
    50
    >>> ts_new = strip_singletons(ts)
    >>> ts_new.num_sites
    44
    >>>

Thus, we have removed 6 singleton sites from the tree sequence.

.. todo::

    Add another example here where we use the array oriented API to edit
    the nodes and edges of a tree sequence. Perhaps decapitating would be a
    good example?

*******************
Working with Tables
*******************


Tables provide a convenient method for viewing, importing and exporting tree
sequences, and are closely tied to the underlying data structures.
There are eight tables that together define a tree sequence,
although some may be empty,
and together they form a :class:`TableCollection`.
The tables are defined in :ref:`Table Definitions <sec_table_definitions>`,
and the :ref:`Tables API <sec_tables_api>` section describes how to work with them.
Here we make some general remarks about what you can, and cannot do with them.


``msprime`` provides direct access to the columns of each table as
``numpy`` arrays: for instance, if ``n`` is a ``NodeTable``, then ``n.time``
will return an array containing the birth times of the individuals whose genomes
are represented by the nodes in the table.
*However*, it is important to note that this is *not* a shallow copy:
modifying ``n.time`` will *not* change the node table ``n``.  This may change in
the future, but currently there are three ways to modify tables: ``.add_row()``,
``.set_columns()``, and ``.append_columns()``
(and also ``.clear()``, which empties the table).

For example, a node table could be constructed using ``.add_row()`` as
follows::

    n = msprime.NodeTable()
    sv = [True, True, True, False, False, False, False]
    tv = [0.0, 0.0, 0.0, 0.4, 0.5, 0.7, 1.0]
    pv = [0, 0, 0, 0, 0, 0, 0]
    for s, t, p in zip(sv, tv, pv):
        n.add_row(flags=s, population=p, time=t)


obtaining::

    >>> print(n)
    id    flags    population    individual    time    metadata
    0    1    0    -1    0.0
    1    1    0    -1    0.0
    2    1    0    -1    0.0
    3    0    0    -1    0.4
    4    0    0    -1    0.5
    5    0    0    -1    0.7
    6    0    0    -1    1.0


The ``.add_row()`` method is natural (and should be reasonably efficient) if
new records appear one-by-one. In the example above it would have been more
natural to use ``.set_columns()`` - equivalently::

    n = msprime.NodeTable()
    n.set_columns(flags=sv, population=pv, time=tv)

Since columns cannot be modified directly as properties of the tables,
they must be extracted, modified, then replaced.
For example, here we add 1.4 to every ``time`` except the first
in the node table constructed above (using ``numpy`` indexing)::

    tn = n.time
    tn[1:] = tn[1:] + 1.4
    n.set_columns(flags=n.flags, population=n.population, time=tn)

The result is::

    >>> print(n)
    id    flags    population    individual    time    metadata
    0    1    0    -1    0.0
    1    1    0    -1    1.4
    2    1    0    -1    1.4
    3    0    0    -1    1.8
    4    0    0    -1    1.9
    5    0    0    -1    2.1
    6    0    0    -1    2.4


*****************************
Overview of the Tables Format
*****************************

The :ref:`Table Definitions <sec_table_definitions>` section gives a precise
definition of how a tree sequence is stored in a collection of tables.
Here we give an overview. Consider the following sequence of trees::

    time ago
    --------
       1.0         6
                 ┏━┻━━┓
                 ┃    ┃
       0.7       ┃    ╋                     5
                 ┃    ┃                   ┏━┻━┓
       0.5       ┃    4         4         ┃   4
                 ┃  ┏━┻━┓     ┏━┻━┓       ┃ ┏━┻━┓
                 ┃  ┃   ┃     ┃   ╋       ┃ ┃   ┃
       0.4       ┃  ┃   ┃     ┃   3       ┃ ┃   ┃
                 ┃  ┃   ┃     ┃ ┏━┻━┓     ┃ ┃   ┃
                 ┃  ┃   ┃     ┃ ┃   ╋     ┃ ┃   ┃
       0.0       0  1   2     1 0   2     0 1   2

    position 0.0          0.2         0.8         1.0

Ancestral recombination events have produced three different trees
that relate the three sampled genomes ``0``, ``1``, and ``2`` to each other
along the chromosome of length 1.0.

Each node in each of the above trees represents a particular ancestral genome
(a *haploid* genome; diploid individuals would be represented by two nodes).
We record when each of nodes lived in a :class:`NodeTable`::

    NodeTable:

    id      flags    population   time
    0       1        0            0
    1       1        0            0
    2       1        0            0
    3       0        0            0.4
    4       0        0            0.5
    5       0        0            0.7
    6       0        0            1.0

Importantly, the first column, ``id``, is not actually recorded, and is
only shown when printing out node tables (as here) for convenience.
The second column, ``flags`` records a ``1`` for the individuals that are *samples*,
i.e., whose entire genealogical history is recorded by these trees.
(Note that the trees above record that node 3 inherited from node 4
on the middle portion of the genome, but not on the ends.)

We next need to record each tree's edges. Since some edges are present
in more than one tree (e.g., node 1 inherits from node 4 across
the entire sequence), we record in the :class:`EdgeTable` each edge
and the genomic region for which it appears in the trees::


    EdgeTable:

    left    right   parent  children
    0.2     0.8     3       0
    0.2     0.8     3       2
    0.0     1.0     4       1
    0.0     0.2     4       2
    0.8     1.0     4       2
    0.2     0.8     4       3
    0.8     1.0     5       0
    0.8     1.0     5       4
    0.0     0.2     6       0
    0.0     0.2     6       4

Since node 3 is most recent, the edge that says that nodes 0 and 2 inherit
from node 3 on the interval between 0.2 and 0.8 comes first.  Next are the
edges from node 4: there are six of these, two for each of the three genomic
intervals over which node 4 is ancestor to a distinct set of nodes.  At this
point, we know the full tree on the middle interval.  Finally, edges
specifying the common ancestor of 0 and 4 on the remaining intervals (parents 6
and 5 respectively) allow us to construct all trees across the entire interval.

There are three mutations in the depiction above,
marked by ``╋``: one above node ``4`` on the first tree,
and the other two above nodes ``2`` and ``3`` on the second tree.
Suppose that the first mutation occurs at position 0.1 and the mutations in the
second tree both occurred at the same position, at 0.5 (with a back mutation).
To record the inheritance patterns of these, we need only record
the positions on the genome at which they occurred,
and on which edge (equivalently, above which node) they occurred.
The positions are recorded in the :class:`SiteTable`::

    SiteTable:

    id    position    ancestral_state
    0    0.1         0
    1    0.5         0

As with node tables, the ``id`` column is **not** actually recorded, but is
implied by the position in the table.  The results of the
actual mutations are then recorded::

    MutationTable:

    site    node    derived_state
    0        4        1
    1        3        1
    1        2        0

This would then result in the following (two-locus) haplotypes for the three
samples::

    sample  haplotype
    ------  ---------
    0       01
    1       10
    2       10


To create these tables, and the corresponding tree sequence, we would
create a :class:`TableCollection`, and then use its
:meth:`TableCollection.tree_sequence` method::

    tables = msprime.TableCollection(sequence_length=1.0)

    # Nodes
    sv = [True, True, True, False, False, False, False]
    tv = [0.0, 0.0, 0.0, 0.4, 0.5, 0.7, 1.0]

    for is_sample, t in zip(sv, tv):
     flags = msprime.NODE_IS_SAMPLE if is_sample else 0
     tables.nodes.add_row(flags=flags, time=t)

    # Edges
    lv = [0.2, 0.2, 0.0, 0.2, 0.8, 0.0, 0.8, 0.2, 0.8, 0.8, 0.0, 0.0]
    rv = [0.8, 0.8, 0.2, 0.8, 1.0, 0.2, 1.0, 0.8, 1.0, 1.0, 0.2, 0.2]
    pv = [3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 6, 6]
    cv = [0, 2, 1, 1, 1, 2, 2, 3, 0, 4, 0, 4]

    for l, r, p, c in zip(lv, rv, pv, cv):
        tables.edges.add_row(left=l, right=r, parent=p, child=c)

    # Sites
    for p, a in zip([0.1, 0.5], ['0', '0']):
        tables.sites.add_row(position=p, ancestral_state=a)

    # Mutations
    for s, n, d in zip([0, 1, 1], [4, 3, 2], ['1', '1', '0']):
        tables.mutations.add_row(site=s, node=n, derived_state=d)

We can then finally obtain the tree sequence::

    ts = tables.tree_sequence()
    for t in ts.trees():
      print(t.draw(format='unicode'))


.. _sec_tutorial_simulate_from:

*******************************
Completing forwards simulations
*******************************

The ``msprime`` simulator generates tree sequences using the backwards in
time coalescent model. But it is also possible to output tree sequences
from `forwards-time <https://www.biorxiv.org/content/early/2018/01/16/248500>`_
simulators such as `SLiM <https://messerlab.org/slim/>`_.
There are many advantages to using forward-time simulators, but they
are usually quite slow compared to similar coalescent simulations. In this
section we show how to combine the best of both approaches by simulating
the recent past using a forwards-time simulator and then complete the
simulation of the ancient past using ``msprime``. (We sometimes refer to this
"recapitation", as we can think of it as adding a "head" onto a tree sequence.)

First, we define a simple Wright-Fisher simulator which returns a tree sequence
with the properties that we require (please see the :ref:`API <sec_api_simulate_from>`
section for a formal description of these properties):

.. code-block:: python

    import random
    import numpy as np

    def wright_fisher(N, T, L=100, random_seed=None):
        """
        Simulate a Wright-Fisher population of N haploid individuals with L
        discrete loci for T generations. Based on Algorithm W from
        https://www.biorxiv.org/content/biorxiv/early/2018/01/16/248500.full.pdf
        """
        random.seed(random_seed)
        tables = msprime.TableCollection(L)
        P = np.arange(N, dtype=int)
        # Mark the initial generation as samples so that we remember these nodes.
        for j in range(N):
            tables.nodes.add_row(time=T, flags=msprime.NODE_IS_SAMPLE)
        t = T
        while t > 0:
            t -= 1
            Pp = P.copy()
            for j in range(N):
                u = tables.nodes.add_row(time=t, flags=0)
                Pp[j] = u
                a = random.randint(0, N - 1)
                b = random.randint(0, N - 1)
                x = random.randint(1, L - 1)
                tables.edges.add_row(0, x, P[a], u)
                tables.edges.add_row(x, L, P[b], u)
            P = Pp

        # Now do some table manipulations to ensure that the tree sequence
        # that we output has the form that msprime needs to finish the
        # simulation. Much of the complexity here is caused by the tables API
        # not allowing direct access to memory, which will change soon.

        # Mark the extant population as samples also
        flags = tables.nodes.flags
        flags[P] = msprime.NODE_IS_SAMPLE
        tables.nodes.set_columns(flags=flags, time=tables.nodes.time)
        tables.sort()
        # Simplify with respect to the current generation, but ensuring we keep the
        # ancient nodes from the initial population.
        tables.simplify()
        # Unmark the initial generation as samples
        flags = tables.nodes.flags
        time = tables.nodes.time
        flags[:] = 0
        flags[time == 0] = msprime.NODE_IS_SAMPLE
        # The final tables must also have at least one population which
        # the samples are assigned to
        tables.populations.add_row()
        tables.nodes.set_columns(
            flags=flags, time=time,
            population=np.zeros_like(tables.nodes.population))
        return tables.tree_sequence()


We then run a tiny forward simulation of 10 two-locus individuals
for 5 generations, and print out the resulting trees:

.. code-block:: python

    num_loci = 2
    N = 10
    wf_ts = wright_fisher(N, 5, L=num_loci, random_seed=3)
    for tree in wf_ts.trees():
        print("interval = ", tree.interval)
        print(tree.draw(format="unicode"))


We get::

    interval =  (0.0, 1.0)
           0                 7
           ┃                 ┃
           25                ┃
      ┏━━━━┻━━━━┓            ┃
      23        24           ┃
    ┏━┻━┓    ┏━━╋━━━┓        ┃
    ┃   21   ┃  ┃   22       20
    ┃  ┏┻━┓  ┃  ┃  ┏┻━┓   ┏━━╋━━┓
    10 14 19 11 18 15 17  12 13 16

    interval =  (1.0, 2.0)
            0          8    4     7
            ┃          ┃   ┏┻━┓   ┃
            21         ┃   ┃  ┃   ┃
    ┏━━┳━━┳━┻┳━━┳━━┓   ┃   ┃  ┃   ┃
    14 19 10 13 16 18  11  15 17  12

Because our Wright Fisher simulation ran for only 5 generations, there has not
been enough time for the trees to fully coalesce. Therefore, instead of having
one root, the trees have several --- the first tree has 2 and the second 4.
Nodes 0 to 9 in this simulation represent the initial population of the
simulation, and so we can see that all samples in the first tree trace back
to one of two individuals from the initial generation.
These unary branches joining samples and coalesced subtrees to the nodes
in the initial generation are essential as they allow use to correctly
assemble the various fragments of ancestral material into chromosomes
when creating the initial conditions for the coalescent simulation.
(Please see the :ref:`API <sec_api_simulate_from>` section for more details on the
required properties of input tree sequences.)

The process of completing this tree sequence using a coalescent simulation
begins by first examining the root segments on the input trees. We get the
following segments::

    [(0, 2, 0), (0, 2, 7), (1, 2, 8), (1, 2, 4)]

where each segment is a ``(left, right, node)`` tuple. As nodes 0 and 7 are
present in both trees, they have segments spanning both loci. Nodes 8 and 4 are
present only in the second tree, and so they have ancestral segments only for
the second locus. Note that this means that we do *not* simulate the ancestry
of the entire initial generation of the simulation, but rather the exact
minimum that we need in order to complete the ancestry of the current
generation. For instance, root ``8`` has not coalesced over the interval from
``1.0`` to ``2.0``, while root ``0`` has not coalesced over the entire segment
from ``0.0`` to ``2.0``.

We run the coalescent simulation to complete this tree sequence using the
``from_ts`` argument to :func:`.simulate`. Because we have simulated a
two locus system with a recombination rate of ``1 / num_loci`` per generation
in the Wright-Fisher model, we want to use the same system in the coalescent simulation.
To do this we create recombination map using the
:meth:`.RecombinationMap.uniform_map` class method to easily create a
discrete map with the required number of loci.
(Please see the :ref:`API <sec_api_simulate_from>` section for more details on the
restrictions on recombination maps when completing an existing simulation.)
We also use a ``Ne`` value of ``N / 2``
since the Wright-Fisher simulation was haploid and ``msprime`` is diploid.

.. code-block:: python

    recomb_map = msprime.RecombinationMap.uniform_map(num_loci, 1 / num_loci, num_loci)
    coalesced_ts = msprime.simulate(
        Ne=N / 2, from_ts=wf_ts, recombination_map=recomb_map, random_seed=5)



After running this simulation we get the following trees::

    interval =  (0.0, 1.0)
                    26
           ┏━━━━━━━━┻━━━━━━━┓
           0                7
           ┃                ┃
           25               ┃
      ┏━━━━┻━━━━┓           ┃
      23        24          ┃
    ┏━┻━┓    ┏━━╋━━━┓       ┃
    ┃   21   ┃  ┃   22      20
    ┃  ┏┻━┓  ┃  ┃  ┏┻━┓  ┏━━╋━━┓
    10 14 19 11 18 15 17 12 13 16

    interval =  (1.0, 2.0)
                      28
                 ┏━━━━┻━━━━━┓
                 ┃          27
                 ┃        ┏━┻━━┓
                 26       ┃    ┃
            ┏━━━━┻━━━━┓   ┃    ┃
            0         7   4    8
            ┃         ┃  ┏┻━┓  ┃
            21        ┃  ┃  ┃  ┃
    ┏━━┳━━┳━┻┳━━┳━━┓  ┃  ┃  ┃  ┃
    14 19 10 13 16 18 12 15 17 11

The trees have fully coalesced and we've successfully combined a forwards-time
Wright-Fisher simulation with a coalescent simulation: hooray!

-------------------------------------
Why record the initial generation?
-------------------------------------

We can now see why it is essential that the forwards simulator records the
*initial* generation in a tree sequence that will later be used as a
``from_ts`` argument to ``msprime.simulate()``. In the example above, if node
``7`` was not in the tree sequence, we would not know that the segment that
node ``20`` inherits from on ``[0.0, 1.0)`` and the segment that node ``12``
inherits from on ``[1.0, 2.0)`` both exist in the same node (here, node ``7``).

However, note that although the intial generation (above, nodes ``0``, ``4``,
``7``, and ``8``) must be in the tree sequence, they do *not* have to be
samples. The easiest way to do this is to
(a) retain the initial generation as samples throughout the forwards simulation
(so they persist through ``simplify()``), but then (b) before we output
the final tree sequence, we remove the flags that mark them as samples,
so that ``simulate()`` does not simulate their entire history as well. This
is the approach taken in the toy simulator provided above (although we skip
the periodic ``simplify()`` steps which are essential in any practical simulation
for simplicity).

-------------------------------------
Topology gotchas
-------------------------------------

The trees that we output from this combined forwards and backwards simulation
process have some slightly odd properties that are important to be aware of.
In the example above, we can see that the old roots are still present in both trees,
even through they have only one child and are clearly redundant.
This is because the tables of ``from_ts`` have been retained, without modification,
at the top of the tables of the output tree sequence. While this
redundancy is not important for many tasks, there are some cases where
they may cause problems:

1. When computing statistics on the number of nodes, edges or trees in a tree
   sequence, having these unary edges and redundant nodes will slightly
   inflate the values.
2. If you are computing the overall tree "height" by taking the time of the
   root node, you may overestimate the height because there is a unary edge
   above the "real" root (this would happen if one of the trees had already
   coalesced in the forwards-time simulation).

For these reasons it is usually better to remove this redundancy from your
computed tree sequence which is easily done using the
:meth:`.TreeSequence.simplify` method:

.. code-block:: python

    final_ts = coalesced_ts.simplify()

    for tree in final_ts.trees():
        print("interval = ", tree.interval)
        print(tree.draw(format="unicode"))

giving us::

    interval =  (0.0, 1.0)
          17
      ┏━━━┻━━━━┓
      ┃        15
      ┃     ┏━━┻━━┓
      ┃     13    14
      ┃   ┏━┻┓  ┏━╋━━┓
      10  ┃  11 ┃ ┃  12
    ┏━╋━┓ ┃ ┏┻┓ ┃ ┃ ┏┻┓
    2 3 6 0 4 9 1 8 5 7

    interval =  (1.0, 2.0)
              19
        ┏━━━━━┻━━━━━┓
        ┃           18
        ┃         ┏━┻┓
        17        ┃  ┃
    ┏━━━┻━━┓      ┃  ┃
    ┃      ┃      ┃  16
    ┃      ┃      ┃ ┏┻┓
    ┃      11     ┃ ┃ ┃
    ┃ ┏━┳━┳┻┳━┳━┓ ┃ ┃ ┃
    2 4 9 0 3 6 8 1 5 7

This final tree sequence is topologically identical to the original tree sequence,
but has the redundant nodes and edges removed. Note also that he node IDs have been
reassigned so that the samples are 0 to 9 --- if you need the IDs from the original
tree sequence, please set ``map_nodes=True`` when calling ``simplify`` to get a
mapping between the two sets of IDs.
