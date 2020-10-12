.. _sec_ancestry:

===================
Simulating ancestry
===================


*****************
API Documentation
*****************

.. autofunction:: msprime.sim_ancestry()


********
Examples
********

.. _sec_ancestry_samples_ploidy:

++++++++++++++++++
Samples and ploidy
++++++++++++++++++

.. jupyter-kernel:: python3


.. jupyter-execute::
    :hide-code:

    import msprime
    from IPython.display import SVG

The samples argument for :func:`sim_ancestry` is flexible, and allows us
to provide samples in a number of different forms. In single-population
models we can use the numeric form, which gives us :math:`n` samples:

.. jupyter-execute::

    ts = msprime.sim_ancestry(3, random_seed=42)
    SVG(ts.first().draw_svg())


It's important to note that the number of samples refers to the number
of :math:`k`-ploid *individuals*, not the number of sample nodes
in the trees. The ``ploidy`` argument determines the number of sample
nodes per individual, and is ``2`` by default; hence, when we asked
for 3 sample individuals in the example above, we got a tree with
six sample *nodes*.

.. jupyter-execute::

    ts = msprime.sim_ancestry(3, ploidy=1, random_seed=42)
    SVG(ts.first().draw_svg())


.. todo examples of drawing samples from a demography.


.. _sec_ancestry_genome_length:

+++++++++++++
Genome length
+++++++++++++

.. jupyter-kernel:: python3

.. jupyter-execute::
    :hide-code:

    import msprime
    from IPython.display import SVG

There are a number of different ways to specify the length of the
chromosome that we want to simulate. In the absence of recombination
and gene conversion, we assume a genome of length 1:


.. jupyter-execute::

    ts = msprime.sim_ancestry(3)
    ts.sequence_length


If a recombination or gene conversion rate is specified, though, we
must define a ``sequence_length``:

.. jupyter-execute::

    ts = msprime.sim_ancestry(
        3, recombination_rate=0.1, sequence_length=10, random_seed=2)
    ts.sequence_length, ts.num_trees

.. jupyter-execute::
    :hide-code:

    assert ts.num_trees > 1

In this example we have a uniform recombination rate between all
positions along the genome. We can also simulate variable
recombination rates along the genome using a :class:`RateMap`.

.. jupyter-execute::

    rate_map = msprime.RateMap(
        position=[0, 10, 12, 20],
        rate=[0.1, 0.5, 0.1])
    ts = msprime.sim_ancestry(3, recombination_rate=rate_map, random_seed=2)
    ts.sequence_length, ts.num_trees

Here we specify varying recombination rates for a genome of length 20,
and there's a hotspot from position 10 to 12. In this case we don't
need to specify the ``sequence_length`` in the call to ``sim_ancestry``
because it's already defined by the :class:`RateMap`.



.. _sec_ancestry_discrete_genome:

+++++++++++++++
Discrete genome
+++++++++++++++

.. jupyter-kernel:: python3


.. jupyter-execute::
    :hide-code:

    import msprime
    from IPython.display import SVG


By default, we assume that the genome we are simulating is *discrete*
so that genome coordinates are at integer positions:

.. jupyter-execute::

    ts = msprime.sim_ancestry(
        3, recombination_rate=0.1, sequence_length=10, random_seed=2)
    ts.sequence_length

.. jupyter-execute::
    :hide-code:

    assert 1 < ts.num_trees < 5


.. jupyter-execute::

    SVG(ts.draw_svg())

We can also simulate a continous genome by setting
``discrete_genome=False``:

.. jupyter-execute::

    ts = msprime.sim_ancestry(
        3, recombination_rate=0.25, sequence_length=1, discrete_genome=False,
        random_seed=33)
    SVG(ts.draw_svg())

.. jupyter-execute::
    :hide-code:

    assert 1 < ts.num_trees < 5

Here we see that the breakpoints along the genome occur at floating point
positions. Simulating a continuous genome sequence can be useful for
theoretical work, but we recommend using discrete coordinates for most
purposes.



.. _sec_ancestry_recombination:

+++++++++++++
Recombination
+++++++++++++

.. jupyter-kernel:: python3


.. jupyter-execute::
    :hide-code:

    import msprime
    from IPython.display import SVG

RECOMB

.. jupyter-execute::

    ts = msprime.sim_ancestry(
        3, recombination_rate=0.1, sequence_length=10, random_seed=2)
    ts.sequence_length

.. jupyter-execute::
    :hide-code:

    assert 1 < ts.num_trees < 5

.. _sec_ancestry_gene_conversion:

+++++++++++++++
Gene conversion
+++++++++++++++

.. jupyter-kernel:: python3


.. jupyter-execute::
    :hide-code:

    import msprime
    from IPython.display import SVG

GC

.. jupyter-execute::

    ts = msprime.sim_ancestry(
        3, recombination_rate=0.1, sequence_length=10, random_seed=2)
    ts.sequence_length

.. jupyter-execute::
    :hide-code:

    assert 1 < ts.num_trees < 5

