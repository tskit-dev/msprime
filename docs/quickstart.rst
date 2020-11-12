.. _sec_quickstart:

==========
Quickstart
==========

.. todo:: This is a rough draft and needs more work.

*************
Installation
*************

If you have `conda <https://docs.conda.io/en/latest/>`_ installed::

    $ conda install -c conda-forge msprime

otherwise::

    $ python3 -m pip install msprime

If these commands don't work, or for more installation options,
please see the :ref:`sec_installation` page for more details.


********
Ancestry
********

Msprime can simulate ancestral histories for a set of sample
genomes under a variety of evolutionary models. The default model
is the `coalescent <https://en.wikipedia.org/wiki/Coalescent_theory>`_,
which assumes a single randomly mating population of a fixed size.
In the simplest case of no
`recombination <https://en.wikipedia.org/wiki/Genetic_recombination>`_
the result of an ancestry simulation is a genealogical `tree
<https://en.wikipedia.org/wiki/Phylogenetic_tree>`_ relating the simulated
samples to each other and their genetic ancestors. Msprime
can efficiently simulate :ref:`recombination <sec_ancestry_recombination>`
and other processes which result in *multiple* trees along the
genome. The output of an ancestry simulation is a therefore
tree *sequence* which we use the :ref:`tskit <tskit:sec_introduction>`
library to represent. Tskit has a rich set of
features for analysing these genealogical histories.

.. jupyter-kernel:: python3


.. jupyter-execute::

    import msprime
    from IPython.display import SVG

    # Simulate an ancestral history for 3 diploid samples under the coalescent
    ts = msprime.sim_ancestry(3)
    # Visualise the simulated ancestral history.
    SVG(ts.draw_svg())


.. todo::
    We want a list of quick pointers here to relevant sections of the
    documentation and tutorials.


*********
Mutations
*********

The :func:`.sim_ancestry` function generates a simulated ancestral
history for some samples. This is often all we need for many purposes.
If we want `genome sequence <https://en.wikipedia.org/wiki/Genome>`_
we must also simulate some
`mutations <https://en.wikipedia.org/wiki/Mutation>`_ on these trees.

.. fixme This should use sim_mutations

.. jupyter-execute::

    # Simulate an ancestral history for 3 diploid samples under the coalescent
    ts = msprime.sim_ancestry(3)
    mutated_ts = msprime.mutate(ts, rate=0.1)
    SVG(mutated_ts.draw_svg())

.. todo:: Some example chunks where we show how to do something simple
    with the sequences and maybe how to export to VCF.


.. todo:: List of pointers to the relevant sections of the documentation.


**********
Demography
**********

By default ancestry simulations assume an extremely simple
population structure in which a single randomly mating population
of a fixed size exists for all time. For most simulations this
is an unrealistic assumption, and so msprime provides a way
to describe more complex demographic histories.

.. todo:: The demography.sample example here isn't a very good one
    and the .sample function needs to updated.

.. jupyter-execute::

    # Create a 1D stepping stone model of demograpy
    demography = msprime.Demography.stepping_stone_1d(10, migration_rate=0.1)
    # Take one diploid sample each from the first two demes
    samples = demography.sample(1, 1)
    # Simulate an ancestral history for this demography and sample.
    ts = msprime.sim_ancestry(samples=samples, demography=demography)
    ts.tables.nodes

.. todo:: This example is rough and needs to be updated once the demography
    API is sharpened up a bit.
