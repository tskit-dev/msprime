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

.. jupyter-execute::

    # Create a 1D stepping stone model of demograpy
    demography = msprime.Demography.stepping_stone_model([100] * 10, migration_rate=0.1)
    # Take one diploid sample each from the first and last demes
    samples = {0: 1, 9: 1}
    # Simulate an ancestral history for this demography and sample.
    ts = msprime.sim_ancestry(samples=samples, demography=demography)
    ts.tables.nodes

.. todo:: Links into more detailed documentation

******************
Upgrading from 0.x
******************

This section is to help 0.x users of the get up to speed quickly, summarising the new
APIs and their main differences to what you are used to.

The main change is that there are two new functions, :func:`.sim_ancestry` and
:func:`.sim_mutations` which correspond to the 0.x functions :func:`.simulate`
and :func:`.mutate`. The 0.x functions are **deprecated** but **will continue
to be supported indefinitely**.

+++++++++++++++++++++++
Backwards compatibility
+++++++++++++++++++++++

All existing simulations should work as before, *except* for simulations relying on
the detailed properties of RecombinationMaps. If your code uses the ``num_loci``
property, then it may need to be updated. The reason for this is that ``msprime``
has changed to simulate directly in physical coordinates internally (which has
greatly simplified the code and solved many thorny issues) and this is fundamentally
incompatible with the approach taken in 0.x. In the vast majority of cases, this
will have no effect.

If you are using ``num_loci`` to simulate a discrete genome, it may be simplest to
convert your code to use the new :func:`.sim_ancestry` method. If you were following
a recipe to simulate multiple chromosomes under the DTWF model, please see
the :ref:`updated recipe <sec_ancestry_multiple_chromosomes>`.

++++++++
Ancestry
++++++++

The new :func:`.sim_ancestry` function replaces the 0.x :func:`.simulate`
function and is very similar. There are some important differences though:

* Samples are now **individuals** rather **nodes** (i.e. monoploid
  genomes). Because the default :ref:`ploidy <sec_ancestry_samples_ploidy>`
  is 2 (see the next point) the upshot is that ``sim_ancestry(2)`` will
  result in a tree sequence with *four* sample nodes, not two. (It is
  possible to override this behaviour using the list of :class:`.SampleSet`
  objects argument to ``samples``.)

* There is now a :ref:`sec_ancestry_samples_ploidy` argument, which has
  two effects:

  #. Sets the default number of sample nodes per *individual*

  #. Changes the timescale of the coalescent process (TODO link to a section
     that explains this effect.) By default ``ploidy`` is 2 and
     time is scaled scaled in units of 4N generations, which is the same as
     msprime 0.x.

* Rather than two arguments ``num_samples`` and ``samples``, the
  :func:`.sim_ancestry` function has a single argument ``samples`` which
  has different behaviour depending on the type of arguments provided.
  See :ref:`sec_ancestry_samples` for details.

  Note in particular that a list of ``Sample`` objects is **not** supported.

* Similarly, there is now one argument ``recombination_rate`` which can
  be either a single value or a :class:`.RateMap` object. Note that the
  0.x :class:`.RecombinationMap` is deprecated and not supported as input
  to :func:`.sim_ancestry`. See :ref:`sec_ancestry_recombination` for more
  details.

* Simulations are peformed on a **discrete** genome by default. To get the
  0.x behaviour of a continuous genome, set ``discrete_genome=False``.
  See :ref:`sec_ancestry_discrete_genome` for more details.

* The ``from_ts`` argument used has been renamed to ``initial_state`` and
  accepts either a :class:`tskit.TableCollection` or :class:`tskit.TreeSequence`
  argument. See :ref:`sec_ancestry_initial_state` for details.

* There is **no** ``mutation_rate`` argument to :func:`.sim_ancestry`: use
  :func:`.sim_mutations` instead.

* The ``population_configurations``, ``migration_matrix`` and ``demographic_events``
  arguments have been replace with a single argument ``demography``, which must take
  a :class:`.Demography` instance. (See the next section for more details.)

++++++++++
Demography
++++++++++

* A new :class:`.Demography` object has been added for version 1.0 which
  encapsulates the functionality needed to define and debug demographic models
  in msprime. Demographic models can only be specified to ``sim_ancestry``
  using an instance of this class.

* It is easy to create a :class:`.Demography` from the 0.x
  ``population_configurations``, ``migration_matrix`` and ``demographic_events``
  values using the :meth:`.Demography.from_old_style` method.

* The :class:`.DemographyDebugger` class should no longer be instantiated
  directly; instead use the :meth:`.Demography.debug` method.

+++++++++
Mutations
+++++++++

* For symmetry with the :func:`.sim_ancestry` function, there is now a :func:`.sim_mutations`
  function. The 0.x :func:`.mutate` function is **deprecated**.

* The :func:`.sim_mutations` function works on a **discrete** genome by default.

* There are now also many new mutation models which support;
  see :ref:`sec_mutations` for details. These are *not* supported in the deprecated
  :func:`.mutate` function.


+++++++++
Utilities
+++++++++

* The 0.x class :class:`.RecombinationMap` has been **deprecated** in favour of the new
  :class:`.RateMap`. This was to (a) generalise the interface to accomodate varying
  rates of mutation and gene conversion along the genome; and (b) convert to a
  more modern numpy-based API.
