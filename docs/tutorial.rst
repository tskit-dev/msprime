.. _sec-tutorial:

========
Tutorial
========

This is the tutorial for the Python interface to the ``msprime``
library. Detailed :ref:`sec-api` is also available for this
library. An :program:`ms`-compatible :ref:`command line interface <sec-cli>`
is also available if you wish to use ``msprime`` directly within
an existing work flow.


****************
Simulating trees
****************

Running simulations is very straightforward in ``msprime``::

    >>> import msprime
    >>> tree = msprime.simulate_tree(5)
    >>> print(tree)
    {1: 6, 2: 8, 3: 6, 4: 8, 5: 7, 6: 7, 7: 9, 8: 9, 9: 0}

Here, we simulate the coalescent process for a sample of size
5 and print out a summary of the resulting tree. Trees are
represented within ``msprime`` in a slightly unusual way. In
the majority of libraries dealing with trees, each node is
represented as an object in memory and the relationship
between nodes as pointers between these objects. In ``msprime``,
however, nodes are *integers*: the leaves (i.e., our sample) are the
integers :math:`1` to :math:`n`, and every internal node is
some positive integer greater than :math:`n`. The result of printing
the tree is a summary of how these nodes relate to each other
in terms of their parents. For example, we can see that the parent
of nodes 1 and 3 is node 6.

This relationship can be seen more clearly in a picture:

.. image:: _static/simple-tree.svg
   :width: 200px
   :alt: A simple coalescent tree

This image shows the same tree as in the example but drawn out in
a more familiar format (images like this can be drawn for any
tree using the :meth:`~msprime.SparseTree.draw` method).
We can see that the leaves of the tree
are labelled with 1 to 5, and all the internal nodes of the tree
are also integers with the root of the tree being 9. Also shown here
are the times for each internal node, in coalescent time units. (The
time for all leaves is 0, and so we don't show this information
to avoid clutter.)

Knowing that our leaves are 1 to 5, we can easily trace our path
back to the root for a particular sample using the
:meth:`~msprime.SparseTree.get_parent` method::

    >>> u = 1
    >>> while u != 0:
    >>>     print("node {}: time = {}".format(u, tree.get_time(u)))
    >>>     u = tree.get_parent(u)
    node 1: time = 0.0
    node 6: time = 0.0269802913256
    node 7: time = 0.251686777821
    node 9: time = 0.446340881302

In this code chunk we iterate up the tree starting from node 1 and
stopping when we get to the root. We know that a node is the root
if its parent is 0, which is a special, reserved node. We also use
the :meth:`~msprime.SparseTree.get_time` method to get the time
for each node, which corresponds to the time at which the coalescence
event happened during the simulation (in coalescent time units).
We can also obtain the length of a branch joining a node to
its parent using the :meth:`~msprime.SparseTree.get_branch_length`
method::

    >>> print(tree.get_branch_length(7))
    0.194654103481

The branch length for node 7 is 0.19 as the time for node 7 is 0.25,
and the time of its parent is 0.44.

*************
Multiple loci
*************

Simulating the history of a single locus is a very useful, but we are most
often interesting in simulating the history of our sample across large genomic
regions under the influence of recombination. The ``msprime`` API is
specifically designed to make this common requirement both easy and efficient.
Throughout ``msprime`` we use the term locus to refer to a non-recombining
stretch of sequence. Recombination occurs between loci but never within loci.
It is often most straightforward to regard a single base pair as a 'locus'. In
``msprime`` we support up to :math:`2^{32}` loci, which is sufficient to
represent all but the very largest of known chromosome sizes. See the
:ref:`sec-api` for a discussion of the precise recombination model used.


We simulate the trees across over a number of loci using the
:func:`msprime.simulate()` function::

    >>> tree_sequence = msprime.simulate(
    ...     5, num_loci=10, scaled_recombination_rate=0.1, random_seed=19)
    >>> for tree in tree_sequence.trees():
    ...     print(tree.get_interval(), str(tree), sep="\t")
    (0, 5)  {1: 7, 2: 6, 3: 7, 4: 10, 5: 6, 6: 8, 7: 8, 8: 10, 10: 0}
    (5, 10) {1: 7, 2: 6, 3: 7, 4: 9, 5: 6, 6: 9, 7: 10, 9: 10, 10: 0}

In this example, we simulate the history of our sample of 5 individuals
over 10 loci, with a scaled recombination rate of 0.1 between adjacent
pairs of loci. (We also provide the ``random_seed`` parameter here
as we wish to use this exact example again later; if we don't provide
a random seed, one is generated automatically.)
Unlike the :func:`msprime.simulate_tree` function which
returns a tree, the :func:`msprime.simulate` function returns a
*tree sequence*, which encapsulates all of the information in the
sequence of correlated trees over the simulated region. The
:class:`msprime.TreeSequence` class provides an array of methods to
simplify working with these trees and some efficient methods for
common tasks that take advantage of the strong correlation structure
of the trees in the sequence.

In this example, we use the :meth:`~msprime.TreeSequence.trees`
method to iterate over the trees in the sequence. For each tree
we print out the interval the tree covers (i.e., the genomic
coordinates which all share precisely this tree) using the
:meth:`~msprime.SparseTree.get_interval` method. We also print
out the summary of each tree in terms of the parent values for
each tree. Again, these differences are best illustrated by
some images:

.. image:: _static/simple-tree-sequence-0.svg
   :width: 200px
   :alt: A simple coalescent tree

.. image:: _static/simple-tree-sequence-1.svg
   :width: 200px
   :alt: A simple coalescent tree

(We have suppressed the node time labels here for clarity.) We can see
that these trees share a great deal of their structure, but that there are
also important differences between the trees.


.. warning:: Do not store the values returned from the
    :meth:`~msprime.TreeSequence.trees` iterator in a list and operate
    on them afterwards! For efficiency reasons ``msprime`` uses the same
    instance of :class:`msprime.SparseTree` for each tree in the sequence
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
is specified using the ``scaled_mutation_rate`` parameter to the
:func:`msprime.simulate` method. For example, to add some mutations
to our example above, we can use::

    >>> tree_sequence = msprime.simulate(
    >>>     5, num_loci=10, scaled_recombination_rate=0.1,
    >>>     scaled_mutation_rate=0.2, random_seed=19)
    >>> print("Total mutations = ", tree_sequence.get_num_mutations())
    >>> for tree in tree_sequence.trees():
    >>>     print(tree.get_interval(), list(tree.mutations()), sep="\t")
    Total mutations =  2
    (0, 5)  [(0.20106735406443477, 8)]
    (5, 10) [(9.032968991668895, 7)]

In this example (which has the same genealogies as our example above because
we use the same random seed), we generate a total of two mutations, which
happen to fall as one on each tree. Mutations are represented as a
tuple ``(position, node)``, where ``position`` is the location of the mutation
in genomic coordinates and ``node`` is the node in the tree above which the
mutation occurs. Positions are given as a floating point value as we are
using the infinite sites model. Every mutation falls on exactly one tree
and we obtain the mutations for a particular tree using the
:meth:`~msprime.TreeSequence.mutations` method. Mutations are always returned
in increasing order of position. The mutations for this example are shown
on the trees here as red boxes:

.. image:: _static/mutations-tree-sequence-0.svg
   :width: 200px
   :alt: A simple coalescent tree with mutations

.. image:: _static/mutations-tree-sequence-1.svg
   :width: 200px
   :alt: A simple coalescent tree with mutations

We can calculate the allele frequency of mutations easily and
efficiently using the :meth:`~msprime.SparseTree.get_num_leaves`
which returns the number of leaves underneath a particular node.
For example,::

    >>> for tree in tree_sequence.trees():
    ...    for position, node in tree.mutations():
    ...        print("Mutation @ position {} has frequency {}".format(
    ...            position, tree.get_num_leaves(node) / tree.get_sample_size()))
    Mutation @ position 0.201067354064 has frequency 0.8
    Mutation @ position 9.03296899167 has frequency 0.4
