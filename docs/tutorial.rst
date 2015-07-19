.. _sec-tutorial:

========
Tutorial
========

This is the tutorial for the Python interface to the ``msprime``
library. Detailed :ref:`sec-api` is also availale for this
library. An ``ms``-compatible :ref:`command line interface <sec-cli>`
is also available if you wish to use ``msprime`` directly within
an existing workflow.


****************
Simulating trees
****************

Running simulations is very straightforward in ``msprime``::

    >>> import msprime
    >>> tree = msprime.simulate(5)
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

Simulating the history of a single locus is a very useful, but we are
most often interesting in simulating the history of our sample across
large genomic regions under the influence of recombination. The ``msprime``
API is specifically designed to make this common requirement both easy
and efficient. Throughout ``msprime`` we use the term locus to
refer to a non-recombining strech of sequence. Recombination occurs
between loci but never within loci. It is often most straightforward
to regard a single base pair as a 'locus'. In ``msprime`` we support
up to :math:`2^{32}` loci, which is sufficient to represent all but
the very largest of known chromosome sizes.

**TODO** explain recombination rate. This is the same as used in
``ms`` but we do not multiply by the number of loci.

We simulate the trees across over a number of loci using the
:func:`msprime.simulate()` function::

    >>> tree_sequence = msprime.simulate(
    >>>     5, num_loci=10, scaled_recombination_rate=0.1)
    >>> for tree in tree_sequence.sparse_trees():
    >>>     print(tree.get_interval(), str(tree), sep="\t")
    (0, 6)  {1: 9, 2: 8, 3: 8, 4: 7, 5: 7, 7: 9, 8: 10, 9: 10, 10: 0}
    (6, 10) {1: 9, 2: 6, 3: 10, 4: 6, 5: 7, 6: 7, 7: 9, 9: 10, 10: 0}

In this example, we simulate the history of our sample of 5 individuals
over 10 loci, with a scaled recombination rate of 0.1 between adjacent
pairs of loci. Unlike the :func:`msprime.simulate_tree` function which
returns a tree, the :func:`msprime.simulate` function returns a
*tree sequence*, which encapsulates all of the information in the
sequence of correlated trees over the simulated region. The
:class:`msprime.TreeSequence` class provides an array of methods to
simplify working with these trees and some efficient methods for
common tasks that take advantage of the strong correlation structure
of the trees in the sequence.

In this example, we use the :meth:`~msprime.TreeSequence.sparse_trees`
method to iterate over the trees in the sequence. For each tree
we print out the interval the tree covers (i.e., the genomic
coordinates which all share precisely this tree) using the
:meth:`~method.TreeSequence.get_interval` method. We also print
out the summary of each tree in terms of the parent values for
each tree. Again, these differences are best illustrated by
some images:

.. image:: _static/simple-tree-sequence-0.svg
   :width: 200px
   :alt: A simple coalescent tree

.. image:: _static/simple-tree-sequence-1.svg
   :width: 200px
   :alt: A simple coalescent tree

(We have supressed the node time labels here for clarity.) We can see
that these trees share a great deal of their structure, but that there are
also important differences between the trees.


.. warning:: Do not store the values returned from the
    :meth:`~msprime.TreeSequence.sparse_trees` iterator in a list and operate
    on them afterwards! For efficiency reasons ``msprime`` uses the same
    instance of :class:`msprime.SparseTree` for each tree in the sequence
    and updates the internal state for each new tree. Therefore, if you store
    the trees returned from the iterator in a list, they will all refer
    to the same tree.


*********
Mutations
*********

Mutations are generated in msprime
