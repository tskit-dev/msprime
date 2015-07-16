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
explicitly represented as an object in memory and the relationship
between nodes as pointers between these objects. In ``msprime``,
however, nodes are *integers*: the leaves (i.e., our sample) are the
integers :math:`1` to :math:`n`, and every internal node is
some positive integer greater than :math:`n`. The result of printing
of printing  the tree is a summary of how these nodes relate to
each other.

This relationship can be seen more clearly in a picture:

.. image:: _static/simple-tree.svg
   :width: 200px
   :alt: A simple coalescent tree

This image shows the same tree as in the example but drawn out in
a more familiar format (images like this can be drawn for any
tree using the :meth:`~msprime.SparseTree.draw` method).
We can see that the leaves of the tree
are labelled with 1 to 5, and all the internal nodes of the tree
are also integers with the root of the tree being 9.

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

Traversing from the leaves upwards using parent pointers is useful
in many instances, but it is also convenient to traverse the tree
downwards in some applications.
