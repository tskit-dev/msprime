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
    (0, 1):{1: 6, 2: 7, 3: 9, 4: 6, 5: 7, 6: 8, 7: 8, 8: 9, 9: 0}

Here, we simulate the coalescent process for a sample of size
5 and print out a summary of the resulting tree. Trees are
represented within ``msprime`` in a slightly unusual way. In
the majority of libraries dealing with trees, each node is
explicitly represented as an object in memory and the relationship
between nodes as pointers between these objects. In ``msprime``,
however, nodes are *integers*: the leaves (i.e., our sample) are the
integers :math:`1` to :math:`n`, and every internal node is
some positive integer greater than :math:`n`.
