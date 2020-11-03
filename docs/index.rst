===================================
Welcome to msprime's documentation!
===================================

This is the documentation for ``msprime``, a reimplementation of Hudson's
classical :command:`ms` simulator.

If you are looking for help on a specific issue or would like to ask a
question to fellow ``msprime`` users, please send an email to the
`mailing list <https://groups.google.com/group/msprime-users>`_. By asking
questions on the mailing list, you can help to build a searchable knowledge
base.

=====================
OLD Introduction text
=====================

.. todo:: refactor this to be more relevant for 1.0. In particular,
    we don't need to talk about ms so much.

The primary goal of ``msprime`` is to efficiently and conveniently
generate coalescent trees for a sample under a range of evolutionary
scenarios. The library is a reimplementation of Hudson's seminal
``ms`` program, and aims to eventually reproduce all its functionality.
``msprime`` differs from ``ms`` in some important ways:

1. ``msprime`` is *much* more efficient than ``ms``, both in terms of
   memory usage and simulation time. In fact, ``msprime`` is also
   much more efficient than simulators based on approximations to the
   coalescent with recombination model, especially for simulations
   with very large sample sizes. ``msprime`` can easily simulate
   chromosome sized regions for hundreds of thousands of samples.

2. ``msprime`` is primarily designed to be used through its
   :ref:`Python API <sec_api>` to simplify the workflow associated with
   running and analysing simulations. (However, we do provide an
   ``ms`` compatible :ref:`command line interface <sec_cli>` to
   plug in to existing workflows.) For many simulations we first
   write a script to generate the command line parameters we
   want to run, then fork shell processes to run the simulations,
   and then parse the results to obtain the genealogies in a form
   we can use. With ``msprime`` all of this can be done directly
   in Python, which is both simpler and far more efficient.

3. ``msprime`` does not use Newick trees for interchange as they
   are extremely inefficient in terms of storage space and the
   time needed to generate and parse them. Instead, we use the
   `tskit library <https://tskit.readthedocs.io/en/stable/index.html>`_
   which allows us to store and process very large scale simulation
   results efficiently.


.. This TOC is a work in progress while we refactor the documentation for 1.0.
   The idea is we refactor all the stuff that's currently in the tutorial and
   API sections into the appropriate high-level sections.

Contents:
=========

.. toctree::
   :maxdepth: 2

   installation
   quickref
   demography
   ancestry
   mutations
   likelihoods
   utilities
   cli
   development
   CITATION
   changelog

   tutorial
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

