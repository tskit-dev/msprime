.. _sec-introduction:

============
Introduction
============

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
   :ref:`Python API <sec-api>` to simplify the workflow associated with
   running and analysing simulations. (However, we do provide an
   ``ms``-compatible :ref:`command line interface <sec-cli>` to
   plug in to existing workflows.) For many simulations we first
   write a script to generate the command line parameters we
   want to run, then fork shell processes to run the simulations,
   and then parse the results to obtain the genealogies in a form
   we can use. With ``msprime`` all of this can be done directly
   in Python, which is both simpler and far more efficient.

3. ``msprime`` does not use Newick trees for interchange as they
   are extremely inefficient in terms of the time required to
   generate and parse, as well as the space required to store them.
   Instead, we use a :ref:`well-defined <sec-file-format>` format using the
   powerful `HDF5 <https://www.hdfgroup.org/HDF5/>`_ standard. This
   format allows us to store genealogical data very concisely,
   particularly for large sample sizes.

