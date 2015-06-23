.. _sec-cli:

======================
Command line interface
======================

This is the documentation for the :command:`mspms` program, an :command:`ms`-compatible
command line interface to the ``msprime`` library. This interface should
be useful for legacy applications, where it can be used as a drop-in
replacement for :command:`ms`. This interface is not recommended for new applications,
particularly if the simulated trees are required as part of the output
as Newick is very inefficient. The :ref:`Python API <sec-api>` is the recommended interface,
providing direct access to the structures used within ``msprime``.

.. note::

   We use :command:`mspms` for the :command:`ms`-compatible command line interface as another
   (fully POSIX compliant) CLI called :command:`msp` is planned. However,
   :command:`mspms` will be supported for the forseeable future.

******************
Supported Features
******************

:command:`mspms` supports a subset of :command:`ms`'s functionality. Please
`open an issue <https://github.com/jeromekelleher/msprime/issues>`_ on
GitHub if there is a feature of :command:`ms` that you would like to see
added. We  currently support:

- Basic functionality (sample size, replicates, tree and haplotype output);
- Recombination (via the ``-r`` option);
- Exponentially growing/shrinking population size (via the ``-G`` option);
- Demographic events (via the ``-eG`` option and ``-eN`` options).

Spatial structure and gene-conversion are not currently supported, but
are planned for future releases.

****************
Argument details
****************

This section provides the detailed listing of the arguments to
:command:`mspms` (also available via ``mspms --help``). See
the `documentation for ms
<http://thirteen-01.stat.iastate.edu/snoweye/phyclust/document/msdoc.pdf>`_
for details on how these values should be interpreted.

.. argparse::
    :module: msprime.cli
    :func: get_parser
    :prog: mspms
    :nodefault:
