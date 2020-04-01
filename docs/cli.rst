.. _sec_cli:

======================
Command line interface
======================

Two command_line applications are provided with ``msprime``: :ref:`sec_msp` and
:ref:`sec_mspms`. The :command:`msp` program is an experimental interface for
interacting with the library, and is a POSIX compliant command line
interface. The :command:`mspms` program is a fully-:command:`ms` compatible
interface. This is useful for those who wish to get started quickly with using
the library, and also as a means of plugging ``msprime`` into existing work
flows. However, there is a substantial overhead involved in translating data
from ``msprime``'s native history file into legacy formats, and so new code
should use the :ref:`Python API <sec_api>` where possible.

.. _sec_msp:

***
msp
***

The ``msp`` program provides a convenient interface to the :ref:`msprime API
<sec_api>`. It is based on subcommands that either generate or consume a
tree sequence file. The ``simulate`` subcommand runs a
simulation storing the results in a file. The other commands are concerned with
converting this file into other formats.

.. warning:: This tool is very new, and the interface may need to change
    over time. This should be considered an alpha feature!

++++++++++++
msp simulate
++++++++++++

:command:`msp simulate` provides a command line interface to the
:func:`msprime.simulate` API function. Using the parameters provided at the
command line, we run a simulation and then save the resulting tree sequence
to the file provided as an argument.

.. argparse::
    :module: msprime.cli
    :func: get_msp_parser
    :prog: msp
    :path: simulate
    :nodefault:

.. note:: The way in which recombination and mutation rates are specified
    is different to :command:`ms`. In :command:`ms` these rates are scaled by the
    length of the simulated region, whereas we use rates per unit distance.
    The rationale for this change is to simplify running simulations on a
    variety of sequence lengths, so that we need to change only one parameter
    and not three simultaneously. See :ref:`sec_api` for more on this point.


.. TODO remove this information and add deprecation notices for the various
.. commands once the tskit CLI has been implemented.

.. .. _sec_msp_upgrade:

.. +++++++++++
.. msp upgrade
.. +++++++++++

.. :command:`msp upgrade` is a command line tool to convert tree sequence
.. files written by older versions of msprime to the latest version.
.. This tool requires `h5py <http://www.h5py.org/>`_, so please ensure that
.. it is installed. The upgrade process involves creating a new tree sequence
.. file from the records stored in the older file and is non-destructive.

.. .. argparse::
..     :module: msprime.cli
..     :func: get_msp_parser
..     :prog: msp
..     :path: upgrade
..     :nodefault:

.. ++++++++++++
.. msp vcf
.. ++++++++++++

.. :command:`msp vcf` is a command line interface to the
.. :meth:`msprime.TreeSequence.write_vcf` method. It prints out the coalescence
.. vcf in a history file in a tab-delimited text format.

.. .. argparse::
..     :module: msprime.cli
..     :func: get_msp_parser
..     :prog: msp
..     :path: vcf
..     :nodefault:

.. ++++++++++
.. msp newick
.. ++++++++++

.. :command:`msp newick` prints out the marginal genealogies in the tree
.. sequence in newick format.

.. .. argparse::
..     :module: msprime.cli
..     :func: get_msp_parser
..     :prog: msp
..     :path: newick
..     :nodefault:

.. ++++++++++++++++++++++++++++++++++++++++++++++++++++
.. msp (nodes, edges, sites, mutations, or provenances)
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++

.. The commands
.. :command:`msp nodes`,
.. :command:`msp edges`,
.. :command:`msp sites`,
.. :command:`msp mutations`, and
.. :command:`msp provenances`
.. each print out the respective table in tabular format from the tree sequence.
.. See :ref:`sec_interchange` for a description of these tables.

.. .. argparse::
..     :module: msprime.cli
..     :func: get_msp_parser
..     :prog: msp
..     :path: nodes
..     :nodefault:


.. ++++++++++++++
.. msp haplotypes
.. ++++++++++++++

.. :command:`msp haplotypes` prints out the haplotypes of each sampled genome
.. described in the tree sequence.
.. This only works with single-character allelic states.

.. .. argparse::
..     :module: msprime.cli
..     :func: get_msp_parser
..     :prog: msp
..     :path: haplotypes
..     :nodefault:


.. .. todo::
..     Provide individuals and populations commands.

.. _sec_mspms:

*****
mspms
*****

The :command:`mspms` program is an :command:`ms`-compatible
command line interface to the ``msprime`` library. This interface should
be useful for legacy applications, where it can be used as a drop-in
replacement for :command:`ms`. This interface is not recommended for new applications,
particularly if the simulated trees are required as part of the output
as Newick is very inefficient. The :ref:`Python API <sec_api>` is the recommended interface,
providing direct access to the structures used within ``msprime``.


++++++++++++++++++
Supported Features
++++++++++++++++++

:command:`mspms` supports a subset of :command:`ms`'s functionality. Please
`open an issue <https://github.com/tskit-dev/msprime/issues>`_ on
GitHub if there is a feature of :command:`ms` that you would like to see
added. We  currently support:

- Basic functionality (sample size, replicates, tree and haplotype output);
- Recombination (via the ``-r`` option);
- Spatial structure with arbitrary migration matrices;
- Support for :command:`ms` demographic events. (The implementation of the
  ``-es`` option is limited, and has restrictions on how it may be
  combined with other options.)

Gene-conversion is not currently supported, but is planned for a future release.

++++++++++++++++
Argument details
++++++++++++++++

This section provides the detailed listing of the arguments to
:command:`mspms` (also available via ``mspms --help``). See
the `documentation for ms
<http://thirteen-01.stat.iastate.edu/snoweye/phyclust/document/msdoc.pdf>`_
for details on how these values should be interpreted.

.. argparse::
    :module: msprime.cli
    :func: get_mspms_parser
    :prog: mspms
    :nodefault:


