=======
msprime
=======

Msprime is a reimplementation of Hudson's classical ms program for modern
datasets. The Python API and storage format are currently under development
and are not fully documented, but the command line interface ``mspms`` is
reliable and ready for use. This program provides a fully ``ms`` compatible
interface, and can be used as a drop-in replacement in existing workflows.

*************
Requirements
*************

Msprime requires Python 2.7+ (Python 3 versions are fully supported from
3.1 onwards), the `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_,
and `HDF5 <https://www.hdfgroup.org/HDF5/>`_ version 1.8 or later. These
packages are available for all major platforms. For example, to install on
Debian/Ubuntu use::

    # apt-get install python-dev libgsl0-dev libhdf5-serial-dev pkg-config



----
TODO
----

- Need installation instructions for an RPM based distribution,
  HomeBrew/MacPorts and possibly a linux distro where we are missing
  root access.

------------
FreeBSD 10.0
------------

To install on FreeBSD, install the gsl and hdf5-18 packages as
follows::

    # pkg install gsl hdf5-18

**TODO** Finish installation instructions.

*************
Installation
*************

The simplest method of installation is to use PyPI and pip::

    $ pip install msprime --user --pre

will install msprime your user Python installation. This should
work in most cases (but see the `Requirements`_ section).

To use the ms-compatible command line program, you must ensure
that the ``~/.local/bin/`` directory is in your ``PATH``.


