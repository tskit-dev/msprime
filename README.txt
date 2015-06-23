=======
msprime
=======

Msprime is a reimplementation of Hudson's classical ms program for modern
datasets. The Python API and storage format are currently under development
and are not fully documented, but the command line interface ``mspms`` is
reliable and ready for use. This program provides a fully ``ms`` compatible
interface, and can be used as a drop-in replacement in existing workflows.

***********
Quick Start
***********

To install and run ``msprime`` on a Ubuntu 14.10 installation do the
following::

    $ sudo apt-get install python-dev libgsl0-dev libhdf5-serial-dev
    $ sudo pip install msprime
    $ mspms 2 1 -t 1
    /usr/local/bin/mspms 2 1 -t 1
    5338 8035 23205

    //
    segsites: 3
    positions: 0.014 0.045 0.573
    100
    011


If you do not wish to install ``msprime`` to your system, you can try
it out in a `virtualenv <http://virtualenv.pypa.io/en/latest/>`_ as
follows::

    $ virtualenv msprime-env
    $ source msprime-env/bin/activate
    (msprime-env) $ pip install msprime
    (msprime-env) $ mspms


*************
Requirements
*************

Msprime requires Python 2.7+ (Python 3 versions are fully supported from
3.1 onwards), the `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_,
and `HDF5 <https://www.hdfgroup.org/HDF5/>`_ version 1.8 or later. These
packages are available for all major platforms. For example, to install on
Debian/Ubuntu use::

    # apt-get install python-dev libgsl0-dev libhdf5-serial-dev pkg-config

For Redhat/Fedora use::

    # yum install gsl-devel hdf5-devel

On FreeBSD we can use ``pkg`` to install the requirements::

    # pkg install gsl hdf5-18

----
TODO
----

- Add instructions for HomeBrew/MacPorts.

*************
Installation
*************

The simplest method of installation is to use PyPI and pip::
    # pip install msprime

will install msprime your Python installation. This will
work in most cases (but see the `Requirements`_ section). See
below for platform specific build instructions when this
fails.

If you do not have root access to your machine, you can install
``msprime`` into your local Python installation as follows::

    $ pip install msprime --user

To use the ``mspms`` you must ensure
that the ``~/.local/bin/`` directory is in your ``PATH``, or
simply run it using::

    $ ~/.local/bin/mspms

To uninstall ``msprime``, simply run::

    $ pip uninstall msprime

------------------------------
Platform specific installation
------------------------------

This section contains instructions to build on specific platforms
that require build time flags for GSL and HDF5.

++++++++++++
FreeBSD 10.0
++++++++++++

Install the prerequisitites, and build ``msprime`` as follows::

    # pkg install gsl hdf5-18
    # CFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib pip install msprime

This assumes that root is logged in using a bash shell. For other shells,
different methods are need to get the ``CFLAGS`` and ``LDFLAGS`` environment
variables.

---
OSX
---

*TODO*
