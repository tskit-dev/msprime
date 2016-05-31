=======
msprime
=======

Msprime is a reimplementation of Hudson's classical ms program for modern
datasets. The Python API and storage format are currently alpha quality
and not fully documented, but the command line interface ``mspms`` is
reliable and ready for use. This program provides a fully ``ms`` compatible
interface, and can be used as a drop-in replacement in existing workflows.

Msprime can simulate the coalescent with recombination much
faster than programs based on the Sequentially Markov Coalescent
for large sample sizes and has very reasonable memory requirements. Msprime
makes it possible to simulate chromosome sized regions with hundreds of
thousands of samples.

If you use ``msprime`` in your work, please cite the `PLOS Computational
Biology <http://dx.doi.org/10.1371/journal.pcbi.1004842>`_ paper.
See `here <https://msprime.readthedocs.org/en/stable/CITATION.html>`_ for
full citation details.

Please see the `documentation <https://msprime.readthedocs.org/en/stable/>`_
for further details.

***********
Quick Start
***********

To install and run ``msprime`` on a fresh Ubuntu 15.10 installation, do the
following::

    $ sudo apt-get install pkg-config python-dev python-pip libgsl0-dev hdf5-tools libhdf5-serial-dev
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

See below for installation instructions for Macs.

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

To install the dependencies on OS X, we can use `Homebrew <http://brew.sh/>`_::

    $ brew update
    $ brew install gsl homebrew/science/hdf5

************
Installation
************

The simplest method of installation is to use PyPI and pip::

    # pip install msprime

This will work in most cases, once the `Requirements`_ have been
satisfied. See below for platform specific build instructions when this
fails.

If you do not have root access to your machine, you can install
``msprime`` into your local Python installation as follows::

    $ pip install msprime --user

To use the ``mspms`` program you must ensure
that the ``~/.local/bin`` directory is in your ``PATH``, or
simply run it using::

    $ ~/.local/bin/mspms

To uninstall ``msprime``, simply run::

    $ pip uninstall msprime

------------------------------
Platform specific installation
------------------------------

This section contains instructions to build on platforms
that require build time flags for GSL and HDF5.

++++++++++++
FreeBSD 10.0
++++++++++++

Install the prerequisitites, and build ``msprime`` as follows::

    # pkg install gsl hdf5-18
    # CFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib pip install msprime

This assumes that root is logged in using a bash shell. For other shells,
different methods are need to set the ``CFLAGS`` and ``LDFLAGS`` environment
variables.

++++
OS X
++++

First, ensure that Homebrew is installed and up-to-date::

    $ brew update

We need to ensure that the version of Python we used is installed via Homebrew
(there can be issues with linking to HDF5 if we use the built-in version of
Python or a version from Anaconda). Therefore, we install Python 3 using
homebrew::

    $ brew install python3
    $ pip3 install --upgrade pip setuptools

The previous step can be skipped if you wish to use your own Python installation,
and already have a working pip.

Now install the dependencies and msprime::

    $ brew install gsl homebrew/science/hdf5
    $ pip3 install msprime

Check if it works::

    $ mspms 10 1 -T

****************
Tested platforms
****************

Msprime is highly portable, and has been successfully built and tested
on the following platforms:

====================    ========        ======          ===========
Operating system        Platform        Python          Compiler
====================    ========        ======          ===========
Debian jessie           x86-64          2.7.9           gcc 4.9.2
Debian jessie           x86-64          3.4.2           gcc 4.9.2
Debian wheezy           i686            2.7.3           gcc 4.7.2
Fedora 20               x86-64          2.7.5           gcc 4.8.3
FreeBSD 10              x86-64          2.7.6           clang 3.4.1
SunOS 5.10              SPARC           3.4.1           gcc 3.4.3
OS X 10.10.5            x86-64          3.4.2           clang 7.0.0
====================    ========        ======          ===========
