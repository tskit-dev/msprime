.. _sec_installation:

############
Installation
############

There are two basic options for installing ``msprime``: either through
pre-built binary packages using :ref:`sec_installation_conda` or
by compiling locally using :ref:`sec_installation_pip`. We recommend using ``conda``
for most users, although ``pip`` can be more convenient in certain cases.

.. _sec_installation_conda:

=====
conda
=====

Pre-built binary packages for ``msprime`` are available through
`conda <https://conda.io/docs/>`_, and built using `conda-forge <https://conda-forge.org/>`_.
Packages for Python 2.7, 3.5 and 3.6 are available for Linux and OSX.
For Windows, only Python 3.5 and 3.6 are supported.

We strongly recommend using Python 3.

***********
Quick Start
***********

1. Install ``conda`` using `miniconda <https://conda.io/miniconda.html>`_.
   Make sure you follow the instructions to fully activate your ``conda``
   installation!
2. Set up the `conda-forge channel <https://conda-forge.org/>`_ using
   ``conda config --add channels conda-forge``.
3. Install msprime: ``conda install msprime``.
4. Try it out: ``msp --version``.


There are several different ways to obtain ``conda``. Please see the
`anaconda installation documentation <https://docs.anaconda.com/anaconda/install/>`_
for full details.


************
Full example
************

In this example we create a
`conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_
and install ``msprime`` into it on an OSX machine.
We assume that ``conda`` has been installed  and bash shell is being used (Windows users will need to modify the
commands slightly).

.. code-block:: none

    $ conda create --name msprime-env
    Solving environment: done

    ## Package Plan ##

      environment location: /Users/jk/miniconda3/envs/msprime-env


    Proceed ([y]/n)? y

    Preparing transaction: done
    Verifying transaction: done
    Executing transaction: done
    #
    # To activate this environment, use
    #
    #     $ conda activate msprime-env
    #
    # To deactivate an active environment, use
    #
    #     $ conda deactivate

    $ source activate msprime-env
    (msprime-env) $ conda install -c conda-forge msprime
    Solving environment: done

    ## Package Plan ##

      environment location: /Users/jk/miniconda3/envs/msprime-env

      added / updated specs:
        - msprime


    The following NEW packages will be INSTALLED:

        ca-certificates: 2018.1.18-0           conda-forge
        certifi:         2018.1.18-py36_0      conda-forge
        gsl:             1.16-0                conda-forge
        hdf5:            1.10.1-2              conda-forge
        intel-openmp:    2018.0.0-h8158457_8
        libgfortran:     3.0.1-h93005f0_2
        mkl:             2018.0.1-hfbd8650_4
        msprime:         0.5.0b2-py36_3        conda-forge
        ncurses:         5.9-10                conda-forge
        numpy:           1.14.1-py36h8a80b8c_1
        openssl:         1.0.2n-0              conda-forge
        pip:             9.0.1-py36_1          conda-forge
        pyparsing:       2.2.0-py36_0          conda-forge
        python:          3.6.4-0               conda-forge
        readline:        7.0-0                 conda-forge
        setuptools:      38.5.1-py36_0         conda-forge
        six:             1.11.0-py36_1         conda-forge
        sqlite:          3.20.1-2              conda-forge
        svgwrite:        1.1.12-py_0           conda-forge
        tk:              8.6.7-0               conda-forge
        wheel:           0.30.0-py36_2         conda-forge
        xz:              5.2.3-0               conda-forge
        zlib:            1.2.11-0              conda-forge

    Proceed ([y]/n)? y

    Preparing transaction: done
    Verifying transaction: done
    Executing transaction: done
    (msprime-env) $ python
    Python 3.6.4 | packaged by conda-forge | (default, Dec 23 2017, 16:54:01)
    [GCC 4.2.1 Compatible Apple LLVM 6.1.0 (clang-602.0.53)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import msprime
    >>> tree = msprime.simulate(5).first()
    >>> print(tree.draw(format="unicode"))
      8
    ┏━┻━━┓
    ┃    7
    ┃  ┏━┻━┓
    ┃  ┃   6
    ┃  ┃  ┏┻┓
    ┃  5  ┃ ┃
    ┃ ┏┻┓ ┃ ┃
    1 0 4 2 3


Please see the `conda documentation <https://conda.io/docs/index.html>`_ for
more details on managing packages and environments.


.. _sec_installation_pip:

===
pip
===

Installing using ``pip`` is more flexible than ``conda`` as it
can support more versions of Python, and the locations of the
various dependencies can be specified.

***********
Quick Start
***********

To install and run ``msprime`` on a fresh Ubuntu 15.10 installation, do the
following:

.. code-block:: bash

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

.. _sec-requirements:

*************
Requirements
*************

Msprime requires Python 2.7+ (Python 3 versions are fully supported from
3.1 onwards), the `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_,
and `HDF5 <https://www.hdfgroup.org/HDF5/>`_ version 1.8 or later. These
packages are available for all major platforms. For example, to install on
Debian/Ubuntu use (as root)::

    $ apt-get install python-dev libgsl0-dev libhdf5-serial-dev pkg-config

For Redhat/Fedora use::

    $ yum install gsl-devel hdf5-devel

On FreeBSD we can use ``pkg`` to install the requirements::

    $ pkg install gsl hdf5-18

To install the dependencies on OS X, we can use `Homebrew <http://brew.sh/>`_::

    $ brew update
    $ brew install gsl homebrew/science/hdf5

************
Installation
************

The simplest method of installation is to use PyPI and pip::

    $ pip install msprime

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

    $ pkg install gsl hdf5-18
    $ CFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib pip install msprime

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
