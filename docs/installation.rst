.. _sec_installation:

############
Installation
############

There are three options for installing ``msprime``: 

1. :ref:`sec_installation_conda`: the recommended options for most users

2. :ref:`sec_installation_pip`: more flexibility and control for advanced users

3. :ref:`sec_linux_container`: to use
   ``msprime`` inside an extremely isolated environment


.. _sec_installation_conda:

=========
Via Conda
=========

Pre-built binary packages for ``msprime`` are available through
`conda <https://conda.io/docs/>`_, and built using `conda-forge <https://conda-forge.org/>`_.
Packages for Python 3.6, 3.7 and 3.8 are available for Linux, OSX and Windows.

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

=======
Via Pip
=======

Installing using ``pip`` is more flexible than ``conda`` as it
can support more versions of Python, and the locations of the
various dependencies can be specified.

.. warning::
    If you installed Python using anaconda/miniconda, please install
    msprime using ``conda``. Link/load time errors will occur if you mix
    a conda installed Python with system libraries.

To install msprime via pip, first
:ref:`sec_installation_system_requirements`, then :ref:`sec_run_pip`.

.. _sec_installation_system_requirements:

*************************
check system requirements
*************************

Msprime has a number of requirements which may or may not already be
installed on your system:

* Python 3.6+ and pip
* `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_ (GSL)
* other requirements depending on the specific platform

To make sure you have Python 3.6+ and pip installed, run

.. code-block:: bash

    $ python3 -m pip --version

to make sure you see Python 3.6 or greater.
If you do not, do installation :ref:`via conda <sec_installation_conda>`.

To install GSL, follow instructions per your platform.

.. glossary::

    Debian/Ubuntu
        ::

            $ apt-get install python-dev libgsl0-dev

    Redhat/Fedora
        ::

            $ yum install gsl-devel

    FreeBSD
        ::

            $ pkg install gsl

    OS X
        We recommend using :ref:`sec_installation_conda` to install ``msprime`` on OS X.
        However, it is also possible to install using `Homebrew <http://brew.sh/>`_:
        ::

            $ brew update
            $ brew install gsl

    Windows
        Use :ref:`sec_installation_conda`, do not install via ``pip`` on Windows.


There may be additional requirements depending on the specific platform. The
next instructions running pip might not succeed. Depending your platform, you
might be able to determine missing requirements. If not, install :ref:`via
conda <sec_installation_conda>`.


.. _sec_run_pip:

***************
run pip install
***************

We can install ``msprime`` easily using pip::

    $ python3 -m pip install msprime

(It is generally better to use ``python3 -m pip`` rather than call ``pip``
directly since this allows you to control which installation of Python the
package is installed to.)

If you do not have root access to your machine, you can install
``msprime`` into your local Python installation as follows::

    $ python3 -m pip install msprime --user

To use the ``mspms`` program you must ensure
that the ``~/.local/bin`` directory is in your ``PATH``, or
simply run it using::

    $ ~/.local/bin/mspms

To uninstall ``msprime``, simply run::

    $ python3 -m pip uninstall msprime


.. _sec_linux_container:

=============
Via Container 
=============

An `open container <https://opencontainers.org/>`_ image (aka docker image) is built on
`Dockerhub <https://hub.docker.com/r/tskit/msprime>`_ for each release of
msprime. Each image is `tagged <https://hub.docker.com/r/tskit/msprime/tags>`_
with the corresponding release. For example, for msprime release 0.7.5, the
corresponding image tag is ``tskit/msprime:0.7.5``.

To run a container, you can use `docker <https://www.docker.com/>`_,
`Singularity <https://sylabs.io/singularity/>`_,
`podman <https://podman.io/>`_ or similar tools supporting docker images.

******
docker
******

`docker` requires root privilege to run a container:

.. code-block:: bash

    $ sudo docker run -it tskit/msprime:<release> mspms 10 1 -T

******
podman
******

podman can run an msprime container without root privilege:

.. code-block:: bash

    $ podman run -it docker.io/tskit/msprime:<release> mspms 10 1 -T

***********
Singularity
***********

A docker image can also be converted to a Singularity container and
then run without root privilege:

.. code-block:: bash

    $ singularity pull docker://tskit/msprime:<release> msprime-<release>.simg
    $ singularity exec msprime-<release>.simg mspms 10 1 -T

.. note::

  It is possible that your current environment may conflict with the environment in the singularity container.
  There are two workarounds:

  1.  Ignore your home with the conflicting environment with ``--contain`` or ``-H </new/path/to/home> -e``

  .. code-block:: bash

      $ singularity shell --contain msprime-release-0.7.3.simg
      Singularity: Invoking an interactive shell within container...

      Singularity msprime-release-0.7.3.simg:~> python3
      Python 3.6.8 (default, Jan 14 2019, 11:02:34)
      [GCC 8.0.1 20180414 (experimental) [trunk revision 259383]] on linux
      Type "help", "copyright", "credits" or "license" for more information.
      >>> import msprime
      >>>

  or use a different path as your home that does not have a conflicting environment

  .. code-block:: bash

      $ singularity shell -H </new/path/to/home> -e msprime-release-0.7.3.simg
      Singularity: Invoking an interactive shell within container...

      Singularity msprime-release-0.7.3.simg:~/cnn_classify_demography> python3
      Python 3.6.8 (default, Jan 14 2019, 11:02:34)
      [GCC 8.0.1 20180414 (experimental) [trunk revision 259383]] on linux
      Type "help", "copyright", "credits" or "license" for more information.
      >>> import msprime
      >>>

  2. In python get rid of your local path

  .. code-block:: bash

      $ singularity shell msprime-release-0.7.3.simg
      Singularity: Invoking an interactive shell within container...

      Singularity msprime-release-0.7.3.simg:~> python3
      Python 3.6.8 (default, Jan 14 2019, 11:02:34)
      [GCC 8.0.1 20180414 (experimental) [trunk revision 259383]] on linux
      Type "help", "copyright", "credits" or "license" for more information.
      >>> import sys
      >>> for _path in sys.path:
      ...     if ".local" in _path:
      ...             sys.path.remove(_path)
      ...
      >>> import msprime
      >>>


For more information on Singularity, see https://www.sylabs.io/guides/3.6/user-guide/

