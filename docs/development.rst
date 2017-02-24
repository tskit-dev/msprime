.. _sec-development:

=======================
Developer documentation
=======================

If you would like to add some features to ``msprime``, please read the
following. If you think there is anything missing,
please open an `issue <http://github.com/jeromekelleher/msprime/issues>`_ or
`pull request <http://github.com/jeromekelleher/msprime/pulls>`_ on GitHub!

**********
Quickstart
**********

- Make a fork of the msprime repo on `GitHub <github.com/jeromekelleher/msprime>`_
- Clone your fork into a local directory.
- Install the :ref:`basic requirements <sec-requirements>`.
- Install the development requirements using ``pip install -r requirements.txt``.
- Build the lowlevel module by running ``make`` in the project root. If you
  are using Python 2.7, run ``make ext2`` and if you are using Python 3.x,
  run ``make ext3``.
- Run the tests to ensure everything has worked: ``nosetests -vs``. These should
  all pass.
- Make your changes in a local branch, and open a pull request on GitHub when you
  are ready. Please make sure that (a) the tests pass before you open the PR; and
  (b) your code passes PEP8 checks (see below for a git commit hook below to ensure this
  happens automatically) the before opening the PR.

********
Overview
********

There are three main parts of ``msprime``, in increasing order of complexity:

1. High-level Python. The Python-API and command line interface tools are all defined
   in the ``msprime`` directory.

2. C library. The underlying high-performance C code is written as a standalone library.
   All of the code for this library is in the ``lib`` directory.

3. Low-level Python-C interface. The interface between the Python and C code is the
   ``_msprimemodule.c`` file, which defines the ``_msprime`` module.


Each of these aspects has its own coding conventions and development tools, which are
documented in the following sections.

*****************
High-level Python
*****************

+++++++++++
Conventions
+++++++++++

All Python code follows the `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ style
guide, and is checked using the `flake8 <http://flake8.pycqa.org/en/latest/>`_ as
part of the continuous integration tests. In particular, lines must be no longer than
89 characters.

To avoid failing CI tests, it's a good idea to install a local `commit hook
<http://git-scm.com/book/gr/v2/Customizing-Git-Git-Hooks>`_ to automatically check
that code conforms to PEP8 before committing. Adding this to your ``.git/hooks/pre-commit``
should do the trick:

.. code-block:: bash

    # Run flake8 to check for lint errors.
    exec flake8 --max-line-length 89 setup.py msprime tests


+++++++++
Packaging
+++++++++

``msprime`` is packaged and distributed as Python module, and follows the current
`best-practices <http://packaging.python.org>`_ advocated by the
`Python Packaging Authority <http://pypa.io/en/latest/>`_. The primary means of
distribution is though `PyPI <http://pypi.python.org/pypi/msprime>`_, which provides the
canonical source for each release.

A package for `conda <http://conda.io/docs/>`_ is also available on
`BioConda <bioconda.github.io/recipes/msprime/README.html>`_.

+++++
Tests
+++++

The tests for the high-level code are in the ``tests`` directory, and run using
`nose <http://nose.readthedocs.io/en/latest/>`_. A lot of the simulation and basic
tests are contained in the `tests/test_highlevel.py`` file, but more recently
smaller test files with more focussed tests are preferred (e.g., ``test_vcf.py``,
``test_demography.py``).

All new code must have high test coverage, which will be checked as part of the
continuous integration tests by `CodeCov <https://codecov.io/gh/jeromekelleher/msprime/>`_.

+++++++++++++++++++++++++++++++++
Interfacing with low-level module
+++++++++++++++++++++++++++++++++

Much of the high-level Python code only exists to provide a simpler interface to
the low-level ``_msprime`` module. As such, many objects (such as ``TreeSequence``)
are really just a shallow layer on top of the corresponding low-level object.
The convention here is to keep a reference to the low-level object via
an instance variable such as ``self._ll_tree_sequence``.

+++++++++++++++++++++++
Command line interfaces
+++++++++++++++++++++++

The command line interfaces for ``msprime`` are defined in the ``msprime/cli.py`` file.
Each CLI has a single entry point (e.g. ``msp_main``) which is invoked to run the
program. These entry points are registered with ``setuptools`` using the
``console_scripts`` argument in ``setup.py``, which allows them to be deployed as
first-class executable programs in a cross-platform manner.

There are simple scripts in the root of the project (e.g., ``msp_dev.py``)
which are used for development. For example, to run the development version of
``mspms``, use ``python mspms_dev.py``.

*********
C Library
*********


******************
Python C Interface
******************

++++++++++++++++++++++++
Testing for memory leaks
++++++++++++++++++++++++

The Python C API can be subtle, and it is easy to get the reference counting wrong.
The ``stress_lowlevel.py`` script makes it easier to track down memory leaks
when they do occur. The script runs the unit tests in a loop, and outputs
memory usage statistics.

*****************
Statistical tests
*****************





