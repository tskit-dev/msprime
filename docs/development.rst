.. _sec-development:

=======================
Developer documentation
=======================

If you would like to add some features to ``msprime``, please read the
following. If you think there is anything missing,
please open an `issue <http://github.com/tskit-dev/msprime/issues>`_ or
`pull request <http://github.com/tskit-dev/msprime/pulls>`_ on GitHub!

**********
Quickstart
**********

- Make a fork of the msprime repo on `GitHub <github.com/tskit-dev/msprime>`_
- Clone your fork into a local directory, making sure that the **submodules
  are correctly initialised**::

  $ git clone git@github.com:tskit-dev/msprime.git --recurse-submodules

  For an already checked out repo, the submodules can be initialised using::

  $ git submodule update --init --recursive

- Install the :ref:`basic requirements <sec_installation_system_requirements>`.
- Install the Python development requirements using ``pip install -r requirements/development.txt``.
- Build the low level module by running ``make`` in the project root. If you
  are using Python 2.7, run ``make ext2`` and if you are using Python 3.x,
  run ``make ext3``. Python 3.x is the default if no arguments are provided.
- Run the tests to ensure everything has worked: ``python -m nose -vs``. These should
  all pass.
- Make your changes in a local branch, and open a pull request on GitHub when you
  are ready. Please make sure that (a) the tests pass before you open the PR; and
  (b) your code passes PEP8 checks (see below for a git commit hook to ensure this
  happens automatically) before opening the PR.

****************************
Continuous integration tests
****************************

Three different continuous integration providers are used, which run different
combinations of tests on different platforms:

1. `Travis CI <https://travis-ci.org/>`_ runs tests on Linux and OSX using the
   `Conda <https://conda.io/docs/>`__ infrastructure for the system level
   requirements. All supported versions of Python are tested here.

2. `CircleCI <https://circleci.com/>`_ Runs all Python tests using the apt-get
   infrastructure for system requirements. Additionally, the low-level tests
   are run, coverage statistics calculated using `CodeCov <https://codecov.io/gh>`__,
   and the documentation built.

3. `AppVeyor <https://www.appveyor.com/>`_ Runs Python tests on 32 and 64 bit
   Windows using conda.

+++++++++++++++++++++++++++++++++++++++++++++++++
Running tests on multiple Python versions locally
+++++++++++++++++++++++++++++++++++++++++++++++++

On `Travis CI <https://travis-ci.org/>`_ all supported Python versions are tested.
If you'd like to test multiple versions locally, you can use `tox`:

.. code-block:: bash

    echo \
    '[tox]
    envlist = py27,py35
    [testenv]
    deps= -rrequirements/development.txt
    commands=nosetests' > tox.ini && tox

Note that if the `requirements/development.txt` have been updated since
initially running `tox`, you may need to `recreate them <http://tox.readthedocs.io/en/latest/example/basic.html#forcing-re-creation-of-virtual-environments>`_:

them:

.. code-block:: bash

    tox --recreate -e py27,py35

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

Throughout this document, we assume that the ``msprime`` package is built and
run locally *within* the project directory. That is, ``msprime`` is *not* installed
into the Python installation using ``pip install -e`` or setuptools `development
mode <http://setuptools.readthedocs.io/en/latest/setuptools.html#id23>`_. Please
ensure that you build the low-level module using (e.g.) ``make ext3`` and that
the shared object file is in the project root.

+++++++++++
Conventions
+++++++++++

All Python code follows the `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ style
guide, and is checked using the `flake8 <http://flake8.pycqa.org/en/latest/>`_  tool as
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
`conda-forge <https://github.com/conda-forge/msprime-feedstock>`_.

+++++
Tests
+++++

The tests for the high-level code are in the ``tests`` directory, and run using
`nose <http://nose.readthedocs.io/en/latest/>`_. A lot of the simulation and basic
tests are contained in the ``tests/test_highlevel.py`` file, but more recently
smaller test files with more focussed tests are preferred (e.g., ``test_vcf.py``,
``test_demography.py``).

All new code must have high test coverage, which will be checked as part of the
continuous integration tests by `CodeCov <https://codecov.io/gh/tskit-dev/msprime/>`_.

+++++++++++++++++++++++++++++++++
Interfacing with low-level module
+++++++++++++++++++++++++++++++++

Much of the high-level Python code only exists to provide a simpler interface to
the low-level ``_msprime`` module. As such, many objects (such as ``RecombinationMap``)
are really just a shallow layer on top of the corresponding low-level object.
The convention here is to keep a reference to the low-level object via
a private instance variable such as ``self._ll_recombination_map``.

+++++++++++++++++++++++
Command line interfaces
+++++++++++++++++++++++

The command line interfaces for ``msprime`` are defined in the ``msprime/cli.py`` file.
Each CLI has a single entry point (e.g. ``msp_main``) which is invoked to run the
program. These entry points are registered with ``setuptools`` using the
``console_scripts`` argument in ``setup.py``, which allows them to be deployed as
first-class executable programs in a cross-platform manner.

There are simple scripts in the root of the project (currently: ``msp_dev.py``,
``mspms_dev.py``) which are used for development. For example, to run the
development version of ``mspms`` use ``python mspms_dev.py``.

*********
C Library
*********

The low-level code for ``msprime`` is written in C, and is structured as a
standalone library. This code is all contained in the ``lib`` directory.
Although the code is structured as a library, it is not intended to be used
outside of the ``msprime`` project! The interfaces at the C level change
considerably over time, and are deliberately undocumented.

++++++
Basics
++++++

To compile and develop the C code, a few extra development libraries are needed.
`Libconfig <http://www.hyperrealm.com/libconfig/>`_ is used for the development CLI
and `CUnit <http://cunit.sourceforge.net>`_ for unit tests. On Debian/Ubuntu, these
can be installed using

.. code-block:: bash

    $ sudo apt-get install libcunit1-dev libconfig-dev

Compile the code locally run ``make`` in the ``lib`` directory.


+++++++++++++++
Development CLI
+++++++++++++++

When developing the C code, it is usually best to use the development CLI to invoke
the code. This is much simpler than going through the Python interface, and allows
tools such as `valgrind <http://valgrind.org>`_ to be used directly. For example,
when developing new simulation functionality, you should get the basic work done
using the CLI and only move over to the Python API once you are reasonably sure
that the code works properly.

The development CLI is written using `libconfig
<http://www.hyperrealm.com/libconfig/>`_ to parse the simulation parameters
file, and `argtable3 <https://github.com/argtable/argtable3>`_ to parse the
command line arguments. The ``argtable3`` code is included in the source (but
not used in the distributed binaries, since this is strictly a development
tool).

The CLI is run as follows:

.. code-block:: bash

    $ ./main <command> <arguments>

Running the ``main`` program without arguments will print out a summary of the
options.

.. warning

    The development CLI is a tool used to develop the msprime API, and not a
    polished artefact intended for users. There is quite a lot of code left
    over from earlier debugging which might not make immediate sense. Some
    commands may not work as expected, or indeed at all. Please feel free to
    tidy it up if you would like to improve it!

The most important command for simulator development is ``simulate``,
which takes a configuration file as a parameter and writes the resulting
simulation to an output file in the native ``.trees`` format. For example,

.. code-block:: bash

    $ ./main simulate dev.cfg -o out.trees

The development configuration file describes the simulation that we want to
run, and uses the
`libconfig syntax <http://www.hyperrealm.com/libconfig/libconfig_manual.html#Configuration-Files>`_.
An example is given in the file ``dev.cfg`` which should have sufficient documentation
to be self-explanatory.

.. warning

    It is important to note that all values in the low-level C code are in
    scaled coalescent units. The high-level Python API defines values in units
    of generations, but for the C code all time is measured in coalescent units.

++++++++++
Unit Tests
++++++++++

The C-library has an extensive suite of unit tests written using
`CUnit <http://cunit.sourceforge.net>`_. These tests aim to establish that the
low-level APIs work correctly over a variety of inputs, and particularly, that
the tests don't result in leaked memory or illegal memory accesses. The tests should be
periodically run under valgrind to make sure of this.

Tests are split into ``simulation_tests`` which covers functionality specific to the
simulation logic, and ``tests`` which covers everything else. To run all the tests
in a given suite, type ``./tests`` or ``./simulation_tests``.
To run a specific test, provide this test name as a command line argument,
e.g.:

.. code-block:: bash

    $ ./simulation_tests fenwick_tree


While 100% test coverage is not feasible for C code, we aim to cover all code
that can be reached. (Some classes of error such as malloc failures
and IO errors are difficult to simulate in C.) Code coverage statistics are
automatically tracked using `CodeCov <https://codecov.io/gh/tskit-dev/msprime/>`_.

++++++++++++++++++
Coding conventions
++++++++++++++++++

The code is written using the `C99 <https://en.wikipedia.org/wiki/C99>`_ standard. All
variable declarations should be done at the start of a function, and functions
kept short and simple where at all possible.

No global or module level variables are used for production code.

The code is organised following object-oriented principles. Each 'class' is defined using
a struct, which encapsulates all the data it requires. Every 'method' on this class
is then a function that takes this struct as its first parameter. Each class has
an ``alloc`` method, which is responsible for allocating memory and a ``free`` method
which frees all memory used by the object. For example, the
`Fenwick tree <https://en.wikipedia.org/wiki/Fenwick_tree>`_ class is defined as
follows:

.. code-block:: C

    typedef struct {
        size_t size;
        size_t log_size;
        int64_t *tree;
        int64_t *values;
    } fenwick_t;

    int fenwick_alloc(fenwick_t *self, size_t initial_size);
    int fenwick_free(fenwick_t *self);
    int64_t fenwick_get_total(fenwick_t *self);

This defines the ``fenwick_t`` struct, and alloc and free methods and a method
to return the total of the tree. Note that we follow the Python convention
and use ``self`` to refer to the current instance.

Most objects also provide a ``print_state`` method, which is useful for
debugging.

This object-oriented structure means that the code is fully thread safe.


++++++++++++++
Error handling
++++++++++++++

A critical element of producing reliable C programs is consistent error handling
and checking of return values. All return values **must** be checked! In msprime,
all functions (except the most trivial accessors) return an integer to indicate
success or failure. Any negative value is an error, and must be handled accordingly.
The following pattern is canonical:

.. code-block:: C

        ret = msp_do_something(self, argument);
        if (ret != 0) {
            goto out;
        }
        // rest of function
    out:
        return ret;

Here we test the return value of ``msp_do_something`` and if it is non-zero,
abort the function and return this same value from the current function. This
is a bit like throwing an exception in higher-level languages, but discipline
is required to ensure that the error codes are propagated back to the original
caller correctly.

Particular care must be taken in functions that allocate memory, because
we must ensure that this memory is freed in all possible success and
failure scenarios. The following pattern is used throughout for this purpose:

.. code-block:: C

        double x = NULL;

        x = malloc(n * sizeof(double));
        if (x == NULL) {
            ret = MSP_ERR_NO_MEMORY;
            goto out;
        }
        // rest of function
    out:
        if (x != NULL) {
            free(x);
        }
        return ret;


It is vital here that ``x`` is initialised to ``NULL`` so that we are guaranteed
correct behaviour in all cases. For this reason, the convention is to declare all
pointer variables on a single line and to initialise them to ``NULL`` as part
of the declaration.

Error codes are defined in ``err.h``, and these can be translated into a
message using ``msp_strerror(err)``.

++++++++++++++++
Running valgrind
++++++++++++++++

Valgrind is an essential development tool, and is used extensively. (Being able
to run valgrind was one of the motivating factors in the C-library architecture.
It is difficult to run valgrind on a Python extension module, and so the simplest
way to ensure that the low-level code is memory-tight is to separate it out
into an independent library.)

Any new C unit tests that are written should be verified using valgrind to
ensure that no memory is leaked. The entire test suite should be run
through valgrind periodically also to detect any leaks or illegal
memory accesses that have been overlooked.

******************
Python C Interface
******************

++++++++
Overview
++++++++

The Python C interface is written using the
`Python C API <https://docs.python.org/3.6/c-api/>`_ and the code is in the
``_msprimemodule.c`` file. When compiled, this produces the ``_msprime`` module,
which is imported by the high-level module. The low-level Python module is
not intended to be used directly and may change arbitrarily over time.

The usual pattern in the low-level Python API is to define a Python class
which corresponds to a given "class" in the C API. For example, we define
a ``RecombinationMap`` class, which is essentially a thin wrapper around the
``recomb_map_t`` type from the C library.

The ``_msprimemodule.c`` file follows the standard conventions given in the
`Python documentation <https://docs.python.org/3.6/extending/index.html>`_.


+++++++++
Compiling
+++++++++

The ``setup.py`` file descibes the requirements for the low-level ``_msprime``
module and how it is built from source. To build the module so that it is available
for use in the current working directory, run

.. code-block:: bash

    $ python setup.py build_ext --inplace

A development Makefile is also provided in the project root, so that running
``make ext2`` or ``make ext3`` should build the extension module for either
Python 2 or Python 3.

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

To ensure that ``msprime`` is simulating the correct process we run many statistical
tests. Since these tests are quite expensive (taking some hours to run) and
difficult to automatically validate, they are not run as part of CI but instead
as a pre-release sanity check. They are also very useful to run when developing
new simulation functionality, as subtle statistical bugs can easily slip in
unnoticed.

The statistical tests are all run via the ``verification.py`` script in the project root.
The script has some extra dependencies listed in the ``requirements/verification.txt``,
which can be installed using ``pip install -r`` or ``conda install --file``. Run
this script using:

.. code-block:: bash

    $ python verification.py


The statistical tests depend on compiled programs in the ``data`` directory.
This includes a customised version of ``ms`` and a locally compiled version of
`scrm <https://scrm.github.io/>`_. These programs must be compiled before
running the statistical tests, and can be built by running ``make`` in the
``data`` directory. If this is successful, there should be several binaries
like ``ms`` and ``ms_summary_stats`` present in the ``data``
directory.

The ``verification.py`` script contains lots of different tests, each one
identified by a particular "key". To run all the tests, run the script without
any arguments. To run some specific tests, provide the required keys as command
line arguments.

Many of the tests involve creating an ``ms`` command line, running it
line on ``ms`` and ``msprime`` and comparing the statistical properties of the
results. The output of each test is a series of plots, written to a directory
named after test. For example, results for the ``admixture-1-pop2`` test are
written in the ``tmp__NOBACKUP__/admixture-1-pop2/`` directory (the prefix is
not important here and can be changed). The majority of the results are
QQ-plots of the statistics in question comparing ``ms`` and ``msprime``.

There are also several "analytical" tests, which compare the distributions of
values from ``msprime`` with analytical expectations.

*************
Documentation
*************

Documentation is written using `Sphinx <http://www.sphinx-doc.org/en/stable/>`_
and contained in the ``docs`` directory. It is written in the
`reStructuredText <http://docutils.sourceforge.net/rst.html>`_ format and
is deployed automatically to `readthedocs <https://readthedocs.org/>`_. To
build the documentation locally run ``make`` in the ``docs`` directory.
This should build the HTML documentation in ``docs/_build/html/``.
