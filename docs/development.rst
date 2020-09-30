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

- Make a fork of the msprime repo on `GitHub <http://github.com/tskit-dev/msprime>`_
- Clone your fork into a local directory, making sure that the **submodules
  are correctly initialised**::

  $ git clone git@github.com:tskit-dev/msprime.git --recurse-submodules

  For an already checked out repo, the submodules can be initialised using::

  $ git submodule update --init --recursive

- Install the :ref:`basic requirements <sec_installation_system_requirements>`.
- Install the Python development requirements using ``pip install -r requirements/development.txt``.
- Build the low level module by running ``make`` in the project root.
- Run the tests to ensure everything has worked: ``python3 -m pytest``. These should
  all pass.
- Install the pre-commit checks: ``pre-commit install``
- Make your changes in a local branch. On each commit a `pre-commit hook
  <https://pre-commit.com/>`_  will run
  checks for code style and common problems.
  Sometimes these will report "files were modified by this hook" ``git add``
  and ``git commit --amend`` will update the commit with the automatically modified
  version.
  The modifications made are for consistency, code readability and designed to
  minimise merge conflicts. They are guaranteed not to modify the functionality of the
  code. To run the checks without committing use ``pre-commit run``. To bypass
  the checks (to save or get feedback on work-in-progress) use ``git commit
  --no-verify``
- If you have modifed the C code then
  ``clang-format -i lib/tests/* lib/!(avl).{c,h}`` will format the code to
  satisfy CI checks.
- When ready open a pull request on GitHub. Please make sure that the tests pass before
  you open the PR, unless you want to ask the community for help with a failing test.
- See the `tskit documentation <https://tskit.readthedocs.io/en/latest/development.html#github-workflow>`_
  for more details on the recommended GitHub workflow.

****************************
Continuous integration tests
****************************

Three different continuous integration providers are used, which run different
combinations of tests on different platforms:

1. A `Github action <https://help.github.com/en/actions>`_ runs `pre-commit
   <https://pre-commit.com/>`_ to run a variety of code style and quality checks.

2. `Travis CI <https://travis-ci.org/>`_ runs tests on Linux and OSX using the
   `Conda <https://conda.io/docs/>`__ infrastructure for the system level
   requirements. All supported versions of Python are tested here.

3. `CircleCI <https://circleci.com/>`_ Runs all Python tests using the apt-get
   infrastructure for system requirements. Additionally, the low-level tests
   are run, coverage statistics calculated using `CodeCov <https://codecov.io/gh>`__,
   and the documentation built.

4. `AppVeyor <https://www.appveyor.com/>`_ Runs Python tests on Windows using conda.

********
Overview
********

There are three main parts of ``msprime``, in increasing order of complexity:

1. High-level Python. The Python-API and command line interface tools are all defined
   in the ``msprime`` directory.

2. C library. The underlying high-performance C code is written as a standalone library.
   All of the code for this library is in the ``lib`` directory.

3. Low-level Python-C interface. The interface between the Python and C code is the
   ``msprime/_msprimemodule.c`` file, which defines the ``msprime._msprime`` module.


Each of these aspects has its own coding conventions and development tools, which are
documented in the following sections.

*****************
High-level Python
*****************

Throughout this document, we assume that the ``msprime`` package is built and
run locally *within* the project directory. That is, ``msprime`` is *not* installed
into the Python installation using ``pip install -e`` or setuptools `development
mode <http://setuptools.readthedocs.io/en/latest/setuptools.html#id23>`_. Please
ensure that you build the low-level module using (e.g.) ``make`` and that
the shared object file is in the ``msprime`` directory. This will have a name
like ``_msprime.cpython-38-x86_64-linux-gnu.so``, depending on your platform
and Python version.

+++++++++++
Conventions
+++++++++++

All Python code follows the `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ style
guide, and is checked using the `flake8 <http://flake8.pycqa.org/en/latest/>`_  tool as
part of the continuous integration tests. `Black <https://github.com/psf/black>`_ is
used as part of the pre-commit hook for python code style and formatting.

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
`pytest <https://docs.pytest.org/en/stable/>`_. A lot of the simulation and basic
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
development version of ``mspms`` use ``python3 mspms_dev.py``.

*********
C Library
*********

The low-level code for ``msprime`` is written in C, and is structured as a
standalone library. This code is all contained in the ``lib`` directory.
Although the code is structured as a library, it is not intended to be used
outside of the ``msprime`` project! The interfaces at the C level change
considerably over time, and are deliberately undocumented.

+++++++++
Toolchain
+++++++++

To compile and develop the C code, a few extra development libraries are needed.
`Libconfig <http://www.hyperrealm.com/libconfig/>`_ is used for the development CLI
and `CUnit <http://cunit.sourceforge.net>`_ for unit tests. We use the
`meson <https://mesonbuild.com>`_ build system in conjunction with `ninja-build
<ninja-build.org>`_ to to compile the unit tests and
development CLI. On Debian/Ubuntu, these can be installed using

.. code-block:: bash

    $ sudo apt-get install libcunit1-dev libconfig-dev ninja-build

Meson is best installed via ``pip``:

.. code-block:: bash

    $ python3 -m pip install meson --user

On macOS rather than use ``apt-get`` for installation of these requirements
a combination of ``homebrew`` and ``pip`` can be used (working as of 2020-01-15).

.. code-block:: bash

    $ brew install cunit
    $ python3 -m pip install meson --user
    $ python3 -m pip install ninja --user

On macOS, conda builds are generally done using ``clang`` packages that are kept up to date:

.. code-block:: bash

    $ conda install clang_osx-64  clangxx_osx-64

In order to make sure that these compilers work correctly (*e.g.*, so that they can find
other dependencies installed via ``conda``), you need to compile ``msprime`` with this command
on versions of macOS older than "Mojave":

.. code-block:: bash

    $ CONDA_BUILD_SYSROOT=/ python3 setup.py build_ext -i

On more recent macOS releases, you may omit the ``CONDA_BUILD_SYSROOT`` prefix.

.. note::

   The use of the C toolchain on macOS is a moving target.  The above advice
   was written on 23 January, 2020 and was validated by a few ``msprime`` contributors.
   Caveat emptor, etc..

+++++++++
Compiling
+++++++++

Meson keeps all compiled binaries in a build directory (this has many advantages
such as allowing multiple builds with different options to coexist). It depends on
a ``meson.build`` file which is in the ``lib`` directory. To set up the initial build
directory, run

.. code-block:: bash

    $ cd lib
    $ meson build

The easiest way to compile the :ref:`sec-development-c-unit-tests`
is to run ``ninja -C build``. (Alternatively,
you can ``cd`` into the ``build`` directory and run ``ninja``). All the
compiled binaries are then in the ``build`` directory, so to run, for example, the
``test_ancestry`` unit tests, use ``./build/test_ancestry``. A handy shortcut
to compile the code and run all the unit tests is:

.. code-block:: bash

    $ ninja -C build test

The `mesonic <www.vim.org/scripts/script.php?script_id=5378>`_ plugin for vim
simplifies this process and allows code to be compiled seamlessly within the
editor.

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
tool). The source code is in ``dev-tools/dev-cli.c``.

After building, the CLI is run as follows:

.. code-block:: bash

    $ ./build/dev-cli <command> <arguments>

Running the ``dev-cli`` program without arguments will print out a summary of the
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

    $ ./build/dev-cli simulate dev-tools/example.cfg -o out.trees

The development configuration file describes the simulation that we want to
run, and uses the
`libconfig syntax <http://www.hyperrealm.com/libconfig/libconfig_manual.html#Configuration-Files>`_.
An example is given in the file ``dev-tools/example.cfg`` which should have sufficient documentation
to be self-explanatory.

.. warning

    It is important to note that all values in the low-level C code are in
    scaled coalescent units. The high-level Python API defines values in units
    of generations, but for the C code all time is measured in coalescent units.


.. _sec-development-c-unit-tests:

++++++++++
Unit Tests
++++++++++

The C-library has an extensive suite of unit tests written using
`CUnit <http://cunit.sourceforge.net>`_. These tests aim to establish that the
low-level APIs work correctly over a variety of inputs, and particularly, that
the tests don't result in leaked memory or illegal memory accesses. The tests should be
periodically run under valgrind to make sure of this.

Tests are defined in the ``tests`` directory, roughly split into suites
defined in different files. For example, the tests associated with Fenwick
trees are defined in the ``tests/tests_fenwick.c`` file. To run all the
tests in this suite, use run using ``./build/test_fenwick``.
To run a specific test in a particular suite, provide the name of the
test name as a command line argument, e.g.:

.. code-block:: bash

    $ ./build/test_fenwick test_fenwick_expand

While 100% test coverage is not feasible for C code, we aim to cover all code
that can be reached. (Some classes of error such as malloc failures
and IO errors are difficult to simulate in C.) Code coverage statistics are
automatically tracked using `CodeCov <https://codecov.io/gh/tskit-dev/msprime/>`_.

++++++++++
Code Style
++++++++++

C code is formatted using
`clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_
with a custom configuration.
To ensure that your code is correctly formatted, you can run

.. code-block:: bash

   make clang-format

in the project root before submitting a pull request. Alternatively,
you can run ``clang-format -i *.[c,h]`` in the ``lib`` directory.

Vim users may find the
`vim-clang-format <https://github.com/rhysd/vim-clang-format>`_
plugin useful for automatically formatting code.


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
        double *tree;
        double *values;
    } fenwick_t;

    int fenwick_alloc(fenwick_t *self, size_t initial_size);
    int fenwick_free(fenwick_t *self);
    double fenwick_get_total(fenwick_t *self);

This defines the ``fenwick_t`` struct, and alloc and free methods and a method
to return the total of the tree. Note that we follow the Python convention
and use ``self`` to refer to the current instance.

Most objects also provide a ``print_state`` method, which is useful for
debugging.

Please see the documentation for the `tskit C API
<https://tskit.readthedocs.io/en/stable/c-api.html#sec-c-api-overview-structure>`_
for more details on the how APIs are structured.

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

        double *x = NULL;

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

The Python C interface is written using the
`Python C API <https://docs.python.org/3.6/c-api/>`_ and the code is in the
``msprime/_msprimemodule.c`` file. When compiled, this produces the
``msprime._msprime`` module,
which is imported by the high-level module. The low-level Python module is
not intended to be used directly and may change arbitrarily over time.

The conventions used within the low-level module here closely follow
those in ``tskit``; please see the
`documentation
<https://tskit.readthedocs.io/en/stable/development.html#python-c-interface>`_
for more information.

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

    $ python3 verification.py


The statistical tests depend on compiled programs in the ``data`` directory.
This includes a customised version of ``ms`` and a locally compiled version of
`scrm <https://scrm.github.io/>`_. These programs must be compiled before
running the statistical tests, and can be built by running ``make`` in the
``data`` directory. If this is successful, there should be several binaries
like ``ms`` and ``ms_summary_stats`` present in the ``data``
directory.

Please the comments at the top of the ``verification.py`` script for details
on how to write and run these tests.


************
Benchmarking
************

Benchmarks to measure performance are in the ``benchmarks`` folder and are run using
`airspeed velocity <https://asv.readthedocs.io/en/stable/index.html>`_.
An automated system runs the benchmarks on each push to the main branch and uploads
the results to `this github pages site` <https://tskit-dev.github.io/msprime-asv>_.
These benchmarks can also be run locally to compare your branch with the main branch.
Your changes must be in a commit to be measured. To run the benchmarks::

    asv run asv run HEAD...main~1

This will run the benchmarks for the latest main branch commit and all commits on
your current branch (the syntax for choosing commits is the same as ``git log``).
The following commands then make a browsable report (link given in output of
the command)::

    asv publish
    asv preview

Note the following tips:

- Specifying the range of commits to run uses the same syntax as git log.
  For example, to run for a single commit, use ``asv run 88fbbc33^!``

- Be careful when running ``asv dev`` or using ``python=same`` as
  this can use the *installed* version of msprime rather than the local
  development version. This can lead to confusing results! When tuning
  benchmarks it's better to commit often and use (e.g.)
  ``asv run HEAD^! --show-stderr -b Hudson.time_large_sample_size``.


****************
Containerization
****************

To run msprime in a container, see the
:ref:`installation instructions as Linux container <sec_linux_container>`.

You can use ``docker`` to locally build an image, but it requires root access:

.. code-block:: bash

    $ sudo docker build -t tskit/msprime .

`podman <https://podman.io/>`_ can build and run images without root privilege.

.. code-block:: bash

    $ podman build -t tskit/msprime .


*************
Documentation
*************

Documentation is written using `Sphinx <http://www.sphinx-doc.org/en/stable/>`_
and contained in the ``docs`` directory. It is written in the
`reStructuredText <http://docutils.sourceforge.net/rst.html>`_ format and
is deployed automatically to `readthedocs <https://readthedocs.org/>`_. To
build the documentation locally run ``make`` in the ``docs`` directory.
This should build the HTML documentation in ``docs/_build/html/``.


***************
Troubleshooting
***************

- If ``make`` is giving you strange errors, or if tests are failing for
  strange reasons, try running ``make clean`` in the project root
  and then rebuilding.
- Beware of multiple versions of the python library installed by different
  programs (e.g., pip versus installing locally from source)! In python,
  ``msprime.__file__`` will tell you the location of the package that is being
  used.
