(sec_development)=

# Development

If you would like to add some features to `msprime`, please read the
following. If you think there is anything missing,
please open an [issue](<http://github.com/tskit-dev/msprime/issues>) or
[pull request](<http://github.com/tskit-dev/msprime/pulls>) on GitHub!

(sec_development_quickstart)=

## Quickstart

- Fork the msprime repo on [GitHub](<http://github.com/tskit-dev/msprime>) and
  clone your fork, making sure that the **submodules are correctly initialised**:

  ```
  $ git clone git@github.com:YOUR_GITHUB_USERNAME/msprime.git --recurse-submodules
  ```

  For an already checked out repo, initialise submodules with:

  ```
  $ git submodule update --init --recursive
  ```
- Install dependencies and build the low-level module: `uv sync`
- Install the prek pre-commit hook: `uv run prek install`
- Run the tests to ensure everything is working: `uv run pytest`
- Make your changes in a local branch and open a pull request when ready.

See the [tskit developer documentation](https://tskit.dev/tskit/docs/stable/development.html)
for details on the recommended git workflow, running lint checks, and building
the documentation.

(sec_development_requirements)=
## Requirements

### System requirements

To develop with msprime you will need to have [GSL](https://www.gnu.org/software/gsl/)
installed and a working compiler. Please see the {ref}`sec_pip_install_source`
section for how to install GSL on some common platforms.

:::{important}
You still need to install GSL and have a working compiler if you are working
on the documentation because it requires a locally built version of the
{ref}`C module<sec_development_c_module>`.
:::

### Python requirements

The packages needed for development are specified as dependency groups in
``python/pyproject.toml``. Install them along with the low-level C extension using:

```
$ uv sync
```

## Overview

There are three main parts of `msprime`, in increasing order of complexity:

1. High-level Python. The Python-API and command line interface tools are all defined
   in the `msprime` directory.
2. C library. The underlying high-performance C code is written as a standalone library.
   All of the code for this library is in the `lib` directory.
3. Low-level Python-C interface. The interface between the Python and C code is the
   `msprime/_msprimemodule.c` file, which defines the `msprime._msprime` module.

Each of these aspects has its own coding conventions and development tools, which are
documented in the following sections.


## Continuous integration tests

CI uses shared GitHub Actions workflows from the tskit-dev ecosystem. See the
[repo administration guide](https://github.com/tskit-dev/.github/blob/main/repo_administration.md)
for details of the available workflows and how they are configured.


(sec_development_documentation)=
## Documentation

The msprime manual exists to provide users with a comprehensive and authoritative
source of information on msprime's interfaces. From a high-level, the documentation
is split into two main sections:

1. The API Reference documentation provides a concise and precise description of
   a particular function or class. See the {ref}`sec_development_documentation_api`
   section for details.
2. Thematically structured sections which discuss the functionality and
   explain features via minimal examples.

Further documentation where features are combined to perform specific
tasks is provided in the [tskit tutorials](https://tskit.dev/tutorials) site.

Documentation can be written on GitHub or locally. If you are new to contributing
to msprime and you will be making minor edits to markdown text, you may find
it easier to make edits on GitHub. To do this, hover your mouse over the GitHub
icon in the top right corner of the documentation, and click "suggest edit".
You can then edit and preview markdown in GitHub's user interface. Clicking
"propose change" at the bottom of the page will commit to a new branch on your fork.
If you do not already have a fork of the msprime repository, GitHub will prompt you
to "fork this repository" - go ahead and create your fork.
You can then create a pull request for your proposed change. On the other hand,
if you are already familiar with contributing to msprime, or have more than simple
markdown edits to add, you can edit and
{ref}`build <sec_development_documentation_building>` the documentation locally.
To do that, first follow the {ref}`Quickstart <sec_development_quickstart>`.
Once you have created and checked out a "topic branch", you are ready to start
editing the documentation.

:::{note}
The documentation requires a working build of the low-level C module. This is
built automatically by `uv sync`, but if you see inscrutable import errors a
mismatch between the installed and local versions of the module is a common cause.
Rebuild by running `make` in the project root.
:::

(sec_development_documentation_building)=
### Building

To build the documentation locally, go to the `docs` directory and run `make`
(ensure that the {ref}`sec_development_requirements` have been installed).
This will build the HTML documentation in  `docs/_build/html/`. You can now
view the local build of the HTML in your local browser (if you do not know how
to do this, try double clicking the HTML file).

:::{note}
If you are having some problems with getting the documentation to build
try running ``make clean`` which will delete all of the HTML and cached
Jupyter notebook content.
:::

(sec_development_documentation_api)=
### API Reference

API docstrings are written in rST. The msprime docstrings predate the switch to
MyST Markdown and converting them would be a significant effort, so rST remains
the format for docstrings for the foreseeable future.

Docstrings should be **concise** and **precise**. Examples should not be
embedded directly in docstrings; instead, each significant parameter should
link to the relevant section in the narrative documentation.

See the [tskit developer documentation](https://tskit.dev/tskit/docs/stable/development.html)
for guidance on markup languages, code examples, and cross referencing.

## High-level Python

The `msprime` package is installed in editable mode by `uv sync`. The low-level
C extension is also built automatically at that point.

### Conventions

Python code is formatted and linted using [ruff](https://docs.astral.sh/ruff/),
run automatically via [prek](https://prek.j178.dev) on each commit.

### Tests

Tests are in the `python/tests` directory and run with
[pytest](<https://docs.pytest.org/en/stable/>). Core simulation and basic tests
are in `python/tests/test_highlevel.py`; more focused tests are in smaller files
(e.g., `test_demography.py`). Run the suite with:

```
$ uv run pytest
```

All new code must have high test coverage, tracked by
[CodeCov](<https://codecov.io/gh/tskit-dev/msprime/>).

### Interfacing with low-level module

Much of the high-level Python code only exists to provide a simpler interface to
the low-level `_msprime` module. As such, many objects (such as `RecombinationMap`)
are really just a shallow layer on top of the corresponding low-level object.
The convention here is to keep a reference to the low-level object via
a private instance variable such as `self._ll_recombination_map`.

## C Library

The low-level code for `msprime` is written in C, and is structured as a
standalone library. This code is all contained in the `lib` directory.
Although the code is structured as a library, it is not intended to be used
outside of the `msprime` project! The interfaces at the C level change
considerably over time, and are deliberately undocumented.

### Toolchain

To compile and develop the C code, a few extra development libraries are needed.
[Libconfig](http://www.hyperrealm.com/libconfig/) is used for the development CLI
and [CUnit](http://cunit.sourceforge.net) for unit tests. We use the
[meson](https://mesonbuild.com) build system in conjunction with [ninja-build](http://ninja-build.org) to compile the unit tests and
development CLI. On Debian/Ubuntu, these can be installed using

```{code-block} bash

$ sudo apt-get install libcunit1-dev libconfig-dev ninja-build

```

Meson can be installed with:

```{code-block} bash

$ uv tool install meson

```

On macOS, `brew install cunit` can be used in place of `apt-get`.

### Compiling

Meson keeps all compiled binaries in a build directory (this has many advantages
such as allowing multiple builds with different options to coexist). It depends on
a `meson.build` file which is in the `lib` directory. To set up the initial build
directory, run

```{code-block} bash

$ cd lib
$ meson build

```

The easiest way to compile the {ref}`sec_development_c_unit_tests`
is to run `ninja -C build`. (Alternatively,
you can `cd` into the `build` directory and run `ninja`). All the
compiled binaries are then in the `build` directory, so to run, for example, the
`test_ancestry` unit tests, use `./build/test_ancestry`. A handy shortcut
to compile the code and run all the unit tests is:

```{code-block} bash

$ ninja -C build test

```

The [mesonic](http://www.vim.org/scripts/script.php?script_id=5378) plugin for vim
simplifies this process and allows code to be compiled seamlessly within the
editor.

### Development CLI

When developing the C code, it is usually best to use the development CLI to invoke
the code. This is much simpler than going through the Python interface, and allows
tools such as [valgrind](<http://valgrind.org>) to be used directly. For example,
when developing new simulation functionality, you should get the basic work done
using the CLI and only move over to the Python API once you are reasonably sure
that the code works properly.

The development CLI is written using [libconfig](http://www.hyperrealm.com/libconfig/) to parse the simulation parameters
file, and [argtable3](https://github.com/argtable/argtable3) to parse the
command line arguments. The `argtable3` code is included in the source (but
not used in the distributed binaries, since this is strictly a development
tool). The source code is in `dev-tools/dev-cli.c`.

After building, the CLI is run as follows:

```{code-block} bash

$ ./build/dev-cli <command> <arguments>

```

Running the `dev-cli` program without arguments will print out a summary of the
options.

<!---
warning

The development CLI is a tool used to develop the msprime API, and not a
polished artefact intended for users. There is quite a lot of code left
over from earlier debugging which might not make immediate sense. Some
commands may not work as expected, or indeed at all. Please feel free to
tidy it up if you would like to improve it!
-->

The most important command for simulator development is `simulate`,
which takes a configuration file as a parameter and writes the resulting
simulation to an output file in the native `.trees` format. For example,

```{code-block} bash

$ ./build/dev-cli simulate dev-tools/example.cfg -o out.trees

```

The development configuration file describes the simulation that we want to
run, and uses the
[libconfig syntax](<http://www.hyperrealm.com/libconfig/libconfig_manual.html#Configuration-Files>).
An example is given in the file `dev-tools/example.cfg` which should have sufficient documentation
to be self-explanatory.

<!---
warning

It is important to note that all values in the low-level C code are in
scaled coalescent units. The high-level Python API defines values in units
of generations, but for the C code all time is measured in coalescent units.
-->

(sec_development_c_unit_tests)=

### Unit Tests

The C-library has an extensive suite of unit tests written using
[CUnit](<http://cunit.sourceforge.net>). These tests aim to establish that the
low-level APIs work correctly over a variety of inputs, and particularly, that
the tests don't result in leaked memory or illegal memory accesses. The tests should be
periodically run under valgrind to make sure of this.

Tests are defined in the `tests` directory, roughly split into suites
defined in different files. For example, the tests associated with Fenwick
trees are defined in the `tests/tests_fenwick.c` file. To run all the
tests in this suite, use run using `./build/test_fenwick`.
To run a specific test in a particular suite, provide the name of the
test name as a command line argument, e.g.:

```{code-block} bash

$ ./build/test_fenwick test_fenwick_expand

```

While 100% test coverage is not feasible for C code, we aim to cover all code
that can be reached. (Some classes of error such as malloc failures
and IO errors are difficult to simulate in C.) Code coverage statistics are
automatically tracked using [CodeCov](<https://codecov.io/gh/tskit-dev/msprime/>).

### Code Style

C code is formatted using
[clang-format](<https://clang.llvm.org/docs/ClangFormat.html>)
with a custom configuration. This is checked automatically by prek. To format
all files manually run:

```{code-block} bash

$ uv run prek --all-files

```

### Coding conventions

The code is written using the [C99](<https://en.wikipedia.org/wiki/C99>) standard. All
variable declarations should be done at the start of a function, and functions
kept short and simple where at all possible.

No global or module level variables are used for production code.

The code is organised following object-oriented principles. Each 'class' is defined using
a struct, which encapsulates all the data it requires. Every 'method' on this class
is then a function that takes this struct as its first parameter. Each class has
an `alloc` method, which is responsible for allocating memory and a `free` method
which frees all memory used by the object. For example, the
[Fenwick tree](<https://en.wikipedia.org/wiki/Fenwick_tree>) class is defined as
follows:

```{code-block} C

typedef struct {
    size_t size;
    size_t log_size;
    double *tree;
    double *values;
} fenwick_t;

int fenwick_alloc(fenwick_t *self, size_t initial_size);
int fenwick_free(fenwick_t *self);
double fenwick_get_total(fenwick_t *self);

```

This defines the `fenwick_t` struct, and alloc and free methods and a method
to return the total of the tree. Note that we follow the Python convention
and use `self` to refer to the current instance.

Most objects also provide a `print_state` method, which is useful for
debugging.

```{eval-rst}
.. todo:: Change to intersphinx mapping for this link.
```

Please see the documentation for the
[tskit C API](https://tskit.readthedocs.io/en/stable/c-api.html#sec-c-api-overview-structure)
for more details on the how APIs are structured.

### Error handling

A critical element of producing reliable C programs is consistent error handling
and checking of return values. All return values **must** be checked! In msprime,
all functions (except the most trivial accessors) return an integer to indicate
success or failure. Any negative value is an error, and must be handled accordingly.
The following pattern is canonical:

```{code-block} C

    ret = msp_do_something(self, argument);
    if (ret != 0) {
        goto out;
    }
    // rest of function
out:
    return ret;

```

Here we test the return value of `msp_do_something` and if it is non-zero,
abort the function and return this same value from the current function. This
is a bit like throwing an exception in higher-level languages, but discipline
is required to ensure that the error codes are propagated back to the original
caller correctly.

Particular care must be taken in functions that allocate memory, because
we must ensure that this memory is freed in all possible success and
failure scenarios. The following pattern is used throughout for this purpose:

```{code-block} C

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

```

It is vital here that `x` is initialised to `NULL` so that we are guaranteed
correct behaviour in all cases. For this reason, the convention is to declare all
pointer variables on a single line and to initialise them to `NULL` as part
of the declaration.

Error codes are defined in `err.h`, and these can be translated into a
message using `msp_strerror(err)`.

### Running valgrind

Valgrind is an essential development tool, and is used extensively. (Being able
to run valgrind was one of the motivating factors in the C-library architecture.
It is difficult to run valgrind on a Python extension module, and so the simplest
way to ensure that the low-level code is memory-tight is to separate it out
into an independent library.)

Any new C unit tests that are written should be verified using valgrind to
ensure that no memory is leaked. The entire test suite should be run
through valgrind periodically also to detect any leaks or illegal
memory accesses that have been overlooked.


(sec_development_c_module)=
## Python C Interface

The Python C interface is written using the
[Python C API](<https://docs.python.org/3.6/c-api/>) and the code is in the
`msprime/_msprimemodule.c` file. When compiled, this produces the
`msprime._msprime` module,
which is imported by the high-level module. The low-level Python module is
not intended to be used directly and may change arbitrarily over time.

The conventions used within the low-level module here closely follow
those in `tskit`; please see the
[tskit developer documentation](https://tskit.dev/tskit/docs/stable/development.html#python-c-interface)
for more information.

## Statistical tests

To ensure that `msprime` is simulating the correct process we run many statistical
tests. Since these tests are quite expensive (taking some hours to run) and
difficult to automatically validate, they are not run as part of CI but instead
as a pre-release sanity check. They are also very useful to run when developing
new simulation functionality, as subtle statistical bugs can easily slip in
unnoticed.

The statistical tests are all run via the `verification.py` script in the project root,
using extra dependencies declared in the `verification` dependency group. Run using:

```{code-block} bash

$ uv run --group verification python verification.py

```

The statistical tests depend on compiled programs in the `data` directory.
This includes a customised version of `ms` and a locally compiled version of
[scrm](<https://scrm.github.io/>). These programs must be compiled before
running the statistical tests, and can be built by running `make` in the
`data` directory. If this is successful, there should be several binaries
like `ms` and `ms_summary_stats` present in the `data`
directory.

Please read the comments at the top of the `verification.py` script for details
on how to write and run these tests.

## Benchmarking

Benchmarks to measure performance are in the `benchmarks` folder and are run using
[airspeed velocity](https://asv.readthedocs.io/en/stable/index.html).
A system that runs the benchmarks on each push to the main branch.
These benchmarks can also be run locally to compare your branch with the main branch.
Your changes must be in a commit to be measured. To run the benchmarks:

```
asv run HEAD...main~1
```

This will run the benchmarks for the latest main branch commit and all commits on
your current branch (the syntax for choosing commits is the same as `git log`).
The following commands then make a browsable report (link given in output of
the command):

```
asv publish
asv preview
```

Note the following tips:

- Specifying the range of commits to run uses the same syntax as git log.
  For example, to run for a single commit, use `asv run 88fbbc33^!`
- Be careful when running `asv dev` or using `python=same` as
  this can use the *installed* version of msprime rather than the local
  development version. This can lead to confusing results! When tuning
  benchmarks it's better to commit often and use (e.g.)
  `asv run HEAD^! --show-stderr -b Hudson.time_large_sample_size`.
- You may want to benchmark a specific list of commits exclusively. To do so, put the commits' hashes in a file and use the command: `asv run --show-stderr --skip-existing HASHFILE:hashestobenchmark.txt`

- There is a script: `benchmarks/check_asv.sh` that can be used to benchmark recent commits.

## Troubleshooting

- If `make` is giving you strange errors, or if tests are failing for
  strange reasons, try running `make clean` in the project root
  and then rebuilding.
- Beware of multiple versions of the python library being visible. In python,
  `msprime.__file__` will tell you the location of the package that is being
  used.
