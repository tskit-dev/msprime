(sec_cli)=

# Command line interface

Two command line applications are provided with `msprime`: {ref}`sec_msp` and
{ref}`sec_mspms`. The {command}`msp` program is a POSIX compliant command line
interface to the library. The {command}`mspms` program is a fully-{command}`ms`
compatible interface. This is useful for those who wish to get started quickly
with using
the library, and also as a means of plugging `msprime` into existing work
flows. However, there is a substantial overhead involved in translating data
from `msprime`'s native history file into legacy formats like `ms`, and so new code
should use the Python API where possible.

(sec_msp)=

## msp

The {command}`msp` program provides a convenient interface to the `msprime` API.
It is based on subcommands that either generate or consume a
tree sequence file. The `ancestry` sub-command simulates tree sequences from the
 coalescent with recombination. The `mutations` sub-command places
mutations onto an existing tree sequence. Several
{ref}`mutation models<sec_mutations_models>` are available.
The deprecated `simulate` sub-command simulates tree sequences from the
coalescent with recombination and mutations.

:::{important}
The {command}`msp ancestry` and {command}`msp mutations` commands write
their output to stdout by default so make sure you redirect this to a
file or use the ``--output`` option.
:::

### Examples

Simulate 10 diploid samples from a population of size 1000
with a genome of 100 base pairs, and write the output to ``ancestry.trees``:

```sh
$ msp ancestry 10 -N 1000 -L 100 > ancestry.trees
```

Simulate mutations from the {class}`Jukes-Cantor<.JC69>` mutation model
at rate 0.01 per-base per-generation, and write this to ``mutations.trees``.

```sh
$ msp mutations 0.01 ancestry.trees > mutations.trees
```

Do the same simulations, but pipe the output of the ancestry simulation
directly to the mutations simulation:
```sh
$ msp ancestry 10 -N 1000 -L 100 | msp mutations 0.01 > combined.trees
```

Show a summary of the properties of the trees file ``combined.trees``:
```
$ tskit info combined.trees
sequence_length:  100.0
trees:            1
samples:          20
individuals:      10
nodes:            39
edges:            38
sites:            100
mutations:        16997
migrations:       0
populations:      1
provenances:      2
```

See the {ref}`tskit documentation<tskit:sec_cli>` for more details on the
{program}`tskit` command line interface.

(sec_msp_ancestry)=

### msp ancestry

{command}`msp ancestry` generates coalescent simulations with recombination
from a constant population size and stores the result as a tree sequence in
an output file. This sub-command is an interface to the
{func}`msprime.sim_ancestry` API function.

```{eval-rst}
.. argparse::
    :module: msprime.cli
    :func: get_msp_parser
    :prog: msp
    :path: ancestry
    :nodefault:
```


(sec_msp_mutate)=

### msp mutations

{command}`msp mutate` can be used to add mutations to a tree sequence and store
a copy of the resulting tree sequence in a second file. This
sub-command is an interface to the
{func}`msprime.sim_mutations` API function.

```{eval-rst}
.. argparse::
    :module: msprime.cli
    :func: get_msp_parser
    :prog: msp
    :path: mutations
    :nodefault:
```

### msp simulate

:::{warning}
The {command}`msp simulate` command is deprecated.
:::

{command}`msp simulate` generates coalescent simulations with recombination
from a constant population size and stores the result as a tree sequence in
an output file.
{command}`msp simulate` is deprecated, but will be supported indefinitely.
{ref}`sec_msp_mutate` provides further mutation models.
{command}`msp ancestry` will provide further demographic scenarios from
which to simulate tree sequences. {command}`msp simulate` is an
interface to the deprecated {func}`msprime.simulate` API function.

```{eval-rst}
.. argparse::
    :module: msprime.cli
    :func: get_msp_parser
    :prog: msp
    :path: simulate
    :nodefault:
```

(sec_mspms)=

## mspms

The {command}`mspms` program is an {command}`ms`-compatible
command line interface to the `msprime` library. This interface should
be useful for legacy applications, where it can be used as a drop-in
replacement for {command}`ms`. This interface is not recommended for new applications,
particularly if the simulated trees are required as part of the output
as Newick is very inefficient. The Python API is the recommended interface,
providing direct access to the structures used within `msprime`.

### Supported Features

{command}`mspms` supports a subset of {command}`ms`'s functionality. Please
[open an issue](<https://github.com/tskit-dev/msprime/issues>) on
GitHub if there is a feature of {command}`ms` that you would like to see
added. We  currently support:

- Basic functionality (sample size, replicates, tree and haplotype output);
- Recombination (via the `-r` option);
- Gene-conversion (via the `-c` option);
- Spatial structure with arbitrary migration matrices;
- Support for {command}`ms` demographic events. (The implementation of the
  `-es` option is limited, and has restrictions on how it may be
  combined with other options.)

### Argument details

This section provides the detailed listing of the arguments to
{command}`mspms` (also available via `mspms --help`). See
the [documentation for ms](<http://thirteen-01.stat.iastate.edu/snoweye/phyclust/document/msdoc.pdf>)
for details on how these values should be interpreted.

```{warning}
Due to quirks in Python's argparse module, negative growth rates
written in exponential form (e.g. `-eG 1.0 -1e-5`) are not recognised as an
option argument. To work around this, specify the argument using quotes
and a leading space, e.g. `-eG 1.0 ' -1e-5'`, or avoid scientific notation,
e.g. `-eG 1.0 -0.00001`.
```

```{eval-rst}
.. argparse::
    :module: msprime.cli
    :func: get_mspms_parser
    :prog: mspms
    :nodefault:

```
