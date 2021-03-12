(sec_cli)=

# Command line interface

Two command line applications are provided with `msprime`: {ref}`sec_msp` and
{ref}`sec_mspms`. The {command}`msp` program is a POSIX compliant command line 
interface to the library. The {command}`mspms` program is a fully-{command}`ms`
compatible interface. This is useful for those who wish to get started quickly
with using
the library, and also as a means of plugging `msprime` into existing work
flows. However, there is a substantial overhead involved in translating data
from `msprime`'s native history file into legacy formats, and so new code
should use the Python API where possible.

(sec_msp)=

## msp

The {command}`msp` program provides a convenient interface to the `msprime` API.
It is based on subcommands that either generate or consume a
tree sequence file. The `ancestry` sub-command simulates tree sequences from the
 coalescent with recombination. The `mutate` sub-command places 
mutations onto an existing tree sequence. Several mutation models are available.
The deprecated `simulate` sub-command simulates tree sequences from the
 coalescent with recombination and mutations.

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

### msp mutate

{command}`msp mutate` can be used to add mutations to a tree sequence and store
a copy of the resulting tree sequence in a second file. This
sub-command is an interface to the
{func}`msprime.sim_mutations` API function. 

```{eval-rst}
.. argparse::
    :module: msprime.cli
    :func: get_msp_parser
    :prog: msp
    :path: mutate
    :nodefault:
```



### msp simulate

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

```{eval-rst}
.. argparse::
    :module: msprime.cli
    :func: get_mspms_parser
    :prog: mspms
    :nodefault:

```
