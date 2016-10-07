********************
[0.4.0] - 2016-10-16
********************

Major release providing new functionality and laying groundwork for
upcoming functionality.

**Breaking changes**:

- The HDF5 file format has been changed to allow for non-binary trees
  and to improve performance. It is now both smaller and faster to
  load. However, msprime cannot directly load tree sequence files
  written by older versions. The ``msp upgrade`` utility has been
  developed to provide an upgrade path for existing users, so that
  files written by older versions of msprime can be converted to the
  newer format and read by version 0.4.x of msprime.

- The tuples returned by the ``mutations`` method contains an element.
  This will break code doing things like

      for pos, node in ts.mutations():
          print(pos, node)

  For better forward compatibility, code should use named attributes
  rather than positional access:

      for mutation in ts.mutations():
          print(mutation.position, mutation.node)

- Similarly, the undocumented ``variants`` method has some major changes:

  1. The returned tuple has two new values, ``node`` and ``index``
    in the middle of the tuple (but see the point above about using
    named attributes).

  2. The returned genotypes are by default numpy arrays. To revert
    to the old behaviour of returning Python bytes objects, use the
    ``as_bytes`` argument to the ``variants()`` method.

**New features**:

- Historical samples. Using the ``samples`` argument to ``simulate``
  users can specify the location and time of all samples explicitly.

- HDF5 file upgrade utility ``msp upgrade``

- Support for non-binary trees in the tree sequence, and relaxation
  of the requirements on input tree sequences using the read_txt()
  function.

- Integration with numpy, with zero-copy access to the low-level C API.

- Documented the variants() method that provides access to the sample
  genotypes as either numpy arrays or Python bytes objects.

- New LdCalculator class that allows very fast calculation of r^2 values.

- Initial support for threading.

- The values returned mutations() method now also contain an ``index``
  attribute. This makes many operations simpler.

- New TreeSequence.get_time() method that returns the time a sample
  was sampled at.

**Performance improvements**:

- File load times substantially reduced by pre-computing and storing
  traversal indexes.

- O(1) implementation of TreeSequence.get_num_trees()

- Improved control of enabled tree features in TreeSequence.trees()
  method using the ``leaf_lists`` and ``leaf_counts`` arguments.

**Bug fixes**:

- Fixed a precision problem with DemographyDebugger. #37

- Segfault on large haplotypes. #29

********************
[0.3.2] - 2016-07-21
********************

Feature release adding new import and export features to the API
and CLI.

- New ``TreeSequence.write_records`` and ``TreeSequence.write_mutations``
  methods to serialise a tree sequence in a human readable text format.

- New ``msprime.load_txt()`` method that parses the above formats, and
  allows msprime to read in data from external sources.

- New ``TreeSequence.write_vcf`` method to write mutation information
  in VCF format.

- Miscellaneous documentation fixes.


********************
[0.3.1] - 2016-06-24
********************

Feature release adding population related methods to the API.

- New ``TreeSequence.get_population(sample_id)`` method.

- New ``TreeSequence.get_samples(population_id)`` method.

- Added the optional ``samples`` argument to the
  ``TreeSequence.get_pairwise_diversity`` method.

- Fixed a potential low-level buffer overrun problem.


********************
[0.3.0] - 2016-05-31
********************

Bugfix release affecting all users of the Python API. Version 0.2.0 contained a
confusing and inconsistent mix of times and rates being expressed in both
coalescent units and generations. This release changes _all_ times and rates
used when describing demographic models to generations, and also changes
all population sizes to be absolute. In the interest of consistency, the
units of the trees output by msprime are also changed to generations. This
is a major breaking change, and will require updates to all scripts using the
API.

This release also include some performance improvements and additional
functionality.

Mspms users are not affected, other than benefiting from performance
improvements.

Breaking changes:

- Time values are now rescaled into generations when a TreeSequence is
  created, and so all times associated with tree nodes are measured in
  generations. The time values in any existing HDF5 file will now be
  interpreted as being in generations, so stored simulations must be
  rerun. To minimise the chance of this happening silently, we have
  incremented the file format major version number, so that attempts
  to read older versions will fail.

- Growth rate values for the PopulationConfiguration class are now
  per generation, and population sizes are absolute. These were in
  coalescent units and relative to Ne previously.

- GrowthRateChangeEvents and SizeChangeEvents have been replaced with
  a single class, PopulationParametersChange. This new class takes
  an initial_size as the absolute population size, and growth_rate
  per generation. Since the change in units was a breaking one,
  potentially leading to subtle and confusing bugs, we decided that
  the name refactoring would at least ensure that users would need
  to be aware that the change had been made. This API should now
  be stable, and will not be changed again without an excellent
  reason.

- MigrationRateChangeEvent has been renamed to MigrationRateChange
  and the migration rates are now per-generation.

- MassMigrationEvent has been renamed to MassMigration, and the
  values of source and destination swapped, fixing the bug in
  issue #14.

- The TreeSequence.records() method now returns an extra value,
  potentially breaking client code.

Improvements:

- Added tutorial for demographic events.

- Added DemographyDebugger class to help view the changes in populations
  over time.

- Added population tracking for coalescent events. We can now determine
  the population associated with every tree node. The relevant information
  has been added to the HDF5 file format.

- Improved performance for replication by reusing the same low-level
  simulator instance. This leads to significant improvements for large
  numbers of replicates of small simulations. Issue #8.

- Changed the TreeSequence.records() method to return named tuples.

- Added get_total_branch_length method. Issue #12.

- Fixed bug in reading Hapmap files. Issue #13.

********************
[0.2.0] - 2016-05-05
********************

Major update release, adding significant new functionality to the Python
API and several breaking changes. All code written for the 0.1.x API
will be affected, unfortunately.

Breaking changes:

- Sample IDs are now zero indexed. In previous versions of msprime, the
  samples were numbered from 1 to n inclusive, which is not Pythonic.
  This change has been made to make the API more usable, but will
  cause issues for existing code.

- There is now an Ne parameter to simulate(), and recombination,
  mutation and migration rates are now all per-generation. The
  keyword arguments have been changed to recombination_rate
  and mutation_rate, which should mean that silent errors will
  be avoided. All rates in existing code will need to be
  divided by 4 as a result of this. This change was made to make
  working with recombination maps and per generation recombination
  rates easier.

- Msprime now uses continuous values to represent coordinates, and
  the num_loci parameter has been replaced with a new length parameter
  to simulate(). Internally, a discrete recombination model is still
  used, but by default the potential number of discrete sites is
  very large and effectively continuous. True discrete recombination
  models can still be specified by using the recombination_map
  argument to simulate.

- The population_models argument to simulate() has been removed, and
  replaced with the population_configuration and demographic_events
  parameters. This was necessary to provide the full demographic
  model.

- The HDF5 file format has been updated to accommodate the continuous
  coordinates, along with other minor changes. As a consequence,
  simulation results will be somewhat larger. Stored simulations will
  need to be re-run and saved.

- Removed the random_seed key from the provenance JSON strings.

- Removed the simulate_tree() function, as it seemed to offer little
  extra value.


New features:

- Simulation of variable recombination rates via arbitrary recombination
  maps.

- Full support for population structure and demographic events.

- API support for replication via the num_replicates argument to simulate().

- Fully reworked random generation mechanisms, so that in the nominal
  case a single instance of gsl_rng is used throughout the entire
  simulation session.

- Addition of several miscellaneous methods to the TreeSequence API.

- Added NULL_NODE constant to make tree traversals more readable.

*********************
[0.1.10] - 2016-04-21
*********************

Bugfix release. Fixes serious issue affecting simulations with small
sample sizes.

https://github.com/jeromekelleher/msprime/issues/7

All users of mspms should update immediately and any analyses using
a small sample size (< 10) with mutations should be repeated.

Many thanks to Konrad Lohse for identifying the issue.

********************
[0.1.9] - 2016-04-01
********************

Bugfix release. Fixes serious issue affecting random seeds in mspms.

https://github.com/jeromekelleher/msprime/issues/6

All users of mspms should update immediately and any analyses using
the ``-seeds`` option in mspms should be repeated.

Many thanks to Derek Setter for identifying the issue.

********************
[0.1.8] - 2016-02-17
********************

Transitional release providing population structure support for the
ms-compatible command line interface. A considerable amount of low-level
plumbing was required to provide the required flexibility. This is currently
not visible from the high-level API, but will shortly be made available in the
forthcoming 0.2.x series.

The current implementation of migration should work well for small numbers of
populations (e.g. < 10), but will not scale well for large numbers of
populations.

+++++++
Changes
+++++++

- Added the -I, -m, -ma, -em, -eM, -ema, -eG, -eg, -eN, -en,
  -ej and -es options to mspms. These should provide full ms
  compatability, except for the -es option which is currently
  limited in scope.

- Added some extra keys to the low-level configuration JSON in
  the HDF5 file format to describe the population structure.
  This will be documented in a future release.

- Added a `get_pairwise_diversity` method to the TreeSequence
  class to efficiently calculate the population genetics
  statistic pi.
