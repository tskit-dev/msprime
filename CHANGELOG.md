# Changelog

## [1.2.1] - 2022-XX-XX

**New features**

- Add a `MicrosatMutationModel` mutation model class, that
  represents a generalized interface for constructing mutational
  models appropriate to STRs. In addition 3 specific microsat models
  are added `SMM`, `TPM`, and `EL2`.  
  ({issue}`2013`, {user}`andrewkern`).

- Raise an error if `log_arg_likelihood` is called on a recombinant tree
  sequence that has no recombination nodes (({issue}`2123`, {pr}`2124`,
  {user}`hyanwong`)

## [1.2.0] - 2022-05-18

**New features**

- Add the `FixedPedigree` ancestry model and various infrastructure for
  importing pedigree information into msprime.

**Bug fixes**:

- Fix rare assertion trip in the single sweep model caused by numerical jitter.
  ({issue}`1966`, {pr}`2038`, {user}`jeromekelleher`, {user}`molpopgen`)

- Fix edge case in `Demography.from_old_style()`
  ({issue}`2047`, {pr}`2048`, {user}`grahamgower`)

**Maintenance**:

- Documentation improvements ({pr}`2054`, {pr}`2033`, {pr}`2011`
  {user}`petrelharp`, {user}`gregorgorjanc`)

## [1.1.1] - 2022-02-10

**Bug fixes**:

- Fix (very) rare assertion trip caused by underlying GSL bug.
  ({issue}`1997`, {pr}`2000`, {user}`chriscrsmith`, {user}`molpopgen`,
  {user}`andrewkern`)

**Maintenance**:

- Various documentation improvements.

## [1.1.0] - 2021-12-14

**New features**

- Add support for tree sequence ``time_units`` field. The ``time_units`` will
  be set to "generations" for the output of ``sim_ancestry`` (and ``simulate``),
  unless the ``initial_state`` argument is used. In this case, the
  ``time_units`` value will be inherited from the input.
  ({pr}`1953`, {issue}`1951`, {issue}`1877`, {issue}`1948`, {user}`jeromekelleher`).

**Bug fixes**:

- Raise an error if `Demography.from_demes()` is passed a model
  with non-zero `selfing_rate` or `cloning_rate` values (which msprime
  does not support).
  ({pr}`1938`, {issue}`1937`, {user}`grahamgower`).

- Do not assume Population metadata schemas contain the ``properties``
  and ``additionalProperties`` attributes ({issue}`1947`, {pr}`1954`,
  {user}`jeromekelleher`).

- Read the population name from PopulationConfiguration ``metadata`` in
  ``Demography.from_old_style`` ({issue}`1950`, {pr}`1954`,
  {user}`jeromekelleher`)

**Maintenance**:

- Update tskit to Python 0.4.0 and C 0.99.15.

## [1.0.4] - 2021-12-01

**New features**:

- Support for Demes 0.2.0, which introduces a change to how pulse
  sources and proportions are specified.
  ({pr}`1936`, {issue}`1930`, {user}`apragsdale`)

## [1.0.3] - 2021-11-12

This is a bugfix release recommended for all users.

**New features**:

- Support for running full ARG simulations with gene conversion
  ({pr}`1801`, {issue}`1773`, {user}`JereKoskela`).

- Improved performance when running many small simulations
  ({pr}`1909`, {user}`jeromekelleher`.)

- Update to tskit C API 0.99.14 ({pr}`1829`).

**Bug fixes**:

- Fix bug in full ARG simulation with missing regions of the genome,
  where ARG nodes were not correctly returned. ({issue}`1893`,
  {user}`jeromekelleher`, {user}`hyl317`)

- Fix memory leak when running ``sim_ancestry`` in a loop
  ({pr}`1904`, {issue}`1899`, {user}`jeromekelleher`, {user}`grahamgower`).

- Fix printing small values in rate maps ({pr}`1906`, {issue}`1905`,
  {user}`petrelharp`).

## [1.0.2] - 2021-06-29

Improved Demes support and minor bugfixes.

**New features**:

- Support for Demes input and logging in the msp simulate CLI
  ({pr}`1716`, {user}`jeromekelleher`).
- Add ``Demography.to_demes`` method for creating a Demes demographic model
  from an msprime demography ({pr}`1724`, {user}`grahamgower`).
- Improved mapping of Demes models to Demography objects
  ({pr}`1758`, {pr}`1757`, {pr}`1756` {user}`apragsdale`).
- Improved numerical algorithms in DemographyDebugger ({pr}`1788`,
  {user}`grahamgower`, {user}`petrelharp`).

**Bugfixes**:

- Raise an error if running full ARG simulations with gene conversion
  ({issue}`1774`).

## [1.0.1] - 2021-05-10

Minor feature release with experimental Demes support.

- Change the semantics of Admixture events slightly so that ancestral
  populations that are inactive, are marked as active ({pr}`1662`,
  {issue}`1657`, {user}`jeromekelleher`, {user}`apragsdale`)

- Initial support for Demes via the ``Demography.from_demes`` method.
  ({pr}`1662`, {issue}`1675`, {user}`jeromekelleher`, {user}`apragsdale`)

## [1.0.0] - 2021-04-14

Msprime 1.0 is a major update, recommended for all users. It introduces
new APIs and many new features, which still retaining compatibility
with nearly all existing code.

**New features**:

- Add new top-level functions `sim_ancestry` and `sim_mutations`, intended
  as long term replacements for `simulate` and `mutate`.
- Add new `Demography` class as a replacement for the 0.x
  ``population_configurations``, ``migration_matrix`` and ``demographic_events``
  parameters to simulate. Many new features to improve ease of use and to
  help avoid errors.
- Change the underlying simulation from working in discrete genetic coordinates
  and translating to physical coordinates, to simulating directly in physical
  coordinates. This fixes some long standing bugs and makes it far simpler
  to add features such as gene conversion. ({user}`daniel-goldstein`)
- Support discrete or continuous genomes in simulations directly.
- Gene conversion ({user}`fbaumdicker`)
- Selective sweeps and low-level infrastructure for the structured coalescent
  ({user}`andrewkern`, {user}`gbinux`)
- Multiple merger coalescent models ({user}`jerekoskela`, {user}`shajoezhu`,
  {user}`TPPSellinger`)
- Different ploidy levels for coalescent models ({user}`jerekoskela`)
- Functions to parse species trees and set up simulation models according
  to the species tree. ({user}`mmatschiner` {issue}`893` {issue}`929` {issue}`931`)
- Add an modern array oriented RateMap class as to replace the RecombinationMap
  class ({user}`jeromekelleher`, {user}`grahamgower`, {user}`hyanwong`).
- Functions to compute the likelihood of a particular ARG realisation
  under the coalescent with recombination ({user}`jerekoskela`).
- Mutations are assigned times ({user}`petrelharp`)
- Finite sites mutations ({user}`petrelharp`)
- Several instances of matrix mutation models: JC69, GTR, HKY, F84, BLOSUM62,
  PAM. ({user}`GertjanBisschop`, {user}`petrelharp`).
- An implementation of the SLiM mutation model ({user}`petrelharp`)
- An infinite alleles mutation model ({user}`jeromekelleher`)
- Methods to track the possible location of lineages, and compute the
  coalescence rates over time ({user}`apragsdale`, {user}`petrelharp`
  {user}`grahamgower`).
- Two new CLI commands, `msp ancestry` and `msp mutations` which correspond
  to the `sim_ancestry` and `sim_mutations` functions. These can be pipelined
  using stdin/stdout. ({user}`winni2k`, {user}`jeromekelleher`)
- Complete provenance recording of all arguments to simulate and mutate.
  Adds argument record_provenance to simulate, which allows recording of
  provenances to be disabled, for example when they are large.
  ({user}`benjeffery` {pr}`914`).
- Add replicate_index to simulate, allowing output of a single tree sequence
  from a set of replicates. ({user}`benjeffery` {pr}`914`).
- Details of the simulation are written to the DEBUG log periodically.
  This can help debug long-running simulations. ({user}`jeromekelleher`,
  {pr}`1080`).
- Binary wheels for PyPI ({user}`benjeffery`)

**Breaking changes**:

- Require Python 3.7+.
- The class form for specifying models (e.g., `msprime.StandardCoalescent()`)
  no longer take a `reference_size` argument. ({user}`jeromekelleher`,
  {pr}`1028`).
- The `simulate` function only takes one positional argument, and all other
  arguments are keyword-only.
- The `msp` CLI has been stripped of all existing sub-commands except for
  `simulate`. The old sub-commands are provided by the `tskit`
  CLI or the `TreeSequence` API in `tskit`.
- The `end_time` option now allows events up to and including the specified
  max time. Previously, events occurred strictly before the max time.
  ({user}`apragsdale`, {issue}`1304`, {pr}`1521`)
- The class msprime.Population now refers to the new demography API. This
  may break old code that depends on the tskit table entities being
  exported from msprime (other classes are now formally deprecated and
  will raise a warning; see the deprecations section).
- The semantics of ``start_time`` and exponential growth models has changed;
  previously, exponential population size changes would be calculated starting
  at the simulation ``start_time`` at the beginning of the simulation. Now
  population growth rates are always computed from time 0 (unless other
  events occur). The old semantics, if needed, can be recovered by setting
  the ``growth_rate`` to 0 in the initial population declaration and adding
  a population parameters change event at the ``start_time`` to change the
  ``growth_rate`` to the desired value.

**Performance improvements**

- Much better simulation performance for models with large numbers
  of populations ({user}`jeromekelleher`, {pr}`1069`).
- Significant performance improvements for simulating from recombination
  maps ({user}`ivan-krukov`, {user}`castedo`)

**Deprecations**:

- Deprecate module attributes that were moved to tskit.
  ({user}`benjeffery`, {issue}`991`, {pr}`1158`)

## [0.7.4] - 2019-12-05

**Bug fixes**:

- Fix error in mspms output of tree spans. In previous versions, the length of
  genome spanned by trees in the newick output was incorrect in certain situations
  (specifically, when "invisible" recombinations are present so that two or more
  identical trees are printed out). Thanks to {user}`fbaumdicker` for spotting
  the problem. ({user}`jeromekelleher`, {pr}`837`, {issue}`836`)
- Fix assertion tripped when we have very low recombination rates in the DTWF
  model. Thanks to {user}`terhorst` for the bug report.
  ({user}`jeromekelleher`, {pr}`833`, {issue}`831`).
- Fix bug in memory allocation when simulating mutations on a tree sequence
  that already contains many mutations. Thanks to {user}`santaci` for the
  bug report. ({user}`jeromekelleher`, {user}`petrelharp`, {pr}`838`,
  {issue}`806`)

**New features**:

- Add the new Census event, which allows us to place nodes on all extant
  branches at a given time ({user}`gtsambos` {pr}`799`).
- Improved error reporting for input parameters, in particular
  demographic events ({pr}`829`).

**Documentation**:

- Improved container documentation ({user}`agladstein`, {pr}`822`, {issue}`809`).
- Improved developer docs for macs ({user}`gtsambos`, {user}`molpopgen`, {pr}`805`).
- Clarify meaning of migration matrix ({user}`petrelharp`, {pr}`830`).

## [0.7.3] - 2019-08-03

**Bug fixes**:

- Support for SMC models coupled with the record_full_arg feature was
  erroneously removed in a previous version ({issue}`795`). The feature
  has been reinstated ({pr}`796`).

## [0.7.2] - 2019-07-30

**Breaking changes**

- The random trajectory has been changed slightly to improve handling
  of ancient sampling events ({pr}`782`). Thus, simulations for a given
  random seed will not be identical to previous versions, if ancient
  samples are used.

**New features**

- Automated Docker builds ({user}`agladstein`; {pr}`661`)
- Add mean coalescence time to DemographyDebugger ({user}`petrelharp`; {pr}`779`).
- Improve MassMigration descriptions in DemographyDebugger
  ({user}`marianne-aspbury`; {pr}`791`).

**Bug fixes**:

- In very, very, very rare cases it was possible to generate a
  zero waiting time until the next coalescent event, leading to
  zero branch lengths in the output tree sequence and an error
  being raised ({user}`molpopgen`, {user}`DL42`, {user}`jeromekelleher`;
  {issue}`783`, {pr}`785`).

## [0.7.1] - 2019-06-08

**New features**

- Discrete Time Wright-Fisher simulation model ({user}`DomNelson`).
- SMC/SMC' simulation models ({user}`jeromekelleher`).
- Mixed simulation models ({user}`jeromekelleher`).
- Specify `end_time` to allow early-finish for simulations ({user}`jeromekelleher`).
- Calculation of historical coalescence rates in the DemographyDebugger
  ({user}`jgallowa07`, {user}`petrelharp`).
- Additional information on population sizes in DemographyDebugger
  ({user}`andrewkern`).
- Remove support for Python 2 ({user}`hugovk`).
- Allow specifying metadata for populations ({user}`jeromekelleher`).

**Bug fixes**:

- Various minor bug and doc fixes from {user}`hyanwong`, {user}`petrelharp`,
  {user}`brianzhang01`, {user}`mufernando` and {user}`andrewkern`.

## [0.7.1b1] - 2019-05-31

Early release making DTWF code available to beta testers.

## [0.7.0] - 2019-02-19

Separation of tskit from msprime. Msprime is now solely dedicated to simulating
the coalescent, and all infrastructure for working with succinct tree sequences
is now provided by tskit. To ensure compatibility, msprime now imports code
from tskit under the old names, which should ensure that all code continues
to work without changes.

**New features**

- Ability to record the full ARG ({user}`jerekoskela`; {issue}`665`)

**Bug fixes**:

- Fix deprecation warning ({issue}`695`).

## [0.7.0a1] - 2019-01-14

Alpha release for testing the tskit/msprime split.

## [0.6.2] - 2018-12-04

Minor bugfix release.

**New features**:

- Add provenance recording option to simplify (#601)
- Minor performance improvement (#598)

**Bug fixes**:

- Fix performance regression in replication (#608)

## [0.6.1] - 2018-08-25

Significant features for integration with forwards-time simulators plus
improvements and bugfixes.

**Breaking changes**:

- Change in the semantics of how populations are treated by simplify. By
  default, populations that are not referenced will now be removed from the
  data model. This can be avoided by setting `filter_populations=False`.
- Simplify now raises an error if called on a set of tables that contain
  one or more migrations.

**New features**:

- The simulate() function now supports a `from_ts` argument allowing
  msprime to complete the ancestry in tree sequences generated by
  forward simulations (#503, #541, #572, #581).
- Add start_time and end_time parameters to the `mutate` function (#508).
- Add `reduce_to_site_topology` argument to simplify. This allows us to
  find the minimal tree sequence that would be visible from a given set
  of sites, and is also a useful compression method if we are only interested
  in the observed sequences. (#545, #307).
- Simplify generalised to support individuals, and the `filter_populations`,
  `filter_individuals` and `filter_sites` parameters added to allow
  filtering of unreferenced objects from the data model. (#567).
- Default random seeds are now generated from a sequence initialised by
  a system source of randomness (#534). Random seeds should also be safely generated
  across multiple processes.
- Full text I/0 support for Individuals and Populations (#498, #555)
- Substantially improved performance in `msprime.load` for large tables
  and significant refactoring of C code (#559, #567, #569).
- Improved performance of generating genotypes (#580).
- Formal schema for tree sequence provenance (#566, #583).
- Many updates to documentation.

**Bug fixes**:

- Throw a more intelligible error during simulation if a topology is produced
  where the time of a parent is equal to the time of the child. (#570, #87).
- Pickle supported in the TableCollection object. (#574, #577).

**Deprecated**:

- The `filter_zero_mutation_sites` parameter for simplify has been deprecated
  in favour of `filter_sites`.

## [0.6.0] - 2018-06-20

This release is focused on ensuring interoperability with the forthcoming SLiM
3.0 release, which has support for outputting tree sequences in msprime's
.trees format. The release represents a substantial step towards the goal of
separating the `tskit` code from `msprime`. It removes the troublesome HDF5
dependency in favour of the much simpler `kastore` library.

The principle new features are the mutate() function which allows us to easily
add mutations to any tree sequence, preliminary support for Individuals and
Populations within the data model, and the addition of the new TableCollection
object as the central structure in the Tables API.

**Breaking changes**:

- Files stored in the HDF5 format will need to upgraded using the
  `msp upgrade` command.

**New features**:

- The mutate function (#507).
- Removed HDF5 library dependency. Now use the embedded kastore library
  for storing data.
- Numpy and h5py are now install time dependencies, solving some installation
  headaches.
- The new TableCollection type  gives much tighter integration with the
  low-level library. Functions like sort_tables and simplify_tables are
  now methods of this class. The load_tables function has been replaced
  by TableCollection.tree_sequence. These functions still work, but are
  deprecated.
- Preliminary support for Individual and Population types in the Tables
  API and for TreeSequences.
- Add 'root' argument to SparseTree.newick and support for arbitrary
  node labels (#510).
- Larger numbers of alleles now supported via 16-bit genotypes (#466).
- Substantially improved simplify performance when there is a large
  number of sites (#453).

**Bug fixes**:

- Fix bug in tree drawing with many roots (#486)
- Fix segfault in accessing trees with zero roots (#515)
- Fix bug where DemographyDebugger was modifying the input sample sizes (#407)

**Deprecated**:

- sort_tables is deprecated in favour of TableCollection.sort().
- simplify_tables is deprecated in favour of TableCollection.simplify().
- load_tables is deprecated in favour of TableCollection.tree_sequence().

## [0.5.0] - 2018-02-26

This is a major update to the underlying data structures in msprime to
generalise the information that can be modelled, and allow
for data from external sources to be efficiently processed. The
new Tables API enables efficient interchange of tree sequence data using
numpy arrays. Many updates have also been made to the tree sequence
API to make it more Pythonic and general. Most changes are backwards
compatible, however.

**Breaking changes**:

- The `SparseTree.mutations()` and `TreeSequence.mutations()` iterators no
  longer support tuple-like access to values. For example, code like

  > > for x, u, j in ts.mutations():
  >
  > : print("mutation at position", x, "node = ", u)
  >
  will no longer work. Code using the old `Mutation.position` and
  `Mutation.index` will still work through deprecated aliases,
  but new code should access these values through `Site.position`
  and `Site.id`, respectively.
- The `TreeSequence.diffs()` method no longer works. Please use
  the `TreeSequence.edge_diffs()` method instead.
- `TreeSequence.get_num_records()` no longer works. Any code using
  this or the `records()` iterator should be rewritten to work with
  the `edges()` iterator and num_edges instead.
- Files stored in the HDF5 format will need to upgraded using the
  `msp upgrade` command.

**New features**:

- The API has been made more Pythonic by replacing (e.g.)
  `tree.get_parent(u)` with `tree.parent(u)`, and
  `tree.get_total_branch_length()` with `tree.total_branch_length`.
  The old forms have been maintained as deprecated aliases. (#64)
- Efficient interchange of tree sequence data using the new Tables
  API. This consists of classes representing the various
  tables (e.g. `NodeTable`) and some utility functions (such
  as `load_tables`, `sort_tables`, etc).
- Support for a much more general class of tree sequence topologies.
  For example, trees with multiple roots are fully supported.
- Substantially generalised mutation model. Mutations now occur at
  specific sites, which can be associated with zero to many mutations.
  Each site has an ancestral state (any character string) and
  each mutation a derived state (any character string).
- Substantially updated documentation to rigorously define the
  underlying data model and requirements for imported data.
- The `variants()` method now returns a list of alleles for each
  site, and genotypes are indexes into this array. This is both
  consistent with existing usage and works with the newly generalised
  mutation model, which allows arbitrary strings of characters as
  mutational states.
- Add the formal concept of a sample, and distinguished from 'leaves'.
  Change `tracked_leaves`, etc. to `tracked_samples` (#225).
  Also rename `sample_size` to `num_samples` for consistency (#227).
- The simplify() method returns subsets of a large tree sequence.
- TreeSequence.first() returns the first tree in sequence.
- Windows support. Msprime is now routinely tested on Windows as
  part of the suite of continuous integration tests.
- Newick output is not supported for more general trees. (#117)
- The `genotype_matrix` method allows efficient access to the
  full genotype matrix. (#306)
- The variants iterator no longer uses a single buffer for
  genotype data, removing a common source of error (#253).
- Unicode and ASCII output formats for `SparseTree.draw()`.
- `SparseTree.draw()` renders tree in the more conventional 'square
  shoulders' format.
- `SparseTree.draw()` by default returns an SVG string, so it can
  be easily displayed in a Jupyter notebook. (#204)
- Preliminary support for a broad class of site-based statistics,
  including Patterson's f-statistics, has been added, through
  the {}`SiteStatCalculator`, and its branch length analog,
  {}`BranchLengthStatCalculator`.  The interface is still in development,
  and is expected may change.

**Bug fixes**:

- Duplicate site no longer possible (#159)
- Fix for incorrect population sizes in DemographyDebugger (#66).

**Deprecated**:

- The `records` iterator has been deprecated, and the underlying data
  model has moved away from the concept of coalescence records. The
  structure of a tree sequence is now defined in terms of a set of nodes
  and edges, essentially a normalised version of coalescence records.
- Changed `population_id` to `population` in various DemographicEvent
  classes for consistency. The old `population_id` argument is kept as a
  deprecated alias.
- Changed `destination` to `dest` in MassMigrationEvent. The old
  `destination` argument is retained as a deprecated alias.
- Changed `sample_size` to `num_samples` in TreeSequence and
  SparseTree. The older versions are retained as deprecated aliases.
- Change `get_num_leaves` to `num_samples` in SparseTree. The
  `get_num_leaves` method (and other related methods) that have
  been retained for backwards compatibility are semantically incorrect,
  in that they now return the number of **samples**. This should have
  no effect on existing code, since samples and leaves were synonymous.
  New code should use the documented `num_samples` form.
- Accessing the `position` attribute on a `Mutation` or
  `Variant` object is now deprecated, as this is a property of a `Site`.
- Accessing the `index` attribute on a `Mutation` or `Variant` object
  is now deprecated. Please use `variant.site.id` instead. In general,
  objects with IDs (i.e., derived from tables) now have an `id` field.
- Various `get_` methods in TreeSequence and SparseTree have been
  replaced by more Pythonic alternatives.

## [0.4.0] - 2016-10-16

Major release providing new functionality and laying groundwork for
upcoming functionality.

**Breaking changes**:

- The HDF5 file format has been changed to allow for non-binary trees
  and to improve performance. It is now both smaller and faster to
  load. However, msprime cannot directly load tree sequence files
  written by older versions. The `msp upgrade` utility has been
  developed to provide an upgrade path for existing users, so that
  files written by older versions of msprime can be converted to the
  newer format and read by version 0.4.x of msprime.
- The tuples returned by the `mutations` method contains an element.
  This will break code doing things like

  > > for pos, node in ts.mutations():
  >
  > : print(pos, node)
  >
  For better forward compatibility, code should use named attributes
  rather than positional access:

  > > for mutation in ts.mutations():
  >
  > : print(mutation.position, mutation.node)
  >
- Similarly, the undocumented `variants` method has some major changes:

  1. The returned tuple has two new values, `node` and `index`
     in the middle of the tuple (but see the point above about using
     named attributes).
  2. The returned genotypes are by default numpy arrays. To revert
     to the old behaviour of returning Python bytes objects, use the
     `as_bytes` argument to the `variants()` method.

**New features**:

- Historical samples. Using the `samples` argument to `simulate`
  users can specify the location and time of all samples explicitly.
- HDF5 file upgrade utility `msp upgrade`
- Support for non-binary trees in the tree sequence, and relaxation
  of the requirements on input tree sequences using the read_txt()
  function.
- Integration with numpy, with zero-copy access to the low-level C API.
- Documented the variants() method that provides access to the sample
  genotypes as either numpy arrays or Python bytes objects.
- New LdCalculator class that allows very fast calculation of r^2 values.
- Initial support for threading.
- The values returned mutations() method now also contain an `index`
  attribute. This makes many operations simpler.
- New TreeSequence.get_time() method that returns the time a sample
  was sampled at.

**Performance improvements**:

- File load times substantially reduced by pre-computing and storing
  traversal indexes.
- O(1) implementation of TreeSequence.get_num_trees()
- Improved control of enabled tree features in TreeSequence.trees()
  method using the `leaf_lists` and `leaf_counts` arguments.

**Bug fixes**:

- Fixed a precision problem with DemographyDebugger. #37
- Segfault on large haplotypes. #29

## [0.3.2] - 2016-07-21

Feature release adding new import and export features to the API
and CLI.

- New `TreeSequence.write_records` and `TreeSequence.write_mutations`
  methods to serialise a tree sequence in a human readable text format.
- New `msprime.load_txt()` method that parses the above formats, and
  allows msprime to read in data from external sources.
- New `TreeSequence.write_vcf` method to write mutation information
  in VCF format.
- Miscellaneous documentation fixes.

## [0.3.1] - 2016-06-24

Feature release adding population related methods to the API.

- New `TreeSequence.get_population(sample_id)` method.
- New `TreeSequence.get_samples(population_id)` method.
- Added the optional `samples` argument to the
  `TreeSequence.get_pairwise_diversity` method.
- Fixed a potential low-level buffer overrun problem.

## [0.3.0] - 2016-05-31

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

## [0.2.0] - 2016-05-05

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

## [0.1.10] - 2016-04-21

Bugfix release. Fixes serious issue affecting simulations with small
sample sizes.

<https://github.com/jeromekelleher/msprime/issues/7>

All users of mspms should update immediately and any analyses using
a small sample size (< 10) with mutations should be repeated.

Many thanks to Konrad Lohse for identifying the issue.

## [0.1.9] - 2016-04-01

Bugfix release. Fixes serious issue affecting random seeds in mspms.

<https://github.com/jeromekelleher/msprime/issues/6>

All users of mspms should update immediately and any analyses using
the `-seeds` option in mspms should be repeated.

Many thanks to Derek Setter for identifying the issue.

## [0.1.8] - 2016-02-17

Transitional release providing population structure support for the
ms-compatible command line interface. A considerable amount of low-level
plumbing was required to provide the required flexibility. This is currently
not visible from the high-level API, but will shortly be made available in the
forthcoming 0.2.x series.

The current implementation of migration should work well for small numbers of
populations (e.g. < 10), but will not scale well for large numbers of
populations.

### Changes

- Added the -I, -m, -ma, -em, -eM, -ema, -eG, -eg, -eN, -en,
  -ej and -es options to mspms. These should provide full ms
  compatibility, except for the -es option which is currently
  limited in scope.
- Added some extra keys to the low-level configuration JSON in
  the HDF5 file format to describe the population structure.
  This will be documented in a future release.
- Added a {}`get_pairwise_diversity` method to the TreeSequence
  class to efficiently calculate the population genetics
  statistic pi.
