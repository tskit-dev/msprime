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
