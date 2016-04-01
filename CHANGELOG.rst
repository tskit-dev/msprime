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
