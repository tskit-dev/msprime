--------------------
[1.0.0] - 2022-05-24
--------------------

This major release marks the point at which the documented API becomes stable and supported.

**Breaking changes**

- Change the type of genotypes to ``int32_t``, removing the TSK_16_BIT_GENOTYPES flag option.
  (:user:`benjeffery`, :issue:`463`, :pr:`2108`)

- ``tsk_variant_t`` now includes its ``tsk_site_t`` rather than pointing to it.
  (:user:`benjeffery`, :issue:`2161`, :pr:`2162`)

- Rename ``TSK_TAKE_TABLES`` to ``TSK_TAKE_OWNERSHIP``.
  (:user:`benjeffery`, :issue:`2221`, :pr:`2222`)

- ``TSK_DEBUG``, ``TSK_NO_INIT``, ``TSK_NO_CHECK_INTEGRITY`` and ``TSK_TAKE_OWNERSHIP`` have moved to ``core.h``
  (:user:`benjeffery`, :issue:`2218`, :pr:`2230`))

- Rename several flags:
     - All flags to ``simplify`` for example ``TSK_KEEP_INPUT_ROOTS`` becomes ``TSK_SIMPLIFY_KEEP_INPUT_ROOTS``.
     - All flags to ``subset`` for example ``TSK_KEEP_UNREFERENCED`` becomes ``TSK_SUBSET_KEEP_UNREFERENCED``.
     - ``TSK_BUILD_INDEXES`` -> ``TSK_TS_INIT_BUILD_INDEXES``
     - ``TSK_NO_METADATA`` -> ``TSK_TABLE_NO_METADATA``
     - ``TSK_NO_EDGE_METADATA`` -> ``TSK_TC_NO_EDGE_METADATA``

  (:user:`benjeffery`, :issue:`1720`, :pr:`2226`, :pr:`2229`, :pr:`2224`)

- Remove the generic ``TSK_ERR_OUT_OF_BOUNDS`` - replacing with specific errors.
  Remove ``TSK_ERR_NON_SINGLE_CHAR_MUTATION`` which was unused.
  (:user:`benjeffery`, :pr:`2260`)

- Reorder stats API methods to place ``result`` as the last argument. (:user:`benjeffery`, :pr:`2292`, :issue:`2285`)

**Features**

- Make dumping of tables and tree sequences to disk a zero-copy operation.
  (:user:`benjeffery`, :issue:`2111`, :pr:`2124`)

- Add ``edge`` attribute to ``mutation_t`` struct and make available in tree sequence.
  (:user:`jeromekelleher`, :issue:`685`, :pr:`2279`)

- Reduce peak memory usage in ``tsk_treeseq_simplify``.
  (:user:`jeromekelleher`, :issue:`2287`, :pr:`2288`)

----------------------
[0.99.15] - 2021-12-07
----------------------

**Breaking changes**

- The ``tables`` argument to ``tsk_treeseq_init`` is no longer ``const``, to allow for future no-copy tree sequence creation.
  (:user:`benjeffery`, :issue:`1718`, :pr:`1719`)
- Additional consistency checks for mutation tables are now run by ``tsk_table_collection_check_integrity``
  even when ``TSK_CHECK_MUTATION_ORDERING`` is not passed in. (:user:`petrelharp`, :issue:`1713`, :pr:`1722`)

- ``num_tracked_samples`` and ``num_samples`` in ``tsk_tree_t`` are now typed as ``tsk_size_t``
  (:user:`benjeffery`, :issue:`1723`, :pr:`1727`)

- The previously deprecated option ``TSK_SAMPLE_COUNTS`` has been removed. (:user:`benjeffery`, :issue:`1744`, :pr:`1761`).
- Individuals are no longer guaranteed or required to be topologically sorted in a tree sequence.
  ``tsk_table_collection_sort`` no longer sorts individuals.
  (:user:`benjeffery`, :issue:`1774`, :pr:`1789`)

- The ``tsk_tree_t.left_root`` member has been removed. Client code can be updated
  most easily by using the equivalent ``tsk_tree_get_left_root`` function. However,
  it may be worth considering updating code to use either the standard traversal
  functions (which automatically iterate over roots) or to use the ``virtual_root``
  member (which may lead to more concise code). (:user:`jeromekelleher`, :issue:`1796`,
  :pr:`1862`)

- Rename ``tsk_tree_t.left`` and ``tsk_tree_t.right`` members to
  ``tsk_tree_t.interval.left`` and ``tsk_tree_t.interval.right`` respectively.
  (:user:`jeromekelleher`, :issue:`1686`, :pr:`1913`)

- ``kastore`` is now vendored into this repo instead of being a git submodule. Developers need to run
  ``git submodule update``. (:user:`jeromekelleher`, :issue:`1687`, :pr:`1973`)

- ``Tree`` arrays such as ``left_sib``, ``right_child`` etc. now have an additional
  "virtual root" node at the end. (:user:`jeromekelleher`, :issue:`1691`, :pr:`1704`)

- ``marked`` and ``mark`` have been removed from ``tsk_tree_t``. (:user:`jeromekelleher`, :pr:`1936`)

**Features**

- Add ``tsk_table_collection_individual_topological_sort`` to sort the individuals as this is no longer done by the
  default sort. (:user:`benjeffery`, :issue:`1774`, :pr:`1789`)

- The default behaviour for table size growth is now to double the current size of the table,
  up to a threshold. To keep the previous behaviour, use (e.g.)
  ``tsk_edge_table_set_max_rows_increment(tables->edges, 1024)``, which results in adding
  space for 1024 additional rows each time we run out of space in the edge table.
  (:user:`benjeffery`, :issue:`5`, :pr:`1683`)
- ``tsk_table_collection_check_integrity`` now has a ``TSK_CHECK_MIGRATION_ORDERING`` flag. (:user:`petrelharp`, :pr:`1722`)

- The default behaviour for ragged column growth is now to double the current size of the column,
  up to a threshold. To keep the previous behaviour, use (e.g.)
  ``tsk_node_table_set_max_metadata_length_increment(tables->nodes, 1024)``, which results in adding
  space for 1024 additional entries each time we run out of space in the ragged column.
  (:user:`benjeffery`, :issue:`1703`, :pr:`1709`)

- Support for compiling the C library on Windows using msys2 (:user:`jeromekelleher`,
  :pr:`1742`).

- Add ``time_units`` to ``tsk_table_collection_t`` to describe the units of the time dimension of the
  tree sequence. This is then used to geerate an error if ``time_units`` is ``uncalibrated`` when
  using the branch lengths in statistics. (:user:`benjeffery`, :issue:`1644`, :pr:`1760`)

- Add the ``TSK_LOAD_SKIP_TABLES`` option to load just the top-level information from a
  file. Also add the ``TSK_CMP_IGNORE_TABLES`` option to compare only the top-level
  information in two table collections. (:user:`clwgg`, :pr:`1882`, :issue:`1854`).

- Add reference sequence.
  (:user:`jeromekelleher`, :user:`benjeffery`, :issue:`146`, :pr:`1911`, :pr:`1944`, :pr:`1911`)

- Add the ``TSK_LOAD_SKIP_REFERENCE_SEQUENCE`` option to load a table collection
  without the reference sequence. Also add the TSK_CMP_IGNORE_REFERENCE_SEQUENCE
  option to compare two table collections without comparing their reference
  sequence. (:user:`clwgg`, :pr:`2019`, :issue:`1971`).

- Add a "virtual root" to ``Tree`` arrays such as ``left_sib``, ``right_child`` etc.
  The virtual root is appended to each array, has all real roots as its children,
  but is not the parent of any node. Simplifies traversal algorithms.
  (:user:`jeromekelleher`, :issue:`1691`, :pr:`1704`)

- Add ``num_edges`` to ``tsk_tree_t`` to count the edges that define the topology of
  the tree. (:user:`jeromekelleher`, :pr:`1704`)

- Add the ``tsk_tree_get_size_bound`` function which returns an upper bound on the number of nodes reachable from
  the roots of a tree. Useful for tree stack allocations (:user:`jeromekelleher`, :pr:`1704`).

- Add ``MetadataSchema.permissive_json`` for an easy way to get the simplest schema.


----------------------
[0.99.14] - 2021-09-03
----------------------

**Breaking changes**

- 64 bits are now used to store the sizes of ragged table columns such as metadata,
  allowing them to hold more data. As such ``tsk_size_t`` is now 64 bits wide.
  This change is fully backwards and forwards compatible for all tree-sequences whose
  ragged column sizes fit into 32 bits. New tree-sequences with
  large offset arrays that require 64 bits will fail to load in previous versions with
  error ``TSK_ERR_BAD_COLUMN_TYPE``.
  (:user:`jeromekelleher`, :issue:`343`, :issue:`1527`, :issue:`1528`, :issue:`1530`,
  :issue:`1554`, :issue:`1573`, :issue:`1589`,:issue:`1598`,:issue:`1628`, :pr:`1571`,
  :pr:`1579`, :pr:`1585`, :pr:`1590`, :pr:`1602`, :pr:`1618`, :pr:`1620`, :pr:`1652`).

**Features**

- Add `tsk_X_table_update_row` methods which allow modifying single rows of tables
  (:user:`jeromekelleher`, :issue:`1545`, :pr:`1552`).

----------------------
[0.99.13] - 2021-07-08
----------------------
**Fixes**

- Fix segfault when very large columns overflow
  (:user:`bhaller`, :user:`benjeffery`, :issue:`1509`, :pr:`1511`).

----------------------
[0.99.12] - 2021-05-14
----------------------

**Breaking changes**

- Removed ``TSK_NO_BUILD_INDEXES``.
  Not building indexes is now the default behaviour of `tsk_table_collection_dump` and related functions.
  (:user:`molpopgen`, :issue:`1327`, :pr:`1337`).

**Features**

- Add ``tsk_*_table_extend`` methods to append to a table from another
  (:user:`benjeffery`, :issue:`1271`, :pr:`1287`).

**Fixes**

----------------------
[0.99.11] - 2021-03-16
----------------------

**Features**

- Add ``parents`` to the individual table to enable recording of pedigrees
  (:user:`ivan-krukov`, :user:`benjeffery`, :issue:`852`, :pr:`1125`, :pr:`866`, :pr:`1153`, :pr:`1177`, :pr:`1199`).

- Added a ``tsk_table_collection_canonicalise`` method, that allows checking for equality between
  tables that are equivalent up to reordering (:user:`petrelharp`, :user:`mufernando`, :pr:`1108`).

- Removed a previous requirement on ``tsk_table_collection_union``, allowing for unioning of
  new information both above and below shared history (:user:`petrelharp`, :user:`mufernando`, :pr:`1108`).

- Support migrations in tsk_table_collection_sort. (:user:`jeromekelleher`,
  :issue:`22`, :issue:`117`, :pr:`1131`).

**Breaking changes**

- Method ``tsk_individual_table_add_row`` has an extra arguments ``parents`` and ``parents_length``.

- Add an ``options`` argument to ``tsk_table_collection_subset`` (:user:`petrelharp`, :pr:`1108`),
  to allow for retaining the order of populations.

- Mutation error codes have changed

**Changes**

- Allow mutations that have the same derived state as their parent mutation.
  (:user:`benjeffery`, :issue:`1180`, :pr:`1233`)

- File minor version change to support individual parents

----------------------
[0.99.10] - 2021-01-25
----------------------

Minor bugfix on internal APIs

---------------------
[0.99.9] - 2021-01-22
---------------------

**Features**

- Add ``TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS`` flag to simplify, which allows the user to
  keep unary nodes only if they belong to a tabled individual. This is useful for
  simplification in forwards simulations (:user:`hyanwong`, :issue:`1113`, :pr:`1119`).


---------------------
[0.99.8] - 2020-11-27
---------------------

**Features**

- Add ``tsk_treeseq_genetic_relatedness`` for calculating genetic relatedness between
  pairs of sets of nodes (:user:`brieuclehmann`, :issue:`1021`, :pr:`1023`, :issue:`974`,
  :issue:`973`, :pr:`898`).

- Exposed ``tsk_table_collection_set_indexes`` to the API
  (:user:`benjeffery`, :issue:`870`, :pr:`921`).

**Breaking changes**

- Added an ``options`` argument to ``tsk_table_collection_equals``
  and table equality methods to allow for more flexible equality criteria
  (e.g., ignore top-level metadata and schema or provenance tables).
  Existing code should add an extra final parameter ``0`` to retain the
  current behaviour (:user:`mufernando`, :user:`jeromekelleher`,
  :issue:`896`, :pr:`897`, :issue:`913`, :pr:`917`).

- Changed default behaviour of ``tsk_table_collection_clear`` to not clear
  provenances and added ``options`` argument to optionally clear provenances
  and schemas (:user:`benjeffery`, :issue:`929`, :pr:`1001`).

- Renamed ``ts.trait_regression`` to ``ts.trait_linear_model``.

---------------------
[0.99.7] - 2020-09-29
---------------------

- Added ``TSK_INCLUDE_TERMINAL`` option to ``tsk_diff_iter_init`` to output the last edges
  at the end of a tree sequence (:user:`hyanwong`, :issue:`783`, :pr:`787`).

- Added ``tsk_bug_assert`` for assertions that should be compiled into release binaries
  (:user:`benjeffery`, :pr:`860`).

---------------------
[0.99.6] - 2020-09-04
---------------------

**Bugfixes**

- :issue:`823` - Fix mutation time error when using
  ``tsk_table_collection_simplify`` with ``TSK_SIMPLIFY_KEEP_INPUT_ROOTS``
  (:user:`petrelharp`, :pr:`823`).

---------------------
[0.99.5] - 2020-08-27
---------------------

**Breaking changes**

- The macro ``TSK_IMPUTE_MISSING_DATA`` is renamed to ``TSK_ISOLATED_NOT_MISSING``
  (:user:`benjeffery`, :issue:`716`, :pr:`794`)

**New features**

- Add a ``TSK_SIMPLIFY_KEEP_INPUT_ROOTS`` option to simplify which, if enabled, adds edges
  from the MRCAs of samples in the simplified tree sequence back to the roots
  in the input tree sequence (:user:`jeromekelleher`, :issue:`775`, :pr:`782`).

**Bugfixes**

- :issue:`777` - Mutations over isolated samples were incorrectly decoded as
  missing data. (:user:`jeromekelleher`, :pr:`778`)

- :issue:`776` - Fix a segfault when a partial list of samples
  was provided to the ``variants`` iterator. (:user:`jeromekelleher`, :pr:`778`)

---------------------
[0.99.4] - 2020-08-12
---------------------

**Note**

- The ``TSK_VERSION_PATCH`` macro was incorrectly set to ``4`` for 0.99.3, so both
  0.99.4 and 0.99.3 have the same value.

**Changes**

- Mutation times can be a mixture of known and unknown as long as for each
  individual site  they are either all known or all unknown (:user:`benjeffery`, :pr:`761`).

**Bugfixes**

- Fix for including core.h under C++ (:user:`petrelharp`, :pr:`755`).

---------------------
[0.99.3] - 2020-07-27
---------------------

**Breaking changes**

- ``tsk_mutation_table_add_row`` has an extra ``time`` argument. If the time
  is unknown ``TSK_UNKNOWN_TIME`` should be passed.
  (:user:`benjeffery`, :pr:`672`)

- Change genotypes from unsigned to signed to accommodate missing data
  (see :issue:`144` for discussion). This only affects users of the
  ``tsk_vargen_t`` class. Genotypes are now stored as int8_t and int16_t
  types rather than the former unsigned types. The field names in the
  genotypes union of the ``tsk_variant_t`` struct returned by ``tsk_vargen_next``
  have been renamed to ``i8`` and ``i16`` accordingly; care should be
  taken when updating client code to ensure that types are correct. The number
  of distinct alleles supported by 8 bit genotypes has therefore dropped
  from 255 to 127, with a similar reduction for 16 bit genotypes.

- Change the ``tsk_vargen_init`` method to take an extra parameter ``alleles``.
  To keep the current behaviour, set this parameter to NULL.

- Edges can now have metadata. Hence edge methods now take two extra arguments:
  metadata and metadata length. The file format has also changed to accommodate this,
  but is backwards compatible. Edge metadata can be disabled for a table collection with
  the TSK_NO_EDGE_METADATA flag.
  (:user:`benjeffery`, :pr:`496`, :pr:`712`)

- Migrations can now have metadata. Hence migration methods now take two extra arguments:
  metadata and metadata length. The file format has also changed to accommodate this,
  but is backwards compatible.
  (:user:`benjeffery`, :pr:`505`)

- The text dump of tables with metadata now includes the metadata schema as a header.
  (:user:`benjeffery`, :pr:`493`)

- Bad tree topologies are detected earlier, so that it is no longer possible
  to create a tsk_treeseq_t object which contains a parent with contradictory
  children on an interval. Previously an error occured when some operation
  building the trees was attempted (:user:`jeromekelleher`, :pr:`709`).

**New features**

- New methods to perform set operations on table collections.
  ``tsk_table_collection_subset`` subsets and reorders table collections by nodes
  (:user:`mufernando`, :user:`petrelharp`, :pr:`663`, :pr:`690`).
  ``tsk_table_collection_union`` forms the node-wise union of two table collections
  (:user:`mufernando`, :user:`petrelharp`, :issue:`381`, :pr:`623`).

- Mutations now have an optional double-precision floating-point ``time`` column.
  If not specified, this defaults to a particular NaN value (``TSK_UNKNOWN_TIME``)
  indicating that the time is unknown. For a tree sequence to be considered valid
  it must meet new criteria for mutation times, see :ref:`sec_mutation_requirements`.
  Add ``tsk_table_collection_compute_mutation_times`` and new flag to
  ``tsk_table_collection_check_integrity``:``TSK_CHECK_MUTATION_TIME``. Table sorting
  orders mutations by non-increasing time per-site, which is also a requirement for a
  valid tree sequence.
  (:user:`benjeffery`, :pr:`672`)

- Add ``metadata`` and ``metadata_schema`` fields to table collection, with accessors on
  tree sequence. These store arbitrary bytes and are optional in the file format.
  (:user: `benjeffery`, :pr:`641`)

- Add the ``TSK_SIMPLIFY_KEEP_UNARY`` option to simplify (:user:`gtsambos`). See :issue:`1`
  and :pr:`143`.

- Add a ``set_root_threshold`` option to tsk_tree_t which allows us to set the
  number of samples a node must be an ancestor of to be considered a root
  (:pr:`462`).

- Change the semantics of tsk_tree_t so that sample counts are always
  computed, and add a new ``TSK_NO_SAMPLE_COUNTS`` option to turn this
  off (:pr:`462`).

- Tables with metadata now have an optional `metadata_schema` field that can contain
  arbitrary bytes. (:user:`benjeffery`, :pr:`493`)

- Tables loaded from a file can now be edited in the same way as any other
  table collection (:user:`jeromekelleher`, :issue:`536`, :pr:`530`.

- Support for reading/writing to arbitrary file streams with the loadf/dumpf
  variants for tree sequence and table collection load/dump
  (:user:`jeromekelleher`, :user:`grahamgower`, :issue:`565`, :pr:`599`).

- Add low-level sorting API and ``TSK_NO_CHECK_INTEGRITY`` flag
  (:user:`jeromekelleher`, :pr:`627`, :issue:`626`).

- Add extension of Kendall-Colijn tree distance metric for tree sequences
  computed by ``tsk_treeseq_kc_distance``
  (:user:`daniel-goldstein`, :pr:`548`)

**Deprecated**

- The ``TSK_SAMPLE_COUNTS`` options is now ignored and  will print out a warning
  if used (:pr:`462`).

---------------------
[0.99.2] - 2019-03-27
---------------------

Bugfix release. Changes:

- Fix incorrect errors on tbl_collection_dump (#132)
- Catch table overflows (#157)

---------------------
[0.99.1] - 2019-01-24
---------------------

Refinements to the C API as we move towards 1.0.0. Changes:

- Change the ``_tbl_`` abbreviation to ``_table_`` to improve readability.
  Hence, we now have, e.g., ``tsk_node_table_t`` etc.
- Change ``tsk_tbl_size_t`` to ``tsk_size_t``.
- Standardise public API to use ``tsk_size_t`` and ``tsk_id_t`` as appropriate.
- Add ``tsk_flags_t`` typedef and consistently use this as the type used to
  encode bitwise flags. To avoid confusion, functions now have an ``options``
  parameter.
- Rename ``tsk_table_collection_position_t`` to ``tsk_bookmark_t``.
- Rename ``tsk_table_collection_reset_position`` to ``tsk_table_collection_truncate``
  and ``tsk_table_collection_record_position`` to ``tsk_table_collection_record_num_rows``.
- Generalise ``tsk_table_collection_sort`` to take a bookmark as start argument.
- Relax restriction that nodes in the ``samples`` argument to simplify must
  currently be marked as samples. (https://github.com/tskit-dev/tskit/issues/72)
- Allow ``tsk_table_collection_simplify`` to take a NULL samples argument to
  specify "all samples in the current tables".
- Add support for building as a meson subproject.

---------------------
[0.99.0] - 2019-01-14
---------------------

Initial alpha version of the tskit C API tagged. Version 0.99.x
represents the series of releases leading to version 1.0.0 which
will be the first stable release. After 1.0.0, semver rules
regarding API/ABI breakage will apply; however, in the 0.99.x
series arbitrary changes may happen.

--------------------
[0.0.0] - 2019-01-10
--------------------

Initial extraction of tskit code from msprime. Relicense to MIT.
Code copied at hash 29921408661d5fe0b1a82b1ca302a8b87510fd23
