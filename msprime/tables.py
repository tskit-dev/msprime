#
# Copyright (C) 2017 University of Oxford
#
# This file is part of msprime.
#
# msprime is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msprime is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msprime.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Tree sequence IO via the tables API.
"""
from __future__ import division
from __future__ import print_function

import base64
import collections
import datetime

from six.moves import copyreg

import _msprime

# We need numpy for most API operations that involve tables. We should make
# numpy a hard requirement soon and get rid of these conditional imports.
try:
    import numpy as np
except ImportError:
    pass


def text_encode_metadata(metadata):
    """
    Returns the specified metadata bytes object encoded as an ASCII-safe
    string.
    """
    return base64.b64encode(metadata).decode('utf8')


def text_decode_metadata(encoded):
    """
    Decodes the specified ASCII encoding of binary metadata and returns the
    resulting bytes.
    """
    return base64.b64decode(encoded.encode('utf8'))


NodeTableRow = collections.namedtuple(
    "NodeTableRow",
    ["flags", "time", "population", "metadata"])


class NodeTable(_msprime.NodeTable):
    """
    A table defining the nodes in a tree sequence. See the
    :ref:`definitions <sec-node-table-definition>` for details on the columns
    in this table and the
    :ref:`tree sequence requirements <sec-valid-tree-sequence-requirements>` section
    for the properties needed for a node table to be a part of a valid tree sequence.

    Example::

        >>> t = msprime.NodeTable()
        >>> t.add_row(flags=0, time=1)
        0
        >>> t.add_row(flags=1, time=2)
        1
        >>> print(t)
        id      flags   population      time    metadata
        0       0       -1      1.00000000000000
        1       1       -1      2.00000000000000
        >>> print(t[0])
        NodeTableRow(flags=0, time=1.0, population=-1, metadata=b'')
        >>> print(t[-1])
        NodeTableRow(flags=1, time=2.0, population=-1, metadata=b'')
        >>> len(t)
        2

    :warning: The numpy arrays returned by table attribute accesses are **copies**
        of the underlying data. In particular, this means that you cannot edit
        the values in the columns by updating the attribute arrays.

        **NOTE:** this behaviour may change in future.

    :ivar time: The array of time values.
    :vartype time: numpy.ndarray, dtype=np.float64
    :ivar flags: The array of flags values.
    :vartype flags: numpy.ndarray, dtype=np.uint32
    :ivar population: The array of population IDs.
    :vartype population: numpy.ndarray, dtype=np.int32
    :ivar metadata: The flattened array of binary metadata values.
    :vartype metadata: numpy.ndarray, dtype=np.int8
    :ivar metadata_offset: The array of offsets into the metadata column.
    :vartype metadata_offset: numpy.ndarray, dtype=np.uint32
    """

    def __str__(self):
        time = self.time
        flags = self.flags
        population = self.population
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tflags\tpopulation\ttime\tmetadata\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\t{}\t{:.14f}\t{}\n".format(
                j, flags[j], population[j], time[j], text_encode_metadata(metadata[j]))
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.flags, other.flags) and
                np.array_equal(self.population, other.population) and
                np.array_equal(self.time, other.time) and
                np.array_equal(self.metadata, other.metadata) and
                np.array_equal(self.metadata_offset, other.metadata_offset))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        return NodeTableRow(*self.get_row(index))

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(
            time=state["time"], flags=state["flags"], population=state["population"],
            metadata=state["metadata"], metadata_offset=state["metadata_offset"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = NodeTable()
        copy.set_columns(
            flags=self.flags, time=self.time, population=self.population,
            metadata=self.metadata, metadata_offset=self.metadata_offset)
        return copy

    def add_row(self, flags=0, time=0, population=-1, metadata=None):
        """
        Adds a new row to this :class:`NodeTable` and returns the ID of the
        corresponding node.

        :param int flags: The bitwise flags for the new node.
        :param float time: The birth time for the new node.
        :param int population: The ID of the population in which the new node was born.
            Defaults to the :const:`.NULL_POPULATION`.
        :param bytes metadata: The binary-encoded metadata for the new node. If not
            specified or None, a zero-length byte string is stored.
        :return: The ID of the newly added node.
        :rtype: int
        """
        return super(NodeTable, self).add_row(flags, time, population, metadata)

    def set_columns(
            self, flags, time, population=None, metadata=None, metadata_offset=None):
        """
        Sets the values for each column in this NodeTable using the values in
        the specified arrays. Overwrites the current state of the table.

        The ``flags``, ``time`` and ``population`` arrays must all be of the same length,
        which is equal to the number of nodes the table will contain. The
        ``metadata`` and ``metadata_offset`` must be supplied together, and
        meet the requirements for :ref:`sec-encoding-ragged-columns`.

        :param flags: The bitwise flags for each node. Required.
        :type flags: numpy.ndarray, dtype=np.uint32
        :param time: The time values for each node. Required.
        :type time: numpy.ndarray, dtype=np.float64
        :param population: The population values for each node. If not specified
            or None, the :const:`.NULL_POPULATION` value is stored for each node.
        :type population: numpy.ndarray, dtype=np.int32
        :param population: The population values for each node. If not specified
            or None, the :const:`.NULL_POPULATION` value is stored for each node.
        :type population: numpy.ndarray, dtype=np.int32
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        super(NodeTable, self).set_columns(
            flags, time, population=population, metadata=metadata,
            metadata_offset=metadata_offset)

    def append_columns(
            self, flags, time, population=None, metadata=None, metadata_offset=None):
        """
        Appends the specified arrays to the columns in this table.

        The ``flags``, ``time`` and ``population`` arrays must all be of the same length,
        which is equal to the number of nodes that will be added to the table. The
        ``metadata`` and ``metadata_offset`` must be supplied together, and
        meet the requirements for :ref:`sec-encoding-ragged-columns`.

        :param flags: The bitwise flags for each node. Required.
        :type flags: numpy.ndarray, dtype=np.uint32
        :param time: The time values for each node. Required.
        :type time: numpy.ndarray, dtype=np.float64
        :param population: The population values for each node. If not specified
            or None, the :const:`.NULL_POPULATION` value is stored for each node.
        :type population: numpy.ndarray, dtype=np.int32
        :param population: The population values for each node. If not specified
            or None, the :const:`.NULL_POPULATION` value is stored for each node.
        :type population: numpy.ndarray, dtype=np.int32
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        super(NodeTable, self).append_columns(
            flags, time, population=population,
            metadata=metadata, metadata_offset=metadata_offset)


# Pickle support. See copyreg registration for this function below.
def _pickle_node_table(table):
    state = {
        "time": table.time,
        "flags": table.flags,
        "population": table.population,
        "metadata": table.metadata,
        "metadata_offset": table.metadata_offset,
    }
    return NodeTable, tuple(), state


class EdgeTable(_msprime.EdgeTable):
    """
    Class for tables describing all edges in a tree sequence, of the form

    left	right	parent	child
    0.0     0.4     3       0
    0.0     0.4     3       2
    0.4     1.0     3       0
    0.4     1.0     3       1
    0.4     1.0     3       2
    0.0     0.4     4       1
    0.0     0.4     4       3

    These describe the half-open genomic interval affected: `[left, right)`,
    the `parent` and the `child` on that interval.

    Requirements: to describe a valid tree sequence, a `EdgeTable` (and
    corresponding `NodeTable`, to provide birth times) must satisfy:

    1. any two edges that share a child must be nonoverlapping, and
    2. the birth times of the `parent` in an edge must be strictly
        greater than the birth times of the `child` in that edge.

    Furthermore, for algorithmic requirements

    4. the smallest `left` coordinate must be 0.0,
    5. the table must be sorted so that birth time of the `parent` increases
        with table row, and
    6. any two edges corresponding to the same `parent` must be
        nonoverlapping.

    It is an additional requirement that the complete ancestry of each sample
    must be specified, but this is harder to verify.

    It is not required that all records corresponding to the same parent be
    adjacent in the table.

    `simplify_tables()` may be used to convert noncontradictory tables
    into tables satisfying the full set of requirements.
    """
    def __str__(self):
        left = self.left
        right = self.right
        parent = self.parent
        child = self.child
        ret = "id\tleft\t\tright\t\tparent\tchild\n"
        for j in range(self.num_rows):
            ret += "{}\t{:.8f}\t{:.8f}\t{}\t{}\n".format(
                j, left[j], right[j], parent[j], child[j])
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.left, other.left) and
                np.array_equal(self.right, other.right) and
                np.array_equal(self.parent, other.parent) and
                np.array_equal(self.child, other.child))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(
            left=state["left"], right=state["right"], parent=state["parent"],
            child=state["child"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = EdgeTable()
        copy.set_columns(
            left=self.left, right=self.right, parent=self.parent, child=self.child)
        return copy


# Pickle support. See copyreg registration for this function below.
def _edge_table_pickle(table):
    state = {
        "left": table.left,
        "right": table.right,
        "parent": table.parent,
        "child": table.child,
    }
    return EdgeTable, tuple(), state


class MigrationTable(_msprime.MigrationTable):
    def __str__(self):
        left = self.left
        right = self.right
        node = self.node
        source = self.source
        dest = self.dest
        time = self.time
        ret = "id\tleft\tright\tnode\tsource\tdest\ttime\n"
        for j in range(self.num_rows):
            ret += "{}\t{:.8f}\t{:.8f}\t{}\t{}\t{}\t{:.8f}\n".format(
                j, left[j], right[j], node[j], source[j], dest[j], time[j])
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.left, other.left) and
                np.array_equal(self.right, other.right) and
                np.array_equal(self.node, other.node) and
                np.array_equal(self.source, other.source) and
                np.array_equal(self.dest, other.dest) and
                np.array_equal(self.time, other.time))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(
            left=state["left"], right=state["right"], node=state["node"],
            source=state["source"], dest=state["dest"], time=state["time"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = MigrationTable()
        copy.set_columns(
            left=self.left, right=self.right, node=self.node, source=self.source,
            dest=self.dest, time=self.time)
        return copy


# Pickle support. See copyreg registration for this function below.
def _migration_table_pickle(table):
    state = {
        "left": table.left,
        "right": table.right,
        "node": table.node,
        "source": table.source,
        "dest": table.dest,
        "time": table.time,
    }
    return MigrationTable, tuple(), state


class SiteTable(_msprime.SiteTable):
    """
    Class for tables describing all sites at which mutations have occurred in a
    tree sequence, of the form

    id	position	ancestral_state
    0	0.1     	0
    1	0.5     	0

    Here ``id`` is not stored directly, but is determined by the row index in
    the table.  ``position`` is the position along the genome, and
    ``ancestral_state`` gives the allele at the root of the tree at that
    position.
    """
    def __str__(self):
        position = self.position
        ancestral_state = unpack_strings(
            self.ancestral_state, self.ancestral_state_offset)
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tposition\tancestral_state\tmetadata\n"
        for j in range(self.num_rows):
            ret += "{}\t{:.8f}\t{}\t{}\n".format(
                j, position[j], ancestral_state[j],
                text_encode_metadata(metadata[j]))
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.position, other.position) and
                np.array_equal(self.ancestral_state, other.ancestral_state) and
                np.array_equal(
                    self.ancestral_state_offset, other.ancestral_state_offset) and
                np.array_equal(self.metadata, other.metadata) and
                np.array_equal(self.metadata_offset, other.metadata_offset))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(
            position=state["position"],
            ancestral_state=state["ancestral_state"],
            ancestral_state_offset=state["ancestral_state_offset"],
            metadata=state["metadata"],
            metadata_offset=state["metadata_offset"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = SiteTable()
        copy.set_columns(
            position=self.position,
            ancestral_state=self.ancestral_state,
            ancestral_state_offset=self.ancestral_state_offset,
            metadata=self.metadata,
            metadata_offset=self.metadata_offset)
        return copy


# Pickle support. See copyreg registration for this function below.
def _site_table_pickle(table):
    state = {
        "position": table.position,
        "ancestral_state": table.ancestral_state,
        "ancestral_state_offset": table.ancestral_state_offset,
        "metadata": table.metadata,
        "metadata_offset": table.metadata_offset,
    }
    return SiteTable, tuple(), state


class MutationTable(_msprime.MutationTable):
    """
    Class for tables describing all mutations that have occurred in a tree
    sequence, of the form

    site	node	derived_state
    0	4	1
    1	3	1
    1	2	0

    Here ``site`` is the index in the SiteTable of the site at which the
    mutation occurred, ``node`` is the index in the NodeTable of the node who
    is the first node inheriting the mutation, and ``derived_state`` is the
    allele resulting from this mutation.

    It is required that mutations occurring at the same node are sorted in
    reverse time order.
    """
    def __str__(self):
        site = self.site
        node = self.node
        parent = self.parent
        derived_state = unpack_strings(self.derived_state, self.derived_state_offset)
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tsite\tnode\tderived_state\tparent\tmetadata\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                j, site[j], node[j], derived_state[j], parent[j],
                text_encode_metadata(metadata[j]))
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.site, other.site) and
                np.array_equal(self.node, other.node) and
                np.array_equal(self.parent, other.parent) and
                np.array_equal(self.derived_state, other.derived_state) and
                np.array_equal(
                    self.derived_state_offset, other.derived_state_offset) and
                np.array_equal(self.metadata, other.metadata) and
                np.array_equal(self.metadata_offset, other.metadata_offset))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(
            site=state["site"], node=state["node"], parent=state["parent"],
            derived_state=state["derived_state"],
            derived_state_offset=state["derived_state_offset"],
            metadata=state["metadata"], metadata_offset=state["metadata_offset"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = MutationTable()
        copy.set_columns(
            site=self.site, node=self.node, parent=self.parent,
            derived_state=self.derived_state,
            derived_state_offset=self.derived_state_offset,
            metadata=self.metadata, metadata_offset=self.metadata_offset)
        return copy


# Pickle support. See copyreg registration for this function below.
def _mutation_table_pickle(table):
    state = {
        "site": table.site,
        "node": table.node,
        "parent": table.parent,
        "derived_state": table.derived_state,
        "derived_state_offset": table.derived_state_offset,
        "metadata": table.metadata,
        "metadata_offset": table.metadata_offset,
    }
    return MutationTable, tuple(), state


class ProvenanceTable(_msprime.ProvenanceTable):
    """
    TODO Document
    """
    def add_row(self, record, timestamp=None):
        """
        Adds a new row to this ProvenanceTable consisting of the specified record and
        timestamp. If timestamp is not specified, it is automatically generated from
        the current time.

        :param str record: A provenance record, describing the parameters and
            environment used to generate the current set of tables.
        :param str timestamp: A string timestamp. This should be in ISO8601 form.
        """
        if timestamp is None:
            timestamp = datetime.datetime.now().isoformat()
        # Note that the order of the positional arguments has been reversed
        # from the low-level module, which is a bit confusing. However, we
        # want the default behaviour here to be to add a row to the table at
        # the current time as simply as possible.
        return super(ProvenanceTable, self).add_row(record=record, timestamp=timestamp)

    def __str__(self):
        timestamp = unpack_strings(self.timestamp, self.timestamp_offset)
        record = unpack_strings(self.record, self.record_offset)
        ret = "id\ttimestamp\trecord\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\t{}\n".format(j, timestamp[j], record[j])
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.timestamp, other.timestamp) and
                np.array_equal(self.timestamp_offset, other.timestamp_offset) and
                np.array_equal(self.record, other.record) and
                np.array_equal(self.record_offset, other.record_offset))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(
            timestamp=state["timestamp"],
            timestamp_offset=state["timestamp_offset"],
            record=state["record"],
            record_offset=state["record_offset"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = ProvenanceTable()
        copy.set_columns(
            timestamp=self.timestamp,
            timestamp_offset=self.timestamp_offset,
            record=self.record,
            record_offset=self.record_offset)
        return copy


# Pickle support. See copyreg registration for this function below.
def _provenance_table_pickle(table):
    state = {
        "timestamp": table.timestamp,
        "timestamp_offset": table.timestamp_offset,
        "record": table.record,
        "record_offset": table.record_offset,
    }
    return ProvenanceTable, tuple(), state


# Pickle support for the various tables. We are forced to use copyreg.pickle
# here to support Python 2. For Python 3, we can just use the __setstate__.
# It would be cleaner to attach the pickle_*_table functions to the classes
# themselves, but this causes issues with Mocking on readthedocs. Sigh.
copyreg.pickle(NodeTable, _pickle_node_table)
copyreg.pickle(EdgeTable, _edge_table_pickle)
copyreg.pickle(MigrationTable, _migration_table_pickle)
copyreg.pickle(SiteTable, _site_table_pickle)
copyreg.pickle(MutationTable, _mutation_table_pickle)
copyreg.pickle(ProvenanceTable, _provenance_table_pickle)


class TableCollection(object):
    """
    A collection of tables. This is a convenience class allowing for convenient
    printing and comparisons of a collection of related tables.
    """
    def __init__(
            self, nodes=None, edges=None, migrations=None, sites=None, mutations=None,
            provenances=None):
        self.nodes = nodes
        self.edges = edges
        self.migrations = migrations
        self.sites = sites
        self.mutations = mutations
        self.provenances = provenances

    def asdict(self):
        """
        Returns this TableCollection as a dictionary mapping the keys "nodes",
        "edges", etc to their respective table objects.
        """
        return {
            "nodes": self.nodes,
            "edges": self.edges,
            "migrations": self.migrations,
            "sites": self.sites,
            "mutations": self.mutations,
            "provenances": self.provenances
        }

    def __banner(self, title):
        width = 60
        line = "#" * width
        title_line = "#   {}".format(title)
        title_line += " " * (width - len(title_line) - 1)
        title_line += "#"
        return line + "\n" + title_line + "\n" + line + "\n"

    def __str__(self):
        s = self.__banner("Nodes")
        s += str(self.nodes) + "\n"
        s += self.__banner("Edges")
        s += str(self.edges) + "\n"
        s += self.__banner("Sites")
        s += str(self.sites) + "\n"
        s += self.__banner("Mutations")
        s += str(self.mutations) + "\n"
        s += self.__banner("Migrations")
        s += str(self.migrations) + "\n"
        s += self.__banner("Provenances")
        s += str(self.provenances)
        return s

    # TODO add support for __eq__ and __ne__


#############################################
# Table functions.
#############################################


def sort_tables(*args, **kwargs):
    """
    Sorts the given tables in place, as follows:

    Edges are ordered by

    - time of parent, then
    - parent node ID, then
    - child node ID, then
    - left endpoint.

    Sites are ordered by position, and Mutations are ordered by site.

    If the ``edge_start`` parameter is provided, this specifies the index
    in the edge table where sorting should start. Only rows with index
    greater than or equal to ``edge_start`` are sorted; rows before this index
    are not affected. This parameter is provided to allow for efficient sorting
    when the user knows that the edges up to a given index are already sorted.

    This function is more general than ``TreeSequence.simplify()``, since it can
    be applied to tables not satisfying the tree sequence
    :ref:`sec-ordering-criteria`
    (and that hence could not be loaded into a TreeSequence).

    :param NodeTable nodes: The tree sequence nodes (required).
    :param EdgeTable edges: The tree sequence edges (required).
    :param MigrationTable migrations: The tree sequence migrations (optional).
    :param SiteTable sites: The tree sequence sites (optional, but required if
         ``mutations`` is provided)
    :param MutationTable mutations: The tree sequence mutations (optional, but
         required if ``sites`` is provided).
    :param int edge_start: The index in the edge table where sorting starts
        (default=0; must be <= len(edges)).
    """
    kwargs_copy = dict(kwargs)
    # If provenances is supplied as a keyword argument just ignore it. This is
    # because we'll often call sort_tables(**t.asdict()), and the provenances
    # entry breaks this pattern.
    kwargs_copy.pop("provenances", None)
    return _msprime.sort_tables(*args, **kwargs_copy)


def simplify_tables(*args, **kwargs):
    """
    Simplifies the tables, in place, to retain only the information necessary
    to reconstruct the tree sequence describing the given ``samples``.  This
    will change the ID of the nodes, so that the individual ``samples[k]]``
    will have ID ``k`` in the result. The resulting NodeTable will have only
    the first ``len(samples)`` individuals marked as samples. The
    ``sequence_length`` can be provided but is otherwise inferred from the
    largest right edge. The mapping from node IDs in the current set of tables
    to their equivalent values in the simplified tables is returned as a numpy
    array. If an array ``a`` is returned by this function and ``u`` is the ID of
    a node in the input table, then ``a[u]`` is the ID of this node in the
    output table. For any node ``u`` that is not mapped into the output tables,
    this mapping will equal ``-1``.

    Tables operated on by this function must: be sorted (see ``sort_tables``),
    have children be born strictly after their parents, and the intervals on
    which any individual is a child must be disjoint; but other than this the
    tables need not satisfy remaining requirements to specify a valid tree
    sequence (but the resulting tables will).

    :param list[int] samples: A list of Node IDs of individuals to retain as samples.
    :param NodeTable nodes: The NodeTable to be simplified.
    :param EdgeTable edges: The EdgeTable to be simplified.
    :param MigrationTable migrations: The MigrationTable to be simplified.
    :param SiteTable sites: The SiteTable to be simplified.
    :param MutationTable mutations: The MutationTable to be simplified.
    :param bool filter_invariant_sites: Whether to remove sites that have no
        mutations from the output (default: True).
    :param float sequence_length: The length of the sequence.
    :return: A numpy array mapping node IDs in the input tables to their
        corresponding node IDs in the output tables.
    :rtype: numpy array (dtype=np.int32).
    """
    return _msprime.simplify_tables(*args, **kwargs)


def pack_bytes(data):
    """
    Packs the specified list of bytes into a flattened numpy array of 8 bit integers
    and corresponding lengths.
    """
    n = len(data)
    offsets = np.zeros(n + 1, dtype=np.uint32)
    for j in range(n):
        offsets[j + 1] = offsets[j] + len(data[j])
    column = np.zeros(offsets[-1], dtype=np.int8)
    for j, value in enumerate(data):
        column[offsets[j]: offsets[j + 1]] = bytearray(value)
    return column, offsets


def unpack_bytes(packed, offset):
    """
    Unpacks a list of bytes from the specified numpy arrays of packed byte
    data and corresponding offsets.
    """
    # This could be done a lot more efficiently...
    ret = []
    for j in range(offset.shape[0] - 1):
        raw = packed[offset[j]: offset[j + 1]].tobytes()
        ret.append(raw)
    return ret


def pack_strings(strings, encoding="utf8"):
    """
    Packs the specified list of strings into a flattened numpy array of 8 bit integers
    and corresponding lengths.
    """
    return pack_bytes([bytearray(s.encode(encoding)) for s in strings])


def unpack_strings(packed, offset, encoding="utf8"):
    """
    Unpacks a list of strings from the specified numpy arrays of packed byte
    data and corresponding offsets.
    """
    return [b.decode(encoding) for b in unpack_bytes(packed, offset)]
