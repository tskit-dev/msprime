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


EdgeTableRow = collections.namedtuple(
    "EdgeTableRow",
    ["left", "right", "parent", "child"])


MigrationTableRow = collections.namedtuple(
    "MigrationTableRow",
    ["left", "right", "node", "source", "dest", "time"])


SiteTableRow = collections.namedtuple(
    "SiteTableRow",
    ["position", "ancestral_state", "metadata"])


MutationTableRow = collections.namedtuple(
    "MutationTableRow",
    ["site", "node", "derived_state", "parent", "metadata"])


ProvenanceTableRow = collections.namedtuple(
    "ProvenanceTableRow",
    ["timestamp", "record"])


class NodeTable(_msprime.NodeTable):
    """
    A table defining the nodes in a tree sequence. See the
    :ref:`definitions <sec_node_table_definition>` for details on the columns
    in this table and the
    :ref:`tree sequence requirements <sec_valid_tree_sequence_requirements>` section
    for the properties needed for a node table to be a part of a valid tree sequence.

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
    :ivar metadata: The flattened array of binary metadata values. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata: numpy.ndarray, dtype=np.int8
    :ivar metadata_offset: The array of offsets into the metadata column. See
        :ref:`sec_tables_api_binary_columns` for more details.
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

    def clear(self):
        """
        Deletes all rows in this table.
        """
        super(NodeTable, self).clear()

    def reset(self):
        # Deprecated alias for clear
        self.clear()

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
        Sets the values for each column in this :class:`.NodeTable` using the values in
        the specified arrays. Overwrites any data currently stored in the table.

        The ``flags``, ``time`` and ``population`` arrays must all be of the same length,
        which is equal to the number of nodes the table will contain. The
        ``metadata`` and ``metadata_offset`` parameters must be supplied together, and
        meet the requirements for :ref:`sec_encoding_ragged_columns`.
        See :ref:`sec_tables_api_binary_columns` for more information.

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
        Appends the specified arrays to the end of the columns in this
        :class:`NodeTable`. This allows many new rows to be added at once.

        The ``flags``, ``time`` and ``population`` arrays must all be of the same length,
        which is equal to the number of nodes that will be added to the table. The
        ``metadata`` and ``metadata_offset`` parameters must be supplied together, and
        meet the requirements for :ref:`sec_encoding_ragged_columns`.
        See :ref:`sec_tables_api_binary_columns` for more information.

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
    A table defining the edges in a tree sequence. See the
    :ref:`definitions <sec_edge_table_definition>` for details on the columns
    in this table and the
    :ref:`tree sequence requirements <sec_valid_tree_sequence_requirements>` section
    for the properties needed for an edge table to be a part of a valid tree sequence.

    :warning: The numpy arrays returned by table attribute accesses are **copies**
        of the underlying data. In particular, this means that you cannot edit
        the values in the columns by updating the attribute arrays.

        **NOTE:** this behaviour may change in future.

    :ivar left: The array of left coordinates.
    :vartype left: numpy.ndarray, dtype=np.float64
    :ivar right: The array of right coordinates.
    :vartype right: numpy.ndarray, dtype=np.float64
    :ivar parent: The array of parent node IDs.
    :vartype parent: numpy.ndarray, dtype=np.int32
    :ivar child: The array of child node IDs.
    :vartype child: numpy.ndarray, dtype=np.int32
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

    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        return EdgeTableRow(*self.get_row(index))

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

    def clear(self):
        """
        Deletes all rows in this table.
        """
        super(EdgeTable, self).clear()

    def reset(self):
        # Deprecated alias for clear
        self.clear()

    def add_row(self, left, right, parent, child):
        """
        Adds a new row to this :class:`EdgeTable` and returns the ID of the
        corresponding edge.

        :param float left: The left coordinate (inclusive).
        :param float right: The right coordinate (exclusive).
        :param int parent: The ID of parent node.
        :param int child: The ID of child node.
        :return: The ID of the newly added edge.
        :rtype: int
        """
        return super(EdgeTable, self).add_row(left, right, parent, child)

    def set_columns(self, left, right, parent, child):
        """
        Sets the values for each column in this :class:`.EdgeTable` using the values
        in the specified arrays. Overwrites any data currently stored in the table.

        All four parameters are mandatory, and must be numpy arrays of the
        same length (which is equal to the number of edges the table will contain).

        :param left: The left coordinates (inclusive).
        :type left: numpy.ndarray, dtype=np.float64
        :param right: The right coordinates (exclusive).
        :type right: numpy.ndarray, dtype=np.float64
        :param parent: The parent node IDs.
        :type parent: numpy.ndarray, dtype=np.int32
        :param child: The child node IDs.
        :type child: numpy.ndarray, dtype=np.int32
        """
        super(EdgeTable, self).set_columns(left, right, parent, child)

    def append_columns(self, left, right, parent, child):
        """
        Appends the specified arrays to the end of the columns of this
        :class:`EdgeTable`. This allows many new rows to be added at once.

        All four parameters are mandatory, and must be numpy arrays of the
        same length (which is equal to the number of additional edges to
        add to the table).

        :param left: The left coordinates (inclusive).
        :type left: numpy.ndarray, dtype=np.float64
        :param right: The right coordinates (exclusive).
        :type right: numpy.ndarray, dtype=np.float64
        :param parent: The parent node IDs.
        :type parent: numpy.ndarray, dtype=np.int32
        :param child: The child node IDs.
        :type child: numpy.ndarray, dtype=np.int32
        """
        super(EdgeTable, self).append_columns(left, right, parent, child)


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
    """
    A table defining the migrations in a tree sequence. See the
    :ref:`definitions <sec_migration_table_definition>` for details on the columns
    in this table and the
    :ref:`tree sequence requirements <sec_valid_tree_sequence_requirements>` section
    for the properties needed for a migration table to be a part of a valid tree
    sequence.

    :warning: The numpy arrays returned by table attribute accesses are **copies**
        of the underlying data. In particular, this means that you cannot edit
        the values in the columns by updating the attribute arrays.

        **NOTE:** this behaviour may change in future.

    :ivar left: The array of left coordinates.
    :vartype left: numpy.ndarray, dtype=np.float64
    :ivar right: The array of right coordinates.
    :vartype right: numpy.ndarray, dtype=np.float64
    :ivar node: The array of node IDs.
    :vartype node: numpy.ndarray, dtype=np.int32
    :ivar source: The array of source population IDs.
    :vartype source: numpy.ndarray, dtype=np.int32
    :ivar dest: The array of destination population IDs.
    :vartype dest: numpy.ndarray, dtype=np.int32
    :ivar time: The array of time values.
    :vartype time: numpy.ndarray, dtype=np.float64
    """

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

    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        return MigrationTableRow(*self.get_row(index))

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

    def clear(self):
        """
        Deletes all rows in this table.
        """
        super(MigrationTable, self).clear()

    def reset(self):
        # Deprecated alias for clear
        self.clear()

    def add_row(self, left, right, node, source, dest, time):
        """
        Adds a new row to this :class:`MigrationTable` and returns the ID of the
        corresponding migration.

        :param float left: The left coordinate (inclusive).
        :param float right: The right coordinate (exclusive).
        :param int node: The node ID.
        :param int source: The ID of the source population.
        :param int dest: The ID of the destination population.
        :param float time: The time of the migration event.
        :return: The ID of the newly added migration.
        :rtype: int
        """
        return super(MigrationTable, self).add_row(left, right, node, source, dest, time)

    def set_columns(self, left, right, node, source, dest, time):
        """
        Sets the values for each column in this :class:`.MigrationTable` using the values
        in the specified arrays. Overwrites any data currently stored in the table.

        All six parameters are mandatory, and must be numpy arrays of the
        same length (which is equal to the number of migrations the table will contain).

        :param left: The left coordinates (inclusive).
        :type left: numpy.ndarray, dtype=np.float64
        :param right: The right coordinates (exclusive).
        :type right: numpy.ndarray, dtype=np.float64
        :param node: The node IDs.
        :type node: numpy.ndarray, dtype=np.int32
        :param source: The source population IDs.
        :type source: numpy.ndarray, dtype=np.int32
        :param dest: The destination population IDs.
        :type dest: numpy.ndarray, dtype=np.int32
        :param time: The time of each migration.
        :type time: numpy.ndarray, dtype=np.int64
        """
        super(MigrationTable, self).set_columns(left, right, node, source, dest, time)

    def append_columns(self, left, right, node, source, dest, time):
        """
        Appends the specified arrays to the end of the columns of this
        :class:`MigrationTable`. This allows many new rows to be added at once.

        All six parameters are mandatory, and must be numpy arrays of the
        same length (which is equal to the number of additional migrations
        to add to the table).

        :param left: The left coordinates (inclusive).
        :type left: numpy.ndarray, dtype=np.float64
        :param right: The right coordinates (exclusive).
        :type right: numpy.ndarray, dtype=np.float64
        :param node: The node IDs.
        :type node: numpy.ndarray, dtype=np.int32
        :param source: The source population IDs.
        :type source: numpy.ndarray, dtype=np.int32
        :param dest: The destination population IDs.
        :type dest: numpy.ndarray, dtype=np.int32
        :param time: The time of each migration.
        :type time: numpy.ndarray, dtype=np.int64
        """
        super(MigrationTable, self).append_columns(
            left, right, node, source, dest, time)


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
    A table defining the sites in a tree sequence. See the
    :ref:`definitions <sec_site_table_definition>` for details on the columns
    in this table and the
    :ref:`tree sequence requirements <sec_valid_tree_sequence_requirements>` section
    for the properties needed for a site table to be a part of a valid tree
    sequence.

    :warning: The numpy arrays returned by table attribute accesses are **copies**
        of the underlying data. In particular, this means that you cannot edit
        the values in the columns by updating the attribute arrays.

        **NOTE:** this behaviour may change in future.

    :ivar position: The array of site position coordinates.
    :vartype position: numpy.ndarray, dtype=np.float64
    :ivar ancestral_state: The flattened array of ancestral state strings.
        See :ref:`sec_tables_api_text_columns` for more details.
    :vartype ancestral_state: numpy.ndarray, dtype=np.int8
    :ivar ancestral_state_offset: The offsets of rows in the ancestral_state
        array. See :ref:`sec_tables_api_text_columns` for more details.
    :vartype ancestral_state_offset: numpy.ndarray, dtype=np.uint32
    :ivar metadata: The flattened array of binary metadata values. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata: numpy.ndarray, dtype=np.int8
    :ivar metadata_offset: The array of offsets into the metadata column. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata_offset: numpy.ndarray, dtype=np.uint32
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

    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        pos, ancestral_state, _, _, metadata = self.get_row(index)
        return SiteTableRow(pos, ancestral_state, metadata)

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

    def clear(self):
        """
        Deletes all rows in this table.
        """
        super(SiteTable, self).clear()

    def reset(self):
        # Deprecated alias for clear
        self.clear()

    def add_row(self, position, ancestral_state, metadata=None):
        """
        Adds a new row to this :class:`SiteTable` and returns the ID of the
        corresponding site.

        :param float position: The position of this site in genome coordinates.
        :param str ancestral_state: The state of this site at the root of the tree.
        :param bytes metadata: The binary-encoded metadata for the new node. If not
            specified or None, a zero-length byte string is stored.
        :return: The ID of the newly added site.
        :rtype: int
        """
        return super(SiteTable, self).add_row(position, ancestral_state, metadata)

    def set_columns(
            self, position, ancestral_state, ancestral_state_offset,
            metadata=None, metadata_offset=None):
        """
        Sets the values for each column in this :class:`.SiteTable` using the values
        in the specified arrays. Overwrites any data currently stored in the table.

        The ``position``, ``ancestral_state`` and ``ancestral_state_offset``
        parameters are mandatory, and must be 1D numpy arrays. The length
        of the ``position`` array determines the number of rows in table.
        The ``ancestral_state`` and ``ancestral_state_offset`` parameters must
        be supplied together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_text_columns` for more information). The
        ``metadata`` and ``metadata_offset`` parameters must be supplied
        together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_binary_columns` for more information).

        :param position: The position of each site in genome coordinates.
        :type position: numpy.ndarray, dtype=np.float64
        :param ancestral_state: The flattened ancestral_state array. Required.
        :type ancestral_state: numpy.ndarray, dtype=np.int8
        :param ancestral_state_offset: The offsets into the ``ancestral_state`` array.
        :type ancestral_state_offset: numpy.ndarray, dtype=np.uint32.
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        super(SiteTable, self).set_columns(
            position, ancestral_state=ancestral_state,
            ancestral_state_offset=ancestral_state_offset,
            metadata=metadata, metadata_offset=metadata_offset)

    def append_columns(
            self, position, ancestral_state, ancestral_state_offset,
            metadata=None, metadata_offset=None):
        """
        Appends the specified arrays to the end of the columns of this
        :class:`SiteTable`. This allows many new rows to be added at once.

        The ``position``, ``ancestral_state`` and ``ancestral_state_offset``
        parameters are mandatory, and must be 1D numpy arrays. The length
        of the ``position`` array determines the number of additional rows
        to add the table.
        The ``ancestral_state`` and ``ancestral_state_offset`` parameters must
        be supplied together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_text_columns` for more information). The
        ``metadata`` and ``metadata_offset`` parameters must be supplied
        together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_binary_columns` for more information).

        :param position: The position of each site in genome coordinates.
        :type position: numpy.ndarray, dtype=np.float64
        :param ancestral_state: The flattened ancestral_state array. Required.
        :type ancestral_state: numpy.ndarray, dtype=np.int8
        :param ancestral_state_offset: The offsets into the ``ancestral_state`` array.
        :type ancestral_state_offset: numpy.ndarray, dtype=np.uint32.
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        super(SiteTable, self).append_columns(
            position, ancestral_state=ancestral_state,
            ancestral_state_offset=ancestral_state_offset,
            metadata=metadata, metadata_offset=metadata_offset)


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
    A table defining the mutations in a tree sequence. See the
    :ref:`definitions <sec_mutation_table_definition>` for details on the columns
    in this table and the
    :ref:`tree sequence requirements <sec_valid_tree_sequence_requirements>` section
    for the properties needed for a mutation table to be a part of a valid tree
    sequence.

    :warning: The numpy arrays returned by table attribute accesses are **copies**
        of the underlying data. In particular, this means that you cannot edit
        the values in the columns by updating the attribute arrays.

        **NOTE:** this behaviour may change in future.

    :ivar site: The array of site IDs.
    :vartype site: numpy.ndarray, dtype=np.int32
    :ivar node: The array of node IDs.
    :vartype node: numpy.ndarray, dtype=np.int32
    :ivar derived_state: The flattened array of derived state strings.
        See :ref:`sec_tables_api_text_columns` for more details.
    :vartype derived_state: numpy.ndarray, dtype=np.int8
    :ivar derived_state_offset: The offsets of rows in the derived_state
        array. See :ref:`sec_tables_api_text_columns` for more details.
    :vartype derived_state_offset: numpy.ndarray, dtype=np.uint32
    :ivar parent: The array of parent mutation IDs.
    :vartype parent: numpy.ndarray, dtype=np.int32
    :ivar metadata: The flattened array of binary metadata values. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata: numpy.ndarray, dtype=np.int8
    :ivar metadata_offset: The array of offsets into the metadata column. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata_offset: numpy.ndarray, dtype=np.uint32
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

    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        site, node, derived_state, parent, metadata = self.get_row(index)
        return MutationTableRow(site, node, derived_state, parent, metadata)

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

    def clear(self):
        """
        Deletes all rows in this table.
        """
        super(MutationTable, self).clear()

    def reset(self):
        # Deprecated alias for clear
        self.clear()

    def add_row(self, site, node, derived_state, parent=-1, metadata=None):
        """
        Adds a new row to this :class:`MutationTable` and returns the ID of the
        corresponding mutation.

        :param int site: The ID of the site that this mutation occurs at.
        :param int node: The ID of the first node inheriting this mutation.
        :param str derived_state: The state of the site at this mutation's node.
        :param int parent: The ID of the parent mutation. If not specified,
            defaults to the :attr:`NULL_MUTATION`.
        :param bytes metadata: The binary-encoded metadata for the new node. If not
            specified or None, a zero-length byte string is stored.
        :return: The ID of the newly added mutation.
        :rtype: int
        """
        return super(MutationTable, self).add_row(
                site, node, derived_state, parent, metadata)

    def set_columns(
            self, site, node, derived_state, derived_state_offset,
            parent=None, metadata=None, metadata_offset=None):
        """
        Sets the values for each column in this :class:`.MutationTable` using the values
        in the specified arrays. Overwrites any data currently stored in the table.

        The ``site``, ``node``, ``derived_state`` and ``derived_state_offset``
        parameters are mandatory, and must be 1D numpy arrays. The
        ``site`` and ``node`` (also ``parent``, if supplied) arrays
        must be of equal length, and determine the number of rows in the table.
        The ``derived_state`` and ``derived_state_offset`` parameters must
        be supplied together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_text_columns` for more information). The
        ``metadata`` and ``metadata_offset`` parameters must be supplied
        together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_binary_columns` for more information).

        :param site: The ID of the site each mutation occurs at.
        :type site: numpy.ndarray, dtype=np.int32
        :param node: The ID of the node each mutation is associated with.
        :type node: numpy.ndarray, dtype=np.int32
        :param derived_state: The flattened derived_state array. Required.
        :type derived_state: numpy.ndarray, dtype=np.int8
        :param derived_state_offset: The offsets into the ``derived_state`` array.
        :type derived_state_offset: numpy.ndarray, dtype=np.uint32.
        :param parent: The ID of the parent mutation for each mutation.
        :type parent: numpy.ndarray, dtype=np.int32
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        super(MutationTable, self).set_columns(
            site=site, node=node, parent=parent,
            derived_state=derived_state, derived_state_offset=derived_state_offset,
            metadata=metadata, metadata_offset=metadata_offset)

    def append_columns(
            self, site, node, derived_state, derived_state_offset,
            parent=None, metadata=None, metadata_offset=None):
        """
        Appends the specified arrays to the end of the columns of this
        :class:`MutationTable`. This allows many new rows to be added at once.

        The ``site``, ``node``, ``derived_state`` and ``derived_state_offset``
        parameters are mandatory, and must be 1D numpy arrays. The
        ``site`` and ``node`` (also ``parent``, if supplied) arrays
        must be of equal length, and determine the number of additional
        rows to add to the table.
        The ``derived_state`` and ``derived_state_offset`` parameters must
        be supplied together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_text_columns` for more information). The
        ``metadata`` and ``metadata_offset`` parameters must be supplied
        together, and meet the requirements for
        :ref:`sec_encoding_ragged_columns` (see
        :ref:`sec_tables_api_binary_columns` for more information).

        :param site: The ID of the site each mutation occurs at.
        :type site: numpy.ndarray, dtype=np.int32
        :param node: The ID of the node each mutation is associated with.
        :type node: numpy.ndarray, dtype=np.int32
        :param derived_state: The flattened derived_state array. Required.
        :type derived_state: numpy.ndarray, dtype=np.int8
        :param derived_state_offset: The offsets into the ``derived_state`` array.
        :type derived_state_offset: numpy.ndarray, dtype=np.uint32.
        :param parent: The ID of the parent mutation for each mutation.
        :type parent: numpy.ndarray, dtype=np.int32
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        super(MutationTable, self).append_columns(
            site=site, node=node, parent=parent,
            derived_state=derived_state, derived_state_offset=derived_state_offset,
            metadata=metadata, metadata_offset=metadata_offset)


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
    .. todo::
        This class is provisional, and the API may change in the future.
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

    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        return ProvenanceTableRow(*self.get_row(index))

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

    def clear(self):
        """
        Deletes all rows in this table.
        """
        super(ProvenanceTable, self).clear()

    def reset(self):
        # Deprecated alias for clear
        self.clear()


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


def sort_tables(
        nodes, edges, migrations=None, sites=None, mutations=None,
        provenances=None, edge_start=0):
    """
    Sorts the given tables **in place**, ensuring that all tree
    sequence ordering requirements are met. See
    the :ref:`sec_valid_tree_sequence_requirements` section for details on these
    requirements.

    If the ``edge_start`` parameter is provided, this specifies the index
    in the edge table where sorting should start. Only rows with index
    greater than or equal to ``edge_start`` are sorted; rows before this index
    are not affected. This parameter is provided to allow for efficient sorting
    when the user knows that the edges up to a given index are already sorted.

    The input node table is not affected by this function.

    Edges are sorted as follows:

    - time of parent, then
    - parent node ID, then
    - child node ID, then
    - left endpoint.

    Note that this sorting order exceeds the
    :ref:`edge sorting requirements <sec_edge_requirements>` for a valid
    tree sequence. For a valid tree sequence, we require that all edges for a
    given parent ID are adjacent, but we do not require that they be listed in
    sorted order.

    Sites are sorted by position.

    Mutations are sorted by site ID.

    Migrations and provenances are not currently affected by this function.
    However, this behaviour is likely to change in the future.

    :param NodeTable nodes: The tree sequence nodes (required).
    :param EdgeTable edges: The tree sequence edges (required).
    :param MigrationTable migrations: The tree sequence migrations (optional).
    :param SiteTable sites: The tree sequence sites (optional, but required if
         ``mutations`` is provided)
    :param MutationTable mutations: The tree sequence mutations (optional, but
         required if ``sites`` is provided).
    :param ProvenanceTable provenances: Ignored. This argument is provided to
        support calling the function like ``sort_tables(**tables.asdict())``.
    :param int edge_start: The index in the edge table where sorting starts
        (default=0; must be <= len(edges)).
    """
    # TODO update the low-level module to accept None and remove this
    kwargs = {"nodes": nodes, "edges": edges, "edge_start": edge_start}
    if migrations is not None:
        kwargs["migrations"] = migrations
    if sites is not None:
        kwargs["sites"] = sites
    if mutations is not None:
        kwargs["mutations"] = mutations
    return _msprime.sort_tables(**kwargs)


def simplify_tables(
        samples, nodes, edges, migrations=None, sites=None, mutations=None,
        sequence_length=0, filter_zero_mutation_sites=True):
    """
    Simplifies the tables, **in place**, to retain only the information necessary
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
    :param bool filter_zero_mutation_sites: Whether to remove sites that have no
        mutations from the output (default: True).
    :param float sequence_length: The length of the sequence.
    :return: A numpy array mapping node IDs in the input tables to their
        corresponding node IDs in the output tables.
    :rtype: numpy array (dtype=np.int32).
    """
    # TODO update the low-level module to accept None and remove this
    kwargs = {
        "samples": samples, "nodes": nodes, "edges": edges,
        "sequence_length": sequence_length,
        "filter_zero_mutation_sites": filter_zero_mutation_sites}
    if migrations is not None:
        kwargs["migrations"] = migrations
    if sites is not None:
        kwargs["sites"] = sites
    if mutations is not None:
        kwargs["mutations"] = mutations
    return _msprime.simplify_tables(**kwargs)


def pack_bytes(data):
    """
    Packs the specified list of bytes into a flattened numpy array of 8 bit integers
    and corresponding offsets. See :ref:`sec_encoding_ragged_columns` for details
    of this encoding.

    :param list[bytes] data: The list of bytes values to encode.
    :return: The tuple (packed, offset) of numpy arrays representing the flattened
        input data and offsets.
    :rtype: numpy.array (dtype=np.int8), numpy.array (dtype=np.uint32).
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
    data and corresponding offsets. See :ref:`sec_encoding_ragged_columns` for details
    of this encoding.

    :param numpy.ndarray packed: The flattened array of byte values.
    :param numpy.ndarray offset: The array of offsets into the ``packed`` array.
    :return: The list of bytes values unpacked from the parameter arrays.
    :rtype: list[bytes]
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
    and corresponding offsets using the specified text encoding.
    See :ref:`sec_encoding_ragged_columns` for details of this encoding of
    columns of variable length data.

    :param list[str] data: The list of strings to encode.
    :param str encoding: The text encoding to use when converting string data
        to bytes. See the :mod:`codecs` module for information on available
        string encodings.
    :return: The tuple (packed, offset) of numpy arrays representing the flattened
        input data and offsets.
    :rtype: numpy.array (dtype=np.int8), numpy.array (dtype=np.uint32).
    """
    return pack_bytes([bytearray(s.encode(encoding)) for s in strings])


def unpack_strings(packed, offset, encoding="utf8"):
    """
    Unpacks a list of strings from the specified numpy arrays of packed byte
    data and corresponding offsets using the specified text encoding.
    See :ref:`sec_encoding_ragged_columns` for details of this encoding of
    columns of variable length data.

    :param numpy.ndarray packed: The flattened array of byte values.
    :param numpy.ndarray offset: The array of offsets into the ``packed`` array.
    :param str encoding: The text encoding to use when converting string data
        to bytes. See the :mod:`codecs` module for information on available
        string encodings.
    :return: The list of strings unpacked from the parameter arrays.
    :rtype: list[str]
    """
    return [b.decode(encoding) for b in unpack_bytes(packed, offset)]
