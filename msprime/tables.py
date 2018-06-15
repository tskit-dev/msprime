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

import numpy as np
from six.moves import copyreg

import _msprime
# This circular import is ugly but it seems hard to avoid it since table collection
# and tree sequence depend on each other. Unless they're in the same module they
# need to import each other. In Py3 at least we can import the modules but we
# can't do this in Py3.
import msprime


IndividualTableRow = collections.namedtuple(
    "IndividualTableRow",
    ["flags", "location", "metadata"])


NodeTableRow = collections.namedtuple(
    "NodeTableRow",
    ["flags", "time", "population", "individual", "metadata"])


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


PopulationTableRow = collections.namedtuple(
    "PopulationTableRow",
    ["metadata"])


ProvenanceTableRow = collections.namedtuple(
    "ProvenanceTableRow",
    ["timestamp", "record"])


# TODO We could abstract quite a lot more functionality up into this baseclass
# if each class kept a list of its columns. Then it would be pretty simple to
# define generic implementation of copy, etc.


class BaseTable(object):
    """
    Superclass of high-level tables. Not intended for direct instantiation.
    """
    def __init__(self, ll_table, row_class):
        self.ll_table = ll_table
        self.row_class = row_class

    @property
    def num_rows(self):
        return self.ll_table.num_rows

    @property
    def max_rows(self):
        return self.ll_table.max_rows

    @property
    def max_rows_increment(self):
        return self.ll_table.max_rows_increment

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = bool(self.ll_table.equals(other.ll_table))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        return self.row_class(*self.ll_table.get_row(index))

    def clear(self):
        """
        Deletes all rows in this table.
        """
        self.ll_table.clear()

    def reset(self):
        # Deprecated alias for clear
        self.clear()

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(**state)


class IndividualTable(BaseTable):
    """
    A table defining the individuals in a tree sequence. See the
    :ref:`definitions <sec_individual_table_definition>` for details on the columns
    in this table and the
    :ref:`tree sequence requirements <sec_valid_tree_sequence_requirements>` section
    for the properties needed for an individual table to be a part of a valid tree
    sequence.

    :warning: The numpy arrays returned by table attribute accesses are **copies**
        of the underlying data. In particular, this means that you cannot edit
        the values in the columns by updating the attribute arrays.

        **NOTE:** this behaviour may change in future.

    :ivar flags: The array of flags values.
    :vartype flags: numpy.ndarray, dtype=np.uint32
    :ivar location: The flattened array of floating point location values. See
        :ref:`sec_encoding_ragged_columns` for more details.
    :vartype location: numpy.ndarray, dtype=np.float64
    :ivar location_offset: The array of offsets into the location column. See
        :ref:`sec_encoding_ragged_columns` for more details.
    :vartype location_offset: numpy.ndarray, dtype=np.uint32
    :ivar metadata: The flattened array of binary metadata values. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata: numpy.ndarray, dtype=np.int8
    :ivar metadata_offset: The array of offsets into the metadata column. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata_offset: numpy.ndarray, dtype=np.uint32
    """
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.IndividualTable(max_rows_increment=max_rows_increment)
        super(IndividualTable, self).__init__(ll_table, IndividualTableRow)

    @property
    def flags(self):
        return self.ll_table.flags

    @property
    def location(self):
        return self.ll_table.location

    @property
    def location_offset(self):
        return self.ll_table.location_offset

    @property
    def metadata(self):
        return self.ll_table.metadata

    @property
    def metadata_offset(self):
        return self.ll_table.metadata_offset

    def __str__(self):
        flags = self.flags
        location = self.location
        location_offset = self.location_offset
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tflags\tlocation\tmetadata\n"
        for j in range(self.num_rows):
            md = base64.b64encode(metadata[j]).decode('utf8')
            location_str = ",".join(map(
                str, location[location_offset[j]: location_offset[j + 1]]))
            ret += "{}\t{}\t{}\t{}\n".format(j, flags[j], location_str, md)
        return ret[:-1]

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = IndividualTable()
        copy.set_columns(
            flags=self.flags,
            location=self.location, location_offset=self.location_offset,
            metadata=self.metadata, metadata_offset=self.metadata_offset)
        return copy

    def add_row(self, flags=0, location=None, metadata=None):
        """
        Adds a new row to this :class:`IndividualTable` and returns the ID of the
        corresponding individual.

        :param int flags: The bitwise flags for the new node.
        :param array-like location: A list of numeric values or one-dimensional numpy
            array describing the location of this individual. If not specified
            or None, a zero-dimensional location is stored.
        :param bytes metadata: The binary-encoded metadata for the new node. If not
            specified or None, a zero-length byte string is stored.
        :return: The ID of the newly added node.
        :rtype: int
        """
        return self.ll_table.add_row(flags=flags, location=location, metadata=metadata)

    def set_columns(
            self, flags, location=None, location_offset=None,
            metadata=None, metadata_offset=None):
        """
        Sets the values for each column in this :class:`.IndividualTable` using the
        values in the specified arrays. Overwrites any data currently stored in
        the table.

        The ``flags`` array is mandatory and defines the number of individuals
        the table will contain.
        The ``location`` and ``location_offset`` parameters must be supplied
        together, and meet the requirements for :ref:`sec_encoding_ragged_columns`.
        The ``metadata`` and ``metadata_offset`` parameters must be supplied
        together, and meet the requirements for :ref:`sec_encoding_ragged_columns`.
        See :ref:`sec_tables_api_binary_columns` for more information.

        :param flags: The bitwise flags for each node. Required.
        :type flags: numpy.ndarray, dtype=np.uint32
        :param location: The flattened location array. Must be specified along
            with ``location_offset``. If not specified or None, an empty location
            value is stored for each node.
        :type location: numpy.ndarray, dtype=np.float64
        :param location_offset: The offsets into the ``location`` array.
        :type location_offset: numpy.ndarray, dtype=np.uint32.
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        self.ll_table.set_columns(
            flags, location=location, location_offset=location_offset,
            metadata=metadata, metadata_offset=metadata_offset)

    def append_columns(
            self, flags, location=None, location_offset=None, metadata=None,
            metadata_offset=None):
        """
        Appends the specified arrays to the end of the columns in this
        :class:`IndividualTable`. This allows many new rows to be added at once.

        The ``flags`` array is mandatory and defines the number of
        extra individuals to add to the table.
        The ``location`` and ``location_offset`` parameters must be supplied
        together, and meet the requirements for :ref:`sec_encoding_ragged_columns`.
        The ``metadata`` and ``metadata_offset`` parameters must be supplied
        together, and meet the requirements for :ref:`sec_encoding_ragged_columns`.
        See :ref:`sec_tables_api_binary_columns` for more information.

        :param flags: The bitwise flags for each node. Required.
        :type flags: numpy.ndarray, dtype=np.uint32
        :param location: The flattened location array. Must be specified along
            with ``location_offset``. If not specified or None, an empty location
            value is stored for each node.
        :type location: numpy.ndarray, dtype=np.float64
        :param location_offset: The offsets into the ``location`` array.
        :type location_offset: numpy.ndarray, dtype=np.uint32.
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        self.ll_table.append_columns(
            flags, location=location, location_offset=location_offset,
            metadata=metadata, metadata_offset=metadata_offset)


# Pickle support. See copyreg registration for this function below.
def _pickle_individual_table(table):
    state = {
        "flags": table.flags,
        "location": table.location,
        "location_offset": table.location_offset,
        "metadata": table.metadata,
        "metadata_offset": table.metadata_offset,
    }
    return IndividualTable, tuple(), state


class NodeTable(BaseTable):
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
    :vartype individual: numpy.ndarray, dtype=np.int32
    :ivar metadata: The flattened array of binary metadata values. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata: numpy.ndarray, dtype=np.int8
    :ivar metadata_offset: The array of offsets into the metadata column. See
        :ref:`sec_tables_api_binary_columns` for more details.
    :vartype metadata_offset: numpy.ndarray, dtype=np.uint32
    """
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.NodeTable(max_rows_increment=max_rows_increment)
        super(NodeTable, self).__init__(ll_table, NodeTableRow)

    @property
    def time(self):
        return self.ll_table.time

    @property
    def flags(self):
        return self.ll_table.flags

    @property
    def population(self):
        return self.ll_table.population

    @property
    def individual(self):
        return self.ll_table.individual

    # EXPERIMENTAL interface for setting a single column. This is done
    # quite a bit in tests. Not part of the public API as yet, but we
    # probably will want to allow something like this in general.
    @individual.setter
    def individual(self, individual):
        self.set_columns(
            flags=self.flags, time=self.time, metadata=self.metadata,
            metadata_offset=self.metadata_offset, individual=individual)

    @property
    def metadata(self):
        return self.ll_table.metadata

    @property
    def metadata_offset(self):
        return self.ll_table.metadata_offset

    def __str__(self):
        time = self.time
        flags = self.flags
        population = self.population
        individual = self.individual
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tflags\tpopulation\tindividual\ttime\tmetadata\n"
        for j in range(self.num_rows):
            md = base64.b64encode(metadata[j]).decode('utf8')
            ret += "{}\t{}\t{}\t{}\t{:.14f}\t{}\n".format(
                j, flags[j], population[j], individual[j], time[j], md)
        return ret[:-1]

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = NodeTable()
        copy.set_columns(
            flags=self.flags, time=self.time, population=self.population,
            individual=self.individual, metadata=self.metadata,
            metadata_offset=self.metadata_offset)
        return copy

    def add_row(self, flags=0, time=0, population=-1, individual=-1, metadata=None):
        """
        Adds a new row to this :class:`NodeTable` and returns the ID of the
        corresponding node.

        :param int flags: The bitwise flags for the new node.
        :param float time: The birth time for the new node.
        :param int population: The ID of the population in which the new node was born.
            Defaults to the :const:`.NULL_POPULATION`.
        :param int individual: The ID of the individual in which the new node was born.
            Defaults to the :const:`.NULL_INDIVIDUAL`.
        :param bytes metadata: The binary-encoded metadata for the new node. If not
            specified or None, a zero-length byte string is stored.
        :return: The ID of the newly added node.
        :rtype: int
        """
        return self.ll_table.add_row(flags, time, population, individual, metadata)

    def set_columns(
            self, flags, time, population=None, individual=None, metadata=None,
            metadata_offset=None):
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
        :param individual: The individual values for each node. If not specified
            or None, the :const:`.NULL_INDIVIDUAL` value is stored for each node.
        :type individual: numpy.ndarray, dtype=np.int32
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        self.ll_table.set_columns(
            flags, time, population=population, individual=individual, metadata=metadata,
            metadata_offset=metadata_offset)

    def append_columns(
            self, flags, time, population=None, individual=None, metadata=None,
            metadata_offset=None):
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
        :param individual: The individual values for each node. If not specified
            or None, the :const:`.NULL_INDIVIDUAL` value is stored for each node.
        :type individual: numpy.ndarray, dtype=np.int32
        :param metadata: The flattened metadata array. Must be specified along
            with ``metadata_offset``. If not specified or None, an empty metadata
            value is stored for each node.
        :type metadata: numpy.ndarray, dtype=np.int8
        :param metadata_offset: The offsets into the ``metadata`` array.
        :type metadata_offset: numpy.ndarray, dtype=np.uint32.
        """
        self.ll_table.append_columns(
            flags, time, population=population, individual=individual,
            metadata=metadata, metadata_offset=metadata_offset)


# Pickle support. See copyreg registration for this function below.
def _pickle_node_table(table):
    state = {
        "time": table.time,
        "flags": table.flags,
        "population": table.population,
        "individual": table.individual,
        "metadata": table.metadata,
        "metadata_offset": table.metadata_offset,
    }
    return NodeTable, tuple(), state


class EdgeTable(BaseTable):
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
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.EdgeTable(max_rows_increment=max_rows_increment)
        super(EdgeTable, self).__init__(ll_table, EdgeTableRow)

    @property
    def left(self):
        return self.ll_table.left

    @property
    def right(self):
        return self.ll_table.right

    @property
    def parent(self):
        return self.ll_table.parent

    @property
    def child(self):
        return self.ll_table.child

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

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = EdgeTable()
        copy.set_columns(
            left=self.left, right=self.right, parent=self.parent, child=self.child)
        return copy

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
        return self.ll_table.add_row(left, right, parent, child)

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
        self.ll_table.set_columns(left, right, parent, child)

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
        self.ll_table.append_columns(left, right, parent, child)


# Pickle support. See copyreg registration for this function below.
def _edge_table_pickle(table):
    state = {
        "left": table.left,
        "right": table.right,
        "parent": table.parent,
        "child": table.child,
    }
    return EdgeTable, tuple(), state


class MigrationTable(BaseTable):
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
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.MigrationTable(max_rows_increment=max_rows_increment)
        super(MigrationTable, self).__init__(ll_table, MigrationTableRow)

    @property
    def left(self):
        return self.ll_table.left

    @property
    def right(self):
        return self.ll_table.right

    @property
    def node(self):
        return self.ll_table.node

    @property
    def source(self):
        return self.ll_table.source

    @property
    def dest(self):
        return self.ll_table.dest

    @property
    def time(self):
        return self.ll_table.time

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

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = MigrationTable()
        copy.set_columns(
            left=self.left, right=self.right, node=self.node, source=self.source,
            dest=self.dest, time=self.time)
        return copy

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
        return self.ll_table.add_row(left, right, node, source, dest, time)

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
        self.ll_table.set_columns(left, right, node, source, dest, time)

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
        self.ll_table.append_columns(left, right, node, source, dest, time)


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


class SiteTable(BaseTable):
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
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.SiteTable(max_rows_increment=max_rows_increment)
        super(SiteTable, self).__init__(ll_table, SiteTableRow)

    @property
    def position(self):
        return self.ll_table.position

    @property
    def ancestral_state(self):
        return self.ll_table.ancestral_state

    @property
    def ancestral_state_offset(self):
        return self.ll_table.ancestral_state_offset

    @property
    def metadata(self):
        return self.ll_table.metadata

    @property
    def metadata_offset(self):
        return self.ll_table.metadata_offset

    def __str__(self):
        position = self.position
        ancestral_state = unpack_strings(
            self.ancestral_state, self.ancestral_state_offset)
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tposition\tancestral_state\tmetadata\n"
        for j in range(self.num_rows):
            md = base64.b64encode(metadata[j]).decode('utf8')
            ret += "{}\t{:.8f}\t{}\t{}\n".format(
                j, position[j], ancestral_state[j], md)
        return ret[:-1]

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
        return self.ll_table.add_row(position, ancestral_state, metadata)

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
        self.ll_table.set_columns(
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
        self.ll_table.append_columns(
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


class MutationTable(BaseTable):
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
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.MutationTable(max_rows_increment=max_rows_increment)
        super(MutationTable, self).__init__(ll_table, MutationTableRow)

    @property
    def site(self):
        return self.ll_table.site

    @property
    def node(self):
        return self.ll_table.node

    @property
    def parent(self):
        return self.ll_table.parent

    @property
    def derived_state(self):
        return self.ll_table.derived_state

    @property
    def derived_state_offset(self):
        return self.ll_table.derived_state_offset

    @property
    def metadata(self):
        return self.ll_table.metadata

    @property
    def metadata_offset(self):
        return self.ll_table.metadata_offset

    def __str__(self):
        site = self.site
        node = self.node
        parent = self.parent
        derived_state = unpack_strings(self.derived_state, self.derived_state_offset)
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tsite\tnode\tderived_state\tparent\tmetadata\n"
        for j in range(self.num_rows):
            md = base64.b64encode(metadata[j]).decode('utf8')
            ret += "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                j, site[j], node[j], derived_state[j], parent[j], md)
        return ret[:-1]

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
        return self.ll_table.add_row(
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
        self.ll_table.set_columns(
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
        self.ll_table.append_columns(
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


class PopulationTable(BaseTable):
    # FIXME
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.PopulationTable(max_rows_increment=max_rows_increment)
        super(PopulationTable, self).__init__(ll_table, PopulationTableRow)

    @property
    def metadata(self):
        return self.ll_table.metadata

    @property
    def metadata_offset(self):
        return self.ll_table.metadata_offset

    def add_row(self, metadata=None):
        """
        .. todo:: document.
        """
        return self.ll_table.add_row(metadata=metadata)

    def __str__(self):
        metadata = unpack_strings(self.metadata, self.metadata_offset)
        ret = "id\tmetadata\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\n".format(j, metadata[j])
        return ret[:-1]

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = PopulationTable()
        copy.set_columns(
            metadata=self.metadata,
            metadata_offset=self.metadata_offset)
        return copy

    def set_columns(self, metadata=None, metadata_offset=None):
        self.ll_table.set_columns(metadata=metadata, metadata_offset=metadata_offset)

    def append_columns(self, metadata=None, metadata_offset=None):
        self.ll_table.append_columns(metadata=metadata, metadata_offset=metadata_offset)


# Pickle support. See copyreg registration for this function below.
def _population_table_pickle(table):
    state = {
        "metadata": table.metadata,
        "metadata_offset": table.metadata_offset,
    }
    return PopulationTable, tuple(), state


class ProvenanceTable(BaseTable):
    """
    .. todo::
        This class is provisional, and the API may change in the future.
    """
    # FIXME
    def __init__(self, max_rows_increment=0, ll_table=None):
        if ll_table is None:
            ll_table = _msprime.ProvenanceTable(max_rows_increment=max_rows_increment)
        super(ProvenanceTable, self).__init__(ll_table, ProvenanceTableRow)

    @property
    def record(self):
        return self.ll_table.record

    @property
    def record_offset(self):
        return self.ll_table.record_offset

    @property
    def timestamp(self):
        return self.ll_table.timestamp

    @property
    def timestamp_offset(self):
        return self.ll_table.timestamp_offset

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
        return self.ll_table.add_row(record=record, timestamp=timestamp)

    def set_columns(
            self, timestamp=None, timestamp_offset=None,
            record=None, record_offset=None):
        self.ll_table.set_columns(
            timestamp=timestamp, timestamp_offset=timestamp_offset,
            record=record, record_offset=record_offset)

    def append_columns(
            self, timestamp=None, timestamp_offset=None,
            record=None, record_offset=None):
        self.ll_table.append_columns(
            timestamp=timestamp, timestamp_offset=timestamp_offset,
            record=record, record_offset=record_offset)

    def __str__(self):
        timestamp = unpack_strings(self.timestamp, self.timestamp_offset)
        record = unpack_strings(self.record, self.record_offset)
        ret = "id\ttimestamp\trecord\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\t{}\n".format(j, timestamp[j], record[j])
        return ret[:-1]

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
copyreg.pickle(IndividualTable, _pickle_individual_table)
copyreg.pickle(NodeTable, _pickle_node_table)
copyreg.pickle(EdgeTable, _edge_table_pickle)
copyreg.pickle(MigrationTable, _migration_table_pickle)
copyreg.pickle(SiteTable, _site_table_pickle)
copyreg.pickle(MutationTable, _mutation_table_pickle)
copyreg.pickle(PopulationTable, _population_table_pickle)
copyreg.pickle(ProvenanceTable, _provenance_table_pickle)


class TableCollection(object):
    """
    A collection of mutable tables defining a tree sequence. See the
    :ref:`sec_data_model` section for definition on the various tables
    and how they together define a :class:`TreeSequence`. Arbitrary
    data can be stored in a TableCollection, but there are certain
    :ref:`requirements <sec_valid_tree_sequence_requirements>` that must be
    satisfied for these tables to be interpreted as a tree sequence.

    To obtain a :class:`TreeSequence` instance corresponding to the current
    state of a ``TableCollection``, please use the :meth:`.tree_sequence`
    method.

    :ivar individuals: The individual table.
    :vartype individuals: IndividualTable
    :ivar edges: The edge table.
    :vartype edges: EdgeTable
    :ivar migrations: The migration table.
    :vartype migrations: MigrationTable
    :ivar sites: The site table.
    :vartype sites: SiteTable
    :ivar mutations: The mutation table.
    :vartype mutations: MutationTable
    :ivar populations: The population table.
    :vartype populations: PopulationTable
    :ivar provenances: The provenance table.
    :vartype provenances: ProvenanceTable
    :ivar sequence_length: The sequence length defining the coordinate
        space.
    :vartype sequence_length: float
    """
    def __init__(self, sequence_length=0, ll_tables=None):
        if ll_tables is None:
            ll_tables = _msprime.TableCollection(
                individuals=_msprime.IndividualTable(),
                nodes=_msprime.NodeTable(),
                edges=_msprime.EdgeTable(),
                migrations=_msprime.MigrationTable(),
                sites=_msprime.SiteTable(),
                mutations=_msprime.MutationTable(),
                populations=_msprime.PopulationTable(),
                provenances=_msprime.ProvenanceTable(),
                sequence_length=sequence_length)
        self.ll_tables = ll_tables
        self.__individuals = IndividualTable(ll_table=self.ll_tables.individuals)
        self.__nodes = NodeTable(ll_table=self.ll_tables.nodes)
        self.__edges = EdgeTable(ll_table=self.ll_tables.edges)
        self.__migrations = MigrationTable(ll_table=self.ll_tables.migrations)
        self.__sites = SiteTable(ll_table=self.ll_tables.sites)
        self.__mutations = MutationTable(ll_table=self.ll_tables.mutations)
        self.__populations = PopulationTable(ll_table=self.ll_tables.populations)
        self.__provenances = ProvenanceTable(ll_table=self.ll_tables.provenances)

    @property
    def individuals(self):
        return self.__individuals

    @property
    def nodes(self):
        return self.__nodes

    @property
    def edges(self):
        return self.__edges

    @property
    def migrations(self):
        return self.__migrations

    @property
    def sites(self):
        return self.__sites

    @property
    def mutations(self):
        return self.__mutations

    @property
    def populations(self):
        return self.__populations

    @property
    def provenances(self):
        return self.__provenances

    @property
    def sequence_length(self):
        return self.ll_tables.sequence_length

    def asdict(self):
        """
        Returns this TableCollection as a dictionary mapping the keys "nodes",
        "edges", etc to their respective table objects.
        """
        return {
            # Leaving individuals out for now to until stuff that depends on it
            # is implemented.
            "individuals": self.individuals,
            "nodes": self.nodes,
            "edges": self.edges,
            "migrations": self.migrations,
            "sites": self.sites,
            "mutations": self.mutations,
            "populations": self.populations,
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
        s = self.__banner("Individuals")
        s += str(self.individuals) + "\n"
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

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = bool(self.ll_tables.equals(other.ll_tables))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def tree_sequence(self):
        """
        Returns a :class:`TreeSequence` instance with the structure defined by the
        tables in this :class:`TableCollection`. If the table collection is not
        in canonical form (i.e., does not meet sorting requirements) or cannot be
        interpreted as a tree sequence an exception is raised. The
        :meth:`.sort` method may be used to ensure that input sorting requirements
        are met.

        :return: A :class:`TreeSequence` instance reflecting the structures
            defined in this set of tables.
        :rtype: .TreeSequence
        """
        return msprime.TreeSequence.load_tables(self)

    def simplify(self, samples, filter_zero_mutation_sites=True):
        """
        Simplifies the tables in place to retain only the information necessary
        to reconstruct the tree sequence describing the given ``samples``.
        This will change the ID of the nodes, so that the individual
        ``samples[k]]`` will have ID ``k`` in the result. The resulting
        NodeTable will have only the first ``len(samples)`` individuals marked
        as samples. The mapping from node IDs in the current set of tables to
        their equivalent values in the simplified tables is returned as a numpy
        array. If an array ``a`` is returned by this function and ``u`` is the
        ID of a node in the input table, then ``a[u]`` is the ID of this node
        in the output table. For any node ``u`` that is not mapped into the
        output tables, this mapping will equal ``-1``.

        Tables operated on by this function must: be sorted (see
        :meth:`TableCollection.sort`)), have children be born strictly after their
        parents, and the intervals on which any individual is a child must be
        disjoint; but other than this the tables need not satisfy remaining
        requirements to specify a valid tree sequence (but the resulting tables
        will).

        :param list[int] samples: A list of Node IDs of individuals to retain
            as samples.
        :param bool filter_zero_mutation_sites: Whether to remove sites that have no
            mutations from the output (default: True).
        :return: A numpy array mapping node IDs in the input tables to their
            corresponding node IDs in the output tables.
        :rtype: numpy array (dtype=np.int32).
        """
        return self.ll_tables.simplify(samples, filter_zero_mutation_sites)

    def sort(self, edge_start=0):
        """
        Sorts the tables in place, ensuring that all tree sequence ordering
        requirements are met. See the
        :ref:`sec_valid_tree_sequence_requirements` section for details on
        these requirements.

        If the ``edge_start`` parameter is provided, this specifies the index
        in the edge table where sorting should start. Only rows with index
        greater than or equal to ``edge_start`` are sorted; rows before this index
        are not affected. This parameter is provided to allow for efficient sorting
        when the user knows that the edges up to a given index are already sorted.

        The individual, node, population and provenance tables are not affected
        by this method.

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

        :param int edge_start: The index in the edge table where sorting starts
            (default=0; must be <= len(edges)).
        """
        self.ll_tables.sort(edge_start)

    def compute_mutation_parents(self):
        # TODO document
        self.ll_tables.compute_mutation_parents()

    def deduplicate_sites(self):
        # TODO document
        self.ll_tables.deduplicate_sites()


#############################################
# Table functions.
#############################################


def sort_tables(
        nodes, edges, migrations=None, sites=None, mutations=None,
        provenances=None, individuals=None, populations=None, edge_start=0):
    """
    **This function is deprecated. Please use TableCollection.sort() instead**

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

    :param NodeTable nodes: The nodes of the tree sequence (required).
    :param EdgeTable edges: The edges of the tree sequence (required).
    :param MigrationTable migrations: The tree sequence's migrations (optional).
    :param SiteTable sites: The tree sequence's sites (optional, but required if
         ``mutations`` is provided)
    :param MutationTable mutations: The tree sequence's mutations (optional, but
         required if ``sites`` is provided).
    :param ProvenanceTable provenances: Ignored. This argument is provided to
        support calling the function like ``sort_tables(**tables.asdict())``.
    :param PopulationTable populations: Ignored. This argument is provided to
        support calling the function like ``sort_tables(**tables.asdict())``.
    :param IndividualTable individuals: Ignored. This argument is provided to
        support calling the function like ``sort_tables(**tables.asdict())``.
    :param int edge_start: The index in the edge table where sorting starts
        (default=0; must be <= len(edges)).
    """
    if migrations is None:
        migrations = MigrationTable()
    if sites is None:
        sites = SiteTable()
    if mutations is None:
        mutations = MutationTable()
    sequence_length = 0
    if len(edges) > 0:
        sequence_length = edges.right.max()
    try:
        ll_tables = _msprime.TableCollection(
            individuals=_msprime.IndividualTable(),
            nodes=nodes.ll_table,
            edges=edges.ll_table,
            migrations=migrations.ll_table,
            sites=sites.ll_table,
            mutations=mutations.ll_table,
            populations=_msprime.PopulationTable(),
            provenances=_msprime.ProvenanceTable(),
            sequence_length=sequence_length)
    except AttributeError as e:
        raise TypeError(str(e))
    return ll_tables.sort(edge_start)


def simplify_tables(
        samples, nodes, edges, migrations=None, sites=None, mutations=None,
        sequence_length=0, filter_zero_mutation_sites=True):
    """
    **This function is deprecated. Please use TableCollection.simplify() instead**

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
    if migrations is None:
        migrations = MigrationTable()
    if sites is None:
        sites = SiteTable()
    if mutations is None:
        mutations = MutationTable()
    if sequence_length == 0 and len(edges) > 0:
        sequence_length = edges.right.max()
    try:
        ll_tables = _msprime.TableCollection(
            individuals=_msprime.IndividualTable(),
            nodes=nodes.ll_table,
            edges=edges.ll_table,
            migrations=migrations.ll_table,
            sites=sites.ll_table,
            mutations=mutations.ll_table,
            populations=_msprime.PopulationTable(),
            provenances=_msprime.ProvenanceTable(),
            sequence_length=sequence_length)
    except AttributeError as e:
        raise TypeError(str(e))
    return ll_tables.simplify(
        samples, filter_zero_mutation_sites=filter_zero_mutation_sites)


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
