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

# import datetime

from six.moves import copyreg

import _msprime

# We need numpy for most API operations that involve tables. We should make
# numpy a hard requirement soon and get rid of these conditional imports.
try:
    import numpy as np
except ImportError:
    pass


class NodeTable(_msprime.NodeTable):
    """
    Class for tables describing all nodes in the tree sequence, of the form
        id     is_sample  population   time
        0      1          0            0.0
        1      1          1            0.0
        2      0          0            0.0
        3      1          0            0.5
        4      0          2            2.1
    Node IDs are *not* recorded; rather the `id` column shows the row index, so
    that the `k`-th row describes the node whose ID is `k`.  `is_sample`
    records whether the node is a sample (=1) or not (=0).  `population` is an
    integer population ID, and `time` is the time since that individual was
    born, as a float.

    Requirements:
        1. All birth times must be greater than or equal to zero.

    It is not required that the `time` column be ordered or that all samples
    must be at the top.
    """
    def __str__(self):
        time = self.time
        flags = self.flags
        population = self.population
        metadata = unpack_bytes(self.metadata, self.metadata_offset)
        ret = "id\tis_sample\tpopulation\ttime\tmetadata\n"
        for j in range(self.num_rows):
            # Not clear how we should deal with printing metadata out here.
            # Probably try to decode as utf8, but printout as raw bytes if
            # this fails??
            ret += "{}\t{}\t{}\t{:.14f}\t{}\n".format(
                j, flags[j], population[j], time[j], metadata[j])
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
    TODO Update docs to relfect EDGES.


    Class for tables describing all edges in a tree sequence, of the form
        left	right	parent	child
        0.0     0.4     3       0
        0.0     0.4     3       2
        0.4     1.0     3       0,
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
        ret = "id\tposition\tancestral_state\n"
        for j in range(self.num_rows):
            ret += "{}\t{:.8f}\t{}\n".format(j, position[j], ancestral_state[j])
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.position, other.position) and
                np.array_equal(self.ancestral_state, other.ancestral_state) and
                np.array_equal(
                    self.ancestral_state_offset, other.ancestral_state_offset))
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return self.num_rows

    # Unpickle support
    def __setstate__(self, state):
        self.set_columns(
            position=state["position"], ancestral_state=state["ancestral_state"],
            ancestral_state_offset=state["ancestral_state_offset"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = SiteTable()
        copy.set_columns(
            position=self.position, ancestral_state=self.ancestral_state,
            ancestral_state_offset=self.ancestral_state_offset)
        return copy


# Pickle support. See copyreg registration for this function below.
def _site_table_pickle(table):
    state = {
        "position": table.position,
        "ancestral_state": table.ancestral_state,
        "ancestral_state_offset": table.ancestral_state_offset,
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
        derived_state = unpack_strings(
            self.derived_state, self.derived_state_offset)
        ret = "id\tsite\tnode\tderived_state\tparent\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\t{}\t{}\t{}\n".format(
                j, site[j], node[j], derived_state[j], parent[j])
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
                    self.derived_state_offset, other.derived_state_offset))
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
            derived_state_offset=state["derived_state_offset"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = MutationTable()
        copy.set_columns(
            site=self.site, node=self.node, parent=self.parent,
            derived_state=self.derived_state,
            derived_state_offset=self.derived_state_offset)
        return copy


# Pickle support. See copyreg registration for this function below.
def _mutation_table_pickle(table):
    state = {
        "site": table.site,
        "node": table.node,
        "parent": table.parent,
        "derived_state": table.derived_state,
        "derived_state_offset": table.derived_state_offset,
    }
    return MutationTable, tuple(), state


class ProvenanceTable(_msprime.ProvenanceTable):
    """
    TODO Document
    """
    # def add_row(self, record, timestamp=None):
    #     if timestamp is None:
    #         timestamp = datetime.datetime.now().isoformat()
    #     print(timestamp)

    #     super(ProvenanceTable, self).add_row(record=record, timestamp=timestamp)

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
            self, nodes=None, edges=None, migrations=None, sites=None, mutations=None):
        self.nodes = nodes
        self.edges = edges
        self.migrations = migrations
        self.sites = sites
        self.mutations = mutations

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
            "mutations": self.mutations
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
        s += str(self.migrations)
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

    .. todo:: Update this documentation to describe the keyword arguments and
       combinations that are allowed.

    :param NodeTable nodes:
    :param EdgeTable edges:
    :param MigrationTable migrations:
    :param SiteTable sites:
    :param MutationTable mutations:
    :param int edge_start: The index in the edge table where sorting starts.
    """
    return _msprime.sort_tables(*args, **kwargs)


def simplify_tables(*args, **kwargs):
    """
    Simplifies the tables, in place, to retain only the information necessary
    to reconstruct the tree sequence describing the given ``samples``.  This
    will change the ID of the nodes, so that the individual ``samples[k]]``
    will have ID ``k`` in the result.  The resulting NodeTable will have only
    the first ``len(samples)`` individuals marked as samples.

    Tables operated on by this function must: be sorted (see ``sort_tables``),
    have children be born strictly after their parents, and the intervals on
    which any individual is a child must be disjoint; but other than this the
    tables need not satisfy remaining requirements to specify a valid tree
    sequence (but the resulting tables will).

    :param list samples: A list of Node IDs of individuals to retain as samples.
    :param NodeTable nodes: The NodeTable to be simplified.
    :param EdgeTable edges: The EdgeTable to be simplified.
    :param MigrationTable migrations: The MigrationTable to be simplified.
    :param SiteTable sites: The SiteTable to be simplified.
    :param MutationTable mutations: The MutationTable to be simplified.
    :param bool filter_invariant_sites: Whether to remove sites that have no
        mutations from the output (default: True).
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
