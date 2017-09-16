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
        id	flags	population	time
        0	1	0		0.0
        1	1	1		0.0
        2	0	0		0.0
        3	1	0		0.5
        4	0	2		2.1
    Node IDs are *not* recorded; rather the `id` column shows the row index, so
    that the `k`-th row describes the node whose ID is `k`.  `flags` currently
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
        ret = "id\tflags\tpopulation\ttime\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\t{}\t\t{:.14f}\n".format(j, flags[j], population[j], time[j])
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.flags, other.flags) and
                np.array_equal(self.population, other.population) and
                np.array_equal(self.time, other.time) and
                np.array_equal(self.name, other.name) and
                np.array_equal(self.name_length, other.name_length))
        return ret

    def __len__(self):
        return self.num_rows

    # Pickle support. See copyreg registration for this function below.
    def _pickle(table):
        state = {
            "time": table.time,
            "flags": table.flags,
            "population": table.population,
            "name": table.name,
            "name_length": table.name_length,
        }
        return NodeTable, tuple(), state

    # Unpickle support
    def __setstate__(self, state):
        self.__init__()
        self.set_columns(
            time=state["time"], flags=state["flags"], population=state["population"],
            name=state["name"], name_length=state["name_length"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = NodeTable()
        copy.set_columns(
            flags=self.flags, time=self.time, population=self.population,
            name=self.name, name_length=self.name_length)
        return copy


class EdgesetTable(_msprime.EdgesetTable):
    """
    Class for tables describing all edgesets in a tree sequence, of the form
        left	right	parent	children
        0.0     0.4     3       0,2
        0.4     1.0     3       0,1,2
        0.0     0.4     4       1,3
    These describe the half-open genomic interval affected: `[left, right)`,
    the `parent` and the `children` on that interval.

    Requirements: to describe a valid tree sequence, a `EdgesetTable` (and
    corresponding `NodeTable`, to provide birth times) must satisfy:
        1. each list of children must be in sorted order,
        2. any two edgesets that share a child must be nonoverlapping, and
        3. the birth times of the `parent` in an edgeset must be strictly
            greater than the birth times of the `children` in that edgeset.
    Furthermore, for algorithmic requirements
        4. the smallest `left` coordinate must be 0.0,
        5. the the table must be sorted by birth time of the `parent`, and
        6. any two edgesets corresponding to the same `parent` must be nonoverlapping.
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
        children = self.children
        children_length = self.children_length
        ret = "id\tleft\t\tright\t\tparent\tchildren\n"
        offset = 0
        for j in range(self.num_rows):
            row_children = children[offset: offset + children_length[j]]
            offset += children_length[j]
            ret += "{}\t{:.8f}\t{:.8f}\t{}\t{}\n".format(
                j, left[j], right[j], parent[j], ",".join(map(str, row_children)))
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.left, other.left) and
                np.array_equal(self.right, other.right) and
                np.array_equal(self.parent, other.parent) and
                np.array_equal(self.children, other.children) and
                np.array_equal(self.children_length, other.children_length))
        return ret

    def __len__(self):
        return self.num_rows

    # Pickle support. See copyreg registration for this function below.
    def _pickle(self):
        state = {
            "left": self.left,
            "right": self.right,
            "parent": self.parent,
            "children": self.children,
            "children_length": self.children_length,
        }
        return EdgesetTable, tuple(), state

    # Unpickle support
    def __setstate__(self, state):
        self.__init__()
        self.set_columns(
            left=state["left"], right=state["right"], parent=state["parent"],
            children=state["children"], children_length=state["children_length"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = EdgesetTable()
        copy.set_columns(
            left=self.left, right=self.right, parent=self.parent,
            children=self.children, children_length=self.children_length)
        return copy


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

    def __len__(self):
        return self.num_rows

    # Pickle support. See copyreg registration for this function below.
    def _pickle(self):
        state = {
            "left": self.left,
            "right": self.right,
            "node": self.node,
            "source": self.source,
            "dest": self.dest,
            "time": self.time,
        }
        return MigrationTable, tuple(), state

    # Unpickle support
    def __setstate__(self, state):
        self.__init__()
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
            self.ancestral_state, self.ancestral_state_length)
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
                    self.ancestral_state_length, other.ancestral_state_length))
        return ret

    def __len__(self):
        return self.num_rows

    # Pickle support. See copyreg registration for this function below.
    def _pickle(self):
        state = {
            "position": self.position,
            "ancestral_state": self.ancestral_state,
            "ancestral_state_length": self.ancestral_state_length,
        }
        return SiteTable, tuple(), state

    # Unpickle support
    def __setstate__(self, state):
        self.__init__()
        self.set_columns(
            position=state["position"], ancestral_state=state["ancestral_state"],
            ancestral_state_length=state["ancestral_state_length"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = SiteTable()
        copy.set_columns(
            position=self.position, ancestral_state=self.ancestral_state,
            ancestral_state_length=self.ancestral_state_length)
        return copy


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
        derived_state = unpack_strings(
            self.derived_state, self.derived_state_length)
        ret = "id\tsite\tnode\tderived_state\n"
        for j in range(self.num_rows):
            ret += "{}\t{}\t{}\t{}\n".format(j, site[j], node[j], derived_state[j])
        return ret[:-1]

    def __eq__(self, other):
        ret = False
        if type(other) is type(self):
            ret = (
                np.array_equal(self.site, other.site) and
                np.array_equal(self.node, other.node) and
                np.array_equal(self.derived_state, other.derived_state) and
                np.array_equal(
                    self.derived_state_length, other.derived_state_length))
        return ret

    def __len__(self):
        return self.num_rows

    # Pickle support. See copyreg registration for this function below.
    def _pickle(self):
        state = {
            "site": self.site,
            "node": self.node,
            "derived_state": self.derived_state,
            "derived_state_length": self.derived_state_length,
        }
        return MutationTable, tuple(), state

    # Unpickle support
    def __setstate__(self, state):
        self.__init__()
        self.set_columns(
            site=state["site"], node=state["node"],
            derived_state=state["derived_state"],
            derived_state_length=state["derived_state_length"])

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = MutationTable()
        copy.set_columns(
            site=self.site, node=self.node, derived_state=self.derived_state,
            derived_state_length=self.derived_state_length)
        return copy


# Pickle support for the various tables. We are forced to use copyreg.pickle
# here to support Python 2. For Python 3, we can just use the __setstate__.
copyreg.pickle(NodeTable, NodeTable._pickle)
copyreg.pickle(EdgesetTable, EdgesetTable._pickle)
copyreg.pickle(MigrationTable, MigrationTable._pickle)
copyreg.pickle(SiteTable, SiteTable._pickle)
copyreg.pickle(MutationTable, MutationTable._pickle)


#############################################
# Table functions.
#############################################


def sort_tables(*args, **kwargs):
    """
    Sorts the given tables in place, as follows:

    Edgesets are ordered by

    - time of parent, then
    - parent node ID, then
    - left endpoint.

    For each edgeset, the ``children`` are sorted by increasing node ID.

    Sites are ordered by position, and Mutations are ordered by site.

    Note: for general edgeset tables this only defines a partial ordering, but
    for strict tables (namely, those for which edgesets belonging to a given
    parent do not overlap) this enforces a complete ordering.

    .. todo:: Update this documentation to describe the keyword arguments and
       combinations that are allowed.

    :param NodeTable nodes:
    :param EdgesetTable edgesets:
    :param MigrationTable migrations:
    :param SiteTable sites:
    :param MutationTable mutations:
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
    :param EdgesetTable edgesets: The NodeTable to be simplified.
    :param MigrationTable migrations: The MigrationTable to be simplified.
    :param SiteTable sites: The SiteTable to be simplified.
    :param MutationTable mutations: The MutationTable to be simplified.
    :param bool filter_invariant_sites: Whether to remove sites that have no
        mutations from the output (default: True).
    """
    return _msprime.simplify_tables(*args, **kwargs)


def pack_strings(strings):
    """
    Packs the specified list of strings into a flattened numpy array of characters
    and corresponding lengths.
    """
    lengths = np.array([len(s) for s in strings], dtype=np.uint32)
    encoded = ("".join(strings)).encode()
    return np.fromstring(encoded, dtype=np.int8), lengths


def unpack_strings(packed, length):
    """
    Unpacks a list of string from the specified numpy arrays of packed character
    data and corresponding lengths.
    """
    # This could be done a lot more efficiently...
    ret = []
    offset = 0
    for l in length:
        raw = packed[offset: offset + l].tostring()
        ret.append(raw.decode())
        offset += l
    return ret
