#
# Copyright (C) 2015-2017 University of Oxford
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
Module responsible to generating and reading tree files.
"""
from __future__ import division
from __future__ import print_function

import collections
import gzip
import json
import math
import random
import sys

try:
    import svgwrite
    _svgwrite_imported = True
except ImportError:
    _svgwrite_imported = False

try:
    import numpy as np
    _numpy_imported = True
except ImportError:
    _numpy_imported = False


import _msprime
import msprime.environment

# Make the low-level generator appear like its from this module
from _msprime import RandomGenerator
from _msprime import MutationGenerator
from _msprime import NODE_IS_SAMPLE
from _msprime import sort_tables  # NOQA
from _msprime import simplify_tables  # NOQA

NULL_NODE = -1

NULL_POPULATION = -1


def check_numpy():
    if not _numpy_imported:
        raise RuntimeError("numpy is required for this operation.")


CoalescenceRecord = collections.namedtuple(
    "CoalescenceRecord",
    ["left", "right", "node", "children", "time", "population"])


Migration = collections.namedtuple(
    "Migration",
    ["left", "right", "node", "source", "dest", "time"])


Site = collections.namedtuple(
    "Site",
    ["position", "ancestral_state", "index", "mutations"])


Mutation = collections.namedtuple(
    "Mutation",
    ["site", "node", "derived_state"])


# This is provided for backwards compatibility with the deprecated mutations()
# iterator.
DeprecatedMutation = collections.namedtuple(
    "DeprecatedMutation",
    ["position", "node", "index"])


# This form was chosen to try to break as little of existing code as possible.
# It seems likely that the `node` member was not used for much in this
# iterator, so that it shouldn't break too much if we replace this with the
# underlying site. The position and index values are retained for compatability,
# although they are not redundant. There is a wider question about what the
# genotypes values should be when we have non binary mutations, and whether
# this is a viable long-term API.
Variant = collections.namedtuple(
    "Variant",
    ["position", "site", "index", "genotypes"])


Sample = collections.namedtuple(
    "Sample",
    ["population", "time"])


TableTuple = collections.namedtuple(
    "TableTuple",
    ["nodes", "edgesets", "migrations", "sites", "mutations"])


class SimpleContainer(object):

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return repr(self.__dict__)


class Node(SimpleContainer):
    def __init__(
            self, time=0, population=NULL_POPULATION, name="", is_sample=False):
        self.time = time
        self.population = population
        self.name = name
        self.flags = 0
        if is_sample:
            self.flags |= NODE_IS_SAMPLE

    def is_sample(self):
        return self.flags & NODE_IS_SAMPLE


class Edgeset(SimpleContainer):
    def __init__(self, parent, children, left=0, right=1):
        self.left = left
        self.right = right
        self.parent = parent
        self.children = children


#############################################
# Table definitions. The classes here are thin wrappers for the low-level classes,
# and basically just exist for debugging purposes.
#############################################


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

    TODO: `TreeSequence.simplify()` will accept edgesets not satisfying the
    fourth requirement, producing a `TreeSequence` whose edgesets are of this form.
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

    def copy(self):
        """
        Returns a deep copy of this table.
        """
        copy = MutationTable()
        copy.set_columns(
            site=self.site, node=self.node, derived_state=self.derived_state,
            derived_state_length=self.derived_state_length)
        return copy


def pack_strings(strings):
    """
    Packs the specified list of strings into a flattened numpy array of characters
    and corresponding lengths.
    """
    check_numpy()
    lengths = np.array([len(s) for s in strings], dtype=np.uint32)
    encoded = ("".join(strings)).encode()
    return np.fromstring(encoded, dtype=np.int8), lengths


def unpack_strings(packed, length):
    """
    Unpacks a list of string from the specified numpy arrays of packed character
    data and corresponding lengths.
    """
    # This could be done a lot more efficiently...
    check_numpy()
    ret = []
    offset = 0
    for l in length:
        raw = packed[offset: offset + l].tostring()
        ret.append(raw.decode())
        offset += l
    return ret


def almost_equal(a, b, rel_tol=1e-9, abs_tol=0.0):
    """
    Returns true if the specified pair of integers are equal to
    within the specified tolerances.

    The signature and implementation are taken from PEP 485,
    https://www.python.org/dev/peps/pep-0485/
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def get_provenance_dict(command, parameters):
    """
    Returns a dictionary encoding an execution of msprime.

    Note: this format is incomplete and provisional.
    """
    document = {
        "software": "msprime",
        "version": msprime.environment.__version__,
        "command": command,
        "parameters": parameters,
        "environment": msprime.environment.get_environment()
    }
    return document


class TreeDrawer(object):
    """
    A class to draw sparse trees in SVG format.
    """
    def __init__(
            self, tree, width=200, height=200, show_times=False,
            show_mutation_locations=True, show_mutation_labels=False,
            show_internal_node_labels=True, show_leaf_node_labels=True):
        self._width = width
        self._height = height
        self._show_times = show_times
        self._show_mutation_locations = show_mutation_locations
        self._show_mutation_labels = show_mutation_labels
        self._show_internal_node_labels = show_internal_node_labels
        self._show_leaf_node_labels = show_leaf_node_labels
        self._x_scale = width / (tree.get_sample_size() + 2)
        t = tree.get_time(tree.get_root())
        # Leave a margin of 20px top and bottom
        y_padding = 20
        self._y_scale = (height - 2 * y_padding) / t
        self._tree = tree
        self._x_coords = {}
        self._y_coords = {}
        for u in tree.nodes():
            scaled_t = tree.get_time(u) * self._y_scale
            self._y_coords[u] = height - scaled_t - y_padding
        self._sample_x = 1
        self._assign_x_coordinates(self._tree.get_root())
        self._mutations = []
        node_mutations = collections.defaultdict(list)
        for site in tree.sites():
            for mutation in site.mutations:
                node_mutations[mutation.node].append(mutation)
        for child, mutations in node_mutations.items():
            n = len(mutations)
            parent = tree.parent(child)
            x = self._x_coords[child]
            y1 = self._y_coords[child]
            y2 = self._y_coords[parent]
            chunk = (y2 - y1) / (n + 1)
            for k, mutation in enumerate(mutations):
                z = x, y1 + (k + 1) * chunk
                self._mutations.append((z, mutation))

    def draw(self):
        """
        Writes the SVG description of this tree and returns the resulting XML
        code as text.
        """
        dwg = svgwrite.Drawing(size=(self._width, self._height), debug=True)
        lines = dwg.add(dwg.g(id='lines', stroke='black'))
        labels = dwg.add(dwg.g(font_size=14, text_anchor="middle"))
        for u in self._tree.nodes():
            v = self._tree.get_parent(u)
            x = self._x_coords[u], self._y_coords[u]
            dwg.add(dwg.circle(center=x, r=3))
            dx = [0]
            dy = None
            if self._tree.is_sample(u):
                dy = [20]
            elif v == msprime.NULL_NODE:
                dy = [-5]
            else:
                dx = [-10]
                dy = [-5]
            condition = (
                (self._tree.is_sample(u) and self._show_leaf_node_labels) or
                (self._tree.is_internal(u) and self._show_internal_node_labels))
            if condition:
                labels.add(dwg.text(str(u), x, dx=dx, dy=dy))
            if self._show_times and self._tree.is_internal(u):
                dx[0] += 25
                labels.add(dwg.text(
                    "t = {:.2f}".format(self._tree.get_time(u)), x, dx=dx,
                    dy=dy))
            if v != NULL_NODE:
                y = self._x_coords[v], self._y_coords[v]
                lines.add(dwg.line(x, (x[0], y[1])))
                lines.add(dwg.line((x[0], y[1]), y))
        for x, mutation in self._mutations:
            r = 3
            if self._show_mutation_locations:
                dwg.add(dwg.rect(
                    insert=(x[0] - r, x[1] - r), size=(2 * r, 2 * r), fill="red"))
            if self._show_mutation_labels:
                dx = [8 * r]
                dy = [-2 * r]
                labels.add(dwg.text("{}".format(mutation.site), x, dx=dx, dy=dy))
        return dwg.tostring()

    def _assign_x_coordinates(self, node):
        """
        Assign x coordinates to all nodes underneath this node.
        """
        if self._tree.is_internal(node):
            children = self._tree.get_children(node)
            for c in children:
                self._assign_x_coordinates(c)
            coords = [self._x_coords[c] for c in children]
            a = min(coords)
            b = max(coords)
            self._x_coords[node] = (a + (b - a) / 2)
        else:
            self._x_coords[node] = self._sample_x * self._x_scale
            self._sample_x += 1


# TODO:
# - Pickle and copy support
class SparseTree(object):
    """
    A SparseTree is a single tree in a :class:`.TreeSequence`. In a sparse tree
    for a sample of size :math:`n`, the samples are nodes :math:`0` to :math:`n
    - 1` inclusive and internal nodes are integers :math:`\geq n`. The value of
    these nodes is strictly increasing as we ascend the tree and the root of
    the tree is the node with the largest value that is reachable from  the
    samples. Each node in the tree has a parent which is obtained using the
    :meth:`.get_parent` method. The parent of the root node is the
    :const:`.NULL_NODE`, :math:`-1`. Similarly, each internal node has a
    pair of children, which are obtained using the :meth:`.get_children`
    method. Each node in the tree has a time associated with it in generations.
    This value is obtained using the :meth:`.SparseTree.get_time` method.

    Sparse trees are not intended to be instantiated directly, and are
    obtained as part of a :class:`.TreeSequence` using the
    :meth:`.trees` method.
    """
    def __init__(self, ll_sparse_tree):
        self._ll_sparse_tree = ll_sparse_tree

    def branch_length(self, u):
        return self.get_branch_length(u)

    def get_branch_length(self, u):
        """
        Returns the length of the branch (in generations) joining the
        specified node to its parent. This is equivalent to

        >>> tree.get_time(tree.get_parent(u)) - tree.get_time(u)

        Note that this is not related to the value returned by
        :meth:`.get_length`, which describes the length of the interval
        covered by the tree in genomic coordinates.

        :param int u: The node of interest.
        :return: The branch length from u to its parent.
        :rtype: float
        """
        return self.get_time(self.get_parent(u)) - self.get_time(u)

    @property
    def total_branch_length(self):
        return self.get_total_branch_length()

    def get_total_branch_length(self):
        """
        Returns the sum of all the branch lengths in this tree (in
        units of generations). This is equivalent to

        >>> sum(
        >>>    tree.get_branch_length(u) for u in tree.nodes()
        >>>    if u != tree.get_root())

        :return: The sum of all the branch lengths in this tree.
        """
        root = self.get_root()
        return sum(
            self.get_branch_length(u) for u in self.nodes() if u != root)

    def mrca(self, u, v):
        return self.get_mrca(u, v)

    def get_mrca(self, u, v):
        """
        Returns the most recent common ancestor of the specified nodes.

        :param int u: The first node.
        :param int v: The second node.
        :return: The most recent common ancestor of u and v.
        :rtype: int
        """
        return self._ll_sparse_tree.get_mrca(u, v)

    def tmrca(self, u, v):
        return self.get_tmrca(u, v)

    def get_tmrca(self, u, v):
        """
        Returns the time of the most recent common ancestor of the specified
        nodes. This is equivalent to::

            tree.get_time(tree.get_mrca(u, v))

        :param int u: The first node.
        :param int v: The second node.
        :return: The time of the most recent common ancestor of u and v.
        :rtype: float
        """
        return self.get_time(self.get_mrca(u, v))

    def parent(self, u):
        return self.get_parent(u)

    def get_parent(self, u):
        """
        Returns the parent of the specified node. Returns
        the :const:`.NULL_NODE` if u is the root or is not a node in
        the current tree.

        :param int u: The node of interest.
        :return: The parent of u.
        :rtype: int
        """
        return self._ll_sparse_tree.get_parent(u)

    def children(self, u):
        return self.get_children(u)

    def get_children(self, u):
        """
        Returns the children of the specified node as a tuple :math:`(v, w)`.
        For internal nodes, this tuple is always in sorted order such that
        :math:`v < w`. If u is a leaf or is not a node in the current tree,
        return the empty tuple.

        :param int u: The node of interest.
        :return: The children of u as a pair of integers
        :rtype: tuple
        """
        return self._ll_sparse_tree.get_children(u)

    def time(self, u):
        return self.get_time(u)

    def get_time(self, u):
        """
        Returns the time of the specified node in generations. Returns 0 if u
        is not a node in the current tree.

        :param int u: The node of interest.
        :return: The time of u.
        :rtype: float
        """
        return self._ll_sparse_tree.get_time(u)

    def population(self, u):
        return self.get_population(u)

    def get_population(self, u):
        """
        Returns the population associated with the specified node. If the
        specified node is not a member of this tree or population level
        information was not stored in the tree sequence,
        :const:`.NULL_POPULATION` is returned.

        :param int u: The node of interest.
        :return: The ID of the population associated with node u.
        :rtype: int
        """
        return self._ll_sparse_tree.get_population(u)

    def is_internal(self, u):
        """
        Returns True if the specified node is not a leaf. A node is internal
        if it has one or more children in the current tree.

        :param int u: The node of interest.
        :return: True if u is not a leaf node.
        :rtype: bool
        """
        return not self.is_leaf(u)

    def is_leaf(self, u):
        """
        Returns True if the specified node is a leaf. A node :math:`u` is a
        leaf if it has zero children.

        :param int u: The node of interest.
        :return: True if u is a leaf node.
        :rtype: bool
        """
        return len(self.children(u)) == 0

    def is_sample(self, u):
        """
        Returns True if the specified node is a sample. A node :math:`u` is a
        sample if it has been marked as a sample in the parent tree sequence.

        :param int u: The node of interest.
        :return: True if u is a sample.
        :rtype: bool
        """
        return bool(self._ll_sparse_tree.is_sample(u))

    @property
    def num_nodes(self):
        # TODO documnent
        return self._ll_sparse_tree.get_num_nodes()

    @property
    def root(self):
        return self.get_root()

    def get_root(self):
        """
        Returns the root of this tree.

        :return: The root node.
        :rtype: int
        """
        return self._ll_sparse_tree.get_root()

    @property
    def index(self):
        return self.get_index()

    def get_index(self):
        """
        Returns the index this tree occupies in the parent tree sequence.
        This index is zero based, so the first tree in the sequence has index
        0.

        :return: The index of this tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_index()

    @property
    def interval(self):
        return self.get_interval()

    def get_interval(self):
        """
        Returns the coordinates of the genomic interval that this tree
        represents the history of. The interval is returned as a tuple
        :math:`(l, r)` and is a half-open interval such that the left
        coordinate is inclusive and the right coordinate is exclusive. This
        tree therefore applies to all genomic locations :math:`x` such that
        :math:`l \leq x < r`.

        :return: A tuple (l, r) representing the left-most (inclusive)
            and right-most (exclusive) coordinates of the genomic region
            covered by this tree.
        :rtype: tuple
        """
        return (
            self._ll_sparse_tree.get_left(), self._ll_sparse_tree.get_right()
        )

    @property
    def length(self):
        return self.get_length()

    def get_length(self):
        """
        Returns the length of the genomic interval that this tree represents.
        This is defined as :math:`r - l`, where :math:`(l, r)` is the genomic
        interval returned by :meth:`.get_interval`.

        :return: The length of the genomic interval covered by this tree.
        :rtype: int
        """
        l, r = self.get_interval()
        return r - l

    @property
    def sample_size(self):
        return self.get_sample_size()

    def get_sample_size(self):
        """
        Returns the sample size for this tree. This is the number of sample
        nodes in the tree.

        :return: The number of sample nodes in the tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_sample_size()

    def draw(
            self, path=None, width=200, height=200, times=False,
            mutation_locations=True, mutation_labels=False,
            internal_node_labels=True, leaf_node_labels=True, show_times=None):
        """
        Returns a representation of this tree in SVG format.

        :param str path: The path to the file to write the SVG. If None, do not
            write to file.
        :param int width: The width of the image in pixels.
        :param int height: The height of the image in pixels.
        :param bool times: If True, show time labels at each internal node.
        :param bool mutation_locations: If True, show mutations as points over nodes.
        :param bool mutation_labels: If True, show labels for mutations.
        :param bool internal_node_labels: If True, show labels for internal nodes.
        :param bool leaf_node_labels: If True, show labels for leaf nodes.
        :param bool show_times: Deprecated alias for ``times``.
        :return: A representation of this tree in SVG format.
        :rtype: str
        """
        # show_times is a deprecated alias for times.
        if show_times is not None:
            times = show_times
        if not _svgwrite_imported:
            raise ImportError(
                "svgwrite is not installed. try `pip install svgwrite`")
        td = TreeDrawer(
                self, width=width, height=height, show_times=times,
                show_mutation_locations=mutation_locations,
                show_mutation_labels=mutation_labels,
                show_internal_node_labels=internal_node_labels,
                show_leaf_node_labels=leaf_node_labels)
        svg = td.draw()
        if path is not None:
            with open(path, "w") as f:
                f.write(svg)
        return svg

    @property
    def num_mutations(self):
        return self.get_num_mutations()

    def get_num_mutations(self):
        """
        Returns the number of mutations on this tree.

        :return: The number of mutations on this tree.
        :rtype: int
        """
        return sum(len(site.mutations) for site in self.sites())

    def sites(self):
        for ll_site in self._ll_sparse_tree.get_sites():
            pos, ancestral_state, mutations, index = ll_site
            yield Site(
                position=pos, ancestral_state=ancestral_state, index=index,
                mutations=[Mutation(*mutation) for mutation in mutations])

    def mutations(self):
        """
        Returns an iterator over the mutations in this tree. Each
        mutation is represented as a tuple :math:`(x, u, j)` where :math:`x`
        is the position of the mutation in the sequence in chromosome
        coordinates, :math:`u` is the node over which the mutation
        occurred and :math:`j` is the zero-based index of the mutation within
        the overall tree sequence. Mutations are returned in non-decreasing
        order of position and increasing index.

        Each mutation returned is an instance of
        :func:`collections.namedtuple`, and may be accessed via the attributes
        ``position``, ``node`` and ``index`` as well as the usual positional
        approach. This is the recommended interface for working with mutations
        as it is both more readable and also ensures that code is forward
        compatible with future extensions.

        :return: An iterator of all :math:`(x, u, j)` tuples defining
            the mutations in this tree.
        :rtype: iter
        """
        # TODO deprecate
        for site in self.sites():
            for mutation in site.mutations:
                yield DeprecatedMutation(
                    position=site.position, node=mutation.node, index=site.index)

    def get_leaves(self, u):
        # Deprecated alias for samples. See the discussion in the get_num_leaves
        # method for why this method is here and why it is semantically incorrect.
        # The 'leaves' iterator below correctly returns the leaves below a given
        # node.
        return self.samples(u)

    def leaves(self, u=None):
        """
        Returns an iterator over all the leaves in this tree that are
        underneath the specified node. If u is not specified, return all leaves
        in the tree.

        :param int u: The node of interest.
        :return: An iterator over all leaves in the subtree rooted at u.
        :rtype: iterator
        """
        if u is None:
            u = self.root
        for v in self.nodes(u):
            if self.is_leaf(v):
                yield v

    def _sample_generator(self, u):
        for v in self.nodes(u):
            if self.is_sample(v):
                yield v

    def samples(self, u=None):
        """
        Returns an iterator over all the samples in this tree that are
        underneath the specified node. If u is a sample, it is included in the
        returned iterator. If u is not specified, return all samples in the tree.

        If the :meth:`.TreeSequence.trees` method is called with
        ``sample_lists=True``, this method uses an efficient algorithm to find
        the samples. If not, a simple traversal based method is used.

        :param int u: The node of interest.
        :return: An iterator over all samples in the subtree rooted at u.
        :rtype: iterator
        """
        if u is None:
            u = self.root
        if self._ll_sparse_tree.get_flags() & _msprime.SAMPLE_LISTS:
            return _msprime.SampleListIterator(self._ll_sparse_tree, u)
        else:
            return self._sample_generator(u)

    def get_num_leaves(self, u):
        # Deprecated alias for num_samples. The method name is inaccurate
        # as this will count the number of tracked _samples_. This is only provided to
        # avoid breaking existing code and should not be used in new code. We could
        # change this method to be semantically correct and just count the
        # number of leaves we hit in the leaves() iterator. However, this would
        # have the undesirable effect of making code that depends on the constant
        # time performance of get_num_leaves many times slower. So, the best option
        # is to leave this method as is, and to slowly deprecate it out. Once this
        # has been removed, we might add in a ``num_leaves`` method that returns the
        # length of the leaves() iterator as one would expect.
        return self.num_samples(u)

    def num_samples(self, u=None):
        return self.get_num_samples(u)

    def get_num_samples(self, u=None):
        """
        Returns the number of samples in this tree underneath the specified
        node (including the node itself). If u is not specified return
        the total number of samples in the tree.

        If the :meth:`.TreeSequence.trees` method is called with
        ``sample_counts=True`` this method is a constant time operation. If not,
        a slower traversal based algorithm is used to count the samples.

        :param int u: The node of interest.
        :return: The number of samples in the subtree rooted at u.
        :rtype: int
        """
        if u is None:
            u = self.root
        return self._ll_sparse_tree.get_num_samples(u)

    def get_num_tracked_leaves(self, u):
        # Deprecated alias for num_tracked_samples. The method name is inaccurate
        # as this will count the number of tracked _samples_. This is only provided to
        # avoid breaking existing code and should not be used in new code.
        return self.num_tracked_samples(u)

    def num_tracked_samples(self, u=None):
        return self.get_num_tracked_samples(u)

    def get_num_tracked_samples(self, u=None):
        """
        Returns the number of samples in the set specified in the
        ``tracked_samples`` parameter of the :meth:`.TreeSequence.trees` method
        underneath the specified node. If the input node is not specified,
        return the total number of tracked samples in the tree.

        This is a constant time operation.

        :param int u: The node of interest.
        :return: The number of samples within the set of tracked samples in
            the subtree rooted at u.
        :rtype: int
        :raises RuntimeError: if the :meth:`.TreeSequence.trees`
            method is not called with ``sample_counts=True``.
        """
        if u is None:
            u = self.root
        if not (self._ll_sparse_tree.get_flags() & _msprime.SAMPLE_COUNTS):
            raise RuntimeError(
                "The get_num_tracked_samples method is only supported "
                "when sample_counts=True.")
        return self._ll_sparse_tree.get_num_tracked_samples(u)

    def _preorder_traversal(self, u):
        stack = [u]
        while len(stack) > 0:
            v = stack.pop()
            if self.is_internal(v):
                stack.extend(reversed(self.get_children(v)))
            yield v

    def _postorder_traversal(self, u):
        stack = [u]
        k = NULL_NODE
        while stack:
            v = stack[-1]
            if self.is_internal(v) and v != k:
                stack.extend(reversed(self.get_children(v)))
            else:
                k = self.get_parent(v)
                yield stack.pop()

    def _inorder_traversal(self, u):
        stack = [u]
        k, j = NULL_NODE, NULL_NODE
        while stack:
            v = stack.pop()
            if self.is_internal(v) and v != k and v != j:
                children = self.get_children(v)
                j = stack[-1] if stack else NULL_NODE
                stack.extend([children[1], v, children[0]])
            else:
                k = self.get_parent(v)
                yield v

    def _levelorder_traversal(self, u):
        queue = collections.deque([u])
        while queue:
            v = queue.popleft()
            if self.is_internal(v):
                queue.extend(self.get_children(v))
            yield v

    def nodes(self, root=None, order="preorder"):
        """
        Returns an iterator over the nodes in this tree. If the root parameter
        is provided, iterate over the nodes in the subtree rooted at this
        node. If this is None, iterate over all nodes. If the order parameter
        is provided, iterate over the nodes in required tree traversal order.

        :param int root: The root of the subtree we are traversing.
        :param str order: The traversal ordering. Currently 'preorder',
            'inorder', 'postorder' and 'levelorder' ('breadthfirst')
            are supported.
        :rtype: iterator
        """
        u = self.get_root() if root is None else root
        if order == "preorder":
            return self._preorder_traversal(u)
        elif order == "postorder":
            return self._postorder_traversal(u)
        elif order == "inorder":
            return self._inorder_traversal(u)
        elif order == "levelorder" or order == "breadthfirst":
            return self._levelorder_traversal(u)
        else:
            raise ValueError(
                "Traversal ordering '{}' not supported".format(order))

    @property
    def parent_dict(self):
        return self.get_parent_dict()

    def get_parent_dict(self):
        pi = {
            u: self.parent(u) for u in range(self.num_nodes)
            if self.parent(u) != NULL_NODE}
        return pi

    @property
    def time_dict(self):
        return self.get_time_dict()

    def get_time_dict(self):
        return {
            u: self.time(u) for u in range(self.num_nodes)
            if len(self.children(u)) != 0 or self.parent(u) != NULL_NODE}

    def __str__(self):
        return str(self.get_parent_dict())

    def __eq__(self, other):
        # TODO this should really use the semantics of the C-level equality definition.
        # This is currently only really used for testing AFAIK.
        return (
            self.get_sample_size() == other.get_sample_size() and
            self.get_parent_dict() == other.get_parent_dict() and
            self.get_time_dict() == other.get_time_dict() and
            self.get_interval() == other.get_interval() and
            self.get_root() == other.get_root() and
            self.get_index() == other.get_index() and
            list(self.sites()) == list(other.sites()))

    def __ne__(self, other):
        return not self.__eq__(other)


def _get_random_seed():
    return random.randint(1, 2**32 - 1)


def _check_population_configurations(population_configurations):
    err = (
        "Population configurations must be a list of "
        "PopulationConfiguration instances")
    if not isinstance(population_configurations, collections.Iterable):
        raise TypeError(err)
    for config in population_configurations:
        if not isinstance(config, PopulationConfiguration):
            raise TypeError(err)


def _replicate_generator(
        sim, mutation_generator, num_replicates, provenance_dict):
    """
    Generator function for the many-replicates case of the simulate
    function.
    """
    # TODO like in the single replicate case, we need to encode the
    # simulation parameters so that particular simulations can be
    # replicated. This will also involve encoding the state of the
    # random generator.
    provenance = [json.dumps(provenance_dict).encode()]
    # Should use range here, but Python 2 makes this awkward...
    j = 0
    while j < num_replicates:
        j += 1
        sim.run()
        tree_sequence = sim.get_tree_sequence(mutation_generator, provenance)
        yield tree_sequence
        sim.reset()


def simulator_factory(
        sample_size=None,
        Ne=1,
        random_generator=None,
        length=None,
        recombination_rate=None,
        recombination_map=None,
        population_configurations=None,
        migration_matrix=None,
        samples=None,
        demographic_events=[],
        model=None,
        record_migrations=False):
    """
    Convenience method to create a simulator instance using the same
    parameters as the `simulate` function. Primarily used for testing.
    """
    condition = (
        sample_size is None and
        population_configurations is None and
        samples is None)
    if condition:
        raise ValueError(
            "Either sample_size, population_configurations or samples "
            "be specified")
    if sample_size is not None:
        if samples is not None:
            raise ValueError(
                "Cannot specify sample size and samples simultaneously.")
        if population_configurations is not None:
            raise ValueError(
                "Cannot specify sample size and population_configurations "
                "simultaneously.")
        s = Sample(population=0, time=0.0)
        the_samples = [s for _ in range(sample_size)]
    # If we have population configurations we may have embedded sample_size
    # values telling us how many samples to take from each population.
    if population_configurations is not None:
        _check_population_configurations(population_configurations)
        if samples is None:
            the_samples = []
            for j, conf in enumerate(population_configurations):
                if conf.sample_size is not None:
                    the_samples += [(j, 0) for _ in range(conf.sample_size)]
        else:
            for conf in population_configurations:
                if conf.sample_size is not None:
                    raise ValueError(
                        "Cannot specify population configuration sample size"
                        "and samples simultaneously")
            the_samples = samples
    elif samples is not None:
        the_samples = samples

    if recombination_map is None:
        the_length = 1 if length is None else length
        the_rate = 0 if recombination_rate is None else recombination_rate
        if the_length <= 0:
            raise ValueError("Cannot provide non-positive sequence length")
        if the_rate < 0:
            raise ValueError("Cannot provide negative recombination rate")
        recomb_map = RecombinationMap.uniform_map(the_length, the_rate)
    else:
        if length is not None or recombination_rate is not None:
            raise ValueError(
                "Cannot specify length/recombination_rate along with "
                "a recombination map")
        recomb_map = recombination_map

    sim = TreeSimulator(the_samples, recomb_map)
    sim.set_store_migrations(record_migrations)
    sim.set_effective_population_size(Ne)
    sim.set_model(model)
    rng = random_generator
    if rng is None:
        rng = RandomGenerator(_get_random_seed())
    sim.set_random_generator(rng)
    if population_configurations is not None:
        sim.set_population_configurations(population_configurations)
    if migration_matrix is not None:
        sim.set_migration_matrix(migration_matrix)
    if demographic_events is not None:
        sim.set_demographic_events(demographic_events)
    return sim


def simulate(
        sample_size=None,
        Ne=1,
        length=None,
        recombination_rate=None,
        recombination_map=None,
        mutation_rate=None,
        population_configurations=None,
        migration_matrix=None,
        demographic_events=[],
        samples=None,
        model=None,
        record_migrations=False,
        random_seed=None,
        mutation_generator=None,
        num_replicates=None):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns the resulting :class:`.TreeSequence`.

    :param int sample_size: The number of individuals in our sample.
        If not specified or None, this defaults to the sum of the
        subpopulation sample sizes. Either ``sample_size``,
        ``population_configurations`` or ``samples`` must be specified.
    :param float Ne: The effective (diploid) population size for the reference
        population. This determines the factor by which the per-generation
        recombination and mutation rates are scaled in the simulation.
        This defaults to 1 if not specified.
    :param float length: The length of the simulated region in bases.
        This parameter cannot be used along with ``recombination_map``.
        Defaults to 1 if not specified.
    :param float recombination_rate: The rate of recombination per base
        per generation. This parameter cannot be used along with
        ``recombination_map``. Defaults to 0 if not specified.
    :param recombination_map: The map
        describing the changing rates of recombination along the simulated
        chromosome. This parameter cannot be used along with the
        ``recombination_rate`` or ``length`` parameters, as these
        values are encoded within the map. Defaults to a uniform rate as
        described in the ``recombination_rate`` parameter if not specified.
    :type recombination_map: :class:`.RecombinationMap`
    :param float mutation_rate: The rate of mutation per base per
        generation. If not specified, no mutations are generated.
    :param list population_configurations: The list of
        :class:`.PopulationConfiguration` instances describing the
        sampling configuration, relative sizes and growth rates of
        the populations to be simulated. If this is not specified,
        a single population with a sample of size ``sample_size``
        is assumed.
    :type population_configurations: list or None.
    :param list migration_matrix: The matrix describing the rates
        of migration between all pairs of populations. If :math:`N`
        populations are defined in the ``population_configurations``
        parameter, then the migration matrix must be an
        :math:`N\\times N` matrix consisting of :math:`N` lists of
        length :math:`N`.
    :param list demographic_events: The list of demographic events to
        simulate. Demographic events describe changes to the populations
        in the past. Events should be supplied in non-decreasing
        order of time. Events with the same time value will be applied
        sequentially in the order that they were supplied before the
        simulation algorithm continues with the next time step.
    :param list samples: The list specifying the location and time of
        all samples. This parameter may be used to specify historical
        samples, and cannot be used in conjunction with the ``sample_size``
        parameter. Each sample is a (``population_id``, ``time``) pair
        such that the sample in position ``j`` in the list of samples
        is drawn in the specified population at the specfied time. Time
        is measured in generations, as elsewhere.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated. Valid random
        seeds must be between 1 and :math:`2^{32} - 1`.
    :param int num_replicates: The number of replicates of the specified
        parameters to simulate. If this is not specified or None,
        no replication is performed and a :class:`.TreeSequence` object
        returned. If :obj:`num_replicates` is provided, the specified
        number of replicates is performed, and an iterator over the
        resulting :class:`.TreeSequence` objects returned.
    :return: The :class:`.TreeSequence` object representing the results
        of the simulation if no replication is performed, or an
        iterator over the independent replicates simulated if the
        :obj:`num_replicates` parameter has been used.
    :rtype: :class:`.TreeSequence` or an iterator over
        :class:`.TreeSequence` replicates.
    :warning: If using replication, do not store the results of the
        iterator in a list! For performance reasons, the same
        underlying object may be used for every TreeSequence
        returned which will most likely lead to unexpected behaviour.
    """
    seed = random_seed
    if random_seed is None:
        seed = _get_random_seed()
    rng = RandomGenerator(seed)
    sim = simulator_factory(
        sample_size=sample_size, random_generator=rng,
        Ne=Ne, length=length,
        recombination_rate=recombination_rate,
        recombination_map=recombination_map,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        samples=samples, model=model, record_migrations=record_migrations)
    # The provenance API is very tentative, and only included now as a
    # pre-alpha feature.
    parameters = {"TODO": "encode simulation parameters"}
    provenance = get_provenance_dict("simulate", parameters)
    if mutation_generator is None:
        mu = 0 if mutation_rate is None else mutation_rate
        mutation_generator = MutationGenerator(rng, mu)
    else:
        if mutation_rate is not None:
            raise ValueError(
                "Cannot specify both mutation_rate and mutation_generator")
    if num_replicates is None:
        return next(_replicate_generator(sim, mutation_generator, 1, provenance))
    else:
        return _replicate_generator(sim, mutation_generator, num_replicates, provenance)


def load(path):
    """
    Loads a tree sequence from the specified file path. This
    file must be in the HDF5 file format produced by the
    :meth:`.TreeSequence.dump` method.

    :param str path: The file path of the HDF5 file containing the
        tree sequence we wish to load.
    :return: The tree sequence object containing the information
        stored in the specified file path.
    :rtype: :class:`msprime.TreeSequence`
    """
    return TreeSequence.load(path)


def load_tables(*args, **kwargs):
    return TreeSequence.load_tables(*args, **kwargs)


def parse_nodes(source):
    """
    Parse the specified file-like object and return a NodeTable instance.  The
    object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``is_sample``, ``time``, and
    optionally, ``population``.  Further requirements are described in
    :class:`NodeTable`.  Note that node ``id`` is not included, but implied by
    order in the file.
    """
    # Read the header and find the indexes of the required fields.
    table = NodeTable()
    header = source.readline().split()
    is_sample_index = header.index("is_sample")
    time_index = header.index("time")
    population_index = None
    try:
        population_index = header.index("population")
    except ValueError:
        pass
    for line in source:
        tokens = line.split()
        if len(tokens) > 0:
            is_sample = int(tokens[is_sample_index])
            time = float(tokens[time_index])
            flags = 0
            if is_sample != 0:
                flags |= NODE_IS_SAMPLE
            population = NULL_POPULATION
            if population_index is not None:
                population = int(tokens[population_index])
            table.add_row(flags=flags, time=time, population=population)
    return table


def parse_edgesets(source):
    """
    Parse the specified file-like object and return a EdgesetTable instance.
    The object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``left``, ``right``, ``parent``,
    and ``children``.  The ``children`` field is a comma-separated list of base
    10 integer values.  Further requirements are described in
    :class:`EdgesetTable`.
    """
    table = EdgesetTable()
    header = source.readline().split()
    left_index = header.index("left")
    right_index = header.index("right")
    parent_index = header.index("parent")
    children_index = header.index("children")
    table = EdgesetTable()
    for line in source:
        tokens = line.split()
        if len(tokens) > 0:
            left = float(tokens[left_index])
            right = float(tokens[right_index])
            parent = int(tokens[parent_index])
            children = tuple(map(int, tokens[children_index].split(",")))
            table.add_row(
                left=left, right=right, parent=parent, children=children)
    return table


def parse_sites(source):
    """
    Parse the specified file-like object and return a SiteTable instance.  The
    object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``position`` and
    ``ancestral_state``.  Further requirements are described in
    :class:`SiteTable`.
    """
    header = source.readline().split()
    position_index = header.index("position")
    ancestral_state_index = header.index("ancestral_state")
    table = SiteTable()
    for line in source:
        tokens = line.split()
        if len(tokens) > 0:
            position = float(tokens[position_index])
            ancestral_state = tokens[ancestral_state_index]
            table.add_row(position=position, ancestral_state=ancestral_state)
    return table


def parse_mutations(source):
    """
    Parse the specified file-like object and return a MutationTable instance.
    The object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``site``, ``node``, and
    ``derived_state``.  Further requirements are described in
    :class:`MutationTable`.
    """
    header = source.readline().split()
    site_index = header.index("site")
    node_index = header.index("node")
    derived_state_index = header.index("derived_state")
    table = MutationTable()
    for line in source:
        tokens = line.split()
        if len(tokens) > 0:
            site = int(tokens[site_index])
            node = int(tokens[node_index])
            derived_state = tokens[derived_state_index]
            table.add_row(site=site, node=node, derived_state=derived_state)
    return table


def load_text(nodes, edgesets, sites=None, mutations=None):
    """
    Loads a tree sequence from the specified file paths. The files input here
    are in a simple whitespace delimited tabular format such as output by the
    :meth:`.TreeSequence.dump_text` method.  This method is intended as a
    convenient interface for importing external data into msprime; the HDF5
    based file format using by :meth:`msprime.load` will be many times more
    efficient that using the text based formats.

    ``nodes`` and ``edgesets`` must be a file-like object containing text with
    whitespace delimited columns,  parsable by :func:`parse_nodes` and
    :func:`parse_edgesets`, respectively.  Further requirements are described
    in :class:`NodeTable` and :class:`EdgesetTable`.

    ``sites`` and ``mutations`` are optional, but if included must be similar,
    parsable by :func:`parse_sites` and :func:`parse_mutations`, respecively.
    Further requirements are described in :class:`SiteTable` and
    :class:`MutationTable`.

    An example of a simple tree sequence for four samples with
    three distinct trees is as follows.

    nodes::

        is_sample   time    population
        1           0.0     0
        1           0.0     0
        1           0.0     0
        1           0.0     0
        0           0.071   0
        0           0.090   0
        0           0.170   0
        0           0.202   0
        0           0.253   0

    edgesets::

        left    right   node    children
        2       10      4       2,3
        0       2       5       1,3
        2       10      5       1,4
        0       7       6       0,5
        7       10      7       0,5
        0       2       8       2,6


    This example is equivalent to the tree sequence illustrated in Figure 4 of
    the `PLoS Computational Biology paper
    <http://dx.doi.org/10.1371/journal.pcbi.1004842>`_. Nodes are given here in
    time order (since this is a backwards-in-time tree sequence), but they may
    be allocated in any order. In particular, left-to-right tree sequences are
    fully supported.

    An example of a ``sites`` and ``mutations`` file for the tree sequence
    defined in the previous example is as follows.

    sites::

        position    ancestral_state
        0.1         0
        8.5         0

    mutations::

        site    node    derived_state
        0       3       1
        1       6       1
        1       0       0


    :param stream nodes: The file-type object containing text describing a NodeTable.
    :param stream edgesets: The file-type object containing text
        describing a EdgesetTable.
    :param stream sites: The file-type object containing text describing a SiteTable.
    :param stream mutations: The file-type object containing text
        describing a MutationTable.
    :return: The tree sequence object containing the information
        stored in the specified file paths.
    :rtype: :class:`msprime.TreeSequence`
    """
    node_table = parse_nodes(nodes)
    edgeset_table = parse_edgesets(edgesets)
    site_table = SiteTable()
    mutation_table = MutationTable()
    if sites is not None:
        site_table = parse_sites(sites)
    if mutations is not None:
        mutation_table = parse_mutations(mutations)
    return load_tables(
        nodes=node_table, edgesets=edgeset_table, sites=site_table,
        mutations=mutation_table)


class TreeSimulator(object):
    """
    Class to simulate trees under the standard neutral coalescent with
    recombination.
    """
    def __init__(self, samples, recombination_map):
        if len(samples) < 2:
            raise ValueError("Sample size must be >= 2")
        if len(samples) >= 2**32:
            raise ValueError("sample_size must be < 2**32")
        self._sample_size = len(samples)
        self._samples = samples
        if not isinstance(recombination_map, RecombinationMap):
            raise TypeError("RecombinationMap instance required")
        self._model = "hudson"
        self._recombination_map = recombination_map
        self._random_generator = None
        self._effective_population_size = 1
        self._population_configurations = [PopulationConfiguration()]
        self._migration_matrix = [[0]]
        self._demographic_events = []
        self._store_migrations = False
        # Set default block sizes to 64K objects.
        # TODO does this give good performance in a range of scenarios?
        block_size = 64 * 1024
        # We always need at least n segments, so no point in making
        # allocation any smaller than this.
        self._segment_block_size = max(block_size, self._sample_size)
        self._avl_node_block_size = block_size
        self._node_mapping_block_size = block_size
        self._coalescence_record_block_size = block_size
        self._migration_block_size = block_size
        # TODO is it useful to bring back the API to set this? Mostly
        # the amount of memory required is tiny.
        self._max_memory = sys.maxsize
        self._ll_sim = None

        # TODO document these public attributes.
        self.node_table = NodeTable(block_size)
        self.edgeset_table = EdgesetTable(block_size)
        self.migration_table = MigrationTable(block_size)
        self.mutation_type_table = SiteTable(1)
        self.mutation_table = MutationTable(block_size)

    def set_random_generator(self, random_generator):
        """
        Sets the random generator instance for this simulator to the specified
        value.
        """
        self._random_generator = random_generator

    def get_random_generator(self):
        return self._random_generator

    def get_sample_size(self):
        return self._sample_size

    def get_samples(self):
        return self._samples

    def get_recombinatation_map(self):
        return self._recombination_map

    def get_effective_population_size(self):
        return self._effective_population_size

    def get_per_locus_scaled_recombination_rate(self):
        """
        Returns the rate of recombination between pair of adjacent
        loci per 4 Ne generations. This is the rate at which recombination
        occurs in the underlying simulation in genetic coordinates. This
        is _not_ the per-base rate of recombination.
        """
        return (
            4 * self.get_effective_population_size() *
            self._recombination_map.get_per_locus_recombination_rate())

    def get_scaled_migration_matrix(self):
        """
        Returns the migratation matrix scaled in coalescent time units.
        """
        return [
            [4 * self.get_effective_population_size() * m for m in row]
            for row in self.get_migration_matrix()]

    def get_model(self):
        return self._model

    def get_migration_matrix(self):
        return self._migration_matrix

    def get_num_loci(self):
        return self._recombination_map.get_num_loci()

    def get_population_configurations(self):
        return self._population_configurations

    def get_sample_configuration(self):
        return [conf.sample_size for conf in self._population_configurations]

    def get_demographic_events(self):
        return self._demographic_events

    def get_num_breakpoints(self):
        return self._ll_sim.get_num_breakpoints()

    def get_breakpoints(self):
        """
        Returns the recombination breakpoints translated into physical
        coordinates.
        """
        return [
            self._recombination_map.genetic_to_physical(x)
            for x in self._ll_sim.get_breakpoints()]

    def get_used_memory(self):
        return self._ll_sim.get_used_memory()

    def get_time(self):
        # The low-level simulation returns time in coalescent units.
        return (
            self.get_effective_population_size() * 4 * self._ll_sim.get_time())

    def get_avl_node_block_size(self):
        return self._avl_node_block_size

    def get_coalescence_record_block_size(self):
        return self._coalescence_record_block_size

    def get_node_mapping_block_size(self):
        return self._node_mapping_block_size

    def get_segment_block_size(self):
        return self._segment_block_size

    def get_num_avl_node_blocks(self):
        return self._ll_sim.get_num_avl_node_blocks()

    def get_num_coalescence_record_blocks(self):
        return self._ll_sim.get_num_coalescence_record_blocks()

    def get_num_node_mapping_blocks(self):
        return self._ll_sim.get_num_node_mapping_blocks()

    def get_num_segment_blocks(self):
        return self._ll_sim.get_num_segment_blocks()

    def get_num_common_ancestor_events(self):
        return self._ll_sim.get_num_common_ancestor_events()

    def get_num_rejected_common_ancestor_events(self):
        return self._ll_sim.get_num_rejected_common_ancestor_events()

    def get_num_recombination_events(self):
        return self._ll_sim.get_num_recombination_events()

    def get_num_populations(self):
        return len(self._population_configurations)

    def get_num_migration_events(self):
        N = self.get_num_populations()
        matrix = [[0 for j in range(N)] for k in range(N)]
        flat = self._ll_sim.get_num_migration_events()
        for j in range(N):
            for k in range(N):
                matrix[j][k] = flat[j * N + k]
        return matrix

    def get_total_num_migration_events(self):
        return sum(self._ll_sim.get_num_migration_events())

    def get_num_multiple_recombination_events(self):
        return self._ll_sim.get_num_multiple_recombination_events()

    def get_max_memory(self):
        return self._ll_sim.get_max_memory()

    def get_configuration(self):
        return json.loads(self._ll_sim.get_configuration_json())

    def set_effective_population_size(self, effective_population_size):
        if effective_population_size <= 0:
            raise ValueError("Cannot set Ne to a non-positive value.")
        self._effective_population_size = effective_population_size

    def set_store_migrations(self, store_migrations):
        self._store_migrations = store_migrations

    def set_migration_matrix(self, migration_matrix):
        err = (
            "migration matrix must be a N x N square matrix encoded "
            "as a list-of-lists, where N is the number of populations "
            "defined in the population_configurations. The diagonal "
            "elements of this matrix must be zero. For example, a "
            "valid matrix for a 3 population system is "
            "[[0, 1, 1], [1, 0, 1], [1, 1, 0]]")
        if self._population_configurations is None:
            raise ValueError(
                "Cannot specify a migration matrix without also providing a "
                "population_configurations argument.")
        N = len(self._population_configurations)
        if not isinstance(migration_matrix, list):
            raise TypeError(err)
        if len(migration_matrix) != N:
            raise ValueError(err)
        for row in migration_matrix:
            if not isinstance(row, list):
                raise TypeError(err)
            if len(row) != N:
                raise ValueError(err)
        self._migration_matrix = migration_matrix

    def set_population_configurations(self, population_configurations):
        _check_population_configurations(population_configurations)
        self._population_configurations = population_configurations
        # Now set the default migration matrix.
        N = len(self._population_configurations)
        self._migration_matrix = [[0 for j in range(N)] for k in range(N)]

    def set_demographic_events(self, demographic_events):
        err = (
            "Demographic events must be a list of DemographicEvent instances "
            "sorted in non-decreasing order of time.")
        if not isinstance(demographic_events, collections.Iterable):
            raise TypeError(err)
        for event in demographic_events:
            if not isinstance(event, DemographicEvent):
                raise TypeError(err)
        self._demographic_events = demographic_events

    def set_model(self, model):
        """
        Sets the simulation model to the specified value. This may be either a string
        or a SimulationModel instance. If None, the default simulation model is used
        (i.e., Hudson's algorithm).
        """
        model_map = {
            "hudson": StandardCoalescent(),
            "smc": SmcApproxCoalescent(),
            "smc_prime": SmcPrimeApproxCoalescent()
        }
        if model is None:
            model_instance = StandardCoalescent()
        elif isinstance(model, str):
            lower_model = model.lower()
            if lower_model not in model_map:
                raise ValueError("Model '{}' unknown. Choose from {}".format(
                    model, list(model_map.keys())))
            model_instance = model_map[lower_model]
        else:
            if not isinstance(model, SimulationModel):
                raise TypeError(
                    "Simulation model must be a string or an instance of "
                    "SimulationModel")
            model_instance = model
        self._model = model_instance

    def set_segment_block_size(self, segment_block_size):
        self._segment_block_size = segment_block_size

    def set_avl_node_block_size(self, avl_node_block_size):
        self._avl_node_block_size = avl_node_block_size

    def set_node_mapping_block_size(self, node_mapping_block_size):
        self._node_mapping_block_size = node_mapping_block_size

    def set_coalescence_record_block_size(self, coalescence_record_block_size):
        self._coalescence_record_block_size = coalescence_record_block_size

    def create_ll_instance(self):
        # Now, convert the high-level values into their low-level
        # counterparts.
        d = len(self._population_configurations)
        Ne = self.get_effective_population_size()
        # The migration matrix must be flattened.
        scaled_migration_matrix = self.get_scaled_migration_matrix()
        ll_migration_matrix = [0 for j in range(d**2)]
        for j in range(d):
            for k in range(d):
                ll_migration_matrix[j * d + k] = scaled_migration_matrix[j][k]
        ll_population_configuration = [
            conf.get_ll_representation(Ne)
            for conf in self._population_configurations]
        ll_demographic_events = [
            event.get_ll_representation(d, Ne)
            for event in self._demographic_events]
        ll_simulation_model = self._model.get_ll_representation()
        ll_recombination_rate = self.get_per_locus_scaled_recombination_rate()
        ll_samples = [(pop, time / (4 * Ne)) for pop, time in self._samples]
        ll_sim = _msprime.Simulator(
            samples=ll_samples,
            random_generator=self._random_generator,
            num_loci=self._recombination_map.get_num_loci(),
            migration_matrix=ll_migration_matrix,
            population_configuration=ll_population_configuration,
            demographic_events=ll_demographic_events,
            model=ll_simulation_model,
            store_migrations=self._store_migrations,
            scaled_recombination_rate=ll_recombination_rate,
            max_memory=self._max_memory,
            segment_block_size=self._segment_block_size,
            avl_node_block_size=self._avl_node_block_size,
            node_mapping_block_size=self._node_mapping_block_size,
            coalescence_record_block_size=self._coalescence_record_block_size,
            migration_block_size=self._migration_block_size)
        return ll_sim

    def run(self):
        """
        Runs the simulation until complete coalescence has occurred.
        """
        if self._random_generator is None:
            raise ValueError("A random generator instance must be set")
        if self._ll_sim is None:
            self._ll_sim = self.create_ll_instance()
        self._ll_sim.run()

    def get_tree_sequence(self, mutation_generator=None, provenance_strings=[]):
        """
        Returns a TreeSequence representing the state of the simulation.
        """
        ll_recomb_map = self._recombination_map.get_ll_recombination_map()
        Ne = self.get_effective_population_size()
        self._ll_sim.populate_tables(
            self.node_table, self.edgeset_table, self.migration_table,
            Ne=Ne, recombination_map=ll_recomb_map)
        if mutation_generator is not None:
            mutation_generator.generate(
                self.node_table, self.edgeset_table, self.mutation_type_table,
                self.mutation_table)
        ll_tree_sequence = _msprime.TreeSequence()
        ll_tree_sequence.load_tables(
            self.node_table, self.edgeset_table, self.migration_table,
            self.mutation_type_table, self.mutation_table,
            provenance_strings=provenance_strings)
        return TreeSequence(ll_tree_sequence)

    def reset(self):
        """
        Resets the simulation so that we can perform another replicate.
        """
        if self._ll_sim is not None:
            self._ll_sim.reset()


class TreeSequence(object):
    """
    A TreeSequence represents the information generated in a coalescent
    simulation. This includes all the trees across the simulated region,
    along with the mutations (if any are present).
    """

    def __init__(self, ll_tree_sequence):
        self._ll_tree_sequence = ll_tree_sequence

    @property
    def ll_tree_sequence(self):
        return self.get_ll_tree_sequence()

    def get_ll_tree_sequence(self):
        return self._ll_tree_sequence

    @classmethod
    def load(cls, path):
        ts = _msprime.TreeSequence()
        ts.load(path)
        return TreeSequence(ts)

    @classmethod
    def load_tables(cls, **kwargs):
        ts = _msprime.TreeSequence()
        ts.load_tables(**kwargs)
        return TreeSequence(ts)

    def copy(self, sites=None):
        # Experimental API. Return a copy of this tree sequence, optionally with
        # the sites set to the specified list.
        node_table = msprime.NodeTable()
        edgeset_table = msprime.EdgesetTable()
        migration_table = msprime.MigrationTable()
        site_table = msprime.SiteTable()
        mutation_table = msprime.MutationTable()
        self._ll_tree_sequence.dump_tables(
            nodes=node_table, edgesets=edgeset_table, migrations=migration_table,
            sites=site_table, mutations=mutation_table)
        if sites is not None:
            site_table.reset()
            mutation_table.reset()
            for j, site in enumerate(sites):
                site_table.add_row(site.position, site.ancestral_state)
                for mutation in site.mutations:
                    mutation_table.add_row(j, mutation.node, mutation.derived_state)
        new_ll_ts = _msprime.TreeSequence()
        new_ll_ts.load_tables(
            nodes=node_table, edgesets=edgeset_table, migrations=migration_table,
            sites=site_table, mutations=mutation_table)
        return TreeSequence(new_ll_ts)

    @property
    def provenance(self):
        return self.get_provenance()

    def get_provenance(self):
        return self._ll_tree_sequence.get_provenance_strings()

    def newick_trees(self, precision=3, breakpoints=None, Ne=1):
        # TODO document this method.
        iterator = _msprime.NewickConverter(
            self._ll_tree_sequence, precision, Ne)
        if breakpoints is None:
            for length, tree in iterator:
                yield length, tree
        else:
            trees_covered = 0
            j = 0
            # TODO this is ugly. Update the alg so we don't need this
            # bracketing.
            bp = [0] + breakpoints + [self.get_sequence_length()]
            for length, tree in iterator:
                trees_covered += length
                while bp[j] < trees_covered:
                    j += 1
                    yield bp[j] - bp[j - 1], tree

    def dump(self, path, zlib_compression=False):
        """
        Writes the tree sequence to the specified file path.

        :param str path: The file path to write the TreeSequence to.
        :param bool zlib_compression: If True, use HDF5's native
            compression when storing the data leading to smaller
            file size. When loading, data will be decompressed
            transparently, but load times will be significantly slower.
        """
        self._ll_tree_sequence.dump(path, zlib_compression)

    def dump_tables(
            self, nodes=None, edgesets=None, migrations=None, sites=None,
            mutations=None):
        """
        Copy the contents of the tables underlying the tree sequence to the
        specified objects.

        :param NodeTable nodes: The NodeTable to load the nodes into.
        :param EdgesetTable edgesets: The EdgesetTable to load the edgesets into.
        :param MigrationTable migrations: The MigrationTable to load the migrations into.
        :param SiteTable sites: The SiteTable to load the sites into.
        :param MutationTable mutations: The NodeTable to load the mutations into.

        :return: A TableTuple containing all tables underlying the tree sequence.
        :rtype: TableTuple
        """
        # TODO document this and test the semantics to passing in new tables
        # as well as returning the updated tables.
        if nodes is None:
            nodes = NodeTable()
        if edgesets is None:
            edgesets = EdgesetTable()
        if migrations is None:
            migrations = MigrationTable()
        if sites is None:
            sites = SiteTable()
        if mutations is None:
            mutations = MutationTable()
        self._ll_tree_sequence.dump_tables(
            nodes=nodes, edgesets=edgesets, migrations=migrations, sites=sites,
            mutations=mutations)
        return TableTuple(
            nodes=nodes, edgesets=edgesets, migrations=migrations, sites=sites,
            mutations=mutations)

    def dump_text(
            self, nodes=None, edgesets=None, sites=None, mutations=None, precision=6):
        """
        Writes a text representation of the tables underlying the tree sequence
        to the specified connections.

        :param stream nodes: The file-like object (having a .write() method) to write
            the NodeTable to.
        :param stream edgesets: The file-like object to write the EdgesetTable to.
        :param stream sites: The file-like object to write the SiteTable to.
        :param stream mutations: The file-like object to write the MutationTable to.
        :param int precision: The number of digits of precision.
        """

        if nodes is not None:
            print("is_sample", "time", "population", sep="\t", file=nodes)
            for node in self.nodes():
                row = (
                    "{is_sample:d}\t"
                    "{time:.{precision}f}\t"
                    "{population:d}\t").format(
                        precision=precision, is_sample=node.is_sample(), time=node.time,
                        population=node.population)
                print(row, file=nodes)

        if edgesets is not None:
            print("left", "right", "parent", "children", sep="\t", file=edgesets)
            for edgeset in self.edgesets():
                children = ",".join(str(u) for u in edgeset.children)
                row = (
                    "{left:.{precision}f}\t"
                    "{right:.{precision}f}\t"
                    "{parent:d}\t"
                    "{children}").format(
                        precision=precision, left=edgeset.left, right=edgeset.right,
                        parent=edgeset.parent, children=children)
                print(row, file=edgesets)

        if sites is not None:
            print("position", "ancestral_state", sep="\t", file=sites)
            for site in self.sites():
                row = (
                    "{position:.{precision}f}\t"
                    "{ancestral_state}").format(
                        precision=precision, position=site.position,
                        ancestral_state=site.ancestral_state)
                print(row, file=sites)

        if mutations is not None:
            print("site", "node", "derived_state", sep="\t", file=mutations)
            for site in self.sites():
                for mutation in site.mutations:
                    row = (
                        "{site}\t"
                        "{node}\t"
                        "{derived_state}").format(
                            site=mutation.site, node=mutation.node,
                            derived_state=mutation.derived_state)
                    print(row, file=mutations)

    def dump_samples_text(self, samples, precision=6):
        """
        Writes a text representation of the entries in the NodeTable
        corresponding to samples to the specified connections.

        :param stream samples: The file-like object to write the subset of the NodeTable
            describing the samples to, with an extra column, `id`.
        :param int precision: The number of digits of precision.
        """

        print("id", "is_sample", "time", "population", sep="\t", file=samples)
        for node_id in self.samples():
            node = self.node(node_id)
            row = (
                "{node_id:d}\t"
                "{is_sample:d}\t"
                "{time:.{precision}f}\t"
                "{population:d}").format(
                    precision=precision, is_sample=node.is_sample(), time=node.time,
                    population=node.population, node_id=node_id)
            print(row, file=samples)

    @property
    def sample_size(self):
        return self.get_sample_size()

    def get_sample_size(self):
        """
        Returns the sample size for this tree sequence. This is the number
        of sample nodes in each tree.

        :return: The number of sample nodes in the tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_sample_size()

    @property
    def sequence_length(self):
        return self.get_sequence_length()

    def get_sequence_length(self):
        """
        Returns the sequence length in this tree sequence. This defines the
        genomic scale over which tree coordinates are defined. Given a
        tree sequence with a sequence length :math:`L`, the constituent
        trees will be defined over the half-closed interval
        :math:`(0, L]`. Each tree then covers some subset of this
        interval --- see :meth:`msprime.SparseTree.get_interval` for details.

        :return: The length of the sequence in this tree sequence in bases.
        :rtype: float
        """
        return self._ll_tree_sequence.get_sequence_length()

    @property
    def num_edgesets(self):
        return self._ll_tree_sequence.get_num_edgesets()

    @property
    def num_records(self):
        return self.get_num_records()

    # TODO deprecate
    def get_num_records(self):
        """
        Returns the number of coalescence records in this tree sequence.
        See the :meth:`.records` method for details on these objects.

        :return: The number of coalescence records defining this tree
            sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_edgesets()

    @property
    def num_trees(self):
        return self.get_num_trees()

    def get_num_trees(self):
        """
        Returns the number of distinct trees in this tree sequence. This
        is equal to the number of trees returned by the :meth:`.trees`
        method.

        :return: The number of trees in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_trees()

    @property
    def num_sites(self):
        return self.get_num_sites()

    def get_num_sites(self):
        return self._ll_tree_sequence.get_num_sites()

    @property
    def num_mutations(self):
        return self.get_num_mutations()

    def get_num_mutations(self):
        """
        Returns the number of mutations in this tree sequence. See
        the :meth:`msprime.TreeSequence.mutations` method for details on how
        mutations are defined.

        :return: The number of mutations in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_mutations()

    @property
    def num_nodes(self):
        return self.get_num_nodes()

    def get_num_nodes(self):
        """
        Returns the number of nodes in this tree sequence. This 1 + the
        largest value :math:`u` such that `u` is a node in any of the
        constituent trees.

        :return: The total number of nodes in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_nodes()

    # TODO deprecate
    def records(self):
        """
        Returns an iterator over the coalescence records in this tree
        sequence in time-sorted order. Each record is a tuple
        :math:`(l, r, u, c, t, d)` defining the assignment of a tree node
        across an interval. The range of this record is the half-open
        genomic interval :math:`[l, r)`, such that it applies to all
        positions :math:`l \leq x < r`. Each record represents the
        assignment of a pair of children :math:`c` to a parent
        parent :math:`u`. This assignment happens at :math:`t` generations
        in the past within the population with ID :math:`d`. If population
        information was not stored for this tree sequence then the
        population ID will be :const:`.NULL_POPULATION`.

        Each record returned is an instance of :func:`collections.namedtuple`,
        and may be accessed via the attributes ``left``, ``right``, ``node``,
        ``children``, ``time`` and ``population``, as well as the usual
        positional approach. For example, if we wished to print out the genomic
        length of each record, we could write::

        >>> for record in tree_sequence.records():
        >>>     print(record.right - record.left)

        :return: An iterator of all :math:`(l, r, u, c, t, d)` tuples defining
            the coalescence records in this tree sequence.
        :rtype: iter
        """
        for j in range(self.get_num_records()):
            yield CoalescenceRecord(*self._ll_tree_sequence.get_record(j))

    def migrations(self):
        for j in range(self._ll_tree_sequence.get_num_migrations()):
            yield Migration(*self._ll_tree_sequence.get_migration(j))

    def nodes(self):
        for j in range(self.num_nodes):
            yield self.node(j)

    def edgesets(self):
        for j in range(self.num_edgesets):
            left, right, parent, children = self._ll_tree_sequence.get_edgeset(j)
            yield Edgeset(parent=parent, children=children, left=left, right=right)

    def diffs(self):
        """
        Returns an iterator over the differences between adjacent trees in this
        tree sequence. Each diff returned by this method is a tuple of the form
        `(length, records_out, records_in)`. The `length` is the length of the
        genomic interval covered by the current tree, and is equivalent to the
        value returned by :meth:`msprime.SparseTree.get_length`. The
        `records_out` value is list of :math:`(u, c, t)` tuples, and
        corresponds to the coalescence records that have been invalidated by
        moving to the current tree.  As in the :meth:`.records` method,
        :math:`u` is a tree node, :math:`c` is a tuple containing its children,
        and :math:`t` is the time the event occurred.  These records are
        returned in time-decreasing order, such that the record affecting the
        highest parts of the tree (i.e., closest to the root) are returned
        first.  The `records_in` value is also a list of :math:`(u, c, t)`
        tuples, and these describe the records that must be applied to create
        the tree covering the current interval. These records are returned in
        time-increasing order, such that the records affecting the lowest parts
        of the tree (i.e., closest to the present day) are returned first.

        :return: An iterator over the diffs between adjacent trees in this
            tree sequence.
        :rtype: iter
        """
        return _msprime.TreeDiffIterator(self._ll_tree_sequence)

    def sites(self):
        for j in range(self.num_sites):
            pos, ancestral_state, mutations, index = self._ll_tree_sequence.get_site(j)
            yield Site(
                position=pos, ancestral_state=ancestral_state, index=index,
                mutations=[Mutation(*mutation) for mutation in mutations])

    def mutations(self):
        """
        Returns an iterator over the mutations in this tree sequence. Each
        mutation is represented as a tuple :math:`(x, u, j)` where :math:`x`
        is the position of the mutation in the sequence in chromosome
        coordinates, :math:`u` is the node over which the mutation
        occurred and :math:`j` is the zero-based index of the mutation within
        the overall tree sequence. Mutations are returned in non-decreasing
        order of position and increasing index.

        Each mutation returned is an instance of
        :func:`collections.namedtuple`, and may be accessed via the attributes
        ``position``, ``node`` and ``index`` as well as the usual positional
        approach. This is the recommended interface for working with mutations
        as it is both more readable and also ensures that code is forward
        compatible with future extensions.

        :return: An iterator of all :math:`(x, u, j)` tuples defining
            the mutations in this tree sequence.
        :rtype: iter
        """
        # TODO deprecate
        for site in self.sites():
            for mutation in site.mutations:
                yield DeprecatedMutation(
                    position=site.position, node=mutation.node, index=site.index)

    def breakpoints(self):
        """
        Returns an iterator over the breakpoints along the chromosome,
        including the two extreme points 0 and L. This is equivalent to

        >>> [0] + [t.get_interval()[1] for t in self.trees()]

        although we do not build an explicit list.

        :return: An iterator over all the breakpoints along the simulated
            sequence.
        :rtype: iter
        """
        yield 0
        for t in self.trees():
            yield t.get_interval()[1]

    def trees(
            self, tracked_samples=None, sample_counts=True, sample_lists=False,
            tracked_leaves=None, leaf_counts=None, leaf_lists=None):
        """
        Returns an iterator over the trees in this tree sequence. Each value
        returned in this iterator is an instance of
        :class:`.SparseTree`.

        The ``sample_counts`` and ``sample_lists`` parameters control the
        features that are enabled for the resulting trees. If ``sample_counts``
        is True, then it is possible to count the number of samples underneath
        a particular node in constant time using the :meth:`.get_num_samples`
        method. If ``sample_lists`` is True a more efficient algorithm is
        used in the :meth:`.SparseTree.samples` method.

        The ``tracked_samples`` parameter can be used to efficiently count the
        number of samples in a given set that exist in a particular subtree
        using the :meth:`.SparseTree.get_num_tracked_samples` method. It is an
        error to use the ``tracked_samples`` parameter when the ``sample_counts``
        flag is False.

        :warning: Do not store the results of this iterator in a list!
           For performance reasons, the same underlying object is used
           for every tree returned which will most likely lead to unexpected
           behaviour.

        :param list tracked_samples: The list of samples to be tracked and
            counted using the :meth:`.SparseTree.get_num_tracked_samples`
            method.
        :param bool sample_counts: If True, support constant time sample counts
            via the :meth:`.SparseTree.get_num_samples` and
            :meth:`.SparseTree.get_num_tracked_samples` methods.
        :param bool sample_lists: If True, provide more efficient access
            to the samples beneath a give node using the
            :meth:`.SparseTree.samples` method.
        :return: An iterator over the sparse trees in this tree sequence.
        :rtype: iter
        """
        # tracked_leaves, leaf_counts and leaf_lists are deprecated aliases
        # for tracked_samples, sample_counts and sample_lists respectively.
        # These are left over from an older version of the API when leaves
        # and samples were synonymous.
        if tracked_leaves is not None:
            tracked_samples = tracked_leaves
        if leaf_counts is not None:
            sample_counts = leaf_counts
        if leaf_lists is not None:
            sample_lists = leaf_lists
        flags = 0
        if sample_counts:
            flags |= _msprime.SAMPLE_COUNTS
        elif tracked_samples is not None:
            raise ValueError("Cannot set tracked_samples without sample_counts")
        if sample_lists:
            flags |= _msprime.SAMPLE_LISTS
        kwargs = {"flags": flags}
        if tracked_samples is not None:
            kwargs["tracked_samples"] = tracked_samples
        ll_sparse_tree = _msprime.SparseTree(self._ll_tree_sequence, **kwargs)
        iterator = _msprime.SparseTreeIterator(ll_sparse_tree)
        sparse_tree = SparseTree(ll_sparse_tree)
        for _ in iterator:
            yield sparse_tree

    def haplotypes(self):
        """
        Returns an iterator over the haplotypes resulting from the trees
        and mutations in this tree sequence as a string of '1's and '0's.
        The iterator returns a total of :math:`n` strings, each of which
        contains :math:`s` characters (:math:`n` is the sample size
        returned by :meth:`msprime.TreeSequence.get_sample_size` and
        :math:`s` is the number of mutations returned by
        :meth:`msprime.TreeSequence.get_num_mutations`). The first
        string returned is the haplotype for sample `0`, and so on.

        :return: An iterator over the haplotype strings for the samples in
            this tree sequence.
        :rtype: iter
        """
        return HaplotypeGenerator(self).haplotypes()

    def variants(self, as_bytes=False):
        """
        Returns an iterator over the variants in this tree sequence. Each
        variant corresponds to a single mutation and is represented as a tuple
        :math:`(x, u, j, g)`. The values of :math:`x`, :math:`u` and :math:`j`
        are identical to the values returned by the
        :meth:`.TreeSequence.mutations` method, and :math:`g` represents the
        sample genotypes for this variant. Thus, :math:`g[k]` is the observed
        state for sample :math:`k` at this site; zero represents the
        ancestral type and one the derived type.

        Each variant returned is an instance of :func:`collections.namedtuple`,
        and may be accessed via the attributes ``position``, ``node``,
        ``index`` and ``genotypes`` as well as the usual positional approach.
        This is the recommended interface for working with variants as it is
        both more readable and also ensures that code is forward compatible
        with future extensions.

        The returned genotypes may be either a numpy array of 1 byte unsigned
        integer 0/1 values, or a Python bytes object of '0'/'1' ASCII
        characters. This behaviour is controller by the ``as_bytes`` parameter.
        The default behaviour is to return a numpy array, which is
        substantially more efficient.

        :warning: The same numpy array is used to represent genotypes between
            iterations, so if you wish the store the results of this
            iterator you **must** take a copy of the array. This warning
            does not apply when ``as_bytes`` is True, as a new bytes object
            is allocated for each variant.

        :param bool as_bytes: If True, the genotype values will be returned
            as a Python bytes object. This is useful in certain situations
            (i.e., directly printing the genotypes) or when numpy is
            not available. Otherwise, genotypes are returned as a numpy
            array (the default).
        :return: An iterator of all :math:`(x, u, j, g)` tuples defining
            the variants in this tree sequence.
        """
        # TODO finalise API and documnent. See comments for the Variant type
        # for discussion on why the present form was chosen.
        n = self.get_sample_size()
        genotypes_buffer = bytearray(n)
        iterator = _msprime.VariantGenerator(
            self._ll_tree_sequence, genotypes_buffer, as_bytes)
        if as_bytes:
            for pos, ancestral_state, mutations, index in iterator:
                site = Site(
                    position=pos, ancestral_state=ancestral_state, index=index,
                    mutations=[Mutation(*mutation) for mutation in mutations])
                g = bytes(genotypes_buffer)
                yield Variant(position=pos, site=site, index=index, genotypes=g)
        else:
            check_numpy()
            g = np.frombuffer(genotypes_buffer, "u1", n)
            for pos, ancestral_state, mutations, index in iterator:
                site = Site(
                    position=pos, ancestral_state=ancestral_state, index=index,
                    mutations=[Mutation(*mutation) for mutation in mutations])
                yield Variant(position=pos, site=site, index=index, genotypes=g)

    def pairwise_diversity(self, samples=None):
        return self.get_pairwise_diversity(samples)

    def get_pairwise_diversity(self, samples=None):
        """
        Returns the value of pi, the pairwise nucleotide site diversity,
        which is the average number of mutations that differ between a randomly
        chosen pair of samples.  If `samples` is specified, calculate the
        diversity within this set.

        :param iterable samples: The set of samples within which we calculate
            the diversity. If None, calculate diversity within the entire
            sample.
        :return: The pairwise nucleotide site diversity.
        :rtype: float
        """
        if samples is None:
            samples = self.samples()
        else:
            samples = list(samples)
        return self._ll_tree_sequence.get_pairwise_diversity(samples)

    def mean_pairwise_tmrca(self, sample_sets, windows):
        """
        Finds the mean time to most recent common ancestor between pairs of samples
        as described in mean_pairwise_tmrca_matrix (which uses this function).
        Returns the upper triangle (including the diagonal) in row-major order,
        so if the output is `x`, then:

        >>> k=0
        >>> for w in range(len(windows)-1):
        >>>     for i in range(len(sample_sets)):
        >>>         for j in range(i,len(sample_sets)):
        >>>             trmca[i,j] = tmrca[j,i] = x[w][k]
        >>>             k += 1

        will fill out the matrix of mean TMRCAs in the `i`th window between (and
        within) each group of samples in `sample_sets` in the matrix `tmrca`.
        Alternatively, if `names` labels the sample_sets, the output labels are:

        >>> [".".join(names[i],names[j]) for i in range(len(names))
        >>>         for j in range(i,len(names))]

        :param list sample_sets: A list of sets of IDs of samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of the upper triangle of mean TMRCA values in row-major
            order, including the diagonal.
        """
        ns = len(sample_sets)
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i]*(n[j]-x[j]) + (n[i]-x[i])*x[j])
                    for i in range(ns) for j in range(i, ns)]

        out = self.branch_stats_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        # corrects the diagonal for self comparisons
        # and note factor of two for tree length -> real time
        for w in range(len(windows)-1):
            k = 0
            for i in range(ns):
                for j in range(i, ns):
                    if i == j:
                        if n[i] == 1:
                            out[w][k] = np.nan
                        else:
                            out[w][k] /= float(2 * n[i] * (n[i] - 1))
                    else:
                        out[w][k] /= float(2 * n[i] * n[j])
                    k += 1

        return out

    def mean_pairwise_tmrca_matrix(self, sample_sets, windows):
        """
        Finds the mean time to most recent common ancestor between pairs of
        samples from each set of samples and in each window. Returns a numpy
        array indexed by (window, sample_set, sample_set).  Diagonal entries are
        corrected so that the value gives the mean pairwise TMRCA for *distinct*
        samples, but it is not checked whether the sample_sets are disjoint
        (so offdiagonals are not corrected).  For this reason, if an element of
        `sample_sets` has only one element, the corresponding diagonal will be
        NaN.

        The mean TMRCA between two samples is defined to be one-half the length
        of all edges separating them in the tree at a uniformly chosen position
        on the genome.

        :param list sample_sets: A list of sets of IDs of samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of the upper triangle of mean TMRCA values in row-major
            order, including the diagonal.
        """
        x = self.mean_pairwise_tmrca(sample_sets, windows)
        ns = len(sample_sets)
        nw = len(windows) - 1
        A = np.ones((nw, ns, ns), dtype=float)
        for w in range(nw):
            k = 0
            for i in range(ns):
                for j in range(i, ns):
                    A[w, i, j] = A[w, j, i] = x[w][k]
                    k += 1
        return A

    def Y_vector(self, sample_sets, windows, indices):
        """
        Finds the 'Y' statistic between three sample_sets.  The sample_sets should
        be disjoint (the computation works fine, but if not the result depends
        on the amount of overlap).  If the sample_sets are A, B, and C, then the
        result gives the mean total length of any edge in the tree between a
        and the most recent common ancestor of b and c, where a, b, and c are
        random draws from A, B, and C respectively.

        The result is, for each window, a vector whose k-th entry is
            Y(sample_sets[indices[k][0]], sample_sets[indices[k][1]],
              sample_sets[indices[k][2]]).

        :param list sample_sets: A list of *three* sets of IDs of samples: (A,B,C).
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of triples of indices of sample_sets.
        :return: A list of numeric vectors of length equal to the length of
            indices, computed separately on each window.
        """
        for u in indices:
            if not len(u) == 3:
                raise ValueError("All indices should be of length 3.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i] * (n[j] - x[j]) * (n[k] - x[k])
                          + (n[i] - x[i]) * x[j] * x[k]) for i, j, k in indices]

        out = self.branch_stats_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        # corrects the diagonal for self comparisons
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                out[w][u] /= float(n[indices[u][0]] * n[indices[u][1]]
                                   * n[indices[u][2]])

        return out

    def Y(self, sample_sets, windows):
        """
        Finds the 'Y' statistic between the three groups of samples
        in sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).
        If the sample_sets are A, B, and C, then the result gives the mean total
        length of any edge in the tree between a and the most recent common
        ancestor of b and c, where a, b, and c are random draws from A, B, and
        C respectively.

        :param list sample_sets: A list of *three* sets of IDs of samples: (A,B,C).
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of numeric values computed separately on each window.
        """
        return self.Y_vector(sample_sets, windows, indices=[(0, 1, 2)])

    def f4_vector(self, sample_sets, windows, indices):
        """
        Finds the Patterson's f4 statistics between multiple subsets of four
        groups of sample_sets. The sample_sets should be disjoint (the computation
        works fine, but if not the result depends on the amount of overlap).

        :param list sample_sets: A list of four sets of IDs of samples: (A,B,C,D)
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of 4-tuples of indices of sample_sets.
        :return: A list of values of f4(A,B;C,D) of length equal to the length of
            indices, computed separately on each window.
        """
        for u in indices:
            if not len(u) == 4:
                raise ValueError("All tuples in indices should be of length 4.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float((x[i] * n[j] - x[j] * n[i]) * (x[k] * n[l] - x[l] * n[k]))
                    for i, j, k, l in indices]

        out = self.branch_stats_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        # corrects the diagonal for self comparisons
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                out[w][u] /= float(n[indices[u][0]] * n[indices[u][1]]
                                   * n[indices[u][2]] * n[indices[u][3]])

        return out

    def f4(self, sample_sets, windows):
        """
        Finds the Patterson's f4 statistics between the four groups of samples
        in sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).

        :param list sample_sets: A list of four sets of IDs of samples: (A,B,C,D)
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of values of f4(A,B;C,D) computed separately on each window.
        """
        if not len(sample_sets) == 4:
            raise ValueError("sample_sets should be of length 4.")
        return self.f4_vector(sample_sets, windows, indices=[(0, 1, 2, 3)])

    def f3_vector(self, sample_sets, windows, indices):
        """
        Finds the Patterson's f3 statistics between multiple subsets of three
        groups of samples in sample_sets. The sample_sets should be disjoint (the
        computation works fine, but if not the result depends on the amount of
        overlap).

        f3(A;B,C) is f4(A,B;A,C) corrected to not include self comparisons.

        If A does not contain at least three samples, the result is NaN.

        :param list sample_sets: A list of sets of IDs of samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of triples of indices of sample_sets.
        :return: A list of values of f3(A,B,C) computed separately on each window.
        """
        for u in indices:
            if not len(u) == 3:
                raise ValueError("All tuples in indices should be of length 3.")
        n = [len(x) for x in sample_sets]

        def f(x):
            return [float(x[i] * (x[i] - 1) * (n[j] - x[j]) * (n[k] - x[k])
                          + (n[i] - x[i]) * (n[i] - x[i] - 1) * x[j] * x[k]
                          - x[i] * (n[i] - x[i]) * (n[j] - x[j]) * x[k]
                          - (n[i] - x[i]) * x[i] * x[j] * (n[k] - x[k]))
                    for i, j, k in indices]

        out = self.branch_stats_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                if n[indices[u][0]] == 1:
                    out[w][u] = np.nan
                else:
                    out[w][u] /= float(n[indices[u][0]] * (n[indices[u][0]]-1)
                                       * n[indices[u][1]] * n[indices[u][2]])

        return out

    def f3(self, sample_sets, windows):
        """
        Finds the Patterson's f3 statistics between the three groups of samples
        in sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).

        f3(A;B,C) is f4(A,B;A,C) corrected to not include self comparisons.

        :param list sample_sets: A list of *three* sets of IDs of samples: (A,B,C),
            with the first set having at least two samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of values of f3(A,B,C) computed separately on each window.
        """
        if not len(sample_sets) == 3:
            raise ValueError("sample_sets should be of length 3.")
        return self.f3_vector(sample_sets, windows, indices=[(0, 1, 2)])

    def f2_vector(self, sample_sets, windows, indices):
        """
        Finds the Patterson's f2 statistics between multiple subsets of pairs
        of samples in sample_sets. The sample_sets should be disjoint (the
        computation works fine, but if not the result depends on the amount of
        overlap).

        f2(A;B) is f4(A,B;A,B) corrected to not include self comparisons.

        :param list sample_sets: A list of sets of IDs of samples, each having at
            least two samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :param list indices: A list of pairs of indices of sample_sets.
        :return: A list of values of f2(A,C) computed separately on each window.
        """
        for u in indices:
            if not len(u) == 2:
                raise ValueError("All tuples in indices should be of length 2.")
        n = [len(x) for x in sample_sets]
        for xlen in n:
            if not xlen > 1:
                raise ValueError("All sample_sets must have at least two samples.")

        def f(x):
            return [float(x[i] * (x[i] - 1) * (n[j] - x[j]) * (n[j] - x[j] - 1)
                          + (n[i] - x[i]) * (n[i] - x[i] - 1) * x[j] * (x[j] - 1)
                          - x[i] * (n[i] - x[i]) * (n[j] - x[j]) * x[j]
                          - (n[i] - x[i]) * x[i] * x[j] * (n[j] - x[j]))
                    for i, j in indices]

        out = self.branch_stats_vector(sample_sets, weight_fun=f, windows=windows)
        # move this division outside of f(x) so it only has to happen once
        for w in range(len(windows)-1):
            for u in range(len(indices)):
                out[w][u] /= float(n[indices[u][0]] * (n[indices[u][0]]-1)
                                   * n[indices[u][1]] * (n[indices[u][1]] - 1))

        return out

    def f2(self, sample_sets, windows):
        """
        Finds the Patterson's f2 statistics between the three groups of samples
        in sample_sets. The sample_sets should be disjoint (the computation works
        fine, but if not the result depends on the amount of overlap).

        f2(A;B) is f4(A,B;A,B) corrected to not include self comparisons.

        :param list sample_sets: A list of *two* sets of IDs of samples: (A,B),
            each having at least two samples.
        :param iterable windows: The breakpoints of the windows (including start
            and end, so has one more entry than number of windows).
        :return: A list of values of f2(A,B) computed separately on each window.
        """
        if not len(sample_sets) == 2:
            raise ValueError("sample_sets should be of length 2.")
        return self.f2_vector(sample_sets, windows, indices=[(0, 1)])

    def branch_stats(self, sample_sets, weight_fun):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a boolean.  A branch in a tree is weighted by weight_fun(x),
        where x[i] is the number of samples in sample_sets[i] below that
        branch.  This finds the sum of all counted branches for each tree,
        and averages this across the tree sequence, weighted by genomic length.
        '''
        out = self.branch_stats_vector(sample_sets, lambda x: [weight_fun(x)])
        assert len(out) == 1 and len(out[0]) == 1
        return out[0][0]

    def branch_stats_windowed(self, sample_sets, weight_fun, windows=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a boolean.  A branch in a tree is weighted by weight_fun(x),
        where x[i] is the number of samples in sample_sets[i] below that
        branch.  This finds the sum of all counted branches for each tree,
        and averages this across the tree sequence, weighted by genomic length.
        '''
        out = self.branch_stats_vector(sample_sets, lambda x: [weight_fun(x)], windows)
        assert len(out[0]) == 1
        return [x[0] for x in out]

    def branch_stats_vector(self, sample_sets, weight_fun, windows=None):
        '''
        Here sample_sets is a list of lists of samples, and weight_fun is a function
        whose argument is a list of integers of the same length as sample_sets
        that returns a boolean.  A branch in a tree is weighted by weight_fun(x),
        where x[i] is the number of samples in sample_sets[i] below that
        branch.  This finds the sum of all counted branches for each tree,
        and averages this across the tree sequence, weighted by genomic length.

        It does this separately for each window [windows[i], windows[i+1])
        and returns the values in a list.
        '''
        if windows is None:
            windows = (0, self.sequence_length)
        for U in sample_sets:
            if len(U) != len(set(U)):
                raise ValueError(
                    "elements of sample_sets cannot contain repeated elements.")
            for u in U:
                if not self.node(u).is_sample():
                    raise ValueError("Not all elements of sample_sets are samples.")
        num_windows = len(windows) - 1
        if windows[0] != 0.0:
            raise ValueError(
                "Windows must start at the start of the sequence (at 0.0).")
        if windows[-1] != self.sequence_length:
            raise ValueError("Windows must extend to the end of the sequence.")
        for k in range(num_windows):
            if windows[k + 1] <= windows[k]:
                raise ValueError("Windows must be increasing.")
        # initialize
        num_sample_sets = len(sample_sets)
        n_out = len(weight_fun([0 for a in range(num_sample_sets)]))
        S = [[0.0 for j in range(n_out)] for _ in range(num_windows)]
        L = [0.0 for j in range(n_out)]
        # print("sample_sets:", sample_sets)
        # print("n_out:",n_out)
        N = self.num_nodes
        X = [[int(u in a) for a in sample_sets] for u in range(N)]
        # we will essentially construct the tree
        pi = [-1 for j in range(N)]
        node_time = [0.0 for u in range(N)]
        # keep track of where we are for the windows
        chrom_pos = 0.0
        # index of *left-hand* end of the current window
        window_num = 0
        for length, records_out, records_in in self.diffs():
            for sign, records in ((-1, records_out), (+1, records_in)):
                for node, children, time in records:
                    # print("Record (",sign,"):",node,children,time)
                    # print("\t",X, "-->", L)
                    if sign == +1:
                        node_time[node] = time
                    dx = [0 for k in range(num_sample_sets)]
                    for child in children:
                        if sign == +1:
                            pi[child] = node
                        for k in range(num_sample_sets):
                            dx[k] += sign * X[child][k]
                        w = weight_fun(X[child])
                        dt = (node_time[pi[child]] - node_time[child])
                        for j in range(n_out):
                            L[j] += sign * dt * w[j]
                        # print("\t\tchild:",child,"+=",sign,"*",weight_fun(X[child]),
                        #    "*(",node_time[pi[child]],"-",node_time[child],")","-->",L)
                        if sign == -1:
                            pi[child] = -1
                    old_w = weight_fun(X[node])
                    for k in range(num_sample_sets):
                        X[node][k] += dx[k]
                    if pi[node] != -1:
                        w = weight_fun(X[node])
                        dt = (node_time[pi[node]] - node_time[node])
                        for j in range(n_out):
                            L[j] += dt * (w[j]-old_w[j])
                        # print("\t\tnode:",node,"+=",dt,"*(",weight_fun(X[node]),"-",
                        #   old_w,") -->",L)
                    # propagate change up the tree
                    u = pi[node]
                    if u != -1:
                        next_u = pi[u]
                        while u != -1:
                            old_w = weight_fun(X[u])
                            for k in range(num_sample_sets):
                                X[u][k] += dx[k]
                            # need to update X for the root,
                            # but the root does not have a branch length
                            if next_u != -1:
                                w = weight_fun(X[u])
                                dt = (node_time[pi[u]] - node_time[u])
                                for j in range(n_out):
                                    L[j] += dt*(w[j] - old_w[j])
                                # print("\t\tanc:",u,"+=",dt,"*(",weight_fun(X[u]),"-",
                                #    old_w,") -->",L)
                            u = next_u
                            next_u = pi[next_u]
                    # print("\t",X, "-->", L)
            # print("next tree:",L,length)
            while chrom_pos + length >= windows[window_num + 1]:
                # wrap up the last window
                this_length = windows[window_num + 1] - chrom_pos
                window_length = windows[window_num + 1] - windows[window_num]
                for j in range(n_out):
                    S[window_num][j] += L[j] * this_length
                    S[window_num][j] /= window_length
                length -= this_length
                # start the next
                if window_num < num_windows - 1:
                    window_num += 1
                    chrom_pos = windows[window_num]
                else:
                    # skips the else statement below
                    break
            else:
                for j in range(n_out):
                    S[window_num][j] += L[j] * length
                chrom_pos += length
        return S

    def node(self, u):
        flags, time, population, name = self._ll_tree_sequence.get_node(u)
        return Node(
            time=time, population=population, name=name,
            is_sample=flags & NODE_IS_SAMPLE)

    def time(self, u):
        return self.get_time(u)

    def get_time(self, u):
        """
        Returns the time that the specified ID was alive at.

        :param int u: The individual ID of interest.
        :return: The time at which the specified individual was alive at.
        :rtype: int
        """
        if u < 0 or u >= self.get_num_nodes():
            raise ValueError("ID out of bounds")
        node = self.node(u)
        return node.time

    def population(self, u):
        return self.get_population(u)

    def get_population(self, u):
        """
        Returns the population ID for the specified sample ID.

        :param int u: The individual  ID of interest.
        :return: The population ID where the specified individual lived.
            Returns :const:`.NULL_POPULATION` if no population information
            is available.
        :rtype: int
        """
        if u < 0 or u >= self.get_num_nodes():
            raise ValueError("ID out of bounds")
        node = self.node(u)
        return node.population

    def samples(self, population_id=None):
        return self.get_samples(population_id)

    def get_samples(self, population_id=None):
        """
        Returns the samples matching the specified population ID.

        :param int population_id: The population of interest. If None,
            return all samples.
        :return: The ID of the population we wish to find samples from.
            If None, return samples from all populations.
        :rtype: list
        """
        samples = self._ll_tree_sequence.get_samples()
        if population_id is not None:
            samples = [
                u for u in samples if self.get_population(u) == population_id]
        return samples

    def write_vcf(self, output, ploidy=1, contig_id="1"):
        """
        Writes a VCF formatted file to the specified file-like object. If a
        ploidy value is supplied, allele values are combined among adjacent
        samples to form a phased genotype of the required ploidy. For example,
        if we have a ploidy of 2 and a sample of size 6, then we will have
        3 diploid samples in the output, consisting of the combined alleles
        for samples [0, 1], [2, 3] and [4, 5]. If we had alleles 011110 at
        a particular variant, then we would output the genotypes 0|1, 1|1
        and 1|0 in VCF. Sample names are generated by appending the index
        to the prefix ``msp_`` such that we would have the sample names
        ``msp_0``, ``msp_1`` and ``msp_2`` in the running example.

        Example usage:

        >>> with open("output.vcf", "w") as vcf_file:
        >>>     tree_sequence.write_vcf(vcf_file, 2)

        :param File output: The file-like object to write the VCF output.
        :param int ploidy: The ploidy of the individual samples in the
            VCF. This sample size must be divisible by ploidy.
        :param str contig_id: The value of the CHROM column in the output VCF.
        """
        if ploidy < 1:
            raise ValueError("Ploidy must be >= sample size")
        if self.get_sample_size() % ploidy != 0:
            raise ValueError("Sample size must be divisible by ploidy")
        converter = _msprime.VcfConverter(
            self._ll_tree_sequence, ploidy=ploidy, contig_id=contig_id)
        output.write(converter.get_header())
        for record in converter:
            output.write(record)

    def simplify(self, samples=None, filter_root_mutations=True):
        if samples is None:
            samples = self.get_samples()
        ll_ts = _msprime.TreeSequence()
        self._ll_tree_sequence.simplify(ll_ts, samples, filter_root_mutations)
        new_ts = msprime.TreeSequence(ll_ts)
        # FIXME provenance
        # for provenance in self.get_provenance():
        #     new_ts.add_provenance(provenance)
        # parameters = {"TODO": "encode subset parameters"}
        # new_ts_provenance = get_provenance_dict("simplify", parameters)
        # new_ts.add_provenance(json.dumps(new_ts_provenance))
        return new_ts


class HaplotypeGenerator(object):

    def __init__(self, tree_sequence):
        self._tree_sequence = tree_sequence
        ts = self._tree_sequence.get_ll_tree_sequence()
        self._ll_haplotype_generator = _msprime.HaplotypeGenerator(ts)

    def get_haplotype(self, sample_id):
        return self._ll_haplotype_generator.get_haplotype(sample_id)

    def haplotypes(self):
        j = 0
        # Would use range here except for Python 2..
        while j < self._tree_sequence.get_sample_size():
            yield self.get_haplotype(j)
            j += 1


class RecombinationMap(object):
    """
    A RecombinationMap represents the changing rates of recombination
    along a chromosome. This is defined via two lists of numbers:
    ``positions`` and ``rates``, which must be of the same length.
    Given an index j in these lists, the rate of recombination
    per base per generation is ``rates[j]`` over the interval
    ``positions[j]`` to ``positions[j + 1]``. Consequently, the first
    position must be zero, and by convention the last rate value
    is also required to be zero (although it does not used).

    :param list positions: The positions (in bases) denoting the
        distinct intervals where recombination rates change. These can
        be floating point values.
    :param list rates: The list of rates corresponding to the supplied
        ``positions``. Recombination rates are specified per base,
        per generation.
    :param int num_loci: The maximum number of non-recombining loci
        in the underlying simulation. By default this is set to
        the largest possible value, allowing the maximum resolution
        in the recombination process. However, for a finite sites
        model this can be set to smaller values.
    """
    DEFAULT_NUM_LOCI = 2**32 - 1
    """
    The default number of non-recombining loci in a RecombinationMap.
    """
    def __init__(self, positions, rates, num_loci=None):
        m = self.DEFAULT_NUM_LOCI
        if num_loci is not None:
            m = num_loci
        self._ll_recombination_map = _msprime.RecombinationMap(
            m, positions, rates)

    @classmethod
    def uniform_map(cls, length, rate, num_loci=None):
        return cls([0, length], [rate, 0], num_loci)

    @classmethod
    def read_hapmap(cls, filename):
        """
        Parses the specified file in HapMap format. These files must contain
        a single header line (which is ignored), and then each subsequent
        line denotes a position/rate pair. Positions are in units of bases,
        and recombination rates in centimorgans/Megabase. The first column
        in this file is ignored, as are subsequence columns after the
        Position and Rate. A sample of this format is as follows::

            Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
            chr1	55550	2.981822	0.000000
            chr1	82571	2.082414	0.080572
            chr1	88169	2.081358	0.092229
            chr1	254996	3.354927	0.439456
            chr1	564598	2.887498	1.478148

        :param str filename: The name of the file to be parsed. This may be
            in plain text or gzipped plain text.
        """
        positions = []
        rates = []
        if filename.endswith(".gz"):
            f = gzip.open(filename)
        else:
            f = open(filename)
        try:
            # Skip the header line
            f.readline()
            for j, line in enumerate(f):
                pos, rate, = map(float, line.split()[1:3])
                if j == 0:
                    if pos != 0:
                        positions.append(0)
                        rates.append(0)
                positions.append(pos)
                # Rate is expressed in centimorgans per megabase, which
                # we convert to per-base rates
                rates.append(rate * 1e-8)
            assert rate == 0
        finally:
            f.close()
        return cls(positions, rates)

    def get_ll_recombination_map(self):
        return self._ll_recombination_map

    def physical_to_genetic(self, physical_x):
        return self._ll_recombination_map.physical_to_genetic(physical_x)

    def genetic_to_physical(self, genetic_x):
        return self._ll_recombination_map.genetic_to_physical(genetic_x)

    def get_total_recombination_rate(self):
        return self._ll_recombination_map.get_total_recombination_rate()

    def get_per_locus_recombination_rate(self):
        return self._ll_recombination_map.get_per_locus_recombination_rate()

    def get_size(self):
        return self._ll_recombination_map.get_size()

    def get_num_loci(self):
        return self._ll_recombination_map.get_num_loci()

    def get_positions(self):
        return self._ll_recombination_map.get_positions()

    def get_length(self):
        # this is a wasteful creation of a list; push the method down
        # into the low-level API.
        return self._ll_recombination_map.get_positions()[-1]

    def get_rates(self):
        return self._ll_recombination_map.get_rates()


class PopulationConfiguration(object):
    """
    The initial configuration of a population (or deme) in a simulation.

    :param int sample_size: The number of initial samples that are drawn
        from this population.
    :param float initial_size: The absolute size of the population at time
        zero. Defaults to the reference population size :math:`N_e`.
    :param float growth_rate: The exponential growth rate of the population
        per generation. Growth rates can be negative. This is zero for a
        constant population size. Defaults to 0.
    """
    def __init__(self, sample_size=None, initial_size=None, growth_rate=0.0):
        self.sample_size = sample_size
        self.initial_size = initial_size
        self.growth_rate = growth_rate

    def get_ll_representation(self, Ne):
        """
        Returns the low-level representation of this PopulationConfiguration.
        """
        initial_size = Ne if self.initial_size is None else self.initial_size
        return {
            "initial_size": initial_size / Ne,
            "growth_rate": self.growth_rate * 4 * Ne
        }


# TODO use this and its inverse throughout the code.
def generations_to_coalescent(time, Ne):
    return time / (4 * Ne)


class DemographicEvent(object):
    def __init__(self, type_, time):
        self.type = type_
        self.time = time

    def _get_scaled_time(self, Ne):
        return generations_to_coalescent(self.time, Ne)


class PopulationParametersChange(DemographicEvent):
    """
    Changes the demographic parameters of a population at a given time.

    This event generalises the ``-eg``, ``-eG``, ``-en`` and ``-eN``
    options from ``ms``. Note that unlike ``ms`` we do not automatically
    set growth rates to zero when the population size is changed.

    :param float time: The time at which this event occurs in generations.
    :param float initial_size: The absolute size of the population
        at the beginning of the time slice starting at ``time``. If None,
        this is calculated according to the initial population size and
        growth rate over the preceding time slice.
    :param float growth_rate: The new per-generation growth rate. If None,
        the growth rate is not changed. Defaults to None.
    :param int population_id: The ID of the population affected. If
        ``population_id`` is None, the changes affect all populations
        simultaneously.
    """
    def __init__(
            self, time, initial_size=None, growth_rate=None,
            population_id=None):
        super(PopulationParametersChange, self).__init__(
            "population_parameters_change", time)
        if growth_rate is None and initial_size is None:
            raise ValueError(
                "Must specify one or more of growth_rate and initial_size")
        self.time = time
        self.growth_rate = growth_rate
        self.initial_size = initial_size
        self.population_id = -1 if population_id is None else population_id

    def get_ll_representation(self, num_populations, Ne):
        ret = {
            "type": self.type,
            "time": self._get_scaled_time(Ne),
            "population_id": self.population_id
        }
        if self.growth_rate is not None:
            ret["growth_rate"] = self.growth_rate * 4 * Ne
        if self.initial_size is not None:
            ret["initial_size"] = self.initial_size / Ne
        return ret

    def __str__(self):
        s = "Population parameter change for {}: ".format(self.population_id)
        if self.initial_size is not None:
            s += "initial_size -> {} ".format(self.initial_size)
        if self.growth_rate is not None:
            s += "growth_rate -> {} ".format(self.growth_rate)
        return s


class MigrationRateChange(DemographicEvent):
    """
    Changes the rate of migration to a new value at a specific time.

    :param float time: The time at which this event occurs in generations.
    :param float rate: The new per-generation migration rate.
    :param tuple matrix_index: A tuple of two population IDs descibing
        the matrix index of interest. If ``matrix_index`` is None, all
        non-diagonal entries of the migration matrix are changed
        simultaneously.
    """
    def __init__(self, time, rate, matrix_index=None):
        super(MigrationRateChange, self).__init__(
            "migration_rate_change", time)
        self.rate = rate
        self.matrix_index = matrix_index

    def get_ll_representation(self, num_populations, Ne):
        matrix_index = -1
        scaled_rate = 4 * Ne * self.rate
        if self.matrix_index is not None:
            matrix_index = (
                self.matrix_index[0] * num_populations + self.matrix_index[1])
        return {
            "type": self.type,
            "time": self._get_scaled_time(Ne),
            "migration_rate": scaled_rate,
            "matrix_index": matrix_index
        }

    def __str__(self):
        if self.matrix_index is None:
            ret = "Migration rate change to {} everywhere".format(self.rate)
        else:
            ret = "Migration rate change for {} to {}".format(
                self.matrix_index, self.rate)
        return ret


class MassMigration(DemographicEvent):
    """
    A mass migration event in which some fraction of the population in
    one deme simultaneously move to another deme, viewed backwards in
    time. For each lineage currently present in the source population,
    they move to the destination population with probability equal to
    ``proportion``.

    This event class generalises the population split (``-ej``) and
    admixture (``-es``) events from ``ms``. Note that MassMigrations
    do *not* have any side effects on the migration matrix.

    :param float time: The time at which this event occurs in generations.
    :param int source: The ID of the source population.
    :param int destination: The ID of the destination population.
    :param float proportion: The probability that any given lineage within
        the source population migrates to the destination population.
    """
    def __init__(self, time, source, destination, proportion=1.0):
        super(MassMigration, self).__init__("mass_migration", time)
        self.source = source
        self.destination = destination
        self.proportion = proportion

    def get_ll_representation(self, num_populations, Ne):
        return {
            "type": self.type,
            "time": self._get_scaled_time(Ne),
            "source": self.source,
            "destination": self.destination,
            "proportion": self.proportion
        }

    def __str__(self):
        return (
            "Mass migration: lineages move from {} to {} with "
            "probability {}".format(
                self.source, self.destination, self.proportion))


class SimpleBottleneck(DemographicEvent):
    # This is an unsupported/undocumented demographic event.
    def __init__(self, time, population_id=0, proportion=1.0):
        super(SimpleBottleneck, self).__init__("simple_bottleneck", time)
        self.population_id = population_id
        self.proportion = proportion

    def get_ll_representation(self, num_populations, Ne):
        return {
            "type": self.type,
            "time": self._get_scaled_time(Ne),
            "population_id": self.population_id,
            "proportion": self.proportion
        }

    def __str__(self):
        return (
            "Simple bottleneck: lineages in population {} coalesce "
            "probability {}".format(self.population_id, self.proportion))


class InstantaneousBottleneck(DemographicEvent):
    # TODO document

    def __init__(self, time, population_id=0, strength=1.0):
        super(InstantaneousBottleneck, self).__init__(
            "instantaneous_bottleneck", time)
        self.population_id = population_id
        self.strength = strength

    def get_ll_representation(self, num_populations, Ne):
        return {
            "type": self.type,
            "time": self._get_scaled_time(Ne),
            "population_id": self.population_id,
            "strength": generations_to_coalescent(self.strength, Ne)
        }

    def __str__(self):
        return (
            "Instantaneous bottleneck in population {}: equivalent to {} "
            "generations of the coalescent".format(
                self.population_id, self.strength))


class SimulationModel(object):
    """
    Superclass of all simulation models.
    """
    name = None

    def get_ll_representation(self):
        return {"name": self.name}


class StandardCoalescent(SimulationModel):
    """
    The classical coalescent with recombination model (i.e., Hudson's algorithm).
    """
    name = "hudson"


class SmcApproxCoalescent(SimulationModel):
    # TODO document
    name = "smc"


class SmcPrimeApproxCoalescent(SimulationModel):
    # TODO document
    name = "smc_prime"


class ParametricSimulationModel(SimulationModel):
    """
    The superclass of simulation models that require parameters.
    """
    def get_ll_representation(self):
        d = super(ParametricSimulationModel, self).get_ll_representation()
        d.update(self.__dict__)
        return d


class BetaCoalescent(ParametricSimulationModel):
    # TODO document.
    name = "beta"

    # TODO what is a meaningful value for this parameter? Ideally, the default
    # would be the equivalent of the Kingman coalescent or something similar.
    def __init__(self, alpha=1, truncation_point=None):
        self.alpha = alpha
        if truncation_point is None:
            truncation_point = sys.float_info.max
        self.truncation_point = truncation_point


class DiracCoalescent(ParametricSimulationModel):
    # TODO document
    name = "dirac"

    # TODO What is a meaningful default for this value? See above.
    def __init__(self, psi=0.5, c=10.0):
        self.psi = psi
        self.c = c


class Population(object):
    """
    Simple class to represent the state of a population in terms of its
    demographic parameters. This is intended to be initialised from the
    corresponding low-level values so that they can be rescaled back into
    input units.
    """
    def __init__(
            self, Ne=None, sample_size=None, initial_size=None,
            growth_rate=None):
        self.Ne = Ne
        self.initial_size = initial_size * self.Ne
        self.growth_rate = growth_rate / (4 * Ne)

    def get_size(self, time):
        """
        Gets the size of the population after the specified amount of
        time.
        """
        size = self.initial_size
        if self.growth_rate != 0:
            size = self.initial_size * math.exp(-self.growth_rate * time)
        return size


class DemographyDebugger(object):
    """
    A class to facilitate debugging of population parameters and migration
    rates in the past.
    """
    def __init__(
            self, Ne=1, population_configurations=None, migration_matrix=None,
            demographic_events=[]):
        self._precision = 3
        self._Ne = Ne
        self._simulator = simulator_factory(
            Ne=Ne,
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        self._demographic_events = sorted(
            demographic_events, key=lambda e: e.time)

    def _scaled_time_to_generations(self, t):
        return 4 * self._Ne * t

    def _print_populations(
            self, start_time, end_time, migration_matrix, populations, output):
        field_width = self._precision + 6
        growth_rate_field_width = 14
        sep_str = " | "
        N = len(migration_matrix)
        fmt = (
            "{id:<2} "
            "{start_size:^{field_width}}"
            "{end_size:^{field_width}}"
            "{growth_rate:>{growth_rate_field_width}}")
        print(fmt.format(
            id="", start_size="start", end_size="end",
            growth_rate="growth_rate", field_width=field_width,
            growth_rate_field_width=growth_rate_field_width), end=sep_str,
            file=output)
        for k in range(N):
            print(
                "{0:^{1}}".format(k, field_width), end="", file=output)
        print(file=output)
        h = "-" * (field_width - 1)
        print(
            fmt.format(
                id="", start_size=h, end_size=h, growth_rate=h,
                field_width=field_width,
                growth_rate_field_width=growth_rate_field_width),
            end=sep_str, file=output)
        for k in range(N):
            s = "-" * (field_width - 1)
            print("{0:<{1}}".format(s, field_width), end="", file=output)
        print(file=output)
        for j, pop in enumerate(populations):
            s = (
                "{id:<2}|"
                "{start_size:^{field_width}.{precision}g}"
                "{end_size:^{field_width}.{precision}g}"
                "{growth_rate:>{growth_rate_field_width}.{precision}g}"
                ).format(
                    id=j, start_size=pop.initial_size,
                    end_size=pop.get_size(end_time - start_time),
                    growth_rate=pop.growth_rate,
                    precision=self._precision, field_width=field_width,
                    growth_rate_field_width=growth_rate_field_width)
            print(s, end=sep_str, file=output)
            for k in range(N):
                x = migration_matrix[j][k]
                print("{0:^{1}.{2}g}".format(
                    x, field_width, self._precision), end="",
                    file=output)
            print(file=output)

    def print_history(self, output=sys.stdout):
        """
        Prints a summary of the history of the populations.
        """
        abs_tol = 1e-9
        ll_sim = self._simulator.create_ll_instance()
        N = self._simulator.get_num_populations()
        Ne = self._Ne
        start_time = 0
        scaled_end_time = 0
        event_index = 0
        while not math.isinf(scaled_end_time):
            events = []
            while (
                    event_index < len(self._demographic_events) and
                    almost_equal(
                        self._demographic_events[event_index].time,
                        start_time, abs_tol=abs_tol)):
                events.append(self._demographic_events[event_index])
                event_index += 1
            if len(events) > 0:
                print(
                    "Events @ generation {}".format(start_time), file=output)
            for event in events:
                assert almost_equal(event.time, start_time, abs_tol=abs_tol)
                print("   -", event, file=output)
            print(file=output)
            scaled_end_time = ll_sim.debug_demography()
            end_time = self._scaled_time_to_generations(scaled_end_time)
            m = ll_sim.get_migration_matrix()
            migration_matrix = [
                [m[j * N + k] / (4 * Ne) for j in range(N)] for k in range(N)]
            populations = [
                Population(Ne=self._Ne, **d)
                for d in ll_sim.get_population_configuration()]
            s = "Epoch: {} -- {} generations".format(start_time, end_time)
            print("=" * len(s), file=output)
            print(s, file=output)
            print("=" * len(s), file=output)
            self._print_populations(
                start_time, end_time, migration_matrix, populations, output)
            print(file=output)
            start_time = end_time


def harmonic_number(n):
    """
    Returns the nth Harmonic number.
    """
    EulerGamma = 0.5772156649
    return math.log(n) + EulerGamma
