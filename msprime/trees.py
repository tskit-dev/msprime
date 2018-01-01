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
import itertools
import json
import sys

try:
    import numpy as np
    _numpy_imported = True
except ImportError:
    _numpy_imported = False


import _msprime
import msprime.drawing as drawing
import msprime.provenance as provenance
import msprime.tables as tables

from _msprime import NODE_IS_SAMPLE

NULL_NODE = -1

NULL_POPULATION = -1

NULL_MUTATION = -1

IS_PY2 = sys.version_info[0] < 3


def check_numpy():
    if not _numpy_imported:
        raise RuntimeError("numpy is required for this operation.")


CoalescenceRecord = collections.namedtuple(
    "CoalescenceRecord",
    ["left", "right", "node", "children", "time", "population"])


Migration = collections.namedtuple(
    "Migration",
    ["left", "right", "node", "source", "dest", "time"])


# TODO We need to get rid of the these namedtuples where possible and
# make proper classes where possible. Also, need to standardise on 'id'
# rather than index throughout.

Site = collections.namedtuple(
    "Site",
    ["position", "ancestral_state", "index", "mutations", "metadata"])


Mutation = collections.namedtuple(
    "Mutation",
    ["site", "node", "derived_state", "parent", "id", "metadata"])


# This is provided for backwards compatibility with the deprecated mutations()
# iterator.
DeprecatedMutation = collections.namedtuple(
    "DeprecatedMutation",
    ["position", "node", "index"])


# This form was chosen to try to break as little of existing code as possible.
# The position and index values are retained for compatability,
# although they are redundant.
Variant = collections.namedtuple(
    "Variant",
    ["position", "site", "index", "genotypes", "alleles"])


# TODO this interface is rubbish. Should have much better printing options.
class SimpleContainer(object):

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return repr(self.__dict__)


class Node(SimpleContainer):
    def __init__(
            self, id_=None, time=0, population=NULL_POPULATION, metadata="",
            is_sample=False):
        self.id = id_
        self.time = time
        self.population = population
        self.metadata = metadata
        self.flags = 0
        if is_sample:
            self.flags |= NODE_IS_SAMPLE

    def is_sample(self):
        return self.flags & NODE_IS_SAMPLE


class Edge(SimpleContainer):
    def __init__(self, left, right, parent, child):
        self.left = left
        self.right = right
        self.parent = parent
        self.child = child

    def __repr__(self):
        return "{{left={:.3f}, right={:.3f}, parent={}, child={}}}".format(
            self.left, self.right, self.parent, self.child)


class Edgeset(SimpleContainer):
    def __init__(self, left, right, parent, children):
        self.left = left
        self.right = right
        self.parent = parent
        self.children = children

    def __repr__(self):
        return "{{left={:.3f}, right={:.3f}, parent={}, children={}}}".format(
            self.left, self.right, self.parent, self.children)


class Provenance(SimpleContainer):
    def __init__(self, id_=None, timestamp=None, record=None):
        self.id = id_
        self.timestamp = timestamp
        self.record = record


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
        roots = self.roots
        return sum(
            self.get_branch_length(u) for u in self.nodes() if u not in roots)

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

    def left_child(self, u):
        return self._ll_sparse_tree.get_left_child(u)

    def right_child(self, u):
        return self._ll_sparse_tree.get_right_child(u)

    def left_sib(self, u):
        return self._ll_sparse_tree.get_left_sib(u)

    def right_sib(self, u):
        return self._ll_sparse_tree.get_right_sib(u)

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
        Returns the time of the specified node in generations.

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
        """
        Returns the number of nodes in the sparse tree.

        :rtype: int
        """
        return self._ll_sparse_tree.get_num_nodes()

    @property
    def num_roots(self):
        """
        The number of roots in this tree, as defined in the :attr:`.roots` attribute.

        Requires O(number of roots) time.

        :rtype: int
        """
        return self._ll_sparse_tree.get_num_roots()

    @property
    def roots(self):
        """
        The list of roots in this tree. A root is defined as a unique endpoint of
        the paths starting at samples. We can define the set of roots as follows:

        .. code-block:: python

            roots = set()
            for u in tree_sequence.samples():
                while tree.parent(u) != msprime.NULL_NODE:
                    u = tree.parent(u)
                roots.add(u)
            # roots is now the set of all roots in this tree.
            assert sorted(roots) == sorted(tree.roots)

        The roots of the tree are returned in a list, in no particular order.

        Requires O(number of roots) time.

        :rtype: list
        """
        roots = []
        u = self.left_root
        while u != NULL_NODE:
            roots.append(u)
            u = self.right_sib(u)
        return roots

    @property
    def root(self):
        """
        The root of this tree. If the tree contains multiple roots, a ValueError is
        raised indicating that the :attr:`.roots` attribute should be used instead.

        :return: The root node.
        :rtype: int
        :raises: ValueError if this tree contains more than one root.
        """
        root = self.left_root
        if root != NULL_NODE and self.right_sib(root) != NULL_NODE:
            raise ValueError("More than one root exists. Use tree.roots instead")
        return root

    def get_root(self):
        # Deprecated alias for self.root
        return self.root

    @property
    def left_root(self):
        return self._ll_sparse_tree.get_left_root()

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
        return self._ll_sparse_tree.get_left(), self._ll_sparse_tree.get_right()

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

    # The sample_size (or num_samples) is really a property of the tree sequence,
    # and so we should provide access to this via a tree.tree_sequence.num_samples
    # property access. However, we can't just remove the method as a lot of code
    # may depend on it. To complicate things a bit more, sample_size has been
    # changed to num_samples elsewhere for consistency. We can't do this here
    # because there is already a num_samples method which returns the number of
    # samples below a particular node. The best thing to do is probably to
    # undocument the sample_size property, but keep it around for ever.

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
            self, path=None, width=None, height=None, times=False,
            mutation_locations=True, mutation_labels=False,
            internal_node_labels=True, leaf_node_labels=True, show_times=None,
            node_label_text=None, format=None):
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
        :param map node_label_text: If specified, show custom labels for the nodes
            that are present in the map. Any nodes not specified in the map will
            have have their default labels.
        :param str format: The format of the returned image. Currently supported
            are 'svg', 'ascii' and 'unicode'.
        :param bool show_times: Deprecated alias for ``times``.
        :return: A representation of this tree in SVG format.
        :rtype: str
        """
        # show_times is a deprecated alias for times.
        if show_times is not None:
            times = show_times
        output = drawing.draw_tree(
            self, format=format, width=width, height=height, times=times,
            mutation_locations=mutation_locations, mutation_labels=mutation_labels,
            internal_node_labels=internal_node_labels, leaf_node_labels=leaf_node_labels,
            node_label_text=node_label_text)
        if path is not None:
            with open(path, "w") as f:
                f.write(output)
        return output

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
            pos, ancestral_state, mutations, index, metadata = ll_site
            yield Site(
                position=pos, ancestral_state=ancestral_state, index=index,
                mutations=[Mutation(*mutation) for mutation in mutations],
                metadata=metadata)

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
        roots = [u]
        if u is None:
            roots = self.roots
        for root in roots:
            for v in self.nodes(root):
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
        roots = [u]
        if u is None:
            roots = self.roots
        for root in roots:
            if self._ll_sparse_tree.get_flags() & _msprime.SAMPLE_LISTS:
                for v in _msprime.SampleListIterator(self._ll_sparse_tree, root):
                    yield v
            else:
                for v in self._sample_generator(root):
                    yield v

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
        roots = [u]
        if u is None:
            roots = self.roots
        return sum(self._ll_sparse_tree.get_num_samples(u) for u in roots)

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
        roots = [u]
        if u is None:
            roots = self.roots
        if not (self._ll_sparse_tree.get_flags() & _msprime.SAMPLE_COUNTS):
            raise RuntimeError(
                "The get_num_tracked_samples method is only supported "
                "when sample_counts=True.")
        return sum(self._ll_sparse_tree.get_num_tracked_samples(root) for root in roots)

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
        # TODO add a nonrecursive version of the inorder traversal.
        children = self.get_children(u)
        mid = len(children) // 2
        for c in children[:mid]:
            for v in self._inorder_traversal(c):
                yield v
        yield u
        for c in children[mid:]:
            for v in self._inorder_traversal(c):
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
        methods = {
            "preorder": self._preorder_traversal,
            "inorder": self._inorder_traversal,
            "postorder": self._postorder_traversal,
            "levelorder": self._levelorder_traversal,
            "breadthfirst": self._levelorder_traversal
        }
        try:
            iterator = methods[order]
        except KeyError:
            raise ValueError("Traversal ordering '{}' not supported".format(order))
        roots = [root]
        if root is None:
            roots = self.roots
        for u in roots:
            for v in iterator(u):
                yield v

    def newick(self, precision=14, time_scale=1):
        s = self._ll_sparse_tree.get_newick(precision=precision, time_scale=time_scale)
        if not IS_PY2:
            s = s.decode()
        return s

    @property
    def parent_dict(self):
        return self.get_parent_dict()

    def get_parent_dict(self):
        pi = {
            u: self.parent(u) for u in range(self.num_nodes)
            if self.parent(u) != NULL_NODE}
        return pi

    def __str__(self):
        return str(self.get_parent_dict())


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
    """
    Loads a tree sequence from the table objects provided.  If sequence_length
    is 0 or not provided, it is inferred to be equal to the largest right value
    in the edges.

    :param NodeTable nodes: The NodeTable.
    :param EdgeTable edges: The EdgeTable.
    :param float sequence_length: The length of the sequence.
    :param MigrationTable migrations: The MigrationTable.
    :param SiteTable sites: The SiteTable.
    :param MutationTable mutations: The MutationTable.
    :param ProvenanceTable provenances: The ProvenanceTable.
    :return: A :class:`msprime.TreeSequence` consistent with the tables.
    :rtype: TreeSequence
    """
    return TreeSequence.load_tables(*args, **kwargs)


def parse_nodes(source, sep=None):
    """
    Parse the specified file-like object and return a NodeTable instance.  The
    object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``is_sample``, ``time``, and
    optionally, ``population``.  Further requirements are described in
    :class:`NodeTable`.  Note that node ``id`` is not included, but implied by
    order in the file.
    """
    # Read the header and find the indexes of the required fields.
    table = tables.NodeTable()
    header = source.readline().strip("\n").split(sep)
    is_sample_index = header.index("is_sample")
    time_index = header.index("time")
    population_index = None
    metadata_index = None
    try:
        population_index = header.index("population")
    except ValueError:
        pass
    try:
        metadata_index = header.index("metadata")
    except ValueError:
        pass
    for line in source:
        tokens = line.split(sep)
        if len(tokens) >= 2:
            is_sample = int(tokens[is_sample_index])
            time = float(tokens[time_index])
            flags = 0
            if is_sample != 0:
                flags |= NODE_IS_SAMPLE
            population = NULL_POPULATION
            if population_index is not None:
                population = int(tokens[population_index])
            metadata = b''
            if metadata_index is not None and metadata_index < len(tokens):
                metadata = tables.text_decode_metadata(tokens[metadata_index])
            table.add_row(
                flags=flags, time=time, population=population, metadata=metadata)
    return table


def parse_edges(source, sep=None):
    """
    Parse the specified file-like object and return a EdgeTable instance.
    The object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``left``, ``right``, ``parent``,
    and ``child``. Several edges may specified at the same time using a single
    line by making the ``child`` field a comma-separated list. Further
    requirements are described in :class:`EdgeTable`.
    """
    table = tables.EdgeTable()
    header = source.readline().strip("\n").split(sep)
    left_index = header.index("left")
    right_index = header.index("right")
    parent_index = header.index("parent")
    children_index = header.index("child")
    table = tables.EdgeTable()
    for line in source:
        tokens = line.split(sep)
        if len(tokens) >= 4:
            left = float(tokens[left_index])
            right = float(tokens[right_index])
            parent = int(tokens[parent_index])
            children = tuple(map(int, tokens[children_index].split(",")))
            for child in children:
                table.add_row(left=left, right=right, parent=parent, child=child)
    return table


def parse_sites(source, sep=None):
    """
    Parse the specified file-like object and return a SiteTable instance.  The
    object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``position`` and
    ``ancestral_state``.  Further requirements are described in
    :class:`SiteTable`.
    """
    header = source.readline().strip("\n").split(sep)
    position_index = header.index("position")
    ancestral_state_index = header.index("ancestral_state")
    metadata_index = None
    try:
        metadata_index = header.index("metadata")
    except ValueError:
        pass
    table = tables.SiteTable()
    for line in source:
        tokens = line.split(sep)
        if len(tokens) >= 2:
            position = float(tokens[position_index])
            ancestral_state = tokens[ancestral_state_index]
            metadata = b''
            if metadata_index is not None and metadata_index < len(tokens):
                metadata = tables.text_decode_metadata(tokens[metadata_index])
            table.add_row(
                position=position, ancestral_state=ancestral_state, metadata=metadata)
    return table


def parse_mutations(source, sep=None):
    """
    Parse the specified file-like object and return a MutationTable instance.
    The object must contain text with whitespace delimited columns, which are
    labeled with headers and contain columns ``site``, ``node``, and
    ``derived_state``. An optional ``parent`` column may also be supplied.
    Further requirements are described in :class:`MutationTable`.
    """
    header = source.readline().strip("\n").split(sep)
    site_index = header.index("site")
    node_index = header.index("node")
    derived_state_index = header.index("derived_state")
    parent_index = None
    parent = NULL_MUTATION
    try:
        parent_index = header.index("parent")
    except ValueError:
        pass
    metadata_index = None
    try:
        metadata_index = header.index("metadata")
    except ValueError:
        pass
    table = tables.MutationTable()
    for line in source:
        tokens = line.split(sep)
        if len(tokens) >= 3:
            site = int(tokens[site_index])
            node = int(tokens[node_index])
            derived_state = tokens[derived_state_index]
            if parent_index is not None:
                parent = int(tokens[parent_index])
            metadata = b''
            if metadata_index is not None and metadata_index < len(tokens):
                metadata = tables.text_decode_metadata(tokens[metadata_index])
            table.add_row(
                site=site, node=node, derived_state=derived_state, parent=parent,
                metadata=metadata)
    return table


def load_text(nodes, edges, sites=None, mutations=None, sequence_length=0, sep=None):
    """
    Loads a tree sequence from the specified file paths. The files input here
    are in a simple whitespace delimited tabular format such as output by the
    :meth:`.TreeSequence.dump_text` method.  This method is intended as a
    convenient interface for importing external data into msprime; the HDF5
    based file format using by :meth:`msprime.load` will be many times more
    efficient that using the text based formats.

    ``nodes`` and ``edges`` must be a file-like object containing text with
    whitespace delimited columns,  parsable by :func:`parse_nodes` and
    :func:`parse_edges`, respectively.  Further requirements are described
    in :class:`NodeTable` and :class:`EdgeTable`.

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

    edges::

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


    TODO: add description of the field separator argument, linking to the builtin
    string.split function and describe when it is useful.

    :param stream nodes: The file-type object containing text describing a NodeTable.
    :param stream edges: The file-type object containing text
        describing a EdgeTable.
    :param stream sites: The file-type object containing text describing a SiteTable.
    :param stream mutations: The file-type object containing text
        describing a MutationTable.
    :param float sequence_length: The sequence length of the returned tree sequence. If
        not supplied or zero this will be inferred from the set of edges.
    :param str sep: The field separator, as defined by the Python str.split function.
    :return: The tree sequence object containing the information
        stored in the specified file paths.
    :rtype: :class:`msprime.TreeSequence`
    """
    node_table = parse_nodes(nodes, sep=sep)
    edge_table = parse_edges(edges, sep=sep)
    site_table = tables.SiteTable()
    mutation_table = tables.MutationTable()
    if sites is not None:
        site_table = parse_sites(sites, sep=sep)
    if mutations is not None:
        mutation_table = parse_mutations(mutations, sep=sep)
    tables.sort_tables(
        nodes=node_table, edges=edge_table, sites=site_table, mutations=mutation_table)
    return load_tables(
        nodes=node_table, edges=edge_table, sites=site_table, mutations=mutation_table,
        sequence_length=sequence_length)


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

    @property
    def tables(self):
        """
        A copy of the tables underlying this tree sequence. See also
        :meth:`.dump_tables`.

        :return: A :class:`.TableCollection` containing all a copy of the
            tables underlying this tree sequence.
        :rtype: TableCollection
        """
        return self.dump_tables()

    def dump_tables(
            self, nodes=None, edges=None, migrations=None, sites=None,
            mutations=None, provenances=None):
        """
        Copy the contents of the tables underlying the tree sequence to the
        specified objects.

        :param NodeTable nodes: The NodeTable to load the nodes into.
        :param EdgeTable edges: The EdgeTable to load the edges into.
        :param MigrationTable migrations: The MigrationTable to load the migrations into.
        :param SiteTable sites: The SiteTable to load the sites into.
        :param MutationTable mutations: The MutationTable to load the mutations into.
        :param ProvenanceTable provenances: The ProvenanceTable to load the provenances
            into.
        :return: A :class:`.TableCollection` containing all tables underlying
            the tree sequence.
        :rtype: TableCollection
        """
        # TODO document this and test the semantics to passing in new tables
        # as well as returning the updated tables.
        if nodes is None:
            nodes = tables.NodeTable()
        if edges is None:
            edges = tables.EdgeTable()
        if migrations is None:
            migrations = tables.MigrationTable()
        if sites is None:
            sites = tables.SiteTable()
        if mutations is None:
            mutations = tables.MutationTable()
        if provenances is None:
            provenances = tables.ProvenanceTable()
        self._ll_tree_sequence.dump_tables(
            nodes=nodes, edges=edges, migrations=migrations, sites=sites,
            mutations=mutations, provenances=provenances)
        return tables.TableCollection(
            nodes=nodes, edges=edges, migrations=migrations, sites=sites,
            mutations=mutations, provenances=provenances)

    def dump_text(
            self, nodes=None, edges=None, sites=None, mutations=None, provenances=None,
            precision=6):
        """
        Writes a text representation of the tables underlying the tree sequence
        to the specified connections.

        :param stream nodes: The file-like object (having a .write() method) to write
            the NodeTable to.
        :param stream edges: The file-like object to write the EdgeTable to.
        :param stream sites: The file-like object to write the SiteTable to.
        :param stream mutations: The file-like object to write the MutationTable to.
        :param stream provenances: The file-like object to write the ProvenanceTable to.
        :param int precision: The number of digits of precision.
        """

        if nodes is not None:
            print(
                "id", "is_sample", "time", "population", "metadata", sep="\t",
                file=nodes)
            for node in self.nodes():
                row = (
                    "{id:d}\t"
                    "{is_sample:d}\t"
                    "{time:.{precision}f}\t"
                    "{population:d}\t"
                    "{metadata}").format(
                        precision=precision, id=node.id,
                        is_sample=node.is_sample(), time=node.time,
                        population=node.population,
                        metadata=tables.text_encode_metadata(node.metadata))
                print(row, file=nodes)

        if edges is not None:
            print("left", "right", "parent", "child", sep="\t", file=edges)
            for edge in self.edges():
                row = (
                    "{left:.{precision}f}\t"
                    "{right:.{precision}f}\t"
                    "{parent:d}\t"
                    "{child:d}").format(
                        precision=precision, left=edge.left, right=edge.right,
                        parent=edge.parent, child=edge.child)
                print(row, file=edges)

        if sites is not None:
            print("position", "ancestral_state", "metadata", sep="\t", file=sites)
            for site in self.sites():
                row = (
                    "{position:.{precision}f}\t"
                    "{ancestral_state}\t"
                    "{metadata}").format(
                        precision=precision, position=site.position,
                        ancestral_state=site.ancestral_state,
                        metadata=tables.text_encode_metadata(site.metadata))
                print(row, file=sites)

        if mutations is not None:
            print(
                "site", "node", "derived_state", "parent", "metadata",
                sep="\t", file=mutations)
            for site in self.sites():
                for mutation in site.mutations:
                    row = (
                        "{site}\t"
                        "{node}\t"
                        "{derived_state}\t"
                        "{parent}\t"
                        "{metadata}").format(
                            site=mutation.site, node=mutation.node,
                            derived_state=mutation.derived_state,
                            parent=mutation.parent,
                            metadata=tables.text_encode_metadata(mutation.metadata))
                    print(row, file=mutations)

        if provenances is not None:
            print("id", "timestamp", "record", sep="\t", file=provenances)
            for provenance in self.provenances():
                row = (
                    "{id}\t"
                    "{timestamp}\t"
                    "{record}\t").format(
                        id=provenance.id,
                        timestamp=provenance.timestamp,
                        record=provenance.record)
                print(row, file=provenances)

    # num_samples was originally called sample_size, and so we must keep sample_size
    # around as a deprecated alias.
    @property
    def num_samples(self):
        """
        Returns the number of samples in this tree sequence. This is the number
        of sample nodes in each tree.

        :return: The number of sample nodes in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_samples()

    @property
    def sample_size(self):
        # Deprecated alias for num_samples
        return self.num_samples

    def get_sample_size(self):
        # Deprecated alias for num_samples
        return self.num_samples

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
    def num_edges(self):
        return self._ll_tree_sequence.get_num_edges()

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
        return self._ll_tree_sequence.get_num_edges()

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

    @property
    def num_provenances(self):
        return self._ll_tree_sequence.get_num_provenances()

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
        t = [node.time for node in self.nodes()]
        pop = [node.population for node in self.nodes()]
        for e in self.edgesets():
            yield CoalescenceRecord(
                e.left, e.right, e.parent, e.children, t[e.parent], pop[e.parent])

    def migrations(self):
        for j in range(self._ll_tree_sequence.get_num_migrations()):
            yield Migration(*self._ll_tree_sequence.get_migration(j))

    def provenances(self):
        for j in range(self.num_provenances):
            yield self.provenance(j)

    def nodes(self):
        for j in range(self.num_nodes):
            yield self.node(j)

    def edges(self):
        for j in range(self.num_edges):
            left, right, parent, child = self._ll_tree_sequence.get_edge(j)
            yield Edge(left=left, right=right, parent=parent, child=child)

    def edgesets(self):
        # TODO the order that these records are returned in is not well specified.
        # Hopefully this does not matter, and we can just state that the ordering
        # should not be depended on.
        children = collections.defaultdict(set)
        active_edgesets = {}
        for (left, right), edges_out, edges_in in self.edge_diffs():
            # Complete and return any edgesets that are affected by this tree
            # transition
            parents = iter(edge.parent for edge in itertools.chain(edges_out, edges_in))
            for parent in parents:
                if parent in active_edgesets:
                    edgeset = active_edgesets.pop(parent)
                    edgeset.right = left
                    edgeset.children = sorted(children[parent])
                    yield edgeset
            for edge in edges_out:
                children[edge.parent].remove(edge.child)
            for edge in edges_in:
                children[edge.parent].add(edge.child)
            # Update the active edgesets
            for edge in itertools.chain(edges_out, edges_in):
                if len(children[edge.parent]) > 0 and edge.parent not in active_edgesets:
                    active_edgesets[edge.parent] = Edgeset(left, right, edge.parent, [])

        for parent in active_edgesets.keys():
            edgeset = active_edgesets[parent]
            edgeset.right = self.sequence_length
            edgeset.children = sorted(children[edgeset.parent])
            yield edgeset

    def edge_diffs(self):
        iterator = _msprime.TreeDiffIterator(self._ll_tree_sequence)
        for interval, edge_tuples_out, edge_tuples_in in iterator:
            edges_out = [Edge(*e) for e in edge_tuples_out]
            edges_in = [Edge(*e) for e in edge_tuples_in]
            yield interval, edges_out, edges_in

    def site(self, index):
        ll_site = self._ll_tree_sequence.get_site(index)
        pos, ancestral_state, mutations, index, metadata = ll_site
        return Site(
            position=pos, ancestral_state=ancestral_state, index=index,
            mutations=[Mutation(*mutation) for mutation in mutations],
            metadata=metadata)

    def sites(self):
        for j in range(self.num_sites):
            yield self.site(j)

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
        hapgen = _msprime.HaplotypeGenerator(self._ll_tree_sequence)
        j = 0
        # Would use range here except for Python 2.
        while j < self.num_samples:
            yield hapgen.get_haplotype(j)
            j += 1

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

        :param bool as_bytes: If True, the genotype values will be returned
            as a Python bytes object. This is useful in certain situations
            (i.e., directly printing the genotypes) or when numpy is
            not available. Otherwise, genotypes are returned as a numpy
            array (the default).
        :return: An iterator of all :math:`(x, u, j, g)` tuples defining
            the variants in this tree sequence.
        """
        # See comments for the Variant type for discussion on why the
        # present form was chosen.
        check_numpy()
        iterator = _msprime.VariantGenerator(self._ll_tree_sequence)
        for ll_site, genotypes, alleles in iterator:
            pos, ancestral_state, mutations, index, metadata = ll_site
            site = Site(
                position=pos, ancestral_state=ancestral_state, index=index,
                mutations=[Mutation(*mutation) for mutation in mutations],
                metadata=metadata)
            if as_bytes:
                if any(len(allele) > 1 for allele in alleles):
                    raise ValueError(
                        "as_bytes only supported for single-letter alleles")
                bytes_genotypes = np.empty(self.num_samples, dtype=np.uint8)
                lookup = np.array([ord(a[0]) for a in alleles], dtype=np.uint8)
                bytes_genotypes[:] = lookup[genotypes]
                genotypes = bytes_genotypes.tobytes()
            v = Variant(
                position=pos, site=site, index=index, genotypes=genotypes,
                alleles=alleles)
            yield v

    def genotype_matrix(self):
        return self._ll_tree_sequence.get_genotype_matrix()

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

    def node(self, u):
        flags, time, population, metadata = self._ll_tree_sequence.get_node(u)
        return Node(
            id_=u, time=time, population=population, metadata=metadata,
            is_sample=flags & NODE_IS_SAMPLE)

    def provenance(self, id_):
        timestamp, record = self._ll_tree_sequence.get_provenance(id_)
        return Provenance(id_=id_, timestamp=timestamp, record=record)

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

    def simplify(self, samples=None, filter_zero_mutation_sites=True, map_nodes=False):
        """
        Returns a simplified tree sequence that retains only the history of
        the nodes given in the list ``samples``. If ``map_nodes`` is true,
        also return a numpy array mapping the node IDs in this tree sequence to
        their node IDs in the simplified tree tree sequence. If a node u is not
        present in the new tree sequence, the value of this mapping will be
        NULL_NODE (-1).

        **TODO** Document the details of how node IDs are transformed.

        If you wish to convert a set of tables that do not satisfy all
        requirements for building a TreeSequence, then use
        ``simplify_tables()``.

        :param list samples: The list of nodes for which to retain information.
        :param bool filter_zero_mutation_sites: If True, remove any sites that have
            no mutations in the simplified tree sequence. Defaults to True.
        :param bool map_nodes: If True, return a tuple containing the resulting
            tree sequence and a numpy array mapping node IDs in the current tree
            sequence to their corresponding node IDs in the returned tree sequence.
            If False (the default), return only the tree sequence object itself.
        :return: The simplified tree sequence, or (if ``map_nodes`` is True)
            a tuple containing the simplified tree sequence and a numpy array
            mapping source node IDs to their corresponding IDs in the new tree
            sequence.
        :rtype: .TreeSequence or a (.TreeSequence, numpy.array) tuple
        """
        check_numpy()
        t = self.dump_tables()
        if samples is None:
            samples = self.get_samples()
        node_map = tables.simplify_tables(
            samples=samples, sequence_length=self.sequence_length,
            nodes=t.nodes, edges=t.edges,
            sites=t.sites, mutations=t.mutations,
            filter_zero_mutation_sites=filter_zero_mutation_sites)
        # TODO add simplify arguments here??
        t.provenances.add_row(record=json.dumps(
            provenance.get_provenance_dict("simplify", [])))
        new_ts = load_tables(
            nodes=t.nodes, edges=t.edges, migrations=t.migrations, sites=t.sites,
            mutations=t.mutations, provenances=t.provenances,
            sequence_length=self.sequence_length)

        if map_nodes:
            return new_ts, node_map
        else:
            return new_ts

    # Unsupported old methods.
    def diffs(self):
        raise NotImplementedError(
            "This method is no longer supported. Please use the "
            "TreeSequence.edge_diffs() method instead")

    def newick_trees(self, precision=3, breakpoints=None, Ne=1):
        raise NotImplementedError(
            "This method is no longer supported. Please use the SparseTree.newick"
            " method instead")
