# -*- coding: utf-8 -*-
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


# TODO this interface is rubbish. Should have much better printing options.
# TODO we should be use __slots__ here probably.
class SimpleContainer(object):

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return repr(self.__dict__)


class Node(SimpleContainer):
    """
    A :ref:`node <sec_node_table_definition>` in a tree sequence.

    Modifying the attributes in this class will have **no effect** on the
    underlying tree sequence data.

    :ivar id: The integer ID of this node. Varies from 0 to
        :attr:`.TreeSequence.num_nodes` - 1.
    :vartype id: int
    :ivar flags: The bitwise flags for this node.
    :vartype flags: int
    :ivar time: The birth time of the individual represented by this node.
    :vartype float: float
    :ivar population: The integer ID of the population that this node was born in.
    :vartype population: int
    :ivar metadata: The :ref:`metadata <sec_metadata_definition>` for this node.
    :vartype metadata: bytes
    """
    def __init__(
            self, id_=None, flags=0, time=0, population=NULL_POPULATION, metadata=""):
        self.id = id_
        self.time = time
        self.population = population
        self.metadata = metadata
        self.flags = flags

    def is_sample(self):
        """
        Returns True if this node is a sample. This value is derived from the
        ``flag`` variable.

        :rtype: bool
        """
        return self.flags & NODE_IS_SAMPLE


class Edge(SimpleContainer):
    """
    An :ref:`edge <sec_edge_table_definition>` in a tree sequence.

    Modifying the attributes in this class will have **no effect** on the
    underlying tree sequence data.

    :ivar left: The left coordinate of this edge.
    :vartype left: float
    :ivar right: The right coordinate of this edge.
    :vartype right: float
    :ivar parent: The integer ID of the parent node for this edge.
        To obtain further information about a node with a given ID, use
        :meth:`.TreeSequence.node`.
    :vartype parent: int
    :ivar child: The integer ID of the child node for this edge.
        To obtain further information about a node with a given ID, use
        :meth:`.TreeSequence.node`.
    :vartype child: int
    """
    def __init__(self, left, right, parent, child):
        self.left = left
        self.right = right
        self.parent = parent
        self.child = child

    def __repr__(self):
        return "{{left={:.3f}, right={:.3f}, parent={}, child={}}}".format(
            self.left, self.right, self.parent, self.child)


class Site(SimpleContainer):
    """
    A :ref:`site <sec_site_table_definition>` in a tree sequence.

    Modifying the attributes in this class will have **no effect** on the
    underlying tree sequence data.

    :ivar id: The integer ID of this site. Varies from 0 to
        :attr:`.TreeSequence.num_sites` - 1.
    :vartype id: int
    :ivar position: The floating point location of this site in genome coordinates.
        Ranges from 0 (inclusive) to :attr:`.TreeSequence.sequence_length`
        (exclusive).
    :vartype position: float
    :ivar ancestral_state: The ancestral state at this site (i.e., the state
        inherited by nodes, unless mutations occur).
    :vartype ancestral_state: str
    :ivar metadata: The :ref:`metadata <sec_metadata_definition>` for this site.
    :vartype metadata: bytes
    :ivar mutations: The list of mutations at this site. Mutations
        within a site are returned in the order they are specified in the
        underlying :class:`.MutationTable`.
    :vartype mutations: list[:class:`.Mutation`]
    """
    def __init__(self, id_, position, ancestral_state, mutations, metadata):
        self.id = id_
        self.position = position
        self.ancestral_state = ancestral_state
        self.mutations = mutations
        self.metadata = metadata


class Mutation(SimpleContainer):
    """
    A :ref:`mutation <sec_mutation_table_definition>` in a tree sequence.

    Modifying the attributes in this class will have **no effect** on the
    underlying tree sequence data.

    :ivar id: The integer ID of this mutation. Varies from 0 to
        :attr:`.TreeSequence.num_mutations` - 1.
    :vartype id: int
    :ivar site: The integer ID of the site that this mutation occurs at. To obtain
        further information about a site with a given ID use
        :meth:`.TreeSequence.site`.
    :vartype site: int
    :ivar node: The integer ID of the first node that inherits this mutation.
        To obtain further information about a node with a given ID, use
        :meth:`.TreeSequence.node`.
    :vartype node: int
    :ivar derived_state: The derived state for this mutation. This is the state
        inherited by nodes in the subtree rooted at this mutation's node, unless
        another mutation occurs.
    :vartype derived_state: str
    :ivar parent: The integer ID of this mutation's parent mutation. When multiple
        mutations occur at a site along a path in the tree, mutations must
        record the mutation that is immediately above them. If the mutation does
        not have a parent, this is equal to the :const:`NULL_MUTATION` (-1).
        To obtain further information about a mutation with a given ID, use
        :meth:`.TreeSequence.mutation`.
    :vartype parent: int
    :ivar metadata: The :ref:`metadata <sec_metadata_definition>` for this site.
    :vartype metadata: bytes
    """
    def __init__(self, id_, site, node, derived_state, parent, metadata):
        self.id = id_
        self.site = site
        self.node = node
        self.derived_state = derived_state
        self.parent = parent
        self.metadata = metadata


class Migration(SimpleContainer):
    """
    A :ref:`migration <sec_migration_table_definition>` in a tree sequence.

    Modifying the attributes in this class will have **no effect** on the
    underlying tree sequence data.

    :ivar left: The left end of the genomic interval covered by this
        migration (inclusive).
    :vartype left: float
    :ivar right: The right end of the genomic interval covered by this migration
        (exclusive).
    :vartype right: float
    :ivar node: The integer ID of the node involved in this migration event.
        To obtain further information about a node with a given ID, use
        :meth:`.TreeSequence.node`.
    :vartype node: int
    :ivar source: The source population ID.
    :vartype source: int
    :ivar dest: The destination population ID.
    :vartype dest: int
    :ivar time: The time at which this migration occured at.
    :vartype time: float
    """
    def __init__(self, left, right, node, source, dest, time):
        self.left = left
        self.right = right
        self.node = node
        self.source = source
        self.dest = dest
        self.time = time


class Variant(SimpleContainer):
    """
    A variant is represents the observed variation among the samples
    for a given site. A variant consists (a) of a reference to the
    :class:`.Site` instance in question; (b) the **alleles** that may be
    observed at the samples for this site; and (c) the **genotypes**
    mapping sample IDs to the observed alleles.

    Each element in the ``alleles`` tuple is a string, representing the
    actual observed state for a given sample. The first element of this
    tuple is guaranteed to be the same as the site's ``ancestral_state`` value.
    The list of alleles is also guaranteed not to contain any duplicates.
    However, allelic values may be listed that are not referred to by any
    samples. For example, if we have a site that is fixed for the derived state
    (i.e., we have a mutation over the tree root), all genotypes will be 1, but
    the alleles list will be equal to ``('0', '1')``. Other than the
    ancestral state being the first allele, the alleles are listed in
    no particular order, and the ordering should not be relied upon.

    The ``genotypes`` represent the observed allelic states for each sample,
    such that ``var.alleles[var.genotypes[j]]`` gives the string allele
    for sample ID ``j``. Thus, the elements of the genotypes array are
    indexes into the ``alleles`` list. The genotypes are provided in this
    way via a numpy array to enable efficient calculations.

    Modifying the attributes in this class will have **no effect** on the
    underlying tree sequence data.

    :ivar site: The site object for this variant.
    :vartype site: :class:`.Site`
    :ivar alleles: A tuple of the allelic values that may be observed at the
        samples at the current site. The first element of this tuple is always
        the sites's ancestral state.
    :vartype alleles: tuple(str)
    :ivar genotypes: An array of indexes into the list ``alleles``, giving the
        state of each sample at the current site.
    :vartype genotypes: numpy.ndarray
    """
    def __init__(self, site, alleles, genotypes):
        self.site = site
        self.alleles = alleles
        self.genotypes = genotypes
        # Deprecated aliases to avoid breaking existing code.
        self.position = site.position
        self.index = site.id


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


def add_deprecated_mutation_attrs(site, mutation):
    """
    Add in attributes for the older deprecated way of defining
    mutations. These attributes will be removed in future releases
    and are deliberately undocumented in version 0.5.0.
    """
    mutation.position = site.position
    mutation.index = site.id
    return mutation


class SparseTree(object):
    """
    A SparseTree is a single tree in a :class:`.TreeSequence`. The SparseTree
    implementation differs from most tree implementations by using **integer
    node IDs** to refer to nodes rather than objects. Thus, when we wish to
    find the parent of the node with ID '0', we use ``tree.parent(0)``, which
    returns another integer. If '0' does not have a parent in the current tree
    (e.g., if it is a root), then the special value :const:`.NULL_NODE`
    (:math:`-1`) is returned. The children of a node are found using the
    :meth:`.children` method. To obtain information about a particular node,
    one may either use ``tree.tree_sequence.node(u)`` to obtain the
    corresponding :class:`Node` instance, or use the :meth:`.time` or
    :meth:`.population` shorthands. Tree traversals in various orders
    is possible using the :meth:`.SparseTree.nodes` iterator.

    Sparse trees are not intended to be instantiated directly, and are
    obtained as part of a :class:`.TreeSequence` using the
    :meth:`.trees` method.
    """
    def __init__(self, ll_sparse_tree, tree_sequence):
        self._ll_sparse_tree = ll_sparse_tree
        self._tree_sequence = tree_sequence

    @property
    def tree_sequence(self):
        """
        Returns the tree sequence that this tree is from.

        :return: The parent tree sequence for this tree.
        :rtype: :class:`.TreeSequence`
        """
        return self._tree_sequence

    def get_branch_length(self, u):
        # Deprecated alias for branch_length
        return self.branch_length(u)

    def branch_length(self, u):
        """
        Returns the length of the branch (in generations) joining the
        specified node to its parent. This is equivalent to

        >>> tree.time(tree.parent(u)) - tree.time(u)

        Note that this is not related to the value returned by
        :attr:`.length`, which describes the length of the interval
        covered by the tree in genomic coordinates.

        :param int u: The node of interest.
        :return: The branch length from u to its parent.
        :rtype: float
        """
        return self.time(self.get_parent(u)) - self.time(u)

    def get_total_branch_length(self):
        # Deprecated alias for total_branch_length
        return self.total_branch_length

    @property
    def total_branch_length(self):
        """
        Returns the sum of all the branch lengths in this tree (in
        units of generations). This is equivalent to

        >>> sum(
        >>>    tree.branch_length(u) for u in tree.nodes()
        >>>    if u not in self.roots)

        :return: The sum of all the branch lengths in this tree.
        :rtype: float
        """
        return sum(
            self.get_branch_length(u) for u in self.nodes() if u not in self.roots)

    def get_mrca(self, u, v):
        # Deprecated alias for mrca
        return self.mrca(u, v)

    def mrca(self, u, v):
        """
        Returns the most recent common ancestor of the specified nodes.

        :param int u: The first node.
        :param int v: The second node.
        :return: The most recent common ancestor of u and v.
        :rtype: int
        """
        return self._ll_sparse_tree.get_mrca(u, v)

    def get_tmrca(self, u, v):
        # Deprecated alias for tmrca
        return self.tmrca(u, v)

    def tmrca(self, u, v):
        """
        Returns the time of the most recent common ancestor of the specified
        nodes. This is equivalent to::

        >>> tree.time(tree.mrca(u, v))

        :param int u: The first node.
        :param int v: The second node.
        :return: The time of the most recent common ancestor of u and v.
        :rtype: float
        """
        return self.get_time(self.get_mrca(u, v))

    def get_parent(self, u):
        # Deprecated alias for parent
        return self.parent(u)

    def parent(self, u):
        """
        Returns the parent of the specified node. Returns
        the :const:`.NULL_NODE` if u is the root or is not a node in
        the current tree.

        :param int u: The node of interest.
        :return: The parent of u.
        :rtype: int
        """
        return self._ll_sparse_tree.get_parent(u)

    # Quintuply linked tree structure.

    def left_child(self, u):
        return self._ll_sparse_tree.get_left_child(u)

    def right_child(self, u):
        return self._ll_sparse_tree.get_right_child(u)

    def left_sib(self, u):
        return self._ll_sparse_tree.get_left_sib(u)

    def right_sib(self, u):
        return self._ll_sparse_tree.get_right_sib(u)

    # TODO do we also have right_root?
    @property
    def left_root(self):
        return self._ll_sparse_tree.get_left_root()

    def get_children(self, u):
        # Deprecated alias for self.children
        return self.children(u)

    def children(self, u):
        """
        Returns the children of the specified node ``u`` as a tuple of integer node IDs.
        If ``u`` is a leaf, return the empty tuple.

        :param int u: The node of interest.
        :return: The children of ``u`` as a tuple of integers
        :rtype: tuple(int)
        """
        return self._ll_sparse_tree.get_children(u)

    def get_time(self, u):
        # Deprecated alias for self.time
        return self.time(u)

    def time(self, u):
        """
        Returns the time of the specified node in generations.
        Equivalent to ``tree.tree_sequence.node(u).time``.

        :param int u: The node of interest.
        :return: The time of u.
        :rtype: float
        """
        return self._ll_sparse_tree.get_time(u)

    def get_population(self, u):
        # Deprecated alias for self.population
        return self.population(u)

    def population(self, u):
        """
        Returns the population associated with the specified node.
        Equivalent to ``tree.tree_sequence.node(u).population``.

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

        :return: The list of roots in this tree.
        :rtype: list
        """
        roots = []
        u = self.left_root
        while u != NULL_NODE:
            roots.append(u)
            u = self.right_sib(u)
        return roots

    def get_root(self):
        # Deprecated alias for self.root
        return self.root

    @property
    def root(self):
        """
        The root of this tree. If the tree contains multiple roots, a ValueError is
        raised indicating that the :attr:`.roots` attribute should be used instead.

        :return: The root node.
        :rtype: int
        :raises: :class:`ValueError` if this tree contains more than one root.
        """
        root = self.left_root
        if root != NULL_NODE and self.right_sib(root) != NULL_NODE:
            raise ValueError("More than one root exists. Use tree.roots instead")
        return root

    def get_index(self):
        # Deprecated alias for self.index
        return self.index

    @property
    def index(self):
        """
        Returns the index this tree occupies in the parent tree sequence.
        This index is zero based, so the first tree in the sequence has index 0.

        :return: The index of this tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_index()

    def get_interval(self):
        # Deprecated alias for self.interval
        return self.interval

    @property
    def interval(self):
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

    def get_length(self):
        # Deprecated alias for self.length
        return self.length

    @property
    def length(self):
        """
        Returns the length of the genomic interval that this tree represents.
        This is defined as :math:`r - l`, where :math:`(l, r)` is the genomic
        interval returned by :attr:`.interval`.

        :return: The length of the genomic interval covered by this tree.
        :rtype: int
        """
        left, right = self.get_interval()
        return right - left

    # The sample_size (or num_samples) is really a property of the tree sequence,
    # and so we should provide access to this via a tree.tree_sequence.num_samples
    # property access. However, we can't just remove the method as a lot of code
    # may depend on it. To complicate things a bit more, sample_size has been
    # changed to num_samples elsewhere for consistency. We can't do this here
    # because there is already a num_samples method which returns the number of
    # samples below a particular node. The best thing to do is probably to
    # undocument the sample_size property, but keep it around for ever.

    def get_sample_size(self):
        # Deprecated alias for self.sample_size
        return self.sample_size

    @property
    def sample_size(self):
        """
        Returns the sample size for this tree. This is the number of sample
        nodes in the tree.

        :return: The number of sample nodes in the tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_sample_size()

    def draw(
            self, path=None, width=None, height=None,
            node_labels=None, node_colours=None,
            mutation_labels=None, mutation_colours=None,
            format=None):
        """
        Returns a drawing of this tree.

        When working in a Jupyter notebook, use the ``IPython.display.SVG``
        function to display the SVG output from this function inline in the notebook::

            >>> SVG(tree.draw())

        The unicode format uses unicode `box drawing characters
        <https://en.wikipedia.org/wiki/Box-drawing_character>`_ to render the tree.
        This allows rendered trees to be printed out to the terminal::

            >>> print(tree.draw(format="unicode"))
              6
            ┏━┻━┓
            ┃   5
            ┃ ┏━┻┓
            ┃ ┃  4
            ┃ ┃ ┏┻┓
            3 0 1 2

        The ``node_labels`` argument allows the user to specify custom labels
        for nodes, or no labels at all::

            >>> print(tree.draw(format="unicode", node_labels={}))
              ┃
            ┏━┻━┓
            ┃   ┃
            ┃ ┏━┻┓
            ┃ ┃  ┃
            ┃ ┃ ┏┻┓
            ┃ ┃ ┃ ┃

        :param str path: The path to the file to write the output. If None, do not
            write to file.
        :param int width: The width of the image in pixels. If not specified, either
            defaults to the minimum size required to depict the tree (text formats)
            or 200 pixels.
        :param int height: The height of the image in pixels. If not specified, either
            defaults to the minimum size required to depict the tree (text formats)
            or 200 pixels.
        :param map node_labels: If specified, show custom labels for the nodes
            that are present in the map. Any nodes not specified in the map will
            not have a node label.
        :param map node_colours: If specified, show custom colours for nodes. (Only
            supported in the SVG format.)
        :param str format: The format of the returned image. Currently supported
            are 'svg', 'ascii' and 'unicode'.
        :return: A representation of this tree in the requested format.
        :rtype: str
        """
        output = drawing.draw_tree(
            self, format=format, width=width, height=height,
            node_labels=node_labels, node_colours=node_colours,
            mutation_labels=mutation_labels, mutation_colours=mutation_colours)
        if path is not None:
            with open(path, "w") as f:
                f.write(output)
        return output

    def get_num_mutations(self):
        return self.num_mutations

    @property
    def num_mutations(self):
        """
        Returns the total number of mutations across all sites on this tree.

        :return: The total number of mutations over all sites on this tree.
        :rtype: int
        """
        return sum(len(site.mutations) for site in self.sites())

    @property
    def num_sites(self):
        """
        Returns the number of sites on this tree.

        :return: The number of sites on this tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_num_sites()

    def sites(self):
        """
        Returns an iterator over all the :ref:`sites <sec_site_table_definition>`
        in this tree. Sites are returned in order of increasing ID
        (and also position). See the :class:`Site` class for details on
        the available fields for each site.

        :return: An iterator over all sites in this tree.
        :rtype: iter(:class:`.Site`)
        """
        # TODO change the low-level API to just return the IDs of the sites.
        for ll_site in self._ll_sparse_tree.get_sites():
            _, _, _, id_, _ = ll_site
            yield self.tree_sequence.site(id_)

    def mutations(self):
        """
        Returns an iterator over all the
        :ref:`mutations <sec_mutation_table_definition>` in this tree.
        Mutations are returned in order of nondecreasing site ID.
        See the :class:`Mutation` class for details on the available fields for
        each mutation.

        The returned iterator is equivalent to iterating over all sites
        and all mutations in each site, i.e.::

            >>> for site in tree.sites():
            >>>     for mutation in site.mutations:
            >>>         yield mutation

        :return: An iterator over all mutations in this tree.
        :rtype: iter(:class:`.Mutation`)
        """
        for site in self.sites():
            for mutation in site.mutations:
                yield add_deprecated_mutation_attrs(site, mutation)

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

    def get_num_samples(self, u=None):
        # Deprecated alias for num_samples.
        return self.num_samples(u)

    def num_samples(self, u=None):
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

    def get_num_tracked_samples(self, u=None):
        # Deprecated alias for num_tracked_samples
        return self.num_tracked_samples(u)

    def num_tracked_samples(self, u=None):
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
        :return: An iterator over the nodes in the tree in some traversal order.
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
        """
        Returns a `newick encoding <https://en.wikipedia.org/wiki/Newick_format>`_
        of this tree. Leaf nodes are labelled with their numerical ID + 1,
        and internal nodes are not labelled.

        This method is currently primarily for ms-compatibility and
        is not intended as a consistent means of data interchange.

        :param int precision: The numerical precision with which branch lengths are
            printed.
        :param float time_scale: A value which all branch lengths are multiplied by.
        :return: A newick representation of this tree.
        :rtype: str
        """
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
    Loads a tree sequence from the specified file path. This file must be in the
    :ref:`HDF5 file format <sec_hdf5_file_format>` produced by the
    :meth:`.TreeSequence.dump` method.

    :param str path: The file path of the HDF5 file containing the
        tree sequence we wish to load.
    :return: The tree sequence object containing the information
        stored in the specified file path.
    :rtype: :class:`msprime.TreeSequence`
    """
    return TreeSequence.load(path)


def load_tables(
        nodes, edges, migrations=None, sites=None, mutations=None,
        provenances=None, sequence_length=0):
    """
    Loads the tree sequence data from the specified table objects, and
    returns the resulting :class:`.TreeSequence` object. These tables
    must fulfil the properties required for an input tree sequence as
    described in the :ref:`sec_valid_tree_sequence_requirements` section.

    The ``sequence_length`` parameter determines the
    :attr:`.TreeSequence.sequence_length` of the returned tree sequence. If it
    is 0 or not specified, the value is taken to be the maximum right
    coordinate of the input edges. This parameter is useful in degenerate
    situations (such as when there are zero edges), but can usually be ignored.

    :param NodeTable nodes: The :ref:`node table <sec_node_table_definition>`
        (required).
    :param EdgeTable edges: The :ref:`edge table <sec_edge_table_definition>`
        (required).
    :param MigrationTable migrations: The :ref:`migration table
        <sec_migration_table_definition>` (optional).
    :param SiteTable sites: The :ref:`site table <sec_site_table_definition>`
        (optional; but if supplied, ``mutations`` must also be specified).
    :param MutationTable mutations: The :ref:`mutation table
        <sec_mutation_table_definition>` (optional; but if supplied, ``sites``
        must also be specified).
    :param ProvenanceTable provenances: The :ref:`provenance table
        <sec_provenance_table_definition>` (optional).
    :param float sequence_length: The sequence length of the returned tree sequence. If
        not supplied or zero this will be inferred from the set of edges.
    :return: A :class:`.TreeSequence` consistent with the specified tables.
    :rtype: TreeSequence
    """
    # TODO update the low-level module to accept None and remove this
    kwargs = {"nodes": nodes, "edges": edges, "sequence_length": sequence_length}
    if migrations is not None:
        kwargs["migrations"] = migrations
    if sites is not None:
        kwargs["sites"] = sites
    if mutations is not None:
        kwargs["mutations"] = mutations
    if provenances is not None:
        kwargs["provenances"] = provenances
    return TreeSequence.load_tables(**kwargs)


def parse_nodes(source, strict=True):
    """
    Parse the specified file-like object containing a whitespace delimited
    description of a node table and returns the corresponding :class:`NodeTable`
    instance. See the :ref:`node text format <sec_node_text_format>` section
    for the details of the required format and the
    :ref:`node table definition <sec_node_table_definition>` section for the
    required properties of the contents.

    See :func:`.load_text` for a detailed explanation of the ``strict`` parameter.

    :param stream source: The file-like object containing the text.
    :param bool strict: If True, require strict tab delimiting (default). If
        False, a relaxed whitespace splitting algorithm is used.
    """
    sep = None
    if strict:
        sep = "\t"
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


def parse_edges(source, strict=True):
    """
    Parse the specified file-like object containing a whitespace delimited
    description of a edge table and returns the corresponding :class:`EdgeTable`
    instance. See the :ref:`edge text format <sec_edge_text_format>` section
    for the details of the required format and the
    :ref:`edge table definition <sec_edge_table_definition>` section for the
    required properties of the contents.

    See :func:`.load_text` for a detailed explanation of the ``strict`` parameter.

    :param stream source: The file-like object containing the text.
    :param bool strict: If True, require strict tab delimiting (default). If
        False, a relaxed whitespace splitting algorithm is used.
    """
    sep = None
    if strict:
        sep = "\t"
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


def parse_sites(source, strict=True):
    """
    Parse the specified file-like object containing a whitespace delimited
    description of a site table and returns the corresponding :class:`SiteTable`
    instance. See the :ref:`site text format <sec_site_text_format>` section
    for the details of the required format and the
    :ref:`site table definition <sec_site_table_definition>` section for the
    required properties of the contents.

    See :func:`.load_text` for a detailed explanation of the ``strict`` parameter.

    :param stream source: The file-like object containing the text.
    :param bool strict: If True, require strict tab delimiting (default). If
        False, a relaxed whitespace splitting algorithm is used.
    """
    sep = None
    if strict:
        sep = "\t"
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


def parse_mutations(source, strict=True):
    """
    Parse the specified file-like object containing a whitespace delimited
    description of a mutation table and returns the corresponding :class:`MutationTable`
    instance. See the :ref:`mutation text format <sec_mutation_text_format>` section
    for the details of the required format and the
    :ref:`mutation table definition <sec_mutation_table_definition>` section for the
    required properties of the contents.

    See :func:`.load_text` for a detailed explanation of the ``strict`` parameter.

    :param stream source: The file-like object containing the text.
    :param bool strict: If True, require strict tab delimiting (default). If
        False, a relaxed whitespace splitting algorithm is used.
    """
    sep = None
    if strict:
        sep = "\t"
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


def load_text(nodes, edges, sites=None, mutations=None, sequence_length=0, strict=True):
    """
    Parses the tree sequence data from the specified file-like objects, and
    returns the resulting :class:`.TreeSequence` object. The format
    for these files is documented in the :ref:`sec_text_file_format` section,
    and is produced by the :meth:`.TreeSequence.dump_text` method. Further
    properties required for an input tree sequence are described in the
    :ref:`sec_valid_tree_sequence_requirements` section. This method is intended as a
    convenient interface for importing external data into msprime; the HDF5
    based file format using by :meth:`msprime.load` is many times more
    efficient than this text format.

    The ``nodes`` and ``edges`` parameters are mandatory and must be file-like
    objects containing text with whitespace delimited columns,  parsable by
    :func:`parse_nodes` and :func:`parse_edges`, respectively. ``sites`` and
    ``mutations`` are optional, and must be parsable by :func:`parse_sites` and
    :func:`parse_mutations`, respectively.

    The ``sequence_length`` parameter determines the
    :attr:`.TreeSequence.sequence_length` of the returned tree sequence. If it
    is 0 or not specified, the value is taken to be the maximum right
    coordinate of the input edges. This parameter is useful in degenerate
    situations (such as when there are zero edges), but can usually be ignored.

    The ``strict`` parameter controls the field delimiting algorithm that
    is used. If ``strict`` is True (the default), we require exactly one
    tab character separating each field. If ``strict`` is False, a more relaxed
    whitespace delimiting algorithm is used, such that any run of whitespace
    is regarded as a field separator. In most situations, ``strict=False``
    is more convenient, but it can lead to error in certain situations. For
    example, if a deletion is encoded in the mutation table this will not
    be parseable when ``strict=False``.

    After parsing the tables, :func:`sort_tables` is called to ensure that
    the loaded tables satisfy the tree sequence :ref:`ordering requirements
    <sec_valid_tree_sequence_requirements>`. Note that this may result in the
    IDs of various entities changing from their positions in the input file.

    :param stream nodes: The file-like object containing text describing a
        :class:`.NodeTable`.
    :param stream edges: The file-like object containing text
        describing an :class:`.EdgeTable`.
    :param stream sites: The file-like object containing text describing a
        :class:`.SiteTable`.
    :param stream mutations: The file-like object containing text
        describing a :class:`MutationTable`.
    :param float sequence_length: The sequence length of the returned tree sequence. If
        not supplied or zero this will be inferred from the set of edges.
    :param bool strict: If True, require strict tab delimiting (default). If
        False, a relaxed whitespace splitting algorithm is used.
    :return: The tree sequence object containing the information
        stored in the specified file paths.
    :rtype: :class:`msprime.TreeSequence`
    """
    node_table = parse_nodes(nodes, strict=strict)
    edge_table = parse_edges(edges, strict=strict)
    site_table = tables.SiteTable()
    mutation_table = tables.MutationTable()
    if sites is not None:
        site_table = parse_sites(sites, strict=strict)
    if mutations is not None:
        mutation_table = parse_mutations(mutations, strict=strict)
    tables.sort_tables(
        nodes=node_table, edges=edge_table, sites=site_table, mutations=mutation_table)
    return load_tables(
        nodes=node_table, edges=edge_table, sites=site_table, mutations=mutation_table,
        sequence_length=sequence_length)


class TreeSequence(object):
    """
    A single tree sequence, as defined by the :ref:`data model <sec_data_model>`.
    A TreeSequence instance can be created from a set of
    :ref:`tables <sec_table_definitions>` using :func:`.load_tables`; or loaded
    from a set of text files using :func:`.load_text`; or, loaded from a
    native file using :func:`load`.

    TreeSequences are immutable. To change the data held in a particular
    tree sequence, first output the informatinn to a set of tables
    (using :meth:`.dump_tables`), edit those tables using the
    :ref:`tables api <sec_tables_api>`, and create a new tree sequence using
    :func:`.load_tables`.

    The :meth:`.trees` method iterates over all trees in a tree sequence, and
    the :meth:`.variants` method iterates over all sites and their genotypes.
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
        return self.get_sequence_length()

    def get_sequence_length(self):
        return self._ll_tree_sequence.get_sequence_length()

    @property
    def num_edges(self):
        """
        Returns the number of :ref:`edges <sec_edge_table_definition>` in this
        tree sequence.

        :return: The number of edges in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_edges()

    def get_num_trees(self):
        # Deprecated alias for self.num_trees
        return self.num_trees

    @property
    def num_trees(self):
        """
        Returns the number of distinct trees in this tree sequence. This
        is equal to the number of trees returned by the :meth:`.trees`
        method.

        :return: The number of trees in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_trees()

    def get_num_sites(self):
        # Deprecated alias for self.num_sites
        return self._ll_tree_sequence.get_num_sites()

    @property
    def num_sites(self):
        """
        Returns the number of sites in this tree sequence.

        :return: The number of sites in this tree sequence.
        :rtype: int
        """
        return self.get_num_sites()

    def get_num_mutations(self):
        # Deprecated alias for self.num_mutations
        return self.num_mutations

    @property
    def num_mutations(self):
        """
        Returns the number of :ref:`mutations <sec_mutation_table_definition>`
        in this tree sequence.

        :return: The number of mutations in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_mutations()

    def get_num_nodes(self):
        # Deprecated alias for self.num_nodes
        return self.num_nodes

    @property
    def num_nodes(self):
        """
        Returns the number of :ref:`nodes <sec_node_table_definition>` in
        this tree sequence.

        :return: The number of nodes in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_nodes()

    @property
    def num_provenances(self):
        """
        Returns the number of :ref:`provenances <sec_provenance_table_definition>`
        in this tree sequence.

        :return: The number of provenances in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_provenances()

    @property
    def num_migrations(self):
        """
        Returns the number of :ref:`migrations <sec_migration_table_definition>`
        in this tree sequence.

        :return: The number of migrations in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_migrations()

    def migrations(self):
        """
        Returns an iterator over all the
        :ref:`migrations <sec_migration_table_definition>` in this tree sequence.

        Migrations are returned in nondecreasing order of the ``time`` value.

        :return: An iterator over all migrations.
        :rtype: iter(:class:`.Migration`)
        """
        for j in range(self._ll_tree_sequence.get_num_migrations()):
            yield Migration(*self._ll_tree_sequence.get_migration(j))

    def provenances(self):
        for j in range(self.num_provenances):
            yield self.provenance(j)

    def nodes(self):
        """
        Returns an iterator over all the :ref:`nodes <sec_node_table_definition>`
        in this tree sequence.

        :return: An iterator over all nodes.
        :rtype: iter(:class:`.Node`)
        """
        for j in range(self.num_nodes):
            yield self.node(j)

    def edges(self):
        """
        Returns an iterator over all the :ref:`edges <sec_edge_table_definition>`
        in this tree sequence. Edges are returned in the order required
        for a :ref:`valid tree sequence <sec_valid_tree_sequence_requirements>`. So,
        edges are guaranteed to be ordered such that (a) all parents with a
        given ID are contiguous; (b) edges are returned in non-descreasing
        order of parent time; (c) within the edges for a given parent, edges
        are sorted first by child ID and then by left coordinate.

        :return: An iterator over all edges.
        :rtype: iter(:class:`.Edge`)
        """
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

    def sites(self):
        """
        Returns an iterator over all the :ref:`sites <sec_site_table_definition>`
        in this tree sequence. Sites are returned in order of increasing ID
        (and also position). See the :class:`Site` class for details on
        the available fields for each site.

        :return: An iterator over all sites.
        :rtype: iter(:class:`.Site`)
        """
        for j in range(self.num_sites):
            yield self.site(j)

    def mutations(self):
        """
        Returns an iterator over all the
        :ref:`mutations <sec_mutation_table_definition>` in this tree sequence.
        Mutations are returned in order of nondecreasing site ID.
        See the :class:`Mutation` class for details on the available fields for
        each mutation.

        The returned iterator is equivalent to iterating over all sites
        and all mutations in each site, i.e.::

            >>> for site in tree_sequence.sites():
            >>>     for mutation in site.mutations:
            >>>         yield mutation

        :return: An iterator over all mutations in this tree sequence.
        :rtype: iter(:class:`.Mutation`)
        """
        for site in self.sites():
            for mutation in site.mutations:
                yield add_deprecated_mutation_attrs(site, mutation)

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

    def first(self):
        """
        Returns the first tree in this :class:`.TreeSequence`. To iterate over all
        trees in the sequence, use the :meth:`.trees` method.

        Currently does not support the extra options for the :meth:`.trees` method.

        :return: The first tree in this tree sequence.
        :rtype: :class:`.SparseTree`.
        """
        return next(self.trees())

    def trees(
            self, tracked_samples=None, sample_counts=True, sample_lists=False,
            tracked_leaves=None, leaf_counts=None, leaf_lists=None):
        """
        Returns an iterator over the trees in this tree sequence. Each value
        returned in this iterator is an instance of :class:`.SparseTree`.

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
            # TODO remove this when we allow numpy arrays in the low-level API.
            kwargs["tracked_samples"] = list(tracked_samples)
        ll_sparse_tree = _msprime.SparseTree(self._ll_tree_sequence, **kwargs)
        iterator = _msprime.SparseTreeIterator(ll_sparse_tree)
        sparse_tree = SparseTree(ll_sparse_tree, self)
        for _ in iterator:
            yield sparse_tree

    def haplotypes(self):
        """
        Returns an iterator over the haplotypes resulting from the trees
        and mutations in this tree sequence as a string.
        The iterator returns a total of :math:`n` strings, each of which
        contains :math:`s` characters (:math:`n` is the sample size
        returned by :attr:`msprime.TreeSequence.num_samples` and
        :math:`s` is the number of sites returned by
        :attr:`msprime.TreeSequence.num_sites`). The first
        string returned is the haplotype for sample `0`, and so on.
        For a given haplotype ``h``, the value of ``h[j]`` is the observed
        allelic state at site ``j``.

        See also the :meth:`variants` iterator for site-centric access
        to sample genotypes.

        This method is only supported for single-letter alleles.

        :return: An iterator over the haplotype strings for the samples in
            this tree sequence.
        :rtype: iter
        :raises: LibraryError if called on a tree sequence containing
            multiletter alleles.
        """
        hapgen = _msprime.HaplotypeGenerator(self._ll_tree_sequence)
        j = 0
        # Would use range here except for Python 2.
        while j < self.num_samples:
            yield hapgen.get_haplotype(j)
            j += 1

    def variants(self, as_bytes=False):
        """
        Returns an iterator over the variants in this tree sequence. See the
        :class:`Variant` class for details on the fields of each returned
        object. By default the ``genotypes`` for the variants are numpy arrays,
        corresponding to indexes into the ``alleles`` array. If the
        ``as_bytes`` parameter is true, these allelic values are recorded
        directly into a bytes array.

        .. note::
            The ``as_bytes`` parameter is kept as a compatibility
            option for older code. It is not the recommended way of
            accessing variant data, and will be deprecated in a later
            release. Another method will be provided to obtain the allelic
            states for each site directly.

        :param bool as_bytes: If True, the genotype values will be returned
            as a Python bytes object. This is useful in certain situations
            (i.e., directly printing the genotypes) or when numpy is
            not available. Otherwise, genotypes are returned as a numpy
            array (the default).
        :return: An iterator of all variants this tree sequence.
        :rtype: iter(:class:`Variant`)
        """
        # See comments for the Variant type for discussion on why the
        # present form was chosen.
        check_numpy()
        iterator = _msprime.VariantGenerator(self._ll_tree_sequence)
        for site_id, genotypes, alleles in iterator:
            site = self.site(site_id)
            if as_bytes:
                if any(len(allele) > 1 for allele in alleles):
                    raise ValueError(
                        "as_bytes only supported for single-letter alleles")
                bytes_genotypes = np.empty(self.num_samples, dtype=np.uint8)
                lookup = np.array([ord(a[0]) for a in alleles], dtype=np.uint8)
                bytes_genotypes[:] = lookup[genotypes]
                genotypes = bytes_genotypes.tobytes()
            yield Variant(site, alleles, genotypes)

    def genotype_matrix(self):
        """
        Returns an :math:`m \\times n` numpy array of the genotypes in this
        tree sequence, where :math:`m` is the number of site and :math:`n`
        the number of samples. The genotypes are the indexes into the array
        of ``alleles``, as described for the :class:`Variant` class. The value
        0 always corresponds to the ancestal state, and values > 0 represent
        distinct derived states.

        .. warning::
            This method can consume a **very large** amount of memory! If
            all genotypes are not needed at once, it is usually better to
            access them sequentially using the :meth:`.variants` iterator.

        :return: The full matrix of genotypes.
        :rtype: numpy.ndarray (dtype=np.uint8)
        """
        return self._ll_tree_sequence.get_genotype_matrix()

    def get_pairwise_diversity(self, samples=None):
        # Deprecated alias for self.pairwise_diversity
        return self.pairwise_diversity(samples)

    def pairwise_diversity(self, samples=None):
        """
        Returns the value of :math:`\pi`, the pairwise nucleotide site diversity,
        which is the average number of mutations that differ between a randomly
        chosen pair of samples.  If `samples` is specified, calculate the
        diversity within this set.

        .. note:: This method does not currently support sites that have more
            than one mutation. Using it on such a tree sequence will raise
            a LibraryError with an "Unsupported operation" message.

        :param iterable samples: The set of samples within which we calculate
            the diversity. If None, calculate diversity within the entire sample.
        :return: The pairwise nucleotide site diversity.
        :rtype: float
        """
        if samples is None:
            samples = self.samples()
        return self._ll_tree_sequence.get_pairwise_diversity(list(samples))

    def node(self, id_):
        """
        Returns the :ref:`node <sec_node_table_definition>` in this tree sequence
        with the specified ID.

        :rtype: :class:`.Node`
        """
        flags, time, population, metadata = self._ll_tree_sequence.get_node(id_)
        return Node(
            id_=id_, flags=flags, time=time, population=population, metadata=metadata)

    def mutation(self, id_):
        """
        Returns the :ref:`mutation <sec_mutation_table_definition>` in this tree sequence
        with the specified ID.

        :rtype: :class:`.Mutation`
        """
        ll_mut = self._ll_tree_sequence.get_mutation(id_)
        return Mutation(
            id_=id_, site=ll_mut[0], node=ll_mut[1], derived_state=ll_mut[2],
            parent=ll_mut[3], metadata=ll_mut[4])

    def site(self, id_):
        """
        Returns the :ref:`site <sec_site_table_definition>` in this tree sequence
        with the specified ID.

        :rtype: :class:`.Site`
        """
        ll_site = self._ll_tree_sequence.get_site(id_)
        pos, ancestral_state, ll_mutations, _, metadata = ll_site
        mutations = [self.mutation(mut_id) for mut_id in ll_mutations]
        return Site(
            id_=id_, position=pos, ancestral_state=ancestral_state,
            mutations=mutations, metadata=metadata)

    def provenance(self, id_):
        timestamp, record = self._ll_tree_sequence.get_provenance(id_)
        return Provenance(id_=id_, timestamp=timestamp, record=record)

    def get_samples(self, population_id=None):
        # Deprecated alias for samples()
        return self.samples(population_id)

    def samples(self, population=None, population_id=None):
        """
        Returns an array of the sample node IDs in this tree sequence. If the
        ``population`` parameter is specified, only return sample IDs from this
        population.

        :param int population: The population of interest. If None,
            return all samples.
        :param int population_id: Deprecated alias for ``population``.
        :return: A numpy array of the node IDs for the samples of interest.
        :rtype: numpy.ndarray (dtype=np.int32)
        """
        if population is not None and population_id is not None:
            raise ValueError(
                "population_id and population are aliases. Cannot specify both")
        if population_id is not None:
            population = population_id
        # TODO the low-level tree sequence should perform this operation natively
        # and return a numpy array.
        samples = self._ll_tree_sequence.get_samples()
        if population is not None:
            samples = [
                u for u in samples if self.get_population(u) == population]
        return np.array(samples, dtype=np.int32)

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
        their node IDs in the simplified tree tree sequence. If a node ``u`` is not
        present in the new tree sequence, the value of this mapping will be
        NULL_NODE (-1).

        In the returned tree sequence, the node with ID ``0`` corresponds to
        ``samples[0]``, node ``1`` corresponds to ``samples[1]``, and so on. Node
        IDs in the returned tree sequence are then allocated sequentially
        in time order. Note that this does **not** necessarily mean that nodes
        in the returned tree sequence will be in strict time order (as we
        may have internal or ancient samples).

        If you wish to convert a set of tables that do not satisfy all
        requirements for building a TreeSequence, then use
        :func:`.simplify_tables()`.

        :param list samples: The list of nodes for which to retain information. This
            may be a numpy array (or array-like) object (dtype=np.int32).
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

    ############################################
    #
    # Deprecated APIs. These are either already unsupported, or will be unsupported in a
    # later release.
    #
    ############################################

    def get_time(self, u):
        # Deprecated. Use ts.node(u).time
        if u < 0 or u >= self.get_num_nodes():
            raise ValueError("ID out of bounds")
        node = self.node(u)
        return node.time

    def get_population(self, u):
        # Deprecated. Use ts.node(u).population
        if u < 0 or u >= self.get_num_nodes():
            raise ValueError("ID out of bounds")
        node = self.node(u)
        return node.population

    def records(self):
        # Deprecated. Use either ts.edges() or ts.edgesets().
        t = [node.time for node in self.nodes()]
        pop = [node.population for node in self.nodes()]
        for e in self.edgesets():
            yield CoalescenceRecord(
                e.left, e.right, e.parent, e.children, t[e.parent], pop[e.parent])

    # Unsupported old methods.

    def get_num_records(self):
        raise NotImplementedError(
            "This method is no longer supported. Please use the "
            "TreeSequence.num_edges if possible to work with edges rather "
            "than coalescence records. If not, please use len(list(ts.edgesets())) "
            "which should return the number of coalescence records, as previously "
            "defined. Please open an issue on GitHub if this is "
            "important for your workflow.")

    def diffs(self):
        raise NotImplementedError(
            "This method is no longer supported. Please use the "
            "TreeSequence.edge_diffs() method instead")

    def newick_trees(self, precision=3, breakpoints=None, Ne=1):
        raise NotImplementedError(
            "This method is no longer supported. Please use the SparseTree.newick"
            " method instead")
