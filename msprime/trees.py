#
# Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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

import _msprime

# Make the low-level generator appear like its from this module
from _msprime import RandomGenerator

NULL_NODE = -1
NULL_POPULATION = -1


CoalescenceRecord = collections.namedtuple(
    "CoalescenceRecord",
    ["left", "right", "node", "children", "time", "population"])


class TreeDrawer(object):
    """
    A class to draw sparse trees in SVG format.
    """
    def __init__(self, tree, width, height, show_times):
        self._width = width
        self._height = height
        self._show_times = show_times
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
        self._leaf_x = 1
        self._assign_x_coordinates(self._tree.get_root())
        self._mutations = []
        for pos, u in tree.mutations():
            x = self._x_coords[u], self._y_coords[u]
            v = tree.get_parent(u)
            y = self._x_coords[v], self._y_coords[v]
            z = (x[0] + y[0]) / 2, (x[1] + y[1]) / 2
            self._mutations.append(z)

    def write(self, path):
        """
        Writes the SVG description of this tree to the specified
        path.
        """
        dwg = svgwrite.Drawing(
            path, size=(self._width, self._height), debug=True)
        lines = dwg.add(dwg.g(id='lines', stroke='black'))
        labels = dwg.add(dwg.g(font_size=14))
        for u in self._tree.nodes():
            v = self._tree.get_parent(u)
            x = self._x_coords[u], self._y_coords[u]
            dwg.add(dwg.circle(center=x, r=3))
            dx = [0]
            dy = None
            if self._tree.is_leaf(u):
                dy = [20]
                dx = [-4]
            elif v == 0:
                dy = [-5]
            else:
                dx = [-20]
            labels.add(dwg.text(str(u), x, dx=dx, dy=dy))
            if self._show_times and self._tree.is_internal(u):
                dx[0] += 25
                labels.add(dwg.text(
                    "t = {:.2f}".format(self._tree.get_time(u)), x, dx=dx,
                    dy=dy)
                )
            if v != NULL_NODE:
                y = self._x_coords[v], self._y_coords[v]
                lines.add(dwg.line(x, y))
        for x in self._mutations:
            dwg.add(dwg.rect(insert=x, size=(6, 6), fill="red"))
        dwg.save()

    def _assign_x_coordinates(self, node):
        """
        Assign x coordinates to all nodes underneath this node.
        """
        if self._tree.is_internal(node):
            children = self._tree.get_children(node)
            for c in children:
                self._assign_x_coordinates(c)
            # We now have x coords for both children
            c1 = self._x_coords[children[0]]
            c2 = self._x_coords[children[1]]
            a = min(c1, c2)
            b = max(c1, c2)
            self._x_coords[node] = (a + (b - a) / 2)
        else:
            self._x_coords[node] = self._leaf_x * self._x_scale
            self._leaf_x += 1


# TODO:
# - Pre, post, and inorder traversals of the nodes as iterators.
# - Pickle and copy support
class SparseTree(object):
    """
    A SparseTree is a single tree in a :class:`.TreeSequence`. In a sparse tree
    for a sample of size :math:`n`, the leaves are nodes :math:`0` to :math:`n
    - 1` inclusive and internal nodes are integers :math:`\geq n`. The value of
    these nodes is strictly increasing as we ascend the tree and the root of
    the tree is the node with the largest value that is reachable from  the
    leaves. Each node in the tree has a parent, which is non-zero for all
    non-root nodes reachable from the leaves. This value is obtained using the
    :meth:`.get_parent` method. The parent of the root node is the
    :const:`.NULL_NODE` node, :math:`-1`. Similarly, each internal node has a
    pair of children, which are obtained using the :meth:`.get_children`
    method. Each node in the tree has a time associated with it in generations.
    This value is obtained using the :meth:`.get_time` method.

    Sparse trees are not intended to be instantiated directly, and are
    obtained as part of a :class:`.TreeSequence` using the
    :meth:`.trees` method.
    """
    def __init__(self, ll_sparse_tree):
        self._ll_sparse_tree = ll_sparse_tree

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
        :rtype: int
        """
        return self.get_time(self.get_parent(u)) - self.get_time(u)

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

    def get_mrca(self, u, v):
        """
        Returns the most recent common ancestor of the specified nodes.

        :param int u: The first node.
        :param int v: The second node.
        :return: The most recent common ancestor of u and v.
        :rtype: int
        """
        return self._ll_sparse_tree.get_mrca(u, v)

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

    def get_parent(self, u):
        """
        Returns the parent of the specified node. Returns 0 if u is the
        root or is not a node in the current tree.

        :param int u: The node of interest.
        :return: The parent of u.
        :rtype: int
        """
        return self._ll_sparse_tree.get_parent(u)

    def get_children(self, u):
        """
        Returns the children of the specified node as a tuple (v, w). Returns
        the tuple (0, 0) if u is a leaf or is not a node in the current tree.

        :param int u: The node of interest.
        :return: The children of u as a pair of integers
        :rtype: tuple
        """
        return self._ll_sparse_tree.get_children(u)

    def get_time(self, u):
        """
        Returns the time of the specified node in generations. Returns 0 if u
        is a leaf or is not a node in the current tree.

        :param int u: The node of interest.
        :return: The time of u.
        :rtype: float
        """
        return self._ll_sparse_tree.get_time(u)

    def get_population(self, u):
        """
        Returns the population associated with the specified node. For leaf
        nodes this is the population of the sample, and for internal nodes this
        is the population where the corresponding coalescence occured. If the
        specified node is not a member of this tree or population level
        information was not stored in the tree sequence,
        :data:`msprime.NULL_POPULATION` is returned.

        :param int u: The node of interest.
        :return: The ID of the population associated with node u.
        :rtype: int
        """
        return self._ll_sparse_tree.get_population(u)

    def is_internal(self, u):
        """
        Returns True if the specified node is not a leaf.

        :param int u: The node of interest.
        :return: True if u is not a leaf node.
        :rtype: bool
        """
        return not self.is_leaf(u)

    def is_leaf(self, u):
        """
        Returns True if the specified node is a leaf. A node :math:`u` is a
        leaf if :math:`1 \leq u \leq n` for a sample size :math:`n`.

        :param int u: The node of interest.
        :return: True if u is a leaf node.
        :rtype: bool
        """
        return 0 <= u < self.get_sample_size()

    def get_root(self):
        """
        Returns the root of this tree.

        :return: The root node.
        :rtype: int
        """
        return self._ll_sparse_tree.get_root()

    def get_index(self):
        """
        Returns the index this tree occupies in the parent tree sequence.
        This index is zero based, so the first tree in the sequence has index
        0.

        :return: The index of this tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_index()

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

    def get_sample_size(self):
        """
        Returns the sample size for this tree. This is the number of leaf
        nodes in the tree.

        :return: The number of leaf nodes in the tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_sample_size()

    def draw(self, path, width=200, height=200, show_times=False):
        """
        Draws a representation of this tree to the specified path in SVG
        format.

        :param str path: The path to the file to write the SVG.
        :param int width: The width of the image in pixels.
        :param int height: The height of the image in pixels.
        :param bool show_times: If True, show time labels at each internal
            node.
        """
        if not _svgwrite_imported:
            raise ImportError(
                "svgwrite is not installed. try `pip install svgwrite`")
        td = TreeDrawer(self, width, height, show_times)
        td.write(path)

    def get_num_mutations(self):
        """
        Returns the number of mutations on this tree.

        :return: The number of mutations on this tree.
        :rtype: int
        """
        return self._ll_sparse_tree.get_num_mutations()

    def mutations(self):
        """
        Returns an iterator over the mutations on this tree. Each mutation
        is represented as a tuple (position, node), and mutations
        returned in increasing order of position.

        :return: The mutations in this tree.
        :rtype: iterator
        """
        return iter(self._ll_sparse_tree.get_mutations())

    def leaves(self, u):
        """
        Returns an iterator over all the leaves in this tree underneath
        the specified node.

        :param int u: The node of interest.
        :return: An iterator over all leaves in the subtree rooted at u.
        :rtype: iterator
        """
        return _msprime.LeafListIterator(self._ll_sparse_tree, u)

    def get_num_leaves(self, u):
        """
        Returns the number of leaves in this tree underneath the specified
        node. This is a constant time operation.

        :param int u: The node of interest.
        :return: The number of leaves in the subtree rooted at u.
        :rtype: int
        """
        return self._ll_sparse_tree.get_num_leaves(u)

    def get_num_tracked_leaves(self, u):
        """
        Returns the number of leaves in the set specified in the
        ``tracked_leaves`` parameter of the :meth:`msprime.TreeSequence.trees`
        method underneath the specified node. This is a constant time
        operation.

        :param int u: The node of interest.
        :return: The number of leaves within the set of tracked leaves in
            the subtree rooted at u.
        :rtype: int
        """
        return self._ll_sparse_tree.get_num_tracked_leaves(u)

    def _preorder_traversal(self, u):
        stack = [u]
        while len(stack) > 0:
            v = stack.pop()
            if self.is_internal(v):
                stack.extend(reversed(self.get_children(v)))
            yield v

    def nodes(self, root=None, order="preorder"):
        """
        Returns an iterator over the nodes in this tree. If the root parameter
        is provided, iterate over the nodes in the subtree rooted at this
        node. If this is None, iterate over all nodes. If the order parameter
        is provided, iterate over the nodes in required tree traversal order.

        :param int root: The root of the subtree we are traversing.
        :param str order: The traversal ordering. Currently only 'preorder'
            is supported.
        :rtype: iterator
        """
        u = self.get_root() if root is None else root
        if order == "preorder":
            return self._preorder_traversal(u)
        else:
            raise ValueError(
                "Traversal ordering '{}' not supported".format(order))

    def get_parent_dict(self):
        pi = {}
        for j in range(self.get_sample_size()):
            u = j
            while u != NULL_NODE and u not in pi:
                pi[u] = self.get_parent(u)
                u = pi[u]
        return pi

    def get_time_dict(self):
        tau = {}
        for j in range(self.get_sample_size()):
            u = j
            while u != NULL_NODE and u not in tau:
                tau[u] = self.get_time(u)
                u = self.get_parent(u)
        return tau

    def __str__(self):
        return str(self.get_parent_dict())

    def __eq__(self, other):
        return (
            self.get_sample_size() == other.get_sample_size() and
            self.get_parent_dict() == other.get_parent_dict() and
            self.get_time_dict() == other.get_time_dict() and
            self.get_interval() == other.get_interval() and
            self.get_root() == other.get_root() and
            self.get_index() == other.get_index() and
            list(self.mutations()) == list(other.mutations()))

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


def _replicate_generator(sim, rng, mutation_rate, num_replicates):
    """
    Generator function for the many-replicates case of the simulate
    function.
    """
    # Should use range here, but Python 2 makes this awkward...
    j = 0
    while j < num_replicates:
        j += 1
        sim.run()
        tree_sequence = sim.get_tree_sequence()
        tree_sequence.generate_mutations(mutation_rate, rng)
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
        demographic_events=[]):
    """
    Convenience method to create a simulator instance using the same
    parameters as the `simulate` function. Primarily used for testing.
    """
    if sample_size is None and population_configurations is None:
        raise ValueError(
            "Either sample_size or population_configurations must "
            "be specified")
    n = sample_size
    if population_configurations is not None:
        _check_population_configurations(population_configurations)
        n = sum(conf.sample_size for conf in population_configurations)
        if sample_size is not None and n != sample_size:
            raise ValueError(
                "Overall sample_size and the sum of population sample sizes "
                "must be equal")
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
    sim = TreeSimulator(n, recomb_map)
    sim.set_effective_population_size(Ne)
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
        random_seed=None,
        num_replicates=None):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns the resulting :class:`.TreeSequence`.

    :param int sample_size: The number of individuals in our sample.
        If not specified or None, this defaults to the sum of the
        subpopulation sample sizes. Either ``sample_size`` or
        ``population_configurations`` must be specified.
    :param float Ne: The effective population size for the reference
        population. This determines the factor by which the
        per-generation recombination and mutation rates are scaled
        in the simulation. This defaults to 1 if not specified.
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
        demographic_events=demographic_events)
    mu = 0 if mutation_rate is None else mutation_rate
    if num_replicates is None:
        sim.run()
        tree_sequence = sim.get_tree_sequence()
        tree_sequence.generate_mutations(mu, rng)
        return tree_sequence
    else:
        return _replicate_generator(sim, rng, mu, num_replicates)


def load(path):
    """
    Loads a tree sequence from the specified file path. This
    file must be in the HDF5 file format produced by the
    :meth:`msprime.TreeSequence.dump` method.

    :param str path: The file path to write the TreeSequence to.
    :return: The tree sequence object containing the information
        stored in the specified file path.
    :rtype: :class:`msprime.TreeSequence`
    """
    return TreeSequence.load(path)


class TreeSimulator(object):
    """
    Class to simulate trees under the standard neutral coalescent with
    recombination.
    """
    def __init__(self, sample_size, recombination_map):
        if not isinstance(sample_size, int):
            raise TypeError("Sample size must be an integer")
        if sample_size < 2:
            raise ValueError("Sample size must be >= 2")
        if sample_size >= 2**32:
            raise ValueError("sample_size must be < 2**32")
        self._sample_size = sample_size
        if not isinstance(recombination_map, RecombinationMap):
            raise TypeError("RecombinationMap instance required")
        self._recombination_map = recombination_map
        self._random_generator = None
        self._effective_population_size = 1
        self._population_configurations = [
            PopulationConfiguration(sample_size)]
        self._migration_matrix = [[0]]
        self._demographic_events = []
        # Set default block sizes to 64K objects.
        # TODO does this give good performance in a range of scenarios?
        block_size = 64 * 1024
        # We always need at least n segments, so no point in making
        # allocation any smaller than this.
        self._segment_block_size = max(block_size, sample_size)
        self._avl_node_block_size = block_size
        self._node_mapping_block_size = block_size
        self._coalescence_record_block_size = block_size
        # TODO is it useful to bring back the API to set this? Mostly
        # the amount of memory required is tiny.
        self._max_memory = sys.maxsize
        self._ll_sim = None

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
        ll_recombination_rate = self.get_per_locus_scaled_recombination_rate()
        ll_sim = _msprime.Simulator(
            sample_size=self._sample_size,
            random_generator=self._random_generator,
            num_loci=self._recombination_map.get_num_loci(),
            migration_matrix=ll_migration_matrix,
            population_configuration=ll_population_configuration,
            demographic_events=ll_demographic_events,
            scaled_recombination_rate=ll_recombination_rate,
            max_memory=self._max_memory,
            segment_block_size=self._segment_block_size,
            avl_node_block_size=self._avl_node_block_size,
            node_mapping_block_size=self._node_mapping_block_size,
            coalescence_record_block_size=self._coalescence_record_block_size)
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

    def get_tree_sequence(self):
        """
        Returns a TreeSequence representing the state of the simulation.
        """
        ll_tree_sequence = _msprime.TreeSequence()
        ll_recomb_map = self._recombination_map.get_ll_recombination_map()
        Ne = self.get_effective_population_size()
        ll_tree_sequence.create(self._ll_sim, ll_recomb_map, Ne)
        return TreeSequence(ll_tree_sequence)

    def reset(self):
        """
        Resets the simulation so that we can perform another replicate.
        """
        if self._ll_sim is not None:
            self._ll_sim.reset()

    def get_ms_command_line(
            self, executable="ms", num_replicates=1, output_trees=True,
            mutation_rate=None):
        """
        Returns an command line for ms that is equivalent to the parameters
        for this simulator.

        :param str executable: The path to the ms-compatible binary.
        :param int num_relicates: The number of replicates to simulate.
        :param float mutation_rate: The rate of mutation per generation.
        :return: A list of command line arguments that can be used to invoke
            ms with equivalent parameters to this simulator.
        :rtype: list
        """
        L = self._recombination_map.get_length()
        m = self._recombination_map.get_num_loci()
        if m != L:
            raise ValueError(
                "Only recombination maps where L = m are supported")
        r = self._recombination_map.get_total_recombination_rate()
        rho = r
        args = [executable, str(self._sample_size), str(num_replicates)]
        if output_trees:
            args += ["-T"]
        if rho > 0 or L > 1:
            args += ["-r", str(rho), str(L)]
        if mutation_rate is not None:
            scaled_mutation_rate = (
                mutation_rate * 4 * self.get_effective_population_size() * L)
            args += ["-t", str(scaled_mutation_rate)]
        for conf in self._population_configurations:
            if conf.growth_rate > 0:
                args.extend(["-G", str(conf.growth_rate)])
        if len(self._population_configurations) > 1:
            # TODO add -I arguments
            assert False
        for event in self._demographic_events:
            args.extend(event.get_ms_arguments())
        return args


class TreeSequence(object):
    """
    A TreeSequence represents the information generated in a coalescent
    simulation. This includes all the trees across the simulated region,
    along with the mutations (if any are present).
    """

    def __init__(self, ll_tree_sequence):
        self._ll_tree_sequence = ll_tree_sequence

    def get_ll_tree_sequence(self):
        return self._ll_tree_sequence

    @classmethod
    def load(cls, path):
        ts = _msprime.TreeSequence()
        ts.load(path)
        return TreeSequence(ts)

    def get_parameters(self):
        return json.loads(self._ll_tree_sequence.get_simulation_parameters())

    def get_mutations(self):
        # TODO should we provide this???
        return self._ll_tree_sequence.get_mutations()

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

    def get_sample_size(self):
        """
        Returns the sample size for this tree sequence. This is the number
        of leaf nodes in each tree.

        :return: The number of leaf nodes in the tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_sample_size()

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

    def get_num_records(self):
        """
        Returns the number of coalescence records in this tree sequence.
        See the :meth:`.records` method for details on these objects.

        :return: The number of coalescence records defining this tree
            sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_records()

    def get_num_trees(self):
        """
        Returns the number of distinct trees in this tree sequence. This
        is equal to the number of trees returned by the :meth:`.trees`
        method.

        :return: The number of trees in this tree sequence.
        :rtype: int
        """
        count = 0
        for _ in self.trees():
            count += 1
        return count

    def get_num_mutations(self):
        """
        Returns the number of mutations in this tree sequence. See
        the :meth:`msprime.TreeSequence.mutations` method for details on how
        mutations are defined.

        :return: The number of mutations in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_mutations()

    def get_num_nodes(self):
        """
        Returns the number of nodes in this tree sequence. This is the
        largest value :math:`u` such that `u` is a node in any of the
        constituent trees.

        :return: The total number of nodes in this tree sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_nodes()

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
        population ID will be :data:`msprime.NULL_POPULATION`.

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
        of the tree (i.e., closest to the leaves) are returned first.

        :return: An iterator over the diffs between adjacent trees in this
            tree sequence.
        :rtype: iter
        """
        return _msprime.TreeDiffIterator(self._ll_tree_sequence)

    def mutations(self):
        """
        Returns an iterator over the mutations in this tree sequence. Each
        mutation is represented as a tuple (position, node), and mutations
        are returned in increasing order of position.

        :return: The mutations in this tree sequence.
        :rtype: iter
        """
        return iter(self._ll_tree_sequence.get_mutations())

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

    def trees(self, tracked_leaves=[]):
        """
        Returns an iterator over the trees in this tree sequence. Each value
        returned in this iterator is an instance of
        :class:`msprime.SparseTree`.  The ``tracked_leaves`` parameter can be
        used to efficiently count the number of leaves in a given set that
        exist in a particular subtree using the
        :meth:`msprime.SparseTree.get_num_tracked_leaves` method.

        :warning: Do not store the results of this iterator in a list!
           For performance reasons, the same underlying object is used
           for every tree returned which will most likely lead to unexpected
           behaviour.

        :param list tracked_leaves: The list of leaves to be tracked and
            counted using the
            :meth:`msprime.SparseTree.get_num_tracked_leaves` method.
        :return: An iterator over the sparse trees in this tree sequence.
        :rtype: iter
        """
        ll_sparse_tree = _msprime.SparseTree(
            self._ll_tree_sequence, tracked_leaves)
        iterator = _msprime.SparseTreeIterator(
            self._ll_tree_sequence, ll_sparse_tree)
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
        string returned is the haplotype for sample `1`, and so on.

        :return: An iterator over the haplotype strings for the samples in
            this tree sequence.
        :rtype: iter
        """
        return HaplotypeGenerator(self).haplotypes()

    def variants(self):
        return _msprime.VariantGenerator(self._ll_tree_sequence)

    def generate_mutations(self, mutation_rate, random_generator):
        # TODO document this function when it's ready to be brought back
        # into the public interface. We would need to document the
        # RandomGenerator as well.
        self._ll_tree_sequence.generate_mutations(
            mutation_rate, random_generator)

    def set_mutations(self, mutations):
        """
        Sets the mutations in this tree sequence to the specified list of
        `(position, node)` tuples.  Each entry in the list must be a tuple of
        the form :math:`(x, u)`, where :math:`x` is a floating point value
        defining a genomic position and :math:`u` is an integer defining a tree
        node. A genomic position :math:`x` must satisfy :math:`0 \leq x < L`
        where :math:`L` is the sequence length (see
        :meth:`.get_sequence_length`). A node :math:`u` must satisfy
        :math:`0 < u \leq N` where :math:`N` is the largest valued node in
        the tree sequence (see :meth:`.get_num_nodes`).

        :param list mutations: The list of mutations to be assigned to this
            tree sequence.
        """
        self._ll_tree_sequence.set_mutations(mutations)

    def get_pairwise_diversity(self, samples=None):
        """
        Returns the value of pi, the pairwise nucleotide site diversity.
        If `samples` is specified, calculate the diversity within this set.

        Note that we do not check for duplicates within the list of samples,
        and if samples are provided multiple times incorrect results will
        be returned.

        :param iterable samples: The set of samples within which we calculate
            the diversity. If None, calculate diversity within the entire
            sample.
        :return: The pairwise nucleotide site diversity.
        :rtype: float
        """
        if samples is None:
            leaves = list(range(self.get_sample_size()))
        else:
            leaves = list(samples)
        return self._ll_tree_sequence.get_pairwise_diversity(leaves)

    def get_population(self, sample):
        """
        Returns the population ID for the specified sample ID.

        :param int sample: The sample ID of interest.
        :return: The population ID where the specified sample was drawn.
            Returns NULL_POPULATION if no population information is
            available.
        :rtype: int
        """
        if sample < 0 or sample >= self.get_sample_size():
            raise ValueError("Sample ID out of bounds")
        return self._ll_tree_sequence.get_population(sample)

    def get_samples(self, population_id=None):
        """
        Returns the samples matching the specified population ID.

        :param int population_id: The population of interest. If None,
            return all samples.
        :return: The ID of the population we wish to find samples from.
            If None, return samples from all populations.
        :rtype: list
        """
        samples = list(range(self.get_sample_size()))
        if population_id is not None:
            samples = [
                u for u in samples if self.get_population(u) == population_id]
        return samples


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
        from this population. Defaults to 0.
    :param float initial_size: The absolute size of the population at time
        zero. Defaults to the reference population size :math:`N_e`.
    :param float growth_rate: The exponential growth rate of the population
        per generation. Growth rates can be negative. This is zero for a
        constant population size. Defaults to 0.
    """
    def __init__(self, sample_size=0, initial_size=None, growth_rate=0.0):
        self.sample_size = sample_size
        self.initial_size = initial_size
        self.growth_rate = growth_rate

    def get_ll_representation(self, Ne):
        """
        Returns the low-level representation of this PopulationConfiguration.
        """
        initial_size = Ne if self.initial_size is None else self.initial_size
        return {
            "sample_size": self.sample_size,
            "initial_size": initial_size / Ne,
            "growth_rate": self.growth_rate * 4 * Ne
        }


class DemographicEvent(object):
    def __init__(self, type_, time):
        self.type = type_
        self.time = time

    def _get_scaled_time(self, Ne):
        return self.time / (4 * Ne)

    def __str__(self):
        raise NotImplementedError()

    def apply(self, populations, migration_matrix):
        raise NotImplementedError()


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

    def get_ms_arguments(self):
        raise NotImplementedError()

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

    def get_ms_arguments(self):
        raise NotImplemented()

    def __str__(self):
        if self.matrix_index is None:
            ret = "Migration rate change to {} everywhere".format(self.rate)
        else:
            ret = "Migration rate change for {} to {}".format(
                self.matrix_index, self.rate)
        return ret

    def apply(self, populations, migration_matrix):
        if self.matrix_index is None:
            # Change all non-diagonal values.
            for j in range(self._num_populations):
                for k in range(self._num_populations):
                    if j != k:
                        migration_matrix[j][k] = self.rate
        else:
            j, k = self.matrix_index
            migration_matrix[j][k] = self.rate


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

    def get_ms_arguments(self):
        raise NotImplemented()

    def __str__(self):
        return (
            "Mass migration: lineages move from {} to {} with "
            "probability {}".format(
                self.source, self.destination, self.proportion))

    def apply(self, populations, migration_matrix):
        pass


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
        self._demographic_events = collections.defaultdict(list)
        for event in demographic_events:
            self._demographic_events[event.time].append(event)

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
        ll_sim = self._simulator.create_ll_instance()
        N = self._simulator.get_num_populations()
        Ne = self._Ne
        start_time = 0
        scaled_end_time = 0
        while not math.isinf(scaled_end_time):
            events = self._demographic_events[start_time]
            if len(events) > 0:
                print(
                    "Events @ generation {}".format(start_time),
                    file=output)
            for event in events:
                assert event.time == start_time
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
