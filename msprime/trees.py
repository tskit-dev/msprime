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
            if v != 0:
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
# - get_total_branch_length
# - Pickle and copy support
class SparseTree(object):
    """
    A SparseTree is a single tree in a :class:`.TreeSequence`. In a sparse tree
    for a sample of size :math:`n`, the leaves are nodes :math:`1` to :math:`n`
    inclusive and internal nodes are integers :math:`> n`. The value of these
    nodes is strictly increasing as we ascend the tree and the root of the tree
    is the node with the largest value that is reachable from  the leaves.
    Each node in the tree has a parent, which is non-zero for all non-root
    nodes reachable from the leaves. This value is obtained using the
    :meth:`.get_parent` method. The parent of the root node is 0. Similarly,
    each internal node has a pair of children, which are obtained using
    the :meth:`.get_children` method. Each node in the tree has a time
    associated with it in coalescent time units. This value is obtained using
    the :meth:`.get_time` method.

    Sparse trees are not intended to be instantiated directly, and are
    obtained as part of a :class:`.TreeSequence` using the
    :meth:`.trees` method.
    """
    def __init__(self, ll_sparse_tree):
        self._ll_sparse_tree = ll_sparse_tree

    def get_branch_length(self, u):
        """
        Returns the length of the branch (in time units) joining the
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
        Returns the time of the specified node. Returns 0 if u is a leaf
        or is not a node in the current tree.

        :param int u: The node of interest.
        :return: The time of u.
        :rtype: float
        """
        return self._ll_sparse_tree.get_time(u)

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
        return 1 <= u <= self.get_sample_size()

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
        for j in range(1, self.get_sample_size() + 1):
            u = j
            while u != 0 and u not in pi:
                pi[u] = self.get_parent(u)
                u = pi[u]
        return pi

    def get_time_dict(self):
        tau = {}
        for j in range(1, self.get_sample_size() + 1):
            u = j
            while u != 0 and u not in tau:
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


def simulate(
        sample_size, num_loci=1, scaled_recombination_rate=0.0,
        scaled_mutation_rate=None,
        population_models=[], random_seed=None, max_memory=None):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns the resulting :class:`.TreeSequence`.

    :param int sample_size: The number of individuals in our sample.
    :param int num_loci: The length of the simulated region in
        discrete non-recombining loci.
    :param float scaled_recombination_rate: The rate of recombination
        between adjacent loci per :math:`4N_0` generations.
    :param float scaled_mutation_rate: The rate of mutation
        per locus per :math:`4N_0` generations.
    :param list population_models: The list of :class:`.PopulationModel`
        instances describing the demographic history of the population.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated.
    :param int,str max_memory: The maximum amount of memory used
        during the simulation. If this is exceeded, the simulation will
        terminate with a :class:`LibraryError` exception. By default
        this is not set, and no memory limits are imposed.
    :return: The :class:`.TreeSequence` object representing the results
        of the simulation.
    :rtype: :class:`.TreeSequence`
    """
    sim = TreeSimulator(sample_size)
    sim.set_num_loci(num_loci)
    sim.set_scaled_recombination_rate(scaled_recombination_rate)
    sim.set_random_seed(random_seed)
    sim.set_max_memory(max_memory)
    # Reinterpret the population models in terms of the new interface.
    # This will be deprecated in later releases in favour of a direct
    # approach.
    initial_size = 1
    initial_growth_rate = 0
    demographic_events = []
    for model in population_models:
        if isinstance(model, ConstantPopulationModel):
            if model.start_time == 0:
                initial_size = model.size
            else:
                demographic_events.append(SizeChangeEvent(
                    model.start_time, model.size))
        elif isinstance(model, ExponentialPopulationModel):
            if model.start_time == 0:
                initial_growth_rate = model.alpha
            else:
                demographic_events.append(GrowthRateChangeEvent(
                    model.start_time, model.alpha))
        else:
            raise TypeError("Arguments must be PopulationModel instances")
    sim.set_population_configurations([PopulationConfiguration(
        sample_size, initial_size, initial_growth_rate)])
    sim.set_demographic_events(demographic_events)
    sim.run()
    tree_sequence = sim.get_tree_sequence()
    if scaled_mutation_rate is not None:
        tree_sequence.generate_mutations(scaled_mutation_rate, random_seed)
    return tree_sequence


def simulate_tree(
        sample_size, scaled_mutation_rate=None, population_models=[],
        random_seed=None, max_memory=None):
    """
    Simulates the coalescent at a single locus for the specified sample size
    under the specified list of population models. Returns a
    :class:`.SparseTree` representing the results.

    :param int sample_size: The number of individuals in our sample.
    :param float scaled_mutation_rate: The rate of mutation
        per :math:`4N_0` generations.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated.
    :param list population_models: The list of :class:`.PopulationModel`
        instances describing the demographic history of the population.
    :param int,str max_memory: The maximum amount of memory used
        during the simulation. If this is exceeded, the simulation will
        terminate with a :class:`LibraryError` exception. By default
        this is not set, and no memory limits are imposed.
    :return: The :class:`.SparseTree` object representing the results
        of the simulation.
    :rtype: :class:`.SparseTree`
    """
    tree_sequence = simulate(
        sample_size, scaled_mutation_rate=scaled_mutation_rate,
        population_models=population_models, random_seed=random_seed,
        max_memory=max_memory)
    return next(tree_sequence.trees())


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
    def __init__(self, sample_size):
        if not isinstance(sample_size, int):
            raise TypeError("Sample size must be an integer")
        if sample_size < 2:
            raise ValueError("Sample size must be >= 2")
        if sample_size >= 2**32:
            raise ValueError("sample_size must be < 2**32")
        self._sample_size = sample_size
        self._scaled_recombination_rate = 0.0
        self._num_loci = 1
        self._migration_matrix = None
        self._population_configurations = None
        self._demographic_events = []
        self._random_seed = None
        self._segment_block_size = None
        self._avl_node_block_size = None
        self._node_mapping_block_size = None
        self._coalescence_record_block_size = None
        self._max_memory = None
        self._ll_sim = None

    def get_sample_size(self):
        return self._sample_size

    def get_scaled_recombination_rate(self):
        return self._scaled_recombination_rate

    def get_migration_matrix(self):
        return self._migration_matrix

    def get_num_loci(self):
        return self._num_loci

    def get_random_seed(self):
        return self._random_seed

    def get_population_configurations(self):
        return self._population_configurations

    def get_sample_configuration(self):
        return [conf.sample_size for conf in self._population_configurations]

    def get_demographic_events(self):
        return self._demographic_events

    def get_num_breakpoints(self):
        return self._ll_sim.get_num_breakpoints()

    def get_breakpoints(self):
        return self._ll_sim.get_breakpoints()

    def get_used_memory(self):
        return self._ll_sim.get_used_memory()

    def get_time(self):
        return self._ll_sim.get_time()

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

    def set_num_loci(self, num_loci):
        if not isinstance(num_loci, int):
            raise TypeError("num_loci must be an integer")
        if num_loci < 1:
            raise ValueError("Positive number of loci required")
        if num_loci >= 2**32:
            raise ValueError("num_loci must be < 2**32")
        self._num_loci = num_loci

    def set_scaled_recombination_rate(self, scaled_recombination_rate):
        if scaled_recombination_rate < 0:
            raise ValueError("Recombination rate cannot be negative")
        self._scaled_recombination_rate = scaled_recombination_rate

    def set_migration_matrix(self, migration_matrix):
        self._migration_matrix = migration_matrix

    def set_population_configurations(self, population_configurations):
        self._population_configurations = population_configurations

    def set_demographic_events(self, demographic_events):
        self._demographic_events = demographic_events

    def set_random_seed(self, random_seed):
        self._random_seed = random_seed

    def set_segment_block_size(self, segment_block_size):
        self._segment_block_size = segment_block_size

    def set_avl_node_block_size(self, avl_node_block_size):
        self._avl_node_block_size = avl_node_block_size

    def set_node_mapping_block_size(self, node_mapping_block_size):
        self._node_mapping_block_size = node_mapping_block_size

    def set_coalescence_record_block_size(self, coalescence_record_block_size):
        self._coalescence_record_block_size = coalescence_record_block_size

    def set_max_memory(self, max_memory):
        """
        Sets the approximate maximum memory used by the simulation
        to the specified value.  This can be suffixed with
        K, M or G to specify units of Kibibytes, Mibibytes or Gibibytes.
        """
        if max_memory is None:
            self._max_memory = max_memory
        else:
            s = max_memory
            d = {"K": 2**10, "M": 2**20, "G": 2**30}
            multiplier = 1
            value = s
            if s.endswith(tuple(d.keys())):
                value = s[:-1]
                multiplier = d[s[-1]]
            n = int(value)
            self._max_memory = n * multiplier

    def _set_environment_defaults(self):
        """
        Sets sensible default values for the memory usage parameters.
        """
        # Set the block sizes using our estimates.
        n = self._sample_size
        m = self._num_loci
        if self._population_configurations is None:
            self._population_configurations = [PopulationConfiguration(n)]
        N = len(self._population_configurations)
        if self._migration_matrix is None:
            self._migration_matrix = [
                [0.0 for j in range(N)] for k in range(N)]
        if self._demographic_events is None:
            self._demographic_events = []
        rho = 4 * self._scaled_recombination_rate * (m - 1)
        num_trees = min(m // 2, rho * harmonic_number(n - 1))
        b = 10  # Baseline maximum
        num_trees = max(b, int(num_trees))
        num_avl_nodes = max(b, 4 * n + num_trees)
        # TODO This is total guesswork. We need to plot this for a range
        # of values and see what a good approximation is.
        num_segments = max(b, int(math.log(n) * rho))
        if self._avl_node_block_size is None:
            self._avl_node_block_size = num_avl_nodes
        if self._segment_block_size is None:
            self._segment_block_size = num_segments
        if self._node_mapping_block_size is None:
            self._node_mapping_block_size = num_trees
        if self._coalescence_record_block_size is None:
            memory = 16 * 2**10  # 16M
            # Each coalescence record is 32bytes
            self._coalescence_record_block_size = memory // 32
        if self._random_seed is None:
            self._random_seed = random.randint(0, 2**31 - 1)
        if self._max_memory is None:
            self._max_memory = sys.maxsize  # Unlimited

    def run(self):
        """
        Runs the simulation until complete coalescence has occurred.
        """
        assert self._ll_sim is None
        self._set_environment_defaults()
        # We flatten the migration matrix for the low-level interface.
        N = len(self._population_configurations)
        ll_migration_matrix = [0 for j in range(N**2)]
        for j in range(N):
            for k in range(N):
                ll_migration_matrix[j * N + k] = self._migration_matrix[j][k]
        ll_population_configuration = [
            conf.get_ll_representation()
            for conf in self._population_configurations]
        ll_demographic_events = [
            event.get_ll_representation(N)
            for event in self._demographic_events]
        self._ll_sim = _msprime.Simulator(
            sample_size=self._sample_size,
            num_loci=self._num_loci,
            migration_matrix=ll_migration_matrix,
            population_configuration=ll_population_configuration,
            demographic_events=ll_demographic_events,
            scaled_recombination_rate=self._scaled_recombination_rate,
            random_seed=self._random_seed,
            max_memory=self._max_memory,
            segment_block_size=self._segment_block_size,
            avl_node_block_size=self._avl_node_block_size,
            node_mapping_block_size=self._node_mapping_block_size,
            coalescence_record_block_size=self._coalescence_record_block_size)
        self._ll_sim.run()

    def get_tree_sequence(self):
        """
        Returns a TreeSequence representing the state of the simulation.
        """
        ll_tree_sequence = _msprime.TreeSequence()
        ll_tree_sequence.create(self._ll_sim)
        return TreeSequence(ll_tree_sequence)

    def reset(self):
        """
        Resets the simulation so that we can perform another replicate.
        """
        self._ll_sim = None

    def get_ms_command_line(
            self, executable="ms", num_replicates=1, output_trees=True,
            scaled_mutation_rate=None):
        """
        Returns an command line for ms that is equivalent to the parameters
        for this simulator.

        :param str executable: The path to the ms-compatible binary.
        :param int num_relicates: The number of replicates to simulate.
        :param float scaled_mutation_rate: The rate of mutation
            per :math:`4N_0` generations.
        :return: A list of command line arguments that can be used to invoke
            ms with equivalent parameters to this simulator.
        :rtype: list
        """
        self._set_environment_defaults()
        m = self._num_loci
        rho = self.get_scaled_recombination_rate() * (m - 1)
        args = [executable, str(self._sample_size), str(num_replicates)]
        if output_trees:
            args += ["-T"]
        if m > 1:
            args += ["-r", str(rho), str(m)]
        if scaled_mutation_rate is not None:
            mu = scaled_mutation_rate * m
            args += ["-t", str(mu)]
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

    def newick_trees(self, precision=3, breakpoints=None):
        # TODO document this method.
        iterator = _msprime.NewickConverter(self._ll_tree_sequence, precision)
        if breakpoints is None:
            for length, tree in iterator:
                yield length, tree
        else:
            trees_covered = 0
            j = 0
            # TODO this is ugly. Update the alg so we don't need this
            # bracketing.
            bp = [0] + breakpoints + [self.get_num_loci()]
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

    def get_num_loci(self):
        """
        Returns the number of loci in this tree sequence. This defined the
        genomic scaler over which tree coordinates are defined. Given a
        tree sequence with :math:`m` loci, the constituent trees will be
        defined over the half-closed interval :math:`(0, m])`. Each tree
        then covers some subset of this interval --- see
        :meth:`msprime.SparseTree.get_interval` for details.

        :return: The number of discrete non-recombining loci in this tree
            sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_loci()

    def get_num_records(self):
        """
        Returns the number of coalescence records in this tree sequence.
        See the :meth:`.records` method for details on these objects.

        :return: The number of coalescence records defining this tree
            sequence.
        :rtype: int
        """
        return self._ll_tree_sequence.get_num_records()

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
        :math:`(l, r, u, c, t)` defining the assignment of a tree node
        across an interval. The range of this record is as the half-open
        genomic interval :math:`[l, r)`, such that it applies to all
        positions :math:`l \leq x < r`. Each record represents the
        assignment of a pair of children :math:`c` to a parent
        :math:`u`. This assignment happens at time :math:`t` in
        coalescent units.

        :return: An iterator of all :math:`(l, r, u, c, t)` tuples defining
            the coalescence records in this tree sequence.
        :rtype: iter
        """
        for j in range(self.get_num_records()):
            yield self._ll_tree_sequence.get_record(j)

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
        returned in increasing order of position.

        :return: The mutations in this tree sequence.
        :rtype: iter
        """
        return iter(self._ll_tree_sequence.get_mutations())

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

    def generate_mutations(self, scaled_mutation_rate, random_seed=None):
        """
        Generates mutation according to the infinite sites model. This
        method over-writes any existing mutations stored in the tree
        sequence.

        :param float scaled_mutation_rate: The rate of mutation
            per locus per :math:`4N_0` generations.
        :param int random_seed: The random seed to use when generating
            mutations.
        """
        seed = random_seed
        if random_seed is None:
            seed = random.randint(0, 2**31)
        self._ll_tree_sequence.generate_mutations(scaled_mutation_rate, seed)

    def set_mutations(self, mutations):
        """
        Sets the mutations in this tree sequence to the specified list of
        `(position, node)` tuples.  Each entry in the list must be a tuple of
        the form :math:`(x, u)`, where :math:`x` is a floating point value
        defining a genomic position and :math:`u` is an integer defining a tree
        node. A genomic position :math:`x` must satisfy :math:`0 \leq x < m`
        where :math:`m` is the number of loci (see :meth:`.get_num_loci`). A
        node :math:`u` must satisfy :math:`0 < u \leq N` where :math:`N` is the
        largest valued node in the tree sequence (see :meth:`.get_num_nodes`).

        :param list mutations: The list of mutations to be assigned to this
            tree sequence.
        """
        self._ll_tree_sequence.set_mutations(mutations)

    def get_pairwise_diversity(self):
        """
        Returns the value of pi, the pairwise nucleotide site diversity.

        :return: The pairwise nucleotide site diversity.
        :rtype: iter
        """
        pi = 0
        n = self.get_sample_size()
        denom = n * (n - 1) / 2
        for t in self.trees():
            for _, node in t.mutations():
                k = t.get_num_leaves(node)
                pi += k * (n - k) / denom
        return pi


class HaplotypeGenerator(object):

    def __init__(self, tree_sequence):
        self._tree_sequence = tree_sequence
        ts = self._tree_sequence.get_ll_tree_sequence()
        self._ll_haplotype_generator = _msprime.HaplotypeGenerator(ts)

    def get_haplotype(self, sample_id):
        return self._ll_haplotype_generator.get_haplotype(sample_id)

    def haplotypes(self):
        for j in range(1, self._tree_sequence.get_sample_size() + 1):
            yield self.get_haplotype(j)


class PopulationConfiguration(object):
    """
    TODO document.
    """
    def __init__(self, sample_size=0, initial_size=1.0, growth_rate=0.0):
        self.sample_size = sample_size
        self.initial_size = initial_size
        self.growth_rate = growth_rate

    def get_ll_representation(self):
        """
        Returns the low-level representation of this PopulationConfiguration.
        """
        return self.__dict__


class DemographicEvent(object):
    def __init__(self, type_, time):
        self.type = type_
        self.time = time


class GrowthRateChangeEvent(DemographicEvent):
    def __init__(self, time, growth_rate, population_id=-1):
        super(GrowthRateChangeEvent, self).__init__(
            "growth_rate_change", time)
        self.time = time
        self.growth_rate = growth_rate
        self.population_id = population_id

    def get_ll_representation(self, num_populations):
        return {
            "type": self.type,
            "time": self.time,
            "growth_rate": self.growth_rate,
            "population_id": self.population_id
        }

    def get_ms_arguments(self):
        if self.population_id == -1:
            return ["-eG", str(self.time), str(self.growth_rate)]
        else:
            return [
                "-eg", str(self.time), str(self.population_id + 1),
                str(self.growth_rate)]


class SizeChangeEvent(DemographicEvent):
    def __init__(self, time, size, population_id=-1):
        super(SizeChangeEvent, self).__init__("size_change", time)
        self.size = size
        self.population_id = population_id

    def get_ll_representation(self, num_populations):
        return {
            "type": self.type,
            "time": self.time,
            "size": self.size,
            "population_id": self.population_id
        }

    def get_ms_arguments(self):
        if self.population_id == -1:
            return ["-eN", str(self.time), str(self.size)]
        else:
            return [
                "-en", str(self.time), str(self.population_id + 1),
                str(self.growth_rate)]


class MigrationRateChangeEvent(DemographicEvent):
    def __init__(self, time, rate, matrix_index=None):
        super(MigrationRateChangeEvent, self).__init__(
            "migration_rate_change", time)
        self.rate = rate
        self.matrix_index = matrix_index

    def _convert_matrix_index(self, matrix_index):
        ret = -1
        if matrix_index is not None:
            ret = matrix_index[0] * self._num_populations + matrix_index[1]
        return ret

    def get_ll_representation(self, num_populations):
        matrix_index = -1
        if self.matrix_index is not None:
            matrix_index = (
                self.matrix_index[0] * num_populations + self.matrix_index[1])
        return {
            "type": self.type,
            "time": self.time,
            "migration_rate": self.rate,
            "matrix_index": matrix_index
        }

    def get_ms_arguments(self):
        raise NotImplemented()


class MassMigrationEvent(DemographicEvent):
    def __init__(self, time, source, destination, proportion):
        super(MassMigrationEvent, self).__init__("mass_migration", time)
        self.source = source
        self.destination = destination
        self.proportion = proportion

    def get_ll_representation(self, num_populations):
        return {
            "type": self.type,
            "time": self.time,
            "source": self.source,
            "destination": self.destination,
            "proportion": self.proportion
        }

    def get_ms_arguments(self):
        raise NotImplemented()


##############################
#
# Deprecated interface for specifying population demography
#
##############################


class PopulationModel(object):
    """
    Superclass of simulation population models.
    """
    def __init__(self, start_time):
        self.start_time = start_time


class ConstantPopulationModel(PopulationModel):
    """
    A population model in which the size of the population
    is a fixed multiple ``size`` of :math:`N_0` (the Wright-Fisher population
    size at time 0), which starts at the specified time.

    :param float start_time: The time (in coalescent units) at which
        this population model begins.
    :param float size: The size of the population under this model
        relative to :math:`N_0`.
    """
    def __init__(self, start_time, size):
        super(ConstantPopulationModel, self).__init__(start_time)
        self.size = size


class ExponentialPopulationModel(PopulationModel):
    """
    A population model in which the size is exponentially growing (or
    shrinking). If we have a :attr:`start_time` of :math:`s`, the population
    size at a time :math:`t` (measured in units of :math:`4N_0` generations)
    is :math:`N_s e^{\\alpha (t - s)}`, where :math:`N_s` is the population
    size at time :math:`s`.

    :param float start_time: The time (in coalescent units) at which
        this population model begins.
    :param float alpha: The exponential growth (or contraction) rate
        :math:`\\alpha`.
    """
    def __init__(self, start_time, alpha):
        super(ExponentialPopulationModel, self).__init__(start_time)
        self.alpha = alpha


def harmonic_number(n):
    """
    Returns the nth Harmonic number.
    """
    EulerGamma = 0.5772156649
    return math.log(n) + EulerGamma
