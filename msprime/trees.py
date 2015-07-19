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

import math
import random

import svgwrite

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
    :meth:`.sparse_trees` method.
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
        :rtype: iter
        """
        return iter(self._ll_sparse_tree.get_mutations())

    def leaves(self, u):
        """
        Returns an iterator over all the leaves in this tree underneath
        the specified node.

        :param int u: The node of interest.
        :return: An iterator over all leaves in the subtree rooted at u.
        :rtype: iter
        """
        stack = [u]
        while len(stack) > 0:
            v = stack.pop()
            if self.is_leaf(v):
                yield v
            else:
                stack.extend(self.get_children(v))

    def nodes(self):
        return self.get_parent_dict().keys()

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
            list(self.mutations()) == list(other.mutations()))

    def __ne__(self, other):
        return not self.__eq__(other)


def simulate(
        sample_size, num_loci=1, scaled_recombination_rate=0.0,
        scaled_mutation_rate=None,
        population_models=[], random_seed=None, max_memory="1G"):
    """
    Simulates the coalescent with recombination under the specified model
    parameters and returns the resulting :class:`.TreeSequence`.

    **TODO** concise description of the model parameters and how we
    can run the simulations we are interested in.

    :param int sample_size: The number of individuals in our sample.
    :param int num_loci: The length of the simulated region in
        bases.
    :param float scaled_recombination_rate: The rate of recombination
        between adjacent bases per :math:`4N` generations.
    :param float scaled_mutation_rate: The rate of mutation
        per site per :math:`4N` generations.
    :param list population_models: The list of :class:`.PopulationModel`
        instances describing the demographic history of the population.
    :param int random_seed: The random seed. If this is `None`, a
        random seed will be automatically generated.
    :param int,str max_memory: The maximum amount of memory used
        during the simulation. If this is exceeded, the simulation will
        terminate with a :class:`LibraryError` exception.
    :return: The :class:`.TreeSequence` object representing the results
        of the simulation.
    :rtype: :class:`.TreeSequence`
    """
    sim = TreeSimulator(sample_size)
    sim.set_num_loci(num_loci)
    sim.set_scaled_recombination_rate(scaled_recombination_rate)
    sim.set_random_seed(random_seed)
    sim.set_max_memory(max_memory)
    for m in population_models:
        sim.add_population_model(m)
    sim.run()
    tree_sequence = sim.get_tree_sequence()
    if scaled_mutation_rate is not None:
        tree_sequence.generate_mutations(scaled_mutation_rate, random_seed)
    return tree_sequence


def simulate_tree(
        sample_size, scaled_mutation_rate=0, population_models=[],
        random_seed=None, max_memory="10M"):
    """
    Simulates the coalescent at a single locus for the specified sample size
    under the specified list of population models. Returns a SparseTree
    representing the results.
    """
    tree_sequence = simulate(
        sample_size, population_models=population_models,
        random_seed=random_seed, max_memory=max_memory)
    tree_sequence.generate_mutations(scaled_mutation_rate, random_seed)
    return next(tree_sequence.sparse_trees())


def load(path):
    """
    Loads a tree sequence file from the specified path.
    """
    return TreeSequence.load(path)


class TreeSimulator(object):
    """
    Class to simulate trees under the standard neutral coalescent with
    recombination.
    """
    def __init__(self, sample_size):
        self._sample_size = sample_size
        self._scaled_recombination_rate = 1.0
        self._num_loci = 1
        self._population_models = []
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

    def get_num_loci(self):
        return self._num_loci

    def get_random_seed(self):
        return self._random_seed

    def get_population_models(self):
        return self._population_models

    def get_num_breakpoints(self):
        return self._ll_sim.get_num_breakpoints()

    def get_breakpoints(self):
        return self._ll_sim.get_breakpoints()

    def get_used_memory(self):
        return self._ll_sim.get_used_memory()

    def get_time(self):
        return self._ll_sim.get_time()

    def get_num_avl_node_blocks(self):
        return self._ll_sim.get_num_avl_node_blocks()

    def get_num_coalescence_record_blocks(self):
        return self._ll_sim.get_num_coalescence_record_blocks()

    def get_num_node_mapping_blocks(self):
        return self._ll_sim.get_num_node_mapping_blocks()

    def get_num_segment_blocks(self):
        return self._ll_sim.get_num_segment_blocks()

    def get_num_coancestry_events(self):
        return self._ll_sim.get_num_coancestry_events()

    def get_num_recombination_events(self):
        return self._ll_sim.get_num_recombination_events()

    def get_num_multiple_recombination_events(self):
        return self._ll_sim.get_num_multiple_recombination_events()

    def get_max_memory(self):
        return self._ll_sim.get_max_memory()

    def add_population_model(self, pop_model):
        self._population_models.append(pop_model)

    def set_num_loci(self, num_loci):
        self._num_loci = num_loci

    def set_scaled_recombination_rate(self, scaled_recombination_rate):
        self._scaled_recombination_rate = scaled_recombination_rate

    def set_effective_population_size(self, effective_population_size):
        self._effective_population_size = effective_population_size

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
        # First check to make sure they are sane.
        if not isinstance(n, int):
            raise TypeError("Sample size must be an integer")
        if not isinstance(m, int):
            raise TypeError("Number of loci must be an integer")
        if n < 2:
            raise ValueError("Sample size must be >= 2")
        if m < 1:
            raise ValueError("Postive number of loci required")
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
            self._max_memory = 10 * 1024 * 1024  # 10MiB by default

    def run(self):
        """
        Runs the simulation until complete coalescence has occured.
        """
        # Sort the models by start time
        models = sorted(self._population_models, key=lambda m: m.start_time)
        models = [m.get_ll_model() for m in models]
        assert self._ll_sim is None
        self._set_environment_defaults()
        self._ll_sim = _msprime.Simulator(
            sample_size=self._sample_size,
            num_loci=self._num_loci,
            population_models=models,
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

    @classmethod
    def load(cls, path):
        ts = _msprime.TreeSequence()
        ts.load(path)
        return TreeSequence(ts)

    def get_sample_size(self):
        return self._ll_tree_sequence.get_sample_size()

    def get_num_loci(self):
        return self._ll_tree_sequence.get_num_loci()

    def get_num_records(self):
        return self._ll_tree_sequence.get_num_records()

    def get_num_mutations(self):
        return self._ll_tree_sequence.get_num_mutations()

    def get_num_nodes(self):
        return self._ll_tree_sequence.get_num_nodes()

    def get_mutations(self):
        return self._ll_tree_sequence.get_mutations()

    def records(self):
        for j in range(self.get_num_records()):
            yield self._ll_tree_sequence.get_record(j)

    def diffs(self):
        return _msprime.TreeDiffIterator(self._ll_tree_sequence)

    def mutations(self):
        """
        Returns an iterator over the mutations on this tree. Each mutation
        is represented as a tuple (position, node), and mutations
        returned in increasing order of position.

        :return: The mutations in this tree.
        :rtype: iter
        """
        return iter(self._ll_tree_sequence.get_mutations())

    def sparse_trees(self):
        """
        Returns an iterator over the sparse trees in this tree
        sequence.
        """
        ll_sparse_tree = _msprime.SparseTree(self._ll_tree_sequence)
        iterator = _msprime.SparseTreeIterator(
            self._ll_tree_sequence, ll_sparse_tree)
        sparse_tree = SparseTree(ll_sparse_tree)
        for _ in iterator:
            yield sparse_tree

    def newick_trees(self, precision=3, breakpoints=None):
        iterator = _msprime.NewickConverter(self._ll_tree_sequence, precision)
        if breakpoints is None:
            for length, tree in iterator:
                yield length, tree
        else:
            trees_covered = 0
            j = 0
            for length, tree in iterator:
                trees_covered += length
                while breakpoints[j] < trees_covered:
                    j += 1
                    yield breakpoints[j] - breakpoints[j - 1], tree

    def generate_mutations(self, scaled_mutation_rate, random_seed=None):
        """
        Generates mutation according to the infinite sites model. This
        method over-writes any existing mutations stored in the tree
        sequence.
        """
        seed = random_seed
        if random_seed is None:
            seed = random.randint(0, 2**31)
        self._ll_tree_sequence.generate_mutations(scaled_mutation_rate, seed)

    def haplotypes(self):
        """
        Returns an iterator over the haplotypes resulting from the trees
        and mutations in this tree sequence as a string of '1's and '0's.
        """
        return HaplotypeGenerator(self).haplotypes()


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


class PopulationModel(object):
    """
    Superclass of simulation population models.
    """
    def __init__(self, start_time):
        self.start_time = start_time

    def get_ll_model(self):
        """
        Returns the low-level model corresponding to this population
        model.
        """
        return self.__dict__


class ConstantPopulationModel(PopulationModel):
    """
    A population model in which the size of the population
    is a fixed multiple ``size`` of the size at time 0, which
    starts at the specified time.

    :param float start_time: The time (in coalescent units) at which
        this population model begins.
    :param float size: The size of the population under this model
        relative to its size at time 0.
    """
    def __init__(self, start_time, size):
        super(ConstantPopulationModel, self).__init__(start_time)
        self.size = size
        self.type = _msprime.POP_MODEL_CONSTANT


class ExponentialPopulationModel(PopulationModel):
    """
    Class representing an exponentially growing or shrinking population.
    TODO document model.
    """
    def __init__(self, start_time, alpha):
        super(ExponentialPopulationModel, self).__init__(start_time)
        self.alpha = alpha
        self.type = _msprime.POP_MODEL_EXPONENTIAL


def harmonic_number(n):
    """
    Returns the nth Harmonic number.
    """
    EulerGamma = 0.5772156649
    return math.log(n) + EulerGamma
