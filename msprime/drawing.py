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
Module responsible for visualisations.
"""
from __future__ import division
from __future__ import print_function

import array
import collections
import sys

try:
    import svgwrite
    _svgwrite_imported = True
except ImportError:
    _svgwrite_imported = False

IS_PY2 = sys.version_info[0] < 3

NULL_NODE = -1


def draw_tree(
        tree, width=None, height=None, node_labels=None, node_colours=None,
        mutation_labels=None, mutation_colours=None, format=None):
    # See tree.draw() for documentation on these arguments.
    if format is None:
        format = "SVG"
    fmt = format.lower()
    supported_formats = ["svg", "ascii", "unicode"]
    if fmt not in supported_formats:
        raise ValueError("Unknown format '{}'. Supported formats are {}".format(
            format, supported_formats))
    if fmt == "svg":
        if not _svgwrite_imported:
            raise ImportError(
                "svgwrite is not installed. try `pip install svgwrite`")
        if width is None:
            width = 200
        if height is None:
            height = 200
        cls = SvgTreeDrawer
    elif fmt == "ascii":
        cls = AsciiTreeDrawer
    elif fmt == "unicode":
        if IS_PY2:
            raise ValueError("Unicode tree drawing not supported on Python 2")
        cls = UnicodeTreeDrawer

    td = cls(
        tree, width=width, height=height,
        node_labels=node_labels, node_colours=node_colours,
        mutation_labels=mutation_labels, mutation_colours=mutation_colours)
    return td.draw()


class TreeDrawer(object):
    """
    A class to draw sparse trees in SVG format.
    """

    discretise_coordinates = False

    def _discretise(self, x):
        """
        Discetises the specified value, if necessary.
        """
        ret = x
        if self.discretise_coordinates:
            ret = int(round(x))
        return ret

    def __init__(
            self, tree, width=None, height=None, node_labels=None, node_colours=None,
            mutation_labels=None, mutation_colours=None):
        self._tree = tree
        self._num_leaves = len(list(tree.leaves()))
        self._width = width
        self._height = height
        self._x_coords = {}
        self._y_coords = {}
        self._node_labels = {}
        self._node_colours = {}
        self._mutation_labels = {}
        self._mutation_colours = {}

        # Set the node labels and colours.
        for u in tree.nodes():
            if node_labels is None:
                self._node_labels[u] = str(u)
            else:
                self._node_labels[u] = None
        if node_labels is not None:
            for node, label in node_labels.items():
                self._node_labels[node] = label
        if node_colours is not None:
            for node, colour in node_colours.items():
                self._node_colours[node] = colour

        # Set the mutation labels.
        for site in tree.sites():
            for mutation in site.mutations:
                if mutation_labels is None:
                    self._mutation_labels[mutation.id] = str(mutation.id)
                else:
                    self._mutation_labels[mutation.id] = None
        if mutation_labels is not None:
            for mutation, label in mutation_labels.items():
                self._mutation_labels[mutation] = label
        if mutation_colours is not None:
            for mutation, colour in mutation_colours.items():
                self._mutation_colours[mutation] = colour

        self._assign_coordinates()


class SvgTreeDrawer(TreeDrawer):
    """
    Draws trees in SVG format using the svgwrite library.
    """

    def _assign_coordinates(self):
        y_padding = 20
        t = 1
        if self._tree.num_roots > 0:
            t = max(self._tree.time(root) for root in self._tree.roots)
        # Do we have any mutations over a root?
        mutations_over_root = any(
            self._tree.parent(mut.node) == NULL_NODE for mut in self._tree.mutations())
        root_branch_length = 0
        if mutations_over_root:
            # Allocate a fixed about of space to show the mutations on the
            # 'root branch'
            root_branch_length = self._height / 10
        self._y_scale = (self._height - root_branch_length - 2 * y_padding) / t
        self._y_coords[-1] = y_padding
        for u in self._tree.nodes():
            scaled_t = self._tree.get_time(u) * self._y_scale
            self._y_coords[u] = self._height - scaled_t - y_padding
        self._x_scale = self._width / (self._num_leaves + 2)
        self._leaf_x = 1
        for root in self._tree.roots:
            self._assign_x_coordinates(root)
        self._mutations = []
        node_mutations = collections.defaultdict(list)
        for site in self._tree.sites():
            for mutation in site.mutations:
                node_mutations[mutation.node].append(mutation)
        for child, mutations in node_mutations.items():
            n = len(mutations)
            parent = self._tree.parent(child)
            # Ignore any mutations that are above non-roots that are
            # not in the current tree.
            if child in self._x_coords:
                x = self._x_coords[child]
                y1 = self._y_coords[child]
                y2 = self._y_coords[parent]
                chunk = (y2 - y1) / (n + 1)
                for k, mutation in enumerate(mutations):
                    z = x, self._discretise(y1 + (k + 1) * chunk)
                    self._mutations.append((z, mutation))

    def _assign_x_coordinates(self, node):
        """
        Assign x coordinates to all nodes underneath this node.
        """
        if self._tree.is_internal(node):
            children = self._tree.children(node)
            for c in children:
                self._assign_x_coordinates(c)
            coords = [self._x_coords[c] for c in children]
            a = min(coords)
            b = max(coords)
            self._x_coords[node] = self._discretise(a + (b - a) / 2)
        else:
            self._x_coords[node] = self._discretise(self._leaf_x * self._x_scale)
            self._leaf_x += 1

    def draw(self):
        """
        Writes the SVG description of this tree and returns the resulting XML
        code as text.
        """
        dwg = svgwrite.Drawing(size=(self._width, self._height), debug=True)
        lines = dwg.add(dwg.g(id='lines', stroke='black'))
        left_labels = dwg.add(dwg.g(font_size=14, text_anchor="start"))
        right_labels = dwg.add(dwg.g(font_size=14, text_anchor="end"))
        mid_labels = dwg.add(dwg.g(font_size=14, text_anchor="middle"))
        for u in self._tree.nodes():
            v = self._tree.get_parent(u)
            x = self._x_coords[u], self._y_coords[u]
            colour = "black"
            if self._node_colours.get(u, None) is not None:
                colour = self._node_colours[u]
            dwg.add(dwg.circle(center=x, r=3, fill=colour))
            dx = 0
            dy = -5
            labels = mid_labels
            if self._tree.is_leaf(u):
                dy = 20
            elif self._tree.parent(u) != NULL_NODE:
                dx = 5
                if self._tree.left_sib(u) == NULL_NODE:
                    dx *= -1
                    labels = right_labels
                else:
                    labels = left_labels
            if self._node_labels[u] is not None:
                labels.add(dwg.text(self._node_labels[u], (x[0] + dx, x[1] + dy)))
            if self._tree.parent(u) != NULL_NODE:
                y = self._x_coords[v], self._y_coords[v]
                lines.add(dwg.line(x, (x[0], y[1])))
                lines.add(dwg.line((x[0], y[1]), y))

        # Experimental stuff to render the mutation labels. Not working very
        # well at the moment.
        left_labels = dwg.add(dwg.g(
            font_size=14, text_anchor="start", font_style="italic",
            alignment_baseline="middle"))
        right_labels = dwg.add(dwg.g(
            font_size=14, text_anchor="end", font_style="italic",
            alignment_baseline="middle"))
        for x, mutation in self._mutations:
            r = 3
            colour = "red"
            if self._mutation_colours.get(mutation.id, None) is not None:
                colour = self._mutation_colours[mutation.id]
            dwg.add(dwg.rect(
                insert=(x[0] - r, x[1] - r), size=(2 * r, 2 * r), fill=colour))
            dx = 5
            if self._tree.left_sib(mutation.node) == NULL_NODE:
                dx *= -1
                labels = right_labels
            else:
                labels = left_labels
            if self._mutation_labels[mutation.id] is not None:
                dy = 1.5 * r
                labels.add(dwg.text(
                    self._mutation_labels[mutation.id], (x[0] + dx, x[1] + dy)))
        return dwg.tostring()


class TextTreeDrawer(TreeDrawer):
    """
    Abstract superclass of TreeDrawers that draw trees in a text buffer.
    """
    discretise_coordinates = True

    array_type = None  # the type used for the array.array canvas
    background_char = None  # The fill char
    eol_char = None  # End of line
    left_down_char = None  # left corner of a horizontal line
    right_down_char = None  # right corner of a horizontal line
    horizontal_line_char = None  # horizontal line fill
    vertical_line_char = None  # vertial line fill
    mid_up_char = None  # char in a horizontal line going up
    mid_down_char = None  # char in a horizontal line going down
    mid_up_down_char = None  # char in a horizontal line going down and up

    def _convert_text(self, text):
        """
        Converts the specified string into an array representation that can be
        filled into the text buffer.
        """
        raise NotImplementedError()

    def _assign_coordinates(self):
        # Get the age of each node and rank them.
        times = set(self._tree.time(u) for u in self._tree.nodes())
        depth = {t: 2 * j for j, t in enumerate(sorted(times, reverse=True))}
        for u in self._tree.nodes():
            self._y_coords[u] = depth[self._tree.time(u)]
        self._height = 0
        if len(self._y_coords) > 0:
            self._height = max(self._y_coords.values()) + 1
        # Get the overall width and assign x coordinates.
        x = 0
        for root in self._tree.roots:
            for u in self._tree.nodes(root, order="postorder"):
                if self._tree.is_leaf(u):
                    label_size = 1
                    if self._node_labels[u] is not None:
                        label_size = len(self._node_labels[u])
                    self._x_coords[u] = x
                    x += label_size + 1
                else:
                    coords = [self._x_coords[c] for c in self._tree.children(u)]
                    if len(coords) == 1:
                        self._x_coords[u] = coords[0]
                    else:
                        a = min(coords)
                        b = max(coords)
                        assert b - a > 1
                        self._x_coords[u] = int(round((a + (b - a) / 2)))
            x += 1
        self._width = x + 1

    def _draw(self):
        w = self._width
        h = self._height

        # Create a width * height canvas of spaces.
        canvas = array.array(self.array_type, (w * h) * [self.background_char])
        for u in self._tree.nodes():
            col = self._x_coords[u]
            row = self._y_coords[u]
            j = row * w + col
            label = self._convert_text(self._node_labels[u])
            n = len(label)
            canvas[j: j + n] = label
            if self._tree.is_internal(u):
                children = self._tree.children(u)
                row += 1
                left = min(self._x_coords[v] for v in children)
                right = max(self._x_coords[v] for v in children)
                for col in range(left + 1, right):
                    canvas[row * w + col] = self.horizontal_line_char
                if len(self._tree.children(u)) == 1:
                    canvas[row * w + self._x_coords[u]] = self.vertical_line_char
                else:
                    canvas[row * w + self._x_coords[u]] = self.mid_up_char
                for v in children:
                    col = self._x_coords[v]
                    canvas[row * w + col] = self.mid_down_char
                    if col == self._x_coords[u]:
                        canvas[row * w + col] = self.mid_up_down_char
                    for j in range(row + 1, self._y_coords[v]):
                        canvas[j * w + col] = self.vertical_line_char
                if left == right:
                    canvas[row * w + left] = self.vertical_line_char
                else:
                    canvas[row * w + left] = self.left_down_char
                    canvas[row * w + right] = self.right_down_char

        # Put in the EOLs last so that if we can't overwrite them.
        for row in range(h):
            canvas[row * w + w - 1] = self.eol_char
        return canvas


class AsciiTreeDrawer(TextTreeDrawer):
    """
    Draws an ASCII rendering of a tree.
    """
    array_type = 'b'
    background_char = ord(' ')
    eol_char = ord('\n')
    left_down_char = ord('+')
    right_down_char = ord('+')
    horizontal_line_char = ord('-')
    vertical_line_char = ord('|')
    mid_up_char = ord('+')
    mid_down_char = ord('+')
    mid_up_down_char = ord('+')

    def _convert_text(self, text):
        if text is None:
            text = "|"  # vertical line char
        return array.array(self.array_type, text.encode())

    def draw(self):
        s = self._draw().tostring()
        if not IS_PY2:
            s = s.decode()
        return s


class UnicodeTreeDrawer(TextTreeDrawer):
    """
    Draws an Unicode rendering of a tree using box drawing characters.
    """
    array_type = 'u'
    background_char = ' '
    eol_char = '\n'
    left_down_char = "\u250F"
    right_down_char = "\u2513"
    horizontal_line_char = "\u2501"
    vertical_line_char = "\u2503"
    mid_up_char = "\u253b"
    mid_down_char = "\u2533"
    mid_up_down_char = "\u254b"

    def _convert_text(self, text):
        if text is None:
            text = self.vertical_line_char
        return array.array(self.array_type, text)

    def draw(self):
        return self._draw().tounicode()
