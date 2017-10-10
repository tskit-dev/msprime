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
        tree, width=None, height=None, times=False,
        mutation_locations=True, mutation_labels=False,
        internal_node_labels=True, leaf_node_labels=True, show_times=None,
        node_label_text=None, format=None):
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
        td = SvgTreeDrawer(
                tree, width=width, height=height, show_times=times,
                show_mutation_locations=mutation_locations,
                show_mutation_labels=mutation_labels,
                show_internal_node_labels=internal_node_labels,
                show_leaf_node_labels=leaf_node_labels,
                node_label_text=node_label_text)
    elif fmt == "ascii":
        td = AsciiTreeDrawer(
                tree, width=width, height=height, show_times=times,
                show_mutation_locations=mutation_locations,
                show_mutation_labels=mutation_labels,
                show_internal_node_labels=internal_node_labels,
                node_label_text=node_label_text,
                show_leaf_node_labels=leaf_node_labels)
    elif fmt == "unicode":
        if IS_PY2:
            raise ValueError("Unicode tree drawing not supported on Python 2")
        td = UnicodeTreeDrawer(
                tree, width=width, height=height, show_times=times,
                show_mutation_locations=mutation_locations,
                show_mutation_labels=mutation_labels,
                show_internal_node_labels=internal_node_labels,
                node_label_text=node_label_text,
                show_leaf_node_labels=leaf_node_labels)
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
            self, tree, width=None, height=None, show_times=False,
            show_mutation_locations=True, show_mutation_labels=False,
            show_internal_node_labels=True, show_leaf_node_labels=True,
            node_label_text=None):
        self._tree = tree
        self._show_times = show_times
        self._show_mutation_locations = show_mutation_locations
        self._show_mutation_labels = show_mutation_labels
        self._show_internal_node_labels = show_internal_node_labels
        self._show_leaf_node_labels = show_leaf_node_labels
        self._num_leaves = len(list(tree.leaves()))
        self._width = width
        self._height = height
        self._x_coords = {}
        self._y_coords = {}
        self._node_label_text = {}
        for u in tree.nodes():
            self._node_label_text[u] = str(u)
        if node_label_text is not None:
            for node, label in node_label_text.items():
                self._node_label_text[node] = label
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
        self._y_scale = (self._height - 2 * y_padding) / t
        for u in self._tree.nodes():
            scaled_t = self._tree.get_time(u) * self._y_scale
            self._y_coords[u] = self._height - scaled_t - y_padding
        self._x_scale = self._width / (self._num_leaves + 2)
        self._sample_x = 1
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
            self._x_coords[node] = self._discretise(self._sample_x * self._x_scale)
            self._sample_x += 1

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
            elif self._tree.parent(u) == NULL_NODE:
                dy = [-5]
            else:
                dx = [-10]
                dy = [-5]
            condition = (
                (self._tree.is_sample(u) and self._show_leaf_node_labels) or
                (self._tree.is_internal(u) and self._show_internal_node_labels))
            if condition:
                labels.add(dwg.text(self._node_label_text[u], x, dx=dx, dy=dy))
            if self._show_times and self._tree.is_internal(u):
                dx[0] += 25
                labels.add(dwg.text(
                    "t = {:.2f}".format(self._tree.get_time(u)), x, dx=dx,
                    dy=dy))
            if self._tree.parent(u) != NULL_NODE:
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


class TextTreeDrawer(TreeDrawer):
    """
    Abstract superclass of TreeDrawers that draw trees in a text buffer.
    """
    discretise_coordinates = True

    array_type = None  # the type used for the array.array canvas
    background_char = None  # The fill char
    eol_char = None  # End of line
    left_shoulder_char = None  # left corner of a horizontal line
    right_shoulder_char = None  # right corner of a horizontal line
    horizontal_line_char = None  # horizontal line fill
    vertical_line_char = None  # vertial line fill
    below_node_char = None  # char just below centre of node label

    def _convert_text(self, text):
        """
        Converts the specified string into an array representation that can be
        filled into the text buffer.
        """
        raise NotImplementedError()

    def _assign_coordinates(self):
        # Get the depth of every node.
        for u in self._tree.nodes():
            v = u
            depth = -1
            while v != NULL_NODE:
                v = self._tree.parent(v)
                depth += 1
            self._y_coords[u] = 3 * depth
        self._height = 0
        if len(self._y_coords) > 0:
            self._height = max(self._y_coords.values()) + 1
        # Get the overall width and assign x coordinates.
        x = 0
        for root in self._tree.roots:
            for u in self._tree.nodes(root, order="postorder"):
                if self._tree.is_leaf(u):
                    self._x_coords[u] = x
                    x += len(self._node_label_text[u]) + 1
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
            label = self._convert_text(self._node_label_text[u])
            n = len(label)
            canvas[j: j + n] = label
            if self._tree.is_internal(u):
                children = self._tree.children(u)
                row += 1
                left = min(self._x_coords[v] for v in children)
                right = max(self._x_coords[v] for v in children)
                canvas[row * w + left] = self.left_shoulder_char
                canvas[row * w + right] = self.right_shoulder_char
                for col in range(left + 1, right):
                    canvas[row * w + col] = self.horizontal_line_char
                if len(self._tree.children(u)) == 1:
                    canvas[row * w + self._x_coords[u]] = self.vertical_line_char
                else:
                    canvas[row * w + self._x_coords[u]] = self.below_node_char
                top = row + 1
                for v in children:
                    col = self._x_coords[v]
                    for row in range(top, self._y_coords[v]):
                        canvas[row * w + col] = self.vertical_line_char
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
    left_shoulder_char = ord('+')
    right_shoulder_char = ord('+')
    horizontal_line_char = ord('-')
    vertical_line_char = ord('|')
    below_node_char = ord('+')

    def _convert_text(self, text):
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
    left_shoulder_char = "\u250F"
    right_shoulder_char = "\u2513"
    horizontal_line_char = "\u2501"
    vertical_line_char = "\u2503"
    below_node_char = "\u2537"

    def _convert_text(self, text):
        return array.array(self.array_type, text)

    def draw(self):
        return self._draw().tounicode()
