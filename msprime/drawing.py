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
from __future__ import unicode_literals

import array
import collections

try:
    import svgwrite
    _svgwrite_imported = True
except ImportError:
    _svgwrite_imported = False


def draw_tree(
        tree, width=200, height=200, times=False,
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
        td = SvgTreeDrawer(
                tree, width=width, height=height, show_times=times,
                show_mutation_locations=mutation_locations,
                show_mutation_labels=mutation_labels,
                show_internal_node_labels=internal_node_labels,
                show_leaf_node_labels=leaf_node_labels,
                node_label_text=node_label_text,
                y_padding=20, x_padding=0)
    elif fmt == "ascii":
        td = AsciiTreeDrawer(
                tree, width=width, height=height, show_times=times,
                show_mutation_locations=mutation_locations,
                show_mutation_labels=mutation_labels,
                show_internal_node_labels=internal_node_labels,
                node_label_text=node_label_text,
                show_leaf_node_labels=leaf_node_labels)
    elif fmt == "unicode":
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
            self, tree, width=200, height=200, show_times=False,
            show_mutation_locations=True, show_mutation_labels=False,
            show_internal_node_labels=True, show_leaf_node_labels=True,
            node_label_text=None, x_padding=0, y_padding=0):
        self._width = width
        self._height = height
        self._show_times = show_times
        self._show_mutation_locations = show_mutation_locations
        self._show_mutation_labels = show_mutation_labels
        self._show_internal_node_labels = show_internal_node_labels
        self._show_leaf_node_labels = show_leaf_node_labels
        self._x_scale = width / (tree.get_sample_size() + 2)
        t = tree.get_time(tree.get_root())
        self._y_scale = (height - 2 * y_padding) / t
        self._tree = tree
        self._x_coords = {}
        self._y_coords = {}
        self._node_label_text = {}
        for u in tree.nodes():
            scaled_t = tree.get_time(u) * self._y_scale
            self._y_coords[u] = self._discretise(height - scaled_t - y_padding)
            self._node_label_text[u] = str(u)
        if node_label_text is not None:
            for node, label in node_label_text.items():
                self._node_label_text[node] = label
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
                z = x, self._discretise(y1 + (k + 1) * chunk)
                self._mutations.append((z, mutation))

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
            self._x_coords[node] = self._discretise(a + (b - a) / 2)
        else:
            self._x_coords[node] = self._discretise(self._sample_x * self._x_scale)
            self._sample_x += 1


class SvgTreeDrawer(TreeDrawer):
    """
    Draws trees in SVG format using the svgwrite library.
    """

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
            elif u == self._tree.root:
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
            if u != self._tree.root:
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


class AsciiTreeDrawer(TreeDrawer):
    """
    Draws an ASCII rendering of a tree.
    """
    discretise_coordinates = True

    def draw(self):
        w = self._width
        h = self._height + 1

        # Create a width * height canvas of spaces.
        canvas = bytearray(w * h)
        for row in range(h):
            for col in range(w - 1):
                canvas[row * w + col] = ord(' ')
            canvas[row * w + w - 1] = ord('\n')
        for u in self._tree.nodes():
            col = self._x_coords[u]
            row = self._y_coords[u]
            j = row * w + col
            label = self._node_label_text[u].encode()
            n = len(label)
            canvas[j - n // 2: j + n // 2 + int(n % 2 == 1)] = label
            if self._tree.is_internal(u):
                children = self._tree.children(u)
                row += 1
                left = min(self._x_coords[v] for v in children)
                right = max(self._x_coords[v] for v in children)
                for col in range(left, right):
                    canvas[row * w + col] = ord('-')
                canvas[row * w + self._x_coords[u]] = ord('+')
                top = row + 1
                for v in children:
                    col = self._x_coords[v]
                    for row in range(top, self._y_coords[v] - 1):
                        canvas[row * w + col] = ord('|')
        return str(bytes(canvas).decode())


class UnicodeTreeDrawer(TreeDrawer):
    """
    Draws an Unicode rendering of a tree using box drawing characters.
    """
    discretise_coordinates = True

    def draw(self):
        w = self._width
        h = self._height + 1

        # Create a width * height canvas of spaces.
        canvas = array.array(b'u', " " * (w * h))
        for row in range(h):
            canvas[row * w + w - 1] = '\n'

        for u in self._tree.nodes():
            col = self._x_coords[u]
            row = self._y_coords[u]
            j = row * w + col
            label = array.array(b'u', self._node_label_text[u])
            n = len(label)
            canvas[j - n // 2: j + n // 2 + int(n % 2 == 1)] = label
            if self._tree.is_internal(u):
                children = self._tree.children(u)
                row += 1
                left = min(self._x_coords[v] for v in children)
                right = max(self._x_coords[v] for v in children)
                canvas[row * w + left] = "\u250F"
                canvas[row * w + right] = "\u2513"
                for col in range(left + 1, right):
                    canvas[row * w + col] = "\u2501"
                canvas[row * w + self._x_coords[u]] = "\u2537"
                top = row + 1
                for v in children:
                    col = self._x_coords[v]
                    for row in range(top, self._y_coords[v] - 1):
                        canvas[row * w + col] = "\u2503"

        return canvas.tounicode()
