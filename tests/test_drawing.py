# -*- coding: utf-8 -*-
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
Test cases for visualisation in msprime.
"""
from __future__ import print_function
from __future__ import division

import os
import sys
import tempfile
import unittest
import xml.etree

import six
import msprime
import tests.tsutil as tsutil

IS_PY2 = sys.version_info[0] < 3


class TestTreeDraw(unittest.TestCase):
    """
    Tests for the tree drawing functionality.
    """
    def get_binary_tree(self):
        ts = msprime.simulate(10, random_seed=1, mutation_rate=1)
        return next(ts.trees())

    def get_nonbinary_tree(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=0.1, population=0, proportion=0.5)]
        ts = msprime.simulate(
            10, recombination_rate=5, mutation_rate=10,
            demographic_events=demographic_events, random_seed=1)
        for t in ts.trees():
            for u in t.nodes():
                if len(t.children(u)) > 2:
                    return t
        assert False

    def get_multiroot_tree(self):
        ts = msprime.simulate(15, random_seed=1)
        # Take off the top quarter of edges
        tables = ts.dump_tables()
        edges = tables.edges
        n = len(edges) - len(edges) // 4
        edges.set_columns(
            left=edges.left[:n], right=edges.right[:n],
            parent=edges.parent[:n], child=edges.child[:n])
        ts = msprime.load_tables(nodes=tables.nodes, edges=edges)
        for t in ts.trees():
            if t.num_roots > 1:
                return t
        assert False

    def get_mutations_over_roots_tree(self):
        ts = msprime.simulate(15, random_seed=1)
        ts = tsutil.decapitate(ts, 20)
        tables = ts.tables
        sites = msprime.SiteTable()
        mutations = msprime.MutationTable()
        delta = 1.0 / (ts.num_nodes + 1)
        x = 0
        for node in range(ts.num_nodes):
            site_id = sites.add_row(x, ancestral_state="0")
            x += delta
            mutations.add_row(site_id, node=node, derived_state="1")
        ts = msprime.load_tables(
            nodes=tables.nodes, edges=tables.edges,
            sites=sites, mutations=mutations)
        tree = ts.first()
        assert any(
            tree.parent(mut.node) == msprime.NULL_NODE
            for mut in tree.mutations())
        return tree

    def get_unary_node_tree(self):
        ts = msprime.simulate(2, random_seed=1)
        tables = ts.dump_tables()
        edges = tables.edges
        # Take out all the edges except 1
        n = 1
        edges.set_columns(
            left=edges.left[:n], right=edges.right[:n],
            parent=edges.parent[:n], child=edges.child[:n])
        ts = msprime.load_tables(nodes=tables.nodes, edges=edges)
        for t in ts.trees():
            for u in t.nodes():
                if len(t.children(u)) == 1:
                    return t
        assert False

    def get_empty_tree(self):
        nodes = msprime.NodeTable()
        edges = msprime.EdgeTable()
        ts = msprime.load_tables(nodes=nodes, edges=edges, sequence_length=1)
        return next(ts.trees())


class TestFormats(TestTreeDraw):
    """
    Tests that formats are recognised correctly.
    """
    def test_svg_variants(self):
        t = self.get_binary_tree()
        for svg in ["svg", "SVG", "sVg"]:
            output = t.draw(format=svg)
            root = xml.etree.ElementTree.fromstring(output)
            self.assertEqual(root.tag, "{http://www.w3.org/2000/svg}svg")

    def test_default(self):
        # Default is SVG
        t = self.get_binary_tree()
        output = t.draw(format=None)
        root = xml.etree.ElementTree.fromstring(output)
        self.assertEqual(root.tag, "{http://www.w3.org/2000/svg}svg")
        output = t.draw()
        root = xml.etree.ElementTree.fromstring(output)
        self.assertEqual(root.tag, "{http://www.w3.org/2000/svg}svg")

    def test_ascii_variants(self):
        t = self.get_binary_tree()
        for fmt in ["ascii", "ASCII", "AScii"]:
            output = t.draw(format=fmt)
            self.assertRaises(
                xml.etree.ElementTree.ParseError, xml.etree.ElementTree.fromstring,
                output)

    def test_unicode_variants(self):
        t = self.get_binary_tree()
        for fmt in ["unicode", "UNICODE", "uniCODE"]:
            if IS_PY2:
                self.assertRaises(ValueError, t.draw, format=fmt)
            else:
                output = t.draw(format=fmt)
                self.assertRaises(
                    xml.etree.ElementTree.ParseError, xml.etree.ElementTree.fromstring,
                    output)

    def test_bad_formats(self):
        t = self.get_binary_tree()
        for bad_format in ["", "ASC", "SV", "jpeg"]:
            self.assertRaises(ValueError, t.draw, format=bad_format)


# TODO we should gather some of these tests into a superclass as they are
# very similar for SVG and ASCII.

class TestDrawText(TestTreeDraw):
    """
    Tests the ASCII tree drawing method.
    """
    drawing_format = "ascii"
    example_label = "XXX"

    def verify_basic_text(self, text):
        self.assertTrue(isinstance(text, str))
        # TODO surely something else we can verify about this...

    def test_draw_defaults(self):
        t = self.get_binary_tree()
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_draw_nonbinary(self):
        t = self.get_nonbinary_tree()
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_draw_multiroot(self):
        t = self.get_multiroot_tree()
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_draw_mutations_over_roots(self):
        t = self.get_mutations_over_roots_tree()
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_draw_unary(self):
        t = self.get_unary_node_tree()
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_draw_empty_tree(self):
        t = self.get_empty_tree()
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_even_num_children_tree(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           1
        2   1           2
        3   1           1
        4   1           4
        5   1           5
        6   1           7
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       6       0
        0       1       6       1
        0       1       6       2
        0       1       6       3
        0       1       6       4
        0       1       6       5
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        t = next(ts.trees())
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_odd_num_children_tree(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           1
        2   1           2
        3   1           1
        4   1           4
        5   1           5
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       5       0
        0       1       5       1
        0       1       5       2
        0       1       5       3
        0       1       5       4
        """)
        ts = msprime.load_text(nodes, edges, strict=False)
        t = next(ts.trees())
        text = t.draw(format=self.drawing_format)
        self.verify_basic_text(text)

    def test_node_labels(self):
        t = self.get_binary_tree()
        labels = {u: self.example_label for u in t.nodes()}
        text = t.draw(format=self.drawing_format, node_labels=labels)
        self.verify_basic_text(text)
        j = 0
        for _ in t.nodes():
            j = text[j:].find(self.example_label)
            self.assertNotEqual(j, -1)

    def test_no_node_labels(self):
        t = self.get_binary_tree()
        labels = {}
        text = t.draw(format=self.drawing_format, node_labels=labels)
        self.verify_basic_text(text)
        for u in t.nodes():
            self.assertEqual(text.find(str(u)), -1)


@unittest.skipIf(IS_PY2, "Unicode tree drawing not supported on Python 2")
class TestDrawUnicode(TestDrawText):
    """
    Tests the Unicode tree drawing method
    """
    drawing_format = "unicode"
    example_label = "\u20ac" * 10  # euro symbol

    def verify_text_rendering(self, drawn, drawn_tree, debug=False):
        if debug:
            print("Drawn:")
            print(drawn)
            print("Expected:")
            print(drawn_tree)
        tree_lines = drawn_tree.splitlines()
        drawn_lines = drawn.splitlines()
        self.assertEqual(len(tree_lines), len(drawn_lines))
        for l1, l2 in zip(tree_lines, drawn_lines):
            # Trailing white space isn't significant.
            self.assertEqual(l1.rstrip(), l2.rstrip())

    def test_simple_tree(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   1           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       2       0
        0       1       2       1
        """)
        tree = (
            " 2 \n"
            "┏┻┓\n"
            "0 1")
        ts = msprime.load_text(nodes, edges, strict=False)
        t = next(ts.trees())
        drawn = t.draw(format="unicode")
        self.verify_text_rendering(drawn, tree)

    def test_trident_tree(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   1           0
        3   1           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       3       0
        0       1       3       1
        0       1       3       2
        """)
        tree = (
            "  3  \n"
            "┏━╋━┓\n"
            "0 1 2\n")
        ts = msprime.load_text(nodes, edges, strict=False)
        t = next(ts.trees())
        drawn = t.draw(format="unicode")
        self.verify_text_rendering(drawn, tree)

    def test_pitchfork_tree(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           0
        2   1           0
        3   1           0
        4   1           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       4       0
        0       1       4       1
        0       1       4       2
        0       1       4       3
        """)
        tree = (
            "   4   \n"
            "┏━┳┻┳━┓\n"
            "0 1 2 3\n")
        ts = msprime.load_text(nodes, edges, strict=False)
        t = next(ts.trees())
        # No labels
        tree = (
            "   ┃   \n"
            "┏━┳┻┳━┓\n"
            "┃ ┃ ┃ ┃\n")
        drawn = t.draw(format="unicode", node_labels={})
        self.verify_text_rendering(drawn, tree)
        # Some lables
        tree = (
            "   ┃   \n"
            "┏━┳┻┳━┓\n"
            "0 ┃ ┃ 3\n")
        drawn = t.draw(format="unicode", node_labels={0: "0", 3: "3"})
        self.verify_text_rendering(drawn, tree)

    def test_stick_tree(self):
        nodes = six.StringIO("""\
        id  is_sample   time
        0   1           0
        1   1           1
        2   1           2
        """)
        edges = six.StringIO("""\
        left    right   parent  child
        0       1       1       0
        0       1       2       1
        """)
        tree = (
            "2\n"
            "┃\n"
            "1\n"
            "┃\n"
            "0\n")
        ts = msprime.load_text(nodes, edges, strict=False)
        t = next(ts.trees())
        drawn = t.draw(format="unicode")
        self.verify_text_rendering(drawn, tree)


class TestDrawSvg(TestTreeDraw):
    """
    Tests the SVG tree drawing.
    """
    def verify_basic_svg(self, svg, width=200, height=200):
        root = xml.etree.ElementTree.fromstring(svg)
        self.assertEqual(root.tag, "{http://www.w3.org/2000/svg}svg")
        self.assertEqual(width, int(root.attrib["width"]))
        self.assertEqual(height, int(root.attrib["height"]))

    def test_draw_file(self):
        t = self.get_binary_tree()
        fd, filename = tempfile.mkstemp(prefix="msprime_viz_")
        try:
            os.close(fd)
            svg = t.draw(path=filename)
            self.assertGreater(os.path.getsize(filename), 0)
            with open(filename) as tmp:
                other_svg = tmp.read()
            self.assertEqual(svg, other_svg)
        finally:
            os.unlink(filename)

    def test_draw_defaults(self):
        t = self.get_binary_tree()
        svg = t.draw()
        self.verify_basic_svg(svg)

    def test_draw_nonbinary(self):
        t = self.get_nonbinary_tree()
        svg = t.draw()
        self.verify_basic_svg(svg)

    def test_draw_multiroot(self):
        t = self.get_multiroot_tree()
        svg = t.draw()
        self.verify_basic_svg(svg)

    def test_draw_mutations_over_roots(self):
        t = self.get_mutations_over_roots_tree()
        svg = t.draw()
        self.verify_basic_svg(svg)

    def test_draw_unary(self):
        t = self.get_unary_node_tree()
        svg = t.draw()
        self.verify_basic_svg(svg)

    def test_draw_empty(self):
        t = self.get_empty_tree()
        svg = t.draw()
        self.verify_basic_svg(svg)

    def test_width_height(self):
        t = self.get_binary_tree()
        w = 123
        h = 456
        svg = t.draw(width=w, height=h)
        self.verify_basic_svg(svg, w, h)

    def test_node_labels(self):
        t = self.get_binary_tree()
        labels = {u: "XXX" for u in t.nodes()}
        svg = t.draw(format="svg", node_labels=labels)
        self.verify_basic_svg(svg)
        self.assertEqual(svg.count("XXX"), t.num_nodes)

    def test_one_node_label(self):
        t = self.get_binary_tree()
        labels = {0: "XXX"}
        svg = t.draw(format="svg", node_labels=labels)
        self.verify_basic_svg(svg)
        self.assertEqual(svg.count("XXX"), 1)

    def test_no_node_labels(self):
        t = self.get_binary_tree()
        labels = {}
        svg = t.draw(format="svg", node_labels=labels)
        self.verify_basic_svg(svg)
        # Can't really test for much here if we don't understand the SVG

    def test_one_node_colour(self):
        t = self.get_binary_tree()
        colour = "rgb(0, 1, 2)"
        colours = {0: colour}
        svg = t.draw(format="svg", node_colours=colours)
        self.verify_basic_svg(svg)
        self.assertEqual(svg.count('fill="{}"'.format(colour)), 1)

    def test_all_nodes_colour(self):
        t = self.get_binary_tree()
        colours = {u: "rgb({}, {}, {})".format(u, u, u) for u in t.nodes()}
        svg = t.draw(format="svg", node_colours=colours)
        self.verify_basic_svg(svg)
        for colour in colours.values():
            self.assertEqual(svg.count('fill="{}"'.format(colour)), 1)

    def test_mutation_labels(self):
        t = self.get_binary_tree()
        labels = {u.id: "XXX" for u in t.mutations()}
        svg = t.draw(format="svg", mutation_labels=labels)
        self.verify_basic_svg(svg)
        self.assertEqual(svg.count("XXX"), t.num_mutations)

    def test_one_mutation_label(self):
        t = self.get_binary_tree()
        labels = {0: "XXX"}
        svg = t.draw(format="svg", mutation_labels=labels)
        self.verify_basic_svg(svg)
        self.assertEqual(svg.count("XXX"), 1)

    def test_no_mutation_labels(self):
        t = self.get_binary_tree()
        labels = {}
        svg = t.draw(format="svg", mutation_labels=labels)
        self.verify_basic_svg(svg)
        # Can't really test for much here if we don't understand the SVG

    def test_one_mutation_colour(self):
        t = self.get_binary_tree()
        colour = "rgb(0, 1, 2)"
        colours = {0: colour}
        svg = t.draw(format="svg", mutation_colours=colours)
        self.verify_basic_svg(svg)
        self.assertEqual(svg.count('fill="{}"'.format(colour)), 1)

    def test_all_mutations_colour(self):
        t = self.get_binary_tree()
        colours = {
            mut.id: "rgb({}, {}, {})".format(mut.id, mut.id, mut.id)
            for mut in t.mutations()}
        svg = t.draw(format="svg", mutation_colours=colours)
        self.verify_basic_svg(svg)
        for colour in colours.values():
            self.assertEqual(svg.count('fill="{}"'.format(colour)), 1)
