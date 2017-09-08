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
import tempfile
import unittest
import xml.etree

import msprime


class TestTreeDraw(unittest.TestCase):
    """
    Tests for the tree drawing functionality.
    """
    def get_binary_tree(self):
        ts = msprime.simulate(10, random_seed=1, mutation_rate=1)
        return next(ts.trees())

    def get_nonbinary_tree(self):
        demographic_events = [
            msprime.SimpleBottleneck(time=0.1, proportion=0.5)]
        ts = msprime.simulate(
            10, recombination_rate=5, mutation_rate=10,
            demographic_events=demographic_events, random_seed=1)
        for t in ts.trees():
            for u in t.nodes():
                if len(t.children(u)) > 2:
                    return t
        assert False


class TestDrawAscii(TestTreeDraw):
    """
    Tests the ASCII tree drawing method.
    """
    def verify_basic_text(self, text):
        self.assertTrue(isinstance(text, str))
        # TODO surely something else we can verify about this...

    def test_draw_defaults(self):
        t = self.get_binary_tree()
        text = t.draw_ascii()
        self.verify_basic_text(text)

    def test_draw_nonbinary(self):
        t = self.get_nonbinary_tree()
        text = t.draw_ascii()
        self.verify_basic_text(text)

    def test_labels(self):
        t = self.get_binary_tree()
        labels = {u: "XXX" for u in t.nodes()}
        text = t.draw_ascii(labels=labels)
        self.verify_basic_text(text)
        j = 0
        for _ in t.nodes():
            j = text[j:].find("XXX")
            self.assertNotEqual(j, -1)


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

    def test_width_height(self):
        t = self.get_binary_tree()
        w = 123
        h = 456
        svg = t.draw(width=w, height=h)
        self.verify_basic_svg(svg, w, h)

    def test_boolean_flags(self):
        flags = [
            "times", "mutation_locations", "mutation_labels", "internal_node_labels",
            "leaf_node_labels", "show_times"]
        t = self.get_binary_tree()
        for flag in flags:
            for boolean in [True, False]:
                svg = t.draw(**{flag: boolean})
                self.verify_basic_svg(svg)
