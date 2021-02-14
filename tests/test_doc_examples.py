#
# Copyright (C) 2020 University of Oxford
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
Test that the documentation examples do what they are supposed to.
"""
import contextlib
import io
import sys
from unittest import mock

from docs import examples


@contextlib.contextmanager
def capture_stdout():
    new_out = io.StringIO()
    old_out = sys.stdout
    try:
        sys.stdout = new_out
        yield sys.stdout
    finally:
        sys.stdout = old_out


class TestDocumentationExamples:
    def test_segregating_sites(self):
        with capture_stdout() as stdout:
            examples.segregating_sites(10, 5, 10)
        output = stdout.getvalue().splitlines()
        assert len(output) == 3
        assert output[1].startswith("Observed")
        assert output[2].startswith("Analytical")

    def test_migration(self):
        with capture_stdout() as stdout:
            examples.migration_example(10)
        output = stdout.getvalue().splitlines()
        assert len(output) == 2
        assert output[0].startswith("Observed")
        assert output[1].startswith("Predicted")

    def test_logging_info(self):
        with mock.patch("daiquiri.setup"):
            examples.logging_info_example()

    def test_logging_debug(self):
        with mock.patch("daiquiri.setup"):
            examples.logging_debug_example()
