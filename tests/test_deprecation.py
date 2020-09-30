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
Test cases for deprecating legacy functionality
"""
import pytest

import msprime


class TestDeprecation:
    def test_moved_to_tskit(self):
        # Note that "Population" is removed here as it now refers
        # to msprime.demography.Population
        for name in [
            "Individual",
            "Node",
            "Edge",
            "Site",
            "Mutation",
            "Migration",
            "Variant",
            "Edgeset",
            "Provenance",
            "Tree",
            "SparseTree",
            "TreeSequence",
            "IndividualTable",
            "NodeTable",
            "EdgeTable",
            "SiteTable",
            "MutationTable",
            "MigrationTable",
            "PopulationTable",
            "ProvenanceTable",
            "TableCollection",
            "LdCalculator",
            "load",
            "load_text",
            "parse_nodes",
            "parse_edges",
            "parse_individuals",
            "parse_sites",
            "parse_mutations",
            "pack_strings",
            "pack_bytes",
            "unpack_bytes",
            "unpack_strings",
            "validate_provenance",
            "NODE_IS_SAMPLE",
            "FORWARD",
            "REVERSE",
        ]:
            with pytest.warns(UserWarning):
                getattr(msprime, name)
