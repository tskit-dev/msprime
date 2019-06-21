#
# Copyright (C) 2015-2018 University of Oxford
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
Msprime is a reimplementation of Hudson's classical ms simulator for
modern datasets.
"""
# flake8: NOQA
import tskit

# Compatibility layer for old code using the tskit API exported by msprime.
NULL_NODE = tskit.NULL
NULL_POPULATION = tskit.NULL
NULL_INDIVIDUAL = tskit.NULL
NULL_MUTATION = tskit.NULL

from tskit import (
    Individual, Node, Edge, Site, Mutation, Migration, Population,
    Variant, Edgeset, Provenance)
from tskit import Tree as SparseTree  # Rename SparseTree to Tree in tskit
from tskit import TreeSequence
from tskit import (
    IndividualTable, NodeTable, EdgeTable, SiteTable, MutationTable,
    MigrationTable, PopulationTable, ProvenanceTable, TableCollection)
from tskit import LdCalculator
from tskit import load, load_text
from tskit import (
    parse_nodes, parse_edges, parse_individuals, parse_sites, parse_mutations)
from tskit import pack_strings, pack_bytes, unpack_bytes, unpack_strings
from tskit import validate_provenance

from tskit import NODE_IS_SAMPLE, FORWARD, REVERSE

from msprime.provenance import __version__
from msprime.simulations import *
from msprime.exceptions import *
from msprime.mutations import *
from msprime.likelihood import *
