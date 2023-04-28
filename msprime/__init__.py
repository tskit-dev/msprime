# Turn off flake8 and reorder-python-imports for this file.
# flake8: NOQA
# noreorder
#
# Copyright (C) 2015-2021 University of Oxford
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
Msprime is a population genetics simulator.
Please see the documentation at https://tskit.dev/msprime/docs/
for more information.
"""

from msprime._msprime import (
    NODE_IS_CEN_EVENT,
    NODE_IS_CA_EVENT,
    NODE_IS_MIG_EVENT,
    NODE_IS_RE_EVENT,
    NODE_IS_GC_EVENT,
    NODE_IS_PASS_THROUGH,
)

from msprime.ancestry import (
    AncestryModel,
    BetaCoalescent,
    DiracCoalescent,
    DiscreteTimeWrightFisher,
    SampleSet,
    sim_ancestry,
    SmcApproxCoalescent,
    SmcPrimeApproxCoalescent,
    StandardCoalescent,
    SweepGenicSelection,
    FixedPedigree,
    TimeUnitsMismatchWarning,
    NodeType,
)

from msprime.core import __version__

from msprime.demography import (
    Demography,
    DemographyDebugger,
    IncompletePopulationMetadataWarning,
    Population,
)

from msprime.pedigrees import PedigreeBuilder, parse_pedigree, write_pedigree
from msprime.intervals import RateMap
from msprime.likelihood import log_arg_likelihood, log_mutation_likelihood

from msprime.mutations import (
    BinaryMutationModel,
    BLOSUM62,
    F84,
    GTR,
    HKY,
    InfiniteAlleles,
    JC69,
    MatrixMutationModel,
    MicrosatMutationModel,
    MutationModel,
    SMM,
    TPM,
    EL2,
    NUCLEOTIDES,
    PAM,
    SLiMMutationModel,
    sim_mutations,
)

# Imports for deprecated 0.x classes and functions. We keep these separate
# here for clarity, but they will be maintained indefinitely.

from msprime.ancestry import (
    Sample,
    simulate,
    SimulationModelChange,
)
from msprime.demography import (
    CensusEvent,
    InstantaneousBottleneck,
    MassMigration,
    MigrationRateChange,
    PopulationConfiguration,
    PopulationParametersChange,
    SimpleBottleneck,
)
from msprime.mutations import (
    BINARY,
    mutate,
    NUCLEOTIDES,
    InfiniteSites,
)
from msprime.intervals import RecombinationMap

__all__ = [
    "AncestryModel",
    "BINARY",
    "BLOSUM62",
    "BetaCoalescent",
    "BinaryMutationModel",
    "CensusEvent",
    "Demography",
    "DemographyDebugger",
    "DiracCoalescent",
    "DiscreteTimeWrightFisher",
    "F84",
    "GTR",
    "HKY",
    "IncompletePopulationMetadataWarning",
    "InfiniteAlleles",
    "InfiniteSites",
    "InstantaneousBottleneck",
    "JC69",
    "MassMigration",
    "MatrixMutationModel",
    "MicrosatMutationModel",
    "SMM",
    "TPM",
    "EL2",
    "MigrationRateChange",
    "MutationModel",
    "NODE_IS_CA_EVENT",
    "NODE_IS_CEN_EVENT",
    "NODE_IS_MIG_EVENT",
    "NODE_IS_RE_EVENT",
    "NODE_IS_GC_EVENT",
    "NodeType",
    "NUCLEOTIDES",
    "PAM",
    "PedigreeBuilder",
    "Population",
    "PopulationConfiguration",
    "PopulationParametersChange",
    "RateMap",
    "RecombinationMap",
    "SLiMMutationModel",
    "Sample",
    "SampleSet",
    "SimpleBottleneck",
    "SimulationModelChange",
    "SmcApproxCoalescent",
    "SmcPrimeApproxCoalescent",
    "StandardCoalescent",
    "SweepGenicSelection",
    "FixedPedigree",
    "log_arg_likelihood",
    "log_mutation_likelihood",
    "mutate",
    "parse_pedigree",
    "sim_ancestry",
    "sim_mutations",
    "simulate",
    "write_pedigree",
]


def __make_tskit_deprecations():
    """
    Update the namespace to include redirects from the old tskit classes
    that used to be in msprime. These will raise a formal warning for
    now and will be removed in 1.2.

    We do this in a function to avoid cluttering up the top-level namespace
    with more stuff.
    """

    import sys
    import tskit
    import warnings
    import types

    old_module = sys.modules["msprime"]

    # Many attributes were moved to tskit, we use a facade here so we don't break old
    # code, but do emit a warning
    msprime_names_now_in_tskit = [
        "Edge",
        "EdgeTable",
        "Edgeset",
        "FORWARD",
        "Individual",
        "IndividualTable",
        "LdCalculator",
        "Migration",
        "MigrationTable",
        "Mutation",
        "MutationTable",
        "NODE_IS_SAMPLE",
        "Node",
        "NodeTable",
        "PopulationTable",
        "Provenance",
        "ProvenanceTable",
        "REVERSE",
        "Site",
        "SiteTable",
        "TableCollection",
        "Tree",
        "TreeSequence",
        "Variant",
        "load",
        "load_text",
        "parse_nodes",
        "pack_bytes",
        "pack_strings",
        "parse_edges",
        "parse_individuals",
        "parse_mutations",
        "parse_sites",
        "unpack_bytes",
        "unpack_strings",
        "validate_provenance",
        "NULL",
    ]

    old_msprime_stuff_renamed = {
        "NULL_NODE": "NULL",
        "NULL_POPULATION": "NULL",
        "NULL_INDIVIDUAL": "NULL",
        "NULL_MUTATION": "NULL",
        "SparseTree": "Tree",
    }

    class DeprecationFacade(types.ModuleType):
        def __getattr__(self, name):
            if name in old_msprime_stuff_renamed:
                new_name = old_msprime_stuff_renamed[name]
                warnings.warn(
                    f"'{__name__}.{name}' is deprecated and will be removed"
                    f" in future versions. Use 'tskit.{new_name}'.",
                    category=FutureWarning,
                )
                return getattr(tskit, new_name)
            if name in msprime_names_now_in_tskit:
                result = getattr(tskit, name)
                warnings.warn(
                    f"'{__name__}.{name}' is deprecated and will be removed"
                    f" in future versions. Use 'tskit.{name}'.",
                    category=FutureWarning,
                )
                return result
            raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    # Patch the facade into the dict of loaded modules
    facade = DeprecationFacade("msprime")
    facade.__dict__.update(old_module.__dict__)
    sys.modules["msprime"] = facade


__make_tskit_deprecations()
