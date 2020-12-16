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

# TODO document these flags
from msprime import _msprime
from msprime._msprime import NODE_IS_CA_EVENT
from msprime._msprime import NODE_IS_CEN_EVENT
from msprime._msprime import NODE_IS_MIG_EVENT
from msprime._msprime import NODE_IS_RE_EVENT

from msprime.core import __version__
from msprime.exceptions import *
from msprime.demography import *
from msprime.ancestry import *
from msprime.pedigrees import *
from msprime.mutations import *
from msprime.likelihood import *
from msprime.intervals import *

from types import ModuleType

old_module = sys.modules["msprime"]

# Many attributes were moved to tskit, we use a facade here so we don't break old
# code, but do emit a warning
class DeprecationFacade(ModuleType):
    def __getattr__(self, name):
        try:
            # Tree used to be SparseTree
            if name == "SparseTree":
                name = "Tree"
            result = getattr(tskit, name)
            warnings.warn(
                f"'{__name__}.{name}' is deprecated and will be removed"
                f" in future versions. Use 'tskit.{name}'."
            )
        except AttributeError:
            raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
        return result


# Patch the facade into the dict of loaded modules
facade = sys.modules["msprime"] = DeprecationFacade("msprime")
facade.__dict__.update(old_module.__dict__)
