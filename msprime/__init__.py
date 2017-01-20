#
# Copyright (C) 2015 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
from __future__ import print_function
from __future__ import division

from _msprime import FORWARD  # NOQA
from _msprime import REVERSE  # NOQA

from msprime.environment import __version__  # NOQA
from msprime.formats import *  # NOQA
from msprime.trees import *  # NOQA
from msprime.stats import *  # NOQA
