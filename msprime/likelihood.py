#
# Copyright (C) 2019 University of Oxford
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
Module responsible for computing likelihoods.
"""
import math

import _msprime


def unnormalised_log_mutation_likelihood(arg, theta):
    # log_likelihood of mutations on a given ARG up to a normalising constant
    # that depends on the pattern of observed mutations, but not on the ARG
    # or the mutation rate
    tables = arg.tables
    time = tables.nodes.time
    total_material = 0
    for e in tables.edges:
        total_material += (e.right - e.left) * (time[e.parent] - time[e.child])
    number_of_mutations = len(tables.mutations)
    if theta == 0:
        if number_of_mutations == 0:
            ret = 0
        else:
            ret = -float("inf")
    else:
        ret = (number_of_mutations * math.log(total_material * theta) -
               total_material * theta)
    for tree in arg.trees():
        for site in tree.sites():
            mutation = site.mutations[0]
            child = mutation.node
            parent = tree.parent(child)
            potential_branch_length = tree.branch_length(child)
            while tree.parent(parent) is not None and len(tree.children(parent)) == 1:
                child = parent
                parent = tree.parent(child)
                potential_branch_length += tree.branch_length(child)
            child = mutation.node
            while len(tree.children(child)) == 1:
                child = tree.children(child)[0]
                potential_branch_length += tree.branch_length(child)
            ret += math.log(potential_branch_length / total_material)
    return ret


def log_arg_likelihood(arg, recombination_rate, Ne=0.25):
    # TODO: Ne should default to 1 for compatability with msprime.simulate. Setting
    # to 1/4 now to keep the tests working.

    # Get the tables into the format we need to interchange with the low-level code.
    lw_tables = _msprime.LightweightTableCollection()
    lw_tables.fromdict(arg.tables.asdict())
    return _msprime.log_likelihood_arg(
        lw_tables, Ne=Ne, recombination_rate=recombination_rate)
