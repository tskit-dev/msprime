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
import collections
import math
import msprime


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


def log_arg_likelihood(arg, rho):
    tables = arg.tables
    number_of_lineages = arg.num_samples
    number_of_links = number_of_lineages * arg.sequence_length
    number_of_edges = len(tables.edges)
    edges_above = collections.defaultdict(list)
    for edge in tables.edges:
        edges_above[edge.child].append(edge)
    edge = 0
    time = 0
    ret = 0
    while edge < number_of_edges:
        parent = tables.edges[edge].parent
        rate = number_of_lineages * (number_of_lineages - 1) / 2 + number_of_links * rho
        ret -= rate * (tables.nodes[parent].time - time)
        time = tables.nodes[parent].time
        child = tables.edges[edge].child
        if tables.nodes[parent].flags == msprime.NODE_IS_RE_EVENT:
            while tables.edges[edge].child == child:
                if tables.edges[edge].parent != tables.edges[edge - 1].parent:
                    gap = tables.edges[edge].left - tables.edges[edge - 1].right
                edge += 1
            number_of_links -= gap
            number_of_lineages += 1
            if gap == 0:
                # we evaluate the density rather than the probability if there is no gap
                gap = 1
            if rho == 0:
                ret -= float("inf")
            else:
                ret += math.log(rho * gap)
        else:
            segment_length_in_children = -tables.edges.left[edge]
            while edge < number_of_edges and tables.edges[edge].child == child:
                edge += 1
            segment_length_in_children += (tables.edges.right[edge - 1] -
                                           tables.edges.left[edge])
            child = tables.edges[edge].child
            while edge < number_of_edges and tables.edges[edge].child == child:
                edge += 1
            segment_length_in_children += tables.edges.right[edge - 1]
            if parent in edges_above:
                segment_length_in_parent = (edges_above[parent][-1].right -
                                            edges_above[parent][0].left)
                number_of_lineages -= 1
                number_of_links -= segment_length_in_children - segment_length_in_parent
            else:
                number_of_lineages -= 2
                number_of_links -= segment_length_in_children
    return ret
