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
    edges_above = collections.defaultdict(list)
    edges_below = collections.defaultdict(list)
    total_material = 0
    time = tables.nodes.time
    for edge in tables.edges:
        total_material += ((edge.right - edge.left) *
                           (time[edge.parent] - time[edge.child]))
        edges_above[edge.child].append(edge)
        edges_below[edge.parent].append(edge)
    number_of_mutations = len(tables.mutations)
    if theta == 0:
        if number_of_mutations == 0:
            ret = 0
        else:
            ret = -float("inf")
    else:
        ret = (number_of_mutations * math.log(total_material * theta) -
               total_material * theta)
    for mut in arg.mutations():
        mutant_location = tables.sites[mut.site].position
        index = 0
        to_check = edges_above[mut.node]
        while to_check[index].right < mutant_location:
            index += 1
        edge = edges_above[mut.node][index]
        potential_branch_length = time[edge.parent] - time[edge.child]
        continue_leafwards = time[edge.child] > 0
        tmp_edge = edge
        while continue_leafwards:
            count = 0
            for e in edges_below[tmp_edge.child]:
                if e.left <= mutant_location and mutant_location < e.right:
                    count += 1
                    tmp_edge = e
            if count == 1:
                potential_branch_length += time[tmp_edge.parent] - time[tmp_edge.child]
                continue_leafwards = time[tmp_edge.child] > 0
            else:
                continue_leafwards = False
        tmp_edge = edge
        continue_rootwards = tmp_edge.parent in edges_above
        while continue_rootwards:
            count = 0
            for e in edges_below[tmp_edge.parent]:
                if e.left <= mutant_location and mutant_location < e.right:
                    count += 1
            if count == 1:
                to_check = edges_above[tmp_edge.parent]
                index = 0
                while to_check[index].right < mutant_location:
                    index += 1
                tmp_edge = to_check[index]
                potential_branch_length += time[tmp_edge.parent] - time[tmp_edge.child]
                continue_rootwards = tmp_edge.parent in edges_above
            else:
                continue_rootwards = False
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
