import math
import msprime
import tskit
import numpy as np

def unnormalised_log_mutation_likelihood(arg, theta):
    # log_likelihood of mutations on a given ARG up to a normalising constant
    # that depends on the pattern of observed mutations, but not on the ARG
    # or the mutation rate
    tables = arg.tables
    number_of_edges = len(tables.edges)
    edges_above = {}
    edge = 0
    total_material = 0
    while edge < number_of_edges:
        child = tables.edges[edge].child
        right = tables.edges[edge].right
        left = tables.edges[edge].left
        parent = tables.edges[edge].parent
        total_material += (right - left) * (tables.nodes[parent].time - tables.nodes[child].time)
        edges_above[child] = [edge]
        edge += 1
        if edge == number_of_edges:
            break
        while tables.edges[edge].child == child:
            right = tables.edges[edge].right
            left = tables.edges[edge].left
            parent = tables.edges[edge].parent
            total_material += (right - left) * (tables.nodes[parent].time - tables.nodes[child].time)
            edges_above[child].append(edge)
            edge += 1
            if edge == number_of_edges:
                break
    edges_below = {}
    edge = 0
    while edge < number_of_edges:
        parent = tables.edges[edge].parent
        edges_below[parent] = [edge]
        edge += 1
        if edge == number_of_edges:
            break
        while tables.edges[edge].parent == parent:
            edges_below[parent].append(edge)
            edge += 1
            if edge == number_of_edges:
                break
    number_of_mutations = len(tables.mutations)
    if number_of_mutations == 0 and theta == 0:
        ret = -float("inf")
    else:
        ret = number_of_mutations * math.log(total_material * theta) - total_material * theta
    for mut in arg.mutations():
        mutant_location = tables.sites[mut.site].position
        ind = 0
        edge = edges_above[mut.node][0]
        while tables.edges[edge].right < mutant_location:
            ind += 1
            edge = edges_above[mut.node][ind]
        parent = tables.edges[edge].parent
        child = tables.edges[edge].child
        potential_branch_length = tables.nodes[parent].time - tables.nodes[child].time
        continue_downwards = True
        if tables.nodes[child].time == 0:
            continue_downwards = False
        tmp_edge = edge
        while continue_downwards:
            count = 0
            to_check = edges_below[child]
            for e in to_check:
                if tables.edges[e].left <= mutant_location and mutant_location < tables.edges[e].right:
                    count += 1
                    tmp_edge = e
            if count == 1:
                parent = tables.edges[tmp_edge].parent
                child = tables.edges[tmp_edge].child
                potential_branch_length += tables.nodes[parent].time - tables.nodes[child].time
                if tables.nodes[child].time == 0:
                    continue_downwards = False
            else:
                continue_downwards = False
        continue_upwards = True
        tmp_edge = edge
        parent = tables.edges[edge].parent
        child = tables.edges[edge].child
        if parent not in edges_above:
            continue_upwards = False
        while continue_upwards:
            count_up = 0
            to_check = edges_above[parent]
            for e in to_check:
                if tables.edges[e].left <= mutant_location and mutant_location < tables.edges[e].right:
                    count_up += 1
                    tmp_edge = e
            count_down = 0
            to_check = edges_below[parent]
            for e in to_check:
                if tables.edges[e].left <= mutant_location and mutant_location < tables.edges[e].right:
                    count_down += 1
            if count_up == 1 and count_down == 1:
                parent = tables.edges[tmp_edge].parent
                child = tables.edges[tmp_edge].child
                potential_branch_length += tables.nodes[parent].time - tables.nodes[child].time
                if parent not in edges_above:
                    continue_upwards = False
            else:
                continue_upwards = False
        ret += math.log(potential_branch_length / total_material)
    return ret

def log_arg_likelihood(arg, rho):
    tables = arg.tables
    number_of_lineages = arg.num_samples
    number_of_links = number_of_lineages * arg.sequence_length
    number_of_edges = len(tables.edges)
    edge = 0
    parent = tables.edges[edge].parent
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
            while tables.edges[edge].child == child:
                edge += 1
                if edge == number_of_edges:
                    break
            segment_length_in_children += tables.edges.right[edge - 1] - tables.edges.left[edge]
            child = tables.edges[edge].child
            while tables.edges[edge].child == child:
                edge += 1
                if edge == number_of_edges:
                    break
            segment_length_in_children += tables.edges.right[edge - 1]
            parent_edges = edge + np.flatnonzero(tables.edges.child[edge:] == parent)
            if len(parent_edges) > 0:
                segment_length_in_parent = tables.edges.right[parent_edges[-1]] - tables.edges.left[parent_edges[0]]
                number_of_lineages -= 1
                number_of_links -= segment_length_in_children - segment_length_in_parent
            else:
                number_of_lineages -= 2
                number_of_links -= segment_length_in_children
    return ret
