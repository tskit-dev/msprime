#
# Copyright (C) 2017 University of Oxford
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
A collection of utilities to edit and construct tree sequences.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import json
import random

import numpy as np

import msprime.provenance as provenance
import msprime


def add_provenance(provenance_table, method_name):
    d = provenance.get_provenance_dict("tsutil.{}".format(method_name))
    provenance_table.add_row(json.dumps(d))


def subsample_sites(ts, num_sites):
    """
    Returns a copy of the specified tree sequence with a random subsample of the
    specified number of sites.
    """
    t = ts.dump_tables()
    t.sites.reset()
    t.mutations.reset()
    sites_to_keep = set(random.sample(list(range(ts.num_sites)), num_sites))
    for site in ts.sites():
        if site.id in sites_to_keep:
            site_id = len(t.sites)
            t.sites.add_row(
                position=site.position, ancestral_state=site.ancestral_state)
            for mutation in site.mutations:
                t.mutations.add_row(
                    site=site_id, derived_state=mutation.derived_state,
                    node=mutation.node, parent=mutation.parent)
    add_provenance(t.provenances, "subsample_sites")
    return msprime.load_tables(**t.asdict())


def decapitate(ts, num_edges):
    """
    Returns a copy of the specified tree sequence in which the specified number of
    edges have been retained.
    """
    t = ts.dump_tables()
    t.edges.set_columns(
        left=t.edges.left[:num_edges], right=t.edges.right[:num_edges],
        parent=t.edges.parent[:num_edges], child=t.edges.child[:num_edges])
    add_provenance(t.provenances, "decapitate")
    return t.tree_sequence()


def insert_branch_mutations(ts, mutations_per_branch=1):
    """
    Returns a copy of the specified tree sequence with a mutation on every branch
    in every tree.
    """
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        site = tables.sites.add_row(position=tree.interval[0], ancestral_state='0')
        for root in tree.roots:
            state = {root: 0}
            mutation = {root: -1}
            stack = [root]
            while len(stack) > 0:
                u = stack.pop()
                stack.extend(tree.children(u))
                v = tree.parent(u)
                if v != msprime.NULL_NODE:
                    state[u] = state[v]
                    parent = mutation[v]
                    for j in range(mutations_per_branch):
                        state[u] = (state[u] + 1) % 2
                        mutation[u] = tables.mutations.add_row(
                            site=site, node=u, derived_state=str(state[u]),
                            parent=parent)
                        parent = mutation[u]
    add_provenance(tables.provenances, "insert_branch_mutations")
    return tables.tree_sequence()


def insert_branch_sites(ts):
    """
    Returns a copy of the specified tree sequence with a site on every branch
    of every tree.
    """
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        left, right = tree.interval
        delta = (right - left) / len(list(tree.nodes()))
        x = left
        for u in tree.nodes():
            if tree.parent(u) != msprime.NULL_NODE:
                site = tables.sites.add_row(position=x, ancestral_state='0')
                tables.mutations.add_row(site=site, node=u, derived_state='1')
                x += delta
    add_provenance(tables.provenances, "insert_branch_sites")
    return tables.tree_sequence()


def insert_multichar_mutations(ts, seed=1, max_len=10):
    """
    Returns a copy of the specified tree sequence with multiple chararacter
    mutations on a randomly chosen branch in every tree.
    """
    rng = random.Random(seed)
    letters = ["A", "C", "T", "G"]
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        ancestral_state = rng.choice(letters) * rng.randint(0, max_len)
        site = tables.sites.add_row(
            position=tree.interval[0], ancestral_state=ancestral_state)
        nodes = list(tree.nodes())
        nodes.remove(tree.root)
        u = rng.choice(nodes)
        derived_state = ancestral_state
        while ancestral_state == derived_state:
            derived_state = rng.choice(letters) * rng.randint(0, max_len)
        tables.mutations.add_row(site=site, node=u, derived_state=derived_state)
    add_provenance(tables.provenances, "insert_multichar_mutations")
    return tables.tree_sequence()


def insert_random_ploidy_individuals(ts, max_ploidy=5, max_dimension=3, seed=1):
    """
    Takes random contiguous subsets of the samples an assigns them to individuals.
    Also creates random locations in variable dimensions in the unit interval.
    """
    rng = random.Random(seed)
    samples = np.array(ts.samples(), dtype=int)
    j = 0
    tables = ts.dump_tables()
    tables.individuals.clear()
    individual = tables.nodes.individual[:]
    individual[:] = msprime.NULL_INDIVIDUAL
    while j < len(samples):
        ploidy = rng.randint(0, max_ploidy)
        nodes = samples[j: min(j + ploidy, len(samples))]
        dimension = rng.randint(0, max_dimension)
        location = [rng.random() for _ in range(dimension)]
        ind_id = tables.individuals.add_row(location=location)
        individual[nodes] = ind_id
        j += ploidy
    tables.nodes.individual = individual
    return msprime.load_tables(**tables.asdict())


def permute_nodes(ts, node_map):
    """
    Returns a copy of the specified tree sequence such that the nodes are
    permuted according to the specified map.
    """
    tables = ts.dump_tables()
    tables.nodes.clear()
    tables.edges.clear()
    tables.mutations.clear()
    # Mapping from nodes in the new tree sequence back to nodes in the original
    reverse_map = [0 for _ in node_map]
    for j in range(ts.num_nodes):
        reverse_map[node_map[j]] = j
    old_nodes = list(ts.nodes())
    for j in range(ts.num_nodes):
        old_node = old_nodes[reverse_map[j]]
        tables.nodes.add_row(
            flags=old_node.flags, metadata=old_node.metadata,
            population=old_node.population, time=old_node.time)
    for edge in ts.edges():
        tables.edges.add_row(
            left=edge.left, right=edge.right, parent=node_map[edge.parent],
            child=node_map[edge.child])
    for site in ts.sites():
        for mutation in site.mutations:
            tables.mutations.add_row(
                site=site.id, derived_state=mutation.derived_state,
                node=node_map[mutation.node], metadata=mutation.metadata)
    tables.sort()
    add_provenance(tables.provenances, "permute_nodes")
    return tables.tree_sequence()


def insert_redundant_breakpoints(ts):
    """
    Builds a new tree sequence containing redundant breakpoints.
    """
    tables = ts.dump_tables()
    tables.edges.reset()
    for r in ts.edges():
        x = r.left + (r.right - r.left) / 2
        tables.edges.add_row(left=r.left, right=x, child=r.child, parent=r.parent)
        tables.edges.add_row(left=x, right=r.right, child=r.child, parent=r.parent)
    add_provenance(tables.provenances, "insert_redundant_breakpoints")
    new_ts = msprime.load_tables(**tables.asdict())
    assert new_ts.num_edges == 2 * ts.num_edges
    return new_ts


def single_childify(ts):
    """
    Builds a new equivalent tree sequence which contains an extra node in the
    middle of all exising branches.
    """
    tables = ts.dump_tables()

    time = tables.nodes.time[:]
    tables.edges.reset()
    for edge in ts.edges():
        # Insert a new node in between the parent and child.
        t = time[edge.child] + (time[edge.parent] - time[edge.child]) / 2
        u = tables.nodes.add_row(time=t)
        tables.edges.add_row(
            left=edge.left, right=edge.right, parent=u, child=edge.child)
        tables.edges.add_row(
            left=edge.left, right=edge.right, parent=edge.parent, child=u)
    tables.sort()
    add_provenance(tables.provenances, "insert_redundant_breakpoints")
    return tables.tree_sequence()


def add_random_metadata(ts, seed=1, max_length=10):
    """
    Returns a copy of the specified tree sequence with random metadata assigned
    to the nodes, sites and mutations.
    """
    tables = ts.dump_tables()
    np.random.seed(seed)

    length = np.random.randint(0, max_length, ts.num_nodes)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    # Older versions of numpy didn't have a dtype argument for randint, so
    # must use astype instead.
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    nodes = tables.nodes
    nodes.set_columns(
        flags=nodes.flags, population=nodes.population, time=nodes.time,
        metadata_offset=offset, metadata=metadata,
        individual=nodes.individual)

    length = np.random.randint(0, max_length, ts.num_sites)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    sites = tables.sites
    sites.set_columns(
        position=sites.position,
        ancestral_state=sites.ancestral_state,
        ancestral_state_offset=sites.ancestral_state_offset,
        metadata_offset=offset, metadata=metadata)

    length = np.random.randint(0, max_length, ts.num_mutations)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    mutations = tables.mutations
    mutations.set_columns(
        site=mutations.site,
        node=mutations.node,
        parent=mutations.parent,
        derived_state=mutations.derived_state,
        derived_state_offset=mutations.derived_state_offset,
        metadata_offset=offset, metadata=metadata)

    length = np.random.randint(0, max_length, ts.num_individuals)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    individuals = tables.individuals
    individuals.set_columns(
        flags=individuals.flags,
        location=individuals.location,
        location_offset=individuals.location_offset,
        metadata_offset=offset, metadata=metadata)

    length = np.random.randint(0, max_length, ts.num_populations)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    populations = tables.populations
    populations.set_columns(metadata_offset=offset, metadata=metadata)

    add_provenance(tables.provenances, "add_random_metadata")
    ts = msprime.load_tables(**tables.asdict())
    return ts


def jiggle_samples(ts):
    """
    Returns a copy of the specified tree sequence with the sample nodes switched
    around. The first n / 2 existing samples become non samples, and the last
    n / 2 node become samples.
    """
    tables = ts.dump_tables()
    nodes = tables.nodes
    flags = nodes.flags
    oldest_parent = tables.edges.parent[-1]
    n = ts.sample_size
    flags[:n // 2] = 0
    flags[oldest_parent - n // 2: oldest_parent] = 1
    nodes.set_columns(flags, nodes.time)
    add_provenance(tables.provenances, "jiggle_samples")
    return msprime.load_tables(**tables.asdict())


def generate_site_mutations(tree, position, mu, site_table, mutation_table,
                            multiple_per_node=True):
    """
    Generates mutations for the site at the specified position on the specified
    tree. Mutations happen at rate mu along each branch. The site and mutation
    information are recorded in the specified tables.  Note that this records
    more than one mutation per edge.
    """
    assert tree.interval[0] <= position < tree.interval[1]
    states = {"A", "C", "G", "T"}
    state = random.choice(list(states))
    site_table.add_row(position, state)
    site = site_table.num_rows - 1
    stack = [(tree.root, state, msprime.NULL_MUTATION)]
    while len(stack) != 0:
        u, state, parent = stack.pop()
        if u != tree.root:
            branch_length = tree.branch_length(u)
            x = random.expovariate(mu)
            new_state = state
            while x < branch_length:
                new_state = random.choice(list(states - set(state)))
                if multiple_per_node and (state != new_state):
                    mutation_table.add_row(site, u, new_state, parent)
                    parent = mutation_table.num_rows - 1
                    state = new_state
                x += random.expovariate(mu)
            else:
                if (not multiple_per_node) and (state != new_state):
                    mutation_table.add_row(site, u, new_state, parent)
                    parent = mutation_table.num_rows - 1
                    state = new_state
        stack.extend(reversed([(v, state, parent) for v in tree.children(u)]))


def jukes_cantor(ts, num_sites, mu, multiple_per_node=True, seed=None):
    """
    Returns a copy of the specified tree sequence with Jukes-Cantor mutations
    applied at the specfied rate at the specifed number of sites. Site positions
    are chosen uniformly.
    """
    random.seed(seed)
    positions = [ts.sequence_length * random.random() for _ in range(num_sites)]
    positions.sort()
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    trees = ts.trees()
    t = next(trees)
    for position in positions:
        while position >= t.interval[1]:
            t = next(trees)
        generate_site_mutations(t, position, mu, tables.sites, tables.mutations,
                                multiple_per_node=multiple_per_node)
    add_provenance(tables.provenances, "jukes_cantor")
    new_ts = msprime.load_tables(**tables.asdict())
    return new_ts


def compute_mutation_parent(ts):
    """
    Compute the `parent` column of a MutationTable. Correct computation uses
    topological information in the nodes and edges, as well as the fact that
    each mutation must be listed after the mutation on whose background it
    occurred (i.e., its parent).

    :param TreeSequence ts: The tree sequence to compute for.  Need not
        have a valid mutation parent column.
    """
    mutation_parent = np.zeros(ts.num_mutations, dtype=np.int32) - 1
    # Maps nodes to the bottom mutation on each branch
    bottom_mutation = np.zeros(ts.num_nodes, dtype=np.int32) - 1
    for tree in ts.trees():
        for site in tree.sites():
            # Go forward through the mutations creating a mapping from the
            # mutations to the nodes. If we see more than one mutation
            # at a node, then these must be parents since we're assuming
            # they are in order.
            for mutation in site.mutations:
                if bottom_mutation[mutation.node] != msprime.NULL_MUTATION:
                    mutation_parent[mutation.id] = bottom_mutation[mutation.node]
                bottom_mutation[mutation.node] = mutation.id
            # There's no point in checking the first mutation, since this cannot
            # have a parent.
            for mutation in site.mutations[1:]:
                if mutation_parent[mutation.id] == msprime.NULL_MUTATION:
                    v = tree.parent(mutation.node)
                    # Traverse upwards until we find a another mutation or root.
                    while v != msprime.NULL_NODE \
                            and bottom_mutation[v] == msprime.NULL_MUTATION:
                        v = tree.parent(v)
                    if v != msprime.NULL_NODE:
                        mutation_parent[mutation.id] = bottom_mutation[v]
            # Reset the maps for the next site.
            for mutation in site.mutations:
                bottom_mutation[mutation.node] = msprime.NULL_MUTATION
            assert np.all(bottom_mutation == -1)
    return mutation_parent


def algorithm_T(ts):
    """
    Simple implementation of algorithm T from the PLOS paper, taking into
    account tree sequences with gaps and other complexities.
    """
    sequence_length = ts.sequence_length
    edges = list(ts.edges())
    M = len(edges)
    time = [ts.node(edge.parent).time for edge in edges]
    in_order = sorted(range(M), key=lambda j: (
        edges[j].left, time[j], edges[j].parent, edges[j].child))
    out_order = sorted(range(M), key=lambda j: (
        edges[j].right, -time[j], -edges[j].parent, -edges[j].child))
    j = 0
    k = 0
    left = 0
    parent = [-1 for _ in range(ts.num_nodes)]
    while j < M or left < sequence_length:
        while k < M and edges[out_order[k]].right == left:
            edge = edges[out_order[k]]
            parent[edge.child] = -1
            k += 1
        while j < M and edges[in_order[j]].left == left:
            edge = edges[in_order[j]]
            parent[edge.child] = edge.parent
            j += 1
        right = sequence_length
        if j < M:
            right = min(right, edges[in_order[j]].left)
        if k < M:
            right = min(right, edges[out_order[k]].right)
        yield (left, right), parent
        left = right
