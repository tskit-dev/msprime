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
    return msprime.load_tables(
        nodes=t.nodes, edges=t.edges, sites=t.sites, mutations=t.mutations,
        provenances=t.provenances, sequence_length=ts.sequence_length)


def insert_branch_mutations(ts, mutations_per_branch=1):
    """
    Returns a copy of the specified tree sequence with a mutation on every branch
    in every tree.
    """
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    for tree in ts.trees():
        site = len(sites)
        sites.add_row(position=tree.interval[0], ancestral_state='0')
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
                        mutation[u] = len(mutations)
                        mutations.add_row(
                            site=site, node=u, derived_state=str(state[u]),
                            parent=parent)
                        parent = mutation[u]
    tables = ts.tables
    add_provenance(tables.provenances, "insert_branch_mutations")
    return msprime.load_tables(
        nodes=tables.nodes, edges=tables.edges, sites=sites, mutations=mutations,
        provenances=tables.provenances)


def insert_multichar_mutations(ts, seed=1, max_len=10):
    """
    Returns a copy of the specified tree sequence with multiple chararacter
    mutations on a randomly chosen branch in every tree.
    """
    rng = random.Random(seed)
    letters = ["A", "C", "T", "G"]
    sites = msprime.SiteTable()
    mutations = msprime.MutationTable()
    for tree in ts.trees():
        site = len(sites)
        ancestral_state = rng.choice(letters) * rng.randint(0, max_len)
        sites.add_row(position=tree.interval[0], ancestral_state=ancestral_state)
        u = rng.choice(list(tree.nodes()))
        derived_state = ancestral_state
        while ancestral_state == derived_state:
            derived_state = rng.choice(letters) * rng.randint(0, max_len)
        mutations.add_row(site=site, node=u, derived_state=derived_state)
    tables = ts.tables
    add_provenance(tables.provenances, "insert_multichar_mutations")
    return msprime.load_tables(
        nodes=tables.nodes, edges=tables.edges, sites=sites, mutations=mutations,
        provenances=tables.provenances)


def permute_nodes(ts, node_map):
    """
    Returns a copy of the specified tree sequence such that the nodes are
    permuted according to the specified map.
    """
    # Mapping from nodes in the new tree sequence back to nodes in the original
    reverse_map = [0 for _ in node_map]
    for j in range(ts.num_nodes):
        reverse_map[node_map[j]] = j
    old_nodes = list(ts.nodes())
    new_nodes = msprime.NodeTable()
    for j in range(ts.num_nodes):
        old_node = old_nodes[reverse_map[j]]
        new_nodes.add_row(
            flags=old_node.flags, metadata=old_node.metadata,
            population=old_node.population, time=old_node.time)
    new_edges = msprime.EdgeTable()
    for edge in ts.edges():
        new_edges.add_row(
            left=edge.left, right=edge.right, parent=node_map[edge.parent],
            child=node_map[edge.child])
    new_sites = msprime.SiteTable()
    new_mutations = msprime.MutationTable()
    for site in ts.sites():
        new_sites.add_row(
            position=site.position, ancestral_state=site.ancestral_state)
        for mutation in site.mutations:
            new_mutations.add_row(
                site=site.id, derived_state=mutation.derived_state,
                node=node_map[mutation.node])
    msprime.sort_tables(
        nodes=new_nodes, edges=new_edges, sites=new_sites, mutations=new_mutations)
    provenances = ts.dump_tables().provenances
    add_provenance(provenances, "permute_nodes")
    return msprime.load_tables(
        nodes=new_nodes, edges=new_edges, sites=new_sites, mutations=new_mutations,
        provenances=provenances)


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
    edges = tables.edges
    nodes = tables.nodes
    sites = tables.sites
    mutations = tables.mutations

    time = nodes.time[:]
    edges.reset()
    for edge in ts.edges():
        # Insert a new node in between the parent and child.
        u = len(nodes)
        t = time[edge.child] + (time[edge.parent] - time[edge.child]) / 2
        nodes.add_row(time=t)
        edges.add_row(
            left=edge.left, right=edge.right, parent=u, child=edge.child)
        edges.add_row(
            left=edge.left, right=edge.right, parent=edge.parent, child=u)
    msprime.sort_tables(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations)
    add_provenance(tables.provenances, "insert_redundant_breakpoints")
    new_ts = msprime.load_tables(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations,
        provenances=tables.provenances)
    return new_ts


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
        metadata_offset=offset, metadata=metadata)

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
    add_provenance(tables.provenances, "add_random_metadata")
    ts = msprime.load_tables(
        nodes=nodes, edges=tables.edges, sites=sites, mutations=mutations,
        provenances=tables.provenances, migrations=tables.migrations)
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
    return msprime.load_tables(
        nodes=nodes, edges=tables.edges, provenances=tables.provenances)


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
    sites = msprime.SiteTable(num_sites)
    mutations = msprime.MutationTable(num_sites)
    trees = ts.trees()
    t = next(trees)
    for position in positions:
        while position >= t.interval[1]:
            t = next(trees)
        generate_site_mutations(t, position, mu, sites, mutations,
                                multiple_per_node=multiple_per_node)
    tables = ts.dump_tables()
    add_provenance(tables.provenances, "jukes_cantor")
    new_ts = msprime.load_tables(
        nodes=tables.nodes, edges=tables.edges, sites=sites, mutations=mutations,
        provenances=tables.provenances)
    return new_ts


def compute_mutation_parent(ts):
    """
    (Re-)compute the `parent` column of a MutationTable. Doing this uses
    topological information in the nodes and edgesets, as well as the fact that
    each mutation must be listed after the mutation on whose background it occurred
    (i.e., its parent).

    :param TreeSequence ts: The tree sequence to compute for.  Need not
        have a valid mutation parent column.
    """
    if ts.num_mutations == 0:
        return []
    # sites are ordered by position,
    #  and mutations by site.
    # mutation_parent = np.repeat(-1, [ts.num_mutations])
    mutation_parent = [-1 for _ in range(ts.num_mutations)]
    for t in ts.trees():
        for site in t.sites():
            # If there is more than one mutation on a given node,
            # they will be in time-increasing order.
            node_map = {}
            for mut in site.mutations:
                u = mut.node
                while u != msprime.NULL_NODE and u not in node_map:
                    u = t.parent(u)
                if u != msprime.NULL_NODE:
                    mutation_parent[mut.id] = node_map[u].id
                    # # for checking, we would
                    # assert node_map[u].id == mut.parent
                node_map[mut.node] = mut
    return mutation_parent
