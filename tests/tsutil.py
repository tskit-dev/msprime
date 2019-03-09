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
import json
import random

import numpy as np

import tskit.provenance as provenance
import msprime


def add_provenance(provenance_table, method_name):
    d = provenance.get_provenance_dict({"command": "tsutil.{}".format(method_name)})
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
    return t.tree_sequence()


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
    return tables.tree_sequence()


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
    new_ts = tables.tree_sequence()
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
    ts = tables.tree_sequence()
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
    return tables.tree_sequence()


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
    new_ts = tables.tree_sequence()
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


class LinkedTree(object):
    """
    Straightforward implementation of the quintuply linked tree for developing
    and testing the sample lists feature.

    NOTE: The interface is pretty awkward; it's not intended for anything other
    than testing.
    """
    def __init__(self, tree_sequence, tracked_samples=None):
        self.tree_sequence = tree_sequence
        num_nodes = tree_sequence.num_nodes
        # Quintuply linked tree.
        self.parent = [-1 for _ in range(num_nodes)]
        self.left_sib = [-1 for _ in range(num_nodes)]
        self.right_sib = [-1 for _ in range(num_nodes)]
        self.left_child = [-1 for _ in range(num_nodes)]
        self.right_child = [-1 for _ in range(num_nodes)]
        self.left_sample = [-1 for _ in range(num_nodes)]
        self.right_sample = [-1 for _ in range(num_nodes)]
        # This is too long, but it's convenient for printing.
        self.next_sample = [-1 for _ in range(num_nodes)]

        self.sample_index_map = [-1 for _ in range(num_nodes)]
        samples = tracked_samples
        if tracked_samples is None:
            samples = list(tree_sequence.samples())
        for j in range(len(samples)):
            u = samples[j]
            self.sample_index_map[u] = j
            self.left_sample[u] = j
            self.right_sample[u] = j

    def __str__(self):
        fmt = "{:<5}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}\n"
        s = fmt.format(
            "node", "parent", "lsib", "rsib", "lchild", "rchild",
            "nsamp", "lsamp", "rsamp")
        for u in range(self.tree_sequence.num_nodes):
            s += fmt.format(
                u, self.parent[u],
                self.left_sib[u], self.right_sib[u],
                self.left_child[u], self.right_child[u],
                self.next_sample[u], self.left_sample[u], self.right_sample[u])
        # Strip off trailing newline
        return s[:-1]

    def remove_edge(self, edge):
        p = edge.parent
        c = edge.child
        lsib = self.left_sib[c]
        rsib = self.right_sib[c]
        if lsib == -1:
            self.left_child[p] = rsib
        else:
            self.right_sib[lsib] = rsib
        if rsib == -1:
            self.right_child[p] = lsib
        else:
            self.left_sib[rsib] = lsib
        self.parent[c] = -1
        self.left_sib[c] = -1
        self.right_sib[c] = -1

    def insert_edge(self, edge):
        p = edge.parent
        c = edge.child
        assert self.parent[c] == -1, "contradictory edges"
        self.parent[c] = p
        u = self.right_child[p]
        if u == -1:
            self.left_child[p] = c
            self.left_sib[c] = -1
            self.right_sib[c] = -1
        else:
            self.right_sib[u] = c
            self.left_sib[c] = u
            self.right_sib[c] = -1
        self.right_child[p] = c

    def update_sample_list(self, parent):
        # This can surely be done more efficiently and elegantly. We are iterating
        # up the tree and iterating over all the siblings of the nodes we visit,
        # rebuilding the links as we go. This results in visiting the same nodes
        # over again, which if we have nodes with many siblings will surely be
        # expensive. Another consequence of the current approach is that the
        # next pointer contains an arbitrary value for the rightmost sample of
        # every root. This should point to NULL ideally, but it's quite tricky
        # to do in practise. It's easier to have a slightly uglier iteration
        # over samples.
        #
        # In the future it would be good have a more efficient version of this
        # algorithm using next and prev pointers that we keep up to date at all
        # times, and which we use to patch the lists together more efficiently.
        u = parent
        while u != -1:
            sample_index = self.sample_index_map[u]
            if sample_index != -1:
                self.right_sample[u] = self.left_sample[u]
            else:
                self.right_sample[u] = -1
                self.left_sample[u] = -1
            v = self.left_child[u]
            while v != -1:
                if self.left_sample[v] != -1:
                    assert self.right_sample[v] != -1
                    if self.left_sample[u] == -1:
                        self.left_sample[u] = self.left_sample[v]
                        self.right_sample[u] = self.right_sample[v]
                    else:
                        self.next_sample[self.right_sample[u]] = self.left_sample[v]
                        self.right_sample[u] = self.right_sample[v]
                v = self.right_sib[v]
            u = self.parent[u]

    def sample_lists(self):
        """
        Iterate over the the trees in this tree sequence, yielding the (left, right)
        interval tuples. The tree state is maintained internally.

        See note above about the cruddiness of this interface.
        """
        ts = self.tree_sequence
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

        while j < M or left < sequence_length:
            while k < M and edges[out_order[k]].right == left:
                edge = edges[out_order[k]]
                self.remove_edge(edge)
                self.update_sample_list(edge.parent)
                k += 1
            while j < M and edges[in_order[j]].left == left:
                edge = edges[in_order[j]]
                self.insert_edge(edge)
                self.update_sample_list(edge.parent)
                j += 1
            right = sequence_length
            if j < M:
                right = min(right, edges[in_order[j]].left)
            if k < M:
                right = min(right, edges[out_order[k]].right)
            yield left, right
            left = right


def mean_descendants(ts, reference_sets):
    """
    Returns the mean number of nodes from the specified reference sets
    where the node is ancestral to at least one of the reference nodes. Returns a
    ``(ts.num_nodes, len(reference_sets))`` dimensional numpy array.
    """
    # Check the inputs (could be done more efficiently here)
    all_reference_nodes = set()
    for reference_set in reference_sets:
        U = set(reference_set)
        if len(U) != len(reference_set):
            raise ValueError("Cannot have duplicate values within set")
        if len(all_reference_nodes & U) != 0:
            raise ValueError("Sample sets must be disjoint")
        all_reference_nodes |= U

    K = len(reference_sets)
    C = np.zeros((ts.num_nodes, K))
    parent = np.zeros(ts.num_nodes, dtype=int) - 1
    # The -1th element of ref_count is for all nodes in the reference set.
    ref_count = np.zeros((ts.num_nodes, K + 1), dtype=int)
    last_update = np.zeros(ts.num_nodes)
    total_length = np.zeros(ts.num_nodes)

    def update_counts(edge, sign):
        # Update the counts and statistics for a given node. Before we change the
        # node counts in the given direction, check to see if we need to update
        # statistics for that node. When a node count changes, we add the
        # accumulated statistic value for the span since that node was last updated.
        v = edge.parent
        while v != -1:
            if last_update[v] != left:
                if ref_count[v, K] > 0:
                    length = left - last_update[v]
                    C[v] += length * ref_count[v, :K]
                    total_length[v] += length
                last_update[v] = left
            ref_count[v] += sign * ref_count[edge.child]
            v = parent[v]

    # Set the intitial conditions.
    for j in range(K):
        ref_count[reference_sets[j], j] = 1
    ref_count[ts.samples(), K] = 1

    for (left, right), edges_out, edges_in in ts.edge_diffs():
        for edge in edges_out:
            parent[edge.child] = -1
            update_counts(edge, -1)
        for edge in edges_in:
            parent[edge.child] = edge.parent
            update_counts(edge, +1)

    # Finally, add the stats for the last tree and divide by the total
    # length that each node was an ancestor to > 0 samples.
    for v in range(ts.num_nodes):
        if ref_count[v, K] > 0:
            length = ts.sequence_length - last_update[v]
            total_length[v] += length
            C[v] += length * ref_count[v, :K]
        if total_length[v] != 0:
            C[v] /= total_length[v]
    return C


def genealogical_nearest_neighbours(ts, focal, reference_sets):

    reference_set_map = np.zeros(ts.num_nodes, dtype=int) - 1
    for k, reference_set in enumerate(reference_sets):
        for u in reference_set:
            if reference_set_map[u] != -1:
                raise ValueError("Duplicate value in reference sets")
            reference_set_map[u] = k

    K = len(reference_sets)
    A = np.zeros((len(focal), K))
    L = np.zeros(len(focal))
    parent = np.zeros(ts.num_nodes, dtype=int) - 1
    sample_count = np.zeros((ts.num_nodes, K), dtype=int)

    # Set the intitial conditions.
    for j in range(K):
        sample_count[reference_sets[j], j] = 1

    for (left, right), edges_out, edges_in in ts.edge_diffs():
        for edge in edges_out:
            parent[edge.child] = -1
            v = edge.parent
            while v != -1:
                sample_count[v] -= sample_count[edge.child]
                v = parent[v]
        for edge in edges_in:
            parent[edge.child] = edge.parent
            v = edge.parent
            while v != -1:
                sample_count[v] += sample_count[edge.child]
                v = parent[v]

        # Process this tree.
        for j, u in enumerate(focal):
            focal_reference_set = reference_set_map[u]
            p = parent[u]
            while p != msprime.NULL_NODE:
                total = np.sum(sample_count[p])
                if total > 1:
                    break
                p = parent[p]
            if p != msprime.NULL_NODE:
                length = right - left
                L[j] += length
                scale = length / (total - int(focal_reference_set != -1))
                for k, reference_set in enumerate(reference_sets):
                    n = sample_count[p, k] - int(focal_reference_set == k)
                    A[j, k] += n * scale

    # Avoid division by zero
    L[L == 0] = 1
    A /= L.reshape((len(focal), 1))
    return A
