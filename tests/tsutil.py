#
# Copyright (C) 2017-2020 University of Oxford
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
import tskit


def add_provenance(provenance_table, method_name):
    d = tskit.provenance.get_provenance_dict({"command": f"tsutil.{method_name}"})
    provenance_table.add_row(json.dumps(d))


def decapitate(ts, num_edges):
    """
    Returns a copy of the specified tree sequence in which the specified number of
    edges have been retained.
    """
    t = ts.dump_tables()
    t.edges.set_columns(
        left=t.edges.left[:num_edges],
        right=t.edges.right[:num_edges],
        parent=t.edges.parent[:num_edges],
        child=t.edges.child[:num_edges],
    )
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
        site = tables.sites.add_row(position=tree.interval[0], ancestral_state="0")
        for root in tree.roots:
            state = {root: 0}
            mutation = {root: -1}
            stack = [root]
            while len(stack) > 0:
                u = stack.pop()
                stack.extend(tree.children(u))
                v = tree.parent(u)
                if v != tskit.NULL:
                    state[u] = state[v]
                    parent = mutation[v]
                    for _ in range(mutations_per_branch):
                        state[u] = (state[u] + 1) % 2
                        mutation[u] = tables.mutations.add_row(
                            site=site,
                            node=u,
                            derived_state=str(state[u]),
                            parent=parent,
                        )
                        parent = mutation[u]
    add_provenance(tables.provenances, "insert_branch_mutations")
    return tables.tree_sequence()


def insert_site(ts):
    """
    Returns a copy of the specified tree sequence with a new site
    and no mutation.
    """
    tables = ts.dump_tables()
    tables.sites.add_row(position=ts.sequence_length / 2, ancestral_state="XX")
    add_provenance(tables.provenances, "insert_site")
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
            position=tree.interval[0], ancestral_state=ancestral_state
        )
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
    individual[:] = tskit.NULL
    while j < len(samples):
        ploidy = rng.randint(0, max_ploidy)
        nodes = samples[j : min(j + ploidy, len(samples))]
        dimension = rng.randint(0, max_dimension)
        location = [rng.random() for _ in range(dimension)]
        ind_id = tables.individuals.add_row(location=location)
        individual[nodes] = ind_id
        j += ploidy
    tables.nodes.individual = individual
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
        flags=nodes.flags,
        population=nodes.population,
        time=nodes.time,
        metadata_offset=offset,
        metadata=metadata,
        individual=nodes.individual,
    )

    length = np.random.randint(0, max_length, ts.num_sites)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    sites = tables.sites
    sites.set_columns(
        position=sites.position,
        ancestral_state=sites.ancestral_state,
        ancestral_state_offset=sites.ancestral_state_offset,
        metadata_offset=offset,
        metadata=metadata,
    )

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
        metadata_offset=offset,
        metadata=metadata,
    )

    length = np.random.randint(0, max_length, ts.num_individuals)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    individuals = tables.individuals
    individuals.set_columns(
        flags=individuals.flags,
        location=individuals.location,
        location_offset=individuals.location_offset,
        metadata_offset=offset,
        metadata=metadata,
    )

    length = np.random.randint(0, max_length, ts.num_populations)
    offset = np.cumsum(np.hstack(([0], length)), dtype=np.uint32)
    metadata = np.random.randint(-127, 127, offset[-1]).astype(np.int8)
    populations = tables.populations
    populations.set_columns(metadata_offset=offset, metadata=metadata)

    add_provenance(tables.provenances, "add_random_metadata")
    ts = tables.tree_sequence()
    return ts
