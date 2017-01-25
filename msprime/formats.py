#
# Copyright (C) 2016 Jerome Kelleher <jerome.kelleher@well.ox.ac.uk>
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
Module responsible for converting tree sequence files from older
formats.
"""
from __future__ import division
from __future__ import print_function

import json
import logging

try:
    import h5py
    _h5py_imported = True
    # Numpy is required by h5py, so we can safely use it.
    import numpy as np
except ImportError:
    _h5py_imported = False


import msprime
import _msprime


def _check_h5py():
    if not _h5py_imported:
        raise RuntimeError("h5py is required for converting HDF5 files.")


def _get_v2_provenance(command, attrs):
    """
    Returns the V2 tree provenance attributes reformatted as a V3
    provenance string.
    """
    environment = {}
    parameters = {}
    # Try to get the provenance strings. Malformed JSON should not prevent us
    # from finishing the conversion.
    try:
        environment = json.loads(attrs["environment"])
    except ValueError:
        logging.warn("Failed to convert environment provenance")
    try:
        parameters = json.loads(attrs["parameters"])
    except ValueError:
        logging.warn("Failed to convert parameters provenance")
    provenance = msprime.get_provenance_dict(command, parameters)
    provenance["version"] = environment.get("msprime_version", "Unknown_version")
    provenance["environment"] = environment
    return json.dumps(provenance).encode()


def _get_upgrade_provenance(root):
    """
    Returns the provenance string from upgrading the specified HDF5 file.
    """
    # TODO add more parameters here like filename, etc.
    parameters = {
        "source_version": list(map(int, root.attrs["format_version"]))
    }
    s = json.dumps(msprime.get_provenance_dict("upgrade", parameters))
    return s.encode()


def _load_legacy_hdf5_v2(root):
    # Get the coalescence records
    trees_group = root["trees"]
    provenance = [
        _get_v2_provenance("generate_trees", trees_group.attrs),
    ]
    left = np.array(trees_group["left"])
    right = np.array(trees_group["right"])
    node = np.array(trees_group["node"])
    children = np.array(trees_group["children"])
    population = np.array(trees_group["population"])
    time = np.array(trees_group["time"])
    num_records = len(left)
    records = num_records * [None]
    for j in range(num_records):
        records[j] = msprime.CoalescenceRecord(
            left=left[j], right=right[j], node=node[j],
            children=tuple(children[j]), time=time[j],
            population=population[j])

    # Get the samples (if present)
    samples = None
    if "samples" in root:
        samples_group = root["samples"]
        population = np.array(samples_group["population"])
        time = None
        if "time" in samples_group:
            time = np.array(samples_group["time"])
        sample_size = len(population)
        samples = sample_size * [None]
        for j in range(sample_size):
            t = 0
            if time is not None:
                t = time[j]
            samples[j] = msprime.Sample(population=population[j], time=t)
    else:
        sample_size = min(record.node for record in records)
        samples = [msprime.Sample(0, 0) for _ in range(sample_size)]

    # Get the mutations (if present)
    mutations = []
    if "mutations" in root:
        mutations_group = root["mutations"]
        provenance.append(
            _get_v2_provenance("generate_mutations", mutations_group.attrs))
        position = np.array(mutations_group["position"])
        node = np.array(mutations_group["node"])
        num_mutations = len(node)
        mutations = num_mutations * [None]
        for j in range(num_mutations):
            mutations[j] = msprime.Mutation(
                position=position[j], nodes=(node[j],), index=j)

    ll_ts = _msprime.TreeSequence()
    provenance.append(_get_upgrade_provenance(root))
    ll_ts.load_records(samples, records, mutations, provenance)
    return ll_ts


def _load_legacy_hdf5_v3(root):
    # get the trees group for the records and samples
    trees_group = root["trees"]
    nodes_group = trees_group["nodes"]
    time = np.array(nodes_group["time"])
    population = np.array(nodes_group["population"])

    breakpoints = np.array(trees_group["breakpoints"])
    records_group = trees_group["records"]
    left = np.array(records_group["left"])
    right = np.array(records_group["right"])
    record_node = np.array(records_group["node"])
    flat_children = np.array(records_group["children"])
    num_children = np.array(records_group["num_children"])
    num_records = left.shape[0]
    records = [None for _ in range(num_records)]
    offset = 0
    for j in range(num_records):
        children = tuple(flat_children[offset: offset + num_children[j]])
        offset += num_children[j]
        u = record_node[j]
        records[j] = msprime.CoalescenceRecord(
            left=breakpoints[left[j]], right=breakpoints[right[j]], node=u,
            children=children, time=time[u], population=population[u])

    sample_size = np.min(record_node)
    samples = [
        msprime.Sample(time=time[v], population=population[v])
        for v in range(sample_size)]

    mutations = []
    if "mutations" in root:
        mutations_group = root["mutations"]
        position = np.array(mutations_group["position"])
        node = np.array(mutations_group["node"])
        num_mutations = len(node)
        mutations = num_mutations * [None]
        for j in range(num_mutations):
            mutations[j] = msprime.Mutation(
                position=position[j], nodes=(node[j],), index=j)

    provenance = []
    if "provenance" in root:
        provenance = list(root["provenance"])
    provenance.append(_get_upgrade_provenance(root))
    ll_ts = _msprime.TreeSequence()
    ll_ts.load_records(samples, records, mutations, provenance)
    return ll_ts


def load_legacy(filename):
    """
    Reads the specified msprime HDF5 file and returns a tree sequence. This
    method is only intended to be used to read old format HDF5 files.
    """
    _check_h5py()
    loaders = {
        2: _load_legacy_hdf5_v2,
        3: _load_legacy_hdf5_v3,
    }
    root = h5py.File(filename, "r")
    if 'format_version' not in root.attrs:
        raise ValueError("HDF5 file not in msprime format")
    format_version = root.attrs['format_version']
    if format_version[0] not in loaders:
        raise ValueError("Version {} not supported for loading".format(format_version))
    try:
        ll_ts = loaders[format_version[0]](root)
    finally:
        root.close()
    return msprime.TreeSequence(ll_ts)


def _dump_legacy_hdf5_v2(tree_sequence, root):
    root.attrs["format_version"] = (2, 999)
    root.attrs["sample_size"] = tree_sequence.get_sample_size()
    root.attrs["sequence_length"] = tree_sequence.get_sequence_length()
    left = []
    right = []
    node = []
    children = []
    time = []
    population = []
    for record in tree_sequence.records():
        left.append(record.left)
        right.append(record.right)
        node.append(record.node)
        if len(record.children) != 2:
            raise ValueError("V2 files only support binary records")
        children.append(record.children)
        time.append(record.time)
        population.append(record.population)
    l = len(time)
    trees = root.create_group("trees")
    trees.attrs["environment"] = json.dumps({"msprime_version": 0})
    trees.attrs["parameters"] = "{}"
    trees.create_dataset("left", (l, ), data=left, dtype=float)
    trees.create_dataset("right", (l, ), data=right, dtype=float)
    trees.create_dataset("time", (l, ), data=time, dtype=float)
    trees.create_dataset("node", (l, ), data=node, dtype="u4")
    trees.create_dataset("population", (l, ), data=population, dtype="u1")
    trees.create_dataset(
        "children", (l, 2), data=children, dtype="u4")
    samples = root.create_group("samples")
    population = []
    time = []
    l = tree_sequence.get_sample_size()
    for u in range(l):
        time.append(tree_sequence.get_time(u))
        population.append(tree_sequence.get_population(u))
    samples.create_dataset("time", (l, ), data=time, dtype=float)
    samples.create_dataset("population", (l, ), data=population, dtype="u1")
    if tree_sequence.get_num_mutations() > 0:
        node = []
        position = []
        for mutation in tree_sequence.mutations():
            if len(mutation.nodes) != 1:
                raise ValueError("v2 does not support recurrent mutations")
            node.append(mutation.nodes[0])
            position.append(mutation.position)
        l = len(node)
        mutations = root.create_group("mutations")
        mutations.attrs["environment"] = json.dumps({"msprime_version": 0})
        mutations.attrs["parameters"] = "{}"
        mutations.create_dataset("position", (l, ), data=position, dtype=float)
        mutations.create_dataset("node", (l, ), data=node, dtype="u4")


def _dump_legacy_hdf5_v3(tree_sequence, root):
    root.attrs["format_version"] = (3, 999)
    root.attrs["sample_size"] = 0
    root.attrs["sequence_length"] = 0
    trees = root.create_group("trees")
    # Get the breakpoints from the records.
    left = [cr.left for cr in tree_sequence.records()]
    breakpoints = np.unique(left + [tree_sequence.sequence_length])
    trees.create_dataset(
        "breakpoints", (len(breakpoints), ), data=breakpoints, dtype=float)

    left = []
    right = []
    node = []
    children = []
    num_children = []
    time = []
    for cr in tree_sequence.records():
        node.append(cr.node)
        left.append(np.searchsorted(breakpoints, cr.left))
        right.append(np.searchsorted(breakpoints, cr.right))
        children.extend(cr.children)
        num_children.append(len(cr.children))
        time.append(cr.time)
    records_group = trees.create_group("records")
    l = len(num_children)
    records_group.create_dataset("left", (l, ), data=left, dtype="u4")
    records_group.create_dataset("right", (l, ), data=right, dtype="u4")
    records_group.create_dataset("node", (l, ), data=node, dtype="u4")
    records_group.create_dataset("num_children", (l, ), data=num_children, dtype="u4")
    records_group.create_dataset(
        "children", (len(children), ), data=children, dtype="u4")

    indexes_group = trees.create_group("indexes")
    I = sorted(range(l), key=lambda j: (left[j], time[j]))
    O = sorted(range(l), key=lambda j: (right[j], -time[j]))
    indexes_group.create_dataset("insertion_order", (l, ), data=I, dtype="u4")
    indexes_group.create_dataset("removal_order", (l, ), data=O, dtype="u4")

    nodes_group = trees.create_group("nodes")
    population = np.zeros(tree_sequence.num_nodes, dtype="u4")
    time = np.zeros(tree_sequence.num_nodes, dtype=float)
    tree = next(tree_sequence.trees())
    for u in range(tree_sequence.sample_size):
        population[u] = tree.population(u)
        time[u] = tree.time(u)
    for cr in tree_sequence.records():
        population[cr.node] = cr.population
        time[cr.node] = cr.time
    l = tree_sequence.num_nodes
    nodes_group.create_dataset("time", (l, ), data=time, dtype=float)
    nodes_group.create_dataset("population", (l, ), data=population, dtype="u4")

    node = []
    position = []
    for mutation in tree_sequence.mutations():
        if len(mutation.nodes) != 1:
            raise ValueError("Recurrent mutations not supported")
        node.append(mutation.nodes[0])
        position.append(mutation.position)
    l = len(position)
    if l > 0:
        mutations = root.create_group("mutations")
        mutations.create_dataset("position", (l, ), data=position, dtype=float)
        mutations.create_dataset("node", (l, ), data=node, dtype="u4")

    provenance = tree_sequence.get_provenance()
    if len(provenance) > 0:
        root.create_dataset("provenance", (len(provenance), ), data=provenance)


def dump_legacy(tree_sequence, filename, version=3):
    """
    Writes the specified tree sequence to a HDF5 file in the specified
    legacy file format version.
    """
    _check_h5py()
    dumpers = {
        2: _dump_legacy_hdf5_v2,
        3: _dump_legacy_hdf5_v3
    }
    if version not in dumpers:
        raise ValueError("Version {} file format is supported".format(version))
    root = h5py.File(filename, "w")
    try:
        dumpers[version](tree_sequence, root)
    finally:
        root.close()
    import shutil
    shutil.copyfile(filename, "tmp.hdf5")
