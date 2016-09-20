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

try:
    import h5py
    _h5py_imported = True
    # Numpy is required by h5py, so we can safely use it.
    import numpy as np
except ImportError:
    _h5py_imported = False


import msprime
import _msprime


def _get_provenance(command, attrs):
    """
    Returns the V2 tree provenance attributes reformatted as a V3
    provenance string.
    """
    environment = json.loads(attrs["environment"])
    parameters = json.loads(attrs["parameters"])
    provenance = msprime.get_provenance_dict(command, parameters)
    provenance["version"] = environment["msprime_version"]
    provenance["environment"] = environment
    return json.dumps(provenance)


def _load_legacy_hdf5(root):
    if 'format_version' not in root.attrs:
        raise ValueError("HDF5 file not in msprime format")
    format_version = root.attrs['format_version']
    if format_version[0] != 2:
        raise ValueError("Only version 2.x are supported")

    # Get the coalescence records
    trees_group = root["trees"]
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

    # Get the mutations (if present)
    mutations = None
    if "mutations" in root:
        mutations_group = root["mutations"]
        position = np.array(mutations_group["position"])
        node = np.array(mutations_group["node"])
        num_mutations = len(node)
        mutations = num_mutations * [None]
        for j in range(num_mutations):
            mutations[j] = msprime.Mutation(
                position=position[j], node=node[j], index=j)

    ll_ts = _msprime.TreeSequence()
    if samples is None:
        ll_ts.load_records(records)
    else:
        ll_ts.load_records(records, samples)
    ll_ts.add_provenance_string(
        _get_provenance("generate_trees", trees_group.attrs))
    if mutations is not None:
        ll_ts.set_mutations(mutations)
        ll_ts.add_provenance_string(
            _get_provenance("generate_mutations", mutations_group.attrs))
    ll_ts.add_provenance_string(
        json.dumps(msprime.get_provenance_dict("upgrade", {})))
    return ll_ts


def load_legacy(filename):
    """
    Reads the specified msprime HDF5 file and returns a tree sequence. This
    method is only intended to be used to read old format HDF5 files.
    """
    if not _h5py_imported:
        raise RuntimeError("h5py is required for converting HDF5 files.")
    root = h5py.File(filename, "r")
    try:
        ll_ts = _load_legacy_hdf5(root)
    finally:
        root.close()
    return msprime.TreeSequence(ll_ts)


def _dump_legacy_hdf5(tree_sequence, root):
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
            node.append(mutation.node)
            position.append(mutation.position)
        l = len(node)
        mutations = root.create_group("mutations")
        mutations.attrs["environment"] = json.dumps({"msprime_version": 0})
        mutations.attrs["parameters"] = "{}"
        mutations.create_dataset("position", (l, ), data=position, dtype=float)
        mutations.create_dataset("node", (l, ), data=node, dtype="u4")


def dump_legacy(tree_sequence, filename, version=2):
    """
    Writes the specified tree sequence to a HDF5 file in the specified
    legacy file format version.
    """
    if version != 2:
        raise ValueError("Only version 2 file format is supported")
    if not _h5py_imported:
        raise RuntimeError("h5py is required for converting HDF5 files.")
    root = h5py.File(filename, "w")
    try:
        _dump_legacy_hdf5(tree_sequence, root)
    finally:
        root.close()
