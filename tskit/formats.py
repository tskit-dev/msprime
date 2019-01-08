"""
Module responsible for converting tree sequence files from older
formats.
"""
from __future__ import division
from __future__ import print_function

import datetime
import json
import logging

import h5py
import numpy as np

import tskit
import tskit.provenance as provenance
import tskit.exceptions as exceptions


def _get_v2_provenance(command, attrs):
    """
    Returns the V2 tree provenance attributes reformatted as a provenance record.
    """
    environment = {}
    parameters = {}
    # Try to get the provenance strings. Malformed JSON should not prevent us
    # from finishing the conversion.
    try:
        environment = json.loads(str(attrs["environment"]))
    except ValueError:
        logging.warn("Failed to convert environment provenance")
    try:
        parameters = json.loads(str(attrs["parameters"]))
    except ValueError:
        logging.warn("Failed to convert parameters provenance")
    parameters["command"] = command
    provenance_dict = provenance.get_provenance_dict(parameters)
    provenance_dict["version"] = environment.get("msprime_version", "Unknown_version")
    provenance_dict["environment"] = environment
    return json.dumps(provenance_dict).encode()


def _get_upgrade_provenance(root):
    """
    Returns the provenance string from upgrading the specified HDF5 file.
    """
    # TODO add more parameters here like filename, etc.
    parameters = {
        "command": "upgrade",
        "source_version": list(map(int, root.attrs["format_version"]))
    }
    s = json.dumps(provenance.get_provenance_dict(parameters))
    return s.encode()


def _convert_hdf5_mutations(
        mutations_group, sites, mutations, remove_duplicate_positions):
    """
    Loads the v2/v3 into the specified tables.
    """
    position = np.array(mutations_group["position"])
    node = np.array(mutations_group["node"], dtype=np.int32)
    unique_position, index = np.unique(position, return_index=True)
    if unique_position.shape != position.shape:
        if remove_duplicate_positions:
            position = position[index]
            node = node[index]
        else:
            # TODO add the number of duplicates so that we can improve the
            # error message.
            raise exceptions.DuplicatePositionsError()
    num_mutations = position.shape[0]
    sites.set_columns(
        position=position,
        ancestral_state=ord("0") * np.ones(num_mutations, dtype=np.int8),
        ancestral_state_offset=np.arange(num_mutations + 1, dtype=np.uint32))
    mutations.set_columns(
        node=node,
        site=np.arange(num_mutations, dtype=np.int32),
        derived_state=ord("1") * np.ones(num_mutations, dtype=np.int8),
        derived_state_offset=np.arange(num_mutations + 1, dtype=np.uint32))


def _set_populations(tables):
    """
    Updates PopulationTable suitable to represent the populations referred to
    in the node table.
    """
    if len(tables.nodes) > 0:
        for _ in range(np.max(tables.nodes.population) + 1):
            tables.populations.add_row()


def _load_legacy_hdf5_v2(root, remove_duplicate_positions):
    # Get the coalescence records
    trees_group = root["trees"]
    old_timestamp = datetime.datetime.min.isoformat()
    provenances = tskit.ProvenanceTable()
    provenances.add_row(
        timestamp=old_timestamp,
        record=_get_v2_provenance("generate_trees", trees_group.attrs))
    num_rows = trees_group["node"].shape[0]
    index = np.arange(num_rows, dtype=int)
    parent = np.zeros(2 * num_rows, dtype=np.int32)
    parent[2 * index] = trees_group["node"]
    parent[2 * index + 1] = trees_group["node"]
    left = np.zeros(2 * num_rows, dtype=np.float64)
    left[2 * index] = trees_group["left"]
    left[2 * index + 1] = trees_group["left"]
    right = np.zeros(2 * num_rows, dtype=np.float64)
    right[2 * index] = trees_group["right"]
    right[2 * index + 1] = trees_group["right"]
    child = np.array(trees_group["children"], dtype=np.int32).flatten()

    tables = tskit.TableCollection(np.max(right))
    tables.edges.set_columns(left=left, right=right, parent=parent, child=child)

    cr_node = np.array(trees_group["node"], dtype=np.int32)
    num_nodes = max(np.max(child), np.max(cr_node)) + 1
    sample_size = np.min(cr_node)
    flags = np.zeros(num_nodes, dtype=np.uint32)
    population = np.zeros(num_nodes, dtype=np.int32)
    time = np.zeros(num_nodes, dtype=np.float64)
    flags[:sample_size] = tskit.NODE_IS_SAMPLE
    cr_population = np.array(trees_group["population"], dtype=np.int32)
    cr_time = np.array(trees_group["time"])
    time[cr_node] = cr_time
    population[cr_node] = cr_population
    if "samples" in root:
        samples_group = root["samples"]
        population[:sample_size] = samples_group["population"]
        if "time" in samples_group:
            time[:sample_size] = samples_group["time"]
    tables.nodes.set_columns(flags=flags, population=population, time=time)
    _set_populations(tables)

    if "mutations" in root:
        mutations_group = root["mutations"]
        _convert_hdf5_mutations(
            mutations_group, tables.sites, tables.mutations, remove_duplicate_positions)
        provenances.add_row(
            timestamp=old_timestamp,
            record=_get_v2_provenance("generate_mutations", mutations_group.attrs))
    tables.provenances.add_row(_get_upgrade_provenance(root))
    tables.sort()
    return tables.tree_sequence()


def _load_legacy_hdf5_v3(root, remove_duplicate_positions):
    # get the trees group for the records and samples
    trees_group = root["trees"]
    nodes_group = trees_group["nodes"]
    time = np.array(nodes_group["time"])

    breakpoints = np.array(trees_group["breakpoints"])
    records_group = trees_group["records"]
    left_indexes = np.array(records_group["left"])
    right_indexes = np.array(records_group["right"])
    record_node = np.array(records_group["node"], dtype=np.int32)
    num_nodes = time.shape[0]
    sample_size = np.min(record_node)
    flags = np.zeros(num_nodes, dtype=np.uint32)
    flags[:sample_size] = tskit.NODE_IS_SAMPLE

    children_length = np.array(records_group["num_children"], dtype=np.uint32)
    total_rows = np.sum(children_length)
    left = np.zeros(total_rows, dtype=np.float64)
    right = np.zeros(total_rows, dtype=np.float64)
    parent = np.zeros(total_rows, dtype=np.int32)
    record_left = breakpoints[left_indexes]
    record_right = breakpoints[right_indexes]
    k = 0
    for j in range(left_indexes.shape[0]):
        for _ in range(children_length[j]):
            left[k] = record_left[j]
            right[k] = record_right[j]
            parent[k] = record_node[j]
            k += 1
    tables = tskit.TableCollection(np.max(right))
    tables.nodes.set_columns(
        flags=flags,
        time=nodes_group["time"],
        population=nodes_group["population"])
    _set_populations(tables)
    tables.edges.set_columns(
        left=left, right=right, parent=parent, child=records_group["children"])
    if "mutations" in root:
        _convert_hdf5_mutations(
            root["mutations"], tables.sites, tables.mutations,
            remove_duplicate_positions)
    old_timestamp = datetime.datetime.min.isoformat()
    if "provenance" in root:
        for record in root["provenance"]:
            tables.provenances.add_row(timestamp=old_timestamp, record=record)
    tables.provenances.add_row(_get_upgrade_provenance(root))
    tables.sort()
    return tables.tree_sequence()


def load_legacy(filename, remove_duplicate_positions=False):
    """
    Reads the specified msprime HDF5 file and returns a tree sequence. This
    method is only intended to be used to read old format HDF5 files.

    If remove_duplicate_positions is True, remove all sites (except the
    first) that contain duplicate positions. If this is False, any input
    files that contain duplicate positions will raise an DuplicatePositionsError.
    """
    loaders = {
        2: _load_legacy_hdf5_v2,
        3: _load_legacy_hdf5_v3,
        10: _load_legacy_hdf5_v10,
    }
    root = h5py.File(filename, "r")
    if 'format_version' not in root.attrs:
        raise ValueError("HDF5 file not in msprime format")
    format_version = root.attrs['format_version']
    if format_version[0] not in loaders:
        raise ValueError("Version {} not supported for loading".format(format_version))
    try:
        ts = loaders[format_version[0]](root, remove_duplicate_positions)
    finally:
        root.close()
    return ts


def raise_hdf5_format_error(filename, original_exception):
    """
    Tries to open the specified file as a legacy HDF5 file. If it looks like
    an msprime format HDF5 file, raise an error advising to run tskit upgrade.
    """
    try:
        with h5py.File(filename, "r") as root:
            version = tuple(root.attrs["format_version"])
            raise exceptions.VersionTooOldError(
                "File format {} is too old. Please use the ``tskit upgrade`` command "
                "to upgrade this file to the latest version".format(version))
    except (IOError, OSError, KeyError):
        raise exceptions.FileFormatError(str(original_exception))


def _dump_legacy_hdf5_v2(tree_sequence, root):
    root.attrs["format_version"] = (2, 999)
    root.attrs["sample_size"] = tree_sequence.get_sample_size()
    root.attrs["sequence_length"] = tree_sequence.get_sequence_length(),
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
    length = len(time)
    trees = root.create_group("trees")
    trees.attrs["environment"] = json.dumps({"msprime_version": 0})
    trees.attrs["parameters"] = "{}"
    trees.create_dataset("left", (length, ), data=left, dtype=float)
    trees.create_dataset("right", (length, ), data=right, dtype=float)
    trees.create_dataset("time", (length, ), data=time, dtype=float)
    trees.create_dataset("node", (length, ), data=node, dtype="u4")
    trees.create_dataset("population", (length, ), data=population, dtype="u1")
    trees.create_dataset(
        "children", (length, 2), data=children, dtype="u4")
    samples = root.create_group("samples")
    population = []
    time = []
    length = tree_sequence.get_sample_size()
    for u in range(length):
        time.append(tree_sequence.get_time(u))
        population.append(tree_sequence.get_population(u))
    samples.create_dataset("time", (length, ), data=time, dtype=float)
    samples.create_dataset("population", (length, ), data=population, dtype="u1")
    if tree_sequence.get_num_mutations() > 0:
        node = []
        position = []
        for site in tree_sequence.sites():
            if len(site.mutations) != 1:
                raise ValueError("v2 does not support recurrent mutations")
            if site.ancestral_state != "0" or site.mutations[0].derived_state != "1":
                raise ValueError("v2 does not support non-binary mutations")
            position.append(site.position)
            node.append(site.mutations[0].node)
        length = len(node)
        mutations = root.create_group("mutations")
        mutations.attrs["environment"] = json.dumps({"msprime_version": 0})
        mutations.attrs["parameters"] = "{}"
        mutations.create_dataset("position", (length, ), data=position, dtype=float)
        mutations.create_dataset("node", (length, ), data=node, dtype="u4")


def _dump_legacy_hdf5_v3(tree_sequence, root):
    root.attrs["format_version"] = (3, 999)
    root.attrs["sample_size"] = 0,
    root.attrs["sequence_length"] = 0,
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
    length = len(num_children)
    records_group.create_dataset("left", (length, ), data=left, dtype="u4")
    records_group.create_dataset("right", (length, ), data=right, dtype="u4")
    records_group.create_dataset("node", (length, ), data=node, dtype="u4")
    records_group.create_dataset(
        "num_children", (length, ), data=num_children, dtype="u4")
    records_group.create_dataset(
        "children", (len(children), ), data=children, dtype="u4")

    indexes_group = trees.create_group("indexes")
    left_index = sorted(range(length), key=lambda j: (left[j], time[j]))
    right_index = sorted(range(length), key=lambda j: (right[j], -time[j]))
    indexes_group.create_dataset(
        "insertion_order", (length, ), data=left_index, dtype="u4")
    indexes_group.create_dataset(
        "removal_order", (length, ), data=right_index, dtype="u4")

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
    length = tree_sequence.num_nodes
    nodes_group.create_dataset("time", (length, ), data=time, dtype=float)
    nodes_group.create_dataset("population", (length, ), data=population, dtype="u4")

    node = []
    position = []
    for site in tree_sequence.sites():
        if len(site.mutations) != 1:
            raise ValueError("v3 does not support recurrent mutations")
        if site.ancestral_state != "0" or site.mutations[0].derived_state != "1":
            raise ValueError("v3 does not support non-binary mutations")
        position.append(site.position)
        node.append(site.mutations[0].node)
    length = len(position)
    if length > 0:
        mutations = root.create_group("mutations")
        mutations.create_dataset("position", (length, ), data=position, dtype=float)
        mutations.create_dataset("node", (length, ), data=node, dtype="u4")


def _add_dataset(group, name, data):
    # In the HDF5 format any zero-d arrays must be excluded.
    if data.shape[0] > 0:
        group.create_dataset(name, data=data)


def _dump_legacy_hdf5_v10(tree_sequence, root):
    root.attrs["format_version"] = (10, 999)
    root.attrs["sample_size"] = 0,
    root.attrs["sequence_length"] = tree_sequence.sequence_length,
    tables = tree_sequence.dump_tables()

    nodes = root.create_group("nodes")
    _add_dataset(nodes, "time", tables.nodes.time)
    _add_dataset(nodes, "flags", tables.nodes.flags)
    _add_dataset(nodes, "population", tables.nodes.population)
    _add_dataset(nodes, "metadata", tables.nodes.metadata)
    _add_dataset(nodes, "metadata_offset", tables.nodes.metadata_offset)

    edges = root.create_group("edges")
    if len(tables.edges) > 0:
        edges.create_dataset("left", data=tables.edges.left)
        edges.create_dataset("right", data=tables.edges.right)
        edges.create_dataset("parent", data=tables.edges.parent)
        edges.create_dataset("child", data=tables.edges.child)

        left = tables.edges.left
        right = tables.edges.right
        time = tables.nodes.time[tables.edges.parent]
        # We can do this more efficiently if we ever need to do it for anything
        # other than testing.
        indexes_group = edges.create_group("indexes")
        length = len(tables.edges)
        left_index = sorted(range(length), key=lambda j: (left[j], time[j]))
        right_index = sorted(range(length), key=lambda j: (right[j], -time[j]))
        indexes_group.create_dataset(
            "insertion_order", data=left_index, dtype="u4")
        indexes_group.create_dataset(
            "removal_order", data=right_index, dtype="u4")

    migrations = root.create_group("migrations")
    if len(tables.migrations) > 0:
        migrations.create_dataset("left", data=tables.migrations.left)
        migrations.create_dataset("right", data=tables.migrations.right)
        migrations.create_dataset("node", data=tables.migrations.node)
        migrations.create_dataset("source", data=tables.migrations.source)
        migrations.create_dataset("dest", data=tables.migrations.dest)
        migrations.create_dataset("time", data=tables.migrations.time)

    sites = root.create_group("sites")
    _add_dataset(sites, "position", tables.sites.position)
    _add_dataset(sites, "ancestral_state", tables.sites.ancestral_state)
    _add_dataset(sites, "ancestral_state_offset", tables.sites.ancestral_state_offset)
    _add_dataset(sites, "metadata", tables.sites.metadata)
    _add_dataset(sites, "metadata_offset", tables.sites.metadata_offset)

    mutations = root.create_group("mutations")
    _add_dataset(mutations, "site", tables.mutations.site)
    _add_dataset(mutations, "node", tables.mutations.node)
    _add_dataset(mutations, "parent", tables.mutations.parent)
    _add_dataset(mutations, "derived_state", tables.mutations.derived_state)
    _add_dataset(
        mutations, "derived_state_offset", tables.mutations.derived_state_offset)
    _add_dataset(mutations, "metadata", tables.mutations.metadata)
    _add_dataset(mutations, "metadata_offset", tables.mutations.metadata_offset)

    provenances = root.create_group("provenances")
    _add_dataset(provenances, "timestamp", tables.provenances.timestamp)
    _add_dataset(provenances, "timestamp_offset", tables.provenances.timestamp_offset)
    _add_dataset(provenances, "record", tables.provenances.record)
    _add_dataset(provenances, "record_offset", tables.provenances.record_offset)


def _load_legacy_hdf5_v10(root, remove_duplicate_positions=False):
    # We cannot have duplicate positions in v10, so this parameter is ignored
    sequence_length = root.attrs["sequence_length"]
    tables = tskit.TableCollection(sequence_length)

    nodes_group = root["nodes"]
    metadata = None
    metadata_offset = None
    if "metadata" in nodes_group:
        metadata = nodes_group["metadata"]
        metadata_offset = nodes_group["metadata_offset"]
    if "flags" in nodes_group:
        tables.nodes.set_columns(
            flags=nodes_group["flags"],
            population=nodes_group["population"],
            time=nodes_group["time"],
            metadata=metadata,
            metadata_offset=metadata_offset)

    edges_group = root["edges"]
    if "left" in edges_group:
        tables.edges.set_columns(
            left=edges_group["left"],
            right=edges_group["right"],
            parent=edges_group["parent"],
            child=edges_group["child"])

    migrations_group = root["migrations"]
    if "left" in migrations_group:
        tables.migrations.set_columns(
            left=migrations_group["left"],
            right=migrations_group["right"],
            node=migrations_group["node"],
            source=migrations_group["source"],
            dest=migrations_group["dest"],
            time=migrations_group["time"])

    sites_group = root["sites"]
    if "position" in sites_group:
        metadata = None
        metadata_offset = None
        if "metadata" in sites_group:
            metadata = sites_group["metadata"]
            metadata_offset = sites_group["metadata_offset"]
        tables.sites.set_columns(
            position=sites_group["position"],
            ancestral_state=sites_group["ancestral_state"],
            ancestral_state_offset=sites_group["ancestral_state_offset"],
            metadata=metadata,
            metadata_offset=metadata_offset)

    mutations_group = root["mutations"]
    if "site" in mutations_group:
        metadata = None
        metadata_offset = None
        if "metadata" in mutations_group:
            metadata = mutations_group["metadata"]
            metadata_offset = mutations_group["metadata_offset"]
        tables.mutations.set_columns(
            site=mutations_group["site"],
            node=mutations_group["node"],
            parent=mutations_group["parent"],
            derived_state=mutations_group["derived_state"],
            derived_state_offset=mutations_group["derived_state_offset"],
            metadata=metadata,
            metadata_offset=metadata_offset)

    provenances_group = root["provenances"]
    if "timestamp" in provenances_group:
        timestamp = provenances_group["timestamp"]
        timestamp_offset = provenances_group["timestamp_offset"]
        if "record" in provenances_group:
            record = provenances_group["record"]
            record_offset = provenances_group["record_offset"]
        else:
            record = np.empty_like(timestamp)
            record_offset = np.zeros_like(timestamp_offset)
        tables.provenances.set_columns(
            timestamp=timestamp,
            timestamp_offset=timestamp_offset,
            record=record,
            record_offset=record_offset)
    tables.provenances.add_row(_get_upgrade_provenance(root))
    _set_populations(tables)
    return tables.tree_sequence()


def dump_legacy(tree_sequence, filename, version=3):
    """
    Writes the specified tree sequence to a HDF5 file in the specified
    legacy file format version.
    """
    dumpers = {
        2: _dump_legacy_hdf5_v2,
        3: _dump_legacy_hdf5_v3,
        10: _dump_legacy_hdf5_v10,
    }
    if version not in dumpers:
        raise ValueError("Version {} file format is supported".format(version))
    root = h5py.File(filename, "w")
    try:
        dumpers[version](tree_sequence, root)
    finally:
        root.close()
