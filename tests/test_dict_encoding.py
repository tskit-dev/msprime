"""
Test cases for the low-level dictionary encoding used to move
data around in C.

The code for doing this was copied from the the _tskitmodule.c file
into the msprime C module. These tests were also copied in to ensure
that no errors are introduced into the code over time.
"""
import numpy as np
import pytest
import tskit

import msprime
import msprime._msprime as c_module


def get_example_tables():
    """
    Return a tree sequence that has data in all fields.
    """
    pop_configs = [msprime.PopulationConfiguration(5) for _ in range(2)]
    migration_matrix = [[0, 1], [1, 0]]
    ts = msprime.simulate(
        population_configurations=pop_configs,
        migration_matrix=migration_matrix,
        mutation_rate=1,
        record_migrations=True,
        random_seed=1,
    )

    tables = ts.dump_tables()
    for j in range(ts.num_samples):
        tables.individuals.add_row(flags=j, location=np.arange(j), metadata=b"x" * j)
    tables.nodes.clear()
    for node in ts.nodes():
        tables.nodes.add_row(
            flags=node.flags,
            time=node.time,
            population=node.population,
            individual=node.id if node.id < ts.num_samples else -1,
            metadata=b"y" * node.id,
        )
    tables.edges.clear()
    for edge in ts.edges():
        tables.edges.add_row(
            left=edge.left,
            right=edge.right,
            child=edge.child,
            parent=edge.parent,
            metadata=b"y" * edge.id,
        )
    tables.sites.clear()
    for site in ts.sites():
        tables.sites.add_row(
            position=site.position,
            ancestral_state="A" * site.id,
            metadata=b"q" * site.id,
        )
    tables.mutations.clear()
    for mutation in ts.mutations():
        mut_id = tables.mutations.add_row(
            site=mutation.site,
            node=mutation.node,
            time=0,
            parent=-1,
            derived_state="C" * mutation.id,
            metadata=b"x" * mutation.id,
        )
        # Add another mutation on the same branch.
        tables.mutations.add_row(
            site=mutation.site,
            node=mutation.node,
            time=0,
            parent=mut_id,
            derived_state="G" * mutation.id,
            metadata=b"y" * mutation.id,
        )
    tables.migrations.clear()
    for migration in ts.migrations():
        tables.migrations.add_row(
            left=migration.left,
            right=migration.right,
            node=migration.node,
            source=migration.source,
            dest=migration.dest,
            time=migration.time,
            metadata=b"y" * migration.id,
        )
    for j in range(10):
        tables.populations.add_row(metadata=b"p" * j)
        tables.provenances.add_row(timestamp="x" * j, record="y" * j)
    tables.metadata_schema = tskit.MetadataSchema(
        {
            "codec": "struct",
            "type": "object",
            "properties": {
                "top-level": {
                    "type": "array",
                    "items": {"type": "integer", "binaryFormat": "B"},
                    "noLengthEncodingExhaustBuffer": True,
                }
            },
        }
    )
    tables.metadata = {"top-level": [1, 2, 3, 4]}
    for table in [
        "individuals",
        "nodes",
        "edges",
        "migrations",
        "sites",
        "mutations",
        "populations",
    ]:
        t = getattr(tables, table)
        t.metadata_schema = tskit.MetadataSchema(
            {
                "codec": "struct",
                "type": "object",
                "properties": {table: {"type": "string", "binaryFormat": "50p"}},
            }
        )
    return tables


class TestEncodingVersion:
    def test_version(self):
        lwt = c_module.LightweightTableCollection()
        assert lwt.asdict()["encoding_version"] == (1, 1)


class TestRoundTrip:
    """
    Tests if we can do a simple round trip on simulated data.
    """

    def verify(self, tables):
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(tables.asdict())
        other_tables = tskit.TableCollection.fromdict(lwt.asdict())
        assert tables == other_tables

    def test_simple(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=2)
        self.verify(ts.tables)

    def test_empty(self):
        tables = tskit.TableCollection(sequence_length=1)
        self.verify(tables)

    def test_individuals(self):
        n = 10
        ts = msprime.simulate(n, mutation_rate=1, random_seed=2)
        tables = ts.dump_tables()
        for j in range(n):
            tables.individuals.add_row(flags=j, location=(j, j), metadata=b"x" * j)
        self.verify(tables)

    def test_sequence_length(self):
        ts = msprime.simulate(
            10, recombination_rate=0.1, mutation_rate=1, length=0.99, random_seed=2
        )
        self.verify(ts.tables)

    def test_migration(self):
        pop_configs = [msprime.PopulationConfiguration(5) for _ in range(2)]
        migration_matrix = [[0, 1], [1, 0]]
        ts = msprime.simulate(
            population_configurations=pop_configs,
            migration_matrix=migration_matrix,
            mutation_rate=1,
            record_migrations=True,
            random_seed=1,
        )
        self.verify(ts.tables)

    def test_example(self):
        tables = get_example_tables()
        tables.metadata_schema = tskit.MetadataSchema(
            {
                "codec": "struct",
                "type": "object",
                "properties": {"top-level": {"type": "string", "binaryFormat": "50p"}},
            }
        )
        tables.metadata = {"top-level": "top-level-metadata"}
        for table in [
            "individuals",
            "nodes",
            "edges",
            "migrations",
            "sites",
            "mutations",
            "populations",
        ]:
            t = getattr(tables, table)
            t.packset_metadata([f"{table}-{i}".encode() for i in range(t.num_rows)])
            t.metadata_schema = tskit.MetadataSchema(
                {
                    "codec": "struct",
                    "type": "object",
                    "properties": {table: {"type": "string", "binaryFormat": "50p"}},
                }
            )

        self.verify(tables)


class TestMissingData:
    """
    Tests what happens when we have missing data in the encoded dict.
    """

    def test_missing_sequence_length(self):
        tables = get_example_tables()
        d = tables.asdict()
        del d["sequence_length"]
        lwt = c_module.LightweightTableCollection()
        with pytest.raises(TypeError):
            lwt.fromdict(d)

    def test_missing_metadata(self):
        tables = get_example_tables()
        assert tables.metadata != b""
        d = tables.asdict()
        del d["metadata"]
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(d)
        tables = tskit.TableCollection.fromdict(lwt.asdict())
        # Empty byte field still gets interpreted by schema
        assert tables.metadata == {"top-level": []}

    def test_missing_metadata_schema(self):
        tables = get_example_tables()
        assert str(tables.metadata_schema) != ""
        d = tables.asdict()
        del d["metadata_schema"]
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(d)
        tables = tskit.TableCollection.fromdict(lwt.asdict())
        assert str(tables.metadata_schema) == ""

    def test_missing_tables(self):
        tables = get_example_tables()
        d = tables.asdict()
        table_names = d.keys() - {
            "sequence_length",
            "metadata",
            "metadata_schema",
            "encoding_version",
        }
        for table_name in table_names:
            d = tables.asdict()
            del d[table_name]
            lwt = c_module.LightweightTableCollection()
            with pytest.raises(TypeError):
                lwt.fromdict(d)


class TestBadTypes:
    """
    Tests for setting each column to a type that can't be converted to 1D numpy array.
    """

    def verify_columns(self, value):
        tables = get_example_tables()
        d = tables.asdict()
        table_names = set(d.keys()) - {
            "sequence_length",
            "metadata",
            "metadata_schema",
            "encoding_version",
        }
        for table_name in table_names:
            table_dict = d[table_name]
            for colname in set(table_dict.keys()) - {"metadata_schema"}:
                copy = dict(table_dict)
                copy[colname] = value
                lwt = c_module.LightweightTableCollection()
                d = tables.asdict()
                d[table_name] = copy
                with pytest.raises(ValueError):
                    lwt.fromdict(d)

    def test_2d_array(self):
        self.verify_columns([[1, 2], [3, 4]])

    def test_str(self):
        self.verify_columns("aserg")

    def test_bad_top_level_types(self):
        tables = get_example_tables()
        d = tables.asdict()
        for key in set(d.keys()) - {"encoding_version"}:
            bad_type_dict = tables.asdict()
            # A list should be a ValueError for both the tables and sequence_length
            bad_type_dict[key] = ["12345"]
            lwt = c_module.LightweightTableCollection()
            with pytest.raises(TypeError):
                lwt.fromdict(bad_type_dict)


class TestBadLengths:
    """
    Tests for setting each column to a length incompatible with the table.
    """

    def verify(self, num_rows):

        tables = get_example_tables()
        d = tables.asdict()
        table_names = set(d.keys()) - {
            "sequence_length",
            "metadata",
            "metadata_schema",
            "encoding_version",
        }
        for table_name in sorted(table_names):
            table_dict = d[table_name]
            for colname in set(table_dict.keys()) - {"metadata_schema"}:
                copy = dict(table_dict)
                copy[colname] = table_dict[colname][:num_rows].copy()
                lwt = c_module.LightweightTableCollection()
                d = tables.asdict()
                d[table_name] = copy
                with pytest.raises(ValueError):
                    lwt.fromdict(d)

    def test_two_rows(self):
        self.verify(2)

    def test_zero_rows(self):
        self.verify(0)


class TestRequiredAndOptionalColumns:
    """
    Tests that specifying None for some columns will give the intended
    outcome.
    """

    def verify_required_columns(self, tables, table_name, required_cols):
        d = tables.asdict()
        table_dict = {col: None for col in d[table_name].keys()}
        for col in required_cols:
            table_dict[col] = d[table_name][col]
        lwt = c_module.LightweightTableCollection()
        d[table_name] = table_dict
        lwt.fromdict(d)
        other = lwt.asdict()
        for col in required_cols:
            assert np.array_equal(other[table_name][col], table_dict[col])

        # Any one of these required columns as None gives an error.
        for col in required_cols:
            d = tables.asdict()
            copy = dict(table_dict)
            copy[col] = None
            d[table_name] = copy
            lwt = c_module.LightweightTableCollection()
            with pytest.raises(TypeError):
                lwt.fromdict(d)

        # Removing any one of these required columns gives an error.
        for col in required_cols:
            d = tables.asdict()
            copy = dict(table_dict)
            del copy[col]
            d[table_name] = copy
            lwt = c_module.LightweightTableCollection()
            with pytest.raises(TypeError):
                lwt.fromdict(d)

    def verify_optional_column(self, tables, table_len, table_name, col_name):
        d = tables.asdict()
        table_dict = d[table_name]
        table_dict[col_name] = None
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(d)
        out = lwt.asdict()
        assert np.array_equal(
            out[table_name][col_name], np.zeros(table_len, dtype=np.int32) - 1
        )

    def verify_offset_pair(
        self, tables, table_len, table_name, col_name, required=False
    ):
        offset_col = col_name + "_offset"

        if not required:
            d = tables.asdict()
            table_dict = d[table_name]
            table_dict[col_name] = None
            table_dict[offset_col] = None
            lwt = c_module.LightweightTableCollection()
            lwt.fromdict(d)
            out = lwt.asdict()
            assert out[table_name][col_name].shape == (0,)
            assert np.array_equal(
                out[table_name][offset_col],
                np.zeros(table_len + 1, dtype=np.uint32),
            )
            d = tables.asdict()
            table_dict = d[table_name]
            del table_dict[col_name]
            del table_dict[offset_col]
            lwt = c_module.LightweightTableCollection()
            lwt.fromdict(d)
            out = lwt.asdict()
            assert out[table_name][col_name].shape == (0,)
            assert np.array_equal(
                out[table_name][offset_col],
                np.zeros(table_len + 1, dtype=np.uint32),
            )

        # Setting one or the other raises a TypeError
        d = tables.asdict()
        table_dict = d[table_name]
        table_dict[col_name] = None
        lwt = c_module.LightweightTableCollection()
        with pytest.raises(TypeError):
            lwt.fromdict(d)

        d = tables.asdict()
        table_dict = d[table_name]
        del table_dict[col_name]
        lwt = c_module.LightweightTableCollection()
        with pytest.raises(TypeError):
            lwt.fromdict(d)

        d = tables.asdict()
        table_dict = d[table_name]
        table_dict[offset_col] = None
        lwt = c_module.LightweightTableCollection()
        with pytest.raises(TypeError):
            lwt.fromdict(d)

        d = tables.asdict()
        table_dict = d[table_name]
        del table_dict[offset_col]
        lwt = c_module.LightweightTableCollection()
        with pytest.raises(TypeError):
            lwt.fromdict(d)

        d = tables.asdict()
        table_dict = d[table_name]
        bad_offset = np.zeros_like(table_dict[offset_col])
        bad_offset[:-1] = table_dict[offset_col][:-1][::-1]
        bad_offset[-1] = table_dict[offset_col][-1]
        table_dict[offset_col] = bad_offset
        lwt = c_module.LightweightTableCollection()
        with pytest.raises(c_module.LibraryError):
            lwt.fromdict(d)

    def verify_metadata_schema(self, tables, table_name):
        d = tables.asdict()
        d[table_name]["metadata_schema"] = None
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(d)
        out = lwt.asdict()
        assert "metadata_schema" not in out[table_name]
        tables = tskit.TableCollection.fromdict(out)
        assert str(getattr(tables, table_name).metadata_schema) == ""

    def test_individuals(self):
        tables = get_example_tables()
        self.verify_required_columns(tables, "individuals", ["flags"])
        self.verify_offset_pair(
            tables, len(tables.individuals), "individuals", "location"
        )
        self.verify_offset_pair(
            tables, len(tables.individuals), "individuals", "metadata"
        )
        self.verify_metadata_schema(tables, "individuals")

    def test_nodes(self):
        tables = get_example_tables()
        self.verify_offset_pair(tables, len(tables.nodes), "nodes", "metadata")
        self.verify_optional_column(tables, len(tables.nodes), "nodes", "population")
        self.verify_optional_column(tables, len(tables.nodes), "nodes", "individual")
        self.verify_required_columns(tables, "nodes", ["flags", "time"])
        self.verify_metadata_schema(tables, "nodes")

    def test_edges(self):
        tables = get_example_tables()
        self.verify_required_columns(
            tables, "edges", ["left", "right", "parent", "child"]
        )
        self.verify_offset_pair(tables, len(tables.edges), "edges", "metadata")
        self.verify_metadata_schema(tables, "edges")

    def test_migrations(self):
        tables = get_example_tables()
        self.verify_required_columns(
            tables, "migrations", ["left", "right", "node", "source", "dest", "time"]
        )
        self.verify_offset_pair(
            tables, len(tables.migrations), "migrations", "metadata"
        )
        self.verify_optional_column(tables, len(tables.nodes), "nodes", "individual")
        self.verify_metadata_schema(tables, "migrations")

    def test_sites(self):
        tables = get_example_tables()
        self.verify_required_columns(
            tables, "sites", ["position", "ancestral_state", "ancestral_state_offset"]
        )
        self.verify_offset_pair(tables, len(tables.sites), "sites", "metadata")
        self.verify_metadata_schema(tables, "sites")

    def test_mutations(self):
        tables = get_example_tables()
        self.verify_required_columns(
            tables,
            "mutations",
            ["site", "node", "derived_state", "derived_state_offset"],
        )
        self.verify_offset_pair(tables, len(tables.mutations), "mutations", "metadata")
        self.verify_metadata_schema(tables, "mutations")
        # Verify optional time column
        d = tables.asdict()
        d["mutations"]["time"] = None
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(d)
        out = lwt.asdict()
        assert all(np.isnan(val) for val in out["mutations"]["time"])

    def test_populations(self):
        tables = get_example_tables()
        self.verify_required_columns(
            tables, "populations", ["metadata", "metadata_offset"]
        )
        self.verify_metadata_schema(tables, "populations")
        self.verify_offset_pair(tables, len(tables.nodes), "nodes", "metadata", True)

    def test_provenances(self):
        tables = get_example_tables()
        self.verify_required_columns(
            tables,
            "provenances",
            ["record", "record_offset", "timestamp", "timestamp_offset"],
        )

    def test_top_level_metadata(self):
        tables = get_example_tables()
        d = tables.asdict()
        # None should give default value
        d["metadata"] = None
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(d)
        out = lwt.asdict()
        assert "metadata" not in out
        tables = tskit.TableCollection.fromdict(out)
        # We only removed the metadata, not the schema. So empty bytefield
        # still gets interpreted
        assert tables.metadata == {"top-level": []}
        # Missing is tested in TestMissingData above

    def test_top_level_metadata_schema(self):
        tables = get_example_tables()
        d = tables.asdict()
        # None should give default value
        d["metadata_schema"] = None
        lwt = c_module.LightweightTableCollection()
        lwt.fromdict(d)
        out = lwt.asdict()
        assert "metadata_schema" not in out
        tables = tskit.TableCollection.fromdict(out)
        assert str(tables.metadata_schema) == ""
        # Missing is tested in TestMissingData above
