/*
 * MIT License
 *
 * Copyright (c) 2019-2022 Tskit Developers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "testlib.h"
#include <tskit/tables.h>

#include <float.h>
#include <unistd.h>
#include <stdlib.h>

static void
reverse_migrations(tsk_table_collection_t *tables)
{
    int ret;
    tsk_migration_table_t migrations;
    tsk_migration_t migration;
    tsk_id_t j, ret_id;

    /* Easy way to copy the metadata schema */
    ret = tsk_migration_table_copy(&tables->migrations, &migrations, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_clear(&migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = (tsk_id_t) tables->migrations.num_rows - 1; j >= 0; j--) {
        ret = tsk_migration_table_get_row(&tables->migrations, j, &migration);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret_id = tsk_migration_table_add_row(&migrations, migration.left,
            migration.right, migration.node, migration.source, migration.dest,
            migration.time, migration.metadata, migration.metadata_length);
        CU_ASSERT_FATAL(ret_id >= 0);
    }

    ret = tsk_migration_table_copy(&migrations, &tables->migrations, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_migration_table_free(&migrations);
}

static void
reverse_edges(tsk_table_collection_t *tables)
{
    int ret;
    tsk_edge_table_t edges;
    tsk_edge_t edge;
    tsk_id_t j, ret_id;

    /* Easy way to copy the metadata schema */
    ret = tsk_edge_table_copy(&tables->edges, &edges, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_clear(&edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = (tsk_id_t) tables->edges.num_rows - 1; j >= 0; j--) {
        ret = tsk_edge_table_get_row(&tables->edges, j, &edge);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret_id = tsk_edge_table_add_row(&edges, edge.left, edge.right, edge.parent,
            edge.child, edge.metadata, edge.metadata_length);
        CU_ASSERT_FATAL(ret_id >= 0);
    }

    ret = tsk_edge_table_copy(&edges, &tables->edges, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_edge_table_free(&edges);
}

static void
reverse_mutations(tsk_table_collection_t *tables)
{
    int ret;
    tsk_mutation_table_t mutations;
    tsk_mutation_t mutation;
    tsk_id_t j, ret_id;
    tsk_id_t new_parent;
    tsk_id_t n = (tsk_id_t) tables->mutations.num_rows;

    ret = tsk_mutation_table_init(&mutations, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = n - 1; j >= 0; j--) {
        ret = tsk_mutation_table_get_row(&tables->mutations, j, &mutation);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        new_parent = (mutation.parent == TSK_NULL) ? TSK_NULL : n - mutation.parent - 1;
        ret_id = tsk_mutation_table_add_row(&mutations, mutation.site, mutation.node,
            new_parent, mutation.time, mutation.derived_state,
            mutation.derived_state_length, mutation.metadata, mutation.metadata_length);
        CU_ASSERT_FATAL(ret_id >= 0);
    }

    ret = tsk_mutation_table_copy(&mutations, &tables->mutations, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_mutation_table_free(&mutations);
}

static void
insert_edge_metadata(tsk_table_collection_t *tables)
{
    int ret;
    tsk_edge_table_t edges;
    tsk_edge_t edge;
    tsk_id_t j, ret_id;
    char metadata[100];

    ret = tsk_edge_table_init(&edges, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) tables->edges.num_rows; j++) {
        ret = tsk_edge_table_get_row(&tables->edges, j, &edge);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        snprintf(metadata, sizeof(metadata), "md_%lld\n", (long long) j);
        ret_id = tsk_edge_table_add_row(&edges, edge.left, edge.right, edge.parent,
            edge.child, metadata, (tsk_size_t) strlen(metadata));
        CU_ASSERT_FATAL(ret_id >= 0);
    }
    ret = tsk_edge_table_copy(&edges, &tables->edges, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_edge_table_free(&edges);
}

static void
test_table_collection_equals_options(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tc1, tc2;

    char example_time_units[100] = "An example of time units with unicode ⏰";
    char example_metadata[100] = "An example of metadata with unicode 🎄🌳🌴🌲🎋";
    char example_metadata_schema[100]
        = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_time_units_length = (tsk_size_t) strlen(example_time_units);
    tsk_size_t example_metadata_length = (tsk_size_t) strlen(example_metadata);
    tsk_size_t example_metadata_schema_length
        = (tsk_size_t) strlen(example_metadata_schema);

    // Test equality empty tables
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_init(&tc2, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_TRUE(ret);

    // Adding some meat to the tables
    ret_id = tsk_node_table_add_row(&tc1.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(&tc1.nodes, TSK_NODE_IS_SAMPLE, 1.0, 0, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id
        = tsk_individual_table_add_row(&tc1.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tc1.populations, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tc1.edges, 0.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tc1.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tc1.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);

    // Equality of empty vs non-empty
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_FALSE(ret);
    ret = tsk_table_collection_copy(&tc1, &tc2, TSK_NO_INIT);
    CU_ASSERT_EQUAL(ret, 0);

    // Equivalent except for time_units
    ret = tsk_table_collection_set_metadata(
        &tc1, example_time_units, example_time_units_length);
    CU_ASSERT_EQUAL(ret, 0);

    // Equivalent except for metadata
    ret = tsk_table_collection_set_metadata(
        &tc1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_TS_METADATA);
    CU_ASSERT_TRUE(ret);
    /* TSK_CMP_IGNORE_METADATA implies TSK_CMP_IGNORE_TS_METADATA */
    ret = tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_METADATA);
    CU_ASSERT_TRUE(ret);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_FALSE(ret);
    ret = tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_PROVENANCE);
    CU_ASSERT_FALSE(ret);
    ret = tsk_table_collection_set_metadata(
        &tc2, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_TRUE(ret);
    ret = tsk_table_collection_set_metadata_schema(
        &tc1, example_metadata_schema, example_metadata_schema_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_TS_METADATA);
    CU_ASSERT_TRUE(ret);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_FALSE(ret);
    ret = tsk_table_collection_set_metadata_schema(
        &tc2, example_metadata_schema, example_metadata_schema_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_TRUE(ret);

    // Ignore provenance
    ret_id = tsk_provenance_table_add_row(&tc1.provenances, "time", 4, "record", 6);
    CU_ASSERT_EQUAL(ret_id, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_PROVENANCE);
    CU_ASSERT_TRUE(ret);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_FALSE(ret);
    ret = tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_TS_METADATA);
    CU_ASSERT_FALSE(ret);
    ret_id = tsk_provenance_table_add_row(&tc2.provenances, "time", 4, "record", 6);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_PROVENANCE);
    CU_ASSERT_TRUE(ret);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_TRUE(ret);

    // Ignore provenance timestamp
    ret_id = tsk_provenance_table_add_row(&tc1.provenances, "time", 4, "record", 6);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_provenance_table_add_row(&tc2.provenances, "other", 5, "record", 6);
    CU_ASSERT_FATAL(ret_id >= 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_PROVENANCE));
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_TIMESTAMPS));

    // Ignore provenance and top-level metadata.
    ret = tsk_provenance_table_clear(&tc1.provenances);
    CU_ASSERT_EQUAL(ret, 0);
    example_metadata[0] = 'J';
    ret = tsk_table_collection_set_metadata(
        &tc1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_equals(&tc1, &tc2, 0);
    CU_ASSERT_FALSE(ret);
    ret = tsk_table_collection_equals(
        &tc1, &tc2, TSK_CMP_IGNORE_TS_METADATA | TSK_CMP_IGNORE_PROVENANCE);
    CU_ASSERT_TRUE(ret);

    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);

    // Check what happens when one of the tables just differs by metadata.
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_population_table_add_row(&tc1.populations, "metadata", 8);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_population_table_add_row(&tc2.populations, "", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_METADATA));

    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);

    // Ignore tables
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata(
        &tc1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_set_metadata(
        &tc2, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    // Add one row for each table we're ignoring
    ret_id
        = tsk_individual_table_add_row(&tc1.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(&tc1.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tc1.edges, 0.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_migration_table_add_row(&tc1.migrations, 0, 0, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tc1.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tc1.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tc1.populations, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_TABLES));

    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);

    // Ignore reference sequence
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata(
        &tc1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_set_metadata(
        &tc2, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_reference_sequence_set_data(&tc1.reference_sequence, "A", 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    CU_ASSERT_TRUE(
        tsk_table_collection_equals(&tc1, &tc2, TSK_CMP_IGNORE_REFERENCE_SEQUENCE));

    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
}

static void
test_table_collection_simplify_errors(void)
{
    int ret;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    tsk_id_t ret_id;
    const char *individuals = "1      0.25     -2\n";
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret_id = tsk_site_table_add_row(&tables.sites, 0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = tsk_table_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SITE_POSITION);

    /* Out of order positions */
    tables.sites.position[0] = 0.5;
    ret = tsk_table_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_SITES);

    /* Position out of bounds */
    tables.sites.position[0] = 1.5;
    ret = tsk_table_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SITE_POSITION);
    tsk_site_table_truncate(&tables.sites, 0);
    tables.sites.position[0] = 0;

    /* Individual out of bounds */
    parse_individuals(individuals, &tables.individuals);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.num_rows, 1);
    ret = tsk_table_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    /* TODO More tests for this: see
     * https://github.com/tskit-dev/msprime/issues/517 */

    tsk_table_collection_free(&tables);
}

static void
test_reference_sequence_state_machine(void)
{

    tsk_reference_sequence_t r1;

    tsk_reference_sequence_init(&r1, 0);
    CU_ASSERT_EQUAL(r1.data, NULL);
    CU_ASSERT_EQUAL(r1.url, NULL);
    CU_ASSERT_EQUAL(r1.metadata, NULL);
    CU_ASSERT_EQUAL(r1.metadata_schema, NULL);
    CU_ASSERT_TRUE(tsk_reference_sequence_is_null(&r1));

    CU_ASSERT_EQUAL(tsk_reference_sequence_set_data(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    /* Setting the value back to NULL makes the reference whole object NULL */
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_data(&r1, NULL, 0), 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_is_null(&r1));
    tsk_reference_sequence_free(&r1);
    CU_ASSERT_TRUE(tsk_reference_sequence_is_null(&r1));

    /* Any empty string is the same thing. */
    tsk_reference_sequence_init(&r1, 0);
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_data(&r1, "", 0), 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_is_null(&r1));
    tsk_reference_sequence_free(&r1);

    tsk_reference_sequence_init(&r1, 0);
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_url(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    tsk_reference_sequence_free(&r1);

    tsk_reference_sequence_init(&r1, 0);
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_metadata(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    tsk_reference_sequence_free(&r1);

    tsk_reference_sequence_init(&r1, 0);
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_metadata_schema(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    tsk_reference_sequence_free(&r1);

    tsk_reference_sequence_init(&r1, 0);
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_metadata(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_metadata_schema(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_url(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_data(&r1, "x", 1), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));

    CU_ASSERT_EQUAL(tsk_reference_sequence_set_metadata(&r1, "", 0), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_metadata_schema(&r1, "", 0), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_url(&r1, "", 0), 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_is_null(&r1));
    CU_ASSERT_EQUAL(tsk_reference_sequence_set_data(&r1, "", 0), 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_is_null(&r1));

    tsk_reference_sequence_free(&r1);
}

static void
test_reference_sequence_take(void)
{
    int ret;
    tsk_reference_sequence_t r1;
    tsk_reference_sequence_t r2;
    const char *const_data = "data";
    const char *const_metadata = "metadata";
    char *takeset_data = strdup(const_data);
    char *takeset_metadata = strdup(const_metadata);

    ret = tsk_reference_sequence_init(&r1, 0);

    ret = tsk_reference_sequence_set_data(&r1, const_data, strlen(const_data));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_set_metadata(
        &r1, const_metadata, strlen(const_metadata));
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_reference_sequence_init(&r2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_equals(&r1, &r2, 0));
    ret = tsk_reference_sequence_takeset_data(&r2, takeset_data, strlen(takeset_data));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_equals(&r1, &r2, 0));
    ret = tsk_reference_sequence_takeset_metadata(
        &r2, takeset_metadata, strlen(takeset_metadata));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    /* Writing over these with copies doesn't lose memory */
    ret = tsk_reference_sequence_set_data(&r2, const_data, strlen(const_data));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_set_metadata(
        &r2, const_metadata, strlen(const_metadata));
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    /* The original copies are gone, make some new ones */
    takeset_data = strdup(const_data);
    takeset_metadata = strdup(const_metadata);

    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_takeset_data(&r1, takeset_data, strlen(takeset_data));
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_takeset_metadata(
        &r1, takeset_metadata, strlen(takeset_metadata));
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    tsk_reference_sequence_free(&r1);
    tsk_reference_sequence_free(&r2);
}

static void
test_reference_sequence(void)
{
    int ret;
    tsk_reference_sequence_t r1;
    tsk_reference_sequence_t r2;

    const char example_data[100] = "An example string with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_data_length = (tsk_size_t) strlen(example_data);
    const char example_url[100] = "An example url with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_url_length = (tsk_size_t) strlen(example_url);
    const char example_metadata[100] = "An example metadata with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_metadata_length = (tsk_size_t) strlen(example_metadata);
    const char example_schema[100] = "An example schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_schema_length = (tsk_size_t) strlen(example_schema);

    tsk_reference_sequence_init(&r1, 0);
    tsk_reference_sequence_init(&r2, 0);

    /* NULL sequences are initially equal */
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_data(&r1, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_data(&r1, "", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_data(&r2, "", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_data(&r1, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_data(&r2, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_url(&r1, example_url, example_url_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_equals(&r1, &r2, 0));
    ret = tsk_reference_sequence_set_url(&r2, example_url, example_url_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_metadata(
        &r1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_equals(&r1, &r2, 0));
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, TSK_CMP_IGNORE_METADATA));
    ret = tsk_reference_sequence_set_metadata(
        &r2, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, TSK_CMP_IGNORE_METADATA));

    ret = tsk_reference_sequence_set_metadata_schema(
        &r1, example_schema, example_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_reference_sequence_equals(&r1, &r2, 0));
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, TSK_CMP_IGNORE_METADATA));
    ret = tsk_reference_sequence_set_metadata_schema(
        &r2, example_schema, example_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, TSK_CMP_IGNORE_METADATA));

    // Test copy
    tsk_reference_sequence_free(&r1);
    tsk_reference_sequence_free(&r2);

    tsk_reference_sequence_init(&r1, 0);
    ret = tsk_reference_sequence_set_data(&r1, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_copy(&r1, &r2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_url(&r1, example_url, example_url_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_copy(&r1, &r2, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_metadata(
        &r1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_copy(&r1, &r2, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    ret = tsk_reference_sequence_set_metadata_schema(
        &r1, example_schema, example_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_copy(&r1, &r2, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_reference_sequence_equals(&r1, &r2, 0));

    tsk_reference_sequence_free(&r1);
    tsk_reference_sequence_free(&r2);
}

static void
test_table_collection_reference_sequence(void)
{
    int ret;
    tsk_table_collection_t tc1, tc2;

    char example_data[100] = "An example string with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_data_length = (tsk_size_t) strlen(example_data);
    char example_url[100] = "An example url with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_url_length = (tsk_size_t) strlen(example_url);
    char example_metadata[100] = "An example metadata with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_metadata_length = (tsk_size_t) strlen(example_metadata);
    char example_schema[100] = "An example schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_schema_length = (tsk_size_t) strlen(example_schema);

    // Test equality
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_data(
        &tc1.reference_sequence, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_data(
        &tc2.reference_sequence, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_url(
        &tc1.reference_sequence, example_url, example_url_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_reference_sequence_set_url(
        &tc2.reference_sequence, example_url, example_url_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_metadata(
        &tc1.reference_sequence, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_reference_sequence_set_metadata(
        &tc2.reference_sequence, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_metadata_schema(
        &tc1.reference_sequence, example_schema, example_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_reference_sequence_set_metadata_schema(
        &tc2.reference_sequence, example_schema, example_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    // Test copy
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_reference_sequence_set_data(
        &tc1.reference_sequence, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_copy(&tc1, &tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_url(
        &tc1.reference_sequence, example_url, example_url_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_copy(&tc1, &tc2, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_metadata(
        &tc1.reference_sequence, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_copy(&tc1, &tc2, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_reference_sequence_set_metadata_schema(
        &tc1.reference_sequence, example_schema, example_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_copy(&tc1, &tc2, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);

    // Test dump and load
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tc1.sequence_length = 1.0;
    ret = tsk_reference_sequence_set_data(
        &tc1.reference_sequence, example_data, example_data_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_set_url(
        &tc1.reference_sequence, example_url, example_url_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_set_metadata(
        &tc1.reference_sequence, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_reference_sequence_set_metadata_schema(
        &tc1.reference_sequence, example_schema, example_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_dump(&tc1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tc2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
}

static void
test_table_collection_has_reference_sequence(void)
{
    int ret;
    tsk_table_collection_t tc;

    ret = tsk_table_collection_init(&tc, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tc.sequence_length = 1.0;

    CU_ASSERT_FALSE(tsk_table_collection_has_reference_sequence(&tc));
    ret = tsk_reference_sequence_set_data(&tc.reference_sequence, "A", 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_reference_sequence(&tc));
    /* Goes back to NULL by setting a empty string. See
     * test_reference_sequence_state_machine for detailed tests. */
    ret = tsk_reference_sequence_set_data(&tc.reference_sequence, "", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_reference_sequence(&tc));

    tsk_table_collection_free(&tc);
}

static void
test_table_collection_metadata(void)
{
    int ret;
    tsk_table_collection_t tc1, tc2;

    char example_metadata[100] = "An example of metadata with unicode 🎄🌳🌴🌲🎋";
    char *takeset_metadata;
    char example_metadata_schema[100]
        = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_metadata_length = (tsk_size_t) strlen(example_metadata);
    tsk_size_t example_metadata_schema_length
        = (tsk_size_t) strlen(example_metadata_schema);

    // Test equality
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_table_collection_set_metadata(
        &tc1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_table_collection_set_metadata(
        &tc2, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_table_collection_set_metadata_schema(
        &tc1, example_metadata_schema, example_metadata_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_table_collection_set_metadata_schema(
        &tc2, example_metadata_schema, example_metadata_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    // Test copy
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata(
        &tc1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_copy(&tc1, &tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    ret = tsk_table_collection_set_metadata_schema(
        &tc1, example_metadata_schema, example_metadata_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_copy(&tc1, &tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    // Test dump and load with empty metadata and schema
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tc1.sequence_length = 1.0;
    ret = tsk_table_collection_dump(&tc1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tc2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    // Test dump and load with set metadata and schema
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tc1.sequence_length = 1.0;
    ret = tsk_table_collection_set_metadata(
        &tc1, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata_schema(
        &tc1, example_metadata_schema, example_metadata_schema_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_dump(&tc1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tc2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);

    takeset_metadata = tsk_malloc(example_metadata_length * sizeof(char));
    CU_ASSERT_FATAL(takeset_metadata != NULL);
    memcpy(takeset_metadata, &example_metadata,
        (size_t)(example_metadata_length * sizeof(char)));

    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_takeset_metadata(
        &tc1, takeset_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(tc1.metadata, &example_metadata, example_metadata_length), 0);
    tsk_table_collection_free(&tc1);
}

static void
test_table_collection_time_units(void)
{
    int ret;
    tsk_table_collection_t tc1, tc2;

    char example_time_units[100] = "An example of time units with unicode ⏰";
    tsk_size_t example_time_units_length = (tsk_size_t) strlen(example_time_units);

    // Test equality
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_init(&tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_table_collection_set_time_units(
        &tc1, example_time_units, example_time_units_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&tc1, &tc2, 0));
    ret = tsk_table_collection_set_time_units(
        &tc2, example_time_units, example_time_units_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    // Test copy
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_time_units(
        &tc1, example_time_units, example_time_units_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_copy(&tc1, &tc2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    // Test dump and load with default time_units
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(ret, strncmp(tc1.time_units, TSK_TIME_UNITS_UNKNOWN, 7));
    tc1.sequence_length = 1.0;
    ret = tsk_table_collection_dump(&tc1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tc2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));

    // Test dump and load with set time_units and schema
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
    ret = tsk_table_collection_init(&tc1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tc1.sequence_length = 1.0;
    ret = tsk_table_collection_set_time_units(
        &tc1, example_time_units, example_time_units_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_dump(&tc1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tc2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tc1, &tc2, 0));
    tsk_table_collection_free(&tc1);
    tsk_table_collection_free(&tc2);
}

static void
test_node_table(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_node_table_t table, table2;
    tsk_node_t node, node2;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    tsk_flags_t *flags;
    tsk_id_t *population;
    double *time;
    tsk_id_t *individual;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    char metadata_copy[test_metadata_length + 1];
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    metadata_copy[test_metadata_length] = '\0';
    ret = tsk_node_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_node_table_set_max_rows_increment(&table, 1);
    tsk_node_table_set_max_metadata_length_increment(&table, 1);
    tsk_node_table_print_state(&table, _devnull);
    ret = tsk_node_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_node_table_add_row(&table, (tsk_flags_t) j, (double) j, j, j,
            test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
        CU_ASSERT_EQUAL(table.flags[j], (tsk_flags_t) j);
        CU_ASSERT_EQUAL(table.time[j], j);
        CU_ASSERT_EQUAL(table.population[j], j);
        CU_ASSERT_EQUAL(table.individual[j], j);
        CU_ASSERT_EQUAL(table.num_rows, (tsk_size_t) j + 1);
        CU_ASSERT_EQUAL(
            table.metadata_length, (tsk_size_t)(j + 1) * test_metadata_length);
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], table.metadata_length);
        /* check the metadata */
        tsk_memcpy(metadata_copy, table.metadata + table.metadata_offset[j],
            test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(metadata_copy, test_metadata, test_metadata_length);
        ret = tsk_node_table_get_row(&table, (tsk_id_t) j, &node);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(node.id, j);
        CU_ASSERT_EQUAL(node.flags, (tsk_size_t) j);
        CU_ASSERT_EQUAL(node.time, j);
        CU_ASSERT_EQUAL(node.population, j);
        CU_ASSERT_EQUAL(node.individual, j);
        CU_ASSERT_EQUAL(node.metadata_length, test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(node.metadata, test_metadata, test_metadata_length);
    }

    /* Test equality with and without metadata */
    tsk_node_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the metadata values */
    table2.metadata[0] = 0;
    CU_ASSERT_FALSE(tsk_node_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the last metadata entry */
    table2.metadata_offset[table2.num_rows]
        = table2.metadata_offset[table2.num_rows - 1];
    CU_ASSERT_FALSE(tsk_node_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Delete all metadata */
    tsk_memset(table2.metadata_offset, 0,
        (table2.num_rows + 1) * sizeof(*table2.metadata_offset));
    CU_ASSERT_FALSE(tsk_node_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    tsk_node_table_free(&table2);

    CU_ASSERT_EQUAL(tsk_node_table_get_row(&table, (tsk_id_t) num_rows, &node),
        TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_node_table_print_state(&table, _devnull);
    ret = tsk_node_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_node_table_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    num_rows *= 2;
    flags = tsk_malloc(num_rows * sizeof(tsk_flags_t));
    CU_ASSERT_FATAL(flags != NULL);
    tsk_memset(flags, 1, num_rows * sizeof(tsk_flags_t));
    population = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(population != NULL);
    tsk_memset(population, 2, num_rows * sizeof(tsk_id_t));
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    tsk_memset(time, 0, num_rows * sizeof(double));
    individual = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(individual != NULL);
    tsk_memset(individual, 3, num_rows * sizeof(tsk_id_t));
    metadata = tsk_malloc(num_rows * sizeof(char));
    tsk_memset(metadata, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        metadata_offset[j] = (tsk_size_t) j;
    }
    ret = tsk_node_table_set_columns(&table, num_rows, flags, time, population,
        individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.population, population, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.individual, individual, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    tsk_node_table_print_state(&table, _devnull);
    ret = tsk_node_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Append another num_rows onto the end */
    ret = tsk_node_table_append_columns(&table, num_rows, flags, time, population,
        individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.flags + num_rows, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.population, population, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.population + num_rows, population, num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.individual, individual, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.individual + num_rows, individual, num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    tsk_node_table_print_state(&table, _devnull);
    ret = tsk_node_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Truncate back to the original number of rows. */
    ret = tsk_node_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.population, population, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.individual, individual, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    ret = tsk_node_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* If population is NULL it should be set to -1. If metadata is NULL all metadatas
     * should be set to the empty string. If individual is NULL it should be set to -1.
     */
    num_rows = 10;
    tsk_memset(population, 0xff, num_rows * sizeof(tsk_id_t));
    tsk_memset(individual, 0xff, num_rows * sizeof(tsk_id_t));
    ret = tsk_node_table_set_columns(
        &table, num_rows, flags, time, NULL, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.population, population, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.individual, individual, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        num_rows * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* flags and time cannot be NULL */
    ret = tsk_node_table_set_columns(
        &table, num_rows, NULL, time, population, individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_table_set_columns(&table, num_rows, flags, NULL, population,
        individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_table_set_columns(
        &table, num_rows, flags, time, population, individual, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_table_set_columns(
        &table, num_rows, flags, time, population, individual, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* if metadata and metadata_offset are both null, all metadatas are zero length */
    num_rows = 10;
    tsk_memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_size_t));
    ret = tsk_node_table_set_columns(
        &table, num_rows, flags, time, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    ret = tsk_node_table_append_columns(
        &table, num_rows, flags, time, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.flags + num_rows, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset + num_rows, metadata_offset,
                        num_rows * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    tsk_node_table_print_state(&table, _devnull);
    ret = tsk_node_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Test extend method */
    ret = tsk_node_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_init(&table2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_node_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, 0));
    ret = tsk_node_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_node_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_node_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    ret = tsk_node_table_set_columns(&table2, num_rows, flags, time, population,
        individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_node_table_equals(&table, &table2, 0));
    ret = tsk_node_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_node_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_node_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, 0));
    ret = tsk_node_table_set_columns(&table2, num_rows, flags, time, population,
        individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Copy a subset */
    ret = tsk_node_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_node_table_equals(&table, &table2, 0));
    ret = tsk_node_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_row_subset; j++) {
        ret = tsk_node_table_get_row(&table, j, &node);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_node_table_get_row(&table2, row_subset[j], &node2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(node.flags, node2.flags);
        CU_ASSERT_EQUAL(node.time, node2.time);
        CU_ASSERT_EQUAL(node.population, node2.population);
        CU_ASSERT_EQUAL(node.individual, node2.individual);
        CU_ASSERT_EQUAL(node.metadata_length, node2.metadata_length);
        CU_ASSERT_EQUAL(tsk_memcmp(node.metadata, node2.metadata,
                            node.metadata_length * sizeof(*node.metadata)),
            0);
    }

    ret = tsk_node_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, 0);
    CU_ASSERT_EQUAL(table.metadata_schema, NULL);
    const char *example = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example);
    const char *example2 = "A different example 🎄🌳🌴🌲🎋";
    tsk_size_t example2_length = (tsk_size_t) strlen(example);
    tsk_node_table_set_metadata_schema(&table, example, example_length);
    CU_ASSERT_EQUAL(table.metadata_schema_length, example_length);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_schema, example, example_length), 0);

    tsk_node_table_copy(&table, &table2, TSK_NO_INIT);
    CU_ASSERT_EQUAL(table.metadata_schema_length, table2.metadata_schema_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_schema, table2.metadata_schema, example_length), 0);
    tsk_node_table_set_metadata_schema(&table2, example, example_length);
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, 0));
    tsk_node_table_set_metadata_schema(&table2, example2, example2_length);
    CU_ASSERT_FALSE(tsk_node_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_node_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));

    tsk_node_table_clear(&table);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_node_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_node_table_free(&table2);
    CU_ASSERT_EQUAL(ret, 0);
    free(flags);
    free(population);
    free(time);
    free(metadata);
    free(metadata_offset);
    free(individual);
}

static void
test_node_table_takeset(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_node_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    tsk_flags_t *flags;
    double *time;
    tsk_id_t *population;
    tsk_id_t *individual;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    tsk_size_t zeros[num_rows + 1];
    tsk_id_t neg_ones[num_rows];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    tsk_memset(neg_ones, 0xff, num_rows * sizeof(tsk_id_t));
    /* Make a table to copy from */
    ret = tsk_node_table_init(&source_table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_node_table_add_row(&source_table, (tsk_flags_t) j, (double) j + 1,
            j + 2, j + 3, test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    flags = tsk_malloc(num_rows * sizeof(tsk_flags_t));
    CU_ASSERT_FATAL(flags != NULL);
    tsk_memcpy(flags, source_table.flags, num_rows * sizeof(tsk_flags_t));
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    tsk_memcpy(time, source_table.time, num_rows * sizeof(double));
    population = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(population != NULL);
    tsk_memcpy(population, source_table.population, num_rows * sizeof(tsk_id_t));
    individual = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(individual != NULL);
    tsk_memcpy(individual, source_table.individual, num_rows * sizeof(tsk_id_t));
    metadata = tsk_malloc(num_rows * test_metadata_length * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    tsk_memcpy(
        metadata, source_table.metadata, num_rows * test_metadata_length * sizeof(char));
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    tsk_memcpy(metadata_offset, source_table.metadata_offset,
        (num_rows + 1) * sizeof(tsk_size_t));

    ret = tsk_node_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_node_table_add_row(
        &table, (tsk_flags_t) 1, 2, 3, 4, test_metadata, test_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_node_table_takeset_columns(&table, num_rows, flags, time, population,
        individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_node_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_node_table_takeset_columns(
        &table, num_rows, NULL, time, population, individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_table_takeset_columns(&table, num_rows, flags, NULL, population,
        individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_table_takeset_columns(
        &table, num_rows, flags, time, population, individual, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_table_takeset_columns(
        &table, num_rows, flags, time, population, individual, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_node_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    flags = tsk_malloc(num_rows * sizeof(tsk_flags_t));
    CU_ASSERT_FATAL(flags != NULL);
    tsk_memcpy(flags, source_table.flags, num_rows * sizeof(tsk_flags_t));
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    tsk_memcpy(time, source_table.time, num_rows * sizeof(double));
    /* if metadata and offset are both null, all entries are zero length,
       individual and population default to -1 */
    num_rows = 10;
    ret = tsk_node_table_takeset_columns(
        &table, num_rows, flags, time, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.population, neg_ones, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.individual, neg_ones, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    ret = tsk_node_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_node_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_node_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_node_table_t table;
    tsk_node_t row;
    const char *metadata = "ABC";

    ret = tsk_node_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_node_table_add_row(&table, 0, 1.0, 2, 3, metadata, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(&table, 1, 2.0, 3, 4, metadata, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(&table, 2, 3.0, 4, 5, metadata, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_node_table_update_row(&table, 0, 1, 2.0, 3, 4, &metadata[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 1);
    CU_ASSERT_EQUAL_FATAL(row.time, 2.0);
    CU_ASSERT_EQUAL_FATAL(row.population, 3);
    CU_ASSERT_EQUAL_FATAL(row.individual, 4);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_node_table_update_row(&table, 0, row.flags + 1, row.time + 1,
        row.population + 1, row.individual + 1, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 2);
    CU_ASSERT_EQUAL_FATAL(row.time, 3.0);
    CU_ASSERT_EQUAL_FATAL(row.population, 4);
    CU_ASSERT_EQUAL_FATAL(row.individual, 5);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_node_table_update_row(&table, 0, 0, 0, 0, 0, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 0);
    CU_ASSERT_EQUAL_FATAL(row.time, 0);
    CU_ASSERT_EQUAL_FATAL(row.population, 0);
    CU_ASSERT_EQUAL_FATAL(row.individual, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_node_table_update_row(&table, 1, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 0);
    CU_ASSERT_EQUAL_FATAL(row.time, 0);
    CU_ASSERT_EQUAL_FATAL(row.population, 0);
    CU_ASSERT_EQUAL_FATAL(row.individual, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_node_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 2);
    CU_ASSERT_EQUAL_FATAL(row.time, 3.0);
    CU_ASSERT_EQUAL_FATAL(row.population, 4);
    CU_ASSERT_EQUAL_FATAL(row.individual, 5);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_node_table_update_row(&table, 3, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    tsk_node_table_free(&table);
}

static void
test_edge_table_with_options(tsk_flags_t options)
{
    int ret;
    tsk_edge_table_t table, table2;
    tsk_size_t num_rows = 100;
    tsk_id_t j, ret_id;
    tsk_edge_t edge, edge2;
    tsk_id_t *parent, *child;
    double *left, *right;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    char metadata_copy[test_metadata_length + 1];
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    metadata_copy[test_metadata_length] = '\0';
    ret = tsk_edge_table_init(&table, options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_edge_table_set_max_rows_increment(&table, 1);
    tsk_edge_table_set_max_metadata_length_increment(&table, 1);
    tsk_edge_table_print_state(&table, _devnull);
    ret = tsk_edge_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        if (options & TSK_TABLE_NO_METADATA) {
            ret_id = tsk_edge_table_add_row(&table, (double) j, (double) j, j, j,
                test_metadata, test_metadata_length);
            CU_ASSERT_EQUAL(ret_id, TSK_ERR_METADATA_DISABLED);
            ret_id
                = tsk_edge_table_add_row(&table, (double) j, (double) j, j, j, NULL, 0);
        } else {
            ret_id = tsk_edge_table_add_row(&table, (double) j, (double) j, j, j,
                test_metadata, test_metadata_length);
        }
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
        CU_ASSERT_EQUAL(table.left[j], j);
        CU_ASSERT_EQUAL(table.right[j], j);
        CU_ASSERT_EQUAL(table.parent[j], j);
        CU_ASSERT_EQUAL(table.child[j], j);
        CU_ASSERT_EQUAL(table.num_rows, (tsk_size_t) j + 1);
        if (options & TSK_TABLE_NO_METADATA) {
            CU_ASSERT_EQUAL(table.metadata_length, 0);
            CU_ASSERT_EQUAL(table.metadata, NULL);
            CU_ASSERT_EQUAL(table.metadata_offset, NULL);
        } else {
            CU_ASSERT_EQUAL(
                table.metadata_length, (tsk_size_t)(j + 1) * test_metadata_length);
            CU_ASSERT_EQUAL(table.metadata_offset[j + 1], table.metadata_length);
            /* check the metadata */
            tsk_memcpy(metadata_copy, table.metadata + table.metadata_offset[j],
                test_metadata_length);
            CU_ASSERT_NSTRING_EQUAL(metadata_copy, test_metadata, test_metadata_length);
        }

        ret = tsk_edge_table_get_row(&table, (tsk_id_t) j, &edge);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(edge.id, j);
        CU_ASSERT_EQUAL(edge.left, j);
        CU_ASSERT_EQUAL(edge.right, j);
        CU_ASSERT_EQUAL(edge.parent, j);
        CU_ASSERT_EQUAL(edge.child, j);
        if (options & TSK_TABLE_NO_METADATA) {
            CU_ASSERT_EQUAL(edge.metadata_length, 0);
            CU_ASSERT_EQUAL(edge.metadata, NULL);
        } else {
            CU_ASSERT_EQUAL(edge.metadata_length, test_metadata_length);
            CU_ASSERT_NSTRING_EQUAL(edge.metadata, test_metadata, test_metadata_length);
        }
    }
    ret = tsk_edge_table_get_row(&table, (tsk_id_t) num_rows, &edge);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    tsk_edge_table_print_state(&table, _devnull);
    ret = tsk_edge_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    num_rows *= 2;
    left = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    tsk_memset(left, 0, num_rows * sizeof(double));
    right = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    tsk_memset(right, 0, num_rows * sizeof(double));
    parent = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parent != NULL);
    tsk_memset(parent, 1, num_rows * sizeof(tsk_id_t));
    child = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(child != NULL);
    tsk_memset(child, 1, num_rows * sizeof(tsk_id_t));
    metadata = tsk_malloc(num_rows * sizeof(char));
    tsk_memset(metadata, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        metadata_offset[j] = (tsk_size_t) j;
    }
    if (options & TSK_TABLE_NO_METADATA) {
        ret = tsk_edge_table_set_columns(
            &table, num_rows, left, right, parent, child, metadata, metadata_offset);
        CU_ASSERT_EQUAL(ret, TSK_ERR_METADATA_DISABLED);
        ret = tsk_edge_table_set_columns(
            &table, num_rows, left, right, parent, child, NULL, NULL);
    } else {
        ret = tsk_edge_table_set_columns(
            &table, num_rows, left, right, parent, child, metadata, metadata_offset);
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    if (options & TSK_TABLE_NO_METADATA) {
        CU_ASSERT_EQUAL(table.metadata, NULL);
        CU_ASSERT_EQUAL(table.metadata_offset, NULL);
        CU_ASSERT_EQUAL(table.metadata_length, 0);
    } else {
        CU_ASSERT_EQUAL(
            tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
        CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                            (num_rows + 1) * sizeof(tsk_size_t)),
            0);
        CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    }

    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    /* Append another num_rows to the end. */
    if (options & TSK_TABLE_NO_METADATA) {
        ret = tsk_edge_table_append_columns(
            &table, num_rows, left, right, parent, child, metadata, metadata_offset);
        CU_ASSERT_EQUAL(ret, TSK_ERR_METADATA_DISABLED);
        ret = tsk_edge_table_append_columns(
            &table, num_rows, left, right, parent, child, NULL, NULL);
    } else {
        ret = tsk_edge_table_append_columns(
            &table, num_rows, left, right, parent, child, metadata, metadata_offset);
    }
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.left + num_rows, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.right + num_rows, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parent + num_rows, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.child + num_rows, child, num_rows * sizeof(tsk_id_t)), 0);
    if (options & TSK_TABLE_NO_METADATA) {
        CU_ASSERT_EQUAL(table.metadata, NULL);
        CU_ASSERT_EQUAL(table.metadata_offset, NULL);
        CU_ASSERT_EQUAL(table.metadata_length, 0);
    } else {
        CU_ASSERT_EQUAL(
            tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
        CU_ASSERT_EQUAL(
            tsk_memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
        CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    }

    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_edge_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    if (options & TSK_TABLE_NO_METADATA) {
        CU_ASSERT_EQUAL(table.metadata, NULL);
        CU_ASSERT_EQUAL(table.metadata_offset, NULL);
        CU_ASSERT_EQUAL(table.metadata_length, 0);
    } else {
        CU_ASSERT_EQUAL(
            tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
        CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                            (num_rows + 1) * sizeof(tsk_size_t)),
            0);
        CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    }
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    ret = tsk_edge_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Test equality with and without metadata */
    tsk_edge_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    if (!(options & TSK_TABLE_NO_METADATA)) {
        /* Change the metadata values */
        table2.metadata[0] = 0;
        CU_ASSERT_FALSE(tsk_edge_table_equals(&table, &table2, 0));
        CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
        /* Change the last metadata entry */
        table2.metadata_offset[table2.num_rows]
            = table2.metadata_offset[table2.num_rows - 1];
        CU_ASSERT_FALSE(tsk_edge_table_equals(&table, &table2, 0));
        CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
        /* Delete all metadata */
        tsk_memset(table2.metadata_offset, 0,
            (table2.num_rows + 1) * sizeof(*table2.metadata_offset));
        CU_ASSERT_FALSE(tsk_edge_table_equals(&table, &table2, 0));
        CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    }
    tsk_edge_table_free(&table2);

    /* Inputs cannot be NULL */
    ret = tsk_edge_table_set_columns(
        &table, num_rows, NULL, right, parent, child, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_set_columns(
        &table, num_rows, left, NULL, parent, child, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_set_columns(
        &table, num_rows, left, right, NULL, child, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_set_columns(
        &table, num_rows, left, right, parent, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_set_columns(
        &table, num_rows, left, right, parent, child, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_set_columns(
        &table, num_rows, left, right, parent, child, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* if metadata and metadata_offset are both null, all metadatas are zero length */
    num_rows = 10;
    tsk_memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_size_t));
    ret = tsk_edge_table_set_columns(
        &table, num_rows, left, right, parent, child, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    if (options & TSK_TABLE_NO_METADATA) {
        CU_ASSERT_EQUAL(table.metadata, NULL);
        CU_ASSERT_EQUAL(table.metadata_offset, NULL);
    } else {
        CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                            (num_rows + 1) * sizeof(tsk_size_t)),
            0);
    }
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    ret = tsk_edge_table_append_columns(
        &table, num_rows, left, right, parent, child, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.left + num_rows, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.right + num_rows, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parent + num_rows, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.child + num_rows, child, num_rows * sizeof(tsk_id_t)), 0);
    if (options & TSK_TABLE_NO_METADATA) {
        CU_ASSERT_EQUAL(table.metadata, NULL);
        CU_ASSERT_EQUAL(table.metadata_offset, NULL);
    } else {
        CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                            (num_rows + 1) * sizeof(tsk_size_t)),
            0);
        CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset + num_rows, metadata_offset,
                            num_rows * sizeof(tsk_size_t)),
            0);
    }
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    tsk_edge_table_print_state(&table, _devnull);
    ret = tsk_edge_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Test extend method */
    ret = tsk_edge_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_init(&table2, options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_edge_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, 0));
    ret = tsk_edge_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_edge_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_edge_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    if (options & TSK_TABLE_NO_METADATA) {
        ret = tsk_edge_table_set_columns(
            &table2, num_rows, left, right, parent, child, NULL, NULL);
    } else {
        ret = tsk_edge_table_set_columns(
            &table2, num_rows, left, right, parent, child, metadata, metadata_offset);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_edge_table_equals(&table, &table2, 0));
    ret = tsk_edge_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_edge_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_edge_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, 0));
    if (options & TSK_TABLE_NO_METADATA) {
        ret = tsk_edge_table_set_columns(
            &table2, num_rows, left, right, parent, child, NULL, NULL);
    } else {
        ret = tsk_edge_table_set_columns(
            &table2, num_rows, left, right, parent, child, metadata, metadata_offset);
    }
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Copy a subset */
    ret = tsk_edge_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_edge_table_equals(&table, &table2, 0));
    ret = tsk_edge_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_row_subset; j++) {
        ret = tsk_edge_table_get_row(&table, j, &edge);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_edge_table_get_row(&table2, row_subset[j], &edge2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(edge.parent, edge2.parent);
        CU_ASSERT_EQUAL(edge.child, edge2.child);
        CU_ASSERT_EQUAL(edge.left, edge2.left);
        CU_ASSERT_EQUAL(edge.right, edge2.right);
        CU_ASSERT_EQUAL(edge.metadata_length, edge2.metadata_length)
        CU_ASSERT_EQUAL(tsk_memcmp(edge.metadata, edge2.metadata,
                            edge.metadata_length * sizeof(*edge.metadata)),
            0);
    }

    ret = tsk_edge_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, 0);
    CU_ASSERT_EQUAL(table.metadata_schema, NULL);
    const char *example = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example);
    const char *example2 = "A different example 🎄🌳🌴🌲🎋";
    tsk_size_t example2_length = (tsk_size_t) strlen(example);
    ret = tsk_edge_table_set_metadata_schema(&table, example, example_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, example_length);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_schema, example, example_length), 0);

    ret = tsk_edge_table_copy(&table, &table2, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, table2.metadata_schema_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_schema, table2.metadata_schema, example_length), 0);
    ret = tsk_edge_table_set_metadata_schema(&table2, example, example_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, 0));
    ret = tsk_edge_table_set_metadata_schema(&table2, example2, example2_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_edge_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_edge_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));

    ret = tsk_edge_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    ret = tsk_edge_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_edge_table_free(&table2);
    CU_ASSERT_EQUAL(ret, 0);
    free(left);
    free(right);
    free(parent);
    free(child);
    free(metadata);
    free(metadata_offset);
}

static void
test_edge_table(void)
{
    test_edge_table_with_options(0);
    test_edge_table_with_options(TSK_TABLE_NO_METADATA);
}

static void
test_edge_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_edge_table_t table;
    tsk_edge_t row;
    const char *metadata = "ABC";

    ret = tsk_edge_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_edge_table_add_row(&table, 0, 1.0, 2, 3, metadata, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&table, 1, 2.0, 3, 4, metadata, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&table, 2, 3.0, 4, 5, metadata, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_edge_table_update_row(&table, 0, 1, 2.0, 3, 4, &metadata[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 1);
    CU_ASSERT_EQUAL_FATAL(row.right, 2.0);
    CU_ASSERT_EQUAL_FATAL(row.parent, 3);
    CU_ASSERT_EQUAL_FATAL(row.child, 4);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_edge_table_update_row(&table, 0, row.left + 1, row.right + 1,
        row.parent + 1, row.child + 1, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 2);
    CU_ASSERT_EQUAL_FATAL(row.right, 3.0);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.child, 5);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_edge_table_update_row(&table, 0, 0, 0, 0, 0, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 0);
    CU_ASSERT_EQUAL_FATAL(row.right, 0);
    CU_ASSERT_EQUAL_FATAL(row.parent, 0);
    CU_ASSERT_EQUAL_FATAL(row.child, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_edge_table_update_row(&table, 1, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 0);
    CU_ASSERT_EQUAL_FATAL(row.right, 0);
    CU_ASSERT_EQUAL_FATAL(row.parent, 0);
    CU_ASSERT_EQUAL_FATAL(row.child, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_edge_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 2);
    CU_ASSERT_EQUAL_FATAL(row.right, 3.0);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.child, 5);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_edge_table_update_row(&table, 3, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);

    tsk_edge_table_free(&table);
}

static void
test_edge_table_update_row_no_metadata(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_edge_table_t table;
    tsk_edge_t row;
    const char *metadata = "ABC";

    ret = tsk_edge_table_init(&table, TSK_TABLE_NO_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_edge_table_add_row(&table, 0, 1.0, 2, 3, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&table, 1, 2.0, 3, 4, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&table, 2, 3.0, 4, 5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_edge_table_update_row(&table, 0, 1, 2.0, 3, 4, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 1);
    CU_ASSERT_EQUAL_FATAL(row.right, 2.0);
    CU_ASSERT_EQUAL_FATAL(row.parent, 3);
    CU_ASSERT_EQUAL_FATAL(row.child, 4);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_edge_table_update_row(&table, 0, row.left + 1, row.right + 1,
        row.parent + 1, row.child + 1, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 2);
    CU_ASSERT_EQUAL_FATAL(row.right, 3.0);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.child, 5);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_edge_table_update_row(&table, 1, 0, 0, 0, 0, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_METADATA_DISABLED);

    tsk_edge_table_free(&table);
}

static void
test_edge_table_takeset_with_options(tsk_flags_t table_options)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_edge_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    double *left;
    double *right;
    tsk_id_t *parent;
    tsk_id_t *child;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    tsk_size_t zeros[num_rows + 1];
    tsk_id_t neg_ones[num_rows];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    tsk_memset(neg_ones, 0xff, num_rows * sizeof(tsk_id_t));
    /* Make a table to copy from */
    ret = tsk_edge_table_init(&source_table, table_options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        if (table_options & TSK_TABLE_NO_METADATA) {
            ret_id = tsk_edge_table_add_row(
                &source_table, (double) j, (double) j + 1, j + 2, j + 3, NULL, 0);

        } else {
            ret_id = tsk_edge_table_add_row(&source_table, (double) j, (double) j + 1,
                j + 2, j + 3, test_metadata, test_metadata_length);
        }
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    left = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    tsk_memcpy(left, source_table.left, num_rows * sizeof(double));
    right = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    tsk_memcpy(right, source_table.right, num_rows * sizeof(double));
    parent = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parent != NULL);
    tsk_memcpy(parent, source_table.parent, num_rows * sizeof(tsk_id_t));
    child = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(child != NULL);
    tsk_memcpy(child, source_table.child, num_rows * sizeof(tsk_id_t));
    if (table_options & TSK_TABLE_NO_METADATA) {
        metadata = NULL;
        metadata_offset = NULL;
        test_metadata = NULL;
        test_metadata_length = 0;
    } else {
        metadata = tsk_malloc(num_rows * test_metadata_length * sizeof(char));
        CU_ASSERT_FATAL(metadata != NULL);
        tsk_memcpy(metadata, source_table.metadata,
            num_rows * test_metadata_length * sizeof(char));
        metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
        CU_ASSERT_FATAL(metadata_offset != NULL);
        tsk_memcpy(metadata_offset, source_table.metadata_offset,
            (num_rows + 1) * sizeof(tsk_size_t));
    }

    ret = tsk_edge_table_init(&table, table_options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_edge_table_add_row(
        &table, 1, 2, 3, 4, test_metadata, test_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_edge_table_takeset_columns(
        &table, num_rows, left, right, parent, child, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_edge_table_takeset_columns(
        &table, num_rows, NULL, right, parent, child, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_takeset_columns(
        &table, num_rows, left, NULL, parent, child, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_takeset_columns(
        &table, num_rows, left, right, NULL, child, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_table_takeset_columns(
        &table, num_rows, left, right, parent, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    if (table_options & TSK_TABLE_NO_METADATA) {
        /* It isn't used, so any pointer does for testing that presence of metadata
            fails */
        ret = tsk_edge_table_takeset_columns(
            &table, num_rows, left, right, parent, child, (char *) child, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_METADATA_DISABLED);
    } else {
        ret = tsk_edge_table_takeset_columns(
            &table, num_rows, left, right, parent, child, NULL, metadata_offset);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
        ret = tsk_edge_table_takeset_columns(
            &table, num_rows, left, right, parent, child, metadata, NULL);
        CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    }

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_edge_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    ret = tsk_edge_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_edge_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_edge_table_takeset(void)
{
    test_edge_table_takeset_with_options(TSK_TABLE_NO_METADATA);
    test_edge_table_takeset_with_options(0);
}

static void
test_edge_table_copy_semantics(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t t1, t2;
    tsk_edge_table_t edges;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_edge_metadata(&t1);

    /* t1 now has metadata. We should be able to copy to another table with metadata */
    ret = tsk_table_collection_copy(&t1, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);

    /* We should not be able to copy into a table with no metadata */
    ret = tsk_table_collection_copy(&t1, &t2, TSK_TC_NO_EDGE_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_METADATA_DISABLED);
    tsk_table_collection_free(&t2);

    tsk_table_collection_free(&t1);
    ret = tsk_treeseq_copy_tables(&ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* t1 has no metadata, but metadata is enabled. We should be able to copy
     * into a table with either metadata enabled or disabled.
     */
    ret = tsk_table_collection_copy(&t1, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);

    ret = tsk_table_collection_copy(&t1, &t2, TSK_TC_NO_EDGE_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);

    /* Try copying into a table directly */
    ret = tsk_edge_table_copy(&t1.edges, &edges, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&t1.edges, &edges, 0));
    tsk_edge_table_free(&edges);

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(&ts);
}

static void
test_edge_table_squash(void)
{
    int ret;
    tsk_table_collection_t tables;

    const char *nodes_ex = "1  0       -1   -1\n"
                           "1  0       -1   -1\n"
                           "0  0.253   -1   -1\n";
    const char *edges_ex = "0  2   2   0\n"
                           "2  10  2   0\n"
                           "0  2   2   1\n"
                           "2  10  2   1\n";

    /*
      2
     / \
    0   1
    */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 10;

    parse_nodes(nodes_ex, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 3);
    parse_edges(edges_ex, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 4);

    ret = tsk_edge_table_squash(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // Check output.
    CU_ASSERT_EQUAL(tables.edges.num_rows, 2);

    // Free things.
    tsk_table_collection_free(&tables);
}

static void
test_edge_table_squash_multiple_parents(void)
{
    int ret;
    tsk_table_collection_t tables;

    const char *nodes_ex = "1  0.000   -1    -1\n"
                           "1  0.000   -1    -1\n"
                           "1  0.000   -1    -1\n"
                           "1  0.000   -1    -1\n"
                           "0  1.000   -1    -1\n"
                           "0  1.000   -1    -1\n";
    const char *edges_ex = "5  10  5   3\n"
                           "5  10  5   2\n"
                           "0  5   5   3\n"
                           "0  5   5   2\n"
                           "4  10  4   1\n"
                           "0  4   4   1\n"
                           "4  10  4   0\n"
                           "0  4   4   0\n";
    /*
                4       5
               / \     / \
              0   1   2   3
    */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 10;

    parse_nodes(nodes_ex, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 6);
    parse_edges(edges_ex, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 8);

    ret = tsk_edge_table_squash(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // Check output.
    CU_ASSERT_EQUAL(tables.edges.num_rows, 4);

    // Free things.
    tsk_table_collection_free(&tables);
}

static void
test_edge_table_squash_empty(void)
{
    int ret;
    tsk_table_collection_t tables;

    const char *nodes_ex = "1  0       -1   -1\n"
                           "1  0       -1   -1\n"
                           "0  0.253   -1   -1\n";
    const char *edges_ex = "";

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 10;

    parse_nodes(nodes_ex, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 3);
    parse_edges(edges_ex, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 0);

    ret = tsk_edge_table_squash(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // Free things.
    tsk_table_collection_free(&tables);
}

static void
test_edge_table_squash_single_edge(void)
{
    int ret;
    tsk_table_collection_t tables;

    const char *nodes_ex = "1  0   -1   -1\n"
                           "0  0   -1   -1\n";
    const char *edges_ex = "0  1   1   0\n";
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    parse_nodes(nodes_ex, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 2);
    parse_edges(edges_ex, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 1);

    ret = tsk_edge_table_squash(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // Free things.
    tsk_table_collection_free(&tables);
}

static void
test_edge_table_squash_bad_intervals(void)
{
    int ret;
    tsk_table_collection_t tables;

    const char *nodes_ex = "1  0   -1   -1\n"
                           "0  0   -1   -1\n";
    const char *edges_ex = "0  0.6   1   0\n"
                           "0.4  1   1   0\n";

    ret = tsk_table_collection_init(&tables, TSK_TC_NO_EDGE_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    parse_nodes(nodes_ex, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 2);
    parse_edges(edges_ex, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 2);

    ret = tsk_edge_table_squash(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_EDGES_CONTRADICTORY_CHILDREN);

    // Free things.
    tsk_table_collection_free(&tables);
}

static void
test_edge_table_squash_metadata(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 10;
    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 0, 1, 1, "metadata", 8);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_edge_table_squash(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANT_PROCESS_EDGES_WITH_METADATA);

    tsk_table_collection_free(&tables);

    ret = tsk_table_collection_init(&tables, TSK_TC_NO_EDGE_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 10;
    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 0, 1, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_edge_table_squash(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_table_collection_free(&tables);
}

static void
test_site_table(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_site_table_t table, table2;
    tsk_size_t num_rows, j;
    char *ancestral_state;
    char *metadata;
    double *position;
    tsk_site_t site, site2;
    tsk_size_t *ancestral_state_offset;
    tsk_size_t *metadata_offset;
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    ret = tsk_site_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_site_table_set_max_rows_increment(&table, 1);
    tsk_site_table_set_max_metadata_length_increment(&table, 1);
    tsk_site_table_set_max_ancestral_state_length_increment(&table, 1);
    tsk_site_table_print_state(&table, _devnull);
    ret = tsk_site_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_site_table_add_row(&table, 0, "A", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    CU_ASSERT_EQUAL(table.position[0], 0);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[0], 0);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 1);
    CU_ASSERT_EQUAL(table.metadata_offset[0], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[1], 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, 1);

    ret = tsk_site_table_get_row(&table, 0, &site);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(site.position, 0);
    CU_ASSERT_EQUAL(site.ancestral_state_length, 1);
    CU_ASSERT_NSTRING_EQUAL(site.ancestral_state, "A", 1);
    CU_ASSERT_EQUAL(site.metadata_length, 0);

    ret_id = tsk_site_table_add_row(&table, 1, "AA", 2, "{}", 2);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    CU_ASSERT_EQUAL(table.position[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[2], 3);
    CU_ASSERT_EQUAL(table.metadata_offset[1], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[2], 2);
    CU_ASSERT_EQUAL(table.metadata_length, 2);
    CU_ASSERT_EQUAL(table.num_rows, 2);

    ret = tsk_site_table_get_row(&table, 1, &site);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(site.position, 1);
    CU_ASSERT_EQUAL(site.ancestral_state_length, 2);
    CU_ASSERT_NSTRING_EQUAL(site.ancestral_state, "AA", 2);
    CU_ASSERT_EQUAL(site.metadata_length, 2);
    CU_ASSERT_NSTRING_EQUAL(site.metadata, "{}", 2);

    ret_id = tsk_site_table_add_row(&table, 2, "A", 1, "metadata", 8);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    CU_ASSERT_EQUAL(table.position[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[3], 4);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 4);
    CU_ASSERT_EQUAL(table.metadata_offset[3], 10);
    CU_ASSERT_EQUAL(table.metadata_length, 10);
    CU_ASSERT_EQUAL(table.num_rows, 3);

    ret = tsk_site_table_get_row(&table, 3, &site);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);

    tsk_site_table_print_state(&table, _devnull);
    ret = tsk_site_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_site_table_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[0], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[0], 0);

    num_rows = 100;
    position = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(position != NULL);
    ancestral_state = tsk_malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(ancestral_state != NULL);
    ancestral_state_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(ancestral_state_offset != NULL);
    metadata = tsk_malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);

    for (j = 0; j < num_rows; j++) {
        position[j] = (double) j;
        ancestral_state[j] = (char) j;
        ancestral_state_offset[j] = (tsk_size_t) j;
        metadata[j] = (char) ('A' + j);
        metadata_offset[j] = (tsk_size_t) j;
    }
    ancestral_state_offset[num_rows] = num_rows;
    metadata_offset[num_rows] = num_rows;

    ret = tsk_site_table_set_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.position, position, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.ancestral_state, ancestral_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, num_rows);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    /* Append another num rows */
    ret = tsk_site_table_append_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.position, position, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.position + num_rows, position, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.ancestral_state, ancestral_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.ancestral_state + num_rows, ancestral_state,
                        num_rows * sizeof(char)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 2 * num_rows);

    /* truncate back to num_rows */
    ret = tsk_site_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.position, position, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.ancestral_state, ancestral_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, num_rows);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    ret = tsk_site_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Test equality with and without metadata */
    tsk_site_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the metadata values */
    table2.metadata[0] = 0;
    CU_ASSERT_FALSE(tsk_site_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the last metadata entry */
    table2.metadata_offset[table2.num_rows]
        = table2.metadata_offset[table2.num_rows - 1];
    CU_ASSERT_FALSE(tsk_site_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Delete all metadata */
    tsk_memset(table2.metadata_offset, 0,
        (table2.num_rows + 1) * sizeof(*table2.metadata_offset));
    CU_ASSERT_FALSE(tsk_site_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    tsk_site_table_free(&table2);

    /* Inputs cannot be NULL */
    ret = tsk_site_table_set_columns(&table, num_rows, NULL, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_table_set_columns(&table, num_rows, position, NULL,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_table_set_columns(
        &table, num_rows, position, ancestral_state, NULL, metadata, metadata_offset);
    /* Metadata and metadata_offset must both be null */
    ret = tsk_site_table_set_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_table_set_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Set metadata to NULL */
    ret = tsk_site_table_set_columns(
        &table, num_rows, position, ancestral_state, ancestral_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_EQUAL(tsk_memcmp(table.position, position, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.ancestral_state, ancestral_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, num_rows);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    /* Test extend method */
    ret = tsk_site_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_init(&table2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_site_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, 0));
    ret = tsk_site_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_site_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_site_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    ret = tsk_site_table_set_columns(&table2, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_site_table_equals(&table, &table2, 0));
    ret = tsk_site_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_site_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_site_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, 0));
    ret = tsk_site_table_set_columns(&table2, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Copy a subset */
    ret = tsk_site_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_site_table_equals(&table, &table2, 0));
    ret = tsk_site_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_row_subset; j++) {
        ret = tsk_site_table_get_row(&table, (tsk_id_t) j, &site);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_site_table_get_row(&table2, row_subset[j], &site2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(site.position, site2.position);
        CU_ASSERT_EQUAL(site.ancestral_state_length, site2.ancestral_state_length);
        CU_ASSERT_EQUAL(site.metadata_length, site2.metadata_length);
        CU_ASSERT_EQUAL(tsk_memcmp(site.ancestral_state, site2.ancestral_state,
                            site.ancestral_state_length * sizeof(*site.ancestral_state)),
            0);
        CU_ASSERT_EQUAL(tsk_memcmp(site.metadata, site2.metadata,
                            site.metadata_length * sizeof(*site.metadata)),
            0);
    }

    /* Test for bad offsets */
    ancestral_state_offset[0] = 1;
    ret = tsk_site_table_set_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    ancestral_state_offset[0] = 0;
    ancestral_state_offset[num_rows] = 0;
    ret = tsk_site_table_set_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    ancestral_state_offset[0] = 0;

    metadata_offset[0] = 0;
    ret = tsk_site_table_set_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    metadata_offset[0] = 0;
    metadata_offset[num_rows] = 0;
    ret = tsk_site_table_set_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    ret = tsk_site_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, 0);
    CU_ASSERT_EQUAL(table.metadata_schema, NULL);
    const char *example = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example);
    const char *example2 = "A different example 🎄🌳🌴🌲🎋";
    tsk_size_t example2_length = (tsk_size_t) strlen(example);
    tsk_site_table_set_metadata_schema(&table, example, example_length);
    CU_ASSERT_EQUAL(table.metadata_schema_length, example_length);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_schema, example, example_length), 0);

    tsk_site_table_copy(&table, &table2, TSK_NO_INIT);
    CU_ASSERT_EQUAL(table.metadata_schema_length, table2.metadata_schema_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_schema, table2.metadata_schema, example_length), 0);
    tsk_site_table_set_metadata_schema(&table2, example, example_length);
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, 0));
    tsk_site_table_set_metadata_schema(&table2, example2, example2_length);
    CU_ASSERT_FALSE(tsk_site_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_site_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));

    ret = tsk_site_table_clear(&table);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_site_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_site_table_free(&table2);
    CU_ASSERT_EQUAL(ret, 0);

    free(position);
    free(ancestral_state);
    free(ancestral_state_offset);
    free(metadata);
    free(metadata_offset);
}

static void
test_site_table_takeset(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_site_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    double *position;
    char *ancestral_state;
    tsk_size_t *ancestral_state_offset;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_ancestral_state = "red";
    tsk_size_t test_ancestral_state_length = 3;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    tsk_size_t zeros[num_rows + 1];
    tsk_id_t neg_ones[num_rows];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    tsk_memset(neg_ones, 0xff, num_rows * sizeof(tsk_id_t));
    /* Make a table to copy from */
    ret = tsk_site_table_init(&source_table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_site_table_add_row(&source_table, (double) j, test_ancestral_state,
            test_ancestral_state_length, test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    position = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(position != NULL);
    tsk_memcpy(position, source_table.position, num_rows * sizeof(double));
    ancestral_state = tsk_malloc(num_rows * test_ancestral_state_length * sizeof(char));
    CU_ASSERT_FATAL(ancestral_state != NULL);
    tsk_memcpy(ancestral_state, source_table.ancestral_state,
        num_rows * test_ancestral_state_length * sizeof(char));
    ancestral_state_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(ancestral_state_offset != NULL);
    tsk_memcpy(ancestral_state_offset, source_table.ancestral_state_offset,
        (num_rows + 1) * sizeof(tsk_size_t));
    metadata = tsk_malloc(num_rows * test_metadata_length * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    tsk_memcpy(
        metadata, source_table.metadata, num_rows * test_metadata_length * sizeof(char));
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    tsk_memcpy(metadata_offset, source_table.metadata_offset,
        (num_rows + 1) * sizeof(tsk_size_t));

    ret = tsk_site_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_site_table_add_row(&table, 1, test_ancestral_state,
        test_ancestral_state_length, test_metadata, test_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_site_table_takeset_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_site_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_site_table_takeset_columns(&table, num_rows, NULL, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_table_takeset_columns(&table, num_rows, position, NULL,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_table_takeset_columns(
        &table, num_rows, position, ancestral_state, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_table_takeset_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_table_takeset_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* Check bad offset in ancestral_state */
    ancestral_state_offset[0] = 1;
    ret = tsk_site_table_takeset_columns(&table, num_rows, position, ancestral_state,
        ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_site_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    position = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(position != NULL);
    tsk_memcpy(position, source_table.position, num_rows * sizeof(double));
    ancestral_state = tsk_malloc(num_rows * test_ancestral_state_length * sizeof(char));
    CU_ASSERT_FATAL(ancestral_state != NULL);
    tsk_memcpy(ancestral_state, source_table.ancestral_state,
        num_rows * test_ancestral_state_length * sizeof(char));
    ancestral_state_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(ancestral_state_offset != NULL);
    tsk_memcpy(ancestral_state_offset, source_table.ancestral_state_offset,
        (num_rows + 1) * sizeof(tsk_size_t));
    /* if metadata and offset are both null, all entries are zero length*/
    num_rows = 10;
    ret = tsk_site_table_takeset_columns(
        &table, num_rows, position, ancestral_state, ancestral_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    ret = tsk_site_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_site_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_site_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_site_table_t table;
    tsk_site_t row;
    const char *ancestral_state = "XYZ";
    const char *metadata = "ABC";

    ret = tsk_site_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_site_table_add_row(&table, 0, ancestral_state, 1, metadata, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&table, 1, ancestral_state, 2, metadata, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&table, 2, ancestral_state, 3, metadata, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_site_table_update_row(
        &table, 0, 1, &ancestral_state[1], 1, &metadata[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.position, 1);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_site_table_update_row(&table, 0, row.position + 1, row.ancestral_state,
        row.ancestral_state_length, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.position, 2);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_site_table_update_row(&table, 0, row.position, row.ancestral_state,
        row.ancestral_state_length, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.position, 2);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_site_table_update_row(
        &table, 0, row.position, NULL, 0, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.position, 2);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_site_table_update_row(&table, 0, 2, ancestral_state, 3, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.position, 2);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[0], 'X');
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[1], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[2], 'Z');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_site_table_update_row(&table, 1, 5, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.position, 5);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_site_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.position, 2);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[0], 'X');
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[1], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.ancestral_state[2], 'Z');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_site_table_update_row(&table, 3, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);

    tsk_site_table_free(&table);
}

static void
test_mutation_table(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_mutation_table_t table, table2;
    tsk_size_t num_rows = 100;
    tsk_size_t max_len = 20;
    tsk_size_t k, len;
    tsk_id_t j;
    tsk_id_t *node;
    tsk_id_t *parent;
    tsk_id_t *site;
    double *time;
    char *derived_state, *metadata;
    char c[max_len + 1];
    tsk_size_t *derived_state_offset, *metadata_offset;
    tsk_mutation_t mutation, mutation2;
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    for (j = 0; j < (tsk_id_t) max_len; j++) {
        c[j] = (char) ('A' + j);
    }

    ret = tsk_mutation_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_mutation_table_set_max_rows_increment(&table, 1);
    tsk_mutation_table_set_max_metadata_length_increment(&table, 1);
    tsk_mutation_table_set_max_derived_state_length_increment(&table, 1);
    tsk_mutation_table_print_state(&table, _devnull);
    ret = tsk_mutation_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    len = 0;
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        k = TSK_MIN((tsk_size_t) j + 1, max_len);
        ret_id = tsk_mutation_table_add_row(&table, j, j, j, (double) j, c, k, c, k);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
        CU_ASSERT_EQUAL(table.site[j], j);
        CU_ASSERT_EQUAL(table.node[j], j);
        CU_ASSERT_EQUAL(table.parent[j], j);
        CU_ASSERT_EQUAL(table.time[j], j);
        CU_ASSERT_EQUAL(table.derived_state_offset[j], len);
        CU_ASSERT_EQUAL(table.metadata_offset[j], len);
        CU_ASSERT_EQUAL(table.num_rows, (tsk_size_t) j + 1);
        len += k;
        CU_ASSERT_EQUAL(table.derived_state_offset[j + 1], len);
        CU_ASSERT_EQUAL(table.derived_state_length, len);
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], len);
        CU_ASSERT_EQUAL(table.metadata_length, len);

        ret = tsk_mutation_table_get_row(&table, (tsk_id_t) j, &mutation);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(mutation.id, j);
        CU_ASSERT_EQUAL(mutation.site, j);
        CU_ASSERT_EQUAL(mutation.node, j);
        CU_ASSERT_EQUAL(mutation.parent, j);
        CU_ASSERT_EQUAL(mutation.time, j);
        CU_ASSERT_EQUAL(mutation.metadata_length, k);
        CU_ASSERT_NSTRING_EQUAL(mutation.metadata, c, k);
        CU_ASSERT_EQUAL(mutation.derived_state_length, k);
        CU_ASSERT_NSTRING_EQUAL(mutation.derived_state, c, k);
    }
    ret = tsk_mutation_table_get_row(&table, (tsk_id_t) num_rows, &mutation);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    tsk_mutation_table_print_state(&table, _devnull);
    ret = tsk_mutation_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    num_rows *= 2;
    site = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(site != NULL);
    node = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    parent = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parent != NULL);
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    derived_state = tsk_malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(derived_state != NULL);
    derived_state_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(derived_state_offset != NULL);
    metadata = tsk_malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        node[j] = j;
        site[j] = j + 1;
        parent[j] = j + 2;
        time[j] = (double) (j + 3);
        derived_state[j] = 'Y';
        derived_state_offset[j] = (tsk_size_t) j;
        metadata[j] = 'M';
        metadata_offset[j] = (tsk_size_t) j;
    }

    derived_state_offset[num_rows] = num_rows;
    metadata_offset[num_rows] = num_rows;
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.derived_state, derived_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* Append another num_rows */
    ret = tsk_mutation_table_append_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.site + num_rows, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.node + num_rows, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parent + num_rows, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.derived_state, derived_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.derived_state, derived_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.derived_state_length, 2 * num_rows);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_mutation_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.derived_state, derived_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* Test equality with and without metadata */
    tsk_mutation_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the metadata values */
    table2.metadata[0] = 0;
    CU_ASSERT_FALSE(tsk_mutation_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the last metadata entry */
    table2.metadata_offset[table2.num_rows]
        = table2.metadata_offset[table2.num_rows - 1];
    CU_ASSERT_FALSE(tsk_mutation_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Delete all metadata */
    tsk_memset(table2.metadata_offset, 0,
        (table2.num_rows + 1) * sizeof(*table2.metadata_offset));
    CU_ASSERT_FALSE(tsk_mutation_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    tsk_mutation_table_free(&table2);

    ret = tsk_mutation_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Check all this again, except with parent == NULL, time == NULL
     * and metadata == NULL. */
    tsk_memset(parent, 0xff, num_rows * sizeof(tsk_id_t));
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        time[j] = TSK_UNKNOWN_TIME;
    }
    tsk_memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_size_t));
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, NULL, NULL,
        derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.derived_state, derived_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.derived_state_offset, derived_state_offset,
                        num_rows * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    /* Append another num_rows */
    ret = tsk_mutation_table_append_columns(&table, num_rows, site, node, NULL, NULL,
        derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.site + num_rows, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.node + num_rows, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parent + num_rows, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.derived_state, derived_state, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.derived_state + num_rows, derived_state,
                        num_rows * sizeof(char)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    /* Inputs except parent, time, metadata and metadata_offset cannot be NULL*/
    ret = tsk_mutation_table_set_columns(&table, num_rows, NULL, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, NULL, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, parent, time,
        NULL, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, parent, time,
        derived_state, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Inputs except parent, time, metadata and metadata_offset cannot be NULL*/
    ret = tsk_mutation_table_append_columns(&table, num_rows, NULL, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_append_columns(&table, num_rows, site, NULL, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_append_columns(&table, num_rows, site, node, parent, time,
        NULL, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_append_columns(&table, num_rows, site, node, parent, time,
        derived_state, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_append_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_append_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Test extend method */
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        parent[j] = j + 2;
        time[j] = (double) (j + 3);
        metadata[j] = (char) ('A' + j);
        metadata_offset[j] = (tsk_size_t) j;
    }
    metadata_offset[num_rows] = num_rows;
    ret = tsk_mutation_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_init(&table2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_mutation_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, 0));
    ret = tsk_mutation_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_mutation_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_mutation_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    ret = tsk_mutation_table_set_columns(&table2, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_mutation_table_equals(&table, &table2, 0));
    ret = tsk_mutation_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_mutation_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_mutation_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, 0));
    ret = tsk_mutation_table_set_columns(&table2, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Copy a subset */
    ret = tsk_mutation_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_mutation_table_equals(&table, &table2, 0));
    ret = tsk_mutation_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (k = 0; k < num_row_subset; k++) {
        ret = tsk_mutation_table_get_row(&table, (tsk_id_t) k, &mutation);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_mutation_table_get_row(&table2, row_subset[k], &mutation2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(mutation.site, mutation2.site);
        CU_ASSERT_EQUAL(mutation.node, mutation2.node);
        CU_ASSERT_EQUAL(mutation.parent, mutation2.parent);
        CU_ASSERT_EQUAL(mutation.time, mutation2.time);
        CU_ASSERT_EQUAL(mutation.derived_state_length, mutation2.derived_state_length);
        CU_ASSERT_EQUAL(mutation.metadata_length, mutation2.metadata_length);
        CU_ASSERT_EQUAL(
            tsk_memcmp(mutation.derived_state, mutation2.derived_state,
                mutation.derived_state_length * sizeof(*mutation.derived_state)),
            0);
        CU_ASSERT_EQUAL(tsk_memcmp(mutation.metadata, mutation2.metadata,
                            mutation.metadata_length * sizeof(*mutation.metadata)),
            0);
    }

    /* Test for bad offsets */
    derived_state_offset[0] = 1;
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    derived_state_offset[0] = 0;
    derived_state_offset[num_rows] = 0;
    ret = tsk_mutation_table_set_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    ret = tsk_mutation_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, 0);
    CU_ASSERT_EQUAL(table.metadata_schema, NULL);
    const char *example = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example);
    const char *example2 = "A different example 🎄🌳🌴🌲🎋";
    tsk_size_t example2_length = (tsk_size_t) strlen(example);
    tsk_mutation_table_set_metadata_schema(&table, example, example_length);
    CU_ASSERT_EQUAL(table.metadata_schema_length, example_length);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_schema, example, example_length), 0);

    tsk_mutation_table_copy(&table, &table2, TSK_NO_INIT);
    CU_ASSERT_EQUAL(table.metadata_schema_length, table2.metadata_schema_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_schema, table2.metadata_schema, example_length), 0);
    tsk_mutation_table_set_metadata_schema(&table2, example, example_length);
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, 0));
    tsk_mutation_table_set_metadata_schema(&table2, example2, example2_length);
    CU_ASSERT_FALSE(tsk_mutation_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));

    tsk_mutation_table_clear(&table);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.derived_state_length, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_mutation_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_mutation_table_free(&table2);
    CU_ASSERT_EQUAL(ret, 0);
    free(site);
    free(node);
    free(parent);
    free(time);
    free(derived_state);
    free(derived_state_offset);
    free(metadata);
    free(metadata_offset);
}

static void
test_mutation_table_takeset(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_mutation_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    tsk_id_t *site;
    tsk_id_t *node;
    tsk_id_t *parent;
    double *time;
    char *derived_state;
    tsk_size_t *derived_state_offset;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_derived_state = "red";
    tsk_size_t test_derived_state_length = 3;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    tsk_size_t zeros[num_rows + 1];
    tsk_id_t neg_ones[num_rows];
    double unknown_times[num_rows];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    tsk_memset(neg_ones, 0xff, num_rows * sizeof(tsk_id_t));
    /* Make a table to copy from */
    ret = tsk_mutation_table_init(&source_table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        unknown_times[j] = TSK_UNKNOWN_TIME;
        ret_id = tsk_mutation_table_add_row(&source_table, j, j + 1, j + 2,
            (double) j + 3, test_derived_state, test_derived_state_length, test_metadata,
            test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    site = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(site != NULL);
    tsk_memcpy(site, source_table.site, num_rows * sizeof(tsk_id_t));
    node = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    tsk_memcpy(node, source_table.node, num_rows * sizeof(tsk_id_t));
    parent = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parent != NULL);
    tsk_memcpy(parent, source_table.parent, num_rows * sizeof(tsk_id_t));
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    tsk_memcpy(time, source_table.time, num_rows * sizeof(double));
    derived_state = tsk_malloc(num_rows * test_derived_state_length * sizeof(char));
    CU_ASSERT_FATAL(derived_state != NULL);
    tsk_memcpy(derived_state, source_table.derived_state,
        num_rows * test_derived_state_length * sizeof(char));
    derived_state_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(derived_state_offset != NULL);
    tsk_memcpy(derived_state_offset, source_table.derived_state_offset,
        (num_rows + 1) * sizeof(tsk_size_t));
    metadata = tsk_malloc(num_rows * test_metadata_length * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    tsk_memcpy(
        metadata, source_table.metadata, num_rows * test_metadata_length * sizeof(char));
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    tsk_memcpy(metadata_offset, source_table.metadata_offset,
        (num_rows + 1) * sizeof(tsk_size_t));

    ret = tsk_mutation_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_mutation_table_add_row(&table, 1, 1, 1, 1, test_derived_state,
        test_derived_state_length, test_metadata, test_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_mutation_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, NULL, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, NULL, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    // Parent and time not tested as they have deafults
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, node, parent, time,
        NULL, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, node, parent, time,
        derived_state, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Check error on bad derived_state offset */
    derived_state_offset[0] = 1;
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, node, parent, time,
        derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_mutation_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    // Re init non-optional arrays
    site = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(site != NULL);
    tsk_memcpy(site, source_table.site, num_rows * sizeof(tsk_id_t));
    node = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    tsk_memcpy(node, source_table.node, num_rows * sizeof(tsk_id_t));
    derived_state = tsk_malloc(num_rows * test_derived_state_length * sizeof(char));
    CU_ASSERT_FATAL(derived_state != NULL);
    tsk_memcpy(derived_state, source_table.derived_state,
        num_rows * test_derived_state_length * sizeof(char));
    derived_state_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(derived_state_offset != NULL);
    tsk_memcpy(derived_state_offset, source_table.derived_state_offset,
        (num_rows + 1) * sizeof(tsk_size_t));
    /* if metadata and offset are both null, all entries are zero length, if parent or
     * time are NULL they default to null values*/
    num_rows = 10;
    ret = tsk_mutation_table_takeset_columns(&table, num_rows, site, node, NULL, NULL,
        derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parent, neg_ones, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.time, unknown_times, num_rows * sizeof(tsk_id_t)), 0);

    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    ret = tsk_mutation_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_mutation_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_mutation_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_mutation_table_t table;
    tsk_mutation_t row;
    const char *derived_state = "XYZ";
    const char *metadata = "ABC";

    ret = tsk_mutation_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id
        = tsk_mutation_table_add_row(&table, 0, 1, 2, 3, derived_state, 1, metadata, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&table, 1, 2, 3, 4, derived_state, 2, metadata, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&table, 2, 3, 4, 5, derived_state, 3, metadata, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_mutation_table_update_row(
        &table, 0, 1, 2, 3, 4, &derived_state[1], 1, &metadata[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.site, 1);
    CU_ASSERT_EQUAL_FATAL(row.node, 2);
    CU_ASSERT_EQUAL_FATAL(row.parent, 3);
    CU_ASSERT_EQUAL_FATAL(row.time, 4);
    CU_ASSERT_EQUAL_FATAL(row.derived_state_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.derived_state[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_mutation_table_update_row(&table, 0, row.site + 1, row.node + 1,
        row.parent + 1, row.time + 1, row.derived_state, row.derived_state_length,
        row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.site, 2);
    CU_ASSERT_EQUAL_FATAL(row.node, 3);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.time, 5);
    CU_ASSERT_EQUAL_FATAL(row.derived_state_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.derived_state[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_mutation_table_update_row(&table, 0, row.site, row.node, row.parent,
        row.time, row.derived_state, row.derived_state_length, row.metadata,
        row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.site, 2);
    CU_ASSERT_EQUAL_FATAL(row.node, 3);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.time, 5);
    CU_ASSERT_EQUAL_FATAL(row.derived_state_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.derived_state[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_mutation_table_update_row(&table, 0, row.site, row.node, row.parent,
        row.time, NULL, 0, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.site, 2);
    CU_ASSERT_EQUAL_FATAL(row.node, 3);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.time, 5);
    CU_ASSERT_EQUAL_FATAL(row.derived_state_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_mutation_table_update_row(
        &table, 0, 2, 3, 4, 5, derived_state, 3, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.site, 2);
    CU_ASSERT_EQUAL_FATAL(row.node, 3);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.time, 5);
    CU_ASSERT_EQUAL_FATAL(row.derived_state_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.derived_state[0], 'X');
    CU_ASSERT_EQUAL_FATAL(row.derived_state[1], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.derived_state[2], 'Z');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_mutation_table_update_row(&table, 1, 5, 6, 7, 8, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.site, 5);
    CU_ASSERT_EQUAL_FATAL(row.node, 6);
    CU_ASSERT_EQUAL_FATAL(row.parent, 7);
    CU_ASSERT_EQUAL_FATAL(row.time, 8);
    CU_ASSERT_EQUAL_FATAL(row.derived_state_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_mutation_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.site, 2);
    CU_ASSERT_EQUAL_FATAL(row.node, 3);
    CU_ASSERT_EQUAL_FATAL(row.parent, 4);
    CU_ASSERT_EQUAL_FATAL(row.time, 5);
    CU_ASSERT_EQUAL_FATAL(row.derived_state_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.derived_state[0], 'X');
    CU_ASSERT_EQUAL_FATAL(row.derived_state[1], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.derived_state[2], 'Z');
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_mutation_table_update_row(&table, 3, 0, 0, 0, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);

    tsk_mutation_table_free(&table);
}

static void
test_migration_table(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_migration_table_t table, table2;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    tsk_id_t *node;
    tsk_id_t *source, *dest;
    double *left, *right, *time;
    tsk_migration_t migration, migration2;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    char metadata_copy[test_metadata_length + 1];
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    metadata_copy[test_metadata_length] = '\0';
    ret = tsk_migration_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_migration_table_set_max_rows_increment(&table, 1);
    tsk_migration_table_print_state(&table, _devnull);
    ret = tsk_migration_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_migration_table_add_row(&table, (double) j, (double) j, j, j, j,
            (double) j, test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
        CU_ASSERT_EQUAL(table.left[j], j);
        CU_ASSERT_EQUAL(table.right[j], j);
        CU_ASSERT_EQUAL(table.node[j], j);
        CU_ASSERT_EQUAL(table.source[j], j);
        CU_ASSERT_EQUAL(table.dest[j], j);
        CU_ASSERT_EQUAL(table.time[j], j);
        CU_ASSERT_EQUAL(table.num_rows, (tsk_size_t) j + 1);
        CU_ASSERT_EQUAL(
            table.metadata_length, (tsk_size_t)(j + 1) * test_metadata_length);
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], table.metadata_length);
        /* check the metadata */
        tsk_memcpy(metadata_copy, table.metadata + table.metadata_offset[j],
            test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(metadata_copy, test_metadata, test_metadata_length);

        ret = tsk_migration_table_get_row(&table, (tsk_id_t) j, &migration);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(migration.id, j);
        CU_ASSERT_EQUAL(migration.left, j);
        CU_ASSERT_EQUAL(migration.right, j);
        CU_ASSERT_EQUAL(migration.node, j);
        CU_ASSERT_EQUAL(migration.source, j);
        CU_ASSERT_EQUAL(migration.dest, j);
        CU_ASSERT_EQUAL(migration.time, j);
        CU_ASSERT_EQUAL(migration.metadata_length, test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(migration.metadata, test_metadata, test_metadata_length);
    }
    ret = tsk_migration_table_get_row(&table, (tsk_id_t) num_rows, &migration);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    tsk_migration_table_print_state(&table, _devnull);
    ret = tsk_migration_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    num_rows *= 2;
    left = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    tsk_memset(left, 1, num_rows * sizeof(double));
    right = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    tsk_memset(right, 2, num_rows * sizeof(double));
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    tsk_memset(time, 3, num_rows * sizeof(double));
    node = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    tsk_memset(node, 4, num_rows * sizeof(tsk_id_t));
    source = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(source != NULL);
    tsk_memset(source, 5, num_rows * sizeof(tsk_id_t));
    dest = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(dest != NULL);
    tsk_memset(dest, 6, num_rows * sizeof(tsk_id_t));
    metadata = tsk_malloc(num_rows * sizeof(char));
    tsk_memset(metadata, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        metadata_offset[j] = (tsk_size_t) j;
    }

    ret = tsk_migration_table_set_columns(&table, num_rows, left, right, node, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* Append another num_rows */
    ret = tsk_migration_table_append_columns(&table, num_rows, left, right, node, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.left + num_rows, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.right + num_rows, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.node + num_rows, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.source + num_rows, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.dest + num_rows, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_migration_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* Test equality with and without metadata */
    tsk_migration_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the metadata values */
    table2.metadata[0] = 0;
    CU_ASSERT_FALSE(tsk_migration_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the last metadata entry */
    table2.metadata_offset[table2.num_rows]
        = table2.metadata_offset[table2.num_rows - 1];
    CU_ASSERT_FALSE(tsk_migration_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Delete all metadata */
    tsk_memset(table2.metadata_offset, 0,
        (table2.num_rows + 1) * sizeof(*table2.metadata_offset));
    CU_ASSERT_FALSE(tsk_migration_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    tsk_migration_table_free(&table2);

    ret = tsk_migration_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* inputs cannot be NULL */
    ret = tsk_migration_table_set_columns(&table, num_rows, NULL, right, node, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_set_columns(&table, num_rows, left, NULL, node, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_set_columns(&table, num_rows, left, right, NULL, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_set_columns(&table, num_rows, left, right, node, NULL,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_set_columns(&table, num_rows, left, right, node, source,
        NULL, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_set_columns(&table, num_rows, left, right, node, source,
        dest, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_set_columns(
        &table, num_rows, left, right, node, source, dest, time, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_set_columns(
        &table, num_rows, left, right, node, source, dest, time, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    tsk_migration_table_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);

    /* if metadata and metadata_offset are both null, all metadatas are zero length */
    num_rows = 10;
    tsk_memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_size_t));
    ret = tsk_migration_table_set_columns(
        &table, num_rows, left, right, node, source, dest, time, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    ret = tsk_migration_table_append_columns(
        &table, num_rows, left, right, node, source, dest, time, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.left + num_rows, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.right + num_rows, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.node + num_rows, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.source + num_rows, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.dest + num_rows, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset + num_rows, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    tsk_migration_table_print_state(&table, _devnull);
    ret = tsk_migration_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Test extend method */
    ret = tsk_migration_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_init(&table2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_migration_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, 0));
    ret = tsk_migration_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_migration_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_migration_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    ret = tsk_migration_table_set_columns(&table2, num_rows, left, right, node, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_migration_table_equals(&table, &table2, 0));
    ret = tsk_migration_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_migration_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_migration_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, 0));
    ret = tsk_migration_table_set_columns(&table2, num_rows, left, right, node, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Copy a subset */
    ret = tsk_migration_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_migration_table_equals(&table, &table2, 0));
    ret = tsk_migration_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_row_subset; j++) {
        ret = tsk_migration_table_get_row(&table, j, &migration);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_migration_table_get_row(&table2, row_subset[j], &migration2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(migration.source, migration2.source);
        CU_ASSERT_EQUAL(migration.dest, migration2.dest);
        CU_ASSERT_EQUAL(migration.node, migration2.node);
        CU_ASSERT_EQUAL(migration.left, migration2.left);
        CU_ASSERT_EQUAL(migration.right, migration2.right);
        CU_ASSERT_EQUAL(migration.time, migration2.time);
        CU_ASSERT_EQUAL(migration.metadata_length, migration2.metadata_length);
        CU_ASSERT_EQUAL(tsk_memcmp(migration.metadata, migration2.metadata,
                            migration.metadata_length * sizeof(*migration.metadata)),
            0);
    }

    ret = tsk_migration_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, 0);
    CU_ASSERT_EQUAL(table.metadata_schema, NULL);
    const char *example = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example);
    const char *example2 = "A different example 🎄🌳🌴🌲🎋";
    tsk_size_t example2_length = (tsk_size_t) strlen(example);
    tsk_migration_table_set_metadata_schema(&table, example, example_length);
    CU_ASSERT_EQUAL(table.metadata_schema_length, example_length);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_schema, example, example_length), 0);

    tsk_migration_table_copy(&table, &table2, TSK_NO_INIT);
    CU_ASSERT_EQUAL(table.metadata_schema_length, table2.metadata_schema_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_schema, table2.metadata_schema, example_length), 0);
    tsk_migration_table_set_metadata_schema(&table2, example, example_length);
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, 0));
    tsk_migration_table_set_metadata_schema(&table2, example2, example2_length);
    CU_ASSERT_FALSE(tsk_migration_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(tsk_migration_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));

    tsk_migration_table_clear(&table);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_migration_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_migration_table_free(&table2);
    CU_ASSERT_EQUAL(ret, 0);

    free(left);
    free(right);
    free(time);
    free(node);
    free(source);
    free(dest);
    free(metadata);
    free(metadata_offset);
}

static void
test_migration_table_takeset(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_migration_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    double *left;
    double *right;
    tsk_id_t *node;
    tsk_id_t *source;
    tsk_id_t *dest;
    double *time;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    tsk_size_t zeros[num_rows + 1];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    /* Make a table to copy from */
    ret = tsk_migration_table_init(&source_table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_migration_table_add_row(&source_table, (double) j, (double) j + 1,
            j + 2, j + 3, j + 4, (double) j + 5, test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    left = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    tsk_memcpy(left, source_table.left, num_rows * sizeof(double));
    right = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    tsk_memcpy(right, source_table.right, num_rows * sizeof(double));
    node = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    tsk_memcpy(node, source_table.node, num_rows * sizeof(tsk_id_t));
    source = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(source != NULL);
    tsk_memcpy(source, source_table.source, num_rows * sizeof(tsk_id_t));
    dest = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(dest != NULL);
    tsk_memcpy(dest, source_table.dest, num_rows * sizeof(tsk_id_t));
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    tsk_memcpy(time, source_table.time, num_rows * sizeof(double));
    metadata = tsk_malloc(num_rows * test_metadata_length * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    tsk_memcpy(
        metadata, source_table.metadata, num_rows * test_metadata_length * sizeof(char));
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    tsk_memcpy(metadata_offset, source_table.metadata_offset,
        (num_rows + 1) * sizeof(tsk_size_t));

    ret = tsk_migration_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_migration_table_add_row(
        &table, 1, 1, 1, 1, 1, 1, test_metadata, test_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_migration_table_takeset_columns(&table, num_rows, left, right, node,
        source, dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_migration_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_migration_table_takeset_columns(&table, num_rows, NULL, right, node,
        source, dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_takeset_columns(&table, num_rows, left, NULL, node, source,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_takeset_columns(&table, num_rows, left, right, NULL,
        source, dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_takeset_columns(&table, num_rows, left, right, node, NULL,
        dest, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_takeset_columns(&table, num_rows, left, right, node,
        source, NULL, time, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_takeset_columns(&table, num_rows, left, right, node,
        source, dest, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_takeset_columns(
        &table, num_rows, left, right, node, source, dest, time, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_table_takeset_columns(
        &table, num_rows, left, right, node, source, dest, time, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_migration_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    // Re init non-optional arrays
    left = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    tsk_memcpy(left, source_table.left, num_rows * sizeof(double));
    right = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    tsk_memcpy(right, source_table.right, num_rows * sizeof(double));
    node = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    tsk_memcpy(node, source_table.node, num_rows * sizeof(tsk_id_t));
    source = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(source != NULL);
    tsk_memcpy(source, source_table.source, num_rows * sizeof(tsk_id_t));
    dest = tsk_malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(dest != NULL);
    tsk_memcpy(dest, source_table.dest, num_rows * sizeof(tsk_id_t));
    time = tsk_malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    tsk_memcpy(time, source_table.time, num_rows * sizeof(double));
    /* if metadata and offset are both null, all entries are zero length */
    num_rows = 10;
    ret = tsk_migration_table_takeset_columns(
        &table, num_rows, left, right, node, source, dest, time, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    ret = tsk_migration_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_migration_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_migration_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_migration_table_t table;
    tsk_migration_t row;
    const char *metadata = "ABC";

    ret = tsk_migration_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_migration_table_add_row(&table, 0, 1.0, 2, 3, 4, 5, metadata, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_migration_table_add_row(&table, 1, 2.0, 3, 4, 5, 6, metadata, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_migration_table_add_row(&table, 2, 3.0, 4, 5, 6, 7, metadata, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_migration_table_update_row(&table, 0, 1, 2.0, 3, 4, 5, 6, &metadata[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 1);
    CU_ASSERT_EQUAL_FATAL(row.right, 2.0);
    CU_ASSERT_EQUAL_FATAL(row.node, 3);
    CU_ASSERT_EQUAL_FATAL(row.source, 4);
    CU_ASSERT_EQUAL_FATAL(row.dest, 5);
    CU_ASSERT_EQUAL_FATAL(row.time, 6);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_migration_table_update_row(&table, 0, row.left + 1, row.right + 1,
        row.node + 1, row.source + 1, row.dest + 1, row.time + 1, row.metadata,
        row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 2);
    CU_ASSERT_EQUAL_FATAL(row.right, 3.0);
    CU_ASSERT_EQUAL_FATAL(row.node, 4);
    CU_ASSERT_EQUAL_FATAL(row.source, 5);
    CU_ASSERT_EQUAL_FATAL(row.dest, 6);
    CU_ASSERT_EQUAL_FATAL(row.time, 7);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_migration_table_update_row(&table, 0, 0, 0, 0, 0, 0, 0, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 0);
    CU_ASSERT_EQUAL_FATAL(row.right, 0);
    CU_ASSERT_EQUAL_FATAL(row.node, 0);
    CU_ASSERT_EQUAL_FATAL(row.source, 0);
    CU_ASSERT_EQUAL_FATAL(row.dest, 0);
    CU_ASSERT_EQUAL_FATAL(row.time, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_migration_table_update_row(&table, 1, 0, 0, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 0);
    CU_ASSERT_EQUAL_FATAL(row.right, 0);
    CU_ASSERT_EQUAL_FATAL(row.node, 0);
    CU_ASSERT_EQUAL_FATAL(row.source, 0);
    CU_ASSERT_EQUAL_FATAL(row.dest, 0);
    CU_ASSERT_EQUAL_FATAL(row.time, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_migration_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.left, 2);
    CU_ASSERT_EQUAL_FATAL(row.right, 3.0);
    CU_ASSERT_EQUAL_FATAL(row.node, 4);
    CU_ASSERT_EQUAL_FATAL(row.source, 5);
    CU_ASSERT_EQUAL_FATAL(row.dest, 6);
    CU_ASSERT_EQUAL_FATAL(row.time, 7);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_migration_table_update_row(&table, 3, 0, 0, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);

    tsk_migration_table_free(&table);
}

static void
test_individual_table(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_individual_table_t table, table2;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    tsk_size_t k;
    tsk_flags_t *flags;
    double *location;
    tsk_id_t *parents;
    char *metadata;
    tsk_size_t *metadata_offset;
    tsk_size_t *parents_offset;
    tsk_size_t *location_offset;
    tsk_individual_t individual;
    tsk_individual_t individual2;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    char metadata_copy[test_metadata_length + 1];
    tsk_size_t spatial_dimension = 2;
    tsk_size_t num_parents = 3;
    double test_location[spatial_dimension];
    tsk_id_t test_parents[num_parents];
    tsk_size_t zeros[num_rows + 1];
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    for (k = 0; k < spatial_dimension; k++) {
        test_location[k] = (double) k;
    }
    for (k = 0; k < num_parents; k++) {
        test_parents[k] = (tsk_id_t) k + 42;
    }
    metadata_copy[test_metadata_length] = '\0';
    ret = tsk_individual_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_individual_table_set_max_rows_increment(&table, 1);
    tsk_individual_table_set_max_metadata_length_increment(&table, 1);
    tsk_individual_table_set_max_location_length_increment(&table, 1);
    tsk_individual_table_set_max_parents_length_increment(&table, 1);

    tsk_individual_table_print_state(&table, _devnull);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_individual_table_add_row(&table, (tsk_flags_t) j, test_location,
            spatial_dimension, test_parents, num_parents, test_metadata,
            test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
        CU_ASSERT_EQUAL(table.flags[j], (tsk_flags_t) j);
        for (k = 0; k < spatial_dimension; k++) {
            test_location[k] = (double) k;
            CU_ASSERT_EQUAL(
                table.location[spatial_dimension * (size_t) j + k], test_location[k]);
        }
        CU_ASSERT_EQUAL(
            table.metadata_length, (tsk_size_t)(j + 1) * test_metadata_length);
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], table.metadata_length);
        /* check the metadata */
        tsk_memcpy(metadata_copy, table.metadata + table.metadata_offset[j],
            test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(metadata_copy, test_metadata, test_metadata_length);

        ret = tsk_individual_table_get_row(&table, (tsk_id_t) j, &individual);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(individual.id, j);
        CU_ASSERT_EQUAL(individual.flags, (tsk_flags_t) j);
        CU_ASSERT_EQUAL(individual.location_length, spatial_dimension);
        CU_ASSERT_NSTRING_EQUAL(
            individual.location, test_location, spatial_dimension * sizeof(double));
        CU_ASSERT_EQUAL(individual.metadata_length, test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(
            individual.metadata, test_metadata, test_metadata_length);
    }

    /* Test equality with and without metadata */
    tsk_individual_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_individual_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_individual_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the metadata values */
    table2.metadata[0] = 0;
    CU_ASSERT_FALSE(tsk_individual_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_individual_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the last metadata entry */
    table2.metadata_offset[table2.num_rows]
        = table2.metadata_offset[table2.num_rows - 1];
    CU_ASSERT_FALSE(tsk_individual_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_individual_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Delete all metadata */
    tsk_memset(table2.metadata_offset, 0,
        (table2.num_rows + 1) * sizeof(*table2.metadata_offset));
    CU_ASSERT_FALSE(tsk_individual_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_individual_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    tsk_individual_table_free(&table2);

    ret = tsk_individual_table_get_row(&table, (tsk_id_t) num_rows, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_individual_table_print_state(&table, _devnull);
    tsk_individual_table_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    num_rows *= 2;
    flags = tsk_malloc(num_rows * sizeof(tsk_flags_t));
    CU_ASSERT_FATAL(flags != NULL);
    for (k = 0; k < num_rows; k++) {
        flags[k] = (tsk_flags_t)(k + num_rows);
    }
    location = tsk_malloc(spatial_dimension * num_rows * sizeof(double));
    CU_ASSERT_FATAL(location != NULL);
    for (k = 0; k < spatial_dimension * num_rows; k++) {
        location[k] = (double) (k + (num_rows * 2));
    }
    location_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(location_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        location_offset[j] = (tsk_size_t) j * spatial_dimension;
    }
    parents = tsk_malloc(num_parents * num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parents != NULL);
    for (k = 0; k < num_parents * num_rows; k++) {
        parents[k] = (tsk_id_t)(k + (num_rows * 4));
    }
    parents_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(parents_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        parents_offset[j] = (tsk_size_t) j * num_parents;
    }
    metadata = tsk_malloc(num_rows * sizeof(char));
    for (k = 0; k < num_rows; k++) {
        metadata[k] = (char) ((k % 58) + 65);
    }
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        metadata_offset[j] = (tsk_size_t) j;
    }
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location, location,
                        spatial_dimension * num_rows * sizeof(double)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location_offset, location_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents, parents, num_parents * num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parents_offset, parents_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.location_length, spatial_dimension * num_rows);
    CU_ASSERT_EQUAL(table.parents_length, num_parents * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    tsk_individual_table_print_state(&table, _devnull);

    /* Append another num_rows onto the end */
    ret = tsk_individual_table_append_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.flags + num_rows, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location, location,
                        spatial_dimension * num_rows * sizeof(double)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location + spatial_dimension * num_rows, location,
                        spatial_dimension * num_rows * sizeof(double)),
        0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents, parents, num_parents * num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parents + num_parents * num_rows, parents,
                        num_parents * num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.parents_length, 2 * num_parents * num_rows);
    CU_ASSERT_EQUAL(table.location_length, 2 * spatial_dimension * num_rows);
    tsk_individual_table_print_state(&table, _devnull);
    ret = tsk_individual_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Truncate back to num_rows */
    ret = tsk_individual_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location, location,
                        spatial_dimension * num_rows * sizeof(double)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location_offset, location_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents, parents, num_parents * num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parents_offset, parents_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset, metadata_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.location_length, spatial_dimension * num_rows);
    CU_ASSERT_EQUAL(table.parents_length, num_parents * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    tsk_individual_table_print_state(&table, _devnull);

    ret = tsk_individual_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* flags can't be NULL */
    ret = tsk_individual_table_set_columns(&table, num_rows, NULL, location,
        location_offset, parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* location and location offset must be simultaneously NULL or not */
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, location, NULL,
        parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, NULL,
        location_offset, NULL, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* parents and parents offset must be simultaneously NULL or not */
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, location,
        location_offset, parents, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, location,
        location_offset, NULL, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* if location and location_offset are both null, all locations are zero length */
    num_rows = 10;
    ret = tsk_individual_table_set_columns(
        &table, num_rows, flags, NULL, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.location_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.location_length, 0);
    ret = tsk_individual_table_append_columns(
        &table, num_rows, flags, NULL, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.location_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location_offset + num_rows, zeros,
                        num_rows * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.location_length, 0);
    tsk_individual_table_print_state(&table, _devnull);
    ret = tsk_individual_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* if parents and parents_offset are both null, all parents are zero length */
    num_rows = 10;
    ret = tsk_individual_table_set_columns(
        &table, num_rows, flags, NULL, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.parents_length, 0);
    ret = tsk_individual_table_append_columns(
        &table, num_rows, flags, NULL, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parents_offset + num_rows, zeros,
                        num_rows * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.parents_length, 0);
    tsk_individual_table_print_state(&table, _devnull);
    ret = tsk_individual_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* if metadata and metadata_offset are both null, all metadatas are zero length */
    num_rows = 10;
    ret = tsk_individual_table_set_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, flags, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location, location,
                        spatial_dimension * num_rows * sizeof(double)),
        0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents, parents, num_parents * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    ret = tsk_individual_table_append_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location, location,
                        spatial_dimension * num_rows * sizeof(double)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.location + spatial_dimension * num_rows, location,
                        spatial_dimension * num_rows * sizeof(double)),
        0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents, parents, num_parents * num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.parents + num_parents * num_rows, parents,
                        num_parents * num_rows * sizeof(tsk_id_t)),
        0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_offset + num_rows, zeros,
                        num_rows * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    tsk_individual_table_print_state(&table, _devnull);
    tsk_individual_table_dump_text(&table, _devnull);

    /* Test extend method */
    ret = tsk_individual_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_init(&table2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_individual_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_individual_table_equals(&table, &table2, 0));
    ret = tsk_individual_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_individual_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_individual_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_individual_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    ret = tsk_individual_table_set_columns(&table2, num_rows, flags, location,
        location_offset, parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_individual_table_equals(&table, &table2, 0));
    ret = tsk_individual_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_individual_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_individual_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_individual_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_individual_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_individual_table_equals(&table, &table2, 0));
    ret = tsk_individual_table_set_columns(&table2, num_rows, flags, location,
        location_offset, parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Copy a subset */
    ret = tsk_individual_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_individual_table_equals(&table, &table2, 0));
    ret = tsk_individual_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (k = 0; k < num_row_subset; k++) {
        ret = tsk_individual_table_get_row(&table, (tsk_id_t) k, &individual);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_individual_table_get_row(&table2, row_subset[k], &individual2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(individual.flags, individual2.flags);
        CU_ASSERT_EQUAL(individual.location_length, individual2.location_length);
        CU_ASSERT_EQUAL(individual.parents_length, individual2.parents_length);
        CU_ASSERT_EQUAL(individual.metadata_length, individual2.metadata_length);
        CU_ASSERT_EQUAL(tsk_memcmp(individual.location, individual2.location,
                            individual.location_length * sizeof(*individual.location)),
            0);
        CU_ASSERT_EQUAL(tsk_memcmp(individual.parents, individual2.parents,
                            individual.parents_length * sizeof(*individual.parents)),
            0);
        CU_ASSERT_EQUAL(tsk_memcmp(individual.metadata, individual2.metadata,
                            individual.metadata_length * sizeof(*individual.metadata)),
            0);
    }

    ret = tsk_individual_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, 0);
    CU_ASSERT_EQUAL(table.metadata_schema, NULL);
    const char *example = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example);
    const char *example2 = "A different example 🎄🌳🌴🌲🎋";
    tsk_size_t example2_length = (tsk_size_t) strlen(example);
    tsk_individual_table_set_metadata_schema(&table, example, example_length);
    CU_ASSERT_EQUAL(table.metadata_schema_length, example_length);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_schema, example, example_length), 0);

    tsk_individual_table_copy(&table, &table2, TSK_NO_INIT);
    CU_ASSERT_EQUAL(table.metadata_schema_length, table2.metadata_schema_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_schema, table2.metadata_schema, example_length), 0);
    tsk_individual_table_set_metadata_schema(&table2, example, example_length);
    CU_ASSERT_TRUE(tsk_individual_table_equals(&table, &table2, 0));
    tsk_individual_table_set_metadata_schema(&table2, example2, example2_length);
    CU_ASSERT_FALSE(tsk_individual_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_individual_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));

    tsk_individual_table_clear(&table);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    ret = tsk_individual_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_individual_table_free(&table2);
    CU_ASSERT_EQUAL(ret, 0);
    free(flags);
    free(location);
    free(location_offset);
    free(parents);
    free(parents_offset);
    free(metadata);
    free(metadata_offset);
}

static void
test_individual_table_takeset(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_individual_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    tsk_size_t k;
    tsk_flags_t *flags;
    double *location;
    tsk_id_t *parents;
    char *metadata;
    tsk_size_t *metadata_offset;
    tsk_size_t *parents_offset;
    tsk_size_t *location_offset;
    tsk_size_t spatial_dimension = 2;
    tsk_size_t num_parents = 3;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    double test_location[spatial_dimension];
    tsk_id_t test_parents[num_parents];
    tsk_size_t zeros[num_rows + 1];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    /* Make a table to copy from */
    ret = tsk_individual_table_init(&source_table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (k = 0; k < spatial_dimension; k++) {
        test_location[k] = (double) k;
    }
    for (k = 0; k < num_parents; k++) {
        test_parents[k] = (tsk_id_t) k + 42;
    }
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_individual_table_add_row(&source_table, (tsk_flags_t) j,
            test_location, spatial_dimension, test_parents, num_parents, test_metadata,
            test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    flags = tsk_malloc(num_rows * sizeof(tsk_flags_t));
    CU_ASSERT_FATAL(flags != NULL);
    tsk_memcpy(flags, source_table.flags, num_rows * sizeof(tsk_flags_t));
    location = tsk_malloc(spatial_dimension * num_rows * sizeof(double));
    CU_ASSERT_FATAL(location != NULL);
    tsk_memcpy(
        location, source_table.location, spatial_dimension * num_rows * sizeof(double));
    location_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(location_offset != NULL);
    tsk_memcpy(location_offset, source_table.location_offset,
        (num_rows + 1) * sizeof(tsk_size_t));
    parents = tsk_malloc(num_parents * num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parents != NULL);
    tsk_memcpy(parents, source_table.parents, num_parents * num_rows * sizeof(tsk_id_t));
    parents_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(parents_offset != NULL);
    tsk_memcpy(parents_offset, source_table.parents_offset,
        (num_rows + 1) * sizeof(tsk_size_t));
    metadata = tsk_malloc(num_rows * test_metadata_length * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    tsk_memcpy(
        metadata, source_table.metadata, num_rows * test_metadata_length * sizeof(char));
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    tsk_memcpy(metadata_offset, source_table.metadata_offset,
        (num_rows + 1) * sizeof(tsk_size_t));

    ret = tsk_individual_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_individual_table_add_row(&table, (tsk_flags_t) 1, test_location,
        spatial_dimension, test_parents, num_parents, test_metadata,
        test_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_individual_table_takeset_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_individual_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    /* location and location offset must be simultaneously NULL or not */
    ret = tsk_individual_table_takeset_columns(&table, num_rows, flags, location, NULL,
        parents, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_table_takeset_columns(&table, num_rows, flags, NULL,
        location_offset, NULL, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* parents and parents offset must be simultaneously NULL or not */
    ret = tsk_individual_table_takeset_columns(&table, num_rows, flags, location,
        location_offset, parents, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_table_takeset_columns(&table, num_rows, flags, location,
        location_offset, NULL, parents_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_individual_table_takeset_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_table_takeset_columns(&table, num_rows, flags, location,
        location_offset, parents, parents_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_individual_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    /* if ragged array and offset are both null, all entries are zero length,
       NULL flags mean all zero entries */
    num_rows = 10;
    ret = tsk_individual_table_takeset_columns(
        &table, num_rows, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(tsk_memcmp(table.flags, zeros, num_rows * sizeof(tsk_flags_t)), 0);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.location_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.location_length, 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.parents_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)), 0);
    CU_ASSERT_EQUAL(table.parents_length, 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_offset, zeros, (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    ret = tsk_individual_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_individual_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_individual_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_individual_table_t table;
    tsk_individual_t row;
    double location[] = { 0, 1, 2 };
    tsk_id_t parents[] = { 0, 1, 2 };
    const char *metadata = "ABC";

    ret = tsk_individual_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id
        = tsk_individual_table_add_row(&table, 0, location, 1, parents, 1, metadata, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_individual_table_add_row(&table, 1, location, 2, parents, 2, metadata, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_individual_table_add_row(&table, 2, location, 3, parents, 3, metadata, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_individual_table_update_row(
        &table, 0, 1, &location[1], 1, &parents[1], 1, &metadata[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 1);
    CU_ASSERT_EQUAL_FATAL(row.location_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.location[0], 1.0);
    CU_ASSERT_EQUAL_FATAL(row.parents_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.parents[0], 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_individual_table_update_row(&table, 0, row.flags + 1, row.location,
        row.location_length, row.parents, row.parents_length, row.metadata,
        row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 2);
    CU_ASSERT_EQUAL_FATAL(row.location_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.location[0], 1.0);
    CU_ASSERT_EQUAL_FATAL(row.parents_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.parents[0], 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_individual_table_update_row(&table, 0, row.flags, location, 1, row.parents,
        row.parents_length, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 2);
    CU_ASSERT_EQUAL_FATAL(row.location_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.location[0], 0.0);
    CU_ASSERT_EQUAL_FATAL(row.parents_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.parents[0], 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_individual_table_update_row(&table, 0, row.flags, NULL, 0, row.parents,
        row.parents_length, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 2);
    CU_ASSERT_EQUAL_FATAL(row.location_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.parents_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.parents[0], 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_individual_table_update_row(
        &table, 0, 2, location, 3, parents, 3, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 2);
    CU_ASSERT_EQUAL_FATAL(row.location_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.location[0], 0);
    CU_ASSERT_EQUAL_FATAL(row.location[1], 1);
    CU_ASSERT_EQUAL_FATAL(row.location[2], 2);
    CU_ASSERT_EQUAL_FATAL(row.parents_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.parents[0], 0);
    CU_ASSERT_EQUAL_FATAL(row.parents[1], 1);
    CU_ASSERT_EQUAL_FATAL(row.parents[2], 2);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_individual_table_update_row(&table, 1, 5, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 5);
    CU_ASSERT_EQUAL_FATAL(row.location_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.parents_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_individual_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.flags, 2);
    CU_ASSERT_EQUAL_FATAL(row.location_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.location[0], 0);
    CU_ASSERT_EQUAL_FATAL(row.location[1], 1);
    CU_ASSERT_EQUAL_FATAL(row.location[2], 2);
    CU_ASSERT_EQUAL_FATAL(row.parents_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.parents[0], 0);
    CU_ASSERT_EQUAL_FATAL(row.parents[1], 1);
    CU_ASSERT_EQUAL_FATAL(row.parents[2], 2);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_individual_table_update_row(&table, 3, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    tsk_individual_table_free(&table);
}

static void
test_population_table(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_population_table_t table, table2;
    tsk_size_t num_rows = 100;
    tsk_size_t max_len = 20;
    tsk_size_t k, len;
    tsk_id_t j;
    char *metadata;
    char c[max_len + 1];
    tsk_size_t *metadata_offset;
    tsk_population_t population, population2;
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    for (j = 0; j < (tsk_id_t) max_len; j++) {
        c[j] = (char) ('A' + j);
    }

    ret = tsk_population_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_population_table_set_max_rows_increment(&table, 1);
    tsk_population_table_set_max_metadata_length_increment(&table, 1);
    tsk_population_table_print_state(&table, _devnull);
    ret = tsk_population_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Adding zero length metadata with NULL should be fine */

    ret_id = tsk_population_table_add_row(&table, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, 1);
    CU_ASSERT_EQUAL(table.metadata_offset[0], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[1], 0);
    tsk_population_table_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);

    len = 0;
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        k = TSK_MIN((tsk_size_t) j + 1, max_len);
        ret_id = tsk_population_table_add_row(&table, c, k);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
        CU_ASSERT_EQUAL(table.metadata_offset[j], len);
        CU_ASSERT_EQUAL(table.num_rows, (tsk_size_t) j + 1);
        len += k;
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], len);
        CU_ASSERT_EQUAL(table.metadata_length, len);

        ret = tsk_population_table_get_row(&table, (tsk_id_t) j, &population);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(population.id, j);
        CU_ASSERT_EQUAL(population.metadata_length, k);
        CU_ASSERT_NSTRING_EQUAL(population.metadata, c, k);
    }

    /* Test equality with and without metadata */
    tsk_population_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_population_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_population_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the metadata values */
    table2.metadata[0] = 0;
    CU_ASSERT_FALSE(tsk_population_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_population_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Change the last metadata entry */
    table2.metadata_offset[table2.num_rows]
        = table2.metadata_offset[table2.num_rows - 1];
    CU_ASSERT_FALSE(tsk_population_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_population_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    /* Delete all metadata */
    tsk_memset(table2.metadata_offset, 0,
        (table2.num_rows + 1) * sizeof(*table2.metadata_offset));
    CU_ASSERT_FALSE(tsk_population_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_population_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));
    tsk_population_table_free(&table2);

    ret = tsk_population_table_get_row(&table, (tsk_id_t) num_rows, &population);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tsk_population_table_print_state(&table, _devnull);
    ret = tsk_population_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    num_rows *= 2;
    metadata = tsk_malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        metadata[j] = 'M';
        metadata_offset[j] = (tsk_size_t) j;
    }

    metadata_offset[num_rows] = num_rows;
    ret = tsk_population_table_set_columns(&table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* Append another num_rows */
    ret = tsk_population_table_append_columns(
        &table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_population_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    ret = tsk_population_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Metadata = NULL gives an error */
    ret = tsk_population_table_set_columns(&table, num_rows, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_population_table_set_columns(&table, num_rows, metadata, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_population_table_set_columns(&table, num_rows, NULL, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Test extend method */
    ret = tsk_population_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_init(&table2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_population_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_population_table_equals(&table, &table2, 0));
    ret = tsk_population_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_population_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_population_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_population_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    ret = tsk_population_table_set_columns(&table2, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_population_table_equals(&table, &table2, 0));
    ret = tsk_population_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_population_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_population_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_population_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_population_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_population_table_equals(&table, &table2, 0));
    ret = tsk_population_table_set_columns(&table2, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Copy a subset */
    ret = tsk_population_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_population_table_equals(&table, &table2, 0));
    ret = tsk_population_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (k = 0; k < num_row_subset; k++) {
        ret = tsk_population_table_get_row(&table, (tsk_id_t) k, &population);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_population_table_get_row(&table2, row_subset[k], &population2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(population.metadata_length, population2.metadata_length);
        CU_ASSERT_EQUAL(tsk_memcmp(population.metadata, population2.metadata,
                            population.metadata_length * sizeof(*population.metadata)),
            0);
    }

    /* Test for bad offsets */
    metadata_offset[0] = 1;
    ret = tsk_population_table_set_columns(&table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    metadata_offset[0] = 0;
    metadata_offset[num_rows] = 0;
    ret = tsk_population_table_set_columns(&table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    ret = tsk_population_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.metadata_schema_length, 0);
    CU_ASSERT_EQUAL(table.metadata_schema, NULL);
    const char *example = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example);
    const char *example2 = "A different example 🎄🌳🌴🌲🎋";
    tsk_size_t example2_length = (tsk_size_t) strlen(example);
    tsk_population_table_set_metadata_schema(&table, example, example_length);
    CU_ASSERT_EQUAL(table.metadata_schema_length, example_length);
    CU_ASSERT_EQUAL(tsk_memcmp(table.metadata_schema, example, example_length), 0);

    tsk_population_table_copy(&table, &table2, TSK_NO_INIT);
    CU_ASSERT_EQUAL(table.metadata_schema_length, table2.metadata_schema_length);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.metadata_schema, table2.metadata_schema, example_length), 0);
    tsk_population_table_set_metadata_schema(&table2, example, example_length);
    CU_ASSERT_TRUE(tsk_population_table_equals(&table, &table2, 0));
    tsk_population_table_set_metadata_schema(&table2, example2, example2_length);
    CU_ASSERT_FALSE(tsk_population_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_population_table_equals(&table, &table2, TSK_CMP_IGNORE_METADATA));

    tsk_population_table_clear(&table);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_population_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    tsk_population_table_free(&table2);
    CU_ASSERT_EQUAL(ret, 0);

    free(metadata);
    free(metadata_offset);
}

static void
test_population_table_takeset(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_population_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    char *metadata;
    tsk_size_t *metadata_offset;
    const char *test_metadata = "test";
    tsk_size_t test_metadata_length = 4;
    tsk_size_t zeros[num_rows + 1];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    /* Make a table to copy from */
    ret = tsk_population_table_init(&source_table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_population_table_add_row(
            &source_table, test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    metadata = tsk_malloc(num_rows * test_metadata_length * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    tsk_memcpy(
        metadata, source_table.metadata, num_rows * test_metadata_length * sizeof(char));
    metadata_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    tsk_memcpy(metadata_offset, source_table.metadata_offset,
        (num_rows + 1) * sizeof(tsk_size_t));

    ret = tsk_population_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_population_table_add_row(&table, test_metadata, test_metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_population_table_takeset_columns(
        &table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_population_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    ret = tsk_population_table_takeset_columns(&table, num_rows, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_population_table_takeset_columns(&table, num_rows, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_population_table_takeset_columns(&table, num_rows, NULL, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Test bad offset */
    metadata_offset[0] = 1;
    ret = tsk_population_table_takeset_columns(
        &table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_population_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    ret = tsk_population_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_population_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_population_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_population_table_t table;
    tsk_population_t row;
    const char *metadata = "ABC";

    ret = tsk_population_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_population_table_add_row(&table, metadata, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&table, metadata, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&table, metadata, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_population_table_update_row(&table, 0, &metadata[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_population_table_update_row(&table, 0, row.metadata, row.metadata_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'B');

    ret = tsk_population_table_update_row(&table, 0, metadata, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_population_table_update_row(&table, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 0);

    ret = tsk_population_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.metadata_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.metadata[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.metadata[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.metadata[2], 'C');

    ret = tsk_population_table_update_row(&table, 3, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    tsk_population_table_free(&table);
}

static void
test_provenance_table(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_provenance_table_t table, table2;
    tsk_size_t num_rows = 100;
    tsk_size_t j;
    char *timestamp;
    tsk_size_t *timestamp_offset;
    const char *test_timestamp = "2017-12-06T20:40:25+00:00";
    tsk_size_t test_timestamp_length = (tsk_size_t) strlen(test_timestamp);
    char timestamp_copy[test_timestamp_length + 1];
    char *record;
    tsk_size_t *record_offset;
    const char *test_record = "{\"json\"=1234}";
    tsk_size_t test_record_length = (tsk_size_t) strlen(test_record);
    char record_copy[test_record_length + 1];
    tsk_provenance_t provenance, provenance2;
    tsk_id_t row_subset[6] = { 1, 9, 1, 0, 2, 2 };
    tsk_size_t num_row_subset = 6;

    timestamp_copy[test_timestamp_length] = '\0';
    record_copy[test_record_length] = '\0';
    ret = tsk_provenance_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_provenance_table_set_max_rows_increment(&table, 1);
    tsk_provenance_table_set_max_timestamp_length_increment(&table, 1);
    tsk_provenance_table_set_max_record_length_increment(&table, 1);
    tsk_provenance_table_print_state(&table, _devnull);
    ret = tsk_provenance_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_rows; j++) {
        ret_id = tsk_provenance_table_add_row(&table, test_timestamp,
            test_timestamp_length, test_record, test_record_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, (tsk_id_t) j);
        CU_ASSERT_EQUAL(table.timestamp_length, (j + 1) * test_timestamp_length);
        CU_ASSERT_EQUAL(table.timestamp_offset[j + 1], table.timestamp_length);
        CU_ASSERT_EQUAL(table.record_length, (j + 1) * test_record_length);
        CU_ASSERT_EQUAL(table.record_offset[j + 1], table.record_length);
        /* check the timestamp */
        tsk_memcpy(timestamp_copy, table.timestamp + table.timestamp_offset[j],
            test_timestamp_length);
        CU_ASSERT_NSTRING_EQUAL(timestamp_copy, test_timestamp, test_timestamp_length);
        /* check the record */
        tsk_memcpy(
            record_copy, table.record + table.record_offset[j], test_record_length);
        CU_ASSERT_NSTRING_EQUAL(record_copy, test_record, test_record_length);

        ret = tsk_provenance_table_get_row(&table, (tsk_id_t) j, &provenance);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(provenance.id, (tsk_id_t) j);
        CU_ASSERT_EQUAL(provenance.timestamp_length, test_timestamp_length);
        CU_ASSERT_NSTRING_EQUAL(
            provenance.timestamp, test_timestamp, test_timestamp_length);
        CU_ASSERT_EQUAL(provenance.record_length, test_record_length);
        CU_ASSERT_NSTRING_EQUAL(provenance.record, test_record, test_record_length);
    }
    ret = tsk_provenance_table_get_row(&table, (tsk_id_t) num_rows, &provenance);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);
    tsk_provenance_table_print_state(&table, _devnull);
    ret = tsk_provenance_table_dump_text(&table, _devnull);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_provenance_table_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.timestamp_length, 0);
    CU_ASSERT_EQUAL(table.record_length, 0);

    num_rows *= 2;
    timestamp = tsk_malloc(num_rows * sizeof(char));
    tsk_memset(timestamp, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(timestamp != NULL);
    timestamp_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(timestamp_offset != NULL);
    record = tsk_malloc(num_rows * sizeof(char));
    tsk_memset(record, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(record != NULL);
    record_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(record_offset != NULL);
    for (j = 0; j < num_rows + 1; j++) {
        timestamp_offset[j] = j;
        record_offset[j] = j;
    }
    ret = tsk_provenance_table_set_columns(
        &table, num_rows, timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.timestamp, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.timestamp_offset, timestamp_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.record, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.record_offset, record_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.timestamp_length, num_rows);
    CU_ASSERT_EQUAL(table.record_length, num_rows);
    tsk_provenance_table_print_state(&table, _devnull);

    /* Append another num_rows onto the end */
    ret = tsk_provenance_table_append_columns(
        &table, num_rows, timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.timestamp, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.timestamp + num_rows, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.record, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(table.record + num_rows, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.timestamp_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.record_length, 2 * num_rows);
    tsk_provenance_table_print_state(&table, _devnull);

    /* Truncate back to num_rows */
    ret = tsk_provenance_table_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.timestamp, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.timestamp_offset, timestamp_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.record, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(tsk_memcmp(table.record_offset, record_offset,
                        (num_rows + 1) * sizeof(tsk_size_t)),
        0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.timestamp_length, num_rows);
    CU_ASSERT_EQUAL(table.record_length, num_rows);
    tsk_provenance_table_print_state(&table, _devnull);

    /* Test equality with and without timestamp */
    tsk_provenance_table_copy(&table, &table2, 0);
    CU_ASSERT_TRUE(tsk_provenance_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_provenance_table_equals(&table, &table2, TSK_CMP_IGNORE_TIMESTAMPS));
    /* Change the timestamp values */
    table2.timestamp[0] = 0;
    CU_ASSERT_FALSE(tsk_provenance_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_provenance_table_equals(&table, &table2, TSK_CMP_IGNORE_TIMESTAMPS));
    /* Change the last timestamp entry */
    table2.timestamp_offset[table2.num_rows]
        = table2.timestamp_offset[table2.num_rows - 1];
    CU_ASSERT_FALSE(tsk_provenance_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_provenance_table_equals(&table, &table2, TSK_CMP_IGNORE_TIMESTAMPS));
    /* Delete all timestamps */
    tsk_memset(table2.timestamp_offset, 0,
        (table2.num_rows + 1) * sizeof(*table2.timestamp_offset));
    CU_ASSERT_FALSE(tsk_provenance_table_equals(&table, &table2, 0));
    CU_ASSERT_TRUE(
        tsk_provenance_table_equals(&table, &table2, TSK_CMP_IGNORE_TIMESTAMPS));
    tsk_provenance_table_free(&table2);

    /* Test equality with and without timestamp */
    tsk_provenance_table_copy(&table, &table2, 0);
    table2.record_length = 0;
    CU_ASSERT_FALSE(tsk_provenance_table_equals(&table, &table2, 0));
    tsk_provenance_table_free(&table2);

    ret = tsk_provenance_table_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* No arguments can be null */
    ret = tsk_provenance_table_set_columns(
        &table, num_rows, NULL, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_table_set_columns(
        &table, num_rows, timestamp, NULL, record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_table_set_columns(
        &table, num_rows, timestamp, timestamp_offset, NULL, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_table_set_columns(
        &table, num_rows, timestamp, timestamp_offset, record, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Test extend method */
    ret = tsk_provenance_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_init(&table2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Can't extend from self */
    ret = tsk_provenance_table_extend(&table, &table, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANNOT_EXTEND_FROM_SELF);

    /* Two empty tables */
    CU_ASSERT_TRUE(tsk_provenance_table_equals(&table, &table2, 0));
    ret = tsk_provenance_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_provenance_table_equals(&table, &table2, 0));

    /* Row out of bounds */
    ret = tsk_provenance_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);

    /* Num rows out of bounds */
    ret = tsk_provenance_table_extend(&table, &table2, num_rows * 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);

    /* Copy rows in order if index NULL */
    ret = tsk_provenance_table_set_columns(
        &table2, num_rows, timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_provenance_table_equals(&table, &table2, 0));
    ret = tsk_provenance_table_extend(&table, &table2, table2.num_rows, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_provenance_table_equals(&table, &table2, 0));

    /* Copy nothing if index not NULL but length zero */
    ret = tsk_provenance_table_extend(&table, &table2, 0, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_provenance_table_equals(&table, &table2, 0));

    /* Copy first N rows in order if index NULL */
    ret = tsk_provenance_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_extend(&table, &table2, num_rows / 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_truncate(&table2, num_rows / 2);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_provenance_table_equals(&table, &table2, 0));
    ret = tsk_provenance_table_set_columns(
        &table2, num_rows, timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Copy a subset */
    ret = tsk_provenance_table_truncate(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_provenance_table_equals(&table, &table2, 0));
    ret = tsk_provenance_table_extend(&table, &table2, num_row_subset, row_subset, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_row_subset; j++) {
        ret = tsk_provenance_table_get_row(&table, (tsk_id_t) j, &provenance);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_provenance_table_get_row(&table2, row_subset[j], &provenance2);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(provenance.timestamp_length, provenance2.timestamp_length);
        CU_ASSERT_EQUAL(provenance.record_length, provenance2.record_length);
        CU_ASSERT_EQUAL(tsk_memcmp(provenance.timestamp, provenance2.timestamp,
                            provenance.timestamp_length * sizeof(*provenance.timestamp)),
            0);
        CU_ASSERT_EQUAL(tsk_memcmp(provenance.record, provenance2.record,
                            provenance.record_length * sizeof(*provenance.record)),
            0);
    }

    tsk_provenance_table_free(&table);
    tsk_provenance_table_free(&table2);
    free(timestamp);
    free(timestamp_offset);
    free(record);
    free(record_offset);
}

static void
test_provenance_table_takeset(void)
{
    int ret = 0;
    tsk_id_t ret_id;
    tsk_provenance_table_t source_table, table;
    tsk_size_t num_rows = 100;
    tsk_id_t j;
    char *timestamp;
    tsk_size_t *timestamp_offset;
    char *record;
    tsk_size_t *record_offset;
    const char *test_timestamp = "red";
    tsk_size_t test_timestamp_length = 3;
    const char *test_record = "test";
    tsk_size_t test_record_length = 4;
    tsk_size_t zeros[num_rows + 1];
    tsk_id_t neg_ones[num_rows];

    tsk_memset(zeros, 0, (num_rows + 1) * sizeof(tsk_size_t));
    tsk_memset(neg_ones, 0xff, num_rows * sizeof(tsk_id_t));
    /* Make a table to copy from */
    ret = tsk_provenance_table_init(&source_table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret_id = tsk_provenance_table_add_row(&source_table, test_timestamp,
            test_timestamp_length, test_record, test_record_length);
        CU_ASSERT_EQUAL_FATAL(ret_id, j);
    }

    /* Prepare arrays to be taken */
    timestamp = tsk_malloc(num_rows * test_timestamp_length * sizeof(char));
    CU_ASSERT_FATAL(timestamp != NULL);
    tsk_memcpy(timestamp, source_table.timestamp,
        num_rows * test_timestamp_length * sizeof(char));
    timestamp_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(timestamp_offset != NULL);
    tsk_memcpy(timestamp_offset, source_table.timestamp_offset,
        (num_rows + 1) * sizeof(tsk_size_t));
    record = tsk_malloc(num_rows * test_record_length * sizeof(char));
    CU_ASSERT_FATAL(record != NULL);
    tsk_memcpy(
        record, source_table.record, num_rows * test_record_length * sizeof(char));
    record_offset = tsk_malloc((num_rows + 1) * sizeof(tsk_size_t));
    CU_ASSERT_FATAL(record_offset != NULL);
    tsk_memcpy(
        record_offset, source_table.record_offset, (num_rows + 1) * sizeof(tsk_size_t));

    ret = tsk_provenance_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add one row so that we can check takeset frees it */
    ret_id = tsk_provenance_table_add_row(
        &table, test_timestamp, test_timestamp_length, test_record, test_record_length);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    ret = tsk_provenance_table_takeset_columns(
        &table, num_rows, timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_provenance_table_equals(&source_table, &table, 0));

    /* Test error states, all of these must not take the array, or free existing */
    ret = tsk_provenance_table_takeset_columns(
        &table, num_rows, NULL, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_table_takeset_columns(
        &table, num_rows, timestamp, NULL, record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_table_takeset_columns(
        &table, num_rows, timestamp, timestamp_offset, NULL, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_table_takeset_columns(
        &table, num_rows, timestamp, timestamp_offset, record, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Bad offsets */
    timestamp_offset[0] = 1;
    ret = tsk_provenance_table_takeset_columns(
        &table, num_rows, timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    timestamp_offset[0] = 0;
    record_offset[0] = 1;
    ret = tsk_provenance_table_takeset_columns(
        &table, num_rows, timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    /* Truncation after takeset keeps memory and max_rows */
    ret = tsk_provenance_table_clear(&table);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(table.max_rows, num_rows);

    ret = tsk_provenance_table_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_provenance_table_free(&source_table);
    CU_ASSERT_EQUAL(ret, 0);
}

static void
test_provenance_table_update_row(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_provenance_table_t table;
    tsk_provenance_t row;
    const char *timestamp = "XYZ";
    const char *record = "ABC";

    ret = tsk_provenance_table_init(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_provenance_table_add_row(&table, timestamp, 1, record, 1);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_provenance_table_add_row(&table, timestamp, 2, record, 2);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_provenance_table_add_row(&table, timestamp, 3, record, 3);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_provenance_table_update_row(&table, 0, &timestamp[1], 1, &record[1], 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.timestamp_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.timestamp[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.record_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.record[0], 'B');

    ret = tsk_provenance_table_update_row(
        &table, 0, row.timestamp, row.timestamp_length, row.record, row.record_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.timestamp_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.timestamp[0], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.record_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.record[0], 'B');

    ret = tsk_provenance_table_update_row(&table, 0, row.timestamp,
        row.timestamp_length - 1, row.record, row.record_length);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.timestamp_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.record_length, 1);
    CU_ASSERT_EQUAL_FATAL(row.record[0], 'B');

    ret = tsk_provenance_table_update_row(&table, 0, timestamp, 3, record, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_get_row(&table, 0, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.timestamp_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.timestamp[0], 'X');
    CU_ASSERT_EQUAL_FATAL(row.timestamp[1], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.timestamp[2], 'Z');
    CU_ASSERT_EQUAL_FATAL(row.record_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.record[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.record[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.record[2], 'C');

    ret = tsk_provenance_table_update_row(&table, 1, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_get_row(&table, 1, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.timestamp_length, 0);
    CU_ASSERT_EQUAL_FATAL(row.record_length, 0);

    ret = tsk_provenance_table_get_row(&table, 2, &row);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(row.timestamp_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.timestamp[0], 'X');
    CU_ASSERT_EQUAL_FATAL(row.timestamp[1], 'Y');
    CU_ASSERT_EQUAL_FATAL(row.timestamp[2], 'Z');
    CU_ASSERT_EQUAL_FATAL(row.record_length, 3);
    CU_ASSERT_EQUAL_FATAL(row.record[0], 'A');
    CU_ASSERT_EQUAL_FATAL(row.record[1], 'B');
    CU_ASSERT_EQUAL_FATAL(row.record[2], 'C');

    ret = tsk_provenance_table_update_row(&table, 3, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);

    tsk_provenance_table_free(&table);
}

static void
test_table_size_increments(void)
{
    int ret;
    tsk_table_collection_t tables;
    tsk_size_t new_size;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_metadata_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_location_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_metadata_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_metadata_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_ancestral_state_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_metadata_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_derived_state_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_metadata_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_metadata_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_timestamp_length_increment, 0);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_record_length_increment, 0);

    /* Setting to non-zero sets to that size */
    new_size = 1;
    ret = tsk_individual_table_set_max_rows_increment(&tables.individuals, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows_increment, new_size);
    ret = tsk_individual_table_set_max_metadata_length_increment(
        &tables.individuals, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_metadata_length_increment, new_size);
    ret = tsk_individual_table_set_max_location_length_increment(
        &tables.individuals, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_location_length_increment, new_size);

    ret = tsk_node_table_set_max_rows_increment(&tables.nodes, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows_increment, new_size);
    ret = tsk_node_table_set_max_metadata_length_increment(&tables.nodes, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length_increment, new_size);

    ret = tsk_edge_table_set_max_rows_increment(&tables.edges, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows_increment, new_size);
    ret = tsk_edge_table_set_max_metadata_length_increment(&tables.edges, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_metadata_length_increment, new_size);

    ret = tsk_site_table_set_max_rows_increment(&tables.sites, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows_increment, new_size);
    ret = tsk_site_table_set_max_metadata_length_increment(&tables.sites, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_metadata_length_increment, new_size);
    ret = tsk_site_table_set_max_ancestral_state_length_increment(
        &tables.sites, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_ancestral_state_length_increment, new_size);

    ret = tsk_mutation_table_set_max_rows_increment(&tables.mutations, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows_increment, new_size);
    ret = tsk_mutation_table_set_max_metadata_length_increment(
        &tables.mutations, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_metadata_length_increment, new_size);
    ret = tsk_mutation_table_set_max_derived_state_length_increment(
        &tables.mutations, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_derived_state_length_increment, new_size);

    ret = tsk_migration_table_set_max_rows_increment(&tables.migrations, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows_increment, new_size);
    ret = tsk_migration_table_set_max_metadata_length_increment(
        &tables.migrations, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_metadata_length_increment, new_size);

    ret = tsk_population_table_set_max_rows_increment(&tables.populations, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows_increment, new_size);
    ret = tsk_population_table_set_max_metadata_length_increment(
        &tables.populations, new_size);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_metadata_length_increment, new_size);

    ret = tsk_provenance_table_set_max_rows_increment(&tables.provenances, new_size);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows_increment, new_size);
    ret = tsk_provenance_table_set_max_timestamp_length_increment(
        &tables.provenances, new_size);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_timestamp_length_increment, new_size);
    ret = tsk_provenance_table_set_max_record_length_increment(
        &tables.provenances, new_size);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_record_length_increment, new_size);

    tsk_table_collection_free(&tables);
}

static void
test_table_expansion(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_table_collection_t tables2;

    ret = tsk_table_collection_init(&tables2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Individual table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 1);

    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 2097153);
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_individual_table_set_max_rows_increment(&tables.individuals, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_individual_table_set_max_rows_increment(&tables.individuals, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_individual_table_extend(
        &tables.individuals, &tables2.individuals, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);

    /* Node table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 1);

    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 2097153);
    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_node_table_set_max_rows_increment(&tables.nodes, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_node_table_extend(&tables.nodes, &tables2.nodes, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_node_table_set_max_rows_increment(&tables.nodes, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_extend(
        &tables.nodes, &tables2.nodes, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);

    /* Edge table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 1);

    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 2097153);
    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_edge_table_set_max_rows_increment(&tables.edges, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_edge_table_extend(&tables.edges, &tables2.edges, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_edge_table_set_max_rows_increment(&tables.edges, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_extend(
        &tables.edges, &tables2.edges, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.edges.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);

    /* Migration table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 1);

    ret_id = tsk_migration_table_add_row(&tables.migrations, 0, 0, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 2097153);
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_migration_table_set_max_rows_increment(&tables.migrations, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_migration_table_set_max_rows_increment(&tables.migrations, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_extend(
        &tables.migrations, &tables2.migrations, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.migrations.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);

    /* Site table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 1);

    ret_id = tsk_site_table_add_row(&tables.sites, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 2097153);
    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_site_table_set_max_rows_increment(&tables.sites, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_site_table_extend(&tables.sites, &tables2.sites, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_site_table_set_max_rows_increment(&tables.sites, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_extend(
        &tables.sites, &tables2.sites, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.sites.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);

    /* Mutation table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 1);

    ret_id = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 0, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_mutation_table_extend(&tables.mutations, &tables2.mutations, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 2097153);
    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_mutation_table_set_max_rows_increment(&tables.mutations, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_mutation_table_set_max_rows_increment(&tables.mutations, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_extend(
        &tables.mutations, &tables2.mutations, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);

    /* Population table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 1);

    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 2097153);
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_population_table_set_max_rows_increment(&tables.populations, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_population_table_set_max_rows_increment(&tables.populations, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_extend(
        &tables.populations, &tables2.populations, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.populations.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);

    /* Provenance table */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 1);

    ret_id = tsk_provenance_table_add_row(&tables.provenances, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);

    /*Extending by a small amount results in 1024 rows in the first case*/
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 1024);

    /*Extending by an amount that fits doesn't grow the table*/
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 1023, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 1024);

    /*Extending by an amount that doesn't fit doubles the table*/
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 1024, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 2048);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 4096, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 4097);

    /*After extending beyond 2^21 subsequent extension doesn't double but adds 2^21*/
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 2097152, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 2097153);
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 2097154, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 4194305);

    /*Extending by more rows than possible results in overflow*/
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, TSK_MAX_ID + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 4194305);

    /*Setting a custom extension uses that*/
    ret = tsk_provenance_table_set_max_rows_increment(&tables.provenances, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 4194305, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 4194305 + 42);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_provenance_table_set_max_rows_increment(&tables.provenances, TSK_MAX_ID);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_provenance_table_extend(
        &tables.provenances, &tables2.provenances, 4194305 + 42 + 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TABLE_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.provenances.max_rows, 4194305 + 42);

    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&tables2);
}

static void
test_ragged_expansion(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    char *data = tsk_malloc(104857600 * sizeof(char));

    /* Test with node table metadata */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 1);

    /*Extending by a small amount results in 65536 bytes in the first case*/
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, 2);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 65536);

    /*Extending by an amount that fits doesn't grow the column*/
    ret_id
        = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, 65534);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 65536);

    /*Extending by an amount that doesn't fit doubles the column*/
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 65536 * 2);

    /*Extending by an amount greater than the next double extends to that amount*/
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data,
        1 + (65536 * 2 * 2 - 2 - 65534 - 1));
    CU_ASSERT_EQUAL_FATAL(ret_id, 3);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 2 + 65534 + 1 + 196608);

    /*After extending beyond 100MB subsequent extension doesn't double but adds 100MB*/
    ret_id = tsk_node_table_add_row(
        &tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, 104857600);
    CU_ASSERT_EQUAL_FATAL(ret_id, 4);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 105119745);
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 5);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 105119745 + 104857600);

    /*Extending by more bytes than possible results in overflow*/
    ret_id = tsk_node_table_add_row(
        &tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, TSK_MAX_SIZE);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 105119745 + 104857600);

    tsk_node_table_free(&tables.nodes);
    ret = tsk_node_table_init(&tables.nodes, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /*Setting a custom extension uses that*/
    ret = tsk_node_table_set_max_metadata_length_increment(&tables.nodes, 42);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, 3);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 43);

    /*Setting a custom extension that overflows errors*/
    ret = tsk_node_table_set_max_metadata_length_increment(&tables.nodes, TSK_MAX_SIZE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, data, 41);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.max_metadata_length, 43);

    tsk_table_collection_free(&tables);
    tsk_safe_free(data);
}

static void
test_link_ancestors_input_errors(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_edge_table_t result;
    tsk_id_t samples[] = { 0, 1 };
    tsk_id_t ancestors[] = { 4, 6 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add an edge with some metadata */
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 7);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 1, 7, 6, "metadata", 8);
    CU_ASSERT_FATAL(ret_id > 0);

    ret = tsk_edge_table_init(&result, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_link_ancestors(
        &tables, NULL, 2, ancestors, 2, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANT_PROCESS_EDGES_WITH_METADATA);
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
    tsk_edge_table_free(&result);

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_init(&result, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_link_ancestors(
        &tables, NULL, 2, ancestors, 2, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Bad sample IDs */
    samples[0] = -1;
    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 2, ancestors, 2, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    /* Bad ancestor IDs */
    samples[0] = 0;
    ancestors[0] = -1;
    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 2, ancestors, 2, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    /* Duplicate sample IDs */
    ancestors[0] = 4;
    samples[0] = 1;
    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 2, ancestors, 2, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    /* Duplicate sample IDs */
    ancestors[0] = 6;
    samples[0] = 0;
    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 2, ancestors, 2, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);

    /* TODO more tests! */

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
    tsk_edge_table_free(&result);
}

static void
test_link_ancestors_single_tree(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_edge_table_t result;
    tsk_id_t samples[] = { 0, 1 };
    tsk_id_t ancestors[] = { 4, 6 };
    size_t i;
    double res_left = 0;
    double res_right = 1;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_init(&result, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 2, ancestors, 2, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // Check we get the right result.
    CU_ASSERT_EQUAL(result.num_rows, 3);
    tsk_id_t res_parent[] = { 4, 4, 6 };
    tsk_id_t res_child[] = { 0, 1, 4 };
    for (i = 0; i < result.num_rows; i++) {
        CU_ASSERT_EQUAL(res_parent[i], result.parent[i]);
        CU_ASSERT_EQUAL(res_child[i], result.child[i]);
        CU_ASSERT_EQUAL(res_left, result.left[i]);
        CU_ASSERT_EQUAL(res_right, result.right[i]);
    }

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
    tsk_edge_table_free(&result);
}

static void
test_link_ancestors_no_edges(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_edge_table_t result;
    tsk_id_t samples[] = { 2 };
    tsk_id_t ancestors[] = { 4 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_init(&result, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 1, ancestors, 1, 0, &result);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_table_collection_free(&tables);
    tsk_edge_table_free(&result);
    tsk_treeseq_free(&ts);
}

static void
test_link_ancestors_samples_and_ancestors_overlap(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_edge_table_t result;
    tsk_id_t samples[] = { 0, 1, 2, 4 };
    tsk_id_t ancestors[] = { 4 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_init(&result, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 4, ancestors, 1, 0, &result);

    // tsk_edge_table_print_state(&result, stdout);

    CU_ASSERT_EQUAL_FATAL(ret, 0);
    // Check we get the right result.
    CU_ASSERT_EQUAL(result.num_rows, 2);
    size_t i;
    tsk_id_t res_parent = 4;
    tsk_id_t res_child[] = { 0, 1 };
    double res_left = 0;
    double res_right = 1;
    for (i = 0; i < result.num_rows; i++) {
        CU_ASSERT_EQUAL(res_parent, result.parent[i]);
        CU_ASSERT_EQUAL(res_child[i], result.child[i]);
        CU_ASSERT_EQUAL(res_left, result.left[i]);
        CU_ASSERT_EQUAL(res_right, result.right[i]);
    }

    tsk_table_collection_free(&tables);
    tsk_edge_table_free(&result);
    tsk_treeseq_free(&ts);
}

static void
test_link_ancestors_paper(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_edge_table_t result;
    tsk_id_t samples[] = { 0, 1, 2 };
    tsk_id_t ancestors[] = { 5, 6, 7 };

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_init(&result, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 3, ancestors, 3, 0, &result);

    // tsk_edge_table_print_state(&result, stdout);

    // Check we get the right result.
    CU_ASSERT_EQUAL(result.num_rows, 6);
    size_t i;
    tsk_id_t res_parent[] = { 5, 5, 6, 6, 7, 7 };
    tsk_id_t res_child[] = { 1, 2, 0, 5, 0, 5 };
    double res_left[] = { 0, 2, 0, 0, 7, 7 };
    double res_right[] = { 10, 10, 7, 7, 10, 10 };
    for (i = 0; i < result.num_rows; i++) {
        CU_ASSERT_EQUAL(res_parent[i], result.parent[i]);
        CU_ASSERT_EQUAL(res_child[i], result.child[i]);
        CU_ASSERT_EQUAL(res_left[i], result.left[i]);
        CU_ASSERT_EQUAL(res_right[i], result.right[i]);
    }

    tsk_table_collection_free(&tables);
    tsk_edge_table_free(&result);
    tsk_treeseq_free(&ts);
}

static void
test_link_ancestors_multiple_to_single_tree(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_edge_table_t result;
    tsk_id_t samples[] = { 1, 3 };
    tsk_id_t ancestors[] = { 5 };
    size_t i;
    tsk_id_t res_parent = 5;
    tsk_id_t res_child[] = { 1, 3 };
    double res_left = 0;
    double res_right = 10;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_init(&result, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_link_ancestors(
        &tables, samples, 2, ancestors, 1, 0, &result);

    CU_ASSERT_EQUAL(result.num_rows, 2);
    for (i = 0; i < result.num_rows; i++) {
        CU_ASSERT_EQUAL(res_parent, result.parent[i]);
        CU_ASSERT_EQUAL(res_child[i], result.child[i]);
        CU_ASSERT_EQUAL(res_left, result.left[i]);
        CU_ASSERT_EQUAL(res_right, result.right[i]);
    }

    tsk_table_collection_free(&tables);
    tsk_edge_table_free(&result);
    tsk_treeseq_free(&ts);
}

static void
verify_ibd_segment_list(tsk_identity_segment_list_t *list, tsk_size_t num_nodes)
{
    tsk_identity_segment_t *seg;
    double total_span = 0;
    tsk_size_t num_segments = 0;
    /* double last_right = 0; */

    for (seg = list->head; seg != NULL; seg = seg->next) {
        CU_ASSERT_FATAL(seg->left < seg->right);
        CU_ASSERT_FATAL(seg->node >= 0);
        CU_ASSERT_FATAL(seg->node < (tsk_id_t) num_nodes);
        total_span += seg->right - seg->left;
        num_segments++;

        /* TODO the segments are not necessarily in order - issue #1682 */
        /* CU_ASSERT_FATAL(seg->left >= last_right); */
        /* last_right = seg->right; */
    }
    CU_ASSERT_EQUAL_FATAL(total_span, list->total_span);
    CU_ASSERT_EQUAL_FATAL(num_segments, list->num_segments);
}

static void
verify_ibd_result(tsk_identity_segments_t *result)
{
    int ret;
    tsk_size_t j;
    tsk_id_t a, b;
    int64_t index;
    tsk_size_t total_segments = 0;
    double total_span = 0;
    tsk_size_t num_pairs = tsk_identity_segments_get_num_pairs(result);
    tsk_id_t *pairs
        = tsk_malloc(2 * tsk_identity_segments_get_num_pairs(result) * sizeof(*pairs));
    tsk_id_t *pairs2
        = tsk_malloc(2 * tsk_identity_segments_get_num_pairs(result) * sizeof(*pairs));
    tsk_identity_segment_list_t **lists
        = tsk_malloc(tsk_identity_segments_get_num_pairs(result) * sizeof(*lists));
    tsk_avl_node_int_t **avl_nodes
        = tsk_malloc(result->pair_map.size * sizeof(*avl_nodes));

    CU_ASSERT_FATAL(pairs != NULL);
    CU_ASSERT_FATAL(pairs2 != NULL);
    CU_ASSERT_FATAL(avl_nodes != NULL);
    CU_ASSERT_FATAL(lists != NULL);
    CU_ASSERT_EQUAL_FATAL(num_pairs, result->pair_map.size);
    tsk_identity_segments_print_state(result, _devnull);

    ret = tsk_identity_segments_get_keys(result, pairs);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_avl_tree_int_ordered_nodes(&result->pair_map, avl_nodes);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_pairs; j++) {
        a = pairs[2 * j];
        b = pairs[2 * j + 1];
        index = a * (int64_t) result->num_nodes + b;
        CU_ASSERT(a < b);
        CU_ASSERT_EQUAL(tsk_avl_tree_int_search(&result->pair_map, index), avl_nodes[j]);
        index = b * (int64_t) result->num_nodes + a;
        CU_ASSERT_EQUAL(tsk_avl_tree_int_search(&result->pair_map, index), NULL);
    }

    ret = tsk_identity_segments_get_items(result, pairs2, lists);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_pairs; j++) {
        CU_ASSERT_EQUAL_FATAL(pairs[2 * j], pairs2[2 * j]);
        CU_ASSERT_EQUAL_FATAL(pairs[2 * j + 1], pairs2[2 * j + 1]);
        verify_ibd_segment_list(lists[j], result->num_nodes);
        total_segments += lists[j]->num_segments;
        total_span += lists[j]->total_span;
    }
    CU_ASSERT_EQUAL_FATAL(result->num_segments, total_segments);
    CU_ASSERT_DOUBLE_EQUAL(result->total_span, total_span, 1e-6);

    free(pairs);
    free(pairs2);
    free(lists);
    free(avl_nodes);
}

static void
test_ibd_segments_debug(void)
{
    tsk_treeseq_t ts;
    int ret;
    tsk_identity_segments_t result;
    tsk_size_t sizes[] = { 2, 2 };
    tsk_id_t samples[] = { 0, 1, 2, 3 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    tsk_set_debug_stream(_devnull);
    /* Run the DEBUG code */
    ret = tsk_table_collection_ibd_within(
        ts.tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_DEBUG);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_between(
        ts.tables, &result, 2, sizes, samples, 0.0, DBL_MAX, TSK_DEBUG);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_within(
        ts.tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_DEBUG | TSK_IBD_STORE_PAIRS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_within(
        ts.tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_DEBUG | TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_identity_segments_free(&result);

    tsk_set_debug_stream(stdout);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_caterpillar_tree(void)
{
    int ret;
    tsk_identity_segments_t result;
    tsk_treeseq_t *ts = caterpillar_tree(100, 1, 5);

    /* We're just testing out the memory expansion in ibd_finder */
    ret = tsk_table_collection_ibd_within(ts->tables, &result, NULL, 0, 0.0, DBL_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_identity_segments_free(&result);

    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_ibd_segments_single_tree(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    tsk_size_t sizes[] = { 1, 1 };
    tsk_identity_segments_t result;
    tsk_identity_segment_list_t *list = NULL;
    tsk_identity_segment_t *seg = NULL;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Only get IBD segs for (0, 1) */
    ret = tsk_table_collection_ibd_within(
        &tables, &result, samples, 2, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_identity_segments_get(&result, samples[0], samples[1], &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(list != NULL);
    seg = list->head;
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 4);

    /* Queries for other sample pairs fail */
    ret = tsk_identity_segments_get(&result, 0, 2, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(list, NULL);
    ret = tsk_identity_segments_get(&result, 1, 3, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(list, NULL);

    tsk_identity_segments_print_state(&result, _devnull);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 1);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_total_span(&result), 1);
    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);

    /* Get IBD segs among all pairs of samples */
    ret = tsk_table_collection_ibd_within(
        &tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* We have 4 samples, so 4 choose 2 sample pairs */
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 6);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_total_span(&result), 6);

    ret = tsk_identity_segments_get(&result, 0, 1, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    seg = list->head;
    CU_ASSERT_FATAL(seg != NULL);
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 4);

    ret = tsk_identity_segments_get(&result, 3, 0, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    seg = list->head;
    CU_ASSERT_FATAL(seg != NULL);
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 6);

    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);

    /* Get segs between {0} and {1} */
    ret = tsk_table_collection_ibd_between(
        ts.tables, &result, 2, sizes, samples, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_ibd_result(&result);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 1);
    ret = tsk_identity_segments_get(&result, 0, 1, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    seg = list->head;
    CU_ASSERT_FATAL(seg != NULL);
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 4);

    tsk_identity_segments_free(&result);

    /* within an empty list gives no segments */
    ret = tsk_table_collection_ibd_within(&tables, &result, samples, 0, 0.0, DBL_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 0);
    tsk_identity_segments_free(&result);

    /* Between an empty list gives no segments */
    ret = tsk_table_collection_ibd_between(
        ts.tables, &result, 0, sizes, samples, 0.0, DBL_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 0);
    tsk_identity_segments_free(&result);

    /* Between one empty list gives no segments*/
    sizes[0] = 0;
    ret = tsk_table_collection_ibd_between(
        ts.tables, &result, 2, sizes, samples, 0.0, DBL_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 0);
    tsk_identity_segments_free(&result);
    sizes[0] = 2;

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_single_tree_options(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_identity_segments_t result;
    tsk_identity_segment_list_t *list = NULL;
    tsk_id_t pairs[12];
    tsk_identity_segment_list_t *lists[6];
    tsk_flags_t options[2];
    int k;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_ibd_within(&tables, &result, NULL, 0, 0.0, DBL_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* We have 4 samples, so 4 choose 2 sample pairs */
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 6);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_total_span(&result), 6);
    /* out-of-bounds is still detected */
    ret = tsk_identity_segments_get(&result, 0, 100, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    /* By default all specific queries fail on the ibd_segments result */
    ret = tsk_identity_segments_get(&result, 0, 1, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IBD_PAIRS_NOT_STORED);
    ret = tsk_identity_segments_get_keys(&result, pairs);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IBD_PAIRS_NOT_STORED);
    ret = tsk_identity_segments_get_items(&result, pairs, lists);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IBD_PAIRS_NOT_STORED);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_within(
        &tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_IBD_STORE_PAIRS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* out-of-bounds is still detected */
    ret = tsk_identity_segments_get(&result, 0, 100, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    /* Getters for the lists now work, but the lists themselves are NULL */
    ret = tsk_identity_segments_get(&result, 0, 1, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(list->head, NULL);
    CU_ASSERT_EQUAL_FATAL(list->total_span, 1);
    CU_ASSERT_EQUAL_FATAL(list->num_segments, 1);
    ret = tsk_identity_segments_get_keys(&result, pairs);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(pairs[0], 0);
    CU_ASSERT_EQUAL_FATAL(pairs[1], 1);
    ret = tsk_identity_segments_get_items(&result, pairs, lists);
    CU_ASSERT_EQUAL_FATAL(pairs[0], 0);
    CU_ASSERT_EQUAL_FATAL(pairs[1], 1);
    CU_ASSERT_EQUAL_FATAL(lists[0]->head, NULL);
    CU_ASSERT_EQUAL_FATAL(lists[0]->total_span, 1);
    CU_ASSERT_EQUAL_FATAL(lists[0]->num_segments, 1);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_identity_segments_free(&result);

    /* store_segments implies store_pairs */
    options[0] = TSK_IBD_STORE_SEGMENTS;
    options[1] = TSK_IBD_STORE_PAIRS | TSK_IBD_STORE_SEGMENTS;
    for (k = 0; k < 2; k++) {

        ret = tsk_table_collection_ibd_within(
            &tables, &result, NULL, 0, 0.0, DBL_MAX, options[k]);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* out-of-bounds is still detected */
        ret = tsk_identity_segments_get(&result, 0, 100, &list);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
        ret = tsk_identity_segments_get(&result, 0, 1, &list);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_FATAL(list->head != NULL);
        CU_ASSERT_EQUAL_FATAL(list->head->left, 0);
        CU_ASSERT_EQUAL_FATAL(list->head->right, 1);
        CU_ASSERT_EQUAL_FATAL(list->head->next, NULL);
        CU_ASSERT_EQUAL_FATAL(list->total_span, 1);
        CU_ASSERT_EQUAL_FATAL(list->num_segments, 1);
        ret = tsk_identity_segments_get_keys(&result, pairs);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(pairs[0], 0);
        CU_ASSERT_EQUAL_FATAL(pairs[1], 1);
        ret = tsk_identity_segments_get_items(&result, pairs, lists);
        CU_ASSERT_EQUAL_FATAL(pairs[0], 0);
        CU_ASSERT_EQUAL_FATAL(pairs[1], 1);
        CU_ASSERT_FATAL(lists[0]->head != NULL);
        CU_ASSERT_EQUAL_FATAL(lists[0]->head->left, 0);
        CU_ASSERT_EQUAL_FATAL(lists[0]->head->right, 1);
        CU_ASSERT_EQUAL_FATAL(lists[0]->head->next, NULL);
        CU_ASSERT_EQUAL_FATAL(lists[0]->total_span, 1);
        CU_ASSERT_EQUAL_FATAL(lists[0]->num_segments, 1);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_identity_segments_free(&result);
    }

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_single_tree_between(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1, 2, 3 };
    tsk_size_t sizes[] = { 2, 2 };
    tsk_identity_segments_t result;
    tsk_identity_segment_list_t *list = NULL;
    tsk_identity_segment_t *seg = NULL;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Get segs between {0, 1} and {2, 3} */
    ret = tsk_table_collection_ibd_between(
        ts.tables, &result, 2, sizes, samples, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_ibd_result(&result);
    CU_ASSERT_EQUAL_FATAL(tsk_identity_segments_get_num_segments(&result), 4);

    ret = tsk_identity_segments_get(&result, 0, 2, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    seg = list->head;
    CU_ASSERT_FATAL(seg != NULL);
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 6);

    ret = tsk_identity_segments_get(&result, 0, 3, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    seg = list->head;
    CU_ASSERT_FATAL(seg != NULL);
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 6);

    ret = tsk_identity_segments_get(&result, 1, 2, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    seg = list->head;
    CU_ASSERT_FATAL(seg != NULL);
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 6);

    ret = tsk_identity_segments_get(&result, 1, 3, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    seg = list->head;
    CU_ASSERT_FATAL(seg != NULL);
    CU_ASSERT_EQUAL_FATAL(seg->next, NULL);
    CU_ASSERT_EQUAL_FATAL(seg->left, 0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 1);
    CU_ASSERT_EQUAL_FATAL(seg->node, 6);

    tsk_identity_segments_free(&result);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_multiple_trees(void)
{
    int ret;
    tsk_size_t j, k;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1, 2 };
    tsk_id_t pairs[][2] = { { 0, 1 }, { 0, 2 } };
    tsk_size_t num_samples = 3;
    tsk_size_t num_pairs = 2;
    tsk_identity_segments_t result;
    double true_left[2][2] = { { 0.0, 0.75 }, { 0.75, 0.0 } };
    double true_right[2][2] = { { 0.75, 1.0 }, { 1.0, 0.75 } };
    double true_node[2][2] = { { 4, 5 }, { 5, 6 } };
    tsk_identity_segment_list_t *list;
    tsk_identity_segment_t *seg;

    tsk_treeseq_from_text(&ts, 2, multiple_tree_ex_nodes, multiple_tree_ex_edges, NULL,
        NULL, NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_ibd_within(
        &tables, &result, samples, num_samples, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_pairs; j++) {
        ret = tsk_identity_segments_get(&result, pairs[j][0], pairs[j][1], &list);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL_FATAL(list->num_segments, 2);
        k = 0;
        for (seg = list->head; seg != NULL; seg = seg->next) {
            CU_ASSERT_EQUAL_FATAL(seg->left, true_left[j][k]);
            CU_ASSERT_EQUAL_FATAL(seg->right, true_right[j][k]);
            CU_ASSERT_EQUAL_FATAL(seg->node, true_node[j][k]);
            k++;
        }
        CU_ASSERT_EQUAL_FATAL(list->num_segments, k);
    }

    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_within(
        &tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_empty_result(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    tsk_identity_segments_t result;
    tsk_identity_segment_list_t *list;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_ibd_within(
        &tables, &result, samples, 1, 0.0, 0.5, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_identity_segments_get(&result, samples[0], samples[1], &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(list == NULL);

    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_min_span_max_time(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_identity_segments_t result;
    tsk_identity_segment_list_t *list;
    tsk_identity_segment_t *seg;

    tsk_treeseq_from_text(&ts, 2, multiple_tree_ex_nodes, multiple_tree_ex_edges, NULL,
        NULL, NULL, NULL, NULL, 0);

    ret = tsk_table_collection_ibd_within(
        ts.tables, &result, NULL, 0, 0.5, 3.0, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_identity_segments_get(&result, 0, 1, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(list->num_segments, 1);
    seg = list->head;
    CU_ASSERT_EQUAL_FATAL(seg->left, 0.0);
    CU_ASSERT_EQUAL_FATAL(seg->right, 0.75);
    CU_ASSERT_EQUAL_FATAL(seg->node, 4);

    ret = tsk_identity_segments_get(&result, 1, 2, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(list, NULL);

    ret = tsk_identity_segments_get(&result, 0, 2, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(list, NULL);

    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);

    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_errors(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1, 2 };
    tsk_id_t duplicate_samples[] = { 0, 1, 0 };
    tsk_id_t samples2[] = { -1, 1 };
    tsk_size_t sample_set_sizes[] = { 3 };
    tsk_identity_segments_t result;
    tsk_identity_segment_list_t *list;

    tsk_treeseq_from_text(&ts, 2, multiple_tree_ex_nodes, multiple_tree_ex_edges, NULL,
        NULL, NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // Invalid sample IDs
    ret = tsk_table_collection_ibd_within(
        &tables, &result, samples2, 1, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_between(&tables, &result, 1, sample_set_sizes,
        samples2, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_identity_segments_free(&result);

    // Bad length or time
    ret = tsk_table_collection_ibd_within(&tables, &result, samples, 2, 0.0, -1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    tsk_identity_segments_free(&result);
    ret = tsk_table_collection_ibd_within(&tables, &result, samples, 2, -1, 0.0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_between(&tables, &result, 1, sample_set_sizes,
        samples, -1, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    tsk_identity_segments_free(&result);
    ret = tsk_table_collection_ibd_between(
        &tables, &result, 1, sample_set_sizes, samples, 0, -1, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    tsk_identity_segments_free(&result);

    // Duplicate samples
    ret = tsk_table_collection_ibd_within(
        &tables, &result, duplicate_samples, 3, 0.0, DBL_MAX, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);
    tsk_identity_segments_free(&result);

    ret = tsk_table_collection_ibd_between(&tables, &result, 1, sample_set_sizes,
        duplicate_samples, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SAMPLE);
    tsk_identity_segments_free(&result);

    // Check for bad inputs to result
    ret = tsk_table_collection_ibd_within(
        &tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_identity_segments_get(&result, 0, -1, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_identity_segments_get(&result, -1, 0, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_identity_segments_get(&result, 0, 100, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    ret = tsk_identity_segments_get(&result, 100, 0, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_identity_segments_get(&result, 0, 5, &list);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(list, NULL);
    /* TODO add more checks here */
    ret = tsk_identity_segments_get(&result, 0, 0, &list);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SAME_NODES_IN_PAIR);

    tsk_identity_segments_free(&result);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_samples_are_descendants(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_id_t samples[] = { 0, 1, 2, 3, 4, 5 };
    tsk_size_t num_samples = 6;
    tsk_identity_segments_t result;
    tsk_id_t pairs[][2] = { { 0, 2 }, { 0, 4 }, { 2, 4 }, { 1, 3 }, { 1, 5 }, { 3, 5 } };
    tsk_size_t num_pairs = 6;
    tsk_id_t true_node[] = { 2, 4, 4, 3, 5, 5 };
    tsk_size_t j;
    tsk_identity_segment_list_t *list;
    tsk_identity_segment_t *seg;

    tsk_treeseq_from_text(&ts, 1, multi_root_tree_ex_nodes, multi_root_tree_ex_edges,
        NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_table_collection_ibd_within(
        ts.tables, &result, samples, num_samples, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < num_pairs; j++) {
        tsk_identity_segments_get(&result, pairs[j][0], pairs[j][1], &list);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_FATAL(list != NULL);
        CU_ASSERT_EQUAL_FATAL(list->num_segments, 1);
        seg = list->head;

        CU_ASSERT_EQUAL_FATAL(seg->left, 0);
        CU_ASSERT_EQUAL_FATAL(seg->right, 1);
        CU_ASSERT_EQUAL_FATAL(seg->node, true_node[j]);
    }

    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_multiple_ibd_paths(void)
{
    int ret;
    tsk_size_t j, k;
    tsk_treeseq_t ts;
    tsk_id_t pairs[][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
    tsk_size_t num_pairs = 3;
    tsk_identity_segments_t result;
    double true_left[3][2] = { { 0.2, 0.0 }, { 0.2, 0.0 }, { 0.0, 0.2 } };
    double true_right[3][2] = { { 1.0, 0.2 }, { 1.0, 0.2 }, { 0.2, 1.0 } };
    double true_node[3][2] = { { 4, 5 }, { 3, 5 }, { 4, 4 } };
    tsk_identity_segment_list_t *list;
    tsk_identity_segment_t *seg;

    tsk_treeseq_from_text(&ts, 2, multi_path_tree_ex_nodes, multi_path_tree_ex_edges,
        NULL, NULL, NULL, NULL, NULL, 0);

    ret = tsk_table_collection_ibd_within(
        ts.tables, &result, NULL, 0, 0.0, DBL_MAX, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < num_pairs; j++) {
        tsk_identity_segments_get(&result, pairs[j][0], pairs[j][1], &list);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        k = 0;
        for (seg = list->head; seg != NULL; seg = seg->next) {
            CU_ASSERT_EQUAL_FATAL(seg->left, true_left[j][k]);
            CU_ASSERT_EQUAL_FATAL(seg->right, true_right[j][k]);
            CU_ASSERT_EQUAL_FATAL(seg->node, true_node[j][k]);
            k++;
        }
        CU_ASSERT_EQUAL_FATAL(k, 2);
    }

    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);
    tsk_treeseq_free(&ts);
}

static void
test_ibd_segments_odd_topologies(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };
    tsk_id_t samples1[] = { 0, 2 };
    tsk_identity_segments_t result;

    tsk_treeseq_from_text(
        &ts, 1, odd_tree1_ex_nodes, odd_tree1_ex_edges, NULL, NULL, NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // Multiple roots.
    ret = tsk_table_collection_ibd_within(
        &tables, &result, samples, 1, 0, 0, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);

    // Parent is a sample.
    ret = tsk_table_collection_ibd_within(
        &tables, &result, samples1, 1, 0, 0, TSK_IBD_STORE_SEGMENTS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    verify_ibd_result(&result);
    tsk_identity_segments_free(&result);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplify_tables_drops_indexes(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_id_t samples[] = { 0, 1 };

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_table_collection_has_index(&tables, 0))
    ret = tsk_table_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0))

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_simplify_empty_tables(void)
{
    int ret;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_table_collection_simplify(&tables, NULL, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 0);

    tsk_table_collection_free(&tables);
}

static void
test_simplify_metadata(void)
{
    int ret;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 10;
    tsk_edge_table_add_row(&tables.edges, 0, 0, 1, 1, "metadata", 8);
    ret = tsk_table_collection_simplify(&tables, NULL, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_CANT_PROCESS_EDGES_WITH_METADATA);

    tsk_table_collection_free(&tables);
}

static void
test_edge_update_invalidates_index(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);

    /* Any operations on the edge table should now invalidate the index */
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_index(&tables, 0))
    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0));
    /* Even though the actual indexes still exist */
    CU_ASSERT_FALSE(tables.indexes.edge_insertion_order == NULL);
    CU_ASSERT_FALSE(tables.indexes.edge_removal_order == NULL);
    CU_ASSERT_EQUAL_FATAL(tables.indexes.num_edges, tsk_treeseq_get_num_edges(&ts));

    ret = tsk_treeseq_copy_tables(&ts, &tables, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_index(&tables, 0))
    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 1, 0, 1, NULL, 0);
    CU_ASSERT_TRUE(ret_id > 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0));
    /* Even though the actual indexes still exist */
    CU_ASSERT_FALSE(tables.indexes.edge_insertion_order == NULL);
    CU_ASSERT_FALSE(tables.indexes.edge_removal_order == NULL);
    CU_ASSERT_EQUAL_FATAL(tables.indexes.num_edges, tsk_treeseq_get_num_edges(&ts));

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_copy_table_collection(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables, tables_copy;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add some migrations, population and provenance */
    ret_id = tsk_migration_table_add_row(&tables.migrations, 0, 1, 2, 3, 4, 5, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_migration_table_add_row(&tables.migrations, 1, 2, 3, 4, 5, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_population_table_add_row(&tables.populations, "metadata", 8);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_population_table_add_row(&tables.populations, "other", 5);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_provenance_table_add_row(&tables.provenances, "time", 4, "record", 6);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_provenance_table_add_row(&tables.provenances, "time ", 5, "record ", 7);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);

    tsk_table_collection_copy(&tables, &tables_copy, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &tables_copy, 0));

    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&tables_copy);
    tsk_treeseq_free(&ts);
}

static void
test_sort_tables_offsets(void)
{
    int ret;
    tsk_treeseq_t *ts;
    tsk_table_collection_t tables, copy;
    tsk_bookmark_t bookmark;

    ts = caterpillar_tree(10, 5, 5);
    ret = tsk_treeseq_copy_tables(ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Check that setting edge offset = len(edges) does nothing */
    reverse_edges(&tables);
    ret = tsk_table_collection_copy(&tables, &copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.edges = tables.edges.num_rows;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &copy, 0));

    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Check that setting migration offset = len(migrations) does nothing */
    reverse_migrations(&tables);
    ret = tsk_table_collection_copy(&tables, &copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.migrations = tables.migrations.num_rows;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &copy, 0));

    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tables.sites.num_rows > 2);
    CU_ASSERT_FATAL(tables.mutations.num_rows > 2);

    /* Check that setting mutation and site offset = to the len
     * of the tables leaves them untouched. */
    reverse_mutations(&tables);
    /* Swap the positions of the first two sites, as a quick way
     * to disorder the site table */
    tables.sites.position[0] = tables.sites.position[1];
    tables.sites.position[1] = 0;
    ret = tsk_table_collection_copy(&tables, &copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.sites = tables.sites.num_rows;
    bookmark.mutations = tables.mutations.num_rows;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &copy, 0));

    /* Anything other than len(table) leads to an error for sites
     * and mutations, and we can't specify one without the other. */
    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.sites = tables.sites.num_rows;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SORT_OFFSET_NOT_SUPPORTED);

    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.mutations = tables.mutations.num_rows;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SORT_OFFSET_NOT_SUPPORTED);

    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.sites = tables.sites.num_rows - 1;
    bookmark.mutations = tables.mutations.num_rows - 1;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SORT_OFFSET_NOT_SUPPORTED);

    /* Individuals must either all be sorted or all skipped */
    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Add a parent relation that unsorts the table */
    tables.individuals.parents[0] = 5;
    ret = tsk_table_collection_copy(&tables, &copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.individuals = tables.individuals.num_rows;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy, 0));

    /* Check that sorting would have had no effect as individuals not in default sort*/
    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&tables, &copy, 0));

    /* Individual bookmark ignored */
    tsk_memset(&bookmark, 0, sizeof(bookmark));
    bookmark.individuals = tables.individuals.num_rows - 1;
    ret = tsk_table_collection_sort(&tables, &bookmark, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&copy);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_sort_tables_drops_indexes_with_options(tsk_flags_t tc_options)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, tc_options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_table_collection_has_index(&tables, 0))
    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0))

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_sort_tables_drops_indexes(void)
{
    test_sort_tables_drops_indexes_with_options(0);
    test_sort_tables_drops_indexes_with_options(TSK_TC_NO_EDGE_METADATA);
}

static void
test_sort_tables_edge_metadata(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t t1, t2;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    insert_edge_metadata(&t1);
    ret = tsk_table_collection_copy(&t1, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    reverse_edges(&t1);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&t1, &t2, 0));
    ret = tsk_table_collection_sort(&t1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
    tsk_treeseq_free(&ts);
}

static void
test_sort_tables_no_edge_metadata(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t t1, t2;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t1, TSK_TC_NO_EDGE_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(t1.edges.options & TSK_TABLE_NO_METADATA);
    ret = tsk_table_collection_copy(&t1, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(t2.edges.options & TSK_TABLE_NO_METADATA);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    reverse_edges(&t1);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&t1, &t2, 0));
    ret = tsk_table_collection_sort(&t1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);

    ret = tsk_table_collection_copy(&t1, &t2, TSK_TC_NO_EDGE_METADATA);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(t1.edges.options & TSK_TABLE_NO_METADATA);
    CU_ASSERT_TRUE(t2.edges.options & TSK_TABLE_NO_METADATA);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    reverse_edges(&t1);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&t1, &t2, 0));
    ret = tsk_table_collection_sort(&t1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(&ts);
}

static void
test_sort_tables_errors(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_bookmark_t pos;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL,
        single_tree_ex_sites, single_tree_ex_mutations, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_memset(&pos, 0, sizeof(pos));
    /* Everything 0 should be fine */
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Everything is sorted already */
    pos.edges = tables.edges.num_rows;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    pos.edges = (tsk_size_t) -1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);

    pos.edges = tables.edges.num_rows + 1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);

    tsk_memset(&pos, 0, sizeof(pos));
    pos.migrations = (tsk_size_t) -1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);

    pos.migrations = tables.migrations.num_rows + 1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);

    /* Node, population and provenance positions are ignored */
    tsk_memset(&pos, 0, sizeof(pos));
    pos.nodes = 1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_memset(&pos, 0, sizeof(pos));
    pos.populations = 1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_memset(&pos, 0, sizeof(pos));
    pos.provenances = 1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Specifying only one of sites or mutations is an error */
    tsk_memset(&pos, 0, sizeof(pos));
    pos.sites = 1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SORT_OFFSET_NOT_SUPPORTED);

    tsk_memset(&pos, 0, sizeof(pos));
    pos.mutations = 1;
    ret = tsk_table_collection_sort(&tables, &pos, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SORT_OFFSET_NOT_SUPPORTED);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_sort_tables_mutation_times(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables, t1, t2;
    const char *sites = "0       0\n"
                        "0.1     0\n"
                        "0.2     0\n"
                        "0.3     0\n";
    const char *mutations = "0   0  1  -1  3\n"
                            "1   1  1  -1  3\n"
                            "2   4  1  -1  8\n"
                            "2   1  0  -1   4\n"
                            "2   2  1  -1  3\n"
                            "2   1  1  -1   2\n"
                            "3   6  1  -1  10\n";

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 7);
    tables.nodes.time[4] = 6;
    tables.nodes.time[5] = 8;
    tables.nodes.time[6] = 10;
    parse_edges(single_tree_ex_edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 6);
    parse_sites(sites, &tables.sites);
    parse_mutations(mutations, &tables.mutations);
    CU_ASSERT_EQUAL_FATAL(tables.sites.num_rows, 4);
    CU_ASSERT_EQUAL_FATAL(tables.mutations.num_rows, 7);
    tables.sequence_length = 1.0;

    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check to make sure we have legal mutations */
    ret = tsk_treeseq_init(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_copy_tables(&ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_copy(&t1, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    reverse_mutations(&t1);
    CU_ASSERT_FALSE(tsk_table_collection_equals(&t1, &t2, 0));
    ret = tsk_table_collection_sort(&t1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_sort_tables_canonical_errors(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    tsk_id_t null_p[] = { -1 };
    tsk_id_t zero_p[] = { 0 };
    tsk_id_t one_p[] = { 1 };

    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.0, "x", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 2, 0.0, "a", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 3, 0.0, "b", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 1, 0.0, "c", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 2, 0.0, "d", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_table_collection_canonicalise(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_PARENT_INCONSISTENT);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_FATAL(ret == 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 2, 0.0, "a", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 3, 0.0, "b", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 1, 0.0, "c", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, -1, 0.0, "d", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_table_collection_canonicalise(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0.0, TSK_NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0.0, TSK_NULL, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, one_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, zero_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_table_collection_canonicalise(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_PARENT_CYCLE);

    ret = tsk_individual_table_clear(&tables.individuals);
    CU_ASSERT_FATAL(ret == 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, zero_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, zero_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_table_collection_canonicalise(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_SELF_PARENT);

    ret = tsk_individual_table_clear(&tables.individuals);
    CU_ASSERT_FATAL(ret == 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, null_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, zero_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_table_collection_canonicalise(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tsk_table_collection_free(&tables);
}

static void
test_sort_tables_canonical(void)
{
    int ret;
    tsk_table_collection_t t1, t2;
    // this is single_tree_ex with individuals and populations
    const char *nodes = "1  0   -1    1\n"
                        "1  0    2    3\n"
                        "1  0    0   -1\n"
                        "1  0   -1    3\n"
                        "0  1    2   -1\n"
                        "0  2   -1    2\n"
                        "0  3   -1   -1\n";
    const char *individuals = "0 0.0 1\n"
                              "0 1.0 -1\n"
                              "0 2.0 1,3\n"
                              "0 3.0 -1,1\n";
    const char *sites = "0       0\n"
                        "0.2     0\n"
                        "0.1     0\n";
    const char *mutations = "0   0  2   3 0.5\n"
                            "2   1  1  -1 0.5\n"
                            "1   4  3  -1   3\n"
                            "0   4  1  -1 2.5\n"
                            "2   2  1  -1   2\n"
                            "1   1  5   7 0.5\n"
                            "1   2  1  -1   2\n"
                            "1   1  4   2 0.5\n"
                            "1   1  6   7 0.5\n";
    const char *nodes_sorted = "1  0   -1    0\n"
                               "1  0    0    1\n"
                               "1  0    1   -1\n"
                               "1  0   -1    1\n"
                               "0  1    0   -1\n"
                               "0  2   -1    2\n"
                               "0  3   -1   -1\n";
    const char *individuals_sorted = "0 1.0 -1\n"
                                     "0 3.0 -1,0\n"
                                     "0 2.0 0,1\n";
    const char *sites_sorted = "0       0\n"
                               "0.1     0\n"
                               "0.2     0\n";
    const char *mutations_sorted = "0   4  1  -1 2.5\n"
                                   "0   0  2   0 0.5\n"
                                   "1   2  1  -1   2\n"
                                   "1   1  1  -1 0.5\n"
                                   "2   4  3  -1   3\n"
                                   "2   2  1  -1   2\n"
                                   "2   1  4   4 0.5\n"
                                   "2   1  5   6 0.5\n"
                                   "2   1  6   6 0.5\n";
    const char *individuals_sorted_kept = "0 1.0 -1\n"
                                          "0 3.0 -1,0\n"
                                          "0 2.0 0,1\n"
                                          "0 0.0 0\n";

    ret = tsk_table_collection_init(&t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    t1.sequence_length = 1.0;
    ret = tsk_table_collection_init(&t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    t2.sequence_length = 1.0;

    parse_nodes(nodes, &t1.nodes);
    CU_ASSERT_EQUAL_FATAL(t1.nodes.num_rows, 7);
    parse_individuals(individuals, &t1.individuals);
    CU_ASSERT_EQUAL_FATAL(t1.individuals.num_rows, 4);
    tsk_population_table_add_row(&t1.populations, "A", 1);
    tsk_population_table_add_row(&t1.populations, "B", 1);
    tsk_population_table_add_row(&t1.populations, "C", 1);
    parse_edges(single_tree_ex_edges, &t1.edges);
    CU_ASSERT_EQUAL_FATAL(t1.edges.num_rows, 6);
    parse_sites(sites, &t1.sites);
    CU_ASSERT_EQUAL_FATAL(t1.sites.num_rows, 3);
    parse_mutations(mutations, &t1.mutations);
    CU_ASSERT_EQUAL_FATAL(t1.mutations.num_rows, 9);

    ret = tsk_table_collection_canonicalise(&t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(nodes_sorted, &t2.nodes);
    tsk_population_table_add_row(&t2.populations, "C", 1);
    tsk_population_table_add_row(&t2.populations, "A", 1);
    CU_ASSERT_EQUAL_FATAL(t2.nodes.num_rows, 7);
    parse_individuals(individuals_sorted, &t2.individuals);
    CU_ASSERT_EQUAL_FATAL(t2.individuals.num_rows, 3);
    parse_edges(single_tree_ex_edges, &t2.edges);
    CU_ASSERT_EQUAL_FATAL(t2.edges.num_rows, 6);
    parse_sites(sites_sorted, &t2.sites);
    parse_mutations(mutations_sorted, &t2.mutations);
    CU_ASSERT_EQUAL_FATAL(t2.sites.num_rows, 3);
    CU_ASSERT_EQUAL_FATAL(t2.mutations.num_rows, 9);

    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    ret = tsk_table_collection_clear(&t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_clear(&t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // now with KEEP_UNREFERENCED
    parse_nodes(nodes, &t1.nodes);
    parse_individuals(individuals, &t1.individuals);
    tsk_population_table_add_row(&t1.populations, "A", 1);
    tsk_population_table_add_row(&t1.populations, "B", 1);
    tsk_population_table_add_row(&t1.populations, "C", 1);
    parse_edges(single_tree_ex_edges, &t1.edges);
    parse_sites(sites, &t1.sites);
    parse_mutations(mutations, &t1.mutations);

    ret = tsk_table_collection_canonicalise(&t1, TSK_SUBSET_KEEP_UNREFERENCED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_nodes(nodes_sorted, &t2.nodes);
    tsk_population_table_add_row(&t2.populations, "C", 1);
    tsk_population_table_add_row(&t2.populations, "A", 1);
    tsk_population_table_add_row(&t2.populations, "B", 1);
    parse_individuals(individuals_sorted_kept, &t2.individuals);
    CU_ASSERT_EQUAL_FATAL(t2.individuals.num_rows, 4);
    parse_edges(single_tree_ex_edges, &t2.edges);
    parse_sites(sites_sorted, &t2.sites);
    parse_mutations(mutations_sorted, &t2.mutations);

    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t2);
    tsk_table_collection_free(&t1);
}

static void
test_sort_tables_migrations(void)
{
    int ret;
    tsk_treeseq_t *ts;
    tsk_table_collection_t tables, copy;

    ts = caterpillar_tree(13, 1, 1);
    ret = tsk_treeseq_copy_tables(ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tables.migrations.num_rows > 0);

    ret = tsk_table_collection_copy(&tables, &copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &copy, 0));

    reverse_migrations(&tables);
    CU_ASSERT_FATAL(!tsk_table_collection_equals(&tables, &copy, 0));
    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_migration_table_equals(&tables.migrations, &copy.migrations, 0));
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &copy, 0));

    /* Make sure we test the deeper comparison keys. The full key is
     * (time, source, dest, left, node) */
    tsk_migration_table_clear(&tables.migrations);

    /* params = left, right, node, source, dest, time */
    tsk_migration_table_add_row(&tables.migrations, 0, 1, 0, 0, 1, 0, NULL, 0);
    tsk_migration_table_add_row(&tables.migrations, 0, 1, 1, 0, 1, 0, NULL, 0);
    ret = tsk_migration_table_copy(&tables.migrations, &copy.migrations, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    reverse_migrations(&tables);
    CU_ASSERT_FATAL(!tsk_table_collection_equals(&tables, &copy, 0));
    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_migration_table_equals(&tables.migrations, &copy.migrations, 0));

    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&copy);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_sort_tables_individuals(void)
{
    int ret;
    tsk_table_collection_t tables, copy;
    const char *individuals = "1      0.25   2,3 0\n"
                              "2      0.5    5,-1  1\n"
                              "3      0.3    -1  2\n"
                              "4      0.3    -1  3\n"
                              "5      0.3    3   4\n"
                              "6      0.3    4   5\n";
    const char *individuals_cycle = "1      0.2    2  0\n"
                                    "2      0.5    0  1\n"
                                    "3      0.3    1  2\n";
    const tsk_id_t bad_parents[] = { 200 };
    tsk_id_t ret_id;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1.0;
    parse_individuals(individuals, &tables.individuals);

    ret = tsk_table_collection_copy(&tables, &copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Table sort doesn't touch individuals by default*/
    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &copy, 0));

    /* Not calling with TSK_CHECK_TREES so casting is safe */
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_INDIVIDUAL_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_INDIVIDUALS);

    ret = tsk_table_collection_individual_topological_sort(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_INDIVIDUAL_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check that the sort is stable */
    tsk_table_collection_free(&copy);
    ret = tsk_table_collection_copy(&tables, &copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_individual_topological_sort(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &copy, 0));

    /* Errors on bad table */
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, bad_parents, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 6);
    ret = tsk_table_collection_individual_topological_sort(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    /* Errors on cycle */
    tsk_individual_table_clear(&tables.individuals);
    parse_individuals(individuals_cycle, &tables.individuals);
    ret = tsk_table_collection_individual_topological_sort(&tables, 0);
    CU_ASSERT_EQUAL(ret, TSK_ERR_INDIVIDUAL_PARENT_CYCLE);

    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&copy);
}

static void
test_sorter_interface(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;
    tsk_table_sorter_t sorter;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, &tables, 0));

    /* Nominal case */
    reverse_edges(&tables);
    CU_ASSERT_FALSE(tsk_table_collection_equals(ts.tables, &tables, 0));
    ret = tsk_table_sorter_init(&sorter, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_sorter_run(&sorter, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, &tables, 0));
    CU_ASSERT_EQUAL(sorter.user_data, NULL);
    tsk_table_sorter_free(&sorter);

    /* If we set the sort_edges function to NULL then we should leave the
     * node table as is. */
    reverse_edges(&tables);
    CU_ASSERT_FALSE(tsk_edge_table_equals(&ts.tables->edges, &tables.edges, 0));
    ret = tsk_table_sorter_init(&sorter, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    sorter.sort_edges = NULL;
    ret = tsk_table_sorter_run(&sorter, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_edge_table_equals(&ts.tables->edges, &tables.edges, 0));
    tsk_table_sorter_free(&sorter);

    /* Reversing again should make them equal */
    reverse_edges(&tables);
    CU_ASSERT_TRUE(tsk_edge_table_equals(&ts.tables->edges, &tables.edges, 0));

    /* Do not check integrity before sorting */
    reverse_edges(&tables);
    CU_ASSERT_FALSE(tsk_table_collection_equals(ts.tables, &tables, 0));
    ret = tsk_table_sorter_init(&sorter, &tables, TSK_NO_CHECK_INTEGRITY);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_sorter_run(&sorter, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, &tables, 0));
    tsk_table_sorter_free(&sorter);

    /* The user_data shouldn't be touched */
    reverse_edges(&tables);
    CU_ASSERT_FALSE(tsk_table_collection_equals(ts.tables, &tables, 0));
    ret = tsk_table_sorter_init(&sorter, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    sorter.user_data = (void *) &ts;
    ret = tsk_table_sorter_run(&sorter, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(ts.tables, &tables, 0));
    CU_ASSERT_EQUAL_FATAL(sorter.user_data, &ts);
    tsk_table_sorter_free(&sorter);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_dump_unindexed_with_options(tsk_flags_t tc_options)
{
    tsk_table_collection_t tables, loaded;
    int ret;

    ret = tsk_table_collection_init(&tables, tc_options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, &tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes.num_rows, 7);
    parse_edges(single_tree_ex_edges, &tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges.num_rows, 6);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0));
    ret = tsk_table_collection_dump(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0));

    ret = tsk_table_collection_load(&loaded, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&loaded, 0));
    CU_ASSERT_TRUE(tsk_node_table_equals(&tables.nodes, &loaded.nodes, 0));
    CU_ASSERT_TRUE(tsk_edge_table_equals(&tables.edges, &loaded.edges, 0));

    tsk_table_collection_free(&loaded);
    tsk_table_collection_free(&tables);
}

static void
test_dump_unindexed(void)
{
    test_dump_unindexed_with_options(0);
    test_dump_unindexed_with_options(TSK_TC_NO_EDGE_METADATA);
}

static void
test_dump_load_empty_with_options(tsk_flags_t tc_options)
{
    int ret;
    tsk_table_collection_t t1, t2;

    ret = tsk_table_collection_init(&t1, tc_options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    t1.sequence_length = 1.0;
    ret = tsk_table_collection_dump(&t1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
}

static void
test_dump_load_empty(void)
{
    test_dump_load_empty_with_options(0);
    test_dump_load_empty_with_options(TSK_TC_NO_EDGE_METADATA);
}

static void
test_dump_load_unsorted_with_options(tsk_flags_t tc_options)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t t1, t2;

    ret = tsk_table_collection_init(&t1, tc_options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    t1.sequence_length = 1.0;

    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 1, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 3);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 2, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 4);

    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 3, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 4, 3, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 3, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 4, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 3);

    /* Verify that it's unsorted */
    ret = (int) tsk_table_collection_check_integrity(&t1, TSK_CHECK_EDGE_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME);

    ret = tsk_table_collection_dump(&t1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&t1, 0));
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&t1, 0));
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
}

static void
test_dump_load_unsorted(void)
{
    test_dump_load_unsorted_with_options(0);
    test_dump_load_unsorted_with_options(TSK_TC_NO_EDGE_METADATA);
}

static void
test_dump_load_metadata_schema(void)
{
    int ret;
    tsk_table_collection_t t1, t2;

    ret = tsk_table_collection_init(&t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    t1.sequence_length = 1.0;
    char example[100] = "An example of metadata schema with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_length = (tsk_size_t) strlen(example) + 4;
    tsk_node_table_set_metadata_schema(
        &t1.nodes, strcat(example, "node"), example_length);
    tsk_edge_table_set_metadata_schema(
        &t1.edges, strcat(example, "edge"), example_length);
    tsk_site_table_set_metadata_schema(
        &t1.sites, strcat(example, "site"), example_length);
    tsk_mutation_table_set_metadata_schema(
        &t1.mutations, strcat(example, "muta"), example_length);
    tsk_migration_table_set_metadata_schema(
        &t1.migrations, strcat(example, "migr"), example_length);
    tsk_individual_table_set_metadata_schema(
        &t1.individuals, strcat(example, "indi"), example_length);
    tsk_population_table_set_metadata_schema(
        &t1.populations, strcat(example, "popu"), example_length);
    ret = tsk_table_collection_dump(&t1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
}

static void
test_dump_fail_no_file(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t t1;

    ret = tsk_table_collection_init(&t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    t1.sequence_length = 1.0;

    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 1, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 3);
    ret_id = tsk_node_table_add_row(
        &t1.nodes, TSK_NODE_IS_SAMPLE, 2, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 4);

    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 3, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 4, 3, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 3, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    ret_id = tsk_edge_table_add_row(&t1.edges, 0, 1, 4, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 3);

    /* Verify that it's unsorted */
    ret = (int) tsk_table_collection_check_integrity(&t1, TSK_CHECK_EDGE_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME);

    /* Make sure the file doesn't exist beforehand. */
    unlink(_tmp_file_name);
    errno = 0;

    CU_ASSERT_EQUAL(access(_tmp_file_name, F_OK), -1);

    tsk_table_collection_free(&t1);
}

static void
test_load_reindex(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_dump(&ts, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_drop_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0));
    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_index(&tables, 0));

    ret = tsk_table_collection_drop_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Dump the unindexed version */
    ret = tsk_table_collection_dump(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0));
    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_index(&tables, 0));

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_table_overflow(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_size_t max_rows = ((tsk_size_t) TSK_MAX_ID) + 1;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Simulate overflows */
    tables.individuals.max_rows = max_rows;
    tables.individuals.num_rows = max_rows;
    ret_id
        = tsk_individual_table_add_row(&tables.individuals, 0, 0, 0, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tables.nodes.max_rows = max_rows;
    tables.nodes.num_rows = max_rows;
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tables.edges.max_rows = max_rows;
    tables.edges.num_rows = max_rows;
    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tables.migrations.max_rows = max_rows;
    tables.migrations.num_rows = max_rows;
    ret_id = tsk_migration_table_add_row(&tables.migrations, 0, 0, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tables.sites.max_rows = max_rows;
    tables.sites.num_rows = max_rows;
    ret_id = tsk_site_table_add_row(&tables.sites, 0, 0, 0, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tables.mutations.max_rows = max_rows;
    tables.mutations.num_rows = max_rows;
    ret_id = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 0, 0, 0, 0, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tables.provenances.max_rows = max_rows;
    tables.provenances.num_rows = max_rows;
    ret_id = tsk_provenance_table_add_row(&tables.provenances, 0, 0, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tables.populations.max_rows = max_rows;
    tables.populations.num_rows = max_rows;
    ret_id = tsk_population_table_add_row(&tables.populations, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_TABLE_OVERFLOW);

    tsk_table_collection_free(&tables);
}

static void
test_column_overflow(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_size_t too_big = TSK_MAX_SIZE;
    double zero = 0;
    char zeros[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    tsk_id_t id_zeros[] = { 0, 0, 0, 0, 0, 0, 0, 0 };

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // location
    /* We can't trigger a column overflow with one element because the parameter
     * value is 32 bit */
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, &zero, 1, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    // Check normal overflow from additional length
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, too_big, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    // Check overflow from minimum increment
    ret = tsk_individual_table_set_max_location_length_increment(
        &tables.individuals, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 1, NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    // parents
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, id_zeros, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 1);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, too_big, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_individual_table_set_max_parents_length_increment(
        &tables.individuals, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    // metadata
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 2);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_individual_table_set_max_metadata_length_increment(
        &tables.individuals, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, 0, 0, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, 0, 0, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_node_table_set_max_metadata_length_increment(&tables.nodes, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 0, 0, 0, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 0, 0, 0, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 0, 0, 0, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_edge_table_set_max_metadata_length_increment(&tables.edges, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0, 0, 0, 0, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    ret_id = tsk_site_table_add_row(&tables.sites, 0, zeros, 1, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    // ancestral state
    ret_id = tsk_site_table_add_row(&tables.sites, 0, NULL, too_big, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_site_table_set_max_ancestral_state_length_increment(
        &tables.sites, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0, NULL, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    // metadata
    ret_id = tsk_site_table_add_row(&tables.sites, 0, NULL, 0, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_site_table_set_max_metadata_length_increment(&tables.sites, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0, NULL, 0, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 0, 0, zeros, 1, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    // derived state
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, 0, 0, NULL, too_big, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_mutation_table_set_max_derived_state_length_increment(
        &tables.mutations, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 0, 0, NULL, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    // metadata
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, 0, 0, NULL, 0, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_mutation_table_set_max_metadata_length_increment(
        &tables.mutations, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(&tables.mutations, 0, 0, 0, 0, NULL, 0, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    ret_id = tsk_provenance_table_add_row(&tables.provenances, zeros, 1, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0)
    // timestamp
    ret_id = tsk_provenance_table_add_row(&tables.provenances, NULL, too_big, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_provenance_table_set_max_timestamp_length_increment(
        &tables.provenances, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_provenance_table_add_row(&tables.provenances, NULL, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    // record
    ret_id = tsk_provenance_table_add_row(&tables.provenances, NULL, 0, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_provenance_table_set_max_record_length_increment(
        &tables.provenances, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_provenance_table_add_row(&tables.provenances, NULL, 0, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    ret_id = tsk_population_table_add_row(&tables.populations, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_population_table_set_max_metadata_length_increment(
        &tables.populations, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    ret_id = tsk_migration_table_add_row(&tables.migrations, 0, 0, 0, 0, 0, 0, zeros, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0, 0, 0, 0, 0, 0, NULL, too_big);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);
    ret = tsk_migration_table_set_max_metadata_length_increment(
        &tables.migrations, too_big);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(&tables.migrations, 0, 0, 0, 0, 0, 0, NULL, 1);
    CU_ASSERT_EQUAL_FATAL(ret_id, TSK_ERR_COLUMN_OVERFLOW);

    tsk_table_collection_free(&tables);
}

static void
test_table_collection_check_integrity_with_options(tsk_flags_t tc_options)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    const char *individuals = "1      0.25     -1\n"
                              "2      0.5,0.25 2\n"
                              "3      0.5,0.25 0\n";

    ret = tsk_table_collection_init(&tables, tc_options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    /* nodes */
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, INFINITY, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret_id);
    /* Not calling with TSK_CHECK_TREES so casting is safe */
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TIME_NONFINITE);

    ret = tsk_node_table_clear(&tables.nodes);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret_id);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_NO_CHECK_POPULATION_REFS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    ret = tsk_node_table_clear(&tables.nodes);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret_id);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    ret = tsk_node_table_clear(&tables.nodes);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 1.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* edges */
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, TSK_NULL, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_PARENT);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 2, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 1, TSK_NULL, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NULL_CHILD);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 1, 2, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, INFINITY, 1, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, -1.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_LEFT_LESS_ZERO);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.1, 1, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_RIGHT_GREATER_SEQ_LENGTH);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.5, 0.1, 1, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 0.5, 0, 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_NODE_TIME_ORDERING);

    ret = tsk_edge_table_clear(&tables.edges);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* sites */
    ret_id = tsk_site_table_add_row(&tables.sites, INFINITY, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SITE_POSITION);

    ret = tsk_site_table_clear(&tables.sites);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, -0.5, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SITE_POSITION);

    ret = tsk_site_table_clear(&tables.sites);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 1.5, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SITE_POSITION);

    ret = tsk_site_table_clear(&tables.sites);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.5, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.5, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, TSK_CHECK_SITE_DUPLICATES);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SITE_POSITION);

    ret = tsk_site_table_clear(&tables.sites);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.5, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.4, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, TSK_CHECK_SITE_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_SITES);

    ret = tsk_site_table_clear(&tables.sites);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.5, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.6, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    /* mutations */
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 2, 0, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 2, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    /* A mixture of known and unknown times on a site fails */
    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_TIME_HAS_BOTH_KNOWN_AND_UNKNOWN);

    /* But on different sites, passes */
    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 0, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(&tables.mutations, 0, 1, 2, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 1, 0, 1.0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_PARENT_EQUAL);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id
        = tsk_mutation_table_add_row(&tables.mutations, 0, 1, 1, 1.0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 1, TSK_NULL, 1.0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_PARENT_AFTER_CHILD);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 1, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, 0, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_PARENT_DIFFERENT_SITE);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 1, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_MUTATIONS);

    /* Unknown times pass */
    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Correctly ordered times pass */
    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 1, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 1, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Incorrectly ordered times fail */
    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 1, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_MUTATIONS);

    /* Putting incorrectly ordered times on diff sites passes */
    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 1, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 0, TSK_NULL, 2, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 0, TSK_NULL, 1, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, NAN, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TIME_NONFINITE);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, INFINITY, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TIME_NONFINITE);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, TSK_NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_TIME_YOUNGER_THAN_NODE);

    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, TSK_NULL, 1, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(&tables.mutations, 1, 1, 0, 2, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_TIME_OLDER_THAN_PARENT_MUTATION);
    ret = tsk_mutation_table_clear(&tables.mutations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MUTATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* migrations */
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, 0.5, 2, 0, 1, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, 0.5, 1, 2, 1, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, 0.5, 1, 0, 2, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, 0.5, 1, 0, 1, INFINITY, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_TIME_NONFINITE);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, 0.5, 1, 0, 1, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, 0.5, 1, 1, 0, 0.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_MIGRATION_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_MIGRATIONS);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, INFINITY, 1, 0, 1, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_GENOME_COORDS_NONFINITE);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, -0.3, 0.5, 1, 0, 1, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_LEFT_LESS_ZERO);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.0, 1.5, 1, 0, 1, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_RIGHT_GREATER_SEQ_LENGTH);

    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.6, 0.5, 1, 0, 1, 1.5, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_EDGE_INTERVAL);
    ret = tsk_migration_table_clear(&tables.migrations);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    parse_individuals(individuals, &tables.individuals);
    CU_ASSERT_EQUAL_FATAL(tables.individuals.num_rows, 3);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_INDIVIDUAL_ORDERING);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_INDIVIDUALS);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Check that an individual can't be its own parent */
    tables.individuals.parents[0] = 0;
    tables.individuals.parents[1] = 1;
    tables.individuals.parents[2] = 2;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_SELF_PARENT);

    tables.individuals.parents[0] = -2;
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);

    tsk_table_collection_free(&tables);
}

static void
test_table_collection_check_integrity_no_populations(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_id_t ret_num_trees;
    tsk_treeseq_t ts;
    tsk_table_collection_t tables;

    tsk_treeseq_from_text(&ts, 10, paper_ex_nodes, paper_ex_edges, NULL, paper_ex_sites,
        paper_ex_mutations, paper_ex_individuals, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Add in some bad population references and check that we can use
     * TSK_NO_CHECK_POPULATION_REFS with TSK_CHECK_TREES */
    tables.nodes.population[0] = 10;
    /* Not calling with TSK_CHECK_TREES so casting is safe */
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    ret_num_trees = tsk_table_collection_check_integrity(&tables, TSK_CHECK_TREES);
    CU_ASSERT_EQUAL_FATAL(ret_num_trees, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_NO_CHECK_POPULATION_REFS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_num_trees = tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_TREES | TSK_NO_CHECK_POPULATION_REFS);
    /* CHECK_TREES returns the number of trees */
    CU_ASSERT_EQUAL_FATAL(ret_num_trees, 3);
    tables.nodes.population[0] = TSK_NULL;

    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_migration_table_add_row(
        &tables.migrations, 0.4, 0.5, 1, 0, 1, 1.5, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret_id, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    ret_num_trees = tsk_table_collection_check_integrity(&tables, TSK_CHECK_TREES);
    CU_ASSERT_EQUAL_FATAL(ret_num_trees, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    ret = (int) tsk_table_collection_check_integrity(
        &tables, TSK_NO_CHECK_POPULATION_REFS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_num_trees = tsk_table_collection_check_integrity(
        &tables, TSK_CHECK_TREES | TSK_NO_CHECK_POPULATION_REFS);
    CU_ASSERT_EQUAL_FATAL(ret_num_trees, 3);

    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_table_collection_check_integrity(void)
{
    test_table_collection_check_integrity_with_options(0);
    test_table_collection_check_integrity_with_options(TSK_TC_NO_EDGE_METADATA);
}

static void
test_table_collection_subset_with_options(tsk_flags_t options)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_table_collection_t tables_copy;
    int k;
    tsk_id_t nodes[4];
    tsk_id_t zero_p[] = { 0 };
    tsk_id_t one_p[] = { 1 };

    ret = tsk_table_collection_init(&tables, options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    ret = tsk_table_collection_init(&tables_copy, options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // does not error on empty tables
    ret = tsk_table_collection_subset(&tables, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // four nodes from two diploids; the first is from pop 0
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 1.0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 2.0, TSK_NULL, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    // unused individual who is the parent of others
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, zero_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, one_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    // unused individual
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, one_p, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    // unused population
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 2, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.4, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    // unused site
    ret_id = tsk_site_table_add_row(&tables.sites, 0.5, "C", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, 0, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    // empty nodes should get empty tables
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(&tables_copy, NULL, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.nodes.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.individuals.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.populations.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.sites.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.mutations.num_rows, 0);

    // unless NO_CHANGE_POPULATIONS is provided
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(
        &tables_copy, NULL, 0, TSK_SUBSET_NO_CHANGE_POPULATIONS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.nodes.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.individuals.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.sites.num_rows, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.mutations.num_rows, 0);
    CU_ASSERT_FATAL(
        tsk_population_table_equals(&tables.populations, &tables_copy.populations, 0));

    // or KEEP_UNREFERENCED
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(
        &tables_copy, NULL, 0, TSK_SUBSET_KEEP_UNREFERENCED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.nodes.num_rows, 0);
    CU_ASSERT_FATAL(
        tsk_individual_table_equals(&tables.individuals, &tables_copy.individuals, 0));
    CU_ASSERT_EQUAL_FATAL(tables_copy.populations.num_rows, 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.mutations.num_rows, 0);
    CU_ASSERT_FATAL(tsk_site_table_equals(&tables.sites, &tables_copy.sites, 0));

    // or both
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(&tables_copy, NULL, 0,
        TSK_SUBSET_KEEP_UNREFERENCED | TSK_SUBSET_NO_CHANGE_POPULATIONS);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.nodes.num_rows, 0);
    CU_ASSERT_FATAL(
        tsk_individual_table_equals(&tables.individuals, &tables_copy.individuals, 0));
    CU_ASSERT_EQUAL_FATAL(tables_copy.mutations.num_rows, 0);
    CU_ASSERT_FATAL(
        tsk_population_table_equals(&tables.populations, &tables_copy.populations, 0));
    CU_ASSERT_FATAL(tsk_site_table_equals(&tables.sites, &tables_copy.sites, 0));

    // the identity transformation, since unused pops are at the end
    for (k = 0; k < 4; k++) {
        nodes[k] = k;
    }
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(
        &tables_copy, nodes, 4, TSK_SUBSET_KEEP_UNREFERENCED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &tables_copy, 0));

    // or, remove unused things:
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(&tables_copy, nodes, 4, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_node_table_equals(&tables.nodes, &tables_copy.nodes, 0));
    CU_ASSERT_EQUAL_FATAL(tables_copy.individuals.num_rows, 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.populations.num_rows, 1);
    CU_ASSERT_EQUAL_FATAL(tables_copy.sites.num_rows, 2);
    CU_ASSERT_FATAL(
        tsk_mutation_table_equals(&tables.mutations, &tables_copy.mutations, 0));

    // reverse twice should get back to the start, since unused pops are at the end
    for (k = 0; k < 4; k++) {
        nodes[k] = 3 - k;
    }
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT | options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(
        &tables_copy, nodes, 4, TSK_SUBSET_KEEP_UNREFERENCED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(
        &tables_copy, nodes, 4, TSK_SUBSET_KEEP_UNREFERENCED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = (int) tsk_table_collection_check_integrity(&tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &tables_copy, 0));

    tsk_table_collection_free(&tables_copy);
    tsk_table_collection_free(&tables);
}

static void
test_table_collection_subset(void)
{
    test_table_collection_subset_with_options(0);
    test_table_collection_subset_with_options(TSK_TC_NO_EDGE_METADATA);
}

static void
test_table_collection_subset_unsorted(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_table_collection_t tables_copy;
    int k;
    tsk_id_t nodes[3];
    tsk_id_t one_p[] = { 1 };

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    ret = tsk_table_collection_init(&tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // these tables are a big mess
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.5, TSK_NULL, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(&tables.nodes, 0, 1.0, TSK_NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, one_p, 1, NULL, 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 0.5, 2, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.5, 1.0, 2, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.4, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, 2, TSK_UNKNOWN_TIME, "B", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    // but still, this should leave them unchanged
    for (k = 0; k < 3; k++) {
        nodes[k] = k;
    }
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(
        &tables_copy, nodes, 3, TSK_SUBSET_KEEP_UNREFERENCED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &tables_copy, 0));

    tsk_table_collection_free(&tables_copy);
    tsk_table_collection_free(&tables);
}

static void
test_table_collection_subset_errors(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_table_collection_t tables_copy;
    tsk_id_t nodes[4] = { 0, 1, 2, 3 };

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    ret = tsk_table_collection_init(&tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // four nodes from two diploids; the first is from pop 0
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 1.0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 2.0, TSK_NULL, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Migrations are not supported */
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_migration_table_add_row(&tables_copy.migrations, 0, 1, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.migrations.num_rows, 1);
    ret = tsk_table_collection_subset(&tables_copy, nodes, 4, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATIONS_NOT_SUPPORTED);

    // test out of bounds nodes
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    nodes[0] = -1;
    ret = tsk_table_collection_subset(&tables_copy, nodes, 4, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);
    nodes[0] = 6;
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_subset(&tables_copy, nodes, 4, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_NODE_OUT_OF_BOUNDS);

    // check integrity
    nodes[0] = 0;
    nodes[1] = 1;
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_truncate(&tables_copy.nodes, 3);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_node_table_add_row(
        &tables_copy.nodes, TSK_NODE_IS_SAMPLE, 0.0, -2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = tsk_table_collection_subset(&tables_copy, nodes, 4, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    tsk_table_collection_free(&tables);
    tsk_table_collection_free(&tables_copy);
}

static void
test_table_collection_union(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_table_collection_t tables_empty;
    tsk_table_collection_t tables_copy;
    tsk_id_t node_mapping[3];
    tsk_id_t parents[2] = { -1, -1 };
    char example_metadata[100] = "An example of metadata with unicode 🎄🌳🌴🌲🎋";
    tsk_size_t example_metadata_length = (tsk_size_t) strlen(example_metadata);

    tsk_memset(node_mapping, 0xff, sizeof(node_mapping));

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    ret = tsk_table_collection_init(&tables_empty, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables_empty.sequence_length = 1;
    ret = tsk_table_collection_init(&tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // does not error on empty tables
    ret = tsk_table_collection_union(&tables, &tables_empty, node_mapping, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // does not error on empty tables but that differ on top level metadata
    ret = tsk_table_collection_set_metadata(
        &tables, example_metadata, example_metadata_length);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_union(&tables, &tables_empty, node_mapping, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // three nodes, two pop, three ind, two edge, two site, two mut
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 1, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.5, 1, 2, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, parents, 2, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    parents[0] = 0;
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, parents, 2, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    parents[1] = 1;
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, parents, 2, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 2, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.4, "T", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 1, 1, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_sort(&tables, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // union with empty should not change
    // other is empty
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_union(
        &tables_copy, &tables_empty, node_mapping, TSK_UNION_NO_CHECK_SHARED);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &tables_copy, 0));
    // self is empty
    ret = tsk_table_collection_clear(&tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_union(
        &tables_copy, &tables, node_mapping, TSK_UNION_NO_CHECK_SHARED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &tables_copy, 0));

    // union all shared nodes + subset original nodes = original table
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_union(
        &tables_copy, &tables, node_mapping, TSK_UNION_NO_CHECK_SHARED);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    node_mapping[0] = 0;
    node_mapping[1] = 1;
    node_mapping[2] = 2;
    ret = tsk_table_collection_subset(&tables_copy, node_mapping, 3, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tables, &tables_copy, 0));

    // union with one shared node
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    node_mapping[0] = TSK_NULL;
    node_mapping[1] = TSK_NULL;
    node_mapping[2] = 2;
    ret = tsk_table_collection_union(&tables_copy, &tables, node_mapping, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(
        tables_copy.populations.num_rows, tables.populations.num_rows + 2);
    CU_ASSERT_EQUAL_FATAL(
        tables_copy.individuals.num_rows, tables.individuals.num_rows + 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.nodes.num_rows, tables.nodes.num_rows + 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.edges.num_rows, tables.edges.num_rows + 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.sites.num_rows, tables.sites.num_rows);
    CU_ASSERT_EQUAL_FATAL(tables_copy.mutations.num_rows, tables.mutations.num_rows + 2);

    // union with one shared node, but no add pop
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    node_mapping[0] = TSK_NULL;
    node_mapping[1] = TSK_NULL;
    node_mapping[2] = 2;
    ret = tsk_table_collection_union(
        &tables_copy, &tables, node_mapping, TSK_UNION_NO_ADD_POP);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.populations.num_rows, tables.populations.num_rows);
    CU_ASSERT_EQUAL_FATAL(
        tables_copy.individuals.num_rows, tables.individuals.num_rows + 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.nodes.num_rows, tables.nodes.num_rows + 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.edges.num_rows, tables.edges.num_rows + 2);
    CU_ASSERT_EQUAL_FATAL(tables_copy.sites.num_rows, tables.sites.num_rows);
    CU_ASSERT_EQUAL_FATAL(tables_copy.mutations.num_rows, tables.mutations.num_rows + 2);

    tsk_table_collection_free(&tables_copy);
    tsk_table_collection_free(&tables_empty);
    tsk_table_collection_free(&tables);
}

static void
test_table_collection_union_middle_merge(void)
{
    /* Test ability to have non-shared history both above and below the
     * shared bits. The full genealogy, in `tu`, is:
     *  3   4
     *   \ /
     *    2
     *   / \
     *  0   1
     * and the left lineage is in `ta` and right in `tb` */
    int ret;
    tsk_id_t ret_id;
    tsk_id_t node_mapping[] = { TSK_NULL, 1, TSK_NULL };
    tsk_id_t node_order[] = { 0, 3, 1, 2, 4 };
    tsk_table_collection_t ta, tb, tu;
    ret = tsk_table_collection_init(&ta, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ta.sequence_length = 1;
    ret = tsk_table_collection_init(&tb, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tb.sequence_length = 1;
    ret = tsk_table_collection_init(&tu, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tu.sequence_length = 1;

    ret_id = tsk_node_table_add_row(
        &tu.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0); // node u0
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &ta.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0); // node a0 = u0
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tu.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0); // node u1
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tb.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, NULL, 0); // node b0 = u1
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tu.nodes, 0, 1, TSK_NULL, TSK_NULL, NULL, 0); // node u2
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tu.edges, 0, 1, 2, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tu.edges, 0, 1, 2, 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &ta.nodes, 0, 1, TSK_NULL, TSK_NULL, NULL, 0); // node a1 = u2
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&ta.edges, 0, 1, 1, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tb.nodes, 0, 1, TSK_NULL, TSK_NULL, NULL, 0); // node b1 = u2
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tb.edges, 0, 1, 1, 0, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tu.nodes, 0, 2, TSK_NULL, TSK_NULL, NULL, 0); // node u3
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tu.edges, 0, 0.5, 3, 2, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &ta.nodes, 0, 2, TSK_NULL, TSK_NULL, NULL, 0); // node a2 = u3
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&ta.edges, 0, 0.5, 2, 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tu.nodes, 0, 2, TSK_NULL, TSK_NULL, NULL, 0); // node u4
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tu.edges, 0.5, 1, 4, 2, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_node_table_add_row(
        &tb.nodes, 0, 2, TSK_NULL, TSK_NULL, NULL, 0); // node b2 = u4
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tb.edges, 0.5, 1, 2, 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);

    ret_id = tsk_site_table_add_row(&ta.sites, 0.25, "A", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&ta.sites, 0.75, "X", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tb.sites, 0.25, "A", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tb.sites, 0.75, "X", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tu.sites, 0.25, "A", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tu.sites, 0.75, "X", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);

    ret_id = tsk_mutation_table_add_row(
        &tu.mutations, 0, 3, TSK_NULL, 3.5, "B", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &ta.mutations, 0, 2, TSK_NULL, 3.5, "B", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tu.mutations, 0, 2, TSK_NULL, 1.5, "D", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &ta.mutations, 0, 1, TSK_NULL, 1.5, "D", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tb.mutations, 0, 1, TSK_NULL, 1.5, "D", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tu.mutations, 0, 2, TSK_NULL, 1.2, "E", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &ta.mutations, 0, 1, TSK_NULL, 1.2, "E", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tb.mutations, 0, 1, TSK_NULL, 1.2, "E", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tu.mutations, 0, 0, TSK_NULL, 0.5, "C", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &ta.mutations, 0, 0, TSK_NULL, 0.5, "C", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tu.mutations, 1, 4, TSK_NULL, 2.4, "Y", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tb.mutations, 1, 2, TSK_NULL, 2.4, "Y", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tu.mutations, 1, 1, TSK_NULL, 0.4, "Z", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tb.mutations, 1, 0, TSK_NULL, 0.4, "Z", 1, NULL, 0);
    CU_ASSERT(ret_id >= 0);

    ret = tsk_table_collection_build_index(&ta, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_compute_mutation_parents(&ta, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_build_index(&tb, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_compute_mutation_parents(&tb, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_build_index(&tu, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_compute_mutation_parents(&tu, 0);
    CU_ASSERT_EQUAL(ret, 0);

    ret = tsk_table_collection_union(&ta, &tb, node_mapping, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_table_collection_subset(&ta, node_order, 5, 0);
    CU_ASSERT_EQUAL(ret, 0);
    ret = tsk_provenance_table_clear(&ta.provenances);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_FATAL(tsk_table_collection_equals(&tu, &ta, 0));

    tsk_table_collection_free(&ta);
    tsk_table_collection_free(&tb);
    tsk_table_collection_free(&tu);
}

static void
test_table_collection_union_errors(void)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    tsk_table_collection_t tables_copy;
    tsk_id_t node_mapping[] = { 0, 1 };

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;
    ret = tsk_table_collection_init(&tables_copy, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    // two nodes, two pop, two ind, one edge, one site, one mut
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.5, 1, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    // trigger diff histories error
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id = tsk_mutation_table_add_row(
        &tables_copy.mutations, 0, 1, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = tsk_table_collection_union(&tables_copy, &tables, node_mapping, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNION_DIFF_HISTORIES);

    // Migrations are not supported
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_migration_table_add_row(&tables_copy.migrations, 0, 1, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(tables_copy.migrations.num_rows, 1);
    ret = tsk_table_collection_union(
        &tables_copy, &tables, node_mapping, TSK_UNION_NO_CHECK_SHARED);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATIONS_NOT_SUPPORTED);

    // test out of bounds node_mapping
    node_mapping[0] = -4;
    node_mapping[1] = 6;
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_union(&tables_copy, &tables, node_mapping, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNION_BAD_MAP);

    // check integrity
    node_mapping[0] = 0;
    node_mapping[1] = 1;
    ret_id = tsk_node_table_add_row(
        &tables_copy.nodes, TSK_NODE_IS_SAMPLE, 0.0, -2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = tsk_table_collection_union(&tables_copy, &tables, node_mapping, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    ret = tsk_table_collection_copy(&tables, &tables_copy, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, -2, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret = tsk_table_collection_union(&tables, &tables_copy, node_mapping, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);

    tsk_table_collection_free(&tables_copy);
    tsk_table_collection_free(&tables);
}

static void
test_table_collection_clear_with_options(tsk_flags_t options)
{
    int ret;
    tsk_id_t ret_id;
    tsk_table_collection_t tables;
    bool clear_provenance = !!(options & TSK_CLEAR_PROVENANCE);
    bool clear_metadata_schemas = !!(options & TSK_CLEAR_METADATA_SCHEMAS);
    bool clear_ts_metadata = !!(options & TSK_CLEAR_TS_METADATA_AND_SCHEMA);
    tsk_bookmark_t num_rows;
    tsk_bookmark_t expected_rows = { .provenances = clear_provenance ? 0 : 1 };
    tsk_size_t expected_len = clear_metadata_schemas ? 0 : 4;
    tsk_size_t expected_len_ts = clear_ts_metadata ? 0 : 4;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id
        = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0.5, 1, 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_individual_table_add_row(
        &tables.individuals, 0, NULL, 0, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_population_table_add_row(&tables.populations, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_edge_table_add_row(&tables.edges, 0.0, 1.0, 1, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_site_table_add_row(&tables.sites, 0.2, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_mutation_table_add_row(
        &tables.mutations, 0, 0, TSK_NULL, TSK_UNKNOWN_TIME, NULL, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);
    ret_id = tsk_migration_table_add_row(&tables.migrations, 0, 1, 0, 0, 0, 0, NULL, 0);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_table_collection_build_index(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_individual_table_set_metadata_schema(&tables.individuals, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_node_table_set_metadata_schema(&tables.nodes, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_edge_table_set_metadata_schema(&tables.edges, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_migration_table_set_metadata_schema(&tables.migrations, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_site_table_set_metadata_schema(&tables.sites, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_mutation_table_set_metadata_schema(&tables.mutations, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_population_table_set_metadata_schema(&tables.populations, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_set_time_units(&tables, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata(&tables, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_set_metadata_schema(&tables, "test", 4);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret_id = tsk_provenance_table_add_row(&tables.provenances, "today", 5, "test", 4);
    CU_ASSERT_FATAL(ret_id >= 0);

    ret = tsk_table_collection_clear(&tables, options);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_table_collection_record_num_rows(&tables, &num_rows);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(num_rows.individuals, expected_rows.individuals);
    CU_ASSERT_EQUAL(num_rows.nodes, expected_rows.nodes);
    CU_ASSERT_EQUAL(num_rows.edges, expected_rows.edges);
    CU_ASSERT_EQUAL(num_rows.migrations, expected_rows.migrations);
    CU_ASSERT_EQUAL(num_rows.sites, expected_rows.sites);
    CU_ASSERT_EQUAL(num_rows.mutations, expected_rows.mutations);
    CU_ASSERT_EQUAL(num_rows.populations, expected_rows.populations);
    CU_ASSERT_EQUAL(num_rows.provenances, expected_rows.provenances);

    CU_ASSERT_FALSE(tsk_table_collection_has_index(&tables, 0));

    CU_ASSERT_EQUAL(tables.individuals.metadata_schema_length, expected_len);
    CU_ASSERT_EQUAL(tables.nodes.metadata_schema_length, expected_len);
    CU_ASSERT_EQUAL(tables.edges.metadata_schema_length, expected_len);
    CU_ASSERT_EQUAL(tables.migrations.metadata_schema_length, expected_len);
    CU_ASSERT_EQUAL(tables.sites.metadata_schema_length, expected_len);
    CU_ASSERT_EQUAL(tables.mutations.metadata_schema_length, expected_len);
    CU_ASSERT_EQUAL(tables.populations.metadata_schema_length, expected_len);
    CU_ASSERT_EQUAL(tables.metadata_schema_length, expected_len_ts);
    CU_ASSERT_EQUAL(tables.metadata_length, expected_len_ts);
    CU_ASSERT_EQUAL(tables.time_units_length, 4);

    tsk_table_collection_free(&tables);
}

static void
test_table_collection_clear(void)
{
    test_table_collection_clear_with_options(0);
    test_table_collection_clear_with_options(TSK_CLEAR_PROVENANCE);
    test_table_collection_clear_with_options(TSK_CLEAR_METADATA_SCHEMAS);
    test_table_collection_clear_with_options(TSK_CLEAR_TS_METADATA_AND_SCHEMA);
    test_table_collection_clear_with_options(
        TSK_CLEAR_PROVENANCE | TSK_CLEAR_METADATA_SCHEMAS);
    test_table_collection_clear_with_options(
        TSK_CLEAR_PROVENANCE | TSK_CLEAR_TS_METADATA_AND_SCHEMA);
    test_table_collection_clear_with_options(
        TSK_CLEAR_METADATA_SCHEMAS | TSK_CLEAR_TS_METADATA_AND_SCHEMA);
    test_table_collection_clear_with_options(TSK_CLEAR_PROVENANCE
                                             | TSK_CLEAR_METADATA_SCHEMAS
                                             | TSK_CLEAR_TS_METADATA_AND_SCHEMA);
}

static void
test_table_collection_takeset_indexes(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_table_collection_t t1, t2;
    tsk_id_t *ins;
    tsk_id_t *rem;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges, NULL, NULL,
        NULL, NULL, NULL, 0);
    ret = tsk_treeseq_copy_tables(&ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ins = tsk_malloc(t1.edges.num_rows * sizeof(*ins));
    CU_ASSERT_FATAL(ins != NULL);
    rem = tsk_malloc(t1.edges.num_rows * sizeof(*rem));
    CU_ASSERT_FATAL(rem != NULL);
    memcpy(ins, t1.indexes.edge_insertion_order,
        (size_t)(t1.edges.num_rows * sizeof(*ins)));
    memcpy(
        rem, t1.indexes.edge_removal_order, (size_t)(t1.edges.num_rows * sizeof(*rem)));

    ret = tsk_table_collection_copy(&t1, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_drop_index(&t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_takeset_indexes(&t2, ins, rem);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(
        tsk_memcmp(t1.indexes.edge_insertion_order, t2.indexes.edge_insertion_order,
            t1.edges.num_rows * sizeof(*ins)),
        0);
    CU_ASSERT_EQUAL(tsk_memcmp(t1.indexes.edge_removal_order,
                        t2.indexes.edge_removal_order, t1.edges.num_rows * sizeof(*rem)),
        0);

    ret = tsk_table_collection_takeset_indexes(&t2, ins, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_table_collection_takeset_indexes(&t2, NULL, rem);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
    tsk_treeseq_free(&ts);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_node_table", test_node_table },
        { "test_node_table_update_row", test_node_table_update_row },
        { "test_node_table_takeset", test_node_table_takeset },
        { "test_edge_table", test_edge_table },
        { "test_edge_table_update_row", test_edge_table_update_row },
        { "test_edge_table_update_row_no_metadata",
            test_edge_table_update_row_no_metadata },
        { "test_edge_table_takeset", test_edge_table_takeset },
        { "test_edge_table_copy_semantics", test_edge_table_copy_semantics },
        { "test_edge_table_squash", test_edge_table_squash },
        { "test_edge_table_squash_multiple_parents",
            test_edge_table_squash_multiple_parents },
        { "test_edge_table_squash_empty", test_edge_table_squash_empty },
        { "test_edge_table_squash_single_edge", test_edge_table_squash_single_edge },
        { "test_edge_table_squash_bad_intervals", test_edge_table_squash_bad_intervals },
        { "test_edge_table_squash_metadata", test_edge_table_squash_metadata },
        { "test_site_table", test_site_table },
        { "test_site_table_update_row", test_site_table_update_row },
        { "test_site_table_takeset", test_site_table_takeset },
        { "test_mutation_table", test_mutation_table },
        { "test_mutation_table_update_row", test_mutation_table_update_row },
        { "test_mutation_table_takeset", test_mutation_table_takeset },
        { "test_migration_table", test_migration_table },
        { "test_migration_table_update_row", test_migration_table_update_row },
        { "test_migration_table_takeset", test_migration_table_takeset },
        { "test_individual_table", test_individual_table },
        { "test_individual_table_takeset", test_individual_table_takeset },
        { "test_individual_table_update_row", test_individual_table_update_row },
        { "test_population_table", test_population_table },
        { "test_population_table_update_row", test_population_table_update_row },
        { "test_population_table_takeset", test_population_table_takeset },
        { "test_provenance_table", test_provenance_table },
        { "test_provenance_table_update_row", test_provenance_table_update_row },
        { "test_provenance_table_takeset", test_provenance_table_takeset },
        { "test_table_size_increments", test_table_size_increments },
        { "test_table_expansion", test_table_expansion },
        { "test_ragged_expansion", test_ragged_expansion },
        { "test_table_collection_equals_options", test_table_collection_equals_options },
        { "test_table_collection_simplify_errors",
            test_table_collection_simplify_errors },
        { "test_table_collection_time_units", test_table_collection_time_units },
        { "test_table_collection_reference_sequence",
            test_table_collection_reference_sequence },
        { "test_table_collection_has_reference_sequence",
            test_table_collection_has_reference_sequence },
        { "test_table_collection_metadata", test_table_collection_metadata },
        { "test_reference_sequence_state_machine",
            test_reference_sequence_state_machine },
        { "test_reference_sequence_take", test_reference_sequence_take },
        { "test_reference_sequence", test_reference_sequence },

        { "test_simplify_tables_drops_indexes", test_simplify_tables_drops_indexes },
        { "test_simplify_empty_tables", test_simplify_empty_tables },
        { "test_simplify_metadata", test_simplify_metadata },
        { "test_link_ancestors_no_edges", test_link_ancestors_no_edges },
        { "test_link_ancestors_input_errors", test_link_ancestors_input_errors },
        { "test_link_ancestors_single_tree", test_link_ancestors_single_tree },
        { "test_link_ancestors_paper", test_link_ancestors_paper },
        { "test_link_ancestors_samples_and_ancestors_overlap",
            test_link_ancestors_samples_and_ancestors_overlap },
        { "test_link_ancestors_multiple_to_single_tree",
            test_link_ancestors_multiple_to_single_tree },
        { "test_ibd_segments_debug", test_ibd_segments_debug },
        { "test_ibd_segments_caterpillar_tree", test_ibd_segments_caterpillar_tree },
        { "test_ibd_segments_single_tree", test_ibd_segments_single_tree },
        { "test_ibd_segments_single_tree_options",
            test_ibd_segments_single_tree_options },
        { "test_ibd_segments_multiple_trees", test_ibd_segments_multiple_trees },
        { "test_ibd_segments_empty_result", test_ibd_segments_empty_result },
        { "test_ibd_segments_min_span_max_time", test_ibd_segments_min_span_max_time },
        { "test_ibd_segments_single_tree_between",
            test_ibd_segments_single_tree_between },
        { "test_ibd_segments_samples_are_descendants",
            test_ibd_segments_samples_are_descendants },
        { "test_ibd_segments_multiple_ibd_paths", test_ibd_segments_multiple_ibd_paths },
        { "test_ibd_segments_odd_topologies", test_ibd_segments_odd_topologies },
        { "test_ibd_segments_errors", test_ibd_segments_errors },
        { "test_sorter_interface", test_sorter_interface },
        { "test_sort_tables_canonical_errors", test_sort_tables_canonical_errors },
        { "test_sort_tables_canonical", test_sort_tables_canonical },
        { "test_sort_tables_drops_indexes", test_sort_tables_drops_indexes },
        { "test_sort_tables_edge_metadata", test_sort_tables_edge_metadata },
        { "test_sort_tables_errors", test_sort_tables_errors },
        { "test_sort_tables_individuals", test_sort_tables_individuals },
        { "test_sort_tables_mutation_times", test_sort_tables_mutation_times },
        { "test_sort_tables_migrations", test_sort_tables_migrations },
        { "test_sort_tables_no_edge_metadata", test_sort_tables_no_edge_metadata },
        { "test_sort_tables_offsets", test_sort_tables_offsets },
        { "test_edge_update_invalidates_index", test_edge_update_invalidates_index },
        { "test_copy_table_collection", test_copy_table_collection },
        { "test_dump_unindexed", test_dump_unindexed },
        { "test_dump_load_empty", test_dump_load_empty },
        { "test_dump_load_unsorted", test_dump_load_unsorted },
        { "test_dump_load_metadata_schema", test_dump_load_metadata_schema },
        { "test_dump_fail_no_file", test_dump_fail_no_file },
        { "test_load_reindex", test_load_reindex },
        { "test_table_overflow", test_table_overflow },
        { "test_column_overflow", test_column_overflow },
        { "test_table_collection_check_integrity",
            test_table_collection_check_integrity },
        { "test_table_collection_check_integrity_no_populations",
            test_table_collection_check_integrity_no_populations },
        { "test_table_collection_subset", test_table_collection_subset },
        { "test_table_collection_subset_unsorted",
            test_table_collection_subset_unsorted },
        { "test_table_collection_subset_errors", test_table_collection_subset_errors },
        { "test_table_collection_union", test_table_collection_union },
        { "test_table_collection_union_middle_merge",
            test_table_collection_union_middle_merge },
        { "test_table_collection_union_errors", test_table_collection_union_errors },
        { "test_table_collection_clear", test_table_collection_clear },
        { "test_table_collection_takeset_indexes",
            test_table_collection_takeset_indexes },
        { NULL, NULL },
    };

    return test_main(tests, argc, argv);
}
