#include "testlib.h"
#include "tsk_tables.h"

#include <unistd.h>
#include <stdlib.h>

typedef struct {
    const char *name;
    void *array;
    tsk_tbl_size_t len;
    int type;
} write_table_col_t;

static void
write_table_cols(kastore_t *store, write_table_col_t *write_cols, size_t num_cols)
{
    size_t j;
    int ret;

    for (j = 0; j < num_cols; j++) {
        ret = kastore_puts(store, write_cols[j].name, write_cols[j].array,
                write_cols[j].len, write_cols[j].type, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
test_format_data_load_errors(void)
{
    size_t uuid_size = 36;
    char uuid[uuid_size];
    char format_name[TSK_FILE_FORMAT_NAME_LENGTH];
    double L[2];
    uint32_t version[2] = {
        TSK_FILE_FORMAT_VERSION_MAJOR, TSK_FILE_FORMAT_VERSION_MINOR};
    write_table_col_t write_cols[] = {
        {"format/name", (void *) format_name, sizeof(format_name), KAS_INT8},
        {"format/version", (void *) version, 2, KAS_UINT32},
        {"sequence_length", (void *) L, 1, KAS_FLOAT64},
        {"uuid", (void *) uuid, (tsk_tbl_size_t) uuid_size, KAS_INT8},
    };
    tsk_tbl_collection_t tables;
    kastore_t store;
    size_t j;
    int ret;

    L[0] = 1;
    L[1] = 0;
    memcpy(format_name, TSK_FILE_FORMAT_NAME, sizeof(format_name));
    /* Note: this will fail if we ever start parsing the form of the UUID */
    memset(uuid, 0, uuid_size);

    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    /* We've only defined the format headers, so we should fail immediately
     * after with key not found */
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_KEY_NOT_FOUND);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Version too old */
    version[0] = TSK_FILE_FORMAT_VERSION_MAJOR - 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_VERSION_TOO_OLD);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Version too new */
    version[0] = TSK_FILE_FORMAT_VERSION_MAJOR + 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_VERSION_TOO_NEW);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    version[0] = TSK_FILE_FORMAT_VERSION_MAJOR;

    /* Bad version length */
    write_cols[1].len = 0;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[1].len = 2;

    /* Bad format name length */
    write_cols[0].len = 0;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[0].len = TSK_FILE_FORMAT_NAME_LENGTH;

    /* Bad format name */
    format_name[0] = 'X';
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    format_name[0] = 't';

    /* Bad type for sequence length. */
    write_cols[2].type = KAS_FLOAT32;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_TYPE_MISMATCH);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[2].type = KAS_FLOAT64;

    /* Bad length for sequence length. */
    write_cols[2].len = 2;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[2].len = 1;

    /* Bad value for sequence length. */
    L[0] = -1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SEQUENCE_LENGTH);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    L[0] = 1;

    /* Wrong length for uuid */
    write_cols[3].len = 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[3].len = (tsk_tbl_size_t) uuid_size;

    /* Missing keys */
    for (j = 0; j < sizeof(write_cols) / sizeof(*write_cols) - 1; j++) {
        ret = kastore_open(&store, _tmp_file_name, "w", 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        write_table_cols(&store, write_cols, j);
        ret = kastore_close(&store);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
        CU_ASSERT_TRUE(tsk_is_kas_error(ret));
        CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_KEY_NOT_FOUND);
        CU_ASSERT_STRING_EQUAL(tsk_strerror(ret), kas_strerror(KAS_ERR_KEY_NOT_FOUND));
        ret = tsk_tbl_collection_free(&tables);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
test_dump_unindexed(void)
{
    tsk_tbl_collection_t tables, loaded;
    int ret;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    tables.sequence_length = 1;
    parse_nodes(single_tree_ex_nodes, tables.nodes);
    CU_ASSERT_EQUAL_FATAL(tables.nodes->num_rows, 7);
    parse_edges(single_tree_ex_edges, tables.edges);
    CU_ASSERT_EQUAL_FATAL(tables.edges->num_rows, 6);
    CU_ASSERT_FALSE(tsk_tbl_collection_is_indexed(&tables));
    ret = tsk_tbl_collection_dump(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tbl_collection_is_indexed(&tables));

    ret = tsk_tbl_collection_load(&loaded, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_tbl_collection_is_indexed(&loaded));
    CU_ASSERT_TRUE(tsk_node_tbl_equals(tables.nodes, loaded.nodes));
    CU_ASSERT_TRUE(tsk_edge_tbl_equals(tables.edges, loaded.edges));

    tsk_tbl_collection_free(&loaded);
    tsk_tbl_collection_free(&tables);
}

static void
test_tbl_collection_load_errors(void)
{
    tsk_tbl_collection_t tables;
    int ret;
    const char *str;

    ret = tsk_tbl_collection_load(&tables, "/", 0);
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);

    tsk_tbl_collection_free(&tables);
}

static void
test_tbl_collection_dump_errors(void)
{
    tsk_tbl_collection_t tables;
    int ret;
    const char *str;

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_dump(&tables, "/", 0);
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);

    tsk_tbl_collection_free(&tables);
}
static void
test_tbl_collection_simplify_errors(void)
{
    int ret;
    tsk_tbl_collection_t tables;
    tsk_id_t samples[] = {0, 1};

    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1;

    ret = tsk_site_tbl_add_row(tables.sites, 0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_site_tbl_add_row(tables.sites, 0, "A", 1, NULL, 0);
    CU_ASSERT_FATAL(ret >= 0);
    ret = tsk_tbl_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_DUPLICATE_SITE_POSITION);

    /* Out of order positions */
    tables.sites->position[0] = 0.5;
    ret = tsk_tbl_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_UNSORTED_SITES);

    /* Position out of bounds */
    tables.sites->position[0] = 1.5;
    ret = tsk_tbl_collection_simplify(&tables, samples, 0, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SITE_POSITION);

    /* TODO More tests for this: see
     * https://github.com/tskit-dev/msprime/issues/517 */

    tsk_tbl_collection_free(&tables);
}

static void
test_load_tsk_node_tbl_errors(void)
{
    char format_name[TSK_FILE_FORMAT_NAME_LENGTH];
    tsk_tbl_size_t uuid_size = 36;
    char uuid[uuid_size];
    double L = 1;
    double time = 0;
    double flags = 0;
    int32_t population = 0;
    int32_t individual = 0;
    int8_t metadata = 0;
    uint32_t metadata_offset[] = {0, 1};
    uint32_t version[2] = {
        TSK_FILE_FORMAT_VERSION_MAJOR, TSK_FILE_FORMAT_VERSION_MINOR};
    write_table_col_t write_cols[] = {
        {"nodes/time", (void *) &time, 1, KAS_FLOAT64},
        {"nodes/flags", (void *) &flags, 1, KAS_UINT32},
        {"nodes/population", (void *) &population, 1, KAS_INT32},
        {"nodes/individual", (void *) &individual, 1, KAS_INT32},
        {"nodes/metadata", (void *) &metadata, 1, KAS_UINT8},
        {"nodes/metadata_offset", (void *) metadata_offset, 2, KAS_UINT32},
        {"format/name", (void *) format_name, sizeof(format_name), KAS_INT8},
        {"format/version", (void *) version, 2, KAS_UINT32},
        {"uuid", (void *) uuid, uuid_size, KAS_INT8},
        {"sequence_length", (void *) &L, 1, KAS_FLOAT64},
    };
    tsk_tbl_collection_t tables;
    kastore_t store;
    int ret;

    memcpy(format_name, TSK_FILE_FORMAT_NAME, sizeof(format_name));
    /* Note: this will fail if we ever start parsing the form of the UUID */
    memset(uuid, 0, uuid_size);

    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    /* We've only defined the format headers and nodes, so we should fail immediately
     * after with key not found */
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_KEY_NOT_FOUND);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Wrong type for time */
    write_cols[0].type = KAS_INT64;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[0].type = KAS_FLOAT64;

    /* Wrong length for flags */
    write_cols[1].len = 0;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[1].len = 1;

    /* Missing key */
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols) - 1);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_KEY_NOT_FOUND);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Wrong length for metadata offset */
    write_cols[5].len = 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_tbl_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_tbl_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[5].len = 2;

}

static void
test_node_table(void)
{
    int ret;
    tsk_node_tbl_t table;
    tsk_node_t node;
    uint32_t num_rows = 100;
    uint32_t j;
    uint32_t *flags;
    tsk_id_t *population;
    double *time;
    tsk_id_t *individual;
    char *metadata;
    uint32_t *metadata_offset;
    const char *test_metadata = "test";
    size_t test_metadata_length = 4;
    char metadata_copy[test_metadata_length + 1];

    metadata_copy[test_metadata_length] = '\0';
    ret = tsk_node_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_node_tbl_set_max_rows_increment(&table, 1);
    tsk_node_tbl_set_max_metadata_length_increment(&table, 1);
    tsk_node_tbl_print_state(&table, _devnull);
    tsk_node_tbl_dump_text(&table, _devnull);

    for (j = 0; j < num_rows; j++) {
        ret = tsk_node_tbl_add_row(&table, j, j, (tsk_id_t) j, (tsk_id_t) j,
                test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        CU_ASSERT_EQUAL(table.flags[j], j);
        CU_ASSERT_EQUAL(table.time[j], j);
        CU_ASSERT_EQUAL(table.population[j], j);
        CU_ASSERT_EQUAL(table.individual[j], j);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
        CU_ASSERT_EQUAL(table.metadata_length, (j + 1) * test_metadata_length);
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], table.metadata_length);
        /* check the metadata */
        memcpy(metadata_copy, table.metadata + table.metadata_offset[j], test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(metadata_copy, test_metadata, test_metadata_length);
        ret = tsk_node_tbl_get_row(&table, j, &node);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(node.id, j);
        CU_ASSERT_EQUAL(node.flags, j);
        CU_ASSERT_EQUAL(node.time, j);
        CU_ASSERT_EQUAL(node.population, j);
        CU_ASSERT_EQUAL(node.individual, j);
        CU_ASSERT_EQUAL(node.metadata_length, test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(node.metadata, test_metadata, test_metadata_length);
    }
    CU_ASSERT_EQUAL(tsk_node_tbl_get_row(&table, num_rows, &node),
            TSK_ERR_NODE_OUT_OF_BOUNDS);
    tsk_node_tbl_print_state(&table, _devnull);
    tsk_node_tbl_dump_text(&table, _devnull);

    tsk_node_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    num_rows *= 2;
    flags = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(flags != NULL);
    memset(flags, 1, num_rows * sizeof(uint32_t));
    population = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(population != NULL);
    memset(population, 2, num_rows * sizeof(uint32_t));
    time = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    memset(time, 0, num_rows * sizeof(double));
    individual = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(individual != NULL);
    memset(individual, 3, num_rows * sizeof(uint32_t));
    metadata = malloc(num_rows * sizeof(char));
    memset(metadata, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    for (j = 0; j < num_rows + 1; j++) {
        metadata_offset[j] = j;
    }
    ret = tsk_node_tbl_set_columns(&table, num_rows, flags, time, population,
            individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.population, population, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.individual, individual, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    tsk_node_tbl_print_state(&table, _devnull);
    tsk_node_tbl_dump_text(&table, _devnull);

    /* Append another num_rows onto the end */
    ret = tsk_node_tbl_append_columns(&table, num_rows, flags, time, population,
            individual, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.flags + num_rows, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.population, population, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.population + num_rows, population,
                num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.individual, individual, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.individual + num_rows, individual,
                num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    tsk_node_tbl_print_state(&table, _devnull);
    tsk_node_tbl_dump_text(&table, _devnull);

    /* Truncate back to the original number of rows. */
    ret = tsk_node_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.population, population, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.individual, individual, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    ret = tsk_node_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* If population is NULL it should be set to -1. If metadata is NULL all metadatas
     * should be set to the empty string. If individual is NULL it should be set to -1. */
    num_rows = 10;
    memset(population, 0xff, num_rows * sizeof(uint32_t));
    memset(individual, 0xff, num_rows * sizeof(uint32_t));
    ret = tsk_node_tbl_set_columns(&table, num_rows, flags, time, NULL, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.population, population, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.individual, individual, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* flags and time cannot be NULL */
    ret = tsk_node_tbl_set_columns(&table, num_rows, NULL, time, population, individual,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_tbl_set_columns(&table, num_rows, flags, NULL, population, individual,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_tbl_set_columns(&table, num_rows, flags, time, population, individual,
            NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_node_tbl_set_columns(&table, num_rows, flags, time, population, individual,
            metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* if metadata and metadata_offset are both null, all metadatas are zero length */
    num_rows = 10;
    memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_tbl_size_t));
    ret = tsk_node_tbl_set_columns(&table, num_rows, flags, time, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    ret = tsk_node_tbl_append_columns(&table, num_rows, flags, time, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.flags + num_rows, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset + num_rows, metadata_offset,
                num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    tsk_node_tbl_print_state(&table, _devnull);
    tsk_node_tbl_dump_text(&table, _devnull);

    tsk_node_tbl_free(&table);
    free(flags);
    free(population);
    free(time);
    free(metadata);
    free(metadata_offset);
    free(individual);
}

static void
test_edge_table(void)
{
    int ret;
    tsk_edge_tbl_t table;
    tsk_tbl_size_t num_rows = 100;
    tsk_id_t j;
    tsk_edge_t edge;
    tsk_id_t *parent, *child;
    double *left, *right;

    ret = tsk_edge_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_edge_tbl_set_max_rows_increment(&table, 1);
    tsk_edge_tbl_print_state(&table, _devnull);
    tsk_edge_tbl_dump_text(&table, _devnull);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret = tsk_edge_tbl_add_row(&table, (double) j, (double) j, j, j);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        CU_ASSERT_EQUAL(table.left[j], j);
        CU_ASSERT_EQUAL(table.right[j], j);
        CU_ASSERT_EQUAL(table.parent[j], j);
        CU_ASSERT_EQUAL(table.child[j], j);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
        ret = tsk_edge_tbl_get_row(&table, (tsk_tbl_size_t) j, &edge);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(edge.id, j);
        CU_ASSERT_EQUAL(edge.left, j);
        CU_ASSERT_EQUAL(edge.right, j);
        CU_ASSERT_EQUAL(edge.parent, j);
        CU_ASSERT_EQUAL(edge.child, j);
    }
    ret = tsk_edge_tbl_get_row(&table, num_rows, &edge);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    tsk_edge_tbl_print_state(&table, _devnull);
    tsk_edge_tbl_dump_text(&table, _devnull);

    num_rows *= 2;
    left = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    memset(left, 0, num_rows * sizeof(double));
    right = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    memset(right, 0, num_rows * sizeof(double));
    parent = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parent != NULL);
    memset(parent, 1, num_rows * sizeof(tsk_id_t));
    child = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(child != NULL);
    memset(child, 1, num_rows * sizeof(tsk_id_t));

    ret = tsk_edge_tbl_set_columns(&table, num_rows, left, right, parent, child);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    /* Append another num_rows to the end. */
    ret = tsk_edge_tbl_append_columns(&table, num_rows, left, right, parent, child);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.left + num_rows, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right + num_rows, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent + num_rows, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.child + num_rows, child, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_edge_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.child, child, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    ret = tsk_edge_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Inputs cannot be NULL */
    ret = tsk_edge_tbl_set_columns(&table, num_rows, NULL, right, parent, child);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_tbl_set_columns(&table, num_rows, left, NULL, parent, child);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_tbl_set_columns(&table, num_rows, left, right, NULL, child);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_edge_tbl_set_columns(&table, num_rows, left, right, parent, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    tsk_edge_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);

    tsk_edge_tbl_free(&table);
    free(left);
    free(right);
    free(parent);
    free(child);
}

static void
test_site_table(void)
{
    int ret;
    tsk_site_tbl_t table;
    tsk_tbl_size_t num_rows, j;
    char *ancestral_state;
    char *metadata;
    double *position;
    tsk_site_t site;
    tsk_tbl_size_t *ancestral_state_offset;
    tsk_tbl_size_t *metadata_offset;

    ret = tsk_site_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_site_tbl_set_max_rows_increment(&table, 1);
    tsk_site_tbl_set_max_metadata_length_increment(&table, 1);
    tsk_site_tbl_set_max_ancestral_state_length_increment(&table, 1);
    tsk_site_tbl_print_state(&table, _devnull);
    tsk_site_tbl_dump_text(&table, _devnull);

    ret = tsk_site_tbl_add_row(&table, 0, "A", 1, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(table.position[0], 0);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[0], 0);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 1);
    CU_ASSERT_EQUAL(table.metadata_offset[0], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[1], 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, 1);

    ret = tsk_site_tbl_get_row(&table, 0, &site);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(site.position, 0);
    CU_ASSERT_EQUAL(site.ancestral_state_length, 1);
    CU_ASSERT_NSTRING_EQUAL(site.ancestral_state, "A", 1);
    CU_ASSERT_EQUAL(site.metadata_length, 0);

    ret = tsk_site_tbl_add_row(&table, 1, "AA", 2, "{}", 2);
    CU_ASSERT_EQUAL_FATAL(ret, 1);
    CU_ASSERT_EQUAL(table.position[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[2], 3);
    CU_ASSERT_EQUAL(table.metadata_offset[1], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[2], 2);
    CU_ASSERT_EQUAL(table.metadata_length, 2);
    CU_ASSERT_EQUAL(table.num_rows, 2);

    ret = tsk_site_tbl_get_row(&table, 1, &site);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(site.position, 1);
    CU_ASSERT_EQUAL(site.ancestral_state_length, 2);
    CU_ASSERT_NSTRING_EQUAL(site.ancestral_state, "AA", 2);
    CU_ASSERT_EQUAL(site.metadata_length, 2);
    CU_ASSERT_NSTRING_EQUAL(site.metadata, "{}", 2);

    ret = tsk_site_tbl_add_row(&table, 2, "A", 1, "metadata", 8);
    CU_ASSERT_EQUAL_FATAL(ret, 2);
    CU_ASSERT_EQUAL(table.position[1], 1);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[3], 4);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 4);
    CU_ASSERT_EQUAL(table.metadata_offset[3], 10);
    CU_ASSERT_EQUAL(table.metadata_length, 10);
    CU_ASSERT_EQUAL(table.num_rows, 3);

    ret = tsk_site_tbl_get_row(&table, 3, &site);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_SITE_OUT_OF_BOUNDS);

    tsk_site_tbl_print_state(&table, _devnull);
    tsk_site_tbl_dump_text(&table, _devnull);
    tsk_site_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.ancestral_state_offset[0], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[0], 0);

    num_rows = 100;
    position = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(position != NULL);
    ancestral_state = malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(ancestral_state != NULL);
    ancestral_state_offset = malloc((num_rows + 1) * sizeof(uint32_t));
    CU_ASSERT_FATAL(ancestral_state_offset != NULL);
    metadata = malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = malloc((num_rows + 1) * sizeof(uint32_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);

    for (j = 0; j < num_rows; j++) {
        position[j] = (double) j;
        ancestral_state[j] = (char) j;
        ancestral_state_offset[j] = (tsk_tbl_size_t) j;
        metadata[j] = (char) ('A' + j);
        metadata_offset[j] = (tsk_tbl_size_t) j;
    }
    ancestral_state_offset[num_rows] = num_rows;
    metadata_offset[num_rows] = num_rows;

    ret = tsk_site_tbl_set_columns(&table, num_rows, position,
            ancestral_state, ancestral_state_offset,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.position, position,
                num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.ancestral_state, ancestral_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, num_rows);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    /* Append another num rows */
    ret = tsk_site_tbl_append_columns(&table, num_rows, position, ancestral_state,
            ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.position, position,
                num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.position + num_rows, position,
                num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.ancestral_state, ancestral_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.ancestral_state + num_rows, ancestral_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata + num_rows, metadata,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 2 * num_rows);

    /* truncate back to num_rows */
    ret = tsk_site_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.position, position,
                num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.ancestral_state, ancestral_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, num_rows);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    ret = tsk_site_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Inputs cannot be NULL */
    ret = tsk_site_tbl_set_columns(&table, num_rows, NULL, ancestral_state,
            ancestral_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_tbl_set_columns(&table, num_rows, position, NULL, ancestral_state_offset,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_tbl_set_columns(&table, num_rows, position, ancestral_state, NULL,
            metadata, metadata_offset);
    /* Metadata and metadata_offset must both be null */
    ret = tsk_site_tbl_set_columns(&table, num_rows, position, ancestral_state,
            ancestral_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_site_tbl_set_columns(&table, num_rows, position, ancestral_state,
            ancestral_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Set metadata to NULL */
    ret = tsk_site_tbl_set_columns(&table, num_rows, position,
            ancestral_state, ancestral_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    memset(metadata_offset, 0, (num_rows + 1) * sizeof(uint32_t));
    CU_ASSERT_EQUAL(memcmp(table.position, position,
                num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.ancestral_state, ancestral_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, num_rows);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    /* Test for bad offsets */
    ancestral_state_offset[0] = 1;
    ret = tsk_site_tbl_set_columns(&table, num_rows, position,
            ancestral_state, ancestral_state_offset,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    ancestral_state_offset[0] = 0;
    ancestral_state_offset[num_rows] = 0;
    ret = tsk_site_tbl_set_columns(&table, num_rows, position,
            ancestral_state, ancestral_state_offset,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    ancestral_state_offset[0] = 0;

    metadata_offset[0] = 0;
    ret = tsk_site_tbl_set_columns(&table, num_rows, position,
            ancestral_state, ancestral_state_offset,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    metadata_offset[0] = 0;
    metadata_offset[num_rows] = 0;
    ret = tsk_site_tbl_set_columns(&table, num_rows, position,
            ancestral_state, ancestral_state_offset,
            metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    ret = tsk_site_tbl_clear(&table);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.ancestral_state_length, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_site_tbl_free(&table);
    free(position);
    free(ancestral_state);
    free(ancestral_state_offset);
    free(metadata);
    free(metadata_offset);
}

static void
test_mutation_table(void)
{
    int ret;
    tsk_mutation_tbl_t table;
    tsk_tbl_size_t num_rows = 100;
    tsk_tbl_size_t max_len = 20;
    tsk_tbl_size_t k, len;
    tsk_id_t j;
    tsk_id_t *node;
    tsk_id_t *parent;
    tsk_id_t *site;
    char *derived_state, *metadata;
    char c[max_len + 1];
    tsk_tbl_size_t *derived_state_offset, *metadata_offset;
    tsk_mutation_t mutation;

    for (j = 0; j < (tsk_id_t) max_len; j++) {
        c[j] = (char) ('A' + j);
    }

    ret = tsk_mutation_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_mutation_tbl_set_max_rows_increment(&table, 1);
    tsk_mutation_tbl_set_max_metadata_length_increment(&table, 1);
    tsk_mutation_tbl_set_max_derived_state_length_increment(&table, 1);
    tsk_mutation_tbl_print_state(&table, _devnull);
    tsk_mutation_tbl_dump_text(&table, _devnull);

    len = 0;
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        k = TSK_MIN((tsk_tbl_size_t) j + 1, max_len);
        ret = tsk_mutation_tbl_add_row(&table, j, j, j, c, k, c, k);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        CU_ASSERT_EQUAL(table.site[j], j);
        CU_ASSERT_EQUAL(table.node[j], j);
        CU_ASSERT_EQUAL(table.parent[j], j);
        CU_ASSERT_EQUAL(table.derived_state_offset[j], len);
        CU_ASSERT_EQUAL(table.metadata_offset[j], len);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
        len += k;
        CU_ASSERT_EQUAL(table.derived_state_offset[j + 1], len);
        CU_ASSERT_EQUAL(table.derived_state_length, len);
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], len);
        CU_ASSERT_EQUAL(table.metadata_length, len);

        ret = tsk_mutation_tbl_get_row(&table, (tsk_tbl_size_t) j, &mutation);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(mutation.id, j);
        CU_ASSERT_EQUAL(mutation.site, j);
        CU_ASSERT_EQUAL(mutation.node, j);
        CU_ASSERT_EQUAL(mutation.parent, j);
        CU_ASSERT_EQUAL(mutation.metadata_length, k);
        CU_ASSERT_NSTRING_EQUAL(mutation.metadata, c, k);
        CU_ASSERT_EQUAL(mutation.derived_state_length, k);
        CU_ASSERT_NSTRING_EQUAL(mutation.derived_state, c, k);
    }
    ret = tsk_mutation_tbl_get_row(&table, num_rows, &mutation);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MUTATION_OUT_OF_BOUNDS);
    tsk_mutation_tbl_print_state(&table, _devnull);
    tsk_mutation_tbl_dump_text(&table, _devnull);

    num_rows *= 2;
    site = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(site != NULL);
    node = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    parent = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(parent != NULL);
    derived_state = malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(derived_state != NULL);
    derived_state_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(derived_state_offset != NULL);
    metadata = malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        node[j] = j;
        site[j] = j + 1;
        parent[j] = j + 2;
        derived_state[j] = 'Y';
        derived_state_offset[j] = (tsk_tbl_size_t) j;
        metadata[j] = 'M';
        metadata_offset[j] = (tsk_tbl_size_t) j;
    }

    derived_state_offset[num_rows] = num_rows;
    metadata_offset[num_rows] = num_rows;
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, parent,
            derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* Append another num_rows */
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, site, node, parent, derived_state,
            derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.site + num_rows, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node + num_rows, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent + num_rows, parent,
                num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.derived_state_length, 2 * num_rows);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_mutation_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    ret = tsk_mutation_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Check all this again, except with parent == NULL and metadata == NULL. */
    memset(parent, 0xff, num_rows * sizeof(tsk_id_t));
    memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_tbl_size_t));
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, NULL,
            derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state_offset, derived_state_offset,
                num_rows * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    /* Append another num_rows */
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, site, node, NULL, derived_state,
            derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.site, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.site + num_rows, site, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node + num_rows, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent, parent, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.parent + num_rows, parent,
                num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.derived_state + num_rows, derived_state,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.derived_state_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);


    /* Inputs except parent, metadata and metadata_offset cannot be NULL*/
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, NULL, node, parent,
            derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, NULL, parent,
            derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, parent,
            NULL, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, parent,
            derived_state, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, parent,
            derived_state, derived_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, parent,
            derived_state, derived_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Inputs except parent, metadata and metadata_offset cannot be NULL*/
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, NULL, node, parent,
            derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, site, NULL, parent,
            derived_state, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, site, node, parent,
            NULL, derived_state_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, site, node, parent,
            derived_state, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, site, node, parent,
            derived_state, derived_state_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_mutation_tbl_append_columns(&table, num_rows, site, node, parent,
            derived_state, derived_state_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Test for bad offsets */
    derived_state_offset[0] = 1;
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, parent,
            derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    derived_state_offset[0] = 0;
    derived_state_offset[num_rows] = 0;
    ret = tsk_mutation_tbl_set_columns(&table, num_rows, site, node, parent,
            derived_state, derived_state_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    tsk_mutation_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.derived_state_length, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_mutation_tbl_free(&table);
    free(site);
    free(node);
    free(parent);
    free(derived_state);
    free(derived_state_offset);
    free(metadata);
    free(metadata_offset);
}

static void
test_migration_table(void)
{
    int ret;
    tsk_migration_tbl_t table;
    tsk_tbl_size_t num_rows = 100;
    tsk_id_t j;
    tsk_id_t *node;
    tsk_id_t *source, *dest;
    double *left, *right, *time;
    tsk_migration_t migration;

    ret = tsk_migration_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_migration_tbl_set_max_rows_increment(&table, 1);
    tsk_migration_tbl_print_state(&table, _devnull);
    tsk_migration_tbl_dump_text(&table, _devnull);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret = tsk_migration_tbl_add_row(&table, j, j, j, j, j, j);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        CU_ASSERT_EQUAL(table.left[j], j);
        CU_ASSERT_EQUAL(table.right[j], j);
        CU_ASSERT_EQUAL(table.node[j], j);
        CU_ASSERT_EQUAL(table.source[j], j);
        CU_ASSERT_EQUAL(table.dest[j], j);
        CU_ASSERT_EQUAL(table.time[j], j);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);

        ret = tsk_migration_tbl_get_row(&table, (tsk_tbl_size_t) j, &migration);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(migration.id, j);
        CU_ASSERT_EQUAL(migration.left, j);
        CU_ASSERT_EQUAL(migration.right, j);
        CU_ASSERT_EQUAL(migration.node, j);
        CU_ASSERT_EQUAL(migration.source, j);
        CU_ASSERT_EQUAL(migration.dest, j);
        CU_ASSERT_EQUAL(migration.time, j);
    }
    ret = tsk_migration_tbl_get_row(&table, num_rows, &migration);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_MIGRATION_OUT_OF_BOUNDS);
    tsk_migration_tbl_print_state(&table, _devnull);
    tsk_migration_tbl_dump_text(&table, _devnull);

    num_rows *= 2;
    left = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(left != NULL);
    memset(left, 1, num_rows * sizeof(double));
    right = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(right != NULL);
    memset(right, 2, num_rows * sizeof(double));
    time = malloc(num_rows * sizeof(double));
    CU_ASSERT_FATAL(time != NULL);
    memset(time, 3, num_rows * sizeof(double));
    node = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(node != NULL);
    memset(node, 4, num_rows * sizeof(tsk_id_t));
    source = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(source != NULL);
    memset(source, 5, num_rows * sizeof(tsk_id_t));
    dest = malloc(num_rows * sizeof(tsk_id_t));
    CU_ASSERT_FATAL(dest != NULL);
    memset(dest, 6, num_rows * sizeof(tsk_id_t));

    ret = tsk_migration_tbl_set_columns(&table, num_rows, left, right, node, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    /* Append another num_rows */
    ret = tsk_migration_tbl_append_columns(&table, num_rows, left, right, node, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.left + num_rows, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right + num_rows, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time + num_rows, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node + num_rows, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.source + num_rows, source,
                num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.dest + num_rows, dest,
                num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_migration_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.left, left, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.right, right, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.time, time, num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.node, node, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.source, source, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.dest, dest, num_rows * sizeof(tsk_id_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);

    ret = tsk_migration_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* inputs cannot be NULL */
    ret = tsk_migration_tbl_set_columns(&table, num_rows, NULL, right, node, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_tbl_set_columns(&table, num_rows, left, NULL, node, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_tbl_set_columns(&table, num_rows, left, right, NULL, source,
            dest, time);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_tbl_set_columns(&table, num_rows, left, right, node, NULL,
            dest, time);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_tbl_set_columns(&table, num_rows, left, right, node, source,
            NULL, time);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_migration_tbl_set_columns(&table, num_rows, left, right, node, source,
            dest, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    tsk_migration_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);

    tsk_migration_tbl_free(&table);
    free(left);
    free(right);
    free(time);
    free(node);
    free(source);
    free(dest);
}

static void
test_individual_table(void)
{
    int ret = 0;
    tsk_individual_tbl_t table;
    /* tsk_tbl_collection_t tables, tables2; */
    tsk_tbl_size_t num_rows = 100;
    tsk_id_t j;
    tsk_tbl_size_t k;
    uint32_t *flags;
    double *location;
    char *metadata;
    tsk_tbl_size_t *metadata_offset;
    tsk_tbl_size_t *location_offset;
    tsk_individual_t individual;
    const char *test_metadata = "test";
    tsk_tbl_size_t test_metadata_length = 4;
    char metadata_copy[test_metadata_length + 1];
    tsk_tbl_size_t spatial_dimension = 2;
    double test_location[spatial_dimension];

    for (k = 0; k < spatial_dimension; k++) {
        test_location[k] = (double) k;
    }
    metadata_copy[test_metadata_length] = '\0';
    ret = tsk_individual_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_individual_tbl_set_max_rows_increment(&table, 1);
    tsk_individual_tbl_set_max_metadata_length_increment(&table, 1);
    tsk_individual_tbl_set_max_location_length_increment(&table, 1);

    tsk_individual_tbl_print_state(&table, _devnull);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        ret = tsk_individual_tbl_add_row(&table, (uint32_t) j, test_location,
                (size_t) spatial_dimension, test_metadata, test_metadata_length);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        CU_ASSERT_EQUAL(table.flags[j], j);
        for (k = 0; k < spatial_dimension; k++) {
            test_location[k] = (double) k;
            CU_ASSERT_EQUAL(table.location[spatial_dimension * (size_t) j + k],
                    test_location[k]);
        }
        CU_ASSERT_EQUAL(table.metadata_length, (tsk_tbl_size_t) (j + 1) * test_metadata_length);
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], table.metadata_length);
        /* check the metadata */
        memcpy(metadata_copy, table.metadata + table.metadata_offset[j],
                test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(metadata_copy, test_metadata, test_metadata_length);

        ret = tsk_individual_tbl_get_row(&table, (tsk_tbl_size_t) j, &individual);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(individual.id, j);
        CU_ASSERT_EQUAL(individual.flags, j);
        CU_ASSERT_EQUAL(individual.location_length, spatial_dimension);
        CU_ASSERT_NSTRING_EQUAL(individual.location, test_location,
                spatial_dimension * sizeof(double));
        CU_ASSERT_EQUAL(individual.metadata_length, test_metadata_length);
        CU_ASSERT_NSTRING_EQUAL(individual.metadata, test_metadata, test_metadata_length);
    }
    ret = tsk_individual_tbl_get_row(&table, num_rows, &individual);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS);
    tsk_individual_tbl_print_state(&table, _devnull);
    tsk_individual_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    num_rows *= 2;
    flags = malloc(num_rows * sizeof(uint32_t));
    CU_ASSERT_FATAL(flags != NULL);
    memset(flags, 1, num_rows * sizeof(uint32_t));
    location = malloc(spatial_dimension * num_rows * sizeof(double));
    CU_ASSERT_FATAL(location != NULL);
    memset(location, 0, spatial_dimension * num_rows * sizeof(double));
    location_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(location_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        location_offset[j] = (tsk_tbl_size_t) j * spatial_dimension;
    }
    metadata = malloc(num_rows * sizeof(char));
    memset(metadata, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);
    for (j = 0; j < (tsk_id_t) num_rows + 1; j++) {
        metadata_offset[j] = (tsk_tbl_size_t) j;
    }
    ret = tsk_individual_tbl_set_columns(&table, num_rows, flags,
            location, location_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location, location,
                spatial_dimension * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location_offset, location_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.location_length, spatial_dimension * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    tsk_individual_tbl_print_state(&table, _devnull);

    /* Append another num_rows onto the end */
    ret = tsk_individual_tbl_append_columns(&table, num_rows, flags, location,
            location_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.flags + num_rows, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata + num_rows, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location, location,
                spatial_dimension * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location + spatial_dimension * num_rows,
                location, spatial_dimension * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    tsk_individual_tbl_print_state(&table, _devnull);
    tsk_individual_tbl_dump_text(&table, _devnull);

    /* Truncate back to num_rows */
    ret = tsk_individual_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location, location,
                spatial_dimension * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location_offset, location_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.location_length, spatial_dimension * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);
    tsk_individual_tbl_print_state(&table, _devnull);

    ret = tsk_individual_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* flags can't be NULL */
    ret = tsk_individual_tbl_set_columns(&table, num_rows, NULL,
            location, location_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* location and location offset must be simultaneously NULL or not */
    ret = tsk_individual_tbl_set_columns(&table, num_rows, flags,
            location, NULL, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_tbl_set_columns(&table, num_rows, flags,
            NULL, location_offset, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    /* metadata and metadata offset must be simultaneously NULL or not */
    ret = tsk_individual_tbl_set_columns(&table, num_rows, flags,
            location, location_offset, NULL, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_individual_tbl_set_columns(&table, num_rows, flags,
            location, location_offset, metadata, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* if location and location_offset are both null, all locations are zero length */
    num_rows = 10;
    memset(location_offset, 0, (num_rows + 1) * sizeof(tsk_tbl_size_t));
    ret = tsk_individual_tbl_set_columns(&table, num_rows, flags,
            NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.location_offset, location_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.location_length, 0);
    ret = tsk_individual_tbl_append_columns(&table, num_rows, flags, NULL, NULL, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.location_offset, location_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location_offset + num_rows, location_offset,
                num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.location_length, 0);
    tsk_individual_tbl_print_state(&table, _devnull);
    tsk_individual_tbl_dump_text(&table, _devnull);

    /* if metadata and metadata_offset are both null, all metadatas are zero length */
    num_rows = 10;
    memset(metadata_offset, 0, (num_rows + 1) * sizeof(tsk_tbl_size_t));
    ret = tsk_individual_tbl_set_columns(&table, num_rows, flags,
            location, location_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.flags, flags, num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location, location,
                spatial_dimension * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    ret = tsk_individual_tbl_append_columns(&table, num_rows, flags, location,
            location_offset, NULL, NULL);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.location, location,
                spatial_dimension * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.location + spatial_dimension * num_rows,
                location, spatial_dimension * num_rows * sizeof(double)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset, metadata_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata_offset + num_rows, metadata_offset,
                num_rows * sizeof(uint32_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, 0);
    tsk_individual_tbl_print_state(&table, _devnull);
    tsk_individual_tbl_dump_text(&table, _devnull);

    ret = tsk_individual_tbl_free(&table);
    CU_ASSERT_EQUAL(ret, 0);
    free(flags);
    free(location);
    free(location_offset);
    free(metadata);
    free(metadata_offset);
}

static void
test_population_table(void)
{
    int ret;
    tsk_population_tbl_t table;
    tsk_tbl_size_t num_rows = 100;
    tsk_tbl_size_t max_len = 20;
    tsk_tbl_size_t k, len;
    tsk_id_t j;
    char *metadata;
    char c[max_len + 1];
    tsk_tbl_size_t *metadata_offset;
    tsk_population_t population;

    for (j = 0; j < (tsk_id_t) max_len; j++) {
        c[j] = (char) ('A' + j);
    }

    ret = tsk_population_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_population_tbl_set_max_rows_increment(&table, 1);
    tsk_population_tbl_set_max_metadata_length_increment(&table, 1);
    tsk_population_tbl_print_state(&table, _devnull);
    tsk_population_tbl_dump_text(&table, _devnull);
    /* Adding zero length metadata with NULL should be fine */

    ret = tsk_population_tbl_add_row(&table, NULL, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_EQUAL(table.metadata_length, 0);
    CU_ASSERT_EQUAL(table.num_rows, 1);
    CU_ASSERT_EQUAL(table.metadata_offset[0], 0);
    CU_ASSERT_EQUAL(table.metadata_offset[1], 0);
    tsk_population_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);

    len = 0;
    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        k = TSK_MIN((tsk_tbl_size_t) j + 1, max_len);
        ret = tsk_population_tbl_add_row(&table, c, k);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        CU_ASSERT_EQUAL(table.metadata_offset[j], len);
        CU_ASSERT_EQUAL(table.num_rows, j + 1);
        len += k;
        CU_ASSERT_EQUAL(table.metadata_offset[j + 1], len);
        CU_ASSERT_EQUAL(table.metadata_length, len);

        ret = tsk_population_tbl_get_row(&table, (tsk_tbl_size_t) j, &population);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(population.id, j);
        CU_ASSERT_EQUAL(population.metadata_length, k);
        CU_ASSERT_NSTRING_EQUAL(population.metadata, c, k);
    }
    ret = tsk_population_tbl_get_row(&table, num_rows, &population);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_POPULATION_OUT_OF_BOUNDS);
    tsk_population_tbl_print_state(&table, _devnull);
    tsk_population_tbl_dump_text(&table, _devnull);

    num_rows *= 2;
    metadata = malloc(num_rows * sizeof(char));
    CU_ASSERT_FATAL(metadata != NULL);
    metadata_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(metadata_offset != NULL);

    for (j = 0; j < (tsk_id_t) num_rows; j++) {
        metadata[j] = 'M';
        metadata_offset[j] = (tsk_tbl_size_t) j;
    }

    metadata_offset[num_rows] = num_rows;
    ret = tsk_population_tbl_set_columns(&table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    /* Append another num_rows */
    ret = tsk_population_tbl_append_columns(&table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata + num_rows, metadata,
                num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.metadata_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);

    /* Truncate back to num_rows */
    ret = tsk_population_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.metadata, metadata, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.metadata_length, num_rows);

    ret = tsk_population_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* Metadata = NULL gives an error */
    ret = tsk_population_tbl_set_columns(&table, num_rows, NULL, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_population_tbl_set_columns(&table, num_rows, metadata, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_population_tbl_set_columns(&table, num_rows, NULL, metadata_offset);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    /* Test for bad offsets */
    metadata_offset[0] = 1;
    ret = tsk_population_tbl_set_columns(&table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);
    metadata_offset[0] = 0;
    metadata_offset[num_rows] = 0;
    ret = tsk_population_tbl_set_columns(&table, num_rows, metadata, metadata_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_OFFSET);

    tsk_population_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.metadata_length, 0);

    tsk_population_tbl_free(&table);
    free(metadata);
    free(metadata_offset);
}

static void
test_provenance_table(void)
{
    int ret;
    tsk_provenance_tbl_t table;
    tsk_tbl_size_t num_rows = 100;
    tsk_tbl_size_t j;
    char *timestamp;
    uint32_t *timestamp_offset;
    const char *test_timestamp = "2017-12-06T20:40:25+00:00";
    size_t test_timestamp_length = strlen(test_timestamp);
    char timestamp_copy[test_timestamp_length + 1];
    char *record;
    uint32_t *record_offset;
    const char *test_record = "{\"json\"=1234}";
    size_t test_record_length = strlen(test_record);
    char record_copy[test_record_length + 1];
    tsk_provenance_t provenance;

    timestamp_copy[test_timestamp_length] = '\0';
    record_copy[test_record_length] = '\0';
    ret = tsk_provenance_tbl_alloc(&table, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_provenance_tbl_set_max_rows_increment(&table, 1);
    tsk_provenance_tbl_set_max_timestamp_length_increment(&table, 1);
    tsk_provenance_tbl_set_max_record_length_increment(&table, 1);
    tsk_provenance_tbl_print_state(&table, _devnull);
    tsk_provenance_tbl_dump_text(&table, _devnull);

    for (j = 0; j < num_rows; j++) {
        ret = tsk_provenance_tbl_add_row(&table, test_timestamp, test_timestamp_length,
                test_record, test_record_length);
        CU_ASSERT_EQUAL_FATAL(ret, j);
        CU_ASSERT_EQUAL(table.timestamp_length, (j + 1) * test_timestamp_length);
        CU_ASSERT_EQUAL(table.timestamp_offset[j + 1], table.timestamp_length);
        CU_ASSERT_EQUAL(table.record_length, (j + 1) * test_record_length);
        CU_ASSERT_EQUAL(table.record_offset[j + 1], table.record_length);
        /* check the timestamp */
        memcpy(timestamp_copy, table.timestamp + table.timestamp_offset[j],
                test_timestamp_length);
        CU_ASSERT_NSTRING_EQUAL(timestamp_copy, test_timestamp, test_timestamp_length);
        /* check the record */
        memcpy(record_copy, table.record + table.record_offset[j],
                test_record_length);
        CU_ASSERT_NSTRING_EQUAL(record_copy, test_record, test_record_length);

        ret = tsk_provenance_tbl_get_row(&table, j, &provenance);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_EQUAL(provenance.id, j);
        CU_ASSERT_EQUAL(provenance.timestamp_length, test_timestamp_length);
        CU_ASSERT_NSTRING_EQUAL(provenance.timestamp, test_timestamp,
                test_timestamp_length);
        CU_ASSERT_EQUAL(provenance.record_length, test_record_length);
        CU_ASSERT_NSTRING_EQUAL(provenance.record, test_record,
                test_record_length);
    }
    ret = tsk_provenance_tbl_get_row(&table, num_rows, &provenance);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_PROVENANCE_OUT_OF_BOUNDS);
    tsk_provenance_tbl_print_state(&table, _devnull);
    tsk_provenance_tbl_dump_text(&table, _devnull);
    tsk_provenance_tbl_clear(&table);
    CU_ASSERT_EQUAL(table.num_rows, 0);
    CU_ASSERT_EQUAL(table.timestamp_length, 0);
    CU_ASSERT_EQUAL(table.record_length, 0);

    num_rows *= 2;
    timestamp = malloc(num_rows * sizeof(char));
    memset(timestamp, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(timestamp != NULL);
    timestamp_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(timestamp_offset != NULL);
    record = malloc(num_rows * sizeof(char));
    memset(record, 'a', num_rows * sizeof(char));
    CU_ASSERT_FATAL(record != NULL);
    record_offset = malloc((num_rows + 1) * sizeof(tsk_tbl_size_t));
    CU_ASSERT_FATAL(record_offset != NULL);
    for (j = 0; j < num_rows + 1; j++) {
        timestamp_offset[j] = j;
        record_offset[j] = j;
    }
    ret = tsk_provenance_tbl_set_columns(&table, num_rows,
            timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.timestamp, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.timestamp_offset, timestamp_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.record, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.record_offset, record_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.timestamp_length, num_rows);
    CU_ASSERT_EQUAL(table.record_length, num_rows);
    tsk_provenance_tbl_print_state(&table, _devnull);

    /* Append another num_rows onto the end */
    ret = tsk_provenance_tbl_append_columns(&table, num_rows,
            timestamp, timestamp_offset, record, record_offset);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.timestamp, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.timestamp + num_rows, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.record, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.record + num_rows, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(table.num_rows, 2 * num_rows);
    CU_ASSERT_EQUAL(table.timestamp_length, 2 * num_rows);
    CU_ASSERT_EQUAL(table.record_length, 2 * num_rows);
    tsk_provenance_tbl_print_state(&table, _devnull);

    /* Truncate back to num_rows */
    ret = tsk_provenance_tbl_truncate(&table, num_rows);
    CU_ASSERT_EQUAL(ret, 0);
    CU_ASSERT_EQUAL(memcmp(table.timestamp, timestamp, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.timestamp_offset, timestamp_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(memcmp(table.record, record, num_rows * sizeof(char)), 0);
    CU_ASSERT_EQUAL(memcmp(table.record_offset, record_offset,
                (num_rows + 1) * sizeof(tsk_tbl_size_t)), 0);
    CU_ASSERT_EQUAL(table.num_rows, num_rows);
    CU_ASSERT_EQUAL(table.timestamp_length, num_rows);
    CU_ASSERT_EQUAL(table.record_length, num_rows);
    tsk_provenance_tbl_print_state(&table, _devnull);

    ret = tsk_provenance_tbl_truncate(&table, num_rows + 1);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_TABLE_POSITION);

    /* No arguments can be null */
    ret = tsk_provenance_tbl_set_columns(&table, num_rows, NULL, timestamp_offset,
            record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_tbl_set_columns(&table, num_rows, timestamp, NULL,
            record, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_tbl_set_columns(&table, num_rows, timestamp, timestamp_offset,
            NULL, record_offset);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);
    ret = tsk_provenance_tbl_set_columns(&table, num_rows, timestamp, timestamp_offset,
            record, NULL);
    CU_ASSERT_EQUAL(ret, TSK_ERR_BAD_PARAM_VALUE);

    tsk_provenance_tbl_free(&table);
    free(timestamp);
    free(timestamp_offset);
    free(record);
    free(record_offset);
}

static void
test_simplify_tables_drops_indexes(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;
    tsk_id_t samples[] = {0, 1};

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges,
            NULL, NULL, NULL, NULL, NULL);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_tbl_collection_is_indexed(&tables))
    ret = tsk_tbl_collection_simplify(&tables, samples, 2, 0, NULL);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_tbl_collection_is_indexed(&tables))

    tsk_tbl_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

static void
test_sort_tables_drops_indexes(void)
{
    int ret;
    tsk_treeseq_t ts;
    tsk_tbl_collection_t tables;

    tsk_treeseq_from_text(&ts, 1, single_tree_ex_nodes, single_tree_ex_edges,
            NULL, NULL, NULL, NULL, NULL);
    ret = tsk_tbl_collection_alloc(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_copy_tables(&ts, &tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_tbl_collection_is_indexed(&tables))
    ret = tsk_tbl_collection_sort(&tables, 0, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_tbl_collection_is_indexed(&tables))

    tsk_tbl_collection_free(&tables);
    tsk_treeseq_free(&ts);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        {"test_node_table", test_node_table},
        {"test_edge_table", test_edge_table},
        {"test_site_table", test_site_table},
        {"test_mutation_table", test_mutation_table},
        {"test_migration_table", test_migration_table},
        {"test_individual_table", test_individual_table},
        {"test_population_table", test_population_table},
        {"test_provenance_table", test_provenance_table},
        {"test_format_data_load_errors", test_format_data_load_errors},
        {"test_dump_unindexed", test_dump_unindexed},
        {"test_tbl_collection_load_errors", test_tbl_collection_load_errors},
        {"test_tbl_collection_dump_errors", test_tbl_collection_dump_errors},
        {"test_tbl_collection_simplify_errors", test_tbl_collection_simplify_errors},
        {"test_load_tsk_node_tbl_errors", test_load_tsk_node_tbl_errors},
        {"test_simplify_tables_drops_indexes", test_simplify_tables_drops_indexes},
        {"test_sort_tables_drops_indexes", test_sort_tables_drops_indexes},
        {NULL},
    };

    return test_main(tests, argc, argv);
}
